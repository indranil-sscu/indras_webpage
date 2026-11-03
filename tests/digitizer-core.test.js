// Unit tests for assets/js/digitizer-core.js — run with: node --test tests/*.test.js
const test = require("node:test");
const assert = require("node:assert");
const D = require("../assets/js/digitizer-core.js");

const near = (got, exp, tol) => assert.ok(Math.abs(got - exp) <= tol, `${got} != ${exp} ± ${tol}`);

// Synthetic "image": data (x, y) -> pixel via rotation, anisotropic scale and offset
function makeImage(angleDeg, sx, sy, ox, oy, log) {
  const th = (angleDeg * Math.PI) / 180, C = Math.cos(th), S = Math.sin(th);
  return (x, y) => {
    const X = log && log.x ? Math.log10(x) : x, Y = log && log.y ? Math.log10(y) : y;
    const u0 = sx * X, v0 = -sy * Y; // image y grows downwards
    return [ox + C * u0 - S * v0, oy + S * u0 + C * v0];
  };
}

test("axis-aligned linear plot", () => {
  const img = makeImage(0, 50, 30, 100, 400);
  const cal = D.calibrate({ X1: img(0, 0), X2: img(10, 0), Y1: img(0, 0), Y2: img(0, 5) }, { x1: 0, x2: 10, y1: 0, y2: 5 });
  const [x, y] = D.toData(cal, ...img(3.3, 2.2));
  near(x, 3.3, 1e-9); near(y, 2.2, 1e-9);
});

test("rotated and sheared scan, axis points away from the origin", () => {
  const img = makeImage(3.5, 42, 37, 120, 610);
  // calibration ticks not at the origin: x-axis drawn at y = 1, y-axis at x = 400
  const cal = D.calibrate(
    { X1: img(400, 1), X2: img(800, 1), Y1: img(400, 1), Y2: img(400, 9) },
    { x1: 400, x2: 800, y1: 1, y2: 9 }
  );
  for (const [x, y] of [[550, 3], [712.5, 7.25], [450, 8.9]]) {
    const [gx, gy] = D.toData(cal, ...img(x, y));
    near(gx, x, 1e-8); near(gy, y, 1e-8);
  }
});

test("log-log axes", () => {
  const log = { x: true, y: true };
  const img = makeImage(0, 120, 90, 80, 500, log);
  const cal = D.calibrate({ X1: img(1, 1e-3), X2: img(1000, 1e-3), Y1: img(1, 1e-3), Y2: img(1, 10) },
    { x1: 1, x2: 1000, y1: 1e-3, y2: 10 }, log);
  const [x, y] = D.toData(cal, ...img(37, 0.042));
  near(x, 37, 1e-9); near(y, 0.042, 1e-12);
});

test("toPixel inverts toData", () => {
  const img = makeImage(-2, 10, 20, 5, 300, { y: true });
  const cal = D.calibrate({ X1: img(0, 1), X2: img(50, 1), Y1: img(0, 1), Y2: img(0, 100) }, { x1: 0, x2: 50, y1: 1, y2: 100 }, { y: true });
  const p = img(17, 33);
  const back = D.toPixel(cal, ...D.toData(cal, ...p));
  near(back[0], p[0], 1e-9); near(back[1], p[1], 1e-9);
});

test("invalid calibrations are rejected", () => {
  const p = { X1: [0, 0], X2: [10, 0], Y1: [0, 0], Y2: [0, -10] };
  assert.throws(() => D.calibrate(p, { x1: 1, x2: 1, y1: 0, y2: 1 }));
  assert.throws(() => D.calibrate(p, { x1: 0, x2: 1, y1: 0, y2: 1 }, { x: true }));
  assert.throws(() => D.calibrate({ X1: [0, 0], X2: [0, 0], Y1: [0, 0], Y2: [0, -10] }, { x1: 0, x2: 1, y1: 0, y2: 1 }));
});

test("CSV export", () => {
  const pts = [{ series: "A", x: 2, y: 0.5 }, { series: "A", x: 1, y: 0.25 }];
  assert.strictEqual(D.toCSV(pts, { sortByX: true }), "x,y\n1,0.25\n2,0.5\n");
  const multi = D.toCSV([...pts, { series: "B, run 2", x: 1.5, y: 3 }], { sortByX: true });
  assert.strictEqual(multi, 'series,x,y\nA,1,0.25\nA,2,0.5\n"B, run 2",1.5,3\n');
});

test("nearest point lookup", () => {
  assert.strictEqual(D.nearest([[0, 0], [10, 10], [12, 9]], 11.5, 9.2, 5), 2);
  assert.strictEqual(D.nearest([[0, 0]], 50, 50, 5), -1);
});
