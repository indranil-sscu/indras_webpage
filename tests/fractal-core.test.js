// Tests for assets/js/fractal-core.js against analytical results. Run: node --test tests/*.test.js
const test = require("node:test");
const assert = require("node:assert");
const F = require("../assets/js/fractal-core.js");

const near = (got, exp, tol) => assert.ok(Math.abs(got - exp) <= tol, `${got} != ${exp} ± ${tol}`);

test("escape counts for simple points", () => {
  assert.strictEqual(F.escape(0, 0, 1, 0, 100), 3); // 0, 1, 2, 5: |z₃| > 2
  assert.strictEqual(F.escape(0, 0, -2, 0, 1000), -1); // 0, −2, 2, 2, …: stays on |z| = 2
  assert.strictEqual(F.escape(0, 0, 0.26, 0, 1000) > 0, true); // just right of the cusp at ¼
  assert.strictEqual(F.mandelbrot(0, 0, 1000), -1);
  assert.strictEqual(F.mandelbrot(-1, 0, 1000), -1); // period-2 cycle 0, −1, 0, −1
  assert.strictEqual(F.escape(0, 0, 0, 1, 1000), -1); // c = i: 0, i, −1 + i, −i, −1 + i, … (pre-periodic)
  assert.strictEqual(F.escape(3, 0, 0, 0, 10), 0); // already outside
});

test("interior shortcut agrees with brute-force iteration", () => {
  let checked = 0;
  for (let i = 0; i < 120; i++) for (let j = 0; j < 80; j++) {
    const x = -2 + (2.6 * i) / 119, y = (1.3 * j) / 79;
    if (F.inMainBulbs(x, y)) { checked++; assert.strictEqual(F.escape(0, 0, x, y, 3000), -1, `${x},${y}`); }
  }
  assert.ok(checked > 300);
  assert.ok(F.inMainBulbs(0.2, 0) && F.inMainBulbs(-1.2, 0) && !F.inMainBulbs(0.3, 0) && !F.inMainBulbs(-0.6, 0.6));
});

test("Julia set of c = 0 is the unit disc, with exact escape times", () => {
  for (const r of [1.01, 1.3, 1.9, 2.5]) {
    // z_n = z0^(2^n): first n with r^(2^n) > R is ceil(log2(ln R / ln r))
    const R = 2, n = Math.ceil(Math.log2(Math.log(R) / Math.log(r)));
    assert.strictEqual(F.julia(r * Math.cos(1), r * Math.sin(1), 0, 0, 1000, R), Math.max(0, n));
    // smooth count: n − log2(ln|z_n| / ln R) = log2(ln R / ln r), whatever the bailout
    for (const RR of [2, 16, 256]) if (r < RR) near(F.julia(r, 0, 0, 0, 1000, RR, true), Math.log2(Math.log(RR) / Math.log(r)), 1e-9);
  }
  assert.strictEqual(F.julia(0.99, 0, 0, 0, 2000), -1);
  assert.strictEqual(F.julia(0, 0.5, 0, 0, 2000), -1);
});

test("smooth count is continuous across iteration bands", () => {
  // along a ray just outside the set, neighbouring points never jump by more than a small amount
  let prev = null, maxJump = 0;
  for (let k = 0; k <= 2000; k++) {
    const x = 0.26 + (0.5 * k) / 2000, v = F.mandelbrot(x, 0, 5000, 256, true);
    if (prev !== null) maxJump = Math.max(maxJump, Math.abs(v - prev));
    prev = v;
  }
  assert.ok(maxJump < 0.5, `max jump ${maxJump}`);
});

test("orbits and rendering", () => {
  const o = F.orbit(0, 0, 1, 0, 50);
  assert.deepStrictEqual(o, [[0, 0], [1, 0], [2, 0], [5, 0]]);
  assert.strictEqual(F.orbit(0, 0, -1, 0, 10).length, 11);
  const w = 40, h = 30, buf = new Float32Array(w * h);
  const esc = F.renderRows(buf, w, h, { cx: -0.5, cy: 0, scale: 3 / w }, { mode: "mandelbrot", maxIter: 200 }, 0, h);
  assert.ok(esc > 0 && esc < w * h);
  // the image is mirror-symmetric about the real axis
  for (let j = 0; j < h; j++) for (let i = 0; i < w; i++) assert.strictEqual(buf[j * w + i], buf[(h - 1 - j) * w + i]);
  // Julia sets have 180° rotational symmetry: J(−z) = J(z)
  const jb = new Float32Array(w * h);
  F.renderRows(jb, w, h, { cx: 0, cy: 0, scale: 3 / w }, { mode: "julia", c: F.JULIA_PRESETS.rabbit.c, maxIter: 200 }, 0, h);
  for (let j = 0; j < h; j++) for (let i = 0; i < w; i++) assert.strictEqual(jb[j * w + i], jb[(h - 1 - j) * w + (w - 1 - i)]);
});

test("area of the Mandelbrot set ≈ 1.5066", () => {
  // best estimates: 1.50659177 ± 8e-8 (Förstemann 2012); a grid count at finite iterations slightly overestimates
  const A = F.mandelbrotArea(800, 1000);
  near(A, 1.5066, 0.02);
});

test("Julia presets lie where they should relative to the Mandelbrot set", () => {
  for (const k of ["rabbit", "dendrite", "sanmarco", "siegel"]) assert.strictEqual(F.mandelbrot(...F.JULIA_PRESETS[k].c, 5000), -1, k); // connected Julia sets
  assert.ok(F.mandelbrot(...F.JULIA_PRESETS.dust.c, 5000) > 0); // outside: the Julia set is dust
});

test("colours are valid and the interior uses the inside colour", () => {
  for (const pal of Object.keys(F.PALETTES)) for (const v of [0, 1.5, 10, 333]) F.colorOf(v, pal, 32).forEach((ch) => assert.ok(ch >= 0 && ch <= 255 && Number.isInteger(ch)));
  assert.deepStrictEqual(F.colorOf(-1, "ocean", 32, [1, 2, 3]), [1, 2, 3]);
});
