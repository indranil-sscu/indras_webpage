/*
 * digitizer-core.js — maths behind the plot digitizer (/tools/digitizer/).
 * No DOM access; unit-tested in Node (tests/digitizer-core.test.js).
 *
 * Calibration: the user clicks two points on the x-axis (X1, X2) and two on
 * the y-axis (Y1, Y2) and types their values. We fit an affine map from
 * pixel (u, v) to data (x, y):
 *     x = a·u + b·v + c,   y = d·u + e·v + f
 * with the conditions
 *     x(X1) = x1, x(X2) = x2, x is constant along the Y1→Y2 line,
 *     y(Y1) = y1, y(Y2) = y2, y is constant along the X1→X2 line.
 * This handles images that are rotated, scaled or sheared (e.g. a scanned
 * page). Logarithmic axes are calibrated in log10 space.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.DigitizerCore = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /** Solve a 3×3 linear system M·s = r (Cramer's rule). */
  function solve3(M, r) {
    const det = (m) =>
      m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
      m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
      m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    const D = det(M);
    if (Math.abs(D) < 1e-12) throw new Error("Calibration points are degenerate (axis points coincide or the axes are parallel).");
    return [0, 1, 2].map((k) => det(M.map((row, i) => row.map((v, j) => (j === k ? r[i] : v)))) / D);
  }

  /**
   * Build a calibration.
   * @param {object} p  pixel positions {X1:[u,v], X2, Y1, Y2}
   * @param {object} v  axis values {x1, x2, y1, y2}
   * @param {object} log {x:boolean, y:boolean}
   */
  function calibrate(p, v, log) {
    log = log || {};
    const tx = (val) => (log.x ? Math.log10(val) : val);
    const ty = (val) => (log.y ? Math.log10(val) : val);
    if (log.x && !(v.x1 > 0 && v.x2 > 0)) throw new Error("Log x-axis values must be positive.");
    if (log.y && !(v.y1 > 0 && v.y2 > 0)) throw new Error("Log y-axis values must be positive.");
    if (v.x1 === v.x2) throw new Error("The two x-axis values must differ.");
    if (v.y1 === v.y2) throw new Error("The two y-axis values must differ.");

    const dY = [p.Y2[0] - p.Y1[0], p.Y2[1] - p.Y1[1]];
    const dX = [p.X2[0] - p.X1[0], p.X2[1] - p.X1[1]];
    const [a, b, c] = solve3(
      [[p.X1[0], p.X1[1], 1], [p.X2[0], p.X2[1], 1], [dY[0], dY[1], 0]],
      [tx(v.x1), tx(v.x2), 0]
    );
    const [d, e, f] = solve3(
      [[p.Y1[0], p.Y1[1], 1], [p.Y2[0], p.Y2[1], 1], [dX[0], dX[1], 0]],
      [ty(v.y1), ty(v.y2), 0]
    );
    return { a, b, c, d, e, f, log: { x: !!log.x, y: !!log.y } };
  }

  /** Pixel (u, v) -> data (x, y). */
  function toData(cal, u, v) {
    const X = cal.a * u + cal.b * v + cal.c;
    const Y = cal.d * u + cal.e * v + cal.f;
    return [cal.log.x ? Math.pow(10, X) : X, cal.log.y ? Math.pow(10, Y) : Y];
  }

  /** Data (x, y) -> pixel (u, v) (inverse of toData). */
  function toPixel(cal, x, y) {
    const X = cal.log.x ? Math.log10(x) : x;
    const Y = cal.log.y ? Math.log10(y) : y;
    const det = cal.a * cal.e - cal.b * cal.d;
    const rx = X - cal.c, ry = Y - cal.f;
    return [(cal.e * rx - cal.b * ry) / det, (-cal.d * rx + cal.a * ry) / det];
  }

  /** Format points as CSV. points: [{series, x, y}] */
  function toCSV(points, opts) {
    opts = opts || {};
    const prec = opts.precision || 6;
    const fmt = (n) => Number(n.toPrecision(prec)).toString();
    let rows = points.slice();
    if (opts.sortByX) rows.sort((p, q) => (p.series < q.series ? -1 : p.series > q.series ? 1 : p.x - q.x));
    const multi = new Set(rows.map((p) => p.series)).size > 1;
    const head = multi ? "series,x,y" : "x,y";
    const quote = (s) => (/[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s);
    return [head, ...rows.map((p) => (multi ? quote(p.series) + "," : "") + fmt(p.x) + "," + fmt(p.y))].join("\n") + "\n";
  }

  /** Index of the point nearest to pixel (u, v) within maxDist, or −1. */
  function nearest(pixelPoints, u, v, maxDist) {
    let best = -1, bd = maxDist * maxDist;
    pixelPoints.forEach((p, i) => {
      const d = (p[0] - u) ** 2 + (p[1] - v) ** 2;
      if (d <= bd) { bd = d; best = i; }
    });
    return best;
  }

  return { calibrate, toData, toPixel, toCSV, nearest };
});
