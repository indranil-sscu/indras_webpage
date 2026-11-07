/*
 * fractal-core.js: the mathematics behind the Mandelbrot and Julia set explorer in
 * "Science for the Curious". Escape-time iteration of z → z² + c, the smooth
 * (continuous) iteration count, orbits, a quick interior test, an area estimate and
 * colour palettes. Pure functions, no DOM; unit-tested in Node (tests/fractal-core.test.js)
 * against analytical results.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.Fractal = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /** Points of the main cardioid and the period-2 disc are in the Mandelbrot set (closed-form test). */
  function inMainBulbs(x, y) {
    const xq = x - 0.25, q = xq * xq + y * y;
    if (q * (q + xq) <= 0.25 * y * y) return true; // main cardioid
    return (x + 1) * (x + 1) + y * y <= 0.0625; // disc of radius 1/4 around −1
  }
  /**
   * Iterate z → z² + c from z0 = (zx, zy) with c = (cx, cy) until |z| > R or maxIter steps.
   * Returns the number of steps n taken when |z_n| first exceeded R (0 if |z0| > R already),
   * or −1 if it never did within maxIter. With `smooth`, returns the continuous count
   *   ν = n − log2( ln|z_n| / ln R ),
   * which is exactly log2(ln R / ln|z0|) for c = 0 and varies continuously between pixels.
   */
  function escape(zx, zy, cx, cy, maxIter, R, smooth) {
    R = R || 2;
    const R2 = R * R;
    let x = zx, y = zy, x2 = x * x, y2 = y * y, n = 0;
    while (x2 + y2 <= R2) {
      if (n >= maxIter) return -1;
      y = 2 * x * y + cy; x = x2 - y2 + cx; x2 = x * x; y2 = y * y; n++;
    }
    if (!smooth) return n;
    return Math.max(0, n - Math.log2((0.5 * Math.log(x2 + y2)) / Math.log(R))); // ≥ 0 so that −1 always means "inside"
  }
  /** Mandelbrot escape count for the point c (z0 = 0), with the interior shortcut. */
  function mandelbrot(cx, cy, maxIter, R, smooth) {
    if (inMainBulbs(cx, cy)) return -1;
    return escape(0, 0, cx, cy, maxIter, R, smooth);
  }
  /** Julia escape count for the starting point z0 with fixed parameter c. */
  const julia = (zx, zy, cx, cy, maxIter, R, smooth) => escape(zx, zy, cx, cy, maxIter, R, smooth);
  /** Orbit z0, z1, … (up to n points, stopping after the first point outside |z| = R). Returns [[x, y], …]. */
  function orbit(zx, zy, cx, cy, n, R) {
    R = R || 2;
    const pts = [[zx, zy]];
    let x = zx, y = zy;
    for (let i = 0; i < n && x * x + y * y <= R * R; i++) { [x, y] = [x * x - y * y + cx, 2 * x * y + cy]; pts.push([x, y]); }
    return pts;
  }
  /**
   * Fill rows j0 … j1−1 of a w × h buffer (Float32Array, row-major, row 0 at the top) with escape values.
   * view = {cx, cy, scale} with scale in units per pixel; o = {mode: "mandelbrot" | "julia", c: [x, y], maxIter, R, smooth}.
   * Returns the number of pixels that escaped.
   */
  function renderRows(buf, w, h, view, o, j0, j1) {
    const julia_ = o.mode === "julia", c = o.c || [0, 0], R = o.R || 2, it = o.maxIter;
    let esc = 0;
    for (let j = j0; j < j1; j++) {
      const y = view.cy + (h / 2 - j - 0.5) * view.scale;
      for (let i = 0; i < w; i++) {
        const x = view.cx + (i + 0.5 - w / 2) * view.scale;
        const v = julia_ ? escape(x, y, c[0], c[1], it, R, o.smooth) : mandelbrot(x, y, it, R, o.smooth);
        buf[j * w + i] = v; if (v >= 0) esc++;
      }
    }
    return esc;
  }
  /** Area of the Mandelbrot set estimated by counting non-escaping grid points (upper half, mirrored). */
  function mandelbrotArea(n, maxIter) {
    const x0 = -2, x1 = 0.5, y1 = 1.25, dx = (x1 - x0) / n, dy = y1 / Math.ceil(n / 2), ny = Math.ceil(n / 2);
    let inside = 0;
    for (let j = 0; j < ny; j++) for (let i = 0; i < n; i++) if (mandelbrot(x0 + (i + 0.5) * dx, (j + 0.5) * dy, maxIter) < 0) inside++;
    return 2 * inside * dx * dy;
  }

  /* ---------------- colour ---------------- */
  /** Cosine gradients (after I. Quilez): colour(t) = a + b cos 2π(c t + d), channels in 0 … 1. */
  const PALETTES = {
    ocean: { name: "Ocean", a: [0.5, 0.5, 0.5], b: [0.5, 0.5, 0.5], c: [1, 1, 1], d: [0.0, 0.1, 0.2] },
    fire: { name: "Fire", a: [0.5, 0.3, 0.1], b: [0.5, 0.4, 0.2], c: [1, 1, 1], d: [0.0, 0.08, 0.18] },
    rainbow: { name: "Rainbow", a: [0.5, 0.5, 0.5], b: [0.5, 0.5, 0.5], c: [1, 1, 1], d: [0.0, 0.33, 0.67] },
    grey: { name: "Greyscale", a: [0.5, 0.5, 0.5], b: [0.5, 0.5, 0.5], c: [0.5, 0.5, 0.5], d: [0.5, 0.5, 0.5] },
  };
  /** Colour of an escape value v (−1 = inside) as [r, g, b] in 0 … 255. `cycle` = iterations per colour cycle. */
  function colorOf(v, palette, cycle, inside) {
    if (v < 0) return inside || [0, 0, 0];
    const P = PALETTES[palette] || PALETTES.ocean, t = Math.sqrt(Math.max(0, v)) / Math.sqrt(cycle || 32);
    return [0, 1, 2].map((k) => Math.round(255 * Math.min(1, Math.max(0, P.a[k] + P.b[k] * Math.cos(2 * Math.PI * (P.c[k] * t + P.d[k]))))));
  }
  /** Well-known Julia parameters. */
  const JULIA_PRESETS = {
    rabbit: { name: "Douady rabbit", c: [-0.123, 0.745] },
    dendrite: { name: "Dendrite (c = i)", c: [0, 1] },
    sanmarco: { name: "San Marco (c = −3/4)", c: [-0.75, 0] },
    siegel: { name: "Siegel disc", c: [-0.390541, -0.586788] },
    spiral: { name: "Spirals (c = −0.8 + 0.156i)", c: [-0.8, 0.156] },
    dust: { name: "Cantor dust (c = 0.3 + 0.6i)", c: [0.3, 0.6] },
  };

  return { inMainBulbs, escape, mandelbrot, julia, orbit, renderRows, mandelbrotArea, PALETTES, colorOf, JULIA_PRESETS };
});
