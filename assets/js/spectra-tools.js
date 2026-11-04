/*
 * spectra-tools.js — spectrum handling for /tools/:
 *   - reading two-column spectra (text / CSV)
 *   - Jacobian conversion between wavelength and wavenumber axes
 * Pure functions, no DOM access (tests/spectra-tools.test.js).
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.SpectraTools = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /** Parse two numeric columns (any of comma, tab, semicolon, whitespace). Lines that
   *  don't start with two numbers (headers, comments) are skipped. */
  function parseTwoColumn(text) {
    const pts = [];
    for (const line of text.split(/\r?\n/)) {
      const f = line.trim().split(/[\s,;]+/);
      if (f.length < 2) continue;
      const x = parseFloat(f[0]), y = parseFloat(f[1]);
      if (isFinite(x) && isFinite(y)) pts.push([x, y]);
    }
    return pts;
  }

  /** Jacobian conversion. direction "nm2cm": I(ν̃) = I(λ)·λ²/10⁷; "cm2nm": I(λ) = I(ν̃)·ν̃²/10⁷.
   *  Conserves the integrated intensity. Output sorted by the new axis. */
  function jacobian(points, direction, normalize) {
    let out = points
      .filter(([x]) => x > 0)
      .map(([x, y]) => [1e7 / x, (y * x * x) / 1e7]); // same algebra both ways
    out.sort((p, q) => p[0] - q[0]);
    if (normalize) {
      const m = Math.max(...out.map((p) => Math.abs(p[1])));
      if (m > 0) out = out.map(([x, y]) => [x, y / m]);
    }
    return out;
  }

  /** Trapezoidal integral of sorted [x, y] points. */
  function trapz(pts) {
    let s = 0;
    for (let i = 1; i < pts.length; i++) s += 0.5 * (pts[i][1] + pts[i - 1][1]) * (pts[i][0] - pts[i - 1][0]);
    return s;
  }

  return { parseTwoColumn, jacobian, trapz };
});
