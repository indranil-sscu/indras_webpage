/*
 * spectra-tools.js — spectrum handling for /tools/:
 *   - reading two-column spectra (text / CSV)
 *   - Jacobian conversion between wavelength and wavenumber axes
 *   - first-order (field) autocorrelation and pulse simulation from a spectrum
 * Pure functions, no DOM access (tests/spectra-tools.test.js).
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.SpectraTools = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  const c = 299792458; // m/s
  const C_NM_PER_FS = c * 1e-6; // 299.792458 nm/fs

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

  /** FWHM of a sampled, single-peaked curve (linear interpolation at half maximum). */
  function fwhm(xs, ys) {
    let im = 0;
    for (let i = 1; i < ys.length; i++) if (ys[i] > ys[im]) im = i;
    const half = ys[im] / 2;
    let i1 = im, i2 = im;
    while (i1 > 0 && ys[i1] > half) i1--;
    while (i2 < ys.length - 1 && ys[i2] > half) i2++;
    const interp = (a, b) => xs[a] + ((half - ys[a]) * (xs[b] - xs[a])) / (ys[b] - ys[a]);
    return interp(i2, i2 - 1) - interp(i1, i1 + 1);
  }

  /**
   * Spectrum on a uniform frequency grid (PHz = fs⁻¹).
   * Gaussian: intensity spectrum Gaussian in frequency, centred at c/λ0, with the
   * frequency FWHM spanned by λ0 ± Δλ/2.
   */
  function gaussianSpectrum(lambda0_nm, dlambda_nm, n) {
    const nu0 = C_NM_PER_FS / lambda0_nm;
    const dnu = C_NM_PER_FS / (lambda0_nm - dlambda_nm / 2) - C_NM_PER_FS / (lambda0_nm + dlambda_nm / 2);
    const sig = dnu / (2 * Math.sqrt(2 * Math.LN2));
    const N = n || 512, span = 16 * sig; // ±8σ in intensity keeps the field amplitude √S untruncated
    const nu = [], S = [];
    for (let i = 0; i < N; i++) {
      const v = nu0 - span / 2 + (span * i) / (N - 1);
      nu.push(v); S.push(Math.exp(-((v - nu0) ** 2) / (2 * sig * sig)));
    }
    return { nu, S, nu0 };
  }

  /** Experimental spectrum I(λ) [nm] → uniform frequency grid with the Jacobian. */
  function spectrumFromWavelength(points, n) {
    const conv = points.filter(([l]) => l > 0).map(([l, I]) => [C_NM_PER_FS / l, Math.max(0, I) * (l * l) / C_NM_PER_FS]).sort((a, b) => a[0] - b[0]);
    if (conv.length < 3) throw new Error("Need at least three data points.");
    const N = n || 512, lo = conv[0][0], hi = conv[conv.length - 1][0];
    const nu = [], S = [];
    let j = 0;
    for (let i = 0; i < N; i++) {
      const v = lo + ((hi - lo) * i) / (N - 1);
      while (j < conv.length - 2 && conv[j + 1][0] < v) j++;
      const [x0, y0] = conv[j], [x1, y1] = conv[j + 1];
      nu.push(v); S.push(y0 + ((y1 - y0) * (v - x0)) / (x1 - x0 || 1));
    }
    let w = 0, ws = 0;
    nu.forEach((v, i) => { w += S[i]; ws += S[i] * v; });
    return { nu, S, nu0: ws / w };
  }

  /**
   * Simulate from a spectrum {nu [fs⁻¹], S, nu0}:
   *  - field autocorrelation G1(τ) = ∫ S(ν) e^{−i2πντ} dν (interferogram = Re G1, envelope = |G1|)
   *  - pulse field E(t) = ∫ √S(ν) e^{iφ(ν)} e^{−i2πνt} dν with φ = ½·GDD·(ω − ω0)²
   *  - spectrum recovered from the interferogram (Wiener–Khinchin)
   * Delays/times in fs, GDD in fs².
   */
  function simulate(spec, gdd_fs2, opts) {
    opts = opts || {};
    const { nu, S, nu0 } = spec;
    const dnu = nu[1] - nu[0];
    const amp = S.map((s) => Math.sqrt(Math.max(0, s)));
    const phase = nu.map((v) => 0.5 * gdd_fs2 * (2 * Math.PI * (v - nu0)) ** 2);

    // envelope durations: TL and chirped pulse intensity FWHM (computed on a coarse grid)
    const pulse = (gdd, T, M) => {
      const ts = [], I = [], re = [], env = [];
      for (let k = 0; k < M; k++) {
        const t = -T + (2 * T * k) / (M - 1);
        let sr = 0, si = 0;
        for (let i = 0; i < nu.length; i++) {
          const ph = (gdd ? phase[i] : 0) - 2 * Math.PI * (nu[i] - nu0) * t; // in the frame of the carrier
          sr += amp[i] * Math.cos(ph); si += amp[i] * Math.sin(ph);
        }
        const e = Math.hypot(sr, si);
        ts.push(t); I.push(e * e); env.push(e);
        re.push((sr * Math.cos(2 * Math.PI * nu0 * t) + si * Math.sin(2 * Math.PI * nu0 * t)));
      }
      return { ts, I, env, re };
    };
    // estimate the TL duration from the spectral FWHM to size the time grids
    const specFwhm = fwhm(nu, S); // fs⁻¹
    const tlEst = (2 * Math.LN2) / Math.PI / specFwhm; // fs
    const tl = pulse(false, 3 * tlEst, 1201);
    const tlFwhm = fwhm(tl.ts, tl.I);
    const stretchEst = tlFwhm * Math.sqrt(1 + ((4 * Math.LN2 * gdd_fs2) / (tlFwhm * tlFwhm)) ** 2);
    const ch = pulse(true, 3 * Math.max(stretchEst, tlFwhm), 1601);
    const chFwhm = fwhm(ch.ts, ch.I);

    // interferogram: delay grid fine enough to sample the carrier
    const T = opts.delayRange_fs || Math.max(4 * tlFwhm, 30);
    const dt = Math.min(0.2, 1 / (8 * nu[nu.length - 1]));
    const M = Math.min(8001, Math.floor((2 * T) / dt) + 1);
    const tau = [], G = [], Genv = [];
    for (let k = 0; k < M; k++) {
      const t = -T + (2 * T * k) / (M - 1);
      let sr = 0, si = 0;
      for (let i = 0; i < nu.length; i++) { const ph = -2 * Math.PI * nu[i] * t; sr += S[i] * Math.cos(ph); si += S[i] * Math.sin(ph); }
      tau.push(t); G.push(sr); Genv.push(Math.hypot(sr, si));
    }
    const g0 = Math.max(...Genv);
    const Gn = G.map((g) => g / g0), Genvn = Genv.map((g) => g / g0);

    // recover the spectrum from the interferogram (cosine transform), back on the ν grid
    const dtau = tau[1] - tau[0];
    const rec = nu.map((v) => { let s = 0; for (let k = 0; k < M; k++) s += Gn[k] * Math.cos(2 * Math.PI * v * tau[k]); return s * dtau; });
    const rmax = Math.max(...rec), smax = Math.max(...S);

    return {
      tlFwhm_fs: tlFwhm,
      stretchedFwhm_fs: chFwhm,
      acEnvelopeFwhm_fs: fwhm(tau, Genvn),
      tau, interferogram: Gn, envelope: Genvn,
      pulse: { t: ch.ts, field: ch.re.map((x) => x / Math.max(...ch.env)), envelope: ch.env.map((x) => x / Math.max(...ch.env)) },
      spectrum: { lambda_nm: nu.map((v) => C_NM_PER_FS / v), original: S.map((s) => s / smax), recovered: rec.map((r) => r / rmax) },
    };
  }

  return { parseTwoColumn, jacobian, trapz, fwhm, gaussianSpectrum, spectrumFromWavelength, simulate, C_NM_PER_FS };
});
