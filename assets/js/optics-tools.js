/*
 * optics-tools.js — physics for the optics calculators on /tools/.
 * Pure functions, no DOM access (unit-tested in Node: tests/optics-tools.test.js).
 *
 * Wavelengths inside the dispersion formulas are in micrometres (Sellmeier
 * convention); public functions take nanometres unless stated otherwise.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.OpticsTools = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  const c = 299792458; // m/s

  /* ------------------------------------------------------------------ *
   * Refractive-index models
   * ------------------------------------------------------------------ */
  // Three-term Sellmeier: n² − 1 = Σ B_i λ² / (λ² − C_i), λ in µm, C_i in µm²
  const sellmeier3 = (B, C) => (l) => {
    const l2 = l * l;
    let s = 1;
    for (let i = 0; i < B.length; i++) s += (B[i] * l2) / (l2 - C[i]);
    return Math.sqrt(s);
  };
  const sq = (x) => x * x;

  const MATERIALS = {
    fused_silica: {
      name: "Fused silica",
      ref: "I. H. Malitson, J. Opt. Soc. Am. 55, 1205 (1965)",
      range: [0.21, 6.7],
      n: sellmeier3([0.6961663, 0.4079426, 0.8974794], [sq(0.0684043), sq(0.1162414), sq(9.896161)]),
    },
    bk7: {
      name: "N-BK7 (Schott)",
      ref: "Schott optical glass data sheet",
      range: [0.3, 2.5],
      n: sellmeier3([1.03961212, 0.231792344, 1.01046945], [0.00600069867, 0.0200179144, 103.560653]),
    },
    sf10: {
      name: "SF10 (Schott)",
      ref: "Schott optical glass data sheet",
      range: [0.38, 2.5],
      n: sellmeier3([1.62153902, 0.256287842, 1.64447552], [0.0122241457, 0.0595736775, 147.468793]),
    },
    caf2: {
      name: "CaF₂",
      ref: "I. H. Malitson, Appl. Opt. 2, 1103 (1963)",
      range: [0.23, 9.7],
      n: sellmeier3([0.5675888, 0.4710914, 3.8484723], [sq(0.050263605), sq(0.1003909), sq(34.649040)]),
    },
    sapphire_o: {
      name: "Sapphire (ordinary ray)",
      ref: "I. H. Malitson & M. J. Dodge, J. Opt. Soc. Am. 62, 1405 (1972)",
      range: [0.2, 5.0],
      n: sellmeier3([1.4313493, 0.65054713, 5.3414021], [sq(0.0726631), sq(0.1193242), sq(18.028251)]),
    },
    bbo_o: {
      name: "β-BBO (ordinary ray)",
      ref: "D. Eimerl et al., J. Appl. Phys. 62, 1968 (1987)",
      range: [0.22, 1.06],
      n: (l) => Math.sqrt(2.7405 + 0.0184 / (l * l - 0.0179) - 0.0155 * l * l),
    },
    bbo_e: {
      name: "β-BBO (extraordinary ray)",
      ref: "D. Eimerl et al., J. Appl. Phys. 62, 1968 (1987)",
      range: [0.22, 1.06],
      n: (l) => Math.sqrt(2.3730 + 0.0128 / (l * l - 0.0156) - 0.0044 * l * l),
    },
  };

  /* ------------------------------------------------------------------ *
   * 1. Dispersion of a material
   *    n_g = n − λ dn/dλ
   *    GVD = λ³/(2πc²) d²n/dλ²
   *    TOD = −λ⁴/(4π²c³) (3 d²n/dλ² + λ d³n/dλ³)
   * ------------------------------------------------------------------ */
  function derivatives(nfun, l_um) {
    // central differences; step sizes chosen to balance truncation and round-off
    const h1 = 1e-4, h2 = 1e-3, h3 = 4e-3;
    const d1 = (nfun(l_um + h1) - nfun(l_um - h1)) / (2 * h1);
    const d2 = (nfun(l_um + h2) - 2 * nfun(l_um) + nfun(l_um - h2)) / (h2 * h2);
    const d3 = (nfun(l_um + 2 * h3) - 2 * nfun(l_um + h3) + 2 * nfun(l_um - h3) - nfun(l_um - 2 * h3)) / (2 * h3 * h3 * h3);
    return { d1, d2, d3 }; // per µm, µm⁻², µm⁻³
  }

  /** Optical properties of `materialKey` at `lambda_nm`. */
  function dispersion(materialKey, lambda_nm) {
    const m = MATERIALS[materialKey];
    if (!m) throw new Error("Unknown material: " + materialKey);
    const l = lambda_nm / 1000; // µm
    const n = m.n(l);
    const { d1, d2, d3 } = derivatives(m.n, l);
    const lm = l * 1e-6; // m
    const d2m = d2 * 1e12, d3m = d3 * 1e18; // per m², per m³
    const gvd = (lm ** 3 / (2 * Math.PI * c * c)) * d2m; // s²/m
    const tod = -(lm ** 4 / (4 * Math.PI * Math.PI * c ** 3)) * (3 * d2m + lm * d3m); // s³/m
    return {
      n,
      dn_dl_um: d1,
      ng: n - l * d1,
      gvd_fs2_mm: gvd * 1e30 / 1e3, // fs²/mm
      tod_fs3_mm: tod * 1e45 / 1e3, // fs³/mm
      inRange: l >= m.range[0] && l <= m.range[1],
    };
  }

  /* ------------------------------------------------------------------ *
   * 2. Pulse compressors (double pass)
   * ------------------------------------------------------------------ */
  /** Treacy grating pair. G = perpendicular grating separation (mm),
   *  grooves per mm, angle of incidence (deg), first order. */
  function gratingCompressor(lambda_nm, grooves_per_mm, incidence_deg, separation_mm) {
    const l = lambda_nm * 1e-9;
    const d = 1e-3 / grooves_per_mm;
    const sinD = l / d - Math.sin((incidence_deg * Math.PI) / 180);
    if (Math.abs(sinD) >= 1) throw new Error("No first diffraction order at this angle");
    const cosD = Math.sqrt(1 - sinD * sinD);
    const G = separation_mm * 1e-3;
    const gdd = -(l ** 3 * G) / (Math.PI * c * c * d * d * cosD ** 3); // s², double pass
    // TOD/GDD ratio for a grating pair: −(3λ / 2πc)(1 + λ sinD / (d cos²D))
    const tod = gdd * (-(3 * l) / (2 * Math.PI * c)) * (1 + (l * sinD) / (d * cosD * cosD));
    return {
      diffraction_deg: (Math.asin(sinD) * 180) / Math.PI,
      gdd_fs2: gdd * 1e30,
      tod_fs3: tod * 1e45,
    };
  }

  /** Brewster-angle prism pair, double pass (leading-order Fork et al. 1984):
   *  GDD ≈ −(4 L λ³ / π c²)(dn/dλ)² + GVD_material × glass path. */
  function prismCompressor(materialKey, lambda_nm, apexSeparation_mm, glassPath_mm) {
    const disp = dispersion(materialKey, lambda_nm);
    const l = lambda_nm * 1e-9;
    const dndl = disp.dn_dl_um * 1e6; // per m
    const L = apexSeparation_mm * 1e-3;
    const geom = (-(4 * L * l ** 3) / (Math.PI * c * c)) * dndl * dndl * 1e30; // fs²
    const material = disp.gvd_fs2_mm * (glassPath_mm || 0);
    const n = disp.n;
    return {
      gdd_fs2: geom + material,
      geometric_fs2: geom,
      material_fs2: material,
      brewster_deg: (Math.atan(n) * 180) / Math.PI,
      apex_deg: 180 - 2 * ((Math.atan(n) * 180) / Math.PI), // minimum-deviation prism cut for Brewster
    };
  }

  /** Separation needed for a target (negative) GDD. */
  function separationForGdd(gddTarget_fs2, gddPerMm_fs2, offset_fs2) {
    return (gddTarget_fs2 - (offset_fs2 || 0)) / gddPerMm_fs2;
  }

  /* ------------------------------------------------------------------ *
   * 3. Birefringent phase matching in negative uniaxial crystals
   * ------------------------------------------------------------------ */
  const CRYSTALS = {
    bbo: { name: "β-BBO", no: MATERIALS.bbo_o.n, ne: MATERIALS.bbo_e.n, ref: MATERIALS.bbo_o.ref, range: [0.22, 1.06] },
    kdp: {
      name: "KDP",
      ref: "F. Zernike, J. Opt. Soc. Am. 54, 1215 (1964)",
      range: [0.21, 1.5],
      no: (l) => Math.sqrt(2.259276 + 0.01008956 / (l * l - 0.012942625) + (13.00522 * l * l) / (l * l - 400)),
      ne: (l) => Math.sqrt(2.132668 + 0.008637494 / (l * l - 0.012281043) + (3.2279924 * l * l) / (l * l - 400)),
    },
  };

  const neTheta = (no, ne, th) => 1 / Math.sqrt(sq(Math.cos(th) / no) + sq(Math.sin(th) / ne));

  /** Walk-off angle of an e-wave at polar angle θ (radians). */
  function walkoff(no, ne, th) {
    return Math.atan((sq(neTheta(no, ne, th)) / 2) * (1 / sq(ne) - 1 / sq(no)) * Math.sin(2 * th));
  }

  /** Sum-frequency phase matching λ1 + λ2 → λ3 (nm). type "I" (o+o→e) or "II" (e+o→e). */
  function phaseMatch(crystalKey, l1_nm, l2_nm, type) {
    const X = CRYSTALS[crystalKey];
    if (!X) throw new Error("Unknown crystal: " + crystalKey);
    const l1 = l1_nm / 1000, l2 = l2_nm / 1000, l3 = 1 / (1 / l1 + 1 / l2);
    const k = (n, l) => n / l; // ∝ wavevector
    const no3 = X.no(l3), ne3 = X.ne(l3);
    const mismatch = (th) => {
      const k1 = type === "II" ? k(neTheta(X.no(l1), X.ne(l1), th), l1) : k(X.no(l1), l1);
      const k2 = k(X.no(l2), l2);
      return k(neTheta(no3, ne3, th), l3) - k1 - k2;
    };
    // Δk decreases monotonically from θ = 0 to 90°: bisection
    let a = 0, b = Math.PI / 2;
    if (mismatch(a) < 0 || mismatch(b) > 0) {
      return { possible: false, lambda3_nm: l3 * 1000 };
    }
    for (let i = 0; i < 100; i++) {
      const m = (a + b) / 2;
      if (mismatch(m) > 0) a = m; else b = m;
    }
    const th = (a + b) / 2;
    return {
      possible: true,
      lambda3_nm: l3 * 1000,
      theta_deg: (th * 180) / Math.PI,
      walkoff_mrad: walkoff(no3, ne3, th) * 1e3,
      inRange: [l1, l2, l3].every((l) => l >= X.range[0] && l <= X.range[1]),
    };
  }

  /* ------------------------------------------------------------------ *
   * 4. Gaussian beams (w = 1/e² intensity radius)
   * ------------------------------------------------------------------ */
  function gaussianBeam(lambda_nm, w0_um, M2, z_mm) {
    const l = lambda_nm * 1e-9, w0 = w0_um * 1e-6, m2 = M2 || 1;
    const zR = (Math.PI * w0 * w0) / (m2 * l);
    const z = (z_mm || 0) * 1e-3;
    return {
      rayleigh_mm: zR * 1e3,
      divergence_mrad: ((m2 * l) / (Math.PI * w0)) * 1e3, // half-angle
      w_at_z_um: w0 * Math.sqrt(1 + sq(z / zR)) * 1e6,
      R_at_z_mm: z === 0 ? Infinity : (z * (1 + sq(zR / z))) * 1e3,
      fwhm_um: w0 * Math.sqrt(2 * Math.LN2) * 1e6,
    };
  }

  /** Focusing a collimated beam of 1/e² diameter D (mm) with focal length f (mm). */
  function focusSpot(lambda_nm, D_mm, f_mm, M2) {
    const l = lambda_nm * 1e-9, w_in = (D_mm * 1e-3) / 2, f = f_mm * 1e-3, m2 = M2 || 1;
    const w0 = (m2 * l * f) / (Math.PI * w_in);
    const zR = (Math.PI * w0 * w0) / (m2 * l);
    return { w0_um: w0 * 1e6, spot_diameter_um: 2 * w0 * 1e6, rayleigh_mm: zR * 1e3, depth_of_focus_mm: 2 * zR * 1e3 };
  }

  /* ------------------------------------------------------------------ *
   * 5. Intensity autocorrelation deconvolution factors (τ_AC / τ_pulse)
   * ------------------------------------------------------------------ */
  const AC_FACTOR = { gaussian: Math.SQRT2, sech2: 1.5427, lorentzian: 2 };
  const pulseFromAutocorrelation = (tauAC_fs, shape) => tauAC_fs / AC_FACTOR[shape];

  /* ------------------------------------------------------------------ *
   * 6. Sampling / FFT axis design for time-domain scans
   * ------------------------------------------------------------------ */
  const CM1_PER_THZ = 1e12 / (c * 100); // 33.356 cm⁻¹ per THz

  /** dt in fs, N points, optional rotating-frame (carrier) frequency in cm⁻¹. */
  function fftAxis(dt_fs, N, frame_cm1) {
    const df_THz = 1e3 / (N * dt_fs); // 1/(N dt), fs⁻¹ → THz
    const fNyq_THz = 1e3 / (2 * dt_fs);
    const kmin = -Math.floor(N / 2), kmax = Math.ceil(N / 2) - 1;
    const f0 = frame_cm1 || 0;
    return {
      window_fs: N * dt_fs,
      resolution_THz: df_THz,
      resolution_cm1: df_THz * CM1_PER_THZ,
      nyquist_THz: fNyq_THz,
      nyquist_cm1: fNyq_THz * CM1_PER_THZ,
      axisMin_cm1: f0 + kmin * df_THz * CM1_PER_THZ,
      axisMax_cm1: f0 + kmax * df_THz * CM1_PER_THZ,
      zeroFrequencyIndex: Math.floor(N / 2), // index of f = 0 after fftshift (0-based)
    };
  }

  /** Apparent (aliased) frequency of a signal at f_cm1 sampled every dt_fs,
   *  measured in a frame rotating at frame_cm1. Returned in absolute cm⁻¹. */
  function aliasedFrequency(f_cm1, dt_fs, frame_cm1) {
    const fs_cm1 = (1e3 / dt_fs) * CM1_PER_THZ; // sampling frequency
    const rel = f_cm1 - (frame_cm1 || 0);
    const folded = rel - fs_cm1 * Math.round(rel / fs_cm1); // into [−fs/2, fs/2]
    return { apparent_cm1: (frame_cm1 || 0) + folded, aliased: Math.abs(rel) > fs_cm1 / 2, sampling_cm1: fs_cm1 };
  }

  return {
    MATERIALS,
    dispersion,
    gratingCompressor,
    prismCompressor,
    separationForGdd,
    CRYSTALS,
    phaseMatch,
    gaussianBeam,
    focusSpot,
    AC_FACTOR,
    pulseFromAutocorrelation,
    CM1_PER_THZ,
    fftAxis,
    aliasedFrequency,
  };
});
