/*
 * sci-tools.js — physics behind the conversion calculators on /tools/.
 * Pure functions, no DOM access, so they can be unit-tested in Node.
 *
 * Constants: exact SI (2019 redefinition) unless noted; Hartree energy from CODATA 2018.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.SciTools = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  const C = {
    h: 6.62607015e-34, // Planck constant, J s (exact)
    c: 299792458, // speed of light, m/s (exact)
    e: 1.602176634e-19, // elementary charge, C (exact)
    kB: 1.380649e-23, // Boltzmann constant, J/K (exact)
    NA: 6.02214076e23, // Avogadro constant, 1/mol (exact)
    Eh: 4.3597447222071e-18, // Hartree energy, J (CODATA 2018)
    cal: 4.184, // thermochemical calorie, J (exact)
  };
  C.hbar = C.h / (2 * Math.PI);
  C.Ry = C.Eh / 2;

  /* ------------------------------------------------------------------ *
   * 1. Energy / frequency / wavelength units
   *    Each unit converts a value <-> energy per particle in joules.
   * ------------------------------------------------------------------ */
  const ENERGY_UNITS = {
    nm: { label: "Wavelength", unit: "nm", toJ: (v) => (C.h * C.c) / (v * 1e-9), fromJ: (J) => (C.h * C.c) / J / 1e-9 },
    um: { label: "Wavelength", unit: "µm", toJ: (v) => (C.h * C.c) / (v * 1e-6), fromJ: (J) => (C.h * C.c) / J / 1e-6 },
    cm1: { label: "Wavenumber", unit: "cm⁻¹", toJ: (v) => C.h * C.c * v * 100, fromJ: (J) => J / (C.h * C.c * 100) },
    THz: { label: "Frequency", unit: "THz", toJ: (v) => C.h * v * 1e12, fromJ: (J) => J / C.h / 1e12 },
    radfs: { label: "Angular frequency", unit: "rad/fs", toJ: (v) => C.hbar * v * 1e15, fromJ: (J) => J / C.hbar / 1e15 },
    fs: { label: "Period", unit: "fs", toJ: (v) => C.h / (v * 1e-15), fromJ: (J) => C.h / J / 1e-15 },
    eV: { label: "Energy", unit: "eV", toJ: (v) => v * C.e, fromJ: (J) => J / C.e },
    meV: { label: "Energy", unit: "meV", toJ: (v) => v * 1e-3 * C.e, fromJ: (J) => J / C.e / 1e-3 },
    Eh: { label: "Energy", unit: "Hartree", toJ: (v) => v * C.Eh, fromJ: (J) => J / C.Eh },
    Ry: { label: "Energy", unit: "Rydberg", toJ: (v) => v * C.Ry, fromJ: (J) => J / C.Ry },
    kJmol: { label: "Molar energy", unit: "kJ/mol", toJ: (v) => (v * 1e3) / C.NA, fromJ: (J) => (J * C.NA) / 1e3 },
    kcalmol: { label: "Molar energy", unit: "kcal/mol", toJ: (v) => (v * 1e3 * C.cal) / C.NA, fromJ: (J) => (J * C.NA) / (1e3 * C.cal) },
    K: { label: "Temperature (E = k<sub>B</sub>T)", unit: "K", toJ: (v) => C.kB * v, fromJ: (J) => J / C.kB },
    J: { label: "Energy", unit: "J", toJ: (v) => v, fromJ: (J) => J },
  };

  /** Convert `value` in unit `from` to every unit. Returns {unitKey: number}. */
  function convertEnergy(value, from) {
    const u = ENERGY_UNITS[from];
    if (!u) throw new Error("Unknown unit: " + from);
    if (!(value > 0)) throw new Error("Value must be positive");
    const J = u.toJ(value);
    const out = {};
    for (const k in ENERGY_UNITS) out[k] = ENERGY_UNITS[k].fromJ(J);
    return out;
  }

  /* ------------------------------------------------------------------ *
   * 2. Spectral bandwidth: Δλ at λ0  <->  Δν, Δν̃, ΔE
   *    Exact for a band [λ0 − Δλ/2, λ0 + Δλ/2] (no small-bandwidth approximation).
   * ------------------------------------------------------------------ */
  function bandwidthFromNm(lambda0_nm, dlambda_nm) {
    const l1 = (lambda0_nm - dlambda_nm / 2) * 1e-9;
    const l2 = (lambda0_nm + dlambda_nm / 2) * 1e-9;
    if (!(l1 > 0)) throw new Error("Δλ must be smaller than 2 λ0");
    const dnu_Hz = C.c / l1 - C.c / l2;
    return {
      dnu_THz: dnu_Hz / 1e12,
      dk_cm1: dnu_Hz / C.c / 100,
      dE_meV: (C.h * dnu_Hz) / C.e * 1e3,
    };
  }

  /** Exact inverse of bandwidthFromNm: frequency width (THz) at λ0 -> Δλ (nm).
   *  Solves Δν = c·Δλ / (λ0² − Δλ²/4) for Δλ. */
  function bandwidthToNm(lambda0_nm, dnu_THz) {
    const l0 = lambda0_nm * 1e-9;
    const dnu = dnu_THz * 1e12;
    if (!(dnu > 0)) throw new Error("Bandwidth must be positive");
    // (√(c² + b²) − c) / (Δν/2) rewritten to avoid cancellation for narrow bands
    return ((2 * dnu * l0 * l0) / (Math.sqrt(C.c * C.c + dnu * dnu * l0 * l0) + C.c)) * 1e9;
  }

  /* ------------------------------------------------------------------ *
   * 3. Transform-limited pulses (FWHM in intensity)
   * ------------------------------------------------------------------ */
  const TBP = {
    gaussian: (2 * Math.LN2) / Math.PI, // 0.4413
    sech2: 0.314833, // (2 acosh(√2) / π)² = 0.3148
  };

  /** Shortest pulse (fs) supported by a spectrum of FWHM Δλ (nm) at λ0 (nm). */
  function transformLimitedDuration(lambda0_nm, dlambda_nm, shape) {
    const k = TBP[shape];
    if (!k) throw new Error("Unknown pulse shape: " + shape);
    const dnu = bandwidthFromNm(lambda0_nm, dlambda_nm).dnu_THz * 1e12;
    return (k / dnu) * 1e15;
  }

  /** Minimum spectral FWHM (nm) needed for a pulse of duration τ (fs) at λ0. */
  function transformLimitedBandwidth(lambda0_nm, tau_fs, shape) {
    const k = TBP[shape];
    if (!k) throw new Error("Unknown pulse shape: " + shape);
    const dnu_THz = (k / (tau_fs * 1e-15)) / 1e12;
    return bandwidthToNm(lambda0_nm, dnu_THz);
  }

  /* ------------------------------------------------------------------ *
   * 4. Optical delay line: stage displacement <-> time delay
   *    t = passes · n_g · L / c   (passes = 2 for a retro-reflector)
   * ------------------------------------------------------------------ */
  function delayFromDistance(distance_mm, passes, groupIndex) {
    return (passes * (groupIndex || 1) * distance_mm * 1e-3) / C.c * 1e15; // fs
  }
  function distanceFromDelay(delay_fs, passes, groupIndex) {
    return (delay_fs * 1e-15 * C.c) / (passes * (groupIndex || 1)) * 1e3; // mm
  }

  /* ------------------------------------------------------------------ *
   * 5. Pulse energy, peak power, fluence and intensity
   *    Gaussian beam (1/e² diameter D), temporal shape Gaussian or sech².
   *    Peak power  P0 = k · E / τ,  k = 0.9394 (Gaussian), 0.8814 (sech²)
   *    Peak fluence F0 = 2E / (π w²), peak intensity I0 = 2 P0 / (π w²), w = D/2
   * ------------------------------------------------------------------ */
  const PEAK_FACTOR = {
    gaussian: 2 * Math.sqrt(Math.LN2 / Math.PI), // 0.9394
    sech2: Math.acosh(Math.SQRT2), // ln(1+√2) = 0.8814
  };

  function pulseParameters({ avgPower_W, repRate_Hz, duration_fs, diameter_um, wavelength_nm, shape }) {
    const E = avgPower_W / repRate_Hz; // J
    const out = { energy_J: E };
    if (duration_fs > 0) {
      out.peakPower_W = (PEAK_FACTOR[shape || "gaussian"] * E) / (duration_fs * 1e-15);
    }
    if (diameter_um > 0) {
      const w = (diameter_um * 1e-6) / 2; // 1/e² radius, m
      const area = (Math.PI * w * w) / 2; // effective area for peak values, m²
      out.peakFluence_Jcm2 = E / area / 1e4;
      if (out.peakPower_W) out.peakIntensity_Wcm2 = out.peakPower_W / area / 1e4;
    }
    if (wavelength_nm > 0) {
      out.photonsPerPulse = E / ((C.h * C.c) / (wavelength_nm * 1e-9));
    }
    return out;
  }

  /* ------------------------------------------------------------------ *
   * 6. Dispersion: broadening of an unchirped Gaussian pulse by GDD
   *    τ_out = τ_in · sqrt(1 + (4 ln2 · GDD / τ_in²)²)
   * ------------------------------------------------------------------ */
  function gddBroadening(tau_in_fs, gdd_fs2) {
    const x = (4 * Math.LN2 * gdd_fs2) / (tau_in_fs * tau_in_fs);
    return tau_in_fs * Math.sqrt(1 + x * x);
  }

  return {
    C,
    ENERGY_UNITS,
    convertEnergy,
    bandwidthFromNm,
    bandwidthToNm,
    TBP,
    transformLimitedDuration,
    transformLimitedBandwidth,
    delayFromDistance,
    distanceFromDelay,
    PEAK_FACTOR,
    pulseParameters,
    gddBroadening,
  };
});
