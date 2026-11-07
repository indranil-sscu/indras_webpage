/*
 * twodes-core.js: the numerical engine of the 2D electronic spectroscopy (2DES) simulator.
 * Pure functions, no DOM; runs in the page, in a Web Worker (twodes-worker.js) and in Node,
 * where tests/twodes-core.test.js checks it against analytical results.
 *
 * CONVENTIONS (impulsive limit, third-order response R(t3, T, t1); times in fs, angular
 * frequencies in rad/fs):
 *   Built-in models use the E ∝ exp(−iωt) convention, in which a rephasing coherence evolves as
 *   exp(+iω_eg t1) during t1 and every emitting coherence as exp(−iω t3) during t3:
 *     R_R  ∝ exp(+iω_a t1) exp(−iω_b t3),   R_NR ∝ exp(−iω_a t1) exp(−iω_b t3).
 *   Spectra are one-sided double Fourier integrals  S(ω1, ω3) = ∫0^∞ dt1 ∫0^∞ dt3 R K,
 *   with kernel K = exp(i s1 ω1 t1) exp(i s3 ω3 t3). The "standard" choice is
 *     rephasing (s1, s3) = (−1, +1),  non-rephasing (s1, s3) = (+1, +1),
 *   which puts both pathways in the (+ω1, +ω3) quadrant. The absorptive spectrum is
 *   Re[S_R + S_NR] (Khalil, Demirdöven & Tokmakoff 2003); with these conventions GSB and SE
 *   are positive and ESA negative. The emitted field's factor of i and the local-oscillator
 *   phase are absorbed into this sign convention.
 *
 * ROTATING FRAME: data may be stored with a carrier ωc removed, defined (for any kernel) as
 *   R_rot(t) = R_lab(t) · K_std(ωc, t),   i.e.  R_rot = R_lab e^{−iωc t1 + iωc t3} (rephasing),
 *   R_rot = R_lab e^{+iωc t1 + iωc t3} (non-rephasing).
 * Then S(ω) = ∫ R_rot exp(i s Ω t) with Ω = ω − σ ωc, σ = s / s_std, so the sampled data only
 * need to resolve frequencies near ωc. Flipping a kernel sign mirrors that frequency axis.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.TwoDES = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";
  const TAU = 2 * Math.PI;
  const HBAR = 0.6582119569; // eV fs
  const C_CM_FS = 2.99792458e-5; // speed of light, cm/fs
  const KB = 8.617333262e-5; // eV/K

  /* ======================= units ======================= */
  const UNITS = {
    cm: { label: "cm⁻¹", to: (w) => w / (TAU * C_CM_FS), from: (v) => v * TAU * C_CM_FS, digits: 0 },
    ev: { label: "eV", to: (w) => w * HBAR, from: (v) => v / HBAR, digits: 4 },
    thz: { label: "THz", to: (w) => (w / TAU) * 1000, from: (v) => (v * TAU) / 1000, digits: 2 },
    radfs: { label: "rad/fs", to: (w) => w, from: (v) => v, digits: 4 },
  };
  const cm = (v) => UNITS.cm.from(v);

  /* ======================= FFT ======================= */
  const tables = new Map();
  function fftTable(n) {
    let t = tables.get(n);
    if (!t) {
      const c = new Float64Array(n / 2), s = new Float64Array(n / 2), rev = new Uint32Array(n), bits = Math.round(Math.log2(n));
      for (let k = 0; k < n / 2; k++) { c[k] = Math.cos((TAU * k) / n); s[k] = Math.sin((TAU * k) / n); }
      for (let i = 0; i < n; i++) { let r = 0; for (let b = 0; b < bits; b++) r |= ((i >> b) & 1) << (bits - 1 - b); rev[i] = r; }
      t = { c, s, rev }; tables.set(n, t);
    }
    return t;
  }
  const isPow2 = (n) => n > 0 && (n & (n - 1)) === 0;
  const nextPow2 = (n) => { let p = 1; while (p < n) p <<= 1; return p; };
  /** In-place radix-2 FFT: X_k = Σ_j x_j exp(sign · 2πi jk/n), sign = ±1, n a power of two. */
  function fft(re, im, sign) {
    const n = re.length;
    if (n <= 1) return;
    if (!isPow2(n)) throw new Error("FFT length must be a power of two");
    const { c, s, rev } = fftTable(n);
    for (let i = 0; i < n; i++) { const j = rev[i]; if (j > i) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; } }
    for (let size = 2; size <= n; size <<= 1) {
      const half = size >> 1, step = n / size;
      for (let start = 0; start < n; start += size) {
        for (let k = 0; k < half; k++) {
          const wr = c[k * step], wi = sign * s[k * step], a = start + k, b = a + half;
          const xr = re[b] * wr - im[b] * wi, xi = re[b] * wi + im[b] * wr;
          re[b] = re[a] - xr; im[b] = im[a] - xi; re[a] += xr; im[a] += xi;
        }
      }
    }
  }
  /** Direct O(n²) DFT with the same definition (used by the tests). */
  function dft(re, im, sign) {
    const n = re.length, or = new Float64Array(n), oi = new Float64Array(n);
    for (let k = 0; k < n; k++) for (let j = 0; j < n; j++) {
      const a = (sign * TAU * j * k) / n, c = Math.cos(a), s = Math.sin(a);
      or[k] += re[j] * c - im[j] * s; oi[k] += re[j] * s + im[j] * c;
    }
    return { re: or, im: oi };
  }

  /* ======================= conventions ======================= */
  const STD = { R: [-1, 1], NR: [1, 1] };
  const CONVENTIONS = {
    standard: { name: "Standard: rephasing ω₁ flipped", R: [-1, 1], NR: [1, 1],
      text: "S_R = ∫∫ R_R e^{−iω₁t₁} e^{+iω₃t₃},  S_NR = ∫∫ R_NR e^{+iω₁t₁} e^{+iω₃t₃}: both pathways appear at (+ω₁, +ω₃)." },
    unflipped: { name: "Rephasing ω₁ not flipped", R: [1, 1], NR: [1, 1],
      text: "S_R = ∫∫ R_R e^{+iω₁t₁} e^{+iω₃t₃},  S_NR = ∫∫ R_NR e^{+iω₁t₁} e^{+iω₃t₃}: the rephasing signal appears at (−ω₁, +ω₃); the absorptive sum uses S_R(−ω₁)." },
    conjugate: { name: "Conjugate time convention (e^{+iωt})", R: [1, -1], NR: [-1, -1],
      text: "S_R = ∫∫ R_R e^{+iω₁t₁} e^{−iω₃t₃},  S_NR = ∫∫ R_NR e^{−iω₁t₁} e^{−iω₃t₃}: for data written with e^{+iωt}; responses written with e^{−iωt} (the built-in models) then appear at (−ω₁, −ω₃)." },
  };

  /* ======================= windows and the 2D transform ======================= */
  const WINDOWS = { none: "None (rectangular)", hann: "Half Hann (cos²)", hamming: "Half Hamming", gauss: "Half Gaussian (1% at t_max)" };
  /** Weights for samples t_i (fs) on [0, tmax]: apodisation × exp(−γ_b t/2) broadening × trapezoid first-point factor. */
  function taper(ts, o) {
    const n = ts.length, tmax = n > 1 ? ts[n - 1] + (ts[n - 1] - ts[n - 2]) : 1, w = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const x = Math.min(1, ts[i] / tmax);
      let v = 1;
      if (o.window === "hann") v = Math.cos((Math.PI * x) / 2) ** 2;
      else if (o.window === "hamming") v = 0.54 + 0.46 * Math.cos(Math.PI * x);
      else if (o.window === "gauss") v = Math.exp(-Math.log(100) * x * x);
      if (o.broaden > 0) v *= Math.exp(-0.5 * o.broaden * ts[i]);
      w[i] = v;
    }
    return w;
  }
  const uniformTimes = (n, dt) => Float64Array.from({ length: n }, (_, i) => i * dt);
  /**
   * 2D spectrum of a sampled response x (row-major, rows = t1, columns = t3) on a uniform grid with
   * step dt, by FFT. o = {sign1, sign3, off1, off3, pad, window, broaden (rad/fs FWHM), halfFirst}.
   * Returns {re, im, w1, w3, N1, N3} with ω ascending along both axes (row index = ω1).
   */
  function spectrum2D(x, n1, n3, dt, o) {
    const pad = o.pad || 1, N1 = nextPow2(n1) * pad, N3 = nextPow2(n3) * pad;
    const re = new Float64Array(N1 * N3), im = new Float64Array(N1 * N3);
    const w1 = taper(uniformTimes(n1, dt), o), w3 = taper(uniformTimes(n3, dt), o);
    if (o.halfFirst !== false) { w1[0] *= 0.5; w3[0] *= 0.5; }
    for (let i = 0; i < n1; i++) for (let k = 0; k < n3; k++) {
      const f = w1[i] * w3[k] * dt * dt, s = i * n3 + k, d = i * N3 + k;
      re[d] = x.re[s] * f; im[d] = x.im[s] * f;
    }
    const rr = new Float64Array(N3), ri = new Float64Array(N3);
    for (let i = 0; i < n1; i++) {
      const b = i * N3;
      for (let k = 0; k < N3; k++) { rr[k] = re[b + k]; ri[k] = im[b + k]; }
      fft(rr, ri, o.sign3);
      for (let k = 0; k < N3; k++) { re[b + k] = rr[k]; im[b + k] = ri[k]; }
    }
    const cr = new Float64Array(N1), ci = new Float64Array(N1);
    for (let k = 0; k < N3; k++) {
      for (let i = 0; i < N1; i++) { cr[i] = re[i * N3 + k]; ci[i] = im[i * N3 + k]; }
      fft(cr, ci, o.sign1);
      for (let i = 0; i < N1; i++) { re[i * N3 + k] = cr[i]; im[i * N3 + k] = ci[i]; }
    }
    // fftshift into ascending frequency order
    const outR = new Float64Array(N1 * N3), outI = new Float64Array(N1 * N3);
    for (let a = 0; a < N1; a++) {
      const ra = ((a + N1 / 2) % N1) * N3;
      for (let b = 0; b < N3; b++) { const s = ra + ((b + N3 / 2) % N3); outR[a * N3 + b] = re[s]; outI[a * N3 + b] = im[s]; }
    }
    return { re: outR, im: outI, w1: freqAxis(N1, dt, o.off1 || 0), w3: freqAxis(N3, dt, o.off3 || 0), N1, N3 };
  }
  /** Ascending angular-frequency axis of an N-point transform with step dt, shifted by off: off + 2π(a − N/2)/(N dt). */
  const freqAxis = (N, dt, off) => Float64Array.from({ length: N }, (_, a) => off + (TAU * (a - N / 2)) / (N * dt));
  /** Trapezoid weights for an arbitrary increasing grid (the last point carries a full step, matching the FFT path). */
  function trapezoid(ts, halfFirst) {
    const n = ts.length, w = new Float64Array(n);
    if (n === 1) { w[0] = 1; return w; }
    w[0] = (ts[1] - ts[0]) * (halfFirst === false ? 1 : 0.5);
    for (let i = 1; i < n - 1; i++) w[i] = (ts[i + 1] - ts[i - 1]) / 2;
    w[n - 1] = ts[n - 1] - ts[n - 2];
    return w;
  }
  /**
   * Direct (nonuniform) double Fourier sum on arbitrary increasing grids t1s, t3s, separable,
   * O(n1·n3·N3 + n1·N3·N1). Frequency grid as for the FFT with the reference step dtRef.
   */
  function spectrum2DDirect(x, t1s, t3s, dtRef, o) {
    const n1 = t1s.length, n3 = t3s.length, pad = o.pad || 1, N1 = nextPow2(n1) * pad, N3 = nextPow2(n3) * pad;
    const a1 = taper(t1s, o), a3 = taper(t3s, o), q1 = trapezoid(t1s, o.halfFirst), q3 = trapezoid(t3s, o.halfFirst);
    const W1 = freqAxis(N1, dtRef, 0), W3 = freqAxis(N3, dtRef, 0);
    // stage 1: along t3 for each t1 row
    const tr = new Float64Array(n1 * N3), ti = new Float64Array(n1 * N3);
    const er = new Float64Array(n3), ei = new Float64Array(n3);
    for (let b = 0; b < N3; b++) {
      for (let k = 0; k < n3; k++) { const ph = o.sign3 * W3[b] * t3s[k], f = a3[k] * q3[k]; er[k] = f * Math.cos(ph); ei[k] = f * Math.sin(ph); }
      for (let i = 0; i < n1; i++) {
        let sr = 0, si = 0; const base = i * n3;
        for (let k = 0; k < n3; k++) { const xr = x.re[base + k], xi = x.im[base + k]; sr += xr * er[k] - xi * ei[k]; si += xr * ei[k] + xi * er[k]; }
        tr[i * N3 + b] = sr; ti[i * N3 + b] = si;
      }
    }
    // stage 2: along t1
    const re = new Float64Array(N1 * N3), im = new Float64Array(N1 * N3), fr = new Float64Array(n1), fi = new Float64Array(n1);
    for (let a = 0; a < N1; a++) {
      for (let i = 0; i < n1; i++) { const ph = o.sign1 * W1[a] * t1s[i], f = a1[i] * q1[i]; fr[i] = f * Math.cos(ph); fi[i] = f * Math.sin(ph); }
      for (let b = 0; b < N3; b++) {
        let sr = 0, si = 0;
        for (let i = 0; i < n1; i++) { const yr = tr[i * N3 + b], yi = ti[i * N3 + b]; sr += yr * fr[i] - yi * fi[i]; si += yr * fi[i] + yi * fr[i]; }
        re[a * N3 + b] = sr; im[a * N3 + b] = si;
      }
    }
    return { re, im, w1: freqAxis(N1, dtRef, o.off1 || 0), w3: freqAxis(N3, dtRef, o.off3 || 0), N1, N3 };
  }
  /** Mirror the ω1 axis of a spectrum (Ω → −Ω on the periodic grid). */
  function mirror1(S) {
    const { N1, N3 } = S, re = new Float64Array(N1 * N3), im = new Float64Array(N1 * N3);
    for (let a = 0; a < N1; a++) { const s = ((N1 - a) % N1) * N3; for (let b = 0; b < N3; b++) { re[a * N3 + b] = S.re[s + b]; im[a * N3 + b] = S.im[s + b]; } }
    return { re, im, w1: Float64Array.from(S.w1, (_, a) => -S.w1[(N1 - a) % N1]), w3: S.w3, N1, N3 };
  }
  /** Total complex spectrum S_R + S_NR on the non-rephasing axes (mirroring the rephasing ω1 axis when the conventions differ). */
  function combine(SR, SNR, mirror) {
    const R = mirror ? mirror1(SR) : SR, n = SNR.re.length, re = new Float64Array(n), im = new Float64Array(n);
    for (let j = 0; j < n; j++) { re[j] = R.re[j] + SNR.re[j]; im[j] = R.im[j] + SNR.im[j]; }
    return { re, im, w1: SNR.w1, w3: SNR.w3, N1: SNR.N1, N3: SNR.N3 };
  }

  /* ======================= model responses ======================= */
  /** Isotropic orientational average ⟨(a·e1)(b·e2)(c·e3)(d·e4)⟩ ×5 for 2D/3D dipole vectors; pol "zzzz" or "zzyy". */
  function orient(a, b, c, d, pol) {
    const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + (u[2] || 0) * (v[2] || 0);
    const ab = dot(a, b), cd = dot(c, d), ac = dot(a, c), bd = dot(b, d), ad = dot(a, d), bc = dot(b, c);
    return pol === "zzyy" ? (4 * ab * cd - ac * bd - ad * bc) / 6 : (ab * cd + ac * bd + ad * bc) / 3;
  }
  /** FWHM → Gaussian standard deviation. */
  const fwhmToSigma = (f) => f / (2 * Math.sqrt(2 * Math.LN2));
  const cexp = (re, im) => { const m = Math.exp(re); return [m * Math.cos(im), m * Math.sin(im)]; };

  /**
   * Built-in models. Energies in rad/fs, times in fs. Each defines its parameters (for the UI)
   * and either paths(p, T) → Liouville pathways (separable), or grid(p, T, ...) → full response.
   * A pathway: {kind: GSB|SE|ESA, path: R|NR, amp, w1, w3, g1, g3, tf: [re, im] (waiting-time factor), label}.
   */
  const P = (k, label, type, def, extra) => Object.assign({ k, label, type, def }, extra || {});
  const MODELS = {
    tls: {
      name: "Two-level system",
      params: [P("w", "Transition energy ω_eg", "energy", cm(15000)), P("A", "Amplitude A (|μ|⁴)", "number", 1, { step: 0.1 }),
        P("tau2", "Dephasing time T₂ (Γ = 1/T₂)", "time", 50, { min: 1 }), P("T1", "Population lifetime T₁", "time", 1000, { min: 1 }),
        P("inh", "Static disorder, FWHM", "energy", 0, { min: 0 })],
      paths(p) {
        const G = 1 / p.tau2, tf = (T) => [Math.exp(-T / p.T1), 0], out = [];
        for (const path of ["R", "NR"]) for (const kind of ["GSB", "SE"]) out.push({ kind, path, amp: p.A, w1: p.w, w3: p.w, g1: G, g3: G, tf, label: `${kind} ${path === "R" ? "rephasing" : "non-rephasing"}` });
        return out;
      },
      transitions: (p) => [p.w],
      notes: "Bloch (Markovian) dephasing Γ = 1/T₂ on both coherence times, population decay e^{−T/T₁} during T (the ground-state bleach recovers as the excited state returns). Gives one diagonal peak; homogeneous linewidth FWHM = 2Γ on both axes. Optional static disorder: a Gaussian distribution of ω_eg.",
    },
    ladder: {
      name: "Three-level ladder: GSB, SE and ESA",
      params: [P("w", "Transition energy ω_eg", "energy", cm(15000)), P("anh", "ESA shift Δ (ω_fe = ω_eg − Δ)", "energy", cm(400)),
        P("kappa", "ESA strength |μ_fe|²/|μ_eg|²", "number", 1.5, { step: 0.1, min: 0 }), P("esaSign", "ESA sign", "sign", -1),
        P("tau2", "Dephasing time T₂ (g–e)", "time", 60, { min: 1 }), P("tauF", "Dephasing time (e–f)", "time", 60, { min: 1 }),
        P("T1", "Population lifetime T₁", "time", 1000, { min: 1 }), P("inh", "Static disorder, FWHM", "energy", 0, { min: 0 })],
      paths(p) {
        const G = 1 / p.tau2, Gf = 1 / p.tauF, wf = p.w - p.anh, tf = (T) => [Math.exp(-T / p.T1), 0], out = [];
        for (const path of ["R", "NR"]) {
          const n = path === "R" ? "rephasing" : "non-rephasing";
          out.push({ kind: "GSB", path, amp: 1, w1: p.w, w3: p.w, g1: G, g3: G, tf, label: `GSB ${n}` });
          out.push({ kind: "SE", path, amp: 1, w1: p.w, w3: p.w, g1: G, g3: G, tf, label: `SE ${n}` });
          out.push({ kind: "ESA", path, amp: p.esaSign * p.kappa, w1: p.w, w3: wf, g1: G, g3: Gf, tf, label: `ESA ${n}` });
        }
        return out;
      },
      transitions: (p) => [p.w, p.w - p.anh],
      notes: "Ground g, excited e and doubly excited f states. GSB and SE are positive at (ω_eg, ω_eg); ESA appears at (ω_eg, ω_fe) with the chosen sign (negative in a real measurement: it is extra absorption). With Δ = 0 and |μ_fe|² = 2|μ_eg|² (a harmonic oscillator) the ESA cancels GSB + SE exactly. Bloch dephasing, single T₁.",
    },
    dimer: {
      name: "Excitonic dimer (two coupled transitions)",
      params: [P("e1", "Site energy ε₁", "energy", cm(15000)), P("e2", "Site energy ε₂", "energy", cm(15600)), P("J", "Coupling J", "energy", cm(300)),
        P("mu1", "Dipole |μ₁|", "number", 1, { step: 0.1, min: 0 }), P("mu2", "Dipole |μ₂|", "number", 1, { step: 0.1, min: 0 }),
        P("phi", "Angle between dipoles", "angle", 60, { min: 0, max: 180 }), P("bx", "Biexciton binding shift Δ_b", "energy", 0),
        P("tau2", "Dephasing time T₂ (g–exciton)", "time", 60, { min: 1 }), P("tauF", "Dephasing time (exciton–f)", "time", 60, { min: 1 }),
        P("tauC", "Inter-exciton coherence time (during T)", "time", 200, { min: 1 }), P("T1", "Population lifetime T₁", "time", 1000, { min: 1 }),
        P("inh", "Static disorder, FWHM (correlated)", "energy", 0, { min: 0 }), P("pol", "Polarization", "pol", "zzzz")],
      exciton(p) {
        // H = [[ε1, J], [J, ε2]]; eigenvectors (cos θ, sin θ) and (−sin θ, cos θ), tan 2θ = 2J/(ε1 − ε2)
        const th = 0.5 * Math.atan2(2 * p.J, p.e1 - p.e2), c = Math.cos(th), s = Math.sin(th);
        const vec = [[c, s], [-s, c]], E = vec.map(([a, b]) => a * a * p.e1 + b * b * p.e2 + 2 * a * b * p.J);
        const m1 = [p.mu1, 0, 0], ang = (p.phi * Math.PI) / 180, m2 = [p.mu2 * Math.cos(ang), p.mu2 * Math.sin(ang), 0];
        const lin = (a, u, b, v) => [a * u[0] + b * v[0], a * u[1] + b * v[1], a * u[2] + b * v[2]];
        const mu = vec.map(([a, b]) => lin(a, m1, b, m2)); // g → exciton α
        const muf = vec.map(([a, b]) => lin(a, m2, b, m1)); // exciton α → f = |12⟩
        return { E, vec, mu, muf, Ef: p.e1 + p.e2 - p.bx, theta: th };
      },
      paths(p) {
        const X = this.exciton(p), G = 1 / p.tau2, Gf = 1 / p.tauF, Gc = 1 / p.tauC, pol = p.pol, out = [];
        const pop = (T) => [Math.exp(-T / p.T1), 0];
        const coh = (w) => (T) => cexp(-Gc * T, -w * T); // exp(−iωT − Γc T)
        const name = ["low", "high"], lo = X.E[0] <= X.E[1] ? 0 : 1, nm = (k) => (k === lo ? "−" : "+");
        for (let a = 0; a < 2; a++) for (let b = 0; b < 2; b++) {
          const A = X.mu[a], B = X.mu[b], Af = X.muf[a], Bf = X.muf[b], same = a === b;
          const tfR = same ? pop : coh(X.E[b] - X.E[a]), tfNR = same ? pop : coh(X.E[a] - X.E[b]);
          const tag = `α=${nm(a)}, β=${nm(b)}`;
          out.push({ kind: "GSB", path: "R", amp: orient(A, A, B, B, pol), w1: X.E[a], w3: X.E[b], g1: G, g3: G, tf: pop, label: `GSB rephasing, ${tag}` });
          out.push({ kind: "SE", path: "R", amp: orient(A, B, A, B, pol), w1: X.E[a], w3: X.E[b], g1: G, g3: G, tf: tfR, label: `SE rephasing, ${tag}` });
          out.push({ kind: "ESA", path: "R", amp: -orient(A, B, Bf, Af, pol), w1: X.E[a], w3: X.Ef - X.E[a], g1: G, g3: Gf, tf: tfR, label: `ESA rephasing, ${tag}` });
          out.push({ kind: "GSB", path: "NR", amp: orient(A, A, B, B, pol), w1: X.E[a], w3: X.E[b], g1: G, g3: G, tf: pop, label: `GSB non-rephasing, ${tag}` });
          out.push({ kind: "SE", path: "NR", amp: orient(A, B, B, A, pol), w1: X.E[a], w3: X.E[a], g1: G, g3: G, tf: tfNR, label: `SE non-rephasing, ${tag}` });
          out.push({ kind: "ESA", path: "NR", amp: -orient(A, B, Af, Bf, pol), w1: X.E[a], w3: X.Ef - X.E[b], g1: G, g3: Gf, tf: tfNR, label: `ESA non-rephasing, ${tag}` });
        }
        void name;
        return out;
      },
      transitions(p) { const X = this.exciton(p); return [X.E[0], X.E[1], X.Ef - X.E[0], X.Ef - X.E[1]]; },
      notes: "Frenkel exciton dimer: two sites with energies ε₁, ε₂ and coupling J form excitons E± by diagonalising [[ε₁, J], [J, ε₂]]; the two-exciton state |12⟩ has energy ε₁ + ε₂ − Δ_b. Transition dipoles are rotated into the exciton basis and averaged over random orientations (all-parallel or ⟨ZZYY⟩ polarization). All GSB, SE and ESA pathways are included. For J = 0 and Δ_b = 0 (with equal dephasing), the cross-peak pathways cancel exactly, as they must for independent transitions. Electronic coherences during T oscillate at |E₊ − E₋| on rephasing cross peaks and non-rephasing diagonal peaks. No population transfer between excitons (secular, uncoupled Bloch model).",
    },
    vibronic: {
      name: "Vibronic two-level system (displaced oscillator)",
      params: [P("w", "0–0 transition energy ω₀₀", "energy", cm(15000)), P("wv", "Vibrational frequency ω_v", "energy", cm(600)),
        P("S", "Huang–Rhys factor S", "number", 0.5, { step: 0.05, min: 0, max: 5 }), P("tauV", "Vibrational damping time", "time", 1000, { min: 1 }),
        P("tau2", "Electronic dephasing time T₂", "time", 60, { min: 1 }), P("T1", "Population lifetime T₁", "time", 1000, { min: 1 }),
        P("temp", "Temperature", "temp", 0, { min: 0, max: 2000 }), P("inh", "Static disorder, FWHM", "energy", 0, { min: 0 })],
      grid: true,
      transitions: (p) => [-2, -1, 0, 1, 2, 3].map((n) => p.w + n * p.wv),
      notes: "One underdamped vibrational mode linearly displaced in the excited state, treated exactly with the second-order cumulant (harmonic) response functions R₁–R₄ and the lineshape function g(t) = S[(n̄+1)(1 − e^{−(iω_v+γ_v)t}) + n̄(1 − e^{(iω_v−γ_v)t})], n̄ the thermal occupation. The damping γ_v is phenomenological. Electronic Bloch dephasing and T₁ are added. Peaks appear at ω₀₀ + nω_v (absorption) and ω₀₀ − nω_v (emission) with Franck–Condon weights; vibrational wavepackets make peaks oscillate in T with period 2π/ω_v. No Duschinsky rotation or frequency change.",
    },
  };
  const DEFAULT_PARAMS = (key) => Object.fromEntries(MODELS[key].params.map((q) => [q.k, q.def]));

  /** Rotating-frame phase per time step: R: (+1 on t1), NR: (−1 on t1); t3 always −1. */
  function pathwayResponse(paths, T, n1, n3, dt, wc, inh, comps) {
    const R = { re: new Float64Array(n1 * n3), im: new Float64Array(n1 * n3) }, NR = { re: new Float64Array(n1 * n3), im: new Float64Array(n1 * n3) };
    const e1r = new Float64Array(n1), e1i = new Float64Array(n1), e3r = new Float64Array(n3), e3i = new Float64Array(n3);
    for (const p of paths) {
      if (comps && comps[p.kind] === false) continue;
      const tf = p.tf(T), cr = p.amp * tf[0], ci = p.amp * tf[1];
      if (cr === 0 && ci === 0) continue;
      const a = p.path === "R" ? 1 : -1, out = p.path === "R" ? R : NR;
      for (let i = 0; i < n1; i++) { const t = i * dt, m = Math.exp(-p.g1 * t), ph = a * (p.w1 - wc) * t; e1r[i] = m * Math.cos(ph); e1i[i] = m * Math.sin(ph); }
      for (let k = 0; k < n3; k++) { const t = k * dt, m = Math.exp(-p.g3 * t), ph = -(p.w3 - wc) * t; e3r[k] = m * Math.cos(ph); e3i[k] = m * Math.sin(ph); }
      for (let i = 0; i < n1; i++) {
        const xr = cr * e1r[i] - ci * e1i[i], xi = cr * e1i[i] + ci * e1r[i], b = i * n3;
        for (let k = 0; k < n3; k++) { out.re[b + k] += xr * e3r[k] - xi * e3i[k]; out.im[b + k] += xr * e3i[k] + xi * e3r[k]; }
      }
    }
    applyDisorder(R, NR, n1, n3, dt, inh);
    return { R, NR };
  }
  /** Fully correlated Gaussian static disorder of all transition energies: R × exp(−σ²(t1 − t3)²/2), NR × exp(−σ²(t1 + t3)²/2). */
  function applyDisorder(R, NR, n1, n3, dt, inhFWHM) {
    if (!(inhFWHM > 0)) return;
    const s2 = fwhmToSigma(inhFWHM) ** 2;
    for (let i = 0; i < n1; i++) for (let k = 0; k < n3; k++) {
      const j = i * n3 + k, fR = Math.exp(-0.5 * s2 * ((i - k) * dt) ** 2), fN = Math.exp(-0.5 * s2 * ((i + k) * dt) ** 2);
      R.re[j] *= fR; R.im[j] *= fR; NR.re[j] *= fN; NR.im[j] *= fN;
    }
  }
  /** Vibronic (cumulant) response on the grid, rotating frame with carrier wc. comps: {GSB, SE}. */
  function vibronicResponse(p, T, n1, n3, dt, wc, comps) {
    const gam = 1 / p.tauV, ev = HBAR * p.wv, nb = p.temp > 0 ? 1 / (Math.exp(ev / (KB * p.temp)) - 1) : 0;
    // g(t) as [re, im]
    const g = (t) => {
      const d = Math.exp(-gam * t), c = Math.cos(p.wv * t), s = Math.sin(p.wv * t);
      // (n+1)(1 − d e^{−iωt}) + n(1 − d e^{+iωt})
      return [p.S * ((nb + 1) * (1 - d * c) + nb * (1 - d * c)), p.S * ((nb + 1) * (d * s) + nb * (-d * s))];
    };
    const m = n1 + n3, gA = [], gT = g(T), g1 = [], g3 = [], g1T = [], gT3 = [];
    for (let i = 0; i < n1; i++) { g1.push(g(i * dt)); g1T.push(g(i * dt + T)); }
    for (let k = 0; k < n3; k++) { g3.push(g(k * dt)); gT3.push(g(T + k * dt)); }
    for (let j = 0; j < m; j++) gA.push(g(j * dt + T));
    const R = { re: new Float64Array(n1 * n3), im: new Float64Array(n1 * n3) }, NR = { re: new Float64Array(n1 * n3), im: new Float64Array(n1 * n3) };
    const G = 1 / p.tau2, pop = Math.exp(-T / p.T1), d = p.w - wc;
    const useSE = !comps || comps.SE !== false, useGSB = !comps || comps.GSB !== false;
    for (let i = 0; i < n1; i++) for (let k = 0; k < n3; k++) {
      const a = g1[i], b = gT, c = g3[k], ab = g1T[i], bc = gT3[k], abc = gA[i + k];
      const el = Math.exp(-G * (i + k) * dt) * pop, j = i * n3 + k;
      // conjugates: x* = [x0, −x1]
      let rr = 0, ri = 0, nr = 0, ni = 0;
      if (useSE) {
        // R2 = exp(−a* + b − c* − ab* − bc + abc*)
        let E = cexp(-a[0] + b[0] - c[0] - ab[0] - bc[0] + abc[0], a[1] + b[1] + c[1] + ab[1] - bc[1] - abc[1]); rr += E[0]; ri += E[1];
        // R1 = exp(−a − b* − c* + ab + bc* − abc)
        E = cexp(-a[0] - b[0] - c[0] + ab[0] + bc[0] - abc[0], -a[1] + b[1] + c[1] + ab[1] - bc[1] - abc[1]); nr += E[0]; ni += E[1];
      }
      if (useGSB) {
        // R3 = exp(−a* + b* − c − ab* − bc* + abc*)
        let E = cexp(-a[0] + b[0] - c[0] - ab[0] - bc[0] + abc[0], a[1] - b[1] - c[1] + ab[1] + bc[1] - abc[1]); rr += E[0]; ri += E[1];
        // R4 = exp(−a − b − c + ab + bc − abc)
        E = cexp(-a[0] - b[0] - c[0] + ab[0] + bc[0] - abc[0], -a[1] - b[1] - c[1] + ab[1] + bc[1] - abc[1]); nr += E[0]; ni += E[1];
      }
      const t1 = i * dt, t3 = k * dt;
      const phR = d * (t1 - t3), phN = -d * (t1 + t3);
      const cR = Math.cos(phR), sR = Math.sin(phR), cN = Math.cos(phN), sN = Math.sin(phN);
      R.re[j] = el * (rr * cR - ri * sR); R.im[j] = el * (rr * sR + ri * cR);
      NR.re[j] = el * (nr * cN - ni * sN); NR.im[j] = el * (nr * sN + ni * cN);
    }
    applyDisorder(R, NR, n1, n3, dt, p.inh);
    return { R, NR };
  }
  /** Response of a built-in model at waiting time T on an n1 × n3 grid (rotating frame, carrier wc). */
  function modelResponse(key, p, T, n1, n3, dt, wc, comps) {
    const M = MODELS[key];
    if (!M) throw new Error("Unknown model " + key);
    if (M.grid) return Object.assign(vibronicResponse(p, T, n1, n3, dt, wc, comps), { paths: null });
    const paths = M.paths(p);
    return Object.assign(pathwayResponse(paths, T, n1, n3, dt, wc, p.inh, comps), { paths });
  }

  /* ======================= restricted expression language ======================= */
  const FUNCS = {
    exp: 1, log: 1, sqrt: 1, sin: 1, cos: 1, tan: 1, sinh: 1, cosh: 1, tanh: 1, abs: 1, arg: 1, re: 1, im: 1, conj: 1,
    cm: 1, ev: 1, thz: 1, pow: 2,
  };
  const CONSTS = { i: [0, 1], pi: [Math.PI, 0], e: [Math.E, 0] };
  const own = (o, k) => Object.prototype.hasOwnProperty.call(o, k); // never consult prototypes ("constructor", "__proto__", …)
  function tokenize(src) {
    const out = []; let i = 0;
    while (i < src.length) {
      const ch = src[i];
      if (/\s/.test(ch)) { i++; continue; }
      const num = /^(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?/.exec(src.slice(i));
      if (num) { out.push({ t: "num", v: parseFloat(num[0]), pos: i }); i += num[0].length; continue; }
      const id = /^[A-Za-z_][A-Za-z0-9_]*/.exec(src.slice(i));
      if (id) { out.push({ t: "id", v: id[0], pos: i }); i += id[0].length; continue; }
      if ("+-*/^(),".includes(ch)) { out.push({ t: ch, pos: i }); i++; continue; }
      throw new Error(`Unexpected character “${ch}” at position ${i + 1}`);
    }
    return out;
  }
  /** Parse to an AST; vars is the set of allowed variable names. Grammar: sum of products of unary powers. */
  function parseExpr(src, vars) {
    if (typeof src !== "string" || !src.trim()) throw new Error("Empty expression");
    if (src.length > 4000) throw new Error("Expression too long");
    const tk = tokenize(src); let p = 0, depth = 0;
    const peek = () => tk[p], take = (t) => { const x = tk[p]; if (!x || (t && x.t !== t)) throw new Error(`Expected “${t}” ${x ? "at position " + (x.pos + 1) : "at the end"}`); p++; return x; };
    function expr() {
      if (++depth > 200) throw new Error("Expression nested too deeply");
      let a = term();
      while (peek() && (peek().t === "+" || peek().t === "-")) { const op = take().t; a = { t: "bin", op, a, b: term() }; }
      depth--; return a;
    }
    function term() {
      let a = unary();
      while (peek() && (peek().t === "*" || peek().t === "/")) { const op = take().t; a = { t: "bin", op, a, b: unary() }; }
      if (peek() && (peek().t === "num" || peek().t === "id" || peek().t === "(")) throw new Error(`Missing operator before position ${peek().pos + 1} (write 2*x, not 2x)`);
      return a;
    }
    function unary() {
      if (peek() && peek().t === "-") { take(); return { t: "neg", a: unary() }; }
      if (peek() && peek().t === "+") { take(); return unary(); }
      return power();
    }
    function power() {
      const a = primary();
      if (peek() && peek().t === "^") { take(); return { t: "bin", op: "^", a, b: unary() }; }
      return a;
    }
    function primary() {
      const x = peek();
      if (!x) throw new Error("Unexpected end of expression");
      if (x.t === "num") { p++; return { t: "num", v: x.v }; }
      if (x.t === "(") { p++; const e = expr(); take(")"); return e; }
      if (x.t === "id") {
        p++;
        if (peek() && peek().t === "(") {
          if (!own(FUNCS, x.v)) throw new Error(`Unknown function “${x.v}”`);
          take("("); const args = [expr()];
          while (peek() && peek().t === ",") { take(","); args.push(expr()); }
          take(")");
          if (args.length !== FUNCS[x.v]) throw new Error(`${x.v}() takes ${FUNCS[x.v]} argument${FUNCS[x.v] > 1 ? "s" : ""}`);
          return { t: "call", f: x.v, args };
        }
        if (own(CONSTS, x.v)) return { t: "const", v: CONSTS[x.v] };
        if (!vars.has(x.v)) throw new Error(`Unknown name “${x.v}”`);
        return { t: "var", n: x.v };
      }
      throw new Error(`Unexpected “${x.t}” at position ${x.pos + 1}`);
    }
    const ast = expr();
    if (p < tk.length) throw new Error(`Unexpected “${tk[p].v || tk[p].t}” at position ${tk[p].pos + 1}`);
    return ast;
  }
  // complex helpers on [re, im]
  const C = {
    add: (a, b) => [a[0] + b[0], a[1] + b[1]], sub: (a, b) => [a[0] - b[0], a[1] - b[1]],
    mul: (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]],
    div: (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; },
    exp: (a) => cexp(a[0], a[1]), log: (a) => [Math.log(Math.hypot(a[0], a[1])), Math.atan2(a[1] || 0, a[0])], // principal branch, cut along the negative real axis
    pow: (a, b) => {
      if (b[1] === 0 && Number.isInteger(b[0]) && Math.abs(b[0]) <= 64) { let r = [1, 0], x = a, n = Math.abs(b[0]); while (n) { if (n & 1) r = C.mul(r, x); x = C.mul(x, x); n >>= 1; } return b[0] < 0 ? C.div([1, 0], r) : r; }
      if (a[0] === 0 && a[1] === 0) return [0, 0];
      return C.exp(C.mul(b, C.log(a)));
    },
  };
  const F = {
    exp: C.exp, log: C.log, sqrt: (a) => C.pow(a, [0.5, 0]),
    sin: (a) => [Math.sin(a[0]) * Math.cosh(a[1]), Math.cos(a[0]) * Math.sinh(a[1])],
    cos: (a) => [Math.cos(a[0]) * Math.cosh(a[1]), -Math.sin(a[0]) * Math.sinh(a[1])],
    tan: (a) => C.div(F.sin(a), F.cos(a)),
    sinh: (a) => [Math.sinh(a[0]) * Math.cos(a[1]), Math.cosh(a[0]) * Math.sin(a[1])],
    cosh: (a) => [Math.cosh(a[0]) * Math.cos(a[1]), Math.sinh(a[0]) * Math.sin(a[1])],
    tanh: (a) => C.div(F.sinh(a), F.cosh(a)),
    abs: (a) => [Math.hypot(a[0], a[1]), 0], arg: (a) => [Math.atan2(a[1], a[0]), 0],
    re: (a) => [a[0], 0], im: (a) => [a[1], 0], conj: (a) => [a[0], -a[1]],
    cm: (a) => [cm(a[0]), cm(a[1])], ev: (a) => [UNITS.ev.from(a[0]), UNITS.ev.from(a[1])], thz: (a) => [UNITS.thz.from(a[0]), UNITS.thz.from(a[1])],
    pow: C.pow,
  };
  /** Compile an AST to a function env → [re, im] (env maps names to [re, im]). */
  function compile(ast) {
    switch (ast.t) {
      case "num": { const v = [ast.v, 0]; return () => v; }
      case "const": { const v = ast.v; return () => v; }
      case "var": { const n = ast.n; return (env) => env[n]; }
      case "neg": { const a = compile(ast.a); return (env) => { const x = a(env); return [0 - x[0], 0 - x[1]]; }; } // 0 − x avoids −0, which would flip branch cuts
      case "call": { const f = F[ast.f], args = ast.args.map(compile); return args.length === 1 ? (env) => f(args[0](env)) : (env) => f(args[0](env), args[1](env)); }
      case "bin": {
        const a = compile(ast.a), b = compile(ast.b);
        const op = { "+": C.add, "-": C.sub, "*": C.mul, "/": C.div, "^": C.pow }[ast.op];
        return (env) => op(a(env), b(env));
      }
    }
    throw new Error("Bad expression");
  }
  const RESERVED = new Set(["t1", "t3", "T", ...Object.keys(CONSTS), ...Object.keys(FUNCS)]);
  /** Parse "name = expression" lines (constants for the response expressions). Returns {name: [re, im]}. */
  function parseParams(text) {
    const env = Object.create(null);
    (text || "").split(/\r?\n|;/).forEach((line, n) => {
      const s = line.replace(/#.*/, "").trim();
      if (!s) return;
      const m = /^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$/.exec(s);
      if (!m) throw new Error(`Parameter line ${n + 1}: write name = value`);
      if (RESERVED.has(m[1]) || m[1].startsWith("__")) throw new Error(`Parameter line ${n + 1}: “${m[1]}” is a reserved name`);
      const v = compile(parseExpr(m[2], new Set(Object.keys(env))))(env);
      if (!isFinite(v[0]) || !isFinite(v[1])) throw new Error(`Parameter line ${n + 1}: value is not finite`);
      env[m[1]] = v;
    });
    return env;
  }
  /** Evaluate lab-frame expressions for R_R and R_NR on the grid and move them to the rotating frame. */
  function exprResponse(src, T, n1, n3, dt, wc) {
    const params = parseParams(src.params), vars = new Set(["t1", "t3", "T", ...Object.keys(params)]);
    const out = {};
    for (const key of ["R", "NR"]) {
      const text = key === "R" ? src.R : src.NR;
      if (!text || !text.trim()) { out[key] = null; continue; }
      let fn;
      try { fn = compile(parseExpr(text, vars)); } catch (e) { throw new Error(`${key === "R" ? "Rephasing" : "Non-rephasing"} expression: ${e.message}`); }
      const env = Object.assign(Object.create(null), params, { T: [T, 0] }), arr = { re: new Float64Array(n1 * n3), im: new Float64Array(n1 * n3) };
      const s1 = STD[key][0];
      for (let i = 0; i < n1; i++) {
        const t1 = i * dt; env.t1 = [t1, 0];
        for (let k = 0; k < n3; k++) {
          const t3 = k * dt; env.t3 = [t3, 0];
          const v = fn(env);
          if (!isFinite(v[0]) || !isFinite(v[1])) throw new Error(`${key === "R" ? "Rephasing" : "Non-rephasing"} expression is not finite at t1 = ${t1} fs, t3 = ${t3} fs`);
          const ph = wc * (s1 * t1 + t3), c = Math.cos(ph), s = Math.sin(ph), j = i * n3 + k;
          arr.re[j] = v[0] * c - v[1] * s; arr.im[j] = v[0] * s + v[1] * c;
        }
      }
      out[key] = arr;
    }
    return out;
  }

  /* ======================= CSV import ======================= */
  const PATH_ALIASES = { r: "R", rephasing: "R", reph: "R", nr: "NR", nonrephasing: "NR", "non-rephasing": "NR", "non_rephasing": "NR" };
  /**
   * Parse a response-function CSV. Columns t1, T, t3, re, im and optionally pathway (R / NR),
   * by header name or in that order. Lines starting with # are comments; "# carrier_rad_per_fs = x"
   * sets the rotating-frame carrier. o = {timeUnit: "fs" | "ps", pathway: "R" | "NR" | "signedR" | "signedNR"}
   * (pathway applies when there is no pathway column; "signedR" means t1 < 0 holds the rephasing
   * response, "signedNR" the non-rephasing one). Returns {errors, warnings, meta, rows}.
   */
  function parseCSV(text, o) {
    o = o || {};
    const errors = [], warnings = [], meta = {}, rows = [];
    const scale = o.timeUnit === "ps" ? 1000 : 1;
    const lines = String(text || "").split(/\r?\n/);
    let cols = null, lineNo = 0, missing = [], bad = [];
    for (const raw of lines) {
      lineNo++;
      const line = raw.trim();
      if (!line) continue;
      if (line[0] === "#") {
        const m = /carrier_rad_per_fs\s*[=:]\s*([-+0-9.eE]+)/.exec(line); if (m) meta.carrier = parseFloat(m[1]);
        const u = /time_unit\s*[=:]\s*(fs|ps)/.exec(line); if (u) meta.timeUnit = u[1];
        continue;
      }
      const cells = line.split(/\s*[,;\t]\s*|\s+/).filter((c, i, a) => !(c === "" && i === a.length - 1));
      if (!cols) {
        if (cells.some((c) => /[A-Za-z]/.test(c) && !/^[-+]?(\d|\.\d)/.test(c) && !/^nan$/i.test(c))) {
          const names = cells.map((c) => c.toLowerCase().replace(/[^a-z0-9]/g, ""));
          const find = (...alts) => names.findIndex((n) => alts.some((a) => n === a || n.startsWith(a)));
          cols = { t1: find("t1", "tau"), T: names.findIndex((n) => n === "t" || n === "t2" || n.startsWith("tw") || n.startsWith("pop") || n === "waiting" || n.startsWith("tfs") || n.startsWith("tps")), t3: find("t3", "tdet"), re: find("re", "real"), im: find("im", "imag"), path: find("pathway", "path", "signal") };
          // "T" vs "t1"/"t3": header "T" lowercases to "t"
          const missingCols = ["t1", "T", "t3", "re", "im"].filter((k) => cols[k] < 0);
          if (missingCols.length) { errors.push(`Header is missing column${missingCols.length > 1 ? "s" : ""}: ${missingCols.join(", ")} (expected t1, T, t3, re, im[, pathway])`); return { errors, warnings, meta, rows }; }
          continue;
        }
        cols = cells.length >= 6 ? { t1: 0, T: 1, t3: 2, re: 3, im: 4, path: 5 } : { t1: 0, T: 1, t3: 2, re: 3, im: 4, path: -1 };
      }
      const need = Math.max(cols.t1, cols.T, cols.t3, cols.re, cols.im) + 1;
      if (cells.length < need || cells.slice(0, need).some((c) => c === "" || /^nan$/i.test(c))) { missing.push(lineNo); continue; }
      const num = (k) => Number(cells[cols[k]]);
      const r = { t1: num("t1") * scale, T: num("T") * scale, t3: num("t3") * scale, re: num("re"), im: num("im"), line: lineNo };
      if (![r.t1, r.T, r.t3, r.re, r.im].every(Number.isFinite)) { bad.push(lineNo); continue; }
      if (cols.path >= 0) {
        const pth = PATH_ALIASES[String(cells[cols.path] || "").toLowerCase()];
        if (!pth) { errors.push(`Line ${lineNo}: unknown pathway “${cells[cols.path]}” (use R or NR)`); continue; }
        r.path = pth;
      }
      rows.push(r);
    }
    const list = (a) => a.slice(0, 5).join(", ") + (a.length > 5 ? ` and ${a.length - 5} more` : "");
    if (missing.length) errors.push(`Missing values on line${missing.length > 1 ? "s" : ""} ${list(missing)}`);
    if (bad.length) errors.push(`Non-numeric or non-finite values on line${bad.length > 1 ? "s" : ""} ${list(bad)}`);
    if (!rows.length && !errors.length) errors.push("No data rows found");
    if (meta.timeUnit && meta.timeUnit !== (o.timeUnit || "fs")) warnings.push(`The file says its times are in ${meta.timeUnit}, but ${o.timeUnit || "fs"} is selected.`);
    // assign pathways
    for (const r of rows) {
      if (r.path) { if (r.t1 < 0) errors.push(`Line ${r.line}: negative t1 in a file with a pathway column`); continue; }
      const mode = o.pathway || "R";
      if (mode === "R" || mode === "NR") { r.path = mode; if (r.t1 < 0) { errors.push(`Line ${r.line}: negative t1; choose a signed-t1 file type if t1 < 0 holds the other pathway`); break; } }
      else { const neg = mode === "signedR" ? "R" : "NR", pos = neg === "R" ? "NR" : "R"; r.path = r.t1 < 0 ? neg : pos; r.t1 = Math.abs(r.t1); }
    }
    // with signed t1, the t1 = 0 row belongs to both pathways
    if (o.pathway === "signedR" || o.pathway === "signedNR") {
      const extra = rows.filter((r) => r.t1 === 0).map((r) => Object.assign({}, r, { path: r.path === "R" ? "NR" : "R" }));
      rows.push(...extra);
    }
    let tmax = 0; for (const r of rows) tmax = Math.max(tmax, Math.abs(r.t1), Math.abs(r.t3));
    if (tmax > 0 && tmax < 1 && (o.timeUnit || "fs") === "fs") warnings.push("All times are below 1 fs: are the times perhaps in ps?");
    return { errors, warnings, meta, rows };
  }
  /** Sorted unique values. */
  const uniq = (a) => [...new Set(a)].sort((x, y) => x - y);
  function checkAxis(ax, name, warnings, errors) {
    if (ax.length < 4) { errors.push(`${name} has only ${ax.length} distinct values; at least 4 are needed`); return { uniform: false, dt: NaN }; }
    if (ax[0] < 0) errors.push(`${name} has negative values`);
    if (ax[0] > 0) warnings.push(`${name} starts at ${+ax[0].toFixed(3)} fs, not 0: the one-sided transform assumes the response starts at time zero, so expect a phase twist`);
    const d = ax.slice(1).map((v, i) => v - ax[i]), mn = Math.min(...d), mx = Math.max(...d), mean = (ax[ax.length - 1] - ax[0]) / (ax.length - 1);
    const uniform = (mx - mn) / mean < 1e-6;
    return { uniform, dt: uniform ? mean : mn, mean, min: mn, max: mx };
  }
  /**
   * Build a validated dataset from parsed rows (one or more files). Every pathway and waiting time must
   * share one complete rectilinear (t1, t3) grid. Returns {errors, warnings, data} with
   * data = {t1, t3, Ts, R: [slices] | null, NR: [slices] | null, uniform, dt, carrier}.
   */
  function buildDataset(rowsList, o) {
    o = o || {};
    const errors = [], warnings = [], rows = [].concat(...rowsList);
    if (!rows.length) return { errors: ["No data"], warnings, data: null };
    const seen = new Set(), dups = [];
    for (const r of rows) { const key = `${r.path}|${r.T}|${r.t1}|${r.t3}`; if (seen.has(key)) dups.push(r.line); else seen.add(key); }
    if (dups.length) errors.push(`Duplicate coordinates (same pathway, T, t1, t3) on line${dups.length > 1 ? "s" : ""} ${dups.slice(0, 5).join(", ")}${dups.length > 5 ? ` and ${dups.length - 5} more` : ""}`);
    const t1 = uniq(rows.map((r) => r.t1)), t3 = uniq(rows.map((r) => r.t3)), Ts = uniq(rows.map((r) => r.T));
    const paths = uniq(rows.map((r) => r.path));
    const a1 = checkAxis(t1, "t1", warnings, errors), a3 = checkAxis(t3, "t3", warnings, errors);
    if (errors.length) return { errors, warnings, data: null };
    const i1 = new Map(t1.map((v, i) => [v, i])), i3 = new Map(t3.map((v, i) => [v, i])), n1 = t1.length, n3 = t3.length;
    const data = { t1: Float64Array.from(t1), t3: Float64Array.from(t3), Ts, R: null, NR: null, carrier: o.carrier || 0 };
    for (const pth of paths) {
      data[pth] = Ts.map(() => ({ re: new Float64Array(n1 * n3), im: new Float64Array(n1 * n3), count: 0 }));
    }
    const Tidx = new Map(Ts.map((v, i) => [v, i]));
    for (const r of rows) { const s = data[r.path][Tidx.get(r.T)], j = i1.get(r.t1) * n3 + i3.get(r.t3); s.re[j] = r.re; s.im[j] = r.im; s.count++; }
    for (const pth of paths) data[pth].forEach((s, k) => {
      if (s.count !== n1 * n3) errors.push(`${pth === "R" ? "Rephasing" : "Non-rephasing"} data at T = ${Ts[k]} fs: incomplete grid, ${s.count} of ${n1} × ${n3} = ${n1 * n3} points (every T slice and pathway must share one rectangular t1 × t3 grid)`);
    });
    if (errors.length) return { errors, warnings, data: null };
    // numerical sanity
    let peak = 0;
    for (const pth of paths) for (const s of data[pth]) for (let j = 0; j < s.re.length; j++) peak = Math.max(peak, Math.hypot(s.re[j], s.im[j]));
    if (!(peak > 0)) errors.push("The response is zero everywhere");
    else if (peak > 1e150 || peak < 1e-150) errors.push(`The response magnitude (${peak.toExponential(1)}) is outside the safe numerical range; rescale the data`);
    if (errors.length) return { errors, warnings, data: null };
    for (const pth of paths) {
      const s = data[pth][0]; let tail = 0;
      for (let i = 0; i < n1; i++) tail = Math.max(tail, Math.hypot(s.re[i * n3 + n3 - 1], s.im[i * n3 + n3 - 1]));
      for (let k = 0; k < n3; k++) tail = Math.max(tail, Math.hypot(s.re[(n1 - 1) * n3 + k], s.im[(n1 - 1) * n3 + k]));
      if (tail > 0.05 * peak) warnings.push(`${pth === "R" ? "Rephasing" : "Non-rephasing"} response is still ${Math.round((100 * tail) / peak)}% of its maximum at the last time point: truncation will cause sinc ripples; consider an apodisation window`);
    }
    if (!data.R) warnings.push("No rephasing data: only the non-rephasing spectrum can be computed (the absorptive spectrum needs both)");
    if (!data.NR) warnings.push("No non-rephasing data: only the rephasing spectrum can be computed (the absorptive spectrum needs both)");
    data.uniform = a1.uniform && a3.uniform;
    data.n1 = n1; data.n3 = n3;
    if (a1.uniform && a3.uniform && Math.abs(a1.dt - a3.dt) > 1e-9 * a1.dt) { data.uniform = false; warnings.push(`t1 and t3 have different steps (${+a1.dt.toFixed(4)} and ${+a3.dt.toFixed(4)} fs); the direct Fourier sum is used`); }
    data.dt = Math.min(a1.dt || a1.min, a3.dt || a3.min);
    if (!a1.uniform) warnings.push(`t1 is not uniformly sampled (steps ${+a1.min.toFixed(3)} to ${+a1.max.toFixed(3)} fs)`);
    if (!a3.uniform) warnings.push(`t3 is not uniformly sampled (steps ${+a3.min.toFixed(3)} to ${+a3.max.toFixed(3)} fs)`);
    return { errors, warnings, data };
  }
  /** Linear interpolation of a dataset slice from its (nonuniform) axes onto a uniform grid with step dt. */
  function resampleSlice(s, t1, t3, dt) {
    const m1 = Math.floor(t1[t1.length - 1] / dt + 1e-9) + 1, m3 = Math.floor(t3[t3.length - 1] / dt + 1e-9) + 1, n3 = t3.length;
    const locate = (ax, t) => { let lo = 0, hi = ax.length - 1; if (t <= ax[0]) return [0, 0]; if (t >= ax[hi]) return [hi - 1, 1]; while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (ax[mid] <= t) lo = mid; else hi = mid; } return [lo, (t - ax[lo]) / (ax[lo + 1] - ax[lo])]; };
    const re = new Float64Array(m1 * m3), im = new Float64Array(m1 * m3);
    for (let i = 0; i < m1; i++) {
      const [a, fa] = locate(t1, i * dt);
      for (let k = 0; k < m3; k++) {
        const [b, fb] = locate(t3, k * dt), j00 = a * n3 + b, j01 = j00 + 1, j10 = j00 + n3, j11 = j10 + 1, w00 = (1 - fa) * (1 - fb), w01 = (1 - fa) * fb, w10 = fa * (1 - fb), w11 = fa * fb, j = i * m3 + k;
        re[j] = w00 * s.re[j00] + w01 * s.re[j01] + w10 * s.re[j10] + w11 * s.re[j11];
        im[j] = w00 * s.im[j00] + w01 * s.im[j01] + w10 * s.im[j10] + w11 * s.im[j11];
      }
    }
    return { re, im, n1: m1, n3: m3 };
  }
  /** Export a response as CSV that parseCSV reads back (lab or rotating frame recorded in the header). */
  function responseCSV(resp, n1, n3, dt, T, carrier) {
    const lines = ["# 2DES response function exported by Indra's Webpage 2DES simulator", `# carrier_rad_per_fs = ${carrier}`, "# time_unit = fs",
      "# rotating frame: R_rot = R_lab exp(-i wc t1 + i wc t3) (rephasing), exp(+i wc t1 + i wc t3) (non-rephasing)", "pathway,t1,T,t3,re,im"];
    for (const key of ["R", "NR"]) {
      const x = resp[key]; if (!x) continue;
      for (let i = 0; i < n1; i++) for (let k = 0; k < n3; k++) { const j = i * n3 + k; lines.push(`${key},${+(i * dt).toPrecision(10)},${T},${+(k * dt).toPrecision(10)},${x.re[j].toPrecision(9)},${x.im[j].toPrecision(9)}`); }
    }
    return lines.join("\n") + "\n";
  }

  /* ======================= the full calculation ======================= */
  /**
   * job = {source: "model" | "expr" | "data", model, params, comps, expr: {params, R, NR}, data (from buildDataset),
   *        T, Tindex, numerics: {dt, n1, n3, pad, window, broaden, halfFirst, convention, carrier, nonuniform: "direct" | "resample"}}
   * Returns {R, NR, total (spectra or null), resp (time-domain, rotating frame), grid: {n1, n3, dt}, paths, warnings, info}.
   */
  function run(job) {
    const nm = job.numerics, conv = CONVENTIONS[nm.convention] || CONVENTIONS.standard, warnings = [];
    let n1 = nm.n1, n3 = nm.n3, dt = nm.dt, wc = nm.carrier || 0, resp, paths = null, direct = null;
    if (job.source === "model") {
      const m = modelResponse(job.model, job.params, job.T, n1, n3, dt, wc, job.comps); paths = m.paths; resp = { R: m.R, NR: m.NR };
      const half = Math.PI / dt, out = MODELS[job.model].transitions.call(MODELS[job.model], job.params).filter((w) => Math.abs(w - wc) > half);
      if (out.length) warnings.push(`Aliasing: ${out.length} transition${out.length > 1 ? "s lie" : " lies"} outside the sampled band ωc ± π/Δt. Reduce the time step or move the carrier.`);
    } else if (job.source === "expr") {
      resp = exprResponse(job.expr, job.T, n1, n3, dt, wc);
    } else {
      const D = job.data, k = Math.max(0, Math.min(D.Ts.length - 1, job.Tindex || 0));
      wc = D.carrier || 0;
      if (D.uniform) { n1 = D.n1; n3 = D.n3; dt = D.dt; resp = { R: D.R ? D.R[k] : null, NR: D.NR ? D.NR[k] : null }; }
      else if (nm.nonuniform === "resample") {
        resp = {}; for (const key of ["R", "NR"]) { if (!D[key]) { resp[key] = null; continue; } const s = resampleSlice(D[key][k], D.t1, D.t3, D.dt); resp[key] = s; n1 = s.n1; n3 = s.n3; }
        dt = D.dt; warnings.push(`Nonuniform data resampled onto a uniform ${+dt.toFixed(3)} fs grid by bilinear interpolation of the real and imaginary parts; interpolation acts as a low-pass filter.`);
      } else { resp = { R: D.R ? D.R[k] : null, NR: D.NR ? D.NR[k] : null }; direct = { t1: D.t1, t3: D.t3 }; n1 = D.n1; n3 = D.n3; dt = (D.t1[D.t1.length - 1] - D.t1[0]) / (D.n1 - 1); warnings.push("Nonuniform data: spectra computed by a direct nonuniform Fourier sum with trapezoid weights (no interpolation). The frequency grid uses the mean step; the alias-free band is set by the smallest step."); }
    }
    const opts = (key) => {
      const [s1, s3] = conv[key], [r1, r3] = STD[key];
      return { sign1: s1, sign3: s3, off1: (s1 / r1) * wc, off3: (s3 / r3) * wc, pad: nm.pad, window: nm.window, broaden: nm.broaden, halfFirst: nm.halfFirst };
    };
    const tr = (key) => {
      if (!resp[key]) return null;
      return direct ? spectrum2DDirect(resp[key], direct.t1, direct.t3, dt, opts(key)) : spectrum2D(resp[key], n1, n3, dt, opts(key));
    };
    const R = tr("R"), NR = tr("NR");
    const mirror = conv.R[0] / STD.R[0] !== conv.NR[0] / STD.NR[0];
    const total = R && NR ? combine(R, NR, mirror) : null;
    // truncation check on the time-domain data
    let peak = 0, tail = 0;
    for (const key of ["R", "NR"]) {
      const x = resp[key]; if (!x) continue;
      for (let j = 0; j < x.re.length; j++) peak = Math.max(peak, Math.hypot(x.re[j], x.im[j]));
      for (let i = 0; i < n1; i++) tail = Math.max(tail, Math.hypot(x.re[i * n3 + n3 - 1], x.im[i * n3 + n3 - 1]));
      for (let k = 0; k < n3; k++) tail = Math.max(tail, Math.hypot(x.re[(n1 - 1) * n3 + k], x.im[(n1 - 1) * n3 + k]));
    }
    if (peak > 0 && tail > 0.01 * peak && nm.window === "none") warnings.push(`Truncation: the response is still ${(100 * tail / peak).toFixed(1)}% of its maximum at t_max. Expect sinc ripples; extend t_max or apply a window.`);
    const N = (R || NR).N1, info = {
      nyquist: Math.PI / dt, gridSpacing1: TAU / ((R || NR).N1 * dt), gridSpacing3: TAU / ((R || NR).N3 * dt),
      resolution1: TAU / (n1 * dt), resolution3: TAU / (n3 * dt), tmax1: (n1 - 1) * dt, tmax3: (n3 - 1) * dt, N1: N, N3: (R || NR).N3, carrier: wc, mirror,
    };
    return { R, NR, total, resp, grid: { n1, n3, dt, direct: !!direct }, paths, warnings, info };
  }
  /** Value of S(ω1, ω3) for each pathway at one frequency point and a list of waiting times (direct sum; same windows). */
  function waitingTrace(job, w1, w3, Ts) {
    const out = { T: [], R: [], NR: [], total: [] }, conv = CONVENTIONS[job.numerics.convention] || CONVENTIONS.standard;
    const mirror = conv.R[0] / STD.R[0] !== conv.NR[0] / STD.NR[0];
    for (let q = 0; q < Ts.length; q++) {
      const j = Object.assign({}, job, { T: Ts[q], Tindex: q });
      const r = run(Object.assign({}, j, { numerics: Object.assign({}, j.numerics, { pad: 1 }) })), dt = r.grid.dt, nm = job.numerics, wc = r.info.carrier;
      const val = (key, ww1) => {
        const x = r.resp[key]; if (!x) return null;
        const [s1, s3] = conv[key], [q1, q3] = STD[key], O1 = ww1 - (s1 / q1) * wc, O3 = w3 - (s3 / q3) * wc;
        const n1 = r.grid.n1, n3 = r.grid.n3, ts1 = r.grid.direct ? job.data.t1 : uniformTimes(n1, dt), ts3 = r.grid.direct ? job.data.t3 : uniformTimes(n3, dt);
        const a1 = taper(ts1, nm), a3 = taper(ts3, nm), q1w = trapezoid(ts1, nm.halfFirst), q3w = trapezoid(ts3, nm.halfFirst);
        let sr = 0, si = 0;
        for (let i = 0; i < n1; i++) {
          let rr = 0, ri = 0;
          for (let k = 0; k < n3; k++) { const ph = s3 * O3 * ts3[k], f = a3[k] * q3w[k], c = f * Math.cos(ph), s = f * Math.sin(ph), jj = i * n3 + k; rr += x.re[jj] * c - x.im[jj] * s; ri += x.re[jj] * s + x.im[jj] * c; }
          const ph = s1 * O1 * ts1[i], f = a1[i] * q1w[i], c = f * Math.cos(ph), s = f * Math.sin(ph);
          sr += rr * c - ri * s; si += rr * s + ri * c;
        }
        return [sr, si];
      };
      const R = val("R", mirror ? -w1 : w1), NR = val("NR", w1);
      out.T.push(Ts[q]); out.R.push(R); out.NR.push(NR); out.total.push(R && NR ? [R[0] + NR[0], R[1] + NR[1]] : null);
    }
    return out;
  }
  /**
   * Convergence check: compare the total (or available) spectrum at (dt, n) with (dt/2, 2n) [sampling]
   * and (dt, 2n) [window length] on exactly coincident frequency points inside `range`
   * ({w1: [lo, hi], w3: [lo, hi]} in rad/fs; default: the central half of the sampled band).
   * Returns the largest changes relative to the largest |S| in that region.
   */
  function convergence(job, range) {
    const nm = job.numerics, base = run(job);
    const pick = (r) => r.total || r.R || r.NR;
    const B = pick(base);
    const half = Math.PI / (2 * base.grid.dt), c1 = (B.w1[0] + B.w1[B.N1 - 1]) / 2, c3 = (B.w3[0] + B.w3[B.N3 - 1]) / 2;
    const r1 = range && range.w1 ? range.w1 : [c1 - half, c1 + half], r3 = range && range.w3 ? range.w3 : [c3 - half, c3 + half];
    const finer = run(Object.assign({}, job, { numerics: Object.assign({}, nm, { dt: nm.dt / 2, n1: nm.n1 * 2, n3: nm.n3 * 2 }) }));
    const longer = run(Object.assign({}, job, { numerics: Object.assign({}, nm, { n1: nm.n1 * 2, n3: nm.n3 * 2 }) }));
    const F_ = pick(finer), L = pick(longer);
    let dF = 0, dL = 0, maxB = 0;
    for (let a = 0; a < B.N1; a++) {
      if (B.w1[a] < r1[0] || B.w1[a] > r1[1]) continue;
      for (let b = 0; b < B.N3; b++) {
        if (B.w3[b] < r3[0] || B.w3[b] > r3[1]) continue;
        const j = a * B.N3 + b; maxB = Math.max(maxB, Math.hypot(B.re[j], B.im[j]));
        const jf = (a + B.N1 / 2) * F_.N3 + (b + B.N3 / 2); dF = Math.max(dF, Math.hypot(F_.re[jf] - B.re[j], F_.im[jf] - B.im[j]));
        const jl = 2 * a * L.N3 + 2 * b; dL = Math.max(dL, Math.hypot(L.re[jl] - B.re[j], L.im[jl] - B.im[j]));
      }
    }
    return { sampling: maxB ? dF / maxB : NaN, window: maxB ? dL / maxB : NaN };
  }

  /* ======================= small helpers for analysis and tests ======================= */
  /** Full width at half maximum of a sampled peak (linear interpolation of the half-maximum crossings). */
  function fwhm(xs, ys) {
    let im = 0; for (let i = 1; i < ys.length; i++) if (ys[i] > ys[im]) im = i;
    const h = ys[im] / 2; let l = im, r = im;
    while (l > 0 && ys[l] > h) l--;
    while (r < ys.length - 1 && ys[r] > h) r++;
    const xl = xs[l] + ((h - ys[l]) * (xs[l + 1] - xs[l])) / (ys[l + 1] - ys[l]);
    const xr = xs[r - 1] + ((h - ys[r - 1]) * (xs[r] - xs[r - 1])) / (ys[r] - ys[r - 1]);
    return xr - xl;
  }
  /** Index of the nearest axis value. */
  const nearest = (ax, v) => { let b = 0; for (let i = 1; i < ax.length; i++) if (Math.abs(ax[i] - v) < Math.abs(ax[b] - v)) b = i; return b; };

  return {
    TAU, HBAR, C_CM_FS, UNITS, cm, fft, dft, nextPow2, CONVENTIONS, STD, WINDOWS, taper, spectrum2D, spectrum2DDirect, freqAxis, trapezoid,
    mirror1, combine, orient, MODELS, DEFAULT_PARAMS, modelResponse, pathwayResponse, vibronicResponse, applyDisorder, fwhmToSigma,
    tokenize, parseExpr, compile, parseParams, exprResponse, parseCSV, buildDataset, resampleSlice, responseCSV, run, waitingTrace, convergence, fwhm, nearest,
  };
});
