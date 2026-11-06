/*
 * tmm-engine.js — transfer-matrix method (TMM) for planar, isotropic multilayers,
 * plus coupled-oscillator tools for polariton analysis. Used by /tools/tmm/.
 * Pure functions, no DOM access; unit-tested in Node (tests/tmm-engine.test.js).
 *
 * Convention: time dependence e^{−iωt}, complex index ñ = n + ik (k ≥ 0 absorbs),
 * permittivity ε = ñ². Light is incident from a transparent medium of index n₀.
 *
 *   in-plane wave vector (conserved):   k_x = k₀ n₀ sin θ₀,          k₀ = 2π/λ
 *   normal wave vector in layer j:      k_z,j = √(k₀² ε_j − k_x²),   Im k_z,j ≥ 0
 *   optical admittance:                 η_j = k_z,j (s, TE),   η_j = ε_j / k_z,j (p, TM)
 *   phase thickness:                    δ_j = k_z,j d_j
 *   characteristic (Abelès) matrix:     M_j = [[cos δ_j, −i sin δ_j / η_j], [−i η_j sin δ_j, cos δ_j]]
 *   [B, C]ᵀ = (Π_j M_j) [1, η_sub]ᵀ  (tangential E and H at the front surface)
 *   r = (η₀B − C)/(η₀B + C),  t = 2η₀/(η₀B + C)
 *   R = |r|²,  T = Re(η_sub)/Re(η₀) · |t|²,  A = 1 − R − T
 * Internally k_z is expressed in units of k₀ (q = k_z/k₀ = √(ε − n₀² sin²θ₀)).
 * References: M. Born & E. Wolf, Principles of Optics, 7th ed. (CUP, 1999), §1.6;
 * P. Yeh, Optical Waves in Layered Media (Wiley, 1988); H. A. Macleod, Thin-Film Optical Filters.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.TMM = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /* ---------------- complex numbers: [re, im] ----------------
   * Small, allocation-light helpers. The hot loop in solve() works on plain
   * scalars and does not allocate per layer beyond the admittance.
   */
  const cx = (re, im) => [re, im || 0];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1]];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1]];
  const mul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
  const div = (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; };
  const conj = (a) => [a[0], -a[1]];
  const abs2 = (a) => a[0] * a[0] + a[1] * a[1];
  /** principal square root (Re ≥ 0) */
  const csqrt = (a) => { const r = Math.hypot(a[0], a[1]); const re = Math.sqrt(Math.max(0, (r + a[0]) / 2)); const im = Math.sqrt(Math.max(0, (r - a[0]) / 2)); return [re, a[1] < 0 ? -im : im]; };
  const ccos = (a) => [Math.cos(a[0]) * Math.cosh(a[1]), -Math.sin(a[0]) * Math.sinh(a[1])];
  const csin = (a) => [Math.sin(a[0]) * Math.cosh(a[1]), Math.cos(a[0]) * Math.sinh(a[1])];
  const cexp = (a) => { const e = Math.exp(a[0]); return [e * Math.cos(a[1]), e * Math.sin(a[1])]; };
  const scale = (a, s) => [a[0] * s, a[1] * s];

  /* ---------------- materials ----------------
   * Each material is { kind, nk(lambda_nm) -> [n, k] } with k ≥ 0 for absorption.
   */
  const HC_EV_NM = 1239.841984; // E[eV] = 1239.84 / λ[nm]

  /** Refractive index from a complex permittivity ε (Im ε ≥ 0). */
  function nkFromEps(eps) {
    const s = csqrt(eps);
    return s[1] < 0 ? [-s[0], -s[1]] : s;
  }
  const epsFromNk = (nk) => [nk[0] * nk[0] - nk[1] * nk[1], 2 * nk[0] * nk[1]];

  const material = {
    constant: (n, k) => ({ kind: "constant", nk: () => [n, k || 0] }),
    /** ε(ν̃) = ε∞ + Σ f·ν₀² / (ν₀² − ν̃² − iγν̃), ν in cm⁻¹ (works for any spectral unit). */
    lorentz: (epsInf, oscillators) => ({
      kind: "lorentz",
      eps: (lambda_nm) => {
        const nu = 1e7 / lambda_nm;
        let e = cx(epsInf, 0);
        for (const o of oscillators) {
          const num = cx(o.f * o.nu0 * o.nu0, 0);
          const den = cx(o.nu0 * o.nu0 - nu * nu, -o.gamma * nu);
          e = add(e, div(num, den));
        }
        return e;
      },
      nk(lambda_nm) { return nkFromEps(this.eps(lambda_nm)); },
    }),
    /** Drude metal: ε(E) = ε∞ − Ep² / (E² + iγE), energies in eV (ħω_p = Ep, ħγ_p = γ). */
    drude: (epsInf, Ep_eV, gamma_eV) => ({
      kind: "drude",
      eps: (lambda_nm) => {
        const E = HC_EV_NM / lambda_nm;
        return sub(cx(epsInf, 0), div(cx(Ep_eV * Ep_eV, 0), cx(E * E, gamma_eV * E)));
      },
      nk(lambda_nm) { return nkFromEps(this.eps(lambda_nm)); },
    }),
    /** Real index from any function n(λ in µm), e.g. a Sellmeier equation. */
    fromFunction: (nOfMicron, kConst) => ({ kind: "function", nk: (lambda_nm) => [nOfMicron(lambda_nm / 1000), kConst || 0] }),
  };

  /* ---------------- core solver ---------------- */
  /**
   * Normalised normal wave vector q = k_z/k₀ = √(ε − β²), β = n₀ sin θ₀,
   * on the physical branch: Im q > 0 (decaying), or Re q ≥ 0 when q is real (outgoing).
   */
  function kzNorm(eps, beta2) {
    let q = csqrt([eps[0] - beta2, eps[1]]);
    if (q[1] < 0 || (q[1] === 0 && q[0] < 0)) q = [-q[0], -q[1]];
    if (q[0] === 0 && q[1] === 0) q = [1e-12, 0]; // grazing / exactly critical: avoid 1/0
    return q;
  }
  /** Admittance η (in units of k₀ for s; dimensionless ratio form for p). */
  function admittance(eps, q, pol) { return pol === "p" ? div(eps, q) : q; }

  /**
   * @param {object} stack {incident: material, layers: [{material, d_nm}], substrate: material}
   * @param {number} lambda_nm vacuum wavelength
   * @param {number} theta_deg angle of incidence in the incident medium
   * @param {"s"|"p"} pol polarisation
   * @returns {R, T, A, r, t}
   */
  function solve(stack, lambda_nm, theta_deg, pol) {
    const nk0 = stack.incident.nk(lambda_nm);
    if (nk0[1] !== 0) throw new Error("The incident medium must be transparent.");
    const beta = nk0[0] * Math.sin((theta_deg * Math.PI) / 180);
    const beta2 = beta * beta;
    const k0 = (2 * Math.PI) / lambda_nm;
    const eps0 = epsFromNk(nk0), epsS = epsFromNk(stack.substrate.nk(lambda_nm));
    const eta0 = admittance(eps0, kzNorm(eps0, beta2), pol);
    const etaS = admittance(epsS, kzNorm(epsS, beta2), pol);
    // M = Π M_j as scalars (m11 … m22, real and imaginary parts)
    let ar = 1, ai = 0, br = 0, bi = 0, cr = 0, ci = 0, dr = 1, di = 0;
    for (const L of stack.layers) {
      const eps = epsFromNk(L.material.nk(lambda_nm));
      const q = kzNorm(eps, beta2), eta = admittance(eps, q, pol);
      const dl = [q[0] * k0 * L.d_nm, q[1] * k0 * L.d_nm];
      const c = ccos(dl), s = csin(dl);
      // e12 = −i s/η, e21 = −i η s
      const sEta = div(s, eta), etaS_ = mul(eta, s);
      const e11r = c[0], e11i = c[1], e12r = sEta[1], e12i = -sEta[0], e21r = etaS_[1], e21i = -etaS_[0];
      const n11r = ar * e11r - ai * e11i + br * e21r - bi * e21i, n11i = ar * e11i + ai * e11r + br * e21i + bi * e21r;
      const n12r = ar * e12r - ai * e12i + br * e11r - bi * e11i, n12i = ar * e12i + ai * e12r + br * e11i + bi * e11r;
      const n21r = cr * e11r - ci * e11i + dr * e21r - di * e21i, n21i = cr * e11i + ci * e11r + dr * e21i + di * e21r;
      const n22r = cr * e12r - ci * e12i + dr * e11r - di * e11i, n22i = cr * e12i + ci * e12r + dr * e11i + di * e11r;
      ar = n11r; ai = n11i; br = n12r; bi = n12i; cr = n21r; ci = n21i; dr = n22r; di = n22i;
    }
    const B = add([ar, ai], mul([br, bi], etaS));
    const C = add([cr, ci], mul([dr, di], etaS));
    const den = add(mul(eta0, B), C);
    const r = div(sub(mul(eta0, B), C), den);
    const t = div(scale(eta0, 2), den);
    const R = abs2(r);
    const T = (etaS[0] / eta0[0]) * abs2(t);
    return { R, T, A: 1 - R - T, r, t };
  }

  /** Unpolarised average of s and p. */
  function solveUnpolarised(stack, lambda_nm, theta_deg) {
    const s = solve(stack, lambda_nm, theta_deg, "s"), p = solve(stack, lambda_nm, theta_deg, "p");
    return { R: (s.R + p.R) / 2, T: (s.T + p.T) / 2, A: (s.A + p.A) / 2 };
  }
  const solveAny = (stack, l, th, pol) => (pol === "u" ? solveUnpolarised(stack, l, th) : solve(stack, l, th, pol));

  /**
   * Tangential electric-field intensity |E(z)|² / |E_inc|² through the stack
   * (z = 0 at the first interface). Sampled every dz_nm; includes `tail_nm` of substrate.
   * Fields are obtained by back-propagating [E, H] = [1, η_sub] from the substrate.
   */
  function fieldProfile(stack, lambda_nm, theta_deg, pol, dz_nm, tail_nm) {
    const n0 = stack.incident.nk(lambda_nm)[0];
    const beta = n0 * Math.sin((theta_deg * Math.PI) / 180), beta2 = beta * beta;
    const k0 = (2 * Math.PI) / lambda_nm;
    const lay = (m) => { const eps = epsFromNk(m.nk(lambda_nm)); const q = kzNorm(eps, beta2); return { q, eta: admittance(eps, q, pol), k0 }; };
    const subA = lay(stack.substrate);
    const layers = stack.layers.map((L) => ({ ...lay(L.material), d: L.d_nm }));
    const back = new Array(layers.length);
    let EH = [cx(1, 0), subA.eta];
    for (let j = layers.length - 1; j >= 0; j--) {
      back[j] = EH;
      EH = propagate(layers[j], layers[j].d, EH);
    }
    const eta0 = lay(stack.incident).eta;
    // E(0) = E_inc(1 + r), H(0) = η₀ E_inc (1 − r)  ⇒  E_inc = (η₀E(0) + H(0)) / (2η₀)
    const Einc = div(add(mul(eta0, EH[0]), EH[1]), scale(eta0, 2));
    const norm = abs2(Einc);
    const z = [], I2 = [], layerIndex = [];
    let z0 = 0;
    layers.forEach((L, j) => {
      const n = Math.max(1, Math.ceil(L.d / dz_nm));
      for (let i = 0; i <= n; i++) {
        const depth = (L.d * i) / n;
        const f = propagate(L, L.d - depth, back[j]);
        z.push(z0 + depth); I2.push(abs2(f[0]) / norm); layerIndex.push(j);
      }
      z0 += L.d;
    });
    // substrate: transmitted wave e^{i k_z z}, |·|² = e^{−2 k₀ Im(q) z}
    const nt = Math.max(1, Math.ceil((tail_nm || 0) / dz_nm));
    for (let i = 1; i <= nt; i++) {
      const dz = ((tail_nm || 0) * i) / nt;
      z.push(z0 + dz); I2.push(Math.exp(-2 * k0 * subA.q[1] * dz) / norm); layerIndex.push(layers.length);
    }
    return { z, I2, layerIndex, total_nm: z0 };
  }
  /** [E, H] a distance `dist` in front of a point where the fields are EH (same characteristic matrix as solve). */
  function propagate(L, dist, EH) {
    const dl = scale(L.q, L.k0 * dist);
    const c = ccos(dl), s = csin(dl);
    const mi = [0, -1];
    return [add(mul(c, EH[0]), mul(div(mul(mi, s), L.eta), EH[1])), add(mul(mul(mi, L.eta), mul(s, EH[0])), mul(c, EH[1]))];
  }

  /* ---------------- structure helpers ---------------- */
  /** Quarter-wave pairs (high/low) designed for lambda0: returns layers [H, L, H, L, …]. */
  function dbr(nH, nL, lambda0_nm, pairs, highFirst) {
    const H = { material: material.constant(nH, 0), d_nm: lambda0_nm / (4 * nH), label: "H" };
    const L = { material: material.constant(nL, 0), d_nm: lambda0_nm / (4 * nL), label: "L" };
    const out = [];
    for (let i = 0; i < pairs; i++) out.push(...(highFirst === false ? [L, H] : [H, L]));
    return out;
  }

  /* ---------------- maps ---------------- */
  /** quantity ∈ "R" | "T" | "A". Returns {x[], E_eV[], Z[x][E]} for x = angle (deg). */
  function angleMap(stack, pol, quantity, thetas, energies_eV) {
    const Z = thetas.map((th) => energies_eV.map((E) => solveAny(stack, HC_EV_NM / E, th, pol)[quantity]));
    return { x: thetas, E: energies_eV, Z };
  }
  /**
   * Map versus in-plane wave vector k∥ (µm⁻¹) instead of angle: k∥ = (2π/λ) n₀ sin θ.
   * Points outside the light cone of the incident medium (sin θ > 1) are NaN.
   */
  function kMap(stack, pol, quantity, kpar_um, energies_eV) {
    const Z = kpar_um.map((k) => energies_eV.map((E) => {
      const lam = HC_EV_NM / E, n0 = stack.incident.nk(lam)[0];
      const st = (k * lam) / (2 * Math.PI * 1000 * n0);
      return Math.abs(st) > 1 ? NaN : solveAny(stack, lam, (Math.asin(st) * 180) / Math.PI, pol)[quantity];
    }));
    return { x: kpar_um, E: energies_eV, Z, axis: "k" };
  }
  /** k∥ (µm⁻¹) for an energy (eV) and angle (deg) in a medium of index n₀. */
  const kParallel = (E_eV, theta_deg, n0) => ((2 * Math.PI * 1000) / (HC_EV_NM / E_eV)) * (n0 || 1) * Math.sin((theta_deg * Math.PI) / 180);
  /** R, T, A against angle at one wavelength (e.g. surface-plasmon resonance in the Kretschmann geometry). */
  function angleScan(stack, lambda_nm, thetas, pol) {
    const R = [], T = [], A = [];
    thetas.forEach((th) => { const r = solveAny(stack, lambda_nm, th, pol); R.push(r.R); T.push(r.T); A.push(r.A); });
    return { theta: thetas, R, T, A };
  }
  /** Repeat a block of layers: (L₁ … L_m) × times. */
  const repeatBlock = (layers, times) => Array.from({ length: Math.max(0, Math.round(times)) }, () => layers.map((L) => ({ ...L }))).flat();
  /** Rabi splitting of a cavity completely filled with a Lorentz medium (ε∞, f) at zero detuning: ħΩ_R ≈ E₀ √(f/ε∞). */
  const filledCavityRabi = (E0, f, epsInf) => E0 * Math.sqrt(f / epsInf);

  /** Map versus the thickness of layer `layerIdx` (normal or fixed angle). */
  function thicknessMap(stack, layerIdx, thicknesses_nm, theta_deg, pol, quantity, energies_eV) {
    const Z = thicknesses_nm.map((d) => {
      const layers = stack.layers.map((L, j) => (j === layerIdx ? { ...L, d_nm: d } : L));
      const st = { ...stack, layers };
      return energies_eV.map((E) => solveAny(st, HC_EV_NM / E, theta_deg, pol)[quantity]);
    });
    return { x: thicknesses_nm, E: energies_eV, Z };
  }

  /* ---------------- polariton analysis ---------------- */
  /** Cavity photon dispersion E_c(θ) = E₀ / √(1 − (n₀ sinθ / n_eff)²), θ in the incident medium of index n₀ (default 1). */
  const cavityEnergy = (E0, neff, theta_deg, n0) => E0 / Math.sqrt(1 - ((n0 || 1) * Math.sin((theta_deg * Math.PI) / 180)) ** 2 / (neff * neff));

  /** Two coupled oscillators (cavity Ec, exciton/vibration Ex, coupling g = ħΩ_R/2). */
  function coupledOscillator(Ec, Ex, g) {
    const mean = (Ec + Ex) / 2, half = Math.sqrt(g * g + ((Ec - Ex) / 2) ** 2);
    const Eup = mean + half, Elp = mean - half;
    // Hopfield coefficients of the lower polariton: |X|² matter, |C|² photon
    const det = Ec - Ex;
    const X2 = 0.5 * (1 + det / (2 * half)); // matter fraction of the LP
    return { Eup, Elp, LP_matter: X2, LP_photon: 1 - X2, UP_matter: 1 - X2, UP_photon: X2 };
  }

  /** Local extrema along one spectrum. kind "max" (T, A) or "min" (R). Returns indices. */
  function findPeaks(y, kind, minProminence) {
    const s = kind === "min" ? y.map((v) => -v) : y.slice();
    const peaks = [];
    for (let i = 1; i < s.length - 1; i++) {
      if (s[i] > s[i - 1] && s[i] >= s[i + 1]) {
        // prominence: drop to the lowest point before a higher peak on either side
        let l = s[i], r = s[i];
        for (let j = i - 1; j >= 0 && s[j] <= s[i]; j--) l = Math.min(l, s[j]);
        for (let j = i + 1; j < s.length && s[j] <= s[i]; j++) r = Math.min(r, s[j]);
        if (s[i] - Math.max(l, r) >= (minProminence || 0)) peaks.push(i);
      }
    }
    return peaks;
  }

  /** Extract resonance positions from a map: [{x, E}] (strongest `perColumn` features per column). */
  function extractResonances(map, quantity, minProminence, perColumn) {
    const kind = quantity === "R" ? "min" : "max";
    const pts = [];
    map.x.forEach((x, i) => {
      const col = map.Z[i];
      let idx = findPeaks(col, kind, minProminence);
      idx.sort((a, b) => (kind === "min" ? col[a] - col[b] : col[b] - col[a]));
      idx.slice(0, perColumn || 2).forEach((k) => pts.push({ x, E: map.E[k] }));
    });
    return pts;
  }

  /**
   * Follow the lower and upper polariton branches across the map by continuity.
   * Seeds at the first column with the nearest features below/above Ex (ignoring
   * |E − Ex| < exclude, the uncoupled absorption line), then at each column takes the
   * feature nearest a linear extrapolation of the branch. A branch stops when no
   * feature lies within the tolerance, so other cavity orders are not picked up.
   * seedWithin: a branch must start within this many columns of the first (default: any).
   */
  function trackPolaritons(map, quantity, minProminence, Ex, exclude, seedWithin) {
    const kind = quantity === "R" ? "min" : "max";
    const dE = Math.abs(map.E[1] - map.E[0]) || 1e-6;
    const cols = map.x.map((x, i) => findPeaks(map.Z[i], kind, minProminence).map((k) => map.E[k]).filter((e) => Math.abs(e - Ex) > (exclude || 0)));
    const branch = (side) => {
      const out = [];
      for (let i = 0; i < cols.length; i++) {
        const cand = cols[i].filter((e) => (side < 0 ? e < Ex : e > Ex));
        if (!out.length) {
          if (i >= (seedWithin || Infinity)) break;
          if (!cand.length) continue;
          out.push({ x: map.x[i], E: side < 0 ? Math.max(...cand) : Math.min(...cand), i });
          continue;
        }
        const a = out[out.length - 1], b = out[out.length - 2];
        const pred = b ? a.E + ((a.E - b.E) * (i - a.i)) / (a.i - b.i) : a.E;
        const step = b ? Math.abs(a.E - b.E) / (a.i - b.i) : 0;
        const tol = Math.max(6 * dE, 3 * step * (i - a.i), 0.01 * Math.abs(Ex));
        let best = null;
        cand.forEach((e) => { if (Math.abs(e - pred) <= tol && (best === null || Math.abs(e - pred) < Math.abs(best - pred))) best = e; });
        if (best === null) { if (i - a.i > 2) break; continue; }
        out.push({ x: map.x[i], E: best, i });
      }
      return out.map(({ x, E }) => ({ x, E, branch: side < 0 ? "LP" : "UP" }));
    };
    return [...branch(-1), ...branch(1)];
  }

  /** Nelder–Mead simplex minimiser. */
  function nelderMead(f, x0, step, iters) {
    const n = x0.length;
    let simplex = [x0.slice()];
    for (let i = 0; i < n; i++) { const v = x0.slice(); v[i] += step[i]; simplex.push(v); }
    let vals = simplex.map(f);
    for (let it = 0; it < (iters || 800); it++) {
      const order = vals.map((v, i) => i).sort((a, b) => vals[a] - vals[b]);
      simplex = order.map((i) => simplex[i]); vals = order.map((i) => vals[i]);
      const cen = new Array(n).fill(0);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) cen[j] += simplex[i][j] / n;
      const pt = (a, b, t) => a.map((v, j) => v + t * (b[j] - v));
      const xr = pt(cen, simplex[n], -1), fr = f(xr);
      if (fr < vals[0]) { const xe = pt(cen, simplex[n], -2), fe = f(xe); [simplex[n], vals[n]] = fe < fr ? [xe, fe] : [xr, fr]; }
      else if (fr < vals[n - 1]) { simplex[n] = xr; vals[n] = fr; }
      else {
        const xc = pt(cen, simplex[n], 0.5), fc = f(xc);
        if (fc < vals[n]) { simplex[n] = xc; vals[n] = fc; }
        else { for (let i = 1; i <= n; i++) { simplex[i] = pt(simplex[0], simplex[i], 0.5); vals[i] = f(simplex[i]); } }
      }
      if (Math.abs(vals[n] - vals[0]) < 1e-14 * (1 + Math.abs(vals[0]))) break;
    }
    const best = vals.indexOf(Math.min(...vals));
    return { x: simplex[best], f: vals[best] };
  }

  /**
   * Fit angle-resolved polariton points with the two-oscillator model.
   * Free parameters: E0 (cavity energy at normal incidence), n_eff, g; Ex fixed.
   * Each point is assigned to the nearer branch. n0 = index of the incident medium.
   */
  function fitPolaritons(points, Ex, guess, n0) {
    const cost = ([E0, neff, g]) => {
      if (!(E0 > 0 && neff > 0.5 && g >= 0)) return 1e9;
      let s = 0;
      for (const p of points) {
        const Ec = cavityEnergy(E0, neff, p.x, n0);
        if (!isFinite(Ec)) return 1e9;
        const { Eup, Elp } = coupledOscillator(Ec, Ex, g);
        s += Math.min((p.E - Eup) ** 2, (p.E - Elp) ** 2);
      }
      return s / points.length;
    };
    const g0 = guess || { E0: Ex, neff: 1.6, g: 0.05 };
    const res = nelderMead(cost, [g0.E0, g0.neff, g0.g], [0.05 * g0.E0, 0.2, Math.max(0.01, 0.5 * g0.g)], 3000);
    const [E0, neff, g] = res.x;
    return { E0, neff, g, rabi: 2 * g, rms: Math.sqrt(res.f) };
  }

  return {
    complex: { cx, add, sub, mul, div, conj, abs2, sqrt: csqrt, cos: ccos, sin: csin, exp: cexp, scale },
    HC_EV_NM, material, nkFromEps, epsFromNk, kzNorm, kParallel, solve, solveUnpolarised, solveAny, fieldProfile, dbr,
    angleMap, kMap, angleScan, repeatBlock, filledCavityRabi, thicknessMap, cavityEnergy, coupledOscillator, findPeaks, extractResonances, trackPolaritons, nelderMead, fitPolaritons,
  };
});
