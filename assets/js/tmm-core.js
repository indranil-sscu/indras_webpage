/*
 * tmm-core.js — transfer-matrix method (TMM) for planar, isotropic multilayers,
 * plus coupled-oscillator tools for polariton analysis. Used by /tools/tmm/.
 * Pure functions, no DOM access; unit-tested in Node (tests/tmm-core.test.js).
 *
 * Formalism: characteristic (Abelès) matrices in the admittance form of
 * H. A. Macleod, "Thin-Film Optical Filters" (n − ik convention internally):
 *   M_j = [[cos δ_j, i sin δ_j / η_j], [i η_j sin δ_j, cos δ_j]],
 *   δ_j = 2π N_j d_j cos θ_j / λ,  η_j = N_j cos θ_j (s),  N_j / cos θ_j (p)
 *   [B, C]ᵀ = Π M_j · [1, η_sub]ᵀ
 *   r = (η₀B − C)/(η₀B + C),  R = |r|²,
 *   T = 4η₀ Re(η_sub) / |η₀B + C|²,  A = 4η₀ Re(BC* − η_sub) / |η₀B + C|².
 * Materials are specified with the usual physics convention ñ = n + ik
 * (k ≥ 0 absorbing) and converted internally.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.TMM = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /* ---------------- complex arithmetic: [re, im] ---------------- */
  const cx = (re, im) => [re, im || 0];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1]];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1]];
  const mul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
  const div = (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; };
  const conj = (a) => [a[0], -a[1]];
  const abs2 = (a) => a[0] * a[0] + a[1] * a[1];
  const csqrt = (a) => { const r = Math.hypot(a[0], a[1]); const re = Math.sqrt((r + a[0]) / 2); const im = Math.sqrt(Math.max(0, (r - a[0]) / 2)); return [re, a[1] < 0 ? -im : im]; };
  const ccos = (a) => [Math.cos(a[0]) * Math.cosh(a[1]), -Math.sin(a[0]) * Math.sinh(a[1])];
  const csin = (a) => [Math.sin(a[0]) * Math.cosh(a[1]), Math.cos(a[0]) * Math.sinh(a[1])];
  const scale = (a, s) => [a[0] * s, a[1] * s];
  const I = [0, 1];

  /* ---------------- materials ----------------
   * Each material is { name, nk(lambda_nm) -> [n, k] } with k ≥ 0 for absorption.
   */
  const HC_EV_NM = 1239.841984; // E[eV] = 1239.84 / λ[nm]

  /** Refractive index from a complex permittivity ε (Im ε ≥ 0). */
  function nkFromEps(eps) {
    const s = csqrt(eps);
    return s[1] < 0 ? [-s[0], -s[1]] : s;
  }

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
    /** Drude metal: ε(E) = ε∞ − Ep² / (E² + iγE), energies in eV. */
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
   * @param {object} stack {incident: material, layers: [{material, d_nm}], substrate: material}
   * @param {number} lambda_nm vacuum wavelength
   * @param {number} theta_deg angle of incidence in the incident medium
   * @param {"s"|"p"} pol polarisation
   * @returns {R, T, A, r}
   */
  function solve(stack, lambda_nm, theta_deg, pol) {
    const th0 = (theta_deg * Math.PI) / 180;
    const nk0 = stack.incident.nk(lambda_nm);
    if (nk0[1] !== 0) throw new Error("The incident medium must be transparent.");
    const n0 = nk0[0];
    const sin0 = n0 * Math.sin(th0); // invariant n·sinθ
    const N = (nk) => [nk[0], -nk[1]]; // physics n+ik  ->  Macleod n−ik
    const admittance = (Nj) => {
      // cosθ_j from Snell's law; choose the decaying branch (Im(N cosθ) ≤ 0)
      let ct = csqrt(sub(cx(1, 0), div(cx(sin0 * sin0, 0), mul(Nj, Nj))));
      let Nc = mul(Nj, ct);
      if (Nc[1] > 0 || (Nc[1] === 0 && Nc[0] < 0)) { ct = scale(ct, -1); Nc = scale(Nc, -1); }
      const eta = pol === "p" ? div(Nj, ct) : Nc;
      return { ct, Nc, eta };
    };
    const inc = admittance(N(nk0));
    const subA = admittance(N(stack.substrate.nk(lambda_nm)));
    // M = Π M_j
    let m11 = cx(1, 0), m12 = cx(0, 0), m21 = cx(0, 0), m22 = cx(1, 0);
    for (const L of stack.layers) {
      const a = admittance(N(L.material.nk(lambda_nm)));
      const delta = scale(a.Nc, (2 * Math.PI * L.d_nm) / lambda_nm);
      const c = ccos(delta), s = csin(delta);
      const a11 = c, a12 = div(mul(I, s), a.eta), a21 = mul(mul(I, a.eta), s), a22 = c;
      [m11, m12, m21, m22] = [add(mul(m11, a11), mul(m12, a21)), add(mul(m11, a12), mul(m12, a22)),
                              add(mul(m21, a11), mul(m22, a21)), add(mul(m21, a12), mul(m22, a22))];
    }
    const B = add(m11, mul(m12, subA.eta));
    const C = add(m21, mul(m22, subA.eta));
    const eta0 = inc.eta;
    const den = add(mul(eta0, B), C);
    const r = div(sub(mul(eta0, B), C), den);
    const d2 = abs2(den);
    const R = abs2(r);
    const T = (4 * eta0[0] * subA.eta[0]) / d2;
    const A = (4 * eta0[0] * (mul(B, conj(C))[0] - subA.eta[0])) / d2;
    return { R, T, A, r };
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
   */
  function fieldProfile(stack, lambda_nm, theta_deg, pol, dz_nm, tail_nm) {
    const th0 = (theta_deg * Math.PI) / 180;
    const n0 = stack.incident.nk(lambda_nm)[0];
    const sin0 = n0 * Math.sin(th0);
    const N = (nk) => [nk[0], -nk[1]];
    const adm = (Nj) => {
      let ct = csqrt(sub(cx(1, 0), div(cx(sin0 * sin0, 0), mul(Nj, Nj))));
      let Nc = mul(Nj, ct);
      if (Nc[1] > 0 || (Nc[1] === 0 && Nc[0] < 0)) { ct = scale(ct, -1); Nc = scale(Nc, -1); }
      return { Nc, eta: pol === "p" ? div(Nj, ct) : Nc };
    };
    const subA = adm(N(stack.substrate.nk(lambda_nm)));
    const layers = stack.layers.map((L) => ({ ...adm(N(L.material.nk(lambda_nm))), d: L.d_nm }));
    // walk back from the substrate: fields at the back of each layer
    const back = new Array(layers.length);
    let EH = [cx(1, 0), subA.eta];
    for (let j = layers.length - 1; j >= 0; j--) {
      back[j] = EH;
      EH = propagate(layers[j], layers[j].d, EH, lambda_nm);
    }
    const eta0 = adm(N(stack.incident.nk(lambda_nm))).eta;
    const Einc = div(add(mul(eta0, EH[0]), EH[1]), scale(eta0, 2));
    const norm = abs2(Einc);
    const z = [], I2 = [], layerIndex = [];
    let z0 = 0;
    layers.forEach((L, j) => {
      const n = Math.max(1, Math.ceil(L.d / dz_nm));
      for (let i = 0; i <= n; i++) {
        const depth = (L.d * i) / n; // from the front of the layer
        const f = propagate(L, L.d - depth, back[j], lambda_nm);
        z.push(z0 + depth); I2.push(abs2(f[0]) / norm); layerIndex.push(j);
      }
      z0 += L.d;
    });
    // substrate: transmitted wave decays (or not) away from the last interface
    const nt = Math.max(1, Math.ceil((tail_nm || 0) / dz_nm));
    for (let i = 1; i <= nt; i++) {
      const dz = ((tail_nm || 0) * i) / nt;
      const ph = scale(subA.Nc, (-2 * Math.PI * dz) / lambda_nm); // e^{-i 2π N cosθ z/λ} (n−ik convention)
      const decay = Math.exp(-2 * ph[1]); // |e^{i·ph}|² with Im(ph) ≥ 0
      z.push(z0 + dz); I2.push(decay / norm); layerIndex.push(layers.length);
    }
    return { z, I2, layerIndex, total_nm: z0 };
  }
  function propagate(L, dist, EH, lambda_nm) {
    // [E, H] at distance `dist` in front of a point where the fields are EH
    const delta = scale(L.Nc, (2 * Math.PI * dist) / lambda_nm);
    const c = ccos(delta), s = csin(delta);
    return [add(mul(c, EH[0]), mul(div(mul(I, s), L.eta), EH[1])), add(mul(mul(I, L.eta), mul(s, EH[0])), mul(c, EH[1]))];
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
   */
  function trackPolaritons(map, quantity, minProminence, Ex, exclude) {
    const kind = quantity === "R" ? "min" : "max";
    const dE = Math.abs(map.E[1] - map.E[0]) || 1e-6;
    const cols = map.x.map((x, i) => findPeaks(map.Z[i], kind, minProminence).map((k) => map.E[k]).filter((e) => Math.abs(e - Ex) > (exclude || 0)));
    const branch = (side) => {
      const out = [];
      for (let i = 0; i < cols.length; i++) {
        const cand = cols[i].filter((e) => (side < 0 ? e < Ex : e > Ex));
        if (!out.length) {
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
    complex: { cx, add, sub, mul, div, abs2, csqrt },
    HC_EV_NM, material, nkFromEps, solve, solveUnpolarised, solveAny, fieldProfile, dbr,
    angleMap, thicknessMap, cavityEnergy, coupledOscillator, findPeaks, extractResonances, trackPolaritons, nelderMead, fitPolaritons,
  };
});
