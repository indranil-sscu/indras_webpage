// Unit tests for assets/js/tmm-core.js — run with: node --test tests/*.test.js
const test = require("node:test");
const assert = require("node:assert");
const T = require("../assets/js/tmm-core.js");
const M = T.material;

const near = (got, exp, tol) => assert.ok(Math.abs(got - exp) <= tol, `${got} != ${exp} ± ${tol}`);
const air = M.constant(1), glass = M.constant(1.5);

test("Fresnel: normal incidence, Brewster angle, total internal reflection", () => {
  const st = { incident: air, layers: [], substrate: glass };
  near(T.solve(st, 600, 0, "s").R, 0.04, 1e-12);
  near(T.solve(st, 600, 0, "p").T, 0.96, 1e-12);
  near(T.solve(st, 600, (Math.atan(1.5) * 180) / Math.PI, "p").R, 0, 1e-15);
  const tir = { incident: glass, layers: [], substrate: air };
  near(T.solve(tir, 600, 60, "s").R, 1, 1e-12);
  near(T.solve(tir, 600, 60, "p").R, 1, 1e-12);
});

test("quarter-wave anti-reflection coating", () => {
  const n = Math.sqrt(1.5);
  const st = { incident: air, layers: [{ material: M.constant(n), d_nm: 600 / (4 * n) }], substrate: glass };
  near(T.solve(st, 600, 0, "s").R, 0, 1e-15);
});

test("Bragg mirror peak reflectance matches the closed form", () => {
  const nH = 2.3, nL = 1.46, ns = 1.52, N = 6;
  const st = { incident: air, layers: T.dbr(nH, nL, 600, N), substrate: M.constant(ns) };
  const Rth = ((nL ** (2 * N) - ns * nH ** (2 * N)) / (nL ** (2 * N) + ns * nH ** (2 * N))) ** 2;
  near(T.solve(st, 600, 0, "s").R, Rth, 1e-12);
});

test("Fabry–Pérot slab follows the Airy formula", () => {
  const st = { incident: air, layers: [{ material: M.constant(2), d_nm: 300 }], substrate: air };
  near(T.solve(st, 600, 0, "s").T, 1, 1e-12); // λ = 2nd/m
  near(T.solve(st, 480, 0, "s").T, 1 / (1 + 0.5625), 1e-12); // anti-resonance, F = 4R/(1−R)²
});

test("energy conservation R + T + A = 1 with absorption, both polarisations", () => {
  const st = { incident: air, layers: [{ material: M.constant(2.0, 0.3), d_nm: 80 }, { material: M.drude(4, 9, 0.05), d_nm: 20 }], substrate: glass };
  for (const pol of ["s", "p"]) for (const th of [0, 30, 60, 80]) {
    const r = T.solve(st, 550, th, pol);
    near(r.R + r.T + r.A, 1, 1e-10);
    assert.ok(r.A > 0 && r.R >= 0 && r.T >= 0);
  }
});

test("Lorentz and Drude media are absorbing (k ≥ 0)", () => {
  const L = M.lorentz(2.25, [{ nu0: 16000, gamma: 200, f: 0.05 }]);
  assert.ok(L.nk(625)[1] > 0.1);
  assert.ok(M.drude(1, 9, 0.02).nk(600)[1] > 1); // metallic
});

test("field profile: antinode at the cavity centre; absorption equals ∫(4πnk/λ)|E|²", () => {
  const lam0 = T.HC_EV_NM / 2.1, nb = 1.5;
  const layers = [...T.dbr(2.3, 1.46, lam0, 6).reverse(), { material: M.constant(nb), d_nm: lam0 / (2 * nb) }, ...T.dbr(2.3, 1.46, lam0, 6)];
  const fp = T.fieldProfile({ incident: air, layers, substrate: air }, lam0, 0, "s", 0.5, 50);
  let im = 0; fp.I2.forEach((v, i) => { if (v > fp.I2[im]) im = i; });
  const centre = layers.slice(0, 12).reduce((s, L) => s + L.d_nm, 0) + lam0 / (4 * nb);
  near(fp.z[im], centre, 0.5);
  near(fp.I2[fp.I2.length - 1], 1, 1e-6); // T = 1 on resonance, n0 = ns

  const st = { incident: air, layers: [{ material: glass, d_nm: 100 }, { material: M.constant(1.8, 0.2), d_nm: 60 }, { material: glass, d_nm: 100 }], substrate: glass };
  const f = T.fieldProfile(st, 550, 0, "s", 0.05, 0);
  let A = 0;
  for (let i = 1; i < f.z.length; i++) if (f.layerIndex[i] === 1 && f.layerIndex[i - 1] === 1) A += 0.5 * (f.I2[i] + f.I2[i - 1]) * (f.z[i] - f.z[i - 1]);
  near(A * ((4 * Math.PI * 1.8 * 0.2) / 550), T.solve(st, 550, 0, "s").A, 1e-5);
});

test("strong coupling: Rabi splitting scales as √f", () => {
  const Ex = 2.1, lam0 = T.HC_EV_NM / Ex, nb = 1.5;
  const split = (f) => {
    const act = M.lorentz(nb * nb, [{ nu0: 1e7 / lam0, gamma: (0.002 / Ex) * (1e7 / lam0), f }]);
    const st = { incident: air, layers: [...T.dbr(2.3, 1.46, lam0, 6).reverse(), { material: act, d_nm: lam0 / (2 * nb) }, ...T.dbr(2.3, 1.46, lam0, 6)], substrate: air };
    const E = [], t = [];
    for (let e = 1.9; e <= 2.3; e += 0.0001) { E.push(e); t.push(T.solve(st, T.HC_EV_NM / e, 0, "s").T); }
    const pk = T.findPeaks(t, "max", 1e-3).map((i) => E[i]);
    return Math.min(...pk.filter((e) => e > Ex)) - Math.max(...pk.filter((e) => e < Ex));
  };
  near(split(0.02) / split(0.005), 2, 0.02);
});

test("coupled-oscillator model and Hopfield coefficients", () => {
  const r = T.coupledOscillator(2.1, 2.1, 0.05);
  near(r.Eup - r.Elp, 0.1, 1e-12);
  near(r.LP_matter, 0.5, 1e-12);
  const far = T.coupledOscillator(1.5, 2.1, 0.05); // cavity far below: LP is photon-like
  assert.ok(far.LP_photon > 0.99);
});

test("fit recovers synthetic polariton dispersion and a TMM anticrossing", () => {
  const syn = [];
  for (let a = 0; a <= 40; a += 2) { const o = T.coupledOscillator(T.cavityEnergy(2.05, 1.7, a), 2.1, 0.06); syn.push({ x: a, E: o.Eup }, { x: a, E: o.Elp }); }
  const f = T.fitPolaritons(syn, 2.1, { E0: 2.0, neff: 1.5, g: 0.03 });
  near(f.E0, 2.05, 1e-4); near(f.neff, 1.7, 1e-3); near(f.g, 0.06, 1e-4);

  const Ex = 2.1, lam0 = T.HC_EV_NM / Ex, nb = 1.5;
  const act = M.lorentz(nb * nb, [{ nu0: 1e7 / lam0, gamma: (0.002 / Ex) * (1e7 / lam0), f: 0.02 }]);
  const st = { incident: air, layers: [...T.dbr(2.3, 1.46, lam0, 6).reverse(), { material: act, d_nm: lam0 / (2 * nb) }, ...T.dbr(2.3, 1.46, lam0, 6)], substrate: air };
  const th = []; for (let a = 0; a <= 40; a += 2) th.push(a);
  const E = []; for (let e = 1.95; e <= 2.3; e += 0.0004) E.push(e);
  const fit = T.fitPolaritons(T.extractResonances(T.angleMap(st, "s", "T", th, E), "T", 0.002, 2), Ex, { E0: 2.08, neff: 1.6, g: 0.04 });
  near(fit.E0, 2.1, 0.01); // resonant cavity
  near(fit.rabi, 0.1193, 0.003); // equals the splitting measured directly at normal incidence
  assert.ok(fit.rms < 0.002);
});

test("polariton fit with light incident from a dense medium (n0 = 1.4, e.g. CaF2 window)", () => {
  const syn = [];
  for (let a = 0; a <= 50; a += 2) { const o = T.coupledOscillator(T.cavityEnergy(1720, 1.45, a, 1.4), 1730, 60); syn.push({ x: a, E: o.Eup }, { x: a, E: o.Elp }); }
  const f = T.fitPolaritons(syn, 1730, { E0: 1700, neff: 2, g: 30 }, 1.4);
  near(f.E0, 1720, 0.01); near(f.neff, 1.45, 1e-4); near(f.g, 60, 0.01);
});

test("branch tracking ignores other cavity orders (mid-IR Fabry–Pérot, VSC-like)", () => {
  const cm = (k) => k / 8065.543937, Ex = cm(1730);
  const au = M.drude(9.5, 9.0, 0.07);
  const st = { incident: M.constant(1.4), layers: [{ material: au, d_nm: 10 }, { material: M.lorentz(2.0, [{ nu0: 1730, gamma: 15, f: 0.025 }]), d_nm: 6060 }, { material: au, d_nm: 10 }], substrate: M.constant(1.43) };
  const th = []; for (let a = 0; a <= 60; a += 2) th.push(a);
  const E = []; for (let k = 1400; k <= 2100; k += 1) E.push(cm(k));
  const pts = T.trackPolaritons(T.angleMap(st, "p", "T", th, E), "T", 0.01, Ex, cm(22.5));
  const fit = T.fitPolaritons(pts, Ex, { E0: cm(1729), neff: 2.2, g: cm(97) }, 1.4);
  near(fit.rabi / cm(1), 195, 3);
  near(fit.E0 / cm(1), 1729, 3);
  near(fit.neff, Math.SQRT2, 0.03);
});
