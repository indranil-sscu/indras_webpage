// Unit tests for assets/js/sci-tools.js — run with: node --test tests/*.test.js
const test = require("node:test");
const assert = require("node:assert");
const T = require("../assets/js/sci-tools.js");

const close = (got, exp, rel = 1e-9) =>
  assert.ok(Math.abs(got - exp) <= rel * Math.abs(exp), `${got} != ${exp}`);

test("energy units against CODATA 2018 reference values", () => {
  const eV = T.convertEnergy(1, "eV");
  close(eV.cm1, 8065.543937, 1e-9);
  close(eV.nm, 1239.841984, 1e-9);
  close(eV.THz, 241.7989242, 1e-9);
  close(eV.kJmol, 96.48533212, 1e-9);
  close(eV.kcalmol, 23.06054783, 1e-9);
  close(eV.K, 11604.51812, 1e-9);
  const Eh = T.convertEnergy(1, "Eh");
  close(Eh.eV, 27.211386245988, 1e-12);
  close(Eh.cm1, 219474.6313632, 1e-12);
  close(Eh.kcalmol, 627.5094740631, 1e-12);
  close(T.convertEnergy(1, "Ry").eV, 13.605693122994, 1e-12);
  close(T.convertEnergy(800, "nm").cm1, 12500, 1e-12);
});

test("every unit round-trips", () => {
  for (const u of Object.keys(T.ENERGY_UNITS)) {
    close(T.convertEnergy(T.convertEnergy(2.5, "eV")[u], u).eV, 2.5, 1e-12);
  }
});

test("spectral bandwidth is exact and invertible", () => {
  const bw = T.bandwidthFromNm(800, 10);
  close(bw.dnu_THz, (299792458 / 795e-9 - 299792458 / 805e-9) / 1e12, 1e-12);
  for (const [l, d] of [[800, 10], [400, 80], [1030, 3], [1550, 0.01]]) {
    close(T.bandwidthToNm(l, T.bandwidthFromNm(l, d).dnu_THz), d, 1e-12);
  }
});

test("transform-limited pulses", () => {
  close(T.TBP.gaussian, 0.441271, 1e-5);
  close(T.TBP.sech2, 0.314833, 1e-5);
  close(T.transformLimitedDuration(800, 10, "gaussian"), 94.199, 1e-4);
  close(T.transformLimitedBandwidth(800, 94.19934654, "gaussian"), 10, 1e-8);
});

test("delay line: 1 mm on a retro-reflector is 6.671 ps", () => {
  close(T.delayFromDistance(1, 2, 1), 6671.281904, 1e-9);
  close(T.distanceFromDelay(6671.281904, 2, 1), 1, 1e-9);
});

test("pulse energy, peak power, fluence", () => {
  const p = T.pulseParameters({ avgPower_W: 1, repRate_Hz: 1000, duration_fs: 100, diameter_um: 100, wavelength_nm: 800, shape: "gaussian" });
  close(p.energy_J, 1e-3);
  close(p.peakPower_W, 9.394372786e9, 1e-9);
  close(p.peakFluence_Jcm2, 25.46479089, 1e-9);
  close(p.photonsPerPulse, 4.027293254e15, 1e-9);
});

test("GDD broadening of a Gaussian pulse", () => {
  close(T.gddBroadening(10, 100), 29.47413819, 1e-9);
  close(T.gddBroadening(50, 0), 50, 1e-12);
});
