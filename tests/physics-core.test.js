// Tests for assets/js/physics-core.js against analytical results. Run: node --test tests/*.test.js
const test = require("node:test");
const assert = require("node:assert");
const P = require("../assets/js/physics-core.js");

const near = (got, exp, tol) => assert.ok(Math.abs(got - exp) <= tol, `${got} != ${exp} ± ${tol}`);
const DT = 0.002; // the default time step of the app (2 ms)

function measuredPeriod(kind, y0, p, tEnd, dt, idx) {
  const m = P.periodMeter(), s = P.simulate(kind, y0, p, tEnd, dt || DT);
  s.t.forEach((t, i) => m.push(t, s.y[i][idx || 0]));
  return m.period;
}

test("RK4 is fourth-order accurate on y' = y", () => {
  const err = (dt) => { let y = [1]; for (let i = 0; i < Math.round(1 / dt); i++) y = P.rk4((v) => [v[0]], y, dt); return Math.abs(y[0] - Math.E); };
  const ratio = err(0.02) / err(0.01);
  near(Math.log2(ratio), 4, 0.1);
});

test("small-angle period of a simple pendulum is 2π√(L/g)", () => {
  for (const [L, g] of [[1, 9.80665], [0.5, 1.625], [2, 24.79]]) {
    const T0 = P.smallAnglePeriod(L, g);
    near(measuredPeriod("simple", [0.01, 0], { L, g }, 3.2 * T0), T0, 2e-4 * T0);
  }
  near(P.smallAnglePeriod(1, 9.80665), 2.0064, 1e-4);
});

test("large-amplitude period matches the exact elliptic-integral result", () => {
  // T/T0 = 2K(sin(θ0/2))/π: 1.0174 at 30°, 1.1803 at 90°, 1.7622 at 150° (standard tables)
  near(P.exactPeriod(1, 1, Math.PI / 6) / (2 * Math.PI), 1.017408797, 1e-8);
  near(P.exactPeriod(1, 1, Math.PI / 2) / (2 * Math.PI), 1.180340599, 1e-8);
  near(P.exactPeriod(1, 1, (5 * Math.PI) / 6) / (2 * Math.PI), 1.762204, 1e-5);
  const p = { L: 1, g: 9.80665 };
  for (const a of [Math.PI / 4, Math.PI / 2, 2.6]) {
    const T = P.exactPeriod(p.L, p.g, a);
    near(measuredPeriod("simple", [a, 0], p, 3.2 * T), T, 2e-4 * T);
  }
});

test("simple pendulum conserves energy without damping and loses it with damping", () => {
  const p = { L: 1, g: 9.80665 }, s = P.simulate("simple", [2.5, 0], p, 30, DT), E0 = s.E[0];
  s.E.forEach((E) => near(E, E0, 1e-6 * E0));
  const pd = { L: 1, g: 9.80665, gamma: 0.3 };
  const sd = P.simulate("simple", [0.02, 0], pd, 20, DT);
  // small oscillations decay as exp(−γt/2), so the energy falls as exp(−γt)
  near(sd.E[sd.E.length - 1] / sd.E[0], Math.exp(-0.3 * 20), 0.15 * Math.exp(-0.3 * 20));
  const Td = P.dampedPeriod(1, 9.80665, 0.3);
  near(measuredPeriod("simple", [0.02, 0], pd, 8, DT), Td, 2e-3 * Td);
});

test("double pendulum: energy drift below 1e-3 over 30 s at the default step", () => {
  const p = { m1: 1, m2: 1, L1: 1, L2: 1, g: 9.80665 };
  for (const y0 of [[2.0, 2.5, 0, 0], [(3 * Math.PI) / 4, -0.3, 0, 0], [3.0, 3.1, 0, 0], [1.2, -1.6, 2, -3]]) {
    const s = P.simulate("double", y0, p, 30, DT), E0 = s.E[0];
    const drift = Math.max(...s.E.map((E) => Math.abs(E - E0) / E0));
    assert.ok(drift < 1e-3, `drift ${drift} for ${y0}`);
  }
  // heavier lower bob and unequal rods
  const q = { m1: 1, m2: 3, L1: 1, L2: 0.6, g: 9.80665 }, s = P.simulate("double", [2.4, 1.0, 0, 0], q, 30, DT);
  assert.ok(Math.max(...s.E.map((E) => Math.abs(E - s.E[0]) / s.E[0])) < 1e-3);
});

test("double pendulum: small-oscillation normal modes", () => {
  // m1 = m2, L1 = L2: θ2 = ±√2 θ1, ω² = (2 ∓ √2) g/L
  const g = 9.80665, p = { m1: 1, m2: 1, L1: 1, L2: 1, g };
  const [w1, w2] = P.doubleNormalModes(1, 1, 1, g);
  near(w1, Math.sqrt((2 - Math.SQRT2) * g), 1e-12);
  near(w2, Math.sqrt((2 + Math.SQRT2) * g), 1e-12);
  const a = 0.002;
  near(measuredPeriod("double", [a, Math.SQRT2 * a, 0, 0], p, 3.2 * (2 * Math.PI) / w1, DT, 0), (2 * Math.PI) / w1, 5e-4 * (2 * Math.PI) / w1);
  near(measuredPeriod("double", [a, -Math.SQRT2 * a, 0, 0], p, 3.2 * (2 * Math.PI) / w2, DT, 0), (2 * Math.PI) / w2, 5e-4 * (2 * Math.PI) / w2);
});

test("double pendulum: geometry and energy zero point", () => {
  const p = { m1: 2, m2: 1, L1: 1.5, L2: 1, g: 10 };
  const pos = P.doublePositions([Math.PI / 2, 0, 0, 0], p);
  near(pos.x1, 1.5, 1e-15); near(pos.y1, 0, 1e-15); near(pos.x2, 1.5, 1e-15); near(pos.y2, -1, 1e-15);
  assert.strictEqual(P.doubleEnergy([0, 0, 0, 0], p), 0);
  near(P.doubleEnergy([Math.PI, Math.PI, 0, 0], p), 2 * (3 * 10 * 1.5) + 2 * 10 * 1, 1e-12); // both upright
});
