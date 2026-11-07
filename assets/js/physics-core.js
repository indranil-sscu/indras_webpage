/*
 * physics-core.js: the physics behind the pendulum playground in "Science for the Curious".
 * Simple (optionally damped) and double pendulums integrated with the classical
 * fourth-order Runge–Kutta method, their energies, the exact large-amplitude period
 * and a period meter. Pure functions, no DOM; unit-tested in Node (tests/physics-core.test.js)
 * against analytical results.
 *
 * Angles are measured from the downward vertical, in radians; SI units throughout.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.Physics = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /** One classical RK4 step for y' = f(y) (autonomous). y is an array; returns a new array. */
  function rk4(f, y, dt) {
    const n = y.length, k1 = f(y), t = new Array(n);
    for (let i = 0; i < n; i++) t[i] = y[i] + 0.5 * dt * k1[i];
    const k2 = f(t);
    for (let i = 0; i < n; i++) t[i] = y[i] + 0.5 * dt * k2[i];
    const k3 = f(t);
    for (let i = 0; i < n; i++) t[i] = y[i] + dt * k3[i];
    const k4 = f(t), out = new Array(n);
    for (let i = 0; i < n; i++) out[i] = y[i] + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    return out;
  }

  /* ---------------- simple pendulum ---------------- */
  /** θ'' = −(g/L) sin θ − γ θ'. State [θ, ω]; p = {L, g, gamma}. */
  const simpleDeriv = (p) => (y) => [y[1], -(p.g / p.L) * Math.sin(y[0]) - (p.gamma || 0) * y[1]];
  const simpleStep = (y, p, dt) => rk4(simpleDeriv(p), y, dt);
  /** Energy per unit mass (J/kg), zero at the bottom: ½L²ω² + gL(1 − cos θ). */
  const simpleEnergy = (y, p) => 0.5 * p.L * p.L * y[1] * y[1] + p.g * p.L * (1 - Math.cos(y[0]));
  /** Small-angle period 2π√(L/g). */
  const smallAnglePeriod = (L, g) => 2 * Math.PI * Math.sqrt(L / g);
  /** Arithmetic–geometric mean. */
  function agm(a, b) {
    for (let i = 0; i < 60 && Math.abs(a - b) > 1e-15 * a; i++) [a, b] = [(a + b) / 2, Math.sqrt(a * b)];
    return (a + b) / 2;
  }
  /** Exact period of an undamped pendulum released from rest at amplitude θ0:
   *  T = 4√(L/g) K(sin(θ0/2)) = 2π√(L/g) / AGM(1, cos(θ0/2)). */
  const exactPeriod = (L, g, theta0) => smallAnglePeriod(L, g) / agm(1, Math.cos(Math.abs(theta0) / 2));
  /** Small-angle period with linear damping γ: 2π / √(g/L − γ²/4) (Infinity if over-damped). */
  function dampedPeriod(L, g, gamma) {
    const w2 = g / L - (gamma * gamma) / 4;
    return w2 > 0 ? (2 * Math.PI) / Math.sqrt(w2) : Infinity;
  }

  /* ---------------- double pendulum ---------------- */
  /**
   * Two point masses on massless rods. State [θ1, θ2, ω1, ω2]; p = {m1, m2, L1, L2, g}.
   * Equations of motion from the Lagrangian (e.g. Shinbrot et al. 1992).
   */
  const doubleDeriv = (p) => (y) => {
    const { m1, m2, L1, L2, g } = p, [a1, a2, w1, w2] = y, d = a1 - a2, sd = Math.sin(d), cd = Math.cos(d);
    const den = 2 * m1 + m2 - m2 * Math.cos(2 * d);
    const al1 = (-g * (2 * m1 + m2) * Math.sin(a1) - m2 * g * Math.sin(a1 - 2 * a2) - 2 * sd * m2 * (w2 * w2 * L2 + w1 * w1 * L1 * cd)) / (L1 * den);
    const al2 = (2 * sd * (w1 * w1 * L1 * (m1 + m2) + g * (m1 + m2) * Math.cos(a1) + w2 * w2 * L2 * m2 * cd)) / (L2 * den);
    return [w1, w2, al1, al2];
  };
  const doubleStep = (y, p, dt) => rk4(doubleDeriv(p), y, dt);
  /** Total energy (J), zero with both rods hanging straight down. */
  function doubleEnergy(y, p) {
    const { m1, m2, L1, L2, g } = p, [a1, a2, w1, w2] = y;
    const T = 0.5 * m1 * L1 * L1 * w1 * w1 + 0.5 * m2 * (L1 * L1 * w1 * w1 + L2 * L2 * w2 * w2 + 2 * L1 * L2 * w1 * w2 * Math.cos(a1 - a2));
    const V = (m1 + m2) * g * L1 * (1 - Math.cos(a1)) + m2 * g * L2 * (1 - Math.cos(a2));
    return T + V;
  }
  /** Cartesian positions of the bobs, pivot at the origin, y up. */
  function doublePositions(y, p) {
    const x1 = p.L1 * Math.sin(y[0]), y1 = -p.L1 * Math.cos(y[0]);
    return { x1, y1, x2: x1 + p.L2 * Math.sin(y[1]), y2: y1 - p.L2 * Math.cos(y[1]) };
  }
  /** Small-oscillation normal-mode angular frequencies of a double pendulum with L1 = L2 = L,
   *  μ = m2/m1: ω² = (g/L)(1 + μ ∓ √(μ(1 + μ))). Returns [slow, fast]. */
  function doubleNormalModes(m1, m2, L, g) {
    const mu = m2 / m1, s = Math.sqrt(mu * (1 + mu));
    return [Math.sqrt((g / L) * (1 + mu - s)), Math.sqrt((g / L) * (1 + mu + s))];
  }

  /* ---------------- measuring a period ---------------- */
  /**
   * A period meter: feed it successive samples (t, x); it detects upward zero crossings of x
   * (with linear interpolation between samples) and reports the time between the last two.
   */
  function periodMeter() {
    let tp = null, xp = null, lastCross = null, period = NaN, count = 0;
    return {
      push(t, x) {
        if (xp !== null && xp < 0 && x >= 0) {
          const tc = tp + ((t - tp) * -xp) / (x - xp);
          if (lastCross !== null) { period = tc - lastCross; count++; }
          lastCross = tc;
        }
        tp = t; xp = x;
        return period;
      },
      get period() { return period; },
      get count() { return count; },
      reset() { tp = xp = lastCross = null; period = NaN; count = 0; },
    };
  }
  /** Integrate a pendulum for tEnd seconds with step dt. kind "simple" | "double". Returns {t[], y[][], E[]}. */
  function simulate(kind, y0, p, tEnd, dt) {
    const step = kind === "double" ? doubleStep : simpleStep, energy = kind === "double" ? doubleEnergy : simpleEnergy;
    const n = Math.round(tEnd / dt), t = [0], ys = [y0.slice()], E = [energy(y0, p)];
    let y = y0.slice();
    for (let i = 1; i <= n; i++) { y = step(y, p, dt); t.push(i * dt); ys.push(y); E.push(energy(y, p)); }
    return { t, y: ys, E };
  }
  /** Gravitational acceleration at the surface (m/s²), standard values. */
  const GRAVITY = { earth: 9.80665, moon: 1.625, mars: 3.721, jupiter: 24.79 };

  return {
    rk4, simpleDeriv, simpleStep, simpleEnergy, smallAnglePeriod, agm, exactPeriod, dampedPeriod,
    doubleDeriv, doubleStep, doubleEnergy, doublePositions, doubleNormalModes, periodMeter, simulate, GRAVITY,
  };
});
