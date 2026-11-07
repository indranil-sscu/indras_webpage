// Tests for assets/js/curious-core.js against exact results. Run: node --test tests/*.test.js
const test = require("node:test");
const assert = require("node:assert");
const C = require("../assets/js/curious-core.js");

const near = (got, exp, tol) => assert.ok(Math.abs(got - exp) <= tol, `${got} != ${exp} ± ${tol}`);

test("seeded generator is reproducible and uniform", () => {
  const a = C.rng(42), b = C.rng(42);
  for (let i = 0; i < 5; i++) assert.strictEqual(a(), b());
  const r = C.rng(7), N = 200000;
  let s = 0, s2 = 0;
  for (let i = 0; i < N; i++) { const u = r(); assert.ok(u >= 0 && u < 1); s += u; s2 += u * u; }
  near(s / N, 0.5, 0.005);
  near(s2 / N - (s / N) ** 2, 1 / 12, 0.002);
});

test("exact dice-sum distributions", () => {
  const two = C.diceSumDistribution(2, 6);
  assert.strictEqual(two.min, 2); assert.strictEqual(two.max, 12);
  near(two.p[7 - 2], 6 / 36, 1e-15);
  near(two.p[2 - 2], 1 / 36, 1e-15);
  near(two.p.reduce((s, v) => s + v, 0), 1, 1e-12);
  const three = C.diceSumDistribution(3, 6);
  near(three.p[10 - 3], 27 / 216, 1e-15);
  near(three.mean, 10.5, 0); near(three.variance, 35 / 4, 1e-12);
  // mean and variance from the distribution itself
  const d = C.diceSumDistribution(4, 20);
  const m = d.p.reduce((s, v, i) => s + v * (i + d.min), 0);
  near(m, d.mean, 1e-9);
  near(d.p.reduce((s, v, i) => s + v * (i + d.min - m) ** 2, 0), d.variance, 1e-8);
});

test("simulated dice agree with theory (law of large numbers)", () => {
  const r = C.rng(1), N = 120000, counts = new Array(11).fill(0);
  for (let i = 0; i < N; i++) { const f = C.rollDice(2, 6, r); counts[f[0] + f[1] - 2]++; }
  const th = C.diceSumDistribution(2, 6).p;
  counts.forEach((c, i) => near(c / N, th[i], 4 * Math.sqrt(th[i] * (1 - th[i]) / N)));
});

test("binomial probabilities and longest runs", () => {
  near(C.binomialPMF(10, 5, 0.5), 252 / 1024, 1e-12);
  near([...Array(21).keys()].reduce((s, k) => s + C.binomialPMF(20, k, 0.3), 0), 1, 1e-12);
  assert.strictEqual(C.longestRun([1, 0, 0, 0, 1, 1]), 3);
  assert.strictEqual(C.longestRun([]), 0);
});

test("Monty Hall: switching wins (n−1)/n, staying 1/n", () => {
  for (const doors of [3, 5, 10]) {
    const r = C.rng(doors), N = 60000;
    let stay = 0, sw = 0;
    for (let i = 0; i < N; i++) {
      const g = C.montyGame(doors, r);
      assert.strictEqual(g.opened.length, doors - 2);
      assert.ok(!g.opened.includes(g.car) && !g.opened.includes(g.pick) && g.offered !== g.pick);
      stay += g.winStay; sw += g.winSwitch;
      assert.ok(g.winStay !== g.winSwitch); // exactly one strategy wins
    }
    const th = C.montyTheory(doors);
    near(stay / N, th.stay, 4 * Math.sqrt(th.stay * (1 - th.stay) / N));
    near(sw / N, th.switch, 4 * Math.sqrt(th.switch * (1 - th.switch) / N));
  }
  near(C.montyTheory(3).switch, 2 / 3, 1e-15);
});

test("Fourier coefficients match numerical integration", () => {
  for (const wave of ["square", "sawtooth", "triangle", "pulse"]) {
    const d = 0.3, M = 8192;
    const samples = Array.from({ length: M }, (_, j) => C.WAVES[wave].f(-Math.PI + (2 * Math.PI * j) / M, d));
    const num = C.dftCoefficients(samples, 9), ex = C.fourierCoefficients(wave, 9, d);
    near(num.a0half, C.WAVES[wave].a0half(d), 2e-3);
    ex.forEach((c, i) => { near(num.coeffs[i].a, c.a, 2e-3); near(num.coeffs[i].b, c.b, 2e-3); });
  }
  near(C.fourierCoefficients("square", 3)[2].b, 4 / (3 * Math.PI), 1e-15);
});

test("partial sums converge, with the Gibbs overshoot at a jump", () => {
  // square wave: S_K(π/2) → 1 (Leibniz series); the peak → (2/π)Si(π) = 1.17898, i.e. ≈ 9 % of the jump of 2
  near(C.partialSum(C.fourierCoefficients("square", 2001), Math.PI / 2), 1, 1e-3);
  near(C.maxOfPartialSum("square", 201, 0, 40000), 1.17898, 2e-3);
  // triangle wave is continuous: no overshoot, and the power fraction approaches 1 quickly
  assert.ok(C.maxOfPartialSum("triangle", 51, 0) <= 1 + 1e-9);
  near(C.powerFraction("triangle", 1), 96 / Math.PI ** 4, 1e-3); // (8/π²)²/2 ÷ 1/3
  near(C.powerFraction("square", 1), 8 / Math.PI ** 2, 1e-3);
  assert.ok(C.powerFraction("square", 99) > 0.99);
});
