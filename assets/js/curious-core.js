/*
 * curious-core.js: the mathematics behind "Science for the Curious"
 * (probability playground, Monty Hall, Fourier series). Pure functions, no DOM;
 * unit-tested in Node (tests/curious-core.test.js) against exact results.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.Curious = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /* ---------------- random numbers ---------------- */
  /** Small, fast, seedable PRNG (mulberry32); returns a function giving uniform numbers in [0, 1). */
  function rng(seed) {
    let a = (seed === undefined ? (Math.random() * 2 ** 32) >>> 0 : seed >>> 0) || 1;
    return function () {
      a = (a + 0x6d2b79f5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  const randInt = (r, n) => Math.floor(r() * n); // 0 … n−1

  /* ---------------- dice and coins ---------------- */
  /**
   * Exact distribution of the sum of `n` fair dice with `sides` faces,
   * by repeated convolution. Returns {min, max, p: [P(min) … P(max)], mean, variance}.
   */
  function diceSumDistribution(n, sides) {
    let p = [1];
    for (let d = 0; d < n; d++) {
      const q = new Array(p.length + sides - 1).fill(0);
      p.forEach((v, i) => { for (let f = 0; f < sides; f++) q[i + f] += v / sides; });
      p = q;
    }
    const mean = (n * (sides + 1)) / 2, variance = (n * (sides * sides - 1)) / 12;
    return { min: n, max: n * sides, p, mean, variance };
  }
  /** Roll `n` dice; returns the faces (1 … sides). */
  const rollDice = (n, sides, r) => Array.from({ length: n }, () => 1 + randInt(r, sides));
  /** Binomial probability C(n, k) p^k (1 − p)^(n − k), computed in log space. */
  function binomialPMF(n, k, p) {
    if (k < 0 || k > n) return 0;
    if (p === 0) return k === 0 ? 1 : 0;
    if (p === 1) return k === n ? 1 : 0;
    let lc = 0;
    for (let i = 1; i <= k; i++) lc += Math.log(n - k + i) - Math.log(i);
    return Math.exp(lc + k * Math.log(p) + (n - k) * Math.log(1 - p));
  }
  /** Longest run of identical outcomes in an array. */
  function longestRun(a) {
    let best = 0, cur = 0;
    a.forEach((v, i) => { cur = i && v === a[i - 1] ? cur + 1 : 1; best = Math.max(best, cur); });
    return best;
  }
  /** Expected longest run of heads-or-tails in N fair tosses ≈ log2(N) (Schilling 1990); exact value from a recursion for small N. */
  const expectedLongestRunApprox = (N) => (N > 0 ? Math.log2(N) : 0);

  /* ---------------- Monty Hall ---------------- */
  /**
   * One game with `doors` doors (≥ 3). The car is placed at random, the player picks a door,
   * the host (who knows where the car is) opens doors − 2 goat doors, leaving one other door closed.
   * Returns {car, pick, offered, opened[], winStay, winSwitch}.
   */
  function montyGame(doors, r, pick) {
    const car = randInt(r, doors);
    if (pick === undefined) pick = randInt(r, doors);
    // the door the host leaves closed: the car if the player missed it, otherwise a random goat door
    let offered;
    if (pick !== car) offered = car;
    else { offered = randInt(r, doors - 1); if (offered >= pick) offered++; }
    const opened = [];
    for (let d = 0; d < doors; d++) if (d !== pick && d !== offered) opened.push(d);
    return { car, pick, offered, opened, winStay: pick === car, winSwitch: offered === car };
  }
  /** Exact win probabilities for the generalised game above. */
  const montyTheory = (doors) => ({ stay: 1 / doors, switch: (doors - 1) / doors });

  /* ---------------- Fourier series ---------------- */
  /**
   * Waveforms of period 2π. Each returns its value f(x) and its real Fourier
   * coefficients a_k, b_k (k ≥ 1) and the mean a0/2:
   *   f(x) = a0/2 + Σ_k [a_k cos(kx) + b_k sin(kx)]
   */
  const wrap = (x) => { const t = (x + Math.PI) % (2 * Math.PI); return (t < 0 ? t + 2 * Math.PI : t) - Math.PI; }; // to (−π, π]
  const WAVES = {
    square: {
      name: "Square wave",
      f: (x) => { const t = wrap(x); return t === 0 || Math.abs(t) === Math.PI ? 0 : t > 0 ? 1 : -1; },
      a0half: () => 0,
      coef: (k) => ({ a: 0, b: k % 2 ? 4 / (Math.PI * k) : 0 }),
      formula: "b_k = 4/(πk) for odd k",
    },
    sawtooth: {
      name: "Sawtooth wave",
      f: (x) => { const t = wrap(x); return Math.abs(t) === Math.PI ? 0 : t / Math.PI; },
      a0half: () => 0,
      coef: (k) => ({ a: 0, b: (2 * (k % 2 ? 1 : -1)) / (Math.PI * k) }),
      formula: "b_k = 2(−1)^(k+1)/(πk)",
    },
    triangle: {
      name: "Triangle wave",
      f: (x) => (2 / Math.PI) * Math.asin(Math.sin(x)),
      a0half: () => 0,
      coef: (k) => ({ a: 0, b: k % 2 ? ((((k - 1) / 2) % 2 ? -1 : 1) * 8) / (Math.PI * Math.PI * k * k) : 0 }),
      formula: "b_k = 8(−1)^((k−1)/2)/(π²k²) for odd k",
    },
    pulse: {
      name: "Pulse train",
      // 1 for |x| < πd, 0 otherwise (duty cycle d)
      f: (x, d) => { const t = Math.abs(wrap(x)), w = Math.PI * d; return t < w ? 1 : t === w ? 0.5 : 0; },
      a0half: (d) => d,
      coef: (k, d) => ({ a: (2 * Math.sin(Math.PI * k * d)) / (Math.PI * k), b: 0 }),
      formula: "a_k = 2 sin(πkd)/(πk), mean d",
    },
  };
  /** First K harmonics: [{k, a, b, amp, phase}], amp·sin(kx + phase) = a cos kx + b sin kx. */
  function fourierCoefficients(wave, K, d) {
    const W = WAVES[wave];
    return Array.from({ length: K }, (_, i) => {
      const k = i + 1, { a, b } = W.coef(k, d);
      return { k, a, b, amp: Math.hypot(a, b), phase: Math.atan2(a, b) };
    });
  }
  /** Partial sum S_K(x) = a0/2 + Σ_{k ≤ K} (a_k cos kx + b_k sin kx). */
  function partialSum(coeffs, x, a0half) {
    let s = a0half || 0;
    for (const c of coeffs) s += c.a * Math.cos(c.k * x) + c.b * Math.sin(c.k * x);
    return s;
  }
  /** Numerical Fourier coefficients of samples f(x_j), x_j = −π + 2πj/M (rectangle rule, exact for trigonometric polynomials). */
  function dftCoefficients(samples, K) {
    const M = samples.length;
    const xs = samples.map((_, j) => -Math.PI + (2 * Math.PI * j) / M);
    const a0half = samples.reduce((s, v) => s + v, 0) / M;
    const out = [];
    for (let k = 1; k <= K; k++) {
      let a = 0, b = 0;
      samples.forEach((v, j) => { a += v * Math.cos(k * xs[j]); b += v * Math.sin(k * xs[j]); });
      a *= 2 / M; b *= 2 / M;
      out.push({ k, a, b, amp: Math.hypot(a, b), phase: Math.atan2(a, b) });
    }
    return { a0half, coeffs: out };
  }
  /** Fraction of the signal's power captured by the first K harmonics (Parseval), estimated on a fine grid. */
  function powerFraction(wave, K, d) {
    const W = WAVES[wave], M = 4096;
    let total = 0;
    for (let j = 0; j < M; j++) { const v = W.f(-Math.PI + ((j + 0.5) * 2 * Math.PI) / M, d); total += v * v; }
    total /= M; // mean of f²
    const m = W.a0half(d);
    let partial = m * m;
    fourierCoefficients(wave, K, d).forEach((c) => (partial += (c.a * c.a + c.b * c.b) / 2));
    return partial / total;
  }
  /** Maximum of the partial sum over one period (used to show the Gibbs overshoot). */
  function maxOfPartialSum(wave, K, d, samples) {
    const co = fourierCoefficients(wave, K, d), a0 = WAVES[wave].a0half(d), n = samples || 4000;
    let m = -Infinity;
    for (let j = 0; j < n; j++) m = Math.max(m, partialSum(co, -Math.PI + (2 * Math.PI * j) / n, a0));
    return m;
  }

  return {
    rng, randInt, diceSumDistribution, rollDice, binomialPMF, longestRun, expectedLongestRunApprox,
    montyGame, montyTheory, WAVES, wrap, fourierCoefficients, partialSum, dftCoefficients, powerFraction, maxOfPartialSum,
  };
});
