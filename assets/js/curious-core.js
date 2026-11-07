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

  /* ---------------- the birthday paradox ---------------- */
  /** Exact probability that at least two of n people share a birthday, with d equally likely days:
   *  1 − ∏_{k=0}^{n−1} (d − k)/d. */
  function birthdayExact(n, d) {
    d = d || 365;
    if (n > d) return 1;
    let q = 1;
    for (let k = 0; k < n; k++) q *= (d - k) / d;
    return 1 - q;
  }
  /** Smallest group with at least probability p of a shared birthday (23 for p = ½, d = 365). */
  function birthdayThreshold(p, d) {
    d = d || 365;
    for (let n = 1; n <= d + 1; n++) if (birthdayExact(n, d) >= p) return n;
    return d + 1;
  }
  /** One room of n people with uniformly random birthdays 0 … d−1.
   *  Returns {days, shared, pairs (number of matching pairs), sharedDays: Set of days held by ≥ 2 people}. */
  function birthdayRoom(n, d, r) {
    d = d || 365;
    const days = Array.from({ length: n }, () => randInt(r, d)), count = new Map();
    days.forEach((x) => count.set(x, (count.get(x) || 0) + 1));
    let pairs = 0; const sharedDays = new Set();
    count.forEach((c, x) => { if (c > 1) { pairs += (c * (c - 1)) / 2; sharedDays.add(x); } });
    return { days, shared: pairs > 0, pairs, sharedDays };
  }
  /** Fast check used for large batches: does a room of n share a birthday? (stops at the first match) */
  function birthdayHasMatch(n, d, r, seen) {
    d = d || 365;
    if (n > d) return true;
    seen = seen && seen.length >= d ? seen.fill(0) : new Uint8Array(d);
    for (let i = 0; i < n; i++) { const x = randInt(r, d); if (seen[x]) return true; seen[x] = 1; }
    return false;
  }
  /** Expected number of matching pairs, C(n, 2)/d (exact, by linearity of expectation). */
  const birthdayExpectedPairs = (n, d) => (n * (n - 1)) / 2 / (d || 365);

  /* ---------------- Simpson's paradox ---------------- */
  /**
   * Two treatments (A, B) in two groups (1, 2). t = {A: [[s1, n1], [s2, n2]], B: [[s1, n1], [s2, n2]]},
   * s = successes, n = cases. Returns the success rates in each group and overall, which treatment wins
   * each comparison, and whether the overall comparison reverses both within-group comparisons.
   */
  function simpson(t) {
    const rate = (s, n) => (n > 0 ? s / n : NaN);
    const out = {};
    ["A", "B"].forEach((k) => {
      const [[s1, n1], [s2, n2]] = t[k];
      out[k] = { g1: rate(s1, n1), g2: rate(s2, n2), all: rate(s1 + s2, n1 + n2), share2: rate(n2, n1 + n2), n: n1 + n2, s: s1 + s2 };
    });
    const cmp = (a, b) => (!isFinite(a) || !isFinite(b) ? 0 : a > b ? 1 : a < b ? -1 : 0);
    const w1 = cmp(out.A.g1, out.B.g1), w2 = cmp(out.A.g2, out.B.g2), wAll = cmp(out.A.all, out.B.all);
    out.winner = { g1: w1, g2: w2, all: wAll }; // +1: A higher, −1: B higher, 0: tie
    out.reversal = w1 !== 0 && w1 === w2 && wAll === -w1;
    return out;
  }
  /**
   * Keep a treatment's group success rates and its total number of cases, and move cases between the
   * groups so that a fraction `share2` of them is in group 2. Successes are rounded to whole people.
   */
  function reallocate(row, share2) {
    const [[s1, n1], [s2, n2]] = row, N = n1 + n2;
    const p1 = n1 ? s1 / n1 : 0, p2 = n2 ? s2 / n2 : 0;
    const m2 = Math.round(N * Math.min(1, Math.max(0, share2))), m1 = N - m2;
    return [[Math.round(p1 * m1), m1], [Math.round(p2 * m2), m2]];
  }
  /** Classic data sets. Kidney stones: Charig et al., BMJ 292, 879 (1986); A = open surgery, B = percutaneous nephrolithotomy. */
  const SIMPSON_PRESETS = {
    kidney: { name: "Kidney-stone treatments (1986)", groups: ["Small stones", "Large stones"], labels: ["Open surgery", "Keyhole surgery"], outcome: "treatment success",
      A: [[81, 87], [192, 263]], B: [[234, 270], [55, 80]] },
    batting: { name: "Batting averages, 1995–96", groups: ["1995", "1996"], labels: ["Derek Jeter", "David Justice"], outcome: "hits per at-bat",
      A: [[12, 48], [183, 582]], B: [[104, 411], [45, 140]] },
    none: { name: "No paradox (balanced groups)", groups: ["Group 1", "Group 2"], labels: ["Treatment A", "Treatment B"], outcome: "success",
      A: [[90, 100], [70, 100]], B: [[85, 100], [65, 100]] },
  };

  return {
    rng, randInt, diceSumDistribution, rollDice, binomialPMF, longestRun, expectedLongestRunApprox,
    montyGame, montyTheory, WAVES, wrap, fourierCoefficients, partialSum, dftCoefficients, powerFraction, maxOfPartialSum,
    birthdayExact, birthdayThreshold, birthdayRoom, birthdayHasMatch, birthdayExpectedPairs, simpson, reallocate, SIMPSON_PRESETS,
  };
});
