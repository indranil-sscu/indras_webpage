/*
 * dsfd-core.js — enumerates double-sided Feynman diagrams (Liouville pathways)
 * for the n-th order response of a multilevel system.
 * No DOM access; unit-tested in Node (tests/dsfd-core.test.js).
 *
 * Conventions (Mukamel, "Principles of Nonlinear Optical Spectroscopy";
 * Hamm & Zanni, "Concepts and Methods of 2D Infrared Spectroscopy"):
 *  - Time runs upward; the left line is the ket, the right line the bra.
 *  - The system starts in the ground-state population |g⟩⟨g|.
 *  - Rotating-wave approximation: each interaction either excites or
 *    de-excites the ket or the bra by one manifold.
 *  - Wave-vector contribution:
 *        ket excitation  → +k     ket de-excitation → −k
 *        bra excitation  → −k     bra de-excitation → +k
 *  - The final (signal-emitting) interaction acts on the ket, de-excites it
 *    and leaves a population (diagonal) state.
 *  - Each diagram carries the sign (−1)^n, n = number of bra interactions
 *    (the emission is not counted).
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.DSFD = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /**
   * Build a level scheme.
   * @param {string[]} singles  names of the singly excited states, e.g. ["e"] or ["a","b"]
   * @param {string[]} doubles  names of doubly excited states, e.g. [] or ["f"]
   */
  function levelScheme(singles, doubles) {
    const manifold = { g: 0 };
    singles.forEach((s) => (manifold[s] = 1));
    (doubles || []).forEach((s) => (manifold[s] = 2));
    const states = Object.keys(manifold);
    // dipole-allowed transitions: between adjacent manifolds only
    const up = (s) => states.filter((t) => manifold[t] === manifold[s] + 1);
    const down = (s) => states.filter((t) => manifold[t] === manifold[s] - 1);
    return { manifold, states, up, down };
  }

  /**
   * Enumerate all pathways of order N = signs.length.
   * @param {number[]} signs  wave-vector sign (+1 / −1) for k_1 … k_N, in time order
   * @param {object}   scheme from levelScheme()
   * @returns {object[]} diagrams
   */
  function enumerate(signs, scheme) {
    const N = signs.length;
    const out = [];
    const step = (j, ket, bra, steps) => {
      if (j === N) {
        // emission: ket de-excites onto the bra state, giving a population
        if (scheme.down(ket).includes(bra)) {
          out.push(makeDiagram(signs, steps, { side: "ket", from: ket, to: bra, kind: "emit" }, scheme));
        }
        return;
      }
      const s = signs[j];
      const opts = [];
      if (s > 0) {
        scheme.up(ket).forEach((t) => opts.push({ side: "ket", from: ket, to: t, kind: "up" }));
        scheme.down(bra).forEach((t) => opts.push({ side: "bra", from: bra, to: t, kind: "down" }));
      } else {
        scheme.down(ket).forEach((t) => opts.push({ side: "ket", from: ket, to: t, kind: "down" }));
        scheme.up(bra).forEach((t) => opts.push({ side: "bra", from: bra, to: t, kind: "up" }));
      }
      for (const o of opts) {
        const nk = o.side === "ket" ? o.to : ket;
        const nb = o.side === "bra" ? o.to : bra;
        step(j + 1, nk, nb, steps.concat([{ ...o, k: j + 1, sign: s }]));
      }
    };
    step(0, "g", "g", []);
    return out;
  }

  function makeDiagram(signs, steps, emission, scheme) {
    // density-matrix element after each interaction (index 0 = initial)
    const rho = [["g", "g"]];
    for (const s of steps) {
      const [k, b] = rho[rho.length - 1];
      rho.push(s.side === "ket" ? [s.to, b] : [k, s.to]);
    }
    const nBra = steps.filter((s) => s.side === "bra").length;
    const d = {
      interactions: steps,
      emission,
      rho, // rho[i] is the state during time interval t_i (i = 1..N)
      sign: nBra % 2 === 0 ? +1 : -1,
      nBra,
    };
    d.type = classify(d, scheme, signs);
    d.rephasing = isRephasing(d, scheme);
    d.label = d.rho.slice(1).map(([a, b]) => `|${a}⟩⟨${b}|`).join(" → ");
    return d;
  }

  /** GSB / SE / ESA / DQC classification (third order); generic otherwise. */
  function classify(d, scheme, signs) {
    const N = signs.length, m = scheme.manifold;
    if (N === 1) return "Linear";
    if (N !== 3) return "—";
    const [k2, b2] = d.rho[2];
    if (Math.abs(m[k2] - m[b2]) === 2) return "DQC"; // |f⟩⟨g| or |g⟩⟨f| during t2
    if (k2 === "g" && b2 === "g") return "GSB";
    const third = d.interactions[2];
    if (third.kind === "up") return "ESA";
    return "SE";
  }

  /** For third order: rephasing if the t1 and t3 coherences oscillate with opposite sign. */
  function isRephasing(d, scheme) {
    if (d.interactions.length !== 3) return null;
    const m = scheme.manifold;
    const sgn = ([a, b]) => Math.sign(m[a] - m[b]);
    return sgn(d.rho[1]) * sgn(d.rho[3]) < 0;
  }

  /** Transition dipoles in time order, lower state first, e.g. ["ge","ge","ge","ge"] or ["ge","ef",…]. */
  function dipoles(d, scheme) {
    const m = scheme.manifold;
    return d.interactions.concat([d.emission]).map((s) => (m[s.from] < m[s.to] ? s.from + s.to : s.to + s.from));
  }

  /** Propagator factors G_{ab}(t_i) = exp(−i ω_ab t_i − Γ_ab t_i) for each interval. */
  function propagators(d) {
    return d.rho.slice(1).map(([a, b], i) => ({ t: i + 1, a, b, population: a === b }));
  }

  /** Phase-matching direction as a string, e.g. "−k₁ + k₂ + k₃". */
  function direction(signs) {
    const sub = "₀₁₂₃₄₅₆₇₈₉";
    return signs.map((s, i) => (i === 0 ? (s < 0 ? "−" : "") : s < 0 ? " − " : " + ") + "k" + sub[i + 1]).join("");
  }

  const PRESETS = {
    linear: { name: "Linear absorption", signs: [+1] },
    rephasing: { name: "Rephasing (photon echo)", signs: [-1, +1, +1] },
    nonrephasing: { name: "Non-rephasing", signs: [+1, -1, +1] },
    dqc: { name: "Double-quantum coherence", signs: [+1, +1, -1] },
  };

  return { levelScheme, enumerate, dipoles, propagators, direction, PRESETS };
});
