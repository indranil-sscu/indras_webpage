/*
 * dsfd-core.js — enumerates double-sided Feynman diagrams (Liouville pathways)
 * for the n-th order response of a multilevel system, lays them out, and
 * renders them as SVG or TikZ. No DOM access; unit-tested in Node
 * (tests/dsfd-core.test.js).
 *
 * Conventions (S. Mukamel, "Principles of Nonlinear Optical Spectroscopy";
 * P. Hamm & M. Zanni, "Concepts and Methods of 2D Infrared Spectroscopy";
 * A. Tokmakoff, lecture notes on nonlinear spectroscopy):
 *  - Time runs upward; the left rail is the ket, the right rail the bra.
 *  - The system starts in the ground-state population |g⟩⟨g|.
 *  - Rotating-wave approximation: each field interaction excites or
 *    de-excites the ket or the bra by one manifold.
 *  - Arrows pointing toward a rail are absorptions (the ket or bra is
 *    excited); arrows pointing away are emissions (de-excitation).
 *  - Wave-vector contribution:
 *        ket excitation  → +k     ket de-excitation → −k
 *        bra excitation  → −k     bra de-excitation → +k
 *    so k_sig = Σ ±k_i with the sign read off each arrow.
 *  - The final (signal) interaction acts on the ket, de-excites it onto the
 *    bra state and leaves a population.
 *  - Each diagram carries the sign (−1)^n, n = number of bra-side
 *    interactions (the signal emission is not counted). The overall minus
 *    sign of excited-state absorption follows from this rule because ESA
 *    pathways have an odd number of bra interactions.
 *  - During interval t_i the coherence |a⟩⟨b| evolves as exp(−iω_ab t_i − Γ_ab t_i).
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.DSFD = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  const SUB = "₀₁₂₃₄₅₆₇₈₉";

  /**
   * Build a level scheme.
   * @param {string[]} singles  names of the singly excited states, e.g. ["e"] or ["a","b"]
   * @param {string[]} doubles  names of doubly excited states, e.g. [] or ["f"]
   * @param {object}  [opts]   {twoQuantum: true} also allows g ↔ doubles transitions, for χ⁽²⁾
   *                           processes through an electronic state (SFG/DFG with g, v, n)
   */
  function levelScheme(singles, doubles, opts) {
    const manifold = { g: 0 };
    singles.forEach((s) => (manifold[s] = 1));
    (doubles || []).forEach((s) => (manifold[s] = 2));
    const states = Object.keys(manifold);
    const two = !!(opts && opts.twoQuantum);
    const allowed = (a, b) => { const dm = manifold[b] - manifold[a]; return dm === 1 || (two && dm === 2 && manifold[a] === 0); };
    // dipole-allowed transitions: between adjacent manifolds (and g ↔ doubles if twoQuantum)
    const up = (s) => states.filter((t) => allowed(s, t));
    const down = (s) => states.filter((t) => allowed(t, s));
    return { manifold, states, up, down, twoQuantum: two };
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
        // signal emission: a dipole-allowed de-excitation of the ket onto the bra state, leaving a population
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
      signs: signs.slice(),
      rho, // rho[i] is the state during time interval t_i (i = 1..N)
      sign: nBra % 2 === 0 ? +1 : -1,
      nBra,
    };
    d.type = classify(d, scheme, signs);
    d.rephasing = isRephasing(d, scheme);
    d.label = d.rho.slice(1).map(([a, b]) => `|${a}⟩⟨${b}|`).join(" → ");
    return d;
  }

  /** GSB / SE / ESA / DQC classification (third order); SFG / DFG at second order. */
  function classify(d, scheme, signs) {
    const N = signs.length, m = scheme.manifold;
    if (N === 1) return "Linear";
    if (N === 2) return signs[0] + signs[1] === 2 ? "SFG" : "DFG";
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

  /**
   * Oscillating phase factor during each interval, written with positive
   * transition frequencies: |a⟩⟨b| with a above b → e^{−iω_ab t},
   * a below b → e^{+iω_ba t}, a population → no oscillation.
   * Returns [{t, sign: −1|+1|0, upper, lower}] and a plain-text form.
   */
  function phaseFactors(d, scheme, intervalNames) {
    const m = scheme.manifold;
    const names = intervalNames || d.rho.slice(1).map((_, i) => "t" + SUB[i + 1]);
    const parts = d.rho.slice(1).map(([a, b], i) => {
      if (a === b) return { t: names[i], sign: 0, upper: a, lower: b };
      // within one manifold (e.g. |a⟩⟨b| of a dimer) keep the ket-first order
      if (m[a] === m[b]) return { t: names[i], sign: -1, upper: a, lower: b };
      return m[a] > m[b] ? { t: names[i], sign: -1, upper: a, lower: b } : { t: names[i], sign: +1, upper: b, lower: a };
    });
    const text = parts.filter((p) => p.sign !== 0).map((p) => `${p.sign < 0 ? "−" : "+"}iω_${p.upper}${p.lower}·${p.t}`).join(" ");
    return { parts, text: text ? `exp(${text})` : "1" };
  }

  /** Phase-matching direction as a string, e.g. "−k₁ + k₂ + k₃". Optional field labels replace the subscripts. */
  function direction(signs, fieldLabels) {
    const lab = (i) => (fieldLabels && fieldLabels[i] ? fieldLabels[i] : SUB[i + 1]);
    const terms = signs.map((s, i) => ({ s, l: lab(i) }));
    if (fieldLabels) {
      // combine repeated labels, e.g. +k_pu − k_pu + k_pr → k_pr
      const net = {};
      terms.forEach((t) => (net[t.l] = (net[t.l] || 0) + t.s));
      const kept = Object.entries(net).filter(([, v]) => v !== 0);
      if (!kept.length) return "0";
      return kept.map(([l, v], i) => (i === 0 ? (v < 0 ? "−" : "") : v < 0 ? " − " : " + ") + (Math.abs(v) > 1 ? Math.abs(v) : "") + "k_" + l).join("");
    }
    return terms.map((t, i) => (i === 0 ? (t.s < 0 ? "−" : "") : t.s < 0 ? " − " : " + ") + "k" + t.l).join("");
  }

  const PRESETS = {
    rephasing: { name: "2D rephasing (photon echo)", signs: [-1, +1, +1], scheme: { singles: "e", doubles: "f" } },
    nonrephasing: { name: "2D non-rephasing", signs: [+1, -1, +1], scheme: { singles: "e", doubles: "f" } },
    pumpprobe: {
      name: "Transient absorption (pump–probe)", signSets: [[-1, +1, +1], [+1, -1, +1]],
      fieldLabels: ["pu", "pu", "pr"], intervals: ["t₁ = 0", "T", "t"], scheme: { singles: "e", doubles: "f" },
    },
    dqc: { name: "Double-quantum coherence", signs: [+1, +1, -1], scheme: { singles: "e", doubles: "f" } },
    sfg: { name: "Sum-frequency generation (IR + vis)", signs: [+1, +1], fieldLabels: ["IR", "vis"], scheme: { singles: "v", doubles: "n", twoQuantum: true } },
    dfg: { name: "Difference-frequency generation", signs: [+1, -1], fieldLabels: ["1", "2"], scheme: { singles: "v", doubles: "n", twoQuantum: true } },
    linear: { name: "Linear absorption", signs: [+1], scheme: { singles: "e", doubles: "" } },
  };

  /* ------------------------------------------------------------------ *
   * Layout: geometry shared by the SVG and TikZ renderers.
   * Coordinates are in px with y pointing down (SVG convention).
   * ------------------------------------------------------------------ */
  const SLOPE = Math.tan((15 * Math.PI) / 180); // interaction arrows rise at 15°

  /**
   * @param {object} d diagram from enumerate()
   * @param {object} o {intervals: [names for t₁…t_n], fieldLabels: [...], timeAxis: bool (default off), header: string}
   */
  function layout(d, o) {
    o = o || {};
    const n = d.interactions.length;
    const withAxis = o.timeAxis === true;
    const gap = 50, axisX = 16, kx = withAxis ? 140 : 96, bx = kx + 60, L = 46, dy = L * SLOPE;
    const top0 = o.header ? 34 : 18;
    const yAt = (j) => top0 + gap * (n + 1 - j) + 38; // j = 0 initial … n+1 signal
    const bottom = yAt(0) + 10, W = bx + 100, H = bottom + 10;
    const SUBS = (i) => (o.fieldLabels && o.fieldLabels[i] ? "_{" + o.fieldLabels[i] + "}" : "_" + (i + 1));
    const lay = { W, H, n, rails: [], arrows: [], rho: [], axis: null, brackets: [], header: o.header || "" };
    const railTop = yAt(n + 1) - 36;
    lay.rails.push({ x: kx, y1: bottom, y2: railTop }, { x: bx, y1: bottom, y2: railTop });
    // density-matrix labels: initial below the first interaction, then between interactions
    const mid = (kx + bx) / 2;
    lay.rho.push({ x: mid, y: (yAt(0) + yAt(1)) / 2 + 4, ket: d.rho[0][0], bra: d.rho[0][1] });
    for (let j = 1; j <= n; j++) lay.rho.push({ x: mid, y: (yAt(j) + yAt(j + 1)) / 2 + 4, ket: d.rho[j][0], bra: d.rho[j][1] });
    lay.rho.push({ x: mid, y: (yAt(n + 1) + railTop) / 2 + 2, ket: d.emission.to, bra: d.emission.to });
    // interaction arrows
    d.interactions.concat([d.emission]).forEach((s, idx) => {
      const y = yAt(idx + 1), sig = s.kind === "emit";
      let a;
      if (s.side === "ket") {
        a = s.kind === "up" ? { x1: kx - L, y1: y + dy, x2: kx - 1.5, y2: y } : { x1: kx - 1.5, y1: y, x2: kx - L, y2: y - dy };
        a.lx = Math.min(a.x1, a.x2) - 4; a.anchor = "end";
      } else {
        a = s.kind === "up" ? { x1: bx + L, y1: y + dy, x2: bx + 1.5, y2: y } : { x1: bx + 1.5, y1: y, x2: bx + L, y2: y - dy };
        a.lx = Math.max(a.x1, a.x2) + 4; a.anchor = "start";
      }
      a.ly = (a.y1 + a.y2) / 2 + 4;
      a.signal = sig;
      a.inward = s.kind === "up";
      a.side = s.side;
      a.text = sig ? "k_{sig}" : (s.sign > 0 ? "+" : "−") + "k" + SUBS(s.k - 1);
      lay.arrows.push(a);
    });
    // time axis with bracketed intervals
    if (withAxis) {
      lay.axis = { x: axisX, y1: bottom, y2: railTop - 4 };
      // interval t_j runs from interaction j to interaction j+1 (where the arrows touch the rails)
      const uni = (t) => String(t).replace(/_(\d)/g, (_, d) => SUB[+d]);
      const names = (o.intervals || Array.from({ length: n }, (_, i) => "t_" + (i + 1))).map(uni);
      lay.ticks = Array.from({ length: n + 1 }, (_, j) => yAt(j + 1));
      for (let j = 1; j <= n; j++) lay.brackets.push({ x: axisX + 8, y1: yAt(j) - 1.5, y2: yAt(j + 1) + 1.5, ym: (yAt(j) + yAt(j + 1)) / 2, text: names[j - 1] });
    }
    return lay;
  }

  /* ---------------- text helpers ---------------- */
  // "k_{sig}", "k_1", "t_2", "t₁ = 0" → SVG tspans / LaTeX
  function svgText(s) {
    const esc = (t) => String(t).replace(/&/g, "&amp;").replace(/</g, "&lt;");
    return esc(s).replace(/_\{([^}]*)\}|_(\w)/g, (_, a, b) => `<tspan baseline-shift="sub" font-size="75%">${a || b}</tspan>`);
  }
  function texText(s) {
    const map = { "₀": "_0", "₁": "_1", "₂": "_2", "₃": "_3", "₄": "_4", "₅": "_5", "−": "-", "τ": "\\tau ", "ω": "\\omega " };
    let t = String(s).replace(/[₀-₅−τω]/g, (c) => map[c]);
    t = t.replace(/k_\{sig\}/g, "\\mathbf{k}_{\\mathrm{s}}").replace(/k_\{([^}]*)\}/g, "\\mathbf{k}_{\\mathrm{$1}}").replace(/k_(\w)/g, "\\mathbf{k}_{$1}");
    return t;
  }
  const texState = (s) => String(s).replace(/[^A-Za-z0-9']/g, "");

  /* ---------------- palettes ---------------- */
  const PATHWAY_COLORS = { GSB: "#2B6CB0", SE: "#2C7A7B", ESA: "#C53030", DQC: "#6B46C1", SFG: "#2B6CB0", DFG: "#2C7A7B", Linear: "#4A5568", "—": "#4A5568" };
  /** palette: {ink, bg, muted, accent(type)} — "print" is pure black on white. */
  function palette(kind, dark) {
    if (kind === "mono") return { ink: "#000000", bg: "#ffffff", muted: "#000000", accent: () => "#000000", mono: true };
    if (dark) return { ink: "#E6E8EB", bg: "#161A20", muted: "#8A919C", accent: (t) => ({ GSB: "#63B3ED", SE: "#4FD1C5", ESA: "#FC8181", DQC: "#B794F4" }[t] || "#CBD5E0") };
    return { ink: "#111111", bg: "#ffffff", muted: "#666666", accent: (t) => PATHWAY_COLORS[t] || "#4A5568" };
  }

  /** Points of a sinusoidal (wavy) arrow shaft from (x1,y1) to (x2,y2), leaving `tail` px straight for the head. */
  function wavyPath(x1, y1, x2, y2, amp, wl, tail) {
    const len = Math.hypot(x2 - x1, y2 - y1), ux = (x2 - x1) / len, uy = (y2 - y1) / len, nx = -uy, ny = ux;
    const L = len - tail, pts = [];
    for (let s = 0; s <= L; s += 1) { const w = amp * Math.sin((2 * Math.PI * s) / wl); pts.push([x1 + ux * s + nx * w, y1 + uy * s + ny * w]); }
    pts.push([x2, y2]);
    return "M" + pts.map((p) => p[0].toFixed(2) + "," + p[1].toFixed(2)).join(" L");
  }

  /**
   * SVG string for one diagram (inline presentation attributes only, so it
   * survives export to Illustrator / Inkscape / LaTeX).
   * opts: {pal, wavy, type, sign}
   */
  function toSVG(lay, opts) {
    const P = opts.pal, accent = P.accent(opts.type), font = "Times New Roman, STIX Two Text, serif";
    const out = [];
    const mk = (id, col) => `<marker id="${id}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse" markerUnits="userSpaceOnUse"><path d="M0,1 L10,5 L0,9 z" fill="${col}"/></marker>`;
    const uid = "m" + Math.random().toString(36).slice(2, 8);
    out.push(`<defs>${mk(uid + "i", P.ink)}${mk(uid + "s", accent)}</defs>`);
    out.push(`<rect width="${lay.W}" height="${lay.H}" fill="${P.bg}"/>`);
    if (lay.header) out.push(`<text x="10" y="18" font-family="${font}" font-size="13" fill="${accent}" font-weight="bold">${svgText(lay.header)}</text>`);
    lay.rails.forEach((r) => out.push(`<line x1="${r.x}" y1="${r.y1}" x2="${r.x}" y2="${r.y2}" stroke="${P.ink}" stroke-width="2" stroke-linecap="square"/>`));
    lay.rho.forEach((r) => out.push(`<text x="${r.x}" y="${r.y}" text-anchor="middle" font-family="${font}" font-size="13" fill="${P.ink}">|${svgText(r.ket)}⟩⟨${svgText(r.bra)}|</text>`));
    lay.arrows.forEach((a) => {
      const col = a.signal ? accent : P.ink, marker = `url(#${uid}${a.signal ? "s" : "i"})`;
      if (a.signal) {
        out.push(`<line x1="${a.x1}" y1="${a.y1}" x2="${a.x2}" y2="${a.y2}" stroke="${col}" stroke-width="1.6" stroke-dasharray="5,3" marker-end="${marker}"/>`);
      } else if (opts.wavy) {
        out.push(`<path d="${wavyPath(a.x1, a.y1, a.x2, a.y2, 2.6, 8, 7)}" fill="none" stroke="${col}" stroke-width="1.4" marker-end="${marker}"/>`);
      } else {
        out.push(`<line x1="${a.x1}" y1="${a.y1}" x2="${a.x2}" y2="${a.y2}" stroke="${col}" stroke-width="1.5" marker-end="${marker}"/>`);
      }
      out.push(`<text x="${a.lx}" y="${a.ly}" font-family="${font}" font-size="13" font-style="italic" fill="${col}" text-anchor="${a.anchor}">${svgText(a.text)}</text>`);
    });
    if (lay.axis) {
      const A = lay.axis;
      out.push(`<line x1="${A.x}" y1="${A.y1}" x2="${A.x}" y2="${A.y2}" stroke="${P.muted}" stroke-width="1" marker-end="url(#${uid}i)"/>`);
      out.push(`<text x="${A.x}" y="${A.y2 - 6}" text-anchor="middle" font-family="${font}" font-size="13" font-style="italic" fill="${P.muted}">t</text>`);
      (lay.ticks || []).forEach((y) => out.push(`<line x1="${A.x - 3}" y1="${y}" x2="${A.x + 3}" y2="${y}" stroke="${P.muted}" stroke-width="1"/>`));
      lay.brackets.forEach((b) => {
        out.push(`<path d="M${A.x + 3},${b.y1} H${b.x} V${b.y2} H${A.x + 3}" fill="none" stroke="${P.muted}" stroke-width="1"/>`);
        out.push(`<text x="${b.x + 5}" y="${b.ym}" dominant-baseline="central" font-family="${font}" font-size="12" font-style="italic" fill="${P.muted}">${svgText(b.text)}</text>`);
      });
    }
    return `<svg xmlns="http://www.w3.org/2000/svg" width="${lay.W}" height="${lay.H}" viewBox="0 0 ${lay.W} ${lay.H}">${out.join("")}</svg>`;
  }

  /**
   * Standalone LaTeX/TikZ document. items: [{lay, type, sign}] → one tikzpicture
   * (one page of the standalone PDF) per diagram.
   */
  function toTikZ(items, opts) {
    opts = opts || {};
    const mono = !!opts.mono, wavy = !!opts.wavy;
    const used = [...new Set(items.map((it) => it.type))];
    const cname = (t) => "path" + String(t).replace(/[^A-Za-z]/g, "");
    const f = (v) => (+v).toFixed(1);
    const lines = [
      "% Double-sided Feynman diagrams — generated by https://indranil-sscu.github.io/indras_webpage/tools/dsfd/",
      "% Compile with pdflatex; each diagram is one page of the PDF.",
      "\\documentclass[tikz,border=4pt]{standalone}",
      "\\usetikzlibrary{arrows.meta, decorations.pathmorphing}",
      ...used.map((t) => `\\definecolor{${cname(t)}}{HTML}{${(mono ? "000000" : (PATHWAY_COLORS[t] || "#4A5568").slice(1))}}`),
      "\\tikzset{",
      "  rail/.style={line width=1.5pt},",
      "  field/.style={-{Stealth[length=5pt,width=4pt]}, line width=0.9pt" + (wavy ? ", decorate, decoration={snake, amplitude=1.6pt, segment length=5pt, post length=5pt}" : "") + "},",
      "  signal/.style={-{Stealth[length=5pt,width=4pt]}, line width=0.9pt, dashed},",
      "  axis/.style={-{Stealth[length=4pt]}, line width=0.5pt, gray},",
      "}",
      "\\begin{document}",
    ];
    items.forEach((it) => {
      const L = it.lay, col = cname(it.type);
      lines.push("\\begin{tikzpicture}[x=0.75pt, y=-0.75pt, every node/.style={font=\\small}]");
      if (L.header) lines.push(`  \\node[anchor=west, text=${col}, font=\\small\\bfseries] at (10,14) {${L.header.replace(/[#%&_]/g, (c) => "\\" + c).replace(/−/g, "$-$").replace(/[^\x00-\x7F]/g, "")}};`);
      L.rails.forEach((r) => lines.push(`  \\draw[rail] (${f(r.x)},${f(r.y1)}) -- (${f(r.x)},${f(r.y2)});`));
      L.rho.forEach((r) => lines.push(`  \\node at (${f(r.x)},${f(r.y - 4)}) {$|${texState(r.ket)}\\rangle\\langle ${texState(r.bra)}|$};`));
      L.arrows.forEach((a) => {
        lines.push(`  \\draw[${a.signal ? "signal, " + col : "field"}] (${f(a.x1)},${f(a.y1)}) -- (${f(a.x2)},${f(a.y2)});`);
        lines.push(`  \\node[anchor=${a.anchor === "end" ? "east" : "west"}${a.signal ? ", text=" + col : ""}] at (${f(a.lx)},${f(a.ly - 4)}) {$${texText(a.text)}$};`);
      });
      if (L.axis) {
        lines.push(`  \\draw[axis] (${f(L.axis.x)},${f(L.axis.y1)}) -- (${f(L.axis.x)},${f(L.axis.y2)}) node[above, text=gray] {$t$};`);
        (L.ticks || []).forEach((y) => lines.push(`  \\draw[gray, line width=0.4pt] (${f(L.axis.x - 3)},${f(y)}) -- (${f(L.axis.x + 3)},${f(y)});`));
        L.brackets.forEach((b) => {
          lines.push(`  \\draw[gray, line width=0.4pt] (${f(L.axis.x + 3)},${f(b.y1)}) -- (${f(b.x)},${f(b.y1)}) -- (${f(b.x)},${f(b.y2)}) -- (${f(L.axis.x + 3)},${f(b.y2)});`);
          lines.push(`  \\node[anchor=west, text=gray, font=\\footnotesize] at (${f(b.x + 2)},${f(b.ym)}) {$${texText(b.text)}$};`);
        });
      }
      lines.push("\\end{tikzpicture}");
    });
    lines.push("\\end{document}");
    return lines.join("\n") + "\n";
  }

  /** Insert a pHYs chunk (pixels per metre) into a PNG byte array, e.g. 300 dpi → 11811 px/m. */
  function pngWithDpi(bytes, dpi) {
    const ppm = Math.round(dpi / 0.0254);
    const crcTable = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
    const crc = (buf) => { let c = 0xffffffff; for (let i = 0; i < buf.length; i++) c = crcTable[(c ^ buf[i]) & 0xff] ^ (c >>> 8); return (c ^ 0xffffffff) >>> 0; };
    const chunk = new Uint8Array(4 + 4 + 9 + 4);
    const dv = new DataView(chunk.buffer);
    dv.setUint32(0, 9); chunk.set([0x70, 0x48, 0x59, 0x73], 4); // "pHYs"
    dv.setUint32(8, ppm); dv.setUint32(12, ppm); chunk[16] = 1; // unit: metre
    dv.setUint32(17, crc(chunk.subarray(4, 17)));
    const ihdrEnd = 8 + 25; // signature + IHDR chunk (13 data bytes + 12)
    const out = new Uint8Array(bytes.length + chunk.length);
    out.set(bytes.subarray(0, ihdrEnd), 0); out.set(chunk, ihdrEnd); out.set(bytes.subarray(ihdrEnd), ihdrEnd + chunk.length);
    return out;
  }

  return { levelScheme, enumerate, dipoles, propagators, phaseFactors, direction, layout, toSVG, toTikZ, palette, pngWithDpi, PATHWAY_COLORS, PRESETS };
});
