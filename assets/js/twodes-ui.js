/*
 * twodes-ui.js: the interface model of the 2DES simulator (v2). Pure data and functions, no DOM:
 * Explore-mode presets and their simplified controls, the parameter categories of Advanced mode,
 * the capability matrix, and the job builder that turns interface state into an engine job.
 * Tested in Node (tests/twodes-ui.test.js) so that the interface provably feeds the engine the same
 * inputs as before the redesign.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory(require("./twodes-core.js"));
  else root.TwoDESUI = factory(root.TwoDES);
})(typeof self !== "undefined" ? self : this, function (D) {
  "use strict";
  const cm = D.cm;

  /** Numerics used by the Explore presets (hidden in Explore mode, editable in the Research Workbench). */
  const NUM_OPTICAL = { dt: 4, n1: 128, n3: 128, pad: 4, window: "none", broaden: 0, halfFirst: true, convention: "standard", carrier: cm(15000), autoCarrier: true, nonuniform: "direct" };
  const NUM_PAPER = Object.assign({}, NUM_OPTICAL, { dt: 10, carrier: cm(12870) });

  const rate = { // dephasing rate Γ = 1/T₂ (rad/fs) ↔ T₂ (fs)
    get: (p) => 1 / p.tau2,
    set: (p, g) => { p.tau2 = 1 / g; },
    min: cm(10), max: cm(300), step: cm(1),
  };
  const centre = { // centre of a pair of transitions; the splitting is kept
    get: (p) => (p.e1 + p.e2) / 2,
    set: (p, v) => { const d = p.e2 - p.e1; p.e1 = v - d / 2; p.e2 = v + d / 2; },
  };

  /**
   * Explore-mode systems. Each maps the four beginner controls (transition energy, coupling, dephasing, waiting
   * time) onto real model parameters; a control is null with a reason when it does not apply to that system.
   */
  const SYSTEMS = {
    tls: {
      name: "Single two-level system", model: "tls", params: {}, num: NUM_OPTICAL, T: 0, range: null,
      about: "One electronic transition (a ground and an excited state). The 2D spectrum has a single peak on the diagonal: the transition is excited at ω₁ and detected at the same frequency ω₃. Its width comes from dephasing, the loss of phase memory of the optical coherence.",
      energy: { label: "Transition energy", tip: "Energy gap between the ground and excited state; sets where the peak sits on both axes.", get: (p) => p.w, set: (p, v) => { p.w = v; } },
      coupling: null, couplingWhy: "A single transition has nothing to couple to.",
      dephasing: Object.assign({ label: "Dephasing rate Γ", tip: "How fast the optical coherence loses phase memory (Γ = 1/T₂). The peak's full width at half maximum is 2Γ on both axes." }, rate),
    },
    uncoupled: {
      name: "Two uncoupled transitions", model: "dimer", params: { J: 0, bx: 0 }, num: NUM_OPTICAL, T: 0, range: null,
      about: "Two independent transitions, 600 cm⁻¹ apart. Each gives its own diagonal peak. Cross-peak pathways (ground-state bleach, stimulated emission and excited-state absorption) exist but cancel exactly, because exciting one transition does not change the other.",
      energy: Object.assign({ label: "Transition energy (centre of the pair)", tip: "Mean energy of the two transitions; their 600 cm⁻¹ splitting is kept (change it in Advanced Models)." }, centre),
      coupling: null, couplingWhy: "Uncoupled by definition: choose “Coupled two-level system” to add coupling.",
      dephasing: Object.assign({ label: "Dephasing rate Γ", tip: "Optical dephasing rate of both transitions (Γ = 1/T₂); peak FWHM = 2Γ." }, rate),
    },
    coupled: {
      name: "Coupled two-level system", model: "dimer", params: { J: cm(300) }, num: NUM_OPTICAL, T: 0, range: null,
      about: "Two transitions coupled electronically (an excitonic dimer). Coupling mixes them into two excitons: the peaks move apart, dipole strength is redistributed, and cross peaks appear because the pathways no longer cancel. During T the excitons can beat against each other (electronic coherence).",
      energy: Object.assign({ label: "Transition energy (centre of the pair)", tip: "Mean energy of the two uncoupled (site) transitions; their 600 cm⁻¹ splitting is kept." }, centre),
      coupling: { label: "Coupling strength J", tip: "Electronic (excitonic) coupling between the two transitions. Larger J mixes them more strongly and splits the excitons further apart.", get: (p) => p.J, set: (p, v) => { p.J = v; }, min: 0, max: cm(800), step: cm(5), unit: "energy" },
      dephasing: Object.assign({ label: "Dephasing rate Γ", tip: "Optical dephasing rate of the exciton transitions (Γ = 1/T₂); peak FWHM = 2Γ." }, rate),
    },
    vibronic: {
      name: "Simple vibronic system", model: "vibronic", params: {}, num: NUM_OPTICAL, T: 0, range: null,
      about: "One electronic transition coupled to one molecular vibration. Absorption shows a progression of peaks at ω₀₀ + nω_v, and emission at ω₀₀ − nω_v. Vibrational wavepackets make the peaks oscillate in waiting time with the vibrational period.",
      energy: { label: "0–0 transition energy", tip: "Energy of the lowest vibronic transition (no vibrational quanta in either state).", get: (p) => p.w, set: (p, v) => { p.w = v; } },
      coupling: { label: "Vibronic coupling (Huang–Rhys S)", tip: "How strongly the vibration is displaced on excitation. S sets the relative intensities of the vibronic peaks (0–1 / 0–0 ≈ S).", get: (p) => p.S, set: (p, v) => { p.S = v; }, min: 0, max: 2, step: 0.05, unit: "number" },
      dephasing: Object.assign({ label: "Dephasing rate Γ", tip: "Electronic (optical) dephasing rate (Γ = 1/T₂); each vibronic peak has FWHM ≈ 2Γ." }, rate),
    },
    brownian: {
      name: "Brownian oscillator (BChl a, Sahu & Tiwari 2023)", model: "brownian", params: {}, num: NUM_PAPER, T: 50, range: [11800, 14200],
      about: "Bacteriochlorophyll a at 77 K, as in Sahu & Tiwari, J. Phys. Chem. Lett. 14, 4617 (2023): one electronic transition coupled to five molecular vibrations and a fluctuating environment (Brownian oscillators), with static disorder. Disorder stretches the peak along the diagonal; the environment sets its width across it.",
      energy: { label: "Electronic gap ω_eg", tip: "Vertical (Franck–Condon) energy gap; the absorption maximum sits close to it.", get: (p) => p.w, set: (p, v) => { p.w = v; } },
      coupling: null, couplingWhy: "Vibrational and environmental couplings are set by the oscillator table (Advanced Models).",
      dephasing: null, dephasingWhy: "Dephasing arises from the Brownian-oscillator bath (anti-diagonal width ≈ 65 cm⁻¹).",
    },
  };
  /** Workbench/Advanced examples, kept from v1 so that existing links and expectations still hold. */
  const EXAMPLES = {
    sahu: { label: "Sahu & Tiwari 2023: BChl a monomer (Fig. 1)", system: "brownian" },
    ex1: { label: "1. Isolated two-level system", system: "tls" },
    ex2: { label: "2. Two uncoupled transitions", system: "uncoupled" },
    ex3: { label: "3. Coupled two-state system (dimer)", system: "coupled" },
    ex4: { label: "4. Vibronic system", system: "vibronic" },
    ex5: { label: "GSB, SE and ESA (three-level ladder)", model: "ladder", params: {}, num: NUM_OPTICAL, T: 0 },
  };

  /** State (model, params, numerics, T, view range) for an Explore system or example. */
  function stateFor(key) {
    const ex = EXAMPLES[key], sys = SYSTEMS[key] || (ex && ex.system ? SYSTEMS[ex.system] : null) || ex;
    if (!sys) throw new Error("Unknown system " + key);
    return {
      model: sys.model, params: Object.assign(D.DEFAULT_PARAMS(sys.model), sys.params || {}), num: Object.assign({}, sys.num || NUM_OPTICAL),
      T: sys.T || 0, range: sys.range ? { w1: sys.range.map(cm), w3: sys.range.map(cm) } : null,
    };
  }

  /** Advanced-mode sections for every model parameter. */
  const SECTIONS = [
    { key: "electronic", title: "Electronic system", keys: ["w", "A", "e1", "e2", "J", "mu1", "mu2", "phi", "bx", "tau2", "tauC", "T1"] },
    { key: "vibrational", title: "Vibrational dynamics", keys: ["wv", "S", "tauV", "modes"] },
    { key: "environment", title: "Environment and disorder", keys: ["inh", "sigma", "temp"] },
    { key: "pathways", title: "Nonlinear pathways", keys: ["anh", "kappa", "esaSign", "tauF", "pol"] },
  ];
  const sectionOf = (k) => (SECTIONS.find((s) => s.keys.includes(k)) || SECTIONS[0]).key;

  /**
   * Which model implements which physical effect. "—" means not implemented anywhere: such effects are listed so
   * that users know they are absent, never shown as controls.
   */
  const CAPABILITIES = [
    ["Transition energies", { tls: "ω_eg", ladder: "ω_eg, ω_fe", dimer: "ε₁, ε₂ → excitons", vibronic: "ω₀₀", brownian: "ω_eg" }],
    ["Transition dipoles", { tls: "amplitude A", ladder: "|μ_fe|²/|μ_eg|²", dimer: "|μ₁|, |μ₂|, angle; orientational average" }],
    ["Excitonic coupling", { dimer: "J" }],
    ["Population relaxation (T₁)", { tls: "✓", ladder: "✓", dimer: "✓ (no transfer between excitons)", vibronic: "✓", brownian: "optional" }],
    ["Homogeneous (Bloch) dephasing", { tls: "✓", ladder: "✓", dimer: "✓", vibronic: "✓", brownian: "optional extra" }],
    ["Vibrational frequency", { vibronic: "ω_v", brownian: "underdamped oscillators" }],
    ["Huang–Rhys factor", { vibronic: "S", brownian: "S = λ/ω per oscillator" }],
    ["Vibrational dephasing", { vibronic: "damping time", brownian: "γ per oscillator" }],
    ["Reorganisation energy", { vibronic: "λ = Sω_v (implicit)", brownian: "λ per oscillator" }],
    ["Static (inhomogeneous) disorder", { tls: "✓", ladder: "✓", dimer: "✓ (fully correlated)", vibronic: "✓", brownian: "σ" }],
    ["Spectral diffusion", { brownian: "✓ via critically or over-damped oscillators" }],
    ["Frequency-correlation time", { brownian: "set by ω and γ of the damped oscillators" }],
    ["Temperature", { vibronic: "thermal vibrational population", brownian: "quantum bath at T" }],
    ["GSB and SE pathways", { tls: "✓", ladder: "✓", dimer: "✓", vibronic: "✓", brownian: "✓" }],
    ["ESA pathways", { ladder: "✓ (sign selectable)", dimer: "✓ (two-exciton state)" }],
    ["Electronic coherence during T", { dimer: "✓" }],
    ["Vibrational coherence during T", { vibronic: "✓", brownian: "✓" }],
    ["Population / coherence transfer", {}],
    ["Double-quantum coherence (DQC)", {}],
    ["Finite pulse durations", {}],
  ];

  /** Carrier at the centre of a model's transitions (the "auto" rotating frame). */
  function autoCarrier(model, p) {
    const M = D.MODELS[model];
    if (model === "vibronic") return p.w + p.S * p.wv;
    if (model === "brownian") return p.w;
    const tr = M.transitions.call(M, p).filter((w) => isFinite(w));
    return (Math.min(...tr) + Math.max(...tr)) / 2;
  }
  /** Engine job from interface state {source, model, params, comps, expr, data, T, Tindex, num}. */
  function buildJob(st) {
    const num = Object.assign({}, st.num);
    if (st.source === "model" && num.autoCarrier) num.carrier = autoCarrier(st.model, st.params);
    const job = { source: st.source, T: st.T, numerics: num };
    if (st.source === "model") Object.assign(job, { model: st.model, params: st.params, comps: st.comps });
    else if (st.source === "expr") job.expr = st.expr;
    else { job.data = st.data; job.Tindex = st.Tindex || 0; }
    return job;
  }

  return { SYSTEMS, EXAMPLES, stateFor, SECTIONS, sectionOf, CAPABILITIES, autoCarrier, buildJob, NUM_OPTICAL, NUM_PAPER };
});
