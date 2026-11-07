// Tests for the v2 interface model (assets/js/twodes-ui.js) and the analysis helpers it relies on.
// Run: node --test tests/*.test.js
const test = require("node:test");
const assert = require("node:assert");
const D = require("../assets/js/twodes-core.js");
const UI = require("../assets/js/twodes-ui.js");

const near = (got, exp, tol, msg) => assert.ok(Math.abs(got - exp) <= tol, `${msg ? msg + ": " : ""}${got} != ${exp} ± ${tol}`);
const jobOf = (key, over) => { const st = UI.stateFor(key); return UI.buildJob(Object.assign({ source: "model", comps: { GSB: true, SE: true, ESA: true } }, st, over || {})); };

test("Explore presets are the v1 examples, with the documented parameters", () => {
  const p = (key) => UI.stateFor(key).params;
  assert.deepStrictEqual(p("tls"), D.DEFAULT_PARAMS("tls"));
  assert.deepStrictEqual(p("uncoupled"), Object.assign(D.DEFAULT_PARAMS("dimer"), { J: 0, bx: 0 }));
  assert.deepStrictEqual(p("coupled"), Object.assign(D.DEFAULT_PARAMS("dimer"), { J: D.cm(300) }));
  assert.deepStrictEqual(p("vibronic"), D.DEFAULT_PARAMS("vibronic"));
  assert.deepStrictEqual(p("brownian"), D.DEFAULT_PARAMS("brownian"));
  // the paper preset keeps the paper's numerics and window
  const b = UI.stateFor("brownian"); assert.strictEqual(b.num.dt, 10); assert.strictEqual(b.T, 50); near(D.UNITS.cm.to(b.range.w1[0]), 11800, 1e-6);
  // examples resolve to the same states as the systems they name
  assert.deepStrictEqual(UI.stateFor("ex3"), UI.stateFor("coupled")); assert.deepStrictEqual(UI.stateFor("sahu"), UI.stateFor("brownian"));
  assert.strictEqual(UI.stateFor("ex5").model, "ladder");
});

test("the interface builds exactly the engine job of v1 for identical inputs", () => {
  // v1 built jobs as {source, T, numerics (carrier from the centre of the transitions), model, params, comps}
  const st = UI.stateFor("coupled"), comps = { GSB: true, SE: true, ESA: true };
  const job = UI.buildJob(Object.assign({ source: "model", comps }, st));
  const X = D.MODELS.dimer.exciton(st.params), tr = [X.E[0], X.E[1], X.Ef - X.E[0], X.Ef - X.E[1]];
  near(job.numerics.carrier, (Math.min(...tr) + Math.max(...tr)) / 2, 1e-15);
  assert.deepStrictEqual(Object.keys(job).sort(), ["T", "comps", "model", "numerics", "params", "source"]);
  // and the spectra are bit-for-bit those of the engine run on the same job
  const a = D.run(job), b = D.run(JSON.parse(JSON.stringify(job)));
  assert.deepStrictEqual(Array.from(a.total.re), Array.from(b.total.re));
  // formula and data sources pass their inputs through unchanged
  const ej = UI.buildJob({ source: "expr", expr: { params: "", R: "1", NR: "1" }, T: 5, num: UI.NUM_OPTICAL });
  assert.deepStrictEqual(ej.expr, { params: "", R: "1", NR: "1" }); assert.strictEqual(ej.numerics.carrier, UI.NUM_OPTICAL.carrier);
});

test("Explore controls map onto model parameters and back", () => {
  for (const key of Object.keys(UI.SYSTEMS)) {
    const sys = UI.SYSTEMS[key], p = UI.stateFor(key).params;
    for (const c of ["energy", "coupling", "dephasing"]) {
      const ctl = sys[c];
      if (!ctl) { assert.ok(sys[c + "Why"], `${key}: ${c} has no reason for being absent`); continue; }
      const q = JSON.parse(JSON.stringify(p)), v = ctl.get(q) * 1.1 + (ctl.get(q) === 0 ? 0.01 : 0);
      ctl.set(q, v); near(ctl.get(q), v, 1e-12 * Math.abs(v) + 1e-15, `${key} ${c}`);
    }
  }
  // the pair keeps its splitting when its centre moves; Γ = 1/T₂
  const p = UI.stateFor("coupled").params; UI.SYSTEMS.coupled.energy.set(p, D.cm(14000));
  near(D.UNITS.cm.to(p.e2 - p.e1), 600, 1e-6); near(D.UNITS.cm.to((p.e1 + p.e2) / 2), 14000, 1e-6);
  const q = UI.stateFor("tls").params; UI.SYSTEMS.tls.dephasing.set(q, D.cm(100)); near(q.tau2, 1 / D.cm(100), 1e-12);
});

test("every model parameter belongs to an Advanced section; capabilities never claim an absent model", () => {
  for (const [k, M] of Object.entries(D.MODELS)) for (const q of M.params) assert.ok(UI.SECTIONS.some((s) => s.keys.includes(q.k)), `${k}.${q.k} has no section`);
  for (const [effect, by] of UI.CAPABILITIES) for (const m of Object.keys(by)) assert.ok(D.MODELS[m], `${effect}: unknown model ${m}`);
  // unsupported effects are listed with no implementing model
  for (const name of ["Double-quantum coherence (DQC)", "Finite pulse durations", "Population / coherence transfer"]) assert.deepStrictEqual(UI.CAPABILITIES.find((r) => r[0] === name)[1], {});
  // ESA is claimed only where the engine has ESA pathways
  assert.deepStrictEqual(Object.keys(UI.CAPABILITIES.find((r) => r[0] === "ESA pathways")[1]).sort(), Object.entries(D.MODEL_PATHWAYS).filter(([, v]) => v.includes("ESA")).map(([k]) => k).sort());
});

test("pathway decomposition: parts add up exactly to the total, with signs and phases intact", () => {
  for (const key of ["tls", "coupled", "vibronic"]) {
    const d = D.decompose(Object.assign(jobOf(key), { T: 37 }));
    assert.ok(d.maxDev < 1e-12, `${key}: Σ parts − total = ${d.maxDev}`);
  }
  const lad = D.decompose(jobOf("ex5")); assert.deepStrictEqual(lad.kinds, ["GSB", "SE", "ESA"]);
  // ESA is negative in the absorptive spectrum, GSB and SE positive (the convention of the engine)
  const S = lad.parts.ESA.total, a = D.nearest(S.w1, D.cm(15000)), b = D.nearest(S.w3, D.cm(14600));
  assert.ok(S.re[a * S.N3 + b] < 0); assert.ok(lad.parts.GSB.total.re[a * S.N3 + a] > 0);
  assert.throws(() => D.decompose({ source: "expr" }), /built-in model/);
});

test("difference spectra need identical grids", () => {
  const A = D.run(jobOf("tls")), B = D.run(Object.assign(jobOf("tls"), { T: 200 }));
  const d = D.difference(A.total, B.total); near(d.re[100], A.total.re[100] - B.total.re[100], 1e-15);
  const C = D.run(Object.assign(jobOf("tls"), { numerics: Object.assign({}, jobOf("tls").numerics, { dt: 5 }) }));
  assert.match(D.incompatibility(A.total, C.total), /Different ω₁ axes|Different grids/);
  assert.throws(() => D.difference(A.total, C.total));
});

test("cuts, direct peak measures and fits recover a known Lorentzian", () => {
  // two-level system, absorptive spectrum: an exact product of Lorentzians with FWHM 2Γ on each axis
  const job = jobOf("tls", { num: Object.assign({}, UI.NUM_OPTICAL, { n1: 256, n3: 256, pad: 4 }) }); job.numerics = Object.assign({}, job.numerics, { n1: 256, n3: 256, pad: 4 });
  const r = D.run(job), S = r.total, w = D.cm(15000), G = 1 / 50;
  for (const kind of ["h", "v"]) {
    const c = D.cut(S, w, w, kind, 0.15), st = D.peakStats(c.x, c.re);
    near(st.fwhm, 2 * G, 0.02 * 2 * G, `${kind} direct FWHM`); near(st.x0, 0, st.dxGrid, `${kind} position`);
    const f = D.fitPeak(c.x, c.re, "lorentz");
    near(f.fwhm, 2 * G, 0.02 * 2 * G, `${kind} fitted FWHM`); near(f.x0, 0, 0.002, `${kind} fitted position`); assert.ok(f.r2 > 0.999);
  }
  // along the anti-diagonal the projected width of a round Lorentzian product is narrower (≈ 0.64 × 2Γ for Lorentz × Lorentz)
  const an = D.peakStats(...((c) => [c.x, c.re])(D.cut(S, w, w, "anti", 0.15)));
  assert.ok(an.fwhm < 2 * G && an.fwhm > 0.5 * 2 * G);
  // synthetic Gaussian: fit returns the exact parameters and the analytic area
  const xs = Array.from({ length: 121 }, (_, i) => -3 + i * 0.05), ys = xs.map((x) => 2 * Math.exp((-4 * Math.LN2 * (x - 0.3) ** 2) / 0.8 ** 2) + 0.1);
  const g = D.fitPeak(xs, ys, "gauss"); near(g.fwhm, 0.8, 1e-4); near(g.x0, 0.3, 1e-4); near(g.amp, 2, 1e-4); near(g.area, 2 * 0.8 * Math.sqrt(Math.PI / (4 * Math.LN2)), 1e-3);
  // warnings: a coarse grid and a peak at the edge are flagged
  assert.ok(D.peakStats([0, 1, 2, 3, 4], [0, 0.4, 1, 0.4, 0]).warnings.some((w) => /grid points across the FWHM/.test(w)));
  assert.ok(D.peakStats([0, 1, 2, 3], [1, 0.8, 0.5, 0.2]).warnings.some((w) => /edge/.test(w)));
});
