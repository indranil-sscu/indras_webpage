// Generates tests/fixtures/twodes-golden.json from the 2DES engine. It was run once on the v1.1 engine
// (main @ 099b204) before the v2 interface redesign; tests/twodes-regression.test.js checks that the engine
// still reproduces these numbers. Re-run only after an intentional, validated change to the numerics.
const D = require("../assets/js/twodes-core.js");
const fs = require("fs");
const cm = D.cm;
const base = { dt: 4, n1: 64, n3: 64, pad: 2, window: "none", broaden: 0, halfFirst: true, convention: "standard", carrier: cm(15000) };
const J = (model, params, T, num, comps) => ({ source: "model", model, params: Object.assign(D.DEFAULT_PARAMS(model), params || {}), T: T || 0, comps, numerics: Object.assign({}, base, num || {}) });
const cases = {
  tls_default: J("tls"),
  tls_T300_inh: J("tls", { inh: cm(300), tau2: 80 }, 300),
  tls_hann_pad4: J("tls", {}, 0, { window: "hann", pad: 4, broaden: cm(20) }),
  tls_unflipped: J("tls", {}, 0, { convention: "unflipped" }),
  tls_conjugate: J("tls", {}, 0, { convention: "conjugate" }),
  tls_nohalf: J("tls", {}, 0, { halfFirst: false }),
  ladder_default: J("ladder", {}, 100),
  dimer_J0: J("dimer", { J: 0 }, 50),
  dimer_J300: J("dimer", { J: cm(300) }, 70),
  dimer_zzyy_noESA: J("dimer", { J: cm(200), pol: "zzyy" }, 40, {}, { GSB: true, SE: true, ESA: false }),
  vibronic_default: J("vibronic", {}, 60, { carrier: cm(15300) }),
  vibronic_hot: J("vibronic", { temp: 300, S: 0.8 }, 25, { carrier: cm(15300) }),
  brownian_paper: J("brownian", {}, 50, { dt: 10, n1: 64, n3: 64, carrier: cm(12870) }),
  brownian_GSBonly: J("brownian", {}, 120, { dt: 10, n1: 64, n3: 64, carrier: cm(12870) }, { GSB: true, SE: false }),
  expr_tls: { source: "expr", T: 150, numerics: Object.assign({}, base), expr: { params: "A = 1\nw = cm(15000)\nG = 1/50\nT1 = 1000", R: "2*A*exp(i*w*t1 - G*t1)*exp(-i*w*t3 - G*t3)*exp(-T/T1)", NR: "2*A*exp(-i*w*t1 - G*t1)*exp(-i*w*t3 - G*t3)*exp(-T/T1)" } },
};
/** A compact fingerprint of a spectrum: sums, norm and a few exact samples. */
function print(S) {
  if (!S) return null;
  let sr = 0, si = 0, n2 = 0;
  for (let j = 0; j < S.re.length; j++) { sr += S.re[j]; si += S.im[j]; n2 += S.re[j] ** 2 + S.im[j] ** 2; }
  const idx = [0, 7, Math.floor(S.re.length / 3), Math.floor(S.re.length / 2) + 5, S.re.length - 1];
  return { N1: S.N1, N3: S.N3, w1: [S.w1[0], S.w1[S.N1 - 1]], w3: [S.w3[0], S.w3[S.N3 - 1]], sr, si, n2, samples: idx.map((j) => [j, S.re[j], S.im[j]]) };
}
const out = { engine: "twodes-core v1.1 (main @ 099b204)", cases: {} };
for (const [k, job] of Object.entries(cases)) { const r = D.run(job); out.cases[k] = { job, R: print(r.R), NR: print(r.NR), total: print(r.total), warnings: r.warnings }; }
// coherence-map fingerprint (small region)
const cmJob = J("brownian", {}, 0, { dt: 10, n1: 64, n3: 64, carrier: cm(12870) });
const info = D.cmPrepare(cmJob, { T0: 10, T1: 400, dT: 10, region: { w1: [cm(12400), cm(13200)], w3: [cm(12400), cm(13200)] }, pad: 2 });
const M = D.cmMaps(cm(350)), sum = (a) => a.reduce((s, x) => s + x, 0);
out.cm = { job: cmJob, opts: { T0: 10, T1: 400, dT: 10, region: { w1: [cm(12400), cm(13200)], w3: [cm(12400), cm(13200)] }, pad: 2 }, wT: cm(350), tau: info.tau, real: sum(M.real), pos: sum(M.pos), neg: sum(M.neg) };
fs.writeFileSync(__dirname + "/fixtures/twodes-golden.json", JSON.stringify(out, null, 1));
console.log("wrote", Object.keys(out.cases).length, "cases");
