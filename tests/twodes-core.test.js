// Validation tests for assets/js/twodes-core.js (2DES simulator). Run: node --test tests/*.test.js
const test = require("node:test");
const assert = require("node:assert");
const D = require("../assets/js/twodes-core.js");

const near = (got, exp, tol, msg) => assert.ok(Math.abs(got - exp) <= tol, `${msg ? msg + ": " : ""}${got} != ${exp} ± ${tol}`);
const NUM = { dt: 4, n1: 128, n3: 128, pad: 2, window: "none", broaden: 0, halfFirst: true, convention: "standard", carrier: D.cm(15000) };
const job = (model, over, num) => ({ source: "model", model, params: Object.assign(D.DEFAULT_PARAMS(model), over || {}), T: 0, numerics: Object.assign({}, NUM, num || {}) });
const absorptive = (r) => Float64Array.from(r.total.re);
function argmax2(S, arr) { let b = 0; for (let j = 1; j < arr.length; j++) if (arr[j] > arr[b]) b = j; return { a: Math.floor(b / S.N3), b: b % S.N3, v: arr[b] }; }
const sliceRow = (S, arr, a) => Array.from({ length: S.N3 }, (_, b) => arr[a * S.N3 + b]);
const sliceCol = (S, arr, b) => Array.from({ length: S.N1 }, (_, a) => arr[a * S.N3 + b]);

test("FFT agrees with the direct DFT for both signs", () => {
  const n = 64, re = Float64Array.from({ length: n }, (_, i) => Math.sin(i * 0.7) + 0.1 * i), im = Float64Array.from({ length: n }, (_, i) => Math.cos(i * 1.3));
  for (const sign of [-1, 1]) {
    const r = re.slice(), m = im.slice(); D.fft(r, m, sign); const d = D.dft(re, im, sign);
    for (let k = 0; k < n; k++) { near(r[k], d.re[k], 1e-10); near(m[k], d.im[k], 1e-10); }
  }
});

test("1. an isolated transition gives one peak at (ω_eg, ω_eg)", () => {
  for (const w of [14500, 15000, 15800]) {
    const r = D.run(job("tls", { w: D.cm(w) })), S = r.total, A = absorptive(r), m = argmax2(S, A);
    near(D.UNITS.cm.to(S.w1[m.a]), w, 20, "ω1"); near(D.UNITS.cm.to(S.w3[m.b]), w, 20, "ω3");
    assert.ok(m.v > 0, "GSB + SE peak is positive");
  }
});

test("2. stronger homogeneous dephasing broadens the peak: FWHM = 2Γ", () => {
  const widths = [];
  for (const tau2 of [30, 60, 120]) {
    const r = D.run(job("tls", { tau2 }, { n1: 256, n3: 256, pad: 4 })), S = r.total, A = absorptive(r), m = argmax2(S, A);
    const wf = D.fwhm(S.w3, sliceRow(S, A, m.a)), wf1 = D.fwhm(S.w1, sliceCol(S, A, m.b));
    near(wf, 2 / tau2, 0.02 * (2 / tau2), "ω3 width"); near(wf1, 2 / tau2, 0.02 * (2 / tau2), "ω1 width");
    widths.push(wf);
  }
  assert.ok(widths[0] > widths[1] && widths[1] > widths[2]);
});

test("3. population relaxation reduces the signal with waiting time as exp(−T/T1)", () => {
  const base = job("tls", { T1: 400 });
  const peak = (T) => { const r = D.run(Object.assign({}, base, { T })); return argmax2(r.total, absorptive(r)).v; };
  const p0 = peak(0);
  for (const T of [200, 400, 800]) near(peak(T) / p0, Math.exp(-T / 400), 1e-9);
});

test("4. rephasing vs non-rephasing frequency axes under each convention", () => {
  const w = D.cm(15200);
  for (const [conv, sR, sN, s3] of [["standard", 1, 1, 1], ["unflipped", -1, 1, 1], ["conjugate", -1, -1, -1]]) {
    const r = D.run(job("tls", { w }, { convention: conv }));
    for (const [S, s1] of [[r.R, sR], [r.NR, sN]]) {
      const mag = Float64Array.from(S.re, (v, j) => Math.hypot(v, S.im[j])), m = argmax2(S, mag);
      near(S.w1[m.a], s1 * w, 0.004, `${conv} ω1`); near(S.w3[m.b], s3 * w, 0.004, `${conv} ω3`);
    }
  }
  // with static disorder the rephasing peak is elongated along the diagonal, the non-rephasing one is not
  const r = D.run(job("tls", { inh: D.cm(400), tau2: 100 }, { n1: 256, n3: 256 }));
  const shape = (S) => {
    const mag = Float64Array.from(S.re, (v, j) => Math.hypot(v, S.im[j])), m = argmax2(S, mag);
    let sxy = 0, sxx = 0, syy = 0;
    for (let a = 0; a < S.N1; a++) for (let b = 0; b < S.N3; b++) { const v = mag[a * S.N3 + b]; if (v < 0.3 * m.v) continue; const x = S.w3[b] - S.w3[m.b], y = S.w1[a] - S.w1[m.a]; sxy += v * x * y; sxx += v * x * x; syy += v * y * y; }
    return sxy / Math.sqrt(sxx * syy);
  };
  assert.ok(shape(r.R) > 0.5, "rephasing: positive correlation (diagonal elongation)");
  assert.ok(shape(r.NR) < shape(r.R) - 0.3, "non-rephasing is much less diagonal");
});

test("5. the absorptive spectrum is the product of absorptive Lorentzians", () => {
  const tau2 = 50, G = 1 / tau2, w = D.cm(15000);
  const r = D.run(job("tls", { tau2 }, { n1: 512, n3: 512, dt: 2, pad: 1 })), S = r.total, A = absorptive(r);
  // analytic: Re[S_R + S_NR] = 2A·2Γ/(Γ² + Δ1²)·Γ/(Γ² + Δ3²) ×2 (GSB + SE), A = 1
  let worst = 0, top = 0;
  for (let a = 0; a < S.N1; a += 7) for (let b = 0; b < S.N3; b += 7) {
    const d1 = S.w1[a] - w, d3 = S.w3[b] - w;
    if (Math.abs(d1) > 0.15 || Math.abs(d3) > 0.15) continue;
    const ex = 2 * ((2 * G) / (G * G + d1 * d1)) * (G / (G * G + d3 * d3));
    worst = Math.max(worst, Math.abs(A[a * S.N3 + b] - ex)); top = Math.max(top, ex);
  }
  assert.ok(worst / top < 2e-3, `max deviation ${(worst / top * 100).toFixed(3)}% of the peak`);
});

test("6. zero padding does not change the physical linewidth", () => {
  const ws = [1, 2, 4, 8].map((pad) => { const r = D.run(job("tls", { tau2: 40 }, { pad })), S = r.total, A = absorptive(r), m = argmax2(S, A); return D.fwhm(S.w3, sliceRow(S, A, m.a)); });
  // pad 1 samples the peak coarsely; from pad 2 on the interpolated width is stable and equals 2Γ
  for (const w of ws.slice(1)) near(w, 2 / 40, 0.03 * (2 / 40));
  near(ws[1], ws[3], 0.01 * ws[3]);
  // the grid gets finer, the line does not get narrower
  const grids = [1, 2, 4, 8].map((pad) => D.run(job("tls", {}, { pad })).info.gridSpacing3);
  near(grids[0] / grids[3], 8, 1e-12);
});

test("7. FFT spectrum equals a direct numerical Fourier integration (small dataset)", () => {
  const n1 = 12, n3 = 10, dt = 3, x = { re: new Float64Array(n1 * n3), im: new Float64Array(n1 * n3) };
  for (let j = 0; j < n1 * n3; j++) { x.re[j] = Math.sin(1.7 * j) + 0.3; x.im[j] = Math.cos(0.9 * j * j); }
  for (const o of [{ sign1: -1, sign3: 1, pad: 2, window: "hann", broaden: 0.01, halfFirst: true }, { sign1: 1, sign3: -1, pad: 1, window: "none", broaden: 0, halfFirst: false }]) {
    const A = D.spectrum2D(x, n1, n3, dt, o);
    const ts1 = Float64Array.from({ length: n1 }, (_, i) => i * dt), ts3 = Float64Array.from({ length: n3 }, (_, i) => i * dt);
    const B = D.spectrum2DDirect(x, ts1, ts3, dt, o);
    // brute force at a few points: S = Σ w x exp(i s1 ω1 t1 + i s3 ω3 t3)
    const w1 = D.taper(ts1, o), w3 = D.taper(ts3, o), q1 = D.trapezoid(ts1, o.halfFirst), q3 = D.trapezoid(ts3, o.halfFirst);
    for (const [a, b] of [[0, 0], [3, 5], [A.N1 - 1, 2], [A.N1 / 2, A.N3 / 2]]) {
      let sr = 0, si = 0;
      for (let i = 0; i < n1; i++) for (let k = 0; k < n3; k++) {
        const ph = o.sign1 * A.w1[a] * ts1[i] + o.sign3 * A.w3[b] * ts3[k], f = w1[i] * w3[k] * q1[i] * q3[k], j = i * n3 + k;
        sr += f * (x.re[j] * Math.cos(ph) - x.im[j] * Math.sin(ph)); si += f * (x.re[j] * Math.sin(ph) + x.im[j] * Math.cos(ph));
      }
      near(A.re[a * A.N3 + b], sr, 1e-10); near(A.im[a * A.N3 + b], si, 1e-10);
    }
    for (let j = 0; j < A.re.length; j++) { near(A.re[j], B.re[j], 1e-10); near(A.im[j], B.im[j], 1e-10); }
  }
});

test("8. flipping a sign convention mirrors the corresponding frequency axis", () => {
  const p = { w: D.cm(15300), tau2: 45 }, st = D.run(job("tls", p)), un = D.run(job("tls", p, { convention: "unflipped" })), cj = D.run(job("tls", p, { convention: "conjugate" }));
  const N1 = st.R.N1, N3 = st.R.N3;
  for (let a = 0; a < N1; a++) {
    const am = (N1 - a) % N1;
    if (a > 0) near(un.R.w1[am], -st.R.w1[a], 1e-12); // a = 0 is the Nyquist point, its own mirror on the periodic grid
    for (let b = 0; b < N3; b += 5) {
      const bm = (N3 - b) % N3;
      near(un.R.re[am * N3 + b], st.R.re[a * N3 + b], 1e-9); near(un.R.im[am * N3 + b], st.R.im[a * N3 + b], 1e-9);
      // conjugate convention (all kernel signs flipped): both axes mirrored, values unchanged: S_conj(−ω1, −ω3) = S_std(ω1, ω3)
      near(cj.NR.re[am * N3 + bm], st.NR.re[a * N3 + b], 1e-9); near(cj.NR.im[am * N3 + bm], st.NR.im[a * N3 + b], 1e-9);
    }
  }
  // the absorptive spectrum is the same whichever way the rephasing ω1 axis is presented
  for (let j = 0; j < st.total.re.length; j += 11) near(un.total.re[j], st.total.re[j], 1e-9);
});

test("convergence check: refining dt and extending t_max change a well-sampled spectrum very little", () => {
  // within ±0.15 rad/fs (±800 cm⁻¹) of the line; the residual sampling change is the aliasing of the 1/Ω tails of the one-sided transform
  const w = D.cm(15000), range = { w1: [w - 0.15, w + 0.15], w3: [w - 0.15, w + 0.15] };
  const c = D.convergence(job("tls", { tau2: 50 }, { n1: 128, n3: 128, dt: 4 }), range);
  assert.ok(c.sampling < 5e-3, `sampling ${c.sampling}`); assert.ok(c.window < 1e-3, `window ${c.window}`);
  const c2 = D.convergence(job("tls", { tau2: 50 }, { n1: 256, n3: 256, dt: 2 }), range);
  assert.ok(c2.sampling < c.sampling / 3, `halving dt reduces the sampling error (${c2.sampling} vs ${c.sampling})`);
  const bad = D.convergence(job("tls", { tau2: 300 }, { n1: 32, n3: 32, dt: 4 }));
  assert.ok(bad.window > 0.05, `a truncated response is flagged (${bad.window})`);
});

test("rotating frame: the spectrum does not depend on the carrier", () => {
  const a = D.run(job("tls", { w: D.cm(15100) }, { carrier: D.cm(15000), dt: 2, n1: 256, n3: 256 }));
  const b = D.run(job("tls", { w: D.cm(15100) }, { carrier: D.cm(15200), dt: 2, n1: 256, n3: 256 }));
  // compare at the peak (grids are offset by the carrier difference, so interpolate via the analytic peak value)
  const pa = argmax2(a.total, absorptive(a)), pb = argmax2(b.total, absorptive(b));
  near(a.total.w3[pa.b], b.total.w3[pb.b], 0.01); near(pa.v, pb.v, 0.02 * pa.v);
  assert.ok(D.run(job("tls", { w: D.cm(15000) }, { carrier: D.cm(15000), dt: 8 })).warnings.length === 0);
  assert.ok(D.run(job("tls", { w: D.cm(19500) }, { carrier: D.cm(15000), dt: 4 })).warnings.some((w) => /Aliasing/.test(w))); // band edge at ±4170 cm⁻¹
});

test("ladder: ESA appears with its sign at (ω_eg, ω_eg − Δ); harmonic limit cancels", () => {
  const p = { anh: D.cm(500), kappa: 1.5, esaSign: -1 };
  const r = D.run(job("ladder", p)), S = r.total, A = absorptive(r);
  const a = D.nearest(S.w1, D.cm(15000)), b = D.nearest(S.w3, D.cm(14500));
  assert.ok(A[a * S.N3 + b] < 0, "ESA negative");
  const pos = D.run(job("ladder", Object.assign({}, p, { esaSign: 1 }))); assert.ok(absorptive(pos)[a * S.N3 + b] > 0, "sign is configurable");
  const h = D.run(job("ladder", { anh: 0, kappa: 2, tauF: 60, tau2: 60 }));
  let m = 0; for (const v of absorptive(h)) m = Math.max(m, Math.abs(v)); assert.ok(m < 1e-9, `harmonic ladder cancels (${m})`);
});

test("dimer: uncoupled transitions have no cross peaks; coupling creates them", () => {
  const cross = (r) => { const S = r.total, A = absorptive(r); const a = D.nearest(S.w1, D.cm(15000)), b = D.nearest(S.w3, D.cm(15600)); return A[a * S.N3 + b]; };
  const diag = (r) => { const S = r.total, A = absorptive(r); const a = D.nearest(S.w1, D.cm(15000)); return A[a * S.N3 + a]; };
  // J = 0: the dimer spectrum is exactly the sum of two independent two-level systems (all cross-peak pathways cancel),
  // at every waiting time and for both polarizations (orientation factor 1 for ZZZZ, 1/3 for ZZYY, for each isolated transition)
  for (const pol of ["zzzz", "zzyy"]) for (const T of [0, 37, 120]) {
    const u = D.run(Object.assign(job("dimer", { J: 0, pol }), { T })), f = pol === "zzzz" ? 1 : 1 / 3;
    const a = D.run(Object.assign(job("tls", { w: D.cm(15000), tau2: 60, A: f }), { T })), b = D.run(Object.assign(job("tls", { w: D.cm(15600), tau2: 60, A: f }), { T }));
    let worst = 0, top = 0;
    for (let j = 0; j < u.total.re.length; j++) { worst = Math.max(worst, Math.abs(u.total.re[j] - a.total.re[j] - b.total.re[j]), Math.abs(u.total.im[j] - a.total.im[j] - b.total.im[j])); top = Math.max(top, Math.abs(u.total.re[j])); }
    assert.ok(worst < 1e-9 * top, `${pol}, T = ${T}: J = 0 dimer differs from two TLS by ${worst / top}`);
  }
  // coupling: the cross peak is no longer just the tails of the diagonal peaks
  const u = D.run(job("dimer", { J: 0 })), c = D.run(job("dimer", { J: D.cm(150) }));
  assert.ok(Math.abs(cross(c) - cross(u)) > 0.05 * diag(c), `coupled cross peak ${cross(c)} vs uncoupled ${cross(u)} (diag ${diag(c)})`);
  // exciton energies: (ε1+ε2)/2 ± sqrt(((ε1−ε2)/2)² + J²)
  const X = D.MODELS.dimer.exciton(D.DEFAULT_PARAMS("dimer")), p = D.DEFAULT_PARAMS("dimer"), m = (p.e1 + p.e2) / 2, d = Math.hypot((p.e1 - p.e2) / 2, p.J);
  near(Math.min(...X.E), m - d, 1e-12); near(Math.max(...X.E), m + d, 1e-12);
  // H-dimer: parallel identical dipoles, ε1 = ε2, J > 0: the lower exciton is dark
  const H = D.MODELS.dimer.exciton(Object.assign({}, p, { e2: p.e1, phi: 0 })), lo = H.E[0] < H.E[1] ? 0 : 1;
  near(Math.hypot(...H.mu[lo]), 0, 1e-12); near(Math.hypot(...H.mu[1 - lo]), Math.SQRT2, 1e-12);
});

test("dimer: electronic coherence beats at |E+ − E−| on the rephasing cross peak", () => {
  const p = Object.assign(D.DEFAULT_PARAMS("dimer"), { tauC: 1e6, T1: 1e6 }), X = D.MODELS.dimer.exciton(p), w = Math.abs(X.E[1] - X.E[0]);
  const j = { source: "model", model: "dimer", params: p, comps: { GSB: false, ESA: false }, T: 0, numerics: NUM };
  const lo = Math.min(...X.E), hi = Math.max(...X.E), Ts = Array.from({ length: 80 }, (_, i) => i * 2);
  const tr = D.waitingTrace(j, lo, hi, Ts).R.map((v) => v[0]);
  // the rephasing SE cross peak (α = low, β = high) oscillates as cos(ωT): find the period from zero crossings
  const zc = []; for (let i = 1; i < tr.length; i++) if (Math.sign(tr[i]) !== Math.sign(tr[i - 1])) zc.push(Ts[i - 1] + (Ts[i] - Ts[i - 1]) * tr[i - 1] / (tr[i - 1] - tr[i]));
  const period = 2 * (zc[zc.length - 1] - zc[0]) / (zc.length - 1);
  near(period, (2 * Math.PI) / w, 0.01 * (2 * Math.PI) / w);
});

test("vibronic: S = 0 reduces to the two-level system; sidebands at ω00 + ω_v; beating period 2π/ω_v", () => {
  const tl = D.run(job("tls", { tau2: 60, T1: 1000, A: 1 })), v0 = D.run(job("vibronic", { S: 0 }));
  for (let j = 0; j < tl.total.re.length; j += 13) near(v0.total.re[j], tl.total.re[j], 1e-9 * Math.abs(tl.total.re[j]) + 1e-12);
  const p = { S: 0.6, wv: D.cm(700), tauV: 2000 };
  const r = D.run(job("vibronic", p)), S = r.total, A = absorptive(r);
  const at = (x, y) => A[D.nearest(S.w1, D.cm(x)) * S.N3 + D.nearest(S.w3, D.cm(y))];
  assert.ok(at(15700, 15700) > 0.1 * at(15000, 15000), "0–1 diagonal sideband present");
  assert.ok(at(15700, 15000) > 0.1 * at(15000, 15000), "cross peak (ω00 + ω_v, ω00) present");
  // GSB-only rephasing cross peak oscillates with the vibrational period
  const j = Object.assign(job("vibronic", p), { comps: { SE: false } }), w = p.wv, Ts = Array.from({ length: 100 }, (_, i) => i * 1.5);
  const tr = D.waitingTrace(j, D.cm(15700), D.cm(15000), Ts).R.map((v) => v[0]), mean = tr.reduce((s, v) => s + v, 0) / tr.length;
  const zc = []; for (let i = 1; i < tr.length; i++) if (Math.sign(tr[i] - mean) !== Math.sign(tr[i - 1] - mean)) zc.push(Ts[i - 1] + 1.5 * (tr[i - 1] - mean) / (tr[i - 1] - tr[i]));
  near(2 * (zc[zc.length - 1] - zc[0]) / (zc.length - 1), (2 * Math.PI) / w, 0.03 * (2 * Math.PI) / w);
});

test("expression mode reproduces the built-in two-level model", () => {
  const src = { params: "A = 1\nw = cm(15000)\nG = 1/50\nT1 = 1000", R: "2*A*exp(i*w*t1 - G*t1)*exp(-i*w*t3 - G*t3)*exp(-T/T1)", NR: "2*A*exp(-i*w*t1 - G*t1)*exp(-i*w*t3 - G*t3)*exp(-T/T1)" };
  const e = D.run({ source: "expr", expr: src, T: 150, numerics: NUM }), m = D.run(Object.assign(job("tls", { tau2: 50 }), { T: 150 }));
  for (let j = 0; j < e.total.re.length; j += 7) { near(e.total.re[j], m.total.re[j], 1e-8); near(e.total.im[j], m.total.im[j], 1e-8); }
});

test("expression parser: precedence, complex arithmetic and safety", () => {
  const ev = (s, env) => D.compile(D.parseExpr(s, new Set(Object.keys(env || {}))))(env || {});
  const eq = (s, v) => { const r = ev(s); near(r[0], v, 1e-12, s); near(r[1], 0, 1e-12, s); };
  eq("2+3*4", 14); eq("-2^2", -4); eq("2^3^2", 512); eq("(1+2)*3-4/2", 7); eq("2*-3", -6);
  const e = ev("exp(i*pi)"); near(e[0], -1, 1e-15); near(e[1], 0, 1e-15);
  const s = ev("sqrt(-4)"); near(s[0], 0, 1e-15); near(s[1], 2, 1e-15);
  const c = ev("conj(1+2*i)*(1+2*i)"); near(c[0], 5, 1e-15); near(c[1], 0, 1e-15);
  near(ev("cm(1)")[0], 2 * Math.PI * 2.99792458e-5, 1e-18);
  near(ev("x/2", { x: [3, 1] })[1], 0.5, 1e-15);
  for (const bad of ["alert(1)", "constructor", "x.y", "2x", "a[0]", "exp(1", "1+", "this", "__proto__", "process.exit()", "`1`", "Function('a')"]) assert.throws(() => D.parseExpr(bad, new Set()), undefined, bad);
  assert.throws(() => D.parseParams("t1 = 3"), /reserved/);
  assert.throws(() => D.parseParams("__proto__ = 3"), /reserved/);
  for (const name of ["toString", "valueOf", "hasOwnProperty"]) assert.throws(() => D.parseExpr(name + "(1)", new Set()), /Unknown function/);
  assert.deepStrictEqual(D.parseParams("a = 2\nb = a*3 # comment").b, [6, 0]);
});

test("CSV import: round trip, validation and nonuniform grids", () => {
  const j = Object.assign(job("dimer", {}, { n1: 32, n3: 32, dt: 4 }), { T: 50 }), r = D.run(j);
  const text = D.responseCSV(r.resp, 32, 32, 4, 50, NUM.carrier), parsed = D.parseCSV(text, { timeUnit: "fs" });
  assert.deepStrictEqual(parsed.errors, []); near(parsed.meta.carrier, NUM.carrier, 1e-9);
  const ds = D.buildDataset([parsed.rows], { carrier: parsed.meta.carrier });
  assert.deepStrictEqual(ds.errors, []); assert.ok(ds.data.uniform); assert.strictEqual(ds.data.n1, 32);
  const back = D.run({ source: "data", data: ds.data, Tindex: 0, numerics: Object.assign({}, NUM, { n1: 32, n3: 32 }) });
  let top = 0; for (const v of r.total.re) top = Math.max(top, Math.abs(v));
  for (let k = 0; k < r.total.re.length; k += 3) near(back.total.re[k], r.total.re[k], 1e-7 * top); // 9 significant digits in the file
  // units: the same data in ps
  const ps = text.split("\n").map((l) => (/^(R|NR),/.test(l) ? l.split(",").map((v, i) => (i >= 1 && i <= 3 ? String(+v / 1000) : v)).join(",") : l)).join("\n");
  const dps = D.buildDataset([D.parseCSV(ps, { timeUnit: "ps" }).rows]); near(dps.data.dt, 4, 1e-9);
  // errors
  const hdr = "t1,T,t3,re,im\n";
  assert.match(D.parseCSV(hdr + "0,0,0,1,0\n0,0,4,,0\n", {}).errors.join(), /Missing values on line 3/);
  assert.match(D.parseCSV(hdr + "0,0,0,1,abc\n", {}).errors.join(), /Non-numeric/);
  assert.match(D.parseCSV("t1,t3,re,im\n0,0,1,0\n", {}).errors.join(), /missing column/);
  const grid = (f) => { let s = hdr; for (const t1 of [0, 2, 4, 6]) for (const t3 of [0, 2, 4, 6]) if (f(t1, t3)) s += `${t1},0,${t3},1,0\n`; return s; };
  assert.match(D.buildDataset([D.parseCSV(grid((a, b) => !(a === 2 && b === 4)), {}).rows]).errors.join(), /incomplete grid, 15 of 4 × 4/);
  assert.match(D.buildDataset([D.parseCSV(grid(() => true) + "0,0,0,2,0\n", {}).rows]).errors.join(), /Duplicate coordinates/);
  // signed t1: negative t1 holds the rephasing pathway; t = 0 shared
  let sg = hdr; for (const t1 of [-6, -4, -2, 0, 2, 4, 6]) for (const t3 of [0, 2, 4, 6]) sg += `${t1},0,${t3},${t1 < 0 ? 1 : 2},0\n`;
  const sd = D.buildDataset([D.parseCSV(sg, { pathway: "signedR" }).rows]);
  assert.deepStrictEqual(sd.errors, []); assert.strictEqual(sd.data.R[0].re[1 * 4], 1); assert.strictEqual(sd.data.NR[0].re[1 * 4], 2);
  // nonuniform t1: direct transform equals the uniform result where the data coincide (exactly sampled function)
  const nonu = [0, 2, 4, 7, 10, 14, 18, 24, 30, 40, 52, 66, 80, 100, 124, 150, 180, 210, 250, 300];
  const f = (t1, t3) => [Math.exp(-t1 / 40 - t3 / 40) * Math.cos(0.05 * (t1 - t3)), Math.exp(-t1 / 40 - t3 / 40) * Math.sin(0.05 * (t1 - t3))];
  let nu = hdr; for (const t1 of nonu) for (let t3 = 0; t3 <= 300; t3 += 2) { const v = f(t1, t3); nu += `${t1},0,${t3},${v[0]},${v[1]}\n`; }
  const dn = D.buildDataset([D.parseCSV(nu, { pathway: "R" }).rows]);
  assert.deepStrictEqual(dn.errors, []); assert.ok(!dn.data.uniform); assert.ok(dn.warnings.some((w) => /t1 is not uniformly sampled/.test(w)));
  const rn = D.run({ source: "data", data: dn.data, Tindex: 0, numerics: Object.assign({}, NUM, { nonuniform: "direct" }) });
  assert.ok(rn.R && !rn.NR && !rn.total); assert.ok(rn.warnings.some((w) => /direct nonuniform/.test(w)));
  const mag = Float64Array.from(rn.R.re, (v, k) => Math.hypot(v, rn.R.im[k])), m = argmax2(rn.R, mag);
  near(rn.R.w1[m.a], 0.05, 0.01); near(rn.R.w3[m.b], 0.05, 0.01); // carrier 0: peak at Ω = +0.05 rad/fs on both axes (standard kernels)
  const rs = D.run({ source: "data", data: dn.data, Tindex: 0, numerics: Object.assign({}, NUM, { nonuniform: "resample" }) });
  assert.ok(rs.warnings.some((w) => /resampled/.test(w)));
});
