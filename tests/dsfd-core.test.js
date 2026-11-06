// Unit tests for assets/js/dsfd-core.js — run with: node --test tests/*.test.js
const test = require("node:test");
const assert = require("node:assert");
const F = require("../assets/js/dsfd-core.js");

const twoLevel = F.levelScheme(["e"], []);
const threeLevel = F.levelScheme(["e"], ["f"]);
const dimer = F.levelScheme(["a", "b"], ["f"]);
const R = F.PRESETS.rephasing.signs, NR = F.PRESETS.nonrephasing.signs, DQ = F.PRESETS.dqc.signs;
const near = (got, exp, tol) => assert.ok(Math.abs(got - exp) <= tol, `${got} != ${exp} ± ${tol}`);
const types = (ds) => ds.map((d) => d.type).sort().join(",");

test("linear absorption: one pathway |e⟩⟨g|", () => {
  const ds = F.enumerate([+1], twoLevel);
  assert.strictEqual(ds.length, 1);
  assert.deepStrictEqual(ds[0].rho[1], ["e", "g"]);
  assert.strictEqual(ds[0].sign, +1);
});

test("two-level system: GSB + SE in rephasing and non-rephasing, no DQC", () => {
  assert.strictEqual(types(F.enumerate(R, twoLevel)), "GSB,SE");
  assert.strictEqual(types(F.enumerate(NR, twoLevel)), "GSB,SE");
  assert.strictEqual(F.enumerate(DQ, twoLevel).length, 0);
});

test("three-level ladder: GSB, SE, ESA; ESA has opposite sign", () => {
  const ds = F.enumerate(R, threeLevel);
  assert.strictEqual(types(ds), "ESA,GSB,SE");
  const sign = Object.fromEntries(ds.map((d) => [d.type, d.sign]));
  assert.deepStrictEqual(sign, { GSB: 1, SE: 1, ESA: -1 });
  assert.ok(ds.every((d) => d.rephasing === true));
  assert.ok(F.enumerate(NR, threeLevel).every((d) => d.rephasing === false));
});

test("double-quantum coherence: two pathways through |f⟩⟨g| with opposite signs", () => {
  const ds = F.enumerate(DQ, threeLevel);
  assert.strictEqual(ds.length, 2);
  assert.ok(ds.every((d) => d.type === "DQC" && d.rho[2][0] === "f" && d.rho[2][1] === "g"));
  assert.deepStrictEqual(ds.map((d) => d.sign).sort(), [-1, 1]);
});

test("excitonic dimer (a, b, f): 4 GSB + 4 SE + 4 ESA, coherences during t2", () => {
  const ds = F.enumerate(R, dimer);
  const count = (t) => ds.filter((d) => d.type === t).length;
  assert.deepStrictEqual([count("GSB"), count("SE"), count("ESA")], [4, 4, 4]);
  const t2coh = ds.filter((d) => d.type === "SE" && d.rho[2][0] !== d.rho[2][1]);
  assert.strictEqual(t2coh.length, 2); // |a⟩⟨b| and |b⟩⟨a| during the waiting time
});

test("every pathway ends in a population and obeys phase matching", () => {
  for (const scheme of [twoLevel, threeLevel, dimer]) {
    for (const signs of [R, NR, DQ, [+1, -1, +1, -1, +1]]) {
      for (const d of F.enumerate(signs, scheme)) {
        const last = d.emission;
        assert.strictEqual(last.to, d.rho[signs.length][1]); // ket lands on the bra state
        d.interactions.forEach((s, j) => {
          const k = (s.side === "ket") === (s.kind === "up") ? +1 : -1;
          assert.strictEqual(k, signs[j]);
        });
      }
    }
  }
});

test("dipole product and direction label", () => {
  const esa = F.enumerate(R, threeLevel).find((d) => d.type === "ESA");
  assert.deepStrictEqual(F.dipoles(esa, threeLevel), ["ge", "ge", "ef", "ef"]);
  assert.strictEqual(F.direction(R), "−k₁ + k₂ + k₃");
});

test("second order: SFG and DFG pathways end in a population and carry χ(2) labels", () => {
  const vib = F.levelScheme(["v"], ["n"], { twoQuantum: true });
  assert.strictEqual(F.enumerate(F.PRESETS.sfg.signs, F.levelScheme(["v"], ["n"])).length, 0); // needs n → g emission
  const sfg = F.enumerate(F.PRESETS.sfg.signs, vib);
  assert.strictEqual(sfg.length, 1);
  assert.deepStrictEqual(sfg[0].rho.map((r) => r.join("")), ["gg", "vg", "ng"]);
  assert.strictEqual(sfg[0].type, "SFG");
  const dfg = F.enumerate(F.PRESETS.dfg.signs, vib);
  assert.ok(dfg.length >= 1 && dfg.every((d) => d.type === "DFG"));
  // third-order results are unchanged by the even-order rule
  assert.strictEqual(F.enumerate(R, threeLevel).length, 3);
});

test("pump–probe preset: rephasing + non-rephasing sets, signal along k_pr", () => {
  const P = F.PRESETS.pumpprobe;
  const ds = P.signSets.flatMap((s) => F.enumerate(s, threeLevel));
  assert.strictEqual(ds.length, 6);
  assert.strictEqual(F.direction(P.signSets[0], P.fieldLabels), "k_pr");
  assert.strictEqual(F.direction(P.signSets[1], P.fieldLabels), "k_pr");
});

test("phase factors: rephasing SE oscillates as e^{+iω_eg t1 − iω_eg t3}", () => {
  const se = F.enumerate(R, twoLevel).find((d) => d.type === "SE");
  const pf = F.phaseFactors(se, twoLevel);
  assert.deepStrictEqual(pf.parts.map((p) => p.sign), [+1, 0, -1]);
  assert.strictEqual(pf.text, "exp(+iω_eg·t₁ −iω_eg·t₃)");
  const nr = F.enumerate(NR, twoLevel).find((d) => d.type === "SE");
  assert.deepStrictEqual(F.phaseFactors(nr, twoLevel).parts.map((p) => p.sign), [-1, 0, -1]);
});

test("layout: rails, n + 1 arrows at 15°, inward arrows for absorption", () => {
  const esa = F.enumerate(R, threeLevel).find((d) => d.type === "ESA");
  const L = F.layout(esa, { intervals: ["t_1", "t_2", "t_3"] });
  assert.strictEqual(L.rails.length, 2);
  assert.strictEqual(L.arrows.length, 4);
  L.arrows.forEach((a) => near(Math.abs(a.y2 - a.y1) / Math.abs(a.x2 - a.x1), Math.tan(Math.PI / 12), 0.03));
  esa.interactions.forEach((s, i) => assert.strictEqual(L.arrows[i].inward, s.kind === "up"));
  assert.ok(L.arrows[3].signal && !L.arrows[3].inward);
  assert.strictEqual(L.brackets.length, 3);
});

test("SVG and TikZ export are well formed", () => {
  const ds = F.enumerate(R, threeLevel);
  const items = ds.map((d) => ({ lay: F.layout(d, { header: d.type }), type: d.type, sign: d.sign }));
  const svg = F.toSVG(items[0].lay, { pal: F.palette("color"), type: items[0].type });
  assert.ok(svg.startsWith("<svg") && svg.endsWith("</svg>"));
  assert.strictEqual((svg.match(/<line /g) || []).length, 2 + 4 + 1); // rails + arrows + time axis
  const wavy = F.toSVG(items[0].lay, { pal: F.palette("mono"), type: items[0].type, wavy: true });
  assert.strictEqual((wavy.match(/<path d="M[^"]*( L[^"]*){5,}"/g) || []).length, 3); // three wavy field arrows (signal stays dashed)
  assert.ok(!/#2B6CB0|#2C7A7B|#C53030/i.test(wavy)); // print mode is black and white
  const tex = F.toTikZ(items, { wavy: true });
  assert.ok(tex.includes("\\usetikzlibrary{arrows.meta, decorations.pathmorphing}"));
  assert.strictEqual((tex.match(/\\begin\{tikzpicture\}/g) || []).length, 3);
  assert.strictEqual((tex.match(/\{/g) || []).length, (tex.match(/\}/g) || []).length);
  assert.ok(tex.includes("\\definecolor{pathESA}{HTML}{C53030}"));
  assert.ok(tex.includes("$|f\\rangle\\langle e|$"));
  assert.ok(!/[^\x00-\x7F]/.test(tex.split("\n").slice(2).join("\n")), "TikZ body must be pure ASCII for pdflatex");
});

test("PNG pHYs chunk for 300 dpi", () => {
  const zlib = require("zlib");
  const sig = [137, 80, 78, 71, 13, 10, 26, 10];
  const ihdr = new Uint8Array(25); ihdr.set([0, 0, 0, 13, 73, 72, 68, 82]);
  const png = new Uint8Array([...sig, ...ihdr, 0, 0, 0, 0, 73, 69, 78, 68, 174, 66, 96, 130]);
  const out = F.pngWithDpi(png, 300);
  assert.strictEqual(out.length, png.length + 21);
  const dv = new DataView(out.buffer, 33);
  assert.strictEqual(String.fromCharCode(...out.subarray(37, 41)), "pHYs");
  assert.strictEqual(dv.getUint32(8), 11811);
  assert.strictEqual(dv.getUint32(17), zlib.crc32(Buffer.from(out.subarray(37, 50))));
});
