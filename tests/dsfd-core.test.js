// Unit tests for assets/js/dsfd-core.js — run with: node --test tests/*.test.js
const test = require("node:test");
const assert = require("node:assert");
const F = require("../assets/js/dsfd-core.js");

const twoLevel = F.levelScheme(["e"], []);
const threeLevel = F.levelScheme(["e"], ["f"]);
const dimer = F.levelScheme(["a", "b"], ["f"]);
const R = F.PRESETS.rephasing.signs, NR = F.PRESETS.nonrephasing.signs, DQ = F.PRESETS.dqc.signs;
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
