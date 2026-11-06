// Tests for spectra-tools.js and the features added for parity with the
// MuSIG unit converter — run with: node --test tests/*.test.js
const test = require("node:test");
const assert = require("node:assert");
const P = require("../assets/js/spectra-tools.js");
const T = require("../assets/js/sci-tools.js");
const O = require("../assets/js/optics-tools.js");

const near = (got, exp, tol) => assert.ok(Math.abs(got - exp) <= tol, `${got} != ${exp} ± ${tol}`);

test("fs⁻¹ (PHz) unit", () => {
  near(T.convertEnergy(800, "nm").PHz, 0.3747405725, 1e-9);
});

test("bandwidth from edges in any unit", () => {
  const a = T.bandwidthFromEdges(795, 805, "nm");
  near(a.dk_cm1, 156.2561038, 1e-6);
  near(a.dnu_PHz, 0.0046844401, 1e-10);
  near(a.domega_radfs, 0.0294332055, 1e-9);
  near(a.tl_gaussian_fs, T.transformLimitedDuration(800, 10, "gaussian"), 1e-9);
  const b = T.bandwidthFromEdges(12000, 12500, "cm1");
  near(b.dl_nm, 33.3333333, 1e-6);
  near(b.dnu_PHz, 0.0149896229, 1e-10);
  near(T.bandwidthFromEdges(0.36, 0.38, "PHz").dnu_THz, 20, 1e-9);
});

test("new materials match catalogue refractive indices", () => {
  near(O.dispersion("sapphire_e", 587.6).n, 1.7600, 2e-4);
  near(O.dispersion("yag", 1064).n, 1.8150, 1e-3);
  near(O.dispersion("calcite_o", 589.3).n, 1.6584, 2e-4);
  near(O.dispersion("znse", 10600).n, 2.403, 1e-3);
});

test("GVD/TOD agree with exact symbolic derivatives", () => {
  // reference values from symbolic differentiation of the Sellmeier equations
  near(O.dispersion("fused_silica", 800).gvd_fs2_mm, 36.16200, 1e-3);
  near(O.dispersion("fused_silica", 800).tod_fs3_mm, 27.49726, 1e-3);
  near(O.dispersion("yag", 800).gvd_fs2_mm, 98.45419, 1e-3);
  near(O.dispersion("bk7", 1030).tod_fs3_mm, 45.52228, 1e-3);
});

test("dispersion stack adds components", () => {
  const s = O.stackDispersion([{ material: "fused_silica", thickness_mm: 2 }, { material: "bk7", thickness_mm: 3 }], 800);
  near(s.gdd_fs2, 2 * 36.162 + 3 * 44.6518, 2e-3);
  assert.strictEqual(s.parts.length, 2);
});

test("two-column parser skips headers and accepts mixed delimiters", () => {
  const p = P.parseTwoColumn("wavelength,intensity\n# comment\n700, 0.1\n701\t0.2\n702;0.3\n703 0.4 extra\n");
  assert.deepStrictEqual(p, [[700, 0.1], [701, 0.2], [702, 0.3], [703, 0.4]]);
});

test("Jacobian conversion conserves area and round-trips", () => {
  const pts = [];
  for (let l = 700; l <= 900; l += 0.5) pts.push([l, Math.exp(-Math.pow((l - 800) / 30, 2))]);
  const k = P.jacobian(pts, "nm2cm");
  near(P.trapz(k), P.trapz(pts), 1e-3);
  assert.ok(k[0][0] < k[k.length - 1][0]); // ascending wavenumber
  const back = P.jacobian(k, "cm2nm");
  back.forEach((p, i) => near(p[1], pts[i][1], 1e-12));
  assert.strictEqual(Math.max(...P.jacobian(pts, "nm2cm", true).map((p) => p[1])), 1);
});

test("thin-film materials for the TMM solver", () => {
  near(O.dispersion("pmma", 587.6).n, 1.4905, 5e-4);
  near(O.dispersion("tio2", 589.3).n, 2.6124, 2e-3);
  near(O.dispersion("ta2o5", 550).n, 2.142, 2e-3);
});
