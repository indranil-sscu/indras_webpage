// Unit tests for assets/js/optics-tools.js — run with: node --test tests/*.test.js
const test = require("node:test");
const assert = require("node:assert");
const O = require("../assets/js/optics-tools.js");

const near = (got, exp, tol) => assert.ok(Math.abs(got - exp) <= tol, `${got} != ${exp} ± ${tol}`);

test("refractive indices match catalogue values", () => {
  near(O.dispersion("bk7", 587.56).n, 1.5168, 1e-4); // N-BK7 n_d
  near(O.dispersion("fused_silica", 587.56).n, 1.4585, 1e-4); // fused silica n_d
  near(O.dispersion("sapphire_o", 587.56).n, 1.7682, 1e-4); // sapphire n_o,d
  near(O.dispersion("fused_silica", 800).n, 1.4533, 1e-4);
  near(O.dispersion("bbo_o", 800).n, 1.6614, 1e-4);
});

test("group index and GVD at 800 nm", () => {
  near(O.dispersion("fused_silica", 800).ng, 1.4671, 1e-4);
  near(O.dispersion("fused_silica", 800).gvd_fs2_mm, 36.16, 0.05);
  near(O.dispersion("bk7", 800).gvd_fs2_mm, 44.65, 0.05);
  near(O.dispersion("caf2", 800).gvd_fs2_mm, 27.9, 0.1);
  near(O.dispersion("sapphire_o", 800).gvd_fs2_mm, 58.0, 0.2);
  near(O.dispersion("sf10", 800).gvd_fs2_mm, 156.5, 0.3);
  near(O.dispersion("fused_silica", 800).tod_fs3_mm, 27.5, 0.2);
});

test("phase-matching angles", () => {
  near(O.phaseMatch("kdp", 1064, 1064, "I").theta_deg, 41.2, 0.2);
  near(O.phaseMatch("kdp", 1064, 1064, "II").theta_deg, 59.1, 0.3);
  near(O.phaseMatch("bbo", 1064, 1064, "I").theta_deg, 22.8, 0.2);
  near(O.phaseMatch("bbo", 800, 800, "I").theta_deg, 29.1, 0.3);
  near(O.phaseMatch("bbo", 1064, 1064, "I").walkoff_mrad, 55.8, 1);
  assert.strictEqual(O.phaseMatch("bbo", 400, 400, "I").possible, false); // 200 nm SHG not phase-matchable
});

test("compressors", () => {
  const p = O.prismCompressor("fused_silica", 800, 1000, 0);
  near(p.apex_deg, 69.1, 0.1); // Brewster prism apex angle for fused silica
  near(O.prismCompressor("sf10", 800, 100, 0).apex_deg, 60.6, 0.1);
  assert.ok(p.gdd_fs2 < 0);
  const g = O.gratingCompressor(800, 1200, 30, 100);
  assert.ok(g.gdd_fs2 < 0 && g.tod_fs3 > 0);
  // GDD scales linearly with separation
  near(O.gratingCompressor(800, 1200, 30, 200).gdd_fs2 / g.gdd_fs2, 2, 1e-12);
});

test("Gaussian beams", () => {
  near(O.gaussianBeam(800, 1000, 1, 0).rayleigh_mm, 3926.99, 0.01);
  near(O.focusSpot(800, 4, 100, 1).w0_um, 12.732, 0.001);
  near(O.gaussianBeam(800, 10, 1, O.gaussianBeam(800, 10, 1, 0).rayleigh_mm).w_at_z_um, 10 * Math.SQRT2, 1e-9);
});

test("autocorrelation factors", () => {
  near(O.pulseFromAutocorrelation(141.42, "gaussian"), 100, 0.01);
  near(O.pulseFromAutocorrelation(154.27, "sech2"), 100, 0.01);
});

test("FFT axis and aliasing", () => {
  const a = O.fftAxis(1, 1000, 0);
  near(a.resolution_cm1, 33.356, 0.001);
  near(a.nyquist_cm1, 16678.2, 0.1);
  const al = O.aliasedFrequency(12500, 2, 0); // 800 nm sampled every 2 fs aliases
  assert.ok(al.aliased);
  near(Math.abs(al.apparent_cm1), 16678.2 - 12500, 0.1);
  assert.ok(!O.aliasedFrequency(12500, 2, 12000).aliased); // rotating frame removes aliasing
});
