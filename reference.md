---
layout: page
title: Constants & Lab Cheat Sheet
permalink: /reference/
crumbs:
  - title: Explore Tools
    url: /explore/
short_title: Cheat sheet
description: Physical constants (CODATA), laser and calibration lines, and optics rules of thumb for the spectroscopy lab.
version: "1.0"
updated: 2026-11-07
source: reference.md
---

A quick reference for the optics and spectroscopy lab. For conversions, use the [calculators](../tools/).

## Physical constants

Values marked *exact* are fixed by the 2019 SI definitions; the rest are CODATA 2018 recommended values.

| Quantity | Symbol | Value | Unit |
|---|---|---|---|
| Speed of light in vacuum | *c* | 299 792 458 (*exact*) | m s⁻¹ |
| Planck constant | *h* | 6.626 070 15 × 10⁻³⁴ (*exact*) | J s |
| Reduced Planck constant | *ħ* | 1.054 571 817 × 10⁻³⁴ | J s |
| Elementary charge | *e* | 1.602 176 634 × 10⁻¹⁹ (*exact*) | C |
| Boltzmann constant | *k*<sub>B</sub> | 1.380 649 × 10⁻²³ (*exact*) | J K⁻¹ |
| Avogadro constant | *N*<sub>A</sub> | 6.022 140 76 × 10²³ (*exact*) | mol⁻¹ |
| Molar gas constant | *R* | 8.314 462 618 | J mol⁻¹ K⁻¹ |
| Faraday constant | *F* | 96 485.332 12 | C mol⁻¹ |
| Vacuum permittivity | *ε*₀ | 8.854 187 8128 × 10⁻¹² | F m⁻¹ |
| Vacuum permeability | *μ*₀ | 1.256 637 062 12 × 10⁻⁶ | N A⁻² |
| Electron mass | *m*<sub>e</sub> | 9.109 383 7015 × 10⁻³¹ | kg |
| Proton mass | *m*<sub>p</sub> | 1.672 621 923 69 × 10⁻²⁷ | kg |
| Atomic mass constant | *m*<sub>u</sub> | 1.660 539 066 60 × 10⁻²⁷ | kg |
| Fine-structure constant | *α* | 7.297 352 5693 × 10⁻³ (1/*α* = 137.035 999 084) | – |
| Bohr radius | *a*₀ | 5.291 772 109 03 × 10⁻¹¹ | m |
| Hartree energy | *E*<sub>h</sub> | 4.359 744 722 2071 × 10⁻¹⁸ (= 27.211 386 246 eV) | J |
| Rydberg constant | *R*<sub>∞</sub> | 10 973 731.568 160 | m⁻¹ |
| Bohr magneton | *μ*<sub>B</sub> | 9.274 010 0783 × 10⁻²⁴ | J T⁻¹ |
| Stefan–Boltzmann constant | *σ* | 5.670 374 419 × 10⁻⁸ | W m⁻² K⁻⁴ |
| Debye (dipole moment) | D | 3.335 64 × 10⁻³⁰ | C m |

### Handy combinations

| | |
|---|---|
| *hc* | 1239.841 98 eV nm = 1.986 446 × 10⁻²⁵ J m |
| 1 eV | 8065.544 cm⁻¹ = 241.799 THz = 96.485 kJ/mol = 23.061 kcal/mol |
| 1 cm⁻¹ | 29.979 GHz = 0.123 984 meV; one period = 33.356 ps |
| 1 THz | 33.356 cm⁻¹ = 4.1357 meV |
| *k*<sub>B</sub>*T* at 298.15 K | 207.22 cm⁻¹ = 25.693 meV = 2.479 kJ/mol |
| 1 mW at 800 nm | 4.03 × 10¹⁵ photons per second |

## Common laser lines

| Source | Wavelength |
|---|---|
| ArF / KrF / XeCl excimer | 193 / 248 / 308 nm |
| Nitrogen laser | 337.1 nm |
| He–Cd | 325.0, 441.6 nm |
| Ar⁺ ion | 488.0, 514.5 nm |
| He–Ne | 632.8 nm |
| Ruby | 694.3 nm |
| Ti:sapphire | 650–1100 nm (typically ~800 nm) |
| Nd:YAG and harmonics | 1064, 532, 355, 266 nm |
| Yb-doped (fibre, KGW, YAG) | 1030 nm |
| Er-doped fibre | 1550–1560 nm |
| CO₂ | 10.6 µm |
| Common laser diodes | 405, 450, 520, 638, 785, 980 nm |

## Wavelength calibration lines (in air)

| Lamp | Lines (nm) |
|---|---|
| Mercury (Hg) | 253.652, 296.728, 313.155, 365.015, 404.656, 435.833, 546.074, 576.960, 579.066 |
| Neon (Ne) | 585.249, 588.190, 594.483, 614.306, 640.225, 650.653, 659.895, 667.828, 692.947, 703.241 |

## Optics & spectroscopy cheat sheet

**Linewidths and beam sizes**
- Gaussian: FWHM = 2√(2 ln 2) σ = 2.3548 σ.
- Gaussian beam: 1/e² radius *w* = 2σ of the intensity profile; intensity FWHM = 1.1774 *w*.
- Width conversion: Δν̃ [cm⁻¹] ≈ 10⁷ Δλ / λ² (λ, Δλ in nm). Coherence length *l*<sub>c</sub> ≈ λ²/Δλ.

**Polarisation**
- Malus's law: *I* = *I*₀ cos²θ.
- Half-wave plate: rotating the plate by θ rotates linear polarisation by 2θ.
- Quarter-wave plate with its axis at 45° to linear polarisation gives circular polarisation (and vice versa).
- Magic angle 54.7° (cos²θ = 1/3) between pump and probe polarisations removes orientational (anisotropy) dynamics from pump–probe signals.
- Anisotropy *r* = (*S*<sub>∥</sub> − *S*<sub>⊥</sub>)/(*S*<sub>∥</sub> + 2*S*<sub>⊥</sub>); isotropic signal ∝ *S*<sub>∥</sub> + 2*S*<sub>⊥</sub>.

**Interfaces**
- Brewster angle: tan θ<sub>B</sub> = *n*₂/*n*₁ (p-polarised reflection vanishes).
- Normal-incidence reflectance: *R* = ((*n*₁ − *n*₂)/(*n*₁ + *n*₂))², about 4 % per uncoated glass surface.
- Grating equation: *m*λ = *d*(sin θ<sub>i</sub> + sin θ<sub>m</sub>).

**Absorption and populations**
- Beer–Lambert: *A* = −log₁₀ *T* = ε *c* *l*. OD 1 = 10 % transmission, OD 2 = 1 %.
- Boltzmann ratio: *N*₁/*N*₀ = (*g*₁/*g*₀) exp(−Δ*E*/*k*<sub>B</sub>*T*). For a 1000 cm⁻¹ mode at 298 K: ≈ 0.8 %.

**Ultrafast pulses**
- Time–bandwidth products (FWHM): Gaussian 0.441, sech² 0.315.
- Delay: 1 mm of stage travel on a retro-reflector = 6.671 ps; 1 µm = 6.67 fs.
- GVD at 800 nm: fused silica 36 fs²/mm, BK7 45 fs²/mm, sapphire 58 fs²/mm, CaF₂ 28 fs²/mm.
