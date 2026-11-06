---
layout: default
title: Home
---

## About

I'm Indranil Roy, a PhD candidate in the Solid State and Structural Chemistry Unit (SSCU) at the Indian Institute of Science, Bengaluru, working with Prof. Vivek Tiwari in the [Multidimensional Spectroscopy and Imaging Group](https://sites.google.com/view/musig/). My research uses ultrafast two-dimensional electronic spectroscopy and vibronic models to study energy flow in molecular aggregates and light–matter coupled systems.

My academic homepage, with research and publications, is at **[sites.google.com/iisc.ac.in/indranilroy](https://sites.google.com/iisc.ac.in/indranilroy)**. You can reach me at [indranilroy@iisc.ac.in](mailto:indranilroy@iisc.ac.in).

---

## About these notes

Short, practical write-ups of scripts and command-line workflows I use for data processing and scientific computing: things that took me a while to figure out, written down so others don't have to. They're working notes rather than polished software, so test before relying on them.

---

## Tools

- [Spectroscopy & optics calculators](./tools/): unit conversion (nm, cm⁻¹, eV, THz, Hartree, kJ/mol, …), bandwidth and transform-limited pulses, delay stages, pulse energy and intensity, GDD broadening, material dispersion (n, GVD, TOD) and multi-component setups, prism and grating compressors, SHG/SFG phase matching, Gaussian beams and focusing, autocorrelation deconvolution, sampling/FFT axis design, and spectral Jacobian (nm ⇄ cm⁻¹) conversion of measured spectra.
- [Plot digitizer](./tools/digitizer/): extract numbers from an image of a published plot (linear or log axes, multiple series, CSV export). The image stays in your browser.
- [Double-sided Feynman diagram generator](./tools/dsfd/): enumerate and draw all Liouville pathways for 2D rephasing and non-rephasing, pump–probe, double-quantum and second-order (SFG/DFG) signals in two-level, three-level and excitonic systems, with signs, GSB/SE/ESA classification, wave-vector sums, phase factors and response-function terms; export SVG, 300 dpi PNG or compilable TikZ/LaTeX, in colour or black and white.
- [Multilayer TMM / polariton solver](./tools/tmm/): transfer-matrix R/T/A of planar stacks (metals, Bragg mirrors, Lorentz/Drude/Sellmeier layers, repeat blocks), spectra against energy, wavelength or angle (surface-plasmon resonance), dispersion maps against angle or k<sub>∥</sub>, thickness (detuning) maps, intracavity field profiles, live oscillator sliders, and coupled-oscillator fits of the anticrossing (Rabi splitting, Hopfield coefficients, strong-coupling check); CSV and PNG export.
- [Constants & lab cheat sheet](./reference/): CODATA constants, laser and calibration lines, and the optics rules of thumb that are easy to forget.

Everything runs in your browser.

---

## Recent posts

- [Batch Plotting in Python](./posts/batch-plotting) — overlay many data files in one publication-style plot from the command line.

See all [posts](./posts/).
