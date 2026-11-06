# Indra's Webpage

Source code of [indras_webpage](https://indranil-sscu.github.io/indras_webpage/), a small Jekyll site on **GitHub Pages** for programming notes, project logs and scientific-computing snippets from my work in spectroscopy, optics and computational chemistry.

> My main academic homepage is **[sites.google.com/iisc.ac.in/indranilroy](https://sites.google.com/iisc.ac.in/indranilroy)**.

## About Me

I am a PhD candidate in the Solid State and Structural Chemistry Unit (SSCU), Indian Institute of Science, Bengaluru, working with Prof. Vivek Tiwari in the [Multidimensional Spectroscopy and Imaging Group](https://sites.google.com/view/musig/). My research uses femtosecond two-dimensional electronic spectroscopy and vibronic model Hamiltonians to study energy flow in molecular aggregates, light-harvesting systems and molecular polaritons.

[Google Scholar](https://scholar.google.com/citations?user=Cz84z-0AAAAJ) · [ORCID](https://orcid.org/0009-0006-3503-9150) · [LinkedIn](https://www.linkedin.com/in/skew-hermitian/) · [indranilroy@iisc.ac.in](mailto:indranilroy@iisc.ac.in)

## Features

- Light and dark themes (toggle in the header; follows the system setting until you choose), defined as CSS custom properties in `assets/css/site.css`, with the theme applied before first paint by an inline script in `_includes/head.html`; canvas helpers in `assets/js/theme.js`
- Blog-style code notes and project logs
- Tips and scripts for scientific computing
- [Spectroscopy & optics calculators](https://indranil-sscu.github.io/indras_webpage/tools/) (`tools.html`, logic in `assets/js/sci-tools.js` and `assets/js/optics-tools.js`; plot digitizer in `digitizer.html` + `assets/js/digitizer-core.js`; Feynman diagram generator in `dsfd.html` + `assets/js/dsfd-core.js`; multilayer TMM / polariton solver in `tmm.html` + `assets/js/tmm-engine.js`; reference sheet in `reference.md`; run the tests with `node --test tests/*.test.js`)

## Structure

```plaintext
.
├── index.md       # Home page (Markdown)
├── tools.html     # Calculators page (each with a theory card: equations, symbols, references)
├── tmm.html       # Multilayer TMM / polariton solver (/tools/tmm/)
├── dsfd.html      # Double-sided Feynman diagram generator (/tools/dsfd/)
├── digitizer.html # Plot digitizer (/tools/digitizer/)
├── reference.md   # Constants & lab cheat sheet
├── assets/css/    # Design tokens and light/dark themes (site.css)
├── assets/js/     # Tool logic, no dependencies (sci-, optics-, spectra-tools, tmm-engine, dsfd-core, digitizer-core, theme)
├── tests/         # Node unit tests (node --test tests/*.test.js)
├── posts.md       # Listing of all posts
├── posts/         # Individual project notes
├── _config.yml    # Jekyll site configuration
├── _layouts/      # Page layouts
├── _includes/     # head (theme script), header (theme toggle) and footer partials
└── README.md      # This file
```
