# Indra's Webpage

Source code of [indras_webpage](https://indranil-sscu.github.io/indras_webpage/), a small Jekyll site on **GitHub Pages** for programming notes, project logs and scientific-computing snippets from my work in spectroscopy, optics and computational chemistry.

> My main academic homepage is **[sites.google.com/iisc.ac.in/indranilroy](https://sites.google.com/iisc.ac.in/indranilroy)**.

## About Me

I am a PhD candidate in the Solid State and Structural Chemistry Unit (SSCU), Indian Institute of Science, Bengaluru, working with Prof. Vivek Tiwari in the [Multidimensional Spectroscopy and Imaging Group](https://sites.google.com/view/musig/). My research uses femtosecond two-dimensional electronic spectroscopy and vibronic model Hamiltonians to study energy flow in molecular aggregates, light-harvesting systems and molecular polaritons.

[Google Scholar](https://scholar.google.com/citations?user=Cz84z-0AAAAJ) · [ORCID](https://orcid.org/0009-0006-3503-9150) · [LinkedIn](https://www.linkedin.com/in/skew-hermitian/) · [indranilroy@iisc.ac.in](mailto:indranilroy@iisc.ac.in)

## What is here

**Indra's Lab** is an interactive science site: laboratory tools for spectroscopy and optics, data utilities, hands-on mathematics experiments for everyone, and tutorials. Everything runs client-side (plain HTML, CSS and JavaScript on GitHub Pages; no frameworks, no tracking).

- **Explore Tools** (`/explore/`): every tool as a card with a live preview, filterable by category and searchable. Cards come from `_data/tools.yml`.
- **Physics & optics:** calculators (`/tools/`), multilayer TMM and polariton solver (`/tools/tmm/`), double-sided Feynman diagram generator (`/tools/dsfd/`), constants cheat sheet (`/reference/`).
- **Data & analysis:** plot digitizer (`/tools/digitizer/`), spectral Jacobian and FFT axis designer (sections of `/tools/`).
- **Science for the Curious** (`/curious/`): probability playground, Monty Hall problem, Fourier series visualiser.
- **Tutorials & Notes** (`/posts/`) and **About** (`/about/`).

Run the tests (scientific calculations against analytical results and benchmarks) with `node --test tests/*.test.js`.

## Structure

```plaintext
.
├── index.html            # Homepage
├── explore.html          # Tool directory (/explore/)
├── curious/              # Science for the Curious: index, probability, monty-hall, fourier
├── tools.html            # Spectroscopy & optics calculators (/tools/)
├── tmm.html              # Multilayer TMM / polariton solver (/tools/tmm/)
├── dsfd.html             # Double-sided Feynman diagram generator (/tools/dsfd/)
├── digitizer.html        # Plot digitizer (/tools/digitizer/)
├── reference.md          # Constants & lab cheat sheet (/reference/)
├── posts.md, posts/      # Tutorials & notes
├── about.md              # About the author and the site
├── _data/                # tools.yml (tool catalogue), categories.yml
├── _includes/            # head (theme script), header, footer, tool card, breadcrumbs, tool metadata
├── _layouts/             # default (skip link, landmarks) and page (breadcrumbs, title, version bar)
├── assets/css/site.css   # Design tokens, light/dark themes and components
├── assets/js/            # Tool engines (sci-, optics-, spectra-tools, tmm-engine, dsfd-core, digitizer-core,
│                         #   curious-core), plotkit (charts), previews (card previews), theme
├── assets/img/           # Social-sharing image
└── tests/                # Node unit tests
```
