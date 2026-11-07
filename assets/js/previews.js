/*
 * previews.js: small live previews on the tool cards, drawn with the tools' own
 * code (TMM engine, Feynman-diagram core, optics data, curious-core) rather than
 * pictures. Each preview is drawn only when its card scrolls into view, and is
 * redrawn when the theme or the card size changes.
 */
(function () {
  "use strict";
  const PK = () => window.PlotKit;
  const svgNS = "http://www.w3.org/2000/svg";
  const canvasIn = (el) => { let c = el.querySelector("canvas"); if (!c) { el.innerHTML = ""; c = document.createElement("canvas"); el.appendChild(c); } return c; };
  const svgIn = (el, markup) => { el.innerHTML = markup; };
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  const VIRIDIS = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]];
  const cmap = (t) => { t = Math.min(1, Math.max(0, t)) * 4; const i = Math.min(3, Math.floor(t)), f = t - i; return VIRIDIS[i].map((c, k) => Math.round(c + f * (VIRIDIS[i + 1][k] - c))); };

  const R = {
    /* Ag / J-aggregate / Ag microcavity: angle-resolved reflectance computed with the TMM engine */
    tmm(el) {
      const T = window.TMM; if (!T) return;
      const M = T.material, ag = M.drude(3.7, 9.1, 0.018), cm = 8065.543937;
      const st = { incident: M.constant(1), substrate: M.constant(1.52), layers: [{ material: ag, d_nm: 30 }, { material: M.lorentz(2.2, [{ nu0: 2.1 * cm, gamma: 0.05 * cm, f: 0.08 }]), d_nm: 155 }, { material: ag, d_nm: 100 }] };
      const nx = 48, ny = 64, th = Array.from({ length: nx }, (_, i) => (60 * i) / (nx - 1)), E = Array.from({ length: ny }, (_, j) => 1.7 + (0.85 * j) / (ny - 1));
      const map = T.angleMap(st, "s", "R", th, E);
      const P = PK().frame(canvasIn(el), { x0: 0, x1: 60, y0: 1.7, y1: 2.55, bare: true });
      const off = document.createElement("canvas"); off.width = nx; off.height = ny;
      const o = off.getContext("2d"), img = o.createImageData(nx, ny);
      const flat = map.Z.flat(), lo = Math.min(...flat), hi = Math.max(...flat);
      for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) { const [r, g, b] = cmap((map.Z[i][j] - lo) / (hi - lo)); const k = ((ny - 1 - j) * nx + i) * 4; img.data.set([r, g, b, 255], k); }
      o.putImageData(img, 0, 0);
      P.ctx.imageSmoothingEnabled = true; P.ctx.drawImage(off, 0, 0, P.W, P.H);
      P.ctx.setLineDash([4, 4]); P.ctx.strokeStyle = "rgba(255,255,255,.75)"; P.ctx.lineWidth = 1;
      const y = P.H - ((2.1 - 1.7) / 0.85) * P.H; P.ctx.beginPath(); P.ctx.moveTo(0, y); P.ctx.lineTo(P.W, y); P.ctx.stroke(); P.ctx.setLineDash([]);
      el.setAttribute("aria-label", "Polariton anticrossing in an angle-resolved reflectance map computed with the TMM solver");
    },
    /* two pathways from the Feynman-diagram generator */
    dsfd(el) {
      const F = window.DSFD; if (!F) return;
      const dark = window.SiteTheme.get() === "dark";
      const pal = dark ? { ink: "#E2E8F0", bg: "transparent", muted: "#94A3B8", accent: () => "#E2E8F0", mono: true } : { ...F.palette("mono"), bg: "transparent" };
      const sc = F.levelScheme(["e"], ["f"]), ds = F.enumerate([-1, 1, 1], sc).sort((a, b) => F.hzLabel(a).localeCompare(F.hzLabel(b)));
      const pick = [ds[0], ds[2]].map((d) => F.toSVG(F.layout(d, { groundZero: true }), { pal, type: d.type }).replace("<svg ", '<svg preserveAspectRatio="xMidYMid meet" '));
      svgIn(el, `<div style="display:flex;justify-content:center;gap:4%;height:100%;padding:6px 0;box-sizing:border-box">${pick.map((s) => `<div style="height:100%;flex:0 1 38%">${s}</div>`).join("")}</div>`);
      el.querySelectorAll("rect").forEach((r) => r.setAttribute("fill", "transparent"));
      el.setAttribute("aria-label", "Two double-sided Feynman diagrams drawn by the generator");
    },
    /* Gaussian pulse in time with its spectrum (transform-limited) */
    pulse(el) {
      const c = 299.792458; // nm/fs
      const l0 = 800, dl = 30, tau = window.SciTools ? window.SciTools.transformLimitedDuration(l0, dl, "gaussian") : (0.4413 * l0 * l0) / (c * dl);
      const w = (2 * Math.PI * c) / l0;
      const P = PK().frame(canvasIn(el), { x0: -2.4 * tau, x1: 2.4 * tau, y0: -1.15, y1: 1.25, bare: true });
      const ts = Array.from({ length: 600 }, (_, i) => -2.4 * tau + (4.8 * tau * i) / 599);
      const env = ts.map((t) => Math.exp(-2 * Math.LN2 * (t / tau) ** 2));
      const pal = PK().palette();
      PK().line(P, ts, ts.map((t, i) => env[i] * Math.cos(w * t)), { color: P.K.muted, width: 1 });
      PK().line(P, ts, env, { color: pal.a, width: 2 });
      PK().line(P, ts, env.map((v) => -v), { color: pal.a, width: 2, dash: [3, 3] });
      P.ctx.fillStyle = P.K.fg; P.ctx.font = "11px " + P.font; P.ctx.textAlign = "left";
      P.ctx.fillText(`Δλ = ${dl} nm at ${l0} nm  →  Δt = ${tau.toFixed(1)} fs`, 10, P.H - 10);
      el.setAttribute("aria-label", `Transform-limited Gaussian pulse of ${tau.toFixed(1)} fs`);
    },
    /* fused-silica index and GVD from the optics module */
    dispersion(el) {
      const O = window.OpticsTools; if (!O) return;
      const ls = Array.from({ length: 160 }, (_, i) => 400 + (1200 * i) / 159);
      const d = ls.map((l) => O.dispersion("fused_silica", l));
      const n = d.map((x) => x.n), g = d.map((x) => x.gvd_fs2_mm);
      const pal = PK().palette();
      const P = PK().frame(canvasIn(el), { x0: 400, x1: 1600, y0: 0, y1: 1, bare: true, pad: {} });
      const norm = (a) => { const lo = Math.min(...a), hi = Math.max(...a); return a.map((v) => 0.12 + (0.76 * (v - lo)) / (hi - lo)); };
      PK().line(P, ls, norm(n), { color: pal.a, width: 2 }); PK().line(P, ls, norm(g), { color: pal.b, width: 2 });
      PK().legend(P, [{ label: "n(λ)", color: pal.a }, { label: "GVD(λ)", color: pal.b }]);
      P.ctx.fillStyle = P.K.muted; P.ctx.font = "11px " + P.font; P.ctx.fillText("fused silica, 400–1600 nm", 10, P.H - 10);
      el.setAttribute("aria-label", "Refractive index and group-velocity dispersion of fused silica");
    },
    /* a schematic of the digitizer workflow: calibrate four axis points, then click the data */
    digitizer(el) {
      const ink = css("--text-secondary"), acc = css("--accent"), mut = css("--border-strong");
      const pts = Array.from({ length: 9 }, (_, i) => { const x = 40 + i * 27; return [x, 120 - 85 * Math.exp(-i / 2.6)]; });
      svgIn(el, `<svg xmlns="${svgNS}" viewBox="0 0 300 150" preserveAspectRatio="xMidYMid meet"><g fill="none" stroke="${mut}"><path d="M36 128V18M36 128H286"/></g>
        <path d="M${pts.map((p) => p.join(",")).join(" L")}" fill="none" stroke="${ink}" stroke-width="1.5"/>
        ${pts.filter((_, i) => i % 2).map((p) => `<circle cx="${p[0]}" cy="${p[1]}" r="3.2" fill="${acc}"/>`).join("")}
        ${[[60, 128, "X₁"], [250, 128, "X₂"], [36, 110, "Y₁"], [36, 30, "Y₂"]].map(([x, y, t]) => `<g stroke="${acc}" stroke-width="1.2"><path d="M${x - 6} ${y}h12M${x} ${y - 6}v12"/></g><text x="${x + 6}" y="${y - 6}" font-size="10" fill="${acc}" style="font-family:var(--font-mono)">${t}</text>`).join("")}
      </svg>`);
      el.setAttribute("aria-label", "Digitizing a plot: four calibration points and picked data points");
    },
    /* wavelength → wavenumber conversion with and without the Jacobian */
    jacobian(el) {
      const S = window.SpectraTools; if (!S) return;
      const pts = []; for (let l = 450; l <= 1100; l += 2) pts.push([l, Math.exp(-(((l - 600) / 50) ** 2)) + 0.8 * Math.exp(-(((l - 900) / 70) ** 2))]);
      const jac = S.jacobian(pts, "nm2cm", true);
      const naive = pts.map(([x, y]) => [1e7 / x, y]).sort((a, b) => a[0] - b[0]);
      const pal = PK().palette();
      const P = PK().frame(canvasIn(el), { x0: 9000, x1: 22300, y0: -0.05, y1: 1.15, bare: true });
      PK().line(P, naive.map((p) => p[0]), naive.map((p) => p[1]), { color: P.K.muted, width: 1.5, dash: [4, 3] });
      PK().line(P, jac.map((p) => p[0]), jac.map((p) => p[1]), { color: pal.b, width: 2 });
      PK().legend(P, [{ label: "with Jacobian", color: pal.b }, { label: "axis relabelled", color: P.K.muted, dash: [4, 3] }]);
      el.setAttribute("aria-label", "A spectrum converted to wavenumber with and without the Jacobian factor");
    },
    /* aliasing: a fast signal sampled too slowly looks slow */
    sampling(el) {
      const pal = PK().palette();
      const P = PK().frame(canvasIn(el), { x0: 0, x1: 1, y0: -1.3, y1: 1.3, bare: true });
      const ts = Array.from({ length: 500 }, (_, i) => i / 499), f = 9, fs = 10;
      PK().line(P, ts, ts.map((t) => Math.sin(2 * Math.PI * f * t)), { color: P.K.muted, width: 1 });
      PK().line(P, ts, ts.map((t) => -Math.sin(2 * Math.PI * (fs - f) * t)), { color: pal.b, width: 2, dash: [5, 3] });
      const tsS = Array.from({ length: fs + 1 }, (_, i) => i / fs);
      PK().dots(P, tsS, tsS.map((t) => Math.sin(2 * Math.PI * f * t)), { color: pal.a, r: 3.4 });
      el.setAttribute("aria-label", "Aliasing: samples of a fast oscillation trace out a slow one");
    },
    /* 2000 rolls of two dice against the exact distribution */
    dice(el) {
      const C = window.Curious; if (!C) return;
      const r = C.rng(2024), N = 2000, counts = new Array(11).fill(0);
      for (let i = 0; i < N; i++) { const f = C.rollDice(2, 6, r); counts[f[0] + f[1] - 2]++; }
      const th = C.diceSumDistribution(2, 6).p, xs = counts.map((_, i) => i + 2), pal = PK().palette();
      const P = PK().frame(canvasIn(el), { x0: 1.3, x1: 12.7, y0: 0, y1: 0.2, bare: true });
      PK().bars(P, xs, counts.map((c) => c / N), { color: pal.a, w: 0.72 });
      PK().dots(P, xs, th, { color: pal.b, r: 3.2 });
      PK().legend(P, [{ label: "2000 rolls", color: pal.a }, { label: "theory", color: pal.b, dot: true }], 8, 8);
      el.setAttribute("aria-label", "Histogram of 2000 simulated rolls of two dice compared with the exact distribution");
    },
    /* Monty Hall: three doors and the long-run win rates of staying and switching */
    montyhall(el) {
      const C = window.Curious; if (!C) return;
      const r = C.rng(3), N = 3000; let s = 0, w = 0;
      for (let i = 0; i < N; i++) { const g = C.montyGame(3, r); s += g.winStay; w += g.winSwitch; }
      const acc = css("--accent"), b = css("--accent-accent"), ink = css("--text-primary"), mut = css("--border-strong"), surf = css("--bg-surface");
      const door = (x, open, label) => `<g><rect x="${x}" y="22" width="44" height="70" rx="3" fill="${open ? "none" : surf}" stroke="${mut}" stroke-width="1.5"/>${open ? `<text x="${x + 22}" y="64" text-anchor="middle" font-size="11" fill="${ink}" style="font-family:var(--font-sans)">goat</text>` : `<circle cx="${x + 36}" cy="58" r="2.2" fill="${mut}"/>`}<text x="${x + 22}" y="106" text-anchor="middle" font-size="10" fill="${ink}" style="font-family:var(--font-mono)">${label}</text></g>`;
      const bar = (y, f, col, lab) => `<text x="184" y="${y + 9}" font-size="10" fill="${ink}" style="font-family:var(--font-sans)">${lab}</text><rect x="228" y="${y}" width="${60 * f}" height="11" rx="2" fill="${col}"/><text x="${232 + 60 * f}" y="${y + 9}" font-size="10" fill="${ink}" style="font-family:var(--font-mono)">${(100 * f).toFixed(0)}%</text>`;
      svgIn(el, `<svg xmlns="${svgNS}" viewBox="0 0 300 130" preserveAspectRatio="xMidYMid meet">${door(14, false, "pick")}${door(66, true, "")}${door(118, false, "?")}
        ${bar(42, s / N, b, "stay")}${bar(66, w / N, acc, "switch")}<text x="184" y="100" font-size="9.5" fill="${ink}" opacity=".7" style="font-family:var(--font-sans)">${N} simulated games</text></svg>`);
      el.setAttribute("aria-label", `Monty Hall: in ${N} simulated games staying won ${(100 * s / N).toFixed(0)} percent, switching ${(100 * w / N).toFixed(0)} percent`);
    },
    /* square wave and its Fourier partial sums */
    fourier(el) {
      const C = window.Curious; if (!C) return;
      const pal = PK().palette(), xs = Array.from({ length: 500 }, (_, i) => -Math.PI + (2 * Math.PI * i) / 499);
      const P = PK().frame(canvasIn(el), { x0: -Math.PI, x1: Math.PI, y0: -1.45, y1: 1.45, bare: true });
      PK().line(P, xs, xs.map((x) => C.WAVES.square.f(x)), { color: P.K.muted, width: 1.2, dash: [4, 3] });
      [[1, pal.c], [3, pal.b], [15, pal.a]].forEach(([K, col]) => { const co = C.fourierCoefficients("square", K); PK().line(P, xs, xs.map((x) => C.partialSum(co, x)), { color: col, width: 1.7 }); });
      PK().legend(P, [{ label: "1 harmonic", color: pal.c }, { label: "3", color: pal.b }, { label: "15", color: pal.a }]);
      el.setAttribute("aria-label", "A square wave built from 1, 3 and 15 Fourier harmonics");
    },
    /* birthday paradox: exact curve with simulated rooms of each size */
    birthday(el) {
      const C = window.Curious; if (!C) return;
      const r = C.rng(23), pal = PK().palette(), xs = Array.from({ length: 70 }, (_, i) => i + 1);
      const P = PK().frame(canvasIn(el), { x0: 0, x1: 70, y0: -0.04, y1: 1.06, bare: true });
      PK().line(P, [0, 70], [0.5, 0.5], { color: P.K.muted, width: 1, dash: [3, 4] });
      PK().line(P, [23, 23], [-0.04, 1.06], { color: P.K.muted, width: 1, dash: [3, 4] });
      PK().line(P, xs, xs.map((n) => C.birthdayExact(n)), { color: pal.a, width: 2 });
      const sx = xs.filter((n) => n % 3 === 2), seen = new Uint8Array(365);
      PK().dots(P, sx, sx.map((n) => { let h = 0; for (let i = 0; i < 300; i++) h += C.birthdayHasMatch(n, 365, r, seen); return h / 300; }), { color: pal.b, r: 2.6 });
      P.ctx.fillStyle = P.K.fg; P.ctx.font = "11px " + P.font; P.ctx.textAlign = "left";
      P.ctx.fillText("23 people → 50.7%", P.X(23) + 6, P.Y(0.5) + 15);
      el.setAttribute("aria-label", "Probability of a shared birthday against group size: exact curve and simulated rooms, crossing 50% at 23 people");
    },
    /* double pendulum: 6 s of the lower bob's path, integrated with RK4 */
    pendulum(el) {
      const Ph = window.Physics; if (!Ph) return;
      const p = { m1: 1, m2: 1, L1: 1, L2: 1, g: 9.81 }, pal = PK().palette();
      const P = PK().frame(canvasIn(el), { x0: -1, x1: 1, y0: -1, y1: 1, bare: true }), c = P.ctx;
      const sc = (P.H / 2 - 10) / 2, ox = P.W / 2, oy = P.H / 2, px = (x, y) => [ox + x * sc, oy - y * sc];
      let y = [2.2, 2.6, 0, 0], z = [2.2 + 1e-4, 2.6, 0, 0];
      const pa = [], pb = [];
      for (let i = 0; i < 3000; i++) { y = Ph.doubleStep(y, p, 0.002); z = Ph.doubleStep(z, p, 0.002); if (i % 4 === 0) { const q = Ph.doublePositions(y, p), u = Ph.doublePositions(z, p); pa.push(px(q.x2, q.y2)); pb.push(px(u.x2, u.y2)); } }
      [[pb, pal.b], [pa, pal.a]].forEach(([pts, col]) => { c.strokeStyle = col; c.lineWidth = 1.3; c.globalAlpha = 0.8; c.beginPath(); pts.forEach((q, i) => (i ? c.lineTo(...q) : c.moveTo(...q))); c.stroke(); });
      c.globalAlpha = 1;
      const q = Ph.doublePositions(y, p), b1 = px(q.x1, q.y1), b2 = px(q.x2, q.y2);
      c.strokeStyle = P.K.fg; c.lineWidth = 1.6; c.beginPath(); c.moveTo(ox, oy); c.lineTo(...b1); c.lineTo(...b2); c.stroke();
      c.fillStyle = pal.a; [b1, b2].forEach((b) => { c.beginPath(); c.arc(b[0], b[1], 4.5, 0, 2 * Math.PI); c.fill(); });
      c.fillStyle = P.K.fg; c.beginPath(); c.arc(ox, oy, 2.5, 0, 2 * Math.PI); c.fill();
      el.setAttribute("aria-label", "Paths of two double pendulums started almost identically, diverging chaotically");
    },
    /* the Mandelbrot set, coloured by smooth escape time */
    fractal(el) {
      const Fr = window.Fractal; if (!Fr) return;
      const cv = canvasIn(el), P = PK().frame(cv, { x0: 0, x1: 1, y0: 0, y1: 1, bare: true });
      const dpr = Math.min(2, window.devicePixelRatio || 1), w = Math.round(P.W * dpr), h = Math.round(P.H * dpr);
      const buf = new Float32Array(w * h), view = { cx: -0.75, cy: 0, scale: Math.max(2.9 / w, 2.5 / h) };
      Fr.renderRows(buf, w, h, view, { mode: "mandelbrot", maxIter: 120, R: 256, smooth: true }, 0, h);
      const off = document.createElement("canvas"); off.width = w; off.height = h;
      const o = off.getContext("2d"), img = o.createImageData(w, h), ins = window.SiteTheme.get() === "dark" ? [6, 9, 18] : [15, 23, 42];
      for (let k = 0; k < w * h; k++) img.data.set([...Fr.colorOf(buf[k], "ocean", 24, ins), 255], 4 * k);
      o.putImageData(img, 0, 0); P.ctx.drawImage(off, 0, 0, P.W, P.H);
      el.setAttribute("aria-label", "The Mandelbrot set coloured by escape time");
    },
    reference(el) {
      const ink = css("--text-primary"), mut = css("--text-muted"), acc = css("--accent");
      const rows = [["c", "299 792 458 m/s", "exact"], ["h", "6.626 070 15e−34 J s", "exact"], ["1 eV", "8065.544 cm⁻¹", ""], ["800 nm", "12 500 cm⁻¹", "1.55 eV"]];
      svgIn(el, `<svg xmlns="${svgNS}" viewBox="0 0 300 130" preserveAspectRatio="xMidYMid meet">${rows.map(([a, b, c], i) => `<text x="22" y="${32 + i * 26}" font-size="12" fill="${acc}" style="font-family:var(--font-mono)">${a}</text><text x="80" y="${32 + i * 26}" font-size="11.5" fill="${ink}" style="font-family:var(--font-mono)">${b}</text><text x="282" y="${32 + i * 26}" text-anchor="end" font-size="9.5" fill="${mut}" style="font-family:var(--font-sans)">${c}</text>`).join("")}</svg>`);
      el.setAttribute("aria-label", "Physical constants and conversions from the cheat sheet");
    },
    python(el) {
      const ink = css("--text-primary"), mut = css("--text-muted"), kw = css("--accent"), str = css("--accent-accent");
      const L = [[["for", kw], [" f ", ink], ["in", kw], [" sorted(files):", ink]], [["    x, y = np.loadtxt(f, unpack=", ink], ["True", kw], [")", ink]], [["    ax.plot(x, y, label=f.stem)", ink]], [["ax.legend(); fig.savefig(", ink], ["\"overlay.pdf\"", str], [")", ink]], [["# one command, every trial in one figure", mut]]];
      svgIn(el, `<svg xmlns="${svgNS}" viewBox="0 0 300 130" preserveAspectRatio="xMidYMid meet">${L.map((parts, i) => `<text x="12" y="${28 + i * 21}" font-size="10.5" style="font-family:var(--font-mono)" xml:space="preserve">${parts.map(([t, c]) => `<tspan fill="${c}">${t.replace(/&/g, "&amp;").replace(/</g, "&lt;")}</tspan>`).join("")}</text>`).join("")}</svg>`);
      el.setAttribute("aria-label", "Python code that overlays many data files in one figure");
    },
    equation(el) {
      const ink = css("--text-primary"), mut = css("--text-muted"), acc = css("--accent");
      const serif = "Cambria Math, STIX Two Math, Latin Modern Math, Cambria, Times New Roman, serif";
      svgIn(el, `<svg xmlns="${svgNS}" viewBox="0 0 300 130" preserveAspectRatio="xMidYMid meet" font-family="${serif}">
        <text x="150" y="46" text-anchor="middle" font-size="16" fill="${ink}"><tspan font-style="italic">M</tspan><tspan font-size="11" dy="4">j</tspan><tspan dy="-4"> = [ cos </tspan><tspan font-style="italic">δ</tspan><tspan>   −</tspan><tspan font-style="italic">i</tspan><tspan> sin </tspan><tspan font-style="italic">δ</tspan><tspan>/</tspan><tspan font-style="italic">η</tspan><tspan> ; … ]</tspan></text>
        <text x="150" y="82" text-anchor="middle" font-size="16" fill="${acc}">ħΩ<tspan font-size="11" dy="4">R</tspan><tspan dy="-4"> = 2</tspan><tspan font-style="italic">g</tspan><tspan> ∝ √</tspan><tspan font-style="italic">f</tspan></text>
        <text x="150" y="112" text-anchor="middle" font-size="11" fill="${mut}" style="font-family:var(--font-sans)">theory notes with references</text></svg>`);
      el.setAttribute("aria-label", "Transfer-matrix and Rabi-splitting equations");
    },
  };

  const drawn = new Set();
  function draw(el) {
    const fn = R[el.dataset.preview]; if (!fn) return;
    try { fn(el); el.setAttribute("role", "img"); drawn.add(el); } catch (e) { console.error("preview", el.dataset.preview, e); }
  }
  function init() {
    const els = [...document.querySelectorAll(".preview[data-preview]")];
    if (!els.length) return;
    if ("IntersectionObserver" in window) {
      const io = new IntersectionObserver((entries) => entries.forEach((e) => { if (e.isIntersecting) { io.unobserve(e.target); draw(e.target); } }), { rootMargin: "200px" });
      els.forEach((el) => io.observe(el));
    } else els.forEach(draw);
    window.SiteTheme.onChange(() => drawn.forEach(draw));
    let t; window.addEventListener("resize", () => { clearTimeout(t); t = setTimeout(() => drawn.forEach(draw), 200); });
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init); else init();
})();
