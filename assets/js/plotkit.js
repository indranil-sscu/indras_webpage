/*
 * plotkit.js: tiny canvas plotting helpers shared by the card previews and the
 * "Science for the Curious" apps. Theme-aware and devicePixelRatio-sharp
 * (via SiteTheme.fit). No dependencies.
 *
 *   const P = PlotKit.frame(canvas, {x0, x1, y0, y1, xlabel, ylabel, bare, xticks, yticks, xfmt, yfmt});
 *   PlotKit.line(P, xs, ys, {color, width, dash});  PlotKit.bars(P, xs, ys, {color, w});
 *   PlotKit.dots(P, xs, ys, {color, r});            PlotKit.legend(P, [{label, color, dash}]);
 */
(function () {
  "use strict";
  const FONT = () => getComputedStyle(document.body).fontFamily;
  const nice = (v) => (Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(1) : String(+v.toPrecision(4)));

  /** Size the canvas to its CSS box (or the given size), draw the frame and return mapping helpers. */
  function frame(canvas, o) {
    o = o || {};
    if (o.height) canvas.style.height = o.height + "px";
    const box = canvas.getBoundingClientRect();
    const w = Math.max(60, Math.round(o.width || box.width || +canvas.getAttribute("width") || 600));
    const h = Math.max(40, Math.round(o.height || box.height || +canvas.getAttribute("height") || 300));
    const { ctx, W, H } = window.SiteTheme.fit(canvas, w, h);
    const K = window.SiteTheme.colors();
    const pad = o.bare ? { l: 6, r: 6, t: 6, b: 6 } : Object.assign({ l: 46, r: 12, t: 12, b: 34 }, o.pad || {});
    ctx.fillStyle = K.bg; ctx.fillRect(0, 0, W, H);
    const x0 = o.x0, x1 = o.x1, y0 = o.y0, y1 = o.y1;
    const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
    const X = (x) => pad.l + ((x - x0) / (x1 - x0 || 1)) * pw;
    const Y = (y) => pad.t + ph - ((y - y0) / (y1 - y0 || 1)) * ph;
    const P = { ctx, W, H, K, X, Y, pad, pw, ph, x0, x1, y0, y1, font: FONT() };
    if (!o.bare) {
      ctx.font = "11px " + P.font; ctx.lineWidth = 1;
      const nx = o.xticks || 5, ny = o.yticks || 4, fx = o.xfmt || nice, fy = o.yfmt || nice;
      ctx.strokeStyle = K.grid;
      for (let i = 0; i <= ny; i++) { const y = pad.t + (ph * i) / ny; ctx.beginPath(); ctx.moveTo(pad.l, Math.round(y) + 0.5); ctx.lineTo(W - pad.r, Math.round(y) + 0.5); ctx.stroke(); }
      ctx.strokeStyle = K.frame; ctx.strokeRect(pad.l + 0.5, pad.t + 0.5, pw, ph);
      ctx.fillStyle = K.muted; ctx.textAlign = "center";
      if (o.xints) {
        // integer-valued axis (counts, harmonics, dice sums): label whole numbers only
        const step = Math.max(1, Math.ceil((x1 - x0) / 14));
        for (let v = Math.ceil(x0); v <= x1; v += step) ctx.fillText(String(v), X(v), H - pad.b + 14);
      } else for (let i = 0; i <= nx; i++) { const v = x0 + ((x1 - x0) * i) / nx; ctx.fillText(fx(v), X(v), H - pad.b + 14); }
      ctx.textAlign = "right";
      for (let i = 0; i <= ny; i++) { const v = y0 + ((y1 - y0) * i) / ny; ctx.fillText(fy(v), pad.l - 5, Y(v) + 4); }
      ctx.fillStyle = K.fg; ctx.textAlign = "center";
      if (o.xlabel) ctx.fillText(o.xlabel, pad.l + pw / 2, H - 4);
      if (o.ylabel) { ctx.save(); ctx.translate(11, pad.t + ph / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(o.ylabel, 0, 0); ctx.restore(); }
    }
    return P;
  }
  function line(P, xs, ys, o) {
    o = o || {};
    const c = P.ctx; c.save(); c.beginPath(); c.rect(P.pad.l, P.pad.t, P.pw, P.ph); c.clip();
    c.strokeStyle = o.color || P.K.accent; c.lineWidth = o.width || 1.6; c.setLineDash(o.dash || []); c.lineJoin = "round";
    c.beginPath(); let on = false;
    xs.forEach((x, i) => { const y = ys[i]; if (!isFinite(y)) { on = false; return; } on ? c.lineTo(P.X(x), P.Y(y)) : c.moveTo(P.X(x), P.Y(y)); on = true; });
    c.stroke(); c.restore();
  }
  function bars(P, xs, ys, o) {
    o = o || {};
    const c = P.ctx, w = (o.w || 0.8) * (P.pw / Math.max(1, (P.x1 - P.x0)));
    c.fillStyle = o.color || P.K.accent;
    xs.forEach((x, i) => { const y0 = P.Y(Math.max(P.y0, 0)), y = P.Y(ys[i]); c.fillRect(P.X(x) - w / 2, Math.min(y, y0), w, Math.abs(y0 - y)); });
  }
  function dots(P, xs, ys, o) {
    o = o || {};
    const c = P.ctx; c.fillStyle = o.color || P.K.fg;
    xs.forEach((x, i) => { c.beginPath(); c.arc(P.X(x), P.Y(ys[i]), o.r || 2.6, 0, 2 * Math.PI); c.fill(); });
  }
  function legend(P, items, x, y) {
    const c = P.ctx; c.font = "11px " + P.font; c.textAlign = "left";
    const w = Math.max(...items.map((it) => c.measureText(it.label).width)) + 32;
    x = x === undefined ? P.W - P.pad.r - w - 4 : x; y = y === undefined ? P.pad.t + 4 : y;
    c.globalAlpha = 0.9; c.fillStyle = P.K.bg; c.fillRect(x, y, w, items.length * 15 + 6); c.globalAlpha = 1;
    items.forEach((it, i) => {
      const yy = y + 12 + i * 15;
      c.strokeStyle = it.color; c.fillStyle = it.color; c.lineWidth = 2; c.setLineDash(it.dash || []);
      if (it.dot) { c.beginPath(); c.arc(x + 13, yy - 3, 2.8, 0, 2 * Math.PI); c.fill(); }
      else { c.beginPath(); c.moveTo(x + 5, yy - 3); c.lineTo(x + 21, yy - 3); c.stroke(); }
      c.setLineDash([]); c.fillStyle = P.K.fg; c.fillText(it.label, x + 26, yy);
    });
  }
  /** Series colours that work in both themes. */
  function palette() {
    const dark = window.SiteTheme.get() === "dark";
    return dark ? { a: "#38bdf8", b: "#fb923c", c: "#34d399", d: "#c4b5fd", e: "#f87171" } : { a: "#0369a1", b: "#c2410c", c: "#047857", d: "#7c3aed", e: "#b91c1c" };
  }
  window.PlotKit = { frame, line, bars, dots, legend, palette, nice };
})();
