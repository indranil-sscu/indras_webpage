/*
 * theme.js — light/dark theme switching and canvas helpers shared by all pages.
 * The initial theme is applied by the inline script in _includes/head.html;
 * this file wires the toggle button and exposes window.SiteTheme:
 *
 *   SiteTheme.get()            "light" | "dark"
 *   SiteTheme.set(t)           switch and remember the choice
 *   SiteTheme.colors()         canvas colours read from the CSS custom properties
 *   SiteTheme.fit(canvas)      size the backing store for devicePixelRatio;
 *                              returns {ctx, W, H} in CSS pixels
 *   SiteTheme.onChange(fn)     call fn() after every theme change
 *   SiteTheme.savePNG(canvas, name)
 */
(function () {
  "use strict";
  var root = document.documentElement;
  var listeners = [];

  function get() { return root.getAttribute("data-theme") === "dark" ? "dark" : "light"; }
  function notify() {
    listeners.forEach(function (fn) { try { fn(get()); } catch (e) { console.error(e); } });
    window.dispatchEvent(new CustomEvent("themechange", { detail: { theme: get() } }));
  }
  function set(t, remember) {
    root.setAttribute("data-theme", t === "dark" ? "dark" : "light");
    if (remember !== false) { try { localStorage.setItem("theme", get()); } catch (e) {} }
    notify();
  }

  // follow the system setting until the visitor chooses explicitly
  if (window.matchMedia) {
    var mq = window.matchMedia("(prefers-color-scheme: dark)");
    var follow = function (e) {
      var saved = null;
      try { saved = localStorage.getItem("theme"); } catch (err) {}
      if (saved !== "light" && saved !== "dark") set(e.matches ? "dark" : "light", false);
    };
    if (mq.addEventListener) mq.addEventListener("change", follow); else if (mq.addListener) mq.addListener(follow);
  }

  function css(name, fallback) {
    var v = getComputedStyle(root).getPropertyValue(name).trim();
    return v || fallback;
  }
  function colors() {
    return {
      bg: css("--canvas-bg", "#fff"),
      fg: css("--canvas-fg", "#222"),
      muted: css("--canvas-muted", "#888"),
      grid: css("--canvas-grid", "#e6e8ec"),
      frame: css("--canvas-frame", "#ccc"),
      accent: css("--accent", "#1f63c6"),
      accent2: css("--accent-accent", "#c2410c"),
      surface: css("--bg-surface", "#fff"),
      text: css("--text-primary", "#111"),
      dark: get() === "dark",
    };
  }

  /** Scale a canvas for crisp rendering on high-DPI screens. The logical size is
   *  the canvas's original width/height attributes (kept in data-w / data-h). */
  function fit(canvas, w, h) {
    if (!canvas.dataset.w) { canvas.dataset.w = canvas.getAttribute("width") || canvas.width; canvas.dataset.h = canvas.getAttribute("height") || canvas.height; }
    var W = w || +canvas.dataset.w, H = h || +canvas.dataset.h;
    if (w) canvas.dataset.w = w;
    if (h) canvas.dataset.h = h;
    var dpr = Math.max(1, Math.min(3, window.devicePixelRatio || 1));
    canvas.width = Math.round(W * dpr);
    canvas.height = Math.round(H * dpr);
    if (!canvas.style.aspectRatio) canvas.style.aspectRatio = W + " / " + H;
    var ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx: ctx, W: W, H: H, dpr: dpr };
  }

  function download(blob, name) {
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(function () { URL.revokeObjectURL(a.href); }, 1000);
  }
  function savePNG(canvas, name) { canvas.toBlob(function (b) { if (b) download(b, name || "plot.png"); }, "image/png"); }

  window.SiteTheme = {
    get: get, set: set, colors: colors, fit: fit, download: download, savePNG: savePNG,
    onChange: function (fn) { listeners.push(fn); },
  };

  function wire() {
    var b = document.getElementById("theme-toggle");
    if (b) b.addEventListener("click", function () { set(get() === "dark" ? "light" : "dark"); });
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", wire); else wire();
})();
