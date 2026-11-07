/*
 * twodes-worker.js: runs the 2DES engine (twodes-core.js) off the main thread so the page stays
 * responsive during large transforms. Messages: {id, kind: "run" | "trace" | "convergence", job, extra}.
 */
/* global importScripts, TwoDES */
importScripts("twodes-core.js");

/** Make a run() result structured-cloneable (pathway waiting-time factors are functions) and list transferable buffers. */
function pack(r, T) {
  const transfer = [];
  const add = (S) => { if (S) { transfer.push(S.re.buffer, S.im.buffer); } return S; };
  [r.R, r.NR, r.total].forEach(add);
  if (r.paths) r.paths = r.paths.map((p) => ({ kind: p.kind, path: p.path, amp: p.amp, w1: p.w1, w3: p.w3, g1: p.g1, g3: p.g3, label: p.label, tfT: p.tf(T) }));
  return { r, transfer };
}

self.onmessage = (e) => {
  const { id, kind, job, extra } = e.data;
  try {
    if (kind === "run") {
      const { r, transfer } = pack(TwoDES.run(job), job.T);
      self.postMessage({ id, ok: true, result: r }, transfer);
    } else if (kind === "trace") {
      self.postMessage({ id, ok: true, result: TwoDES.waitingTrace(job, extra.w1, extra.w3, extra.Ts) });
    } else if (kind === "cmPrepare") {
      const info = TwoDES.cmPrepare(job, extra);
      self.postMessage({ id, ok: true, result: info });
    } else if (kind === "cmMaps") {
      const m = TwoDES.cmMaps(extra.wT);
      self.postMessage({ id, ok: true, result: { real: m.real, pos: m.pos, neg: m.neg, w1: m.w1.slice(), w3: m.w3.slice(), na: m.na, nb: m.nb, tau: m.tau } }, [m.real.buffer, m.pos.buffer, m.neg.buffer]);
    } else if (kind === "convergence") {
      self.postMessage({ id, ok: true, result: TwoDES.convergence(job, extra.range) });
    }
  } catch (err) {
    self.postMessage({ id, ok: false, error: String((err && err.message) || err) });
  }
};
