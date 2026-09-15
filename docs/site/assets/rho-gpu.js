/* rho-gpu.js -- the page side of the browser Pollard rho engine.
 *
 * Owns the consent and the readout, nothing else: the walk runs in
 * ./rho-gpu-worker.js on a Web Worker, over WebGPU.
 *
 * IT IS OFF UNTIL ASKED. A page that quietly spends a visitor's GPU is
 * indistinguishable from the thing browsers ship blocklists for, and a
 * throughput figure taken from an unaware visitor's laptop is not worth
 * having either. So: the control starts off, the choice is remembered in
 * localStorage, the device and the job are named before anything starts, and a
 * hidden tab pauses by default -- pausing holds the trails and the
 * distinguished-point table, so nothing is lost by it.
 */
(function () {
  "use strict";

  var root = document.getElementById("rho-gpu");
  if (!root) return;

  var el = function (id) { return document.getElementById(id); };
  var toggle = el("rho-run"), bg = el("rho-bg"), select = el("rho-job");
  var STORE = "rho-gpu.v1";

  var fmt = function (n) { return isFinite(n) ? Number(n).toLocaleString("en-US") : "—"; };
  var rate = function (n) {
    if (!isFinite(n) || n <= 0) return "—";
    if (n >= 1e9) return (n / 1e9).toFixed(2) + "G steps/s";
    if (n >= 1e6) return (n / 1e6).toFixed(2) + "M steps/s";
    if (n >= 1e3) return (n / 1e3).toFixed(1) + "k steps/s";
    return n.toFixed(0) + " steps/s";
  };
  var duration = function (s) {
    if (!isFinite(s) || s < 0) return "—";
    if (s < 90) return s.toFixed(0) + " s";
    if (s < 5400) return (s / 60).toFixed(1) + " min";
    if (s < 172800) return (s / 3600).toFixed(1) + " h";
    if (s < 3.15e9) return (s / 86400).toFixed(1) + " d";
    return (s / 3.156e7).toPrecision(3) + " y";
  };
  var set = function (id, text, cls) {
    var node = el(id);
    if (!node) return;
    node.textContent = text;
    if (cls !== undefined) node.className = "mono " + cls;
  };

  var prefs = {};
  try { prefs = JSON.parse(localStorage.getItem(STORE) || "{}") || {}; } catch (e) { prefs = {}; }
  var save = function () {
    try {
      localStorage.setItem(STORE, JSON.stringify({ on: !!toggle.checked, bg: !!bg.checked, job: select.value }));
    } catch (e) { /* private mode: the page works, the choice just is not kept */ }
  };

  if (!("gpu" in navigator)) {
    root.classList.add("unavailable");
    toggle.disabled = true;
    set("rho-state", "WebGPU not available in this browser", "s-bad");
    return;
  }

  var worker = null, job = null, expected = null, finished = false;

  var startWorker = function () {
    worker = new Worker("./assets/rho-gpu-worker.js", { type: "module" });
    worker.onmessage = function (ev) { handle(ev.data || {}); };
    worker.onerror = function (e) {
      set("rho-state", "worker failed: " + (e.message || "unknown"), "s-bad");
      toggle.checked = false;
      save();
    };
    worker.postMessage({ type: "jobs" });
  };

  var handle = function (m) {
    if (m.type === "jobs") {
      select.innerHTML = "";
      m.jobs.forEach(function (j) {
        var opt = document.createElement("option");
        opt.value = j.id;
        opt.textContent = j.label + " — " + j.bits + "-bit curve, rho ≈ 2^" + j.rho_log2;
        select.appendChild(opt);
      });
      if (prefs.job) select.value = prefs.job;
      if (!select.value && m.jobs.length) select.value = m.jobs[0].id;
      if (toggle.checked) start();
      return;
    }
    if (m.type === "error") { set("rho-state", m.message, "s-bad"); return; }
    if (m.type === "selftest") { set("rho-check", m.ok ? "passed — " + m.detail : "FAILED", m.ok ? "s-good" : "s-bad"); return; }
    if (m.type === "ready") {
      job = m.job;
      expected = Math.pow(2, m.job.rho_log2);
      set("rho-device", m.device);
      set("rho-shape", fmt(m.slots) + " walks, " + m.batch + " per thread, 2^" + m.r_bits + " jump table");
      set("rho-target", m.job.label + " — expected 2^" + m.job.rho_log2 + " steps, one point per 2^" + m.job.dp_bits);
      set("rho-state", "running", "s-good");
      return;
    }
    if (m.type === "stat") {
      set("rho-rate", m.paused ? "paused" : rate(m.rate));
      set("rho-steps", fmt(m.steps) + (expected ? "  (" + (100 * m.steps / expected).toFixed(2) + "% of 2^" + job.rho_log2 + ")" : ""));
      set("rho-dps", fmt(m.dps) + " reported, " + fmt(m.distinct) + " distinct");
      set("rho-eta", m.rate > 0 && expected ? duration(Math.max(0, expected - m.steps) / m.rate) : "—");
      set("rho-elapsed", duration(m.seconds));
      if (m.paused) set("rho-state", "paused (tab hidden)", "s-warn");
      else if (!m.solved) set("rho-state", "running", "s-good");
      if (m.dropped) set("rho-note", fmt(m.dropped) + " reports dropped by a full buffer; the engine shrinks its dispatch to recover", "s-warn");
      else if (m.table_full) set("rho-note", "distinguished-point table full — restart on a smaller curve to finish a solve", "s-warn");
      return;
    }
    if (m.type === "solved") {
      finished = true;
      set("rho-state", "solved: k = " + m.k + " (verified kG = Q)", "s-good");
      /* The last periodic stat is up to half a second stale; the solve carries
       * the step count the answer actually cost, so show that one. */
      set("rho-steps", fmt(m.steps) + (expected ? "  (" + (100 * m.steps / expected).toFixed(2) + "% of 2^" + job.rho_log2 + ")" : ""));
      set("rho-elapsed", duration(m.seconds));
      set("rho-note", "in " + fmt(m.steps) + " steps and " + duration(m.seconds)
        + "; the two trails that met were " + m.trails.join(" and ") + " steps long", "s-good");
      toggle.checked = false;
      save();
      return;
    }
    if (m.type === "stopped") {
      if (toggle.checked) return;   // a restart is already on its way
      set("rho-rate", m.rate > 0 ? rate(m.rate) + " at the end" : "—");
      /* A solved instance stops too; "stopped" must not overwrite the answer. */
      if (!finished) set("rho-state", "stopped", "");
    }
  };

  var start = function () {
    if (!worker) { startWorker(); return; }   // jobs message will start it
    finished = false;
    set("rho-state", "requesting a GPU adapter…", "");
    set("rho-check", "—", "");
    set("rho-note", "", "");
    worker.postMessage({ type: "start", job: select.value, paused: shouldPause() });
  };

  var shouldPause = function () { return document.hidden && !bg.checked; };

  toggle.checked = !!prefs.on;
  bg.checked = !!prefs.bg;

  toggle.addEventListener("change", function () {
    save();
    if (toggle.checked) start();
    else if (worker) worker.postMessage({ type: "stop" });
  });

  bg.addEventListener("change", function () {
    save();
    if (worker) worker.postMessage({ type: "pause", on: shouldPause() });
  });

  select.addEventListener("change", function () {
    save();
    if (toggle.checked && worker) {
      /* A job change is a different curve, table and instance: stop, then
       * start clean rather than carrying trails across curves. */
      worker.postMessage({ type: "stop" });
      setTimeout(start, 250);
    }
  });

  document.addEventListener("visibilitychange", function () {
    if (worker && toggle.checked) worker.postMessage({ type: "pause", on: shouldPause() });
  });

  window.addEventListener("pagehide", function () { if (worker) worker.postMessage({ type: "stop" }); });

  startWorker();
})();
