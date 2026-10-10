/* rho-gpu-worker.js -- the browser Pollard rho engine.
 *
 * Runs in a Web Worker so the walk never competes with the page for the main
 * thread, and drives ./rho-gpu.wgsl over WebGPU. The host arithmetic it needs
 * -- jump table, walk seeds, trail replay, collision-to-logarithm -- is in
 * ./rho-gpu-host.js; what is here is the device: buffers, dispatch pacing, the
 * self-test that establishes the device agrees with the host before a single
 * measurement is reported, and the reseeding of retired walk slots.
 *
 * WHAT THIS IS FOR. It measures how fast a browser GPU steps an r-adding walk,
 * on the scale-down curves this repository calibrates attacks against, and it
 * solves them end to end so the number is a measurement of a working solver
 * rather than of an inner loop. It is not aimed at a deployed curve, and the
 * ladder stops where a browser can still finish. A step rate here is a
 * throughput figure; the cost ledger's S = ops / sqrt(n) is what turns it into
 * a claim about an attack, and that conversion is done on the page.
 *
 * Protocol with the page (postMessage):
 *   in : {type:"jobs"} | {type:"start", job, workers} | {type:"pause", on}
 *      | {type:"stop"}
 *   out: {type:"jobs", jobs} | {type:"ready", device, job} | {type:"stat", ...}
 *      | {type:"dp", ...} | {type:"solved", ...} | {type:"error", message}
 *      | {type:"selftest", ok, detail} | {type:"stopped"}
 */
"use strict";

import {
  mod, invMod, limbs, fromLimbs, n0inv16, replay, solveCollision, prepare, JOBS,
  WG_SIZE, BATCH, R_BITS, R, DP_TABLE_CAP, TARGET_MS,
} from "./rho-gpu-host.js";

let stop = false;
let paused = false;
let running = false;

const post = (m) => self.postMessage(m);

self.onmessage = (ev) => {
  const msg = ev.data || {};
  if (msg.type === "jobs") {
    post({ type: "jobs", jobs: JOBS.map((j) => ({ id: j.id, label: j.label, bits: j.bits, rho_log2: j.rho_log2, dp_bits: j.dp_bits })) });
  } else if (msg.type === "start") {
    if (running) return;
    stop = false;
    paused = !!msg.paused;
    run(msg).catch((e) => { running = false; post({ type: "error", message: String(e && e.message || e) }); post({ type: "stopped" }); });
  } else if (msg.type === "pause") {
    /* Pausing holds the distinguished-point table and every live trail: a
     * hidden tab should cost nothing and lose nothing. */
    paused = !!msg.on;
  } else if (msg.type === "stop") {
    stop = true;
  }
};

/* -------------------------------------------------------------- the loop -- */

async function run(msg) {
  running = true;
  const job = JOBS.find((j) => j.id === msg.job) || JOBS[0];
  const st = prepare(job);
  if (st.L < 2) throw new Error("curve too small for the 32-bit distinguishing test");

  if (!navigator.gpu) throw new Error("WebGPU is not available in this browser");
  const adapter = await navigator.gpu.requestAdapter();
  if (!adapter) throw new Error("no WebGPU adapter (the GPU may be blocklisted)");
  const device = await adapter.requestDevice();
  device.lost.then((info) => { stop = true; post({ type: "error", message: "GPU device lost: " + (info && info.reason || "unknown") }); });

  const info = (adapter.info || {});
  const deviceLabel = [info.vendor, info.architecture || info.device].filter(Boolean).join(" ") || "WebGPU device";

  /* Walk count. More walks is not better: every walk is abandoned at its
   * first distinguished point, so the trails alive at any moment carry a
   * combined tail of slots * 2^dp_bits steps that contributes nothing to a
   * collision. Filling a big GPU with a small instance therefore spends most
   * of its work on that tail -- 8192 walks on toy32 overshoot the expected
   * 2^15.7 steps by nearly an order of magnitude. Cap the occupancy at the
   * expected collision cost divided by the distinguishing interval, which
   * keeps the tail at or below the search itself, then clamp to one workgroup
   * at the bottom and the device's real limit at the top. */
  const maxThreads = Math.min(adapter.limits.maxComputeWorkgroupsPerDimension * WG_SIZE, 1 << 16);
  const useful = Math.pow(2, job.rho_log2 - job.dp_bits);
  const want = Math.max(WG_SIZE * BATCH, Math.min((msg.workers | 0) || 8192, useful));
  const threads = Math.min(maxThreads, Math.ceil(want / BATCH / WG_SIZE) * WG_SIZE);
  const slots = threads * BATCH;

  const source = await (await fetch("./rho-gpu.wgsl")).text();
  const wgsl = source.replace("//__CONFIG__", [
    "const LIMBS      = " + st.L + "u;",
    "const ACC_WORDS  = " + (st.L + 2) + "u;",
    "const R          = " + R + "u;",
    "const BATCH      = " + BATCH + "u;",
    "const WG_SIZE    = " + WG_SIZE + "u;",
  ].join("\n"));
  const module = device.createShaderModule({ code: wgsl, label: "pollard-rho" });
  const compilation = await module.getCompilationInfo();
  const fatal = compilation.messages.filter((m) => m.type === "error");
  if (fatal.length) throw new Error("shader: " + fatal[0].message);
  const pipeline = await device.createComputePipelineAsync({ layout: "auto", compute: { module, entryPoint: "main" } });

  /* ---- buffers ---- */
  const L = st.L;
  const REPORT_WORDS = 4 + L;
  /* A slot can report at most one distinguished point per dispatch -- it goes
   * dead on the first -- so `slots` records is exactly enough and a report can
   * never be dropped. */
  const REPORT_CAP = slots;
  const PARAM_WORDS = 3 * L + 5;
  const paramsHost = new Uint32Array(PARAM_WORDS);
  paramsHost.set(limbs(st.p, L), 0);
  paramsHost.set(limbs(st.p - 2n, L), L);
  paramsHost.set(limbs(mod(st.Rmont, st.p), L), 2 * L);
  paramsHost[3 * L + 0] = n0inv16(st.p);
  paramsHost[3 * L + 1] = (st.p - 2n).toString(2).length;
  paramsHost[3 * L + 2] = 256;              // iters, retuned below
  paramsHost[3 * L + 3] = REPORT_CAP;
  paramsHost[3 * L + 4] = 0xffffffff;       // dp_mask: self-test first
  const ITERS_OFF = (3 * L + 2) * 4, DPMASK_OFF = (3 * L + 4) * 4;

  const S = (n) => n * 4;
  const usage = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC;
  const paramsBuf = device.createBuffer({ size: S(PARAM_WORDS), usage });
  const tableBuf = device.createBuffer({ size: S(2 * R * L), usage });
  const WALK_WORDS = 2 * L + 4;
  const walksBuf = device.createBuffer({ size: S(slots * WALK_WORDS), usage });
  const reportsBuf = device.createBuffer({ size: S(REPORT_CAP * REPORT_WORDS), usage });
  const countersBuf = device.createBuffer({ size: S(4), usage });
  const readCounters = device.createBuffer({ size: S(4), usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
  const readReports = device.createBuffer({ size: S(REPORT_CAP * REPORT_WORDS), usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
  const readWalks = device.createBuffer({ size: S(slots * WALK_WORDS), usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });

  const bind = device.createBindGroup({
    layout: pipeline.getBindGroupLayout(0),
    entries: [paramsBuf, tableBuf, walksBuf, reportsBuf, countersBuf].map((buffer, binding) => ({ binding, resource: { buffer } })),
  });

  device.queue.writeBuffer(paramsBuf, 0, paramsHost);

  const tableHost = new Uint32Array(2 * R * L);
  for (let j = 0; j < R; j++) {
    tableHost.set(limbs(st.toMont(st.table[j][0]), L), 2 * j * L);
    tableHost.set(limbs(st.toMont(st.table[j][1]), L), (2 * j + 1) * L);
  }
  device.queue.writeBuffer(tableBuf, 0, tableHost);

  /* Host mirror: what each slot was seeded with. The point itself lives on the
   * device; only the seed indices are needed to replay it. */
  const restart = new Uint32Array(slots);
  const walksHost = new Uint32Array(slots * WALK_WORDS);
  const seedSlot = (slot) => {
    const { P } = st.seedTrail(slot, restart[slot]);
    const off = slot * WALK_WORDS;
    walksHost.set(limbs(st.toMont(P[0]), L), off);
    walksHost.set(limbs(st.toMont(P[1]), L), off + L);
    walksHost[off + 2 * L + 0] = 1;             // live
    walksHost[off + 2 * L + 1] = 0;             // steps
    walksHost[off + 2 * L + 2] = restart[slot];
    walksHost[off + 2 * L + 3] = 0;
  };
  for (let s = 0; s < slots; s++) seedSlot(s);
  device.queue.writeBuffer(walksBuf, 0, walksHost);

  const dispatch = async (iters) => {
    device.queue.writeBuffer(paramsBuf, ITERS_OFF, new Uint32Array([iters]));
    device.queue.writeBuffer(countersBuf, 0, new Uint32Array(4));
    const enc = device.createCommandEncoder();
    const pass = enc.beginComputePass();
    pass.setPipeline(pipeline);
    pass.setBindGroup(0, bind);
    pass.dispatchWorkgroups(Math.ceil(threads / WG_SIZE));
    pass.end();
    enc.copyBufferToBuffer(countersBuf, 0, readCounters, 0, S(4));
    enc.copyBufferToBuffer(reportsBuf, 0, readReports, 0, readReports.size);
    enc.copyBufferToBuffer(walksBuf, 0, readWalks, 0, readWalks.size);
    device.queue.submit([enc.finish()]);
    await Promise.all([
      readCounters.mapAsync(GPUMapMode.READ),
      readReports.mapAsync(GPUMapMode.READ),
      readWalks.mapAsync(GPUMapMode.READ),
    ]);
    const counters = new Uint32Array(readCounters.getMappedRange().slice(0));
    const reports = new Uint32Array(readReports.getMappedRange().slice(0));
    const walks = new Uint32Array(readWalks.getMappedRange().slice(0));
    readCounters.unmap(); readReports.unmap(); readWalks.unmap();
    return { counters, reports, walks };
  };

  /* ---- self-test: 32 steps on the device against 32 steps on the host ----
   * Distinguishing is off (dp_mask all ones needs the top 24 bits of x clear),
   * so every slot survives and every slot is comparable. A mismatch stops the
   * engine: a walk whose device and host steps differ can still produce
   * collisions, and every one of them would be meaningless. */
  const TEST_STEPS = 32;
  {
    const { walks } = await dispatch(TEST_STEPS);
    let bad = null;
    for (let s = 0; s < Math.min(slots, 16) && bad === null; s++) {
      const off = s * WALK_WORDS;
      const got = [fromLimbs(walks, off, L), fromLimbs(walks, off + L, L)];
      let { P } = st.seedTrail(s, 0);
      for (let i = 0; i < walks[off + 2 * L + 1]; i++) P = st.E.add(P, st.table[st.partition(P)]);
      if (st.toMont(P[0]) !== got[0] || st.toMont(P[1]) !== got[1]) bad = s;
    }
    if (bad !== null) throw new Error("self-test failed: device and host disagree on slot " + bad);
    post({ type: "selftest", ok: true, detail: TEST_STEPS + " steps x " + Math.min(slots, 16) + " walks reproduced in BigInt" });
  }

  /* Reseed everything, switch distinguishing on, and start measuring. */
  for (let s = 0; s < slots; s++) { restart[s] = 1; seedSlot(s); }
  device.queue.writeBuffer(walksBuf, 0, walksHost);
  device.queue.writeBuffer(paramsBuf, DPMASK_OFF, new Uint32Array([st.dpMask]));

  post({
    type: "ready", device: deviceLabel,
    job: { id: job.id, label: job.label, bits: job.bits, rho_log2: job.rho_log2, dp_bits: job.dp_bits },
    slots, batch: BATCH, r_bits: R_BITS,
  });

  const dps = new Map();      // x (hex) -> {slot, restart, steps}
  const ITERS_CAP = Math.max(64, Math.min(1 << 20, Math.pow(2, job.dp_bits - 1)));
  let iters = 64, steps = 0, reported = 0, dropped = 0, solved = null, full = false;
  const t0 = performance.now();
  let lastStat = 0, rate = 0;

  while (!stop) {
    if (paused) {
      post({ type: "stat", steps, rate: 0, dps: reported, distinct: dps.size, dropped, seconds: (performance.now() - t0) / 1000, iters, reseeded: 0, table_full: full, paused: true, solved: false });
      await new Promise((r) => setTimeout(r, 250));
      continue;
    }
    const t1 = performance.now();
    const { counters, reports, walks } = await dispatch(iters);
    const dt = performance.now() - t1;

    /* dt spans the submit, the walk and the readback, so `rate` is a sustained
     * rate rather than a kernel time: it is what the page reports, and what an
     * expected-time-left figure has to be divided by. */
    steps += counters[1];
    dropped += counters[2];
    const n = Math.min(counters[0], REPORT_CAP);
    reported += n;
    rate = counters[1] / Math.max(dt, 1) * 1000;

    /* Retune the dispatch so one submit is about TARGET_MS of work: the tab
     * stays responsive and a stop request is honoured inside one dispatch.
     *
     * The cap is not about latency. A slot retires on its first distinguished
     * point and only the host can reseed it, so a dispatch longer than the
     * expected 2^dp_bits steps between points leaves most slots idle for most
     * of it. Half that expectation keeps about 80% of the walks stepping.
     */
    const scale = TARGET_MS / Math.max(dt, 1);
    iters = Math.max(16, Math.min(ITERS_CAP, Math.round(iters * Math.min(4, Math.max(0.25, scale)))));

    for (let i = 0; i < n && !solved; i++) {
      const off = i * REPORT_WORDS;
      const slot = reports[off], rst = reports[off + 1], stp = reports[off + 2];
      const xHex = fromLimbs(reports, off + 4, L).toString(16);
      const prev = dps.get(xHex);
      if (prev && !(prev.slot === slot && prev.restart === rst)) {
        solved = await solveCollision(st, prev, { slot, restart: rst, steps: stp }, xHex);
        if (solved) post({ type: "solved", ...solved, job: job.id, steps, seconds: (performance.now() - t0) / 1000 });
        else post({ type: "dp", note: "collision did not replay; discarded", x: xHex });
      } else if (!prev && dps.size < DP_TABLE_CAP) {
        dps.set(xHex, { slot, restart: rst, steps: stp });
      } else if (!prev) {
        full = true;      // the table, not the search, is what ran out
      }
    }

    /* Reseed the slots the dispatch retired (a distinguished point, or the
     * ~2/p exceptional addition). */
    let reseeded = 0;
    for (let s = 0; s < slots; s++) {
      if (walks[s * WALK_WORDS + 2 * L] !== 0) continue;
      restart[s] += 1;
      seedSlot(s);
      device.queue.writeBuffer(walksBuf, S(s * WALK_WORDS), walksHost, s * WALK_WORDS, WALK_WORDS);
      reseeded++;
    }

    if (performance.now() - lastStat > 500) {
      lastStat = performance.now();
      post({
        type: "stat", steps, rate, dps: reported, distinct: dps.size, dropped,
        seconds: (performance.now() - t0) / 1000, iters, reseeded,
        table_full: full, solved: !!solved,
      });
    }
    if (solved) break;
    /* Yield to the event loop so a stop lands promptly. */
    await new Promise((r) => setTimeout(r, 0));
  }

  running = false;
  post({ type: "stopped", steps, seconds: (performance.now() - t0) / 1000, rate });
}

