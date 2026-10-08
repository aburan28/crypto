/* rho_gpu_emulate.mjs -- tests for the browser WebGPU Pollard rho engine
 * (docs/site/assets/rho-gpu.{wgsl,js}, rho-gpu-host.js, rho-gpu-worker.js).
 *
 *   node scripts/rho_gpu_emulate.mjs
 *
 * WHAT THIS ESTABLISHES, AND WHAT IT DOES NOT. The shader's arithmetic is
 * re-implemented below in u32 JavaScript, mirroring rho-gpu.wgsl function by
 * function, and checked against BigInt on every job curve: the 16-bit-limb
 * CIOS Montgomery multiply and its carry bounds, the conditional add/subtract,
 * Fermat inversion, the batched inverse, the partition and distinguishing
 * predicates, and the full step. Then a complete ECDLP is solved on toy32
 * through the emulated device and the engine's own replay/solve path, which
 * covers the part no unit test reaches: that a collision between two device
 * trails reported as (slot, restart, steps, x) recovers the logarithm.
 *
 * What it cannot establish is that rho-gpu.wgsl -- the text that actually runs
 * on a GPU -- matches this emulation, because nothing here compiles WGSL. That
 * gap is closed at runtime instead: the worker replays its first 32 device
 * steps in BigInt before it reports a single measurement, and refuses to run
 * if they disagree. No measurement from this engine is published from an
 * unverified device.
 */
import {
  mod, invMod, limbs, fromLimbs, n0inv16, prepare, replay, solveCollision,
  JOBS, BATCH, R,
} from "../docs/site/assets/rho-gpu-host.js";

let failures = 0, checks = 0;
const ok = (cond, what) => {
  checks++;
  if (!cond) { failures++; console.log("  FAIL  " + what); }
};

/* ----------------------------------------------- the shader, in u32 JS ----
 * Every operation below is the WGSL function of the same name. u32 wrap is
 * explicit (>>> 0); the products are of 16-bit values, so they are exact in a
 * double and the >>> 0 is a formality that also documents the bound.
 */
function Device(p, L) {
  const P = limbs(p, L);
  const PM2 = limbs(p - 2n, L);
  const N0 = n0inv16(p);
  const PM2_BITS = (p - 2n).toString(2).length;
  const ONE = limbs(mod(1n << (16n * BigInt(L)), p), L);

  const ge = (a, b) => {
    for (let i = L - 1; i >= 0; i--) if (a[i] !== b[i]) return a[i] > b[i];
    return true;
  };
  const subRaw = (a, b) => {
    const r = new Uint32Array(L);
    let borrow = 0;
    for (let i = 0; i < L; i++) {
      const t = (a[i] + 0x10000 - b[i] - borrow) >>> 0;
      r[i] = t & 0xffff;
      borrow = 1 - (t >>> 16);
    }
    return r;
  };
  const addRawP = (a) => {
    const r = new Uint32Array(L);
    let carry = 0;
    for (let i = 0; i < L; i++) {
      const t = a[i] + P[i] + carry;
      ok(t <= 0xffffffff, "fe_add_raw_p stays in u32");
      r[i] = t & 0xffff;
      carry = t >>> 16;
    }
    ok(carry === 0, "a + p does not overflow the representation (p < 2^(16L-1))");
    return r;
  };
  const sub = (a, b) => (ge(a, b) ? subRaw(a, b) : subRaw(addRawP(a), b));
  const isZero = (a) => a.every((v) => v === 0);

  const montMul = (a, b) => {
    const acc = new Uint32Array(L + 2);
    for (let i = 0; i < L; i++) {
      let c = 0;
      for (let j = 0; j < L; j++) {
        const t = acc[j] + a[j] * b[i] + c;
        ok(t <= 0xffffffff, "CIOS accumulation stays in u32");
        acc[j] = t & 0xffff;
        c = t >>> 16;
      }
      let t = acc[L] + c;
      acc[L] = t & 0xffff;
      acc[L + 1] = t >>> 16;

      const m = (acc[0] * N0) & 0xffff;
      t = acc[0] + m * P[0];
      c = t >>> 16;
      for (let j = 1; j < L; j++) {
        const u = acc[j] + m * P[j] + c;
        ok(u <= 0xffffffff, "CIOS reduction stays in u32");
        acc[j - 1] = u & 0xffff;
        c = u >>> 16;
      }
      t = acc[L] + c;
      acc[L - 1] = t & 0xffff;
      acc[L] = acc[L + 1] + (t >>> 16);
    }
    let r = acc.slice(0, L);
    if (acc[L] !== 0 || ge(r, P)) r = subRaw(r, P);
    return r;
  };

  const inv = (a) => {
    let r = Uint32Array.from(ONE);
    for (let i = PM2_BITS; i > 0; i--) {
      const bit = i - 1;
      r = montMul(r, r);
      if (((PM2[bit >> 4] >>> (bit & 15)) & 1) === 1) r = montMul(r, a);
    }
    return r;
  };

  /* The shader's batched inverse, over the same BATCH slots. */
  const batchInv = (d) => {
    const pref = [];
    let run = d[0];
    pref[0] = run;
    for (let s = 1; s < BATCH; s++) { run = montMul(run, d[s]); pref[s] = run; }
    let acc = inv(run);
    for (let i = BATCH; i > 1; i--) {
      const s = i - 1;
      const invS = montMul(acc, pref[s - 1]);
      acc = montMul(acc, d[s]);
      d[s] = invS;
    }
    d[0] = acc;
    return d;
  };

  const partition = (x) => x[0] & (R - 1);
  const isDp = (x, mask) => ((((x[0] | (x[1] << 16)) >>> 0) >>> 8) & mask) === 0;

  return { P, ONE, ge, sub, subRaw, isZero, montMul, inv, batchInv, partition, isDp, L };
}

/* -------------------------------------------------------------- the runs -- */

function fieldChecks(job) {
  const st = prepare(job);
  const dev = Device(st.p, st.L);
  const L = st.L;
  const toL = (x) => limbs(x, L);
  const back = (arr) => fromLimbs(arr, 0, L);
  let seed = 12345n;
  const nextRand = () => { seed = (seed * 6364136223846793005n + 1442695040888963407n) & ((1n << 64n) - 1n); return mod(seed, st.p); };

  for (let t = 0; t < 40; t++) {
    const x = nextRand(), y = nextRand();
    const xm = st.toMont(x), ym = st.toMont(y);
    ok(back(dev.montMul(toL(xm), toL(ym))) === st.toMont(mod(x * y, st.p)),
      job.id + ": mont_mul agrees with BigInt");
    ok(back(dev.sub(toL(xm), toL(ym))) === st.toMont(mod(x - y, st.p)),
      job.id + ": fe_sub agrees with BigInt");
    if (x !== 0n) {
      ok(back(dev.inv(toL(xm))) === st.toMont(invMod(x, st.p)),
        job.id + ": fe_inv agrees with BigInt");
    }
  }

  /* Batched inverse over BATCH random denominators. */
  const xs = [], d = [];
  for (let s = 0; s < BATCH; s++) { const v = nextRand() || 1n; xs.push(v); d.push(toL(st.toMont(v))); }
  dev.batchInv(d);
  for (let s = 0; s < BATCH; s++) {
    ok(back(d[s]) === st.toMont(invMod(xs[s], st.p)), job.id + ": batched inverse slot " + s);
  }

  /* One emulated device step against one BigInt step, on the real table. */
  let Pt = st.seedTrail(0, 0).P;
  let cur = [toL(st.toMont(Pt[0])), toL(st.toMont(Pt[1]))];
  for (let i = 0; i < 64; i++) {
    const j = dev.partition(cur[0]);
    ok(j === st.partition(Pt), job.id + ": partition agrees at step " + i);
    ok(dev.isDp(cur[0], st.dpMask) === st.isDp(Pt), job.id + ": dp predicate agrees at step " + i);
    const tx = toL(st.toMont(st.table[j][0])), ty = toL(st.toMont(st.table[j][1]));
    const den = dev.sub(tx, cur[0]);
    ok(!dev.isZero(den), job.id + ": step " + i + " is a generic addition");
    const lam = dev.montMul(dev.sub(ty, cur[1]), dev.inv(den));
    const x3 = dev.sub(dev.sub(dev.montMul(lam, lam), cur[0]), tx);
    const y3 = dev.sub(dev.montMul(lam, dev.sub(cur[0], x3)), cur[1]);
    cur = [x3, y3];
    Pt = st.E.add(Pt, st.table[j]);
    ok(back(x3) === st.toMont(Pt[0]) && back(y3) === st.toMont(Pt[1]),
      job.id + ": emulated step " + i + " matches the BigInt walk");
  }
}

/* End-to-end: run the emulated device until two trails meet at a distinguished
 * point, then solve through the engine's own replay path. */
async function endToEnd(job, slots = 64, maxSteps = 1 << 22) {
  const st = prepare(job);
  const dev = Device(st.p, st.L);
  const L = st.L, toL = (x) => limbs(x, L), back = (a) => fromLimbs(a, 0, L);

  const restart = new Uint32Array(slots);
  const state = [], steps = new Uint32Array(slots);
  for (let s = 0; s < slots; s++) {
    const P = st.seedTrail(s, 0).P;
    state.push([toL(st.toMont(P[0])), toL(st.toMont(P[1]))]);
  }

  const dps = new Map();
  let total = 0, solved = null;
  while (!solved && total < maxSteps) {
    for (let s = 0; s < slots && !solved; s++) {
      const j = dev.partition(state[s][0]);
      const tx = toL(st.toMont(st.table[j][0])), ty = toL(st.toMont(st.table[j][1]));
      const den = dev.sub(tx, state[s][0]);
      if (dev.isZero(den)) { restart[s]++; steps[s] = 0; const P = st.seedTrail(s, restart[s]).P; state[s] = [toL(st.toMont(P[0])), toL(st.toMont(P[1]))]; continue; }
      const lam = dev.montMul(dev.sub(ty, state[s][1]), dev.inv(den));
      const x3 = dev.sub(dev.sub(dev.montMul(lam, lam), state[s][0]), tx);
      const y3 = dev.sub(dev.montMul(lam, dev.sub(state[s][0], x3)), state[s][1]);
      state[s] = [x3, y3];
      steps[s]++; total++;
      if (!dev.isDp(x3, st.dpMask)) continue;

      const xHex = back(x3).toString(16);
      const rec = { slot: s, restart: restart[s], steps: steps[s] };
      const prev = dps.get(xHex);
      if (prev && !(prev.slot === rec.slot && prev.restart === rec.restart)) {
        solved = await solveCollision(st, prev, rec, xHex);
        ok(solved !== null, job.id + ": the collision at " + xHex + " solved");
      } else if (!prev) {
        dps.set(xHex, rec);
      }
      restart[s]++; steps[s] = 0;
      const P = st.seedTrail(s, restart[s]).P;
      state[s] = [toL(st.toMont(P[0])), toL(st.toMont(P[1]))];
    }
  }
  ok(solved !== null, job.id + ": solved end to end within " + maxSteps + " steps");
  if (solved) {
    const k = BigInt(solved.k);
    const check = st.E.mul(k, st.G);
    ok(check[0] === st.Q[0] && check[1] === st.Q[1], job.id + ": kG == Q");
    console.log("  " + job.id + ": k = " + k + "  (" + total + " emulated steps, "
      + dps.size + " distinguished points, trails " + solved.trails.join(" and ") + ")");
  }
}

console.log("field and walk agreement, per job curve:");
for (const job of JOBS) fieldChecks(job);
console.log("end-to-end solve through the emulated device:");
await endToEnd(JOBS[0]);

console.log(checks + " checks, " + failures + " failures");
process.exit(failures === 0 ? 0 : 1);
