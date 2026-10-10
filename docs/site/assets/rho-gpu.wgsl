// rho-gpu.wgsl -- Pollard rho r-adding walk for WebGPU.
//
// The walk definition is the one in gpu/ecc/rho.cuh and gpu/ecc/ecref.py, so a
// browser step and a CUDA step are the same step:
//
//   partition(P)  = limb0(internal(x)) & (R - 1)
//   step          P <- P + M[partition(P)],  M[j] = c_j P + d_j Q
//   distinguished ((low32(internal(x)) >> 8) & DP_MASK) == 0
//
// `internal` is Montgomery form, as on the CUDA generic path: hashing the
// internal representation saves a conversion per step and is still a
// deterministic function of the point, which is all the walk needs. The host
// side (rho-gpu-worker.js) reproduces both predicates bit-for-bit in BigInt,
// which is what makes a device trail replayable and a collision solvable.
//
// The negation map is NOT implemented. It would buy sqrt(2) steps and cost the
// fruitless-cycle machinery that rho.cuh needs for it; this engine measures
// browser throughput, so the simpler walk is the honest one to measure.
//
// Field arithmetic uses 16-BIT limbs. WGSL has no 64-bit integer type and no
// widening multiply, so a 32x32 product cannot be formed at all; with 16-bit
// limbs every partial product, plus an accumulator word and a carry, fits in
// u32 exactly (see mont_mul). LIMBS, R_BITS, DP_MASK, BATCH and the modulus
// are substituted into this source before compilation.

//__CONFIG__

alias Fe = array<u32, LIMBS>;

struct Params {
  p        : Fe,          // modulus, 16-bit limbs, little-endian
  pm2      : Fe,          // p - 2, the Fermat inversion exponent
  one_mont : Fe,          // R mod p
  n0inv    : u32,         // -p^-1 mod 2^16
  pm2_bits : u32,         // bit length of p - 2
  iters    : u32,         // walk steps per dispatch
  out_cap  : u32,         // capacity of the report buffer, in records
  dp_mask  : u32,         // distinguished iff ((low32(x) >> 8) & dp_mask) == 0
};

// One walk slot: the current point, plus a liveness flag. A slot goes dead on
// a reported distinguished point or an exceptional addition; the host reseeds
// dead slots from fresh (a, b) coefficients before the next dispatch.
struct Walk {
  x      : Fe,
  y      : Fe,
  live   : u32,
  steps  : u32,           // steps since this slot was seeded
  restart: u32,           // how many times the host has reseeded it
  pad    : u32,
};

// A report is {slot, restart, steps, x}: enough for the host to replay the
// trail and recover its (a, b) coefficients. Coordinates stay on the device
// otherwise.
struct Report {
  slot    : u32,
  restart : u32,
  steps   : u32,
  pad     : u32,
  x       : Fe,
};

@group(0) @binding(0) var<storage, read>       params  : Params;
@group(0) @binding(1) var<storage, read>       table   : array<Fe>;   // 2*R entries: x0,y0,x1,y1,...
@group(0) @binding(2) var<storage, read_write> walks   : array<Walk>;
@group(0) @binding(3) var<storage, read_write> reports : array<Report>;
@group(0) @binding(4) var<storage, read_write> counters: array<atomic<u32>>; // 0: reports, 1: steps done, 2: reports dropped

// ---------------------------------------------------------------- field ----

fn fe_is_zero(a: Fe) -> bool {
  var acc = 0u;
  for (var i = 0u; i < LIMBS; i++) { acc |= a[i]; }
  return acc == 0u;
}

// a >= b, unsigned, limb-wise from the top.
fn fe_ge(a: Fe, b: Fe) -> bool {
  for (var i = LIMBS; i > 0u; i--) {
    let k = i - 1u;
    if (a[k] != b[k]) { return a[k] > b[k]; }
  }
  return true;
}

fn fe_sub_raw(a: Fe, b: Fe) -> Fe {
  var r: Fe;
  var borrow = 0u;
  for (var i = 0u; i < LIMBS; i++) {
    let t = a[i] + 0x10000u - b[i] - borrow;
    r[i] = t & 0xffffu;
    borrow = 1u - (t >> 16u);
  }
  return r;
}

fn fe_add_raw_p(a: Fe) -> Fe {
  var r: Fe;
  var carry = 0u;
  for (var i = 0u; i < LIMBS; i++) {
    let t = a[i] + params.p[i] + carry;
    r[i] = t & 0xffffu;
    carry = t >> 16u;
  }
  return r;   // a + p < 2^(16*LIMBS) because p < 2^(16*LIMBS - 1)
}

fn fe_add(a: Fe, b: Fe) -> Fe {
  var r: Fe;
  var carry = 0u;
  for (var i = 0u; i < LIMBS; i++) {
    let t = a[i] + b[i] + carry;
    r[i] = t & 0xffffu;
    carry = t >> 16u;
  }
  if (carry != 0u || fe_ge(r, params.p)) { r = fe_sub_raw(r, params.p); }
  return r;
}

fn fe_sub(a: Fe, b: Fe) -> Fe {
  if (fe_ge(a, b)) { return fe_sub_raw(a, b); }
  return fe_sub_raw(fe_add_raw_p(a), b);
}

fn mont_mul(a: Fe, b: Fe) -> Fe {
  var acc: array<u32, ACC_WORDS>;
  for (var i = 0u; i < LIMBS + 2u; i++) { acc[i] = 0u; }

  for (var i = 0u; i < LIMBS; i++) {
    var c = 0u;
    for (var j = 0u; j < LIMBS; j++) {
      let t = acc[j] + a[j] * b[i] + c;
      acc[j] = t & 0xffffu;
      c = t >> 16u;
    }
    var t = acc[LIMBS] + c;
    acc[LIMBS] = t & 0xffffu;
    acc[LIMBS + 1u] = t >> 16u;

    let m = (acc[0] * params.n0inv) & 0xffffu;
    t = acc[0] + m * params.p[0];
    c = t >> 16u;
    for (var j = 1u; j < LIMBS; j++) {
      let u = acc[j] + m * params.p[j] + c;
      acc[j - 1u] = u & 0xffffu;
      c = u >> 16u;
    }
    t = acc[LIMBS] + c;
    acc[LIMBS - 1u] = t & 0xffffu;
    acc[LIMBS] = acc[LIMBS + 1u] + (t >> 16u);
  }

  var r: Fe;
  for (var i = 0u; i < LIMBS; i++) { r[i] = acc[i]; }
  if (acc[LIMBS] != 0u || fe_ge(r, params.p)) { r = fe_sub_raw(r, params.p); }
  return r;
}

fn mont_sqr(a: Fe) -> Fe { return mont_mul(a, a); }

// a^(p-2) mod p: one inversion per BATCH walks, so its cost is amortised by
// the batched-inverse trick below rather than paid per step.
fn fe_inv(a: Fe) -> Fe {
  var r = params.one_mont;
  var i = params.pm2_bits;
  loop {
    if (i == 0u) { break; }
    i--;
    r = mont_sqr(r);
    if (((params.pm2[i >> 4u] >> (i & 15u)) & 1u) == 1u) { r = mont_mul(r, a); }
  }
  return r;
}

// ------------------------------------------------------------- the walk ----

fn walk_partition(x: Fe) -> u32 { return x[0] & (R - 1u); }

fn is_dp(x: Fe) -> bool {
  let low32 = x[0] | (x[1] << 16u);
  return ((low32 >> 8u) & params.dp_mask) == 0u;
}

@compute @workgroup_size(WG_SIZE)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let base = gid.x * BATCH;
  if (base >= arrayLength(&walks)) { return; }

  // Registers for the BATCH slots this thread owns.
  var px: array<Fe, BATCH>;
  var py: array<Fe, BATCH>;
  var live: array<u32, BATCH>;
  var steps: array<u32, BATCH>;

  var steps_in = 0u;
  for (var s = 0u; s < BATCH; s++) {
    let w = walks[base + s];
    px[s] = w.x; py[s] = w.y; live[s] = w.live; steps[s] = w.steps;
    steps_in += w.steps;
  }

  var d: array<Fe, BATCH>;      // x2 - x1 per slot, the denominators
  var num: array<Fe, BATCH>;    // y2 - y1 per slot
  var tx: array<Fe, BATCH>;     // x2 per slot
  var pref: array<Fe, BATCH>;   // running products for the batched inverse

  for (var it = 0u; it < params.iters; it++) {
    var any = 0u;
    for (var s = 0u; s < BATCH; s++) {
      if (live[s] == 0u) { d[s] = params.one_mont; continue; }
      any = 1u;
      let j = walk_partition(px[s]);
      tx[s] = table[2u * j];
      let ty = table[2u * j + 1u];
      d[s] = fe_sub(tx[s], px[s]);
      num[s] = fe_sub(ty, py[s]);
      // x1 == x2 is a doubling or an inverse pair: probability ~2/p per step,
      // and handling it needs a different formula. Retire the slot instead and
      // let the host reseed it; the trail is abandoned, not miscomputed.
      if (fe_is_zero(d[s])) { live[s] = 0u; d[s] = params.one_mont; }
    }
    if (any == 0u) { break; }

    // Montgomery's batched inverse: one fe_inv for BATCH denominators.
    var run = d[0];
    pref[0] = run;
    for (var s = 1u; s < BATCH; s++) { run = mont_mul(run, d[s]); pref[s] = run; }
    var acc = fe_inv(run);
    for (var i = BATCH; i > 1u; i--) {
      let s = i - 1u;
      let inv_s = mont_mul(acc, pref[s - 1u]);
      acc = mont_mul(acc, d[s]);
      d[s] = inv_s;
    }
    d[0] = acc;

    for (var s = 0u; s < BATCH; s++) {
      if (live[s] == 0u) { continue; }
      let lam = mont_mul(num[s], d[s]);
      let x3 = fe_sub(fe_sub(mont_sqr(lam), px[s]), tx[s]);
      let y3 = fe_sub(mont_mul(lam, fe_sub(px[s], x3)), py[s]);
      px[s] = x3; py[s] = y3;
      steps[s] = steps[s] + 1u;
      if (is_dp(x3)) {
        let slot = base + s;
        let idx = atomicAdd(&counters[0], 1u);
        if (idx < params.out_cap) {
          reports[idx].slot = slot;
          reports[idx].restart = walks[slot].restart;
          reports[idx].steps = steps[s];
          reports[idx].x = x3;
        } else {
          atomicAdd(&counters[2], 1u);   // dropped: host shrinks `iters`
        }
        live[s] = 0u;                    // host reseeds this slot
      }
    }
  }

  // counters[1] accumulates the steps this dispatch actually performed, which
  // is what a steps/second figure has to be measured from: a slot that went
  // dead early did fewer steps than `iters`.
  var steps_out = 0u;
  for (var s = 0u; s < BATCH; s++) { steps_out += steps[s]; }
  atomicAdd(&counters[1], steps_out - steps_in);

  for (var s = 0u; s < BATCH; s++) {
    walks[base + s].x = px[s];
    walks[base + s].y = py[s];
    walks[base + s].live = live[s];
    walks[base + s].steps = steps[s];
  }
}
