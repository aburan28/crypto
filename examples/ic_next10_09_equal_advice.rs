//! IC next-10 item 9: equal-advice audit of the Koblitz online IC wins and the
//! E-3 free-relation census.
//!
//! Public synthetic Koblitz fixtures only (n <= 63, single-word fields).
//!
//! Modes:
//!
//! * `targets <n> <a> <count> <label> <points_out> <scalars_out>`:
//!   deterministic unseen public targets `Q = [d]G`, `d = H(label, i) mod r`.
//!   The points file holds one `[x,y]` per line (what the measured arms read);
//!   the scalar file is kept only for the verifier.
//! * `bl <n> <a> <table_bytes> <targets.jsonl-or-points> <out.jsonl>`:
//!   Bernstein–Lange preprocessing rho on signed-Frobenius classes. The
//!   precomputation stores `T = table_bytes / 16` distinguished points
//!   `(canonical x, log)` found by `G`-only walks; the online walk from
//!   `Q + [c]G` uses the *same* iteration map. Distinguished-point density
//!   `2^-t` with `t = round(log2 sqrt(N / T))`, `N = r / (2n)`. Prefix tables
//!   (the first `T/4^j` distinct DPs) are evaluated on the same targets from
//!   the same precomputation. Env: `BL_ONLINE_MAP_SEED` (negative control: a
//!   different jump table online), `BL_PREFIXES` (default 3), `BL_SEED`.
//! * `census <n> <a> <K> <out.jsonl>`: builds the compact-orbit factor base
//!   exactly as `koblitz_orbit_dlp_fast_online.rs` `construct:<n>:<a>:<K>`,
//!   enumerates every unordered pair state `(l <= r, relative)` and both S3
//!   roots, sorts canonical (Frobenius-rotation-minimal) roots, verifies every
//!   collision as a four-point group identity, and reports the rank mod r of
//!   the free homogeneous relation rows.
//!
//! The base construction, normal basis and S3 solver are copied from
//! `examples/koblitz_orbit_dlp_fast_online.rs` (private items there).

#![recursion_limit = "256"]

use crypto_lib::cryptanalysis::koblitz_fast_arith::{s3_x_roots, FastBinaryCurve, FastPoint};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::io::Write;
use std::time::Instant;

// ---------------------------------------------------------------- field helpers

struct Linear {
    tables: Vec<[u64; 256]>,
}

impl Linear {
    fn from_images(images: &[u64]) -> Self {
        let chunks = images.len().div_ceil(8);
        let mut tables = vec![[0u64; 256]; chunks];
        for (chunk, table) in tables.iter_mut().enumerate() {
            for byte in 1usize..256 {
                let low = byte & (byte - 1);
                let bit = chunk * 8 + byte.trailing_zeros() as usize;
                let image = images.get(bit).copied().unwrap_or(0);
                table[byte] = table[low] ^ image;
            }
        }
        Self { tables }
    }

    #[inline(always)]
    fn apply(&self, v: u64) -> u64 {
        let mut acc = 0u64;
        for (chunk, table) in self.tables.iter().enumerate() {
            acc ^= table[((v >> (8 * chunk)) & 0xff) as usize];
        }
        acc
    }
}

struct NormalBasis {
    n: u32,
    mask: u64,
    to_normal: Linear,
    to_poly: Linear,
}

impl NormalBasis {
    fn new(gf: &Gf2) -> Self {
        let n = gf.n as usize;
        let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        for candidate in 2u64.. {
            let mut conjugates = Vec::with_capacity(n);
            let mut value = candidate & mask;
            for _ in 0..n {
                conjugates.push(value);
                value = gf.sqr(value);
            }
            let Some(inverse_columns) = invert_columns(&conjugates, n) else {
                continue;
            };
            return Self {
                n: gf.n,
                mask,
                to_normal: Linear::from_images(&inverse_columns),
                to_poly: Linear::from_images(&conjugates),
            };
        }
        unreachable!()
    }

    #[inline(always)]
    fn rotate(&self, v: u64, k: u32) -> u64 {
        if k == 0 {
            return v;
        }
        ((v << k) | (v >> (self.n - k))) & self.mask
    }

    #[inline(always)]
    fn canonical(&self, v: u64) -> (u64, u32) {
        let mut best = v;
        let mut best_shift = 0;
        let mut current = v;
        for shift in 1..self.n {
            current = ((current << 1) | (current >> (self.n - 1))) & self.mask;
            if current < best {
                best = current;
                best_shift = shift;
            }
        }
        (best, best_shift)
    }
}

fn invert_columns(columns: &[u64], n: usize) -> Option<Vec<u64>> {
    let mut rows: Vec<u128> = (0..n)
        .map(|r| {
            let mut row = 0u128;
            for (i, &column) in columns.iter().enumerate() {
                if (column >> r) & 1 == 1 {
                    row |= 1u128 << i;
                }
            }
            row | (1u128 << (64 + r))
        })
        .collect();
    for pivot in 0..n {
        let found = (pivot..n).find(|&r| (rows[r] >> pivot) & 1 == 1)?;
        rows.swap(pivot, found);
        for r in 0..n {
            if r != pivot && (rows[r] >> pivot) & 1 == 1 {
                rows[r] ^= rows[pivot];
            }
        }
    }
    let mut images = vec![0u64; n];
    for (i, row) in rows.iter().enumerate() {
        let high = (row >> 64) as u64;
        for (r, image) in images.iter_mut().enumerate() {
            if (high >> r) & 1 == 1 {
                *image |= 1u64 << i;
            }
        }
    }
    Some(images)
}

fn mulmod(a: u64, b: u64, r: u64) -> u64 {
    ((a as u128 * b as u128) % r as u128) as u64
}

fn powmod(mut a: u64, mut e: u64, r: u64) -> u64 {
    let mut acc = 1u64;
    while e > 0 {
        if e & 1 == 1 {
            acc = mulmod(acc, a, r);
        }
        a = mulmod(a, a, r);
        e >>= 1;
    }
    acc
}

#[inline(always)]
fn mix(mut value: u64) -> u64 {
    value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    value ^ (value >> 31)
}

/// Deterministic scalar stream from a label.
struct Stream {
    state: [u8; 32],
    counter: u64,
}

impl Stream {
    fn new(label: &str) -> Self {
        Self {
            state: *blake3::hash(label.as_bytes()).as_bytes(),
            counter: 0,
        }
    }

    fn next_u128(&mut self) -> u128 {
        let mut hasher = blake3::Hasher::new();
        hasher.update(&self.state);
        hasher.update(&self.counter.to_le_bytes());
        self.counter += 1;
        let digest = hasher.finalize();
        u128::from_le_bytes(digest.as_bytes()[..16].try_into().unwrap())
    }

    fn scalar(&mut self, r: u64) -> u64 {
        loop {
            let v = (self.next_u128() % r as u128) as u64;
            if v != 0 {
                return v;
            }
        }
    }
}

struct Setup {
    n: u32,
    a: u8,
    r: u64,
    lambda: u64,
    curve: KoblitzCurve,
    fast: FastBinaryCurve,
    gf: Gf2,
    b: u64,
    basis: NormalBasis,
    generator: FastPoint,
}

fn setup(n: u32, a: u8) -> Setup {
    let curve = KoblitzCurve::new(a, n).expect("admitted Koblitz rung");
    let r = curve.subgroup_order.to_u64().expect("r < 2^64");
    let lambda = curve.lambda.to_u64().unwrap();
    let fast = FastBinaryCurve::new(&curve.curve.irreducible, a as u64).expect("n <= 63");
    let gf = Gf2::new(&curve.curve.irreducible);
    let b = gf.from_element(&curve.curve.b);
    let basis = NormalBasis::new(&gf);
    let generator: FastPoint = match curve.generator() {
        crypto_lib::binary_ecc::BinaryPoint::Affine { x, y } => Some((fast.word(x), fast.word(y))),
        crypto_lib::binary_ecc::BinaryPoint::Infinity => panic!("generator must be affine"),
    };
    // Frobenius acts as lambda on the subgroup and as a left rotation by one
    // in the normal basis.
    let (gx, gy) = generator.unwrap();
    assert_eq!(
        fast.scalar_mul(generator, &BigUint::from(lambda)),
        Some((gf.sqr(gx), gf.sqr(gy)))
    );
    for probe in 1..512u64 {
        let x = probe.wrapping_mul(0x2545_F491_4F6C_DD1D) & basis.mask;
        assert_eq!(
            basis.to_poly.apply(basis.rotate(basis.to_normal.apply(x), 1)),
            gf.sqr(x)
        );
    }
    Setup {
        n,
        a,
        r,
        lambda,
        curve,
        fast,
        gf,
        b,
        basis,
        generator,
    }
}

fn peak_rss_bytes() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            let kb: u64 = rest.trim().trim_end_matches("kB").trim().parse().ok()?;
            return Some(kb * 1024);
        }
    }
    None
}

// ---------------------------------------------------------------- targets mode

fn mode_targets(args: &[String]) {
    let n: u32 = args[2].parse().unwrap();
    let a: u8 = args[3].parse().unwrap();
    let count: usize = args[4].parse().unwrap();
    let label = &args[5];
    let s = setup(n, a);
    let mut stream = Stream::new(&format!("ic-next10-09-targets|{n}|{a}|{label}"));
    let mut points = std::fs::File::create(&args[6]).unwrap();
    let mut scalars = std::fs::File::create(&args[7]).unwrap();
    for _ in 0..count {
        let d = stream.scalar(s.r);
        let (x, y) = s.fast.scalar_mul(s.generator, &BigUint::from(d)).unwrap();
        writeln!(points, "[{x},{y}]").unwrap();
        writeln!(scalars, "{d}").unwrap();
    }
}

fn read_points(path: &str) -> Vec<(u64, u64)> {
    std::fs::read_to_string(path)
        .expect("read targets")
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let [x, y]: [u64; 2] = serde_json::from_str(l.trim()).expect("[x,y]");
            (x, y)
        })
        .collect()
}

// ---------------------------------------------------------------- BL rho mode

/// Walk state: canonical point `P = [u]G + [v]Q` (online) or `[u]G` (offline).
#[derive(Clone, Copy)]
struct Lane {
    x: u64,
    y: u64,
    /// Normal-basis canonical x (class key).
    key: u64,
    u: u64,
    v: u64,
    len: u64,
}

const JUMPS: usize = 2048;
/// Look-ahead retries against fruitless 2-cycles.
const LOOKAHEAD: usize = 4;

struct Walk<'a> {
    s: &'a Setup,
    jump_points: Vec<(u64, u64)>,
    jump_logs: Vec<u64>,
    lambda_pow: Vec<u64>,
    dp_mask: u64,
    map_salt: u64,
}

impl<'a> Walk<'a> {
    fn new(s: &'a Setup, map_seed: &str, dp_bits: u32) -> Self {
        let mut stream = Stream::new(&format!("ic-next10-09-jumps|{}|{}|{map_seed}", s.n, s.a));
        let mut jump_points = Vec::with_capacity(JUMPS);
        let mut jump_logs = Vec::with_capacity(JUMPS);
        for _ in 0..JUMPS {
            let m = stream.scalar(s.r);
            jump_points.push(s.fast.scalar_mul(s.generator, &BigUint::from(m)).unwrap());
            jump_logs.push(m);
        }
        let mut lambda_pow = Vec::with_capacity(s.n as usize);
        let mut acc = 1u64;
        for _ in 0..s.n {
            lambda_pow.push(acc);
            acc = mulmod(acc, s.lambda, s.r);
        }
        Self {
            s,
            jump_points,
            jump_logs,
            lambda_pow,
            dp_mask: (1u64 << dp_bits) - 1,
            map_salt: mix(stream.next_u128() as u64),
        }
    }

    #[inline(always)]
    fn partition(&self, key: u64) -> usize {
        (mix(key ^ self.map_salt) >> 40) as usize % JUMPS
    }

    #[inline(always)]
    fn is_dp(&self, key: u64) -> bool {
        (mix(key.wrapping_add(0x5851_F42D_4C95_7F2D)) & self.dp_mask) == 0
    }

    /// Canonical class representative of `(x, y)` with coefficients `(u, v)`.
    #[inline(always)]
    fn canon(&self, x: u64, y: u64, u: u64, v: u64) -> Lane {
        let b = &self.s.basis;
        let r = self.s.r;
        let xn = b.to_normal.apply(x);
        let (key, k) = b.canonical(xn);
        let (cx, mut cy) = if k == 0 {
            (x, y)
        } else {
            (
                b.to_poly.apply(key),
                b.to_poly.apply(b.rotate(b.to_normal.apply(y), k)),
            )
        };
        let mut m = self.lambda_pow[k as usize];
        let other = cx ^ cy;
        if other < cy {
            cy = other;
            m = (r - m) % r;
        }
        Lane {
            x: cx,
            y: cy,
            key,
            u: mulmod(u, m, r),
            v: mulmod(v, m, r),
            len: 0,
        }
    }

    /// One iteration with look-ahead, scalar path (one inversion).
    fn step_single(&self, lane: &Lane) -> Option<Lane> {
        let mut j = self.partition(lane.key);
        let mut candidate = None;
        for attempt in 0..LOOKAHEAD {
            let (jx, jy) = self.jump_points[j];
            let sum = self.s.fast.add(Some((lane.x, lane.y)), Some((jx, jy)))?;
            let next = self.canon(sum.0, sum.1, (lane.u + self.jump_logs[j]) % self.s.r, lane.v);
            let nj = self.partition(next.key);
            candidate = Some(next);
            if nj != j || attempt + 1 == LOOKAHEAD {
                break;
            }
            j = (j + 1) % JUMPS;
        }
        candidate.map(|mut c| {
            c.len = lane.len + 1;
            c
        })
    }
}

/// Batched lanes: one shared inversion per round; look-ahead retries and the
/// rare exceptional sums fall back to `step_single`, so every lane follows
/// exactly the map of `step_single`.
fn step_batch(walk: &Walk, lanes: &mut [Lane], dxs: &mut Vec<u64>, scratch: &mut Vec<u64>) {
    let gf = &walk.s.gf;
    let a = walk.s.fast.a;
    dxs.clear();
    for lane in lanes.iter() {
        let j = walk.partition(lane.key);
        let (jx, _) = walk.jump_points[j];
        let d = lane.x ^ jx;
        dxs.push(if d == 0 { 1 } else { d });
    }
    gf.batch_inv(dxs, scratch);
    for (lane, &inv) in lanes.iter_mut().zip(dxs.iter()) {
        let j = walk.partition(lane.key);
        let (jx, jy) = walk.jump_points[j];
        if lane.x == jx {
            *lane = walk.step_single(lane).unwrap_or(Lane { len: u64::MAX, ..*lane });
            continue;
        }
        let lambda = gf.mul(lane.y ^ jy, inv);
        let x3 = gf.sqr(lambda) ^ lambda ^ lane.x ^ jx ^ a;
        let y3 = gf.mul(lambda, lane.x ^ x3) ^ x3 ^ lane.y;
        let next = walk.canon(x3, y3, (lane.u + walk.jump_logs[j]) % walk.s.r, lane.v);
        if walk.partition(next.key) == j {
            // Look-ahead branch: recompute this step on the scalar path.
            *lane = walk.step_single(lane).unwrap_or(Lane { len: u64::MAX, ..*lane });
            continue;
        }
        let len = lane.len + 1;
        *lane = Lane { len, ..next };
    }
}

fn mode_bl(args: &[String]) {
    let n: u32 = args[2].parse().unwrap();
    let a: u8 = args[3].parse().unwrap();
    let table_bytes: u64 = args[4].parse().unwrap();
    let targets = read_points(&args[5]);
    let mut out = std::fs::File::create(&args[6]).unwrap();
    let seed = std::env::var("BL_SEED").unwrap_or_else(|_| "v1".to_owned());
    let online_map_seed = std::env::var("BL_ONLINE_MAP_SEED").unwrap_or_else(|_| seed.clone());
    let prefixes: u32 = std::env::var("BL_PREFIXES")
        .ok()
        .map(|v| v.parse().unwrap())
        .unwrap_or(3);
    let precompute_lanes: usize = 256;

    let process_started = Instant::now();
    let s = setup(n, a);
    let r = s.r;
    let classes = r as f64 / (2.0 * n as f64);
    let entries = (table_bytes / 16) as usize;
    let dp_bits = ((classes / entries as f64).sqrt().log2()).round().max(0.0) as u32;
    let walk_len = 1u64 << dp_bits;
    let cap = 20 * walk_len;
    let walk = Walk::new(&s, &seed, dp_bits);
    let setup_ms = process_started.elapsed().as_secs_f64() * 1000.0;

    // ---- precomputation: G-only walks, first T distinct DPs in order.
    let pre_started = Instant::now();
    let mut stream = Stream::new(&format!("ic-next10-09-pre|{n}|{a}|{seed}"));
    let fresh = |stream: &mut Stream| -> Lane {
        let u = stream.scalar(r);
        let (x, y) = s.fast.scalar_mul(s.generator, &BigUint::from(u)).unwrap();
        walk.canon(x, y, u, 0)
    };
    let mut lanes: Vec<Lane> = (0..precompute_lanes).map(|_| fresh(&mut stream)).collect();
    // The table is kept sorted by (key, epoch<<48 | u); u < r < 2^48. A DP
    // found again later is dropped (the earliest epoch is kept), so prefix
    // tables are "the first T/4^j distinct DPs" up to flush granularity.
    assert!(r < (1u64 << 48), "epoch packing needs r < 2^48");
    let u_mask = (1u64 << 48) - 1;
    let mut table: Vec<(u64, u64)> = Vec::with_capacity(entries + entries / 32 + 1024);
    let mut pre_steps = 0u64;
    let mut pre_walks = 0u64;
    let mut pre_abandoned = 0u64;
    let mut dxs = Vec::with_capacity(precompute_lanes);
    let mut scratch = Vec::with_capacity(precompute_lanes);
    let mut duplicates = 0u64;
    let mut epoch = 0u64;
    let flush_every = (entries / 64).max(1024);
    let mut pending = 0usize;
    while table.len() < entries {
        step_batch(&walk, &mut lanes, &mut dxs, &mut scratch);
        pre_steps += lanes.len() as u64;
        for lane in lanes.iter_mut() {
            if lane.len == u64::MAX || lane.len > cap {
                pre_abandoned += 1;
                *lane = fresh(&mut stream);
                continue;
            }
            if walk.is_dp(lane.key) {
                table.push((lane.key, (epoch << 48) | lane.u));
                pending += 1;
                pre_walks += 1;
                *lane = fresh(&mut stream);
            }
        }
        if pending >= flush_every || table.len() >= entries {
            table.sort_unstable();
            let before = table.len();
            table.dedup_by(|later, earlier| later.0 == earlier.0);
            duplicates += (before - table.len()) as u64;
            pending = 0;
            epoch += 1;
        }
    }
    // Trim to exactly `entries`, dropping the newest epoch's surplus.
    if table.len() > entries {
        table.sort_unstable_by_key(|e| (e.1 >> 48, e.0));
        table.truncate(entries);
        table.sort_unstable();
    }
    let pre_ms = pre_started.elapsed().as_secs_f64() * 1000.0;
    // Spot-verify table logs.
    for &(key, packed) in table.iter().step_by((entries / 64).max(1)) {
        let p = s.fast.scalar_mul(s.generator, &BigUint::from(packed & u_mask)).unwrap();
        let c = walk.canon(p.0, p.1, 1, 0);
        assert_eq!(c.key, key, "table entry log does not reproduce its class");
    }
    let summary_pre = json!({
        "kind":"bl_precompute", "n":n, "a":a, "r":r, "classes_N":classes,
        "table_bytes_budget":table_bytes, "entries":entries, "bytes_per_entry":16,
        "dp_bits":dp_bits, "walk_len_2t":walk_len, "walk_cap":cap,
        "precompute_steps":pre_steps, "precompute_walks":pre_walks,
        "precompute_abandoned":pre_abandoned, "precompute_duplicate_dps":duplicates,
        "precompute_ms":pre_ms, "precompute_lanes":precompute_lanes,
        "ns_per_precompute_step":pre_ms * 1e6 / pre_steps as f64,
        "setup_ms":setup_ms, "map_seed":seed, "online_map_seed":online_map_seed,
        "jumps":JUMPS, "lookahead":LOOKAHEAD,
    });
    writeln!(out, "{summary_pre}").unwrap();
    println!("{summary_pre}");

    // ---- online: per prefix table size, every target, single walk.
    let online_walk = Walk::new(&s, &online_map_seed, dp_bits);
    for level in 0..prefixes {
        let size = entries >> (2 * level);
        if size == 0 {
            break;
        }
        let prefix: Option<Vec<(u64, u64)>> = (level > 0).then(|| {
            // Earliest epochs first, without copying the full table.
            let mut histogram = vec![0usize; (epoch + 2) as usize];
            for e in table.iter() {
                histogram[(e.1 >> 48) as usize] += 1;
            }
            let mut cumulative = 0usize;
            let mut boundary = 0u64;
            for (index, &count) in histogram.iter().enumerate() {
                if cumulative + count >= size {
                    boundary = index as u64;
                    break;
                }
                cumulative += count;
            }
            let mut partial = size - cumulative;
            let mut chosen = Vec::with_capacity(size);
            for &e in table.iter() {
                let ep = e.1 >> 48;
                if ep < boundary {
                    chosen.push(e);
                } else if ep == boundary && partial > 0 {
                    chosen.push(e);
                    partial -= 1;
                }
            }
            chosen
        });
        let table: &[(u64, u64)] = prefix.as_deref().unwrap_or(&table);
        let mut steps_all = Vec::new();
        let mut ms_all = Vec::new();
        for (index, &(qx, qy)) in targets.iter().enumerate() {
            let mut tstream = Stream::new(&format!("ic-next10-09-online|{n}|{a}|{seed}|{index}"));
            let started = Instant::now();
            let mut steps = 0u64;
            let mut walks = 0u64;
            let mut misses = 0u64;
            let mut abandoned = 0u64;
            let mut result = None;
            let budget = 400 * walk_len.max(1);
            'outer: while steps < budget {
                walks += 1;
                let c = tstream.scalar(r);
                let cg = s.fast.scalar_mul(s.generator, &BigUint::from(c));
                let Some((sx, sy)) = s.fast.add(cg, Some((qx, qy))) else {
                    continue;
                };
                let mut lane = online_walk.canon(sx, sy, c, 1);
                loop {
                    if online_walk.is_dp(lane.key) {
                        if let Ok(pos) = table.binary_search_by_key(&lane.key, |e| e.0) {
                            let ut = table[pos].1 & u_mask;
                            // [u]G + [v]Q = [ut]G  =>  d = (ut - u) / v.
                            if lane.v != 0 {
                                let d = mulmod((ut + r - lane.u) % r, powmod(lane.v, r - 2, r), r);
                                if s.fast.scalar_mul(s.generator, &BigUint::from(d)) == Some((qx, qy)) {
                                    result = Some(d);
                                    break 'outer;
                                }
                            }
                        }
                        misses += 1;
                        break;
                    }
                    if lane.len > cap || steps >= budget {
                        abandoned += 1;
                        break;
                    }
                    match online_walk.step_single(&lane) {
                        Some(next) => lane = next,
                        None => {
                            abandoned += 1;
                            break;
                        }
                    }
                    steps += 1;
                }
            }
            let elapsed = started.elapsed().as_secs_f64() * 1000.0;
            steps_all.push(steps);
            ms_all.push(elapsed);
            let record = json!({
                "kind":"bl_online_target", "n":n, "a":a, "target_index":index,
                "table_entries":size, "table_bytes":size * 16, "prefix_level":level,
                "dp_bits":dp_bits, "online_steps":steps, "walks":walks,
                "dp_misses":misses, "abandoned":abandoned,
                "solved":result.is_some(), "recovered_scalar":result,
                "online_ms":elapsed, "online_map_seed":online_map_seed,
                "budget_steps":budget,
            });
            writeln!(out, "{record}").unwrap();
        }
        let solved = steps_all.len();
        let mean_steps = steps_all.iter().sum::<u64>() as f64 / solved.max(1) as f64;
        let mean_ms = ms_all.iter().sum::<f64>() / solved.max(1) as f64;
        println!(
            "{}",
            json!({"kind":"bl_online_summary","prefix_level":level,"table_entries":size,
                   "mean_steps":mean_steps,"mean_ms":mean_ms,"targets":solved})
        );
    }
    let summary = json!({
        "kind":"bl_summary", "process_ms":process_started.elapsed().as_secs_f64() * 1000.0,
        "peak_rss_bytes":peak_rss_bytes(),
    });
    writeln!(out, "{summary}").unwrap();
    println!("{summary}");
}

// ---------------------------------------------------------------- census mode

/// Copied from `koblitz_orbit_dlp_fast_online.rs::construct_base`.
fn construct_base(
    s: &Setup,
    columns: usize,
) -> (Vec<(u64, u64)>, Vec<(usize, u64)>, Vec<[u64; 2]>) {
    let fast = &s.fast;
    let gf = &s.gf;
    let n = s.n as usize;
    let r = s.r;
    let mut seen = HashSet::new();
    let mut points = Vec::with_capacity(columns * 2 * n);
    let mut labels = Vec::with_capacity(columns * 2 * n);
    let mut reps = Vec::with_capacity(columns);
    let mut raw_x = 1u64;
    while reps.len() < columns {
        for point in fast.points_with_x(s.b, raw_x) {
            let Some((px, py)) = fast.scalar_mul(point, &s.curve.cofactor) else {
                continue;
            };
            if px <= 1 {
                continue;
            }
            let mut key = px;
            let mut x = px;
            for _ in 1..n {
                x = gf.sqr(x);
                key = key.min(x);
            }
            if !seen.insert(key) {
                continue;
            }
            let column = reps.len();
            reps.push([px, py]);
            let (mut cx, mut cy) = (px, py);
            let mut coefficient = 1u64;
            for _ in 0..n {
                points.push((cx, cy));
                labels.push((column, coefficient));
                points.push((cx, cx ^ cy));
                labels.push((column, (r - coefficient) % r));
                cx = gf.sqr(cx);
                cy = gf.sqr(cy);
                coefficient = mulmod(coefficient, s.lambda, r);
            }
            if reps.len() == columns {
                break;
            }
        }
        raw_x += 1;
    }
    (points, labels, reps)
}

fn mode_census(args: &[String]) {
    let n: u32 = args[2].parse().unwrap();
    let a: u8 = args[3].parse().unwrap();
    let k: usize = args[4].parse().unwrap();
    let mut out = std::fs::File::create(&args[5]).unwrap();
    let started = Instant::now();
    let s = setup(n, a);
    let r = s.r;
    let nn = n as usize;
    let (points, labels, reps) = construct_base(&s, k);
    let mut by_x: std::collections::HashMap<u64, Vec<usize>> = std::collections::HashMap::new();
    for (i, p) in points.iter().enumerate() {
        by_x.entry(p.0).or_default().push(i);
    }
    let shifted: Vec<Vec<u64>> = reps
        .iter()
        .map(|rep| {
            let mut row = Vec::with_capacity(nn);
            let mut v = rep[0];
            for _ in 0..nn {
                row.push(v);
                v = s.gf.sqr(v);
            }
            row
        })
        .collect();
    let base_ms = started.elapsed().as_secs_f64() * 1000.0;

    // Enumerate unordered pair-state orbits: l <= r; for l == r skip the mate
    // relative (n - rel) and rel == 0 (P + P / P - P are not pair states of
    // distinct points). Payload: l (16) | r (16) | rel (16) | root bit (1).
    let enum_started = Instant::now();
    let mut keys: Vec<(u64, u64)> = Vec::with_capacity(k * (k + 1) * nn + 16);
    let mut pair_candidates = 0u64;
    let mut root_states = 0u64;
    for l in 0..k {
        for rr in l..k {
            for rel in 0..nn {
                if l == rr && (rel == 0 || rel > nn - rel) {
                    continue;
                }
                pair_candidates += 1;
                let left = shifted[l][0];
                let right = shifted[rr][rel];
                if let Some(roots) = s3_x_roots(&s.gf, s.b, left, right) {
                    root_states += 1;
                    let distinct = if roots[0] == roots[1] { 1 } else { 2 };
                    for (bit, &root) in roots.iter().take(distinct).enumerate() {
                        let (canon, _) = s.basis.canonical(s.basis.to_normal.apply(root));
                        let payload =
                            (l as u64) | ((rr as u64) << 16) | ((rel as u64) << 32) | ((bit as u64) << 48);
                        keys.push((canon, payload));
                    }
                }
            }
        }
    }
    let enum_ms = enum_started.elapsed().as_secs_f64() * 1000.0;
    let sort_started = Instant::now();
    keys.sort_unstable();
    let sort_ms = sort_started.elapsed().as_secs_f64() * 1000.0;

    // Collisions: every pair inside an equal-key run.
    let verify_started = Instant::now();
    let mut collisions_pairs = 0u64;
    let mut verified = 0u64;
    let mut trivial = 0u64;
    let mut failed = 0u64;
    let mut max_run = 1usize;
    let mut rows: Vec<Vec<u64>> = Vec::new();
    let mut records: Vec<Value> = Vec::new();
    let unpack = |p: u64| -> (usize, usize, usize, usize) {
        (
            (p & 0xffff) as usize,
            ((p >> 16) & 0xffff) as usize,
            ((p >> 32) & 0xffff) as usize,
            ((p >> 48) & 1) as usize,
        )
    };
    let mut i = 0;
    while i < keys.len() {
        let mut j = i + 1;
        while j < keys.len() && keys[j].0 == keys[i].0 {
            j += 1;
        }
        max_run = max_run.max(j - i);
        for p in i..j {
            for q in p + 1..j {
                collisions_pairs += 1;
                let (l1, r1, rel1, _) = unpack(keys[p].1);
                let (l2, r2, rel2, _) = unpack(keys[q].1);
                // Find a signed four-point identity among the x-codes, with a
                // Frobenius rotation on the second pair: rotate pair 2 by every
                // shift t (the canonical key equality fixes one t, try all).
                let x1 = shifted[l1][0];
                let x2 = shifted[r1][rel1];
                // Prefer a non-trivial identity: a collision is trivial only
                // when every signed four-point identity it carries has a zero
                // log row (the same pair orbit seen twice).
                let mut found: Option<([usize; 4], Vec<u64>)> = None;
                let mut any_identity = false;
                'shift: for t in 0..nn {
                    let x3 = shifted[l2][t % nn];
                    let x4 = shifted[r2][(rel2 + t) % nn];
                    let (Some(c1), Some(c2), Some(c3), Some(c4)) =
                        (by_x.get(&x1), by_x.get(&x2), by_x.get(&x3), by_x.get(&x4))
                    else {
                        continue;
                    };
                    for &p1 in c1 {
                        for &p2 in c2 {
                            let s12 = s.fast.add(Some(points[p1]), Some(points[p2]));
                            for &p3 in c3 {
                                let s123 = s.fast.add(s12, Some(points[p3]));
                                for &p4 in c4 {
                                    if s.fast.add(s123, Some(points[p4])).is_none() {
                                        any_identity = true;
                                        let indices = [p1, p2, p3, p4];
                                        let mut row = vec![0u64; k];
                                        for &idx in &indices {
                                            let (column, coefficient) = labels[idx];
                                            row[column] = (row[column] + coefficient) % r;
                                        }
                                        if row.iter().any(|&v| v != 0) {
                                            found = Some((indices, row));
                                            break 'shift;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                let Some((indices, row)) = found else {
                    if any_identity {
                        trivial += 1;
                    } else {
                        failed += 1;
                    }
                    continue;
                };
                verified += 1;
                if records.len() < 2000 {
                    records.push(json!({
                        "kind":"free_relation", "states":[[l1,r1,rel1],[l2,r2,rel2]],
                        "point_indices":indices,
                        "points":indices.iter().map(|&x| [points[x].0, points[x].1]).collect::<Vec<_>>(),
                        "support":row.iter().enumerate().filter(|(_, v)| **v != 0).map(|(c, v)| [c as u64, *v]).collect::<Vec<_>>(),
                    }));
                }
                rows.push(row);
            }
        }
        i = j;
    }
    let verify_ms = verify_started.elapsed().as_secs_f64() * 1000.0;

    // Rank mod r of the homogeneous rows (columns = orbit logs).
    let rank_started = Instant::now();
    let mut pivots: Vec<Option<Vec<u64>>> = vec![None; k];
    let mut rank = 0usize;
    for row in &rows {
        let mut row = row.clone();
        for column in 0..k {
            if row[column] == 0 {
                continue;
            }
            if let Some(pivot) = &pivots[column] {
                let factor = row[column];
                for (value, &pv) in row.iter_mut().zip(pivot.iter()).skip(column) {
                    *value = (*value + r - mulmod(factor, pv, r)) % r;
                }
            } else {
                let inverse = powmod(row[column], r - 2, r);
                for value in row.iter_mut().skip(column) {
                    *value = mulmod(*value, inverse, r);
                }
                pivots[column] = Some(row);
                rank += 1;
                break;
            }
        }
    }
    let rank_ms = rank_started.elapsed().as_secs_f64() * 1000.0;
    let columns_touched = {
        let mut seen = vec![false; k];
        for row in &rows {
            for (c, v) in row.iter().enumerate() {
                if *v != 0 {
                    seen[c] = true;
                }
            }
        }
        seen.iter().filter(|&&x| x).count()
    };
    for rec in &records {
        writeln!(out, "{rec}").unwrap();
    }
    let m = keys.len() as f64;
    let classes = r as f64 / (2.0 * n as f64);
    let summary = json!({
        "kind":"free_relation_census", "n":n, "a":a, "r":r, "K":k,
        "factor_base_points":points.len(), "pair_state_orbits":pair_candidates,
        "root_states":root_states, "root_keys":keys.len(),
        "classes_r_over_2n":classes,
        "expected_pairs_uniform":m * (m - 1.0) / 2.0 / classes,
        "e3_prediction_K4n3_over_r":(k as f64).powi(4) * (n as f64).powi(3) / r as f64,
        "leads_prediction_n4K4_over_r":(k as f64).powi(4) * (n as f64).powi(4) / r as f64,
        "collision_pairs":collisions_pairs, "verified_nontrivial":verified,
        "trivial_rows":trivial, "verification_failed":failed, "max_run":max_run,
        "rank_mod_r":rank, "columns_touched":columns_touched,
        "timing_ms":{"base":base_ms,"enumerate":enum_ms,"sort":sort_ms,"verify":verify_ms,"rank":rank_ms,
                     "total":started.elapsed().as_secs_f64()*1000.0},
        "peak_rss_bytes":peak_rss_bytes(),
    });
    writeln!(out, "{summary}").unwrap();
    println!("{summary}");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("targets") => mode_targets(&args),
        Some("bl") => mode_bl(&args),
        Some("census") => mode_census(&args),
        _ => panic!("usage: targets|bl|census ..."),
    }
}
