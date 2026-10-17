//! Precomputed distinguished-point Pollard rho (Bernstein–Lange) on the exact
//! legacy n=83/a=1 Koblitz curve: the matched-precomputation control for the
//! N83 compact-orbit index calculus.
//!
//! Usage: n83_precomputed_rho <config.json> <targets.jsonl> <records.jsonl>
//!
//! The targets file holds one JSON array of two decimal coordinates per line
//! (numbers or strings). No scalar or answer sidecar is read.
//!
//! Setup (target independent, timed and reported, outside the online
//! interval exactly like the IC's base/index/rank setup): curve, normal basis,
//! the `G`-only jump table, then `G`-only distinguished-point walks whose
//! endpoints `(fingerprint, a)` fill a table, sorted and deduplicated. The
//! table is bounded either by wall time (`wall_seconds`) or by its byte size
//! (`table_bytes`, 16 bytes per entry).
//!
//! Online, per target (one clock, start before the first operation on the
//! target point, stop after the `[d]G = Q` recovery check): target
//! on-curve and subgroup checks, `Q`-dependent walk starts `c·G + Q`, walks
//! on the same step function to a distinguished point, table lookups (plus
//! a per-target table of the target's own distinguished points, so an
//! ordinary rho collision also counts), and the recovery check. A separate
//! generic library replay of every recovered scalar runs after the clock
//! stops and is recorded as verification.
//!
//! The step is the repository's fastest strong-rho step: the fold field
//! backend of `examples/n83_rho_optimization_eval.rs` (`Gf2_128`), the
//! production jumps, partition, distinguished-point hash and signed-Frobenius
//! canonicalization of `koblitz_strong_rho`, and batched slope inversion.
#![recursion_limit = "256"]

use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::koblitz_strong_rho::{
    dp_fingerprint, DpWalkStats, JumpG, RawPointG, RhoField, StrongRhoCharges, StrongRhoG,
    WalkStateG,
};
use crypto_lib::cryptanalysis::semaev_decomp::Gf2_128;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde::Deserialize;
use serde_json::{json, Value};

const LEGACY_R: u64 = 8_569_786_107_849_059;
const LEGACY_G: [u128; 2] = [
    6_614_244_651_398_184_860_292_015,
    8_868_260_410_029_249_391_539_529,
];
const ENTRY_BYTES: u64 = 16;

/// The fold backend of `n83_rho_optimization_eval.rs`, unchanged.
struct FoldField(Gf2_128);

impl RhoField for FoldField {
    type E = u128;
    fn degree(&self) -> u32 {
        self.0.n
    }
    #[inline(always)]
    fn mul(&self, a: u128, b: u128) -> u128 {
        self.0.mul(a, b)
    }
    #[inline(always)]
    fn sqr(&self, a: u128) -> u128 {
        self.0.sqr(a)
    }
    fn batch_inv(&self, xs: &mut [u128], scratch: &mut Vec<u128>) {
        self.0.batch_inv(xs, scratch);
    }
}

type Rho = StrongRhoG<FoldField, u64>;

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Config {
    arm: String,
    backend: String,
    dp_bits: u32,
    walk_cap_factor: u64,
    precompute_lanes: usize,
    precompute_threads: usize,
    /// Stop the table walks at this wall time from the start of setup.
    budget_wall_seconds: Option<f64>,
    /// Stop when the table holds this many bytes (16 per entry).
    budget_table_bytes: Option<u64>,
    online_lanes: usize,
    online_step_cap: u64,
    table_seed: u64,
    online_seed: u64,
}

impl Config {
    fn parse(bytes: &[u8]) -> Result<Self, String> {
        let c: Self =
            serde_json::from_slice(bytes).map_err(|e| format!("invalid strict config: {e}"))?;
        let wall_ok = c.budget_wall_seconds.is_none_or(|s| s.is_finite() && s > 0.0);
        let bytes_ok = c.budget_table_bytes.is_none_or(|b| b >= ENTRY_BYTES);
        let some_budget = c.budget_wall_seconds.is_some() || c.budget_table_bytes.is_some();
        if c.backend != "fold"
            || !(1..=20).contains(&c.dp_bits)
            || !(2..=64).contains(&c.walk_cap_factor)
            || !(1..=256).contains(&c.precompute_lanes)
            || !(1..=64).contains(&c.precompute_threads)
            || !(1..=256).contains(&c.online_lanes)
            || c.online_step_cap < 1_000
            || !wall_ok
            || !bytes_ok
            || !some_budget
            || c.arm.is_empty()
        {
            return Err("config out of range (backend must be fold, one budget required)".into());
        }
        Ok(c)
    }
}

fn sha256(bytes: &[u8]) -> String {
    hex::encode(crypto_lib::hash::sha256(bytes))
}

fn parse_targets(text: &str) -> Result<Vec<[u128; 2]>, String> {
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let inner = line
            .strip_prefix('[')
            .and_then(|l| l.strip_suffix(']'))
            .ok_or(format!("line {i}: expected a JSON array"))?;
        let parts: Vec<&str> = inner.split(',').map(|p| p.trim().trim_matches('"')).collect();
        if parts.len() != 2 {
            return Err(format!("line {i}: expected two coordinates"));
        }
        let mut point = [0u128; 2];
        for (slot, part) in point.iter_mut().zip(&parts) {
            if part.is_empty() || !part.bytes().all(|b| b.is_ascii_digit()) {
                return Err(format!("line {i}: coordinates must be decimal"));
            }
            *slot = part.parse::<u128>().map_err(|_| format!("line {i}: bad coordinate"))?;
            if *slot >> 83 != 0 {
                return Err(format!("line {i}: coordinate is not a canonical field word"));
            }
        }
        out.push(point);
    }
    if out.is_empty() {
        return Err("no targets".into());
    }
    Ok(out)
}

fn point_strings(point: &BinaryPoint) -> Result<[String; 2], String> {
    match point {
        BinaryPoint::Infinity => Err("affine point required".into()),
        BinaryPoint::Affine { x, y } => {
            Ok([x.to_biguint().to_string(), y.to_biguint().to_string()])
        }
    }
}

fn assert_legacy_curve(curve: &KoblitzCurve) -> Result<(), String> {
    if curve.n != 83
        || curve.a != 1
        || curve.curve.irreducible.low_terms != [0, 2, 4, 7]
        || curve.subgroup_order.to_u64() != Some(LEGACY_R)
        || curve.cofactor.to_u64() != Some(1_128_547_018)
        || point_strings(curve.generator())? != [LEGACY_G[0].to_string(), LEGACY_G[1].to_string()]
    {
        return Err("constructor does not match the exact legacy n83 a1 curve and generator".into());
    }
    Ok(())
}

fn peak_rss_bytes() -> Option<u64> {
    // SAFETY: getrusage fills the initialized plain-data structure for this process.
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    if unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage) } != 0 {
        return None;
    }
    let raw = u64::try_from(usage.ru_maxrss).ok()?;
    #[cfg(target_os = "macos")]
    {
        Some(raw)
    }
    #[cfg(not(target_os = "macos"))]
    {
        Some(raw.saturating_mul(1024))
    }
}

fn cpu_seconds() -> f64 {
    // SAFETY: as above.
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    if unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage) } != 0 {
        return f64::NAN;
    }
    let tv = |t: libc::timeval| t.tv_sec as f64 + t.tv_usec as f64 * 1e-6;
    tv(usage.ru_utime) + tv(usage.ru_stime)
}

fn derive_seed(domain: &[u8], seed: u64, extra: &[u8]) -> u64 {
    let mut material = domain.to_vec();
    material.extend_from_slice(&seed.to_le_bytes());
    material.extend_from_slice(extra);
    u64::from_le_bytes(blake3::hash(&material).as_bytes()[..8].try_into().unwrap())
}

/// `G`-only walks from `c·G` (`b = 0`), `c` stepped by a random stride,
/// writing `(fingerprint, a)` into `slots` until it is full or the deadline
/// passes. Returns the slots written.
fn precompute_thread(
    rho: &Rho,
    jumps: &[JumpG<u128, u64>],
    config: &Config,
    seed: u64,
    deadline: Option<Instant>,
    slots: &mut [(u64, u64)],
    stats: &mut DpWalkStats,
) -> usize {
    let r = rho.modulus();
    let mut rng = StdRng::seed_from_u64(seed);
    let stride_a = rng.gen_range(1..r);
    let stride = rho.scalar_mul(rho.generator(), stride_a);
    let mut cursor_a = rng.gen_range(0..r);
    let mut cursor = rho.scalar_mul(rho.generator(), cursor_a);
    let walk_cap = config.walk_cap_factor << config.dp_bits;
    let mut starts = 0u64;
    let mut written = 0usize;
    if slots.is_empty() {
        return 0;
    }
    rho.dp_walks(
        jumps,
        config.precompute_lanes,
        config.dp_bits,
        walk_cap,
        || {
            starts += 1;
            if let Some(deadline) = deadline {
                if starts % 64 == 0 && Instant::now() >= deadline {
                    return None;
                }
            }
            let start = WalkStateG {
                point: cursor,
                a: cursor_a,
                b: 0,
            };
            cursor = rho.add(cursor, stride);
            cursor_a = (cursor_a + stride_a) % r;
            Some(start)
        },
        |state, _| {
            slots[written] = (dp_fingerprint(state.point), state.a);
            written += 1;
            written >= slots.len()
        },
        stats,
    );
    written
}

struct Table {
    entries: Vec<(u64, u64)>,
    stats: DpWalkStats,
    raw_entries: u64,
    stopped_by: &'static str,
    walk_ns: u128,
    sort_ns: u128,
}

/// Fill one preallocated slot array (split evenly across threads), then
/// compact, sort and deduplicate it in place: peak table memory is the slot
/// array itself. Under a wall budget alone the array is sized from a
/// conservative step-rate bound and grows only if a thread fills its share.
fn build_table(rho: &Rho, jumps: &[JumpG<u128, u64>], config: &Config, start: Instant) -> Table {
    let threads = config.precompute_threads;
    let deadline = config
        .budget_wall_seconds
        .map(|seconds| start + std::time::Duration::from_secs_f64(seconds));
    let capacity = match config.budget_table_bytes {
        Some(bytes) => (bytes / ENTRY_BYTES) as usize,
        // Wall budget only: 2e7 steps per second per thread bounds any core.
        None => {
            let seconds = config.budget_wall_seconds.unwrap();
            ((seconds * 2e7 * threads as f64) / (1u64 << config.dp_bits) as f64) as usize + 1024
        }
    };
    let per_thread = capacity.div_ceil(threads);
    let walk_started = Instant::now();
    let mut entries: Vec<(u64, u64)> = vec![(0, 0); capacity];
    let mut stats = DpWalkStats::default();
    let mut filled = Vec::with_capacity(threads);
    std::thread::scope(|scope| {
        let handles: Vec<_> = entries
            .chunks_mut(per_thread)
            .enumerate()
            .map(|(t, chunk)| {
                let seed = derive_seed(b"n83-bl-table-walks-v1\0", config.table_seed, &t.to_le_bytes());
                scope.spawn(move || {
                    let mut stats = DpWalkStats::default();
                    let written = precompute_thread(rho, jumps, config, seed, deadline, chunk, &mut stats);
                    (written, stats)
                })
            })
            .collect();
        for handle in handles {
            let (written, s) = handle.join().unwrap();
            stats.steps += s.steps;
            stats.walks += s.walks;
            stats.distinguished += s.distinguished;
            stats.escapes += s.escapes;
            stats.capped += s.capped;
            filled.push(written);
        }
    });
    let walk_ns = walk_started.elapsed().as_nanos();
    let full = filled
        .iter()
        .zip(entries.chunks(per_thread))
        .all(|(&w, chunk)| w == chunk.len());
    let stopped_by = if full { "table_bytes" } else { "wall_seconds" };
    let sort_started = Instant::now();
    let mut end = 0usize;
    for (t, &written) in filled.iter().enumerate() {
        let from = t * per_thread;
        entries.copy_within(from..from + written, end);
        end += written;
    }
    entries.truncate(end);
    let raw_entries = entries.len() as u64;
    entries.sort_unstable_by_key(|e| e.0);
    entries.dedup_by_key(|e| e.0);
    if config.budget_table_bytes.is_none() {
        entries.shrink_to_fit();
    }
    Table {
        entries,
        stats,
        raw_entries,
        stopped_by,
        walk_ns,
        sort_ns: sort_started.elapsed().as_nanos(),
    }
}

fn lookup(table: &[(u64, u64)], key: u64) -> Option<u64> {
    table
        .binary_search_by_key(&key, |e| e.0)
        .ok()
        .map(|i| table[i].1)
}

fn inverse_mod(v: u64, m: u64) -> Option<u64> {
    let (mut t, mut new_t) = (0i128, 1i128);
    let (mut r, mut new_r) = (m as i128, v as i128);
    while new_r != 0 {
        let q = r / new_r;
        (t, new_t) = (new_t, t - q * new_t);
        (r, new_r) = (new_r, r - q * new_r);
    }
    (r == 1).then(|| t.rem_euclid(m as i128) as u64)
}

fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

fn sub_mod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 + m as u128 - b as u128) % m as u128) as u64
}

#[derive(Default)]
struct OnlineCounts {
    table_lookups: u64,
    table_hits: u64,
    own_hits: u64,
    failed_candidates: u64,
    zero_denominators: u64,
}

struct Online {
    scalar: Option<u64>,
    source: &'static str,
    stats: DpWalkStats,
    counts: OnlineCounts,
    on_curve: bool,
    in_subgroup: bool,
    validation_ns: u128,
    walk_ns: u128,
    online_ns: u128,
}

fn on_curve(rho: &Rho, gf: &Gf2_128, q: RawPointG<u128>) -> bool {
    let _ = rho;
    match q {
        RawPointG::Infinity => false,
        RawPointG::Affine { x, y } => {
            // y^2 + xy = x^3 + x^2 + 1 (a = b = 1).
            let x2 = gf.sqr(x);
            gf.sqr(y) ^ gf.mul(x, y) == gf.mul(x2, x) ^ x2 ^ 1
        }
    }
}

fn solve_online(
    rho: &Rho,
    gf: &Gf2_128,
    jumps: &[JumpG<u128, u64>],
    table: &[(u64, u64)],
    config: &Config,
    q: RawPointG<u128>,
    seed: u64,
) -> Online {
    let r = rho.modulus();
    let online_started = Instant::now();
    let validation_started = Instant::now();
    let on_curve = on_curve(rho, gf, q);
    let in_subgroup = on_curve && rho.scalar_mul(q, r) == RawPointG::Infinity;
    let validation_ns = validation_started.elapsed().as_nanos();
    let walk_started = Instant::now();
    let mut stats = DpWalkStats::default();
    let mut counts = OnlineCounts::default();
    let mut scalar = None;
    let mut source = "none";
    if in_subgroup {
        let mut rng = StdRng::seed_from_u64(seed);
        let stride_a = rng.gen_range(1..r);
        let stride = rho.scalar_mul(rho.generator(), stride_a);
        let mut cursor_a = rng.gen_range(0..r);
        let mut cursor = rho.add(rho.scalar_mul(rho.generator(), cursor_a), q);
        let walk_cap = config.walk_cap_factor << config.dp_bits;
        let max_walks = (config.online_step_cap >> config.dp_bits).max(1);
        let mut walks = 0u64;
        let mut own: HashMap<u64, (u64, u64)> = HashMap::new();
        rho.dp_walks(
            jumps,
            config.online_lanes,
            config.dp_bits,
            walk_cap,
            || {
                walks += 1;
                if walks > max_walks {
                    return None;
                }
                let start = WalkStateG {
                    point: cursor,
                    a: cursor_a,
                    b: 1,
                };
                cursor = rho.add(cursor, stride);
                cursor_a = (cursor_a + stride_a) % r;
                Some(start)
            },
            |state, _| {
                let key = dp_fingerprint(state.point);
                counts.table_lookups += 1;
                // P = a·G + b·Q meets a'·G (+ b'·Q): d = (a' − a) / (b − b').
                let mut try_pair = |hit_a: u64, hit_b: u64, from_table: bool| -> bool {
                    let denominator = sub_mod(state.b, hit_b, r);
                    let Some(inv) = inverse_mod(denominator, r) else {
                        counts.zero_denominators += 1;
                        return false;
                    };
                    let candidate = mul_mod(sub_mod(hit_a, state.a, r), inv, r);
                    if rho.scalar_mul(rho.generator(), candidate) == q {
                        scalar = Some(candidate);
                        source = if from_table { "precomputed_table" } else { "own_walks" };
                        return true;
                    }
                    counts.failed_candidates += 1;
                    false
                };
                if let Some(hit_a) = lookup(table, key) {
                    counts.table_hits += 1;
                    if try_pair(hit_a, 0, true) {
                        return true;
                    }
                }
                if let Some(&(hit_a, hit_b)) = own.get(&key) {
                    counts.own_hits += 1;
                    if try_pair(hit_a, hit_b, false) {
                        return true;
                    }
                } else {
                    own.insert(key, (state.a, state.b));
                }
                false
            },
            &mut stats,
        );
    }
    let walk_ns = walk_started.elapsed().as_nanos();
    Online {
        scalar,
        source,
        stats,
        counts,
        on_curve,
        in_subgroup,
        validation_ns,
        walk_ns,
        online_ns: online_started.elapsed().as_nanos(),
    }
}

fn stats_json(s: &DpWalkStats) -> Value {
    json!({
        "steps": s.steps,
        "walks": s.walks,
        "distinguished": s.distinguished,
        "escapes": s.escapes,
        "capped": s.capped,
    })
}

fn run(config_path: &str, targets_path: &str, out_path: &str) -> Result<bool, String> {
    let config_bytes = fs::read(config_path).map_err(|e| format!("read config: {e}"))?;
    let config = Config::parse(&config_bytes)?;
    let targets_bytes = fs::read(targets_path).map_err(|e| format!("read targets: {e}"))?;
    let targets = parse_targets(std::str::from_utf8(&targets_bytes).map_err(|_| "targets not UTF-8")?)?;
    let mut out = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(out_path)
        .map_err(|e| format!("open records: {e}"))?;
    let mut emit = |value: Value| -> Result<(), String> {
        writeln!(out, "{}", serde_json::to_string(&value).unwrap()).map_err(|e| e.to_string())?;
        out.flush().map_err(|e| e.to_string())
    };
    let source = include_bytes!("n83_precomputed_rho.rs");
    let executable = std::env::current_exe().ok().and_then(|p| fs::read(p).ok());
    let started_unix_ms = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| e.to_string())?
        .as_millis();
    #[cfg(target_arch = "x86_64")]
    let clmul = std::arch::is_x86_feature_detected!("pclmulqdq");
    #[cfg(not(target_arch = "x86_64"))]
    let clmul = false;

    // ---- target-independent setup -------------------------------------
    let setup_started = Instant::now();
    let setup_cpu_before = cpu_seconds();
    let curve = KoblitzCurve::new(1, 83).ok_or("legacy curve construction failed")?;
    assert_legacy_curve(&curve)?;
    let gf = Gf2_128::new(&curve.curve.irreducible);
    let rho: Rho = StrongRhoG::from_parts(
        FoldField(Gf2_128::new(&curve.curve.irreducible)),
        1u128,
        LEGACY_R,
        curve.lambda.to_u64().ok_or("lambda does not fit u64")?,
        RawPointG::<u128>::from_binary_wide(curve.generator()),
    );
    let mut charges = StrongRhoCharges::default();
    let jump_seed = derive_seed(b"n83-bl-jumps-v1\0", config.table_seed, &[]);
    let jumps = rho.jumps(jump_seed, &mut charges);
    let prepare_ns = setup_started.elapsed().as_nanos();
    let table = build_table(&rho, &jumps, &config, setup_started);
    let setup_ns = setup_started.elapsed().as_nanos();
    let setup_cpu_s = cpu_seconds() - setup_cpu_before;
    let classes = LEGACY_R as f64 / rho.automorphisms() as f64;
    let w = (1u64 << config.dp_bits) as f64;
    let distinct = table.entries.len() as f64;
    emit(json!({
        "schema": "n83-precomputed-rho/v1",
        "record": "setup",
        "arm": config.arm,
        "started_unix_ms": started_unix_ms,
        "config": serde_json::from_slice::<Value>(&config_bytes).map_err(|e| e.to_string())?,
        "config_sha256": sha256(&config_bytes),
        "targets_sha256": sha256(&targets_bytes),
        "target_count": targets.len(),
        "worker_source_sha256": sha256(source),
        "compiled_source_sha256": {
            "examples/n83_precomputed_rho.rs": sha256(source),
            "src/cryptanalysis/koblitz_strong_rho.rs": sha256(include_bytes!("../src/cryptanalysis/koblitz_strong_rho.rs")),
            "src/cryptanalysis/semaev_decomp.rs": sha256(include_bytes!("../src/cryptanalysis/semaev_decomp.rs")),
        },
        "compiled_source_revision": option_env!("N83_SOURCE_COMMIT"),
        "executable_sha256": executable.as_deref().map(sha256),
        "arithmetic_kernel": {"backend": "fold", "carryless_multiply": if clmul {"pclmulqdq"} else {"portable"}},
        "curve_id": curve.curve_id().icv1,
        "subgroup_order": LEGACY_R.to_string(),
        "automorphisms": rho.automorphisms(),
        "orbit_classes": classes,
        "jump_seed": jump_seed,
        "jump_count": jumps.len(),
        "setup_ns": setup_ns,
        "setup_cpu_seconds": setup_cpu_s,
        "setup_phase_ns": {"curve_basis_jumps": prepare_ns, "table_walks": table.walk_ns, "table_sort_dedup": table.sort_ns},
        "table_walks": stats_json(&table.stats),
        "table_raw_entries": table.raw_entries,
        "table_stopped_by": table.stopped_by,
        "table_slot_bytes_allocated": table.entries.capacity() as u64 * ENTRY_BYTES,
        "table_distinct_entries": table.entries.len(),
        "table_entry_bytes": ENTRY_BYTES,
        "table_bytes": table.entries.len() as u64 * ENTRY_BYTES,
        "table_walk_rate_steps_per_second": table.stats.steps as f64 / (table.walk_ns as f64 * 1e-9),
        "predicted_online_steps": classes / (distinct * w) + w,
        "peak_rss_bytes_after_setup": peak_rss_bytes(),
    }))?;

    // ---- per-target online ---------------------------------------------
    let mut all_ok = true;
    for (index, coords) in targets.iter().enumerate() {
        let q = RawPointG::Affine { x: coords[0], y: coords[1] };
        let seed = derive_seed(
            b"n83-bl-online-v1\0",
            config.online_seed,
            &[coords[0].to_le_bytes(), coords[1].to_le_bytes()].concat(),
        );
        let online = solve_online(&rho, &gf, &jumps, &table.entries, &config, q, seed);
        // Independent generic-library replay, outside the online clock.
        let replay_started = Instant::now();
        let q_lib = BinaryPoint::Affine {
            x: F2mElement::from_biguint(&BigUint::from(coords[0]), 83),
            y: F2mElement::from_biguint(&BigUint::from(coords[1]), 83),
        };
        let replay_ok = online
            .scalar
            .is_some_and(|d| curve.mul(curve.generator(), &BigUint::from(d)) == q_lib);
        let replay_ns = replay_started.elapsed().as_nanos();
        let status = if !online.in_subgroup {
            "REJECTED_PUBLIC_POINT"
        } else if online.scalar.is_none() {
            "UNKNOWN_STEP_CAP"
        } else if !replay_ok {
            "FAILED_SCALAR_REPLAY"
        } else {
            "PASS_VERIFIED_PUBLIC_TARGET"
        };
        all_ok &= status == "PASS_VERIFIED_PUBLIC_TARGET";
        emit(json!({
            "schema": "n83-precomputed-rho/v1",
            "record": "target",
            "arm": config.arm,
            "index": index,
            "status": status,
            "public_target_q": [coords[0].to_string(), coords[1].to_string()],
            "target_on_curve": online.on_curve,
            "target_in_subgroup": online.in_subgroup,
            "recovered_scalar": online.scalar.map(|d| d.to_string()),
            "recovery_source": online.source,
            "generic_scalar_replay_verified": replay_ok,
            "generic_replay_ns_excluded": replay_ns,
            "online_ns": online.online_ns,
            "online_ms": online.online_ns as f64 / 1e6,
            "online_phase_ns": {"target_validation": online.validation_ns, "walks_lookup_and_recovery_check": online.walk_ns},
            "online_interval": {
                "start_event": "before_target_on_curve_and_subgroup_check",
                "stop_event": "after_fast_[d]G=Q_recovery_check",
                "excluded": ["table setup (reported in the setup record)", "generic library replay", "record I/O"],
            },
            "online_walks": stats_json(&online.stats),
            "online_counts": {
                "table_lookups": online.counts.table_lookups,
                "table_hits": online.counts.table_hits,
                "own_walk_hits": online.counts.own_hits,
                "failed_candidates": online.counts.failed_candidates,
                "zero_denominators": online.counts.zero_denominators,
            },
            "walk_seed": seed,
        }))?;
    }
    emit(json!({
        "schema": "n83-precomputed-rho/v1",
        "record": "end",
        "arm": config.arm,
        "all_verified": all_ok,
        "peak_rss_bytes": peak_rss_bytes(),
        "process_cpu_seconds": cpu_seconds(),
    }))?;
    Ok(all_ok)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: n83_precomputed_rho <config.json> <targets.jsonl> <records.jsonl>");
        std::process::exit(2);
    }
    match run(&args[1], &args[2], &args[3]) {
        Ok(true) => {}
        Ok(false) => std::process::exit(1),
        Err(message) => {
            eprintln!("error: {message}");
            std::process::exit(2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn targets_parse_numbers_and_strings_and_refuse_bad_rows() {
        let t = parse_targets("[6274346871052084705058109, 7777137791885914859810330]\n[\"1\",\"2\"]\n").unwrap();
        assert_eq!(t.len(), 2);
        assert_eq!(t[1], [1, 2]);
        assert!(parse_targets("").is_err());
        assert!(parse_targets("[1]").is_err());
        assert!(parse_targets("[1,-2]").is_err());
        assert!(parse_targets(&format!("[{},1]", 1u128 << 83)).is_err());
        assert!(parse_targets("{\"scalar\":5}").is_err());
    }

    #[test]
    fn strict_config_rejects_unknown_fields_and_bad_budgets() {
        let ok = br#"{"arm":"a","backend":"fold","dp_bits":7,"walk_cap_factor":16,"precompute_lanes":64,"precompute_threads":1,"budget_wall_seconds":640.0,"budget_table_bytes":null,"online_lanes":32,"online_step_cap":100000000,"table_seed":1,"online_seed":2}"#;
        assert!(Config::parse(ok).is_ok());
        let bytes = String::from_utf8(ok.to_vec()).unwrap();
        for bad in [
            bytes.replace("fold", "wide"),
            bytes.replace("\"budget_wall_seconds\":640.0", "\"budget_wall_seconds\":-1.0"),
            bytes.replace("\"budget_wall_seconds\":640.0", "\"budget_wall_seconds\":null"),
            bytes.replace("\"budget_table_bytes\":null", "\"budget_table_bytes\":3"),
            bytes.replace("\"online_seed\":2", "\"online_seed\":2,\"scalar\":3"),
        ] {
            assert!(Config::parse(bad.as_bytes()).is_err(), "{bad}");
        }
    }

    #[test]
    fn modular_helpers_match_definitions() {
        let r = LEGACY_R;
        for v in [1u64, 2, 12345, r - 1] {
            assert_eq!(mul_mod(v, inverse_mod(v, r).unwrap(), r), 1);
        }
        assert_eq!(inverse_mod(0, r), None);
        assert_eq!(sub_mod(1, 2, r), r - 1);
    }

    /// Small end-to-end solve on the exact curve: a time-bounded table, then
    /// three known-scalar targets (correctness-only seeds, not a corpus).
    #[test]
    fn exact_curve_small_table_recovers_known_scalars() {
        let curve = KoblitzCurve::new(1, 83).unwrap();
        assert_legacy_curve(&curve).unwrap();
        let gf = Gf2_128::new(&curve.curve.irreducible);
        let rho: Rho = StrongRhoG::from_parts(
            FoldField(Gf2_128::new(&curve.curve.irreducible)),
            1,
            LEGACY_R,
            curve.lambda.to_u64().unwrap(),
            RawPointG::<u128>::from_binary_wide(curve.generator()),
        );
        let mut charges = StrongRhoCharges::default();
        let jumps = rho.jumps(7, &mut charges);
        let config = Config::parse(br#"{"arm":"test","backend":"fold","dp_bits":12,"walk_cap_factor":16,"precompute_lanes":64,"precompute_threads":2,"budget_wall_seconds":20.0,"budget_table_bytes":null,"online_lanes":64,"online_step_cap":400000000,"table_seed":1,"online_seed":2}"#).unwrap();
        let table = build_table(&rho, &jumps, &config, Instant::now());
        assert!(table.entries.windows(2).all(|w| w[0].0 < w[1].0));
        for d in [3u64, 1 << 40, LEGACY_R - 5] {
            let q = rho.scalar_mul(rho.generator(), d);
            let online = solve_online(&rho, &gf, &jumps, &table.entries, &config, q, d);
            assert!(online.in_subgroup);
            assert_eq!(online.scalar, Some(d));
        }
        let bad = RawPointG::Affine { x: 5u128, y: 7u128 };
        let online = solve_online(&rho, &gf, &jumps, &table.entries, &config, bad, 1);
        assert!(!online.on_curve && online.scalar.is_none());
    }
}
