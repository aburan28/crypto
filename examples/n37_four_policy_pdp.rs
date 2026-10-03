//! Complete n37 four-policy, natural-target m3 PDP and rank producer.
//! Protocol: research/notes/ecc2k130/n37_four_policy_pdp_20261003/PROTOCOL.md
//! It reads point-only targets and never opens fixture-scalar files.

#![recursion_limit = "256"]

#[path = "support/n37_policy_common.rs"]
mod common;

use common::{point_from_value, Context, K, N, R};
use crypto_lib::cryptanalysis::binary_velu::{velu_point_map, Curve};
use crypto_lib::cryptanalysis::ec_index_calculus::gaussian_eliminate_mod_n;
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::cryptanalysis::koblitz_relation_solver::{IncrementalRelationSolver, RowStatus};
use crypto_lib::hash::sha256::sha256;
use flate2::read::GzDecoder;
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::fs::OpenOptions;
use std::io::{Read, Write};
use std::time::Instant;

const NOTE: &str = "research/notes/ecc2k130/n37_four_policy_pdp_20261003";
const SUPPORT_GZIP: &str =
    "research/notes/ecc2k130/n37_four_policy_support_20261003/RESULT.json.gz";
const SUPPORT_GZIP_SHA: &str = "8eac2ae4b8d8f0fc452b7cd7c8edbd3183558e55454abbd67349665bd64610d5";
const SUPPORT_RAW_SHA: &str = "05664103a6dab090a8b4298f34c0964e6f253f625cd73d661562d85a2fce87af";
const B03: &str =
    "research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b03.points.jsonl";
const B03_SHA: &str = "84400a2914f06e4a001d0f113f0195952f692634d4ff8285e23ade2599e2bde2";
const B04: &str =
    "research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b04.points.jsonl";
const B04_SHA: &str = "72c5b5361a0b846c9626269ce3928820cb7d420aae479d31dc26f5cbf7e78aa6";
const P: usize = 2 * K * N;
const RAW_TABLE: usize = 1 + P + P * (P + 1) / 2;
const PROBE_CAP: usize = 256;
const PROBE_SEED: u64 = 0x6e33_375f_3462_7064;

#[derive(Clone)]
struct Factor {
    point: FastPoint,
    column: usize,
    coefficient: u64,
}

struct Inputs {
    policies: Vec<Vec<Factor>>,
    source_targets: Vec<FastPoint>,
    leaf_targets: Vec<FastPoint>,
    source_generator: FastPoint,
    leaf_generator: FastPoint,
    validation_scalar_muls: u64,
    validation_map_calls: u64,
}

#[derive(Clone)]
struct Query {
    indices: Option<Vec<usize>>,
    lookups: u32,
    subtractions: u32,
    verification_adds: u32,
}

impl Query {
    fn value(&self) -> Value {
        match &self.indices {
            Some(indices) => json!({
                "status":"hit", "indices":indices, "arity":indices.len(),
                "lookups":self.lookups, "subtractions":self.subtractions,
                "verification_adds":self.verification_adds,
            }),
            None => json!({
                "status":"proved_miss", "indices":Value::Null, "arity":Value::Null,
                "lookups":self.lookups, "subtractions":self.subtractions,
                "verification_adds":self.verification_adds,
            }),
        }
    }
}

struct Table {
    pairs: Vec<(u64, u32)>,
    singles: Vec<(u64, u32)>,
    raw_entries: usize,
    pair_additions: usize,
    retained_bytes_lower_bound: usize,
}

impl Table {
    fn build(curve: &Curve, factors: &[Factor]) -> Result<Self, String> {
        if factors.len() != P {
            return Err("factor array is not 3108 points".into());
        }
        let mut pairs = Vec::new();
        pairs
            .try_reserve_exact(RAW_TABLE)
            .map_err(|error| format!("reserve complete pair table: {error}"))?;
        pairs.push((FastPoint::INFINITY.pack(), 0));
        for (index, factor) in factors.iter().enumerate() {
            pairs.push((factor.point.pack(), (index + 1) as u32));
        }
        for i in 0..P {
            for j in i..P {
                let sum = curve.fast.add(factors[i].point, factors[j].point);
                pairs.push((sum.pack(), (P + 1 + i * P + j) as u32));
            }
        }
        if pairs.len() != RAW_TABLE {
            return Err("complete pair enumeration count changed".into());
        }
        pairs.sort_unstable(); // (full-point key, enumeration code): first witness survives.
        pairs.dedup_by_key(|entry| entry.0);
        let mut singles: Vec<(u64, u32)> = std::iter::once((FastPoint::INFINITY.pack(), 0))
            .chain(
                factors
                    .iter()
                    .enumerate()
                    .map(|(index, factor)| (factor.point.pack(), (index + 1) as u32)),
            )
            .collect();
        singles.sort_unstable();
        singles.dedup_by_key(|entry| entry.0);
        let retained_bytes_lower_bound =
            (pairs.capacity() + singles.capacity()) * std::mem::size_of::<(u64, u32)>();
        Ok(Self {
            pairs,
            singles,
            raw_entries: RAW_TABLE,
            pair_additions: P * (P + 1) / 2,
            retained_bytes_lower_bound,
        })
    }

    fn lookup(&self, key: u64, arity: usize) -> Option<u32> {
        let entries = if arity == 2 {
            &self.singles
        } else {
            &self.pairs
        };
        entries
            .binary_search_by_key(&key, |entry| entry.0)
            .ok()
            .map(|index| entries[index].1)
    }

    fn summary(&self) -> Value {
        json!({
            "raw_candidates":self.raw_entries,
            "distinct_sums":self.pairs.len(),
            "collisions":self.raw_entries-self.pairs.len(),
            "identity_singleton_distinct_sums":self.singles.len(),
            "pair_group_additions":self.pair_additions,
            "retained_bytes_lower_bound":self.retained_bytes_lower_bound,
        })
    }
}

fn decode(code: u32) -> Vec<usize> {
    let code = code as usize;
    if code == 0 {
        Vec::new()
    } else if code <= P {
        vec![code - 1]
    } else {
        let offset = code - P - 1;
        vec![offset / P, offset % P]
    }
}

fn query(
    table: &Table,
    curve: &Curve,
    factors: &[Factor],
    target: FastPoint,
    arity: usize,
) -> Result<Query, String> {
    let mut lookups = 0u32;
    let mut subtractions = 0u32;
    for scan in 0..=P {
        let (residual, extra) = if scan == 0 {
            (target, None)
        } else {
            subtractions += 1;
            let index = scan - 1;
            (
                curve.fast.add(target, curve.fast.neg(factors[index].point)),
                Some(index),
            )
        };
        lookups += 1;
        if let Some(code) = table.lookup(residual.pack(), arity) {
            let mut indices = decode(code);
            if let Some(index) = extra {
                indices.push(index);
            }
            if indices.len() > arity || indices.iter().any(|&index| index >= P) {
                return Err("table returned invalid factor indices".into());
            }
            let mut sum = FastPoint::INFINITY;
            for &index in &indices {
                sum = curve.fast.add(sum, factors[index].point);
            }
            if sum != target {
                return Err("first full-point table hit failed group replay".into());
            }
            return Ok(Query {
                verification_adds: indices.len() as u32,
                indices: Some(indices),
                lookups,
                subtractions,
            });
        }
    }
    Ok(Query {
        indices: None,
        lookups,
        subtractions,
        verification_adds: 0,
    })
}

fn row_from_indices(factors: &[Factor], indices: &[usize]) -> Vec<u64> {
    let mut row = vec![0u64; K];
    for &index in indices {
        let factor = &factors[index];
        row[factor.column] = (row[factor.column] + factor.coefficient) % R;
    }
    row
}

fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut value = *state;
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn probe_scalars() -> Vec<u64> {
    let mut state = PROBE_SEED;
    let mut seen = HashSet::new();
    let mut scalars = Vec::with_capacity(PROBE_CAP);
    while scalars.len() < PROBE_CAP {
        let scalar = 1 + splitmix64(&mut state) % (R - 1);
        if seen.insert(scalar) {
            scalars.push(scalar);
        }
    }
    scalars
}

fn read_support() -> Result<Value, String> {
    let compressed = common::verified_bytes(SUPPORT_GZIP, SUPPORT_GZIP_SHA)?;
    let mut raw = Vec::new();
    GzDecoder::new(&compressed[..])
        .read_to_end(&mut raw)
        .map_err(|error| format!("support gzip: {error}"))?;
    if hex::encode(sha256(&raw)) != SUPPORT_RAW_SHA {
        return Err("expanded support SHA-256 differs".into());
    }
    serde_json::from_slice(&raw).map_err(|error| format!("support JSON: {error}"))
}

fn map_to_leaf(ctx: &Context, point: FastPoint) -> Result<FastPoint, String> {
    velu_point_map(
        ctx.bridge.map_point(point).ok_or("field bridge point")?,
        &ctx.isogeny.kernel,
        37,
        &ctx.archived_source.irr,
    )
    .ok_or("undefined degree-73 map".into())
}

fn read_targets(path: &str, sha: &str, curve: &Curve) -> Result<Vec<FastPoint>, String> {
    let bytes = common::verified_bytes(path, sha)?;
    let text = String::from_utf8(bytes).map_err(|error| format!("{path}: UTF-8: {error}"))?;
    let mut points = Vec::new();
    for (index, line) in text.lines().enumerate() {
        let words: [u64; 2] =
            serde_json::from_str(line).map_err(|error| format!("{path}:{index}: {error}"))?;
        let point = FastPoint::affine(words[0], words[1]);
        if !curve.fast.is_on_curve(point) || !curve.fast.mul_u64(point, R).infinity {
            return Err(format!("{path}:{index}: invalid subgroup target"));
        }
        points.push(point);
    }
    if points.len() != 1024 {
        return Err(format!("{path}: expected 1024 public targets"));
    }
    Ok(points)
}

fn load_inputs(ctx: &Context) -> Result<Inputs, String> {
    let support = read_support()?;
    if support["schema"] != "n37-four-policy-support-v1"
        || support["status"] != "PASS"
        || support["r"].as_u64() != Some(R)
        || support["lambda"].as_u64() != Some(ctx.lambda)
    {
        return Err("support manifest identity changed".into());
    }
    let policies = support["policies"].as_array().ok_or("support policies")?;
    if policies.len() != 4 {
        return Err("support manifest does not have four policies".into());
    }
    let powers = ctx.lambda_powers();
    let expected_names = ["original", "transported", "descendant_native", "pullback"];
    let mut parsed = Vec::with_capacity(4);
    let mut scalar_muls = 0u64;
    for (policy_index, policy) in policies.iter().enumerate() {
        let entries = policy["entries"].as_array().ok_or("support entries")?;
        if policy["name"] != expected_names[policy_index]
            || policy["columns"].as_u64() != Some(K as u64)
            || policy["signed_classes"].as_u64() != Some((K * N) as u64)
            || policy["physical_points"].as_u64() != Some(P as u64)
            || entries.len() != P
        {
            return Err(format!("support policy {policy_index} metadata differs"));
        }
        let curve = if policy_index == 0 || policy_index == 3 {
            &ctx.control
        } else {
            &ctx.leaf
        };
        let mut factors = Vec::with_capacity(P);
        let mut packed = HashSet::with_capacity(P);
        let mut classes = HashSet::with_capacity(K * N);
        for (index, entry) in entries.iter().enumerate() {
            let point = point_from_value(&entry["point"])?;
            let column = index / (2 * N);
            let exponent = (index / 2) % N;
            let sign = if index % 2 == 0 { 1i64 } else { -1i64 };
            let coefficient = if sign == 1 {
                powers[exponent]
            } else {
                (R - powers[exponent]) % R
            };
            scalar_muls += 1;
            if entry["column"].as_u64() != Some(column as u64)
                || entry["exponent"].as_u64() != Some(exponent as u64)
                || entry["sign"].as_i64() != Some(sign)
                || entry["coefficient"].as_u64() != Some(coefficient)
                || point.infinity
                || point.x == 0
                || !curve.fast.is_on_curve(point)
                || !curve.fast.mul_u64(point, R).infinity
                || !packed.insert(point.pack())
            {
                return Err(format!("invalid policy {policy_index} factor {index}"));
            }
            if sign == 1 && !classes.insert(point.x) {
                return Err(format!("policy {policy_index} signed class collision"));
            }
            factors.push(Factor {
                point,
                column,
                coefficient,
            });
        }
        if classes.len() != K * N {
            return Err(format!("policy {policy_index} class shortfall"));
        }
        parsed.push(factors);
    }
    let mut map_calls = 0u64;
    for index in 0..P {
        map_calls += 2;
        if map_to_leaf(ctx, parsed[0][index].point)? != parsed[1][index].point
            || map_to_leaf(ctx, parsed[3][index].point)? != parsed[2][index].point
        {
            return Err(format!("support transport pairing failed at {index}"));
        }
    }
    let source_generator = ctx.control.fast.lift(&ctx.kc.generator());
    let leaf_generator = map_to_leaf(ctx, source_generator)?;
    map_calls += 1;
    if !ctx.control.fast.mul_u64(source_generator, R).infinity
        || !ctx.leaf.fast.mul_u64(leaf_generator, R).infinity
    {
        return Err("source or leaf generator has wrong order".into());
    }
    scalar_muls += 2;
    let mut source_targets = read_targets(B03, B03_SHA, &ctx.control)?;
    source_targets.extend(read_targets(B04, B04_SHA, &ctx.control)?);
    scalar_muls += source_targets.len() as u64;
    let mut leaf_targets = Vec::with_capacity(source_targets.len());
    for &target in &source_targets {
        let image = map_to_leaf(ctx, target)?;
        map_calls += 1;
        scalar_muls += 1;
        if image.infinity
            || !ctx.leaf.fast.is_on_curve(image)
            || !ctx.leaf.fast.mul_u64(image, R).infinity
        {
            return Err("mapped target failed leaf subgroup check".into());
        }
        leaf_targets.push(image);
    }
    Ok(Inputs {
        policies: parsed,
        source_targets,
        leaf_targets,
        source_generator,
        leaf_generator,
        validation_scalar_muls: scalar_muls,
        validation_map_calls: map_calls,
    })
}

fn log_from_witness(factors: &[Factor], indices: &[usize], base_logs: &[u64]) -> u64 {
    indices
        .iter()
        .map(|&index| factors[index].coefficient as u128 * base_logs[factors[index].column] as u128)
        .sum::<u128>()
        .rem_euclid(R as u128) as u64
}

fn run_policy(
    name: &str,
    ctx: &Context,
    curve: &Curve,
    factors: &[Factor],
    generator: FastPoint,
    inputs: &Inputs,
    scalars: &[u64],
    source_policy: bool,
) -> Result<Value, String> {
    let started = Instant::now();
    let build_started = Instant::now();
    let table = Table::build(curve, factors)?;
    let build_ms = build_started.elapsed().as_secs_f64() * 1_000.0;
    let mut ranker =
        IncrementalRelationSolver::new(K, &BigUint::from(R)).ok_or("ranker modulus")?;
    let mut independent = Vec::<(Vec<u64>, u64)>::new();
    let mut rank_trace = Vec::new();
    let mut rank_lookups = 0u64;
    let mut rank_subtractions = 0u64;
    let mut rank_witness_adds = 0u64;
    let rank_started = Instant::now();
    for &scalar in scalars {
        let point = curve.fast.mul_u64(generator, scalar);
        let witness = query(&table, curve, factors, point, 3)?;
        rank_lookups += witness.lookups as u64;
        rank_subtractions += witness.subtractions as u64;
        rank_witness_adds += witness.verification_adds as u64;
        let mut row_value = Value::Null;
        let status = if let Some(indices) = &witness.indices {
            let row = row_from_indices(factors, indices);
            let mut augmented = row.clone();
            augmented.push(0);
            augmented.push(scalar);
            let status = ranker.add_row(augmented);
            if status == RowStatus::Inconsistent {
                return Err(format!("{name}: verified relation contradicts rank system"));
            }
            if status == RowStatus::Independent {
                independent.push((row.clone(), scalar));
            }
            row_value = json!(row);
            format!("{status:?}")
        } else {
            "ProvedMiss".to_string()
        };
        rank_trace.push(json!({
            "scalar":scalar, "m3":witness.value(), "row":row_value,
            "row_status":status, "rank_after":ranker.rank(),
        }));
        if ranker.rank() == K {
            break;
        }
    }
    let rank_ms = rank_started.elapsed().as_secs_f64() * 1_000.0;

    let solve_started = Instant::now();
    let mut base_logs: Option<Vec<u64>> = None;
    if ranker.rank() == K {
        let mut matrix: Vec<Vec<BigUint>> = independent
            .iter()
            .map(|(row, _)| row.iter().copied().map(BigUint::from).collect())
            .collect();
        let mut rhs: Vec<BigUint> = independent
            .iter()
            .map(|(_, scalar)| BigUint::from(*scalar))
            .collect();
        let values = gaussian_eliminate_mod_n(&mut matrix, &mut rhs, &BigUint::from(R))
            .ok_or("full-rank dense solve failed")?;
        let logs: Vec<u64> = values
            .iter()
            .map(|value| value.to_u64_digits().first().copied().unwrap_or(0))
            .collect();
        if logs.len() != K {
            return Err("dense solve did not return 42 logs".into());
        }
        for (column, &log) in logs.iter().enumerate() {
            if curve.fast.mul_u64(generator, log) != factors[column * 2 * N].point {
                return Err(format!("{name}: base log {column} failed full-point check"));
            }
        }
        base_logs = Some(logs);
    }
    let solve_ms = solve_started.elapsed().as_secs_f64() * 1_000.0;

    let targets_started = Instant::now();
    let mut target_results = Vec::with_capacity(inputs.source_targets.len());
    let mut m2_hits = 0usize;
    let mut m3_hits = 0usize;
    let mut logs_verified = 0usize;
    let mut m2_lookups = 0u64;
    let mut m3_lookups = 0u64;
    let mut m2_subtractions = 0u64;
    let mut m3_subtractions = 0u64;
    let mut witness_adds = 0u64;
    for index in 0..inputs.source_targets.len() {
        let target = if source_policy {
            inputs.source_targets[index]
        } else {
            inputs.leaf_targets[index]
        };
        let m2 = query(&table, curve, factors, target, 2)?;
        let m3 = query(&table, curve, factors, target, 3)?;
        if m2.indices.is_some() && m3.indices.is_none() {
            return Err(format!(
                "{name}: m2 hit outside m3 support at target {index}"
            ));
        }
        m2_hits += usize::from(m2.indices.is_some());
        m3_hits += usize::from(m3.indices.is_some());
        m2_lookups += m2.lookups as u64;
        m3_lookups += m3.lookups as u64;
        m2_subtractions += m2.subtractions as u64;
        m3_subtractions += m3.subtractions as u64;
        witness_adds += (m2.verification_adds + m3.verification_adds) as u64;
        let recovered = match (&m3.indices, &base_logs) {
            (Some(indices), Some(logs)) => {
                let log = log_from_witness(factors, indices, logs);
                // Every answer is checked in both coupled groups, regardless
                // of the policy on which its witness was found.
                if ctx.control.fast.mul_u64(inputs.source_generator, log)
                    != inputs.source_targets[index]
                    || ctx.leaf.fast.mul_u64(inputs.leaf_generator, log)
                        != inputs.leaf_targets[index]
                {
                    return Err(format!(
                        "{name}: recovered target {index} failed full-point checks"
                    ));
                }
                logs_verified += 1;
                Some(log)
            }
            _ => None,
        };
        target_results.push(json!({
            "index":index, "m2":m2.value(), "m3":m3.value(),
            "recovered_log":recovered,
        }));
    }
    let targets_ms = targets_started.elapsed().as_secs_f64() * 1_000.0;
    let rank_probe_count = rank_trace.len();
    let base_logs_verified = base_logs.is_some();
    let b03_m2_hits = target_results[..1024]
        .iter()
        .filter(|row| row["m2"]["status"] == "hit")
        .count();
    let b04_m2_hits = target_results[1024..]
        .iter()
        .filter(|row| row["m2"]["status"] == "hit")
        .count();
    let b03_m3_hits = target_results[..1024]
        .iter()
        .filter(|row| row["m3"]["status"] == "hit")
        .count();
    let b04_m3_hits = target_results[1024..]
        .iter()
        .filter(|row| row["m3"]["status"] == "hit")
        .count();
    Ok(json!({
        "name":name, "curve":if source_policy { "source" } else { "leaf" },
        "table":table.summary(),
        "rank":ranker.rank(), "rank_probe_count":rank_probe_count,
        "rank_misses":rank_trace.iter().filter(|row| row["row_status"]=="ProvedMiss").count(),
        "dependent_rows":ranker.dependent_rows(),
        "rank_trace":rank_trace, "base_logs":base_logs,
        "targets":target_results,
        "summary":{
            "m2_hits":m2_hits, "m3_hits":m3_hits,
            "b03_m2_hits":b03_m2_hits, "b04_m2_hits":b04_m2_hits,
            "b03_m3_hits":b03_m3_hits, "b04_m3_hits":b04_m3_hits,
            "verified_target_logs":logs_verified,
        },
        "operation_counts":{
            "rank_probe_scalar_multiplications":rank_probe_count,
            "rank_lookups":rank_lookups, "rank_subtractions":rank_subtractions,
            "rank_witness_verification_adds":rank_witness_adds,
            "base_log_verification_scalar_multiplications":if base_logs_verified { K } else { 0 },
            "target_log_verification_scalar_multiplications":2*logs_verified,
            "m2_lookups":m2_lookups, "m3_lookups":m3_lookups,
            "m2_subtractions":m2_subtractions, "m3_subtractions":m3_subtractions,
            "target_witness_verification_adds":witness_adds,
        },
        "phase_ms_descriptive":{
            "build":build_ms, "rank":rank_ms, "solve_verify_base":solve_ms,
            "targets":targets_ms, "total":started.elapsed().as_secs_f64()*1000.0,
        },
    }))
}

fn paired_discordance(a: &[Value], b: &[Value], arity: &str) -> (usize, usize) {
    let mut a_only = 0usize;
    let mut b_only = 0usize;
    for (left, right) in a.iter().zip(b) {
        let hit_a = left[arity]["status"] == "hit";
        let hit_b = right[arity]["status"] == "hit";
        a_only += usize::from(hit_a && !hit_b);
        b_only += usize::from(hit_b && !hit_a);
    }
    (a_only, b_only)
}

fn mcnemar_exact(a_only: usize, b_only: usize) -> Value {
    let n = a_only + b_only;
    let k = a_only.min(b_only);
    let mut choose = BigUint::from(1u32);
    let mut tail = BigUint::from(0u32);
    for i in 0..=k {
        tail += &choose;
        if i < k {
            choose = choose * BigUint::from(n - i) / BigUint::from(i + 1);
        }
    }
    let denominator = BigUint::from(1u32) << n;
    let numerator = (tail * BigUint::from(2u32)).min(denominator.clone());
    json!({
        "a_only":a_only, "b_only":b_only,
        "two_sided_p_numerator":numerator.to_string(),
        "two_sided_p_denominator_power_of_two":n,
        "p_lt_0_01":numerator * BigUint::from(100u32) < denominator,
    })
}

#[cfg(unix)]
fn peak_rss_bytes() -> Option<u64> {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
    if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } != 0 {
        return None;
    }
    let rss = unsafe { usage.assume_init() }.ru_maxrss;
    if rss < 0 {
        return None;
    }
    #[cfg(target_os = "macos")]
    {
        Some(rss as u64)
    }
    #[cfg(not(target_os = "macos"))]
    {
        Some((rss as u64).saturating_mul(1024))
    }
}

#[cfg(not(unix))]
fn peak_rss_bytes() -> Option<u64> {
    None
}

fn run() -> Result<Value, String> {
    let started = Instant::now();
    let ctx = Context::load()?;
    let inputs = load_inputs(&ctx)?;
    let input_ms = started.elapsed().as_secs_f64() * 1_000.0;
    let scalars = probe_scalars();
    let names = ["original", "transported", "descendant_native", "pullback"];
    let mut outcomes = Vec::with_capacity(4);
    for (index, name) in names.into_iter().enumerate() {
        let source = index == 0 || index == 3;
        let curve = if source { &ctx.control } else { &ctx.leaf };
        let generator = if source {
            inputs.source_generator
        } else {
            inputs.leaf_generator
        };
        outcomes.push(run_policy(
            name,
            &ctx,
            curve,
            &inputs.policies[index],
            generator,
            &inputs,
            &scalars,
            source,
        )?);
    }
    for &(source, leaf) in &[(0usize, 1usize), (3, 2)] {
        if outcomes[source]["targets"] != outcomes[leaf]["targets"]
            || outcomes[source]["rank_trace"] != outcomes[leaf]["rank_trace"]
            || outcomes[source]["base_logs"] != outcomes[leaf]["base_logs"]
            || outcomes[source]["rank"] != outcomes[leaf]["rank"]
            || outcomes[source]["table"]["distinct_sums"]
                != outcomes[leaf]["table"]["distinct_sums"]
        {
            return Err(format!("paired policies {source}/{leaf} differ"));
        }
    }
    let original_targets = outcomes[0]["targets"]
        .as_array()
        .ok_or("original targets")?;
    let native_targets = outcomes[3]["targets"]
        .as_array()
        .ok_or("pullback targets")?;
    let (m2_a, m2_b) = paired_discordance(original_targets, native_targets, "m2");
    let (m3_a, m3_b) = paired_discordance(original_targets, native_targets, "m3");
    let m3_test = mcnemar_exact(m3_a, m3_b);
    let m3_gap = m3_a.abs_diff(m3_b);
    let m3_saturated =
        outcomes[0]["summary"]["m3_hits"] == 2048 && outcomes[3]["summary"]["m3_hits"] == 2048;
    let rank_complete = outcomes.iter().all(|row| row["rank"] == K);
    let all_logs_verified = outcomes
        .iter()
        .all(|row| row["summary"]["verified_target_logs"] == 2048);
    let decision = if !rank_complete {
        "RANK_INCOMPLETE"
    } else if m3_saturated {
        "M3_SATURATED_NONDISCRIMINATING"
    } else if m3_gap >= 21 && m3_test["p_lt_0_01"] == true && rank_complete {
        "FIXED_BLOCK_M3_SELECTION_LEAD"
    } else {
        "NO_FIXED_BLOCK_M3_SELECTION_LEAD"
    };
    Ok(json!({
        "schema":"n37-four-policy-pdp-v1", "status":"PASS",
        "protocol_commit":"10043b86",
        "protocol":format!("{NOTE}/PROTOCOL.md"),
        "support_gzip_sha256":SUPPORT_GZIP_SHA,
        "support_raw_sha256":SUPPORT_RAW_SHA,
        "point_only_target_sha256":{"b03":B03_SHA,"b04":B04_SHA},
        "target_count":2048, "columns":K, "physical_points_each":P,
        "rank_probe_seed":format!("0x{PROBE_SEED:016x}"), "rank_probe_cap":PROBE_CAP,
        "shared_input_validation_counts":{
            "scalar_multiplications":inputs.validation_scalar_muls,
            "degree73_map_calls":inputs.validation_map_calls,
        },
        "shared_input_ms_descriptive":input_ms,
        "policies":outcomes,
        "selection_comparison":{
            "m2":mcnemar_exact(m2_a,m2_b),
            "m3":m3_test,
            "m3_absolute_hit_gap":m3_gap,
            "m3_saturated":m3_saturated,
            "rank_42_each":rank_complete,
            "verified_2048_logs_each":all_logs_verified,
            "decision":decision,
        },
        "cost_status":"stage diagnostic; archived support excludes cold source/leaf selection and transport construction; S and rho ratio unset",
        "s":Value::Null, "rho_ratio":Value::Null,
        "peak_rss_before_output_bytes":peak_rss_bytes(),
        "total_ms_descriptive":started.elapsed().as_secs_f64()*1_000.0,
        "host":{"os":std::env::consts::OS,"arch":std::env::consts::ARCH},
    }))
}

fn main() {
    let output = std::env::args()
        .nth(1)
        .expect("usage: n37_four_policy_pdp OUTPUT.json");
    assert!(
        std::env::args().nth(2).is_none(),
        "unexpected extra argument"
    );
    let value = match run() {
        Ok(value) => value,
        Err(error) => json!({"schema":"n37-four-policy-pdp-v1","status":"FAIL","error":error}),
    };
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&output)
        .expect("refuse to overwrite PDP result");
    file.write_all(&serde_json::to_vec_pretty(&value).expect("serialize PDP result"))
        .expect("write PDP result");
    file.write_all(b"\n").expect("terminate PDP result");
    println!("{} {}", value["status"], output);
    if value["status"] != "PASS" {
        std::process::exit(1);
    }
}
