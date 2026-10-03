//! Independent full-point hash-table and modular-rank replay for the n37 PDP gate.
//! Unlike the producer's sorted packed table, this recomputes complete source
//! supports in hash maps keyed by affine full points; every hit is re-added
//! with the general binary-curve law before any row or logarithm is accepted.

#![recursion_limit = "256"]

#[path = "support/n37_policy_common.rs"]
mod common;

use common::{point_from_value, Context, K, N, R};
use crypto_lib::binary_ecc::{curve::scalar_mul, BinaryPoint};
use crypto_lib::cryptanalysis::binary_velu::velu_point_map;
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::hash::sha256::sha256;
use flate2::read::GzDecoder;
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};
use std::fs::OpenOptions;
use std::io::{Read, Write};

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

struct HashTable {
    pairs: HashMap<FastPoint, u32>,
    singles: HashMap<FastPoint, u32>,
}

impl HashTable {
    fn build(ctx: &Context, factors: &[Factor]) -> Result<Self, String> {
        let mut pairs = HashMap::new();
        pairs
            .try_reserve(RAW_TABLE)
            .map_err(|error| format!("replay pair-table reservation: {error}"))?;
        pairs.insert(FastPoint::INFINITY, 0);
        for (index, factor) in factors.iter().enumerate() {
            pairs.entry(factor.point).or_insert((index + 1) as u32);
        }
        for i in 0..P {
            for j in i..P {
                let sum = ctx.control.fast.add(factors[i].point, factors[j].point);
                pairs.entry(sum).or_insert((P + 1 + i * P + j) as u32);
            }
        }
        let mut singles = HashMap::with_capacity(P + 1);
        singles.insert(FastPoint::INFINITY, 0);
        for (index, factor) in factors.iter().enumerate() {
            singles.entry(factor.point).or_insert((index + 1) as u32);
        }
        Ok(Self { pairs, singles })
    }

    fn get(&self, point: FastPoint, arity: usize) -> Option<u32> {
        if arity == 2 {
            self.singles.get(&point).copied()
        } else {
            self.pairs.get(&point).copied()
        }
    }
}

fn decode(code: u32) -> Vec<usize> {
    let value = code as usize;
    if value == 0 {
        Vec::new()
    } else if value <= P {
        vec![value - 1]
    } else {
        let pair = value - (P + 1);
        vec![pair / P, pair % P]
    }
}

fn general_sum(ctx: &Context, factors: &[Factor], indices: &[usize]) -> BinaryPoint {
    let mut sum = BinaryPoint::Infinity;
    for &index in indices {
        let point = ctx.control.fast.lower(factors[index].point);
        sum = ctx.kc.add(&sum, &point);
    }
    sum
}

fn query(
    ctx: &Context,
    table: &HashTable,
    factors: &[Factor],
    target: FastPoint,
    arity: usize,
) -> Result<Value, String> {
    let mut lookups = 0u32;
    let mut subtractions = 0u32;
    for scan in 0..=P {
        let residual = if scan == 0 {
            target
        } else {
            let index = scan - 1;
            subtractions += 1;
            ctx.control
                .fast
                .add(target, ctx.control.fast.neg(factors[index].point))
        };
        lookups += 1;
        if let Some(code) = table.get(residual, arity) {
            let mut indices = decode(code);
            if scan > 0 {
                indices.push(scan - 1);
            }
            if indices.len() > arity || indices.iter().any(|&index| index >= P) {
                return Err("replay hash table produced invalid factor indices".into());
            }
            if general_sum(ctx, factors, &indices) != ctx.control.fast.lower(target) {
                return Err("general group law rejects replay witness".into());
            }
            return Ok(json!({
                "status":"hit", "indices":indices, "arity":indices.len(),
                "lookups":lookups, "subtractions":subtractions,
                "verification_adds":indices.len(),
            }));
        }
    }
    Ok(json!({
        "status":"proved_miss", "indices":Value::Null, "arity":Value::Null,
        "lookups":lookups, "subtractions":subtractions,
        "verification_adds":0,
    }))
}

struct IndependentRank {
    pivots: Vec<Option<Vec<u64>>>,
    rank: usize,
    dependent: usize,
}

fn addmod(a: u64, b: u64) -> u64 {
    ((a as u128 + b as u128) % R as u128) as u64
}

fn submod(a: u64, b: u64) -> u64 {
    ((a as u128 + R as u128 - b as u128) % R as u128) as u64
}

fn mulmod(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % R as u128) as u64
}

fn powmod(mut base: u64, mut exponent: u64) -> u64 {
    let mut result = 1u64;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = mulmod(result, base);
        }
        base = mulmod(base, base);
        exponent >>= 1;
    }
    result
}

impl IndependentRank {
    fn new() -> Self {
        Self {
            pivots: vec![None; K],
            rank: 0,
            dependent: 0,
        }
    }

    fn add(&mut self, row: &[u64], rhs: u64) -> Result<&'static str, String> {
        if row.len() != K {
            return Err("replay rank row has wrong width".into());
        }
        let mut work = row.to_vec();
        work.push(rhs);
        for col in 0..K {
            if work[col] == 0 {
                continue;
            }
            if let Some(pivot) = &self.pivots[col] {
                let factor = work[col];
                for j in col..=K {
                    work[j] = submod(work[j], mulmod(factor, pivot[j]));
                }
            } else {
                let inverse = powmod(work[col], R - 2);
                if mulmod(work[col], inverse) != 1 {
                    return Err("non-invertible replay pivot".into());
                }
                for value in &mut work[col..=K] {
                    *value = mulmod(*value, inverse);
                }
                self.pivots[col] = Some(work);
                self.rank += 1;
                return Ok("Independent");
            }
        }
        if work[K] != 0 {
            return Err("independent replay found inconsistent verified row".into());
        }
        self.dependent += 1;
        Ok("Dependent")
    }

    fn solution(&self) -> Option<Vec<u64>> {
        if self.rank != K {
            return None;
        }
        let mut logs = vec![0u64; K];
        for col in (0..K).rev() {
            let row = self.pivots[col].as_ref()?;
            let mut value = row[K];
            for j in col + 1..K {
                value = submod(value, mulmod(row[j], logs[j]));
            }
            logs[col] = value;
        }
        Some(logs)
    }
}

fn row_from_indices(factors: &[Factor], indices: &[usize]) -> Vec<u64> {
    let mut row = vec![0u64; K];
    for &index in indices {
        let factor = &factors[index];
        row[factor.column] = addmod(row[factor.column], factor.coefficient);
    }
    row
}

fn log_from_indices(factors: &[Factor], indices: &[usize], logs: &[u64]) -> u64 {
    indices.iter().fold(0u64, |value, &index| {
        let factor = &factors[index];
        addmod(value, mulmod(factor.coefficient, logs[factor.column]))
    })
}

fn map_to_leaf(ctx: &Context, point: FastPoint) -> Result<FastPoint, String> {
    velu_point_map(
        ctx.bridge.map_point(point).ok_or("replay field bridge")?,
        &ctx.isogeny.kernel,
        37,
        &ctx.archived_source.irr,
    )
    .ok_or("replay degree-73 map undefined".into())
}

fn support_manifest() -> Result<Value, String> {
    let bytes = common::verified_bytes(SUPPORT_GZIP, SUPPORT_GZIP_SHA)?;
    let mut raw = Vec::new();
    GzDecoder::new(&bytes[..])
        .read_to_end(&mut raw)
        .map_err(|error| format!("replay support gzip: {error}"))?;
    if hex::encode(sha256(&raw)) != SUPPORT_RAW_SHA {
        return Err("replay support expanded SHA-256 differs".into());
    }
    serde_json::from_slice(&raw).map_err(|error| format!("replay support JSON: {error}"))
}

fn public_targets() -> Result<Vec<FastPoint>, String> {
    let mut targets = Vec::with_capacity(2048);
    for (path, sha) in [(B03, B03_SHA), (B04, B04_SHA)] {
        let bytes = common::verified_bytes(path, sha)?;
        let text = String::from_utf8(bytes).map_err(|error| format!("{path}: {error}"))?;
        let mut count = 0usize;
        for line in text.lines() {
            let words: [u64; 2] =
                serde_json::from_str(line).map_err(|error| format!("{path}: {error}"))?;
            targets.push(FastPoint::affine(words[0], words[1]));
            count += 1;
        }
        if count != 1024 {
            return Err(format!("{path}: expected 1024 targets"));
        }
    }
    Ok(targets)
}

fn policy_factors(ctx: &Context, support: &Value) -> Result<Vec<Vec<Factor>>, String> {
    let policies = support["policies"]
        .as_array()
        .ok_or("support policy array")?;
    if policies.len() != 4 || support["status"] != "PASS" {
        return Err("replay support policy identity differs".into());
    }
    let names = ["original", "transported", "descendant_native", "pullback"];
    let mut powers = Vec::with_capacity(N);
    let mut power = 1u64;
    for _ in 0..N {
        powers.push(power);
        power = mulmod(power, ctx.lambda);
    }
    if power != 1 {
        return Err("lambda action does not close at 37".into());
    }
    let mut output = Vec::with_capacity(4);
    for (policy_index, policy) in policies.iter().enumerate() {
        let entries = policy["entries"].as_array().ok_or("support entries")?;
        if policy["name"] != names[policy_index] || entries.len() != P {
            return Err("support policy name or size differs".into());
        }
        let curve = if policy_index == 0 || policy_index == 3 {
            &ctx.control
        } else {
            &ctx.leaf
        };
        let mut seen = HashSet::new();
        let mut factors = Vec::with_capacity(P);
        for (index, entry) in entries.iter().enumerate() {
            let point = point_from_value(&entry["point"])?;
            let column = index / (2 * N);
            let exponent = (index / 2) % N;
            let coefficient = if index % 2 == 0 {
                powers[exponent]
            } else {
                submod(0, powers[exponent])
            };
            let sign = if index % 2 == 0 { 1i64 } else { -1i64 };
            if entry["column"].as_u64() != Some(column as u64)
                || entry["exponent"].as_u64() != Some(exponent as u64)
                || entry["sign"].as_i64() != Some(sign)
                || entry["coefficient"].as_u64() != Some(coefficient)
                || point.infinity
                || !curve.fast.is_on_curve(point)
                || !curve.fast.mul_u64(point, R).infinity
                || !seen.insert(point)
            {
                return Err(format!(
                    "replay support factor {policy_index}/{index} differs"
                ));
            }
            factors.push(Factor {
                point,
                column,
                coefficient,
            });
        }
        output.push(factors);
    }
    for index in 0..P {
        if map_to_leaf(ctx, output[0][index].point)? != output[1][index].point
            || map_to_leaf(ctx, output[3][index].point)? != output[2][index].point
        {
            return Err(format!("replay paired factor differs at {index}"));
        }
    }
    Ok(output)
}

fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut x = *state;
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

fn probe_scalars() -> Vec<u64> {
    let mut state = PROBE_SEED;
    let mut seen = HashSet::new();
    let mut values = Vec::with_capacity(PROBE_CAP);
    while values.len() < PROBE_CAP {
        let scalar = 1 + splitmix64(&mut state) % (R - 1);
        if seen.insert(scalar) {
            values.push(scalar);
        }
    }
    values
}

fn exact_mcnemar(a_only: usize, b_only: usize) -> Value {
    let trials = a_only + b_only;
    let tail_end = a_only.min(b_only);
    let mut coefficient = BigUint::from(1u32);
    let mut tail = BigUint::from(1u32);
    for i in 1..=tail_end {
        coefficient = coefficient * BigUint::from(trials - i + 1) / BigUint::from(i);
        tail += &coefficient;
    }
    let denominator = BigUint::from(1u32) << trials;
    let numerator = (tail * BigUint::from(2u32)).min(denominator.clone());
    json!({
        "a_only":a_only, "b_only":b_only,
        "two_sided_p_numerator":numerator.to_string(),
        "two_sided_p_denominator_power_of_two":trials,
        "p_lt_0_01":numerator * BigUint::from(100u32) < denominator,
    })
}

fn query_sum(rows: &[Value], query_name: &str, field: &str) -> Result<u64, String> {
    rows.iter()
        .map(|row| {
            row[query_name][field]
                .as_u64()
                .ok_or_else(|| format!("replay query counter {query_name}/{field}"))
        })
        .sum()
}

fn verified_source_scalar(ctx: &Context, log: u64, point: FastPoint) -> bool {
    scalar_mul(&ctx.kc.curve, ctx.kc.generator(), &BigUint::from(log))
        == ctx.control.fast.lower(point)
}

fn replay_source_policy(
    ctx: &Context,
    manifest_policy: &Value,
    leaf_policy: &Value,
    factors: &[Factor],
    leaf_factors: &[Factor],
    targets: &[FastPoint],
    leaf_targets: &[FastPoint],
    scalars: &[u64],
) -> Result<(usize, usize, usize), String> {
    let table = HashTable::build(ctx, factors)?;
    for policy in [manifest_policy, leaf_policy] {
        if policy["table"]["raw_candidates"].as_u64() != Some(RAW_TABLE as u64)
            || policy["table"]["distinct_sums"].as_u64() != Some(table.pairs.len() as u64)
            || policy["table"]["collisions"].as_u64()
                != Some((RAW_TABLE - table.pairs.len()) as u64)
            || policy["table"]["identity_singleton_distinct_sums"].as_u64()
                != Some(table.singles.len() as u64)
            || policy["table"]["pair_group_additions"].as_u64() != Some((P * (P + 1) / 2) as u64)
            || policy["table"]["retained_bytes_lower_bound"]
                .as_u64()
                .is_none_or(|bytes| {
                    bytes
                        < ((table.pairs.len() + table.singles.len())
                            * std::mem::size_of::<(u64, u32)>()) as u64
                })
        {
            return Err("replay complete-table counts differ".into());
        }
    }
    let source_generator = ctx.control.fast.lift(ctx.kc.generator());
    let mut rank = IndependentRank::new();
    let trace = manifest_policy["rank_trace"]
        .as_array()
        .ok_or("rank trace")?;
    let leaf_trace = leaf_policy["rank_trace"]
        .as_array()
        .ok_or("leaf rank trace")?;
    if trace != leaf_trace {
        return Err("paired rank traces differ".into());
    }
    let mut replay_probe_count = 0usize;
    for (index, &scalar) in scalars.iter().enumerate() {
        if rank.rank == K {
            break;
        }
        replay_probe_count += 1;
        let target = ctx.control.fast.mul_u64(source_generator, scalar);
        let witness = query(ctx, &table, factors, target, 3)?;
        let mut row = Value::Null;
        let status = if let Some(indices) = witness["indices"].as_array() {
            let indices: Vec<usize> = indices
                .iter()
                .map(|v| v.as_u64().map(|n| n as usize).ok_or("witness index"))
                .collect::<Result<_, _>>()?;
            let coefficients = row_from_indices(factors, &indices);
            let kind = rank.add(&coefficients, scalar)?;
            row = json!(coefficients);
            kind
        } else {
            "ProvedMiss"
        };
        let expected = json!({
            "scalar":scalar, "m3":witness, "row":row,
            "row_status":status, "rank_after":rank.rank,
        });
        if trace.get(index) != Some(&expected) {
            return Err(format!("rank probe {index} differs in independent replay"));
        }
    }
    if trace.len() != replay_probe_count {
        return Err("rank trace length differs from independently replayed probes".into());
    }
    if rank.rank < K && trace.len() != PROBE_CAP {
        return Err("incomplete rank stopped before 256 probes".into());
    }
    let solved = rank.solution();
    if let Some(logs) = &solved {
        for (column, &log) in logs.iter().enumerate() {
            if !verified_source_scalar(ctx, log, factors[column * 2 * N].point)
                || map_to_leaf(ctx, factors[column * 2 * N].point)?
                    != leaf_factors[column * 2 * N].point
            {
                return Err(format!("replayed base log {column} failed general law"));
            }
        }
    }
    for policy in [manifest_policy, leaf_policy] {
        if policy["rank"].as_u64() != Some(rank.rank as u64)
            || policy["rank_probe_count"].as_u64() != Some(trace.len() as u64)
            || policy["dependent_rows"].as_u64() != Some(rank.dependent as u64)
            || policy["base_logs"] != json!(solved)
        {
            return Err("replay rank or base logs differ".into());
        }
    }

    let rows = manifest_policy["targets"]
        .as_array()
        .ok_or("source target rows")?;
    let paired_rows = leaf_policy["targets"]
        .as_array()
        .ok_or("leaf target rows")?;
    if rows.len() != 2048 || rows != paired_rows {
        return Err("paired target rows or count differ".into());
    }
    let mut m2_hits = 0usize;
    let mut m3_hits = 0usize;
    let mut verified_logs = 0usize;
    for (index, (&target, &leaf_target)) in targets.iter().zip(leaf_targets).enumerate() {
        let m2 = query(ctx, &table, factors, target, 2)?;
        let m3 = query(ctx, &table, factors, target, 3)?;
        m2_hits += usize::from(m2["status"] == "hit");
        m3_hits += usize::from(m3["status"] == "hit");
        let recovered = match (&solved, m3["indices"].as_array()) {
            (Some(logs), Some(indices)) => {
                let indices: Vec<usize> = indices
                    .iter()
                    .map(|v| v.as_u64().map(|n| n as usize).ok_or("target witness index"))
                    .collect::<Result<_, _>>()?;
                for &factor_index in &indices {
                    if map_to_leaf(ctx, factors[factor_index].point)?
                        != leaf_factors[factor_index].point
                    {
                        return Err(format!("target {index} leaf factor differs"));
                    }
                }
                if map_to_leaf(ctx, target)? != leaf_target {
                    return Err(format!("target {index} image differs"));
                }
                let log = log_from_indices(factors, &indices, logs);
                if !verified_source_scalar(ctx, log, target) {
                    return Err(format!("target {index} logarithm failed general law"));
                }
                verified_logs += 1;
                Some(log)
            }
            _ => None,
        };
        let expected = json!({
            "index":index, "m2":m2, "m3":m3, "recovered_log":recovered,
        });
        if rows[index] != expected {
            return Err(format!("target {index} exact PDP result differs"));
        }
    }
    for policy in [manifest_policy, leaf_policy] {
        let b03_m2_hits = rows[..1024]
            .iter()
            .filter(|row| row["m2"]["status"] == "hit")
            .count();
        let b04_m2_hits = rows[1024..]
            .iter()
            .filter(|row| row["m2"]["status"] == "hit")
            .count();
        let b03_m3_hits = rows[..1024]
            .iter()
            .filter(|row| row["m3"]["status"] == "hit")
            .count();
        let b04_m3_hits = rows[1024..]
            .iter()
            .filter(|row| row["m3"]["status"] == "hit")
            .count();
        if policy["summary"]["m2_hits"].as_u64() != Some(m2_hits as u64)
            || policy["summary"]["m3_hits"].as_u64() != Some(m3_hits as u64)
            || policy["summary"]["verified_target_logs"].as_u64() != Some(verified_logs as u64)
            || policy["summary"]["b03_m2_hits"].as_u64() != Some(b03_m2_hits as u64)
            || policy["summary"]["b04_m2_hits"].as_u64() != Some(b04_m2_hits as u64)
            || policy["summary"]["b03_m3_hits"].as_u64() != Some(b03_m3_hits as u64)
            || policy["summary"]["b04_m3_hits"].as_u64() != Some(b04_m3_hits as u64)
            || policy["rank_misses"].as_u64()
                != Some(
                    trace
                        .iter()
                        .filter(|row| row["row_status"] == "ProvedMiss")
                        .count() as u64,
                )
        {
            return Err("replay target yield or log counts differ".into());
        }
        let expected_counts = json!({
            "rank_probe_scalar_multiplications":trace.len(),
            "rank_lookups":query_sum(trace,"m3","lookups")?,
            "rank_subtractions":query_sum(trace,"m3","subtractions")?,
            "rank_witness_verification_adds":query_sum(trace,"m3","verification_adds")?,
            "base_log_verification_scalar_multiplications":if solved.is_some(){K}else{0},
            "target_log_verification_scalar_multiplications":2*verified_logs,
            "m2_lookups":query_sum(rows,"m2","lookups")?,
            "m3_lookups":query_sum(rows,"m3","lookups")?,
            "m2_subtractions":query_sum(rows,"m2","subtractions")?,
            "m3_subtractions":query_sum(rows,"m3","subtractions")?,
            "target_witness_verification_adds":query_sum(rows,"m2","verification_adds")?
                +query_sum(rows,"m3","verification_adds")?,
        });
        if policy["operation_counts"] != expected_counts {
            return Err("replay operation counters differ".into());
        }
    }
    Ok((m2_hits, m3_hits, verified_logs))
}

fn audit(manifest: &Value) -> Result<Value, String> {
    if manifest["schema"] != "n37-four-policy-pdp-v1"
        || manifest["status"] != "PASS"
        || manifest["protocol_commit"] != "10043b86"
        || manifest["support_gzip_sha256"] != SUPPORT_GZIP_SHA
        || manifest["support_raw_sha256"] != SUPPORT_RAW_SHA
        || manifest["point_only_target_sha256"]["b03"] != B03_SHA
        || manifest["point_only_target_sha256"]["b04"] != B04_SHA
        || manifest["target_count"].as_u64() != Some(2048)
        || manifest["columns"].as_u64() != Some(K as u64)
        || manifest["physical_points_each"].as_u64() != Some(P as u64)
        || manifest["rank_probe_seed"] != format!("0x{PROBE_SEED:016x}")
        || manifest["rank_probe_cap"].as_u64() != Some(PROBE_CAP as u64)
        || manifest["shared_input_validation_counts"]
            != json!({
                "scalar_multiplications":4*P+2+2*2048,
                "degree73_map_calls":2*P+1+2048,
            })
        || manifest["cost_status"]
            != "stage diagnostic; archived support excludes cold source/leaf selection and transport construction; S and rho ratio unset"
        || !manifest["s"].is_null()
        || !manifest["rho_ratio"].is_null()
    {
        return Err("PDP manifest or claim boundary differs".into());
    }
    let ctx = Context::load()?;
    let support = support_manifest()?;
    let factors = policy_factors(&ctx, &support)?;
    let targets = public_targets()?;
    let source_generator = ctx.control.fast.lift(ctx.kc.generator());
    let leaf_generator = map_to_leaf(&ctx, source_generator)?;
    if !ctx.control.fast.mul_u64(source_generator, R).infinity
        || !ctx.leaf.fast.mul_u64(leaf_generator, R).infinity
    {
        return Err("replay generator order differs".into());
    }
    let mut leaf_targets = Vec::with_capacity(targets.len());
    for (index, &point) in targets.iter().enumerate() {
        if !ctx.control.fast.is_on_curve(point) || !ctx.control.fast.mul_u64(point, R).infinity {
            return Err(format!("public point {index} is invalid"));
        }
        let image = map_to_leaf(&ctx, point)?;
        if image.infinity
            || !ctx.leaf.fast.is_on_curve(image)
            || !ctx.leaf.fast.mul_u64(image, R).infinity
        {
            return Err(format!("public image {index} is invalid"));
        }
        leaf_targets.push(image);
    }
    let policies = manifest["policies"].as_array().ok_or("PDP policy table")?;
    if policies.len() != 4 {
        return Err("PDP result lacks four policies".into());
    }
    let names = ["original", "transported", "descendant_native", "pullback"];
    for (index, name) in names.iter().enumerate() {
        if policies[index]["name"] != *name {
            return Err(format!("PDP policy {index} name differs"));
        }
    }
    let scalars = probe_scalars();
    let source = replay_source_policy(
        &ctx,
        &policies[0],
        &policies[1],
        &factors[0],
        &factors[1],
        &targets,
        &leaf_targets,
        &scalars,
    )?;
    let pullback = replay_source_policy(
        &ctx,
        &policies[3],
        &policies[2],
        &factors[3],
        &factors[2],
        &targets,
        &leaf_targets,
        &scalars,
    )?;
    let source_rows = policies[0]["targets"].as_array().ok_or("source rows")?;
    let pullback_rows = policies[3]["targets"].as_array().ok_or("pullback rows")?;
    let mut m2_source_only = 0usize;
    let mut m2_native_only = 0usize;
    let mut m3_source_only = 0usize;
    let mut m3_native_only = 0usize;
    for (left, right) in source_rows.iter().zip(pullback_rows) {
        for (arity, a, b) in [
            ("m2", &mut m2_source_only, &mut m2_native_only),
            ("m3", &mut m3_source_only, &mut m3_native_only),
        ] {
            let l = left[arity]["status"] == "hit";
            let r = right[arity]["status"] == "hit";
            *a += usize::from(l && !r);
            *b += usize::from(r && !l);
        }
    }
    let m2_test = exact_mcnemar(m2_source_only, m2_native_only);
    let m3_test = exact_mcnemar(m3_source_only, m3_native_only);
    let m3_gap = m3_source_only.abs_diff(m3_native_only);
    let saturated = source.1 == 2048 && pullback.1 == 2048;
    let full_rank = policies
        .iter()
        .all(|policy| policy["rank"].as_u64() == Some(K as u64));
    let all_logs = source.2 == 2048 && pullback.2 == 2048;
    let decision = if !full_rank {
        "RANK_INCOMPLETE"
    } else if saturated {
        "M3_SATURATED_NONDISCRIMINATING"
    } else if m3_gap >= 21 && m3_test["p_lt_0_01"] == true {
        "FIXED_BLOCK_M3_SELECTION_LEAD"
    } else {
        "NO_FIXED_BLOCK_M3_SELECTION_LEAD"
    };
    let expected_comparison = json!({
        "m2":m2_test, "m3":m3_test,
        "m3_absolute_hit_gap":m3_gap,
        "m3_saturated":saturated,
        "rank_42_each":full_rank,
        "verified_2048_logs_each":all_logs,
        "decision":decision,
    });
    if manifest["selection_comparison"] != expected_comparison {
        return Err("replayed source-native decision differs".into());
    }
    Ok(json!({
        "schema":"n37-four-policy-pdp-replay-v1", "status":"PASS",
        "recomputed_source_supports":2,
        "inferred_and_checked_isogenous_arms":2,
        "replayed_target_decisions":4*2048,
        "replayed_rank_trajectories":4,
        "source_policy_m2_hits":source.0, "source_policy_m3_hits":source.1,
        "source_policy_verified_logs":source.2,
        "native_policy_m2_hits":pullback.0, "native_policy_m3_hits":pullback.1,
        "native_policy_verified_logs":pullback.2,
    }))
}

fn main() {
    let mut args = std::env::args().skip(1);
    let input = args
        .next()
        .expect("usage: n37_four_policy_pdp_replay RESULT.json RECEIPT.json");
    let output = args.next().expect("missing replay receipt path");
    assert!(args.next().is_none(), "unexpected extra argument");
    let bytes = std::fs::read(&input).expect("read PDP result");
    let digest = hex::encode(sha256(&bytes));
    let manifest: Value = serde_json::from_slice(&bytes).expect("parse PDP result");
    let receipt = match audit(&manifest) {
        Ok(mut value) => {
            value["manifest_sha256"] = json!(digest);
            value
        }
        Err(error) => json!({
            "schema":"n37-four-policy-pdp-replay-v1", "status":"FAIL",
            "manifest_sha256":digest, "error":error,
        }),
    };
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&output)
        .expect("refuse to overwrite replay receipt");
    file.write_all(&serde_json::to_vec_pretty(&receipt).expect("serialize replay receipt"))
        .expect("write replay receipt");
    file.write_all(b"\n").expect("terminate replay receipt");
    println!("{} {}", receipt["status"], output);
    if receipt["status"] != "PASS" {
        std::process::exit(1);
    }
}
