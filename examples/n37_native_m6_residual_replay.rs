//! Independent source-curve replay of the descendant-native n37 residual run.
//! Builds at-most-3 reachable sets by iterative group-law closure, rather
//! than using the producer's ordered multiset table or its query routine.
//!
//! `cargo run --release --example n37_native_m6_residual_replay -- RAW.json RECEIPT.json`

use crypto_lib::binary_ecc::IrreduciblePoly;
use crypto_lib::cryptanalysis::binary_velu::Curve;
use crypto_lib::cryptanalysis::ec_index_calculus::gaussian_eliminate_mod_n;
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::fs;
use std::path::Path;

const BASE_PATH: &str = "research/notes/ecc2k130/n37_native_basis_bridge_20261002/NATIVE42.json";
const BASE_SHA: &str = "bd4bd8af982bcc65234ae1fae1ad52eb0ed41a1c203bf4d2556fd8e87f788a4c";
const FROZEN_PATH: &str = "research/notes/ecc2k130/disjoint_cold_v2_20261001/FROZEN.json";
const FROZEN_SHA: &str = "da958a3f1117dd1b88703fed2055c5c9f64255516e1e920f35d193bf53c6a88d";
const POINTS_PATH: &str =
    "research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b01.points.jsonl";
const POINTS_SHA: &str = "78553fdff5ae66521d3a1052978962258e78d48c92285df6e1df962034b43717";
const FIXTURE_PATH: &str =
    "research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b01.fixture.jsonl";
const FIXTURE_SHA: &str = "3453995c0423e4911ad4a6afa7cfe50ec1d6b0e6c67a771bda590f1d6ed680df";
const R: u64 = 230_603_167;
const K: usize = 42;
const SHIFT_SEED: u64 = 0x6e33_375f_7265_7331;
const SHIFT_CAP: usize = 16;

fn read_verified(path: &str, expected: &str) -> Result<Vec<u8>, String> {
    let bytes = fs::read(path).map_err(|error| format!("read {path}: {error}"))?;
    let actual = hex::encode(sha256(&bytes));
    if actual != expected {
        return Err(format!("SHA-256 mismatch for {path}: {actual}"));
    }
    Ok(bytes)
}

fn words(value: &Value) -> Result<[u64; 2], String> {
    serde_json::from_value(value.clone()).map_err(|error| format!("point words: {error}"))
}

fn point(value: &Value) -> Result<FastPoint, String> {
    let [x, y] = words(value)?;
    Ok(FastPoint::affine(x, y))
}

fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut value = *state;
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn frozen_shifts() -> Vec<u64> {
    let mut state = SHIFT_SEED;
    let mut seen = HashSet::new();
    let mut shifts = Vec::with_capacity(SHIFT_CAP);
    while shifts.len() < SHIFT_CAP {
        let candidate = 1 + splitmix64(&mut state) % (R - 1);
        if seen.insert(candidate) {
            shifts.push(candidate);
        }
    }
    shifts
}

fn coefficients(value: &Value) -> Result<Vec<i8>, String> {
    let vector: Vec<i8> = serde_json::from_value(value["coefficients"].clone())
        .map_err(|error| format!("witness coefficient vector: {error}"))?;
    if vector.len() != K {
        return Err("witness row width differs from K=42".into());
    }
    Ok(vector)
}

fn sum_coefficients(curve: &Curve, base: &[FastPoint], vector: &[i8]) -> FastPoint {
    let mut sum = FastPoint::INFINITY;
    for (&base_point, &coefficient) in base.iter().zip(vector) {
        let signed = if coefficient < 0 {
            curve.fast.neg(base_point)
        } else {
            base_point
        };
        for _ in 0..coefficient.unsigned_abs() {
            sum = curve.fast.add(sum, signed);
        }
    }
    sum
}

fn net_l1(vector: &[i8]) -> usize {
    vector.iter().map(|&c| c.unsigned_abs() as usize).sum()
}

fn reachable_sets(curve: &Curve, base: &[FastPoint]) -> Vec<HashSet<FastPoint>> {
    let signed: Vec<_> = base.iter().flat_map(|&p| [p, curve.fast.neg(p)]).collect();
    let mut reach = vec![HashSet::from([FastPoint::INFINITY])];
    for _ in 1..=3 {
        let previous = reach.last().expect("prior reachable set");
        let mut next = previous.clone();
        for &sum in previous {
            for &term in &signed {
                next.insert(curve.fast.add(sum, term));
            }
        }
        reach.push(next);
    }
    reach
}

fn member(
    curve: &Curve,
    target: FastPoint,
    left: &HashSet<FastPoint>,
    right: &HashSet<FastPoint>,
) -> bool {
    left.iter().any(|&sum| {
        let complement = curve.fast.add(target, curve.fast.neg(sum));
        right.contains(&complement)
    })
}

fn verify_witness(
    curve: &Curve,
    base: &[FastPoint],
    target: FastPoint,
    value: &Value,
    arity: usize,
) -> Result<Vec<i8>, String> {
    let vector = coefficients(value)?;
    if net_l1(&vector) > arity || sum_coefficients(curve, base, &vector) != target {
        return Err(format!("invalid {arity}-summand witness"));
    }
    Ok(vector)
}

fn run(result_path: &Path, receipt_path: &Path) -> Result<(), String> {
    let result_bytes = fs::read(result_path)
        .map_err(|error| format!("read {}: {error}", result_path.display()))?;
    let result_sha = hex::encode(sha256(&result_bytes));
    let result: Value =
        serde_json::from_slice(&result_bytes).map_err(|error| format!("raw JSON: {error}"))?;
    if result["schema"] != "n37-native-m6-residual-v1"
        || result["mode"] != "full"
        || result["target_count"].as_u64() != Some(1024)
    {
        return Err("replay requires the full frozen n37 result".into());
    }
    let base_manifest: Value = serde_json::from_slice(&read_verified(BASE_PATH, BASE_SHA)?)
        .map_err(|error| format!("base JSON: {error}"))?;
    let frozen: Value = serde_json::from_slice(&read_verified(FROZEN_PATH, FROZEN_SHA)?)
        .map_err(|error| format!("frozen JSON: {error}"))?;
    let spec = &frozen["specs"]["n37_L1024"];
    let modulus = IrreduciblePoly {
        degree: 37,
        low_terms: serde_json::from_value(spec["field_modulus_low_terms"].clone())
            .map_err(|error| format!("control modulus: {error}"))?,
    };
    let curve = Curve::with_order(37, &modulus, 0, 1, 137_439_487_532).ok_or("control curve")?;
    let generator = point(&spec["generator"])?;
    if !curve.fast.is_on_curve(generator) || !curve.fast.mul_u64(generator, R).infinity {
        return Err("invalid frozen generator".into());
    }
    let rows = base_manifest["rows"].as_array().ok_or("base rows")?;
    if rows.len() != K {
        return Err("base row count differs from K=42".into());
    }
    let mut base = Vec::new();
    for row in rows {
        let p = point(&row["control_pullback"])?;
        if !curve.fast.is_on_curve(p) || !curve.fast.mul_u64(p, R).infinity {
            return Err("invalid source pullback point".into());
        }
        base.push(p);
    }
    let reach = reachable_sets(&curve, &base);
    let relations = result["relations"].as_array().ok_or("relation records")?;
    let mut matrix = Vec::new();
    let mut rhs = Vec::new();
    let mut relation_misses = 0usize;
    for (index, row) in relations.iter().enumerate() {
        if row["probe_index"].as_u64() != Some(index as u64) {
            return Err("relation probe index mismatch".into());
        }
        let scalar = row["scalar"].as_u64().ok_or("relation scalar")?;
        if scalar == 0 || scalar >= R {
            return Err("relation scalar out of range".into());
        }
        let target = curve.fast.mul_u64(generator, scalar);
        let exact_member = member(&curve, target, &reach[3], &reach[3]);
        if row["witness"].is_null() {
            relation_misses += 1;
            if exact_member {
                return Err(format!("relation probe {index} falsely reported a miss"));
            }
            continue;
        }
        if !exact_member {
            return Err(format!("relation probe {index} falsely reported support"));
        }
        let vector = verify_witness(&curve, &base, target, &row["witness"], 6)?;
        matrix.push(
            vector
                .iter()
                .map(|&c| BigUint::from((c as i64).rem_euclid(R as i64) as u64))
                .collect(),
        );
        rhs.push(BigUint::from(scalar));
    }
    let solved = gaussian_eliminate_mod_n(&mut matrix, &mut rhs, &BigUint::from(R))
        .ok_or("independent dense relation matrix did not have full rank")?;
    if solved.len() != K || result["rank"].as_u64() != Some(K as u64) {
        return Err("independent rank differs from K=42".into());
    }
    let logs: Vec<u64> = solved
        .iter()
        .map(|value| value.to_u64_digits().first().copied().unwrap_or(0))
        .collect();
    let reported_logs: Vec<u64> = serde_json::from_value(result["base_logs"].clone())
        .map_err(|error| format!("reported base logs: {error}"))?;
    if logs != reported_logs {
        return Err("independently solved base logs differ".into());
    }
    for (&base_point, &log) in base.iter().zip(&logs) {
        if curve.fast.mul_u64(generator, log) != base_point {
            return Err("independent base logarithm check failed".into());
        }
    }

    let public: Vec<[u64; 2]> = String::from_utf8(read_verified(POINTS_PATH, POINTS_SHA)?)
        .map_err(|error| format!("point text: {error}"))?
        .lines()
        .map(|line| serde_json::from_str(line).map_err(|error| error.to_string()))
        .collect::<Result<_, _>>()?;
    let fixture: Vec<Value> = String::from_utf8(read_verified(FIXTURE_PATH, FIXTURE_SHA)?)
        .map_err(|error| format!("fixture text: {error}"))?
        .lines()
        .map(|line| serde_json::from_str(line).map_err(|error| error.to_string()))
        .collect::<Result<_, _>>()?;
    let targets = result["targets"].as_array().ok_or("target records")?;
    if public.len() != 1024 || fixture.len() != 1024 || targets.len() != 1024 {
        return Err("full target cardinality differs from 1024".into());
    }
    if result["shift_seed"] != format!("0x{SHIFT_SEED:016x}")
        || result["shift_cap"].as_u64() != Some(SHIFT_CAP as u64)
        || result["scalar_multiplications"]["shift_precomputation"].as_u64()
            != Some(SHIFT_CAP as u64)
    {
        return Err("frozen shift configuration differs".into());
    }
    let shifts = frozen_shifts();
    let reported_shifts = result["shifts"].as_array().ok_or("shift records")?;
    if reported_shifts.len() != SHIFT_CAP {
        return Err("shift list cardinality differs".into());
    }
    for (index, (&shift, row)) in shifts.iter().zip(reported_shifts).enumerate() {
        if row["scalar"].as_u64() != Some(shift) {
            return Err(format!("shift {index} differs from frozen schedule"));
        }
        let _ = words(&row["leaf_point"])?;
    }
    let mut direct_hits = 0usize;
    let mut shifted_hits = 0usize;
    let mut verified_logs = 0usize;
    let mut total_attempts = 0usize;
    let mut maximum_attempts = 0usize;
    for (index, target_row) in targets.iter().enumerate() {
        let source = FastPoint::affine(public[index][0], public[index][1]);
        if target_row["index"].as_u64() != Some(index as u64)
            || words(&target_row["source"])? != public[index]
            || words(&fixture[index]["published_q"])? != public[index]
        {
            return Err(format!("target {index} does not match frozen public input"));
        }
        let attempts = target_row["attempts"].as_array().ok_or("target attempts")?;
        if attempts.is_empty() || attempts.len() > SHIFT_CAP + 1 {
            return Err(format!("target {index} has invalid attempt count"));
        }
        total_attempts += attempts.len();
        maximum_attempts = maximum_attempts.max(attempts.len());
        let mut recovered = false;
        for (attempt_index, attempt) in attempts.iter().enumerate() {
            let shift = if attempt_index == 0 {
                0
            } else {
                shifts[attempt_index - 1]
            };
            if attempt["attempt_index"].as_u64() != Some(attempt_index as u64)
                || attempt["shift"].as_u64() != Some(shift)
            {
                return Err(format!("target {index} attempt schedule differs"));
            }
            let _ = words(&attempt["query_leaf"])?;
            if attempt_index == 0 && attempt["query_leaf"] != target_row["leaf"] {
                return Err(format!("target {index} direct leaf point differs"));
            }
            let query = if shift == 0 {
                source
            } else {
                curve.fast.add(source, curve.fast.mul_u64(generator, shift))
            };
            let supported = member(&curve, query, &reach[3], &reach[3]);
            if attempt["witness"].is_null() == supported {
                return Err(format!(
                    "target {index} attempt {attempt_index} support decision differs"
                ));
            }
            if supported {
                if attempt_index + 1 != attempts.len() {
                    return Err(format!("target {index} continued after a verified hit"));
                }
                let vector = verify_witness(&curve, &base, query, &attempt["witness"], 6)?;
                let shifted_log = vector
                    .iter()
                    .zip(&logs)
                    .map(|(&c, &log)| c as i128 * log as i128)
                    .sum::<i128>()
                    .rem_euclid(R as i128);
                let independently_derived =
                    (shifted_log - shift as i128).rem_euclid(R as i128) as u64;
                let reported = target_row["recovered_log"]
                    .as_u64()
                    .ok_or("supported target lacks reported logarithm")?;
                if reported != independently_derived
                    || curve.fast.mul_u64(generator, reported) != source
                    || fixture[index]["published_fixture_scalar"].as_u64() != Some(reported)
                {
                    return Err(format!("target {index} scalar failed independent check"));
                }
                if attempt_index == 0 {
                    direct_hits += 1;
                } else {
                    shifted_hits += 1;
                }
                verified_logs += 1;
                recovered = true;
            } else if attempt_index + 1 == attempts.len() && attempts.len() < SHIFT_CAP + 1 {
                return Err(format!("target {index} stopped before the shift cap"));
            }
        }
        if !recovered && !target_row["recovered_log"].is_null() {
            return Err(format!("missed target {index} reports a logarithm"));
        }
    }
    if result["direct_hits"].as_u64() != Some(direct_hits as u64)
        || result["shifted_hits"].as_u64() != Some(shifted_hits as u64)
        || result["unresolved_targets"].as_u64() != Some((1024 - verified_logs) as u64)
        || result["verified_target_logs"].as_u64() != Some(verified_logs as u64)
        || result["total_target_attempts"].as_u64() != Some(total_attempts as u64)
        || result["maximum_target_attempts"].as_u64() != Some(maximum_attempts as u64)
        || result["residual_shift_adds"].as_u64() != Some((total_attempts - 1024) as u64)
        || result["relation_misses"].as_u64() != Some(relation_misses as u64)
    {
        return Err("aggregate counts differ from independent replay".into());
    }
    let receipt = json!({
        "schema":"n37-native-m6-residual-source-replay-v1",
        "raw_result_sha256":result_sha,
        "input_sha256":{
            "base":BASE_SHA,"frozen":FROZEN_SHA,"points":POINTS_SHA,"fixture":FIXTURE_SHA
        },
        "method":"Iterative source-curve reachable sets through 3 factors; exact 3+3 membership for every direct and frozen-shift query; full source pullback witness replay; independent dense modular solve; source scalar checks and post-solve fixture equality",
        "reachable_source_sums_at_most":[reach[0].len(),reach[1].len(),reach[2].len(),reach[3].len()],
        "relation_probes":relations.len(),"relation_misses":relation_misses,
        "independent_rank":K,"base_logs_verified":K,
        "shift_seed":format!("0x{SHIFT_SEED:016x}"),"shift_cap":SHIFT_CAP,
        "direct_hits":direct_hits,"shifted_hits":shifted_hits,
        "unresolved_targets":1024-verified_logs,
        "total_target_attempts":total_attempts,"maximum_target_attempts":maximum_attempts,
        "verified_target_logs":verified_logs,"complete_batch":verified_logs==1024
    });
    let bytes = serde_json::to_vec_pretty(&receipt).map_err(|error| error.to_string())?;
    fs::write(receipt_path, bytes)
        .map_err(|error| format!("write {}: {error}", receipt_path.display()))?;
    println!(
        "replay rank={} direct={} shifted={} logs={} receipt={}",
        K,
        direct_hits,
        shifted_hits,
        verified_logs,
        receipt_path.display()
    );
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: n37_native_m6_residual_replay RAW.json RECEIPT.json");
        std::process::exit(2);
    }
    if let Err(error) = run(Path::new(&args[1]), Path::new(&args[2])) {
        eprintln!("n37_native_m6_residual_replay: {error}");
        std::process::exit(1);
    }
}
