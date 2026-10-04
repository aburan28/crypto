//! Independent general-group-law replay for the n37 point-only target gate.
//! The producer's fast curve, folded table and target recovery are never used.

use crypto_lib::binary_ecc::curve::{point_add, scalar_mul};
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ic_boundary::koblitz_instance;
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::collections::HashSet;
use std::fs;
use std::path::Path;

const SOURCE: &str = "research/notes/ecc2k130/n37_four_policy_support_20261003/SOURCE42.jsonl";
const SOURCE_SHA: &str = "0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75";
const RANK: &str = "research/notes/ecc2k130/shared_rank_folded_table_20261003/RAW.json";
const RANK_SHA: &str = "293d273c071c280326825c5a4f6d12bce00e46349735e5475b16bcd40264dd3b";
const R: u64 = 230_603_167;
const K: usize = 42;

fn read(path: &Path) -> Result<Vec<u8>, String> {
    fs::read(path).map_err(|error| format!("{}: {error}", path.display()))
}

fn pinned(path: &str, digest: &str) -> Result<Vec<u8>, String> {
    let bytes = read(Path::new(path))?;
    if sha256_hex(&bytes) != digest {
        return Err(format!("{path}: SHA-256 differs"));
    }
    Ok(bytes)
}

fn xy(value: &Value) -> Result<[u64; 2], String> {
    serde_json::from_value(value.clone()).map_err(|error| format!("point: {error}"))
}

fn general_point(value: &Value) -> Result<BinaryPoint, String> {
    let [x, y] = xy(value)?;
    Ok(BinaryPoint::Affine {
        x: F2mElement::from_biguint(&BigUint::from(x), 37),
        y: F2mElement::from_biguint(&BigUint::from(y), 37),
    })
}

fn key(
    inst: &crypto_lib::cryptanalysis::ic_boundary::BinaryInstance,
    value: &Value,
) -> Result<u64, String> {
    let [x, y] = xy(value)?;
    let p = crypto_lib::cryptanalysis::koblitz_fast::FastPoint::affine(x, y);
    if !inst.fast.is_on_curve(p) {
        return Err("off-curve point in frozen inventory".into());
    }
    let mut cursor = x;
    let mut minimum = x;
    for _ in 1..37 {
        cursor = inst.gf.sqr(cursor);
        minimum = minimum.min(cursor);
    }
    Ok(minimum)
}

fn u64_at(value: &Value, name: &str) -> Result<u64, String> {
    value[name]
        .as_u64()
        .ok_or_else(|| format!("missing {name}"))
}

fn fixture_replay(points: &[u8], fixture: &Value) -> Result<Value, String> {
    let inst = koblitz_instance(0, 37).ok_or("n37 curve unavailable")?;
    let kc = KoblitzCurve::new(0, 37).ok_or("general n37 curve unavailable")?;
    if fixture["schema"] != "n37-shared-rank-target-fixture-v1"
        || fixture["curve"] != inst.name
        || u64_at(fixture, "subgroup_order")? != R
        || fixture["generator"] != json!([inst.generator.x, inst.generator.y])
        || fixture["domain"] != "n37-shared-rank-target-v1"
        || fixture["source_sha256"] != SOURCE_SHA
        || fixture["rank_sha256"] != RANK_SHA
        || fixture["points_sha256"] != sha256_hex(points)
    {
        return Err("fixture identity or point-file digest mismatch".into());
    }
    let inventory = fixture["inventory"].as_array().ok_or("inventory absent")?;
    if inventory.is_empty() {
        return Err("empty inventory".into());
    }
    if fixture["inventory_sha256"]
        != sha256_hex(&serde_json::to_vec(inventory).map_err(|e| e.to_string())?)
    {
        return Err("inventory list digest mismatch".into());
    }
    let mut previous_path = "";
    let mut excluded = HashSet::new();
    let mut historical_rows = 0usize;
    for entry in inventory {
        let path = entry["path"].as_str().ok_or("inventory path absent")?;
        if !path.starts_with("research/notes/ecc2k130/")
            || !path.ends_with(".points.jsonl")
            || path <= previous_path
        {
            return Err("inventory paths not sorted and scoped".into());
        }
        previous_path = path;
        let bytes = read(Path::new(path))?;
        if entry["sha256"] != sha256_hex(&bytes) {
            return Err(format!("inventory digest mismatch: {path}"));
        }
        let lines: Vec<_> = bytes
            .split(|&b| b == b'\n')
            .filter(|line| !line.is_empty())
            .collect();
        if u64_at(entry, "rows")? != lines.len() as u64 {
            return Err(format!("inventory row count mismatch: {path}"));
        }
        historical_rows += lines.len();
        let is_n37 = Path::new(path)
            .file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.starts_with("n37_"));
        if is_n37 {
            for line in lines {
                let point: Value =
                    serde_json::from_slice(line).map_err(|e| format!("{path}: {e}"))?;
                excluded.insert(key(&inst, &point)?);
            }
        }
    }
    let source = pinned(SOURCE, SOURCE_SHA)?;
    let header: Value = serde_json::from_slice(source.split(|&b| b == b'\n').next().unwrap())
        .map_err(|e| format!("source header: {e}"))?;
    let base = header["factor_base_point_coordinates"]
        .as_array()
        .ok_or("source base absent")?;
    if base.len() != 3_108 {
        return Err("source point count mismatch".into());
    }
    for point in base {
        excluded.insert(key(&inst, point)?);
    }
    let rank: Value = serde_json::from_slice(&pinned(RANK, RANK_SHA)?)
        .map_err(|e| format!("rank transcript: {e}"))?;
    let probes = rank["relations"].as_array().ok_or("rank probes absent")?;
    if probes.len() != 55 || rank["verified"] != true {
        return Err("rank setup is not verified".into());
    }
    for probe in probes {
        excluded.insert(key(
            &inst,
            &json!([u64_at(probe, "point_x")?, u64_at(probe, "point_y")?]),
        )?);
    }
    if u64_at(fixture, "excluded_orbits_before_selection")? != excluded.len() as u64 {
        return Err("initial exclusion count mismatch".into());
    }
    let target_rows = fixture["targets"].as_array().ok_or("targets absent")?;
    if target_rows.len() != 16 || u64_at(fixture, "target_count")? != 16 {
        return Err("target count mismatch".into());
    }
    let point_rows: Vec<_> = points
        .split(|&b| b == b'\n')
        .filter(|line| !line.is_empty())
        .collect();
    if point_rows.len() != target_rows.len() {
        return Err("point-only file length mismatch".into());
    }
    let mut accepted = 0usize;
    let mut zero_rejections = 0u64;
    let mut orbit_rejections = 0u64;
    for candidate in 0..u64_at(fixture, "candidate_count")? {
        let digest = sha256(format!("n37-shared-rank-target-v1|{candidate}").as_bytes());
        let scalar = (BigUint::from_bytes_be(&digest) % BigUint::from(R))
            .to_u64()
            .ok_or("candidate scalar outside u64")?;
        if scalar == 0 {
            zero_rejections += 1;
            continue;
        }
        let q = scalar_mul(&kc.curve, kc.generator(), &BigUint::from(scalar));
        let fast = inst.fast.lift(&q);
        let q_value = json!([fast.x, fast.y]);
        if !excluded.insert(key(&inst, &q_value)?) {
            orbit_rejections += 1;
            continue;
        }
        let row = target_rows
            .get(accepted)
            .ok_or("more accepted candidates than target count")?;
        let point_only: Value = serde_json::from_slice(point_rows[accepted])
            .map_err(|e| format!("point-only row {accepted}: {e}"))?;
        if u64_at(row, "candidate_index")? != candidate
            || u64_at(row, "scalar")? != scalar
            || row["candidate_sha256"] != hex::encode(digest)
            || row["point"] != q_value
            || point_only != q_value
            || general_point(&point_only)? != q
        {
            return Err(format!(
                "candidate or full-point replay mismatch at accepted {accepted}"
            ));
        }
        accepted += 1;
    }
    if accepted != 16
        || u64_at(fixture, "zero_rejections")? != zero_rejections
        || u64_at(fixture, "orbit_rejections")? != orbit_rejections
    {
        return Err("selection or rejection counts differ".into());
    }
    Ok(json!({
        "status":"PASS",
        "points_sha256":sha256_hex(points),
        "inventory_sha256":fixture["inventory_sha256"],
        "inventory_files":inventory.len(),
        "historical_rows":historical_rows,
        "excluded_orbits":fixture["excluded_orbits_before_selection"],
        "targets":accepted,
        "candidate_count":fixture["candidate_count"],
        "zero_rejections":zero_rejections,
        "orbit_rejections":orbit_rejections,
    }))
}

fn mulmod(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % R as u128) as u64
}

fn modpow(mut base: u64, mut exponent: u64) -> u64 {
    let mut result = 1;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = mulmod(result, base);
        }
        base = mulmod(base, base);
        exponent >>= 1;
    }
    result
}

fn strip_timing(value: &mut Value) {
    match value {
        Value::Object(map) => {
            map.remove("wall_ns");
            for value in map.values_mut() {
                strip_timing(value);
            }
        }
        Value::Array(values) => {
            for value in values {
                strip_timing(value);
            }
        }
        _ => {}
    }
}

fn phase_cost(phase: &Value, name: &str) -> Result<(u64, f64), String> {
    let wall = u64_at(phase, "wall_ns")?;
    let adds = phase["group_ops"]["adds"]
        .as_u64()
        .ok_or(format!("{name} adds absent"))?;
    let doubles = phase["group_ops"]["doubles"]
        .as_u64()
        .ok_or(format!("{name} doubles absent"))?;
    let gae = phase["gae"].as_f64().ok_or(format!("{name} GAE absent"))?;
    // Default calibration prices these three mandatory native units at
    // one GAE each, even though optional field conversions remain unset.
    let native = &phase["native"];
    let expected = adds
        + doubles
        + native["lookups"].as_u64().unwrap_or(0)
        + native["row_ops"].as_u64().unwrap_or(0)
        + native["target_guard_probes"].as_u64().unwrap_or(0);
    if gae != expected as f64 {
        return Err(format!(
            "{name} default-calibration GAE differs from priced ledger"
        ));
    }
    Ok((wall, gae))
}

fn target_replay(raw: &Value, fixture: &Value, points: &[u8]) -> Result<Value, String> {
    let kc = KoblitzCurve::new(0, 37).ok_or("general n37 curve unavailable")?;
    let h = kc.cofactor.to_u64().ok_or("cofactor outside u64")?;
    let inverse_h = modpow(h % R, R - 2);
    if mulmod(h % R, inverse_h) != 1 {
        return Err("cofactor inverse check failed".into());
    }
    if raw["curve"] != "icv1-f2m37-tm534059-32aad96b"
        || raw["points_sha256"] != sha256_hex(points)
        || raw["verified"] != true
        || u64_at(&raw["spec"], "residual_seed")? != 202610031649
        || u64_at(&raw["spec"], "max_attempts")? != 64
    {
        return Err("target report identity, policy or verified flag differs".into());
    }
    let source = pinned(SOURCE, SOURCE_SHA)?;
    let header: Value = serde_json::from_slice(source.split(|&b| b == b'\n').next().unwrap())
        .map_err(|e| format!("source header: {e}"))?;
    let point_values = header["factor_base_point_coordinates"]
        .as_array()
        .ok_or("base points absent")?;
    let label_values = header["factor_base_point_labels"]
        .as_array()
        .ok_or("base labels absent")?;
    if point_values.len() != 3_108 || label_values.len() != point_values.len() {
        return Err("source base length differs".into());
    }
    let base_points: Vec<_> = point_values
        .iter()
        .map(general_point)
        .collect::<Result<_, _>>()?;
    let labels: Vec<[u64; 2]> = label_values.iter().map(xy).collect::<Result<_, _>>()?;
    let archived_rank: Value = serde_json::from_slice(&pinned(RANK, RANK_SHA)?)
        .map_err(|e| format!("archived rank: {e}"))?;
    let mut expected_rank = archived_rank.clone();
    let mut actual_rank = raw["rank"].clone();
    strip_timing(&mut expected_rank);
    strip_timing(&mut actual_rank);
    if expected_rank != actual_rank {
        return Err("rank setup differs from independently verified #1302 outside timing".into());
    }
    let logs: Vec<u64> = serde_json::from_value(archived_rank["column_logs"].clone())
        .map_err(|e| format!("archived column logs: {e}"))?;
    if logs.len() != K {
        return Err("archived rank is not 42 columns".into());
    }
    let target_rows = raw["targets"]
        .as_array()
        .ok_or("target transcript absent")?;
    let frozen_rows = fixture["targets"]
        .as_array()
        .ok_or("fixture targets absent")?;
    if target_rows.len() != 16 || frozen_rows.len() != 16 {
        return Err("target transcript length differs".into());
    }
    let mut total_gae = archived_rank["total_gae"]
        .as_f64()
        .ok_or("rank GAE absent")?;
    let mut total_attempts = 0usize;
    let mut direct_hits = 0usize;
    for (index, (target, frozen)) in target_rows.iter().zip(frozen_rows).enumerate() {
        let q = general_point(&frozen["point"])?;
        let reported = json!([u64_at(target, "point_x")?, u64_at(target, "point_y")?]);
        if reported != frozen["point"]
            || u64_at(target, "index")? != index as u64
            || !kc.curve.is_on_curve(&q)
            || scalar_mul(&kc.curve, &q, &BigUint::from(R)) != BinaryPoint::Infinity
        {
            return Err(format!("invalid target identity or subgroup at {index}"));
        }
        let attempts = target["attempts"]
            .as_array()
            .ok_or("target attempts absent")?;
        if attempts.is_empty() || attempts.len() > 64 {
            return Err(format!("target {index} has invalid attempt count"));
        }
        if target["query"]["native"]["attempts"].as_u64() != Some(attempts.len() as u64) {
            return Err(format!("target {index} attempt ledger differs"));
        }
        let mut rng = StdRng::seed_from_u64(202610031649 ^ index as u64);
        let mut last_log = None;
        let mut hits = 0u64;
        for (number, attempt) in attempts.iter().enumerate() {
            let expected_a = if number == 0 { 0 } else { rng.gen_range(1..R) };
            let shifted = if expected_a == 0 {
                q.clone()
            } else {
                let multiple = scalar_mul(&kc.curve, kc.generator(), &BigUint::from(expected_a));
                point_add(&kc.curve, &q, &multiple)
            };
            let residual = general_point(&json!([
                u64_at(attempt, "residual_x")?,
                u64_at(attempt, "residual_y")?
            ]))?;
            if u64_at(attempt, "number")? != number as u64
                || u64_at(attempt, "residual_scalar")? != expected_a
                || shifted != residual
            {
                return Err(format!(
                    "target {index} residual differs at attempt {number}"
                ));
            }
            if attempt["witness"].is_null() {
                if !attempt["row"].is_null() || number + 1 == attempts.len() {
                    return Err(format!("target {index} miss row or final status differs"));
                }
                continue;
            }
            hits += 1;
            if number + 1 != attempts.len() {
                return Err(format!("target {index} continued after hit"));
            }
            let witness: Vec<usize> = serde_json::from_value(attempt["witness"].clone())
                .map_err(|e| format!("target {index} witness: {e}"))?;
            if witness.is_empty() || witness.len() > 3 {
                return Err(format!("target {index} witness length differs"));
            }
            let mut sum = BinaryPoint::Infinity;
            let mut dense = vec![0u64; K];
            for &w in &witness {
                let point = base_points
                    .get(w)
                    .ok_or("target witness index outside base")?;
                sum = point_add(&kc.curve, &sum, point);
                let [column, coefficient] = labels[w];
                if column >= K as u64 || coefficient >= R {
                    return Err("target witness label outside base".into());
                }
                let col = column as usize;
                dense[col] = (dense[col] + coefficient) % R;
            }
            if sum != residual || attempt["row"] != json!(dense) {
                return Err(format!("target {index} witness sum or folded row differs"));
            }
            let projected = dense
                .iter()
                .zip(&logs)
                .fold(0u64, |sum, (&coefficient, &log)| {
                    (sum + mulmod(coefficient, log)) % R
                });
            let recovered = (mulmod(inverse_h, projected) + R - expected_a) % R;
            if target["recovered_log"].as_u64() != Some(recovered)
                || frozen["scalar"].as_u64() != Some(recovered)
                || scalar_mul(&kc.curve, kc.generator(), &BigUint::from(recovered)) != q
            {
                return Err(format!("target {index} recovered scalar differs"));
            }
            last_log = Some(recovered);
        }
        if hits != 1
            || last_log.is_none()
            || target["verified"] != true
            || target["relation_check"]["native"]["witnesses_checked"] != 1
            || target["recovery_check"]["native"]["scalar_replays"] != 1
        {
            return Err(format!(
                "target {index} recovery or verification count differs"
            ));
        }
        let mut wall_sum = 0u64;
        let mut gae_sum = 0f64;
        for name in [
            "query",
            "pdp",
            "relation_check",
            "descent",
            "recovery_check",
        ] {
            let (wall, gae) = phase_cost(&target[name], name)?;
            wall_sum += wall;
            gae_sum += gae;
        }
        if u64_at(target, "online_wall_ns")? != wall_sum
            || target["online_gae"].as_f64() != Some(gae_sum)
        {
            return Err(format!("target {index} exclusive phase sum differs"));
        }
        total_gae += gae_sum;
        total_attempts += attempts.len();
        direct_hits += usize::from(attempts.len() == 1);
    }
    if raw["total_gae"].as_f64() != Some(total_gae) {
        return Err("shared setup plus all target GAE differs".into());
    }
    Ok(json!({
        "status":"PASS",
        "rank_raw_sha256":RANK_SHA,
        "source_sha256":SOURCE_SHA,
        "targets_verified":target_rows.len(),
        "total_attempts":total_attempts,
        "direct_hits":direct_hits,
        "shared_setup_gae":archived_rank["total_gae"],
        "cold_gae_lower_bound":total_gae,
    }))
}

fn run(args: &[String]) -> Result<Value, String> {
    if args.len() != 4 && args.len() != 5 {
        return Err("usage: n37_shared_rank_target_replay POINTS FIXTURE [RAW] RECEIPT".into());
    }
    let points = read(Path::new(&args[1]))?;
    let fixture_bytes = read(Path::new(&args[2]))?;
    let fixture: Value = serde_json::from_slice(&fixture_bytes).map_err(|e| e.to_string())?;
    let input = fixture_replay(&points, &fixture)?;
    if args.len() == 4 {
        return Ok(
            json!({"input":input,"fixture_sha256":sha256_hex(&fixture_bytes),"status":"PASS"}),
        );
    }
    let raw_bytes = read(Path::new(&args[3]))?;
    let raw: Value = serde_json::from_slice(&raw_bytes).map_err(|e| e.to_string())?;
    let target = target_replay(&raw, &fixture, &points)?;
    Ok(json!({
        "input":input,
        "target":target,
        "fixture_sha256":sha256_hex(&fixture_bytes),
        "raw_sha256":sha256_hex(&raw_bytes),
        "status":"PASS",
    }))
}

fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().collect();
    let receipt_path = args.last().ok_or("receipt output absent")?;
    let result = run(&args);
    let receipt = match &result {
        Ok(value) => value.clone(),
        Err(error) => json!({"status":"FAIL","error":error}),
    };
    serde_json::to_writer_pretty(
        fs::File::create(receipt_path).map_err(|e| e.to_string())?,
        &receipt,
    )
    .map_err(|e| e.to_string())?;
    result.map(|_| ())
}
