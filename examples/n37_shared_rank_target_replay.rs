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

fn key(inst: &crypto_lib::cryptanalysis::ic_boundary::BinaryInstance, value: &Value) -> Result<u64, String> {
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
    value[name].as_u64().ok_or_else(|| format!("missing {name}"))
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
    if fixture["inventory_sha256"] != sha256_hex(&serde_json::to_vec(inventory).map_err(|e| e.to_string())?) {
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
        let lines: Vec<_> = bytes.split(|&b| b == b'\n').filter(|line| !line.is_empty()).collect();
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
                let point: Value = serde_json::from_slice(line).map_err(|e| format!("{path}: {e}"))?;
                excluded.insert(key(&inst, &point)?);
            }
        }
    }
    let source = pinned(SOURCE, SOURCE_SHA)?;
    let header: Value = serde_json::from_slice(source.split(|&b| b == b'\n').next().unwrap())
        .map_err(|e| format!("source header: {e}"))?;
    let base = header["factor_base_point_coordinates"].as_array().ok_or("source base absent")?;
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
        excluded.insert(key(&inst, &json!([u64_at(probe, "point_x")?, u64_at(probe, "point_y")?]))?);
    }
    if u64_at(fixture, "excluded_orbits_before_selection")? != excluded.len() as u64 {
        return Err("initial exclusion count mismatch".into());
    }
    let target_rows = fixture["targets"].as_array().ok_or("targets absent")?;
    if target_rows.len() != 16 || u64_at(fixture, "target_count")? != 16 {
        return Err("target count mismatch".into());
    }
    let point_rows: Vec<_> = points.split(|&b| b == b'\n').filter(|line| !line.is_empty()).collect();
    if point_rows.len() != target_rows.len() {
        return Err("point-only file length mismatch".into());
    }
    let mut accepted = 0usize;
    let mut zero_rejections = 0u64;
    let mut orbit_rejections = 0u64;
    for candidate in 0..u64_at(fixture, "candidate_count")? {
        let digest = sha256(format!("n37-shared-rank-target-v1|{candidate}").as_bytes());
        let scalar = (BigUint::from_bytes_be(&digest) % BigUint::from(R))
            .to_u64().ok_or("candidate scalar outside u64")?;
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
        let row = target_rows.get(accepted).ok_or("more accepted candidates than target count")?;
        let point_only: Value = serde_json::from_slice(point_rows[accepted])
            .map_err(|e| format!("point-only row {accepted}: {e}"))?;
        if u64_at(row, "candidate_index")? != candidate
            || u64_at(row, "scalar")? != scalar
            || row["candidate_sha256"] != hex::encode(digest)
            || row["point"] != q_value
            || point_only != q_value
            || general_point(&point_only)? != q
        {
            return Err(format!("candidate or full-point replay mismatch at accepted {accepted}"));
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

fn target_replay(_raw: &Value, _fixture: &Value, _points: &[u8]) -> Result<Value, String> {
    Err("target transcript replay is pending implementation".into())
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
        return Ok(json!({"input":input,"fixture_sha256":sha256_hex(&fixture_bytes),"status":"PASS"}));
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
