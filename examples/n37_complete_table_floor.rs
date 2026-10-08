//! Count the mandatory pair-sum construction of the four frozen n37 policies.
//! This is a cold setup lower bound for the exact complete m<=3 oracle, not a
//! logarithm solver or an IC/rho speed claim. The protocol was committed in
//! PR #1280 before this producer or the new rho session was run.

#[path = "support/n37_policy_common.rs"]
#[allow(dead_code)]
mod common;

use common::{point_from_value, verified_bytes, Context, K, N, R};
use crypto_lib::cryptanalysis::ic_boundary::{BinaryGroup, CountedGroup, GroupOps};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::hash::sha256::sha256;
use flate2::read::GzDecoder;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::fs::OpenOptions;
use std::io::{Read, Write};

const SUPPORT_PATH: &str =
    "research/notes/ecc2k130/n37_four_policy_support_20261003/RESULT.json.gz";
const SUPPORT_SHA: &str = "8eac2ae4b8d8f0fc452b7cd7c8edbd3183558e55454abbd67349665bd64610d5";
const P: usize = 2 * K * N;
const PAIRS: usize = P * (P + 1) / 2;

fn support() -> Result<Value, String> {
    let compressed = verified_bytes(SUPPORT_PATH, SUPPORT_SHA)?;
    let mut reader = GzDecoder::new(&compressed[..]);
    let mut raw = Vec::new();
    reader
        .read_to_end(&mut raw)
        .map_err(|error| format!("support gzip: {error}"))?;
    let value: Value =
        serde_json::from_slice(&raw).map_err(|error| format!("support JSON: {error}"))?;
    if value["schema"] != "n37-four-policy-support-v1"
        || value["status"] != "PASS"
        || value["r"].as_u64() != Some(R)
        || value["source_base_hash"].as_str() != Some(common::SOURCE_BASE_HASH)
        || value["source_sha256"].as_str() != Some(common::SOURCE_SHA)
        || value["native_sha256"].as_str() != Some(common::NATIVE_SHA)
    {
        return Err("frozen support identity changed".into());
    }
    Ok(value)
}

fn count_policy(
    policy: &Value,
    expected_name: &str,
    expected_curve: &str,
    group: BinaryGroup<'_>,
) -> Result<Value, String> {
    if policy["name"] != expected_name
        || policy["curve"] != expected_curve
        || policy["columns"].as_u64() != Some(K as u64)
        || policy["signed_classes"].as_u64() != Some((K * N) as u64)
        || policy["physical_points"].as_u64() != Some(P as u64)
    {
        return Err(format!("{expected_name}: wrong policy metadata"));
    }
    let entries = policy["entries"]
        .as_array()
        .ok_or_else(|| format!("{expected_name}: no entries"))?;
    if entries.len() != P {
        return Err(format!("{expected_name}: expected {P} physical points"));
    }
    let mut seen = HashSet::with_capacity(P);
    let mut points = Vec::with_capacity(P);
    for (index, entry) in entries.iter().enumerate() {
        let point = point_from_value(&entry["point"])?;
        if point.infinity || !group.0.is_on_curve(point) || !seen.insert(point.pack()) {
            return Err(format!(
                "{expected_name}: invalid or duplicate point {index}"
            ));
        }
        points.push(point);
    }

    // The frozen PDP producer builds every i<=j pair using fast.add. We
    // replay the same complete enumeration through the ecbench group ledger.
    // Hashing every result makes both the arithmetic and its order auditable.
    let mut packed_sums = Vec::with_capacity(PAIRS * 8);
    let mut ops = GroupOps::default();
    for (i, &left) in points.iter().enumerate() {
        for &right in &points[i..] {
            let sum: FastPoint = group.add(&mut ops, left, right);
            packed_sums.extend_from_slice(&sum.pack().to_le_bytes());
        }
    }
    if ops.adds != PAIRS as u64
        || ops.doubles != 0
        || ops.scalar_mults != 0
        || packed_sums.len() != PAIRS * 8
    {
        return Err(format!("{expected_name}: counted pair ledger changed"));
    }
    Ok(json!({
        "policy": expected_name,
        "curve": expected_curve,
        "physical_points": P,
        "pair_additions": ops.adds,
        "pair_doublings": ops.doubles,
        "pair_scalar_mults": ops.scalar_mults,
        "pair_gae": ops.gae(),
        "ordered_pair_sums_sha256": hex::encode(sha256(&packed_sums)),
        "unpriced": ["point parsing", "curve and uniqueness checks", "table sorting", "table memory", "field arithmetic", "factor-base and isogeny setup", "relations", "linear algebra", "descent", "recovery"],
    }))
}

fn run() -> Result<Value, String> {
    let ctx = Context::load()?;
    let frozen = support()?;
    let policies = frozen["policies"]
        .as_array()
        .ok_or("support policies missing")?;
    if policies.len() != 4 {
        return Err("expected four policies".into());
    }
    let expected = [
        ("original", "source"),
        ("transported", "leaf"),
        ("descendant_native", "leaf"),
        ("pullback", "source"),
    ];
    let mut results = Vec::with_capacity(4);
    for (policy, (name, curve)) in policies.iter().zip(expected) {
        let group = match curve {
            "source" => BinaryGroup(&ctx.control.fast),
            "leaf" => BinaryGroup(&ctx.leaf.fast),
            _ => return Err("unknown policy curve".into()),
        };
        results.push(count_policy(policy, name, curve, group)?);
    }
    let r_sqrt = (R as f64).sqrt();
    Ok(json!({
        "schema": "n37-complete-table-floor/v1",
        "status": "PASS",
        "support_gzip_sha256": SUPPORT_SHA,
        "curve_slug": "icv1-f2m37-tm534059-32aad96b",
        "r": R,
        "signed_automorphism_order": 74,
        "columns_each": K,
        "signed_classes_each": K * N,
        "physical_points_each": P,
        "pairs_each": PAIRS,
        "cold_setup_gae_floor_each": PAIRS,
        "cold_setup_s_floor_each": PAIRS as f64 / r_sqrt,
        "generic_rho_s_floor": (std::f64::consts::PI / (2.0 * 74.0)).sqrt(),
        "scope": "mandatory table-construction group operations for the frozen complete at-most-three-summand oracle; no complete IC cost or speedup",
        "policies": results,
    }))
}

fn main() {
    let mut args = std::env::args().skip(1);
    let output = args
        .next()
        .expect("usage: n37_complete_table_floor OUTPUT.json");
    assert!(args.next().is_none(), "unexpected extra argument");
    let result = run().expect("count frozen n37 complete tables");
    let mut file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&output)
        .expect("refuse to overwrite count receipt");
    serde_json::to_writer_pretty(&mut file, &result).expect("write count receipt");
    file.write_all(b"\n").expect("newline");
    println!("{output}");
}
