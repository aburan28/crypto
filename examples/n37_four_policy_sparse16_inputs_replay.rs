//! Independent scalar and acceptance replay of the sparse n37 target freeze.

#![recursion_limit = "256"]

#[path = "support/n37_policy_common.rs"]
#[allow(dead_code)]
mod common;

use common::{Context, R};
use crypto_lib::binary_ecc::{curve::scalar_mul, BinaryPoint};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;

const PRIOR: [(&str, &str); 2] = [
    (
        "research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b03.points.jsonl",
        "84400a2914f06e4a001d0f113f0195952f692634d4ff8285e23ade2599e2bde2",
    ),
    (
        "research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b04.points.jsonl",
        "72c5b5361a0b846c9626269ce3928820cb7d420aae479d31dc26f5cbf7e78aa6",
    ),
];

fn key(ctx: &Context, point: FastPoint) -> Result<u64, String> {
    if point.infinity || !ctx.control.fast.is_on_curve(point) {
        return Err("bad finite point in orbit key".into());
    }
    let mut value = point.x;
    let mut minimum = value;
    for _ in 0..37 {
        minimum = minimum.min(value);
        value = ctx.control.gf.sqr(value);
    }
    if value != point.x {
        return Err("orbit did not close".into());
    }
    Ok(minimum)
}

fn checked_file(path: &Path, expected: &str) -> Result<Vec<u8>, String> {
    let bytes = fs::read(path).map_err(|error| format!("read {}: {error}", path.display()))?;
    if hex::encode(sha256(&bytes)) != expected {
        return Err(format!("SHA-256 mismatch: {}", path.display()));
    }
    Ok(bytes)
}

fn point_lines(bytes: &[u8]) -> Result<Vec<[u64; 2]>, String> {
    bytes
        .split(|&byte| byte == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| serde_json::from_slice(line).map_err(|error| error.to_string()))
        .collect()
}

fn replay(directory: &Path) -> Result<Value, String> {
    let ctx = Context::load()?;
    let freeze: Value = serde_json::from_slice(
        &fs::read(directory.join("FROZEN.json"))
            .map_err(|error| format!("read input freeze: {error}"))?,
    )
    .map_err(|error| format!("parse input freeze: {error}"))?;
    if freeze["schema"] != "n37-four-policy-sparse16-input-v1"
        || freeze["status"] != "FROZEN"
        || freeze["subgroup_order"].as_u64() != Some(R)
        || freeze["prior_rows"].as_u64() != Some(2048)
        || freeze["prior_point_sha256"]["b03"] != PRIOR[0].1
        || freeze["prior_point_sha256"]["b04"] != PRIOR[1].1
    {
        return Err("freeze identity differs".into());
    }
    let generator = ctx.control.fast.lift(ctx.kc.generator());
    if freeze["generator"] != json!([generator.x, generator.y]) {
        return Err("freeze generator differs".into());
    }
    let mut seen = HashSet::new();
    let mut prior_rows = 0usize;
    for (path, sha) in PRIOR {
        let data = common::verified_bytes(path, sha)?;
        for [x, y] in point_lines(&data)? {
            let point = FastPoint::affine(x, y);
            if scalar_mul(
                &ctx.kc.curve,
                &ctx.control.fast.lower(point),
                &BigUint::from(R),
            ) != BinaryPoint::Infinity
            {
                return Err("prior point failed general subgroup check".into());
            }
            seen.insert(key(&ctx, point)?);
            prior_rows += 1;
        }
    }
    if prior_rows != 2048 {
        return Err("prior row count differs".into());
    }
    let blocks = freeze["blocks"].as_array().ok_or("freeze blocks")?;
    if blocks.len() != 2 {
        return Err("freeze block count differs".into());
    }
    let mut verified = 0usize;
    for (block, spec) in blocks.iter().enumerate() {
        let point_name = format!("points-b{block}.jsonl");
        let label_name = format!("labels-b{block}.jsonl");
        if spec["block"].as_u64() != Some(block as u64)
            || spec["rows"].as_u64() != Some(1024)
            || spec["point_file"] != point_name
            || spec["label_file"] != label_name
        {
            return Err(format!("block {block} metadata differs"));
        }
        let point_sha = spec["point_sha256"].as_str().ok_or("point digest")?;
        let label_sha = spec["label_sha256"].as_str().ok_or("label digest")?;
        let points = point_lines(&checked_file(&directory.join(point_name), point_sha)?)?;
        let labels: Vec<Value> = checked_file(&directory.join(label_name), label_sha)?
            .split(|&byte| byte == b'\n')
            .filter(|line| !line.is_empty())
            .map(|line| serde_json::from_slice(line).map_err(|error| error.to_string()))
            .collect::<Result<_, _>>()?;
        if points.len() != 1024 || labels.len() != 1024 {
            return Err(format!("block {block} rows differ"));
        }
        let candidates = spec["candidates_scanned"]
            .as_u64()
            .ok_or("candidate count")?;
        let mut accepted = 0usize;
        let mut zero_rejects = 0u64;
        let mut orbit_rejects = 0u64;
        for candidate in 0..candidates {
            let input = format!("n37-four-policy-sparse16-20261004|{block}|{candidate}");
            let scalar = (BigUint::from_bytes_be(&sha256(input.as_bytes())) % BigUint::from(R))
                .to_u64()
                .ok_or("reduced scalar")?;
            if scalar == 0 {
                zero_rejects += 1;
                continue;
            }
            let general = scalar_mul(&ctx.kc.curve, ctx.kc.generator(), &BigUint::from(scalar));
            let point = ctx.control.fast.lift(&general);
            if point.infinity || !ctx.control.fast.is_on_curve(point) {
                return Err(format!("candidate {candidate} general point failed"));
            }
            if !seen.insert(key(&ctx, point)?) {
                orbit_rejects += 1;
                continue;
            }
            if accepted >= points.len()
                || labels[accepted]["candidate"].as_u64() != Some(candidate)
                || labels[accepted]["scalar"].as_u64() != Some(scalar)
                || labels[accepted]["point"] != json!([point.x, point.y])
                || points[accepted] != [point.x, point.y]
            {
                return Err(format!(
                    "block {block} accepted candidate {candidate} differs"
                ));
            }
            accepted += 1;
            verified += 1;
        }
        if accepted != 1024
            || zero_rejects != spec["zero_rejects"].as_u64().ok_or("zero rejects")?
            || orbit_rejects != spec["orbit_rejects"].as_u64().ok_or("orbit rejects")?
        {
            return Err(format!("block {block} acceptance counts differ"));
        }
    }
    Ok(json!({
        "schema":"n37-four-policy-sparse16-input-replay-v1",
        "status":"PASS","verified_points":verified,
        "prior_rows":prior_rows,"distinct_orbits_including_prior":seen.len(),
        "frozen_sha256":hex::encode(sha256(&fs::read(directory.join("FROZEN.json")).map_err(|error| error.to_string())?)),
    }))
}

fn main() {
    let directory = std::env::args()
        .nth(1)
        .expect("usage: n37_four_policy_sparse16_inputs_replay INPUT_DIR OUTPUT.json");
    let output = std::env::args().nth(2).expect("output path");
    assert!(std::env::args().nth(3).is_none(), "unexpected argument");
    let value = replay(Path::new(&directory)).expect("input replay");
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&output)
        .expect("refuse to overwrite input replay");
    file.write_all(format!("{}\n", serde_json::to_string_pretty(&value).unwrap()).as_bytes())
        .expect("write input replay");
    println!("PASS {output}");
}
