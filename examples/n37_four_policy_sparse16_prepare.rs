//! Freeze new point-only n37 targets for the sparse four-policy PDP gate.
//! The scalar labels are written separately and are never producer inputs.

#![recursion_limit = "256"]

#[path = "support/n37_policy_common.rs"]
#[allow(dead_code)]
mod common;

use common::{Context, R};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::fs::OpenOptions;
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

fn orbit_key(ctx: &Context, point: FastPoint) -> Result<u64, String> {
    if point.infinity || !ctx.control.fast.is_on_curve(point) {
        return Err("orbit key requires a finite source point".into());
    }
    let mut x = point.x;
    let mut key = x;
    for _ in 1..37 {
        x = ctx.control.gf.sqr(x);
        key = key.min(x);
    }
    if ctx.control.gf.sqr(x) != point.x {
        return Err("Frobenius closure failed".into());
    }
    Ok(key)
}

fn file_new(path: &Path, bytes: &[u8]) -> Result<String, String> {
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|error| format!("create {}: {error}", path.display()))?;
    file.write_all(bytes)
        .map_err(|error| format!("write {}: {error}", path.display()))?;
    Ok(hex::encode(sha256(bytes)))
}

fn prepare(out_dir: &Path) -> Result<Value, String> {
    let ctx = Context::load()?;
    let generator = ctx.control.fast.lift(ctx.kc.generator());
    if generator.infinity || !ctx.control.fast.mul_u64(generator, R).infinity {
        return Err("registered generator failed subgroup check".into());
    }
    let mut seen = HashSet::new();
    let mut prior_rows = 0usize;
    for (path, expected_sha) in PRIOR {
        let data = common::verified_bytes(path, expected_sha)?;
        for line in data
            .split(|&byte| byte == b'\n')
            .filter(|line| !line.is_empty())
        {
            let [x, y]: [u64; 2] = serde_json::from_slice(line)
                .map_err(|error| format!("parse prior point: {error}"))?;
            let point = FastPoint::affine(x, y);
            if !ctx.control.fast.mul_u64(point, R).infinity {
                return Err("prior point outside subgroup".into());
            }
            seen.insert(orbit_key(&ctx, point)?);
            prior_rows += 1;
        }
    }
    if prior_rows != 2048 {
        return Err("prior point files did not contain 2048 rows".into());
    }
    let mut blocks = Vec::new();
    for block in 0..2usize {
        let mut point_bytes = Vec::new();
        let mut label_bytes = Vec::new();
        let mut accepted = 0usize;
        let mut candidate = 0u64;
        let mut zero_rejects = 0u64;
        let mut orbit_rejects = 0u64;
        while accepted < 1024 {
            let input = format!("n37-four-policy-sparse16-20261004|{block}|{candidate}");
            let digest = sha256(input.as_bytes());
            let scalar = (BigUint::from_bytes_be(&digest) % BigUint::from(R))
                .to_u64()
                .ok_or("reduced scalar does not fit u64")?;
            if scalar == 0 {
                zero_rejects += 1;
                candidate += 1;
                continue;
            }
            let point = ctx.control.fast.mul_u64(generator, scalar);
            if point.infinity || !ctx.control.fast.is_on_curve(point) {
                return Err("generated point invalid".into());
            }
            if !seen.insert(orbit_key(&ctx, point)?) {
                orbit_rejects += 1;
                candidate += 1;
                continue;
            }
            point_bytes.extend_from_slice(format!("{}\n", json!([point.x, point.y])).as_bytes());
            label_bytes.extend_from_slice(
                format!(
                    "{}\n",
                    json!({"candidate":candidate,"scalar":scalar,"point":[point.x,point.y]})
                )
                .as_bytes(),
            );
            accepted += 1;
            candidate += 1;
        }
        let point_name = format!("points-b{block}.jsonl");
        let label_name = format!("labels-b{block}.jsonl");
        let point_sha = file_new(&out_dir.join(&point_name), &point_bytes)?;
        let label_sha = file_new(&out_dir.join(&label_name), &label_bytes)?;
        blocks.push(json!({
            "block":block,"rows":accepted,"candidates_scanned":candidate,
            "zero_rejects":zero_rejects,"orbit_rejects":orbit_rejects,
            "point_file":point_name,"point_sha256":point_sha,
            "label_file":label_name,"label_sha256":label_sha,
        }));
    }
    Ok(json!({
        "schema":"n37-four-policy-sparse16-input-v1",
        "status":"FROZEN","curve":"n37-koblitz-a0-b1",
        "subgroup_order":R,"generator":[generator.x,generator.y],
        "prior_point_sha256":{"b03":PRIOR[0].1,"b04":PRIOR[1].1},
        "prior_rows":prior_rows,"blocks":blocks,
    }))
}

fn main() {
    let directory = std::env::args()
        .nth(1)
        .expect("usage: n37_four_policy_sparse16_prepare OUTPUT_DIR");
    assert!(std::env::args().nth(2).is_none(), "unexpected argument");
    let out_dir = Path::new(&directory);
    std::fs::create_dir_all(out_dir).expect("create output directory");
    let freeze = prepare(out_dir).expect("freeze inputs");
    let bytes = format!("{}\n", serde_json::to_string_pretty(&freeze).unwrap());
    file_new(&out_dir.join("FROZEN.json"), bytes.as_bytes()).expect("save freeze");
    println!("FROZEN {}", out_dir.display());
}
