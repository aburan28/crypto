//! Freeze point-only n37 targets independently of the target solver.
//! Protocol: research/notes/ecc2k130/n37_shared_rank_targets_20261003/PROTOCOL.md

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ic_boundary::{koblitz_instance, BinaryInstance};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde::Serialize;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::process::Command;

const SOURCE: &str = "research/notes/ecc2k130/n37_four_policy_support_20261003/SOURCE42.jsonl";
const SOURCE_SHA: &str = "0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75";
const RANK: &str = "research/notes/ecc2k130/shared_rank_folded_table_20261003/RAW.json";
const RANK_SHA: &str = "293d273c071c280326825c5a4f6d12bce00e46349735e5475b16bcd40264dd3b";
const DOMAIN: &str = "n37-shared-rank-target-v1";
const R: u64 = 230_603_167;
const COUNT: usize = 16;

#[derive(Serialize)]
struct InventoryEntry {
    path: String,
    sha256: String,
    rows: usize,
}

fn point(value: &Value) -> Result<FastPoint, String> {
    let [x, y]: [u64; 2] = serde_json::from_value(value.clone())
        .map_err(|error| format!("point coordinates: {error}"))?;
    Ok(FastPoint::affine(x, y))
}

fn key(inst: &BinaryInstance, p: FastPoint) -> Result<u64, String> {
    if p.infinity || !inst.fast.is_on_curve(p) {
        return Err("invalid n37 point in exclusion inventory".into());
    }
    let mut x = p.x;
    let mut minimum = x;
    for _ in 1..inst.n {
        x = inst.gf.sqr(x);
        minimum = minimum.min(x);
    }
    Ok(minimum)
}

fn pinned(path: &str, digest: &str) -> Result<Vec<u8>, String> {
    let bytes = fs::read(path).map_err(|error| format!("read {path}: {error}"))?;
    if sha256_hex(&bytes) != digest {
        return Err(format!("{path}: SHA-256 differs from pinned source"));
    }
    Ok(bytes)
}

fn run(points_path: &Path, fixture_path: &Path) -> Result<(), String> {
    let inst = koblitz_instance(0, 37).ok_or("n37 curve unavailable")?;
    if inst.name != "icv1-f2m37-tm534059-32aad96b" || inst.r != R {
        return Err("curve identity differs from protocol".into());
    }
    let output = Command::new("git")
        .args(["ls-files", "research/notes/ecc2k130"])
        .output()
        .map_err(|error| format!("git ls-files: {error}"))?;
    if !output.status.success() {
        return Err("git ls-files failed".into());
    }
    let tracked = String::from_utf8(output.stdout).map_err(|error| error.to_string())?;
    let mut paths: Vec<_> = tracked
        .lines()
        .filter(|path| path.ends_with(".points.jsonl"))
        .collect();
    paths.sort_unstable();
    if paths.is_empty() {
        return Err("no tracked point-only inventory".into());
    }
    let mut excluded = HashSet::new();
    let mut inventory = Vec::with_capacity(paths.len());
    for path in paths {
        let bytes = fs::read(path).map_err(|error| format!("read inventory {path}: {error}"))?;
        let is_n37 = Path::new(path)
            .file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.starts_with("n37_"));
        let mut rows = 0usize;
        for line in bytes
            .split(|&byte| byte == b'\n')
            .filter(|line| !line.is_empty())
        {
            rows += 1;
            if is_n37 {
                let value: Value = serde_json::from_slice(line)
                    .map_err(|error| format!("parse {path}:{rows}: {error}"))?;
                excluded.insert(key(&inst, point(&value)?)?);
            }
        }
        inventory.push(InventoryEntry {
            path: path.to_owned(),
            sha256: sha256_hex(&bytes),
            rows,
        });
    }
    let source = pinned(SOURCE, SOURCE_SHA)?;
    let header: Value = serde_json::from_slice(source.split(|&b| b == b'\n').next().unwrap())
        .map_err(|error| format!("source header: {error}"))?;
    let base = header["factor_base_point_coordinates"]
        .as_array()
        .ok_or("source base points absent")?;
    if base.len() != 3_108 {
        return Err("source base point count differs".into());
    }
    for value in base {
        excluded.insert(key(&inst, point(value)?)?);
    }
    let rank: Value = serde_json::from_slice(&pinned(RANK, RANK_SHA)?)
        .map_err(|error| format!("rank record: {error}"))?;
    let probes = rank["relations"].as_array().ok_or("rank probes absent")?;
    if probes.len() != 55 {
        return Err("rank probe count differs".into());
    }
    for probe in probes {
        excluded.insert(key(
            &inst,
            FastPoint::affine(
                probe["point_x"].as_u64().ok_or("probe x")?,
                probe["point_y"].as_u64().ok_or("probe y")?,
            ),
        )?);
    }
    let excluded_before_selection = excluded.len();
    let inventory_sha256 =
        sha256_hex(&serde_json::to_vec(&json!(&inventory)).map_err(|error| error.to_string())?);
    let mut targets = Vec::with_capacity(COUNT);
    let mut zero_rejections = 0u64;
    let mut orbit_rejections = 0u64;
    let mut j = 0u64;
    while targets.len() < COUNT {
        let digest = sha256(format!("{DOMAIN}|{j}").as_bytes());
        let scalar = (BigUint::from_bytes_be(&digest) % BigUint::from(R))
            .to_u64()
            .ok_or("candidate scalar outside u64")?;
        j += 1;
        if scalar == 0 {
            zero_rejections += 1;
            continue;
        }
        let q = inst.fast.mul_u64(inst.generator, scalar);
        if !excluded.insert(key(&inst, q)?) {
            orbit_rejections += 1;
            continue;
        }
        targets.push(json!({
            "candidate_index":j - 1,
            "candidate_sha256":hex::encode(digest),
            "scalar":scalar,
            "point":[q.x,q.y],
        }));
    }
    let mut file = fs::File::create(points_path).map_err(|error| error.to_string())?;
    for target in &targets {
        writeln!(file, "{}", target["point"]).map_err(|error| error.to_string())?;
    }
    let point_bytes = fs::read(points_path).map_err(|error| error.to_string())?;
    let fixture = json!({
        "schema":"n37-shared-rank-target-fixture-v1",
        "curve":inst.name,
        "subgroup_order":R,
        "generator":[inst.generator.x,inst.generator.y],
        "domain":DOMAIN,
        "target_count":COUNT,
        "candidate_count":j,
        "zero_rejections":zero_rejections,
        "orbit_rejections":orbit_rejections,
        "excluded_orbits_before_selection":excluded_before_selection,
        "inventory_sha256":inventory_sha256,
        "inventory":inventory,
        "source_sha256":SOURCE_SHA,
        "rank_sha256":RANK_SHA,
        "points_sha256":sha256_hex(&point_bytes),
        "targets":targets,
    });
    serde_json::to_writer_pretty(
        fs::File::create(fixture_path).map_err(|error| error.to_string())?,
        &fixture,
    )
    .map_err(|error| error.to_string())
}

fn main() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let points = args
        .next()
        .ok_or("usage: n37_shared_rank_target_fixture POINTS FIXTURE")?;
    let fixture = args.next().ok_or("missing FIXTURE output")?;
    if args.next().is_some() {
        return Err("too many arguments".into());
    }
    run(Path::new(&points), Path::new(&fixture))
}
