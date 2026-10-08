//! Check the preregistered ecbench public point against the frozen n37 orbit
//! exclusions without computing a logarithm or invoking either solver.

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ic_boundary::{koblitz_instance, BinaryInstance};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use serde_json::{json, Value};
use std::collections::{BTreeSet, HashSet};
use std::fs;
use std::path::Path;
use std::process::Command;

const PRIOR: &str = "research/notes/ecc2k130/n37_shared_rank_targets_20261003";
const SOURCE: &str = "research/notes/ecc2k130/n37_four_policy_support_20261003/SOURCE42.jsonl";
const SOURCE_SHA: &str = "0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75";
const RANK: &str = "research/notes/ecc2k130/shared_rank_folded_table_20261003/RAW.json";
const RANK_SHA: &str = "293d273c071c280326825c5a4f6d12bce00e46349735e5475b16bcd40264dd3b";

fn hex(s: &str) -> Result<u64, String> {
    u64::from_str_radix(s.trim_start_matches("0x"), 16).map_err(|e| e.to_string())
}

fn point(v: &Value) -> Result<FastPoint, String> {
    let [x, y]: [u64; 2] = serde_json::from_value(v.clone()).map_err(|e| e.to_string())?;
    Ok(FastPoint::affine(x, y))
}

fn orbit_key(inst: &BinaryInstance, p: FastPoint) -> Result<u64, String> {
    if p.infinity || !inst.fast.is_on_curve(p) {
        return Err("invalid point in orbit inventory".into());
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
    let bytes = fs::read(path).map_err(|e| format!("{path}: {e}"))?;
    if sha256_hex(&bytes) != digest {
        return Err(format!("{path}: SHA-256 differs from frozen inventory"));
    }
    Ok(bytes)
}

fn add_points(
    inst: &BinaryInstance,
    excluded: &mut HashSet<u64>,
    bytes: &[u8],
) -> Result<usize, String> {
    let mut rows = 0;
    for line in bytes.split(|&b| b == b'\n').filter(|line| !line.is_empty()) {
        rows += 1;
        let value: Value = serde_json::from_slice(line).map_err(|e| e.to_string())?;
        excluded.insert(orbit_key(inst, point(&value)?)?);
    }
    Ok(rows)
}

fn run(plan_path: &Path, out_path: &Path) -> Result<(), String> {
    let plan_bytes = fs::read(plan_path).map_err(|e| e.to_string())?;
    let plan: Value = serde_json::from_slice(&plan_bytes).map_err(|e| e.to_string())?;
    let workloads = plan["workloads"]
        .as_array()
        .ok_or("plan has no workloads")?;
    if workloads.len() != 1 || plan["spec_id"] != "ECS1h7421a101192f" {
        return Err("plan is not the frozen one-target spec".into());
    }
    let w = &workloads[0];
    if w["target_law"] != "hash_to_subgroup_v1"
        || w["target_seed"] != 202610041751u64
        || w["target_index"] != 0
        || !w["planted"].is_null()
    {
        return Err("plan does not contain the frozen public-point law".into());
    }
    let coordinates = w["target"].as_array().ok_or("target not an array")?;
    let q = FastPoint::affine(
        hex(coordinates[0].as_str().ok_or("target x not hex")?)?,
        hex(coordinates[1].as_str().ok_or("target y not hex")?)?,
    );
    let inst = koblitz_instance(0, 37).ok_or("n37 Koblitz curve unavailable")?;
    if inst.name != "icv1-f2m37-tm534059-32aad96b" || inst.r != 230603167 {
        return Err("curve differs from preregistration".into());
    }
    let q_key = orbit_key(&inst, q)?;
    if !inst.fast.mul_u64(q, inst.r).infinity {
        return Err("public target is not in the specified subgroup".into());
    }
    let fixture: Value = serde_json::from_slice(
        &fs::read(format!("{PRIOR}/FIXTURE.json")).map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;
    let inventory = fixture["inventory"]
        .as_array()
        .ok_or("prior inventory absent")?;
    let mut paths = BTreeSet::new();
    let mut excluded = HashSet::new();
    let mut point_rows = 0;
    for entry in inventory {
        let path = entry["path"].as_str().ok_or("inventory path absent")?;
        let sha = entry["sha256"].as_str().ok_or("inventory hash absent")?;
        let bytes = pinned(path, sha)?;
        paths.insert(path.to_owned());
        let rows = if Path::new(path)
            .file_name()
            .and_then(|x| x.to_str())
            .is_some_and(|x| x.starts_with("n37_"))
        {
            add_points(&inst, &mut excluded, &bytes)?
        } else {
            bytes
                .split(|&b| b == b'\n')
                .filter(|l| !l.is_empty())
                .count()
        };
        if entry["rows"].as_u64() != Some(rows as u64) {
            return Err(format!("{path}: row count differs from prior inventory"));
        }
        point_rows += rows;
    }
    let prior_points = format!("{PRIOR}/TARGETS.points.jsonl");
    let bytes = pinned(
        &prior_points,
        fixture["points_sha256"].as_str().ok_or("point hash")?,
    )?;
    let prior_target_rows = add_points(&inst, &mut excluded, &bytes)?;
    paths.insert(prior_points);
    let tracked = Command::new("git")
        .args(["ls-files", "research/notes/ecc2k130"])
        .output()
        .map_err(|e| e.to_string())?;
    if !tracked.status.success() {
        return Err("git ls-files failed".into());
    }
    let tracked_paths: BTreeSet<String> = String::from_utf8(tracked.stdout)
        .map_err(|e| e.to_string())?
        .lines()
        .filter(|x| x.ends_with(".points.jsonl"))
        .map(str::to_owned)
        .collect();
    if tracked_paths != paths {
        return Err("tracked point-file set differs from frozen exclusions".into());
    }
    let source: Value = serde_json::from_slice(
        pinned(SOURCE, SOURCE_SHA)?
            .split(|&b| b == b'\n')
            .next()
            .ok_or("empty source")?,
    )
    .map_err(|e| e.to_string())?;
    let base = source["factor_base_point_coordinates"]
        .as_array()
        .ok_or("source base absent")?;
    for value in base {
        excluded.insert(orbit_key(&inst, point(value)?)?);
    }
    let rank: Value =
        serde_json::from_slice(&pinned(RANK, RANK_SHA)?).map_err(|e| e.to_string())?;
    let probes = rank["relations"].as_array().ok_or("rank probes absent")?;
    for probe in probes {
        excluded.insert(orbit_key(
            &inst,
            FastPoint::affine(
                probe["point_x"].as_u64().ok_or("probe x")?,
                probe["point_y"].as_u64().ok_or("probe y")?,
            ),
        )?);
    }
    let mut sorted: Vec<u64> = excluded.iter().copied().collect();
    sorted.sort_unstable();
    let mut orbit_bytes = b"n37-excluded-orbits-v1\0".to_vec();
    for x in &sorted {
        orbit_bytes.extend_from_slice(&x.to_le_bytes());
    }
    let fresh = !excluded.contains(&q_key);
    let result = json!({
        "schema": "n37-shared-rank-ecbench-input-check/v1",
        "spec_id": plan["spec_id"],
        "workload_id": w["workload_id"],
        "plan_sha256": sha256_hex(&plan_bytes),
        "target": [q.x, q.y],
        "target_orbit_key": q_key,
        "subgroup_checked": true,
        "fresh_orbit": fresh,
        "inventory_files": paths.len(),
        "inventory_rows": point_rows + prior_target_rows,
        "base_points": base.len(),
        "rank_probes": probes.len(),
        "excluded_orbits": sorted.len(),
        "excluded_orbits_sha256": sha256_hex(&orbit_bytes),
        "source_sha256": SOURCE_SHA,
        "rank_sha256": RANK_SHA,
    });
    serde_json::to_writer_pretty(
        fs::File::create(out_path).map_err(|e| e.to_string())?,
        &result,
    )
    .map_err(|e| e.to_string())?;
    if !fresh {
        return Err("preregistered point collides with an excluded orbit".into());
    }
    Ok(())
}

fn main() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let plan = args
        .next()
        .ok_or("usage: input_check PLAN_JSON OUTPUT_JSON")?;
    let out = args.next().ok_or("missing output path")?;
    if args.next().is_some() {
        return Err("too many arguments".into());
    }
    run(Path::new(&plan), Path::new(&out))
}
