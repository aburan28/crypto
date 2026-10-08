//! Independently rederive the n37 target-only instruction gate from raw profiles.
//! Usage: ecbench_n37_online_ir_analyze ARTIFACT_DIR JOBS.json RECORDS.jsonl OUT.json

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use crypto_lib::cryptanalysis::ecbench::{callgrind, canonical};
use serde::Deserialize;
use serde_json::{json, Value};

const ARCHIVED_RECORDS_SHA256: &str =
    "60838cb661b6df9ba801cca089975ce0d6bca259e1bb18f148461340fbe0f169";
const FROZEN_JOBS_SHA256: &str = "a2d63e9063c98955e7fdf49b137d1c9baafc9912b903f8aa79cd69a43035a32f";
const SESSION_ID: &str = "ECBS1he4e70185e6e2";
const ORIGINAL_ENV_CLASS: &str = "ECBENV2hd82681268e96";
const R: f64 = 230_603_167.0;
const BOOTSTRAP_SAMPLES: usize = 20_000;
const BOOTSTRAP_SEED: u64 = 202_610_044_184;
const ARMS: [&str; 4] = ["ic-k8", "ic-k16", "ic-k16-control", "rho-strong"];

#[derive(Deserialize)]
struct Job {
    seq: u64,
    arm: String,
    workload_id: String,
    algorithm_seed: String,
    expected_method_id: String,
    expected_recovered: String,
    archived_status: String,
}

#[derive(Deserialize)]
struct Row {
    seq: u64,
    arm: String,
    workload_id: String,
    algorithm_seed: String,
    expected_method_id: String,
    expected_recovered: String,
    status: String,
    pre_online_ir: Option<u64>,
    online_ir: Option<u64>,
    post_online_ir: Option<u64>,
    solve_ir: Option<u64>,
}

#[derive(Deserialize)]
struct Census {
    schema: String,
    unit: String,
    rows: Vec<Row>,
}

fn read(path: &Path) -> Result<Vec<u8>, String> {
    std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))
}

fn parse_json<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T, String> {
    serde_json::from_slice(&read(path)?).map_err(|e| format!("{}: {e}", path.display()))
}

fn field<'a>(value: &'a Value, key: &str) -> Result<&'a str, String> {
    value[key]
        .as_str()
        .ok_or_else(|| format!("missing string {key}"))
}

fn check_sums(dir: &Path) -> Result<(), String> {
    let sums = std::fs::read_to_string(dir.join("SHA256SUMS"))
        .map_err(|e| format!("{}: {e}", dir.display()))?;
    let mut names = BTreeSet::new();
    for line in sums.lines() {
        let mut fields = line.split_whitespace();
        let expected = fields.next().ok_or("empty SHA256SUMS row")?;
        let named = fields.next().ok_or("SHA256SUMS row has no path")?;
        let name = Path::new(named)
            .file_name()
            .and_then(|s| s.to_str())
            .ok_or("SHA256SUMS path has no file name")?;
        if !names.insert(name.to_string()) {
            return Err(format!("duplicate SHA256SUMS entry {name}"));
        }
        let actual = canonical::sha256_hex(&read(&dir.join(name))?);
        if actual != expected {
            return Err(format!("{name}: SHA-256 differs from the CI receipt"));
        }
    }
    for needed in ["input.json", "child.json", "ir.json", "valgrind.log"] {
        if !names.contains(needed) {
            return Err(format!("SHA256SUMS misses {needed}"));
        }
    }
    Ok(())
}

fn next(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

fn comparison(groups: &[BTreeMap<String, u64>], a: &str, b: &str) -> Value {
    let pairs: Vec<(f64, f64)> = groups
        .iter()
        .map(|group| (group[a] as f64, group[b] as f64))
        .collect();
    let ratio =
        pairs.iter().map(|pair| pair.0).sum::<f64>() / pairs.iter().map(|pair| pair.1).sum::<f64>();
    let ratios: Vec<f64> = pairs.iter().map(|pair| pair.0 / pair.1).collect();
    let mut state = BOOTSTRAP_SEED;
    let mut samples = Vec::with_capacity(BOOTSTRAP_SAMPLES);
    for _ in 0..BOOTSTRAP_SAMPLES {
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for _ in 0..pairs.len() {
            let index = (next(&mut state) % pairs.len() as u64) as usize;
            numerator += pairs[index].0;
            denominator += pairs[index].1;
        }
        samples.push(numerator / denominator);
    }
    samples.sort_by(f64::total_cmp);
    json!({
        "numerator": a,
        "denominator": b,
        "ratio_of_sums": ratio,
        "paired_ratio_range": [ratios.iter().copied().fold(f64::INFINITY, f64::min),
                               ratios.iter().copied().fold(f64::NEG_INFINITY, f64::max)],
        "bootstrap_95": [samples[499], samples[19_499]],
        "bootstrap_samples": BOOTSTRAP_SAMPLES,
        "bootstrap_seed": BOOTSTRAP_SEED,
    })
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        return Err(
            "usage: ecbench_n37_online_ir_analyze ARTIFACT_DIR JOBS.json RECORDS.jsonl OUT.json"
                .into(),
        );
    }
    let (artifact, jobs_path, records_path, output) = (
        PathBuf::from(&args[1]),
        PathBuf::from(&args[2]),
        PathBuf::from(&args[3]),
        PathBuf::from(&args[4]),
    );
    let records_bytes = read(&records_path)?;
    if canonical::sha256_hex(&records_bytes) != ARCHIVED_RECORDS_SHA256 {
        return Err("archived records SHA-256 differs from the frozen input".into());
    }
    let archived: BTreeMap<u64, Value> = std::str::from_utf8(&records_bytes)
        .map_err(|e| e.to_string())?
        .lines()
        .map(|line| {
            let record: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
            let seq = record["seq"].as_u64().ok_or("archived row misses seq")?;
            Ok((seq, record))
        })
        .collect::<Result<_, String>>()?;
    if archived.len() != 384 {
        return Err("archived session is not the frozen 384-record session".into());
    }
    if canonical::sha256_hex(&read(&jobs_path)?) != FROZEN_JOBS_SHA256 {
        return Err("frozen JOBS.json SHA-256 differs".into());
    }
    let jobs: Vec<Job> = parse_json(&jobs_path)?;
    let independent: Value = parse_json(&artifact.join("INDEPENDENT-AUDIT.json"))?;
    if independent["ok"] != true
        || field(&independent, "session_id")? != SESSION_ID
        || independent["records"].as_u64() != Some(384)
        || independent["verified_records"].as_u64() != Some(384)
        || independent["replays"].as_array().map(Vec::len) != Some(320)
        || independent["replays"]
            .as_array()
            .ok_or("independent replay list missing")?
            .iter()
            .any(|replay| replay["reproduced"] != true || replay["detail"] != "identical")
        || field(&independent, "auditor_env_class_id")? == ORIGINAL_ENV_CLASS
        || field(&independent["files"], "records.jsonl")? != ARCHIVED_RECORDS_SHA256
    {
        return Err("independent Linux receipt does not reproduce the frozen session".into());
    }
    let census_path = artifact.join("CENSUS.json");
    let census: Census = parse_json(&census_path)?;
    if census.schema != "ecbench.callgrind_online_census/v1"
        || census.unit != "callgrind.Ir"
        || jobs.len() != 64
        || census.rows.len() != 64
    {
        return Err("wrong census schema, unit, or frozen job count".into());
    }
    let jobs_by_seq: BTreeMap<u64, &Job> = jobs.iter().map(|j| (j.seq, j)).collect();
    if jobs_by_seq.len() != jobs.len() {
        return Err("duplicate job seq".into());
    }
    if jobs_by_seq.keys().copied().collect::<Vec<_>>() != (64..128).collect::<Vec<_>>() {
        return Err("frozen jobs do not cover the first measured round".into());
    }
    let mut online_groups: BTreeMap<String, BTreeMap<String, u64>> = BTreeMap::new();
    let mut solve_groups: BTreeMap<String, BTreeMap<String, u64>> = BTreeMap::new();
    let mut per_target_parts: BTreeMap<String, BTreeMap<String, Value>> = BTreeMap::new();
    let mut seen = BTreeSet::new();
    for row in &census.rows {
        let job = jobs_by_seq
            .get(&row.seq)
            .ok_or_else(|| format!("unplanned seq {}", row.seq))?;
        if !seen.insert(row.seq)
            || row.arm != job.arm
            || row.workload_id != job.workload_id
            || row.algorithm_seed != job.algorithm_seed
            || row.expected_method_id != job.expected_method_id
            || row.expected_recovered != job.expected_recovered
            || row.status != "verified"
            || job.archived_status != "verified"
        {
            return Err(format!("seq {} differs from the frozen job", row.seq));
        }
        let archived_row = archived
            .get(&row.seq)
            .ok_or_else(|| format!("archived seq {} absent", row.seq))?;
        if archived_row["round"] != 1
            || archived_row["warmup"] != false
            || field(archived_row, "arm")? != job.arm
            || field(&archived_row["workload"], "workload_id")? != job.workload_id
            || field(archived_row, "algorithm_seed")? != job.algorithm_seed
            || field(&archived_row["method"], "method_id")? != job.expected_method_id
            || field(&archived_row["outcome"], "recovered")? != job.expected_recovered
            || field(&archived_row["outcome"], "status")? != "verified"
        {
            return Err(format!("seq {} differs from the archived run", row.seq));
        }
        let dir = artifact.join(format!("seq-{}", row.seq));
        check_sums(&dir)?;
        let input: Value = parse_json(&dir.join("input.json"))?;
        let child: Value = parse_json(&dir.join("child.json"))?;
        let ir: Value = parse_json(&dir.join("ir.json"))?;
        if field(&input, "expected_workload_id")? != job.workload_id
            || field(&input, "expected_method_id")? != job.expected_method_id
            || input["algorithm_seed"].as_u64().map(|v| v.to_string())
                != Some(job.algorithm_seed.clone())
            || field(&child, "workload_id")? != job.workload_id
            || child["error"] != Value::Null
            || child["report"]["exhausted"] != false
            || child["report"]["recovered"].as_u64().map(|v| v.to_string())
                != Some(job.expected_recovered.clone())
            || child["method_id"] != input["expected_method_id"]
        {
            return Err(format!("seq {} child failed exact scalar replay", row.seq));
        }
        let raw = callgrind::online_ir(&dir.join("profile"))?;
        if row.pre_online_ir != Some(raw.pre_online_ir)
            || row.online_ir != Some(raw.online_ir)
            || row.post_online_ir != Some(raw.post_online_ir)
            || row.solve_ir != Some(raw.solve_ir)
            || ir["pre_online_ir"].as_u64() != Some(raw.pre_online_ir)
            || ir["online_ir"].as_u64() != Some(raw.online_ir)
            || ir["post_online_ir"].as_u64() != Some(raw.post_online_ir)
            || ir["solve_ir"].as_u64() != Some(raw.solve_ir)
            || ir["schema"] != "ecbench.callgrind_online_ir/v1"
            || raw
                .pre_online_ir
                .checked_add(raw.online_ir)
                .and_then(|n| n.checked_add(raw.post_online_ir))
                != Some(raw.solve_ir)
        {
            return Err(format!("seq {} raw Callgrind parts disagree", row.seq));
        }
        let old = online_groups
            .entry(job.workload_id.clone())
            .or_default()
            .insert(job.arm.clone(), raw.online_ir);
        if old.is_some() {
            return Err(format!("duplicate arm {} on {}", job.arm, job.workload_id));
        }
        solve_groups
            .entry(job.workload_id.clone())
            .or_default()
            .insert(job.arm.clone(), raw.solve_ir);
        per_target_parts
            .entry(job.workload_id.clone())
            .or_default()
            .insert(
                job.arm.clone(),
                json!({"pre_online_ir": raw.pre_online_ir, "online_ir": raw.online_ir,
                       "post_online_ir": raw.post_online_ir, "solve_ir": raw.solve_ir}),
            );
    }
    if seen.len() != 64 || online_groups.len() != 16 || solve_groups.len() != 16 {
        return Err("missing job or workload".into());
    }
    let per_target: Vec<Value> = per_target_parts
        .iter()
        .map(|(workload_id, counts)| json!({"workload_id": workload_id, "parts": counts}))
        .collect();
    let online_groups: Vec<BTreeMap<String, u64>> = online_groups.into_values().collect();
    let solve_groups: Vec<BTreeMap<String, u64>> = solve_groups.into_values().collect();
    if online_groups
        .iter()
        .chain(&solve_groups)
        .any(|group| ARMS.iter().any(|arm| !group.contains_key(*arm)) || group.len() != 4)
    {
        return Err("a target has missing or unexpected arms".into());
    }
    let summaries = |groups: &[BTreeMap<String, u64>]| -> Vec<Value> {
        ARMS.iter()
            .map(|arm| {
                let counts: Vec<u64> = groups.iter().map(|group| group[*arm]).collect();
                let mean = counts.iter().map(|&n| n as f64).sum::<f64>() / counts.len() as f64;
                json!({
                    "arm": arm,
                    "runs": counts.len(),
                    "mean_ir": mean,
                    "s_ir_per_target": mean / R.sqrt(),
                    "min_ir": counts.iter().min(),
                    "max_ir": counts.iter().max(),
                })
            })
            .collect()
    };
    let comparisons = |groups: &[BTreeMap<String, u64>]| -> Vec<Value> {
        vec![
            comparison(groups, "ic-k8", "ic-k16"),
            comparison(groups, "ic-k8", "rho-strong"),
            comparison(groups, "ic-k16", "rho-strong"),
            comparison(groups, "ic-k16", "ic-k16-control"),
        ]
    };
    let online_comparisons = comparisons(&online_groups);
    let solve_comparisons = comparisons(&solve_groups);
    let max_aa_deviation = online_groups
        .iter()
        .map(|group| (group["ic-k16"] as f64 / group["ic-k16-control"] as f64 - 1.0).abs())
        .fold(0.0, f64::max);
    let solve_max_aa_deviation = solve_groups
        .iter()
        .map(|group| (group["ic-k16"] as f64 / group["ic-k16-control"] as f64 - 1.0).abs())
        .fold(0.0, f64::max);
    let decision = if max_aa_deviation > 0.02 {
        "no_selection_aa_gate_failed"
    } else if online_comparisons[0]["bootstrap_95"][0]
        .as_f64()
        .ok_or("bootstrap low missing")?
        > 1.10
    {
        "prioritize_k16_for_isolated_n37_online_wall_gate"
    } else if online_comparisons[0]["bootstrap_95"][1]
        .as_f64()
        .ok_or("bootstrap high missing")?
        < 0.90
    {
        "prioritize_k8_for_isolated_n37_online_wall_gate"
    } else {
        "retain_both_for_isolated_n37_online_wall_gate"
    };
    let result = json!({
        "schema": "ecbench.n37_online_ir_decision/v1",
        "unit": "callgrind.Ir",
        "normalised_unit": "crypto.S.callgrind_ir",
        "subgroup_order": R as u64,
        "workloads": online_groups.len(),
        "profiles": seen.len(),
        "archived_records_sha256": ARCHIVED_RECORDS_SHA256,
        "jobs_sha256": canonical::sha256_hex(&read(&jobs_path)?),
        "census_sha256": canonical::sha256_hex(&read(&census_path)?),
        "provenance_sha256": canonical::sha256_hex(&read(&artifact.join("PROVENANCE.txt"))?),
        "independent_audit_sha256": canonical::sha256_hex(&read(&artifact.join("INDEPENDENT-AUDIT.json"))?),
        "independent_auditor_env_class_id": field(&independent, "auditor_env_class_id")?,
        "all_archived_and_profiled_scalars_verified": true,
        "online": {"arms": summaries(&online_groups), "comparisons": online_comparisons,
                   "aa_max_absolute_relative_deviation": max_aa_deviation},
        "complete_solve": {"arms": summaries(&solve_groups), "comparisons": solve_comparisons,
                          "aa_max_absolute_relative_deviation": solve_max_aa_deviation},
        "per_target": per_target,
        "decision": decision,
        "claim_limit": "same-binary Callgrind target-only and complete-solve user-space instruction diagnostics; no isolated native online wall speedup or n41/n53 transfer",
    });
    let bytes = serde_json::to_vec_pretty(&result).map_err(|e| e.to_string())?;
    let mut bytes = bytes;
    bytes.push(b'\n');
    std::fs::write(&output, bytes).map_err(|e| format!("{}: {e}", output.display()))?;
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("ecbench_n37_online_ir_analyze: {e}");
        std::process::exit(1);
    }
}
