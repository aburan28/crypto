//! Independently rederive the frozen n37 Callgrind census from its raw parts.
//! Usage: ecbench_callgrind_analyze ARTIFACT_DIR JOBS.json RECORDS.jsonl OUT.json

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use crypto_lib::cryptanalysis::ecbench::{callgrind, canonical};
use serde::Deserialize;
use serde_json::{json, Value};

const ARCHIVED_RECORDS_SHA256: &str =
    "82c8b31640074b2f829d88e546bf487ee4b598edf4ef26f779984d625ed7de3c";
const R: f64 = 230_603_167.0;
const BOOTSTRAP_SAMPLES: usize = 20_000;
const BOOTSTRAP_SEED: u64 = 20_261_004;
const ARMS: [&str; 4] = ["ic-k16", "ic-k42", "ic-k42-control", "rho-strong"];

#[derive(Deserialize)]
struct Job {
    seq: u64,
    arm: String,
    workload_id: String,
    algorithm_seed: String,
    expected_recovered: String,
    archived_status: String,
}

#[derive(Deserialize)]
struct Row {
    seq: u64,
    arm: String,
    workload_id: String,
    algorithm_seed: String,
    expected_recovered: String,
    status: String,
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
            "usage: ecbench_callgrind_analyze ARTIFACT_DIR JOBS.json RECORDS.jsonl OUT.json".into(),
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
    let jobs: Vec<Job> = parse_json(&jobs_path)?;
    let census_path = artifact.join("CENSUS.json");
    let census: Census = parse_json(&census_path)?;
    if census.schema != "ecbench.callgrind_census/v1"
        || census.unit != "callgrind.Ir"
        || jobs.len() != 32
        || census.rows.len() != 32
    {
        return Err("wrong census schema, unit, or frozen job count".into());
    }
    let jobs_by_seq: BTreeMap<u64, &Job> = jobs.iter().map(|j| (j.seq, j)).collect();
    if jobs_by_seq.len() != jobs.len() {
        return Err("duplicate job seq".into());
    }
    let mut groups: BTreeMap<String, BTreeMap<String, u64>> = BTreeMap::new();
    let mut seen = BTreeSet::new();
    for row in &census.rows {
        let job = jobs_by_seq
            .get(&row.seq)
            .ok_or_else(|| format!("unplanned seq {}", row.seq))?;
        if !seen.insert(row.seq)
            || row.arm != job.arm
            || row.workload_id != job.workload_id
            || row.algorithm_seed != job.algorithm_seed
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
            || field(archived_row, "arm")? != job.arm
            || field(&archived_row["workload"], "workload_id")? != job.workload_id
            || field(archived_row, "algorithm_seed")? != job.algorithm_seed
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
            || field(&input, "expected_method_id")? != field(&archived_row["method"], "method_id")?
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
        let raw = callgrind::solve_ir(&dir.join("profile"))?;
        if row.solve_ir != Some(raw.solve_ir)
            || ir["solve_ir"].as_u64() != Some(raw.solve_ir)
            || ir["solve_parts"].as_u64() != Some(raw.solve_parts as u64)
        {
            return Err(format!("seq {} raw Callgrind parts disagree", row.seq));
        }
        let old = groups
            .entry(job.workload_id.clone())
            .or_default()
            .insert(job.arm.clone(), raw.solve_ir);
        if old.is_some() {
            return Err(format!("duplicate arm {} on {}", job.arm, job.workload_id));
        }
    }
    if seen.len() != 32 || groups.len() != 8 {
        return Err("missing job or workload".into());
    }
    let groups: Vec<BTreeMap<String, u64>> = groups.into_values().collect();
    if groups
        .iter()
        .any(|group| ARMS.iter().any(|arm| !group.contains_key(*arm)) || group.len() != 4)
    {
        return Err("a target has missing or unexpected arms".into());
    }
    let summaries: Vec<Value> = ARMS
        .iter()
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
        .collect();
    let k16_k42 = comparison(&groups, "ic-k16", "ic-k42");
    let k16_rho = comparison(&groups, "ic-k16", "rho-strong");
    let aa = comparison(&groups, "ic-k42", "ic-k42-control");
    let max_aa_deviation = groups
        .iter()
        .map(|group| (group["ic-k42"] as f64 / group["ic-k42-control"] as f64 - 1.0).abs())
        .fold(0.0, f64::max);
    let retain_k16 = max_aa_deviation <= 0.02
        && k16_k42["bootstrap_95"][1]
            .as_f64()
            .ok_or("bootstrap high missing")?
            < 1.0;
    let result = json!({
        "schema": "ecbench.callgrind_decision/v1",
        "unit": "callgrind.Ir",
        "normalised_unit": "crypto.S.callgrind_ir",
        "subgroup_order": R as u64,
        "workloads": groups.len(),
        "profiles": seen.len(),
        "archived_records_sha256": ARCHIVED_RECORDS_SHA256,
        "jobs_sha256": canonical::sha256_hex(&read(&jobs_path)?),
        "census_sha256": canonical::sha256_hex(&read(&census_path)?),
        "provenance_sha256": canonical::sha256_hex(&read(&artifact.join("PROVENANCE.txt"))?),
        "all_archived_and_profiled_scalars_verified": true,
        "arms": summaries,
        "comparisons": [k16_k42, k16_rho, aa],
        "aa_max_absolute_relative_deviation": max_aa_deviation,
        "decision": if retain_k16 { "retain_k16_instruction_lead" }
                    else { "stop_column_lead_pending_cost_or_aa_investigation" },
        "claim_limit": "same-binary Callgrind user-space instruction diagnostic; no isolated native wall speedup or cross-unit floor ratio",
    });
    let bytes = serde_json::to_vec_pretty(&result).map_err(|e| e.to_string())?;
    let mut bytes = bytes;
    bytes.push(b'\n');
    std::fs::write(&output, bytes).map_err(|e| format!("{}: {e}", output.display()))?;
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("ecbench_callgrind_analyze: {e}");
        std::process::exit(1);
    }
}
