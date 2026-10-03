//! Audit and summarize the frozen n37 b02 one-target native IC/rho panel.
//! Reads only recorded results and independent replay receipts; no solving.
//! `cargo run --release --example n37_native_one_target_analyze -- SUMMARY.json`

use crypto_lib::hash::sha256::sha256;
use serde_json::{json, Value};
use std::fs;
use std::path::Path;

const NOTE: &str = "research/notes/ecc2k130/n37_native_one_target_20261003";

fn read(path: &Path) -> Result<Vec<u8>, String> {
    fs::read(path).map_err(|error| format!("read {}: {error}", path.display()))
}

fn digest(bytes: &[u8]) -> String {
    hex::encode(sha256(bytes))
}

fn parse(path: &Path) -> Result<(Value, String), String> {
    let bytes = read(path)?;
    let value = serde_json::from_slice(&bytes)
        .map_err(|error| format!("parse {}: {error}", path.display()))?;
    Ok((value, digest(&bytes)))
}

fn required_u64(value: &Value, key: &str) -> Result<u64, String> {
    value[key]
        .as_u64()
        .ok_or_else(|| format!("missing integer {key}"))
}

fn median(values: &[u64]) -> u64 {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    sorted[sorted.len() / 2]
}

fn real_time(path: &Path) -> Result<f64, String> {
    let text = String::from_utf8(read(path)?)
        .map_err(|error| format!("timer text {}: {error}", path.display()))?;
    let line = text
        .lines()
        .find(|line| line.starts_with("real "))
        .ok_or_else(|| format!("missing process real time in {}", path.display()))?;
    line[5..]
        .parse::<f64>()
        .map_err(|error| format!("process real time: {error}"))
}

fn phase_sum(value: &Value, keys: &[&str]) -> Result<u64, String> {
    keys.iter().try_fold(0u64, |sum, key| {
        value[*key]
            .as_u64()
            .and_then(|duration| sum.checked_add(duration))
            .ok_or_else(|| format!("missing or overflowing phase {key}"))
    })
}

fn run(output: &Path) -> Result<(), String> {
    let root = Path::new(NOTE);
    let source_hashes = json!({
        "ic":digest(&read(Path::new("examples/n37_native_m6_one_target.rs"))?),
        "rho":digest(&read(Path::new("examples/koblitz_rho_batch_ks_strong_online.rs"))?),
        "replay":digest(&read(Path::new("examples/n37_native_m6_one_target_replay.rs"))?),
        "runner":digest(&read(&root.join("run_panel.sh"))?),
        "host":digest(&read(&root.join("HOST.json"))?),
        "protocol":digest(&read(&root.join("PROTOCOL.md"))?)
    });
    let mut rows = Vec::new();
    let mut per_point = Vec::new();
    let mut direct = 0usize;
    let mut shifted = 0usize;
    let mut all_attempts = 0usize;
    let mut max_attempts = 0usize;
    let mut min_ic_online = u64::MAX;
    let mut max_ic_online = 0u64;
    let mut min_rho_online = u64::MAX;
    let mut max_rho_online = 0u64;

    for index in 0..32usize {
        let mut answer: Option<u64> = None;
        let mut frozen_ic: Option<Value> = None;
        let mut ic_times = Vec::new();
        let mut rho_times = Vec::new();
        let mut rho_steps = Vec::new();
        let mut ic_real = Vec::new();
        let mut rho_real = Vec::new();
        for repeat in 0..5usize {
            let stem = format!("q{index:02}_r{repeat}");
            let prefix = root.join("raw").join(&stem);
            let ic_path = prefix.with_extension("ic.json");
            let rho_path = prefix.with_extension("rho.jsonl");
            let replay_path = prefix.with_extension("replay.json");
            let (ic, ic_sha) = parse(&ic_path)?;
            let (replay, replay_sha) = parse(&replay_path)?;
            let rho_bytes = read(&rho_path)?;
            let rho_sha = digest(&rho_bytes);
            let rho_text = String::from_utf8(rho_bytes)
                .map_err(|error| format!("rho text {stem}: {error}"))?;
            let rho_rows: Vec<Value> = rho_text
                .lines()
                .map(|line| serde_json::from_str(line).map_err(|error| error.to_string()))
                .collect::<Result<_, _>>()?;
            if rho_rows.len() != 2 {
                return Err(format!("{stem}: rho row count"));
            }
            let (rho, rho_summary) = (&rho_rows[0], &rho_rows[1]);
            let scalar = required_u64(&ic["target"], "recovered_log")?;
            let ic_online = required_u64(&ic, "online_ns")?;
            let rho_online = required_u64(rho, "online_ns")?;
            let expected_seed = 2026100110102u64 + 1000 * index as u64 + repeat as u64;
            if ic["schema"] != "n37-native-m6-one-target-v1"
                || ic["target_count"] != 1
                || ic["target_index"] != index
                || ic["target"]["index"] != index
                || ic["rank"] != 42
                || ic["complete_one_target"] != true
                || replay["schema"] != "n37-native-m6-one-target-source-replay-v1"
                || replay["complete_one_target"] != true
                || replay["online_phase_sums_exact"] != true
                || replay["rho_source_scalar_verified"] != true
                || replay["target_index"] != index
                || replay["raw_result_sha256"] != ic_sha
                || replay["rho_result_sha256"] != rho_sha
                || rho["published_q"] != ic["target"]["source"]
                || rho["recovered_fixture_scalar"] != scalar
                || rho["verified"] != true
                || rho_summary["all_verified"] != true
                || rho_summary["fixtures"] != 1
                || rho_summary["batch_seed"] != expected_seed
                || replay["rho_seed"] != expected_seed
                || phase_sum(
                    &ic["online_phase_ns"],
                    &[
                        "target_query",
                        "target_pdp",
                        "target_relation_check",
                        "target_descent",
                        "recovery_check",
                    ],
                )? != ic_online
                || phase_sum(rho, &["walk_ns", "collision_ns", "recovery_check_ns"])? != rho_online
            {
                return Err(format!(
                    "{stem}: identity, correctness, or phase audit failed"
                ));
            }
            if let Some(prior) = answer {
                if scalar != prior
                    || frozen_ic.as_ref()
                        != Some(&json!({
                            "base_logs":ic["base_logs"],
                            "relations":ic["relations"],
                            "shifts":ic["shifts"],
                            "attempts":ic["target"]["attempts"],
                            "build_counts":ic["build_counts"],
                            "relation_oracle_counts":ic["relation_oracle_counts"],
                            "direct_oracle_counts":ic["direct_oracle_counts"],
                            "residual_oracle_counts":ic["residual_oracle_counts"],
                            "scalar_multiplications":ic["scalar_multiplications"]
                        }))
                {
                    return Err(format!("{stem}: IC answer or deterministic work changed"));
                }
            } else {
                answer = Some(scalar);
                frozen_ic = Some(json!({
                    "base_logs":ic["base_logs"],
                    "relations":ic["relations"],
                    "shifts":ic["shifts"],
                    "attempts":ic["target"]["attempts"],
                    "build_counts":ic["build_counts"],
                    "relation_oracle_counts":ic["relation_oracle_counts"],
                    "direct_oracle_counts":ic["direct_oracle_counts"],
                    "residual_oracle_counts":ic["residual_oracle_counts"],
                    "scalar_multiplications":ic["scalar_multiplications"]
                }));
                direct += ic["direct_hits"].as_u64().unwrap_or(0) as usize;
                shifted += ic["shifted_hits"].as_u64().unwrap_or(0) as usize;
                let attempts = required_u64(&ic, "total_target_attempts")? as usize;
                all_attempts += attempts;
                max_attempts = max_attempts.max(attempts);
            }
            let ic_real_s = real_time(&prefix.with_extension("ic.stderr"))?;
            let rho_real_s = real_time(&prefix.with_extension("rho.stderr"))?;
            ic_times.push(ic_online);
            rho_times.push(rho_online);
            rho_steps.push(required_u64(rho, "walk_steps")?);
            ic_real.push(ic_real_s);
            rho_real.push(rho_real_s);
            min_ic_online = min_ic_online.min(ic_online);
            max_ic_online = max_ic_online.max(ic_online);
            min_rho_online = min_rho_online.min(rho_online);
            max_rho_online = max_rho_online.max(rho_online);
            rows.push(json!({
                "candidate_id":"n37-native-m6-42-shift16",
                "workload_id":format!("n37-b02-q{index:02}"),
                "run_id":format!("r{repeat}"),
                "ic_sha256":ic_sha,"rho_sha256":rho_sha,"replay_sha256":replay_sha,
                "scalar":scalar,"ic_online_ns":ic_online,"rho_online_ns":rho_online,
                "ic_process_real_s":ic_real_s,"rho_process_real_s":rho_real_s,
                "rho_seed":expected_seed,"rho_walk_steps":rho["walk_steps"]
            }));
        }
        per_point.push(json!({
            "workload_id":format!("n37-b02-q{index:02}"),
            "verified_scalar":answer,
            "attempts":frozen_ic.as_ref().unwrap()["attempts"].as_array().unwrap().len(),
            "ic_online_ns_median":median(&ic_times),
            "rho_online_ns_median":median(&rho_times),
            "rho_walk_steps_median":median(&rho_steps),
            "ic_online_ns_by_repeat":ic_times,
            "rho_online_ns_by_repeat":rho_times,
            "ic_process_real_s_by_repeat":ic_real,
            "rho_process_real_s_by_repeat":rho_real
        }));
    }
    if direct + shifted != 32 || rows.len() != 160 {
        return Err("incomplete frozen panel".into());
    }
    let report = json!({
        "schema":"n37-native-one-target-panel-analysis-v1",
        "protocol":format!("{NOTE}/PROTOCOL.md"),
        "source_sha256":source_hashes,
        "host_resource_envelope":serde_json::from_slice::<Value>(&read(&root.join("HOST.json"))?)
            .map_err(|error| format!("host JSON: {error}"))?["resource_envelope"],
        "point_count":32,"repeats_per_point":5,"pair_count":160,
        "all_ic_rho_source_replays_pass":true,
        "all_ic_deterministic_outputs_and_counts_match_by_point":true,
        "direct_points":direct,"shifted_points":shifted,
        "total_queries_first_repetition":all_attempts,
        "maximum_queries_per_point":max_attempts,
        "online_ns_range_all_runs":{
            "ic":[min_ic_online,max_ic_online],"rho":[min_rho_online,max_rho_online]
        },
        "measurement_status":"complete one-target correctness and raw counts; unisolated descriptive time; common calibrated operation unit and S unset",
        "S":Value::Null,"speedup":Value::Null,
        "per_point":per_point,"runs":rows
    });
    fs::write(output, serde_json::to_vec_pretty(&report).unwrap())
        .map_err(|error| format!("write {}: {error}", output.display()))?;
    println!("32/32 Q, 160/160 replayed pairs; direct={direct}, shifted={shifted}, queries={all_attempts}, max={max_attempts}");
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        eprintln!("usage: n37_native_one_target_analyze SUMMARY.json");
        std::process::exit(2);
    }
    if let Err(error) = run(Path::new(&args[1])) {
        eprintln!("n37_native_one_target_analyze: {error}");
        std::process::exit(1);
    }
}
