//! Recompute the frozen compact-orbit held-out K comparison from raw rows.

use serde_json::{json, Value};
use std::path::Path;

const STUDY: &str = "experiments/koblitz-base-size-cold-panel-20261004";
const RSS_LIMIT: u64 = 16 * 1024 * 1024 * 1024;

fn one(path: impl AsRef<Path>) -> Result<Value, String> {
    let path = path.as_ref();
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let rows: Vec<&[u8]> = bytes
        .split(|&byte| byte == b'\n')
        .filter(|row| !row.is_empty())
        .collect();
    if rows.len() != 1 {
        return Err(format!("{}: expected one JSON row", path.display()));
    }
    serde_json::from_slice(rows[0]).map_err(|e| format!("{}: {e}", path.display()))
}

fn whole(value: &Value, key: &str) -> Result<u64, String> {
    value[key].as_u64().ok_or_else(|| format!("missing {key}"))
}

fn measured(value: &Value, key: &str) -> Result<f64, String> {
    let result = value[key]
        .as_f64()
        .ok_or_else(|| format!("missing {key}"))?;
    if !result.is_finite() || result < 0.0 {
        return Err(format!("invalid {key}"));
    }
    Ok(result)
}

fn stats(values: &[f64]) -> Value {
    if values.is_empty() {
        return Value::Null;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let middle = sorted.len() / 2;
    let median = if sorted.len() % 2 == 0 {
        (sorted[middle - 1] + sorted[middle]) / 2.0
    } else {
        sorted[middle]
    };
    json!({"median":median,"min":sorted[0],"max":sorted[sorted.len()-1]})
}

struct Arm {
    report: Value,
    cold_by_repeat: Vec<Option<f64>>,
    online_by_repeat: Vec<Option<f64>>,
}

fn arm(n: u64, k: u64, target: &Value, hash_seed: u64, rho_seed: u64) -> Result<Arm, String> {
    let mut rows = Vec::new();
    let mut cold_by_repeat = vec![None; 6];
    let mut online_by_repeat = vec![None; 6];
    let mut ic_cold = Vec::new();
    let mut rho_cold = Vec::new();
    let mut ic_online = Vec::new();
    let mut rho_online = Vec::new();
    let mut cold_ratio = Vec::new();
    let mut online_ratio = Vec::new();
    let mut rho_over_ic_cold = Vec::new();
    let mut rho_over_ic_online = Vec::new();
    let mut rank_pdp = Vec::new();
    let mut index_build = Vec::new();
    let mut rank_matrix = Vec::new();
    let mut rank_probes_mean = Vec::new();
    let mut actual_base = None;
    let mut base_hash = None;
    let mut root_entries = None;

    for repeat in 1..=6 {
        let dir = format!("{STUDY}/runs/holdout/n{n}/k{k}/r{repeat}");
        let status = one(format!("{dir}/status.jsonl"))?;
        if status["phase"] != "holdout"
            || status["n"] != n
            || status["K"] != k
            || status["repeat"] != repeat
        {
            return Err(format!("n{n} K{k} r{repeat}: status identity differs"));
        }
        let rho_status = whole(&status, "rho_status")?;
        let ic_status = whole(&status, "ic_status")?;
        let replay_status = status["replay_status"].as_u64();
        if rho_status != 0 || ic_status != 0 || replay_status != Some(0) {
            let reason = std::fs::read_to_string(format!("{dir}/replay.stderr"))
                .ok()
                .map(|text| text.trim().to_owned())
                .filter(|text| !text.is_empty());
            rows.push(json!({
                "repeat":repeat,"verified":false,"rho_status":rho_status,
                "ic_status":ic_status,"replay_status":replay_status,
                "failure_reason":reason
            }));
            continue;
        }

        let replay = one(format!("{dir}/replay.jsonl"))?;
        let summary = one(format!("{dir}/ic_summary.jsonl"))?;
        let ic = one(format!("{dir}/ic_targets.jsonl"))?;
        let rho = one(format!("{dir}/rho.jsonl"))?;
        let base = one(format!("{dir}/base.jsonl"))?;
        let b = whole(&summary, "factor_base_points")?;
        let hash = summary["base_hash"]
            .as_str()
            .ok_or_else(|| format!("n{n} K{k}: base hash missing"))?;
        let entries = whole(&summary, "root_table_entries")?;
        if actual_base.is_some_and(|prior| prior != b)
            || base_hash.as_ref().is_some_and(|prior| prior != hash)
            || root_entries.is_some_and(|prior| prior != entries)
        {
            return Err(format!("n{n} K{k}: base or index drift"));
        }
        actual_base = Some(b);
        base_hash = Some(hash.to_owned());
        root_entries = Some(entries);
        if replay["verified"] != true
            || replay["cold_phase_verified"] != true
            || replay["n"] != n
            || replay["orbit_columns"] != k
            || replay["public_hash_seed"] != hash_seed
            || replay["rho_batch_seed"] != rho_seed
            || replay["replayed_base_logs"] != b
            || summary["n"] != n
            || summary["orbit_columns"] != k
            || summary["rank"] != k
            || summary["targets_solved"] != 1
            || summary["targets_failed"] != 0
            || summary["base_hash"] != base["base_hash"]
            || ic["published_q"] != *target
            || rho["published_q"] != *target
            || rho["verified"] != true
            || rho["reference_grade"] != "strong"
            || rho["quotient_mode"] != "signed_frobenius"
            || rho["target_kind"] != "public_hash_to_curve_cofactor"
            || rho["public_hash_seed"] != hash_seed
            || rho["batch_seed"] != rho_seed
            || ic["group_verified"] != true
            || ic["recovered_scalar"] != rho["recovered_fixture_scalar"]
            || whole(&replay, "peak_rss_bytes")? > RSS_LIMIT
            || whole(&replay, "rho_peak_rss_bytes")? > RSS_LIMIT
        {
            return Err(format!(
                "n{n} K{k} r{repeat}: replay or frozen-input mismatch"
            ));
        }
        let i_cold = measured(&replay, "ic_cold_in_process_ms")?;
        let r_cold = measured(&replay, "rho_cold_in_process_ms")?;
        let i_online = measured(&replay, "ic_online_ms")?;
        let r_online = measured(&replay, "rho_online_ms")?;
        if i_cold == 0.0 || r_cold == 0.0 || i_online == 0.0 || r_online == 0.0 {
            return Err(format!("n{n} K{k} r{repeat}: zero timing denominator"));
        }
        if (i_cold - measured(&summary, "cold_in_process_ms")?).abs() > 0.000001 {
            return Err(format!("n{n} K{k} r{repeat}: cold interval differs"));
        }
        let pdp = measured(&summary["cold_phase_ms"], "rank_pdp")?;
        let index = measured(&summary["cold_phase_ms"], "precompute_index")?;
        let matrix = measured(&summary["cold_phase_ms"], "rank_matrix_build")?;
        let probes = measured(&summary, "rank_probes_mean")?;
        ic_cold.push(i_cold);
        rho_cold.push(r_cold);
        ic_online.push(i_online);
        rho_online.push(r_online);
        cold_ratio.push(i_cold / r_cold);
        online_ratio.push(i_online / r_online);
        rho_over_ic_cold.push(r_cold / i_cold);
        rho_over_ic_online.push(r_online / i_online);
        rank_pdp.push(pdp);
        index_build.push(index);
        rank_matrix.push(matrix);
        rank_probes_mean.push(probes);
        cold_by_repeat[repeat as usize - 1] = Some(i_cold);
        online_by_repeat[repeat as usize - 1] = Some(i_online);
        rows.push(json!({
            "repeat":repeat,"verified":true,"rho_status":0,"ic_status":0,
            "replay_status":0,"ic_cold_ms":i_cold,"rho_cold_ms":r_cold,
            "ic_online_ms":i_online,"rho_online_ms":r_online,
            "ic_over_rho_cold":i_cold/r_cold,
            "ic_over_rho_online":i_online/r_online,
            "rho_over_ic_cold":r_cold/i_cold,
            "rho_over_ic_online":r_online/i_online,
            "rank_pdp_ms":pdp,"precompute_index_ms":index,
            "rank_matrix_build_ms":matrix,
            "rank_probes_mean":probes,
            "rank_attempts":whole(&summary,"rank_attempts")?,
            "rank_failures":whole(&summary,"rank_failures")?,
            "peak_rss_bytes":whole(&replay,"peak_rss_bytes")?,
            "rho_peak_rss_bytes":whole(&replay,"rho_peak_rss_bytes")?
        }));
    }
    Ok(Arm {
        report: json!({
            "K":k,"actual_base_points":actual_base,"base_hash":base_hash,
            "root_table_entries":root_entries,"verified_pairs":ic_cold.len(),
            "rows":rows,
            "spread":{
                "ic_cold_ms":stats(&ic_cold),"rho_cold_ms":stats(&rho_cold),
                "ic_online_ms":stats(&ic_online),"rho_online_ms":stats(&rho_online),
                "ic_over_rho_cold":stats(&cold_ratio),
                "ic_over_rho_online":stats(&online_ratio),
                "rho_over_ic_cold":stats(&rho_over_ic_cold),
                "rho_over_ic_online":stats(&rho_over_ic_online),
                "rank_pdp_ms":stats(&rank_pdp),
                "precompute_index_ms":stats(&index_build),
                "rank_matrix_build_ms":stats(&rank_matrix),
                "rank_probes_mean":stats(&rank_probes_mean)
            }
        }),
        cold_by_repeat,
        online_by_repeat,
    })
}

fn cell(
    n: u64,
    expected_baseline: u64,
    expected_selected: u64,
    hash_seed: u64,
    rho_seed: u64,
) -> Result<Value, String> {
    let pilot = one(format!("{STUDY}/PILOT_ANALYSIS.json"))?;
    let selection = pilot["cells"]
        .as_array()
        .ok_or("pilot cells missing")?
        .iter()
        .find(|cell| cell["n"] == n)
        .ok_or_else(|| format!("n{n} pilot selection missing"))?;
    if selection["baseline_K"] != expected_baseline
        || selection["selected_smaller_K"] != expected_selected
    {
        return Err(format!("n{n}: published selection differs"));
    }
    let target = one(format!("{STUDY}/targets/n{n}_holdout.jsonl"))?;
    if !target.is_array() || target.as_array().unwrap().len() != 2 {
        return Err(format!("n{n}: malformed held-out point"));
    }
    let baseline = arm(n, expected_baseline, &target, hash_seed, rho_seed)?;
    let selected = arm(n, expected_selected, &target, hash_seed, rho_seed)?;
    let mut paired_cold = Vec::new();
    let mut paired_online = Vec::new();
    let mut paired_rows = Vec::new();
    for repeat in 1..=6 {
        let index = repeat - 1;
        let b_cold = baseline.cold_by_repeat[index];
        let s_cold = selected.cold_by_repeat[index];
        let b_online = baseline.online_by_repeat[index];
        let s_online = selected.online_by_repeat[index];
        let cold_ratio = b_cold.zip(s_cold).map(|(b, s)| s / b);
        let online_ratio = b_online.zip(s_online).map(|(b, s)| s / b);
        if let Some(ratio) = cold_ratio {
            paired_cold.push(ratio);
        }
        if let Some(ratio) = online_ratio {
            paired_online.push(ratio);
        }
        paired_rows.push(json!({
            "repeat":repeat,"selected_over_baseline_ic_cold":cold_ratio,
            "selected_over_baseline_ic_online":online_ratio
        }));
    }
    Ok(json!({
        "n":n,"heldout_target":target,"public_hash_seed":hash_seed,
        "rank_and_rho_seed":rho_seed,"baseline_K":expected_baseline,
        "selected_smaller_K":expected_selected,
        "baseline":baseline.report,"selected":selected.report,
        "paired_rows":paired_rows,
        "paired_selected_over_baseline_ic_cold":stats(&paired_cold),
        "paired_selected_over_baseline_ic_online":stats(&paired_online)
    }))
}

fn main() {
    let result = (|| -> Result<Value, String> {
        let n41 = cell(41, 85, 64, 41261207, 410041)?;
        let n53 = cell(53, 220, 160, 53261207, 530053)?;
        Ok(json!({
            "schema":"koblitz_base_size_holdout/v1",
            "source_commit":"e37c50e9dbbc88148a12be648735894c37ba0b1a",
            "producer_evidence_class":"exploratory_unisolated_native_example_diagnostic",
            "formal_ecbench_vs_rho_claim":false,
            "controlled_wall_speedup":null,
            "selection_from_pilot_only":true,
            "cells":[n41,n53]
        }))
    })();
    match result {
        Ok(value) => println!("{value}"),
        Err(error) => {
            eprintln!("held-out analysis failed: {error}");
            std::process::exit(1);
        }
    }
}
