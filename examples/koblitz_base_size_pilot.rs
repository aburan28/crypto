//! Recompute the frozen compact-orbit K pilot and select without held-out data.
//! This analyzer reads only pilot paths from the preregistered experiment.

use serde_json::{json, Value};
use std::path::Path;

const STUDY: &str = "experiments/koblitz-base-size-cold-panel-20261004";

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

fn nonnegative(value: &Value, key: &str) -> Result<f64, String> {
    let result = value[key]
        .as_f64()
        .ok_or_else(|| format!("missing {key}"))?;
    if !result.is_finite() || result < 0.0 {
        return Err(format!("invalid {key}"));
    }
    Ok(result)
}

fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    sorted[sorted.len() / 2]
}

fn cell(n: u64, ks: &[u64], hash_seed: u64, rho_seed: u64) -> Result<Value, String> {
    let target = one(format!("{STUDY}/targets/n{n}_pilot.jsonl"))?;
    if !target.is_array() || target.as_array().unwrap().len() != 2 {
        return Err(format!("n{n}: malformed frozen pilot point"));
    }
    let baseline_k = *ks.last().ok_or("empty K grid")?;
    let mut rows = Vec::new();
    let mut summaries = Vec::new();
    let mut eligible_smaller = Vec::new();
    let mut baseline_median = None;

    for &k in ks {
        let mut verified_cold = Vec::new();
        let mut actual_base = None;
        let mut base_hash = None;
        let mut entries = None;
        let mut rank_attempts = Vec::new();
        let mut rank_failures = Vec::new();
        let mut solved_targets = Vec::new();
        for repeat in 1..=3 {
            let dir = format!("{STUDY}/runs/pilot/n{n}/k{k}/r{repeat}");
            let status = one(format!("{dir}/status.jsonl"))?;
            if status["phase"] != "pilot"
                || status["n"] != n
                || status["K"] != k
                || status["repeat"] != repeat
            {
                return Err(format!("n{n} K{k} r{repeat}: status identity differs"));
            }
            let rho_status = whole(&status, "rho_status")?;
            let ic_status = whole(&status, "ic_status")?;
            let replay_status = status["replay_status"].as_u64();
            let mut diagnostics = Value::Null;
            let mut ic_cold = None;
            let mut rho_cold = None;
            let mut ic_online = None;
            let mut rho_online = None;
            if ic_status == 0 {
                let summary = one(format!("{dir}/ic_summary.jsonl"))?;
                if summary["n"] != n || summary["orbit_columns"] != k {
                    return Err(format!("n{n} K{k} r{repeat}: producer identity differs"));
                }
                let b = whole(&summary, "factor_base_points")?;
                let hash = summary["base_hash"]
                    .as_str()
                    .ok_or_else(|| format!("n{n} K{k}: base hash missing"))?;
                let root_entries = whole(&summary, "root_table_entries")?;
                if actual_base.is_some_and(|prior| prior != b)
                    || base_hash.as_ref().is_some_and(|prior| prior != hash)
                    || entries.is_some_and(|prior| prior != root_entries)
                {
                    return Err(format!("n{n} K{k}: factor-base or index drift"));
                }
                actual_base = Some(b);
                base_hash = Some(hash.to_owned());
                entries = Some(root_entries);
                rank_attempts.push(whole(&summary, "rank_attempts")?);
                rank_failures.push(whole(&summary, "rank_failures")?);
                solved_targets.push(whole(&summary, "targets_solved")?);
                ic_cold = Some(nonnegative(&summary, "cold_in_process_ms")?);
                diagnostics = json!({
                    "factor_base_points": b,
                    "orbit_columns": k,
                    "root_table_entries": root_entries,
                    "rank": whole(&summary, "rank")?,
                    "rank_attempts": whole(&summary, "rank_attempts")?,
                    "rank_failures": whole(&summary, "rank_failures")?,
                    "targets_solved": whole(&summary, "targets_solved")?,
                    "targets_failed": whole(&summary, "targets_failed")?,
                    "rank_pdp_ms": nonnegative(&summary["cold_phase_ms"], "rank_pdp")?,
                    "precompute_index_ms": nonnegative(&summary["cold_phase_ms"], "precompute_index")?,
                    "rank_matrix_build_ms": nonnegative(&summary["cold_phase_ms"], "rank_matrix_build")?,
                    "peak_rss_bytes": whole(&summary, "peak_rss_bytes")?
                });
            }
            if replay_status == Some(0) {
                if rho_status != 0 || ic_status != 0 {
                    return Err(format!("n{n} K{k} r{repeat}: replay passed failed arm"));
                }
                let replay = one(format!("{dir}/replay.jsonl"))?;
                let summary = one(format!("{dir}/ic_summary.jsonl"))?;
                let ic = one(format!("{dir}/ic_targets.jsonl"))?;
                let rho = one(format!("{dir}/rho.jsonl"))?;
                if replay["verified"] != true
                    || replay["cold_phase_verified"] != true
                    || replay["n"] != n
                    || replay["orbit_columns"] != k
                    || replay["public_hash_seed"] != hash_seed
                    || replay["rho_batch_seed"] != rho_seed
                    || replay["replayed_base_logs"] != summary["factor_base_points"]
                    || summary["rank"] != k
                    || summary["targets_solved"] != 1
                    || summary["targets_failed"] != 0
                    || ic["published_q"] != target
                    || rho["published_q"] != target
                    || rho["public_hash_seed"] != hash_seed
                    || rho["batch_seed"] != rho_seed
                    || rho["verified"] != true
                    || ic["group_verified"] != true
                    || ic["recovered_scalar"] != rho["recovered_fixture_scalar"]
                {
                    return Err(format!("n{n} K{k} r{repeat}: independent replay mismatch"));
                }
                let measured_cold = nonnegative(&replay, "ic_cold_in_process_ms")?;
                if ic_cold.is_none_or(|reported| (reported - measured_cold).abs() > 0.000001) {
                    return Err(format!("n{n} K{k} r{repeat}: cold interval differs"));
                }
                ic_online = Some(nonnegative(&replay, "ic_online_ms")?);
                rho_online = Some(nonnegative(&replay, "rho_online_ms")?);
                rho_cold = Some(nonnegative(&replay, "rho_cold_in_process_ms")?);
                verified_cold.push(measured_cold);
            }
            let reason = if replay_status != Some(0) {
                std::fs::read_to_string(format!("{dir}/replay.stderr"))
                    .ok()
                    .map(|text| text.trim().to_owned())
                    .filter(|text| !text.is_empty())
            } else {
                None
            };
            rows.push(json!({
                "n":n,"K":k,"repeat":repeat,
                "rho_status":rho_status,"ic_status":ic_status,"replay_status":replay_status,
                "verified":replay_status==Some(0),
                "failure_reason":reason,
                "ic_cold_ms":ic_cold,"rho_cold_ms":rho_cold,
                "ic_online_ms":ic_online,"rho_online_ms":rho_online,
                "diagnostics":diagnostics
            }));
        }
        let eligible = verified_cold.len() == 3;
        let median_cold = eligible.then(|| median(&verified_cold));
        if k == baseline_k {
            baseline_median = median_cold;
        } else if let Some(cost) = median_cold {
            eligible_smaller.push((k, cost));
        }
        summaries.push(json!({
            "K":k,"actual_base_points":actual_base,"base_hash":base_hash,
            "root_table_entries":entries,"verified_pairs":verified_cold.len(),
            "eligible":eligible,"median_verified_ic_cold_ms":median_cold,
            "rank_attempts":rank_attempts,"rank_failures":rank_failures,
            "solved_targets":solved_targets
        }));
    }
    eligible_smaller.sort_by(|(ka, ca), (kb, cb)| ca.total_cmp(cb).then(ka.cmp(kb)));
    let selected = eligible_smaller.first().map(|(k, _)| *k);
    Ok(json!({
        "n":n,"pilot_target":target,"public_hash_seed":hash_seed,
        "rank_and_rho_seed":rho_seed,"baseline_K":baseline_k,
        "baseline_median_ic_cold_ms":baseline_median,
        "selected_smaller_K":selected,
        "selection_rule":"all three pairs verified; lowest median IC cold among smaller K; exact tie to smaller K",
        "K_summaries":summaries,"pair_rows":rows
    }))
}

fn main() {
    let result = (|| -> Result<Value, String> {
        let n41 = cell(41, &[20, 32, 48, 64, 85], 41261107, 410041)?;
        let n53 = cell(53, &[80, 100, 128, 160, 220], 53261107, 530053)?;
        Ok(json!({
            "schema":"koblitz_base_size_pilot/v1",
            "source_commit":"e37c50e9dbbc88148a12be648735894c37ba0b1a",
            "producer_evidence_class":"exploratory_unisolated_native_example_diagnostic",
            "heldout_data_read":false,
            "formal_ecbench_vs_rho_claim":false,
            "cells":[n41,n53]
        }))
    })();
    match result {
        Ok(value) => println!("{value}"),
        Err(error) => {
            eprintln!("pilot analysis failed: {error}");
            std::process::exit(1);
        }
    }
}
