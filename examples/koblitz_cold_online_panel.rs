//! Deterministic analysis of the frozen n41/n53 one-target diagnostic panel.
//! Run from the repository root after every pair has independently replayed.

use serde_json::{json, Value};
use std::path::Path;

const STUDY: &str = "experiments/koblitz-n41-n53-single-target-20261004";

fn read_one(path: impl AsRef<Path>) -> Result<Value, String> {
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

fn measured(value: &Value, key: &str) -> Result<f64, String> {
    let v = value[key]
        .as_f64()
        .ok_or_else(|| format!("missing {key}"))?;
    if !v.is_finite() || v < 0.0 {
        return Err(format!("invalid {key}"));
    }
    Ok(v)
}

fn whole(value: &Value, key: &str) -> Result<u64, String> {
    value[key].as_u64().ok_or_else(|| format!("missing {key}"))
}

fn stats(values: &[f64]) -> Value {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let len = sorted.len();
    let median = if len.is_multiple_of(2) {
        (sorted[len / 2 - 1] + sorted[len / 2]) / 2.0
    } else {
        sorted[len / 2]
    };
    json!({"median": median, "min": sorted[0], "max": sorted[len - 1]})
}

fn analyze_cell(n: u32, columns: u64, hash_seed: u64, rho_seed: u64) -> Result<Value, String> {
    let (curve_slug, curve_ec1) = match n {
        41 => ("icv1-f2m41-tm2308219-7f48b14a", "EC1N41Ce0he09550ab560a"),
        53 => ("icv1-f2m53-tm56619371-dac20a85", "EC1N53Ce0hb097de99be9a"),
        _ => return Err(format!("unsupported frozen curve n={n}")),
    };
    let target = read_one(format!("{STUDY}/targets/n{n}.jsonl"))?;
    let mut base_hash = None;
    let mut point_count = None;
    let mut root_table_entries = None;
    let mut scalar = None;
    let mut pairs = Vec::new();
    let mut ic_online = Vec::new();
    let mut rho_online = Vec::new();
    let mut ic_cold = Vec::new();
    let mut rho_cold = Vec::new();
    let mut online_ratio = Vec::new();
    let mut cold_ratio = Vec::new();
    let mut ic_over_rho_online = Vec::new();
    let mut ic_over_rho_cold = Vec::new();
    let mut rank_pdp = Vec::new();
    let mut index_build = Vec::new();
    let mut rank_matrix = Vec::new();
    let mut cold_without_rank_pdp = Vec::new();
    let mut cold_without_rank_pdp_and_index = Vec::new();
    let mut rank_attempts = Vec::new();
    let mut rank_failures = Vec::new();
    let mut max_ic_rss = 0;
    let mut max_rho_rss = 0;

    for repeat in 1..=6 {
        let dir = format!("{STUDY}/runs/n{n}/r{repeat}");
        let receipt = read_one(format!("{dir}/replay.jsonl"))?;
        let summary = read_one(format!("{dir}/ic_summary.jsonl"))?;
        let ic = read_one(format!("{dir}/ic_targets.jsonl"))?;
        let rho = read_one(format!("{dir}/rho.jsonl"))?;
        let base = read_one(format!("{dir}/base.jsonl"))?;

        if receipt["verified"] != true
            || receipt["cold_phase_verified"] != true
            || receipt["n"] != n
            || receipt["public_hash_seed"] != hash_seed
            || receipt["rho_batch_seed"] != rho_seed
            || receipt["orbit_columns"] != columns
            || summary["rank"] != columns
            || summary["targets_solved"] != 1
            || summary["targets_failed"] != 0
            || summary["base_hash"] != base["base_hash"]
            || ic["published_q"] != target
            || rho["published_q"] != target
            || rho["verified"] != true
            || ic["group_verified"] != true
            || ic["recovered_scalar"] != rho["recovered_fixture_scalar"]
        {
            return Err(format!("n{n} r{repeat}: frozen input or replay mismatch"));
        }
        if let Some(previous) = &base_hash {
            if previous != &base["base_hash"] {
                return Err(format!("n{n} r{repeat}: base digest drift"));
            }
        } else {
            base_hash = Some(base["base_hash"].clone());
        }
        if let Some(previous) = point_count {
            if previous != whole(&receipt, "factor_base_points")? {
                return Err(format!("n{n} r{repeat}: actual base size drift"));
            }
        } else {
            point_count = Some(whole(&receipt, "factor_base_points")?);
        }
        if let Some(previous) = root_table_entries {
            if previous != whole(&summary, "root_table_entries")? {
                return Err(format!("n{n} r{repeat}: root-index size drift"));
            }
        } else {
            root_table_entries = Some(whole(&summary, "root_table_entries")?);
        }
        if let Some(previous) = scalar {
            if previous != whole(&receipt, "target_scalar")? {
                return Err(format!("n{n} r{repeat}: recovered scalar drift"));
            }
        } else {
            scalar = Some(whole(&receipt, "target_scalar")?);
        }

        let i_online = measured(&receipt, "ic_online_ms")?;
        let r_online = measured(&receipt, "rho_online_ms")?;
        let i_cold = measured(&receipt, "ic_cold_in_process_ms")?;
        let r_cold = measured(&receipt, "rho_cold_in_process_ms")?;
        if i_online == 0.0 || i_cold == 0.0 || r_online == 0.0 || r_cold == 0.0 {
            return Err(format!("n{n} r{repeat}: zero timing denominator"));
        }
        let o_ratio = r_online / i_online;
        let c_ratio = r_cold / i_cold;
        let i_over_r_online = i_online / r_online;
        let i_over_r_cold = i_cold / r_cold;
        let rank_pdp_ms = measured(&summary["cold_phase_ms"], "rank_pdp")?;
        let index_ms = measured(&summary["cold_phase_ms"], "precompute_index")?;
        let matrix_ms = measured(&summary["cold_phase_ms"], "rank_matrix_build")?;
        let attempts = whole(&summary, "rank_attempts")?;
        let failures = whole(&summary, "rank_failures")?;
        let i_rss = whole(&receipt, "peak_rss_bytes")?;
        let r_rss = whole(&receipt, "rho_peak_rss_bytes")?;
        max_ic_rss = max_ic_rss.max(i_rss);
        max_rho_rss = max_rho_rss.max(r_rss);
        ic_online.push(i_online);
        rho_online.push(r_online);
        ic_cold.push(i_cold);
        rho_cold.push(r_cold);
        online_ratio.push(o_ratio);
        cold_ratio.push(c_ratio);
        ic_over_rho_online.push(i_over_r_online);
        ic_over_rho_cold.push(i_over_r_cold);
        rank_pdp.push(rank_pdp_ms);
        index_build.push(index_ms);
        rank_matrix.push(matrix_ms);
        cold_without_rank_pdp.push(i_cold - rank_pdp_ms);
        cold_without_rank_pdp_and_index.push(i_cold - rank_pdp_ms - index_ms);
        rank_attempts.push(attempts);
        rank_failures.push(failures);
        pairs.push(json!({
            "repeat": repeat,
            "ic_online_ms": i_online,
            "rho_online_ms": r_online,
            "rho_over_ic_online": o_ratio,
            "ic_over_rho_online": i_over_r_online,
            "ic_cold_ms": i_cold,
            "rho_cold_ms": r_cold,
            "rho_over_ic_cold": c_ratio,
            "ic_over_rho_cold": i_over_r_cold,
            "rank_pdp_ms": rank_pdp_ms,
            "precompute_index_ms": index_ms,
            "rank_matrix_build_ms": matrix_ms,
            "counterfactual_zero_rank_pdp_ms": i_cold - rank_pdp_ms,
            "counterfactual_zero_rank_pdp_and_index_ms": i_cold - rank_pdp_ms - index_ms,
            "rank_attempts": attempts,
            "rank_failures": failures,
            "ic_peak_rss_bytes": i_rss,
            "rho_peak_rss_bytes": r_rss,
            "replay": "pass"
        }));
    }

    Ok(json!({
        "n": n,
        "a": 0,
        "curve_slug": curve_slug,
        "curve_ec1": curve_ec1,
        "candidate_id": null,
        "workload_id": null,
        "identity_status": "producer diagnostic; no formal IC1 candidate/workload manifest",
        "target": target,
        "public_hash_seed": hash_seed,
        "rho_batch_seed": rho_seed,
        "factor_base_points": point_count,
        "orbit_columns": columns,
        "root_table_entries": root_table_entries,
        "base_hash": base_hash,
        "recovered_scalar": scalar,
        "verified_pairs": pairs.len(),
        "rank_attempts": rank_attempts,
        "rank_failures": rank_failures,
        "max_ic_rss_bytes": max_ic_rss,
        "max_rho_rss_bytes": max_rho_rss,
        "pair_rows": pairs,
        "timing_spread": {
            "ic_online_ms": stats(&ic_online),
            "rho_online_ms": stats(&rho_online),
            "rho_over_ic_online": stats(&online_ratio),
            "ic_over_rho_online": stats(&ic_over_rho_online),
            "ic_cold_ms": stats(&ic_cold),
            "rho_cold_ms": stats(&rho_cold),
            "rho_over_ic_cold": stats(&cold_ratio),
            "ic_over_rho_cold": stats(&ic_over_rho_cold),
            "rank_pdp_ms": stats(&rank_pdp),
            "precompute_index_ms": stats(&index_build),
            "rank_matrix_build_ms": stats(&rank_matrix),
            "counterfactual_zero_rank_pdp_ms": stats(&cold_without_rank_pdp),
            "counterfactual_zero_rank_pdp_and_index_ms": stats(&cold_without_rank_pdp_and_index)
        }
    }))
}

fn main() {
    let result = (|| -> Result<Value, String> {
        let n41 = analyze_cell(41, 85, 41_261_004, 410_041)?;
        let n53 = analyze_cell(53, 220, 53_261_004, 530_053)?;
        Ok(json!({
            "schema": "koblitz_cold_online_panel/v1",
            "evidence_class": "exploratory_unisolated_native_example_diagnostic",
            "formal_ecbench_vs_rho_claim": false,
            "controlled_wall_speedup": null,
            "same_point_pairs_per_cell": 6,
            "cells": [n41, n53]
        }))
    })();
    match result {
        Ok(value) => println!("{value}"),
        Err(error) => {
            eprintln!("panel analysis failed: {error}");
            std::process::exit(1);
        }
    }
}
