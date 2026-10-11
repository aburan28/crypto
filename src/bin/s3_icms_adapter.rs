//! ICMS command-adapter producer for the frozen S3 baseline/candidate binaries.
//!
//! Usage: s3_icms_adapter SOLVER N A COLUMNS RELATION_SEED PUBLIC_POINT.json
//!        WORKERS METRICS.json
//! The solver receives only the public point, never its fixture scalar.

use crypto_lib::hash::sha256::sha256;
use serde_json::{json, Value};
use std::env;
use std::fs;
use std::io::{self, Write};
use std::path::PathBuf;
use std::process::Command;

const PHASES: [(&str, &str); 5] = [
    ("target_query", "target_query"),
    ("target_pdp", "target_pdp_charged"),
    ("target_relation_check", "target_relation_check"),
    ("target_descent", "target_descent"),
    ("recovery_check", "target_recovery_check"),
];

fn as_ms(timing: &Value, key: &str) -> Result<f64, String> {
    let value = timing[key]
        .as_f64()
        .ok_or_else(|| format!("missing timing_ms.{key}"))?;
    if !value.is_finite() || value < 0.0 {
        return Err(format!("invalid timing_ms.{key}"));
    }
    Ok(value)
}

fn metrics(
    report: &Value,
    target: &Value,
    n: u64,
    a: u64,
    columns: u64,
    relation_seed: u64,
    workers: u64,
    binary_sha256: &str,
) -> Result<Value, String> {
    if report["n"] != n
        || report["a"] != a
        || report["orbit_columns"] != columns
        || report["relation_seed"] != relation_seed
        || report["target"] != *target
        || report["target_count"] != 1
        || report["relation_collection_workers"] != workers
        || report["group_verified"] != true
        || report["recovered_scalar"].is_null()
    {
        return Err("solver report fails instance, target, worker or replay check".into());
    }
    let timing = &report["timing_ms"];
    let online = as_ms(timing, "target_online_after_reusable_setup")?;
    if online <= 0.0 {
        return Err("online interval must be positive".into());
    }
    let mut phase_map = serde_json::Map::new();
    let mut phase_sum = 0.0;
    for (standard, native) in PHASES {
        let ms = as_ms(timing, native)?;
        phase_sum += ms;
        phase_map.insert(
            standard.into(),
            json!({"wall_ns": (ms * 1e6).round() as u64,
                                                 "ops": null, "ops_unit": null}),
        );
    }
    if (phase_sum - online).abs() > 0.02 {
        return Err("five exclusive phases do not sum to the online interval".into());
    }
    let online_ns = (online * 1e6).round() as u64;
    let s3_calls = report["target_s3_calls"]
        .as_u64()
        .ok_or("missing target S3 count")?
        + report["rank_target_s3_calls"]
            .as_u64()
            .ok_or("missing rank S3 count")?;
    Ok(json!({
        "outcome": {"status": "complete", "verified": true,
                    "recovered": report["recovered_scalar"], "target": target},
        "units": {
            "time.wall_ns": {"total": online_ns, "deterministic": false,
                             "host_dependent_because": ["target-online wall clock on the measured host"]},
            "count.s3_solves": {"total": s3_calls, "deterministic": true,
                                "host_dependent_because": []}
        },
        "reference": null,
        "metrics": {
            "instance": {"degree": n, "koblitz_a": a,
                         "r": report["subgroup_order"].to_string()},
            "factor_base": {"family": "explicit", "usable_points": report["factor_base_points"],
                            "columns": report["orbit_columns"],
                            "orbit_representatives": report["orbit_columns"],
                            "native": {"enumerated_set_digest": report["factor_base_digest"],
                                       "index_entries": report["root_table_entries"]}},
            "decomposition": {"oracle": "PDP4root", "arity": 4,
                              "targets_tried": report["rank_attempts"],
                              "relations_found": report["rank_verified_relations"],
                              "native": {"s3_calls": s3_calls,
                                         "target_s3_calls": report["target_s3_calls"],
                                         "rank_s3_calls": report["rank_target_s3_calls"]}},
            "linear_algebra": {"method": "incremental_gauss",
                               "rows": report["rank_new_rows"],
                               "rank": report["rank"],
                               "dependent": report["rank_dependent_rows"]}
        },
        "phases": phase_map,
        "windows": {
            "online_one_target": {"conformance": "exact", "wall_ns": online_ns,
                                  "ops": null, "ops_unit": null,
                                  "composition": {
                                      "target_query": "timing_ms.target_query",
                                      "target_pdp": "timing_ms.target_pdp_charged",
                                      "target_relation_check": "timing_ms.target_relation_check",
                                      "target_descent": "timing_ms.target_descent",
                                      "recovery_check": "timing_ms.target_recovery_check"
                                  }},
            "whole_process": {"conformance": "runner", "wall_ns": null,
                              "ops": null, "ops_unit": null}
        },
        "consistency": [],
        "producer": {"kind": "s3_pair_root_solver", "binary_sha256": binary_sha256,
                     "online_boundary": "after reusable base/index setup through recovered scalar point replay",
                     "native_report_kind": report["kind"]}
    }))
}

fn write_metrics(path: &PathBuf, value: &Value) -> Result<(), String> {
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    let mut bytes = serde_json::to_vec_pretty(value).map_err(|e| e.to_string())?;
    bytes.push(b'\n');
    fs::write(path, bytes).map_err(|e| e.to_string())
}

fn run() -> Result<(), String> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 9 {
        return Err("usage: s3_icms_adapter SOLVER N A COLUMNS RELATION_SEED PUBLIC_POINT.json WORKERS METRICS.json".into());
    }
    let solver = PathBuf::from(&args[1]);
    let n = args[2].parse::<u64>().map_err(|e| e.to_string())?;
    let a = args[3].parse::<u64>().map_err(|e| e.to_string())?;
    let columns = args[4].parse::<u64>().map_err(|e| e.to_string())?;
    let relation_seed = args[5].parse::<u64>().map_err(|e| e.to_string())?;
    let workers = args[7].parse::<u64>().map_err(|e| e.to_string())?;
    let metrics_path = PathBuf::from(&args[8]);
    let raw_path = metrics_path.with_extension("solver.jsonl");
    let binary_sha256 = hex::encode(sha256(&fs::read(&solver).map_err(|e| e.to_string())?));
    let target: Value = serde_json::from_slice(&fs::read(&args[6]).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    let output = Command::new(&solver)
        .args([&args[2], &args[3], &args[4], &args[5], &args[6]])
        .arg(&raw_path)
        .arg(&args[7])
        .output()
        .map_err(|e| e.to_string())?;
    io::stdout()
        .write_all(&output.stdout)
        .map_err(|e| e.to_string())?;
    io::stderr()
        .write_all(&output.stderr)
        .map_err(|e| e.to_string())?;
    if !output.status.success() {
        write_metrics(
            &metrics_path,
            &json!({"outcome": {"status": "error",
            "verified": false, "reason": format!("solver exit {}", output.status)}}),
        )?;
        return Err(format!("solver exit {}", output.status));
    }
    let report: Value = serde_json::from_slice(&fs::read(&raw_path).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    let binary_after = hex::encode(sha256(&fs::read(&solver).map_err(|e| e.to_string())?));
    if binary_after != binary_sha256 {
        return Err("solver binary changed during run".into());
    }
    let value = metrics(
        &report,
        &target,
        n,
        a,
        columns,
        relation_seed,
        workers,
        &binary_sha256,
    )?;
    write_metrics(&metrics_path, &value)
}

fn main() {
    if let Err(error) = run() {
        eprintln!("s3_icms_adapter: {error}");
        std::process::exit(2);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn report() -> Value {
        json!({
            "n": 41, "a": 0, "orbit_columns": 244, "relation_seed": 20260928,
            "target": [1, 2], "target_count": 1,
            "relation_collection_workers": 14, "group_verified": true,
            "recovered_scalar": 3, "target_s3_calls": 4,
            "rank_target_s3_calls": 5,
            "timing_ms": {
                "target_query": 1.0, "target_pdp_charged": 2.0,
                "target_relation_check": 3.0, "target_descent": 4.0,
                "target_recovery_check": 5.0,
                "target_online_after_reusable_setup": 15.0
            }
        })
    }

    #[test]
    fn exact_online_window_accepts_complete_verified_run() {
        let row = metrics(&report(), &json!([1, 2]), 41, 0, 244, 20260928, 14, "hash").unwrap();
        assert_eq!(row["windows"]["online_one_target"]["wall_ns"], 15_000_000);
        assert_eq!(row["phases"].as_object().unwrap().len(), 5);
    }

    #[test]
    fn changed_target_or_incomplete_phase_total_is_rejected() {
        let mut row = report();
        assert!(metrics(&row, &json!([9, 9]), 41, 0, 244, 20260928, 14, "hash").is_err());
        row["timing_ms"]["target_pdp_charged"] = json!(100.0);
        assert!(metrics(&row, &json!([1, 2]), 41, 0, 244, 20260928, 14, "hash").is_err());
    }
}
