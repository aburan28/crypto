//! Rebuild the preregistered m=83 step-rate summary from sealed ecbench sessions.
//! This is a stage diagnostic, never a complete-DLP comparison.

use std::collections::BTreeMap;
use std::path::Path;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ecbench::record::Record;
use crypto_lib::cryptanalysis::ecbench::runner::read_records;
use serde_json::json;

fn measured(dir: &Path) -> Result<Vec<Record>, String> {
    let records = read_records(dir)?;
    let measured: Vec<_> = records.into_iter().filter(|r| !r.warmup).collect();
    if measured.len() != 6 {
        return Err(format!("{}: expected six measured records", dir.display()));
    }
    Ok(measured)
}

fn rate(r: &Record) -> Result<f64, String> {
    let ns = r.time.solve_wall_ns.ok_or("missing solve wall time")?;
    let steps = *r
        .counters
        .get("walk_operations")
        .ok_or("missing walk_operations")?;
    if steps == 0 {
        return Err("zero walk operations".into());
    }
    Ok(ns as f64 / steps as f64)
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_by(f64::total_cmp);
    (values[values.len() / 2 - 1] + values[values.len() / 2]) / 2.0
}

fn key(r: &Record) -> (String, u32) {
    (r.workload.workload_id.clone(), r.round)
}

fn digest(dir: &Path) -> Result<String, String> {
    let bytes = std::fs::read(dir.join("records.jsonl")).map_err(|e| e.to_string())?;
    Ok(sha256_hex(&bytes))
}

fn run(root: &Path) -> Result<(), String> {
    let path = |s: &str| root.join("sessions").join(s);
    let narrow = measured(&path("n61-narrow"))?;
    let wide = measured(&path("n61-wide"))?;
    let gate = measured(&path("m83"))?;
    let by_key: BTreeMap<_, _> = wide.iter().map(|r| (key(r), r)).collect();
    let mut narrow_ns = Vec::new();
    let mut wide_ns = Vec::new();
    for a in &narrow {
        let b = by_key.get(&key(a)).ok_or("n61 pair missing")?;
        if a.algorithm_seed != b.algorithm_seed
            || a.outcome.status != "verified"
            || b.outcome.status != "verified"
            || a.outcome.recovered != b.outcome.recovered
            || a.cost.total_gae != b.cost.total_gae
            || a.counters != b.counters
            || a.phases
                .iter()
                .map(|p| (p.adds, p.doubles, p.scalar_mults))
                .collect::<Vec<_>>()
                != b.phases
                    .iter()
                    .map(|p| (p.adds, p.doubles, p.scalar_mults))
                    .collect::<Vec<_>>()
        {
            return Err(format!("n61 narrow/wide differ at {:?}", key(a)));
        }
        narrow_ns.push(rate(a)?);
        wide_ns.push(rate(b)?);
    }
    let mut gate_ns = Vec::new();
    for r in &gate {
        if r.outcome.status != "exhausted" || r.counters.get("walk_operations") != Some(&20_000_000)
        {
            return Err(format!("m83 record {} did not reach the budget", r.seq));
        }
        gate_ns.push(rate(r)?);
    }
    let n_med = median(narrow_ns.clone());
    let w_med = median(wide_ns.clone());
    let g_med = median(gate_ns.clone());
    let floor_steps =
        (std::f64::consts::PI * 2_417_851_639_230_796_216_685_689_f64 / (2.0 * 166.0)).sqrt();
    let projection_hours = |overhead: f64| floor_steps * overhead * g_med / 1e9 / 3600.0;
    let report = json!({
        "schema": "ecbench.gate_table/v1",
        "source_records_sha256": {
            "n61_narrow": digest(&path("n61-narrow"))?,
            "n61_wide": digest(&path("n61-wide"))?,
            "m83": digest(&path("m83"))?
        },
        "measured_runs_per_arm": 6,
        "n61_step_for_step_identical": true,
        "median_ns_per_walk_operation": {
            "n61_narrow": n_med,
            "n61_wide": w_med,
            "m83": g_med
        },
        "per_run_ns_per_walk_operation": {
            "n61_narrow": narrow_ns,
            "n61_wide": wide_ns,
            "m83": gate_ns
        },
        "wide_over_narrow_n61": w_med / n_med,
        "m83_over_n61_wide": g_med / w_med,
        "gate_floor_steps": floor_steps,
        "illustrative_single_core_hours_at_m83": {
            "floor": projection_hours(1.0),
            "overhead_3_percent": projection_hours(1.03),
            "overhead_8_percent": projection_hours(1.08)
        },
        "wall_level": "L0",
        "classification": "stage_diagnostic",
        "note": "L0 wall times and the hours projection are descriptive; no full m83 solve or wall-speed claim follows"
    });
    println!(
        "{}",
        serde_json::to_string_pretty(&report).map_err(|e| e.to_string())?
    );
    Ok(())
}

fn main() {
    let root = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "research/ecbench_m83_gate_20261003".into());
    if let Err(e) = run(Path::new(&root)) {
        eprintln!("ecbench_gate_table: {e}");
        std::process::exit(1);
    }
}
