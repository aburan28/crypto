//! Round-2 analysis for the drastic prime conductor-gap study
//! (research/drastic_prime_conductor_gap_round2_20261012). Reads sealed
//! ecbench records; never launches a solver. Same statistical core as
//! round 1 (20 000 resamples, frozen seed 20261012). Two quotients per
//! pair: the same-method control (rho.negation on both curves) and the
//! best-method discrepancy (rho.negation on the floor over rho.cubic on
//! the j = 0 crater), plus the within-curve crater quotient and the
//! round-1 lead replication. Wall-time quotients are descriptive only.
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::Path;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ecbench::runner::read_records;
use crypto_lib::cryptanalysis::ecbench::stats::{mean, quantile, Resampler};
use serde_json::{json, Value};

const SAMPLES: usize = 20_000;
const SEED: u64 = 20261012;

fn ratio_of(a: &[f64], b: &[f64]) -> Value {
    // quotient mean(b) / mean(a)
    let ratio = mean(b).unwrap() / mean(a).unwrap();
    let mut rng = Resampler::new(SEED);
    let mut boot = Vec::with_capacity(SAMPLES);
    for _ in 0..SAMPLES {
        let sa: f64 = (0..a.len()).map(|_| a[rng.below(a.len())]).sum();
        let sb: f64 = (0..b.len()).map(|_| b[rng.below(b.len())]).sum();
        boot.push((sb / b.len() as f64) / (sa / a.len() as f64));
    }
    json!({
        "quotient": ratio,
        "ci95": [quantile(&boot, 0.025), quantile(&boot, 0.975)],
        "ci99": [quantile(&boot, 0.005), quantile(&boot, 0.995)],
        "resamples": SAMPLES, "seed": SEED,
        "method": "independent workload bootstrap on each fixed curve",
    })
}

fn arg(args: &[String], name: &str) -> String {
    let i = args.iter().position(|s| s == name).expect("required argument");
    args.get(i + 1).expect("argument value").clone()
}

/// (arm, per-workload S means, per-workload wall ms, slug) for one
/// session, split by curve when the session holds two.
fn samples(
    dir: &Path,
    arms: &[&str],
) -> Result<Vec<Value>, String> {
    let spec: Value = serde_json::from_slice(
        &std::fs::read(dir.join("spec.json")).map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;
    let curves = spec["workloads"]["curves"].as_array().ok_or("curves")?;
    let records = read_records(dir)?;
    let measured: Vec<_> = records.iter().filter(|r| !r.warmup).collect();
    let failed = measured
        .iter()
        .filter(|r| r.outcome.status != "verified")
        .map(|r| json!({"run_id": r.run_id, "arm": r.arm, "status": r.outcome.status}))
        .collect::<Vec<_>>();
    let mut out = Vec::new();
    for (ci, c) in curves.iter().enumerate() {
        for arm in arms {
            let mut workloads = BTreeMap::<String, Vec<f64>>::new();
            let mut walls = BTreeMap::<String, Vec<f64>>::new();
            let mut statuses = BTreeMap::<String, usize>::new();
            let mut count = 0usize;
            for r in measured.iter() {
                if r.arm != *arm
                    || serde_json::to_value(&r.workload.curve_spec).unwrap() != *c
                {
                    continue;
                }
                count += 1;
                *statuses.entry(r.outcome.status.clone()).or_default() += 1;
                if r.outcome.status == "verified" {
                    let s = r.cost.s.ok_or("verified run has no S")?;
                    workloads
                        .entry(r.workload.workload_id.clone())
                        .or_default()
                        .push(s);
                    let wall = r.time.solve_wall_ns.unwrap_or(0) as f64 / 1e6;
                    walls.entry(r.workload.workload_id.clone()).or_default().push(wall);
                }
            }
            let s_vals = workloads.values().map(|v| mean(v).unwrap()).collect::<Vec<f64>>();
            let w_vals = walls.values().map(|v| mean(v).unwrap()).collect::<Vec<f64>>();
            out.push(json!({
                "curve_index": ci, "arm": arm,
                "slug": records.iter().find(|r| serde_json::to_value(&r.workload.curve_spec).unwrap() == *c)
                    .map(|r| r.workload.curve.slug.clone()),
                "measured": count, "statuses": statuses,
                "valid": count == 32 && workloads.len() == 32
                    && statuses.get("verified") == Some(&count),
                "s_per_workload": s_vals, "wall_ms_per_workload": w_vals,
            }));
        }
    }
    let source = std::fs::read(dir.join("records.jsonl")).map_err(|e| e.to_string())?;
    let _ = sha256_hex(&source);
    Ok(out
        .into_iter()
        .chain(std::iter::once(json!({"failed": failed})))
        .collect())
}

fn aa_check(dir: &Path) -> Result<Value, String> {
    let records = read_records(dir)?;
    let measured: Vec<_> = records.iter().filter(|r| !r.warmup).collect();
    let mut pairs = 0usize;
    let mut mismatch = Vec::new();
    for r in measured.iter().filter(|r| r.arm == "rho-neg") {
        let controls: Vec<_> = measured
            .iter()
            .filter(|x| {
                x.arm == "rho-neg-aa"
                    && x.workload.workload_id == r.workload.workload_id
                    && x.round == r.round
            })
            .collect();
        pairs += 1;
        if controls.len() != 1
            || controls[0].cost != r.cost
            || controls[0].counters != r.counters
            || controls[0].outcome.recovered != r.outcome.recovered
        {
            mismatch.push(r.run_id.clone());
        }
    }
    Ok(json!({"pairs": pairs, "mismatches": mismatch}))
}

fn pick(rows: &[Value], curve_index: usize, arm: &str) -> Option<Vec<f64>> {
    rows.iter()
        .find(|r| r["curve_index"] == curve_index && r["arm"] == arm)
        .and_then(|r| {
            (r["valid"] == true)
                .then(|| r["s_per_workload"].as_array().cloned())
                .flatten()
        })
        .map(|v| v.iter().filter_map(|x| x.as_f64()).collect::<Vec<f64>>())
}

fn pick_wall(rows: &[Value], curve_index: usize, arm: &str) -> Option<Vec<f64>> {
    rows.iter()
        .find(|r| r["curve_index"] == curve_index && r["arm"] == arm)
        .and_then(|r| r["wall_ms_per_workload"].as_array().cloned())
        .map(|v| v.iter().filter_map(|x| x.as_f64()).collect::<Vec<f64>>())
}

fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    let root = arg(&args, "--sessions");
    let out = arg(&args, "--out");
    let cells = [
        ("d3-ell1016839-rho", 1016839u64),
        ("d3-ell1016839-cubic", 1016839),
        ("d3-ell2033609-rho", 2033609),
        ("d3-ell2033609-cubic", 2033609),
        ("replicate-d43-ell1016839", 0),
    ];
    let mut sessions = BTreeMap::new();
    for (cell, _ell) in cells {
        let dir = Path::new(&root).join(cell);
        let arms: &[&str] = if cell.ends_with("cubic") {
            &["rho-cubic"]
        } else {
            &["rho-neg", "rho-neg-aa"]
        };
        let rows = samples(&dir, arms)?;
        let aa = if cell.ends_with("rho") || cell.starts_with("replicate") {
            aa_check(&dir)?
        } else {
            json!({"pairs": 0, "mismatches": []})
        };
        let failed = rows
            .iter()
            .find(|r| r.get("failed").is_some())
            .and_then(|r| r["failed"].as_array().cloned())
            .unwrap_or_default();
        sessions.insert(
            cell.to_string(),
            json!({"rows": rows, "aa": aa, "failed": failed,
                   "records_dir": cell}),
        );
    }

    let mut pairs = Vec::new();
    for ell in [1016839u64, 2033609u64] {
        let rho = sessions[&format!("d3-ell{ell}-rho")]["rows"].as_array().unwrap();
        let cubic = sessions[&format!("d3-ell{ell}-cubic")]["rows"].as_array().unwrap();
        let crater_neg = pick(rho, 0, "rho-neg");
        let floor_neg = pick(rho, 1, "rho-neg");
        let crater_cubic = pick(cubic, 0, "rho-cubic");
        let crater_neg_wall = pick_wall(rho, 0, "rho-neg");
        let floor_neg_wall = pick_wall(rho, 1, "rho-neg");
        let crater_cubic_wall = pick_wall(cubic, 0, "rho-cubic");
        let same = match (&crater_neg, &floor_neg) {
            (Some(a), Some(b)) => ratio_of(a, b),
            _ => json!({"status": "inadmissible"}),
        };
        let best = match (&crater_cubic, &floor_neg) {
            (Some(c), Some(f)) => ratio_of(c, f),
            _ => json!({"status": "inadmissible"}),
        };
        let within = match (&crater_cubic, &crater_neg) {
            (Some(c), Some(n)) => ratio_of(c, n),
            _ => json!({"status": "inadmissible"}),
        };
        let best_wall = match (&crater_cubic_wall, &floor_neg_wall) {
            (Some(c), Some(f)) => ratio_of(c, f),
            _ => json!({"status": "inadmissible"}),
        };
        pairs.push(json!({
            "ell": ell, "d": 3,
            "same_method_rho_neg": same,
            "best_method_floor_neg_over_crater_cubic": best,
            "within_crater_neg_over_cubic": within,
            "best_method_wall_descriptive": best_wall,
            "crater_wall_neg_descriptive": crater_neg_wall.map(|w| mean(&w)),
            "crater_wall_cubic_descriptive": crater_cubic_wall.map(|w| mean(&w)),
            "aa": sessions[&format!("d3-ell{ell}-rho")]["aa"],
            "failed": sessions[&format!("d3-ell{ell}-rho")]["failed"],
        }));
    }

    // round-1 lead replication
    let repl = sessions["replicate-d43-ell1016839"]["rows"].as_array().unwrap();
    let replication = match (pick(repl, 0, "rho-neg"), pick(repl, 1, "rho-neg")) {
        (Some(a), Some(b)) => json!({
            "same_method_rho_neg": ratio_of(&a, &b),
            "aa": sessions["replicate-d43-ell1016839"]["aa"],
            "failed": sessions["replicate-d43-ell1016839"]["failed"],
        }),
        _ => json!({"status": "inadmissible"}),
    };

    let doc = json!({
        "schema": "drastic-prime-conductor-gap-round2-analysis/v1",
        "seed": SEED, "resamples": SAMPLES,
        "pairs": pairs,
        "replication_d43_ell1016839": replication,
        "status": "count diagnostics on fixed pairs; frozen rules in PROTOCOL.md",
    });
    std::fs::create_dir_all(&out).unwrap();
    std::fs::write(
        Path::new(&out).join("analysis.json"),
        serde_json::to_vec_pretty(&doc).unwrap(),
    )
    .unwrap();

    let mut md = String::from(
        "# Drastic prime conductor gap, round 2: the j = 0 crater\n\n\
`S` = counted GAE / sqrt(r). The best-method quotient is `S(rho.negation, floor) / S(rho.cubic, crater)`; the predicted value is sqrt(3) = 1.7321. Wall quotients are descriptive (L0 host).\n\n\
| pair | same-method (99%) | best-method (99%) | within crater (99%) | wall best-method (95%) |\n|---|---|---|---|---|\n",
    );
    for p in &pairs {
        let cell = |v: &Value, k: &str| -> String {
            let q = &v[k];
            if q.get("quotient").is_none() {
                return "inadmissible".into();
            }
            format!(
                "{:.4} [{:.4} to {:.4}]",
                q["quotient"].as_f64().unwrap(),
                q["ci99"][0].as_f64().unwrap(),
                q["ci99"][1].as_f64().unwrap()
            )
        };
        let _ = writeln!(
            md,
            "| D=-3 ell={} | {} | {} | {} | {} |",
            p["ell"],
            cell(p, "same_method_rho_neg"),
            cell(p, "best_method_floor_neg_over_crater_cubic"),
            cell(p, "within_crater_neg_over_cubic"),
            cell(p, "best_method_wall_descriptive"),
        );
    }
    let repl_q = &replication["same_method_rho_neg"];
    if repl_q.get("quotient").is_some() {
        let _ = writeln!(
            md,
            "\nRound-1 lead replication (d43-ell1016839, 32 fresh targets): floor/crater {:.6}, 99% [{:.6} to {:.6}].",
            repl_q["quotient"].as_f64().unwrap(),
            repl_q["ci99"][0].as_f64().unwrap(),
            repl_q["ci99"][1].as_f64().unwrap()
        );
    }
    std::fs::write(Path::new(&out).join("analysis.md"), md).unwrap();
    println!("wrote {}/analysis.json and analysis.md", out);
    Ok(())
}
