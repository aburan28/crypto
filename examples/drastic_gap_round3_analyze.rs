//! Round-3 paired analysis (research/drastic_prime_conductor_gap_round3_20261012).
//! Reads sealed ecbench records; never launches a solver. The paired
//! statistic divides per-index costs of the SAME planted scalar on both
//! curves (target law uniform_scalar_sha256_shared_v1) and bootstraps
//! indices jointly, removing between-target variance. The independent
//! quotient of the prior rounds is reported beside it.
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::Path;

use crypto_lib::cryptanalysis::ecbench::runner::read_records;
use crypto_lib::cryptanalysis::ecbench::stats::{mean, quantile, Resampler};
use serde_json::{json, Value};

const SAMPLES: usize = 20_000;
const SEED: u64 = 20261012;

fn log_paired(a: &[f64], b: &[f64]) -> Value {
    // b/a per index (same target index = same scalar), mean of logs,
    // bootstrap over indices jointly.
    let n = a.len().min(b.len());
    let d: Vec<f64> = (0..n).map(|i| (b[i] / a[i]).ln()).collect();
    let stat = d.iter().sum::<f64>() / n as f64;
    let mut rng = Resampler::new(SEED);
    let mut boot = Vec::with_capacity(SAMPLES);
    for _ in 0..SAMPLES {
        let mut s = 0.0;
        for _ in 0..n {
            s += d[rng.below(n)];
        }
        boot.push((s / n as f64).exp());
    }
    json!({
        "quotient": stat.exp(),
        "ci95": [quantile(&boot, 0.025), quantile(&boot, 0.975)],
        "ci99": [quantile(&boot, 0.005), quantile(&boot, 0.995)],
        "pairs": n, "resamples": SAMPLES, "seed": SEED,
        "method": "paired per-index log quotient, joint bootstrap over indices",
    })
}

fn independent(a: &[f64], b: &[f64]) -> Value {
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
        "ci99": [quantile(&boot, 0.005), quantile(&boot, 0.995)],
        "method": "independent-workload bootstrap (the prior rounds' law)",
    })
}

fn arg(args: &[String], name: &str) -> String {
    let i = args
        .iter()
        .position(|s| s == name)
        .expect("required argument");
    args.get(i + 1).expect("argument value").clone()
}

/// Per-workload S keyed by target index (the planted k's index), for one
/// (session, arm, curve). Only measured, verified, shared-law records.
fn per_index(dir: &Path, arm: &str, curve_index: usize) -> Result<BTreeMap<u64, f64>, String> {
    let spec: Value =
        serde_json::from_slice(&std::fs::read(dir.join("spec.json")).map_err(|e| e.to_string())?)
            .map_err(|e| e.to_string())?;
    let curves = spec["workloads"]["curves"].as_array().ok_or("curves")?;
    if curves.len() <= curve_index {
        return Err("curve index out of range".into());
    }
    let shared = spec["workloads"]["target_kind"] == json!("planted_shared");
    if !shared {
        return Err("session is not planted_shared".into());
    }
    let records = read_records(dir)?;
    let mut out = BTreeMap::new();
    for r in records.iter().filter(|r| !r.warmup) {
        if r.arm != arm
            || serde_json::to_value(&r.workload.curve_spec).unwrap() != curves[curve_index]
        {
            continue;
        }
        if r.outcome.status != "verified" {
            return Err(format!("measured failure in {} arm {}", dir.display(), arm));
        }
        let s = r.cost.s.ok_or("verified run has no S")?;
        out.insert(r.workload.target_index, s);
    }
    Ok(out)
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

fn vec_of(m: &BTreeMap<u64, f64>) -> Vec<f64> {
    m.values().copied().collect()
}

fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    let root = arg(&args, "--sessions");
    let out = arg(&args, "--out");
    let mut results = Vec::new();

    // D = -3 pairs, re-measured paired: rho-neg both, rho-cubic crater
    for ell in [1016839u64, 2033609u64] {
        let rho_dir = Path::new(&root).join(format!("d3-ell{ell}-rho"));
        let cubic_dir = Path::new(&root).join(format!("d3-ell{ell}-cubic"));
        let crater = per_index(&rho_dir, "rho-neg", 0)?;
        let floor = per_index(&rho_dir, "rho-neg", 1)?;
        let cubic = per_index(&cubic_dir, "rho-cubic", 0)?;
        let (c, f, q) = (vec_of(&crater), vec_of(&floor), vec_of(&cubic));
        results.push(json!({
            "family": "d3", "ell": ell,
            "same_method_paired": log_paired(&c, &f),
            "same_method_independent": independent(&c, &f),
            "best_method_paired": log_paired(&q, &f),
            "within_crater_paired": log_paired(&q, &c),
            "aa": aa_check(&rho_dir)?,
            "prediction": "sqrt(3) = 1.7321",
        }));
    }

    // D = -4 pair: rho-neg both, rho-quartic crater
    {
        let doc: Value = serde_json::from_slice(
            &std::fs::read(Path::new(&root).join("..").join("pairs").join("pairs.json"))
                .map_err(|e| e.to_string())?,
        )
        .map_err(|e| e.to_string())?;
        let ell = doc["pairs"]
            .as_array()
            .and_then(|p| p.iter().find(|x| x["family"] == "d4"))
            .and_then(|x| x["ell"].as_u64())
            .ok_or("no d4 pair")?;
        let rho_dir = Path::new(&root).join(format!("d4-ell{ell}-rho"));
        let auto_dir = Path::new(&root).join(format!("d4-ell{ell}-auto"));
        let crater = per_index(&rho_dir, "rho-neg", 0)?;
        let floor = per_index(&rho_dir, "rho-neg", 1)?;
        let quartic = per_index(&auto_dir, "rho-quartic", 0)?;
        let (c, f, qq) = (vec_of(&crater), vec_of(&floor), vec_of(&quartic));
        results.push(json!({
            "family": "d4", "ell": ell,
            "same_method_paired": log_paired(&c, &f),
            "same_method_independent": independent(&c, &f),
            "best_method_paired": log_paired(&qq, &f),
            "within_crater_paired": log_paired(&qq, &c),
            "aa": aa_check(&rho_dir)?,
            "prediction": "sqrt(2) = 1.4142",
        }));
    }

    // D = -19 gap sweep
    {
        let doc: Value = serde_json::from_slice(
            &std::fs::read(Path::new(&root).join("..").join("pairs").join("pairs.json"))
                .map_err(|e| e.to_string())?,
        )
        .map_err(|e| e.to_string())?;
        let ell = doc["pairs"]
            .as_array()
            .and_then(|p| p.iter().find(|x| x["family"] == "d19"))
            .and_then(|x| x["ell"].as_u64())
            .ok_or("no d19 pair")?;
        let rho_dir = Path::new(&root).join(format!("d19-ell{ell}-rho"));
        let crater = per_index(&rho_dir, "rho-neg", 0)?;
        let floor = per_index(&rho_dir, "rho-neg", 1)?;
        let (c, f) = (vec_of(&crater), vec_of(&floor));
        results.push(json!({
            "family": "d19", "ell": ell,
            "same_method_paired": log_paired(&c, &f),
            "same_method_independent": independent(&c, &f),
            "aa": aa_check(&rho_dir)?,
            "prediction": "1 (the generic theorem at a 30x gap)",
        }));
    }

    let doc = json!({
        "schema": "drastic-prime-conductor-gap-round3-analysis/v1",
        "seed": SEED, "resamples": SAMPLES,
        "target_law": "uniform_scalar_sha256_shared_v1",
        "results": results,
        "status": "paired count diagnostics; frozen rules in PROTOCOL.md",
    });
    std::fs::create_dir_all(&out).unwrap();
    std::fs::write(
        Path::new(&out).join("analysis.json"),
        serde_json::to_vec_pretty(&doc).unwrap(),
    )
    .unwrap();
    let mut md = String::from(
        "# Drastic prime conductor gap, round 3: paired targets\n\n\
Paired planted-shared scalars (the same `k` on both curves); per-index log quotients, joint bootstrap over indices (20,000 resamples, frozen seed). The independent-law quotient is reported beside each paired one.\n\n\
| family | ell | prediction | same-method paired (99%) | best-method paired (99%) | same-method independent (99%) |\n|---|---|---|---|---|---|\n",
    );
    for r in &results {
        let q = |k: &str| -> String {
            let v = &r[k];
            if v.get("quotient").is_none() {
                return "inadmissible".into();
            }
            format!(
                "{:.4} [{:.4} to {:.4}]",
                v["quotient"].as_f64().unwrap(),
                v["ci99"][0].as_f64().unwrap(),
                v["ci99"][1].as_f64().unwrap()
            )
        };
        let _ = writeln!(
            md,
            "| {} | {} | {} | {} | {} | {} |",
            r["family"],
            r["ell"],
            r["prediction"],
            q("same_method_paired"),
            q("best_method_paired"),
            q("same_method_independent"),
        );
    }
    std::fs::write(Path::new(&out).join("analysis.md"), md).unwrap();
    println!("wrote {}/analysis.json and analysis.md", out);
    Ok(())
}
