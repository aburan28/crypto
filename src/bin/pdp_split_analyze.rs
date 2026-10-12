//! `pdp_split_analyze` — the analysis of `research/pdp_split_dreg_20261011`
//! (`PROTOCOL.md`, "Pricing the PDP variants in S").
//!
//! ```sh
//! target/release/pdp_split_analyze \
//!   --runs research/pdp_split_dreg_20261011/runs \
//!   --rho-session research/pdp_split_dreg_20261011/sessions/rho-01 \
//!   --calibration docs/ic/calibration.json \
//!   [--replica research/pdp_split_dreg_20261011/runs/aborted-launch-1] \
//!   --out research/pdp_split_dreg_20261011/ANALYSIS.json
//! ```
//!
//! Reads only recorded values.  Per cell and (plan, g): completion, the
//! exhaustive cross-check, solving degrees, F4 word XORs and wall time.
//! Per complete (plan, g): the relation-phase projection in `S`, priced
//! at the pinned XOR ratio, against both rho arms of the matched session.
//! Prints a Markdown table on stdout.

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use serde_json::{json, Value};

fn read_jsonl(path: &Path) -> Vec<Value> {
    fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).expect("JSON line"))
        .collect()
}

fn jsonl_files(dir: &Path) -> Vec<std::path::PathBuf> {
    let mut v: Vec<_> = fs::read_dir(dir)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "jsonl"))
        .collect();
    v.sort();
    v
}

fn mean(xs: &[f64]) -> Option<f64> {
    (!xs.is_empty()).then(|| xs.iter().sum::<f64>() / xs.len() as f64)
}

fn u(v: &Value, k: &str) -> u64 {
    v[k].as_u64()
        .unwrap_or_else(|| panic!("missing integer {k}"))
}

/// Mean `S` and mean online wall of each (arm, field degree) in a session.
fn rho_summary(session: &Path) -> BTreeMap<(String, u64), (f64, f64, usize)> {
    let mut acc: BTreeMap<(String, u64), (Vec<f64>, Vec<f64>)> = BTreeMap::new();
    for r in read_jsonl(&session.join("records.jsonl")) {
        if r["warmup"].as_bool() == Some(true) || r["outcome"]["status"] != "verified" {
            continue;
        }
        let key = (
            r["arm"].as_str().expect("arm").to_string(),
            r["workload"]["curve"]["field_degree"]
                .as_u64()
                .expect("field_degree"),
        );
        let e = acc.entry(key).or_default();
        e.0.push(r["cost"]["s"].as_f64().expect("cost.s"));
        e.1.push(r["online"]["wall_ns"].as_f64().unwrap_or(f64::NAN));
    }
    acc.into_iter()
        .map(|(k, (s, w))| {
            let n = s.len();
            (k, (mean(&s).unwrap(), mean(&w).unwrap(), n))
        })
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let flag = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let runs_dir = flag("--runs").expect("--runs DIR");
    let session = flag("--rho-session").expect("--rho-session DIR");
    let calibration_path = flag("--calibration").expect("--calibration FILE");
    let out_path = flag("--out").expect("--out FILE");
    let calibration: Value =
        serde_json::from_str(&fs::read_to_string(&calibration_path).expect("calibration"))
            .expect("calibration JSON");
    let rho = rho_summary(Path::new(&session));

    // Deterministic fields of a replica launch, keyed like the main runs.
    let mut replica: BTreeMap<String, (u64, u64, usize)> = BTreeMap::new();
    if let Some(dir) = flag("--replica") {
        for f in jsonl_files(Path::new(&dir)) {
            for r in read_jsonl(&f).into_iter().filter(|r| r["kind"] == "run") {
                let key = format!(
                    "{}:{}:{}:{}:{}:{}",
                    r["m"], r["n"], r["ell"], r["target"], r["plan"], r["g"]
                );
                replica.insert(
                    key,
                    (
                        u(&r, "solving_degree_max"),
                        u(&r, "word_xors"),
                        r["verified"].as_array().map_or(0, |v| v.len()),
                    ),
                );
            }
        }
    }
    let (mut replica_compared, mut replica_mismatch) = (0usize, 0usize);

    let mut cells_out = Vec::new();
    let mut table = String::from(
        "| cell m:a:n:ℓ | plan | g | unknowns | targets complete | exact | solving degree (per target) | mean word XORs | mean F4 wall s | λ | projected S (relations) | ÷ strong rho S | ÷ signed-Frobenius rho S |\n\
         |---|---|--:|--:|--:|---|---|--:|--:|--:|--:|--:|--:|\n",
    );
    for f in jsonl_files(Path::new(&runs_dir)) {
        let lines = read_jsonl(&f);
        let Some(cell) = lines.iter().find(|l| l["kind"] == "cell") else {
            continue;
        };
        let (m, a, n, ell) = (u(cell, "m"), u(cell, "a"), u(cell, "n"), u(cell, "ell"));
        let r = u(cell, "r") as f64;
        let fb = u(cell, "factor_base_liftable_abscissae") as f64;
        let lambda =
            u(cell, "lambda_decompositions") as f64 / u(cell, "lambda_targets").max(1) as f64;
        let cal_key = format!("koblitz/K_{a} / GF(2^{n})");
        let ratio = calibration["instances"][&cal_key]["ns_per_word_xor"].as_f64();
        let strong = rho.get(&("rho-strong".to_string(), n)).copied();
        let frob = rho.get(&("rho-frob".to_string(), n)).copied();

        let mut groups: BTreeMap<(String, u64), Vec<&Value>> = BTreeMap::new();
        for run in lines.iter().filter(|l| l["kind"] == "run") {
            let key = format!(
                "{}:{}:{}:{}:{}:{}",
                run["m"], run["n"], run["ell"], run["target"], run["plan"], run["g"]
            );
            if let Some(&(sd, x, v)) = replica.get(&key) {
                replica_compared += 1;
                let same = sd == u(run, "solving_degree_max")
                    && x == u(run, "word_xors")
                    && v == run["verified"].as_array().map_or(0, |v| v.len());
                replica_mismatch += usize::from(!same);
            }
            groups
                .entry((run["plan"].as_str().unwrap().to_string(), u(run, "g")))
                .or_default()
                .push(run);
        }
        let mut variants = Vec::new();
        for ((plan, g), runs) in &groups {
            let complete: Vec<&&Value> = runs
                .iter()
                .filter(|r| !r["complete_and_exact"].is_null())
                .collect();
            let exact_fail = complete
                .iter()
                .filter(|r| r["complete_and_exact"] == false)
                .count();
            let sds: Vec<u64> = runs.iter().map(|r| u(r, "solving_degree_max")).collect();
            let xors: Vec<f64> = runs.iter().map(|r| u(r, "word_xors") as f64).collect();
            let walls: Vec<f64> = runs
                .iter()
                .map(|r| u(r, "f4_wall_ns") as f64 * 1e-9)
                .collect();
            let all_complete = complete.len() == runs.len() && exact_fail == 0;
            let mean_x = mean(&xors).unwrap_or(0.0);
            let mean_w = mean(&walls).unwrap_or(0.0);
            // PROTOCOL.md: calls ≈ (F + 1)/λ; S ≈ calls · X̄ · ratio / √r.
            let projection = match (all_complete, ratio, lambda > 0.0) {
                (true, Some(q), true) => {
                    let calls = (fb + 1.0) / lambda;
                    let s = calls * mean_x * q / r.sqrt();
                    Some(json!({
                        "pdp_calls": calls,
                        "word_xors": calls * mean_x,
                        "gae": calls * mean_x * q,
                        "s_relations": s,
                        "over_strong_rho_s": strong.map(|v| s / v.0),
                        "over_signed_frobenius_rho_s": frob.map(|v| s / v.0),
                        "relation_phase_wall_s_descriptive": calls * mean_w,
                    }))
                }
                _ => None,
            };
            let n_vars = u(runs[0], "n_vars");
            let fmt = |x: Option<f64>| x.map_or("—".into(), |v| format!("{v:.3}"));
            let pj = projection.as_ref();
            let s_proj = pj.and_then(|p| p["s_relations"].as_f64());
            table.push_str(&format!(
                "| {m}:{a}:{n}:{ell} | {plan} | {g} | {n_vars} | {}/{} | {} | {} | {:.3e} | {:.3} | {:.4} | {} | {} | {} |\n",
                complete.len(),
                runs.len(),
                if exact_fail == 0 { "yes" } else { "**NO**" },
                sds.iter().map(|d| d.to_string()).collect::<Vec<_>>().join(", "),
                mean_x,
                mean_w,
                lambda,
                s_proj.map_or("—".into(), |s| format!("{s:.4e}")),
                fmt(pj.and_then(|p| p["over_strong_rho_s"].as_f64())),
                fmt(pj.and_then(|p| p["over_signed_frobenius_rho_s"].as_f64())),
            ));
            variants.push(json!({
                "plan": plan, "g": g, "n_vars": n_vars,
                "input_degree": u(runs[0], "input_degree"),
                "targets": runs.len(), "complete": complete.len(), "exact_failures": exact_fail,
                "solving_degree_per_target": sds,
                "word_xors_per_target": xors,
                "f4_wall_s_per_target_descriptive": walls,
                "verified_per_target": runs.iter().map(|r| r["verified"].as_array().map_or(0, |v| v.len())).collect::<Vec<_>>(),
                "truth_per_target": runs.iter().map(|r| u(r, "truth_count")).collect::<Vec<_>>(),
                "unliftable_solutions": runs.iter().map(|r| u(r, "unliftable_solutions")).sum::<u64>(),
                "false_solutions": runs.iter().map(|r| u(r, "false_solutions")).sum::<u64>(),
                "branches_timed_out": runs.iter().map(|r| u(r, "branches_timed_out")).sum::<u64>(),
                "mean_word_xors": mean_x,
                "mean_f4_wall_s_descriptive": mean_w,
                "projection": projection,
            }));
        }
        cells_out.push(json!({
            "file": f.file_name().unwrap().to_string_lossy(),
            "m": m, "a": a, "n": n, "ell": ell,
            "curve_slug": cell["curve_slug"], "r": r, "factor_base": fb,
            "surplus": cell["surplus_n_minus_m_ell"],
            "lambda": lambda, "lambda_targets": cell["lambda_targets"],
            "pinned_ns_per_word_xor": ratio, "calibration_key": cal_key,
            "rho_strong": strong.map(|(s, w, k)| json!({"mean_s": s, "mean_online_wall_ns": w, "runs": k})),
            "rho_signed_frobenius": frob.map(|(s, w, k)| json!({"mean_s": s, "mean_online_wall_ns": w, "runs": k})),
            "variants": variants,
        }));
    }
    let analysis = json!({
        "schema": "pdp_split_dreg.analysis/v1",
        "inputs": {"runs": runs_dir, "rho_session": session, "calibration": calibration_path, "replica": flag("--replica")},
        "projection_rule": "S_relations = (F+1)/lambda * mean_word_xors * ns_per_word_xor / sqrt(r); relation phase only; pinned ratio from another host",
        "replica_runs_compared": replica_compared,
        "replica_runs_mismatched": replica_mismatch,
        "cells": cells_out,
    });
    fs::write(
        &out_path,
        serde_json::to_string_pretty(&analysis).unwrap() + "\n",
    )
    .expect("write");
    print!("{table}");
    println!("\nreplica: {replica_compared} runs compared, {replica_mismatch} mismatched (solving degree, word XORs, verified count)");
}
