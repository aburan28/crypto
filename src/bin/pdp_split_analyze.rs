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
    let runs_dir = flag("--runs").unwrap_or_default();
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
    let run_files_empty = runs_dir.is_empty();
    let run_files = if runs_dir.is_empty() {
        Vec::new()
    } else {
        jsonl_files(Path::new(&runs_dir))
    };
    for f in run_files {
        let lines = read_jsonl(&f);
        let Some(cell) = lines.iter().find(|l| l["kind"] == "cell") else {
            continue;
        };
        let (m, a, n, ell) = (u(cell, "m"), u(cell, "a"), u(cell, "n"), u(cell, "ell"));
        let r = u(cell, "r") as f64;
        // A closure cell (round 2, `--base closure`) needs its folded
        // columns and its per-system yield; a standard cell its base size
        // and per-target yield (one system per target).
        let closure = cell["base"] == "closure";
        let fb = if closure {
            u(cell, "closure_columns") as f64
        } else {
            u(cell, "factor_base_liftable_abscissae") as f64
        };
        let lambda = if closure {
            u(cell, "lambda_system_decompositions") as f64 / u(cell, "lambda_systems").max(1) as f64
        } else {
            u(cell, "lambda_decompositions") as f64 / u(cell, "lambda_targets").max(1) as f64
        };
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
                "| {m}:{a}:{n}:{ell}{} | {plan} | {g} | {n_vars} | {}/{} | {} | {} | {:.3e} | {:.3} | {:.4} | {} | {} | {} |\n",
                if closure { " closure" } else { "" },
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
            "m": m, "a": a, "n": n, "ell": ell, "base": cell["base"].as_str().unwrap_or("standard"),
            "curve_slug": cell["curve_slug"], "r": r, "factor_base": fb,
            "surplus": cell["surplus_n_minus_m_ell"],
            "lambda": lambda, "lambda_targets": cell["lambda_targets"],
            "pinned_ns_per_word_xor": ratio, "calibration_key": cal_key,
            "rho_strong": strong.map(|(s, w, k)| json!({"mean_s": s, "mean_online_wall_ns": w, "runs": k})),
            "rho_signed_frobenius": frob.map(|(s, w, k)| json!({"mean_s": s, "mean_online_wall_ns": w, "runs": k})),
            "variants": variants,
        }));
    }
    let ecbench_dirs: Vec<String> = args
        .iter()
        .enumerate()
        .filter(|(_, a)| *a == "--ecbench")
        .filter_map(|(i, _)| args.get(i + 1).cloned())
        .collect();
    let (ecbench_json, ecbench_table) = ecbench_summary(&ecbench_dirs, &calibration);
    let analysis = json!({
        "schema": "pdp_split_dreg.analysis/v1",
        "ecbench": ecbench_json,
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
    if !run_files_empty {
        print!("{table}");
    }
    print!("{ecbench_table}");
    println!("\nreplica: {replica_compared} runs compared, {replica_mismatch} mismatched (solving degree, word XORs, verified count)");
}

fn median(xs: &mut [f64]) -> Option<f64> {
    if xs.is_empty() {
        return None;
    }
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let k = xs.len();
    Some(if k % 2 == 1 {
        xs[k / 2]
    } else {
        (xs[k / 2 - 1] + xs[k / 2]) / 2.0
    })
}

/// Per (curve, arm) of ecbench sessions: S as recorded (a lower bound for
/// the IC arms: F4 work counted, not charged), F4 word XORs and systems
/// per run, and the pinned-ratio sensitivity `S + XORs · ns_per_word_xor /
/// √r` per run; per IC arm, the per-workload ratios of that priced S to
/// each rho arm's S on the same workload, and closure arms to their
/// standard twin (`clo-X` ↔ `std-X`).  Everything is read from records.
fn ecbench_summary(dirs: &[String], calibration: &Value) -> (Value, String) {
    // (n, arm) -> workload -> (s, priced, xors, systems, verified)
    type Row = (f64, f64, Option<f64>, Option<f64>, bool);
    let mut by: BTreeMap<(u64, String), BTreeMap<String, Row>> = BTreeMap::new();
    let mut meta: BTreeMap<u64, (String, f64, Option<f64>)> = BTreeMap::new();
    for d in dirs {
        for r in read_jsonl(&Path::new(d).join("records.jsonl")) {
            if r["warmup"].as_bool() == Some(true) {
                continue;
            }
            let w = &r["workload"];
            let n = w["curve"]["field_degree"].as_u64().expect("field_degree");
            let a = w["curve_spec"]["a"].as_u64().unwrap_or(0);
            let rr = w["curve"]["r"].as_f64().expect("r");
            let key = format!("koblitz/K_{a} / GF(2^{n})");
            let ratio = calibration["instances"][&key]["ns_per_word_xor"].as_f64();
            meta.insert(
                n,
                (
                    w["curve"]["slug"].as_str().unwrap_or("?").to_string(),
                    rr,
                    ratio,
                ),
            );
            let verified = r["outcome"]["status"] == "verified";
            let s = r["cost"]["s"].as_f64().unwrap_or(f64::NAN);
            let xors = r["solver"]["ops"].as_f64();
            let systems = r["solver"]["extra"]["systems"].as_f64();
            let priced = match (xors, ratio) {
                (Some(x), Some(q)) => s + x * q / rr.sqrt(),
                (None, _) => s,
                _ => f64::NAN,
            };
            let arm = r["arm"].as_str().expect("arm").to_string();
            let wid = w["workload_id"].as_str().expect("workload_id").to_string();
            by.entry((n, arm))
                .or_default()
                .insert(wid, (s, priced, xors, systems, verified));
        }
    }
    let mut out = Vec::new();
    let mut table = String::from(
        "\n| curve | arm | verified | median S (recorded) | median F4 XORs | median systems | median priced S | priced ÷ strong rho S, median (min–max) | priced ÷ signed-Frobenius rho S, median | priced ÷ standard twin, median |\n\
         |---|---|--:|--:|--:|--:|--:|---|--:|--:|\n",
    );
    for ((n, arm), rows) in &by {
        let (slug, _, ratio) = meta[n].clone();
        let col = |f: &dyn Fn(&Row) -> Option<f64>| -> Vec<f64> {
            rows.values()
                .filter_map(f)
                .filter(|v| v.is_finite())
                .collect()
        };
        let verified = rows.values().filter(|r| r.4).count();
        let paired = |other: &str| -> Vec<f64> {
            let Some(o) = by.get(&(*n, other.to_string())) else {
                return Vec::new();
            };
            rows.iter()
                .filter_map(|(w, r)| o.get(w).map(|q| r.1 / q.1))
                .filter(|v| v.is_finite())
                .collect()
        };
        let mut vs_strong = paired("rho-strong");
        let mut vs_frob = paired("rho-frob");
        let twin = arm.strip_prefix("clo-").map(|t| format!("std-{t}"));
        let mut vs_twin = twin.as_deref().map(paired).unwrap_or_default();
        let range = |v: &[f64]| {
            let lo = v.iter().cloned().fold(f64::INFINITY, f64::min);
            let hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            (lo, hi)
        };
        let (lo, hi) = range(&vs_strong);
        let ms = median(&mut col(&|r| Some(r.0)));
        let mx = median(&mut col(&|r| r.2));
        let msys = median(&mut col(&|r| r.3));
        let mp = median(&mut col(&|r| Some(r.1)));
        let m_strong = median(&mut vs_strong.clone());
        let m_frob = median(&mut vs_frob);
        let m_twin = median(&mut vs_twin);
        let f = |x: Option<f64>| x.map_or("—".to_string(), |v| format!("{v:.4}"));
        let e = |x: Option<f64>| x.map_or("—".to_string(), |v| format!("{v:.3e}"));
        table.push_str(&format!(
            "| {slug} | {arm} | {verified}/{} | {} | {} | {} | {} | {} ({}) | {} | {} |\n",
            rows.len(),
            f(ms),
            e(mx),
            f(msys),
            f(mp),
            f(m_strong),
            if vs_strong.is_empty() {
                "—".to_string()
            } else {
                format!("{lo:.3}–{hi:.3}")
            },
            f(m_frob),
            f(m_twin),
        ));
        vs_strong.sort_by(|a, b| a.partial_cmp(b).unwrap());
        out.push(json!({
            "n": n, "curve_slug": slug, "arm": arm, "pinned_ns_per_word_xor": ratio,
            "runs": rows.len(), "verified": verified,
            "median_s_recorded": ms, "median_f4_word_xors": mx, "median_systems": msys,
            "median_priced_s": mp,
            "priced_over_strong_rho_per_workload": vs_strong,
            "median_priced_over_strong_rho": m_strong,
            "median_priced_over_signed_frobenius_rho": m_frob,
            "median_priced_over_standard_twin": m_twin,
            "per_workload": rows.iter().map(|(w, r)| json!({
                "workload": w, "s": r.0, "priced_s": r.1, "f4_word_xors": r.2,
                "systems": r.3, "verified": r.4})).collect::<Vec<_>>(),
        }));
    }
    (
        json!({
            "sessions": dirs,
            "pricing": "priced S = recorded S + F4 word XORs * pinned ns_per_word_xor (docs/ic/calibration.json, another host) / sqrt(r); a sensitivity, not the recorded S",
            "rows": out,
        }),
        if dirs.is_empty() {
            String::new()
        } else {
            table
        },
    )
}
