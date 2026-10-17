//! Analysis of the F4/F5/F6 cross-family round
//! (`research/f4_f5_f6_cross_family_20261011/PROTOCOL.md`).
//!
//! Reads only the sessions' `records.jsonl` and applies the protocol's
//! decision rules as written:
//!
//! 1. F5 against F4 (matrix) at each cell: `W(matrix-F5)/W(matrix-F4)`.
//! 2. Inherited against matrix F4: `W(matrix-F4)/W(inherited-F4)`.
//! 3. F6-IC against inherited F4: `log2(W(inherited)/W(F6-IC))`.
//! 4. Every IC arm's lower-bound `S` against each same-target rho arm.
//! 5. The capped matrix arms at m = 19: outcome and counts per run.
//!
//! `W` is the record's counted solver work (`solver.ops`, word XORs). A
//! cell whose IC run fails verification is reported as failed and enters
//! no ratio; a cell with more than half its targets timed out gets no
//! median.
//!
//! Usage: `f4f5f6_cross_family_analyze OUT.json SESSION_DIR...`

use serde_json::{json, Value};
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
struct Run {
    arm: String,
    workload: String,
    status: String,
    slug: String,
    family: String,
    m: u64,
    r: f64,
    s: Option<f64>,
    word_xors: Option<u64>,
    calls: Option<u64>,
    reductions: Option<u64>,
    geometry: Option<u64>,
    relations: Option<u64>,
    trials: Option<u64>,
    relations_wall_ns: Option<u64>,
    base_points: Option<u64>,
    isolation: String,
}

fn read(dir: &str) -> Result<(String, Vec<Run>), String> {
    let path = format!("{dir}/records.jsonl");
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{path}: {e}"))?;
    let session: Value = std::fs::read_to_string(format!("{dir}/session.json"))
        .ok()
        .and_then(|t| serde_json::from_str(&t).ok())
        .unwrap_or(Value::Null);
    let id = session["session_id"].as_str().unwrap_or("?").to_string();
    let mut runs = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let r: Value = serde_json::from_str(line).map_err(|e| format!("{path}: {e}"))?;
        if r["warmup"].as_bool().unwrap_or(false) {
            continue;
        }
        let curve = &r["workload"]["curve"];
        let phase = |name: &str| {
            r["phases"]
                .as_array()
                .and_then(|ps| ps.iter().find(|p| p["name"] == name).cloned())
        };
        runs.push(Run {
            arm: r["arm"].as_str().unwrap_or("?").to_string(),
            workload: r["workload"]["workload_id"]
                .as_str()
                .unwrap_or("?")
                .to_string(),
            status: r["outcome"]["status"].as_str().unwrap_or("?").to_string(),
            slug: curve["slug"].as_str().unwrap_or("?").to_string(),
            family: curve["family"].as_str().unwrap_or("?").to_string(),
            m: curve["field_degree"].as_u64().unwrap_or(0),
            r: curve["r"].as_f64().unwrap_or(f64::NAN),
            s: r["cost"]["s"].as_f64(),
            word_xors: r["solver"]["ops"].as_u64(),
            calls: r["solver"]["calls"].as_u64(),
            reductions: r["solver"]["extra"]["reductions"].as_u64(),
            geometry: r["solver"]["extra"]["geometric_group_additions"].as_u64(),
            relations: r["counters"]["relations_found"].as_u64(),
            trials: r["counters"]["targets_tried"].as_u64(),
            relations_wall_ns: phase("relations").and_then(|p| p["wall_ns"].as_u64()),
            base_points: r["factor_base"]["signed_points"].as_u64(),
            isolation: r["isolation"]["level"].as_str().unwrap_or("?").to_string(),
        });
    }
    Ok((id, runs))
}

fn median(v: &mut [f64]) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    let n = v.len();
    Some(if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    })
}

fn summary(v: &[f64]) -> Value {
    let mut w = v.to_vec();
    let med = median(&mut w);
    json!({
        "n": v.len(),
        "median": med,
        "min": w.first(),
        "max": w.last(),
        "values": v,
    })
}

const IC_ARMS: [&str; 4] = ["matrix-f4", "matrix-f5", "inherited-f4", "f6-ic"];
const RHO_ARMS: [&str; 2] = ["rho-strong", "rho-negation"];

fn cell(slug: &str, runs: &[&Run]) -> Value {
    let first = runs[0];
    // arm -> workload -> run
    let mut by: BTreeMap<&str, BTreeMap<&str, &Run>> = BTreeMap::new();
    for r in runs {
        by.entry(r.arm.as_str())
            .or_default()
            .insert(r.workload.as_str(), r);
    }
    let mut arms = serde_json::Map::new();
    let mut failed = Vec::new();
    for (arm, rs) in &by {
        let n = rs.len();
        let verified = rs.values().filter(|r| r.status == "verified").count();
        let timed_out = rs.values().filter(|r| r.status == "timeout").count();
        if IC_ARMS.contains(arm) && verified + timed_out < n {
            failed.push(arm.to_string());
        }
        let ok: Vec<&&Run> = rs.values().filter(|r| r.status == "verified").collect();
        let field = |f: &dyn Fn(&Run) -> Option<f64>| -> Value {
            let v: Vec<f64> = ok.iter().filter_map(|r| f(r)).collect();
            if v.is_empty() || timed_out * 2 > n {
                Value::Null
            } else {
                summary(&v)
            }
        };
        let mut statuses: BTreeMap<String, usize> = BTreeMap::new();
        let mut isolation: BTreeMap<String, usize> = BTreeMap::new();
        for r in rs.values() {
            *statuses.entry(r.status.clone()).or_default() += 1;
            *isolation.entry(r.isolation.clone()).or_default() += 1;
        }
        arms.insert(
            arm.to_string(),
            json!({
                "runs": n,
                "statuses": statuses,
                "isolation": isolation,
                "s": field(&|r| r.s),
                "word_xors": field(&|r| r.word_xors.map(|x| x as f64)),
                "solver_calls": field(&|r| r.calls.map(|x| x as f64)),
                "reductions": field(&|r| r.reductions.map(|x| x as f64)),
                "geometric_additions": field(&|r| r.geometry.map(|x| x as f64)),
                "relations": field(&|r| r.relations.map(|x| x as f64)),
                "trials": field(&|r| r.trials.map(|x| x as f64)),
                "relations_wall_s": field(&|r| r.relations_wall_ns.map(|x| x as f64 / 1e9)),
                "base_points": ok.iter().find_map(|r| r.base_points),
            }),
        );
    }
    // per-target ratios over workloads both arms verified
    let paired = |a: &str, b: &str, f: &dyn Fn(&Run, &Run) -> Option<f64>| -> Vec<f64> {
        let (Some(ra), Some(rb)) = (by.get(a), by.get(b)) else {
            return Vec::new();
        };
        ra.iter()
            .filter_map(|(w, x)| {
                let y = rb.get(w)?;
                (x.status == "verified" && y.status == "verified")
                    .then(|| f(x, y))
                    .flatten()
            })
            .collect()
    };
    let ratio = |a: &str, b: &str| {
        paired(a, b, &|x, y| {
            Some(x.word_xors? as f64 / y.word_xors? as f64)
        })
    };
    let failed_cell = !failed.is_empty();
    let mut q = serde_json::Map::new();

    // Q1: F5 against F4
    let f5f4 = ratio("matrix-f5", "matrix-f4");
    if !f5f4.is_empty() {
        let same_tree = paired("matrix-f5", "matrix-f4", &|x, y| {
            Some(
                (x.relations == y.relations && x.trials == y.trials && x.reductions == y.reductions)
                    as u8 as f64,
            )
        });
        let all_below = f5f4.iter().all(|&v| v < 1.0);
        let all_same = same_tree.iter().all(|&v| v == 1.0);
        q.insert(
            "q1_f5_over_f4_word_xors".into(),
            json!({
                "ratio": summary(&f5f4),
                "same_tree_on_every_target": all_same,
                "verdict": if failed_cell { "failed cell" }
                    else if all_below && all_same { "F5 cuts F4's work" }
                    else { "reported ratio (rule not met)" },
            }),
        );
    }
    // Q2: matrix F4 against inherited
    let f4inh = ratio("matrix-f4", "inherited-f4");
    if !f4inh.is_empty() {
        q.insert(
            "q2_matrix_f4_over_inherited_word_xors".into(),
            json!({"ratio": summary(&f4inh)}),
        );
    }
    // Q3: F6-IC against inherited
    let l = paired("inherited-f4", "f6-ic", &|x, y| {
        Some((x.word_xors? as f64 / y.word_xors? as f64).log2())
    });
    if !l.is_empty() {
        let mut w = l.clone();
        let med = median(&mut w).unwrap();
        let all_pos = l.iter().all(|&v| v > 0.0);
        let verdict = if failed_cell {
            "failed cell"
        } else if med >= 0.1 && all_pos {
            "F6-IC saves solver work"
        } else if med.abs() <= 0.05 {
            "no saving"
        } else {
            "mixed"
        };
        let s_ratio = paired("f6-ic", "inherited-f4", &|x, y| Some(x.s? / y.s?));
        // Derived, not a preregistered rule: the arms' trees can diverge
        // (each engine returns its own decomposition, so runs stop after
        // different numbers of trials), so also report work per solver call.
        let per_call = paired("inherited-f4", "f6-ic", &|x, y| {
            let a = x.word_xors? as f64 / x.calls? as f64;
            let b = y.word_xors? as f64 / y.calls? as f64;
            Some((a / b).log2())
        });
        let same_calls = paired("inherited-f4", "f6-ic", &|x, y| {
            Some((x.calls == y.calls) as u8 as f64)
        });
        let wall = paired("inherited-f4", "f6-ic", &|x, y| {
            Some(x.relations_wall_ns? as f64 / y.relations_wall_ns? as f64)
        });
        q.insert(
            "q3_log2_inherited_over_f6_word_xors".into(),
            json!({
                "log2_ratio": summary(&l),
                "verdict": verdict,
                "derived_log2_per_call": summary(&per_call),
                "targets_with_equal_solver_calls": same_calls.iter().filter(|&&v| v == 1.0).count(),
                "s_f6_over_inherited": summary(&s_ratio),
                "relations_wall_inherited_over_f6": summary(&wall),
            }),
        );
    }
    // Q4: every IC arm against each rho arm
    let mut q4 = serde_json::Map::new();
    for ic in IC_ARMS {
        for rho in RHO_ARMS {
            let v = paired(ic, rho, &|x, y| Some(x.s? / y.s?));
            if v.is_empty() {
                continue;
            }
            let below = v.iter().all(|&x| x < 1.0);
            q4.insert(
                format!("{ic}/{rho}"),
                json!({
                    "s_ratio": summary(&v),
                    "verdict": if below { "IC lower bound below rho on every target (proves nothing; word XORs unpriced)" }
                        else if v.iter().all(|&x| x > 1.0) { "IC slower than rho on every target (lower bound above rho)" }
                        else { "mixed" },
                }),
            );
        }
    }
    if !q4.is_empty() {
        q.insert("q4_s_over_rho".into(), Value::Object(q4));
    }
    let floor = if first.family == "koblitz" {
        (std::f64::consts::PI / (4.0 * first.m as f64)).sqrt()
    } else {
        (std::f64::consts::PI / 4.0).sqrt()
    };
    json!({
        "slug": slug,
        "family": first.family,
        "m": first.m,
        "r": first.r,
        "floor_s": floor,
        "failed_ic_arms": failed,
        "arms": arms,
        "questions": q,
    })
}

fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        return Err("usage: f4f5f6_cross_family_analyze OUT.json SESSION_DIR...".into());
    }
    let mut sessions = serde_json::Map::new();
    for dir in &args[2..] {
        let (id, runs) = read(dir)?;
        let mut cells: BTreeMap<String, Vec<&Run>> = BTreeMap::new();
        for r in &runs {
            cells.entry(r.slug.clone()).or_default().push(r);
        }
        let cells: Vec<Value> = cells.iter().map(|(slug, rs)| cell(slug, rs)).collect();
        let name = std::path::Path::new(dir)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or(dir)
            .to_string();
        sessions.insert(name, json!({"session_id": id, "dir": dir, "cells": cells}));
    }
    let out = json!({
        "schema": "f4f5f6-cross-family-analysis/v1",
        "protocol": "research/f4_f5_f6_cross_family_20261011/PROTOCOL.md",
        "unit": {"s": "ecbench.gae / sqrt(r), lower bound (solver word XORs unpriced; F6-IC geometric additions charged)",
                 "word_xors": "solver.ops, counted, deterministic"},
        "sessions": sessions,
    });
    let text = serde_json::to_string_pretty(&out).map_err(|e| e.to_string())?;
    std::fs::write(&args[1], text + "\n").map_err(|e| e.to_string())?;
    println!("wrote {}", args[1]);
    Ok(())
}
