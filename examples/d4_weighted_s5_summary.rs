//! Summarise `research/d4_symmetrised_s5_20261015/raw/*.jsonl` (written by
//! `d4_weighted_s5`) into the study's tables and `summary.json`: per size the
//! median and range of each arm's cost, the paired ratio, the ratio to the
//! quoted plain `C₄`, the lift-back share, exactness counts, and the
//! handover in E18's `τ_T` frame.  Arithmetic on recorded values only.
//!
//! ```text
//! d4_weighted_s5_summary RAW_DIR > summary.txt   # also writes RAW_DIR/../summary.json
//! ```

use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::Path;

use serde_json::{json, Value};

/// Quoted plain `C₄` (§11.17, `experiments/24_gaudry_quartic_c4.json`).
const C4_PLAIN: f64 = 1.213e12;
/// E18 `τ_T`-folded handover at `C₄_PLAIN`, bits of `n` (§9.3).
const NSTAR_TAU_PLAIN: f64 = 157.9;

fn read(path: &Path) -> Vec<Value> {
    fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).expect("jsonl"))
        .collect()
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn range(v: &[f64]) -> (f64, f64) {
    v.iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &x| {
            (a.min(x), b.max(x))
        })
}

fn nstar(c: f64) -> f64 {
    NSTAR_TAU_PLAIN - 4.0 * (C4_PLAIN / c).log2()
}

fn f(v: &Value, k: &str) -> f64 {
    v[k].as_f64().unwrap_or(f64::NAN)
}

fn main() {
    let dir = env::args().nth(1).expect("RAW_DIR");
    let dir = Path::new(&dir);
    // `--prefix v2_` reads `v2_tiny.jsonl` etc. (the re-run after the
    // eigenspace-resolution fix); the default reads the first run.
    let args: Vec<String> = env::args().collect();
    let prefix = args
        .iter()
        .position(|a| a == "--prefix")
        .map_or(String::new(), |i| args[i + 1].clone());
    let out_name = format!("{prefix}summary.json");
    let mut rows = read(&dir.join(format!("{prefix}tiny.jsonl")));
    rows.extend(read(&dir.join(format!("{prefix}measured.jsonl"))));
    rows.extend(read(&dir.join(format!("{prefix}w_only.jsonl"))));
    let controls = read(&dir.join("controls.jsonl"));

    // Exactness and structure over every target.
    let mut by_p: BTreeMap<(u64, u64), Vec<&Value>> = BTreeMap::new();
    for r in &rows {
        by_p.entry((r["p"].as_u64().unwrap(), r["seed"].as_u64().unwrap()))
            .or_default()
            .push(r);
    }
    let mut cells = Vec::new();
    println!(
        "| p | seed | targets | w solved | dim 512 | deg-27 reruns | C_w median (min–max) | L_w median (max) | L_w / C_w max | C_c8 median (min–max) | R_c8 median (min–max) | R_plain median (min–max) | rational sols (w) | split / rooted / candidates | decompositions (w) | w = brute | w = mitm | c8 = w | planted found w / c8 |"
    );
    println!("|--:|--:|--:|--:|--:|--:|:--|:--|--:|:--|:--|:--|--:|:--|--:|:--|:--|:--|:--|");
    for ((p, seed), rs) in &by_p {
        let solved: Vec<&&Value> = rs.iter().filter(|r| r["w_solved"] == true).collect();
        let dim512 = solved
            .iter()
            .filter(|r| r["w_solve"]["dim"].as_u64() == Some(512))
            .count();
        let reruns = solved
            .iter()
            .filter(|r| r["w_degree_used"].as_u64().is_some_and(|d| d > 26))
            .count();
        let mut cw: Vec<f64> = solved
            .iter()
            .map(|r| f(&r["w_solve"], "total_muls"))
            .collect();
        let mut lw: Vec<f64> = solved.iter().map(|r| f(&r["w_lift"], "muls")).collect();
        let share = solved
            .iter()
            .map(|r| f(&r["w_lift"], "muls") / f(&r["w_solve"], "total_muls"))
            .fold(0.0f64, f64::max);
        let c8: Vec<&&Value> = rs.iter().filter(|r| r["c8_solved"] == true).collect();
        let mut cc8: Vec<f64> = c8.iter().map(|r| f(&r["c8_solve"], "total_muls")).collect();
        let mut rc8: Vec<f64> = c8
            .iter()
            .filter(|r| r["w_solved"] == true)
            .map(|r| f(&r["c8_solve"], "total_muls") / f(&r["w_solve"], "total_muls"))
            .collect();
        let mut rpl: Vec<f64> = cw.iter().map(|c| C4_PLAIN / c).collect();
        let sum = |k: &str| -> u64 { rs.iter().map(|r| r[k].as_u64().unwrap_or(0)).sum() };
        let lsum = |k: &str| -> u64 {
            rs.iter()
                .map(|r| r["w_lift"][k].as_u64().unwrap_or(0))
                .sum()
        };
        let count = |k: &str| -> (usize, usize) {
            let have: Vec<bool> = rs.iter().filter_map(|r| r[k].as_bool()).collect();
            (have.iter().filter(|&&b| b).count(), have.len())
        };
        let (cwr, lwr, c8r, rc8r, rplr) = (
            range(&cw),
            range(&lw),
            range(&cc8),
            range(&rc8),
            range(&rpl),
        );
        let (cwm, lwm, c8m, rc8m, rplm) = (
            median(&mut cw),
            median(&mut lw),
            median(&mut cc8),
            median(&mut rc8),
            median(&mut rpl),
        );
        let (wb, wbn) = count("w_equals_brute");
        let (wm, wmn) = count("w_equals_mitm");
        let (cw8, cw8n) = count("c8_equals_w");
        let (pw, pwn) = count("planted_found_w");
        let (pc, pcn) = count("planted_found_c8");
        println!(
            "| {p} | {seed} | {} | {} | {dim512} | {reruns} | {cwm:.3e} ({:.3e}–{:.3e}) | {lwm:.2e} ({:.2e}) | {share:.2e} | {} | {} | {rplm:.0} ({:.0}–{:.0}) | {} | {} / {} / {} | {} | {wb}/{wbn} | {wm}/{wmn} | {cw8}/{cw8n} | {pw}/{pwn} / {pc}/{pcn} |",
            rs.len(),
            solved.len(),
            cwr.0,
            cwr.1,
            lwr.1,
            if cc8.is_empty() { "—".into() } else { format!("{c8m:.3e} ({:.3e}–{:.3e})", c8r.0, c8r.1) },
            if rc8.is_empty() { "—".into() } else { format!("{rc8m:.1} ({:.1}–{:.1})", rc8r.0, rc8r.1) },
            rplr.0,
            rplr.1,
            sum("w_rational_solutions"),
            lsum("split"),
            lsum("rooted"),
            lsum("candidates"),
            sum("w_decompositions"),
        );
        cells.push(json!({
            "p": p, "seed": seed, "targets": rs.len(), "w_solved": solved.len(),
            "dim512": dim512, "reruns_deg27": reruns,
            "c_w_median": cwm, "c_w_min": cwr.0, "c_w_max": cwr.1,
            "l_w_median": lwm, "l_w_max": lwr.1, "l_over_c_max": share,
            "c_c8_median": if cc8.is_empty() { Value::Null } else { json!(c8m) },
            "c_c8_min": if cc8.is_empty() { Value::Null } else { json!(c8r.0) },
            "c_c8_max": if cc8.is_empty() { Value::Null } else { json!(c8r.1) },
            "r_c8_median": if rc8.is_empty() { Value::Null } else { json!(rc8m) },
            "r_c8_min": if rc8.is_empty() { Value::Null } else { json!(rc8r.0) },
            "r_c8_max": if rc8.is_empty() { Value::Null } else { json!(rc8r.1) },
            "r_plain_median": rplm, "r_plain_min": rplr.0, "r_plain_max": rplr.1,
            "w_equals_brute": [wb, wbn], "w_equals_mitm": [wm, wmn], "c8_equals_w": [cw8, cw8n],
            "planted_found_w": [pw, pwn], "planted_found_c8": [pc, pcn],
            "rational_solutions": sum("w_rational_solutions"),
            "split": lsum("split"), "rooted": lsum("rooted"), "candidates": lsum("candidates"),
            "decompositions": sum("w_decompositions"),
        }));
    }

    // The handover, E18's τ_T frame (derived), on one paired target set:
    // the measured sizes, seed 1, targets solved by both arms.
    let measured = [269u64, 1033, 4093];
    let paired: Vec<&Value> = rows
        .iter()
        .filter(|r| measured.contains(&r["p"].as_u64().unwrap()))
        .filter(|r| r["seed"].as_u64() == Some(1))
        .filter(|r| r["w_solved"] == true && r["c8_solved"] == true)
        .collect();
    let mut cw_all: Vec<f64> = paired
        .iter()
        .map(|r| f(&r["w_solve"], "total_muls"))
        .collect();
    let mut cwl_all: Vec<f64> = paired
        .iter()
        .map(|r| f(&r["w_solve"], "total_muls") + f(&r["w_lift"], "muls"))
        .collect();
    let mut c8_all: Vec<f64> = paired
        .iter()
        .map(|r| f(&r["c8_solve"], "total_muls"))
        .collect();
    let n_paired = paired.len();
    let (cwm, cwlm, c8m) = (
        median(&mut cw_all),
        median(&mut cwl_all),
        median(&mut c8_all),
    );
    println!();
    println!("Handover in E18's τ_T frame (derived; extrapolation): n* = {NSTAR_TAU_PLAIN} − 4·log2({C4_PLAIN:.3e} / C)");
    println!(
        "Paired set: p ∈ {{269, 1033, 4093}}, seed 1, solved by both arms: {n_paired} targets."
    );
    println!("| solve | C (median, paired set) | C4_plain / C | n* (bits) |");
    println!("|:--|--:|--:|--:|");
    let hand = [
        ("plain S₄ form (quoted)", C4_PLAIN),
        ("C8, five unknowns", c8m),
        ("w, four weighted unknowns", cwm),
        ("w + lift-back", cwlm),
    ];
    let mut hv = Vec::new();
    for (name, c) in hand {
        println!(
            "| {name} | {c:.3e} | {:.1} | {:.1} |",
            C4_PLAIN / c,
            nstar(c)
        );
        hv.push(json!({"solve": name, "c": c, "ratio": C4_PLAIN / c, "nstar_bits": nstar(c)}));
    }
    println!();
    println!("Controls:");
    for c in &controls {
        println!("{c}");
    }
    let out = json!({
        "c4_plain_quoted": C4_PLAIN,
        "nstar_tau_plain_bits": NSTAR_TAU_PLAIN,
        "cells": cells,
        "handover_tau_frame": hv,
        "handover_paired_targets": n_paired,
        "controls": controls,
    });
    let path = dir.parent().unwrap().join(out_name);
    fs::write(&path, serde_json::to_string_pretty(&out).unwrap()).expect("write summary");
    eprintln!("wrote {}", path.display());
}
