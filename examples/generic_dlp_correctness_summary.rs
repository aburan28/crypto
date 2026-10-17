//! Summary for `research/generic_dlp_correctness_20261017` (PROTOCOL.md):
//! reads the study's `raw/` counts and `transcripts/`, writes
//! `summary.json` beside them, and prints a Markdown summary.  Every number
//! it reports is a count; nothing is timed.
//!
//! ```sh
//! cargo run --release --example generic_dlp_correctness_summary -- research/generic_dlp_correctness_20261017
//! ```

use std::collections::BTreeMap;
use std::path::Path;

use serde_json::{json, Value};

fn jsonl(path: &Path) -> Vec<Value> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).expect("json line"))
        .collect()
}

fn text(path: &Path) -> String {
    std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The sections of a collab transcript, keyed by their header
/// (`<curve> neg=<bool> dp=<bits>`), each with its full text and its counts
/// of `dp` and `dead` walker lines.
fn collab_sections(t: &str) -> BTreeMap<String, (String, u64, u64)> {
    let mut out: BTreeMap<String, (String, u64, u64)> = BTreeMap::new();
    let mut key = String::new();
    for line in t.lines() {
        if line.contains(" neg=") && !line.starts_with(' ') {
            key = line.split(':').next().unwrap().to_string();
            out.entry(key.clone()).or_default();
        }
        let e = out.entry(key.clone()).or_default();
        e.0.push_str(line);
        e.0.push('\n');
        let mut words = line.split_whitespace();
        let _walker = words.next();
        match words.next() {
            Some("dp") => e.1 += 1,
            Some("dead") => e.2 += 1,
            _ => {}
        }
    }
    out
}

/// A Pohlig–Hellman transcript line split into (everything but the trailing
/// step count, the step count); Cheon lines and headers have no count.
fn ph_line(line: &str) -> (String, Option<u64>) {
    let is_ph = !line.starts_with("cheon")
        && (line.contains(" d=") || line.trim_start().starts_with("wrong-order:"));
    if !is_ph {
        return (line.to_string(), None);
    }
    let (head, last) = line.rsplit_once(' ').expect("step count");
    (
        head.to_string(),
        Some(last.parse().expect("integer step count")),
    )
}

fn main() {
    let dir = std::path::PathBuf::from(std::env::args().nth(1).expect("study directory"));
    let before = jsonl(&dir.join("raw/counts_before.jsonl"));
    let after = jsonl(&dir.join("raw/counts_after.jsonl"));

    // ── collab walkers ──
    let key = |v: &Value| {
        (
            v["curve"].as_str().unwrap().to_string(),
            v["negation_map"].as_bool().unwrap(),
        )
    };
    let mut collab = Vec::new();
    println!("| curve | negation map | walkers | DPs before → after | dead trails before → after | DPs verified before → after |");
    println!("|:--|:--|--:|:--|:--|:--|");
    for b in before.iter().filter(|v| v["kind"] == "collab") {
        let a = after
            .iter()
            .find(|v| v["kind"] == "collab" && key(v) == key(b))
            .expect("matching after row");
        println!(
            "| {} | {} | {} | {} → {} | {} → {} | {} → {} |",
            b["curve"].as_str().unwrap(),
            b["negation_map"],
            b["walkers"],
            b["dps"],
            a["dps"],
            b["dead_trails"],
            a["dead_trails"],
            b["dps_verified"],
            a["dps_verified"]
        );
        collab.push(json!({
            "curve": b["curve"], "n_bits": b["n_bits"], "negation_map": b["negation_map"],
            "dp_bits": b["dp_bits"], "walkers": b["walkers"],
            "before": {"dps": b["dps"], "dead_trails": b["dead_trails"], "dps_verified": b["dps_verified"], "dp_steps": b["dp_steps"], "dead_steps": b["dead_steps"]},
            "after": {"dps": a["dps"], "dead_trails": a["dead_trails"], "dps_verified": a["dps_verified"], "dp_steps": a["dp_steps"], "dead_steps": a["dead_steps"]},
            "all_dps_verified": b["dps"] == b["dps_verified"] && a["dps"] == a["dps_verified"],
            "unchanged": b["dps"] == a["dps"] && b["dead_trails"] == a["dead_trails"] && b["dp_steps"] == a["dp_steps"] && b["dead_steps"] == a["dead_steps"],
        }));
    }

    // ── Pohlig–Hellman, driver panel ──
    let mut ph = Vec::new();
    let (mut sb, mut sa, mut correct, mut same) = (0u64, 0u64, true, true);
    for b in before.iter().filter(|v| v["kind"] == "ph") {
        let a = after
            .iter()
            .find(|v| v["kind"] == "ph" && v["target"] == b["target"] && v["bound"] == b["bound"])
            .expect("matching after row");
        let (x, y) = (
            b["total_steps"].as_u64().unwrap(),
            a["total_steps"].as_u64().unwrap(),
        );
        sb += x;
        sa += y;
        correct &= b["correct"] == json!(true) && a["correct"] == json!(true);
        same &= b["recovered"] == a["recovered"];
        ph.push(json!({"target": b["target"], "bound": b["bound"], "steps_before": x, "steps_after": y, "ratio": x as f64 / y as f64}));
    }
    let ratios: Vec<f64> = ph.iter().map(|r| r["ratio"].as_f64().unwrap()).collect();
    let (rmin, rmax) = (
        ratios.iter().cloned().fold(f64::INFINITY, f64::min),
        ratios.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    println!();
    println!(
        "Pohlig–Hellman driver panel: {} solves, all correct before and after: {correct}, recovered logs identical: {same}; total_steps {sb} → {sa} ({:.2}×), per solve {rmin:.2}–{rmax:.2}×",
        ph.len(),
        sb as f64 / sa as f64
    );

    // ── pinned transcripts ──
    let (cb, ca) = (
        collab_sections(&text(&dir.join("transcripts/collab_before.txt"))),
        collab_sections(&text(&dir.join("transcripts/collab_after.txt"))),
    );
    let mut sections = Vec::new();
    let mut noneg_identical = true;
    for (k, (tb, pb, db)) in &cb {
        let (ta, pa, da) = ca.get(k).expect("same sections after");
        if k.contains("neg=false") {
            noneg_identical &= tb == ta;
        }
        sections.push(json!({"section": k, "dp_before": pb, "dead_before": db, "dp_after": pa, "dead_after": da, "identical": tb == ta}));
    }
    let neg_dead = |at: &BTreeMap<String, (String, u64, u64)>| -> u64 {
        at.iter()
            .filter(|(k, _)| k.contains("neg=true"))
            .map(|(_, v)| v.2)
            .sum()
    };
    let (pb, pa) = (
        text(&dir.join("transcripts/ph_cheon_before.txt")),
        text(&dir.join("transcripts/ph_cheon_after.txt")),
    );
    let (lb, la): (Vec<_>, Vec<_>) = (
        pb.lines().map(ph_line).collect(),
        pa.lines().map(ph_line).collect(),
    );
    let ph_same_but_steps = lb.len() == la.len() && lb.iter().zip(&la).all(|(x, y)| x.0 == y.0);
    let steps = |v: &[(String, Option<u64>)]| -> (u64, usize) {
        (
            v.iter().filter_map(|x| x.1).sum(),
            v.iter().filter(|x| x.1.is_some()).count(),
        )
    };
    let ((tb, nb), (ta, _)) = (steps(&lb), steps(&la));
    println!(
        "pinned collab transcripts: neg=false sections byte-identical: {noneg_identical}; neg=true dead trails {} → {}",
        neg_dead(&cb),
        neg_dead(&ca)
    );
    println!(
        "pinned Pohlig–Hellman/Cheon transcript: identical apart from step counts: {ph_same_but_steps}; {nb} Pohlig–Hellman lines, total_steps {tb} → {ta}"
    );

    let summary = json!({
        "schema": "generic-dlp-correctness-summary/v1",
        "collab": collab,
        "ph_panel": {"solves": ph.len(), "all_correct": correct, "recovered_identical": same,
                     "total_steps_before": sb, "total_steps_after": sa, "ratio": sb as f64 / sa as f64,
                     "per_solve_ratio_min": rmin, "per_solve_ratio_max": rmax, "rows": ph},
        "pinned": {
            "collab_sections": sections,
            "collab_neg_false_identical": noneg_identical,
            "collab_neg_true_dead_before": neg_dead(&cb),
            "collab_neg_true_dead_after": neg_dead(&ca),
            "ph_identical_apart_from_steps": ph_same_but_steps,
            "ph_lines": nb,
            "ph_total_steps_before": tb,
            "ph_total_steps_after": ta,
        },
    });
    std::fs::write(
        dir.join("summary.json"),
        serde_json::to_string_pretty(&summary).unwrap() + "\n",
    )
    .expect("write summary.json");
    eprintln!("wrote {}", dir.join("summary.json").display());
}
