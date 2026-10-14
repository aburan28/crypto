//! The Shantz–Teske delta method, measured net: for fixed `(n, m)` on `K_0`
//! sweep the subspace dimension `ℓ` over `m·ℓ = n − δ` and report, per
//! `δ`, the solving degree, F4 wall and word XORs per system, the fraction
//! of targets that decompose (measured against the `2^{−δ}` prediction),
//! and the net cost per found relation.
//!
//! ```bash
//! CARGO_TARGET_DIR=target cargo run --release --example delta_method_sweep -- \
//!     --out research/ic_techniques_20261012/a2_cyclotomic_and_delta \
//!     --targets 32 --cell-timeout-secs 120 --budget-secs 1200
//! ```
//!
//! Cells run in child processes (`--cell JSON`) with a wall budget, one
//! JSON line per target, so a cell too slow for its budget is truncated
//! and reported truncated rather than stalling the sweep.  Deterministic
//! seeds; public synthetic targets `[k]G`; no logarithm computed.

use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use crypto_lib::cryptanalysis::delta_method::{
    delta_dimensions, run_delta_cell, seeded_targets, standard_base, summarise_cell,
    DeltaCellResult, DeltaChoice, DeltaTargetRecord,
};
use crypto_lib::cryptanalysis::koblitz_groebner::SolverEngine;
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

struct Args {
    out: PathBuf,
    targets: usize,
    seed: u64,
    cell_timeout: Duration,
    budget: Duration,
    node_budget: usize,
    delta_min: i32,
    delta_max: i32,
    degrees: Vec<u32>,
    summands: Vec<usize>,
    a: u8,
    engine: String,
    cell: Option<String>,
    /// Suffix of the output files (`delta_results_<label>.json`), so
    /// several bounded runs can share one directory.
    label: Option<String>,
}

fn parse_args() -> Args {
    let mut args = Args {
        out: PathBuf::from("research/ic_techniques_20261012/a2_cyclotomic_and_delta"),
        targets: 32,
        seed: 20261012,
        cell_timeout: Duration::from_secs(120),
        budget: Duration::from_secs(1200),
        node_budget: 4096,
        delta_min: -2,
        delta_max: 4,
        degrees: vec![17, 23, 31],
        summands: vec![2, 3],
        a: 0,
        engine: "inherited".into(),
        cell: None,
        label: None,
    };
    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut value = || it.next().expect("flag value");
        match flag.as_str() {
            "--out" => args.out = PathBuf::from(value()),
            "--targets" => args.targets = value().parse().expect("targets"),
            "--seed" => args.seed = value().parse().expect("seed"),
            "--cell-timeout-secs" => {
                args.cell_timeout = Duration::from_secs(value().parse().expect("secs"))
            }
            "--budget-secs" => args.budget = Duration::from_secs(value().parse().expect("secs")),
            "--node-budget" => args.node_budget = value().parse().expect("budget"),
            "--delta-min" => args.delta_min = value().parse().expect("delta"),
            "--delta-max" => args.delta_max = value().parse().expect("delta"),
            "--degrees" => {
                args.degrees = value()
                    .split(',')
                    .map(|s| s.parse().expect("degree"))
                    .collect()
            }
            "--summands" => {
                args.summands = value().split(',').map(|s| s.parse().expect("m")).collect()
            }
            "--a" => args.a = value().parse().expect("a"),
            "--engine" => args.engine = value(),
            "--cell" => args.cell = Some(value()),
            "--label" => args.label = Some(value()),
            other => panic!("unknown flag {other}"),
        }
    }
    args
}

fn engine_of(name: &str) -> SolverEngine {
    match name {
        "matrix" | "matrix_f4" => SolverEngine::MatrixF4 { max_degree: 3 },
        "inherited" | "inherited_f4" => SolverEngine::InheritedF4 { max_degree: 3 },
        "f5" | "matrix_f5" => SolverEngine::MatrixF5 { max_degree: 3 },
        other => panic!("unknown engine {other}"),
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct CellSpec {
    a: u8,
    n: u32,
    m: usize,
    ell: u32,
    targets: usize,
    seed: u64,
    node_budget: usize,
    engine: String,
}

fn run_cell_child(spec: &CellSpec) {
    let kc = KoblitzCurve::new(spec.a, spec.n).expect("usable Koblitz curve");
    let targets = seeded_targets(&kc, spec.targets, spec.seed);
    let stdout = std::io::stdout();
    let result = run_delta_cell(
        &kc,
        spec.m,
        spec.ell,
        &targets,
        engine_of(&spec.engine),
        spec.node_budget,
        None,
        |record| {
            let mut lock = stdout.lock();
            writeln!(lock, "{}", serde_json::to_string(record).unwrap()).unwrap();
            lock.flush().unwrap();
        },
    )
    .expect("cell runs");
    let mut v = serde_json::to_value(&result).unwrap();
    v.as_object_mut().unwrap().remove("records");
    println!("{}", json!({"summary": v}));
}

/// Parent: run the child under a wall budget, collecting records.
fn run_cell_parent(spec: &CellSpec, timeout: Duration) -> DeltaCellResult {
    let started = Instant::now();
    let exe = std::env::current_exe().expect("current exe");
    let mut child = Command::new(exe)
        .arg("--cell")
        .arg(serde_json::to_string(spec).unwrap())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .expect("spawn cell child");
    let stdout = child.stdout.take().expect("child stdout");
    let (tx, rx) = std::sync::mpsc::channel::<String>();
    std::thread::spawn(move || {
        for line in BufReader::new(stdout).lines().map_while(Result::ok) {
            if tx.send(line).is_err() {
                break;
            }
        }
    });
    let deadline = started + timeout;
    let mut records: Vec<DeltaTargetRecord> = Vec::new();
    let mut truncated = false;
    let take = |line: &str, records: &mut Vec<DeltaTargetRecord>| {
        if let Ok(r) = serde_json::from_str::<DeltaTargetRecord>(line) {
            eprintln!(
                "    n{} m{} ℓ{} target k={} {} {:.0} ms deg {} splits {}",
                spec.n,
                spec.m,
                spec.ell,
                r.scalar,
                if r.decomposed {
                    "DECOMPOSED"
                } else {
                    "refuted"
                },
                r.wall_ms,
                r.max_degree_built,
                r.splits
            );
            records.push(r);
        }
    };
    loop {
        match rx.recv_timeout(Duration::from_millis(200)) {
            Ok(line) => take(&line, &mut records),
            Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {}
            Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => break,
        }
        if let Ok(Some(_)) = child.try_wait() {
            while let Ok(line) = rx.try_recv() {
                take(&line, &mut records);
            }
            break;
        }
        if Instant::now() >= deadline {
            let _ = child.kill();
            let _ = child.wait();
            truncated = true;
            break;
        }
    }
    let _ = child.wait();
    let kc = KoblitzCurve::new(spec.a, spec.n).expect("usable Koblitz curve");
    let (base, _) = standard_base(&kc, spec.ell).expect("base");
    summarise_cell(
        &kc,
        spec.m,
        spec.ell,
        &base,
        spec.targets,
        records,
        truncated,
    )
}

fn main() {
    let args = parse_args();
    if let Some(cell) = &args.cell {
        let spec: CellSpec = serde_json::from_str(cell).expect("cell spec");
        run_cell_child(&spec);
        return;
    }
    let sweep_started = Instant::now();
    let deadline = sweep_started + args.budget;
    std::fs::create_dir_all(&args.out).expect("output dir");

    // Plan: every admissible ℓ per (n, m), cheapest (smallest system) first
    // within each (n, m) so a budget stop loses the largest cells.
    let mut plan: Vec<(CellSpec, DeltaChoice)> = Vec::new();
    for &n in &args.degrees {
        if KoblitzCurve::new(args.a, n).is_none() {
            eprintln!("K_{}/F_2^{n}: no usable subgroup, skipped", args.a);
            continue;
        }
        for &m in &args.summands {
            for choice in delta_dimensions(n, m, args.delta_min, args.delta_max) {
                plan.push((
                    CellSpec {
                        a: args.a,
                        n,
                        m,
                        ell: choice.ell,
                        targets: args.targets,
                        seed: args.seed,
                        node_budget: args.node_budget,
                        engine: args.engine.clone(),
                    },
                    choice,
                ));
            }
        }
    }
    eprintln!(
        "{} cells, {} targets each, {} s per cell, {} s total budget, engine {}",
        plan.len(),
        args.targets,
        args.cell_timeout.as_secs(),
        args.budget.as_secs(),
        args.engine
    );
    let mut results: Vec<DeltaCellResult> = Vec::new();
    let mut skipped: Vec<Value> = Vec::new();
    for (spec, choice) in &plan {
        let remaining = deadline.saturating_duration_since(Instant::now());
        if remaining < Duration::from_secs(5) {
            skipped.push(json!({"n": spec.n, "m": spec.m, "ell": spec.ell, "delta": choice.delta, "reason": "sweep budget exhausted"}));
            continue;
        }
        eprintln!(
            "-- n={} m={} ℓ={} δ={} vars={}",
            spec.n, spec.m, spec.ell, choice.delta, choice.n_vars
        );
        let result = run_cell_parent(spec, args.cell_timeout.min(remaining));
        eprintln!(
            "   measured {}/{} decomposed {} ({:.3} vs predicted {:.3}, 2^-δ = {}), mean {:.1} ms, deg {}, split_free {}, net {:.1} ms/relation{}",
            result.targets_measured,
            result.targets_requested,
            result.decomposed,
            result.measured_fraction,
            result.predicted_fraction,
            result.predicted_delta_factor,
            result.mean_wall_ms,
            result.max_degree_built,
            result.split_free,
            result.net_ms_per_relation,
            if result.truncated { " [TRUNCATED]" } else { "" }
        );
        results.push(result);
    }

    let report = json!({
        "lead": "Shantz–Teske delta method, net cost per relation",
        "generated_utc": format!("unix:{}", std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0)),
        "args": {
            "a": args.a, "degrees": args.degrees, "summands": args.summands,
            "delta_min": args.delta_min, "delta_max": args.delta_max,
            "targets": args.targets, "seed": args.seed,
            "cell_timeout_secs": args.cell_timeout.as_secs(), "budget_secs": args.budget.as_secs(),
            "node_budget": args.node_budget, "engine": args.engine,
            "factor_base": "standard polynomial-basis subspace <1, z, …, z^(ℓ−1)>, every rational point above it",
            "oracle": "koblitz_index_calculus::groebner_decompose (chained S3 Weil descent, m·ℓ + (m−2)·n unknowns)",
        },
        "host": {"os": std::env::consts::OS, "arch": std::env::consts::ARCH},
        "cells": results,
        "skipped": skipped,
        "wall_ms": sweep_started.elapsed().as_secs_f64() * 1e3,
    });
    let suffix = args
        .label
        .as_ref()
        .map(|l| format!("_{l}"))
        .unwrap_or_default();
    let path = args.out.join(format!("delta_results{suffix}.json"));
    std::fs::write(&path, serde_json::to_string_pretty(&report).unwrap()).expect("write results");
    write_table(
        &args.out.join(format!("delta_summary{suffix}.md")),
        &results,
        &skipped,
    );
    eprintln!(
        "wrote {} after {:.1} s",
        path.display(),
        sweep_started.elapsed().as_secs_f64()
    );
}

fn write_table(path: &Path, results: &[DeltaCellResult], skipped: &[Value]) {
    let mut s = String::from("# Delta method summary (generated)\n\n| n | m | δ | ℓ | vars | |F| | targets | decomposed | measured frac | predicted frac | 2^-δ | F4 max deg | split-free | mean splits | mean ms/system | median ms | mean word XORs | net ms/relation | net XORs/relation | truncated |\n|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|:--|\n");
    for r in results {
        s.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {}/{} | {} | {:.3} | {:.3} | {} | {} | {} | {:.2} | {:.1} | {:.1} | {:.3e} | {:.1} | {:.3e} | {} |\n",
            r.n, r.m, r.delta, r.ell, r.n_vars, r.base_points, r.targets_measured, r.targets_requested,
            r.decomposed, r.measured_fraction, r.predicted_fraction, r.predicted_delta_factor,
            r.max_degree_built, r.split_free, r.mean_splits, r.mean_wall_ms, r.median_wall_ms, r.mean_word_xors,
            r.net_ms_per_relation, r.net_word_xors_per_relation, r.truncated
        ));
    }
    if !skipped.is_empty() {
        s.push_str("\nSkipped cells:\n");
        for v in skipped {
            s.push_str(&format!("- {v}\n"));
        }
    }
    std::fs::write(path, s).expect("write summary");
}
