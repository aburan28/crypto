//! Matched pair-table IC vs signed-Frobenius rho in exclusive group ops.
//!
//! Frozen protocol: `RESEARCH_ECC2K130_RHO_PARITY.md`. Do not change the
//! cells, seeds, row formula, window rule, or conversion after seeing a
//! cell. `--quick` is a smoke path, not a result.
//!
//! ```bash
//! cargo run --release --example rho_parity_e2e
//! cargo run --release --example rho_parity_e2e -- --quick
//! ```

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    binary_method_group_ops, build_subgroup_orbit_factor_base, expected_binary_method_group_ops,
    koblitz_signed_frobenius_rho_with_progress, point_key, points_with_x, probe_scalar,
    rho_exclusive_group_ops, rho_expected_steps, DecompositionStrategy, FactorBaseLogSolver,
    IndividualLogSolver, KoblitzCurve, KoblitzIcOptions, KoblitzSignedRhoOptions, PairSumTable,
    RelationCollector, RelationWorkUnit, PRECOMPUTE_BATCH_TRIALS,
};
use num_bigint::BigUint;
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde_json::json;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

const ALGORITHM_SEED: u64 = 0x5EED_0007;
const FB_SEED: u64 = 43;
const TARGET_SEEDS: [u64; 3] = [0x51EE_D001, 0x51EE_D002, 0x51EE_D003];
const REPS: u32 = 3;
const MAX_TRIALS: usize = 4_096;
const WALKS: u64 = 64;
const PROBE_RUN: u64 = 64;

struct Cell {
    a: u8,
    n: u32,
}

const CELLS: [Cell; 5] = [
    Cell { a: 0, n: 13 },
    Cell { a: 1, n: 17 },
    Cell { a: 0, n: 19 },
    Cell { a: 1, n: 19 },
    Cell { a: 0, n: 23 },
];

fn cell_name(cell: &Cell) -> String {
    format!("n{}a{}", cell.n, cell.a)
}

fn recipe_floor(n: u32) -> usize {
    6 * n as usize
}

fn hash_target(curve: &KoblitzCurve, seed: u64) -> BinaryPoint {
    for counter in 0u64..1_000_000 {
        let mut h = blake3::Hasher::new();
        h.update(b"ic-workflow-public-target-v1\0");
        h.update(&curve.n.to_le_bytes());
        h.update(&[curve.a]);
        h.update(&curve.k.to_le_bytes());
        h.update(&curve.b_index.to_le_bytes());
        h.update(&seed.to_le_bytes());
        h.update(&counter.to_le_bytes());
        let digest = h.finalize();
        let x = u64::from_le_bytes(digest.as_bytes()[..8].try_into().unwrap())
            & ((1u64 << curve.n) - 1);
        let x = F2mElement::from_biguint(&BigUint::from(x), curve.n);
        let mut lifts = points_with_x(&curve.curve, &x);
        lifts.sort_by_key(point_key);
        if lifts.is_empty() {
            continue;
        }
        let q = curve.mul(
            &lifts[usize::from(digest.as_bytes()[8] & 1) % lifts.len()],
            &curve.cofactor,
        );
        if q != BinaryPoint::Infinity {
            return q;
        }
    }
    panic!("hash-to-curve exhausted");
}

fn log_ops(log: &BigUint) -> u64 {
    match log.to_u64_digits().first().copied() {
        Some(s) => binary_method_group_ops(s),
        None => 0,
    }
}

fn walked_collection_ops(seed: u64, trials: u64, r: u64, batches: u64) -> u64 {
    if trials == 0 {
        return 0;
    }
    let stride = StdRng::seed_from_u64(seed ^ 0x5354_5249_4445_5f30).gen_range(1..r.max(2));
    let mut ops = batches.saturating_mul(binary_method_group_ops(stride));
    let last = trials - 1;
    let first_run = 0u64;
    let last_run = last / PROBE_RUN;
    for run in first_run..=last_run {
        let run_start = run * PROBE_RUN;
        let run_end = ((run + 1) * PROBE_RUN).min(trials);
        let anchor = probe_scalar(seed, run_start / PROBE_RUN, r);
        ops = ops
            .saturating_add(binary_method_group_ops(anchor))
            .saturating_add(run_end - run_start);
    }
    ops
}

fn descent_ops(trials: u64, r_bits: u32) -> u64 {
    let mul = expected_binary_method_group_ops(r_bits);
    let placement = 3 * mul + (WALKS - 1);
    if trials == 0 {
        return placement;
    }
    let steps = trials.saturating_sub(WALKS).div_ceil(WALKS);
    placement
        .saturating_add(steps.saturating_mul(WALKS))
        .saturating_add(mul)
}

fn ic_then_rho(n: u32, a: u8, target_index: usize, rep: u32) -> bool {
    (n + u32::from(a) + target_index as u32 + rep) % 2 == 0
}

fn run_rho(curve: &KoblitzCurve, target: &BinaryPoint, seed: u64) -> (bool, u64, u64, u64, u128) {
    let options = KoblitzSignedRhoOptions {
        seed,
        max_iterations_per_restart: 1 << 20,
        ..KoblitzSignedRhoOptions::default()
    };
    let started = Instant::now();
    let report = koblitz_signed_frobenius_rho_with_progress(curve, target, &options, &mut |_| {});
    let ns = started.elapsed().as_nanos();
    let r_bits = curve.subgroup_order.bits() as u32;
    let ops = rho_exclusive_group_ops(&report.charges, r_bits);
    let verified = report.verified
        && report.recovered_log.as_ref().is_some_and(|d| {
            d < &curve.subgroup_order && curve.mul(curve.generator(), d) == *target
        });
    (
        verified,
        ops,
        report.iterations,
        report.charges.walk_group_additions,
        ns,
    )
}

fn run_ic(
    curve: &KoblitzCurve,
    target: &BinaryPoint,
) -> Result<(bool, u64, serde_json::Value, u128), String> {
    let started = Instant::now();
    let floor = recipe_floor(curve.n);
    let fb = build_subgroup_orbit_factor_base(curve, FB_SEED, floor)
        .map_err(|e| format!("factor base: {e}"))?;
    let r = curve.subgroup_order.to_u64_digits()[0];
    let r_bits = curve.subgroup_order.bits() as u32;
    let k = fb.unknowns();
    let rows = PairSumTable::optimal_folded_rows(k, fb.points.len(), r);
    let window = fb.points.len().saturating_sub(1).max(1);
    let (pair, table_adds) =
        PairSumTable::build_folded_rows(curve, &fb, rows).ok_or("folded pair table")?;
    let opts = KoblitzIcOptions {
        m: 3,
        seed: ALGORITHM_SEED,
        max_trials: MAX_TRIALS,
        strategy: DecompositionStrategy::PairTable,
        allow_direct_relation: false,
        collection_window: Some(window),
        ..KoblitzIcOptions::default()
    };
    let collector = RelationCollector::with_pair_table(curve, &fb, &opts, Some(&pair))
        .ok_or("collector (base cannot decompose)")?;
    let mut solver = FactorBaseLogSolver::new(curve, &fb, &opts).ok_or("no projected columns")?;
    let mut trials = 0u64;
    let mut batches = 0u64;
    let mut outcome = None;
    while trials < MAX_TRIALS as u64 && outcome.is_none() {
        let count = (PRECOMPUTE_BATCH_TRIALS as u64).min(MAX_TRIALS as u64 - trials);
        let (rels, report) = collector.collect(RelationWorkUnit {
            seed: ALGORITHM_SEED,
            start: trials,
            count,
        });
        trials += report.trials as u64;
        batches += 1;
        solver.push(&rels);
        outcome = solver.try_solve();
    }
    let (table, log_report) = outcome.ok_or("log database incomplete")?;
    if !table.verify(curve) {
        return Err("column logs failed [x]G = R".into());
    }
    let descent =
        IndividualLogSolver::new(curve, &fb, &table, &opts, Some(&pair)).ok_or("descent setup")?;
    let (found, descent_report) = descent.solve(target).ok_or("descent failed")?;
    let verified = found < curve.subgroup_order && curve.mul(curve.generator(), &found) == *target;
    let cofactor = curve.cofactor.to_u64_digits().first().copied().unwrap_or(1);
    let cert: u64 = table.columns.iter().map(|(_, log)| log_ops(log)).sum();
    let ops = table_adds
        .saturating_add(
            u64::try_from(k)
                .unwrap_or(0)
                .saturating_mul(binary_method_group_ops(cofactor)),
        )
        .saturating_add(walked_collection_ops(ALGORITHM_SEED, trials, r, batches))
        .saturating_add(cert)
        .saturating_add(descent_ops(descent_report.trials as u64, r_bits));
    let ns = started.elapsed().as_nanos();
    let detail = json!({
        "factor_base_size": fb.points.len(),
        "orbits": k,
        "folded_rows": rows,
        "pair_table_additions": table_adds,
        "pair_table_entries": pair.len(),
        "collection_trials": trials,
        "relations": log_report.relations,
        "descent_trials": descent_report.trials,
        "window": window,
        "recipe_floor": floor,
    });
    Ok((verified, ops, detail, ns))
}

fn main() {
    let quick = env::args().any(|a| a == "--quick");
    let targets = if quick { 1usize } else { TARGET_SEEDS.len() };
    let reps = if quick { 1u32 } else { REPS };
    let cells: &[Cell] = if quick { &CELLS[..1] } else { &CELLS };
    let out = env::args()
        .position(|a| a == "--out")
        .and_then(|i| env::args().nth(i + 1))
        .map(PathBuf::from);

    println!("rho-parity protocol seed={ALGORITHM_SEED:#x} fb_seed={FB_SEED} quick={quick}");
    println!(
        "{:>8} {:>4} {:>4} {:>8} {:>10} {:>10} {:>8} {:>6} {:>6} {}",
        "cell", "tgt", "rep", "order", "G_ic", "G_rho", "alpha", "S_ic", "S_rho", "ok"
    );

    let mut rows = Vec::new();
    let mut all_ok = true;
    for cell in cells {
        let Some(curve) = KoblitzCurve::new(cell.a, cell.n) else {
            println!("{:>8} skip (no curve)", cell_name(cell));
            all_ok = false;
            continue;
        };
        let r = curve.subgroup_order.to_u64_digits()[0];
        let sqrt_r = (r as f64).sqrt();
        let expected_rho = rho_expected_steps(r, curve.n);
        for t in 0..targets {
            let target = hash_target(&curve, TARGET_SEEDS[t]);
            for rep in 0..reps {
                let ic_first = ic_then_rho(cell.n, cell.a, t, rep);
                let rho_seed = ALGORITHM_SEED ^ (t as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
                let (ic, rho) = if ic_first {
                    let ic = run_ic(&curve, &target);
                    let rho = run_rho(&curve, &target, rho_seed);
                    (ic, rho)
                } else {
                    let rho = run_rho(&curve, &target, rho_seed);
                    let ic = run_ic(&curve, &target);
                    (ic, rho)
                };
                let (rho_ok, g_rho, rho_iters, rho_walk, rho_ns) = rho;
                match ic {
                    Ok((ic_ok, g_ic, detail, ic_ns)) => {
                        let alpha = if g_rho == 0 {
                            f64::INFINITY
                        } else {
                            g_ic as f64 / g_rho as f64
                        };
                        let ok = ic_ok && rho_ok && alpha <= 1.0;
                        all_ok &= ic_ok && rho_ok;
                        println!(
                            "{:>8} {:>4} {:>4} {:>8} {:>10} {:>10} {:>8.4} {:>6.3} {:>6.3} {}",
                            cell_name(cell),
                            t,
                            rep,
                            if ic_first { "ic-rho" } else { "rho-ic" },
                            g_ic,
                            g_rho,
                            alpha,
                            g_ic as f64 / sqrt_r,
                            g_rho as f64 / sqrt_r,
                            if ok { "yes" } else { "NO" }
                        );
                        rows.push(json!({
                            "cell": cell_name(cell),
                            "a": cell.a,
                            "n": cell.n,
                            "subgroup_order": curve.subgroup_order.to_string(),
                            "target_index": t,
                            "target_seed": TARGET_SEEDS[t],
                            "rep": rep,
                            "order": if ic_first { "ic-rho" } else { "rho-ic" },
                            "ic_verified": ic_ok,
                            "rho_verified": rho_ok,
                            "g_ic": g_ic,
                            "g_rho": g_rho,
                            "alpha": alpha,
                            "s_ic": g_ic as f64 / sqrt_r,
                            "s_rho": g_rho as f64 / sqrt_r,
                            "rho_iterations": rho_iters,
                            "rho_walk_additions": rho_walk,
                            "rho_expected_steps": expected_rho,
                            "ic_ns": ic_ns,
                            "rho_ns": rho_ns,
                            "ic": detail,
                            "gate": ok,
                        }));
                    }
                    Err(reason) => {
                        all_ok = false;
                        println!(
                            "{:>8} {:>4} {:>4} {:>8} ic failed: {reason}",
                            cell_name(cell),
                            t,
                            rep,
                            if ic_first { "ic-rho" } else { "rho-ic" }
                        );
                        rows.push(json!({
                            "cell": cell_name(cell),
                            "a": cell.a,
                            "n": cell.n,
                            "target_index": t,
                            "rep": rep,
                            "ic_verified": false,
                            "rho_verified": rho_ok,
                            "g_rho": g_rho,
                            "error": reason,
                            "gate": false,
                        }));
                    }
                }
            }
        }
    }

    let useful: Vec<&serde_json::Value> = rows
        .iter()
        .filter(|r| r["ic_verified"] == true && r["rho_verified"] == true)
        .collect();
    let mut cell_means = Vec::new();
    for cell in cells {
        let name = cell_name(cell);
        let xs: Vec<f64> = useful
            .iter()
            .filter(|r| r["cell"] == name)
            .filter_map(|r| r["alpha"].as_f64())
            .collect();
        if xs.is_empty() {
            cell_means.push(json!({"cell": name, "pairs": 0, "mean_alpha": null, "max_alpha": null, "gate": false}));
            continue;
        }
        let mean = xs.iter().sum::<f64>() / xs.len() as f64;
        let max = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        cell_means.push(json!({
            "cell": name,
            "pairs": xs.len(),
            "mean_alpha": mean,
            "max_alpha": max,
            "gate": mean <= 1.0 && max <= 1.0,
        }));
    }
    let gate = cell_means.iter().all(|c| c["gate"] == true) && all_ok && !quick;
    let summary = json!({
        "protocol": "RESEARCH_ECC2K130_RHO_PARITY.md",
        "algorithm_seed": ALGORITHM_SEED,
        "fb_seed": FB_SEED,
        "target_seeds": TARGET_SEEDS,
        "quick": quick,
        "unit": "exclusive group operations (adds + binary-method scalar muls)",
        "class": "engineering",
        "n131_claim": false,
        "all_cases_gate": gate,
        "cells": cell_means,
        "rows": rows,
    });
    println!("{}", serde_json::to_string_pretty(&summary).unwrap());
    if let Some(dir) = out {
        fs::create_dir_all(&dir).expect("out dir");
        fs::write(
            dir.join("summary.json"),
            serde_json::to_vec_pretty(&summary).unwrap(),
        )
        .expect("summary");
    }
}
