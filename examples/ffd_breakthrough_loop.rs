//! # FFD breakthrough loop — the executable shell of `RESEARCH_FFD_WORKFLOW.md`.
//!
//! One invocation = one snapshot of the falsification-driven loop:
//!
//! 1. Runs the currently-decidable experiments (EXP-A curve-independence,
//!    the operating-point scaling sweep for P1, the regime sweep for P6).
//! 2. Applies the **pre-registered decision gates** (G-P5, G-P1, G-P6) to
//!    the fresh data and prints a verdict per prediction.
//! 3. Writes a reproducible JSON snapshot to `experiments/ffd_loop_<ts>.json`
//!    so the numbers survive the ephemeral container and feed the ledger.
//! 4. Prints the **next experiment** to run (the queue head whose gate is
//!    closest to decidable but still blocked).
//!
//! This does not *replace* judgement — it mechanises steps 0, 3, and the
//! arithmetic of step 4 in the workflow, so a human (or the autolab) starts
//! each iteration from a verdict, not from raw numbers.
//!
//! ```bash
//! cargo run --release --example ffd_breakthrough_loop
//! ```

use crypto_lib::cryptanalysis::pc_degree_avg::{
    mean_dstar_spread, run_pc_curve_independence_sweep, run_pc_operating_point_sweep,
    run_pc_regime_sweep, AvgRow,
};
use std::time::{SystemTime, UNIX_EPOCH};

/// Pre-registered gate thresholds (mirror RESEARCH_FFD_WORKFLOW.md §4).
const G_P5_SPREAD_SUPPORT: f64 = 0.5; // spread below this ⇒ curve-independence supported
const G_P6_OVERDET_RATIO: f64 = 2.0; // ρ ≥ this should collapse D* to ≈ 2
const G_P6_COLLAPSE_MEAN: f64 = 2.05; // mean D* below this counts as "collapsed"

fn main() {
    let seed = 0x_FFD_10_09; // fixed for reproducibility of this snapshot
    let targets = 64;

    println!("════════════════════════════════════════════════════════════════");
    println!(" FFD breakthrough loop — snapshot");
    println!(" (gates: RESEARCH_FFD_WORKFLOW.md §4)");
    println!("════════════════════════════════════════════════════════════════");

    // ── EXP-A: curve-independence (P5) ──────────────────────────────
    let n_a = 8;
    let nsub_a = 3;
    let n_curves = 8;
    let exp_a = run_pc_curve_independence_sweep(n_a, nsub_a, n_curves, 12, targets, seed);
    let spread = mean_dstar_spread(&exp_a);
    let p5 = judge_p5(spread);
    println!(
        "\n[EXP-A] curve-independence  n={n_a} n'={nsub_a}  ({n_curves} curves)"
    );
    for (i, r) in exp_a.iter().enumerate() {
        println!(
            "   curve #{i}: mean D*={}  hist {}",
            fmt_mean(r),
            r.refutation.histogram_str()
        );
    }
    println!(
        "   → spread(mean D*) = {}   GATE G-P5: {}",
        spread.map(|s| format!("{s:.3}")).unwrap_or_else(|| "—".into()),
        p5
    );

    // ── P1: operating-point scaling, split by parity ────────────────
    let opp = run_pc_operating_point_sweep(6..=12, 1, 16, targets, seed.wrapping_add(1));
    println!("\n[EXP/scale] operating-point (margin=1, ρ≳1), split by parity:");
    print_table(&opp);
    let odd: Vec<&AvgRow> = opp.iter().filter(|r| r.n % 2 == 1).collect();
    let even: Vec<&AvgRow> = opp.iter().filter(|r| r.n % 2 == 0).collect();
    let p1_odd = judge_p1(&odd, "odd");
    let p1_even = judge_p1(&even, "even");
    println!("   GATE G-P1 (odd  parity): {p1_odd}");
    println!("   GATE G-P1 (even parity): {p1_even}");

    // ── P6: over-determination collapse ─────────────────────────────
    let regime = run_pc_regime_sweep(10, 14, targets, seed.wrapping_add(2));
    println!("\n[EXP/regime] determination-ratio at fixed n=10:");
    print_table(&regime);
    let p6 = judge_p6(&regime);
    println!("   GATE G-P6 (over-determination ⇒ D*→2): {p6}");

    // ── Snapshot to experiments/ ────────────────────────────────────
    // Fixed seed ⇒ reproducible run ⇒ a single canonical file (overwritten
    // each run) rather than a timestamped pile. `generated_at` records the
    // wall clock inside the JSON for provenance.
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let path = "experiments/ffd_loop_latest.json".to_string();
    let json = build_json(seed, targets, ts, &exp_a, spread, &p5, &opp, &p1_odd, &p1_even, &regime, &p6);
    match std::fs::write(&path, &json) {
        Ok(_) => println!("\n[snapshot] wrote {path}"),
        Err(e) => println!("\n[snapshot] FAILED to write {path}: {e}"),
    }

    // ── Next experiment (queue head still blocked) ──────────────────
    println!("\n[next] {}", next_experiment(&p5, &p1_odd, &p1_even));
    println!("════════════════════════════════════════════════════════════════");
}

// ── Gate judgements ─────────────────────────────────────────────────

#[derive(Clone)]
struct Verdict {
    status: &'static str, // supported | killed | blocked | inconclusive
    detail: String,
}
impl std::fmt::Display for Verdict {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ({})", self.status.to_uppercase(), self.detail)
    }
}

/// G-P5: small mean-D* spread across curves ⇒ curve-independence supported.
fn judge_p5(spread: Option<f64>) -> Verdict {
    match spread {
        None => Verdict {
            status: "blocked",
            detail: "fewer than 2 defined means".into(),
        },
        Some(s) if s < G_P5_SPREAD_SUPPORT => Verdict {
            status: "supported",
            detail: format!("spread {s:.3} < {G_P5_SPREAD_SUPPORT}"),
        },
        Some(s) => Verdict {
            status: "killed",
            detail: format!("spread {s:.3} ≥ {G_P5_SPREAD_SUPPORT}: D* depends on curve"),
        },
    }
}

/// G-P1: positive slope of mean D* vs n (within a parity class), ≥ 3 points.
fn judge_p1(rows: &[&AvgRow], label: &str) -> Verdict {
    let pts: Vec<(f64, f64)> = rows
        .iter()
        .filter_map(|r| r.refutation.mean.map(|m| (r.n as f64, m)))
        .collect();
    if pts.len() < 3 {
        return Verdict {
            status: "blocked",
            detail: format!("{label}: only {} points (need ≥3)", pts.len()),
        };
    }
    let slope = ols_slope(&pts);
    if slope > 0.02 {
        Verdict {
            status: "supported",
            detail: format!("{label}: slope {slope:+.3}/n over {} pts", pts.len()),
        }
    } else if slope.abs() <= 0.02 {
        Verdict {
            status: "inconclusive",
            detail: format!("{label}: flat slope {slope:+.3} (reach-limited)"),
        }
    } else {
        Verdict {
            status: "killed",
            detail: format!("{label}: negative slope {slope:+.3}"),
        }
    }
}

/// G-P6: rows with ρ ≥ 2 should have mean D* ≈ 2.
fn judge_p6(rows: &[AvgRow]) -> Verdict {
    let overdet: Vec<&AvgRow> = rows.iter().filter(|r| r.ratio() >= G_P6_OVERDET_RATIO).collect();
    if overdet.is_empty() {
        return Verdict {
            status: "blocked",
            detail: "no ρ ≥ 2 cells in sweep".into(),
        };
    }
    let all_collapsed = overdet
        .iter()
        .all(|r| r.refutation.mean.map_or(false, |m| m <= G_P6_COLLAPSE_MEAN));
    if all_collapsed {
        Verdict {
            status: "supported",
            detail: format!("all {} over-determined cells have mean D* ≤ {G_P6_COLLAPSE_MEAN}", overdet.len()),
        }
    } else {
        Verdict {
            status: "killed",
            detail: "an over-determined cell did NOT collapse to D*≈2".into(),
        }
    }
}

/// Ordinary-least-squares slope of y on x.
fn ols_slope(pts: &[(f64, f64)]) -> f64 {
    let n = pts.len() as f64;
    let sx: f64 = pts.iter().map(|p| p.0).sum();
    let sy: f64 = pts.iter().map(|p| p.1).sum();
    let sxx: f64 = pts.iter().map(|p| p.0 * p.0).sum();
    let sxy: f64 = pts.iter().map(|p| p.0 * p.1).sum();
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-12 {
        0.0
    } else {
        (n * sxy - sx * sy) / denom
    }
}

/// Pick the next experiment to run: the highest-priority gate that is not
/// yet supported (mirrors the §5 queue, re-prioritised from verdicts).
fn next_experiment(p5: &Verdict, p1_odd: &Verdict, p1_even: &Verdict) -> String {
    if p5.status == "blocked" {
        return "EXP-A: widen curve-independence batch (need ≥2 defined means)".into();
    }
    if p5.status == "killed" {
        return "EXP-C: build descent_expansion KEYED ON THE CURVE (P5 killed: D* is curve-dependent)".into();
    }
    // P5 supported → expansion predictor can be field/basis-keyed.
    if p1_odd.status == "blocked" || p1_even.status == "blocked" {
        return "EXP-D: build sparse-F4 backend to reach 2n'≳14 (G-P1 reach-limited)".into();
    }
    if p1_odd.status == "inconclusive" || p1_even.status == "inconclusive" {
        return "EXP-D: extend operating-point reach; current slope flat within noise".into();
    }
    // P1 has a verdict on both parities → the central test is next.
    "EXP-C/E: build descent_expansion (field/basis-keyed) + basis sweep, then run G-P2".into()
}

// ── Output helpers ──────────────────────────────────────────────────

fn fmt_mean(r: &AvgRow) -> String {
    r.refutation
        .mean
        .map(|m| format!("{m:.3}"))
        .unwrap_or_else(|| "—".into())
}

fn print_table(rows: &[AvgRow]) {
    println!("    {:>3} {:>3} {:>5} {:>5} {:>9} {:>12}", "n", "n'", "ρ", "N", "mean D*", "hist");
    for r in rows {
        println!(
            "    {:>3} {:>3} {:>5.2} {:>5} {:>9} {:>12}",
            r.n,
            r.n_sub,
            r.ratio(),
            r.n_targets,
            fmt_mean(r),
            r.refutation.histogram_str()
        );
    }
}

/// Hand-rolled JSON (avoids requiring serde derives on the harness types;
/// the snapshot is flat and stable).
#[allow(clippy::too_many_arguments)]
fn build_json(
    seed: u64,
    targets: u32,
    generated_at: u64,
    exp_a: &[AvgRow],
    spread: Option<f64>,
    p5: &Verdict,
    opp: &[AvgRow],
    p1_odd: &Verdict,
    p1_even: &Verdict,
    regime: &[AvgRow],
    p6: &Verdict,
) -> String {
    let row_json = |r: &AvgRow| {
        format!(
            "{{\"n\":{},\"n_sub\":{},\"vars\":{},\"eqs\":{},\"ratio\":{:.4},\"n_targets\":{},\"mean_dstar\":{},\"min_dstar\":{},\"max_dstar\":{},\"none\":{},\"hist\":\"{}\"}}",
            r.n,
            r.n_sub,
            r.num_vars,
            r.num_eqs,
            r.ratio(),
            r.n_targets,
            r.refutation.mean.map(|m| format!("{m:.4}")).unwrap_or_else(|| "null".into()),
            r.refutation.min.map(|v| v.to_string()).unwrap_or_else(|| "null".into()),
            r.refutation.max.map(|v| v.to_string()).unwrap_or_else(|| "null".into()),
            r.refutation.n_none,
            r.refutation.histogram_str(),
        )
    };
    let arr = |rows: &[AvgRow]| {
        rows.iter().map(row_json).collect::<Vec<_>>().join(",")
    };
    let verdict = |v: &Verdict| format!("{{\"status\":\"{}\",\"detail\":\"{}\"}}", v.status, v.detail.replace('"', "'"));
    format!(
        "{{\n  \"schema\": \"ffd_breakthrough_loop/v1\",\n  \"seed\": {seed},\n  \"generated_at\": {generated_at},\n  \"targets_per_cell\": {targets},\n  \"exp_a_curve_independence\": {{\n    \"spread_mean_dstar\": {},\n    \"gate_p5\": {},\n    \"rows\": [{}]\n  }},\n  \"scaling_operating_point\": {{\n    \"gate_p1_odd\": {},\n    \"gate_p1_even\": {},\n    \"rows\": [{}]\n  }},\n  \"regime_n10\": {{\n    \"gate_p6\": {},\n    \"rows\": [{}]\n  }}\n}}\n",
        spread.map(|s| format!("{s:.4}")).unwrap_or_else(|| "null".into()),
        verdict(p5),
        arr(exp_a),
        verdict(p1_odd),
        verdict(p1_even),
        arr(opp),
        verdict(p6),
        arr(regime),
    )
}
