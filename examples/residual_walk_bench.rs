//! **Partial-decomposition residual walks — measurement panel.**
//!
//! Companion to `crypto_lib::cryptanalysis::residual_walk` and
//! `RESEARCH_RESIDUAL_WALKS.md`.  Compares five ways of turning residual
//! collisions into index-calculus relations on toy prime-order curves,
//! all sharing one factor base, one operation counter, one verifier and
//! one rank tracker:
//!
//! - `A`  independent random partial sums `(a, b, k-tuple)`;
//! - `B`  local-mutation walk (swap one factor-base point per step);
//! - `C1` collision-preserving r-adding walk on residuals;
//! - `C2` collision-preserving fresh-hash walk `s_{t+1} = H(L(s_t))`;
//! - `R`  plain Pollard rho, the reference.
//!
//! ```text
//! cargo run --release --example residual_walk_bench -- --bits 24 --fb 256 --trials 3
//! cargo run --release --example residual_walk_bench -- --panel --json experiments/20_residual_walk_panel.json
//! cargo run --release --example residual_walk_bench -- --panel --quick
//! cargo run --release --example residual_walk_bench -- --baseline --json experiments/20_residual_walk_baseline.json
//! ```
//!
//! `--baseline` is the fixed optimisation protocol (`n ≈ 2^24` and
//! `2^28`, `B = 256`, `k = 3`, seeds 1–3, every strategy, plus `C1` and
//! `R` with 8 distinguished-point bits; about a minute).  Score it, or
//! compare it against the frozen baseline, with
//! `scripts/residual_walk_scoreboard.py`.
//!
//! Every run verifies each relation by scalar multiplication and scores
//! the recovered logarithm against the planted one.

use std::env;
use std::fs;

use crypto_lib::cryptanalysis::residual_walk::{
    generate_instance, mitm_four_decomposition, run_strategy, FactorBase, Instance, MitmReport,
    Strategy, StrategyReport, WalkOptions,
};
use serde::Serialize;

#[derive(Serialize)]
struct Panel {
    generated_by: &'static str,
    recovery: Vec<StrategyReport>,
    fixed_budget: Vec<StrategyReport>,
    rank_cap: Vec<StrategyReport>,
    filter_trap: Vec<StrategyReport>,
    distinguished_points: Vec<StrategyReport>,
    mitm: Vec<MitmReport>,
}

fn header() {
    println!(
        "{:<3} {:>4} {:>5} {:>4} {:>3} {:>10} {:>9} {:>11} {:>11} {:>6} {:>6} {:>5} {:>5} {:>4} {:>8} {:>7}",
        "tag", "bits", "B", "seed", "dp", "samples", "table", "total_ops", "ops/rel", "indep",
        "rank", "full", "triv", "ok", "ops/rho", "ms"
    );
}

fn row(r: &StrategyReport) {
    println!(
        "{:<3} {:>4} {:>5} {:>4} {:>3} {:>10} {:>9} {:>11} {:>11} {:>6} {:>6} {:>5} {:>5} {:>4} {:>8.2} {:>7.0}",
        r.tag,
        r.bits,
        r.factor_base,
        "-",
        r.dp_bits,
        r.samples,
        r.table_entries,
        r.total_ops,
        r.ops_per_independent_relation
            .map(|v| format!("{v:.0}"))
            .unwrap_or_else(|| "-".into()),
        r.relations_independent,
        r.rank,
        r.full_decompositions,
        r.collisions_trivial,
        match r.correct {
            Some(true) => "yes",
            Some(false) => "WRONG",
            None => "no",
        },
        r.total_ops as f64 / r.predicted_rho_steps,
        r.wall_ms,
    );
}

fn mitm_row(r: &MitmReport) {
    println!(
        "MITM bits={} B={} pairs={} coincidences={} targets={} rel={} indep={} rank={} ops={} ops/rel={} pred_matches/target={:.2} ok={:?} ms={:.0}",
        r.bits,
        r.factor_base,
        r.pair_table,
        r.pair_coincidences,
        r.targets,
        r.relations_verified,
        r.relations_independent,
        r.rank,
        r.total_ops,
        r.ops_per_independent_relation
            .map(|v| format!("{v:.0}"))
            .unwrap_or_else(|| "-".into()),
        r.predicted_matches_per_target,
        r.correct,
        r.wall_ms
    );
}

fn setup(bits: u32, fb: usize, seed: u64) -> (Instance, FactorBase) {
    let inst = generate_instance(bits, seed);
    let fb = FactorBase::build(&inst.curve, fb);
    (inst, fb)
}

fn run_all(
    inst: &Instance,
    fb: &FactorBase,
    strategies: &[Strategy],
    opts: &WalkOptions,
    out: &mut Vec<StrategyReport>,
) {
    for &s in strategies {
        let r = run_strategy(inst, fb, s, opts);
        row(&r);
        out.push(r);
    }
}

fn panel(quick: bool) -> Panel {
    let mut p = Panel {
        generated_by: "examples/residual_walk_bench.rs --panel",
        recovery: Vec::new(),
        fixed_budget: Vec::new(),
        rank_cap: Vec::new(),
        filter_trap: Vec::new(),
        distinguished_points: Vec::new(),
        mitm: Vec::new(),
    };
    let budget = 1u64 << 34;

    // ── P1: full target recovery ───────────────────────────────────
    println!("\n== P1: operations to recover d (stop at first solve) ==");
    header();
    let mut cases: Vec<(u32, usize, u64)> = Vec::new();
    let sizes: &[u32] = if quick { &[20, 24] } else { &[20, 24, 28] };
    let seeds: u64 = if quick { 1 } else { 3 };
    for &bits in sizes {
        for &fbsize in &[64usize, 256] {
            for seed in 1..=seeds {
                cases.push((bits, fbsize, seed));
            }
        }
    }
    if !quick {
        cases.push((28, 1024, 1));
        cases.push((32, 256, 1));
        cases.push((32, 1024, 1));
    }
    for (bits, fbsize, seed) in cases {
        let (inst, fb) = setup(bits, fbsize, seed);
        println!(
            "-- n = {} (2^{:.1}), B = {}, seed {}",
            inst.curve.n,
            (inst.curve.n as f64).log2(),
            fb.len(),
            seed
        );
        let opts = WalkOptions {
            k: 3,
            max_ops: budget,
            seed,
            ..WalkOptions::default()
        };
        run_all(&inst, &fb, &Strategy::ALL, &opts, &mut p.recovery);
    }

    // ── P2: relation yield at a fixed budget ───────────────────────
    println!("\n== P2: relations at a fixed budget (no early stop) ==");
    header();
    let (bits, fbsize) = if quick { (22, 128) } else { (26, 256) };
    let (inst, fb) = setup(bits, fbsize, 1);
    let fixed = (2.0 * inst.curve.n as f64 * (fb.len() as f64 + 1.0)).sqrt() as u64 * 8;
    println!(
        "-- n = {}, B = {}, budget {} ops",
        inst.curve.n,
        fb.len(),
        fixed
    );
    let opts = WalkOptions {
        k: 3,
        max_ops: fixed,
        seed: 1,
        stop_when_solved: false,
        ..WalkOptions::default()
    };
    run_all(&inst, &fb, &Strategy::ALL, &opts, &mut p.fixed_budget);

    // ── P3: rank cap of the r-adding residual walk ─────────────────
    println!("\n== P3: r-adding residual walk, rank vs number of multipliers ==");
    header();
    let (inst, fb) = setup(if quick { 20 } else { 24 }, 256, 2);
    for r in [8usize, 32, 128, 256] {
        let opts = WalkOptions {
            multipliers: r,
            max_ops: 1 << 24,
            seed: 2,
            stop_when_solved: false,
            ..WalkOptions::default()
        };
        let rep = run_strategy(&inst, &fb, Strategy::RAddingResidualWalk, &opts);
        println!("r = {r}:");
        row(&rep);
        p.rank_cap.push(rep);
    }

    // ── P4: the rejection-filter trap ──────────────────────────────
    println!("\n== P4: independent samples with residuals filtered to x < M ==");
    header();
    let (inst, fb) = setup(if quick { 20 } else { 24 }, 64, 3);
    let pp = inst.curve.p;
    for shift in [0u32, 2, 4, 6, 8] {
        let opts = WalkOptions {
            k: 3,
            max_ops: budget,
            seed: 3,
            filter_bound: if shift == 0 { None } else { Some(pp >> shift) },
            ..WalkOptions::default()
        };
        let rep = run_strategy(&inst, &fb, Strategy::IndependentSamples, &opts);
        println!("M = p / 2^{shift}:");
        row(&rep);
        p.filter_trap.push(rep);
    }

    // ── P5: distinguished points on the collision-preserving walks ─
    println!("\n== P5: distinguished-point storage (memory vs operations) ==");
    header();
    let (inst, fb) = setup(if quick { 22 } else { 28 }, 256, 4);
    for dp in [0u32, 4, 8] {
        for s in [
            Strategy::RAddingResidualWalk,
            Strategy::FreshHashWalk,
            Strategy::PlainRho,
        ] {
            let opts = WalkOptions {
                k: 3,
                dp_bits: dp,
                max_ops: budget,
                seed: 4,
                ..WalkOptions::default()
            };
            let rep = run_strategy(&inst, &fb, s, &opts);
            row(&rep);
            p.distinguished_points.push(rep);
        }
    }

    // ── P6: meet in the middle on 4-decompositions ─────────────────
    println!("\n== P6: meet-in-the-middle 4-decompositions ==");
    let sizes: &[u32] = if quick { &[20] } else { &[20, 24, 28] };
    for &bits in sizes {
        let inst = generate_instance(bits, 5);
        let base = (4.0 * inst.curve.n as f64).powf(0.25).ceil() as usize;
        for mult in [1usize, 2] {
            let fb = FactorBase::build(&inst.curve, base * mult);
            let rep = mitm_four_decomposition(&inst, &fb, 100_000, 5);
            mitm_row(&rep);
            p.mitm.push(rep);
        }
    }
    p
}

/// The fixed optimisation protocol: two sizes, one factor base, three
/// seeds, every strategy, plus the distinguished-point variants of the
/// two one-op-per-step walks.  Deterministic for a given build.
fn baseline() -> Vec<StrategyReport> {
    let mut out = Vec::new();
    header();
    for bits in [24u32, 28] {
        for seed in 1..=3u64 {
            let (inst, fb) = setup(bits, 256, seed);
            println!(
                "-- n = {} (2^{:.1}), B = {}, seed {}",
                inst.curve.n,
                (inst.curve.n as f64).log2(),
                fb.len(),
                seed
            );
            let opts = WalkOptions {
                k: 3,
                max_ops: 1 << 34,
                seed,
                ..WalkOptions::default()
            };
            run_all(&inst, &fb, &Strategy::ALL, &opts, &mut out);
            let dp = WalkOptions {
                dp_bits: 8,
                ..opts.clone()
            };
            run_all(
                &inst,
                &fb,
                &[Strategy::RAddingResidualWalk, Strategy::PlainRho],
                &dp,
                &mut out,
            );
        }
    }
    out
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut bits = 24u32;
    let mut fbsize = 256usize;
    let mut k = 3usize;
    let mut trials = 1u64;
    let mut seed = 1u64;
    let mut dp = 0u32;
    let mut budget = 1u64 << 34;
    let mut strategies: Vec<Strategy> = Strategy::ALL.to_vec();
    let mut json: Option<String> = None;
    let mut do_panel = false;
    let mut do_baseline = false;
    let mut quick = false;
    let mut i = 0;
    let value = |i: &mut usize, args: &[String]| -> String {
        *i += 1;
        args.get(*i).cloned().unwrap_or_else(|| {
            eprintln!("missing value for {}", args[*i - 1]);
            std::process::exit(2);
        })
    };
    while i < args.len() {
        match args[i].as_str() {
            "--bits" => bits = value(&mut i, &args).parse().expect("--bits"),
            "--fb" => fbsize = value(&mut i, &args).parse().expect("--fb"),
            "--k" => k = value(&mut i, &args).parse().expect("--k"),
            "--trials" => trials = value(&mut i, &args).parse().expect("--trials"),
            "--seed" => seed = value(&mut i, &args).parse().expect("--seed"),
            "--dp" => dp = value(&mut i, &args).parse().expect("--dp"),
            "--budget" => budget = value(&mut i, &args).parse().expect("--budget"),
            "--strategies" => {
                strategies = value(&mut i, &args)
                    .split(',')
                    .map(|s| {
                        Strategy::parse(s.trim()).unwrap_or_else(|| {
                            eprintln!("unknown strategy {s}; use A,B,C1,C2,R");
                            std::process::exit(2);
                        })
                    })
                    .collect();
            }
            "--json" => json = Some(value(&mut i, &args)),
            "--panel" => do_panel = true,
            "--baseline" => do_baseline = true,
            "--quick" => quick = true,
            "--help" | "-h" => {
                println!(
                    "usage: residual_walk_bench [--bits N] [--fb B] [--k K] [--trials T] [--seed S] \
                     [--dp BITS] [--budget OPS] [--strategies A,B,C1,C2,R] [--json FILE] [--panel [--quick]] [--baseline]"
                );
                return;
            }
            other => {
                eprintln!("unknown argument {other}");
                std::process::exit(2);
            }
        }
        i += 1;
    }

    if do_baseline {
        let reports = baseline();
        if let Some(path) = json {
            fs::write(&path, serde_json::to_string_pretty(&reports).unwrap()).expect("write json");
            println!("\nwrote {path}");
        }
        return;
    }

    if do_panel {
        let p = panel(quick);
        if let Some(path) = json {
            fs::write(&path, serde_json::to_string_pretty(&p).unwrap()).expect("write json");
            println!("\nwrote {path}");
        }
        return;
    }

    let mut reports = Vec::new();
    header();
    for t in 0..trials {
        let (inst, fb) = setup(bits, fbsize, seed + t);
        println!(
            "-- p = {}, n = {} (2^{:.1}), B = {}, d = {}",
            inst.curve.p,
            inst.curve.n,
            (inst.curve.n as f64).log2(),
            fb.len(),
            inst.d
        );
        let opts = WalkOptions {
            k,
            max_ops: budget,
            dp_bits: dp,
            seed: seed + t,
            ..WalkOptions::default()
        };
        for &s in &strategies {
            let o = if s.collision_preserving() {
                opts.clone()
            } else {
                WalkOptions {
                    dp_bits: 0,
                    ..opts.clone()
                }
            };
            let r = run_strategy(&inst, &fb, s, &o);
            row(&r);
            reports.push(r);
        }
    }
    if let Some(path) = json {
        fs::write(&path, serde_json::to_string_pretty(&reports).unwrap()).expect("write json");
        println!("wrote {path}");
    }
}
