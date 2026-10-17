//! **E19 of `RESEARCH_GLV_INVARIANT_FACTOR_BASES.md` (§9.5): the first
//! end-to-end `k = 4` index calculus with the weighted `D₄` solve.**
//!
//! On E18's 2-torsion instances and `τ_T`-folded base, every target is
//! decomposed by `d4_weighted_s5`'s weighted solve at degree 26, one
//! verified decomposition per target becomes a folded row, and Wiedemann
//! runs until `[d]G = Q` verifies.  Every phase is charged in `F_p`
//! multiplications against `16` negation-folded walks on the same curve.
//!
//! ```text
//! # one curve per invocation survives a restart of the host
//! cargo run --release --example quartic_pipeline_d4 -- --sizes 53 --first-seed 1 --seeds 1 --json experiments/27_quartic_pipeline_d4_53s1.json
//! ```
//!
//! Tables: `cargo run --release --example glv_invariant_experiment_tables
//! -- experiments/27_quartic_pipeline_d4_*.json`.

use std::env;
use std::fs;

use crypto_lib::cryptanalysis::quartic_folds::run_k4_d4_pipeline;
use serde_json::Value;

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut sizes: Vec<u64> = vec![53, 101, 149, 197, 269];
    let mut seeds = 2u64;
    let mut first_seed = 1u64;
    let mut rho_runs = 16usize;
    let mut max_residuals = 200_000u64;
    let mut json: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--sizes" => {
                i += 1;
                sizes = args[i]
                    .split(',')
                    .map(|v| v.parse().expect("--sizes"))
                    .collect();
            }
            "--first-seed" => {
                i += 1;
                first_seed = args[i].parse().expect("--first-seed");
            }
            "--seeds" => {
                i += 1;
                seeds = args[i].parse().expect("--seeds");
            }
            "--rho-runs" => {
                i += 1;
                rho_runs = args[i].parse().expect("--rho-runs");
            }
            "--max-residuals" => {
                i += 1;
                max_residuals = args[i].parse().expect("--max-residuals");
            }
            "--json" => {
                i += 1;
                json = Some(args[i].clone());
            }
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    let mut rows: Vec<Value> = Vec::new();
    for &p in &sizes {
        for seed in first_seed..first_seed + seeds {
            let r = match run_k4_d4_pipeline(p, seed, rho_runs, max_residuals) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("e19 p={p} seed={seed}: {e}");
                    continue;
                }
            };
            let rel = r.folded.relations.max(1) as f64;
            let pct = |v: u64| 100.0 * v as f64 / r.total_muls as f64;
            println!(
                "e19 p={} seed={} n=2^{:.1} h={} |F|={} cols {} | residuals {} unsolved {} reruns {} | decs {} off-base {} failures {} | rel {} y {:.1} core {} φ {:.3} | C median {:.3e} mean {:.3e} | solve {:.1}% LA {:.2}% other {:.2}% | S {:.3e} rho(neg) {:.3e} ratio {:.3e} | correct {} rho {} | {:.0} s",
                r.p, r.seed, r.bits, r.cofactor, r.base_points, r.folded.columns,
                r.residuals, r.unsolved, r.reruns,
                r.decompositions, r.off_base, r.check_failures,
                r.folded.relations, r.residuals as f64 / rel, r.folded.unknowns, r.folded.phi,
                r.call.median, r.call.mean,
                pct(r.weil_muls + r.solve_muls), pct(r.la_fp_muls),
                pct(r.precompute_muls + r.base_muls + r.generation_muls + r.lift_muls + r.check_muls),
                r.total_muls as f64, r.rho_matched_muls, r.ratio,
                r.folded.correct, r.rho.iter().all(|w| w.correct),
                r.wall_ms / 1e3
            );
            let mut v = serde_json::to_value(&r).unwrap();
            v["experiment"] = Value::from("e19");
            rows.push(v);
            if let Some(path) = &json {
                fs::write(path, serde_json::to_string_pretty(&rows).unwrap()).expect("write json");
            }
        }
    }
}
