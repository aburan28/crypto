//! **Joux–Vitse three-point decompositions at `k = 4`: `C′` and the method
//! end to end.**
//!
//! Companion to `research/notes/index-calculus/RESEARCH_RHO_PARITY_PROGRAMME.md`.
//!
//! ```text
//! cargo run --release --example jv_quartic -- --exp cprime --sizes 269,521,769,1033 --seeds 2 --residuals 200 --constructed 40 --json experiments/26_jv_quartic_cprime.json
//! cargo run --release --example jv_quartic -- --exp dlp --sizes 269,521,769,1033 --seeds 2 --rho-runs 16 --check-every 256 --json experiments/26_jv_quartic_dlp.json
//! ```

use std::env;
use std::fs;

use crypto_lib::cryptanalysis::jv_quartic::{
    run_jv4_cprime, run_jv4_dlp, CprimeReport, Jv4DlpReport,
};
use serde::Serialize;

#[derive(Serialize)]
#[serde(untagged)]
#[allow(clippy::large_enum_variant)]
enum Row {
    Cprime(CprimeReport),
    Dlp(Jv4DlpReport),
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut exp = "cprime".to_string();
    let mut sizes: Vec<u64> = vec![269, 521, 769, 1033];
    let mut seeds = 2u64;
    let mut first_seed = 1u64;
    let mut residuals = 200usize;
    let mut constructed = 40usize;
    let mut rho_runs = 16usize;
    let mut check_every = 256u64;
    let mut json: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--exp" => {
                i += 1;
                exp = args[i].clone();
            }
            "--sizes" => {
                i += 1;
                sizes = args[i]
                    .split(',')
                    .map(|v| v.parse().expect("--sizes"))
                    .collect();
            }
            "--seeds" => {
                i += 1;
                seeds = args[i].parse().expect("--seeds");
            }
            "--first-seed" => {
                i += 1;
                first_seed = args[i].parse().expect("--first-seed");
            }
            "--residuals" => {
                i += 1;
                residuals = args[i].parse().expect("--residuals");
            }
            "--constructed" => {
                i += 1;
                constructed = args[i].parse().expect("--constructed");
            }
            "--rho-runs" => {
                i += 1;
                rho_runs = args[i].parse().expect("--rho-runs");
            }
            "--check-every" => {
                i += 1;
                check_every = args[i].parse().expect("--check-every");
            }
            "--json" => {
                i += 1;
                json = Some(args[i].clone());
            }
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    let mut rows: Vec<Row> = Vec::new();
    for &p in &sizes {
        for seed in first_seed..first_seed + seeds {
            match exp.as_str() {
                "cprime" => {
                    let r = run_jv4_cprime(p, seed, residuals, constructed);
                    println!(
                        "p={:>5} seed={} n=2^{:.1} |F|={:>4} | random {} decomposable {} rate {:.5} (1/6p = {:.5}) | planted {}/{} | mismatches {} undetermined {} | C' {:.0} [{}, {}] = weil {:.0} + f4 {:.0} + roots {:.0} + signs {:.1} ops | F4 degree {:.2} [{}, {}] rows {:.0} cols {:.0} {:.3} ms | decomposable C' {:.0} degree {:.2} | mitm3 {:.0} ops | {:.1} s",
                        r.p, r.seed, r.bits, r.base, r.random_residuals, r.random_decomposable,
                        r.decomposition_rate, r.expected_rate, r.planted_found, r.constructed_residuals,
                        r.mismatches, r.undetermined, r.c_prime.mean, r.c_prime.min, r.c_prime.max,
                        r.weil_muls.mean, r.f4_muls.mean, r.roots_muls.mean, r.sign_group_ops.mean,
                        r.solving_degree.mean, r.solving_degree.min, r.solving_degree.max,
                        r.max_rows.mean, r.max_cols.mean, r.f4_ms.mean, r.c_prime_decomposable.mean,
                        r.solving_degree_decomposable.mean, r.mitm3_group_ops, r.wall_ms / 1e3
                    );
                    rows.push(Row::Cprime(r));
                }
                "dlp" => {
                    let r = run_jv4_dlp(p, seed, rho_runs, check_every);
                    println!(
                        "p={:>5} seed={} n=2^{:.1} |F|={:>4} | residuals {:>8} relations {:>5} rate {:.5} (1/6p = {:.5}) res/floor {:.3} | checked {} mismatches {} undetermined {} | C' {:.0} | LA {:.3e} mod-n over {} attempts, core {} φ {:.3} /N² {:.2} | S {:.1} = walk {:.1}% oracle {:.1}% LA {:.1}% setup {:.1}% | rho S {:.3} ± {:.3} ({} runs, all correct: {}) | S/rho {:.1} (relations {:.1}, r {:.3}; formula {:.1}) | solved {} correct {} | {:.0} s",
                        r.p, r.seed, r.bits, r.base, r.residuals, r.relations, r.decomposition_rate,
                        r.expected_rate, r.residual_ratio, r.cross_checked, r.mismatches, r.undetermined,
                        r.c_prime, r.la_ops as f64, r.la_attempts, r.unknowns, r.phi, r.wiedemann_constant,
                        r.s, 100.0 * r.walk_muls as f64 / r.total_muls as f64,
                        100.0 * r.oracle_muls as f64 / r.total_muls as f64,
                        100.0 * r.la_muls as f64 / r.total_muls as f64,
                        100.0 * (r.setup_muls + r.precompute_muls + r.verify_muls) as f64 / r.total_muls as f64,
                        r.rho_s_mean, r.rho_s_sd, r.rho.len(), r.rho.iter().all(|x| x.correct),
                        r.s_over_rho, r.relation_over_rho, r.r, r.predicted_s_over_rho,
                        r.solved, r.correct, r.wall_ms / 1e3
                    );
                    rows.push(Row::Dlp(r));
                }
                other => panic!("unknown experiment {other}"),
            }
            if let Some(path) = &json {
                fs::write(path, serde_json::to_string_pretty(&rows).unwrap()).expect("write json");
            }
        }
    }
}
