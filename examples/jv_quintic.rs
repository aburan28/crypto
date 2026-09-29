//! **Joux–Vitse four-point decompositions at `k = 5`: `C″`.**
//!
//! Companion to `research/notes/index-calculus/RESEARCH_RHO_PARITY_PROGRAMME.md` §6.
//!
//! ```text
//! cargo run --release --example jv_quintic -- --sizes 271,521,761,1031 --seeds 1 --residuals 4 --constructed 4 --max-degree 24 --budget 1800 --json experiments/27_jv_quintic_csecond.json
//! ```

use std::env;
use std::fs;

use crypto_lib::cryptanalysis::jv_quintic::{run_jv5_csecond, CsecondReport};

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut sizes: Vec<u64> = vec![271, 521, 761, 1031];
    let mut seeds = 1u64;
    let mut first_seed = 1u64;
    let mut residuals = 4usize;
    let mut constructed = 4usize;
    let mut max_degree = 24u32;
    let mut budget = 1800f64;
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
            "--max-degree" => {
                i += 1;
                max_degree = args[i].parse().expect("--max-degree");
            }
            "--budget" => {
                i += 1;
                budget = args[i].parse().expect("--budget");
            }
            "--json" => {
                i += 1;
                json = Some(args[i].clone());
            }
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    let mut rows: Vec<CsecondReport> = Vec::new();
    for &p in &sizes {
        for seed in first_seed..first_seed + seeds {
            let r = run_jv5_csecond(p, seed, residuals, constructed, max_degree, budget);
            println!(
                "p={:>5} seed={} n=2^{:.1} |F|={:>4} c_add {:.0} | precompute {:.2e} muls {:.1} s | random {} decomposable {} (1/24p = {:.5}) | planted {}/{} | mismatches {} unverified {} undetermined {} timed out {} | C'' {:.3e} [{:.2e}, {:.2e}] = weil {:.0} + f4 {:.3e} + roots {:.0} + signs {:.1} ops | F4 degree reached {:.1} rows {:.0} cols {:.0} {:.1} s | decomposable C'' {:.3e} degree {:.1} cols {:.0} {:.1} s | mitm4 {:.0} ops | {:.0} s",
                r.p, r.seed, r.bits, r.base, r.fp_muls_per_add, r.precompute_muls as f64,
                r.precompute_ms / 1e3, r.random_residuals, r.random_decomposable, r.expected_rate,
                r.planted_found, r.constructed_residuals, r.mismatches, r.unverified, r.undetermined,
                r.timed_out, r.c_second.mean, r.c_second.min as f64, r.c_second.max as f64,
                r.weil_muls.mean, r.f4_muls.mean, r.roots_muls.mean, r.sign_group_ops.mean,
                r.degree_reached.mean, r.max_rows.mean, r.max_cols.mean, r.f4_ms.mean / 1e3,
                r.c_second_decomposable.mean, r.degree_reached_decomposable.mean,
                r.max_cols_decomposable.mean, r.f4_ms_decomposable.mean / 1e3, r.mitm4_group_ops,
                r.wall_ms / 1e3
            );
            rows.push(r);
            if let Some(path) = &json {
                fs::write(path, serde_json::to_string_pretty(&rows).unwrap()).expect("write json");
            }
        }
    }
}
