//! **Gaudry-style index calculus on `E(F_{p³})` versus rho — measurement.**
//!
//! Companion to `crypto_lib::cryptanalysis::gaudry_cubic` and §11 of
//! `RESEARCH_RESIDUAL_WALKS.md`.  For each prime `p ≡ 1 (mod 3)` it
//! generates a prime-order curve over `F_{p³}`, the subspace factor
//! base `{P : x(P) ∈ F_p}`, collects triple decompositions with the
//! Weil-restricted `S₃` meet-in-the-middle oracle until `d` is
//! determined, and runs Pollard rho on the same group.
//!
//! ```text
//! cargo run --release --example gaudry_cubic_bench -- --p 271 --seed 1
//! cargo run --release --example gaudry_cubic_bench -- --protocol --json experiments/21_gaudry_cubic.json
//! ```

use std::env;
use std::fs;

use crypto_lib::cryptanalysis::gaudry_cubic::{
    generate_instance3, run_gaudry, run_rho3, GaudryReport, RhoReport3, SubspaceBase,
};
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde::Serialize;

#[derive(Serialize)]
struct Row {
    gaudry: GaudryReport,
    rho: RhoReport3,
}

fn one(p: u64, seed: u64, max_residuals: u64) -> Row {
    let inst = generate_instance3(p, seed);
    let mut rng = StdRng::seed_from_u64(seed);
    let base = SubspaceBase::build(&inst, &mut rng);
    let g = run_gaudry(&inst, &base, seed, max_residuals);
    let r = run_rho3(&inst, seed);
    println!(
        "p={:>5} n=2^{:5.1} B={:>4} | gaudry: residuals={:>6} decomp={:>5} rate={:.3} indep={:>4} pair_tests={:>9} Fp_muls/test={:>7.0} Fp/add={:>5.1} group_ops={:>9} total_ops={:>12.0} S={:>10.1} ok={:?} {:>7.0} ms | rho: steps={:>8} S={:>6.1} ok={:?}",
        g.p, g.bits, g.base, g.residuals, g.decompositions, g.decomposition_rate, g.relations_independent,
        g.pair_tests, g.fp_muls_per_pair_test, g.fp_muls_per_add, g.group_ops, g.total_ops, g.s, g.correct, g.wall_ms,
        r.steps, r.s, r.correct
    );
    Row { gaudry: g, rho: r }
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut p = 271u64;
    let mut seed = 1u64;
    let mut protocol = false;
    let mut json: Option<String> = None;
    let mut max_residuals = 1_000_000u64;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--p" => {
                i += 1;
                p = args[i].parse().expect("--p");
            }
            "--seed" => {
                i += 1;
                seed = args[i].parse().expect("--seed");
            }
            "--max-residuals" => {
                i += 1;
                max_residuals = args[i].parse().expect("--max-residuals");
            }
            "--json" => {
                i += 1;
                json = Some(args[i].clone());
            }
            "--protocol" => protocol = true,
            other => {
                eprintln!("unknown argument {other}");
                std::process::exit(2);
            }
        }
        i += 1;
    }
    let rows: Vec<Row> = if protocol {
        // Primes ≡ 1 mod 3 giving n ≈ 2^24, 2^27, 2^30, 2^33.
        let mut out = Vec::new();
        for &pp in &[271u64, 523, 1039, 2083] {
            for s in 1..=2u64 {
                out.push(one(pp, s, max_residuals));
            }
        }
        out
    } else {
        vec![one(p, seed, max_residuals)]
    };
    if let Some(path) = json {
        fs::write(&path, serde_json::to_string_pretty(&rows).unwrap()).expect("write json");
        println!("wrote {path}");
    }
}
