//! Known-answer j=0 orbit IC vs Pollard ρ on the research-bench curve ladder.
//!
//! ```bash
//! cargo run --release --example j0_known_answer_bits -- 14
//! cargo run --release --example j0_known_answer_bits -- 16
//! ```
//!
//! Public synthetic only. Not a ledger promotion by itself.

use crypto_lib::cryptanalysis::ec_index_calculus::pollard_rho_ecdlp;
use crypto_lib::cryptanalysis::ec_index_calculus_j0::j0_index_calculus_dlp;
use crypto_lib::cryptanalysis::research_bench::bench_curves_j0;
use crypto_lib::utils::random::random_scalar;
use num_traits::Zero;
use std::time::Instant;

fn main() {
    let bits: u32 = std::env::args()
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(16);

    let Some((_, curve)) = bench_curves_j0().into_iter().find(|(b, _)| *b == bits) else {
        eprintln!("no j=0 bench curve at {bits} bits");
        std::process::exit(2);
    };

    let g = curve.generator();
    let a_fe = curve.a_fe();
    let x_truth = random_scalar(&curve.n);
    if x_truth.is_zero() {
        eprintln!("x_truth = 0");
        std::process::exit(2);
    }
    let q = g.scalar_mul(&x_truth, &a_fe);

    let target_orbits = (1usize << (bits / 2 + 2)).min(40);
    let extra = 8;
    let max_trials = (1usize << (bits + 4)).min(200_000);
    let rho_steps = (1usize << ((bits / 2) + 8)).min(5_000_000);

    println!("curve={} p_bits≈{bits} n={} orbits={target_orbits} trials={max_trials} rho_steps={rho_steps}",
        curve.name, curve.n);

    let t0 = Instant::now();
    let ic = j0_index_calculus_dlp(&curve, &g, &q, target_orbits, extra, max_trials);
    let ic_ms = t0.elapsed().as_secs_f64() * 1e3;
    let ic_ok = ic.as_ref().map(|r| g.scalar_mul(r, &a_fe) == q).unwrap_or(false);
    let ic_match = ic.as_ref().map(|r| r == &x_truth).unwrap_or(false);

    let t1 = Instant::now();
    let rho = pollard_rho_ecdlp(&curve, &g, &q, rho_steps);
    let rho_ms = t1.elapsed().as_secs_f64() * 1e3;
    let rho_ok = rho.as_ref().map(|r| g.scalar_mul(r, &a_fe) == q).unwrap_or(false);

    let agree = match (&ic, &rho) {
        (Some(a), Some(b)) => a == b,
        _ => false,
    };

    println!(
        "{}",
        serde_json::json!({
            "bits": bits,
            "curve": curve.name,
            "x_truth": x_truth.to_string(),
            "ic_ms": ic_ms,
            "ic_ok": ic_ok,
            "ic_matches_truth": ic_match,
            "ic_recovered": ic.as_ref().map(|v| v.to_string()),
            "rho_ms": rho_ms,
            "rho_ok": rho_ok,
            "rho_recovered": rho.as_ref().map(|v| v.to_string()),
            "ic_agrees_rho": agree,
            "target_orbits": target_orbits,
            "max_trials": max_trials,
            "rho_steps": rho_steps,
            "claim_boundary": "synthetic known-answer j=0 IC vs rho",
        })
    );

    if !ic_ok || !rho_ok || !agree {
        std::process::exit(1);
    }
}
