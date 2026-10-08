//! Split-stage j=0 orbit IC ladder rung with a paired rho control.
//!
//! ```bash
//! cargo run --release --example j0_stage_ladder -- 16 20261005
//! cargo run --release --example j0_stage_ladder -- 20 20261005
//! ```
//!
//! Public synthetic known-answer only.  One deterministic seeded target
//! per invocation; the scalar is a validation-only sidecar and is never
//! supplied to either solver.  IC stage timers are split
//! (factor_base -> relations -> LA -> verify) per the prime
//! `end_to_end_dlp` ledger next-target; rho runs on the identical
//! public point as the agreement control (this is not a vs_rho claim —
//! at these sizes IC wall exceeds rho wall).

use crypto_lib::cryptanalysis::ec_index_calculus::pollard_rho_ecdlp;
use crypto_lib::cryptanalysis::ec_index_calculus_j0::j0_index_calculus_dlp_staged;
use crypto_lib::cryptanalysis::research_bench::bench_curves_j0;
use num_traits::Zero;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::time::Instant;

fn main() {
    let arguments: Vec<String> = std::env::args().collect();
    let bits: u32 = arguments
        .get(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(20);
    let target_seed: u64 = arguments
        .get(2)
        .and_then(|a| a.parse().ok())
        .unwrap_or(20261005);

    let Some((_, curve)) = bench_curves_j0().into_iter().find(|(b, _)| *b == bits) else {
        eprintln!("no j=0 bench curve at {bits} bits");
        std::process::exit(2);
    };

    let g = curve.generator();
    let a_fe = curve.a_fe();

    // Deterministic target: validation-only sidecar, never solver input.
    let mut rng = StdRng::seed_from_u64(target_seed);
    let limbs = (curve.n.bits() / 64 + 1) as usize;
    let mut x_truth = num_bigint::BigUint::zero();
    for i in 0..limbs {
        let word: u64 = rng.gen();
        x_truth = x_truth | (num_bigint::BigUint::from(word) << (64 * i));
    }
    x_truth = (x_truth % (&curve.n - 1u32)) + 1u32;
    if x_truth.is_zero() || x_truth >= curve.n {
        eprintln!("x_truth out of range");
        std::process::exit(2);
    }
    let target_generation_started = Instant::now();
    let q = g.scalar_mul(&x_truth, &a_fe);
    let target_generation_ms = target_generation_started.elapsed().as_secs_f64() * 1e3;
    let (q_x, q_y) = match &q {
        crypto_lib::ecc::point::Point::Affine { x, y } => (x.value.clone(), y.value.clone()),
        _ => std::process::exit(2),
    };

    // Frozen parameter policy: identical to the 16-bit rung policy in
    // j0_known_answer_bits (orbits capped at 40, extra 8), with the
    // 20-bit+ rung's per-relation trial cap raised to 1M and a bounded
    // retry budget so a heavy relation cannot abort the rung.
    let target_orbits = (1usize << (bits / 2 + 2)).min(40);
    let extra = 8;
    let max_trials = (1usize << (bits + 4)).min(1_000_000);
    let max_relation_attempts = 64;
    let rho_steps = (1usize << ((bits / 2) + 8)).min(5_000_000);

    println!(
        "curve={} p_bits≈{bits} n={} orbits={target_orbits} trials={max_trials} attempts={max_relation_attempts} rho_steps={rho_steps} seed={target_seed}",
        curve.name, curve.n
    );

    let t0 = Instant::now();
    let staged = j0_index_calculus_dlp_staged(
        &curve,
        &g,
        &q,
        target_orbits,
        extra,
        max_trials,
        max_relation_attempts,
    );
    let ic_ms = t0.elapsed().as_secs_f64() * 1e3;
    let (ic_recovered, stage) = match staged {
        Some((x, report)) => (Some(x), Some(report)),
        None => (None, None),
    };
    let ic_ok = ic_recovered
        .as_ref()
        .map(|r| g.scalar_mul(r, &a_fe) == q)
        .unwrap_or(false);
    let ic_match = ic_recovered
        .as_ref()
        .map(|r| r == &x_truth)
        .unwrap_or(false);

    let t1 = Instant::now();
    let rho = pollard_rho_ecdlp(&curve, &g, &q, rho_steps);
    let rho_ms = t1.elapsed().as_secs_f64() * 1e3;
    let rho_ok = rho
        .as_ref()
        .map(|r| g.scalar_mul(r, &a_fe) == q)
        .unwrap_or(false);

    let agree = match (&ic_recovered, &rho) {
        (Some(a), Some(b)) => a == b,
        _ => false,
    };

    let stage_json = stage
        .map(|s| serde_json::json!(s))
        .unwrap_or(serde_json::json!(null));

    println!(
        "{}",
        serde_json::json!({
            "bits": bits,
            "curve": curve.name,
            "target_seed": target_seed,
            "x_truth": x_truth.to_string(),
            "public_target_q": [q_x.to_string(), q_y.to_string()],
            "target_generation_ms_excluded": target_generation_ms,
            "target_orbits": target_orbits,
            "extra_relations": extra,
            "max_trials_per_relation": max_trials,
            "max_relation_attempts": max_relation_attempts,
            "ic_ms": ic_ms,
            "ic_ok": ic_ok,
            "ic_matches_truth": ic_match,
            "ic_recovered": ic_recovered.as_ref().map(|v| v.to_string()),
            "ic_stages": stage_json,
            "rho_ms": rho_ms,
            "rho_ok": rho_ok,
            "rho_recovered": rho.as_ref().map(|v| v.to_string()),
            "rho_steps_cap": rho_steps,
            "ic_agrees_rho": agree,
            "claim_boundary": "synthetic known-answer j=0 IC stage ladder with paired rho agreement control; not a vs_rho claim (IC wall exceeds rho wall at these sizes)",
        })
    );

    if !ic_ok || !rho_ok || !agree {
        std::process::exit(1);
    }
}
