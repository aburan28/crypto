//! Semaev-S₄ 3-decomposition IC ladder on the a=−3 deployed-shape curves —
//! the prime-regime `decomposition` stage lever.
//!
//! ```bash
//! cargo run --release --example a3_s4_ladder -- 16 20261005
//! cargo run --release --example a3_s4_ladder -- 20 20261005
//! cargo run --release --example a3_s4_ladder -- 24 20261005
//! ```
//!
//! Relations decompose `R = aG + bQ` as a signed sum of **three**
//! factor-base points via the S₄ quartic and Cantor–Zassenhaus root
//! finding (`find_roots_fp_fast`).  This is the first 3-decomposition
//! relation family past toy sizes on the P-256 / GOST CryptoPro-B curve
//! shape (prime order, a = p − 3, generic j).  Known-answer only; the
//! scalar is a validation-only sidecar, never a solver input.
//!
//! Oracle class: direct quartic root finding (no Gröbner / SAT system,
//! so FFD is inapplicable — recorded per the decomposition schema).

use crypto_lib::cryptanalysis::ec_index_calculus::{
    ec_index_calculus_dlp_s4_staged, pollard_rho_ecdlp,
};
use crypto_lib::cryptanalysis::research_bench::bench_curves_a_minus_3;
use crypto_lib::cryptanalysis::residual_walk::is_prime_u64;
use crypto_lib::ecc::point::Point;
use num_bigint::BigUint;
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
        .unwrap_or(20261006);
    let class_filter: Option<String> = arguments.get(3).cloned();

    let Some((_, curve)) = bench_curves_a_minus_3()
        .into_iter()
        .find(|(b, c)| *b == bits && class_filter.as_deref().map_or(true, |f| c.name.contains(f)))
    else {
        eprintln!("no a=-3 bench curve at {bits} bits matching {class_filter:?}");
        std::process::exit(2);
    };

    let a_fe = curve.a_fe();
    let g = curve.generator();
    let p_u64 = curve.p.to_u64_digits().first().copied().unwrap_or(0);
    let n_u64 = curve.n.to_u64_digits().first().copied().unwrap_or(0);
    if p_u64 == 0 || n_u64 == 0 || !is_prime_u64(p_u64) || !is_prime_u64(n_u64) {
        eprintln!("fixture primality check failed");
        std::process::exit(2);
    }
    if g.scalar_mul(&curve.n, &a_fe) != Point::Infinity {
        eprintln!("[n]G is not infinity");
        std::process::exit(2);
    }
    if curve.a != (&curve.p - 3u32) {
        eprintln!("fixture is not an a = p - 3 curve");
        std::process::exit(2);
    }

    // Deterministic target: validation-only sidecar, never solver input.
    let mut rng = StdRng::seed_from_u64(target_seed);
    let limbs = (curve.n.bits() / 64 + 1) as usize;
    let mut x_truth = BigUint::zero();
    for i in 0..limbs {
        let word: u64 = rng.gen();
        x_truth = x_truth | (BigUint::from(word) << (64 * i));
    }
    x_truth = (x_truth % (&curve.n - 1u32)) + 1u32;
    let target_generation_started = Instant::now();
    let q = g.scalar_mul(&x_truth, &a_fe);
    let target_generation_ms = target_generation_started.elapsed().as_secs_f64() * 1e3;
    let (q_x, q_y) = match &q {
        Point::Affine { x, y } => (x.value.clone(), y.value.clone()),
        _ => std::process::exit(2),
    };

    let fb_size = 120;
    let extra = 8;
    let max_trials = 1_000_000;
    let max_relation_attempts = 64;
    let rho_steps = (1usize << ((bits / 2) + 8)).min(5_000_000);

    println!(
        "curve={} p_bits≈{bits} n={} fb={fb_size} s4-trials={max_trials} attempts={max_relation_attempts} seed={target_seed}",
        curve.name, curve.n
    );

    let t0 = Instant::now();
    let staged =
        ec_index_calculus_dlp_s4_staged(&curve, &g, &q, fb_size, extra, max_trials, max_relation_attempts);
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
            "curve_class": if curve.name.contains("cryptoproclass") { "cryptopro_b_shape" } else { "p256_shape" },
            "decomposition": {
                "m_summands": 3,
                "system_degree": 4,
                "oracle_class": "quartic_root_find_cantor_zassenhaus",
                "ffd": {"status": "inapplicable", "reason": "direct univariate quartic root finding; no Boolean/Gröbner system solved"},
                "unknowns": fb_size,
            },
            "target_seed": target_seed,
            "x_truth": x_truth.to_string(),
            "public_target_q": [q_x.to_string(), q_y.to_string()],
            "target_generation_ms_excluded": target_generation_ms,
            "fb_size": fb_size,
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
            "ic_agrees_rho": agree,
            "claim_boundary": "synthetic known-answer S4 3-decomposition IC ladder on a=-3 (P-256/CryptoPro-B shape) curves; not a vs_rho claim; no deployed-curve claim",
        })
    );

    if !ic_ok || !rho_ok || !agree {
        std::process::exit(1);
    }
}
