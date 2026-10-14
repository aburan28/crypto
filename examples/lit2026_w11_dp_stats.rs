//! lit2026 W11 (d): distinguished-point statistics of the repository's
//! signed-Frobenius (negation × Frobenius class) rho walk
//! ([`crypto_lib::cryptanalysis::koblitz_strong_rho`]) against the geometric
//! model, after the exploration-mode analysis of Panny et al. (ePrint
//! 2026/1823, HyperSolver).
//!
//! ```sh
//! cargo run --release --example lit2026_w11_dp_stats -- <a> <n> <dp_bits,...> <solves> [lanes]
//! ```
//!
//! For each `dp_bits` it runs `solves` independent solves (seeds `1..`, one
//! per rayon task, each walk single-threaded) of `K_a / F_{2^n}` with the
//! strong reference parameters (`lanes` default 32, fruitless walks
//! abandoned), records the length of every walk that ends on a
//! distinguished point and of every abandoned walk, and prints:
//!
//! - the DP hit rate per step against `θ = 2^{-dp_bits}`;
//! - the trail-length mean against the geometric `(1−θ)/θ`, and against the
//!   competing-risk geometric with the measured abandonment hazard
//!   `f = abandoned / steps`, mean `(1−θ)(1−f)/(1 − (1−θ)(1−f))`;
//! - the Kolmogorov–Smirnov distance of the DP trail lengths from both
//!   models' laws, conditioned on ending at a DP;
//! - a histogram in units of `1/θ`.
//!
//! All statistics are computed here, in Rust; the plotting script only draws
//! the printed values.

use std::time::Instant;

use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::koblitz_strong_rho::{
    RawPoint, StrongRho, StrongRhoCharges, StrongRhoParams,
};

/// KS distance of the sample from `F(k) = 1 − q^{k+1}` (a geometric law on
/// `0, 1, …` with continuation probability `q`).
fn ks_geometric(sorted: &[u32], q: f64) -> f64 {
    let n = sorted.len() as f64;
    let mut d: f64 = 0.0;
    let mut i = 0usize;
    while i < sorted.len() {
        let k = sorted[i];
        let mut j = i;
        while j < sorted.len() && sorted[j] == k {
            j += 1;
        }
        let model = 1.0 - q.powf(k as f64 + 1.0);
        let below = if k == 0 { 0.0 } else { 1.0 - q.powf(k as f64) };
        d = d
            .max((j as f64 / n - model).abs())
            .max((i as f64 / n - below).abs());
        i = j;
    }
    d
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    assert!(a.len() >= 4, "usage: a n dp_bits,... solves [lanes]");
    let curve_a: u8 = a[0].parse().unwrap();
    let n: u32 = a[1].parse().unwrap();
    let dps: Vec<u32> = a[2].split(',').map(|s| s.parse().unwrap()).collect();
    let solves: u64 = a[3].parse().unwrap();
    let lanes: usize = a.get(4).map(|s| s.parse().unwrap()).unwrap_or(32);
    let curve = KoblitzCurve::new(curve_a, n).expect("Koblitz curve");
    let rho = StrongRho::new(&curve);
    let modulus = rho.modulus();
    for &dp_bits in &dps {
        let t0 = Instant::now();
        let params = StrongRhoParams {
            lanes,
            dp_bits,
            ..StrongRhoParams::default()
        };
        let runs: Vec<_> = (1..=solves)
            .into_par_iter()
            .map(|seed| {
                let mut rng = StdRng::seed_from_u64(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15));
                let secret = rng.gen_range(1..modulus);
                let qb = curve.mul(curve.generator(), &BigUint::from(secret));
                let q = RawPoint::from_binary(&qb);
                let mut charges = StrongRhoCharges::default();
                let jumps = rho.jumps(seed ^ 0xD15C_0000, &mut charges);
                let out = rho
                    .solve(q, &jumps, &mut rng, &params, charges)
                    .expect("solved within the cap");
                assert_eq!(out.scalar, secret, "verified logarithm");
                out
            })
            .collect();
        let wall = t0.elapsed().as_secs_f64();
        let steps: u64 = runs.iter().map(|r| r.walk_steps).sum();
        let mut lens: Vec<u32> = runs
            .iter()
            .flat_map(|r| r.dp_trail_lengths.iter().copied())
            .collect();
        let abandoned: Vec<u32> = runs
            .iter()
            .flat_map(|r| r.abandoned_trail_lengths.iter().copied())
            .collect();
        let fruitless: u64 = runs.iter().map(|r| r.fruitless).sum();
        let capped: u64 = runs.iter().map(|r| r.capped).sum();
        let walks: u64 = runs.iter().map(|r| r.walks).sum();
        lens.sort_unstable();
        let theta = 0.5f64.powi(dp_bits as i32);
        let dpn = lens.len() as f64;
        let hit_rate = dpn / steps as f64;
        let mean = lens.iter().map(|&x| x as f64).sum::<f64>() / dpn;
        let var = lens.iter().map(|&x| (x as f64 - mean).powi(2)).sum::<f64>() / (dpn - 1.0);
        let f = abandoned.len() as f64 / steps as f64;
        let q_geo = 1.0 - theta;
        let q_cr = (1.0 - theta) * (1.0 - f);
        let mean_geo = q_geo / (1.0 - q_geo);
        let mean_cr = q_cr / (1.0 - q_cr);
        let ks_geo = ks_geometric(&lens, q_geo);
        let ks_cr = ks_geometric(&lens, q_cr);
        // histogram in units of 1/θ: 32 bins of width 0.25/θ up to 8/θ
        let w = 0.25 / theta;
        let mut hist = vec![0u64; 33];
        for &l in &lens {
            let b = ((l as f64) / w) as usize;
            hist[b.min(32)] += 1;
        }
        let model_hist: Vec<f64> = (0..33)
            .map(|b| {
                let lo = (b as f64 * w).ceil();
                let hi = if b == 32 {
                    f64::INFINITY
                } else {
                    ((b + 1) as f64 * w).ceil()
                };
                let p = q_geo.powf(lo) - if hi.is_finite() { q_geo.powf(hi) } else { 0.0 };
                p * dpn
            })
            .collect();
        let ideal: f64 = runs[0].ideal_steps;
        let ratios: Vec<f64> = runs.iter().map(|r| r.walk_steps as f64 / ideal).collect();
        let mean_ratio = ratios.iter().sum::<f64>() / ratios.len() as f64;
        let mean_ab = if abandoned.is_empty() {
            0.0
        } else {
            abandoned.iter().map(|&x| x as f64).sum::<f64>() / abandoned.len() as f64
        };
        println!(
            "{{\"curve\":\"K{curve_a}/F2^{n}\",\"a\":{curve_a},\"n\":{n},\"r\":{modulus},\"automorphisms\":{},\"lanes\":{lanes},\"dp_bits\":{dp_bits},\"theta\":{theta},\"solves\":{solves},\"steps\":{steps},\"walks\":{walks},\"dps\":{},\"hit_rate\":{hit_rate:.6e},\"hit_rate_over_theta\":{:.5},\"trail_mean\":{mean:.3},\"trail_var\":{var:.3},\"geo_mean\":{mean_geo:.3},\"trail_mean_over_geo\":{:.5},\"abandon_hazard\":{f:.4e},\"fruitless\":{fruitless},\"capped\":{capped},\"abandoned_mean_len\":{mean_ab:.3},\"cr_mean\":{mean_cr:.3},\"trail_mean_over_cr\":{:.5},\"ks_geo\":{ks_geo:.5},\"ks_cr\":{ks_cr:.5},\"ks_crit_95\":{:.5},\"ideal_steps\":{ideal:.1},\"steps_over_ideal_mean\":{mean_ratio:.4},\"steps_over_ideal\":{:?},\"hist_bin_width_over_inv_theta\":0.25,\"hist\":{:?},\"hist_model_geo\":{:?},\"wall_s\":{wall:.2}}}",
            runs[0].automorphisms,
            lens.len(),
            hit_rate / theta,
            mean / mean_geo,
            mean / mean_cr,
            1.358 / dpn.sqrt(),
            ratios.iter().map(|x| (x * 1e4).round() / 1e4).collect::<Vec<_>>(),
            hist,
            model_hist.iter().map(|x| (x * 10.0).round() / 10.0).collect::<Vec<_>>()
        );
    }
}
