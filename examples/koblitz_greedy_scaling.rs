//! Does the greedy coverage selector keep its advantage as the subgroup
//! grows past what the dense selector can enumerate?
//!
//! For each Koblitz curve, at each column count `U`, two bases of `U`
//! signed Frobenius orbits are built by the sampled selector
//! (`koblitz_greedy_sampled`), which needs no enumeration of the subgroup:
//!
//!   * `greedy` — each step scores 64 freshly drawn orbits by the number of
//!     new sum abscissae and adds the best; one base per seed;
//!   * `random` — `U` orbits drawn uniformly by the same sampler (the
//!     control; same granularity, same points, no coverage information).
//!
//! Coverage of every base is then recomputed from the base alone by
//! `pair_sum_coverage` (all sums of two points, hashed), which is exact and
//! does not use the selector's own bookkeeping; for each greedy base the
//! selector's count is compared with it.  On the subgroups the dense
//! selector can handle, the same quantity is what the production census
//! reports (`koblitz_greedy_sampled` tests).
//!
//! Output: `rows.json` (one row per base) and `summary.json` (per curve and
//! `U`, means over seeds and a 95% bootstrap interval on the ratio of mean
//! coverage, resampling seeds).
//!
//!     cargo run --release --example koblitz_greedy_scaling -- \
//!         <out-dir> <degrees,comma-separated> <seeds> <max-points> [pool=64] [columns,comma-separated]

use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use crypto_lib::cryptanalysis::koblitz_greedy_base::Budget;
use crypto_lib::cryptanalysis::koblitz_greedy_sampled::{
    pair_sum_coverage, sampled_random, sampled_select,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_traits::ToPrimitive;
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde::Serialize;

/// Candidate orbits scored per greedy step, unless given on the command line.
const DEFAULT_POOL: usize = 64;
/// Column counts swept, capped by the point budget.
const GRID: [usize; 8] = [2, 4, 8, 16, 32, 48, 64, 96];
const RESAMPLES: usize = 4000;

#[derive(Clone, Serialize)]
struct Row {
    a: u8,
    n: u32,
    subgroup_bits: f64,
    arm: &'static str,
    seed: u64,
    pool: usize,
    columns: usize,
    points: usize,
    covered: u64,
    targets: u64,
    coverage: f64,
    build_ms: f64,
    /// Greedy only: the selector's own count, and whether it equals the
    /// independent pair-sum count.
    model_covered: Option<u64>,
    model_match: Option<bool>,
}

#[derive(Serialize)]
struct Summary {
    a: u8,
    n: u32,
    subgroup_bits: f64,
    columns: usize,
    pool: usize,
    points: usize,
    seeds: usize,
    greedy_mean_coverage: f64,
    random_mean_coverage: f64,
    random_sd_coverage: f64,
    /// Greedy mean coverage over random mean coverage.
    ratio: f64,
    ratio_ci95: [f64; 2],
    greedy_build_ms_mean: f64,
    random_build_ms_mean: f64,
    all_models_match: bool,
}

fn mean(xs: &[f64]) -> f64 {
    xs.iter().sum::<f64>() / xs.len() as f64
}

fn sd(xs: &[f64]) -> f64 {
    if xs.len() < 2 {
        return 0.0;
    }
    let m = mean(xs);
    (xs.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (xs.len() - 1) as f64).sqrt()
}

fn bootstrap_ratio(num: &[f64], den: &[f64], rng: &mut StdRng) -> (f64, [f64; 2]) {
    let point = mean(num) / mean(den);
    let mut draws: Vec<f64> = (0..RESAMPLES)
        .map(|_| {
            let mn = (0..num.len())
                .map(|_| num[rng.gen_range(0..num.len())])
                .sum::<f64>()
                / num.len() as f64;
            let md = (0..den.len())
                .map(|_| den[rng.gen_range(0..den.len())])
                .sum::<f64>()
                / den.len() as f64;
            mn / md
        })
        .collect();
    draws.sort_by(|a, b| a.total_cmp(b));
    let lo = draws[(0.025 * RESAMPLES as f64) as usize];
    let hi = draws[((0.975 * RESAMPLES as f64) as usize).min(RESAMPLES - 1)];
    (point, [lo, hi])
}

fn main() {
    let mut args = std::env::args().skip(1);
    let out_dir = PathBuf::from(args.next().expect("out-dir"));
    let degrees: Vec<u32> = args
        .next()
        .expect("degrees")
        .split(',')
        .map(|s| s.parse().expect("degree"))
        .collect();
    let seeds: u64 = args.next().and_then(|s| s.parse().ok()).expect("seeds");
    let max_points: usize = args
        .next()
        .and_then(|s| s.parse().ok())
        .expect("max-points");
    let pool: usize = args
        .next()
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_POOL);
    let grid: Vec<usize> = match args.next() {
        Some(list) => list
            .split(',')
            .map(|s| s.parse().expect("columns"))
            .collect(),
        None => GRID.to_vec(),
    };
    fs::create_dir_all(&out_dir).expect("out-dir");

    let mut rows: Vec<Row> = Vec::new();
    let mut summaries: Vec<Summary> = Vec::new();
    let mut boot = StdRng::seed_from_u64(0x5343_414c_494e_4701);

    for &n in &degrees {
        for a in [0u8, 1u8] {
            let Some(kc) = KoblitzCurve::new(a, n) else {
                eprintln!("skip K_{a}/2^{n}: no curve");
                continue;
            };
            let bits = kc.subgroup_order.to_f64().map_or(f64::NAN, f64::log2);
            eprintln!("K_{a}/2^{n}: r = 2^{bits:.1}");
            for &u in grid.iter().filter(|&&u| u * 2 * n as usize <= max_points) {
                let mut g_cov = Vec::new();
                let mut r_cov = Vec::new();
                let mut g_ms = Vec::new();
                let mut r_ms = Vec::new();
                let mut all_match = true;
                let mut points = 0usize;
                for seed in 1..=seeds {
                    let start = Instant::now();
                    let Ok((fb, rep)) = sampled_select(&kc, Budget::Columns(u), pool, seed) else {
                        continue;
                    };
                    let build_ms = start.elapsed().as_secs_f64() * 1000.0;
                    let Ok((covered, targets)) = pair_sum_coverage(&kc, &fb) else {
                        continue;
                    };
                    let matches = rep.covered_model as u64 == covered;
                    all_match &= matches;
                    points = fb.points.len();
                    g_cov.push(covered as f64 / targets as f64);
                    g_ms.push(build_ms);
                    rows.push(Row {
                        a,
                        n,
                        subgroup_bits: bits,
                        arm: "greedy",
                        seed,
                        pool,
                        columns: u,
                        points,
                        covered,
                        targets,
                        coverage: covered as f64 / targets as f64,
                        build_ms,
                        model_covered: Some(rep.covered_model as u64),
                        model_match: Some(matches),
                    });
                    let start = Instant::now();
                    let Ok(fb) = sampled_random(&kc, u, seed) else {
                        continue;
                    };
                    let build_ms = start.elapsed().as_secs_f64() * 1000.0;
                    let Ok((covered, targets)) = pair_sum_coverage(&kc, &fb) else {
                        continue;
                    };
                    r_cov.push(covered as f64 / targets as f64);
                    r_ms.push(build_ms);
                    rows.push(Row {
                        a,
                        n,
                        subgroup_bits: bits,
                        arm: "random",
                        seed,
                        pool: 0,
                        columns: u,
                        points: fb.points.len(),
                        covered,
                        targets,
                        coverage: covered as f64 / targets as f64,
                        build_ms,
                        model_covered: None,
                        model_match: None,
                    });
                }
                if g_cov.is_empty() || r_cov.is_empty() {
                    continue;
                }
                let (ratio, ci) = bootstrap_ratio(&g_cov, &r_cov, &mut boot);
                eprintln!(
                    "K_{a}/2^{n} U={u:>3} ({points} pts): greedy {:.3e} random {:.3e}±{:.1e} \
                     ratio {ratio:.4} [{:.4}, {:.4}]  greedy {:.0} ms  {}",
                    mean(&g_cov),
                    mean(&r_cov),
                    sd(&r_cov),
                    ci[0],
                    ci[1],
                    mean(&g_ms),
                    if all_match {
                        "model agrees"
                    } else {
                        "MODEL DISAGREES"
                    },
                );
                summaries.push(Summary {
                    a,
                    n,
                    subgroup_bits: bits,
                    columns: u,
                    pool,
                    points,
                    seeds: g_cov.len(),
                    greedy_mean_coverage: mean(&g_cov),
                    random_mean_coverage: mean(&r_cov),
                    random_sd_coverage: sd(&r_cov),
                    ratio,
                    ratio_ci95: ci,
                    greedy_build_ms_mean: mean(&g_ms),
                    random_build_ms_mean: mean(&r_ms),
                    all_models_match: all_match,
                });
            }
        }
    }
    fs::write(
        out_dir.join("rows.json"),
        serde_json::to_string_pretty(&rows).expect("rows"),
    )
    .expect("write rows");
    fs::write(
        out_dir.join("summary.json"),
        serde_json::to_string_pretty(&summaries).expect("summary"),
    )
    .expect("write summary");
    eprintln!("wrote {} rows, {} summaries", rows.len(), summaries.len());
}
