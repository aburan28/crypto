//! Relation yield of greedy coverage-driven factor bases at equal column
//! count, scored by the production census.
//!
//! The expected-trial formula charges one unknown per signed Frobenius
//! orbit ("column"): `(U + 1 + extra) / coverage`.  So two bases are
//! compared at equal `U`, and the comparison isolates how well the orbits
//! were chosen:
//!
//!   * `greedy` — `build_greedy_column_factor_base`, every candidate orbit
//!     scored at each step; deterministic, one base per `U`;
//!   * `greedy_pool64` — the same rule scoring a seeded sample of 64
//!     orbits per step, one base per seed;
//!   * `random` — `build_random_column_factor_base`, `U` orbits drawn
//!     uniformly from the same orbit list, one base per seed (the control).
//!
//! Every base is scored by `witness_list` and `census_from_witnesses` on
//! the exhaustive target set `⟨G⟩ ∖ {O}`.  For each greedy base the
//! builder's scalar-model coverage is also compared with the census.
//!
//! The shipped sampler (`build_subgroup_orbit_factor_base_with_cost`) is
//! recorded as a reference only, at the point budgets the sampler is
//! specified in: its batches overshoot those budgets, so its column count
//! is reported beside each row rather than matched.
//!
//! Output: `rows.json` (one row per base) and `summary.json` (per curve
//! and `U`: means over seeds, and 95% bootstrap intervals on the coverage
//! ratio, resampling seeds).  Nothing is computed beyond these.
//!
//!     cargo run --release --example koblitz_greedy_factor_base_yield -- \
//!         <out-dir> <degrees,comma-separated> <seeds>

use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use crypto_lib::cryptanalysis::koblitz_factor_base_search::{
    census_from_witnesses, witness_list, TargetSet,
};
use crypto_lib::cryptanalysis::koblitz_greedy_base::{
    build_greedy_column_factor_base, build_random_column_factor_base,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_subgroup_orbit_factor_base_with_cost, FrobeniusFactorBase, KoblitzCurve, PairSumTable,
};
use num_traits::ToPrimitive;
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde::Serialize;

/// Candidate orbits scored per step by the pooled greedy arm.
const POOL: usize = 64;
/// Column counts swept, capped by the orbits available.
const SAMPLE_U: [usize; 15] = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48];
/// Reference budgets for the shipped sampler, as multiples of `√r`.
const PROD_GRID: [f64; 4] = [0.5, 1.0, 2.0, 4.0];
/// Bootstrap resamples per interval.
const RESAMPLES: usize = 4000;

#[derive(Clone, Serialize)]
struct Row {
    a: u8,
    n: u32,
    subgroup_order: u64,
    /// `greedy`, `greedy_pool64`, `random` or `production` (reference).
    arm: &'static str,
    seed: u64,
    /// Columns requested (greedy and random) or implied by the budget.
    columns: usize,
    points: usize,
    covered: usize,
    targets: usize,
    coverage: f64,
    /// `(columns + 1) · targets / covered`, the census's own formula.
    expected_trials: f64,
    build_ms: f64,
    census_ms: f64,
    /// Greedy arms: the builder's scalar-model coverage at this column count.
    model_covered: Option<usize>,
    /// Greedy arms: whether the model count equals the census.
    model_match: Option<bool>,
    /// Production reference: the point budget requested.
    requested_points: Option<usize>,
}

#[derive(Serialize)]
struct Summary {
    a: u8,
    n: u32,
    subgroup_order: u64,
    columns: usize,
    targets: usize,
    random_seeds: usize,
    random_mean_coverage: f64,
    random_sd_coverage: f64,
    greedy_coverage: f64,
    /// Greedy coverage over the mean random coverage, with a 95% bootstrap
    /// interval from resampling the random seeds.
    coverage_ratio: f64,
    coverage_ratio_ci95: [f64; 2],
    /// Trial formula at equal `U`: random mean trials over greedy trials.
    trial_ratio: f64,
    pool_mean_coverage: f64,
    /// Coverage ratio of the pooled greedy over the mean random coverage.
    pool_coverage_ratio: f64,
    pool_coverage_ratio_ci95: [f64; 2],
    greedy_build_ms: f64,
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

/// Ratio of two means with a percentile bootstrap over seeds of each arm.
/// The greedy side is one deterministic base when `den` has one entry, so
/// its interval reflects the random seeds alone.
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

/// Production census of one base on the shared exhaustive target set:
/// `(covered, expected trials, census wall time in ms)`.
fn census(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    targets: &TargetSet,
) -> Option<(usize, f64, f64)> {
    let start = Instant::now();
    let table = PairSumTable::build(kc, fb)?;
    let witnesses = witness_list(kc, fb, &table, targets, 2);
    let census_ms = start.elapsed().as_secs_f64() * 1000.0;
    let c = census_from_witnesses(
        &witnesses,
        targets,
        2,
        fb.signed_orbits.len(),
        fb.signed_orbits.len(),
        0,
        census_ms,
    );
    Some((c.covered, c.expected_trials, census_ms))
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
    fs::create_dir_all(&out_dir).expect("out-dir");

    let mut rows: Vec<Row> = Vec::new();
    let mut summaries: Vec<Summary> = Vec::new();
    let mut boot_rng = StdRng::seed_from_u64(0x4752_4545_4459_0001);

    for &n in &degrees {
        for a in [0u8, 1u8] {
            let Some(kc) = KoblitzCurve::new(a, n) else {
                eprintln!("skip K_{a}/2^{n}: no curve");
                continue;
            };
            let Some(r) = kc.subgroup_order.to_u64() else {
                continue;
            };
            let targets = TargetSet::new(&kc, 0, u64::MAX, 0);
            let t = targets.len();
            let units_available = build_greedy_column_factor_base(&kc, 1, 0, 0)
                .map(|(_, rep)| rep.units_available)
                .unwrap_or(0);
            if units_available == 0 {
                continue;
            }
            eprintln!("K_{a}/2^{n}: r = {r} ({t} targets), {units_available} orbits");

            // Reference: the shipped sampler, at its own point budgets.
            let sqrt_r = (r as f64).sqrt();
            for &c in &PROD_GRID {
                let b = ((c * sqrt_r).round() as usize).max(2);
                for seed in 1..=seeds {
                    let start = Instant::now();
                    let Ok((fb, _)) = build_subgroup_orbit_factor_base_with_cost(&kc, seed, b)
                    else {
                        continue;
                    };
                    let build_ms = start.elapsed().as_secs_f64() * 1000.0;
                    let Some((covered, trials, census_ms)) = census(&kc, &fb, &targets) else {
                        continue;
                    };
                    rows.push(Row {
                        a,
                        n,
                        subgroup_order: r,
                        arm: "production",
                        seed,
                        columns: fb.signed_orbits.len(),
                        points: fb.points.len(),
                        covered,
                        targets: t,
                        coverage: covered as f64 / t as f64,
                        expected_trials: trials,
                        build_ms,
                        census_ms,
                        model_covered: None,
                        model_match: None,
                        requested_points: Some(b),
                    });
                }
            }

            let ucap = units_available.min(*SAMPLE_U.last().expect("grid"));
            for &u in SAMPLE_U.iter().filter(|&&u| u <= ucap) {
                // Deterministic greedy at U columns.
                let start = Instant::now();
                let Ok((g_fb, g_rep)) = build_greedy_column_factor_base(&kc, u, 0, 0) else {
                    continue;
                };
                let g_build_ms = start.elapsed().as_secs_f64() * 1000.0;
                let Some((g_cov, g_trials, g_census_ms)) = census(&kc, &g_fb, &targets) else {
                    continue;
                };
                let g_row = Row {
                    a,
                    n,
                    subgroup_order: r,
                    arm: "greedy",
                    seed: 0,
                    columns: u,
                    points: g_fb.points.len(),
                    covered: g_cov,
                    targets: t,
                    coverage: g_cov as f64 / t as f64,
                    expected_trials: g_trials,
                    build_ms: g_build_ms,
                    census_ms: g_census_ms,
                    model_covered: Some(g_rep.covered_model),
                    model_match: Some(g_rep.covered_model == g_cov),
                    requested_points: None,
                };

                // Pooled greedy and uniform random, one base per seed.
                let mut pool_rows = Vec::new();
                let mut rand_rows = Vec::new();
                for seed in 1..=seeds {
                    let start = Instant::now();
                    if let Ok((fb, _)) = build_greedy_column_factor_base(&kc, u, POOL, seed) {
                        let build_ms = start.elapsed().as_secs_f64() * 1000.0;
                        if let Some((covered, trials, census_ms)) = census(&kc, &fb, &targets) {
                            pool_rows.push(Row {
                                a,
                                n,
                                subgroup_order: r,
                                arm: "greedy_pool64",
                                seed,
                                columns: u,
                                points: fb.points.len(),
                                covered,
                                targets: t,
                                coverage: covered as f64 / t as f64,
                                expected_trials: trials,
                                build_ms,
                                census_ms,
                                model_covered: None,
                                model_match: None,
                                requested_points: None,
                            });
                        }
                    }
                    let start = Instant::now();
                    if let Ok(fb) = build_random_column_factor_base(&kc, u, seed) {
                        let build_ms = start.elapsed().as_secs_f64() * 1000.0;
                        if let Some((covered, trials, census_ms)) = census(&kc, &fb, &targets) {
                            rand_rows.push(Row {
                                a,
                                n,
                                subgroup_order: r,
                                arm: "random",
                                seed,
                                columns: u,
                                points: fb.points.len(),
                                covered,
                                targets: t,
                                coverage: covered as f64 / t as f64,
                                expected_trials: trials,
                                build_ms,
                                census_ms,
                                model_covered: None,
                                model_match: None,
                                requested_points: None,
                            });
                        }
                    }
                }
                if rand_rows.is_empty() {
                    rows.push(g_row);
                    rows.extend(pool_rows);
                    continue;
                }

                let r_cov: Vec<f64> = rand_rows.iter().map(|x| x.coverage).collect();
                let r_trials: Vec<f64> = rand_rows.iter().map(|x| x.expected_trials).collect();
                let r_build: Vec<f64> = rand_rows.iter().map(|x| x.build_ms).collect();
                let p_cov: Vec<f64> = pool_rows.iter().map(|x| x.coverage).collect();
                let all_match = g_row.model_match == Some(true);
                let (ratio, ci) = bootstrap_ratio(&[g_row.coverage], &r_cov, &mut boot_rng);
                let (pool_ratio, pool_ci) = if p_cov.is_empty() {
                    (f64::NAN, [f64::NAN, f64::NAN])
                } else {
                    bootstrap_ratio(&p_cov, &r_cov, &mut boot_rng)
                };
                let s = Summary {
                    a,
                    n,
                    subgroup_order: r,
                    columns: u,
                    targets: t,
                    random_seeds: r_cov.len(),
                    random_mean_coverage: mean(&r_cov),
                    random_sd_coverage: sd(&r_cov),
                    greedy_coverage: g_row.coverage,
                    coverage_ratio: ratio,
                    coverage_ratio_ci95: ci,
                    trial_ratio: mean(&r_trials) / g_row.expected_trials,
                    pool_mean_coverage: if p_cov.is_empty() {
                        f64::NAN
                    } else {
                        mean(&p_cov)
                    },
                    pool_coverage_ratio: pool_ratio,
                    pool_coverage_ratio_ci95: pool_ci,
                    greedy_build_ms: g_row.build_ms,
                    random_build_ms_mean: mean(&r_build),
                    all_models_match: all_match,
                };
                eprintln!(
                    "K_{a}/2^{n} U={u:>2}: greedy {}/{t} ({:.3}), random {:.1}±{:.1} ({:.3}) \
                     ratio {:.3} [{:.3}, {:.3}]  pool {:.3} ratio {:.3}  model {}",
                    g_row.covered,
                    g_row.coverage,
                    s.random_mean_coverage * t as f64,
                    s.random_sd_coverage * t as f64,
                    s.random_mean_coverage,
                    s.coverage_ratio,
                    s.coverage_ratio_ci95[0],
                    s.coverage_ratio_ci95[1],
                    s.pool_mean_coverage,
                    s.pool_coverage_ratio,
                    if all_match { "agrees" } else { "DISAGREES" },
                );
                summaries.push(s);
                rows.push(g_row);
                rows.extend(pool_rows);
                rows.extend(rand_rows);
            }
        }
    }

    fs::write(
        out_dir.join("rows.json"),
        serde_json::to_string_pretty(&rows).expect("rows json"),
    )
    .expect("write rows");
    fs::write(
        out_dir.join("summary.json"),
        serde_json::to_string_pretty(&summaries).expect("summary json"),
    )
    .expect("write summary");
    eprintln!("wrote {} rows, {} summaries", rows.len(), summaries.len());
}
