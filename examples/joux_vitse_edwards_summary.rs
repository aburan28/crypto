//! Summarise `PROTOCOL_EDWARDS.md`'s runs (`data/edwards/`): the replay
//! factor per residual, the complete DLPs against rho and against the
//! counting floor, and the crossover model (labelled extrapolation).
//!
//! ```text
//! cargo run --release --example joux_vitse_edwards_summary -- \
//!     research/joux_vitse_n1_20261015/data/edwards > research/joux_vitse_n1_20261015/summary_edwards.json
//! ```

use std::env;
use std::fs;

use serde_json::{json, Value};

fn rows(path: &std::path::Path) -> Vec<Value> {
    fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter_map(|l| serde_json::from_str(l).ok())
        .collect()
}

fn u(v: &Value) -> f64 {
    v.as_f64()
        .or_else(|| v.as_str().and_then(|s| s.parse().ok()))
        .unwrap_or(f64::NAN)
}

/// Exact (Garwood) 95 % interval for a Poisson mean from `k` events, by
/// bisection on the Poisson tails.
fn poisson95(k: u64) -> (f64, f64) {
    let cdf = |lam: f64, k: u64| -> f64 {
        // P(X ≤ k)
        let mut term = (-lam).exp();
        let mut sum = term;
        for i in 1..=k {
            term *= lam / i as f64;
            sum += term;
        }
        sum
    };
    let solve = |f: &dyn Fn(f64) -> f64| -> f64 {
        let (mut lo, mut hi) = (0.0f64, 1e4f64);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if f(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };
    let lower = if k == 0 {
        0.0
    } else {
        solve(&|lam| 0.025 - (1.0 - cdf(lam, k - 1)))
    };
    let upper = solve(&|lam| cdf(lam, k) - 0.025);
    (lower, upper)
}

fn choose(n: f64, k: u32) -> f64 {
    (0..k).fold(1.0, |acc, i| acc * (n - i as f64) / (i as f64 + 1.0))
}

fn main() {
    let dir = std::path::PathBuf::from(env::args().nth(1).expect("data/edwards directory"));
    let mut files: Vec<_> = fs::read_dir(&dir)
        .expect("read dir")
        .filter_map(Result::ok)
        .map(|e| e.path())
        .collect();
    files.sort();
    let (mut cost_out, mut dlp_out) = (Vec::new(), Vec::new());
    let (mut replay_all, mut c_add) = (Vec::new(), Vec::new());
    let mut checks_failed = 0u64;
    for path in &files {
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        let rs = rows(path);
        let setup = rs.iter().find(|r| r["kind"] == "setup").cloned();
        if name.starts_with("cost_") {
            let ts: Vec<&Value> = rs.iter().filter(|r| r["kind"] == "cost_trial").collect();
            if ts.is_empty() {
                continue;
            }
            let full: Vec<f64> = ts
                .iter()
                .map(|t| u(&t["full"]["solve_muls"]) + u(&t["full"]["weil_muls"]))
                .collect();
            let rep: Vec<f64> = ts
                .iter()
                .map(|t| {
                    u(&t["replay"]["replay_muls"])
                        + u(&t["replay"]["solve_muls"])
                        + u(&t["replay"]["weil_muls"])
                })
                .collect();
            let ratio: Vec<f64> = full.iter().zip(&rep).map(|(a, b)| a / b).collect();
            let bad = ts
                .iter()
                .filter(|t| {
                    t["same_decompositions"] != true
                        || t["oracle_agree"] != true
                        || t["verified"] != true
                })
                .count();
            checks_failed += bad as u64;
            let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
            let min = |v: &[f64]| v.iter().cloned().fold(f64::INFINITY, f64::min);
            let max = |v: &[f64]| v.iter().cloned().fold(0.0, f64::max);
            replay_all.extend(rep.iter().cloned());
            cost_out.push(json!({
                "file": name, "p": setup.as_ref().map(|s| s["p"].clone()),
                "seed": setup.as_ref().map(|s| s["seed"].clone()),
                "base_points": setup.as_ref().map(|s| s["base_points"].clone()),
                "residuals": ts.len(),
                "full_mean": mean(&full), "full_min": min(&full), "full_max": max(&full),
                "replay_mean": mean(&rep), "replay_min": min(&rep), "replay_max": max(&rep),
                "factor_mean": mean(&ratio), "factor_min": min(&ratio), "factor_max": max(&ratio),
                "replay_mismatches": ts.iter().filter(|t| t["replay"]["replay_mismatch"] == true).count(),
                "checks_failed": bad,
            }));
        } else if name.starts_with("dlp_") {
            let (Some(setup), Some(rho)) = (setup, rs.iter().find(|r| r["kind"] == "rho")) else {
                continue;
            };
            c_add.push(u(&rho["fp_muls_per_op"]));
            // the last `ic` row: a resumed run may carry an earlier one from a
            // segment that stopped at `--max-trials`
            let ic = rs.iter().rev().find(|r| r["kind"] == "ic").cloned();
            let n = u(&setup["n"]);
            let base = u(&setup["base_points"]);
            // the census expectation (RESEARCH_K5 §3): 8 subgroup points per
            // class-quadruple, so a random residual decomposes at rate
            // 8·C(|F|, 4)/n on average — an expectation, not a bound
            let floor_trials = n / (8.0 * choose(base, 4));
            let relations = rs.iter().filter(|r| r["kind"] == "relation").count() as f64;
            let row = match &ic {
                Some(ic) => {
                    checks_failed += (u(&ic["oracle_disagree"]) + u(&ic["relations_unverified"]))
                        as u64
                        + u64::from(ic["verified"] != true)
                        + u64::from(rho["verified"] != true);
                    let trials = u(&ic["trials"]);
                    let rels = u(&ic["relations"]);
                    let tpr = trials / rels;
                    // relations ~ Poisson(trials / expected): ratio of the
                    // measured to the census trials per relation, 95 %
                    let (lo, hi) = poisson95(rels as u64);
                    let ratio_ci = [rels / hi, rels / lo.max(1e-300)];
                    json!({
                        "file": name, "p": setup["p"], "seed": setup["seed"], "n": setup["n"],
                        "log2_n": setup["log2_n"], "base_points": base,
                        "complete": true, "verified": ic["verified"], "rho_verified": rho["verified"],
                        "trials": trials, "relations": ic["relations"], "unknowns": ic["unknowns"],
                        "trials_per_relation": tpr, "census_trials_per_relation": floor_trials,
                        "trials_per_relation_over_census": tpr / floor_trials,
                        "trials_per_relation_over_census_95": [ratio_ci[0] * tpr / floor_trials, ratio_ci[1] * tpr / floor_trials],
                        "c_trial": (u(&ic["online_muls"]) - u(&ic["linear_algebra_mulmods"])) / trials,
                        "online_muls": ic["online_muls"], "rho_group_ops": rho["group_ops"],
                        "rho_fp_muls": rho["fp_muls"], "rho_expected_walk_ops": rho["expected_walk_ops"],
                        "ratio_online_ic_over_rho": ic["ratio_online_ic_over_rho"],
                        "ratio_online_ic_over_expected_rho_walk": ic["ratio_online_ic_over_expected_rho_walk"],
                        "agg": ic["agg"], "oracle_disagree": ic["oracle_disagree"],
                    })
                }
                None => json!({
                    "file": name, "p": setup["p"], "seed": setup["seed"], "n": setup["n"],
                    "base_points": base, "complete": false, "relations_recorded": relations,
                    "census_trials_per_relation": floor_trials,
                    "rho_group_ops": rho["group_ops"], "rho_verified": rho["verified"],
                }),
            };
            dlp_out.push(row);
        }
    }
    // Model (extrapolation): IC ≈ (p/4) · ρ · 192 p · C, rho ≈ √(πn/4) c_add,
    // n = p⁵/4, with C the mean replayed trial at the cost-mode sizes and
    // ρ the measured trials per relation over the counting floor.
    let c_trial = replay_all.iter().sum::<f64>() / replay_all.len().max(1) as f64;
    let c_add_m = c_add.iter().sum::<f64>() / c_add.len().max(1) as f64;
    let rhos: Vec<f64> = dlp_out
        .iter()
        .filter(|d| d["complete"] == true)
        .map(|d| u(&d["trials_per_relation_over_census"]))
        .collect();
    let rho_floor = if rhos.is_empty() {
        f64::NAN
    } else {
        (rhos.iter().map(|x| x.ln()).sum::<f64>() / rhos.len() as f64).exp()
    };
    let k_rho = (std::f64::consts::PI / 16.0).sqrt(); // √(π n /4) = k_rho p^{5/2}
    let p_star_at = |rho: f64| (48.0 * rho * c_trial / (k_rho * c_add_m)).powi(2);
    let p_star = p_star_at(1.0);
    let p_star_measured = p_star_at(rho_floor);
    let model_rows: Vec<Value> = [8.0, 10.0, 16.0, 20.0, 30.0, 40.0, 45.0, 50.0]
        .iter()
        .map(|&b: &f64| {
            let p = 2f64.powf(b);
            let ic = 48.0 * p * p * c_trial;
            let rho = k_rho * p.powf(2.5) * c_add_m;
            json!({"log2_p": b, "log2_n": 5.0 * b - 2.0, "ic_over_rho": ic / rho})
        })
        .collect();
    let out = json!({
        "protocol": "research/joux_vitse_n1_20261015/PROTOCOL_EDWARDS.md",
        "checks_failed": checks_failed,
        "cost": cost_out,
        "dlp": dlp_out,
        "model_extrapolation": {
            "label": "extrapolation: IC = (p/4)·192p·C_trial (relations |F|+1 ≈ p/4 at the census rate 1/(192p), verified to three digits in RESEARCH_K5 §3), rho = √(πn/4)·c_add with the negation map, n = p⁵/4; C_trial the mean replayed trial over the cost runs. Sensitivity: the same at the measured trials-per-relation ratio",
            "c_trial_replay": c_trial, "c_add": c_add_m,
            "p_star": p_star, "log2_p_star": p_star.log2(),
            "log2_n_star": 5.0 * p_star.log2() - 2.0,
            "measured_ratio_to_census": rho_floor,
            "log2_n_star_at_measured_ratio": 5.0 * p_star_measured.log2() - 2.0,
            "rows": model_rows,
        },
        "committed_reference": {
            "source": "research/notes/index-calculus/RESEARCH_K5_TORSION_JOUX_VITSE.md §3",
            "c_second": 7.3e7, "log2_n_star": 233,
        },
    });
    println!("{}", serde_json::to_string_pretty(&out).unwrap());
}
