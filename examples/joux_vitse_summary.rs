//! Summarise the frozen Joux–Vitse `n = 5` study
//! (`research/joux_vitse_n1_20261015/PROTOCOL.md`, "Analysis plan"):
//! per-rung trial cost, decomposition rate with Wilson intervals, the
//! counting floor, cost per relation, log–log fits, the DLP arms against
//! rho, and the crossover model (labelled extrapolation).
//!
//! ```text
//! cargo run --release --example joux_vitse_summary -- research/joux_vitse_n1_20261015/data \
//!     > research/joux_vitse_n1_20261015/summary.json
//! ```

use std::collections::BTreeMap;
use std::env;
use std::fs;

use serde_json::{json, Value};

fn read_jsonl(path: &std::path::Path) -> Vec<Value> {
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

/// Wilson 95 % interval for `k` successes in `n`.
fn wilson(k: f64, n: f64) -> (f64, f64) {
    if n == 0.0 {
        return (f64::NAN, f64::NAN);
    }
    let z = 1.959964f64;
    let ph = k / n;
    let d = 1.0 + z * z / n;
    let c = (ph + z * z / (2.0 * n)) / d;
    let h = z * (ph * (1.0 - ph) / n + z * z / (4.0 * n * n)).sqrt() / d;
    ((c - h).max(0.0), c + h)
}

/// Least squares `y = a + b x`: `(a, b, se_b)`.
fn ols(xs: &[f64], ys: &[f64]) -> (f64, f64, f64) {
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let sxx: f64 = xs.iter().map(|x| (x - mx).powi(2)).sum();
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let b = sxy / sxx;
    let a = my - b * mx;
    let rss: f64 = xs
        .iter()
        .zip(ys)
        .map(|(x, y)| (y - a - b * x).powi(2))
        .sum();
    let se = if n > 2.0 {
        (rss / (n - 2.0) / sxx).sqrt()
    } else {
        f64::NAN
    };
    (a, b, se)
}

fn choose(n: f64, k: u32) -> f64 {
    (0..k).fold(1.0, |acc, i| acc * (n - i as f64) / (i as f64 + 1.0))
}

fn main() {
    let dir = std::path::PathBuf::from(env::args().nth(1).expect("data directory"));
    // ── ladder ────────────────────────────────────────────────────────────
    let mut rungs: BTreeMap<u64, Vec<Value>> = BTreeMap::new();
    let mut entries: Vec<_> = fs::read_dir(dir.join("ladder"))
        .map(|d| d.filter_map(Result::ok).map(|e| e.path()).collect())
        .unwrap_or_default();
    entries.sort();
    for path in entries {
        for row in read_jsonl(&path) {
            let p = row["p"].as_u64().unwrap();
            rungs.entry(p).or_default().push(row);
        }
    }
    let mut rung_out = Vec::new();
    let mut checks_failed = 0u64;
    for (&p, rows) in &rungs {
        let f4: Vec<&Value> = rows.iter().filter(|r| r["kind"] == "f4_trial").collect();
        let sums: Vec<&Value> = rows.iter().filter(|r| r["kind"] == "summary").collect();
        let setups: Vec<&Value> = rows.iter().filter(|r| r["kind"] == "setup").collect();
        if sums.len() < 2 {
            continue; // rung incomplete: reported by its absence
        }
        let cost: Vec<f64> = f4
            .iter()
            .map(|r| {
                let x = &r["report"];
                u(&x["muls_system"]) + u(&x["muls_f4"]) + u(&x["muls_solve"])
            })
            .collect();
        let n = cost.len() as f64;
        let mean = cost.iter().sum::<f64>() / n;
        let sd = (cost.iter().map(|c| (c - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
        let min = cost.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = cost.iter().cloned().fold(0.0, f64::max);
        let replayed = f4
            .iter()
            .filter(|r| r["report"]["replayed"] == true)
            .count();
        let mismatch = f4
            .iter()
            .filter(|r| r["report"]["replay_mismatch"] == true)
            .count();
        let consistent = f4
            .iter()
            .filter(|r| r["report"]["consistent"] == true)
            .count();
        let disagree = f4.iter().filter(|r| r["agree"] != true).count();
        let unverified = f4
            .iter()
            .filter(|r| r["relations_verified"] != true)
            .count();
        checks_failed += (disagree + unverified) as u64;
        let f4_hits = f4
            .iter()
            .filter(|r| u(&r["mitm_decompositions"]) > 0.0)
            .count() as f64;
        let f4_decomps: f64 = f4.iter().map(|r| u(&r["mitm_decompositions"])).sum();
        let f4_rel_trials = f4
            .iter()
            .filter(|r| {
                r["report"]["relations"]
                    .as_array()
                    .is_some_and(|a| !a.is_empty())
            })
            .count();
        let solutions: f64 = f4
            .iter()
            .map(|r| u(&r["report"]["split"]["solutions"]))
            .sum();
        let not_split: f64 = f4
            .iter()
            .map(|r| u(&r["report"]["split"]["not_split"]))
            .sum();
        let mitm_trials: f64 = sums.iter().map(|s| u(&s["mitm_trials"])).sum();
        let mitm_hits: f64 = sums.iter().map(|s| u(&s["mitm_hit_trials"])).sum();
        let mitm_decomps: f64 = sums.iter().map(|s| u(&s["mitm_decompositions"])).sum();
        let trials = mitm_trials + n;
        let hits = mitm_hits + f4_hits;
        let rate = hits / trials;
        let (lo, hi) = wilson(hits, trials);
        let base: f64 =
            setups.iter().map(|s| u(&s["base_points"])).sum::<f64>() / setups.len() as f64;
        let e = (p as f64).powi(5); // #E = p⁵ + O(p^{5/2})
        let expected_decomps = 16.0 * choose(base, 4) / e;
        let rho_max_distinct = expected_decomps; // ≥ P(at least one)
        let rho_max_multiset = 16.0 * choose(base + 3.0, 4) / e;
        let decomps_per_trial = (mitm_decomps + f4_decomps) / trials;
        let c_rel = mean / rate;
        let sym: f64 = setups
            .iter()
            .map(|s| u(&s["sym_precompute_muls"]))
            .sum::<f64>()
            / setups.len() as f64;
        let trace: f64 = setups
            .iter()
            .map(|s| u(&s["trace_record_muls"]))
            .sum::<f64>()
            / setups.len() as f64;
        rung_out.push(json!({
            "p": p, "base_points_mean": base, "f4_trials": n,
            "c_trial_mean": mean, "c_trial_sd": sd, "c_trial_min": min, "c_trial_max": max,
            "replayed": replayed, "replay_mismatch": mismatch, "consistent": consistent,
            "solutions": solutions, "not_split": not_split,
            "f4_trials_with_relation": f4_rel_trials,
            "mitm_disagree": disagree, "relations_unverified": unverified,
            "rate_trials": trials, "rate_hits": hits, "rate": rate,
            "rate_wilson95": [lo, hi],
            "rate_times_24p": rate * 24.0 * p as f64,
            "decompositions_per_trial": decomps_per_trial,
            "expected_decompositions_per_trial_16C": expected_decomps,
            "rate_over_floor_distinct": rate / rho_max_distinct,
            "rate_floor_multiset": rho_max_multiset,
            "trials_per_relation": 1.0 / rate,
            "trials_per_relation_wilson95": [1.0 / hi, 1.0 / lo],
            "c_rel": c_rel, "c_rel_wilson95": [mean / hi, mean / lo],
            "setup_sym_precompute_muls": sym, "setup_trace_record_muls": trace,
        }));
    }
    let lp: Vec<f64> = rung_out.iter().map(|r| u(&r["p"]).ln()).collect();
    let lc: Vec<f64> = rung_out
        .iter()
        .map(|r| u(&r["c_trial_mean"]).ln())
        .collect();
    let fit_trial = if lp.len() >= 2 {
        let (a, b, se) = ols(&lp, &lc);
        json!({"intercept_ln": a, "slope": b, "slope_se": se, "rungs": lp.len()})
    } else {
        Value::Null
    };
    let good: Vec<&Value> = rung_out
        .iter()
        .filter(|r| u(&r["rate_hits"]) >= 5.0)
        .collect();
    let gp: Vec<f64> = good.iter().map(|r| u(&r["p"]).ln()).collect();
    let gc: Vec<f64> = good.iter().map(|r| u(&r["c_rel"]).ln()).collect();
    let gt: Vec<f64> = good
        .iter()
        .map(|r| u(&r["trials_per_relation"]).ln())
        .collect();
    let (fit_rel, fit_tpr) = if gp.len() >= 2 {
        let (a, b, se) = ols(&gp, &gc);
        let (a2, b2, se2) = ols(&gp, &gt);
        (
            json!({"intercept_ln": a, "slope": b, "slope_se": se, "rungs": gp.len()}),
            json!({"intercept_ln": a2, "slope": b2, "slope_se": se2, "rungs": gp.len()}),
        )
    } else {
        (Value::Null, Value::Null)
    };
    // ── DLP arms ──────────────────────────────────────────────────────────
    let mut dlp_out = Vec::new();
    let mut dlp_entries: Vec<_> = fs::read_dir(dir.join("dlp"))
        .map(|d| d.filter_map(Result::ok).map(|e| e.path()).collect())
        .unwrap_or_default();
    dlp_entries.sort();
    let mut c_add = Vec::new();
    for path in dlp_entries {
        let rows = read_jsonl(&path);
        let get = |k: &str| rows.iter().find(|r| r["kind"] == k).cloned();
        let (Some(inst), Some(rho)) = (get("instance"), get("rho")) else {
            continue;
        };
        c_add.push(u(&rho["fp_muls_per_op"]));
        let ic = get("ic");
        let relations = rows.iter().filter(|r| r["kind"] == "relation").count();
        dlp_out.push(json!({
            "file": path.file_name().unwrap().to_string_lossy(),
            "p": inst["p"], "seed": inst["seed"], "order": inst["order"], "ell": inst["ell"],
            "cofactor": inst["cofactor"], "k_secret": inst["k_secret"],
            "rho_k": rho["k"], "rho_verified": rho["verified"],
            "rho_group_ops": rho["group_ops"], "rho_fp_muls": rho["fp_muls"],
            "rho_expected_ops_no_negation": (std::f64::consts::PI * u(&inst["ell"]) / 2.0).sqrt(),
            "ic_complete": ic.is_some(),
            "ic": ic.clone().unwrap_or(Value::Null),
            "relations_recorded": relations,
            "ratio_online_ic_over_rho": ic.as_ref().map(|i| i["ratio_online_ic_over_rho"].clone()),
            // against the expected walk alone (no multiplier set-up), at the
            // measured cost of one group operation
            "ratio_online_ic_over_expected_rho_walk": ic.as_ref().map(|i| {
                u(&i["online_muls"])
                    / ((std::f64::consts::PI * u(&inst["ell"]) / 2.0).sqrt()
                        * u(&rho["fp_muls_per_op"]))
            }),
        }));
    }
    // ── Model (extrapolation) ─────────────────────────────────────────────
    let c_add_mean = if c_add.is_empty() {
        f64::NAN
    } else {
        c_add.iter().sum::<f64>() / c_add.len() as f64
    };
    let c_trial_all =
        rung_out.iter().map(|r| u(&r["c_trial_mean"])).sum::<f64>() / rung_out.len().max(1) as f64;
    // trials per relation = κ q with κ from the fit at slope 1 (geometric mean)
    let kappa = if good.is_empty() {
        f64::NAN
    } else {
        (good
            .iter()
            .map(|r| (u(&r["trials_per_relation"]) / u(&r["p"])).ln())
            .sum::<f64>()
            / good.len() as f64)
            .exp()
    };
    let model = |q: f64| -> Value {
        let base = q / 2.0;
        let relation = (base + 1.0) * kappa * q * c_trial_all;
        let la = 4.0 * base * base * 25.0; // 4 entries/row, Wiedemann, 25 F_p-mult-equiv per mod-ℓ mult
        let rho_strong = (std::f64::consts::PI * q.powi(5) / 4.0).sqrt() * c_add_mean;
        json!({"log2_q": q.log2(), "log2_group": 5.0 * q.log2(),
               "ic_relation_muls": relation, "ic_la_muls": la,
               "rho_negation_muls": rho_strong,
               "ic_over_rho": (relation + la) / rho_strong})
    };
    let qstar =
        (kappa * c_trial_all / (2.0 * (std::f64::consts::PI / 4.0).sqrt() * c_add_mean)).powi(2);
    let ladder_model: Vec<Value> = [10.0, 16.0, 20.0, 30.0, 40.0, 50.0, 60.0, 64.0, 70.0]
        .iter()
        .map(|&b: &f64| model(2f64.powf(b)))
        .collect();
    let out = json!({
        "protocol": "research/joux_vitse_n1_20261015/PROTOCOL.md",
        "checks_failed": checks_failed,
        "rungs": rung_out,
        "fit_c_trial_vs_p": fit_trial,
        "fit_c_rel_vs_p": fit_rel,
        "fit_trials_per_relation_vs_p": fit_tpr,
        "dlp": dlp_out,
        "model_extrapolation": {
            "label": "extrapolation from the measured C_trial (mean over rungs), trials per relation = κ q (κ from the rungs with ≥ 5 hits), c_add from the DLP rho arms, rho with the negation map on a prime-order group of p⁵ points; |F| = q/2",
            "c_trial": c_trial_all, "kappa": kappa, "c_add": c_add_mean,
            "q_star": qstar, "log2_q_star": qstar.log2(), "log2_group_at_q_star": 5.0 * qstar.log2(),
            "rows": ladder_model,
        },
        "thesis_reference": {
            "source": "Vitse, PhD thesis (UVSQ, 2011), §3.3.1 Table 3.1, §7.3.3, Remarque 7.3.6",
            "trial_muls_f4": 1e10, "trial_muls_f4remake_approx": 3.3e9,
            "trial_muls_stated_7_3_3": 3e9,
            "crossover_base_field_bits_odd_p": 60,
        },
        "d4_n4_reference": {
            "source": "PR #1863 (open, stacked on #1792): C4 = 5.5–8.5e7 F_p mults, p = 29…4049; IC 1.39e6× rho at p = 29 (ecbench)",
            "c4": 6.9e7,
        },
    });
    println!("{}", serde_json::to_string_pretty(&out).unwrap());
}
