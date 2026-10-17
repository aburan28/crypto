//! Summarise `research/gaudry_quartic_d4_char2_20261015` (PROTOCOL.md
//! Parts A and C) from its recorded rows and sessions: per configuration
//! the mean `C₄` with a 95 % bootstrap interval over targets (10,000
//! resamples, seed 1; with replay over the replayed solves), phase means,
//! closing and exactness counts; paired ratios; the trend of `S_IC / S_rho`
//! against `log₂ n` over the ecbench sessions and its extrapolation to 1;
//! and §11.16's model with the measured binary inputs.  Writes
//! `raw/summary.json`.
//!
//! ```sh
//! cargo run --release --example gaudry_quartic_bin_summary -- research/gaudry_quartic_d4_char2_20261015
//! ```

use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::Path;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

const PHASES: [&str; 9] = [
    "weil_muls",
    "triangularise_muls",
    "echelon_muls",
    "nf_muls",
    "charpoly_muls",
    "roots_muls",
    "eigenvector_muls",
    "verify_muls",
    "split_muls",
];

/// rho with the negation map, `S` in theory (`ecbench table`'s column).
const RHO_THEORY_S: f64 = 0.886;

fn lines(path: &Path) -> Vec<Value> {
    fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).expect("a JSON row"))
        .collect()
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

fn bootstrap(n: usize, stat: impl Fn(&[usize]) -> f64) -> [f64; 2] {
    let mut rng = StdRng::seed_from_u64(1);
    let mut s: Vec<f64> = (0..10_000)
        .map(|_| {
            let idx: Vec<usize> = (0..n).map(|_| rng.gen_range(0..n)).collect();
            stat(&idx)
        })
        .collect();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    [s[249], s[9_749]]
}

fn total(r: &Value) -> f64 {
    r["solve"]["total_muls"].as_f64().unwrap()
}

fn steady(rows: &[Value]) -> Vec<&Value> {
    rows.iter()
        .filter(|r| r["solve"]["replay"] != "recorded")
        .collect()
}

fn config(raw: &Path, name: &str) -> Option<Value> {
    let rows = lines(&raw.join(format!("{name}.jsonl")));
    if rows.is_empty() {
        return None;
    }
    let log = fs::read_to_string(raw.join(format!("{name}.log"))).unwrap_or_default();
    let unsolved = log.lines().filter(|l| l.contains("did not close")).count();
    let unsolved_paired = log
        .lines()
        .filter(|l| l.contains("with a pair"))
        .filter(|l| {
            // "mitm N decompositions, K with a pair": all of them paired.
            let n: Vec<&str> = l.split_whitespace().collect();
            n.len() > 5 && n[3] == n[5]
        })
        .count();
    let st = steady(&rows);
    let v: Vec<f64> = st.iter().map(|r| total(r)).collect();
    let ci = bootstrap(v.len(), |idx| {
        mean(&idx.iter().map(|&i| v[i]).collect::<Vec<_>>())
    });
    let first = rows
        .iter()
        .find(|r| r["solve"]["replay"] == "recorded")
        .map(total);
    let phases: BTreeMap<&str, f64> = PHASES
        .iter()
        .map(|&k| {
            (
                k,
                mean(
                    &st.iter()
                        .map(|r| r["solve"][k].as_f64().unwrap())
                        .collect::<Vec<_>>(),
                ),
            )
        })
        .collect();
    Some(json!({
        "name": name,
        "m": rows[0]["m"],
        "planted": rows[0]["planted"],
        "form": rows[0]["form"],
        "macaulay_degree": rows[0]["macaulay_degree"],
        "rows": rows.last().unwrap()["solve"]["rows"],
        "cols": rows[0]["solve"]["cols"],
        "dim": rows[0]["solve"]["dim"],
        "chains": rows[0]["solve"]["chains"],
        "closed": rows.len(),
        "unsolved": unsolved,
        "unsolved_all_decompositions_paired": unsolved_paired,
        "decompositions_agree": rows.iter().filter(|r| r["decompositions_agree"].as_bool() == Some(true)).count(),
        "mitm_decompositions": rows.iter().map(|r| r["mitm_decompositions"].as_u64().unwrap()).sum::<u64>(),
        "degenerate_eigenspaces": rows.iter().map(|r| r["solve"]["degenerate_eigenspaces"].as_u64().unwrap()).sum::<u64>(),
        "replays": rows.iter().filter(|r| r["solve"]["replay"] == "replayed").count(),
        "fallbacks": rows.iter().filter(|r| r["solve"]["replay"] == "fell back").count(),
        "first_solve_muls": first,
        "n_steady": v.len(),
        "mean_muls": mean(&v),
        "ci95_muls": ci,
        "phase_means": phases,
        "mean_wall_s": mean(&st.iter().map(|r| r["solve"]["total_ms"].as_f64().unwrap() / 1e3).collect::<Vec<_>>()),
    }))
}

fn paired(raw: &Path, a: &str, b: &str) -> Option<Value> {
    let ra = lines(&raw.join(format!("{a}.jsonl")));
    let rb = lines(&raw.join(format!("{b}.jsonl")));
    let mut pairs = Vec::new();
    for x in steady(&ra) {
        let t = &x["target"];
        if let Some(y) = steady(&rb).into_iter().find(|y| &y["target"] == t) {
            pairs.push((
                total(x),
                total(y),
                x["solver_decompositions"] == y["solver_decompositions"],
            ));
        }
    }
    if pairs.is_empty() {
        return None;
    }
    let ratio = |idx: &[usize]| {
        idx.iter().map(|&i| pairs[i].0).sum::<f64>() / idx.iter().map(|&i| pairs[i].1).sum::<f64>()
    };
    let all: Vec<usize> = (0..pairs.len()).collect();
    Some(json!({
        "a": a, "b": b, "pairs": pairs.len(),
        "same_decompositions": pairs.iter().filter(|p| p.2).count(),
        "ratio": ratio(&all), "ci95": bootstrap(pairs.len(), ratio),
    }))
}

/// The sessions' IC runs: per run `(m, log₂ n, h, S_IC, trials,
/// relations, ops, per_add)` and the rho runs' `S`.
fn session_points(dir: &Path) -> (Vec<Value>, Vec<Value>) {
    let mut ic = Vec::new();
    let mut rho = Vec::new();
    for s in ["char2-forms-m6", "char2-ladder", "char2-m14"] {
        let sd = dir.join("sessions").join(s);
        let plan: Value = match fs::read_to_string(sd.join("plan.json")) {
            Ok(t) => serde_json::from_str(&t).unwrap(),
            Err(_) => continue,
        };
        let facts: BTreeMap<String, Value> = plan["workloads"]
            .as_array()
            .unwrap()
            .iter()
            .map(|w| {
                (
                    w["workload_id"].as_str().unwrap().to_string(),
                    w["curve"].clone(),
                )
            })
            .collect();
        for r in lines(&sd.join("records.jsonl")) {
            let wid = r["workload"]["workload_id"].as_str().unwrap().to_string();
            let c = &facts[&wid];
            let order = c["group_order"]
                .as_f64()
                .unwrap_or_else(|| c["group_order"].as_str().unwrap().parse().unwrap());
            let h = c["cofactor"]
                .as_f64()
                .unwrap_or_else(|| c["cofactor"].as_str().unwrap().parse().unwrap());
            let m = r["workload"]["curve_spec"]["m"].as_u64().unwrap();
            let s_val = r["cost"]["s"].as_f64().unwrap();
            let status = r["outcome"]["status"].as_str().unwrap().to_string();
            let arm = r["arm"].as_str().unwrap().to_string();
            if arm.starts_with("rho") {
                rho.push(json!({"session": s, "arm": arm, "m": m, "workload": wid, "s": s_val, "status": status}));
                continue;
            }
            let native = &r["phases"][2]["native"];
            ic.push(json!({
                "session": s, "arm": arm, "m": m, "workload": wid,
                "log2_n": order.log2(), "h": h, "status": status,
                "s": s_val, "ratio_to_rho_theory": s_val / RHO_THEORY_S,
                "trials": native["trials"], "relations": native["relations"],
                "solver_ops": r["solver"]["ops"], "extra": r["solver"]["extra"],
            }));
        }
    }
    (ic, rho)
}

fn fit(xs: &[f64], ys: &[f64]) -> (f64, f64) {
    let (mx, my) = (mean(xs), mean(ys));
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    let b = sxy / sxx;
    (my - b * mx, b)
}

fn main() {
    let dir = env::args().nth(1).expect("the study directory");
    let dir = Path::new(&dir);
    let raw = dir.join("raw");
    let mut names: Vec<String> = [
        "abl_plain_s",
        "abl_plain_sr",
        "abl_weighted",
        "abl_weighted_s",
        "abl_weighted_sr",
        "abl_frob",
        "abl_frob_s",
        "abl_frob_sr",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    for m in [5, 6, 7, 8, 10, 12] {
        names.push(format!("m{m}_planted"));
        names.push(format!("m{m}_uniform"));
    }
    names.push("m15_uniform_cost_only".into());
    let configs: Vec<Value> = names.iter().filter_map(|n| config(&raw, n)).collect();
    let pairs: Vec<Value> = [
        ("abl_plain_sr", "abl_weighted_sr"),
        ("abl_weighted_sr", "abl_frob_sr"),
        ("abl_plain_sr", "abl_frob_sr"),
        ("abl_frob", "abl_frob_s"),
        ("abl_frob_s", "abl_frob_sr"),
        ("abl_weighted", "abl_frob"),
    ]
    .iter()
    .filter_map(|(a, b)| paired(&raw, a, b))
    .collect();
    // Part C1: the measured trend over Frobenius runs with m ≥ 8.
    let (ic, rho) = session_points(dir);
    let pts: Vec<(f64, f64)> = ic
        .iter()
        .filter(|p| {
            p["arm"] == "ic-frobenius" && p["status"] == "verified" && p["m"].as_u64().unwrap() >= 8
        })
        .map(|p| {
            (
                p["log2_n"].as_f64().unwrap(),
                p["ratio_to_rho_theory"].as_f64().unwrap().log2(),
            )
        })
        .collect();
    let trend = if pts.len() >= 3 {
        let (xs, ys): (Vec<f64>, Vec<f64>) = pts.iter().cloned().unzip();
        let (a, b) = fit(&xs, &ys);
        // Leave one out: the spread of the extrapolation.
        let loo: Vec<(f64, f64)> = (0..pts.len())
            .map(|k| {
                let (x2, y2): (Vec<f64>, Vec<f64>) = pts
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != k)
                    .map(|(_, p)| *p)
                    .unzip();
                let (a2, b2) = fit(&x2, &y2);
                (b2, -a2 / b2)
            })
            .collect();
        json!({
            "points": pts, "intercept": a, "slope": b, "log2_n_star": -a / b,
            "leave_one_out": loo,
            "note": "log2(S_IC / 0.886) against log2 n, Frobenius arm, verified runs, m >= 8; extrapolated to 0",
        })
    } else {
        json!({"points": pts, "note": "fewer than three points"})
    };
    // Part C2: the model with measured inputs.
    let cfg = |n: &str| configs.iter().find(|c| c["name"] == n).cloned();
    let largest = ic
        .iter()
        .filter(|p| p["arm"] == "ic-frobenius" && p["status"] == "verified")
        .max_by_key(|p| p["m"].as_u64().unwrap());
    let model = match (cfg("m15_uniform_cost_only"), cfg("m12_uniform"), largest) {
        (Some(c15), Some(c12), Some(big)) => {
            let roots = |c: &Value| c["phase_means"]["roots_muls"].as_f64().unwrap();
            let base = c15["mean_muls"].as_f64().unwrap() - roots(&c15);
            let (r12, r15) = (roots(&c12), roots(&c15));
            let c_of = |m: f64| base + r15 + (r15 - r12) / 3.0 * (m - 15.0);
            let spr = big["trials"].as_f64().unwrap() / big["relations"].as_f64().unwrap();
            let c_add = big["extra"]["fp_muls_per_add"].as_f64().unwrap();
            let mut rows = Vec::new();
            for (label, rinf) in [("phi=0.73", 0.34 / 4.0), ("phi=1", 0.63 / 4.0)] {
                for h in [1.0f64, 2.0, 4.0] {
                    // √h (spr (q/4) C(m) / (1.3 q² c_add) + r∞) = 1, by bisection in m.
                    let ratio = |m: f64| {
                        let q = 2f64.powf(m);
                        h.sqrt() * (spr * (q / 4.0) * c_of(m) / (1.3 * q * q * c_add) + rinf)
                    };
                    let (mut lo, mut hi) = (4.0, 64.0);
                    for _ in 0..200 {
                        let mid = 0.5 * (lo + hi);
                        if ratio(mid) > 1.0 {
                            lo = mid;
                        } else {
                            hi = mid;
                        }
                    }
                    rows.push(
                        json!({"r_inf": label, "h": h, "m_star": lo, "log2_n_star": 4.0 * lo}),
                    );
                }
            }
            json!({
                "spr_trials_per_relation": spr, "spr_from_m": big["m"], "c_add_fq_muls": c_add,
                "c4_at_m15": c15["mean_muls"], "roots_at_m12": r12, "roots_at_m15": r15,
                "crossover": rows,
                "model": "S/rho = sqrt(h) (spr (q/4) C(m) / (1.3 q^2 c_add) + r_inf); C(m) = C4(15) with root finding linear in m; r_inf carried over from the odd study (an assumption)",
            })
        }
        _ => json!({"note": "inputs missing"}),
    };
    let out = json!({
        "schema": "gaudry_quartic_d4_char2.summary/v1",
        "configs": configs, "paired_ratios": pairs,
        "sessions": {"ic": ic, "rho": rho},
        "trend": trend, "model": model,
    });
    let path = raw.join("summary.json");
    fs::write(&path, serde_json::to_string_pretty(&out).unwrap() + "\n").unwrap();
    eprintln!("wrote {}", path.display());
}
