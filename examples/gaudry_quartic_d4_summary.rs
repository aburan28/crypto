//! Summarise the `C₄` ladder of `research/gaudry_quartic_d4_cost_20261014`
//! (PROTOCOL.md, Parts A and C) from the rows `gaudry_quartic_d4` wrote:
//! per configuration, the mean `C₄` with a 95 % bootstrap interval over
//! targets (10,000 resamples, seed 1; with replay, over the replayed
//! solves, the recorded first solve reported apart), the phase means, the
//! unsolved count from the log, and exactness against the enumeration; the
//! paired committed/new ratios; and §11.16's crossover with the measured
//! `C₄`.  Writes `summary.json` in the directory it reads.
//!
//! ```sh
//! cargo run --release --example gaudry_quartic_d4_summary -- research/gaudry_quartic_d4_cost_20261014/raw
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

fn rows(dir: &Path, name: &str) -> Vec<Value> {
    let path = dir.join(format!("{name}.jsonl"));
    fs::read_to_string(&path)
        .unwrap_or_default()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).expect("a row"))
        .collect()
}

fn unsolved(dir: &Path, name: &str) -> usize {
    fs::read_to_string(dir.join(format!("{name}.log")))
        .unwrap_or_default()
        .lines()
        .filter(|l| l.contains("did not close"))
        .count()
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

/// 95 % percentile interval of `stat` over resamples of `0..n`.
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

fn config(dir: &Path, name: &str) -> Option<Value> {
    let all = rows(dir, name);
    if all.is_empty() {
        return None;
    }
    // With replay the first solve of the file records the shape.
    let replay = fs::read_to_string(dir.join(format!("{name}.log")))
        .unwrap_or_default()
        .contains("shape replays:");
    let (first, steady): (Option<&Value>, Vec<&Value>) = if replay {
        (all.first(), all.iter().skip(1).collect())
    } else {
        (None, all.iter().collect())
    };
    let v: Vec<f64> = steady.iter().map(|r| total(r)).collect();
    let ci = bootstrap(v.len(), |idx| {
        mean(&idx.iter().map(|&i| v[i]).collect::<Vec<_>>())
    });
    let phases: BTreeMap<&str, f64> = PHASES
        .iter()
        .map(|&k| {
            (
                k,
                mean(
                    &steady
                        .iter()
                        .map(|r| r["solve"][k].as_f64().unwrap())
                        .collect::<Vec<_>>(),
                ),
            )
        })
        .collect();
    let agree = all
        .iter()
        .filter(|r| r["decompositions_agree"].as_bool() == Some(true))
        .count();
    let topped: Vec<u64> = steady
        .iter()
        .map(|r| r["solve"]["topped_up"].as_u64().unwrap_or(0))
        .collect();
    Some(json!({
        "name": name,
        "p": all[0]["p"],
        "planted": all[0]["planted"],
        "form": all[0]["form"],
        "macaulay_degree": all[0]["macaulay_degree"],
        "rows": all.last().unwrap()["solve"]["rows"],
        "cols": all[0]["solve"]["cols"],
        "solved": all.len(),
        "unsolved": unsolved(dir, name),
        "decompositions_agree": agree,
        "mitm_decompositions": all.iter().map(|r| r["mitm_decompositions"].as_u64().unwrap()).sum::<u64>(),
        "targets_with_decompositions": all.iter().filter(|r| r["mitm_decompositions"].as_u64().unwrap() > 0).count(),
        "degenerate_eigenspaces": all.iter().map(|r| r["solve"]["degenerate_eigenspaces"].as_u64().unwrap()).sum::<u64>(),
        "replay": replay,
        "first_solve_muls": first.map(total),
        "mean_muls": mean(&v),
        "ci95_muls": ci,
        "min_muls": v.iter().cloned().fold(f64::INFINITY, f64::min),
        "max_muls": v.iter().cloned().fold(0.0, f64::max),
        "n_steady": v.len(),
        "phase_means": phases,
        "topped_up_solves": topped.iter().filter(|&&t| t > 0).count(),
        "mean_wall_s": mean(&steady.iter().map(|r| r["solve"]["total_ms"].as_f64().unwrap() / 1e3).collect::<Vec<_>>()),
    }))
}

/// Ratio of mean `C₄` (`a` over `b`) on the targets both solved, with a
/// paired bootstrap; replayed rows of `b` only when `b` replays (the first,
/// recorded, solve is left out of both).
fn paired(dir: &Path, a: &str, b: &str) -> Option<Value> {
    let ra = rows(dir, a);
    let rb = rows(dir, b);
    let skip_first = |name: &str| {
        fs::read_to_string(dir.join(format!("{name}.log")))
            .unwrap_or_default()
            .contains("shape replays:")
    };
    let first_b = if skip_first(b) {
        rb.first().map(|r| r["target"].as_u64().unwrap())
    } else {
        None
    };
    let mut pairs = Vec::new();
    for x in &ra {
        let t = x["target"].as_u64().unwrap();
        if Some(t) == first_b {
            continue;
        }
        if let Some(y) = rb.iter().find(|y| y["target"].as_u64().unwrap() == t) {
            pairs.push((
                total(x),
                total(y),
                x["rational_solutions"] == y["rational_solutions"]
                    && x["solver_decompositions"] == y["solver_decompositions"],
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
        "same_solutions": pairs.iter().filter(|p| p.2).count(),
        "ratio": ratio(&all),
        "ci95": bootstrap(pairs.len(), ratio),
    }))
}

/// §11.16: `S/rho = √h (coef · C(p)/p + r∞)`, `coef = 181/(4·126)`; the
/// `n* = p⁴` where it is 1, by fixed point with `C(p) = base + roots(p)`,
/// `roots` linear in `log₂ p` through the two measured points.
fn crossover(base: f64, roots: Option<(f64, f64, f64, f64)>, h: f64, rinf: f64) -> f64 {
    let coef = 181.0 / (4.0 * 126.0);
    let c = |p: f64| match roots {
        Some((l1, r1, l2, r2)) => base + r1 + (r2 - r1) / (l2 - l1) * (p.log2() - l1),
        None => base,
    };
    let target = 1.0 / h.sqrt() - rinf;
    let mut p = 1e6;
    for _ in 0..200 {
        p = coef * c(p) / target;
    }
    4.0 * p.log2()
}

fn main() {
    let dir = env::args().nth(1).expect("the raw directory");
    let dir = Path::new(&dir);
    let names = [
        "abl_committed",
        "abl_s",
        "abl_sc",
        "abl_scr",
        "abl_w",
        "abl_ws",
        "abl_wc",
        "abl_wcs",
        "abl_wcsr",
        "p29_new_planted",
        "p29_new_uniform",
        "p29_committed_planted",
        "p29_committed_uniform",
        "p269_new_planted",
        "p269_new_uniform",
        "p269_committed_planted",
        "p1009_new_planted",
        "p1009_new_uniform",
        "p1009_committed_planted",
        "p4049_new_planted",
        "p4049_new_uniform",
        "p4049_committed_planted",
    ];
    let configs: Vec<Value> = names.iter().filter_map(|n| config(dir, n)).collect();
    let pairs: Vec<Value> = [
        ("abl_committed", "abl_wcsr"),
        ("abl_committed", "abl_s"),
        ("abl_scr", "abl_wcsr"),
        ("abl_sc", "abl_wcs"),
        ("p29_committed_planted", "p29_new_planted"),
        ("p29_committed_uniform", "p29_new_uniform"),
        ("p269_committed_planted", "p269_new_planted"),
        ("p1009_committed_planted", "p1009_new_planted"),
        ("p4049_committed_planted", "p4049_new_planted"),
    ]
    .iter()
    .filter_map(|(a, b)| paired(dir, a, b))
    .collect();
    // Crossover inputs: the largest p's replayed new solves (planted and
    // uniform pooled), its root-finding term at p = 1009 and 4049.
    let pooled = |p: u64| -> Option<(Vec<f64>, f64)> {
        let mut v = Vec::new();
        let mut roots = Vec::new();
        for n in [format!("p{p}_new_planted"), format!("p{p}_new_uniform")] {
            let rs = rows(dir, &n);
            for r in rs.iter().skip(1) {
                v.push(total(r));
                roots.push(r["solve"]["roots_muls"].as_f64().unwrap());
            }
        }
        (!v.is_empty()).then(|| (v, mean(&roots)))
    };
    let mut cross = Vec::new();
    if let (Some((v4049, r4049)), Some((_, r1009))) = (pooled(4049), pooled(1009)) {
        let c = mean(&v4049);
        let ci = bootstrap(v4049.len(), |idx| {
            mean(&idx.iter().map(|&i| v4049[i]).collect::<Vec<_>>())
        });
        let roots = Some((1009f64.log2(), r1009, 4049f64.log2(), r4049));
        for (label, rinf) in [("phi=0.73", 0.34 / 4.0), ("phi=1", 0.63 / 4.0)] {
            for h in [1.0, 2.0, 4.0] {
                let n_mid = crossover(c - r4049, roots, h, rinf);
                let n_lo = crossover(ci[0] - r4049, roots, h, rinf);
                let n_hi = crossover(ci[1] - r4049, roots, h, rinf);
                let n_committed = crossover(2.54e10, None, h, rinf);
                cross.push(json!({
                    "r_inf": label, "h": h,
                    "log2_n_star_new": n_mid, "log2_n_star_new_ci95": [n_lo, n_hi],
                    "log2_n_star_committed_c4_2.54e10": n_committed,
                }));
            }
        }
        cross.push(json!({
            "c4_new_at_4049": c, "c4_new_at_4049_ci95": ci,
            "roots_muls_at_1009": r1009, "roots_muls_at_4049": r4049,
            "model": "S/rho = sqrt(h) * (181/(4*126) * C(p)/p + r_inf); C(p) = (C4 at 4049 - roots at 4049) + roots linear in log2 p through 1009 and 4049; n* = p^4",
        }));
    }
    let out = json!({
        "schema": "gaudry_quartic_d4_cost.summary/v1",
        "configs": configs,
        "paired_ratios": pairs,
        "crossover": cross,
    });
    let path = dir.join("summary.json");
    fs::write(&path, serde_json::to_string_pretty(&out).unwrap() + "\n").unwrap();
    eprintln!("wrote {}", path.display());
}
