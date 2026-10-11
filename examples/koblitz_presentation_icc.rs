//! **Curve-vs-presentation variance split** for `koblitz_presentation_sweep`
//! output: the Rust port of `research/notes/koblitz-isogeny/presentation_icc.py`
//! (same model, estimator, metrics and percentile CI; the bootstrap draws
//! come from a different RNG, so CI endpoints differ in the last digits).
//!
//! Two-way crossed random effects without replication,
//! `Y[E,V] = μ + u_E + w_V + e`, by ANOVA mean squares on the complete
//! curves × subspaces grid; `ICC_E = σ²_u / (σ²_u + σ²_w + σ²_e)`; 95 % CI
//! by 2000 bootstrap resamples of curves.
//!
//! Metrics (per cell): yield per `|F|^m`; `ln µs` per call; work per call
//! (Gröbner reductions for m = 2, `S₄` roots for the m = 3 e₃ engine);
//! `ln` projected solve `(unknowns + 4) / yield × time per call`.
//!
//! Run: `cargo run --release --example koblitz_presentation_icc <pres.jsonl>...`

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};

type Metric = (&'static str, fn(&Value) -> f64);

fn f(r: &Value, k: &str) -> f64 {
    r[k].as_f64().unwrap_or(0.0)
}

fn m_of(r: &Value) -> i32 {
    r["m"].as_i64().unwrap_or(2) as i32
}

const METRICS: [Metric; 4] = [
    ("yield_per_Fm", |r| {
        (f(r, "relations") / f(r, "probes")) / f(r, "factor_base_points").powi(m_of(r))
    }),
    ("log_us_per_call", |r| {
        (f(r, "groebner_ns") / f(r, "groebner_calls") / 1e3).ln()
    }),
    ("work_per_call", |r| {
        let work = if m_of(r) == 3 {
            f(r, "sym_roots")
        } else {
            f(r, "reductions")
        };
        work / f(r, "groebner_calls")
    }),
    ("log_projected_solve", |r| {
        let y = (f(r, "relations") / f(r, "probes")).max(1e-9);
        ((f(r, "unknowns") + 4.0) / y * f(r, "groebner_ns") / f(r, "groebner_calls")).ln()
    }),
];

/// `(ICC_E, σ²_u, σ²_w, σ²_e)` of a curves × subspaces grid.
fn icc(grid: &[Vec<f64>]) -> (f64, f64, f64, f64) {
    let (a, b) = (grid.len() as f64, grid[0].len() as f64);
    let gm = grid.iter().flatten().sum::<f64>() / (a * b);
    let rm: Vec<f64> = grid.iter().map(|r| r.iter().sum::<f64>() / b).collect();
    let cm: Vec<f64> = (0..grid[0].len())
        .map(|j| grid.iter().map(|r| r[j]).sum::<f64>() / a)
        .collect();
    let ss_r = b * rm.iter().map(|x| (x - gm).powi(2)).sum::<f64>();
    let ss_c = a * cm.iter().map(|x| (x - gm).powi(2)).sum::<f64>();
    let ss_t: f64 = grid.iter().flatten().map(|x| (x - gm).powi(2)).sum();
    let (ms_r, ms_c) = (ss_r / (a - 1.0), ss_c / (b - 1.0));
    let ms_e = (ss_t - ss_r - ss_c) / ((a - 1.0) * (b - 1.0));
    let su = ((ms_r - ms_e) / b).max(0.0);
    let sw = ((ms_c - ms_e) / a).max(0.0);
    let tot = su + sw + ms_e;
    (if tot > 0.0 { su / tot } else { 0.0 }, su, sw, ms_e)
}

fn main() {
    let mut groups: BTreeMap<(u64, u64, u64, i64), Vec<Value>> = BTreeMap::new();
    for path in std::env::args().skip(1) {
        for line in std::fs::read_to_string(&path).expect("read").lines() {
            let r: Value = serde_json::from_str(line).expect("json");
            if r.get("skipped").is_some()
                || r["inconsistent"].as_u64().unwrap_or(0) > 0
                || !r["bsgs_ok"].as_bool().unwrap_or(true)
            {
                continue;
            }
            let key = (
                r["n"].as_u64().unwrap(),
                r["a2"].as_u64().unwrap(),
                r["l"].as_u64().unwrap(),
                r["m"].as_i64().unwrap_or(2),
            );
            groups.entry(key).or_default().push(r);
        }
    }
    let mut out = serde_json::Map::new();
    for ((n, a2, l, m), rs) in groups {
        let seeds: BTreeSet<u64> = rs
            .iter()
            .map(|r| r["basis_seed"].as_u64().unwrap())
            .collect();
        let mut by: BTreeMap<u64, BTreeMap<u64, &Value>> = BTreeMap::new();
        for r in &rs {
            by.entry(r["a6"].as_u64().unwrap())
                .or_default()
                .insert(r["basis_seed"].as_u64().unwrap(), r);
        }
        let curves: Vec<u64> = by
            .iter()
            .filter(|(_, v)| seeds.iter().all(|s| v.contains_key(s)))
            .map(|(c, _)| *c)
            .collect();
        let mut res = serde_json::Map::new();
        res.insert("curves".into(), json!(curves.len()));
        res.insert("subspaces".into(), json!(seeds.len()));
        for (name, metric) in METRICS {
            let grid: Vec<Vec<f64>> = curves
                .iter()
                .map(|c| seeds.iter().map(|s| metric(by[c][s])).collect())
                .collect();
            let (point, su, sw, se) = icc(&grid);
            let mut rng = StdRng::seed_from_u64(7);
            let mut boots: Vec<f64> = (0..2000)
                .map(|_| {
                    let g: Vec<Vec<f64>> = (0..grid.len())
                        .map(|_| grid[rng.gen_range(0..grid.len())].clone())
                        .collect();
                    icc(&g).0
                })
                .collect();
            boots.sort_by(f64::total_cmp);
            let tot = su + sw + se;
            res.insert(
                name.into(),
                json!({"icc_curve": point, "ci95": [boots[50], boots[1949]],
                       "var_curve": su, "var_subspace": sw, "var_residual": se,
                       "share_subspace": if tot > 0.0 { sw / tot } else { 0.0 }}),
            );
        }
        out.insert(format!("n={n} a2={a2} l={l} m={m}"), Value::Object(res));
    }
    println!(
        "{}",
        serde_json::to_string_pretty(&Value::Object(out)).unwrap()
    );
}
