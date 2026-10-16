//! IC improvement item 02: paired-ratio statistics for
//! `ic_improve_02_equivariant_oracle` output.
//!
//! For each `--ratio label:num_field/den_field`, per size `n`: the geometric
//! mean of `num/den` over targets with a 95 % percentile bootstrap CI
//! (10,000 resamples of targets, seed 1), then the OLS slope of
//! `ln ratio` on `ln n` over all targets with a bootstrap CI (targets
//! resampled within each size).  Only rows with `"gate": true` are used;
//! the count of failed gates is printed.  Output: one JSON line per ratio.
//!
//! ```bash
//! cargo run --release --example ic_improve_02_fit -- --file runs/oracle.jsonl \
//!     --ratio I/C:i_cpu/c_cpu --ratio AA:i_cpu/i2_cpu
//! ```

use std::collections::BTreeMap;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

const B: usize = 10_000;

fn pct(v: &mut [f64], q: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let i = ((v.len() as f64 - 1.0) * q).round() as usize;
    v[i]
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

fn slope(points: &[(f64, f64)]) -> f64 {
    let mx = points.iter().map(|p| p.0).sum::<f64>() / points.len() as f64;
    let my = points.iter().map(|p| p.1).sum::<f64>() / points.len() as f64;
    let sxy: f64 = points.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let sxx: f64 = points.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum();
    sxy / sxx
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let mut files = Vec::new();
    let mut ratios = Vec::new();
    let mut gate_field = "gate".to_string();
    let mut i = 1;
    while i < argv.len() {
        match argv[i].as_str() {
            "--file" => {
                i += 1;
                files.push(argv[i].clone());
            }
            "--ratio" => {
                i += 1;
                ratios.push(argv[i].clone());
            }
            "--gate" => {
                i += 1;
                gate_field = argv[i].clone();
            }
            _ => {}
        }
        i += 1;
    }
    let mut rows: Vec<serde_json::Value> = Vec::new();
    for f in &files {
        for line in std::fs::read_to_string(f).expect("read").lines() {
            if let Ok(v) = serde_json::from_str::<serde_json::Value>(line) {
                rows.push(v);
            }
        }
    }
    let failed = rows
        .iter()
        .filter(|r| r.get(&gate_field).and_then(|g| g.as_bool()) == Some(false))
        .count();
    for spec in &ratios {
        let (label, fields) = spec.split_once(':').expect("label:num/den");
        let (num, den) = fields.split_once('/').expect("num/den");
        let mut by_n: BTreeMap<u64, Vec<f64>> = BTreeMap::new();
        for r in &rows {
            let (Some(a), Some(b)) = (r[num].as_f64(), r[den].as_f64()) else {
                continue;
            };
            if a <= 0.0 || b <= 0.0 {
                continue;
            }
            by_n.entry(r["n"].as_u64().unwrap())
                .or_default()
                .push((a / b).ln());
        }
        let mut rng = StdRng::seed_from_u64(1);
        let mut per_n = Vec::new();
        for (n, logs) in &by_n {
            let est = mean(logs).exp();
            let mut boots: Vec<f64> = (0..B)
                .map(|_| {
                    let s: f64 = (0..logs.len())
                        .map(|_| logs[rng.gen_range(0..logs.len())])
                        .sum();
                    (s / logs.len() as f64).exp()
                })
                .collect();
            let lo = pct(&mut boots, 0.025);
            let hi = pct(&mut boots, 0.975);
            per_n.push(serde_json::json!({
                "n": n, "k": logs.len(), "geomean": est, "lo": lo, "hi": hi,
                "min": logs.iter().cloned().fold(f64::INFINITY, f64::min).exp(),
                "max": logs.iter().cloned().fold(f64::NEG_INFINITY, f64::max).exp(),
            }));
        }
        let pts: Vec<(f64, f64)> = by_n
            .iter()
            .flat_map(|(n, logs)| logs.iter().map(move |&y| ((*n as f64).ln(), y)))
            .collect();
        let (b_est, b_lo, b_hi) = if by_n.len() >= 2 {
            let est = slope(&pts);
            let mut boots: Vec<f64> = (0..B)
                .map(|_| {
                    let mut p = Vec::with_capacity(pts.len());
                    for (n, logs) in &by_n {
                        for _ in 0..logs.len() {
                            p.push(((*n as f64).ln(), logs[rng.gen_range(0..logs.len())]));
                        }
                    }
                    slope(&p)
                })
                .collect();
            (est, pct(&mut boots, 0.025), pct(&mut boots, 0.975))
        } else {
            (f64::NAN, f64::NAN, f64::NAN)
        };
        println!(
            "{}",
            serde_json::json!({
                "ratio": label, "num": num, "den": den, "gate_failures": failed,
                "per_n": per_n,
                "slope_vs_ln_n": if b_est.is_nan() { serde_json::Value::Null } else { serde_json::json!({"b": b_est, "lo": b_lo, "hi": b_hi}) },
            })
        );
    }
}
