//! IC improvement round 2, item 03: paired-ratio statistics and gate
//! training for `ic_improve2_03_fallrow_gate` output.
//!
//! `--ratio label:num/den` (repeatable): per `cell`, the geometric mean of
//! `num/den` over targets with a 95 % percentile bootstrap CI (10,000
//! resamples of targets, seed 1), and the verdict (ADVANCES if lo > 1,
//! REGRESSES if hi < 1, else MATCHES).  Only rows with `"gate_ok": true`
//! are used; failed rows are counted.  One JSON line per (ratio, cell).
//!
//! `--train`: the frozen procedure of PROTOCOL.md §5 on train-mode rows;
//! prints the chosen `θ` and `π₀` with the table it was chosen from.
//!
//! ```bash
//! cargo run --release --example ic_improve2_03_fit -- --file runs/test.jsonl --ratio I/G:i_cpu/g_cpu
//! cargo run --release --example ic_improve2_03_fit -- --file runs/train.jsonl --train
//! ```

use std::collections::BTreeMap;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::Value;

const B: usize = 10_000;

fn pct(v: &mut [f64], q: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let i = ((v.len() as f64 - 1.0) * q).round() as usize;
    v[i]
}

/// Geometric mean of `xs` (given as logs) and its 95 % bootstrap CI.
fn gm_ci(logs: &[f64]) -> (f64, f64, f64) {
    let n = logs.len();
    let m = logs.iter().sum::<f64>() / n as f64;
    let mut rng = StdRng::seed_from_u64(1);
    let mut boots: Vec<f64> = (0..B)
        .map(|_| (0..n).map(|_| logs[rng.gen_range(0..n)]).sum::<f64>() / n as f64)
        .collect();
    let lo = pct(&mut boots, 0.025);
    let hi = pct(&mut boots, 0.975);
    (m.exp(), lo.exp(), hi.exp())
}

fn verdict(lo: f64, hi: f64) -> &'static str {
    if lo > 1.0 {
        "ADVANCES"
    } else if hi < 1.0 {
        "REGRESSES"
    } else {
        "MATCHES"
    }
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let k = v.len();
    if k % 2 == 1 {
        v[k / 2]
    } else {
        0.5 * (v[k / 2 - 1] + v[k / 2])
    }
}

fn f(r: &Value, k: &str) -> f64 {
    r.get(k).and_then(Value::as_f64).unwrap_or_else(|| panic!("field {k}"))
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let mut files = Vec::new();
    let mut ratios = Vec::new();
    let mut train = false;
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
            "--train" => train = true,
            _ => {}
        }
        i += 1;
    }
    let mut rows: Vec<Value> = Vec::new();
    for file in &files {
        for line in std::fs::read_to_string(file).expect("read").lines() {
            if let Ok(v) = serde_json::from_str::<Value>(line) {
                rows.push(v);
            }
        }
    }
    let failed = rows
        .iter()
        .filter(|r| r.get("gate_ok").and_then(Value::as_bool) != Some(true))
        .count();
    let good: Vec<&Value> = rows
        .iter()
        .filter(|r| r.get("gate_ok").and_then(Value::as_bool) == Some(true))
        .collect();
    let mut cells: BTreeMap<String, Vec<&Value>> = BTreeMap::new();
    for r in &good {
        cells
            .entry(r["cell"].as_str().unwrap().to_string())
            .or_default()
            .push(r);
    }
    eprintln!("rows {} ok {} failed {failed}", rows.len(), good.len());

    for spec in &ratios {
        let (label, fields) = spec.split_once(':').expect("label:num/den");
        let (num, den) = fields.split_once('/').expect("num/den");
        for (cell, rs) in &cells {
            let logs: Vec<f64> = rs.iter().map(|r| (f(r, num) / f(r, den)).ln()).collect();
            let (gm, lo, hi) = gm_ci(&logs);
            let fires = rs
                .iter()
                .filter(|r| r.get("fires").and_then(Value::as_bool) == Some(true))
                .count();
            let gpv = median(rs.iter().map(|r| f(r, "gens_per_var")).collect());
            let i_med = median(rs.iter().map(|r| f(r, "i_cpu")).collect());
            println!(
                r#"{{"ratio":"{label}","cell":"{cell}","targets":{},"gm":{gm:.4},"lo":{lo:.4},"hi":{hi:.4},"verdict":"{}","fires":{fires},"gpv_median":{gpv:.4},"i_cpu_median":{i_med:.6},"failed_total":{failed}}}"#,
                rs.len(),
                verdict(lo, hi)
            );
        }
    }

    if train {
        // §5.1: the non-firing path.
        let mut pooled = Vec::new();
        let mut lr_ok = true;
        let mut table: Vec<(String, f64, (f64, f64, f64), (f64, f64, f64))> = Vec::new();
        for (cell, rs) in &cells {
            let lr: Vec<f64> = rs.iter().map(|r| (f(r, "i_cpu") / f(r, "lr_cpu")).ln()).collect();
            let fl: Vec<f64> = rs.iter().map(|r| (f(r, "i_cpu") / f(r, "f_cpu")).ln()).collect();
            pooled.extend(lr.iter().copied());
            let lr_s = gm_ci(&lr);
            let f_s = gm_ci(&fl);
            if lr_s.2 < 1.0 {
                lr_ok = false;
            }
            let gpv = median(rs.iter().map(|r| f(r, "gens_per_var")).collect());
            table.push((cell.clone(), gpv, f_s, lr_s));
        }
        let pooled_gm = (pooled.iter().sum::<f64>() / pooled.len() as f64).exp();
        let pi0 = if lr_ok && pooled_gm >= 1.0 { "lr" } else { "i" };
        // §5.2: the threshold.
        let mut gpvs: Vec<f64> = table.iter().map(|t| t.1).collect();
        gpvs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        gpvs.dedup();
        let mut candidates: Vec<f64> = gpvs.windows(2).map(|w| 0.5 * (w[0] + w[1])).collect();
        candidates.push(f64::INFINITY);
        let mut best: Option<(f64, f64)> = None;
        let mut scan = Vec::new();
        for &theta in &candidates {
            let mut admissible = true;
            let mut score = 0.0;
            for (_, gpv, f_s, lr_s) in &table {
                let s = if *gpv >= theta {
                    *f_s
                } else if pi0 == "lr" {
                    *lr_s
                } else {
                    (1.0, 1.0, 1.0)
                };
                if s.2 < 1.0 {
                    admissible = false;
                }
                score += s.0.ln();
            }
            scan.push(format!(
                r#"{{"theta":{},"admissible":{admissible},"score":{score:.4}}}"#,
                if theta.is_finite() { format!("{theta:.4}") } else { "null".into() }
            ));
            if admissible {
                let better = match best {
                    None => true,
                    Some((bt, bs)) => score > bs + 1e-12 || ((score - bs).abs() <= 1e-12 && theta > bt),
                };
                if better {
                    best = Some((theta, score));
                }
            }
        }
        let cells_json: Vec<String> = table
            .iter()
            .map(|(c, g, fs, ls)| {
                format!(
                    r#"{{"cell":"{c}","gpv_median":{g:.4},"I/F":[{:.4},{:.4},{:.4}],"I/LR":[{:.4},{:.4},{:.4}]}}"#,
                    fs.0, fs.1, fs.2, ls.0, ls.1, ls.2
                )
            })
            .collect();
        let (theta, score) = best.expect("+inf with pi0 admissible or none");
        println!(
            r#"{{"pi0":"{pi0}","lr_pooled_gm":{pooled_gm:.4},"lr_no_regression":{lr_ok},"theta":{},"score":{score:.4},"scan":[{}],"cells":[{}]}}"#,
            if theta.is_finite() { format!("{theta:.4}") } else { "null".into() },
            scan.join(","),
            cells_json.join(",")
        );
    }
}
