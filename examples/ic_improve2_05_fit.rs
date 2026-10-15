//! IC improvement round 2, item 05: paired-ratio statistics and the frozen
//! θ-training rule for `ic_improve2_05_linear_harvest` output.
//!
//! `--ratio label:num/den` (repeatable): per `cell`, the geometric mean of
//! `num/den` over targets, with a 95 % percentile bootstrap CI (10,000
//! resamples of targets, seed 1) and the verdict (ADVANCES if lo > 1,
//! REGRESSES if hi < 1, else MATCHES).  Only rows with `"gate_ok": true`
//! are used.  `--medians a,b` prints per-cell medians of fields.
//!
//! `--train`: PROTOCOL.md "Training".  G_θ = r where gpv = gens/occ_vars < θ,
//! else i.  The candidates are +∞ and the midpoints of the per-cell median
//! gpv values.  Admissible means no cell's I/G_θ CI upper bound is below 1.
//! The rule picks the largest pooled geometric mean, with ties going to the
//! larger θ.
//!
//! ```bash
//! cargo run --release --example ic_improve2_05_fit -- --file runs/train_m3.jsonl --train
//! ```
use std::collections::BTreeMap;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::Value;

const B: usize = 10_000;

fn pct(v: &mut [f64], q: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[((v.len() as f64 - 1.0) * q).round() as usize]
}

fn gm_ci(logs: &[f64]) -> (f64, f64, f64) {
    let n = logs.len();
    let m = logs.iter().sum::<f64>() / n as f64;
    let mut rng = StdRng::seed_from_u64(1);
    let mut boots: Vec<f64> = (0..B)
        .map(|_| (0..n).map(|_| logs[rng.gen_range(0..n)]).sum::<f64>() / n as f64)
        .collect();
    (m.exp(), pct(&mut boots, 0.025).exp(), pct(&mut boots, 0.975).exp())
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
    if k == "gpv" {
        return f(r, "gens") / f(r, "occ_vars");
    }
    r.get(k)
        .and_then(Value::as_f64)
        .unwrap_or_else(|| panic!("field {k}"))
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let (mut files, mut ratios, mut medians, mut train) = (vec![], vec![], vec![], false);
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
            "--medians" => {
                i += 1;
                medians.extend(argv[i].split(',').map(str::to_string));
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
    let ok = |r: &Value| r.get("gate_ok").and_then(Value::as_bool) == Some(true);
    let failed = rows.iter().filter(|r| !ok(r)).count();
    let mut cells: BTreeMap<String, Vec<&Value>> = BTreeMap::new();
    for r in rows.iter().filter(|r| ok(r)) {
        cells
            .entry(r["cell"].as_str().unwrap().to_string())
            .or_default()
            .push(r);
    }
    eprintln!("rows {} failed {failed}", rows.len());
    for spec in &ratios {
        let (label, fields) = spec.split_once(':').expect("label:num/den");
        let (num, den) = fields.split_once('/').expect("num/den");
        for (cell, rs) in &cells {
            let logs: Vec<f64> = rs.iter().map(|r| (f(r, num) / f(r, den)).ln()).collect();
            let (gm, lo, hi) = gm_ci(&logs);
            println!(
                r#"{{"ratio":"{label}","cell":"{cell}","targets":{},"gm":{gm:.4},"lo":{lo:.4},"hi":{hi:.4},"verdict":"{}","gpv_median":{:.4},"i_cpu_median":{:.6},"failed_total":{failed}}}"#,
                rs.len(),
                verdict(lo, hi),
                median(rs.iter().map(|r| f(r, "gpv")).collect()),
                median(rs.iter().map(|r| f(r, "i_cpu")).collect()),
            );
        }
    }
    if !medians.is_empty() {
        for (cell, rs) in &cells {
            let parts: Vec<String> = medians
                .iter()
                .map(|k| format!(r#""{k}":{:.6}"#, median(rs.iter().map(|r| f(r, k)).collect())))
                .collect();
            println!(r#"{{"medians":"{cell}","targets":{},{}}}"#, rs.len(), parts.join(","));
        }
    }
    if train {
        let mut gpvs: Vec<f64> = cells
            .values()
            .map(|rs| median(rs.iter().map(|r| f(r, "gpv")).collect()))
            .collect();
        gpvs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        gpvs.dedup();
        let mut candidates: Vec<f64> = gpvs.windows(2).map(|w| 0.5 * (w[0] + w[1])).collect();
        candidates.push(f64::INFINITY);
        let mut best: Option<(f64, f64)> = None;
        let mut scan = Vec::new();
        for &theta in &candidates {
            let mut admissible = true;
            let mut pooled = Vec::new();
            let mut per = Vec::new();
            for (cell, rs) in &cells {
                let logs: Vec<f64> = rs
                    .iter()
                    .map(|r| {
                        let g = if f(r, "gpv") < theta { f(r, "r_cpu") } else { f(r, "i_cpu") };
                        (f(r, "i_cpu") / g).ln()
                    })
                    .collect();
                let (gm, lo, hi) = gm_ci(&logs);
                if hi < 1.0 {
                    admissible = false;
                }
                pooled.extend(logs);
                per.push(format!(r#"{{"cell":"{cell}","gm":{gm:.4},"lo":{lo:.4},"hi":{hi:.4}}}"#));
            }
            let score = (pooled.iter().sum::<f64>() / pooled.len() as f64).exp();
            let th = if theta.is_finite() { format!("{theta:.4}") } else { "null".into() };
            scan.push(format!(
                r#"{{"theta":{th},"admissible":{admissible},"pooled_gm":{score:.4},"cells":[{}]}}"#,
                per.join(",")
            ));
            if admissible
                && best.is_none_or(|(bt, bs)| score > bs + 1e-12 || ((score - bs).abs() <= 1e-12 && theta > bt))
            {
                best = Some((theta, score));
            }
        }
        let chosen = match best {
            Some((t, s)) if t.is_finite() => format!(r#""theta":{t:.4},"pooled_gm":{s:.4}"#),
            Some((_, s)) => format!(r#""theta":null,"pooled_gm":{s:.4}"#),
            None => r#""theta":"none admissible""#.to_string(),
        };
        println!(r#"{{"train":true,{chosen},"scan":[{}]}}"#, scan.join(","));
    }
}
