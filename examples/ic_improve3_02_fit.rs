//! IC improvement round 3, item 02: paired-ratio statistics, the frozen
//! gate-training rule and counter-identity checks for
//! `ic_improve3_02_harvest_gate` output.
//!
//! * `--ratio label:num/den` (repeatable): per cell, the geometric mean of
//!   `num/den` over targets with a 95 % percentile bootstrap CI (10,000
//!   resamples of targets, seed 1) and the verdict (ADVANCES if lo > 1,
//!   REGRESSES if hi < 1, else MATCHES).  Rows with `"gate_ok": false` are
//!   counted and dropped.
//! * `--medians a,b`: per-cell medians of fields (`gpv` = gens/occ_vars).
//! * `--train`: PROTOCOL.md "Training".  C(f, θ) fires where G0 fires
//!   (gpv < 0.9607) or f ≥ θ, f ∈ {rank (`root_rank`), occ (`occ_vars`)}.
//!   Per target the ratio G0/C is 1 where C takes G0's path and `g/r`
//!   where C adds firing.  Admissible: no cell's CPU or word-op CI upper
//!   bound below 1.  Choice: the largest pooled CPU geometric mean; ties to
//!   the larger θ, then `rank` before `occ`.
//! * `--same a=b` (repeatable): rows whose arms agree on every counter.
//! * `--ref FILE --ref-arm r --arm a`: rows whose arm `a` agrees with FILE's
//!   arm `r` (matched by cell without `pipe-` and target) on every counter
//!   FILE records.
//! * `--fires a:SPEC`: check arm `a` harvested exactly where the frozen
//!   rule SPEC (`rank=θ` or `occ=θ`) predicts from the logged features.
//!
//! ```bash
//! cargo run --release --example ic_improve3_02_fit -- --file runs/train_a.jsonl --train
//! ```
use std::collections::BTreeMap;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::Value;

const B: usize = 10_000;
const GPV_THETA: f64 = 0.9607;

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
    (
        m.exp(),
        pct(&mut boots, 0.025).exp(),
        pct(&mut boots, 0.975).exp(),
    )
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
    match k {
        "gpv" => f(r, "gens") / f(r, "occ_vars"),
        "rank" => f(r, "root_rank"),
        "occ" => f(r, "occ_vars"),
        _ => r
            .get(k)
            .and_then(Value::as_f64)
            .unwrap_or_else(|| panic!("field {k}")),
    }
}

/// Does row `r` carry field `k` (derived fields count as present)?
fn has(r: &Value, k: &str) -> bool {
    matches!(k, "gpv" | "rank" | "occ") || r.get(k).is_some()
}

fn g0_fires(r: &Value) -> bool {
    f(r, "gpv") < GPV_THETA
}

/// Does C(feature, θ) add firing on this target (fires where G0 does not)?
fn adds(r: &Value, feature: &str, theta: f64) -> bool {
    !g0_fires(r) && f(r, feature) >= theta
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let (mut files, mut ratios, mut medians, mut train) = (vec![], vec![], vec![], false);
    let (mut same, mut refs, mut cmp_arms, mut fires): (
        Vec<String>,
        Vec<String>,
        Vec<String>,
        Vec<String>,
    ) = (vec![], vec![], vec![], vec![]);
    let mut ref_arm = String::from("i");
    let mut i = 1;
    while i < argv.len() {
        let next = |i: &mut usize| {
            *i += 1;
            argv[*i].clone()
        };
        match argv[i].as_str() {
            "--file" => files.push(next(&mut i)),
            "--ratio" => ratios.push(next(&mut i)),
            "--medians" => medians.extend(next(&mut i).split(',').map(str::to_string)),
            "--train" => train = true,
            "--same" => same.push(next(&mut i)),
            "--ref" => refs.push(next(&mut i)),
            "--ref-arm" => ref_arm = next(&mut i),
            "--arm" => cmp_arms.push(next(&mut i)),
            "--fires" => fires.push(next(&mut i)),
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
        let key = format!(
            "{}{}",
            r["cell"].as_str().unwrap(),
            r.get("mode")
                .and_then(Value::as_str)
                .map(|m| format!("/{m}"))
                .unwrap_or_default()
        );
        cells.entry(key).or_default().push(r);
    }
    eprintln!("rows {} failed {failed}", rows.len());
    for spec in &ratios {
        let (label, fields) = spec.split_once(':').expect("label:num/den");
        let (num, den) = fields.split_once('/').expect("num/den");
        for (cell, rs) in &cells {
            let rs: Vec<&Value> = rs.iter().copied().filter(|r| has(r, num) && has(r, den)).collect();
            if rs.is_empty() {
                continue;
            }
            let logs: Vec<f64> = rs.iter().map(|r| (f(r, num) / f(r, den)).ln()).collect();
            let (gm, lo, hi) = gm_ci(&logs);
            println!(
                r#"{{"ratio":"{label}","cell":"{cell}","targets":{},"gm":{gm:.4},"lo":{lo:.4},"hi":{hi:.4},"verdict":"{}","i_cpu_median":{:.6},"failed_total":{failed}}}"#,
                rs.len(),
                verdict(lo, hi),
                median(rs.iter().map(|r| f(r, "i_cpu")).collect()),
            );
        }
    }
    if !medians.is_empty() {
        for (cell, rs) in &cells {
            let parts: Vec<String> = medians
                .iter()
                .filter(|k| rs.iter().all(|r| has(r, k)))
                .map(|k| {
                    format!(
                        r#""{k}":{:.6}"#,
                        median(rs.iter().map(|r| f(r, k)).collect())
                    )
                })
                .collect();
            println!(
                r#"{{"medians":"{cell}","targets":{},{}}}"#,
                rs.len(),
                parts.join(",")
            );
        }
    }
    const COUNTERS: [&str; 8] = [
        "splits",
        "reductions",
        "elim",
        "word_ops",
        "calls",
        "spec_ops",
        "roots",
        "root_hash",
    ];
    for spec in &same {
        let (a, b) = spec.split_once('=').expect("a=b");
        let agree = rows
            .iter()
            .filter(|r| {
                COUNTERS.iter().all(|c| {
                    let (x, y) = (r.get(format!("{a}_{c}")), r.get(format!("{b}_{c}")));
                    x == y
                }) && r.get(format!("{a}_found")) == r.get(format!("{b}_found"))
            })
            .count();
        println!(
            r#"{{"same":"{a}={b}","rows":{},"identical":{agree}}}"#,
            rows.len()
        );
    }
    if !refs.is_empty() {
        let mut reference: BTreeMap<(String, u64), Value> = BTreeMap::new();
        for file in &refs {
            for line in std::fs::read_to_string(file).expect("read ref").lines() {
                if let Ok(v) = serde_json::from_str::<Value>(line) {
                    if let (Some(c), Some(t)) = (v["cell"].as_str(), v["target"].as_u64()) {
                        reference.insert((c.trim_start_matches("pipe-").to_string(), t), v.clone());
                    }
                }
            }
        }
        for a in &cmp_arms {
            let (mut matched, mut agree) = (0, 0);
            for r in &rows {
                let key = (
                    r["cell"]
                        .as_str()
                        .unwrap()
                        .trim_start_matches("pipe-")
                        .to_string(),
                    r["target"].as_u64().unwrap(),
                );
                let Some(o) = reference.get(&key) else {
                    continue;
                };
                matched += 1;
                if COUNTERS.iter().all(|c| {
                    let rk = format!("{ref_arm}_{c}");
                    o.get(&rk).is_none_or(|ov| Some(ov) == r.get(format!("{a}_{c}")))
                }) {
                    agree += 1;
                }
            }
            println!(
                r#"{{"ref_arm":"{ref_arm}","arm":"{a}","matched":{matched},"identical":{agree}}}"#
            );
        }
    }
    for spec in &fires {
        let (arm, rule) = spec.split_once(':').expect("arm:rule");
        let (feature, theta) = rule.split_once('=').expect("feature=θ");
        let theta: f64 = theta.parse().expect("θ");
        let mut agree = 0;
        let rows: Vec<&Value> = rows.iter().filter(|r| has(r, "r_harvested_rows")).collect();
        for r in &rows {
            let predicted = g0_fires(r) || f(r, feature) >= theta;
            let fired = f(r, &format!("{arm}_harvests")) > 0.0;
            // A firing gate harvests only when the root's tail has rows.
            let has_rows = f(r, "r_harvested_rows") > 0.0;
            if fired == (predicted && has_rows) {
                agree += 1;
            }
        }
        println!(
            r#"{{"fires":"{spec}","rows":{},"agree":{agree}}}"#,
            rows.len()
        );
    }
    if train {
        let mut scan = Vec::new();
        // (feature, θ, score) of the best admissible candidate.
        let mut best: Option<(String, f64, f64)> = None;
        for feature in ["rank", "occ"] {
            let mut values: Vec<f64> = cells
                .values()
                .map(|rs| median(rs.iter().map(|r| f(r, feature)).collect()))
                .collect();
            values.sort_by(|a, b| a.partial_cmp(b).unwrap());
            values.dedup();
            let mut candidates: Vec<f64> = vec![values[0] - 1.0];
            candidates.extend(values.windows(2).map(|w| (0.5 * (w[0] + w[1])).floor()));
            candidates.push(f64::INFINITY);
            for &theta in &candidates {
                let mut admissible = true;
                let mut pooled = Vec::new();
                let mut per = Vec::new();
                for (cell, rs) in &cells {
                    let ratio = |r: &&Value, num: &str, den: &str| {
                        if adds(r, feature, theta) {
                            (f(r, num) / f(r, den)).ln()
                        } else {
                            0.0
                        }
                    };
                    let cpu: Vec<f64> = rs.iter().map(|r| ratio(r, "g_cpu", "r_cpu")).collect();
                    let wo: Vec<f64> = rs
                        .iter()
                        .map(|r| ratio(r, "g_word_ops", "r_word_ops"))
                        .collect();
                    let (gm, lo, hi) = gm_ci(&cpu);
                    let (wgm, wlo, whi) = gm_ci(&wo);
                    if hi < 1.0 || whi < 1.0 {
                        admissible = false;
                    }
                    let added = rs.iter().filter(|r| adds(r, feature, theta)).count();
                    pooled.extend(cpu);
                    per.push(format!(
                        r#"{{"cell":"{cell}","added":{added},"gm":{gm:.4},"lo":{lo:.4},"hi":{hi:.4},"w_gm":{wgm:.4},"w_lo":{wlo:.4},"w_hi":{whi:.4}}}"#
                    ));
                }
                let score = (pooled.iter().sum::<f64>() / pooled.len() as f64).exp();
                let th = if theta.is_finite() {
                    format!("{theta}")
                } else {
                    "null".into()
                };
                scan.push(format!(
                    r#"{{"feature":"{feature}","theta":{th},"admissible":{admissible},"pooled_gm":{score:.4},"cells":[{}]}}"#,
                    per.join(",")
                ));
                if admissible && theta.is_finite() {
                    let better = match &best {
                        None => true,
                        Some((bf, bt, bs)) => {
                            score > bs + 1e-12
                                || ((score - bs).abs() <= 1e-12
                                    && (theta > *bt
                                        || (theta == *bt && feature == "rank" && bf == "occ")))
                        }
                    };
                    if better {
                        best = Some((feature.to_string(), theta, score));
                    }
                }
            }
        }
        let chosen = match best {
            Some((feat, t, s)) => {
                format!(r#""feature":"{feat}","theta":{t},"pooled_gm":{s:.4},"spec":"root,gpv=0.9607,{feat}={t}""#)
            }
            None => r#""feature":null,"theta":"none admissible""#.to_string(),
        };
        println!(r#"{{"train":true,{chosen}}}"#);
        for line in scan {
            println!("{line}");
        }
    }
}
