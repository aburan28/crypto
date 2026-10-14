//! # Tables for research/fiber_ic_scaling_20261011
//!
//! Reads the scaling sessions (one per cofactor) and, per curve, takes:
//! - rho's mean `S`;
//! - every IC arm's mean `S` over verified measured runs;
//! - the **envelope**: the cheapest closed-table arm on that curve.
//!
//! It then fits `log₂(S_envelope / S_rho)` against `log₂ #E` by least
//! squares, with a 95 percent interval from the residuals (Student t,
//! `n − 2` degrees of freedom), and prints the protocol's three tests.
//!
//! ```bash
//! cargo run --release --example fiber_scaling_tables -- \
//!   --dir research/fiber_ic_scaling_20261011 \
//!   --out research/fiber_ic_scaling_20261011/TABLES.md \
//!   --json research/fiber_ic_scaling_20261011/tables.json
//! ```

use std::collections::BTreeMap;
use std::fmt::Write as _;

use serde_json::{json, Value};

fn arg(args: &[String], name: &str) -> String {
    let i = args
        .iter()
        .position(|a| a == name)
        .unwrap_or_else(|| panic!("missing {name}"));
    args[i + 1].clone()
}

/// Two-sided 97.5 % Student t quantile for small degrees of freedom.
fn t975(df: usize) -> f64 {
    const T: [f64; 30] = [
        12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160,
        2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056,
        2.052, 2.048, 2.045, 2.042,
    ];
    if df == 0 {
        f64::NAN
    } else if df <= 30 {
        T[df - 1]
    } else {
        1.96
    }
}

/// Least squares `y = a + b·x`: `(b, half-width of b's 95 % interval, a)`.
fn fit(xs: &[f64], ys: &[f64]) -> (f64, f64, f64) {
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
    let se = (rss / (n - 2.0) / sxx).sqrt();
    (b, t975(xs.len() - 2) * se, a)
}

#[derive(Default)]
struct Arm {
    runs: u64,
    verified: u64,
    s_sum: f64,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let dir = arg(&args, "--dir");
    let mut md = String::new();
    let mut js = serde_json::Map::new();
    let mut by_h_size: BTreeMap<(String, u64), Vec<f64>> = BTreeMap::new();
    let mut sessions: Vec<String> = std::fs::read_dir(format!("{dir}/sessions"))
        .map(|d| {
            d.filter_map(|e| e.ok())
                .filter(|e| e.path().join("records.jsonl").exists())
                .map(|e| e.file_name().to_string_lossy().to_string())
                .collect()
        })
        .unwrap_or_default();
    sessions.sort();
    // `--combine NAME=a+b` also analyses sessions a and b as one (the same
    // curves and targets, more arms): the envelope then spans both.
    let mut groups: Vec<(String, Vec<String>)> = sessions
        .iter()
        .map(|s| (s.clone(), vec![s.clone()]))
        .collect();
    if let Some(i) = args.iter().position(|a| a == "--combine") {
        let (n, parts) = args[i + 1].split_once('=').expect("--combine NAME=a+b");
        groups.push((n.to_string(), parts.split('+').map(String::from).collect()));
    }
    for (name, members) in groups {
        let text: String = members
            .iter()
            .map(|m| std::fs::read_to_string(format!("{dir}/sessions/{m}/records.jsonl")).unwrap())
            .collect();
        // curve slug → (log2 N, bits, cofactor, arms)
        let mut curves: BTreeMap<(u64, String), (f64, u64, u64, BTreeMap<String, Arm>)> =
            BTreeMap::new();
        for line in text.lines() {
            let r: Value = serde_json::from_str(line).unwrap();
            if r["warmup"] == true {
                continue;
            }
            let c = &r["workload"]["curve"];
            let n = c["group_order"].as_u64().unwrap();
            let entry = curves
                .entry((n, c["slug"].as_str().unwrap().to_string()))
                .or_insert_with(|| {
                    (
                        (n as f64).log2(),
                        c["field_bits"].as_u64().unwrap_or(0),
                        c["cofactor"].as_u64().unwrap_or(0),
                        BTreeMap::new(),
                    )
                });
            let a = entry
                .3
                .entry(r["arm"].as_str().unwrap().to_string())
                .or_default();
            a.runs += 1;
            if r["outcome"]["status"] == "verified" {
                a.verified += 1;
                a.s_sum += r["cost"]["s"].as_f64().unwrap_or(0.0);
            }
        }
        let mean = |a: &Arm| (a.verified > 0).then(|| a.s_sum / a.verified as f64);
        let arm_names: Vec<String> = curves
            .values()
            .flat_map(|c| c.3.keys().cloned())
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect();
        writeln!(md, "## `{name}`\n").unwrap();
        write!(md, "| curve | bits | h | log₂ #E |").unwrap();
        for a in &arm_names {
            write!(md, " {a} |").unwrap();
        }
        writeln!(md, " envelope arm | envelope ÷ rho | blind ÷ envelope |").unwrap();
        write!(md, "|:--|--:|--:|--:|").unwrap();
        for _ in &arm_names {
            write!(md, "--:|").unwrap();
        }
        writeln!(md, ":--|--:|--:|").unwrap();
        let (mut xs, mut ys) = (Vec::new(), Vec::new());
        let mut rows = Vec::new();
        let mut h_seen = 0;
        let mut incomplete = 0;
        for ((_, slug), (lg, bits, h, arms)) in &curves {
            h_seen = *h;
            write!(md, "| `{slug}` | {bits} | {h} | {lg:.2} |").unwrap();
            for a in &arm_names {
                let cell = arms.get(a).map_or("–".to_string(), |x| {
                    let m = mean(x).map_or("–".into(), |s| format!("{s:.1}"));
                    if x.verified < x.runs {
                        format!("{m} ({}/{})", x.verified, x.runs)
                    } else {
                        m
                    }
                });
                write!(md, " {cell} |").unwrap();
            }
            let rho = arms.get("rho-neg").and_then(mean);
            let env = arms
                .iter()
                .filter(|(k, x)| k.starts_with("closed-") && x.verified == x.runs && x.runs > 0)
                .filter_map(|(k, x)| mean(x).map(|s| (s, k.clone())))
                .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            let blind = arms
                .get("blind-b128")
                .and_then(|x| (x.verified == x.runs).then(|| mean(x)).flatten());
            match (rho, &env) {
                (Some(rho), Some((s, k))) => {
                    let ratio = s / rho;
                    xs.push(*lg);
                    ys.push(ratio.log2());
                    by_h_size
                        .entry((name.clone(), *bits))
                        .or_default()
                        .push(ratio);
                    let bl = blind.map_or("–".into(), |b| format!("{:.2}", b / s));
                    writeln!(md, " {k} | {ratio:.2} | {bl} |").unwrap();
                    rows.push(json!({"slug": slug, "bits": bits, "cofactor": h, "log2_order": lg,
                        "rho_s": rho, "envelope_arm": k, "envelope_s": s, "ratio": ratio,
                        "blind_over_envelope": blind.map(|b| b / s),
                        "arms": arms.iter().map(|(k, x)| (k.clone(), json!({"runs": x.runs, "verified": x.verified, "s_mean": mean(x)}))).collect::<BTreeMap<_, _>>()}));
                }
                _ => {
                    incomplete += 1;
                    writeln!(md, " – | – | – |").unwrap();
                }
            }
        }
        let fitted = if xs.len() >= 3 {
            let (b, hw, a) = fit(&xs, &ys);
            writeln!(
                md,
                "\nEnvelope ÷ rho against size, h = {h_seen}: slope of log₂(ratio) on log₂ #E = **{b:.3}** (95 % interval {:.3} to {:.3}; prediction 1/6 = 0.167), over {} curves. Intercept {a:.2}: the fitted ratio would reach 1 at log₂ #E = {:.1}.{}\n",
                b - hw,
                b + hw,
                xs.len(),
                if b > 0.0 { -a / b } else { f64::NAN },
                if incomplete > 0 { format!(" {incomplete} curve(s) had no complete closed arm or no rho and are left out.") } else { String::new() }
            )
            .unwrap();
            Some(
                json!({"slope": b, "slope_low": b - hw, "slope_high": b + hw, "intercept": a, "curves": xs.len()}),
            )
        } else {
            writeln!(
                md,
                "\nFewer than three curves with complete arms: no fit.\n"
            )
            .unwrap();
            None
        };
        js.insert(
            name,
            json!({"cofactor": h_seen, "curves": rows, "fit": fitted}),
        );
    }
    // P-S2: cofactor 72 over cofactor 4 at matched size, for the
    // pre-registered h = 72 session and for any combined h = 72 analysis.
    let mut ps2 = Vec::new();
    let names72: std::collections::BTreeSet<String> = by_h_size
        .keys()
        .map(|(n, _)| n.clone())
        .filter(|n| n.starts_with("scaling_h72") && n != "scaling_h72_small_b")
        .collect();
    for n72 in names72 {
        writeln!(md, "## Cofactor 72 (`{n72}`) against cofactor 4 at matched size\n\n| bits | mean envelope ÷ rho, h = 4 | h = 72 | ratio |\n|--:|--:|--:|--:|").unwrap();
        let sizes: std::collections::BTreeSet<u64> = by_h_size.keys().map(|(_, b)| *b).collect();
        for b in sizes {
            let m = |n: &str| {
                by_h_size
                    .get(&(n.to_string(), b))
                    .map(|v| v.iter().sum::<f64>() / v.len() as f64)
            };
            if let (Some(a), Some(c)) = (m("scaling_h4"), m(&n72)) {
                writeln!(md, "| {b} | {a:.2} | {c:.2} | {:.2} |", c / a).unwrap();
                ps2.push(
                    json!({"h72_analysis": n72, "bits": b, "h4": a, "h72": c, "ratio": c / a}),
                );
            }
        }
        writeln!(md).unwrap();
    }
    js.insert("cofactor_comparison".into(), Value::Array(ps2));
    std::fs::write(arg(&args, "--out"), &md).expect("write");
    std::fs::write(
        arg(&args, "--json"),
        serde_json::to_string_pretty(&Value::Object(js)).unwrap() + "\n",
    )
    .expect("write");
    print!("{md}");
}
