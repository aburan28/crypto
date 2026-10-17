//! lit2026 W11: summary statistics of a run directory (medians, min–max,
//! per-round ratios, identity checks), in Rust so the plotting script only
//! draws recorded values.
//!
//! ```sh
//! cargo run --release --example lit2026_w11_summarize -- <run_dir> [classwalk_dir]
//! ```
//!
//! Reads `<run_dir>/raw/{bsgs_time.jsonl, binary_sect163.jsonl}` and, when
//! given, the class-walk session directory (`samples.tsv` and each run's
//! `experiments/koblitz_isogeny_class_walk.json`); writes
//! `<run_dir>/raw/summary.json`.

use std::collections::BTreeMap;
use std::path::Path;

use serde_json::{json, Value};

fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 {
        s[n / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}

fn mm(v: &[f64]) -> (f64, f64) {
    (
        v.iter().cloned().fold(f64::INFINITY, f64::min),
        v.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    )
}

fn jsonl(p: &Path) -> Vec<Value> {
    std::fs::read_to_string(p)
        .map(|s| {
            s.lines()
                .filter(|l| l.starts_with('{'))
                .map(|l| serde_json::from_str(l).expect("json line"))
                .collect()
        })
        .unwrap_or_default()
}

fn strip_seconds(v: &mut Value) {
    match v {
        Value::Object(m) => {
            m.remove("seconds");
            for x in m.values_mut() {
                strip_seconds(x);
            }
        }
        Value::Array(a) => {
            for x in a {
                strip_seconds(x);
            }
        }
        _ => {}
    }
}

/// Per-round ratios `num/den` matched by round.
fn round_ratios(num: &BTreeMap<u64, f64>, den: &BTreeMap<u64, f64>) -> Vec<f64> {
    num.iter()
        .filter_map(|(r, a)| den.get(r).map(|b| a / b))
        .collect()
}

fn ratio_block(num: &BTreeMap<u64, f64>, den: &BTreeMap<u64, f64>) -> Value {
    let a: Vec<f64> = num.values().cloned().collect();
    let b: Vec<f64> = den.values().cloned().collect();
    let rr = round_ratios(num, den);
    let (lo, hi) = mm(&rr);
    json!({
        "ratio_of_medians": median(&a) / median(&b),
        "median_round_ratio": median(&rr),
        "round_ratio_min": lo,
        "round_ratio_max": hi,
        "rounds": rr.len(),
        "interval_kind": "min-max of rounds",
    })
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let run = Path::new(&a[0]);
    let raw = run.join("raw");
    let mut out = serde_json::Map::new();

    // ── BSGS wall time ──
    let rows = jsonl(&raw.join("bsgs_time.jsonl"));
    let mut by: BTreeMap<(u64, String), BTreeMap<u64, f64>> = BTreeMap::new();
    let mut digests: BTreeMap<u64, Vec<(String, String, u64)>> = BTreeMap::new();
    let mut loads: Vec<f64> = Vec::new();
    for r in &rows {
        let bits = r["bits"].as_u64().unwrap();
        let arm = r["arm"].as_str().unwrap().to_string();
        let round = r["round"].as_u64().unwrap();
        by.entry((bits, arm.clone()))
            .or_default()
            .insert(round, r["total_s"].as_f64().unwrap());
        digests.entry(bits).or_default().push((
            arm,
            r["answers_digest"].as_str().unwrap().to_string(),
            r["table_entries"].as_u64().unwrap(),
        ));
        if let Some(l) = r["uptime_load"].as_str() {
            if let Some(x) = l.split_whitespace().next() {
                loads.push(x.trim_end_matches(',').parse().unwrap_or(f64::NAN));
            }
        }
    }
    let mut bsgs = Vec::new();
    for (&bits, ds) in &digests {
        let same = ds.iter().all(|d| d.1 == ds[0].1 && d.2 == ds[0].2);
        let get = |arm: &str| {
            by.get(&(bits, arm.to_string()))
                .cloned()
                .unwrap_or_default()
        };
        let (o, p, pa, cz) = (get("original"), get("previous"), get("paired"), get("coz"));
        let arm_stats = |m: &BTreeMap<u64, f64>| {
            let v: Vec<f64> = m.values().cloned().collect();
            let (lo, hi) = mm(&v);
            json!({"median_s": median(&v), "min_s": lo, "max_s": hi, "samples": v.len()})
        };
        bsgs.push(json!({
            "bits": bits,
            "outputs_identical": same,
            "original": arm_stats(&o), "previous": arm_stats(&p),
            "paired": arm_stats(&pa), "coz": arm_stats(&cz),
            "original_over_paired": ratio_block(&o, &pa),
            "previous_over_paired": ratio_block(&p, &pa),
            "original_over_coz": ratio_block(&o, &cz),
            "previous_over_coz": ratio_block(&p, &cz),
            "original_over_previous": ratio_block(&o, &p),
        }));
    }
    if !bsgs.is_empty() {
        let (lo, hi) = mm(&loads);
        out.insert(
            "bsgs_time".into(),
            json!({"rows": bsgs, "load1_min": lo, "load1_max": hi}),
        );
    }

    // ── binary sect163 evaluators ──
    let rows = jsonl(&raw.join("binary_sect163.jsonl"));
    if !rows.is_empty() {
        let mut per: BTreeMap<u64, Value> = BTreeMap::new();
        for ell in [3u64, 5, 7] {
            let rs: Vec<&Value> = rows.iter().filter(|r| r["ell"] == ell).collect();
            if rs.is_empty() {
                continue;
            }
            let arms: Vec<String> = rs[0]["time_two_points"]
                .as_object()
                .unwrap()
                .keys()
                .cloned()
                .collect();
            let mut m = serde_json::Map::new();
            let base: Vec<f64> = rs
                .iter()
                .map(|r| {
                    r["time_two_points"]["velu_x_map_original"]["median_us"]
                        .as_f64()
                        .unwrap()
                })
                .collect();
            for arm in &arms {
                let v: Vec<f64> = rs
                    .iter()
                    .map(|r| r["time_two_points"][arm]["median_us"].as_f64().unwrap())
                    .collect();
                let ratios: Vec<f64> = base.iter().zip(&v).map(|(b, x)| b / x).collect();
                let (lo, hi) = mm(&ratios);
                m.insert(
                    arm.clone(),
                    json!({
                        "median_us_over_steps": median(&v),
                        "velu_x_map_over_this_median": median(&ratios),
                        "ratio_min": lo, "ratio_max": hi, "steps": v.len(),
                        "counts_two_points": rs[0]["counts_two_points"].get(arm).cloned().unwrap_or(Value::Null),
                    }),
                );
            }
            per.insert(ell, Value::Object(m));
        }
        out.insert("binary_sect163".into(), json!(per));
    }

    // ── challenge ladder: iterations / √(πn/4) per level ──
    let rows = jsonl(&raw.join("ladder_solve.jsonl"));
    if !rows.is_empty() {
        let mut by: BTreeMap<u64, Vec<f64>> = BTreeMap::new();
        let mut verified = true;
        for r in &rows {
            by.entry(r["bits"].as_u64().unwrap())
                .or_default()
                .push(r["ratio"].as_f64().unwrap());
            verified &= r["verified"].as_bool().unwrap_or(false);
        }
        let mut levels = Vec::new();
        let mut all = Vec::new();
        for (bits, v) in &by {
            let n = v.len() as f64;
            let mean = v.iter().sum::<f64>() / n;
            let sd = (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
            all.extend_from_slice(v);
            levels.push(
                json!({"bits": bits, "seeds": v.len(), "mean_ratio": mean, "sd": sd,
                "stderr": sd / n.sqrt(), "median_ratio": median(v)}),
            );
        }
        let n = all.len() as f64;
        let mean = all.iter().sum::<f64>() / n;
        let sd = (all.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
        out.insert(
            "challenge_ladder".into(),
            json!({"all_verified": verified, "levels": levels, "pooled": {"runs": all.len(),
                "mean_ratio": mean, "sd": sd, "stderr": sd / n.sqrt(), "median_ratio": median(&all)}}),
        );
    }

    // ── class walk ──
    if let Some(dir) = a.get(1) {
        let dir = Path::new(dir);
        let tsv = std::fs::read_to_string(dir.join("samples.tsv")).expect("samples.tsv");
        let mut wall: BTreeMap<String, BTreeMap<u64, f64>> = BTreeMap::new();
        let mut inner: BTreeMap<String, BTreeMap<u64, f64>> = BTreeMap::new();
        let mut snaps: Vec<(String, Value)> = Vec::new();
        for line in tsv.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            let (arm, round, secs) = (
                f[0].to_string(),
                f[1].parse::<u64>().unwrap(),
                f[2].parse::<f64>().unwrap(),
            );
            wall.entry(arm.clone()).or_default().insert(round, secs);
            let snap_path = dir
                .join(format!("wd-{arm}-{round}"))
                .join("experiments/koblitz_isogeny_class_walk.json");
            let mut snap: Value =
                serde_json::from_str(&std::fs::read_to_string(&snap_path).expect("snapshot"))
                    .unwrap();
            let s = snap["cases"][0]["seconds"].as_f64().unwrap();
            inner.entry(arm.clone()).or_default().insert(round, s);
            strip_seconds(&mut snap);
            snaps.push((arm, snap));
        }
        let identical = snaps.iter().all(|(_, s)| *s == snaps[0].1);
        let stats = |m: &BTreeMap<u64, f64>| {
            let v: Vec<f64> = m.values().cloned().collect();
            let (lo, hi) = mm(&v);
            json!({"median_s": median(&v), "min_s": lo, "max_s": hi, "samples": v.len()})
        };
        let g = |m: &BTreeMap<String, BTreeMap<u64, f64>>, k: &str| {
            m.get(k).cloned().unwrap_or_default()
        };
        out.insert(
            "classwalk_19_1".into(),
            json!({
                "snapshots_identical_without_seconds": identical,
                "process_wall": {
                    "original": stats(&g(&wall, "original")),
                    "previous": stats(&g(&wall, "previous")),
                    "candidate": stats(&g(&wall, "candidate")),
                    "original_over_candidate": ratio_block(&g(&wall, "original"), &g(&wall, "candidate")),
                    "previous_over_candidate": ratio_block(&g(&wall, "previous"), &g(&wall, "candidate")),
                },
                "walk_seconds_reported": {
                    "original": stats(&g(&inner, "original")),
                    "previous": stats(&g(&inner, "previous")),
                    "candidate": stats(&g(&inner, "candidate")),
                    "original_over_candidate": ratio_block(&g(&inner, "original"), &g(&inner, "candidate")),
                    "previous_over_candidate": ratio_block(&g(&inner, "previous"), &g(&inner, "candidate")),
                },
            }),
        );
    }

    let path = raw.join("summary.json");
    std::fs::write(
        &path,
        serde_json::to_string_pretty(&Value::Object(out)).unwrap(),
    )
    .unwrap();
    println!("wrote {}", path.display());
}
