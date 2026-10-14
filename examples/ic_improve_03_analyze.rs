//! IC improvement item 03: analysis of the paired arm runs.
//!
//! Reads `data/n<n>_K<K>_s<seed>_<arm>.{summary.json,jsonl}` written by
//! `scripts/run_arms.sh`, gates correctness (every target verified, recovered
//! scalar equal to item 09's held-back scalar, base logs identical across arms
//! and seeds), and reports per-replicate costs and paired ratios against arm
//! `A` with 95% percentile bootstrap intervals over seeds (ratio of means).
//!
//! Usage: `ic_improve_03_analyze <data_dir> <scalars_dir> <out_results.json> <out_rows.csv>`
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_traits::ToPrimitive;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::io::Write;

const SIZES: [(u32, usize); 3] = [(41, 255), (53, 220), (61, 500)];
const ARMS: [&str; 5] = ["inc", "A", "F", "R", "FR"];
const TARGETS: usize = 8;

#[derive(Clone, Debug, Default)]
struct Rep {
    // counted unit (S3 root-solve equivalents), mean over targets
    c_index: f64,
    c_rank: f64,
    c_target: f64,
    c_filter: f64,
    c_harvest: f64,
    c_total: f64,
    // wall ms
    w_setup_other: f64,
    w_index: f64,
    w_harvest: f64,
    w_rank: f64,
    w_la: f64,
    w_target: f64,
    w_total: f64,
    rank_attempts: f64,
    free_rank: f64,
    filter_rejections: f64,
}

fn get(v: &Value, path: &[&str]) -> f64 {
    let mut cur = v;
    for p in path {
        cur = &cur[*p];
    }
    cur.as_f64().unwrap_or(0.0)
}

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

/// Ratio of means sum(cand)/sum(base) and its 95% percentile bootstrap over pairs.
fn paired_ratio(cand: &[f64], base: &[f64], seed: u64) -> (f64, f64, f64) {
    let k = cand.len();
    let point = cand.iter().sum::<f64>() / base.iter().sum::<f64>();
    let mut rng = Lcg(seed);
    let mut draws = Vec::with_capacity(10_000);
    for _ in 0..10_000 {
        let (mut c, mut b) = (0.0, 0.0);
        for _ in 0..k {
            let i = (rng.next() % k as u64) as usize;
            c += cand[i];
            b += base[i];
        }
        draws.push(c / b);
    }
    draws.sort_by(|x, y| x.partial_cmp(y).unwrap());
    (point, draws[249], draws[9_749])
}

fn verdict(lo: f64, hi: f64) -> &'static str {
    if hi < 1.0 {
        "ADVANCES"
    } else if lo > 1.0 {
        "REGRESSES"
    } else {
        "MATCHES"
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 5, "usage: <data_dir> <scalars_dir> <results.json> <rows.csv>");
    let data = &args[1];
    let scalars_dir = &args[2];
    let mut csv = std::fs::File::create(&args[4]).unwrap();
    writeln!(csv, "n,K,seed,arm,c_index,c_rank,c_target,c_filter,c_harvest,c_total,w_setup_other,w_index,w_harvest,w_rank,w_la,w_target,w_total,rank_attempts,free_rank,filter_rejections").unwrap();
    let mut results = serde_json::Map::new();
    let mut failures: Vec<String> = Vec::new();
    let mut pooled: BTreeMap<(&str, &str), (Vec<f64>, Vec<f64>)> = BTreeMap::new();
    for (n, k) in SIZES {
        let r = KoblitzCurve::new(0, n)
            .unwrap()
            .subgroup_order
            .to_u64()
            .unwrap();
        let scalars: Vec<u64> = std::fs::read_to_string(format!("{scalars_dir}/scalars_n{n}.txt"))
            .unwrap()
            .lines()
            .take(TARGETS)
            .map(|l| l.trim().parse().unwrap())
            .collect();
        let mut reps: BTreeMap<&str, BTreeMap<u64, Rep>> = BTreeMap::new();
        let mut logs_hashes: BTreeMap<String, usize> = BTreeMap::new();
        for seed in 1..=64u64 {
            for arm in ARMS {
                let stem = format!("{data}/n{n}_K{k}_s{seed}_{arm}");
                let Ok(summary_text) = std::fs::read_to_string(format!("{stem}.summary.json")) else {
                    continue;
                };
                let Ok(summary) = serde_json::from_str::<Value>(summary_text.trim()) else {
                    failures.push(format!("{stem}: no summary"));
                    continue;
                };
                let targets: Vec<Value> = std::fs::read_to_string(format!("{stem}.jsonl"))
                    .unwrap()
                    .lines()
                    .map(|l| serde_json::from_str(l).unwrap())
                    .collect();
                if targets.len() != TARGETS {
                    failures.push(format!("{stem}: {} target rows", targets.len()));
                    continue;
                }
                let mut ok = true;
                for (i, t) in targets.iter().enumerate() {
                    let verified = t["group_verified"].as_bool() == Some(true);
                    let recovered = t["recovered_scalar"].as_u64();
                    if !verified || recovered != Some(scalars[i] % r) {
                        ok = false;
                        failures.push(format!("{stem}: target {i} verified={verified} recovered={recovered:?}"));
                    }
                }
                if !ok {
                    continue;
                }
                if let Some(h) = summary["item03"]["logs_hash"].as_str() {
                    *logs_hashes.entry(h.to_owned()).or_default() += 1;
                }
                let t_root: f64 = targets.iter().map(|t| get(t, &["root_calls"])).sum::<f64>()
                    / TARGETS as f64;
                let t_wall: f64 =
                    targets.iter().map(|t| get(t, &["online_ms"])).sum::<f64>() / TARGETS as f64;
                let filter_checks = get(&summary, &["item03", "filter_checks"]);
                // Filter checks: rank-stage ones are setup, one per target is online.
                let per_target_filter = if filter_checks > 0.0 { 1.0 } else { 0.0 };
                let c_filter = filter_checks - per_target_filter * TARGETS as f64 + per_target_filter;
                let c_index = get(&summary, &["index_pair_candidates"]);
                let c_rank = get(&summary, &["rank_root_calls_total"]);
                let c_harvest = get(&summary, &["item03", "harvest_additions"]);
                let tm = &summary["timing_ms"];
                let w_index = get(tm, &["index_build"]);
                let w_harvest = get(&summary, &["item03", "harvest_ms"]);
                let w_setup_total = get(tm, &["setup_total"]);
                let w_rank = get(tm, &["rank_stage"]);
                let w_la = get(tm, &["linear_algebra"]);
                let rep = Rep {
                    c_index,
                    c_rank,
                    c_target: t_root,
                    c_filter,
                    c_harvest,
                    c_total: c_index + c_rank + t_root + c_filter + c_harvest,
                    w_setup_other: w_setup_total - w_index - w_harvest,
                    w_index,
                    w_harvest,
                    w_rank,
                    w_la,
                    w_target: t_wall,
                    w_total: w_setup_total + w_rank + w_la + t_wall,
                    rank_attempts: get(&summary, &["rank_attempts"]),
                    free_rank: get(&summary, &["item03", "free_rank"]),
                    filter_rejections: get(&summary, &["item03", "filter_rejections"]),
                };
                writeln!(
                    csv,
                    "{n},{k},{seed},{arm},{},{},{},{},{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{},{},{}",
                    rep.c_index, rep.c_rank, rep.c_target, rep.c_filter, rep.c_harvest, rep.c_total,
                    rep.w_setup_other, rep.w_index, rep.w_harvest, rep.w_rank, rep.w_la, rep.w_target,
                    rep.w_total, rep.rank_attempts, rep.free_rank, rep.filter_rejections
                )
                .unwrap();
                reps.entry(arm).or_default().insert(seed, rep);
            }
        }
        if logs_hashes.len() > 1 {
            failures.push(format!("n={n}: {} distinct base-log hashes", logs_hashes.len()));
        }
        // Paired comparisons against A (and inc against A as the A/A null).
        let base = reps.get("A").cloned().unwrap_or_default();
        let mut size_out = serde_json::Map::new();
        for arm in ["inc", "F", "R", "FR"] {
            let Some(cand) = reps.get(arm) else { continue };
            let seeds: Vec<u64> = cand.keys().filter(|s| base.contains_key(s)).copied().collect();
            if seeds.is_empty() {
                continue;
            }
            let pick = |f: fn(&Rep) -> f64, m: &BTreeMap<u64, Rep>| -> Vec<f64> {
                seeds.iter().map(|s| f(&m[s])).collect()
            };
            let metrics: [(&str, fn(&Rep) -> f64); 12] = [
                ("counted_total", |r| r.c_total),
                ("counted_index", |r| r.c_index),
                ("counted_rank", |r| r.c_rank),
                ("counted_target", |r| r.c_target),
                ("counted_filter_plus_harvest", |r| r.c_filter + r.c_harvest),
                ("wall_total", |r| r.w_total),
                ("wall_index", |r| r.w_index),
                ("wall_harvest_plus_other_setup", |r| r.w_harvest + r.w_setup_other),
                ("wall_rank", |r| r.w_rank),
                ("wall_la", |r| r.w_la),
                ("wall_target", |r| r.w_target),
                ("rank_attempts", |r| r.rank_attempts),
            ];
            let mut arm_out = serde_json::Map::new();
            for (i, (name, f)) in metrics.iter().enumerate() {
                let c = pick(*f, cand);
                let b = pick(*f, &base);
                let (point, lo, hi) = if b.iter().sum::<f64>() > 0.0 {
                    paired_ratio(&c, &b, 0xC0FFEE ^ (n as u64 * 1000 + i as u64))
                } else {
                    (f64::NAN, f64::NAN, f64::NAN)
                };
                arm_out.insert(
                    (*name).into(),
                    json!({"ratio":point, "ci95":[lo, hi],
                           "mean_candidate":c.iter().sum::<f64>()/c.len() as f64,
                           "mean_A":b.iter().sum::<f64>()/b.len() as f64}),
                );
                if *name == "counted_total" || *name == "wall_total" {
                    let e = pooled.entry((arm, name)).or_default();
                    // Pool as per-size normalised values (divide by the size's A mean).
                    let bm = b.iter().sum::<f64>() / b.len() as f64;
                    e.0.extend(c.iter().map(|x| x / bm));
                    e.1.extend(b.iter().map(|x| x / bm));
                }
            }
            let ct = &arm_out["counted_total"]["ci95"];
            arm_out.insert(
                "verdict_counted_total".into(),
                json!(verdict(ct[0].as_f64().unwrap(), ct[1].as_f64().unwrap())),
            );
            arm_out.insert("replicates".into(), json!(seeds.len()));
            let fr: Vec<f64> = pick(|r| r.free_rank, cand);
            let rej: Vec<f64> = pick(|r| r.filter_rejections, cand);
            arm_out.insert("free_rank_mean".into(), json!(fr.iter().sum::<f64>() / fr.len() as f64));
            arm_out.insert("filter_rejections_total".into(), json!(rej.iter().sum::<f64>()));
            arm_out.insert(
                "expected_rank_ratio_K_minus_F_over_K".into(),
                json!((k as f64 - fr.iter().sum::<f64>() / fr.len() as f64) / k as f64),
            );
            size_out.insert(arm.into(), Value::Object(arm_out));
        }
        size_out.insert("K".into(), json!(k));
        size_out.insert("r".into(), json!(r));
        size_out.insert("A_replicates".into(), json!(base.len()));
        size_out.insert("distinct_logs_hashes".into(), json!(logs_hashes.len()));
        results.insert(format!("n{n}"), Value::Object(size_out));
    }
    let mut pooled_out = serde_json::Map::new();
    for ((arm, name), (c, b)) in &pooled {
        let (point, lo, hi) = paired_ratio(c, b, 0xBEEF);
        pooled_out.insert(
            format!("{arm}.{name}"),
            json!({"ratio":point, "ci95":[lo, hi], "pairs":c.len(), "verdict":verdict(lo, hi)}),
        );
    }
    results.insert("pooled_size_normalised".into(), Value::Object(pooled_out));
    results.insert("correctness_failures".into(), json!(failures));
    let out = Value::Object(results);
    std::fs::write(&args[3], serde_json::to_string_pretty(&out).unwrap()).unwrap();
    println!("{}", serde_json::to_string_pretty(&out).unwrap());
}
