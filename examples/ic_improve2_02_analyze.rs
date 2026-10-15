//! IC improvement round 2, item 02: analysis of the free-relation K sweep.
//!
//! Reads `data/runs/n<n>_K<K>_s<seed>_<arm>.{summary.json,jsonl}` written by
//! `research/ic_improve_20261015/02-freerel-production/scripts/run_job.sh`.
//!
//! Modes:
//! - `aa <runs> <traces> <out.json>`: A/A identity. For every `aa` run (new
//!   binary, flag unset) and its `incd` partner (frozen incumbent), the summary
//!   line, every target record and the rank trace must be byte-identical after
//!   deleting wall-clock, RSS and hardware-counter fields. Also records the
//!   blake3 hash of the base logs per (n, K) from the incumbent's rank trace.
//! - `train <runs> <scalars> <aa.json> <frozen_K.json>`: training table over
//!   seeds 101-104 and the frozen per-arm K* (argmin of mean counted cost).
//! - `test <runs> <scalars> <aa.json> <frozen_K.json> <results.json> <rows.csv>`:
//!   held-out paired comparisons with 95% percentile bootstrap intervals.
//!
//! Counted unit per run: index_pair_candidates + rank_root_calls_total +
//! mean target root_calls + harvest point additions (complete solve per target).
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_traits::ToPrimitive;
use serde_json::{json, Map, Value};
use std::collections::BTreeMap;
use std::io::Write;

const TARGETS: usize = 8;
const SIZES: [u32; 3] = [41, 53, 61];
const TRAIN_SEEDS: std::ops::RangeInclusive<u64> = 101..=104;
const K_DEF: [(u32, usize); 3] = [(41, 255), (53, 220), (61, 500)];

#[derive(Clone, Debug, Default)]
struct Rep {
    c_index: f64,
    c_rank: f64,
    c_target: f64,
    c_harvest: f64,
    c_total: f64,
    w_index: f64,
    w_harvest: f64,
    w_rank: f64,
    w_target: f64,
    w_total: f64,
    rank_attempts: f64,
    free_rank: f64,
    rss: f64,
}

fn num(v: &Value, path: &[&str]) -> f64 {
    let mut cur = v;
    for p in path {
        cur = &cur[*p];
    }
    cur.as_f64().unwrap_or(0.0)
}

fn read_lines(path: &str) -> Option<Vec<Value>> {
    let text = std::fs::read_to_string(path).ok()?;
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).ok())
        .collect()
}

/// Delete wall-clock, RSS and hardware-counter fields, recursively.
fn strip(v: &Value) -> Value {
    match v {
        Value::Object(map) => {
            let mut out = Map::new();
            for (k, x) in map {
                let timing = k.ends_with("_ms")
                    || k.ends_with("_ns")
                    || k == "peak_rss_bytes"
                    || k == "online_hw_counts";
                if !timing {
                    out.insert(k.clone(), strip(x));
                }
            }
            Value::Object(out)
        }
        Value::Array(xs) => Value::Array(xs.iter().map(strip).collect()),
        other => other.clone(),
    }
}

fn logs_hash(logs: &[u64]) -> String {
    let mut hasher = blake3::Hasher::new();
    for &log in logs {
        hasher.update(&log.to_le_bytes());
    }
    hasher.finalize().to_hex().to_string()
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

/// Ratio of summed means over groups; bootstrap resamples seeds within each group.
fn paired_ratio(groups: &[(Vec<f64>, Vec<f64>)], seed: u64) -> Value {
    let mean = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
    let point = groups.iter().map(|(c, _)| mean(c)).sum::<f64>()
        / groups.iter().map(|(_, b)| mean(b)).sum::<f64>();
    let mut rng = Lcg(seed);
    let mut draws = Vec::with_capacity(10_000);
    for _ in 0..10_000 {
        let (mut c, mut b) = (0.0, 0.0);
        for (cand, base) in groups {
            let k = cand.len();
            let (mut cs, mut bs) = (0.0, 0.0);
            for _ in 0..k {
                let i = (rng.next() % k as u64) as usize;
                cs += cand[i];
                bs += base[i];
            }
            c += cs / k as f64;
            b += bs / k as f64;
        }
        draws.push(c / b);
    }
    draws.sort_by(|x, y| x.partial_cmp(y).unwrap());
    json!({"ratio":point, "ci95":[draws[249], draws[9_749]],
           "pairs":groups.iter().map(|(c, _)| c.len()).sum::<usize>()})
}

fn side(r: &Value) -> i32 {
    let lo = r["ci95"][0].as_f64().unwrap();
    let hi = r["ci95"][1].as_f64().unwrap();
    if hi < 1.0 {
        -1
    } else if lo > 1.0 {
        1
    } else {
        0
    }
}

fn verdict(total: &Value, wall: &Value, phases: &[&Value]) -> String {
    let base = match side(total) {
        -1 => "ADVANCES",
        1 => "REGRESSES",
        _ => {
            let sides: Vec<i32> = phases.iter().map(|p| side(p)).collect();
            if sides.contains(&-1) && sides.contains(&1) {
                "TRADES"
            } else {
                "MATCHES"
            }
        }
    };
    let s = side(total);
    if s != 0 && side(wall) == -s {
        format!("{base} (contradicted by wall)")
    } else {
        base.to_owned()
    }
}

struct Ctx {
    runs: String,
    scalars: BTreeMap<u32, Vec<u64>>,
    order: BTreeMap<u32, u64>,
    logs: BTreeMap<String, String>,
}

impl Ctx {
    fn new(runs: &str, scalars_dir: &str, aa_path: Option<&str>) -> Self {
        let mut scalars = BTreeMap::new();
        let mut order = BTreeMap::new();
        for n in SIZES {
            let r = KoblitzCurve::new(0, n)
                .unwrap()
                .subgroup_order
                .to_u64()
                .unwrap();
            order.insert(n, r);
            let s: Vec<u64> = std::fs::read_to_string(format!("{scalars_dir}/scalars_n{n}.txt"))
                .unwrap()
                .lines()
                .take(TARGETS)
                .map(|l| l.trim().parse::<u64>().unwrap() % r)
                .collect();
            scalars.insert(n, s);
        }
        let mut logs = BTreeMap::new();
        if let Some(p) = aa_path {
            let aa: Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
            for (k, v) in aa["logs_hash"].as_object().unwrap() {
                logs.insert(k.clone(), v.as_str().unwrap().to_owned());
            }
        }
        Ctx {
            runs: runs.to_owned(),
            scalars,
            order,
            logs,
        }
    }

    /// Load and gate one run. Err(reason) on any correctness failure.
    fn load(&self, n: u32, k: usize, seed: u64, arm: &str) -> Option<Result<Rep, String>> {
        let stem = format!("{}/n{n}_K{k}_s{seed}_{arm}", self.runs);
        let summary_text = std::fs::read_to_string(format!("{stem}.summary.json")).ok()?;
        let Ok(summary) = serde_json::from_str::<Value>(summary_text.trim()) else {
            return Some(Err(format!("{stem}: no summary (killed?)")));
        };
        let targets = read_lines(&format!("{stem}.jsonl")).unwrap_or_default();
        if targets.len() != TARGETS {
            return Some(Err(format!("{stem}: {} target rows", targets.len())));
        }
        for (i, t) in targets.iter().enumerate() {
            let verified = t["group_verified"].as_bool() == Some(true);
            let recovered = t["recovered_scalar"].as_u64();
            if !verified || recovered != Some(self.scalars[&n][i]) {
                return Some(Err(format!(
                    "{stem}: target {i} verified={verified} recovered={recovered:?}"
                )));
            }
        }
        if summary["rank"].as_u64() != Some(k as u64) {
            return Some(Err(format!("{stem}: rank {}", summary["rank"])));
        }
        let fr = &summary["free_relations"];
        if arm == "free" {
            if fr["enabled"].as_bool() != Some(true) || fr["failed"].as_u64() != Some(0) {
                return Some(Err(format!("{stem}: harvest failed or flag off")));
            }
            match self.logs.get(&format!("n{n}_K{k}")) {
                Some(h) if fr["logs_hash"].as_str() == Some(h.as_str()) => {}
                Some(_) => return Some(Err(format!("{stem}: base logs differ from incumbent"))),
                None => return Some(Err(format!("{stem}: no incumbent logs for n{n}_K{k}"))),
            }
        } else if !fr.is_null() {
            return Some(Err(format!("{stem}: free_relations present on {arm}")));
        }
        let _ = self.order[&n];
        let mean_t =
            |key: &str| targets.iter().map(|t| num(t, &[key])).sum::<f64>() / TARGETS as f64;
        let c_index = num(&summary, &["index_pair_candidates"]);
        let c_rank = num(&summary, &["rank_root_calls_total"]);
        let c_target = mean_t("root_calls");
        let c_harvest = num(fr, &["harvest_additions"]);
        let tm = &summary["timing_ms"];
        let w_setup = num(tm, &["setup_total"]);
        let w_rank = num(tm, &["rank_stage"]);
        let w_la = num(tm, &["linear_algebra"]);
        let w_target = mean_t("online_ms");
        Some(Ok(Rep {
            c_index,
            c_rank,
            c_target,
            c_harvest,
            c_total: c_index + c_rank + c_target + c_harvest,
            w_index: num(tm, &["index_build"]),
            w_harvest: num(fr, &["harvest_ms"]),
            w_rank,
            w_target,
            w_total: w_setup + w_rank + w_la + w_target,
            rank_attempts: num(&summary, &["rank_attempts"]),
            free_rank: num(fr, &["free_rank"]),
            rss: num(&summary, &["peak_rss_bytes"]),
        }))
    }
}

fn mode_aa(runs: &str, traces: &str, out: &str) {
    let mut pairs = Vec::new();
    let mut logs = Map::new();
    let mut failures = Vec::new();
    let mut entries: Vec<String> = std::fs::read_dir(runs)
        .unwrap()
        .filter_map(|e| e.ok()?.file_name().into_string().ok())
        .filter(|f| f.ends_with("_aa.summary.json"))
        .collect();
    entries.sort();
    for f in entries {
        let stem = f.trim_end_matches("_aa.summary.json").to_owned();
        let mut identical = true;
        let mut detail = Vec::new();
        for (suffix, a_path, b_path) in [
            (
                "summary",
                format!("{runs}/{stem}_aa.summary.json"),
                format!("{runs}/{stem}_incd.summary.json"),
            ),
            (
                "targets",
                format!("{runs}/{stem}_aa.jsonl"),
                format!("{runs}/{stem}_incd.jsonl"),
            ),
            (
                "rank_trace",
                format!("{traces}/{stem}_aa.rank.jsonl"),
                format!("{traces}/{stem}_incd.rank.jsonl"),
            ),
        ] {
            let (Some(a), Some(b)) = (read_lines(&a_path), read_lines(&b_path)) else {
                identical = false;
                detail.push(format!("{suffix}: missing"));
                continue;
            };
            let sa: Vec<String> = a.iter().map(|v| strip(v).to_string()).collect();
            let sb: Vec<String> = b.iter().map(|v| strip(v).to_string()).collect();
            let same = sa == sb && !sa.is_empty();
            if !same {
                identical = false;
                detail.push(format!("{suffix}: differs"));
            }
            if suffix == "rank_trace" {
                if let Some(sol) = b
                    .iter()
                    .find(|v| v["kind"] == "compact_orbit_rank_solution")
                {
                    let l: Vec<u64> = serde_json::from_value(sol["logs"].clone()).unwrap();
                    let key: String = stem.split("_s").next().unwrap().to_owned();
                    let h = logs_hash(&l);
                    if let Some(prev) = logs.get(&key) {
                        if prev.as_str() != Some(h.as_str()) {
                            failures.push(format!("{stem}: logs differ across seeds"));
                        }
                    }
                    logs.insert(key, json!(h));
                }
                detail.push(format!("{suffix}: {} lines", a.len()));
            } else {
                detail.push(format!("{suffix}: {} lines", a.len()));
            }
        }
        if !identical {
            failures.push(format!("{stem}: {}", detail.join("; ")));
        }
        pairs.push(json!({"run":stem, "identical":identical, "detail":detail}));
    }
    let result = json!({
        "kind":"ic_improve2_02_aa",
        "rule":"summary, target records and rank trace byte-identical after deleting *_ms, *_ns, peak_rss_bytes, online_hw_counts",
        "pairs_checked":pairs.len(),
        "pairs_identical":pairs.iter().filter(|p| p["identical"] == true).count(),
        "failures":failures,
        "pairs":pairs,
        "logs_hash":logs,
    });
    std::fs::write(
        out,
        format!("{}\n", serde_json::to_string_pretty(&result).unwrap()),
    )
    .unwrap();
    println!(
        "A/A: {}/{} identical; failures: {:?}",
        result["pairs_identical"], result["pairs_checked"], result["failures"]
    );
}

fn ks_in(runs: &str, n: u32, arm: &str, seeds: &[u64]) -> Vec<usize> {
    let mut ks: Vec<usize> = std::fs::read_dir(runs)
        .unwrap()
        .filter_map(|e| e.ok()?.file_name().into_string().ok())
        .filter_map(|f| {
            let rest = f.strip_prefix(&format!("n{n}_K"))?;
            let (k, rest) = rest.split_once("_s")?;
            let (seed, rest) = rest.split_once('_')?;
            (rest == format!("{arm}.summary.json") && seeds.contains(&seed.parse::<u64>().ok()?))
                .then(|| k.parse().ok())?
        })
        .collect();
    ks.sort();
    ks.dedup();
    ks
}

fn rep_json(reps: &[Rep]) -> Value {
    let m = |f: fn(&Rep) -> f64| reps.iter().map(f).sum::<f64>() / reps.len() as f64;
    let sd = |f: fn(&Rep) -> f64| {
        let mu = m(f);
        (reps.iter().map(|r| (f(r) - mu).powi(2)).sum::<f64>() / (reps.len() as f64 - 1.0).max(1.0))
            .sqrt()
    };
    json!({
        "seeds":reps.len(),
        "c_total_mean":m(|r| r.c_total), "c_total_sd":sd(|r| r.c_total),
        "c_index":m(|r| r.c_index), "c_rank_mean":m(|r| r.c_rank),
        "c_target_mean":m(|r| r.c_target), "c_harvest_mean":m(|r| r.c_harvest),
        "w_total_mean_ms":m(|r| r.w_total), "w_total_sd_ms":sd(|r| r.w_total),
        "w_index_mean_ms":m(|r| r.w_index), "w_rank_mean_ms":m(|r| r.w_rank),
        "w_harvest_mean_ms":m(|r| r.w_harvest), "w_target_mean_ms":m(|r| r.w_target),
        "rank_attempts_mean":m(|r| r.rank_attempts), "free_rank_mean":m(|r| r.free_rank),
        "peak_rss_max_bytes":reps.iter().map(|r| r.rss).fold(0.0, f64::max),
    })
}

fn mode_train(ctx: &Ctx, out: &str) {
    let seeds: Vec<u64> = TRAIN_SEEDS.collect();
    let mut table = Map::new();
    let mut frozen = Map::new();
    let mut failures = Vec::new();
    for n in SIZES {
        let mut size_out = Map::new();
        let mut best: BTreeMap<&str, (usize, f64)> = BTreeMap::new();
        for arm in ["inc", "free"] {
            let mut arm_out = Map::new();
            for k in ks_in(&ctx.runs, n, arm, &seeds) {
                let mut reps = Vec::new();
                let mut eligible = true;
                for &seed in &seeds {
                    match ctx.load(n, k, seed, arm) {
                        Some(Ok(r)) => reps.push(r),
                        Some(Err(e)) => {
                            eligible = false;
                            failures.push(e);
                        }
                        None => eligible = false,
                    }
                }
                let mut v = if reps.is_empty() {
                    json!({})
                } else {
                    rep_json(&reps)
                };
                v["eligible"] = json!(eligible);
                if eligible {
                    let c = v["c_total_mean"].as_f64().unwrap();
                    if best.get(arm).is_none_or(|&(_, b)| c < b) {
                        best.insert(arm, (k, c));
                    }
                }
                arm_out.insert(k.to_string(), v);
            }
            size_out.insert(arm.to_owned(), Value::Object(arm_out));
        }
        table.insert(format!("n{n}"), Value::Object(size_out));
        frozen.insert(
            format!("n{n}"),
            json!({
                "K_inc":best.get("inc").map(|b| b.0), "K_free":best.get("free").map(|b| b.0),
                "C_inc_train":best.get("inc").map(|b| b.1), "C_free_train":best.get("free").map(|b| b.1),
                "K_def":K_DEF.iter().find(|(m, _)| *m == n).unwrap().1,
            }),
        );
    }
    let result = json!({
        "kind":"ic_improve2_02_frozen_K",
        "rule":"per (n, arm) K* = argmin over the grid of the mean counted cost over training seeds 101-104; a K with any failed run is ineligible",
        "frozen":frozen, "failures":failures, "training":table,
    });
    std::fs::write(
        out,
        format!("{}\n", serde_json::to_string_pretty(&result).unwrap()),
    )
    .unwrap();
    println!(
        "{}",
        serde_json::to_string_pretty(&result["frozen"]).unwrap()
    );
    println!("failures: {:?}", result["failures"]);
}

fn mode_test(ctx: &Ctx, frozen_path: &str, out: &str, csv_path: &str) {
    let frozen: Value =
        serde_json::from_str(&std::fs::read_to_string(frozen_path).unwrap()).unwrap();
    let mut csv = std::fs::File::create(csv_path).unwrap();
    writeln!(csv, "n,K,seed,arm,split,c_index,c_rank,c_target,c_harvest,c_total,w_index,w_harvest,w_rank,w_target,w_total,rank_attempts,free_rank,peak_rss_bytes").unwrap();
    let mut failures = Vec::new();
    // Every run, training and test, into rows.csv.
    let mut all: BTreeMap<(u32, usize, String), BTreeMap<u64, Rep>> = BTreeMap::new();
    for n in SIZES {
        for arm in ["inc", "free"] {
            let seeds: Vec<u64> = (1..=200).collect();
            for k in ks_in(&ctx.runs, n, arm, &seeds) {
                for &seed in &seeds {
                    match ctx.load(n, k, seed, arm) {
                        Some(Ok(r)) => {
                            let split = if TRAIN_SEEDS.contains(&seed) {
                                "train"
                            } else {
                                "test"
                            };
                            writeln!(csv, "{n},{k},{seed},{arm},{split},{},{},{},{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{},{},{}",
                                r.c_index, r.c_rank, r.c_target, r.c_harvest, r.c_total,
                                r.w_index, r.w_harvest, r.w_rank, r.w_target, r.w_total,
                                r.rank_attempts, r.free_rank, r.rss).unwrap();
                            all.entry((n, k, arm.to_owned()))
                                .or_default()
                                .insert(seed, r);
                        }
                        Some(Err(e)) => failures.push(e),
                        None => {}
                    }
                }
            }
        }
    }
    let test_seeds = |n: u32| -> Vec<u64> {
        if n == 61 {
            (1..=12).collect()
        } else {
            (1..=16).collect()
        }
    };
    // A comparison: candidate (arm, K) over base (arm, K), paired on test seeds.
    type Pick = fn(&Rep) -> f64;
    let metrics: [(&str, Pick); 6] = [
        ("counted_total", |r| r.c_total),
        ("counted_index", |r| r.c_index),
        ("counted_rank", |r| r.c_rank),
        ("counted_target", |r| r.c_target),
        ("wall_total", |r| r.w_total),
        ("wall_rank", |r| r.w_rank),
    ];
    let pairs_for = |n: u32, ck: usize, ca: &str, bk: usize, ba: &str| -> Vec<(Rep, Rep)> {
        let empty = BTreeMap::new();
        let c = all.get(&(n, ck, ca.to_owned())).unwrap_or(&empty);
        let b = all.get(&(n, bk, ba.to_owned())).unwrap_or(&empty);
        test_seeds(n)
            .into_iter()
            .filter_map(|s| Some((c.get(&s)?.clone(), b.get(&s)?.clone())))
            .collect()
    };
    let compare = |spec: &[(u32, usize, &str, usize, &str)], seed: u64| -> Value {
        let mut o = Map::new();
        for (name, f) in metrics {
            let groups: Vec<(Vec<f64>, Vec<f64>)> = spec
                .iter()
                .map(|&(n, ck, ca, bk, ba)| {
                    let p = pairs_for(n, ck, ca, bk, ba);
                    (
                        p.iter().map(|(c, _)| f(c)).collect(),
                        p.iter().map(|(_, b)| f(b)).collect(),
                    )
                })
                .collect();
            if groups.iter().any(|(c, _)| c.is_empty()) {
                return json!({"missing":true});
            }
            o.insert(name.to_owned(), paired_ratio(&groups, seed));
        }
        let v = verdict(
            &o["counted_total"],
            &o["wall_total"],
            &[&o["counted_index"], &o["counted_rank"]],
        );
        o.insert("verdict".to_owned(), json!(v));
        // Means on the paired seeds, for the tables.
        if spec.len() == 1 {
            let (n, ck, ca, bk, ba) = spec[0];
            let p = pairs_for(n, ck, ca, bk, ba);
            let c: Vec<Rep> = p.iter().map(|x| x.0.clone()).collect();
            let b: Vec<Rep> = p.iter().map(|x| x.1.clone()).collect();
            o.insert(
                "candidate".to_owned(),
                json!({"arm":ca, "K":ck, "means":rep_json(&c)}),
            );
            o.insert(
                "base".to_owned(),
                json!({"arm":ba, "K":bk, "means":rep_json(&b)}),
            );
        }
        Value::Object(o)
    };
    let fk = |n: u32, key: &str| frozen["frozen"][format!("n{n}")][key].as_u64().unwrap() as usize;
    let mut sizes = Map::new();
    let (mut pooled_b, mut pooled_c) = (Vec::new(), Vec::new());
    let mut pooled_a: BTreeMap<&str, Vec<(u32, usize, &str, usize, &str)>> = BTreeMap::new();
    for n in SIZES {
        let (ki, kf, kd) = (fk(n, "K_inc"), fk(n, "K_free"), fk(n, "K_def"));
        let mut s = Map::new();
        s.insert("K".to_owned(), json!({"K_inc":ki, "K_free":kf, "K_def":kd}));
        let mut same = Map::new();
        for (label, k) in [("K_inc", ki), ("K_free", kf), ("K_def", kd)] {
            same.insert(
                format!("{label}={k}"),
                compare(&[(n, k, "free", k, "inc")], 7 + n as u64),
            );
            pooled_a
                .entry(label)
                .or_default()
                .push((n, k, "free", k, "inc"));
        }
        s.insert("a_same_K".to_owned(), Value::Object(same));
        s.insert(
            "b_tuned_vs_tuned".to_owned(),
            compare(&[(n, kf, "free", ki, "inc")], 11 + n as u64),
        );
        s.insert(
            "c_tuned_free_vs_inc_K_def".to_owned(),
            compare(&[(n, kf, "free", kd, "inc")], 13 + n as u64),
        );
        s.insert(
            "inc_tuned_vs_inc_K_def".to_owned(),
            compare(&[(n, ki, "inc", kd, "inc")], 17 + n as u64),
        );
        pooled_b.push((n, kf, "free", ki, "inc"));
        pooled_c.push((n, kf, "free", kd, "inc"));
        sizes.insert(format!("n{n}"), Value::Object(s));
    }
    let mut pooled = Map::new();
    pooled.insert("b_tuned_vs_tuned".to_owned(), compare(&pooled_b, 101));
    pooled.insert(
        "c_tuned_free_vs_inc_K_def".to_owned(),
        compare(&pooled_c, 103),
    );
    for (label, spec) in &pooled_a {
        pooled.insert(format!("a_same_K_{label}"), compare(spec, 107));
    }
    let result = json!({
        "kind":"ic_improve2_02_results",
        "statistic":"ratio of means over held-out seeds (pooled: ratio of summed per-size means), 95% percentile bootstrap over seeds within size, 10,000 resamples",
        "counted_unit":"index_pair_candidates + rank_root_calls_total + mean target root_calls + harvest point additions",
        "frozen":frozen["frozen"], "sizes":sizes, "pooled":pooled,
        "correctness_failures":failures,
    });
    std::fs::write(
        out,
        format!("{}\n", serde_json::to_string_pretty(&result).unwrap()),
    )
    .unwrap();
    for n in SIZES {
        let s = &result["sizes"][format!("n{n}")];
        let b = &s["b_tuned_vs_tuned"];
        println!(
            "n{n} K={} (b) C {} wall {} -> {}",
            s["K"], b["counted_total"], b["wall_total"], b["verdict"]
        );
    }
    let b = &result["pooled"]["b_tuned_vs_tuned"];
    println!(
        "pooled (b) C {} wall {} -> {}",
        b["counted_total"], b["wall_total"], b["verdict"]
    );
    println!("failures: {:?}", result["correctness_failures"]);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("aa") if args.len() == 5 => mode_aa(&args[2], &args[3], &args[4]),
        Some("train") if args.len() == 6 => {
            mode_train(&Ctx::new(&args[2], &args[3], Some(&args[4])), &args[5])
        }
        Some("test") if args.len() == 8 => mode_test(
            &Ctx::new(&args[2], &args[3], Some(&args[4])),
            &args[5],
            &args[6],
            &args[7],
        ),
        _ => panic!(
            "usage: aa <runs> <traces> <out> | train <runs> <scalars> <aa.json> <frozen.json> | test <runs> <scalars> <aa.json> <frozen.json> <results.json> <rows.csv>"
        ),
    }
}
