//! Analysis for `research/ic_improve_20261014/04b-genus3-2lp-rebalance/` (PROTOCOL.md §4, §6–§8).
//!
//!   tune   --dir RESULTS_TUNING --out tuning/kappa_choice.json
//!   report --dir RESULTS_TEST --rho RESULTS_RHO --choice tuning/kappa_choice.json --out analysis.json --md tables.md
//!
//! Reads only the tool's own JSON-lines records. Computes everything the report states;
//! the figure script only draws `analysis.json`.
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write as _;

const GRID: [f64; 5] = [0.5, 0.7, 1.0, 1.4, 2.0];
const BOOT_B: usize = 10000;
const BOOT_SEED: u64 = 20261016;

fn splitmix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        splitmix(self.0)
    }
    fn below(&mut self, n: usize) -> usize {
        ((self.next() as u128 * n as u128) >> 64) as usize
    }
}

#[derive(Clone, Debug)]
struct La {
    total: f64,
    la: f64,
    lanczos: f64,
    rows: f64,
    cols: f64,
    nnz: f64,
    w_row: f64,
    w_col: f64,
    iters: u64,
    attempts: u64,
    verified: bool,
    planted: bool,
}
#[derive(Clone, Debug)]
struct Run {
    k: u32,
    q: u64,
    curve: u32,
    target: u32,
    kappa: f64,
    ell: u64,
    rel: f64,
    las: HashMap<String, La>,
    rss: f64,
    wall: f64,
}

fn read_records(dir: &str) -> Vec<Value> {
    let mut out = vec![];
    let mut stack = vec![std::path::PathBuf::from(dir)];
    while let Some(p) = stack.pop() {
        if p.is_dir() {
            for e in std::fs::read_dir(&p).unwrap() {
                stack.push(e.unwrap().path());
            }
        } else if p.extension().map(|e| e == "jsonl").unwrap_or(false)
            && !p
                .file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with("failures")
        {
            for line in std::fs::read_to_string(&p).unwrap().lines() {
                if line.trim().is_empty() {
                    continue;
                }
                out.push(
                    serde_json::from_str::<Value>(line)
                        .unwrap_or_else(|e| panic!("{}: {}", p.display(), e)),
                );
            }
        }
    }
    out
}
fn f(v: &Value, k: &str) -> f64 {
    v[k].as_f64().unwrap_or_else(|| {
        v[k].as_str()
            .map(|s| s.parse().unwrap())
            .unwrap_or(f64::NAN)
    })
}
fn u(v: &Value, k: &str) -> u64 {
    v[k].as_u64().unwrap_or_else(|| panic!("field {}", k))
}
fn runs(dir: &str) -> Vec<Run> {
    read_records(dir)
        .into_iter()
        .filter(|v| v["cmd"] == "solve")
        .map(|v| {
            let mut las = HashMap::new();
            for l in v["la"].as_array().unwrap() {
                las.insert(
                    l["method"].as_str().unwrap().to_string(),
                    La {
                        total: f(l, "total_mults"),
                        la: f(l, "elim_mults") + f(l, "lanczos_mults"),
                        lanczos: f(l, "lanczos_mults"),
                        rows: f(l, "rows"),
                        cols: f(l, "cols"),
                        nnz: f(l, "nnz"),
                        w_row: f(l, "w_row"),
                        w_col: f(l, "w_col"),
                        iters: u(l, "lanczos_iters"),
                        attempts: u(l, "attempts"),
                        verified: l["log_verified"].as_bool().unwrap(),
                        planted: l["log_equals_planted"].as_bool().unwrap(),
                    },
                );
            }
            Run {
                k: u(&v, "k") as u32,
                q: u(&v, "q"),
                curve: u(&v, "curve") as u32,
                target: u(&v, "target") as u32,
                kappa: f(&v, "kappa"),
                ell: u(&v, "ell"),
                rel: f(&v["rel"], "rel_mults"),
                las,
                rss: f(&v, "rss_mb"),
                wall: f(&v, "wall_s"),
            }
        })
        .collect()
}
fn kap_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}
fn ok(r: &Run, m: &str) -> bool {
    r.las
        .get(m)
        .map(|l| l.verified && l.planted)
        .unwrap_or(false)
}

fn cmd_tune(dir: &str, out: &str) {
    let rs = runs(dir);
    let mut choice = serde_json::Map::new();
    let mut md = String::new();
    for k in 14..=20u32 {
        let mut table = vec![];
        let (mut best_b, mut best_s) = ((f64::INFINITY, 0.0), (f64::INFINITY, 0.0));
        for &kap in GRID.iter() {
            let sel: Vec<&Run> = rs
                .iter()
                .filter(|r| {
                    r.k == k
                        && kap_eq(r.kappa, kap)
                        && (5..=7).contains(&r.curve)
                        && (1..=2).contains(&r.target)
                })
                .collect();
            let complete = sel.len() == 6
                && sel
                    .iter()
                    .all(|r| ok(r, "bfs_lanczos") && ok(r, "sge_lanczos"));
            let mb = sel.iter().map(|r| r.las["bfs_lanczos"].total).sum::<f64>()
                / sel.len().max(1) as f64;
            let ms = sel.iter().map(|r| r.las["sge_lanczos"].total).sum::<f64>()
                / sel.len().max(1) as f64;
            if complete {
                if mb <= best_b.0 {
                    best_b = (mb, kap);
                }
                if ms <= best_s.0 {
                    best_s = (ms, kap);
                }
            } // ties: later (larger) kappa wins
            table.push(json!({"kappa": kap, "n": sel.len(), "complete": complete, "mean_total_bl_pipeline": mb, "mean_total_a2_pipeline": ms,
                "mean_rel": sel.iter().map(|r| r.rel).sum::<f64>() / sel.len().max(1) as f64}));
        }
        let ks = |b: (f64, f64)| {
            if b.0.is_finite() {
                json!(b.1)
            } else {
                Value::Null
            }
        }; // no complete kappa: no choice
        let _ = writeln!(md, "k={} kappa*={} kappa**={}", k, ks(best_b), ks(best_s));
        choice.insert(
            k.to_string(),
            json!({"kappa_star": ks(best_b), "kappa_2star": ks(best_s), "table": table}),
        );
    }
    let doc = json!({"protocol": "research/ic_improve_20261014/04b-genus3-2lp-rebalance/PROTOCOL.md §4", "grid": GRID, "rule": "arg min of mean total_mults over training curves 5-7 x targets 1-2; ties to the larger kappa", "choice": choice});
    std::fs::write(out, serde_json::to_string_pretty(&doc).unwrap() + "\n").unwrap();
    print!("{}", md);
}

// ------------------------------------------------------------------ statistics
fn ols(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = x.iter().map(|a| (a - mx) * (a - mx)).sum();
    sxy / sxx
}
fn pct(v: &mut [f64], lo: f64, hi: f64) -> (f64, f64) {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    let i = ((lo * n as f64).floor() as usize).min(n - 1);
    let j = (((hi * n as f64).ceil() as usize).max(1) - 1).min(n - 1);
    (v[i], v[j])
}

const ARMS: [&str; 4] = ["BL", "A1", "A2", "A3"];
const PHASES: [&str; 3] = ["total", "rel", "la"];

/// cost[arm][phase][k] = per-instance vector (instances in fixed order); None where the arm is not defined
type Costs = BTreeMap<(String, String, u32), Vec<f64>>;

fn arm_record<'a>(
    runs: &'a [Run],
    choice: &Value,
    arm: &str,
    k: u32,
    c: u32,
    t: u32,
) -> Option<(&'a Run, &'static str)> {
    let (kap, m) = match arm {
        "BL" => (2.0, "bfs_lanczos"),
        "A2" => (2.0, "sge_lanczos"),
        "A1" => {
            if k < 14 {
                return None;
            }
            (
                choice["choice"][k.to_string()]["kappa_star"].as_f64()?,
                "bfs_lanczos",
            )
        }
        "A3" => {
            if k < 14 {
                return None;
            }
            (
                choice["choice"][k.to_string()]["kappa_2star"].as_f64()?,
                "sge_lanczos",
            )
        }
        _ => unreachable!(),
    };
    runs.iter()
        .find(|r| r.k == k && r.curve == c && r.target == t && kap_eq(r.kappa, kap))
        .map(|r| (r, m))
}

fn stats_for(
    costs: &Costs,
    qs: &BTreeMap<u32, f64>,
    idx: &BTreeMap<u32, Vec<usize>>,
    arm: &str,
    phase: &str,
) -> (
    BTreeMap<u32, f64>,
    BTreeMap<u32, f64>,
    Option<f64>,
    Option<f64>,
) {
    let mut mean = BTreeMap::new();
    for (&k, ix) in idx {
        if let Some(v) = costs.get(&(arm.to_string(), phase.to_string(), k)) {
            mean.insert(k, ix.iter().map(|&i| v[i]).sum::<f64>() / ix.len() as f64);
        }
    }
    let mut loc = BTreeMap::new();
    for (&k, &c) in &mean {
        if let Some(&c2) = mean.get(&(k + 1)) {
            loc.insert(k, (c2 / c).ln() / (qs[&(k + 1)] / qs[&k]).ln());
        }
    }
    let fit_ks: Vec<u32> = (14..=20).filter(|k| mean.contains_key(k)).collect();
    let fit = if fit_ks.len() == 7 {
        Some(ols(
            &fit_ks.iter().map(|k| qs[k].ln()).collect::<Vec<_>>(),
            &fit_ks.iter().map(|k| mean[k].ln()).collect::<Vec<_>>(),
        ))
    } else {
        None
    };
    let tr_ks: Vec<u32> = (14..=19).filter(|k| loc.contains_key(k)).collect();
    let trend = if tr_ks.len() == 6 {
        Some(ols(
            &tr_ks
                .iter()
                .map(|k| (qs[k].log2() + qs[&(k + 1)].log2()) / 2.0)
                .collect::<Vec<_>>(),
            &tr_ks.iter().map(|k| loc[k]).collect::<Vec<_>>(),
        ))
    } else {
        None
    };
    (mean, loc, fit, trend)
}

fn cmd_report(dir: &str, rho_dir: &str, choice_path: &str, out: &str, md_path: &str) {
    let rs = runs(dir);
    let choice: Value =
        serde_json::from_str(&std::fs::read_to_string(choice_path).unwrap()).unwrap();
    let inst: Vec<(u32, u32)> = (1..=4).flat_map(|c| (1..=3).map(move |t| (c, t))).collect();
    let mut qs = BTreeMap::new();
    let mut ells: BTreeMap<(u32, u32), u64> = BTreeMap::new();
    for r in &rs {
        qs.insert(r.k, r.q as f64);
        ells.insert((r.k, r.curve), r.ell);
    }
    // ---- completeness and verification
    let mut costs: Costs = BTreeMap::new();
    let mut weights: BTreeMap<(String, u32), Vec<(f64, f64, f64, f64, f64)>> = BTreeMap::new();
    let mut missing = vec![];
    let mut failed = vec![];
    let mut acct_bad = vec![];
    let mut acct_checked = 0u64;
    let mut kappas_used: BTreeMap<(String, u32), f64> = BTreeMap::new();
    let mut rss_max = 0f64;
    let mut wall_max = 0f64;
    for r in &rs {
        rss_max = rss_max.max(r.rss);
        wall_max = wall_max.max(r.wall);
        for (m, l) in &r.las {
            if !(l.verified && l.planted) {
                failed.push(format!(
                    "k={} curve={} target={} kappa={} {}",
                    r.k, r.curve, r.target, r.kappa, m
                ));
            }
            // exact Lanczos accounting (PROTOCOL §8): attempts = 1
            if l.attempts == 1 {
                acct_checked += 1;
                let e = r.ell - 2;
                let inv = (64 - e.leading_zeros() as u64 - 1) + (e.count_ones() as u64 - 1);
                let (n, nnz, it) = (l.rows, l.nnz, l.iters as f64);
                let expect = it * (2.0 * nnz + 7.0 * n + 3.0 + inv as f64) - (2.0 * n + 1.0)
                    + 3.0 * nnz
                    + 2.0 * n
                    + inv as f64
                    + 1.0;
                if (l.lanczos - expect).abs() > 0.5 {
                    acct_bad.push(format!(
                        "k={} c={} t={} kappa={} {}: counted {} expected {}",
                        r.k, r.curve, r.target, r.kappa, m, l.lanczos, expect
                    ));
                }
            }
        }
    }
    let ks: Vec<u32> = qs.keys().cloned().collect();
    for arm in ARMS {
        for &k in &ks {
            let mut v: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
            let mut w = vec![];
            let mut complete = true;
            for &(c, t) in &inst {
                match arm_record(&rs, &choice, arm, k, c, t) {
                    Some((r, m)) if ok(r, m) => {
                        let l = &r.las[m];
                        v.entry("total").or_default().push(l.total);
                        v.entry("rel").or_default().push(r.rel);
                        v.entry("la").or_default().push(l.la);
                        w.push((l.w_row, l.w_col, l.rows, l.cols, l.nnz));
                        kappas_used.insert((arm.to_string(), k), r.kappa);
                    }
                    Some(_) => {
                        complete = false;
                    }
                    None => {
                        if arm == "BL" || arm == "A2" || k >= 14 {
                            missing.push(format!("{} k={} curve={} target={}", arm, k, c, t));
                            complete = false;
                        } else {
                            complete = false;
                        }
                    }
                }
            }
            if complete {
                for p in PHASES {
                    costs.insert((arm.to_string(), p.to_string(), k), v[p].clone());
                }
                weights.insert((arm.to_string(), k), w);
            }
        }
    }
    // ---- point estimates
    let full_idx: BTreeMap<u32, Vec<usize>> =
        ks.iter().map(|&k| (k, (0..inst.len()).collect())).collect();
    let mut point = serde_json::Map::new();
    for arm in ARMS {
        for p in PHASES {
            let (mean, loc, fit, trend) = stats_for(&costs, &qs, &full_idx, arm, p);
            point.insert(format!("{}/{}", arm, p), json!({"mean": mean.iter().map(|(k, v)| (k.to_string(), json!(v))).collect::<serde_json::Map<_, _>>(),
            "local": loc.iter().map(|(k, v)| (k.to_string(), json!(v))).collect::<serde_json::Map<_, _>>(), "fit_14_20": fit, "trend_14_19": trend}));
        }
    }
    // paired ratios (geometric mean) per size, per phase
    let ratio = |idx: &BTreeMap<u32, Vec<usize>>, arm: &str, p: &str, k: u32| -> Option<f64> {
        let a = costs.get(&(arm.to_string(), p.to_string(), k))?;
        let b = costs.get(&("BL".to_string(), p.to_string(), k))?;
        let ix = &idx[&k];
        Some((ix.iter().map(|&i| (a[i] / b[i]).ln()).sum::<f64>() / ix.len() as f64).exp())
    };
    // ---- bootstrap (paired: one resample of instances per size shared by every arm)
    let mut rng = Rng(BOOT_SEED);
    let mut bs: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for _ in 0..BOOT_B {
        let idx: BTreeMap<u32, Vec<usize>> = ks
            .iter()
            .map(|&k| (k, (0..inst.len()).map(|_| rng.below(inst.len())).collect()))
            .collect();
        let mut fits = HashMap::new();
        for arm in ARMS {
            for p in PHASES {
                let (mean, loc, fit, trend) = stats_for(&costs, &qs, &idx, arm, p);
                for (k, v) in &mean {
                    bs.entry(format!("{}/{}/mean/{}", arm, p, k))
                        .or_default()
                        .push(*v);
                }
                for (k, v) in &loc {
                    bs.entry(format!("{}/{}/local/{}", arm, p, k))
                        .or_default()
                        .push(*v);
                }
                if let Some(x) = fit {
                    bs.entry(format!("{}/{}/fit", arm, p)).or_default().push(x);
                    fits.insert((arm, p), x);
                }
                if let Some(x) = trend {
                    bs.entry(format!("{}/{}/trend", arm, p))
                        .or_default()
                        .push(x);
                }
                if arm != "BL" {
                    for &k in &ks {
                        if let Some(r) = ratio(&idx, arm, p, k) {
                            bs.entry(format!("{}/{}/ratio/{}", arm, p, k))
                                .or_default()
                                .push(r);
                        }
                    }
                }
            }
        }
        for arm in ["A1", "A2", "A3"] {
            for p in PHASES {
                if let (Some(a), Some(b)) = (fits.get(&(arm, p)), fits.get(&("BL", p))) {
                    bs.entry(format!("{}/{}/fitdiff", arm, p))
                        .or_default()
                        .push(a - b);
                }
            }
        }
    }
    let mut ci = serde_json::Map::new();
    for (key, v) in bs.iter_mut() {
        let (l95, h95) = pct(v, 0.025, 0.975);
        let (l99, h99) = pct(v, 0.005, 0.995);
        ci.insert(
            key.clone(),
            json!({"lo95": l95, "hi95": h95, "lo99": l99, "hi99": h99}),
        );
    }
    let pr = |key: &str| -> serde_json::Map<String, Value> {
        ci.get(key)
            .and_then(|v| v.as_object().cloned())
            .unwrap_or_default()
    };
    // point ratios
    let mut ratios = serde_json::Map::new();
    for arm in ["A1", "A2", "A3"] {
        for p in PHASES {
            for &k in &ks {
                if let Some(r) = ratio(&full_idx, arm, p, k) {
                    ratios.insert(
                        format!("{}/{}/{}", arm, p, k),
                        json!({"geo_mean": r, "ci": pr(&format!("{}/{}/ratio/{}", arm, p, k))}),
                    );
                }
            }
        }
    }
    // ---- A/A control: BL fit on target 1 vs targets 2-3
    let sub = |ts: &[u32]| -> Vec<usize> {
        inst.iter()
            .enumerate()
            .filter(|(_, x)| ts.contains(&x.1))
            .map(|(i, _)| i)
            .collect()
    };
    let (s1, s23) = (sub(&[1]), sub(&[2, 3]));
    let mut rng2 = Rng(BOOT_SEED ^ 0xAA);
    let mut aa = vec![];
    for _ in 0..BOOT_B {
        let i1: BTreeMap<u32, Vec<usize>> = ks
            .iter()
            .map(|&k| (k, (0..s1.len()).map(|_| s1[rng2.below(s1.len())]).collect()))
            .collect();
        let i2: BTreeMap<u32, Vec<usize>> = ks
            .iter()
            .map(|&k| {
                (
                    k,
                    (0..s23.len()).map(|_| s23[rng2.below(s23.len())]).collect(),
                )
            })
            .collect();
        if let (Some(a), Some(b)) = (
            stats_for(&costs, &qs, &i1, "BL", "total").2,
            stats_for(&costs, &qs, &i2, "BL", "total").2,
        ) {
            aa.push(a - b);
        }
    }
    let aa_point = match (
        stats_for(
            &costs,
            &qs,
            &ks.iter().map(|&k| (k, s1.clone())).collect(),
            "BL",
            "total",
        )
        .2,
        stats_for(
            &costs,
            &qs,
            &ks.iter().map(|&k| (k, s23.clone())).collect(),
            "BL",
            "total",
        )
        .2,
    ) {
        (Some(a), Some(b)) => Some(a - b),
        _ => None,
    };
    let aa_ci = if aa.is_empty() {
        (f64::NAN, f64::NAN)
    } else {
        pct(&mut aa, 0.025, 0.975)
    };
    // ---- negative control: A2 relation phase == BL relation phase
    let mut pairing_ok = true;
    for &k in &ks {
        if let (Some(a), Some(b)) = (
            costs.get(&("A2".into(), "rel".into(), k)),
            costs.get(&("BL".into(), "rel".into(), k)),
        ) {
            if a != b {
                pairing_ok = false;
            }
        }
    }
    // ---- decisions (PROTOCOL §7)
    let g =
        |key: &str, f: &str| -> f64 { ci.get(key).and_then(|v| v[f].as_f64()).unwrap_or(f64::NAN) };
    let four_thirds = 4.0 / 3.0;
    let bl_fit = point["BL/total"]["fit_14_20"].as_f64();
    let half95 = (g("BL/total/fit", "hi95") - g("BL/total/fit", "lo95")) / 2.0;
    let d0 = g("BL/total/fit", "lo99") > four_thirds || g("BL/total/trend", "lo99") > 0.0;
    let mut verdicts = serde_json::Map::new();
    for arm in ["A1", "A2", "A3"] {
        let e_lo = g(&format!("{}/total/fit", arm), "lo99");
        let e_hi = g(&format!("{}/total/fit", arm), "hi99");
        let t_lo = g(&format!("{}/total/trend", arm), "lo99");
        let _t_hi = g(&format!("{}/total/trend", arm), "hi99");
        let d_lo = g(&format!("{}/total/fitdiff", arm), "lo99");
        let d_hi = g(&format!("{}/total/fitdiff", arm), "hi99");
        let r_lo = g(&format!("{}/total/ratio/20", arm), "lo99");
        let r_hi = g(&format!("{}/total/ratio/20", arm), "hi99");
        let removes = e_lo <= four_thirds && t_lo <= 0.0;
        let v = if e_lo.is_nan() {
            "incomplete".to_string()
        } else if removes {
            "removes the departure".into()
        } else if d_hi < 0.0 {
            "reduces the exponent".into()
        } else if d_lo > 0.0 {
            "raises the exponent".into()
        } else if r_hi < 1.0 || r_lo > 1.0 {
            "constant only".into()
        } else {
            "no measured effect".into()
        };
        let _ = e_hi;
        verdicts.insert(arm.to_string(), json!({"verdict": v, "E99": [e_lo, e_hi], "T99_lo": t_lo, "E_minus_BL_99": [d_lo, d_hi], "ratio20_99": [r_lo, r_hi], "secondary": arm == "A3"}));
    }
    // ---- rho reference
    let rho: Vec<Value> = read_records(rho_dir)
        .into_iter()
        .filter(|v| v["cmd"] == "rho")
        .collect();
    let mut rho_tab = serde_json::Map::new();
    for &k in &ks {
        let sel: Vec<&Value> = rho.iter().filter(|v| u(v, "k") as u32 == k).collect();
        if sel.is_empty() {
            continue;
        }
        let okn = sel
            .iter()
            .filter(|v| {
                v["log_verified"].as_bool().unwrap() && v["log_equals_planted"].as_bool().unwrap()
            })
            .count();
        rho_tab.insert(k.to_string(), json!({"n": sel.len(), "verified": okn, "mean_S": sel.iter().map(|v| f(v, "S")).sum::<f64>() / sel.len() as f64,
            "mean_steps_per_sqrt_ell": sel.iter().map(|v| f(v, "steps_per_sqrt_ell")).sum::<f64>() / sel.len() as f64}));
    }
    // S per arm (mean over instances of total / sqrt(ell))
    let mut s_tab = serde_json::Map::new();
    for arm in ARMS {
        for &k in &ks {
            if let Some(v) = costs.get(&(arm.to_string(), "total".to_string(), k)) {
                let s: f64 = inst
                    .iter()
                    .enumerate()
                    .map(|(i, &(c, _))| v[i] / (ells[&(k, c)] as f64).sqrt())
                    .sum::<f64>()
                    / inst.len() as f64;
                s_tab.insert(format!("{}/{}", arm, k), json!(s));
            }
        }
    }
    // IC / rho in S (descriptive, PROTOCOL §8): raw, and with rho's measured walk renormalised to the ideal sqrt(pi/2)
    let mut icrho = serde_json::Map::new();
    for (k, r) in &rho_tab {
        for arm in ["BL", "A2"] {
            if let Some(sv) = s_tab.get(&format!("{}/{}", arm, k)) {
                let (srho, steps) = (
                    r["mean_S"].as_f64().unwrap(),
                    r["mean_steps_per_sqrt_ell"].as_f64().unwrap(),
                );
                let ideal = (std::f64::consts::PI / 2.0).sqrt();
                icrho.insert(
                    format!("{}/{}", arm, k),
                    json!({"raw": sv.as_f64().unwrap() / srho, "rho_renormalised": sv.as_f64().unwrap() / (srho * ideal / steps), "rho_steps_per_sqrt_ell": steps, "rho_n": r["n"]}),
                );
            }
        }
    }
    let mut wt = serde_json::Map::new();
    for ((arm, k), w) in &weights {
        let n = w.len() as f64;
        wt.insert(format!("{}/{}", arm, k), json!({"w_row": w.iter().map(|x| x.0).sum::<f64>() / n, "w_col": w.iter().map(|x| x.1).sum::<f64>() / n,
        "rows": w.iter().map(|x| x.2).sum::<f64>() / n, "cols": w.iter().map(|x| x.3).sum::<f64>() / n, "nnz": w.iter().map(|x| x.4).sum::<f64>() / n}));
    }
    let doc = json!({
        "protocol": "research/ic_improve_20261014/04b-genus3-2lp-rebalance/PROTOCOL.md", "bootstrap": {"B": BOOT_B, "seed": BOOT_SEED, "unit": "12 test instances resampled within each size, shared by all arms"},
        "q": qs.iter().map(|(k, v)| (k.to_string(), json!(v))).collect::<serde_json::Map<_, _>>(),
        "kappa_used": kappas_used.iter().map(|((a, k), v)| (format!("{}/{}", a, k), json!(v))).collect::<serde_json::Map<_, _>>(),
        "point": point, "ci": ci, "ratios": ratios, "S": s_tab, "weights": wt, "rho": rho_tab, "ic_over_rho": icrho,
        "controls": {"pairing_rel_identical": pairing_ok, "aa_fit_diff_point": aa_point, "aa_fit_diff_95": [aa_ci.0, aa_ci.1], "lanczos_accounting_checked": acct_checked, "lanczos_accounting_mismatches": acct_bad,
                     "power_half_width_95_E_BL": half95, "power_ok": half95 <= 0.03},
        "decision": {"D0_baseline_departs": d0, "E_BL": bl_fit, "E_BL_99": [g("BL/total/fit", "lo99"), g("BL/total/fit", "hi99")], "T_BL_99": [g("BL/total/trend", "lo99"), g("BL/total/trend", "hi99")], "arms": verdicts},
        "integrity": {"runs": rs.len(), "failed_logs": failed, "missing": missing, "max_rss_mb": rss_max, "max_wall_s": wall_max}
    });
    std::fs::write(out, serde_json::to_string_pretty(&doc).unwrap() + "\n").unwrap();
    // ---- markdown tables
    let mut md = String::new();
    let fmt_ci = |key: &str| -> String {
        let lo = g(key, "lo95");
        let hi = g(key, "hi95");
        if lo.is_nan() {
            "—".into()
        } else {
            format!("[{:.3}, {:.3}]", lo, hi)
        }
    };
    let pt = |arm: &str, p: &str, what: &str, k: u32| -> Option<f64> {
        point[&format!("{}/{}", arm, p)][what][&k.to_string()].as_f64()
    };
    let _ = writeln!(md, "### Mean total cost per size (multiplications; 12 test instances; 95% bootstrap interval of the mean)\n");
    let _ = writeln!(md, "| k | q | BL κ=2 | A1 (κ*) | A2 SGE κ=2 | A3 (κ**) | S_BL | S_A1 | S_A2 | S_A3 | S_rho (measured) |\n|---|---|---|---|---|---|---|---|---|---|---|");
    for &k in &ks {
        let mut row = format!("| {} | {} |", k, qs[&k]);
        for arm in ARMS {
            row += &match pt(arm, "total", "mean", k) {
                Some(m) => format!(
                    " {:.4e} (κ={}) |",
                    m,
                    kappas_used
                        .get(&(arm.to_string(), k))
                        .map(|x| x.to_string())
                        .unwrap_or("?".into())
                ),
                None => " — |".into(),
            };
        }
        for arm in ARMS {
            row += &match s_tab.get(&format!("{}/{}", arm, k)) {
                Some(s) => format!(" {:.1} |", s.as_f64().unwrap()),
                None => " — |".into(),
            };
        }
        row += &match rho_tab.get(&k.to_string()) {
            Some(r) => format!(
                " {:.1} (n={}, {} verified) |",
                r["mean_S"].as_f64().unwrap(),
                r["n"],
                r["verified"]
            ),
            None => " — |".into(),
        };
        let _ = writeln!(md, "{}", row);
    }
    for p in PHASES {
        let _ = writeln!(
            md,
            "\n### Local exponents, phase `{}` (point [95% interval])\n",
            p
        );
        let _ = writeln!(md, "| k→k+1 | BL | A1 | A2 | A3 |\n|---|---|---|---|---|");
        for &k in &ks {
            if k == 20 {
                continue;
            }
            let mut row = format!("| {}→{} |", k, k + 1);
            for arm in ARMS {
                row += &match pt(arm, p, "local", k) {
                    Some(e) => format!(
                        " {:.3} {} |",
                        e,
                        fmt_ci(&format!("{}/{}/local/{}", arm, p, k))
                    ),
                    None => " — |".into(),
                };
            }
            let _ = writeln!(md, "{}", row);
        }
    }
    let _ = writeln!(
        md,
        "\n### Fit over 2^14–2^20 and trend of local exponents (point [95%] [99%])\n"
    );
    let _ = writeln!(md, "| arm | phase | E (fit) | E − 4/3 | E − E_BL | trend T (per log₂ q) |\n|---|---|---|---|---|---|");
    for arm in ARMS {
        for p in PHASES {
            let key = format!("{}/{}", arm, p);
            let e = point[&key]["fit_14_20"].as_f64();
            let t = point[&key]["trend_14_19"].as_f64();
            let e_bl = point[&format!("BL/{}", p)]["fit_14_20"].as_f64();
            let ff = |x: Option<f64>, kk: &str| -> String {
                match x {
                    Some(v) => format!(
                        "{:.3} [{:.3}, {:.3}] [{:.3}, {:.3}]",
                        v,
                        g(kk, "lo95"),
                        g(kk, "hi95"),
                        g(kk, "lo99"),
                        g(kk, "hi99")
                    ),
                    None => "—".into(),
                }
            };
            let diff = match (e, e_bl) {
                (Some(a), Some(b)) if arm != "BL" => format!(
                    "{:+.3} [{:+.3}, {:+.3}] [{:+.3}, {:+.3}]",
                    a - b,
                    g(&format!("{}/fitdiff", key), "lo95"),
                    g(&format!("{}/fitdiff", key), "hi95"),
                    g(&format!("{}/fitdiff", key), "lo99"),
                    g(&format!("{}/fitdiff", key), "hi99")
                ),
                _ => "—".into(),
            };
            let _ = writeln!(
                md,
                "| {} | {} | {} | {} | {} | {} |",
                arm,
                p,
                ff(e, &format!("{}/fit", key)),
                e.map(|v| format!("{:+.3}", v - four_thirds))
                    .unwrap_or("—".into()),
                diff,
                ff(t, &format!("{}/trend", key))
            );
        }
    }
    let _ = writeln!(
        md,
        "\n### Paired cost ratios against BL (geometric mean over the 12 instances [95%] [99%])\n"
    );
    let _ = writeln!(
        md,
        "| k | A1/BL total | A2/BL total | A2/BL la | A3/BL total |\n|---|---|---|---|---|"
    );
    for &k in &ks {
        let cell = |arm: &str, p: &str| -> String {
            match ratios.get(&format!("{}/{}/{}", arm, p, k)) {
                Some(r) => {
                    let c = &r["ci"];
                    format!(
                        "{:.3} [{:.3}, {:.3}] [{:.3}, {:.3}]",
                        r["geo_mean"].as_f64().unwrap(),
                        c["lo95"].as_f64().unwrap_or(f64::NAN),
                        c["hi95"].as_f64().unwrap_or(f64::NAN),
                        c["lo99"].as_f64().unwrap_or(f64::NAN),
                        c["hi99"].as_f64().unwrap_or(f64::NAN)
                    )
                }
                None => "—".into(),
            }
        };
        let _ = writeln!(
            md,
            "| {} | {} | {} | {} | {} |",
            k,
            cell("A1", "total"),
            cell("A2", "total"),
            cell("A2", "la"),
            cell("A3", "total")
        );
    }
    let _ = writeln!(md, "\n### Matrix entering Lanczos (mean over instances)\n");
    let _ = writeln!(md, "| k | BL rows | BL nnz/row | BL nnz/col | A2 rows | A2 nnz/row | A2 nnz/col | A1 κ* rows | A1 nnz/row |\n|---|---|---|---|---|---|---|---|---|");
    for &k in &ks {
        let w = |arm: &str, f: &str| -> String {
            wt.get(&format!("{}/{}", arm, k))
                .map(|x| format!("{:.2}", x[f].as_f64().unwrap()))
                .unwrap_or("—".into())
        };
        let _ = writeln!(
            md,
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} |",
            k,
            w("BL", "rows"),
            w("BL", "w_row"),
            w("BL", "w_col"),
            w("A2", "rows"),
            w("A2", "w_row"),
            w("A2", "w_col"),
            w("A1", "rows"),
            w("A1", "w_row")
        );
    }
    let _ = writeln!(md, "\n### Decision (PROTOCOL §7)\n");
    let _ = writeln!(
        md,
        "- D0, baseline departs: **{}** (E_BL 99% [{:.3}, {:.3}], T_BL 99% [{:.4}, {:.4}])",
        d0,
        g("BL/total/fit", "lo99"),
        g("BL/total/fit", "hi99"),
        g("BL/total/trend", "lo99"),
        g("BL/total/trend", "hi99")
    );
    for arm in ["A1", "A2", "A3"] {
        let v = &verdicts[arm];
        let _ = writeln!(
            md,
            "- {}{}: **{}**",
            arm,
            if arm == "A3" { " (secondary)" } else { "" },
            v["verdict"].as_str().unwrap()
        );
    }
    let _ = writeln!(
        md,
        "- Power: half-width of the 95% interval of E_BL = {:.4} ({})",
        half95,
        if half95 <= 0.03 {
            "≤ 0.03, resolves the smallest departure of interest"
        } else {
            "> 0.03: straddling verdicts are inconclusive at this resolution"
        }
    );
    let _ = writeln!(md, "- Controls: A2 relation phase identical to BL on every instance: {}; A/A difference of E_BL (target 1 vs targets 2–3): {} [95% {:.3}, {:.3}]; Lanczos accounting checked on {} solves, mismatches {}",
        pairing_ok, aa_point.map(|x| format!("{:+.3}", x)).unwrap_or("—".into()), aa_ci.0, aa_ci.1, acct_checked, acct_bad.len());
    let _ = writeln!(md, "- Integrity: {} solve records; failed logs {}; missing arm instances {}; max RSS {:.0} MB; max job wall {:.0} s", rs.len(), failed.len(), missing.len(), rss_max, wall_max);
    std::fs::write(md_path, md).unwrap();
    println!("wrote {} and {}", out, md_path);
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let get = |k: &str| -> String {
        a.iter()
            .position(|x| x == k)
            .map(|i| a[i + 1].clone())
            .unwrap_or_else(|| panic!("missing {}", k))
    };
    match a.get(1).map(|s| s.as_str()) {
        Some("tune") => cmd_tune(&get("--dir"), &get("--out")),
        Some("report") => cmd_report(
            &get("--dir"),
            &get("--rho"),
            &get("--choice"),
            &get("--out"),
            &get("--md"),
        ),
        _ => {
            eprintln!(
                "usage: tune --dir D --out F | report --dir D --rho D --choice F --out F --md F"
            );
            std::process::exit(2);
        }
    }
}
