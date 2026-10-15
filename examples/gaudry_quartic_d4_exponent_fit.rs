//! Fit the measured `k = 4` exponents of `research/gaudry_quartic_d4_exponent_20261015`
//! (PROTOCOL.md, "Analysis") from the committed ecbench records, and give
//! the crossover with rho from those fits.
//!
//! For every IC arm and rung: the mean, over the rung's targets, of the
//! relations phase, the linear algebra (Wiedemann, priced), the other phases
//! (factor base, set-up, verification) and the total, in GAE; and, derived,
//! the dense elimination that ran on the same rows as a shadow, priced at
//! the same `row_op` price.  Each is fitted as `ln(mean) = c + β ln p` by
//! OLS over the rungs, with a Student-`t` interval and a bootstrap interval
//! (10,000 replicates resampling targets within each rung, seed 1); the
//! headline interval is their union.  Rho is fitted the same way, and
//! `S_rho = gae/√r` is pooled.  The crossover solves
//! `Σ_k A_k p^{β_k} = S̄_rho √(p⁴/h)` per replicate.  Writes
//! `raw/summary.json`.
//!
//! ```sh
//! cargo run --release --example gaudry_quartic_d4_exponent_fit -- research/gaudry_quartic_d4_exponent_20261015
//! ```

use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;
use std::path::Path;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

const REPLICATES: usize = 10_000;
const IC_ARMS: [&str; 3] = ["ic-d4-w", "ic-d4-lp2", "ic-d4-lp1"];
/// Quantities fitted per IC arm; the first three sum to `total`.
const QUANTITIES: [&str; 6] = [
    "relations",
    "linear_algebra",
    "other",
    "total",
    "shadow_gauss",
    "trials",
];
const HS: [f64; 3] = [2.0, 4.0, 16.0];

#[derive(Clone, Debug)]
struct IcRec {
    p: u64,
    q: BTreeMap<&'static str, f64>,
    s: f64,
    workload: String,
    native: BTreeMap<String, f64>,
}

#[derive(Clone, Debug)]
struct RhoRec {
    p: u64,
    gae: f64,
    s: f64,
}

fn t975(df: usize) -> f64 {
    const T: [f64; 30] = [
        12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160,
        2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056,
        2.052, 2.048, 2.045, 2.042,
    ];
    T[(df.max(1) - 1).min(29)]
}

/// OLS of `y` on `x`: `(c, β, se(β), se(c))`.
fn ols(x: &[f64], y: &[f64]) -> (f64, f64, f64, f64) {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let sxx: f64 = x.iter().map(|a| (a - mx).powi(2)).sum();
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let b = sxy / sxx;
    let c = my - b * mx;
    let rss: f64 = x.iter().zip(y).map(|(a, v)| (v - c - b * a).powi(2)).sum();
    let s2 = if x.len() > 2 { rss / (n - 2.0) } else { 0.0 };
    let se_b = (s2 / sxx).sqrt();
    let se_c = (s2 * (1.0 / n + mx * mx / sxx)).sqrt();
    (c, b, se_b, se_c)
}

fn pct(v: &mut [f64], q: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let i = ((v.len() - 1) as f64 * q).round() as usize;
    v[i]
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

fn num(v: &Value) -> f64 {
    v.as_f64()
        .or_else(|| v.as_str().and_then(|s| s.parse().ok()))
        .unwrap_or(0.0)
}

fn load(dir: &Path) -> (Vec<(String, IcRec)>, Vec<RhoRec>, Value) {
    let mut ic = Vec::new();
    let mut rho = Vec::new();
    let mut seen_rho = BTreeSet::new();
    let mut status: BTreeMap<String, BTreeMap<String, u64>> = BTreeMap::new();
    let mut sessions = Vec::new();
    let mut entries: Vec<_> = fs::read_dir(dir.join("sessions"))
        .expect("sessions/")
        .flatten()
        .filter(|e| e.path().join("records.jsonl").exists())
        .collect();
    entries.sort_by_key(|e| e.file_name());
    for e in entries {
        let text = fs::read_to_string(e.path().join("records.jsonl")).unwrap();
        let mut sid = String::new();
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            let r: Value = serde_json::from_str(line).expect("record");
            sid = r["session_id"].as_str().unwrap_or("").to_string();
            let arm = r["arm"].as_str().unwrap().to_string();
            let st = r["outcome"]["status"].as_str().unwrap_or("?").to_string();
            *status
                .entry(arm.clone())
                .or_default()
                .entry(st.clone())
                .or_default() += 1;
            if st != "verified" {
                continue;
            }
            let p = r["workload"]["curve_spec"]["p"].as_u64().unwrap();
            let r_ord = num(&r["workload"]["curve"]["r"]);
            let gae = num(&r["cost"]["total_gae"]);
            let s = num(&r["cost"]["s"]);
            let workload = r["workload"]["workload_id"].as_str().unwrap().to_string();
            if arm == "rho-neg" {
                let key = (workload.clone(), r["algorithm_seed"].to_string());
                if seen_rho.insert(key) {
                    rho.push(RhoRec { p, gae, s });
                }
                let _ = r_ord;
                continue;
            }
            let phase = |name: &str| -> &Value {
                r["phases"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .find(|ph| ph["name"] == name)
                    .unwrap_or(&Value::Null)
            };
            let rel = num(&phase("relations")["gae"]);
            let la = num(&phase("linear_algebra")["gae"]);
            let la_ops = num(&phase("linear_algebra")["native"]["row_ops"]);
            let shadow_ops = num(&phase("linear_algebra")["native"]["shadow_gauss_row_ops"]);
            let shadow = if la_ops > 0.0 {
                shadow_ops * la / la_ops
            } else {
                0.0
            };
            let mut native = BTreeMap::new();
            for ph in ["relations", "linear_algebra"] {
                if let Some(m) = phase(ph)["native"].as_object() {
                    for (k, v) in m {
                        native.insert(format!("{ph}.{k}"), num(v));
                    }
                }
            }
            let trials = native.get("relations.trials").copied().unwrap_or(0.0);
            let mut q = BTreeMap::new();
            q.insert("relations", rel);
            q.insert("linear_algebra", la);
            q.insert("other", gae - rel - la);
            q.insert("total", gae);
            q.insert("shadow_gauss", shadow);
            q.insert("trials", trials);
            ic.push((
                arm,
                IcRec {
                    p,
                    q,
                    s,
                    workload,
                    native,
                },
            ));
        }
        sessions.push(json!({"dir": e.file_name().to_string_lossy(), "session_id": sid}));
    }
    (ic, rho, json!({"status": status, "sessions": sessions}))
}

/// Group record indices by rung.
fn by_rung<T>(recs: &[T], p_of: impl Fn(&T) -> u64) -> BTreeMap<u64, Vec<usize>> {
    let mut m: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    for (i, r) in recs.iter().enumerate() {
        m.entry(p_of(r)).or_default().push(i);
    }
    m
}

/// Per-rung means of a value under a choice of record indices per rung.
fn rung_means(
    rungs: &BTreeMap<u64, Vec<usize>>,
    pick: &BTreeMap<u64, Vec<usize>>,
    val: impl Fn(usize) -> f64,
) -> Vec<(f64, f64)> {
    rungs
        .keys()
        .map(|p| {
            let v: Vec<f64> = pick[p].iter().map(|&i| val(i)).collect();
            ((*p as f64).ln(), mean(&v))
        })
        .collect()
}

fn fit_points(pts: &[(f64, f64)]) -> Option<(f64, f64, f64, f64)> {
    let pts: Vec<_> = pts.iter().filter(|(_, y)| *y > 0.0).collect();
    if pts.len() < 3 {
        return None;
    }
    let x: Vec<f64> = pts.iter().map(|(a, _)| *a).collect();
    let y: Vec<f64> = pts.iter().map(|(_, b)| b.ln()).collect();
    Some(ols(&x, &y))
}

fn resample(rungs: &BTreeMap<u64, Vec<usize>>, rng: &mut StdRng) -> BTreeMap<u64, Vec<usize>> {
    rungs
        .iter()
        .map(|(p, idx)| {
            (
                *p,
                (0..idx.len())
                    .map(|_| idx[rng.gen_range(0..idx.len())])
                    .collect(),
            )
        })
        .collect()
}

/// ln of the IC cost minus ln of rho's at `x = ln p`.
fn gap(ic: &[(f64, f64)], rho: &dyn Fn(f64) -> f64, x: f64) -> f64 {
    let s: f64 = ic.iter().map(|(c, b)| (c + b * x).exp()).sum();
    s.ln() - rho(x)
}

/// Crossings of the IC cost with rho for `ln p ∈ [ln 29, 64 ln 2]`:
/// the first `p` where IC falls below rho and, if it rises above again,
/// where.
fn crossings(ic: &[(f64, f64)], rho: &dyn Fn(f64) -> f64) -> (Option<f64>, Option<f64>) {
    let (lo, hi) = ((29f64).ln(), 64.0 * 2f64.ln());
    let steps = 4000;
    let mut out = (None, None);
    let mut prev = gap(ic, rho, lo);
    for i in 1..=steps {
        let x = lo + (hi - lo) * i as f64 / steps as f64;
        let g = gap(ic, rho, x);
        if prev.signum() != g.signum() {
            let (mut a, mut b) = (x - (hi - lo) / steps as f64, x);
            for _ in 0..60 {
                let m = 0.5 * (a + b);
                if gap(ic, rho, m).signum() == gap(ic, rho, a).signum() {
                    a = m;
                } else {
                    b = m;
                }
            }
            let root = 0.5 * (a + b) / 2f64.ln();
            if g < 0.0 && out.0.is_none() {
                out.0 = Some(root);
            } else if g > 0.0 && out.0.is_some() && out.1.is_none() {
                out.1 = Some(root);
            }
        }
        prev = g;
    }
    out
}

fn interval(v: &mut Vec<f64>) -> Value {
    if v.is_empty() {
        return Value::Null;
    }
    json!([pct(v, 0.025), pct(v, 0.975)])
}

fn main() {
    let dir = env::args().nth(1).expect("study directory");
    let dir = Path::new(&dir);
    let (ic_all, rho, meta) = load(dir);
    let rho_rungs = by_rung(&rho, |r| r.p);

    // Rho.
    let rho_pts = rung_means(&rho_rungs, &rho_rungs, |i| rho[i].gae);
    let (rc, rb, rse, _) = fit_points(&rho_pts).expect("rho fit");
    let s_rho_pts = rung_means(&rho_rungs, &rho_rungs, |i| rho[i].s);
    let (_, s_slope, s_slope_se, _) = fit_points(&s_rho_pts).expect("S_rho fit");
    let s_rho_mean = mean(&rho.iter().map(|r| r.s).collect::<Vec<_>>());

    let mut rng = StdRng::seed_from_u64(1);
    let mut out_arms = serde_json::Map::new();
    let mut rho_boot: Vec<(f64, f64, f64)> = Vec::with_capacity(REPLICATES);
    let mut rho_slope_boot = Vec::new();
    let mut s_rho_boot = Vec::new();
    for _ in 0..REPLICATES {
        let pick = resample(&rho_rungs, &mut rng);
        let pts = rung_means(&rho_rungs, &pick, |i| rho[i].gae);
        let (c, b, _, _) = fit_points(&pts).unwrap();
        let all: Vec<f64> = pick.values().flatten().map(|&i| rho[i].s).collect();
        let sm = mean(&all);
        rho_boot.push((c, b, sm));
        rho_slope_boot.push(b);
        s_rho_boot.push(sm);
    }

    for arm in IC_ARMS {
        let recs: Vec<IcRec> = ic_all
            .iter()
            .filter(|(a, _)| a == arm)
            .map(|(_, r)| r.clone())
            .collect();
        if recs.is_empty() {
            continue;
        }
        let rungs = by_rung(&recs, |r| r.p);
        let mut per_rung = Vec::new();
        for (p, idx) in &rungs {
            let mut m = serde_json::Map::new();
            m.insert("p".into(), json!(p));
            m.insert("n".into(), json!(idx.len()));
            for q in QUANTITIES {
                let v: Vec<f64> = idx.iter().map(|&i| recs[i].q[q]).collect();
                m.insert(q.into(), json!(mean(&v)));
                m.insert(format!("{q}_values"), json!(v));
            }
            let s: Vec<f64> = idx.iter().map(|&i| recs[i].s).collect();
            m.insert("s".into(), json!(mean(&s)));
            m.insert("s_values".into(), json!(s));
            m.insert(
                "workloads".into(),
                json!(idx.iter().map(|&i| recs[i].workload.clone()).collect::<Vec<_>>()),
            );
            for key in [
                "relations.relations",
                "relations.solver_calls",
                "relations.solver_ops",
                "relations.lp_merged_full",
                "relations.lp_full_direct",
                "relations.lp_partials_one_large",
                "relations.lp_partials_two_large",
                "relations.lp_discarded",
                "relations.lp_small_columns",
                "relations.lp_merge_row_ops",
                "linear_algebra.columns",
                "linear_algebra.wiedemann_dimension",
                "linear_algebra.wiedemann_touched_columns",
                "linear_algebra.wiedemann_attempts",
                "linear_algebra.wiedemann_nonzeros",
                "linear_algebra.row_ops",
                "linear_algebra.shadow_gauss_row_ops",
                "linear_algebra.shadow_gauss_pin_rank",
                "linear_algebra.shadow_gauss_agrees",
            ] {
                let v: Vec<f64> = idx
                    .iter()
                    .map(|&i| recs[i].native.get(key).copied().unwrap_or(0.0))
                    .collect();
                m.insert(key.into(), json!(v));
            }
            per_rung.push(Value::Object(m));
        }
        // Point fits.
        let mut fits = serde_json::Map::new();
        let mut point_coef: BTreeMap<&str, (f64, f64)> = BTreeMap::new();
        for q in QUANTITIES {
            let pts = rung_means(&rungs, &rungs, |i| recs[i].q[q]);
            let Some((c, b, se, _)) = fit_points(&pts) else {
                continue;
            };
            point_coef.insert(q, (c, b));
            let df = pts.iter().filter(|(_, y)| *y > 0.0).count() - 2;
            let local: Vec<Value> = pts
                .windows(2)
                .filter(|w| w[0].1 > 0.0 && w[1].1 > 0.0)
                .map(|w| {
                    json!({"p": [w[0].0.exp().round(), w[1].0.exp().round()],
                           "slope": (w[1].1.ln() - w[0].1.ln()) / (w[1].0 - w[0].0)})
                })
                .collect();
            fits.insert(
                q.into(),
                json!({"slope": b, "const_ln": c, "slope_se": se, "df": df,
                       "slope_t95": [b - t975(df) * se, b + t975(df) * se],
                       "local_slopes": local}),
            );
        }
        // Bootstrap: refit every quantity and the crossovers per replicate.
        let mut rng = StdRng::seed_from_u64(1);
        let mut slopes: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
        let mut cross: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        let mut cross_none: BTreeMap<String, u64> = BTreeMap::new();
        for (rep, &(rc_b, rb_b, sm_b)) in rho_boot.iter().enumerate() {
            let pick = resample(&rungs, &mut rng);
            let mut coef: BTreeMap<&str, (f64, f64)> = BTreeMap::new();
            for q in QUANTITIES {
                let pts = rung_means(&rungs, &pick, |i| recs[i].q[q]);
                if let Some((c, b, _, _)) = fit_points(&pts) {
                    slopes.entry(q).or_default().push(b);
                    coef.insert(q, (c, b));
                }
            }
            let _ = rep;
            for (variant, parts) in [
                ("wiedemann", vec!["relations", "linear_algebra", "other"]),
                ("dense_gauss", vec!["relations", "shadow_gauss", "other"]),
            ] {
                let ic: Vec<(f64, f64)> = parts.iter().filter_map(|q| coef.get(q).copied()).collect();
                for h in HS {
                    let rho_f = move |x: f64| (sm_b * ((4.0 * x).exp() / h).sqrt()).ln();
                    let key = format!("{variant}/h={h}");
                    match crossings(&ic, &rho_f).0 {
                        Some(v) => cross.entry(key).or_default().push(v),
                        None => *cross_none.entry(key).or_default() += 1,
                    }
                }
                let rho_pl = move |x: f64| rc_b + rb_b * x;
                let key = format!("{variant}/rho_power_law");
                match crossings(&ic, &rho_pl).0 {
                    Some(v) => cross.entry(key).or_default().push(v),
                    None => *cross_none.entry(key).or_default() += 1,
                }
            }
        }
        for q in QUANTITIES {
            if let (Some(f), Some(v)) = (fits.get_mut(q), slopes.get_mut(q)) {
                let boot = [pct(v, 0.025), pct(v, 0.975)];
                let t = f["slope_t95"].as_array().unwrap().clone();
                let union = [boot[0].min(num(&t[0])), boot[1].max(num(&t[1]))];
                f["slope_boot95"] = json!(boot);
                f["slope_ci95_union"] = json!(union);
            }
        }
        // Point crossovers.
        let mut crossovers = Vec::new();
        for (variant, parts) in [
            ("wiedemann", vec!["relations", "linear_algebra", "other"]),
            ("dense_gauss", vec!["relations", "shadow_gauss", "other"]),
        ] {
            let ic: Vec<(f64, f64)> = parts
                .iter()
                .filter_map(|q| point_coef.get(q).copied())
                .collect();
            let mut push = |label: String, rho_f: &dyn Fn(f64) -> f64| {
                let (first, back) = crossings(&ic, rho_f);
                let key = format!("{variant}/{label}");
                let mut v = cross.remove(&key).unwrap_or_default();
                let n_none = cross_none.get(&key).copied().unwrap_or(0);
                // A percentile interval only when at least 97.5 % of the
                // replicates cross; otherwise the count says it.
                let ci = if (n_none as f64) <= 0.025 * REPLICATES as f64 {
                    interval(&mut v)
                } else {
                    Value::Null
                };
                crossovers.push(json!({
                    "variant": variant, "rho": label,
                    "log2_p_star": first, "log2_group_star": first.map(|x| 4.0 * x),
                    "log2_p_back_above": back,
                    "log2_p_star_boot95": ci,
                    "log2_group_star_boot95": ci.as_array().map(|a| json!([4.0 * num(&a[0]), 4.0 * num(&a[1])])),
                    "replicates_without_crossing": n_none,
                    "top_rung_p": *rungs.keys().last().unwrap(),
                    "log2_reach_beyond_top_rung": first.map(|x| x - (*rungs.keys().last().unwrap() as f64).log2()),
                }));
            };
            for h in HS {
                push(format!("h={h}"), &move |x: f64| {
                    (s_rho_mean * ((4.0 * x).exp() / h).sqrt()).ln()
                });
            }
            push("rho_power_law".into(), &move |x: f64| rc + rb * x);
        }
        out_arms.insert(
            arm.into(),
            json!({"per_rung": per_rung, "fits": fits, "crossovers": crossovers}),
        );
    }

    // Paired: large-prime arms against the plain arm on the same workload.
    let mut paired = Vec::new();
    for lp in ["ic-d4-lp2", "ic-d4-lp1"] {
        for (a, r) in ic_all.iter().filter(|(a, _)| a == lp) {
            if let Some((_, w)) = ic_all
                .iter()
                .find(|(b, x)| b == "ic-d4-w" && x.workload == r.workload)
            {
                paired.push(json!({
                    "arm": a, "p": r.p, "workload": r.workload,
                    "s_lp": r.s, "s_plain": w.s, "ratio": r.s / w.s,
                    "la_lp": r.q["linear_algebra"], "la_plain": w.q["linear_algebra"],
                    "trials_lp": r.q["trials"], "trials_plain": w.q["trials"],
                }));
            }
        }
    }

    let rho_per_rung: Vec<Value> = rho_rungs
        .iter()
        .map(|(p, idx)| {
            json!({"p": p, "n": idx.len(),
                   "gae": mean(&idx.iter().map(|&i| rho[i].gae).collect::<Vec<_>>()),
                   "s": mean(&idx.iter().map(|&i| rho[i].s).collect::<Vec<_>>())})
        })
        .collect();
    let df = rho_pts.len() - 2;
    let summary = json!({
        "schema": "gaudry_quartic_d4_exponent_fit/v1",
        "replicates": REPLICATES,
        "records": meta,
        "rho": {
            "per_rung": rho_per_rung,
            "fit_gae": {"slope": rb, "const_ln": rc, "slope_t95": [rb - t975(df) * rse, rb + t975(df) * rse],
                        "slope_boot95": interval(&mut rho_slope_boot)},
            "s_rho_mean": s_rho_mean, "s_rho_boot95": interval(&mut s_rho_boot),
            "s_rho_slope": s_slope, "s_rho_slope_t95": [s_slope - t975(df) * s_slope_se, s_slope + t975(df) * s_slope_se],
            "records": rho.len(),
        },
        "arms": out_arms,
        "paired_large_prime": paired,
        "model_11_16_from_pr1863": {
            "formula": "S/rho = sqrt(h) (181 (p/4) C4 / (126 p^2) + r_inf/4), C4 = 8.34e7",
            "log2_group_star": {"h=1": [100.7, 101.2], "h=2": [103.0, 103.7], "h=4": [105.4, 106.5]},
            "source": "research/gaudry_quartic_d4_cost_20261014/raw/summary.json (phi = 0.73 .. 1)"
        }
    });
    let path = dir.join("raw/summary.json");
    fs::write(&path, serde_json::to_string_pretty(&summary).unwrap()).unwrap();
    println!("wrote {}", path.display());
}
