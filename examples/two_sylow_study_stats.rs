//! Derived statistics of `research/pkm_isogeny_tower_p192_20261015`
//! (PROTOCOL.md: H1, H2, S1–S4, the draws-per-success table and the
//! `t = 48` extrapolation), from the recorded cell outputs.
//!
//! ```text
//! cargo run --release --example two_sylow_study_stats -- DATA_DIR > derived.json
//! ```

use serde_json::{json, Value};

fn load(p: &str) -> Value {
    serde_json::from_str(&std::fs::read_to_string(p).unwrap_or_else(|e| panic!("{p}: {e}")))
        .unwrap()
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

/// Two-sided Poisson interval `[lo, hi]` with mass `1 − alpha` around mean
/// `lambda`: exact quantiles up to `lambda = 500` (`e^−λ` stays representable), normal above.
fn poisson_interval(lambda: f64, alpha: f64) -> (u64, u64) {
    if lambda > 500.0 {
        let z = if alpha < 0.002 {
            3.290_526_731
        } else {
            1.959_963_985
        };
        let s = lambda.sqrt();
        return (
            (lambda - z * s).floor().max(0.0) as u64,
            (lambda + z * s).ceil() as u64,
        );
    }
    let mut pmf = (-lambda).exp();
    let mut cdf = pmf;
    let mut k = 0u64;
    let mut lo = None;
    while cdf < 1.0 - alpha / 2.0 {
        if lo.is_none() && cdf >= alpha / 2.0 {
            lo = Some(k);
        }
        k += 1;
        pmf *= lambda / k as f64;
        cdf += pmf;
    }
    (lo.unwrap_or(0).min(k), k)
}

fn main() {
    let dir = std::env::args().nth(1).expect("DATA_DIR");
    // ── C1: per-draw cost ──
    let mut c1 = serde_json::Map::new();
    let mut med = std::collections::HashMap::new();
    for fam in ["z2xz8", "z2xz2"] {
        for arm in ["baseline", "fast-chi", "fast"] {
            let mut v: Vec<f64> = (1..=7)
                .map(|r| {
                    load(&format!("{dir}/c1/{fam}-{arm}-rep{r}.json"))["us_per_draw"]
                        .as_f64()
                        .unwrap()
                })
                .collect();
            let (mn, mx) = (
                v.iter().cloned().fold(f64::MAX, f64::min),
                v.iter().cloned().fold(0.0, f64::max),
            );
            let m = median(&mut v);
            med.insert((fam, arm), (m, mn, mx));
            c1.insert(
                format!("{fam}/{arm}"),
                json!({"median_us": m, "min_us": mn, "max_us": mx, "reps": 7, "draws_per_rep": 65536}),
            );
        }
    }
    let r = |fam, a, b| med[&(fam, a)].0 / med[&(fam, b)].0;
    let s1 = r("z2xz8", "fast", "baseline");
    let s2 = r("z2xz8", "fast", "fast-chi");
    let s3 = med[&("z2xz8", "fast")].0;
    let s4 = r("z2xz2", "fast", "baseline");
    let claims = json!({
        "S1": {"ratio": s1, "threshold": 0.40, "pass": s1 <= 0.40,
               "spread_overlap": med[&("z2xz8","fast")].2 >= med[&("z2xz8","baseline")].1},
        "S2": {"ratio": s2, "threshold": 0.85, "pass": s2 <= 0.85,
               "spread_overlap": med[&("z2xz8","fast")].2 >= med[&("z2xz8","fast-chi")].1},
        "S3": {"median_us": s3, "threshold": 12.0, "pass": s3 <= 12.0},
        "S4": {"ratio": s4, "threshold": 0.65, "pass": s4 <= 0.65,
               "spread_overlap": med[&("z2xz2","fast")].2 >= med[&("z2xz2","baseline")].1},
    });
    // ── C3/C4: depth distributions against the Haar model ──
    let mut dist = serde_json::Map::new();
    let mut verdicts = serde_json::Map::new();
    for (cell, fam) in [("c3", "z2xz2"), ("c4", "z2xz8")] {
        let s = load(&format!("{dir}/{cell}/summary.json"));
        let hist: Vec<u64> = s["hist_depth"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap())
            .collect();
        let n: u64 = hist.iter().sum::<u64>() - hist[0];
        let singular = hist[0];
        let maxb = hist.iter().rposition(|&c| c > 0).unwrap();
        let mut rows = Vec::new();
        let mut consistent = true;
        let tmin = if fam == "z2xz2" { 1 } else { 4 };
        for t in tmin..=maxb.max(tmin) {
            let ge: u64 = hist[t..].iter().sum();
            let p_model = match (fam, t) {
                ("z2xz2", 1) => 1.0,
                ("z2xz2", t) => 3.0 * 2f64.powi(-(t as i32)),
                (_, t) => 2f64.powi(3 - t as i32),
            };
            let lambda = n as f64 * p_model;
            let (lo, hi) = poisson_interval(lambda, 0.001);
            let (lo95, hi95) = poisson_interval(ge as f64, 0.05);
            let tested = lambda >= 20.0 && !(fam == "z2xz2" && t == 1);
            let inside = ge >= lo && ge <= hi;
            if tested && !inside {
                consistent = false;
            }
            // Draws per success: 1 / P̂ with a 95 % interval from the count.
            let dps = if ge > 0 {
                n as f64 / ge as f64
            } else {
                f64::NAN
            };
            rows.push(json!({
                "t": t, "count_ge": ge, "frac_ge": ge as f64 / n as f64,
                "model_frac_ge": p_model, "model_count": lambda,
                "model_interval_999": [lo, hi], "tested": tested, "inside": inside,
                "draws_per_success": dps,
                "draws_per_success_95": [n as f64 / hi95.max(1) as f64,
                                         if lo95 > 0 { n as f64 / lo95 as f64 } else { f64::INFINITY }],
                "model_draws_per_success": 1.0 / p_model,
                "ratio_obs_model": ge as f64 / lambda,
            }));
        }
        // P(b = 1) for z2xz2.
        if fam == "z2xz2" {
            let lambda = n as f64 / 4.0;
            let (lo, hi) = poisson_interval(lambda, 0.001);
            let inside = hist[1] >= lo && hist[1] <= hi;
            consistent &= inside;
            verdicts.insert("H1_b_eq_1".into(), json!({"count": hist[1], "model": lambda, "interval_999": [lo, hi], "inside": inside}));
        }
        // Gaps between consecutive draws with b ≥ t (record-min 16).
        let recs: Vec<(u64, u64)> = s["records"]
            .as_array()
            .unwrap()
            .iter()
            .map(|r| (r[0].as_u64().unwrap(), r[1].as_u64().unwrap()))
            .collect();
        let mut gaps = serde_json::Map::new();
        for t in [16u64, 20, 24, 28] {
            let seeds: Vec<u64> = recs.iter().filter(|r| r.1 >= t).map(|r| r.0).collect();
            if seeds.len() >= 2 {
                let g: Vec<f64> = seeds.windows(2).map(|w| (w[1] - w[0]) as f64).collect();
                let mean = g.iter().sum::<f64>() / g.len() as f64;
                let var = g.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
                    / (g.len() as f64 - 1.0).max(1.0);
                gaps.insert(t.to_string(), json!({"successes": seeds.len(), "first_seed": seeds[0],
                    "mean_gap": mean, "se": (var / g.len() as f64).sqrt(),
                    "model": if fam == "z2xz8" { 2f64.powi(t as i32 - 3) } else { 2f64.powi(t as i32) / 3.0 }}));
            } else {
                gaps.insert(
                    t.to_string(),
                    json!({"successes": seeds.len(), "first_seed": seeds.first()}),
                );
            }
        }
        let deepest: Vec<Value> = recs
            .iter()
            .filter(|r| r.1 as usize == maxb)
            .map(|r| json!([r.0, r.1]))
            .collect();
        dist.insert(
            fam.into(),
            json!({"cell": cell, "draws": n, "singular": singular,
            "elapsed_s": s["elapsed_s"], "us_per_draw_long_run": s["us_per_draw"],
            "hist_depth": hist, "max_b": maxb, "deepest": deepest,
            "rows": rows, "gaps": gaps, "model_consistent": consistent}),
        );
        verdicts.insert(
            if fam == "z2xz2" {
                "H1".into()
            } else {
                "H2".into()
            },
            json!({"consistent": consistent}),
        );
    }
    // ── Derived: t = 48 on z2xz8 ──
    let (m, mn, mx) = med[&("z2xz8", "fast")];
    let draws48 = 2f64.powi(45);
    let ch = |us: f64| draws48 * us / 3.6e9;
    let t48 = json!({
        "label": "derived (H2 model draws x C1 fast median)",
        "draws": draws48, "core_hours": ch(m), "core_hours_min": ch(mn), "core_hours_max": ch(mx),
        "core_years": ch(m) / 8766.0,
        "z2xz2_draws": 2f64.powi(48) / 3.0,
        "z2xz2_core_hours": 2f64.powi(48) / 3.0 * med[&("z2xz2", "fast")].0 / 3.6e9,
        "card_previous_est_core_hours": 39094,
    });
    let out = json!({"schema": "pkm_isogeny_tower_p192.derived/v1", "c1": c1, "claims": claims,
        "distributions": dist, "verdicts": verdicts, "t48": t48});
    println!("{}", serde_json::to_string_pretty(&out).unwrap());
}
