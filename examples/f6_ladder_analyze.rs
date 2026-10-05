//! Derive the F6-IC size-ladder decision from committed ecbench sessions,
//! exactly as `research/f6_ic_ecbench_ladder_20261005/PROTOCOL.md` and its
//! `AMENDMENT_1.md` fix it.
//!
//! Usage: f6_ladder_analyze CALIBRATION.json OUT.json SESSION_DIR...
//!
//! Reads only each session's `records.jsonl`. The paired unit is a
//! (target, round): within a round the three arms share one algorithm seed,
//! which is checked, as is the F4 and F6 arms' query stream. A target's
//! value is the mean over its measured rounds. For every size the analysis
//! reports both IC arms' lower-bound `S` against the floor and same-target
//! strong rho, the Boolean word XORs `W4` and `W6`, and F6-IC's geometric
//! additions `G6`. It fits the slope of `log2(W4/W6)` against `m` with a
//! within-size target bootstrap at the frozen seed, the per-relation
//! geometry check, the first-round-only sensitivity, and applies the frozen
//! rules. It computes nothing the protocol does not name.

use std::collections::BTreeMap;
use std::path::Path;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

const BOOTSTRAP_SAMPLES: usize = 10_000;
const BOOTSTRAP_SEED: u64 = 202_610_050_001;
const MIN_TARGETS: usize = 6;
const RHO_ARM: &str = "rho-strong";
const PHASES: [&str; 5] = [
    "factor_base",
    "oracle_setup",
    "relations",
    "linear_algebra",
    "verify",
];

/// The calibration entries the protocol's sensitivity row may use, by
/// field degree: `docs/ic/calibration.json`'s `koblitz` keys for the
/// `a = 0` curves `icv1-f2m13-t181-515ee569` and
/// `icv1-f2m31-tm90707-c95f16f5`.
const PINNED_WORD_XOR: [(u64, &str); 2] = [
    (13, "koblitz/K_0 / GF(2^13)"),
    (31, "koblitz/K_0 / GF(2^31)"),
];

#[derive(Clone, Debug)]
struct Run {
    arm: String,
    workload: String,
    round: u64,
    seed: String,
    m: u64,
    r: f64,
    slug: String,
    verified: bool,
    status: String,
    s: Option<f64>,
    total_gae: Option<f64>,
    word_xors: Option<u64>,
    reductions: Option<u64>,
    geometry: Option<u64>,
    relations: Option<u64>,
    trials: Option<u64>,
    phases: BTreeMap<String, f64>,
    isolation: String,
    solve_wall_ns: Option<u64>,
}

fn read_runs(dir: &Path) -> Result<Vec<Run>, String> {
    let path = dir.join("records.jsonl");
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut out = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let r: Value =
            serde_json::from_str(line).map_err(|e| format!("{}: {e}", path.display()))?;
        if r["warmup"].as_bool().unwrap_or(false) {
            continue;
        }
        let curve = &r["workload"]["curve"];
        let status = r["outcome"]["status"].as_str().unwrap_or("?").to_string();
        let phases = r["phases"]
            .as_array()
            .map(|ps| {
                ps.iter()
                    .filter_map(|p| Some((p["name"].as_str()?.to_string(), p["gae"].as_f64()?)))
                    .collect()
            })
            .unwrap_or_default();
        out.push(Run {
            arm: r["arm"].as_str().unwrap_or("?").to_string(),
            workload: r["workload"]["workload_id"]
                .as_str()
                .unwrap_or("?")
                .to_string(),
            round: r["round"].as_u64().ok_or("record without round")?,
            seed: r["algorithm_seed"].as_str().unwrap_or("?").to_string(),
            m: curve["field_degree"]
                .as_u64()
                .ok_or("record without field_degree")?,
            r: curve["r"].as_f64().ok_or("record without r")?,
            slug: curve["slug"].as_str().unwrap_or("?").to_string(),
            verified: status == "verified",
            status,
            s: r["cost"]["s"].as_f64(),
            total_gae: r["cost"]["total_gae"].as_f64(),
            word_xors: r["solver"]["ops"].as_u64(),
            reductions: r["solver"]["extra"]["reductions"].as_u64(),
            geometry: r["solver"]["extra"]["geometric_group_additions"].as_u64(),
            relations: r["counters"]["relations_found"].as_u64(),
            trials: r["counters"]["targets_tried"].as_u64(),
            phases,
            isolation: r["isolation"]["level"].as_str().unwrap_or("?").to_string(),
            solve_wall_ns: r["time"]["solve_wall_ns"].as_u64(),
        });
    }
    Ok(out)
}

fn median(mut v: Vec<f64>) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    Some(if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    })
}

fn mean(v: &[f64]) -> Option<f64> {
    (!v.is_empty()).then(|| v.iter().sum::<f64>() / v.len() as f64)
}

fn min_max(v: &[f64]) -> Option<(f64, f64)> {
    let lo = v.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    (!v.is_empty()).then_some((lo, hi))
}

/// Ordinary least squares slope and intercept of `y` on `x`.
fn ols(points: &[(f64, f64)]) -> Option<(f64, f64)> {
    if points.len() < 2 {
        return None;
    }
    let n = points.len() as f64;
    let mx = points.iter().map(|p| p.0).sum::<f64>() / n;
    let my = points.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx: f64 = points.iter().map(|p| (p.0 - mx).powi(2)).sum();
    if sxx == 0.0 {
        return None;
    }
    let sxy: f64 = points.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let b = sxy / sxx;
    Some((b, my - b * mx))
}

/// The slope's 95% percentile interval, resampling targets within each size.
fn bootstrap(points: &[(f64, f64)]) -> (usize, Option<f64>, Option<f64>) {
    let mut by_size: BTreeMap<u64, Vec<(f64, f64)>> = BTreeMap::new();
    for p in points {
        by_size.entry(p.0 as u64).or_default().push(*p);
    }
    let mut rng = StdRng::seed_from_u64(BOOTSTRAP_SEED);
    let mut slopes = Vec::with_capacity(BOOTSTRAP_SAMPLES);
    for _ in 0..BOOTSTRAP_SAMPLES {
        let mut sample = Vec::with_capacity(points.len());
        for pts in by_size.values() {
            for _ in 0..pts.len() {
                sample.push(pts[rng.gen_range(0..pts.len())]);
            }
        }
        if let Some((b, _)) = ols(&sample) {
            slopes.push(b);
        }
    }
    slopes.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let pct = |q: f64| -> Option<f64> {
        (!slopes.is_empty()).then(|| slopes[(q * (slopes.len() - 1) as f64).round() as usize])
    };
    (slopes.len(), pct(0.025), pct(0.975))
}

/// One (target, round): the two IC arms and rho on one seed.
struct Pair {
    round: u64,
    f4: Run,
    f6: Run,
    rho: Option<Run>,
}

impl Pair {
    fn log2_ratio(&self) -> Option<f64> {
        Some((self.f4.word_xors? as f64 / self.f6.word_xors? as f64).log2())
    }
}

struct Target {
    m: u64,
    workload: String,
    pairs: Vec<Pair>,
}

impl Target {
    /// The mean over this target's rounds of a per-pair quantity.
    fn mean(&self, f: impl Fn(&Pair) -> Option<f64>) -> Option<f64> {
        let v: Vec<f64> = self.pairs.iter().filter_map(&f).collect();
        (v.len() == self.pairs.len()).then(|| mean(&v)).flatten()
    }

    /// The same quantity on the first measured round only.
    fn first(&self, f: impl Fn(&Pair) -> Option<f64>) -> Option<f64> {
        self.pairs.iter().min_by_key(|p| p.round).and_then(f)
    }
}

fn arm_of<'p>(p: &'p Pair, arm: &str) -> &'p Run {
    if arm == "ic-f4" {
        &p.f4
    } else {
        &p.f6
    }
}

fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        return Err("usage: f6_ladder_analyze CALIBRATION.json OUT.json SESSION_DIR...".into());
    }
    let calibration: Value = serde_json::from_str(
        &std::fs::read_to_string(&args[1]).map_err(|e| format!("{}: {e}", args[1]))?,
    )
    .map_err(|e| e.to_string())?;
    let mut all = Vec::new();
    let mut sessions = Vec::new();
    for d in &args[3..] {
        let runs = read_runs(Path::new(d))?;
        sessions.push(json!({"dir": d, "measured_records": runs.len()}));
        all.extend(runs);
    }

    // (m, workload, round) -> arm -> run.
    let mut cells: BTreeMap<(u64, String, u64), BTreeMap<String, Run>> = BTreeMap::new();
    for r in &all {
        cells
            .entry((r.m, r.workload.clone(), r.round))
            .or_default()
            .insert(r.arm.clone(), r.clone());
    }

    let mut seed_mismatches = Vec::new();
    let mut divergent = Vec::new();
    let mut incomplete = Vec::new();
    let mut by_target: BTreeMap<(u64, String), (Vec<Pair>, bool)> = BTreeMap::new();
    for ((m, w, round), arms) in &cells {
        let entry = by_target
            .entry((*m, w.clone()))
            .or_insert_with(|| (Vec::new(), true));
        let seeds: Vec<&String> = arms.values().map(|r| &r.seed).collect();
        if seeds.windows(2).any(|s| s[0] != s[1]) {
            seed_mismatches.push(json!({"m": m, "workload": w, "round": round}));
        }
        match (arms.get("ic-f4"), arms.get("ic-f6")) {
            (Some(f4), Some(f6)) if f4.verified && f6.verified => {
                if f4.trials != f6.trials || f4.relations != f6.relations {
                    divergent.push(json!({"m": m, "workload": w, "round": round,
                        "f4": {"trials": f4.trials, "relations": f4.relations},
                        "f6": {"trials": f6.trials, "relations": f6.relations}}));
                }
                entry.0.push(Pair {
                    round: *round,
                    f4: f4.clone(),
                    f6: f6.clone(),
                    rho: arms.get(RHO_ARM).filter(|r| r.verified).cloned(),
                });
            }
            (f4, f6) => {
                entry.1 = false;
                incomplete.push(json!({"m": m, "workload": w, "round": round,
                    "f4": f4.map(|r| r.status.clone()), "f6": f6.map(|r| r.status.clone())}));
            }
        }
    }
    // A target completes when every measured round verified in both IC arms.
    let targets: Vec<Target> = by_target
        .into_iter()
        .filter(|(_, (pairs, complete))| *complete && !pairs.is_empty())
        .map(|((m, workload), (pairs, _))| Target { m, workload, pairs })
        .collect();

    let sizes: Vec<u64> = {
        let mut v: Vec<u64> = all.iter().map(|r| r.m).collect();
        v.sort();
        v.dedup();
        v
    };

    let mut per_size = Vec::new();
    let mut completed_by_size: BTreeMap<u64, usize> = BTreeMap::new();
    for &m in &sizes {
        let ts: Vec<&Target> = targets.iter().filter(|t| t.m == m).collect();
        completed_by_size.insert(m, ts.len());
        let any = all.iter().find(|r| r.m == m).unwrap();
        let floor = (std::f64::consts::PI / (4.0 * m as f64)).sqrt();
        let med =
            |f: &dyn Fn(&Target) -> Option<f64>| median(ts.iter().filter_map(|t| f(t)).collect());
        let mut arms = serde_json::Map::new();
        for arm in ["ic-f4", "ic-f6"] {
            let s_of = |t: &Target| t.mean(|p| arm_of(p, arm).s);
            let vs_rho: Vec<f64> = ts
                .iter()
                .filter_map(|t| t.mean(|p| Some(arm_of(p, arm).s? / p.rho.as_ref()?.s?)))
                .collect();
            let mut phase_median = serde_json::Map::new();
            for name in PHASES {
                phase_median.insert(
                    name.into(),
                    json!(med(
                        &|t| t.mean(|p| arm_of(p, arm).phases.get(name).copied())
                    )),
                );
            }
            let sensitivity =
                PINNED_WORD_XOR
                    .iter()
                    .find(|(deg, _)| *deg == m)
                    .and_then(|(_, key)| {
                        let ratio = calibration["instances"][key]["ns_per_word_xor"].as_f64()?;
                        let priced = med(&|t| {
                            t.mean(|p| {
                                let a = arm_of(p, arm);
                                Some((a.total_gae? + a.word_xors? as f64 * ratio) / a.r.sqrt())
                            })
                        });
                        Some(json!({"calibration_key": key,
                        "ns_per_word_xor_over_ns_per_add": ratio,
                        "median_s_with_word_xors_priced": priced}))
                    });
            arms.insert(
                arm.into(),
                json!({
                    "median_s_lower": med(&s_of),
                    "median_s_lower_over_floor": med(&|t| Some(s_of(t)? / floor)),
                    "median_s_lower_over_rho": median(vs_rho.clone()),
                    "min_max_s_lower_over_rho": min_max(&vs_rho),
                    "median_word_xors": med(&|t| t.mean(|p| Some(arm_of(p, arm).word_xors? as f64))),
                    "median_reductions": med(&|t| t.mean(|p| Some(arm_of(p, arm).reductions? as f64))),
                    "median_geometric_additions": med(&|t| t.mean(|p| Some(arm_of(p, arm).geometry? as f64))),
                    "median_relations": med(&|t| t.mean(|p| Some(arm_of(p, arm).relations? as f64))),
                    "median_phase_gae": phase_median,
                    "median_solve_wall_ns_practicality_only": med(&|t| t.mean(|p| Some(arm_of(p, arm).solve_wall_ns? as f64))),
                    "pinned_sensitivity": sensitivity,
                }),
            );
        }
        let ys: Vec<f64> = ts.iter().filter_map(|t| t.mean(Pair::log2_ratio)).collect();
        let mut levels: BTreeMap<String, usize> = BTreeMap::new();
        for r in all.iter().filter(|r| r.m == m) {
            *levels.entry(r.isolation.clone()).or_default() += 1;
        }
        per_size.push(json!({
            "m": m,
            "slug": any.slug,
            "r": any.r,
            "s_floor": floor,
            "rounds_per_target": ts.first().map(|t| t.pairs.len()),
            "targets_completed_both_ic_arms": ts.len(),
            "rho_median_s": med(&|t| t.mean(|p| p.rho.as_ref()?.s)),
            "arms": arms,
            "median_log2_w4_over_w6": median(ys.clone()),
            "min_max_log2_w4_over_w6": min_max(&ys),
            "abandon_rule_all_w6_ge_w4": (m == 19).then(|| ts.iter().all(|t| {
                t.pairs.iter().all(|p| p.f6.word_xors.unwrap_or(u64::MAX) >= p.f4.word_xors.unwrap_or(0))
            })),
            "isolation_levels": levels,
        }));
    }

    // Primary fit (per-target mean over rounds) and first-round sensitivity.
    let fit_of = |f: &dyn Fn(&Target) -> Option<f64>| -> Value {
        let points: Vec<(f64, f64)> = targets
            .iter()
            .filter_map(|t| Some((t.m as f64, f(t)?)))
            .collect();
        let fit = ols(&points);
        let (n, lo, hi) = bootstrap(&points);
        json!({"points": points.len(), "slope_per_unit_m": fit.map(|f| f.0),
               "intercept": fit.map(|f| f.1),
               "bootstrap": {"samples": n, "seed": BOOTSTRAP_SEED, "ci95": [lo, hi],
                             "scheme": "targets resampled within each size"}})
    };
    let primary = fit_of(&|t| t.mean(Pair::log2_ratio));
    let first_round = fit_of(&|t| t.first(Pair::log2_ratio));
    let slope_of = |f: &dyn Fn(&Target) -> Option<f64>| -> Option<f64> {
        let pts: Vec<(f64, f64)> = targets
            .iter()
            .filter_map(|t| Some((t.m as f64, f(t)?.log2())))
            .collect();
        ols(&pts).map(|(b, _)| b)
    };
    let geometry_slope =
        slope_of(&|t| t.mean(|p| Some(p.f6.geometry? as f64 / p.f6.relations?.max(1) as f64)));
    let w4_slope =
        slope_of(&|t| t.mean(|p| Some(p.f4.word_xors? as f64 / p.f4.relations?.max(1) as f64)));

    let decide = |fit: &Value| -> &'static str {
        let undecidable = sizes.len() < 4
            || sizes
                .iter()
                .any(|m| completed_by_size.get(m).copied().unwrap_or(0) < MIN_TARGETS);
        let lo = fit["bootstrap"]["ci95"][0].as_f64();
        if undecidable {
            "not decidable: fewer than six targets completed at some size, or fewer than four sizes"
        } else if lo.is_some_and(|lo| lo > 0.0)
            && matches!((geometry_slope, w4_slope), (Some(g), Some(w)) if g <= w)
        {
            "supported: a stage-level exponent lead (reductions unpriced; no whole-method claim)"
        } else {
            "rejected: the F6-IC gain is a constant factor; class engineering; close the F6 thread"
        }
    };
    let decision = decide(&primary);
    let sensitivity_decision = decide(&first_round);

    let rows: Vec<Value> = targets
        .iter()
        .map(|t| {
            json!({
                "m": t.m,
                "workload": t.workload,
                "rounds": t.pairs.len(),
                "mean_log2_w4_over_w6": t.mean(Pair::log2_ratio),
                "first_round_log2_w4_over_w6": t.first(Pair::log2_ratio),
                "mean_w4_word_xors": t.mean(|p| Some(p.f4.word_xors? as f64)),
                "mean_w6_word_xors": t.mean(|p| Some(p.f6.word_xors? as f64)),
                "mean_f4_reductions": t.mean(|p| Some(p.f4.reductions? as f64)),
                "mean_f6_reductions": t.mean(|p| Some(p.f6.reductions? as f64)),
                "mean_g6_geometric_additions": t.mean(|p| Some(p.f6.geometry? as f64)),
                "mean_relations": t.mean(|p| Some(p.f4.relations? as f64)),
                "mean_trials": t.mean(|p| Some(p.f4.trials? as f64)),
                "mean_s_lower_f4": t.mean(|p| p.f4.s),
                "mean_s_lower_f6": t.mean(|p| p.f6.s),
                "mean_s_rho": t.mean(|p| p.rho.as_ref()?.s),
            })
        })
        .collect();

    let out = json!({
        "schema": "f6-ic-ladder-analysis/v2",
        "protocol": ["research/f6_ic_ecbench_ladder_20261005/PROTOCOL.md",
                     "research/f6_ic_ecbench_ladder_20261005/AMENDMENT_1.md"],
        "sessions": sessions,
        "seed_mismatches": seed_mismatches,
        "incomplete_pairs": incomplete,
        "divergent_query_streams": divergent,
        "targets": rows,
        "per_size": per_size,
        "fit": {
            "statistic": "log2(W4/W6), W = Boolean solver word XORs over all PDP calls; per-target mean over rounds",
            "primary": primary,
            "first_round_only_sensitivity": first_round,
            "slope_log2_geometry_per_relation_f6": geometry_slope,
            "slope_log2_word_xors_per_relation_f4": w4_slope,
            "slope_log2_s_lower_ic_f4": slope_of(&|t| t.mean(|p| p.f4.s)),
            "slope_log2_s_lower_ic_f6": slope_of(&|t| t.mean(|p| p.f6.s)),
            "slope_log2_s_rho": slope_of(&|t| t.mean(|p| p.rho.as_ref()?.s)),
        },
        "decision": decision,
        "first_round_only_decision": sensitivity_decision,
        "decisions_agree": decision == sensitivity_decision,
    });
    let text = serde_json::to_string_pretty(&out).map_err(|e| e.to_string())?;
    std::fs::write(&args[2], format!("{text}\n")).map_err(|e| format!("{}: {e}", args[2]))?;
    println!("{text}");
    Ok(())
}
