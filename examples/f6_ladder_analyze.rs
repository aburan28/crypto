//! Derive the F6-IC size-ladder decision from committed ecbench sessions,
//! exactly as `research/f6_ic_ecbench_ladder_20261005/PROTOCOL.md` fixes it.
//!
//! Usage: f6_ladder_analyze CALIBRATION.json OUT.json SESSION_DIR...
//!
//! Reads only each session's `records.jsonl`. For every size it reports the
//! IC arms' lower-bound `S` against the floor and the same-target strong-rho
//! `S`, the Boolean word XORs `W4`, `W6`, and F6-IC's geometric additions
//! `G6`. It fits the protocol's slope of `log2(W4/W6)` against `m`, with a
//! within-size bootstrap at the frozen seed, and applies the frozen rules.
//! It computes nothing the protocol does not name.

use std::collections::BTreeMap;
use std::path::Path;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

const BOOTSTRAP_SAMPLES: usize = 10_000;
const BOOTSTRAP_SEED: u64 = 202_610_050_001;
const MIN_TARGETS: usize = 6;
const IC_ARMS: [&str; 2] = ["ic-f4", "ic-f6"];
const RHO_ARM: &str = "rho-strong";

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
}

fn u(v: &Value) -> Option<u64> {
    v.as_u64()
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
            m: u(&curve["field_degree"]).ok_or("record without field_degree")?,
            r: curve["r"].as_f64().ok_or("record without r")?,
            slug: curve["slug"].as_str().unwrap_or("?").to_string(),
            verified: status == "verified",
            status,
            s: r["cost"]["s"].as_f64(),
            total_gae: r["cost"]["total_gae"].as_f64(),
            word_xors: u(&r["solver"]["ops"]),
            reductions: u(&r["solver"]["extra"]["reductions"]),
            geometry: u(&r["solver"]["extra"]["geometric_group_additions"]),
            relations: u(&r["counters"]["relations_found"]),
            trials: u(&r["counters"]["targets_tried"]),
            phases,
            isolation: r["isolation"]["level"].as_str().unwrap_or("?").to_string(),
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

fn min_max(v: &[f64]) -> Option<(f64, f64)> {
    let lo = v.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    (!v.is_empty()).then_some((lo, hi))
}

/// Ordinary least squares slope and intercept of `y` on `x`.
fn ols(points: &[(f64, f64)]) -> Option<(f64, f64)> {
    let n = points.len() as f64;
    if points.len() < 2 {
        return None;
    }
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

/// Every measured round of a deterministic arm must carry the same counts.
fn check_rounds(runs: &[&Run]) -> Result<(), String> {
    let first = runs[0];
    for r in &runs[1..] {
        if r.status != first.status
            || r.word_xors != first.word_xors
            || r.geometry != first.geometry
            || r.relations != first.relations
            || r.total_gae.map(f64::to_bits) != first.total_gae.map(f64::to_bits)
        {
            return Err(format!(
                "{} {}: measured rounds disagree on counts",
                first.arm, first.workload
            ));
        }
    }
    Ok(())
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

    // (m, workload, arm) -> first measured run, after checking rounds agree.
    let mut grouped: BTreeMap<(u64, String, String), Vec<&Run>> = BTreeMap::new();
    for r in &all {
        grouped
            .entry((r.m, r.workload.clone(), r.arm.clone()))
            .or_default()
            .push(r);
    }
    let mut cell: BTreeMap<(u64, String, String), &Run> = BTreeMap::new();
    let mut round_disagreements = Vec::new();
    for (k, rs) in &grouped {
        if let Err(e) = check_rounds(rs) {
            round_disagreements.push(e);
        }
        cell.insert(k.clone(), rs[0]);
    }

    let sizes: Vec<u64> = {
        let mut v: Vec<u64> = all.iter().map(|r| r.m).collect();
        v.sort();
        v.dedup();
        v
    };

    // Per-target paired rows.
    struct Target {
        m: u64,
        workload: String,
        f4: Run,
        f6: Run,
        rho: Option<Run>,
    }
    let mut targets: Vec<Target> = Vec::new();
    let mut incomplete = Vec::new();
    let mut divergent = Vec::new();
    for (k, f4) in &cell {
        if k.2 != "ic-f4" {
            continue;
        }
        let f6 = cell.get(&(k.0, k.1.clone(), "ic-f6".to_string()));
        let rho = cell
            .get(&(k.0, k.1.clone(), RHO_ARM.to_string()))
            .map(|r| (*r).clone());
        match f6 {
            Some(f6) if f4.verified && f6.verified => {
                if f4.trials != f6.trials || f4.relations != f6.relations {
                    divergent.push(json!({"m": k.0, "workload": k.1,
                        "f4": {"trials": f4.trials, "relations": f4.relations},
                        "f6": {"trials": f6.trials, "relations": f6.relations}}));
                }
                targets.push(Target {
                    m: k.0,
                    workload: k.1.clone(),
                    f4: (*f4).clone(),
                    f6: (*f6).clone(),
                    rho,
                });
            }
            _ => incomplete.push(json!({"m": k.0, "workload": k.1,
                "f4": f4.status, "f6": f6.map(|r| r.status.clone())})),
        }
    }

    let log2_ratio = |t: &Target| -> Option<f64> {
        Some((t.f4.word_xors? as f64 / t.f6.word_xors? as f64).log2())
    };

    // Per size.
    let mut per_size = Vec::new();
    let mut completed_by_size: BTreeMap<u64, usize> = BTreeMap::new();
    for &m in &sizes {
        let ts: Vec<&Target> = targets.iter().filter(|t| t.m == m).collect();
        completed_by_size.insert(m, ts.len());
        let any = all.iter().find(|r| r.m == m).unwrap();
        let floor = (std::f64::consts::PI / (4.0 * m as f64)).sqrt();
        let rho_s: Vec<f64> = cell
            .iter()
            .filter(|(k, r)| k.0 == m && k.2 == RHO_ARM && r.verified)
            .filter_map(|(_, r)| r.s)
            .collect();
        let mut arms = serde_json::Map::new();
        for arm in IC_ARMS {
            let pick = |t: &Target| {
                if arm == "ic-f4" {
                    t.f4.clone()
                } else {
                    t.f6.clone()
                }
            };
            let s: Vec<f64> = ts.iter().filter_map(|t| pick(t).s).collect();
            let vs_rho: Vec<f64> = ts
                .iter()
                .filter_map(|t| Some(pick(t).s? / t.rho.as_ref().filter(|r| r.verified)?.s?))
                .collect();
            let mut phase_median = serde_json::Map::new();
            for name in [
                "factor_base",
                "oracle_setup",
                "relations",
                "linear_algebra",
                "verify",
            ] {
                let v: Vec<f64> = ts
                    .iter()
                    .filter_map(|t| pick(t).phases.get(name).copied())
                    .collect();
                phase_median.insert(name.into(), json!(median(v)));
            }
            let sensitivity =
                PINNED_WORD_XOR
                    .iter()
                    .find(|(deg, _)| *deg == m)
                    .and_then(|(_, key)| {
                        let ratio = calibration["instances"][key]["ns_per_word_xor"].as_f64()?;
                        let priced: Vec<f64> = ts
                            .iter()
                            .filter_map(|t| {
                                let p = pick(t);
                                Some((p.total_gae? + p.word_xors? as f64 * ratio) / p.r.sqrt())
                            })
                            .collect();
                        Some(
                            json!({"calibration_key": key, "ns_per_word_xor_over_ns_per_add": ratio,
                        "median_s_with_word_xors_priced": median(priced)}),
                        )
                    });
            arms.insert(
                arm.into(),
                json!({
                    "median_s_lower": median(s.clone()),
                    "median_s_lower_over_floor": median(s.iter().map(|x| x / floor).collect()),
                    "median_s_lower_over_rho": median(vs_rho.clone()),
                    "min_max_s_lower_over_rho": min_max(&vs_rho),
                    "median_word_xors": median(ts.iter().filter_map(|t| pick(t).word_xors.map(|x| x as f64)).collect()),
                    "median_reductions": median(ts.iter().filter_map(|t| pick(t).reductions.map(|x| x as f64)).collect()),
                    "median_geometric_additions": median(ts.iter().filter_map(|t| pick(t).geometry.map(|x| x as f64)).collect()),
                    "median_relations": median(ts.iter().filter_map(|t| pick(t).relations.map(|x| x as f64)).collect()),
                    "median_phase_gae": phase_median,
                    "pinned_sensitivity": sensitivity,
                }),
            );
        }
        let ys: Vec<f64> = ts.iter().filter_map(|t| log2_ratio(t)).collect();
        let mut levels: BTreeMap<String, usize> = BTreeMap::new();
        for r in all.iter().filter(|r| r.m == m) {
            *levels.entry(r.isolation.clone()).or_default() += 1;
        }
        per_size.push(json!({
            "m": m,
            "slug": any.slug,
            "r": any.r,
            "s_floor": floor,
            "targets_completed_both_ic_arms": ts.len(),
            "rho_median_s": median(rho_s),
            "arms": arms,
            "median_log2_w4_over_w6": median(ys.clone()),
            "min_max_log2_w4_over_w6": min_max(&ys),
            "abandon_rule_all_w6_ge_w4": (m == 19).then(|| ts.iter().all(|t|
                t.f6.word_xors.unwrap_or(u64::MAX) >= t.f4.word_xors.unwrap_or(0))),
            "isolation_levels": levels,
        }));
    }

    // The slope and its within-size bootstrap.
    let points: Vec<(f64, f64)> = targets
        .iter()
        .filter_map(|t| Some((t.m as f64, log2_ratio(t)?)))
        .collect();
    let fit = ols(&points);
    let mut by_size: BTreeMap<u64, Vec<(f64, f64)>> = BTreeMap::new();
    for t in &targets {
        if let Some(y) = log2_ratio(t) {
            by_size.entry(t.m).or_default().push((t.m as f64, y));
        }
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
        (!slopes.is_empty()).then(|| slopes[((q * (slopes.len() - 1) as f64).round()) as usize])
    };
    let (ci_lo, ci_hi) = (pct(0.025), pct(0.975));

    // Per-relation growth of F6's charged geometry against F4's word XORs.
    let per_rel = |f: &dyn Fn(&Target) -> Option<f64>| -> Option<f64> {
        let pts: Vec<(f64, f64)> = targets
            .iter()
            .filter_map(|t| Some((t.m as f64, f(t)?.log2())))
            .collect();
        ols(&pts).map(|(b, _)| b)
    };
    let geometry_slope = per_rel(&|t| Some(t.f6.geometry? as f64 / t.f6.relations?.max(1) as f64));
    let w4_slope = per_rel(&|t| Some(t.f4.word_xors? as f64 / t.f4.relations?.max(1) as f64));

    // Exponents of the lower-bound S against m (log2 S per unit m).
    let s_slope = |arm: &str| -> Option<f64> {
        let pts: Vec<(f64, f64)> = targets
            .iter()
            .filter_map(|t| {
                let s = if arm == "ic-f4" { t.f4.s? } else { t.f6.s? };
                Some((t.m as f64, s.log2()))
            })
            .collect();
        ols(&pts).map(|(b, _)| b)
    };

    let undecidable = sizes
        .iter()
        .any(|m| completed_by_size.get(m).copied().unwrap_or(0) < MIN_TARGETS)
        || sizes.len() < 4;
    let supported = !undecidable
        && ci_lo.is_some_and(|lo| lo > 0.0)
        && matches!((geometry_slope, w4_slope), (Some(g), Some(w)) if g <= w);
    let decision = if undecidable {
        "not decidable: fewer than six targets completed at some size, or fewer than four sizes"
    } else if supported {
        "supported: a stage-level exponent lead (reductions unpriced; no whole-method claim)"
    } else {
        "rejected: the F6-IC gain is a constant factor; class engineering; close the F6 thread"
    };

    let rows: Vec<Value> = targets
        .iter()
        .map(|t| {
            json!({
                "m": t.m,
                "workload": t.workload,
                "w4_word_xors": t.f4.word_xors,
                "w6_word_xors": t.f6.word_xors,
                "log2_w4_over_w6": log2_ratio(t),
                "f4_reductions": t.f4.reductions,
                "f6_reductions": t.f6.reductions,
                "g6_geometric_additions": t.f6.geometry,
                "relations": t.f4.relations,
                "trials": t.f4.trials,
                "s_lower_f4": t.f4.s,
                "s_lower_f6": t.f6.s,
                "s_rho": t.rho.as_ref().filter(|r| r.verified).and_then(|r| r.s),
            })
        })
        .collect();
    let out = json!({
        "schema": "f6-ic-ladder-analysis/v1",
        "targets": rows,
        "protocol": "research/f6_ic_ecbench_ladder_20261005/PROTOCOL.md",
        "sessions": sessions,
        "round_disagreements": round_disagreements,
        "incomplete_targets": incomplete,
        "divergent_query_streams": divergent,
        "per_size": per_size,
        "fit": {
            "statistic": "log2(W4/W6), W = Boolean solver word XORs over all PDP calls",
            "points": points.len(),
            "slope_per_unit_m": fit.map(|f| f.0),
            "intercept": fit.map(|f| f.1),
            "bootstrap": {"samples": slopes.len(), "seed": BOOTSTRAP_SEED,
                          "ci95": [ci_lo, ci_hi], "scheme": "targets resampled within each size"},
            "slope_log2_geometry_per_relation_f6": geometry_slope,
            "slope_log2_word_xors_per_relation_f4": w4_slope,
            "slope_log2_s_lower_ic_f4": s_slope("ic-f4"),
            "slope_log2_s_lower_ic_f6": s_slope("ic-f6"),
        },
        "decision": decision,
    });
    let text = serde_json::to_string_pretty(&out).map_err(|e| e.to_string())?;
    std::fs::write(&args[2], format!("{text}\n")).map_err(|e| format!("{}: {e}", args[2]))?;
    println!("{text}");
    Ok(())
}
