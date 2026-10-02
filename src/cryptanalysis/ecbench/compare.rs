//! Comparing two arms.
//!
//! **Operations** are the metric (AGENTS.md §2, §6).  For each workload
//! the mean `S` of each arm's verified measured runs is taken, and the
//! ratio `B / A` is the geometric mean of the per-workload ratios, with a
//! stratified bootstrap interval.  Within one session the runs are paired
//! by round (both arms ran under the same algorithm seed) and resampled
//! as pairs.  The mean, not the median, is the statistic: the expected
//! cost `√(πr/2A)` the floor is stated in is a mean.
//!
//! **Wall time** is a practicality note.  It is compared only inside one
//! session (same host class, interleaved), only over pairs where both
//! runs earned the spec's `isolation_required`, only with at least five
//! such rounds, and it is set beside the A/A interval of the session's
//! control arm when there is one.  Anything less is printed as
//! descriptive and is never a result.
//!
//! A comparison is refused, not qualified, when the arms ran different
//! workloads or units.  It is reported `incomplete` when any measured run
//! of either arm did not verify, and `bounded` when either arm left work
//! unpriced.  Classifying a change as advance, engineering, relabelling
//! or accounting (AGENTS.md §3) stays with the author.

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::json;

use crate::cryptanalysis::ecbench::canonical::short_id;
use crate::cryptanalysis::ecbench::record::Record;
use crate::cryptanalysis::ecbench::runner::{read_records, read_session};
use crate::cryptanalysis::ecbench::spec::{Level, Spec};
use crate::cryptanalysis::ecbench::stats::{
    bootstrap_ci, cluster_bootstrap_ci, geomean, mean, median,
};

pub const COMPARISON_SCHEMA: &str = "ecbench.comparison/v1";

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ArmSummary {
    pub session_id: String,
    pub arm: String,
    pub method_id: String,
    pub method: String,
    pub measured: u64,
    pub verified: u64,
    pub statuses: BTreeMap<String, u64>,
    pub mean_s: Option<f64>,
    pub mean_ratio_to_floor: Option<f64>,
    pub lower_bound: bool,
    pub deterministic: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WorkloadRow {
    pub workload_id: String,
    pub slug: String,
    pub r: u64,
    pub floor_s: f64,
    pub a_runs: u64,
    pub b_runs: u64,
    pub a_mean_s: Option<f64>,
    pub b_mean_s: Option<f64>,
    pub ratio_b_over_a: Option<f64>,
}

/// One curve's ratio: the per-size figure a claim about scaling reads.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CurveRatio {
    pub slug: String,
    pub log2_r: f64,
    pub workloads: u64,
    pub ratio_b_over_a: Option<f64>,
    pub ci95: Option<(f64, f64)>,
    pub ci_method: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct OpsResult {
    /// `ok`, `incomplete` (a measured run did not verify) or `empty`.
    pub status: String,
    pub unit: String,
    /// Geometric mean over workloads of mean-S ratios.
    pub ratio_b_over_a: Option<f64>,
    pub ci95: Option<(f64, f64)>,
    /// `cluster` (workloads and runs resampled, the default with two or
    /// more workloads) or `within` (runs only: one workload, so the
    /// interval says nothing about other targets).
    pub ci_method: String,
    pub paired: bool,
    /// Pooled over curves of different sizes; read `curves` for scaling.
    pub pooled_curves: u64,
    /// Either arm left work unpriced: the ratio is between bounds.
    pub bounded: bool,
    pub resamples: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WallResult {
    /// `admitted`, `descriptive` or `refused`.
    pub status: String,
    pub reasons: Vec<String>,
    pub required_level: Level,
    pub pairs: u64,
    pub pairs_below_level: u64,
    /// Median over pairs of `solve_wall(B) / solve_wall(A)`.
    pub median_ratio: Option<f64>,
    pub ci95: Option<(f64, f64)>,
    /// The same statistic between the session's two arms that run the
    /// same method, when it has them: the noise floor.
    pub aa_arms: Option<(String, String)>,
    pub aa_median_ratio: Option<f64>,
    pub aa_ci95: Option<(f64, f64)>,
    /// The B/A interval excludes 1 and does not overlap the A/A interval.
    pub outside_noise: Option<bool>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Comparison {
    pub schema: String,
    pub comparison_id: String,
    pub a: ArmSummary,
    pub b: ArmSummary,
    pub workloads: Vec<WorkloadRow>,
    pub curves: Vec<CurveRatio>,
    pub ops: OpsResult,
    pub wall: WallResult,
    pub verdict: String,
}

fn summarize(session_id: &str, arm: &str, recs: &[&Record]) -> Result<ArmSummary, String> {
    let first = recs
        .first()
        .ok_or_else(|| format!("arm `{arm}` has no measured runs"))?;
    let mut statuses = BTreeMap::new();
    for r in recs {
        *statuses.entry(r.outcome.status.clone()).or_insert(0) += 1;
    }
    let ok: Vec<&&Record> = recs.iter().filter(|r| r.counts()).collect();
    let s: Vec<f64> = ok.iter().filter_map(|r| r.cost.s).collect();
    let f: Vec<f64> = ok
        .iter()
        .filter_map(|r| r.boundaries.ratio_to_floor)
        .collect();
    Ok(ArmSummary {
        session_id: session_id.into(),
        arm: arm.into(),
        method_id: first.method.method_id.clone(),
        method: first.method.id.clone(),
        measured: recs.len() as u64,
        verified: ok.len() as u64,
        statuses,
        mean_s: mean(&s),
        mean_ratio_to_floor: mean(&f),
        lower_bound: recs.iter().any(|r| r.cost.lower_bound),
        deterministic: recs.iter().all(|r| r.cost.deterministic),
    })
}

fn measured<'a>(recs: &'a [Record], arm: &str) -> Vec<&'a Record> {
    recs.iter().filter(|r| r.arm == arm && !r.warmup).collect()
}

/// Compare arm `b` against arm `a`.  `b_dir` names a second session for
/// a cross-session operations comparison; wall time is then refused.
pub fn compare(
    a_dir: &Path,
    a_arm: &str,
    b_dir: Option<&Path>,
    b_arm: &str,
    resamples: usize,
    seed: u64,
) -> Result<Comparison, String> {
    let sa = read_session(a_dir)?;
    let ra = read_records(a_dir)?;
    let (sb, rb) = match b_dir {
        Some(d) => (read_session(d)?, read_records(d)?),
        None => (sa.clone(), ra.clone()),
    };
    let same_session = sa.session_id == sb.session_id;
    if same_session && a_arm == b_arm {
        return Err("compare two different arms".into());
    }
    let ma = measured(&ra, a_arm);
    let mb = measured(&rb, b_arm);
    let a = summarize(&sa.session_id, a_arm, &ma)?;
    let b = summarize(&sb.session_id, b_arm, &mb)?;

    // Refusals: different workloads, different units.
    let wa: BTreeSet<&str> = ma.iter().map(|r| r.workload.workload_id.as_str()).collect();
    let wb: BTreeSet<&str> = mb.iter().map(|r| r.workload.workload_id.as_str()).collect();
    if wa != wb {
        return Err(format!(
            "the arms ran different workloads ({} vs {}); a comparison pairs identical workloads",
            wa.len(),
            wb.len()
        ));
    }
    if ma
        .iter()
        .chain(mb.iter())
        .any(|r| r.cost.unit != ma[0].cost.unit)
    {
        return Err("the arms report different units".into());
    }

    // Operations, per workload.
    let mut rows = Vec::new();
    let mut strata_pairs: Vec<Vec<(f64, f64)>> = Vec::new();
    let mut strata_a: Vec<Vec<f64>> = Vec::new();
    let mut strata_b: Vec<Vec<f64>> = Vec::new();
    for wid in &wa {
        let sa_runs: Vec<&&Record> = ma
            .iter()
            .filter(|r| r.workload.workload_id == *wid && r.counts())
            .collect();
        let sb_runs: Vec<&&Record> = mb
            .iter()
            .filter(|r| r.workload.workload_id == *wid && r.counts())
            .collect();
        let xa: Vec<f64> = sa_runs.iter().filter_map(|r| r.cost.s).collect();
        let xb: Vec<f64> = sb_runs.iter().filter_map(|r| r.cost.s).collect();
        let any = ma
            .iter()
            .find(|r| r.workload.workload_id == *wid)
            .expect("workload is in a");
        let (am, bm) = (mean(&xa), mean(&xb));
        rows.push(WorkloadRow {
            workload_id: wid.to_string(),
            slug: any.workload.curve.slug.clone(),
            r: any.workload.curve.r,
            floor_s: any.boundaries.floor_s,
            a_runs: xa.len() as u64,
            b_runs: xb.len() as u64,
            a_mean_s: am,
            b_mean_s: bm,
            ratio_b_over_a: match (am, bm) {
                (Some(a), Some(b)) if a > 0.0 => Some(b / a),
                _ => None,
            },
        });
        if same_session {
            let by_round_a: BTreeMap<u32, f64> = sa_runs
                .iter()
                .filter_map(|r| Some((r.round, r.cost.s?)))
                .collect();
            let pairs: Vec<(f64, f64)> = sb_runs
                .iter()
                .filter_map(|r| Some((*by_round_a.get(&r.round)?, r.cost.s?)))
                .collect();
            strata_pairs.push(pairs);
        }
        strata_a.push(xa);
        strata_b.push(xb);
    }
    let ratios: Vec<f64> = rows.iter().filter_map(|r| r.ratio_b_over_a).collect();
    let ratio = if ratios.len() == rows.len() {
        geomean(&ratios)
    } else {
        None
    };
    // One stratum per workload: (a, b) pairs by round in one session,
    // tagged runs of either arm across sessions.
    let strata: Vec<Vec<(u8, f64, f64)>> = if same_session {
        strata_pairs
            .iter()
            .map(|p| p.iter().map(|(x, y)| (2u8, *x, *y)).collect())
            .collect()
    } else {
        strata_a
            .iter()
            .zip(&strata_b)
            .map(|(x, y)| {
                x.iter()
                    .map(|v| (0u8, *v, 0.0))
                    .chain(y.iter().map(|v| (1u8, 0.0, *v)))
                    .collect()
            })
            .collect()
    };
    let stat = |s: &[Vec<(u8, f64, f64)>]| -> Option<f64> {
        let per: Vec<f64> = s
            .iter()
            .map(|p| {
                let a = mean(
                    &p.iter()
                        .filter(|x| x.0 != 1)
                        .map(|x| x.1)
                        .collect::<Vec<_>>(),
                )?;
                let b = mean(
                    &p.iter()
                        .filter(|x| x.0 != 0)
                        .map(|x| x.2)
                        .collect::<Vec<_>>(),
                )?;
                (a > 0.0).then(|| b / a)
            })
            .collect::<Option<Vec<f64>>>()?;
        geomean(&per)
    };
    let interval = |idx: &[usize], salt: u64| -> (Option<(f64, f64)>, String) {
        let sub: Vec<Vec<(u8, f64, f64)>> = idx.iter().map(|&i| strata[i].clone()).collect();
        if sub.len() >= 2 {
            (
                cluster_bootstrap_ci(&sub, resamples, seed ^ salt, stat),
                "cluster".into(),
            )
        } else {
            (
                bootstrap_ci(&sub, resamples, seed ^ salt, stat),
                "within".into(),
            )
        }
    };
    let all: Vec<usize> = (0..strata.len()).collect();
    let (ci, ci_method) = interval(&all, 0);
    let mut slugs: Vec<String> = rows.iter().map(|r| r.slug.clone()).collect();
    slugs.sort();
    slugs.dedup();
    let curves: Vec<CurveRatio> = slugs
        .iter()
        .map(|slug| {
            let idx: Vec<usize> = rows
                .iter()
                .enumerate()
                .filter(|(_, r)| &r.slug == slug)
                .map(|(i, _)| i)
                .collect();
            let rs: Vec<f64> = idx.iter().filter_map(|&i| rows[i].ratio_b_over_a).collect();
            let (ci95, ci_method) = interval(&idx, 0);
            CurveRatio {
                slug: slug.clone(),
                log2_r: (rows[idx[0]].r as f64).log2(),
                workloads: idx.len() as u64,
                ratio_b_over_a: if rs.len() == idx.len() {
                    geomean(&rs)
                } else {
                    None
                },
                ci95,
                ci_method,
            }
        })
        .collect();
    let incomplete = a.verified < a.measured || b.verified < b.measured;
    let ops = OpsResult {
        status: if ratio.is_none() {
            "empty".into()
        } else if incomplete {
            "incomplete".into()
        } else {
            "ok".into()
        },
        unit: ma[0].cost.unit.clone(),
        ratio_b_over_a: ratio,
        ci95: ci,
        ci_method,
        paired: same_session,
        pooled_curves: curves.len() as u64,
        bounded: a.lower_bound || b.lower_bound,
        resamples,
    };

    // Wall time.
    let spec: Option<Spec> = std::fs::read_to_string(a_dir.join("spec.json"))
        .ok()
        .and_then(|t| Spec::from_json(&t).ok());
    let required = spec
        .as_ref()
        .map(|s| s.measurement.isolation_required)
        .unwrap_or(Level::L2);
    let wall = wall_compare(&ra, a_arm, b_arm, same_session, required, resamples, seed);

    let fmt_ci = |c: Option<(f64, f64)>| {
        c.map(|(l, h)| format!(" [{l:.3}, {h:.3}]"))
            .unwrap_or_default()
    };
    let verdict = format!(
        "ops ({}): {} / {} = {}{} ({}) over {} workloads on {} curve(s){}{}; wall: {}{}",
        ops.status,
        b_arm,
        a_arm,
        ops.ratio_b_over_a
            .map(|r| format!("{r:.4}"))
            .unwrap_or_else(|| "unknown".into()),
        fmt_ci(ops.ci95),
        ops.ci_method,
        rows.len(),
        ops.pooled_curves,
        if ops.bounded {
            ", bounded (unpriced work)"
        } else {
            ""
        },
        if incomplete {
            ", NOT admissible: a measured run did not verify"
        } else {
            ""
        },
        wall.status,
        match (wall.status.as_str(), wall.median_ratio) {
            ("admitted", Some(m)) => format!(" median {m:.4}{}", fmt_ci(wall.ci95)),
            ("descriptive", Some(m)) => format!(
                " (median {m:.4}, not a result: {})",
                wall.reasons.join("; ")
            ),
            _ => format!(" ({})", wall.reasons.join("; ")),
        },
    );
    let (comparison_id, _) = short_id(
        "ECC1",
        &json!({
            "a": [a.session_id, a.arm, a.method_id],
            "b": [b.session_id, b.arm, b.method_id],
            "resamples": resamples,
            "seed": seed.to_string(),
        }),
    )?;
    Ok(Comparison {
        schema: COMPARISON_SCHEMA.into(),
        comparison_id,
        a,
        b,
        workloads: rows,
        curves,
        ops,
        wall,
        verdict,
    })
}

/// Per-pair wall ratios `B / A` over (workload, round) where both
/// verified and both earned `required`; `(ratios by workload, pairs
/// below level)`.
fn wall_pairs(recs: &[Record], a: &str, b: &str, required: Level) -> (Vec<Vec<f64>>, u64, u64) {
    let mut by: BTreeMap<(String, u32), (Option<&Record>, Option<&Record>)> = BTreeMap::new();
    for r in recs.iter().filter(|r| r.counts()) {
        let key = (r.workload.workload_id.clone(), r.round);
        if r.arm == a {
            by.entry(key).or_default().0 = Some(r);
        } else if r.arm == b {
            by.entry(key).or_default().1 = Some(r);
        }
    }
    let mut strata: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let (mut pairs, mut below) = (0u64, 0u64);
    for ((wid, _), (ra, rb)) in by {
        let (Some(ra), Some(rb)) = (ra, rb) else {
            continue;
        };
        let (Some(ta), Some(tb)) = (ra.time.solve_wall_ns, rb.time.solve_wall_ns) else {
            continue;
        };
        pairs += 1;
        if ra.isolation.level < required || rb.isolation.level < required {
            below += 1;
            continue;
        }
        if ta > 0 {
            strata.entry(wid).or_default().push(tb as f64 / ta as f64);
        }
    }
    (strata.into_values().collect(), pairs, below)
}

fn wall_compare(
    recs: &[Record],
    a: &str,
    b: &str,
    same_session: bool,
    required: Level,
    resamples: usize,
    seed: u64,
) -> WallResult {
    let mut w = WallResult {
        status: "refused".into(),
        reasons: vec![],
        required_level: required,
        pairs: 0,
        pairs_below_level: 0,
        median_ratio: None,
        ci95: None,
        aa_arms: None,
        aa_median_ratio: None,
        aa_ci95: None,
        outside_noise: None,
    };
    if !same_session {
        w.reasons
            .push("the arms ran in different sessions, so no pair shares a host window".into());
        return w;
    }
    let med = |s: &[Vec<f64>]| median(&s.concat());
    let (strata, pairs, below) = wall_pairs(recs, a, b, required);
    w.pairs = pairs;
    w.pairs_below_level = below;
    let admitted: usize = strata.iter().map(|s| s.len()).sum();
    // Descriptive figure over every verified pair, whatever its level.
    let (all, _, _) = wall_pairs(recs, a, b, Level::L0);
    w.median_ratio = med(&all);
    if below > 0 {
        w.reasons.push(format!(
            "{below} of {pairs} pairs did not reach {}",
            required.name()
        ));
    }
    if admitted < 5 {
        w.reasons
            .push(format!("{admitted} admitted pairs; at least 5 are needed"));
        w.status = "descriptive".into();
        return w;
    }
    w.median_ratio = med(&strata);
    w.ci95 = bootstrap_ci(&strata, resamples, seed ^ 0xA11, med);
    w.status = "admitted".into();
    // The noise floor: two arms of the session that run one method.
    let mut arms: BTreeMap<String, String> = BTreeMap::new();
    for r in recs {
        arms.entry(r.arm.clone())
            .or_insert_with(|| r.method.method_id.clone());
    }
    let names: Vec<(&String, &String)> = arms.iter().collect();
    'outer: for i in 0..names.len() {
        for j in (i + 1)..names.len() {
            if names[i].1 == names[j].1 {
                let (s, _, _) = wall_pairs(recs, names[i].0, names[j].0, required);
                w.aa_arms = Some((names[i].0.clone(), names[j].0.clone()));
                w.aa_median_ratio = med(&s);
                w.aa_ci95 = bootstrap_ci(&s, resamples, seed ^ 0xAA, med);
                break 'outer;
            }
        }
    }
    if let (Some((lo, hi)), Some((alo, ahi))) = (w.ci95, w.aa_ci95) {
        w.outside_noise = Some((hi < 1.0 || lo > 1.0) && (hi < alo || lo > ahi));
    }
    w
}
