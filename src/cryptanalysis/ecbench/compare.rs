//! Comparing two arms.
//!
//! **Operations** are the metric (AGENTS.md §2, §6).  Runs are matched by
//! `(workload, round)`: within one session the two arms of a round ran
//! the same workload under the same algorithm seed, and across two
//! sessions of one spec they did too.  The ratio is
//!
//! ```text
//!   B / A  =  Σ S_B / Σ S_A        over the matched, verified pairs
//! ```
//!
//! which is AGENTS.md §8's `baseline_total / candidate_total` in `S`: a
//! ratio of totals, never a mean of ratios, which is biased whenever the
//! arms spread differently (rho's cost varies with its seed, BSGS's does
//! not).  Its interval is a two-stage bootstrap of the same statistic:
//! workloads resampled, then pairs within each.  With one workload only
//! the second stage exists and the comparison says so.
//!
//! **Wall time** is a practicality note.  It is compared only inside one
//! session (same host class, interleaved), only over pairs where both
//! runs earned the spec's `isolation_required`, only with at least five
//! such pairs, and it is set beside the A/A interval of the session's two
//! arms that run one method, when it has five admitted pairs too.
//! Anything less is printed as descriptive and is never a result.
//!
//! A comparison is refused, not qualified, when the arms ran different
//! workloads or units.  It is `incomplete` when a measured run of either
//! arm did not verify, `partial` when some runs have no partner, and
//! `bounded` when either arm left work unpriced.  Classifying a change as
//! advance, engineering, relabelling or accounting (AGENTS.md §3) stays
//! with the author.

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::json;

use crate::cryptanalysis::ecbench::canonical::short_id;
use crate::cryptanalysis::ecbench::record::Record;
use crate::cryptanalysis::ecbench::runner::{read_records, read_session};
use crate::cryptanalysis::ecbench::spec::{Level, Spec};
use crate::cryptanalysis::ecbench::stats::{bootstrap_ci, cluster_bootstrap_ci, mean, median};

pub const COMPARISON_SCHEMA: &str = "ecbench.comparison/v1";

/// Wall-time and A/A figures need at least this many admitted pairs.
pub const MIN_WALL_PAIRS: usize = 5;

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
    /// Matched, verified pairs on this workload.
    pub pairs: u64,
    pub a_mean_s: Option<f64>,
    pub b_mean_s: Option<f64>,
    /// `Σ S_B / Σ S_A` over the pairs.
    pub ratio_b_over_a: Option<f64>,
}

/// One curve's ratio: the per-size figure a claim about scaling reads.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CurveRatio {
    pub slug: String,
    pub log2_r: f64,
    pub workloads: u64,
    pub pairs: u64,
    pub ratio_b_over_a: Option<f64>,
    pub ci95: Option<(f64, f64)>,
    pub ci_method: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct OpsResult {
    /// `ok`; `incomplete` (a measured run did not verify); `partial`
    /// (some verified runs have no partner); or `empty`.
    pub status: String,
    pub unit: String,
    /// `Σ S_B / Σ S_A` over every matched pair, all curves pooled.
    pub ratio_b_over_a: Option<f64>,
    pub ci95: Option<(f64, f64)>,
    /// `cluster` (workloads, then pairs, resampled; two or more
    /// workloads) or `within` (pairs only: one workload, so the interval
    /// says nothing about other targets).
    pub ci_method: String,
    pub pairs: u64,
    /// Verified measured runs of either arm with no partner.
    pub unpaired_runs: u64,
    /// The pairs ran under the same algorithm seeds.
    pub same_seeds: bool,
    pub same_session: bool,
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
    /// same method, when it has at least five admitted pairs: the noise
    /// floor.
    pub aa_arms: Option<(String, String)>,
    pub aa_pairs: u64,
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

/// `Σ b / Σ a` over pairs, `None` when undefined.
fn ratio_of_sums(pairs: impl Iterator<Item = (f64, f64)>) -> Option<f64> {
    let (mut sa, mut sb, mut n) = (0.0, 0.0, 0usize);
    for (a, b) in pairs {
        sa += a;
        sb += b;
        n += 1;
    }
    (n > 0 && sa > 0.0).then(|| sb / sa)
}

/// Matched pairs `(S_A, S_B)` per workload, with the seeds' agreement.
struct Matched {
    by_workload: BTreeMap<String, Vec<(f64, f64)>>,
    unpaired: u64,
    same_seeds: bool,
}

fn match_pairs(ma: &[&Record], mb: &[&Record]) -> Matched {
    let key = |r: &Record| (r.workload.workload_id.clone(), r.round);
    let a: BTreeMap<(String, u32), &Record> = ma
        .iter()
        .filter(|r| r.counts() && r.cost.s.is_some())
        .map(|r| (key(r), *r))
        .collect();
    let b: BTreeMap<(String, u32), &Record> = mb
        .iter()
        .filter(|r| r.counts() && r.cost.s.is_some())
        .map(|r| (key(r), *r))
        .collect();
    let mut by_workload: BTreeMap<String, Vec<(f64, f64)>> = BTreeMap::new();
    let mut same_seeds = true;
    for (k, ra) in &a {
        if let Some(rb) = b.get(k) {
            same_seeds &= ra.algorithm_seed == rb.algorithm_seed;
            by_workload
                .entry(k.0.clone())
                .or_default()
                .push((ra.cost.s.unwrap_or(0.0), rb.cost.s.unwrap_or(0.0)));
        }
    }
    let paired: usize = by_workload.values().map(|v| v.len()).sum();
    let unpaired = (a.len() - paired) as u64 + (b.len() - paired) as u64;
    Matched {
        by_workload,
        unpaired,
        same_seeds,
    }
}

/// The ratio of sums over `strata` with its interval: two-stage when
/// there are two or more strata, within-stratum otherwise.
fn ratio_with_interval(
    strata: &[Vec<(f64, f64)>],
    resamples: usize,
    seed: u64,
) -> (Option<f64>, Option<(f64, f64)>, String) {
    let stat = |s: &[Vec<(f64, f64)>]| ratio_of_sums(s.iter().flatten().copied());
    let ratio = stat(strata);
    if strata.len() >= 2 {
        (
            ratio,
            cluster_bootstrap_ci(strata, resamples, seed, stat),
            "cluster".into(),
        )
    } else {
        (
            ratio,
            bootstrap_ci(strata, resamples, seed, stat),
            "within".into(),
        )
    }
}

/// Compare arm `b` against arm `a`.  `b_dir` names a second session (a
/// candidate binary's run of the same spec); wall time is then refused.
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

    let matched = match_pairs(&ma, &mb);
    let workload_of: BTreeMap<&str, &Record> = ma
        .iter()
        .map(|r| (r.workload.workload_id.as_str(), *r))
        .collect();
    let empty = Vec::new();
    let rows: Vec<WorkloadRow> = wa
        .iter()
        .map(|wid| {
            let pairs = matched.by_workload.get(*wid).unwrap_or(&empty);
            let any = workload_of[wid];
            WorkloadRow {
                workload_id: wid.to_string(),
                slug: any.workload.curve.slug.clone(),
                r: any.workload.curve.r,
                floor_s: any.boundaries.floor_s,
                pairs: pairs.len() as u64,
                a_mean_s: mean(&pairs.iter().map(|p| p.0).collect::<Vec<_>>()),
                b_mean_s: mean(&pairs.iter().map(|p| p.1).collect::<Vec<_>>()),
                ratio_b_over_a: ratio_of_sums(pairs.iter().copied()),
            }
        })
        .collect();
    let strata: Vec<Vec<(f64, f64)>> = matched
        .by_workload
        .values()
        .filter(|v| !v.is_empty())
        .cloned()
        .collect();
    let (ratio, ci, ci_method) = ratio_with_interval(&strata, resamples, seed);
    let mut slugs: Vec<String> = rows.iter().map(|r| r.slug.clone()).collect();
    slugs.sort();
    slugs.dedup();
    let curves: Vec<CurveRatio> = slugs
        .iter()
        .map(|slug| {
            let mine: Vec<&WorkloadRow> = rows.iter().filter(|r| &r.slug == slug).collect();
            let strata: Vec<Vec<(f64, f64)>> = mine
                .iter()
                .filter_map(|r| matched.by_workload.get(&r.workload_id))
                .filter(|v| !v.is_empty())
                .cloned()
                .collect();
            let (ratio_b_over_a, ci95, ci_method) = ratio_with_interval(&strata, resamples, seed);
            CurveRatio {
                slug: slug.clone(),
                log2_r: (mine[0].r as f64).log2(),
                workloads: mine.len() as u64,
                pairs: strata.iter().map(|s| s.len() as u64).sum(),
                ratio_b_over_a,
                ci95,
                ci_method,
            }
        })
        .collect();
    let pairs: u64 = strata.iter().map(|s| s.len() as u64).sum();
    let incomplete = a.verified < a.measured || b.verified < b.measured;
    let ops = OpsResult {
        status: if ratio.is_none() {
            "empty".into()
        } else if incomplete {
            "incomplete".into()
        } else if matched.unpaired > 0 {
            "partial".into()
        } else {
            "ok".into()
        },
        unit: ma[0].cost.unit.clone(),
        ratio_b_over_a: ratio,
        ci95: ci,
        ci_method,
        pairs,
        unpaired_runs: matched.unpaired,
        same_seeds: matched.same_seeds,
        same_session,
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
        "ops ({}): {} / {} = {}{} ({}) over {} pairs on {} workloads, {} curve(s){}{}; wall: {}{}",
        ops.status,
        b_arm,
        a_arm,
        ops.ratio_b_over_a
            .map(|r| format!("{r:.4}"))
            .unwrap_or_else(|| "unknown".into()),
        fmt_ci(ops.ci95),
        ops.ci_method,
        ops.pairs,
        rows.len(),
        ops.pooled_curves,
        if ops.bounded {
            ", bounded (unpriced work)"
        } else {
            ""
        },
        match ops.status.as_str() {
            "incomplete" => ", NOT admissible: a measured run did not verify",
            "partial" => ", partial: some runs have no partner",
            _ => "",
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
/// verified and both earned `required`; `(ratios by workload, pairs,
/// pairs below level)`.
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
        aa_pairs: 0,
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
    if admitted < MIN_WALL_PAIRS {
        w.reasons.push(format!(
            "{admitted} admitted pairs; at least {MIN_WALL_PAIRS} are needed"
        ));
        w.status = "descriptive".into();
        return w;
    }
    w.median_ratio = med(&strata);
    w.ci95 = bootstrap_ci(&strata, resamples, seed ^ 0xA11, med);
    w.status = "admitted".into();
    // The noise floor: two arms of the session that run one method, with
    // enough admitted pairs to have an interval at all.
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
                let n: usize = s.iter().map(|v| v.len()).sum();
                w.aa_arms = Some((names[i].0.clone(), names[j].0.clone()));
                w.aa_pairs = n as u64;
                if n >= MIN_WALL_PAIRS {
                    w.aa_median_ratio = med(&s);
                    w.aa_ci95 = bootstrap_ci(&s, resamples, seed ^ 0xAA, med);
                } else {
                    w.reasons.push(format!(
                        "the A/A arms have {n} admitted pairs; at least {MIN_WALL_PAIRS} are needed for a noise floor"
                    ));
                }
                break 'outer;
            }
        }
    }
    if let (Some((lo, hi)), Some((alo, ahi))) = (w.ci95, w.aa_ci95) {
        w.outside_noise = Some((hi < 1.0 || lo > 1.0) && (hi < alo || lo > ahi));
    }
    w
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ratio_of_sums_is_the_totals_ratio() {
        // Two workloads, unequal spread: the mean of per-workload ratios
        // (1.5) is not the totals ratio (11/9).
        let pairs = [(1.0, 2.0), (8.0, 9.0)];
        assert_eq!(ratio_of_sums(pairs.iter().copied()), Some(11.0 / 9.0));
        assert_eq!(ratio_of_sums(std::iter::empty()), None);
    }

    #[test]
    fn identical_arms_give_a_unit_ratio_and_a_point_interval() {
        let strata = vec![vec![(1.0, 1.0), (2.0, 2.0)], vec![(3.0, 3.0)]];
        let (r, ci, method) = ratio_with_interval(&strata, 500, 1);
        assert_eq!(r, Some(1.0));
        assert_eq!(ci, Some((1.0, 1.0)));
        assert_eq!(method, "cluster");
    }
}
