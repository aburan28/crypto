//! The pricing band of AGENTS.md §5a.
//!
//! Almost no ECDLP cost in a record is complete: rho, BSGS and the kangaroo
//! leave canonicalisations, partition hashes and table work unpriced, and an
//! index-calculus arm is complete only where `docs/ic/calibration.json` has a
//! pinned row for its curve.  Comparing a calibrated IC cost with an
//! uncalibrated rho cost is lopsided.  The band prices every arm alike: the
//! native work a record leaves unpriced *and* the native work a pinned
//! calibration priced are both charged at one weight `w` (group-addition
//! equivalents per native operation), and two arms are ranked only when the
//! order holds at every weight of the band.
//!
//! For a record, with `G` its charged GAE, `C` the GAE its calibration put on
//! native work, `n_C` the operations behind `C`, and `U` the count behind its
//! unpriced counters,
//!
//! ```text
//! priced_gae(w) = G − C + w · (U + n_C),   priced_S(w) = priced_gae(w) / √r.
//! ```
//!
//! A record that names an unpriced counter without a count has no priced
//! cost above `w = 0`; it stays a lower bound and cannot be ranked.

use std::collections::BTreeMap;

use serde::Serialize;

use super::methods::PRICED_NATIVE;
use super::record::Record;
use crate::cryptanalysis::ic_boundary::pinned_calibration_ratios;

/// The band every comparison uses unless it states a wider one.
pub const DEFAULT_WEIGHTS: [f64; 3] = [0.0, 0.1, 1.0];

/// The count behind each unpriced counter of `r`, `None` where the record
/// names the counter but carries no count for it.  A counter `X_uncharged`
/// is looked up in the top-level counters, then in each phase's native map
/// as `X_uncharged` or `X`, and a `solver_*` counter in the solver block.
pub fn unpriced_counts(r: &Record) -> BTreeMap<String, Option<u64>> {
    let mut out = BTreeMap::new();
    for name in &r.cost.unpriced {
        let base = name.strip_suffix("_uncharged").unwrap_or(name);
        if let Some(&v) = r.counters.get(name) {
            out.insert(name.clone(), Some(v));
            continue;
        }
        let mut total = 0u64;
        let mut found = false;
        for p in &r.phases {
            if let Some(&v) = p.native.get(name).or_else(|| p.native.get(base)) {
                total += v;
                found = true;
            }
        }
        if !found && base.starts_with("solver_") {
            if let Some(s) = &r.solver {
                out.insert(name.clone(), Some(s.ops));
                continue;
            }
        }
        out.insert(name.clone(), found.then_some(total));
    }
    out
}

/// Total unpriced operations of `r` (counted ones only) and whether any
/// unpriced counter has no count.
pub fn unpriced_total(r: &Record) -> (u64, bool) {
    let counts = unpriced_counts(r);
    let missing = counts.values().any(Option::is_none);
    (counts.values().flatten().sum(), missing)
}

fn regime_of(r: &Record) -> &'static str {
    match r.workload.curve.family.as_str() {
        "koblitz" => "koblitz",
        "prime" => "prime",
        _ => "char2",
    }
}

/// The GAE a pinned calibration put on native work in `r`, and the number
/// of operations behind it.  `Err` when the record says it was priced by a
/// pinned calibration that this binary's table cannot reproduce.
pub fn calibrated_native(r: &Record) -> Result<(f64, u64), String> {
    let pinned: Vec<String> = r
        .detail
        .get("calibration_pinned_units")
        .and_then(|v| v.as_array())
        .map(|a| {
            a.iter()
                .filter_map(|x| x.as_str().map(str::to_string))
                .collect()
        })
        .unwrap_or_default();
    if pinned.is_empty() {
        return Ok((0.0, 0));
    }
    let table =
        pinned_calibration_ratios(regime_of(r), &r.workload.curve.slug).ok_or_else(|| {
            format!(
                "{}: priced by pinned units {:?} but no calibration row for {}",
                r.record_id, pinned, r.workload.curve.slug
            )
        })?;
    let ratio = |unit: &str| table.get(unit).and_then(|v| v.as_f64());
    let mut gae = 0.0;
    let mut ops = 0u64;
    for p in &r.phases {
        for (counter, unit) in PRICED_NATIVE {
            let n = p.native.get(*counter).copied().unwrap_or(0);
            if n > 0 && pinned.iter().any(|u| u == unit) {
                let k = ratio(unit).ok_or_else(|| format!("{}: no ratio {unit}", r.record_id))?;
                gae += n as f64 * k;
                ops += n;
            }
        }
    }
    if let Some(s) = &r.solver {
        if s.ops > 0
            && s.op_unit.contains("word XOR")
            && pinned.iter().any(|u| u == "ns_per_word_xor")
            && !r.cost.unpriced.iter().any(|u| u.starts_with("solver_"))
        {
            let k = ratio("ns_per_word_xor")
                .ok_or_else(|| format!("{}: no ratio ns_per_word_xor", r.record_id))?;
            gae += s.ops as f64 * k;
            ops += s.ops;
        }
    }
    Ok((gae, ops))
}

/// `priced_gae(w)` of the module comment; `None` for a record with no
/// recorded cost, or one with an uncounted unpriced counter at `w > 0`.
pub fn priced_gae(r: &Record, w: f64) -> Result<Option<f64>, String> {
    let Some(g) = r.cost.total_gae else {
        return Ok(None);
    };
    let (u, missing) = unpriced_total(r);
    if missing && w > 0.0 {
        return Ok(None);
    }
    let (c, n_c) = calibrated_native(r)?;
    Ok(Some(g - c + w * (u + n_c) as f64))
}

/// `priced_S(w)`.
pub fn priced_s(r: &Record, w: f64) -> Result<Option<f64>, String> {
    let sq = (r.workload.curve.r as f64).sqrt();
    Ok(priced_gae(r, w)?.map(|g| g / sq))
}

/// One arm on one curve across the band.
#[derive(Clone, Debug, Serialize)]
pub struct BandRow {
    pub curve: String,
    pub log2_r: f64,
    pub arm: String,
    pub runs: u64,
    pub verified: u64,
    /// Mean priced `S` over verified runs at each weight, `None` when any
    /// verified run has no priced cost there.
    pub mean_s: Vec<Option<f64>>,
    /// Unpriced counters that carry no count on some run.
    pub uncounted: Vec<String>,
}

/// The verdict for a pair of arms on one curve.
#[derive(Clone, Debug, Serialize)]
pub struct BandRank {
    pub curve: String,
    pub a: String,
    pub b: String,
    pub verdict: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct BandReport {
    pub schema: String,
    pub weights: Vec<f64>,
    pub rows: Vec<BandRow>,
    pub ranks: Vec<BandRank>,
}

/// The band report over `records` (warm-ups ignored).
pub fn band(records: &[Record], weights: &[f64]) -> Result<BandReport, String> {
    let mut by: BTreeMap<(String, String), Vec<&Record>> = BTreeMap::new();
    for r in records.iter().filter(|r| !r.warmup) {
        by.entry((r.workload.curve.slug.clone(), r.arm.clone()))
            .or_default()
            .push(r);
    }
    let mut rows = Vec::new();
    for ((curve, arm), rs) in &by {
        let ok: Vec<&&Record> = rs.iter().filter(|r| r.counts()).collect();
        let mut mean_s = Vec::new();
        for &w in weights {
            let mut vals = Vec::new();
            let mut all = !ok.is_empty();
            for r in &ok {
                match priced_s(r, w)? {
                    Some(v) => vals.push(v),
                    None => all = false,
                }
            }
            mean_s.push(all.then(|| vals.iter().sum::<f64>() / vals.len() as f64));
        }
        let mut uncounted: Vec<String> = rs
            .iter()
            .flat_map(|r| {
                unpriced_counts(r)
                    .into_iter()
                    .filter(|(_, v)| v.is_none())
                    .map(|(k, _)| k)
            })
            .collect();
        uncounted.sort();
        uncounted.dedup();
        rows.push(BandRow {
            curve: curve.clone(),
            log2_r: (rs[0].workload.curve.r as f64).log2(),
            arm: arm.clone(),
            runs: rs.len() as u64,
            verified: ok.len() as u64,
            mean_s,
            uncounted,
        });
    }
    rows.sort_by(|a, b| {
        a.log2_r
            .partial_cmp(&b.log2_r)
            .unwrap()
            .then(a.arm.cmp(&b.arm))
    });
    let mut ranks = Vec::new();
    for (i, a) in rows.iter().enumerate() {
        for b in rows[i + 1..].iter().filter(|b| b.curve == a.curve) {
            let verdict = if a.verified < a.runs || b.verified < b.runs {
                "incomplete".to_string()
            } else if a.mean_s.iter().chain(&b.mean_s).any(Option::is_none) {
                "unrankable (uncounted unpriced work)".to_string()
            } else {
                let lt = a
                    .mean_s
                    .iter()
                    .zip(&b.mean_s)
                    .all(|(x, y)| x.unwrap() < y.unwrap());
                let gt = a
                    .mean_s
                    .iter()
                    .zip(&b.mean_s)
                    .all(|(x, y)| x.unwrap() > y.unwrap());
                match (lt, gt) {
                    (true, _) => format!("{} < {}", a.arm, b.arm),
                    (_, true) => format!("{} < {}", b.arm, a.arm),
                    _ => "unresolved inside the band".to_string(),
                }
            };
            ranks.push(BandRank {
                curve: a.curve.clone(),
                a: a.arm.clone(),
                b: b.arm.clone(),
                verdict,
            });
        }
    }
    Ok(BandReport {
        schema: "ecbench.band/v1".into(),
        weights: weights.to_vec(),
        rows,
        ranks,
    })
}

/// The band as a Markdown table.
pub fn markdown(rep: &BandReport) -> String {
    let mut s = String::from("| curve | log2 r | arm | verified/runs |");
    for w in &rep.weights {
        s.push_str(&format!(" S(w={w}) |"));
    }
    s.push_str(" uncounted |\n|:--|--:|:--|--:|");
    for _ in &rep.weights {
        s.push_str("--:|");
    }
    s.push_str(":--|\n");
    for r in &rep.rows {
        s.push_str(&format!(
            "| {} | {:.2} | {} | {}/{} |",
            r.curve, r.log2_r, r.arm, r.verified, r.runs
        ));
        for v in &r.mean_s {
            match v {
                Some(v) => s.push_str(&format!(" {v:.3} |")),
                None => s.push_str(" — |"),
            }
        }
        let u = if r.uncounted.is_empty() {
            "—".to_string()
        } else {
            r.uncounted.join(", ")
        };
        s.push_str(&format!(" {u} |\n"));
    }
    s
}
