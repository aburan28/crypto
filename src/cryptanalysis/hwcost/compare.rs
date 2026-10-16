//! Two builds measured in alternation (A1, B1, A2, B2): per arm, the
//! candidate-to-baseline ratios `r_k = B_k / A_k`, the same-build drift
//! `d = |A2 / A1 − 1|`, and a classification that calls an effect only when
//! both ratios clear `max(2d, floor)` on the same side.  Also the
//! `fp64bf / fp64` kernel ratio inside each run (round-2 claim F1).

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};

use super::measure::Measurements;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ArmCompare {
    pub arm: String,
    pub a1: f64,
    pub a2: f64,
    pub b1: f64,
    pub b2: f64,
    pub r1: f64,
    pub r2: f64,
    pub drift: f64,
    pub margin: f64,
    /// `candidate faster`, `candidate slower`, `no detected effect`.  For a
    /// bandwidth or clock arm (higher is better) the ratio is inverted
    /// first, so "faster" always means the candidate does better.
    pub class: String,
}

/// Arms where a larger value is better (rates rather than times).
fn higher_is_better(arm: &str) -> bool {
    arm.starts_with("mem.") || arm.starts_with("clock.")
}

pub fn compare(
    a1: &Measurements,
    a2: &Measurements,
    b1: &Measurements,
    b2: &Measurements,
    floor: f64,
) -> Vec<ArmCompare> {
    let arms: BTreeSet<String> = a1.samples.iter().map(|s| s.arm.clone()).collect();
    let mut out = Vec::new();
    for arm in arms {
        let (Some(xa1), Some(xa2), Some(xb1), Some(xb2)) = (
            a1.median(&arm),
            a2.median(&arm),
            b1.median(&arm),
            b2.median(&arm),
        ) else {
            continue;
        };
        let inv = higher_is_better(&arm);
        let ratio = |b: f64, a: f64| if inv { a / b } else { b / a };
        let (r1, r2) = (ratio(xb1, xa1), ratio(xb2, xa2));
        let drift = (xa2 / xa1 - 1.0).abs();
        let margin = (2.0 * drift).max(floor);
        let class = if r1 < 1.0 - margin && r2 < 1.0 - margin {
            "candidate faster"
        } else if r1 > 1.0 + margin && r2 > 1.0 + margin {
            "candidate slower"
        } else {
            "no detected effect"
        };
        out.push(ArmCompare {
            arm,
            a1: xa1,
            a2: xa2,
            b1: xb1,
            b2: xb2,
            r1,
            r2,
            drift,
            margin,
            class: class.into(),
        });
    }
    out
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PairRow {
    pub run: String,
    pub batch: usize,
    pub fp64_ns: f64,
    pub fp64bf_ns: f64,
    pub ratio: f64,
}

/// `fp64bf / fp64` kernel medians at every batch size of every run.
pub fn fp64_pairs(runs: &[(String, &Measurements)]) -> Vec<PairRow> {
    field_pairs(runs, "fp64", "fp64bf")
}

/// `b / a` kernel medians (two fields of the same runs) at every batch size.
pub fn field_pairs(runs: &[(String, &Measurements)], a: &str, b_field: &str) -> Vec<PairRow> {
    let mut out = Vec::new();
    for (name, m) in runs {
        for &b in &m.params.batches {
            if let (Some(x), Some(y)) = (
                m.median(&format!("kernel.{a}.{b}")),
                m.median(&format!("kernel.{b_field}.{b}")),
            ) {
                out.push(PairRow {
                    run: name.clone(),
                    batch: b,
                    fp64_ns: x,
                    fp64bf_ns: y,
                    ratio: y / x,
                });
            }
        }
    }
    out
}

pub fn compare_csv(rows: &[ArmCompare]) -> String {
    let mut s = String::from("arm,a1,a2,b1,b2,r1,r2,drift,margin,class\n");
    for r in rows {
        s.push_str(&format!(
            "{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{}\n",
            r.arm, r.a1, r.a2, r.b1, r.b2, r.r1, r.r2, r.drift, r.margin, r.class
        ));
    }
    s
}

/// The rows as CSV, the two kernel columns named after fields `a` and `b`.
pub fn pairs_csv(rows: &[PairRow], a: &str, b: &str) -> String {
    let mut s = format!("run,batch,{a}_ns,{b}_ns,ratio\n");
    for r in rows {
        s.push_str(&format!(
            "{},{},{:.4},{:.4},{:.6}\n",
            r.run, r.batch, r.fp64_ns, r.fp64bf_ns, r.ratio
        ));
    }
    s
}
