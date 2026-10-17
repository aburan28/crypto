//! Two builds measured in alternation (A1, B1, A2, B2): per arm, the
//! candidate-to-baseline ratios `r_k = B_k / A_k`, the same-build drift
//! `d = |A2 / A1 − 1|`, and a classification that calls an effect only when
//! both ratios clear `max(2d, floor)` on the same side.  Also the
//! `fp64bf / fp64` kernel ratio inside each run (round-2 claim F1), and
//! threshold claims over a set of arms ([`ratio_claim`]).

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};

use super::fit::Claim;
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

/// The arms named `prefix…` whose batch size (the last dotted field, where
/// it is a number) is at least `min_batch`; arms without one always count.
fn select<'a>(rows: &'a [ArmCompare], prefix: &str, min_batch: usize) -> Vec<&'a ArmCompare> {
    rows.iter()
        .filter(|r| {
            r.arm.starts_with(prefix)
                && r.arm
                    .rsplit('.')
                    .next()
                    .and_then(|b| b.parse::<usize>().ok())
                    .is_none_or(|b| b >= min_batch)
        })
        .collect()
}

/// A claim over a set of arms from `ID:PREFIX:MIN_BATCH:T`: the largest
/// ratio max(r1, r2) is at most `T`, or with `dev`, the largest deviation
/// |r_k − 1| is.  No matching arm is an error, not a pass.
pub fn ratio_claim(rows: &[ArmCompare], spec: &str, dev: bool) -> Result<Claim, String> {
    let parts: Vec<&str> = spec.split(':').collect();
    let [id, prefix, min_batch, t] = parts[..] else {
        return Err(format!("claim {spec:?}: want ID:PREFIX:MIN_BATCH:T"));
    };
    let min_batch: usize = min_batch.parse().map_err(|e| format!("{spec}: {e}"))?;
    let t: f64 = t.parse().map_err(|e| format!("{spec}: {e}"))?;
    let arms = select(rows, prefix, min_batch);
    if arms.is_empty() {
        return Err(format!(
            "claim {id}: no arm matches {prefix:?} at batch ≥ {min_batch}"
        ));
    }
    let value = arms
        .iter()
        .flat_map(|r| [r.r1, r.r2])
        .map(|r| if dev { (r - 1.0).abs() } else { r })
        .fold(f64::MIN, f64::max);
    let what = if dev {
        "max |r_k − 1|"
    } else {
        "max(r1, r2)"
    };
    Ok(Claim::at_most(
        id,
        &format!(
            "{what} over {} arms {prefix}* (batch ≥ {min_batch})",
            arms.len()
        ),
        value,
        t,
    ))
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

#[cfg(test)]
mod tests {
    use super::*;

    fn row(arm: &str, r1: f64, r2: f64) -> ArmCompare {
        ArmCompare {
            arm: arm.into(),
            a1: 1.0,
            a2: 1.0,
            b1: r1,
            b2: r2,
            r1,
            r2,
            drift: 0.0,
            margin: 0.02,
            class: String::new(),
        }
    }

    #[test]
    fn ratio_claims_select_by_prefix_and_batch() {
        let rows = [
            row("kernel.fp64.16", 0.9, 0.95),
            row("kernel.fp64.4096", 0.5, 0.45),
            row("kernel.fp64.65536", 0.4, 0.6),
            row("kernel.fp64cm.4096", 1.03, 0.92),
            row("fp64.mul.latency", 1.04, 1.02),
        ];
        // B = 16 is below the cut; `kernel.fp64cm.` does not match `kernel.fp64.`.
        let c = ratio_claim(&rows, "P1:kernel.fp64.:4096:0.85", false).unwrap();
        assert_eq!((c.value, c.verdict.as_str()), (0.6, "PASS"));
        let c = ratio_claim(&rows, "P2:fp64.mul.latency:0:1.03", false).unwrap();
        assert_eq!((c.value, c.verdict.as_str()), (1.04, "FAIL"));
        let c = ratio_claim(&rows, "C1:kernel.fp64cm.:4096:0.10", true).unwrap();
        assert!((c.value - 0.08).abs() < 1e-12 && c.verdict == "PASS");
        assert!(ratio_claim(&rows, "X:kernel.w1.:0:1", false).is_err());
        assert!(ratio_claim(&rows, "X:kernel.fp64.:0", false).is_err());
    }
}
