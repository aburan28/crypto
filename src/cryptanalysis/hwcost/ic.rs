//! An index-calculus pipeline priced in wall time and dollars next to its
//! matched rho, from a frozen `ecbench` campaign (one session directory per
//! curve, `records.jsonl` each, and the campaign's AGENTS.md §5a band file).
//!
//! Per curve, for the IC arm and the matched rho arm:
//!
//! - `P(w)`, the GAE with unpriced native work priced at `w ∈ {0, 0.1, 1}`
//!   GAE per operation, is the band's `mean_s[w] · √r`;
//! - the measured wall time per solve is the mean `time.solve_wall_ns` of
//!   the counted (non-warm-up), verified runs, per phase as well;
//! - `t(w) = wall_rho / P_rho(w)` is what one priced GAE costs in wall time
//!   on the campaign host, from the matched rho, and `P_IC(w) · t(w)` is the
//!   IC wall time the band predicts at that weight;
//! - `w*` solves `P_IC(w) / P_rho(w) = wall_IC / wall_rho`: the weight per
//!   unpriced native operation that wall time actually charged.
//!
//! Dollars are machine-hours at a stated price; each solve ran on one core
//! of the campaign's 4-vCPU host.  Exponents fitted over a declared size tier
//! extrapolate wall time and `P(w)` to larger groups, and IC is ranked
//! against rho at a size only where wall time and every band weight agree.

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use super::fit::linfit;

pub const WEIGHTS: [f64; 3] = [0.0, 0.1, 1.0];

#[derive(Deserialize)]
struct BandFile {
    rows: Vec<BandRow>,
}

#[derive(Clone, Debug, Deserialize)]
struct BandRow {
    curve: String,
    log2_r: f64,
    arm: String,
    verified: u64,
    /// `null` where the arm is unrankable at that weight (an unpriced
    /// counter without a count).
    mean_s: Vec<Option<f64>>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PhaseTime {
    pub name: String,
    pub gae: f64,
    pub wall_ms: f64,
}

/// One arm on one curve.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ArmCell {
    pub arm: String,
    pub runs: usize,
    /// `P(w)` at [`WEIGHTS`].
    pub p: [f64; 3],
    pub wall_ms_mean: f64,
    pub wall_ms_median: f64,
    pub phases: Vec<PhaseTime>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Pair {
    pub curve: String,
    pub session: String,
    pub log2_r: f64,
    pub ic: ArmCell,
    pub rho: ArmCell,
    /// Wall ns per priced GAE of the matched rho, at each weight.
    pub ns_per_gae: [f64; 3],
    /// IC wall time the band predicts at each weight, ms.
    pub ic_wall_pred_ms: [f64; 3],
    pub err_pred: [f64; 3],
    pub wall_ratio: f64,
    pub band_ratio: [f64; 3],
    /// Inside the band's ratio range.
    pub wall_ratio_in_band: bool,
    pub w_star: f64,
    pub usd_per_solve_ic: f64,
    pub usd_per_solve_rho: f64,
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len().max(1) as f64
}

fn median(mut v: Vec<f64>) -> f64 {
    super::measure::median(&mut v).unwrap_or(f64::NAN)
}

/// The counted, verified runs of `arm` in a session's records.
fn arm_runs(
    records: &Path,
    arm: &str,
) -> Result<(Vec<f64>, BTreeMap<String, (Vec<f64>, Vec<f64>)>, String), String> {
    let text =
        std::fs::read_to_string(records).map_err(|e| format!("{}: {e}", records.display()))?;
    let mut walls = Vec::new();
    let mut phases: BTreeMap<String, (Vec<f64>, Vec<f64>)> = BTreeMap::new();
    let mut curve = String::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let v: serde_json::Value =
            serde_json::from_str(line).map_err(|e| format!("{}: {e}", records.display()))?;
        if v["arm"] != arm
            || v["warmup"].as_bool().unwrap_or(false)
            || v["outcome"]["status"] != "verified"
        {
            continue;
        }
        let wall = v["time"]["solve_wall_ns"]
            .as_f64()
            .ok_or_else(|| format!("{}: a run of {arm} has no solve_wall_ns", records.display()))?;
        walls.push(wall / 1e6);
        if curve.is_empty() {
            curve = v["workload"]["curve"]["slug"]
                .as_str()
                .unwrap_or("")
                .to_string();
        }
        for ph in v["phases"].as_array().into_iter().flatten() {
            let e = phases
                .entry(ph["name"].as_str().unwrap_or("?").to_string())
                .or_default();
            e.0.push(ph["gae"].as_f64().unwrap_or(0.0));
            e.1.push(ph["wall_ns"].as_f64().unwrap_or(0.0) / 1e6);
        }
    }
    Ok((walls, phases, curve))
}

fn cell(
    band: &[BandRow],
    curve: &str,
    arm: &str,
    records: &Path,
) -> Result<Option<(ArmCell, f64)>, String> {
    let Some(b) = band.iter().find(|r| r.curve == curve && r.arm == arm) else {
        return Ok(None);
    };
    let s: Vec<f64> = b.mean_s.iter().flatten().copied().collect();
    if s.len() != 3 {
        return Err(format!(
            "{curve} {arm}: band row has {} priced weights of 3 (unrankable at some weight)",
            s.len()
        ));
    }
    let (walls, phases, _) = arm_runs(records, arm)?;
    if walls.len() as u64 != b.verified {
        return Err(format!(
            "{curve} {arm}: {} counted verified runs in the records, {} in the band",
            walls.len(),
            b.verified
        ));
    }
    let root = (b.log2_r / 2.0).exp2();
    Ok(Some((
        ArmCell {
            arm: arm.into(),
            runs: walls.len(),
            p: [s[0] * root, s[1] * root, s[2] * root],
            wall_ms_mean: mean(&walls),
            wall_ms_median: median(walls.clone()),
            phases: phases
                .into_iter()
                .map(|(name, (g, w))| PhaseTime {
                    name,
                    gae: mean(&g),
                    wall_ms: mean(&w),
                })
                .collect(),
        },
        b.log2_r,
    )))
}

/// Pair the IC arm with its matched rho on every session's curve.
pub fn pairs(
    band_path: &Path,
    sessions: &Path,
    ic_arm: &str,
    rho_arm: &str,
    usd_per_core_hour: f64,
) -> Result<Vec<Pair>, String> {
    let s =
        std::fs::read_to_string(band_path).map_err(|e| format!("{}: {e}", band_path.display()))?;
    let band: BandFile =
        serde_json::from_str(&s).map_err(|e| format!("{}: {e}", band_path.display()))?;
    let mut dirs: Vec<_> = std::fs::read_dir(sessions)
        .map_err(|e| format!("{}: {e}", sessions.display()))?
        .filter_map(|d| d.ok().map(|d| d.path()))
        .filter(|p| p.join("records.jsonl").exists())
        .collect();
    dirs.sort();
    let mut out = Vec::new();
    for dir in dirs {
        let records = dir.join("records.jsonl");
        let (_, _, curve) = arm_runs(&records, ic_arm)?;
        if curve.is_empty() {
            continue;
        }
        let (Some((ic, log2_r)), Some((rho, _))) = (
            cell(&band.rows, &curve, ic_arm, &records)?,
            cell(&band.rows, &curve, rho_arm, &records)?,
        ) else {
            continue;
        };
        let ns_per_gae = [0, 1, 2].map(|i| rho.wall_ms_mean * 1e6 / rho.p[i]);
        let pred = [0, 1, 2].map(|i| ic.p[i] * ns_per_gae[i] / 1e6);
        let wall_ratio = ic.wall_ms_mean / rho.wall_ms_mean;
        let band_ratio = [0, 1, 2].map(|i| ic.p[i] / rho.p[i]);
        let (lo, hi) = band_ratio
            .iter()
            .fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)));
        // P(w) = P(0) + w·Q, linear in w.
        let (q_ic, q_rho) = (ic.p[2] - ic.p[0], rho.p[2] - rho.p[0]);
        let denom = q_ic - wall_ratio * q_rho;
        let w_star = if denom.abs() > 0.0 {
            (wall_ratio * rho.p[0] - ic.p[0]) / denom
        } else {
            f64::NAN
        };
        let usd = |ms: f64| ms / 3.6e6 * usd_per_core_hour;
        out.push(Pair {
            curve: curve.clone(),
            session: dir.file_name().unwrap().to_string_lossy().into_owned(),
            log2_r,
            ns_per_gae,
            ic_wall_pred_ms: pred,
            err_pred: [0, 1, 2].map(|i| pred[i] / ic.wall_ms_mean - 1.0),
            wall_ratio,
            band_ratio,
            wall_ratio_in_band: wall_ratio >= lo && wall_ratio <= hi,
            w_star,
            usd_per_solve_ic: usd(ic.wall_ms_mean),
            usd_per_solve_rho: usd(rho.wall_ms_mean),
            ic,
            rho,
        });
    }
    out.sort_by(|a, b| a.log2_r.partial_cmp(&b.log2_r).unwrap());
    Ok(out)
}

/// `log₂ y = a + α · log₂ r` over the tier.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Fit {
    pub quantity: String,
    pub arm: String,
    pub log2_c: f64,
    pub alpha: f64,
    pub points: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Projection {
    pub log2_r: f64,
    pub ic_wall_core_years: f64,
    pub rho_wall_core_years: f64,
    pub ic_usd: f64,
    pub rho_usd: f64,
    pub wall_ratio: f64,
    pub band_ratio: [f64; 3],
    /// `IC > rho`, `IC < rho`, or `unranked` when wall time and the band
    /// weights disagree on the order.
    pub rank: String,
}

pub fn fits(pairs: &[Pair], tier: (f64, f64)) -> Result<Vec<Fit>, String> {
    let sel: Vec<&Pair> = pairs
        .iter()
        .filter(|p| p.log2_r >= tier.0 && p.log2_r <= tier.1)
        .collect();
    if sel.len() < 2 {
        return Err(format!("{} sizes in the tier, need at least 2", sel.len()));
    }
    let xs: Vec<f64> = sel.iter().map(|p| p.log2_r).collect();
    let mut out = Vec::new();
    let mut push = |quantity: &str, arm: &str, ys: Vec<f64>| -> Result<(), String> {
        let (a, b) = linfit(&xs, &ys).ok_or("degenerate fit")?;
        out.push(Fit {
            quantity: quantity.into(),
            arm: arm.into(),
            log2_c: a,
            alpha: b,
            points: ys.len(),
        });
        Ok(())
    };
    for (arm, get) in [("ic", 0usize), ("rho", 1)] {
        let c = |p: &&Pair| {
            if get == 0 {
                p.ic.clone()
            } else {
                p.rho.clone()
            }
        };
        push(
            "wall_ms",
            arm,
            sel.iter().map(|p| c(p).wall_ms_mean.log2()).collect(),
        )?;
        for (i, w) in WEIGHTS.iter().enumerate() {
            push(
                &format!("P(w={w})"),
                arm,
                sel.iter().map(|p| c(p).p[i].log2()).collect(),
            )?;
        }
    }
    Ok(out)
}

pub fn project(fits: &[Fit], targets: &[f64], usd_per_core_hour: f64) -> Vec<Projection> {
    let at = |q: &str, arm: &str, x: f64| {
        let f = fits
            .iter()
            .find(|f| f.quantity == q && f.arm == arm)
            .expect("fit");
        (f.log2_c + f.alpha * x).exp2()
    };
    targets
        .iter()
        .map(|&x| {
            let (ic_ms, rho_ms) = (at("wall_ms", "ic", x), at("wall_ms", "rho", x));
            let band = [0, 1, 2].map(|i| {
                let q = format!("P(w={})", WEIGHTS[i]);
                at(&q, "ic", x) / at(&q, "rho", x)
            });
            let all = [ic_ms / rho_ms, band[0], band[1], band[2]];
            let rank = if all.iter().all(|&r| r > 1.0) {
                "IC > rho"
            } else if all.iter().all(|&r| r < 1.0) {
                "IC < rho"
            } else {
                "unranked"
            };
            let years = |ms: f64| ms / 1e3 / (365.25 * 86400.0);
            Projection {
                log2_r: x,
                ic_wall_core_years: years(ic_ms),
                rho_wall_core_years: years(rho_ms),
                ic_usd: ic_ms / 3.6e6 * usd_per_core_hour,
                rho_usd: rho_ms / 3.6e6 * usd_per_core_hour,
                wall_ratio: ic_ms / rho_ms,
                band_ratio: band,
                rank: rank.into(),
            }
        })
        .collect()
}

/// The protocol's D1 quantities over the paired curves.
pub fn summary(pairs: &[Pair]) -> Vec<super::fit::Claim> {
    use super::fit::Claim;
    let inside = pairs.iter().filter(|p| p.wall_ratio_in_band).count() as f64;
    let mut ws: Vec<f64> = pairs
        .iter()
        .map(|p| p.w_star)
        .filter(|w| w.is_finite())
        .collect();
    let w_med = super::measure::median(&mut ws).unwrap_or(f64::NAN);
    vec![
        Claim::report("D1.in_band", &format!("curves (of {}) whose measured IC/rho wall ratio lies inside the band's ratio range", pairs.len()), inside),
        Claim::report("D1.w_star_median", "median implied weight w* (GAE per unpriced native operation)", w_med),
    ]
}

/// Projection in GAE (the fitted `P(w)` exponents) converted to time at
/// `ns_per_gae[w]`, the matched rho's wall ns per priced GAE at the largest
/// measured size.  Unlike the fitted wall-time exponents this carries no
/// per-solve floor into the slope.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GaeProjection {
    pub log2_r: f64,
    pub w: f64,
    pub ns_per_gae: f64,
    pub ic_core_years: f64,
    pub rho_core_years: f64,
    pub ic_usd: f64,
    pub rho_usd: f64,
    pub ratio: f64,
}

pub fn project_gae(
    fits: &[Fit],
    ns_per_gae: [f64; 3],
    targets: &[f64],
    usd_per_core_hour: f64,
) -> Vec<GaeProjection> {
    let at = |q: &str, arm: &str, x: f64| {
        let f = fits
            .iter()
            .find(|f| f.quantity == q && f.arm == arm)
            .expect("fit");
        (f.log2_c + f.alpha * x).exp2()
    };
    let mut out = Vec::new();
    for &x in targets {
        for (i, &w) in WEIGHTS.iter().enumerate() {
            let q = format!("P(w={w})");
            let (ic, rho) = (at(&q, "ic", x), at(&q, "rho", x));
            let s = |g: f64| g * ns_per_gae[i] * 1e-9;
            out.push(GaeProjection {
                log2_r: x,
                w,
                ns_per_gae: ns_per_gae[i],
                ic_core_years: s(ic) / (365.25 * 86400.0),
                rho_core_years: s(rho) / (365.25 * 86400.0),
                ic_usd: s(ic) / 3600.0 * usd_per_core_hour,
                rho_usd: s(rho) / 3600.0 * usd_per_core_hour,
                ratio: ic / rho,
            });
        }
    }
    out
}

pub fn pairs_csv(pairs: &[Pair]) -> String {
    let mut s = String::from("session,curve,log2_r,ic_arm,rho_arm,ic_runs,rho_runs,ic_wall_ms,rho_wall_ms,ic_p_w0,ic_p_w01,ic_p_w1,rho_p_w0,rho_p_w01,rho_p_w1,ns_per_gae_w0,ns_per_gae_w01,ns_per_gae_w1,ic_pred_ms_w0,ic_pred_ms_w01,ic_pred_ms_w1,err_w0,err_w01,err_w1,wall_ratio,band_ratio_w0,band_ratio_w01,band_ratio_w1,in_band,w_star,usd_ic,usd_rho\n");
    for p in pairs {
        s.push_str(&format!(
            "{},{},{:.4},{},{},{},{},{:.4},{:.4},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.3},{:.3},{:.3},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{},{:.4},{:.3e},{:.3e}\n",
            p.session, p.curve, p.log2_r, p.ic.arm, p.rho.arm, p.ic.runs, p.rho.runs, p.ic.wall_ms_mean, p.rho.wall_ms_mean,
            p.ic.p[0], p.ic.p[1], p.ic.p[2], p.rho.p[0], p.rho.p[1], p.rho.p[2],
            p.ns_per_gae[0], p.ns_per_gae[1], p.ns_per_gae[2], p.ic_wall_pred_ms[0], p.ic_wall_pred_ms[1], p.ic_wall_pred_ms[2],
            p.err_pred[0], p.err_pred[1], p.err_pred[2], p.wall_ratio, p.band_ratio[0], p.band_ratio[1], p.band_ratio[2],
            p.wall_ratio_in_band, p.w_star, p.usd_per_solve_ic, p.usd_per_solve_rho
        ));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn projection_ranks_only_when_wall_and_band_agree() {
        let f = |q: &str, arm: &str, c: f64, a: f64| Fit {
            quantity: q.into(),
            arm: arm.into(),
            log2_c: c,
            alpha: a,
            points: 3,
        };
        let mut fits = vec![f("wall_ms", "ic", 0.0, 0.6), f("wall_ms", "rho", 0.0, 0.5)];
        for w in WEIGHTS {
            fits.push(f(&format!("P(w={w})"), "ic", 1.0, 0.6));
            fits.push(f(&format!("P(w={w})"), "rho", 0.0, 0.5));
        }
        let p = project(&fits, &[40.0], 1.0);
        assert_eq!(p[0].rank, "IC > rho");
        assert!((p[0].wall_ratio - 4.0f64.exp2()).abs() < 1e-9);
        fits[0].log2_c = -10.0;
        assert_eq!(project(&fits, &[40.0], 1.0)[0].rank, "unranked");
    }
}
