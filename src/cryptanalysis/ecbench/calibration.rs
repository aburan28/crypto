//! Empirical check of the wall interval and admission rule using repeated A/A
//! sessions. The existing quality score remains an operational rubric.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use super::compare::{Comparison, COMPARISON_SCHEMA};

pub const CALIBRATION_SCHEMA: &str = "ecbench.calibration/v1";

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ControlClass {
    pub sessions: u64,
    /// Sessions with both measured contention channels present.
    pub telemetry_complete_sessions: u64,
    /// Sessions with at least 5% foreign CPU busy or run-queue delay.
    pub contention_evidenced_sessions: u64,
    pub intervals: u64,
    pub intervals_covering_one: u64,
    pub empirical_coverage: Option<f64>,
    pub coverage_wilson95: Option<(f64, f64)>,
    pub admitted: u64,
    pub admission_rate: Option<f64>,
    pub admission_wilson95: Option<(f64, f64)>,
    pub false_admissions: u64,
    pub false_admission_rate: Option<f64>,
    pub false_admission_wilson95: Option<(f64, f64)>,
    pub missing_interval: u64,
    pub refused: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Calibration {
    pub schema: String,
    pub quiet: ControlClass,
    pub contended: ControlClass,
    pub env_class_ids: Vec<String>,
    pub missing_session_manifests: u64,
    pub input_sha256: Vec<(String, String)>,
    pub status: String,
    pub notes: Vec<String>,
}

fn wilson95(success: u64, n: u64) -> Option<(f64, f64)> {
    if n == 0 || success > n {
        return None;
    }
    let p = success as f64 / n as f64;
    let z = 1.959963984540054;
    let z2 = z * z;
    let nf = n as f64;
    let center = (p + z2 / (2.0 * nf)) / (1.0 + z2 / nf);
    let half = z * ((p * (1.0 - p) / nf + z2 / (4.0 * nf * nf)).sqrt()) / (1.0 + z2 / nf);
    Some(((center - half).max(0.0), (center + half).min(1.0)))
}

fn class(rows: &[Comparison]) -> ControlClass {
    let mut c = ControlClass {
        sessions: rows.len() as u64,
        telemetry_complete_sessions: 0,
        contention_evidenced_sessions: 0,
        intervals: 0,
        intervals_covering_one: 0,
        empirical_coverage: None,
        coverage_wilson95: None,
        admitted: 0,
        admission_rate: None,
        admission_wilson95: None,
        false_admissions: 0,
        false_admission_rate: None,
        false_admission_wilson95: None,
        missing_interval: 0,
        refused: 0,
    };
    for row in rows {
        if let Some(q) = &row.quality {
            if q.p95_foreign_busy_percent.is_some()
                && q.p95_run_queue_percent.is_some()
                && q.contention_observed_runs >= q.verified_pairs.saturating_mul(2)
            {
                c.telemetry_complete_sessions += 1;
            }
            if q.p95_foreign_busy_percent.is_some_and(|v| v >= 5.0)
                || q.p95_run_queue_percent.is_some_and(|v| v >= 5.0)
            {
                c.contention_evidenced_sessions += 1;
            }
        }
        if let Some((lo, hi)) = row.wall.ci95 {
            c.intervals += 1;
            if lo <= 1.0 && hi >= 1.0 {
                c.intervals_covering_one += 1;
            }
            if row.wall.status == "admitted" {
                c.admitted += 1;
                if lo > 1.0 || hi < 1.0 {
                    c.false_admissions += 1;
                }
            }
        } else {
            c.missing_interval += 1;
        }
        if row.wall.status == "refused" {
            c.refused += 1;
        }
    }
    if c.intervals > 0 {
        c.empirical_coverage = Some(c.intervals_covering_one as f64 / c.intervals as f64);
        c.coverage_wilson95 = wilson95(c.intervals_covering_one, c.intervals);
    }
    if c.sessions > 0 {
        c.admission_rate = Some(c.admitted as f64 / c.sessions as f64);
        c.admission_wilson95 = wilson95(c.admitted, c.sessions);
    }
    if c.sessions > 0 {
        c.false_admission_rate = Some(c.false_admissions as f64 / c.sessions as f64);
        c.false_admission_wilson95 = wilson95(c.false_admissions, c.sessions);
    }
    c
}

fn read(
    path: &Path,
    seen: &mut BTreeSet<String>,
) -> Result<(Comparison, String, Option<String>), String> {
    let raw = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let row: Comparison = serde_json::from_slice(&raw).map_err(|e| e.to_string())?;
    if row.schema != COMPARISON_SCHEMA
        || row.a.method_id != row.b.method_id
        || row.a.method != row.b.method
    {
        return Err(format!("{}: expected an A/A comparison", path.display()));
    }
    if !row.ops.same_session || row.a.session_id != row.b.session_id {
        return Err(format!(
            "{}: A/A arms must share one session",
            path.display()
        ));
    }
    if !seen.insert(row.a.session_id.clone()) {
        return Err(format!("{}: duplicate session", path.display()));
    }
    if row.wall.status == "admitted" && row.wall.ci95.is_none() {
        return Err(format!("{}: admitted A/A has no interval", path.display()));
    }
    if let Some((lo, hi)) = row.wall.ci95 {
        if !lo.is_finite() || !hi.is_finite() || lo <= 0.0 || lo > hi {
            return Err(format!("{}: invalid wall interval", path.display()));
        }
    }
    let host_class = if let Some(dir) = path
        .parent()
        .and_then(Path::parent)
        .filter(|dir| dir.join("session.json").exists())
    {
        let session = super::runner::read_session(dir)?;
        if session.session_id != row.a.session_id {
            return Err(format!(
                "{}: source session identity mismatch",
                path.display()
            ));
        }
        if session.status != "complete" || session.records_sha256.is_none() {
            return Err(format!(
                "{}: source session is not sealed and complete",
                path.display()
            ));
        }
        let audit = super::audit::audit(dir, 0)?;
        if !audit.ok {
            return Err(format!(
                "{}: source session audit failed: {}",
                path.display(),
                audit.problems.join("; ")
            ));
        }
        let seed = row.resample_seed.ok_or_else(|| {
            format!(
                "{}: bootstrap seed is absent; regenerate the comparison",
                path.display()
            )
        })?;
        let replay =
            super::compare::compare(dir, &row.a.arm, None, &row.b.arm, row.ops.resamples, seed)?;
        let same_ci = |a: Option<(f64, f64)>, b: Option<(f64, f64)>| match (a, b) {
            (Some((al, ah)), Some((bl, bh))) => {
                (al - bl).abs() <= 1e-12 && (ah - bh).abs() <= 1e-12
            }
            (None, None) => true,
            _ => false,
        };
        if replay.comparison_id != row.comparison_id
            || replay.wall.status != row.wall.status
            || replay.wall.pairs != row.wall.pairs
            || !same_ci(replay.wall.ci95, row.wall.ci95)
            || !same_ci(replay.wall.aa_ci95, row.wall.aa_ci95)
            || replay.wall.outside_noise != row.wall.outside_noise
        {
            return Err(format!(
                "{}: comparison does not reproduce from its source session",
                path.display()
            ));
        }
        Some(session.env_class_id)
    } else {
        None
    };
    Ok((row, super::canonical::sha256_hex(&raw), host_class))
}

/// One comparison per independent session in each class. Inputs are
/// content-bound and duplicate sessions are rejected across both classes.
pub fn evaluate(quiet: &[PathBuf], contended: &[PathBuf]) -> Result<Calibration, String> {
    let mut seen = BTreeSet::new();
    let mut hashes = Vec::new();
    let mut host_classes = BTreeSet::new();
    let mut missing_session_manifests = 0;
    let mut load = |paths: &[PathBuf]| -> Result<Vec<Comparison>, String> {
        paths
            .iter()
            .map(|p| {
                let (r, hash, host_class) = read(p, &mut seen)?;
                hashes.push((p.display().to_string(), hash));
                if let Some(host_class) = host_class {
                    host_classes.insert(host_class);
                } else {
                    missing_session_manifests += 1;
                }
                Ok(r)
            })
            .collect()
    };
    let q = class(&load(quiet)?);
    let c = class(&load(contended)?);
    let mut notes = Vec::new();
    if q.sessions < 100 || c.sessions < 100 {
        notes.push("at least 100 independent sessions in each class are required before calibration is qualified".into());
    }
    if missing_session_manifests > 0 || host_classes.len() != 1 {
        notes.push("calibration requires source session manifests from one host class".into());
    }
    if q.telemetry_complete_sessions != q.sessions || c.telemetry_complete_sessions != c.sessions {
        notes.push("CPU contention channels are incomplete in some controls".into());
    }
    if c.contention_evidenced_sessions != c.sessions {
        notes.push("some contended controls lack at least 5 percent measured contention".into());
    }
    if q.contention_evidenced_sessions != 0 {
        notes.push("some quiet controls exceed the 5 percent contention screen".into());
    }
    if q.missing_interval > 0 {
        notes.push("quiet controls need an admitted interval in every session".into());
    }
    if q.coverage_wilson95.is_none_or(|(lo, _)| lo < 0.90) {
        notes.push("quiet A/A coverage has not cleared the 90 percent Wilson lower bound".into());
    }
    if q.false_admission_wilson95.is_none_or(|(_, hi)| hi > 0.10) {
        notes.push("quiet A/A false-positive upper bound is not below 10 percent".into());
    }
    if c.admission_wilson95.is_none_or(|(_, hi)| hi > 0.10) {
        notes.push("contended-session admission upper bound is not below 10 percent".into());
    }
    let status = if notes.is_empty() {
        "qualified"
    } else {
        "unqualified"
    };
    Ok(Calibration {
        schema: CALIBRATION_SCHEMA.into(),
        quiet: q,
        contended: c,
        env_class_ids: host_classes.into_iter().collect(),
        missing_session_manifests,
        input_sha256: hashes,
        status: status.into(),
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wilson_bounds_and_censoring() {
        assert_eq!(wilson95(0, 0), None);
        assert!(wilson95(95, 100).unwrap().0 < 0.95);
        assert!(wilson95(0, 100).unwrap().1 < 0.04);
        let empty = evaluate(&[], &[]).unwrap();
        assert_eq!(empty.status, "unqualified");
        assert!(empty.quiet.empirical_coverage.is_none());
        assert!(empty.contended.admission_rate.is_none());
    }
}
