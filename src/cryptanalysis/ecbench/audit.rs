//! `ecbench verify`: re-derive a session from its own files.
//!
//! The audit trusts nothing a session says about itself that it can
//! recompute:
//!
//! - every file hashes to what `session.json` says it does;
//! - the spec expands to the same spec id, workloads, execution order
//!   and seeds as `plan.json`;
//! - the host capsule's class id recomputes from its stable facts;
//! - every record is sealed, in sequence, matches its planned execution,
//!   arm, method and workload, and keeps the stderr it hashed;
//! - every recovered logarithm is re-checked here, `[k]G = Q`, and the
//!   record's status must agree with that check;
//! - with `--replay N`, N measured deterministic runs are re-executed and
//!   must reproduce their answer, total, phases, counters and factor base
//!   exactly, and their field-operation block when both the record and
//!   the replay carry one (a record written before the block existed has
//!   none, and unknown is not a disagreement).  Operation counts do not
//!   depend on the host, so a replay on any machine is an independent
//!   check of the figure.
//!
//! The receipt names every file by hash and every check by result; its
//! own SHA-256 is the replay certificate a claim cites.

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::cryptanalysis::ecbench::canonical::sha256_hex;
use crate::cryptanalysis::ecbench::host::{recompute_class, HostCapsule};
use crate::cryptanalysis::ecbench::methods::MethodSpec;
use crate::cryptanalysis::ecbench::record::{
    check_line_seal, child_main, floor_s, grade, json_roundtrip, ChildInput, GradeInput, Record,
    CHILD_INPUT_SCHEMA, GRADING_VERSION,
};
use crate::cryptanalysis::ecbench::runner::{read_record_lines, read_session, PlanDoc};
use crate::cryptanalysis::ecbench::spec::{plan, Arm, Spec};
use crate::cryptanalysis::ecbench::workload::{CurveSpec, Workload};

/// Replay every measured deterministic run.
pub const REPLAY_ALL: usize = usize::MAX;

/// What makes a workload the workload: its id and every field the id is
/// a hash of, plus the planted answer.  The registry-derived facts
/// (`registered`, `ec1`, `curve_uid`) and the constructor's text are not
/// in it, so registering a curve later never breaks an honest session.
fn workload_identity(w: &Workload) -> impl PartialEq + std::fmt::Debug + '_ {
    (
        (&w.workload_id, &w.workload_sha256, &w.target, w.planted),
        (
            &w.curve.icv1,
            w.curve.r,
            w.curve.cofactor,
            &w.curve.generator,
        ),
        (&w.target_law, w.target_seed, w.target_index),
    )
}

/// What makes an arm the arm: its name, role and method identity, not the
/// method's descriptive `entry` text.
fn arm_identity(a: &Arm) -> impl PartialEq + std::fmt::Debug + '_ {
    (
        &a.name,
        a.role,
        &a.method.id,
        &a.method.params,
        &a.method.method_id,
        &a.method.method_sha256,
    )
}

/// Equal to within a few ulps of a parsed float (records are parsed with a
/// float reader that is not always correctly rounded).
fn close(a: f64, b: f64) -> bool {
    a == b || (a - b).abs() <= 1e-12 * a.abs().max(b.abs())
}

/// Replay detail for a record whose gae carries the pre-fix rounding.
pub const LEGACY_IDENTICAL: &str =
    "identical (legacy gae rounding, reproduced from the record's solver term)";

/// The gae figures an `ic.pipeline` record written before the
/// solver-rounding fix carries, from the replay's exact phases.
///
/// Those binaries removed a wall-priced solver term `w` by subtraction,
/// so the relations phase was `(x + w) - w` and the total
/// `(fb + prep + (x + w) + la + ver) - w`: `x` rounded to `w`'s binade.
/// The rounding is a function of `w`, which the record carries
/// (`detail.decomposition.solver.gae`), so it is reproduced here from the
/// record rather than from the replay host's clock.  `w` is accepted only
/// when it also reproduces the record's own pre-removal total
/// (`detail.framework_total_gae_before_solver_removal`), so the record
/// cannot pick a `w` that absorbs a real difference; what stays unchecked
/// is below one ulp of `x + w`.
pub fn legacy_solver_rounding(
    exact: &[f64],
    detail: &serde_json::Value,
) -> Option<(Vec<f64>, f64)> {
    let solver = detail.pointer("/decomposition/solver")?;
    let w = solver.get("gae")?.as_f64()?;
    if exact.len() != 5
        || w.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater)
        || solver.get("priced_by")?.as_str()? == "pinned"
    {
        return None;
    }
    let with_w = exact[2] + w;
    let before = exact[0] + exact[1] + with_w + exact[3] + exact[4];
    let recorded = detail
        .get("framework_total_gae_before_solver_removal")?
        .as_f64()?;
    if json_roundtrip(before).to_bits() != recorded.to_bits() {
        return None;
    }
    let mut phases = exact.to_vec();
    phases[2] = with_w - w;
    Some((phases, before - w))
}

pub const AUDIT_SCHEMA: &str = "ecbench.audit/v1";

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Replay {
    pub seq: u64,
    pub record_id: String,
    pub reproduced: bool,
    pub detail: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuditReport {
    pub schema: String,
    pub session_id: String,
    /// The session's own status: `complete`, or `interrupted` when that
    /// was accepted.
    #[serde(default)]
    pub session_status: String,
    pub ok: bool,
    pub problems: Vec<String>,
    pub checks: BTreeMap<String, bool>,
    pub records: u64,
    pub verified_records: u64,
    /// Every run's level was recomputed from its observations (the
    /// session was graded under this binary's rules).
    #[serde(default)]
    pub regraded: bool,
    pub replays: Vec<Replay>,
    /// SHA-256 of every file the audit read.
    pub files: BTreeMap<String, String>,
    pub auditor_binary_sha256: Option<String>,
    /// The host class the audit ran on: a replay on a class other than
    /// the session's is an independent one.
    #[serde(default)]
    pub auditor_env_class_id: Option<String>,
    pub audited_unix_ms: u128,
}

fn file_sha(dir: &Path, name: &str, files: &mut BTreeMap<String, String>) -> Option<String> {
    let bytes = std::fs::read(dir.join(name)).ok()?;
    let h = sha256_hex(&bytes);
    files.insert(name.to_string(), h.clone());
    Some(h)
}

/// Audit `dir`; replay up to `replay` deterministic measured runs
/// ([`REPLAY_ALL`] for every one).
pub fn audit(dir: &Path, replay: usize) -> Result<AuditReport, String> {
    audit_with(dir, replay, false)
}

/// [`audit`], optionally accepting an `interrupted` session: it is held
/// to every integrity check (hashes, seals, plan, answers, derived
/// figures), with "every planned execution recorded" replaced by "every
/// record it says it wrote".  A session that claims `complete` gets no
/// such relaxation.
pub fn audit_with(
    dir: &Path,
    replay: usize,
    allow_interrupted: bool,
) -> Result<AuditReport, String> {
    let mut problems: Vec<String> = Vec::new();
    let mut checks: BTreeMap<String, bool> = BTreeMap::new();
    let mut files = BTreeMap::new();
    let session = read_session(dir)?;
    file_sha(dir, "session.json", &mut files);
    file_sha(dir, "spec.json", &mut files);

    let mut check = |name: &str, ok: bool, why: String, problems: &mut Vec<String>| {
        checks.insert(name.to_string(), ok);
        if !ok {
            problems.push(format!("{name}: {why}"));
        }
    };

    // Files against the session's hashes.
    let host_sha = file_sha(dir, "host.json", &mut files);
    check(
        "host_hash",
        host_sha.as_deref() == Some(&session.host_sha256[..]),
        "host.json does not hash to session.host_sha256".into(),
        &mut problems,
    );
    let plan_sha = file_sha(dir, "plan.json", &mut files);
    check(
        "plan_hash",
        plan_sha.as_deref() == Some(&session.plan_sha256[..]),
        "plan.json does not hash to session.plan_sha256".into(),
        &mut problems,
    );
    let rec_sha = file_sha(dir, "records.jsonl", &mut files);
    let complete = session.status == "complete";
    let accepted_interruption = allow_interrupted && session.status == "interrupted";
    check(
        "records_hash",
        !(complete || accepted_interruption) || rec_sha == session.records_sha256,
        "records.jsonl does not hash to session.records_sha256".into(),
        &mut problems,
    );
    check(
        "session_complete",
        complete || accepted_interruption,
        format!("session status is `{}`", session.status),
        &mut problems,
    );

    // The host class.
    let host: HostCapsule = serde_json::from_str(
        &std::fs::read_to_string(dir.join("host.json")).map_err(|e| e.to_string())?,
    )
    .map_err(|e| format!("host.json: {e}"))?;
    let class = recompute_class(&host)?;
    check(
        "env_class",
        class == host.env_class_id && class == session.env_class_id,
        format!(
            "recomputed {class}, recorded {} / {}",
            host.env_class_id, session.env_class_id
        ),
        &mut problems,
    );

    // The spec re-expands to the plan.
    let spec_text = std::fs::read_to_string(dir.join("spec.json")).map_err(|e| e.to_string())?;
    let p = plan(Spec::from_json(&spec_text)?)?;
    let doc: PlanDoc = serde_json::from_str(
        &std::fs::read_to_string(dir.join("plan.json")).map_err(|e| e.to_string())?,
    )
    .map_err(|e| format!("plan.json: {e}"))?;
    check(
        "spec_id",
        p.spec_id == session.spec_id && p.spec_id == doc.spec_id,
        format!("spec re-expands to {}", p.spec_id),
        &mut problems,
    );
    check(
        "workloads",
        p.workloads.len() == doc.workloads.len()
            && p.workloads
                .iter()
                .zip(&doc.workloads)
                .all(|(a, b)| workload_identity(a) == workload_identity(b))
            && doc
                .workloads
                .iter()
                .all(|w| w.recompute_id().ok().as_ref() == Some(&w.workload_id)),
        "the spec rebuilds different workloads, or a workload id does not recompute".into(),
        &mut problems,
    );
    check(
        "executions",
        p.executions == doc.executions,
        "the spec expands to a different execution order or seeds".into(),
        &mut problems,
    );
    check(
        "arms",
        p.arms.len() == doc.arms.len()
            && p.arms
                .iter()
                .zip(&doc.arms)
                .all(|(a, b)| arm_identity(a) == arm_identity(b)),
        "the spec resolves to different methods".into(),
        &mut problems,
    );

    // Records.
    let lines = read_record_lines(dir)?;
    let records: Vec<Record> = lines.iter().map(|(_, r)| r.clone()).collect();
    check(
        "record_count",
        if complete {
            records.len() as u64 == session.records_written && records.len() == p.executions.len()
        } else {
            records.len() as u64 == session.records_written && records.len() <= p.executions.len()
        },
        format!(
            "{} records, {} written, {} planned",
            records.len(),
            session.records_written,
            p.executions.len()
        ),
        &mut problems,
    );
    let mut record_problems = Vec::new();
    let mut verified = 0u64;
    for (i, r) in records.iter().enumerate() {
        let tag = format!("record seq {}", r.seq);
        if r.seq != i as u64 {
            record_problems.push(format!("{tag}: out of sequence at line {}", i + 1));
        }
        if !check_line_seal(&lines[i].0, &r.record_id) {
            record_problems.push(format!("{tag}: content does not hash to {}", r.record_id));
        }
        if r.session_id != session.session_id {
            record_problems.push(format!(
                "{tag}: session {} is not {}",
                r.session_id, session.session_id
            ));
        }
        let Some(ex) = p.executions.get(i) else {
            record_problems.push(format!("{tag}: no planned execution"));
            continue;
        };
        let arm = &p.arms[ex.arm];
        let w = &p.workloads[ex.workload];
        if r.round != ex.round
            || r.warmup != ex.warmup
            || r.algorithm_seed != ex.algorithm_seed.to_string()
        {
            record_problems.push(format!(
                "{tag}: round, warm-up or seed differs from the plan"
            ));
        }
        if r.arm != arm.name
            || r.method.method_id != arm.method.method_id
            || r.method.method_sha256 != arm.method.method_sha256
        {
            record_problems.push(format!("{tag}: arm or method differs from the plan"));
        }
        if workload_identity(&r.workload) != workload_identity(w) {
            record_problems.push(format!("{tag}: workload differs from the plan"));
        }
        // Every derived figure must follow from the record's own counts:
        // a hand-edited S, ratio or flag fails here even if resealed.
        let floor = floor_s(r.boundaries.automorphisms_available);
        let sqrt_r = (w.curve.r as f64).sqrt();
        let derived_ok = r.boundaries.automorphisms_available == w.curve.automorphisms_available
            && close(r.boundaries.floor_s, floor)
            && match (r.cost.total_gae, r.cost.s, r.boundaries.ratio_to_floor) {
                (Some(t), Some(sv), Some(q)) => close(sv, t / sqrt_r) && close(q, sv / floor),
                (None, None, None) => true,
                _ => false,
            }
            && r.cost.lower_bound == !r.cost.unpriced.is_empty()
            && (r.cost.total_gae.is_none()
                || r.cost.deterministic == r.cost.nondeterminism.is_empty())
            && r.run_id
                .starts_with(&format!("{}{}R", r.method.method_id, w.workload_id));
        if !derived_ok {
            record_problems.push(format!(
                "{tag}: S, floor, ratio, lower-bound flag or run id does not follow from the record"
            ));
        }
        // Re-check the answer here.
        let inst = p.instance(ex.workload);
        let rec: Option<u128> = r.outcome.recovered.as_deref().and_then(|s| s.parse().ok());
        let target_ok = rec.map(|k| inst.mul_generator_hex(k).as_ref() == Some(&w.target));
        let planted_ok = w.planted.and_then(|p| rec.map(|k| k == p));
        if target_ok != r.outcome.matches_target || planted_ok != r.outcome.matches_planted {
            record_problems.push(format!("{tag}: verification does not reproduce"));
        }
        let should_verify = target_ok == Some(true) && planted_ok != Some(false);
        if (r.outcome.status == "verified") != should_verify {
            record_problems.push(format!(
                "{tag}: status `{}` contradicts the check",
                r.outcome.status
            ));
        }
        if r.outcome.status == "verified" {
            verified += 1;
        }
        if let Some(h) = &r.outcome.stderr_sha256 {
            let name = format!("exec/{:06}.stderr", r.seq);
            if file_sha(dir, &name, &mut files).as_deref() != Some(h) {
                record_problems.push(format!("{tag}: {name} missing or altered"));
            }
        }
    }
    // Regrade every run from its own observations, when the session was
    // graded under the rules this binary applies.
    let regraded = session.grading_version == GRADING_VERSION;
    if regraded {
        for r in &records {
            let g = grade(&GradeInput {
                plan: session.cpu_plan.as_ref(),
                preflight: session.preflight.as_ref(),
                eviction: session.eviction.as_ref(),
                runner_on_run_cpu: session.runner_on_run_cpu,
                capsule: &host,
                thresholds: &session.options.thresholds,
                placement: &r.isolation.placement,
                before: &r.isolation.before,
                after: &r.isolation.after,
                timing: &r.time,
                hz: session.ticks_per_second as f64,
            });
            if g.level != r.isolation.level || g.blockers != r.isolation.blockers {
                record_problems.push(format!(
                    "record seq {}: graded {} with {} blockers, its observations give {} with {}",
                    r.seq,
                    r.isolation.level.name(),
                    r.isolation.blockers.len(),
                    g.level.name(),
                    g.blockers.len()
                ));
            }
        }
    }
    let mut seen = std::collections::BTreeSet::new();
    for r in &records {
        if !seen.insert(r.run_id.as_str()) {
            record_problems.push(format!("record seq {}: run id {} repeats", r.seq, r.run_id));
        }
    }
    check(
        "records",
        record_problems.is_empty(),
        format!("{} record problems", record_problems.len()),
        &mut problems,
    );
    problems.extend(record_problems);

    // Replays: evenly spaced over the measured deterministic verified runs.
    let candidates: Vec<&Record> = records
        .iter()
        .filter(|r| r.counts() && r.cost.deterministic)
        .collect();
    let mut replays = Vec::new();
    if replay > 0 && !candidates.is_empty() {
        let n = replay.min(candidates.len());
        for j in 0..n {
            let r = candidates[j * candidates.len() / n];
            let ex = &p.executions[r.seq as usize];
            let w = &p.workloads[ex.workload];
            let input = ChildInput {
                schema: CHILD_INPUT_SCHEMA.into(),
                curve: CurveSpec::explicit(p.instance(ex.workload))
                    .unwrap_or_else(|| w.curve_spec.clone()),
                target_seed: w.target_seed,
                target_index: w.target_index,
                target_kind: w.kind(),
                expected_workload_id: w.workload_id.clone(),
                method: MethodSpec {
                    id: r.method.id.clone(),
                    params: r.method.params.clone(),
                },
                expected_method_id: r.method.method_id.clone(),
                algorithm_seed: ex.algorithm_seed,
            };
            let out = child_main(&serde_json::to_string(&input).map_err(|e| e.to_string())?);
            let (ok, detail) = match out.report {
                None => (
                    false,
                    format!("replay failed: {}", out.error.unwrap_or_default()),
                ),
                Some(rep) => {
                    let mut diffs = Vec::new();
                    if rep.recovered.map(|k| k.to_string()) != r.outcome.recovered {
                        diffs.push("recovered");
                    }
                    // Integer counts compare exactly; a float compares after
                    // the same serialise-parse round trip the record's took.
                    let same_gae = |phases: &[f64], total: f64| {
                        (
                            Some(json_roundtrip(total).to_bits())
                                == r.cost.total_gae.map(f64::to_bits),
                            phases
                                .iter()
                                .map(|g| json_roundtrip(*g).to_bits())
                                .eq(r.phases.iter().map(|p| p.gae.to_bits())),
                        )
                    };
                    let exact: Vec<f64> = rep.phases.iter().map(|p| p.gae).collect();
                    let (mut total_ok, mut phases_ok) = same_gae(&exact, rep.total_gae);
                    let mut legacy = false;
                    if !(total_ok && phases_ok) {
                        if let Some((ph, tot)) = legacy_solver_rounding(&exact, &r.detail) {
                            if same_gae(&ph, tot) == (true, true) {
                                (total_ok, phases_ok, legacy) = (true, true, true);
                            }
                        }
                    }
                    if !total_ok {
                        diffs.push("total_gae");
                    }
                    let ints = |ps: &[crate::cryptanalysis::ecbench::methods::PhaseRecord]| {
                        ps.iter()
                            .map(|p| {
                                (
                                    p.name.clone(),
                                    p.adds,
                                    p.doubles,
                                    p.scalar_mults,
                                    p.native.clone(),
                                )
                            })
                            .collect::<Vec<_>>()
                    };
                    if ints(&rep.phases) != ints(&r.phases) {
                        diffs.push("phase counts");
                    }
                    if !phases_ok {
                        diffs.push("phase gae");
                    }
                    if rep.counters != r.counters {
                        diffs.push("counters");
                    }
                    if rep.factor_base != r.factor_base {
                        diffs.push("factor_base");
                    }
                    // Field operations compare only when both sides carry
                    // them: a committed record written before they were
                    // counted has none, and unknown is not a disagreement.
                    if let (Some(a), Some(b)) = (rep.field_ops, r.field_ops) {
                        if a != b {
                            diffs.push("field_ops");
                        }
                    }
                    if rep.unpriced != r.cost.unpriced
                        || rep.deterministic != r.cost.deterministic
                        || Some(rep.automorphisms_used) != r.boundaries.automorphisms_used
                    {
                        diffs.push("unpriced, determinism or automorphisms");
                    }
                    (
                        diffs.is_empty(),
                        if diffs.is_empty() && legacy {
                            LEGACY_IDENTICAL.into()
                        } else if diffs.is_empty() {
                            "identical".into()
                        } else {
                            format!("differs in {}", diffs.join(", "))
                        },
                    )
                }
            };
            if !ok {
                problems.push(format!("replay of seq {}: {detail}", r.seq));
            }
            replays.push(Replay {
                seq: r.seq,
                record_id: r.record_id.clone(),
                reproduced: ok,
                detail,
            });
        }
        checks.insert("replays".into(), replays.iter().all(|r| r.reproduced));
    }

    let auditor_binary_sha256 = std::env::current_exe()
        .ok()
        .and_then(|p| std::fs::read(p).ok())
        .map(|b| sha256_hex(&b));
    Ok(AuditReport {
        schema: AUDIT_SCHEMA.into(),
        session_status: session.status.clone(),
        session_id: session.session_id,
        ok: problems.is_empty(),
        problems,
        checks,
        records: records.len() as u64,
        verified_records: verified,
        regraded,
        replays,
        files,
        auditor_binary_sha256,
        auditor_env_class_id: crate::cryptanalysis::ecbench::host::capture()
            .ok()
            .map(|c| c.env_class_id),
        audited_unix_ms: std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_millis())
            .unwrap_or(0),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn detail(w: f64, before: f64, priced_by: &str) -> serde_json::Value {
        json!({
            "decomposition": {"solver": {"gae": w, "priced_by": priced_by}},
            "framework_total_gae_before_solver_removal": before,
        })
    }

    /// The figures of the committed koblitz session's seq 4
    /// (`sub6-da-buch`): the pre-fix rounding is a function of the
    /// record's own solver term, so it reproduces exactly.
    #[test]
    fn legacy_rounding_is_reproduced_from_the_record_alone() {
        let x = 611.5963;
        let exact = [3.223395, 0.0, x, 2.721705, 15.0];
        let w = 91449423.0;
        let before = exact[0] + exact[1] + (x + w) + exact[3] + exact[4];
        let (ph, tot) = legacy_solver_rounding(&exact, &detail(w, before, "measured")).unwrap();
        assert_eq!(ph[2], (x + w) - w);
        assert_eq!(tot, before - w);
        assert_eq!(ph[..2], exact[..2]);
        assert_eq!(ph[3..], exact[3..]);
        // Another host's clock lands `w` in another binade: the old
        // subtraction then disagrees with itself, the exact figure does not.
        let w2 = 29295605.0;
        assert_ne!((x + w).to_bits(), x.to_bits());
        assert_ne!(((x + w) - w).to_bits(), ((x + w2) - w2).to_bits());
    }

    #[test]
    fn legacy_rounding_needs_the_records_own_total_and_a_wall_price() {
        let exact = [3.223395, 0.0, 611.5963, 2.721705, 15.0];
        let w = 91449423.0;
        let before = exact[0] + exact[1] + (exact[2] + w) + exact[3] + exact[4];
        // A solver term that does not reproduce the recorded pre-removal
        // total is not the one the record was written with.
        assert!(legacy_solver_rounding(&exact, &detail(w * 2.0, before, "measured")).is_none());
        assert!(legacy_solver_rounding(&exact, &detail(w, before + 1.0, "measured")).is_none());
        // Pinned solver work was never subtracted.
        assert!(legacy_solver_rounding(&exact, &detail(w, before, "pinned")).is_none());
        assert!(legacy_solver_rounding(&exact, &json!({})).is_none());
    }
}
