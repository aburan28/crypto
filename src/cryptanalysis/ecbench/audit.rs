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
//!   exactly.  Operation counts do not depend on the host, so a replay on
//!   any machine is an independent check of the figure.
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
    check_line_seal, child_main, json_roundtrip, ChildInput, Record, CHILD_INPUT_SCHEMA,
};
use crate::cryptanalysis::ecbench::runner::{read_record_lines, read_session, PlanDoc};
use crate::cryptanalysis::ecbench::spec::{plan, Spec};
use crate::cryptanalysis::ecbench::workload::CurveSpec;

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
    pub ok: bool,
    pub problems: Vec<String>,
    pub checks: BTreeMap<String, bool>,
    pub records: u64,
    pub verified_records: u64,
    pub replays: Vec<Replay>,
    /// SHA-256 of every file the audit read.
    pub files: BTreeMap<String, String>,
    pub auditor_binary_sha256: Option<String>,
    pub audited_unix_ms: u128,
}

fn file_sha(dir: &Path, name: &str, files: &mut BTreeMap<String, String>) -> Option<String> {
    let bytes = std::fs::read(dir.join(name)).ok()?;
    let h = sha256_hex(&bytes);
    files.insert(name.to_string(), h.clone());
    Some(h)
}

/// Audit `dir`; replay up to `replay` deterministic measured runs.
pub fn audit(dir: &Path, replay: usize) -> Result<AuditReport, String> {
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
    check(
        "records_hash",
        !complete || rec_sha == session.records_sha256,
        "records.jsonl does not hash to session.records_sha256".into(),
        &mut problems,
    );
    check(
        "session_complete",
        complete,
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
        p.workloads == doc.workloads,
        "the spec rebuilds different workloads".into(),
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
        p.arms == doc.arms,
        "the spec resolves to different methods".into(),
        &mut problems,
    );

    // Records.
    let lines = read_record_lines(dir)?;
    let records: Vec<Record> = lines.iter().map(|(_, r)| r.clone()).collect();
    check(
        "record_count",
        !complete
            || (records.len() as u64 == session.records_written
                && records.len() == p.executions.len()),
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
        if r.arm != arm.name || r.method != arm.method {
            record_problems.push(format!("{tag}: arm or method differs from the plan"));
        }
        if &r.workload != w {
            record_problems.push(format!("{tag}: workload differs from the plan"));
        }
        // Re-check the answer here.
        let inst = p.instance(ex.workload);
        let rec: Option<u64> = r.outcome.recovered.as_deref().and_then(|s| s.parse().ok());
        let target_ok = rec.map(|k| inst.mul_generator_hex(k).as_ref() == Some(&w.target));
        let planted_ok = rec.map(|k| k == w.planted);
        if target_ok != r.outcome.matches_target || planted_ok != r.outcome.matches_planted {
            record_problems.push(format!("{tag}: verification does not reproduce"));
        }
        let should_verify = target_ok == Some(true) && planted_ok == Some(true);
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
                    if Some(json_roundtrip(rep.total_gae).to_bits())
                        != r.cost.total_gae.map(f64::to_bits)
                    {
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
                    if rep
                        .phases
                        .iter()
                        .map(|p| json_roundtrip(p.gae).to_bits())
                        .ne(r.phases.iter().map(|p| p.gae.to_bits()))
                    {
                        diffs.push("phase gae");
                    }
                    if rep.counters != r.counters {
                        diffs.push("counters");
                    }
                    if rep.factor_base != r.factor_base {
                        diffs.push("factor_base");
                    }
                    (
                        diffs.is_empty(),
                        if diffs.is_empty() {
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
        session_id: session.session_id,
        ok: problems.is_empty(),
        problems,
        checks,
        records: records.len() as u64,
        verified_records: verified,
        replays,
        files,
        auditor_binary_sha256,
        audited_unix_ms: std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_millis())
            .unwrap_or(0),
    })
}
