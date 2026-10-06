//! Content-addressed shared-compute wrapper for the bounded P-256 walk.
//!
//! This module intentionally describes one deterministic prefix, not a
//! scalable shard.  Labels do not enter [`P256IsogenyTask::work_sha256`], so
//! schedulers cannot manufacture independent work by renaming the same walk.

use crate::cryptanalysis::p256_isogeny_campaign::P256_ICV1;
use crate::cryptanalysis::p256_isogeny_walk::{
    to_json_lines, verify, CycleEvent, WalkCertificate, DEFAULT_STEPS, PHI11_SHA256, WALK_DEGREE,
};
use crate::hash::sha256::sha256;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value as Json};

pub const TASK_SCHEMA: &str = "p256.isogeny-walk-task/v1";
pub const WORK_SCHEMA: &str = "p256.isogeny-walk-work/v1";
pub const RESULT_SCHEMA: &str = "p256.isogeny-walk-task-result/v1";
pub const CAIRN_RECEIPT_SCHEMA: &str = "p256.isogeny.cairn.walk-task-receipt/v1";
pub const DIRECTION_RULE: &str = "edge0=lower-j-root;later=exclude-predecessor";
pub const CERTIFICATE_FILE: &str = "certificate.jsonl.gz";
pub const TASK_FILE: &str = "task.json";
pub const RESULT_FILE: &str = "result.json";
pub const MAX_STEPS: u64 = DEFAULT_STEPS;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct P256IsogenyTask {
    pub schema: String,
    pub run_id: String,
    pub task_id: String,
    pub source_commit: String,
    pub curve: String,
    pub degree: u64,
    pub steps: u64,
    pub direction_rule: String,
    pub phi11_sha256: String,
}

#[derive(Serialize)]
struct WorkIdentity<'a> {
    schema: &'static str,
    source_commit: &'a str,
    curve: &'a str,
    degree: u64,
    steps: u64,
    direction_rule: &'a str,
    phi11_sha256: &'a str,
}

impl P256IsogenyTask {
    pub fn new(
        run_id: impl Into<String>,
        task_id: impl Into<String>,
        source_commit: impl Into<String>,
        steps: u64,
    ) -> Self {
        Self {
            schema: TASK_SCHEMA.into(),
            run_id: run_id.into(),
            task_id: task_id.into(),
            source_commit: source_commit.into(),
            curve: P256_ICV1.into(),
            degree: WALK_DEGREE,
            steps,
            direction_rule: DIRECTION_RULE.into(),
            phi11_sha256: PHI11_SHA256.into(),
        }
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.schema != TASK_SCHEMA {
            return Err(format!("schema must be {TASK_SCHEMA}"));
        }
        validate_segment("run_id", &self.run_id)?;
        validate_segment("task_id", &self.task_id)?;
        validate_commit(&self.source_commit)?;
        if self.curve != P256_ICV1
            || self.degree != WALK_DEGREE
            || self.direction_rule != DIRECTION_RULE
            || self.phi11_sha256 != PHI11_SHA256
        {
            return Err("task does not name the frozen P-256 degree-11 walk".into());
        }
        if !(1..=MAX_STEPS).contains(&self.steps) {
            return Err(format!("steps must be between 1 and {MAX_STEPS}"));
        }
        Ok(())
    }

    /// Hash of semantic work.  Scheduler labels are deliberately absent.
    pub fn work_sha256(&self) -> Result<String, String> {
        self.validate()?;
        json_sha256(&WorkIdentity {
            schema: WORK_SCHEMA,
            source_commit: &self.source_commit,
            curve: &self.curve,
            degree: self.degree,
            steps: self.steps,
            direction_rule: &self.direction_rule,
            phi11_sha256: &self.phi11_sha256,
        })
    }

    /// Hash of the full labelled task, useful for binding a result to its
    /// exact input file.
    pub fn task_sha256(&self) -> Result<String, String> {
        self.validate()?;
        json_sha256(self)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VerifiedSummary {
    pub emitted_steps: u64,
    pub verified_steps: u64,
    pub distinct_j: u64,
    pub first_cycle: Option<CycleEvent>,
    pub serialized_edge_bytes: u64,
    pub edge_stream_sha256: String,
}

impl VerifiedSummary {
    fn from_certificate(certificate: &WalkCertificate) -> Self {
        Self {
            emitted_steps: certificate.summary.emitted_steps,
            verified_steps: certificate.summary.verified_steps,
            distinct_j: certificate.summary.distinct_j,
            first_cycle: certificate.summary.first_cycle.clone(),
            serialized_edge_bytes: certificate.summary.serialized_edge_bytes,
            edge_stream_sha256: certificate.summary.edge_stream_sha256.clone(),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct P256IsogenyTaskResult {
    pub schema: String,
    pub task_sha256: String,
    pub work_sha256: String,
    pub certificate_sha256: String,
    pub certificate_bytes: u64,
    pub certificate_semantic_sha256: String,
    pub summary: VerifiedSummary,
}

/// Build a result only after a complete independent replay.  `artifact_bytes`
/// are the exact stored bytes (normally deterministic gzip), while the
/// semantic hash zeros generation timing before serializing the certificate.
pub fn build_result(
    task: &P256IsogenyTask,
    artifact_bytes: &[u8],
    certificate: &WalkCertificate,
) -> Result<P256IsogenyTaskResult, String> {
    task.validate()?;
    let mut checked = certificate.clone();
    verify(&mut checked)?;
    if checked.header.requested_steps != task.steps {
        return Err("certificate step count does not match the task".into());
    }
    let certificate_bytes = u64::try_from(artifact_bytes.len())
        .map_err(|_| "certificate length does not fit u64".to_string())?;
    Ok(P256IsogenyTaskResult {
        schema: RESULT_SCHEMA.into(),
        task_sha256: task.task_sha256()?,
        work_sha256: task.work_sha256()?,
        certificate_sha256: sha256_hex(artifact_bytes),
        certificate_bytes,
        certificate_semantic_sha256: semantic_certificate_sha256(&checked)?,
        summary: VerifiedSummary::from_certificate(&checked),
    })
}

/// Replay a result and require every stored field to equal a fresh derivation.
pub fn verify_result(
    task: &P256IsogenyTask,
    result: &P256IsogenyTaskResult,
    artifact_bytes: &[u8],
    certificate: &WalkCertificate,
) -> Result<(), String> {
    if result.schema != RESULT_SCHEMA {
        return Err(format!("result schema must be {RESULT_SCHEMA}"));
    }
    let expected = build_result(task, artifact_bytes, certificate)?;
    if result != &expected {
        return Err("result metadata does not match the replayed certificate".into());
    }
    Ok(())
}

/// Deterministic digest of the certificate with machine timing removed.
pub fn semantic_certificate_sha256(certificate: &WalkCertificate) -> Result<String, String> {
    let mut normalized = certificate.clone();
    normalized.summary.generation_ms = 0;
    normalized.summary.edges_per_second = 0.0;
    Ok(sha256_hex(&to_json_lines(&normalized)?))
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CairnWalkReceipt {
    pub schema: String,
    pub claim_class: String,
    pub verification: String,
    pub source_commit: String,
    pub curve: String,
    pub degree: u64,
    pub steps: u64,
    pub direction_rule: String,
    pub phi11_sha256: String,
    pub work_sha256: String,
    pub certificate_semantic_sha256: String,
    pub summary: VerifiedSummary,
}

impl CairnWalkReceipt {
    pub fn from_verified(
        task: &P256IsogenyTask,
        result: &P256IsogenyTaskResult,
        artifact_bytes: &[u8],
        certificate: &WalkCertificate,
    ) -> Result<Self, String> {
        verify_result(task, result, artifact_bytes, certificate)?;
        Ok(Self {
            schema: CAIRN_RECEIPT_SCHEMA.into(),
            claim_class: "coordination-receipt".into(),
            verification: "publisher-local-full-replay".into(),
            source_commit: task.source_commit.clone(),
            curve: task.curve.clone(),
            degree: task.degree,
            steps: task.steps,
            direction_rule: task.direction_rule.clone(),
            phi11_sha256: task.phi11_sha256.clone(),
            work_sha256: result.work_sha256.clone(),
            certificate_semantic_sha256: result.certificate_semantic_sha256.clone(),
            summary: result.summary.clone(),
        })
    }

    pub fn validate(
        &self,
        task: &P256IsogenyTask,
        result: &P256IsogenyTaskResult,
    ) -> Result<(), String> {
        task.validate()?;
        if result.schema != RESULT_SCHEMA
            || result.task_sha256 != task.task_sha256()?
            || result.work_sha256 != task.work_sha256()?
        {
            return Err("result identity does not match the task".into());
        }
        let expected = Self {
            schema: CAIRN_RECEIPT_SCHEMA.into(),
            claim_class: "coordination-receipt".into(),
            verification: "publisher-local-full-replay".into(),
            source_commit: task.source_commit.clone(),
            curve: task.curve.clone(),
            degree: task.degree,
            steps: task.steps,
            direction_rule: task.direction_rule.clone(),
            phi11_sha256: task.phi11_sha256.clone(),
            work_sha256: task.work_sha256()?,
            certificate_semantic_sha256: result.certificate_semantic_sha256.clone(),
            summary: result.summary.clone(),
        };
        if self != &expected {
            return Err("Cairn receipt does not match the semantic task result".into());
        }
        Ok(())
    }

    pub fn artifact(&self) -> Result<Json, String> {
        serde_json::to_value(self).map_err(|error| error.to_string())
    }

    pub fn claim_key(&self) -> String {
        format!("p256-isogeny-walk:{}", self.work_sha256)
    }
}

/// Emit, but never submit, a taskq command spec for this one work unit.
pub fn plan_taskq(
    task: &P256IsogenyTask,
    repo: &str,
    queue: &str,
    timeout_seconds: u64,
) -> Result<Json, String> {
    task.validate()?;
    validate_segment("repo", repo)?;
    validate_segment("queue", queue)?;
    if timeout_seconds == 0 {
        return Err("timeout_seconds must be positive".into());
    }
    let task_json = serde_json::to_string(task).map_err(|error| error.to_string())?;
    let work = task.work_sha256()?;
    Ok(json!({
        "schema": "taskq.task-spec/v1",
        "queue": queue,
        "kind": "command",
        "source": {
            "repo": repo,
            "commit": task.source_commit.clone(),
        },
        "command": {
            "argv": [
                "target/release/p256_isogeny_task",
                "run",
                "--task-json",
                task_json,
            ],
            "setup": [[
                "cargo",
                "build",
                "--locked",
                "--release",
                "--bin",
                "p256_isogeny_task",
            ]],
            "cwd": ".",
        },
        "limits": {
            "timeout_seconds": timeout_seconds,
            "setup_timeout_seconds": timeout_seconds.max(3600),
        },
        "retry": {"max_attempts": 3},
        "idempotency_key": format!("p256-isogeny-walk/{work}"),
        "labels": {
            "experiment": "p256-isogeny-walk",
            "run": task.run_id.clone(),
            "task": task.task_id.clone(),
            "work": &work[..16],
        },
    }))
}

fn json_sha256(value: &impl Serialize) -> Result<String, String> {
    let bytes = serde_json::to_vec(value).map_err(|error| error.to_string())?;
    Ok(sha256_hex(&bytes))
}

fn sha256_hex(bytes: &[u8]) -> String {
    hex::encode(sha256(bytes))
}

fn validate_commit(value: &str) -> Result<(), String> {
    if value.len() != 40
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err("source_commit must be a full lower-case 40-hex SHA".into());
    }
    Ok(())
}

fn validate_segment(name: &str, value: &str) -> Result<(), String> {
    let valid = !value.is_empty()
        && value.len() <= 128
        && value != "."
        && value != ".."
        && value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-'));
    if !valid {
        return Err(format!("{name} must be a safe nonempty identifier"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::p256_isogeny_walk::{generate, to_json_lines};

    const COMMIT: &str = "0123456789abcdef0123456789abcdef01234567";

    fn task(run: &str, id: &str, steps: u64) -> P256IsogenyTask {
        P256IsogenyTask::new(run, id, COMMIT, steps)
    }

    #[test]
    fn task_validation_is_fail_closed() {
        let valid = task("pilot", "prefix", 4);
        valid.validate().unwrap();

        let mut bad = valid.clone();
        bad.steps = MAX_STEPS + 1;
        assert!(bad.validate().is_err());
        let mut bad = valid.clone();
        bad.curve.push('0');
        assert!(bad.validate().is_err());
        let mut bad = valid.clone();
        bad.source_commit = "ABCDEF".repeat(6) + "ABCD";
        assert!(bad.validate().is_err());
        let mut bad = valid;
        bad.task_id = "../task".into();
        assert!(bad.validate().is_err());
    }

    #[test]
    fn labels_do_not_create_new_semantic_work() {
        let first = task("run-a", "task-a", 4);
        let second = task("run-b", "task-b", 4);
        assert_eq!(first.work_sha256().unwrap(), second.work_sha256().unwrap());
        assert_ne!(first.task_sha256().unwrap(), second.task_sha256().unwrap());

        let changed = task("run-a", "task-a", 5);
        assert_ne!(first.work_sha256().unwrap(), changed.work_sha256().unwrap());
    }

    #[test]
    fn taskq_plan_pins_source_and_semantic_idempotency() {
        let task = task("pilot", "prefix", 4);
        let plan = plan_taskq(&task, "crypto", "cpu", 1800).unwrap();
        assert_eq!(plan["schema"], "taskq.task-spec/v1");
        assert_eq!(plan["source"]["commit"], COMMIT);
        assert_eq!(
            plan["idempotency_key"],
            format!("p256-isogeny-walk/{}", task.work_sha256().unwrap())
        );
        assert_eq!(plan["retry"]["max_attempts"], 3);
    }

    #[test]
    fn short_result_replays_and_tampering_fails() {
        let pilot = task("pilot", "prefix", 2);
        let certificate = generate(pilot.steps).unwrap();
        let artifact = to_json_lines(&certificate).unwrap();
        let result = build_result(&pilot, &artifact, &certificate).unwrap();
        verify_result(&pilot, &result, &artifact, &certificate).unwrap();

        let receipt =
            CairnWalkReceipt::from_verified(&pilot, &result, &artifact, &certificate).unwrap();
        receipt.validate(&pilot, &result).unwrap();
        assert_eq!(
            receipt.claim_key(),
            format!("p256-isogeny-walk:{}", result.work_sha256)
        );

        let relabelled = task("another-run", "another-task", 2);
        let relabelled_result = build_result(&relabelled, &artifact, &certificate).unwrap();
        let relabelled_receipt = CairnWalkReceipt::from_verified(
            &relabelled,
            &relabelled_result,
            &artifact,
            &certificate,
        )
        .unwrap();
        assert_eq!(receipt, relabelled_receipt);

        let mut tampered = result.clone();
        let replacement = if tampered.certificate_sha256.starts_with('f') {
            "e"
        } else {
            "f"
        };
        tampered.certificate_sha256.replace_range(0..1, replacement);
        assert!(verify_result(&pilot, &tampered, &artifact, &certificate).is_err());

        let mut tampered = certificate.clone();
        tampered.header.record = "not-a-header".into();
        assert!(verify_result(&pilot, &result, &artifact, &tampered).is_err());

        let mut tampered = certificate;
        tampered.summary.bytes_per_edge += 1.0;
        assert!(verify_result(&pilot, &result, &artifact, &tampered).is_err());
    }

    #[test]
    fn unknown_task_fields_are_rejected() {
        let mut value = serde_json::to_value(task("pilot", "prefix", 1)).unwrap();
        value["offset"] = json!(1);
        assert!(serde_json::from_value::<P256IsogenyTask>(value).is_err());
    }
}
