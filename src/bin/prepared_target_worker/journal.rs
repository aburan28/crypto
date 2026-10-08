//! Durable attempt files; inspection proves structure, never native admission.
use super::native;
use crypto_lib::cryptanalysis::{
    koblitz_index_calculus::QueryAttempt,
    prepared_n17_target::F5TargetObserver,
    prepared_sat_control::{canonical_sha, sha256},
};
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
};

pub const QUESTION: &str = "disclosed-n17-target-attempt-files-v1";
pub fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Family {
    MatrixF5,
    Cryptominisat,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Binding {
    pub schema_version: u32,
    pub question: String,
    pub registration_sha256: String,
    pub worker_sha256: String,
    pub family: Family,
    pub target: [u64; 2],
    pub algorithm_seed: u64,
    pub max_queries: usize,
}
impl Binding {
    pub fn validate(&self) -> Result<(), String> {
        native::require(
            self.schema_version == 1
                && self.question == QUESTION
                && digest(&self.registration_sha256)
                && digest(&self.worker_sha256)
                && self.target.iter().all(|&v| v < 1 << 17)
                && (1..=8).contains(&self.max_queries),
            "outside n17 target journal contract",
        )
    }
    fn coefficients(&self) -> Vec<[u64; 2]> {
        let mut rng = StdRng::seed_from_u64(self.algorithm_seed ^ 0x44455343454e5400);
        (0..self.max_queries)
            .map(|_| [rng.gen_range(1..65587), rng.gen_range(1..65587)])
            .collect()
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Header {
    pub trial: usize,
    pub a: u64,
    pub b: u64,
    pub target: [u64; 2],
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Envelope {
    schema_version: u32,
    binding_sha256: String,
    header: Header,
    start_sha256: Option<String>,
    body: Value,
}
fn keys(value: &Value, expected: &[&str]) -> Result<(), String> {
    native::require(
        value.as_object().is_some_and(|v| {
            v.len() == expected.len() && expected.iter().all(|k| v.contains_key(*k))
        }),
        "target record fields differ",
    )
}
fn header(
    binding: &Binding,
    coefficients: &[[u64; 2]],
    body: &Value,
    completion: bool,
) -> Result<Header, String> {
    let row = if binding.family == Family::MatrixF5 && completion {
        &body["attempt"]
    } else {
        body
    };
    let h = Header {
        trial: row["trial"]
            .as_u64()
            .and_then(|v| v.try_into().ok())
            .ok_or("invalid target trial")?,
        a: row["a"].as_u64().ok_or("invalid target coefficient a")?,
        b: row["b"].as_u64().ok_or("invalid target coefficient b")?,
        target: binding.target,
    };
    native::require(
        h.trial < binding.max_queries && [h.a, h.b] == coefficients[h.trial],
        "target record differs from frozen seeded input law",
    )?;
    if binding.family == Family::MatrixF5 {
        if completion {
            keys(body, &["target", "attempt", "record_kind"])?;
            native::require(
                body["target"] == json!(binding.target)
                    && body["record_kind"] == "completed-pdp-before-scalar-recovery",
                "F5 completion target/kind differs",
            )?;
            let _: QueryAttempt =
                serde_json::from_value(body["attempt"].clone()).map_err(|e| e.to_string())?;
        } else {
            keys(
                body,
                &[
                    "trial",
                    "a",
                    "b",
                    "target",
                    "query_rule",
                    "decomposition_started",
                ],
            )?;
            native::require(
                body["target"] == json!(binding.target)
                    && body["query_rule"] == "seeded-sample-aG-plus-bQ"
                    && body["decomposition_started"] == false,
                "F5 start target/kind differs",
            )?;
        }
    } else if !completion {
        keys(
            body,
            &[
                "trial",
                "a",
                "b",
                "public_query",
                "outcome",
                "source_model_valid",
                "witness_indices",
                "candidate_scalar",
                "backend_called",
            ],
        )?;
        native::require(
            matches!(
                body["outcome"].as_str(),
                Some("query_planned" | "identity_query")
            ) && body["source_model_valid"] == false
                && body["backend_called"] == false
                && body["witness_indices"].is_null()
                && body["candidate_scalar"].is_null(),
            "SAT start already claims computation or recovery",
        )?;
    } else {
        native::require(
            body["outcome"].is_string() && body["backend_called"].is_boolean(),
            "SAT completion lacks disposition",
        )?;
    }
    Ok(h)
}
fn write(root: &Path, name: &str, record: &Envelope) -> Result<(), String> {
    native::require(
        serde_json::to_vec(record).map_err(|e| e.to_string())?.len() <= 8 * 1024 * 1024,
        "target record exceeds bound",
    )?;
    native::save(&root.join(name), &json!(record))
}
fn read(root: &Path, name: &str) -> Result<(Envelope, String), String> {
    let bytes = native::read(&root.join(name), 8 * 1024 * 1024)?;
    let record = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    native::require(
        native::read(&root.join(name), 8 * 1024 * 1024)? == bytes,
        "target record changed during read",
    )?;
    Ok((record, sha256(&bytes)))
}
pub struct Journal {
    root: PathBuf,
    binding: Binding,
    binding_sha256: String,
    coefficients: Vec<[u64; 2]>,
    next: usize,
    pending: Option<Envelope>,
    faulted: bool,
}
impl Journal {
    pub fn create(execution: &Path, binding: Binding) -> Result<Self, String> {
        binding.validate()?;
        let root = execution.join("attempts");
        fs::create_dir(&root).map_err(|e| e.to_string())?;
        #[cfg(unix)]
        fs::File::open(execution)
            .and_then(|parent| parent.sync_all())
            .map_err(|e| e.to_string())?;
        let root = root.canonicalize().map_err(|e| e.to_string())?;
        native::save(&root.join("binding.json"), &json!(binding))?;
        let binding_sha256 = canonical_sha(&json!(binding))?;
        let coefficients = binding.coefficients();
        Ok(Self {
            root,
            binding,
            binding_sha256,
            coefficients,
            next: 0,
            pending: None,
            faulted: false,
        })
    }
    pub fn start(&mut self, body: &Value) -> Result<(), String> {
        let result = self.start_inner(body);
        if result.is_err() {
            self.faulted = true;
        }
        result
    }
    fn start_inner(&mut self, body: &Value) -> Result<(), String> {
        native::require(
            !self.faulted && self.pending.is_none(),
            "target writer stopped or has pending query",
        )?;
        let h = header(&self.binding, &self.coefficients, body, false)?;
        native::require(h.trial == self.next, "target start chronology differs")?;
        let record = Envelope {
            schema_version: 1,
            binding_sha256: self.binding_sha256.clone(),
            header: h,
            start_sha256: None,
            body: body.clone(),
        };
        write(&self.root, &format!("start-{:03}.json", self.next), &record)?;
        self.pending = Some(record);
        Ok(())
    }
    pub fn complete(&mut self, body: &Value) -> Result<(), String> {
        let result = self.complete_inner(body);
        if result.is_err() {
            self.faulted = true;
        }
        result
    }
    fn complete_inner(&mut self, body: &Value) -> Result<(), String> {
        native::require(!self.faulted, "target writer stopped after failure")?;
        let start = self
            .pending
            .as_ref()
            .ok_or("target completion lacks start")?;
        let h = header(&self.binding, &self.coefficients, body, true)?;
        native::require(h == start.header, "target completion differs from start")?;
        let record = Envelope {
            schema_version: 1,
            binding_sha256: self.binding_sha256.clone(),
            header: h,
            start_sha256: Some(canonical_sha(&json!(start))?),
            body: body.clone(),
        };
        write(
            &self.root,
            &format!("completed-{:03}.json", self.next),
            &record,
        )?;
        self.pending = None;
        self.next += 1;
        Ok(())
    }
}
impl F5TargetObserver for Journal {
    fn started(&mut self, row: &Value) -> Result<(), String> {
        self.start(row)
    }
    fn completed(&mut self, row: &Value) -> Result<(), String> {
        self.complete(row)
    }
}

pub fn inspect(execution: &Path, expected_registration: &str) -> Result<Value, String> {
    native::require(
        digest(expected_registration),
        "external registration digest required",
    )?;
    let root = execution.join("attempts");
    let binding: Binding =
        serde_json::from_slice(&native::read(&root.join("binding.json"), 65536)?)
            .map_err(|e| e.to_string())?;
    binding.validate()?;
    native::require(
        binding.registration_sha256 == expected_registration,
        "target prefix registration differs",
    )?;
    let bound = canonical_sha(&json!(binding))?;
    let coefficients = binding.coefficients();
    let mut names = BTreeSet::new();
    for entry in fs::read_dir(&root).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        native::require(
            entry.file_type().map_err(|e| e.to_string())?.is_file(),
            "target journal member is not regular",
        )?;
        names.insert(
            entry
                .file_name()
                .into_string()
                .map_err(|_| "non-UTF8 target journal member")?,
        );
        native::require(
            names.len() <= 1 + 2 * binding.max_queries,
            "too many target journal files",
        )?;
    }
    let mut expected = BTreeSet::from(["binding.json".to_string()]);
    let mut complete = Vec::new();
    let mut pending = None;
    for trial in 0..binding.max_queries {
        let name = format!("start-{trial:03}.json");
        if !names.contains(&name) {
            break;
        }
        expected.insert(name.clone());
        let (start, start_bytes_sha256) = read(&root, &name)?;
        native::require(
            start.schema_version == 1
                && start.binding_sha256 == bound
                && start.start_sha256.is_none()
                && start.header.trial == trial
                && header(&binding, &coefficients, &start.body, false)? == start.header,
            "retained target start differs",
        )?;
        let completion_name = format!("completed-{trial:03}.json");
        if !names.contains(&completion_name) {
            pending = Some(json!(start.header));
            break;
        }
        expected.insert(completion_name.clone());
        let (completion, completion_bytes_sha256) = read(&root, &completion_name)?;
        native::require(
            completion.schema_version == 1
                && completion.binding_sha256 == bound
                && completion.header == start.header
                && header(&binding, &coefficients, &completion.body, true)? == start.header
                && completion.start_sha256.as_deref()
                    == Some(canonical_sha(&json!(start))?.as_str()),
            "retained target completion differs",
        )?;
        complete.push(json!({"header":completion.header,
            "start_sha256":start_bytes_sha256,"completed_sha256":completion_bytes_sha256}));
    }
    native::require(
        names == expected,
        "target prefix contains gaps, orphan completions or unknown files",
    )?;
    Ok(
        json!({"schema_version":1,"status":"RETAINED_TARGET_PREFIX_STRUCTURE",
        "binding":binding,"completed_attempts":complete,"pending_start":pending,
        "searches_executed_by_inspection":0,"source_bound_execution_admitted":false,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,
        "full_goal_complete":false,"online_speedup":null}),
    )
}

#[cfg(test)]
#[path = "journal_tests.rs"]
mod tests;
