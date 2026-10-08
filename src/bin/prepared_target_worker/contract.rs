//! New one-target F5 capsule; no old registration can authorize this worker.
use super::{
    journal::{self, Binding, Family},
    native, preparation,
};
use crypto_lib::cryptanalysis::{
    prepared_n17_target::F5TargetPlan,
    prepared_sat_control::{canonical_sha, sha256},
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};
pub const SCOPE: &str = "disclosed-n17-native-one-target-f5-v1";
pub const WORKER: &str = "prepared_target_worker";
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub schema_version: u32,
    pub question: String,
    pub target: [u64; 2],
    pub plan: F5TargetPlan,
    pub worker_timeout_ms: u64,
}
impl Config {
    pub fn validate(&self) -> Result<(), String> {
        self.plan.validate()?;
        native::require(
            self.schema_version == 1
                && self.question == SCOPE
                && self.target.iter().all(|&v| v < 1 << 17)
                && (30_000..=900_000).contains(&self.worker_timeout_ms),
            "outside bounded n17 target worker configuration",
        )
    }
    pub fn binding(&self, hash: &str, worker: &str) -> Binding {
        Binding {
            schema_version: 1,
            question: journal::QUESTION.into(),
            registration_sha256: hash.into(),
            worker_sha256: worker.into(),
            family: Family::MatrixF5,
            target: self.target,
            algorithm_seed: self.plan.algorithm_seed,
            max_queries: self.plan.max_queries,
        }
    }
}
pub fn config(capsule: &Path) -> Result<Config, String> {
    let cfg: Config = serde_json::from_slice(&native::read(&capsule.join("config.json"), 65536)?)
        .map_err(|e| e.to_string())?;
    cfg.validate()?;
    Ok(cfg)
}
pub fn environment() -> BTreeMap<String, String> {
    [
        ("LC_ALL".into(), "C".into()),
        ("RAYON_NUM_THREADS".into(), "1".into()),
    ]
    .into()
}
pub fn identity() -> Value {
    json!({"schema_version":1,"source_manifest_sha256":option_env!("IC_TARGET_SOURCE_MANIFEST_SHA256"),
        "build_sha256":option_env!("IC_TARGET_BUILD_SHA256"),
        "target_os":std::env::consts::OS,"target_arch":std::env::consts::ARCH})
}
pub fn require_identity(compiled: &Value, expected: &Value) -> Result<(), String> {
    native::require(
        compiled == expected
            && compiled["schema_version"] == 1
            && compiled["source_manifest_sha256"]
                .as_str()
                .is_some_and(journal::digest)
            && compiled["build_sha256"]
                .as_str()
                .is_some_and(journal::digest)
            && compiled["target_os"] == "macos"
            && compiled["target_arch"] == "aarch64",
        "target worker lacks matching compiled source/build identity",
    )
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreparationBinding {
    pub capsule: PathBuf,
    pub execution: PathBuf,
    pub registration_sha256: String,
    pub auditor_sha256: String,
    pub producer_sha256: String,
    pub mathematics_sha256: String,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Registration {
    pub schema_version: u32,
    pub scope: String,
    pub source_commit: String,
    pub hardware: Value,
    pub immutable_files: Value,
    pub source_manifest_sha256: String,
    pub config_sha256: String,
    pub host_context_sha256: String,
    pub worker_sha256: String,
    pub auditor_sha256: String,
    pub worker_build_identity: Value,
    pub runtime_environment: BTreeMap<String, String>,
    pub preparation: PreparationBinding,
    pub validation_only: bool,
}
fn pins(record: &Registration) -> Result<(), String> {
    let prep = &record.preparation;
    native::require(
        record.schema_version == 1
            && record.scope == SCOPE
            && record.hardware == json!({"os":"macos","architecture":"aarch64"})
            && record.source_commit.len() == 40
            && record
                .source_commit
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            && [
                &record.source_manifest_sha256,
                &record.config_sha256,
                &record.host_context_sha256,
                &record.worker_sha256,
                &record.auditor_sha256,
                &prep.registration_sha256,
                &prep.auditor_sha256,
                &prep.producer_sha256,
                &prep.mathematics_sha256,
            ]
            .into_iter()
            .all(|s| journal::digest(s))
            && record.runtime_environment == environment()
            && prep.capsule.is_absolute()
            && prep.execution.is_absolute(),
        "target registration scope/source/preparation contract differs",
    )?;
    require_identity(&record.worker_build_identity, &record.worker_build_identity)?;
    native::require(
        record.worker_build_identity["source_manifest_sha256"] == record.source_manifest_sha256,
        "target compiled source pin differs",
    )
}
/// Validate the sealed registration and its small sidecars as data. This is
/// also used for portable build-custody replay, where the immutable tree is
/// checked inside a bounded archive instead of extracted or executed.
pub fn registration(capsule: &Path, expected: &str) -> Result<Registration, String> {
    native::require(
        journal::digest(expected),
        "external target registration seal required",
    )?;
    let bytes = native::read(&capsule.join("registration.json"), 16 * 1024 * 1024)?;
    let record: Registration = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    pins(&record)?;
    let hash = canonical_sha(&json!(record))?;
    native::require(
        hash == expected
            && native::load(&capsule.join("seal.json"))?["registration_sha256"] == hash,
        "target registration differs from external seal",
    )?;
    native::require(
        record.config_sha256 == sha256(&native::read(&capsule.join("config.json"), 65536)?)
            && record.host_context_sha256
                == sha256(&native::read(&capsule.join("host-context.json"), 65536)?),
        "target configuration changed",
    )?;
    config(capsule)?;
    Ok(record)
}
pub fn check_capsule(capsule: &Path, expected: &str) -> Result<Registration, String> {
    let record = registration(capsule, expected)?;
    native::check_tree(&capsule.join("immutable"), &record.immutable_files)?;
    native::require(
        canonical_sha(&native::inventory(&capsule.join("immutable/source"))?)?
            == record.source_manifest_sha256
            && sha256(&native::read(
                &capsule.join("immutable/bin").join(WORKER),
                128 * 1024 * 1024,
            )?) == record.worker_sha256
            && sha256(&native::read(
                &capsule.join("immutable/bin/icprog"),
                128 * 1024 * 1024,
            )?) == record.auditor_sha256,
        "target source or executable binding differs",
    )?;
    Ok(record)
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Claim {
    pub schema_version: u32,
    pub scope: String,
    pub registration_sha256: String,
    pub execution: PathBuf,
    pub worker_sha256: String,
}
pub fn claim(
    capsule: &Path,
    execution: &Path,
    record: &Registration,
    expected: &str,
) -> Result<(), String> {
    let c: Claim = serde_json::from_slice(&native::read(&capsule.join("consumed.json"), 65536)?)
        .map_err(|e| e.to_string())?;
    native::require(
        !record.validation_only
            && c.schema_version == 1
            && c.scope == SCOPE
            && c.registration_sha256 == expected
            && c.worker_sha256 == record.worker_sha256
            && execution.is_absolute()
            && c.execution == execution,
        "target claim differs, is missing, or registration is validation-only",
    )
}
pub fn preparation_paths(
    binding: &PreparationBinding,
) -> Result<(PathBuf, PathBuf, PathBuf), String> {
    let capsule = binding.capsule.canonicalize().map_err(|e| e.to_string())?;
    let execution = binding
        .execution
        .canonicalize()
        .map_err(|e| e.to_string())?;
    native::require(
        capsule == binding.capsule && execution == binding.execution,
        "preparation paths are not canonical original locations",
    )?;
    let (r, hash) = preparation::check_capsule(&capsule)?;
    preparation::claim(&capsule, &execution, &r, &hash)?;
    native::require(
        hash == binding.registration_sha256
            && r.auditor_sha256 == binding.auditor_sha256
            && preparation::config(&capsule)?.plan.family
                == crypto_lib::cryptanalysis::prepared_ordinary::Family::MatrixF5,
        "target preparation family/registration/auditor differs",
    )?;
    let checker = capsule.join("immutable/bin/icprog");
    native::require(
        sha256(&native::read(&checker, 128 * 1024 * 1024)?) == binding.auditor_sha256,
        "original preparation auditor changed",
    )?;
    Ok((capsule, execution, checker))
}
pub fn preparation_receipt(
    binding: &PreparationBinding,
    receipt: &Value,
    producer: &Value,
    producer_bytes: &[u8],
) -> Result<Value, String> {
    let math = &producer["mathematical_input"];
    native::require(
        receipt["schema_version"] == 1
            && receipt["status"] == "PASS_NATIVE_SOURCE_BOUND_ORDINARY_PREPARATION_AUDIT"
            && receipt["registration_sha256"] == binding.registration_sha256
            && receipt["checker_sha256"] == binding.auditor_sha256
            && receipt["source_bound_execution_admitted"] == true
            && receipt["mathematical_preparation_complete"] == true
            && receipt["producer_sha256"] == binding.producer_sha256
            && sha256(producer_bytes) == binding.producer_sha256
            && receipt["mathematics"]["input_sha256"] == binding.mathematics_sha256
            && canonical_sha(math)? == binding.mathematics_sha256
            && receipt["mathematics"]["declared_family"] == "matrix_f5"
            && receipt["mathematics"]["rank"] == 29
            && receipt["mathematics"]["usable_points"] == 62
            && receipt["mathematics"]["panel_complete"] == true,
        "preparation recheck lacks bound native execution, complete rank or exact input",
    )?;
    Ok(math.clone())
}

#[cfg(test)]
#[path = "contract_tests.rs"]
mod tests;
