//! Non-arithmetic supervisor for the frozen Role-1 P-192 weighted-CM gate.
//!
//! This binary performs process sequencing, byte hashing, source-dependency
//! auditing, and receipt construction.  It deliberately contains no curve,
//! ideal, factorization, or certificate arithmetic.

use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::ffi::OsString;
use std::fs::{self, OpenOptions};
use std::io::{Read, Seek, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

type Result<T> = std::result::Result<T, String>;

const EXPERIMENT_ID: &str = "EXP-SCURVE-1a8daf";
const PROTOCOL_VERSION: u64 = 2;
const NATIVE_REPOSITORY: &str = "https://github.com/aburan28/crypto";
const PROTOCOL_REPOSITORY: &str = "https://github.com/aburan28/crypto-autoresearcher";
const ROLE1_ADDENDUM_PATH: &str =
    "experiments/EXP-SCURVE-1a8daf/amendments/v2_addendum_role1_interface.yaml";
const ROLE1_INTERFACE_MANIFEST_PATH: &str =
    "research/p192-weighted-cm-20261007-v2-role1-interface/SHA256SUMS";
const EMBEDDED_PROTOCOL_COMMIT: &str = env!("P192_WCM_EMBEDDED_PROTOCOL_COMMIT");
const EMBEDDED_SUPERVISOR_COMMIT: &str = env!("P192_WCM_EMBEDDED_SUPERVISOR_COMMIT");
const EMBEDDED_VERIFIER_COMMIT: &str = env!("P192_WCM_EMBEDDED_VERIFIER_COMMIT");
const EMBEDDED_BUILD_PROFILE: &str = env!("P192_WCM_EMBEDDED_BUILD_PROFILE");
const EMBEDDED_SUPERVISOR_DIRTY: &str = env!("P192_WCM_EMBEDDED_SUPERVISOR_DIRTY");
const EMBEDDED_GIT_METADATA_PRESENT: &str = env!("P192_WCM_EMBEDDED_GIT_METADATA_PRESENT");
const PRODUCER_PATHS: [&str; 10] = [
    "identity.json",
    "maximal-order.json",
    "factor-base.json",
    "factor-base.sha256",
    "controls.json",
    "reference-box.json",
    "reference-box/candidate-records.bin",
    "reference-box/disposition.bin",
    "reference-box/certificates.bin",
    "verification.json",
];
const COMPARISON_FIELDS: [&str; 8] = [
    "curve_tuple_sha256",
    "cm_order_tuple_sha256",
    "factor_base_sha256",
    "control_results_sha256",
    "reference_candidate_shard_sha256",
    "reference_candidate_shell_sha256",
    "reference_disposition_shard_sha256",
    "reference_disposition_shell_sha256",
];
const PRODUCER_CHECKS: [&str; 11] = [
    "identity_pass",
    "maximal_order_pass",
    "factor_base_pass",
    "controls_pass",
    "ref0_count_1025",
    "ref0_candidate_encoding",
    "ref0_candidate_root",
    "ref0_disposition_encoding",
    "ref0_disposition_root",
    "ref0_certificate_store",
    "source_commit_consistent",
];
const VERIFIER_CHECKS: [&str; 15] = [
    "source_inventory_complete",
    "source_hashes_stable",
    "identity_rederived",
    "maximal_order_rederived",
    "primality_proofs_replayed",
    "factor_base_regenerated",
    "controls_replayed",
    "ref0_order_regenerated",
    "ref0_candidate_bytes_equal",
    "ref0_candidate_roots_equal",
    "ref0_dispositions_regenerated",
    "ref0_disposition_bytes_equal",
    "ref0_disposition_roots_equal",
    "ref0_certificates_replayed",
    "source_commit_consistent",
];
const RUN_MANIFEST_PATHS: [&str; 19] = [
    "identity.json",
    "maximal-order.json",
    "factor-base.json",
    "factor-base.sha256",
    "controls.json",
    "reference-box.json",
    "verification.json",
    "independent-verification.json",
    "independent-verifier-receipt.json",
    "reference-box/candidate-records.bin",
    "reference-box/disposition.bin",
    "reference-box/certificates.bin",
    "independent-agreement.json",
    "dependency-audit.json",
    "command.txt",
    "environment.json",
    "stdout.log",
    "stderr.log",
    "raw-result.json",
];
const STDOUT_LOG_DOMAIN: &[u8] = b"P192-WCM-SUPERVISOR-STDOUT-v1\0";
const STDERR_LOG_DOMAIN: &[u8] = b"P192-WCM-SUPERVISOR-STDERR-v1\0";
const FORBIDDEN_VERIFIER_SOURCE_MARKERS: [&str; 6] = [
    "extern crate crypto_lib",
    "use crypto_lib",
    "crypto_lib::",
    "include!",
    "#[path",
    "../producer/",
];
const PROTECTED_PROTOCOL_PATHS: [&str; 11] = [
    "experiments/EXP-SCURVE-1a8daf/amendments/protocol-amendment.v2.yaml",
    "experiments/EXP-SCURVE-1a8daf/amendments/specification.v2.yaml",
    "experiments/EXP-SCURVE-1a8daf/amendments/v2_addendum_role1_interface.yaml",
    "research/p192-weighted-cm-20261007-v2/run-family.yaml",
    "research/p192-weighted-cm-20261007-v2-role1-interface/README.md",
    "research/p192-weighted-cm-20261007-v2-role1-interface/SHA256SUMS",
    "research/p192-weighted-cm-20261007-v2-role1-interface/check_role1_interface.py",
    "research/p192-weighted-cm-20261007-v2-role1-interface/run-family.role1-interface.yaml",
    "research/p192-weighted-cm-20261007-v2-role1-interface/schemas/role1-artifacts.schema.json",
    "research/p192-weighted-cm-20261007-v2-role1-interface/schemas/role1-dispatch-decision.schema.json",
    "tests/test_p192_wcm_role1_interface.py",
];
const UNAUTHORIZED_STAGE_LABELS: [&str; 7] = [
    "P192-WCM-BOX0",
    "P192-WCM-BOX1-SHELL",
    "P192-WCM-BOX2-SHELL",
    "P192-WCM-RELATIONS",
    "P192-WCM-MAPS",
    "P192-WCM-HOLDOUT",
    "P192-WCM-REPLICATION",
];
const CUSTODY_OUTPUTS: [&str; 6] = [
    "manifest.yaml",
    "command.txt",
    "environment.json",
    "stdout.log",
    "stderr.log",
    "raw-result.json",
];

#[derive(Clone, Debug)]
struct Arguments {
    producer: PathBuf,
    verifier: PathBuf,
    producer_source: PathBuf,
    verifier_source: PathBuf,
    protocol_repository: PathBuf,
    decision: PathBuf,
    run_dir: PathBuf,
    protocol_commit: String,
    source_commit: String,
    verifier_commit: String,
    audit_out: PathBuf,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Artifact {
    path: String,
    byte_length: u64,
    sha256: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct DispatchDocument {
    coordinator_decision: CoordinatorDecision,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct CoordinatorDecision {
    id: String,
    context: String,
    decision: String,
    target_ids: Vec<String>,
    rationale: Vec<String>,
    evidence_refs: Vec<String>,
    limitations: Vec<String>,
    next_actions: Vec<String>,
    authorization: Role1Authorization,
    knowledge_promotion: KnowledgePromotion,
    decided_by: String,
    decided_at: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct KnowledgePromotion {
    promoted: Vec<String>,
    not_warranted: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Role1Authorization {
    schema: String,
    execution_authorized: bool,
    scientific_execution_authorized: bool,
    scientific_result_recording_authorized: bool,
    experiment_id: String,
    protocol_version: u64,
    decision_path: String,
    stage: StageAuthorization,
    run: RunAuthorization,
    protocol: ProtocolAuthorization,
    executables: ExecutableAuthorization,
    dependency_audit: DependencyAuditAuthorization,
    supervisor_argv: Vec<String>,
    producer_argv: Vec<String>,
    verifier_argv: Vec<String>,
    launch: LaunchAuthorization,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct StageAuthorization {
    label: String,
    ordinal: u64,
    role: String,
    composite_invocations: u64,
    retry_authorized: bool,
    resume_authorized: bool,
    advance_to_later_role_authorized: bool,
    unauthorized_stage_labels: Vec<String>,
    custody_outputs: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct RunAuthorization {
    run_id: String,
    run_dir: String,
    must_be_fresh: bool,
    existing_path_refused: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ProtocolAuthorization {
    repository: String,
    repository_path: String,
    commit: String,
    must_be_ancestor_of_decision_commit: bool,
    protected_paths_unchanged_through_decision_commit: bool,
    protected_paths: Vec<String>,
    addendum: BoundProtocolFile,
    interface_manifest: BoundProtocolFile,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct BoundProtocolFile {
    path: String,
    byte_length: u64,
    sha256: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ExecutableAuthorization {
    supervisor: SupervisorAuthorization,
    producer: ProducerAuthorization,
    verifier: VerifierAuthorization,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct SupervisorAuthorization {
    repository: String,
    source_commit: String,
    repository_path: String,
    package_path: String,
    binary_path: String,
    byte_length: u64,
    binary_sha256: String,
    build_profile: String,
    worktree_clean: bool,
    contains_curve_arithmetic: bool,
    embedded_protocol_commit: String,
    embedded_supervisor_commit: String,
    embedded_verifier_commit: String,
    embedded_build_profile: String,
    embedded_supervisor_dirty: bool,
    embedded_git_metadata_present: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ProducerAuthorization {
    repository: String,
    source_commit: String,
    embedded_protocol_commit: String,
    embedded_source_commit: String,
    repository_path: String,
    package_path: String,
    binary_path: String,
    byte_length: u64,
    binary_sha256: String,
    build_profile: String,
    worktree_clean: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct VerifierAuthorization {
    repository: String,
    verifier_commit: String,
    embedded_verifier_commit: String,
    repository_path: String,
    package_path: String,
    binary_path: String,
    byte_length: u64,
    binary_sha256: String,
    build_profile: String,
    worktree_clean: bool,
    crypto_lib_found: bool,
    shared_implementation_components: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct DependencyAuditAuthorization {
    path: String,
    schema: String,
    must_be_fresh: bool,
    receipt_hash_field: String,
    hash_scope: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct LaunchAuthorization {
    no_path_search_for_any_executable: bool,
    runtime_commit_override_forbidden: bool,
    producer_before_verifier: bool,
    close_and_rehash_producer_outputs_before_verifier: bool,
    supervisor_writes_only_run_directory_artifacts: bool,
    verify_protocol_checkout_before_run_directory: bool,
    derive_and_retain_decision_binding: bool,
}

#[derive(Clone, Debug)]
struct DecisionRuntimeBinding {
    supervisor_binary: PathBuf,
    supervisor_binary_byte_length: u64,
    supervisor_binary_sha256: String,
    supervisor_argv: Vec<String>,
    producer_binary_byte_length: u64,
    producer_binary_sha256: String,
    producer_argv: Vec<String>,
    verifier_binary_byte_length: u64,
    verifier_binary_sha256: String,
    verifier_argv: Vec<String>,
}

#[derive(Debug)]
struct RetainedExecutable {
    path: PathBuf,
    source_file: fs::File,
    sealed_image: fs::File,
    byte_length: u64,
    sha256: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ProtocolBindingPhase {
    PreRun,
    PostRun,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
struct SourceHashRecord {
    path: String,
    sha256: String,
}

#[derive(Clone, Debug)]
struct TokenizedSource {
    path: String,
    bytes: Vec<u8>,
    tokens: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
struct SourceCloneMatch {
    producer_path: String,
    producer_token_start: usize,
    producer_token_end_exclusive: usize,
    verifier_path: String,
    verifier_token_start: usize,
    verifier_token_end_exclusive: usize,
    token_count: usize,
    normalized_tokens_sha256: String,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct VerificationEnvelope {
    schema: String,
    experiment_id: String,
    protocol_version: u64,
    protocol_commit: String,
    source_commit: String,
    #[serde(default)]
    verifier_commit: Option<String>,
    role: String,
    #[serde(default)]
    producer_binary_sha256: Option<String>,
    #[serde(default)]
    verifier_binary_sha256: Option<String>,
    #[serde(default)]
    artifacts: Option<Vec<Artifact>>,
    #[serde(default)]
    source_artifacts: Option<Vec<Artifact>>,
    bindings: BTreeMap<String, String>,
    checks: Vec<Value>,
    overall_status: String,
    #[serde(default)]
    independence: Option<Value>,
}

fn usage() -> &'static str {
    "p192_wcm_preflight_supervisor --producer ABS_BIN --verifier ABS_BIN --producer-source ABS_DIR --verifier-source ABS_DIR --protocol-repository ABS_DIR --decision ABS_FILE --run-dir ABS_DIR --protocol-commit HEX40 --source-commit HEX40 --verifier-commit HEX40 --audit-out ABS_FILE"
}

fn parse_args<I, S>(values: I) -> Result<Arguments>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let values = values
        .into_iter()
        .map(|value| {
            value
                .into()
                .into_string()
                .map_err(|_| "arguments must be UTF-8".to_owned())
        })
        .collect::<Result<Vec<_>>>()?;
    let expected = [
        "--producer",
        "--verifier",
        "--producer-source",
        "--verifier-source",
        "--protocol-repository",
        "--decision",
        "--run-dir",
        "--protocol-commit",
        "--source-commit",
        "--verifier-commit",
        "--audit-out",
    ];
    if values.len() != expected.len() * 2 {
        return Err(usage().to_owned());
    }
    for (index, option) in expected.iter().enumerate() {
        if values[index * 2] != *option {
            return Err(format!(
                "expected {option} at argv position {}\n{}",
                index * 2 + 1,
                usage()
            ));
        }
    }
    let arguments = Arguments {
        producer: PathBuf::from(&values[1]),
        verifier: PathBuf::from(&values[3]),
        producer_source: PathBuf::from(&values[5]),
        verifier_source: PathBuf::from(&values[7]),
        protocol_repository: PathBuf::from(&values[9]),
        decision: PathBuf::from(&values[11]),
        run_dir: PathBuf::from(&values[13]),
        protocol_commit: values[15].clone(),
        source_commit: values[17].clone(),
        verifier_commit: values[19].clone(),
        audit_out: PathBuf::from(&values[21]),
    };
    for (label, commit) in [
        ("protocol", &arguments.protocol_commit),
        ("source", &arguments.source_commit),
        ("verifier", &arguments.verifier_commit),
    ] {
        if !lower_hex(commit, 40) {
            return Err(format!(
                "{label} commit must be exactly 40 lowercase hexadecimal characters"
            ));
        }
    }
    Ok(arguments)
}

fn lower_hex(value: &str, length: usize) -> bool {
    lower_hex_digits(value, length) && value.bytes().any(|byte| byte != b'0')
}

fn lower_hex_digits(value: &str, length: usize) -> bool {
    value.len() == length
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn validate_embedded_provenance(arguments: &Arguments) -> Result<()> {
    if EMBEDDED_BUILD_PROFILE != "release"
        || EMBEDDED_SUPERVISOR_DIRTY != "false"
        || EMBEDDED_GIT_METADATA_PRESENT != "true"
    {
        return Err(
            "supervisor is not a clean, Git-bound release build and is not admission eligible"
                .to_owned(),
        );
    }
    if arguments.protocol_commit != EMBEDDED_PROTOCOL_COMMIT
        || arguments.source_commit != EMBEDDED_SUPERVISOR_COMMIT
        || arguments.verifier_commit != EMBEDDED_VERIFIER_COMMIT
    {
        return Err("runtime commits differ from the supervisor's embedded provenance".to_owned());
    }
    Ok(())
}

fn normalized_absolute(path: &Path) -> bool {
    let Some(value) = path.to_str() else {
        return false;
    };
    if !path.is_absolute()
        || !value.starts_with('/')
        || value.starts_with("//")
        || (value.len() > 1 && value.ends_with('/'))
        || value
            .chars()
            .any(|character| character < ' ' || character == '<' || character == '>')
    {
        return false;
    }
    value == "/"
        || value[1..]
            .split('/')
            .all(|component| !component.is_empty() && component != "." && component != "..")
}

fn regular_nonsymlink(path: &Path, label: &str) -> Result<()> {
    if !normalized_absolute(path) {
        return Err(format!("{label} must be a normalized absolute path"));
    }
    let metadata = fs::symlink_metadata(path)
        .map_err(|error| format!("inspect {label} {}: {error}", path.display()))?;
    if !metadata.file_type().is_file() || metadata.file_type().is_symlink() {
        return Err(format!("{label} must be a regular non-symlink file"));
    }
    let canonical = fs::canonicalize(path)
        .map_err(|error| format!("canonicalize {label} {}: {error}", path.display()))?;
    if canonical != path {
        return Err(format!(
            "{label} must not traverse a symlink: {}",
            path.display()
        ));
    }
    Ok(())
}

fn directory_nonsymlink(path: &Path, label: &str) -> Result<()> {
    if !normalized_absolute(path) {
        return Err(format!("{label} must be a normalized absolute path"));
    }
    let metadata = fs::symlink_metadata(path)
        .map_err(|error| format!("inspect {label} {}: {error}", path.display()))?;
    if !metadata.file_type().is_dir() || metadata.file_type().is_symlink() {
        return Err(format!("{label} must be a regular non-symlink directory"));
    }
    let canonical = fs::canonicalize(path)
        .map_err(|error| format!("canonicalize {label} {}: {error}", path.display()))?;
    if canonical != path {
        return Err(format!(
            "{label} must not traverse a symlink: {}",
            path.display()
        ));
    }
    Ok(())
}

fn sha256(bytes: &[u8]) -> [u8; 32] {
    Sha256::digest(bytes).into()
}

fn hash_file(path: &Path) -> Result<(u64, String)> {
    let bytes = fs::read(path).map_err(|error| format!("read {}: {error}", path.display()))?;
    Ok((bytes.len() as u64, hex::encode(sha256(&bytes))))
}

fn hash_open_file(file: &fs::File, label: &str) -> Result<(u64, String)> {
    let mut reader = file
        .try_clone()
        .map_err(|error| format!("clone retained {label} descriptor: {error}"))?;
    reader
        .rewind()
        .map_err(|error| format!("rewind retained {label} descriptor: {error}"))?;
    let mut digest = Sha256::new();
    let mut byte_length = 0u64;
    let mut buffer = [0u8; 64 * 1024];
    loop {
        let count = reader
            .read(&mut buffer)
            .map_err(|error| format!("read retained {label} descriptor: {error}"))?;
        if count == 0 {
            break;
        }
        byte_length = byte_length
            .checked_add(u64::try_from(count).map_err(|_| "read length does not fit u64")?)
            .ok_or_else(|| format!("retained {label} byte length overflow"))?;
        digest.update(&buffer[..count]);
    }
    Ok((byte_length, hex::encode(digest.finalize())))
}

#[cfg(target_os = "linux")]
fn sealed_executable_image(
    source: &fs::File,
    expected_length: u64,
    expected_sha256: &str,
    label: &str,
) -> Result<fs::File> {
    use memfd::{FileSeal, MemfdOptions};
    use std::os::unix::fs::PermissionsExt;

    let image = MemfdOptions::default()
        .allow_sealing(true)
        .create(format!("p192-wcm-{label}"))
        .map_err(|error| format!("create sealable {label} image: {error}"))?;
    image
        .as_file()
        .set_permissions(fs::Permissions::from_mode(0o500))
        .map_err(|error| format!("mark sealed {label} image executable: {error}"))?;
    let mut reader = source
        .try_clone()
        .map_err(|error| format!("clone source {label} descriptor: {error}"))?;
    reader
        .rewind()
        .map_err(|error| format!("rewind source {label} descriptor: {error}"))?;
    let mut writer = image.as_file();
    let copied = std::io::copy(&mut reader, &mut writer)
        .map_err(|error| format!("copy {label} into sealed image: {error}"))?;
    image
        .as_file()
        .sync_all()
        .map_err(|error| format!("sync sealed {label} image: {error}"))?;
    let before_seal = hash_open_file(image.as_file(), &format!("sealed {label} image"))?;
    if copied != expected_length
        || before_seal.0 != expected_length
        || before_seal.1 != expected_sha256
    {
        return Err(format!("sealed {label} image differs from admitted bytes"));
    }
    let required_seals = [
        FileSeal::SealShrink,
        FileSeal::SealGrow,
        FileSeal::SealWrite,
        FileSeal::SealSeal,
    ];
    image
        .add_seals(&required_seals)
        .map_err(|error| format!("seal immutable {label} image: {error}"))?;
    let observed_seals = image
        .seals()
        .map_err(|error| format!("inspect immutable {label} image seals: {error}"))?;
    if required_seals
        .iter()
        .any(|seal| !observed_seals.contains(seal))
    {
        return Err(format!("immutable {label} image lacks a required seal"));
    }
    let sealed = image.into_file();
    let after_seal = hash_open_file(&sealed, &format!("sealed {label} image"))?;
    if after_seal != before_seal {
        return Err(format!("sealed {label} image changed while sealing"));
    }
    Ok(sealed)
}

#[cfg(not(target_os = "linux"))]
fn sealed_executable_image(
    _source: &fs::File,
    _expected_length: u64,
    _expected_sha256: &str,
    label: &str,
) -> Result<fs::File> {
    Err(format!(
        "sealed executable image for {label} requires Linux"
    ))
}

#[cfg(unix)]
fn open_retained_executable(path: &Path, label: &str) -> Result<RetainedExecutable> {
    use std::os::unix::fs::{MetadataExt, OpenOptionsExt, PermissionsExt};

    regular_nonsymlink(path, label)?;
    let file = OpenOptions::new()
        .read(true)
        .custom_flags(libc::O_NOFOLLOW)
        .open(path)
        .map_err(|error| format!("open retained {label} {}: {error}", path.display()))?;
    let metadata = file
        .metadata()
        .map_err(|error| format!("inspect retained {label} descriptor: {error}"))?;
    if !metadata.file_type().is_file() || metadata.permissions().mode() & 0o111 == 0 {
        return Err(format!(
            "{label} must be a regular executable file: {}",
            path.display()
        ));
    }
    let path_metadata =
        fs::metadata(path).map_err(|error| format!("reinspect retained {label} path: {error}"))?;
    if metadata.dev() != path_metadata.dev() || metadata.ino() != path_metadata.ino() {
        return Err(format!(
            "{label} path changed while its descriptor was opened"
        ));
    }
    let (byte_length, sha256) = hash_open_file(&file, label)?;
    if byte_length == 0 || byte_length != metadata.len() {
        return Err(format!("{label} has an empty or unstable byte length"));
    }
    let sealed_image = sealed_executable_image(&file, byte_length, &sha256, label)?;
    Ok(RetainedExecutable {
        path: path.to_owned(),
        source_file: file,
        sealed_image,
        byte_length,
        sha256,
    })
}

#[cfg(not(unix))]
fn open_retained_executable(_path: &Path, label: &str) -> Result<RetainedExecutable> {
    Err(format!(
        "retained-descriptor execution for {label} requires Unix"
    ))
}

#[cfg(target_os = "linux")]
fn require_running_executable(retained: &RetainedExecutable) -> Result<()> {
    use std::os::unix::fs::MetadataExt;

    let running = fs::File::open("/proc/self/exe")
        .map_err(|error| format!("open running supervisor image: {error}"))?;
    let running_metadata = running
        .metadata()
        .map_err(|error| format!("inspect running supervisor image: {error}"))?;
    let retained_metadata = retained
        .source_file
        .metadata()
        .map_err(|error| format!("inspect retained supervisor descriptor: {error}"))?;
    if running_metadata.dev() != retained_metadata.dev()
        || running_metadata.ino() != retained_metadata.ino()
    {
        return Err(
            "current executable path does not identify the running supervisor image".to_owned(),
        );
    }
    Ok(())
}

#[cfg(not(target_os = "linux"))]
fn require_running_executable(_retained: &RetainedExecutable) -> Result<()> {
    Err("running-image identity verification requires Linux /proc/self/exe".to_owned())
}

impl RetainedExecutable {
    fn rehash(&self, label: &str) -> Result<String> {
        let (byte_length, sha256) = hash_open_file(&self.sealed_image, label)?;
        if byte_length != self.byte_length {
            return Err(format!("retained {label} byte length changed during run"));
        }
        Ok(sha256)
    }

    #[cfg(unix)]
    fn verify_path_identity(&self, label: &str) -> Result<()> {
        use std::os::unix::fs::MetadataExt;

        let current = open_retained_executable(&self.path, label)?;
        let retained_metadata = self
            .source_file
            .metadata()
            .map_err(|error| format!("inspect retained {label} descriptor: {error}"))?;
        let current_metadata = current
            .source_file
            .metadata()
            .map_err(|error| format!("inspect current {label} descriptor: {error}"))?;
        if retained_metadata.dev() != current_metadata.dev()
            || retained_metadata.ino() != current_metadata.ino()
            || self.byte_length != current.byte_length
            || self.sha256 != current.sha256
        {
            return Err(format!("{label} path identity changed during run"));
        }
        Ok(())
    }

    #[cfg(not(unix))]
    fn verify_path_identity(&self, label: &str) -> Result<()> {
        Err(format!("{label} path identity verification requires Unix"))
    }
}

fn path_text(path: &Path, label: &str) -> Result<String> {
    if !normalized_absolute(path) {
        return Err(format!("{label} must be a normalized absolute path"));
    }
    path.to_str()
        .map(str::to_owned)
        .ok_or_else(|| format!("{label} must be UTF-8"))
}

fn require_absent(path: &Path, label: &str) -> Result<()> {
    match fs::symlink_metadata(path) {
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Ok(_) => Err(format!("{label} already exists: {}", path.display())),
        Err(error) => Err(format!(
            "inspect absent {label} {}: {error}",
            path.display()
        )),
    }
}

fn run_id(run_dir: &Path) -> Result<String> {
    let id = run_dir
        .file_name()
        .and_then(|value| value.to_str())
        .ok_or_else(|| "RUN_DIR must have a UTF-8 basename".to_owned())?;
    let suffix = id
        .strip_prefix("RUN-SCURVE-")
        .ok_or_else(|| "RUN_DIR basename must be RUN-SCURVE-<6hex>".to_owned())?;
    if !lower_hex_digits(suffix, 6) {
        return Err("RUN_DIR basename must be RUN-SCURVE-<6 lowercase hex>".to_owned());
    }
    Ok(id.to_owned())
}

fn resolved_producer_argv(run_dir: &str) -> Vec<String> {
    vec![
        "p192_weighted_cm".to_owned(),
        "preflight".to_owned(),
        "--curve".to_owned(),
        "p192".to_owned(),
        "--sieve-bound".to_owned(),
        "65521".to_owned(),
        "--map-bound".to_owned(),
        "113".to_owned(),
        "--large-prime-bound".to_owned(),
        "2147483647".to_owned(),
        "--positive-control".to_owned(),
        "D=-23".to_owned(),
        "--reference-shell".to_owned(),
        "REF-0".to_owned(),
        "--reference-v-max".to_owned(),
        "2".to_owned(),
        "--reference-x-max-inclusive".to_owned(),
        "1024".to_owned(),
        "--out".to_owned(),
        run_dir.to_owned(),
    ]
}

fn resolved_verifier_argv(run_dir: &str) -> Vec<String> {
    vec![
        "p192_weighted_cm_verify".to_owned(),
        "preflight".to_owned(),
        "--source".to_owned(),
        run_dir.to_owned(),
        "--out".to_owned(),
        format!("{run_dir}/independent-verification.json"),
    ]
}

fn resolved_supervisor_argv(
    executable: &Path,
    arguments: &Arguments,
    run_dir: &str,
) -> Result<Vec<String>> {
    Ok(vec![
        path_text(executable, "supervisor binary")?,
        "--producer".to_owned(),
        path_text(&arguments.producer, "producer binary")?,
        "--verifier".to_owned(),
        path_text(&arguments.verifier, "verifier binary")?,
        "--producer-source".to_owned(),
        path_text(&arguments.producer_source, "producer source")?,
        "--verifier-source".to_owned(),
        path_text(&arguments.verifier_source, "verifier source")?,
        "--protocol-repository".to_owned(),
        path_text(&arguments.protocol_repository, "protocol repository")?,
        "--decision".to_owned(),
        path_text(&arguments.decision, "dispatch decision")?,
        "--run-dir".to_owned(),
        run_dir.to_owned(),
        "--protocol-commit".to_owned(),
        arguments.protocol_commit.clone(),
        "--source-commit".to_owned(),
        arguments.source_commit.clone(),
        "--verifier-commit".to_owned(),
        arguments.verifier_commit.clone(),
        "--audit-out".to_owned(),
        path_text(&arguments.audit_out, "dependency audit output")?,
    ])
}

fn framed_stream_log(domain: &[u8], first: &[u8], second: &[u8]) -> Result<Vec<u8>> {
    let first_length = u64::try_from(first.len())
        .map_err(|_| "first child stream length does not fit u64".to_owned())?;
    let second_length = u64::try_from(second.len())
        .map_err(|_| "second child stream length does not fit u64".to_owned())?;
    let capacity = domain
        .len()
        .checked_add(16)
        .and_then(|value| value.checked_add(first.len()))
        .and_then(|value| value.checked_add(second.len()))
        .ok_or_else(|| "framed child stream length overflow".to_owned())?;
    let mut bytes = Vec::with_capacity(capacity);
    bytes.extend_from_slice(domain);
    bytes.extend_from_slice(&first_length.to_be_bytes());
    bytes.extend_from_slice(first);
    bytes.extend_from_slice(&second_length.to_be_bytes());
    bytes.extend_from_slice(second);
    Ok(bytes)
}

fn canonical_json(value: &Value) -> Result<Vec<u8>> {
    fn write(value: &Value, output: &mut String) -> Result<()> {
        match value {
            Value::Null => output.push_str("null"),
            Value::Bool(value) => output.push_str(if *value { "true" } else { "false" }),
            Value::Number(number) if number.is_i64() || number.is_u64() => {
                output.push_str(&number.to_string());
            }
            Value::Number(_) => return Err("floating-point JSON is forbidden".to_owned()),
            Value::String(text) => output.push_str(
                &serde_json::to_string(text)
                    .map_err(|error| format!("encode JSON string: {error}"))?,
            ),
            Value::Array(values) => {
                output.push('[');
                for (index, item) in values.iter().enumerate() {
                    if index != 0 {
                        output.push(',');
                    }
                    write(item, output)?;
                }
                output.push(']');
            }
            Value::Object(values) => {
                let mut entries = values.iter().collect::<Vec<_>>();
                entries.sort_by(|left, right| left.0.encode_utf16().cmp(right.0.encode_utf16()));
                output.push('{');
                for (index, (key, item)) in entries.into_iter().enumerate() {
                    if index != 0 {
                        output.push(',');
                    }
                    output.push_str(
                        &serde_json::to_string(key)
                            .map_err(|error| format!("encode JSON key: {error}"))?,
                    );
                    output.push(':');
                    write(item, output)?;
                }
                output.push('}');
            }
        }
        Ok(())
    }
    let mut output = String::new();
    write(value, &mut output)?;
    Ok(output.into_bytes())
}

fn write_atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| format!("invalid output filename {}", path.display()))?;
    let temporary = path.with_file_name(format!(".{name}.tmp"));
    let mut file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&temporary)
        .map_err(|error| format!("create {}: {error}", temporary.display()))?;
    file.write_all(bytes)
        .map_err(|error| format!("write {}: {error}", temporary.display()))?;
    file.sync_all()
        .map_err(|error| format!("fsync {}: {error}", temporary.display()))?;
    drop(file);
    fs::hard_link(&temporary, path).map_err(|error| {
        format!(
            "publish {} as {} without replacement: {error}",
            temporary.display(),
            path.display()
        )
    })?;
    fs::remove_file(&temporary).map_err(|error| {
        format!(
            "remove published temporary {}: {error}",
            temporary.display()
        )
    })?;
    let parent_path = path
        .parent()
        .ok_or_else(|| format!("output has no parent: {}", path.display()))?;
    let parent = OpenOptions::new()
        .read(true)
        .open(parent_path)
        .map_err(|error| format!("open {}: {error}", parent_path.display()))?;
    parent
        .sync_all()
        .map_err(|error| format!("fsync {}: {error}", parent_path.display()))
}

fn inventory(run_dir: &Path, paths: &[&str]) -> Result<Vec<Artifact>> {
    paths
        .iter()
        .map(|relative| {
            let path = run_dir.join(relative);
            regular_nonsymlink(&path, relative)?;
            let (byte_length, sha256) = hash_file(&path)?;
            Ok(Artifact {
                path: (*relative).to_owned(),
                byte_length,
                sha256,
            })
        })
        .collect()
}

fn snapshot(run_dir: &Path) -> Result<Vec<Artifact>> {
    inventory(run_dir, &PRODUCER_PATHS)
}

fn exact_terminal_inventory(run_dir: &Path) -> Result<()> {
    fn visit(
        root: &Path,
        current: &Path,
        files: &mut BTreeSet<String>,
        directories: &mut BTreeSet<String>,
    ) -> Result<()> {
        let mut entries = fs::read_dir(current)
            .map_err(|error| format!("read run directory {}: {error}", current.display()))?
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(|error| format!("read run entry: {error}"))?;
        entries.sort_by_key(|entry| entry.file_name());
        for entry in entries {
            let path = entry.path();
            let metadata = fs::symlink_metadata(&path)
                .map_err(|error| format!("inspect run path {}: {error}", path.display()))?;
            if metadata.file_type().is_symlink() {
                return Err(format!("run inventory contains symlink {}", path.display()));
            }
            if metadata.is_dir() {
                let relative = path
                    .strip_prefix(root)
                    .map_err(|error| format!("run directory prefix: {error}"))?
                    .to_str()
                    .ok_or_else(|| "run inventory contains a non-UTF-8 directory".to_owned())?
                    .replace('\\', "/");
                if !directories.insert(relative.clone()) {
                    return Err(format!("duplicate run inventory directory {relative}"));
                }
                visit(root, &path, files, directories)?;
            } else if metadata.is_file() {
                let relative = path
                    .strip_prefix(root)
                    .map_err(|error| format!("run path prefix: {error}"))?
                    .to_str()
                    .ok_or_else(|| "run inventory contains a non-UTF-8 path".to_owned())?
                    .replace('\\', "/");
                if !files.insert(relative.clone()) {
                    return Err(format!("duplicate run inventory path {relative}"));
                }
            } else {
                return Err(format!(
                    "run inventory contains a non-file object {}",
                    path.display()
                ));
            }
        }
        Ok(())
    }

    let mut observed_files = BTreeSet::new();
    let mut observed_directories = BTreeSet::new();
    visit(
        run_dir,
        run_dir,
        &mut observed_files,
        &mut observed_directories,
    )?;
    let mut expected_files = RUN_MANIFEST_PATHS
        .iter()
        .map(|path| (*path).to_owned())
        .collect::<BTreeSet<_>>();
    expected_files.insert("manifest.yaml".to_owned());
    let expected_directories = BTreeSet::from(["reference-box".to_owned()]);
    if observed_files != expected_files || observed_directories != expected_directories {
        return Err(format!(
            "terminal run inventory differs: observed_files={observed_files:?} expected_files={expected_files:?} observed_directories={observed_directories:?} expected_directories={expected_directories:?}"
        ));
    }
    Ok(())
}

fn collect_rs(root: &Path) -> Result<Vec<SourceHashRecord>> {
    fn visit(root: &Path, current: &Path, output: &mut Vec<SourceHashRecord>) -> Result<()> {
        let mut entries = fs::read_dir(current)
            .map_err(|error| format!("read source directory {}: {error}", current.display()))?
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(|error| format!("read source entry: {error}"))?;
        entries.sort_by_key(|entry| entry.file_name());
        for entry in entries {
            let path = entry.path();
            let metadata = fs::symlink_metadata(&path)
                .map_err(|error| format!("inspect source {}: {error}", path.display()))?;
            if metadata.file_type().is_symlink() {
                return Err(format!("source tree contains symlink {}", path.display()));
            }
            if metadata.is_dir() {
                if current == root && entry.file_name() == "target" {
                    continue;
                }
                visit(root, &path, output)?;
            } else if path.extension().and_then(|value| value.to_str()) == Some("rs") {
                let relative = path
                    .strip_prefix(root)
                    .map_err(|error| format!("source prefix: {error}"))?
                    .to_str()
                    .ok_or_else(|| "source path is not UTF-8".to_owned())?
                    .replace('\\', "/");
                let bytes = fs::read(&path)
                    .map_err(|error| format!("read source {}: {error}", path.display()))?;
                output.push(SourceHashRecord {
                    path: relative,
                    sha256: hex::encode(sha256(&bytes)),
                });
            }
        }
        Ok(())
    }
    let mut output = Vec::new();
    visit(root, root, &mut output)?;
    Ok(output)
}

fn rust_identifier_start(character: char) -> bool {
    character == '_' || character.is_alphabetic()
}

fn rust_identifier_continue(character: char) -> bool {
    rust_identifier_start(character) || character.is_alphanumeric()
}

fn raw_literal_end(bytes: &[u8], start: usize, prefix_length: usize) -> Option<usize> {
    let mut cursor = start.checked_add(prefix_length)?;
    let mut hashes = 0usize;
    while bytes.get(cursor) == Some(&b'#') {
        hashes += 1;
        cursor += 1;
    }
    if bytes.get(cursor) != Some(&b'"') {
        return None;
    }
    cursor += 1;
    while cursor < bytes.len() {
        if bytes[cursor] == b'"'
            && cursor
                .checked_add(1 + hashes)
                .is_some_and(|end| end <= bytes.len())
            && bytes[cursor + 1..cursor + 1 + hashes]
                .iter()
                .all(|byte| *byte == b'#')
        {
            return Some(cursor + 1 + hashes);
        }
        cursor += 1;
    }
    None
}

fn cooked_literal_end(bytes: &[u8], quote_at: usize) -> Result<usize> {
    let quote = *bytes
        .get(quote_at)
        .ok_or_else(|| "literal quote offset is out of range".to_owned())?;
    let mut cursor = quote_at + 1;
    let mut escaped = false;
    while cursor < bytes.len() {
        let byte = bytes[cursor];
        if escaped {
            escaped = false;
        } else if byte == b'\\' {
            escaped = true;
        } else if byte == quote {
            return Ok(cursor + 1);
        } else if quote == b'\'' && (byte == b'\n' || byte == b'\r') {
            return Err("unterminated cooked Rust literal".to_owned());
        }
        cursor += 1;
    }
    Err("unterminated cooked Rust literal".to_owned())
}

fn rust_structural_tokens(bytes: &[u8]) -> Result<Vec<String>> {
    let text = std::str::from_utf8(bytes).map_err(|_| "Rust source is not UTF-8".to_owned())?;
    const KEYWORDS: [&str; 52] = [
        "as", "async", "await", "break", "const", "continue", "crate", "dyn", "else", "enum",
        "extern", "false", "fn", "for", "if", "impl", "in", "let", "loop", "match", "mod", "move",
        "mut", "pub", "ref", "return", "self", "Self", "static", "struct", "super", "trait",
        "true", "type", "union", "unsafe", "use", "where", "while", "abstract", "become", "box",
        "do", "final", "macro", "override", "priv", "typeof", "unsized", "virtual", "yield", "try",
    ];
    const OPERATORS: [&str; 30] = [
        "<<=", ">>=", "..=", "...", "::", "->", "=>", "==", "!=", "<=", ">=", "&&", "||", "+=",
        "-=", "*=", "/=", "%=", "^=", "&=", "|=", "<<", ">>", "..", "+", "-", "*", "/", "%", "^",
    ];

    let mut tokens = Vec::new();
    let mut cursor = 0usize;
    let mut block_comment_depth = 0usize;
    while cursor < bytes.len() {
        if block_comment_depth != 0 {
            if bytes[cursor..].starts_with(b"/*") {
                block_comment_depth += 1;
                cursor += 2;
            } else if bytes[cursor..].starts_with(b"*/") {
                block_comment_depth -= 1;
                cursor += 2;
            } else {
                let character = text[cursor..]
                    .chars()
                    .next()
                    .ok_or_else(|| "invalid comment cursor".to_owned())?;
                cursor += character.len_utf8();
            }
            continue;
        }
        if bytes[cursor..].starts_with(b"//") {
            cursor += 2;
            while cursor < bytes.len() && !matches!(bytes[cursor], b'\n' | b'\r') {
                cursor += 1;
            }
            continue;
        }
        if bytes[cursor..].starts_with(b"/*") {
            block_comment_depth = 1;
            cursor += 2;
            continue;
        }
        if bytes[cursor].is_ascii_whitespace() {
            cursor += 1;
            continue;
        }

        let raw_prefix = if bytes[cursor..].starts_with(b"br") || bytes[cursor..].starts_with(b"cr")
        {
            Some(2)
        } else if bytes[cursor] == b'r' {
            Some(1)
        } else {
            None
        };
        if let Some(prefix_length) = raw_prefix {
            if let Some(end) = raw_literal_end(bytes, cursor, prefix_length) {
                tokens.push("LITERAL".to_owned());
                cursor = end;
                continue;
            }
        }

        let cooked_prefix = if matches!(bytes[cursor], b'b' | b'c')
            && matches!(bytes.get(cursor + 1), Some(b'"' | b'\''))
        {
            Some(1)
        } else {
            None
        };
        if let Some(prefix_length) = cooked_prefix {
            let end = cooked_literal_end(bytes, cursor + prefix_length)?;
            tokens.push("LITERAL".to_owned());
            cursor = end;
            continue;
        }
        if bytes[cursor] == b'"' {
            let end = cooked_literal_end(bytes, cursor)?;
            tokens.push("LITERAL".to_owned());
            cursor = end;
            continue;
        }
        if bytes[cursor] == b'\'' {
            if let Ok(end) = cooked_literal_end(bytes, cursor) {
                tokens.push("LITERAL".to_owned());
                cursor = end;
                continue;
            }
        }

        if bytes[cursor].is_ascii_digit() {
            let start = cursor;
            cursor += 1;
            while cursor < bytes.len() {
                let byte = bytes[cursor];
                let numeric_body = byte.is_ascii_alphanumeric() || byte == b'_';
                let decimal_point = byte == b'.'
                    && bytes.get(cursor + 1) != Some(&b'.')
                    && bytes.get(cursor + 1).is_some_and(u8::is_ascii_digit);
                let exponent_sign = matches!(byte, b'+' | b'-')
                    && matches!(bytes.get(cursor.wrapping_sub(1)), Some(b'e' | b'E'));
                if numeric_body || decimal_point || exponent_sign {
                    cursor += 1;
                } else {
                    break;
                }
            }
            if cursor == start {
                return Err("numeric lexer made no progress".to_owned());
            }
            tokens.push("NUMBER".to_owned());
            continue;
        }

        let mut identifier_start = cursor;
        if bytes[cursor..].starts_with(b"r#") {
            identifier_start += 2;
        }
        let first = text[identifier_start..]
            .chars()
            .next()
            .ok_or_else(|| "invalid identifier cursor".to_owned())?;
        if rust_identifier_start(first) {
            let mut end = identifier_start + first.len_utf8();
            while end < bytes.len() {
                let character = text[end..]
                    .chars()
                    .next()
                    .ok_or_else(|| "invalid identifier continuation".to_owned())?;
                if !rust_identifier_continue(character) {
                    break;
                }
                end += character.len_utf8();
            }
            let identifier = &text[identifier_start..end];
            tokens.push(if KEYWORDS.contains(&identifier) {
                identifier.to_owned()
            } else {
                "IDENT".to_owned()
            });
            cursor = end;
            continue;
        }

        if let Some(operator) = OPERATORS
            .iter()
            .find(|operator| text[cursor..].starts_with(**operator))
        {
            tokens.push((*operator).to_owned());
            cursor += operator.len();
            continue;
        }
        let character = text[cursor..]
            .chars()
            .next()
            .ok_or_else(|| "Rust lexer made no progress".to_owned())?;
        tokens.push(character.to_string());
        cursor += character.len_utf8();
    }
    if block_comment_depth != 0 {
        return Err("unterminated nested Rust block comment".to_owned());
    }
    Ok(tokens)
}

fn token_digest(tokens: &[String]) -> Result<[u8; 32]> {
    let mut digest = Sha256::new();
    digest.update(b"P192-WCM-RUST-TOKENS-v1\0");
    for token in tokens {
        let length = u32::try_from(token.len())
            .map_err(|_| "normalized Rust token exceeds u32 bytes".to_owned())?;
        digest.update(length.to_be_bytes());
        digest.update(token.as_bytes());
    }
    Ok(digest.finalize().into())
}

fn clone_source_files(package_root: &Path) -> Result<(PathBuf, Vec<TokenizedSource>, String)> {
    fn visit(root: &Path, current: &Path, output: &mut Vec<TokenizedSource>) -> Result<()> {
        let mut entries = fs::read_dir(current)
            .map_err(|error| format!("read clone-screen root {}: {error}", current.display()))?
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(|error| format!("read clone-screen entry: {error}"))?;
        entries.sort_by_key(|entry| entry.file_name());
        for entry in entries {
            let path = entry.path();
            let metadata = fs::symlink_metadata(&path).map_err(|error| {
                format!("inspect clone-screen path {}: {error}", path.display())
            })?;
            if metadata.file_type().is_symlink() {
                return Err(format!(
                    "clone-screen source tree contains symlink {}",
                    path.display()
                ));
            }
            if metadata.is_dir() {
                visit(root, &path, output)?;
            } else if metadata.is_file()
                && path.extension().and_then(|value| value.to_str()) == Some("rs")
            {
                let relative = path
                    .strip_prefix(root)
                    .map_err(|error| format!("clone-screen source prefix: {error}"))?
                    .to_str()
                    .ok_or_else(|| "clone-screen source path is not UTF-8".to_owned())?
                    .replace('\\', "/");
                let bytes = fs::read(&path).map_err(|error| {
                    format!("read clone-screen source {}: {error}", path.display())
                })?;
                let tokens = rust_structural_tokens(&bytes)
                    .map_err(|error| format!("tokenize clone-screen source {relative}: {error}"))?;
                output.push(TokenizedSource {
                    path: relative,
                    bytes,
                    tokens,
                });
            }
        }
        Ok(())
    }

    let root = package_root.join("src");
    directory_nonsymlink(&root, "clone-screen source root")?;
    let mut sources = Vec::new();
    visit(&root, &root, &mut sources)?;
    sources.sort_by(|left, right| left.path.cmp(&right.path));
    if sources.is_empty() {
        return Err("clone-screen source root contains no Rust files".to_owned());
    }
    let mut tree = Sha256::new();
    tree.update(b"P192-WCM-RUST-SOURCE-TREE-v1\0");
    for source in &sources {
        let path_length = u32::try_from(source.path.len())
            .map_err(|_| "clone-screen path exceeds u32 bytes".to_owned())?;
        let file_length = u64::try_from(source.bytes.len())
            .map_err(|_| "clone-screen file exceeds u64 bytes".to_owned())?;
        tree.update(path_length.to_be_bytes());
        tree.update(source.path.as_bytes());
        tree.update(file_length.to_be_bytes());
        tree.update(&source.bytes);
    }
    Ok((root, sources, hex::encode(tree.finalize())))
}

fn source_clone_screen(producer_source: &Path, verifier_source: &Path) -> Result<Value> {
    const WINDOW: usize = 64;
    const REJECT_AT: usize = 128;
    let (producer_root, producer, producer_tree_sha256) = clone_source_files(producer_source)?;
    let (verifier_root, verifier, verifier_tree_sha256) = clone_source_files(verifier_source)?;

    let mut verifier_windows: BTreeMap<[u8; 32], Vec<(usize, usize)>> = BTreeMap::new();
    for (file_index, source) in verifier.iter().enumerate() {
        if source.tokens.len() < WINDOW {
            continue;
        }
        for start in 0..=source.tokens.len() - WINDOW {
            verifier_windows
                .entry(token_digest(&source.tokens[start..start + WINDOW])?)
                .or_default()
                .push((file_index, start));
        }
    }

    let mut endpoints = BTreeSet::new();
    for (producer_file, source) in producer.iter().enumerate() {
        if source.tokens.len() < WINDOW {
            continue;
        }
        for producer_start in 0..=source.tokens.len() - WINDOW {
            let digest = token_digest(&source.tokens[producer_start..producer_start + WINDOW])?;
            let Some(candidates) = verifier_windows.get(&digest) else {
                continue;
            };
            for &(verifier_file, verifier_start) in candidates {
                let other = &verifier[verifier_file];
                if source.tokens[producer_start..producer_start + WINDOW]
                    != other.tokens[verifier_start..verifier_start + WINDOW]
                {
                    continue;
                }
                let mut producer_left = producer_start;
                let mut verifier_left = verifier_start;
                while producer_left != 0
                    && verifier_left != 0
                    && source.tokens[producer_left - 1] == other.tokens[verifier_left - 1]
                {
                    producer_left -= 1;
                    verifier_left -= 1;
                }
                let mut producer_right = producer_start + WINDOW;
                let mut verifier_right = verifier_start + WINDOW;
                while producer_right < source.tokens.len()
                    && verifier_right < other.tokens.len()
                    && source.tokens[producer_right] == other.tokens[verifier_right]
                {
                    producer_right += 1;
                    verifier_right += 1;
                }
                endpoints.insert((
                    producer_file,
                    producer_left,
                    producer_right,
                    verifier_file,
                    verifier_left,
                    verifier_right,
                ));
            }
        }
    }

    let endpoint_rows = endpoints.into_iter().collect::<Vec<_>>();
    let mut matches = Vec::new();
    for (index, row) in endpoint_rows.iter().enumerate() {
        let contained = endpoint_rows
            .iter()
            .enumerate()
            .any(|(other_index, other)| {
                index != other_index
                    && row.0 == other.0
                    && row.3 == other.3
                    && other.1 <= row.1
                    && other.2 >= row.2
                    && other.4 <= row.4
                    && other.5 >= row.5
                    && (other.1 < row.1 || other.2 > row.2 || other.4 < row.4 || other.5 > row.5)
            });
        if contained {
            continue;
        }
        let token_count = row.2 - row.1;
        if token_count < WINDOW || row.5 - row.4 != token_count {
            return Err("clone-screen endpoint length invariant failed".to_owned());
        }
        matches.push(SourceCloneMatch {
            producer_path: producer[row.0].path.clone(),
            producer_token_start: row.1,
            producer_token_end_exclusive: row.2,
            verifier_path: verifier[row.3].path.clone(),
            verifier_token_start: row.4,
            verifier_token_end_exclusive: row.5,
            token_count,
            normalized_tokens_sha256: hex::encode(token_digest(
                &producer[row.0].tokens[row.1..row.2],
            )?),
        });
    }
    matches.sort_by(|left, right| {
        right
            .token_count
            .cmp(&left.token_count)
            .then_with(|| left.producer_path.cmp(&right.producer_path))
            .then_with(|| left.producer_token_start.cmp(&right.producer_token_start))
            .then_with(|| {
                left.producer_token_end_exclusive
                    .cmp(&right.producer_token_end_exclusive)
            })
            .then_with(|| left.verifier_path.cmp(&right.verifier_path))
            .then_with(|| left.verifier_token_start.cmp(&right.verifier_token_start))
            .then_with(|| {
                left.verifier_token_end_exclusive
                    .cmp(&right.verifier_token_end_exclusive)
            })
    });
    let max_match_tokens = matches.first().map_or(0, |record| record.token_count);
    Ok(json!({
        "claim_limit": "high_similarity_screen_not_proof_of_zero_conceptual_or_common_lineage",
        "excluded_components": ["target"],
        "include": "**/*.rs",
        "matches": matches,
        "max_match_tokens": max_match_tokens,
        "normalization": "rust-structural-v1",
        "overall_status": if max_match_tokens < REJECT_AT { "PASS" } else { "FAIL" },
        "producer_root": path_text(&producer_root, "clone-screen producer root")?,
        "producer_tree_sha256": producer_tree_sha256,
        "reject_at_tokens": REJECT_AT,
        "report_min_tokens": WINDOW,
        "schema": "p192-wcm-source-clone-screen-v1",
        "verifier_root": path_text(&verifier_root, "clone-screen verifier root")?,
        "verifier_tree_sha256": verifier_tree_sha256,
    }))
}

fn git_raw(root: &Path, arguments: &[&str]) -> Result<Output> {
    Command::new("/usr/bin/git")
        .env_clear()
        .env("LC_ALL", "C")
        .env("GIT_CONFIG_NOSYSTEM", "1")
        .env("GIT_CONFIG_GLOBAL", "/dev/null")
        .env("GIT_NO_REPLACE_OBJECTS", "1")
        .env("GIT_GRAFT_FILE", "/dev/null")
        .env("GIT_NO_LAZY_FETCH", "1")
        .env("GIT_OPTIONAL_LOCKS", "0")
        .env("GIT_TERMINAL_PROMPT", "0")
        .env("GIT_PAGER", "cat")
        .env("GIT_LITERAL_PATHSPECS", "1")
        .arg("--no-pager")
        .arg("--no-replace-objects")
        .args(["-c", "core.fsmonitor=false"])
        .args(["-c", "core.commitGraph=false"])
        .args(["-c", "protocol.ext.allow=never"])
        .args(["-c", "core.hooksPath=/dev/null"])
        .args(["-c", "core.attributesFile=/dev/null"])
        .args(["-c", "diff.external="])
        .arg("-C")
        .arg(root)
        .args(arguments)
        .output()
        .map_err(|error| format!("execute git in {}: {error}", root.display()))
}

fn git_text(root: &Path, arguments: &[&str]) -> Result<String> {
    let output = git_raw(root, arguments)?;
    if !output.status.success() {
        return Err(format!(
            "git {:?} failed in {}: {}",
            arguments,
            root.display(),
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    String::from_utf8(output.stdout)
        .map_err(|_| format!("git {:?} output is not UTF-8", arguments))
        .map(|text| text.trim().to_owned())
}

fn git_tree_entry(
    root: &Path,
    commit: &str,
    relative: &str,
    label: &str,
) -> Result<(String, String)> {
    let output = git_raw(root, &["ls-tree", "-z", commit, "--", relative])?;
    if !output.status.success() {
        return Err(format!(
            "inspect {label} at {commit}: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    let records = output
        .stdout
        .split(|byte| *byte == 0)
        .filter(|record| !record.is_empty())
        .collect::<Vec<_>>();
    if records.len() != 1 {
        return Err(format!("{label} must have exactly one Git tree entry"));
    }
    let separator = records[0]
        .iter()
        .position(|byte| *byte == b'\t')
        .ok_or_else(|| format!("{label} has a malformed Git tree entry"))?;
    let (metadata, encoded_path_with_tab) = records[0].split_at(separator);
    let encoded_path = &encoded_path_with_tab[1..];
    let metadata =
        std::str::from_utf8(metadata).map_err(|_| format!("{label} Git metadata is not UTF-8"))?;
    let fields = metadata.split(' ').collect::<Vec<_>>();
    if fields.len() != 3 || fields[1] != "blob" || !lower_hex_digits(fields[2], 40) {
        return Err(format!("{label} is not one full-hash Git blob entry"));
    }
    let observed_path =
        std::str::from_utf8(encoded_path).map_err(|_| format!("{label} Git path is not UTF-8"))?;
    if observed_path != relative {
        return Err(format!("{label} Git path differs from {relative}"));
    }
    Ok((fields[0].to_owned(), fields[2].to_owned()))
}

fn git_blob(root: &Path, object_id: &str, label: &str) -> Result<Vec<u8>> {
    let output = git_raw(root, &["cat-file", "blob", object_id])?;
    if !output.status.success() {
        return Err(format!(
            "read {label} Git blob: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    Ok(output.stdout)
}

fn raw_commit_parents(root: &Path, commit: &str, label: &str) -> Result<Vec<String>> {
    let output = git_raw(root, &["cat-file", "commit", commit])?;
    if !output.status.success() {
        return Err(format!(
            "read raw {label} commit {commit}: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    let header_end = output
        .stdout
        .windows(2)
        .position(|window| window == b"\n\n")
        .ok_or_else(|| format!("raw {label} commit {commit} lacks a header terminator"))?;
    let mut tree_count = 0usize;
    let mut parents = Vec::new();
    for line in output.stdout[..header_end].split(|byte| *byte == b'\n') {
        let Some(separator) = line.iter().position(|byte| *byte == b' ') else {
            continue;
        };
        let kind = &line[..separator];
        let value = &line[separator + 1..];
        if kind == b"tree" {
            let object = std::str::from_utf8(value)
                .map_err(|_| format!("raw {label} commit has a non-UTF-8 tree id"))?;
            if !lower_hex_digits(object, 40) {
                return Err(format!("raw {label} commit has an invalid tree id"));
            }
            tree_count += 1;
        } else if kind == b"parent" {
            let parent = std::str::from_utf8(value)
                .map_err(|_| format!("raw {label} commit has a non-UTF-8 parent id"))?;
            if !lower_hex_digits(parent, 40) {
                return Err(format!("raw {label} commit has an invalid parent id"));
            }
            parents.push(parent.to_owned());
        }
    }
    if tree_count != 1 || parents.iter().collect::<BTreeSet<_>>().len() != parents.len() {
        return Err(format!(
            "raw {label} commit must have one tree and unique parents"
        ));
    }
    Ok(parents)
}

fn raw_commit_is_strict_ancestor(root: &Path, ancestor: &str, descendant: &str) -> Result<bool> {
    if ancestor == descendant {
        return Ok(false);
    }
    let mut pending = raw_commit_parents(root, descendant, "descendant")?;
    let mut visited = BTreeSet::new();
    while let Some(commit) = pending.pop() {
        if commit == ancestor {
            return Ok(true);
        }
        if !visited.insert(commit.clone()) {
            continue;
        }
        if visited.len() > 1_000_000 {
            return Err("raw commit ancestry exceeds the one-million-commit safety cap".to_owned());
        }
        pending.extend(raw_commit_parents(root, &commit, "ancestor traversal")?);
    }
    Ok(false)
}

fn decision_identifier(path: &Path) -> Result<String> {
    let name = path
        .file_name()
        .and_then(|value| value.to_str())
        .ok_or_else(|| "decision path must have a UTF-8 basename".to_owned())?;
    let identifier = name
        .strip_suffix(".yaml")
        .ok_or_else(|| "decision basename must end in .yaml".to_owned())?;
    let Some(suffix) = identifier.strip_prefix("DEC-") else {
        return Err("decision basename must be DEC-YYYYMMDD-<6hex>.yaml".to_owned());
    };
    let Some((date, nonce)) = suffix.split_once('-') else {
        return Err("decision basename must be DEC-YYYYMMDD-<6hex>.yaml".to_owned());
    };
    if date.len() != 8
        || !date.bytes().all(|byte| byte.is_ascii_digit())
        || !lower_hex_digits(nonce, 6)
    {
        return Err("decision basename must be DEC-YYYYMMDD-<6hex>.yaml".to_owned());
    }
    Ok(identifier.to_owned())
}

fn bound_protocol_file(repository: &Path, relative: &str) -> Result<BoundProtocolFile> {
    let path = repository.join(relative);
    regular_nonsymlink(&path, &format!("bound protocol file {relative}"))?;
    let (byte_length, sha256) = hash_file(&path)?;
    if byte_length == 0 {
        return Err(format!("bound protocol file is empty: {relative}"));
    }
    Ok(BoundProtocolFile {
        path: relative.to_owned(),
        byte_length,
        sha256,
    })
}

fn repository_top_level(package: &Path, label: &str) -> Result<String> {
    directory_nonsymlink(package, &format!("{label} package"))?;
    let top_level = PathBuf::from(git_text(package, &["rev-parse", "--show-toplevel"])?);
    directory_nonsymlink(&top_level, &format!("{label} repository"))?;
    path_text(&top_level, &format!("{label} repository"))
}

fn canonical_role1_run_dir(arguments: &Arguments) -> Result<(String, String)> {
    let id = run_id(&arguments.run_dir)?;
    let expected = arguments
        .protocol_repository
        .join("experiments")
        .join(EXPERIMENT_ID)
        .join("runs")
        .join(&id);
    if arguments.run_dir != expected {
        return Err(format!(
            "RUN_DIR must equal the canonical protocol experiment path {}",
            expected.display()
        ));
    }
    Ok((id, path_text(&arguments.run_dir, "RUN_DIR")?))
}

fn create_canonical_run_directory(arguments: &Arguments) -> Result<()> {
    canonical_role1_run_dir(arguments)?;
    require_absent(&arguments.run_dir, "RUN_DIR")?;
    let run_parent = arguments
        .run_dir
        .parent()
        .ok_or_else(|| "RUN_DIR must have a parent".to_owned())?;
    match fs::symlink_metadata(run_parent) {
        Ok(_) => directory_nonsymlink(run_parent, "RUN_DIR parent")?,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            let parent = run_parent
                .parent()
                .ok_or_else(|| "RUN_DIR parent must have an existing parent".to_owned())?;
            directory_nonsymlink(parent, "RUN_DIR parent ancestor")?;
            fs::create_dir(run_parent).map_err(|error| {
                format!("create RUN_DIR parent {}: {error}", run_parent.display())
            })?;
            directory_nonsymlink(run_parent, "created RUN_DIR parent")?;
        }
        Err(error) => {
            return Err(format!(
                "inspect RUN_DIR parent {}: {error}",
                run_parent.display()
            ));
        }
    }
    fs::create_dir(&arguments.run_dir)
        .map_err(|error| format!("create RUN_DIR {}: {error}", arguments.run_dir.display()))?;
    directory_nonsymlink(&arguments.run_dir, "RUN_DIR")
}

fn expected_role1_authorization(
    arguments: &Arguments,
    runtime: &DecisionRuntimeBinding,
) -> Result<Role1Authorization> {
    let (run_id, run_dir) = canonical_role1_run_dir(arguments)?;
    let supervisor_package = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let supervisor_repository = repository_top_level(&supervisor_package, "supervisor")?;
    let producer_repository = repository_top_level(&arguments.producer_source, "producer")?;
    let verifier_repository = repository_top_level(&arguments.verifier_source, "verifier")?;
    if supervisor_package
        != PathBuf::from(&supervisor_repository).join("tools/p192-weighted-cm-supervisor")
        || arguments.producer_source
            != PathBuf::from(&producer_repository).join("tools/p192-weighted-cm")
        || arguments.verifier_source
            != PathBuf::from(&verifier_repository).join("tools/p192-weighted-cm-verify")
    {
        return Err("supervisor, producer, or verifier package path is not canonical".to_owned());
    }
    let distinct_binary_paths = [
        &runtime.supervisor_binary,
        &arguments.producer,
        &arguments.verifier,
    ]
    .into_iter()
    .collect::<BTreeSet<_>>();
    let distinct_binary_hashes = [
        runtime.supervisor_binary_sha256.as_str(),
        runtime.producer_binary_sha256.as_str(),
        runtime.verifier_binary_sha256.as_str(),
    ]
    .into_iter()
    .collect::<BTreeSet<_>>();
    if distinct_binary_paths.len() != 3
        || distinct_binary_hashes.len() != 3
        || arguments.producer_source == arguments.verifier_source
    {
        return Err(
            "supervisor, producer, and verifier paths, hashes, and package roles must be distinct"
                .to_owned(),
        );
    }
    let protocol_repository_path =
        path_text(&arguments.protocol_repository, "protocol repository")?;
    let decision_path = path_text(&arguments.decision, "dispatch decision")?;
    let producer_binary_path = path_text(&arguments.producer, "producer binary")?;
    let verifier_binary_path = path_text(&arguments.verifier, "verifier binary")?;
    let producer_package_path = path_text(&arguments.producer_source, "producer source")?;
    let verifier_package_path = path_text(&arguments.verifier_source, "verifier source")?;
    let supervisor_package_path = path_text(&supervisor_package, "supervisor source")?;
    let supervisor_binary_path = path_text(&runtime.supervisor_binary, "supervisor binary")?;
    let dependency_audit_path = path_text(&arguments.audit_out, "dependency audit output")?;

    Ok(Role1Authorization {
        schema: "p192-wcm-role1-preflight-authorization-v1".to_owned(),
        execution_authorized: true,
        scientific_execution_authorized: false,
        scientific_result_recording_authorized: false,
        experiment_id: EXPERIMENT_ID.to_owned(),
        protocol_version: PROTOCOL_VERSION,
        decision_path,
        stage: StageAuthorization {
            label: "P192-WCM-PREFLIGHT".to_owned(),
            ordinal: 1,
            role: "composite_producer_and_isolated_verifier".to_owned(),
            composite_invocations: 1,
            retry_authorized: false,
            resume_authorized: false,
            advance_to_later_role_authorized: false,
            unauthorized_stage_labels: UNAUTHORIZED_STAGE_LABELS
                .iter()
                .map(|value| (*value).to_owned())
                .collect(),
            custody_outputs: CUSTODY_OUTPUTS
                .iter()
                .map(|value| (*value).to_owned())
                .collect(),
        },
        run: RunAuthorization {
            run_id,
            run_dir,
            must_be_fresh: true,
            existing_path_refused: true,
        },
        protocol: ProtocolAuthorization {
            repository: PROTOCOL_REPOSITORY.to_owned(),
            repository_path: protocol_repository_path,
            commit: arguments.protocol_commit.clone(),
            must_be_ancestor_of_decision_commit: true,
            protected_paths_unchanged_through_decision_commit: true,
            protected_paths: PROTECTED_PROTOCOL_PATHS
                .iter()
                .map(|value| (*value).to_owned())
                .collect(),
            addendum: bound_protocol_file(&arguments.protocol_repository, ROLE1_ADDENDUM_PATH)?,
            interface_manifest: bound_protocol_file(
                &arguments.protocol_repository,
                ROLE1_INTERFACE_MANIFEST_PATH,
            )?,
        },
        executables: ExecutableAuthorization {
            supervisor: SupervisorAuthorization {
                repository: NATIVE_REPOSITORY.to_owned(),
                source_commit: arguments.source_commit.clone(),
                repository_path: supervisor_repository,
                package_path: supervisor_package_path,
                binary_path: supervisor_binary_path,
                byte_length: runtime.supervisor_binary_byte_length,
                binary_sha256: runtime.supervisor_binary_sha256.clone(),
                build_profile: "release".to_owned(),
                worktree_clean: true,
                contains_curve_arithmetic: false,
                embedded_protocol_commit: arguments.protocol_commit.clone(),
                embedded_supervisor_commit: arguments.source_commit.clone(),
                embedded_verifier_commit: arguments.verifier_commit.clone(),
                embedded_build_profile: "release".to_owned(),
                embedded_supervisor_dirty: false,
                embedded_git_metadata_present: true,
            },
            producer: ProducerAuthorization {
                repository: NATIVE_REPOSITORY.to_owned(),
                source_commit: arguments.source_commit.clone(),
                embedded_protocol_commit: arguments.protocol_commit.clone(),
                embedded_source_commit: arguments.source_commit.clone(),
                repository_path: producer_repository,
                package_path: producer_package_path,
                binary_path: producer_binary_path,
                byte_length: runtime.producer_binary_byte_length,
                binary_sha256: runtime.producer_binary_sha256.clone(),
                build_profile: "release".to_owned(),
                worktree_clean: true,
            },
            verifier: VerifierAuthorization {
                repository: NATIVE_REPOSITORY.to_owned(),
                verifier_commit: arguments.verifier_commit.clone(),
                embedded_verifier_commit: arguments.verifier_commit.clone(),
                repository_path: verifier_repository,
                package_path: verifier_package_path,
                binary_path: verifier_binary_path,
                byte_length: runtime.verifier_binary_byte_length,
                binary_sha256: runtime.verifier_binary_sha256.clone(),
                build_profile: "release".to_owned(),
                worktree_clean: true,
                crypto_lib_found: false,
                shared_implementation_components: Vec::new(),
            },
        },
        dependency_audit: DependencyAuditAuthorization {
            path: dependency_audit_path,
            schema: "p192-wcm-dependency-audit-v1".to_owned(),
            must_be_fresh: true,
            receipt_hash_field: "dependency_audit.audit_output_sha256".to_owned(),
            hash_scope: "exact_raw_file_bytes".to_owned(),
        },
        supervisor_argv: runtime.supervisor_argv.clone(),
        producer_argv: runtime.producer_argv.clone(),
        verifier_argv: runtime.verifier_argv.clone(),
        launch: LaunchAuthorization {
            no_path_search_for_any_executable: true,
            runtime_commit_override_forbidden: true,
            producer_before_verifier: true,
            close_and_rehash_producer_outputs_before_verifier: true,
            supervisor_writes_only_run_directory_artifacts: true,
            verify_protocol_checkout_before_run_directory: true,
            derive_and_retain_decision_binding: true,
        },
    })
}

fn nonempty_strings(values: &[String], label: &str) -> Result<()> {
    if values.is_empty() || values.iter().any(String::is_empty) {
        return Err(format!("{label} must contain one or more nonempty strings"));
    }
    Ok(())
}

fn decision_date(value: &str) -> bool {
    value.len() == 10
        && value.as_bytes()[4] == b'-'
        && value.as_bytes()[7] == b'-'
        && value
            .bytes()
            .enumerate()
            .all(|(index, byte)| index == 4 || index == 7 || byte.is_ascii_digit())
}

fn validate_dispatch_decision(
    decision_bytes: &[u8],
    decision_id: &str,
    arguments: &Arguments,
    runtime: &DecisionRuntimeBinding,
) -> Result<()> {
    let document: DispatchDocument = serde_yaml::from_slice(decision_bytes)
        .map_err(|error| format!("parse closed dispatch decision YAML: {error}"))?;
    let decision = document.coordinator_decision;
    if decision.id != decision_id
        || decision.decision != "approve"
        || decision.target_ids != [EXPERIMENT_ID]
        || decision.decided_by != "coordinator"
        || decision.context.is_empty()
        || !decision_date(&decision.decided_at)
        || decision.id.get(4..12) != Some(&decision.decided_at.replace('-', ""))
    {
        return Err("dispatch decision identity, approval, target, context, coordinator, or date binding failed".to_owned());
    }
    nonempty_strings(&decision.rationale, "dispatch rationale")?;
    nonempty_strings(&decision.evidence_refs, "dispatch evidence_refs")?;
    nonempty_strings(&decision.limitations, "dispatch limitations")?;
    nonempty_strings(&decision.next_actions, "dispatch next_actions")?;
    if !decision
        .evidence_refs
        .iter()
        .any(|value| value == ROLE1_ADDENDUM_PATH)
        || !decision
            .evidence_refs
            .iter()
            .any(|value| value == ROLE1_INTERFACE_MANIFEST_PATH)
        || !decision.knowledge_promotion.promoted.is_empty()
        || decision.knowledge_promotion.not_warranted.is_empty()
    {
        return Err("dispatch evidence or knowledge-promotion boundary failed".to_owned());
    }
    let expected = expected_role1_authorization(arguments, runtime)?;
    if decision.authorization != expected {
        return Err(
            "dispatch authorization does not exactly bind the admitted Role-1 runtime".to_owned(),
        );
    }
    Ok(())
}

fn validate_protocol_worktree(arguments: &Arguments, phase: ProtocolBindingPhase) -> Result<()> {
    match phase {
        ProtocolBindingPhase::PreRun => {
            let status = git_text(
                &arguments.protocol_repository,
                &["status", "--porcelain=v1", "--untracked-files=all"],
            )?;
            if !status.is_empty() {
                return Err("protocol repository worktree is dirty before run".to_owned());
            }
        }
        ProtocolBindingPhase::PostRun => {
            let tracked_status = git_text(
                &arguments.protocol_repository,
                &["status", "--porcelain=v1", "--untracked-files=no"],
            )?;
            if !tracked_status.is_empty() {
                return Err("tracked protocol repository state changed during run".to_owned());
            }
            let relative_run_dir = arguments
                .run_dir
                .strip_prefix(&arguments.protocol_repository)
                .map_err(|_| "canonical RUN_DIR is not inside the protocol repository".to_owned())?
                .to_str()
                .ok_or_else(|| "canonical RUN_DIR relative path must be UTF-8".to_owned())?;
            let allowed_prefix = format!("{relative_run_dir}/");
            let untracked = git_raw(
                &arguments.protocol_repository,
                &["ls-files", "--others", "--exclude-standard", "-z"],
            )?;
            if !untracked.status.success() {
                return Err(format!(
                    "inspect post-run untracked protocol paths: {}",
                    String::from_utf8_lossy(&untracked.stderr).trim()
                ));
            }
            for encoded in untracked.stdout.split(|byte| *byte == 0) {
                if encoded.is_empty() {
                    continue;
                }
                let relative = std::str::from_utf8(encoded)
                    .map_err(|_| "post-run untracked protocol path is not UTF-8".to_owned())?;
                if !relative.starts_with(&allowed_prefix) {
                    return Err(format!(
                        "untracked protocol path outside the canonical RUN_DIR: {relative}"
                    ));
                }
            }
        }
    }
    Ok(())
}

fn protocol_decision_binding(
    arguments: &Arguments,
    runtime: &DecisionRuntimeBinding,
    phase: ProtocolBindingPhase,
) -> Result<Value> {
    directory_nonsymlink(&arguments.protocol_repository, "protocol repository")?;
    regular_nonsymlink(&arguments.decision, "dispatch decision")?;
    let top_level = PathBuf::from(git_text(
        &arguments.protocol_repository,
        &["rev-parse", "--show-toplevel"],
    )?);
    if top_level != arguments.protocol_repository {
        return Err("protocol repository path must equal its Git top level".to_owned());
    }
    validate_protocol_worktree(arguments, phase)?;
    let decision_commit = git_text(&arguments.protocol_repository, &["rev-parse", "HEAD"])?;
    if !lower_hex(&decision_commit, 40) {
        return Err("protocol repository HEAD is not a full nonzero lowercase commit".to_owned());
    }
    for (commit, label) in [
        (&arguments.protocol_commit, "protocol commit"),
        (&decision_commit, "decision commit"),
    ] {
        let output = git_raw(
            &arguments.protocol_repository,
            &["cat-file", "-e", &format!("{commit}^{{commit}}")],
        )?;
        if !output.status.success() {
            return Err(format!("{label} does not resolve to an exact Git commit"));
        }
    }
    if arguments.protocol_commit == decision_commit {
        return Err("protocol commit must be a strict ancestor of decision commit".to_owned());
    }
    if !raw_commit_is_strict_ancestor(
        &arguments.protocol_repository,
        &arguments.protocol_commit,
        &decision_commit,
    )? {
        return Err("protocol commit is not a strict ancestor of decision commit".to_owned());
    }

    let decision_id = decision_identifier(&arguments.decision)?;
    let expected_decision = arguments
        .protocol_repository
        .join("ledger")
        .join("decisions")
        .join(format!("{decision_id}.yaml"));
    if arguments.decision != expected_decision {
        return Err("decision path is not the canonical protocol ledger path".to_owned());
    }
    let decision_relative = format!("ledger/decisions/{decision_id}.yaml");
    let (decision_mode, decision_object) = git_tree_entry(
        &arguments.protocol_repository,
        &decision_commit,
        &decision_relative,
        "dispatch decision",
    )?;
    if decision_mode != "100644" {
        return Err("dispatch decision Git mode must be 100644".to_owned());
    }
    let absent = git_raw(
        &arguments.protocol_repository,
        &[
            "ls-tree",
            "-z",
            &arguments.protocol_commit,
            "--",
            &decision_relative,
        ],
    )?;
    if !absent.status.success() || !absent.stdout.is_empty() {
        return Err("dispatch decision must be absent at the protocol commit".to_owned());
    }
    let decision_bytes = fs::read(&arguments.decision)
        .map_err(|error| format!("read dispatch decision: {error}"))?;
    if git_blob(
        &arguments.protocol_repository,
        &decision_object,
        "dispatch decision",
    )? != decision_bytes
    {
        return Err("working dispatch decision differs from its committed blob".to_owned());
    }

    for relative in PROTECTED_PROTOCOL_PATHS {
        let at_protocol = git_tree_entry(
            &arguments.protocol_repository,
            &arguments.protocol_commit,
            relative,
            &format!("protected path {relative} at protocol commit"),
        )?;
        let at_decision = git_tree_entry(
            &arguments.protocol_repository,
            &decision_commit,
            relative,
            &format!("protected path {relative} at decision commit"),
        )?;
        if at_protocol != at_decision {
            return Err(format!(
                "protected protocol Git entry changed through decision commit: {relative}"
            ));
        }
        let current_path = arguments.protocol_repository.join(relative);
        regular_nonsymlink(
            &current_path,
            &format!("protected protocol path {relative}"),
        )?;
        let current = fs::read(&current_path)
            .map_err(|error| format!("read protected protocol path {relative}: {error}"))?;
        if git_blob(
            &arguments.protocol_repository,
            &at_protocol.1,
            &format!("protected path {relative}"),
        )? != current
        {
            return Err(format!(
                "current protected protocol bytes differ from protocol commit: {relative}"
            ));
        }
    }

    validate_dispatch_decision(&decision_bytes, &decision_id, arguments, runtime)?;

    Ok(json!({
        "decision_commit": decision_commit,
        "decision_git_mode": decision_mode,
        "decision_id": decision_id,
        "decision_path": path_text(&arguments.decision, "dispatch decision")?,
        "decision_sha256": hex::encode(sha256(&decision_bytes)),
        "protected_protocol_paths_unchanged": true,
        "protocol_checkout_head": decision_commit,
        "protocol_commit_is_strict_ancestor": true,
        "protocol_repository_path": path_text(&arguments.protocol_repository, "protocol repository")?,
        "protocol_worktree_clean_before_run": true,
    }))
}

fn verify_package_tree_at_commit(
    repository: &Path,
    package: &Path,
    commit: &str,
    label: &str,
) -> Result<()> {
    let relative_package = package
        .strip_prefix(repository)
        .map_err(|_| format!("{label} package is outside its Git repository"))?
        .to_str()
        .ok_or_else(|| format!("{label} package relative path is not UTF-8"))?;
    if relative_package.is_empty() {
        return Err(format!("{label} package may not equal the repository root"));
    }
    let index = git_raw(
        repository,
        &["ls-files", "-v", "-z", "--", relative_package],
    )?;
    if !index.status.success() {
        return Err(format!(
            "inspect {label} package index flags: {}",
            String::from_utf8_lossy(&index.stderr).trim()
        ));
    }
    for record in index.stdout.split(|byte| *byte == 0) {
        if record.is_empty() {
            continue;
        }
        if record.len() < 3 || record[0] != b'H' || record[1] != b' ' {
            return Err(format!(
                "{label} package index contains assume-unchanged, skip-worktree, or noncanonical state"
            ));
        }
    }

    let tree = git_raw(
        repository,
        &["ls-tree", "-r", "-z", commit, "--", relative_package],
    )?;
    if !tree.status.success() {
        return Err(format!(
            "read {label} package tree at {commit}: {}",
            String::from_utf8_lossy(&tree.stderr).trim()
        ));
    }
    let mut expected_paths = BTreeSet::new();
    for record in tree.stdout.split(|byte| *byte == 0) {
        if record.is_empty() {
            continue;
        }
        let separator = record
            .iter()
            .position(|byte| *byte == b'\t')
            .ok_or_else(|| format!("{label} package tree record lacks a path separator"))?;
        let metadata = std::str::from_utf8(&record[..separator])
            .map_err(|_| format!("{label} package tree metadata is not UTF-8"))?;
        let path = std::str::from_utf8(&record[separator + 1..])
            .map_err(|_| format!("{label} package tree path is not UTF-8"))?;
        let fields = metadata.split(' ').collect::<Vec<_>>();
        if fields.len() != 3
            || fields[1] != "blob"
            || !matches!(fields[0], "100644" | "100755")
            || !lower_hex_digits(fields[2], 40)
            || !path.starts_with(&format!("{relative_package}/"))
            || !expected_paths.insert(path.to_owned())
        {
            return Err(format!("{label} package tree contains an invalid entry"));
        }
        let current = repository.join(path);
        let current_metadata = fs::symlink_metadata(&current)
            .map_err(|error| format!("inspect {label} source {path}: {error}"))?;
        if !current_metadata.file_type().is_file() || current_metadata.file_type().is_symlink() {
            return Err(format!("{label} source is not a regular file: {path}"));
        }
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let executable = current_metadata.permissions().mode() & 0o111 != 0;
            if executable != (fields[0] == "100755") {
                return Err(format!("{label} source mode differs from {commit}: {path}"));
            }
        }
        let bytes =
            fs::read(&current).map_err(|error| format!("read {label} source {path}: {error}"))?;
        if bytes != git_blob(repository, fields[2], &format!("{label} source {path}"))? {
            return Err(format!("{label} source bytes differ from {commit}: {path}"));
        }
    }
    if expected_paths.is_empty() {
        return Err(format!("{label} package tree is empty at {commit}"));
    }

    fn collect_current(
        repository: &Path,
        package: &Path,
        current: &Path,
        paths: &mut BTreeSet<String>,
        label: &str,
    ) -> Result<()> {
        let mut entries = fs::read_dir(current)
            .map_err(|error| format!("read {label} package {}: {error}", current.display()))?
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(|error| format!("read {label} package entry: {error}"))?;
        entries.sort_by_key(|entry| entry.file_name());
        for entry in entries {
            let path = entry.path();
            let metadata = fs::symlink_metadata(&path)
                .map_err(|error| format!("inspect {label} package {}: {error}", path.display()))?;
            if metadata.file_type().is_symlink() {
                return Err(format!(
                    "{label} package contains a symlink: {}",
                    path.display()
                ));
            }
            if current == package && entry.file_name() == "target" {
                if !metadata.is_dir() {
                    return Err(format!(
                        "{label} target exclusion is not a directory: {}",
                        path.display()
                    ));
                }
                continue;
            }
            if metadata.is_dir() {
                collect_current(repository, package, &path, paths, label)?;
            } else if metadata.is_file() {
                let relative = path
                    .strip_prefix(repository)
                    .map_err(|error| format!("{label} package prefix: {error}"))?
                    .to_str()
                    .ok_or_else(|| format!("{label} package path is not UTF-8"))?
                    .replace('\\', "/");
                paths.insert(relative);
            } else {
                return Err(format!(
                    "{label} package contains a non-file object: {}",
                    path.display()
                ));
            }
        }
        Ok(())
    }

    let mut current_paths = BTreeSet::new();
    collect_current(repository, package, package, &mut current_paths, label)?;
    if current_paths != expected_paths {
        return Err(format!(
            "{label} package file inventory differs from {commit}"
        ));
    }
    Ok(())
}

fn clean_git_binding(root: &Path, expected_commit: &str, label: &str) -> Result<Value> {
    directory_nonsymlink(root, &format!("{label} package root"))?;
    let canonical_root = fs::canonicalize(root)
        .map_err(|error| format!("canonicalize {label} package root: {error}"))?;
    if canonical_root != root {
        return Err(format!(
            "{label} package root must be canonical and contain no symlink traversal"
        ));
    }
    let top_level = PathBuf::from(git_text(root, &["rev-parse", "--show-toplevel"])?);
    directory_nonsymlink(&top_level, &format!("{label} Git top level"))?;
    let head = git_text(root, &["rev-parse", "HEAD"])?;
    if head != expected_commit {
        return Err(format!(
            "{label} source HEAD {head} differs from bound commit {expected_commit}"
        ));
    }
    let status = git_text(
        root,
        &["status", "--porcelain=v1", "--untracked-files=normal"],
    )?;
    if !status.is_empty() {
        return Err(format!("{label} source worktree is dirty"));
    }
    verify_package_tree_at_commit(&top_level, root, expected_commit, label)?;
    Ok(json!({
        "git_head": head,
        "git_top_level": top_level,
        "package_root": root,
        "worktree_clean": true,
    }))
}

fn dependency_audit(
    producer_source: &Path,
    verifier_source: &Path,
    source_commit: &str,
    verifier_commit: &str,
) -> Result<Vec<u8>> {
    directory_nonsymlink(producer_source, "producer source")?;
    directory_nonsymlink(verifier_source, "verifier source")?;
    if producer_source.file_name().and_then(|name| name.to_str()) != Some("p192-weighted-cm")
        || verifier_source.file_name().and_then(|name| name.to_str())
            != Some("p192-weighted-cm-verify")
        || producer_source == verifier_source
    {
        return Err(
            "producer/verifier source package roots are not the frozen distinct crates".to_owned(),
        );
    }
    let producer_git = clean_git_binding(producer_source, source_commit, "producer")?;
    let verifier_git = clean_git_binding(verifier_source, verifier_commit, "verifier")?;
    let producer_manifest = fs::read(producer_source.join("Cargo.toml"))
        .map_err(|error| format!("read producer Cargo.toml: {error}"))?;
    let verifier_manifest = fs::read(verifier_source.join("Cargo.toml"))
        .map_err(|error| format!("read verifier Cargo.toml: {error}"))?;
    let verifier_lock = fs::read(verifier_source.join("Cargo.lock"))
        .map_err(|error| format!("read verifier Cargo.lock: {error}"))?;
    let verifier_text = String::from_utf8(verifier_manifest.clone())
        .map_err(|_| "verifier Cargo.toml is not UTF-8".to_owned())?;
    for forbidden in ["crypto_lib", "p192-weighted-cm =", "../p192-weighted-cm"] {
        if verifier_text.contains(forbidden) {
            return Err(format!(
                "verifier manifest contains forbidden dependency marker {forbidden:?}"
            ));
        }
    }
    let producer_files = collect_rs(producer_source)?;
    let verifier_files = collect_rs(verifier_source)?;
    let producer_hashes: BTreeSet<_> = producer_files
        .iter()
        .map(|record| record.sha256.as_str())
        .collect();
    let identical_source_hashes = verifier_files
        .iter()
        .filter(|record| producer_hashes.contains(record.sha256.as_str()))
        .map(|record| json!({"path":record.path,"sha256":record.sha256}))
        .collect::<Vec<_>>();
    if !identical_source_hashes.is_empty() {
        return Err(
            "producer and verifier contain byte-identical Rust source components".to_owned(),
        );
    }
    for record in &verifier_files {
        let text = fs::read_to_string(verifier_source.join(&record.path))
            .map_err(|error| format!("read verifier source {}: {error}", record.path))?;
        for forbidden in FORBIDDEN_VERIFIER_SOURCE_MARKERS {
            if text.contains(forbidden) {
                return Err(format!(
                    "verifier source {} contains forbidden marker {forbidden:?}",
                    record.path
                ));
            }
        }
    }
    let source_clone_screen = source_clone_screen(producer_source, verifier_source)?;
    if source_clone_screen["overall_status"] != "PASS" {
        return Err(format!(
            "producer/verifier structural clone screen failed with max_match_tokens={}",
            source_clone_screen["max_match_tokens"]
        ));
    }
    canonical_json(&json!({
        "identical_source_hashes": identical_source_hashes,
        "producer_cargo_toml_sha256": hex::encode(sha256(&producer_manifest)),
        "producer_git": producer_git,
        "producer_rust_sources": producer_files,
        "schema": "p192-wcm-dependency-audit-v1",
        "source_clone_screen": source_clone_screen,
        "shared_implementation_components": [],
        "verifier_cargo_lock_sha256": hex::encode(sha256(&verifier_lock)),
        "verifier_cargo_toml_sha256": hex::encode(sha256(&verifier_manifest)),
        "verifier_crypto_lib_found": false,
        "verifier_git": verifier_git,
        "verifier_rust_sources": verifier_files,
    }))
}

fn validate_runtime_environment() -> Result<()> {
    let mut forbidden = Vec::new();
    for (name, _) in env::vars_os() {
        let name = name
            .to_str()
            .ok_or_else(|| "runtime environment contains a non-UTF-8 variable name".to_owned())?;
        if name.starts_with("LD_") || name.starts_with("DYLD_") {
            forbidden.push(name.to_owned());
        }
    }
    forbidden.sort();
    if !forbidden.is_empty() {
        return Err(format!(
            "forbidden loader override variables are present: {}",
            forbidden.join(",")
        ));
    }
    Ok(())
}

fn run_child(binary: &RetainedExecutable, arguments: &[String]) -> Result<Output> {
    let (argument_zero, remaining) = arguments
        .split_first()
        .ok_or_else(|| "child argv must contain argv[0]".to_owned())?;
    binary.verify_path_identity("child executable")?;
    if binary.rehash("child executable")? != binary.sha256 {
        return Err(format!(
            "retained child executable changed before launch: {}",
            binary.path.display()
        ));
    }
    #[cfg(target_os = "linux")]
    let descriptor_path = {
        use std::os::fd::AsRawFd;
        PathBuf::from(format!("/proc/self/fd/{}", binary.sealed_image.as_raw_fd()))
    };
    #[cfg(not(target_os = "linux"))]
    return Err("retained-descriptor child execution requires Linux /proc/self/fd".to_owned());
    #[cfg(target_os = "linux")]
    let mut command = Command::new(&descriptor_path);
    #[cfg(unix)]
    {
        use std::os::unix::process::CommandExt;
        command.arg0(argument_zero);
    }
    #[cfg(not(unix))]
    let _ = argument_zero;
    command
        .args(remaining)
        .env_remove("P192_WCM_PROTOCOL_COMMIT")
        .env_remove("P192_WCM_SOURCE_COMMIT")
        .env_remove("P192_WCM_VERIFIER_COMMIT")
        .output()
        .map_err(|error| format!("execute retained {}: {error}", binary.path.display()))
}

fn read_envelope(path: &Path) -> Result<(Vec<u8>, VerificationEnvelope)> {
    regular_nonsymlink(path, "verification JSON")?;
    let bytes = fs::read(path).map_err(|error| format!("read {}: {error}", path.display()))?;
    let value: Value = serde_json::from_slice(&bytes)
        .map_err(|error| format!("parse {}: {error}", path.display()))?;
    if canonical_json(&value)? != bytes {
        return Err(format!("{} is not canonical JSON", path.display()));
    }
    let envelope = serde_json::from_value(value)
        .map_err(|error| format!("verification envelope {}: {error}", path.display()))?;
    Ok((bytes, envelope))
}

fn validate_envelope(
    envelope: &VerificationEnvelope,
    expected_schema: &str,
    expected_role: &str,
    arguments: &Arguments,
) -> Result<()> {
    if envelope.schema != expected_schema
        || envelope.experiment_id != EXPERIMENT_ID
        || envelope.protocol_version != PROTOCOL_VERSION
        || envelope.protocol_commit != arguments.protocol_commit
        || envelope.source_commit != arguments.source_commit
        || envelope.role != expected_role
        || envelope.overall_status != "PASS"
    {
        return Err(format!("{expected_role} verification envelope failed"));
    }
    if envelope.bindings.len() != COMPARISON_FIELDS.len()
        || COMPARISON_FIELDS.iter().any(|field| {
            envelope
                .bindings
                .get(*field)
                .is_none_or(|value| !lower_hex(value, 64))
        })
    {
        return Err(format!("{expected_role} binding inventory failed"));
    }
    let expected_checks: &[&str] = match expected_role {
        "producer_self_check" => &PRODUCER_CHECKS,
        "isolated_independent_verifier" => &VERIFIER_CHECKS,
        _ => return Err(format!("unknown verification role {expected_role}")),
    };
    if envelope.checks.len() != expected_checks.len()
        || envelope
            .checks
            .iter()
            .zip(expected_checks)
            .any(|(check, expected_id)| {
                check.as_object().is_none_or(|object| {
                    object.len() != 2
                        || object.get("check_id").and_then(Value::as_str) != Some(*expected_id)
                        || object.get("status").and_then(Value::as_str) != Some("PASS")
                })
            })
    {
        return Err(format!("{expected_role} ordered checks failed"));
    }
    Ok(())
}

fn validate_role_specific_fields(
    producer: &VerificationEnvelope,
    verifier: &VerificationEnvelope,
) -> Result<()> {
    if producer.verifier_commit.is_some()
        || producer.verifier_binary_sha256.is_some()
        || producer.source_artifacts.is_some()
        || producer.independence.is_some()
        || producer.producer_binary_sha256.is_none()
        || producer.artifacts.is_none()
    {
        return Err("producer verification contains the wrong role-specific fields".to_owned());
    }
    if verifier.producer_binary_sha256.is_some()
        || verifier.artifacts.is_some()
        || verifier.verifier_commit.is_none()
        || verifier.verifier_binary_sha256.is_none()
        || verifier.source_artifacts.is_none()
        || verifier.independence.as_ref()
            != Some(&json!({
                "crypto_lib_dependency": false,
                "no_in_memory_producer_state": true,
                "producer_outputs_closed_before_start": true,
                "separate_process": true,
                "shared_implementation_components": [],
            }))
    {
        return Err("independent verification role or independence fields failed".to_owned());
    }
    Ok(())
}

fn child_receipt(
    binary_sha256: &str,
    commit: &str,
    argv: Vec<String>,
    output: &Output,
    started_sequence: u64,
    completed_sequence: u64,
) -> Value {
    json!({
        "argv": argv,
        "binary_sha256": binary_sha256,
        "build_profile": "release",
        "commit": commit,
        "completed_sequence": completed_sequence,
        "exit_code": output.status.code().unwrap_or(-1),
        "started_sequence": started_sequence,
        "stderr_sha256": hex::encode(sha256(&output.stderr)),
        "stdout_sha256": hex::encode(sha256(&output.stdout)),
        "worktree_clean": true,
    })
}

fn execute(arguments: &Arguments, actual_argv0: &Path) -> Result<Value> {
    validate_embedded_provenance(arguments)?;
    validate_runtime_environment()?;
    let executable = env::current_exe().map_err(|error| format!("resolve supervisor: {error}"))?;
    let actual_argv0_text = actual_argv0
        .to_str()
        .ok_or_else(|| "actual argv[0] must be UTF-8".to_owned())?;
    let executable_text = path_text(&executable, "supervisor binary")?;
    if !normalized_absolute(actual_argv0) || actual_argv0_text != executable_text {
        return Err(
            "actual argv[0] must be the normalized absolute current supervisor executable"
                .to_owned(),
        );
    }
    let supervisor_executable = open_retained_executable(&executable, "supervisor binary")?;
    require_running_executable(&supervisor_executable)?;
    let producer_executable = open_retained_executable(&arguments.producer, "producer binary")?;
    let verifier_executable = open_retained_executable(&arguments.verifier, "verifier binary")?;
    let wrapper_byte_length = supervisor_executable.byte_length;
    let wrapper_sha256_before = supervisor_executable.sha256.clone();
    let producer_binary_byte_length = producer_executable.byte_length;
    let producer_binary_sha256_before = producer_executable.sha256.clone();
    let verifier_binary_byte_length = verifier_executable.byte_length;
    let verifier_binary_sha256_before = verifier_executable.sha256.clone();
    if !normalized_absolute(&arguments.run_dir) {
        return Err("RUN_DIR must be a normalized absolute path".to_owned());
    }
    require_absent(&arguments.run_dir, "RUN_DIR")?;
    let (run_id, run_dir) = canonical_role1_run_dir(arguments)?;
    let expected_audit_out = arguments.run_dir.join("dependency-audit.json");
    if !normalized_absolute(&arguments.audit_out) || arguments.audit_out != expected_audit_out {
        return Err(
            "audit output must be the absent normalized RUN_DIR/dependency-audit.json path"
                .to_owned(),
        );
    }
    require_absent(&arguments.audit_out, "dependency audit output")?;
    let supervisor_git = clean_git_binding(
        Path::new(env!("CARGO_MANIFEST_DIR")),
        &arguments.source_commit,
        "supervisor",
    )?;
    let audit_bytes = dependency_audit(
        &arguments.producer_source,
        &arguments.verifier_source,
        &arguments.source_commit,
        &arguments.verifier_commit,
    )?;
    let supervisor_argv = resolved_supervisor_argv(&executable, arguments, &run_dir)?;
    let producer_argv = resolved_producer_argv(&run_dir);
    let verifier_argv = resolved_verifier_argv(&run_dir);
    let decision_runtime = DecisionRuntimeBinding {
        supervisor_binary: executable.clone(),
        supervisor_binary_byte_length: wrapper_byte_length,
        supervisor_binary_sha256: wrapper_sha256_before.clone(),
        supervisor_argv: supervisor_argv.clone(),
        producer_binary_byte_length,
        producer_binary_sha256: producer_binary_sha256_before.clone(),
        producer_argv: producer_argv.clone(),
        verifier_binary_byte_length,
        verifier_binary_sha256: verifier_binary_sha256_before.clone(),
        verifier_argv: verifier_argv.clone(),
    };
    let command_bytes = canonical_json(&json!(supervisor_argv))?;
    let environment = json!({
        "child_environment_policy": "inherit_except_removed_commit_variables_and_forbid_loader_overrides",
        "forbidden_loader_override_names_present": [],
        "producer_binary_sha256": producer_binary_sha256_before,
        "protocol_commit": arguments.protocol_commit,
        "removed_variable_names": [
            "P192_WCM_PROTOCOL_COMMIT",
            "P192_WCM_SOURCE_COMMIT",
            "P192_WCM_VERIFIER_COMMIT",
        ],
        "schema": "p192-wcm-role1-environment-v1",
        "secret_values_recorded": false,
        "source_commit": arguments.source_commit,
        "supervisor_binary_sha256": wrapper_sha256_before,
        "supervisor_build": {
            "embedded_build_profile": EMBEDDED_BUILD_PROFILE,
            "embedded_git_metadata_present": true,
            "embedded_protocol_commit": EMBEDDED_PROTOCOL_COMMIT,
            "embedded_supervisor_commit": EMBEDDED_SUPERVISOR_COMMIT,
            "embedded_supervisor_dirty": false,
            "embedded_verifier_commit": EMBEDDED_VERIFIER_COMMIT,
        },
        "verifier_binary_sha256": verifier_binary_sha256_before,
        "verifier_commit": arguments.verifier_commit,
    });
    let environment_bytes = canonical_json(&environment)?;
    let decision_binding =
        protocol_decision_binding(arguments, &decision_runtime, ProtocolBindingPhase::PreRun)?;
    create_canonical_run_directory(arguments)?;
    let audit_sha256 = hex::encode(sha256(&audit_bytes));
    let producer_output = run_child(&producer_executable, &producer_argv)?;
    if !producer_output.status.success() {
        return Err(format!(
            "producer failed with exit code {:?}; stdout_sha256={} stderr_sha256={}",
            producer_output.status.code(),
            hex::encode(sha256(&producer_output.stdout)),
            hex::encode(sha256(&producer_output.stderr)),
        ));
    }
    let before = snapshot(&arguments.run_dir)?;

    let verifier_output = run_child(&verifier_executable, &verifier_argv)?;
    if !verifier_output.status.success() {
        return Err(format!(
            "verifier failed with exit code {:?}; stdout_sha256={} stderr_sha256={}",
            verifier_output.status.code(),
            hex::encode(sha256(&verifier_output.stdout)),
            hex::encode(sha256(&verifier_output.stderr)),
        ));
    }
    let after = snapshot(&arguments.run_dir)?;
    if before != after {
        return Err("producer artifacts changed while verifier ran".to_owned());
    }

    supervisor_executable.verify_path_identity("supervisor binary")?;
    require_running_executable(&supervisor_executable)?;
    producer_executable.verify_path_identity("producer binary")?;
    verifier_executable.verify_path_identity("verifier binary")?;
    let producer_binary_sha256 = producer_executable.rehash("producer binary")?;
    let verifier_binary_sha256 = verifier_executable.rehash("verifier binary")?;
    let (producer_verification_bytes, producer_verification) =
        read_envelope(&arguments.run_dir.join("verification.json"))?;
    let (independent_verification_bytes, independent_verification) =
        read_envelope(&arguments.run_dir.join("independent-verification.json"))?;
    validate_envelope(
        &producer_verification,
        "p192-wcm-verification-v1",
        "producer_self_check",
        arguments,
    )?;
    validate_envelope(
        &independent_verification,
        "p192-wcm-independent-verification-v1",
        "isolated_independent_verifier",
        arguments,
    )?;
    validate_role_specific_fields(&producer_verification, &independent_verification)?;
    let audit_bytes_after = dependency_audit(
        &arguments.producer_source,
        &arguments.verifier_source,
        &arguments.source_commit,
        &arguments.verifier_commit,
    )?;
    if audit_bytes_after != audit_bytes {
        return Err("producer or verifier source changed during the run".to_owned());
    }
    if clean_git_binding(
        Path::new(env!("CARGO_MANIFEST_DIR")),
        &arguments.source_commit,
        "supervisor",
    )? != supervisor_git
    {
        return Err("supervisor source binding changed during the run".to_owned());
    }
    write_atomic(&arguments.run_dir.join("command.txt"), &command_bytes)?;
    write_atomic(
        &arguments.run_dir.join("environment.json"),
        &environment_bytes,
    )?;
    write_atomic(&arguments.audit_out, &audit_bytes)?;
    if producer_verification.producer_binary_sha256.as_deref() != Some(&producer_binary_sha256)
        || independent_verification.verifier_binary_sha256.as_deref()
            != Some(&verifier_binary_sha256)
        || independent_verification.verifier_commit.as_deref() != Some(&arguments.verifier_commit)
        || producer_verification.artifacts.as_deref() != Some(&before[..9])
        || independent_verification.source_artifacts.as_ref() != Some(&before)
        || producer_binary_sha256_before != producer_binary_sha256
        || verifier_binary_sha256_before != verifier_binary_sha256
    {
        return Err("child binary, commit, inventory, or independence binding failed".to_owned());
    }

    let comparisons = COMPARISON_FIELDS
        .iter()
        .map(|field| {
            let producer = producer_verification.bindings[*field].clone();
            let verifier = independent_verification.bindings[*field].clone();
            json!({"equal":producer == verifier,"field":field,"producer":producer,"verifier":verifier})
        })
        .collect::<Vec<_>>();
    if comparisons.iter().any(|record| record["equal"] != true) {
        return Err("producer/verifier binding comparison failed".to_owned());
    }
    let agreement = json!({
        "comparisons": comparisons,
        "experiment_id": EXPERIMENT_ID,
        "independent_verification_sha256": hex::encode(sha256(&independent_verification_bytes)),
        "overall_status": "PASS",
        "producer_verification_sha256": hex::encode(sha256(&producer_verification_bytes)),
        "protocol_commit": arguments.protocol_commit,
        "protocol_version": PROTOCOL_VERSION,
        "schema": "p192-wcm-independent-agreement-v1",
        "source_commit": arguments.source_commit,
    });
    let agreement_bytes = canonical_json(&agreement)?;
    let agreement_path = arguments.run_dir.join("independent-agreement.json");
    write_atomic(&agreement_path, &agreement_bytes)?;
    let agreement_sha256 = hex::encode(sha256(&agreement_bytes));

    let wrapper_sha256 = supervisor_executable.rehash("supervisor binary")?;
    if wrapper_sha256 != wrapper_sha256_before {
        return Err("supervisor binary changed during the run".to_owned());
    }
    let receipt = json!({
        "agreement_path": "independent-agreement.json",
        "agreement_sha256": agreement_sha256,
        "dependency_audit": {
            "audit_method": "native manifest/source hash comparison with forbidden-import scan",
            "audit_output_sha256": audit_sha256,
            "shared_implementation_components": [],
            "verifier_crypto_lib_found": false,
        },
        "experiment_id": EXPERIMENT_ID,
        "overall_status": "PASS",
        "decision_binding": decision_binding.clone(),
        "producer": child_receipt(
            &producer_binary_sha256,
            &arguments.source_commit,
            producer_argv,
            &producer_output,
            0,
            1,
        ),
        "producer_artifacts_after": after,
        "producer_artifacts_before": before,
        "protocol_commit": arguments.protocol_commit,
        "protocol_version": PROTOCOL_VERSION,
        "run_dir": run_dir,
        "schema": "p192-wcm-independent-verifier-receipt-v1",
        "source_commit": arguments.source_commit,
        "verifier": child_receipt(
            &verifier_binary_sha256,
            &arguments.verifier_commit,
            verifier_argv,
            &verifier_output,
            2,
            3,
        ),
        "verifier_commit": arguments.verifier_commit,
        "wrapper": {"binary_sha256":wrapper_sha256,"contains_curve_arithmetic":false},
    });
    let receipt_bytes = canonical_json(&receipt)?;
    write_atomic(
        &arguments.run_dir.join("independent-verifier-receipt.json"),
        &receipt_bytes,
    )?;
    let stdout_bytes = framed_stream_log(
        STDOUT_LOG_DOMAIN,
        &producer_output.stdout,
        &verifier_output.stdout,
    )?;
    let stderr_bytes = framed_stream_log(
        STDERR_LOG_DOMAIN,
        &producer_output.stderr,
        &verifier_output.stderr,
    )?;
    write_atomic(&arguments.run_dir.join("stdout.log"), &stdout_bytes)?;
    write_atomic(&arguments.run_dir.join("stderr.log"), &stderr_bytes)?;

    let raw_result = json!({
        "agreement_sha256": agreement_sha256,
        "audit_sha256": audit_sha256,
        "experiment_id": EXPERIMENT_ID,
        "protocol_commit": arguments.protocol_commit,
        "receipt_sha256": hex::encode(sha256(&receipt_bytes)),
        "run_dir": run_dir,
        "run_id": run_id,
        "schema": "p192-wcm-role1-raw-result-v1",
        "scientific_result": false,
        "scientific_results": [],
        "source_commit": arguments.source_commit,
        "stage": "P192-WCM-PREFLIGHT",
        "status": "PASS",
        "verifier_commit": arguments.verifier_commit,
    });
    let raw_result_bytes = canonical_json(&raw_result)?;
    write_atomic(
        &arguments.run_dir.join("raw-result.json"),
        &raw_result_bytes,
    )?;

    let manifest_artifacts = inventory(&arguments.run_dir, &RUN_MANIFEST_PATHS)?;
    let manifest = json!({
        "run": {
            "artifacts": manifest_artifacts,
            "code": {
                "command": "command.txt",
                "commit": arguments.source_commit,
                "repository": NATIVE_REPOSITORY,
                "supervisor_binary_sha256": wrapper_sha256,
            },
            "environment": {
                "path": "environment.json",
                "sha256": hex::encode(sha256(&environment_bytes)),
            },
            "experiment_id": EXPERIMENT_ID,
            "id": run_id,
            "inputs": {
                "parameters": {
                    "candidate_count": 1025,
                    "reference_shell": "REF-0",
                    "stage": "P192-WCM-PREFLIGHT",
                },
            },
            "result": {
                "certificate": {"kind": "none"},
                "raw_result": "raw-result.json",
                "scientific_result": false,
            },
            "status": "completed",
            "timing": {"status": "not_recorded_non_scientific_preflight"},
        },
    });
    let manifest_bytes = canonical_json(&manifest)?;
    write_atomic(&arguments.run_dir.join("manifest.yaml"), &manifest_bytes)?;
    exact_terminal_inventory(&arguments.run_dir)?;
    supervisor_executable.verify_path_identity("supervisor binary")?;
    require_running_executable(&supervisor_executable)?;
    producer_executable.verify_path_identity("producer binary")?;
    verifier_executable.verify_path_identity("verifier binary")?;
    if supervisor_executable.rehash("supervisor binary")? != wrapper_sha256_before
        || producer_executable.rehash("producer binary")? != producer_binary_sha256_before
        || verifier_executable.rehash("verifier binary")? != verifier_binary_sha256_before
        || clean_git_binding(
            Path::new(env!("CARGO_MANIFEST_DIR")),
            &arguments.source_commit,
            "supervisor",
        )? != supervisor_git
    {
        return Err("final executable or supervisor-source custody changed".to_owned());
    }
    if protocol_decision_binding(arguments, &decision_runtime, ProtocolBindingPhase::PostRun)?
        != decision_binding
    {
        return Err("protocol repository or dispatch decision changed during the run".to_owned());
    }
    Ok(raw_result)
}

fn main() -> std::process::ExitCode {
    let result = (|| {
        let mut process_argv = env::args_os();
        let actual_argv0 = process_argv
            .next()
            .ok_or_else(|| "process argv must contain argv[0]".to_owned())?
            .into_string()
            .map_err(|_| "actual argv[0] must be UTF-8".to_owned())?;
        let arguments = parse_args(process_argv)?;
        execute(&arguments, Path::new(&actual_argv0))
    })();
    match result {
        Ok(summary) => match canonical_json(&summary) {
            Ok(bytes) => {
                println!("{}", String::from_utf8_lossy(&bytes));
                std::process::ExitCode::SUCCESS
            }
            Err(error) => {
                eprintln!("p192_wcm_preflight_supervisor: {error}");
                std::process::ExitCode::FAILURE
            }
        },
        Err(error) => {
            eprintln!("p192_wcm_preflight_supervisor: {error}");
            std::process::ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

    fn source_package(label: &str, file: &str, source: &str) -> PathBuf {
        let serial = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let root = env::temp_dir().join(format!(
            "p192-wcm-supervisor-{label}-{}-{serial}",
            std::process::id()
        ));
        let path = root.join("src").join(file);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, source).unwrap();
        root
    }

    fn fixture_git(root: &Path, arguments: &[&str]) -> String {
        let output = git_raw(root, arguments).unwrap();
        assert!(
            output.status.success(),
            "git {arguments:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        String::from_utf8(output.stdout).unwrap().trim().to_owned()
    }

    fn dispatch_fixture_document(
        arguments: &Arguments,
        runtime: &DecisionRuntimeBinding,
        decision_id: &str,
    ) -> DispatchDocument {
        DispatchDocument {
            coordinator_decision: CoordinatorDecision {
                id: decision_id.to_owned(),
                context: "Authorize only the frozen non-scientific Role-1 preflight.".to_owned(),
                decision: "approve".to_owned(),
                target_ids: vec![EXPERIMENT_ID.to_owned()],
                rationale: vec!["All closed pre-dispatch bindings are present.".to_owned()],
                evidence_refs: vec![
                    ROLE1_ADDENDUM_PATH.to_owned(),
                    ROLE1_INTERFACE_MANIFEST_PATH.to_owned(),
                ],
                limitations: vec!["No scientific search or later role is authorized.".to_owned()],
                next_actions: vec!["Run the closed pre-dispatch checker once.".to_owned()],
                authorization: expected_role1_authorization(arguments, runtime).unwrap(),
                knowledge_promotion: KnowledgePromotion {
                    promoted: Vec::new(),
                    not_warranted: "A non-scientific preflight is not research evidence."
                        .to_owned(),
                },
                decided_by: "coordinator".to_owned(),
                decided_at: "2026-10-07".to_owned(),
            },
        }
    }

    fn protocol_binding_fixture() -> (PathBuf, Arguments, String, DecisionRuntimeBinding) {
        let serial = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let root = env::temp_dir().join(format!(
            "p192-wcm-protocol-binding-{}-{serial}",
            std::process::id()
        ));
        fs::create_dir(&root).unwrap();
        fixture_git(&root, &["init", "-q"]);
        fixture_git(&root, &["config", "user.email", "fixture@example.invalid"]);
        fixture_git(&root, &["config", "user.name", "Fixture"]);
        for (index, relative) in PROTECTED_PROTOCOL_PATHS.iter().enumerate() {
            let path = root.join(relative);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, format!("protected-{index}\n")).unwrap();
        }
        let binary_dir = root.join("fixture-bin");
        fs::create_dir(&binary_dir).unwrap();
        let supervisor = binary_dir.join("supervisor");
        let producer = binary_dir.join("producer");
        let verifier = binary_dir.join("verifier");
        fs::write(&supervisor, b"supervisor-fixture").unwrap();
        fs::write(&producer, b"producer-fixture").unwrap();
        fs::write(&verifier, b"verifier-fixture").unwrap();
        fixture_git(&root, &["add", "."]);
        fixture_git(&root, &["commit", "-q", "-m", "protocol"]);
        let protocol_commit = fixture_git(&root, &["rev-parse", "HEAD"]);
        let decision_id = "DEC-20261007-abcdef";
        let decision = root
            .join("ledger")
            .join("decisions")
            .join(format!("{decision_id}.yaml"));
        fs::create_dir_all(decision.parent().unwrap()).unwrap();
        let run_dir = root
            .join("experiments")
            .join(EXPERIMENT_ID)
            .join("runs")
            .join("RUN-SCURVE-abcdef");
        let native_root = PathBuf::from(fixture_git(
            Path::new(env!("CARGO_MANIFEST_DIR")),
            &["rev-parse", "--show-toplevel"],
        ));
        let arguments = Arguments {
            producer: producer.clone(),
            verifier: verifier.clone(),
            producer_source: native_root.join("tools/p192-weighted-cm"),
            verifier_source: native_root.join("tools/p192-weighted-cm-verify"),
            protocol_repository: root.clone(),
            decision,
            run_dir: run_dir.clone(),
            protocol_commit: protocol_commit.clone(),
            source_commit: "a".repeat(40),
            verifier_commit: "b".repeat(40),
            audit_out: run_dir.join("dependency-audit.json"),
        };
        let (supervisor_binary_byte_length, supervisor_binary_sha256) =
            hash_file(&supervisor).unwrap();
        let (producer_binary_byte_length, producer_binary_sha256) = hash_file(&producer).unwrap();
        let (verifier_binary_byte_length, verifier_binary_sha256) = hash_file(&verifier).unwrap();
        let run_dir_text = path_text(&arguments.run_dir, "fixture run dir").unwrap();
        let runtime = DecisionRuntimeBinding {
            supervisor_binary: supervisor.clone(),
            supervisor_binary_byte_length,
            supervisor_binary_sha256,
            supervisor_argv: resolved_supervisor_argv(&supervisor, &arguments, &run_dir_text)
                .unwrap(),
            producer_binary_byte_length,
            producer_binary_sha256,
            producer_argv: resolved_producer_argv(&run_dir_text),
            verifier_binary_byte_length,
            verifier_binary_sha256,
            verifier_argv: resolved_verifier_argv(&run_dir_text),
        };
        let document = dispatch_fixture_document(&arguments, &runtime, decision_id);
        fs::write(
            &arguments.decision,
            serde_yaml::to_string(&document).unwrap(),
        )
        .unwrap();
        fixture_git(&root, &["add", "ledger/decisions/DEC-20261007-abcdef.yaml"]);
        fixture_git(&root, &["commit", "-q", "-m", "decision"]);
        (root, arguments, protocol_commit, runtime)
    }

    #[test]
    fn frozen_comparison_fields_are_unique() {
        assert_eq!(COMPARISON_FIELDS.len(), 8);
        assert_eq!(
            COMPARISON_FIELDS.into_iter().collect::<BTreeSet<_>>().len(),
            8
        );
    }

    #[test]
    fn canonical_json_is_compact_and_sorted() {
        let bytes = canonical_json(&json!({"z":1,"a":[true,"x"]})).unwrap();
        assert_eq!(bytes, br#"{"a":[true,"x"],"z":1}"#);
        assert!(!bytes.ends_with(b"\n"));
    }

    #[test]
    fn parser_requires_exact_order_and_full_commits() {
        let commit = "a".repeat(40);
        let values = vec![
            "--producer".to_owned(),
            "/tmp/p".to_owned(),
            "--verifier".to_owned(),
            "/tmp/v".to_owned(),
            "--producer-source".to_owned(),
            "/tmp/ps".to_owned(),
            "--verifier-source".to_owned(),
            "/tmp/vs".to_owned(),
            "--protocol-repository".to_owned(),
            "/tmp/protocol".to_owned(),
            "--decision".to_owned(),
            "/tmp/protocol/ledger/decisions/DEC-20261007-abcdef.yaml".to_owned(),
            "--run-dir".to_owned(),
            "/tmp/run".to_owned(),
            "--protocol-commit".to_owned(),
            commit.clone(),
            "--source-commit".to_owned(),
            commit.clone(),
            "--verifier-commit".to_owned(),
            commit.clone(),
            "--audit-out".to_owned(),
            "/tmp/audit".to_owned(),
        ];
        assert!(parse_args(values.clone()).is_ok());
        let mut bad = values;
        bad[0] = "--verifier".to_owned();
        assert!(parse_args(bad).is_err());
    }

    #[test]
    fn normalized_absolute_paths_are_lexically_exact() {
        assert!(normalized_absolute(Path::new("/absolute/path")));
        for path in [
            "//absolute/path",
            "/absolute//path",
            "/absolute/./path",
            "/absolute/../path",
            "/absolute/path/",
            "/absolute/<placeholder>",
            "/absolute/control\npath",
            "relative/path",
        ] {
            assert!(!normalized_absolute(Path::new(path)), "accepted {path:?}");
        }
    }

    #[test]
    fn decision_identifier_is_exact() {
        assert_eq!(
            decision_identifier(Path::new(
                "/protocol/ledger/decisions/DEC-20261007-abcdef.yaml"
            ))
            .unwrap(),
            "DEC-20261007-abcdef"
        );
        assert!(decision_identifier(Path::new("DEC-20261007-ABCDEF.yaml")).is_err());
        assert!(decision_identifier(Path::new("DEC-2026107-abcdef.yaml")).is_err());
        assert!(decision_identifier(Path::new("DEC-20261007-abcdef.json")).is_err());
    }

    #[test]
    fn canonical_run_directory_is_inside_the_protocol_experiment() {
        let (root, mut arguments, _protocol_commit, _runtime) = protocol_binding_fixture();
        let expected = root
            .join("experiments")
            .join(EXPERIMENT_ID)
            .join("runs")
            .join("RUN-SCURVE-abcdef");
        assert_eq!(
            canonical_role1_run_dir(&arguments).unwrap().1,
            expected.to_str().unwrap()
        );
        assert!(!expected.parent().unwrap().exists());
        create_canonical_run_directory(&arguments).unwrap();
        assert!(expected.is_dir());
        arguments.run_dir = root.join("RUN-SCURVE-abcdef");
        assert!(canonical_role1_run_dir(&arguments)
            .unwrap_err()
            .contains("canonical protocol experiment path"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn post_run_binding_allows_only_the_canonical_run_tree() {
        let (root, arguments, _protocol_commit, runtime) = protocol_binding_fixture();
        let before =
            protocol_decision_binding(&arguments, &runtime, ProtocolBindingPhase::PreRun).unwrap();
        create_canonical_run_directory(&arguments).unwrap();
        fs::write(arguments.run_dir.join("artifact.json"), b"{}\n").unwrap();
        let after =
            protocol_decision_binding(&arguments, &runtime, ProtocolBindingPhase::PostRun).unwrap();
        assert_eq!(before, after);
        fs::write(root.join("outside-run.txt"), b"unexpected\n").unwrap();
        assert!(
            protocol_decision_binding(&arguments, &runtime, ProtocolBindingPhase::PostRun,)
                .unwrap_err()
                .contains("outside the canonical RUN_DIR")
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn retained_descriptor_executes_the_hashed_inode() {
        let path = Path::new("/usr/bin/true");
        if !path.exists() {
            return;
        }
        let retained = open_retained_executable(path, "true fixture").unwrap();
        let mut attempted_writer = retained.sealed_image.try_clone().unwrap();
        attempted_writer.rewind().unwrap();
        assert!(attempted_writer.write_all(b"mutate").is_err());
        let output = run_child(&retained, &["true".to_owned()]).unwrap();
        assert!(output.status.success());
        assert_eq!(retained.rehash("true fixture").unwrap(), retained.sha256);
    }

    #[test]
    fn terminal_inventory_rejects_extra_empty_directories() {
        let serial = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let root = env::temp_dir().join(format!(
            "p192-wcm-terminal-inventory-{}-{serial}",
            std::process::id()
        ));
        fs::create_dir(&root).unwrap();
        for relative in RUN_MANIFEST_PATHS
            .iter()
            .copied()
            .chain(std::iter::once("manifest.yaml"))
        {
            let path = root.join(relative);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, b"fixture").unwrap();
        }
        exact_terminal_inventory(&root).unwrap();
        fs::create_dir(root.join("extra-empty-directory")).unwrap();
        assert!(exact_terminal_inventory(&root)
            .unwrap_err()
            .contains("observed_directories"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn dispatch_authorization_rejects_false_execution_and_runtime_mutations() {
        let (root, arguments, _protocol_commit, runtime) = protocol_binding_fixture();
        let decision_id = decision_identifier(&arguments.decision).unwrap();
        let valid = dispatch_fixture_document(&arguments, &runtime, &decision_id);
        let valid_bytes = serde_yaml::to_string(&valid).unwrap();
        validate_dispatch_decision(valid_bytes.as_bytes(), &decision_id, &arguments, &runtime)
            .unwrap();

        let mut denied = valid.clone();
        denied
            .coordinator_decision
            .authorization
            .execution_authorized = false;
        assert!(validate_dispatch_decision(
            serde_yaml::to_string(&denied).unwrap().as_bytes(),
            &decision_id,
            &arguments,
            &runtime,
        )
        .unwrap_err()
        .contains("does not exactly bind"));

        let mut argv_mutation = valid.clone();
        argv_mutation
            .coordinator_decision
            .authorization
            .supervisor_argv[0] = "/tmp/unbound-supervisor".to_owned();
        assert!(validate_dispatch_decision(
            serde_yaml::to_string(&argv_mutation).unwrap().as_bytes(),
            &decision_id,
            &arguments,
            &runtime,
        )
        .is_err());

        let mut hash_mutation = valid;
        hash_mutation
            .coordinator_decision
            .authorization
            .executables
            .producer
            .binary_sha256 = "0".repeat(64);
        assert!(validate_dispatch_decision(
            serde_yaml::to_string(&hash_mutation).unwrap().as_bytes(),
            &decision_id,
            &arguments,
            &runtime,
        )
        .is_err());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn dispatch_yaml_rejects_duplicate_and_unknown_fields() {
        let (root, arguments, _protocol_commit, runtime) = protocol_binding_fixture();
        let decision_id = decision_identifier(&arguments.decision).unwrap();
        let document = dispatch_fixture_document(&arguments, &runtime, &decision_id);
        let valid = serde_yaml::to_string(&document).unwrap();
        let duplicate = valid.replacen(
            "  decision: approve\n",
            "  decision: approve\n  decision: approve\n",
            1,
        );
        assert!(validate_dispatch_decision(
            duplicate.as_bytes(),
            &decision_id,
            &arguments,
            &runtime,
        )
        .unwrap_err()
        .contains("parse closed dispatch decision YAML"));
        let unknown = valid.replacen(
            "coordinator_decision:\n",
            "coordinator_decision:\n  unexpected: true\n",
            1,
        );
        assert!(
            validate_dispatch_decision(unknown.as_bytes(), &decision_id, &arguments, &runtime,)
                .unwrap_err()
                .contains("parse closed dispatch decision YAML")
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn protocol_decision_binding_requires_exact_committed_custody() {
        let (root, arguments, protocol_commit, runtime) = protocol_binding_fixture();
        let binding =
            protocol_decision_binding(&arguments, &runtime, ProtocolBindingPhase::PreRun).unwrap();
        assert_eq!(binding["decision_id"], "DEC-20261007-abcdef");
        assert_eq!(binding["decision_git_mode"], "100644");
        assert_eq!(binding["protocol_repository_path"], root.to_str().unwrap());
        assert_eq!(binding["protocol_commit_is_strict_ancestor"], true);
        assert_ne!(binding["decision_commit"], protocol_commit);
        assert_eq!(
            binding["protocol_checkout_head"],
            binding["decision_commit"]
        );
        assert_eq!(binding["decision_sha256"].as_str().unwrap().len(), 64);

        let protected = root.join(PROTECTED_PROTOCOL_PATHS[0]);
        let original = fs::read(&protected).unwrap();
        fs::write(&protected, b"tampered\n").unwrap();
        assert!(
            protocol_decision_binding(&arguments, &runtime, ProtocolBindingPhase::PreRun,)
                .unwrap_err()
                .contains("worktree is dirty before run")
        );
        fs::write(&protected, original).unwrap();
        assert!(
            protocol_decision_binding(&arguments, &runtime, ProtocolBindingPhase::PreRun,).is_ok()
        );
        fs::write(&protected, b"committed drift\n").unwrap();
        fixture_git(&root, &["add", PROTECTED_PROTOCOL_PATHS[0]]);
        fixture_git(&root, &["commit", "-q", "-m", "protected drift"]);
        assert!(
            protocol_decision_binding(&arguments, &runtime, ProtocolBindingPhase::PreRun,)
                .unwrap_err()
                .contains("protected protocol Git entry changed")
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn raw_ancestry_rejects_legacy_graft_forgery() {
        let (root, mut arguments, protocol_commit, runtime) = protocol_binding_fixture();
        let protocol_tree = fixture_git(
            &root,
            &["rev-parse", &format!("{protocol_commit}^{{tree}}")],
        );
        let unrelated = fixture_git(
            &root,
            &[
                "commit-tree",
                &protocol_tree,
                "-m",
                "unrelated protocol root",
            ],
        );
        let decision_commit = fixture_git(&root, &["rev-parse", "HEAD"]);
        fs::write(
            root.join(".git/info/grafts"),
            format!("{decision_commit} {unrelated}\n"),
        )
        .unwrap();
        assert!(
            Command::new("/usr/bin/git")
                .arg("-C")
                .arg(&root)
                .args(["merge-base", "--is-ancestor", &unrelated, &decision_commit,])
                .status()
                .unwrap()
                .success(),
            "fixture must demonstrate that Git's legacy graft overlay forges ancestry"
        );
        assert!(!raw_commit_is_strict_ancestor(&root, &unrelated, &decision_commit).unwrap());
        arguments.protocol_commit = unrelated;
        assert!(
            protocol_decision_binding(&arguments, &runtime, ProtocolBindingPhase::PreRun,)
                .unwrap_err()
                .contains("not a strict ancestor")
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn package_snapshot_rejects_assume_unchanged_source_mutation() {
        let serial = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let root = env::temp_dir().join(format!(
            "p192-wcm-package-custody-{}-{serial}",
            std::process::id()
        ));
        let package = root.join("tools/fixture");
        fs::create_dir_all(package.join("src")).unwrap();
        fs::write(package.join("Cargo.toml"), b"[package]\nname='fixture'\n").unwrap();
        fs::write(package.join("src/lib.rs"), b"pub fn admitted() {}\n").unwrap();
        fixture_git(&root, &["init", "-q"]);
        fixture_git(&root, &["config", "user.email", "fixture@example.invalid"]);
        fixture_git(&root, &["config", "user.name", "Fixture"]);
        fixture_git(&root, &["add", "."]);
        fixture_git(&root, &["commit", "-q", "-m", "fixture"]);
        let commit = fixture_git(&root, &["rev-parse", "HEAD"]);
        verify_package_tree_at_commit(&root, &package, &commit, "fixture").unwrap();
        fixture_git(
            &root,
            &[
                "update-index",
                "--assume-unchanged",
                "tools/fixture/src/lib.rs",
            ],
        );
        fs::write(package.join("src/lib.rs"), b"pub fn malicious() {}\n").unwrap();
        assert!(
            verify_package_tree_at_commit(&root, &package, &commit, "fixture")
                .unwrap_err()
                .contains("assume-unchanged")
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn dependency_scan_allows_schema_names_but_rejects_root_crate_use() {
        let allowed = "pub crypto_lib_dependency: bool; // no `crypto_lib` dependency";
        for marker in FORBIDDEN_VERIFIER_SOURCE_MARKERS {
            assert!(!allowed.contains(marker));
        }
        assert!("use crypto_lib::cryptanalysis;".contains("use crypto_lib"));
        assert!("let x = crypto_lib::foo();".contains("crypto_lib::"));
    }

    #[test]
    fn framed_child_logs_have_unambiguous_lengths() {
        let bytes = framed_stream_log(b"domain\0", b"first\n", b"second").unwrap();
        assert_eq!(&bytes[..7], b"domain\0");
        assert_eq!(u64::from_be_bytes(bytes[7..15].try_into().unwrap()), 6);
        assert_eq!(&bytes[15..21], b"first\n");
        assert_eq!(u64::from_be_bytes(bytes[21..29].try_into().unwrap()), 6);
        assert_eq!(&bytes[29..], b"second");
    }

    #[test]
    fn source_rows_are_objects_not_tuple_arrays() {
        let value = serde_json::to_value(SourceHashRecord {
            path: "src/lib.rs".to_owned(),
            sha256: "a".repeat(64),
        })
        .unwrap();
        assert_eq!(value, json!({"path":"src/lib.rs","sha256":"a".repeat(64)}));
    }

    #[test]
    fn structural_clone_screen_catches_renamed_commented_and_moved_code() {
        let mut producer_text = "pub fn producer() {\n".to_owned();
        let mut verifier_text = "pub fn verifier() { /* moved and renamed */\n".to_owned();
        for index in 0..40 {
            producer_text.push_str(&format!("let producer_{index} = {index};\n"));
            verifier_text.push_str(&format!("let verifier_{index} = {}; // noise\n", index + 1));
        }
        producer_text.push_str("}\n");
        verifier_text.push_str("}\n");
        let producer = source_package("clone-producer", "nested/first.rs", &producer_text);
        let verifier = source_package("clone-verifier", "moved/second.rs", &verifier_text);
        let screen = source_clone_screen(&producer, &verifier).unwrap();
        assert_eq!(screen["overall_status"], "FAIL");
        assert!(screen["max_match_tokens"].as_u64().unwrap() >= 128);
        assert_eq!(screen["matches"][0]["producer_path"], "nested/first.rs");
        assert_eq!(screen["matches"][0]["verifier_path"], "moved/second.rs");
        fs::remove_dir_all(producer).unwrap();
        fs::remove_dir_all(verifier).unwrap();
    }

    #[test]
    fn actual_producer_and_verifier_sources_pass_structural_clone_screen() {
        let supervisor_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let tools_root = supervisor_root
            .parent()
            .expect("supervisor package is nested under tools");
        let producer = tools_root.join("p192-weighted-cm");
        let verifier = tools_root.join("p192-weighted-cm-verify");
        let screen = source_clone_screen(&producer, &verifier).unwrap();
        assert_eq!(
            screen["overall_status"], "PASS",
            "actual source clone screen failed: {screen}"
        );
        assert!(
            screen["max_match_tokens"].as_u64().unwrap() < 128,
            "actual source clone screen exceeded the frozen limit: {screen}"
        );
    }

    #[test]
    fn rust_structural_lexer_normalizes_names_literals_and_nested_comments() {
        let left =
            br###"fn alpha() { /* outer /* inner */ done */ let r#value = br##"left"##; }"###;
        let right = br###"fn beta() { let other = br##"right"##; // trailing
        }"###;
        assert_eq!(
            rust_structural_tokens(left).unwrap(),
            rust_structural_tokens(right).unwrap()
        );
    }

    #[test]
    fn run_identifier_is_exact() {
        assert_eq!(
            run_id(Path::new("/tmp/RUN-SCURVE-000000")).unwrap(),
            "RUN-SCURVE-000000"
        );
        assert!(run_id(Path::new("/tmp/RUN-SCURVE-ABCDEF")).is_err());
        assert!(run_id(Path::new("/tmp/not-a-run")).is_err());
    }
}
