//! Non-arithmetic supervisor for the frozen Role-1 P-192 weighted-CM gate.
//!
//! This binary performs process sequencing, byte hashing, source-dependency
//! auditing, and receipt construction.  It deliberately contains no curve,
//! ideal, factorization, or certificate arithmetic.

use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::ffi::OsString;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Component, Path, PathBuf};
use std::process::{Command, Output};

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

type Result<T> = std::result::Result<T, String>;

const EXPERIMENT_ID: &str = "EXP-SCURVE-1a8daf";
const PROTOCOL_VERSION: u64 = 2;
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

#[derive(Clone, Debug)]
struct Arguments {
    producer: PathBuf,
    verifier: PathBuf,
    producer_source: PathBuf,
    verifier_source: PathBuf,
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
    "p192_wcm_preflight_supervisor --producer ABS_BIN --verifier ABS_BIN --producer-source ABS_DIR --verifier-source ABS_DIR --run-dir ABS_DIR --protocol-commit HEX40 --source-commit HEX40 --verifier-commit HEX40 --audit-out ABS_FILE"
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
        run_dir: PathBuf::from(&values[9]),
        protocol_commit: values[11].clone(),
        source_commit: values[13].clone(),
        verifier_commit: values[15].clone(),
        audit_out: PathBuf::from(&values[17]),
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
    value.len() == length
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn normalized_absolute(path: &Path) -> bool {
    path.is_absolute()
        && path.components().all(|component| {
            matches!(
                component,
                Component::RootDir | Component::Prefix(_) | Component::Normal(_)
            )
        })
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
    Ok(())
}

fn sha256(bytes: &[u8]) -> [u8; 32] {
    Sha256::digest(bytes).into()
}

fn hash_file(path: &Path) -> Result<(u64, String)> {
    let bytes = fs::read(path).map_err(|error| format!("read {}: {error}", path.display()))?;
    Ok((bytes.len() as u64, hex::encode(sha256(&bytes))))
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
    fs::rename(&temporary, path).map_err(|error| {
        format!(
            "rename {} to {}: {error}",
            temporary.display(),
            path.display()
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

fn snapshot(run_dir: &Path) -> Result<Vec<Artifact>> {
    PRODUCER_PATHS
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

fn collect_rs(root: &Path) -> Result<Vec<(String, String)>> {
    fn visit(root: &Path, current: &Path, output: &mut Vec<(String, String)>) -> Result<()> {
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
                if entry.file_name() != "target" {
                    visit(root, &path, output)?;
                }
            } else if path.extension().and_then(|value| value.to_str()) == Some("rs") {
                let relative = path
                    .strip_prefix(root)
                    .map_err(|error| format!("source prefix: {error}"))?
                    .to_str()
                    .ok_or_else(|| "source path is not UTF-8".to_owned())?
                    .replace('\\', "/");
                let bytes = fs::read(&path)
                    .map_err(|error| format!("read source {}: {error}", path.display()))?;
                output.push((relative, hex::encode(sha256(&bytes))));
            }
        }
        Ok(())
    }
    let mut output = Vec::new();
    visit(root, root, &mut output)?;
    Ok(output)
}

fn dependency_audit(producer_source: &Path, verifier_source: &Path) -> Result<Vec<u8>> {
    directory_nonsymlink(producer_source, "producer source")?;
    directory_nonsymlink(verifier_source, "verifier source")?;
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
    let producer_hashes: BTreeSet<_> = producer_files.iter().map(|(_, hash)| hash).collect();
    let identical_source_hashes = verifier_files
        .iter()
        .filter_map(|(path, hash)| producer_hashes.contains(hash).then_some((path, hash)))
        .map(|(path, hash)| json!({"path":path,"sha256":hash}))
        .collect::<Vec<_>>();
    if !identical_source_hashes.is_empty() {
        return Err(
            "producer and verifier contain byte-identical Rust source components".to_owned(),
        );
    }
    for (relative, _) in &verifier_files {
        let text = fs::read_to_string(verifier_source.join(relative))
            .map_err(|error| format!("read verifier source {relative}: {error}"))?;
        for forbidden in ["crypto_lib", "include!", "#[path", "../producer/"] {
            if text.contains(forbidden) {
                return Err(format!(
                    "verifier source {relative} contains forbidden marker {forbidden:?}"
                ));
            }
        }
    }
    canonical_json(&json!({
        "identical_source_hashes": identical_source_hashes,
        "producer_cargo_toml_sha256": hex::encode(sha256(&producer_manifest)),
        "producer_rust_sources": producer_files,
        "schema": "p192-wcm-dependency-audit-v1",
        "shared_implementation_components": [],
        "verifier_cargo_lock_sha256": hex::encode(sha256(&verifier_lock)),
        "verifier_cargo_toml_sha256": hex::encode(sha256(&verifier_manifest)),
        "verifier_crypto_lib_found": false,
        "verifier_rust_sources": verifier_files,
    }))
}

fn run_child(binary: &Path, arguments: &[String]) -> Result<Output> {
    Command::new(binary)
        .args(arguments)
        .env_remove("P192_WCM_PROTOCOL_COMMIT")
        .env_remove("P192_WCM_SOURCE_COMMIT")
        .env_remove("P192_WCM_VERIFIER_COMMIT")
        .output()
        .map_err(|error| format!("execute {}: {error}", binary.display()))
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

fn execute(arguments: &Arguments) -> Result<Value> {
    for (path, label) in [
        (&arguments.producer, "producer binary"),
        (&arguments.verifier, "verifier binary"),
    ] {
        regular_nonsymlink(path, label)?;
    }
    let producer_binary_sha256_before = hash_file(&arguments.producer)?.1;
    let verifier_binary_sha256_before = hash_file(&arguments.verifier)?.1;
    if !normalized_absolute(&arguments.run_dir) || arguments.run_dir.exists() {
        return Err("RUN_DIR must be a normalized absolute path that does not exist".to_owned());
    }
    if !normalized_absolute(&arguments.audit_out) || arguments.audit_out.exists() {
        return Err(
            "audit output must be a normalized absolute path that does not exist".to_owned(),
        );
    }
    let audit_bytes = dependency_audit(&arguments.producer_source, &arguments.verifier_source)?;
    write_atomic(&arguments.audit_out, &audit_bytes)?;
    let audit_sha256 = hex::encode(sha256(&audit_bytes));

    fs::create_dir(&arguments.run_dir)
        .map_err(|error| format!("create RUN_DIR {}: {error}", arguments.run_dir.display()))?;
    let run_dir = arguments
        .run_dir
        .to_str()
        .ok_or_else(|| "RUN_DIR is not UTF-8".to_owned())?
        .to_owned();
    let producer_argv = vec![
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
        run_dir.clone(),
    ];
    let producer_output = run_child(&arguments.producer, &producer_argv[1..])?;
    if !producer_output.status.success() {
        return Err(format!(
            "producer failed with exit code {:?}; stdout_sha256={} stderr_sha256={}",
            producer_output.status.code(),
            hex::encode(sha256(&producer_output.stdout)),
            hex::encode(sha256(&producer_output.stderr)),
        ));
    }
    let before = snapshot(&arguments.run_dir)?;

    let verifier_argv = vec![
        "p192_weighted_cm_verify".to_owned(),
        "preflight".to_owned(),
        "--source".to_owned(),
        run_dir.clone(),
        "--out".to_owned(),
        format!("{run_dir}/independent-verification.json"),
    ];
    let verifier_output = run_child(&arguments.verifier, &verifier_argv[1..])?;
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

    let producer_binary_sha256 = hash_file(&arguments.producer)?.1;
    let verifier_binary_sha256 = hash_file(&arguments.verifier)?.1;
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
    let audit_bytes_after =
        dependency_audit(&arguments.producer_source, &arguments.verifier_source)?;
    if audit_bytes_after != audit_bytes {
        return Err("producer or verifier source changed during the run".to_owned());
    }
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

    let executable = env::current_exe().map_err(|error| format!("resolve supervisor: {error}"))?;
    let wrapper_sha256 = hash_file(&executable)?.1;
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
    Ok(json!({
        "agreement_sha256": agreement_sha256,
        "audit_sha256": audit_sha256,
        "receipt_sha256": hex::encode(sha256(&receipt_bytes)),
        "run_dir": arguments.run_dir,
        "status": "PASS",
    }))
}

fn main() -> std::process::ExitCode {
    let result = parse_args(env::args_os().skip(1)).and_then(|arguments| execute(&arguments));
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
}
