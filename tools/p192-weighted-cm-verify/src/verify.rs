use std::{
    collections::BTreeSet,
    fs::{self, OpenOptions},
    io::{Read, Write},
    path::{Component, Path, PathBuf},
};

use serde::{Deserialize, Serialize};

use crate::{
    controls::verify_controls_bytes,
    encoding::{canonical_json, parse_canonical_json, sha256},
    factor_base::verify_factor_base_bytes,
    identity::{verify_identity_and_maximal_order_bytes, Verdict},
    ref0::{verify_reference_box_bytes, ArtifactRecord},
    Result,
};

pub const SOURCE_PATHS: [&str; 10] = [
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

const PRODUCER_ARTIFACT_PATHS: [&str; 9] = [
    "identity.json",
    "maximal-order.json",
    "factor-base.json",
    "factor-base.sha256",
    "controls.json",
    "reference-box.json",
    "reference-box/candidate-records.bin",
    "reference-box/disposition.bin",
    "reference-box/certificates.bin",
];

const PRODUCER_CHECK_IDS: [&str; 11] = [
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

fn running_verifier_digest() -> Result<[u8; 32]> {
    #[cfg(target_os = "linux")]
    {
        let mut image = fs::File::open("/proc/self/exe")
            .map_err(|error| format!("open running verifier image: {error}"))?;
        let mut bytes = Vec::new();
        image
            .read_to_end(&mut bytes)
            .map_err(|error| format!("stream running verifier image: {error}"))?;
        Ok(sha256(&bytes))
    }
    #[cfg(not(target_os = "linux"))]
    {
        Err("running verifier image hashing requires Linux /proc/self/exe".to_owned())
    }
}

const INDEPENDENT_CHECK_IDS: [&str; 15] = [
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

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BindingValues {
    pub curve_tuple_sha256: String,
    pub cm_order_tuple_sha256: String,
    pub factor_base_sha256: String,
    pub control_results_sha256: String,
    pub reference_candidate_shard_sha256: String,
    pub reference_candidate_shell_sha256: String,
    pub reference_disposition_shard_sha256: String,
    pub reference_disposition_shell_sha256: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ProducerVerification {
    schema: String,
    experiment_id: String,
    protocol_version: u64,
    protocol_commit: String,
    source_commit: String,
    role: String,
    producer_binary_sha256: String,
    artifacts: Vec<ArtifactRecord>,
    bindings: BindingValues,
    checks: Vec<Verdict>,
    overall_status: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Independence {
    pub separate_process: bool,
    pub producer_outputs_closed_before_start: bool,
    pub no_in_memory_producer_state: bool,
    pub crypto_lib_dependency: bool,
    pub shared_implementation_components: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IndependentVerification {
    pub schema: String,
    pub experiment_id: String,
    pub protocol_version: u64,
    pub protocol_commit: String,
    pub source_commit: String,
    pub verifier_commit: String,
    pub role: String,
    pub verifier_binary_sha256: String,
    pub source_artifacts: Vec<ArtifactRecord>,
    pub bindings: BindingValues,
    pub independence: Independence,
    pub checks: Vec<Verdict>,
    pub overall_status: String,
}

#[derive(Clone, Debug)]
struct SourceSnapshot {
    path: &'static str,
    bytes: Vec<u8>,
}

fn lower_hex(value: &str, length: usize) -> bool {
    value.len() == length
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn commit_hex(value: &str) -> bool {
    lower_hex(value, 40) && value.bytes().any(|byte| byte != b'0')
}

fn read_regular(source: &Path, relative: &'static str) -> Result<SourceSnapshot> {
    let path = source.join(relative);
    let metadata = fs::symlink_metadata(&path)
        .map_err(|error| format!("metadata {}: {error}", path.display()))?;
    if !metadata.file_type().is_file() || metadata.file_type().is_symlink() {
        return Err(format!(
            "{} is not a regular non-symlink file",
            path.display()
        ));
    }
    let bytes = fs::read(&path).map_err(|error| format!("read {}: {error}", path.display()))?;
    Ok(SourceSnapshot {
        path: relative,
        bytes,
    })
}

fn exact_directory_entries(path: &Path, expected: &[&str]) -> Result<()> {
    let mut actual = BTreeSet::new();
    for entry in
        fs::read_dir(path).map_err(|error| format!("read_dir {}: {error}", path.display()))?
    {
        let entry = entry.map_err(|error| format!("read_dir entry {}: {error}", path.display()))?;
        let name = entry
            .file_name()
            .into_string()
            .map_err(|_| format!("{} contains a non-UTF8 filename", path.display()))?;
        actual.insert(name);
    }
    let expected: BTreeSet<_> = expected.iter().map(|value| (*value).to_owned()).collect();
    if actual != expected {
        return Err(format!(
            "{} has a noncanonical artifact inventory",
            path.display()
        ));
    }
    Ok(())
}

fn verify_inventory(source: &Path) -> Result<()> {
    let metadata = fs::symlink_metadata(source)
        .map_err(|error| format!("metadata {}: {error}", source.display()))?;
    if !metadata.file_type().is_dir() || metadata.file_type().is_symlink() {
        return Err("source is not a regular non-symlink directory".to_owned());
    }
    exact_directory_entries(
        source,
        &[
            "identity.json",
            "maximal-order.json",
            "factor-base.json",
            "factor-base.sha256",
            "controls.json",
            "reference-box.json",
            "reference-box",
            "verification.json",
        ],
    )?;
    let reference = source.join("reference-box");
    let metadata = fs::symlink_metadata(&reference)
        .map_err(|error| format!("metadata {}: {error}", reference.display()))?;
    if !metadata.file_type().is_dir() || metadata.file_type().is_symlink() {
        return Err("reference-box is not a regular non-symlink directory".to_owned());
    }
    exact_directory_entries(
        &reference,
        &[
            "candidate-records.bin",
            "disposition.bin",
            "certificates.bin",
        ],
    )
}

fn artifact(snapshot: &SourceSnapshot) -> ArtifactRecord {
    ArtifactRecord {
        path: snapshot.path.to_owned(),
        byte_length: snapshot.bytes.len() as u64,
        sha256: hex::encode(sha256(&snapshot.bytes)),
    }
}

fn require_passes(checks: &[Verdict], ids: &[&str], overall: &str) -> Result<()> {
    if overall != "PASS" || checks.len() != ids.len() {
        return Err("verification verdict array is not an overall PASS".to_owned());
    }
    for (check, expected) in checks.iter().zip(ids) {
        if check.check_id != *expected || check.status != "PASS" {
            return Err(format!(
                "verification check {expected} is missing, reordered, or failed"
            ));
        }
    }
    Ok(())
}

fn parse_producer_verification(
    bytes: &[u8],
    protocol_commit: &str,
    source_commit: &str,
    snapshots: &[SourceSnapshot],
    bindings: &BindingValues,
) -> Result<()> {
    let value = parse_canonical_json(bytes)?;
    let manifest: ProducerVerification = serde_json::from_value(value)
        .map_err(|error| format!("producer verification schema: {error}"))?;
    if manifest.schema != "p192-wcm-verification-v1"
        || manifest.experiment_id != "EXP-SCURVE-1a8daf"
        || manifest.protocol_version != 2
        || manifest.protocol_commit != protocol_commit
        || manifest.source_commit != source_commit
        || manifest.role != "producer_self_check"
        || !lower_hex(&manifest.producer_binary_sha256, 64)
        || manifest.artifacts.len() != PRODUCER_ARTIFACT_PATHS.len()
        || manifest.bindings != *bindings
    {
        return Err("producer verification envelope/cross-binding mismatch".to_owned());
    }
    for ((record, expected_path), snapshot) in manifest
        .artifacts
        .iter()
        .zip(PRODUCER_ARTIFACT_PATHS)
        .zip(snapshots.iter())
    {
        if record != &artifact(snapshot) || record.path != expected_path {
            return Err(format!(
                "producer artifact binding mismatch for {expected_path}"
            ));
        }
    }
    require_passes(
        &manifest.checks,
        &PRODUCER_CHECK_IDS,
        &manifest.overall_status,
    )
}

fn snapshots_stable(source: &Path, before: &[SourceSnapshot]) -> Result<()> {
    for original in before {
        let after = read_regular(source, original.path)?;
        if after.bytes != original.bytes {
            return Err(format!(
                "source artifact changed during verification: {}",
                original.path
            ));
        }
    }
    Ok(())
}

pub fn verify_source(
    source: &Path,
    verifier_commit: &str,
    verifier_binary_sha256: [u8; 32],
) -> Result<IndependentVerification> {
    if !commit_hex(verifier_commit) {
        return Err("embedded verifier commit is not lowercase 40-hex".to_owned());
    }
    verify_inventory(source)?;
    let snapshots: Vec<_> = SOURCE_PATHS
        .iter()
        .map(|path| read_regular(source, path))
        .collect::<Result<_>>()?;

    let (identity, maximal_order) =
        verify_identity_and_maximal_order_bytes(&snapshots[0].bytes, &snapshots[1].bytes)?;
    let protocol_commit = identity.manifest.protocol_commit.clone();
    let source_commit = identity.manifest.source_commit.clone();
    if !commit_hex(&protocol_commit) || !commit_hex(&source_commit) {
        return Err("producer protocol/source commit is not lowercase 40-hex".to_owned());
    }
    let factor_base = verify_factor_base_bytes(&snapshots[2].bytes, &snapshots[3].bytes)?;
    if factor_base.manifest.source_commit != source_commit {
        return Err("factor-base source_commit mismatch".to_owned());
    }
    let reference = verify_reference_box_bytes(
        snapshots[5].bytes.clone(),
        &snapshots[6].bytes,
        &snapshots[7].bytes,
        &snapshots[8].bytes,
        &protocol_commit,
        &source_commit,
        &factor_base,
    )?;
    let controls = verify_controls_bytes(
        &snapshots[4].bytes,
        &protocol_commit,
        &source_commit,
        &factor_base,
        &reference.regenerated,
    )?;

    let bindings = BindingValues {
        curve_tuple_sha256: hex::encode(identity.curve_tuple_sha256),
        cm_order_tuple_sha256: hex::encode(maximal_order.cm_order_tuple_sha256),
        factor_base_sha256: hex::encode(factor_base.sha256),
        control_results_sha256: hex::encode(controls.results_sha256),
        reference_candidate_shard_sha256: hex::encode(reference.regenerated.candidate_shard_sha256),
        reference_candidate_shell_sha256: hex::encode(reference.regenerated.candidate_shell_sha256),
        reference_disposition_shard_sha256: hex::encode(
            reference.regenerated.disposition_shard_sha256,
        ),
        reference_disposition_shell_sha256: hex::encode(
            reference.regenerated.disposition_shell_sha256,
        ),
    };
    parse_producer_verification(
        &snapshots[9].bytes,
        &protocol_commit,
        &source_commit,
        &snapshots,
        &bindings,
    )?;
    snapshots_stable(source, &snapshots)?;

    Ok(IndependentVerification {
        schema: "p192-wcm-independent-verification-v1".to_owned(),
        experiment_id: "EXP-SCURVE-1a8daf".to_owned(),
        protocol_version: 2,
        protocol_commit,
        source_commit,
        verifier_commit: verifier_commit.to_owned(),
        role: "isolated_independent_verifier".to_owned(),
        verifier_binary_sha256: hex::encode(verifier_binary_sha256),
        source_artifacts: snapshots.iter().map(artifact).collect(),
        bindings,
        independence: Independence {
            separate_process: true,
            producer_outputs_closed_before_start: true,
            no_in_memory_producer_state: true,
            crypto_lib_dependency: false,
            shared_implementation_components: Vec::new(),
        },
        checks: INDEPENDENT_CHECK_IDS
            .iter()
            .map(|check_id| Verdict {
                check_id: (*check_id).to_owned(),
                status: "PASS".to_owned(),
            })
            .collect(),
        overall_status: "PASS".to_owned(),
    })
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

pub fn run_preflight(source: &Path, output: &Path) -> Result<()> {
    if !build_is_admission_eligible(
        env!("P192_WCM_EMBEDDED_BUILD_PROFILE"),
        env!("P192_WCM_EMBEDDED_VERIFIER_DIRTY"),
    ) {
        return Err(
            "verifier preflight requires a clean release build with inspectable Git metadata"
                .to_owned(),
        );
    }
    if !normalized_absolute(source) || !normalized_absolute(output) {
        return Err("--source and --out must be normalized absolute paths".to_owned());
    }
    let canonical_source = fs::canonicalize(source)
        .map_err(|error| format!("canonicalize --source {}: {error}", source.display()))?;
    if canonical_source != source {
        return Err("--source may not traverse a symlink or alternate path spelling".to_owned());
    }
    if output != source.join("independent-verification.json") {
        return Err("--out must be exactly RUN_DIR/independent-verification.json".to_owned());
    }
    if output.exists() {
        return Err("independent-verification.json already exists".to_owned());
    }
    let verifier_binary_sha256 = running_verifier_digest()?;
    let result = verify_source(
        source,
        env!("P192_WCM_EMBEDDED_VERIFIER_COMMIT"),
        verifier_binary_sha256,
    )?;
    let value = serde_json::to_value(result).map_err(|error| error.to_string())?;
    let bytes = canonical_json(&value)?;
    let temporary = source.join(".independent-verification.json.tmp");
    let mut file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&temporary)
        .map_err(|error| format!("create {}: {error}", temporary.display()))?;
    file.write_all(&bytes)
        .map_err(|error| format!("write {}: {error}", temporary.display()))?;
    file.sync_all()
        .map_err(|error| format!("fsync {}: {error}", temporary.display()))?;
    drop(file);
    fs::rename(&temporary, output).map_err(|error| {
        format!(
            "rename {} to {}: {error}",
            temporary.display(),
            output.display()
        )
    })?;
    let parent = OpenOptions::new()
        .read(true)
        .open(source)
        .map_err(|error| format!("open source directory for fsync: {error}"))?;
    parent
        .sync_all()
        .map_err(|error| format!("fsync source directory: {error}"))?;
    Ok(())
}

fn build_is_admission_eligible(profile: &str, dirty: &str) -> bool {
    profile == "release" && dirty == "false"
}

pub fn expected_output_path(source: &Path) -> PathBuf {
    source.join("independent-verification.json")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_path_is_frozen() {
        let source = Path::new("/tmp/example-run");
        assert_eq!(
            expected_output_path(source),
            Path::new("/tmp/example-run/independent-verification.json")
        );
    }

    #[test]
    fn frozen_check_order_has_fifteen_unique_ids() {
        let unique: BTreeSet<_> = INDEPENDENT_CHECK_IDS.into_iter().collect();
        assert_eq!(unique.len(), 15);
        assert_eq!(INDEPENDENT_CHECK_IDS.len(), 15);
    }

    #[test]
    fn only_clean_release_builds_are_admission_eligible() {
        assert!(build_is_admission_eligible("release", "false"));
        assert!(!build_is_admission_eligible("debug", "false"));
        assert!(!build_is_admission_eligible("release", "true"));
        assert!(!build_is_admission_eligible("unknown", "true"));
        assert!(!commit_hex(&"0".repeat(40)));
        assert!(!commit_hex(&"A".repeat(40)));
        assert!(commit_hex(&"a".repeat(40)));
    }
}
