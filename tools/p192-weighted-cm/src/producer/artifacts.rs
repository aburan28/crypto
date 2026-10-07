//! Closed Role-1 artifact construction and file emission.

use std::collections::BTreeMap;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Component, Path};

use serde_json::{json, Value};

use super::controls;
use super::digest::{sha256, sha256_hex};
use super::evidence;
use super::factor_base;
use super::provenance::BuildProvenance;
use super::ref0::{self, ReferenceOutput, RECORD_COUNT, SHARD_COUNT};
use super::schema::{
    canonical_json_bytes, verdicts, CURVE_UID, EXPERIMENT_ID, PROTOCOL_VERSION,
    SCHEMA_REFERENCE_BOX,
};
use super::{PreflightArgs, Result};

const PRODUCER_PATHS: [&str; 9] = [
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

const WRITE_PATHS: [&str; 9] = [
    "identity.json",
    "maximal-order.json",
    "factor-base.json",
    "factor-base.sha256",
    "controls.json",
    "reference-box/candidate-records.bin",
    "reference-box/disposition.bin",
    "reference-box/certificates.bin",
    "reference-box.json",
];

fn artifact(path: &str, bytes: &[u8]) -> Value {
    json!({
        "byte_length": bytes.len() as u64,
        "path": path,
        "sha256": sha256_hex(bytes),
    })
}

fn reference_document(
    provenance: &BuildProvenance,
    factor_base_sha256: &str,
    output: &ReferenceOutput,
) -> Value {
    let candidate_hash = ref0::hex(&output.candidate_records_sha256);
    let candidate_shard = ref0::hex(&output.candidate_shard_sha256);
    let candidate_shell = ref0::hex(&output.candidate_shell_sha256);
    let disposition_hash = ref0::hex(&output.disposition_records_sha256);
    let disposition_shard = ref0::hex(&output.disposition_shard_sha256);
    let disposition_shell = ref0::hex(&output.disposition_shell_sha256);
    json!({
        "bounds": {
            "v_max": 2,
            "v_min": 1,
            "x_max_inclusive": 1024,
            "x_min": 0,
        },
        "candidate_count": RECORD_COUNT,
        "candidate_records": {
            "byte_length": output.candidate_records.len() as u64,
            "path": "reference-box/candidate-records.bin",
            "sha256": candidate_hash,
        },
        "candidate_shards": [{
            "candidate_shard_sha256": candidate_shard,
            "first_candidate_index": 0,
            "record_count": RECORD_COUNT,
            "records_sha256": candidate_hash,
            "shard_index": 0,
        }],
        "candidate_shell_sha256": candidate_shell,
        "certificates": {
            "byte_length": output.certificates.len() as u64,
            "path": "reference-box/certificates.bin",
            "record_count": output.status_counts.retained(),
            "sha256": sha256_hex(&output.certificates),
        },
        "complete": output.status_counts.invalid == 0
            && output.status_counts.unresolved == 0
            && output.status_counts.total() == RECORD_COUNT,
        "curve_uid": CURVE_UID,
        "disposition_records": {
            "byte_length": output.disposition_records.len() as u64,
            "path": "reference-box/disposition.bin",
            "sha256": disposition_hash,
        },
        "disposition_shards": [{
            "byte_length": output.disposition_records.len() as u64,
            "byte_offset": 0,
            "disposition_shard_sha256": disposition_shard,
            "first_candidate_index": 0,
            "record_count": RECORD_COUNT,
            "records_sha256": disposition_hash,
            "shard_index": 0,
        }],
        "disposition_shell_sha256": disposition_shell,
        "experiment_id": EXPERIMENT_ID,
        "factor_base_sha256": factor_base_sha256,
        "protocol_commit": provenance.protocol_commit,
        "protocol_version": PROTOCOL_VERSION,
        "schema": SCHEMA_REFERENCE_BOX,
        "shard_count": SHARD_COUNT,
        "shell_id": ref0::REFERENCE_SHELL,
        "source_commit": provenance.source_commit,
        "status_counts": {
            "complete": output.status_counts.complete,
            "invalid": output.status_counts.invalid,
            "nonprimitive_duplicate": output.status_counts.nonprimitive_duplicate,
            "one_large_prime": output.status_counts.one_large_prime,
            "rejected": output.status_counts.rejected,
            "two_large_prime": output.status_counts.two_large_prime,
            "unresolved": output.status_counts.unresolved,
        },
    })
}

fn digest_projection(domain: &[u8], value: &Value) -> Result<String> {
    let mut preimage = domain.to_vec();
    preimage.extend_from_slice(&canonical_json_bytes(value)?);
    Ok(sha256_hex(&preimage))
}

fn bindings(
    identity: &Value,
    maximal_order: &Value,
    controls: &Value,
    factor_base_sha256: &str,
    reference: &ReferenceOutput,
) -> Result<Value> {
    Ok(json!({
        "cm_order_tuple_sha256": digest_projection(
            b"P192-WCM-CM-ORDER-TUPLE-v1\0",
            &maximal_order["order"],
        )?,
        "control_results_sha256": controls["control_results_sha256"].clone(),
        "curve_tuple_sha256": digest_projection(
            b"P192-WCM-CURVE-TUPLE-v1\0",
            &identity["curve"],
        )?,
        "factor_base_sha256": factor_base_sha256,
        "reference_candidate_shard_sha256": ref0::hex(&reference.candidate_shard_sha256),
        "reference_candidate_shell_sha256": ref0::hex(&reference.candidate_shell_sha256),
        "reference_disposition_shard_sha256": ref0::hex(&reference.disposition_shard_sha256),
        "reference_disposition_shell_sha256": ref0::hex(&reference.disposition_shell_sha256),
    }))
}

fn verification_document(
    provenance: &BuildProvenance,
    producer_binary_sha256: &str,
    files: &BTreeMap<&'static str, Vec<u8>>,
    binding_values: Value,
) -> Result<Value> {
    let artifacts = PRODUCER_PATHS
        .iter()
        .map(|path| {
            files
                .get(path)
                .map(|bytes| artifact(path, bytes))
                .ok_or_else(|| format!("missing producer artifact bytes for {path}"))
        })
        .collect::<Result<Vec<_>>>()?;
    let checks = [
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
    Ok(json!({
        "artifacts": artifacts,
        "bindings": binding_values,
        "checks": verdicts(&checks),
        "experiment_id": EXPERIMENT_ID,
        "overall_status": "PASS",
        "producer_binary_sha256": producer_binary_sha256,
        "protocol_commit": provenance.protocol_commit,
        "protocol_version": PROTOCOL_VERSION,
        "role": "producer_self_check",
        "schema": "p192-wcm-verification-v1",
        "source_commit": provenance.source_commit,
    }))
}

fn validate_run_dir(path: &Path) -> Result<()> {
    let text = path
        .to_str()
        .ok_or_else(|| "RUN_DIR must be valid UTF-8".to_owned())?;
    let pieces = text.split('/').collect::<Vec<_>>();
    if pieces.iter().any(|piece| matches!(*piece, "." | ".."))
        || pieces
            .iter()
            .enumerate()
            .any(|(index, piece)| piece.is_empty() && index != 0)
    {
        return Err("RUN_DIR contains an empty or dot path segment".to_owned());
    }
    if path
        .components()
        .any(|component| matches!(component, Component::CurDir | Component::ParentDir))
    {
        return Err("RUN_DIR contains a dot path component".to_owned());
    }
    let metadata = fs::symlink_metadata(path)
        .map_err(|error| format!("inspect RUN_DIR {}: {error}", path.display()))?;
    if metadata.file_type().is_symlink() || !metadata.is_dir() {
        return Err("RUN_DIR must be an existing nonsymlink directory".to_owned());
    }
    for ancestor in path.ancestors() {
        let metadata = fs::symlink_metadata(ancestor)
            .map_err(|error| format!("inspect RUN_DIR ancestor {}: {error}", ancestor.display()))?;
        if metadata.file_type().is_symlink() {
            return Err(format!(
                "RUN_DIR ancestor is a symlink: {}",
                ancestor.display()
            ));
        }
    }
    if fs::read_dir(path)
        .map_err(|error| format!("read RUN_DIR {}: {error}", path.display()))?
        .next()
        .is_some()
    {
        return Err("RUN_DIR must be fresh and empty".to_owned());
    }
    Ok(())
}

fn sync_directory(path: &Path) -> Result<()> {
    fs::File::open(path)
        .and_then(|directory| directory.sync_all())
        .map_err(|error| format!("sync directory {}: {error}", path.display()))
}

fn write_new(path: &Path, bytes: &[u8]) -> Result<()> {
    let parent = path
        .parent()
        .ok_or_else(|| format!("artifact has no parent directory: {}", path.display()))?;
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| format!("artifact filename is not UTF-8: {}", path.display()))?;
    let temporary = parent.join(format!(".{name}.p192-wcm.tmp"));
    if path.exists() || temporary.exists() {
        return Err(format!(
            "refusing existing final or temporary artifact for {}",
            path.display()
        ));
    }
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .map_err(|error| format!("create {}: {error}", temporary.display()))?;
    file.write_all(bytes)
        .map_err(|error| format!("write {}: {error}", temporary.display()))?;
    file.sync_all()
        .map_err(|error| format!("sync {}: {error}", temporary.display()))?;
    drop(file);
    fs::rename(&temporary, path).map_err(|error| {
        format!(
            "atomically rename {} to {}: {error}",
            temporary.display(),
            path.display()
        )
    })?;
    sync_directory(parent)
}

pub fn write_preflight(args: &PreflightArgs, provenance: &BuildProvenance) -> Result<()> {
    // Construct and self-check every byte before changing the fresh RUN_DIR.
    let identity = evidence::identity(provenance)?;
    let maximal_order = evidence::maximal_order(provenance)?;
    let factor_base = factor_base::build(args.sieve_bound)?;
    let factor_base_document =
        factor_base::manifest(provenance.source_commit, args.sieve_bound, args.map_bound)?;
    let identity_bytes = canonical_json_bytes(&identity)?;
    let maximal_order_bytes = canonical_json_bytes(&maximal_order)?;
    let factor_base_bytes = canonical_json_bytes(&factor_base_document)?;
    let factor_base_digest = sha256(&factor_base_bytes);
    let factor_base_sha256 = ref0::hex(&factor_base_digest);
    let factor_base_sidecar = format!("{factor_base_sha256}  factor-base.json\n").into_bytes();
    let reference = ref0::build_reference(
        &factor_base,
        &factor_base_digest,
        args.sieve_bound,
        args.large_prime_bound,
    )?;
    if reference.status_counts.invalid != 0
        || reference.status_counts.unresolved != 0
        || reference.status_counts.total() != RECORD_COUNT
    {
        return Err("REF-0 is incomplete or invalid before artifact emission".to_owned());
    }
    let controls = controls::build(
        provenance,
        &factor_base_sha256,
        &factor_base_digest,
        &factor_base,
        &reference,
    )?;
    let controls_bytes = canonical_json_bytes(&controls)?;
    let reference_document = reference_document(provenance, &factor_base_sha256, &reference);
    let reference_bytes = canonical_json_bytes(&reference_document)?;

    let mut files = BTreeMap::new();
    files.insert("identity.json", identity_bytes);
    files.insert("maximal-order.json", maximal_order_bytes);
    files.insert("factor-base.json", factor_base_bytes);
    files.insert("factor-base.sha256", factor_base_sidecar);
    files.insert("controls.json", controls_bytes);
    files.insert("reference-box.json", reference_bytes);
    files.insert(
        "reference-box/candidate-records.bin",
        reference.candidate_records.clone(),
    );
    files.insert(
        "reference-box/disposition.bin",
        reference.disposition_records.clone(),
    );
    files.insert(
        "reference-box/certificates.bin",
        reference.certificates.clone(),
    );
    let executable =
        std::env::current_exe().map_err(|error| format!("resolve producer executable: {error}"))?;
    let executable_bytes = fs::read(&executable)
        .map_err(|error| format!("read producer executable {}: {error}", executable.display()))?;
    let producer_binary_sha256 = sha256_hex(&executable_bytes);
    let binding_values = bindings(
        &identity,
        &maximal_order,
        &controls,
        &factor_base_sha256,
        &reference,
    )?;
    let verification =
        verification_document(provenance, &producer_binary_sha256, &files, binding_values)?;
    let verification_bytes = canonical_json_bytes(&verification)?;

    validate_run_dir(&args.out)?;
    let reference_directory = args.out.join("reference-box");
    fs::create_dir(&reference_directory)
        .map_err(|error| format!("create {}: {error}", reference_directory.display()))?;
    sync_directory(&args.out)?;
    for path in WRITE_PATHS {
        write_new(
            &args.out.join(path),
            files
                .get(path)
                .ok_or_else(|| format!("internal missing bytes for {path}"))?,
        )?;
    }
    write_new(&args.out.join("verification.json"), &verification_bytes)?;
    sync_directory(&args.out)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn producer_inventory_excludes_verifier_and_supervisor_outputs() {
        assert_eq!(PRODUCER_PATHS.len(), 9);
        assert!(!PRODUCER_PATHS.contains(&"verification.json"));
        assert!(!PRODUCER_PATHS.contains(&"independent-verification.json"));
        assert!(!PRODUCER_PATHS.contains(&"independent-agreement.json"));
    }

    #[test]
    fn terminal_dot_segments_are_rejected_lexically() {
        assert!(validate_run_dir(Path::new("not-present/.")).is_err());
        assert!(validate_run_dir(Path::new("not-present/../elsewhere")).is_err());
    }

    #[test]
    fn writer_emits_only_the_ten_closed_producer_artifacts() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let run_dir = std::env::temp_dir().join(format!(
            "p192-wcm-artifact-test-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir(&run_dir).unwrap();
        let args = PreflightArgs {
            curve: "p192".to_owned(),
            sieve_bound: 65_521,
            map_bound: 113,
            large_prime_bound: 2_147_483_647,
            positive_control: "D=-23".to_owned(),
            reference_shell: "REF-0".to_owned(),
            reference_v_max: 2,
            reference_x_max_inclusive: 1024,
            out: run_dir.clone(),
        };
        let provenance = BuildProvenance {
            protocol_commit: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            source_commit: "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
            source_dirty: false,
            git_metadata_present: true,
            build_profile: "release",
        };
        write_preflight(&args, &provenance).unwrap();
        for path in PRODUCER_PATHS
            .iter()
            .copied()
            .chain(std::iter::once("verification.json"))
        {
            assert!(run_dir.join(path).is_file(), "missing {path}");
        }
        assert_eq!(
            fs::metadata(run_dir.join("reference-box/candidate-records.bin"))
                .unwrap()
                .len(),
            16_400
        );
        assert!(!run_dir.join("independent-verification.json").exists());
        assert!(!run_dir.join("independent-agreement.json").exists());
        assert!(!run_dir.join("independent-verifier-receipt.json").exists());
        for path in [
            "identity.json",
            "maximal-order.json",
            "factor-base.json",
            "controls.json",
            "reference-box.json",
            "verification.json",
        ] {
            let bytes = fs::read(run_dir.join(path)).unwrap();
            assert!(!bytes.ends_with(b"\n"));
            let value: Value = serde_json::from_slice(&bytes).unwrap();
            assert_eq!(canonical_json_bytes(&value).unwrap(), bytes);
        }
        fs::remove_dir_all(&run_dir).unwrap();
    }
}
