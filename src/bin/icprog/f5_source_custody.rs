//! Target-free F5 source custody from a validation-only n17 target build.
//! The archive contains the immutable build, never its disclosed-point config.
use super::{
    ordinary_build,
    sat_control::native,
    sat_target_build, target_build,
    target_control::{self, capsule},
    target_math,
};
use crypto_lib::cryptanalysis::{
    prepared_control_archive::verify_tar_inventory,
    prepared_sat_control::{canonical_sha, sha256},
};
use flate2::read::MultiGzDecoder;
use native::{load, read, require, save};
use serde_json::{json, Map, Value};
use std::{fs, io::Read, path::Path};

const QUESTION: &str = "f5-target-free-source-publication-n17-v1";
const STATUS: &str = "PASS_DATA_ONLY_F5_TARGET_FREE_SOURCE_PUBLICATION";
const CARD_SOURCE_QUESTION: &str = "fresh-paired-n17-source-publication-v1";
const CURVE_ID: &str = "EC1N17Ckb1hbbe2b5b6b1e6";
const SIDECARS: [&str; 3] = ["registration.json", "seal.json", "host-context.json"];

fn desc(bytes: &[u8]) -> Value {
    json!({"bytes":bytes.len(),"sha256":sha256(bytes)})
}
fn outside(root: &Path, publication: &Path, out: &Path) -> Result<(), String> {
    let parent = out
        .parent()
        .ok_or("F5 source output lacks parent")?
        .canonicalize()
        .map_err(|e| e.to_string())?;
    require(
        !parent.starts_with(root) && !parent.starts_with(publication),
        "F5 source result must stay outside original and publication",
    )
}
fn full_files(record: &capsule::Registration) -> Result<Map<String, Value>, String> {
    Ok(record
        .immutable_files
        .as_object()
        .ok_or("F5 source immutable inventory missing")?
        .iter()
        .map(|(name, value)| (format!("immutable/{name}"), value.clone()))
        .collect())
}
fn retained(root: &Path) -> Result<Map<String, Value>, String> {
    let mut files = Map::new();
    for name in SIDECARS {
        files.insert(
            name.into(),
            desc(&read(&root.join(name), 16 * 1024 * 1024)?),
        );
    }
    for (name, value) in native::inventory(&root.join("immutable/build-receipts"))?
        .as_object()
        .ok_or("F5 source build receipts missing")?
    {
        files.insert(format!("immutable/build-receipts/{name}"), value.clone());
    }
    Ok(files)
}
fn receipt_subset(record: &capsule::Registration) -> Result<Value, String> {
    Ok(json!(record
        .immutable_files
        .as_object()
        .ok_or("F5 source immutable inventory missing")?
        .iter()
        .filter_map(|(name, value)| name
            .strip_prefix("build-receipts/")
            .map(|relative| (relative.to_string(), value.clone())))
        .collect::<Map<String, Value>>()))
}
fn source_subset(record: &capsule::Registration) -> Result<Value, String> {
    Ok(json!(record
        .immutable_files
        .as_object()
        .ok_or("F5 source immutable inventory missing")?
        .iter()
        .filter_map(|(name, value)| name
            .strip_prefix("source/")
            .map(|relative| (relative.to_string(), value.clone())))
        .collect::<Map<String, Value>>()))
}
fn registration_sidecar(root: &Path, seal: &str) -> Result<capsule::Registration, String> {
    let bytes = read(&root.join("registration.json"), 16 * 1024 * 1024)?;
    let record: capsule::Registration =
        serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    require(
        record.validation_only
            && canonical_sha(&json!(record))? == seal
            && load(&root.join("seal.json"))?["registration_sha256"] == seal,
        "F5 source sidecar is not the sealed validation-only registration",
    )?;
    Ok(record)
}
fn claims(header: &Value) -> Result<(), String> {
    require(
        header["schema_version"] == 1
            && header["question"] == QUESTION
            && header["stage"] == "validation-source-only-before-fresh-card"
            && header["target_config_archived"] == false
            && header["target_card_created"] == false
            && header["target_worker_calls"] == 0
            && header["source_bound_execution_admitted"] == false
            && header["fresh_paired_qualification"] == false
            && header["online_wall_ns"].is_null()
            && header["online_speedup"].is_null(),
        "F5 source publication overstates target-free custody",
    )
}
fn frozen_self(root: &Path, record: &capsule::Registration) -> Result<String, String> {
    capsule::require_identity(&capsule::identity(), &record.worker_build_identity)?;
    let own = std::env::current_exe().map_err(|e| e.to_string())?;
    let hash = sha256(&read(&own, 128 * 1024 * 1024)?);
    require(
        own.canonicalize().map_err(|e| e.to_string())?
            == root
                .join("immutable/bin/icprog")
                .canonicalize()
                .map_err(|e| e.to_string())?
            && hash == record.auditor_sha256,
        "F5 source adoption requires its original frozen controller",
    )?;
    Ok(hash)
}
fn source_descriptor(path: &Path, source_seal: &str, source: &Value) -> Result<Vec<u8>, String> {
    let bytes = read(path, 65536)?;
    let value = target_math::parse(&bytes)?;
    let fields = [
        "schema_version",
        "question",
        "role",
        "curve_id",
        "registration_sha256",
        "source_manifest_sha256",
        "archive_sha256",
        "target_free",
        "source_bound_execution_admitted",
    ];
    require(
        value.as_object().is_some_and(|obj| {
            obj.len() == fields.len() && fields.iter().all(|&field| obj.contains_key(field))
        }) && value["schema_version"] == 1
            && value["question"] == CARD_SOURCE_QUESTION
            && value["role"] == "f5"
            && value["curve_id"] == CURVE_ID
            && value["registration_sha256"] == source_seal
            && value["source_manifest_sha256"] == source["source_manifest_sha256"]
            && value["archive_sha256"] == source["archive_sha256"]
            && value["target_free"] == true
            && value["source_bound_execution_admitted"] == false,
        "F5 source descriptor differs from complete target-free publication",
    )?;
    Ok(bytes)
}

pub(super) fn publish(
    capsule_root: &Path,
    publication: &Path,
    validation_seal: &str,
    out: &Path,
) -> Result<String, String> {
    let root = capsule_root.canonicalize().map_err(|e| e.to_string())?;
    require(
        !root
            .join("consumed.json")
            .try_exists()
            .map_err(|e| e.to_string())?,
        "F5 source validation capsule was consumed",
    )?;
    let record = capsule::check_capsule(&root, validation_seal)?;
    require(
        record.validation_only,
        "F5 source requires validation-only build",
    )?;
    target_build::verify(&root, &record)?;
    target_control::original_preparation(&record.preparation)?;
    let publication_parent = publication
        .parent()
        .ok_or("F5 publication lacks parent")?
        .canonicalize()
        .map_err(|e| e.to_string())?;
    require(
        !publication_parent.starts_with(&root),
        "F5 source publication must stay outside its original capsule",
    )?;
    fs::create_dir(publication).map_err(|e| e.to_string())?;
    let publication = publication.canonicalize().map_err(|e| e.to_string())?;
    outside(&root, &publication, out)?;
    let retained = retained(&root)?;
    for name in retained.keys() {
        let target = publication.join(name);
        fs::create_dir_all(target.parent().ok_or("F5 retained parent missing")?)
            .map_err(|e| e.to_string())?;
        native::create(&target, &read(&root.join(name), 16 * 1024 * 1024)?)?;
    }
    let files = full_files(&record)?;
    ordinary_build::archive(&root, &publication.join("capsule.tar.gz"), &files)?;
    let archive = read(&publication.join("capsule.tar.gz"), 96 * 1024 * 1024)?;
    let header = json!({"schema_version":1,"question":QUESTION,
        "stage":"validation-source-only-before-fresh-card",
        "original_capsule_root":root,"validation_registration_sha256":validation_seal,
        "source_commit":record.source_commit,
        "source_manifest_sha256":record.source_manifest_sha256,
        "worker_sha256":record.worker_sha256,"auditor_sha256":record.auditor_sha256,
        "worker_build_identity":record.worker_build_identity,
        "immutable_files":files,"sidecars":retained,
        "archive":{"path":"capsule.tar.gz","bytes":archive.len(),"sha256":sha256(&archive)},
        "target_config_archived":false,"target_card_created":false,"target_worker_calls":0,
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,
        "online_wall_ns":null,"online_speedup":null});
    claims(&header)?;
    save(&publication.join("SOURCE_PUBLICATION.json"), &header)?;
    let source_seal = canonical_sha(&header)?;
    save(
        &publication.join("SOURCE_SEAL.json"),
        &json!({"source_registration_sha256":source_seal}),
    )?;
    replay(&publication, &source_seal, out)
}

fn verify(publication: &Path, source_seal: &str) -> Result<Value, String> {
    let header = load(&publication.join("SOURCE_PUBLICATION.json"))?;
    claims(&header)?;
    require(
        canonical_sha(&header)? == source_seal
            && load(&publication.join("SOURCE_SEAL.json"))?["source_registration_sha256"]
                == source_seal,
        "F5 external source publication seal differs",
    )?;
    let validation_seal = header["validation_registration_sha256"]
        .as_str()
        .ok_or("F5 original validation seal missing")?;
    let record = registration_sidecar(publication, validation_seal)?;
    let original = Path::new(
        header["original_capsule_root"]
            .as_str()
            .ok_or("F5 original capsule root missing")?,
    );
    require(
        original.is_absolute()
            && header["source_commit"] == record.source_commit
            && header["source_manifest_sha256"] == record.source_manifest_sha256
            && header["worker_sha256"] == record.worker_sha256
            && header["auditor_sha256"] == record.auditor_sha256
            && header["worker_build_identity"] == record.worker_build_identity
            && header["immutable_files"] == json!(full_files(&record)?)
            && header["sidecars"] == json!(retained(publication)?),
        "F5 source publication differs from original sealed build",
    )?;
    require(
        receipt_subset(&record)?
            == native::inventory(&publication.join("immutable/build-receipts"))?
            && canonical_sha(&source_subset(&record)?)? == record.source_manifest_sha256
            && record.immutable_files["bin/prepared_target_worker"]["sha256"]
                == record.worker_sha256
            && record.immutable_files["bin/icprog"]["sha256"] == record.auditor_sha256,
        "F5 source build receipt or executable inventory differs",
    )?;
    target_build::verify_from(publication, &record, original)?;
    require(
        header["archive"]["path"] == "capsule.tar.gz",
        "F5 source archive path differs",
    )?;
    let archive = read(&publication.join("capsule.tar.gz"), 96 * 1024 * 1024)?;
    require(
        desc(&archive)
            == json!({"bytes":header["archive"]["bytes"],"sha256":header["archive"]["sha256"]}),
        "F5 source archive bytes differ",
    )?;
    let mut tar = Vec::new();
    MultiGzDecoder::new(archive.as_slice())
        .take(512 * 1024 * 1024 + 1)
        .read_to_end(&mut tar)
        .map_err(|e| e.to_string())?;
    require(
        tar.len() <= 512 * 1024 * 1024,
        "F5 source archive exceeds decompression bound",
    )?;
    let members = verify_tar_inventory(&tar, &header["immutable_files"])?;
    require(
        header["immutable_files"]
            .as_object()
            .is_some_and(|files| files.keys().all(|name| name.starts_with("immutable/"))),
        "F5 source archive includes a target-bearing sidecar",
    )?;
    Ok(json!({"schema_version":1,"status":STATUS,
        "source_registration_sha256":source_seal,
        "validation_registration_sha256":validation_seal,
        "source_commit":record.source_commit,
        "source_manifest_sha256":record.source_manifest_sha256,
        "worker_sha256":record.worker_sha256,"auditor_sha256":record.auditor_sha256,
        "archive_sha256":sha256(&archive),"archive_regular_files":members,
        "target_config_archived":false,"target_card_created":false,
        "archived_binaries_executed":0,"target_worker_calls":0,
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,
        "online_wall_ns":null,"online_speedup":null}))
}
pub(super) fn replay(publication: &Path, source_seal: &str, out: &Path) -> Result<String, String> {
    let publication = publication.canonicalize().map_err(|e| e.to_string())?;
    outside(&publication, &publication, out)?;
    let result = verify(&publication, source_seal)?;
    save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

/// Bind one later scalar-free card to the exact prepublished F5 executable.
/// Once adoption-start is written, an interruption is terminal for this build.
pub(super) fn adopt_card(
    root: &Path,
    publication: &Path,
    source_seal: &str,
    descriptor_path: &Path,
    card_path: &Path,
    validation_seal: &str,
    out: &Path,
) -> Result<String, String> {
    native::enforce_hardware()?;
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let publication = publication.canonicalize().map_err(|e| e.to_string())?;
    let descriptor_path = descriptor_path.canonicalize().map_err(|e| e.to_string())?;
    let card_path = card_path.canonicalize().map_err(|e| e.to_string())?;
    outside(&root, &publication, out)?;
    require(
        !descriptor_path.starts_with(&root)
            && !card_path.starts_with(&root)
            && !root
                .join("adoption-start.json")
                .try_exists()
                .map_err(|e| e.to_string())?
            && !root
                .join("consumed.json")
                .try_exists()
                .map_err(|e| e.to_string())?,
        "F5 source capsule was already adopted/consumed or card is internal",
    )?;
    let record = capsule::check_capsule(&root, validation_seal)?;
    require(
        record.validation_only,
        "F5 card adoption requires validation-only build",
    )?;
    let own = frozen_self(&root, &record)?;
    target_build::verify(&root, &record)?;
    target_control::original_preparation(&record.preparation)?;
    let source = verify(&publication, source_seal)?;
    let header = load(&publication.join("SOURCE_PUBLICATION.json"))?;
    require(
        header["original_capsule_root"] == json!(root)
            && source["validation_registration_sha256"] == validation_seal
            && source["source_manifest_sha256"] == record.source_manifest_sha256
            && source["worker_sha256"] == record.worker_sha256
            && source["auditor_sha256"] == record.auditor_sha256
            && header["immutable_files"] == json!(full_files(&record)?),
        "F5 card adoption differs from original published source/build",
    )?;
    let descriptor = source_descriptor(&descriptor_path, source_seal, &source)?;
    let descriptor_sha = sha256(&descriptor);
    let (card, card_sha) = sat_target_build::capsule::card(&card_path)?;
    require(
        card.source_publications["f5"] == descriptor_sha,
        "F5 public point card does not bind the original source descriptor",
    )?;
    let mut cfg = capsule::config(&root)?;
    require(
        cfg.target != card.target,
        "F5 fresh card repeats the disclosed validation target",
    )?;
    cfg.target = card.target;
    cfg.validate()?;
    let old_config_sha = sha256(&read(&root.join("config.json"), 65536)?);
    let old_record_sha = sha256(&read(&root.join("registration.json"), 16 * 1024 * 1024)?);
    let card_bytes = read(&card_path, 65536)?;
    require(
        sha256(&card_bytes) == card_sha
            && read(&descriptor_path, 65536)? == descriptor
            && frozen_self(&root, &record)? == own,
        "F5 card/source/controller changed before irreversible adoption",
    )?;
    save(
        &root.join("adoption-start.json"),
        &json!({
        "schema_version":1,"status":"ONE_USE_F5_CARD_ADOPTION_STARTED",
        "validation_registration_sha256":validation_seal,
        "source_registration_sha256":source_seal,"source_publication_root":publication,
        "source_descriptor_path":descriptor_path,"source_descriptor_sha256":descriptor_sha,
        "card_path":card_path,"card_sha256":card_sha,"target":card.target,
        "old_config_sha256":old_config_sha,"old_registration_file_sha256":old_record_sha,
        "controller_sha256":own,"source_bound_execution_admitted":false,
        "target_worker_calls":0,"online_speedup":null}),
    )?;
    let prior = root.join("precard-validation");
    fs::create_dir(&prior).map_err(|e| e.to_string())?;
    for name in ["config.json", "registration.json", "seal.json"] {
        fs::rename(root.join(name), prior.join(name)).map_err(|e| e.to_string())?;
    }
    native::create(&root.join("source-descriptor.json"), &descriptor)?;
    native::create(&root.join("card.json"), &card_bytes)?;
    save(&root.join("config.json"), &json!(cfg))?;
    let mut scientific = record.clone();
    scientific.validation_only = false;
    scientific.config_sha256 = sha256(&read(&root.join("config.json"), 65536)?);
    let scientific_seal = canonical_sha(&json!(scientific))?;
    save(&root.join("registration.json"), &json!(scientific))?;
    save(
        &root.join("seal.json"),
        &json!({"registration_sha256":scientific_seal}),
    )?;
    let checked = capsule::check_capsule(&root, &scientific_seal)?;
    target_build::verify(&root, &checked)?;
    save(
        &root.join("adoption-terminal.json"),
        &json!({
        "schema_version":1,"status":"F5_CARD_REGISTERED_NOT_DISPATCHED",
        "validation_registration_sha256":validation_seal,
        "registration_sha256":scientific_seal,
        "source_registration_sha256":source_seal,
        "source_manifest_sha256":checked.source_manifest_sha256,
        "worker_sha256":checked.worker_sha256,"auditor_sha256":checked.auditor_sha256,
        "source_descriptor_sha256":descriptor_sha,"card_sha256":card_sha,
        "target":card.target,"controller_sha256":own,
        "target_worker_calls":0,"source_bound_execution_admitted":false,
        "online_speedup":null}),
    )?;
    let result = verify_adoption(&root, &scientific_seal)?;
    save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

pub(super) fn verify_adoption(root: &Path, scientific_seal: &str) -> Result<Value, String> {
    let record = capsule::check_capsule(root, scientific_seal)?;
    require(
        !record.validation_only,
        "F5 card was not scientifically registered",
    )?;
    target_build::verify(root, &record)?;
    let own = frozen_self(root, &record)?;
    let start = load(&root.join("adoption-start.json"))?;
    let terminal = load(&root.join("adoption-terminal.json"))?;
    let source_root = Path::new(
        start["source_publication_root"]
            .as_str()
            .ok_or("F5 original source publication missing")?,
    );
    let source_seal = start["source_registration_sha256"]
        .as_str()
        .ok_or("F5 source seal missing")?;
    let source = verify(source_root, source_seal)?;
    let header = load(&source_root.join("SOURCE_PUBLICATION.json"))?;
    let validation_seal = start["validation_registration_sha256"]
        .as_str()
        .ok_or("F5 prior validation seal missing")?;
    let prior = root.join("precard-validation");
    let prior_record = registration_sidecar(&prior, validation_seal)?;
    require(
        header["original_capsule_root"] == json!(root)
            && source["validation_registration_sha256"] == validation_seal
            && sha256(&read(&prior.join("registration.json"), 16 * 1024 * 1024)?)
                == start["old_registration_file_sha256"]
            && sha256(&read(&prior.join("config.json"), 65536)?) == start["old_config_sha256"]
            && prior_record.config_sha256 == start["old_config_sha256"]
            && prior_record.immutable_files == record.immutable_files
            && prior_record.source_commit == record.source_commit
            && prior_record.source_manifest_sha256 == record.source_manifest_sha256
            && prior_record.worker_sha256 == record.worker_sha256
            && prior_record.auditor_sha256 == record.auditor_sha256
            && prior_record.worker_build_identity == record.worker_build_identity
            && json!(prior_record.preparation) == json!(record.preparation),
        "F5 post-card registration differs from original pre-card source",
    )?;
    let descriptor_path = Path::new(
        start["source_descriptor_path"]
            .as_str()
            .ok_or("F5 original descriptor path missing")?,
    );
    let descriptor = source_descriptor(descriptor_path, source_seal, &source)?;
    let card_path = Path::new(
        start["card_path"]
            .as_str()
            .ok_or("F5 original card path missing")?,
    );
    let (card, card_sha) = sat_target_build::capsule::card(card_path)?;
    let descriptor_sha = sha256(&descriptor);
    let cfg = capsule::config(root)?;
    require(
        start["schema_version"] == 1
            && start["status"] == "ONE_USE_F5_CARD_ADOPTION_STARTED"
            && start["controller_sha256"] == own
            && start["source_bound_execution_admitted"] == false
            && start["target_worker_calls"] == 0
            && start["online_speedup"].is_null()
            && start["source_descriptor_sha256"] == descriptor_sha
            && start["card_sha256"] == card_sha
            && start["target"] == json!(card.target)
            && read(&root.join("source-descriptor.json"), 65536)? == descriptor
            && read(&root.join("card.json"), 65536)? == read(card_path, 65536)?
            && card.source_publications["f5"] == descriptor_sha
            && cfg.target == card.target
            && record.config_sha256 == sha256(&read(&root.join("config.json"), 65536)?)
            && terminal
                == json!({
                "schema_version":1,"status":"F5_CARD_REGISTERED_NOT_DISPATCHED",
                "validation_registration_sha256":validation_seal,
                "registration_sha256":scientific_seal,
                "source_registration_sha256":source_seal,
                "source_manifest_sha256":record.source_manifest_sha256,
                "worker_sha256":record.worker_sha256,
                "auditor_sha256":record.auditor_sha256,
                "source_descriptor_sha256":descriptor_sha,"card_sha256":card_sha,
                "target":card.target,"controller_sha256":own,
                "target_worker_calls":0,"source_bound_execution_admitted":false,
                "online_speedup":null}),
        "F5 original source/card adoption record differs",
    )?;
    Ok(
        json!({"schema_version":1,"status":"PASS_SOURCE_BOUND_F5_CARD_ADOPTION",
        "registration_sha256":scientific_seal,
        "validation_registration_sha256":validation_seal,
        "source_registration_sha256":source_seal,
        "source_manifest_sha256":record.source_manifest_sha256,
        "source_archive_sha256":source["archive_sha256"],
        "worker_sha256":record.worker_sha256,"controller_sha256":own,
        "source_descriptor_sha256":descriptor_sha,"card_sha256":card_sha,
        "target":card.target,"target_worker_calls":0,
        "source_bound_execution_admitted":false,"online_speedup":null}),
    )
}

pub(super) fn audit_adoption(
    root: &Path,
    scientific_seal: &str,
    out: &Path,
) -> Result<String, String> {
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    outside(&root, &root, out)?;
    require(
        !root
            .join("consumed.json")
            .try_exists()
            .map_err(|e| e.to_string())?,
        "F5 standalone adoption audit must precede target dispatch",
    )?;
    let result = verify_adoption(&root, scientific_seal)?;
    save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn source_claims_reject_target_or_dispatch() {
        let mut header = json!({"schema_version":1,"question":QUESTION,
            "stage":"validation-source-only-before-fresh-card",
            "target_config_archived":false,"target_card_created":false,"target_worker_calls":0,
            "source_bound_execution_admitted":false,"fresh_paired_qualification":false,
            "online_wall_ns":null,"online_speedup":null});
        claims(&header).unwrap();
        header["target_config_archived"] = json!(true);
        assert!(claims(&header).is_err());
        header["target_config_archived"] = json!(false);
        header["target_worker_calls"] = json!(1);
        assert!(claims(&header).is_err());
    }
}
