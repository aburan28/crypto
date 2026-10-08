//! Portable, data-only custody for a validation-only n17 target build.
//! Archived programs are never extracted or executed by this module.
use super::{ordinary_build, sat_control::native, target_build, target_control::capsule};
use crypto_lib::cryptanalysis::{
    prepared_control_archive::verify_tar_inventory,
    prepared_sat_control::{canonical_sha, sha256},
};
use flate2::read::MultiGzDecoder;
use native::{load, read, require, save};
use serde_json::{json, Value};
use std::{fs, io::Read, path::Path};

const QUESTION: &str = "target-build-custody-n17-v1";
const FILES: [&str; 4] = [
    "config.json",
    "host-context.json",
    "registration.json",
    "seal.json",
];

fn desc(bytes: &[u8]) -> Value {
    json!({"bytes":bytes.len(),"sha256":sha256(bytes)})
}
fn immutable_subset(record: &capsule::Registration, prefix: &str) -> Result<Value, String> {
    Ok(json!(record
        .immutable_files
        .as_object()
        .ok_or("missing target inventory")?
        .iter()
        .filter_map(|(name, value)| name
            .strip_prefix(prefix)
            .map(|relative| (relative.to_string(), value.clone())))
        .collect::<serde_json::Map<String, Value>>()))
}
fn unconsumed(root: &Path) -> Result<(), String> {
    require(
        !root
            .join("consumed.json")
            .try_exists()
            .map_err(|e| e.to_string())?,
        "target validation capsule was consumed",
    )
}
fn validation(record: &capsule::Registration) -> Result<(), String> {
    require(
        record.validation_only,
        "target custody requires validation-only registration",
    )
}
fn sidecars(root: &Path) -> Result<serde_json::Map<String, Value>, String> {
    let mut files = serde_json::Map::new();
    for name in FILES {
        files.insert(
            name.into(),
            desc(&read(&root.join(name), 16 * 1024 * 1024)?),
        );
    }
    for (name, value) in native::inventory(&root.join("immutable/build-receipts"))?
        .as_object()
        .ok_or("missing target build receipts")?
    {
        files.insert(format!("immutable/build-receipts/{name}"), value.clone());
    }
    Ok(files)
}
fn full_files(
    root: &Path,
    record: &capsule::Registration,
) -> Result<serde_json::Map<String, Value>, String> {
    let mut files = serde_json::Map::new();
    for name in FILES {
        files.insert(
            name.into(),
            desc(&read(&root.join(name), 16 * 1024 * 1024)?),
        );
    }
    for (name, value) in record
        .immutable_files
        .as_object()
        .ok_or("missing target inventory")?
    {
        files.insert(format!("immutable/{name}"), value.clone());
    }
    Ok(files)
}
fn claims(header: &Value) -> Result<(), String> {
    require(
        header["schema_version"] == 1
            && header["question"] == QUESTION
            && header["stage"] == "validation-only-not-dispatched"
            && header["scientific_worker_calls"] == 0
            && header["execution_admitted"] == false
            && header["source_bound_execution_admitted"] == false
            && header["full_goal_complete"] == false
            && header["online_wall_ns"].is_null()
            && header["online_speedup"].is_null(),
        "target custody scope or claims differ",
    )
}

pub(super) fn publish(
    root: &Path,
    publication: &Path,
    expected: &str,
    out: &Path,
) -> Result<String, String> {
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let record = capsule::check_capsule(&root, expected)?;
    validation(&record)?;
    unconsumed(&root)?;
    target_build::verify(&root, &record)?;
    fs::create_dir(publication).map_err(|e| e.to_string())?;
    let retained = sidecars(&root)?;
    for name in retained.keys() {
        let target = publication.join(name);
        fs::create_dir_all(target.parent().ok_or("missing target sidecar parent")?)
            .map_err(|e| e.to_string())?;
        native::create(&target, &read(&root.join(name), 16 * 1024 * 1024)?)?;
    }
    ordinary_build::archive(
        &root,
        &publication.join("capsule.tar.gz"),
        &full_files(&root, &record)?,
    )?;
    let archive = read(&publication.join("capsule.tar.gz"), 96 * 1024 * 1024)?;
    let header = json!({"schema_version":1,"question":QUESTION,"stage":"validation-only-not-dispatched",
        "registration_sha256":expected,"capsule_root":root,"source_commit":record.source_commit,
        "immutable_file_count":record.immutable_files.as_object().ok_or("missing inventory")?.len(),
        "files":retained,"archive":{"path":"capsule.tar.gz","bytes":archive.len(),"sha256":sha256(&archive)},
        "scientific_worker_calls":0,"execution_admitted":false,"source_bound_execution_admitted":false,
        "full_goal_complete":false,"online_wall_ns":null,"online_speedup":null});
    save(&publication.join("PUBLICATION.json"), &header)?;
    capsule::check_capsule(&root, expected)?;
    unconsumed(&root)?;
    replay(publication, expected, out)
}
fn verify(publication: &Path, expected: &str) -> Result<Value, String> {
    let header = load(&publication.join("PUBLICATION.json"))?;
    claims(&header)?;
    let record = capsule::registration(publication, expected)?;
    validation(&record)?;
    unconsumed(publication)?;
    require(
        header["registration_sha256"] == expected
            && header["source_commit"] == record.source_commit
            && header["immutable_file_count"]
                == record
                    .immutable_files
                    .as_object()
                    .ok_or("missing inventory")?
                    .len(),
        "target custody descriptor differs",
    )?;
    let original = Path::new(
        header["capsule_root"]
            .as_str()
            .ok_or("missing original target capsule root")?,
    );
    require(
        original.is_absolute(),
        "original target capsule path is not absolute",
    )?;
    let retained = sidecars(publication)?;
    require(
        header["files"] == json!(retained),
        "target custody sidecars changed",
    )?;
    require(
        immutable_subset(&record, "build-receipts/")?
            == native::inventory(&publication.join("immutable/build-receipts"))?,
        "target receipt set differs from original sealed inventory",
    )?;
    require(
        canonical_sha(&immutable_subset(&record, "source/")?)? == record.source_manifest_sha256
            && record.immutable_files[format!("bin/{}", capsule::WORKER)]["sha256"]
                == record.worker_sha256
            && record.immutable_files["bin/icprog"]["sha256"] == record.auditor_sha256,
        "archived target source or binary pins differ",
    )?;
    target_build::verify_from(publication, &record, original)?;
    require(
        header["archive"]["path"] == "capsule.tar.gz",
        "target custody archive path differs",
    )?;
    let archive = read(&publication.join("capsule.tar.gz"), 96 * 1024 * 1024)?;
    require(
        desc(&archive)
            == json!({"bytes":header["archive"]["bytes"],"sha256":header["archive"]["sha256"]}),
        "target custody archive bytes differ",
    )?;
    let mut tar = Vec::new();
    MultiGzDecoder::new(archive.as_slice())
        .take(512 * 1024 * 1024 + 1)
        .read_to_end(&mut tar)
        .map_err(|e| e.to_string())?;
    require(
        tar.len() <= 512 * 1024 * 1024,
        "target custody decompression exceeded bound",
    )?;
    let members = verify_tar_inventory(&tar, &json!(full_files(publication, &record)?))?;
    Ok(
        json!({"schema_version":1,"status":"PASS_DATA_ONLY_TARGET_BUILD_CUSTODY",
        "registration_sha256":expected,"source_commit":record.source_commit,
        "source_manifest_sha256":record.source_manifest_sha256,"archive_sha256":sha256(&archive),
        "archive_regular_files":members,"complete_published_capsule_bytes_verified":true,
        "archived_binaries_executed":0,"scientific_worker_calls":0,"execution_admitted":false,
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,
        "headline_eligible":false,"promotion_eligible":false,"full_goal_complete":false,
        "online_wall_ns":null,"online_speedup":null}),
    )
}
pub(super) fn replay(publication: &Path, expected: &str, out: &Path) -> Result<String, String> {
    let result = verify(publication, expected)?;
    save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn custody_claims_cannot_promote_a_validation_build() {
        let good = json!({"schema_version":1,"question":QUESTION,"stage":"validation-only-not-dispatched",
            "scientific_worker_calls":0,"execution_admitted":false,"source_bound_execution_admitted":false,
            "full_goal_complete":false,"online_wall_ns":null,"online_speedup":null});
        claims(&good).unwrap();
        for (key, value) in [
            ("execution_admitted", json!(true)),
            ("source_bound_execution_admitted", json!(true)),
            ("scientific_worker_calls", json!(1)),
            ("online_wall_ns", json!(0)),
            ("question", json!("scientific-qualification")),
        ] {
            let mut wrong = good.clone();
            wrong[key] = value;
            assert!(claims(&wrong).is_err(), "{key}");
        }
    }
    #[test]
    fn registered_sidecars_reject_mutated_host_context() {
        let fixture = Path::new("research/ic_candidate_tournament_20260915/goal_20260924/native-target-controller-v1/control-v2");
        let seal = "a69894962cfd5d8c71014beee9747f58de8a647abcf4591cbdae4d470a0bba68";
        let root = std::env::temp_dir().join(format!("target-sidecars-{}", std::process::id()));
        fs::create_dir(&root).unwrap();
        for name in FILES {
            fs::copy(fixture.join(name), root.join(name)).unwrap();
        }
        assert!(capsule::registration(&root, seal).is_ok());
        fs::write(root.join("host-context.json"), b"{}").unwrap();
        assert!(capsule::registration(&root, seal).is_err());
        fs::remove_dir_all(root).unwrap();
    }
}
