//! Data-only custody for a target-free SAT source/build capsule.
//! The replay path never extracts or executes an archived program.
use super::{ordinary_build, sat_control::native, sat_target_build, sat_target_build::capsule};
use crypto_lib::cryptanalysis::{
    prepared_control_archive::verify_tar_inventory,
    prepared_sat_control::{canonical_sha, sha256},
};
use flate2::read::MultiGzDecoder;
use native::{load, read, require, save};
use serde_json::{json, Map, Value};
use std::{fs, io::Read, path::Path};

const BUILD_QUESTION: &str = "sat-target-build-custody-n17-v1";
const SCIENTIFIC_QUESTION: &str = "sat-target-scientific-registration-custody-n17-v1";
const SIDECARS: [&str; 4] = [
    "config.json",
    "host-context.json",
    "registration.json",
    "seal.json",
];

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Kind {
    Validation,
    Scientific,
}
impl Kind {
    fn question(self) -> &'static str {
        match self {
            Self::Validation => BUILD_QUESTION,
            Self::Scientific => SCIENTIFIC_QUESTION,
        }
    }
    fn status(self) -> &'static str {
        match self {
            Self::Validation => "PASS_DATA_ONLY_SAT_TARGET_BUILD_CUSTODY",
            Self::Scientific => "PASS_DATA_ONLY_SAT_TARGET_SCIENTIFIC_REGISTRATION_CUSTODY",
        }
    }
    fn stage(self) -> &'static str {
        match self {
            Self::Validation => "validation-only-not-dispatched",
            Self::Scientific => "scientific-registered-not-dispatched",
        }
    }
    fn check(self, record: &capsule::Registration) -> Result<(), String> {
        require(
            record.validation_only == (self == Self::Validation),
            "SAT publication kind differs from sealed registration",
        )
    }
}

fn desc(bytes: &[u8]) -> Value {
    json!({"bytes":bytes.len(),"sha256":sha256(bytes)})
}
fn unconsumed(root: &Path) -> Result<(), String> {
    require(
        !root
            .join("consumed.json")
            .try_exists()
            .map_err(|e| e.to_string())?,
        "SAT target registration was consumed",
    )
}
fn sidecars(root: &Path) -> Result<Map<String, Value>, String> {
    let mut files = Map::new();
    for name in SIDECARS {
        files.insert(
            name.into(),
            desc(&read(&root.join(name), 16 * 1024 * 1024)?),
        );
    }
    for (name, value) in native::inventory(&root.join("immutable/build-receipts"))?
        .as_object()
        .ok_or("SAT build receipt inventory missing")?
    {
        files.insert(format!("immutable/build-receipts/{name}"), value.clone());
    }
    for name in [
        "marked-audit.json",
        "preparation-audit.json",
        "exporter-parity-audit.json",
    ] {
        let path = format!("immutable/assets/controls/{name}");
        if root.join(&path).try_exists().map_err(|e| e.to_string())? {
            files.insert(
                path.clone(),
                desc(&read(&root.join(path), 16 * 1024 * 1024)?),
            );
        }
    }
    Ok(files)
}
fn full_files(root: &Path, record: &capsule::Registration) -> Result<Map<String, Value>, String> {
    let mut files = Map::new();
    for name in SIDECARS {
        files.insert(
            name.into(),
            desc(&read(&root.join(name), 16 * 1024 * 1024)?),
        );
    }
    for (name, value) in record
        .immutable_files
        .as_object()
        .ok_or("SAT immutable inventory missing")?
    {
        files.insert(format!("immutable/{name}"), value.clone());
    }
    Ok(files)
}
fn subset(record: &capsule::Registration, prefix: &str) -> Result<Value, String> {
    Ok(json!(record
        .immutable_files
        .as_object()
        .ok_or("SAT immutable inventory missing")?
        .iter()
        .filter_map(|(name, value)| name
            .strip_prefix(prefix)
            .map(|relative| (relative.to_string(), value.clone())))
        .collect::<Map<String, Value>>()))
}
fn claims(header: &Value, kind: Kind) -> Result<(), String> {
    require(
        header["schema_version"] == 1
            && header["question"] == kind.question()
            && header["stage"] == kind.stage()
            && header["scientific_worker_calls"] == 0
            && header["execution_admitted"] == false
            && header["source_bound_execution_admitted"] == false
            && header["fresh_paired_qualification"] == false
            && header["online_wall_ns"].is_null()
            && header["online_speedup"].is_null(),
        "SAT publication claims exceed target-free custody",
    )
}

pub(super) fn publish(
    capsule_root: &Path,
    publication: &Path,
    expected: &str,
    out: &Path,
    kind: Kind,
) -> Result<String, String> {
    let root = capsule_root.canonicalize().map_err(|e| e.to_string())?;
    let record = capsule::check_capsule(&root, expected)?;
    kind.check(&record)?;
    unconsumed(&root)?;
    sat_target_build::verify(&root, &record)?;
    if kind == Kind::Scientific {
        sat_target_build::original_preparation(&record.preparation)?;
    }
    fs::create_dir(publication).map_err(|e| e.to_string())?;
    let retained = sidecars(&root)?;
    for name in retained.keys() {
        let target = publication.join(name);
        fs::create_dir_all(target.parent().ok_or("SAT sidecar parent missing")?)
            .map_err(|e| e.to_string())?;
        native::create(&target, &read(&root.join(name), 16 * 1024 * 1024)?)?;
    }
    ordinary_build::archive(
        &root,
        &publication.join("capsule.tar.gz"),
        &full_files(&root, &record)?,
    )?;
    let archive = read(&publication.join("capsule.tar.gz"), 96 * 1024 * 1024)?;
    let header = json!({"schema_version":1,"question":kind.question(),"stage":kind.stage(),
        "registration_sha256":expected,"capsule_root":root,"source_commit":record.source_commit,
        "immutable_file_count":record.immutable_files.as_object().ok_or("SAT inventory missing")?.len(),
        "files":retained,"archive":{"path":"capsule.tar.gz","bytes":archive.len(),"sha256":sha256(&archive)},
        "scientific_worker_calls":0,"execution_admitted":false,"source_bound_execution_admitted":false,
        "fresh_paired_qualification":false,"online_wall_ns":null,"online_speedup":null});
    save(&publication.join("PUBLICATION.json"), &header)?;
    capsule::check_capsule(&root, expected)?;
    unconsumed(&root)?;
    replay(publication, expected, out, kind)
}

pub(super) fn verify(publication: &Path, expected: &str, kind: Kind) -> Result<Value, String> {
    let header = load(&publication.join("PUBLICATION.json"))?;
    claims(&header, kind)?;
    let record = capsule::registration(publication, expected)?;
    kind.check(&record)?;
    unconsumed(publication)?;
    require(
        header["registration_sha256"] == expected
            && header["source_commit"] == record.source_commit
            && header["immutable_file_count"]
                == record
                    .immutable_files
                    .as_object()
                    .ok_or("SAT inventory missing")?
                    .len(),
        "SAT publication descriptor differs",
    )?;
    let original = Path::new(
        header["capsule_root"]
            .as_str()
            .ok_or("SAT original root missing")?,
    );
    require(original.is_absolute(), "SAT original root must be absolute")?;
    let retained = sidecars(publication)?;
    require(
        header["files"] == json!(retained)
            && subset(&record, "build-receipts/")?
                == native::inventory(&publication.join("immutable/build-receipts"))?
            && subset(&record, "assets/controls/")?
                == native::inventory(&publication.join("immutable/assets/controls"))?,
        "SAT publication sidecar or receipt inventory differs",
    )?;
    require(
        canonical_sha(&subset(&record, "source/")?)? == record.source_manifest_sha256
            && record.immutable_files[format!("bin/{}", capsule::WORKER)]["sha256"]
                == record.worker_sha256
            && record.immutable_files["bin/icprog"]["sha256"] == record.auditor_sha256
            && record.immutable_files["assets/bin/prepared-exporter"]["sha256"]
                == record.prepared_exporter_sha256
            && record.immutable_files["assets/bin/prepared-cms"]["sha256"]
                == record.prepared_cms_sha256,
        "SAT archived source or executable pins differ",
    )?;
    sat_target_build::verify_from(publication, &record, original)?;
    require(
        header["archive"]["path"] == "capsule.tar.gz",
        "SAT archive path differs",
    )?;
    let archive = read(&publication.join("capsule.tar.gz"), 96 * 1024 * 1024)?;
    require(
        desc(&archive)
            == json!({"bytes":header["archive"]["bytes"],"sha256":header["archive"]["sha256"]}),
        "SAT published archive bytes differ",
    )?;
    let mut tar = Vec::new();
    MultiGzDecoder::new(archive.as_slice())
        .take(512 * 1024 * 1024 + 1)
        .read_to_end(&mut tar)
        .map_err(|e| e.to_string())?;
    require(
        tar.len() <= 512 * 1024 * 1024,
        "SAT archive exceeds decompression bound",
    )?;
    let members = verify_tar_inventory(&tar, &json!(full_files(publication, &record)?))?;
    Ok(json!({"schema_version":1,"status":kind.status(),
        "registration_sha256":expected,"source_commit":record.source_commit,
        "source_manifest_sha256":record.source_manifest_sha256,
        "prepared_exporter_sha256":record.prepared_exporter_sha256,
        "prepared_cms_sha256":record.prepared_cms_sha256,
        "archive_sha256":sha256(&archive),"archive_regular_files":members,
        "complete_published_capsule_bytes_verified":true,"archived_binaries_executed":0,
        "scientific_worker_calls":0,"execution_admitted":false,
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,
        "online_wall_ns":null,"online_speedup":null}))
}

pub(super) fn replay(
    publication: &Path,
    expected: &str,
    out: &Path,
    kind: Kind,
) -> Result<String, String> {
    let result = verify(publication, expected, kind)?;
    save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

/// Recheck the entire scientific publication before a one-use claim.
/// The archived programs remain data and are never extracted.
pub(super) fn scientific_preflight(
    publication: &Path,
    capsule_root: &Path,
    seal: &str,
) -> Result<Value, String> {
    let publication = publication.canonicalize().map_err(|e| e.to_string())?;
    let result = verify(&publication, seal, Kind::Scientific)?;
    let bytes = read(&publication.join("PUBLICATION.json"), 16 * 1024 * 1024)?;
    let header: Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    require(
        header["capsule_root"] == json!(capsule_root),
        "SAT scientific publication names a different original capsule",
    )?;
    Ok(json!({"schema_version":1,"question":SCIENTIFIC_QUESTION,
        "stage":"scientific-registered-not-dispatched",
        "publication_root":publication,"publication_descriptor_sha256":sha256(&bytes),
        "registration_sha256":seal,"source_manifest_sha256":result["source_manifest_sha256"],
        "archive_sha256":result["archive_sha256"],"scientific_worker_calls":0,
        "execution_admitted":false}))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn custody_rejects_promotion_and_kind_swap() {
        let good = json!({"schema_version":1,"question":BUILD_QUESTION,
            "stage":"validation-only-not-dispatched","scientific_worker_calls":0,
            "execution_admitted":false,"source_bound_execution_admitted":false,
            "fresh_paired_qualification":false,"online_wall_ns":null,"online_speedup":null});
        claims(&good, Kind::Validation).unwrap();
        assert!(claims(&good, Kind::Scientific).is_err());
        for (key, value) in [
            ("scientific_worker_calls", json!(1)),
            ("execution_admitted", json!(true)),
            ("source_bound_execution_admitted", json!(true)),
            ("fresh_paired_qualification", json!(true)),
            ("online_wall_ns", json!(0)),
        ] {
            let mut wrong = good.clone();
            wrong[key] = value;
            assert!(claims(&wrong, Kind::Validation).is_err(), "{key}");
        }
    }
}
