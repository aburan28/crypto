//! Portable, data-only custody of an unconsumed ordinary build or preregistration.
//! Publication and replay never run the solver or any archived executable.
use super::{
    ordinary_control::{self, capsule},
    sat_control::native,
};
use crypto_lib::cryptanalysis::{
    prepared_control_archive::{safe_path, verify_tar_inventory},
    prepared_ordinary::Family,
    prepared_sat_control::{canonical_sha, sha256},
};
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression};
use native::{load, read, require, save};
use serde_json::{json, Value};
use std::{
    fs,
    io::{Read, Write},
    path::Path,
};
const FILES: [&str; 4] = [
    "config.json",
    "host-context.json",
    "registration.json",
    "seal.json",
];
const QUESTION: &str = "ordinary-build-custody-n17-v1";
const SCIENTIFIC_QUESTION: &str = "ordinary-scientific-registration-custody-n17-v1";
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Validation,
    Scientific,
}
impl Kind {
    fn question(self) -> &'static str {
        match self {
            Self::Validation => QUESTION,
            Self::Scientific => SCIENTIFIC_QUESTION,
        }
    }
    fn stage(self) -> &'static str {
        match self {
            Self::Validation => "validation-only-not-dispatched",
            Self::Scientific => "scientific-registered-not-dispatched",
        }
    }
    fn status(self) -> &'static str {
        match self {
            Self::Validation => "PASS_DATA_ONLY_ORDINARY_BUILD_CUSTODY",
            Self::Scientific => "PASS_DATA_ONLY_ORDINARY_SCIENTIFIC_REGISTRATION_CUSTODY",
        }
    }
    fn check(
        self,
        record: &capsule::Registration,
        hash: &str,
        expected: &str,
    ) -> Result<(), String> {
        match self {
            Self::Validation => validation(record, hash, expected),
            Self::Scientific => require(
                !record.validation_only && hash == expected,
                "scientific ordinary registration requires its external executable seal",
            ),
        }
    }
}
fn desc(bytes: &[u8]) -> Value {
    json!({"bytes":bytes.len(),"sha256":sha256(bytes)})
}
fn immutable_subset(record: &capsule::Registration, prefix: &str) -> Result<Value, String> {
    Ok(json!(record
        .immutable_files
        .as_object()
        .ok_or("missing immutable files")?
        .iter()
        .filter_map(|(name, desc)| name
            .strip_prefix(prefix)
            .map(|n| (n.to_string(), desc.clone())))
        .collect::<serde_json::Map<String, Value>>()))
}
fn unconsumed(root: &Path) -> Result<(), String> {
    require(
        !root
            .join("consumed.json")
            .try_exists()
            .map_err(|e| e.to_string())?,
        "ordinary capsule was consumed",
    )
}
fn validation(record: &capsule::Registration, hash: &str, expected: &str) -> Result<(), String> {
    require(
        record.validation_only && hash == expected,
        "build custody requires the external validation-only seal",
    )
}
// Shared strict USTAR writer for ordinary and target build custody. It streams
// only registered regular files and never extracts or executes archived data.
pub(super) fn archive(
    root: &Path,
    out: &Path,
    files: &serde_json::Map<String, Value>,
) -> Result<(), String> {
    let file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(out)
        .map_err(|e| e.to_string())?;
    let mut gz = GzEncoder::new(file, Compression::default());
    let mut size = 1024usize;
    for (index, (name, expected)) in files.iter().enumerate() {
        if index.is_multiple_of(512) {
            eprintln!(
                "native build custody: serializing member {index} of {}",
                files.len()
            );
        }
        require(safe_path(name), "unsafe custody member")?;
        let bytes = read(&root.join(name), 256 * 1024 * 1024)?;
        require(
            desc(&bytes) == *expected,
            "custody member changed during serialization",
        )?;
        let (prefix, leaf) = if name.len() <= 100 {
            ("", name.as_str())
        } else {
            name.match_indices('/')
                .map(|(i, _)| (&name[..i], &name[i + 1..]))
                .find(|(p, n)| p.len() <= 155 && n.len() <= 100)
                .ok_or("custody path exceeds USTAR")?
        };
        let mut h = [0u8; 512];
        h[..leaf.len()].copy_from_slice(leaf.as_bytes());
        h[345..345 + prefix.len()].copy_from_slice(prefix.as_bytes());
        h[100..108].copy_from_slice(b"0000644\0");
        h[108..116].copy_from_slice(b"0000000\0");
        h[116..124].copy_from_slice(b"0000000\0");
        h[124..136].copy_from_slice(format!("{:011o}\0", bytes.len()).as_bytes());
        h[136..148].copy_from_slice(b"00000000000\0");
        h[148..156].fill(b' ');
        h[156] = b'0';
        h[257..265].copy_from_slice(b"ustar\x0000");
        let sum = h.iter().map(|&b| b as usize).sum::<usize>();
        h[148..156].copy_from_slice(format!("{sum:06o}\0 ").as_bytes());
        let padding = (512 - bytes.len() % 512) % 512;
        size = size
            .checked_add(512 + bytes.len() + padding)
            .ok_or("custody size overflow")?;
        require(
            size <= 512 * 1024 * 1024,
            "custody exceeds decompression bound",
        )?;
        gz.write_all(&h)
            .and_then(|_| gz.write_all(&bytes))
            .and_then(|_| gz.write_all(&[0u8; 512][..padding]))
            .map_err(|e| e.to_string())?;
    }
    gz.write_all(&[0u8; 1024]).map_err(|e| e.to_string())?;
    gz.finish()
        .map_err(|e| e.to_string())?
        .sync_all()
        .map_err(|e| e.to_string())
}
fn sidecars(root: &Path) -> Result<serde_json::Map<String, Value>, String> {
    let mut files = serde_json::Map::new();
    for name in FILES {
        files.insert(
            name.into(),
            desc(&read(&root.join(name), 16 * 1024 * 1024)?),
        );
    }
    for (name, desc) in native::inventory(&root.join("immutable/build-receipts"))?
        .as_object()
        .ok_or("missing build receipts")?
    {
        files.insert(format!("immutable/build-receipts/{name}"), desc.clone());
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
    for (name, desc) in record
        .immutable_files
        .as_object()
        .ok_or("missing immutable inventory")?
    {
        files.insert(format!("immutable/{name}"), desc.clone());
    }
    Ok(files)
}
#[cfg(test)]
fn claims(header: &Value) -> Result<(), String> {
    claims_for(header, Kind::Validation)
}
fn claims_for(header: &Value, kind: Kind) -> Result<(), String> {
    require(
        header["schema_version"] == 1
            && header["question"] == kind.question()
            && header["stage"] == kind.stage()
            && header["scientific_worker_calls"] == 0
            && header["execution_admitted"] == false
            && header["full_goal_complete"] == false
            && (if kind == Kind::Scientific {
                header["source_bound_execution_admitted"] == false
            } else {
                header
                    .get("source_bound_execution_admitted")
                    .is_none_or(|v| v == false)
            })
            && header["online_wall_ns"].is_null()
            && header["online_speedup"].is_null(),
        "build custody scope or claims differ",
    )
}
pub fn publish(
    root: &Path,
    publication: &Path,
    expected: &str,
    out: &Path,
) -> Result<String, String> {
    publish_for(root, publication, expected, out, Kind::Validation)
}
pub fn publish_scientific(
    root: &Path,
    publication: &Path,
    expected: &str,
    out: &Path,
) -> Result<String, String> {
    publish_for(root, publication, expected, out, Kind::Scientific)
}
fn publish_for(
    root: &Path,
    publication: &Path,
    expected: &str,
    out: &Path,
    kind: Kind,
) -> Result<String, String> {
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let (record, hash) = capsule::check_capsule(&root)?;
    kind.check(&record, &hash, expected)?;
    unconsumed(&root)?;
    ordinary_control::verify_build_from(&root, &record, &root)?;
    fs::create_dir(publication).map_err(|e| e.to_string())?;
    let retained = sidecars(&root)?;
    for name in retained.keys() {
        let target = publication.join(name);
        fs::create_dir_all(target.parent().ok_or("missing sidecar parent")?)
            .map_err(|e| e.to_string())?;
        native::create(&target, &read(&root.join(name), 16 * 1024 * 1024)?)?;
    }
    eprintln!("ordinary build custody: original source/build checks passed; retaining archive");
    archive(
        &root,
        &publication.join("capsule.tar.gz"),
        &full_files(&root, &record)?,
    )?;
    eprintln!("ordinary build custody: archive stream completed; checking descriptor");
    let bytes = read(&publication.join("capsule.tar.gz"), 96 * 1024 * 1024)?;
    eprintln!("ordinary build custody: bounded archive read completed");
    let header = json!({"schema_version":1,"question":kind.question(),"stage":kind.stage(),
        "registration_sha256":hash,"capsule_root":root,"source_commit":record.source_commit,
        "immutable_file_count":record.immutable_files.as_object().ok_or("missing inventory")?.len(),
        "files":retained,"archive":{"path":"capsule.tar.gz","bytes":bytes.len(),"sha256":sha256(&bytes)},
        "scientific_worker_calls":0,"execution_admitted":false,
        "source_bound_execution_admitted":false,"full_goal_complete":false,
        "online_wall_ns":null,"online_speedup":null});
    save(&publication.join("PUBLICATION.json"), &header)?;
    eprintln!(
        "ordinary build custody: publication descriptor retained; rechecking original capsule"
    );
    require(
        capsule::check_capsule(&root)?.1 == hash,
        "ordinary capsule changed during publication",
    )?;
    unconsumed(&root)?;
    eprintln!(
        "ordinary build custody: original capsule unchanged; verifying archived bytes as data"
    );
    replay_for(publication, expected, out, kind)
}
#[cfg(test)]
fn verify(publication: &Path, expected: &str) -> Result<Value, String> {
    verify_for(publication, expected, Kind::Validation)
}
fn verify_for(publication: &Path, expected: &str, kind: Kind) -> Result<Value, String> {
    let header = load(&publication.join("PUBLICATION.json"))?;
    claims_for(&header, kind)?;
    let (record, hash) = capsule::registration(publication)?;
    kind.check(&record, &hash, expected)?;
    unconsumed(publication)?;
    require(
        header["registration_sha256"] == hash
            && header["source_commit"] == record.source_commit
            && header["immutable_file_count"]
                == record
                    .immutable_files
                    .as_object()
                    .ok_or("missing inventory")?
                    .len(),
        "custody descriptor differs",
    )?;
    let original = Path::new(
        header["capsule_root"]
            .as_str()
            .ok_or("missing original capsule root")?,
    );
    require(
        original.is_absolute(),
        "original capsule path is not absolute",
    )?;
    let retained = sidecars(publication)?;
    require(
        header["files"] == json!(retained),
        "build custody sidecars changed",
    )?;
    require(
        immutable_subset(&record, "build-receipts/")?
            == native::inventory(&publication.join("immutable/build-receipts"))?,
        "build receipt set differs from original sealed inventory",
    )?;
    require(
        canonical_sha(&immutable_subset(&record, "source/")?)? == record.source_manifest_sha256
            && record.immutable_files["bin/prepared_ordinary_worker"]["sha256"]
                == record.worker_sha256
            && record.immutable_files["bin/icprog"]["sha256"] == record.auditor_sha256,
        "archived source/binary pins differ",
    )?;
    ordinary_control::verify_build_from(publication, &record, original)?;
    let cfg = capsule::config(publication)?;
    if cfg.plan.family == Family::Cryptominisat {
        let assets = load(&publication.join("immutable/build-receipts/native-assets.json"))?;
        require(
            record.immutable_files["native-archive/assets.tar.gz"]["sha256"] == native::ARCHIVE_SHA
                && record.immutable_files["assets/bin/cms"]["sha256"] == native::CMS_SHA
                && record.immutable_files["assets/bin/exporter"]["sha256"] == native::EXPORTER_SHA
                && assets["archive_sha256"] == native::ARCHIVE_SHA
                && assets["files"] == immutable_subset(&record, "assets/")?,
            "retained native asset pins/receipt differ",
        )?;
    }
    require(
        header["archive"]["path"] == "capsule.tar.gz",
        "custody archive path differs",
    )?;
    let bytes = read(&publication.join("capsule.tar.gz"), 96 * 1024 * 1024)?;
    require(
        desc(&bytes)
            == json!({"bytes":header["archive"]["bytes"],"sha256":header["archive"]["sha256"]}),
        "custody archive bytes differ",
    )?;
    let mut tar = Vec::new();
    MultiGzDecoder::new(bytes.as_slice())
        .take(512 * 1024 * 1024 + 1)
        .read_to_end(&mut tar)
        .map_err(|e| e.to_string())?;
    require(
        tar.len() <= 512 * 1024 * 1024,
        "custody decompression exceeded bound",
    )?;
    let members = verify_tar_inventory(&tar, &json!(full_files(publication, &record)?))?;
    let mut result = json!({"schema_version":1,"status":kind.status(),
        "registration_sha256":hash,"source_commit":record.source_commit,"source_manifest_sha256":record.source_manifest_sha256,
        "archive_sha256":sha256(&bytes),"archive_regular_files":members,"complete_published_capsule_bytes_verified":true,
        "archived_binaries_executed":0,"scientific_worker_calls":0,"execution_admitted":false,
        "source_bound_execution_admitted":false,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,
        "full_goal_complete":false,"online_wall_ns":null,"online_speedup":null});
    if kind == Kind::Scientific {
        result["scientific_registration_published"] = json!(true);
    }
    Ok(result)
}
pub fn replay(publication: &Path, expected: &str, out: &Path) -> Result<String, String> {
    replay_for(publication, expected, out, Kind::Validation)
}
pub fn replay_scientific(publication: &Path, expected: &str, out: &Path) -> Result<String, String> {
    replay_for(publication, expected, out, Kind::Scientific)
}
/// A complete data-only prerequisite for the original one-use controller.
/// The absolute original capsule path must match the one in the publication.
pub(super) fn scientific_preflight(
    publication: &Path,
    capsule: &Path,
    expected: &str,
) -> Result<Value, String> {
    let publication = publication.canonicalize().map_err(|e| e.to_string())?;
    let receipt = verify_for(&publication, expected, Kind::Scientific)?;
    let descriptor = read(&publication.join("PUBLICATION.json"), 16 * 1024 * 1024)?;
    let header: Value = serde_json::from_slice(&descriptor).map_err(|e| e.to_string())?;
    require(
        header["capsule_root"] == json!(capsule),
        "published ordinary registration names a different original capsule",
    )?;
    Ok(json!({"schema_version":1,"question":SCIENTIFIC_QUESTION,
        "stage":"scientific-registered-not-dispatched",
        "publication_root":publication,"publication_descriptor_sha256":sha256(&descriptor),
        "registration_sha256":expected,"source_manifest_sha256":receipt["source_manifest_sha256"],
        "archive_sha256":receipt["archive_sha256"],"scientific_worker_calls":0,
        "execution_admitted":false}))
}
fn replay_for(
    publication: &Path,
    expected: &str,
    out: &Path,
    kind: Kind,
) -> Result<String, String> {
    let receipt = verify_for(publication, expected, kind)?;
    save(out, &receipt)?;
    serde_json::to_string_pretty(&receipt).map_err(|e| e.to_string())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn archive_checks_actual_bytes_and_round_trips_long_registered_names() {
        let root =
            std::env::temp_dir().join(format!("ordinary-build-archive-{}", std::process::id()));
        fs::create_dir(&root).unwrap();
        let name = format!("{}/source", "long".repeat(30));
        fs::create_dir(root.join("long".repeat(30))).unwrap();
        fs::write(root.join(&name), b"control").unwrap();
        let files = [(name.clone(), desc(b"control"))].into_iter().collect();
        archive(&root, &root.join("good.gz"), &files).unwrap();
        let mut tar = Vec::new();
        MultiGzDecoder::new(read(&root.join("good.gz"), 1024 * 1024).unwrap().as_slice())
            .read_to_end(&mut tar)
            .unwrap();
        assert_eq!(verify_tar_inventory(&tar, &json!(files)).unwrap(), 1);
        fs::write(root.join(name), b"changed").unwrap();
        assert!(archive(&root, &root.join("bad.gz"), &files).is_err());
        fs::remove_dir_all(root).unwrap();
    }
    #[test]
    fn custody_cannot_claim_execution_or_change_its_question() {
        let header = json!({"schema_version":1,"question":QUESTION,"stage":"validation-only-not-dispatched",
            "scientific_worker_calls":0,"execution_admitted":false,"full_goal_complete":false,"online_wall_ns":null,"online_speedup":null});
        claims(&header).unwrap();
        let mut wrong = header.clone();
        wrong["execution_admitted"] = json!(true);
        assert!(claims(&wrong).is_err());
        let mut wrong = header.clone();
        wrong["question"] = json!("scientific-qualification");
        assert!(claims(&wrong).is_err());
        let mut wrong = header;
        wrong["online_wall_ns"] = json!(0);
        assert!(claims(&wrong).is_err());
    }
    #[test]
    fn scientific_registration_is_distinct_from_validation_and_execution() {
        let publication = Path::new(
            "research/ic_candidate_tournament_20260915/goal_20260924/native-ordinary-controller-v1/build-validation-v2",
        );
        let (record, seal) = capsule::registration(publication).unwrap();
        assert!(Kind::Validation.check(&record, &seal, &seal).is_ok());
        assert!(Kind::Scientific.check(&record, &seal, &seal).is_err());
        let mut scientific: capsule::Registration = serde_json::from_value(json!(record)).unwrap();
        scientific.validation_only = false;
        assert!(Kind::Scientific.check(&scientific, &seal, &seal).is_ok());
        assert!(Kind::Scientific
            .check(&scientific, &seal, &"0".repeat(64))
            .is_err());
        let changed = std::env::temp_dir().join(format!(
            "ordinary-scientific-seal-control-{}",
            std::process::id()
        ));
        fs::create_dir(&changed).unwrap();
        for name in FILES {
            fs::copy(publication.join(name), changed.join(name)).unwrap();
        }
        let mut changed_registration = load(&changed.join("registration.json")).unwrap();
        changed_registration["validation_only"] = json!(false);
        fs::write(
            changed.join("registration.json"),
            serde_json::to_vec(&changed_registration).unwrap(),
        )
        .unwrap();
        assert!(capsule::registration(&changed).is_err());
        fs::remove_dir_all(changed).unwrap();
        let good = json!({"schema_version":1,"question":SCIENTIFIC_QUESTION,
            "stage":"scientific-registered-not-dispatched","scientific_worker_calls":0,
            "execution_admitted":false,"source_bound_execution_admitted":false,
            "full_goal_complete":false,"online_wall_ns":null,"online_speedup":null});
        claims_for(&good, Kind::Scientific).unwrap();
        assert!(claims_for(&good, Kind::Validation).is_err());
        for (key, value) in [
            ("scientific_worker_calls", json!(1)),
            ("execution_admitted", json!(true)),
            ("source_bound_execution_admitted", json!(true)),
            ("online_wall_ns", json!(0)),
        ] {
            let mut changed = good.clone();
            changed[key] = value;
            assert!(claims_for(&changed, Kind::Scientific).is_err(), "{key}");
        }
    }
    fn retained_control(label: &str) -> std::path::PathBuf {
        let original=Path::new("research/ic_candidate_tournament_20260915/goal_20260924/native-ordinary-controller-v1/build-validation-v2");
        let root =
            std::env::temp_dir().join(format!("ordinary-build-{label}-{}", std::process::id()));
        fs::create_dir(&root).unwrap();
        for name in sidecars(original)
            .unwrap()
            .keys()
            .chain(std::iter::once(&"PUBLICATION.json".to_string()))
        {
            let to = root.join(name);
            fs::create_dir_all(to.parent().unwrap()).unwrap();
            native::create(&to, &read(&original.join(name), 16 * 1024 * 1024).unwrap()).unwrap();
        }
        root
    }
    const SEAL: &str = "02e12b456947c8034f096b47c16dd7b73eb6c15cedb6c0b79164f887ec1406f3";
    #[test]
    fn original_sealed_receipt_set_cannot_be_weakened_by_refreshing_outer_manifest() {
        let root = retained_control("missing");
        fs::remove_file(root.join("immutable/build-receipts/worker-identity.receipt.json"))
            .unwrap();
        let mut header = load(&root.join("PUBLICATION.json")).unwrap();
        header["files"] = json!(sidecars(&root).unwrap());
        fs::write(
            root.join("PUBLICATION.json"),
            serde_json::to_vec(&header).unwrap(),
        )
        .unwrap();
        assert_eq!(
            verify(&root, SEAL).unwrap_err(),
            "build receipt set differs from original sealed inventory"
        );
        fs::remove_dir_all(root).unwrap();
        let root = retained_control("changed");
        let path = root.join("immutable/build-receipts/build.receipt.json");
        let mut receipt = load(&path).unwrap();
        receipt["argv"] = json!(["build", "--release", "--features", "other"]);
        fs::write(path, serde_json::to_vec(&receipt).unwrap()).unwrap();
        let mut header = load(&root.join("PUBLICATION.json")).unwrap();
        header["files"] = json!(sidecars(&root).unwrap());
        fs::write(
            root.join("PUBLICATION.json"),
            serde_json::to_vec(&header).unwrap(),
        )
        .unwrap();
        assert_eq!(
            verify(&root, SEAL).unwrap_err(),
            "build receipt set differs from original sealed inventory"
        );
        fs::remove_dir_all(root).unwrap();
    }
    #[test]
    fn retained_publication_requires_external_seal_and_rejects_runtime_claims() {
        let root = retained_control("claims");
        assert_eq!(
            verify(&root, &"0".repeat(64)).unwrap_err(),
            "build custody requires the external validation-only seal"
        );
        let mut header = load(&root.join("PUBLICATION.json")).unwrap();
        header["source_bound_execution_admitted"] = json!(true);
        fs::write(
            root.join("PUBLICATION.json"),
            serde_json::to_vec(&header).unwrap(),
        )
        .unwrap();
        assert_eq!(
            verify(&root, SEAL).unwrap_err(),
            "build custody scope or claims differ"
        );
        fs::remove_dir_all(root).unwrap();
    }
}
