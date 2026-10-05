//! Complete offline source/build custody for the new n17 target scope.
use super::{
    ordinary_control as files,
    sat_control::{native, native_step},
    target_control::{self, capsule},
};
use crypto_lib::cryptanalysis::prepared_sat_control::{canonical_sha, sha256};
use native::{load, read, require, save};
use serde_json::{json, Value};
use std::{fs, path::Path, process::Command};

pub(super) struct FreezeRequest<'a> {
    pub root: &'a Path,
    pub out: &'a Path,
    pub config: &'a Path,
    pub preparation: &'a Path,
    pub cargo: &'a Path,
    pub rustc: &'a Path,
    pub host: &'a Path,
    pub validation_only: bool,
}
pub(super) fn freeze(request: FreezeRequest<'_>) -> Result<String, String> {
    let FreezeRequest {
        root,
        out,
        config,
        preparation,
        cargo,
        rustc,
        host,
        validation_only,
    } = request;
    native::enforce_hardware()?;
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let cargo = cargo.canonicalize().map_err(|e| e.to_string())?;
    let rustc = rustc.canonicalize().map_err(|e| e.to_string())?;
    let cfg: capsule::Config =
        serde_json::from_slice(&read(config, 65536)?).map_err(|e| e.to_string())?;
    cfg.validate()?;
    let prep: capsule::PreparationBinding =
        serde_json::from_slice(&read(preparation, 65536)?).map_err(|e| e.to_string())?;
    for hash in [
        &prep.registration_sha256,
        &prep.auditor_sha256,
        &prep.producer_sha256,
        &prep.mathematics_sha256,
    ] {
        require(
            target_control::journal::digest(hash),
            "invalid original preparation pin",
        )?;
    }
    require(
        prep.capsule.is_absolute() && prep.execution.is_absolute(),
        "original preparation paths must be absolute",
    )?;
    if !validation_only {
        target_control::original_preparation(&prep)?;
    }
    let host = load(host)?;
    require(
        host["os"] == "macos"
            && host["architecture"] == "aarch64"
            && host["physical_hardware"] == true
            && host["calibrated_performance_environment"] == false,
        "target host context cannot claim calibrated performance",
    )?;
    let head = Command::new("git")
        .args(["-c", "gc.auto=0", "rev-parse", "HEAD"])
        .current_dir(&root)
        .output()
        .map_err(|e| e.to_string())?;
    require(head.status.success(), "target source lacks committed HEAD")?;
    let commit = String::from_utf8(head.stdout)
        .map_err(|e| e.to_string())?
        .trim()
        .to_string();
    let dirty = Command::new("git")
        .args([
            "-c",
            "gc.auto=0",
            "status",
            "--porcelain",
            "--untracked-files=all",
            "--",
            "src",
            "examples",
            "benches",
            "tests",
            "docs",
            "Cargo.toml",
            "build.rs",
        ])
        .current_dir(&root)
        .output()
        .map_err(|e| e.to_string())?;
    require(
        dirty.status.success() && dirty.stdout.is_empty(),
        "commit all target source and embedded data before freeze",
    )?;
    fs::create_dir(out).map_err(|e| e.to_string())?;
    let out = out.canonicalize().map_err(|e| e.to_string())?;
    let immutable = out.join("immutable");
    fs::create_dir(&immutable).map_err(|e| e.to_string())?;
    let source = immutable.join("source");
    fs::create_dir(&source).map_err(|e| e.to_string())?;
    for name in ["src", "examples", "benches", "tests", "docs"] {
        native::copy_tree(&root.join(name), &source.join(name))?;
    }
    for name in ["Cargo.toml", "Cargo.lock"] {
        files::copy_file(&root, &source, name)?;
    }
    if root
        .join("build.rs")
        .try_exists()
        .map_err(|e| e.to_string())?
    {
        files::copy_file(&root, &source, "build.rs")?;
    }
    let receipts = immutable.join("build-receipts");
    fs::create_dir(&receipts).map_err(|e| e.to_string())?;
    let path = std::env::var("PATH").map_err(|e| e.to_string())?;
    native_step(
        &cargo,
        &[
            "vendor".into(),
            "--locked".into(),
            "--offline".into(),
            source.join("vendor").to_string_lossy().into_owned(),
        ],
        &root,
        &receipts.join("vendor.log"),
        &[("PATH".into(), path.clone()), ("LC_ALL".into(), "C".into())],
    )?;
    fs::create_dir(source.join(".cargo")).map_err(|e| e.to_string())?;
    native::create(&source.join(".cargo/config.toml"),b"[source.crates-io]\nreplace-with = \"vendored-sources\"\n[source.vendored-sources]\ndirectory = \"vendor\"\n")?;
    files::reject_ancestor_config(&source)?;
    let snapshot = native::inventory(&source)?;
    let source_sha = canonical_sha(&snapshot)?;
    let cargo_sha = sha256(&read(&cargo, 128 * 1024 * 1024)?);
    let rustc_sha = sha256(&read(&rustc, 128 * 1024 * 1024)?);
    let build = json!({"schema_version":1,"scope":capsule::SCOPE,"source_commit":commit,"checkout_root":root,
        "source_manifest_sha256":source_sha,"cargo_sha256":cargo_sha,"rustc_sha256":rustc_sha,"profile":"release","features":"default",
        "hardware":{"os":"macos","architecture":"aarch64"},"algorithm_flags":"none","native_archive_sha256":null});
    let build_sha = canonical_sha(&build)?;
    save(&receipts.join("build-identity.json"), &build)?;
    let home = out.join("build-home");
    fs::create_dir(&home).map_err(|e| e.to_string())?;
    let target = out.join("build-target");
    let env = vec![
        ("PATH".into(), path),
        ("LC_ALL".into(), "C".into()),
        ("RUSTC".into(), rustc.to_string_lossy().into_owned()),
        ("CARGO_HOME".into(), home.to_string_lossy().into_owned()),
        (
            "CARGO_TARGET_DIR".into(),
            target.to_string_lossy().into_owned(),
        ),
        ("CARGO_INCREMENTAL".into(), "0".into()),
        (
            "IC_TARGET_SOURCE_MANIFEST_SHA256".into(),
            source_sha.clone(),
        ),
        ("IC_TARGET_BUILD_SHA256".into(), build_sha.clone()),
    ];
    for (tool, name) in [(&rustc, "rustc"), (&cargo, "cargo")] {
        native_step(
            tool,
            &["--version".into(), "--verbose".into()],
            &source,
            &receipts.join(format!("{name}.log")),
            &env,
        )?;
    }
    native_step(
        &cargo,
        &[
            "build".into(),
            "--locked".into(),
            "--offline".into(),
            "--release".into(),
            "--bin".into(),
            capsule::WORKER.into(),
            "--bin".into(),
            "icprog".into(),
        ],
        &source,
        &receipts.join("build.log"),
        &env,
    )?;
    native::check_tree(&source, &snapshot)?;
    require(
        sha256(&read(&cargo, 128 * 1024 * 1024)?) == cargo_sha
            && sha256(&read(&rustc, 128 * 1024 * 1024)?) == rustc_sha,
        "target build tools changed during freeze",
    )?;
    let bin = immutable.join("bin");
    fs::create_dir(&bin).map_err(|e| e.to_string())?;
    for name in [capsule::WORKER, "icprog"] {
        files::executable_copy(&target.join("release").join(name), &bin.join(name))?;
    }
    let identity = json!({"schema_version":1,"source_manifest_sha256":source_sha,"build_sha256":build_sha,"target_os":"macos","target_arch":"aarch64"});
    native_step(
        &bin.join(capsule::WORKER),
        &["build-identity".into()],
        &out,
        &receipts.join("worker-identity.log"),
        &capsule::environment().into_iter().collect::<Vec<_>>(),
    )?;
    capsule::require_identity(&load_log(&receipts.join("worker-identity.log"))?, &identity)?;
    save(&out.join("config.json"), &json!(cfg))?;
    save(&out.join("host-context.json"), &host)?;
    let record = capsule::Registration {
        schema_version: 1,
        scope: capsule::SCOPE.into(),
        source_commit: commit,
        hardware: json!({"os":"macos","architecture":"aarch64"}),
        immutable_files: native::inventory(&immutable)?,
        source_manifest_sha256: source_sha,
        config_sha256: sha256(&read(&out.join("config.json"), 65536)?),
        host_context_sha256: sha256(&read(&out.join("host-context.json"), 65536)?),
        worker_sha256: sha256(&read(&bin.join(capsule::WORKER), 128 * 1024 * 1024)?),
        auditor_sha256: sha256(&read(&bin.join("icprog"), 128 * 1024 * 1024)?),
        worker_build_identity: identity,
        runtime_environment: capsule::environment(),
        preparation: prep,
        validation_only,
    };
    save(&out.join("registration.json"), &json!(record))?;
    let seal = canonical_sha(&json!(record))?;
    save(&out.join("seal.json"), &json!({"registration_sha256":seal}))?;
    let bound = capsule::check_capsule(&out, &seal)?;
    verify(&out, &bound)?;
    serde_json::to_string_pretty(&json!({"registration_sha256":seal,"validation_only":validation_only,"scientific_worker_calls":0,"source_bound_execution_admitted":false,"full_goal_complete":false,"online_speedup":null})).map_err(|e|e.to_string())
}
fn load_log(path: &Path) -> Result<Value, String> {
    super::target_math::parse(&read(path, 65536)?)
}
pub(super) fn verify(root: &Path, record: &capsule::Registration) -> Result<(), String> {
    files::reject_ancestor_config(&root.join("immutable/source"))?;
    verify_from(root, record, root)
}
/// Check retained original invocation receipts against archived sidecars.
/// `root` holds the bytes being replayed; `original` names the build's actual
/// capsule location recorded at freeze time. No archived program is run.
pub(super) fn verify_from(
    root: &Path,
    record: &capsule::Registration,
    original: &Path,
) -> Result<(), String> {
    let dir = root.join("immutable/build-receipts");
    let build = load(&dir.join("build-identity.json"))?;
    require(
        build["schema_version"] == 1
            && build["scope"] == capsule::SCOPE
            && build["source_commit"] == record.source_commit
            && canonical_sha(&build)? == record.worker_build_identity["build_sha256"]
            && build["source_manifest_sha256"] == record.source_manifest_sha256
            && build["profile"] == "release"
            && build["features"] == "default"
            && build["algorithm_flags"] == "none"
            && build["hardware"] == record.hardware
            && build["native_archive_sha256"].is_null(),
        "target compiled source/build flags differ",
    )?;
    let source = original.join("immutable/source");
    let build_env = load(&dir.join("build.receipt.json"))?["environment"].clone();
    for name in ["vendor", "rustc", "cargo", "build", "worker-identity"] {
        let receipt = load(&dir.join(format!("{name}.receipt.json")))?;
        let bytes = read(&dir.join(format!("{name}.log")), 16 * 1024 * 1024)?;
        require(
            receipt["wall_ns"]
                .as_str()
                .and_then(|s| s.parse::<u128>().ok())
                .is_some_and(|n| n > 0),
            "target tool receipt lacks measured interval",
        )?;
        let pin = if name == "worker-identity" {
            json!(record.worker_sha256)
        } else if name == "rustc" {
            build["rustc_sha256"].clone()
        } else {
            build["cargo_sha256"].clone()
        };
        require(
            receipt["exit_code"] == 0
                && receipt["program_sha256_before"] == pin
                && receipt["program_sha256_after"] == pin
                && receipt["log_sha256"] == sha256(&bytes),
            "target tool/output/build receipt differs",
        )?;
        if name == "worker-identity" {
            require(
                load_log(&dir.join("worker-identity.log"))? == record.worker_build_identity
                    && receipt["argv"] == json!(["build-identity"])
                    && receipt["cwd"] == json!(original)
                    && receipt["program"]
                        == json!(original.join("immutable/bin").join(capsule::WORKER))
                    && receipt["environment"]
                        == json!(capsule::environment().into_iter().collect::<Vec<_>>()),
                "target observed worker identity differs",
            )?;
        } else if name == "vendor" {
            let env: Vec<(String, String)> = serde_json::from_value(receipt["environment"].clone())
                .map_err(|e| e.to_string())?;
            require(
                receipt["argv"]
                    == json!(["vendor", "--locked", "--offline", source.join("vendor")])
                    && receipt["cwd"] == build["checkout_root"]
                    && env.len() == 2
                    && env[0].0 == "PATH"
                    && !env[0].1.is_empty()
                    && env[1] == ("LC_ALL".into(), "C".into())
                    && receipt["program"] == load(&dir.join("cargo.receipt.json"))?["program"],
                "target vendor invocation differs",
            )?;
        } else {
            require(
                receipt["cwd"] == json!(source) && receipt["environment"] == build_env,
                "target build cwd/environment differs",
            )?;
            if name == "build" {
                require(
                    receipt["argv"]
                        == json!([
                            "build",
                            "--locked",
                            "--offline",
                            "--release",
                            "--bin",
                            capsule::WORKER,
                            "--bin",
                            "icprog"
                        ]),
                    "target build argv differs",
                )?;
                let pairs: Vec<(String, String)> =
                    serde_json::from_value(build_env.clone()).map_err(|e| e.to_string())?;
                let env: std::collections::BTreeMap<String, String> =
                    pairs.iter().cloned().collect();
                require(
                    pairs.len() == 8
                        && env.len() == 8
                        && env.get("IC_TARGET_SOURCE_MANIFEST_SHA256")
                            == Some(&record.source_manifest_sha256)
                        && env.get("IC_TARGET_BUILD_SHA256").map(String::as_str)
                            == record.worker_build_identity["build_sha256"].as_str()
                        && env.get("CARGO_INCREMENTAL").map(String::as_str) == Some("0")
                        && env.get("LC_ALL").map(String::as_str) == Some("C")
                        && env.get("CARGO_HOME")
                            == Some(&original.join("build-home").to_string_lossy().into_owned())
                        && env.get("CARGO_TARGET_DIR")
                            == Some(&original.join("build-target").to_string_lossy().into_owned())
                        && env.get("RUSTC").map(String::as_str)
                            == load(&dir.join("rustc.receipt.json"))?["program"].as_str()
                        && env.contains_key("PATH"),
                    "target frozen build environment differs",
                )?;
                require(
                    receipt["program"] == load(&dir.join("cargo.receipt.json"))?["program"],
                    "target build tool path differs",
                )?;
            } else {
                require(
                    receipt["argv"] == json!(["--version", "--verbose"]),
                    "target tool version argv differs",
                )?;
            }
        }
    }
    Ok(())
}
