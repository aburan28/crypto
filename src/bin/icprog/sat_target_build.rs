//! Target-free source/build freeze for the external-CMS n17 target worker.
//! This creates a sealed build, never a target card or permission to dispatch.
use super::{
    ordinary_control::{self as files, capsule as preparation},
    ordinary_preparation,
    sat_control::{native, native_step},
    target_math,
};
#[path = "../prepared_sat_target_worker/contract.rs"]
#[allow(dead_code)]
pub(super) mod capsule;
#[path = "../prepared_target_worker/journal.rs"]
#[allow(dead_code)]
mod journal;

use crypto_lib::cryptanalysis::prepared_sat_control::{canonical_sha, sha256};
use native::{load, read, require, save};
use serde_json::{json, Value};
use std::{fs, path::Path, process::Command};

const MARKED_CMS_SHA: &str = "128f074850bbd9d973bd50be26c6ffc161e97c7a365b9234f8af171545fe7263";
const MARKED_AUDIT_SHA: &str = "1cecc6161179445317644a92a4ecc9e13fa714f3d2dc015ea5862818f79e3da5";
const ORIGINAL_PREPARATION_AUDIT_SHA: &str =
    "e2e9caa59ffa61713e4ffd47d98a05850fe4e22bef9594d9431c2c7b24d283be";
const MARKED_AUDIT: &str = "research/ic_candidate_tournament_20260915/goal_20260924/native-sat-target-transport-v1/result-v1/original-audit.json";
const PREPARATION_AUDIT: &str = "research/ic_candidate_tournament_20260915/goal_20260924/native-sat-million-registration-v1/result-v1/original-audit.json";

pub(super) struct FreezeRequest<'a> {
    pub root: &'a Path,
    pub out: &'a Path,
    pub config: &'a Path,
    pub preparation_binding: &'a Path,
    pub marked_cms: &'a Path,
    pub cargo: &'a Path,
    pub rustc: &'a Path,
    pub host_context: &'a Path,
    pub validation_only: bool,
}

fn source_commit(root: &Path) -> Result<String, String> {
    let head = Command::new("git")
        .args(["-c", "gc.auto=0", "rev-parse", "HEAD"])
        .current_dir(root)
        .output()
        .map_err(|e| e.to_string())?;
    require(
        head.status.success(),
        "SAT target source lacks committed HEAD",
    )?;
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
            "Cargo.lock",
            "build.rs",
            MARKED_AUDIT,
            PREPARATION_AUDIT,
        ])
        .current_dir(root)
        .output()
        .map_err(|e| e.to_string())?;
    require(
        dirty.status.success() && dirty.stdout.is_empty(),
        "commit the SAT target source and audited controls before freeze",
    )?;
    Ok(commit)
}

/// Reconstruct the original SAT table without executing a solver. The worker
/// reruns the frozen preparation auditor before it opens its online interval.
pub(super) fn original_preparation(binding: &capsule::PreparationBinding) -> Result<Value, String> {
    let (root, execution, _) = capsule::preparation_paths(binding)?;
    let (record, seal) = preparation::check_capsule(&root)?;
    require(
        seal == binding.registration_sha256,
        "SAT preparation seal differs",
    )?;
    files::verify_build_from(&root, &record, &root)?;
    let producer_bytes = read(&execution.join("producer.json"), 16 * 1024 * 1024)?;
    let producer = target_math::parse_ordinary_producer(&producer_bytes)?;
    let math = producer["mathematical_input"].clone();
    let checked = ordinary_preparation::audit(&math)?;
    require(
        sha256(&producer_bytes) == binding.producer_sha256
            && canonical_sha(&math)? == binding.mathematics_sha256
            && checked["declared_family"] == "cryptominisat"
            && checked["panel_complete"] == true
            && checked["mathematical_preparation_complete"] == true
            && checked["rank"] == 29
            && checked["folded_columns"] == 29
            && checked["usable_points"] == 62
            && checked["independently_recovered_column_logs"]
                .as_array()
                .is_some_and(|logs| logs.len() == 29),
        "SAT original preparation lacks source-bound complete logs",
    )?;
    require(
        read(&execution.join("producer.json"), 16 * 1024 * 1024)? == producer_bytes,
        "SAT preparation producer changed during source freeze",
    )?;
    Ok(math)
}

fn bound_control(path: &Path, expected_sha: &str, status: &str) -> Result<Vec<u8>, String> {
    let bytes = read(path, 16 * 1024 * 1024)?;
    let value: Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    require(
        sha256(&bytes) == expected_sha
            && value["status"] == status
            && value["online_speedup"].is_null(),
        "SAT prerequisite audit differs from the accepted original bytes",
    )?;
    Ok(bytes)
}

fn check_prerequisites(
    marked_bytes: &[u8],
    preparation_bytes: &[u8],
    binding: &capsule::PreparationBinding,
) -> Result<(), String> {
    let marked: Value = serde_json::from_slice(marked_bytes).map_err(|e| e.to_string())?;
    let prep: Value = serde_json::from_slice(preparation_bytes).map_err(|e| e.to_string())?;
    require(
        marked["audited_roles"] == 9
            && marked["disclosed_transport_parity_admitted"] == true
            && marked["source_binary_pins_checked"] == true
            && marked["source_bound_target_admitted"] == false
            && marked["fresh_paired_qualification"] == false
            && marked["native_solver_calls_by_auditor"] == 0
            && prep["registration_sha256"] == binding.registration_sha256
            && prep["checker_sha256"] == binding.auditor_sha256
            && prep["producer_sha256"] == binding.producer_sha256
            && prep["mathematics"]["input_sha256"] == binding.mathematics_sha256
            && prep["mathematics"]["rank"] == 29
            && prep["mathematics"]["folded_columns"] == 29
            && prep["source_bound_execution_admitted"] == true
            && prep["mathematical_preparation_complete"] == true,
        "SAT original preparation or marked transport receipt is not bound",
    )
}

pub(super) fn freeze(request: FreezeRequest<'_>) -> Result<String, String> {
    require(
        request.validation_only,
        "SAT scientific freeze awaits exact-built prepared-exporter parity and one-use controller",
    )?;
    native::enforce_hardware()?;
    let root = request.root.canonicalize().map_err(|e| e.to_string())?;
    let cargo = request.cargo.canonicalize().map_err(|e| e.to_string())?;
    let rustc = request.rustc.canonicalize().map_err(|e| e.to_string())?;
    let marked = request
        .marked_cms
        .canonicalize()
        .map_err(|e| e.to_string())?;
    let cfg: capsule::Config =
        serde_json::from_slice(&read(request.config, 65536)?).map_err(|e| e.to_string())?;
    cfg.validate()?;
    require(
        cfg.plan.max_queries == 8
            && cfg.plan.conflict_budget == 1_000_000
            && cfg.plan.solver_timeout_ms == 120_000
            && cfg.plan.controller_timeout_ms == 900_000
            && cfg.worker_timeout_ms == 900_000,
        "SAT scientific target budget differs from frozen four-arm protocol",
    )?;
    let prep: capsule::PreparationBinding =
        serde_json::from_slice(&read(request.preparation_binding, 65536)?)
            .map_err(|e| e.to_string())?;
    original_preparation(&prep)?;
    let host = load(request.host_context)?;
    require(
        host["os"] == "macos"
            && host["architecture"] == "aarch64"
            && host["physical_hardware"] == true
            && host["calibrated_performance_environment"] == false,
        "SAT target source freeze cannot claim calibrated performance",
    )?;
    let marked_bytes = read(&marked, 16 * 1024 * 1024)?;
    require(
        sha256(&marked_bytes) == MARKED_CMS_SHA,
        "marked CMS differs from audited disclosed control",
    )?;
    let marked_audit = bound_control(
        &root.join(MARKED_AUDIT),
        MARKED_AUDIT_SHA,
        "PASS_DISCLOSED_CMS_TRANSPORT_PARITY",
    )?;
    let preparation_audit = bound_control(
        &root.join(PREPARATION_AUDIT),
        ORIGINAL_PREPARATION_AUDIT_SHA,
        "PASS_NATIVE_SOURCE_BOUND_ORDINARY_PREPARATION_AUDIT",
    )?;
    check_prerequisites(&marked_audit, &preparation_audit, &prep)?;
    let commit = source_commit(&root)?;
    fs::create_dir(request.out).map_err(|e| e.to_string())?;
    let out = request.out.canonicalize().map_err(|e| e.to_string())?;
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
    let assets = immutable.join("assets/bin");
    fs::create_dir_all(&assets).map_err(|e| e.to_string())?;
    files::executable_copy(&marked, &assets.join("prepared-cms"))?;
    let control_dir = immutable.join("assets/controls");
    fs::create_dir_all(&control_dir).map_err(|e| e.to_string())?;
    native::create(&control_dir.join("marked-audit.json"), &marked_audit)?;
    native::create(
        &control_dir.join("preparation-audit.json"),
        &preparation_audit,
    )?;
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
    native::create(
        &source.join(".cargo/config.toml"),
        b"[source.crates-io]\nreplace-with = \"vendored-sources\"\n[source.vendored-sources]\ndirectory = \"vendor\"\n",
    )?;
    files::reject_ancestor_config(&source)?;
    let source_inventory = native::inventory(&source)?;
    let source_sha = canonical_sha(&source_inventory)?;
    let cargo_sha = sha256(&read(&cargo, 128 * 1024 * 1024)?);
    let rustc_sha = sha256(&read(&rustc, 128 * 1024 * 1024)?);
    let build = json!({"schema_version":1,"scope":capsule::SCOPE,
        "source_commit":commit,"checkout_root":root,
        "source_manifest_sha256":source_sha,"cargo_sha256":cargo_sha,
        "rustc_sha256":rustc_sha,"profile":"release","features":"default",
        "algorithm_flags":"none","hardware":{"os":"macos","architecture":"aarch64"},
        "marked_cms_sha256":MARKED_CMS_SHA,"marked_control_sha256":MARKED_AUDIT_SHA,
        "preparation_audit_sha256":ORIGINAL_PREPARATION_AUDIT_SHA});
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
            "IC_SAT_TARGET_SOURCE_MANIFEST_SHA256".into(),
            source_sha.clone(),
        ),
        ("IC_SAT_TARGET_BUILD_SHA256".into(), build_sha.clone()),
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
            "--example".into(),
            "koblitz_pdp_export".into(),
        ],
        &source,
        &receipts.join("build.log"),
        &env,
    )?;
    native::check_tree(&source, &source_inventory)?;
    require(
        sha256(&read(&cargo, 128 * 1024 * 1024)?) == cargo_sha
            && sha256(&read(&rustc, 128 * 1024 * 1024)?) == rustc_sha,
        "SAT frozen build tools changed",
    )?;
    let bin = immutable.join("bin");
    fs::create_dir(&bin).map_err(|e| e.to_string())?;
    for name in [capsule::WORKER, "icprog"] {
        files::executable_copy(&target.join("release").join(name), &bin.join(name))?;
    }
    files::executable_copy(
        &target.join("release/examples/koblitz_pdp_export"),
        &assets.join("prepared-exporter"),
    )?;
    save(
        &receipts.join("build-artifacts.json"),
        &json!({"schema_version":1,"scope":capsule::SCOPE,
            "worker_sha256":sha256(&read(&bin.join(capsule::WORKER),128*1024*1024)?),
            "icprog_sha256":sha256(&read(&bin.join("icprog"),128*1024*1024)?),
            "prepared_exporter_sha256":sha256(&read(&assets.join("prepared-exporter"),8*1024*1024)?),
            "prepared_cms_sha256":MARKED_CMS_SHA,
            "build_command_receipt_sha256":sha256(&read(&receipts.join("build.receipt.json"),65536)?)}),
    )?;
    let identity = json!({"schema_version":1,"source_manifest_sha256":source_sha,
        "build_sha256":build_sha,"target_os":"macos","target_arch":"aarch64"});
    native_step(
        &bin.join(capsule::WORKER),
        &["build-identity".into()],
        &out,
        &receipts.join("worker-identity.log"),
        &capsule::environment().into_iter().collect::<Vec<_>>(),
    )?;
    capsule::require_identity(&load(&receipts.join("worker-identity.log"))?, &identity)?;
    save(&out.join("config.json"), &json!(cfg))?;
    save(&out.join("host-context.json"), &host)?;
    let record = capsule::Registration {
        schema_version: 1,
        scope: capsule::SCOPE.into(),
        source_commit: commit.clone(),
        hardware: json!({"os":"macos","architecture":"aarch64"}),
        immutable_files: native::inventory(&immutable)?,
        source_manifest_sha256: source_sha.clone(),
        config_sha256: sha256(&read(&out.join("config.json"), 65536)?),
        host_context_sha256: sha256(&read(&out.join("host-context.json"), 65536)?),
        worker_sha256: sha256(&read(&bin.join(capsule::WORKER), 128 * 1024 * 1024)?),
        auditor_sha256: sha256(&read(&bin.join("icprog"), 128 * 1024 * 1024)?),
        prepared_exporter_sha256: sha256(&read(
            &assets.join("prepared-exporter"),
            8 * 1024 * 1024,
        )?),
        prepared_cms_sha256: MARKED_CMS_SHA.into(),
        worker_build_identity: identity,
        runtime_environment: capsule::environment(),
        preparation: prep,
        validation_only: request.validation_only,
    };
    save(&out.join("registration.json"), &json!(record))?;
    let seal = canonical_sha(&json!(record))?;
    save(&out.join("seal.json"), &json!({"registration_sha256":seal}))?;
    let bound = capsule::check_capsule(&out, &seal)?;
    verify(&out, &bound)?;
    serde_json::to_string_pretty(&json!({"registration_sha256":seal,
        "validation_only":request.validation_only,"source_commit":commit,
        "source_manifest_sha256":source_sha,"prepared_exporter_sha256":bound.prepared_exporter_sha256,
        "prepared_cms_sha256":bound.prepared_cms_sha256,"scientific_worker_calls":0,
        "source_bound_execution_admitted":false,"online_speedup":null}))
        .map_err(|e| e.to_string())
}

pub(super) fn verify(root: &Path, record: &capsule::Registration) -> Result<(), String> {
    let source = root.join("immutable/source");
    files::reject_ancestor_config(&source)?;
    require(
        canonical_sha(&native::inventory(&source)?)? == record.source_manifest_sha256,
        "SAT target source inventory differs",
    )?;
    capsule::check_roles(root, record)?;
    verify_from(root, record, root)
}

/// Recheck original build invocations from copied receipts. `root` may hold
/// only registration sidecars and receipts; `original` is the historical cwd.
/// A publication checker verifies every archived file against the inventory.
pub(super) fn verify_from(
    root: &Path,
    record: &capsule::Registration,
    original: &Path,
) -> Result<(), String> {
    let receipts = root.join("immutable/build-receipts");
    let build = load(&receipts.join("build-identity.json"))?;
    require(
        build["schema_version"] == 1
            && build["scope"] == capsule::SCOPE
            && build["source_commit"] == record.source_commit
            && build["source_manifest_sha256"] == record.source_manifest_sha256
            && build["hardware"] == record.hardware
            && build["profile"] == "release"
            && build["features"] == "default"
            && build["algorithm_flags"] == "none"
            && build["marked_cms_sha256"] == MARKED_CMS_SHA
            && build["marked_control_sha256"] == MARKED_AUDIT_SHA
            && build["preparation_audit_sha256"] == ORIGINAL_PREPARATION_AUDIT_SHA
            && canonical_sha(&build)? == record.worker_build_identity["build_sha256"],
        "SAT target build identity differs",
    )?;
    let source = original.join("immutable/source");
    let mut controls = Vec::new();
    for (name, pin) in [
        ("marked-audit.json", MARKED_AUDIT_SHA),
        ("preparation-audit.json", ORIGINAL_PREPARATION_AUDIT_SHA),
    ] {
        let path = root.join("immutable/assets/controls").join(name);
        let expected = if name == "marked-audit.json" {
            "PASS_DISCLOSED_CMS_TRANSPORT_PARITY"
        } else {
            "PASS_NATIVE_SOURCE_BOUND_ORDINARY_PREPARATION_AUDIT"
        };
        controls.push(bound_control(&path, pin, expected)?);
    }
    check_prerequisites(&controls[0], &controls[1], &record.preparation)?;
    require(
        record.immutable_files["assets/bin/prepared-cms"]["sha256"] == MARKED_CMS_SHA
            && record.prepared_cms_sha256 == MARKED_CMS_SHA,
        "SAT marked CMS binary differs",
    )?;
    let build_receipt = load(&receipts.join("build.receipt.json"))?;
    let build_env = build_receipt["environment"].clone();
    let environment: Vec<(String, String)> =
        serde_json::from_value(build_env.clone()).map_err(|e| e.to_string())?;
    let map = environment
        .into_iter()
        .collect::<std::collections::BTreeMap<_, _>>();
    require(
        map.len() == 8
            && map.get("IC_SAT_TARGET_SOURCE_MANIFEST_SHA256")
                == Some(&record.source_manifest_sha256)
            && map.get("IC_SAT_TARGET_BUILD_SHA256").map(String::as_str)
                == record.worker_build_identity["build_sha256"].as_str()
            && map.get("CARGO_INCREMENTAL").map(String::as_str) == Some("0")
            && map.get("LC_ALL").map(String::as_str) == Some("C")
            && map.get("CARGO_HOME")
                == Some(&original.join("build-home").to_string_lossy().into_owned())
            && map.get("CARGO_TARGET_DIR")
                == Some(&original.join("build-target").to_string_lossy().into_owned())
            && map.get("RUSTC").map(String::as_str)
                == load(&receipts.join("rustc.receipt.json"))?["program"].as_str()
            && map.contains_key("PATH"),
        "SAT compiled identity or offline build environment differs",
    )?;
    let artifacts = load(&receipts.join("build-artifacts.json"))?;
    require(
        artifacts["schema_version"] == 1
            && artifacts["scope"] == capsule::SCOPE
            && artifacts["worker_sha256"] == record.worker_sha256
            && artifacts["icprog_sha256"] == record.auditor_sha256
            && artifacts["prepared_exporter_sha256"] == record.prepared_exporter_sha256
            && artifacts["prepared_cms_sha256"] == record.prepared_cms_sha256
            && artifacts["build_command_receipt_sha256"]
                == sha256(&read(&receipts.join("build.receipt.json"), 65536)?),
        "SAT copied artifacts differ from frozen build receipt",
    )?;
    for name in ["vendor", "rustc", "cargo", "build", "worker-identity"] {
        let receipt = load(&receipts.join(format!("{name}.receipt.json")))?;
        let bytes = read(&receipts.join(format!("{name}.log")), 16 * 1024 * 1024)?;
        let pin = match name {
            "worker-identity" => &record.worker_sha256,
            "rustc" => build["rustc_sha256"].as_str().ok_or("rustc pin missing")?,
            _ => build["cargo_sha256"].as_str().ok_or("cargo pin missing")?,
        };
        require(
            receipt["exit_code"] == 0
                && receipt["program_sha256_before"] == pin
                && receipt["program_sha256_after"] == pin
                && receipt["log_sha256"] == sha256(&bytes)
                && receipt["wall_ns"]
                    .as_str()
                    .and_then(|s| s.parse::<u128>().ok())
                    .is_some_and(|n| n > 0),
            "SAT frozen tool receipt or output differs",
        )?;
        match name {
            "build" => require(
                receipt["argv"]
                    == json!([
                        "build",
                        "--locked",
                        "--offline",
                        "--release",
                        "--bin",
                        capsule::WORKER,
                        "--bin",
                        "icprog",
                        "--example",
                        "koblitz_pdp_export"
                    ])
                    && receipt["cwd"] == json!(source)
                    && receipt["environment"] == build_env
                    && receipt["program"] == load(&receipts.join("cargo.receipt.json"))?["program"],
                "SAT worker/exporter build invocation differs",
            )?,
            "worker-identity" => require(
                receipt["argv"] == json!(["build-identity"])
                    && receipt["cwd"] == json!(original)
                    && receipt["program"]
                        == json!(original.join("immutable/bin").join(capsule::WORKER))
                    && receipt["environment"]
                        == json!(capsule::environment().into_iter().collect::<Vec<_>>())
                    && load(&receipts.join("worker-identity.log"))? == record.worker_build_identity,
                "SAT compiled worker identity invocation differs",
            )?,
            "vendor" => {
                let env: Vec<(String, String)> =
                    serde_json::from_value(receipt["environment"].clone())
                        .map_err(|e| e.to_string())?;
                require(
                    receipt["argv"]
                        == json!(["vendor", "--locked", "--offline", source.join("vendor")])
                        && receipt["cwd"] == build["checkout_root"]
                        && env.len() == 2
                        && env[0].0 == "PATH"
                        && !env[0].1.is_empty()
                        && env[1] == ("LC_ALL".into(), "C".into())
                        && receipt["program"]
                            == load(&receipts.join("cargo.receipt.json"))?["program"],
                    "SAT offline vendor invocation differs",
                )?;
            }
            "rustc" | "cargo" => require(
                receipt["argv"] == json!(["--version", "--verbose"])
                    && receipt["cwd"] == json!(source)
                    && receipt["environment"] == build_env,
                "SAT tool version invocation differs",
            )?,
            _ => unreachable!(),
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    const MARKED: &[u8] = include_bytes!("../../../research/ic_candidate_tournament_20260915/goal_20260924/native-sat-target-transport-v1/result-v1/original-audit.json");
    const PREP: &[u8] = include_bytes!("../../../research/ic_candidate_tournament_20260915/goal_20260924/native-sat-million-registration-v1/result-v1/original-audit.json");

    #[test]
    fn frozen_control_receipts_bind_the_exact_sat_preparation() {
        assert_eq!(sha256(MARKED), MARKED_AUDIT_SHA);
        assert_eq!(sha256(PREP), ORIGINAL_PREPARATION_AUDIT_SHA);
        let binding = capsule::PreparationBinding {
            capsule: PathBuf::from("/not-executed"),
            execution: PathBuf::from("/not-executed"),
            registration_sha256: "3b768aa18127c2aaafb261770cfafe8eed1c7631e9aaf2a42db211b678c3e470"
                .into(),
            auditor_sha256: "0719685743ffe99b0de308503fafc3451cad97b38c1dfb912599f0960f737bfd"
                .into(),
            producer_sha256: "c7d01efe3f0fa9b2ffe33daaa0d8a59edf1bdff2ab0be00b7c12707e196e451a"
                .into(),
            mathematics_sha256: "d5a428b182ef32f9f7e2b7408befa0bf41c4d194bc3a4bd13cf34d2cd38f6c83"
                .into(),
        };
        check_prerequisites(MARKED, PREP, &binding).unwrap();
        let mut wrong = binding.clone();
        wrong.mathematics_sha256 = "0".repeat(64);
        assert!(check_prerequisites(MARKED, PREP, &wrong).is_err());
        let mut marked: Value = serde_json::from_slice(MARKED).unwrap();
        marked["audited_roles"] = json!(8);
        assert!(
            check_prerequisites(&serde_json::to_vec(&marked).unwrap(), PREP, &binding).is_err()
        );
    }

    #[test]
    fn scientific_freeze_cannot_bypass_exporter_parity_gate() {
        let absent = Path::new("/no-sat-source");
        let request = FreezeRequest {
            root: absent,
            out: absent,
            config: absent,
            preparation_binding: absent,
            marked_cms: absent,
            cargo: absent,
            rustc: absent,
            host_context: absent,
            validation_only: false,
        };
        assert!(freeze(request)
            .unwrap_err()
            .contains("exact-built prepared-exporter parity"));
    }
}
