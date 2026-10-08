//! Synthetic negative controls only: no F5/SAT search or runtime admission.
use super::*;
use std::{
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};
struct Directory(PathBuf);
impl Directory {
    fn new(label: &str) -> Self {
        let path = std::env::temp_dir().join(format!(
            "target-controller-{label}-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&path).unwrap();
        Self(path.canonicalize().unwrap())
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}
fn cfg() -> capsule::Config {
    serde_json::from_value(json!({"schema_version":1,"question":capsule::SCOPE,"target":[52411,72106],
        "plan":{"schema_version":1,"question":"prepared-public-target-n17-f5-v1","algorithm_seed":2026100310,"max_queries":8},
        "worker_timeout_ms":900000})).unwrap()
}
fn record() -> capsule::Registration {
    serde_json::from_value(json!({"schema_version":1,"scope":capsule::SCOPE,"source_commit":"1".repeat(40),
        "hardware":{"os":"macos","architecture":"aarch64"},"immutable_files":{},"source_manifest_sha256":"2".repeat(64),
        "config_sha256":"3".repeat(64),"host_context_sha256":"4".repeat(64),"worker_sha256":"5".repeat(64),"auditor_sha256":"6".repeat(64),
        "worker_build_identity":{"schema_version":1,"source_manifest_sha256":"2".repeat(64),"build_sha256":"7".repeat(64),
        "target_os":"macos","target_arch":"aarch64"},"runtime_environment":capsule::environment(),
        "preparation":{"capsule":"/synthetic-no-execution/preparation","execution":"/synthetic-no-execution/results",
            "registration_sha256":"8".repeat(64),"auditor_sha256":"9".repeat(64),"producer_sha256":"a".repeat(64),"mathematics_sha256":"b".repeat(64)},
        "validation_only":false})).unwrap()
}
fn rewrite(path: &Path, value: &Value) {
    fs::write(path, serde_json::to_vec(value).unwrap()).unwrap();
}
#[test]
fn consumed_claim_and_exposure_are_one_use_and_not_freshness() {
    let dir = Directory::new("claim");
    let root = dir.0.join("capsule");
    let execution = dir.0.join("execution");
    fs::create_dir(&root).unwrap();
    fs::create_dir(&execution).unwrap();
    let record = record();
    let seal = "c".repeat(64);
    consume(&root, &execution, &record, &cfg(), &seal).unwrap();
    capsule::claim(&root, &execution, &record, &seal).unwrap();
    let original = read(&root.join("consumed.json"), 65536).unwrap();
    assert!(consume(&root, &execution, &record, &cfg(), &seal).is_err());
    assert_eq!(read(&root.join("consumed.json"), 65536).unwrap(), original);
    assert_eq!(
        load(&execution.join("target-exposure.json")).unwrap(),
        exposure(&cfg(), &record, &seal)
    );
    assert_eq!(
        exposure(&cfg(), &record, &seal)["fresh_target_certified"],
        false
    );
    assert!(capsule::claim(&root, &execution, &record, &"d".repeat(64)).is_err());
    assert!(capsule::claim(&root, &dir.0, &record, &seal).is_err());
}
#[test]
fn validation_development_and_unsealed_capsules_cannot_dispatch() {
    let mut record = record();
    let seal = "c".repeat(64);
    admissible(&record, &seal).unwrap();
    record.validation_only = true;
    assert!(admissible(&record, &seal).is_err());
    record.validation_only = false;
    assert!(admissible(&record, "c").is_err());
    let dir = Directory::new("dev");
    assert!(frozen_self(&dir.0, &record).is_err());
    assert!(original_preparation(&record.preparation).is_err());
}
#[test]
fn interrupted_prefix_without_files_is_not_admitted_or_assumed_no_search() {
    let dir = Directory::new("empty");
    let checked = prefix(&dir.0, &cfg(), &record(), &"c".repeat(64)).unwrap();
    assert_eq!(checked["status"], "NO_DURABLE_TARGET_ATTEMPT_RECORD");
    assert_eq!(checked["source_bound_execution_admitted"], false);
    assert_eq!(checked["online_speedup"], Value::Null);
    assert!(checked.get("solver_calls").is_none());
}
#[test]
fn prefix_worker_and_target_binding_cannot_be_swapped() {
    let dir = Directory::new("binding");
    let manifest = record();
    let seal = "c".repeat(64);
    let _journal =
        journal::Journal::create(&dir.0, cfg().binding(&seal, &manifest.worker_sha256)).unwrap();
    prefix(&dir.0, &cfg(), &manifest, &seal).unwrap();
    let mut changed = manifest;
    changed.worker_sha256 = "d".repeat(64);
    assert!(prefix(&dir.0, &cfg(), &changed, &seal).is_err());
    let mut changed_cfg = cfg();
    changed_cfg.target[0] ^= 1;
    assert!(prefix(&dir.0, &changed_cfg, &record(), &seal).is_err());
}
#[test]
fn whole_interval_requires_nonoverlapping_components_and_checked_sum() {
    let mut producer =
        json!({"costs":{"observed_wall_ns":50},"reusable_validation_ns":10,"solver_setup_ns":20});
    let prep = json!({"child_wall_ns":30});
    let outer = json!({"child_wall_ns":120});
    assert_eq!(interval(&producer, &prep, &outer).unwrap(), 10);
    assert!(interval(&producer, &prep, &json!({"child_wall_ns":109})).is_err());
    producer["solver_setup_ns"] = Value::Null;
    assert!(interval(&producer, &prep, &outer).is_err());
    producer["solver_setup_ns"] = json!(u64::MAX);
    assert!(interval(&producer, &prep, &outer).is_err());
}
#[test]
fn verified_budget_failure_is_not_a_successful_exit_or_win() {
    status(&json!({"exit_code":0,"timed_out":false}), true).unwrap();
    status(&json!({"exit_code":1,"timed_out":false}), false).unwrap();
    assert!(status(&json!({"exit_code":0,"timed_out":false}), false).is_err());
    assert!(status(&json!({"exit_code":1,"timed_out":false}), true).is_err());
    assert!(status(&json!({"exit_code":0,"timed_out":true}), true).is_err());
    assert!(status(&json!({"exit_code":null,"timed_out":true}), false).is_err());
}
#[test]
fn ledger_requires_the_exact_completed_original_role() {
    let dir = Directory::new("ledger");
    let path = dir.0.join("pids");
    fs::write(&path, b"start 123\ndone 123\n").unwrap();
    ledger(&path, 123).unwrap();
    assert!(ledger(&path, 124).is_err());
    for bytes in [
        "start 123\n",
        "start 123\ndone 124\n",
        "start 123\ndone 123\nstart 456\n",
    ] {
        fs::write(&path, bytes).unwrap();
        assert!(ledger(&path, 123).is_err());
    }
}
#[test]
fn audit_output_cannot_mutate_original_capsule_or_execution() {
    let dir = Directory::new("output");
    let root = dir.0.join("capsule");
    let execution = dir.0.join("execution");
    fs::create_dir(&root).unwrap();
    fs::create_dir(&execution).unwrap();
    output_outside(&execution, &root, &dir.0.join("audit.json")).unwrap();
    assert!(output_outside(&execution, &root, &root.join("audit.json")).is_err());
    assert!(output_outside(&execution, &root, &execution.join("audit.json")).is_err());
}
fn terminal() -> Value {
    json!({"schema_version":1,"scope":capsule::SCOPE,"registration_sha256":"c".repeat(64),"controller_sha256":"d".repeat(64),
        "controller_unchanged":true,"controller_read_error":null,"source_gate_passed":true,"source_gate_error":null,
        "worker_drain_passed":true,"worker_drain_error":null,"preparation_auditor_drain_passed":true,"preparation_auditor_drain_error":null,
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,"headline_eligible":false,
        "promotion_eligible":false,"full_goal_complete":false,"online_speedup":null})
}
#[test]
fn controller_source_and_nested_drain_gates_cannot_be_omitted() {
    terminal_binding(&terminal(), &"c".repeat(64), &"d".repeat(64)).unwrap();
    for key in [
        "controller_unchanged",
        "source_gate_passed",
        "worker_drain_passed",
        "preparation_auditor_drain_passed",
    ] {
        let mut changed = terminal();
        changed[key] = json!(false);
        assert!(terminal_binding(&changed, &"c".repeat(64), &"d".repeat(64)).is_err());
        changed.as_object_mut().unwrap().remove(key);
        assert!(terminal_binding(&changed, &"c".repeat(64), &"d".repeat(64)).is_err());
    }
    for key in [
        "source_bound_execution_admitted",
        "fresh_paired_qualification",
        "headline_eligible",
        "promotion_eligible",
        "full_goal_complete",
    ] {
        let mut changed = terminal();
        changed[key] = json!(true);
        assert!(terminal_binding(&changed, &"c".repeat(64), &"d".repeat(64)).is_err());
    }
}
#[test]
fn original_worker_and_nested_preparation_receipts_reject_substitution() {
    let dir = Directory::new("receipts");
    let seal = "c".repeat(64);
    let pin = "d".repeat(64);
    for (name, args, cwd) in [
        ("worker", worker_args(&dir.0, &dir.0, &seal), dir.0.clone()),
        (
            "preparation-auditor",
            preparation_args(&dir.0, &dir.0, &seal, &dir.0.join("audit.json")),
            dir.0.clone(),
        ),
    ] {
        let stem = dir.0.join(name);
        native::create(&stem.with_extension("stdout"), b"synthetic stdout only").unwrap();
        native::create(&stem.with_extension("stderr"), b"").unwrap();
        let argv = std::iter::once("/synthetic-no-execution/bin/role".to_string())
            .chain(args)
            .collect::<Vec<_>>();
        let receipt = json!({"argv":argv,"cwd":cwd,"pid":123,"deadline_ms":900000,"child_wall_ns":120,
            "executable_sha256_before":pin,"executable_sha256_after":pin,"stdout_sha256":sha256(b"synthetic stdout only"),
            "stdout_bytes":21,"stderr_sha256":sha256(b""),"stderr_bytes":0,"output_limit":false,
            "process_group_drain_requested":true,"process_group_drain_confirmed":true,"environment":capsule::environment()});
        let verify = |r: &Value| {
            ordinary_control::verify_receipt(
                &cwd,
                &stem,
                r,
                &argv,
                &pin,
                900000,
                &json!(capsule::environment()),
            )
        };
        assert_eq!(verify(&receipt).unwrap(), 123);
        for (key, value) in [
            ("executable_sha256_after", json!("e".repeat(64))),
            ("cwd", json!(dir.0.join("other"))),
            ("argv", json!(["different-role"])),
            ("environment", json!({"RAYON_NUM_THREADS":"2","LC_ALL":"C"})),
            ("deadline_ms", json!(900001)),
            ("process_group_drain_confirmed", json!(false)),
        ] {
            let mut changed = receipt.clone();
            changed[key] = value;
            assert!(verify(&changed).is_err());
        }
        fs::write(stem.with_extension("stdout"), b"modified bytes").unwrap();
        assert!(verify(&receipt).is_err());
    }
}
fn synthetic_build() -> (Directory, capsule::Registration) {
    let dir = Directory::new("build");
    let root = &dir.0;
    let source = root.join("immutable/source");
    let receipts = root.join("immutable/build-receipts");
    fs::create_dir_all(&source).unwrap();
    fs::create_dir_all(&receipts).unwrap();
    let mut record = record();
    let build = json!({"schema_version":1,"scope":capsule::SCOPE,"source_commit":record.source_commit,
        "checkout_root":"/synthetic-no-execution/source","source_manifest_sha256":record.source_manifest_sha256,
        "cargo_sha256":"a".repeat(64),"rustc_sha256":"b".repeat(64),"profile":"release","features":"default",
        "hardware":record.hardware,"algorithm_flags":"none","native_archive_sha256":null});
    record.worker_build_identity["build_sha256"] = json!(canonical_sha(&build).unwrap());
    save(&receipts.join("build-identity.json"), &build).unwrap();
    let env = vec![
        ("PATH", "/synthetic-no-execution/bin".to_string()),
        ("LC_ALL", "C".into()),
        ("RUSTC", "/synthetic-no-execution/bin/rustc".into()),
        (
            "CARGO_HOME",
            root.join("build-home").to_string_lossy().into_owned(),
        ),
        (
            "CARGO_TARGET_DIR",
            root.join("build-target").to_string_lossy().into_owned(),
        ),
        ("CARGO_INCREMENTAL", "0".into()),
        (
            "IC_TARGET_SOURCE_MANIFEST_SHA256",
            record.source_manifest_sha256.clone(),
        ),
        (
            "IC_TARGET_BUILD_SHA256",
            record.worker_build_identity["build_sha256"]
                .as_str()
                .unwrap()
                .into(),
        ),
    ];
    for name in ["vendor", "rustc", "cargo", "build", "worker-identity"] {
        let (program, args, cwd, environment, pin, bytes) = match name {
            "vendor" => (
                json!("/synthetic-no-execution/bin/cargo"),
                json!(["vendor", "--locked", "--offline", source.join("vendor")]),
                build["checkout_root"].clone(),
                json!([env[0].clone(), env[1].clone()]),
                build["cargo_sha256"].clone(),
                b"synthetic vendor control only".to_vec(),
            ),
            "worker-identity" => (
                json!(root.join("immutable/bin").join(capsule::WORKER)),
                json!(["build-identity"]),
                json!(root),
                json!(capsule::environment().into_iter().collect::<Vec<_>>()),
                json!(record.worker_sha256),
                serde_json::to_vec(&record.worker_build_identity).unwrap(),
            ),
            _ => (
                json!(format!(
                    "/synthetic-no-execution/bin/{}",
                    if name == "rustc" { "rustc" } else { "cargo" }
                )),
                if name == "build" {
                    json!([
                        "build",
                        "--locked",
                        "--offline",
                        "--release",
                        "--bin",
                        capsule::WORKER,
                        "--bin",
                        "icprog"
                    ])
                } else {
                    json!(["--version", "--verbose"])
                },
                json!(source),
                json!(env),
                build[if name == "rustc" {
                    "rustc_sha256"
                } else {
                    "cargo_sha256"
                }]
                .clone(),
                b"synthetic tool control only".to_vec(),
            ),
        };
        native::create(&receipts.join(format!("{name}.log")), &bytes).unwrap();
        save(&receipts.join(format!("{name}.receipt.json")),&json!({"program":program,"argv":args,"cwd":cwd,"environment":environment,
            "exit_code":0,"wall_ns":"10","program_sha256_before":pin,"program_sha256_after":pin,"log_sha256":sha256(&bytes)})).unwrap();
    }
    (dir, record)
}
#[test]
fn frozen_build_receipt_contract_rejects_changed_flags_tools_and_identity() {
    let (dir, record) = synthetic_build();
    target_build::verify(&dir.0, &record).unwrap(); // Receipt control only; this is not a capsule or executable.
    let path = dir.0.join("immutable/build-receipts/build.receipt.json");
    let original = load(&path).unwrap();
    for (key, value) in [
        ("argv", json!(["build", "--release"])),
        ("cwd", json!(dir.0)),
        ("exit_code", json!(1)),
        ("program_sha256_after", json!("f".repeat(64))),
        ("environment", json!([])),
        ("wall_ns", json!("0")),
    ] {
        let mut changed = original.clone();
        changed[key] = value;
        rewrite(&path, &changed);
        assert!(
            target_build::verify(&dir.0, &record).is_err(),
            "accepted changed {key}"
        );
    }
    rewrite(&path, &original);
    let path = dir.0.join("immutable/build-receipts/worker-identity.log");
    let original = fs::read(&path).unwrap();
    fs::write(&path, b"{}").unwrap();
    assert!(target_build::verify(&dir.0, &record).is_err());
    fs::write(&path, &original).unwrap();
    let path = dir.0.join("immutable/build-receipts/build-identity.json");
    let mut build = load(&path).unwrap();
    build["features"] = json!("alternate");
    rewrite(&path, &build);
    assert!(target_build::verify(&dir.0, &record).is_err());
    assert!(frozen_self(&dir.0, &record).is_err());
}

#[test]
fn host_bytes_and_external_seal_are_bound_without_granting_execution() {
    let dir = Directory::new("capsule-data");
    let root = &dir.0;
    fs::create_dir_all(root.join("immutable/source")).unwrap();
    fs::create_dir_all(root.join("immutable/bin")).unwrap();
    let marker = b"synthetic non-executable capsule control only\n";
    // Git cannot retain an empty source directory. Keep the data-only marker
    // in the inventory so this rejection fixture survives a clean checkout.
    native::create(&root.join("immutable/source/CONTROL.txt"), marker).unwrap();
    native::create(&root.join("immutable/bin").join(capsule::WORKER), marker).unwrap();
    native::create(&root.join("immutable/bin/icprog"), marker).unwrap();
    save(&root.join("config.json"), &json!(cfg())).unwrap();
    save(&root.join("host-context.json"),&json!({"scope":"synthetic-non-executable-control-only","os":"macos",
        "architecture":"aarch64","physical_hardware":false,"calibrated_performance_environment":false})).unwrap();
    let mut record = record();
    record.source_manifest_sha256 =
        canonical_sha(&native::inventory(&root.join("immutable/source")).unwrap()).unwrap();
    record.worker_build_identity["source_manifest_sha256"] = json!(record.source_manifest_sha256);
    record.worker_sha256 = sha256(marker);
    record.auditor_sha256 = sha256(marker);
    record.config_sha256 = sha256(&read(&root.join("config.json"), 65536).unwrap());
    record.host_context_sha256 = sha256(&read(&root.join("host-context.json"), 65536).unwrap());
    record.immutable_files = native::inventory(&root.join("immutable")).unwrap();
    record.validation_only = true;
    save(&root.join("registration.json"), &json!(record)).unwrap();
    let seal = canonical_sha(&json!(record)).unwrap();
    save(
        &root.join("seal.json"),
        &json!({"registration_sha256":seal}),
    )
    .unwrap();
    capsule::check_capsule(root, &seal).unwrap();
    assert!(admissible(&record, &seal).is_err());
    assert!(frozen_self(root, &record).is_err());
    assert!(capsule::check_capsule(root, &"c".repeat(64)).is_err());
    if let Some(out) = std::env::var_os("IC_TARGET_CONTROLLER_CONTROL_OUT") {
        let out = PathBuf::from(out);
        native::copy_tree(root, &out).unwrap();
        save(&out.join("CONTROL.json"),&json!({"schema_version":1,"scope":"synthetic-capsule-data-and-dispatch-rejection-only",
            "registration_sha256":seal,"scientific_worker_calls":0,"preparation_panels_executed":0,
            "source_bound_execution_admitted":false,"full_goal_complete":false,"online_speedup":null})).unwrap();
    }
    fs::write(root.join("host-context.json"), b"{}").unwrap();
    assert!(capsule::check_capsule(root, &seal).is_err());
}
