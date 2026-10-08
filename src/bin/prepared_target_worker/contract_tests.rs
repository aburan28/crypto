use super::*;
use std::{
    fs,
    time::{SystemTime, UNIX_EPOCH},
};
fn cfg() -> Config {
    Config {
        schema_version: 1,
        question: SCOPE.into(),
        target: [52411, 72106],
        plan: F5TargetPlan {
            schema_version: 1,
            question: "prepared-public-target-n17-f5-v1".into(),
            algorithm_seed: 2026100310,
            max_queries: 8,
        },
        worker_timeout_ms: 900000,
    }
}
fn record() -> Registration {
    Registration {
        schema_version: 1,
        scope: SCOPE.into(),
        source_commit: "1".repeat(40),
        hardware: json!({"os":"macos","architecture":"aarch64"}),
        immutable_files: json!({}),
        source_manifest_sha256: "2".repeat(64),
        config_sha256: "3".repeat(64),
        host_context_sha256: "6".repeat(64),
        worker_sha256: "4".repeat(64),
        auditor_sha256: "5".repeat(64),
        worker_build_identity: json!({"schema_version":1,"source_manifest_sha256":"2".repeat(64),
            "build_sha256":"6".repeat(64),"target_os":"macos","target_arch":"aarch64"}),
        runtime_environment: environment(),
        preparation: PreparationBinding {
            capsule: "/synthetic-no-execution/preparation".into(),
            execution: "/synthetic-no-execution/results".into(),
            registration_sha256: "7".repeat(64),
            auditor_sha256: "8".repeat(64),
            producer_sha256: "9".repeat(64),
            mathematics_sha256: "a".repeat(64),
        },
        validation_only: false,
    }
}
#[test]
fn config_rejects_unbounded_and_known_scalar_fields() {
    assert!(cfg().validate().is_ok());
    let mut bad = cfg();
    bad.plan.max_queries = 9;
    assert!(bad.validate().is_err());
    bad = cfg();
    bad.target[0] = 1 << 17;
    assert!(bad.validate().is_err());
    let mut value = json!(cfg());
    value["known_scalar"] = json!(24886);
    assert!(serde_json::from_value::<Config>(value).is_err());
    value = json!(cfg());
    value["family"] = json!("cryptominisat");
    assert!(serde_json::from_value::<Config>(value).is_err());
}
#[test]
fn old_capsule_scope_and_mismatched_compiled_identity_are_rejected() {
    let mut r = record();
    assert!(pins(&r).is_ok());
    r.scope = preparation::SCOPE.into();
    assert!(pins(&r).is_err());
    let expected = record().worker_build_identity;
    assert!(require_identity(&json!({}), &expected).is_err());
    let mut changed = expected.clone();
    changed["build_sha256"] = json!("b".repeat(64));
    assert!(require_identity(&changed, &expected).is_err());
    changed = expected.clone();
    changed["source_manifest_sha256"] = Value::Null;
    assert!(require_identity(&changed, &changed).is_err());
}
#[test]
fn claim_rejects_validation_capsule_wrong_seal_and_wrong_execution() {
    let p = std::env::temp_dir().join(format!(
        "ic-target-claim-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&p).unwrap();
    let p = p.canonicalize().unwrap();
    let execution = p.join("execution");
    fs::create_dir(&execution).unwrap();
    let mut r = record();
    let seal = "c".repeat(64);
    native::save(
        &p.join("consumed.json"),
        &json!(Claim {
            schema_version: 1,
            scope: SCOPE.into(),
            registration_sha256: seal.clone(),
            execution: execution.clone(),
            worker_sha256: r.worker_sha256.clone()
        }),
    )
    .unwrap();
    assert!(claim(&p, &execution, &r, &seal).is_ok());
    r.validation_only = true;
    assert!(claim(&p, &execution, &r, &seal).is_err());
    r.validation_only = false;
    assert!(claim(&p, &p, &r, &seal).is_err());
    assert!(claim(&p, &execution, &r, &"d".repeat(64)).is_err());
    fs::remove_dir_all(p).unwrap();
}
#[test]
fn preparation_receipt_pins_reject_missing_native_admission_and_changed_input() {
    // Synthetic receipt syntax only; no capsule/native execution is admitted.
    let producer =
        json!({"mathematical_input":{"kind":"synthetic-no-mathematical-or-native-admission"}});
    let bytes = serde_json::to_vec(&producer).unwrap();
    let mut prep = record().preparation;
    prep.producer_sha256 = sha256(&bytes);
    prep.mathematics_sha256 = canonical_sha(&producer["mathematical_input"]).unwrap();
    let receipt = json!({"schema_version":1,"status":"PASS_NATIVE_SOURCE_BOUND_ORDINARY_PREPARATION_AUDIT",
        "registration_sha256":prep.registration_sha256,"checker_sha256":prep.auditor_sha256,"producer_sha256":prep.producer_sha256,
        "source_bound_execution_admitted":true,"mathematical_preparation_complete":true,
        "mathematics":{"input_sha256":prep.mathematics_sha256,"declared_family":"matrix_f5","rank":29,"usable_points":62,"panel_complete":true}});
    assert!(preparation_receipt(&prep, &receipt, &producer, &bytes).is_ok());
    for field in [
        "source_bound_execution_admitted",
        "mathematical_preparation_complete",
    ] {
        let mut bad = receipt.clone();
        bad[field] = json!(false);
        assert!(preparation_receipt(&prep, &bad, &producer, &bytes).is_err());
    }
    let mut bad = receipt.clone();
    bad["mathematics"]["declared_family"] = json!("cryptominisat");
    assert!(preparation_receipt(&prep, &bad, &producer, &bytes).is_err());
    assert!(preparation_receipt(&prep, &receipt, &producer, b"changed").is_err());
}
