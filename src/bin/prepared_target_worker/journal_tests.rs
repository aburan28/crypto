use super::*;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{PdpAttempt, PdpOutcome, PdpSolverStats};
use std::time::{SystemTime, UNIX_EPOCH};
fn temp(label: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!(
        "ic-target-journal-{label}-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&p).unwrap();
    p.canonicalize().unwrap()
}
fn binding(family: Family) -> Binding {
    Binding {
        schema_version: 1,
        question: QUESTION.into(),
        registration_sha256: "1".repeat(64),
        worker_sha256: "2".repeat(64),
        family,
        target: [52411, 72106],
        algorithm_seed: 2026100310,
        max_queries: 3,
    }
}
fn start(b: &Binding, trial: usize) -> Value {
    let [a, c] = b.coefficients()[trial];
    match b.family {
        Family::MatrixF5 => json!({"trial":trial,"a":a,"b":c,"target":b.target,
            "query_rule":"seeded-sample-aG-plus-bQ","decomposition_started":false}),
        Family::Cryptominisat => json!({"trial":trial,"a":a,"b":c,"public_query":[1,2],
            "outcome":"query_planned","source_model_valid":false,"witness_indices":null,
            "candidate_scalar":null,"backend_called":false}),
    }
}
fn completion(b: &Binding, trial: usize) -> Value {
    let [a, c] = b.coefficients()[trial];
    match b.family {
        Family::MatrixF5 => {
            json!({"target":b.target,"record_kind":"completed-pdp-before-scalar-recovery",
            "attempt":QueryAttempt {trial:trial as u64,a,b:c,pdp:PdpAttempt {outcome:PdpOutcome::ProvedUnsat,points:None,stats:PdpSolverStats::None}}})
        }
        Family::Cryptominisat => {
            json!({"trial":trial,"a":a,"b":c,"outcome":"timeout","backend_called":true})
        }
    }
}
#[test]
fn durable_start_survives_drop_as_unadmitted_pending_prefix() {
    let p = temp("pending");
    let b = binding(Family::MatrixF5);
    let mut writer = Journal::create(&p, b.clone()).unwrap();
    writer.start(&start(&b, 0)).unwrap();
    drop(writer);
    let receipt = inspect(&p, &b.registration_sha256).unwrap();
    assert_eq!(receipt["pending_start"]["trial"], 0);
    assert_eq!(receipt["completed_attempts"].as_array().unwrap().len(), 0);
    assert_eq!(receipt["source_bound_execution_admitted"], false);
    assert!(Journal::create(&p, b).is_err());
    fs::remove_dir_all(p).unwrap();
}
#[test]
fn both_record_shapes_preserve_complete_prefix_with_no_runtime_admission() {
    for family in [Family::MatrixF5, Family::Cryptominisat] {
        let p = temp("complete");
        let b = binding(family);
        let mut writer = Journal::create(&p, b.clone()).unwrap();
        for trial in 0..3 {
            writer.start(&start(&b, trial)).unwrap();
            writer.complete(&completion(&b, trial)).unwrap();
        }
        let receipt = inspect(&p, &b.registration_sha256).unwrap();
        assert_eq!(receipt["completed_attempts"].as_array().unwrap().len(), 3);
        assert!(receipt["pending_start"].is_null());
        assert_eq!(receipt["searches_executed_by_inspection"], 0);
        assert_eq!(receipt["source_bound_execution_admitted"], false);
        assert!(receipt["online_speedup"].is_null());
        fs::remove_dir_all(p).unwrap();
    }
}
#[test]
fn failed_write_poisoning_prevents_overwrite_and_resume() {
    let p = temp("failed");
    let b = binding(Family::MatrixF5);
    let mut writer = Journal::create(&p, b.clone()).unwrap();
    native::create(
        &p.join("attempts/start-000.json"),
        b"retained earlier failed bytes",
    )
    .unwrap();
    assert!(writer.start(&start(&b, 0)).is_err());
    assert!(writer.start(&start(&b, 1)).is_err());
    assert!(writer.complete(&completion(&b, 0)).is_err());
    assert_eq!(
        native::read(&p.join("attempts/start-000.json"), 65536).unwrap(),
        b"retained earlier failed bytes"
    );
    assert!(inspect(&p, &b.registration_sha256).is_err());
    fs::remove_dir_all(p).unwrap();
}
#[test]
fn wrong_target_coefficients_family_and_completion_are_rejected() {
    for mutation in 0..4 {
        let p = temp("changed");
        let b = binding(Family::MatrixF5);
        let mut writer = Journal::create(&p, b.clone()).unwrap();
        if mutation < 3 {
            let mut row = start(&b, 0);
            match mutation {
                0 => row["target"] = json!([1, 2]),
                1 => row["a"] = json!(1),
                _ => row = start(&binding(Family::Cryptominisat), 0),
            };
            assert!(writer.start(&row).is_err());
        } else {
            writer.start(&start(&b, 0)).unwrap();
            assert!(writer.complete(&completion(&b, 1)).is_err());
            assert!(writer.complete(&completion(&b, 0)).is_err());
        }
        fs::remove_dir_all(p).unwrap();
    }
}
#[test]
fn inspector_rejects_orphans_gaps_changed_links_and_wrong_seal() {
    for mutation in 0..4 {
        let p = temp("inspect");
        let b = binding(Family::MatrixF5);
        let mut writer = Journal::create(&p, b.clone()).unwrap();
        writer.start(&start(&b, 0)).unwrap();
        writer.complete(&completion(&b, 0)).unwrap();
        match mutation {
            0 => fs::remove_file(p.join("attempts/start-000.json")).unwrap(),
            1 => native::create(&p.join("attempts/start-002.json"), b"{}").unwrap(),
            2 => {
                let path = p.join("attempts/completed-000.json");
                let mut row = native::load(&path).unwrap();
                row["start_sha256"] = json!("3".repeat(64));
                fs::write(path, serde_json::to_vec(&row).unwrap()).unwrap();
            }
            _ => {
                assert!(inspect(&p, &"4".repeat(64)).is_err());
                fs::remove_dir_all(p).unwrap();
                continue;
            }
        }
        assert!(inspect(&p, &b.registration_sha256).is_err());
        fs::remove_dir_all(p).unwrap();
    }
}
#[test]
fn binding_and_record_limits_are_fail_closed() {
    let mut b = binding(Family::MatrixF5);
    b.max_queries = 9;
    assert!(b.validate().is_err());
    b.max_queries = 3;
    b.target[0] = 1 << 17;
    assert!(b.validate().is_err());
    b.target = [52411, 72106];
    b.registration_sha256 = "Z".repeat(64);
    assert!(b.validate().is_err());
}
#[test]
fn abrupt_record_helper() {
    let Some(root) = std::env::var_os("IC_TARGET_RECORD_CONTROL") else {
        return;
    };
    let p = PathBuf::from(root);
    let b = binding(Family::MatrixF5);
    let mut writer = Journal::create(&p, b.clone()).unwrap();
    writer.start(&start(&b, 0)).unwrap();
    println!("DURABLE_START_READY");
    use std::io::Write;
    std::io::stdout().flush().unwrap();
    // The parent watchdog kills this control process. It runs no solver.
    std::thread::sleep(std::time::Duration::from_secs(30));
    panic!("control watchdog must terminate helper");
}
#[cfg(unix)]
#[test]
fn watchdog_killed_writer_retains_start_and_original_drain_receipt() {
    let p = temp("watchdog");
    let own = std::env::current_exe().unwrap();
    let args = vec![
        "--exact".into(),
        // The journal is also compiled below icprog::target_control. Test names
        // omit the crate prefix, so retain the actual module path in both bins.
        format!(
            "{}::abrupt_record_helper",
            module_path!().split_once("::").unwrap().1
        ),
        "--nocapture".into(),
    ];
    let env = vec![(
        "IC_TARGET_RECORD_CONTROL".into(),
        p.to_string_lossy().into_owned(),
    )];
    let ledger = p.join("control-pids");
    let result = native::measured_child_request(native::ChildRequest {
        program: &own,
        args: &args,
        cwd: &p,
        stem: &p.join("control"),
        deadline_ms: 2000,
        ledger: &ledger,
        helper: false,
        input: None,
        environment: &env,
    })
    .unwrap();
    assert!(result.timed_out);
    assert!(result.stdout.contains("DURABLE_START_READY"));
    native::drain_ledger(&ledger).unwrap();
    let receipt = inspect(&p, &binding(Family::MatrixF5).registration_sha256).unwrap();
    assert_eq!(receipt["pending_start"]["trial"], 0);
    assert_eq!(receipt["completed_attempts"].as_array().unwrap().len(), 0);
    assert_eq!(receipt["source_bound_execution_admitted"], false);
    // Preserve actual control-process receipts separately from scientific data.
    if let Some(dest) = std::env::var_os("IC_TARGET_RECORD_RECEIPTS") {
        let dest = PathBuf::from(dest);
        fs::create_dir(&dest).unwrap();
        for name in ["control.receipt.json", "control.stdout", "control.stderr"] {
            let bytes = native::read(&p.join(name), 16 * 1024 * 1024).unwrap();
            native::create(&dest.join(name), &bytes).unwrap();
        }
        native::save(&dest.join("prefix.json"), &receipt).unwrap();
        native::save(
            &dest.join("attempt-files.json"),
            &native::inventory(&p.join("attempts")).unwrap(),
        )
        .unwrap();
        let records = dest.join("attempts");
        fs::create_dir(&records).unwrap();
        for entry in fs::read_dir(p.join("attempts")).unwrap() {
            let entry = entry.unwrap();
            native::create(
                &records.join(entry.file_name()),
                &native::read(&entry.path(), 8 * 1024 * 1024).unwrap(),
            )
            .unwrap();
        }
        native::save(
            &dest.join("CONTROL.json"),
            &json!({"schema_version":1,
            "kind":"durable-file-only-watchdog-control","cryptographic_searches_executed":0,
            "scientific_registration":false,"source_bound_execution_admitted":false,
            "full_goal_complete":false,"online_speedup":null}),
        )
        .unwrap();
    }
    fs::remove_dir_all(p).unwrap();
}
