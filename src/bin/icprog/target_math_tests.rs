//! Historical bytes + synthetic new envelopes; no worker, F5 or SAT execution.
use super::*;
use std::{
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};
const PREP: &str = include_str!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json");
const TARGET: &str = include_str!("../../../research/ic_candidate_tournament_20260915/goal_20260924/native-f5-control-registration-v1/result-v1/data/execution/worker.stdout");
fn retained() -> (Value, Config, Value) {
    let old: Value = parse(PREP.as_bytes()).unwrap();
    let original: Value = serde_json::from_str(TARGET).unwrap(); // Old elapsed seconds are historical floats, not new input schema.
    let inputs = &old["certificate"]["inputs"];
    let geometry = inputs["base"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| {
            p.as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_str().unwrap().parse::<u64>().unwrap())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let prep = json!({"schema_version":1,"question":"ordinary-preparation-mathematics-n17-v1",
        "plan":{"family":"matrix_f5","algorithm_seed":2026093032u64,"planned_queries":216,"query_law":"independent-probe-scalar-rand08-v1"},
        "fixture":inputs["fixture"],"geometric_base":geometry,"attempts":inputs["attempts"],"stop":"panel_complete",
        "claimed_column_logs":inputs["logs"].as_array().unwrap().iter().map(|v|v["log"].as_str().unwrap().parse::<u64>().unwrap()).collect::<Vec<_>>()});
    let cfg = Config {
        schema_version: 1,
        question: SCOPE.into(),
        target: [52411, 72106],
        plan: Plan {
            schema_version: 1,
            question: "prepared-public-target-n17-f5-v1".into(),
            algorithm_seed: 2026100310,
            max_queries: 8,
        },
        worker_timeout_ms: 900000,
    };
    let mut attempts = original["solutions"][0]["attempts"].clone();
    for row in attempts.as_array_mut().unwrap() {
        if let Some(indices) = row["pdp"]["points"].as_array_mut() {
            for index in indices {
                let p = &original["factor_base"][index.as_u64().unwrap() as usize];
                let point = p
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_str().unwrap().parse::<u64>().unwrap())
                    .collect::<Vec<_>>();
                *index = json!(geometry.iter().position(|v| v == &point).unwrap());
            }
        }
    }
    let mut phases = Map::new();
    let mut online = Map::new();
    for name in PHASES {
        phases.insert(
            name.into(),
            if name == "setup" || (ONLINE.contains(&name) && name != "rho_solve") {
                json!(1)
            } else {
                Value::Null
            },
        );
    }
    for name in ONLINE {
        online.insert(
            name.into(),
            if name != "rho_solve" {
                json!(1)
            } else {
                Value::Null
            },
        );
    }
    let producer = json!({"schema_version":1,"question":"prepared-public-target-n17-f5-producer-v1","plan":cfg.plan,
        "target":cfg.target,"target_count":1,"preparation_family":"matrix_f5","preparation_math_sha256":canonical_sha(&prep).unwrap(),
        "actual_usable_base_points":62,"geometric_base_points":63,"effective_columns":29,
        "dispatch":{"strategy":"Groebner","field_kernel":"portable","pair_table":false,"query_rule":"seeded-sample","summands":3,
            "direct_collision":false,"engine":{"MatrixF5":{"max_degree":3}},"node_budget":8192,"collapse_negation":true,
            "collapse_projected_orbits":false,"attempt_record_contract":"start-and-completed-pdp-before-recovery-v1"},
        "reusable_validation_ns":0,"solver_setup_ns":0,
        "report":{"trials":3,"claimed_scalar":"24886","relation":{"a":45669,"b":63073,"points":attempts[2]["pdp"]["points"]},"attempts":attempts},
        "verified_scalar":24886,"independent_recovery_certificate":{"check":"independent-shift-and-add-scalar-replay",
            "generator":[43693,23339],"scalar":24886,"target":cfg.target,"inside_online_interval":true},"failure":null,"interruption":null,
        "costs":{"schema_version":1,"phases_ns":phases,"observed_wall_ns":6,"online_phases_ns":online,"online_wall_ns":5},
        "exclusive_ic_online_phases_ns":{"target_query":1,"target_pdp":1,"target_relation_check":1,"target_descent":1,"target_recovery_check":1},
        "timing_class":"ordinary-host-diagnostic","candidate_id":null,"workload_id":null,"run_id":null,"source_bound_execution_admitted":false,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,"full_goal_complete":false,"online_speedup":null,"operation_counts":null});
    (prep, cfg, producer)
}
fn temp() -> PathBuf {
    let p = std::env::temp_dir().join(format!(
        "ic-target-math-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&p).unwrap();
    p.canonicalize().unwrap()
}
fn records(p: &Path, cfg: &Config, report: &Value) {
    fs::create_dir(p.join("attempts")).unwrap();
    let root = p.join("attempts");
    let binding = json!({"schema_version":1,"question":RECORD_SCOPE,"registration_sha256":"1".repeat(64),"worker_sha256":"2".repeat(64),
        "family":"matrix_f5","target":cfg.target,"algorithm_seed":cfg.plan.algorithm_seed,"max_queries":cfg.plan.max_queries});
    native::save(&root.join("binding.json"), &binding).unwrap();
    for (trial, attempt) in report["report"]["attempts"]
        .as_array()
        .unwrap()
        .iter()
        .enumerate()
    {
        let a = &attempt["a"];
        let b = &attempt["b"];
        let header = json!({"trial":trial,"a":a,"b":b,"target":cfg.target});
        let start = json!({"schema_version":1,"binding_sha256":canonical_sha(&binding).unwrap(),"header":header,"start_sha256":null,
            "body":{"trial":trial,"a":a,"b":b,"target":cfg.target,"query_rule":"seeded-sample-aG-plus-bQ","decomposition_started":false}});
        native::save(&root.join(format!("start-{trial:03}.json")), &start).unwrap();
        native::save(&root.join(format!("completed-{trial:03}.json")), &json!({"schema_version":1,"binding_sha256":canonical_sha(&binding).unwrap(),
            "header":header,"start_sha256":canonical_sha(&start).unwrap(),"body":{"target":cfg.target,"attempt":attempt,"record_kind":"completed-pdp-before-scalar-recovery"}})).unwrap();
    }
}
fn unadmitted(result: &Value) {
    for k in [
        "source_bound_execution_admitted",
        "fresh_paired_qualification",
        "headline_eligible",
        "promotion_eligible",
        "full_goal_complete",
    ] {
        assert_eq!(result[k], false);
    }
    assert!(result["online_wall_ns"].is_null());
    assert!(result["online_speedup"].is_null());
}
#[test]
fn retained_witness_and_negatives_verify_without_runtime_or_timing_admission() {
    let (prep, cfg, report) = retained();
    let result = verify(&prep, &cfg, &report).unwrap();
    assert_eq!(result["verified_scalar"], 24886);
    assert_eq!(result["outcome_mix"]["proved_unsat"], 2);
    assert_eq!(result["outcome_mix"]["witness"], 1);
    assert_eq!(result["solver_calls_by_auditor"], 0);
    unadmitted(&result);
    if let Some(path) = std::env::var_os("IC_TARGET_MATH_CONTROLS") {
        let p = PathBuf::from(path);
        fs::create_dir(&p).unwrap();
        let execution = p.join("execution");
        fs::create_dir(&execution).unwrap();
        records(&execution, &cfg, &report);
        native::save(&p.join("preparation.json"), &prep).unwrap();
        native::save(&p.join("config.json"), &json!(cfg)).unwrap();
        native::save(&p.join("producer.json"), &report).unwrap();
        native::save(&p.join("CONTROL.json"),&json!({"schema_version":1,"kind":"historical-mathematics-synthetic-producer-and-journal-control",
            "synthetic_timing":true,"synthetic_dispatch_envelope":true,"synthetic_registration_seal":true,"new_solver_calls":0,
            "historical_preparation_sha256":sha256(PREP.as_bytes()),"historical_target_report_sha256":sha256(TARGET.as_bytes()),
            "source_bound_execution_admitted":false,"full_goal_complete":false,"online_speedup":null})).unwrap();
        run(
            &p.join("preparation.json"),
            &p.join("config.json"),
            &p.join("producer.json"),
            &execution,
            &"1".repeat(64),
            &"2".repeat(64),
            &p.join("audit.json"),
        )
        .unwrap();
    }
}
#[test]
fn scalar_relation_witness_and_certificate_tampering_are_rejected() {
    let (prep, cfg, report) = retained();
    for (pointer, value) in [
        ("/verified_scalar", json!(24885)),
        ("/report/claimed_scalar", json!("24885")),
        ("/report/relation/a", json!(1)),
        ("/report/attempts/2/pdp/points/0", json!(0)),
        (
            "/independent_recovery_certificate/inside_online_interval",
            json!(false),
        ),
        ("/independent_recovery_certificate/target", json!([1, 2])),
    ] {
        let mut bad = report.clone();
        *bad.pointer_mut(pointer).unwrap() = value;
        assert!(verify(&prep, &cfg, &bad).is_err(), "accepted {pointer}");
    }
}
#[test]
fn false_negative_with_true_decomposition_is_independently_rejected() {
    let (prep, mut cfg, mut report) = retained();
    cfg.plan.max_queries = 3;
    report["plan"] = json!(cfg.plan);
    report["report"]["attempts"][2]["pdp"]["outcome"] = json!("proved_unsat");
    report["report"]["attempts"][2]["pdp"]["points"] = Value::Null;
    let error = verify(&prep, &cfg, &report).unwrap_err();
    assert!(error.contains("false geometric negative"), "{error}");
}
#[test]
fn inconclusive_budget_exhaustion_is_retained_and_never_a_win() {
    let (prep, mut cfg, report) = retained();
    cfg.plan.max_queries = 3;
    for outcome in ["incomplete", "unsupported"] {
        let mut r = report.clone();
        r["plan"] = json!(cfg.plan);
        let pdp = &mut r["report"]["attempts"][2]["pdp"];
        pdp["outcome"] = json!(outcome);
        pdp["points"] = Value::Null;
        pdp["stats"]["stats"]["exhausted"] = json!(true);
        pdp["stats"]["stats"]["unsupported"] = json!(outcome == "unsupported");
        r["report"]["claimed_scalar"] = Value::Null;
        r["report"]["relation"] = Value::Null;
        r["verified_scalar"] = Value::Null;
        r["independent_recovery_certificate"] = Value::Null;
        r["failure"] = json!("query budget exhausted without scalar recovery");
        let result = verify(&prep, &cfg, &r).unwrap();
        assert_eq!(result["verified_recovery"], false);
        assert_eq!(result["outcome_mix"][outcome], 1);
        unadmitted(&result);
    }
}
#[test]
fn wrong_engine_limits_dispatch_and_identity_are_rejected() {
    let (prep, cfg, report) = retained();
    for (pointer, v) in [
        ("/report/attempts/0/pdp/stats/stats/reductions", json!(8193)),
        (
            "/report/attempts/0/pdp/stats/stats/max_degree_built",
            json!(4),
        ),
        (
            "/report/attempts/0/pdp/stats/engine",
            json!({"MatrixF4":{"max_degree":3}}),
        ),
        ("/dispatch/direct_collision", json!(true)),
        ("/dispatch/pair_table", json!(true)),
        ("/report/attempts/0/pdp/outcome", json!("identity")),
    ] {
        let mut bad = report.clone();
        *bad.pointer_mut(pointer).unwrap() = v;
        assert!(verify(&prep, &cfg, &bad).is_err(), "accepted {pointer}");
    }
}
#[test]
fn query_law_stop_count_preparation_and_schema_changes_are_rejected() {
    let (prep, cfg, report) = retained();
    for (pointer, v) in [
        ("/report/attempts/0/a", json!(1)),
        ("/report/trials", json!(2)),
        ("/target_count", json!(2)),
        ("/target", json!([0, 1])),
        ("/preparation_math_sha256", json!("a".repeat(64))),
        ("/source_bound_execution_admitted", json!(true)),
        ("/interruption", json!({"stage":"completion_record"})),
    ] {
        let mut bad = report.clone();
        *bad.pointer_mut(pointer).unwrap() = v;
        assert!(verify(&prep, &cfg, &bad).is_err());
    }
    let mut bad = report.clone();
    bad.as_object_mut().unwrap().remove("failure");
    assert!(verify(&prep, &cfg, &bad).is_err());
    bad = report.clone();
    bad["known_scalar"] = json!(24886);
    assert!(verify(&prep, &cfg, &bad).is_err());
    let mut bad_cfg = json!(cfg);
    bad_cfg["known_scalar"] = json!(24886);
    assert!(serde_json::from_value::<Config>(bad_cfg).is_err());
    let mut bad_prep = prep.clone();
    bad_prep["claimed_column_logs"][0] = json!(0);
    assert!(verify(&bad_prep, &cfg, &report).is_err());
}
#[test]
fn no_query_may_follow_recovery_and_early_budget_stop_is_rejected() {
    let (prep, cfg, report) = retained();
    let mut bad = report.clone();
    let extra = bad["report"]["attempts"][0].clone();
    bad["report"]["attempts"]
        .as_array_mut()
        .unwrap()
        .push(extra);
    bad["report"]["trials"] = json!(4);
    assert!(verify(&prep, &cfg, &bad).unwrap_err().contains("continued"));
    bad = report.clone();
    bad["report"]["attempts"].as_array_mut().unwrap().pop();
    bad["report"]["trials"] = json!(2);
    bad["report"]["claimed_scalar"] = Value::Null;
    bad["report"]["relation"] = Value::Null;
    bad["verified_scalar"] = Value::Null;
    bad["independent_recovery_certificate"] = Value::Null;
    bad["failure"] = json!("query budget exhausted without scalar recovery");
    assert!(verify(&prep, &cfg, &bad).is_err());
}
#[test]
fn clock_mapping_overflow_extra_phases_and_reusable_work_are_rejected() {
    let (_, _, r) = retained();
    for (pointer, v) in [
        ("/costs/online_wall_ns", json!(4)),
        ("/costs/observed_wall_ns", json!(5)),
        ("/exclusive_ic_online_phases_ns/target_pdp", json!(2)),
        ("/costs/phases_ns/setup", json!(u64::MAX)),
        ("/costs/phases_ns/factor_base", json!(1)),
        ("/costs/online_phases_ns/rho_solve", json!(1)),
    ] {
        let mut bad = r.clone();
        *bad.pointer_mut(pointer).unwrap() = v;
        assert!(verify_costs(&bad).is_err(), "accepted {pointer}");
    }
    let mut bad = r.clone();
    bad["costs"]["online_phases_ns"]["extra"] = json!(0);
    assert!(verify_costs(&bad).is_err());
}
#[test]
fn missing_phase_stays_unknown_even_with_correct_scalar() {
    let (prep, cfg, mut r) = retained();
    r["costs"]["phases_ns"]["target_descent"] = Value::Null;
    r["costs"]["online_phases_ns"]["target_descent"] = Value::Null;
    r["exclusive_ic_online_phases_ns"]["target_descent"] = Value::Null;
    r["costs"]["online_wall_ns"] = json!(4);
    r["costs"]["observed_wall_ns"] = json!(5);
    let result = verify(&prep, &cfg, &r).unwrap();
    assert_eq!(result["verified_recovery"], true);
    assert!(result["timing"]["five_phase_recorded_total_ns"].is_null());
    assert_eq!(result["timing"]["missing_costs"], json!(["target_descent"]));
    unadmitted(&result);
}
#[test]
fn duplicate_keys_nested_floats_and_missing_schema_are_rejected() {
    for bytes in [
        br#"{"x":1,"x":2}"#.as_slice(),
        br#"{"nested":{"x":1,"x":2}}"#.as_slice(),
        br#"{"x":1.0}"#.as_slice(),
        br#"[1,{"x":1e2}]"#.as_slice(),
    ] {
        assert!(parse(bytes).is_err());
    }
    let ordinary =
        br#"{"mathematical_input":{"rank":29},"log_report":{"collection_seconds":0.007025581}}"#;
    assert_eq!(
        parse_ordinary_producer(ordinary).unwrap()["log_report"]["collection_seconds"],
        json!(0.007025581)
    );
    assert!(parse_ordinary_producer(
        br#"{"log_report":{"collection_seconds":0.1,"collection_seconds":0.2}}"#
    )
    .is_err());
    let (_, cfg, report) = retained();
    assert_eq!(
        parse(&serde_json::to_vec(&report).unwrap()).unwrap(),
        report
    );
    let mut c = json!(cfg);
    c.as_object_mut().unwrap().remove("target");
    assert!(serde_json::from_value::<Config>(c).is_err());
}
#[test]
fn original_records_binding_links_pending_orphans_and_body_changes_are_rejected() {
    let (_, cfg, r) = retained();
    for mutation in 0..6 {
        let p = temp();
        records(&p, &cfg, &r);
        assert!(verify_records(&p, &cfg, &r, &"1".repeat(64), &"2".repeat(64)).is_ok());
        match mutation {
            0 => assert!(verify_records(&p, &cfg, &r, &"3".repeat(64), &"2".repeat(64)).is_err()),
            1 => assert!(verify_records(&p, &cfg, &r, &"1".repeat(64), &"3".repeat(64)).is_err()),
            2 => {
                fs::remove_file(p.join("attempts/completed-002.json")).unwrap();
                assert!(verify_records(&p, &cfg, &r, &"1".repeat(64), &"2".repeat(64)).is_err());
            }
            3 => {
                native::save(&p.join("attempts/start-003.json"), &json!({})).unwrap();
                assert!(verify_records(&p, &cfg, &r, &"1".repeat(64), &"2".repeat(64)).is_err());
            }
            _ => {
                let path = p.join("attempts/completed-000.json");
                let mut v = native::load(&path).unwrap();
                if mutation == 4 {
                    v["start_sha256"] = json!("0".repeat(64));
                } else {
                    v["body"]["attempt"]["a"] = json!(1);
                }
                fs::write(path, serde_json::to_vec(&v).unwrap()).unwrap();
                assert!(verify_records(&p, &cfg, &r, &"1".repeat(64), &"2".repeat(64)).is_err());
            }
        }
        fs::remove_dir_all(p).unwrap();
    }
}
#[test]
fn failed_audit_retains_receipt_and_original_input_and_cannot_write_into_evidence() {
    let (prep, cfg, mut r) = retained();
    let p = temp();
    let execution = p.join("execution");
    fs::create_dir(&execution).unwrap();
    records(&execution, &cfg, &r);
    r["verified_scalar"] = json!(24885);
    native::save(&p.join("prep.json"), &prep).unwrap();
    native::save(&p.join("cfg.json"), &json!(cfg)).unwrap();
    native::save(&p.join("producer.json"), &r).unwrap();
    let before = native::read(&p.join("producer.json"), 16 * 1024 * 1024).unwrap();
    let bad_out = execution.join("audit.json");
    assert!(run(
        &p.join("prep.json"),
        &p.join("cfg.json"),
        &p.join("producer.json"),
        &execution,
        &"1".repeat(64),
        &"2".repeat(64),
        &bad_out
    )
    .is_err());
    assert!(!bad_out.exists());
    let out = p.join("failed-audit.json");
    assert!(run(
        &p.join("prep.json"),
        &p.join("cfg.json"),
        &p.join("producer.json"),
        &execution,
        &"1".repeat(64),
        &"2".repeat(64),
        &out
    )
    .is_err());
    let receipt = native::load(&out).unwrap();
    assert_eq!(receipt["status"], "TARGET_MATHEMATICS_NOT_VERIFIED");
    unadmitted(&receipt);
    assert_eq!(
        native::read(&p.join("producer.json"), 16 * 1024 * 1024).unwrap(),
        before
    );
    fs::remove_dir_all(p).unwrap();
}
