//! Retained-data controls. No F5 search, old worker or registration executes.
use super::boundary_tests::{plan, retained_math, retained_report, CLOCK};
use super::*;
use std::{cell::RefCell, rc::Rc};

#[derive(Default)]
struct Records {
    events: Rc<RefCell<Vec<String>>>,
    starts: Vec<Value>,
    completions: Vec<Value>,
    fail_start: Option<usize>,
    fail_complete: Option<usize>,
}
impl F5TargetObserver for Records {
    fn started(&mut self, row: &Value) -> Result<(), String> {
        let trial = row["trial"].as_u64().unwrap() as usize;
        self.events.borrow_mut().push(format!("start:{trial}"));
        self.starts.push(row.clone());
        if self.fail_start == Some(trial) {
            Err("synthetic start write failure".into())
        } else {
            Ok(())
        }
    }
    fn completed(&mut self, row: &Value) -> Result<(), String> {
        let trial = row["attempt"]["trial"].as_u64().unwrap() as usize;
        self.events.borrow_mut().push(format!("complete:{trial}"));
        self.completions.push(row.clone());
        if self.fail_complete == Some(trial) {
            Err("synthetic completion write failure".into())
        } else {
            Ok(())
        }
    }
}

fn replay(
    state: &PreparedN17Target,
    cap: usize,
    records: &mut Records,
    script: IndividualLogReport,
) -> Value {
    let opts = options(2026100310, cap);
    let solver =
        IndividualLogSolver::new(&state.curve, &state.base, &state.table, &opts, None).unwrap();
    let attempts = script.attempts.unwrap();
    let events = records.events.clone();
    let mut bridge = F5RecordBridge {
        public: [52411, 72106],
        observer: records,
    };
    let mut calls = 0;
    state
        .run_observed_result(
            [52411, 72106],
            &plan(cap),
            json!({"kind":"retained-pdp-data-control-no-f5-search"}),
            0,
            |q| {
                solver.solve_n17_durable_control(q, &mut bridge, &mut |a, b| {
                    let original = &attempts[calls];
                    assert_eq!((a, b), (original.a, original.b));
                    assert_eq!(events.borrow().last().unwrap(), &format!("start:{calls}"));
                    events.borrow_mut().push(format!("probe:{calls}"));
                    calls += 1;
                    clock::mark(Phase::TargetPdp);
                    original.pdp.clone()
                })
            },
        )
        .unwrap()
}
fn closure(result: &Value) {
    let phases = result["exclusive_ic_online_phases_ns"].as_object().unwrap();
    assert_eq!(phases.len(), 5);
    assert_eq!(
        phases.values().filter_map(Value::as_u64).sum::<u64>(),
        result["costs"]["online_wall_ns"].as_u64().unwrap()
    );
    for flag in [
        "source_bound_execution_admitted",
        "fresh_paired_qualification",
        "headline_eligible",
        "promotion_eligible",
        "full_goal_complete",
    ] {
        assert_eq!(result[flag], false);
    }
    for field in ["candidate_id", "workload_id", "run_id", "online_speedup"] {
        assert!(result[field].is_null());
    }
}

#[test]
fn records_start_probe_completion_before_recovery_and_preserves_seeded_law() {
    let _guard = CLOCK.lock().unwrap();
    let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
    let original = retained_report(&state);
    let mut records = Records::default();
    let result = replay(&state, 8, &mut records, original.clone());
    assert_eq!(
        *records.events.borrow(),
        [
            "start:0",
            "probe:0",
            "complete:0",
            "start:1",
            "probe:1",
            "complete:1",
            "start:2",
            "probe:2",
            "complete:2"
        ]
    );
    assert_eq!(result["report"]["attempts"], json!(original.attempts));
    assert_eq!(result["verified_scalar"], 24886);
    assert!(result["failure"].is_null());
    assert!(result["interruption"].is_null());
    assert_eq!(records.starts[2]["target"], json!([52411, 72106]));
    assert_eq!(records.starts[2]["decomposition_started"], false);
    assert_eq!(
        result["independent_recovery_certificate"]["inside_online_interval"],
        true
    );
    closure(&result);
}

#[test]
fn failed_start_preserves_prefix_and_never_decomposes_that_query() {
    let _guard = CLOCK.lock().unwrap();
    let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
    let mut records = Records {
        fail_start: Some(1),
        ..Default::default()
    };
    let result = replay(&state, 8, &mut records, retained_report(&state));
    assert_eq!(
        *records.events.borrow(),
        ["start:0", "probe:0", "complete:0", "start:1"]
    );
    assert_eq!(result["report"]["trials"], 2);
    assert_eq!(result["report"]["attempts"].as_array().unwrap().len(), 1);
    assert_eq!(result["interruption"]["stage"], "start_record");
    assert_eq!(result["interruption"]["trial"], 1);
    assert_eq!(
        result["interruption"]["coefficients"],
        json!([43768, 33418])
    );
    assert!(result["verified_scalar"].is_null());
    assert!(result["exclusive_ic_online_phases_ns"]["target_recovery_check"].is_null());
    closure(&result);
}

#[test]
fn failed_witness_completion_never_recovers_or_admits_the_scalar() {
    let _guard = CLOCK.lock().unwrap();
    let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
    let mut records = Records {
        fail_complete: Some(2),
        ..Default::default()
    };
    let result = replay(&state, 8, &mut records, retained_report(&state));
    assert_eq!(records.events.borrow().len(), 9);
    assert_eq!(result["report"]["trials"], 3);
    assert_eq!(result["report"]["attempts"].as_array().unwrap().len(), 3);
    assert_eq!(result["interruption"]["stage"], "completion_record");
    assert!(result["report"]["claimed_scalar"].is_null());
    assert!(result["report"]["relation"].is_null());
    assert!(result["verified_scalar"].is_null());
    assert!(result["exclusive_ic_online_phases_ns"]["target_descent"].is_null());
    assert!(result["exclusive_ic_online_phases_ns"]["target_recovery_check"].is_null());
    closure(&result);
}

#[test]
fn budget_exhaustion_retains_every_unsuccessful_query() {
    let _guard = CLOCK.lock().unwrap();
    let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
    let mut records = Records::default();
    let result = replay(&state, 2, &mut records, retained_report(&state));
    assert_eq!(records.starts.len(), 2);
    assert_eq!(records.completions.len(), 2);
    assert_eq!(result["report"]["trials"], 2);
    assert_eq!(result["report"]["attempts"].as_array().unwrap().len(), 2);
    assert!(result["verified_scalar"].is_null());
    assert!(result["failure"]
        .as_str()
        .unwrap()
        .contains("budget exhausted"));
    assert!(result["interruption"].is_null());
    closure(&result);
}

#[test]
fn invalid_public_point_never_reaches_recording_or_solver() {
    let _guard = CLOCK.lock().unwrap();
    let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
    let mut records = Records::default();
    let result = state
        .solve_f5_recorded([1 << 17, 1], &plan(8), &mut records)
        .unwrap();
    assert!(records.events.borrow().is_empty());
    assert!(result["report"].is_null());
    assert!(result["verified_scalar"].is_null());
    assert!(result["failure"].is_string());
    closure(&result);
}

#[test]
fn internal_recording_entry_rejects_unbounded_or_different_dispatch() {
    let _guard = CLOCK.lock().unwrap();
    let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
    for change in 0..4 {
        let mut opts = options(2026100310, 8);
        match change {
            0 => opts.max_trials = 9,
            1 => opts.allow_direct_relation = true,
            2 => opts.node_budget = 8193,
            _ => opts.engine = SolverEngine::MatrixF5 { max_degree: 4 },
        }
        let solver =
            IndividualLogSolver::new(&state.curve, &state.base, &state.table, &opts, None).unwrap();
        let mut records = Records::default();
        let mut bridge = F5RecordBridge {
            public: [52411, 72106],
            observer: &mut records,
        };
        let stopped = solver
            .solve_n17_durable_control(&point([52411, 72106]), &mut bridge, &mut |_, _| {
                panic!("source gate must precede decomposition")
            })
            .unwrap_err();
        assert_eq!(stopped.stage, "source_gate");
        assert_eq!(stopped.report.trials, 0);
        assert!(records.events.borrow().is_empty());
    }
}

#[test]
fn invalid_model_is_retained_and_stops_before_next_query() {
    let _guard = CLOCK.lock().unwrap();
    let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
    let mut script = retained_report(&state);
    script.attempts.as_mut().unwrap()[0].pdp.outcome = PdpOutcome::InvalidModel;
    let mut records = Records::default();
    let result = replay(&state, 8, &mut records, script);
    assert_eq!(records.events.borrow().len(), 3);
    assert_eq!(result["interruption"]["stage"], "invalid_model");
    assert_eq!(result["report"]["trials"], 1);
    assert!(result["verified_scalar"].is_null());
    closure(&result);
}

#[test]
fn interrupted_transcript_cannot_smuggle_a_claimed_scalar() {
    let _guard = CLOCK.lock().unwrap();
    let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
    let result = state
        .run_observed_result(
            [52411, 72106],
            &plan(8),
            json!({"kind":"corrupt-data-control"}),
            0,
            |_| {
                Err(IndividualLogInterruption {
                    report: Box::new(retained_report(&state)),
                    stage: "completion_record",
                    trial: Some(2),
                    coefficients: Some([45669, 63073]),
                    reason: "synthetic write failure".into(),
                })
            },
        )
        .unwrap();
    assert!(result["verified_scalar"].is_null());
    assert!(result["failure"]
        .as_str()
        .unwrap()
        .contains("cannot retain"));
    closure(&result);
}

#[test]
fn diagnostic_sampling_preserves_observed_and_unobserved_reports() {
    let _guard = CLOCK.lock().unwrap();
    let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
    let original = retained_report(&state);
    let opts = options(2026100310, 8);
    let solver =
        IndividualLogSolver::new(&state.curve, &state.base, &state.table, &opts, None).unwrap();
    for observe in [false, true] {
        let mut calls = 0;
        let result = solver
            .solve_n17_diagnostic_control(&point([52411, 72106]), observe, &mut |a, b| {
                let old = &original.attempts.as_ref().unwrap()[calls];
                assert_eq!((a, b), (old.a, old.b));
                calls += 1;
                old.pdp.clone()
            })
            .unwrap();
        assert_eq!(result.trials, 3);
        assert_eq!(result.log, original.log);
        assert_eq!(json!(result.relation), json!(original.relation));
        assert_eq!(result.attempts.is_some(), observe);
        if observe {
            assert_eq!(result.attempts, original.attempts);
        }
    }
}

#[test]
fn interrupted_stage_and_prefix_shape_are_checked() {
    let _guard = CLOCK.lock().unwrap();
    let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
    let mut report = retained_report(&state);
    report.log = None;
    report.relation = None;
    let mut stopped = IndividualLogInterruption {
        report: Box::new(report),
        stage: "completion_record",
        trial: Some(2),
        coefficients: Some([45669, 63073]),
        reason: "synthetic failure".into(),
    };
    assert!(interruption_shape(&stopped).is_ok());
    stopped.trial = Some(1);
    assert!(interruption_shape(&stopped).is_err());
    stopped.trial = Some(2);
    stopped.stage = "start_record";
    assert!(interruption_shape(&stopped).is_err());
    stopped.report.attempts.as_mut().unwrap().pop();
    assert!(interruption_shape(&stopped).is_ok());
    stopped.stage = "unknown";
    assert!(interruption_shape(&stopped).is_err());
}
