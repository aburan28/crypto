//! Target-free preparation on one fixed, disclosed educational n17 instance.
//! No imported logs or target. Runtime custody is an external frozen gate.
use super::{
    ic_measurement::{self as clock, Phase, Session, Snapshot},
    koblitz_factor_base_search::FactorBaseSpec,
    koblitz_fast::{FastCurve, FastPoint},
    koblitz_groebner::SolverEngine,
    koblitz_index_calculus::{
        probe_scalar, CollectedRelation, DecompositionStrategy, FactorBaseLogSolver,
        FrobeniusFactorBase, KoblitzCurve, KoblitzIcOptions, LinearAlgebra, PdpOutcome,
        RelationCollector, RelationWorkUnit,
    },
    prepared_sat_control::{self as sat, QueryOutput},
};
use num_traits::ToPrimitive;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{collections::BTreeMap, time::Instant};

const ORDER: u64 = 65587;
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Family {
    MatrixF5,
    Cryptominisat,
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Plan {
    pub schema_version: u32,
    pub question: String,
    pub family: Family,
    pub algorithm_seed: u64,
    pub planned_queries: usize,
}
impl Plan {
    pub fn validate(&self) -> Result<(), String> {
        require(
            self.schema_version == 1
                && self.question == "native-target-free-preparation-n17-v1"
                && (1..=512).contains(&self.planned_queries),
            "outside bounded target-free n17 preparation",
        )
    }
}

/// The sealed native CMS executor supplies original source/model receipts.
/// A callback's output cannot itself attest execution custody.
pub trait SatBackend {
    fn query(&mut self, trial: usize, public: [u64; 2]) -> Result<QueryOutput, String>;
}
/// The frozen worker must durably create each start and completion record.
pub trait Observer {
    fn started(&mut self, trial: usize, scalar: u64) -> Result<(), String>;
    fn completed(&mut self, record: &Value) -> Result<(), String>;
}

fn require(ok: bool, message: &str) -> Result<(), String> {
    if ok {
        Ok(())
    } else {
        Err(message.into())
    }
}
fn ns(start: Instant) -> Result<u64, String> {
    start
        .elapsed()
        .as_nanos()
        .try_into()
        .map_err(|_| "duration overflow".into())
}
fn add_window(
    phases: &mut BTreeMap<String, u64>,
    total: &mut u64,
    window: &Snapshot,
) -> Result<(), String> {
    require(
        window.online_wall_ns.is_none(),
        "preparation entered online timing",
    )?;
    require(
        window.phases_ns.values().flatten().sum::<u64>() == window.observed_wall_ns,
        "phase window does not close",
    )?;
    *total = total
        .checked_add(window.observed_wall_ns)
        .ok_or("duration overflow")?;
    for (&name, &cost) in &window.phases_ns {
        if let Some(cost) = cost {
            *phases.entry(name.into()).or_default() += cost;
        }
    }
    Ok(())
}
struct Context {
    curve: KoblitzCurve,
    base: FrobeniusFactorBase,
    options: KoblitzIcOptions,
    fast: FastCurve,
    generator: FastPoint,
    geometry: Vec<FastPoint>,
}
impl Context {
    fn new(plan: &Plan) -> Result<Self, String> {
        clock::mark(Phase::Setup);
        let curve = KoblitzCurve::new(1, 17).ok_or("n17 curve unavailable")?;
        let fast = FastCurve::new(&curve.curve).ok_or("portable native curve unavailable")?;
        let generator = fast.lift(curve.generator());
        require(
            curve.n == 17
                && curve.a == 1
                && curve.subgroup_order.to_u64() == Some(ORDER)
                && curve.group_order.to_u64() == Some(131174)
                && curve.cofactor.to_u64() == Some(2)
                && curve.lambda.to_u64() == Some(17184)
                && curve.curve.irreducible.low_terms == vec![0, 3]
                && sat::encoded(generator) == Some([43693, 23339]),
            "exact disclosed curve differs",
        )?;
        clock::mark(Phase::FactorBase);
        let base = FactorBaseSpec::StandardSubspace { dimension: 6 }.materialize(&curve)?;
        let geometry = base.points.iter().map(|p| fast.lift(p)).collect::<Vec<_>>();
        require(
            geometry.len() == 63 && geometry.iter().all(|p| !p.infinity && p.x < 64),
            "fresh geometric base differs",
        )?;
        clock::mark(Phase::Precompute);
        let options = KoblitzIcOptions {
            m: 3,
            seed: plan.algorithm_seed,
            max_trials: plan.planned_queries,
            strategy: DecompositionStrategy::Groebner,
            engine: SolverEngine::MatrixF5 { max_degree: 3 },
            node_budget: 8192,
            collection_window: None,
            allow_direct_relation: false,
            linear_algebra: LinearAlgebra::Dense,
            ..Default::default()
        };
        Ok(Self {
            curve,
            base,
            options,
            fast,
            generator,
            geometry,
        })
    }
    fn fixture(&self) -> Value {
        json!({"degree":17,"curve_a":1,"subgroup_order":ORDER,"group_order":131174,
            "cofactor":2,"generator":[43693,23339],"lambda":17184,
            "irreducible":{"degree":17,"low_terms":[0,3]},"targets":[],
            "target_seeds":[],"target_scalar_constructed":false})
    }
    fn lift_model(&self, model: &[bool], public: FastPoint) -> Option<[usize; 3]> {
        if model.len() < 51 {
            return None;
        }
        let xs: [u64; 3] =
            std::array::from_fn(|i| (0..6).fold(0, |v, b| v | ((model[6 * i + b] as u64) << b)));
        let candidates = xs.map(|x| {
            self.geometry
                .iter()
                .enumerate()
                .filter(|(_, p)| p.x == x)
                .map(|(i, _)| i)
                .collect::<Vec<_>>()
        });
        for &a in &candidates[0] {
            for &b in &candidates[1] {
                for &c in &candidates[2] {
                    if self.fast.add(
                        self.fast.add(self.geometry[a], self.geometry[b]),
                        self.geometry[c],
                    ) == public
                    {
                        return Some([a, b, c]);
                    }
                }
            }
        }
        None
    }
}

struct Frontend {
    scalar: u64,
    outcome: String,
    indices: Option<[usize; 3]>,
    evidence: Value,
    fatal: bool,
}
fn sat_frontend(
    context: &Context,
    trial: usize,
    scalar: u64,
    backend: &mut impl SatBackend,
) -> Result<Frontend, String> {
    clock::mark(Phase::Queries);
    let point = context.fast.mul_u64(context.generator, scalar);
    let public = sat::encoded(point).ok_or("nonzero subgroup query became identity")?;
    clock::mark(Phase::Pdp);
    let output = backend.query(trial, public)?;
    clock::mark(Phase::RelationCheck);
    let source_check = sat::validate_manifest(&output.manifest, public).and_then(|_| {
        require(
            output.anf.len() <= 16 * 1024 * 1024
                && output.cnf.len() <= 16 * 1024 * 1024
                && output.native.stdout.len() <= 4 * 1024 * 1024,
            "source/model output exceeded bound",
        )
    });
    if let Err(reason) = source_check {
        return Ok(Frontend {
            scalar,
            outcome: "invalid_model".into(),
            indices: None,
            fatal: true,
            evidence: json!({"public_point":public,"reason":reason,
                "manifest":output.manifest,"source_receipt":output.source_receipt,
                "native":{"exit_code":output.native.exit_code,"timed_out":output.native.timed_out,
                    "stdout_bytes":output.native.stdout.len(),"stdout_sha256":sat::sha256(output.native.stdout.as_bytes()),
                    "receipt":output.native.receipt},
                "anf_sha256":sat::sha256(output.anf.as_bytes()),"anf_bytes":output.anf.len(),
                "cnf_sha256":sat::sha256(output.cnf.as_bytes()),"cnf_bytes":output.cnf.len()}),
        });
    }
    let status = sat::native_status(&output.native);
    let mut indices = None;
    let mut reason = None;
    let outcome = match status {
        "SAT_MODEL" => {
            let checked = (|| {
                let count = output.manifest["exports"]["cryptominisat_xor_dimacs"]["variables"]
                    .as_u64()
                    .filter(|&n| (51..=4096).contains(&n))
                    .ok_or("invalid CNF variable count")? as usize;
                let model = sat::parse_model(&output.native.stdout, count)?;
                sat::verify_anf(&output.anf, &model[..51])?;
                sat::verify_cnf(&output.cnf, &model)?;
                context
                    .lift_model(&model, point)
                    .ok_or_else(|| "source model does not readd".to_string())
            })();
            match checked {
                Ok(witness) => {
                    indices = Some(witness);
                    "witness"
                }
                Err(e) => {
                    reason = Some(e);
                    "invalid_model"
                }
            }
        }
        "SOURCE_UNSAT" => "proved_unsat",
        "TIMEOUT" => "timeout",
        "CONFLICT_BUDGET_INCONCLUSIVE" | "UNKNOWN_INCONCLUSIVE" => "incomplete",
        _ => "transport_failure",
    };
    Ok(Frontend {
        scalar,
        outcome: outcome.into(),
        indices,
        fatal: matches!(outcome, "invalid_model" | "transport_failure"),
        evidence: json!({"public_point":public,"native_status":status,"reason":reason,
            "native":output.native,"manifest":output.manifest,"anf":output.anf,"cnf":output.cnf,
            "source_receipt":output.source_receipt}),
    })
}

/// Real MatrixF5 dispatch, with one collector built and no target/log import.
/// Call only from a separately frozen and accepted one-use native worker.
pub fn prepare_f5(plan: &Plan, observer: &mut impl Observer) -> Result<Value, String> {
    plan.validate()?;
    require(plan.family == Family::MatrixF5, "F5 family differs")?;
    let start = Instant::now();
    let init = Session::begin_strict().map_err(str::to_string)?;
    let context = Context::new(plan)?;
    let collector = RelationCollector::new(&context.curve, &context.base, &context.options)
        .ok_or("F5 collector unavailable")?;
    let dispatch =
        serde_json::to_value(collector.admission_dispatch()).map_err(|e| e.to_string())?;
    run_panel(plan, &context, start, init, dispatch, observer, |trial| {
        let (_, report) = collector.collect_observed(RelationWorkUnit {
            seed: plan.algorithm_seed,
            start: trial as u64,
            count: 1,
        });
        let mut attempts = report.attempts.ok_or("F5 omitted attempt")?;
        require(attempts.len() == 1, "F5 attempt count differs")?;
        let attempt = attempts.remove(0);
        require(
            attempt.trial == trial as u64 && attempt.b == 0,
            "F5 chronology or ordinary query differs",
        )?;
        let outcome = serde_json::to_value(attempt.pdp.outcome).map_err(|e| e.to_string())?;
        let outcome = outcome.as_str().ok_or("invalid F5 status")?.to_string();
        let indices = attempt
            .pdp
            .points
            .map(|v| {
                v.try_into()
                    .map_err(|_| "F5 witness arity differs".to_string())
            })
            .transpose()?;
        Ok(Frontend {
            scalar: attempt.a,
            indices,
            evidence: json!({"solver_stats":attempt.pdp.stats}),
            fatal: matches!(
                attempt.pdp.outcome,
                PdpOutcome::InvalidModel | PdpOutcome::Identity
            ),
            outcome,
        })
    })
}

/// External CMS source/model path. It never substitutes the Rust SAT engine.
pub fn prepare_sat(
    plan: &Plan,
    backend: &mut impl SatBackend,
    observer: &mut impl Observer,
) -> Result<Value, String> {
    plan.validate()?;
    require(plan.family == Family::Cryptominisat, "CMS family differs")?;
    let start = Instant::now();
    let init = Session::begin_strict().map_err(str::to_string)?;
    let context = Context::new(plan)?;
    run_panel(
        plan,
        &context,
        start,
        init,
        json!({"family":"external_cryptominisat","source_model_checks":"anf-and-cnf-and-full-point"}),
        observer,
        |trial| {
            sat_frontend(
                &context,
                trial,
                probe_scalar(plan.algorithm_seed, trial as u64, ORDER),
                backend,
            )
        },
    )
}

#[allow(clippy::too_many_arguments)]
fn run_panel(
    plan: &Plan,
    context: &Context,
    start: Instant,
    init: Session,
    dispatch: Value,
    observer: &mut impl Observer,
    mut frontend: impl FnMut(usize) -> Result<Frontend, String>,
) -> Result<Value, String> {
    let mut matrix = FactorBaseLogSolver::new(&context.curve, &context.base, &context.options)
        .ok_or("relation matrix unavailable")?;
    require(matrix.columns() == 29, "folded columns differ")?;
    let mut phases = BTreeMap::<String, u64>::new();
    let mut window_wall = 0u64;
    let initial = init.finish().map_err(str::to_string)?;
    add_window(&mut phases, &mut window_wall, &initial)?;
    let mut windows = vec![json!({"stage":"reusable_initialization","clock":initial})];
    let mut attempts = Vec::new();
    let mut records = Vec::new();
    let mut stopped = None;
    for trial in 0..plan.planned_queries {
        let scalar = probe_scalar(plan.algorithm_seed, trial as u64, ORDER);
        if let Err(e) = observer.started(trial, scalar) {
            stopped = Some(e);
            break;
        }
        let session = Session::begin_strict().map_err(str::to_string)?;
        let mut result = match frontend(trial) {
            Ok(r) => r,
            Err(e) => Frontend {
                scalar,
                outcome: "transport_failure".into(),
                indices: None,
                evidence: json!({"reason":e}),
                fatal: true,
            },
        };
        result.evidence = json!({"original_claim":{"scalar":result.scalar,"outcome":result.outcome,"indices":result.indices},
            "frontend":result.evidence});
        clock::mark(Phase::RelationCheck);
        if result.scalar != scalar {
            result.fatal = true;
            result.outcome = "invalid_model".into();
            result.indices = None;
        }
        if result.outcome == "witness" {
            if let Some(indices) = result.indices.filter(|v| v.iter().all(|&i| i < 63)) {
                let before = matrix.report().rejected_relations;
                matrix.push(&[CollectedRelation {
                    trial: trial as u64,
                    a: result.scalar,
                    points: indices.to_vec(),
                }]);
                if matrix.report().rejected_relations != before {
                    result.outcome = "invalid_model".into();
                    result.indices = None;
                    result.fatal = true;
                }
            } else {
                result.outcome = "invalid_model".into();
                result.indices = None;
                result.fatal = true;
            }
        } else if result.indices.is_some() {
            result.outcome = "invalid_model".into();
            result.indices = None;
            result.fatal = true;
        }
        let costs = session.finish().map_err(str::to_string)?;
        add_window(&mut phases, &mut window_wall, &costs)?;
        let attempt = json!({"trial":trial,"scalar":result.scalar,"outcome":result.outcome,"indices":result.indices});
        let record = json!({"attempt":attempt,"costs":costs,"evidence":result.evidence,
            "accepted_rows":matrix.relations(),"duplicate_relations":matrix.report().duplicate_relations});
        let progress = observer.completed(&record);
        attempts.push(attempt);
        records.push(record);
        if result.fatal || progress.is_err() {
            stopped = Some(
                progress
                    .err()
                    .unwrap_or_else(|| "fatal frontend outcome".into()),
            );
            break;
        }
    }
    let final_clock = Session::begin_strict().map_err(str::to_string)?;
    clock::mark(Phase::RelationLa);
    let solved = matrix.try_solve();
    let logs = solved
        .as_ref()
        .map(|(table, _)| {
            table
                .columns
                .iter()
                .map(|(_, l)| l.to_u64().ok_or("column log overflow"))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()?;
    let final_costs = final_clock.finish().map_err(str::to_string)?;
    add_window(&mut phases, &mut window_wall, &final_costs)?;
    windows.push(json!({"stage":"final_relation_la_and_group_replay","clock":final_costs}));
    let wall_ns = ns(start)?;
    let controller_ns = wall_ns
        .checked_sub(window_wall)
        .ok_or("overlapping preparation windows")?;
    *phases.entry("setup".into()).or_default() += controller_ns;
    require(
        phases.values().sum::<u64>() == wall_ns,
        "preparation accounting does not close",
    )?;
    let complete = attempts.len() == plan.planned_queries;
    let math = json!({"schema_version":1,"question":"ordinary-preparation-mathematics-n17-v1",
        "plan":{"family":plan.family,"algorithm_seed":plan.algorithm_seed,"planned_queries":plan.planned_queries,
            "query_law":"independent-probe-scalar-rand08-v1"},"fixture":context.fixture(),
        "geometric_base":context.geometry.iter().map(|&p|sat::encoded(p)).collect::<Vec<_>>(),
        "attempts":attempts,"stop":if complete{"panel_complete"}else{"interrupted"},"claimed_column_logs":logs});
    Ok(
        json!({"schema_version":1,"question":"native-target-free-preparation-producer-n17-v1",
        "plan":plan,"panel_complete":complete,"stop_reason":stopped,"mathematical_input":math,
        "frontend_dispatch":dispatch,"query_records":records,"phase_windows":windows,
        "relation_matrix":matrix.matrix_snapshot(),"log_report":matrix.report(),
        "preparation_wall_ns":wall_ns,"exclusive_preparation_phases_ns":phases,
        "controller_setup_ns":controller_ns,"wall_timing_class":"ordinary-host-stage-diagnostic",
        "source_bound_execution_admitted":false,"native_custody_requires_external_audit":true,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,
        "full_goal_complete":false,"online_wall_ns":null,"online_speedup":null,"operation_counts":null}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[derive(Default)]
    struct Observed {
        starts: Vec<(usize, u64)>,
        records: Vec<Value>,
        fail_on: Option<usize>,
    }
    impl Observer for Observed {
        fn started(&mut self, t: usize, a: u64) -> Result<(), String> {
            self.starts.push((t, a));
            Ok(())
        }
        fn completed(&mut self, v: &Value) -> Result<(), String> {
            self.records.push(v.clone());
            if self.fail_on == v["attempt"]["trial"].as_u64().map(|t| t as usize) {
                Err("durable progress failure".into())
            } else {
                Ok(())
            }
        }
    }
    fn plan(count: usize) -> Plan {
        Plan {
            schema_version: 1,
            question: "native-target-free-preparation-n17-v1".into(),
            family: Family::MatrixF5,
            algorithm_seed: 2026093032,
            planned_queries: count,
        }
    }
    fn retained() -> Value {
        serde_json::from_str(include_str!(
        "../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json")).unwrap()
    }
    fn fixture_frontend(old: &Value, trial: usize) -> Frontend {
        let attempt = &old["certificate"]["inputs"]["attempts"][trial];
        Frontend {
            scalar: attempt["scalar"].as_u64().unwrap(),
            outcome: attempt["outcome"].as_str().unwrap().into(),
            indices: serde_json::from_value(attempt["indices"].clone()).unwrap(),
            fatal: false,
            evidence: json!({"kind":"retained-correctness-control-no-native-solver"}),
        }
    }
    fn control(
        plan: &Plan,
        observer: &mut Observed,
        f: impl FnMut(usize) -> Result<Frontend, String>,
    ) -> Value {
        let start = Instant::now();
        let init = Session::begin_strict().unwrap();
        let context = Context::new(plan).unwrap();
        run_panel(
            plan,
            &context,
            start,
            init,
            json!({"kind":"retained-correctness-control-no-native-solver"}),
            observer,
            f,
        )
        .unwrap()
    }
    #[test]
    fn retained_panel_derives_logs_without_import_and_solves_once() {
        let old = retained();
        let mut observer = Observed::default();
        let report = control(&plan(216), &mut observer, |t| Ok(fixture_frontend(&old, t)));
        assert_eq!(observer.starts.len(), 216);
        assert_eq!(observer.records.len(), 216);
        assert_eq!(report["panel_complete"], true);
        assert_eq!(report["log_report"]["solve_attempts"], 1);
        assert_eq!(report["log_report"]["relations"], 61);
        assert_eq!(report["log_report"]["rejected_relations"], 0);
        let expected = old["certificate"]["inputs"]["logs"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v["log"].as_str().unwrap().parse::<u64>().unwrap())
            .collect::<Vec<_>>();
        assert_eq!(
            report["mathematical_input"]["claimed_column_logs"],
            json!(expected)
        );
        let costs = report["exclusive_preparation_phases_ns"]
            .as_object()
            .unwrap()
            .values()
            .map(|v| v.as_u64().unwrap())
            .sum::<u64>();
        assert_eq!(costs, report["preparation_wall_ns"].as_u64().unwrap());
        assert_eq!(report["source_bound_execution_admitted"], false);
        assert_eq!(report["full_goal_complete"], false);
        assert!(report["online_speedup"].is_null());
    }
    #[test]
    fn failed_query_and_durable_progress_failure_preserve_completed_prefix() {
        let old = retained();
        let mut observer = Observed::default();
        let report = control(&plan(216), &mut observer, |t| {
            if t == 7 {
                Err("retained transport failure".into())
            } else {
                Ok(fixture_frontend(&old, t))
            }
        });
        assert_eq!(report["query_records"].as_array().unwrap().len(), 8);
        assert_eq!(
            report["mathematical_input"]["attempts"][7]["outcome"],
            "transport_failure"
        );
        assert_eq!(report["mathematical_input"]["stop"], "interrupted");
        assert_eq!(report["log_report"]["solve_attempts"], 0);
        let mut observer = Observed {
            fail_on: Some(6),
            ..Default::default()
        };
        let report = control(&plan(216), &mut observer, |t| Ok(fixture_frontend(&old, t)));
        assert_eq!(report["query_records"].as_array().unwrap().len(), 7);
        assert_eq!(report["stop_reason"], "durable progress failure");
    }
    #[test]
    fn false_witness_is_rejected_with_original_claim_retained() {
        let old = retained();
        let mut observer = Observed::default();
        let report = control(&plan(216), &mut observer, |t| {
            let mut attempt = fixture_frontend(&old, t);
            if t == 6 {
                attempt.indices = Some([0, 0, 0]);
            }
            Ok(attempt)
        });
        assert_eq!(
            report["query_records"][6]["attempt"]["outcome"],
            "invalid_model"
        );
        assert!(report["query_records"][6]["attempt"]["indices"].is_null());
        assert_eq!(
            report["query_records"][6]["evidence"]["original_claim"]["indices"],
            json!([0, 0, 0])
        );
        assert_eq!(report["log_report"]["rejected_relations"], 1);
        assert_eq!(report["log_report"]["relations"], 0);
    }
    #[test]
    fn full_panel_failure_is_count_complete_but_never_runtime_admitted() {
        let mut observer = Observed::default();
        let report = control(&plan(1), &mut observer, |_| Err("frontend failed".into()));
        assert_eq!(report["mathematical_input"]["stop"], "panel_complete");
        assert!(!report["stop_reason"].is_null());
        assert_eq!(report["source_bound_execution_admitted"], false);
        assert!(report["mathematical_input"]["claimed_column_logs"].is_null());
        let old = retained();
        let mut observer = Observed::default();
        let report = control(&plan(1), &mut observer, |t| {
            let mut a = fixture_frontend(&old, t);
            a.scalar += 1;
            Ok(a)
        });
        assert_eq!(
            report["mathematical_input"]["attempts"][0]["outcome"],
            "invalid_model"
        );
        assert_eq!(report["mathematical_input"]["attempts"][0]["scalar"], 6728);
    }
    #[test]
    fn unentered_online_and_relation_phases_remain_null() {
        let old = retained();
        let mut observer = Observed::default();
        let report = control(&plan(1), &mut observer, |t| Ok(fixture_frontend(&old, t)));
        assert!(report["query_records"][0]["costs"]["phases_ns"]["pdp"].is_null());
        assert!(report["query_records"][0]["costs"]["online_wall_ns"].is_null());
        let mut bad = plan(513);
        assert!(bad.validate().is_err());
        bad = plan(1);
        bad.question = "another curve".into();
        assert!(bad.validate().is_err());
        let mut value = serde_json::to_value(plan(1)).unwrap();
        value["target"] = json!([52411, 72106]);
        assert!(serde_json::from_value::<Plan>(value).is_err());
    }

    struct SatFailure {
        wrong_source: bool,
        calls: usize,
    }
    impl SatBackend for SatFailure {
        fn query(&mut self, _: usize, p: [u64; 2]) -> Result<QueryOutput, String> {
            self.calls += 1;
            Ok(QueryOutput {
                native: sat::NativeOutput {
                    exit_code: None,
                    timed_out: true,
                    stdout: String::new(),
                    receipt: json!({"retained_failure_control":true}),
                },
                manifest: json!({"n":17,"curve_a":1,"ell":6,"m":3,"source_variables":51,"source_equations":50,
                    "representation":"symmetrised_s4","factor_base_basis_bitmasks":["1","2","4","8","16","32"],
                    "target":{"x":if self.wrong_source{p[0]^1}else{p[0]},"y":p[1]}}),
                anf: String::new(),
                cnf: String::new(),
                source_receipt: json!({"retained_failure_control":true}),
            })
        }
    }
    #[test]
    fn sat_timeouts_stay_attempts_and_changed_source_is_retained_and_rejected() {
        let mut p = plan(2);
        p.family = Family::Cryptominisat;
        let mut backend = SatFailure {
            wrong_source: false,
            calls: 0,
        };
        let mut observer = Observed::default();
        let report = prepare_sat(&p, &mut backend, &mut observer).unwrap();
        assert_eq!(backend.calls, 2);
        assert_eq!(report["panel_complete"], true);
        for row in report["query_records"].as_array().unwrap() {
            assert_eq!(row["attempt"]["outcome"], "timeout");
            assert!(row["attempt"]["indices"].is_null());
            assert!(row["costs"]["phases_ns"]["pdp"].as_u64().is_some());
        }
        assert_eq!(report["log_report"]["relations"], 0);
        backend = SatFailure {
            wrong_source: true,
            calls: 0,
        };
        observer = Observed::default();
        let report = prepare_sat(&p, &mut backend, &mut observer).unwrap();
        assert_eq!(backend.calls, 1);
        assert_eq!(
            report["query_records"][0]["attempt"]["outcome"],
            "invalid_model"
        );
        assert_eq!(report["mathematical_input"]["stop"], "interrupted");
        assert!(!report["query_records"][0]["evidence"]["frontend"]["manifest"].is_null());
    }
    struct SourceControl {
        flip: bool,
    }
    impl SatBackend for SourceControl {
        fn query(&mut self, _: usize, p: [u64; 2]) -> Result<QueryOutput, String> {
            assert_eq!(p, [62577, 27783]);
            let result:Value=serde_json::from_str(include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/result.json")).unwrap();
            let mut model = result["cnf_assignment"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_bool().unwrap())
                .collect::<Vec<_>>();
            if self.flip {
                model[0] = !model[0];
            }
            let literals = model
                .iter()
                .enumerate()
                .map(|(i, &b)| {
                    if b {
                        (i as i64 + 1).to_string()
                    } else {
                        (-(i as i64 + 1)).to_string()
                    }
                })
                .collect::<Vec<_>>()
                .join(" ");
            Ok(QueryOutput{native:sat::NativeOutput{exit_code:Some(10),timed_out:false,stdout:format!("s SATISFIABLE\nv {literals} 0\n"),receipt:json!({"retained_source_control":true})},
                manifest:serde_json::from_str(include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/manifest.json")).unwrap(),
                anf:include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/instance.anf").into(),
                cnf:include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/instance.xor.cnf").into(),
                source_receipt:json!({"retained_source_control":true})})
        }
    }
    #[test]
    fn retained_sat_source_model_lifts_and_altered_assignment_does_not() {
        let session = Session::begin_strict().unwrap();
        let context = Context::new(&plan(1)).unwrap();
        // Exact historical disclosed control, not the new ordinary query law.
        let scalar = (32326 + 42888 * 24886u64) % ORDER;
        let good = sat_frontend(&context, 1, scalar, &mut SourceControl { flip: false }).unwrap();
        assert_eq!(good.outcome, "witness");
        assert!(good.indices.is_some());
        let bad = sat_frontend(&context, 1, scalar, &mut SourceControl { flip: true }).unwrap();
        assert_eq!(bad.outcome, "invalid_model");
        assert!(bad.indices.is_none());
        assert!(bad.fatal);
        assert!(session.finish().unwrap().online_wall_ns.is_none());
    }
}
