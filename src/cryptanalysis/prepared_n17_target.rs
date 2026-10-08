//! Prepared public-point boundary for the disclosed educational n17 fixture.
//!
//! Mathematical table checks do not attest the preparation's source or yield.
//! A frozen outer controller must bind its original preparation audit, this
//! source, target exposure, limits and execution receipts before admission.
use super::{
    ic_measurement::{self as clock, Phase, Session},
    koblitz_factor_base_search::FactorBaseSpec,
    koblitz_groebner::SolverEngine,
    koblitz_index_calculus::{
        DecompositionStrategy, FactorBaseLogSolver, FactorBaseLogTable, FrobeniusFactorBase,
        IndividualLogInterruption, IndividualLogObserver, IndividualLogReport, IndividualLogSolver,
        KoblitzCurve, KoblitzIcOptions, LinearAlgebra, PdpOutcome, QueryAttempt,
    },
};
use crate::binary_ecc::{BinaryPoint, F2mElement};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{collections::BTreeSet, time::Instant};

// Reuse the audited independent arithmetic source unchanged. It does not call
// the producer's curve/field operations. Both files must enter a future freeze.
#[allow(dead_code)]
#[path = "../bin/icprog/json.rs"]
mod json;
#[allow(dead_code)]
#[path = "../bin/icprog/oracle.rs"]
mod oracle;

const ORDER: u64 = 65587;
const GENERATOR: [u64; 2] = [43693, 23339];

#[path = "prepared_n17_target_sat.rs"]
mod external_sat;
pub use external_sat::{SatTargetBackend, SatTargetPlan};

/// A new frozen worker must durably publish these records before returning.
/// Callback success cannot establish native execution or evidence custody.
pub trait F5TargetObserver {
    fn started(&mut self, query: &Value) -> Result<(), String>;
    fn completed(&mut self, attempt: &Value) -> Result<(), String>;
}

struct F5RecordBridge<'a, O> {
    public: [u64; 2],
    observer: &'a mut O,
}
impl<O: F5TargetObserver> IndividualLogObserver for F5RecordBridge<'_, O> {
    fn started(&mut self, trial: u64, a: u64, b: u64) -> Result<(), String> {
        self.observer.started(&json!({"trial":trial,"a":a,"b":b,
            "target":self.public,"query_rule":"seeded-sample-aG-plus-bQ",
            "decomposition_started":false}))
    }
    fn completed(&mut self, attempt: &QueryAttempt) -> Result<(), String> {
        self.observer
            .completed(&json!({"target":self.public,"attempt":attempt,
            "record_kind":"completed-pdp-before-scalar-recovery"}))
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct F5TargetPlan {
    pub schema_version: u32,
    pub question: String,
    pub algorithm_seed: u64,
    pub max_queries: usize,
}

impl F5TargetPlan {
    pub fn validate(&self) -> Result<(), String> {
        require(
            self.schema_version == 1
                && self.question == "prepared-public-target-n17-f5-v1"
                && (1..=8).contains(&self.max_queries),
            "outside the bounded disclosed n17 target interface",
        )
    }
}

fn require(ok: bool, why: &str) -> Result<(), String> {
    if ok {
        Ok(())
    } else {
        Err(why.into())
    }
}
fn interruption_shape(stopped: &IndividualLogInterruption) -> Result<(), String> {
    let r = &stopped.report;
    let count = r
        .attempts
        .as_ref()
        .ok_or("interruption lacks prefix")?
        .len();
    let chronology = if stopped.stage == "source_gate" {
        stopped.trial.is_none() && stopped.coefficients.is_none() && r.trials == 0 && count == 0
    } else {
        stopped.trial.and_then(|t| usize::try_from(t).ok()) == r.trials.checked_sub(1)
            && stopped
                .coefficients
                .is_some_and(|[a, b]| (1..ORDER).contains(&a) && (1..ORDER).contains(&b))
            && match stopped.stage {
                "start_record" => count.checked_add(1) == Some(r.trials),
                "completion_record" | "invalid_model" => {
                    count == r.trials
                        && r.attempts
                            .as_ref()
                            .unwrap()
                            .last()
                            .is_some_and(|last| Some([last.a, last.b]) == stopped.coefficients)
                }
                _ => false,
            }
    };
    require(chronology, "invalid interrupted stage or prefix chronology")
}
fn fixture() -> Value {
    json!({"degree":17,"curve_a":1,"subgroup_order":ORDER,"group_order":131174,
        "cofactor":2,"generator":GENERATOR,"lambda":17184,
        "irreducible":{"degree":17,"low_terms":[0,3]},"targets":[],
        "target_seeds":[],"target_scalar_constructed":false})
}
fn encoded(p: &BinaryPoint) -> Result<[u64; 2], String> {
    match p {
        BinaryPoint::Affine { x, y } => Ok([
            x.to_biguint().to_u64().ok_or("coordinate overflow")?,
            y.to_biguint().to_u64().ok_or("coordinate overflow")?,
        ]),
        BinaryPoint::Infinity => Err("identity has no affine encoding".into()),
    }
}
fn point(p: [u64; 2]) -> BinaryPoint {
    BinaryPoint::Affine {
        x: F2mElement::from_biguint(&BigUint::from(p[0]), 17),
        y: F2mElement::from_biguint(&BigUint::from(p[1]), 17),
    }
}
fn independent_point(c: &oracle::Curve, p: [u64; 2]) -> Result<oracle::Point, String> {
    c.decode(&json::J::Arr(
        p.into_iter().map(|v| json::J::Int(v.into())).collect(),
    ))
}
fn options(seed: u64, max_queries: usize) -> KoblitzIcOptions {
    KoblitzIcOptions {
        m: 3,
        seed,
        max_trials: max_queries,
        strategy: DecompositionStrategy::Groebner,
        engine: SolverEngine::MatrixF5 { max_degree: 3 },
        node_budget: 8192,
        collection_window: None,
        allow_direct_relation: false,
        linear_algebra: LinearAlgebra::Dense,
        ..Default::default()
    }
}

/// Immutable reusable state derived from an ordinary preparation's math input.
/// No target or known answer is accepted by this constructor.
pub struct PreparedN17Target {
    curve: KoblitzCurve,
    base: FrobeniusFactorBase,
    table: FactorBaseLogTable,
    independent: oracle::Curve,
    geometry: Vec<oracle::Point>,
    preparation_family: String,
    preparation_math_sha256: String,
    reusable_validation_ns: u64,
}
impl PreparedN17Target {
    /// Reconstruct geometry/column order and independently replay every log.
    /// The caller must separately verify original preparation provenance.
    pub fn from_ordinary_math(input: &Value) -> Result<Self, String> {
        let started = Instant::now();
        let preparation_math_sha256 = super::prepared_sat_control::canonical_sha(input)?;
        require(
            input.as_object().is_some_and(|o| {
                o.len() == 8
                    && [
                        "schema_version",
                        "question",
                        "plan",
                        "fixture",
                        "geometric_base",
                        "attempts",
                        "stop",
                        "claimed_column_logs",
                    ]
                    .iter()
                    .all(|&k| o.contains_key(k))
            }) && input["plan"].as_object().is_some_and(|o| {
                o.len() == 4
                    && ["family", "algorithm_seed", "planned_queries", "query_law"]
                        .iter()
                        .all(|&k| o.contains_key(k))
            }) && input["plan"]["algorithm_seed"].as_u64().is_some()
                && input["plan"]["query_law"] == "independent-probe-scalar-rand08-v1",
            "unexpected target or preparation fields/query law",
        )?;
        require(
            input["schema_version"] == 1
                && input["question"] == "ordinary-preparation-mathematics-n17-v1"
                && input["fixture"] == fixture()
                && input["stop"] == "panel_complete",
            "not a complete target-free exact n17 mathematical preparation",
        )?;
        let family = input["plan"]["family"]
            .as_str()
            .ok_or("missing preparation family")?;
        require(
            matches!(family, "matrix_f5" | "cryptominisat"),
            "preparation family differs",
        )?;
        let count = input["plan"]["planned_queries"]
            .as_u64()
            .ok_or("missing panel count")?;
        require(
            (1..=512).contains(&count)
                && input["attempts"]
                    .as_array()
                    .is_some_and(|v| v.len() as u64 == count),
            "incomplete preparation panel",
        )?;
        let independent = oracle::Curve::new(&json::parse(&fixture().to_string())?)?;
        let curve = KoblitzCurve::new(1, 17).ok_or("n17 curve unavailable")?;
        require(
            curve.n == 17
                && curve.a == 1
                && curve.k == 1
                && curve.subgroup_order.to_u64() == Some(ORDER)
                && curve.group_order.to_u64() == Some(131174)
                && curve.cofactor.to_u64() == Some(2)
                && curve.lambda.to_u64() == Some(17184)
                && curve.curve.irreducible.low_terms == vec![0, 3]
                && encoded(curve.generator())? == GENERATOR,
            "installed producer curve differs from exact disclosed fixture",
        )?;
        let base = FactorBaseSpec::StandardSubspace { dimension: 6 }.materialize(&curve)?;
        let actual = base
            .points
            .iter()
            .map(encoded)
            .collect::<Result<Vec<_>, _>>()?;
        require(
            actual.len() == 63 && input["geometric_base"] == json!(actual),
            "geometric order differs",
        )?;
        let geometry = actual
            .iter()
            .map(|&p| independent_point(&independent, p))
            .collect::<Result<Vec<_>, _>>()?;
        require(
            geometry.iter().copied().collect::<BTreeSet<_>>().len() == 63
                && geometry
                    .iter()
                    .filter_map(|&p| independent.mul(p, 2))
                    .collect::<BTreeSet<_>>()
                    .len()
                    == 62,
            "actual usable base count differs",
        )?;
        let opts = options(0, 1);
        let matrix = FactorBaseLogSolver::new(&curve, &base, &opts).ok_or("columns unavailable")?;
        let columns = matrix.matrix_snapshot().column_points;
        let logs = input["claimed_column_logs"]
            .as_array()
            .ok_or("complete logs unavailable")?;
        require(
            columns.len() == 29 && logs.len() == columns.len(),
            "folded columns differ",
        )?;
        let table = FactorBaseLogTable {
            columns: columns
                .iter()
                .zip(logs)
                .map(|((x, y), l)| {
                    let p = [
                        x.parse::<u64>().map_err(|_| "column x overflow")?,
                        y.parse::<u64>().map_err(|_| "column y overflow")?,
                    ];
                    let log = l.as_u64().ok_or("noninteger column log")?;
                    let decoded = independent_point(&independent, p)?;
                    require(
                        log < ORDER
                            && decoded.is_some()
                            && independent.mul(Some(independent.g), log.into()) == decoded,
                        "independent column scalar replay failed",
                    )?;
                    Ok((point(p), BigUint::from(log)))
                })
                .collect::<Result<Vec<_>, String>>()?,
        };
        require(table.verify(&curve), "producer column replay failed")?;
        require(
            IndividualLogSolver::new(&curve, &base, &table, &opts, None).is_some(),
            "table does not cover the reconstructed signed/Frobenius columns",
        )?;
        Ok(Self {
            curve,
            base,
            table,
            independent,
            geometry,
            preparation_family: family.into(),
            preparation_math_sha256,
            reusable_validation_ns: started
                .elapsed()
                .as_nanos()
                .try_into()
                .map_err(|_| "duration overflow")?,
        })
    }

    /// Public input is a point only. This source API is not a dispatch authority.
    /// A new frozen, one-use controller is required before a scientific solve.
    pub fn solve_f5(&self, public: [u64; 2], plan: &F5TargetPlan) -> Result<Value, String> {
        plan.validate()?;
        let setup = Instant::now();
        let opts = options(plan.algorithm_seed, plan.max_queries);
        let solver = IndividualLogSolver::new(&self.curve, &self.base, &self.table, &opts, None)
            .ok_or("prepared F5 dispatch unavailable")?;
        let mut dispatch =
            serde_json::to_value(solver.admission_dispatch()).map_err(|e| e.to_string())?;
        dispatch["engine"] = json!(opts.engine.effective());
        dispatch["node_budget"] = json!(opts.node_budget);
        dispatch["collapse_negation"] = json!(opts.collapse_negation);
        dispatch["collapse_projected_orbits"] = json!(opts.collapse_projected_orbits);
        let solver_setup_ns: u64 = setup
            .elapsed()
            .as_nanos()
            .try_into()
            .map_err(|_| "duration overflow")?;
        self.run_observed(public, plan, dispatch, solver_setup_ns, |q| {
            solver.solve_observed(q)
        })
    }

    /// Same bounded F5 query law, with a start/completion callback per attempt.
    /// This entry still requires a separately frozen one-use native executor.
    pub fn solve_f5_recorded(
        &self,
        public: [u64; 2],
        plan: &F5TargetPlan,
        observer: &mut impl F5TargetObserver,
    ) -> Result<Value, String> {
        plan.validate()?;
        let setup = Instant::now();
        let opts = options(plan.algorithm_seed, plan.max_queries);
        let solver = IndividualLogSolver::new(&self.curve, &self.base, &self.table, &opts, None)
            .ok_or("prepared F5 recording dispatch unavailable")?;
        let mut dispatch =
            serde_json::to_value(solver.admission_dispatch()).map_err(|e| e.to_string())?;
        dispatch["engine"] = json!(opts.engine.effective());
        dispatch["node_budget"] = json!(opts.node_budget);
        dispatch["collapse_negation"] = json!(opts.collapse_negation);
        dispatch["collapse_projected_orbits"] = json!(opts.collapse_projected_orbits);
        dispatch["attempt_record_contract"] = json!("start-and-completed-pdp-before-recovery-v1");
        let solver_setup_ns = setup
            .elapsed()
            .as_nanos()
            .try_into()
            .map_err(|_| "duration overflow")?;
        let mut bridge = F5RecordBridge { public, observer };
        self.run_observed_result(public, plan, dispatch, solver_setup_ns, |q| {
            solver.solve_n17_durable(q, &mut bridge)
        })
    }

    fn run_observed(
        &self,
        public: [u64; 2],
        plan: &F5TargetPlan,
        dispatch: Value,
        solver_setup_ns: u64,
        solve: impl FnOnce(&BinaryPoint) -> IndividualLogReport,
    ) -> Result<Value, String> {
        self.run_observed_result(public, plan, dispatch, solver_setup_ns, |q| Ok(solve(q)))
    }

    fn run_observed_result(
        &self,
        public: [u64; 2],
        plan: &F5TargetPlan,
        dispatch: Value,
        solver_setup_ns: u64,
        solve: impl FnOnce(&BinaryPoint) -> Result<IndividualLogReport, IndividualLogInterruption>,
    ) -> Result<Value, String> {
        let session = Session::begin_strict().map_err(str::to_string)?;
        clock::begin_online(Phase::TargetQuery);
        // Range/on-curve/subgroup work is target dependent and stays charged.
        let decoded = independent_point(&self.independent, public).and_then(|q| {
            require(
                q.is_some() && self.independent.mul(q, ORDER.into()).is_none(),
                "public target is outside the nonidentity prime subgroup",
            )?;
            Ok(q)
        });
        let mut report = None;
        let mut failure = None;
        let mut scalar = None;
        let mut certificate = None;
        let mut interruption = None;
        match decoded {
            Err(e) => failure = Some(e),
            Ok(q) => {
                let (observed, interrupted) = match solve(&point(public)) {
                    Ok(report) => (report, false),
                    Err(stopped) => {
                        let shape = interruption_shape(&stopped);
                        let prefix_shape_valid = shape.is_ok();
                        failure =
                            Some(shape.err().unwrap_or_else(|| {
                                format!("{}: {}", stopped.stage, stopped.reason)
                            }));
                        interruption = Some(json!({"stage":stopped.stage,"trial":stopped.trial,
                            "coefficients":stopped.coefficients,
                            "reason":stopped.reason,"prefix_shape_valid":prefix_shape_valid,
                            "runtime_custody_admitted":false}));
                        (*stopped.report, true)
                    }
                };
                clock::mark(Phase::TargetRelationCheck);
                match self.check_attempt_prefix(q, &observed, plan.max_queries, interrupted) {
                    Err(e) => failure = Some(e),
                    Ok(()) if !interrupted => {
                        clock::mark(Phase::RecoveryCheck);
                        if let Some(k) = observed.log.as_ref() {
                            match k.to_u64().filter(|&k| k < ORDER) {
                                Some(k)
                                    if self.independent.mul(Some(self.independent.g), k.into())
                                        == q =>
                                {
                                    scalar = Some(k);
                                    certificate = Some(
                                        json!({"check":"independent-shift-and-add-scalar-replay",
                                        "generator":GENERATOR,"scalar":k,"target":public,"inside_online_interval":true}),
                                    );
                                }
                                _ => {
                                    failure = Some("independent target scalar replay failed".into())
                                }
                            }
                        } else {
                            failure = Some("query budget exhausted without scalar recovery".into());
                        }
                    }
                    Ok(()) => {} // A partial transcript is never a successful solve.
                }
                report = Some(observed);
            }
        }
        clock::end_online();
        let costs = session.finish().map_err(str::to_string)?;
        require(
            costs
                .online_wall_ns
                .is_some_and(|n| costs.online_phases_ns.values().flatten().sum::<u64>() == n),
            "online phases do not close",
        )?;
        let report = report.as_ref().map(|r| {
            let relation = r
                .relation
                .as_ref()
                .map(|rel| json!({"a":rel.a,"b":rel.b,"points":rel.points}));
            json!({"trials":r.trials,"claimed_scalar":r.log.as_ref().map(ToString::to_string),
                "relation":relation,"attempts":r.attempts})
        });
        // Keep the raw clock namespace and export the five measurement-contract
        // names explicitly. Unentered failure phases remain unknown.
        let exclusive_ic_online_phases_ns = json!({
            "target_query":costs.online_phases_ns["target_query"],
            "target_pdp":costs.online_phases_ns["target_pdp"],
            "target_relation_check":costs.online_phases_ns["target_relation_check"],
            "target_descent":costs.online_phases_ns["target_descent"],
            "target_recovery_check":costs.online_phases_ns["recovery_check"]});
        Ok(
            json!({"schema_version":1,"question":"prepared-public-target-n17-f5-producer-v1",
            "plan":plan,"target":public,"target_count":1,"preparation_family":self.preparation_family,
            "preparation_math_sha256":self.preparation_math_sha256,
            "actual_usable_base_points":62,"geometric_base_points":63,"effective_columns":29,
            "dispatch":dispatch,"reusable_validation_ns":self.reusable_validation_ns,
            "solver_setup_ns":solver_setup_ns,"report":report,"verified_scalar":scalar,
            "independent_recovery_certificate":certificate,"failure":failure,"interruption":interruption,"costs":costs,
            "exclusive_ic_online_phases_ns":exclusive_ic_online_phases_ns,
            "timing_class":"ordinary-host-diagnostic","candidate_id":null,"workload_id":null,"run_id":null,
            "source_bound_execution_admitted":false,"fresh_paired_qualification":false,
            "headline_eligible":false,"promotion_eligible":false,"full_goal_complete":false,
            "online_speedup":null,"operation_counts":null}),
        )
    }

    fn check_attempt_prefix(
        &self,
        q: oracle::Point,
        report: &IndividualLogReport,
        cap: usize,
        interrupted: bool,
    ) -> Result<(), String> {
        let attempts = report
            .attempts
            .as_ref()
            .ok_or("missing failed-query transcript")?;
        require(
            report.trials <= cap
                && (attempts.len() == report.trials
                    || (interrupted && attempts.len().checked_add(1) == Some(report.trials))),
            "attempt chronology/cap differs",
        )?;
        for (i, attempt) in attempts.iter().enumerate() {
            require(
                attempt.trial == i as u64
                    && (1..ORDER).contains(&attempt.a)
                    && (1..ORDER).contains(&attempt.b),
                "invalid probe chronology or coefficients",
            )?;
            if let Some(indices) = &attempt.pdp.points {
                require(
                    attempt.pdp.outcome == PdpOutcome::Witness && indices.len() == 3,
                    "invalid witness disposition or summand count",
                )?;
                let mut sum = None;
                for &index in indices {
                    let p = *self
                        .geometry
                        .get(index)
                        .ok_or("witness index outside geometric base")?;
                    sum = self.independent.add(sum, p);
                }
                let probe = self.independent.add(
                    self.independent
                        .mul(Some(self.independent.g), attempt.a.into()),
                    self.independent.mul(q, attempt.b.into()),
                );
                require(
                    sum == probe,
                    "independent full-point relation replay failed",
                )?;
            } else {
                require(
                    attempt.pdp.outcome != PdpOutcome::Witness,
                    "witness lacks points",
                )?;
            }
        }
        if interrupted {
            return require(
                report.log.is_none() && report.relation.is_none(),
                "interrupted solve cannot retain a claimed scalar or recovery relation",
            );
        }
        match (&report.log, &report.relation) {
            (Some(_), Some(rel)) => require(
                attempts.last().is_some_and(|last| {
                    last.a == rel.a
                        && last.b == rel.b
                        && last.pdp.points.as_ref() == Some(&rel.points)
                }),
                "recovered scalar lacks the final observed relation",
            ),
            (None, None) => require(report.trials == cap, "unexplained early solver termination"),
            _ => Err("scalar/relation completeness differs".into()),
        }
    }
}

#[cfg(test)]
#[path = "prepared_n17_target_f5_tests.rs"]
mod recording_tests;

#[cfg(test)]
mod boundary_tests {
    use super::*;
    use crate::cryptanalysis::koblitz_index_calculus::{
        DescentRelation, PdpAttempt, PdpSolverStats, QueryAttempt,
    };
    use std::sync::Mutex;

    pub(super) static CLOCK: Mutex<()> = Mutex::new(());

    // Retained mathematics only. No old controller, capsule, worker or solver
    // is executed. These deliberately disclosed controls are never fresh data.
    pub(super) fn retained_math() -> Value {
        let old: Value = serde_json::from_str(include_str!(
            "../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json"
        )).unwrap();
        let inputs = &old["certificate"]["inputs"];
        // Historical certificate coordinates are decimal strings. Its record
        // has the same coordinates as integers. Convert this fixture explicitly;
        // never relax the constructor or rewrite the original artifact.
        let curve = KoblitzCurve::new(1, 17).unwrap();
        let base = FactorBaseSpec::StandardSubspace { dimension: 6 }
            .materialize(&curve)
            .unwrap();
        let geometry = base
            .points
            .iter()
            .map(|p| json!(encoded(p).unwrap()))
            .collect::<Vec<_>>();
        assert_eq!(
            json!(geometry),
            old["record"]["factor_base"]["points"],
            "ordered geometric points must match after encoding conversion"
        );
        let opts = options(0, 1);
        let matrix = FactorBaseLogSolver::new(&curve, &base, &opts).unwrap();
        let logs = matrix
            .matrix_snapshot()
            .column_points
            .iter()
            .map(|(x, y)| {
                let p = json!([x.parse::<u64>().unwrap(), y.parse::<u64>().unwrap()]);
                old["record"]["column_logs"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .find(|v| v["point"] == p)
                    .expect("control column point missing")["log"]
                    .clone()
            })
            .collect::<Vec<_>>();
        let mut attempts = inputs["attempts"].clone();
        for a in attempts.as_array_mut().unwrap() {
            if let Some(indices) = a["indices"].as_array_mut() {
                for i in indices {
                    let p = &old["record"]["factor_base"]["points"][i.as_u64().unwrap() as usize];
                    *i = json!(geometry.iter().position(|v| v == p).unwrap());
                }
            }
        }
        json!({"schema_version":1,"question":"ordinary-preparation-mathematics-n17-v1",
            "fixture":fixture(),"plan":{"family":"matrix_f5","planned_queries":216,
                "algorithm_seed":2026093032u64,"query_law":"independent-probe-scalar-rand08-v1"},
            "geometric_base":geometry,"attempts":attempts,
            "stop":"panel_complete","claimed_column_logs":logs})
    }
    pub(super) fn plan(cap: usize) -> F5TargetPlan {
        F5TargetPlan {
            schema_version: 1,
            question: "prepared-public-target-n17-f5-v1".into(),
            algorithm_seed: 2026100310,
            max_queries: cap,
        }
    }
    fn attempt(t: u64, a: u64, b: u64, points: Option<Vec<usize>>) -> QueryAttempt {
        QueryAttempt {
            trial: t,
            a,
            b,
            pdp: PdpAttempt {
                outcome: if points.is_some() {
                    PdpOutcome::Witness
                } else {
                    PdpOutcome::ProvedUnsat
                },
                points,
                stats: PdpSolverStats::None,
            },
        }
    }
    pub(super) fn retained_report(state: &PreparedN17Target) -> IndividualLogReport {
        let old: Value = serde_json::from_str(include_str!(
            "../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json"
        )).unwrap();
        let indices = [26, 36, 46]
            .map(|i| {
                let p: [u64; 2] =
                    serde_json::from_value(old["record"]["factor_base"]["points"][i].clone())
                        .unwrap();
                let decoded = independent_point(&state.independent, p).unwrap();
                state.geometry.iter().position(|&v| v == decoded).unwrap()
            })
            .to_vec();
        IndividualLogReport {
            trials: 3,
            log: Some(BigUint::from(24886u64)),
            relation: Some(DescentRelation {
                a: 45669,
                b: 63073,
                points: indices.clone(),
            }),
            attempts: Some(vec![
                attempt(0, 57707, 37266, None),
                attempt(1, 43768, 33418, None),
                attempt(2, 45669, 63073, Some(indices)),
            ]),
        }
    }
    fn replay(state: &PreparedN17Target, report: IndividualLogReport) -> Value {
        state
            .run_observed(
                [52411, 72106],
                &plan(8),
                json!({"kind":"retained-transcript-control-no-solver-execution"}),
                0,
                |_| {
                    clock::mark(Phase::TargetPdp);
                    clock::mark(Phase::TargetDescent);
                    report
                },
            )
            .unwrap()
    }

    #[test]
    fn reconstructs_exact_geometry_and_rejects_mutated_logs() {
        let _guard = CLOCK.lock().unwrap();
        let original = retained_math();
        let state = PreparedN17Target::from_ordinary_math(&original).unwrap();
        assert_eq!(state.geometry.len(), 63);
        assert_eq!(state.table.len(), 29);
        let mut bad = original.clone();
        bad["claimed_column_logs"][0] = json!(1);
        assert!(PreparedN17Target::from_ordinary_math(&bad).is_err());
        bad = original.clone();
        bad["geometric_base"].as_array_mut().unwrap().swap(1, 2);
        assert!(PreparedN17Target::from_ordinary_math(&bad).is_err());
        bad = original.clone();
        bad["fixture"]["degree"] = json!(19);
        assert!(PreparedN17Target::from_ordinary_math(&bad).is_err());
        bad = original.clone();
        bad["target"] = json!([52411, 72106]);
        assert!(PreparedN17Target::from_ordinary_math(&bad).is_err());
        bad = original.clone();
        bad["plan"]["known_scalar"] = json!(24886);
        assert!(PreparedN17Target::from_ordinary_math(&bad).is_err());
        bad = original.clone();
        bad["claimed_column_logs"][0] = json!(63339.0);
        assert!(PreparedN17Target::from_ordinary_math(&bad).is_err());
        bad = original;
        bad["attempts"].as_array_mut().unwrap().pop();
        assert!(PreparedN17Target::from_ordinary_math(&bad).is_err());
    }

    #[test]
    fn retained_success_charges_independent_replay_and_preserves_failures() {
        let _guard = CLOCK.lock().unwrap();
        let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
        let result = replay(&state, retained_report(&state));
        assert_eq!(result["verified_scalar"], 24886);
        assert!(result["failure"].is_null());
        assert_eq!(result["report"]["attempts"].as_array().unwrap().len(), 3);
        assert_eq!(
            result["independent_recovery_certificate"]["inside_online_interval"],
            true
        );
        let phases = result["costs"]["online_phases_ns"].as_object().unwrap();
        assert_eq!(
            phases.values().filter_map(Value::as_u64).sum::<u64>(),
            result["costs"]["online_wall_ns"].as_u64().unwrap()
        );
        assert!(phases["recovery_check"].as_u64().unwrap() > 0);
        let exclusive = result["exclusive_ic_online_phases_ns"].as_object().unwrap();
        assert_eq!(exclusive.len(), 5);
        assert_eq!(
            exclusive.values().map(|v| v.as_u64().unwrap()).sum::<u64>(),
            result["costs"]["online_wall_ns"].as_u64().unwrap()
        );
        assert_eq!(result["source_bound_execution_admitted"], false);
        assert_eq!(result["fresh_paired_qualification"], false);
        assert!(result["online_speedup"].is_null());
        assert!(result["candidate_id"].is_null());
    }

    #[test]
    fn rejects_corrupt_witness_scalar_and_missing_failed_attempts() {
        let _guard = CLOCK.lock().unwrap();
        let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
        let mut bad = retained_report(&state);
        bad.attempts.as_mut().unwrap()[2]
            .pdp
            .points
            .as_mut()
            .unwrap()[0] = 1;
        let result = replay(&state, bad);
        assert!(result["verified_scalar"].is_null());
        assert!(result["failure"]
            .as_str()
            .unwrap()
            .contains("relation replay failed"));
        bad = retained_report(&state);
        bad.log = Some(BigUint::from(1u64));
        let result = replay(&state, bad);
        assert!(result["verified_scalar"].is_null());
        assert!(result["failure"]
            .as_str()
            .unwrap()
            .contains("scalar replay failed"));
        bad = retained_report(&state);
        bad.attempts.as_mut().unwrap().remove(0);
        assert!(replay(&state, bad)["failure"]
            .as_str()
            .unwrap()
            .contains("chronology"));
    }

    #[test]
    fn retains_budget_exhaustion_and_rejects_unexplained_early_stop() {
        let _guard = CLOCK.lock().unwrap();
        let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
        let report = IndividualLogReport {
            trials: 2,
            log: None,
            relation: None,
            attempts: Some(vec![
                attempt(0, 57707, 37266, None),
                attempt(1, 43768, 33418, None),
            ]),
        };
        let result = state
            .run_observed(
                [52411, 72106],
                &plan(2),
                json!({"kind":"no-search-control"}),
                0,
                |_| report.clone(),
            )
            .unwrap();
        assert_eq!(result["report"]["trials"], 2);
        assert!(result["failure"]
            .as_str()
            .unwrap()
            .contains("budget exhausted"));
        assert!(result["verified_scalar"].is_null());
        assert!(replay(&state, report)["failure"]
            .as_str()
            .unwrap()
            .contains("early solver termination"));
    }

    #[test]
    fn checks_public_point_inside_interval_before_any_solver_call() {
        let _guard = CLOCK.lock().unwrap();
        let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
        for public in [[1 << 17, 0], [0, 0], [0, 1]] {
            let result = state
                .run_observed(public, &plan(1), json!({}), 0, |_| {
                    panic!("invalid target reached solver")
                })
                .unwrap();
            assert!(result["failure"].as_str().is_some());
            assert!(result["report"].is_null());
            assert!(result["costs"]["online_wall_ns"].as_u64().is_some());
            assert!(result["exclusive_ic_online_phases_ns"]["target_pdp"].is_null());
        }
    }

    #[test]
    fn plan_cannot_select_another_curve_backend_or_unbounded_search() {
        assert!(plan(1).validate().is_ok());
        assert!(plan(0).validate().is_err());
        assert!(plan(9).validate().is_err());
        let mut value = serde_json::to_value(plan(1)).unwrap();
        value["known_scalar"] = json!(24886);
        assert!(serde_json::from_value::<F5TargetPlan>(value).is_err());
    }
}
