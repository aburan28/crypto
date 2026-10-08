//! External SAT source/model boundary on the fixed educational n17 fixture.
//! This callback contract does not attest actual native execution or custody.
use super::*;
use crate::cryptanalysis::prepared_sat_control::{self as sat, QueryOutput};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SatTargetPlan {
    pub schema_version: u32,
    pub question: String,
    pub algorithm_seed: u64,
    pub max_queries: usize,
    pub export_nonce: u64,
    pub conflict_budget: u64,
    pub exporter_timeout_ms: u64,
    pub solver_timeout_ms: u64,
    pub controller_timeout_ms: u64,
}
impl SatTargetPlan {
    pub fn validate(&self) -> Result<(), String> {
        require(
            self.schema_version == 1
                && self.question == "prepared-public-target-n17-cms-v1"
                && (1..=8).contains(&self.max_queries)
                && (100_000..=1_000_000).contains(&self.conflict_budget)
                && (1..=30_000).contains(&self.exporter_timeout_ms)
                // Match the registered one-million-conflict natural panel's
                // 120-second per-query watchdog. The separate 900-second
                // one-target controller cap still bounds total work.
                && (1..=120_000).contains(&self.solver_timeout_ms)
                && (30_000..=900_000).contains(&self.controller_timeout_ms),
            "outside the bounded disclosed n17 CMS target interface",
        )
    }
}

/// A separately sealed native worker must implement these durable boundaries,
/// retain original ANF/CNF/stdout and drain native roles before returning.
/// Callback output and booleans cannot establish any execution admission.
pub trait SatTargetBackend {
    fn started(&mut self, query: &Value) -> Result<(), String>;
    fn query(
        &mut self,
        trial: usize,
        public: [u64; 2],
        plan: &SatTargetPlan,
    ) -> Result<QueryOutput, String>;
    fn completed(&mut self, attempt: &Value) -> Result<(), String>;
}

fn mul_mod(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % ORDER as u128) as u64
}
fn inverse(mut a: u64) -> u64 {
    let mut exponent = ORDER - 2;
    let mut value = 1;
    while exponent != 0 {
        if exponent & 1 == 1 {
            value = mul_mod(value, a);
        }
        a = mul_mod(a, a);
        exponent >>= 1;
    }
    value
}

impl PreparedN17Target {
    /// Derive cofactor-image logs from the replayed columns and public
    /// Frobenius/sign action; never load a historical projection table.
    fn projected_logs(&self) -> Result<Vec<Option<u64>>, String> {
        let mut by_point = BTreeMap::new();
        for (p, l) in &self.table.columns {
            let mut current = independent_point(&self.independent, encoded(p)?)?;
            let mut log = l.to_u64().ok_or("column log overflow")?;
            for _ in 0..17 {
                for (image, value) in [
                    (current, log),
                    (self.independent.neg(current), (ORDER - log) % ORDER),
                ] {
                    let key = image.ok_or("column orbit contains identity")?;
                    if let Some(old) = by_point.insert(key, value) {
                        require(old == value, "inconsistent signed/Frobenius logs")?;
                    }
                }
                current = self.independent.frob(current);
                log = mul_mod(log, 17184);
            }
        }
        self.geometry
            .iter()
            .map(|&p| {
                let image = self.independent.mul(p, 2);
                if let Some(point) = image {
                    let log = *by_point
                        .get(&point)
                        .ok_or("usable geometric image missing from columns")?;
                    require(
                        self.independent.mul(Some(self.independent.g), log.into()) == image,
                        "projected log fails independent scalar replay",
                    )?;
                    Ok(Some(log))
                } else {
                    Ok(None)
                }
            })
            .collect()
    }

    fn lift_sat_model(&self, model: &[bool], query: oracle::Point) -> Option<[usize; 3]> {
        if model.len() < 51 {
            return None;
        }
        let xs: [u64; 3] =
            std::array::from_fn(|i| (0..6).fold(0, |v, b| v | ((model[6 * i + b] as u64) << b)));
        let choices = xs.map(|x| {
            self.geometry
                .iter()
                .enumerate()
                .filter(|(_, p)| p.is_some_and(|p| p.0 == x))
                .map(|(i, _)| i)
                .collect::<Vec<_>>()
        });
        for &a in &choices[0] {
            for &b in &choices[1] {
                for &c in &choices[2] {
                    if self.independent.add(
                        self.independent.add(self.geometry[a], self.geometry[b]),
                        self.geometry[c],
                    ) == query
                    {
                        return Some([a, b, c]);
                    }
                }
            }
        }
        None
    }

    /// No native executable is selected here. Actual executor/asset/source
    /// identity and watchdog/drain receipts are mandatory outer admission gates.
    pub fn solve_external_sat(
        &self,
        public: [u64; 2],
        plan: &SatTargetPlan,
        backend: &mut impl SatTargetBackend,
    ) -> Result<Value, String> {
        plan.validate()?;
        let setup = Instant::now();
        let projected_logs = self.projected_logs()?;
        let fast = super::super::koblitz_fast::FastCurve::new(&self.curve.curve)
            .ok_or("portable native curve unavailable")?;
        let generator = fast.lift(self.curve.generator());
        let mut rng = StdRng::seed_from_u64(plan.algorithm_seed ^ 0x44455343454e5400);
        let solver_setup_ns: u64 = setup
            .elapsed()
            .as_nanos()
            .try_into()
            .map_err(|_| "duration overflow")?;
        let session = Session::begin_strict().map_err(str::to_string)?;
        clock::begin_online(Phase::TargetQuery);
        let decoded = independent_point(&self.independent, public).and_then(|q| {
            require(
                q.is_some() && self.independent.mul(q, ORDER.into()).is_none(),
                "public target is outside the nonidentity prime subgroup",
            )?;
            Ok(q)
        });
        let mut attempts = Vec::new();
        let mut scalar = None;
        let mut certificate = None;
        let mut failure = None;
        if let Ok(q) = decoded {
            let target = fast.lift(&point(public));
            for trial in 0..plan.max_queries {
                clock::mark(Phase::TargetQuery);
                let (a, b) = (rng.gen_range(1..ORDER), rng.gen_range(1..ORDER));
                let query = fast.add(fast.mul_u64(generator, a), fast.mul_u64(target, b));
                let public_query = sat::encoded(query);
                let mut row = json!({"trial":trial,"a":a,"b":b,"public_query":public_query,
                    "outcome":if public_query.is_some(){"query_planned"}else{"identity_query"},
                    "source_model_valid":false,"witness_indices":null,
                    "candidate_scalar":null,"backend_called":false});
                let mut fatal = false;
                if let Err(e) = backend.started(&row) {
                    row["outcome"] = json!("start_record_failure");
                    row["reason"] = json!(e);
                    fatal = true;
                } else if let Some(probe) = public_query {
                    clock::mark(Phase::TargetPdp);
                    row["backend_called"] = json!(true);
                    match backend.query(trial, probe, plan) {
                        Err(e) => {
                            row["outcome"] = json!("query_transport_failure");
                            row["reason"] = json!(e);
                            fatal = true;
                        }
                        Ok(output) => {
                            clock::mark(Phase::TargetRelationCheck);
                            let status = sat::native_status(&output.native);
                            row["native_status"] = json!(status);
                            row["source_evidence"] = json!({"manifest":output.manifest,
                                "source_receipt":output.source_receipt,
                                "anf_bytes":output.anf.len(),"anf_sha256":sat::sha256(output.anf.as_bytes()),
                                "cnf_bytes":output.cnf.len(),"cnf_sha256":sat::sha256(output.cnf.as_bytes()),
                                "stdout_bytes":output.native.stdout.len(),"stdout_sha256":sat::sha256(output.native.stdout.as_bytes()),
                                "native":{"exit_code":output.native.exit_code,"timed_out":output.native.timed_out,
                                    "receipt":output.native.receipt},"payload_policy":"outer executor retains original files"});
                            let checked = require(
                                output.anf.len() <= 16 * 1024 * 1024
                                    && output.cnf.len() <= 16 * 1024 * 1024
                                    && output.native.stdout.len() <= 4 * 1024 * 1024,
                                "source/model output exceeded bound",
                            )
                            .and_then(|_| {
                                if matches!(status, "SAT_MODEL" | "SOURCE_UNSAT") {
                                    sat::validate_manifest(&output.manifest, probe)
                                } else {
                                    Ok(())
                                }
                            });
                            if let Err(e) = checked {
                                row["outcome"] = json!("invalid_source_output");
                                row["reason"] = json!(e);
                                fatal = true;
                            } else {
                                match status {
                                    "SAT_MODEL" => {
                                        let model = (|| {
                                            let count = output.manifest["exports"]
                                                ["cryptominisat_xor_dimacs"]["variables"]
                                                .as_u64()
                                                .filter(|n| (51..=4096).contains(n))
                                                .ok_or("invalid expanded variable count")?
                                                as usize;
                                            let model =
                                                sat::parse_model(&output.native.stdout, count)?;
                                            sat::verify_anf(&output.anf, &model[..51])?;
                                            sat::verify_cnf(&output.cnf, &model)?;
                                            let probe_point =
                                                independent_point(&self.independent, probe)?;
                                            Ok::<_, String>((model, probe_point))
                                        })();
                                        match model {
                                            Err(e) => {
                                                row["outcome"] = json!("invalid_source_model");
                                                row["reason"] = json!(e);
                                                fatal = true;
                                            }
                                            Ok((model, probe_point)) => {
                                                row["source_model_valid"] = json!(true);
                                                row["expanded_model_sha256"] = json!(sat::sha256(
                                                    &model
                                                        .iter()
                                                        .map(|&b| u8::from(b))
                                                        .collect::<Vec<_>>()
                                                ));
                                                if let Some(indices) =
                                                    self.lift_sat_model(&model, probe_point)
                                                {
                                                    row["witness_indices"] = json!(indices);
                                                    clock::mark(Phase::TargetDescent);
                                                    let sum = indices.iter().fold(0, |s, &i| {
                                                        (s + projected_logs[i].unwrap_or(0)) % ORDER
                                                    });
                                                    let k = mul_mod(
                                                        (sum + ORDER - mul_mod(2, a)) % ORDER,
                                                        inverse(mul_mod(2, b)),
                                                    );
                                                    row["candidate_scalar"] = json!(k);
                                                    clock::mark(Phase::RecoveryCheck);
                                                    if self
                                                        .independent
                                                        .mul(Some(self.independent.g), k.into())
                                                        == q
                                                    {
                                                        scalar = Some(k);
                                                        row["outcome"] =
                                                            json!("verified_point_witness");
                                                        certificate = Some(
                                                            json!({"check":"independent-shift-and-add-scalar-replay",
                                                            "generator":GENERATOR,"scalar":k,"target":public,"inside_online_interval":true}),
                                                        );
                                                    } else {
                                                        row["outcome"] =
                                                            json!("independent_recovery_failure");
                                                        fatal = true;
                                                    }
                                                } else {
                                                    row["outcome"] = json!("nonlifting_model");
                                                }
                                            }
                                        }
                                    }
                                    // A native source UNSAT is a frontend claim. The frozen
                                    // independent audit must prove the geometric negative.
                                    "SOURCE_UNSAT" => row["outcome"] = json!("source_unsat_claim"),
                                    "TIMEOUT" => row["outcome"] = json!("timeout"),
                                    "CONFLICT_BUDGET_INCONCLUSIVE" | "UNKNOWN_INCONCLUSIVE" => {
                                        row["outcome"] = json!("incomplete")
                                    }
                                    _ => {
                                        row["outcome"] = json!("native_error");
                                        fatal = true;
                                    }
                                }
                            }
                        }
                    }
                }
                if let Err(e) = backend.completed(&row) {
                    row["completion_record_failure"] = json!(e);
                    scalar = None;
                    fatal = true;
                }
                if fatal {
                    failure =
                        Some("fatal source, transport, recovery or durability failure".into());
                }
                attempts.push(row);
                if scalar.is_some() || fatal {
                    break;
                }
            }
            if scalar.is_none() && failure.is_none() {
                failure =
                    Some("query budget exhausted without independently verified recovery".into());
            }
        } else {
            failure = decoded.err();
        }
        clock::end_online();
        let costs = session.finish().map_err(str::to_string)?;
        require(
            costs
                .online_wall_ns
                .is_some_and(|n| costs.online_phases_ns.values().flatten().sum::<u64>() == n),
            "online phases do not close",
        )?;
        let exclusive = json!({"target_query":costs.online_phases_ns["target_query"],
            "target_pdp":costs.online_phases_ns["target_pdp"],
            "target_relation_check":costs.online_phases_ns["target_relation_check"],
            "target_descent":costs.online_phases_ns["target_descent"],
            "target_recovery_check":costs.online_phases_ns["recovery_check"]});
        Ok(
            json!({"schema_version":1,"question":"prepared-public-target-n17-cms-producer-v1",
            "frontend_contract":"external-cryptominisat-source-output","plan":plan,"target":public,"target_count":1,
            "preparation_family":self.preparation_family,"preparation_math_sha256":self.preparation_math_sha256,
            "actual_usable_base_points":62,"geometric_base_points":63,"effective_columns":29,
            "reusable_validation_ns":self.reusable_validation_ns,"solver_setup_ns":solver_setup_ns,
            "attempts":attempts,"verified_scalar":scalar,"independent_recovery_certificate":certificate,
            "failure":failure,"costs":costs,"exclusive_ic_online_phases_ns":exclusive,
            "timing_class":"ordinary-host-diagnostic","candidate_id":null,"workload_id":null,"run_id":null,
            "source_bound_execution_admitted":false,"fresh_paired_qualification":false,"headline_eligible":false,
            "promotion_eligible":false,"full_goal_complete":false,"online_speedup":null,"operation_counts":null}),
        )
    }
}

#[cfg(test)]
mod sat_boundary_tests {
    use super::super::boundary_tests::{retained_math, CLOCK};
    use super::*;
    use crate::cryptanalysis::prepared_sat_control::NativeOutput;
    const MANIFEST: &str = include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/manifest.json");
    const ANF: &str = include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/instance.anf");
    const CNF: &str = include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/instance.xor.cnf");
    const MODEL: &str = include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/result.json");
    const PUBLIC: [u64; 2] = [52411, 72106];

    fn plan() -> SatTargetPlan {
        SatTargetPlan {
            schema_version: 1,
            question: "prepared-public-target-n17-cms-v1".into(),
            algorithm_seed: 2026100102,
            max_queries: 2,
            export_nonce: 2026100103,
            conflict_budget: 100_000,
            exporter_timeout_ms: 30_000,
            solver_timeout_ms: 60_000,
            controller_timeout_ms: 180_000,
        }
    }
    fn model() -> Vec<bool> {
        serde_json::from_str::<Value>(MODEL).unwrap()["cnf_assignment"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_bool().unwrap())
            .collect()
    }
    fn source(corrupt: bool) -> QueryOutput {
        let mut bits = model();
        if corrupt {
            bits[0] = !bits[0];
        }
        let mut stdout = "s SATISFIABLE\nv ".to_string();
        for (i, &bit) in bits.iter().enumerate() {
            stdout.push_str(&format!(
                "{} ",
                if bit { i as i32 + 1 } else { -(i as i32 + 1) }
            ));
        }
        stdout.push_str("0\n");
        QueryOutput {
            native: NativeOutput {
                exit_code: Some(10),
                timed_out: false,
                stdout,
                receipt: json!({"kind":"retained-data-control-no-native-execution"}),
            },
            manifest: serde_json::from_str(MANIFEST).unwrap(),
            anf: ANF.into(),
            cnf: CNF.into(),
            source_receipt: json!({"kind":"retained-source-control"}),
        }
    }
    #[derive(Clone, Copy, PartialEq)]
    enum Mode {
        Good,
        Corrupt,
        WrongManifest,
        Timeout,
        Transport,
        FailStart,
        FailComplete,
    }
    struct Script {
        mode: Mode,
        starts: Vec<Value>,
        completions: Vec<Value>,
        calls: Vec<[u64; 2]>,
    }
    impl Script {
        fn new(mode: Mode) -> Self {
            Self {
                mode,
                starts: Vec::new(),
                completions: Vec::new(),
                calls: Vec::new(),
            }
        }
    }
    impl SatTargetBackend for Script {
        fn started(&mut self, row: &Value) -> Result<(), String> {
            self.starts.push(row.clone());
            if self.mode == Mode::FailStart {
                Err("synthetic start durability failure".into())
            } else {
                Ok(())
            }
        }
        fn query(
            &mut self,
            trial: usize,
            public: [u64; 2],
            p: &SatTargetPlan,
        ) -> Result<QueryOutput, String> {
            assert_eq!(p.algorithm_seed, 2026100102);
            self.calls.push(public);
            if self.mode == Mode::Transport {
                return Err("synthetic native transport failure".into());
            }
            if trial == 0 || self.mode == Mode::Timeout {
                return Ok(QueryOutput {
                    native: NativeOutput {
                        exit_code: if self.mode == Mode::Timeout {
                            None
                        } else {
                            Some(15)
                        },
                        timed_out: self.mode == Mode::Timeout,
                        stdout: if self.mode == Mode::Timeout {
                            String::new()
                        } else {
                            "s INDETERMINATE\n".into()
                        },
                        receipt: json!({"kind":"synthetic-inconclusive-control-no-native-execution"}),
                    },
                    manifest: Value::Null,
                    anf: String::new(),
                    cnf: String::new(),
                    source_receipt: Value::Null,
                });
            }
            assert_eq!(public, [62577, 27783]);
            let mut output = source(self.mode == Mode::Corrupt);
            if self.mode == Mode::WrongManifest {
                output.manifest["target"]["x"] = json!("1");
            }
            Ok(output)
        }
        fn completed(&mut self, row: &Value) -> Result<(), String> {
            self.completions.push(row.clone());
            if self.mode == Mode::FailComplete && row["outcome"] == "verified_point_witness" {
                Err("synthetic completion durability failure".into())
            } else {
                Ok(())
            }
        }
    }
    fn replay(mode: Mode) -> (Value, Script) {
        let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
        let mut backend = Script::new(mode);
        let result = state
            .solve_external_sat(PUBLIC, &plan(), &mut backend)
            .unwrap();
        (result, backend)
    }

    #[test]
    fn derives_all_projected_logs_from_columns_without_imported_projection() {
        let _guard = CLOCK.lock().unwrap();
        let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
        let logs = state.projected_logs().unwrap();
        assert_eq!(logs.len(), 63);
        assert_eq!(logs.iter().flatten().count(), 62);
        for (&p, log) in state.geometry.iter().zip(logs) {
            assert_eq!(
                state.independent.mul(p, 2),
                log.and_then(|k| state.independent.mul(Some(state.independent.g), k.into()))
            );
        }
    }
    #[test]
    fn retained_source_recovers_and_charges_independent_verification() {
        let _guard = CLOCK.lock().unwrap();
        let (result, backend) = replay(Mode::Good);
        assert_eq!(result["verified_scalar"], 24886);
        assert!(result["failure"].is_null());
        assert_eq!(backend.calls.len(), 2);
        assert_eq!(backend.starts.len(), 2);
        assert_eq!(backend.completions.len(), 2);
        assert_eq!(result["attempts"][0]["outcome"], "incomplete");
        assert_eq!(result["attempts"][1]["outcome"], "verified_point_witness");
        assert_eq!(
            result["independent_recovery_certificate"]["inside_online_interval"],
            true
        );
        let phases = result["exclusive_ic_online_phases_ns"].as_object().unwrap();
        assert_eq!(phases.len(), 5);
        assert_eq!(
            phases.values().map(|v| v.as_u64().unwrap()).sum::<u64>(),
            result["costs"]["online_wall_ns"].as_u64().unwrap()
        );
        assert_eq!(result["source_bound_execution_admitted"], false);
        assert_eq!(result["fresh_paired_qualification"], false);
        assert!(result["online_speedup"].is_null());
        assert!(result["candidate_id"].is_null());
    }
    #[test]
    fn rejects_assignment_and_manifest_corruption_with_all_attempts_retained() {
        let _guard = CLOCK.lock().unwrap();
        for (mode, outcome) in [
            (Mode::Corrupt, "invalid_source_model"),
            (Mode::WrongManifest, "invalid_source_output"),
        ] {
            let (result, backend) = replay(mode);
            assert!(result["verified_scalar"].is_null());
            assert_eq!(result["attempts"][1]["outcome"], outcome);
            assert_eq!(backend.completions.len(), 2);
            assert!(result["failure"].as_str().is_some());
        }
    }
    #[test]
    fn timeouts_remain_rows_and_exhaustion_never_becomes_success() {
        let _guard = CLOCK.lock().unwrap();
        let (result, backend) = replay(Mode::Timeout);
        assert!(result["verified_scalar"].is_null());
        assert_eq!(backend.completions.len(), 2);
        for row in result["attempts"].as_array().unwrap() {
            assert_eq!(row["outcome"], "timeout");
            assert_eq!(row["source_model_valid"], false);
        }
        assert!(result["failure"]
            .as_str()
            .unwrap()
            .contains("budget exhausted"));
        assert!(result["exclusive_ic_online_phases_ns"]["target_recovery_check"].is_null());
    }
    #[test]
    fn durable_start_completion_and_transport_failures_prevent_success() {
        let _guard = CLOCK.lock().unwrap();
        for (mode, calls, rows) in [
            (Mode::FailStart, 0, 1),
            (Mode::Transport, 1, 1),
            (Mode::FailComplete, 2, 2),
        ] {
            let (result, backend) = replay(mode);
            assert!(result["verified_scalar"].is_null());
            assert_eq!(backend.calls.len(), calls);
            assert_eq!(result["attempts"].as_array().unwrap().len(), rows);
            assert!(result["failure"].as_str().is_some());
            assert!(result["costs"]["online_wall_ns"].as_u64().is_some());
        }
    }
    #[test]
    fn lift_helper_rejects_a_different_group_query_without_source_unsat_claim() {
        let _guard = CLOCK.lock().unwrap();
        let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
        let bits = model();
        sat::verify_anf(ANF, &bits[..51]).unwrap();
        sat::verify_cnf(CNF, &bits).unwrap();
        assert!(state
            .lift_sat_model(
                &bits,
                independent_point(&state.independent, [62577, 27783]).unwrap()
            )
            .is_some());
        assert!(state
            .lift_sat_model(&bits, Some(state.independent.g))
            .is_none());
        // The retained assignment belongs only to its original source instance;
        // this helper control does not assert a new valid nonlifting source.
    }
    #[test]
    fn invalid_public_targets_never_call_backend_and_keep_failure_timing() {
        let _guard = CLOCK.lock().unwrap();
        let state = PreparedN17Target::from_ordinary_math(&retained_math()).unwrap();
        for public in [[1 << 17, 0], [0, 0], [0, 1]] {
            let mut backend = Script::new(Mode::Good);
            let result = state
                .solve_external_sat(public, &plan(), &mut backend)
                .unwrap();
            assert!(backend.calls.is_empty());
            assert!(backend.starts.is_empty());
            assert!(result["verified_scalar"].is_null());
            assert!(result["failure"].as_str().is_some());
        }
    }
    #[test]
    fn plan_rejects_unbounded_search_and_known_scalar_inputs() {
        assert!(plan().validate().is_ok());
        let mut bad = plan();
        bad.max_queries = 9;
        assert!(bad.validate().is_err());
        bad = plan();
        bad.solver_timeout_ms = 0;
        assert!(bad.validate().is_err());
        let mut matched_panel = plan();
        matched_panel.solver_timeout_ms = 120_000;
        assert!(matched_panel.validate().is_ok());
        matched_panel.solver_timeout_ms = 120_001;
        assert!(matched_panel.validate().is_err());
        let mut value = serde_json::to_value(plan()).unwrap();
        value["known_scalar"] = json!(24886);
        assert!(serde_json::from_value::<SatTargetPlan>(value).is_err());
    }
}
