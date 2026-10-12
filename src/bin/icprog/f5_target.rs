//! Independent mathematics for the bounded, disclosed prepared F5 target.
//! Transport/source admission is a separate gate. No solver is called here.
use super::{
    identity, json as records, oracle::Curve, prepared_f5, prepared_sat, sat_control::native,
    sat_query_law::QueryLaw,
};
use serde_json::{json, Value};
use std::path::Path;

const TARGET: [u64; 2] = [52411, 72106];
const STATE: &str = "edbff76da6442b9f2e5e8235682c9bf1052465f310c765ba8a37d018a60bf107";
const PREPARATION: &str = "research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json";
const RETAINED_REPORT: &str = "research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-v3-control-v1/controls/pipeline.stdout";
const RETAINED_REPORT_SHA: &str =
    "57cb2e8e885c5da4912abd32b767b26ac6263175c5ec98f8f20b93fbc8d250ca";
// Linux debug/test ELF files can exceed the native release-child limit.
// Self-identification remains bounded and uses the regular-file/symlink gate.
const CHECKER_BYTE_LIMIT: u64 = 512 * 1024 * 1024;
const PHASES: [&str; 5] = [
    "target_query",
    "target_pdp",
    "target_relation_check",
    "target_descent",
    "recovery_check",
];

fn number(value: &Value) -> Result<u64, String> {
    value
        .as_u64()
        .ok_or_else(|| "expected a nonnegative integer".into())
}
fn require(condition: bool, message: &str) -> Result<(), String> {
    native::require(condition, message)
}
fn power(mut value: u64, mut exponent: u64, r: u64) -> u64 {
    let mut result = 1;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = result * value % r;
        }
        value = value * value % r;
        exponent >>= 1;
    }
    result
}
fn checker_digest(path: &Path) -> Result<String, String> {
    native::read(path, CHECKER_BYTE_LIMIT)
        .map(|bytes| identity::sha256_hex(&bytes))
        .map_err(|error| {
            format!(
                "F5 replay checker {} (file bytes {:?}, limit {CHECKER_BYTE_LIMIT}): {error}",
                path.display(),
                std::fs::symlink_metadata(path).map(|m| m.len()).ok()
            )
        })
}

/// Postexecution mathematics only. The original registration stays consumed;
/// this command never starts the old worker or substitutes for its frozen audit.
pub fn replay(root: &Path, out: &Path) -> Result<String, String> {
    let own = std::env::current_exe().map_err(|e| e.to_string())?;
    let own_before = checker_digest(&own)?;
    let report_path = root.join(RETAINED_REPORT);
    let raw = native::read(&report_path, 16 * 1024 * 1024)?;
    require(
        identity::sha256_hex(&raw) == RETAINED_REPORT_SHA,
        "retained F5 target report seal differs",
    )?;
    let prep = native::load(&root.join(PREPARATION))?;
    let report: Value = serde_json::from_slice(&raw).map_err(|e| e.to_string())?;
    let mut result = verify(&prep, 2026093032, 8, &report)?;
    require(
        checker_digest(&own)? == own_before
            && identity::sha256_hex(&native::read(&report_path, 16 * 1024 * 1024)?)
                == RETAINED_REPORT_SHA
            && native::load(&root.join(PREPARATION))? == prep,
        "replay inputs or checker changed during verification",
    )?;
    result["schema_version"] = json!(1);
    result["scope"] = json!("retained-disclosed-n17-f5-postexecution-mathematics-only");
    result["producer_report_sha256"] = json!(RETAINED_REPORT_SHA);
    result["historical_producer_build"] = report["generic_build"].clone();
    result["historical_provenance"] = json!({"producer":"retained native Rust worker", "controller":"historical Python orchestration", "timings":"unchanged historical uncalibrated diagnostic; no new performance measurement"});
    result["original_registration_consumed_and_closed"] = json!(true);
    result["checker_binary_sha256"] = json!(own_before);
    result["checker_binary_byte_limit"] = json!(CHECKER_BYTE_LIMIT);
    result["checker_binary_unchanged_before_after"] = json!(true);
    result["fresh_targets_generated"] = json!(0);
    native::save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

/// Full target mathematics for a registered seed and cap; source/drain checks
/// must surround this function before any native execution is admitted.
pub fn verify(prep: &Value, seed: u64, cap: u64, report: &Value) -> Result<Value, String> {
    require(
        (1..=8).contains(&cap),
        "target cap outside disclosed control envelope",
    )?;
    let preparation = prepared_f5::verify(prep)?;
    let inputs = &prep["certificate"]["inputs"];
    let curve = Curve::new(&records::parse(&inputs["fixture"].to_string())?)?;
    let target = curve.decode(&records::parse(&json!(TARGET).to_string())?)?;
    require(
        target.is_some() && curve.mul(target, curve.r as u128).is_none(),
        "bad disclosed target",
    )?;
    let base = inputs["base"]
        .as_array()
        .ok_or("missing prepared geometry")?
        .iter()
        .map(|p| curve.decode(&records::parse(&p.to_string())?))
        .collect::<Result<Vec<_>, _>>()?;
    let (_, projection) = prepared_sat::projected_columns(&curve, &base)?;
    let logs = prep["record"]["column_logs"]
        .as_array()
        .ok_or("missing prepared logs")?
        .iter()
        .map(|item| number(&item["log"]))
        .collect::<Result<Vec<_>, _>>()?;
    let string_points = |points: &[Value]| -> Result<Value, String> {
        Ok(json!(points
            .iter()
            .map(|p| {
                let p = curve
                    .decode(&records::parse(&p.to_string())?)?
                    .ok_or("identity in prepared table")?;
                Ok([p.0.to_string(), p.1.to_string()])
            })
            .collect::<Result<Vec<_>, String>>()?))
    };
    let supplied_logs = json!(prep["record"]["column_logs"]
        .as_array()
        .ok_or("missing state logs")?
        .iter()
        .map(|item| Ok(
            json!({"point":string_points(std::slice::from_ref(&item["point"]))?[0],
            "log":number(&item["log"])?.to_string()})
        ))
        .collect::<Result<Vec<Value>, String>>()?);
    let mut fixture = inputs["fixture"].clone();
    fixture["targets"] = json!([[TARGET[0].to_string(), TARGET[1].to_string()]]);
    fixture["target_seeds"] = json!([null]); // Output provenance; input seeds remain [].
    let expected_config = json!({"solver":"f5","linear_algebra":"dense","summands":3,
        "groebner_degree":3,"node_budget":8192,"conflict_budget":100000,"batch_trials":1,
        "max_trials":cap,"collection_window":null,"factor_base":null,
        "factor_base_cube_root":false,"factor_base_orbits":null,"rho_parallel_walks":32,
        "sparse":crypto_lib::cryptanalysis::koblitz_sparse_la::SparseSolveOptions::default()});
    require(
        report["schema_version"] == 1
            && report["query_schema_version"] == 1
            && report["mode"] == "ic"
            && report["fixture"] == fixture
            && report["effective_config"] == expected_config
            && report["effective_factor_base"] == json!({"kind":"standard_subspace","dimension":6})
            && report["factor_base"]
                == string_points(inputs["base"].as_array().ok_or("missing base")?)?
            && report["column_logs"] == supplied_logs
            && report["columns"] == 29
            && report["summands"] == 3
            && report["trials"] == 0
            && report["relations"] == json!([])
            && report["collection_reports"] == json!([])
            && report["solve_attempts"] == 0
            && report["preparation_mode"] == "imported-certified-log-table-v1"
            && report["preparation_mathematical_state_sha256"] == STATE
            && report["reusable_symbolic_template_prepared"] == true
            && report["target_input"] == "supplied_public_point"
            && report["generic_runtime_policy"] == "default-environment-one-rayon-v1",
        "target-only method, input or preparation differs",
    )?;
    let dispatch = &report["descent_dispatch"];
    require(
        dispatch
            == &json!({"strategy":"Groebner","field_kernel":dispatch["field_kernel"],
        "pair_table":false,"query_rule":"seeded-sample","summands":3,"direct_collision":false})
            && matches!(
                dispatch["field_kernel"].as_str(),
                Some("pmull" | "pclmulqdq" | "portable")
            ),
        "target solver dispatch differs",
    )?;
    let solutions = report["solutions"]
        .as_array()
        .filter(|s| s.len() == 1)
        .ok_or("not exactly one target")?;
    let solution = &solutions[0];
    require(solution["index"] == 0, "target solution index differs")?;
    let attempts = solution["attempts"]
        .as_array()
        .ok_or("missing target attempts")?;
    require(
        !attempts.is_empty()
            && attempts.len() as u64 <= cap
            && number(&solution["trials"])? == attempts.len() as u64,
        "target attempt count exceeds cap or chronology",
    )?;
    let mut law = QueryLaw::new(seed);
    let pairs = prepared_sat::pair_sums(&curve, &base);
    let mut audited = Vec::new();
    let mut recovered = None;
    let mut last_relation = Value::Null;
    for (trial, attempt) in attempts.iter().enumerate() {
        require(
            recovered.is_none() && attempt["trial"] == trial as u64,
            "attempts continued after recovery or changed chronology",
        )?;
        let (a, b) = (law.scalar(), law.scalar());
        require(
            attempt["a"] == a && attempt["b"] == b,
            "target query law differs",
        )?;
        let query = curve.add(
            curve.mul(Some(curve.g), a as u128),
            curve.mul(target, b as u128),
        );
        let pdp = &attempt["pdp"];
        let outcome = pdp["outcome"].as_str().ok_or("missing PDP outcome")?;
        let stats = &pdp["stats"]["stats"];
        if outcome != "identity" {
            require(
                pdp["stats"]["family"] == "groebner"
                    && pdp["stats"]["engine"] == json!({"MatrixF5":{"max_degree":3}})
                    && number(&stats["reductions"])? <= 8192
                    && number(&stats["max_degree_built"])? <= 3,
                "undeclared F5 engine or node/degree limit",
            )?;
            for field in [
                "eliminated",
                "infeasible_branches",
                "oversize",
                "propagations",
                "splits",
            ] {
                number(&stats[field])?;
            }
            require(
                stats["exhausted"].is_boolean() && stats["unsupported"].is_boolean(),
                "missing F5 completion flags",
            )?;
        }
        match outcome {
            "witness" => {
                let indices = pdp["points"]
                    .as_array()
                    .filter(|p| p.len() == 3)
                    .ok_or("missing target witness")?
                    .iter()
                    .map(|i| {
                        i.as_u64()
                            .filter(|&i| i < base.len() as u64)
                            .map(|i| i as usize)
                            .ok_or("invalid target witness index")
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                require(
                    indices.iter().fold(None, |sum, &i| curve.add(sum, base[i])) == query,
                    "target witness does not group-readd",
                )?;
                // A witness can survive budget exhaustion, but never an unsupported frontend.
                require(
                    stats["unsupported"] == false,
                    "unsupported frontend claims witness",
                )?;
                let sum = indices
                    .iter()
                    .filter_map(|&i| projection[i])
                    .fold(0, |sum, (c, k)| (sum + k * logs[c] % curve.r) % curve.r);
                let d = (sum + curve.r - curve.h * a % curve.r) % curve.r
                    * power(curve.h * b % curve.r, curve.r - 2, curve.r)
                    % curve.r;
                require(
                    curve.mul(Some(curve.g), d as u128) == target,
                    "independent target scalar replay failed",
                )?;
                recovered = Some(d);
                last_relation = json!({"a":a,"b":b,"points":indices});
            }
            "identity" => {
                require(
                    query.is_none() && pdp["points"].is_null(),
                    "false identity query",
                )?;
                let d = (curve.r - a) * power(b, curve.r - 2, curve.r) % curve.r;
                require(
                    curve.mul(Some(curve.g), d as u128) == target,
                    "identity recovery failed",
                )?;
                recovered = Some(d);
                last_relation = json!({"a":a,"b":b,"points":null});
            }
            "proved_unsat" => {
                require(
                    pdp["points"].is_null()
                        && stats["unsupported"] == false
                        && stats["exhausted"] == false,
                    "incomplete or unsupported PDP became proved negative",
                )?;
                require(
                    !prepared_sat::has_three_sum(&curve, &base, &pairs, query),
                    "claimed target negative has a geometric decomposition",
                )?;
            }
            "incomplete" | "unsupported" => {
                require(
                    pdp["points"].is_null()
                        && stats["exhausted"] == true
                        && stats["unsupported"] == (outcome == "unsupported"),
                    "inconclusive PDP status/flags differ",
                )?;
            }
            _ => return Err("unaccepted F5 PDP outcome".into()),
        }
        audited.push(json!({"trial":trial,"a":a,"b":b,"public_query":query,"outcome":outcome}));
    }
    let complete = recovered.is_some();
    require(
        report["status"] == if complete { "complete" } else { "incomplete" }
            && report["scalar_verified"] == complete
            && report["scalar_replay_included"] == complete
            && solution["relation"] == last_relation
            && solution["recovered"]
                == recovered
                    .map(|d| json!(d.to_string()))
                    .unwrap_or(Value::Null)
            && (complete || attempts.len() as u64 == cap),
        "target recovery or stop rule differs",
    )?;
    let trace = &report["generic_phase_timing"];
    let phases = trace["online_phases_ns"]
        .as_object()
        .ok_or("missing online phases")?;
    require(
        report["generic_phase_policy"] == "exclusive-owner-thread-v1"
            && report["online_timing_schema"] == 2
            && report["reusable_setup_excluded"] == true
            && phases.len() == 6
            && phases.get("rho_solve") == Some(&Value::Null)
            && PHASES.iter().all(|phase| phases.contains_key(*phase)),
        "online phase schema differs",
    )?;
    let mut sum = 0u64;
    let mut missing = Vec::new();
    for phase in PHASES {
        if phases[phase].is_null() && !complete {
            missing.push(phase);
        } else {
            sum = sum
                .checked_add(number(&phases[phase])?)
                .ok_or("online phase overflow")?;
        }
    }
    let attempted = number(&trace["online_wall_ns"])?;
    require(
        attempted > 0
            && sum == attempted
            && report["online_wall_ns"] == attempted
            && number(&report["outer_online_wall_ns"])? <= attempted
            && attempted <= number(&trace["observed_wall_ns"])?,
        "online interval fails phase closure",
    )?;
    Ok(
        json!({"status":"PASS_NATIVE_PREPARED_F5_TARGET_MATHEMATICS","preparation_admission":preparation,
        "target_complete":complete,"recovered_scalar":recovered,"scalar_independently_verified":complete,
        "audited_attempts":audited,"attempted_online_diagnostic_ns":attempted,
        "complete_online_diagnostic_ns":recovered.map(|_|attempted),"online_phases_ns":trace["online_phases_ns"],
        "unmeasured_phase_slots":missing,"ordinary_queries_executed":0,"native_children_executed":0,
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,
        "headline_eligible":false,"promotion_eligible":false,"online_speedup":null,
        "checker_source_sha256":identity::sha256_hex(include_bytes!("f5_target.rs"))}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    fn inputs() -> (Value, Value) {
        (serde_json::from_str(include_str!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json")).unwrap(),
        serde_json::from_str(include_str!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-v3-control-v1/controls/pipeline.stdout")).unwrap())
    }
    #[test]
    fn portable_replay_creates_once_and_never_reopens_the_control() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let dir = std::env::temp_dir().join(format!("f5-target-replay-{}", std::process::id()));
        std::fs::create_dir(&dir).unwrap();
        let out = dir.join("receipt.json");
        let receipt: Value = serde_json::from_str(&replay(root, &out).unwrap()).unwrap();
        assert_eq!(receipt["source_bound_execution_admitted"], false);
        assert_eq!(receipt["original_registration_consumed_and_closed"], true);
        assert_eq!(receipt["fresh_targets_generated"], 0);
        assert!(receipt["online_speedup"].is_null());
        let original = std::fs::read(&out).unwrap();
        assert!(replay(root, &out).is_err());
        assert_eq!(std::fs::read(&out).unwrap(), original);
        std::fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn retained_disclosed_math_reconstructs_all_queries_and_scalar_without_search() {
        let (prep, report) = inputs();
        let result = verify(&prep, 2026093032, 8, &report).unwrap();
        assert_eq!(result["recovered_scalar"], 24886);
        assert_eq!(result["audited_attempts"].as_array().unwrap().len(), 3);
        assert_eq!(result["native_children_executed"], 0);
        assert_eq!(result["source_bound_execution_admitted"], false);
    }
    #[test]
    fn changed_query_witness_log_recovery_engine_and_phase_fail() {
        for (pointer, replacement) in [
            ("/solutions/0/attempts/0/a", json!(1)),
            ("/solutions/0/attempts/2/pdp/points/0", json!(1)),
            ("/solutions/0/recovered", json!("0")),
            (
                "/solutions/0/attempts/0/pdp/stats/engine/MatrixF5/max_degree",
                json!(4),
            ),
            (
                "/solutions/0/attempts/0/pdp/stats/stats/exhausted",
                json!(true),
            ),
            (
                "/generic_phase_timing/online_phases_ns/target_pdp",
                json!(0),
            ),
            ("/column_logs/0/log", json!("0")),
        ] {
            let (prep, mut report) = inputs();
            *report.pointer_mut(pointer).unwrap() = replacement;
            assert!(
                verify(&prep, 2026093032, 8, &report).is_err(),
                "accepted {pointer}"
            );
        }
    }
    #[test]
    fn witness_cannot_be_reclassified_as_proved_negative() {
        let (prep, mut report) = inputs();
        let pdp = &mut report["solutions"][0]["attempts"][2]["pdp"];
        pdp["outcome"] = json!("proved_unsat");
        pdp["points"] = Value::Null;
        assert!(verify(&prep, 2026093032, 8, &report)
            .unwrap_err()
            .contains("has a geometric decomposition"));
    }
    #[test]
    fn failed_target_keeps_unknown_unperformed_phases_and_cannot_claim_scalar() {
        let (prep, mut report) = inputs();
        report["test_fixture"] = json!(true);
        report["solutions"][0]["attempts"]
            .as_array_mut()
            .unwrap()
            .truncate(1);
        report["solutions"][0]["trials"] = json!(1);
        report["solutions"][0]["relation"] = Value::Null;
        report["solutions"][0]["recovered"] = Value::Null;
        report["effective_config"]["max_trials"] = json!(1);
        report["status"] = json!("incomplete");
        report["scalar_verified"] = json!(false);
        report["scalar_replay_included"] = json!(false);
        report["generic_phase_timing"]["online_phases_ns"] = json!({"target_query":1,"target_pdp":2,
            "target_relation_check":3,"target_descent":null,"recovery_check":null,"rho_solve":null});
        report["generic_phase_timing"]["online_wall_ns"] = json!(6);
        report["generic_phase_timing"]["observed_wall_ns"] = json!(7);
        report["online_wall_ns"] = json!(6);
        report["outer_online_wall_ns"] = json!(5);
        let result = verify(&prep, 2026093032, 1, &report).unwrap();
        assert!(result["recovered_scalar"].is_null());
        assert!(result["complete_online_diagnostic_ns"].is_null());
        assert_eq!(
            result["unmeasured_phase_slots"],
            json!(["target_descent", "recovery_check"])
        );
        assert_eq!(result["target_complete"], false);
        report["scalar_verified"] = json!(true);
        assert!(verify(&prep, 2026093032, 1, &report).is_err());
    }
    #[test]
    fn missing_phase_key_rejects_without_panicking() {
        let (prep, mut report) = inputs();
        let phases = report["generic_phase_timing"]["online_phases_ns"]
            .as_object_mut()
            .unwrap();
        phases.remove("rho_solve");
        phases.insert("unknown_phase".into(), Value::Null);
        assert!(verify(&prep, 2026093032, 8, &report)
            .unwrap_err()
            .contains("phase schema differs"));
    }
}
