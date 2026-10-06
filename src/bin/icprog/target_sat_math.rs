//! Data-only, independent group replay of a prepared n17 SAT target report.
//! Native ANF/CNF/model bytes and execution custody require a separate frozen
//! controller audit; this module never invokes a solver or admits runtime.
use super::{
    json, oracle::Curve, ordinary_preparation, prepared_sat, sat_control::native,
    sat_query_law::QueryLaw, target_math,
};
use crypto_lib::cryptanalysis::{
    prepared_n17_target::SatTargetPlan,
    prepared_sat_control::{canonical_sha, sha256},
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{path::Path, time::Instant};

const SCOPE: &str = "disclosed-n17-native-one-target-cms-v1";

#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Config {
    schema_version: u32,
    question: String,
    target: [u64; 2],
    plan: SatTargetPlan,
    worker_timeout_ms: u64,
}
impl Config {
    fn validate(&self) -> Result<(), String> {
        self.plan.validate()?;
        native::require(
            self.schema_version == 1
                && self.question == SCOPE
                && self.target.iter().all(|&v| v < 1 << 17)
                && (30_000..=900_000).contains(&self.worker_timeout_ms)
                && self.plan.controller_timeout_ms <= self.worker_timeout_ms,
            "outside bounded n17 SAT target mathematical scope",
        )
    }
}
fn keys(value: &Value, fields: &[&str]) -> Result<(), String> {
    native::require(
        value.as_object().is_some_and(|object| {
            object.len() == fields.len() && fields.iter().all(|key| object.contains_key(*key))
        }),
        "SAT target report has missing or undeclared fields",
    )
}
fn number(value: &Value) -> Result<u64, String> {
    value.as_u64().ok_or("missing unsigned integer".into())
}
fn decode(curve: &Curve, point: &Value) -> Result<super::oracle::Point, String> {
    curve.decode(&json::parse(&point.to_string())?)
}
fn inverse(mut value: u64, modulus: u64) -> u64 {
    let mut power = modulus - 2;
    let mut result = 1;
    while power != 0 {
        if power & 1 != 0 {
            result = result * value % modulus;
        }
        value = value * value % modulus;
        power >>= 1;
    }
    result
}
fn provisional(producer: &Value) -> Result<(), String> {
    for key in [
        "source_bound_execution_admitted",
        "fresh_paired_qualification",
        "headline_eligible",
        "promotion_eligible",
        "full_goal_complete",
    ] {
        native::require(producer[key] == false, "SAT producer self-admits a claim")?;
    }
    for key in [
        "candidate_id",
        "workload_id",
        "run_id",
        "online_speedup",
        "operation_counts",
    ] {
        native::require(
            producer[key].is_null(),
            "SAT producer invents an identity or ratio",
        )?;
    }
    Ok(())
}

pub(super) fn verify(preparation: &Value, cfg: &Config, producer: &Value) -> Result<Value, String> {
    cfg.validate()?;
    let prep = ordinary_preparation::audit(preparation)?;
    native::require(
        prep["declared_family"] == "cryptominisat"
            && prep["panel_complete"] == true
            && prep["mathematical_preparation_complete"] == true
            && prep["rank"] == 29
            && prep["usable_points"] == 62
            && prep["folded_columns"] == 29,
        "SAT target lacks a complete independently replayed SAT preparation",
    )?;
    keys(
        producer,
        &[
            "schema_version",
            "question",
            "frontend_contract",
            "plan",
            "target",
            "target_count",
            "preparation_family",
            "preparation_math_sha256",
            "actual_usable_base_points",
            "geometric_base_points",
            "effective_columns",
            "reusable_validation_ns",
            "solver_setup_ns",
            "attempts",
            "verified_scalar",
            "independent_recovery_certificate",
            "failure",
            "costs",
            "exclusive_ic_online_phases_ns",
            "timing_class",
            "candidate_id",
            "workload_id",
            "run_id",
            "source_bound_execution_admitted",
            "fresh_paired_qualification",
            "headline_eligible",
            "promotion_eligible",
            "full_goal_complete",
            "online_speedup",
            "operation_counts",
        ],
    )?;
    native::require(
        producer["schema_version"] == 1
            && producer["question"] == "prepared-public-target-n17-cms-producer-v1"
            && producer["frontend_contract"] == "external-cryptominisat-source-output"
            && producer["plan"] == json!(cfg.plan)
            && producer["target"] == json!(cfg.target)
            && producer["target_count"] == 1
            && producer["preparation_family"] == "cryptominisat"
            && producer["preparation_math_sha256"] == canonical_sha(preparation)?
            && producer["actual_usable_base_points"] == 62
            && producer["geometric_base_points"] == 63
            && producer["effective_columns"] == 29
            && producer["timing_class"] == "ordinary-host-diagnostic",
        "SAT target producer differs from exact preparation, point or plan",
    )?;
    provisional(producer)?;
    number(&producer["reusable_validation_ns"])?;
    number(&producer["solver_setup_ns"])?;

    let curve = Curve::new(&json::parse(&preparation["fixture"].to_string())?)?;
    let target = decode(&curve, &json!(cfg.target))?;
    native::require(
        target.is_some() && curve.mul(target, curve.r as u128).is_none(),
        "SAT target is outside the nonidentity subgroup",
    )?;
    let base = preparation["geometric_base"]
        .as_array()
        .ok_or("missing geometric base")?
        .iter()
        .map(|point| decode(&curve, point))
        .collect::<Result<Vec<_>, _>>()?;
    let (_, projection) = prepared_sat::projected_columns(&curve, &base)?;
    let logs = prep["independently_recovered_column_logs"]
        .as_array()
        .ok_or("missing independently recovered SAT logs")?
        .iter()
        .map(number)
        .collect::<Result<Vec<_>, _>>()?;
    let pairs = prepared_sat::pair_sums(&curve, &base);
    let attempts = producer["attempts"]
        .as_array()
        .filter(|rows| !rows.is_empty() && rows.len() <= cfg.plan.max_queries)
        .ok_or("missing or over-cap SAT target attempts")?;
    let mut law = QueryLaw::new(cfg.plan.algorithm_seed);
    let mut recovered = None;
    let mut fatal = false;
    let mut audited = Vec::new();
    for (trial, row) in attempts.iter().enumerate() {
        native::require(
            recovered.is_none() && !fatal,
            "SAT target continued after stop",
        )?;
        let (a, b) = (law.scalar(), law.scalar());
        let query = curve.add(
            curve.mul(Some(curve.g), a as u128),
            curve.mul(target, b as u128),
        );
        native::require(
            row["trial"] == trial
                && row["a"] == a
                && row["b"] == b
                && row["public_query"] == json!(query.map(|(x, y)| [x, y]))
                && row["source_model_valid"].is_boolean()
                && row["backend_called"].is_boolean(),
            "SAT target attempt differs from independent seeded point law",
        )?;
        let outcome = row["outcome"]
            .as_str()
            .ok_or("missing SAT target outcome")?;
        if query.is_none() {
            native::require(
                outcome == "identity_query"
                    && row["backend_called"] == false
                    && row["source_model_valid"] == false
                    && row["witness_indices"].is_null()
                    && row["candidate_scalar"].is_null(),
                "SAT target claimed native work on an identity query",
            )?;
        } else {
            native::require(
                row["backend_called"] == true,
                "SAT target omitted native query",
            )?;
            match outcome {
                "verified_point_witness" => {
                    native::require(
                        row["native_status"] == "SAT_MODEL" && row["source_model_valid"] == true,
                        "SAT witness lacks source-valid model status",
                    )?;
                    let indices = row["witness_indices"]
                        .as_array()
                        .filter(|v| v.len() == 3)
                        .ok_or("SAT witness lacks three base indices")?
                        .iter()
                        .map(|v| {
                            number(v).and_then(|i| {
                                usize::try_from(i)
                                    .ok()
                                    .filter(|&i| i < base.len())
                                    .ok_or("SAT witness index outside base".into())
                            })
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    native::require(
                        indices.iter().fold(None, |sum, &i| curve.add(sum, base[i])) == query,
                        "SAT witness fails independent full-point addition",
                    )?;
                    let sum = indices.iter().filter_map(|&i| projection[i]).fold(
                        0,
                        |sum, (column, coefficient)| {
                            (sum + coefficient * logs[column] % curve.r) % curve.r
                        },
                    );
                    let denominator = curve.h * b % curve.r;
                    native::require(denominator != 0, "SAT recovery has no inverse")?;
                    let scalar = (sum + curve.r - curve.h * a % curve.r) % curve.r
                        * inverse(denominator, curve.r)
                        % curve.r;
                    native::require(
                        row["candidate_scalar"] == scalar
                            && curve.mul(Some(curve.g), scalar as u128) == target,
                        "SAT scalar fails independent projection or point replay",
                    )?;
                    recovered = Some(scalar);
                }
                "source_unsat_claim" => native::require(
                    row["native_status"] == "SOURCE_UNSAT"
                        && row["source_model_valid"] == false
                        && row["witness_indices"].is_null()
                        && row["candidate_scalar"].is_null()
                        && !prepared_sat::has_three_sum(&curve, &base, &pairs, query),
                    "SAT source UNSAT claim is not a geometric negative",
                )?,
                "nonlifting_model" => native::require(
                    row["native_status"] == "SAT_MODEL"
                        && row["source_model_valid"] == true
                        && row["witness_indices"].is_null()
                        && row["candidate_scalar"].is_null(),
                    "SAT nonlifting model disposition differs",
                )?,
                "timeout" | "incomplete" => native::require(
                    (outcome == "timeout" && row["native_status"] == "TIMEOUT"
                        || outcome == "incomplete"
                            && matches!(
                                row["native_status"].as_str(),
                                Some("CONFLICT_BUDGET_INCONCLUSIVE" | "UNKNOWN_INCONCLUSIVE")
                            ))
                        && row["source_model_valid"] == false
                        && row["witness_indices"].is_null()
                        && row["candidate_scalar"].is_null(),
                    "SAT incomplete attempt was promoted or mislabeled",
                )?,
                "query_transport_failure"
                | "invalid_source_output"
                | "invalid_source_model"
                | "native_error"
                | "independent_recovery_failure" => fatal = true,
                _ => return Err("unknown SAT target disposition".into()),
            }
        }
        audited.push(json!({"trial":trial,"a":a,"b":b,"public_query":query,
            "outcome":outcome,"geometric_negative_checked":outcome=="source_unsat_claim"}));
    }
    native::require(
        producer["verified_scalar"] == json!(recovered),
        "SAT producer scalar differs from independent recovery",
    )?;
    if let Some(scalar) = recovered {
        native::require(
            producer["failure"].is_null()
                && producer["independent_recovery_certificate"]
                    == json!({"check":"independent-shift-and-add-scalar-replay",
                        "generator":[curve.g.0,curve.g.1],"scalar":scalar,
                        "target":cfg.target,"inside_online_interval":true}),
            "SAT producer recovery certificate differs",
        )?;
    } else {
        native::require(
            producer["independent_recovery_certificate"].is_null()
                && producer["failure"]
                    == if fatal {
                        "fatal source, transport, recovery or durability failure"
                    } else {
                        "query budget exhausted without independently verified recovery"
                    }
                && (fatal || attempts.len() == cfg.plan.max_queries),
            "SAT target stopped early or hid a failure",
        )?;
    }
    let timing = target_math::verify_costs(producer)?;
    Ok(
        json!({"schema_version":1,"status":"PASS_SAT_TARGET_MATHEMATICS_RUNTIME_NOT_ADMITTED",
        "target":cfg.target,"preparation_math_sha256":canonical_sha(preparation)?,
        "preparation":prep,"audited_attempts":audited,"verified_recovery":recovered.is_some(),
        "verified_scalar":recovered,"timing":timing,"solver_calls_by_auditor":0,
        "native_source_bytes_verified":false,"source_bound_execution_admitted":false,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,
        "full_goal_complete":false,"online_wall_ns":null,"online_speedup":null}),
    )
}

pub(super) fn run(
    preparation: &Path,
    config: &Path,
    producer: &Path,
    out: &Path,
) -> Result<String, String> {
    let started = Instant::now();
    let inputs = [preparation, config, producer]
        .map(|path| native::read(path, 16 * 1024 * 1024))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
    let result = (|| {
        let cfg: Config =
            serde_json::from_value(target_math::parse(&inputs[1])?).map_err(|e| e.to_string())?;
        verify(
            &target_math::parse(&inputs[0])?,
            &cfg,
            &target_math::parse(&inputs[2])?,
        )
    })();
    for (path, bytes) in [preparation, config, producer].into_iter().zip(&inputs) {
        native::require(
            native::read(path, 16 * 1024 * 1024)? == *bytes,
            "SAT target audit input changed during replay",
        )?;
    }
    let error = result.as_ref().err().cloned();
    let mut receipt = result.unwrap_or_else(|reason| {
        json!({"schema_version":1,"status":"SAT_TARGET_MATHEMATICS_NOT_VERIFIED",
            "error":reason,"verified_recovery":false,"solver_calls_by_auditor":0,
            "native_source_bytes_verified":false,"source_bound_execution_admitted":false,
            "fresh_paired_qualification":false,"headline_eligible":false,
            "promotion_eligible":false,"full_goal_complete":false,
            "online_wall_ns":null,"online_speedup":null})
    });
    receipt["input_sha256"] = json!({"preparation":sha256(&inputs[0]),
        "config":sha256(&inputs[1]),"producer":sha256(&inputs[2])});
    receipt["independent_audit_wall_ns"] = json!(started.elapsed().as_nanos().to_string());
    receipt["audit_scope"] = json!("mathematics only; original native source/model files, worker and execution custody require a frozen audit");
    native::save(out, &receipt)?;
    if let Some(reason) = error {
        return Err(format!(
            "SAT target mathematical audit failed; receipt retained: {reason}"
        ));
    }
    serde_json::to_string_pretty(&receipt).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    const HISTORICAL_F5: &str = include_str!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json");

    #[test]
    fn full_rank_f5_preparation_cannot_be_relabelled_as_sat() {
        let old: Value = serde_json::from_str(HISTORICAL_F5).unwrap();
        let inputs = &old["certificate"]["inputs"];
        let geometry = inputs["base"]
            .as_array()
            .unwrap()
            .iter()
            .map(|point| {
                point
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_str().unwrap().parse::<u64>().unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let preparation = json!({"schema_version":1,"question":"ordinary-preparation-mathematics-n17-v1",
            "plan":{"family":"matrix_f5","algorithm_seed":2026093032u64,"planned_queries":216,
                "query_law":"independent-probe-scalar-rand08-v1"},
            "fixture":inputs["fixture"],"geometric_base":geometry,"attempts":inputs["attempts"],
            "stop":"panel_complete","claimed_column_logs":inputs["logs"].as_array().unwrap()
                .iter().map(|v|v["log"].as_str().unwrap().parse::<u64>().unwrap()).collect::<Vec<_>>()});
        let independently_checked = ordinary_preparation::audit(&preparation).unwrap();
        assert_eq!(independently_checked["rank"], 29);
        assert_eq!(independently_checked["declared_family"], "matrix_f5");
        let cfg = Config {
            schema_version: 1,
            question: SCOPE.into(),
            target: [52411, 72106],
            plan: SatTargetPlan {
                schema_version: 1,
                question: "prepared-public-target-n17-cms-v1".into(),
                algorithm_seed: 2026100302,
                max_queries: 8,
                export_nonce: 2026100521,
                conflict_budget: 1_000_000,
                exporter_timeout_ms: 30_000,
                solver_timeout_ms: 60_000,
                controller_timeout_ms: 900_000,
            },
            worker_timeout_ms: 900_000,
        };
        let err = verify(&preparation, &cfg, &json!({})).unwrap_err();
        assert!(err.contains("complete independently replayed SAT preparation"));
    }
}
