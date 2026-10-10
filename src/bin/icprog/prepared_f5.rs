//! Native reconstruction of the retained F5 preparation; no search or new timing.
//! Keep its original outcome vocabulary and Python provenance. The SAT replay
//! source stays unchanged because its exact bytes belong to a frozen audit.
use super::sat_control::native;
use super::{identity, json as records, oracle::Curve, prepared_sat};
use serde_json::{json, Value};
use std::{collections::BTreeSet, fs, path::Path};

pub const CERTIFICATE_SHA: &str =
    "d8c5d6679fe89561606154785fdfba763835bf261dceb08fae193eff9f5636ae";
const STATE_SHA: &str = "edbff76da6442b9f2e5e8235682c9bf1052465f310c765ba8a37d018a60bf107";
const INPUTS_SHA: &str = "6678815b78b9414d81f94b595d3e10ab38f388511fa3a6831891f685e96f833d";
const CERTIFICATE: &str =
    "research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json";

fn require(condition: bool, message: &str) -> Result<(), String> {
    native::require(condition, message)
}
fn digest(value: &Value) -> Result<String, String> {
    identity::sha256(&records::parse(&value.to_string())?)
}
fn keys(value: &Value, expected: &[&str]) -> Result<(), String> {
    let actual = value.as_object().ok_or("expected object")?;
    require(
        actual.len() == expected.len() && expected.iter().all(|key| actual.contains_key(*key)),
        "preparation fields differ",
    )
}
fn scalar(value: &Value) -> Result<u64, String> {
    value
        .as_u64()
        .or_else(|| value.as_str().and_then(|s| s.parse().ok()))
        .ok_or_else(|| "invalid scalar".into())
}
fn power(mut value: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = result * value % modulus;
        }
        value = value * value % modulus;
        exponent >>= 1;
    }
    result
}
fn insert(pivots: &mut [Option<Vec<u64>>], mut row: Vec<u64>, r: u64) -> Result<bool, String> {
    for column in 0..pivots.len() {
        let coefficient = row[column];
        if coefficient == 0 {
            continue;
        }
        if let Some(pivot) = &pivots[column] {
            for (value, &p) in row.iter_mut().zip(pivot) {
                *value = (*value + r - coefficient * p % r) % r;
            }
        } else {
            let inverse = power(coefficient, r - 2, r);
            for value in &mut row {
                *value = *value * inverse % r;
            }
            pivots[column] = Some(row);
            return Ok(true);
        }
    }
    require(
        row.last() == Some(&0),
        "inconsistent ordinary relation matrix",
    )?;
    Ok(false)
}

/// Sealed entry point used by the controller; refuses arbitrary/resealed inputs.
pub fn verify(document: &Value) -> Result<Value, String> {
    require(
        digest(document)? == CERTIFICATE_SHA,
        "external F5 preparation seal differs",
    )?;
    require(
        digest(&document["certificate"]["inputs"])? == INPUTS_SHA,
        "accepted ordinary input seal differs",
    )?;
    verify_mathematics(document)
}

/// Separate mathematics entry point lets tests challenge retained proof claims.
fn verify_mathematics(document: &Value) -> Result<Value, String> {
    require(
        document["schema_version"] == 1
            && document["native_execution"] == false
            && document["fresh_targets_generated"] == 0
            && document["promotion_eligible"] == false
            && document["online_wall_ns"].is_null()
            && document["online_speedup"].is_null()
            && document["provenance"]["family"] == "f5"
            && document["record_sha256"] == STATE_SHA
            && digest(&document["record"])? == STATE_SHA
            && document["state_id"] == "ICP1hedbff76da644",
        "preparation state or measurement scope differs",
    )?;
    let inputs = &document["certificate"]["inputs"];
    keys(inputs, &["fixture", "base", "logs", "attempts"])?;
    let fixture = &inputs["fixture"];
    require(
        fixture["targets"] == json!([])
            && fixture["target_seeds"] == json!([])
            && fixture["target_scalar_constructed"] == false,
        "preparation contains a target input",
    )?;
    let curve = Curve::new(&records::parse(&fixture.to_string())?)?;
    require(
        (
            curve.n,
            curve.a,
            curve.r,
            curve.h,
            curve.modulus,
            curve.g,
            curve.lam,
        ) == (17, 1, 65587, 2, 131081, (43693, 23339), 17184),
        "outside the disclosed synthetic preparation envelope",
    )?;
    let base = inputs["base"]
        .as_array()
        .ok_or("missing geometric base")?
        .iter()
        .map(|p| curve.decode(&records::parse(&p.to_string())?))
        .collect::<Result<Vec<_>, _>>()?;
    require(
        base.len() == 63
            && !base.contains(&None)
            && base.iter().collect::<BTreeSet<_>>().len() == 63
            && document["record"]["factor_base"]["points"] == json!(base),
        "ordered geometry differs",
    )?;
    let (columns, projection) = prepared_sat::projected_columns(&curve, &base)?;
    let usable = base
        .iter()
        .map(|&p| curve.mul(p, curve.h as u128))
        .filter(|p| p.is_some())
        .collect::<BTreeSet<_>>()
        .len();
    require(
        usable == 62
            && columns.len() == 29
            && document["record"]["projection"] == json!(projection),
        "usable points or folded columns differ",
    )?;
    let pairs = prepared_sat::pair_sums(&curve, &base);
    let attempts = inputs["attempts"]
        .as_array()
        .ok_or("missing ordinary attempts")?;
    require(attempts.len() == 216, "ordinary query count differs")?;
    let mut pivots = vec![None; columns.len()];
    let mut rows = Vec::new();
    let mut trajectory = Vec::new();
    let mut seen = BTreeSet::new();
    let (mut dependent, mut duplicates, mut negatives) = (0, 0, 0);
    for (trial, attempt) in attempts.iter().enumerate() {
        keys(attempt, &["trial", "scalar", "outcome", "indices"])?;
        require(
            attempt["trial"].as_u64() == Some(trial as u64),
            "ordinary chronology differs",
        )?;
        let query_scalar = attempt["scalar"]
            .as_u64()
            .ok_or("ordinary scalar is not an integer")?;
        require(query_scalar < curve.r, "ordinary scalar outside subgroup")?;
        let query = curve.mul(Some(curve.g), query_scalar as u128);
        match attempt["outcome"].as_str() {
            Some("witness") => {
                let indices = attempt["indices"]
                    .as_array()
                    .filter(|i| i.len() == 3)
                    .ok_or("missing three-point witness")?
                    .iter()
                    .map(|i| {
                        i.as_u64()
                            .filter(|&i| i < base.len() as u64)
                            .map(|i| i as usize)
                            .ok_or("invalid witness index")
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                require(
                    indices.iter().fold(None, |sum, &i| curve.add(sum, base[i])) == query,
                    "ordinary witness does not readd",
                )?;
                let mut sorted = indices.clone();
                sorted.sort_unstable();
                if !seen.insert((query_scalar, sorted)) {
                    duplicates += 1;
                } else {
                    let mut row = vec![0; columns.len()];
                    for &i in &indices {
                        if let Some((column, coefficient)) = projection[i] {
                            row[column] = (row[column] + coefficient) % curve.r;
                        }
                    }
                    let rhs = curve.h * query_scalar % curve.r;
                    rows.push(
                        json!({"entries":row.iter().enumerate().filter(|(_,v)|**v!=0)
                        .map(|(i,v)|json!([i,v.to_string()])).collect::<Vec<_>>(),
                        "rhs":rhs.to_string(),"scalar":query_scalar,"indices":indices}),
                    );
                    row.push(rhs);
                    if !insert(&mut pivots, row, curve.r)? {
                        dependent += 1;
                    }
                }
            }
            Some("proved_unsat") => {
                require(
                    attempt["indices"].is_null(),
                    "failed ordinary query contributes a relation",
                )?;
                require(
                    !prepared_sat::has_three_sum(&curve, &base, &pairs, query),
                    "claimed ordinary negative has a decomposition",
                )?;
                negatives += 1;
            }
            _ => return Err("unaccepted original F5 outcome".into()),
        }
        trajectory.push(pivots.iter().filter(|p| p.is_some()).count());
    }
    let rank = pivots.iter().filter(|p| p.is_some()).count();
    require(
        rank == 29 && rows.len() == 61 && negatives == 155 && dependent == 32 && duplicates == 0,
        "ordinary yield or full-rank criterion differs",
    )?;
    let mut logs = vec![0; columns.len()];
    for column in (0..columns.len()).rev() {
        let row = pivots[column].as_ref().ok_or("missing pivot")?;
        let sum = (column + 1..columns.len())
            .map(|i| row[i] * logs[i] % curve.r)
            .fold(0, |sum, value| (sum + value) % curve.r);
        logs[column] = (row[columns.len()] + curve.r - sum) % curve.r;
    }
    let retained = inputs["logs"].as_array().ok_or("missing column logs")?;
    require(retained.len() == columns.len(), "column log count differs")?;
    for (i, (&point, &log)) in columns.iter().zip(&logs).enumerate() {
        keys(&retained[i], &["point", "log"])?;
        require(
            curve.decode(&records::parse(&retained[i]["point"].to_string())?)? == point
                && scalar(&retained[i]["log"])? == log
                && curve.mul(Some(curve.g), log as u128) == point,
            "column log failed independent scalar replay",
        )?;
    }
    require(
        document["record"]["column_logs"]
            == json!(columns
                .iter()
                .zip(&logs)
                .map(|(&p, &l)| json!({"point":p,"log":l}))
                .collect::<Vec<_>>()),
        "state column logs differ",
    )?;
    let matrix = json!({"modulus":curve.r.to_string(),"column_points":columns,"columns":columns.len(),
        "rank":rank,"accepted_rows":rows.len(),"duplicate_relations":duplicates,"dependent_relations":dependent,
        "rows_sha256":digest(&json!(rows))?,"rows":rows});
    let proof = json!({"ordinary_query_count":attempts.len(),"ordinary_outcome_mix":{"proved_unsat":negatives,"witness":61},
        "rank_trajectory":trajectory,"matrix":matrix,"all_column_logs_independently_verified":true,
        "preparation_target_count":0,"previous_target_evidence_retained_in_state":false});
    require(
        document["certificate"]["proof"] == proof,
        "ordinary proof snapshot differs",
    )?;
    Ok(
        json!({"status":"PASS_NATIVE_PREPARED_F5_MATHEMATICS","certificate_sha256":CERTIFICATE_SHA,
        "mathematical_state_sha256":STATE_SHA,"ordinary_queries":attempts.len(),
        "original_outcome_mix":proof["ordinary_outcome_mix"],"verified_relations":rows.len(),
        "independently_proved_geometric_negatives":negatives,"dependent_relations":dependent,
        "duplicate_relations":duplicates,"geometric_points":base.len(),"usable_points":usable,
        "folded_columns":columns.len(),"rank":rank,"verified_column_logs":logs.len(),
        "rows_sha256":matrix["rows_sha256"],"rank_trajectory_sha256":digest(&json!(trajectory))?,
        "target_input_present":false,"new_queries":0,"native_solvers_executed":0,
        "fresh_targets_generated":0,"online_wall_ns":null,"online_speedup":null,"promotion_eligible":false,
        "accounting":"historical ordinary preparation; original Python provenance retained; no new native yield or timing"}),
    )
}

/// This replay never opens the old execution archive or starts a native child.
pub fn run(root: &Path, out: &Path) -> Result<String, String> {
    let path = root.join(CERTIFICATE);
    let bytes = native::read(&path, 2 * 1024 * 1024)?;
    let document: Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    let executable = std::env::current_exe().map_err(|e| e.to_string())?;
    let binary_sha = identity::sha256_hex(&fs::read(&executable).map_err(|e| e.to_string())?);
    let mathematics = verify(&document)?;
    require(
        native::read(&path, 2 * 1024 * 1024)? == bytes
            && identity::sha256_hex(&fs::read(&executable).map_err(|e| e.to_string())?)
                == binary_sha,
        "replay inputs or executable changed",
    )?;
    let result = json!({"schema_version":1,"status":"PASS_NATIVE_RETAINED_F5_PREPARATION_REPLAY",
        "preparation_admission":mathematics,"historical_provenance":document["provenance"],
        "certificate_file_sha256":identity::sha256_hex(&bytes),"checker_binary_sha256":binary_sha,
        "checker_binary_unchanged_before_after":true,"checker_source_sha256":identity::sha256_hex(include_bytes!("prepared_f5.rs")),
        "oracle_source_sha256":identity::sha256_hex(include_bytes!("oracle.rs")),
        "projection_source_sha256":identity::sha256_hex(include_bytes!("prepared_sat.rs")),
        "host_architecture":std::env::consts::ARCH,"host_os":std::env::consts::OS,
        "native_children_executed":0,"original_archive_reopened":false,"source_bound_scientific_runtime_admitted":false,
        "fresh_paired_qualification":false,"headline_eligible":false,"candidate_id":null,
        "online_wall_ns":null,"online_speedup":null,"promotion_eligible":false});
    native::save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn retained() -> Value {
        serde_json::from_str(include_str!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json")).unwrap()
    }
    #[test]
    fn reconstructs_all_original_queries_and_full_rank_without_search() {
        let result = verify(&retained()).unwrap();
        assert_eq!(result["verified_relations"], 61);
        assert_eq!(result["independently_proved_geometric_negatives"], 155);
        assert_eq!(result["dependent_relations"], 32);
        assert_eq!(result["verified_column_logs"], 29);
        assert_eq!(result["new_queries"], 0);
    }
    #[test]
    fn witness_as_negative_and_false_witness_fail_group_checks() {
        let mut value = retained();
        value["certificate"]["inputs"]["attempts"][6]["outcome"] = json!("proved_unsat");
        value["certificate"]["inputs"]["attempts"][6]["indices"] = Value::Null;
        assert!(verify_mathematics(&value)
            .unwrap_err()
            .contains("has a decomposition"));
        let mut value = retained();
        value["certificate"]["inputs"]["attempts"][6]["indices"][0] = json!(0);
        assert!(verify_mathematics(&value)
            .unwrap_err()
            .contains("does not readd"));
    }
    #[test]
    fn log_rank_target_chronology_and_outcome_corruption_are_rejected() {
        for (pointer, replacement) in [
            ("/certificate/inputs/logs/0/log", json!(0)),
            ("/certificate/proof/rank_trajectory/0", json!(1)),
            (
                "/certificate/inputs/fixture/targets",
                json!([[52411, 72106]]),
            ),
            ("/certificate/inputs/attempts/1/trial", json!(0)),
            (
                "/certificate/inputs/attempts/0/outcome",
                json!("CONFLICT_BUDGET_INCONCLUSIVE"),
            ),
            ("/certificate/inputs/attempts/0/indices", json!([0, 1, 2])),
        ] {
            let mut value = retained();
            *value.pointer_mut(pointer).unwrap() = replacement;
            assert!(verify_mathematics(&value).is_err(), "accepted {pointer}");
        }
    }
    #[test]
    fn external_seal_rejects_provenance_swaps_and_promotion() {
        for (pointer, replacement) in [
            ("/provenance/archive_sha256", json!("changed")),
            ("/provenance/family", json!("sat")),
            ("/promotion_eligible", json!(true)),
            ("/online_wall_ns", json!(0)),
            ("/record/projection/1/1", json!(0)),
        ] {
            let mut value = retained();
            *value.pointer_mut(pointer).unwrap() = replacement;
            assert_eq!(
                verify(&value).unwrap_err(),
                "external F5 preparation seal differs"
            );
        }
    }
    #[test]
    fn missing_torsion_reordered_geometry_and_omitted_query_fail() {
        let mut value = retained();
        value["certificate"]["inputs"]["base"]
            .as_array_mut()
            .unwrap()
            .remove(0);
        assert!(verify_mathematics(&value)
            .unwrap_err()
            .contains("geometry differs"));
        let mut value = retained();
        value["certificate"]["inputs"]["base"]
            .as_array_mut()
            .unwrap()
            .swap(0, 1);
        assert!(verify_mathematics(&value)
            .unwrap_err()
            .contains("geometry differs"));
        let mut value = retained();
        value["certificate"]["inputs"]["attempts"]
            .as_array_mut()
            .unwrap()
            .pop();
        assert!(verify_mathematics(&value)
            .unwrap_err()
            .contains("query count differs"));
    }
}
