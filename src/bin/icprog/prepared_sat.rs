//! Independent native reconstruction of the historical SAT preparation.
//! No new relation collection, target, solver dispatch, or performance claim.
use std::collections::{BTreeMap, BTreeSet};

use serde_json::{json, Value};

use super::{
    identity, json as records,
    oracle::{Curve, Point},
};

fn require(condition: bool, message: &str) -> Result<(), String> {
    if condition {
        Ok(())
    } else {
        Err(message.into())
    }
}

fn natural(value: &Value) -> Result<u64, String> {
    value
        .as_u64()
        .or_else(|| value.as_str().and_then(|s| s.parse().ok()))
        .ok_or_else(|| "expected nonnegative integer or decimal integer string".into())
}

fn canonical_sha(value: &Value) -> Result<String, String> {
    identity::sha256(&records::parse(&value.to_string())?)
}

fn power(mut value: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1;
    while exponent != 0 {
        if exponent & 1 == 1 {
            result = result * value % modulus;
        }
        value = value * value % modulus;
        exponent >>= 1;
    }
    result
}

type Projection = Vec<Option<(usize, u64)>>;

pub fn projected_columns(
    curve: &Curve,
    base: &[Point],
) -> Result<(Vec<Point>, Projection), String> {
    let images = base
        .iter()
        .map(|&p| curve.mul(p, curve.h as u128))
        .collect::<Vec<_>>();
    let mut representatives = BTreeSet::new();
    for &image in images.iter().filter(|p| p.is_some()) {
        let mut point = image;
        let mut smallest = image;
        for _ in 0..curve.n {
            smallest = smallest.min(point).min(curve.neg(point));
            point = curve.frob(point);
        }
        representatives.insert(smallest);
    }
    let columns = representatives.into_iter().collect::<Vec<_>>();
    let mut locations = BTreeMap::new();
    for (column, &representative) in columns.iter().enumerate() {
        let (mut point, mut coefficient) = (representative, 1);
        for _ in 0..curve.n {
            for (p, value) in [
                (point, coefficient),
                (curve.neg(point), (curve.r - coefficient) % curve.r),
            ] {
                if let Some(previous) = locations.insert(p, (column, value)) {
                    require(previous == (column, value), "ambiguous orbit coefficient")?;
                }
            }
            point = curve.frob(point);
            coefficient = coefficient * curve.lam % curve.r;
        }
    }
    let projection = images
        .iter()
        .map(|p| {
            if p.is_none() {
                Ok(None)
            } else {
                locations
                    .get(p)
                    .copied()
                    .map(Some)
                    .ok_or_else(|| "missing projected orbit".into())
            }
        })
        .collect::<Result<Vec<_>, String>>()?;
    Ok((columns, projection))
}

/// All unordered pair sums, with replacement and identity sums retained.
pub fn pair_sums(curve: &Curve, base: &[Point]) -> BTreeSet<Point> {
    let mut sums = BTreeSet::new();
    for (index, &left) in base.iter().enumerate() {
        for &right in &base[index..] {
            sums.insert(curve.add(left, right));
        }
    }
    sums
}

pub fn has_three_sum(
    curve: &Curve,
    base: &[Point],
    pairs: &BTreeSet<Point>,
    target: Point,
) -> bool {
    base.iter()
        .any(|&p| pairs.contains(&curve.add(target, curve.neg(p))))
}

fn insert_row(
    pivots: &mut [Option<Vec<u64>>],
    row: Vec<u64>,
    modulus: u64,
) -> Result<bool, String> {
    let mut reduced = row;
    for column in 0..pivots.len() {
        let coefficient = reduced[column];
        if coefficient == 0 {
            continue;
        }
        if let Some(pivot) = &pivots[column] {
            for (value, &p) in reduced.iter_mut().zip(pivot) {
                *value = (*value + modulus - coefficient * p % modulus) % modulus;
            }
        } else {
            let inverse = power(coefficient, modulus - 2, modulus);
            for value in &mut reduced {
                *value = *value * inverse % modulus;
            }
            pivots[column] = Some(reduced);
            return Ok(true);
        }
    }
    require(reduced.last() == Some(&0), "inconsistent relation matrix")?;
    Ok(false)
}

/// Mathematical verification separated from the external seal for mutation controls.
pub fn verify(preparation: &records::J) -> Result<Value, String> {
    let prep: Value = serde_json::from_str(&records::dumps_compact(preparation, false, false))
        .map_err(|e| e.to_string())?;
    let inputs = &prep["certificate"]["inputs"];
    let curve = Curve::new(preparation.at("certificate")?.at("inputs")?.at("fixture")?)?;
    require(
        curve.n == 17 && curve.a == 1 && curve.r == 65587,
        "preparation is outside toy replay envelope",
    )?;
    let base = preparation
        .at("certificate")?
        .at("inputs")?
        .at("base")?
        .as_arr()
        .ok_or("geometric base is not an array")?
        .iter()
        .map(|p| curve.decode(p))
        .collect::<Result<Vec<_>, _>>()?;
    require(
        base.len() == 63
            && !base.contains(&None)
            && base.iter().collect::<BTreeSet<_>>().len() == 63,
        "geometric base is missing, repeated or identity",
    )?;
    let (columns, projection) = projected_columns(&curve, &base)?;
    require(
        columns.len() == 29 && prep["record"]["projection"] == json!(projection),
        "projection differs",
    )?;
    let usable = base
        .iter()
        .map(|&p| curve.mul(p, curve.h as u128))
        .filter(|p| p.is_some())
        .collect::<BTreeSet<_>>()
        .len();
    require(usable == 62, "usable point count differs")?;
    let pairs = pair_sums(&curve, &base);
    let mut pivots = vec![None; columns.len()];
    let mut seen = BTreeSet::new();
    let mut rows = Vec::new();
    let mut trajectory = Vec::new();
    let mut outcomes = BTreeMap::<String, usize>::new();
    let (mut duplicates, mut dependent, mut negatives) = (0, 0, 0);
    let attempts = inputs["attempts"]
        .as_array()
        .ok_or("ordinary attempts are not an array")?;
    require(attempts.len() == 149, "ordinary query count differs")?;
    for (trial, attempt) in attempts.iter().enumerate() {
        require(
            attempt["trial"].as_u64() == Some(trial as u64),
            "ordinary chronology differs",
        )?;
        let scalar = attempt["scalar"]
            .as_u64()
            .ok_or("ordinary scalar is not an integer")?;
        require(scalar < curve.r, "ordinary scalar outside subgroup")?;
        let outcome = attempt["outcome"]
            .as_str()
            .ok_or("missing ordinary outcome")?;
        *outcomes.entry(outcome.into()).or_default() += 1;
        let query = curve.mul(Some(curve.g), scalar as u128);
        match outcome {
            "VALID_POINT_WITNESS" => {
                let indices = attempt["indices"]
                    .as_array()
                    .filter(|i| i.len() == 3)
                    .ok_or("missing ordinary witness")?
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
                if !seen.insert((scalar, sorted)) {
                    duplicates += 1;
                } else {
                    let mut row = vec![0; columns.len()];
                    for &i in &indices {
                        if let Some((column, coefficient)) = projection[i] {
                            row[column] = (row[column] + coefficient) % curve.r;
                        }
                    }
                    let rhs = curve.h * scalar % curve.r;
                    rows.push(
                        json!({"entries":row.iter().enumerate().filter(|(_,v)|**v!=0)
                        .map(|(i,v)|json!([i,v.to_string()])).collect::<Vec<_>>(),
                        "rhs":rhs.to_string(),"scalar":scalar,"indices":indices}),
                    );
                    row.push(rhs);
                    if !insert_row(&mut pivots, row, curve.r)? {
                        dependent += 1;
                    }
                }
            }
            "SOURCE_UNSAT" | "CONFLICT_BUDGET_INCONCLUSIVE" => {
                require(
                    attempt["indices"].is_null(),
                    "failed query contributes a relation",
                )?;
                if outcome == "SOURCE_UNSAT" {
                    require(
                        !has_three_sum(&curve, &base, &pairs, query),
                        "claimed ordinary negative has a decomposition",
                    )?;
                    negatives += 1;
                }
            }
            _ => return Err("unaccepted ordinary outcome".into()),
        }
        trajectory.push(pivots.iter().filter(|p| p.is_some()).count());
    }
    let rank = pivots.iter().filter(|p| p.is_some()).count();
    require(rank == columns.len(), "preparation matrix lacks full rank")?;
    let mut logs = vec![0; columns.len()];
    for column in (0..columns.len()).rev() {
        let pivot = pivots[column].as_ref().ok_or("missing pivot")?;
        let sum = (column + 1..columns.len())
            .map(|i| pivot[i] * logs[i] % curve.r)
            .fold(0, |sum, v| (sum + v) % curve.r);
        logs[column] = (pivot[columns.len()] + curve.r - sum) % curve.r;
    }
    let retained_logs = inputs["logs"]
        .as_array()
        .ok_or("column logs are not an array")?;
    require(
        retained_logs.len() == columns.len(),
        "column log count differs",
    )?;
    for (index, (&point, &log)) in columns.iter().zip(&logs).enumerate() {
        require(
            retained_logs[index]["point"] == json!(point)
                && natural(&retained_logs[index]["log"])? == log
                && curve.mul(Some(curve.g), log as u128) == point,
            "column log failed independent scalar replay",
        )?;
    }
    let matrix = json!({"modulus":curve.r.to_string(),"column_points":columns,"columns":columns.len(),
        "rank":rank,"accepted_rows":rows.len(),"duplicate_relations":duplicates,"dependent_relations":dependent,
        "rows_sha256":canonical_sha(&json!(rows))?,"rows":rows});
    let proof = &prep["certificate"]["proof"];
    require(
        proof["matrix"] == matrix
            && proof["rank_trajectory"] == json!(trajectory)
            && proof["ordinary_outcome_mix"] == json!(outcomes)
            && proof["ordinary_query_count"] == attempts.len(),
        "preparation proof snapshot differs",
    )?;
    require(
        prep["record"]["column_logs"]
            == json!(columns
                .iter()
                .zip(logs)
                .map(|(&p, l)| json!({"point":p,"log":l}))
                .collect::<Vec<_>>()),
        "mathematical state column logs differ",
    )?;
    Ok(
        json!({"status":"PASS_NATIVE_PREPARED_SAT_MATHEMATICS","geometric_points":base.len(),
        "usable_points":usable,"folded_columns":columns.len(),"rank":rank,"verified_column_logs":columns.len(),
        "ordinary_queries":attempts.len(),"ordinary_outcome_mix":outcomes,"verified_relations":rows.len(),
        "independently_proved_ordinary_negatives":negatives,"dependent_relations":dependent,"duplicate_relations":duplicates,
        "rows_sha256":matrix["rows_sha256"],"rank_trajectory_sha256":canonical_sha(&json!(trajectory))?,
        "preparation_source_sha256":identity::sha256_hex(include_bytes!("prepared_sat.rs")),
        "new_queries":0,"native_solvers_executed":0,"fresh_targets_generated":0,
        "accounting":"historical preparation only; original Python controller provenance retained; no native timings measured"}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    fn retained() -> Value {
        serde_json::from_str(include_str!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/sat-preparation.json")).unwrap()
    }
    fn check(value: &Value) -> Result<Value, String> {
        verify(&records::parse(&value.to_string()).unwrap())
    }
    #[test]
    fn reconstructs_all_retained_rows_rank_logs_and_negatives() {
        let report = check(&retained()).unwrap();
        assert_eq!(report["verified_relations"], 37);
        assert_eq!(report["independently_proved_ordinary_negatives"], 106);
        assert_eq!(
            report["ordinary_outcome_mix"]["CONFLICT_BUDGET_INCONCLUSIVE"],
            6
        );
    }
    #[test]
    fn rejects_changed_witness_projection_log_and_rank_trajectory() {
        for pointer in [
            "/certificate/inputs/attempts/0/indices/0",
            "/record/projection/1/1",
            "/certificate/inputs/logs/0/log",
            "/certificate/proof/rank_trajectory/0",
        ] {
            let mut value = retained();
            *value.pointer_mut(pointer).unwrap() = json!(0);
            assert!(check(&value).is_err(), "accepted mutation at {pointer}");
        }
    }
    #[test]
    fn does_not_accept_a_witness_as_unsat_or_budget_as_a_relation() {
        let mut value = retained();
        value["certificate"]["inputs"]["attempts"][0]["outcome"] = json!("SOURCE_UNSAT");
        value["certificate"]["inputs"]["attempts"][0]["indices"] = Value::Null;
        assert!(check(&value).unwrap_err().contains("has a decomposition"));
        let mut value = retained();
        let attempt = value["certificate"]["inputs"]["attempts"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|a| a["outcome"] == "CONFLICT_BUDGET_INCONCLUSIVE")
            .unwrap();
        attempt["indices"] = json!([0, 1, 2]);
        assert!(check(&value)
            .unwrap_err()
            .contains("contributes a relation"));
    }
    #[test]
    fn pair_complement_proof_retains_repetition_torsion_and_identity_sums() {
        let curve = Curve::new(
            &records::parse(&retained()["certificate"]["inputs"]["fixture"].to_string()).unwrap(),
        )
        .unwrap();
        let p = Some(curve.g);
        let base = vec![p, curve.neg(p), Some((0, 1))];
        let sums = pair_sums(&curve, &base);
        assert!(sums.contains(&None));
        assert!(has_three_sum(
            &curve,
            &base,
            &sums,
            curve.add(curve.add(p, p), p)
        ));
        assert!(has_three_sum(&curve, &base, &sums, Some((0, 1))));
    }
}
