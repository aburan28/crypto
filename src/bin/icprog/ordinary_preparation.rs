//! Independent mathematics for target-free synthetic n17 ordinary preparation.
//! Data only: no native invocation admission, solver, timing or promotion.
use super::{
    json as records, oracle::Curve, prepared_sat, sat_control::native, sat_query_law::QueryLaw,
};
use crypto_lib::cryptanalysis::prepared_sat_control::canonical_sha;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
enum Family {
    MatrixF5,
    Cryptominisat,
}
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
#[serde(rename_all = "snake_case")]
enum Stop {
    PanelComplete,
    Interrupted,
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Plan {
    family: Family,
    algorithm_seed: u64,
    planned_queries: usize,
    query_law: String,
}
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
#[serde(rename_all = "snake_case")]
enum Outcome {
    Witness,
    ProvedUnsat,
    Incomplete,
    Unsupported,
    InvalidModel,
    Unresolved,
    Timeout,
    TransportFailure,
    Oom,
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Attempt {
    trial: usize,
    scalar: u64,
    outcome: Outcome,
    #[serde(deserialize_with = "required_nullable")]
    indices: Option<[usize; 3]>,
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Document {
    schema_version: u32,
    question: String,
    plan: Plan,
    fixture: Value,
    geometric_base: Vec<[u64; 2]>,
    attempts: Vec<Attempt>,
    stop: Stop,
    #[serde(deserialize_with = "required_nullable")]
    claimed_column_logs: Option<Vec<u64>>,
}

// Nullable fields must still be present. Serde otherwise treats a missing
// Option as null, erasing whether the producer retained its explicit claim.
fn required_nullable<'de, D, T>(deserializer: D) -> Result<Option<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    Option::<T>::deserialize(deserializer)
}

/// Reproduce the collector's keyed rand-0.8 StdRng without calling its sampler.
fn scalar(seed: u64, trial: usize) -> u64 {
    let key = seed
        ^ 0x5052_4f42_4553_4551
        ^ (trial as u64)
            .wrapping_mul(0x9e37_79b9_7f4a_7c15)
            .rotate_left(17);
    // QueryLaw's public constructor applies the descent domain. Cancel it to
    // obtain the collector's raw keyed PCG/ChaCha12 state, not descent queries.
    QueryLaw::new(key ^ 0x4445_5343_454e_5400).scalar()
}
fn power(mut a: u64, mut e: u64, r: u64) -> u64 {
    let mut out = 1;
    while e != 0 {
        if e & 1 != 0 {
            out = out * a % r;
        }
        a = a * a % r;
        e >>= 1;
    }
    out
}
fn exact_integer(value: &Value, expected: u64) -> bool {
    value
        .as_u64()
        .or_else(|| value.as_str().and_then(|s| s.parse::<u64>().ok()))
        == Some(expected)
}
fn complete_geometry(curve: &Curve) -> Result<BTreeSet<super::oracle::Point>, String> {
    let mut out = BTreeSet::from([Some((0, 1))]);
    for x in 1..64 {
        let inverse = curve.inv(x)?;
        let u = x ^ curve.a ^ curve.fm(inverse, inverse);
        // In odd degree 17, H(u)=sum(u^(2^(2j)),j=0..8) satisfies
        // H(u)^2+H(u)=u+Tr(u). Its equality to u therefore tests whether
        // y=x*H(u), y+x are the complete two roots over this abscissa.
        let mut z = u;
        for _ in 0..8 {
            let squared = curve.fm(z, z);
            z = curve.fm(squared, squared) ^ u;
        }
        if curve.fm(z, z) ^ z == u {
            let y = curve.fm(x, z);
            out.insert(Some((x, y)));
            out.insert(Some((x, y ^ x)));
        }
    }
    Ok(out)
}
struct Matrix {
    pivots: Vec<Option<Vec<u64>>>,
    seen_relations: BTreeSet<(u64, Vec<usize>)>,
    seen_rows: BTreeSet<Vec<u64>>,
    accepted: Vec<Value>,
    duplicates: usize,
    repeated_rows: usize,
    dependent: usize,
}
impl Matrix {
    fn new(columns: usize) -> Self {
        Self {
            pivots: vec![None; columns],
            seen_relations: BTreeSet::new(),
            seen_rows: BTreeSet::new(),
            accepted: vec![],
            duplicates: 0,
            repeated_rows: 0,
            dependent: 0,
        }
    }
    fn rank(&self) -> usize {
        self.pivots.iter().filter(|p| p.is_some()).count()
    }
    fn push(
        &mut self,
        a: u64,
        mut indices: Vec<usize>,
        mut row: Vec<u64>,
        r: u64,
    ) -> Result<(), String> {
        indices.sort_unstable();
        if !self.seen_relations.insert((a, indices)) {
            self.duplicates += 1;
            return Ok(());
        }
        if !self.seen_rows.insert(row.clone()) {
            self.repeated_rows += 1;
        }
        self.accepted
            .push(json!({"coefficients": &row[..row.len()-1], "rhs": row.last()}));
        for column in 0..self.pivots.len() {
            let coefficient = row[column];
            if coefficient == 0 {
                continue;
            }
            if let Some(pivot) = &self.pivots[column] {
                for (value, &p) in row.iter_mut().zip(pivot) {
                    *value = (*value + r - coefficient * p % r) % r;
                }
            } else {
                let inverse = power(coefficient, r - 2, r);
                for value in &mut row {
                    *value = *value * inverse % r;
                }
                self.pivots[column] = Some(row);
                return Ok(());
            }
        }
        native::require(
            row.last() == Some(&0),
            "inconsistent projected relation matrix",
        )?;
        self.dependent += 1;
        Ok(())
    }
    fn solve(&self, r: u64) -> Option<Vec<u64>> {
        if self.rank() != self.pivots.len() {
            return None;
        }
        let mut logs = vec![0; self.pivots.len()];
        for i in (0..logs.len()).rev() {
            let row = self.pivots[i].as_ref()?;
            let sum = (i + 1..logs.len()).fold(0, |v, j| (v + row[j] * logs[j] % r) % r);
            logs[i] = (row[logs.len()] + r - sum) % r;
        }
        Some(logs)
    }
}

pub fn audit(value: &Value) -> Result<Value, String> {
    let input_digest = canonical_sha(value)?; // reject floats before interpretation
    let document: Document = serde_json::from_value(value.clone()).map_err(|e| e.to_string())?;
    native::require(
        document.schema_version == 1
            && document.question == "ordinary-preparation-mathematics-n17-v1",
        "ordinary preparation schema/question differs",
    )?;
    let plan = &document.plan;
    native::require(
        (1..=512).contains(&plan.planned_queries)
            && plan.query_law == "independent-probe-scalar-rand08-v1"
            && document.attempts.len() <= plan.planned_queries,
        "ordinary query law or bounded panel differs",
    )?;
    native::require(
        (document.stop == Stop::PanelComplete) == (document.attempts.len() == plan.planned_queries),
        "panel completion differs from retained chronology",
    )?;
    let fixture = &document.fixture;
    let fields = [
        "degree",
        "curve_a",
        "subgroup_order",
        "group_order",
        "cofactor",
        "generator",
        "lambda",
        "irreducible",
        "targets",
        "target_seeds",
        "target_scalar_constructed",
    ];
    native::require(
        fixture
            .as_object()
            .is_some_and(|o| o.len() == fields.len() && fields.iter().all(|k| o.contains_key(*k)))
            && fixture["degree"] == 17
            && fixture["curve_a"] == 1
            && fixture["targets"] == json!([])
            && fixture["target_seeds"] == json!([])
            && fixture["target_scalar_constructed"] == false,
        "preparation has a target, undeclared fields or another curve",
    )?;
    // Bound all arithmetic inputs before invoking the generic oracle: malformed
    // huge group/cofactor values must reject rather than overflow its checks.
    native::require(
        exact_integer(&fixture["subgroup_order"], 65587)
            && exact_integer(&fixture["group_order"], 131174)
            && exact_integer(&fixture["cofactor"], 2)
            && exact_integer(&fixture["lambda"], 17184)
            && fixture["generator"].as_array().is_some_and(|g| {
                g.len() == 2 && exact_integer(&g[0], 43693) && exact_integer(&g[1], 23339)
            })
            && fixture["irreducible"] == json!({"degree":17,"low_terms":[0,3]}),
        "fixture differs from bounded exact synthetic n17 record",
    )?;
    let curve = Curve::new(&records::parse(&fixture.to_string())?)?;
    native::require(
        (
            curve.n,
            curve.a,
            curve.r,
            curve.h,
            curve.modulus,
            curve.g,
            curve.lam,
        ) == (17, 1, 65587, 2, 131081, (43693, 23339), 17184),
        "outside exact synthetic n17 mathematics",
    )?;
    let base = document
        .geometric_base
        .iter()
        .map(|p| {
            native::require(
                p[0] < 64 && p[1] < 1 << 17,
                "point outside geometric subspace",
            )?;
            curve.decode(&records::parse(&json!(p).to_string())?)
        })
        .collect::<Result<Vec<_>, String>>()?;
    // Reconstruct every fibre independently; a claimed count alone cannot
    // certify completeness. Ordering remains bound by its separate digest.
    native::require(
        base.len() == 63
            && !base.contains(&None)
            && base.iter().copied().collect::<BTreeSet<_>>() == complete_geometry(&curve)?,
        "geometric base is incomplete or duplicated",
    )?;
    let (columns, projection) = prepared_sat::projected_columns(&curve, &base)?;
    let usable = base
        .iter()
        .map(|&p| curve.mul(p, 2))
        .filter(|p| p.is_some())
        .collect::<BTreeSet<_>>()
        .len();
    native::require(
        usable == 62 && columns.len() == 29,
        "usable points or folded columns differ",
    )?;
    let pairs = prepared_sat::pair_sums(&curve, &base);
    let mut matrix = Matrix::new(columns.len());
    let mut mix = BTreeMap::<String, usize>::new();
    let mut trajectory = Vec::new();
    let mut witnesses = 0;
    let mut feasible_inconclusive = 0;
    for (trial, attempt) in document.attempts.iter().enumerate() {
        native::require(
            attempt.trial == trial && attempt.scalar == scalar(plan.algorithm_seed, trial),
            "ordinary chronology or frozen input law differs",
        )?;
        let point = curve.mul(Some(curve.g), attempt.scalar as u128);
        let outcome = serde_json::to_value(attempt.outcome).map_err(|e| e.to_string())?;
        *mix.entry(outcome.as_str().ok_or("invalid outcome")?.into())
            .or_default() += 1;
        match attempt.outcome {
            Outcome::Witness => {
                let indices = attempt.indices.ok_or("witness has no three indices")?;
                native::require(
                    indices.iter().all(|&i| i < base.len()),
                    "witness index outside base",
                )?;
                native::require(
                    indices.iter().fold(None, |s, &i| curve.add(s, base[i])) == point,
                    "ordinary witness does not readd",
                )?;
                witnesses += 1;
                let mut row = vec![0; columns.len()];
                for &i in &indices {
                    if let Some((column, coefficient)) = projection[i] {
                        row[column] = (row[column] + coefficient) % curve.r;
                    }
                }
                row.push(curve.h * attempt.scalar % curve.r);
                matrix.push(attempt.scalar, indices.to_vec(), row, curve.r)?;
            }
            Outcome::ProvedUnsat => {
                native::require(
                    attempt.indices.is_none(),
                    "negative outcome contributes a relation",
                )?;
                native::require(
                    !prepared_sat::has_three_sum(&curve, &base, &pairs, point),
                    "claimed negative has a geometric decomposition",
                )?;
            }
            _ => {
                native::require(
                    attempt.indices.is_none(),
                    "inconclusive outcome contributes a relation",
                )?;
                feasible_inconclusive +=
                    usize::from(prepared_sat::has_three_sum(&curve, &base, &pairs, point));
            }
        }
        trajectory.push(matrix.rank());
    }
    let logs = matrix.solve(curve.r);
    if let Some(values) = &logs {
        for (&point, &log) in columns.iter().zip(values) {
            native::require(
                curve.mul(Some(curve.g), log as u128) == point,
                "recovered factor log fails scalar replay",
            )?;
        }
    }
    if let Some(claimed) = &document.claimed_column_logs {
        native::require(
            logs.as_ref() == Some(claimed),
            "producer log claim lacks full rank or differs from independent solution",
        )?;
    }
    let complete = document.stop == Stop::PanelComplete;
    Ok(
        json!({"schema_version":1,"status":"PASS_ORDINARY_PREPARATION_MATHEMATICS",
        "input_sha256":input_digest,"declared_family":plan.family,"query_law":plan.query_law,
        "planned_queries":plan.planned_queries,"audited_queries":document.attempts.len(),"panel_complete":complete,
        "ordinary_outcome_mix":mix,"verified_witnesses":witnesses,"feasible_inconclusive_queries":feasible_inconclusive,
        "geometric_points":63,"usable_points":usable,"folded_columns":columns.len(),
        "geometric_order_sha256":canonical_sha(&json!(document.geometric_base))?,
        "column_points":columns,"projection":projection,"accepted_rows":matrix.accepted.len(),
        "duplicate_relations":matrix.duplicates,"repeated_projected_rows":matrix.repeated_rows,
        "dependent_relations":matrix.dependent,"rank":matrix.rank(),"rank_trajectory":trajectory,
        "rows_sha256":canonical_sha(&json!(matrix.accepted))?,"independently_recovered_column_logs":logs,
        "mathematical_preparation_complete":complete && logs.is_some() && document.claimed_column_logs.is_some(),
        "source_bound_execution_admitted":false,"native_solvers_executed":0,"new_queries":0,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,
        "full_goal_complete":false,"online_wall_ns":null,"online_speedup":null,
        "accounting":"mathematical input replay only; native execution, source-model custody, costs and rates need separate evidence"}),
    )
}

pub fn run(input: &Path, out: &Path) -> Result<String, String> {
    let before = native::read(input, 8 * 1024 * 1024)?;
    // Deserialize the original bytes as well as the Value. The struct parser
    // rejects duplicate document/plan/attempt keys before Value can erase them.
    let _: Document = serde_json::from_slice(&before).map_err(|e| e.to_string())?;
    let value = serde_json::from_slice(&before).map_err(|e| e.to_string())?;
    let result = audit(&value)?;
    native::require(
        native::read(input, 8 * 1024 * 1024)? == before,
        "ordinary proof input changed during audit",
    )?;
    native::save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn retained() -> Value {
        let old: Value = serde_json::from_str(include_str!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json")).unwrap();
        let input = &old["certificate"]["inputs"];
        let base = input["base"]
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
        json!({"schema_version":1,"question":"ordinary-preparation-mathematics-n17-v1",
            "plan":{"family":"matrix_f5","algorithm_seed":2026093032u64,"planned_queries":216,"query_law":"independent-probe-scalar-rand08-v1"},
            "fixture":input["fixture"],"geometric_base":base,"attempts":input["attempts"],"stop":"panel_complete",
            "claimed_column_logs":input["logs"].as_array().unwrap().iter().map(|v|v["log"].as_str().unwrap().parse::<u64>().unwrap()).collect::<Vec<_>>()})
    }
    #[test]
    fn keyed_query_law_matches_production_without_calling_it_in_audit() {
        for seed in [0, 1, 2026093032, 2026100311, u64::MAX] {
            for trial in [0, 1, 63, 64, 128, 511] {
                assert_eq!(
                    scalar(seed, trial),
                    crypto_lib::cryptanalysis::koblitz_index_calculus::probe_scalar(
                        seed,
                        trial as u64,
                        65587
                    )
                );
            }
        }
    }
    #[test]
    fn retained_rows_logs_and_negative_controls_replay_without_new_yield() {
        let result = audit(&retained()).unwrap();
        assert_eq!(result["verified_witnesses"], 61);
        assert_eq!(result["ordinary_outcome_mix"]["proved_unsat"], 155);
        assert_eq!(result["rank"], 29);
        assert_eq!(result["dependent_relations"], 32);
        assert_eq!(result["duplicate_relations"], 0);
        assert_eq!(result["mathematical_preparation_complete"], true);
        assert_eq!(result["source_bound_execution_admitted"], false);
        assert_eq!(result["new_queries"], 0);
        assert_eq!(result["full_goal_complete"], false);
        assert!(result["online_speedup"].is_null());
    }
    #[test]
    fn changed_law_targets_geometry_witnesses_and_logs_are_rejected() {
        for (pointer, value) in [
            ("/plan/algorithm_seed", json!(2026100311u64)),
            ("/attempts/0/trial", json!(1)),
            ("/fixture/targets", json!([[52411, 72106]])),
            ("/geometric_base/1", json!([0, 1])),
            ("/attempts/6/indices/0", json!(0)),
            ("/claimed_column_logs/0", json!(0)),
            ("/plan/planned_queries", json!(513)),
        ] {
            let mut changed = retained();
            *changed.pointer_mut(pointer).unwrap() = value;
            assert!(audit(&changed).is_err(), "accepted {pointer}");
        }
        let mut false_negative = retained();
        false_negative["attempts"][6]["outcome"] = json!("proved_unsat");
        false_negative["attempts"][6]["indices"] = Value::Null;
        assert!(audit(&false_negative)
            .unwrap_err()
            .contains("geometric decomposition"));
    }
    #[test]
    fn partial_full_rank_panel_stays_incomplete_and_failed_attempts_stay_rows() {
        let mut value = retained();
        value["plan"]["planned_queries"] = json!(512);
        value["stop"] = json!("interrupted");
        let partial = audit(&value).unwrap();
        assert_eq!(partial["rank"], 29);
        assert_eq!(partial["mathematical_preparation_complete"], false);
        value["attempts"][6]["outcome"] = json!("timeout");
        value["attempts"][6]["indices"] = Value::Null;
        value["claimed_column_logs"] = Value::Null;
        let result = audit(&value).unwrap();
        assert_eq!(result["audited_queries"], 216);
        assert_eq!(result["ordinary_outcome_mix"]["timeout"], 1);
        assert_eq!(result["feasible_inconclusive_queries"], 1);
        assert_eq!(result["panel_complete"], false);
        assert_eq!(result["mathematical_preparation_complete"], false);
        value["stop"] = json!("panel_complete");
        assert!(audit(&value).is_err());
    }
    #[test]
    fn raw_duplicates_do_not_create_novel_rank_and_inconsistent_rows_fail() {
        let mut matrix = Matrix::new(2);
        matrix.push(1, vec![0, 1, 2], vec![1, 0, 2], 65587).unwrap();
        matrix.push(1, vec![2, 0, 1], vec![1, 0, 2], 65587).unwrap();
        assert_eq!(matrix.duplicates, 1);
        assert_eq!(matrix.rank(), 1);
        matrix.push(1, vec![3, 4, 5], vec![1, 0, 2], 65587).unwrap();
        assert_eq!(matrix.repeated_rows, 1);
        assert_eq!(matrix.dependent, 1);
        assert!(matrix.push(2, vec![0, 1, 2], vec![1, 0, 3], 65587).is_err());
        assert!(matrix.solve(65587).is_none());
    }
    #[test]
    fn floats_unknown_fields_and_inconclusive_witnesses_do_not_admit() {
        let mut value = retained();
        value["attempts"][0]["scalar"] = json!(6727.0);
        assert!(audit(&value).is_err());
        let mut value = retained();
        value["candidate_id"] = json!("invented");
        assert!(audit(&value).is_err());
        let mut value = retained();
        value["attempts"][6]["outcome"] = json!("incomplete");
        assert!(audit(&value).is_err());
        let mut value = retained();
        value["fixture"]["cofactor"] = json!("170141183460469231731687303715884105727");
        assert!(audit(&value).is_err());
        for pointer in ["/plan", "/attempts/0", "/fixture/irreducible"] {
            let mut value = retained();
            value.pointer_mut(pointer).unwrap()["unknown"] = json!(true);
            assert!(audit(&value).is_err());
        }
        let mut value = retained();
        value.as_object_mut().unwrap().remove("claimed_column_logs");
        assert!(audit(&value).is_err());
        let mut value = retained();
        value["attempts"][0]
            .as_object_mut()
            .unwrap()
            .remove("indices");
        assert!(audit(&value).is_err());
    }

    #[test]
    fn cli_preserves_output_and_rejects_duplicate_struct_keys() {
        let root = std::env::temp_dir().join(format!(
            "ic-ordinary-preparation-cli-{}",
            std::process::id()
        ));
        std::fs::create_dir(&root).unwrap();
        let input = root.join("input.json");
        let out = root.join("audit.json");
        let value = retained();
        std::fs::write(&input, serde_json::to_vec(&value).unwrap()).unwrap();
        run(&input, &out).unwrap();
        let saved = std::fs::read(&out).unwrap();
        assert!(run(&input, &out).is_err());
        assert_eq!(std::fs::read(&out).unwrap(), saved);
        let duplicate = format!(
            "{{\"schema_version\":1,{}",
            &serde_json::to_string(&value).unwrap()[1..]
        );
        std::fs::write(&input, duplicate).unwrap();
        let next = root.join("invalid.json");
        assert!(run(&input, &next).unwrap_err().contains("duplicate field"));
        assert!(!next.exists());
        std::fs::remove_dir_all(&root).unwrap();
    }
}
