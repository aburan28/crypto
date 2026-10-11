//! Independent general-curve replay of the frozen S3-batch status control.
//! Usage: koblitz_s3_batch_status_replay <experiment-dir> <receipt.json>

use crypto_lib::binary_ecc::curve::{point_add, point_neg, scalar_mul};
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::{json, Value};
use std::fs;
use std::path::Path;

fn read_one(path: &Path) -> Value {
    let body = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut lines = body.lines().filter(|line| !line.is_empty());
    let first = lines.next().expect("empty input");
    assert!(
        lines.next().is_none(),
        "expected one row: {}",
        path.display()
    );
    serde_json::from_str(first).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn read_rows(path: &Path) -> Vec<Value> {
    let body = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    body.lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            serde_json::from_str(line).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        })
        .collect()
}

fn number(row: &Value, key: &str) -> u64 {
    row[key].as_u64().unwrap_or_else(|| panic!("missing {key}"))
}

fn point(value: &Value, n: u32) -> BinaryPoint {
    if value.is_null() {
        return BinaryPoint::Infinity;
    }
    let [x, y]: [u64; 2] = serde_json::from_value(value.clone()).expect("point coordinates");
    BinaryPoint::Affine {
        x: F2mElement::from_biguint(&BigUint::from(x), n),
        y: F2mElement::from_biguint(&BigUint::from(y), n),
    }
}

fn mul(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 * b as u128) % modulus as u128) as u64
}

fn power(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1;
    while exponent > 0 {
        if exponent & 1 != 0 {
            result = mul(result, base, modulus);
        }
        base = mul(base, base, modulus);
        exponent >>= 1;
    }
    result
}

fn solve(rows: &[Vec<u64>], columns: usize, modulus: u64) -> (usize, Vec<u64>) {
    let mut matrix = rows.to_vec();
    assert!(matrix.iter().all(|row| row.len() == columns + 1));
    let mut rank = 0;
    for column in 0..columns {
        let Some(pivot) = (rank..matrix.len()).find(|&i| matrix[i][column] != 0) else {
            continue;
        };
        matrix.swap(rank, pivot);
        let inverse = power(matrix[rank][column], modulus - 2, modulus);
        assert_eq!(mul(matrix[rank][column], inverse, modulus), 1);
        for offset in column..=columns {
            matrix[rank][offset] = mul(matrix[rank][offset], inverse, modulus);
        }
        for index in 0..matrix.len() {
            if index == rank || matrix[index][column] == 0 {
                continue;
            }
            let factor = matrix[index][column];
            for offset in column..=columns {
                let product = mul(factor, matrix[rank][offset], modulus);
                matrix[index][offset] = ((matrix[index][offset] as u128 + modulus as u128
                    - product as u128)
                    % modulus as u128) as u64;
            }
        }
        rank += 1;
    }
    assert!(matrix
        .iter()
        .all(|row| { row[..columns].iter().any(|&value| value != 0) || row[columns] == 0 }));
    let logs = if rank == columns {
        matrix[..columns].iter().map(|row| row[columns]).collect()
    } else {
        Vec::new()
    };
    (rank, logs)
}

struct Replayed {
    base_points: usize,
    rank_attempts: usize,
    rank_failures: usize,
    target_scalar: Option<u64>,
}

fn replay_run(dir: &Path, expected_q: &Value, expected_success: bool) -> Replayed {
    let base = read_one(&dir.join("base.jsonl"));
    let trace = read_rows(&dir.join("rank.jsonl"));
    let target = read_one(&dir.join("targets.jsonl"));
    let summary = read_one(&dir.join("summary.jsonl"));
    let n = number(&base, "n") as u32;
    assert_eq!(n, 41);
    assert_eq!(number(&base, "a"), 0);
    let curve = KoblitzCurve::new(0, n).expect("n41 Koblitz curve");
    let modulus = curve.subgroup_order.to_u64().expect("u64 subgroup");
    let generator = curve.generator();
    let columns = number(&base, "orbit_columns") as usize;
    assert_eq!(number(&base, "subgroup_order"), modulus);
    assert_eq!(
        base["field_modulus_low_terms"],
        json!(curve.curve.irreducible.low_terms)
    );
    assert_eq!(number(&summary, "rank"), columns as u64);
    assert_eq!(number(&summary, "orbit_columns"), columns as u64);
    assert_eq!(summary["query_backend"], "point_sum");
    assert_eq!(summary["base_hash"], base["base_hash"]);
    assert_eq!(number(&summary, "targets"), 1);
    let points_json = base["factor_base_point_coordinates"]
        .as_array()
        .expect("base points");
    let labels: Vec<[u64; 2]> =
        serde_json::from_value(base["factor_base_point_labels"].clone()).expect("base labels");
    let reps = base["factor_base_representatives"]
        .as_array()
        .expect("representatives");
    assert_eq!(reps.len(), columns);
    assert_eq!(labels.len(), points_json.len());
    assert_eq!(number(&base, "factor_base_points") as usize, labels.len());
    let points: Vec<BinaryPoint> = points_json.iter().map(|p| point(p, n)).collect();
    assert_eq!(trace[0]["kind"], "compact_orbit_rank_header");
    assert_eq!(trace[0]["base_hash"], base["base_hash"]);
    assert_eq!(number(&trace[0], "subgroup_order"), modulus);
    let solution = trace.last().expect("rank solution");
    assert_eq!(solution["kind"], "compact_orbit_rank_solution");
    let mut rows = Vec::new();
    let mut failures = 0;
    for (index, attempt) in trace[1..trace.len() - 1].iter().enumerate() {
        assert_eq!(attempt["kind"], "compact_orbit_rank_attempt");
        assert_eq!(number(attempt, "attempt_index") as usize, index);
        let before = solve(&rows, columns, modulus).0;
        assert_eq!(number(attempt, "rank_before") as usize, before);
        let pivot = number(attempt, "pivotless_column") as usize;
        let scalar = number(attempt, "scalar");
        assert!(pivot < columns && scalar > 0 && scalar < modulus);
        let expected = point_add(
            &curve.curve,
            &scalar_mul(&curve.curve, generator, &BigUint::from(scalar)),
            &point_neg(&point(&reps[pivot], n)),
        );
        assert_eq!(point(&attempt["target"], n), expected);
        if attempt["found"] == true {
            let witness: Vec<usize> =
                serde_json::from_value(attempt["point_indices"].clone()).expect("rank witness");
            assert_eq!(witness.len(), 4);
            let sum = witness.iter().fold(BinaryPoint::Infinity, |acc, &i| {
                point_add(&curve.curve, &acc, &points[i])
            });
            assert_eq!(sum, expected);
            let mut row = vec![0; columns + 1];
            for i in witness {
                let [column, coefficient] = labels[i];
                assert!((column as usize) < columns && coefficient < modulus);
                row[column as usize] = (row[column as usize] + coefficient) % modulus;
            }
            row[pivot] = (row[pivot] + 1) % modulus;
            row[columns] = scalar;
            let recorded: Vec<u64> =
                serde_json::from_value(attempt["row"].clone()).expect("rank row");
            assert_eq!(row, recorded);
            rows.push(row);
        } else {
            assert_eq!(attempt["found"], false);
            failures += 1;
        }
        assert_eq!(
            number(attempt, "rank_after") as usize,
            solve(&rows, columns, modulus).0
        );
    }
    let (rank, logs) = solve(&rows, columns, modulus);
    assert_eq!(rank, columns);
    assert_eq!(solution["logs"], json!(logs));
    assert_eq!(number(solution, "attempts") as usize, trace.len() - 2);
    assert_eq!(number(solution, "failures") as usize, failures);
    assert_eq!(number(&summary, "rank_attempts") as usize, trace.len() - 2);
    assert_eq!(number(&summary, "rank_failures") as usize, failures);
    for (base_point, &[column, coefficient]) in points.iter().zip(&labels) {
        assert!(curve.curve.is_on_curve(base_point));
        assert_eq!(
            scalar_mul(
                &curve.curve,
                generator,
                &BigUint::from(mul(logs[column as usize], coefficient, modulus))
            ),
            *base_point
        );
    }
    let q = point(expected_q, n);
    assert!(curve.curve.is_on_curve(&q));
    assert_eq!(point(&target["target"], n), q);
    if expected_success {
        assert_eq!(target["group_verified"], true);
    } else {
        assert_ne!(target["group_verified"], true);
    }
    assert_eq!(number(&target, "exit_code"), u64::from(!expected_success));
    assert_eq!(
        number(&summary, "targets_solved"),
        u64::from(expected_success)
    );
    assert_eq!(
        number(&summary, "targets_failed"),
        u64::from(!expected_success)
    );
    let target_scalar = if expected_success {
        let witness: Vec<usize> =
            serde_json::from_value(target["point_indices"].clone()).expect("target witness");
        assert_eq!(witness.len(), 4);
        let sum = witness.iter().fold(BinaryPoint::Infinity, |acc, &i| {
            point_add(&curve.curve, &acc, &points[i])
        });
        assert_eq!(sum, q);
        let recovered = witness.iter().fold(0, |acc, &i| {
            let [column, coefficient] = labels[i];
            (acc + mul(logs[column as usize], coefficient, modulus)) % modulus
        });
        assert_eq!(recovered, number(&target, "recovered_scalar"));
        assert_eq!(
            scalar_mul(&curve.curve, generator, &BigUint::from(recovered)),
            q
        );
        Some(recovered)
    } else {
        assert!(target["point_indices"].is_null());
        assert!(target["recovered_scalar"].is_null());
        None
    };
    Replayed {
        base_points: points.len(),
        rank_attempts: trace.len() - 2,
        rank_failures: failures,
        target_scalar,
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 3, "usage: <experiment-dir> <receipt.json>");
    let dir = Path::new(&args[1]);
    let failed_q = read_one(&dir.join("failed_pilot_q.jsonl"));
    let success_q = read_one(&dir.join("success_holdout_q.jsonl"));
    let old = dir.join("runs/old_failure");
    let fixed = dir.join("runs/fixed_failure");
    let good = dir.join("runs/fixed_success");
    assert_eq!(
        fs::read(old.join("base.jsonl")).unwrap(),
        fs::read(fixed.join("base.jsonl")).unwrap()
    );
    assert_eq!(
        fs::read(old.join("rank.jsonl")).unwrap(),
        fs::read(fixed.join("rank.jsonl")).unwrap()
    );
    let old_target = read_one(&old.join("targets.jsonl"));
    let fixed_target = read_one(&fixed.join("targets.jsonl"));
    for field in [
        "target",
        "group_verified",
        "point_indices",
        "x_codes",
        "pinned_intermediates",
        "recovered_scalar",
        "probes",
    ] {
        assert_eq!(old_target[field], fixed_target[field], "changed {field}");
    }
    assert_eq!(number(&old_target, "exit_code"), 0);
    assert_eq!(number(&fixed_target, "exit_code"), 1);
    assert_eq!(
        fs::read_to_string(old.join("process_exit.txt"))
            .unwrap()
            .trim(),
        "0"
    );
    assert_eq!(
        fs::read_to_string(fixed.join("process_exit.txt"))
            .unwrap()
            .trim(),
        "1"
    );
    assert_eq!(
        fs::read_to_string(good.join("process_exit.txt"))
            .unwrap()
            .trim(),
        "0"
    );
    let failure = replay_run(&fixed, &failed_q, false);
    let success = replay_run(&good, &success_q, true);
    assert_eq!(failure.base_points, 1_640);
    assert_eq!(failure.rank_attempts, 43);
    assert_eq!(failure.rank_failures, 23);
    assert_eq!(success.base_points, 6_970);
    assert_eq!(success.rank_attempts, 85);
    assert_eq!(success.rank_failures, 0);
    let receipt = json!({
        "kind":"koblitz_s3_batch_target_status_replay",
        "verified":true,
        "old_failed_process_exit":0,
        "fixed_failed_process_exit":1,
        "fixed_success_process_exit":0,
        "old_and_fixed_failure_base_and_rank_byte_equal":true,
        "failed_target_recovered_scalar":failure.target_scalar,
        "success_target_recovered_scalar":success.target_scalar,
        "failed_rank_attempts":failure.rank_attempts,
        "failed_rank_failures":failure.rank_failures,
        "success_rank_attempts":success.rank_attempts,
        "independent_general_curve_base_and_relation_replay":true
    });
    fs::write(&args[2], format!("{receipt}\n")).expect("write receipt");
    println!("status replay PASS");
}
