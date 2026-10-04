//! Independently replay a compact-orbit one-target run with the general
//! binary-curve group law and a fresh modular Gaussian elimination.
//!
//! Usage: <base.jsonl> <rank.jsonl> <ic_targets.jsonl>
//!        <ic_summary.jsonl> <rho.jsonl> <target_points.jsonl>

use crypto_lib::binary_ecc::curve::{point_add, point_neg, scalar_mul};
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::{json, Value};

fn number(value: &Value, key: &str) -> Result<u64, String> {
    value[key]
        .as_u64()
        .ok_or_else(|| format!("missing or invalid {key}"))
}

fn point(value: &Value, n: u32) -> Result<BinaryPoint, String> {
    if value.is_null() {
        return Ok(BinaryPoint::Infinity);
    }
    let [x, y]: [u64; 2] =
        serde_json::from_value(value.clone()).map_err(|error| error.to_string())?;
    Ok(BinaryPoint::Affine {
        x: F2mElement::from_biguint(&BigUint::from(x), n),
        y: F2mElement::from_biguint(&BigUint::from(y), n),
    })
}

fn read_one(path: &str) -> Result<Value, String> {
    let bytes = std::fs::read(path).map_err(|error| format!("{path}: {error}"))?;
    let rows: Vec<&[u8]> = bytes
        .split(|&byte| byte == b'\n')
        .filter(|line| !line.is_empty())
        .collect();
    if rows.len() != 1 {
        return Err(format!("{path}: expected exactly one JSON row"));
    }
    serde_json::from_slice(rows[0]).map_err(|error| format!("{path}: {error}"))
}

fn read_rows(path: &str) -> Result<Vec<Value>, String> {
    let bytes = std::fs::read(path).map_err(|error| format!("{path}: {error}"))?;
    bytes
        .split(|&byte| byte == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| serde_json::from_slice(line).map_err(|error| format!("{path}: {error}")))
        .collect()
}

fn mulmod(left: u64, right: u64, modulus: u64) -> u64 {
    ((left as u128 * right as u128) % modulus as u128) as u64
}

fn powmod(mut value: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1u64;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = mulmod(result, value, modulus);
        }
        value = mulmod(value, value, modulus);
        exponent >>= 1;
    }
    result
}

fn reduce(rows: &[Vec<u64>], columns: usize, modulus: u64) -> Result<(usize, Vec<u64>), String> {
    let mut matrix = rows.to_vec();
    if matrix.iter().any(|row| row.len() != columns + 1) {
        return Err("rank row width differs from base columns".into());
    }
    let mut rank = 0usize;
    for column in 0..columns {
        let Some(pivot) = (rank..matrix.len()).find(|&i| matrix[i][column] != 0) else {
            continue;
        };
        matrix.swap(rank, pivot);
        let inverse = powmod(matrix[rank][column], modulus - 2, modulus);
        if mulmod(matrix[rank][column], inverse, modulus) != 1 {
            return Err(format!("noninvertible pivot in column {column}"));
        }
        for offset in column..=columns {
            matrix[rank][offset] = mulmod(matrix[rank][offset], inverse, modulus);
        }
        for row in 0..matrix.len() {
            if row == rank || matrix[row][column] == 0 {
                continue;
            }
            let factor = matrix[row][column];
            for offset in column..=columns {
                let subtrahend = mulmod(factor, matrix[rank][offset], modulus);
                matrix[row][offset] = ((matrix[row][offset] as u128 + modulus as u128
                    - subtrahend as u128)
                    % modulus as u128) as u64;
            }
        }
        rank += 1;
    }
    if matrix
        .iter()
        .any(|row| row[..columns].iter().all(|&v| v == 0) && row[columns] != 0)
    {
        return Err("inconsistent rank relation".into());
    }
    let solution = if rank == columns {
        matrix[..columns].iter().map(|row| row[columns]).collect()
    } else {
        Vec::new()
    };
    Ok((rank, solution))
}

fn replay(args: &[String]) -> Result<Value, String> {
    let base = read_one(&args[1])?;
    let trace = read_rows(&args[2])?;
    let ic = read_one(&args[3])?;
    let summary = read_one(&args[4])?;
    let rho = read_one(&args[5])?;
    let target_input = read_one(&args[6])?;
    let n = number(&base, "n")? as u32;
    let a = number(&base, "a")? as u8;
    let curve = KoblitzCurve::new(a, n).ok_or("curve unavailable")?;
    let modulus = curve
        .subgroup_order
        .to_u64()
        .ok_or("subgroup order outside u64")?;
    if trace.len() < 2
        || trace[0]["kind"] != "compact_orbit_rank_header"
        || number(&trace[0], "n")? != n as u64
        || number(&trace[0], "a")? != a as u64
        || number(&trace[0], "subgroup_order")? != modulus
        || trace[0]["base_hash"] != base["base_hash"]
    {
        return Err("rank trace header differs from the frozen base".into());
    }
    let columns = number(&base, "orbit_columns")? as usize;
    let coordinates = base["factor_base_point_coordinates"]
        .as_array()
        .ok_or("missing factor-base coordinates")?;
    let labels = base["factor_base_point_labels"]
        .as_array()
        .ok_or("missing factor-base labels")?;
    let representatives = base["factor_base_representatives"]
        .as_array()
        .ok_or("missing factor-base representatives")?;
    if number(&base, "subgroup_order")? != modulus
        || base["field_modulus_low_terms"] != json!(curve.curve.irreducible.low_terms)
        || rho["field_modulus_low_terms"] != base["field_modulus_low_terms"]
        || representatives.len() != columns
        || coordinates.len() != labels.len()
        || number(&base, "factor_base_points")? as usize != coordinates.len()
        || number(&summary, "rank")? as usize != columns
        || number(&summary, "factor_base_points")? as usize != coordinates.len()
        || number(&summary, "orbit_columns")? as usize != columns
        || number(&summary, "targets")? != 1
        || number(&summary, "targets_solved")? != 1
        || number(&summary, "targets_failed")? != 0
        || summary["base_hash"] != base["base_hash"]
    {
        return Err("curve, base, or final-rank metadata differs".into());
    }
    let points: Vec<BinaryPoint> = coordinates
        .iter()
        .map(|value| point(value, n))
        .collect::<Result<_, _>>()?;
    let parsed_labels: Vec<[u64; 2]> = labels
        .iter()
        .map(|value| serde_json::from_value(value.clone()).map_err(|error| error.to_string()))
        .collect::<Result<_, _>>()?;
    if parsed_labels
        .iter()
        .any(|&[column, coefficient]| column as usize >= columns || coefficient >= modulus)
    {
        return Err("factor-base label outside matrix or subgroup".into());
    }
    let generator = curve.generator();
    let mut rows = Vec::<Vec<u64>>::new();
    let mut failures = 0usize;
    let mut verified_relations = 0usize;
    for (index, attempt) in trace[1..trace.len() - 1].iter().enumerate() {
        if attempt["kind"] != "compact_orbit_rank_attempt"
            || number(attempt, "attempt_index")? != index as u64
        {
            return Err("rank attempt order or kind differs".into());
        }
        let rank_before = reduce(&rows, columns, modulus)?.0;
        if number(attempt, "rank_before")? as usize != rank_before {
            return Err(format!("rank before attempt {index} differs"));
        }
        let scalar = number(attempt, "scalar")?;
        let pivot = number(attempt, "pivotless_column")? as usize;
        if pivot >= columns || scalar == 0 || scalar >= modulus {
            return Err(format!("rank query {index} outside frozen group or base"));
        }
        let expected = point_add(
            &curve.curve,
            &scalar_mul(&curve.curve, generator, &BigUint::from(scalar)),
            &point_neg(&point(&representatives[pivot], n)?),
        );
        let target = point(&attempt["target"], n)?;
        if expected != target || !curve.curve.is_on_curve(&target) {
            return Err(format!("rank target {index} failed general-law replay"));
        }
        if attempt["found"] == true {
            let witness: Vec<usize> = serde_json::from_value(attempt["point_indices"].clone())
                .map_err(|error| error.to_string())?;
            if witness.len() != 4 || witness.iter().any(|&i| i >= points.len()) {
                return Err(format!("rank witness {index} invalid"));
            }
            let sum = witness.iter().fold(BinaryPoint::Infinity, |acc, &i| {
                point_add(&curve.curve, &acc, &points[i])
            });
            if sum != target {
                return Err(format!("rank witness {index} does not sum to target"));
            }
            let mut row = vec![0u64; columns + 1];
            for &i in &witness {
                let [column, coefficient] = parsed_labels[i];
                row[column as usize] = (row[column as usize] + coefficient) % modulus;
            }
            row[pivot] = (row[pivot] + 1) % modulus;
            row[columns] = scalar;
            let recorded: Vec<u64> = serde_json::from_value(attempt["row"].clone())
                .map_err(|error| error.to_string())?;
            if row != recorded {
                return Err(format!("rank row {index} differs"));
            }
            rows.push(row);
            verified_relations += 1;
        } else if attempt["found"] == false {
            failures += 1;
        } else {
            return Err(format!("rank status {index} missing"));
        }
        let rank_after = reduce(&rows, columns, modulus)?.0;
        if number(attempt, "rank_after")? as usize != rank_after
            || attempt["gained"]
                .as_bool()
                .is_some_and(|gained| gained != (rank_after > rank_before))
        {
            return Err(format!("rank after attempt {index} differs"));
        }
    }
    let solution = trace.last().ok_or("empty rank trace")?;
    let (rank, logs) = reduce(&rows, columns, modulus)?;
    let recorded_logs: Vec<u64> =
        serde_json::from_value(solution["logs"].clone()).map_err(|error| error.to_string())?;
    if rank != columns
        || logs != recorded_logs
        || number(solution, "rank")? as usize != rank
        || number(solution, "attempts")? as usize != verified_relations + failures
        || number(&summary, "rank_attempts")? as usize != verified_relations + failures
        || number(&summary, "rank_relations")? as usize != verified_relations
        || number(&summary, "rank_failures")? as usize != failures
    {
        return Err("independent matrix solution or counts differ".into());
    }
    for (index, (base_point, &[column, coefficient])) in
        points.iter().zip(&parsed_labels).enumerate()
    {
        if !curve.curve.is_on_curve(base_point)
            || scalar_mul(
                &curve.curve,
                generator,
                &BigUint::from(mulmod(logs[column as usize], coefficient, modulus)),
            ) != *base_point
        {
            return Err(format!("base log {index} failed independent scalar replay"));
        }
    }
    let q = point(&target_input, n)?;
    if q != point(&ic["target"], n)?
        || q != point(&rho["published_q"], n)?
        || number(&ic, "n")? != n as u64
        || number(&ic, "a")? != a as u64
        || number(&rho, "n")? != n as u64
        || number(&rho, "a")? != a as u64
        || number(&rho, "subgroup_order")? != modulus
        || ic["published_fixture_scalar"] != Value::Null
        || rho["published_fixture_scalar"] != Value::Null
        || rho["reference_grade"] != "strong"
        || rho["quotient_mode"] != "signed_frobenius"
        || rho["target_kind"] != "public_hash_to_curve_cofactor"
        || number(&rho, "public_hash_seed")? != 370413
        || ic["group_verified"] != true
        || rho["verified"] != true
    {
        return Err("the two arms do not solve the same unknown-scalar public point".into());
    }
    let witness: Vec<usize> =
        serde_json::from_value(ic["point_indices"].clone()).map_err(|error| error.to_string())?;
    if witness.len() != 4 || witness.iter().any(|&i| i >= points.len()) {
        return Err("target witness invalid".into());
    }
    let sum = witness.iter().fold(BinaryPoint::Infinity, |acc, &i| {
        point_add(&curve.curve, &acc, &points[i])
    });
    let recovered = witness.iter().fold(0u64, |acc, &i| {
        let [column, coefficient] = parsed_labels[i];
        (acc + mulmod(logs[column as usize], coefficient, modulus)) % modulus
    });
    if sum != q
        || recovered != number(&ic, "recovered_scalar")?
        || recovered != number(&rho, "recovered_fixture_scalar")?
        || scalar_mul(&curve.curve, generator, &BigUint::from(recovered)) != q
    {
        return Err("target relation or recovered scalar failed replay".into());
    }
    let online_ms = ic["online_ms"]
        .as_f64()
        .ok_or("IC online interval missing")?;
    let phases_ms = ic["target_phase_sum_ms"]
        .as_f64()
        .ok_or("IC phase sum missing")?;
    let rss = number(&summary, "peak_rss_bytes")?;
    if (online_ms - phases_ms).abs() > 0.000001 || rss > 16 * 1024 * 1024 * 1024 {
        return Err("online phase sum or memory acceptance gate failed".into());
    }
    Ok(json!({
        "verified":true,
        "n":n,
        "a":a,
        "subgroup_order":modulus,
        "factor_base_points":points.len(),
        "orbit_columns":columns,
        "rank":rank,
        "rank_attempts":verified_relations + failures,
        "rank_relations":verified_relations,
        "rank_failures":failures,
        "replayed_base_logs":points.len(),
        "target_relation_replayed":true,
        "target_scalar":recovered,
        "ic_online_ms":online_ms,
        "rho_online_ms":rho["walk_ms"].as_f64().ok_or("rho walk time missing")?
            + rho["validation_ms"].as_f64().ok_or("rho validation time missing")?,
        "peak_rss_bytes":rss
    }))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 7 {
        eprintln!("usage: <base> <rank> <ic-target> <ic-summary> <rho> <target-point>");
        std::process::exit(2);
    }
    match replay(&args) {
        Ok(receipt) => println!("{receipt}"),
        Err(error) => {
            eprintln!("replay failed: {error}");
            std::process::exit(1);
        }
    }
}
