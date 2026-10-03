//! Independent replay of the target-blind n37 rank transcript. This uses
//! the frozen point-defined base and the general binary-curve group law;
//! it never calls the producer's base builder, folded oracle or matrix.

use crypto_lib::binary_ecc::curve::{point_add, scalar_mul};
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

const R: u64 = 230_603_167;
const K: usize = 42;
const SOURCE: &str = "research/notes/ecc2k130/n37_four_policy_support_20261003/SOURCE42.jsonl";
const SOURCE_SHA: &str = "0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75";

fn field_point(value: &Value) -> Result<BinaryPoint, String> {
    let words: [u64; 2] = serde_json::from_value(value.clone())
        .map_err(|error| format!("point coordinates: {error}"))?;
    Ok(BinaryPoint::Affine {
        x: F2mElement::from_biguint(&BigUint::from(words[0]), 37),
        y: F2mElement::from_biguint(&BigUint::from(words[1]), 37),
    })
}

fn u64_field(value: &Value, name: &str) -> Result<u64, String> {
    value[name]
        .as_u64()
        .ok_or_else(|| format!("missing {name}"))
}

fn mulmod(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % R as u128) as u64
}

fn submod(a: u64, b: u64) -> u64 {
    ((a as u128 + R as u128 - b as u128) % R as u128) as u64
}

fn powmod(mut x: u64, mut e: u64) -> u64 {
    let mut out = 1u64;
    while e > 0 {
        if e & 1 == 1 {
            out = mulmod(out, x);
        }
        x = mulmod(x, x);
        e >>= 1;
    }
    out
}

/// Rebuild reduced echelon form from the entire accepted-row prefix.
/// This is deliberately separate from the producer's incremental solver.
fn solve_prefix(rows: &[(Vec<u64>, u64)]) -> Result<(usize, Option<Vec<u64>>), String> {
    let mut matrix: Vec<Vec<u64>> = rows
        .iter()
        .map(|(row, rhs)| {
            let mut augmented = row.clone();
            augmented.push(*rhs);
            augmented
        })
        .collect();
    let mut rank = 0usize;
    for col in 0..K {
        let Some(pivot) = (rank..matrix.len()).find(|&i| matrix[i][col] != 0) else {
            continue;
        };
        matrix.swap(rank, pivot);
        let inverse = powmod(matrix[rank][col], R - 2);
        if mulmod(matrix[rank][col], inverse) != 1 {
            return Err(format!("noninvertible replay pivot at column {col}"));
        }
        for j in col..=K {
            matrix[rank][j] = mulmod(matrix[rank][j], inverse);
        }
        for i in 0..matrix.len() {
            if i == rank || matrix[i][col] == 0 {
                continue;
            }
            let factor = matrix[i][col];
            for j in col..=K {
                matrix[i][j] = submod(matrix[i][j], mulmod(factor, matrix[rank][j]));
            }
        }
        rank += 1;
        if rank == matrix.len() {
            break;
        }
    }
    if matrix
        .iter()
        .any(|row| row[..K].iter().all(|&value| value == 0) && row[K] != 0)
    {
        return Err("independent replay found an inconsistent row".into());
    }
    let logs = if rank == K {
        Some(matrix[..K].iter().map(|row| row[K]).collect())
    } else {
        None
    };
    Ok((rank, logs))
}

fn replay(raw: &[u8]) -> Result<Value, String> {
    let report: Value = serde_json::from_slice(raw).map_err(|error| error.to_string())?;
    let source_bytes = std::fs::read(SOURCE).map_err(|error| format!("source base: {error}"))?;
    if sha256_hex(&source_bytes) != SOURCE_SHA {
        return Err("source42 SHA-256 differs from frozen input".into());
    }
    let header: Value = serde_json::from_slice(
        source_bytes
            .split(|&byte| byte == b'\n')
            .next()
            .ok_or("empty source file")?,
    )
    .map_err(|error| error.to_string())?;
    if header["kind"] != "point_defined_factor_base"
        || u64_field(&header, "n")? != 37
        || u64_field(&header, "a")? != 0
        || u64_field(&header, "subgroup_order")? != R
        || u64_field(&header, "orbit_columns")? != K as u64
        || report["curve"] != "icv1-f2m37-tm534059-32aad96b"
    {
        return Err("curve or frozen source metadata mismatch".into());
    }
    let spec = &report["spec"];
    if u64_field(spec, "columns")? != K as u64
        || u64_field(spec, "raw_x_cap")? != 1_000_000
        || u64_field(spec, "rank_seed")? != 202610031137
        || u64_field(spec, "max_trials")? != 1_000_000
    {
        return Err("producer spec differs from preregistration".into());
    }
    let coordinates = header["factor_base_point_coordinates"]
        .as_array()
        .ok_or("source coordinates missing")?;
    let labels = header["factor_base_point_labels"]
        .as_array()
        .ok_or("source labels missing")?;
    if coordinates.len() != 2 * 37 * K || labels.len() != coordinates.len() {
        return Err("source point or label count mismatch".into());
    }
    let mut digest_bytes = b"ic-shared-rank-base-v1\0".to_vec();
    digest_bytes.extend_from_slice(&(K as u64).to_le_bytes());
    let mut points = Vec::with_capacity(coordinates.len());
    let mut parsed_labels = Vec::with_capacity(labels.len());
    for (coordinates, label) in coordinates.iter().zip(labels) {
        let xy: [u64; 2] = serde_json::from_value(coordinates.clone())
            .map_err(|error| format!("source point: {error}"))?;
        let pair: [u64; 2] = serde_json::from_value(label.clone())
            .map_err(|error| format!("source label: {error}"))?;
        if pair[0] >= K as u64 || pair[1] >= R {
            return Err("invalid frozen source label".into());
        }
        digest_bytes.extend_from_slice(&xy[0].to_le_bytes());
        digest_bytes.extend_from_slice(&xy[1].to_le_bytes());
        digest_bytes.push(0);
        digest_bytes.extend_from_slice(&pair[0].to_le_bytes());
        digest_bytes.extend_from_slice(&pair[1].to_le_bytes());
        points.push(field_point(coordinates)?);
        parsed_labels.push(pair);
    }
    let base_sha = sha256_hex(&digest_bytes);
    if report["base_sha256"].as_str() != Some(&base_sha) {
        return Err("producer base differs from frozen source coordinates and labels".into());
    }

    let kc = KoblitzCurve::new(0, 37).ok_or("source curve unavailable")?;
    let h = kc.cofactor.to_u64().ok_or("cofactor outside u64")?;
    let trials = u64_field(&report, "trials")?;
    let transcript = report["relations"].as_array().ok_or("relations missing")?;
    if transcript.len() as u64 != u64_field(&report, "hits")?
        || u64_field(&report, "signed_points")? != points.len() as u64
        || u64_field(&report, "columns")? != K as u64
    {
        return Err("producer hit or base counts differ".into());
    }
    let mut scalars = StdRng::seed_from_u64(202610031137 ^ 0x5348_4152_4544_524b);
    let mut previous_trial = None;
    let mut expected_scalar = 0u64;
    let mut drawn = 0u64;
    let mut rows = Vec::with_capacity(transcript.len());
    for relation in transcript {
        let trial = u64_field(relation, "trial")?;
        if trial >= trials || previous_trial.is_some_and(|previous| trial <= previous) {
            return Err("relation trials are out of order or beyond the cap".into());
        }
        while drawn <= trial {
            expected_scalar = scalars.gen_range(1..R);
            drawn += 1;
        }
        let scalar = u64_field(relation, "scalar")?;
        if scalar != expected_scalar {
            return Err(format!("probe scalar differs at trial {trial}"));
        }
        let probe = scalar_mul(&kc.curve, kc.generator(), &BigUint::from(scalar));
        let reported_probe = field_point(&json!([
            u64_field(relation, "point_x")?,
            u64_field(relation, "point_y")?
        ]))?;
        if probe != reported_probe || !kc.curve.is_on_curve(&probe) {
            return Err(format!("full-point probe differs at trial {trial}"));
        }
        let witness: Vec<usize> = serde_json::from_value(relation["witness"].clone())
            .map_err(|error| format!("witness at trial {trial}: {error}"))?;
        let mut sum = BinaryPoint::Infinity;
        let mut row = vec![0u64; K];
        for &index in &witness {
            let point = points
                .get(index)
                .ok_or("witness index outside source base")?;
            sum = point_add(&kc.curve, &sum, point);
            let [col, coef] = parsed_labels[index];
            row[col as usize] = ((row[col as usize] as u128 + coef as u128) % R as u128) as u64;
        }
        let rhs = mulmod(h % R, scalar);
        let reported_row: Vec<u64> = serde_json::from_value(relation["row"].clone())
            .map_err(|error| format!("row at trial {trial}: {error}"))?;
        if sum != probe || row != reported_row || rhs != u64_field(relation, "rhs")? {
            return Err(format!("witness or matrix row differs at trial {trial}"));
        }
        let previous_rank = solve_prefix(&rows)?.0;
        rows.push((row, rhs));
        let current_rank = solve_prefix(&rows)?.0;
        if u64_field(relation, "rank_after")? != current_rank as u64
            || relation["independent"].as_bool() != Some(current_rank > previous_rank)
        {
            return Err(format!("rank transition differs at trial {trial}"));
        }
        previous_trial = Some(trial);
    }
    let (rank, solved) = solve_prefix(&rows)?;
    if rank != K
        || u64_field(&report, "rank")? != K as u64
        || report["verified"] != true
        || report["exhausted"] != false
        || previous_trial != Some(trials - 1)
    {
        return Err("target-blind rank gate is incomplete".into());
    }
    let logs = solved.ok_or("full-rank matrix has no solution")?;
    let reported_logs: Vec<u64> = serde_json::from_value(report["column_logs"].clone())
        .map_err(|error| format!("column logs: {error}"))?;
    if logs != reported_logs {
        return Err("independent column logs differ".into());
    }
    let mut checked = 0usize;
    for (col, &log) in logs.iter().enumerate() {
        let index = parsed_labels
            .iter()
            .position(|&[c, coef]| c == col as u64 && coef != 0)
            .ok_or("column has no nonzero-coefficient representative")?;
        let coefficient = parsed_labels[index][1];
        let lhs = scalar_mul(
            &kc.curve,
            kc.generator(),
            &BigUint::from(mulmod(coefficient, log)),
        );
        let rhs = scalar_mul(&kc.curve, &points[index], &BigUint::from(h));
        if lhs != rhs {
            return Err(format!("full-point column logarithm {col} differs"));
        }
        checked += 1;
    }
    if report["verification"]["native"]["base_log_columns_checked"].as_u64() != Some(K as u64)
        || report["verification"]["native"]["relation_witnesses_checked"].as_u64()
            != Some(transcript.len() as u64)
    {
        return Err("producer verification counters differ".into());
    }
    if report["rank_search"]["native"]["trials"].as_u64() != Some(trials)
        || report["rank_search"]["native"]["hits"].as_u64()
            != Some(transcript.len() as u64)
        || report["linear_algebra"]["native"]["rows"].as_u64()
            != Some(transcript.len() as u64)
        || report["linear_algebra"]["native"]["rank"].as_u64() != Some(K as u64)
        || report["linear_algebra"]["native"]["dependent_rows"].as_u64()
            != Some(transcript.len() as u64 - K as u64)
    {
        return Err("producer search or rank counters differ".into());
    }
    // The prior independent folded-table gate established these exact
    // source-support counts. They are a frozen ledger identity here, not a
    // second implementation of the table cost model.
    if report["table"]["group_ops"]["adds"].as_u64() != Some(66_822)
        || report["table"]["native"]["pair_table_entries"].as_u64() != Some(64_467)
    {
        return Err("frozen folded-table ledger differs".into());
    }
    let mut phase_sum = 0.0;
    for name in ["base", "table", "rank_search", "linear_algebra", "verification"] {
        phase_sum += report[name]["gae"]
            .as_f64()
            .ok_or_else(|| format!("missing {name} GAE"))?;
    }
    if report["total_gae"].as_f64() != Some(phase_sum) {
        return Err("phase GAE sum differs from reported setup total".into());
    }
    Ok(json!({
        "status":"PASS",
        "raw_sha256":sha256_hex(raw),
        "source_sha256":SOURCE_SHA,
        "base_sha256":base_sha,
        "curve":report["curve"],
        "trials":trials,
        "relations_verified":transcript.len(),
        "rank":rank,
        "column_logs_verified":checked,
    }))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let raw_path = args
        .next()
        .ok_or("usage: n37_shared_rank_replay RAW.json RECEIPT.json")?;
    let receipt_path = args.next().ok_or("receipt path missing")?;
    let raw = std::fs::read(raw_path)?;
    let receipt = match replay(&raw) {
        Ok(receipt) => receipt,
        Err(reason) => json!({"status":"FAIL","raw_sha256":sha256_hex(&raw),"reason":reason}),
    };
    let file = std::fs::File::create(receipt_path)?;
    serde_json::to_writer_pretty(file, &receipt)?;
    if receipt["status"] != "PASS" {
        return Err(receipt["reason"]
            .as_str()
            .unwrap_or("replay failed")
            .to_owned()
            .into());
    }
    Ok(())
}
