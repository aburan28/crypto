//! Independent generic-curve replay and accounting for 75% root-table runs.
//! Usage: n83-load75-analyzer <study-dir> <primary_v2|panel>

use std::collections::HashSet;
use std::fs;
use std::path::Path;

use crypto_lib::binary_ecc::curve::point_neg;
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_bigint::BigUint;
use serde_json::{json, Value};

const ARMS: [&str; 3] = ["k1798_dense", "k2000_load75", "k2202_load75"];

fn read_jsonl(path: &Path) -> Result<Vec<Value>, Box<dyn std::error::Error>> {
    fs::read_to_string(path)?
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| serde_json::from_str(line).map_err(Into::into))
        .collect()
}

fn coordinates(value: &Value) -> Result<[BigUint; 2], Box<dyn std::error::Error>> {
    let pair = value.as_array().ok_or("point must be an array")?;
    if pair.len() != 2 {
        return Err("point needs two coordinates".into());
    }
    Ok([
        pair[0]
            .as_str()
            .ok_or("x must be decimal string")?
            .parse()?,
        pair[1]
            .as_str()
            .ok_or("y must be decimal string")?
            .parse()?,
    ])
}

fn point(value: &Value) -> Result<BinaryPoint, Box<dyn std::error::Error>> {
    if value.is_null() {
        return Ok(BinaryPoint::Infinity);
    }
    let [x, y] = coordinates(value)?;
    if x.bits() > 83 || y.bits() > 83 {
        return Err("noncanonical N83 point coordinate".into());
    }
    Ok(BinaryPoint::Affine {
        x: F2mElement::from_biguint(&x, 83),
        y: F2mElement::from_biguint(&y, 83),
    })
}

fn replay_base(
    dump: &Value,
    curve: &KoblitzCurve,
) -> Result<Vec<BinaryPoint>, Box<dyn std::error::Error>> {
    const R: u64 = 8_569_786_107_849_059;
    const LAMBDA: u64 = 2_998_004_823_950_457;
    let coordinates = dump["factor_base_point_coordinates"]
        .as_array()
        .ok_or("missing base point coordinates")?;
    let labels = dump["factor_base_point_labels"]
        .as_array()
        .ok_or("missing base labels")?;
    let representatives = dump["factor_base_representatives"]
        .as_array()
        .ok_or("missing base representatives")?;
    if coordinates.len() != labels.len()
        || coordinates.len() != 166 * representatives.len()
        || dump["factor_base_points"].as_u64() != Some(coordinates.len() as u64)
        || dump["orbit_columns"].as_u64() != Some(representatives.len() as u64)
    {
        return Err("base dimensions do not match complete signed Frobenius orbits".into());
    }
    let points = coordinates
        .iter()
        .map(point)
        .collect::<Result<Vec<_>, _>>()?;
    let mut distinct = HashSet::with_capacity(points.len());
    for (column, rep_value) in representatives.iter().enumerate() {
        let representative = point(rep_value)?;
        if representative == BinaryPoint::Infinity
            || !curve.curve.is_on_curve(&representative)
            || curve.mul(&representative, &curve.subgroup_order) != BinaryPoint::Infinity
        {
            return Err(format!("base representative {column} fails subgroup replay").into());
        }
        let mut current = representative.clone();
        let mut coefficient = 1u64;
        for rotation in 0..83usize {
            let positive = 166 * column + 2 * rotation;
            let negative = positive + 1;
            if points[positive] != current || points[negative] != point_neg(&current) {
                return Err(format!("base orbit mismatch at {column}/{rotation}").into());
            }
            for (index, expected_coefficient) in
                [(positive, coefficient), (negative, (R - coefficient) % R)]
            {
                let label = labels[index].as_array().ok_or("bad base label")?;
                if label.len() != 2
                    || label[0].as_u64() != Some(column as u64)
                    || label[1].as_u64() != Some(expected_coefficient)
                {
                    return Err(format!("base label mismatch at index {index}").into());
                }
                let BinaryPoint::Affine { x, y } = &points[index] else {
                    return Err("base contains infinity".into());
                };
                if !distinct.insert((x.to_biguint(), y.to_biguint())) {
                    return Err(format!("duplicate base point at index {index}").into());
                }
            }
            current = curve.frobenius(&current);
            coefficient = ((u128::from(coefficient) * u128::from(LAMBDA)) % u128::from(R)) as u64;
        }
        if current != representative || coefficient != 1 {
            return Err(format!("base orbit does not close at column {column}").into());
        }
    }
    Ok(points)
}

fn verify_arm(
    study: &Path,
    class: &str,
    arm: &str,
    answers: &[Value],
    curve: &KoblitzCurve,
) -> Result<Value, Box<dyn std::error::Error>> {
    let dir = study.join("runs").join(class).join(arm);
    if !dir.exists() {
        return Ok(json!({"arm": arm, "status": "MISSING"}));
    }
    let exit: i32 = fs::read_to_string(dir.join("exit_status.txt"))?
        .trim()
        .parse()?;
    if exit != 0 {
        return Ok(
            json!({"arm": arm, "status": if exit == 124 { "TIMEOUT" } else { "FAILED" }, "exit_status":exit}),
        );
    }
    let summaries = read_jsonl(&dir.join("stdout.jsonl"))?;
    if summaries.len() != 1 || summaries[0]["kind"] != "compact_orbit_dlp_summary" {
        return Err(format!("{arm}: missing exact summary").into());
    }
    let summary = &summaries[0];
    let expected_load75 = arm.ends_with("_load75");
    if summary["root_table_load75_enabled"] != expected_load75
        || summary["root_table_dense_enabled"] != !expected_load75
        || summary["root_table_allocated_bytes"].as_u64() != Some(8 * 1024 * 1024 * 1024)
    {
        return Err(format!("{arm}: table mode or 8 GiB capacity mismatch").into());
    }
    let base_points = if class == "panel" {
        let dumps = read_jsonl(
            &study
                .join("runs")
                .join("panel")
                .join(arm)
                .join("base_dump.jsonl"),
        )?;
        if dumps.len() != 1
            || dumps[0]["kind"] != "point_defined_factor_base"
            || dumps[0]["base_hash"] != summary["base_hash"]
            || dumps[0]["factor_base_points"] != summary["factor_base_points"]
        {
            return Err(format!("{arm}: base dump identity mismatch").into());
        }
        Some(replay_base(&dumps[0], curve)?)
    } else {
        None
    };
    let targets = read_jsonl(&dir.join("targets_result.jsonl"))?;
    if targets.len() != answers.len() || summary["targets"].as_u64() != Some(answers.len() as u64) {
        return Err(format!("{arm}: target count mismatch").into());
    }
    let columns = summary["orbit_columns"].as_u64().ok_or("missing columns")?;
    if summary["rank"].as_u64() != Some(columns)
        || summary["targets_solved"].as_u64() != Some(answers.len() as u64)
        || summary["targets_failed"].as_u64() != Some(0)
    {
        return Err(format!("{arm}: rank or solve failure").into());
    }
    let mut probe_sum: u128 = 0;
    let mut online_sum = 0.0f64;
    let mut max_phase_residual = 0.0f64;
    for (index, (target, answer)) in targets.iter().zip(answers).enumerate() {
        let expected_point = coordinates(&answer["point"])?;
        if coordinates(&target["target"])? != expected_point
            || coordinates(&target["published_q"])? != expected_point
            || target["fixture_index"].as_u64() != Some(index as u64)
        {
            return Err(format!("{arm}: public target mismatch at row {index}").into());
        }
        let recovered = target["recovered_scalar"]
            .as_u64()
            .ok_or("missing recovered scalar")?;
        if recovered.to_string() != answer["scalar"].as_str().ok_or("missing sidecar scalar")?
            || target["group_verified"] != true
        {
            return Err(format!("{arm}: scalar or producer replay mismatch at row {index}").into());
        }
        let q = point(&answer["point"])?;
        if curve.mul(curve.generator(), &BigUint::from(recovered)) != q {
            return Err(format!("{arm}: independent generic replay failed at row {index}").into());
        }
        if let Some(base_points) = &base_points {
            let indices = target["point_indices"]
                .as_array()
                .ok_or("missing relation point indices")?;
            if indices.len() != 4 {
                return Err(
                    format!("{arm}: relation does not have four points at row {index}").into(),
                );
            }
            let mut sum = BinaryPoint::Infinity;
            for entry in indices {
                let base_point = base_points
                    .get(entry.as_u64().ok_or("nonintegral relation index")? as usize)
                    .ok_or("relation index outside archived base")?;
                if !curve.curve.is_on_curve(base_point) {
                    return Err(format!("{arm}: relation point is off curve at row {index}").into());
                }
                sum = curve.add(&sum, base_point);
            }
            if sum != q {
                return Err(format!("{arm}: four-point relation failed at row {index}").into());
            }
        }
        probe_sum += u128::from(target["probes"].as_u64().ok_or("missing probes")?);
        let online = target["target_ms"].as_f64().ok_or("missing online time")?;
        let phase_sum = target["target_phase_sum_ms"]
            .as_f64()
            .ok_or("missing phase sum")?;
        let residual = (online - phase_sum).abs();
        if !online.is_finite() || !phase_sum.is_finite() || residual > 1.0 {
            return Err(format!("{arm}: invalid online phase timing at row {index}").into());
        }
        max_phase_residual = max_phase_residual.max(residual);
        online_sum += online;
    }
    let n = targets.len() as f64;
    Ok(json!({
        "arm": arm,
        "status": "PASS_INDEPENDENT_REPLAY",
        "targets": targets.len(),
        "factor_base_points": summary["factor_base_points"],
        "orbit_columns": columns,
        "base_hash": summary["base_hash"],
        "regular_states": summary["regular_states"],
        "root_table_entries": summary["root_table_entries"],
        "root_table_allocated_bytes": summary["root_table_allocated_bytes"],
        "root_table_dense_enabled": summary["root_table_dense_enabled"],
        "root_table_load75_enabled": summary["root_table_load75_enabled"],
        "mean_target_probes": probe_sum as f64 / n,
        "mean_online_ms_exploratory": online_sum / n,
        "max_phase_residual_ms": max_phase_residual,
        "four_point_group_replay": if class == "panel" { "PASS" } else { "NOT_ARCHIVED_FOR_PRIMARY" },
        "complete_base_subgroup_orbit_label_replay": if class == "panel" { "PASS" } else { "NOT_ARCHIVED_FOR_PRIMARY" },
        "setup_ms": summary["timing_ms"]["setup_total"],
        "rank_ms": summary["timing_ms"]["rank_stage"],
        "linear_algebra_ms": summary["timing_ms"]["linear_algebra"],
        "process_total_ms": summary["timing_ms"]["process_total"],
        "peak_rss_bytes": summary["peak_rss_bytes"],
    }))
}

fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut x = *state;
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

fn paired_probe_interval(
    study: &Path,
    left: &str,
    right: &str,
) -> Result<Value, Box<dyn std::error::Error>> {
    let rows = |arm: &str| {
        read_jsonl(
            &study
                .join("runs")
                .join("panel")
                .join(arm)
                .join("targets_result.jsonl"),
        )
    };
    let left_rows = rows(left)?;
    let right_rows = rows(right)?;
    if left_rows.len() != right_rows.len() || left_rows.is_empty() {
        return Err(format!("{left}/{right}: unpaired target rows").into());
    }
    let differences: Vec<i64> = left_rows
        .iter()
        .zip(&right_rows)
        .map(|(a, b)| {
            Ok(
                i64::try_from(a["probes"].as_u64().ok_or("left probes missing")?)?
                    - i64::try_from(b["probes"].as_u64().ok_or("right probes missing")?)?,
            )
        })
        .collect::<Result<_, Box<dyn std::error::Error>>>()?;
    let n = differences.len();
    let mean = differences.iter().map(|&x| x as f64).sum::<f64>() / n as f64;
    let mut state = 20261013256u64;
    let mut resamples = Vec::with_capacity(10_000);
    for _ in 0..10_000 {
        let mut sum: i128 = 0;
        for _ in 0..n {
            let index = (splitmix64(&mut state) % n as u64) as usize;
            sum += i128::from(differences[index]);
        }
        resamples.push(sum as f64 / n as f64);
    }
    resamples.sort_by(f64::total_cmp);
    Ok(json!({
        "left": left,
        "right": right,
        "status": "PASS_PAIRED_PROBE_BOOTSTRAP",
        "direction": "left_minus_right_probes; positive means right has fewer",
        "targets": n,
        "mean_difference_probes": mean,
        "paired_bootstrap_95_percentile_probes": [resamples[250], resamples[9750]],
        "bootstrap_seed": 20261013256u64,
        "bootstrap_resamples": 10_000,
        "bootstrap_index_rule": "SplitMix64 modulo target count, paired resampling with replacement; sorted sample positions 250 and 9750",
    }))
}

fn rho_comparison(
    study: &Path,
    primary_answer: &Value,
    arms: &[Value],
) -> Result<Value, Box<dyn std::error::Error>> {
    let receipt: Value = serde_json::from_slice(&fs::read(
        study.join("../n83_dense_root_20261013/runs/rho_primary/receipt.json"),
    )?)?;
    if receipt["status"] != "PASS_VERIFIED_PUBLIC_TARGET"
        || receipt["generic_scalar_replay_verified"] != true
        || receipt["recovered_scalar"] != primary_answer["scalar"]
        || coordinates(&receipt["public_target_q"])? != coordinates(&primary_answer["point"])?
        || receipt["public_input_sha256"]
            != "21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8"
    {
        return Err("same-point rho receipt did not independently validate".into());
    }
    let rho_ms = receipt["online_ms"]
        .as_f64()
        .ok_or("rho online time missing")?;
    let comparisons: Vec<Value> = arms
        .iter()
        .filter(|arm| arm["status"] == "PASS_INDEPENDENT_REPLAY")
        .map(|arm| {
            let ic_ms = arm["mean_online_ms_exploratory"].as_f64().unwrap();
            json!({
                "arm": arm["arm"],
                "rho_online_ms": rho_ms,
                "ic_online_ms": ic_ms,
                "rho_over_ic_online_ratio_exploratory": rho_ms / ic_ms,
            })
        })
        .collect();
    Ok(json!({
        "status": "PASS_SAME_PUBLIC_POINT",
        "rho_run_id": receipt["run_id"],
        "rho_online_ms": rho_ms,
        "online_interval": receipt["online_interval"],
        "controlled_speedup": null,
        "comparisons": comparisons,
    }))
}

fn prior_k1798_parity(study: &Path, class: &str) -> Result<Value, Box<dyn std::error::Error>> {
    let old_class = if class == "primary_v2" {
        "primary"
    } else {
        class
    };
    let old_dir = study
        .join("../n83_dense_root_20261013/runs")
        .join(old_class)
        .join("k1798_dense");
    let new_dir = study.join("runs").join(class).join("k1798_dense");
    let old_summary = read_jsonl(&old_dir.join("stdout.jsonl"))?;
    let new_summary = read_jsonl(&new_dir.join("stdout.jsonl"))?;
    if old_summary.len() != 1 || new_summary.len() != 1 {
        return Err("K1798 parity summaries missing".into());
    }
    for field in [
        "base_hash",
        "factor_base_points",
        "orbit_columns",
        "regular_states",
        "root_table_entries",
        "root_table_allocated_bytes",
    ] {
        if old_summary[0][field] != new_summary[0][field] {
            return Err(format!("K1798 summary parity mismatch: {field}").into());
        }
    }
    let old_rows = read_jsonl(&old_dir.join("targets_result.jsonl"))?;
    let new_rows = read_jsonl(&new_dir.join("targets_result.jsonl"))?;
    if old_rows.len() != new_rows.len() {
        return Err("K1798 parity target count mismatch".into());
    }
    for (i, (old, new)) in old_rows.iter().zip(&new_rows).enumerate() {
        for field in [
            "target",
            "x_codes",
            "point_indices",
            "probes",
            "recovered_scalar",
        ] {
            if old[field] != new[field] {
                return Err(format!("K1798 first-witness parity mismatch at {i}/{field}").into());
            }
        }
    }
    Ok(json!({"status":"PASS_FIRST_WITNESS_PARITY", "targets":new_rows.len()}))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 || !["primary_v2", "panel"].contains(&args[2].as_str()) {
        return Err("usage: n83-load75-analyzer <study-dir> <primary_v2|panel>".into());
    }
    let study = Path::new(&args[1]);
    let answers = read_jsonl(&study.join("answers.jsonl"))?;
    let selected_answers = if args[2] == "primary_v2" {
        &answers[..1]
    } else {
        &answers[..]
    };
    let curve = KoblitzCurve::new(1, 83).ok_or("N83 curve unavailable")?;
    let arms: Vec<Value> = ARMS
        .iter()
        .map(|arm| verify_arm(study, &args[2], arm, selected_answers, &curve))
        .collect::<Result<_, _>>()?;
    let comparisons = if args[2] == "panel" {
        [
            ("k1798_dense", "k2000_load75"),
            ("k2000_load75", "k2202_load75"),
            ("k1798_dense", "k2202_load75"),
        ]
        .into_iter()
        .map(|(left, right)| {
            let complete = [left, right].iter().all(|name| {
                arms.iter()
                    .any(|arm| arm["arm"] == *name && arm["status"] == "PASS_INDEPENDENT_REPLAY")
            });
            if complete {
                paired_probe_interval(study, left, right)
            } else {
                Ok(json!({"left":left, "right":right, "status":"CENSORED_BY_ARM_EXIT"}))
            }
        })
        .collect::<Result<Vec<_>, _>>()?
    } else {
        Vec::new()
    };
    let rho = if args[2] == "primary_v2" {
        rho_comparison(study, &answers[0], &arms)?
    } else {
        Value::Null
    };
    let k1798_prior_parity = if arms[0]["status"] == "PASS_INDEPENDENT_REPLAY" {
        prior_k1798_parity(study, &args[2])?
    } else {
        json!({"status":"CENSORED_BY_ARM_EXIT"})
    };
    let selected_by_mean_probes = if args[2] == "panel" {
        arms.iter()
            .filter(|arm| arm["status"] == "PASS_INDEPENDENT_REPLAY")
            .min_by(|a, b| {
                a["mean_target_probes"]
                    .as_f64()
                    .unwrap()
                    .total_cmp(&b["mean_target_probes"].as_f64().unwrap())
            })
            .map(|arm| arm["arm"].clone())
            .unwrap_or(Value::Null)
    } else {
        Value::Null
    };
    println!(
        "{}",
        json!({"schema":"n83.load75-analysis/v1", "class":args[2], "arms":arms, "probe_comparisons":comparisons, "rho_same_point":rho, "k1798_prior_parity":k1798_prior_parity, "selected_by_mean_probes":selected_by_mean_probes})
    );
    Ok(())
}
