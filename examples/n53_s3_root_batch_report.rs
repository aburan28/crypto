//! Native, read-only aggregation of the verified N53 root-batch pair records.
//! Usage: n53_s3_root_batch_report <ANALYSIS.json> <RESULT_DATA.json>

use serde_json::{json, Value};
use std::fs;
use std::path::Path;

fn number(row: &Value, key: &str) -> Result<f64, String> {
    row[key].as_f64().ok_or_else(|| format!("missing {key}"))
}

fn count(row: &Value, key: &str) -> Result<u64, String> {
    row[key].as_u64().ok_or_else(|| format!("missing {key}"))
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_by(|a, b| a.total_cmp(b));
    values[values.len() / 2]
}

fn target_block_interval(values: &[f64]) -> [f64; 2] {
    let mut state = 0x53_BA_7C_12_2026u64;
    let mut draws = Vec::with_capacity(10_000);
    for _ in 0..10_000 {
        let mut sample = Vec::with_capacity(9);
        for _ in 0..3 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let target = (state % 3) as usize;
            sample.extend_from_slice(&values[target * 3..target * 3 + 3]);
        }
        draws.push(median(sample));
    }
    draws.sort_by(|a, b| a.total_cmp(b));
    [draws[250], draws[9750]]
}

fn summarize(rows: &[Value], key: &str) -> Result<Value, String> {
    let values = rows
        .iter()
        .map(|row| number(row, key))
        .collect::<Result<Vec<_>, _>>()?;
    if values.len() != 9 || values.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(format!("{key}: expected nine finite nonnegative values"));
    }
    Ok(json!({
        "median":median(values.clone()),
        "target_block_95_interval":target_block_interval(&values),
        "values":values
    }))
}

fn run(input: &Path, output: &Path) -> Result<(), String> {
    if output.exists() {
        return Err(format!("refusing to overwrite {}", output.display()));
    }
    let analysis: Value = serde_json::from_slice(&fs::read(input).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    if analysis["schema"] != "n53-s3-root-batch-analysis-v1"
        || analysis["verified_pairs"] != 9
        || analysis["exact_gate"] != true
        || analysis["inversion_gate"] != true
        || !analysis["missing_or_invalid"]
            .as_array()
            .is_some_and(Vec::is_empty)
    {
        return Err("analysis is incomplete or failed an exact gate".into());
    }
    let rows = analysis["pairs"].as_array().ok_or("pairs absent")?;
    if rows.len() != 9 {
        return Err("pair count differs".into());
    }
    let mut inversion_control = 0u64;
    let mut inversion_candidate = 0u64;
    let mut multiplication_control = 0u64;
    let mut multiplication_candidate = 0u64;
    let mut squaring_control = 0u64;
    let mut squaring_candidate = 0u64;
    for (index, row) in rows.iter().enumerate() {
        if row["status"] != "VERIFIED"
            || row["target_index"] != (index / 3) as u64
            || row["seed_index"] != (index % 3) as u64
            || row["control_replay"]["verified"] != true
            || row["candidate_replay"]["verified"] != true
        {
            return Err(format!("pair {index} order or replay differs"));
        }
        inversion_control += count(&row["control_root_ops"], "inversion_primitives")?;
        inversion_candidate += count(&row["candidate_root_ops"], "inversion_primitives")?;
        multiplication_control += count(&row["control_root_ops"], "field_multiplications")?;
        multiplication_candidate += count(&row["candidate_root_ops"], "field_multiplications")?;
        squaring_control += count(&row["control_root_ops"], "field_squarings")?;
        squaring_candidate += count(&row["candidate_root_ops"], "field_squarings")?;
    }
    let result = json!({
        "schema":"n53-s3-root-batch-result-data-v1",
        "source_analysis":input,
        "verified_pairs":9,
        "cold_candidate_over_control":summarize(rows,"candidate_over_control_cold")?,
        "online_candidate_over_control":summarize(rows,"candidate_over_control_online")?,
        "cold_control_over_rho":summarize(rows,"control_over_rho_cold")?,
        "cold_candidate_over_rho":summarize(rows,"candidate_over_rho_cold")?,
        "online_control_over_rho":summarize(rows,"control_over_rho_online")?,
        "online_candidate_over_rho":summarize(rows,"candidate_over_rho_online")?,
        "rank_root_ops_nine_cells":{
            "control_inversions":inversion_control,
            "candidate_inversions":inversion_candidate,
            "candidate_over_control_inversions":inversion_candidate as f64 / inversion_control as f64,
            "control_field_multiplications":multiplication_control,
            "candidate_field_multiplications":multiplication_candidate,
            "candidate_over_control_field_multiplications":multiplication_candidate as f64 / multiplication_control as f64,
            "control_field_squarings":squaring_control,
            "candidate_field_squarings":squaring_candidate,
            "candidate_over_control_field_squarings":squaring_candidate as f64 / squaring_control as f64
        },
        "timing_class":"shared_macOS_L0_descriptive",
        "common_operation_unit_S":Value::Null
    });
    fs::write(output, format!("{result:#}\n")).map_err(|e| e.to_string())?;
    println!(
        "{}",
        json!({
            "cold_candidate_over_control":result["cold_candidate_over_control"]["median"],
            "online_candidate_over_control":result["online_candidate_over_control"]["median"],
            "cold_candidate_over_rho":result["cold_candidate_over_rho"]["median"],
            "online_candidate_over_rho":result["online_candidate_over_rho"]["median"]
        })
    );
    Ok(())
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: n53_s3_root_batch_report <ANALYSIS.json> <RESULT_DATA.json>");
        std::process::exit(2);
    }
    if let Err(error) = run(Path::new(&args[1]), Path::new(&args[2])) {
        eprintln!("n53 S3 root-batch report: {error}");
        std::process::exit(1);
    }
}
