//! Post-run descriptive supplement. Frozen campaign statistics remain untouched.
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::os::unix::fs::PermissionsExt;
use std::path::Path;
use std::process::{Command, Stdio};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const METRICS: [(&str, &str); 6] = [
    ("elapsed_process_wall_s", "seconds"),
    ("cpu_total_s", "seconds"),
    ("peak_rss_bytes", "bytes"),
    ("producer_preparation_ms", "milliseconds"),
    ("producer_online_ms", "milliseconds"),
    ("walk_steps", "uncalibrated rho walk steps"),
];
fn hash(path: &Path) -> Result<String> {
    let output = Command::new("/usr/bin/sha256sum")
        .arg("--")
        .arg(path)
        .output()?;
    if !output.status.success() {
        return Err("SHA256 command failed".into());
    }
    Ok(String::from_utf8(output.stdout)?
        .split_whitespace()
        .next()
        .ok_or("SHA256 missing")?
        .into())
}
fn hash_bytes(bytes: &[u8]) -> Result<String> {
    let mut child = Command::new("/usr/bin/sha256sum")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;
    child
        .stdin
        .take()
        .ok_or("SHA256 stdin missing")?
        .write_all(bytes)?;
    let output = child.wait_with_output()?;
    if !output.status.success() {
        return Err("source SHA256 failed".into());
    }
    Ok(String::from_utf8(output.stdout)?
        .split_whitespace()
        .next()
        .ok_or("SHA256 missing")?
        .into())
}
fn value(row: &Value, key: &str) -> Result<f64> {
    let raw = if key == "walk_steps" {
        &row["trajectory"]["outcome"]["walk_steps"]
    } else {
        &row[key]
    };
    raw.as_f64()
        .filter(|v| v.is_finite() && *v > 0.)
        .ok_or_else(|| format!("required metric {key} is missing/nonfinite/nonpositive").into())
}
fn summary(values: &[f64]) -> Value {
    assert!(!values.is_empty());
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let n = sorted.len();
    let median = if n % 2 == 1 {
        sorted[n / 2]
    } else {
        sorted[n / 2 - 1] / 2. + sorted[n / 2] / 2.
    };
    json!({"n":n,"arithmetic_mean":values.iter().sum::<f64>()/n as f64,"median":median,"minimum":sorted[0],"maximum":sorted[n-1]})
}
fn summarize_metrics(rows: &[&Value]) -> Result<Value> {
    let mut result = json!({});
    for (key, unit) in METRICS {
        let values = rows
            .iter()
            .map(|r| value(r, key))
            .collect::<Result<Vec<_>>>()?;
        let mut metric = summary(&values);
        metric["unit"] = json!(unit);
        metric["raw_values_in_run_order"] = json!(values);
        result[key] = metric;
    }
    Ok(result)
}
fn analyze(rows: &[Value]) -> Result<Value> {
    if rows.len() != 30 {
        return Err("all thirty verified held-out rows required".into());
    }
    let mut targets = BTreeMap::<String, BTreeMap<String, BTreeMap<u64, &Value>>>::new();
    let mut runs = BTreeSet::new();
    let mut policies = BTreeMap::<String, Value>::new();
    for row in rows {
        if row["schema"] != "n83.speed-observation/v1"
            || row["status"] != "COMPLETED_VERIFIED"
            || row["verified"] != true
            || row["wrong_answer"] != false
            || row["parent_known_answer_scalar_match"] != true
        {
            return Err(
                "unverified or invalid row rejects descriptive complete-panel summary".into(),
            );
        }
        let run = row["run_id"].as_str().ok_or("run id missing")?;
        if !runs.insert(run) {
            return Err("duplicate run id".into());
        }
        let target = row["target_id"].as_str().ok_or("target missing")?;
        let target_number = target
            .strip_prefix("holdout_")
            .and_then(|s| s.parse::<u64>().ok())
            .filter(|q| (1..=5).contains(q))
            .ok_or("not one of frozen five holdout target IDs")?;
        let arm = row["arm"]
            .as_str()
            .filter(|a| ["current", "selected"].contains(a))
            .ok_or("unknown arm")?;
        let seed = row["seed"].as_u64().ok_or("numeric seed missing")?;
        let first_seed = 202610126241 + (target_number - 1) * 3;
        if !(first_seed..first_seed + 3).contains(&seed) {
            return Err("seed differs from frozen target assignment".into());
        }
        if let Some(previous) = policies.insert(arm.into(), row["config"].clone()) {
            if previous != row["config"] {
                return Err("arm policy differs within holdout".into());
            }
        }
        for (key, _) in METRICS {
            value(row, key)?;
        }
        if targets
            .entry(target.into())
            .or_default()
            .entry(arm.into())
            .or_default()
            .insert(seed, row)
            .is_some()
        {
            return Err("duplicate target/arm/seed".into());
        }
    }
    if targets.len() != 5 {
        return Err("five targets required".into());
    }
    let mut arm_rows = BTreeMap::<String, Vec<&Value>>::new();
    for row in rows {
        arm_rows
            .entry(row["arm"].as_str().unwrap().into())
            .or_default()
            .push(row);
    }
    let mut arms = json!({});
    for (arm, values) in &arm_rows {
        if values.len() != 15 {
            return Err("fifteen draws per arm required".into());
        }
        arms[arm] = json!({"config":policies[arm],"metrics":summarize_metrics(values)?});
    }
    let mut per_target = Vec::new();
    let mut ratios = Vec::new();
    let mut ratio_values = Vec::new();
    for (target, block) in &targets {
        let current = block.get("current").ok_or("current arm missing")?;
        let selected = block.get("selected").ok_or("selected arm missing")?;
        if current.len() != 3 || selected.len() != 3 || !current.keys().eq(selected.keys()) {
            return Err("three matching distinct seeds per target required".into());
        }
        let mut summaries = json!({});
        for (arm, draws) in block {
            summaries[arm] = summarize_metrics(&draws.values().copied().collect::<Vec<_>>())?;
        }
        per_target.push(json!({"target_id":target,"arms":summaries}));
        for (seed, c) in current {
            let s = selected[seed];
            if c["recovered_scalar"] != s["recovered_scalar"] {
                return Err("paired recovered scalars differ".into());
            }
            let cw = value(c, "elapsed_process_wall_s")?;
            let sw = value(s, "elapsed_process_wall_s")?;
            let ratio = cw / sw;
            if !ratio.is_finite() {
                return Err("paired ratio nonfinite".into());
            }
            ratio_values.push(ratio);
            ratios.push(json!({"target_id":target,"seed":seed,"current_process_wall_s":cw,"selected_process_wall_s":sw,"current_over_selected_process_wall_ratio":ratio}));
        }
    }
    Ok(
        json!({"schema":"n83.speed-descriptive-summary/v1","status":"ALL30_VERIFIED_COMPLETE_PANEL","observed_targets":5,"draws_per_target_arm":3,"observed_rows":30,
        "scope":"Post-run descriptive supplement for the frozen thirty-row holdout. Fifteen paired launches are descriptive repeated draws, not fifteen independent target clusters. Frozen primary statistics and target-cluster bootstrap remain unchanged. No new experiment or confidence interval.",
        "arms":arms,"per_target":per_target,"paired_current_over_selected_process_wall":{"descriptive_summary":summary(&ratio_values),"all15_pairs":ratios,"confidence_interval":Value::Null,"estimand":"Distribution of fifteen matched launch ratios; descriptive only, distinct from frozen five-target ratio estimand"},
        "missing_value_policy":"All six metrics must be finite and positive on all thirty verified rows. No missing values, failed draws, zero imputation, omission or replacement. Descriptions include every held-out draw."}),
    )
}
fn main_run() -> Result<()> {
    let args = std::env::args().collect::<Vec<_>>();
    if args.len() != 3 {
        return Err("usage: n83-speed-descriptive INPUT_JSONL NEW_OUTPUT_JSON".into());
    }
    let input = Path::new(&args[1]);
    let output = Path::new(&args[2]);
    let rows = BufReader::new(File::open(input)?)
        .lines()
        .filter_map(|l| match l {
            Ok(l) if l.trim().is_empty() => None,
            other => Some(other),
        })
        .map(|l| Ok(serde_json::from_str(&l?)?))
        .collect::<Result<Vec<Value>>>()?;
    let mut result = analyze(&rows)?;
    let executable = std::env::current_exe()?;
    result["provenance"] = json!({"command":args,"input_file":input,"input_sha256":hash(input)?,"executable":executable,"executable_sha256":hash(&executable)?,"compiled_source_sha256":hash_bytes(include_bytes!("main.rs"))?,"compiled_cargo_manifest_sha256":hash_bytes(include_bytes!("../Cargo.toml"))?,"compiled_cargo_lock_sha256":hash_bytes(include_bytes!("../Cargo.lock"))?});
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)?;
    serde_json::to_writer_pretty(&mut file, &result)?;
    file.write_all(b"\n")?;
    file.sync_all()?;
    fs::set_permissions(output, fs::Permissions::from_mode(0o444))?;
    println!(
        "{}",
        json!({"status":result["status"],"output":output,"output_sha256":hash(output)?})
    );
    Ok(())
}
fn main() {
    if let Err(e) = main_run() {
        eprintln!("descriptive analysis rejected: {e}");
        std::process::exit(2);
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn odd_median_ignores_input_order() {
        let s = summary(&[9., 1., 3.]);
        assert_eq!(s["median"], 3.);
        assert_eq!(s["minimum"], 1.);
        assert_eq!(s["maximum"], 9.);
    }
    #[test]
    fn even_median_averages_middle_pair() {
        let s = summary(&[8., 2., 4., 6.]);
        assert_eq!(s["median"], 5.);
        assert_eq!(s["arithmetic_mean"], 5.);
    }
    #[test]
    fn singleton_summary_retains_value() {
        let s = summary(&[7.]);
        assert_eq!(s["median"], 7.);
        assert_eq!(s["n"], 1);
    }
    #[test]
    fn incomplete_panel_rejected() {
        assert!(analyze(&[]).is_err());
    }
    #[test]
    fn absent_metric_rejected() {
        assert!(value(&json!({}), "cpu_total_s").is_err());
        assert!(value(&json!({"cpu_total_s":0}), "cpu_total_s").is_err());
    }
}
