//! Planning sample size for the paired S3 inversion experiment.
//!
//! Each input is a JSON object with `pairs`, one entry per independent target.
//! A pair supplies either `exploratory_online_ratio`, `ratio`, or positive
//! `baseline_online_ms` and `candidate_online_ms`. Repeats of one target must
//! be combined before this program: counting them as independent would make
//! the required sample size too small.

use serde_json::{json, Value};
use std::collections::HashSet;
use std::env;
use std::fs;

const Z_TWO_SIDED_95: f64 = 1.959_963_984_540_054;
const Z_POWER_80: f64 = 0.841_621_233_572_914_3;

fn pair_ratio(pair: &Value) -> Result<f64, String> {
    let ratio = pair
        .get("ratio")
        .or_else(|| pair.get("exploratory_online_ratio"))
        .and_then(Value::as_f64)
        .or_else(|| {
            let baseline = pair.get("baseline_online_ms")?.as_f64()?;
            let candidate = pair.get("candidate_online_ms")?.as_f64()?;
            (candidate > 0.0).then_some(baseline / candidate)
        })
        .ok_or("each pair needs a ratio or both online times")?;
    if !ratio.is_finite() || ratio <= 0.0 {
        return Err("a pair ratio must be finite and positive".into());
    }
    Ok(ratio)
}

fn log_ratios(doc: &Value) -> Result<Vec<f64>, String> {
    let pairs = doc
        .get("pairs")
        .and_then(Value::as_array)
        .ok_or("input needs a pairs array")?;
    if pairs.len() < 2 {
        return Err("at least two independent target pairs are needed".into());
    }
    let mut workloads = HashSet::new();
    let mut logs = Vec::with_capacity(pairs.len());
    for pair in pairs {
        let id = pair
            .get("workload_id")
            .and_then(Value::as_str)
            .filter(|id| !id.is_empty())
            .ok_or("every pair needs a canonical, nonempty workload_id")?;
        if !workloads.insert(id) {
            return Err(format!(
                "duplicate workload_id {id}; aggregate its repeats first"
            ));
        }
        logs.push(pair_ratio(pair)?.ln());
    }
    Ok(logs)
}

fn sample_sd(values: &[f64]) -> f64 {
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance =
        values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64;
    variance.sqrt()
}

fn needed(variance: f64, minimum_ratio: f64) -> usize {
    let z = Z_TWO_SIDED_95 + Z_POWER_80;
    let raw = (z * z * variance / minimum_ratio.ln().powi(2)).ceil();
    raw.max(5.0) as usize
}

fn confirmatory_size(interaction_targets_per_curve: usize) -> usize {
    // Twenty percent pilot-variance allowance, frozen before new target seeds.
    ((interaction_targets_per_curve as f64 * 1.20).ceil() as usize).max(24)
}

fn plan(first: &Value, second: Option<&Value>, minimum_ratio: f64) -> Result<Value, String> {
    if !minimum_ratio.is_finite() || minimum_ratio <= 1.0 {
        return Err("minimum detectable ratio must exceed 1".into());
    }
    let a = log_ratios(first)?;
    let b = second.map(log_ratios).transpose()?;
    let sd_a = sample_sd(&a);
    let sd_b = b.as_ref().map(|x| sample_sd(x)).unwrap_or(sd_a);
    let interaction_targets_per_curve = needed(sd_a * sd_a + sd_b * sd_b, minimum_ratio);
    let mean = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
    Ok(json!({
        "kind": "paired_s3_sample_size_planning_only",
        "alpha_two_sided": 0.05,
        "power": 0.8,
        "minimum_detectable_ratio": minimum_ratio,
        "unit_of_replication": "independent one-target workload; repeat timings for a target must be aggregated first",
        "curve_a": {
            "pilot_targets": a.len(),
            "geometric_mean_baseline_over_candidate": mean(&a).exp(),
            "sd_log_ratio": sd_a,
            "targets_for_within_curve_effect": needed(sd_a * sd_a, minimum_ratio)
        },
        "curve_b": b.as_ref().map(|x| json!({
            "pilot_targets": x.len(),
            "geometric_mean_baseline_over_candidate": mean(x).exp(),
            "sd_log_ratio": sd_b,
            "targets_for_within_curve_effect": needed(sd_b * sd_b, minimum_ratio)
        })),
        "interaction_targets_per_curve": interaction_targets_per_curve,
        "confirmatory_targets_per_curve": confirmatory_size(interaction_targets_per_curve),
        "confirmatory_rule": "max(24, ceil(1.20 * interaction_targets_per_curve)); fresh targets, pilot excluded",
        "second_curve_variance": if b.is_some() { "measured_pilot" } else { "assumed_equal_to_curve_a" },
        "method": "normal-approximation two-sided test of the difference of mean paired log ratios; 80% power; pilot data excluded from the eventual fixed-size confirmatory panel",
        "caution": "A planning calculation is not a confidence interval or a speedup claim. Recalculate on clean, isolated pilot data and freeze the confirmatory size before generating its targets."
    }))
}

fn run() -> Result<(), String> {
    let args: Vec<String> = env::args().collect();
    if !(2..=4).contains(&args.len()) {
        return Err("usage: s3_pair_power CURVE_A.json [CURVE_B.json] [MIN_RATIO=1.10]".into());
    }
    let a: Value = serde_json::from_slice(&fs::read(&args[1]).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    if args.len() == 4 && args[2].parse::<f64>().is_ok() {
        return Err("with three arguments, the second must be CURVE_B.json".into());
    }
    let second_path = args.get(2).filter(|s| s.parse::<f64>().is_err());
    let b: Option<Value> = second_path
        .map(|p| fs::read(p).map_err(|e| e.to_string()))
        .transpose()?
        .map(|data| serde_json::from_slice(&data).map_err(|e| e.to_string()))
        .transpose()?;
    let ratio_arg = if second_path.is_some() {
        args.get(3)
    } else {
        args.get(2)
    };
    let minimum_ratio = ratio_arg
        .map(|s| s.parse::<f64>().map_err(|e| e.to_string()))
        .transpose()?
        .unwrap_or(1.10);
    let result = plan(&a, b.as_ref(), minimum_ratio)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&result).map_err(|e| e.to_string())?
    );
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("s3_pair_power: {error}");
        std::process::exit(2);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interaction_uses_two_curve_variances() {
        let a = json!({"pairs": [{"workload_id": "a", "ratio": 0.8},
                                   {"workload_id": "b", "ratio": 1.2},
                                   {"workload_id": "c", "ratio": 1.0}]});
        let b = json!({"pairs": [{"workload_id": "d", "ratio": 0.9},
                                   {"workload_id": "e", "ratio": 1.1},
                                   {"workload_id": "f", "ratio": 1.0}]});
        let one = plan(&a, None, 1.10).unwrap();
        let two = plan(&a, Some(&b), 1.10).unwrap();
        assert_ne!(
            one["interaction_targets_per_curve"],
            two["interaction_targets_per_curve"]
        );
        assert_eq!(one["second_curve_variance"], "assumed_equal_to_curve_a");
        assert_eq!(two["second_curve_variance"], "measured_pilot");
    }

    #[test]
    fn repeated_target_is_not_a_new_sample() {
        let repeated = json!({"pairs": [{"workload_id": "x", "ratio": 1.1},
                                          {"workload_id": "x", "ratio": 1.2}]});
        assert!(log_ratios(&repeated)
            .unwrap_err()
            .contains("duplicate workload_id"));
        let missing = json!({"pairs": [{"ratio": 1.1}, {"ratio": 1.2}]});
        assert!(log_ratios(&missing).unwrap_err().contains("workload_id"));
    }

    #[test]
    fn invalid_timing_is_rejected() {
        assert!(pair_ratio(&json!({"baseline_online_ms": 4.0,
                                   "candidate_online_ms": 0.0}))
        .is_err());
        assert!(plan(
            &json!({"pairs": [{"ratio": 1.0}, {"ratio": 1.0}]}),
            None,
            1.0
        )
        .is_err());
    }

    #[test]
    fn confirmation_applies_prespecified_inflation_and_floor() {
        assert_eq!(confirmatory_size(5), 24);
        assert_eq!(confirmatory_size(76), 92);
    }
}
