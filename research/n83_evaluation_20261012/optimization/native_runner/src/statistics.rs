//! Frozen, native summaries for the degree-83 holdout panel.
//!
//! The resampling unit is a target, never an individual rho seed. Three rho
//! draws estimate each target's arithmetic cost. The frozen primary estimand
//! is the geometric mean of the five default/selected target cost ratios.
//! Capped arithmetic expected cost is a separate secondary description.

use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};

const EXPECTED_TARGETS: usize = 5;
const EXPECTED_RHO_DRAWS: usize = 3;
const FAILURE_CHARGE_S: f64 = 60.0;
const BOOTSTRAP_REPLICATES: usize = 20_000;
const BOOTSTRAP_SEED: u64 = 202_610_128_399;
const VERIFIED: &str = "COMPLETED_VERIFIED";

#[derive(Default)]
struct TargetRows<'a> {
    default: Vec<&'a Value>,
    selected: Vec<&'a Value>,
    ic: Vec<&'a Value>,
}

struct Arm {
    output: Value,
    seed_ids: BTreeSet<u64>,
    actual_process_mean: Option<f64>,
    process_mean: Option<f64>,
    online_mean: Option<f64>,
}

/// Summarize the frozen five-target, three-seed holdout without omitting failed
/// draws. All input rows are returned unchanged under `raw_rows`.
pub fn summarize(rows: &[Value], ic_rows: &[Value]) -> Value {
    let mut targets: BTreeMap<String, TargetRows<'_>> = BTreeMap::new();
    let mut issues = Vec::<String>::new();
    let mut status_counts = BTreeMap::<String, usize>::new();
    for (index, row) in rows.iter().enumerate() {
        let status = row["status"].as_str().unwrap_or("MISSING_STATUS");
        *status_counts.entry(status.to_owned()).or_default() += 1;
        let Some(target) = target_id(row) else {
            issues.push(format!("rho row {index}: missing or invalid target_id"));
            continue;
        };
        let target_rows = targets.entry(target.to_owned()).or_default();
        match row["arm"].as_str() {
            Some("default") => target_rows.default.push(row),
            Some("selected") => target_rows.selected.push(row),
            _ => issues.push(format!("rho row {index}: unknown or missing arm")),
        }
    }
    let rho_target_ids: BTreeSet<String> = targets.keys().cloned().collect();
    for (index, row) in ic_rows.iter().enumerate() {
        let Some(target) = target_id(row) else {
            issues.push(format!("IC row {index}: missing or invalid target_id"));
            continue;
        };
        targets.entry(target.to_owned()).or_default().ic.push(row);
    }
    let mut per_target = Vec::new();
    let mut process_pairs = Vec::<(f64, f64)>::new();
    let mut charged_pairs = Vec::<(f64, f64)>::new();
    let mut online_pairs = Vec::<(f64, f64)>::new();
    let mut ic_pairs = Vec::<(f64, f64)>::new();
    let mut process_pair_ids = Vec::<String>::new();
    let mut charged_pair_ids = Vec::<String>::new();
    let mut online_pair_ids = Vec::<String>::new();
    let mut ic_pair_ids = Vec::<String>::new();
    let mut process_default = Vec::<f64>::new();
    let mut process_selected = Vec::<f64>::new();
    let mut charge_default = Vec::<f64>::new();
    let mut charge_selected = Vec::<f64>::new();
    let mut online_default = Vec::<f64>::new();
    let mut online_selected = Vec::<f64>::new();
    let mut ic_online = Vec::<f64>::new();

    for (target, target_rows) in &targets {
        let default = summarize_arm(&target_rows.default);
        let selected = summarize_arm(&target_rows.selected);
        let paired_seeds_match =
            default.seed_ids.len() == EXPECTED_RHO_DRAWS && default.seed_ids == selected.seed_ids;
        if !paired_seeds_match {
            issues.push(format!("rho row panel for {target}: default and selected must share the same three distinct numeric seeds"));
        }
        let (ic, ic_mean) = summarize_ic(&target_rows.ic);
        let process_ratio = paired_ratio(default.actual_process_mean, selected.actual_process_mean);
        let online_ratio = paired_ratio(default.online_mean, selected.online_mean);
        let ic_ratio = paired_ratio(ic_mean, selected.online_mean);
        if let (Some(n), Some(d)) = (default.actual_process_mean, selected.actual_process_mean) {
            process_pairs.push((n, d));
            process_pair_ids.push(target.clone());
            process_default.push(n);
            process_selected.push(d);
        }
        if let (Some(n), Some(d)) = (selected.process_mean, default.process_mean) {
            charged_pairs.push((n, d));
            charged_pair_ids.push(target.clone());
            charge_selected.push(n);
            charge_default.push(d);
        }
        if let (Some(n), Some(d)) = (default.online_mean, selected.online_mean) {
            online_pairs.push((n, d));
            online_pair_ids.push(target.clone());
            online_default.push(n);
            online_selected.push(d);
        }
        if let (Some(n), Some(d)) = (ic_mean, selected.online_mean) {
            ic_pairs.push((n, d));
            ic_pair_ids.push(target.clone());
            ic_online.push(n);
        }
        per_target.push(json!({
            "target_id": target,
            "default": default.output,
            "selected": selected.output,
            "paired_rho_seed_sets_match": paired_seeds_match,
            "prepared_ic": ic,
            "default_over_selected_process_wall_ratio": process_ratio,
            "default_over_selected_online_ratio": online_ratio,
            "prepared_ic_over_selected_rho_online_ratio": ic_ratio,
            "selected_rho_over_prepared_ic_online_ratio": paired_ratio(selected.online_mean, ic_mean),
            "selected_over_default_capped_expected_cost_ratio": paired_ratio(selected.process_mean, default.process_mean)
        }));
    }

    let wrong_answer_count = rows
        .iter()
        .filter(|row| row["status"] == "WRONG_ANSWER")
        .count();
    let rho_panel_valid = rho_target_ids.len() == EXPECTED_TARGETS
        && wrong_answer_count == 0
        && issues.iter().all(|issue| !issue.starts_with("rho row"));
    let ic_panel_valid = rho_panel_valid
        && process_pairs.len() == EXPECTED_TARGETS
        && targets.len() == EXPECTED_TARGETS
        && ic_rows.len() == EXPECTED_TARGETS
        && issues.iter().all(|issue| !issue.starts_with("IC row"));
    if rho_target_ids.len() != EXPECTED_TARGETS {
        issues.push(format!(
            "expected {EXPECTED_TARGETS} rho targets; observed {}",
            rho_target_ids.len()
        ));
    }
    if ic_rows.len() != EXPECTED_TARGETS {
        issues.push(format!(
            "expected {EXPECTED_TARGETS} prepared IC rows; observed {}",
            ic_rows.len()
        ));
    }
    if wrong_answer_count > 0 {
        issues.push("WRONG_ANSWER rejects holdout soundness; every aggregate and confidence interval is withheld.".to_owned());
    }

    json!({
        "schema": "n83_optimization_statistics_v1",
        "statistical_scope": "Five disjoint known-answer public synthetic holdout targets on one host; three rho seed draws per target and arm. Targets are sampled without replacement conditional on exclusions of shared signed-Frobenius orbits with the original target, retained base, and other generated targets; they are not unconditional independent draws. Three training Qs used for fitting are outside this holdout. No deployed-curve, asymptotic, or cross-host inference.",
        "expected_targets": EXPECTED_TARGETS,
        "observed_rho_targets": rho_target_ids.len(),
        "observed_all_targets": targets.len(),
        "expected_rho_draws_per_target_arm": EXPECTED_RHO_DRAWS,
        "input_identity_contract": "Caller must supply the already-frozen holdout manifest only, with its declared target IDs and seed assignments. This module verifies five distinct nonempty target IDs, three distinct numeric seeds per arm, and matching arm seed sets within each target; it does not reconstruct the target manifest or prove that rows came from held-out rather than training Qs. Caller must provide IC rows only after a verified uncensored batch with successful process disposition, rank600/count5/basehash semantic validation, and independent scalar/decomposition replay; group_verified must not remain true for a rejected batch.",
        "observed_rho_rows": rows.len(),
        "observed_ic_rows": ic_rows.len(),
        "rho_status_counts": status_counts,
        "wrong_answer_count": wrong_answer_count,
        "holdout_soundness_rejected": wrong_answer_count > 0,
        "failure_charge_wall_s": FAILURE_CHARGE_S,
        "failure_accounting": "Primary ratios require all verified uncensored draws and use actual measured process wall. For the separate capped expected-cost secondary, a verified draw is charged min(actual process wall, 60 seconds), and every TIMEOUT, PROCESS_FAILED, or SEMANTIC_VALIDATION_FAILED draw is charged 60 seconds. Missing draws are not charged or imputed. WRONG_ANSWER rejects holdout soundness.",
        "estimand": "Frozen primary: geometric mean of five default/selected per-target ratios, after averaging the three actual full-process rho wall times within each target and arm. All draws must complete and verify. Arithmetic capped expected cost is explicitly secondary.",
        "primary_default_over_selected_process_wall": comparison(
            &process_pairs, &process_pair_ids, rho_panel_valid, BOOTSTRAP_SEED,
            "seconds", "verified uncensored actual full-process wall", "default", "selected", RatioKind::GeometricTargetRatios
        ),
        "secondary_default_over_selected_online": comparison(
            &online_pairs, &online_pair_ids, rho_panel_valid, BOOTSTRAP_SEED,
            "milliseconds", "verified online only", "default", "selected", RatioKind::GeometricTargetRatios
        ),
        "secondary_prepared_ic_over_selected_rho_online": comparison(
            &ic_pairs, &ic_pair_ids, ic_panel_valid, BOOTSTRAP_SEED + 1,
            "milliseconds", "verified online only, prepared shared IC batch; unequal process caps: IC 1800 seconds versus rho 60 seconds", "prepared_ic", "selected_rho", RatioKind::GeometricTargetRatios
        ),
        "secondary_capped_expected_cost_selected_over_default_process_wall": comparison(
            &charged_pairs, &charged_pair_ids, rho_panel_valid, BOOTSTRAP_SEED,
            "seconds", "capped expected-cost description retaining all attempted draws, including failures; no unqualified full-workload winner", "selected", "default", RatioKind::ArithmeticMeans
        ),
        "prepared_ic_accounting": "One verified prepared IC measurement per holdout target from one shared-preparation batch. The ratio excludes reusable IC preparation and is not a cold-start or end-to-end IC comparison. The target bootstrap does not measure IC setup variation or IC stochastic run-to-run variation.",
        "per_target": per_target,
        "complete_pair_counts": {"verified_process_wall": process_pairs.len(), "capped_expected_cost": charged_pairs.len(), "rho_online": online_pairs.len(), "rho_vs_ic_online": ic_pairs.len()},
        "descriptive_complete_target_means": {
            "default_actual_process_wall_s": descriptive(&process_default),
            "selected_actual_process_wall_s": descriptive(&process_selected),
            "default_charged_process_wall_s": descriptive(&charge_default),
            "selected_charged_process_wall_s": descriptive(&charge_selected),
            "default_verified_online_ms": descriptive(&online_default),
            "selected_verified_online_ms": descriptive(&online_selected),
            "prepared_ic_verified_online_ms": descriptive(&ic_online)
        },
        "input_issues": issues,
        "raw_rows": {"rho": rows, "prepared_ic": ic_rows},
        "missing_value_policy": "Null stays null; no zero imputation. Failed or incomplete draws do not enter verified online or actual process-wall means. Aggregate primary and online ratios and confidence intervals require all five target pairs and all three verified rho draws per relevant arm. IC comparison additionally requires all held-out rho completions and all five verified IC targets."
    })
}

fn target_id(row: &Value) -> Option<&str> {
    row["target_id"].as_str().filter(|id| !id.is_empty())
}

fn nonnegative(value: &Value) -> Option<f64> {
    value.as_f64().filter(|v| v.is_finite() && *v >= 0.0)
}

fn seed_id(row: &Value) -> Option<u64> {
    row["seed"]
        .as_u64()
        .or_else(|| row["seed"].as_str()?.parse::<u64>().ok())
}

fn summarize_arm(rows: &[&Value]) -> Arm {
    let mut issues = Vec::<String>::new();
    let mut counts = BTreeMap::<String, usize>::new();
    let mut seeds = BTreeSet::<u64>::new();
    let mut charge_values = Vec::<f64>::new();
    let mut verified_online_values = Vec::<f64>::new();
    let mut observed_wall_values = Vec::<f64>::new();
    let mut verified_count = 0;
    let mut valid_seed_count = 0;
    for (index, row) in rows.iter().enumerate() {
        let status = row["status"].as_str().unwrap_or("MISSING_STATUS");
        *counts.entry(status.to_owned()).or_default() += 1;
        if let Some(seed) = seed_id(row) {
            valid_seed_count += 1;
            if !seeds.insert(seed) {
                issues.push(format!("draw {index}: duplicate seed"));
            }
        } else {
            issues.push(format!("draw {index}: missing or invalid seed"));
        }
        let wall = nonnegative(&row["elapsed_process_wall_s"]);
        if let Some(wall) = wall {
            observed_wall_values.push(wall);
        }
        let charge = nonnegative(&row["capped_charge_wall_s"]);
        match status {
            VERIFIED => {
                verified_count += 1;
                if let Some(online) = nonnegative(&row["producer_online_ms"]) {
                    verified_online_values.push(online);
                } else {
                    issues.push(format!(
                        "draw {index}: verified row has missing or invalid online time"
                    ));
                }
                match (charge, wall) {
                    (Some(charge), Some(wall)) if charge == wall.min(FAILURE_CHARGE_S) => charge_values.push(charge),
                    _ => issues.push(format!("draw {index}: verified charge must equal min(measured process wall, 60 seconds)")),
                }
            }
            "TIMEOUT" | "PROCESS_FAILED" | "SEMANTIC_VALIDATION_FAILED" => {
                if charge == Some(FAILURE_CHARGE_S) {
                    charge_values.push(FAILURE_CHARGE_S);
                } else {
                    issues.push(format!("draw {index}: failure charge must be 60 seconds"));
                }
            }
            _ => issues.push(format!("draw {index}: unknown or missing status")),
        }
    }
    let complete_draws = rows.len() == EXPECTED_RHO_DRAWS
        && valid_seed_count == EXPECTED_RHO_DRAWS
        && seeds.len() == EXPECTED_RHO_DRAWS;
    if rows.len() != EXPECTED_RHO_DRAWS {
        issues.push(format!(
            "expected {EXPECTED_RHO_DRAWS} draws; observed {}",
            rows.len()
        ));
    }
    let process_mean = if complete_draws && charge_values.len() == EXPECTED_RHO_DRAWS {
        mean(&charge_values)
    } else {
        None
    };
    let actual_process_mean = if complete_draws
        && verified_count == EXPECTED_RHO_DRAWS
        && observed_wall_values.len() == EXPECTED_RHO_DRAWS
    {
        mean(&observed_wall_values)
    } else {
        None
    };
    let online_mean = if complete_draws
        && verified_count == EXPECTED_RHO_DRAWS
        && verified_online_values.len() == EXPECTED_RHO_DRAWS
    {
        mean(&verified_online_values)
    } else {
        None
    };
    Arm {
        output: json!({
            "draw_count": rows.len(),
            "expected_draw_count": EXPECTED_RHO_DRAWS,
            "verified_draw_count": verified_count,
            "failure_draw_count": rows.len().saturating_sub(verified_count),
            "unique_seed_count": seeds.len(),
            "status_counts": counts,
            "complete_charged_process_draws": process_mean.is_some(),
            "complete_verified_actual_process_draws": actual_process_mean.is_some(),
            "complete_verified_online_draws": online_mean.is_some(),
            "mean_charged_process_wall_s": process_mean,
            "mean_verified_actual_process_wall_s": actual_process_mean,
            "mean_verified_online_ms": online_mean,
            "descriptive_observed_process_wall_s": descriptive(&observed_wall_values),
            "descriptive_valid_charge_wall_s": descriptive(&charge_values),
            "descriptive_verified_online_ms": descriptive(&verified_online_values),
            "issues": issues
        }),
        seed_ids: seeds,
        actual_process_mean,
        process_mean,
        online_mean,
    }
}

fn summarize_ic(rows: &[&Value]) -> (Value, Option<f64>) {
    let mut issues = Vec::<String>::new();
    if rows.len() != 1 {
        issues.push(format!(
            "expected one prepared IC row; observed {}",
            rows.len()
        ));
    }
    let verified = rows.len() == 1 && rows[0]["group_verified"].as_bool() == Some(true);
    let time = if verified {
        nonnegative(&rows[0]["target_ms"])
    } else {
        None
    };
    if rows.len() == 1 && !verified {
        issues.push("prepared IC row is not group-verified".to_owned());
    }
    if verified && time.is_none() {
        issues.push("verified prepared IC row has missing or invalid target_ms".to_owned());
    }
    (
        json!({
            "draw_count": rows.len(),
            "group_verified": verified,
            "target_ms": time,
            "issues": issues
        }),
        time,
    )
}

fn mean(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        None
    } else {
        let result = values.iter().sum::<f64>() / values.len() as f64;
        result.is_finite().then_some(result)
    }
}

fn descriptive(values: &[f64]) -> Value {
    json!({
        "n": values.len(),
        "arithmetic_mean": mean(values),
        "minimum": values.iter().copied().reduce(f64::min),
        "maximum": values.iter().copied().reduce(f64::max)
    })
}

fn paired_ratio(numerator: Option<f64>, denominator: Option<f64>) -> Option<f64> {
    match (numerator, denominator) {
        (Some(n), Some(d)) if d > 0.0 => {
            let ratio = n / d;
            ratio.is_finite().then_some(ratio)
        }
        _ => None,
    }
}

#[derive(Clone, Copy)]
enum RatioKind {
    GeometricTargetRatios,
    ArithmeticMeans,
}

impl RatioKind {
    fn name(self) -> &'static str {
        match self {
            Self::GeometricTargetRatios => "geometric mean of paired per-target cost ratios; each rho cost is its three-seed arithmetic mean and each IC cost is its one verified target interval",
            Self::ArithmeticMeans => "ratio of arithmetic means of per-target three-seed capped expected costs",
        }
    }

    fn estimate(self, pairs: &[(f64, f64)]) -> Option<f64> {
        if pairs.is_empty() {
            return None;
        }
        match self {
            Self::ArithmeticMeans => paired_ratio(
                Some(pairs.iter().map(|pair| pair.0).sum()),
                Some(pairs.iter().map(|pair| pair.1).sum()),
            ),
            Self::GeometricTargetRatios => {
                let mut sum_logs = 0.0;
                for &(numerator, denominator) in pairs {
                    if numerator <= 0.0 || denominator <= 0.0 {
                        return None;
                    }
                    sum_logs += numerator.ln() - denominator.ln();
                }
                let result = (sum_logs / pairs.len() as f64).exp();
                (result.is_finite() && result > 0.0).then_some(result)
            }
        }
    }
}

fn comparison(
    pairs: &[(f64, f64)],
    target_ids: &[String],
    panel_valid: bool,
    seed: u64,
    unit: &str,
    accounting: &str,
    numerator_arm: &str,
    denominator_arm: &str,
    kind: RatioKind,
) -> Value {
    let complete = panel_valid && pairs.len() == EXPECTED_TARGETS;
    let numerator_mean = if complete {
        mean(&pairs.iter().map(|pair| pair.0).collect::<Vec<_>>())
    } else {
        None
    };
    let denominator_mean = if complete {
        mean(&pairs.iter().map(|pair| pair.1).collect::<Vec<_>>())
    } else {
        None
    };
    let ratio = if complete { kind.estimate(pairs) } else { None };
    let interval = if ratio.is_some() {
        bootstrap_interval(pairs, seed, BOOTSTRAP_REPLICATES, kind)
    } else {
        None
    };
    let reason = if !panel_valid {
        Some("Input target panel is incomplete or contains unassignable rows; aggregate and CI withheld.")
    } else if pairs.len() != EXPECTED_TARGETS {
        Some("One or more target pairs lack the required valid draws; aggregate and CI withheld. Missing and failed online draws are not imputed.")
    } else if ratio.is_none() {
        Some("A required ratio has a nonpositive cost or is nonfinite; aggregate and CI withheld.")
    } else if interval.is_none() {
        Some("At least one resample has a nonpositive denominator or nonfinite ratio; CI withheld.")
    } else {
        None
    };
    json!({
        "numerator_arm": numerator_arm,
        "denominator_arm": denominator_arm,
        "accounting": accounting,
        "unit": unit,
        "estimand": kind.name(),
        "ratio": ratio,
        "ratio_of_arithmetic_means_descriptive": paired_ratio(numerator_mean, denominator_mean),
        "numerator_arithmetic_mean": numerator_mean,
        "denominator_arithmetic_mean": denominator_mean,
        "ci95_percentile": interval,
        "reciprocal": {
            "numerator_arm": denominator_arm,
            "denominator_arm": numerator_arm,
            "ratio": paired_ratio(Some(1.0), ratio),
            "ci95_percentile": interval.and_then(|[low, high]| {
                if low > 0.0 && high > 0.0 { Some([1.0 / high, 1.0 / low]) } else { None }
            })
        },
        "complete_target_pairs": pairs.len(),
        "expected_target_pairs": EXPECTED_TARGETS,
        "missing_target_pairs": EXPECTED_TARGETS.saturating_sub(pairs.len()),
        "complete_target_ids": target_ids,
        "complete_panel": complete,
        "no_ci_reason": reason,
        "bootstrap": {
            "unit": "paired target cluster, preserving all seed draws through each target arithmetic mean",
            "replicates": BOOTSTRAP_REPLICATES,
            "seed": seed,
            "prng": "SplitMix64 with unbiased rejection sampling for uniform target indices",
            "pairing": "Identical sampled target indices are used in numerator and denominator.",
            "percentiles": [0.025, 0.975],
            "quantile": "linear interpolation at p*(replicates-1)",
            "scope": "Conditional on the frozen selection, the orbit-exclusion target law, this five-target panel, and this host; three draws do not become fifteen independent target clusters. Target-block resampling is an exchangeability approximation for these five conditionally sampled targets."
        }
    })
}

#[derive(Clone)]
struct SplitMix64(u64);

impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn index(&mut self, len: usize) -> usize {
        let bound = len as u64;
        let threshold = bound.wrapping_neg() % bound;
        loop {
            let draw = self.next();
            if draw >= threshold {
                return (draw % bound) as usize;
            }
        }
    }
}

fn bootstrap_interval(
    pairs: &[(f64, f64)],
    seed: u64,
    replicates: usize,
    kind: RatioKind,
) -> Option<[f64; 2]> {
    if pairs.is_empty() || replicates < 2 {
        return None;
    }
    let mut rng = SplitMix64(seed);
    let mut ratios = Vec::with_capacity(replicates);
    for _ in 0..replicates {
        let mut sampled_pairs = Vec::with_capacity(pairs.len());
        for _ in 0..pairs.len() {
            sampled_pairs.push(pairs[rng.index(pairs.len())]);
        }
        ratios.push(kind.estimate(&sampled_pairs)?);
    }
    ratios.sort_by(f64::total_cmp);
    Some([quantile(&ratios, 0.025), quantile(&ratios, 0.975)])
}

fn quantile(sorted: &[f64], probability: f64) -> f64 {
    let position = probability * (sorted.len() - 1) as f64;
    let low = position.floor() as usize;
    let high = position.ceil() as usize;
    sorted[low] + (sorted[high] - sorted[low]) * (position - low as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn panel() -> (Vec<Value>, Vec<Value>) {
        let mut rows = Vec::new();
        let mut ic = Vec::new();
        for target in 0..EXPECTED_TARGETS {
            let target_id = format!("holdout_{target}");
            for (arm, multiplier) in [("default", 1.0), ("selected", 2.0)] {
                for seed in 0..EXPECTED_RHO_DRAWS {
                    let wall = multiplier * (target + 1) as f64;
                    rows.push(json!({
                        "target_id": target_id, "arm": arm, "seed": seed,
                        "status": VERIFIED, "producer_online_ms": wall * 100.0,
                        "elapsed_process_wall_s": wall, "capped_charge_wall_s": wall,
                        "recovered_scalar": "public-known-answer"
                    }));
                }
            }
            ic.push(json!({
                "target_id": target_id, "target_ms": (target + 1) as f64 * 50.0,
                "group_verified": true, "recovered_scalar": "public-known-answer"
            }));
        }
        (rows, ic)
    }

    #[test]
    fn paired_cluster_resampling_preserves_exact_proportionality() {
        let (rows, ic) = panel();
        let result = summarize(&rows, &ic);
        let primary = &result["primary_default_over_selected_process_wall"];
        assert!((primary["ratio"].as_f64().unwrap() - 0.5).abs() < 1e-12);
        for bound in primary["ci95_percentile"].as_array().unwrap() {
            assert!((bound.as_f64().unwrap() - 0.5).abs() < 1e-12);
        }
        assert_eq!(primary["complete_target_pairs"], json!(5));
        let online = &result["secondary_prepared_ic_over_selected_rho_online"];
        assert!((online["ratio"].as_f64().unwrap() - 0.25).abs() < 1e-12);
        for bound in online["ci95_percentile"].as_array().unwrap() {
            assert!((bound.as_f64().unwrap() - 0.25).abs() < 1e-12);
        }
        assert!((online["reciprocal"]["ratio"].as_f64().unwrap() - 4.0).abs() < 1e-12);
        assert_eq!(result["raw_rows"]["rho"], json!(rows));
    }

    #[test]
    fn failure_is_charged_only_in_secondary_and_disqualifies_full_workload_ratios() {
        let (mut rows, ic) = panel();
        let failed = rows
            .iter_mut()
            .find(|row| row["target_id"] == "holdout_0" && row["arm"] == "selected")
            .unwrap();
        failed["status"] = json!("TIMEOUT");
        failed["producer_online_ms"] = Value::Null;
        failed["elapsed_process_wall_s"] = json!(60.1);
        failed["capped_charge_wall_s"] = json!(60.0);
        let result = summarize(&rows, &ic);
        // Selected sums: target zero (60 + 2 + 2)/3, then 4+6+8+10.
        let expected = ((64.0 / 3.0) + 28.0) / 15.0;
        let capped = &result["secondary_capped_expected_cost_selected_over_default_process_wall"];
        let observed = capped["ratio"].as_f64().unwrap();
        assert!((observed - expected).abs() < 1e-12);
        assert!(!capped["ci95_percentile"].is_null());
        assert!(result["primary_default_over_selected_process_wall"]["ratio"].is_null());
        assert!(result["primary_default_over_selected_process_wall"]["ci95_percentile"].is_null());
        assert!(result["secondary_default_over_selected_online"]["ratio"].is_null());
        assert!(result["secondary_default_over_selected_online"]["ci95_percentile"].is_null());
        assert_eq!(
            result["secondary_default_over_selected_online"]["complete_target_pairs"],
            json!(4)
        );
        assert!(
            result["secondary_prepared_ic_over_selected_rho_online"]["ci95_percentile"].is_null()
        );
        assert_eq!(
            result["per_target"][0]["selected"]["descriptive_verified_online_ms"]["n"],
            json!(2)
        );
    }

    #[test]
    fn missing_seed_draw_is_never_zero_imputed() {
        let (mut rows, ic) = panel();
        rows.remove(0);
        let result = summarize(&rows, &ic);
        let primary = &result["primary_default_over_selected_process_wall"];
        assert_eq!(primary["complete_target_pairs"], json!(4));
        assert!(primary["ratio"].is_null());
        assert!(primary["ci95_percentile"].is_null());
        assert!(result["per_target"][0]["default"]["mean_charged_process_wall_s"].is_null());
        assert_eq!(result["raw_rows"]["rho"].as_array().unwrap().len(), 29);
    }

    #[test]
    fn bad_failure_charge_and_duplicate_seed_fail_closed() {
        let (mut rows, ic) = panel();
        rows[0]["status"] = json!("PROCESS_FAILED");
        rows[0]["capped_charge_wall_s"] = json!(0.0);
        rows[1]["seed"] = rows[2]["seed"].clone();
        let result = summarize(&rows, &ic);
        assert!(result["primary_default_over_selected_process_wall"]["ci95_percentile"].is_null());
        assert!(result["per_target"][0]["default"]["mean_charged_process_wall_s"].is_null());
        assert_eq!(
            result["per_target"][0]["default"]["unique_seed_count"],
            json!(2)
        );
    }

    #[test]
    fn bootstrap_is_reproducible_and_uses_geometric_paired_target_ratios() {
        let pairs = [(1.0, 1.0), (10.0, 2.0), (3.0, 3.0), (4.0, 4.0), (5.0, 5.0)];
        let a = bootstrap_interval(
            &pairs,
            BOOTSTRAP_SEED,
            BOOTSTRAP_REPLICATES,
            RatioKind::GeometricTargetRatios,
        )
        .unwrap();
        let b = bootstrap_interval(
            &pairs,
            BOOTSTRAP_SEED,
            BOOTSTRAP_REPLICATES,
            RatioKind::GeometricTargetRatios,
        )
        .unwrap();
        assert_eq!(a, b);
        let ratio = RatioKind::GeometricTargetRatios.estimate(&pairs).unwrap();
        let ratio_of_means = RatioKind::ArithmeticMeans.estimate(&pairs).unwrap();
        let average_ratio = pairs.iter().map(|p| p.0 / p.1).sum::<f64>() / pairs.len() as f64;
        assert!((ratio - 5.0_f64.powf(0.2)).abs() < 1e-12);
        assert!((ratio - ratio_of_means).abs() > 0.1);
        assert!((ratio - average_ratio).abs() > 0.1);
        // The special target is sampled k times with k~Binomial(5, 1/5).
        // Its population percentile endpoints are k=0 and k=3.
        assert!((a[0] - 1.0).abs() < 1e-12);
        assert!((a[1] - 5.0_f64.powf(0.6)).abs() < 1e-12);
        assert!(bootstrap_interval(
            &[(1.0, 0.0)],
            BOOTSTRAP_SEED,
            100,
            RatioKind::GeometricTargetRatios
        )
        .is_none());
    }

    #[test]
    fn prepared_ic_requires_exactly_one_verified_draw_per_matching_target() {
        let (rows, mut ic) = panel();
        ic[0]["group_verified"] = json!(false);
        let result = summarize(&rows, &ic);
        assert!(
            result["secondary_prepared_ic_over_selected_rho_online"]["ci95_percentile"].is_null()
        );
        assert_eq!(
            result["secondary_prepared_ic_over_selected_rho_online"]["complete_target_pairs"],
            json!(4)
        );
        assert!(!result["secondary_default_over_selected_online"]["ci95_percentile"].is_null());
    }

    #[test]
    fn actual_verified_wall_is_preserved_while_secondary_charge_is_capped() {
        let (mut rows, ic) = panel();
        rows[0]["elapsed_process_wall_s"] = json!(61.0);
        rows[0]["capped_charge_wall_s"] = json!(60.0);
        let result = summarize(&rows, &ic);
        assert_eq!(
            result["per_target"][0]["default"]["mean_verified_actual_process_wall_s"],
            json!(21.0)
        );
        let capped = result["per_target"][0]["default"]["mean_charged_process_wall_s"]
            .as_f64()
            .unwrap();
        assert!((capped - 62.0 / 3.0).abs() < 1e-12);
        assert!(!result["primary_default_over_selected_process_wall"]["ci95_percentile"].is_null());
        assert_eq!(
            result["raw_rows"]["rho"][0]["elapsed_process_wall_s"],
            json!(61.0)
        );
    }

    #[test]
    fn wrong_answer_rejects_entire_holdout_including_capped_secondary() {
        let (mut rows, ic) = panel();
        rows[0]["status"] = json!("WRONG_ANSWER");
        let result = summarize(&rows, &ic);
        assert_eq!(result["wrong_answer_count"], json!(1));
        assert_eq!(result["holdout_soundness_rejected"], json!(true));
        for key in [
            "primary_default_over_selected_process_wall",
            "secondary_default_over_selected_online",
            "secondary_prepared_ic_over_selected_rho_online",
            "secondary_capped_expected_cost_selected_over_default_process_wall",
        ] {
            assert!(result[key]["ratio"].is_null());
            assert!(result[key]["ci95_percentile"].is_null());
        }
    }

    #[test]
    fn paired_arm_seed_sets_must_match_and_text_numbers_are_normalized() {
        let (mut rows, ic) = panel();
        // The same seed with another JSON representation is still a duplicate.
        rows[0]["seed"] = json!("1");
        let duplicate = summarize(&rows, &ic);
        assert_eq!(
            duplicate["per_target"][0]["default"]["unique_seed_count"],
            json!(2)
        );
        assert!(
            duplicate["primary_default_over_selected_process_wall"]["ci95_percentile"].is_null()
        );

        let (mut rows, ic) = panel();
        // Each arm now has three unique seeds, but the paired arm sets differ.
        rows[3]["seed"] = json!(999);
        let mismatch = summarize(&rows, &ic);
        assert_eq!(
            mismatch["per_target"][0]["paired_rho_seed_sets_match"],
            json!(false)
        );
        assert!(mismatch["primary_default_over_selected_process_wall"]["ratio"].is_null());
        assert!(
            mismatch["secondary_prepared_ic_over_selected_rho_online"]["ci95_percentile"].is_null()
        );
    }
}
