//! Native target-cluster statistics for the frozen further-speed campaign.
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
const EXPECTED_TARGETS: usize = 5;
const EXPECTED_RHO_DRAWS: usize = 3;
const FAILURE_CHARGE_S: f64 = 60.0;
const BOOTSTRAP_REPLICATES: usize = 20_000;
const BOOTSTRAP_SEED: u64 = 202_610_126_299;
const VERIFIED: &str = "COMPLETED_VERIFIED";
pub fn summarize_aa(rows: &[Value]) -> Value {
    let mut pairs = BTreeMap::<String, BTreeMap<String, &Value>>::new();
    for row in rows {
        let id = row["run_id"].as_str().unwrap_or("");
        let pair = id.split('-').next().unwrap_or("").to_string();
        pairs
            .entry(pair)
            .or_default()
            .insert(row["arm"].as_str().unwrap_or("").into(), row);
    }
    let mut checks = Vec::new();
    let mut wall = Vec::new();
    let mut online = Vec::new();
    for (pair, arms) in &pairs {
        let a = arms.get("a");
        let b = arms.get("b");
        let equal = a.zip(b).is_some_and(|(a, b)| {
            a["status"] == VERIFIED
                && b["status"] == VERIFIED
                && a["verified"] == true
                && b["verified"] == true
                && a["trajectory"].is_object()
                && a["trajectory"] == b["trajectory"]
                && a["seed"] == b["seed"]
                && a["target_id"] == b["target_id"]
        });
        let ratio = a.zip(b).and_then(|(a, b)| {
            paired_ratio(
                nonnegative(&a["elapsed_process_wall_s"]),
                nonnegative(&b["elapsed_process_wall_s"]),
            )
        });
        let on = a.zip(b).and_then(|(a, b)| {
            paired_ratio(
                nonnegative(&a["producer_online_ms"]),
                nonnegative(&b["producer_online_ms"]),
            )
        });
        if equal {
            if let Some(r) = ratio {
                wall.push(r);
            }
            if let Some(r) = on {
                online.push(r);
            }
        }
        checks.push(json!({"pair_id":pair,"both_verified_matching_exported_counters_and_charges":equal,"a_over_b_process_wall_ratio":ratio,"a_over_b_online_ratio":on}));
    }
    let valid = rows.len() == 10
        && pairs.len() == 5
        && checks
            .iter()
            .all(|c| c["both_verified_matching_exported_counters_and_charges"] == true)
        && wall.len() == 5
        && online.len() == 5;
    json!({"schema":"n83.speed-aa-noise/v1","all_complete_verified_and_counter_equal":valid,"expected_pairs":5,"observed_pairs":pairs.len(),"observed_rows":rows.len(),"pairs":checks,
        "process_wall_a_over_b":if valid {noise_description(&wall)}else{Value::Null},"online_a_over_b":if valid {noise_description(&online)}else{Value::Null},"raw_rows":rows,
        "scope":"Five A/A same-target same-seed fresh-process pairs on training_1,referenceFold64/DP6. Runtime noise diagnostic only; excluded from selection and holdout. No target-population confidence interval."})
}
fn noise_description(v: &[f64]) -> Value {
    let mut sorted = v.to_vec();
    sorted.sort_by(f64::total_cmp);
    let avg = mean(v).unwrap();
    let min = sorted[0];
    let max = sorted[sorted.len() - 1];
    json!({"ratios":v,"arithmetic_mean":avg,"geometric_mean":(v.iter().map(|x|x.ln()).sum::<f64>()/v.len() as f64).exp(),"median":quantile(&sorted,0.5),"minimum":min,"maximum":max,"range":max-min,"sample_standard_deviation":(v.iter().map(|x|(x-avg).powi(2)).sum::<f64>()/(v.len()-1) as f64).sqrt()})
}
#[derive(Default)]
struct TargetRows<'a> {
    current: Vec<&'a Value>,
    selected: Vec<&'a Value>,
}
struct Arm {
    output: Value,
    seed_ids: BTreeSet<u64>,
    actual_process_mean: Option<f64>,
    process_mean: Option<f64>,
    online_mean: Option<f64>,
}
pub fn summarize(rows: &[Value]) -> Value {
    let mut targets = BTreeMap::<String, TargetRows<'_>>::new();
    let mut issues = Vec::<String>::new();
    let mut counts = BTreeMap::<String, usize>::new();
    for (index, row) in rows.iter().enumerate() {
        *counts
            .entry(row["status"].as_str().unwrap_or("MISSING_STATUS").into())
            .or_default() += 1;
        let Some(id) = target_id(row) else {
            issues.push(format!("row{index}:missing target"));
            continue;
        };
        let block = targets.entry(id.into()).or_default();
        match row["arm"].as_str() {
            Some("current") => block.current.push(row),
            Some("selected") => block.selected.push(row),
            _ => issues.push(format!("row{index}:unknown arm")),
        }
    }
    let mut per_target = Vec::new();
    let mut process = Vec::new();
    let mut online = Vec::new();
    let mut charged = Vec::new();
    let mut process_ids = Vec::new();
    let mut online_ids = Vec::new();
    let mut charged_ids = Vec::new();
    for (id, block) in &targets {
        let current = summarize_arm(&block.current);
        let selected = summarize_arm(&block.selected);
        let paired = current.seed_ids.len() == 3 && current.seed_ids == selected.seed_ids;
        if !paired {
            issues.push(format!("{id}:requires same three distinct seeds per arm"));
        }
        for (name, arm) in [("current", &current), ("selected", &selected)] {
            if let Some(errors) = arm.output["issues"].as_array() {
                for e in errors {
                    issues.push(format!("{id}/{name}:{e}"));
                }
            }
        }
        if let (Some(c), Some(s)) = (current.actual_process_mean, selected.actual_process_mean) {
            process.push((c, s));
            process_ids.push(id.clone());
        }
        if let (Some(c), Some(s)) = (current.online_mean, selected.online_mean) {
            online.push((c, s));
            online_ids.push(id.clone());
        }
        if let (Some(c), Some(s)) = (current.process_mean, selected.process_mean) {
            charged.push((c, s));
            charged_ids.push(id.clone());
        }
        per_target.push(json!({"target_id":id,"current":current.output,"selected":selected.output,"paired_seed_sets_match":paired,
            "current_over_selected_process_wall_ratio":paired_ratio(current.actual_process_mean,selected.actual_process_mean),
            "current_over_selected_online_ratio":paired_ratio(current.online_mean,selected.online_mean)}));
    }
    let wrong = rows
        .iter()
        .filter(|r| r["wrong_answer"] == true || r["status"] == "WRONG_ANSWER")
        .count();
    if wrong > 0 {
        issues.push("Wrong answer rejects entire holdout".into());
    }
    if rows.len() != 30 || targets.len() != 5 {
        issues.push("Expected five target blocks and all30rows".into());
    }
    let panel_valid = issues.is_empty() && wrong == 0;
    let all_verified = panel_valid
        && rows
            .iter()
            .all(|r| r["status"] == VERIFIED && r["verified"] == true);
    json!({
        "schema":"n83_speed_statistics_v1",
        "statistical_scope":"Five fresh disjoint public synthetic known-answer target clusters conditional on prior-orbit/base exclusions, on one host. Three paired rho seeds per target. No IC ratio or asymptotic claim.",
        "estimand":"Primary geometric mean of five current/selected ratios of per-target three-seed arithmetic full-process wall means. Current is Fold/reference scalar,64 lanes,DP6,escape off. Selected is frozen on separate training targets.",
        "expected_targets":5,"observed_targets":targets.len(),"observed_rows":rows.len(),"expected_rho_draws_per_target_arm":3,
        "all_verified":all_verified,"status_counts":counts,"wrong_answer_count":wrong,"input_issues":issues,
        "primary_current_over_selected_process_wall":comparison(&process,&process_ids,all_verified,BOOTSTRAP_SEED,"seconds","verified uncensored full-process wall","current","selected",RatioKind::GeometricTargetRatios),
        "secondary_current_over_selected_online":comparison(&online,&online_ids,all_verified,BOOTSTRAP_SEED,"milliseconds","verified target-dependent online interval","current","selected",RatioKind::GeometricTargetRatios),
        "secondary_current_over_selected_arithmetic_process_wall":comparison(&process,&process_ids,all_verified,BOOTSTRAP_SEED,"seconds","ratio of arithmetic means of complete actual wall costs","current","selected",RatioKind::ArithmeticMeans),
        "secondary_current_over_selected_capped_expected_cost":comparison(&charged,&charged_ids,panel_valid,BOOTSTRAP_SEED,"seconds","all attempts: verified min(actual wall,60); each failed attempt60; no missing-draw imputation","current","selected",RatioKind::ArithmeticMeans),
        "per_target":per_target,"raw_rows":rows,
        "accounting":"Whole process includes launch, input, target-independent preparation, target-dependent walk, collision, recovery and producer verification. Polling5ms nominal plus scheduling overhead applied equally. Independent external replay is a separate correctness gate, outside process time. CPU,RSS,outcome counters retained in raw rows.",
        "missing_value_policy":"Every attempted row remains; failures charged60 for secondary capped expected cost. Primary withheld unless all30uncensored verified completions. No zero imputation."
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
            "descriptive_observed_cpu_total_s":descriptive(&rows.iter().filter_map(|r|nonnegative(&r["cpu_total_s"])).collect::<Vec<_>>()),
            "descriptive_observed_peak_rss_bytes":descriptive(&rows.iter().filter_map(|r|nonnegative(&r["peak_rss_bytes"])).collect::<Vec<_>>()),
            "descriptive_completed_walk_steps":descriptive(&rows.iter().filter_map(|r|nonnegative(&r["trajectory"]["outcome"]["walk_steps"])).collect::<Vec<_>>()),
            "descriptive_preparation_ms":descriptive(&rows.iter().filter_map(|r|nonnegative(&r["producer_preparation_ms"])).collect::<Vec<_>>()),
            "issues": issues
        }),
        seed_ids: seeds,
        actual_process_mean,
        process_mean,
        online_mean,
    }
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
            Self::GeometricTargetRatios => "geometric mean of paired per-target cost ratios; each rho cost is its three-seed arithmetic mean",
            Self::ArithmeticMeans => "ratio of arithmetic means of per-target three-seed costs",
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
    fn panel() -> Vec<Value> {
        let mut rows = Vec::new();
        for q in 1..=5 {
            for (arm, multiplier) in [("current", 2.0), ("selected", 1.0)] {
                for seed in 0..3 {
                    let wall = multiplier * q as f64;
                    rows.push(json!({"target_id":format!("holdout_{q}"),"arm":arm,"seed":seed,"status":VERIFIED,"verified":true,"wrong_answer":false,
                        "producer_online_ms":wall*100.,"elapsed_process_wall_s":wall,"capped_charge_wall_s":wall}));
                }
            }
        }
        rows
    }
    #[test]
    fn correct_ratio_and_target_clusters() {
        let r = summarize(&panel());
        assert_eq!(r["all_verified"], true);
        assert!(
            (r["primary_current_over_selected_process_wall"]["ratio"]
                .as_f64()
                .unwrap()
                - 2.)
                .abs()
                < 1e-12
        );
        assert_eq!(
            r["primary_current_over_selected_process_wall"]["bootstrap"]["seed"],
            BOOTSTRAP_SEED
        );
        assert_eq!(
            r["primary_current_over_selected_process_wall"]["complete_target_pairs"],
            5
        );
    }
    #[test]
    fn failure_preserved_charged_and_primary_withheld() {
        let mut p = panel();
        p[0]["status"] = json!("TIMEOUT");
        p[0]["verified"] = json!(false);
        p[0]["capped_charge_wall_s"] = json!(60.);
        let r = summarize(&p);
        assert_eq!(r["all_verified"], false);
        assert!(r["primary_current_over_selected_process_wall"]["ratio"].is_null());
        assert!(r["secondary_current_over_selected_capped_expected_cost"]["ratio"].is_number());
        assert_eq!(r["raw_rows"].as_array().unwrap().len(), 30);
    }
    #[test]
    fn wrong_answer_rejects_every_comparison() {
        let mut p = panel();
        p[0]["status"] = json!("WRONG_ANSWER");
        p[0]["wrong_answer"] = json!(true);
        let r = summarize(&p);
        assert!(r["secondary_current_over_selected_capped_expected_cost"]["ratio"].is_null());
    }
    #[test]
    fn duplicate_seed_rejects_panel() {
        let mut p = panel();
        p[0]["seed"] = p[1]["seed"].clone();
        let r = summarize(&p);
        assert_eq!(r["all_verified"], false);
        assert!(r["primary_current_over_selected_process_wall"]["ratio"].is_null());
    }
    #[test]
    fn missing_draw_is_not_imputed() {
        let mut p = panel();
        p.pop();
        let r = summarize(&p);
        assert!(r["secondary_current_over_selected_capped_expected_cost"]["ratio"].is_null());
    }
    #[test]
    fn false_verified_flag_rejects_primary() {
        let mut p = panel();
        p[0]["verified"] = json!(false);
        let r = summarize(&p);
        assert_eq!(r["all_verified"], false);
    }
    #[test]
    fn aa_rejects_counter_mismatch_and_preserves_all_ten_rows() {
        let mut rows = Vec::new();
        for pair in 1..=5 {
            for arm in ["a", "b"] {
                rows.push(json!({"run_id":format!("AA{pair:02}-{arm}"),"arm":arm,"target_id":"training_1","seed":202610126211u64,"status":VERIFIED,"verified":true,"elapsed_process_wall_s":2.,"producer_online_ms":1000.,"trajectory":{"outcome":{"steps":100}}}));
            }
        }
        let good = summarize_aa(&rows);
        assert_eq!(good["all_complete_verified_and_counter_equal"], true);
        assert_eq!(good["process_wall_a_over_b"]["median"], 1.);
        rows[0]["trajectory"]["outcome"]["steps"] = json!(101);
        let bad = summarize_aa(&rows);
        assert_eq!(bad["all_complete_verified_and_counter_equal"], false);
        assert!(bad["process_wall_a_over_b"].is_null());
        assert_eq!(bad["raw_rows"].as_array().unwrap().len(), 10);
    }
    #[test]
    fn bootstrap_deterministic_and_bounded_by_target_ratios() {
        let pairs = vec![(1., 1.), (2., 1.), (3., 1.), (4., 1.), (5., 1.)];
        let a = bootstrap_interval(
            &pairs,
            BOOTSTRAP_SEED,
            20000,
            RatioKind::GeometricTargetRatios,
        )
        .unwrap();
        let b = bootstrap_interval(
            &pairs,
            BOOTSTRAP_SEED,
            20000,
            RatioKind::GeometricTargetRatios,
        )
        .unwrap();
        assert_eq!(a, b);
        assert!(a[0] >= 1. && a[1] <= 5.);
    }
}
