use serde_json::{json, Map, Value};
use std::collections::BTreeSet;

pub const STATES: [&str; 4] = ["unknown", "supported", "refuted", "not_applicable"];

const SCREEN_KEYS: [&str; 6] = [
    "dimension",
    "group_order",
    "frobenius_or_weil_polynomial",
    "field_of_definition",
    "riemann_hurwitz_or_different",
    "polarization_or_jacobian_locus",
];

const OBLIGATION_KEYS: [&str; 9] = [
    "existence",
    "explicit_construction",
    "map_correctness",
    "subgroup_preservation",
    "scalar_recovery",
    "destination_solver_applicability",
    "computational_advantage",
    "scaling_support",
    "novelty_within_boundary",
];

const VERDICT_KEYS: [&str; 5] = [
    "valid_transfer",
    "algorithmically_relevant",
    "modeled_advantage",
    "measured_advantage",
    "novel_reduction_within_boundary",
];

const SOLVER_RECORD_KEYS: [&str; 4] = [
    "applicability",
    "image_subgroup_access",
    "relation_or_search_compatibility",
    "individual_log_and_recovery",
];

const SOLVER_FAMILIES: [&str; 7] = [
    "none",
    "generic_group",
    "index_calculus",
    "weil_descent",
    "pairing_transfer",
    "lattice_or_linear_algebra",
    "other",
];

const PHASE_KEYS: [&str; 9] = [
    "route_search",
    "construction",
    "forward_map",
    "solver_setup",
    "relation_or_walk",
    "linear_algebra",
    "individual_log_or_descent",
    "source_recovery",
    "verification",
];

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Validation {
    pub case_id: String,
    pub errors: Vec<String>,
}

impl Validation {
    pub fn is_valid(&self) -> bool {
        self.errors.is_empty()
    }

    pub fn to_json(&self, assessment: &Value) -> Value {
        json!({
            "case_id": self.case_id,
            "assessment_format": assessment.get("format").cloned().unwrap_or(Value::Null),
            "valid": self.is_valid(),
            "solver": {
                "name": assessment.pointer("/solver_assessment/name").cloned().unwrap_or(Value::Null),
                "family": assessment.pointer("/solver_assessment/family").cloned().unwrap_or(Value::Null),
            },
            "advantage_measurement": {
                "claim_scope": assessment.pointer("/advantage_measurement/claim_scope").cloned().unwrap_or(Value::Null),
                "basis": assessment.pointer("/advantage_measurement/basis").cloned().unwrap_or(Value::Null),
                "net_speedup": assessment.pointer("/advantage_measurement/net_speedup").cloned().unwrap_or(Value::Null),
            },
            "verdicts": assessment.get("verdicts").cloned().unwrap_or(Value::Null),
            "errors": self.errors,
        })
    }
}

fn object_at<'a>(
    value: &'a Value,
    key: &str,
    errors: &mut Vec<String>,
) -> Option<&'a Map<String, Value>> {
    match value.get(key).and_then(Value::as_object) {
        Some(object) => Some(object),
        None => {
            errors.push(format!("{key} must be an object"));
            None
        }
    }
}

fn state_at<'a>(object: &'a Map<String, Value>, key: &str) -> Option<&'a str> {
    object.get(key).and_then(Value::as_str)
}

fn validate_state_record(prefix: &str, key: &str, value: Option<&Value>, errors: &mut Vec<String>) {
    let path = format!("{prefix}.{key}");
    let Some(record) = value.and_then(Value::as_object) else {
        errors.push(format!("{path} must be an object"));
        return;
    };
    let state = state_at(record, "state");
    if !state.is_some_and(|candidate| STATES.contains(&candidate)) {
        errors.push(format!("{path}.state must be one of {}", STATES.join(", ")));
        return;
    }
    let evidence = record.get("evidence").and_then(Value::as_array);
    if evidence.is_none() {
        errors.push(format!("{path}.evidence must be an array"));
        return;
    }
    if matches!(state, Some("supported" | "refuted")) && evidence.is_some_and(Vec::is_empty) {
        errors.push(format!(
            "{path} needs evidence when state is supported or refuted"
        ));
    }
    if !record.get("scope").is_some_and(Value::is_string) {
        errors.push(format!("{path}.scope must be a string"));
    }
}

fn obligation_state<'a>(assessment: &'a Value, key: &str) -> Option<&'a str> {
    assessment
        .get("obligations")?
        .get(key)?
        .get("state")?
        .as_str()
}

fn verdict_state<'a>(assessment: &'a Value, key: &str) -> Option<&'a str> {
    assessment.get("verdicts")?.get(key)?.as_str()
}

fn require_supported(
    assessment: &Value,
    verdict: &str,
    obligations: &[&str],
    prior_verdicts: &[&str],
    errors: &mut Vec<String>,
) {
    if verdict_state(assessment, verdict) != Some("supported") {
        return;
    }
    for obligation in obligations {
        if obligation_state(assessment, obligation) != Some("supported") {
            errors.push(format!(
                "verdicts.{verdict}=supported requires obligations.{obligation}.state=supported"
            ));
        }
    }
    for prior in prior_verdicts {
        if verdict_state(assessment, prior) != Some("supported") {
            errors.push(format!(
                "verdicts.{verdict}=supported requires verdicts.{prior}=supported"
            ));
        }
    }
}

fn finite_positive(value: Option<&Value>) -> Option<f64> {
    value
        .and_then(Value::as_f64)
        .filter(|number| number.is_finite() && *number > 0.0)
}

fn finite_nonnegative(value: Option<&Value>) -> Option<f64> {
    value
        .and_then(Value::as_f64)
        .filter(|number| number.is_finite() && *number >= 0.0)
}

fn validate_cost_arm(prefix: &str, arm: &Value, errors: &mut Vec<String>) {
    let Some(arm) = arm.as_object() else {
        errors.push(format!("{prefix} must be an object or null"));
        return;
    };
    if !arm
        .get("method")
        .and_then(Value::as_str)
        .is_some_and(|method| !method.is_empty())
    {
        errors.push(format!("{prefix}.method must be a nonempty string"));
    }
    if !matches!(
        arm.get("status").and_then(Value::as_str),
        Some("modeled" | "measured")
    ) {
        errors.push(format!("{prefix}.status must be modeled or measured"));
    }
    let total = finite_positive(arm.get("total"));
    if total.is_none() {
        errors.push(format!("{prefix}.total must be a finite positive number"));
    }
    if !arm
        .get("correctness")
        .and_then(Value::as_str)
        .is_some_and(|state| STATES.contains(&state))
    {
        errors.push(format!(
            "{prefix}.correctness must be one of {}",
            STATES.join(", ")
        ));
    }
    if !arm.get("evidence").is_some_and(Value::is_array) {
        errors.push(format!("{prefix}.evidence must be an array"));
    }

    let Some(phases) = arm.get("phases").and_then(Value::as_object) else {
        errors.push(format!("{prefix}.phases must be an object"));
        return;
    };
    let mut phase_sum = 0.0;
    let mut phases_complete = true;
    for phase in PHASE_KEYS {
        match finite_nonnegative(phases.get(phase)) {
            Some(cost) => phase_sum += cost,
            None => {
                phases_complete = false;
                errors.push(format!(
                    "{prefix}.phases.{phase} must be a finite nonnegative number"
                ));
            }
        }
    }
    if let (Some(total), true) = (total, phases_complete) {
        let tolerance = 1e-9_f64.max(total.abs() * 1e-9);
        if (total - phase_sum).abs() > tolerance {
            errors.push(format!(
                "{prefix}.total must equal the sum of the nine phase costs"
            ));
        }
    }
}

fn validate_solver_assessment(assessment: &Value, errors: &mut Vec<String>) {
    let Some(solver) = object_at(assessment, "solver_assessment", errors) else {
        return;
    };
    if !solver
        .get("name")
        .is_some_and(|value| value.is_null() || value.is_string())
    {
        errors.push("solver_assessment.name must be a string or null".into());
    }
    if !solver
        .get("family")
        .and_then(Value::as_str)
        .is_some_and(|family| SOLVER_FAMILIES.contains(&family))
    {
        errors.push(format!(
            "solver_assessment.family must be one of {}",
            SOLVER_FAMILIES.join(", ")
        ));
    }
    if !solver
        .get("implementation")
        .is_some_and(|value| value.is_null() || value.is_string())
    {
        errors.push("solver_assessment.implementation must be a string or null".into());
    }
    for key in SOLVER_RECORD_KEYS {
        validate_state_record("solver_assessment", key, solver.get(key), errors);
    }
}

fn validate_advantage_measurement(assessment: &Value, errors: &mut Vec<String>) {
    let Some(measurement) = object_at(assessment, "advantage_measurement", errors) else {
        return;
    };
    let scope = measurement.get("claim_scope").and_then(Value::as_str);
    if !matches!(
        scope,
        Some("unmeasured" | "stage_only" | "cold_end_to_end" | "amortized")
    ) {
        errors.push(
            "advantage_measurement.claim_scope must be unmeasured, stage_only, cold_end_to_end, or amortized"
                .into(),
        );
    }
    let basis = measurement.get("basis").and_then(Value::as_str);
    if !matches!(basis, Some("unmeasured" | "modeled" | "measured")) {
        errors.push("advantage_measurement.basis must be unmeasured, modeled, or measured".into());
    }
    if matches!(scope, Some("unmeasured")) && basis != Some("unmeasured") {
        errors.push("unmeasured claim_scope requires unmeasured basis".into());
    }
    if matches!(basis, Some("unmeasured")) && scope != Some("unmeasured") {
        errors.push("unmeasured basis requires unmeasured claim_scope".into());
    }
    if !measurement
        .get("unit")
        .is_some_and(|value| value.is_null() || value.is_string())
    {
        errors.push("advantage_measurement.unit must be a string or null".into());
    }
    if !measurement
        .get("target_count")
        .and_then(Value::as_u64)
        .is_some_and(|count| count >= 1)
    {
        errors.push("advantage_measurement.target_count must be a positive integer".into());
    }
    for key in [
        "source_reference",
        "destination_generic_control",
        "destination_candidate",
    ] {
        match measurement.get(key) {
            Some(value) if value.is_null() => {}
            Some(value) => {
                validate_cost_arm(&format!("advantage_measurement.{key}"), value, errors)
            }
            None => errors.push(format!("advantage_measurement.{key} is required")),
        }
    }
    if !measurement
        .get("net_speedup")
        .is_some_and(|value| value.is_null() || finite_positive(Some(value)).is_some())
    {
        errors.push(
            "advantage_measurement.net_speedup must be a finite positive number or null".into(),
        );
    }
    match measurement.get("net_speedup_ci95") {
        Some(value) if value.is_null() => {}
        Some(value) => {
            let bounds = value.as_array();
            let low = bounds.and_then(|values| finite_positive(values.first()));
            let high = bounds.and_then(|values| finite_positive(values.get(1)));
            if bounds.map_or(true, |values| values.len() != 2)
                || low.is_none()
                || high.is_none()
                || low > high
            {
                errors.push(
                    "advantage_measurement.net_speedup_ci95 must be null or [positive low, positive high] with low <= high"
                        .into(),
                );
            }
        }
        None => errors.push("advantage_measurement.net_speedup_ci95 is required".into()),
    }
    for key in ["verified_pairs", "scaling_sizes"] {
        if measurement.get(key).and_then(Value::as_u64).is_none() {
            errors.push(format!(
                "advantage_measurement.{key} must be a nonnegative integer"
            ));
        }
    }
    for key in ["replay_receipts", "unpriced_work", "evidence"] {
        if !measurement.get(key).is_some_and(Value::is_array) {
            errors.push(format!("advantage_measurement.{key} must be an array"));
        }
    }
    if !measurement
        .get("resource_envelope")
        .is_some_and(|value| value.is_null() || value.is_object())
    {
        errors.push("advantage_measurement.resource_envelope must be an object or null".into());
    }
}

fn require_solver_promotion(assessment: &Value, errors: &mut Vec<String>) {
    if verdict_state(assessment, "algorithmically_relevant") != Some("supported") {
        return;
    }
    let solver = assessment.get("solver_assessment");
    if !solver
        .and_then(|value| value.get("name"))
        .and_then(Value::as_str)
        .is_some_and(|name| !name.is_empty())
    {
        errors.push(
            "verdicts.algorithmically_relevant=supported requires a named destination solver"
                .into(),
        );
    }
    if solver
        .and_then(|value| value.get("family"))
        .and_then(Value::as_str)
        == Some("none")
    {
        errors.push(
            "verdicts.algorithmically_relevant=supported requires a non-none solver family".into(),
        );
    }
    for key in SOLVER_RECORD_KEYS {
        if solver
            .and_then(|value| value.get(key))
            .and_then(|value| value.get("state"))
            .and_then(Value::as_str)
            != Some("supported")
        {
            errors.push(format!(
                "verdicts.algorithmically_relevant=supported requires solver_assessment.{key}.state=supported"
            ));
        }
    }
}

fn require_modeled_advantage(assessment: &Value, errors: &mut Vec<String>) {
    if verdict_state(assessment, "modeled_advantage") != Some("supported") {
        return;
    }
    let Some(measurement) = assessment.get("advantage_measurement") else {
        return;
    };
    if measurement.get("claim_scope").and_then(Value::as_str) != Some("cold_end_to_end") {
        errors.push(
            "verdicts.modeled_advantage=supported requires advantage_measurement.claim_scope=cold_end_to_end"
                .into(),
        );
    }
    let basis = measurement.get("basis").and_then(Value::as_str);
    if !matches!(basis, Some("modeled" | "measured")) {
        errors.push(
            "verdicts.modeled_advantage=supported requires a modeled or measured advantage basis"
                .into(),
        );
    }
    if assessment
        .get("cost_comparison")
        .and_then(|value| value.get("modeled_or_measured"))
        .and_then(Value::as_str)
        != basis
    {
        errors.push(
            "verdicts.modeled_advantage=supported requires cost_comparison.modeled_or_measured to match advantage_measurement.basis"
                .into(),
        );
    }
    if !measurement
        .get("unit")
        .and_then(Value::as_str)
        .is_some_and(|unit| !unit.is_empty())
    {
        errors.push("verdicts.modeled_advantage=supported requires a nonempty common unit".into());
    }
    if measurement.get("target_count").and_then(Value::as_u64) != Some(1) {
        errors.push(
            "verdicts.modeled_advantage=supported requires the primary cold one-target workload"
                .into(),
        );
    }
    for key in [
        "source_reference",
        "destination_generic_control",
        "destination_candidate",
    ] {
        let arm = measurement.get(key).and_then(Value::as_object);
        if arm.is_none() {
            errors.push(format!(
                "verdicts.modeled_advantage=supported requires advantage_measurement.{key}"
            ));
        } else {
            if arm.and_then(|value| state_at(value, "correctness")) != Some("supported") {
                errors.push(format!(
                    "verdicts.modeled_advantage=supported requires advantage_measurement.{key}.correctness=supported"
                ));
            }
            if arm
                .and_then(|value| value.get("evidence"))
                .and_then(Value::as_array)
                .map_or(true, Vec::is_empty)
            {
                errors.push(format!(
                    "verdicts.modeled_advantage=supported requires advantage_measurement.{key}.evidence"
                ));
            }
        }
    }
    if measurement
        .get("unpriced_work")
        .and_then(Value::as_array)
        .map_or(true, |work| !work.is_empty())
    {
        errors.push("verdicts.modeled_advantage=supported requires no unpriced work".into());
    }
    if !measurement
        .get("resource_envelope")
        .and_then(Value::as_object)
        .is_some_and(|envelope| !envelope.is_empty())
    {
        errors.push(
            "verdicts.modeled_advantage=supported requires a nonempty resource envelope".into(),
        );
    }
    let source_total = measurement
        .get("source_reference")
        .and_then(|arm| finite_positive(arm.get("total")));
    let candidate_total = measurement
        .get("destination_candidate")
        .and_then(|arm| finite_positive(arm.get("total")));
    let speedup = finite_positive(measurement.get("net_speedup"));
    match (source_total, candidate_total, speedup) {
        (Some(source), Some(candidate), Some(recorded)) => {
            let expected = source / candidate;
            let tolerance = 1e-9_f64.max(expected.abs() * 1e-9);
            if (recorded - expected).abs() > tolerance {
                errors.push(
                    "advantage_measurement.net_speedup must equal source_reference.total / destination_candidate.total"
                        .into(),
                );
            }
            if recorded <= 1.0 {
                errors.push(
                    "verdicts.modeled_advantage=supported requires net_speedup > 1".into(),
                );
            }
        }
        _ => errors.push(
            "verdicts.modeled_advantage=supported requires positive source, candidate, and net-speedup totals"
                .into(),
        ),
    }
    if measurement
        .get("evidence")
        .and_then(Value::as_array)
        .map_or(true, Vec::is_empty)
    {
        errors.push(
            "verdicts.modeled_advantage=supported requires advantage measurement evidence".into(),
        );
    }
}

fn require_measured_advantage(assessment: &Value, errors: &mut Vec<String>) {
    if verdict_state(assessment, "measured_advantage") != Some("supported") {
        return;
    }
    let Some(measurement) = assessment.get("advantage_measurement") else {
        return;
    };
    if measurement.get("basis").and_then(Value::as_str) != Some("measured") {
        errors.push(
            "verdicts.measured_advantage=supported requires advantage_measurement.basis=measured"
                .into(),
        );
    }
    if !assessment
        .pointer("/solver_assessment/implementation")
        .and_then(Value::as_str)
        .is_some_and(|implementation| !implementation.is_empty())
    {
        errors.push(
            "verdicts.measured_advantage=supported requires a nonempty solver implementation revision"
                .into(),
        );
    }
    for key in [
        "source_reference",
        "destination_generic_control",
        "destination_candidate",
    ] {
        if measurement
            .get(key)
            .and_then(|arm| arm.get("status"))
            .and_then(Value::as_str)
            != Some("measured")
        {
            errors.push(format!(
                "verdicts.measured_advantage=supported requires advantage_measurement.{key}.status=measured"
            ));
        }
    }
    if measurement
        .get("verified_pairs")
        .and_then(Value::as_u64)
        .map_or(true, |pairs| pairs < 5)
    {
        errors.push(
            "verdicts.measured_advantage=supported requires at least five verified pairs".into(),
        );
    }
    if measurement
        .get("scaling_sizes")
        .and_then(Value::as_u64)
        .map_or(true, |sizes| sizes < 4)
    {
        errors.push(
            "verdicts.measured_advantage=supported requires at least four measured sizes".into(),
        );
    }
    if measurement
        .get("replay_receipts")
        .and_then(Value::as_array)
        .map_or(true, Vec::is_empty)
    {
        errors.push("verdicts.measured_advantage=supported requires a replay receipt".into());
    }
    let ci_low = measurement
        .get("net_speedup_ci95")
        .and_then(Value::as_array)
        .and_then(|bounds| finite_positive(bounds.first()));
    if ci_low.map_or(true, |low| low <= 1.0) {
        errors.push(
            "verdicts.measured_advantage=supported requires net_speedup_ci95 lower bound > 1"
                .into(),
        );
    }
}

pub fn validate_assessment(assessment: &Value) -> Validation {
    let case_id = assessment
        .get("case_id")
        .and_then(Value::as_str)
        .unwrap_or("<missing-case-id>")
        .to_owned();
    let mut errors = Vec::new();

    if assessment.get("format").and_then(Value::as_str) != Some("correspondence-assessment-v3") {
        errors.push("format must equal correspondence-assessment-v3".into());
    }
    if case_id == "<missing-case-id>" {
        errors.push("case_id must be a string".into());
    }
    for key in ["source", "destination", "transformation", "typed_path"] {
        object_at(assessment, key, &mut errors);
    }

    if let Some(screens) = object_at(assessment, "algebraic_screen", &mut errors) {
        for key in SCREEN_KEYS {
            validate_state_record("algebraic_screen", key, screens.get(key), &mut errors);
        }
    }
    if let Some(obligations) = object_at(assessment, "obligations", &mut errors) {
        for key in OBLIGATION_KEYS {
            validate_state_record("obligations", key, obligations.get(key), &mut errors);
        }
    }
    if let Some(verdicts) = object_at(assessment, "verdicts", &mut errors) {
        for key in VERDICT_KEYS {
            match state_at(verdicts, key) {
                Some(state) if STATES.contains(&state) => {}
                _ => errors.push(format!(
                    "verdicts.{key} must be one of {}",
                    STATES.join(", ")
                )),
            }
        }
    }
    validate_solver_assessment(assessment, &mut errors);
    validate_advantage_measurement(assessment, &mut errors);

    require_supported(
        assessment,
        "valid_transfer",
        &[
            "map_correctness",
            "subgroup_preservation",
            "scalar_recovery",
        ],
        &[],
        &mut errors,
    );
    require_solver_promotion(assessment, &mut errors);
    require_supported(
        assessment,
        "algorithmically_relevant",
        &["destination_solver_applicability"],
        &["valid_transfer"],
        &mut errors,
    );
    require_supported(
        assessment,
        "modeled_advantage",
        &["computational_advantage"],
        &["algorithmically_relevant"],
        &mut errors,
    );
    require_supported(
        assessment,
        "measured_advantage",
        &["scaling_support"],
        &["modeled_advantage"],
        &mut errors,
    );
    require_supported(
        assessment,
        "novel_reduction_within_boundary",
        &["novelty_within_boundary"],
        &["measured_advantage"],
        &mut errors,
    );

    if verdict_state(assessment, "modeled_advantage") == Some("supported") {
        let accounting = assessment.get("cost_comparison").and_then(Value::as_object);
        if accounting.and_then(|v| state_at(v, "accounting_complete")) != Some("supported") {
            errors.push(
                "verdicts.modeled_advantage=supported requires cost_comparison.accounting_complete=supported"
                    .into(),
            );
        }
    }
    if verdict_state(assessment, "measured_advantage") == Some("supported")
        && assessment
            .get("cost_comparison")
            .and_then(|v| v.get("modeled_or_measured"))
            .and_then(Value::as_str)
            != Some("measured")
    {
        errors.push(
            "verdicts.measured_advantage=supported requires cost_comparison.modeled_or_measured=measured"
                .into(),
        );
    }
    require_modeled_advantage(assessment, &mut errors);
    require_measured_advantage(assessment, &mut errors);

    let source_ids: BTreeSet<&str> = assessment
        .get("sources")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(|source| source.get("id").and_then(Value::as_str))
        .collect();
    if source_ids.is_empty() {
        errors.push("sources must contain at least one object with a string id".into());
    }
    for prefix in ["algebraic_screen", "obligations", "solver_assessment"] {
        if let Some(records) = assessment.get(prefix).and_then(Value::as_object) {
            for (key, record) in records {
                let state = record.get("state").and_then(Value::as_str);
                if matches!(state, Some("supported" | "refuted")) {
                    for cited in record
                        .get("evidence")
                        .and_then(Value::as_array)
                        .into_iter()
                        .flatten()
                        .filter_map(Value::as_str)
                    {
                        if !source_ids.contains(cited) {
                            errors.push(format!(
                                "{prefix}.{key}.evidence cites unknown source id {cited}"
                            ));
                        }
                    }
                }
            }
        }
    }
    if let Some(measurement) = assessment
        .get("advantage_measurement")
        .and_then(Value::as_object)
    {
        for cited in measurement
            .get("evidence")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(Value::as_str)
        {
            if !source_ids.contains(cited) {
                errors.push(format!(
                    "advantage_measurement.evidence cites unknown source id {cited}"
                ));
            }
        }
        for key in [
            "source_reference",
            "destination_generic_control",
            "destination_candidate",
        ] {
            for cited in measurement
                .get(key)
                .and_then(Value::as_object)
                .and_then(|arm| arm.get("evidence"))
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
                .filter_map(Value::as_str)
            {
                if !source_ids.contains(cited) {
                    errors.push(format!(
                        "advantage_measurement.{key}.evidence cites unknown source id {cited}"
                    ));
                }
            }
        }
    }

    if assessment
        .get("official_status_change")
        .and_then(Value::as_str)
        != Some("none")
    {
        errors.push("official_status_change must be none for this diagnostic profile".into());
    }

    Validation { case_id, errors }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(state: &str) -> Value {
        json!({"state": state, "evidence": if matches!(state, "supported" | "refuted") { json!(["s1"]) } else { json!([]) }, "scope": "test"})
    }

    fn assessment() -> Value {
        let mut screens = Map::new();
        for key in SCREEN_KEYS {
            screens.insert(key.into(), record("unknown"));
        }
        let mut obligations = Map::new();
        for key in OBLIGATION_KEYS {
            obligations.insert(key.into(), record("unknown"));
        }
        let verdicts = VERDICT_KEYS
            .into_iter()
            .map(|key| (key.to_string(), Value::String("unknown".into())))
            .collect::<Map<String, Value>>();
        json!({
            "format": "correspondence-assessment-v3",
            "case_id": "control",
            "source": {}, "destination": {}, "transformation": {}, "typed_path": {},
            "algebraic_screen": screens,
            "obligations": obligations,
            "solver_assessment": {
                "name": null,
                "family": "none",
                "implementation": null,
                "applicability": record("unknown"),
                "image_subgroup_access": record("unknown"),
                "relation_or_search_compatibility": record("unknown"),
                "individual_log_and_recovery": record("unknown")
            },
            "cost_comparison": {"accounting_complete": "unknown", "modeled_or_measured": "unmeasured"},
            "advantage_measurement": {
                "claim_scope": "unmeasured",
                "basis": "unmeasured",
                "unit": null,
                "target_count": 1,
                "source_reference": null,
                "destination_generic_control": null,
                "destination_candidate": null,
                "net_speedup": null,
                "net_speedup_ci95": null,
                "verified_pairs": 0,
                "scaling_sizes": 0,
                "replay_receipts": [],
                "unpriced_work": [],
                "resource_envelope": null,
                "evidence": []
            },
            "sources": [{"id": "s1"}],
            "verdicts": verdicts,
            "official_status_change": "none"
        })
    }

    fn phases(total: f64) -> Value {
        json!({
            "route_search": 0.0,
            "construction": 0.0,
            "forward_map": 0.0,
            "solver_setup": 0.0,
            "relation_or_walk": total,
            "linear_algebra": 0.0,
            "individual_log_or_descent": 0.0,
            "source_recovery": 0.0,
            "verification": 0.0
        })
    }

    fn arm(method: &str, status: &str, total: f64) -> Value {
        json!({
            "method": method,
            "status": status,
            "total": total,
            "phases": phases(total),
            "correctness": "supported",
            "evidence": ["s1"]
        })
    }

    fn support_transfer_and_solver(value: &mut Value) {
        for key in [
            "map_correctness",
            "subgroup_preservation",
            "scalar_recovery",
            "destination_solver_applicability",
        ] {
            value["obligations"][key] = record("supported");
        }
        value["verdicts"]["valid_transfer"] = json!("supported");
        value["verdicts"]["algorithmically_relevant"] = json!("supported");
        value["solver_assessment"] = json!({
            "name": "control solver",
            "family": "index_calculus",
            "implementation": "control-revision",
            "applicability": record("supported"),
            "image_subgroup_access": record("supported"),
            "relation_or_search_compatibility": record("supported"),
            "individual_log_and_recovery": record("supported")
        });
    }

    fn support_modeled_advantage(value: &mut Value, basis: &str, arm_status: &str) {
        support_transfer_and_solver(value);
        value["obligations"]["computational_advantage"] = record("supported");
        value["verdicts"]["modeled_advantage"] = json!("supported");
        value["cost_comparison"] = json!({
            "accounting_complete": "supported",
            "modeled_or_measured": basis
        });
        value["advantage_measurement"] = json!({
            "claim_scope": "cold_end_to_end",
            "basis": basis,
            "unit": "control_operations",
            "target_count": 1,
            "source_reference": arm("source", arm_status, 120.0),
            "destination_generic_control": arm("destination generic", arm_status, 110.0),
            "destination_candidate": arm("destination candidate", arm_status, 80.0),
            "net_speedup": 1.5,
            "net_speedup_ci95": null,
            "verified_pairs": 0,
            "scaling_sizes": 0,
            "replay_receipts": [],
            "unpriced_work": [],
            "resource_envelope": {"threads": 1, "memory_bytes": 1024},
            "evidence": ["s1"]
        });
    }

    #[test]
    fn accepts_fail_closed_unknown_assessment() {
        assert!(validate_assessment(&assessment()).is_valid());
    }

    #[test]
    fn accepts_supported_transfer_with_all_hard_gates() {
        let mut value = assessment();
        for key in [
            "map_correctness",
            "subgroup_preservation",
            "scalar_recovery",
        ] {
            value["obligations"][key] = record("supported");
        }
        value["verdicts"]["valid_transfer"] = json!("supported");
        assert!(validate_assessment(&value).is_valid());
    }

    #[test]
    fn rejects_promoted_transfer_without_hard_gates() {
        let mut value = assessment();
        value["verdicts"]["valid_transfer"] = json!("supported");
        let result = validate_assessment(&value);
        assert!(!result.is_valid());
        assert_eq!(result.errors.len(), 3);
    }

    #[test]
    fn rejects_unbound_evidence_and_unsupported_measurement_promotion() {
        let mut value = assessment();
        value["obligations"]["map_correctness"] =
            json!({"state":"supported", "evidence":["missing"], "scope":"test"});
        value["verdicts"]["measured_advantage"] = json!("supported");
        let result = validate_assessment(&value);
        assert!(!result.is_valid());
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("unknown source id")));
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("modeled_or_measured=measured")));
    }

    #[test]
    fn rejects_algorithmic_relevance_without_named_solver_gates() {
        let mut value = assessment();
        for key in [
            "map_correctness",
            "subgroup_preservation",
            "scalar_recovery",
            "destination_solver_applicability",
        ] {
            value["obligations"][key] = record("supported");
        }
        value["verdicts"]["valid_transfer"] = json!("supported");
        value["verdicts"]["algorithmically_relevant"] = json!("supported");
        let result = validate_assessment(&value);
        assert!(!result.is_valid());
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("named destination solver")));
    }

    #[test]
    fn accepts_complete_modeled_net_advantage() {
        let mut value = assessment();
        support_modeled_advantage(&mut value, "modeled", "modeled");
        assert!(validate_assessment(&value).is_valid());
    }

    #[test]
    fn rejects_stage_only_or_unpriced_modeled_promotion() {
        let mut value = assessment();
        support_modeled_advantage(&mut value, "modeled", "modeled");
        value["advantage_measurement"]["claim_scope"] = json!("stage_only");
        value["advantage_measurement"]["unpriced_work"] = json!(["linear algebra"]);
        value["advantage_measurement"]["destination_candidate"]["evidence"] = json!([]);
        value["cost_comparison"]["modeled_or_measured"] = json!("unmeasured");
        let result = validate_assessment(&value);
        assert!(!result.is_valid());
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("claim_scope=cold_end_to_end")));
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("no unpriced work")));
        assert!(result.errors.iter().any(|error| error.contains("to match")));
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("destination_candidate.evidence")));
    }

    #[test]
    fn rejects_inconsistent_phase_sum_and_net_ratio() {
        let mut value = assessment();
        support_modeled_advantage(&mut value, "modeled", "modeled");
        value["advantage_measurement"]["destination_candidate"]["total"] = json!(70.0);
        let result = validate_assessment(&value);
        assert!(!result.is_valid());
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("sum of the nine phase costs")));
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("source_reference.total")));
    }

    #[test]
    fn accepts_complete_measured_advantage_control() {
        let mut value = assessment();
        support_modeled_advantage(&mut value, "measured", "measured");
        value["obligations"]["scaling_support"] = record("supported");
        value["verdicts"]["measured_advantage"] = json!("supported");
        value["advantage_measurement"]["verified_pairs"] = json!(8);
        value["advantage_measurement"]["scaling_sizes"] = json!(4);
        value["advantage_measurement"]["replay_receipts"] = json!(["receipt-sha256"]);
        value["advantage_measurement"]["net_speedup_ci95"] = json!([1.1, 1.9]);
        assert!(validate_assessment(&value).is_valid());
    }

    #[test]
    fn rejects_measured_promotion_without_replay_scaling_or_interval() {
        let mut value = assessment();
        support_modeled_advantage(&mut value, "measured", "measured");
        value["obligations"]["scaling_support"] = record("supported");
        value["verdicts"]["measured_advantage"] = json!("supported");
        value["advantage_measurement"]["verified_pairs"] = json!(4);
        value["advantage_measurement"]["scaling_sizes"] = json!(3);
        value["advantage_measurement"]["net_speedup_ci95"] = json!([0.9, 1.9]);
        value["solver_assessment"]["implementation"] = Value::Null;
        let result = validate_assessment(&value);
        assert!(!result.is_valid());
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("five verified pairs")));
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("four measured sizes")));
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("replay receipt")));
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("lower bound > 1")));
        assert!(result
            .errors
            .iter()
            .any(|error| error.contains("implementation revision")));
    }
}
