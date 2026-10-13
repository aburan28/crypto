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
            "valid": self.is_valid(),
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

pub fn validate_assessment(assessment: &Value) -> Validation {
    let case_id = assessment
        .get("case_id")
        .and_then(Value::as_str)
        .unwrap_or("<missing-case-id>")
        .to_owned();
    let mut errors = Vec::new();

    if assessment.get("format").and_then(Value::as_str) != Some("correspondence-assessment-v2") {
        errors.push("format must equal correspondence-assessment-v2".into());
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
    for prefix in ["algebraic_screen", "obligations"] {
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
            "format": "correspondence-assessment-v2",
            "case_id": "control",
            "source": {}, "destination": {}, "transformation": {}, "typed_path": {},
            "algebraic_screen": screens,
            "obligations": obligations,
            "cost_comparison": {"accounting_complete": "unknown", "modeled_or_measured": "unmeasured"},
            "sources": [{"id": "s1"}],
            "verdicts": verdicts,
            "official_status_change": "none"
        })
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
}
