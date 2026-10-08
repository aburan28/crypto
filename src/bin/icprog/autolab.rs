//! The repository's claim checker, natively: `validate_claim` from the
//! autolab (`research/sat_factor_base_review_20260908/autolab/
//! boundary_autolab.py`), against the ledger its protocol names
//! (`docs/ic/boundary_targets.json`, schema 2), with the same checks in the
//! same order and the same messages.

use std::path::Path;

use super::json::{self, J};
use super::stats;

const TASK_ID: &str = "TASK-IC-BOUNDARY-AUTOLAB-20260910";
const PROTOCOL: &str = "research/sat_factor_base_review_20260908/autolab/protocol.json";

/// Abstract `measurement_schema` keys and the concrete names a report may
/// use for them.
const FIELD_ALIASES: &[(&str, &[&str])] = &[
    ("n_or_bits", &["n_or_bits", "n", "bits"]),
    (
        "factor_base_size_F",
        &["factor_base_size_F", "factor_base_size", "F", "|F|"],
    ),
    (
        "orbit_count_K",
        &["orbit_count_K", "orbit_columns", "K", "orbit_columns_K"],
    ),
    (
        "dimension_l_or_dim",
        &["dimension_l_or_dim", "dimension", "l", "ell", "dim"],
    ),
    ("construction_method", &["construction_method"]),
    ("materialized", &["materialized"]),
    ("construction_wall_ms", &["construction_wall_ms"]),
    ("retained_bytes", &["retained_bytes"]),
    ("m_summands", &["m_summands", "m", "summands"]),
    ("unknowns", &["unknowns"]),
    ("system_degree", &["system_degree"]),
    ("eq_var_ratio", &["eq_var_ratio"]),
    (
        "ffd_or_degree_of_regularity",
        &[
            "ffd_or_degree_of_regularity",
            "ffd",
            "degree_of_regularity",
            "DoR",
            "dor",
        ],
    ),
    ("oracle_class", &["oracle_class"]),
    ("median_ms_per_target", &["median_ms_per_target"]),
    ("largest_solvable", &["largest_solvable"]),
    (
        "base_id_or_hash",
        &["base_id_or_hash", "base_id", "base_hash"],
    ),
    (
        "eta_or_coverage_policy",
        &["eta_or_coverage_policy", "eta", "coverage_policy"],
    ),
    (
        "pr_decomposition_or_hit_rate_with_ci",
        &[
            "pr_decomposition_or_hit_rate_with_ci",
            "pr_decomposition",
            "hit_rate",
            "hit_rate_with_ci",
        ],
    ),
    ("trials_per_relation", &["trials_per_relation"]),
    ("target_mix", &["target_mix"]),
    (
        "orbit_columns_K",
        &["orbit_columns_K", "orbit_columns", "K"],
    ),
    ("relations_collected", &["relations_collected"]),
    ("relations_needed", &["relations_needed"]),
    ("surplus", &["surplus"]),
    ("matrix_dims", &["matrix_dims"]),
    ("sparse_or_dense", &["sparse_or_dense"]),
    ("rank_accumulation", &["rank_accumulation"]),
    (
        "la_wall_ms_or_la_charged_ms",
        &["la_wall_ms_or_la_charged_ms", "la_wall_ms", "la_charged_ms"],
    ),
    (
        "recovered_d_verified",
        &["recovered_d_verified", "recovered_d", "d_verified"],
    ),
    ("stage_timers", &["stage_timers"]),
    ("claim_boundary", &["claim_boundary"]),
    ("timing_class", &["timing_class"]),
    ("target_count", &["target_count"]),
    ("ic_target_hash", &["ic_target_hash"]),
    ("rho_target_hash", &["rho_target_hash"]),
    ("ic_online_wall_ms", &["ic_online_wall_ms"]),
    ("rho_online_wall_ms", &["rho_online_wall_ms"]),
    ("online_speedup", &["online_speedup"]),
    ("online_interval", &["online_interval"]),
    ("same_resource_envelope", &["same_resource_envelope"]),
    ("ic_scalar_verified", &["ic_scalar_verified"]),
    ("rho_scalar_verified", &["rho_scalar_verified"]),
    ("independent_validation", &["independent_validation"]),
    (
        "ic_replay_certificate_sha256",
        &["ic_replay_certificate_sha256"],
    ),
    (
        "rho_replay_certificate_sha256",
        &["rho_replay_certificate_sha256"],
    ),
    ("ic_resource_envelope", &["ic_resource_envelope"]),
    ("rho_resource_envelope", &["rho_resource_envelope"]),
    ("rho_policy", &["rho_policy"]),
    ("ic_cost", &["ic_cost"]),
    ("rho_cost", &["rho_cost"]),
    ("automorphism_discount", &["automorphism_discount"]),
    (
        "all_stages_charged_same_series",
        &["all_stages_charged_same_series"],
    ),
    ("verdict", &["verdict"]),
    (
        "independent_replay_pointer",
        &["independent_replay_pointer"],
    ),
    ("fixture_hash", &["fixture_hash"]),
    (
        "executable_or_source_hash",
        &[
            "executable_or_source_hash",
            "executable_hash",
            "source_hash",
        ],
    ),
    ("host_id", &["host_id"]),
    ("resource_caps", &["resource_caps"]),
    ("seeds", &["seeds", "seed"]),
    (
        "claim_boundary_non_claims",
        &["claim_boundary_non_claims", "non_claims", "claim_boundary"],
    ),
];

/// The stages each arm's online interval must name.
const ONLINE_REQUIRED_STAGES: [(&str, &[&str]); 2] = [
    (
        "ic_included_stages",
        &[
            "target_query",
            "target_PDP",
            "target_relation_check",
            "target_descent",
            "target_recovery_check",
        ],
    ),
    (
        "rho_included_stages",
        &["walk", "collision", "recovery_check"],
    ),
];

fn require(cond: bool, message: String) -> Result<(), String> {
    if cond {
        Ok(())
    } else {
        Err(message)
    }
}

/// `str(v)` as Python prints the values these messages hold.
fn py_str(v: Option<&J>) -> String {
    match v {
        None | Some(J::Null) => "None".into(),
        Some(J::Bool(b)) => if *b { "True" } else { "False" }.into(),
        Some(J::Int(i)) => i.to_string(),
        Some(J::Float(x)) => json::py_float(*x),
        Some(J::Str(s)) => s.clone(),
        Some(other) => json::dumps_line(other, false),
    }
}

/// Python's `str.strip()` left something.
fn nonblank(s: &str) -> bool {
    !s.trim().is_empty()
}

/// The protocol's ledger, refused unless it is the schema the protocol
/// requires and fails closed.
pub fn load_ledger(root: &Path) -> Result<J, String> {
    let protocol = json::read(&root.join(PROTOCOL))?;
    require(
        protocol.get("task_id").and_then(J::as_str) == Some(TASK_ID),
        "protocol task_id mismatch".into(),
    )?;
    let cfg = protocol.at("ledger")?;
    let path = root.join(cfg.at("path")?.as_str().ok_or("the ledger path")?);
    require(
        path.is_file(),
        format!("boundary ledger missing: {}", path.display()),
    )?;
    let ledger = json::read(&path)?;
    let required = cfg.at("required_schema_version")?;
    require(
        ledger
            .get("schema_version")
            .is_some_and(|v| json::py_eq(v, required)),
        format!(
            "ledger schema_version must be {}, got {}",
            py_str(Some(required)),
            py_str(ledger.get("schema_version"))
        ),
    )?;
    let schema = ledger.get("measurement_schema");
    require(
        schema.and_then(J::as_obj).is_some(),
        "ledger missing measurement_schema (fail closed)".into(),
    )?;
    require(
        schema.and_then(|s| s.get("fail_closed")) == Some(&J::Bool(true)),
        "measurement_schema.fail_closed must be true".into(),
    )?;
    Ok(ledger)
}

fn field_present(report: &J, key: &str) -> bool {
    let own = [key];
    let aliases = FIELD_ALIASES
        .iter()
        .find(|(k, _)| *k == key)
        .map_or(&own[..], |(_, a)| a);
    aliases.iter().any(|alias| match report.get(alias) {
        None | Some(J::Null) => false,
        Some(J::Str(s)) => nonblank(s),
        Some(_) => true,
    })
}

fn missing_fields(report: &J, keys: &[J]) -> Vec<J> {
    keys.iter()
        .filter(|k| !field_present(report, k.as_str().unwrap_or("")))
        .cloned()
        .collect()
}

/// Python's `math.isclose`.
fn isclose(a: f64, b: f64, rel_tol: f64, abs_tol: f64) -> bool {
    if a == b {
        return true;
    }
    if a.is_infinite() || b.is_infinite() {
        return false;
    }
    let diff = (b - a).abs();
    diff <= (rel_tol * b).abs() || diff <= (rel_tol * a).abs() || diff <= abs_tol
}

/// The autolab's `positive_cost`: a positive finite number, or nothing.
fn positive_cost(v: Option<&J>) -> Option<f64> {
    match v {
        Some(J::Int(i)) if *i > 0 => Some(*i as f64),
        Some(J::Float(x)) if x.is_finite() && *x > 0.0 => Some(*x),
        _ => None,
    }
}

fn hex64(v: Option<&J>) -> bool {
    v.and_then(J::as_str).is_some_and(|s| {
        s.chars().count() == 64 && s.chars().all(|c| "0123456789abcdef".contains(c))
    })
}

/// The digest an IC1 label ends with: `h([0-9a-f]{12,64})$`, searched.
fn label_digest(id: &str) -> Option<&str> {
    // `$` also matches before one final newline.
    let body = id.strip_suffix('\n').unwrap_or(id);
    let (_, hex) = body.rsplit_once('h')?;
    ((12..=64).contains(&hex.len()) && hex.bytes().all(|b| b"0123456789abcdef".contains(&b)))
        .then_some(hex)
}

/// `validate_claim(report, stage=stage, ledger=ledger)`.
pub fn validate_claim(report: &J, stage: &str, ledger: &J) -> Result<J, String> {
    let schema = ledger.at("measurement_schema")?;
    let stage_schema = schema
        .get(stage)
        .ok_or_else(|| format!("unknown stage for measurement schema: {stage}"))?;
    let list = |v: Option<&J>| v.and_then(J::as_arr).map(<[J]>::to_vec).unwrap_or_default();
    let required = list(stage_schema.get("required"));
    let global_required = list(schema.get("global_provenance_required"));
    let missing_stage = missing_fields(report, &required);
    let missing_global = missing_fields(report, &global_required);
    let mut errors: Vec<String> = Vec::new();

    if stage == "vs_rho" {
        let text = |key: &str| report.get(key).and_then(J::as_str);
        for key in ["candidate_id", "workload_id", "run_id"] {
            if !text(key).is_some_and(nonblank) {
                errors.push(format!("{key} must be a nonempty string"));
            }
        }
        let (candidate_id, workload_id, run_id) =
            (text("candidate_id"), text("workload_id"), text("run_id"));
        let candidate_hash = text("candidate_manifest_sha256");
        let workload_hash = text("workload_manifest_sha256");
        for (key, h) in [
            ("candidate_manifest_sha256", candidate_hash),
            ("workload_manifest_sha256", workload_hash),
        ] {
            match h {
                Some(s) if s.chars().count() == 64 => {
                    if !s.chars().all(|c| "0123456789abcdef".contains(c)) {
                        errors.push(format!("{key} must be lowercase hex"));
                    }
                }
                _ => errors.push(format!("{key} must be a full SHA-256 hex digest")),
            }
        }
        if let (Some(id), Some(h)) = (candidate_id, candidate_hash) {
            match label_digest(id) {
                Some(digest) if id.starts_with("IC1") => {
                    if !h.starts_with(digest) {
                        errors.push(
                            "candidate_id digest must match candidate_manifest_sha256".into(),
                        );
                    }
                }
                _ => errors.push(
                    "candidate_id must use the IC1 identity format with a digest suffix".into(),
                ),
            }
        }
        if let (Some(id), Some(h)) = (workload_id, workload_hash) {
            if id.chars().count() != 12 || !id.chars().all(|c| "0123456789abcdef".contains(c)) {
                errors.push("workload_id must be 12 lowercase hex digits".into());
            } else if !h.starts_with(id) {
                errors.push("workload_id must match workload_manifest_sha256".into());
            }
        }
        if let (Some(c), Some(w), Some(r)) = (candidate_id, workload_id, run_id) {
            if !c.is_empty() && !w.is_empty() && !r.is_empty() {
                let number = r
                    .strip_prefix(c)
                    .and_then(|s| s.strip_prefix('W'))
                    .and_then(|s| s.strip_prefix(w))
                    .and_then(|s| s.strip_prefix('R'));
                let ok = number.is_some_and(|d| {
                    d.bytes().next().is_some_and(|b| (b'1'..=b'9').contains(&b))
                        && d.bytes().all(|b| b.is_ascii_digit())
                });
                if !ok {
                    errors.push("run_id must be <candidate_id>W<workload_id>R<run-number>".into());
                }
            }
        }
        if report.get("target_count") != Some(&J::Int(1)) {
            errors.push("target_count must equal 1".into());
        }
        match (text("ic_target_hash"), text("rho_target_hash")) {
            (Some(a), Some(b)) if nonblank(a) && nonblank(b) => {
                if a != b {
                    errors.push("IC and rho target hashes must match".into());
                }
            }
            _ => errors.push("IC and rho target hashes must be nonempty strings".into()),
        }
        let timing = report.get("timing_class").unwrap_or(&J::Null);
        if !list(stage_schema.get("timing_class_enum"))
            .iter()
            .any(|t| json::py_eq(t, timing))
        {
            errors.push("timing_class must be single_target_online_wall".into());
        }
        let is_true = |key: &str| report.get(key) == Some(&J::Bool(true));
        if !is_true("same_resource_envelope") {
            errors.push("same_resource_envelope must be true".into());
        }
        for key in ["ic_scalar_verified", "rho_scalar_verified"] {
            if !is_true(key) {
                errors.push(format!("{key} must be true"));
            }
        }
        if !is_true("independent_validation") {
            errors.push("independent_validation must be true".into());
        }
        for key in [
            "ic_replay_certificate_sha256",
            "rho_replay_certificate_sha256",
        ] {
            if !hex64(report.get(key)) {
                errors.push(format!("{key} must be a full lowercase SHA-256 digest"));
            }
        }
        let ic_res = report.get("ic_resource_envelope").and_then(J::as_obj);
        let rho_res = report.get("rho_resource_envelope").and_then(J::as_obj);
        if ic_res.is_none_or(|kv| kv.is_empty()) {
            errors.push("ic_resource_envelope must be a nonempty object".into());
        }
        if rho_res.is_none_or(|kv| kv.is_empty()) {
            errors.push("rho_resource_envelope must be a nonempty object".into());
        }
        if let (Some(a), Some(b)) = (ic_res, rho_res) {
            if !json::py_eq(&J::Obj(a.to_vec()), &J::Obj(b.to_vec())) {
                errors.push("IC and rho resource envelopes must match exactly".into());
            }
        }

        let ic_ms = positive_cost(report.get("ic_online_wall_ms"));
        let rho_ms = positive_cost(report.get("rho_online_wall_ms"));
        let speedup = positive_cost(report.get("online_speedup"));
        match (ic_ms, rho_ms, speedup) {
            (Some(ic), Some(rho), Some(s)) => {
                if !isclose(s, rho / ic, 1e-9, 1e-12) {
                    errors.push(
                        "online_speedup must equal rho_online_wall_ms / ic_online_wall_ms".into(),
                    );
                }
            }
            _ => errors.push("online times and speedup must be positive finite numbers".into()),
        }

        match report.get("ic_online_phase_ms") {
            Some(J::Obj(_)) => {
                let costs = report.at("ic_online_phase_ms")?;
                let fields = list(stage_schema.get("ic_online_phase_fields"));
                let mut values = Vec::new();
                for key in &fields {
                    let key = key.as_str().unwrap_or("");
                    match costs.get(key) {
                        Some(J::Int(i)) if *i >= 0 => values.push(*i as f64),
                        Some(J::Float(x)) if x.is_finite() && *x >= 0.0 => values.push(*x),
                        _ => errors.push(format!(
                            "ic_online_phase_ms.{key} must be a finite nonnegative number"
                        )),
                    }
                }
                if let (true, Some(ic)) = (values.len() == fields.len(), ic_ms) {
                    if !isclose(stats::fsum(&values), ic, 1e-6, 1e-3) {
                        errors
                            .push("IC exclusive phase costs must sum to ic_online_wall_ms".into());
                    }
                }
            }
            _ => errors.push("ic_online_phase_ms must be an object".into()),
        }

        match report.get("online_interval") {
            Some(interval @ J::Obj(_)) => {
                let mut keys = list(stage_schema.get("online_interval_event_fields"));
                keys.extend(list(stage_schema.get("online_interval_stage_fields")));
                for key in &keys {
                    let key = key.as_str().unwrap_or("");
                    let value = interval.get(key);
                    let valid = if key.ends_with("_event") {
                        value.and_then(J::as_str).is_some_and(nonblank)
                    } else {
                        value.and_then(J::as_arr).is_some_and(|items| {
                            !items.is_empty()
                                && items.iter().all(|item| item.as_str().is_some_and(nonblank))
                        })
                    };
                    if !valid {
                        errors.push(format!("online_interval.{key} is missing or invalid"));
                    }
                }
                for (included_key, stages) in ONLINE_REQUIRED_STAGES {
                    if let Some(included) = interval.get(included_key).and_then(J::as_arr) {
                        let mut missing: Vec<&str> = stages
                            .iter()
                            .filter(|s| !included.iter().any(|i| i.as_str() == Some(**s)))
                            .copied()
                            .collect();
                        missing.sort_unstable();
                        if !missing.is_empty() {
                            errors.push(format!(
                                "online_interval.{included_key} is missing required stages: {}",
                                missing.join(", ")
                            ));
                        }
                    }
                }
            }
            _ => errors.push("online_interval must be an object".into()),
        }

        match report.get("rho_policy") {
            Some(policy @ J::Obj(_)) => {
                let minimums = stage_schema.get("rho_policy_minimums");
                for key in list(stage_schema.get("rho_policy_integer_fields")) {
                    let key = key.as_str().unwrap_or("");
                    let minimum = minimums.and_then(|m| m.get(key));
                    let ok = match policy.get(key) {
                        Some(J::Int(v)) => match minimum.and_then(J::as_f64) {
                            Some(m) => (*v as f64) >= m,
                            None => true,
                        },
                        _ => false,
                    };
                    if !ok {
                        errors.push(format!(
                            "rho_policy.{key} must be an integer >= {}",
                            py_str(minimum)
                        ));
                    }
                }
                for key in list(stage_schema.get("rho_policy_string_fields")) {
                    let key = key.as_str().unwrap_or("");
                    if !policy.get(key).and_then(J::as_str).is_some_and(nonblank) {
                        errors.push(format!("rho_policy.{key} must be a nonempty string"));
                    }
                }
            }
            _ => errors.push("rho_policy must be an object".into()),
        }
    }

    let ok = missing_stage.is_empty() && missing_global.is_empty() && errors.is_empty();
    Ok(J::Obj(vec![
        (
            "schema_version".into(),
            ledger.get("schema_version").cloned().unwrap_or(J::Null),
        ),
        ("stage".into(), J::Str(stage.into())),
        ("fail_closed".into(), J::Bool(true)),
        (
            "status".into(),
            J::Str(if ok { "PASS" } else { "FAIL" }.into()),
        ),
        ("missing_stage_fields".into(), J::Arr(missing_stage)),
        ("missing_global_provenance".into(), J::Arr(missing_global)),
        (
            "validation_errors".into(),
            J::Arr(errors.into_iter().map(J::Str).collect()),
        ),
        ("required_stage_fields".into(), J::Arr(required)),
        ("required_global_provenance".into(), J::Arr(global_required)),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn label_digests_are_found_as_the_pattern_finds_them() {
        assert_eq!(
            label_digest("IC1N41Ckb0fb6560PDP3pairtableRCaimedLAsparseTDwalkISO0he13eaebe8e18"),
            Some("e13eaebe8e18")
        );
        assert_eq!(label_digest("IC1Nh0123456789a"), None);
        assert_eq!(label_digest("IC1h0123456789ab\n"), Some("0123456789ab"));
        assert!(isclose(1.0, 1.0 + 1e-10, 1e-9, 0.0));
        assert!(!isclose(1.0, 1.0 + 1e-8, 1e-9, 0.0));
    }

    /// A claim §23 checked passes here with the same record, and each
    /// tampering is refused with the autolab's own message.
    #[test]
    fn a_frozen_claim_passes_and_each_tampering_is_named() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let ledger = load_ledger(root).unwrap();
        let doc = json::read(
            &root.join("research/ic_single_target_20260930/claims/k0n41/T01-R1.claim.json"),
        )
        .unwrap();
        let claim = doc.at("claim").unwrap();
        let checked = validate_claim(claim, "vs_rho", &ledger).unwrap();
        assert_eq!(
            json::dumps(&checked, 1),
            json::dumps(doc.at("validation").unwrap(), 1)
        );
        let errors_with = |key: &str, value: J| {
            let mut c = claim.clone();
            if let J::Obj(kv) = &mut c {
                match kv.iter_mut().find(|(k, _)| k == key) {
                    Some((_, v)) => *v = value,
                    None => kv.push((key.into(), value)),
                }
            }
            let r = validate_claim(&c, "vs_rho", &ledger).unwrap();
            assert_eq!(r.at("status").unwrap(), &J::Str("FAIL".into()));
            r.at("validation_errors").unwrap().clone()
        };
        let one = |s: &str| J::Arr(vec![J::Str(s.into())]);
        assert_eq!(
            errors_with("target_count", J::Int(2)),
            one("target_count must equal 1")
        );
        let mut envelope = claim.at("rho_resource_envelope").unwrap().clone();
        if let J::Obj(kv) = &mut envelope {
            kv[0].1 = J::Int(4);
        }
        assert_eq!(
            errors_with("rho_resource_envelope", envelope),
            one("IC and rho resource envelopes must match exactly")
        );
        assert_eq!(
            errors_with("online_speedup", J::Float(9.2)),
            one("online_speedup must equal rho_online_wall_ms / ic_online_wall_ms")
        );
        assert_eq!(
            errors_with("ic_online_wall_ms", J::Float(2.5)),
            J::Arr(vec![
                J::Str("online_speedup must equal rho_online_wall_ms / ic_online_wall_ms".into()),
                J::Str("IC exclusive phase costs must sum to ic_online_wall_ms".into()),
            ])
        );
        assert_eq!(
            errors_with("independent_validation", J::Str("yes".into())),
            one("independent_validation must be true")
        );
        assert_eq!(
            errors_with("run_id", J::Str(format!("{}R1", "IC1Nx"))),
            one("run_id must be <candidate_id>W<workload_id>R<run-number>")
        );
    }
}
