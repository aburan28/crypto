//! Native, fail-closed checks for the non-authoritative ECDLP audit plan.
//!
//! This command verifies the planning contract. It does not promote a
//! mathematical, implementation, performance, historical, or security claim.

use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde::de::{self, MapAccess, SeqAccess, Visitor};
use serde::{Deserialize, Deserializer};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Component, Path, PathBuf};
use std::process::Command as ProcessCommand;

#[path = "ecdlp_audit/fixed_f8_order.rs"]
mod fixed_f8_order;

const EXPECTED_SCHEMA_VERSION: &str = "ecdlp-audit-plan/v1";
const EXPECTED_RECORD_STATE: &str = "mock_non_authoritative";
const EXPECTED_WORKSTREAMS: [&str; 8] = ["WS1", "WS2", "WS3", "WS4", "WS5", "WS6", "WS7", "WS8"];
const EXPECTED_LANES: [&str; 5] = ["A", "B", "C", "D", "E"];
const EXPECTED_SCOPE_CLAIM_BOUNDARY: &str =
    "This is a planning mock, not an authoritative ledger, attack result, completion claim, authorization instrument, or full Jacobian index-calculus solver.";
const EXPECTED_REPORT_CLAIM_BOUNDARY: &str =
    "The mock verifies the planning contract only. It cannot promote an attack, safety, completeness, arithmetic, transfer, solver, performance, or historical-superlative claim.";
const PLAN_SCHEMA_JSON: &str = include_str!("../../docs/ecdlp-audit/ecdlp-audit-plan.schema.json");
const MAX_PLAN_BYTES: u64 = 4 * 1024 * 1024;
const MAX_SOURCE_BYTES: u64 = 64 * 1024 * 1024;
const EXPECTED_RECORD_PATH: &str = "docs/ecdlp-audit/ECDLP_AUDIT_PLAN.mock.json";
const EXPECTED_PLAN_PATH: &str =
    "research/notes/ecdlp-general/RESEARCH_COMPREHENSIVE_ECDLP_AUDIT_PLAN.md";
const EXPECTED_PDF_PATH: &str =
    "research/notes/ecdlp-general/RESEARCH_COMPREHENSIVE_ECDLP_AUDIT_PLAN.pdf";
const EXPECTED_DIAGRAM_PATH: &str = "docs/ecdlp-audit/program-dependencies.svg";
const EXPECTED_VALIDATOR_PATH: &str = "src/bin/ecdlp_audit.rs";
const EXPECTED_DELIVERABLE_BINDINGS: [(&str, &str, Option<&str>); 19] = [
    ("WS1-D1", "artifact_replay", None),
    (
        "WS1-D2",
        "review_receipt",
        Some("research/cover_transfer_assessment_20261013/README.md"),
    ),
    ("WS2-D1", "implementation", None),
    (
        "WS2-D2",
        "implementation",
        Some("src/binary_ecc/cover_transfer.rs"),
    ),
    (
        "WS2-D3",
        "certificate_corpus",
        Some("problems/catalog-cover-subgroup-certificates-incomplete.md"),
    ),
    ("WS2-D4", "independent_oracle", None),
    ("WS3-D1", "schema_and_ledger", None),
    ("WS3-D2", "semantic_verifier", None),
    (
        "WS3-D3",
        "identity_audit",
        Some("docs/ecdlp-audit/ledger/IDENTITY_CANONICALIZATION_AUDIT.md"),
    ),
    (
        "WS4-D1",
        "evidence_profiles",
        Some("docs/ecdlp-audit/evidence/evidence-profile-registry.draft.json"),
    ),
    (
        "WS4-D2",
        "transition_validator",
        Some("docs/ecdlp-audit/evidence/README.md"),
    ),
    (
        "WS5-D1",
        "route_registry",
        Some("docs/ecdlp-audit/routes/route-registry.seed.json"),
    ),
    (
        "WS5-D2",
        "coverage_matrix",
        Some("docs/ecdlp-audit/routes/fixtures/valid/coverage-matrix.seed.json"),
    ),
    (
        "WS6-D1",
        "search_protocol",
        Some("docs/ecdlp-audit/corpus/PROTOCOL.md"),
    ),
    (
        "WS6-D2",
        "source_corpus",
        Some("docs/ecdlp-audit/corpus/seed-manifest.json"),
    ),
    (
        "WS7-D1",
        "run_schema",
        Some("docs/ecdlp-audit/accounting/cross-route-cost-receipt.schema.json"),
    ),
    (
        "WS7-D2",
        "claim_checker",
        Some("docs/ecdlp-audit/accounting/README.md"),
    ),
    (
        "WS8-D1",
        "governance_policy",
        Some("docs/ecdlp-audit/governance/POLICY.md"),
    ),
    (
        "WS8-D2",
        "governance_drill",
        Some("docs/ecdlp-audit/governance/README.md"),
    ),
];
const EXPECTED_PINNED_SOURCES: [(&str, &str, &str, &str); 3] = [
    (
        "source.cover-transfer-ledger",
        "docs/curves/cover-transfer-assessments.json",
        "da6ae6d1b030f6b7455bfb042634520d2e93b7aa",
        "9bdf27316c0370ef4e7013d8d36186ec1c7d186858b76c6ac885475df1a456e9",
    ),
    (
        "source.ic-boundary-ledger",
        "docs/ic/boundary_targets.json",
        "2906ecf3f73936eaa14f4b87166d908c710ed2ef",
        "1e062aa4fed838279028871b60a3b23a018d18d4d60d07a0e89700d2ae4aba64",
    ),
    (
        "source.attack-matrix-seed",
        "docs/ECDLP_ATTACK_MATRIX.md",
        "da6ae6d1b030f6b7455bfb042634520d2e93b7aa",
        "42cd637a56339f96cc12ba9e0c7d0c7c28d08bcb488fe69ced07fe561bef39e4",
    ),
];

#[derive(Parser)]
#[command(name = "ecdlp-audit")]
#[command(about = "Validate ECDLP audit control-plane records")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Verify the non-authoritative eight-workstream planning record.
    VerifyPlan {
        /// Planning JSON, relative to --repo-root unless absolute.
        #[arg(long, default_value = "docs/ecdlp-audit/ECDLP_AUDIT_PLAN.mock.json")]
        input: PathBuf,
        /// Repository root used to resolve bound evidence and report artifacts.
        #[arg(long, default_value = ".")]
        repo_root: PathBuf,
    },
    /// Independently replay the exact fixed-GF(8) genus-3 order certificate.
    VerifyFixedF8Order {
        /// Certificate JSON, relative to --repo-root unless absolute.
        #[arg(
            long,
            default_value = "docs/curves/fixed-f8-genus3-order-certificate.json"
        )]
        input: PathBuf,
        /// Repository root used to confine and resolve the certificate path.
        #[arg(long, default_value = ".")]
        repo_root: PathBuf,
    },
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct Plan {
    schema_version: String,
    record_state: String,
    plan_uid: String,
    title: String,
    as_of: String,
    authoritative: bool,
    promotion_allowed: bool,
    scope: Scope,
    source_artifacts: Vec<SourceArtifact>,
    vocabularies: Vocabularies,
    inventory_freeze: InventoryFreeze,
    workstreams: Vec<Workstream>,
    attack_routes: Vec<AttackRoute>,
    subjects: Vec<Subject>,
    assessments: Vec<Assessment>,
    dependencies: Vec<Dependency>,
    completion_metrics: CompletionMetrics,
    governance: Governance,
    open_obligations: Vec<Obligation>,
    report: Report,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct Scope {
    objective: String,
    claim_boundary: String,
    execution_fixture_classes: Vec<String>,
    default_authorization_tiers: Vec<String>,
    live_target_execution: bool,
    secret_key_material: bool,
    production_exploitation: bool,
    outcome_classes: Vec<String>,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct Vocabularies {
    claim_status: Vec<String>,
    requirement_status: Vec<String>,
    maturity_levels: Vec<String>,
    independence_levels: Vec<String>,
    authorization_tiers: Vec<String>,
    cost_status: Vec<String>,
    execution_outcome: Vec<String>,
    review_decision: Vec<String>,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct SourceArtifact {
    id: String,
    path: String,
    role: String,
    revision: Option<String>,
    sha256: Option<String>,
    status: SourceStatus,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq)]
#[serde(rename_all = "snake_case")]
enum SourceStatus {
    Pinned,
    CurrentWorktree,
    Planned,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct InventoryFreeze {
    status: String,
    frozen_at: Option<String>,
    unit_expression: String,
    denominator_manifest_sha256: Option<String>,
    expected_units: Option<u64>,
    note: String,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct Workstream {
    id: String,
    slug: String,
    title: String,
    requirement_status: RequirementStatus,
    depends_on: Vec<String>,
    summary: String,
    deliverables: Vec<Deliverable>,
    acceptance_gates: Vec<AcceptanceGate>,
    metrics: WorkstreamMetrics,
    non_claims: Vec<String>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq)]
#[serde(rename_all = "snake_case")]
enum RequirementStatus {
    VerifiedComplete,
    ImplementedUnverified,
    Partial,
    Blocked,
    NotAttempted,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct Deliverable {
    id: String,
    kind: String,
    description: String,
    path: Option<String>,
    artifact_sha256: Option<String>,
    status: RequirementStatus,
    acceptance_gate_ids: Vec<String>,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct AcceptanceGate {
    id: String,
    required: bool,
    claim_status: ClaimStatus,
    scope: String,
    rationale: String,
    evidence_refs: Vec<String>,
    obligation_id: Option<String>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq)]
#[serde(rename_all = "snake_case")]
enum ClaimStatus {
    Unknown,
    Supported,
    Refuted,
    NotApplicable,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
#[serde(deny_unknown_fields)]
struct WorkstreamMetrics {
    required_gates: usize,
    classified_gates: usize,
    supported_gates: usize,
    open_required_gates: usize,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct AttackRoute {
    id: String,
    lane: String,
    title: String,
    claim_status: ClaimStatus,
    scope: String,
    obligation_id: String,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct Subject {
    unit_uid: String,
    unit_key: Value,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct Assessment {
    assessment_uid: String,
    unit_uid: String,
    revision: u64,
    supersedes: Option<String>,
    claim_status: ClaimStatus,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct Dependency {
    from: String,
    to: String,
    relation: DependencyRelation,
    status: DependencyStatus,
    rationale: String,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd)]
#[serde(rename_all = "snake_case")]
enum DependencyRelation {
    Requires,
    Blocks,
    Validates,
    Governs,
    Supplies,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd)]
#[serde(rename_all = "snake_case")]
enum DependencyStatus {
    Active,
    Resolved,
    Planned,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct CompletionMetrics {
    inventory: CountBlock,
    classification: ClassificationMetrics,
    evidence: EvidenceMetrics,
    authorization: AuthorizationMetrics,
    governance: GovernanceMetrics,
    obligations: ObligationMetrics,
    derived: DerivedMetrics,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct CountBlock {
    expected: Option<u64>,
    present: u64,
    missing: Option<u64>,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct ClassificationMetrics {
    required_claim_cells: Option<u64>,
    classified_claim_cells: u64,
    unknown: u64,
    supported: u64,
    refuted: u64,
    not_applicable: u64,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct EvidenceMetrics {
    independently_replayed: u64,
    matched_controlled: u64,
    end_to_end_certified: u64,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct AuthorizationMetrics {
    execution_units: u64,
    documented: u64,
    blocked_or_unknown: u64,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct GovernanceMetrics {
    promotion_candidates: u64,
    fully_reviewed: u64,
    conflicted_reviews: u64,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ObligationMetrics {
    open_required: usize,
    blocked: usize,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct DerivedMetrics {
    inventory_complete: bool,
    classification_complete: bool,
    authorization_complete: bool,
    promotion_ready: bool,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct Governance {
    status: String,
    required_roles: Vec<String>,
    self_promotion_allowed: bool,
    immutable_evidence: bool,
    promotion_rule: String,
    correction_rule: String,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct Obligation {
    id: String,
    workstream_id: String,
    status: ObligationStatus,
    required: bool,
    description: String,
    exit_evidence: String,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq)]
#[serde(rename_all = "snake_case")]
enum ObligationStatus {
    Open,
    Blocked,
    Resolved,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct Report {
    record_path: String,
    plan_path: String,
    pdf_path: String,
    diagram_path: String,
    validator_path: String,
    claim_boundary: String,
}

fn main() {
    let result = match Cli::parse().command {
        Command::VerifyPlan { input, repo_root } => verify_plan_file(&input, &repo_root),
        Command::VerifyFixedF8Order { input, repo_root } => {
            fixed_f8_order::verify_certificate_file(&input, &repo_root)
                .map(|result| println!("{result}"))
                .map_err(|error| vec![error])
        }
    };

    if let Err(errors) = result {
        for error in &errors {
            eprintln!("ERROR: {error}");
        }
        eprintln!(
            "ECDLP audit verification failed with {} error(s).",
            errors.len()
        );
        std::process::exit(1);
    }
}

fn verify_plan_file(input: &Path, repo_root: &Path) -> Result<(), Vec<String>> {
    let root = match fs::canonicalize(repo_root) {
        Ok(root) => root,
        Err(error) => return Err(vec![format!("cannot resolve repository root: {error}")]),
    };
    let unresolved_input = if input.is_absolute() {
        input.to_path_buf()
    } else {
        root.join(input)
    };
    let input_path = match fs::canonicalize(&unresolved_input) {
        Ok(path) if path.starts_with(&root) && path.is_file() => path,
        Ok(path) if !path.starts_with(&root) => {
            return Err(vec![format!(
                "input plan resolves outside repository: {}",
                path.display()
            )])
        }
        Ok(path) => {
            return Err(vec![format!(
                "input plan is not a file: {}",
                path.display()
            )])
        }
        Err(error) => {
            return Err(vec![format!(
                "cannot resolve input plan {}: {error}",
                unresolved_input.display()
            )])
        }
    };
    let bytes = match read_bounded(&input_path, MAX_PLAN_BYTES, "input plan") {
        Ok(bytes) => bytes,
        Err(error) => return Err(vec![error]),
    };
    let raw_plan = match parse_strict_json(&bytes) {
        Ok(value) => value,
        Err(error) => {
            return Err(vec![format!(
                "cannot parse plan identity preimage: {error}"
            )])
        }
    };
    let schema_errors = validate_plan_schema(&raw_plan);
    if !schema_errors.is_empty() {
        return Err(schema_errors);
    }
    let plan: Plan = match serde_json::from_value(raw_plan.clone()) {
        Ok(plan) => plan,
        Err(error) => return Err(vec![format!("strict typed plan parse failed: {error}")]),
    };
    let errors = validate_plan(&plan, &raw_plan, &root, &input_path);
    if !errors.is_empty() {
        return Err(errors);
    }

    let required_gates = plan
        .workstreams
        .iter()
        .flat_map(|workstream| &workstream.acceptance_gates)
        .filter(|gate| gate.required)
        .count();
    let open_required = plan.completion_metrics.obligations.open_required;
    println!("ECDLP audit plan verification: PASS");
    println!("plan_uid={}", plan.plan_uid);
    println!(
        "workstreams={} required_gates={required_gates} open_required={open_required}",
        plan.workstreams.len()
    );
    println!("authoritative=false promotion_allowed=false subjects=0 assessments=0");
    Ok(())
}

fn validate_plan_schema(instance: &Value) -> Vec<String> {
    let schema = match parse_strict_json(PLAN_SCHEMA_JSON.as_bytes()) {
        Ok(schema) => schema,
        Err(error) => return vec![format!("embedded plan schema is not JSON: {error}")],
    };
    let schema_errors = validate_supported_plan_schema(&schema);
    if !schema_errors.is_empty() {
        return schema_errors;
    }

    let mut violations = Vec::new();
    validate_schema_instance(
        instance,
        &schema,
        &schema,
        "",
        "#",
        &mut Vec::new(),
        &mut violations,
    );
    let truncated = violations.len() == MAX_SCHEMA_VIOLATIONS;
    let mut errors: Vec<String> = violations
        .into_iter()
        .map(|violation| {
            format!(
                "schema violation at {} (schema {}): {}",
                if violation.instance_path.is_empty() {
                    "<root>"
                } else {
                    &violation.instance_path
                },
                violation.schema_path,
                violation.message
            )
        })
        .collect();
    if truncated {
        errors.push(format!(
            "schema validation stopped after {MAX_SCHEMA_VIOLATIONS} violations; additional violations may exist"
        ));
    }
    errors.sort();
    errors.dedup();
    errors
}

struct DuplicateRejectingJson;

impl<'de> Deserialize<'de> for DuplicateRejectingJson {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_any(DuplicateRejectingJsonVisitor)
    }
}

struct DuplicateRejectingJsonVisitor;

impl<'de> Visitor<'de> for DuplicateRejectingJsonVisitor {
    type Value = DuplicateRejectingJson;

    fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str("JSON without duplicate object keys")
    }

    fn visit_unit<E>(self) -> Result<Self::Value, E> {
        Ok(DuplicateRejectingJson)
    }

    fn visit_bool<E>(self, _value: bool) -> Result<Self::Value, E> {
        Ok(DuplicateRejectingJson)
    }

    fn visit_i64<E>(self, _value: i64) -> Result<Self::Value, E> {
        Ok(DuplicateRejectingJson)
    }

    fn visit_u64<E>(self, _value: u64) -> Result<Self::Value, E> {
        Ok(DuplicateRejectingJson)
    }

    fn visit_i128<E>(self, _value: i128) -> Result<Self::Value, E> {
        Ok(DuplicateRejectingJson)
    }

    fn visit_u128<E>(self, _value: u128) -> Result<Self::Value, E> {
        Ok(DuplicateRejectingJson)
    }

    fn visit_f64<E>(self, _value: f64) -> Result<Self::Value, E> {
        Ok(DuplicateRejectingJson)
    }

    fn visit_str<E>(self, _value: &str) -> Result<Self::Value, E> {
        Ok(DuplicateRejectingJson)
    }

    fn visit_string<E>(self, _value: String) -> Result<Self::Value, E> {
        Ok(DuplicateRejectingJson)
    }

    fn visit_seq<A>(self, mut sequence: A) -> Result<Self::Value, A::Error>
    where
        A: SeqAccess<'de>,
    {
        while sequence.next_element::<DuplicateRejectingJson>()?.is_some() {}
        Ok(DuplicateRejectingJson)
    }

    fn visit_map<A>(self, mut map: A) -> Result<Self::Value, A::Error>
    where
        A: MapAccess<'de>,
    {
        let mut keys = BTreeSet::new();
        while let Some(key) = map.next_key::<String>()? {
            if !keys.insert(key.clone()) {
                return Err(de::Error::custom(format!(
                    "duplicate JSON object key {key:?}"
                )));
            }
            map.next_value::<DuplicateRejectingJson>()?;
        }
        Ok(DuplicateRejectingJson)
    }
}

fn parse_strict_json(bytes: &[u8]) -> Result<Value, String> {
    let mut duplicate_deserializer = serde_json::Deserializer::from_slice(bytes);
    DuplicateRejectingJson::deserialize(&mut duplicate_deserializer)
        .map_err(|error| error.to_string())?;
    duplicate_deserializer
        .end()
        .map_err(|error| error.to_string())?;

    let mut value_deserializer = serde_json::Deserializer::from_slice(bytes);
    let value = Value::deserialize(&mut value_deserializer).map_err(|error| error.to_string())?;
    value_deserializer
        .end()
        .map_err(|error| error.to_string())?;
    reject_non_integer_json_numbers(&value)?;
    Ok(value)
}

fn reject_non_integer_json_numbers(value: &Value) -> Result<(), String> {
    match value {
        Value::Number(number) => {
            let rendered = number.to_string();
            if rendered.contains(['.', 'e', 'E']) {
                Err(format!(
                    "floating-point JSON number {rendered} is unsupported"
                ))
            } else {
                Ok(())
            }
        }
        Value::Array(items) => items.iter().try_for_each(reject_non_integer_json_numbers),
        Value::Object(object) => object
            .values()
            .try_for_each(reject_non_integer_json_numbers),
        Value::Null | Value::Bool(_) | Value::String(_) => Ok(()),
    }
}

#[derive(Debug)]
struct SchemaViolation {
    instance_path: String,
    schema_path: String,
    message: String,
}

const PLAN_SCHEMA_DRAFT: &str = "https://json-schema.org/draft/2020-12/schema";
const MAX_SCHEMA_VIOLATIONS: usize = 256;
const SUPPORTED_SCHEMA_KEYWORDS: [&str; 21] = [
    "$schema",
    "$id",
    "$defs",
    "$ref",
    "title",
    "description",
    "type",
    "const",
    "enum",
    "required",
    "properties",
    "additionalProperties",
    "items",
    "minItems",
    "maxItems",
    "minLength",
    "minimum",
    "uniqueItems",
    "pattern",
    "format",
    "allOf",
];
const SUPPORTED_SCHEMA_CONDITIONAL_KEYWORDS: [&str; 2] = ["if", "then"];

fn validate_supported_plan_schema(schema: &Value) -> Vec<String> {
    let mut problems = Vec::new();
    inspect_supported_schema(schema, schema, "#", true, &mut problems);
    inspect_ref_cycles(schema, &mut problems);
    problems.sort();
    problems.dedup();
    problems
        .into_iter()
        .map(|problem| format!("embedded plan schema is not valid Draft 2020-12: {problem}"))
        .collect()
}

fn inspect_supported_schema(
    schema: &Value,
    root: &Value,
    schema_path: &str,
    is_root: bool,
    problems: &mut Vec<String>,
) {
    let Some(object) = schema.as_object() else {
        problems.push(format!(
            "{schema_path} is a boolean or non-object schema; only object schemas are supported"
        ));
        return;
    };

    for keyword in object.keys() {
        if !SUPPORTED_SCHEMA_KEYWORDS.contains(&keyword.as_str())
            && !SUPPORTED_SCHEMA_CONDITIONAL_KEYWORDS.contains(&keyword.as_str())
        {
            problems.push(format!(
                "{} uses unsupported keyword {keyword:?}",
                schema_pointer(schema_path, keyword)
            ));
        }
    }

    if let Some(draft) = object.get("$schema") {
        if !is_root {
            problems.push(format!(
                "{} is only supported at the schema root",
                schema_pointer(schema_path, "$schema")
            ));
        }
        if draft.as_str() != Some(PLAN_SCHEMA_DRAFT) {
            problems.push(format!(
                "{} must identify Draft 2020-12 exactly",
                schema_pointer(schema_path, "$schema")
            ));
        }
    } else if is_root {
        problems.push("# is missing its required $schema declaration".to_owned());
    }

    if let Some(id) = object.get("$id") {
        if !is_root {
            problems.push(format!(
                "{} is only supported at the schema root",
                schema_pointer(schema_path, "$id")
            ));
        }
        if id.as_str().is_none() {
            problems.push(format!(
                "{} must be a string",
                schema_pointer(schema_path, "$id")
            ));
        }
    }
    for annotation in ["title", "description"] {
        if object
            .get(annotation)
            .is_some_and(|value| !value.is_string())
        {
            problems.push(format!(
                "{} must be a string",
                schema_pointer(schema_path, annotation)
            ));
        }
    }

    if let Some(reference) = object.get("$ref") {
        match reference.as_str() {
            Some(reference) => {
                if resolve_plan_schema_ref(root, reference).is_none() {
                    problems.push(format!(
                        "{} must be an existing local #/$defs/<name> reference; remote, nested, escaped, and missing references are unsupported",
                        schema_pointer(schema_path, "$ref")
                    ));
                }
            }
            None => problems.push(format!(
                "{} must be a string",
                schema_pointer(schema_path, "$ref")
            )),
        }
    }

    if let Some(schema_type) = object.get("type") {
        let mut types = Vec::new();
        match schema_type {
            Value::String(schema_type) => types.push(schema_type.as_str()),
            Value::Array(items) if !items.is_empty() => {
                for item in items {
                    if let Some(schema_type) = item.as_str() {
                        types.push(schema_type);
                    } else {
                        problems.push(format!(
                            "{} must contain only strings",
                            schema_pointer(schema_path, "type")
                        ));
                    }
                }
            }
            _ => problems.push(format!(
                "{} must be a type string or a nonempty array of type strings",
                schema_pointer(schema_path, "type")
            )),
        }
        let allowed = ["object", "array", "string", "integer", "boolean", "null"];
        let mut seen = BTreeSet::new();
        for schema_type in types {
            if !allowed.contains(&schema_type) {
                problems.push(format!(
                    "{} uses unsupported type {schema_type:?}",
                    schema_pointer(schema_path, "type")
                ));
            }
            if !seen.insert(schema_type) {
                problems.push(format!(
                    "{} repeats type {schema_type:?}",
                    schema_pointer(schema_path, "type")
                ));
            }
        }
    }

    if let Some(values) = object.get("enum") {
        match values.as_array() {
            Some(values) if !values.is_empty() => {
                for (index, value) in values.iter().enumerate() {
                    if values[..index].contains(value) {
                        problems.push(format!(
                            "{} contains duplicate values",
                            schema_pointer(schema_path, "enum")
                        ));
                        break;
                    }
                }
            }
            _ => problems.push(format!(
                "{} must be a nonempty array",
                schema_pointer(schema_path, "enum")
            )),
        }
    }

    if let Some(required) = object.get("required") {
        inspect_unique_string_array(required, schema_path, "required", problems);
    }
    if let Some(additional) = object.get("additionalProperties") {
        if additional != &Value::Bool(false) {
            problems.push(format!(
                "{} supports only false",
                schema_pointer(schema_path, "additionalProperties")
            ));
        }
    }
    if let Some(unique) = object.get("uniqueItems") {
        if !unique.is_boolean() {
            problems.push(format!(
                "{} must be a boolean",
                schema_pointer(schema_path, "uniqueItems")
            ));
        }
    }

    let min_items = inspect_nonnegative_schema_integer(
        object.get("minItems"),
        schema_path,
        "minItems",
        problems,
    );
    let max_items = inspect_nonnegative_schema_integer(
        object.get("maxItems"),
        schema_path,
        "maxItems",
        problems,
    );
    if min_items
        .zip(max_items)
        .is_some_and(|(minimum, maximum)| minimum > maximum)
    {
        problems.push(format!("{schema_path} has minItems greater than maxItems"));
    }
    inspect_nonnegative_schema_integer(object.get("minLength"), schema_path, "minLength", problems);
    if let Some(minimum) = object.get("minimum") {
        if !matches!(minimum.as_u64(), Some(0 | 1)) {
            problems.push(format!(
                "{} supports only the exact nonnegative integer bounds 0 and 1",
                schema_pointer(schema_path, "minimum")
            ));
        }
    }

    if let Some(pattern) = object.get("pattern") {
        match pattern.as_str() {
            Some(pattern) if is_supported_plan_pattern(pattern) => {}
            Some(pattern) => problems.push(format!(
                "{} uses unsupported pattern {pattern:?}",
                schema_pointer(schema_path, "pattern")
            )),
            None => problems.push(format!(
                "{} must be a string",
                schema_pointer(schema_path, "pattern")
            )),
        }
    }
    if let Some(format) = object.get("format") {
        match format.as_str() {
            Some("date" | "date-time") => {}
            Some(format) => problems.push(format!(
                "{} uses unsupported format {format:?}",
                schema_pointer(schema_path, "format")
            )),
            None => problems.push(format!(
                "{} must be a string",
                schema_pointer(schema_path, "format")
            )),
        }
    }

    let has_if = object.contains_key("if");
    let has_then = object.contains_key("then");
    if has_if != has_then {
        problems.push(format!(
            "{schema_path} must pair if and then; else is intentionally unsupported"
        ));
    }

    if let Some(definitions) = object.get("$defs") {
        if !is_root {
            problems.push(format!(
                "{} is only supported at the schema root",
                schema_pointer(schema_path, "$defs")
            ));
        }
        match definitions.as_object() {
            Some(definitions) => {
                for (name, definition) in definitions {
                    inspect_supported_schema(
                        definition,
                        root,
                        &schema_pointer(&schema_pointer(schema_path, "$defs"), name),
                        false,
                        problems,
                    );
                }
            }
            None => problems.push(format!(
                "{} must be an object",
                schema_pointer(schema_path, "$defs")
            )),
        }
    }
    if let Some(properties) = object.get("properties") {
        match properties.as_object() {
            Some(properties) => {
                for (name, property_schema) in properties {
                    inspect_supported_schema(
                        property_schema,
                        root,
                        &schema_pointer(&schema_pointer(schema_path, "properties"), name),
                        false,
                        problems,
                    );
                }
            }
            None => problems.push(format!(
                "{} must be an object",
                schema_pointer(schema_path, "properties")
            )),
        }
    }
    for keyword in ["items", "if", "then"] {
        if let Some(subschema) = object.get(keyword) {
            inspect_supported_schema(
                subschema,
                root,
                &schema_pointer(schema_path, keyword),
                false,
                problems,
            );
        }
    }
    if let Some(all_of) = object.get("allOf") {
        match all_of.as_array() {
            Some(all_of) if !all_of.is_empty() => {
                for (index, subschema) in all_of.iter().enumerate() {
                    inspect_supported_schema(
                        subschema,
                        root,
                        &schema_pointer(&schema_pointer(schema_path, "allOf"), &index.to_string()),
                        false,
                        problems,
                    );
                }
            }
            _ => problems.push(format!(
                "{} must be a nonempty array",
                schema_pointer(schema_path, "allOf")
            )),
        }
    }
}

fn inspect_unique_string_array(
    value: &Value,
    schema_path: &str,
    keyword: &str,
    problems: &mut Vec<String>,
) {
    let Some(items) = value.as_array() else {
        problems.push(format!(
            "{} must be an array of unique strings",
            schema_pointer(schema_path, keyword)
        ));
        return;
    };
    let mut seen = BTreeSet::new();
    for item in items {
        let Some(item) = item.as_str() else {
            problems.push(format!(
                "{} must contain only strings",
                schema_pointer(schema_path, keyword)
            ));
            continue;
        };
        if !seen.insert(item) {
            problems.push(format!(
                "{} repeats {item:?}",
                schema_pointer(schema_path, keyword)
            ));
        }
    }
}

fn inspect_nonnegative_schema_integer(
    value: Option<&Value>,
    schema_path: &str,
    keyword: &str,
    problems: &mut Vec<String>,
) -> Option<u64> {
    let value = value?;
    match value.as_u64() {
        Some(value) => Some(value),
        None => {
            problems.push(format!(
                "{} must be a nonnegative integer",
                schema_pointer(schema_path, keyword)
            ));
            None
        }
    }
}

fn plan_schema_ref_name(reference: &str) -> Option<&str> {
    let name = reference.strip_prefix("#/$defs/")?;
    if name.is_empty() || name.contains(['/', '~', '#']) {
        return None;
    }
    Some(name)
}

fn resolve_plan_schema_ref<'a>(root: &'a Value, reference: &str) -> Option<&'a Value> {
    let name = plan_schema_ref_name(reference)?;
    root.get("$defs")?.as_object()?.get(name)
}

fn inspect_ref_cycles(root: &Value, problems: &mut Vec<String>) {
    let Some(definitions) = root.get("$defs").and_then(Value::as_object) else {
        return;
    };
    let mut states: BTreeMap<String, u8> = BTreeMap::new();
    let mut stack = Vec::new();
    for name in definitions.keys() {
        inspect_ref_cycle_from(name, root, &mut states, &mut stack, problems);
    }
}

fn inspect_ref_cycle_from(
    name: &str,
    root: &Value,
    states: &mut BTreeMap<String, u8>,
    stack: &mut Vec<String>,
    problems: &mut Vec<String>,
) {
    match states.get(name).copied() {
        Some(2) => return,
        Some(1) => {
            let mut cycle = stack.clone();
            cycle.push(name.to_owned());
            problems.push(format!(
                "#/$defs/{name} participates in unsupported recursive reference cycle {}",
                cycle.join(" -> ")
            ));
            return;
        }
        _ => {}
    }
    let Some(definition) = root
        .get("$defs")
        .and_then(Value::as_object)
        .and_then(|definitions| definitions.get(name))
    else {
        return;
    };

    states.insert(name.to_owned(), 1);
    stack.push(name.to_owned());
    let mut references = Vec::new();
    collect_schema_refs(definition, &mut references);
    for reference in references {
        if let Some(target) = plan_schema_ref_name(&reference) {
            inspect_ref_cycle_from(target, root, states, stack, problems);
        }
    }
    stack.pop();
    states.insert(name.to_owned(), 2);
}

fn collect_schema_refs(schema: &Value, references: &mut Vec<String>) {
    let Some(object) = schema.as_object() else {
        return;
    };
    if let Some(reference) = object.get("$ref").and_then(Value::as_str) {
        references.push(reference.to_owned());
    }
    if let Some(properties) = object.get("properties").and_then(Value::as_object) {
        for subschema in properties.values() {
            collect_schema_refs(subschema, references);
        }
    }
    for keyword in ["items", "if", "then"] {
        if let Some(subschema) = object.get(keyword) {
            collect_schema_refs(subschema, references);
        }
    }
    if let Some(all_of) = object.get("allOf").and_then(Value::as_array) {
        for subschema in all_of {
            collect_schema_refs(subschema, references);
        }
    }
}

fn schema_pointer(base: &str, segment: &str) -> String {
    format!("{base}/{}", json_pointer_segment(segment))
}

fn instance_pointer(base: &str, segment: &str) -> String {
    format!("{base}/{}", json_pointer_segment(segment))
}

fn json_pointer_segment(segment: &str) -> String {
    segment.replace('~', "~0").replace('/', "~1")
}

fn validate_schema_instance(
    instance: &Value,
    schema: &Value,
    root: &Value,
    instance_path: &str,
    schema_path: &str,
    ref_stack: &mut Vec<String>,
    violations: &mut Vec<SchemaViolation>,
) {
    if violations.len() >= MAX_SCHEMA_VIOLATIONS {
        return;
    }
    let object = schema
        .as_object()
        .expect("supported schema inspection guarantees object schemas");

    if let Some(reference) = object.get("$ref").and_then(Value::as_str) {
        let name = plan_schema_ref_name(reference)
            .expect("supported schema inspection guarantees a local definition reference");
        if ref_stack.iter().any(|active| active == name) {
            push_schema_violation(
                violations,
                instance_path,
                &schema_pointer(schema_path, "$ref"),
                "recursive schema references are unsupported",
            );
        } else {
            let target = resolve_plan_schema_ref(root, reference)
                .expect("supported schema inspection guarantees the definition exists");
            ref_stack.push(name.to_owned());
            validate_schema_instance(
                instance,
                target,
                root,
                instance_path,
                &schema_pointer("#/$defs", name),
                ref_stack,
                violations,
            );
            ref_stack.pop();
        }
    }

    if let Some(schema_type) = object.get("type") {
        let matches = match schema_type {
            Value::String(schema_type) => instance_has_schema_type(instance, schema_type),
            Value::Array(types) => types.iter().any(|schema_type| {
                instance_has_schema_type(
                    instance,
                    schema_type
                        .as_str()
                        .expect("supported schema inspection guarantees string types"),
                )
            }),
            _ => unreachable!("supported schema inspection guarantees the type shape"),
        };
        if !matches {
            push_schema_violation(
                violations,
                instance_path,
                &schema_pointer(schema_path, "type"),
                &format!("value does not match declared type {schema_type}"),
            );
        }
    }

    if let Some(expected) = object.get("const") {
        if instance != expected {
            push_schema_violation(
                violations,
                instance_path,
                &schema_pointer(schema_path, "const"),
                "value does not match const",
            );
        }
    }
    if let Some(choices) = object.get("enum").and_then(Value::as_array) {
        if !choices.contains(instance) {
            push_schema_violation(
                violations,
                instance_path,
                &schema_pointer(schema_path, "enum"),
                "value is not one of the enumerated choices",
            );
        }
    }

    if let Some(instance_object) = instance.as_object() {
        if let Some(required) = object.get("required").and_then(Value::as_array) {
            for property in required {
                let property = property
                    .as_str()
                    .expect("supported schema inspection guarantees string required names");
                if !instance_object.contains_key(property) {
                    push_schema_violation(
                        violations,
                        instance_path,
                        &schema_pointer(schema_path, "required"),
                        &format!("required property {property:?} is missing"),
                    );
                }
            }
        }

        let properties = object.get("properties").and_then(Value::as_object);
        if let Some(properties) = properties {
            for (property, property_schema) in properties {
                if let Some(property_value) = instance_object.get(property) {
                    validate_schema_instance(
                        property_value,
                        property_schema,
                        root,
                        &instance_pointer(instance_path, property),
                        &schema_pointer(&schema_pointer(schema_path, "properties"), property),
                        ref_stack,
                        violations,
                    );
                }
            }
        }
        if object.get("additionalProperties") == Some(&Value::Bool(false)) {
            for property in instance_object.keys() {
                if !properties.is_some_and(|properties| properties.contains_key(property)) {
                    push_schema_violation(
                        violations,
                        &instance_pointer(instance_path, property),
                        &schema_pointer(schema_path, "additionalProperties"),
                        &format!("additional property {property:?} is not allowed"),
                    );
                }
            }
        }
    }

    if let Some(items) = instance.as_array() {
        if let Some(item_schema) = object.get("items") {
            for (index, item) in items.iter().enumerate() {
                if violations.len() >= MAX_SCHEMA_VIOLATIONS {
                    break;
                }
                validate_schema_instance(
                    item,
                    item_schema,
                    root,
                    &instance_pointer(instance_path, &index.to_string()),
                    &schema_pointer(schema_path, "items"),
                    ref_stack,
                    violations,
                );
            }
        }
        if let Some(minimum) = object.get("minItems").and_then(Value::as_u64) {
            if (items.len() as u64) < minimum {
                push_schema_violation(
                    violations,
                    instance_path,
                    &schema_pointer(schema_path, "minItems"),
                    &format!("array has fewer than {minimum} items"),
                );
            }
        }
        if let Some(maximum) = object.get("maxItems").and_then(Value::as_u64) {
            if (items.len() as u64) > maximum {
                push_schema_violation(
                    violations,
                    instance_path,
                    &schema_pointer(schema_path, "maxItems"),
                    &format!("array has more than {maximum} items"),
                );
            }
        }
        if object.get("uniqueItems") == Some(&Value::Bool(true)) {
            let mut seen = BTreeSet::new();
            for (index, item) in items.iter().enumerate() {
                if violations.len() >= MAX_SCHEMA_VIOLATIONS {
                    break;
                }
                let encoded = serde_json::to_string(item)
                    .expect("a parsed serde_json::Value always serializes");
                if !seen.insert(encoded) {
                    push_schema_violation(
                        violations,
                        instance_path,
                        &schema_pointer(schema_path, "uniqueItems"),
                        &format!("array item {index} duplicates an earlier item"),
                    );
                }
            }
        }
    }

    if let Some(text) = instance.as_str() {
        if let Some(minimum) = object.get("minLength").and_then(Value::as_u64) {
            if (text.chars().count() as u64) < minimum {
                push_schema_violation(
                    violations,
                    instance_path,
                    &schema_pointer(schema_path, "minLength"),
                    &format!("string has fewer than {minimum} Unicode scalar values"),
                );
            }
        }
        if let Some(pattern) = object.get("pattern").and_then(Value::as_str) {
            if !matches_plan_pattern(pattern, text) {
                push_schema_violation(
                    violations,
                    instance_path,
                    &schema_pointer(schema_path, "pattern"),
                    &format!("string does not match supported pattern {pattern:?}"),
                );
            }
        }
        if let Some(format) = object.get("format").and_then(Value::as_str) {
            let valid = match format {
                "date" => is_rfc3339_full_date(text),
                "date-time" => is_rfc3339_date_time(text),
                _ => unreachable!("supported schema inspection guarantees the format"),
            };
            if !valid {
                push_schema_violation(
                    violations,
                    instance_path,
                    &schema_pointer(schema_path, "format"),
                    &format!("string is not a valid {format}"),
                );
            }
        }
    }

    if let Some(number) = instance.as_number() {
        if let Some(minimum) = object.get("minimum").and_then(Value::as_u64) {
            if !json_number_at_least_supported_minimum(number, minimum) {
                push_schema_violation(
                    violations,
                    instance_path,
                    &schema_pointer(schema_path, "minimum"),
                    &format!("number is less than {minimum}"),
                );
            }
        }
    }

    if let Some(all_of) = object.get("allOf").and_then(Value::as_array) {
        for (index, subschema) in all_of.iter().enumerate() {
            validate_schema_instance(
                instance,
                subschema,
                root,
                instance_path,
                &schema_pointer(&schema_pointer(schema_path, "allOf"), &index.to_string()),
                ref_stack,
                violations,
            );
        }
    }

    if let (Some(condition), Some(consequence)) = (object.get("if"), object.get("then")) {
        let mut condition_violations = Vec::new();
        validate_schema_instance(
            instance,
            condition,
            root,
            instance_path,
            &schema_pointer(schema_path, "if"),
            ref_stack,
            &mut condition_violations,
        );
        if condition_violations.is_empty() {
            validate_schema_instance(
                instance,
                consequence,
                root,
                instance_path,
                &schema_pointer(schema_path, "then"),
                ref_stack,
                violations,
            );
        }
    }
}

fn push_schema_violation(
    violations: &mut Vec<SchemaViolation>,
    instance_path: &str,
    schema_path: &str,
    message: &str,
) {
    if violations.len() < MAX_SCHEMA_VIOLATIONS {
        violations.push(SchemaViolation {
            instance_path: instance_path.to_owned(),
            schema_path: schema_path.to_owned(),
            message: message.to_owned(),
        });
    }
}

fn instance_has_schema_type(instance: &Value, schema_type: &str) -> bool {
    match schema_type {
        "object" => instance.is_object(),
        "array" => instance.is_array(),
        "string" => instance.is_string(),
        "integer" => instance
            .as_number()
            .is_some_and(json_number_is_mathematical_integer),
        "boolean" => instance.is_boolean(),
        "null" => instance.is_null(),
        _ => false,
    }
}

const SUPPORTED_PLAN_PATTERNS: [&str; 13] = [
    r"^urn:ecdlp-audit-plan:1:sha256:[0-9a-f]{64}$",
    r"^[0-9a-f]{64}$",
    r"^source\.[a-z0-9._-]+$",
    r"^[0-9a-f]{40}$",
    r"^WS[1-8]-D[0-9]+$",
    r"^[a-z][a-z0-9_]*$",
    r"^WS[1-8]-G[0-9]+$",
    r"^[a-z][a-z0-9._-]+$",
    r"^WS[1-8]$",
    r"^[a-z0-9]+(?:-[a-z0-9]+)*$",
    r"^route\.[a-e]\.[a-z0-9._-]+$",
    r"^urn:ecdlp-audit-unit:1:sha256:[0-9a-f]{64}$",
    r"^urn:ecdlp-audit-assessment:1:sha256:[0-9a-f]{64}$",
];

fn is_supported_plan_pattern(pattern: &str) -> bool {
    SUPPORTED_PLAN_PATTERNS.contains(&pattern)
}

fn matches_plan_pattern(pattern: &str, value: &str) -> bool {
    match pattern {
        r"^urn:ecdlp-audit-plan:1:sha256:[0-9a-f]{64}$" => value
            .strip_prefix("urn:ecdlp-audit-plan:1:sha256:")
            .is_some_and(|suffix| is_lower_hex(suffix, 64)),
        r"^[0-9a-f]{64}$" => is_lower_hex(value, 64),
        r"^source\.[a-z0-9._-]+$" => value
            .strip_prefix("source.")
            .is_some_and(is_nonempty_lower_identifier_tail),
        r"^[0-9a-f]{40}$" => is_lower_hex(value, 40),
        r"^WS[1-8]-D[0-9]+$" => matches_ws_numbered_id(value, 'D'),
        r"^[a-z][a-z0-9_]*$" => {
            let bytes = value.as_bytes();
            bytes.first().is_some_and(u8::is_ascii_lowercase)
                && bytes[1..]
                    .iter()
                    .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || *byte == b'_')
        }
        r"^WS[1-8]-G[0-9]+$" => matches_ws_numbered_id(value, 'G'),
        r"^[a-z][a-z0-9._-]+$" => {
            let bytes = value.as_bytes();
            bytes.len() >= 2
                && bytes.first().is_some_and(u8::is_ascii_lowercase)
                && bytes[1..].iter().all(|byte| {
                    byte.is_ascii_lowercase()
                        || byte.is_ascii_digit()
                        || matches!(*byte, b'.' | b'_' | b'-')
                })
        }
        r"^WS[1-8]$" => matches_ws_id(value),
        r"^[a-z0-9]+(?:-[a-z0-9]+)*$" => {
            !value.is_empty()
                && value.split('-').all(|segment| {
                    !segment.is_empty()
                        && segment
                            .bytes()
                            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit())
                })
        }
        r"^route\.[a-e]\.[a-z0-9._-]+$" => value
            .strip_prefix("route.")
            .and_then(|tail| {
                let bytes = tail.as_bytes();
                (bytes.len() >= 3 && (b'a'..=b'e').contains(&bytes[0]) && bytes[1] == b'.')
                    .then_some(&tail[2..])
            })
            .is_some_and(is_nonempty_lower_identifier_tail),
        r"^urn:ecdlp-audit-unit:1:sha256:[0-9a-f]{64}$" => value
            .strip_prefix("urn:ecdlp-audit-unit:1:sha256:")
            .is_some_and(|suffix| is_lower_hex(suffix, 64)),
        r"^urn:ecdlp-audit-assessment:1:sha256:[0-9a-f]{64}$" => value
            .strip_prefix("urn:ecdlp-audit-assessment:1:sha256:")
            .is_some_and(|suffix| is_lower_hex(suffix, 64)),
        _ => false,
    }
}

fn is_nonempty_lower_identifier_tail(value: &str) -> bool {
    !value.is_empty()
        && value.bytes().all(|byte| {
            byte.is_ascii_lowercase() || byte.is_ascii_digit() || matches!(byte, b'.' | b'_' | b'-')
        })
}

fn matches_ws_id(value: &str) -> bool {
    let bytes = value.as_bytes();
    bytes.len() == 3 && bytes[0..2] == *b"WS" && (b'1'..=b'8').contains(&bytes[2])
}

fn matches_ws_numbered_id(value: &str, kind: char) -> bool {
    let bytes = value.as_bytes();
    bytes.len() >= 6
        && bytes[0..2] == *b"WS"
        && (b'1'..=b'8').contains(&bytes[2])
        && bytes[3] == b'-'
        && bytes[4] == kind as u8
        && bytes[5..].iter().all(u8::is_ascii_digit)
}

fn decimal_number_parts(number: &serde_json::Number) -> (bool, String, usize, i64) {
    let rendered = number.to_string();
    let (negative, unsigned) = match rendered.strip_prefix('-') {
        Some(unsigned) => (true, unsigned),
        None => (false, rendered.as_str()),
    };
    let (mantissa, exponent_text) = unsigned
        .split_once(['e', 'E'])
        .map_or((unsigned, None), |(mantissa, exponent)| {
            (mantissa, Some(exponent))
        });
    let exponent = exponent_text.map_or(0, |exponent| {
        exponent.parse::<i64>().unwrap_or_else(|_| {
            if exponent.starts_with('-') {
                i64::MIN / 4
            } else {
                i64::MAX / 4
            }
        })
    });
    let (integer, fraction) = mantissa
        .split_once('.')
        .map_or((mantissa, ""), |(integer, fraction)| (integer, fraction));
    let mut digits = String::with_capacity(integer.len() + fraction.len());
    digits.push_str(integer);
    digits.push_str(fraction);
    (negative, digits, fraction.len(), exponent)
}

fn json_number_is_mathematical_integer(number: &serde_json::Number) -> bool {
    let (_, digits, fractional_digits, exponent) = decimal_number_parts(number);
    let fractional_places = fractional_digits as i128 - exponent as i128;
    if fractional_places <= 0 {
        return true;
    }
    let Ok(fractional_places) = usize::try_from(fractional_places) else {
        return digits.bytes().all(|digit| digit == b'0');
    };
    if fractional_places >= digits.len() {
        digits.bytes().all(|digit| digit == b'0')
    } else {
        digits.as_bytes()[digits.len() - fractional_places..]
            .iter()
            .all(|digit| *digit == b'0')
    }
}

fn json_number_at_least_supported_minimum(number: &serde_json::Number, minimum: u64) -> bool {
    let (negative, digits, fractional_digits, exponent) = decimal_number_parts(number);
    let nonzero = digits.bytes().any(|digit| digit != b'0');
    if negative && nonzero {
        return false;
    }
    match minimum {
        0 => true,
        1 => {
            if !nonzero {
                return false;
            }
            let significant_digits = digits.trim_start_matches('0').len() as i128;
            significant_digits + exponent as i128 - fractional_digits as i128 >= 1
        }
        _ => false,
    }
}

fn is_rfc3339_full_date(value: &str) -> bool {
    is_rfc3339_full_date_bytes(value.as_bytes())
}

fn is_rfc3339_full_date_bytes(bytes: &[u8]) -> bool {
    if bytes.len() != 10 || bytes[4] != b'-' || bytes[7] != b'-' {
        return false;
    }
    let Some(year) = parse_ascii_decimal(&bytes[0..4]) else {
        return false;
    };
    let Some(month) = parse_ascii_decimal(&bytes[5..7]) else {
        return false;
    };
    let Some(day) = parse_ascii_decimal(&bytes[8..10]) else {
        return false;
    };
    let Some(days) = days_in_month(year, month) else {
        return false;
    };
    (1..=days).contains(&day)
}

fn is_rfc3339_date_time(value: &str) -> bool {
    let bytes = value.as_bytes();
    if bytes.len() < 20
        || !is_rfc3339_full_date_bytes(&bytes[..10])
        || !matches!(bytes[10], b'T' | b't')
        || bytes[13] != b':'
        || bytes[16] != b':'
    {
        return false;
    }
    let Some(hour) = parse_ascii_decimal(&bytes[11..13]) else {
        return false;
    };
    let Some(minute) = parse_ascii_decimal(&bytes[14..16]) else {
        return false;
    };
    let Some(second) = parse_ascii_decimal(&bytes[17..19]) else {
        return false;
    };
    if hour > 23 || minute > 59 || second > 60 {
        return false;
    }

    let mut cursor = 19;
    if bytes.get(cursor) == Some(&b'.') {
        cursor += 1;
        let fraction_start = cursor;
        while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
            cursor += 1;
        }
        if cursor == fraction_start {
            return false;
        }
    }
    let offset_minutes = match bytes.get(cursor) {
        Some(b'Z' | b'z') if cursor + 1 == bytes.len() => 0,
        Some(b'+' | b'-') => {
            if cursor + 6 != bytes.len() || bytes[cursor + 3] != b':' {
                return false;
            }
            let Some(offset_hour) = parse_ascii_decimal(&bytes[cursor + 1..cursor + 3]) else {
                return false;
            };
            let Some(offset_minute) = parse_ascii_decimal(&bytes[cursor + 4..cursor + 6]) else {
                return false;
            };
            if offset_hour > 23 || offset_minute > 59 {
                return false;
            }
            let magnitude = (offset_hour * 60 + offset_minute) as i32;
            if bytes[cursor] == b'+' {
                magnitude
            } else {
                -magnitude
            }
        }
        _ => return false,
    };
    if second == 60 {
        // RFC 3339 permits :60 only at the UTC end of a day. This is the
        // deterministic format-level rule; it does not claim an IERS table.
        ((hour * 60 + minute) as i32 - offset_minutes).rem_euclid(24 * 60) == 23 * 60 + 59
    } else {
        true
    }
}

fn days_in_month(year: u32, month: u32) -> Option<u32> {
    match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => Some(31),
        4 | 6 | 9 | 11 => Some(30),
        2 if is_leap_year(year) => Some(29),
        2 => Some(28),
        _ => None,
    }
}

fn parse_ascii_decimal(bytes: &[u8]) -> Option<u32> {
    if bytes.is_empty() || !bytes.iter().all(u8::is_ascii_digit) {
        return None;
    }
    bytes.iter().try_fold(0u32, |value, digit| {
        value.checked_mul(10)?.checked_add(u32::from(*digit - b'0'))
    })
}

fn is_leap_year(year: u32) -> bool {
    year.is_multiple_of(4) && (!year.is_multiple_of(100) || year.is_multiple_of(400))
}

fn read_bounded(path: &Path, limit: u64, kind: &str) -> Result<Vec<u8>, String> {
    let metadata = fs::metadata(path)
        .map_err(|error| format!("cannot inspect {kind} {}: {error}", path.display()))?;
    if metadata.len() > limit {
        return Err(format!(
            "{kind} exceeds the {limit}-byte limit: {}",
            path.display()
        ));
    }
    fs::read(path).map_err(|error| format!("cannot read {kind} {}: {error}", path.display()))
}

fn validate_plan(
    plan: &Plan,
    raw_plan: &Value,
    repo_root: &Path,
    input_path: &Path,
) -> Vec<String> {
    let mut errors = Vec::new();
    require(
        plan.schema_version == EXPECTED_SCHEMA_VERSION,
        "unexpected schema_version",
        &mut errors,
    );
    require(
        plan.record_state == EXPECTED_RECORD_STATE,
        "record must remain a non-authoritative mock",
        &mut errors,
    );
    require(
        !plan.authoritative,
        "mock cannot be authoritative",
        &mut errors,
    );
    require(
        !plan.promotion_allowed,
        "mock cannot allow promotion",
        &mut errors,
    );
    require(
        !plan.scope.live_target_execution
            && !plan.scope.secret_key_material
            && !plan.scope.production_exploitation,
        "planning scope cannot authorize live targets, secrets, or production exploitation",
        &mut errors,
    );
    require(
        plan.scope.claim_boundary == EXPECTED_SCOPE_CLAIM_BOUNDARY,
        "planning scope claim boundary drifted",
        &mut errors,
    );
    require(
        strings_match(
            &plan.scope.execution_fixture_classes,
            &["passive_public", "public_synthetic", "known_answer"],
        ) && strings_match(&plan.scope.default_authorization_tiers, &["T0", "T1"]),
        "planning mock execution fixtures or default authorization tiers drifted",
        &mut errors,
    );
    require(
        strings_match(
            &plan.scope.outcome_classes,
            &[
                "exact_scalar_recovery",
                "equivalent_capability",
                "partial_secret_leakage",
                "property_or_provenance_finding",
                "no_recovery_within_boundary",
            ],
        ),
        "planning mock outcome classes drifted",
        &mut errors,
    );
    require(
        strings_match(
            &plan.vocabularies.claim_status,
            &["unknown", "supported", "refuted", "not_applicable"],
        ) && strings_match(
            &plan.vocabularies.requirement_status,
            &[
                "verified_complete",
                "implemented_unverified",
                "partial",
                "blocked",
                "not_attempted",
            ],
        ) && strings_match(
            &plan.vocabularies.maturity_levels,
            &["E0", "E1", "E2", "E3", "E4", "E5", "E6", "E7"],
        ) && strings_match(
            &plan.vocabularies.independence_levels,
            &["R0", "R1", "R2", "R3", "R4", "R5"],
        ) && strings_match(
            &plan.vocabularies.authorization_tiers,
            &["T0", "T1", "T2", "T3", "T4"],
        ) && strings_match(
            &plan.vocabularies.cost_status,
            &[
                "measured",
                "modeled",
                "extrapolated",
                "unknown",
                "not_applicable",
            ],
        ) && strings_match(
            &plan.vocabularies.execution_outcome,
            &[
                "complete",
                "insufficient_relations",
                "budget",
                "timeout",
                "oom",
                "error",
            ],
        ) && strings_match(
            &plan.vocabularies.review_decision,
            &["approved", "changes_requested", "pending", "conflicted"],
        ),
        "planning mock vocabularies drifted from the versioned contract",
        &mut errors,
    );
    require(
        plan.governance.status == "proposed",
        "planning mock governance must remain proposed",
        &mut errors,
    );
    require(
        !plan.governance.self_promotion_allowed,
        "planning mock cannot allow self-promotion",
        &mut errors,
    );
    require(
        plan.governance.immutable_evidence,
        "planning mock governance must require immutable evidence",
        &mut errors,
    );
    require(
        strings_match(
            &plan.governance.required_roles,
            &[
                "author",
                "mathematical_reviewer",
                "implementation_reviewer",
                "evidence_replay_reviewer",
                "release_approver",
            ],
        ),
        "planning mock governance roles drifted",
        &mut errors,
    );

    let workstream_ids: Vec<&str> = plan
        .workstreams
        .iter()
        .map(|item| item.id.as_str())
        .collect();
    require(
        workstream_ids == EXPECTED_WORKSTREAMS,
        "workstreams must be ordered exactly as WS1 through WS8",
        &mut errors,
    );
    require(
        unique(workstream_ids.iter().copied()),
        "workstream ids must be unique",
        &mut errors,
    );

    match expected_plan_uid(raw_plan) {
        Ok(expected) => {
            if plan.plan_uid != expected {
                errors.push(format!(
                    "plan_uid does not match the canonical full-plan preimage; expected {expected}"
                ));
            }
        }
        Err(error) => errors.push(format!("cannot canonicalize plan_uid preimage: {error}")),
    }

    let obligation_by_id: BTreeMap<&str, &Obligation> = plan
        .open_obligations
        .iter()
        .map(|obligation| (obligation.id.as_str(), obligation))
        .collect();
    require(
        obligation_by_id.len() == plan.open_obligations.len(),
        "obligation ids must be unique",
        &mut errors,
    );

    let mut all_gate_ids = BTreeSet::new();
    let mut all_deliverable_ids = BTreeSet::new();
    let mut required_gate_total = 0usize;
    let mut open_required_gate_total = 0usize;
    let mut required_gate_obligation_ids = BTreeSet::new();

    for workstream in &plan.workstreams {
        require(
            workstream.requirement_status != RequirementStatus::VerifiedComplete,
            &format!(
                "{} cannot be verified complete while the planning mock has no assessments",
                workstream.id
            ),
            &mut errors,
        );
        let derived_status = derive_workstream_status(workstream);
        require(
            workstream.requirement_status == derived_status,
            &format!(
                "{} status {:?} does not match deliverable-derived status {:?}",
                workstream.id, workstream.requirement_status, derived_status
            ),
            &mut errors,
        );
        let local_gate_ids: BTreeSet<&str> = workstream
            .acceptance_gates
            .iter()
            .map(|gate| gate.id.as_str())
            .collect();
        require(
            local_gate_ids.len() == workstream.acceptance_gates.len(),
            &format!("{} gate ids must be unique", workstream.id),
            &mut errors,
        );
        for gate_id in &local_gate_ids {
            require(
                all_gate_ids.insert((*gate_id).to_owned()),
                &format!("{} reuses gate id {gate_id}", workstream.id),
                &mut errors,
            );
        }

        let local_deliverable_ids: BTreeSet<&str> = workstream
            .deliverables
            .iter()
            .map(|deliverable| deliverable.id.as_str())
            .collect();
        require(
            local_deliverable_ids.len() == workstream.deliverables.len(),
            &format!("{} deliverable ids must be unique", workstream.id),
            &mut errors,
        );
        for deliverable_id in &local_deliverable_ids {
            require(
                all_deliverable_ids.insert((*deliverable_id).to_owned()),
                &format!("{} reuses deliverable id {deliverable_id}", workstream.id),
                &mut errors,
            );
        }

        let mut referenced_gate_ids = BTreeSet::new();
        for deliverable in &workstream.deliverables {
            require(
                deliverable.status != RequirementStatus::VerifiedComplete,
                &format!(
                    "{} cannot be verified complete while its gates remain unknown",
                    deliverable.id
                ),
                &mut errors,
            );
            require(
                deliverable.id.starts_with(&format!("{}-D", workstream.id)),
                &format!("{} is not namespaced to {}", deliverable.id, workstream.id),
                &mut errors,
            );
            let expected_binding = EXPECTED_DELIVERABLE_BINDINGS
                .iter()
                .find(|(id, _, _)| *id == deliverable.id);
            match expected_binding {
                Some((_, expected_kind, expected_path)) => {
                    require(
                        deliverable.kind == *expected_kind,
                        &format!("{} artifact kind drifted", deliverable.id),
                        &mut errors,
                    );
                    require(
                        deliverable.path.as_deref() == *expected_path,
                        &format!("{} artifact path drifted", deliverable.id),
                        &mut errors,
                    );
                }
                None => errors.push(format!(
                    "{} is not an expected planning-mock deliverable",
                    deliverable.id
                )),
            }
            validate_deliverable_artifact(repo_root, deliverable, &mut errors);
            if matches!(
                deliverable.status,
                RequirementStatus::Partial
                    | RequirementStatus::ImplementedUnverified
                    | RequirementStatus::VerifiedComplete
            ) {
                require(
                    deliverable.path.is_some(),
                    &format!(
                        "{} has status {:?} but no artifact path",
                        deliverable.id, deliverable.status
                    ),
                    &mut errors,
                );
                require(
                    deliverable.artifact_sha256.is_some(),
                    &format!(
                        "{} has status {:?} but no artifact SHA-256",
                        deliverable.id, deliverable.status
                    ),
                    &mut errors,
                );
            }
            if deliverable.status == RequirementStatus::NotAttempted {
                require(
                    deliverable.path.is_none() && deliverable.artifact_sha256.is_none(),
                    &format!("{} is not attempted but binds an artifact", deliverable.id),
                    &mut errors,
                );
            }
            require(
                !deliverable.acceptance_gate_ids.is_empty(),
                &format!("{} has no acceptance gates", deliverable.id),
                &mut errors,
            );
            require(
                unique(deliverable.acceptance_gate_ids.iter().map(String::as_str)),
                &format!("{} repeats an acceptance gate", deliverable.id),
                &mut errors,
            );
            for gate_id in &deliverable.acceptance_gate_ids {
                referenced_gate_ids.insert(gate_id.as_str());
                require(
                    local_gate_ids.contains(gate_id.as_str()),
                    &format!("{} references non-local gate {gate_id}", deliverable.id),
                    &mut errors,
                );
            }
        }

        for gate in workstream
            .acceptance_gates
            .iter()
            .filter(|gate| gate.required)
        {
            require(
                referenced_gate_ids.contains(gate.id.as_str()),
                &format!(
                    "{} required gate {} is not referenced by a deliverable",
                    workstream.id, gate.id
                ),
                &mut errors,
            );
        }

        for gate in &workstream.acceptance_gates {
            require(
                gate.claim_status == ClaimStatus::Unknown,
                &format!("{} must remain unknown in the empty planning mock", gate.id),
                &mut errors,
            );
            require(
                gate.id.starts_with(&format!("{}-G", workstream.id)),
                &format!("{} is not namespaced to {}", gate.id, workstream.id),
                &mut errors,
            );
            if matches!(
                gate.claim_status,
                ClaimStatus::Supported | ClaimStatus::Refuted
            ) {
                require(
                    !gate.evidence_refs.is_empty(),
                    &format!("{} promoted status lacks evidence", gate.id),
                    &mut errors,
                );
            }
            if gate.claim_status == ClaimStatus::Unknown {
                require(
                    gate.obligation_id.is_some(),
                    &format!("{} unknown status lacks an obligation", gate.id),
                    &mut errors,
                );
                require(
                    gate.evidence_refs.is_empty(),
                    &format!("{} unknown status must not cite promoted evidence", gate.id),
                    &mut errors,
                );
            }
            if let Some(obligation_id) = gate.obligation_id.as_deref() {
                if gate.required {
                    require(
                        required_gate_obligation_ids.insert(obligation_id.to_owned()),
                        &format!(
                            "{} reuses required-gate obligation {obligation_id}",
                            gate.id
                        ),
                        &mut errors,
                    );
                }
                require(
                    obligation_by_id.contains_key(obligation_id),
                    &format!("{} references missing obligation {obligation_id}", gate.id),
                    &mut errors,
                );
                if let Some(obligation) = obligation_by_id.get(obligation_id) {
                    require(
                        obligation.workstream_id == workstream.id,
                        &format!("{} obligation belongs to another workstream", gate.id),
                        &mut errors,
                    );
                    if gate.required {
                        require(
                            obligation.required,
                            &format!(
                                "{} required gate references optional obligation {obligation_id}",
                                gate.id
                            ),
                            &mut errors,
                        );
                        require(
                            matches!(
                                obligation.status,
                                ObligationStatus::Open | ObligationStatus::Blocked
                            ),
                            &format!(
                                "{} required gate references closed obligation {obligation_id}",
                                gate.id
                            ),
                            &mut errors,
                        );
                    }
                }
            }
        }

        let required = workstream
            .acceptance_gates
            .iter()
            .filter(|gate| gate.required)
            .count();
        let supported = workstream
            .acceptance_gates
            .iter()
            .filter(|gate| gate.claim_status == ClaimStatus::Supported)
            .count();
        let open_required = workstream
            .acceptance_gates
            .iter()
            .filter(|gate| gate.required && gate.claim_status != ClaimStatus::Supported)
            .count();
        let expected_metrics = WorkstreamMetrics {
            required_gates: required,
            classified_gates: workstream.acceptance_gates.len(),
            supported_gates: supported,
            open_required_gates: open_required,
        };
        require(
            workstream.metrics == expected_metrics,
            &format!("{} metrics do not match its gates", workstream.id),
            &mut errors,
        );
        required_gate_total += required;
        open_required_gate_total += open_required;
    }

    let expected_deliverable_ids: BTreeSet<String> = EXPECTED_DELIVERABLE_BINDINGS
        .iter()
        .map(|(id, _, _)| (*id).to_owned())
        .collect();
    require(
        all_deliverable_ids == expected_deliverable_ids,
        "planning-mock deliverable inventory drifted",
        &mut errors,
    );
    require(
        required_gate_total == 20,
        &format!("expected 20 required planning gates, got {required_gate_total}"),
        &mut errors,
    );
    let required_open_obligation_ids: BTreeSet<String> = plan
        .open_obligations
        .iter()
        .filter(|obligation| {
            obligation.required
                && matches!(
                    obligation.status,
                    ObligationStatus::Open | ObligationStatus::Blocked
                )
        })
        .map(|obligation| obligation.id.clone())
        .collect();
    require(
        plan.completion_metrics.obligations.open_required == required_open_obligation_ids.len(),
        "global open-required count does not match required unresolved obligations",
        &mut errors,
    );
    require(
        required_gate_obligation_ids == required_open_obligation_ids,
        "required gates and required unresolved obligations are not one-to-one",
        &mut errors,
    );
    require(
        open_required_gate_total == required_open_obligation_ids.len(),
        "open required-gate and required-obligation totals differ",
        &mut errors,
    );

    let blocked_obligations = plan
        .open_obligations
        .iter()
        .filter(|obligation| obligation.status == ObligationStatus::Blocked)
        .count();
    require(
        plan.completion_metrics.obligations.blocked == blocked_obligations,
        "global blocked-obligation count does not match obligations",
        &mut errors,
    );
    for obligation in &plan.open_obligations {
        require(
            EXPECTED_WORKSTREAMS.contains(&obligation.workstream_id.as_str()),
            &format!("{} has unknown workstream", obligation.id),
            &mut errors,
        );
        require(
            obligation.status != ObligationStatus::Resolved,
            &format!("resolved obligation {} must not remain open", obligation.id),
            &mut errors,
        );
    }

    validate_dependencies(plan, &workstream_ids, &mut errors);

    let lanes: Vec<&str> = plan
        .attack_routes
        .iter()
        .map(|route| route.lane.as_str())
        .collect();
    require(
        lanes == EXPECTED_LANES,
        "mock must contain one ordered skeleton record for lanes A-E",
        &mut errors,
    );
    require(
        unique(plan.attack_routes.iter().map(|route| route.id.as_str())),
        "attack-route ids must be unique",
        &mut errors,
    );
    for route in &plan.attack_routes {
        require(
            route.claim_status == ClaimStatus::Unknown,
            &format!("{} cannot be promoted in the planning mock", route.id),
            &mut errors,
        );
        require(
            obligation_by_id.contains_key(route.obligation_id.as_str()),
            &format!(
                "{} references missing obligation {}",
                route.id, route.obligation_id
            ),
            &mut errors,
        );
        require(
            route
                .id
                .starts_with(&format!("route.{}.", route.lane.to_ascii_lowercase())),
            &format!(
                "{} is not namespaced to attack-route lane {}",
                route.id, route.lane
            ),
            &mut errors,
        );
        if let Some(obligation) = obligation_by_id.get(route.obligation_id.as_str()) {
            require(
                obligation.workstream_id == "WS5",
                &format!(
                    "{} route obligation {} is not owned by WS5",
                    route.id, route.obligation_id
                ),
                &mut errors,
            );
            require(
                obligation.required
                    && matches!(
                        obligation.status,
                        ObligationStatus::Open | ObligationStatus::Blocked
                    ),
                &format!(
                    "{} route obligation {} is not required and unresolved",
                    route.id, route.obligation_id
                ),
                &mut errors,
            );
        }
    }
    require(
        plan.subjects.is_empty(),
        "planning mock cannot contain authoritative audit subjects",
        &mut errors,
    );
    require(
        plan.assessments.is_empty(),
        "planning mock cannot contain authoritative assessments",
        &mut errors,
    );
    require(
        plan.inventory_freeze.status == "not_frozen",
        "planning mock cannot freeze the authoritative denominator",
        &mut errors,
    );
    require(
        plan.inventory_freeze.frozen_at.is_none()
            && plan.inventory_freeze.denominator_manifest_sha256.is_none()
            && plan.inventory_freeze.expected_units.is_none(),
        "unfrozen planning inventory cannot claim a freeze time, denominator, or expected units",
        &mut errors,
    );

    let inventory = &plan.completion_metrics.inventory;
    require(
        inventory.expected.is_none() && inventory.present == 0 && inventory.missing.is_none(),
        "empty planning mock cannot claim inventory population or completion counts",
        &mut errors,
    );
    let classification = &plan.completion_metrics.classification;
    require(
        classification.required_claim_cells.is_none()
            && classification.classified_claim_cells == 0
            && classification.unknown == 0
            && classification.supported == 0
            && classification.refuted == 0
            && classification.not_applicable == 0,
        "empty planning mock cannot claim classified assessment cells",
        &mut errors,
    );
    let evidence = &plan.completion_metrics.evidence;
    require(
        evidence.independently_replayed == 0
            && evidence.matched_controlled == 0
            && evidence.end_to_end_certified == 0,
        "empty planning mock cannot claim evidence certification",
        &mut errors,
    );
    let authorization = &plan.completion_metrics.authorization;
    require(
        authorization.execution_units == 0
            && authorization.documented == 0
            && authorization.blocked_or_unknown == 0,
        "empty planning mock cannot claim authorization coverage",
        &mut errors,
    );
    let governance = &plan.completion_metrics.governance;
    require(
        governance.promotion_candidates == 0
            && governance.fully_reviewed == 0
            && governance.conflicted_reviews == 0,
        "empty planning mock cannot claim reviewed promotion state",
        &mut errors,
    );

    validate_source_artifacts(plan, repo_root, &mut errors);
    require(
        plan.report.record_path == EXPECTED_RECORD_PATH
            && plan.report.plan_path == EXPECTED_PLAN_PATH
            && plan.report.pdf_path == EXPECTED_PDF_PATH
            && plan.report.diagram_path == EXPECTED_DIAGRAM_PATH
            && plan.report.validator_path == EXPECTED_VALIDATOR_PATH,
        "planning mock report-artifact paths drifted",
        &mut errors,
    );
    require(
        plan.report.claim_boundary == EXPECTED_REPORT_CLAIM_BOUNDARY,
        "planning report claim boundary drifted",
        &mut errors,
    );
    let reported_record = validate_repo_file(
        repo_root,
        &plan.report.record_path,
        "machine planning record",
        &mut errors,
    );
    if let Some(reported_record) = reported_record {
        require(
            input_path.starts_with(repo_root),
            "input plan resolves outside repository",
            &mut errors,
        );
        require(
            reported_record == input_path,
            "input plan does not match report.record_path",
            &mut errors,
        );
    }
    for path in [
        plan.report.plan_path.as_str(),
        plan.report.pdf_path.as_str(),
        plan.report.diagram_path.as_str(),
        plan.report.validator_path.as_str(),
    ] {
        validate_repo_file(repo_root, path, "report artifact", &mut errors);
    }

    let derived = &plan.completion_metrics.derived;
    require(
        !derived.inventory_complete
            && !derived.classification_complete
            && !derived.authorization_complete
            && !derived.promotion_ready,
        "planning mock completion flags must all remain false",
        &mut errors,
    );
    errors
}

fn derive_workstream_status(workstream: &Workstream) -> RequirementStatus {
    let deliverables = &workstream.deliverables;
    if deliverables.is_empty() {
        return RequirementStatus::NotAttempted;
    }
    if deliverables
        .iter()
        .any(|item| item.status == RequirementStatus::Blocked)
        || workstream
            .acceptance_gates
            .iter()
            .any(|gate| gate.required && gate.claim_status == ClaimStatus::Refuted)
    {
        RequirementStatus::Blocked
    } else if deliverables
        .iter()
        .all(|item| item.status == RequirementStatus::VerifiedComplete)
        && workstream
            .acceptance_gates
            .iter()
            .filter(|gate| gate.required)
            .all(|gate| gate.claim_status == ClaimStatus::Supported)
    {
        RequirementStatus::VerifiedComplete
    } else if deliverables.iter().all(|item| {
        matches!(
            item.status,
            RequirementStatus::VerifiedComplete | RequirementStatus::ImplementedUnverified
        )
    }) {
        RequirementStatus::ImplementedUnverified
    } else if deliverables.iter().any(|item| {
        matches!(
            item.status,
            RequirementStatus::Partial
                | RequirementStatus::ImplementedUnverified
                | RequirementStatus::VerifiedComplete
        )
    }) {
        RequirementStatus::Partial
    } else {
        RequirementStatus::NotAttempted
    }
}

fn validate_dependencies(plan: &Plan, workstream_ids: &[&str], errors: &mut Vec<String>) {
    let mut triples = BTreeSet::new();
    let mut blocking: BTreeMap<&str, BTreeSet<&str>> = workstream_ids
        .iter()
        .map(|id| (*id, BTreeSet::new()))
        .collect();
    for dependency in &plan.dependencies {
        require(
            triples.insert((
                dependency.from.as_str(),
                dependency.to.as_str(),
                dependency.relation,
            )),
            "dependencies must be unique",
            errors,
        );
        require(
            dependency.from != dependency.to,
            &format!("self-dependency {}", dependency.from),
            errors,
        );
        require(
            workstream_ids.contains(&dependency.from.as_str())
                && workstream_ids.contains(&dependency.to.as_str()),
            &format!(
                "dependency has unknown endpoint: {}->{}",
                dependency.from, dependency.to
            ),
            errors,
        );
        if dependency.status == DependencyStatus::Active
            && matches!(
                dependency.relation,
                DependencyRelation::Requires | DependencyRelation::Blocks
            )
        {
            if let Some(targets) = blocking.get_mut(dependency.from.as_str()) {
                targets.insert(dependency.to.as_str());
            }
        }
    }

    let mut temporary = BTreeSet::new();
    let mut permanent = BTreeSet::new();
    for id in workstream_ids {
        visit_dependency(
            id,
            &blocking,
            &mut temporary,
            &mut permanent,
            &mut Vec::new(),
            errors,
        );
    }

    for workstream in &plan.workstreams {
        require(
            unique(workstream.depends_on.iter().map(String::as_str)),
            &format!("{} repeats a depends_on workstream", workstream.id),
            errors,
        );
        let declared: BTreeSet<&str> = workstream.depends_on.iter().map(String::as_str).collect();
        let listed: BTreeSet<&str> = plan
            .dependencies
            .iter()
            .filter(|dependency| {
                dependency.from == workstream.id
                    && dependency.relation == DependencyRelation::Requires
                    && dependency.status == DependencyStatus::Active
            })
            .map(|dependency| dependency.to.as_str())
            .collect();
        require(
            declared == listed,
            &format!("{} depends_on and requires edges differ", workstream.id),
            errors,
        );
    }
}

fn visit_dependency<'a>(
    node: &'a str,
    graph: &BTreeMap<&'a str, BTreeSet<&'a str>>,
    temporary: &mut BTreeSet<&'a str>,
    permanent: &mut BTreeSet<&'a str>,
    trail: &mut Vec<&'a str>,
    errors: &mut Vec<String>,
) {
    if permanent.contains(node) {
        return;
    }
    if !temporary.insert(node) {
        let mut cycle = trail.clone();
        cycle.push(node);
        errors.push(format!(
            "cyclic blocking dependency: {}",
            cycle.join(" -> ")
        ));
        return;
    }
    trail.push(node);
    if let Some(neighbors) = graph.get(node) {
        for neighbor in neighbors {
            visit_dependency(neighbor, graph, temporary, permanent, trail, errors);
        }
    }
    trail.pop();
    temporary.remove(node);
    permanent.insert(node);
}

fn validate_source_artifacts(plan: &Plan, repo_root: &Path, errors: &mut Vec<String>) {
    require(
        !plan.source_artifacts.is_empty(),
        "planning mock must retain at least one source artifact",
        errors,
    );
    require(
        unique(
            plan.source_artifacts
                .iter()
                .map(|source| source.id.as_str()),
        ),
        "source artifact ids must be unique",
        errors,
    );
    require(
        plan.source_artifacts.len() == EXPECTED_PINNED_SOURCES.len(),
        "planning mock source-artifact set drifted",
        errors,
    );
    for (id, expected_path, expected_revision, expected_digest) in EXPECTED_PINNED_SOURCES {
        match plan.source_artifacts.iter().find(|source| source.id == id) {
            Some(source) => require(
                source.path == expected_path
                    && source.status == SourceStatus::Pinned
                    && source.revision.as_deref() == Some(expected_revision)
                    && source.sha256.as_deref() == Some(expected_digest),
                &format!("pinned source provenance drifted: {id}"),
                errors,
            ),
            None => errors.push(format!("required pinned source is missing: {id}")),
        }
    }
    for source in &plan.source_artifacts {
        // A pinned source names immutable historical bytes via
        // (commit, repository path, SHA-256).  Its mutable worktree path may
        // legitimately have advanced since that baseline was frozen.
        let worktree_path = if source.status == SourceStatus::Pinned {
            None
        } else {
            validate_repo_file(repo_root, &source.path, "source artifact", errors)
        };
        match source.status {
            SourceStatus::Pinned => {
                let revision = source.revision.as_deref().unwrap_or_default();
                let revision_valid = is_lower_hex(revision, 40);
                require(
                    revision_valid,
                    &format!("pinned source lacks a 40-hex revision: {}", source.id),
                    errors,
                );
                let expected_digest = source.sha256.as_deref();
                let digest_valid = expected_digest.is_some_and(|hash| is_lower_hex(hash, 64));
                require(
                    digest_valid,
                    &format!("pinned source lacks a 64-hex digest: {}", source.id),
                    errors,
                );
                if revision_valid && digest_valid {
                    validate_pinned_revision(
                        repo_root,
                        source,
                        revision,
                        expected_digest.expect("digest validity checked"),
                        errors,
                    );
                }
            }
            SourceStatus::CurrentWorktree => {
                require(
                    source.revision.is_none(),
                    &format!(
                        "current-worktree source must not claim a committed revision: {}",
                        source.id
                    ),
                    errors,
                );
                require(
                    source
                        .sha256
                        .as_deref()
                        .is_some_and(|hash| is_lower_hex(hash, 64)),
                    &format!(
                        "current-worktree source lacks a 64-hex digest: {}",
                        source.id
                    ),
                    errors,
                );
            }
            SourceStatus::Planned => require(
                source.revision.is_none() && source.sha256.is_none(),
                &format!(
                    "planned source must not claim frozen provenance: {}",
                    source.id
                ),
                errors,
            ),
        }
        if let (Some(path), Some(expected)) = (worktree_path, source.sha256.as_deref()) {
            match read_bounded(&path, MAX_SOURCE_BYTES, "source artifact") {
                Ok(bytes) => require(
                    sha256_hex(&bytes) == expected,
                    &format!("source artifact hash mismatch: {}", source.path),
                    errors,
                ),
                Err(error) => errors.push(error),
            }
        }
    }
}

fn validate_pinned_revision(
    repo_root: &Path,
    source: &SourceArtifact,
    revision: &str,
    expected_digest: &str,
    errors: &mut Vec<String>,
) {
    let object = format!("{revision}:{}", source.path);
    match ProcessCommand::new("git")
        .arg("-C")
        .arg(repo_root)
        .args(["cat-file", "-t"])
        .arg(revision)
        .output()
    {
        Ok(output)
            if output.status.success()
                && String::from_utf8_lossy(&output.stdout).trim() == "commit" => {}
        Ok(output) => {
            errors.push(format!(
                "pinned source revision is not a commit: {} at {revision} (git exited {})",
                source.path, output.status
            ));
            return;
        }
        Err(error) => {
            errors.push(format!(
                "cannot inspect pinned source revision {} at {revision}: {error}",
                source.path
            ));
            return;
        }
    }
    match ProcessCommand::new("git")
        .arg("-C")
        .arg(repo_root)
        .args(["cat-file", "-s"])
        .arg(&object)
        .output()
    {
        Ok(output) if output.status.success() => {
            let size = String::from_utf8_lossy(&output.stdout)
                .trim()
                .parse::<u64>();
            match size {
                Ok(size) if size <= MAX_SOURCE_BYTES => {}
                Ok(size) => {
                    errors.push(format!(
                        "pinned source blob exceeds the {MAX_SOURCE_BYTES}-byte limit: {} at {revision} ({size} bytes)",
                        source.path
                    ));
                    return;
                }
                Err(error) => {
                    errors.push(format!(
                        "cannot parse pinned source blob size {} at {revision}: {error}",
                        source.path
                    ));
                    return;
                }
            }
        }
        Ok(output) => {
            errors.push(format!(
                "cannot inspect pinned source blob {} at {revision}: git cat-file exited {}",
                source.path, output.status
            ));
            return;
        }
        Err(error) => {
            errors.push(format!(
                "cannot inspect pinned source blob {} at {revision}: {error}",
                source.path
            ));
            return;
        }
    }
    match ProcessCommand::new("git")
        .arg("-C")
        .arg(repo_root)
        .args(["cat-file", "blob"])
        .arg(&object)
        .output()
    {
        Ok(output) if output.status.success() => require(
            sha256_hex(&output.stdout) == expected_digest,
            &format!(
                "pinned revision blob hash mismatch: {} at {revision}",
                source.path
            ),
            errors,
        ),
        Ok(output) => errors.push(format!(
            "cannot verify pinned source revision {} at {revision}: git cat-file exited {}",
            source.path, output.status
        )),
        Err(error) => errors.push(format!(
            "cannot verify pinned source revision {} at {revision}: {error}",
            source.path
        )),
    }
}

fn validate_repo_file(
    repo_root: &Path,
    raw_path: &str,
    kind: &str,
    errors: &mut Vec<String>,
) -> Option<PathBuf> {
    let relative = Path::new(raw_path);
    let safe = !relative.is_absolute()
        && relative
            .components()
            .all(|component| matches!(component, Component::Normal(_) | Component::CurDir));
    if !safe {
        errors.push(format!("{kind} path escapes repository: {raw_path}"));
        return None;
    }
    let joined = repo_root.join(relative);
    if !joined.is_file() {
        errors.push(format!("{kind} missing: {raw_path}"));
        return None;
    }
    match fs::canonicalize(&joined) {
        Ok(canonical_path) if canonical_path.starts_with(repo_root) => Some(canonical_path),
        Ok(_) => {
            errors.push(format!("{kind} resolves outside repository: {raw_path}"));
            None
        }
        Err(error) => {
            errors.push(format!("cannot resolve {kind} {raw_path}: {error}"));
            None
        }
    }
}

fn validate_deliverable_artifact(
    repo_root: &Path,
    deliverable: &Deliverable,
    errors: &mut Vec<String>,
) {
    match (
        deliverable.path.as_deref(),
        deliverable.artifact_sha256.as_deref(),
    ) {
        (None, None) => {}
        (None, Some(_)) => errors.push(format!(
            "{} has an artifact SHA-256 but no artifact path",
            deliverable.id
        )),
        (Some(path), None) => {
            validate_repo_file(repo_root, path, "deliverable artifact", errors);
            errors.push(format!(
                "{} artifact path is not bound by a SHA-256",
                deliverable.id
            ));
        }
        (Some(path), Some(expected_digest)) => {
            if !is_lower_hex(expected_digest, 64) {
                errors.push(format!(
                    "{} artifact SHA-256 is not 64 lowercase hexadecimal characters",
                    deliverable.id
                ));
                return;
            }
            let Some(canonical_path) =
                validate_repo_file(repo_root, path, "deliverable artifact", errors)
            else {
                return;
            };
            match read_bounded(&canonical_path, MAX_SOURCE_BYTES, "deliverable artifact") {
                Ok(bytes) => require(
                    sha256_hex(&bytes) == expected_digest,
                    &format!("{} artifact SHA-256 mismatch", deliverable.id),
                    errors,
                ),
                Err(error) => errors.push(error),
            }
        }
    }
}

fn unique<'a>(values: impl IntoIterator<Item = &'a str>) -> bool {
    let mut seen = BTreeSet::new();
    values.into_iter().all(|value| seen.insert(value))
}

fn strings_match(actual: &[String], expected: &[&str]) -> bool {
    actual
        .iter()
        .map(String::as_str)
        .eq(expected.iter().copied())
}

fn expected_plan_uid(raw_plan: &Value) -> Result<String, String> {
    let mut preimage = raw_plan.clone();
    let object = preimage
        .as_object_mut()
        .ok_or_else(|| "plan identity preimage is not an object".to_owned())?;
    if object.remove("plan_uid").is_none() {
        return Err("plan identity preimage lacks plan_uid".to_owned());
    }
    let encoded = canonical_utf8(&preimage)?;
    Ok(format!(
        "urn:ecdlp-audit-plan:1:sha256:{}",
        sha256_hex(encoded.as_bytes())
    ))
}

/// Canonical JSON for audit identities: sorted keys, compact encoding, raw
/// UTF-8, and no floats. This deliberately follows `docs/curve-identities.md`;
/// the ecbench identity codec escapes non-ASCII and is therefore not reused.
fn canonical_utf8(value: &Value) -> Result<String, String> {
    fn write(value: &Value, output: &mut String) -> Result<(), String> {
        match value {
            Value::Null => output.push_str("null"),
            Value::Bool(flag) => output.push_str(if *flag { "true" } else { "false" }),
            Value::Number(number) => {
                if number.is_f64() {
                    return Err(format!(
                        "float {number} in an identity; encode it as a string"
                    ));
                }
                output.push_str(&number.to_string());
            }
            Value::String(text) => output.push_str(
                &serde_json::to_string(text)
                    .map_err(|error| format!("cannot encode identity string: {error}"))?,
            ),
            Value::Array(items) => {
                output.push('[');
                for (index, item) in items.iter().enumerate() {
                    if index != 0 {
                        output.push(',');
                    }
                    write(item, output)?;
                }
                output.push(']');
            }
            Value::Object(map) => {
                output.push('{');
                let mut keys: Vec<&String> = map.keys().collect();
                keys.sort();
                for (index, key) in keys.iter().enumerate() {
                    if index != 0 {
                        output.push(',');
                    }
                    output.push_str(
                        &serde_json::to_string(key)
                            .map_err(|error| format!("cannot encode identity key: {error}"))?,
                    );
                    output.push(':');
                    write(&map[*key], output)?;
                }
                output.push('}');
            }
        }
        Ok(())
    }

    let mut output = String::new();
    write(value, &mut output)?;
    Ok(output)
}

fn is_lower_hex(value: &str, length: usize) -> bool {
    value.len() == length
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn require(condition: bool, message: &str, errors: &mut Vec<String>) {
    if !condition {
        errors.push(message.to_owned());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn root() -> PathBuf {
        fs::canonicalize(env!("CARGO_MANIFEST_DIR")).expect("repository root must resolve")
    }

    fn fixture_path() -> PathBuf {
        fs::canonicalize(root().join("docs/ecdlp-audit/ECDLP_AUDIT_PLAN.mock.json"))
            .expect("checked-in planning fixture must resolve")
    }

    fn fixture_value() -> Value {
        serde_json::from_str(include_str!(
            "../../docs/ecdlp-audit/ECDLP_AUDIT_PLAN.mock.json"
        ))
        .expect("checked-in planning fixture must parse as JSON")
    }

    fn fixture() -> Plan {
        serde_json::from_value(fixture_value()).expect("checked-in planning fixture must parse")
    }

    fn validate_fixture(plan: &Plan) -> Vec<String> {
        validate_plan(plan, &fixture_value(), &root(), &fixture_path())
    }

    #[test]
    fn checked_in_plan_passes() {
        assert!(is_rfc3339_date_time("1998-12-31T23:59:60Z"));
        assert!(is_rfc3339_date_time("1998-12-31T15:59:60.123-08:00"));
        assert!(is_rfc3339_date_time("2016-12-31T18:59:60-05:00"));
        assert!(is_rfc3339_date_time("2017-01-01T00:59:60+01:00"));
        assert!(!is_rfc3339_date_time("1998-12-31T23:59:60+01:00"));
        assert!(!is_rfc3339_date_time("1998-12-31T12:30:60+01:00"));
        assert!(!is_rfc3339_date_time("1998-12-31T23:59:61Z"));
        assert!(!is_rfc3339_date_time("2016-12-31T23:59:60-05:00"));
        assert!(parse_strict_json(PLAN_SCHEMA_JSON.as_bytes()).is_ok());
        assert_eq!(validate_plan_schema(&fixture_value()), Vec::<String>::new());
        assert_eq!(
            verify_plan_file(Path::new(EXPECTED_RECORD_PATH), &root()),
            Ok(())
        );
        assert_eq!(validate_fixture(&fixture()), Vec::<String>::new());
    }

    #[test]
    fn strict_json_parser_rejects_ambiguous_inputs() {
        assert_eq!(
            parse_strict_json(br#"{"a":[null,true,-1,2]}"#).unwrap(),
            json!({"a": [null, true, -1, 2]})
        );
        assert!(parse_strict_json(br#"{"count":18446744073709551616}   "#).is_ok());

        let cases: [(&[u8], &str); 6] = [
            (
                br#"{"schema_version":"a","schema_version":"b"}"#,
                "duplicate",
            ),
            (br#"{"properties":{"x":{},"x":{}}}"#, "duplicate"),
            (br#"{"count":1.0}"#, "number"),
            (br#"{"count":1e3}"#, "number"),
            (br#"{"count":1E+3}"#, "number"),
            (br#"{} false"#, "trailing"),
        ];
        for (input, expected) in cases {
            let error = parse_strict_json(input).expect_err("ambiguous JSON must fail closed");
            assert!(
                error.to_ascii_lowercase().contains(expected),
                "expected {expected:?} in strict parser error {error:?}"
            );
        }
    }

    #[test]
    fn schema_inspection_rejects_unsupported_constructs() {
        let cases = [
            (
                json!({
                    "$schema": PLAN_SCHEMA_DRAFT,
                    "type": "object",
                    "unevaluatedProperties": false
                }),
                "unsupported keyword",
            ),
            (
                json!({
                    "$schema": PLAN_SCHEMA_DRAFT,
                    "$ref": "https://example.invalid/schema.json"
                }),
                "remote, nested, escaped, and missing references are unsupported",
            ),
            (
                json!({
                    "$schema": PLAN_SCHEMA_DRAFT,
                    "type": "string",
                    "pattern": "^.*$"
                }),
                "unsupported pattern",
            ),
        ];
        for (schema, expected) in cases {
            let errors = validate_supported_plan_schema(&schema);
            assert!(
                errors.iter().any(|error| error.contains(expected)),
                "missing {expected:?} in {errors:?}"
            );
        }

        let mut excessive = serde_json::Map::new();
        for index in 0..300 {
            excessive.insert(format!("unexpected_{index}"), Value::Null);
        }
        let errors = validate_plan_schema(&Value::Object(excessive));
        assert!(errors.iter().any(|error| error.contains(
            "schema validation stopped after 256 violations; additional violations may exist"
        )));
    }

    #[test]
    fn schema_rejects_serde_valid_contract_violations() {
        let mut cases: Vec<(Value, &str)> = Vec::new();

        let mut missing_frozen_at = fixture_value();
        missing_frozen_at["inventory_freeze"]
            .as_object_mut()
            .expect("inventory freeze is an object")
            .remove("frozen_at");
        cases.push((missing_frozen_at, "/inventory_freeze"));

        let mut missing_nullable_path = fixture_value();
        missing_nullable_path["workstreams"][1]["deliverables"][1]
            .as_object_mut()
            .expect("deliverable is an object")
            .remove("path");
        cases.push((missing_nullable_path, "/workstreams/1/deliverables/1"));

        let mut invalid_date = fixture_value();
        invalid_date["as_of"] = Value::String("2026-99-99".to_owned());
        cases.push((invalid_date, "/as_of"));

        let mut duplicate_dependency = fixture_value();
        duplicate_dependency["workstreams"][1]["depends_on"]
            .as_array_mut()
            .expect("depends_on is an array")
            .push(Value::String("WS1".to_owned()));
        cases.push((duplicate_dependency, "/workstreams/1/depends_on"));

        let mut empty_deliverables = fixture_value();
        empty_deliverables["workstreams"][0]["deliverables"] = json!([]);
        cases.push((empty_deliverables, "/workstreams/0/deliverables"));

        let mut invalid_route_id = fixture_value();
        invalid_route_id["attack_routes"][0]["id"] = Value::String("not-a-route".to_owned());
        cases.push((invalid_route_id, "/attack_routes/0/id"));

        let mut empty_non_claims = fixture_value();
        empty_non_claims["workstreams"][0]["non_claims"] = json!([]);
        cases.push((empty_non_claims, "/workstreams/0/non_claims"));

        let mut unbound_partial_artifact = fixture_value();
        unbound_partial_artifact["workstreams"][0]["deliverables"][1]["artifact_sha256"] =
            Value::Null;
        cases.push((
            unbound_partial_artifact,
            "/workstreams/0/deliverables/1/artifact_sha256",
        ));

        let mut attempted_not_attempted = fixture_value();
        attempted_not_attempted["workstreams"][0]["deliverables"][0]["path"] =
            Value::String("Cargo.toml".to_owned());
        attempted_not_attempted["workstreams"][0]["deliverables"][0]["artifact_sha256"] =
            Value::String("0".repeat(64));
        cases.push((
            attempted_not_attempted,
            "/workstreams/0/deliverables/0/path",
        ));

        let mut weakened_scope_boundary = fixture_value();
        weakened_scope_boundary["scope"]["claim_boundary"] =
            Value::String("Operational recovery is authorized.".to_owned());
        cases.push((weakened_scope_boundary, "/scope/claim_boundary"));

        let mut weakened_report_boundary = fixture_value();
        weakened_report_boundary["report"]["claim_boundary"] =
            Value::String("The report certifies completion.".to_owned());
        cases.push((weakened_report_boundary, "/report/claim_boundary"));

        let mut unpinned_source = fixture_value();
        unpinned_source["source_artifacts"][0]["revision"] = Value::Null;
        cases.push((unpinned_source, "/source_artifacts/0/revision"));

        for (value, expected_path) in cases {
            assert!(
                serde_json::from_value::<Plan>(value.clone()).is_ok(),
                "negative schema case must remain Serde-valid: {expected_path}"
            );
            let errors = validate_plan_schema(&value);
            assert!(
                errors.iter().any(|error| error.contains(expected_path)),
                "missing schema failure at {expected_path}: {errors:?}"
            );
        }
    }

    #[test]
    fn file_entrypoint_enforces_embedded_schema() {
        let mut value = fixture_value();
        value["as_of"] = Value::String("2026-99-99".to_owned());
        let path = root().join(format!(
            "target/ecdlp-audit-schema-negative-{}.json",
            std::process::id()
        ));
        fs::create_dir_all(path.parent().expect("temporary path has a parent"))
            .expect("temporary directory can be created");
        fs::write(
            &path,
            serde_json::to_vec_pretty(&value).expect("negative fixture serializes"),
        )
        .expect("negative fixture can be written");
        let result = verify_plan_file(&path, &root());
        fs::remove_file(&path).expect("negative fixture can be removed");
        let errors = result.expect_err("schema-invalid plan must fail at the file entrypoint");
        assert!(errors
            .iter()
            .any(|error| error.contains("schema violation at /as_of")));
    }

    #[test]
    fn schema_requires_explicit_deliverable_artifact_identity() {
        let mut value = fixture_value();
        value["workstreams"][0]["deliverables"][1]
            .as_object_mut()
            .expect("deliverable is an object")
            .remove("artifact_sha256");
        let errors = validate_plan_schema(&value);
        assert!(errors
            .iter()
            .any(|error| error.contains("/workstreams/0/deliverables/1")));
        assert!(serde_json::from_value::<Plan>(value).is_ok());
    }

    #[test]
    fn promotion_fails_closed() {
        let mut plan = fixture();
        plan.promotion_allowed = true;
        let errors = validate_fixture(&plan);
        assert!(errors
            .iter()
            .any(|error| error.contains("cannot allow promotion")));
    }

    #[test]
    fn missing_unknown_obligation_fails_closed() {
        let mut plan = fixture();
        plan.workstreams[0].acceptance_gates[0].obligation_id = None;
        let errors = validate_fixture(&plan);
        assert!(errors
            .iter()
            .any(|error| error.contains("unknown status lacks an obligation")));
    }

    #[test]
    fn dependency_cycle_is_rejected() {
        let mut plan = fixture();
        plan.dependencies.push(Dependency {
            from: "WS1".to_owned(),
            to: "WS2".to_owned(),
            relation: DependencyRelation::Requires,
            status: DependencyStatus::Active,
            rationale: "negative fixture".to_owned(),
        });
        plan.workstreams[0].depends_on.push("WS2".to_owned());
        let errors = validate_fixture(&plan);
        assert!(errors
            .iter()
            .any(|error| error.contains("cyclic blocking dependency")));
    }

    #[test]
    fn repository_escape_is_rejected() {
        let mut plan = fixture();
        plan.report.validator_path = "../outside".to_owned();
        let errors = validate_fixture(&plan);
        assert!(errors
            .iter()
            .any(|error| error.contains("path escapes repository")));
    }

    #[test]
    fn self_promotion_fails_closed() {
        let mut plan = fixture();
        plan.governance.self_promotion_allowed = true;
        let errors = validate_fixture(&plan);
        assert!(errors
            .iter()
            .any(|error| error.contains("cannot allow self-promotion")));
    }

    #[test]
    fn mock_contract_vocabulary_and_sources_fail_closed() {
        let mut plan = fixture();
        plan.vocabularies.claim_status.clear();
        plan.source_artifacts.clear();
        let errors = validate_fixture(&plan);
        assert!(errors
            .iter()
            .any(|error| error.contains("vocabularies drifted")));
        assert!(errors
            .iter()
            .any(|error| error.contains("at least one source artifact")));
    }

    #[test]
    fn positive_completion_claims_fail_closed() {
        let mut plan = fixture();
        plan.completion_metrics.evidence.end_to_end_certified = 1;
        plan.completion_metrics.governance.fully_reviewed = 1;
        plan.workstreams[0].requirement_status = RequirementStatus::VerifiedComplete;
        let errors = validate_fixture(&plan);
        assert!(errors
            .iter()
            .any(|error| error.contains("cannot claim evidence certification")));
        assert!(errors
            .iter()
            .any(|error| error.contains("cannot claim reviewed promotion state")));
        assert!(errors
            .iter()
            .any(|error| error.contains("cannot be verified complete")));
    }

    #[test]
    fn required_gate_cannot_use_optional_obligation() {
        let mut plan = fixture();
        let obligation_id = plan.workstreams[0].acceptance_gates[0]
            .obligation_id
            .clone()
            .expect("fixture gate has an obligation");
        plan.open_obligations
            .iter_mut()
            .find(|obligation| obligation.id == obligation_id)
            .expect("fixture obligation exists")
            .required = false;
        let errors = validate_fixture(&plan);
        assert!(errors
            .iter()
            .any(|error| error.contains("required gate references optional obligation")));
    }

    #[test]
    fn progressed_deliverable_requires_artifact() {
        for status in [
            RequirementStatus::Partial,
            RequirementStatus::ImplementedUnverified,
            RequirementStatus::VerifiedComplete,
        ] {
            let mut plan = fixture();
            plan.workstreams[0].deliverables[0].status = status;
            plan.workstreams[0].deliverables[0].path = None;
            plan.workstreams[0].deliverables[0].artifact_sha256 = None;
            let errors = validate_fixture(&plan);
            assert!(errors
                .iter()
                .any(|error| { error.contains("WS1-D1") && error.contains("no artifact path") }));
        }
    }

    #[test]
    fn not_attempted_deliverable_cannot_bind_an_artifact() {
        let mut plan = fixture();
        plan.workstreams[0].deliverables[0].path = Some("Cargo.toml".to_owned());
        plan.workstreams[0].deliverables[0].artifact_sha256 = Some("0".repeat(64));
        let errors = validate_fixture(&plan);
        assert!(errors
            .iter()
            .any(|error| error.contains("WS1-D1 is not attempted but binds an artifact")));
    }

    #[test]
    fn workstream_status_is_derived_from_deliverables() {
        let mut plan = fixture();
        plan.workstreams[0].requirement_status = RequirementStatus::ImplementedUnverified;
        let errors = validate_fixture(&plan);
        assert!(errors.iter().any(|error| {
            error.contains("WS1 status ImplementedUnverified")
                && error.contains("deliverable-derived status Partial")
        }));
    }

    #[test]
    fn workstream_status_derivation_covers_progress_and_blockers() {
        let mut plan = fixture();
        let mut workstream = plan.workstreams.remove(0);
        assert_eq!(
            derive_workstream_status(&workstream),
            RequirementStatus::Partial
        );

        for deliverable in &mut workstream.deliverables {
            deliverable.status = RequirementStatus::NotAttempted;
        }
        assert_eq!(
            derive_workstream_status(&workstream),
            RequirementStatus::NotAttempted
        );

        workstream.deliverables[0].status = RequirementStatus::Blocked;
        assert_eq!(
            derive_workstream_status(&workstream),
            RequirementStatus::Blocked
        );

        for deliverable in &mut workstream.deliverables {
            deliverable.status = RequirementStatus::ImplementedUnverified;
        }
        assert_eq!(
            derive_workstream_status(&workstream),
            RequirementStatus::ImplementedUnverified
        );

        for deliverable in &mut workstream.deliverables {
            deliverable.status = RequirementStatus::VerifiedComplete;
        }
        assert_eq!(
            derive_workstream_status(&workstream),
            RequirementStatus::ImplementedUnverified
        );
        for gate in workstream
            .acceptance_gates
            .iter_mut()
            .filter(|gate| gate.required)
        {
            gate.claim_status = ClaimStatus::Supported;
        }
        assert_eq!(
            derive_workstream_status(&workstream),
            RequirementStatus::VerifiedComplete
        );
    }

    #[test]
    fn deliverable_artifacts_are_exactly_path_and_content_bound() {
        let mut unbound = fixture();
        unbound.workstreams[0].deliverables[1].artifact_sha256 = None;
        let errors = validate_fixture(&unbound);
        assert!(errors
            .iter()
            .any(|error| error.contains("WS1-D2 artifact path is not bound by a SHA-256")));

        let mut wrong_path = fixture();
        wrong_path.workstreams[0].deliverables[1].path = Some("Cargo.toml".to_owned());
        let errors = validate_fixture(&wrong_path);
        assert!(errors
            .iter()
            .any(|error| error.contains("WS1-D2 artifact path drifted")));

        let mut wrong_kind = fixture();
        wrong_kind.workstreams[0].deliverables[1].kind = "implementation".to_owned();
        let errors = validate_fixture(&wrong_kind);
        assert!(errors
            .iter()
            .any(|error| error.contains("WS1-D2 artifact kind drifted")));

        let mut wrong_digest = fixture();
        wrong_digest.workstreams[0].deliverables[1].artifact_sha256 = Some("0".repeat(64));
        let errors = validate_fixture(&wrong_digest);
        assert!(errors
            .iter()
            .any(|error| error.contains("WS1-D2 artifact SHA-256 mismatch")));
    }

    #[test]
    fn every_required_gate_must_be_owned_by_a_deliverable() {
        let mut plan = fixture();
        plan.workstreams[0].deliverables[0]
            .acceptance_gate_ids
            .retain(|gate_id| gate_id != "WS1-G2");
        let errors = validate_fixture(&plan);
        assert!(errors.iter().any(|error| {
            error.contains("WS1 required gate WS1-G2 is not referenced by a deliverable")
        }));
    }

    #[test]
    fn planning_claim_boundaries_are_pinned() {
        let mut plan = fixture();
        plan.scope.claim_boundary = "Operational solver authorization.".to_owned();
        plan.report.claim_boundary = "Complete and safe.".to_owned();
        let errors = validate_fixture(&plan);
        assert!(errors
            .iter()
            .any(|error| error.contains("planning scope claim boundary drifted")));
        assert!(errors
            .iter()
            .any(|error| error.contains("planning report claim boundary drifted")));
    }

    #[test]
    fn pinned_source_cannot_be_downgraded() {
        let mut plan = fixture();
        plan.source_artifacts[0].status = SourceStatus::CurrentWorktree;
        plan.source_artifacts[0].revision = None;
        let errors = validate_fixture(&plan);
        assert!(errors
            .iter()
            .any(|error| error.contains("pinned source provenance drifted")));
    }

    #[test]
    fn pinned_source_replays_historical_blob_despite_worktree_divergence() {
        let plan = fixture();
        let source = plan
            .source_artifacts
            .iter()
            .find(|source| source.id == "source.cover-transfer-ledger")
            .expect("fixture contains the historical cover-ledger baseline");
        let revision = source
            .revision
            .as_deref()
            .expect("pinned source has a revision");
        let digest = source
            .sha256
            .as_deref()
            .expect("pinned source has a digest");
        let worktree_bytes = fs::read(root().join(&source.path))
            .expect("current generated cover ledger can be read");
        assert_ne!(
            sha256_hex(&worktree_bytes),
            digest,
            "regression requires the current generated ledger to differ from its historical baseline"
        );

        let mut errors = Vec::new();
        validate_pinned_revision(&root(), source, revision, digest, &mut errors);
        assert_eq!(errors, Vec::<String>::new());

        let mut wrong_digest_errors = Vec::new();
        validate_pinned_revision(
            &root(),
            source,
            revision,
            &"0".repeat(64),
            &mut wrong_digest_errors,
        );
        assert!(wrong_digest_errors
            .iter()
            .any(|error| error.contains("pinned revision blob hash mismatch")));

        let mut wrong_path_source = source.clone();
        wrong_path_source.path = "Cargo.toml".to_owned();
        let mut wrong_path_errors = Vec::new();
        validate_pinned_revision(
            &root(),
            &wrong_path_source,
            revision,
            digest,
            &mut wrong_path_errors,
        );
        assert!(wrong_path_errors
            .iter()
            .any(|error| error.contains("pinned revision blob hash mismatch")));

        let mut wrong_revision_errors = Vec::new();
        validate_pinned_revision(
            &root(),
            source,
            &"0".repeat(40),
            digest,
            &mut wrong_revision_errors,
        );
        assert!(wrong_revision_errors
            .iter()
            .any(|error| error.contains("pinned source revision is not a commit")));
    }

    #[test]
    fn attack_routes_require_unique_lane_namespaces_and_ws5_obligations() {
        let mut duplicate = fixture();
        duplicate.attack_routes[1].id = duplicate.attack_routes[0].id.clone();
        assert!(validate_fixture(&duplicate)
            .iter()
            .any(|error| error.contains("attack-route ids must be unique")));

        let mut wrong_lane = fixture();
        wrong_lane.attack_routes[0].id = "route.b.wrong-lane".to_owned();
        assert!(validate_fixture(&wrong_lane)
            .iter()
            .any(|error| error.contains("not namespaced to attack-route lane A")));

        let mut wrong_owner = fixture();
        wrong_owner.attack_routes[0].obligation_id = "evidence-profile-registry".to_owned();
        assert!(validate_fixture(&wrong_owner)
            .iter()
            .any(|error| error.contains("is not owned by WS5")));
    }

    #[test]
    fn dependency_status_controls_active_graph_and_declarations() {
        let mut downgraded = fixture();
        downgraded.dependencies[0].status = DependencyStatus::Planned;
        assert!(validate_fixture(&downgraded)
            .iter()
            .any(|error| error.contains("depends_on and requires edges differ")));

        let mut planned_reciprocal = fixture();
        planned_reciprocal.dependencies.push(Dependency {
            from: "WS1".to_owned(),
            to: "WS2".to_owned(),
            relation: DependencyRelation::Requires,
            status: DependencyStatus::Planned,
            rationale: "negative fixture".to_owned(),
        });
        let errors = validate_fixture(&planned_reciprocal);
        assert!(!errors
            .iter()
            .any(|error| error.contains("cyclic blocking dependency")));
    }

    #[test]
    fn duplicate_depends_on_and_gate_obligations_are_rejected() {
        let mut duplicate_dependency = fixture();
        duplicate_dependency.workstreams[1]
            .depends_on
            .push("WS1".to_owned());
        assert!(validate_fixture(&duplicate_dependency)
            .iter()
            .any(|error| error.contains("repeats a depends_on workstream")));

        let mut duplicate_obligation = fixture();
        duplicate_obligation.workstreams[0].acceptance_gates[1].obligation_id =
            duplicate_obligation.workstreams[0].acceptance_gates[0]
                .obligation_id
                .clone();
        let errors = validate_fixture(&duplicate_obligation);
        assert!(errors
            .iter()
            .any(|error| error.contains("reuses required-gate obligation")));
        assert!(errors
            .iter()
            .any(|error| error.contains("are not one-to-one")));
    }

    #[test]
    fn open_required_metric_counts_obligations() {
        let mut plan = fixture();
        plan.open_obligations[0].required = false;
        let errors = validate_fixture(&plan);
        assert!(errors.iter().any(|error| {
            error.contains("open-required count does not match required unresolved obligations")
        }));
    }

    #[test]
    fn input_plan_must_match_reported_plan() {
        let plan = fixture();
        let errors = validate_plan(&plan, &fixture_value(), &root(), &root().join("Cargo.toml"));
        assert!(errors
            .iter()
            .any(|error| error.contains("does not match report.record_path")));
    }

    #[test]
    fn input_plan_must_resolve_inside_repository() {
        let plan = fixture();
        let outside = root()
            .parent()
            .expect("repository has a parent")
            .join("outside-plan.json");
        let errors = validate_plan(&plan, &fixture_value(), &root(), &outside);
        assert!(errors
            .iter()
            .any(|error| error.contains("input plan resolves outside repository")));
    }

    #[test]
    fn plan_uid_binds_full_record() {
        let mut value = fixture_value();
        value["scope"]["objective"] = Value::String("mutated objective".to_owned());
        let plan: Plan = serde_json::from_value(value.clone()).expect("mutation stays typed");
        let errors = validate_plan(&plan, &value, &root(), &fixture_path());
        assert!(errors
            .iter()
            .any(|error| error.contains("canonical full-plan preimage")));
    }

    #[test]
    fn unknown_top_level_field_is_rejected() {
        let mut value: Value = serde_json::from_str(include_str!(
            "../../docs/ecdlp-audit/ECDLP_AUDIT_PLAN.mock.json"
        ))
        .unwrap();
        value["unexpected"] = Value::Bool(true);
        assert!(serde_json::from_value::<Plan>(value).is_err());
    }

    #[test]
    fn malformed_json_is_rejected() {
        assert!(serde_json::from_str::<Plan>("{not json").is_err());
    }

    #[test]
    fn identity_codec_uses_sorted_raw_utf8() {
        let value = json!({"é": "雪", "a": [2, 1]});
        let encoded = canonical_utf8(&value).unwrap();
        assert_eq!(encoded, "{\"a\":[2,1],\"é\":\"雪\"}");
        assert_eq!(
            sha256_hex(encoded.as_bytes()),
            "47eda5b8dc14a7e4bb584c33e9a4acbc92c199b000d84dbea2035171cbd03661"
        );
        assert!(canonical_utf8(&json!({"float": 1.5})).is_err());
    }
}
