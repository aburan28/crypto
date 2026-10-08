//! Round 301: exact P-256 standard algebraic-transfer census.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::ecc::CurveParams;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::Serialize;
use serde_json::{json, Value};

const ROUND25_SHA256: &str = "dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf";
const ROUND28_SHA256: &str = "950d4b7f4997e146ee30fc7c61174070595fa5d344b9a30cbc0ca578ae7e60f4";
const ROUND297_SHA256: &str = "53949881daf4a48769c9ef8f869fad633c3dc69d11b0d341a5b7416de2e43982";
const ROUND300_SHA256: &str = "7d65f3e015c6d6c6d24b64ea2e59a6e9e8653c88a1900835838b1bc0a1ec36e5";
const FACTOR_BASE_ID: &str = "FB1hc72514a2a8d3";
const REGISTERED_S17_FACTOR_BASE_ID: &str = "FB1h2f8621cda105";
const ROUND295_RATIO_TO_RHO: f64 = 13.920_747_397_073_491;
const REGISTERED_S17_RATIO_TO_RHO: f64 = 394.425_280;

#[derive(Parser)]
#[command(about = "Certify standard algebraic transfer boundaries for P-256")]
struct Cli {
    #[arg(long)]
    round25: PathBuf,
    #[arg(long)]
    round28: PathBuf,
    #[arg(long)]
    round297: PathBuf,
    #[arg(long)]
    round300: PathBuf,
    #[arg(long)]
    out: PathBuf,
    #[arg(long)]
    assessment: PathBuf,
}

#[derive(Clone, Serialize)]
struct Dependency {
    path: String,
    bytes: u64,
    sha256: String,
    schema: String,
}

#[derive(Clone, Default, Serialize)]
struct OperationCounts {
    modular_exponentiations: u64,
    modular_multiplications: u64,
    modular_squarings: u64,
    integer_multiplications: u64,
    exact_divisions: u64,
    factor_divisibility_tests: u64,
}

#[derive(Clone, Serialize)]
struct PrimePower {
    prime: String,
    exponent: u64,
}

#[derive(Serialize)]
struct OrderReductionStep {
    prime: String,
    occurrence: u64,
    candidate_order: String,
    residue: String,
    divided: bool,
}

#[derive(Serialize)]
struct MinimalityCheck {
    prime: String,
    exponent_in_order: u64,
    tested_exponent: String,
    residue: String,
    residue_is_one: bool,
}

#[derive(Serialize)]
struct CurveInvariants {
    p: String,
    a: String,
    b: String,
    n: String,
    cofactor: String,
    trace: String,
    frobenius_discriminant: String,
    frobenius_discriminant_abs: String,
    anomalous: bool,
    supersingular: bool,
    ordinary: bool,
    j_invariant: String,
    j_is_zero: bool,
    j_is_1728: bool,
    exceptional_automorphism_case: bool,
    prime_field_has_proper_smaller_subfield: bool,
}

#[derive(Serialize)]
struct EmbeddingDegreeCertificate {
    definition: String,
    n_minus_one_factorization: Vec<PrimePower>,
    factorization_product_verified: bool,
    reduction_steps: Vec<OrderReductionStep>,
    exact_embedding_degree: String,
    exact_embedding_degree_bits: u64,
    embedding_degree_factorization: Vec<PrimePower>,
    final_residue: String,
    final_residue_is_one: bool,
    minimality_checks: Vec<MinimalityCheck>,
    all_minimality_checks_pass: bool,
    low_degree_at_most_64: bool,
    extension_field_bit_dimension: String,
    extension_is_smaller_representation: bool,
}

#[derive(Serialize)]
struct ImportedCertificate {
    family: String,
    source_round: u64,
    checked_facts: Vec<String>,
    applicable_non_generic_quotient: bool,
    status: String,
}

#[derive(Serialize)]
struct BoundaryRow {
    variant: String,
    optimistic_ratio_to_rho: Option<f64>,
    complete_cost_measured: bool,
    status: String,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_and_schemas_checked: bool,
    curve_invariants_recomputed: bool,
    n_minus_one_factorization_reverified: bool,
    exact_embedding_degree_certified: bool,
    imported_transfer_certificates_reconciled: bool,
    applicable_non_generic_transfer_found: bool,
    explicit_forward_inverse_and_recovery_implemented: bool,
    structured_residual_degree_at_most_5: bool,
    parity_at_or_below_rho: bool,
    per_usable_relation_below_2_103: bool,
    complete_collection_below_2_120: bool,
    projected_storage_below_2_50: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultReceipt {
    schema: String,
    curve: String,
    screening_round: u64,
    execution_status: String,
    factor_base: String,
    registered_s17_factor_base: String,
    dependencies: Vec<Dependency>,
    curve_invariants: CurveInvariants,
    embedding_degree: EmbeddingDegreeCertificate,
    imported_certificates: Vec<ImportedCertificate>,
    operation_counts: OperationCounts,
    maximum_integer_bits: u64,
    boundary_table: Vec<BoundaryRow>,
    named_families_complete: bool,
    exploration_boundary: String,
    relations_reported_on_p256: u64,
    full_depth_unplanted_p256_relation_attempted: bool,
    transfer_assessment_semantic_sha256: String,
    gates: Gates,
    classification: String,
    dominant_obstruction: String,
    decision: String,
    semantic_evidence_sha256: String,
    result_json_bytes: u64,
}

#[derive(Serialize)]
struct Obligation {
    name: String,
    status: String,
    evidence: String,
    scope: String,
}

#[derive(Serialize)]
struct TransferFamily {
    name: String,
    source: String,
    destination: String,
    transformation_type: String,
    field_of_definition: String,
    source_dimension: Option<u64>,
    destination_dimension: Option<u64>,
    relevant_subgroup: String,
    obligations: Vec<Obligation>,
    conclusion: String,
    exploration_boundary: String,
}

#[derive(Serialize)]
struct TransferAssessment {
    schema: String,
    curve: String,
    screening_round: u64,
    skill_profile: String,
    required_companion_resources_available: bool,
    methodology_resource: String,
    template_resource: String,
    typed_correspondence_graph: Vec<String>,
    families: Vec<TransferFamily>,
    controls: Vec<String>,
    cost_accounting: Vec<String>,
    weakest_open_obligation: String,
    narrowest_supported_finding: String,
    semantic_evidence_sha256: String,
    assessment_json_bytes: u64,
}

fn checked_dependency(
    path: &Path,
    expected_hash: &str,
    expected_schema: &str,
) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = sha256_hex(&bytes);
    if digest != expected_hash {
        return Err(format!(
            "{} SHA-256 is {digest}, expected {expected_hash}",
            path.display()
        ));
    }
    let value: Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("{}: {error}", path.display()))?;
    if value.pointer("/schema").and_then(Value::as_str) != Some(expected_schema) {
        return Err(format!("{} schema mismatch", path.display()));
    }
    Ok((
        Dependency {
            path: path.display().to_string(),
            bytes: bytes.len() as u64,
            sha256: digest,
            schema: expected_schema.into(),
        },
        value,
    ))
}

fn parse_big(value: &str, label: &str) -> Result<BigUint, String> {
    BigUint::parse_bytes(value.as_bytes(), 10).ok_or_else(|| format!("invalid {label}"))
}

fn modpow_counted(
    base: &BigUint,
    exponent: &BigUint,
    modulus: &BigUint,
    operations: &mut OperationCounts,
) -> BigUint {
    operations.modular_exponentiations += 1;
    let mut power = base % modulus;
    let mut exponent = exponent.clone();
    let mut out = BigUint::one();
    while !exponent.is_zero() {
        if exponent.bit(0) {
            out = (&out * &power) % modulus;
            operations.modular_multiplications += 1;
        }
        power = (&power * &power) % modulus;
        operations.modular_squarings += 1;
        exponent >>= 1usize;
    }
    out
}

fn mod_inverse_prime(
    value: &BigUint,
    modulus: &BigUint,
    operations: &mut OperationCounts,
) -> BigUint {
    modpow_counted(value, &(modulus - BigUint::from(2u8)), modulus, operations)
}

fn modular_j(curve: &CurveParams, operations: &mut OperationCounts) -> Result<BigUint, String> {
    let p = &curve.p;
    let a2 = (&curve.a * &curve.a) % p;
    let a3 = (&a2 * &curve.a) % p;
    let b2 = (&curve.b * &curve.b) % p;
    operations.modular_multiplications += 3;
    let four_a3 = (BigUint::from(4u8) * a3) % p;
    let twenty_seven_b2 = (BigUint::from(27u8) * b2) % p;
    operations.modular_multiplications += 2;
    let denominator = (&four_a3 + twenty_seven_b2) % p;
    if denominator.is_zero() {
        return Err("singular P-256 curve".into());
    }
    let inverse = mod_inverse_prime(&denominator, p, operations);
    operations.modular_multiplications += 2;
    Ok((BigUint::from(1728u16) * four_a3 * inverse) % p)
}

fn factors_from_round25(value: &Value) -> Result<Vec<(BigUint, u64)>, String> {
    value
        .pointer("/scalar_order_screen/n_minus_one_factorization")
        .and_then(Value::as_array)
        .ok_or("Round 25 lacks n-1 factorization")?
        .iter()
        .map(|item| {
            let prime = item
                .pointer("/prime")
                .and_then(Value::as_str)
                .ok_or("factor prime missing")?;
            let exponent = item
                .pointer("/exponent")
                .and_then(Value::as_u64)
                .ok_or("factor exponent missing")?;
            Ok((parse_big(prime, "factor prime")?, exponent))
        })
        .collect()
}

fn factorization_product(factors: &[(BigUint, u64)], operations: &mut OperationCounts) -> BigUint {
    let mut product = BigUint::one();
    for (prime, exponent) in factors {
        for _ in 0..*exponent {
            product *= prime;
            operations.integer_multiplications += 1;
        }
    }
    product
}

fn embedding_degree(
    p: &BigUint,
    n: &BigUint,
    factors: &[(BigUint, u64)],
    operations: &mut OperationCounts,
) -> Result<EmbeddingDegreeCertificate, String> {
    let n_minus_one = n - BigUint::one();
    let product = factorization_product(factors, operations);
    if product != n_minus_one {
        return Err("Round 25 factorization does not multiply to n-1".into());
    }
    let mut order = n_minus_one;
    let mut remaining = factors.to_vec();
    let mut steps = Vec::new();
    for (factor_index, (prime, exponent)) in factors.iter().enumerate() {
        for occurrence in 1..=*exponent {
            operations.factor_divisibility_tests += 1;
            if &order % prime != BigUint::zero() {
                break;
            }
            let candidate = &order / prime;
            operations.exact_divisions += 1;
            let residue = modpow_counted(p, &candidate, n, operations);
            let divided = residue == BigUint::one();
            steps.push(OrderReductionStep {
                prime: prime.to_string(),
                occurrence,
                candidate_order: candidate.to_string(),
                residue: residue.to_string(),
                divided,
            });
            if divided {
                order = candidate;
                remaining[factor_index].1 -= 1;
            } else {
                break;
            }
        }
    }
    let final_residue = modpow_counted(p, &order, n, operations);
    let remaining = remaining
        .into_iter()
        .filter(|(_, exponent)| *exponent != 0)
        .collect::<Vec<_>>();
    let mut minimality = Vec::new();
    for (prime, exponent) in &remaining {
        let tested = &order / prime;
        operations.exact_divisions += 1;
        let residue = modpow_counted(p, &tested, n, operations);
        minimality.push(MinimalityCheck {
            prime: prime.to_string(),
            exponent_in_order: *exponent,
            tested_exponent: tested.to_string(),
            residue_is_one: residue == BigUint::one(),
            residue: residue.to_string(),
        });
    }
    let final_is_one = final_residue == BigUint::one();
    let minimal = minimality.iter().all(|check| !check.residue_is_one);
    if !final_is_one || !minimal {
        return Err("embedding-degree certificate failed".into());
    }
    let extension_bits = &order * BigUint::from(curve_field_bits(p));
    Ok(EmbeddingDegreeCertificate {
        definition: "k=ord_n(p), the least positive k with n | p^k-1".into(),
        n_minus_one_factorization: factors
            .iter()
            .map(|(prime, exponent)| PrimePower {
                prime: prime.to_string(),
                exponent: *exponent,
            })
            .collect(),
        factorization_product_verified: true,
        reduction_steps: steps,
        exact_embedding_degree: order.to_string(),
        exact_embedding_degree_bits: order.bits(),
        embedding_degree_factorization: remaining
            .iter()
            .map(|(prime, exponent)| PrimePower {
                prime: prime.to_string(),
                exponent: *exponent,
            })
            .collect(),
        final_residue: final_residue.to_string(),
        final_residue_is_one: final_is_one,
        minimality_checks: minimality,
        all_minimality_checks_pass: minimal,
        low_degree_at_most_64: order <= BigUint::from(64u8),
        extension_field_bit_dimension: extension_bits.to_string(),
        extension_is_smaller_representation: false,
    })
}

fn curve_field_bits(p: &BigUint) -> u64 {
    p.bits()
}

fn curve_invariants(
    curve: &CurveParams,
    round25: &Value,
    operations: &mut OperationCounts,
) -> Result<CurveInvariants, String> {
    if curve.h != 1 {
        return Err("P-256 cofactor is not one".into());
    }
    let p_plus_one = &curve.p + BigUint::one();
    if p_plus_one < curve.n {
        return Err("unexpected negative P-256 trace".into());
    }
    let trace = p_plus_one - &curve.n;
    let trace_squared = &trace * &trace;
    operations.integer_multiplications += 1;
    let four_p = BigUint::from(4u8) * &curve.p;
    operations.integer_multiplications += 1;
    if trace_squared >= four_p {
        return Err("invalid Frobenius discriminant sign".into());
    }
    let discriminant_abs = four_p - trace_squared;
    let imported_trace = round25
        .pointer("/cm_screen/trace")
        .and_then(Value::as_str)
        .ok_or("Round 25 trace missing")?;
    let imported_discriminant_abs = round25
        .pointer("/cm_screen/discriminant_abs")
        .and_then(Value::as_str)
        .ok_or("Round 25 discriminant missing")?;
    if trace.to_string() != imported_trace
        || discriminant_abs.to_string() != imported_discriminant_abs
    {
        return Err("Round 25 trace/discriminant mismatch".into());
    }
    let j = modular_j(curve, operations)?;
    let j_1728 = BigUint::from(1728u16) % &curve.p;
    let supersingular = trace.is_zero();
    Ok(CurveInvariants {
        p: curve.p.to_string(),
        a: curve.a.to_string(),
        b: curve.b.to_string(),
        n: curve.n.to_string(),
        cofactor: curve.h.to_string(),
        trace: trace.to_string(),
        frobenius_discriminant: format!("-{}", discriminant_abs),
        frobenius_discriminant_abs: discriminant_abs.to_string(),
        anomalous: curve.n == curve.p,
        supersingular,
        ordinary: !supersingular,
        j_is_zero: j.is_zero(),
        j_is_1728: j == j_1728,
        exceptional_automorphism_case: j.is_zero() || j == j_1728,
        j_invariant: j.to_string(),
        prime_field_has_proper_smaller_subfield: false,
    })
}

fn imported_certificates(
    round25: &Value,
    round28: &Value,
    round297: &Value,
    round300: &Value,
) -> Result<Vec<ImportedCertificate>, String> {
    let min_degree = round25
        .pointer("/cm_screen/minimum_noninteger_endomorphism_degree")
        .and_then(Value::as_str)
        .ok_or("Round 25 minimum degree missing")?;
    if round25
        .pointer("/cm_screen/frobenius_conductor")
        .and_then(Value::as_str)
        != Some("1")
        || round25
            .pointer("/cm_screen/low_degree_noninteger_transport_exists")
            .and_then(Value::as_bool)
            != Some(false)
    {
        return Err("Round 25 CM certificate mismatch".into());
    }
    let edges = round28
        .pointer("/certificate_audit/edges")
        .and_then(Value::as_u64)
        .ok_or("Round 28 edge count missing")?;
    let models = round28
        .pointer("/certificate_audit/models")
        .and_then(Value::as_u64)
        .ok_or("Round 28 model count missing")?;
    if edges != 4096
        || models != 4097
        || round28
            .pointer("/certificate_audit/all_checks_passed")
            .and_then(Value::as_bool)
            != Some(true)
        || round28
            .pointer("/promotion_gates/exact_non_generic_log_transport")
            .and_then(Value::as_bool)
            != Some(false)
    {
        return Err("Round 28 isogeny certificate mismatch".into());
    }
    let roots = round297
        .pointer("/stabilizer/distinct_roots")
        .and_then(Value::as_u64)
        .ok_or("Round 297 roots missing")?;
    let folded_action = round297
        .pointer("/stabilizer/effective_folded_action_order")
        .and_then(Value::as_u64)
        .ok_or("Round 297 folded action missing")?;
    let quotient = round297
        .pointer("/stabilizer/independent_log_quotient_dimension")
        .and_then(Value::as_u64)
        .ok_or("Round 297 quotient dimension missing")?;
    if roots != 8 || folded_action != 1 || quotient != 164 {
        return Err("Round 297 stabilizer certificate mismatch".into());
    }
    let lifted_rows = round300
        .pointer("/total_lifted_rows_replayed")
        .and_then(Value::as_u64)
        .ok_or("Round 300 replay count missing")?;
    let quotient_gain = round300
        .pointer("/maximum_observed_quotient_rank_gain")
        .and_then(Value::as_u64)
        .ok_or("Round 300 quotient gain missing")?;
    if lifted_rows != 374_400
        || quotient_gain != 0
        || round300
            .pointer("/gates/quotient_identity_holds_everywhere")
            .and_then(Value::as_bool)
            != Some(true)
    {
        return Err("Round 300 cover certificate mismatch".into());
    }
    Ok(vec![
        ImportedCertificate {
            family: "CM/GLV endomorphisms".into(),
            source_round: 25,
            checked_facts: vec![
                "End(E)=Z[pi] is the maximal order (Frobenius conductor one)".into(),
                format!("minimum noninteger endomorphism degree {min_degree}"),
                "every rational-point action restricts to a known scalar".into(),
            ],
            applicable_non_generic_quotient: false,
            status: "closed within exact rational endomorphism ring".into(),
        },
        ImportedCertificate {
            family: "degree-11 isogeny class".into(),
            source_round: 28,
            checked_facts: vec![
                format!("{edges} certified edges and {models} models"),
                "all models generic-j; subgroup transport changes coordinates without log quotient".into(),
            ],
            applicable_non_generic_quotient: false,
            status: "complete registered walk; no quotient".into(),
        },
        ImportedCertificate {
            family: "scalar factor-base stabilizer".into(),
            source_round: 297,
            checked_facts: vec![
                format!("all {roots} possible eighth roots screened"),
                format!("effective folded action order {folded_action}"),
                format!("independent-log quotient dimension {quotient}"),
            ],
            applicable_non_generic_quotient: false,
            status: "only already-folded negation survives".into(),
        },
        ImportedCertificate {
            family: "cover/Jacobian fibre multiplicity".into(),
            source_round: 300,
            checked_facts: vec![
                format!("{lifted_rows} complete lifted-row replays"),
                "sheet-collapse quotient identity holds in every bucket".into(),
                format!("maximum quotient-rank gain {quotient_gain}"),
            ],
            applicable_non_generic_quotient: false,
            status: "pure kernel-sheet multiplicity closed; distinct collapsed decompositions remain open".into(),
        },
    ])
}

fn obligation(name: &str, status: &str, evidence: &str, scope: &str) -> Obligation {
    Obligation {
        name: name.into(),
        status: status.into(),
        evidence: evidence.into(),
        scope: scope.into(),
    }
}

fn build_assessment(
    invariants: &CurveInvariants,
    embedding: &EmbeddingDegreeCertificate,
    imported: &[ImportedCertificate],
) -> TransferAssessment {
    let exact_k = &embedding.exact_embedding_degree;
    let mut families = vec![
        TransferFamily {
            name: "Smart anomalous-curve transfer".into(),
            source: format!("{CURVE_SLUG}/E(F_p)"),
            destination: "additive p-adic residue representation".into(),
            transformation_type: "anomalous-curve logarithm map".into(),
            field_of_definition: "F_p with p-adic lift".into(),
            source_dimension: Some(1),
            destination_dimension: Some(1),
            relevant_subgroup: "requires #E(F_p)=p".into(),
            obligations: vec![obligation(
                "anomalous hypothesis",
                if invariants.anomalous { "supported" } else { "refuted" },
                "Exact p and n were recomputed from CurveParams::p256().",
                "P-256 rational group",
            )],
            conclusion: "P-256 is not anomalous; the required hypothesis fails.".into(),
            exploration_boundary: "Complete applicability check for the anomalous-curve family.".into(),
        },
        TransferFamily {
            name: "MOV/Frey--Ruck pairing transfer".into(),
            source: format!("{CURVE_SLUG}/E(F_p)[n]"),
            destination: format!("F_(p^k)^*, k={exact_k}"),
            transformation_type: "bilinear pairing into the n-torsion embedding field".into(),
            field_of_definition: format!("F_(p^{exact_k})"),
            source_dimension: Some(1),
            destination_dimension: None,
            relevant_subgroup: "prime order n".into(),
            obligations: vec![
                obligation("embedding degree", "supported", "Exact ord_n(p) certificate with prime-divisor minimality checks.", "P-256 n-torsion"),
                obligation("low-degree destination", if embedding.low_degree_at_most_64 { "supported" } else { "refuted" }, "The exact k is compared with the preregistered low-degree bound 64.", "registered pairing family"),
                obligation("measured advantage", "unknown", "No extension-field DLP implementation is run; the exact destination dimension already fails the low-degree applicability gate.", "end-to-end algorithm"),
            ],
            conclusion: "The exact embedding field is not a low-degree transfer and is not a smaller representation.".into(),
            exploration_boundary: "Complete embedding-degree certificate; no runtime extrapolation for finite-field DLP.".into(),
        },
        TransferFamily {
            name: "prime-field subfield descent".into(),
            source: "F_p".into(),
            destination: "proper smaller finite subfield".into(),
            transformation_type: "restriction/descent of field of definition".into(),
            field_of_definition: "F_p, p prime".into(),
            source_dimension: Some(1),
            destination_dimension: None,
            relevant_subgroup: "P-256 rational points".into(),
            obligations: vec![obligation("proper smaller subfield", "refuted", "A prime field has no proper finite subfield.", "base-field descent")],
            conclusion: "There is no smaller base field through which P-256 descends.".into(),
            exploration_boundary: "Does not reject extension-field covers or Weil restriction; those increase dimension and require separate formulas and cost.".into(),
        },
    ];
    for certificate in imported {
        families.push(TransferFamily {
            name: certificate.family.clone(),
            source: format!("{CURVE_SLUG}/prime-order subgroup"),
            destination: "family-specific representation".into(),
            transformation_type: "pinned exact certificate imported by hash".into(),
            field_of_definition: "as recorded in source round".into(),
            source_dimension: Some(1),
            destination_dimension: None,
            relevant_subgroup: "P-256 order-n subgroup".into(),
            obligations: vec![obligation(
                "non-generic quotient",
                if certificate.applicable_non_generic_quotient {
                    "supported"
                } else {
                    "refuted"
                },
                &certificate.checked_facts.join("; "),
                &format!("Round {} registered family", certificate.source_round),
            )],
            conclusion: certificate.status.clone(),
            exploration_boundary: format!(
                "Exactly the construction family and search boundary of Round {}.",
                certificate.source_round
            ),
        });
    }
    let mut assessment = TransferAssessment {
        schema: "p256-algebraic-escape-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 301,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        methodology_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/references/methodology.md (unavailable)".into(),
        template_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/assets/assessment-template.json (unavailable)".into(),
        typed_correspondence_graph: vec![
            format!("{CURVE_SLUG}/E(F_p) --[Smart; requires n=p]--> additive lift: hypothesis false"),
            format!("{CURVE_SLUG}/E(F_p)[n] --[pairing]--> F_(p^{exact_k})^*: exact high-degree destination"),
            format!("{CURVE_SLUG}/E(F_p) --[CM/isogeny/scalar action]--> cyclic image: scalar preserved"),
            format!("{CURVE_SLUG}/relations --[cover sheet pushforward]--> destination rows: kernel directions vanish"),
            "F_p --[subfield descent]--> proper smaller finite field: no such field".into(),
            "future nonstandard construction --[explicit maps plus complete cost]--> possible parity claim: open".into(),
        ],
        families,
        controls: vec![
            "All curve invariants are recomputed from CurveParams::p256().".into(),
            "The complete Round-25 factorization must multiply to n-1 before order reduction.".into(),
            "Embedding order is certified at k and falsified at k/q for every prime q dividing k.".into(),
            "Rounds 25, 28, 297, and 300 are imported only after exact byte-hash and fact checks.".into(),
        ],
        cost_accounting: vec![
            "Measured: exact modular exponentiation/multiplication/squaring counts, integer products/divisions, process telemetry, and artifact bytes.".into(),
            "Exact structural quantities: p, n, trace, discriminant, j, embedding degree, extension dimension, and imported certificate facts.".into(),
            "Unset: any unregistered future transfer, destination DLP implementation, inverse recovery, structured solver, relation yield, sparse linear algebra, and complete P-256 attack cost.".into(),
        ],
        weakest_open_obligation: "Construct a nonstandard P-256-specific map or relation mechanism outside the named families that yields distinct collapsed equations or an easier destination, then provide executable recovery and complete cost below rho.".into(),
        narrowest_supported_finding: "No applicable rho-parity route exists within the registered anomalous, supersingular/pairing, exceptional-automorphism, CM/GLV, degree-11 isogeny, scalar-stabilizer, pure cover-fibre, or base-field subfield families.".into(),
        semantic_evidence_sha256: String::new(),
        assessment_json_bytes: 0,
    };
    let semantic = json!({
        "curve": assessment.curve,
        "typed_correspondence_graph": assessment.typed_correspondence_graph,
        "families": assessment.families,
        "controls": assessment.controls,
        "cost_accounting": assessment.cost_accounting,
        "weakest_open_obligation": assessment.weakest_open_obligation,
        "narrowest_supported_finding": assessment.narrowest_supported_finding,
    });
    assessment.semantic_evidence_sha256 =
        sha256_hex(&serde_json::to_vec(&semantic).expect("assessment semantic JSON"));
    assessment
}

fn boundary_table() -> Vec<BoundaryRow> {
    vec![
        BoundaryRow {
            variant: "Pollard rho".into(),
            optimistic_ratio_to_rho: Some(1.0),
            complete_cost_measured: true,
            status: "reference".into(),
        },
        BoundaryRow {
            variant: format!("minimum-width independent-log {FACTOR_BASE_ID}"),
            optimistic_ratio_to_rho: Some(ROUND295_RATIO_TO_RHO),
            complete_cost_measured: false,
            status: "free-oracle lower boundary; named algebraic escapes do not move it".into(),
        },
        BoundaryRow {
            variant: format!("registered 17-term {REGISTERED_S17_FACTOR_BASE_ID}"),
            optimistic_ratio_to_rho: Some(REGISTERED_S17_RATIO_TO_RHO),
            complete_cost_measured: false,
            status: "free-oracle comparison; unchanged".into(),
        },
        BoundaryRow {
            variant: "unregistered future non-generic construction".into(),
            optimistic_ratio_to_rho: None,
            complete_cost_measured: false,
            status: "open / no algorithm or projection".into(),
        },
    ]
}

fn build_result(cli: &Cli) -> Result<(ResultReceipt, TransferAssessment), String> {
    let (dep25, round25) = checked_dependency(
        &cli.round25,
        ROUND25_SHA256,
        "p256.scalar_orbit_factor_base_screen/v1",
    )?;
    let (dep28, round28) = checked_dependency(
        &cli.round28,
        ROUND28_SHA256,
        "p256.isogenous_dickson_factor_base_screen/v1",
    )?;
    let (dep297, round297) = checked_dependency(
        &cli.round297,
        ROUND297_SHA256,
        "p256-union-scalar-stabilizer/v1",
    )?;
    let (dep300, round300) =
        checked_dependency(&cli.round300, ROUND300_SHA256, "p256-cover-fiber-screen/v1")?;
    for value in [&round25, &round297, &round300] {
        if value.pointer("/curve").and_then(Value::as_str) != Some(CURVE_SLUG) {
            return Err("dependency curve mismatch".into());
        }
    }
    if round28.pointer("/root_curve").and_then(Value::as_str) != Some(CURVE_SLUG) {
        return Err("Round 28 root curve mismatch".into());
    }
    let curve = CurveParams::p256();
    let mut operations = OperationCounts::default();
    let invariants = curve_invariants(&curve, &round25, &mut operations)?;
    let factors = factors_from_round25(&round25)?;
    let embedding = embedding_degree(&curve.p, &curve.n, &factors, &mut operations)?;
    let imported = imported_certificates(&round25, &round28, &round297, &round300)?;
    let named_complete = imported.len() == 4
        && !invariants.anomalous
        && !invariants.supersingular
        && !invariants.exceptional_automorphism_case
        && !embedding.low_degree_at_most_64
        && embedding.final_residue_is_one
        && embedding.all_minimality_checks_pass;
    if !named_complete {
        return Err(
            "named algebraic family census did not reach its expected exact closure".into(),
        );
    }
    let assessment = build_assessment(&invariants, &embedding, &imported);
    let gates = Gates {
        dependency_hashes_and_schemas_checked: true,
        curve_invariants_recomputed: true,
        n_minus_one_factorization_reverified: embedding.factorization_product_verified,
        exact_embedding_degree_certified: embedding.final_residue_is_one
            && embedding.all_minimality_checks_pass,
        imported_transfer_certificates_reconciled: imported.len() == 4,
        applicable_non_generic_transfer_found: false,
        explicit_forward_inverse_and_recovery_implemented: false,
        structured_residual_degree_at_most_5: false,
        parity_at_or_below_rho: false,
        per_usable_relation_below_2_103: false,
        complete_collection_below_2_120: false,
        projected_storage_below_2_50: false,
        promoted: false,
    };
    let classification =
        "standard-algebraic-transfer-census/named-families-negative/parity-blocked";
    let obstruction = "Every registered standard P-256 escape either fails its applicability hypothesis, preserves the prime-order logarithm scalar, maps only to kernel directions, or requires the exact high-degree pairing field. No named family supplies an executable easier destination or multiple distinct collapsed equations, so the 13.920747-times-rho free-oracle boundary remains unchanged.";
    let decision = "Close only the named anomalous, supersingular/pairing, exceptional-automorphism, CM/GLV, registered isogeny, scalar-stabilizer, pure cover-fibre, and base-field descent families. Do not claim universal impossibility. Require a new P-256-specific construction outside this census with explicit recovery and complete cost below rho before another full-depth relation attempt.";
    let maximum_bits = [
        curve.p.bits(),
        curve.n.bits(),
        parse_big(&embedding.exact_embedding_degree, "embedding degree")?.bits(),
        parse_big(
            &embedding.extension_field_bit_dimension,
            "extension dimension",
        )?
        .bits(),
    ]
    .into_iter()
    .max()
    .unwrap_or(0);
    let semantic = json!({
        "curve": CURVE_SLUG,
        "dependencies": [&dep25, &dep28, &dep297, &dep300],
        "curve_invariants": &invariants,
        "embedding_degree": &embedding,
        "imported_certificates": &imported,
        "operation_counts": &operations,
        "boundary_table": boundary_table(),
        "assessment": assessment.semantic_evidence_sha256,
        "gates": &gates,
        "classification": classification,
        "dominant_obstruction": obstruction,
        "decision": decision,
    });
    let semantic_evidence_sha256 =
        sha256_hex(&serde_json::to_vec(&semantic).map_err(|error| error.to_string())?);
    Ok((
        ResultReceipt {
            schema: "p256-algebraic-escape-census/v1".into(),
            curve: CURVE_SLUG.into(),
            screening_round: 301,
            execution_status: "complete".into(),
            factor_base: FACTOR_BASE_ID.into(),
            registered_s17_factor_base: REGISTERED_S17_FACTOR_BASE_ID.into(),
            dependencies: vec![dep25, dep28, dep297, dep300],
            curve_invariants: invariants,
            embedding_degree: embedding,
            imported_certificates: imported,
            operation_counts: operations,
            maximum_integer_bits: maximum_bits,
            boundary_table: boundary_table(),
            named_families_complete: named_complete,
            exploration_boundary: "Exact for the named standard algebraic families and pinned construction boundaries; no universal nonexistence claim for future algorithms.".into(),
            relations_reported_on_p256: 0,
            full_depth_unplanted_p256_relation_attempted: false,
            transfer_assessment_semantic_sha256: assessment.semantic_evidence_sha256.clone(),
            gates,
            classification: classification.into(),
            dominant_obstruction: obstruction.into(),
            decision: decision.into(),
            semantic_evidence_sha256,
            result_json_bytes: 0,
        },
        assessment,
    ))
}

fn write_result(path: &Path, result: &mut ResultReceipt) -> Result<(), String> {
    loop {
        let text = serde_json::to_string_pretty(result).map_err(|error| error.to_string())? + "\n";
        let bytes = text.len() as u64;
        if bytes == result.result_json_bytes {
            fs::write(path, text).map_err(|error| format!("{}: {error}", path.display()))?;
            return Ok(());
        }
        result.result_json_bytes = bytes;
    }
}

fn write_assessment(path: &Path, assessment: &mut TransferAssessment) -> Result<(), String> {
    loop {
        let text =
            serde_json::to_string_pretty(assessment).map_err(|error| error.to_string())? + "\n";
        let bytes = text.len() as u64;
        if bytes == assessment.assessment_json_bytes {
            fs::write(path, text).map_err(|error| format!("{}: {error}", path.display()))?;
            return Ok(());
        }
        assessment.assessment_json_bytes = bytes;
    }
}

fn run(cli: Cli) -> Result<(), String> {
    let (mut result, mut assessment) = build_result(&cli)?;
    write_assessment(&cli.assessment, &mut assessment)?;
    write_result(&cli.out, &mut result)?;
    eprintln!(
        "round 301: embedding degree {} ({} bits), named families complete={}, promoted={}",
        result.embedding_degree.exact_embedding_degree,
        result.embedding_degree.exact_embedding_degree_bits,
        result.named_families_complete,
        result.gates.promoted
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_algebraic_escape_census: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn p256_invariants_are_nonexceptional() {
        let curve = CurveParams::p256();
        let mut operations = OperationCounts::default();
        let j = modular_j(&curve, &mut operations).expect("j");
        assert!(!j.is_zero());
        assert_ne!(j, BigUint::from(1728u16));
        assert_ne!(curve.p, curve.n);
        assert_eq!(curve.h, 1);
    }

    #[test]
    fn counted_modpow_matches_library() {
        let modulus = BigUint::from(149u16);
        let base = BigUint::from(37u8);
        let exponent = BigUint::from(123u8);
        let mut operations = OperationCounts::default();
        assert_eq!(
            modpow_counted(&base, &exponent, &modulus, &mut operations),
            base.modpow(&exponent, &modulus)
        );
        assert_eq!(operations.modular_exponentiations, 1);
        assert_ne!(operations.modular_squarings, 0);
    }

    #[test]
    fn small_order_certificate_reduces_exactly() {
        let p = BigUint::from(2u8);
        let n = BigUint::from(31u8);
        let factors = vec![
            (BigUint::from(2u8), 1),
            (BigUint::from(3u8), 1),
            (BigUint::from(5u8), 1),
        ];
        let mut operations = OperationCounts::default();
        let certificate = embedding_degree(&p, &n, &factors, &mut operations).expect("order");
        assert_eq!(certificate.exact_embedding_degree, "5");
        assert!(certificate.final_residue_is_one);
        assert!(certificate.all_minimality_checks_pass);
    }
}
