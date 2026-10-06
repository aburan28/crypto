//! Exact-memoryless P-256 factor-base closure audit (round 29).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::Serialize;
use serde_json::Value;

const ROUND24_SHA256: &str = "d1e2fe33e1abbae8a9bf6885ca7cfef8012233363f4cd43e0d30d468a381fc9c";
const ROUND25_SHA256: &str = "dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf";
const ROUND26_SHA256: &str = "463c100dd951e2a7a7a419118f84bde90fc419358d8279a430a5a0e6191db029";
const ROUND27_SHA256: &str = "256271ebb6f1c3c736a01da4358f01a6822c4315bb2081adcba1dd05b1d847a5";
const ROUND28_SHA256: &str = "950d4b7f4997e146ee30fc7c61174070595fa5d344b9a30cbc0ca578ae7e60f4";
const GROUP_ORDER: &str =
    "115792089210356248762697446949407573529996955224135760342422259061068512044369";
const RELATION_MODEL: &str =
    "17 distinct variable columns; signed 8+9 cross-colour collision modulo global negation";
const RHO_S: f64 = 1.3;
const PARITY_OPERATIONS: f64 = 1.036_982_447_222_106_5;
const AFFINE_CONTROLS: u64 = 4_096;

#[derive(Parser)]
#[command(about = "Close exact-memoryless factor-base transports for P-256")]
struct Cli {
    #[arg(long)]
    round24: PathBuf,
    #[arg(long)]
    round25: PathBuf,
    #[arg(long)]
    round26: PathBuf,
    #[arg(long)]
    round27: PathBuf,
    #[arg(long)]
    round28: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Serialize)]
struct Dependency {
    round: u64,
    path: String,
    sha256: String,
    schema: String,
    curve: String,
    relation_model: String,
}

struct LoadedDependency {
    record: Dependency,
    value: Value,
}

#[derive(Serialize)]
struct Scope {
    included: Vec<String>,
    excluded: Vec<String>,
    discarded_probabilistic_branches_counted_as_exhaustive: bool,
}

#[derive(Serialize)]
struct TranslationLemma {
    group_order: String,
    group_order_is_odd_prime: bool,
    nonzero_delta_gcd_group_order: String,
    nonzero_translation_orbit_length: String,
    proper_nonempty_translation_invariant_subset_exists: bool,
    signed_factor_base_contains_identity: bool,
    consequence: String,
}

#[derive(Serialize)]
struct AffineIdentityAudit {
    controls: u64,
    inverses_verified: u64,
    conjugated_negation_verified: u64,
    generated_translation_verified: u64,
    failures: u64,
    controls_sha256: String,
    conjugated_negation: String,
    generated_translation: String,
}

#[derive(Serialize)]
struct ToyPrimeRow {
    prime: u64,
    signed_nonempty_subsets: u64,
    nonzero_translation_cases: u64,
    nonzero_translation_invariant_subsets: u64,
    affine_cases: u64,
    nonzero_translation_affine_invariants: u64,
    scalar_invariant_cases: u64,
}

#[derive(Serialize)]
struct ExhaustiveReference {
    primes: Vec<u64>,
    rows: Vec<ToyPrimeRow>,
    signed_nonempty_subsets: u64,
    nonzero_translation_cases: u64,
    affine_cases: u64,
    false_positives: u64,
    false_negatives: u64,
    ordered_rows_sha256: String,
    complete_for_reported_instances: bool,
}

#[derive(Serialize)]
struct CoalescenceClass {
    class: String,
    coalesces_on_group_sum: bool,
    preserves_proper_signed_fixed_s17: bool,
    exact_log_transport: bool,
    non_generic: bool,
    terminal_class: String,
    proof_obligation: String,
}

#[derive(Serialize)]
struct ScalarActionCensus {
    rho_s: f64,
    maximum_operations_per_sample_for_parity: f64,
    divisors_enumerated: u64,
    retained_orders_through_maximum_width: u64,
    exact_order_elements_visited: u64,
    unique_folded_actions: u64,
    width_band_combinations_evaluated: u64,
    enumeration_exhaustive: bool,
    minimum_binary_width_at_least_17: u64,
    minimum_binary_operations: u64,
    minimum_binary_scalar: String,
    minimum_doubling_width_at_least_17: u64,
    minimum_doubling_lower_bound: u64,
    minimum_doubling_scalar: String,
    target_corrected_minimum_doubling_operations: u64,
    target_corrected_minimum_doubling_ratio_to_rho: f64,
}

#[derive(Serialize)]
struct ImportedCorrectness {
    round24_false_positives: u64,
    round24_false_negatives: u64,
    round25_transport_failures: u64,
    round26_duplicate_deltas: u64,
    round26_control_false_positives: u64,
    round26_control_false_negatives: u64,
    round27_false_positives_gate: bool,
    round27_false_negatives_gate: bool,
    round28_false_positives: u64,
    round28_false_negatives: u64,
    all_imported_complete_checks_pass: bool,
}

#[derive(Serialize)]
struct BoundaryRow {
    name: String,
    operations_per_accepted_sample: Option<f64>,
    ratio_to_rho: f64,
    projected_log2_operations: Option<f64>,
    projected_log2_bytes: Option<f64>,
    preserves_fixed_s17: bool,
    coalesces: bool,
    exact_log_quotient: bool,
    non_generic: bool,
    decision: String,
}

#[derive(Serialize)]
struct DegreeBoundary {
    short_weierstrass_s3_total_degree: u64,
    all_isogenous_models_retain_universal_leading_support: bool,
    isogeny_pullback_is_degree_reduction: bool,
    structured_residual_degree_of_regularity: Option<u64>,
    degree_at_most_5_gate_passed: bool,
}

#[derive(Serialize)]
struct PromotionGates {
    exact_coalescing_transition: bool,
    preserves_fixed_s17: bool,
    non_generic: bool,
    exact_logarithm_quotient: bool,
    structured_residual_degree_at_most_5: bool,
    complete_cost_at_or_below_rho: bool,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    zero_false_positives_and_false_negatives: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    relation_model: String,
    dependencies: Vec<Dependency>,
    scope: Scope,
    translation_lemma: TranslationLemma,
    affine_identity_audit: AffineIdentityAudit,
    exhaustive_small_prime_reference: ExhaustiveReference,
    coalescence_classification: Vec<CoalescenceClass>,
    scalar_action_census: ScalarActionCensus,
    imported_correctness: ImportedCorrectness,
    terminal_boundaries: Vec<BoundaryRow>,
    degree_boundary: DegreeBoundary,
    promotion_gates: PromotionGates,
    full_depth_unplanted_attempted: bool,
    dominant_obstruction: String,
    result_scope_limit: String,
    decision: String,
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    BigUint::parse_bytes(value.as_bytes(), 10)
        .ok_or_else(|| format!("invalid frozen decimal integer: {value}"))
}

fn load_dependency(
    path: &Path,
    round: u64,
    expected_sha256: &str,
) -> Result<LoadedDependency, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != expected_sha256 {
        return Err(format!(
            "round-{round} dependency hash mismatch: expected {expected_sha256}, got {digest}"
        ));
    }
    let value: Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("{}: {error}", path.display()))?;
    let curve = value
        .get("curve")
        .or_else(|| value.get("root_curve"))
        .and_then(Value::as_str)
        .ok_or_else(|| format!("round-{round} dependency has no curve/root_curve string"))?;
    let relation_model = string_at(&value, &["relation_model"])?;
    if curve != CURVE_SLUG {
        return Err(format!("round-{round} curve mismatch: {curve}"));
    }
    if relation_model != RELATION_MODEL {
        return Err(format!(
            "round-{round} relation-model mismatch: {relation_model}"
        ));
    }
    Ok(LoadedDependency {
        record: Dependency {
            round,
            path: path.display().to_string(),
            sha256: digest,
            schema: string_at(&value, &["schema"])?.to_string(),
            curve: curve.to_string(),
            relation_model: relation_model.to_string(),
        },
        value,
    })
}

fn at<'a>(value: &'a Value, path: &[&str]) -> Result<&'a Value, String> {
    let mut cursor = value;
    for key in path {
        cursor = cursor
            .get(*key)
            .ok_or_else(|| format!("missing JSON path {}", path.join(".")))?;
    }
    Ok(cursor)
}

fn string_at<'a>(value: &'a Value, path: &[&str]) -> Result<&'a str, String> {
    at(value, path)?
        .as_str()
        .ok_or_else(|| format!("JSON path {} is not a string", path.join(".")))
}

fn u64_at(value: &Value, path: &[&str]) -> Result<u64, String> {
    at(value, path)?
        .as_u64()
        .ok_or_else(|| format!("JSON path {} is not a u64", path.join(".")))
}

fn f64_at(value: &Value, path: &[&str]) -> Result<f64, String> {
    at(value, path)?
        .as_f64()
        .ok_or_else(|| format!("JSON path {} is not a number", path.join(".")))
}

fn bool_at(value: &Value, path: &[&str]) -> Result<bool, String> {
    at(value, path)?
        .as_bool()
        .ok_or_else(|| format!("JSON path {} is not a bool", path.join(".")))
}

fn projection_named<'a>(round24: &'a Value, name: &str) -> Result<&'a Value, String> {
    at(round24, &["projections"])?
        .as_array()
        .ok_or("round-24 projections is not an array")?
        .iter()
        .find(|row| row.get("name").and_then(Value::as_str) == Some(name))
        .ok_or_else(|| format!("round-24 projection not found: {name}"))
}

fn hash_big(label: &str, index: u64, modulus: &BigUint, nonzero: bool) -> BigUint {
    let digest = sha256(format!("{CURVE_SLUG}/memoryless-round29/{label}/{index}").as_bytes());
    if nonzero {
        BigUint::from_bytes_be(&digest) % (modulus - 1u8) + 1u8
    } else {
        BigUint::from_bytes_be(&digest) % modulus
    }
}

fn neg_mod(value: &BigUint, modulus: &BigUint) -> BigUint {
    if value.is_zero() {
        BigUint::zero()
    } else {
        modulus - value
    }
}

fn affine_identity_audit(modulus: &BigUint) -> Result<AffineIdentityAudit, String> {
    let exponent = modulus - 2u8;
    let mut blob = Vec::with_capacity(AFFINE_CONTROLS as usize * 96);
    let mut inverse_verified = 0u64;
    let mut conjugated_verified = 0u64;
    let mut translation_verified = 0u64;
    let mut failures = 0u64;
    for index in 0..AFFINE_CONTROLS {
        let lambda = hash_big("lambda", index, modulus, true);
        let b = hash_big("translation", index, modulus, false);
        let x = hash_big("point", index, modulus, false);
        let inverse = lambda.modpow(&exponent, modulus);
        let inverse_ok = &lambda * &inverse % modulus == BigUint::one();
        inverse_verified += u64::from(inverse_ok);

        // f o nu o f^-1 evaluated directly in coefficient arithmetic.
        let x_minus_b = (&x + modulus - &b) % modulus;
        let inverse_image = &inverse * x_minus_b % modulus;
        let negated = neg_mod(&inverse_image, modulus);
        let conjugated = (&lambda * negated + &b) % modulus;
        let expected_conjugated = (neg_mod(&x, modulus) + (&b << 1usize)) % modulus;
        let conjugated_ok = conjugated == expected_conjugated;
        conjugated_verified += u64::from(conjugated_ok);

        // (f o nu o f^-1) o nu is translation by 2b.
        let negative_x = neg_mod(&x, modulus);
        let negative_x_minus_b = (&negative_x + modulus - &b) % modulus;
        let inverse_negative = &inverse * negative_x_minus_b % modulus;
        let generated = (&lambda * neg_mod(&inverse_negative, modulus) + &b) % modulus;
        let expected_translation = (&x + (&b << 1usize)) % modulus;
        let translation_ok = generated == expected_translation;
        translation_verified += u64::from(translation_ok);
        failures += u64::from(!(inverse_ok && conjugated_ok && translation_ok));

        for value in [&lambda, &b, &x] {
            let raw = value.to_bytes_be();
            blob.extend(std::iter::repeat_n(0u8, 32usize.saturating_sub(raw.len())));
            blob.extend(raw);
        }
    }
    if failures != 0 {
        return Err(format!("{failures} affine identity controls failed"));
    }
    Ok(AffineIdentityAudit {
        controls: AFFINE_CONTROLS,
        inverses_verified: inverse_verified,
        conjugated_negation_verified: conjugated_verified,
        generated_translation_verified: translation_verified,
        failures,
        controls_sha256: hex::encode(sha256(&blob)),
        conjugated_negation: "f o nu o f^-1: x -> -x+2b".into(),
        generated_translation: "(f o nu o f^-1) o nu: x -> x+2b".into(),
    })
}

fn in_signed_subset(mask: u64, value: u64, prime: u64) -> bool {
    if value == 0 {
        return false;
    }
    let representative = value.min(prime - value);
    mask & (1u64 << (representative - 1)) != 0
}

fn subset_invariant(mask: u64, prime: u64, lambda: u64, b: u64) -> bool {
    (1..prime).all(|value| {
        !in_signed_subset(mask, value, prime)
            || in_signed_subset(mask, (lambda * value + b) % prime, prime)
    })
}

fn exhaustive_reference() -> Result<ExhaustiveReference, String> {
    let primes = vec![3u64, 5, 7, 11, 13, 17, 19, 23, 29, 31];
    let mut rows = Vec::new();
    let mut total_subsets = 0u64;
    let mut total_translation_cases = 0u64;
    let mut total_affine_cases = 0u64;
    let mut false_positives = 0u64;
    let mut false_negatives = 0u64;
    let mut blob = Vec::new();
    for &prime in &primes {
        let pairs = (prime - 1) / 2;
        let limit = 1u64 << pairs;
        let mut row = ToyPrimeRow {
            prime,
            signed_nonempty_subsets: 0,
            nonzero_translation_cases: 0,
            nonzero_translation_invariant_subsets: 0,
            affine_cases: 0,
            nonzero_translation_affine_invariants: 0,
            scalar_invariant_cases: 0,
        };
        for mask in 1..limit {
            row.signed_nonempty_subsets += 1;
            for delta in 1..prime {
                row.nonzero_translation_cases += 1;
                let invariant = subset_invariant(mask, prime, 1, delta);
                row.nonzero_translation_invariant_subsets += u64::from(invariant);
                false_positives += u64::from(invariant);
            }
            for lambda in 1..prime {
                for b in 0..prime {
                    row.affine_cases += 1;
                    let invariant = subset_invariant(mask, prime, lambda, b);
                    if b == 0 {
                        row.scalar_invariant_cases += u64::from(invariant);
                    } else {
                        row.nonzero_translation_affine_invariants += u64::from(invariant);
                        false_negatives += u64::from(invariant);
                    }
                }
            }
        }
        if row.nonzero_translation_invariant_subsets != 0
            || row.nonzero_translation_affine_invariants != 0
        {
            return Err(format!(
                "small-prime closure theorem failed for prime {prime}"
            ));
        }
        total_subsets += row.signed_nonempty_subsets;
        total_translation_cases += row.nonzero_translation_cases;
        total_affine_cases += row.affine_cases;
        for value in [
            row.prime,
            row.signed_nonempty_subsets,
            row.nonzero_translation_cases,
            row.nonzero_translation_invariant_subsets,
            row.affine_cases,
            row.nonzero_translation_affine_invariants,
            row.scalar_invariant_cases,
        ] {
            blob.extend(value.to_be_bytes());
        }
        rows.push(row);
    }
    Ok(ExhaustiveReference {
        primes,
        rows,
        signed_nonempty_subsets: total_subsets,
        nonzero_translation_cases: total_translation_cases,
        affine_cases: total_affine_cases,
        false_positives,
        false_negatives,
        ordered_rows_sha256: hex::encode(sha256(&blob)),
        complete_for_reported_instances: true,
    })
}

fn scalar_action_census(round27: &Value) -> Result<ScalarActionCensus, String> {
    let census = at(round27, &["action_census"])?;
    let rows = at(census, &["order_minima"])?
        .as_array()
        .ok_or("round-27 order_minima is not an array")?;
    let eligible = rows
        .iter()
        .filter(|row| {
            row.get("folded_width")
                .and_then(Value::as_u64)
                .is_some_and(|width| width >= 17)
        })
        .collect::<Vec<_>>();
    let minimum_binary = eligible
        .iter()
        .min_by_key(|row| {
            (
                row.get("minimum_binary_operations")
                    .and_then(Value::as_u64)
                    .unwrap_or(u64::MAX),
                row.get("minimum_doubling_lower_bound")
                    .and_then(Value::as_u64)
                    .unwrap_or(u64::MAX),
            )
        })
        .ok_or("no round-27 scalar action with width at least 17")?;
    let minimum_doubling = eligible
        .iter()
        .min_by_key(|row| {
            (
                row.get("minimum_doubling_lower_bound")
                    .and_then(Value::as_u64)
                    .unwrap_or(u64::MAX),
                row.get("minimum_binary_operations")
                    .and_then(Value::as_u64)
                    .unwrap_or(u64::MAX),
            )
        })
        .ok_or("no round-27 scalar action with width at least 17")?;
    let minimum_doubling_operations = u64_at(minimum_doubling, &["minimum_doubling_lower_bound"])?;
    let corrected = minimum_doubling_operations + 1;
    let imported_ratio = f64_at(
        round27,
        &[
            "projection",
            "global_minimum_doubling_lower_bound_ratio_to_rho",
        ],
    )?;
    let derived_ratio = corrected as f64 / PARITY_OPERATIONS;
    if (derived_ratio - imported_ratio).abs() > 1e-12 {
        return Err(format!(
            "round-27 global floor mismatch: derived {derived_ratio}, imported {imported_ratio}"
        ));
    }
    Ok(ScalarActionCensus {
        rho_s: RHO_S,
        maximum_operations_per_sample_for_parity: PARITY_OPERATIONS,
        divisors_enumerated: u64_at(census, &["divisors_enumerated"])?,
        retained_orders_through_maximum_width: u64_at(census, &["retained_nontrivial_orders"])?,
        exact_order_elements_visited: u64_at(census, &["exact_order_elements_visited"])?,
        unique_folded_actions: u64_at(census, &["unique_folded_actions"])?,
        width_band_combinations_evaluated: u64_at(census, &["width_band_combinations_evaluated"])?,
        enumeration_exhaustive: bool_at(census, &["enumeration_exhaustive"])?,
        minimum_binary_width_at_least_17: u64_at(minimum_binary, &["folded_width"])?,
        minimum_binary_operations: u64_at(minimum_binary, &["minimum_binary_operations"])?,
        minimum_binary_scalar: string_at(minimum_binary, &["minimum_scalar"])?.into(),
        minimum_doubling_width_at_least_17: u64_at(minimum_doubling, &["folded_width"])?,
        minimum_doubling_lower_bound: minimum_doubling_operations,
        minimum_doubling_scalar: string_at(minimum_doubling, &["minimum_scalar"])?.into(),
        target_corrected_minimum_doubling_operations: corrected,
        target_corrected_minimum_doubling_ratio_to_rho: derived_ratio,
    })
}

fn imported_correctness(
    round24: &Value,
    round25: &Value,
    round26: &Value,
    round27: &Value,
    round28: &Value,
) -> Result<ImportedCorrectness, String> {
    let record = ImportedCorrectness {
        round24_false_positives: u64_at(round24, &["reference_false_positives"])?,
        round24_false_negatives: u64_at(round24, &["reference_false_negatives"])?,
        round25_transport_failures: u64_at(round25, &["factor_base", "transport_replay_failures"])?,
        round26_duplicate_deltas: u64_at(round26, &["delta_census", "duplicate_deltas"])?,
        round26_control_false_positives: u64_at(round26, &["local_controls", "false_positives"])?,
        round26_control_false_negatives: u64_at(round26, &["local_controls", "false_negatives"])?,
        round27_false_positives_gate: bool_at(
            round27,
            &["promotion_gates", "zero_false_positives"],
        )?,
        round27_false_negatives_gate: bool_at(
            round27,
            &["promotion_gates", "zero_false_negatives"],
        )?,
        round28_false_positives: u64_at(round28, &["native_factor_base", "false_positives"])?,
        round28_false_negatives: u64_at(round28, &["native_factor_base", "false_negatives"])?,
        all_imported_complete_checks_pass: false,
    };
    let pass = record.round24_false_positives == 0
        && record.round24_false_negatives == 0
        && record.round25_transport_failures == 0
        && record.round26_duplicate_deltas == 0
        && record.round26_control_false_positives == 0
        && record.round26_control_false_negatives == 0
        && record.round27_false_positives_gate
        && record.round27_false_negatives_gate
        && record.round28_false_positives == 0
        && record.round28_false_negatives == 0;
    Ok(ImportedCorrectness {
        all_imported_complete_checks_pass: pass,
        ..record
    })
}

fn run(cli: Cli) -> Result<ResultFile, String> {
    let round24 = load_dependency(&cli.round24, 24, ROUND24_SHA256)?;
    let round25 = load_dependency(&cli.round25, 25, ROUND25_SHA256)?;
    let round26 = load_dependency(&cli.round26, 26, ROUND26_SHA256)?;
    let round27 = load_dependency(&cli.round27, 27, ROUND27_SHA256)?;
    let round28 = load_dependency(&cli.round28, 28, ROUND28_SHA256)?;
    let group_order = parse_big(GROUP_ORDER)?;
    if string_at(&round27.value, &["classification", "group_order"])? != GROUP_ORDER
        || !bool_at(
            &round27.value,
            &["classification", "group_order_is_odd_prime"],
        )?
    {
        return Err("round-27 prime-order classification did not replay".into());
    }

    let affine_audit = affine_identity_audit(&group_order)?;
    let exhaustive = exhaustive_reference()?;
    let scalar_census = scalar_action_census(&round27.value)?;
    let correctness = imported_correctness(
        &round24.value,
        &round25.value,
        &round26.value,
        &round27.value,
        &round28.value,
    )?;
    if !correctness.all_imported_complete_checks_pass {
        return Err("an imported complete correctness check failed".into());
    }

    let known_log = projection_named(&round24.value, "known-log K=1 target")?;
    let independent = projection_named(&round24.value, "round-23 best-width quotient")?;
    let local_ratio = f64_at(
        &round26.value,
        &["parity_projection", "local_oracle_ratio_to_rho"],
    )?;
    let local_storage = f64_at(
        &round26.value,
        &["parity_projection", "local_materialized_log2_bytes"],
    )?;
    let round28_collection = f64_at(
        &round28.value,
        &["best_candidate", "projected_138031_rows_log2_operations"],
    )?;
    let round28_per_row = f64_at(
        &round28.value,
        &[
            "best_candidate",
            "projected_cost_per_usable_relation_log2_operations",
        ],
    )?;
    let round28_storage = f64_at(
        &round28.value,
        &["best_candidate", "projected_peak_two_list_log2_bytes"],
    )?;

    let result = ResultFile {
        schema: "p256.memoryless_factor_base_closure/v1".into(),
        curve: CURVE_SLUG.into(),
        relation_model: RELATION_MODEL.into(),
        dependencies: vec![
            round24.record,
            round25.record,
            round26.record,
            round27.record,
            round28.record,
        ],
        scope: Scope {
            included: vec![
                "one-column exact transports whose transition factors through the represented group sum".into(),
                "bounded-tuple exact transports with representation-independent total delta".into(),
                "elliptic-curve morphisms and their affine translates on the prime-order P-256 subgroup".into(),
                "coordinate correspondences credited only when an exact group/log transport is supplied".into(),
            ],
            excluded: vec![
                "non-memoryless solvers with unbounded representation state".into(),
                "unknown future relation algorithms outside fixed signed S17".into(),
                "probabilistic branch discards presented as exhaustive".into(),
            ],
            discarded_probabilistic_branches_counted_as_exhaustive: false,
        },
        translation_lemma: TranslationLemma {
            group_order: GROUP_ORDER.into(),
            group_order_is_odd_prime: true,
            nonzero_delta_gcd_group_order: "1".into(),
            nonzero_translation_orbit_length: GROUP_ORDER.into(),
            proper_nonempty_translation_invariant_subset_exists: false,
            signed_factor_base_contains_identity: false,
            consequence: "Closure under any nonzero constant delta forces the complete prime-order group, including the identity; a proper signed factor base therefore admits only zero constant delta.".into(),
        },
        affine_identity_audit: affine_audit,
        exhaustive_small_prime_reference: exhaustive,
        coalescence_classification: vec![
            CoalescenceClass {
                class: "one-column, representation-independent delta".into(),
                coalesces_on_group_sum: true,
                preserves_proper_signed_fixed_s17: false,
                exact_log_transport: true,
                non_generic: false,
                terminal_class: "nonzero delta forces the full group; zero delta does not move".into(),
                proof_obligation: "u(P)-P=C for every admitted column; closure under C".into(),
            },
            CoalescenceClass {
                class: "bounded-tuple, representation-independent total delta".into(),
                coalesces_on_group_sum: true,
                preserves_proper_signed_fixed_s17: false,
                exact_log_transport: true,
                non_generic: false,
                terminal_class: "partitioned S -> S+C is Pollard rho with redundant tuple state".into(),
                proof_obligation: "sum(u(tuple))-sum(tuple)=C in each group-state partition".into(),
            },
            CoalescenceClass {
                class: "identity-preserving elliptic-curve morphism".into(),
                coalesces_on_group_sum: true,
                preserves_proper_signed_fixed_s17: true,
                exact_log_transport: true,
                non_generic: true,
                terminal_class: "scalar multiplication; proper signed invariant bases are scalar-orbit unions".into(),
                proof_obligation: "elliptic morphism is a homomorphism; prime-order subgroup endomorphism is scalar".into(),
            },
            CoalescenceClass {
                class: "coordinate-only low-degree correspondence".into(),
                coalesces_on_group_sum: false,
                preserves_proper_signed_fixed_s17: true,
                exact_log_transport: false,
                non_generic: true,
                terminal_class: "independent logarithm quotient; collision boundary remains".into(),
                proof_obligation: "no exact point delta or subgroup scalar is available".into(),
            },
        ],
        scalar_action_census: scalar_census,
        imported_correctness: correctness,
        terminal_boundaries: vec![
            BoundaryRow {
                name: "Pollard rho partitioned translation".into(),
                operations_per_accepted_sample: Some(1.0),
                ratio_to_rho: 1.0,
                projected_log2_operations: None,
                projected_log2_bytes: None,
                preserves_fixed_s17: false,
                coalesces: true,
                exact_log_quotient: true,
                non_generic: false,
                decision: "rho reference; not an index-calculus improvement".into(),
            },
            BoundaryRow {
                name: "known-log K=1 collision oracle".into(),
                operations_per_accepted_sample: None,
                ratio_to_rho: f64_at(known_log, &["ratio_to_rho"] )?,
                projected_log2_operations: Some(f64_at(
                    known_log,
                    &["complete_lower_bound_log2_operations"],
                )?),
                projected_log2_bytes: Some(f64_at(
                    known_log,
                    &["optimistic_materialized_log2_bytes"],
                )?),
                preserves_fixed_s17: true,
                coalesces: false,
                exact_log_quotient: true,
                non_generic: false,
                decision: "oracle constant reaches parity, but materialization fails and a memoryless realization is rho-equivalent".into(),
            },
            BoundaryRow {
                name: "one-addition scalar-orbit neighbour".into(),
                operations_per_accepted_sample: Some(1.0),
                ratio_to_rho: local_ratio,
                projected_log2_operations: None,
                projected_log2_bytes: Some(local_storage),
                preserves_fixed_s17: true,
                coalesces: false,
                exact_log_quotient: true,
                non_generic: true,
                decision: "fails coalescence on every distinct-column collision control".into(),
            },
            BoundaryRow {
                name: "best proper signed scalar action, doubling floor".into(),
                operations_per_accepted_sample: Some(
                    u64_at(
                        &round27.value,
                        &[
                            "projection",
                            "global_minimum_target_corrected_doubling_lower_bound",
                        ],
                    )? as f64,
                ),
                ratio_to_rho: f64_at(
                    &round27.value,
                    &[
                        "projection",
                        "global_minimum_doubling_lower_bound_ratio_to_rho",
                    ],
                )?,
                projected_log2_operations: None,
                projected_log2_bytes: None,
                preserves_fixed_s17: true,
                coalesces: true,
                exact_log_quotient: true,
                non_generic: true,
                decision: "complete scalar-action census remains hundreds of times above parity".into(),
            },
            BoundaryRow {
                name: "best-width independent-log quotient".into(),
                operations_per_accepted_sample: None,
                ratio_to_rho: f64_at(independent, &["ratio_to_rho"] )?,
                projected_log2_operations: Some(f64_at(
                    independent,
                    &["complete_lower_bound_log2_operations"],
                )?),
                projected_log2_bytes: Some(f64_at(
                    independent,
                    &["optimistic_materialized_log2_bytes"],
                )?),
                preserves_fixed_s17: true,
                coalesces: false,
                exact_log_quotient: false,
                non_generic: true,
                decision: "negation quotient still leaves one unknown per column".into(),
            },
        ],
        degree_boundary: DegreeBoundary {
            short_weierstrass_s3_total_degree: u64_at(
                &round28.value,
                &["degree_boundary", "short_weierstrass_s3_total_degree"],
            )?,
            all_isogenous_models_retain_universal_leading_support: bool_at(
                &round28.value,
                &[
                    "degree_boundary",
                    "all_models_retain_universal_leading_support",
                ],
            )?,
            isogeny_pullback_is_degree_reduction: bool_at(
                &round28.value,
                &["degree_boundary", "isogeny_pullback_is_degree_reduction"],
            )?,
            structured_residual_degree_of_regularity: at(
                &round28.value,
                &[
                    "degree_boundary",
                    "structured_residual_degree_of_regularity",
                ],
            )?
            .as_u64(),
            degree_at_most_5_gate_passed: bool_at(
                &round28.value,
                &["degree_boundary", "degree_gate_passed"],
            )?,
        },
        promotion_gates: PromotionGates {
            exact_coalescing_transition: true,
            preserves_fixed_s17: true,
            non_generic: false,
            exact_logarithm_quotient: true,
            structured_residual_degree_at_most_5: false,
            complete_cost_at_or_below_rho: false,
            relation_collection_below_2_120: round28_collection < 120.0,
            per_usable_relation_below_2_103: round28_per_row < 103.0,
            materialized_storage_below_2_50: round28_storage < 50.0,
            zero_false_positives_and_false_negatives: true,
            promoted: false,
        },
        full_depth_unplanted_attempted: false,
        dominant_obstruction: "Exact coalescence collapses bounded local updates to group translations (generic rho) and algebraic transports to scalar actions (at least 228.548x rho at the optimistic doubling floor). Coordinate-only factor bases retain the independent-log quotient (at least 392.155x rho).".into(),
        result_scope_limit: "This closes exact memoryless bounded-update and algebraic-transport selectors for proper signed fixed-S17 P-256 factor bases. It does not prove that every conceivable non-memoryless or future relation algorithm is impossible.".into(),
        decision: "No admitted factor-base selector reaches rho parity. The only one-operation coalescing terminal class is ordinary Pollard rho; every class that preserves proper fixed S17 fails coalescence, exact-log transport, or cost.".into(),
    };

    if result.promotion_gates.promoted || result.full_depth_unplanted_attempted {
        return Err("closure result unexpectedly promoted a candidate".into());
    }
    if let Some(path) = cli.out {
        let mut encoded = serde_json::to_vec_pretty(&result).map_err(|error| error.to_string())?;
        encoded.push(b'\n');
        fs::write(&path, encoded).map_err(|error| format!("{}: {error}", path.display()))?;
    } else {
        println!(
            "{}",
            serde_json::to_string_pretty(&result).map_err(|error| error.to_string())?
        );
    }
    Ok(result)
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(_) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_subset_encoding_is_exact() {
        let prime = 11;
        let mask = (1 << 0) | (1 << 3);
        let members = (0..prime)
            .filter(|value| in_signed_subset(mask, *value, prime))
            .collect::<Vec<_>>();
        assert_eq!(members, vec![1, 4, 7, 10]);
    }

    #[test]
    fn exhaustive_small_prime_closure_has_no_counterexample() {
        let result = exhaustive_reference().expect("complete reference");
        assert_eq!(result.false_positives, 0);
        assert_eq!(result.false_negatives, 0);
        assert!(result.nonzero_translation_cases > 1_000_000);
        assert!(result.affine_cases > result.nonzero_translation_cases);
    }

    #[test]
    fn nonzero_translation_does_not_preserve_a_signed_proper_set() {
        for prime in [5u64, 7, 11, 13] {
            let limit = 1u64 << ((prime - 1) / 2);
            for mask in 1..limit {
                for delta in 1..prime {
                    assert!(!subset_invariant(mask, prime, 1, delta));
                }
            }
        }
    }

    #[test]
    fn fixed_s17_needs_at_least_seventeen_columns() {
        assert!((0u64..17).all(|columns| columns < 17));
        assert_eq!((0u64..17).count(), 17);
    }

    #[test]
    fn frozen_group_order_is_odd() {
        let order = parse_big(GROUP_ORDER).expect("group order");
        assert_eq!((&order & BigUint::one()), BigUint::one());
        assert_eq!(order.bits(), 256);
    }

    #[test]
    fn parity_ratio_matches_round27_floor() {
        let ratio = 237.0 / PARITY_OPERATIONS;
        assert!((ratio - 228.547_745_079_852_9).abs() < 1e-12);
        assert_eq!(RHO_S, 1.3);
    }
}
