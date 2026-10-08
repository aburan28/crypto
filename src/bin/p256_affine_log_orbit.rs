//! Exact affine-log-orbit factor-base screen for P-256 (round 39).

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::ct_bignum::U256;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::Serialize;
use serde_json::{json, Value};

const CURVE_SLUG: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const COMPARISON_FB: &str = "FB1h2f8621cda105";
const ROUND19_SHA256: &str = "3096540621408e4a48cfa18963ad01da9686d3527ee26776c8b6cf6f45a71114";
const ROUND25_SHA256: &str = "dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf";
const COLUMNS: u64 = 131_458;
const RELATION_TERMS: usize = 17;
const POSITIVE_TERMS: usize = 8;
const P256_INSTANCES: u64 = 4_096;
const P256_PAIR_CONTROLS: u64 = 4_096;

#[derive(Parser)]
#[command(about = "Screen affine-log-orbit P-256 factor bases")]
struct Cli {
    #[arg(long)]
    round19: PathBuf,
    #[arg(long)]
    round25: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Clone, Serialize)]
struct Dependency {
    round: u64,
    path: String,
    sha256: String,
    schema: String,
}

#[derive(Clone, Serialize)]
struct ImportedDegreeBoundary {
    factor_base: String,
    structured_residual_maximum: u64,
    every_component_complete_and_correct: bool,
    applies_to_affine_log_family: bool,
    unsplit_s17_degree_of_regularity: Option<u64>,
}

#[derive(Clone, Serialize)]
struct ToyCell {
    modulus: u64,
    family: String,
    width: u64,
    signed_terms: u64,
    positive_terms: u64,
    coefficient_patterns: u64,
    candidates_checked: u64,
    identity_true: u64,
    impossible: u64,
    informative_true: u64,
    informative_false: u64,
    modular_inversions: u64,
    anchor_recovery_failures: u64,
    pair_difference_checks: u64,
    pair_difference_failures: u64,
    false_positives: u64,
    false_negatives: u64,
}

#[derive(Clone, Default, Serialize)]
struct GroupOperations {
    additions: u64,
    doublings: u64,
    scalar_multiplications: u64,
}

#[derive(Clone, Serialize)]
struct P256Controls {
    deterministic_instances: u64,
    identity_relations: u64,
    informative_relations: u64,
    impossible_negative_controls: u64,
    exact_positive_group_replays: u64,
    exact_negative_group_rejections: u64,
    recovered_anchor_scalars: u64,
    anchor_recovery_failures: u64,
    distinct_column_failures: u64,
    termwise_aggregate_failures: u64,
    common_translate_pair_checks: u64,
    common_translate_pair_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    instance_stream_sha256: String,
    operations: GroupOperations,
}

#[derive(Clone, Serialize)]
struct ExactTrichotomy {
    factor_base_law: String,
    target_law: String,
    delta_a: String,
    delta_b: String,
    identity_case: String,
    impossible_case: String,
    informative_case: String,
    common_translate_pair_law: String,
    balanced_s17_anchor_coefficient: i64,
    informative_relation_is_direct_dlp_witness: bool,
    identity_relation_supplies_target_information: bool,
}

#[derive(Clone, Serialize)]
struct DegreeBoundary {
    comparison_factor_base: ImportedDegreeBoundary,
    affine_log_structured_residual_degree_of_regularity: Option<u64>,
    affine_log_degree_at_most_5_gate_passed: bool,
    common_translate_support_polynomial_degree: u64,
    unsplit_s17_degree_of_regularity: Option<u64>,
    classification: String,
}

#[derive(Clone, Serialize)]
struct ResourceAccounting {
    columns: u64,
    point_key_bytes: u64,
    materialized_point_bytes: u64,
    coefficient_bytes_per_column: u64,
    materialized_point_and_coefficient_bytes: u64,
    storage_below_2_50: bool,
    ordinary_independent_relation_rows: u64,
    projected_relation_collection_operations: Option<String>,
    projected_cost_per_usable_relation: Option<String>,
    rho_runtime_ratio_claimed: Option<f64>,
    rho_classification: String,
}

#[derive(Clone, Serialize)]
struct Gates {
    zero_false_positives_and_false_negatives: bool,
    exact_group_replay: bool,
    structured_residual_degree_at_most_5: bool,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    informative_selector_non_generic: bool,
    no_discarded_branch_counted_as_exhaustive: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    relation_model: String,
    family: String,
    dependencies: Vec<Dependency>,
    exact_trichotomy: ExactTrichotomy,
    toy_controls: Vec<ToyCell>,
    p256_controls: P256Controls,
    degree_boundary: DegreeBoundary,
    resource_accounting: ResourceAccounting,
    semantic_evidence_sha256: String,
    gates: Gates,
    full_depth_unplanted_attempted: bool,
    selector_implemented_or_benchmarked: bool,
    dominant_obstruction: String,
    decision: String,
}

fn json_at<'a>(value: &'a Value, path: &[&str]) -> Result<&'a Value, String> {
    let mut cursor = value;
    for key in path {
        cursor = cursor
            .get(*key)
            .ok_or_else(|| format!("missing JSON path {}", path.join(".")))?;
    }
    Ok(cursor)
}

fn json_str<'a>(value: &'a Value, path: &[&str]) -> Result<&'a str, String> {
    json_at(value, path)?
        .as_str()
        .ok_or_else(|| format!("JSON path {} is not a string", path.join(".")))
}

fn json_u64(value: &Value, path: &[&str]) -> Result<u64, String> {
    json_at(value, path)?
        .as_u64()
        .ok_or_else(|| format!("JSON path {} is not a u64", path.join(".")))
}

fn json_bool(value: &Value, path: &[&str]) -> Result<bool, String> {
    json_at(value, path)?
        .as_bool()
        .ok_or_else(|| format!("JSON path {} is not a bool", path.join(".")))
}

fn load_dependency(
    path: &Path,
    round: u64,
    expected_hash: &str,
    expected_schema: &str,
) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != expected_hash {
        return Err(format!(
            "round-{round} dependency hash mismatch: expected {expected_hash}, got {digest}"
        ));
    }
    let value: Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("invalid JSON: {error}"))?;
    if json_str(&value, &["schema"])? != expected_schema
        || json_str(&value, &["curve"])? != CURVE_SLUG
    {
        return Err(format!("round-{round} dependency identity mismatch"));
    }
    Ok((
        Dependency {
            round,
            path: path.display().to_string(),
            sha256: digest,
            schema: expected_schema.into(),
        },
        value,
    ))
}

fn add_mod(left: u64, right: u64, modulus: u64) -> u64 {
    (left + right) % modulus
}

fn sub_mod(left: u64, right: u64, modulus: u64) -> u64 {
    (left + modulus - right) % modulus
}

fn mul_mod(left: u64, right: u64, modulus: u64) -> u64 {
    ((u128::from(left) * u128::from(right)) % u128::from(modulus)) as u64
}

fn pow_mod(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1u64;
    while exponent != 0 {
        if exponent & 1 == 1 {
            result = mul_mod(result, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exponent >>= 1;
    }
    result
}

fn toy_coefficients(modulus: u64, common_translate: bool) -> (Vec<u64>, Vec<u64>) {
    let mut a = Vec::with_capacity(5);
    let mut b = Vec::with_capacity(5);
    for index in 0..5u64 {
        a.push(if common_translate {
            1
        } else {
            (index * index + 3 * index + 1) % modulus
        });
        b.push((index * index * index + 5 * index + 2) % modulus);
    }
    (a, b)
}

fn toy_cell(modulus: u64, common_translate: bool) -> ToyCell {
    let (a, b) = toy_coefficients(modulus, common_translate);
    let masks = (0u64..(1u64 << 5))
        .filter(|mask| mask.count_ones() == 2)
        .collect::<Vec<_>>();
    let mut row = ToyCell {
        modulus,
        family: if common_translate {
            "common-translate a_i=1".into()
        } else {
            "general affine-log".into()
        },
        width: 5,
        signed_terms: 5,
        positive_terms: 2,
        coefficient_patterns: masks.len() as u64,
        candidates_checked: 0,
        identity_true: 0,
        impossible: 0,
        informative_true: 0,
        informative_false: 0,
        modular_inversions: 0,
        anchor_recovery_failures: 0,
        pair_difference_checks: 0,
        pair_difference_failures: 0,
        false_positives: 0,
        false_negatives: 0,
    };
    for h in 0..modulus {
        for u in 0..modulus {
            for v in 0..modulus {
                for &mask in &masks {
                    let mut sum_a = 0u64;
                    let mut sum_b = 0u64;
                    for index in 0..5 {
                        if mask & (1u64 << index) != 0 {
                            sum_a = add_mod(sum_a, a[index], modulus);
                            sum_b = add_mod(sum_b, b[index], modulus);
                        } else {
                            sum_a = sub_mod(sum_a, a[index], modulus);
                            sum_b = sub_mod(sum_b, b[index], modulus);
                        }
                    }
                    let target = add_mod(u, mul_mod(v, h, modulus), modulus);
                    let rhs = add_mod(sum_b, mul_mod(sum_a, h, modulus), modulus);
                    let direct = target == rhs;
                    let delta_a = sub_mod(v, sum_a, modulus);
                    let delta_b = sub_mod(sum_b, u, modulus);
                    let classified = if delta_a == 0 {
                        if delta_b == 0 {
                            row.identity_true += 1;
                            true
                        } else {
                            row.impossible += 1;
                            false
                        }
                    } else if mul_mod(delta_a, h, modulus) == delta_b {
                        row.informative_true += 1;
                        row.modular_inversions += 1;
                        let recovered =
                            mul_mod(delta_b, pow_mod(delta_a, modulus - 2, modulus), modulus);
                        if recovered != h {
                            row.anchor_recovery_failures += 1;
                        }
                        true
                    } else {
                        row.informative_false += 1;
                        false
                    };
                    row.candidates_checked += 1;
                    if classified && !direct {
                        row.false_positives += 1;
                    }
                    if !classified && direct {
                        row.false_negatives += 1;
                    }
                }
            }
        }
    }
    if common_translate {
        for h in 0..modulus {
            for left in 0..5 {
                for right in 0..5 {
                    if left == right {
                        continue;
                    }
                    let point_left = add_mod(h, b[left], modulus);
                    let point_right = add_mod(h, b[right], modulus);
                    let observed = sub_mod(point_left, point_right, modulus);
                    let expected = sub_mod(b[left], b[right], modulus);
                    row.pair_difference_checks += 1;
                    if observed != expected {
                        row.pair_difference_failures += 1;
                    }
                }
            }
        }
    }
    row
}

fn big_add_mod(left: &BigUint, right: &BigUint, modulus: &BigUint) -> BigUint {
    (left + right) % modulus
}

fn big_sub_mod(left: &BigUint, right: &BigUint, modulus: &BigUint) -> BigUint {
    if left >= right {
        (left - right) % modulus
    } else {
        (modulus - ((right - left) % modulus)) % modulus
    }
}

fn signed_accumulate(
    accumulator: &BigUint,
    value: &BigUint,
    positive: bool,
    modulus: &BigUint,
) -> BigUint {
    if positive {
        big_add_mod(accumulator, value, modulus)
    } else {
        big_sub_mod(accumulator, value, modulus)
    }
}

fn hash_scalar(domain: &str, modulus: &BigUint, nonzero: bool) -> BigUint {
    let mut value = BigUint::from_bytes_be(&sha256(domain.as_bytes())) % modulus;
    if nonzero && value.is_zero() {
        value = BigUint::one();
    }
    value
}

fn coefficient(index: u64, kind: char, modulus: &BigUint) -> BigUint {
    hash_scalar(
        &format!("{CURVE_SLUG}/affine-log-round39/{kind}/{index}"),
        modulus,
        kind == 'a',
    )
}

fn selected_indices(instance: u64) -> Vec<u64> {
    let mut seen = BTreeSet::new();
    let mut counter = 0u64;
    while seen.len() < RELATION_TERMS {
        let digest = sha256(
            format!("{CURVE_SLUG}/affine-log-round39/indices/{instance}/{counter}").as_bytes(),
        );
        let mut word = [0u8; 8];
        word.copy_from_slice(&digest[..8]);
        seen.insert(u64::from_be_bytes(word) % COLUMNS);
        counter += 1;
    }
    seen.into_iter().collect()
}

fn scalar_mul(
    point: &Projective,
    scalar: &BigUint,
    operations: &mut GroupOperations,
) -> Projective {
    operations.scalar_multiplications += 1;
    let mut result = Projective::IDENTITY;
    for bit in (0..scalar.bits()).rev() {
        result = result.double();
        operations.doublings += 1;
        if scalar.bit(bit) {
            result = result.add(point);
            operations.additions += 1;
        }
    }
    result
}

fn linear_point(
    generator: &Projective,
    anchor: &Projective,
    generator_scalar: &BigUint,
    anchor_scalar: &BigUint,
    operations: &mut GroupOperations,
) -> Projective {
    let left = scalar_mul(generator, generator_scalar, operations);
    let right = scalar_mul(anchor, anchor_scalar, operations);
    operations.additions += 1;
    left.add(&right)
}

fn points_equal(left: &Projective, right: &Projective, operations: &mut GroupOperations) -> bool {
    operations.additions += 1;
    bool::from(left.add(&right.neg()).is_identity())
}

fn append_scalar(stream: &mut Vec<u8>, scalar: &BigUint) {
    let bytes = scalar.to_bytes_be();
    stream.extend(std::iter::repeat_n(0, 32usize.saturating_sub(bytes.len())));
    stream.extend(bytes);
}

fn p256_controls(curve: &CurveParams) -> P256Controls {
    let modulus = &curve.n;
    let generator = Projective::from_affine(
        &U256::from_biguint(&curve.gx),
        &U256::from_biguint(&curve.gy),
    );
    let anchor_scalar = hash_scalar(
        &format!("{CURVE_SLUG}/affine-log-round39/control-anchor"),
        modulus,
        true,
    );
    let mut operations = GroupOperations::default();
    let anchor = scalar_mul(&generator, &anchor_scalar, &mut operations);
    let mut identity_relations = 0u64;
    let mut informative_relations = 0u64;
    let mut exact_positive_group_replays = 0u64;
    let mut exact_negative_group_rejections = 0u64;
    let mut recovered_anchor_scalars = 0u64;
    let mut anchor_recovery_failures = 0u64;
    let mut distinct_column_failures = 0u64;
    let mut termwise_aggregate_failures = 0u64;
    let mut false_positives = 0u64;
    let mut false_negatives = 0u64;
    let mut stream = Vec::new();

    for instance in 0..P256_INSTANCES {
        let indices = selected_indices(instance);
        if indices.len() != RELATION_TERMS {
            distinct_column_failures += 1;
        }
        let mut sum_a = BigUint::zero();
        let mut sum_b = BigUint::zero();
        let mut termwise_rhs = Projective::IDENTITY;
        stream.extend(instance.to_be_bytes());
        for (position, index) in indices.iter().enumerate() {
            let a = coefficient(*index, 'a', modulus);
            let b = coefficient(*index, 'b', modulus);
            let positive = position < POSITIVE_TERMS;
            sum_a = signed_accumulate(&sum_a, &a, positive, modulus);
            sum_b = signed_accumulate(&sum_b, &b, positive, modulus);
            let point = linear_point(&generator, &anchor, &b, &a, &mut operations);
            operations.additions += 1;
            termwise_rhs = if positive {
                termwise_rhs.add(&point)
            } else {
                termwise_rhs.add(&point.neg())
            };
            stream.extend(index.to_be_bytes());
            append_scalar(&mut stream, &a);
            append_scalar(&mut stream, &b);
        }

        let (u, v, informative) = if instance % 2 == 0 {
            identity_relations += 1;
            (sum_b.clone(), sum_a.clone(), false)
        } else {
            informative_relations += 1;
            let delta_a = hash_scalar(
                &format!("{CURVE_SLUG}/affine-log-round39/delta/{instance}"),
                modulus,
                true,
            );
            let delta_h = &delta_a * &anchor_scalar % modulus;
            (
                big_sub_mod(&sum_b, &delta_h, modulus),
                big_add_mod(&sum_a, &delta_a, modulus),
                true,
            )
        };
        append_scalar(&mut stream, &u);
        append_scalar(&mut stream, &v);

        let target = linear_point(&generator, &anchor, &u, &v, &mut operations);
        let aggregate_rhs = linear_point(&generator, &anchor, &sum_b, &sum_a, &mut operations);
        if !points_equal(&termwise_rhs, &aggregate_rhs, &mut operations) {
            termwise_aggregate_failures += 1;
        }
        let positive_replay = points_equal(&target, &termwise_rhs, &mut operations);
        if positive_replay {
            exact_positive_group_replays += 1;
        } else {
            false_negatives += 1;
        }

        let negative_u = big_add_mod(&u, &BigUint::one(), modulus);
        let negative_target = linear_point(&generator, &anchor, &negative_u, &v, &mut operations);
        let negative_replay = points_equal(&negative_target, &termwise_rhs, &mut operations);
        if negative_replay {
            false_positives += 1;
        } else {
            exact_negative_group_rejections += 1;
        }

        let delta_a = big_sub_mod(&v, &sum_a, modulus);
        let delta_b = big_sub_mod(&sum_b, &u, modulus);
        if informative {
            let inverse = delta_a.modpow(&(modulus - 2u8), modulus);
            let recovered = &delta_b * inverse % modulus;
            if recovered == anchor_scalar {
                recovered_anchor_scalars += 1;
            } else {
                anchor_recovery_failures += 1;
            }
        } else if !delta_a.is_zero() || !delta_b.is_zero() {
            false_positives += 1;
        }
    }

    let mut common_translate_pair_failures = 0u64;
    for control in 0..P256_PAIR_CONTROLS {
        let left_index = control * 2 % COLUMNS;
        let mut right_index = (control * 2 + 1) % COLUMNS;
        if right_index == left_index {
            right_index = (right_index + 1) % COLUMNS;
        }
        let left_b = coefficient(left_index, 'b', modulus);
        let right_b = coefficient(right_index, 'b', modulus);
        let left = linear_point(
            &generator,
            &anchor,
            &left_b,
            &BigUint::one(),
            &mut operations,
        );
        let right = linear_point(
            &generator,
            &anchor,
            &right_b,
            &BigUint::one(),
            &mut operations,
        );
        operations.additions += 1;
        let observed = left.add(&right.neg());
        let difference = big_sub_mod(&left_b, &right_b, modulus);
        let expected = scalar_mul(&generator, &difference, &mut operations);
        if !points_equal(&observed, &expected, &mut operations) {
            common_translate_pair_failures += 1;
        }
        stream.extend(control.to_be_bytes());
        append_scalar(&mut stream, &left_b);
        append_scalar(&mut stream, &right_b);
    }

    P256Controls {
        deterministic_instances: P256_INSTANCES,
        identity_relations,
        informative_relations,
        impossible_negative_controls: P256_INSTANCES,
        exact_positive_group_replays,
        exact_negative_group_rejections,
        recovered_anchor_scalars,
        anchor_recovery_failures,
        distinct_column_failures,
        termwise_aggregate_failures,
        common_translate_pair_checks: P256_PAIR_CONTROLS,
        common_translate_pair_failures,
        false_positives,
        false_negatives,
        instance_stream_sha256: hex::encode(sha256(&stream)),
        operations,
    }
}

fn run(cli: &Cli) -> Result<ResultFile, String> {
    let (round19, round19_value) = load_dependency(
        &cli.round19,
        19,
        ROUND19_SHA256,
        "p256.s17_multilevel_selector/v1",
    )?;
    let (round25, round25_value) = load_dependency(
        &cli.round25,
        25,
        ROUND25_SHA256,
        "p256.scalar_orbit_factor_base_screen/v1",
    )?;
    if json_str(&round19_value, &["factor_base", "fb_id"])? != COMPARISON_FB
        || json_u64(&round19_value, &["factor_base", "columns"])? != COLUMNS
        || json_u64(
            &round19_value,
            &["residual_degree_evidence", "structured_residual_maximum"],
        )? != 4
        || json_str(&round25_value, &["factor_base", "fb_id"])? != "FB1he6b6b6e25de6"
        || json_u64(
            &round25_value,
            &["factor_base", "independent_log_quotient_dimension"],
        )? != 1
    {
        return Err("frozen dependency facts changed".into());
    }

    let mut toy_controls = Vec::new();
    for modulus in [19, 23, 29, 31] {
        toy_controls.push(toy_cell(modulus, false));
        toy_controls.push(toy_cell(modulus, true));
    }
    let p256 = p256_controls(&CurveParams::p256());
    let imported = ImportedDegreeBoundary {
        factor_base: COMPARISON_FB.into(),
        structured_residual_maximum: 4,
        every_component_complete_and_correct: json_bool(
            &round19_value,
            &[
                "residual_degree_evidence",
                "every_component_complete_and_correct",
            ],
        )?,
        applies_to_affine_log_family: false,
        unsplit_s17_degree_of_regularity: None,
    };
    let exact_trichotomy = ExactTrichotomy {
        factor_base_law: "P_i=[a_i]H+[b_i]G".into(),
        target_law: "R=[u]G+[v]H".into(),
        delta_a: "v-sum(epsilon_i*a_i) mod n".into(),
        delta_b: "sum(epsilon_i*b_i)-u mod n".into(),
        identity_case: "Delta_a=0 and Delta_b=0: coefficient identity, no anchor information"
            .into(),
        impossible_case: "Delta_a=0 and Delta_b!=0: impossible in the prime-order subgroup".into(),
        informative_case: "Delta_a!=0: a true equality returns h=Delta_b/Delta_a".into(),
        common_translate_pair_law: "(H+[b_i]G)-(H+[b_j]G)=[b_i-b_j]G".into(),
        balanced_s17_anchor_coefficient: -1,
        informative_relation_is_direct_dlp_witness: true,
        identity_relation_supplies_target_information: false,
    };
    let degree_boundary = DegreeBoundary {
        comparison_factor_base: imported,
        affine_log_structured_residual_degree_of_regularity: None,
        affine_log_degree_at_most_5_gate_passed: false,
        common_translate_support_polynomial_degree: COLUMNS,
        unsplit_s17_degree_of_regularity: None,
        classification: "The completed degree-4 structured residual belongs only to FB1h2f8621cda105. Affine-log and common-translate bases have no completed residual-degree measurement; a B-point coordinate support polynomial has degree B."
            .into(),
    };
    let materialized_point_bytes = COLUMNS * 33;
    let materialized_point_and_coefficient_bytes = materialized_point_bytes + COLUMNS * 64;
    let resources = ResourceAccounting {
        columns: COLUMNS,
        point_key_bytes: 33,
        materialized_point_bytes,
        coefficient_bytes_per_column: 64,
        materialized_point_and_coefficient_bytes,
        storage_below_2_50: materialized_point_and_coefficient_bytes < (1u64 << 50),
        ordinary_independent_relation_rows: 0,
        projected_relation_collection_operations: None,
        projected_cost_per_usable_relation: None,
        rho_runtime_ratio_claimed: None,
        rho_classification: "Every informative equality is the original one-dimensional DLP equation. This is an exact reduction, not a measured runtime ratio or a proved generic-group lower bound."
            .into(),
    };
    let toy_failures = toy_controls
        .iter()
        .map(|row| {
            row.false_positives
                + row.false_negatives
                + row.anchor_recovery_failures
                + row.pair_difference_failures
        })
        .sum::<u64>();
    let p256_failures = p256.false_positives
        + p256.false_negatives
        + p256.anchor_recovery_failures
        + p256.distinct_column_failures
        + p256.termwise_aggregate_failures
        + p256.common_translate_pair_failures;
    let exact = toy_failures == 0 && p256_failures == 0;
    let semantic = json!({
        "curve": CURVE_SLUG,
        "round19_sha256": ROUND19_SHA256,
        "round25_sha256": ROUND25_SHA256,
        "trichotomy": &exact_trichotomy,
        "toy_controls": &toy_controls,
        "p256_controls": &p256,
        "degree_boundary": &degree_boundary,
        "resource_accounting": &resources,
    });
    let semantic_evidence_sha256 = hex::encode(sha256(
        &serde_json::to_vec(&semantic).map_err(|error| error.to_string())?,
    ));
    let gates = Gates {
        zero_false_positives_and_false_negatives: exact,
        exact_group_replay: exact,
        structured_residual_degree_at_most_5: false,
        relation_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        materialized_storage_below_2_50: resources.storage_below_2_50,
        informative_selector_non_generic: false,
        no_discarded_branch_counted_as_exhaustive: true,
        promoted: false,
    };
    Ok(ResultFile {
        schema: "p256.affine_log_orbit_screen/v1".into(),
        curve: CURVE_SLUG.into(),
        relation_model: "17 distinct variable columns; signed balanced 8+9".into(),
        family: "affine log orbits P_i=[a_i]H+[b_i]G, including common translates"
            .into(),
        dependencies: vec![round19, round25],
        exact_trichotomy,
        toy_controls,
        p256_controls: p256,
        degree_boundary,
        resource_accounting: resources,
        semantic_evidence_sha256,
        gates,
        full_depth_unplanted_attempted: false,
        selector_implemented_or_benchmarked: false,
        dominant_obstruction: "Known affine log transport collapses every true relation into either a target-independent coefficient identity or the original DLP equation h=Delta_b/Delta_a. Common translations make pair differences free only in the identity branch."
            .into(),
        decision: "Reject affine-log-orbit and common-translate factor bases as an index-calculus route to rho parity. Cheap relations contain no target information; informative relations are direct DLP witnesses, and the family has no measured degree-at-most-five residual system."
            .into(),
    })
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(result) => {
            let mut bytes = match serde_json::to_vec_pretty(&result) {
                Ok(bytes) => bytes,
                Err(error) => {
                    eprintln!("serialization failed: {error}");
                    return ExitCode::FAILURE;
                }
            };
            bytes.push(b'\n');
            if let Some(path) = cli.out {
                if let Err(error) = fs::write(&path, &bytes) {
                    eprintln!("{}: {error}", path.display());
                    return ExitCode::FAILURE;
                }
            } else {
                print!("{}", String::from_utf8_lossy(&bytes));
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("P-256 affine-log-orbit screen failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn toy_trichotomy_is_exact() {
        for common_translate in [false, true] {
            let row = toy_cell(19, common_translate);
            assert_eq!(row.false_positives, 0);
            assert_eq!(row.false_negatives, 0);
            assert_eq!(row.anchor_recovery_failures, 0);
            assert_eq!(row.pair_difference_failures, 0);
            assert!(row.identity_true > 0);
            assert!(row.informative_true > 0);
            assert!(row.impossible > 0);
        }
    }

    #[test]
    fn selected_columns_are_distinct_and_deterministic() {
        let first = selected_indices(7);
        assert_eq!(first.len(), RELATION_TERMS);
        assert_eq!(first, selected_indices(7));
        assert!(first.windows(2).all(|pair| pair[0] < pair[1]));
    }

    #[test]
    fn balanced_s17_common_translate_coefficient_is_minus_one() {
        assert_eq!(
            POSITIVE_TERMS as i64 - (RELATION_TERMS - POSITIVE_TERMS) as i64,
            -1
        );
    }
}
