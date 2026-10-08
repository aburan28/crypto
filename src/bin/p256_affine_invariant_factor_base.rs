//! Complete affine-invariant P-256 factor-base screen (round 27).

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ec_index_calculus::sqrt_mod_p;
use crypto_lib::cryptanalysis::ecbench::canonical::short_id;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{
    WideCurveFacts, WideFactorBaseDump, WideFactorBaseFacts, WideFbPoint, CURVE_EC1, CURVE_ICV1,
    CURVE_SLUG, CURVE_UID, DUMP_SCHEMA, FACTOR_BASE_SCHEMA, POINT_KEY_BYTES, POINT_KEY_ENCODING,
};
use crypto_lib::ct_bignum::U256;
use crypto_lib::ecc::p256_field::P256FieldElement as Fe;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;
use serde_json::{json, Value};

const ROUND25_SHA256: &str = "dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf";
const ROUND26_SHA256: &str = "463c100dd951e2a7a7a419118f84bde90fc419358d8279a430a5a0e6191db029";
const TARGET_COLUMNS: u64 = 131_458;
const MAX_COLUMNS: u64 = 139_592;
const USEFUL_ROWS: u64 = 138_031;
const CONTROLS: u64 = 4_096;
const RHO_S: f64 = 1.3;
const MAX_PARITY_OPS: f64 = 1.036_982_447_222_106_5;

#[derive(Parser)]
#[command(
    about = "Screen every proper affine-invariant P-256 factor base in the frozen width band"
)]
struct Cli {
    #[arg(long)]
    round25: PathBuf,
    #[arg(long)]
    round26: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
    /// Optional native `ecbench.factor_base_dump/v1-wide` materialisation.
    #[arg(long)]
    dump: Option<PathBuf>,
}

#[derive(Clone)]
struct FactorPower {
    prime: BigUint,
    exponent: u32,
}

#[derive(Serialize)]
struct Dependency {
    round: u64,
    path: String,
    sha256: String,
    schema: String,
}

#[derive(Serialize)]
struct ClassificationProof {
    group_order: String,
    group_order_is_odd_prime: bool,
    signed_set_condition: String,
    conjugated_negation: String,
    generated_translation: String,
    nonzero_translation_orbit_length: String,
    proper_affine_invariant_base_requires_b_zero: bool,
    remaining_family: String,
}

#[derive(Default, Serialize)]
struct GroupOperations {
    orbit_generation_additions: u64,
    orbit_generation_doublings: u64,
    transport_replay_additions: u64,
    transport_replay_doublings: u64,
    control_additions: u64,
    control_doublings: u64,
    progression_generation_additions: u64,
    progression_replay_additions: u64,
    normalization_field_multiplications: u64,
    normalization_inversions: u64,
}

#[derive(Serialize)]
struct ProgressionControl {
    columns: u64,
    coefficient_u: String,
    coefficient_v: String,
    internal_edges_replayed: u64,
    internal_edge_failures: u64,
    wrap_replayed: bool,
    wraps_to_first_column: bool,
    wraps_to_negated_first_column: bool,
    wrap_defect_coefficient: String,
    wrap_defect_nonzero: bool,
    cyclic_constant_translation_possible: bool,
    known_logs: bool,
    classification: String,
}

#[derive(Clone)]
struct Action {
    scalar: BigUint,
    source_order: u64,
    width: u64,
    binary_operations: u64,
    doubling_lower_bound: u64,
}

#[derive(Serialize)]
struct OrderMinimum {
    order: u64,
    folded_width: u64,
    exact_order_elements: u64,
    unique_folded_actions: u64,
    minimum_scalar: String,
    minimum_binary_operations: u64,
    minimum_doubling_lower_bound: u64,
}

#[derive(Clone, Serialize)]
struct FrontierRow {
    scalar: String,
    source_order: u64,
    orbit_width: u64,
    orbit_count: u64,
    columns: u64,
    available_orbits: String,
    binary_operations: u64,
    doubling_lower_bound: u64,
    target_corrected_binary_operations: u64,
    target_corrected_doubling_lower_bound: u64,
}

#[derive(Serialize)]
struct ActionCensus {
    n_minus_one_factorization: Vec<String>,
    factorization_product_verified: bool,
    least_primitive_root: String,
    primitive_root_replayed: bool,
    divisors_enumerated: u64,
    retained_nontrivial_orders: u64,
    exact_order_elements_visited: u64,
    unique_folded_actions: u64,
    duplicate_folded_actions_removed: u64,
    width_band_combinations_evaluated: u64,
    ordered_census_sha256: String,
    order_minima: Vec<OrderMinimum>,
    pareto_frontier: Vec<FrontierRow>,
    minimum_cost_candidate: FrontierRow,
    minimum_doubling_candidate: FrontierRow,
    minimum_anchor_candidate: FrontierRow,
    enumeration_exhaustive: bool,
}

#[derive(Clone)]
struct OrbitPoint {
    orbit: u64,
    position: u64,
    point: Projective,
    x: BigUint,
    low_y: BigUint,
}

#[derive(Serialize)]
struct NativeFactorBase {
    fb_id: String,
    fb_sha256: String,
    family: String,
    columns: u64,
    signed_points: u64,
    scalar: String,
    source_order: u64,
    orbit_width: u64,
    independent_anchor_classes: u64,
    anchor_candidates: u64,
    rejected_intersecting_anchors: u64,
    points_sha256: String,
    point_key_encoding: String,
    all_points_nonidentity: bool,
    all_points_on_curve: bool,
    unique_up_to_sign: bool,
    transport_edges_replayed: u64,
    transport_replay_failures: u64,
    controls: u64,
    control_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    dump_schema: String,
    dump_bytes: u64,
    dump_sha256: String,
    rebuild_identity_verified: bool,
    operations: GroupOperations,
}

#[derive(Serialize)]
struct EndToEndProjection {
    rho_s: f64,
    maximum_operations_per_sample_for_parity: f64,
    selected_binary_operations_per_sample: u64,
    selected_doubling_lower_bound_per_sample: u64,
    selected_target_corrected_binary_operations_per_sample: u64,
    selected_target_corrected_doubling_lower_bound_per_sample: u64,
    target_corrected_binary_ratio_to_rho: f64,
    target_corrected_doubling_lower_bound_ratio_to_rho: f64,
    global_minimum_target_corrected_doubling_lower_bound: u64,
    global_minimum_doubling_lower_bound_ratio_to_rho: f64,
    independent_unknown_anchors: u64,
    exact_support_polynomial_degree: u64,
    structured_residual_degree_of_regularity: Option<u32>,
    known_anchor_classification: String,
    unknown_anchor_classification: String,
}

#[derive(Serialize)]
struct PromotionGates {
    zero_false_positives: bool,
    zero_false_negatives: bool,
    exact_transport_replay: bool,
    complete_coalescing_action: bool,
    operations_per_sample_at_most_1_036983: bool,
    structured_residual_degree_at_most_5: bool,
    relation_collection_below_2_120: bool,
    per_usable_row_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    complete_cost_below_rho: bool,
    selector_and_anchor_non_generic: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    relation_model: String,
    target_columns: u64,
    maximum_columns: u64,
    useful_rows: u64,
    dependencies: Vec<Dependency>,
    classification: ClassificationProof,
    arithmetic_progression: ProgressionControl,
    action_census: ActionCensus,
    factor_base: NativeFactorBase,
    projection: EndToEndProjection,
    promotion_gates: PromotionGates,
    full_depth_unplanted_attempted: bool,
    dominant_obstruction: String,
    decision: String,
}

fn parse(value: &str) -> BigUint {
    BigUint::parse_bytes(value.as_bytes(), 10).expect("frozen decimal integer")
}

fn factors() -> Vec<FactorPower> {
    [
        ("2", 4),
        ("3", 1),
        ("71", 1),
        ("131", 1),
        ("373", 1),
        ("3407", 1),
        ("17449", 1),
        ("38189", 1),
        ("187019741", 1),
        ("622491383", 1),
        ("1002328039319", 1),
        ("2624747550333869278416773953", 1),
    ]
    .into_iter()
    .map(|(prime, exponent)| FactorPower {
        prime: parse(prime),
        exponent,
    })
    .collect()
}

fn factor_product(entries: &[FactorPower]) -> BigUint {
    entries.iter().fold(BigUint::one(), |product, entry| {
        product * entry.prime.pow(entry.exponent)
    })
}

fn divisors(entries: &[FactorPower]) -> Vec<BigUint> {
    let mut out = vec![BigUint::one()];
    for entry in entries {
        let old = out.clone();
        let mut power = BigUint::one();
        for _ in 0..entry.exponent {
            power *= &entry.prime;
            out.extend(old.iter().map(|value| value * &power));
        }
    }
    out
}

fn folded_width(order: u64) -> u64 {
    if order.is_even() {
        order / 2
    } else {
        order
    }
}

fn bit_weight(value: &BigUint) -> u64 {
    value
        .to_bytes_le()
        .iter()
        .map(|byte| u64::from(byte.count_ones()))
        .sum()
}

fn action_for(scalar: BigUint, source_order: u64, width: u64) -> Action {
    let bits = scalar.bits();
    let weight = bit_weight(&scalar);
    Action {
        scalar,
        source_order,
        width,
        binary_operations: bits + weight - 2,
        doubling_lower_bound: bits - 1,
    }
}

fn action_is_better(left: &Action, right: &Action) -> bool {
    left.binary_operations < right.binary_operations
        || (left.binary_operations == right.binary_operations
            && (left.doubling_lower_bound < right.doubling_lower_bound
                || (left.doubling_lower_bound == right.doubling_lower_bound
                    && left.scalar < right.scalar)))
}

fn fixed_32(value: &BigUint) -> Result<[u8; 32], String> {
    let raw = value.to_bytes_be();
    if raw.len() > 32 {
        return Err("coefficient exceeded 32 bytes".into());
    }
    let mut out = [0u8; 32];
    out[32 - raw.len()..].copy_from_slice(&raw);
    Ok(out)
}

fn dependency(path: &Path, round: u64, expected: &str) -> Result<Dependency, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != expected {
        return Err(format!(
            "round-{round} dependency hash mismatch: expected {expected}, got {digest}"
        ));
    }
    let value: Value = serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if value.get("curve").and_then(Value::as_str) != Some(CURVE_SLUG) {
        return Err(format!("round-{round} dependency curve mismatch"));
    }
    Ok(Dependency {
        round,
        path: path.display().to_string(),
        sha256: digest,
        schema: value
            .get("schema")
            .and_then(Value::as_str)
            .unwrap_or("unknown")
            .into(),
    })
}

fn primitive_root_replayed(root: &BigUint, modulus: &BigUint, entries: &[FactorPower]) -> bool {
    root.modpow(&(modulus - 1u8), modulus).is_one()
        && entries
            .iter()
            .all(|entry| root.modpow(&((modulus - 1u8) / &entry.prime), modulus) != BigUint::one())
}

fn row(action: &Action, orbit_count: u64, n: &BigUint) -> FrontierRow {
    FrontierRow {
        scalar: action.scalar.to_string(),
        source_order: action.source_order,
        orbit_width: action.width,
        orbit_count,
        columns: action.width * orbit_count,
        available_orbits: ((n - 1u8) / BigUint::from(2 * action.width)).to_string(),
        binary_operations: action.binary_operations,
        doubling_lower_bound: action.doubling_lower_bound,
        target_corrected_binary_operations: action.binary_operations + 1,
        target_corrected_doubling_lower_bound: action.doubling_lower_bound + 1,
    }
}

fn dominates(left: &FrontierRow, right: &FrontierRow) -> bool {
    let no_worse = left.columns <= right.columns
        && left.orbit_count <= right.orbit_count
        && left.binary_operations <= right.binary_operations
        && left.doubling_lower_bound <= right.doubling_lower_bound;
    let strict = left.columns < right.columns
        || left.orbit_count < right.orbit_count
        || left.binary_operations < right.binary_operations
        || left.doubling_lower_bound < right.doubling_lower_bound;
    no_worse && strict
}

fn action_census(curve: &CurveParams) -> Result<(ActionCensus, Action), String> {
    let factorization = factors();
    if factor_product(&factorization) != &curve.n - 1u8 {
        return Err("n-1 factorization did not replay".into());
    }
    let primitive_root = BigUint::from(7u8);
    if !primitive_root_replayed(&primitive_root, &curve.n, &factorization) {
        return Err("primitive-root replay failed".into());
    }
    let all_divisors = divisors(&factorization);
    if all_divisors.len() != 10_240 {
        return Err(format!(
            "expected 10240 divisors, got {}",
            all_divisors.len()
        ));
    }
    let mut orders = all_divisors
        .iter()
        .filter_map(|value| value.to_u64())
        .filter(|order| *order > 2 && folded_width(*order) <= MAX_COLUMNS)
        .collect::<Vec<_>>();
    orders.sort_unstable();

    let mut unique_actions = BTreeMap::<BigUint, Action>::new();
    let mut order_minima = Vec::with_capacity(orders.len());
    let mut exact_order_elements_visited = 0u64;
    let mut duplicate_folded_actions_removed = 0u64;
    for order in &orders {
        let width = folded_width(*order);
        let root = primitive_root.modpow(&((&curve.n - 1u8) / *order), &curve.n);
        let mut current = BigUint::one();
        let mut exact = 0u64;
        let mut local_folded_actions = BTreeSet::new();
        let mut minimum: Option<Action> = None;
        for exponent in 1..*order {
            current = current * &root % &curve.n;
            if exponent.gcd(order) != 1 {
                continue;
            }
            exact += 1;
            exact_order_elements_visited += 1;
            let negative = &curve.n - &current;
            let signed = current.clone().min(negative);
            local_folded_actions.insert(signed.clone());
            let candidate = action_for(signed.clone(), *order, width);
            if minimum
                .as_ref()
                .map(|best| action_is_better(&candidate, best))
                .unwrap_or(true)
            {
                minimum = Some(candidate.clone());
            }
            match unique_actions.entry(signed) {
                std::collections::btree_map::Entry::Vacant(slot) => {
                    slot.insert(candidate);
                }
                std::collections::btree_map::Entry::Occupied(mut slot) => {
                    duplicate_folded_actions_removed += 1;
                    if slot.get().width != width {
                        return Err("one folded action acquired two effective widths".into());
                    }
                    if action_is_better(&candidate, slot.get()) {
                        slot.insert(candidate);
                    }
                }
            }
        }
        let best = minimum.ok_or_else(|| format!("order {order} had no generators"))?;
        order_minima.push(OrderMinimum {
            order: *order,
            folded_width: width,
            exact_order_elements: exact,
            unique_folded_actions: local_folded_actions.len() as u64,
            minimum_scalar: best.scalar.to_string(),
            minimum_binary_operations: best.binary_operations,
            minimum_doubling_lower_bound: best.doubling_lower_bound,
        });
    }

    let mut census_blob = Vec::with_capacity(unique_actions.len() * 56);
    for action in unique_actions.values() {
        census_blob.extend(fixed_32(&action.scalar)?);
        census_blob.extend(action.source_order.to_be_bytes());
        census_blob.extend(action.width.to_be_bytes());
        census_blob.extend(action.binary_operations.to_be_bytes());
    }

    let mut width_minima = BTreeMap::<u64, Action>::new();
    for action in unique_actions.values() {
        let replace = width_minima
            .get(&action.width)
            .map(|best| action_is_better(action, best))
            .unwrap_or(true);
        if replace {
            width_minima.insert(action.width, action.clone());
        }
    }

    let mut combinations = 0u64;
    let mut candidates = Vec::new();
    for action in width_minima.values() {
        let first = TARGET_COLUMNS.div_ceil(action.width);
        let last = MAX_COLUMNS / action.width;
        if first > last {
            continue;
        }
        let available = (&curve.n - 1u8) / BigUint::from(2 * action.width);
        for orbit_count in first..=last {
            combinations += 1;
            if BigUint::from(orbit_count) > available {
                return Err("requested more scalar orbits than exist".into());
            }
        }
        candidates.push(row(action, first, &curve.n));
    }
    if candidates.is_empty() {
        return Err("no invariant union reached the width band".into());
    }
    let mut frontier = candidates
        .iter()
        .enumerate()
        .filter(|(index, candidate)| {
            !candidates
                .iter()
                .enumerate()
                .any(|(other, row)| other != *index && dominates(row, candidate))
        })
        .map(|(_, candidate)| candidate)
        .collect::<Vec<_>>();
    frontier.sort_by(|left, right| {
        left.columns
            .cmp(&right.columns)
            .then_with(|| left.orbit_count.cmp(&right.orbit_count))
            .then_with(|| left.binary_operations.cmp(&right.binary_operations))
            .then_with(|| left.scalar.cmp(&right.scalar))
    });
    let minimum_cost = candidates
        .iter()
        .min_by(|left, right| {
            left.binary_operations
                .cmp(&right.binary_operations)
                .then_with(|| left.doubling_lower_bound.cmp(&right.doubling_lower_bound))
                .then_with(|| left.orbit_count.cmp(&right.orbit_count))
                .then_with(|| left.columns.cmp(&right.columns))
                .then_with(|| left.scalar.cmp(&right.scalar))
        })
        .ok_or("missing minimum-cost candidate")?;
    let minimum_anchor = candidates
        .iter()
        .min_by(|left, right| {
            left.orbit_count
                .cmp(&right.orbit_count)
                .then_with(|| left.binary_operations.cmp(&right.binary_operations))
                .then_with(|| left.columns.cmp(&right.columns))
                .then_with(|| left.scalar.cmp(&right.scalar))
        })
        .ok_or("missing minimum-anchor candidate")?;
    let minimum_doubling = candidates
        .iter()
        .min_by(|left, right| {
            left.doubling_lower_bound
                .cmp(&right.doubling_lower_bound)
                .then_with(|| left.binary_operations.cmp(&right.binary_operations))
                .then_with(|| left.orbit_count.cmp(&right.orbit_count))
                .then_with(|| left.columns.cmp(&right.columns))
                .then_with(|| left.scalar.cmp(&right.scalar))
        })
        .ok_or("missing minimum-doubling candidate")?;
    let selected_action = unique_actions
        .get(&parse(&minimum_cost.scalar))
        .ok_or("selected action disappeared from census")?
        .clone();
    let census = ActionCensus {
        n_minus_one_factorization: factorization
            .iter()
            .map(|entry| format!("{}^{}", entry.prime, entry.exponent))
            .collect(),
        factorization_product_verified: true,
        least_primitive_root: primitive_root.to_string(),
        primitive_root_replayed: true,
        divisors_enumerated: all_divisors.len() as u64,
        retained_nontrivial_orders: orders.len() as u64,
        exact_order_elements_visited,
        unique_folded_actions: unique_actions.len() as u64,
        duplicate_folded_actions_removed,
        width_band_combinations_evaluated: combinations,
        ordered_census_sha256: hex::encode(sha256(&census_blob)),
        order_minima,
        pareto_frontier: frontier.into_iter().cloned().collect(),
        minimum_cost_candidate: row(&selected_action, minimum_cost.orbit_count, &curve.n),
        minimum_doubling_candidate: minimum_doubling.clone(),
        minimum_anchor_candidate: minimum_anchor.clone(),
        enumeration_exhaustive: true,
    };
    Ok((census, selected_action))
}

fn binary_scalar_mul(
    point: &Projective,
    scalar: &BigUint,
    additions: &mut u64,
    doublings: &mut u64,
) -> Projective {
    let mut result = Projective::IDENTITY;
    let mut addend = *point;
    for bit in 0..scalar.bits() {
        if scalar.bit(bit) {
            result = result.add(&addend);
            *additions += 1;
        }
        if bit + 1 < scalar.bits() {
            addend = addend.double();
            *doublings += 1;
        }
    }
    result
}

fn normalize_points(
    points: &[Projective],
    operations: &mut GroupOperations,
) -> Result<Vec<(BigUint, BigUint)>, String> {
    if points.iter().any(|point| bool::from(point.is_identity())) {
        return Err("factor-base orbit contains identity".into());
    }
    let mut prefixes = Vec::with_capacity(points.len());
    let mut product = Fe::ONE;
    for point in points {
        prefixes.push(product);
        product = product.mul(&point.z);
        operations.normalization_field_multiplications += 1;
    }
    let mut inverse = product.inv();
    operations.normalization_inversions += 1;
    let mut inverses = vec![Fe::ZERO; points.len()];
    for index in (0..points.len()).rev() {
        inverses[index] = inverse.mul(&prefixes[index]);
        operations.normalization_field_multiplications += 1;
        if index != 0 {
            inverse = inverse.mul(&points[index].z);
            operations.normalization_field_multiplications += 1;
        }
    }
    Ok(points
        .iter()
        .zip(inverses)
        .map(|(point, inverse)| {
            operations.normalization_field_multiplications += 2;
            (
                point.x.mul(&inverse).to_canonical().to_biguint(),
                point.y.mul(&inverse).to_canonical().to_biguint(),
            )
        })
        .collect())
}

fn derive_anchor(curve: &CurveParams, orbit: u64, attempt: u64) -> Result<Projective, String> {
    for counter in 0u64.. {
        let digest = sha256(
            format!("{CURVE_SLUG}/affine-invariant-round27/anchor/{orbit}/{attempt}/{counter}")
                .as_bytes(),
        );
        let x = BigUint::from_bytes_be(&digest) % &curve.p;
        let rhs = ((&x * &x % &curve.p) * &x + &curve.a * &x + &curve.b) % &curve.p;
        let Some(y) = sqrt_mod_p(&rhs, &curve.p) else {
            continue;
        };
        if y.is_zero() {
            continue;
        }
        let low_y = y.clone().min(&curve.p - &y);
        return Ok(Projective::from_affine(
            &U256::from_biguint(&x),
            &U256::from_biguint(&low_y),
        ));
    }
    unreachable!()
}

fn lower_hex(value: &BigUint) -> String {
    format!("0x{}", value.to_str_radix(16))
}

fn point_key_bytes(x: &BigUint, high_sign: bool) -> Result<[u8; POINT_KEY_BYTES], String> {
    let key = ((x + BigUint::one()) << 1usize) + BigUint::from(high_sign);
    let raw = key.to_bytes_be();
    if raw.len() > POINT_KEY_BYTES {
        return Err("P-256 point key exceeded 33 bytes".into());
    }
    let mut out = [0u8; POINT_KEY_BYTES];
    out[POINT_KEY_BYTES - raw.len()..].copy_from_slice(&raw);
    Ok(out)
}

fn curve_facts(curve: &CurveParams) -> WideCurveFacts {
    WideCurveFacts {
        slug: CURVE_SLUG.into(),
        icv1: CURVE_ICV1.into(),
        family: "prime".into(),
        construction: "CurveParams::p256()".into(),
        field_degree: None,
        field_bits: curve.p.bits() as u32,
        group_order: curve.n.to_string(),
        r: curve.n.to_string(),
        cofactor: curve.h.to_string(),
        generator: [lower_hex(&curve.gx), lower_hex(&curve.gy)],
        automorphisms_available: 2,
        registered: true,
        ec1: CURVE_EC1.into(),
        curve_uid: CURVE_UID.into(),
    }
}

fn progression_control(
    curve: &CurveParams,
    operations: &mut GroupOperations,
) -> Result<ProgressionControl, String> {
    let generator = Projective::from_affine(
        &U256::from_biguint(&curve.gx),
        &U256::from_biguint(&curve.gy),
    );
    let mut points = Vec::with_capacity(TARGET_COLUMNS as usize);
    let mut point = generator;
    points.push(point);
    for _ in 1..TARGET_COLUMNS {
        point = point.add(&generator);
        operations.progression_generation_additions += 1;
        points.push(point);
    }
    let mut failures = 0u64;
    for index in 0..points.len() - 1 {
        let image = points[index].add(&generator);
        operations.progression_replay_additions += 1;
        if !bool::from(image.add(&points[index + 1].neg()).is_identity()) {
            failures += 1;
        }
    }
    let wrap_image = points
        .last()
        .ok_or("empty arithmetic progression")?
        .add(&generator);
    operations.progression_replay_additions += 1;
    let wraps_first = bool::from(wrap_image.add(&points[0].neg()).is_identity());
    let wraps_negative = bool::from(wrap_image.add(&points[0]).is_identity());
    let defect = BigUint::from(TARGET_COLUMNS) % &curve.n;
    Ok(ProgressionControl {
        columns: TARGET_COLUMNS,
        coefficient_u: "1".into(),
        coefficient_v: "1".into(),
        internal_edges_replayed: TARGET_COLUMNS - 1,
        internal_edge_failures: failures,
        wrap_replayed: true,
        wraps_to_first_column: wraps_first,
        wraps_to_negated_first_column: wraps_negative,
        wrap_defect_coefficient: defect.to_string(),
        wrap_defect_nonzero: !defect.is_zero(),
        cyclic_constant_translation_possible: defect.is_zero(),
        known_logs: true,
        classification: "All internal +G edges replay, but wrapping an index subtracts [B]G. The correction depends on representation state, so the proper progression is not a cyclic memoryless walk; with exposed coefficients it is also a generic representation encoding.".into(),
    })
}

fn build_factor_base(
    curve: &CurveParams,
    action: &Action,
    orbit_count: u64,
    operations: &mut GroupOperations,
    dump_path: Option<&Path>,
) -> Result<NativeFactorBase, String> {
    let columns = action.width * orbit_count;
    let wrap_coefficient = action.scalar.modpow(&BigUint::from(action.width), &curve.n);
    if wrap_coefficient != BigUint::one() && wrap_coefficient != &curve.n - 1u8 {
        return Err("selected action did not close modulo sign at its declared width".into());
    }
    let mut used_x = BTreeSet::new();
    let mut orbit_points = Vec::with_capacity(columns as usize);
    let mut anchor_candidates = 0u64;
    let mut rejected = 0u64;
    for orbit in 0..orbit_count {
        let mut attempt = 0u64;
        loop {
            anchor_candidates += 1;
            let anchor = derive_anchor(curve, orbit, attempt)?;
            let mut points = Vec::with_capacity(action.width as usize);
            let mut point = anchor;
            for _ in 0..action.width {
                points.push(point);
                point = binary_scalar_mul(
                    &point,
                    &action.scalar,
                    &mut operations.orbit_generation_additions,
                    &mut operations.orbit_generation_doublings,
                );
            }
            let affine = normalize_points(&points, operations)?;
            let xs = affine
                .iter()
                .map(|(x, _)| x.clone())
                .collect::<BTreeSet<_>>();
            if xs.len() != action.width as usize || xs.iter().any(|x| used_x.contains(x)) {
                rejected += 1;
                attempt += 1;
                continue;
            }
            for x in xs {
                used_x.insert(x);
            }
            for (position, (point, (x, y))) in points.into_iter().zip(affine).enumerate() {
                let low_y = y.clone().min(&curve.p - &y);
                orbit_points.push(OrbitPoint {
                    orbit,
                    position: position as u64,
                    point,
                    x,
                    low_y,
                });
            }
            break;
        }
    }
    if orbit_points.len() != columns as usize {
        return Err("native factor base has the wrong number of columns".into());
    }

    let mut replay_failures = 0u64;
    for start in (0..orbit_points.len()).step_by(action.width as usize) {
        for offset in 0..action.width as usize {
            let index = start + offset;
            let image = binary_scalar_mul(
                &orbit_points[index].point,
                &action.scalar,
                &mut operations.transport_replay_additions,
                &mut operations.transport_replay_doublings,
            );
            let expected = if offset + 1 == action.width as usize {
                orbit_points[start].point
            } else {
                orbit_points[index + 1].point
            };
            let equal = bool::from(image.add(&expected.neg()).is_identity());
            let opposite = bool::from(image.add(&expected).is_identity());
            if !equal && !opposite {
                replay_failures += 1;
            }
        }
    }

    let mut control_failures = 0u64;
    for control in 0..CONTROLS {
        let digest =
            sha256(format!("{CURVE_SLUG}/affine-invariant-round27/control/{control}").as_bytes());
        let index = (BigUint::from_bytes_be(&digest) % BigUint::from(columns))
            .to_usize()
            .ok_or("control index did not fit usize")?;
        let position = orbit_points[index].position;
        let start = index - position as usize;
        let expected = if position + 1 == action.width {
            orbit_points[start].point
        } else {
            orbit_points[index + 1].point
        };
        let image = binary_scalar_mul(
            &orbit_points[index].point,
            &action.scalar,
            &mut operations.control_additions,
            &mut operations.control_doublings,
        );
        let equal = bool::from(image.add(&expected.neg()).is_identity());
        let opposite = bool::from(image.add(&expected).is_identity());
        if !equal && !opposite {
            control_failures += 1;
        }
    }

    let all_points_on_curve = orbit_points.iter().all(|row| {
        (&row.low_y * &row.low_y) % &curve.p
            == ((&row.x * &row.x % &curve.p) * &row.x + &curve.a * &row.x + &curve.b) % &curve.p
    });
    orbit_points.sort_by(|left, right| {
        left.x
            .cmp(&right.x)
            .then_with(|| left.orbit.cmp(&right.orbit))
            .then_with(|| left.position.cmp(&right.position))
    });
    let unique_up_to_sign = orbit_points.windows(2).all(|pair| pair[0].x < pair[1].x);
    if replay_failures != 0 || control_failures != 0 || !all_points_on_curve || !unique_up_to_sign {
        return Err("native affine-invariant factor-base replay failed".into());
    }

    let mut point_blob = Vec::with_capacity(orbit_points.len() * 2 * POINT_KEY_BYTES);
    let mut dump_points = Vec::with_capacity(orbit_points.len() * 2);
    for (column, point) in orbit_points.iter().enumerate() {
        point_blob.extend(point_key_bytes(&point.x, false)?);
        point_blob.extend(point_key_bytes(&point.x, true)?);
        dump_points.push(WideFbPoint {
            x: lower_hex(&point.x),
            y: lower_hex(&point.low_y),
            col: column as u64,
            coef: "1".into(),
        });
        dump_points.push(WideFbPoint {
            x: lower_hex(&point.x),
            y: lower_hex(&(&curve.p - &point.low_y)),
            col: column as u64,
            coef: (&curve.n - 1u8).to_string(),
        });
    }
    let points_sha256 = hex::encode(sha256(&point_blob));
    let mut params = BTreeMap::new();
    params.insert(
        "anchor_derivation".into(),
        "sha256-try-and-increment-lower-root-per-orbit".into(),
    );
    params.insert("independent_anchor_classes".into(), orbit_count.to_string());
    params.insert("negation_folded".into(), "true".into());
    params.insert("orbit_width".into(), action.width.to_string());
    params.insert("point_key_encoding".into(), POINT_KEY_ENCODING.into());
    params.insert("scalar".into(), action.scalar.to_string());
    params.insert("source_order".into(), action.source_order.to_string());
    let identity = json!({
        "schema": FACTOR_BASE_SCHEMA,
        "curve": CURVE_SLUG,
        "family": "affine-invariant-scalar-union",
        "params": params,
        "columns": columns,
        "signed_points": 2 * columns,
        "points_sha256": points_sha256,
    });
    let (fb_id, fb_sha256) = short_id("FB1", &identity)?;
    let factor_base = WideFactorBaseFacts {
        fb_id: fb_id.clone(),
        fb_sha256: fb_sha256.clone(),
        family: "affine-invariant-scalar-union".into(),
        params,
        description: format!(
            "{orbit_count} negation-folded scalar orbits of width {} on {CURVE_SLUG}",
            action.width
        ),
        signed_points: 2 * columns,
        abscissae: columns,
        columns,
        dimension: None,
        points_sha256: points_sha256.clone(),
    };
    let dump = WideFactorBaseDump {
        schema: DUMP_SCHEMA.into(),
        curve: curve_facts(curve),
        factor_base,
        build_adds: operations.orbit_generation_additions,
        build_doubles: operations.orbit_generation_doublings,
        points: dump_points,
    };
    let mut dump_bytes = serde_json::to_vec_pretty(&dump).map_err(|error| error.to_string())?;
    dump_bytes.push(b'\n');
    let dump_sha256 = hex::encode(sha256(&dump_bytes));
    let rebuilt: WideFactorBaseDump =
        serde_json::from_slice(&dump_bytes).map_err(|error| error.to_string())?;
    if rebuilt != dump {
        return Err("native wide dump serialization was not stable".into());
    }
    if let Some(path) = dump_path {
        fs::write(path, &dump_bytes).map_err(|error| format!("{}: {error}", path.display()))?;
    }
    let rebuilt_identity = json!({
        "schema": FACTOR_BASE_SCHEMA,
        "curve": dump.curve.slug,
        "family": dump.factor_base.family,
        "params": dump.factor_base.params,
        "columns": dump.factor_base.columns,
        "signed_points": dump.factor_base.signed_points,
        "points_sha256": dump.factor_base.points_sha256,
    });
    let (rebuilt_id, rebuilt_sha) = short_id("FB1", &rebuilt_identity)?;
    let rebuild_identity_verified = rebuilt_id == fb_id && rebuilt_sha == fb_sha256;
    if !rebuild_identity_verified {
        return Err("native FB1 identity did not rebuild".into());
    }
    Ok(NativeFactorBase {
        fb_id,
        fb_sha256,
        family: "affine-invariant-scalar-union".into(),
        columns,
        signed_points: 2 * columns,
        scalar: action.scalar.to_string(),
        source_order: action.source_order,
        orbit_width: action.width,
        independent_anchor_classes: orbit_count,
        anchor_candidates,
        rejected_intersecting_anchors: rejected,
        points_sha256,
        point_key_encoding: POINT_KEY_ENCODING.into(),
        all_points_nonidentity: true,
        all_points_on_curve,
        unique_up_to_sign,
        transport_edges_replayed: columns,
        transport_replay_failures: replay_failures,
        controls: CONTROLS,
        control_failures,
        false_positives: 0,
        false_negatives: 0,
        dump_schema: DUMP_SCHEMA.into(),
        dump_bytes: dump_bytes.len() as u64,
        dump_sha256,
        rebuild_identity_verified,
        operations: std::mem::take(operations),
    })
}

fn run(cli: &Cli) -> Result<ResultFile, String> {
    let dependencies = vec![
        dependency(&cli.round25, 25, ROUND25_SHA256)?,
        dependency(&cli.round26, 26, ROUND26_SHA256)?,
    ];
    let curve = CurveParams::p256();
    let classification = ClassificationProof {
        group_order: curve.n.to_string(),
        group_order_is_odd_prime: true,
        signed_set_condition: "S=-S, 0 not in S, and f(S)=S".into(),
        conjugated_negation: "f o nu o f^-1: x -> -x+2b".into(),
        generated_translation: "(f o nu o f^-1) o nu: x -> x+2b".into(),
        nonzero_translation_orbit_length: curve.n.to_string(),
        proper_affine_invariant_base_requires_b_zero: true,
        remaining_family: "Unions of scalar orbits on F_n^*/{+-1}".into(),
    };
    let (action_census, selected_action) = action_census(&curve)?;
    let selected_orbits = action_census.minimum_cost_candidate.orbit_count;
    let global_doubling_floor = action_census
        .minimum_doubling_candidate
        .target_corrected_doubling_lower_bound;
    let mut operations = GroupOperations::default();
    let arithmetic_progression = progression_control(&curve, &mut operations)?;
    let factor_base = build_factor_base(
        &curve,
        &selected_action,
        selected_orbits,
        &mut operations,
        cli.dump.as_deref(),
    )?;
    let projection = EndToEndProjection {
        rho_s: RHO_S,
        maximum_operations_per_sample_for_parity: MAX_PARITY_OPS,
        selected_binary_operations_per_sample: selected_action.binary_operations,
        selected_doubling_lower_bound_per_sample: selected_action.doubling_lower_bound,
        selected_target_corrected_binary_operations_per_sample: selected_action.binary_operations
            + 1,
        selected_target_corrected_doubling_lower_bound_per_sample: selected_action
            .doubling_lower_bound
            + 1,
        target_corrected_binary_ratio_to_rho: (selected_action.binary_operations + 1) as f64
            / MAX_PARITY_OPS,
        target_corrected_doubling_lower_bound_ratio_to_rho: (selected_action
            .doubling_lower_bound
            + 1) as f64
            / MAX_PARITY_OPS,
        global_minimum_target_corrected_doubling_lower_bound: global_doubling_floor,
        global_minimum_doubling_lower_bound_ratio_to_rho: global_doubling_floor as f64
            / MAX_PARITY_OPS,
        independent_unknown_anchors: selected_orbits,
        exact_support_polynomial_degree: factor_base.columns,
        structured_residual_degree_of_regularity: None,
        known_anchor_classification: "Publishing every anchor logarithm turns the relation search into a generic representation/claw search; it is not an index-calculus speedup.".into(),
        unknown_anchor_classification: "The native union retains one independent logarithm per orbit before relations. Cheap local transport does not remove these global unknowns.".into(),
    };
    let correctness = factor_base.transport_replay_failures == 0
        && factor_base.control_failures == 0
        && arithmetic_progression.internal_edge_failures == 0;
    let operations_gate = (selected_action.binary_operations + 1) as f64 <= MAX_PARITY_OPS;
    let promotion_gates = PromotionGates {
        zero_false_positives: factor_base.false_positives == 0,
        zero_false_negatives: factor_base.false_negatives == 0,
        exact_transport_replay: correctness,
        complete_coalescing_action: true,
        operations_per_sample_at_most_1_036983: operations_gate,
        structured_residual_degree_at_most_5: false,
        relation_collection_below_2_120: false,
        per_usable_row_below_2_103: false,
        materialized_storage_below_2_50: false,
        complete_cost_below_rho: false,
        selector_and_anchor_non_generic: false,
        promoted: false,
    };
    Ok(ResultFile {
        schema: "p256.affine_invariant_factor_base_screen/v1".into(),
        curve: CURVE_SLUG.into(),
        relation_model: "17 distinct variable columns; signed 8+9 cross-colour collision modulo global negation".into(),
        target_columns: TARGET_COLUMNS,
        maximum_columns: MAX_COLUMNS,
        useful_rows: USEFUL_ROWS,
        dependencies,
        classification,
        arithmetic_progression,
        action_census,
        factor_base,
        projection,
        promotion_gates,
        full_depth_unplanted_attempted: false,
        dominant_obstruction: "Negation forces every proper affine-invariant base to have b=0. The remaining scalar-union actions either cost far more than one group operation per sample or expose/retain many anchor logarithms; progression wrapping reintroduces representation-dependent state.".into(),
        decision: "No proper affine-invariant P-256 factor base passes. Additive progressions are not cyclic at a sub-n width, affine translations force the full group, and the complete remaining scalar-action census misses the rho operation gate before summation-polynomial solving, storage, or sparse linear algebra.".into(),
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
            eprintln!("P-256 affine-invariant screen failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factorization_and_divisor_count_replay() {
        let curve = CurveParams::p256();
        let entries = factors();
        assert_eq!(factor_product(&entries), curve.n - 1u8);
        assert_eq!(divisors(&entries).len(), 10_240);
    }

    #[test]
    fn conjugated_negation_generates_translation() {
        let modulus = 101u64;
        let a = 7u64;
        let b = 13u64;
        let inverse = 29u64;
        assert_eq!(a * inverse % modulus, 1);
        for x in 0..modulus {
            let f_inverse = inverse * ((x + modulus - b) % modulus) % modulus;
            let conjugate = (a * ((modulus - f_inverse) % modulus) + b) % modulus;
            let composed_with_negation = {
                let negative_x = (modulus - x) % modulus;
                let inv = inverse * ((negative_x + modulus - b) % modulus) % modulus;
                (a * ((modulus - inv) % modulus) + b) % modulus
            };
            assert_eq!(conjugate, (2 * b + modulus - x) % modulus);
            assert_eq!(composed_with_negation, (x + 2 * b) % modulus);
        }
    }

    #[test]
    fn folded_width_handles_odd_and_even_orders() {
        assert_eq!(folded_width(3), 3);
        assert_eq!(folded_width(6), 3);
        assert_eq!(folded_width(279_184), 139_592);
    }

    #[test]
    fn proper_progression_cannot_wrap_by_constant_translation() {
        let n = CurveParams::p256().n;
        assert_ne!(BigUint::from(TARGET_COLUMNS) % n, BigUint::zero());
    }
}
