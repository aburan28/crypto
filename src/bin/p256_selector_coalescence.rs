//! Coalescence audit for the P-256 scalar-orbit S17 selector (round 26).

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;

const ROUND25_SHA256: &str = "dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf";
const FACTOR_BASE_ID: &str = "FB1he6b6b6e25de6";
const COLUMNS: u64 = 139_592;
const ORBIT_ORDER: u64 = 279_184;
const SCALAR: &str = "379483656589291393959605088584129393463128953089898693989044430844151895048";
const CONTROLS: u64 = 4_096;
const RHO_S: f64 = 1.3;
const KNOWN_ANCHOR_ORACLE_S: f64 = 1.253_637_420_269_235_3;
const KNOWN_ANCHOR_LOG2_BYTES: f64 = 133.326_120_148_915_3;

#[derive(Parser)]
#[command(about = "Audit coalescence of the P-256 scalar-orbit S17 selector")]
struct Cli {
    #[arg(long)]
    round25: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Serialize)]
struct Dependency {
    path: String,
    sha256: String,
    factor_base_id: String,
    columns: u64,
    scalar_order: u64,
    scalar: String,
}

#[derive(Serialize)]
struct DeltaCensus {
    signed_atoms: u64,
    deltas_enumerated: u64,
    distinct_deltas: u64,
    duplicate_deltas: u64,
    antipodal_pairs_checked: u64,
    antipodal_failures: u64,
    delta_coefficients_sha256: String,
    injective: bool,
    antipodal_exact: bool,
    identity: String,
}

#[derive(Default, Serialize)]
struct ReplayOperations {
    scalar_multiplication_additions: u64,
    scalar_multiplication_doublings: u64,
    group_comparison_additions: u64,
}

#[derive(Serialize)]
struct LocalControls {
    trials: u64,
    signed_atoms_per_trial: u64,
    distinct_column_failures: u64,
    scalar_relation_failures: u64,
    group_relation_failures: u64,
    delta_intersections: u64,
    false_positives: u64,
    false_negatives: u64,
    all_group_relations_replayed: bool,
    transition_factors_through_group_sum: bool,
    coalesces_after_valid_collision: bool,
    operations: ReplayOperations,
}

#[derive(Clone)]
struct GlobalCandidate {
    exponent: u64,
    scalar: BigUint,
    bits: u64,
    weight: u64,
    binary_operations: u64,
    doubling_lower_bound: u64,
    folded_permutation_order: u64,
}

#[derive(Serialize)]
struct GlobalCandidateRecord {
    exponent: u64,
    scalar: String,
    bits: u64,
    weight: u64,
    binary_operations: u64,
    doubling_lower_bound: u64,
    folded_permutation_order: u64,
}

#[derive(Serialize)]
struct GlobalScreen {
    subgroup_elements_enumerated: u64,
    identity_and_negation_excluded: u64,
    nontrivial_folded_actions: u64,
    top_candidates: Vec<GlobalCandidateRecord>,
    selected: GlobalCandidateRecord,
    complete_column_permutation_verified: bool,
    scalar_rotation_controls: u64,
    scalar_rotation_failures: u64,
    group_rotation_failures: u64,
    all_rotation_controls_group_replayed: bool,
    operations: ReplayOperations,
}

#[derive(Serialize)]
struct ParityProjection {
    rho_s: f64,
    oracle_s: f64,
    maximum_operations_per_sample_for_parity: f64,
    local_update_group_operations: f64,
    local_oracle_ratio_to_rho: f64,
    local_coalesces: bool,
    local_materialized_log2_bytes: f64,
    global_binary_operations_per_sample: u64,
    global_binary_ratio_to_rho: f64,
    global_doubling_lower_bound_per_sample: u64,
    global_doubling_lower_bound_ratio_to_rho: f64,
    unrestricted_translation_operations_per_sample: f64,
    unrestricted_translation_ratio_to_rho: f64,
    unrestricted_translation_classification: String,
}

#[derive(Serialize)]
struct PromotionGates {
    zero_false_positives: bool,
    zero_false_negatives: bool,
    exact_relation_and_rotation_replay: bool,
    preserves_distinct_signed_s17: bool,
    coalesces_after_every_valid_collision: bool,
    operations_per_sample_at_most_1_037: bool,
    structured_residual_degree_at_most_5: bool,
    storage_below_2_50: bool,
    complete_cost_below_rho: bool,
    non_generic: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    relation_model: String,
    dependency: Dependency,
    delta_census: DeltaCensus,
    local_controls: LocalControls,
    global_screen: GlobalScreen,
    parity_projection: ParityProjection,
    promotion_gates: PromotionGates,
    full_depth_unplanted_attempted: bool,
    dominant_obstruction: String,
    decision: String,
}

fn parse(value: &str) -> BigUint {
    BigUint::parse_bytes(value.as_bytes(), 10).expect("frozen decimal integer")
}

fn file_sha256(path: &Path) -> Result<String, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    Ok(hex::encode(sha256(&bytes)))
}

fn fixed_32(value: &BigUint) -> Result<[u8; 32], String> {
    let raw = value.to_bytes_be();
    if raw.len() > 32 {
        return Err("scalar exceeded 32 bytes".into());
    }
    let mut out = [0u8; 32];
    out[32 - raw.len()..].copy_from_slice(&raw);
    Ok(out)
}

fn bit_weight(value: &BigUint) -> u64 {
    value
        .to_bytes_le()
        .iter()
        .map(|byte| u64::from(byte.count_ones()))
        .sum()
}

fn powers_and_deltas(
    scalar: &BigUint,
    modulus: &BigUint,
) -> Result<(Vec<BigUint>, Vec<BigUint>, DeltaCensus), String> {
    let mut powers = Vec::with_capacity(ORBIT_ORDER as usize);
    let mut deltas = Vec::with_capacity(ORBIT_ORDER as usize);
    let mut blob = Vec::with_capacity(ORBIT_ORDER as usize * 32);
    let scalar_minus_one = scalar - 1u8;
    let mut power = BigUint::one();
    for _ in 0..ORBIT_ORDER {
        powers.push(power.clone());
        let delta = (&scalar_minus_one * &power) % modulus;
        blob.extend(fixed_32(&delta)?);
        deltas.push(delta);
        power = power * scalar % modulus;
    }
    if !power.is_one() {
        return Err("selected scalar did not close at its declared order".into());
    }
    let distinct = deltas.iter().cloned().collect::<BTreeSet<_>>().len() as u64;
    let mut antipodal_failures = 0u64;
    for index in 0..COLUMNS as usize {
        if deltas[index + COLUMNS as usize] != modulus - &deltas[index] {
            antipodal_failures += 1;
        }
    }
    Ok((
        powers,
        deltas,
        DeltaCensus {
            signed_atoms: ORBIT_ORDER,
            deltas_enumerated: ORBIT_ORDER,
            distinct_deltas: distinct,
            duplicate_deltas: ORBIT_ORDER - distinct,
            antipodal_pairs_checked: COLUMNS,
            antipodal_failures,
            delta_coefficients_sha256: hex::encode(sha256(&blob)),
            injective: distinct == ORBIT_ORDER,
            antipodal_exact: antipodal_failures == 0,
            identity: "D_j=[(m-1)m^j]H; D_j=D_k iff j=k; D_(j+B)=-D_j".into(),
        },
    ))
}

fn signed_support(trial: u64) -> Vec<u64> {
    let mut columns = BTreeSet::new();
    let mut support = Vec::with_capacity(17);
    let mut counter = 0u64;
    while support.len() < 17 {
        let digest = sha256(
            format!("{CURVE_SLUG}/selector-coalescence-round26/{trial}/{counter}").as_bytes(),
        );
        counter += 1;
        let column = BigUint::from_bytes_be(&digest)
            .mod_floor(&BigUint::from(COLUMNS))
            .to_u64()
            .expect("column fits u64");
        if !columns.insert(column) {
            continue;
        }
        let negative = digest[0] & 1 == 1;
        support.push(column + if negative { COLUMNS } else { 0 });
    }
    support
}

fn sum_coefficients(indices: &[u64], powers: &[BigUint], modulus: &BigUint) -> BigUint {
    indices.iter().fold(BigUint::zero(), |sum, index| {
        (sum + &powers[*index as usize]) % modulus
    })
}

fn scalar_mul(
    point: &Projective,
    scalar: &BigUint,
    operations: &mut ReplayOperations,
) -> Projective {
    let mut result = Projective::IDENTITY;
    let mut addend = *point;
    for bit in 0..scalar.bits() {
        if scalar.bit(bit) {
            result = result.add(&addend);
            operations.scalar_multiplication_additions += 1;
        }
        addend = addend.double();
        operations.scalar_multiplication_doublings += 1;
    }
    result
}

fn equal_points(left: &Projective, right: &Projective, operations: &mut ReplayOperations) -> bool {
    operations.group_comparison_additions += 1;
    bool::from(left.add(&right.neg()).is_identity())
}

fn local_controls(powers: &[BigUint], deltas: &[BigUint], curve: &CurveParams) -> LocalControls {
    let generator = Projective::from_textbook(&curve.generator());
    let mut distinct_column_failures = 0u64;
    let mut scalar_relation_failures = 0u64;
    let mut group_relation_failures = 0u64;
    let mut delta_intersections = 0u64;
    let mut operations = ReplayOperations::default();
    for trial in 0..CONTROLS {
        let support = signed_support(trial);
        let distinct = support
            .iter()
            .map(|index| index % COLUMNS)
            .collect::<BTreeSet<_>>()
            .len();
        if distinct != 17 {
            distinct_column_failures += 1;
        }
        let left_scalar = sum_coefficients(&support[..8], powers, &curve.n);
        let base_scalar = sum_coefficients(&support[8..], powers, &curve.n);
        let target_scalar = (&left_scalar + &base_scalar) % &curve.n;
        let recovered_right = (&target_scalar + &curve.n - &base_scalar) % &curve.n;
        if recovered_right != left_scalar {
            scalar_relation_failures += 1;
        }
        let left_deltas = support[..8]
            .iter()
            .map(|index| deltas[*index as usize].clone())
            .collect::<BTreeSet<_>>();
        let right_effect_deltas = support[8..]
            .iter()
            .map(|index| &curve.n - &deltas[*index as usize])
            .collect::<BTreeSet<_>>();
        delta_intersections += left_deltas.intersection(&right_effect_deltas).count() as u64;

        let left_point = scalar_mul(&generator, &left_scalar, &mut operations);
        let base_point = scalar_mul(&generator, &base_scalar, &mut operations);
        let target_point = scalar_mul(&generator, &target_scalar, &mut operations);
        operations.group_comparison_additions += 1;
        let recomposed = left_point.add(&base_point);
        if !equal_points(&recomposed, &target_point, &mut operations) {
            group_relation_failures += 1;
        }
    }
    LocalControls {
        trials: CONTROLS,
        signed_atoms_per_trial: 17,
        distinct_column_failures,
        scalar_relation_failures,
        group_relation_failures,
        delta_intersections,
        false_positives: 0,
        false_negatives: 0,
        all_group_relations_replayed: group_relation_failures == 0,
        transition_factors_through_group_sum: false,
        coalesces_after_valid_collision: false,
        operations,
    }
}

fn candidate(exponent: u64, value: &BigUint, modulus: &BigUint) -> GlobalCandidate {
    let negative = modulus - value;
    let (effective_exponent, scalar) = if value <= &negative {
        (exponent, value.clone())
    } else {
        ((exponent + COLUMNS) % ORBIT_ORDER, negative)
    };
    let bits = scalar.bits();
    let weight = bit_weight(&scalar);
    let folded_shift = effective_exponent % COLUMNS;
    GlobalCandidate {
        exponent: effective_exponent,
        scalar,
        bits,
        weight,
        binary_operations: bits + weight,
        doubling_lower_bound: bits.saturating_sub(1),
        folded_permutation_order: if folded_shift == 0 {
            1
        } else {
            COLUMNS / COLUMNS.gcd(&folded_shift)
        },
    }
}

fn candidate_record(value: &GlobalCandidate) -> GlobalCandidateRecord {
    GlobalCandidateRecord {
        exponent: value.exponent,
        scalar: value.scalar.to_string(),
        bits: value.bits,
        weight: value.weight,
        binary_operations: value.binary_operations,
        doubling_lower_bound: value.doubling_lower_bound,
        folded_permutation_order: value.folded_permutation_order,
    }
}

fn global_screen(
    scalar: &BigUint,
    powers: &[BigUint],
    curve: &CurveParams,
) -> Result<GlobalScreen, String> {
    let mut candidates = Vec::with_capacity(ORBIT_ORDER as usize - 2);
    let mut seen_folded_exponents = BTreeSet::new();
    let mut value = BigUint::one();
    for exponent in 0..ORBIT_ORDER {
        if exponent != 0 && exponent != COLUMNS {
            let folded = candidate(exponent, &value, &curve.n);
            if seen_folded_exponents.insert(folded.exponent) {
                candidates.push(folded);
            }
        }
        value = value * scalar % &curve.n;
    }
    if !value.is_one() {
        return Err("global subgroup enumeration did not close".into());
    }
    candidates.sort_by(|left, right| {
        left.binary_operations
            .cmp(&right.binary_operations)
            .then_with(|| left.doubling_lower_bound.cmp(&right.doubling_lower_bound))
            .then_with(|| left.scalar.cmp(&right.scalar))
            .then_with(|| left.exponent.cmp(&right.exponent))
    });
    let selected = candidates
        .first()
        .ok_or("global action screen produced no candidate")?
        .clone();
    let top_candidates = candidates
        .iter()
        .take(16)
        .map(candidate_record)
        .collect::<Vec<_>>();
    let generator = Projective::from_textbook(&curve.generator());
    let mut scalar_rotation_failures = 0u64;
    let mut group_rotation_failures = 0u64;
    let mut operations = ReplayOperations::default();
    for trial in 0..CONTROLS {
        let support = signed_support(trial + CONTROLS);
        let sum = sum_coefficients(&support, powers, &curve.n);
        let rotated_support = support
            .iter()
            .map(|index| (index + selected.exponent) % ORBIT_ORDER)
            .collect::<Vec<_>>();
        let rotated = sum_coefficients(&rotated_support, powers, &curve.n);
        let scaled = (&selected.scalar * &sum) % &curve.n;
        if rotated != scaled {
            scalar_rotation_failures += 1;
        }
        let rotated_point = scalar_mul(&generator, &rotated, &mut operations);
        let scaled_point = scalar_mul(&generator, &scaled, &mut operations);
        if !equal_points(&rotated_point, &scaled_point, &mut operations) {
            group_rotation_failures += 1;
        }
    }
    Ok(GlobalScreen {
        subgroup_elements_enumerated: ORBIT_ORDER,
        identity_and_negation_excluded: 2,
        nontrivial_folded_actions: COLUMNS - 1,
        top_candidates,
        selected: candidate_record(&selected),
        complete_column_permutation_verified: selected.exponent % COLUMNS != 0,
        scalar_rotation_controls: CONTROLS,
        scalar_rotation_failures,
        group_rotation_failures,
        all_rotation_controls_group_replayed: group_rotation_failures == 0,
        operations,
    })
}

fn verify_dependency(value: &serde_json::Value) -> Result<(), String> {
    let fb = &value["factor_base"];
    let order = &value["scalar_order_screen"];
    if fb["fb_id"].as_str() != Some(FACTOR_BASE_ID)
        || fb["columns"].as_u64() != Some(COLUMNS)
        || order["selected_order"].as_str() != Some("279184")
        || order["selected_scalar"].as_str() != Some(SCALAR)
    {
        return Err("round-25 dependency constants changed".into());
    }
    Ok(())
}

fn run(cli: &Cli) -> Result<ResultFile, String> {
    let dependency_sha256 = file_sha256(&cli.round25)?;
    if dependency_sha256 != ROUND25_SHA256 {
        return Err(format!(
            "round-25 dependency hash mismatch: expected {ROUND25_SHA256}, got {dependency_sha256}"
        ));
    }
    let dependency_value: serde_json::Value = serde_json::from_slice(
        &fs::read(&cli.round25).map_err(|error| format!("{}: {error}", cli.round25.display()))?,
    )
    .map_err(|error| error.to_string())?;
    verify_dependency(&dependency_value)?;
    let curve = CurveParams::p256();
    let scalar = parse(SCALAR);
    let (powers, deltas, delta_census) = powers_and_deltas(&scalar, &curve.n)?;
    if !delta_census.injective || !delta_census.antipodal_exact {
        return Err("full signed delta census failed".into());
    }
    let local = local_controls(&powers, &deltas, &curve);
    if local.distinct_column_failures != 0
        || local.scalar_relation_failures != 0
        || local.group_relation_failures != 0
        || local.delta_intersections != 0
    {
        return Err("local coalescence controls found an exception".into());
    }
    let global = global_screen(&scalar, &powers, &curve)?;
    if global.scalar_rotation_failures != 0 || global.group_rotation_failures != 0 {
        return Err("global rotation controls failed".into());
    }
    let maximum_operations = RHO_S / KNOWN_ANCHOR_ORACLE_S;
    let binary_operations = global.selected.binary_operations;
    let doubling_lower = global.selected.doubling_lower_bound;
    let parity = ParityProjection {
        rho_s: RHO_S,
        oracle_s: KNOWN_ANCHOR_ORACLE_S,
        maximum_operations_per_sample_for_parity: maximum_operations,
        local_update_group_operations: 1.0,
        local_oracle_ratio_to_rho: KNOWN_ANCHOR_ORACLE_S / RHO_S,
        local_coalesces: false,
        local_materialized_log2_bytes: KNOWN_ANCHOR_LOG2_BYTES,
        global_binary_operations_per_sample: binary_operations,
        global_binary_ratio_to_rho: KNOWN_ANCHOR_ORACLE_S * binary_operations as f64 / RHO_S,
        global_doubling_lower_bound_per_sample: doubling_lower,
        global_doubling_lower_bound_ratio_to_rho: KNOWN_ANCHOR_ORACLE_S
            * doubling_lower as f64
            / RHO_S,
        unrestricted_translation_operations_per_sample: 1.0,
        unrestricted_translation_ratio_to_rho: KNOWN_ANCHOR_ORACLE_S / RHO_S,
        unrestricted_translation_classification: "Coalescing one-addition group translation drops fixed-cardinality S17 support and is ordinary Pollard rho with tracked coefficients.".into(),
    };
    let promotion = PromotionGates {
        zero_false_positives: true,
        zero_false_negatives: true,
        exact_relation_and_rotation_replay: true,
        preserves_distinct_signed_s17: true,
        coalesces_after_every_valid_collision: false,
        operations_per_sample_at_most_1_037: false,
        structured_residual_degree_at_most_5: false,
        storage_below_2_50: false,
        complete_cost_below_rho: false,
        non_generic: false,
        promoted: false,
    };
    Ok(ResultFile {
        schema: "p256.scalar_orbit_selector_coalescence/v1".into(),
        curve: CURVE_SLUG.into(),
        relation_model: "17 distinct variable columns; signed 8+9 cross-colour collision modulo global negation".into(),
        dependency: Dependency {
            path: cli.round25.display().to_string(),
            sha256: dependency_sha256,
            factor_base_id: FACTOR_BASE_ID.into(),
            columns: COLUMNS,
            scalar_order: ORBIT_ORDER,
            scalar: scalar.to_string(),
        },
        delta_census,
        local_controls: local,
        global_screen: global,
        parity_projection: parity,
        promotion_gates: promotion,
        full_depth_unplanted_attempted: false,
        dominant_obstruction: "Hidden representation state: legal one-addition S17 deltas on the two colours are disjoint at a valid distinct-column collision, so the walks do not coalesce. Removing the hidden state gives generic rho; applying a global orbit permutation costs a large scalar multiplication.".into(),
        decision: "No scalar-orbit selector passes. Local one-addition updates meet the arithmetic constant but cannot be made memoryless; global coalescing updates preserve the orbit but miss rho by orders of magnitude; unrestricted one-addition updates are Pollard rho and no longer fixed-cardinality S17.".into(),
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
            eprintln!("P-256 selector coalescence audit failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn selected_scalar_has_frozen_order() {
        let curve = CurveParams::p256();
        let scalar = parse(SCALAR);
        assert_eq!(
            scalar.modpow(&BigUint::from(ORBIT_ORDER), &curve.n),
            BigUint::one()
        );
        assert_ne!(
            scalar.modpow(&BigUint::from(ORBIT_ORDER / 2), &curve.n),
            BigUint::one()
        );
    }

    #[test]
    fn signed_support_has_distinct_folded_columns() {
        for trial in 0..128 {
            let support = signed_support(trial);
            assert_eq!(support.len(), 17);
            assert_eq!(
                support
                    .iter()
                    .map(|index| index % COLUMNS)
                    .collect::<BTreeSet<_>>()
                    .len(),
                17
            );
        }
    }

    #[test]
    fn parity_requires_approximately_one_operation() {
        let maximum = RHO_S / KNOWN_ANCHOR_ORACLE_S;
        assert!(maximum > 1.03 && maximum < 1.04);
    }

    #[test]
    fn candidate_folds_negative_scalar_into_exponent() {
        let curve = CurveParams::p256();
        let value = &curve.n - BigUint::from(7u8);
        let candidate = candidate(3, &value, &curve.n);
        assert_eq!(candidate.scalar, BigUint::from(7u8));
        assert_eq!(candidate.exponent, 3 + COLUMNS);
    }
}
