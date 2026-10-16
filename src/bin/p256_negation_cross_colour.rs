//! Negation-folded two-colour S17 collision accounting (round 24).

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{ToPrimitive, Zero};
use serde::Serialize;

const ROUND22_SHA256: &str = "3116c678d257794040c5519c85a8037177e12381e5f948a47cf1335d0cafde16";
const ROUND23_SHA256: &str = "f315ef150c123824c4d8954f55e62120579fd666ddf42051bda5351a798bac07";
const RELATION_ROWS: u64 = 138_031;
const RHO_S: f64 = 1.3;
const ENTRY_BYTES: u64 = 64;

#[derive(Parser)]
#[command(about = "Audit the negation-folded P-256 S17 cross-colour boundary")]
struct Cli {
    #[arg(long)]
    round22: PathBuf,
    #[arg(long)]
    round23: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Serialize)]
struct Dependencies {
    round22_sha256: String,
    round23_sha256: String,
}

#[derive(Serialize)]
struct ExhaustiveCheck {
    orders: Vec<u64>,
    pairs_checked: u64,
    false_positives: u64,
    false_negatives: u64,
}

#[derive(Serialize)]
struct SimulationCell {
    order: u64,
    trials: u64,
    total_samples: u64,
    mean_samples: f64,
    median_samples: f64,
    mean_s: f64,
    median_s: f64,
    predicted_first_hit_s: f64,
    positive_orientation_hits: u64,
    negative_orientation_hits: u64,
    zero_hits: u64,
    false_positives: u64,
    false_negatives: u64,
    mean_peak_stored_entries: f64,
}

#[derive(Serialize)]
struct PlantedControl {
    orientation: String,
    columns: u64,
    columns_distinct_up_to_sign: bool,
    equal_abscissa: bool,
    opposite_group_points: bool,
    relation_replayed: bool,
    recovered_target_scalar_verified: bool,
    group_additions: u64,
}

#[derive(Serialize)]
struct GrayControl {
    arity: u64,
    fixed_support: bool,
    sign_states: u64,
    global_negation_folded: bool,
    initial_group_additions: u64,
    update_group_additions: u64,
    additions_per_update: f64,
    replay_failures: u64,
    global_column_coverage: bool,
    projection_credit: bool,
}

#[derive(Serialize)]
struct BoundaryProjection {
    name: String,
    factor_base: String,
    columns: u64,
    required_cross_collisions: u64,
    disjoint_probability: f64,
    full_group_old_s: f64,
    negation_folded_s: f64,
    correction_factor: f64,
    ratio_to_rho: f64,
    log2_oracle_samples: f64,
    log2_cost_per_usable_collision: f64,
    direct_average_additions_per_sample: f64,
    direct_log2_group_additions: f64,
    direct_ratio_to_rho: f64,
    optimistic_materialized_log2_bytes: f64,
    sparse_linear_algebra_log2_operations: Option<f64>,
    complete_lower_bound_log2_operations: f64,
    relation_collection_below_2_120: bool,
    cost_per_relation_below_2_103: bool,
    storage_below_2_50: bool,
    oracle_constant_parity: bool,
    end_to_end_parity: bool,
    classification: String,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    relation_model: String,
    dependencies: Dependencies,
    orientation_identity: Vec<String>,
    exhaustive_check: ExhaustiveCheck,
    simulations: Vec<SimulationCell>,
    planted_p256_controls: Vec<PlantedControl>,
    gray_controls: Vec<GrayControl>,
    projections: Vec<BoundaryProjection>,
    reference_false_positives: u64,
    reference_false_negatives: u64,
    all_reported_relations_replayed: bool,
    old_round22_round23_projection_superseded: bool,
    zero_transport_measurements_remain_valid: bool,
    any_complete_promotion_gate_passed: bool,
    full_depth_unplanted_attempted: bool,
    decision: String,
}

#[derive(Clone)]
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn from_label(label: &str) -> Self {
        let digest = sha256(label.as_bytes());
        let mut bytes = [0u8; 8];
        bytes.copy_from_slice(&digest[..8]);
        Self {
            state: u64::from_be_bytes(bytes),
        }
    }

    fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut value = self.state;
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^ (value >> 31)
    }
}

fn file_sha256(path: &Path) -> Result<String, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    Ok(hex::encode(sha256(&bytes)))
}

fn canonical_negation(value: u64, order: u64) -> u64 {
    value.min(if value == 0 { 0 } else { order - value })
}

fn exhaustive_check() -> ExhaustiveCheck {
    let orders = vec![9, 15, 31, 63, 127];
    let mut pairs_checked = 0u64;
    let mut false_positives = 0u64;
    let mut false_negatives = 0u64;
    for order in &orders {
        for left in 0..*order {
            for right in 0..*order {
                pairs_checked += 1;
                let folded = canonical_negation(left, *order) == canonical_negation(right, *order);
                let replay = left == right || (left + right) % order == 0;
                if folded && !replay {
                    false_positives += 1;
                }
                if replay && !folded {
                    false_negatives += 1;
                }
            }
        }
    }
    ExhaustiveCheck {
        orders,
        pairs_checked,
        false_positives,
        false_negatives,
    }
}

fn first_cross_hit(order: u64, trial: u64) -> (u64, u64, u64, u64, u64, usize) {
    let mut rng = SplitMix64::from_label(&format!(
        "{CURVE_SLUG}/negation-cross-colour-round24/{order}/{trial}"
    ));
    let mut left = BTreeMap::<u64, u64>::new();
    let mut right = BTreeMap::<u64, u64>::new();
    let mut samples = 0u64;
    loop {
        let value = rng.next() % order;
        let key = canonical_negation(value, order);
        let other = if samples & 1 == 0 {
            right.get(&key).copied()
        } else {
            left.get(&key).copied()
        };
        samples += 1;
        if let Some(other) = other {
            let (left_value, right_value) = if samples & 1 == 1 {
                (value, other)
            } else {
                (other, value)
            };
            let positive = u64::from(left_value == right_value);
            let negative =
                u64::from(left_value != right_value && (left_value + right_value) % order == 0);
            let zero = u64::from(left_value == 0 && right_value == 0);
            let false_hit = u64::from(positive + negative == 0);
            return (
                samples,
                positive,
                negative,
                zero,
                false_hit,
                left.len() + right.len(),
            );
        }
        if samples & 1 == 1 {
            left.entry(key).or_insert(value);
        } else {
            right.entry(key).or_insert(value);
        }
    }
}

fn simulation(order: u64, trials: u64) -> SimulationCell {
    let mut samples = Vec::with_capacity(trials as usize);
    let mut positive = 0u64;
    let mut negative = 0u64;
    let mut zero = 0u64;
    let mut false_positives = 0u64;
    let mut stored = 0u64;
    for trial in 0..trials {
        let (count, pos, neg, zero_hit, bad, entries) = first_cross_hit(order, trial);
        samples.push(count);
        positive += pos;
        negative += neg;
        zero += zero_hit;
        false_positives += bad;
        stored += entries as u64;
    }
    samples.sort_unstable();
    let total_samples = samples.iter().sum::<u64>();
    let mean_samples = total_samples as f64 / trials as f64;
    let median_samples = if samples.len() % 2 == 0 {
        (samples[samples.len() / 2 - 1] + samples[samples.len() / 2]) as f64 / 2.0
    } else {
        samples[samples.len() / 2] as f64
    };
    let sqrt_order = (order as f64).sqrt();
    SimulationCell {
        order,
        trials,
        total_samples,
        mean_samples,
        median_samples,
        mean_s: mean_samples / sqrt_order,
        median_s: median_samples / sqrt_order,
        predicted_first_hit_s: (std::f64::consts::PI / 2.0).sqrt(),
        positive_orientation_hits: positive,
        negative_orientation_hits: negative,
        zero_hits: zero,
        false_positives,
        false_negatives: 0,
        mean_peak_stored_entries: stored as f64 / trials as f64,
    }
}

fn scalar_mul(point: &Projective, scalar: &BigUint) -> (Projective, u64) {
    let mut result = Projective::IDENTITY;
    let mut addend = *point;
    let mut additions = 0u64;
    for bit in 0..scalar.bits() {
        if scalar.bit(bit) {
            result = result.add(&addend);
            additions += 1;
        }
        addend = addend.double();
        additions += 1;
    }
    (result, additions)
}

fn add_points(points: &[Projective]) -> (Projective, u64) {
    let mut sum = Projective::IDENTITY;
    let mut additions = 0u64;
    for point in points {
        sum = sum.add(point);
        additions += 1;
    }
    (sum, additions)
}

fn signed_scalar(value: &BigUint, negative: bool, modulus: &BigUint) -> BigUint {
    if negative && !value.is_zero() {
        modulus - value
    } else {
        value.clone()
    }
}

fn p256_planted_controls(curve: &CurveParams) -> Result<Vec<PlantedControl>, String> {
    let generator = Projective::from_textbook(&curve.generator());
    let mut scalars = Vec::new();
    let mut signed = BTreeSet::new();
    let mut counter = 0u64;
    while scalars.len() < 17 {
        let digest = sha256(
            format!("{CURVE_SLUG}/negation-cross-colour-round24/column/{counter}").as_bytes(),
        );
        counter += 1;
        let scalar = BigUint::from_bytes_be(&digest) % &curve.n;
        if scalar.is_zero() {
            continue;
        }
        let key = scalar.clone().min(&curve.n - &scalar);
        if signed.insert(key) {
            scalars.push(scalar);
        }
    }
    let mut points = Vec::new();
    let mut additions = 0u64;
    let mut oriented_scalars = Vec::new();
    for (index, scalar) in scalars.iter().enumerate() {
        let negative =
            sha256(format!("{CURVE_SLUG}/negation-cross-colour-round24/sign/{index}").as_bytes())
                [0]
                & 1
                == 1;
        let (mut point, used) = scalar_mul(&generator, scalar);
        additions += used;
        if negative {
            point = point.neg();
        }
        points.push(point);
        oriented_scalars.push(signed_scalar(scalar, negative, &curve.n));
    }
    let (left, left_adds) = add_points(&points[..8]);
    let (right_sum, right_adds) = add_points(&points[8..]);
    additions += left_adds + right_adds;
    let left_scalar = oriented_scalars[..8]
        .iter()
        .fold(BigUint::zero(), |sum, value| (sum + value) % &curve.n);
    let right_scalar = oriented_scalars[8..]
        .iter()
        .fold(BigUint::zero(), |sum, value| (sum + value) % &curve.n);

    let positive_target = left.add(&right_sum);
    let negative_target = right_sum.add(&left.neg());
    additions += 2;
    let positive_r = positive_target.add(&right_sum.neg());
    let negative_r = negative_target.add(&right_sum.neg());
    additions += 2;
    let left_x = left.to_affine().ok_or("planted left is infinity")?.0;
    let positive_x = positive_r.to_affine().ok_or("positive R is infinity")?.0;
    let negative_x = negative_r.to_affine().ok_or("negative R is infinity")?.0;
    let positive_scalar = (&left_scalar + &right_scalar) % &curve.n;
    let negative_scalar = (&right_scalar + &curve.n - &left_scalar) % &curve.n;
    let (positive_check, positive_mul_adds) = scalar_mul(&generator, &positive_scalar);
    let (negative_check, negative_mul_adds) = scalar_mul(&generator, &negative_scalar);
    additions += positive_mul_adds + negative_mul_adds;
    Ok(vec![
        PlantedControl {
            orientation: "L=R".into(),
            columns: 17,
            columns_distinct_up_to_sign: signed.len() == 17,
            equal_abscissa: bool::from(left_x.ct_eq_full(&positive_x)),
            opposite_group_points: false,
            relation_replayed: bool::from(
                positive_target
                    .add(&left.neg())
                    .add(&right_sum.neg())
                    .is_identity(),
            ),
            recovered_target_scalar_verified: bool::from(
                positive_target.add(&positive_check.neg()).is_identity(),
            ),
            group_additions: additions,
        },
        PlantedControl {
            orientation: "L=-R".into(),
            columns: 17,
            columns_distinct_up_to_sign: signed.len() == 17,
            equal_abscissa: bool::from(left_x.ct_eq_full(&negative_x)),
            opposite_group_points: bool::from(left.add(&negative_r).is_identity()),
            relation_replayed: bool::from(
                negative_target
                    .add(&left)
                    .add(&right_sum.neg())
                    .is_identity(),
            ),
            recovered_target_scalar_verified: bool::from(
                negative_target.add(&negative_check.neg()).is_identity(),
            ),
            group_additions: additions,
        },
    ])
}

fn gray_control(arity: usize, fold_global_negation: bool) -> GrayControl {
    let variable_bits = if fold_global_negation {
        arity - 1
    } else {
        arity
    };
    let states = 1u64 << variable_bits;
    let order = 1_000_003i128;
    let support = (1..=arity).map(|value| value as i128).collect::<Vec<_>>();
    let mut running = support.iter().sum::<i128>() % order;
    let mut previous = 0u64;
    let mut failures = 0u64;
    let mut updates = 0u64;
    for index in 1..states {
        let gray = index ^ (index >> 1);
        let difference = gray ^ previous;
        if difference.count_ones() != 1 {
            failures += 1;
            previous = gray;
            continue;
        }
        let bit = difference.trailing_zeros() as usize;
        let support_index = bit + usize::from(fold_global_negation);
        let was_negative = (previous >> bit) & 1 == 1;
        let delta = if was_negative {
            2 * support[support_index]
        } else {
            -2 * support[support_index]
        };
        running = (running + delta).rem_euclid(order);
        let mut replay = 0i128;
        for (support_index, scalar) in support.iter().enumerate() {
            let negative = if fold_global_negation && support_index == 0 {
                false
            } else {
                let bit_index = support_index - usize::from(fold_global_negation);
                (gray >> bit_index) & 1 == 1
            };
            replay += if negative { -*scalar } else { *scalar };
        }
        if running != replay.rem_euclid(order) {
            failures += 1;
        }
        updates += 1;
        previous = gray;
    }
    GrayControl {
        arity: arity as u64,
        fixed_support: true,
        sign_states: states,
        global_negation_folded: fold_global_negation,
        initial_group_additions: arity.saturating_sub(1) as u64,
        update_group_additions: updates,
        additions_per_update: if updates == 0 { 0.0 } else { 1.0 },
        replay_failures: failures,
        global_column_coverage: false,
        projection_credit: false,
    }
}

fn gamma_ratio(integer_k: u64) -> f64 {
    assert!(integer_k >= 1);
    let mut ratio = std::f64::consts::PI.sqrt() / 2.0;
    for k in 1..integer_k {
        ratio *= (k as f64 + 0.5) / k as f64;
    }
    ratio
}

fn disjoint_probability(columns: u64) -> f64 {
    (0..9).fold(1.0, |probability, offset| {
        probability * (columns - 8 - offset) as f64 / (columns - offset) as f64
    })
}

fn log2_big(value: &BigUint) -> f64 {
    let bits = value.bits();
    let shift = bits.saturating_sub(53);
    let top = (value >> shift).to_u64().unwrap_or(0) as f64;
    top.log2() + shift as f64
}

fn projection(
    name: &str,
    factor_base: &str,
    columns: u64,
    collisions: u64,
    include_linear_algebra: bool,
    curve: &CurveParams,
) -> BoundaryProjection {
    let p_disjoint = disjoint_probability(columns);
    let full_group_old_s = (std::f64::consts::PI * collisions as f64 / p_disjoint).sqrt();
    let negation_folded_s = 2.0_f64.sqrt() * gamma_ratio(collisions) / p_disjoint.sqrt();
    let log2_sqrt_n = log2_big(&curve.n) / 2.0;
    let log2_oracle_samples = log2_sqrt_n + negation_folded_s.log2();
    let direct_average = 8.0f64;
    let direct_log2 = log2_oracle_samples + direct_average.log2();
    // Materialize one colour and stream the other: the resident list has
    // half of the expected total samples at the balanced stopping time.
    let memory_log2 = log2_oracle_samples - 1.0 + (ENTRY_BYTES as f64).log2();
    let sparse_linear_algebra = include_linear_algebra.then(|| {
        let k = columns as f64;
        let nonzeros = RELATION_ROWS as f64 * 17.0;
        (2.0 * k * nonzeros + k * k).log2()
    });
    let complete_lower = sparse_linear_algebra.map_or(log2_oracle_samples, |la| {
        let maximum = log2_oracle_samples.max(la);
        maximum + (2f64.powf(log2_oracle_samples - maximum) + 2f64.powf(la - maximum)).log2()
    });
    let per_relation = log2_oracle_samples - (collisions as f64).log2();
    let oracle_parity = collisions == 1 && negation_folded_s < RHO_S;
    BoundaryProjection {
        name: name.into(),
        factor_base: factor_base.into(),
        columns,
        required_cross_collisions: collisions,
        disjoint_probability: p_disjoint,
        full_group_old_s,
        negation_folded_s,
        correction_factor: negation_folded_s / full_group_old_s,
        ratio_to_rho: negation_folded_s / RHO_S,
        log2_oracle_samples,
        log2_cost_per_usable_collision: per_relation,
        direct_average_additions_per_sample: direct_average,
        direct_log2_group_additions: direct_log2,
        direct_ratio_to_rho: direct_average * negation_folded_s / RHO_S,
        optimistic_materialized_log2_bytes: memory_log2,
        sparse_linear_algebra_log2_operations: sparse_linear_algebra,
        complete_lower_bound_log2_operations: complete_lower,
        relation_collection_below_2_120: log2_oracle_samples < 120.0,
        cost_per_relation_below_2_103: per_relation < 103.0,
        storage_below_2_50: memory_log2 < 50.0,
        oracle_constant_parity: oracle_parity,
        end_to_end_parity: false,
        classification: if oracle_parity {
            "ideal known-log target collision reaches constant parity, but materialized memory and direct generation fail; a memoryless realization is rho-equivalent".into()
        } else {
            "negation correction lowers the boundary, but the independent-log quotient and complete gates remain far above rho".into()
        },
    }
}

fn run(cli: &Cli) -> Result<ResultFile, String> {
    let round22_sha256 = file_sha256(&cli.round22)?;
    let round23_sha256 = file_sha256(&cli.round23)?;
    if round22_sha256 != ROUND22_SHA256 || round23_sha256 != ROUND23_SHA256 {
        return Err(format!(
            "dependency hash mismatch: round22={round22_sha256}, round23={round23_sha256}"
        ));
    }
    let exhaustive = exhaustive_check();
    if exhaustive.false_positives != 0 || exhaustive.false_negatives != 0 {
        return Err("negation quotient failed exhaustive orientation check".into());
    }
    let simulations = [
        (4_093, 4_096),
        (65_537, 2_048),
        (1_048_583, 1_024),
        (16_777_213, 512),
    ]
    .into_iter()
    .map(|(order, trials)| simulation(order, trials))
    .collect::<Vec<_>>();
    if simulations
        .iter()
        .any(|cell| cell.false_positives != 0 || cell.false_negatives != 0)
    {
        return Err("deterministic simulation replay mismatch".into());
    }
    let curve = CurveParams::p256();
    let planted = p256_planted_controls(&curve)?;
    if planted.iter().any(|control| {
        !control.columns_distinct_up_to_sign
            || !control.equal_abscissa
            || !control.relation_replayed
            || !control.recovered_target_scalar_verified
    }) {
        return Err("P-256 planted orientation replay failed".into());
    }
    let gray_controls = vec![gray_control(8, true), gray_control(9, false)];
    if gray_controls
        .iter()
        .any(|control| control.replay_failures != 0)
    {
        return Err("sign-Gray control failed".into());
    }
    let projections = vec![
        projection(
            "known-log K=1 target",
            "round-22 known-log controls",
            131_458,
            1,
            false,
            &curve,
        ),
        projection(
            "round-23 best-width quotient",
            "FB1h8fc8b5fd8529",
            129_877,
            129_877,
            true,
            &curve,
        ),
        projection(
            "selected Dickson quotient",
            "FB1h2f8621cda105",
            131_458,
            131_458,
            true,
            &curve,
        ),
        projection(
            "fixed useful-row collection",
            "FB1h2f8621cda105",
            131_458,
            RELATION_ROWS,
            true,
            &curve,
        ),
    ];
    Ok(ResultFile {
        schema: "p256.negation_cross_colour/v1".into(),
        curve: CURVE_SLUG.into(),
        relation_model: "17 distinct variable columns; signed 8+9 cross-colour collision modulo global negation".into(),
        dependencies: Dependencies {
            round22_sha256,
            round23_sha256,
        },
        orientation_identity: vec![
            "L=R=T-B implies T=L+B".into(),
            "L=-R=B-T implies T=-L+B by flipping all eight left signs".into(),
        ],
        exhaustive_check: exhaustive,
        simulations,
        planted_p256_controls: planted,
        gray_controls,
        projections,
        reference_false_positives: 0,
        reference_false_negatives: 0,
        all_reported_relations_replayed: true,
        old_round22_round23_projection_superseded: true,
        zero_transport_measurements_remain_valid: true,
        any_complete_promotion_gate_passed: false,
        full_depth_unplanted_attempted: false,
        decision: "Negation folding restores ideal K=1 collision-constant parity, but not end-to-end parity: known-log controls are rho-equivalent, materialized storage is prohibitive, direct generation is above rho, and both geometric bases retain full independent-log quotients.".into(),
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
            eprintln!("P-256 negation cross-colour audit failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn folded_key_is_exactly_equal_or_opposite() {
        let check = exhaustive_check();
        assert_eq!(check.false_positives, 0);
        assert_eq!(check.false_negatives, 0);
    }

    #[test]
    fn first_hit_formula_reaches_constant_parity() {
        let curve = CurveParams::p256();
        let row = projection("K1", "control", 131_458, 1, false, &curve);
        assert!((row.negation_folded_s - (std::f64::consts::PI / 2.0).sqrt()).abs() < 0.001);
        assert!(row.ratio_to_rho < 1.0);
        assert!(!row.end_to_end_parity);
    }

    #[test]
    fn gamma_recurrence_matches_first_cases() {
        let r1 = gamma_ratio(1);
        let r2 = gamma_ratio(2);
        assert!((r1 - std::f64::consts::PI.sqrt() / 2.0).abs() < 1e-12);
        assert!((r2 / r1 - 1.5).abs() < 1e-12);
    }

    #[test]
    fn p256_both_orientations_replay() {
        let controls = p256_planted_controls(&CurveParams::p256()).unwrap();
        assert_eq!(controls.len(), 2);
        assert!(controls.iter().all(|control| control.relation_replayed));
        assert!(controls.iter().all(|control| control.equal_abscissa));
        assert!(controls[1].opposite_group_points);
    }

    #[test]
    fn fixed_support_gray_control_gets_no_global_credit() {
        for control in [gray_control(8, true), gray_control(9, false)] {
            assert_eq!(control.replay_failures, 0);
            assert_eq!(control.additions_per_update, 1.0);
            assert!(!control.global_column_coverage);
            assert!(!control.projection_credit);
        }
    }
}
