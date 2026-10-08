//! Round 305: homogeneous factor-base relations conserve one logarithm scale.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::ecc::{CurveParams, Point};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::Serialize;
use serde_json::{json, Value};

const ROUND295_RESULT_SHA256: &str =
    "9550f2bbaceb9e297fb35c12e79480581ca80bd58b03a8493e86ea7c1eda91a5";
const ROUND295_FACTOR_BASE_SHA256: &str =
    "d27516ca40a612ecf3ebabfa8ae04776084c7948110e20f221da438e1e64d8f8";
const ROUND298_SHA256: &str = "8718bbe26751ff0191164b0665618226a97556c64c6e2feeef72468a38206965";
const ROUND299_SHA256: &str = "e6b78d729e785b1831a061cd39e6e57f10c46da8b8e98418321dc72f56752f35";
const ROUND304_SHA256: &str = "c94e339eca07b96ea0ef890592fb13b80815f0baf52af0e6bd57de17dd46ac9c";
const FACTOR_BASE_ID: &str = "FB1hc72514a2a8d3";
const REGISTERED_S17_FACTOR_BASE_ID: &str = "FB1h2f8621cda105";
const COLUMNS: usize = 164;
const TOY_MODULUS: u8 = 7;
const TOY_WIDTH: usize = 4;
const TOY_PROJECTIVE_PROFILES: u64 = 400;
const UNQUOTIENTED_RATIO: f64 = 13.920_747_397_073_491;
const REGISTERED_S17_RATIO: f64 = 394.425_280;
const UNKNOWN_ANCHOR_ORACLE_RATIO: f64 = 1.446_504_715_695_271_7;
const UNKNOWN_ANCHOR_DIRECT_RATIO: f64 = 11.572_037_725_562_172;
const KNOWN_ANCHOR_ORACLE_RATIO: f64 = 0.964_336_477_130_181;
const KNOWN_ANCHOR_DIRECT_RATIO: f64 = 7.714_691_817_041_447;

#[derive(Parser)]
#[command(about = "Certify homogeneous factor-base anchor conservation on P-256")]
struct Cli {
    #[arg(long)]
    round295_result: PathBuf,
    #[arg(long)]
    round295_factor_base: PathBuf,
    #[arg(long)]
    round298: PathBuf,
    #[arg(long)]
    round299: PathBuf,
    #[arg(long)]
    round304: PathBuf,
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
struct RankOperations {
    matrices: u64,
    rows_materialized: u64,
    entries_materialized: u64,
    pivot_row_scans: u64,
    pivot_swaps: u64,
    modular_inversions: u64,
    modular_multiplications: u64,
    modular_subtractions: u64,
}

impl RankOperations {
    fn add_assign(&mut self, other: &Self) {
        self.matrices += other.matrices;
        self.rows_materialized += other.rows_materialized;
        self.entries_materialized += other.entries_materialized;
        self.pivot_row_scans += other.pivot_row_scans;
        self.pivot_swaps += other.pivot_swaps;
        self.modular_inversions += other.modular_inversions;
        self.modular_multiplications += other.modular_multiplications;
        self.modular_subtractions += other.modular_subtractions;
    }
}

#[derive(Clone, Default, Serialize)]
struct GroupOperations {
    points_constructed: u64,
    equations_replayed: u64,
    scalar_multiplications: u64,
    scalar_bits_processed: u64,
    scalar_one_bits: u64,
    point_additions: u64,
}

#[derive(Clone, Serialize)]
struct ToyProfileExample {
    logs: Vec<u8>,
    target_log: u8,
    static_rows: usize,
    fixed_target_rows: usize,
    static_rank: usize,
    fixed_target_augmented_rank: usize,
    static_plus_one_anchored_rank: usize,
    static_plus_two_anchored_rank: usize,
    fixed_target_plus_one_generator_anchor_rank: usize,
    anchored_right_hand_sides: Vec<u8>,
}

#[derive(Clone, Serialize)]
struct ToyCensus {
    modulus: u8,
    base_columns: usize,
    raw_nonzero_log_vectors: u64,
    projective_profiles: u64,
    coefficient_vectors_tested: u64,
    static_rows_retained: u64,
    fixed_target_rows_retained: u64,
    fixed_target_difference_rows: u64,
    anchored_rows_replayed: u64,
    cyclic_group_equations_replayed: u64,
    replay_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    rank_discrepancies: u64,
    minimum_static_rank: usize,
    maximum_static_rank: usize,
    minimum_fixed_target_augmented_rank: usize,
    maximum_fixed_target_augmented_rank: usize,
    minimum_one_anchored_rank: usize,
    maximum_one_anchored_rank: usize,
    minimum_two_anchored_rank: usize,
    maximum_two_anchored_rank: usize,
    minimum_fixed_target_plus_one_generator_anchor_rank: usize,
    maximum_fixed_target_plus_one_generator_anchor_rank: usize,
    rank_operations: RankOperations,
    first_profile: ToyProfileExample,
    last_profile: ToyProfileExample,
}

#[derive(Clone)]
struct Equation {
    coefficients: Vec<BigUint>,
    rhs: BigUint,
}

#[derive(Clone)]
struct ControlSystems {
    logs: Vec<BigUint>,
    target_log: BigUint,
    systems: Vec<(String, Vec<Equation>, usize, usize)>,
}

#[derive(Clone, Serialize)]
struct StageProfile {
    variant: String,
    rows: usize,
    columns: usize,
    expected_rank: usize,
    rank_mod_2_61_minus_1: usize,
    reversed_rank_mod_2_61_minus_1: usize,
    rank_mod_p256_n: usize,
    reversed_rank_mod_p256_n: usize,
    declared_nullity: usize,
    scalar_equation_replays: u64,
    p256_group_equation_replays: u64,
    replay_failures: u64,
    rank_discrepancies: u64,
    materialized_bytes_mod_p256_n: u64,
    operations: RankOperations,
}

#[derive(Clone, Serialize)]
struct PointReceipt {
    infinity: bool,
    x: Option<String>,
    y: Option<String>,
}

#[derive(Clone, Serialize)]
struct P256Control {
    scalar_labelled_control: bool,
    base_columns: usize,
    augmented_columns: usize,
    witness_log_vector_sha256: String,
    witness_point_vector_sha256: String,
    target_log_sha256: String,
    target_point: PointReceipt,
    stage_profiles: Vec<StageProfile>,
    aggregate_rank_operations: RankOperations,
    group_operations: GroupOperations,
    group_replay_failures: u64,
    all_points_on_curve: bool,
}

#[derive(Clone, Serialize)]
struct BoundaryRow {
    variant: String,
    free_oracle_ratio_to_rho: f64,
    direct_measured_ratio_to_rho: Option<f64>,
    complete_attack_measured: bool,
    status: String,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_and_schemas_checked: bool,
    complete_projective_toy_census: bool,
    zero_false_positives_and_false_negatives: bool,
    static_rank_at_most_k_minus_1: bool,
    fixed_target_retains_global_scale: bool,
    two_post_static_equations_necessary_and_sufficient_in_controls: bool,
    fixed_target_plus_one_generator_anchor_reaches_full_rank: bool,
    dual_modulus_p256_ranks_and_group_replays_passed: bool,
    non_scalar_labelled_two_row_p256_event: bool,
    complete_relation_decomposition_and_recovery: bool,
    structured_residual_degree_at_most_5: bool,
    parity_at_or_below_rho: bool,
    complete_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    projected_storage_below_2_50: bool,
    discarded_probabilistic_branches_counted_as_exhaustive: bool,
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
    columns: usize,
    dependencies: Vec<Dependency>,
    theorem: String,
    target_extension: String,
    exhaustive_toy_census: ToyCensus,
    p256_control: P256Control,
    boundary_table: Vec<BoundaryRow>,
    structured_residual_degree_of_regularity: Option<u64>,
    relations_reported_on_actual_factor_base: u64,
    full_depth_unplanted_p256_relation_attempted: bool,
    exploration_boundary: String,
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
struct TransferAssessment {
    schema: String,
    curve: String,
    screening_round: u64,
    skill_profile: String,
    required_companion_resources_available: bool,
    methodology_resource: String,
    template_resource: String,
    typed_correspondence_graph: Vec<String>,
    obligations: Vec<Obligation>,
    controls: Vec<String>,
    cost_accounting: Vec<String>,
    exploration_boundary: String,
    weakest_open_obligation: String,
    narrowest_supported_finding: String,
    semantic_evidence_sha256: String,
    assessment_json_bytes: u64,
}

fn checked_json(
    path: &Path,
    expected_hash: &str,
    schema: &str,
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
    if value.get("schema").and_then(Value::as_str) != Some(schema) {
        return Err(format!("{} schema mismatch", path.display()));
    }
    Ok((
        Dependency {
            path: path.display().to_string(),
            bytes: bytes.len() as u64,
            sha256: digest,
            schema: schema.into(),
        },
        value,
    ))
}

fn require_ratio(value: &Value, pointer: &str, expected: f64, label: &str) -> Result<(), String> {
    let actual = value
        .pointer(pointer)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("{label} missing at {pointer}"))?;
    if (actual - expected).abs() > 1e-12 {
        return Err(format!("{label} is {actual}, expected {expected}"));
    }
    Ok(())
}

fn dependency_checks(cli: &Cli) -> Result<Vec<Dependency>, String> {
    let (dep295, round295) = checked_json(
        &cli.round295_result,
        ROUND295_RESULT_SHA256,
        "p256-dickson-union-screen/v1",
    )?;
    let (dep_fb, factor_base) = checked_json(
        &cli.round295_factor_base,
        ROUND295_FACTOR_BASE_SHA256,
        "ecbench.factor_base_dump/v1-wide",
    )?;
    let (dep298, round298) = checked_json(
        &cli.round298,
        ROUND298_SHA256,
        "p256-parity-escape-screen/v2",
    )?;
    let (dep299, round299) = checked_json(
        &cli.round299,
        ROUND299_SHA256,
        "p256-transfer-anchor-screen/v1",
    )?;
    let (dep304, round304) = checked_json(
        &cli.round304,
        ROUND304_SHA256,
        "p256-algebraic-quotient-rigidity-screen/v1",
    )?;

    for (label, value) in [
        ("Round 295", &round295),
        ("Round 298", &round298),
        ("Round 299", &round299),
        ("Round 304", &round304),
    ] {
        if value.get("curve").and_then(Value::as_str) != Some(CURVE_SLUG) {
            return Err(format!("{label} curve mismatch"));
        }
    }
    if round295.pointer("/selected/fb_id").and_then(Value::as_str) != Some(FACTOR_BASE_ID)
        || round295
            .pointer("/selected/selection/union_columns")
            .and_then(Value::as_u64)
            != Some(COLUMNS as u64)
        || round295
            .pointer("/selected/factor_base_json_sha256")
            .and_then(Value::as_str)
            != Some(ROUND295_FACTOR_BASE_SHA256)
    {
        return Err("Round 295 selected factor-base mismatch".into());
    }
    if factor_base.pointer("/curve/slug").and_then(Value::as_str) != Some(CURVE_SLUG)
        || factor_base
            .pointer("/factor_base/fb_id")
            .and_then(Value::as_str)
            != Some(FACTOR_BASE_ID)
        || factor_base
            .pointer("/factor_base/columns")
            .and_then(Value::as_u64)
            != Some(COLUMNS as u64)
        || factor_base
            .pointer("/factor_base/signed_points")
            .and_then(Value::as_u64)
            != Some((2 * COLUMNS) as u64)
        || factor_base
            .get("points")
            .and_then(Value::as_array)
            .map(Vec::len)
            != Some(2 * COLUMNS)
    {
        return Err("Round 295 factor-base artifact mismatch".into());
    }
    if round298
        .pointer("/p256_boundary/columns")
        .and_then(Value::as_u64)
        != Some(COLUMNS as u64)
        || round298
            .pointer("/p256_boundary/minimum_independent_rows_per_event")
            .and_then(Value::as_u64)
            != Some(COLUMNS as u64)
        || round298.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
    {
        return Err("Round 298 parity boundary mismatch".into());
    }
    if round299.get("factor_base").and_then(Value::as_str) != Some(FACTOR_BASE_ID)
        || round299
            .pointer("/gates/independent_anchor_eliminated_by_homomorphic_transfer")
            .and_then(Value::as_bool)
            != Some(false)
        || round299.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
    {
        return Err("Round 299 anchor conclusion mismatch".into());
    }
    if round304.get("factor_base").and_then(Value::as_str) != Some(FACTOR_BASE_ID)
        || round304.get("columns").and_then(Value::as_u64) != Some(COLUMNS as u64)
        || round304.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
    {
        return Err("Round 304 quotient conclusion mismatch".into());
    }
    require_ratio(
        &round304,
        "/boundary_table/1/free_oracle_ratio_to_rho",
        UNQUOTIENTED_RATIO,
        "Round 304 unquotiented boundary",
    )?;
    require_ratio(
        &round304,
        "/boundary_table/2/free_oracle_ratio_to_rho",
        UNKNOWN_ANCHOR_ORACLE_RATIO,
        "Round 304 unknown-anchor boundary",
    )?;
    require_ratio(
        &round304,
        "/boundary_table/3/free_oracle_ratio_to_rho",
        KNOWN_ANCHOR_ORACLE_RATIO,
        "Round 304 known-anchor boundary",
    )?;
    Ok(vec![dep295, dep_fb, dep298, dep299, dep304])
}

fn mod7_dot(left: &[u8], right: &[u8]) -> u8 {
    left.iter().zip(right).fold(0u16, |acc, (a, b)| {
        (acc + u16::from(*a) * u16::from(*b)) % u16::from(TOY_MODULUS)
    }) as u8
}

fn mod7_sub(left: u8, right: u8) -> u8 {
    (left + TOY_MODULUS - right) % TOY_MODULUS
}

fn decode_base7(mut code: u64, width: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(width);
    for _ in 0..width {
        out.push((code % u64::from(TOY_MODULUS)) as u8);
        code /= u64::from(TOY_MODULUS);
    }
    out
}

fn inverse_mod7(value: u8) -> u8 {
    (1..TOY_MODULUS)
        .find(|candidate| (u16::from(*candidate) * u16::from(value)) % 7 == 1)
        .expect("nonzero modulo seven")
}

fn projective_toy_vectors() -> Vec<Vec<u8>> {
    let mut unique = BTreeSet::new();
    for code in 1..u64::from(TOY_MODULUS).pow(TOY_WIDTH as u32) {
        let mut values = decode_base7(code, TOY_WIDTH);
        let first = values.iter().copied().find(|value| *value != 0).unwrap();
        let inverse = inverse_mod7(first);
        for value in &mut values {
            *value = ((u16::from(*value) * u16::from(inverse)) % 7) as u8;
        }
        unique.insert([values[0], values[1], values[2], values[3]]);
    }
    unique.into_iter().map(|values| values.to_vec()).collect()
}

fn toy_target(logs: &[u8]) -> u8 {
    let weighted = logs.iter().enumerate().fold(0u64, |acc, (index, value)| {
        acc + (index as u64 + 11) * u64::from(*value)
    });
    (weighted % 6 + 1) as u8
}

fn rank_mod7(rows: &[Vec<u8>], width: usize, reverse: bool) -> (usize, RankOperations) {
    let mut matrix = rows.to_vec();
    if reverse {
        matrix.reverse();
    }
    let mut ops = RankOperations {
        matrices: 1,
        rows_materialized: matrix.len() as u64,
        entries_materialized: (matrix.len() * width) as u64,
        ..RankOperations::default()
    };
    let mut rank = 0usize;
    for column in 0..width {
        let mut pivot = None;
        for row in rank..matrix.len() {
            ops.pivot_row_scans += 1;
            if matrix[row][column] != 0 {
                pivot = Some(row);
                break;
            }
        }
        let Some(pivot) = pivot else {
            continue;
        };
        if pivot != rank {
            matrix.swap(pivot, rank);
            ops.pivot_swaps += 1;
        }
        let inverse = inverse_mod7(matrix[rank][column]);
        ops.modular_inversions += 1;
        for entry in &mut matrix[rank][column..width] {
            *entry = ((u16::from(*entry) * u16::from(inverse)) % 7) as u8;
            ops.modular_multiplications += 1;
        }
        let pivot_row = matrix[rank].clone();
        for row in (rank + 1)..matrix.len() {
            let factor = matrix[row][column];
            if factor == 0 {
                continue;
            }
            for entry_column in column..width {
                let product = ((u16::from(factor) * u16::from(pivot_row[entry_column])) % 7) as u8;
                matrix[row][entry_column] = mod7_sub(matrix[row][entry_column], product);
                ops.modular_multiplications += 1;
                ops.modular_subtractions += 1;
            }
        }
        rank += 1;
        if rank == matrix.len() {
            break;
        }
    }
    (rank, ops)
}

fn toy_census() -> Result<ToyCensus, String> {
    let profiles = projective_toy_vectors();
    if profiles.len() as u64 != TOY_PROJECTIVE_PROFILES {
        return Err("toy projective profile count mismatch".into());
    }
    let mut coefficient_tests = 0u64;
    let mut static_rows_total = 0u64;
    let mut fixed_rows_total = 0u64;
    let mut difference_rows_total = 0u64;
    let mut anchored_rows_total = 0u64;
    let mut replays = 0u64;
    let mut replay_failures = 0u64;
    let mut false_positives = 0u64;
    let mut false_negatives = 0u64;
    let mut rank_discrepancies = 0u64;
    let mut aggregate_ops = RankOperations::default();
    let mut static_ranks = Vec::new();
    let mut fixed_ranks = Vec::new();
    let mut one_ranks = Vec::new();
    let mut two_ranks = Vec::new();
    let mut target_plus_anchor_ranks = Vec::new();
    let mut examples = Vec::new();

    for logs in &profiles {
        let target = toy_target(logs);
        let mut static_rows = Vec::new();
        let mut fixed_base_rows = Vec::new();
        let mut fixed_augmented_rows = Vec::new();
        for code in 0..u64::from(TOY_MODULUS).pow(TOY_WIDTH as u32) {
            let coefficients = decode_base7(code, TOY_WIDTH);
            let value = mod7_dot(&coefficients, logs);
            coefficient_tests += 1;
            if value == 0 {
                replays += 1;
                static_rows.push(coefficients.clone());
            }
            if value == target {
                replays += 1;
                fixed_base_rows.push(coefficients.clone());
                let mut augmented = coefficients;
                augmented.push(TOY_MODULUS - 1);
                if !(mod7_dot(&augmented[..TOY_WIDTH], logs) + 6 * target).is_multiple_of(7) {
                    replay_failures += 1;
                    false_positives += 1;
                }
                fixed_augmented_rows.push(augmented);
            }
        }
        if static_rows.len() != 343 || fixed_base_rows.len() != 343 {
            false_negatives += 1;
        }
        static_rows_total += static_rows.len() as u64;
        fixed_rows_total += fixed_base_rows.len() as u64;

        let first_fixed = fixed_base_rows.first().ok_or("no toy fixed-target row")?;
        let differences: Vec<Vec<u8>> = fixed_base_rows
            .iter()
            .skip(1)
            .map(|row| {
                row.iter()
                    .zip(first_fixed)
                    .map(|(left, right)| mod7_sub(*left, *right))
                    .collect()
            })
            .collect();
        for row in &differences {
            replays += 1;
            if mod7_dot(row, logs) != 0 {
                replay_failures += 1;
                false_positives += 1;
            }
        }
        difference_rows_total += differences.len() as u64;

        let static_lifted: Vec<Vec<u8>> = static_rows
            .iter()
            .map(|row| {
                let mut lifted = row.clone();
                lifted.push(0);
                lifted
            })
            .collect();
        let mut anchors: Vec<(Vec<u8>, u8, (u8, u8))> = Vec::new();
        'candidate: for b in 1..TOY_MODULUS {
            for code in 1..u64::from(TOY_MODULUS).pow(TOY_WIDTH as u32) {
                let coefficients = decode_base7(code, TOY_WIDTH);
                let projection = (mod7_dot(&coefficients, logs), b);
                if let Some((_, _, first_projection)) = anchors.first() {
                    let determinant = mod7_sub(
                        ((u16::from(first_projection.0) * u16::from(projection.1)) % 7) as u8,
                        ((u16::from(first_projection.1) * u16::from(projection.0)) % 7) as u8,
                    );
                    if determinant == 0 {
                        continue;
                    }
                }
                let rhs = (u16::from(projection.0) + u16::from(b) * u16::from(target)) % 7;
                let mut augmented = coefficients;
                augmented.push(b);
                anchors.push((augmented, rhs as u8, projection));
                if anchors.len() == 2 {
                    break 'candidate;
                }
            }
        }
        if anchors.len() != 2 {
            return Err("failed to construct two toy anchored rows".into());
        }
        for (row, rhs, _) in &anchors {
            replays += 1;
            anchored_rows_total += 1;
            let left = (u16::from(mod7_dot(&row[..TOY_WIDTH], logs))
                + u16::from(row[TOY_WIDTH]) * u16::from(target))
                % 7;
            if left as u8 != *rhs {
                replay_failures += 1;
                false_positives += 1;
            }
        }

        let mut static_plus_one = static_lifted.clone();
        static_plus_one.push(anchors[0].0.clone());
        let mut static_plus_two = static_plus_one.clone();
        static_plus_two.push(anchors[1].0.clone());
        let mut scale_anchor = None;
        'scale_anchor: for b in 0..TOY_MODULUS {
            for code in 1..u64::from(TOY_MODULUS).pow(TOY_WIDTH as u32) {
                let coefficients = decode_base7(code, TOY_WIDTH);
                let rhs = (u16::from(mod7_dot(&coefficients, logs))
                    + u16::from(b) * u16::from(target))
                    % 7;
                if rhs == 0 {
                    continue;
                }
                let mut augmented = coefficients;
                augmented.push(b);
                scale_anchor = Some((augmented, rhs as u8));
                break 'scale_anchor;
            }
        }
        let scale_anchor = scale_anchor.ok_or("no toy scale anchor")?;
        replays += 1;
        anchored_rows_total += 1;
        let scale_anchor_left = (u16::from(mod7_dot(&scale_anchor.0[..TOY_WIDTH], logs))
            + u16::from(scale_anchor.0[TOY_WIDTH]) * u16::from(target))
            % 7;
        if scale_anchor_left as u8 != scale_anchor.1 {
            replay_failures += 1;
            false_positives += 1;
        }
        let mut fixed_target_plus_anchor = fixed_augmented_rows.clone();
        fixed_target_plus_anchor.push(scale_anchor.0);
        let systems = [
            (&static_rows, TOY_WIDTH, 3usize),
            (&fixed_augmented_rows, TOY_WIDTH + 1, 4),
            (&differences, TOY_WIDTH, 3),
            (&static_plus_one, TOY_WIDTH + 1, 4),
            (&static_plus_two, TOY_WIDTH + 1, 5),
            (&fixed_target_plus_anchor, TOY_WIDTH + 1, 5),
        ];
        let mut ranks = Vec::new();
        for (rows, width, expected) in systems {
            let (rank, ops) = rank_mod7(rows, width, false);
            let (reverse_rank, reverse_ops) = rank_mod7(rows, width, true);
            aggregate_ops.add_assign(&ops);
            aggregate_ops.add_assign(&reverse_ops);
            if rank != expected || reverse_rank != expected {
                rank_discrepancies += 1;
            }
            ranks.push(rank);
        }
        static_ranks.push(ranks[0]);
        fixed_ranks.push(ranks[1]);
        one_ranks.push(ranks[3]);
        two_ranks.push(ranks[4]);
        target_plus_anchor_ranks.push(ranks[5]);
        examples.push(ToyProfileExample {
            logs: logs.clone(),
            target_log: target,
            static_rows: static_rows.len(),
            fixed_target_rows: fixed_base_rows.len(),
            static_rank: ranks[0],
            fixed_target_augmented_rank: ranks[1],
            static_plus_one_anchored_rank: ranks[3],
            static_plus_two_anchored_rank: ranks[4],
            fixed_target_plus_one_generator_anchor_rank: ranks[5],
            anchored_right_hand_sides: anchors.iter().map(|entry| entry.1).collect(),
        });
    }

    if replay_failures != 0
        || false_positives != 0
        || false_negatives != 0
        || rank_discrepancies != 0
    {
        return Err("toy census replay or rank failure".into());
    }
    Ok(ToyCensus {
        modulus: TOY_MODULUS,
        base_columns: TOY_WIDTH,
        raw_nonzero_log_vectors: u64::from(TOY_MODULUS).pow(TOY_WIDTH as u32) - 1,
        projective_profiles: profiles.len() as u64,
        coefficient_vectors_tested: coefficient_tests,
        static_rows_retained: static_rows_total,
        fixed_target_rows_retained: fixed_rows_total,
        fixed_target_difference_rows: difference_rows_total,
        anchored_rows_replayed: anchored_rows_total,
        cyclic_group_equations_replayed: replays,
        replay_failures,
        false_positives,
        false_negatives,
        rank_discrepancies,
        minimum_static_rank: *static_ranks.iter().min().ok_or("no static ranks")?,
        maximum_static_rank: *static_ranks.iter().max().ok_or("no static ranks")?,
        minimum_fixed_target_augmented_rank: *fixed_ranks.iter().min().ok_or("no fixed ranks")?,
        maximum_fixed_target_augmented_rank: *fixed_ranks.iter().max().ok_or("no fixed ranks")?,
        minimum_one_anchored_rank: *one_ranks.iter().min().ok_or("no one-anchor ranks")?,
        maximum_one_anchored_rank: *one_ranks.iter().max().ok_or("no one-anchor ranks")?,
        minimum_two_anchored_rank: *two_ranks.iter().min().ok_or("no two-anchor ranks")?,
        maximum_two_anchored_rank: *two_ranks.iter().max().ok_or("no two-anchor ranks")?,
        minimum_fixed_target_plus_one_generator_anchor_rank: *target_plus_anchor_ranks
            .iter()
            .min()
            .ok_or("no target-plus-anchor ranks")?,
        maximum_fixed_target_plus_one_generator_anchor_rank: *target_plus_anchor_ranks
            .iter()
            .max()
            .ok_or("no target-plus-anchor ranks")?,
        rank_operations: aggregate_ops,
        first_profile: examples.first().ok_or("no first toy profile")?.clone(),
        last_profile: examples.last().ok_or("no last toy profile")?.clone(),
    })
}

fn mod_sub(left: &BigUint, right: &BigUint, modulus: &BigUint) -> BigUint {
    if left >= right {
        left - right
    } else {
        left + modulus - right
    }
}

fn hash_scalar(label: &str, index: usize, modulus: &BigUint) -> BigUint {
    let digest = sha256_hex(format!("round305/{label}/{index}").as_bytes());
    let value = BigUint::parse_bytes(digest.as_bytes(), 16).expect("hex digest");
    (value % (modulus - BigUint::one())) + BigUint::one()
}

fn big_dot(coefficients: &[BigUint], values: &[BigUint], modulus: &BigUint) -> BigUint {
    coefficients
        .iter()
        .zip(values)
        .fold(BigUint::zero(), |acc, (coefficient, value)| {
            (acc + coefficient * value) % modulus
        })
}

fn control_systems(modulus: &BigUint) -> Result<ControlSystems, String> {
    let logs: Vec<BigUint> = (0..COLUMNS)
        .map(|index| hash_scalar("base-log", index, modulus))
        .collect();
    let target_log = hash_scalar("target-log", 0, modulus);
    let pivot = logs[0].clone();
    let pivot_inverse = pivot.modpow(&(modulus - BigUint::from(2u8)), modulus);
    if (&pivot * &pivot_inverse) % modulus != BigUint::one() {
        return Err("control pivot inverse failed".into());
    }

    let mut static_rows = Vec::with_capacity(COLUMNS - 1);
    for index in 1..COLUMNS {
        let mut coefficients = vec![BigUint::zero(); COLUMNS];
        coefficients[0] = logs[index].clone();
        coefficients[index] = mod_sub(&BigUint::zero(), &pivot, modulus);
        static_rows.push(Equation {
            coefficients,
            rhs: BigUint::zero(),
        });
    }

    let mut fixed_base = vec![BigUint::zero(); COLUMNS + 1];
    fixed_base[0] = (&target_log * &pivot_inverse) % modulus;
    fixed_base[COLUMNS] = modulus - BigUint::one();
    let mut fixed_rows = vec![Equation {
        coefficients: fixed_base.clone(),
        rhs: BigUint::zero(),
    }];
    for relation in &static_rows {
        let mut coefficients = fixed_base.clone();
        for (entry, addend) in coefficients[..COLUMNS]
            .iter_mut()
            .zip(&relation.coefficients)
        {
            *entry = (&*entry + addend) % modulus;
        }
        fixed_rows.push(Equation {
            coefficients,
            rhs: BigUint::zero(),
        });
    }

    let static_lifted: Vec<Equation> = static_rows
        .iter()
        .map(|equation| {
            let mut coefficients = equation.coefficients.clone();
            coefficients.push(BigUint::zero());
            Equation {
                coefficients,
                rhs: BigUint::zero(),
            }
        })
        .collect();
    let mut anchor_one_coefficients = vec![BigUint::zero(); COLUMNS + 1];
    anchor_one_coefficients[0] = BigUint::one();
    anchor_one_coefficients[COLUMNS] = BigUint::one();
    let anchor_one = Equation {
        coefficients: anchor_one_coefficients,
        rhs: (&pivot + &target_log) % modulus,
    };
    let mut anchor_two_coefficients = vec![BigUint::zero(); COLUMNS + 1];
    anchor_two_coefficients[0] = BigUint::from(2u8);
    anchor_two_coefficients[COLUMNS] = BigUint::from(3u8);
    let anchor_two = Equation {
        coefficients: anchor_two_coefficients,
        rhs: (BigUint::from(2u8) * &pivot + BigUint::from(3u8) * &target_log) % modulus,
    };
    let mut static_plus_one = static_lifted.clone();
    static_plus_one.push(anchor_one.clone());
    let mut static_plus_two = static_plus_one.clone();
    static_plus_two.push(anchor_two);
    let mut base_scale_anchor_coefficients = vec![BigUint::zero(); COLUMNS + 1];
    base_scale_anchor_coefficients[0] = BigUint::one();
    let base_scale_anchor = Equation {
        coefficients: base_scale_anchor_coefficients,
        rhs: pivot,
    };
    let mut fixed_target_plus_anchor = fixed_rows.clone();
    fixed_target_plus_anchor.push(base_scale_anchor);

    Ok(ControlSystems {
        logs,
        target_log,
        systems: vec![
            (
                "static-homogeneous".into(),
                static_rows.clone(),
                COLUMNS,
                COLUMNS - 1,
            ),
            (
                "fixed-target-augmented".into(),
                fixed_rows,
                COLUMNS + 1,
                COLUMNS,
            ),
            (
                "fixed-target-differences".into(),
                static_rows,
                COLUMNS,
                COLUMNS - 1,
            ),
            (
                "static-plus-one-anchored".into(),
                static_plus_one,
                COLUMNS + 1,
                COLUMNS,
            ),
            (
                "static-plus-two-anchored".into(),
                static_plus_two,
                COLUMNS + 1,
                COLUMNS + 1,
            ),
            (
                "fixed-target-plus-one-generator-anchor".into(),
                fixed_target_plus_anchor,
                COLUMNS + 1,
                COLUMNS + 1,
            ),
        ],
    })
}

fn modular_rank(
    equations: &[Equation],
    width: usize,
    modulus: &BigUint,
    reverse: bool,
) -> Result<(usize, RankOperations), String> {
    let mut matrix: Vec<Vec<BigUint>> = equations
        .iter()
        .map(|equation| equation.coefficients.clone())
        .collect();
    if matrix.iter().any(|row| row.len() != width) {
        return Err("matrix width mismatch".into());
    }
    if reverse {
        matrix.reverse();
    }
    let rows = matrix.len();
    let mut ops = RankOperations {
        matrices: 1,
        rows_materialized: rows as u64,
        entries_materialized: (rows * width) as u64,
        ..RankOperations::default()
    };
    let exponent = modulus - BigUint::from(2u8);
    let mut rank = 0usize;
    for column in 0..width {
        let mut pivot = None;
        for row in rank..rows {
            ops.pivot_row_scans += 1;
            if !matrix[row][column].is_zero() {
                pivot = Some(row);
                break;
            }
        }
        let Some(pivot) = pivot else {
            continue;
        };
        if pivot != rank {
            matrix.swap(pivot, rank);
            ops.pivot_swaps += 1;
        }
        let inverse = matrix[rank][column].modpow(&exponent, modulus);
        ops.modular_inversions += 1;
        if (&matrix[rank][column] * &inverse) % modulus != BigUint::one() {
            return Err("rank pivot inversion failed".into());
        }
        for entry in &mut matrix[rank][column..] {
            *entry = (&*entry * &inverse) % modulus;
            ops.modular_multiplications += 1;
        }
        let pivot_row = matrix[rank].clone();
        for row in (rank + 1)..rows {
            let factor = matrix[row][column].clone();
            if factor.is_zero() {
                continue;
            }
            for entry_column in column..width {
                let product = (&factor * &pivot_row[entry_column]) % modulus;
                matrix[row][entry_column] = mod_sub(&matrix[row][entry_column], &product, modulus);
                ops.modular_multiplications += 1;
                ops.modular_subtractions += 1;
            }
        }
        rank += 1;
        if rank == rows {
            break;
        }
    }
    Ok((rank, ops))
}

fn scalar_replay_failures(system: &ControlSystems, modulus: &BigUint) -> u64 {
    let mut augmented_witness = system.logs.clone();
    augmented_witness.push(system.target_log.clone());
    system
        .systems
        .iter()
        .flat_map(|(_, equations, width, _)| {
            equations.iter().map(|equation| {
                let witness = if *width == COLUMNS {
                    &system.logs
                } else {
                    &augmented_witness
                };
                u64::from(big_dot(&equation.coefficients, witness, modulus) != equation.rhs)
            })
        })
        .sum()
}

fn point_receipt(point: &Point) -> PointReceipt {
    match point {
        Point::Infinity => PointReceipt {
            infinity: true,
            x: None,
            y: None,
        },
        Point::Affine { x, y } => PointReceipt {
            infinity: false,
            x: Some(x.value.to_string()),
            y: Some(y.value.to_string()),
        },
    }
}

fn count_one_bits(value: &BigUint) -> u64 {
    (0..value.bits()).filter(|bit| value.bit(*bit)).count() as u64
}

fn counted_scalar_mul(
    point: &Point,
    scalar: &BigUint,
    curve: &CurveParams,
    ops: &mut GroupOperations,
) -> Point {
    ops.scalar_multiplications += 1;
    ops.scalar_bits_processed += scalar.bits();
    ops.scalar_one_bits += count_one_bits(scalar);
    point.scalar_mul(scalar, &curve.a_fe())
}

fn replay_group_equation(
    equation: &Equation,
    points: &[Point],
    target: &Point,
    curve: &CurveParams,
    ops: &mut GroupOperations,
) -> bool {
    let mut sum = Point::Infinity;
    for (index, coefficient) in equation.coefficients.iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        let point = if index < COLUMNS {
            &points[index]
        } else {
            target
        };
        let term = counted_scalar_mul(point, coefficient, curve, ops);
        sum = sum.add(&term, &curve.a_fe());
        ops.point_additions += 1;
    }
    let rhs = if equation.rhs.is_zero() {
        Point::Infinity
    } else {
        counted_scalar_mul(&curve.generator(), &equation.rhs, curve, ops)
    };
    ops.equations_replayed += 1;
    sum == rhs
}

fn p256_control(curve: &CurveParams) -> Result<P256Control, String> {
    let p61 = (BigUint::one() << 61usize) - BigUint::one();
    let control61 = control_systems(&p61)?;
    let controln = control_systems(&curve.n)?;
    let scalar_failures =
        scalar_replay_failures(&control61, &p61) + scalar_replay_failures(&controln, &curve.n);
    if scalar_failures != 0 {
        return Err("dual-modulus scalar replay failure".into());
    }

    let generator = curve.generator();
    let mut group_ops = GroupOperations::default();
    let points: Vec<Point> = controln
        .logs
        .iter()
        .map(|log| {
            group_ops.points_constructed += 1;
            counted_scalar_mul(&generator, log, curve, &mut group_ops)
        })
        .collect();
    group_ops.points_constructed += 1;
    let target = counted_scalar_mul(&generator, &controln.target_log, curve, &mut group_ops);
    let all_points_on_curve = points.iter().all(|point| curve.is_on_curve(point))
        && curve.is_on_curve(&target)
        && points.iter().all(|point| !matches!(point, Point::Infinity))
        && !matches!(target, Point::Infinity);
    if !all_points_on_curve {
        return Err("P-256 witness point construction failed".into());
    }

    let mut profiles = Vec::new();
    let mut aggregate = RankOperations::default();
    let mut group_replay_failures = 0u64;
    for index in 0..controln.systems.len() {
        let (name61, equations61, width61, expected61) = &control61.systems[index];
        let (namen, equationsn, widthn, expectedn) = &controln.systems[index];
        if name61 != namen || width61 != widthn || expected61 != expectedn {
            return Err("dual-modulus control system mismatch".into());
        }
        let (rank61, ops61) = modular_rank(equations61, *width61, &p61, false)?;
        let (reverse61, reverse_ops61) = modular_rank(equations61, *width61, &p61, true)?;
        let (rankn, opsn) = modular_rank(equationsn, *widthn, &curve.n, false)?;
        let (reversen, reverse_opsn) = modular_rank(equationsn, *widthn, &curve.n, true)?;
        let rank_discrepancies = [rank61, reverse61, rankn, reversen]
            .into_iter()
            .filter(|rank| *rank != *expectedn)
            .count() as u64;
        let before_replays = group_ops.equations_replayed;
        for equation in equationsn {
            if !replay_group_equation(equation, &points, &target, curve, &mut group_ops) {
                group_replay_failures += 1;
            }
        }
        let mut operations = RankOperations::default();
        for counts in [&ops61, &reverse_ops61, &opsn, &reverse_opsn] {
            operations.add_assign(counts);
            aggregate.add_assign(counts);
        }
        if rank_discrepancies != 0 {
            return Err(format!("{namen} rank discrepancy"));
        }
        profiles.push(StageProfile {
            variant: namen.clone(),
            rows: equationsn.len(),
            columns: *widthn,
            expected_rank: *expectedn,
            rank_mod_2_61_minus_1: rank61,
            reversed_rank_mod_2_61_minus_1: reverse61,
            rank_mod_p256_n: rankn,
            reversed_rank_mod_p256_n: reversen,
            declared_nullity: *widthn - *expectedn,
            scalar_equation_replays: (equations61.len() + equationsn.len()) as u64,
            p256_group_equation_replays: group_ops.equations_replayed - before_replays,
            replay_failures: 0,
            rank_discrepancies,
            materialized_bytes_mod_p256_n: (equationsn.len() * widthn * 32) as u64,
            operations,
        });
    }
    if group_replay_failures != 0 {
        return Err("P-256 group replay failure".into());
    }

    let log_strings: Vec<String> = controln.logs.iter().map(ToString::to_string).collect();
    let point_receipts: Vec<PointReceipt> = points.iter().map(point_receipt).collect();
    Ok(P256Control {
        scalar_labelled_control: true,
        base_columns: COLUMNS,
        augmented_columns: COLUMNS + 1,
        witness_log_vector_sha256: sha256_hex(
            &serde_json::to_vec(&log_strings).map_err(|error| error.to_string())?,
        ),
        witness_point_vector_sha256: sha256_hex(
            &serde_json::to_vec(&point_receipts).map_err(|error| error.to_string())?,
        ),
        target_log_sha256: sha256_hex(controln.target_log.to_string().as_bytes()),
        target_point: point_receipt(&target),
        stage_profiles: profiles,
        aggregate_rank_operations: aggregate,
        group_operations: group_ops,
        group_replay_failures,
        all_points_on_curve,
    })
}

fn boundary_table() -> Vec<BoundaryRow> {
    vec![
        BoundaryRow {
            variant: "Pollard rho".into(),
            free_oracle_ratio_to_rho: 1.0,
            direct_measured_ratio_to_rho: Some(1.0),
            complete_attack_measured: true,
            status: "complete reference".into(),
        },
        BoundaryRow {
            variant: format!("unquotiented geometry-defined {FACTOR_BASE_ID}"),
            free_oracle_ratio_to_rho: UNQUOTIENTED_RATIO,
            direct_measured_ratio_to_rho: None,
            complete_attack_measured: false,
            status: "164 independent log classes".into(),
        },
        BoundaryRow {
            variant: "free maximal homogeneous kernel; one scale remains".into(),
            free_oracle_ratio_to_rho: UNKNOWN_ANCHOR_ORACLE_RATIO,
            direct_measured_ratio_to_rho: None,
            complete_attack_measured: false,
            status: "target coupling plus an independent generator anchor still required; construction omitted".into(),
        },
        BoundaryRow {
            variant: "Round-25 unknown-anchor scalar-labelled control".into(),
            free_oracle_ratio_to_rho: UNKNOWN_ANCHOR_ORACLE_RATIO,
            direct_measured_ratio_to_rho: Some(UNKNOWN_ANCHOR_DIRECT_RATIO),
            complete_attack_measured: false,
            status: "registered direct control remains above rho".into(),
        },
        BoundaryRow {
            variant: "known-anchor scalar-labelled control".into(),
            free_oracle_ratio_to_rho: KNOWN_ANCHOR_ORACLE_RATIO,
            direct_measured_ratio_to_rho: Some(KNOWN_ANCHOR_DIRECT_RATIO),
            complete_attack_measured: false,
            status: "ideal oracle below rho; direct construction above rho".into(),
        },
        BoundaryRow {
            variant: format!("registered 17-term {REGISTERED_S17_FACTOR_BASE_ID}"),
            free_oracle_ratio_to_rho: REGISTERED_S17_RATIO,
            direct_measured_ratio_to_rho: None,
            complete_attack_measured: false,
            status: "unchanged free-perfect-oracle comparison".into(),
        },
    ]
}

fn obligation(name: &str, status: &str, evidence: &str, scope: &str) -> Obligation {
    Obligation {
        name: name.into(),
        status: status.into(),
        evidence: evidence.into(),
        scope: scope.into(),
    }
}

fn build_assessment(toy: &ToyCensus, p256: &P256Control) -> TransferAssessment {
    let mut assessment = TransferAssessment {
        schema: "p256-homogeneous-anchor-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 305,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        methodology_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/references/methodology.md (unavailable)".into(),
        template_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/assets/assessment-template.json (unavailable)".into(),
        typed_correspondence_graph: vec![
            format!("{CURVE_SLUG}/factor-base points --[static homogeneous identities]--> ell^perp"),
            "rank-K-1 static kernel --[quotient]--> one projective log scale".into(),
            "fixed target Q --[unanchored decompositions]--> (ell,d)^perp / one global scale".into(),
            "one generator-anchored row --[known a,b]--> one residual unknown".into(),
            "second independent anchored row --[known a,b]--> full scalar-labelled control".into(),
            "fixed-target coupling plus one independent generator anchor --> full scalar-labelled control".into(),
        ],
        obligations: vec![
            obligation("homogeneous static rank bound", "supported", "Every replayable static row annihilates the nonzero factor-base log vector, so rank is at most K-1.", "Any order-n factor base in the cyclic P-256 subgroup."),
            obligation("fixed-target projective bound", "supported", "Every augmented fixed-target row annihilates (ell,d), so rank is at most K and one global scale remains.", "Unanchored decompositions of one P-256 target."),
            obligation("complete finite control", "supported", &format!("All {} projective Z/7Z log profiles and all {} coefficient vectors were checked with zero discrepancies.", toy.projective_profiles, toy.coefficient_vectors_tested), "Cyclic order-seven group with four base logs and one target log."),
            obligation("P-256 algebra replay", "supported as scalar-labelled control", &format!("{} materialized P-256 equations replayed exactly; dual-modulus ranks reach 163, 164, and 165 at the registered stages.", p256.group_operations.equations_replayed), "Synthetic deterministic log-labelled points, not the geometry-defined factor base."),
            obligation("two-equation P-256 event", "unknown", "No non-scalar-labelled event supplying both target coupling and an independent generator scale anchor is constructed.", "Future non-homomorphic relation mechanisms."),
            obligation("complete below-rho recovery", "unknown", "The free one-anchor boundary is above rho and no new relation/decomposition/recovery pipeline exists.", "End-to-end P-256 attack."),
        ],
        controls: vec![
            format!("The toy census replayed {} cyclic-group equations with zero false positives, false negatives, or rank discrepancies.", toy.cyclic_group_equations_replayed),
            format!("Six P-256 stages were ranked over two fields in both row orders; {} group equations replayed with zero failures.", p256.group_operations.equations_replayed),
            "Round 295, 298, 299, and 304 inputs were imported only after exact byte-hash, schema, curve, factor-base, width, and conclusion checks.".into(),
        ],
        cost_accounting: vec![
            "Measured: complete toy coefficient enumeration, retained rows, rank operations, P-256 scalar multiplications, group additions, equation replays, process telemetry, and artifact sizes.".into(),
            "Imported without extrapolation: Round-304 unquotiented, unknown-anchor, and known-anchor oracle/direct ratios.".into(),
            "Unset: discovery of static identities on FB1hc72514a2a8d3, a target-coupling-plus-scale-anchor event, relation yield, structured solver, sparse linear algebra, target recovery, and complete attack cost.".into(),
        ],
        exploration_boundary: "Exact for static homogeneous group identities and unanchored fixed-target decompositions in a prime-order cyclic subgroup; not a lower bound on a mechanism that supplies target coupling and an independent generator scale anchor in one event.".into(),
        weakest_open_obligation: "Construct and replay one P-256 event that supplies both target coupling and an independent known generator right-hand side without secret scalar labels, then price its complete recovery pipeline below rho.".into(),
        narrowest_supported_finding: "Any free homogeneous structure on a 164-column P-256 factor base leaves at least one global log scale, and unanchored fixed-target decompositions remain projective. Full rank needs two post-static equations: target coupling plus one independent generator anchor, or two independent generator-anchored target rows.".into(),
        semantic_evidence_sha256: String::new(),
        assessment_json_bytes: 0,
    };
    let semantic = json!({
        "curve": assessment.curve,
        "typed_correspondence_graph": assessment.typed_correspondence_graph,
        "obligations": assessment.obligations,
        "controls": assessment.controls,
        "cost_accounting": assessment.cost_accounting,
        "exploration_boundary": assessment.exploration_boundary,
        "weakest_open_obligation": assessment.weakest_open_obligation,
        "narrowest_supported_finding": assessment.narrowest_supported_finding,
    });
    assessment.semantic_evidence_sha256 =
        sha256_hex(&serde_json::to_vec(&semantic).expect("assessment semantic JSON"));
    assessment
}

fn build_result(cli: &Cli) -> Result<(ResultReceipt, TransferAssessment), String> {
    let dependencies = dependency_checks(cli)?;
    let curve = CurveParams::p256();
    if curve.h != 1 || curve.n.bits() != 256 || !curve.is_valid_public_point(&curve.generator()) {
        return Err("P-256 subgroup prerequisites failed".into());
    }
    let toy = toy_census()?;
    let p256 = p256_control(&curve)?;
    let boundaries = boundary_table();
    let gates = Gates {
        dependency_hashes_and_schemas_checked: true,
        complete_projective_toy_census: true,
        zero_false_positives_and_false_negatives: true,
        static_rank_at_most_k_minus_1: true,
        fixed_target_retains_global_scale: true,
        two_post_static_equations_necessary_and_sufficient_in_controls: true,
        fixed_target_plus_one_generator_anchor_reaches_full_rank: true,
        dual_modulus_p256_ranks_and_group_replays_passed: true,
        non_scalar_labelled_two_row_p256_event: false,
        complete_relation_decomposition_and_recovery: false,
        structured_residual_degree_at_most_5: false,
        parity_at_or_below_rho: false,
        complete_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        projected_storage_below_2_50: false,
        discarded_probabilistic_branches_counted_as_exhaustive: false,
        promoted: false,
    };
    let assessment = build_assessment(&toy, &p256);
    let classification =
        "homogeneous-factor-base-relations/projective-scale-conservation/target-plus-anchor-required/parity-blocked";
    let obstruction = "Every homogeneous static relation lies in the orthogonal complement of the nonzero factor-base log vector, so even a free maximal kernel has rank only 163 and leaves one scale. Unanchored decompositions of a fixed target likewise determine only the projective class of (ell,d). Full rank after the static quotient needs two post-static equations: one target coupling plus one independent known generator anchor, or two independent generator-anchored target rows. The registered free one-anchor boundary is already 1.446505 times rho, and no non-scalar-labelled event supplying the required pair is constructed.";
    let decision = "Reject static homogeneous factor-base structure as a complete logarithm quotient. Credit at most rank K-1, require both target coupling and an independent generator scale anchor after that quotient, and charge their discovery and full recovery pipeline. Keep a structured two-equation event open; do not attempt an unplanted full-depth relation.";
    let semantic = json!({
        "curve": CURVE_SLUG,
        "dependencies": &dependencies,
        "theorem": "homogeneous relation rows annihilate the nonzero log vector",
        "toy": &toy,
        "p256": &p256,
        "boundary_table": &boundaries,
        "assessment": assessment.semantic_evidence_sha256,
        "gates": &gates,
        "classification": classification,
        "dominant_obstruction": obstruction,
        "decision": decision,
    });
    let semantic_hash =
        sha256_hex(&serde_json::to_vec(&semantic).map_err(|error| error.to_string())?);
    Ok((
        ResultReceipt {
            schema: "p256-homogeneous-anchor-screen/v1".into(),
            curve: CURVE_SLUG.into(),
            screening_round: 305,
            execution_status: "complete".into(),
            factor_base: FACTOR_BASE_ID.into(),
            registered_s17_factor_base: REGISTERED_S17_FACTOR_BASE_ID.into(),
            columns: COLUMNS,
            dependencies,
            theorem: "For P_i=[ell_i]G with nonzero ell, every replayable homogeneous factor-base row h satisfies h dot ell=0; therefore static rank is at most K-1 and one log scale remains.".into(),
            target_extension: "For Q=[d]G, every unanchored decomposition row (c,-1) annihilates (ell,d), so augmented rank is at most K. Starting from static rank K-1, full rank needs two post-static equations: a target coupling plus one independent generator anchor, or two independent generator-anchored target rows.".into(),
            exhaustive_toy_census: toy,
            p256_control: p256,
            boundary_table: boundaries,
            structured_residual_degree_of_regularity: None,
            relations_reported_on_actual_factor_base: 0,
            full_depth_unplanted_p256_relation_attempted: false,
            exploration_boundary: "Exact for homogeneous static identities and unanchored fixed-target decompositions; not universal over an event carrying both target coupling and an independent known generator coefficient.".into(),
            transfer_assessment_semantic_sha256: assessment.semantic_evidence_sha256.clone(),
            gates,
            classification: classification.into(),
            dominant_obstruction: obstruction.into(),
            decision: decision.into(),
            semantic_evidence_sha256: semantic_hash,
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
        "round 305: toy_profiles={}, p256_replays={}, promoted={}",
        result.exhaustive_toy_census.projective_profiles,
        result.p256_control.group_operations.equations_replayed,
        result.gates.promoted
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_homogeneous_anchor: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complete_toy_projective_census_has_registered_ranks() {
        let census = toy_census().expect("toy census");
        assert_eq!(census.projective_profiles, 400);
        assert_eq!(census.minimum_static_rank, 3);
        assert_eq!(census.maximum_static_rank, 3);
        assert_eq!(census.minimum_fixed_target_augmented_rank, 4);
        assert_eq!(census.maximum_two_anchored_rank, 5);
        assert_eq!(
            census.minimum_fixed_target_plus_one_generator_anchor_rank,
            5
        );
        assert_eq!(census.replay_failures, 0);
    }

    #[test]
    fn p61_control_has_one_scale_until_second_anchor() {
        let modulus = (BigUint::one() << 61usize) - BigUint::one();
        let control = control_systems(&modulus).expect("control");
        assert_eq!(scalar_replay_failures(&control, &modulus), 0);
        for (_, equations, width, expected) in &control.systems {
            assert_eq!(
                modular_rank(equations, *width, &modulus, false)
                    .expect("rank")
                    .0,
                *expected
            );
        }
        assert_eq!(control.systems.last().expect("last").3, COLUMNS + 1);
    }

    #[test]
    fn p256_group_equation_replay_is_exact() {
        let curve = CurveParams::p256();
        let generator = curve.generator();
        let left = generator
            .scalar_mul(&BigUint::from(17u8), &curve.a_fe())
            .add(
                &generator.scalar_mul(&BigUint::from(29u8), &curve.a_fe()),
                &curve.a_fe(),
            );
        let right = generator.scalar_mul(&BigUint::from(46u8), &curve.a_fe());
        assert_eq!(left, right);
        assert!(curve.is_on_curve(&left));
    }
}
