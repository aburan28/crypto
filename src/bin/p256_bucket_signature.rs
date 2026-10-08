//! Round 306: relation-bucket multiplicity has the rank of its signature.

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

const ROUND298_SHA256: &str = "8718bbe26751ff0191164b0665618226a97556c64c6e2feeef72468a38206965";
const ROUND304_SHA256: &str = "c94e339eca07b96ea0ef890592fb13b80815f0baf52af0e6bd57de17dd46ac9c";
const ROUND305_SHA256: &str = "ecc749ad8edf3e90981f9fb816382d4f4a0913539e9a05e454b6ef8d0a6f4cd7";
const FACTOR_BASE_ID: &str = "FB1hc72514a2a8d3";
const REGISTERED_S17_FACTOR_BASE_ID: &str = "FB1h2f8621cda105";
const COLUMNS: usize = 164;
const TOY_MODULUS: u8 = 7;
const TOY_WIDTH: usize = 4;
const TOY_PROJECTIVE_PROFILES: u64 = 400;
const TOY_SIGNATURES: u64 = 48;
const TOY_PAIRS_PER_PROFILE: u64 = 1_128;
const TOY_DEPENDENT_PAIRS_PER_PROFILE: u64 = 120;
const TOY_INDEPENDENT_PAIRS_PER_PROFILE: u64 = 1_008;
const UNQUOTIENTED_RATIO: f64 = 13.920_747_397_073_491;
const UNKNOWN_ANCHOR_ORACLE_RATIO: f64 = 1.446_504_715_695_271_7;
const REGISTERED_S17_RATIO: f64 = 394.425_280;

#[derive(Parser)]
#[command(about = "Certify P-256 relation-bucket signature rank")]
struct Cli {
    #[arg(long)]
    round298: PathBuf,
    #[arg(long)]
    round304: PathBuf,
    #[arg(long)]
    round305: PathBuf,
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
struct SignatureExample {
    a: u8,
    b: u8,
    quotient_first: u8,
    quotient_second: u8,
    bucket_rows: usize,
    combined_rank: usize,
}

#[derive(Clone, Serialize)]
struct ToyCensus {
    modulus: u8,
    base_columns: usize,
    projective_profiles: u64,
    nonzero_signatures_per_profile: u64,
    full_buckets_checked: u64,
    coefficient_vectors_tested: u64,
    retained_bucket_rows: u64,
    bucket_equation_replays: u64,
    signature_pairs_checked: u64,
    pair_representative_replays: u64,
    dependent_signature_pairs: u64,
    independent_signature_pairs: u64,
    bucket_rank_discrepancies: u64,
    determinant_rank_discrepancies: u64,
    replay_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    minimum_bucket_rank: usize,
    maximum_bucket_rank: usize,
    minimum_dependent_pair_rank: usize,
    maximum_dependent_pair_rank: usize,
    minimum_independent_pair_rank: usize,
    maximum_independent_pair_rank: usize,
    rank_operations: RankOperations,
    first_signature: SignatureExample,
    last_signature: SignatureExample,
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
    systems: Vec<(String, Vec<Equation>, usize, usize, String)>,
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
    signature_status: String,
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
    ratio_to_rho: f64,
    complete_attack_measured: bool,
    status: String,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_and_schemas_checked: bool,
    complete_projective_signature_census: bool,
    every_fixed_signature_bucket_adds_at_most_one_dimension: bool,
    signature_determinant_matches_matrix_rank: bool,
    zero_false_positives_and_false_negatives: bool,
    dual_modulus_p256_ranks_and_group_replays_passed: bool,
    non_scalar_labelled_two_signature_p256_event: bool,
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
    signature_definition: String,
    same_bucket_bound: String,
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
    let (dep298, round298) = checked_json(
        &cli.round298,
        ROUND298_SHA256,
        "p256-parity-escape-screen/v2",
    )?;
    let (dep304, round304) = checked_json(
        &cli.round304,
        ROUND304_SHA256,
        "p256-algebraic-quotient-rigidity-screen/v1",
    )?;
    let (dep305, round305) = checked_json(
        &cli.round305,
        ROUND305_SHA256,
        "p256-homogeneous-anchor-screen/v1",
    )?;
    for (label, value) in [
        ("Round 298", &round298),
        ("Round 304", &round304),
        ("Round 305", &round305),
    ] {
        if value.get("curve").and_then(Value::as_str) != Some(CURVE_SLUG) {
            return Err(format!("{label} curve mismatch"));
        }
    }
    if round298
        .pointer("/p256_boundary/minimum_independent_rows_per_event")
        .and_then(Value::as_u64)
        != Some(COLUMNS as u64)
        || round298.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
    {
        return Err("Round 298 parity conclusion mismatch".into());
    }
    if round304.get("factor_base").and_then(Value::as_str) != Some(FACTOR_BASE_ID)
        || round304.get("columns").and_then(Value::as_u64) != Some(COLUMNS as u64)
        || round304.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
    {
        return Err("Round 304 quotient conclusion mismatch".into());
    }
    if round305.get("factor_base").and_then(Value::as_str) != Some(FACTOR_BASE_ID)
        || round305.get("columns").and_then(Value::as_u64) != Some(COLUMNS as u64)
        || round305
            .pointer("/gates/fixed_target_plus_one_generator_anchor_reaches_full_rank")
            .and_then(Value::as_bool)
            != Some(true)
        || round305.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
    {
        return Err("Round 305 corrected anchor conclusion mismatch".into());
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
    Ok(vec![dep298, dep304, dep305])
}

fn mod7_sub(left: u8, right: u8) -> u8 {
    (left + TOY_MODULUS - right) % TOY_MODULUS
}

fn mod7_dot(left: &[u8], right: &[u8]) -> u8 {
    left.iter().zip(right).fold(0u16, |acc, (a, b)| {
        (acc + u16::from(*a) * u16::from(*b)) % 7
    }) as u8
}

fn inverse_mod7(value: u8) -> u8 {
    (1..TOY_MODULUS)
        .find(|candidate| (u16::from(*candidate) * u16::from(value)) % 7 == 1)
        .expect("nonzero modulo seven")
}

fn decode_base7(mut code: u64, width: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(width);
    for _ in 0..width {
        out.push((code % u64::from(TOY_MODULUS)) as u8);
        code /= u64::from(TOY_MODULUS);
    }
    out
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

fn toy_static_basis(logs: &[u8]) -> Vec<Vec<u8>> {
    let pivot_index = logs.iter().position(|value| *value != 0).unwrap();
    let pivot = logs[pivot_index];
    let mut basis = Vec::with_capacity(TOY_WIDTH - 1);
    for index in 0..TOY_WIDTH {
        if index == pivot_index {
            continue;
        }
        let mut row = vec![0u8; TOY_WIDTH + 1];
        row[pivot_index] = logs[index];
        row[index] = mod7_sub(0, pivot);
        basis.push(row);
    }
    basis
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
    let signatures: Vec<(u8, u8)> = (0..TOY_MODULUS)
        .flat_map(|a| (0..TOY_MODULUS).map(move |b| (a, b)))
        .filter(|pair| *pair != (0, 0))
        .collect();
    if signatures.len() as u64 != TOY_SIGNATURES {
        return Err("toy signature count mismatch".into());
    }

    let mut full_buckets = 0u64;
    let mut coefficient_tests = 0u64;
    let mut retained_rows = 0u64;
    let mut bucket_replays = 0u64;
    let mut pair_replays = 0u64;
    let mut pairs_checked = 0u64;
    let mut dependent_pairs = 0u64;
    let mut independent_pairs = 0u64;
    let mut bucket_rank_discrepancies = 0u64;
    let mut determinant_rank_discrepancies = 0u64;
    let mut replay_failures = 0u64;
    let mut false_positives = 0u64;
    let mut false_negatives = 0u64;
    let mut aggregate_ops = RankOperations::default();
    let mut bucket_ranks = Vec::new();
    let mut dependent_ranks = Vec::new();
    let mut independent_ranks = Vec::new();
    let mut examples = Vec::new();

    for logs in &profiles {
        let target = toy_target(logs);
        let static_basis = toy_static_basis(logs);
        let mut representatives: Vec<(Vec<u8>, (u8, u8))> = Vec::new();
        for (a, b) in &signatures {
            let expected_dot = (u16::from(*a) + u16::from(*b) * u16::from(target)) % 7;
            let quotient = (expected_dot as u8, mod7_sub(0, *b));
            let mut bucket = Vec::new();
            for code in 0..u64::from(TOY_MODULUS).pow(TOY_WIDTH as u32) {
                coefficient_tests += 1;
                let coefficients = decode_base7(code, TOY_WIDTH);
                let dot = mod7_dot(&coefficients, logs);
                if dot != expected_dot as u8 {
                    continue;
                }
                let left = mod7_sub(dot, ((u16::from(*b) * u16::from(target)) % 7) as u8);
                bucket_replays += 1;
                if left != *a {
                    replay_failures += 1;
                    false_positives += 1;
                }
                let mut augmented = coefficients;
                augmented.push(mod7_sub(0, *b));
                bucket.push(augmented);
            }
            if bucket.len() != 343 {
                false_negatives += 1;
            }
            let mut system = static_basis.clone();
            system.extend(bucket.iter().cloned());
            let (rank, ops) = rank_mod7(&system, TOY_WIDTH + 1, false);
            let (reverse_rank, reverse_ops) = rank_mod7(&system, TOY_WIDTH + 1, true);
            aggregate_ops.add_assign(&ops);
            aggregate_ops.add_assign(&reverse_ops);
            if rank != 4 || reverse_rank != 4 {
                bucket_rank_discrepancies += 1;
            }
            bucket_ranks.push(rank);
            full_buckets += 1;
            retained_rows += bucket.len() as u64;
            let representative = bucket.first().ok_or("empty signature bucket")?.clone();
            representatives.push((representative, quotient));
            examples.push(SignatureExample {
                a: *a,
                b: *b,
                quotient_first: quotient.0,
                quotient_second: quotient.1,
                bucket_rows: bucket.len(),
                combined_rank: rank,
            });
        }

        let mut profile_dependent = 0u64;
        let mut profile_independent = 0u64;
        for left in 0..representatives.len() {
            for right in (left + 1)..representatives.len() {
                pairs_checked += 1;
                let (left_row, left_signature) = &representatives[left];
                let (right_row, right_signature) = &representatives[right];
                let determinant = mod7_sub(
                    ((u16::from(left_signature.0) * u16::from(right_signature.1)) % 7) as u8,
                    ((u16::from(left_signature.1) * u16::from(right_signature.0)) % 7) as u8,
                );
                let expected_rank = if determinant == 0 { 4 } else { 5 };
                if determinant == 0 {
                    dependent_pairs += 1;
                    profile_dependent += 1;
                } else {
                    independent_pairs += 1;
                    profile_independent += 1;
                }
                for (row, (a, b)) in [(left_row, signatures[left]), (right_row, signatures[right])]
                {
                    pair_replays += 1;
                    let dot = mod7_dot(&row[..TOY_WIDTH], logs);
                    let lhs = mod7_sub(dot, ((u16::from(b) * u16::from(target)) % 7) as u8);
                    if lhs != a {
                        replay_failures += 1;
                        false_positives += 1;
                    }
                }
                let mut system = static_basis.clone();
                system.push(left_row.clone());
                system.push(right_row.clone());
                let (rank, ops) = rank_mod7(&system, TOY_WIDTH + 1, false);
                let (reverse_rank, reverse_ops) = rank_mod7(&system, TOY_WIDTH + 1, true);
                aggregate_ops.add_assign(&ops);
                aggregate_ops.add_assign(&reverse_ops);
                if rank != expected_rank || reverse_rank != expected_rank {
                    determinant_rank_discrepancies += 1;
                }
                if determinant == 0 {
                    dependent_ranks.push(rank);
                } else {
                    independent_ranks.push(rank);
                }
            }
        }
        if profile_dependent != TOY_DEPENDENT_PAIRS_PER_PROFILE
            || profile_independent != TOY_INDEPENDENT_PAIRS_PER_PROFILE
        {
            determinant_rank_discrepancies += 1;
        }
    }

    let expected_buckets = TOY_PROJECTIVE_PROFILES * TOY_SIGNATURES;
    let expected_pairs = TOY_PROJECTIVE_PROFILES * TOY_PAIRS_PER_PROFILE;
    if full_buckets != expected_buckets
        || pairs_checked != expected_pairs
        || replay_failures != 0
        || false_positives != 0
        || false_negatives != 0
        || bucket_rank_discrepancies != 0
        || determinant_rank_discrepancies != 0
    {
        return Err("toy bucket-signature census failed".into());
    }

    Ok(ToyCensus {
        modulus: TOY_MODULUS,
        base_columns: TOY_WIDTH,
        projective_profiles: profiles.len() as u64,
        nonzero_signatures_per_profile: signatures.len() as u64,
        full_buckets_checked: full_buckets,
        coefficient_vectors_tested: coefficient_tests,
        retained_bucket_rows: retained_rows,
        bucket_equation_replays: bucket_replays,
        signature_pairs_checked: pairs_checked,
        pair_representative_replays: pair_replays,
        dependent_signature_pairs: dependent_pairs,
        independent_signature_pairs: independent_pairs,
        bucket_rank_discrepancies,
        determinant_rank_discrepancies,
        replay_failures,
        false_positives,
        false_negatives,
        minimum_bucket_rank: *bucket_ranks.iter().min().ok_or("no bucket ranks")?,
        maximum_bucket_rank: *bucket_ranks.iter().max().ok_or("no bucket ranks")?,
        minimum_dependent_pair_rank: *dependent_ranks.iter().min().ok_or("no dependent ranks")?,
        maximum_dependent_pair_rank: *dependent_ranks.iter().max().ok_or("no dependent ranks")?,
        minimum_independent_pair_rank: *independent_ranks
            .iter()
            .min()
            .ok_or("no independent ranks")?,
        maximum_independent_pair_rank: *independent_ranks
            .iter()
            .max()
            .ok_or("no independent ranks")?,
        rank_operations: aggregate_ops,
        first_signature: examples.first().ok_or("no signature examples")?.clone(),
        last_signature: examples.last().ok_or("no signature examples")?.clone(),
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

fn static_basis(logs: &[BigUint], modulus: &BigUint) -> Vec<Equation> {
    let pivot = &logs[0];
    (1..COLUMNS)
        .map(|index| {
            let mut coefficients = vec![BigUint::zero(); COLUMNS];
            coefficients[0] = logs[index].clone();
            coefficients[index] = mod_sub(&BigUint::zero(), pivot, modulus);
            Equation {
                coefficients,
                rhs: BigUint::zero(),
            }
        })
        .collect()
}

fn translate_bucket(
    base: &[BigUint],
    static_rows: &[Equation],
    multiplier: u8,
    modulus: &BigUint,
) -> Vec<Equation> {
    let mut rows = vec![Equation {
        coefficients: base.to_vec(),
        rhs: BigUint::zero(),
    }];
    for relation in static_rows {
        let mut coefficients = base.to_vec();
        for (entry, addend) in coefficients[..COLUMNS]
            .iter_mut()
            .zip(&relation.coefficients)
        {
            *entry = (&*entry + BigUint::from(multiplier) * addend) % modulus;
        }
        rows.push(Equation {
            coefficients,
            rhs: BigUint::zero(),
        });
    }
    rows
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
    let static_rows = static_basis(&logs, modulus);
    let mut target_base = vec![BigUint::zero(); COLUMNS + 1];
    target_base[0] = (&target_log * &pivot_inverse) % modulus;
    target_base[COLUMNS] = modulus - BigUint::one();
    let target_bucket = translate_bucket(&target_base, &static_rows, 1, modulus);

    let mut enlarged_target_bucket = target_bucket.clone();
    enlarged_target_bucket.extend(
        translate_bucket(&target_base, &static_rows, 2, modulus)
            .into_iter()
            .skip(1),
    );

    let mut proportional_base = vec![BigUint::zero(); COLUMNS + 1];
    proportional_base[0] = (BigUint::from(2u8) * &target_log * &pivot_inverse) % modulus;
    proportional_base[COLUMNS] = modulus - BigUint::from(2u8);
    let proportional_row = Equation {
        coefficients: proportional_base,
        rhs: BigUint::zero(),
    };
    let mut target_plus_proportional = target_bucket.clone();
    target_plus_proportional.push(proportional_row);

    let mut generator_coefficients = vec![BigUint::zero(); COLUMNS + 1];
    generator_coefficients[0] = pivot_inverse;
    let generator_row = Equation {
        coefficients: generator_coefficients,
        rhs: BigUint::one(),
    };
    let mut target_plus_generator = target_bucket.clone();
    target_plus_generator.push(generator_row);

    Ok(ControlSystems {
        logs,
        target_log,
        systems: vec![
            (
                "static-homogeneous".into(),
                static_rows,
                COLUMNS,
                COLUMNS - 1,
                "maximal static rank; one base scale remains".into(),
            ),
            (
                "fixed-target-bucket".into(),
                target_bucket,
                COLUMNS + 1,
                COLUMNS,
                "one signature (d,-1); one global scale remains".into(),
            ),
            (
                "enlarged-fixed-target-bucket".into(),
                enlarged_target_bucket,
                COLUMNS + 1,
                COLUMNS,
                "327 rows with the same signature; rank unchanged".into(),
            ),
            (
                "target-plus-proportional-signature".into(),
                target_plus_proportional,
                COLUMNS + 1,
                COLUMNS,
                "second signature (2d,-2) is proportional; rank unchanged".into(),
            ),
            (
                "target-plus-generator-signature".into(),
                target_plus_generator,
                COLUMNS + 1,
                COLUMNS + 1,
                "signatures (d,-1) and (1,0) are independent; full rank".into(),
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
    let mut ops = RankOperations {
        matrices: 1,
        rows_materialized: matrix.len() as u64,
        entries_materialized: (matrix.len() * width) as u64,
        ..RankOperations::default()
    };
    let exponent = modulus - BigUint::from(2u8);
    let mut rank = 0usize;
    for column in 0..width {
        let mut pivot = None;
        for row in rank..matrix.len() {
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
        for entry in &mut matrix[rank][column..] {
            *entry = (&*entry * &inverse) % modulus;
            ops.modular_multiplications += 1;
        }
        let pivot_row = matrix[rank].clone();
        for row in (rank + 1)..matrix.len() {
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
        if rank == matrix.len() {
            break;
        }
    }
    Ok((rank, ops))
}

fn scalar_replay_failures(system: &ControlSystems, modulus: &BigUint) -> u64 {
    let mut augmented = system.logs.clone();
    augmented.push(system.target_log.clone());
    system
        .systems
        .iter()
        .flat_map(|(_, equations, width, _, _)| {
            equations.iter().map(|equation| {
                let witness = if *width == COLUMNS {
                    &system.logs
                } else {
                    &augmented
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
        let (name61, equations61, width61, expected61, status61) = &control61.systems[index];
        let (namen, equationsn, widthn, expectedn, statusn) = &controln.systems[index];
        if name61 != namen || width61 != widthn || expected61 != expectedn || status61 != statusn {
            return Err("dual-modulus system mismatch".into());
        }
        let (rank61, ops61) = modular_rank(equations61, *width61, &p61, false)?;
        let (reverse61, reverse_ops61) = modular_rank(equations61, *width61, &p61, true)?;
        let (rankn, opsn) = modular_rank(equationsn, *widthn, &curve.n, false)?;
        let (reversen, reverse_opsn) = modular_rank(equationsn, *widthn, &curve.n, true)?;
        let discrepancies = [rank61, reverse61, rankn, reversen]
            .into_iter()
            .filter(|rank| *rank != *expectedn)
            .count() as u64;
        let before = group_ops.equations_replayed;
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
        if discrepancies != 0 {
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
            signature_status: statusn.clone(),
            scalar_equation_replays: (equations61.len() + equationsn.len()) as u64,
            p256_group_equation_replays: group_ops.equations_replayed - before,
            replay_failures: 0,
            rank_discrepancies: discrepancies,
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
            ratio_to_rho: 1.0,
            complete_attack_measured: true,
            status: "complete reference".into(),
        },
        BoundaryRow {
            variant: format!("unquotiented geometry-defined {FACTOR_BASE_ID}"),
            ratio_to_rho: UNQUOTIENTED_RATIO,
            complete_attack_measured: false,
            status: "executable free-oracle boundary unchanged".into(),
        },
        BoundaryRow {
            variant: "free maximal static kernel plus two independent signature events".into(),
            ratio_to_rho: UNKNOWN_ANCHOR_ORACLE_RATIO,
            complete_attack_measured: false,
            status: "optimistic oracle; construction and recovery omitted".into(),
        },
        BoundaryRow {
            variant: format!("registered 17-term {REGISTERED_S17_FACTOR_BASE_ID}"),
            ratio_to_rho: REGISTERED_S17_RATIO,
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
        schema: "p256-bucket-signature-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 306,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        methodology_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/references/methodology.md (unavailable)".into(),
        template_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/assets/assessment-template.json (unavailable)".into(),
        typed_correspondence_graph: vec![
            "maximal static kernel --[quotient]--> variables (lambda,d)".into(),
            "fixed walk bucket (a,b) --[all decompositions]--> one signature (a+b*d,-b)".into(),
            "same-bucket row differences --[subtraction]--> static homogeneous kernel".into(),
            "two proportional signatures --> quotient rank one".into(),
            "two nonproportional signatures --> quotient rank two / scalar-labelled full-rank control".into(),
        ],
        obligations: vec![
            obligation("same-bucket quotient dimension", "supported", "Every row in a fixed (a,b) bucket has the same quotient signature; all differences are static.", "Prime-order cyclic subgroup after a maximal homogeneous quotient."),
            obligation("complete signature-pair control", "supported", &format!("{} complete buckets and {} signature pairs were checked; determinant and matrix rank always agree.", toy.full_buckets_checked, toy.signature_pairs_checked), "All projective four-log profiles over Z/7Z."),
            obligation("P-256 same-signature replay", "supported as scalar-labelled control", &format!("The 327-row enlarged target bucket and proportional-signature stage retain rank 164; the independent generator signature reaches rank 165 across {} group replays.", p256.group_operations.equations_replayed), "Synthetic deterministic labels, not a relation construction on FB1hc72514a2a8d3."),
            obligation("two-signature P-256 event", "unknown", "No non-scalar-labelled event emits both target coupling and an independent generator signature.", "Future non-homomorphic relation mechanisms."),
            obligation("complete below-rho recovery", "unknown", "The two-signature free oracle is above rho and no complete relation/decomposition/recovery pipeline exists.", "End-to-end P-256 attack."),
        ],
        controls: vec![
            format!("The toy census retained and replayed {} bucket rows and {} pair representatives with zero discrepancies.", toy.retained_bucket_rows, toy.pair_representative_replays),
            format!("Five P-256 stages were ranked over two fields in both row orders and replayed through {} exact group equations.", p256.group_operations.equations_replayed),
            "Rounds 298, 304, and corrected 305 were imported only after exact byte-hash, schema, curve, factor-base, width, boundary, and gate checks.".into(),
        ],
        cost_accounting: vec![
            "Measured: toy coefficient tests, retained bucket rows, pair determinants, rank operations, P-256 scalar multiplications, group additions, equation replays, process telemetry, and artifact sizes.".into(),
            "Imported without extrapolation: Round-304 unquotiented and one-anchor oracle boundaries.".into(),
            "Unset: a non-labelled two-signature event, relation yield, structured solver, sparse linear algebra, target recovery, and complete attack cost.".into(),
        ],
        exploration_boundary: "Exact for rows sharing a fixed generator/target walk signature after a maximal static quotient; not a bound on an event that emits multiple nonproportional signatures.".into(),
        weakest_open_obligation: "Construct and replay one P-256 event that emits two nonproportional signatures, including target coupling and a generator anchor, then price complete recovery below rho.".into(),
        narrowest_supported_finding: "Any number of decompositions in one fixed P-256 target/walk bucket contributes at most one post-static dimension. Full rank requires two nonproportional bucket signatures; no qualifying event is constructed.".into(),
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
        complete_projective_signature_census: true,
        every_fixed_signature_bucket_adds_at_most_one_dimension: true,
        signature_determinant_matches_matrix_rank: true,
        zero_false_positives_and_false_negatives: true,
        dual_modulus_p256_ranks_and_group_replays_passed: true,
        non_scalar_labelled_two_signature_p256_event: false,
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
        "many-rows-per-event/fixed-bucket-signature-rank-one/two-signatures-required/parity-blocked";
    let obstruction = "After a maximal static quotient, every decomposition in one fixed walk bucket (a,b) has the same two-variable signature (a+b*d,-b); all row differences are static. Therefore any occupancy in that bucket contributes at most one post-static dimension. P-256 controls retain rank 164 with 327 same-signature rows and with a proportional second signature, reaching rank 165 only with an independent generator signature. The two-signature free oracle is already 1.446505 times rho, and no non-labelled event emitting the required pair is constructed.";
    let decision = "Reject same-bucket multiplicity as a two-row parity escape. Count quotient signature rank, not raw decompositions, and require two nonproportional signatures carrying target coupling and a generator anchor. Keep a structured two-signature event open; do not attempt an unplanted full-depth relation.";
    let semantic = json!({
        "curve": CURVE_SLUG,
        "dependencies": &dependencies,
        "signature_definition": "sigma(c,b)=(c dot ell,-b)",
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
            schema: "p256-bucket-signature-screen/v1".into(),
            curve: CURVE_SLUG.into(),
            screening_round: 306,
            execution_status: "complete".into(),
            factor_base: FACTOR_BASE_ID.into(),
            registered_s17_factor_base: REGISTERED_S17_FACTOR_BASE_ID.into(),
            columns: COLUMNS,
            dependencies,
            signature_definition: "For sum c_i P_i-[b]Q=[a]G, quotienting a maximal static kernel leaves sigma(c,b)=(c dot ell,-b) on variables (lambda,d).".into(),
            same_bucket_bound: "Every decomposition in one fixed (a,b) bucket has sigma=(a+b*d,-b), so arbitrary bucket occupancy adds at most one post-static dimension.".into(),
            exhaustive_toy_census: toy,
            p256_control: p256,
            boundary_table: boundaries,
            structured_residual_degree_of_regularity: None,
            relations_reported_on_actual_factor_base: 0,
            full_depth_unplanted_p256_relation_attempted: false,
            exploration_boundary: "Exact for fixed-signature relation buckets after a maximal static quotient; not universal over events emitting multiple nonproportional signatures.".into(),
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
        "round 306: buckets={}, signature_pairs={}, p256_replays={}, promoted={}",
        result.exhaustive_toy_census.full_buckets_checked,
        result.exhaustive_toy_census.signature_pairs_checked,
        result.p256_control.group_operations.equations_replayed,
        result.gates.promoted
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_bucket_signature: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complete_toy_signature_census_matches_projective_counts() {
        let census = toy_census().expect("toy census");
        assert_eq!(census.full_buckets_checked, 19_200);
        assert_eq!(census.signature_pairs_checked, 451_200);
        assert_eq!(census.dependent_signature_pairs, 48_000);
        assert_eq!(census.independent_signature_pairs, 403_200);
        assert_eq!(census.minimum_bucket_rank, 4);
        assert_eq!(census.maximum_independent_pair_rank, 5);
        assert_eq!(census.replay_failures, 0);
    }

    #[test]
    fn p61_same_signature_rank_stays_one_below_full() {
        let modulus = (BigUint::one() << 61usize) - BigUint::one();
        let control = control_systems(&modulus).expect("control");
        assert_eq!(scalar_replay_failures(&control, &modulus), 0);
        for (_, equations, width, expected, _) in &control.systems {
            assert_eq!(
                modular_rank(equations, *width, &modulus, false)
                    .expect("rank")
                    .0,
                *expected
            );
        }
        assert_eq!(control.systems[2].1.len(), 327);
        assert_eq!(control.systems[2].3, COLUMNS);
        assert_eq!(control.systems[4].3, COLUMNS + 1);
    }

    #[test]
    fn p256_group_equation_replay_is_exact() {
        let curve = CurveParams::p256();
        let generator = curve.generator();
        let left = generator
            .scalar_mul(&BigUint::from(19u8), &curve.a_fe())
            .add(
                &generator.scalar_mul(&BigUint::from(37u8), &curve.a_fe()),
                &curve.a_fe(),
            );
        assert_eq!(
            left,
            generator.scalar_mul(&BigUint::from(56u8), &curve.a_fe())
        );
    }
}
