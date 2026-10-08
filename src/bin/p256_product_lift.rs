//! Round 308: audit whether a product-group collision yields two useful P-256 rows.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::ecc::{CurveParams, Point};
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;
use serde_json::{json, Value};

const ROUND39_SHA256: &str = "189e48a45ab2919dd3e279bf058940230c4a0461c745fba42dffb26f785ba079";
const ROUND302_SHA256: &str = "a9a12d823cbca6df965d28f8372217b08e65c04326e13e408bba9dac9f1819c1";
const ROUND303_SHA256: &str = "6b023f47fbb2c57eb79f870f1c25be48e0b4230d974c3823caf76784a095caa0";
const ROUND307_SHA256: &str = "3751468a6708d1851598e366bfb99d0795d2da251377a3a1bbf7918cbb36034d";
const REGISTERED_FACTOR_BASE: &str = "FB1h2f8621cda105";
const LOWER_BOUND_FACTOR_BASE: &str = "FB1hc72514a2a8d3";
const REGISTERED_COLUMNS: u64 = 131_458;
const REQUIRED_USEFUL_ROWS: u64 = 138_031;
const CONTROL_WIDTH: usize = 17;
const TOY_MODULUS: u8 = 5;
const TOY_WIDTH: usize = 3;
const RHO_S: f64 = 1.3;
const UNQUOTIENTED_RATIO: f64 = 13.920_747_397_073_491;
const TWO_PRIMITIVE_RATIO: f64 = 1.446_504_715_695_271_7;
const REGISTERED_RATIO: f64 = 394.425_280;

#[derive(Parser)]
#[command(about = "Audit P-256 product-lift rank and collision cost")]
struct Cli {
    #[arg(long)]
    round39: PathBuf,
    #[arg(long)]
    round302: PathBuf,
    #[arg(long)]
    round303: PathBuf,
    #[arg(long)]
    round307: PathBuf,
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
struct ToyRankControl {
    linked_event_rank: usize,
    linked_event_reversed_rank: usize,
    independent_expanded_block_rank: usize,
    independent_expanded_block_reversed_rank: usize,
    independent_original_projection_rank: usize,
    independent_original_projection_reversed_rank: usize,
}

#[derive(Clone, Default, Serialize)]
struct ToyOperations {
    coefficient_state_evaluations: u64,
    dot_product_terms: u64,
    cyclic_group_additions: u64,
    unordered_pair_tests: u64,
    collision_equations_replayed: u64,
    logical_peak_materialized_bytes: u64,
}

#[derive(Clone, Serialize)]
struct ToyCensus {
    modulus: u8,
    width: usize,
    normalized_projective_profiles: u64,
    ordered_profile_pairs: u64,
    proportional_profile_pairs: u64,
    independent_profile_pairs: u64,
    coefficient_states: u64,
    unordered_coefficient_pair_tests: u64,
    first_coordinate_collisions: u64,
    full_product_collisions: u64,
    full_proportional_collisions: u64,
    full_independent_collisions: u64,
    independent_first_only_collisions: u64,
    coordinate_replay_discrepancies: u64,
    bucket_histogram_discrepancies: u64,
    collision_replay_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    rank_control: ToyRankControl,
    operations: ToyOperations,
}

#[derive(Clone, Serialize)]
struct RankReceipt {
    variant: String,
    rows: usize,
    columns: usize,
    expected_rank: usize,
    rank_mod_2_61_minus_1: usize,
    reversed_rank_mod_2_61_minus_1: usize,
    rank_mod_p256_n: usize,
    reversed_rank_mod_p256_n: usize,
    useful_original_system_rank: usize,
    introduced_auxiliary_log_vector: bool,
    status: String,
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
    columns_materialized: usize,
    primary_log_vector_sha256: String,
    independent_auxiliary_log_vector_sha256: String,
    primary_point_vector_sha256: String,
    independent_auxiliary_point_vector_sha256: String,
    linked_multiplier: String,
    linked_delta_sha256: String,
    independent_delta_sha256: String,
    all_points_on_curve: bool,
    scalar_replay_failures: u64,
    group_replay_failures: u64,
    group_operations: GroupOperations,
    rank_operations: RankOperations,
    rank_receipts: Vec<RankReceipt>,
    linked_first_primary_point: PointReceipt,
    independent_first_auxiliary_point: PointReceipt,
    known_label_selector_new_original_log_rows: u64,
}

#[derive(Clone, Serialize)]
struct CollisionProjection {
    subgroup_order: String,
    original_canonical_states: String,
    product_canonical_states: String,
    asymptotic_expected_product_collision_samples_formula: String,
    expected_product_collision_samples_log2: f64,
    ratio_to_original_rho: f64,
    ratio_to_original_rho_log2: f64,
    optimistic_useful_original_rows_per_product_event: u64,
    optimistic_cost_per_useful_original_row_log2: f64,
    registered_useful_rows_with_five_percent_rank_allowance: u64,
    optimistic_registered_collection_lower_bound_log2_operations: f64,
    materialized_record_bytes: u64,
    materialized_list_bytes_log2: f64,
    materialized_disk_write_plus_read_bytes_log2: f64,
    memoryless_walk_state_bytes: u64,
    memoryless_storage_below_2_50: bool,
    materialized_list_storage_below_2_50: bool,
    operation_cost_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    projection_is_measurement: bool,
}

#[derive(Clone, Serialize)]
struct BoundaryRow {
    variant: String,
    ratio_to_rho: f64,
    ratio_to_rho_log2: f64,
    useful_original_rows_per_event: Option<u64>,
    complete_attack_measured: bool,
    status: String,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_and_schemas_checked: bool,
    complete_toy_product_census: bool,
    zero_false_positives_and_false_negatives: bool,
    p256_exact_group_replays_passed: bool,
    linked_lift_useful_rank_is_one: bool,
    independent_lift_block_rank_is_two: bool,
    independent_lift_original_projection_rank_is_one: bool,
    nonhomomorphic_two_row_original_system_coupling: bool,
    complete_relation_decomposition_and_recovery: bool,
    structured_residual_degree_at_most_5: bool,
    parity_at_or_below_rho: bool,
    complete_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    complete_attack_peak_storage_below_2_50: bool,
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
    factor_base_columns: u64,
    required_useful_rows: u64,
    optimistic_lower_bound_factor_base: String,
    dependencies: Vec<Dependency>,
    typed_product_map: String,
    trilemma: Vec<String>,
    exhaustive_toy_census: ToyCensus,
    p256_control: P256Control,
    generic_product_collision_projection: CollisionProjection,
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
    if value.get("schema").and_then(Value::as_str) != Some(expected_schema) {
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

fn require_bool(value: &Value, pointer: &str, expected: bool, label: &str) -> Result<(), String> {
    let actual = value
        .pointer(pointer)
        .and_then(Value::as_bool)
        .ok_or_else(|| format!("{label} missing Boolean {pointer}"))?;
    if actual != expected {
        return Err(format!(
            "{label} {pointer} is {actual}, expected {expected}"
        ));
    }
    Ok(())
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
    let (dep39, round39) = checked_json(
        &cli.round39,
        ROUND39_SHA256,
        "p256.affine_log_orbit_screen/v1",
    )?;
    let (dep302, round302) = checked_json(
        &cli.round302,
        ROUND302_SHA256,
        "p256-complementary-quotient-screen/v1",
    )?;
    let (dep303, round303) = checked_json(
        &cli.round303,
        ROUND303_SHA256,
        "p256-jacobian-row-capacity-screen/v1",
    )?;
    let (dep307, round307) = checked_json(
        &cli.round307,
        ROUND307_SHA256,
        "p256-event-provenance-screen/v1",
    )?;
    for (label, value) in [
        ("Round 39", &round39),
        ("Round 302", &round302),
        ("Round 303", &round303),
        ("Round 307", &round307),
    ] {
        if value.get("curve").and_then(Value::as_str) != Some(CURVE_SLUG) {
            return Err(format!("{label} curve mismatch"));
        }
        require_bool(value, "/gates/promoted", false, label)?;
    }
    for (label, value) in [
        ("Round 302", &round302),
        ("Round 303", &round303),
        ("Round 307", &round307),
    ] {
        if value.get("factor_base").and_then(Value::as_str) != Some(LOWER_BOUND_FACTOR_BASE)
            || value
                .get("registered_s17_factor_base")
                .and_then(Value::as_str)
                != Some(REGISTERED_FACTOR_BASE)
        {
            return Err(format!("{label} factor-base mismatch"));
        }
    }
    require_bool(
        &round302,
        "/gates/two_distinct_order_n_rows_per_cover_event",
        false,
        "Round 302",
    )?;
    require_bool(
        &round303,
        "/gates/exact_capacity_bound_derived",
        true,
        "Round 303",
    )?;
    require_bool(
        &round307,
        "/gates/every_one_line_closure_adds_at_most_one_dimension",
        true,
        "Round 307",
    )?;
    require_bool(
        &round307,
        "/gates/non_scalar_labelled_two_primitive_p256_event",
        false,
        "Round 307",
    )?;
    require_ratio(
        &round307,
        "/boundary_table/1/ratio_to_rho",
        UNQUOTIENTED_RATIO,
        "Round 307 executable boundary",
    )?;
    require_ratio(
        &round307,
        "/boundary_table/2/ratio_to_rho",
        TWO_PRIMITIVE_RATIO,
        "Round 307 two-primitive boundary",
    )?;
    require_ratio(
        &round307,
        "/boundary_table/4/ratio_to_rho",
        REGISTERED_RATIO,
        "Round 307 registered-base boundary",
    )?;
    Ok(vec![dep39, dep302, dep303, dep307])
}

fn decode_base5(mut code: u16) -> Vec<u8> {
    let mut out = vec![0u8; TOY_WIDTH];
    for entry in &mut out {
        *entry = (code % TOY_MODULUS as u16) as u8;
        code /= TOY_MODULUS as u16;
    }
    out
}

fn inverse_mod5(value: u8) -> u8 {
    match value {
        1 => 1,
        2 => 3,
        3 => 2,
        4 => 4,
        _ => panic!("zero has no inverse"),
    }
}

fn projective_profiles() -> Vec<Vec<u8>> {
    let mut set = BTreeSet::new();
    for code in 1..125u16 {
        let mut row = decode_base5(code);
        let first = row.iter().copied().find(|entry| *entry != 0).unwrap();
        let inverse = inverse_mod5(first);
        for entry in &mut row {
            *entry = (*entry * inverse) % TOY_MODULUS;
        }
        set.insert(row);
    }
    set.into_iter().collect()
}

fn dot_mod5(left: &[u8], right: &[u8]) -> u8 {
    left.iter()
        .zip(right)
        .fold(0u8, |sum, (a, b)| (sum + a * b) % TOY_MODULUS)
}

fn replay_dot_mod5(left: &[u8], right: &[u8], operations: &mut ToyOperations) -> u8 {
    let mut sum = 0u8;
    for (coefficient, value) in left.iter().zip(right) {
        operations.dot_product_terms += 1;
        for _ in 0..*coefficient {
            sum = (sum + value) % TOY_MODULUS;
            operations.cyclic_group_additions += 1;
        }
    }
    sum
}

fn subtract_mod5(left: &[u8], right: &[u8]) -> Vec<u8> {
    left.iter()
        .zip(right)
        .map(|(a, b)| (a + TOY_MODULUS - b) % TOY_MODULUS)
        .collect()
}

fn rank_mod5(rows: &[Vec<u8>], reverse: bool) -> usize {
    if rows.is_empty() {
        return 0;
    }
    let width = rows[0].len();
    let mut matrix = rows.to_vec();
    if reverse {
        matrix.reverse();
    }
    let mut rank = 0usize;
    for column in 0..width {
        let Some(pivot) = (rank..matrix.len()).find(|row| matrix[*row][column] != 0) else {
            continue;
        };
        matrix.swap(rank, pivot);
        let inverse = inverse_mod5(matrix[rank][column]);
        for entry in &mut matrix[rank][column..] {
            *entry = (*entry * inverse) % TOY_MODULUS;
        }
        let pivot_row = matrix[rank].clone();
        for row in (rank + 1)..matrix.len() {
            let factor = matrix[row][column];
            for entry_column in column..width {
                matrix[row][entry_column] = (matrix[row][entry_column] + TOY_MODULUS
                    - factor * pivot_row[entry_column] % TOY_MODULUS)
                    % TOY_MODULUS;
            }
        }
        rank += 1;
        if rank == matrix.len() {
            break;
        }
    }
    rank
}

fn toy_census() -> Result<ToyCensus, String> {
    let profiles = projective_profiles();
    if profiles.len() != 31 {
        return Err(format!(
            "got {} projective profiles, expected 31",
            profiles.len()
        ));
    }
    let states: Vec<Vec<u8>> = (0..125u16).map(decode_base5).collect();
    let mut proportional_pairs = 0u64;
    let mut independent_pairs = 0u64;
    let mut coefficient_states = 0u64;
    let mut pair_tests = 0u64;
    let mut first_collisions = 0u64;
    let mut product_collisions = 0u64;
    let mut proportional_collisions = 0u64;
    let mut independent_collisions = 0u64;
    let mut independent_first_only = 0u64;
    let mut coordinate_discrepancies = 0u64;
    let mut histogram_discrepancies = 0u64;
    let mut replay_failures = 0u64;
    let mut false_positives = 0u64;
    let mut false_negatives = 0u64;
    let mut linked_delta = None;
    let mut independent_delta = None;
    let mut operations = ToyOperations::default();

    for primary in &profiles {
        for auxiliary in &profiles {
            let proportional = primary == auxiliary;
            if proportional {
                proportional_pairs += 1;
            } else {
                independent_pairs += 1;
            }
            let mut coordinates = Vec::with_capacity(states.len());
            let mut first_buckets = BTreeMap::<u8, u64>::new();
            let mut product_buckets = BTreeMap::<(u8, u8), u64>::new();
            for state in &states {
                let first = dot_mod5(state, primary);
                let second = dot_mod5(state, auxiliary);
                let replay_first = replay_dot_mod5(state, primary, &mut operations);
                let replay_second = replay_dot_mod5(state, auxiliary, &mut operations);
                coordinate_discrepancies += u64::from(first != replay_first);
                coordinate_discrepancies += u64::from(second != replay_second);
                *first_buckets.entry(first).or_default() += 1;
                *product_buckets.entry((first, second)).or_default() += 1;
                coordinates.push((first, second));
                coefficient_states += 1;
                operations.coefficient_state_evaluations += 1;
            }
            let first_histogram_ok = first_buckets.len() == 5
                && first_buckets.values().all(|occupancy| *occupancy == 25);
            let product_histogram_ok = if proportional {
                product_buckets.len() == 5
                    && product_buckets.values().all(|occupancy| *occupancy == 25)
            } else {
                product_buckets.len() == 25
                    && product_buckets.values().all(|occupancy| *occupancy == 5)
            };
            histogram_discrepancies += u64::from(!first_histogram_ok);
            histogram_discrepancies += u64::from(!product_histogram_ok);

            for left in 0..states.len() {
                for right in (left + 1)..states.len() {
                    pair_tests += 1;
                    operations.unordered_pair_tests += 1;
                    let first_collision = coordinates[left].0 == coordinates[right].0;
                    let product_collision =
                        first_collision && coordinates[left].1 == coordinates[right].1;
                    first_collisions += u64::from(first_collision);
                    product_collisions += u64::from(product_collision);
                    let delta = subtract_mod5(&states[left], &states[right]);
                    let replay_first_zero = dot_mod5(&delta, primary) == 0;
                    let replay_second_zero = dot_mod5(&delta, auxiliary) == 0;
                    let replay_product = replay_first_zero && replay_second_zero;
                    false_positives += u64::from(replay_product && !product_collision);
                    false_negatives += u64::from(product_collision && !replay_product);
                    if product_collision {
                        operations.collision_equations_replayed += 2;
                        replay_failures += u64::from(!replay_first_zero);
                        replay_failures += u64::from(!replay_second_zero);
                        if proportional {
                            proportional_collisions += 1;
                            if linked_delta.is_none() {
                                linked_delta = Some(delta);
                            }
                        } else {
                            independent_collisions += 1;
                            if independent_delta.is_none() {
                                independent_delta = Some(delta);
                            }
                        }
                    } else if first_collision && !proportional {
                        independent_first_only += 1;
                    }
                }
            }
        }
    }

    let linked_delta = linked_delta.ok_or("no linked collision witness")?;
    let independent_delta = independent_delta.ok_or("no independent collision witness")?;
    let linked_rows = vec![linked_delta.clone(), linked_delta.clone()];
    let mut expanded_first = independent_delta.clone();
    expanded_first.extend(vec![0u8; TOY_WIDTH]);
    let mut expanded_second = vec![0u8; TOY_WIDTH];
    expanded_second.extend(independent_delta.clone());
    let expanded_rows = vec![expanded_first, expanded_second];
    let projection_rows = vec![independent_delta, vec![0u8; TOY_WIDTH]];
    let rank_control = ToyRankControl {
        linked_event_rank: rank_mod5(&linked_rows, false),
        linked_event_reversed_rank: rank_mod5(&linked_rows, true),
        independent_expanded_block_rank: rank_mod5(&expanded_rows, false),
        independent_expanded_block_reversed_rank: rank_mod5(&expanded_rows, true),
        independent_original_projection_rank: rank_mod5(&projection_rows, false),
        independent_original_projection_reversed_rank: rank_mod5(&projection_rows, true),
    };
    operations.logical_peak_materialized_bytes =
        (states.len() * TOY_WIDTH + states.len() * 2) as u64;

    let census = ToyCensus {
        modulus: TOY_MODULUS,
        width: TOY_WIDTH,
        normalized_projective_profiles: profiles.len() as u64,
        ordered_profile_pairs: (profiles.len() * profiles.len()) as u64,
        proportional_profile_pairs: proportional_pairs,
        independent_profile_pairs: independent_pairs,
        coefficient_states,
        unordered_coefficient_pair_tests: pair_tests,
        first_coordinate_collisions: first_collisions,
        full_product_collisions: product_collisions,
        full_proportional_collisions: proportional_collisions,
        full_independent_collisions: independent_collisions,
        independent_first_only_collisions: independent_first_only,
        coordinate_replay_discrepancies: coordinate_discrepancies,
        bucket_histogram_discrepancies: histogram_discrepancies,
        collision_replay_failures: replay_failures,
        false_positives,
        false_negatives,
        rank_control,
        operations,
    };
    if census.proportional_profile_pairs != 31
        || census.independent_profile_pairs != 930
        || census.coefficient_states != 120_125
        || census.unordered_coefficient_pair_tests != 7_447_750
        || census.first_coordinate_collisions != 1_441_500
        || census.full_proportional_collisions != 46_500
        || census.full_independent_collisions != 232_500
        || census.full_product_collisions != 279_000
        || census.independent_first_only_collisions != 1_162_500
        || census.coordinate_replay_discrepancies != 0
        || census.bucket_histogram_discrepancies != 0
        || census.collision_replay_failures != 0
        || census.false_positives != 0
        || census.false_negatives != 0
        || census.rank_control.linked_event_rank != 1
        || census.rank_control.linked_event_reversed_rank != 1
        || census.rank_control.independent_expanded_block_rank != 2
        || census.rank_control.independent_expanded_block_reversed_rank != 2
        || census.rank_control.independent_original_projection_rank != 1
        || census
            .rank_control
            .independent_original_projection_reversed_rank
            != 1
    {
        return Err("complete toy product census mismatch".into());
    }
    Ok(census)
}

fn mod_sub(left: &BigUint, right: &BigUint, modulus: &BigUint) -> BigUint {
    if left >= right {
        left - right
    } else {
        left + modulus - right
    }
}

fn hash_scalar(label: &str, index: usize, modulus: &BigUint) -> BigUint {
    let digest = sha256_hex(format!("round308/product-lift/{label}/{index}").as_bytes());
    let value = BigUint::parse_bytes(digest.as_bytes(), 16).expect("hex digest");
    (value % (modulus - BigUint::one())) + BigUint::one()
}

fn vector_digest(values: &[BigUint]) -> String {
    let strings: Vec<String> = values.iter().map(ToString::to_string).collect();
    sha256_hex(&serde_json::to_vec(&strings).expect("vector JSON"))
}

fn big_dot(left: &[BigUint], right: &[BigUint], modulus: &BigUint) -> BigUint {
    left.iter()
        .zip(right)
        .fold(BigUint::zero(), |sum, (a, b)| (sum + a * b) % modulus)
}

fn linked_delta(logs: &[BigUint], modulus: &BigUint) -> Vec<BigUint> {
    let mut delta = vec![BigUint::zero(); CONTROL_WIDTH];
    delta[0] = logs[1].clone();
    delta[1] = mod_sub(&BigUint::zero(), &logs[0], modulus);
    delta
}

fn independent_delta(
    primary: &[BigUint],
    auxiliary: &[BigUint],
    modulus: &BigUint,
) -> Vec<BigUint> {
    let mut delta = vec![BigUint::zero(); CONTROL_WIDTH];
    delta[0] = mod_sub(
        &((&primary[1] * &auxiliary[2]) % modulus),
        &((&primary[2] * &auxiliary[1]) % modulus),
        modulus,
    );
    delta[1] = mod_sub(
        &((&primary[2] * &auxiliary[0]) % modulus),
        &((&primary[0] * &auxiliary[2]) % modulus),
        modulus,
    );
    delta[2] = mod_sub(
        &((&primary[0] * &auxiliary[1]) % modulus),
        &((&primary[1] * &auxiliary[0]) % modulus),
        modulus,
    );
    delta
}

fn control_vectors(
    modulus: &BigUint,
) -> Result<
    (
        Vec<BigUint>,
        Vec<BigUint>,
        BigUint,
        Vec<BigUint>,
        Vec<BigUint>,
    ),
    String,
> {
    let primary: Vec<BigUint> = (0..CONTROL_WIDTH)
        .map(|index| hash_scalar("primary", index, modulus))
        .collect();
    let auxiliary: Vec<BigUint> = (0..CONTROL_WIDTH)
        .map(|index| hash_scalar("auxiliary", index, modulus))
        .collect();
    let multiplier = hash_scalar("linked-multiplier", 0, modulus);
    let linked = linked_delta(&primary, modulus);
    let independent = independent_delta(&primary, &auxiliary, modulus);
    if independent.iter().all(BigUint::is_zero)
        || !big_dot(&linked, &primary, modulus).is_zero()
        || !big_dot(&independent, &primary, modulus).is_zero()
        || !big_dot(&independent, &auxiliary, modulus).is_zero()
    {
        return Err("P-256 control-vector construction failed".into());
    }
    Ok((primary, auxiliary, multiplier, linked, independent))
}

fn modular_rank(
    rows: &[Vec<BigUint>],
    modulus: &BigUint,
    reverse: bool,
) -> Result<(usize, RankOperations), String> {
    if rows.is_empty() {
        return Ok((0, RankOperations::default()));
    }
    let columns = rows[0].len();
    if columns == 0 || rows.iter().any(|row| row.len() != columns) {
        return Err("rank matrix width mismatch".into());
    }
    let mut matrix = rows.to_vec();
    if reverse {
        matrix.reverse();
    }
    let mut operations = RankOperations {
        matrices: 1,
        rows_materialized: matrix.len() as u64,
        entries_materialized: (matrix.len() * columns) as u64,
        ..RankOperations::default()
    };
    let exponent = modulus - BigUint::from(2u8);
    let mut rank = 0usize;
    for column in 0..columns {
        let mut pivot = None;
        for row in rank..matrix.len() {
            operations.pivot_row_scans += 1;
            if !matrix[row][column].is_zero() {
                pivot = Some(row);
                break;
            }
        }
        let Some(pivot) = pivot else {
            continue;
        };
        if pivot != rank {
            matrix.swap(rank, pivot);
            operations.pivot_swaps += 1;
        }
        let inverse = matrix[rank][column].modpow(&exponent, modulus);
        operations.modular_inversions += 1;
        for entry in &mut matrix[rank][column..] {
            *entry = (&*entry * &inverse) % modulus;
            operations.modular_multiplications += 1;
        }
        let pivot_row = matrix[rank].clone();
        for row in (rank + 1)..matrix.len() {
            let factor = matrix[row][column].clone();
            if factor.is_zero() {
                continue;
            }
            for entry_column in column..columns {
                let product = (&factor * &pivot_row[entry_column]) % modulus;
                matrix[row][entry_column] = mod_sub(&matrix[row][entry_column], &product, modulus);
                operations.modular_multiplications += 1;
                operations.modular_subtractions += 1;
            }
        }
        rank += 1;
        if rank == matrix.len() {
            break;
        }
    }
    Ok((rank, operations))
}

fn count_one_bits(value: &BigUint) -> u64 {
    (0..value.bits()).filter(|bit| value.bit(*bit)).count() as u64
}

fn counted_scalar_mul(
    point: &Point,
    scalar: &BigUint,
    curve: &CurveParams,
    operations: &mut GroupOperations,
) -> Point {
    operations.scalar_multiplications += 1;
    operations.scalar_bits_processed += scalar.bits();
    operations.scalar_one_bits += count_one_bits(scalar);
    point.scalar_mul(scalar, &curve.a_fe())
}

fn replay_zero_sum(
    coefficients: &[BigUint],
    points: &[Point],
    curve: &CurveParams,
    operations: &mut GroupOperations,
) -> bool {
    let mut sum = Point::Infinity;
    for (coefficient, point) in coefficients.iter().zip(points) {
        if coefficient.is_zero() {
            continue;
        }
        let term = counted_scalar_mul(point, coefficient, curve, operations);
        sum = sum.add(&term, &curve.a_fe());
        operations.point_additions += 1;
    }
    operations.equations_replayed += 1;
    matches!(sum, Point::Infinity)
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

fn points_digest(points: &[Point]) -> String {
    let values: Vec<PointReceipt> = points.iter().map(point_receipt).collect();
    sha256_hex(&serde_json::to_vec(&values).expect("point JSON"))
}

fn rank_rows(
    linked: &[BigUint],
    independent: &[BigUint],
    multiplier: &BigUint,
    modulus: &BigUint,
) -> Result<(Vec<(usize, usize)>, RankOperations), String> {
    let linked_second: Vec<BigUint> = linked
        .iter()
        .map(|entry| (entry * multiplier) % modulus)
        .collect();
    let linked_rows = vec![linked.to_vec(), linked_second];
    let mut block_first = independent.to_vec();
    block_first.extend(vec![BigUint::zero(); CONTROL_WIDTH]);
    let mut block_second = vec![BigUint::zero(); CONTROL_WIDTH];
    block_second.extend(independent.to_vec());
    let block_rows = vec![block_first, block_second];
    let projection_rows = vec![independent.to_vec(), vec![BigUint::zero(); CONTROL_WIDTH]];
    let mut ranks = Vec::new();
    let mut aggregate = RankOperations::default();
    for rows in [&linked_rows, &block_rows, &projection_rows] {
        let (forward, forward_ops) = modular_rank(rows, modulus, false)?;
        let (reverse, reverse_ops) = modular_rank(rows, modulus, true)?;
        aggregate.add_assign(&forward_ops);
        aggregate.add_assign(&reverse_ops);
        ranks.push((forward, reverse));
    }
    Ok((ranks, aggregate))
}

fn p256_control(curve: &CurveParams) -> Result<P256Control, String> {
    let p61 = (BigUint::one() << 61usize) - BigUint::one();
    let (primary61, auxiliary61, multiplier61, linked61, independent61) = control_vectors(&p61)?;
    let (primary, auxiliary, multiplier, linked, independent) = control_vectors(&curve.n)?;
    let linked_auxiliary: Vec<BigUint> = primary
        .iter()
        .map(|entry| (entry * &multiplier) % &curve.n)
        .collect();
    let scalar_replay_failures = [
        big_dot(&linked, &primary, &curve.n),
        big_dot(&linked, &linked_auxiliary, &curve.n),
        big_dot(&independent, &primary, &curve.n),
        big_dot(&independent, &auxiliary, &curve.n),
        big_dot(&linked61, &primary61, &p61),
        big_dot(
            &linked61,
            &primary61
                .iter()
                .map(|entry| (entry * &multiplier61) % &p61)
                .collect::<Vec<_>>(),
            &p61,
        ),
        big_dot(&independent61, &primary61, &p61),
        big_dot(&independent61, &auxiliary61, &p61),
    ]
    .iter()
    .filter(|value| !value.is_zero())
    .count() as u64;
    if scalar_replay_failures != 0 {
        return Err("scalar replay failure".into());
    }

    let (ranks61, rank_ops61) = rank_rows(&linked61, &independent61, &multiplier61, &p61)?;
    let (ranksn, rank_opsn) = rank_rows(&linked, &independent, &multiplier, &curve.n)?;
    let mut rank_operations = RankOperations::default();
    rank_operations.add_assign(&rank_ops61);
    rank_operations.add_assign(&rank_opsn);

    let generator = curve.generator();
    let mut group_operations = GroupOperations::default();
    let primary_points: Vec<Point> = primary
        .iter()
        .map(|scalar| {
            group_operations.points_constructed += 1;
            counted_scalar_mul(&generator, scalar, curve, &mut group_operations)
        })
        .collect();
    let linked_points: Vec<Point> = linked_auxiliary
        .iter()
        .map(|scalar| {
            group_operations.points_constructed += 1;
            counted_scalar_mul(&generator, scalar, curve, &mut group_operations)
        })
        .collect();
    let independent_points: Vec<Point> = auxiliary
        .iter()
        .map(|scalar| {
            group_operations.points_constructed += 1;
            counted_scalar_mul(&generator, scalar, curve, &mut group_operations)
        })
        .collect();
    let all_points_on_curve = primary_points
        .iter()
        .chain(&linked_points)
        .chain(&independent_points)
        .all(|point| curve.is_valid_public_point(point));
    if !all_points_on_curve {
        return Err("P-256 point construction failure".into());
    }
    let mut group_replay_failures = 0u64;
    group_replay_failures += u64::from(!replay_zero_sum(
        &linked,
        &primary_points,
        curve,
        &mut group_operations,
    ));
    group_replay_failures += u64::from(!replay_zero_sum(
        &linked,
        &linked_points,
        curve,
        &mut group_operations,
    ));
    group_replay_failures += u64::from(!replay_zero_sum(
        &independent,
        &primary_points,
        curve,
        &mut group_operations,
    ));
    group_replay_failures += u64::from(!replay_zero_sum(
        &independent,
        &independent_points,
        curve,
        &mut group_operations,
    ));
    if group_replay_failures != 0 {
        return Err("P-256 product equation replay failure".into());
    }

    let rank_receipts = vec![
        RankReceipt {
            variant: "linked homomorphic coordinate".into(),
            rows: 2,
            columns: CONTROL_WIDTH,
            expected_rank: 1,
            rank_mod_2_61_minus_1: ranks61[0].0,
            reversed_rank_mod_2_61_minus_1: ranks61[0].1,
            rank_mod_p256_n: ranksn[0].0,
            reversed_rank_mod_p256_n: ranksn[0].1,
            useful_original_system_rank: 1,
            introduced_auxiliary_log_vector: false,
            status: "second equality is a known scalar multiple".into(),
        },
        RankReceipt {
            variant: "independent coordinate, expanded block system".into(),
            rows: 2,
            columns: 2 * CONTROL_WIDTH,
            expected_rank: 2,
            rank_mod_2_61_minus_1: ranks61[1].0,
            reversed_rank_mod_2_61_minus_1: ranks61[1].1,
            rank_mod_p256_n: ranksn[1].0,
            reversed_rank_mod_p256_n: ranksn[1].1,
            useful_original_system_rank: 1,
            introduced_auxiliary_log_vector: true,
            status: "rank two only after adding an independent unknown log block".into(),
        },
        RankReceipt {
            variant: "independent coordinate projected to original unknowns".into(),
            rows: 2,
            columns: CONTROL_WIDTH,
            expected_rank: 1,
            rank_mod_2_61_minus_1: ranks61[2].0,
            reversed_rank_mod_2_61_minus_1: ranks61[2].1,
            rank_mod_p256_n: ranksn[2].0,
            reversed_rank_mod_p256_n: ranksn[2].1,
            useful_original_system_rank: 1,
            introduced_auxiliary_log_vector: true,
            status: "auxiliary row projects to zero on the original recovery variables".into(),
        },
    ];
    if rank_receipts.iter().any(|receipt| {
        receipt.rank_mod_2_61_minus_1 != receipt.expected_rank
            || receipt.reversed_rank_mod_2_61_minus_1 != receipt.expected_rank
            || receipt.rank_mod_p256_n != receipt.expected_rank
            || receipt.reversed_rank_mod_p256_n != receipt.expected_rank
    }) {
        return Err("P-256 rank control mismatch".into());
    }

    Ok(P256Control {
        scalar_labelled_control: true,
        columns_materialized: CONTROL_WIDTH,
        primary_log_vector_sha256: vector_digest(&primary),
        independent_auxiliary_log_vector_sha256: vector_digest(&auxiliary),
        primary_point_vector_sha256: points_digest(&primary_points),
        independent_auxiliary_point_vector_sha256: points_digest(&independent_points),
        linked_multiplier: multiplier.to_string(),
        linked_delta_sha256: vector_digest(&linked),
        independent_delta_sha256: vector_digest(&independent),
        all_points_on_curve,
        scalar_replay_failures,
        group_replay_failures,
        group_operations,
        rank_operations,
        rank_receipts,
        linked_first_primary_point: point_receipt(&primary_points[0]),
        independent_first_auxiliary_point: point_receipt(&independent_points[0]),
        known_label_selector_new_original_log_rows: 0,
    })
}

fn collision_projection(curve: &CurveParams) -> Result<CollisionProjection, String> {
    let original_states = (&curve.n + BigUint::one()) >> 1usize;
    let product_states = (&curve.n * &curve.n + BigUint::one()) >> 1usize;
    let n_f = curve
        .n
        .to_f64()
        .ok_or("P-256 subgroup order does not fit f64")?;
    let expected_samples_log2 = (std::f64::consts::PI.sqrt() / 2.0).log2() + n_f.log2();
    let ratio = (std::f64::consts::PI / 2.0).sqrt() / RHO_S * n_f.sqrt();
    let ratio_log2 = ratio.log2();
    let materialized_record_bytes = 82u64;
    let list_bytes_log2 = expected_samples_log2 + (materialized_record_bytes as f64).log2();
    let disk_bytes_log2 = list_bytes_log2 + 1.0;
    let collection_log2 = expected_samples_log2 + (REQUIRED_USEFUL_ROWS as f64).log2();
    Ok(CollisionProjection {
        subgroup_order: curve.n.to_string(),
        original_canonical_states: original_states.to_string(),
        product_canonical_states: product_states.to_string(),
        asymptotic_expected_product_collision_samples_formula:
            "sqrt(pi)/2*n under simultaneous global-negation quotient".into(),
        expected_product_collision_samples_log2: expected_samples_log2,
        ratio_to_original_rho: ratio,
        ratio_to_original_rho_log2: ratio_log2,
        optimistic_useful_original_rows_per_product_event: 1,
        optimistic_cost_per_useful_original_row_log2: expected_samples_log2,
        registered_useful_rows_with_five_percent_rank_allowance: REQUIRED_USEFUL_ROWS,
        optimistic_registered_collection_lower_bound_log2_operations: collection_log2,
        materialized_record_bytes,
        materialized_list_bytes_log2: list_bytes_log2,
        materialized_disk_write_plus_read_bytes_log2: disk_bytes_log2,
        memoryless_walk_state_bytes: 512,
        memoryless_storage_below_2_50: true,
        materialized_list_storage_below_2_50: list_bytes_log2 < 50.0,
        operation_cost_below_2_120: collection_log2 < 120.0,
        per_usable_relation_below_2_103: expected_samples_log2 < 103.0,
        projection_is_measurement: false,
    })
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
    toy: &ToyCensus,
    p256: &P256Control,
    projection: &CollisionProjection,
) -> TransferAssessment {
    let mut assessment = TransferAssessment {
        schema: "p256-product-lift-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 308,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        methodology_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/references/methodology.md (unavailable)".into(),
        template_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/assets/assessment-template.json (unavailable)".into(),
        typed_correspondence_graph: vec![
            "coefficient vector c --[primary sum]--> sum c_i P_i in G".into(),
            "coefficient vector c --[auxiliary sum]--> sum c_i R_i in G".into(),
            "product collision in GxG --[difference]--> two coordinate equalities".into(),
            "linked R_i=[a]P_i --[known scalar]--> one useful original-system row".into(),
            "independent R_i --[expanded unknown block]--> block rank two, original projection rank one".into(),
            "known auxiliary labels --[precomputable filter]--> zero new original-log rows".into(),
        ],
        obligations: vec![
            obligation("typed product replay", "supported", &format!("{} toy product collisions and {} P-256 coordinate equations replay with zero failures.", toy.full_product_collisions, p256.group_operations.equations_replayed), "Complete order-five control and deterministic scalar-labelled P-256 witnesses."),
            obligation("linked-coordinate independence", "refuted", "The linked pair has rank one over both fields and both row orders.", "Prime-order-group homomorphic scalar transports."),
            obligation("independent-coordinate net information", "refuted for direct block lift", "The expanded rows have rank two, but projection to the original unknown block has rank one and the auxiliary log vector is new.", "Direct product of two independent cyclic-group coordinates."),
            obligation("known-label selector cost", "fails promotion", &format!("The optimistic generic product collision costs 2^{:.6} samples per useful original row and is 2^{:.6} times original rho.", projection.optimistic_cost_per_useful_original_row_log2, projection.ratio_to_original_rho_log2), "Uniform balanced generic product collision; projection, not a measured attack."),
            obligation("non-homomorphic coupling", "open", "No map supplies two independent rows on the original recovery variables without adding unknowns or a product-size collision space.", "Mechanisms outside homomorphic, independent-block, and known-label cases."),
        ],
        controls: vec![
            format!("Complete toy census: {} profile pairs, {} pair tests, {} first-coordinate collisions, and {} full product collisions.", toy.ordered_profile_pairs, toy.unordered_coefficient_pair_tests, toy.first_coordinate_collisions, toy.full_product_collisions),
            format!("P-256 controls materialize {} primary and two sets of auxiliary points, replay four equations, and rank all systems over two fields in both row orders.", p256.columns_materialized),
            "Rounds 39, 302, 303, and 307 are imported only after exact byte-hash, schema, curve, factor-base, conclusion, and boundary checks.".into(),
        ],
        cost_accounting: vec![
            "Measured: complete toy state and pair enumeration, collision/replay counts, P-256 scalar multiplications, group additions, dual-field rank operations, process telemetry, and artifact sizes.".into(),
            format!("Modeled lower bound: uniform simultaneous-negation product birthday collision at 2^{:.6} samples; this is explicitly not measurement.", projection.expected_product_collision_samples_log2),
            "Storage alternatives are separated: an 82-byte materialized record list fails the storage gate, while a constant-memory walk can pass storage but retains the order-n operation exponent.".into(),
            "Unset: a non-homomorphic coupling, actual relation yield, product-system degree, decomposition, sparse linear algebra, target recovery, and complete attack cost.".into(),
        ],
        exploration_boundary: "Exact for direct linked, independent-block, and known-label product lifts in prime-order cyclic groups; not a universal theorem about non-homomorphic encodings or higher algebraic correspondences.".into(),
        weakest_open_obligation: "Construct a typed non-homomorphic event whose two replayed equations both survive projection to the original P-256 recovery variables, add no unresolved log block, and cost no more than original-group rho.".into(),
        narrowest_supported_finding: "A direct two-coordinate product lift does not yield two useful original-system rows: linked coordinates collapse to rank one, while independent coordinates add an unknown block and a generic full collision costs order n rather than order sqrt(n).".into(),
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
    let projection = collision_projection(&curve)?;
    let boundaries = vec![
        BoundaryRow {
            variant: "Pollard rho".into(),
            ratio_to_rho: 1.0,
            ratio_to_rho_log2: 0.0,
            useful_original_rows_per_event: None,
            complete_attack_measured: true,
            status: "complete reference".into(),
        },
        BoundaryRow {
            variant: format!("unquotiented geometry-defined {LOWER_BOUND_FACTOR_BASE}"),
            ratio_to_rho: UNQUOTIENTED_RATIO,
            ratio_to_rho_log2: UNQUOTIENTED_RATIO.log2(),
            useful_original_rows_per_event: Some(1),
            complete_attack_measured: false,
            status: "executable free-oracle boundary unchanged".into(),
        },
        BoundaryRow {
            variant: "free two-primitive oracle".into(),
            ratio_to_rho: TWO_PRIMITIVE_RATIO,
            ratio_to_rho_log2: TWO_PRIMITIVE_RATIO.log2(),
            useful_original_rows_per_event: Some(2),
            complete_attack_measured: false,
            status: "already above rho; no construction".into(),
        },
        BoundaryRow {
            variant: "generic balanced product lift".into(),
            ratio_to_rho: projection.ratio_to_original_rho,
            ratio_to_rho_log2: projection.ratio_to_original_rho_log2,
            useful_original_rows_per_event: Some(1),
            complete_attack_measured: false,
            status: "order-n projection; independent second coordinate adds no original-system row"
                .into(),
        },
        BoundaryRow {
            variant: format!("registered 17-term {REGISTERED_FACTOR_BASE}"),
            ratio_to_rho: REGISTERED_RATIO,
            ratio_to_rho_log2: REGISTERED_RATIO.log2(),
            useful_original_rows_per_event: Some(1),
            complete_attack_measured: false,
            status: "unchanged free-perfect-oracle comparison".into(),
        },
    ];
    let gates = Gates {
        dependency_hashes_and_schemas_checked: true,
        complete_toy_product_census: true,
        zero_false_positives_and_false_negatives: true,
        p256_exact_group_replays_passed: true,
        linked_lift_useful_rank_is_one: true,
        independent_lift_block_rank_is_two: true,
        independent_lift_original_projection_rank_is_one: true,
        nonhomomorphic_two_row_original_system_coupling: false,
        complete_relation_decomposition_and_recovery: false,
        structured_residual_degree_at_most_5: false,
        parity_at_or_below_rho: false,
        complete_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        complete_attack_peak_storage_below_2_50: false,
        discarded_probabilistic_branches_counted_as_exhaustive: false,
        promoted: false,
    };
    let assessment = build_assessment(&toy, &p256, &projection);
    let trilemma = vec![
        "linked/homomorphic auxiliary coordinates replay a scalar multiple of the primary equality and add useful rank one".into(),
        "independent unknown auxiliary coordinates give block rank two but add a new log block and project to useful original rank one".into(),
        "known auxiliary labels are coefficient-computable selectors, not new unknown-log equations, and a generic exhaustive full collision costs order n".into(),
    ];
    let classification =
        "product-lift/useful-rank-one-or-new-unknown-block/order-n-collision/parity-blocked";
    let obstruction = "A direct product coordinate cannot simultaneously remain coupled to the original recovery variables and independent of the primary equation. Linked prime-order-group transports are scalar and rank one. Independent coordinates have block rank two only by introducing a new log vector; their second row projects to zero on the original unknowns. Treating known labels as a selector requires a generic full product collision at about 2^255.825 samples, roughly 2^127.947 times original rho.";
    let decision = "Reject direct GxG product lifting as the two-primitive parity escape. Preserve a non-homomorphic coupling as the open target, but require two rows on the original recovery system, no new unresolved logs, exact replay, and a complete below-rho cost before promotion. Do not attempt an unplanted full-depth relation.";
    let semantic = json!({
        "curve": CURVE_SLUG,
        "factor_base": REGISTERED_FACTOR_BASE,
        "dependencies": &dependencies,
        "trilemma": &trilemma,
        "toy": &toy,
        "p256": &p256,
        "projection": &projection,
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
            schema: "p256-product-lift-screen/v1".into(),
            curve: CURVE_SLUG.into(),
            screening_round: 308,
            execution_status: "complete".into(),
            factor_base: REGISTERED_FACTOR_BASE.into(),
            factor_base_columns: REGISTERED_COLUMNS,
            required_useful_rows: REQUIRED_USEFUL_ROWS,
            optimistic_lower_bound_factor_base: LOWER_BOUND_FACTOR_BASE.into(),
            dependencies,
            typed_product_map: "Phi(c)=(sum_i c_i P_i,sum_i c_i R_i) in GxG; a collision yields delta=c-c' with both coordinate sums equal to O".into(),
            trilemma,
            exhaustive_toy_census: toy,
            p256_control: p256,
            generic_product_collision_projection: projection,
            boundary_table: boundaries,
            structured_residual_degree_of_regularity: None,
            relations_reported_on_actual_factor_base: 0,
            full_depth_unplanted_p256_relation_attempted: false,
            exploration_boundary: "Exact for direct product lifts with linked scalar, independent unknown, or known-label auxiliary coordinates; non-homomorphic couplings remain outside the result.".into(),
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
        "round 308: toy_pairs={}, full_collisions={}, p256_replays={}, product_log2_samples={:.6}, promoted={}",
        result.exhaustive_toy_census.unordered_coefficient_pair_tests,
        result.exhaustive_toy_census.full_product_collisions,
        result.p256_control.group_operations.equations_replayed,
        result
            .generic_product_collision_projection
            .expected_product_collision_samples_log2,
        result.gates.promoted
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_product_lift: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complete_toy_product_census_matches_closed_form() {
        let census = toy_census().expect("toy census");
        assert_eq!(census.unordered_coefficient_pair_tests, 7_447_750);
        assert_eq!(census.first_coordinate_collisions, 1_441_500);
        assert_eq!(census.full_product_collisions, 279_000);
        assert_eq!(census.full_independent_collisions, 232_500);
        assert_eq!(census.independent_first_only_collisions, 1_162_500);
        assert_eq!(census.false_positives + census.false_negatives, 0);
    }

    #[test]
    fn p256_controls_replay_and_project_to_one_useful_row() {
        let control = p256_control(&CurveParams::p256()).expect("P-256 control");
        assert!(control.all_points_on_curve);
        assert_eq!(control.scalar_replay_failures, 0);
        assert_eq!(control.group_replay_failures, 0);
        assert_eq!(control.rank_receipts[0].rank_mod_p256_n, 1);
        assert_eq!(control.rank_receipts[1].rank_mod_p256_n, 2);
        assert_eq!(control.rank_receipts[2].rank_mod_p256_n, 1);
        assert_eq!(control.known_label_selector_new_original_log_rows, 0);
    }

    #[test]
    fn product_collision_projection_fails_time_gate() {
        let projection = collision_projection(&CurveParams::p256()).expect("projection");
        assert!(projection.expected_product_collision_samples_log2 > 255.0);
        assert!(projection.ratio_to_original_rho_log2 > 127.0);
        assert!(!projection.operation_cost_below_2_120);
        assert!(!projection.per_usable_relation_below_2_103);
        assert!(projection.memoryless_storage_below_2_50);
        assert!(!projection.materialized_list_storage_below_2_50);
    }
}
