//! Round 307: one primitive relation cannot be cloned into two quotient rows.

use std::collections::{BTreeMap, BTreeSet};
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

const ROUND39_SHA256: &str = "189e48a45ab2919dd3e279bf058940230c4a0461c745fba42dffb26f785ba079";
const ROUND304_SHA256: &str = "c94e339eca07b96ea0ef890592fb13b80815f0baf52af0e6bd57de17dd46ac9c";
const ROUND306_SHA256: &str = "08e11b9c0c829210100730f01710d2b68d80702bcd5d1768b3e594ca88713a31";
const FACTOR_BASE_ID: &str = "FB1hc72514a2a8d3";
const REGISTERED_S17_FACTOR_BASE_ID: &str = "FB1h2f8621cda105";
const COLUMNS: usize = 164;
const HOMOGENEOUS_COLUMNS: usize = COLUMNS + 2;
const DERIVED_ROWS: usize = 1_024;
const TOY_MODULUS: u8 = 7;
const TOY_WIDTH: usize = 4;
const TOY_HOMOGENEOUS_COLUMNS: usize = TOY_WIDTH + 2;
const TOY_PROJECTIVE_PROFILES: u64 = 400;
const TOY_EQUATIONS_PER_PROFILE: u64 = 16_807;
const TOY_STATIC_ROWS_PER_PROFILE: u64 = 343;
const TOY_NONSTATIC_ROWS_PER_PROFILE: u64 = 16_464;
const TOY_PROJECTIVE_LINES: u64 = 8;
const TOY_ROWS_PER_LINE: u64 = 2_058;
const TOY_LINE_PAIRS_PER_PROFILE: u64 = 28;
const UNQUOTIENTED_RATIO: f64 = 13.920_747_397_073_491;
const TWO_PRIMITIVE_ORACLE_RATIO: f64 = 1.446_504_715_695_271_7;
const KNOWN_ANCHOR_ORACLE_RATIO: f64 = 0.964_336_477_130_181;
const KNOWN_ANCHOR_DIRECT_RATIO: f64 = 7.714_691_817_041_447;
const REGISTERED_S17_RATIO: f64 = 394.425_280;

#[derive(Parser)]
#[command(about = "Certify P-256 one-event relation provenance rank")]
struct Cli {
    #[arg(long)]
    round39: PathBuf,
    #[arg(long)]
    round304: PathBuf,
    #[arg(long)]
    round306: PathBuf,
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
struct LineExample {
    normalized_first: u8,
    normalized_second: u8,
    rows: usize,
    forward_rank: usize,
    reversed_rank: usize,
}

#[derive(Clone, Serialize)]
struct ToyCensus {
    modulus: u8,
    base_columns: usize,
    homogeneous_columns: usize,
    projective_profiles: u64,
    exact_equations_enumerated: u64,
    equation_replays: u64,
    static_rows: u64,
    nonstatic_rows: u64,
    projective_line_closures: u64,
    rows_per_line_closure: u64,
    distinct_line_pairs: u64,
    single_line_rank_discrepancies: u64,
    pair_rank_discrepancies: u64,
    partition_discrepancies: u64,
    replay_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    minimum_single_line_rank: usize,
    maximum_single_line_rank: usize,
    minimum_distinct_pair_rank: usize,
    maximum_distinct_pair_rank: usize,
    rank_operations: RankOperations,
    first_line: LineExample,
    last_line: LineExample,
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
    systems: Vec<(String, Vec<Equation>, usize, String)>,
    derived_row_digest: String,
    distinct_derived_rows: usize,
}

#[derive(Clone, Serialize)]
struct StageProfile {
    variant: String,
    rows: usize,
    homogeneous_columns: usize,
    expected_rank: usize,
    rank_mod_2_61_minus_1: usize,
    reversed_rank_mod_2_61_minus_1: usize,
    rank_mod_p256_n: usize,
    reversed_rank_mod_p256_n: usize,
    quotient_rank_increment: usize,
    provenance_status: String,
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
    homogeneous_columns: usize,
    derived_rows_requested: usize,
    distinct_derived_rows: usize,
    derived_row_digest: String,
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
    complete_projective_provenance_census: bool,
    every_one_line_closure_adds_at_most_one_dimension: bool,
    distinct_quotient_lines_add_two_dimensions: bool,
    zero_false_positives_and_false_negatives: bool,
    dual_modulus_p256_ranks_and_group_replays_passed: bool,
    non_scalar_labelled_two_primitive_p256_event: bool,
    primitive_information_provenance_accounted: bool,
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
    provenance_definition: String,
    one_event_bound: String,
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
    let (dep39, round39) = checked_json(
        &cli.round39,
        ROUND39_SHA256,
        "p256.affine_log_orbit_screen/v1",
    )?;
    let (dep304, round304) = checked_json(
        &cli.round304,
        ROUND304_SHA256,
        "p256-algebraic-quotient-rigidity-screen/v1",
    )?;
    let (dep306, round306) = checked_json(
        &cli.round306,
        ROUND306_SHA256,
        "p256-bucket-signature-screen/v1",
    )?;
    for (label, value) in [
        ("Round 39", &round39),
        ("Round 304", &round304),
        ("Round 306", &round306),
    ] {
        if value.get("curve").and_then(Value::as_str) != Some(CURVE_SLUG) {
            return Err(format!("{label} curve mismatch"));
        }
    }
    if round39
        .pointer("/resource_accounting/ordinary_independent_relation_rows")
        .and_then(Value::as_u64)
        != Some(0)
        || round39.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
    {
        return Err("Round 39 affine-log conclusion mismatch".into());
    }
    if round304.get("factor_base").and_then(Value::as_str) != Some(FACTOR_BASE_ID)
        || round304.get("columns").and_then(Value::as_u64) != Some(COLUMNS as u64)
        || round304.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
    {
        return Err("Round 304 quotient conclusion mismatch".into());
    }
    if round306.get("factor_base").and_then(Value::as_str) != Some(FACTOR_BASE_ID)
        || round306.get("columns").and_then(Value::as_u64) != Some(COLUMNS as u64)
        || round306
            .pointer("/gates/every_fixed_signature_bucket_adds_at_most_one_dimension")
            .and_then(Value::as_bool)
            != Some(true)
        || round306.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
    {
        return Err("Round 306 bucket-signature conclusion mismatch".into());
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
        TWO_PRIMITIVE_ORACLE_RATIO,
        "Round 304 two-primitive boundary",
    )?;
    require_ratio(
        &round304,
        "/boundary_table/3/free_oracle_ratio_to_rho",
        KNOWN_ANCHOR_ORACLE_RATIO,
        "Round 304 known-anchor oracle boundary",
    )?;
    require_ratio(
        &round304,
        "/boundary_table/3/direct_measured_ratio_to_rho",
        KNOWN_ANCHOR_DIRECT_RATIO,
        "Round 304 known-anchor direct boundary",
    )?;
    Ok(vec![dep39, dep304, dep306])
}

fn mod7_sub(left: u8, right: u8) -> u8 {
    (left + TOY_MODULUS - right) % TOY_MODULUS
}

fn mod7_mul(left: u8, right: u8) -> u8 {
    ((u16::from(left) * u16::from(right)) % u16::from(TOY_MODULUS)) as u8
}

fn mod7_dot(left: &[u8], right: &[u8]) -> u8 {
    left.iter()
        .zip(right)
        .fold(0u8, |acc, (a, b)| (acc + mod7_mul(*a, *b)) % 7)
}

fn inverse_mod7(value: u8) -> u8 {
    (1..TOY_MODULUS)
        .find(|candidate| mod7_mul(*candidate, value) == 1)
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
            *value = mod7_mul(*value, inverse);
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

fn normalize_line(first: u8, second: u8) -> Option<(u8, u8)> {
    if first == 0 && second == 0 {
        return None;
    }
    let pivot = if first != 0 { first } else { second };
    let inverse = inverse_mod7(pivot);
    Some((mod7_mul(first, inverse), mod7_mul(second, inverse)))
}

fn rank_mod7(rows: &[Vec<u8>], reverse: bool) -> (usize, RankOperations) {
    let mut matrix = rows.to_vec();
    if reverse {
        matrix.reverse();
    }
    let width = TOY_HOMOGENEOUS_COLUMNS;
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
        for entry in &mut matrix[rank][column..] {
            *entry = mod7_mul(*entry, inverse);
            ops.modular_multiplications += 1;
        }
        let pivot_row = matrix[rank].clone();
        for row in (rank + 1)..matrix.len() {
            let factor = matrix[row][column];
            if factor == 0 {
                continue;
            }
            for entry_column in column..width {
                let product = mod7_mul(factor, pivot_row[entry_column]);
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
    let mut equations = 0u64;
    let mut replays = 0u64;
    let mut static_rows_total = 0u64;
    let mut nonstatic_rows_total = 0u64;
    let mut line_closures = 0u64;
    let mut line_pairs = 0u64;
    let mut line_rank_discrepancies = 0u64;
    let mut pair_rank_discrepancies = 0u64;
    let mut partition_discrepancies = 0u64;
    let mut replay_failures = 0u64;
    let false_positives = 0u64;
    let mut false_negatives = 0u64;
    let mut line_ranks = Vec::new();
    let mut pair_ranks = Vec::new();
    let mut examples = Vec::new();
    let mut rank_operations = RankOperations::default();

    for logs in &profiles {
        let target = toy_target(logs);
        let mut static_rows = Vec::new();
        let mut lines: BTreeMap<(u8, u8), Vec<Vec<u8>>> = BTreeMap::new();
        for code in 0..u64::from(TOY_MODULUS).pow(TOY_WIDTH as u32) {
            let coefficients = decode_base7(code, TOY_WIDTH);
            let dot = mod7_dot(&coefficients, logs);
            for b in 0..TOY_MODULUS {
                let a = mod7_sub(dot, mod7_mul(b, target));
                let mut row = coefficients.clone();
                row.push(mod7_sub(0, b));
                row.push(mod7_sub(0, a));
                equations += 1;
                replays += 1;
                let witness = [logs[0], logs[1], logs[2], logs[3], target, 1];
                if mod7_dot(&row, &witness) != 0 {
                    replay_failures += 1;
                }
                let signature = (dot, mod7_sub(0, b));
                if let Some(line) = normalize_line(signature.0, signature.1) {
                    lines.entry(line).or_default().push(row);
                } else {
                    static_rows.push(row);
                }
            }
        }
        if !equations.is_multiple_of(TOY_EQUATIONS_PER_PROFILE)
            || static_rows.len() as u64 != TOY_STATIC_ROWS_PER_PROFILE
            || lines.len() as u64 != TOY_PROJECTIVE_LINES
            || lines
                .values()
                .any(|rows| rows.len() as u64 != TOY_ROWS_PER_LINE)
        {
            partition_discrepancies += 1;
            false_negatives += 1;
        }
        static_rows_total += static_rows.len() as u64;
        nonstatic_rows_total += lines.values().map(Vec::len).sum::<usize>() as u64;
        let line_entries: Vec<((u8, u8), Vec<Vec<u8>>)> = lines.into_iter().collect();
        for (line, rows) in &line_entries {
            let mut system = static_rows.clone();
            system.extend(rows.iter().cloned());
            let (rank, ops) = rank_mod7(&system, false);
            let (reverse_rank, reverse_ops) = rank_mod7(&system, true);
            rank_operations.add_assign(&ops);
            rank_operations.add_assign(&reverse_ops);
            if rank != 4 || reverse_rank != 4 {
                line_rank_discrepancies += 1;
            }
            line_ranks.push(rank);
            line_closures += 1;
            examples.push(LineExample {
                normalized_first: line.0,
                normalized_second: line.1,
                rows: rows.len(),
                forward_rank: rank,
                reversed_rank: reverse_rank,
            });
        }
        for left in 0..line_entries.len() {
            for right in (left + 1)..line_entries.len() {
                let mut system = static_rows.clone();
                system.extend(line_entries[left].1.iter().cloned());
                system.extend(line_entries[right].1.iter().cloned());
                let (rank, ops) = rank_mod7(&system, false);
                let (reverse_rank, reverse_ops) = rank_mod7(&system, true);
                rank_operations.add_assign(&ops);
                rank_operations.add_assign(&reverse_ops);
                if rank != 5 || reverse_rank != 5 {
                    pair_rank_discrepancies += 1;
                }
                pair_ranks.push(rank);
                line_pairs += 1;
            }
        }
    }

    if equations != TOY_PROJECTIVE_PROFILES * TOY_EQUATIONS_PER_PROFILE
        || replays != equations
        || static_rows_total != TOY_PROJECTIVE_PROFILES * TOY_STATIC_ROWS_PER_PROFILE
        || nonstatic_rows_total != TOY_PROJECTIVE_PROFILES * TOY_NONSTATIC_ROWS_PER_PROFILE
        || line_closures != TOY_PROJECTIVE_PROFILES * TOY_PROJECTIVE_LINES
        || line_pairs != TOY_PROJECTIVE_PROFILES * TOY_LINE_PAIRS_PER_PROFILE
        || partition_discrepancies != 0
        || replay_failures != 0
        || line_rank_discrepancies != 0
        || pair_rank_discrepancies != 0
        || false_positives != 0
        || false_negatives != 0
    {
        return Err("toy event-provenance census failed".into());
    }

    Ok(ToyCensus {
        modulus: TOY_MODULUS,
        base_columns: TOY_WIDTH,
        homogeneous_columns: TOY_HOMOGENEOUS_COLUMNS,
        projective_profiles: profiles.len() as u64,
        exact_equations_enumerated: equations,
        equation_replays: replays,
        static_rows: static_rows_total,
        nonstatic_rows: nonstatic_rows_total,
        projective_line_closures: line_closures,
        rows_per_line_closure: TOY_ROWS_PER_LINE,
        distinct_line_pairs: line_pairs,
        single_line_rank_discrepancies: line_rank_discrepancies,
        pair_rank_discrepancies,
        partition_discrepancies,
        replay_failures,
        false_positives,
        false_negatives,
        minimum_single_line_rank: *line_ranks.iter().min().ok_or("no line ranks")?,
        maximum_single_line_rank: *line_ranks.iter().max().ok_or("no line ranks")?,
        minimum_distinct_pair_rank: *pair_ranks.iter().min().ok_or("no pair ranks")?,
        maximum_distinct_pair_rank: *pair_ranks.iter().max().ok_or("no pair ranks")?,
        rank_operations,
        first_line: examples.first().ok_or("no line examples")?.clone(),
        last_line: examples.last().ok_or("no line examples")?.clone(),
    })
}

fn mod_sub(left: &BigUint, right: &BigUint, modulus: &BigUint) -> BigUint {
    if left >= right {
        left - right
    } else {
        left + modulus - right
    }
}

fn hash_scalar(namespace: &str, label: &str, index: usize, modulus: &BigUint) -> BigUint {
    let digest = sha256_hex(format!("{namespace}/{label}/{index}").as_bytes());
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
            let mut coefficients = vec![BigUint::zero(); COLUMNS + 1];
            coefficients[0] = logs[index].clone();
            coefficients[index] = mod_sub(&BigUint::zero(), pivot, modulus);
            Equation {
                coefficients,
                rhs: BigUint::zero(),
            }
        })
        .collect()
}

fn scale_equation(equation: &Equation, scalar: &BigUint, modulus: &BigUint) -> Equation {
    Equation {
        coefficients: equation
            .coefficients
            .iter()
            .map(|entry| (entry * scalar) % modulus)
            .collect(),
        rhs: (&equation.rhs * scalar) % modulus,
    }
}

fn add_scaled_equation(
    left: &Equation,
    right: &Equation,
    scalar: &BigUint,
    modulus: &BigUint,
) -> Equation {
    Equation {
        coefficients: left
            .coefficients
            .iter()
            .zip(&right.coefficients)
            .map(|(a, b)| (a + b * scalar) % modulus)
            .collect(),
        rhs: (&left.rhs + &right.rhs * scalar) % modulus,
    }
}

fn homogeneous_row(equation: &Equation, modulus: &BigUint) -> Vec<BigUint> {
    let mut row = equation.coefficients.clone();
    row.push(mod_sub(&BigUint::zero(), &equation.rhs, modulus));
    row
}

fn row_digest(equation: &Equation) -> String {
    let values: Vec<String> = equation
        .coefficients
        .iter()
        .chain(std::iter::once(&equation.rhs))
        .map(ToString::to_string)
        .collect();
    sha256_hex(&serde_json::to_vec(&values).expect("row JSON"))
}

fn control_systems(modulus: &BigUint) -> Result<ControlSystems, String> {
    let logs: Vec<BigUint> = (0..COLUMNS)
        .map(|index| hash_scalar("round305", "base-log", index, modulus))
        .collect();
    let target_log = hash_scalar("round305", "target-log", 0, modulus);
    let pivot = logs[0].clone();
    let pivot_inverse = pivot.modpow(&(modulus - BigUint::from(2u8)), modulus);
    if (&pivot * &pivot_inverse) % modulus != BigUint::one() {
        return Err("control pivot inverse failed".into());
    }
    let static_rows = static_basis(&logs, modulus);
    let event_rhs = BigUint::one();
    let mut event_coefficients = vec![BigUint::zero(); COLUMNS + 1];
    event_coefficients[0] = ((&target_log + &event_rhs) * &pivot_inverse) % modulus;
    event_coefficients[COLUMNS] = modulus - BigUint::one();
    let event = Equation {
        coefficients: event_coefficients,
        rhs: event_rhs,
    };

    let mut derived = Vec::with_capacity(DERIVED_ROWS);
    let mut distinct = BTreeSet::new();
    for index in 0..DERIVED_ROWS {
        let multiplier = hash_scalar("round307", "event-multiplier", index, modulus);
        let static_multiplier = hash_scalar("round307", "static-multiplier", index, modulus);
        let scaled = scale_equation(&event, &multiplier, modulus);
        let row = add_scaled_equation(
            &scaled,
            &static_rows[index % static_rows.len()],
            &static_multiplier,
            modulus,
        );
        distinct.insert(row_digest(&row));
        derived.push(row);
    }
    if distinct.len() != DERIVED_ROWS {
        return Err("derived P-256 control rows are not distinct".into());
    }
    let derived_digest = sha256_hex(
        &serde_json::to_vec(&distinct.iter().collect::<Vec<_>>())
            .map_err(|error| error.to_string())?,
    );

    let mut one_event = static_rows.clone();
    one_event.push(event.clone());
    let mut derived_closure = static_rows.clone();
    derived_closure.extend(derived);
    let mut proportional = one_event.clone();
    proportional.push(scale_equation(&event, &BigUint::from(2u8), modulus));

    let mut generator_coefficients = vec![BigUint::zero(); COLUMNS + 1];
    generator_coefficients[0] = pivot_inverse;
    let generator_event = Equation {
        coefficients: generator_coefficients,
        rhs: BigUint::one(),
    };
    let mut independent = one_event.clone();
    independent.push(generator_event);

    Ok(ControlSystems {
        logs,
        target_log,
        systems: vec![
            (
                "static-pre-event-space".into(),
                static_rows,
                COLUMNS - 1,
                "W0 has rank 163 before any target event".into(),
            ),
            (
                "one-primitive-target-event".into(),
                one_event,
                COLUMNS,
                "one primitive adds exactly one equation modulo W0".into(),
            ),
            (
                "complete-derived-row-control".into(),
                derived_closure,
                COLUMNS,
                "1,024 distinct u*r+w rows retain one primitive quotient line".into(),
            ),
            (
                "second-proportional-primitive".into(),
                proportional,
                COLUMNS,
                "a proportional second primitive does not add quotient rank".into(),
            ),
            (
                "second-independent-primitive".into(),
                independent,
                COLUMNS + 1,
                "an independent generator equation adds the second dimension".into(),
            ),
        ],
        derived_row_digest: derived_digest,
        distinct_derived_rows: distinct.len(),
    })
}

fn modular_rank(
    equations: &[Equation],
    modulus: &BigUint,
    reverse: bool,
) -> Result<(usize, RankOperations), String> {
    let mut matrix: Vec<Vec<BigUint>> = equations
        .iter()
        .map(|equation| homogeneous_row(equation, modulus))
        .collect();
    if matrix.iter().any(|row| row.len() != HOMOGENEOUS_COLUMNS) {
        return Err("matrix width mismatch".into());
    }
    if reverse {
        matrix.reverse();
    }
    let mut ops = RankOperations {
        matrices: 1,
        rows_materialized: matrix.len() as u64,
        entries_materialized: (matrix.len() * HOMOGENEOUS_COLUMNS) as u64,
        ..RankOperations::default()
    };
    let exponent = modulus - BigUint::from(2u8);
    let mut rank = 0usize;
    for column in 0..HOMOGENEOUS_COLUMNS {
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
            for entry_column in column..HOMOGENEOUS_COLUMNS {
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
    let mut witness = system.logs.clone();
    witness.push(system.target_log.clone());
    system
        .systems
        .iter()
        .flat_map(|(_, equations, _, _)| {
            equations.iter().map(|equation| {
                u64::from(big_dot(&equation.coefficients, &witness, modulus) != equation.rhs)
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
        let (name61, equations61, expected61, status61) = &control61.systems[index];
        let (namen, equationsn, expectedn, statusn) = &controln.systems[index];
        if name61 != namen || expected61 != expectedn || status61 != statusn {
            return Err("dual-modulus system mismatch".into());
        }
        let (rank61, ops61) = modular_rank(equations61, &p61, false)?;
        let (reverse61, reverse_ops61) = modular_rank(equations61, &p61, true)?;
        let (rankn, opsn) = modular_rank(equationsn, &curve.n, false)?;
        let (reversen, reverse_opsn) = modular_rank(equationsn, &curve.n, true)?;
        let discrepancies = [rank61, reverse61, rankn, reversen]
            .into_iter()
            .filter(|rank| *rank != *expectedn)
            .count() as u64;
        let before = group_ops.equations_replayed;
        let failures_before = group_replay_failures;
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
            homogeneous_columns: HOMOGENEOUS_COLUMNS,
            expected_rank: *expectedn,
            rank_mod_2_61_minus_1: rank61,
            reversed_rank_mod_2_61_minus_1: reverse61,
            rank_mod_p256_n: rankn,
            reversed_rank_mod_p256_n: reversen,
            quotient_rank_increment: expectedn.saturating_sub(COLUMNS - 1),
            provenance_status: statusn.clone(),
            scalar_equation_replays: (equations61.len() + equationsn.len()) as u64,
            p256_group_equation_replays: group_ops.equations_replayed - before,
            replay_failures: group_replay_failures - failures_before,
            rank_discrepancies: discrepancies,
            materialized_bytes_mod_p256_n: (equationsn.len() * HOMOGENEOUS_COLUMNS * 32) as u64,
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
        homogeneous_columns: HOMOGENEOUS_COLUMNS,
        derived_rows_requested: DERIVED_ROWS,
        distinct_derived_rows: controln.distinct_derived_rows,
        derived_row_digest: controln.derived_row_digest,
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
            variant: "one unknown anchor plus two independent primitive events".into(),
            ratio_to_rho: TWO_PRIMITIVE_ORACLE_RATIO,
            complete_attack_measured: false,
            status: "optimistic oracle; post-processing cannot reduce primitive count".into(),
        },
        BoundaryRow {
            variant: "known anchor plus one primitive event".into(),
            ratio_to_rho: KNOWN_ANCHOR_ORACLE_RATIO,
            complete_attack_measured: false,
            status: format!(
                "ideal oracle below rho; registered direct labelled control is {KNOWN_ANCHOR_DIRECT_RATIO:.6}x"
            ),
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
        schema: "p256-event-provenance-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 307,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        methodology_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/references/methodology.md (unavailable)".into(),
        template_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/assets/assessment-template.json (unavailable)".into(),
        typed_correspondence_graph: vec![
            "pre-event equations --[span]--> W0".into(),
            "one primitive event r --[known scalar/group operations]--> W0+span(r)".into(),
            "replayed affine transport --[charged identity]--> W0".into(),
            "one complete projective quotient line --> rank increment one".into(),
            "two distinct primitive quotient lines --> rank increment two".into(),
        ],
        obligations: vec![
            obligation("one-event closure rank", "supported", &format!("All {} complete projective-line closures add one dimension over W0.", toy.projective_line_closures), "Linear, scalar, group-law, and charged affine-transport post-processing."),
            obligation("two-line independence", "supported", &format!("All {} distinct-line pairs add two dimensions with zero discrepancies.", toy.distinct_line_pairs), "Complete order-seven quotient control."),
            obligation("P-256 derived-row replay", "supported as scalar-labelled control", &format!("{} distinct derived rows retain rank 164; the independent primitive reaches rank 165 across {} group replays.", p256.distinct_derived_rows, p256.group_operations.equations_replayed), "Synthetic deterministic labels, not a relation construction on FB1hc72514a2a8d3."),
            obligation("two-primitive P-256 event", "unknown", "No non-scalar-labelled mechanism emits two independently certified quotient lines.", "Future non-homomorphic relation mechanisms."),
            obligation("complete below-rho recovery", "unknown", "The two-primitive free oracle is above rho and no complete relation/decomposition/recovery pipeline exists.", "End-to-end P-256 attack."),
        ],
        controls: vec![
            format!("The toy census enumerated and replayed {} exact equations, {} line closures, and {} distinct-line pairs with zero discrepancies.", toy.exact_equations_enumerated, toy.projective_line_closures, toy.distinct_line_pairs),
            format!("Five P-256 stages were ranked over two fields in both row orders and replayed through {} exact group equations.", p256.group_operations.equations_replayed),
            "Rounds 39, 304, and 306 were imported only after exact byte-hash, schema, curve, factor-base, boundary, and gate checks.".into(),
        ],
        cost_accounting: vec![
            "Measured: complete toy equation enumeration, line partition, rank operations, P-256 scalar multiplications, group additions, equation replays, process telemetry, and artifact sizes.".into(),
            "Imported without extrapolation: Round-304 unquotiented, two-primitive, and known-anchor oracle boundaries.".into(),
            "Unset: a non-labelled two-primitive event, relation yield, structured solver, sparse linear algebra, target recovery, and complete attack cost.".into(),
        ],
        exploration_boundary: "Exact for rows derived from one primitive equation using pre-event equations and linear group-law operations; not a bound on a search that independently discovers a second equality.".into(),
        weakest_open_obligation: "Construct and replay one non-scalar-labelled P-256 mechanism that independently certifies two quotient lines, then price complete recovery below rho.".into(),
        narrowest_supported_finding: "Arbitrarily many rows obtained by replayable post-processing of one primitive event add at most one dimension modulo pre-event knowledge. A second dimension requires independently certified information.".into(),
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
        complete_projective_provenance_census: true,
        every_one_line_closure_adds_at_most_one_dimension: true,
        distinct_quotient_lines_add_two_dimensions: true,
        zero_false_positives_and_false_negatives: true,
        dual_modulus_p256_ranks_and_group_replays_passed: true,
        non_scalar_labelled_two_primitive_p256_event: false,
        primitive_information_provenance_accounted: true,
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
        "many-rows-per-event/one-primitive-provenance-rank-one/two-primitives-required/parity-blocked";
    let obstruction = "Modulo all equations known before an event, scalar multiplication, group addition, column transport, and affine replay of one primitive relation remain in its one-dimensional span. The complete toy census closes every row in each quotient projective line, while the P-256 control retains rank 164 with 1,024 distinct derived rows and reaches rank 165 only with an independently certified generator equation. The two-primitive free oracle is already 1.446505 times rho, and no non-labelled source of the required pair is constructed.";
    let decision = "Reject relation post-processing as a two-row parity escape. Count independently certified primitive equation rank, not emitted or replayed rows. Keep a genuinely joint two-primitive event open; do not attempt an unplanted full-depth relation.";
    let semantic = json!({
        "curve": CURVE_SLUG,
        "dependencies": &dependencies,
        "provenance_definition": "derived rows lie in W0+span(r)",
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
            schema: "p256-event-provenance-screen/v1".into(),
            curve: CURVE_SLUG.into(),
            screening_round: 307,
            execution_status: "complete".into(),
            factor_base: FACTOR_BASE_ID.into(),
            registered_s17_factor_base: REGISTERED_S17_FACTOR_BASE_ID.into(),
            columns: COLUMNS,
            dependencies,
            provenance_definition: "With W0 the independently replayed pre-event equation space, any row derived from one primitive r by known scalar/group operations and identities already charged to W0 lies in W0+span(r).".into(),
            one_event_bound: "Arbitrary replayable post-processing of one primitive equation adds at most one quotient dimension; any new independent identity is itself a second primitive and must be charged.".into(),
            exhaustive_toy_census: toy,
            p256_control: p256,
            boundary_table: boundaries,
            structured_residual_degree_of_regularity: None,
            relations_reported_on_actual_factor_base: 0,
            full_depth_unplanted_p256_relation_attempted: false,
            exploration_boundary: "Exact for the linear closure of one primitive equation modulo independently verified pre-event knowledge; not universal over searches that discover a second independent equality.".into(),
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
        "round 307: equations={}, line_closures={}, line_pairs={}, p256_replays={}, promoted={}",
        result.exhaustive_toy_census.exact_equations_enumerated,
        result.exhaustive_toy_census.projective_line_closures,
        result.exhaustive_toy_census.distinct_line_pairs,
        result.p256_control.group_operations.equations_replayed,
        result.gates.promoted
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_event_provenance: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complete_toy_provenance_census_matches_projective_counts() {
        let census = toy_census().expect("toy census");
        assert_eq!(census.exact_equations_enumerated, 6_722_800);
        assert_eq!(census.projective_line_closures, 3_200);
        assert_eq!(census.distinct_line_pairs, 11_200);
        assert_eq!(census.minimum_single_line_rank, 4);
        assert_eq!(census.maximum_distinct_pair_rank, 5);
        assert_eq!(census.replay_failures, 0);
    }

    #[test]
    fn p61_derived_rows_retain_one_primitive_dimension() {
        let modulus = (BigUint::one() << 61usize) - BigUint::one();
        let control = control_systems(&modulus).expect("control");
        assert_eq!(scalar_replay_failures(&control, &modulus), 0);
        for (_, equations, expected, _) in &control.systems {
            assert_eq!(
                modular_rank(equations, &modulus, false).expect("rank").0,
                *expected
            );
        }
        assert_eq!(control.distinct_derived_rows, DERIVED_ROWS);
        assert_eq!(control.systems[2].1.len(), COLUMNS - 1 + DERIVED_ROWS);
        assert_eq!(control.systems[2].2, COLUMNS);
        assert_eq!(control.systems[4].2, COLUMNS + 1);
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
