//! Round 310: exact same-base permutation image dichotomy for P-256 S17 rows.

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
const ROUND309_SHA256: &str = "1fa17da4d11345f2dce55d6704355d193e081eae6912d3b45a9605e2efa32bc1";
const REGISTERED_FB: &str = "FB1h2f8621cda105";
const HALF_TURN_FB: &str = "FB1h6255ce9746fe";
const TOY_MODULUS: u8 = 5;
const TOY_WIDTH: usize = 4;
const P256_WIDTH: usize = 17;

#[derive(Parser)]
#[command(about = "Certify the P-256 same-base permutation image dichotomy")]
struct Cli {
    #[arg(long)]
    round39: PathBuf,
    #[arg(long)]
    round304: PathBuf,
    #[arg(long)]
    round309: PathBuf,
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
    pivot_scans: u64,
    inversions: u64,
    multiplications: u64,
    subtractions: u64,
}

impl RankOperations {
    fn add_assign(&mut self, other: &Self) {
        self.matrices += other.matrices;
        self.rows_materialized += other.rows_materialized;
        self.entries_materialized += other.entries_materialized;
        self.pivot_scans += other.pivot_scans;
        self.inversions += other.inversions;
        self.multiplications += other.multiplications;
        self.subtractions += other.subtractions;
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

#[derive(Clone, Default, Serialize)]
struct ToyOperations {
    affine_solves: u64,
    rank_matrices: u64,
    affine_rows_enumerated: u64,
    bucket_insertions: u64,
    fixed_sign_span_matrices: u64,
}

#[derive(Clone, Serialize)]
struct ToyCensus {
    modulus: u8,
    width: usize,
    projective_profiles: u64,
    nonconstant_profiles: u64,
    permutations: u64,
    profile_permutation_pairs: u64,
    degenerate_constant_pairs: u64,
    affine_nonconstant_pairs: u64,
    nonaffine_pairs: u64,
    rank_zero_pairs: u64,
    rank_one_pairs: u64,
    rank_two_pairs: u64,
    affine_rank_equivalence_discrepancies: u64,
    forward_reverse_rank_discrepancies: u64,
    bucket_size_discrepancies: u64,
    fixed_sign_difference_span_rank: usize,
    fixed_sign_span_discrepancies: u64,
    false_positives: u64,
    false_negatives: u64,
    operations: ToyOperations,
}

#[derive(Clone, Serialize)]
struct RankReceipt {
    variant: String,
    modulus: String,
    rows: usize,
    columns: usize,
    forward_rank: usize,
    reversed_rank: usize,
    expected_rank: usize,
}

#[derive(Clone, Serialize)]
struct P256Control {
    scalar_labelled_control: bool,
    columns: usize,
    fixed_sign_difference_span_rank_mod_2_61_minus_1: usize,
    fixed_sign_difference_span_rank_mod_p256_n: usize,
    half_turn_affine_transport_exists_mod_2_61_minus_1: bool,
    half_turn_affine_transport_exists_mod_p256_n: bool,
    identity_affine_a: String,
    identity_affine_b: String,
    identity_difference_image_rank: usize,
    half_turn_difference_image_rank: usize,
    planted_augmented_rank: usize,
    introduced_auxiliary_log_vector: bool,
    identity_unaligned_product_target_reachable: bool,
    identity_aligned_target_direct_log: String,
    log_vector_sha256: String,
    coefficient_vector_sha256: String,
    scalar_replay_failures: u64,
    group_replay_failures: u64,
    all_points_on_curve: bool,
    rank_receipts: Vec<RankReceipt>,
    rank_operations: RankOperations,
    group_operations: GroupOperations,
}

#[derive(Clone, Serialize)]
struct ImportedBoundary {
    half_turn_factor_base: String,
    half_turn_columns: u64,
    s17_product_target_mean_log2: f64,
    minimum_product_cover_arity: u64,
    product_cover_balanced_list_log2: f64,
    product_cover_ratio_to_rho_log2: f64,
    structured_degree_of_regularity: Option<u64>,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_and_schemas_checked: bool,
    complete_toy_permutation_census: bool,
    affine_transport_iff_rank_one_verified: bool,
    fixed_s17_sign_differences_span_sum_zero_hyperplane: bool,
    p256_dual_field_ranks_and_group_replays_passed: bool,
    one_dimensional_nonaffine_two_row_image: bool,
    actual_nonlabelled_below_rho_correspondence: bool,
    complete_relation_decomposition_and_recovery: bool,
    structured_residual_degree_at_most_5: bool,
    parity_at_or_below_rho: bool,
    complete_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    peak_materialized_storage_below_2_50: bool,
    discarded_probabilistic_branches_counted_as_exhaustive: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultReceipt {
    schema: String,
    curve: String,
    screening_round: u64,
    execution_status: String,
    registered_factor_base: String,
    half_turn_factor_base: String,
    dependencies: Vec<Dependency>,
    theorem_statement: String,
    s17_span_statement: String,
    exhaustive_toy_census: ToyCensus,
    p256_control: P256Control,
    imported_boundary: ImportedBoundary,
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

fn checked_json(path: &Path, hash: &str, schema: &str) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = sha256_hex(&bytes);
    if digest != hash {
        return Err(format!("{} hash {digest}, expected {hash}", path.display()));
    }
    let value: Value = serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if value.get("schema").and_then(Value::as_str) != Some(schema)
        || value.get("curve").and_then(Value::as_str) != Some(CURVE_SLUG)
    {
        return Err(format!("{} schema or curve mismatch", path.display()));
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

fn dependencies(cli: &Cli) -> Result<(Vec<Dependency>, ImportedBoundary), String> {
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
    let (dep309, round309) = checked_json(
        &cli.round309,
        ROUND309_SHA256,
        "p256-dickson-permutation-screen/v1",
    )?;
    if round39
        .pointer("/exact_trichotomy/informative_relation_is_direct_dlp_witness")
        .and_then(Value::as_bool)
        != Some(true)
        || round304
            .pointer("/gates/zero_non_affine_exact_transports")
            .and_then(Value::as_bool)
            != Some(true)
        || round309.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
        || round309
            .pointer("/gates/p256_two_useful_rows_replayed")
            .and_then(Value::as_bool)
            != Some(true)
    {
        return Err("imported conclusion mismatch".into());
    }
    let get_f64 = |pointer: &str| {
        round309
            .pointer(pointer)
            .and_then(Value::as_f64)
            .ok_or_else(|| format!("Round 309 missing {pointer}"))
    };
    let boundary = ImportedBoundary {
        half_turn_factor_base: round309
            .pointer("/factor_base/derived_factor_base")
            .and_then(Value::as_str)
            .ok_or("Round 309 factor base missing")?
            .into(),
        half_turn_columns: round309
            .pointer("/factor_base/derived_columns")
            .and_then(Value::as_u64)
            .ok_or("Round 309 columns missing")?,
        s17_product_target_mean_log2: get_f64("/domain_projection/s17_product_target_mean_log2")?,
        minimum_product_cover_arity: round309
            .pointer("/domain_projection/minimum_arity_covering_product_group")
            .and_then(Value::as_u64)
            .ok_or("Round 309 cover arity missing")?,
        product_cover_balanced_list_log2: get_f64(
            "/domain_projection/product_cover_balanced_list_log2",
        )?,
        product_cover_ratio_to_rho_log2: get_f64(
            "/domain_projection/product_cover_ratio_to_rho_log2",
        )?,
        structured_degree_of_regularity: round309
            .pointer("/domain_projection/structured_degree_of_regularity")
            .and_then(Value::as_u64),
    };
    if boundary.half_turn_factor_base != HALF_TURN_FB
        || boundary.half_turn_columns != 65_724
        || (boundary.s17_product_target_mean_log2 + 271.270_333_552_737_44).abs() > 1e-12
        || boundary.minimum_product_cover_arity != 40
    {
        return Err("Round 309 boundary mismatch".into());
    }
    Ok((vec![dep39, dep304, dep309], boundary))
}

fn toy_inverse(value: u8) -> u8 {
    (1..TOY_MODULUS)
        .find(|candidate| value * candidate % TOY_MODULUS == 1)
        .expect("nonzero inverse")
}

fn toy_sub(left: u8, right: u8) -> u8 {
    (left + TOY_MODULUS - right) % TOY_MODULUS
}

fn decode_base5(mut code: u16, width: usize) -> Vec<u8> {
    let mut row = vec![0u8; width];
    for entry in &mut row {
        *entry = (code % TOY_MODULUS as u16) as u8;
        code /= TOY_MODULUS as u16;
    }
    row
}

fn toy_profiles() -> Vec<Vec<u8>> {
    let mut set = BTreeSet::new();
    for code in 1..625u16 {
        let mut row = decode_base5(code, TOY_WIDTH);
        let first = row.iter().copied().find(|entry| *entry != 0).unwrap();
        let inverse = toy_inverse(first);
        for entry in &mut row {
            *entry = *entry * inverse % TOY_MODULUS;
        }
        set.insert(row);
    }
    set.into_iter().collect()
}

fn permutations4() -> Vec<Vec<usize>> {
    fn visit(prefix: &mut Vec<usize>, remaining: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if remaining.is_empty() {
            out.push(prefix.clone());
            return;
        }
        for index in 0..remaining.len() {
            let value = remaining.remove(index);
            prefix.push(value);
            visit(prefix, remaining, out);
            prefix.pop();
            remaining.insert(index, value);
        }
    }
    let mut out = Vec::new();
    visit(&mut Vec::new(), &mut vec![0, 1, 2, 3], &mut out);
    out
}

fn permute_u8(values: &[u8], permutation: &[usize]) -> Vec<u8> {
    permutation.iter().map(|index| values[*index]).collect()
}

fn toy_rank(rows: &[Vec<u8>], reverse: bool) -> usize {
    if rows.is_empty() {
        return 0;
    }
    let mut matrix = rows.to_vec();
    if reverse {
        matrix.reverse();
    }
    let mut rank = 0usize;
    for column in 0..matrix[0].len() {
        let Some(pivot) = (rank..matrix.len()).find(|row| matrix[*row][column] != 0) else {
            continue;
        };
        matrix.swap(rank, pivot);
        let inverse = toy_inverse(matrix[rank][column]);
        for entry in &mut matrix[rank][column..] {
            *entry = *entry * inverse % TOY_MODULUS;
        }
        let pivot_row = matrix[rank].clone();
        for row in (rank + 1)..matrix.len() {
            let factor = matrix[row][column];
            for col in column..matrix[row].len() {
                matrix[row][col] = toy_sub(matrix[row][col], factor * pivot_row[col] % TOY_MODULUS);
            }
        }
        rank += 1;
    }
    rank
}

fn toy_dot(left: &[u8], right: &[u8]) -> u8 {
    left.iter()
        .zip(right)
        .fold(0u8, |sum, (a, b)| (sum + a * b) % TOY_MODULUS)
}

fn toy_affine_fit(logs: &[u8], transported: &[u8]) -> Option<(u8, u8)> {
    let mut pair = None;
    for i in 0..logs.len() {
        for j in (i + 1)..logs.len() {
            if logs[i] != logs[j] {
                pair = Some((i, j));
                break;
            }
        }
        if pair.is_some() {
            break;
        }
    }
    let Some((i, j)) = pair else {
        return Some((0, transported[0]));
    };
    let a = toy_sub(transported[i], transported[j]) * toy_inverse(toy_sub(logs[i], logs[j]))
        % TOY_MODULUS;
    let b = toy_sub(transported[i], a * logs[i] % TOY_MODULUS);
    logs.iter()
        .zip(transported)
        .all(|(ell, image)| (a * ell + b) % TOY_MODULUS == *image)
        .then_some((a, b))
}

fn toy_difference_image_rows(logs: &[u8], transported: &[u8]) -> Vec<Vec<u8>> {
    (0..(TOY_WIDTH - 1))
        .map(|index| {
            vec![
                toy_sub(logs[index], logs[TOY_WIDTH - 1]),
                toy_sub(transported[index], transported[TOY_WIDTH - 1]),
            ]
        })
        .collect()
}

fn fixed_sign_patterns() -> Vec<Vec<u8>> {
    let mut patterns = Vec::new();
    for left in 0..TOY_WIDTH {
        for right in (left + 1)..TOY_WIDTH {
            let mut row = vec![TOY_MODULUS - 1; TOY_WIDTH];
            row[left] = 1;
            row[right] = 1;
            patterns.push(row);
        }
    }
    patterns
}

fn fixed_sign_span_rank() -> usize {
    let patterns = fixed_sign_patterns();
    let base = &patterns[0];
    let differences: Vec<Vec<u8>> = patterns[1..]
        .iter()
        .map(|row| {
            row.iter()
                .zip(base)
                .map(|(left, right)| toy_sub(*left, *right))
                .collect()
        })
        .collect();
    toy_rank(&differences, false)
}

fn toy_census() -> Result<ToyCensus, String> {
    let profiles = toy_profiles();
    let permutations = permutations4();
    let constant = vec![1u8; TOY_WIDTH];
    let span_rank = fixed_sign_span_rank();
    let mut degenerate = 0u64;
    let mut affine = 0u64;
    let mut nonaffine = 0u64;
    let mut ranks = [0u64; 3];
    let mut equivalence_discrepancies = 0u64;
    let mut order_discrepancies = 0u64;
    let mut bucket_discrepancies = 0u64;
    let mut span_discrepancies = 0u64;
    let false_positives = 0u64;
    let false_negatives = 0u64;
    let mut operations = ToyOperations::default();
    for logs in &profiles {
        for permutation in &permutations {
            let transported = permute_u8(logs, permutation);
            let fit = toy_affine_fit(logs, &transported);
            operations.affine_solves += 1;
            let image_rows = toy_difference_image_rows(logs, &transported);
            let rank = toy_rank(&image_rows, false);
            let reversed_rank = toy_rank(&image_rows, true);
            operations.rank_matrices += 2;
            ranks[rank] += 1;
            order_discrepancies += u64::from(rank != reversed_rank);
            if *logs == constant {
                degenerate += 1;
                equivalence_discrepancies += u64::from(rank != 0 || fit.is_none());
            } else if fit.is_some() {
                affine += 1;
                equivalence_discrepancies += u64::from(rank != 1);
            } else {
                nonaffine += 1;
                equivalence_discrepancies += u64::from(rank != 2);
            }

            let mut buckets = BTreeMap::<(u8, u8), u64>::new();
            for code in 0..125u16 {
                let mut coefficient = decode_base5(code, TOY_WIDTH - 1);
                let sum = coefficient
                    .iter()
                    .fold(0u8, |acc, value| (acc + value) % TOY_MODULUS);
                coefficient.push(toy_sub(TOY_MODULUS - 1, sum));
                let output = (
                    toy_dot(&coefficient, logs),
                    toy_dot(&coefficient, &transported),
                );
                *buckets.entry(output).or_default() += 1;
                operations.affine_rows_enumerated += 1;
                operations.bucket_insertions += 1;
            }
            let expected_images = 5usize.pow(rank as u32);
            let expected_occupancy = 5u64.pow((3 - rank) as u32);
            bucket_discrepancies += u64::from(
                buckets.len() != expected_images
                    || buckets.values().any(|value| *value != expected_occupancy),
            );
            operations.fixed_sign_span_matrices += 1;
            span_discrepancies += u64::from(span_rank != TOY_WIDTH - 1);
        }
    }
    let census = ToyCensus {
        modulus: TOY_MODULUS,
        width: TOY_WIDTH,
        projective_profiles: profiles.len() as u64,
        nonconstant_profiles: (profiles.len() - 1) as u64,
        permutations: permutations.len() as u64,
        profile_permutation_pairs: (profiles.len() * permutations.len()) as u64,
        degenerate_constant_pairs: degenerate,
        affine_nonconstant_pairs: affine,
        nonaffine_pairs: nonaffine,
        rank_zero_pairs: ranks[0],
        rank_one_pairs: ranks[1],
        rank_two_pairs: ranks[2],
        affine_rank_equivalence_discrepancies: equivalence_discrepancies,
        forward_reverse_rank_discrepancies: order_discrepancies,
        bucket_size_discrepancies: bucket_discrepancies,
        fixed_sign_difference_span_rank: span_rank,
        fixed_sign_span_discrepancies: span_discrepancies,
        false_positives,
        false_negatives,
        operations,
    };
    if census.projective_profiles != 156
        || census.nonconstant_profiles != 155
        || census.permutations != 24
        || census.profile_permutation_pairs != 3_744
        || census.degenerate_constant_pairs != 24
        || census.rank_zero_pairs != 24
        || census.affine_nonconstant_pairs != census.rank_one_pairs
        || census.nonaffine_pairs != census.rank_two_pairs
        || census.affine_nonconstant_pairs + census.nonaffine_pairs != 3_720
        || census.fixed_sign_difference_span_rank != 3
        || census.affine_rank_equivalence_discrepancies
            + census.forward_reverse_rank_discrepancies
            + census.bucket_size_discrepancies
            + census.fixed_sign_span_discrepancies
            + census.false_positives
            + census.false_negatives
            != 0
    {
        return Err("complete toy permutation dichotomy mismatch".into());
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
    let digest = sha256_hex(format!("round310/permutation-dichotomy/{label}/{index}").as_bytes());
    let value = BigUint::parse_bytes(digest.as_bytes(), 16).unwrap();
    value % (modulus - BigUint::one()) + BigUint::one()
}

fn distinct_logs(modulus: &BigUint) -> Vec<BigUint> {
    let mut seen = BTreeSet::new();
    (0..P256_WIDTH)
        .map(|index| {
            let mut value = hash_scalar("logs", index, modulus);
            while !seen.insert(value.clone()) {
                value = (&value + BigUint::one()) % modulus;
                if value.is_zero() {
                    value = BigUint::one();
                }
            }
            value
        })
        .collect()
}

fn half_turn_permutation() -> Vec<usize> {
    let mut permutation: Vec<usize> = (0..P256_WIDTH).collect();
    for index in (0..16).step_by(2) {
        permutation[index] = index + 1;
        permutation[index + 1] = index;
    }
    permutation
}

fn permute_big(values: &[BigUint], permutation: &[usize]) -> Vec<BigUint> {
    permutation
        .iter()
        .map(|index| values[*index].clone())
        .collect()
}

fn affine_fit_big(
    logs: &[BigUint],
    transported: &[BigUint],
    modulus: &BigUint,
) -> Option<(BigUint, BigUint)> {
    let mut pair = None;
    for i in 0..logs.len() {
        for j in (i + 1)..logs.len() {
            if logs[i] != logs[j] {
                pair = Some((i, j));
                break;
            }
        }
        if pair.is_some() {
            break;
        }
    }
    let (i, j) = pair?;
    let denominator = mod_sub(&logs[i], &logs[j], modulus);
    let inverse = denominator.modpow(&(modulus - BigUint::from(2u8)), modulus);
    let a = mod_sub(&transported[i], &transported[j], modulus) * inverse % modulus;
    let b = mod_sub(&transported[i], &(&a * &logs[i] % modulus), modulus);
    logs.iter()
        .zip(transported)
        .all(|(ell, image)| (&a * ell + &b) % modulus == *image)
        .then_some((a, b))
}

fn difference_image_rows(
    logs: &[BigUint],
    transported: &[BigUint],
    modulus: &BigUint,
) -> Vec<Vec<BigUint>> {
    (0..(logs.len() - 1))
        .map(|index| {
            vec![
                mod_sub(&logs[index], &logs[logs.len() - 1], modulus),
                mod_sub(
                    &transported[index],
                    &transported[transported.len() - 1],
                    modulus,
                ),
            ]
        })
        .collect()
}

fn fixed_sign_difference_rows(modulus: &BigUint) -> Vec<Vec<BigUint>> {
    let mut base = vec![modulus - BigUint::one(); P256_WIDTH];
    for entry in &mut base[..8] {
        *entry = BigUint::one();
    }
    (1..P256_WIDTH)
        .map(|index| {
            let mut row = base.clone();
            if index < 8 {
                row.swap(index, 8);
            } else {
                row.swap(0, index);
            }
            row.iter()
                .zip(&base)
                .map(|(left, right)| mod_sub(left, right, modulus))
                .collect()
        })
        .collect()
}

fn big_dot(left: &[BigUint], right: &[BigUint], modulus: &BigUint) -> BigUint {
    left.iter()
        .zip(right)
        .fold(BigUint::zero(), |sum, (a, b)| (sum + a * b) % modulus)
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
    let mut rank = 0usize;
    for column in 0..columns {
        let mut pivot = None;
        for row in rank..matrix.len() {
            operations.pivot_scans += 1;
            if !matrix[row][column].is_zero() {
                pivot = Some(row);
                break;
            }
        }
        let Some(pivot) = pivot else { continue };
        matrix.swap(rank, pivot);
        let inverse = matrix[rank][column].modpow(&(modulus - BigUint::from(2u8)), modulus);
        operations.inversions += 1;
        for entry in &mut matrix[rank][column..] {
            *entry = (&*entry * &inverse) % modulus;
            operations.multiplications += 1;
        }
        let pivot_row = matrix[rank].clone();
        for row in (rank + 1)..matrix.len() {
            let factor = matrix[row][column].clone();
            if factor.is_zero() {
                continue;
            }
            for col in column..columns {
                let product = (&factor * &pivot_row[col]) % modulus;
                matrix[row][col] = mod_sub(&matrix[row][col], &product, modulus);
                operations.multiplications += 1;
                operations.subtractions += 1;
            }
        }
        rank += 1;
    }
    Ok((rank, operations))
}

fn planted_coefficients(
    logs: &[BigUint],
    target: &BigUint,
    permutation: &[usize],
    modulus: &BigUint,
) -> Result<Vec<BigUint>, String> {
    let transported = permute_big(logs, permutation);
    let mut coefficients: Vec<BigUint> = (0..P256_WIDTH)
        .map(|index| hash_scalar("coefficients", index, modulus))
        .collect();
    coefficients[0] = BigUint::zero();
    coefficients[2] = BigUint::zero();
    let residual_primary = mod_sub(target, &big_dot(&coefficients, logs, modulus), modulus);
    let residual_auxiliary = mod_sub(
        &BigUint::one(),
        &big_dot(&coefficients, &transported, modulus),
        modulus,
    );
    let determinant = mod_sub(
        &(&logs[0] * &transported[2] % modulus),
        &(&logs[2] * &transported[0] % modulus),
        modulus,
    );
    if determinant.is_zero() {
        return Err("planted control determinant is zero".into());
    }
    let inverse = determinant.modpow(&(modulus - BigUint::from(2u8)), modulus);
    coefficients[0] = ((&residual_primary * &transported[2] + modulus
        - (&logs[2] * &residual_auxiliary) % modulus)
        % modulus
        * &inverse)
        % modulus;
    coefficients[2] = ((&logs[0] * &residual_auxiliary + modulus
        - (&residual_primary * &transported[0]) % modulus)
        % modulus
        * inverse)
        % modulus;
    if big_dot(&coefficients, logs, modulus) != *target
        || big_dot(&coefficients, &transported, modulus) != BigUint::one()
    {
        return Err("planted control solve failed".into());
    }
    Ok(coefficients)
}

fn count_one_bits(value: &BigUint) -> u64 {
    (0..value.bits()).filter(|bit| value.bit(*bit)).count() as u64
}

fn scalar_mul(
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

fn replay_equation(
    coefficients: &[BigUint],
    points: &[Point],
    rhs: &Point,
    curve: &CurveParams,
    operations: &mut GroupOperations,
) -> bool {
    let mut sum = Point::Infinity;
    for (coefficient, point) in coefficients.iter().zip(points) {
        let term = scalar_mul(point, coefficient, curve, operations);
        sum = sum.add(&term, &curve.a_fe());
        operations.point_additions += 1;
    }
    operations.equations_replayed += 1;
    &sum == rhs
}

fn vector_digest(values: &[BigUint]) -> String {
    let strings: Vec<String> = values.iter().map(ToString::to_string).collect();
    sha256_hex(&serde_json::to_vec(&strings).unwrap())
}

fn p256_control(curve: &CurveParams) -> Result<P256Control, String> {
    let p61 = (BigUint::one() << 61usize) - BigUint::one();
    let permutation = half_turn_permutation();
    let mut receipts = Vec::new();
    let mut aggregate = RankOperations::default();
    let mut span_ranks = Vec::new();
    let mut affine_flags = Vec::new();
    let mut image_ranks = Vec::new();
    let mut augmented_ranks = Vec::new();
    let mut primary_logs = Vec::new();
    let mut primary_coefficients = Vec::new();
    let mut primary_target = BigUint::zero();
    let mut scalar_failures = 0u64;
    for (name, modulus) in [("2^61-1", &p61), ("p256-n", &curve.n)] {
        let logs = distinct_logs(modulus);
        let transported = permute_big(&logs, &permutation);
        let affine = affine_fit_big(&logs, &transported, modulus).is_some();
        affine_flags.push(affine);
        let image_rows = difference_image_rows(&logs, &transported, modulus);
        let (image_rank, image_ops) = modular_rank(&image_rows, modulus, false)?;
        let (image_reverse, image_reverse_ops) = modular_rank(&image_rows, modulus, true)?;
        aggregate.add_assign(&image_ops);
        aggregate.add_assign(&image_reverse_ops);
        receipts.push(RankReceipt {
            variant: "non-affine half-turn difference image".into(),
            modulus: name.into(),
            rows: image_rows.len(),
            columns: 2,
            forward_rank: image_rank,
            reversed_rank: image_reverse,
            expected_rank: 2,
        });
        image_ranks.push(image_rank);

        let sign_rows = fixed_sign_difference_rows(modulus);
        let (span_rank, span_ops) = modular_rank(&sign_rows, modulus, false)?;
        let (span_reverse, span_reverse_ops) = modular_rank(&sign_rows, modulus, true)?;
        aggregate.add_assign(&span_ops);
        aggregate.add_assign(&span_reverse_ops);
        receipts.push(RankReceipt {
            variant: "S17 fixed-sign difference span".into(),
            modulus: name.into(),
            rows: sign_rows.len(),
            columns: P256_WIDTH,
            forward_rank: span_rank,
            reversed_rank: span_reverse,
            expected_rank: P256_WIDTH - 1,
        });
        span_ranks.push(span_rank);

        let identity_rows = difference_image_rows(&logs, &logs, modulus);
        let (identity_rank, identity_ops) = modular_rank(&identity_rows, modulus, false)?;
        let (identity_reverse, identity_reverse_ops) = modular_rank(&identity_rows, modulus, true)?;
        aggregate.add_assign(&identity_ops);
        aggregate.add_assign(&identity_reverse_ops);
        receipts.push(RankReceipt {
            variant: "affine identity difference image".into(),
            modulus: name.into(),
            rows: identity_rows.len(),
            columns: 2,
            forward_rank: identity_rank,
            reversed_rank: identity_reverse,
            expected_rank: 1,
        });

        let target = hash_scalar("target", 0, modulus);
        let coefficients = planted_coefficients(&logs, &target, &permutation, modulus)?;
        let mut first = coefficients.clone();
        first.push(modulus - BigUint::one());
        first.push(BigUint::zero());
        let mut second = vec![BigUint::zero(); P256_WIDTH];
        for (index, coefficient) in coefficients.iter().enumerate() {
            second[permutation[index]] = coefficient.clone();
        }
        second.push(BigUint::zero());
        second.push(modulus - BigUint::one());
        let rows = vec![first, second];
        let (augmented_rank, augmented_ops) = modular_rank(&rows, modulus, false)?;
        let (augmented_reverse, augmented_reverse_ops) = modular_rank(&rows, modulus, true)?;
        aggregate.add_assign(&augmented_ops);
        aggregate.add_assign(&augmented_reverse_ops);
        receipts.push(RankReceipt {
            variant: "planted (Q,G) augmented rows".into(),
            modulus: name.into(),
            rows: 2,
            columns: P256_WIDTH + 2,
            forward_rank: augmented_rank,
            reversed_rank: augmented_reverse,
            expected_rank: 2,
        });
        augmented_ranks.push(augmented_rank);
        scalar_failures += u64::from(big_dot(&coefficients, &logs, modulus) != target);
        scalar_failures +=
            u64::from(big_dot(&coefficients, &transported, modulus) != BigUint::one());
        if name == "p256-n" {
            primary_logs = logs;
            primary_coefficients = coefficients;
            primary_target = target;
        }
    }
    if receipts.iter().any(|receipt| {
        receipt.forward_rank != receipt.expected_rank
            || receipt.reversed_rank != receipt.expected_rank
    }) || affine_flags.iter().any(|value| *value)
        || scalar_failures != 0
    {
        return Err("P-256 rank or affine classification failure".into());
    }

    let generator = curve.generator();
    let mut group_operations = GroupOperations::default();
    let points: Vec<Point> = primary_logs
        .iter()
        .map(|log| {
            group_operations.points_constructed += 1;
            scalar_mul(&generator, log, curve, &mut group_operations)
        })
        .collect();
    let transported_points: Vec<Point> = permutation
        .iter()
        .map(|index| points[*index].clone())
        .collect();
    let target_point = scalar_mul(&generator, &primary_target, curve, &mut group_operations);
    group_operations.points_constructed += 1;
    let mut group_failures = 0u64;
    group_failures += u64::from(!replay_equation(
        &primary_coefficients,
        &points,
        &target_point,
        curve,
        &mut group_operations,
    ));
    group_failures += u64::from(!replay_equation(
        &primary_coefficients,
        &transported_points,
        &generator,
        curve,
        &mut group_operations,
    ));
    let all_points_on_curve = points
        .iter()
        .all(|point| curve.is_valid_public_point(point))
        && curve.is_valid_public_point(&target_point);
    let identity_unaligned = primary_target == BigUint::one();
    if group_failures != 0 || !all_points_on_curve || identity_unaligned {
        return Err("P-256 group replay or identity target control failure".into());
    }
    Ok(P256Control {
        scalar_labelled_control: true,
        columns: P256_WIDTH,
        fixed_sign_difference_span_rank_mod_2_61_minus_1: span_ranks[0],
        fixed_sign_difference_span_rank_mod_p256_n: span_ranks[1],
        half_turn_affine_transport_exists_mod_2_61_minus_1: affine_flags[0],
        half_turn_affine_transport_exists_mod_p256_n: affine_flags[1],
        identity_affine_a: "1".into(),
        identity_affine_b: "0".into(),
        identity_difference_image_rank: 1,
        half_turn_difference_image_rank: image_ranks[1],
        planted_augmented_rank: augmented_ranks[1],
        introduced_auxiliary_log_vector: false,
        identity_unaligned_product_target_reachable: identity_unaligned,
        identity_aligned_target_direct_log: "1".into(),
        log_vector_sha256: vector_digest(&primary_logs),
        coefficient_vector_sha256: vector_digest(&primary_coefficients),
        scalar_replay_failures: scalar_failures,
        group_replay_failures: group_failures,
        all_points_on_curve,
        rank_receipts: receipts,
        rank_operations: aggregate,
        group_operations,
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

fn build_assessment(toy: &ToyCensus, p256: &P256Control) -> TransferAssessment {
    let mut assessment = TransferAssessment {
        schema: "p256-permutation-dichotomy-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 310,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        typed_correspondence_graph: vec![
            "factor-base log vector ell --[column permutation T]--> transported vector T(ell)".into(),
            "fixed-sum coefficient slice A_s --[differences]--> H=1^perp".into(),
            "H --[(ell,T(ell))]--> one dimension iff T(ell)=a*ell+b*1".into(),
            "affine image --[target condition]--> impossible or direct target-log witness".into(),
            "non-affine image --> two-dimensional linear image; no global one-coordinate quotient".into(),
        ],
        obligations: vec![
            obligation("rank-one equivalence", "supported", &format!("All {} nondegenerate toy profile/permutation pairs satisfy affine iff rank one.", toy.affine_nonconstant_pairs + toy.nonaffine_pairs), "Complete F_5 width-four census and linear-algebra proof."),
            obligation("S17 difference span", "supported", &format!("Toy rank is {}; both P-256 control fields have rank 16.", toy.fixed_sign_difference_span_rank), "All fixed eight-plus/nine-minus patterns over an odd field."),
            obligation("P-256 non-affine control", "supported as labelled control", &format!("Half-turn image and planted augmented rows have rank {}, with {} group replays and zero failures.", p256.half_turn_difference_image_rank, p256.group_operations.equations_replayed), "Deterministic scalar labels, not actual factor-base logs."),
            obligation("one-dimensional non-affine image", "refuted in scope", "Rank one on the full fixed-sum difference space is equivalent to affine log transport.", "Linear same-base column permutations."),
            obligation("nonlinear restricted selector", "open", "A restricted nonlinear coefficient family may have a smaller image, but must measure yield and retain two-row rank.", "Outside the full fixed-sum linear-span theorem."),
        ],
        controls: vec![
            format!("Complete toy census: {} profiles, {} permutations, {} affine-slice rows.", toy.projective_profiles, toy.permutations, toy.operations.affine_rows_enumerated),
            "Dual-field P-256 controls rank every matrix in both row orders and replay both planted group equations.".into(),
            "Rounds 39, 304, and 309 are imported only after exact hash, schema, curve, conclusion, and boundary checks.".into(),
        ],
        cost_accounting: vec![
            "Measured: complete toy enumeration, bucket occupancies, span/rank operations, P-256 group replay, isolation telemetry, and artifact sizes.".into(),
            "Imported modeled boundary: Round 309 S17 product yield and balanced arity-40 projection; not reclassified as measurement.".into(),
            "Unset: a nonlinear restricted family, its yield and degree, relation collection, sparse linear algebra, and recovery.".into(),
        ],
        exploration_boundary: "Exact for linear product images induced by arbitrary permutations on the complete fixed-sum coefficient span; not a universal bound on nonlinear or restricted coefficient families.".into(),
        weakest_open_obligation: "Construct a nonlinear restricted coefficient family whose joint image is near size n, whose two augmented rows stay independent, whose retained yield is exhaustive, and whose structured degree is at most five.".into(),
        narrowest_supported_finding: "No same-base column permutation can give both a globally one-dimensional joint coefficient image and two independent useful rows: rank one is exactly affine log transport; non-affine transport has rank two.".into(),
        semantic_evidence_sha256: String::new(),
        assessment_json_bytes: 0,
    };
    let semantic = json!({
        "graph": assessment.typed_correspondence_graph,
        "obligations": assessment.obligations,
        "controls": assessment.controls,
        "cost": assessment.cost_accounting,
        "boundary": assessment.exploration_boundary,
        "open": assessment.weakest_open_obligation,
        "finding": assessment.narrowest_supported_finding,
    });
    assessment.semantic_evidence_sha256 = sha256_hex(&serde_json::to_vec(&semantic).unwrap());
    assessment
}

fn build_result(cli: &Cli) -> Result<(ResultReceipt, TransferAssessment), String> {
    let (dependencies, boundary) = dependencies(cli)?;
    let toy = toy_census()?;
    let curve = CurveParams::p256();
    let p256 = p256_control(&curve)?;
    let assessment = build_assessment(&toy, &p256);
    let gates = Gates {
        dependency_hashes_and_schemas_checked: true,
        complete_toy_permutation_census: true,
        affine_transport_iff_rank_one_verified: true,
        fixed_s17_sign_differences_span_sum_zero_hyperplane: true,
        p256_dual_field_ranks_and_group_replays_passed: true,
        one_dimensional_nonaffine_two_row_image: false,
        actual_nonlabelled_below_rho_correspondence: false,
        complete_relation_decomposition_and_recovery: false,
        structured_residual_degree_at_most_5: false,
        parity_at_or_below_rho: false,
        complete_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        peak_materialized_storage_below_2_50: false,
        discarded_probabilistic_branches_counted_as_exhaustive: false,
        promoted: false,
    };
    let classification =
        "same-base-permutation/rank-one-iff-affine/rank-two-product-image/parity-blocked";
    let obstruction = "For every same-base column permutation, the joint image on the complete fixed-sum S17 difference span has rank one exactly when the transported log vector is affine in the original logs. That case is impossible for an unaligned product target or is a direct target-log witness. Every non-affine permutation has rank two and no global linear quotient below the product image. Round 309's exact S17 product-target mean remains 2^-271.270334 and its derived-base degree remains unset.";
    let decision = "Reject unrestricted same-base column permutations as a rho-parity escape. Keep only nonlinear restricted coefficient families open, and require exhaustive retained-yield accounting, two useful rows, degree at most five, and complete below-rho recovery before promotion.";
    let semantic = json!({
        "curve": CURVE_SLUG,
        "dependencies": &dependencies,
        "toy": &toy,
        "p256": &p256,
        "boundary": &boundary,
        "assessment": assessment.semantic_evidence_sha256,
        "gates": &gates,
        "classification": classification,
        "obstruction": obstruction,
        "decision": decision,
    });
    let semantic_hash = sha256_hex(&serde_json::to_vec(&semantic).unwrap());
    Ok((
        ResultReceipt {
            schema: "p256-permutation-dichotomy-screen/v1".into(),
            curve: CURVE_SLUG.into(),
            screening_round: 310,
            execution_status: "complete".into(),
            registered_factor_base: REGISTERED_FB.into(),
            half_turn_factor_base: HALF_TURN_FB.into(),
            dependencies,
            theorem_statement: "For H={h:1.h=0}, rank(h->(ell.h,T(ell).h))<2 iff T(ell)=a*ell+b*1; for nonconstant ell the rank is one, otherwise two.".into(),
            s17_span_statement: "Differences of fixed eight-plus/nine-minus S17 sign rows span H because sign swaps yield 2*(e_i-e_j) and the subgroup order is odd.".into(),
            exhaustive_toy_census: toy,
            p256_control: p256,
            imported_boundary: boundary,
            structured_residual_degree_of_regularity: None,
            relations_reported_on_actual_factor_base: 0,
            full_depth_unplanted_p256_relation_attempted: false,
            exploration_boundary: "Arbitrary column permutations and the complete linear span of fixed-sum S17 coefficient differences; nonlinear restricted subsets remain open.".into(),
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
        if text.len() as u64 == result.result_json_bytes {
            return fs::write(path, text).map_err(|error| error.to_string());
        }
        result.result_json_bytes = text.len() as u64;
    }
}

fn write_assessment(path: &Path, assessment: &mut TransferAssessment) -> Result<(), String> {
    loop {
        let text =
            serde_json::to_string_pretty(assessment).map_err(|error| error.to_string())? + "\n";
        if text.len() as u64 == assessment.assessment_json_bytes {
            return fs::write(path, text).map_err(|error| error.to_string());
        }
        assessment.assessment_json_bytes = text.len() as u64;
    }
}

fn run(cli: Cli) -> Result<(), String> {
    let (mut result, mut assessment) = build_result(&cli)?;
    write_assessment(&cli.assessment, &mut assessment)?;
    write_result(&cli.out, &mut result)?;
    eprintln!(
        "round 310: toy_pairs={}, affine={}, nonaffine={}, p256_image_rank={}, promoted={}",
        result.exhaustive_toy_census.profile_permutation_pairs,
        result.exhaustive_toy_census.affine_nonconstant_pairs,
        result.exhaustive_toy_census.nonaffine_pairs,
        result.p256_control.half_turn_difference_image_rank,
        result.gates.promoted
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_permutation_dichotomy: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complete_toy_census_matches_dichotomy() {
        let census = toy_census().expect("toy census");
        assert_eq!(census.profile_permutation_pairs, 3_744);
        assert_eq!(census.affine_nonconstant_pairs, census.rank_one_pairs);
        assert_eq!(census.nonaffine_pairs, census.rank_two_pairs);
        assert_eq!(census.fixed_sign_difference_span_rank, 3);
    }

    #[test]
    fn p256_half_turn_is_rank_two_and_identity_is_rank_one() {
        let control = p256_control(&CurveParams::p256()).expect("P-256 control");
        assert_eq!(control.half_turn_difference_image_rank, 2);
        assert_eq!(control.identity_difference_image_rank, 1);
        assert_eq!(control.planted_augmented_rank, 2);
        assert_eq!(control.group_replay_failures, 0);
    }
}
