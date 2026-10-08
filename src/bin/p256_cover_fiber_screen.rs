//! Round 300: quotient relation multiplicity created by cover fibres.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use blake3::Hasher;
use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::ecc::CurveParams;
use serde::Serialize;
use serde_json::{json, Value};

const ROUND299_SHA256: &str = "e6b78d729e785b1831a061cd39e6e57f10c46da8b8e98418321dc72f56752f35";
const ROUND295_FB_SHA256: &str = "d27516ca40a612ecf3ebabfa8ae04776084c7948110e20f221da438e1e64d8f8";
const FACTOR_BASE_ID: &str = "FB1hc72514a2a8d3";
const REGISTERED_S17_FACTOR_BASE_ID: &str = "FB1h2f8621cda105";
const P: u64 = 1151;
const A: u64 = P - 3;
const B: u64 = 241;
const TOY_GROUP_ORDER: usize = 1192;
const SUBGROUP_ORDER: u64 = 149;
const COFACTOR: u64 = 8;
const COLUMNS: usize = 6;
const ARITY: usize = 3;
const MAP_SCALARS: [u64; 4] = [1, 2, 4, 8];
const ROUND295_RATIO_TO_RHO: f64 = 13.920_747_397_073_491;
const REGISTERED_S17_RATIO_TO_RHO: f64 = 394.425_280;

#[derive(Parser)]
#[command(about = "Measure cover-fibre row rank before and after pushforward quotienting")]
struct Cli {
    #[arg(long)]
    round299: PathBuf,
    #[arg(long)]
    factor_base: PathBuf,
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
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct Affine {
    x: u64,
    y: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
enum ToyPoint {
    Infinity,
    Affine(Affine),
}

#[derive(Clone, Copy, Default, Serialize)]
struct RankOperations {
    additions_or_subtractions: u64,
    multiplications: u64,
    inversions: u64,
}

impl RankOperations {
    fn add_assign(&mut self, other: Self) {
        self.additions_or_subtractions += other.additions_or_subtractions;
        self.multiplications += other.multiplications;
        self.inversions += other.inversions;
    }
}

#[derive(Clone, Copy)]
struct RankMeasurement {
    rank: usize,
    operations: RankOperations,
}

#[derive(Clone, Serialize)]
struct HistogramRow {
    value: u64,
    count: u64,
}

#[derive(Clone, Serialize)]
struct BaseSpec {
    family: String,
    selection_uses_discrete_log_labels: bool,
    folded_scalar_labels: Vec<u64>,
    points: Vec<String>,
}

#[derive(Serialize)]
struct CellResult {
    family: String,
    selection_uses_discrete_log_labels: bool,
    map_scalar: u64,
    algebraic_degree: u64,
    rational_kernel_size: u64,
    rational_kernel_points: Vec<String>,
    fibre_size_histogram_for_six_destinations: Vec<HistogramRow>,
    destination_columns: u64,
    lifted_columns: u64,
    arity: u64,
    destination_raw_rows: u64,
    expected_lifted_rows: u64,
    enumerated_lifted_rows: u64,
    unique_lifted_rows: u64,
    lift_multiplier_per_destination_row: u64,
    nonempty_destination_targets: u64,
    maximum_destination_row_occupancy: u64,
    maximum_lifted_row_occupancy: u64,
    maximum_distinct_source_sums_per_destination_target: u64,
    deck_kernel_rank: u64,
    maximum_raw_lifted_coefficient_rank: u64,
    maximum_raw_lifted_homogeneous_rank: u64,
    maximum_union_with_kernel_rank: u64,
    maximum_homogeneous_union_with_kernel_rank: u64,
    maximum_quotient_coefficient_rank: u64,
    maximum_quotient_homogeneous_rank: u64,
    maximum_collapsed_destination_coefficient_rank: u64,
    maximum_collapsed_destination_homogeneous_rank: u64,
    maximum_raw_rank_removed_by_quotient: u64,
    buckets_with_quotient_identity_failure: u64,
    buckets_with_cover_created_destination_rank: u64,
    missing_fibres: u64,
    duplicate_lift_rows: u64,
    replay_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    fibre_scan_pushforward_doublings: u64,
    destination_sum_group_additions: u64,
    source_sum_group_additions: u64,
    row_pushforward_doublings: u64,
    rank_field_operations_mod_149: RankOperations,
    peak_logical_coefficient_bytes: u64,
    disk_bytes_written_during_enumeration: u64,
    coefficient_rank_histogram: Vec<HistogramRow>,
    quotient_rank_histogram: Vec<HistogramRow>,
    enumeration_blake3: String,
}

#[derive(Serialize)]
struct BoundaryRow {
    variant: String,
    effective_log_classes: Option<u64>,
    cover_sheet_credit: u64,
    optimistic_ratio_to_rho: f64,
    scope: String,
    result: String,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_checked: bool,
    toy_curve_and_prime_subgroup_exact: bool,
    all_sixteen_cells_complete: bool,
    every_lifted_row_group_replayed: bool,
    zero_false_positives_and_false_negatives: bool,
    zero_missing_fibres_and_duplicates: bool,
    quotient_identity_holds_everywhere: bool,
    cover_creates_two_destination_equations_from_one_destination_row: bool,
    explicit_p256_cover_and_recovery_implemented: bool,
    structured_residual_degree_at_most_5: bool,
    parity_at_or_below_rho: bool,
    per_usable_relation_below_2_103: bool,
    complete_collection_below_2_120: bool,
    projected_storage_below_2_50: bool,
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
    dependencies: Vec<Dependency>,
    theorem: String,
    toy_curve: String,
    toy_group_order: u64,
    cofactor: u64,
    prime_subgroup_order: u64,
    subgroup_generator: String,
    ambient_points_enumerated: u64,
    factor_base_controls: Vec<BaseSpec>,
    cells: Vec<CellResult>,
    total_destination_rows: u64,
    total_lifted_rows_replayed: u64,
    total_fibre_scan_pushforward_doublings: u64,
    total_destination_sum_group_additions: u64,
    total_source_sum_group_additions: u64,
    total_row_pushforward_doublings: u64,
    total_rank_field_operations_mod_149: RankOperations,
    maximum_raw_lift_multiplier: u64,
    maximum_observed_raw_rank_gain: u64,
    maximum_observed_quotient_rank_gain: u64,
    p256_prime_order_rational_fibre_statement: String,
    boundary_table: Vec<BoundaryRow>,
    relations_reported_on_p256: u64,
    full_depth_unplanted_p256_relation_attempted: bool,
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
struct TransferPath {
    name: String,
    source: String,
    destination: String,
    transformation_type: String,
    field_of_definition: String,
    source_dimension: Option<u64>,
    destination_dimension: Option<u64>,
    relevant_subgroup: String,
    obligations: Vec<Obligation>,
    conclusion: String,
    exploration_boundary: String,
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
    paths: Vec<TransferPath>,
    controls: Vec<String>,
    cost_accounting: Vec<String>,
    weakest_open_obligation: String,
    narrowest_supported_finding: String,
    semantic_evidence_sha256: String,
    assessment_json_bytes: u64,
}

fn checked_dependency(path: &Path, expected: &str) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = sha256_hex(&bytes);
    if digest != expected {
        return Err(format!(
            "{} SHA-256 is {digest}, expected {expected}",
            path.display()
        ));
    }
    let value =
        serde_json::from_slice(&bytes).map_err(|error| format!("{}: {error}", path.display()))?;
    Ok((
        Dependency {
            path: path.display().to_string(),
            bytes: bytes.len() as u64,
            sha256: digest,
        },
        value,
    ))
}

fn addm(left: u64, right: u64) -> u64 {
    ((left as u128 + right as u128) % P as u128) as u64
}

fn negm(value: u64) -> u64 {
    if value == 0 {
        0
    } else {
        P - value
    }
}

fn subm(left: u64, right: u64) -> u64 {
    addm(left, negm(right))
}

fn mulm(left: u64, right: u64) -> u64 {
    ((left as u128 * right as u128) % P as u128) as u64
}

fn powm(mut base: u64, mut exponent: u64) -> u64 {
    let mut out = 1u64;
    while exponent != 0 {
        if exponent & 1 == 1 {
            out = mulm(out, base);
        }
        base = mulm(base, base);
        exponent >>= 1;
    }
    out
}

fn invm(value: u64) -> u64 {
    assert_ne!(value, 0);
    powm(value, P - 2)
}

fn toy_rhs(x: u64) -> u64 {
    addm(addm(mulm(mulm(x, x), x), mulm(A, x)), B)
}

fn square_roots() -> Vec<Vec<u64>> {
    let mut table = vec![Vec::new(); P as usize];
    for value in 0..P {
        table[mulm(value, value) as usize].push(value);
    }
    table
}

fn negate(point: ToyPoint) -> ToyPoint {
    match point {
        ToyPoint::Infinity => ToyPoint::Infinity,
        ToyPoint::Affine(point) => ToyPoint::Affine(Affine {
            x: point.x,
            y: negm(point.y),
        }),
    }
}

fn add_points(left: ToyPoint, right: ToyPoint) -> ToyPoint {
    match (left, right) {
        (ToyPoint::Infinity, point) | (point, ToyPoint::Infinity) => point,
        (ToyPoint::Affine(left), ToyPoint::Affine(right)) => {
            if left.x == right.x && addm(left.y, right.y) == 0 {
                return ToyPoint::Infinity;
            }
            let slope = if left == right {
                if left.y == 0 {
                    return ToyPoint::Infinity;
                }
                mulm(
                    addm(mulm(3, mulm(left.x, left.x)), A),
                    invm(mulm(2, left.y)),
                )
            } else {
                mulm(subm(right.y, left.y), invm(subm(right.x, left.x)))
            };
            let x = subm(subm(mulm(slope, slope), left.x), right.x);
            let y = subm(mulm(slope, subm(left.x, x)), left.y);
            ToyPoint::Affine(Affine { x, y })
        }
    }
}

fn scalar_mul(mut point: ToyPoint, mut scalar: u64) -> ToyPoint {
    let mut out = ToyPoint::Infinity;
    while scalar != 0 {
        if scalar & 1 == 1 {
            out = add_points(out, point);
        }
        point = add_points(point, point);
        scalar >>= 1;
    }
    out
}

fn power_two_map(mut point: ToyPoint, scalar: u64) -> ToyPoint {
    debug_assert!(scalar.is_power_of_two());
    for _ in 0..scalar.trailing_zeros() {
        point = add_points(point, point);
    }
    point
}

fn all_curve_points(roots: &[Vec<u64>]) -> Vec<ToyPoint> {
    let mut points = vec![ToyPoint::Infinity];
    for x in 0..P {
        for &y in &roots[toy_rhs(x) as usize] {
            points.push(ToyPoint::Affine(Affine { x, y }));
        }
    }
    points.sort_unstable();
    points
}

fn point_label(point: ToyPoint) -> String {
    match point {
        ToyPoint::Infinity => "infinity".into(),
        ToyPoint::Affine(point) => format!("{},{}", point.x, point.y),
    }
}

struct ToyGroup {
    points: Vec<ToyPoint>,
    subgroup_generator: ToyPoint,
    subgroup_points: Vec<ToyPoint>,
}

impl ToyGroup {
    fn build() -> Result<Self, String> {
        let roots = square_roots();
        let points = all_curve_points(&roots);
        if points.len() != TOY_GROUP_ORDER {
            return Err(format!(
                "toy curve has {} points, expected {TOY_GROUP_ORDER}",
                points.len()
            ));
        }
        let generator = points
            .iter()
            .copied()
            .filter(|point| *point != ToyPoint::Infinity)
            .map(|point| scalar_mul(point, COFACTOR))
            .find(|point| {
                *point != ToyPoint::Infinity
                    && scalar_mul(*point, SUBGROUP_ORDER) == ToyPoint::Infinity
            })
            .ok_or("failed to find order-149 generator")?;
        let mut subgroup_points = Vec::with_capacity(SUBGROUP_ORDER as usize);
        let mut point = ToyPoint::Infinity;
        for _ in 0..SUBGROUP_ORDER {
            subgroup_points.push(point);
            point = add_points(point, generator);
        }
        if point != ToyPoint::Infinity
            || subgroup_points
                .iter()
                .copied()
                .collect::<BTreeSet<_>>()
                .len()
                != SUBGROUP_ORDER as usize
        {
            return Err("order-149 subgroup enumeration failed".into());
        }
        Ok(Self {
            points,
            subgroup_generator: generator,
            subgroup_points,
        })
    }
}

fn mod_value(value: i16) -> u64 {
    if value >= 0 {
        value as u64 % SUBGROUP_ORDER
    } else {
        let magnitude = (-value) as u64 % SUBGROUP_ORDER;
        if magnitude == 0 {
            0
        } else {
            SUBGROUP_ORDER - magnitude
        }
    }
}

fn mod_mul(left: u64, right: u64, operations: &mut RankOperations) -> u64 {
    operations.multiplications += 1;
    ((left as u128 * right as u128) % SUBGROUP_ORDER as u128) as u64
}

fn mod_pow(mut base: u64, mut exponent: u64, operations: &mut RankOperations) -> u64 {
    let mut out = 1u64;
    while exponent != 0 {
        if exponent & 1 == 1 {
            out = mod_mul(out, base, operations);
        }
        base = mod_mul(base, base, operations);
        exponent >>= 1;
    }
    out
}

fn matrix_rank(rows: &[Vec<i16>]) -> RankMeasurement {
    if rows.is_empty() {
        return RankMeasurement {
            rank: 0,
            operations: RankOperations::default(),
        };
    }
    let columns = rows[0].len();
    let mut operations = RankOperations::default();
    let mut basis = BTreeMap::<usize, Vec<u64>>::new();
    for row in rows {
        let mut current = row
            .iter()
            .map(|&value| mod_value(value))
            .collect::<Vec<_>>();
        for (&pivot, pivot_row) in &basis {
            let factor = current[pivot];
            if factor == 0 {
                continue;
            }
            for column in pivot..columns {
                let product = mod_mul(factor, pivot_row[column], &mut operations);
                operations.additions_or_subtractions += 1;
                current[column] = (current[column] + SUBGROUP_ORDER - product) % SUBGROUP_ORDER;
            }
        }
        let Some(pivot) = current.iter().position(|&value| value != 0) else {
            continue;
        };
        operations.inversions += 1;
        let inverse = mod_pow(current[pivot], SUBGROUP_ORDER - 2, &mut operations);
        for value in &mut current[pivot..] {
            *value = mod_mul(*value, inverse, &mut operations);
        }
        basis.insert(pivot, current);
    }
    RankMeasurement {
        rank: basis.len(),
        operations,
    }
}

fn histogram_rows(histogram: BTreeMap<usize, u64>) -> Vec<HistogramRow> {
    histogram
        .into_iter()
        .map(|(value, count)| HistogramRow {
            value: value as u64,
            count,
        })
        .collect()
}

fn next_combination(combination: &mut [usize], limit: usize) -> bool {
    let width = combination.len();
    for position in (0..width).rev() {
        if combination[position] < limit - width + position {
            combination[position] += 1;
            for next in position + 1..width {
                combination[next] = combination[next - 1] + 1;
            }
            return true;
        }
    }
    false
}

#[derive(Clone)]
struct DestinationRow {
    coefficients: Vec<i16>,
    selected: Vec<(usize, i16)>,
    target: usize,
    point: ToyPoint,
}

fn destination_rows(base: &[u64], group: &ToyGroup) -> Result<Vec<DestinationRow>, String> {
    let mut rows = Vec::new();
    let mut combination = (0..ARITY).collect::<Vec<_>>();
    loop {
        for mask in 0..1usize << ARITY {
            let mut coefficients = vec![0i16; COLUMNS];
            let mut selected = Vec::with_capacity(ARITY);
            let mut scalar = 0i64;
            let mut point = ToyPoint::Infinity;
            for (position, &column) in combination.iter().enumerate() {
                let coefficient = if mask >> position & 1 == 0 { 1 } else { -1 };
                coefficients[column] = coefficient;
                selected.push((column, coefficient));
                scalar += coefficient as i64 * base[column] as i64;
                let addend = group.subgroup_points[base[column] as usize];
                point = add_points(
                    point,
                    if coefficient == 1 {
                        addend
                    } else {
                        negate(addend)
                    },
                );
            }
            let target = scalar.rem_euclid(SUBGROUP_ORDER as i64) as usize;
            if point != group.subgroup_points[target] {
                return Err("destination scalar/group replay mismatch".into());
            }
            rows.push(DestinationRow {
                coefficients,
                selected,
                target,
                point,
            });
        }
        if !next_combination(&mut combination, COLUMNS) {
            break;
        }
    }
    if rows.len() != 160 {
        return Err(format!(
            "destination row count {}, expected 160",
            rows.len()
        ));
    }
    Ok(rows)
}

fn deck_kernel_basis(fibre_size: usize) -> Vec<Vec<i16>> {
    let columns = COLUMNS * fibre_size;
    let mut rows = Vec::with_capacity(COLUMNS * fibre_size.saturating_sub(1));
    for column in 0..COLUMNS {
        for sheet in 1..fibre_size {
            let mut row = vec![0i16; columns];
            row[column * fibre_size] = -1;
            row[column * fibre_size + sheet] = 1;
            rows.push(row);
        }
    }
    rows
}

fn collapse_row(row: &[i16], fibre_size: usize) -> Vec<i16> {
    (0..COLUMNS)
        .map(|column| {
            row[column * fibre_size..(column + 1) * fibre_size]
                .iter()
                .sum()
        })
        .collect()
}

struct Bucket {
    lifted_rows: Vec<Vec<i16>>,
    unique_lifted_rows: BTreeSet<Vec<i16>>,
    collapsed_rows: BTreeSet<Vec<i16>>,
    source_sums: BTreeSet<ToyPoint>,
}

impl Bucket {
    fn new() -> Self {
        Self {
            lifted_rows: Vec::new(),
            unique_lifted_rows: BTreeSet::new(),
            collapsed_rows: BTreeSet::new(),
            source_sums: BTreeSet::new(),
        }
    }
}

fn evaluate_cell(base: &BaseSpec, scalar: u64, group: &ToyGroup) -> Result<CellResult, String> {
    let labels = base.folded_scalar_labels.clone();
    let destination = destination_rows(&labels, group)?;
    let kernel = group
        .points
        .iter()
        .copied()
        .filter(|&point| power_two_map(point, scalar) == ToyPoint::Infinity)
        .collect::<Vec<_>>();
    let fibre_size = kernel.len();
    if fibre_size == 0 {
        return Err("empty rational kernel".into());
    }
    let mut fibres = Vec::<Vec<ToyPoint>>::new();
    let mut fibre_histogram = BTreeMap::<usize, u64>::new();
    let mut missing_fibres = 0u64;
    for &label in &labels {
        let target = group.subgroup_points[label as usize];
        let fibre = group
            .points
            .iter()
            .copied()
            .filter(|&point| power_two_map(point, scalar) == target)
            .collect::<Vec<_>>();
        *fibre_histogram.entry(fibre.len()).or_default() += 1;
        if fibre.len() != fibre_size {
            missing_fibres += fibre_size.abs_diff(fibre.len()) as u64;
        }
        fibres.push(fibre);
    }
    let lifted_columns = COLUMNS * fibre_size;
    let kernel_basis = deck_kernel_basis(fibre_size);
    let kernel_measurement = matrix_rank(&kernel_basis);
    let expected_kernel_rank = COLUMNS * fibre_size.saturating_sub(1);
    if kernel_measurement.rank != expected_kernel_rank {
        return Err(format!(
            "deck kernel rank {}, expected {expected_kernel_rank}",
            kernel_measurement.rank
        ));
    }

    let mut buckets = BTreeMap::<usize, Bucket>::new();
    let mut replay_failures = 0u64;
    let mut duplicate_lift_rows = 0u64;
    let mut enumerated = 0u64;
    let mut hasher = Hasher::new();
    for row in &destination {
        let choices = fibre_size.pow(ARITY as u32);
        for ordinal in 0..choices {
            let mut quotient = ordinal;
            let mut lifted = vec![0i16; lifted_columns];
            let mut source_sum = ToyPoint::Infinity;
            for &(column, coefficient) in &row.selected {
                let sheet = quotient % fibre_size;
                quotient /= fibre_size;
                lifted[column * fibre_size + sheet] = coefficient;
                let source = fibres[column][sheet];
                source_sum = add_points(
                    source_sum,
                    if coefficient == 1 {
                        source
                    } else {
                        negate(source)
                    },
                );
            }
            let mapped = power_two_map(source_sum, scalar);
            if mapped != row.point || mapped != group.subgroup_points[row.target] {
                replay_failures += 1;
            }
            let collapsed = collapse_row(&lifted, fibre_size);
            if collapsed != row.coefficients {
                replay_failures += 1;
            }
            let bucket = buckets.entry(row.target).or_insert_with(Bucket::new);
            if !bucket.unique_lifted_rows.insert(lifted.clone()) {
                duplicate_lift_rows += 1;
            }
            bucket.collapsed_rows.insert(collapsed);
            bucket.source_sums.insert(source_sum);
            bucket.lifted_rows.push(lifted.clone());
            hasher.update(&(row.target as u64).to_be_bytes());
            for value in lifted {
                hasher.update(&value.to_be_bytes());
            }
            enumerated += 1;
        }
    }

    let expected = destination.len() as u64 * fibre_size.pow(ARITY as u32) as u64;
    let false_negatives = expected.saturating_sub(enumerated);
    let mut rank_operations = kernel_measurement.operations;
    let mut max_destination_occupancy = 0usize;
    let mut max_lifted_occupancy = 0usize;
    let mut max_source_sums = 0usize;
    let mut max_raw_rank = 0usize;
    let mut max_raw_homogeneous_rank = 0usize;
    let mut max_union_rank = 0usize;
    let mut max_homogeneous_union_rank = 0usize;
    let mut max_quotient_rank = 0usize;
    let mut max_quotient_homogeneous_rank = 0usize;
    let mut max_collapsed_rank = 0usize;
    let mut max_collapsed_homogeneous_rank = 0usize;
    let mut max_raw_removed = 0usize;
    let mut quotient_failures = 0u64;
    let mut cover_created_rank = 0u64;
    let mut coefficient_rank_histogram = BTreeMap::<usize, u64>::new();
    let mut quotient_rank_histogram = BTreeMap::<usize, u64>::new();
    let mut peak_logical_bytes = 0u64;
    for bucket in buckets.values() {
        let rows = &bucket.lifted_rows;
        let first = &rows[0];
        let homogeneous = rows
            .iter()
            .skip(1)
            .map(|row| {
                row.iter()
                    .zip(first)
                    .map(|(&value, &origin)| value - origin)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let collapsed = bucket.collapsed_rows.iter().cloned().collect::<Vec<_>>();
        let collapsed_first = &collapsed[0];
        let collapsed_homogeneous = collapsed
            .iter()
            .skip(1)
            .map(|row| {
                row.iter()
                    .zip(collapsed_first)
                    .map(|(&value, &origin)| value - origin)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut union = kernel_basis.clone();
        union.extend(rows.iter().cloned());
        let mut homogeneous_union = kernel_basis.clone();
        homogeneous_union.extend(homogeneous.iter().cloned());
        let raw = matrix_rank(rows);
        let raw_homogeneous = matrix_rank(&homogeneous);
        let union_measurement = matrix_rank(&union);
        let homogeneous_union_measurement = matrix_rank(&homogeneous_union);
        let collapsed_measurement = matrix_rank(&collapsed);
        let collapsed_homogeneous_measurement = matrix_rank(&collapsed_homogeneous);
        for measurement in [
            raw,
            raw_homogeneous,
            union_measurement,
            homogeneous_union_measurement,
            collapsed_measurement,
            collapsed_homogeneous_measurement,
        ] {
            rank_operations.add_assign(measurement.operations);
        }
        let quotient_rank = union_measurement.rank - kernel_measurement.rank;
        let quotient_homogeneous_rank =
            homogeneous_union_measurement.rank - kernel_measurement.rank;
        if quotient_rank != collapsed_measurement.rank
            || quotient_homogeneous_rank != collapsed_homogeneous_measurement.rank
        {
            quotient_failures += 1;
        }
        if quotient_rank > collapsed_measurement.rank {
            cover_created_rank += 1;
        }
        *coefficient_rank_histogram.entry(raw.rank).or_default() += 1;
        *quotient_rank_histogram.entry(quotient_rank).or_default() += 1;
        max_destination_occupancy = max_destination_occupancy.max(collapsed.len());
        max_lifted_occupancy = max_lifted_occupancy.max(rows.len());
        max_source_sums = max_source_sums.max(bucket.source_sums.len());
        max_raw_rank = max_raw_rank.max(raw.rank);
        max_raw_homogeneous_rank = max_raw_homogeneous_rank.max(raw_homogeneous.rank);
        max_union_rank = max_union_rank.max(union_measurement.rank);
        max_homogeneous_union_rank =
            max_homogeneous_union_rank.max(homogeneous_union_measurement.rank);
        max_quotient_rank = max_quotient_rank.max(quotient_rank);
        max_quotient_homogeneous_rank =
            max_quotient_homogeneous_rank.max(quotient_homogeneous_rank);
        max_collapsed_rank = max_collapsed_rank.max(collapsed_measurement.rank);
        max_collapsed_homogeneous_rank =
            max_collapsed_homogeneous_rank.max(collapsed_homogeneous_measurement.rank);
        max_raw_removed = max_raw_removed.max(raw.rank.saturating_sub(quotient_rank));
        let logical_bytes = (rows.len() * lifted_columns * std::mem::size_of::<i16>()
            + kernel_basis.len() * lifted_columns * std::mem::size_of::<i16>())
            as u64;
        peak_logical_bytes = peak_logical_bytes.max(logical_bytes);
    }

    Ok(CellResult {
        family: base.family.clone(),
        selection_uses_discrete_log_labels: base.selection_uses_discrete_log_labels,
        map_scalar: scalar,
        algebraic_degree: scalar * scalar,
        rational_kernel_size: kernel.len() as u64,
        rational_kernel_points: kernel.into_iter().map(point_label).collect(),
        fibre_size_histogram_for_six_destinations: histogram_rows(fibre_histogram),
        destination_columns: COLUMNS as u64,
        lifted_columns: lifted_columns as u64,
        arity: ARITY as u64,
        destination_raw_rows: destination.len() as u64,
        expected_lifted_rows: expected,
        enumerated_lifted_rows: enumerated,
        unique_lifted_rows: buckets
            .values()
            .map(|bucket| bucket.unique_lifted_rows.len() as u64)
            .sum(),
        lift_multiplier_per_destination_row: fibre_size.pow(ARITY as u32) as u64,
        nonempty_destination_targets: buckets.len() as u64,
        maximum_destination_row_occupancy: max_destination_occupancy as u64,
        maximum_lifted_row_occupancy: max_lifted_occupancy as u64,
        maximum_distinct_source_sums_per_destination_target: max_source_sums as u64,
        deck_kernel_rank: kernel_measurement.rank as u64,
        maximum_raw_lifted_coefficient_rank: max_raw_rank as u64,
        maximum_raw_lifted_homogeneous_rank: max_raw_homogeneous_rank as u64,
        maximum_union_with_kernel_rank: max_union_rank as u64,
        maximum_homogeneous_union_with_kernel_rank: max_homogeneous_union_rank as u64,
        maximum_quotient_coefficient_rank: max_quotient_rank as u64,
        maximum_quotient_homogeneous_rank: max_quotient_homogeneous_rank as u64,
        maximum_collapsed_destination_coefficient_rank: max_collapsed_rank as u64,
        maximum_collapsed_destination_homogeneous_rank: max_collapsed_homogeneous_rank as u64,
        maximum_raw_rank_removed_by_quotient: max_raw_removed as u64,
        buckets_with_quotient_identity_failure: quotient_failures,
        buckets_with_cover_created_destination_rank: cover_created_rank,
        missing_fibres,
        duplicate_lift_rows,
        replay_failures,
        false_positives: replay_failures,
        false_negatives,
        fibre_scan_pushforward_doublings: group.points.len() as u64
            * scalar.trailing_zeros() as u64
            * (COLUMNS as u64 + 1),
        destination_sum_group_additions: destination.len() as u64 * ARITY as u64,
        source_sum_group_additions: enumerated * ARITY as u64,
        row_pushforward_doublings: enumerated * scalar.trailing_zeros() as u64,
        rank_field_operations_mod_149: rank_operations,
        peak_logical_coefficient_bytes: peak_logical_bytes,
        disk_bytes_written_during_enumeration: 0,
        coefficient_rank_histogram: histogram_rows(coefficient_rank_histogram),
        quotient_rank_histogram: histogram_rows(quotient_rank_histogram),
        enumeration_blake3: hasher.finalize().to_hex().to_string(),
    })
}

fn base_specs(group: &ToyGroup) -> Vec<BaseSpec> {
    [
        ("affine-progression", true, [1, 3, 5, 7, 9, 11]),
        ("pair-sum-closure", true, [3, 32, 35, 38, 41, 70]),
        ("low-coordinate-geometry", false, [8, 41, 46, 49, 54, 59]),
        ("hash-coordinate-geometry", false, [4, 13, 14, 22, 40, 41]),
    ]
    .into_iter()
    .map(|(family, scalar_defined, labels)| BaseSpec {
        family: family.into(),
        selection_uses_discrete_log_labels: scalar_defined,
        folded_scalar_labels: labels.to_vec(),
        points: labels
            .iter()
            .map(|&label| point_label(group.subgroup_points[label as usize]))
            .collect(),
    })
    .collect()
}

fn obligation(name: &str, status: &str, evidence: &str, scope: &str) -> Obligation {
    Obligation {
        name: name.into(),
        status: status.into(),
        evidence: evidence.into(),
        scope: scope.into(),
    }
}

fn build_assessment(cells: &[CellResult]) -> TransferAssessment {
    let complete = cells.len() == 16;
    let exact = complete
        && cells.iter().all(|cell| {
            cell.replay_failures == 0
                && cell.missing_fibres == 0
                && cell.buckets_with_quotient_identity_failure == 0
        });
    let toy_path = TransferPath {
        name: "complete rational multiplication-map fibres".into(),
        source: "E(F_1151), #E=1192=8*149".into(),
        destination: "E(F_1151), order-149 subgroup H".into(),
        transformation_type: "curve endomorphisms [k], k in {1,2,4,8}, with induced group pushforward".into(),
        field_of_definition: "F_1151".into(),
        source_dimension: Some(1),
        destination_dimension: Some(1),
        relevant_subgroup: "H=< (360,129) >, |H|=149; [k]|H is injective".into(),
        obligations: vec![
            obligation("existence and formulas", "supported", "The standard multiplication maps are evaluated with complete affine group arithmetic.", "registered toy maps"),
            obligation("rational kernel and fibres", if exact { "supported" } else { "unknown" }, "Every one of 1,192 rational points is enumerated for every map and destination column.", "registered toy maps"),
            obligation("subgroup preservation", "supported", "gcd(k,149)=1, and every replayed source row pushes to its exact order-149 destination sum.", "registered toy maps"),
            obligation("quotient rank", if exact { "supported" } else { "unknown" }, "Each bucket checks dim((span(R)+ker(C))/ker(C))=rank(C(R)) over F_149.", "all registered buckets"),
            obligation("exceptional points", "not_applicable", "Complete group-law enumeration includes infinity; multiplication maps have no omitted rational exceptional branch.", "registered toy maps"),
            obligation("measured advantage", "refuted", "Raw sheet rows are removed by the collapse kernel and create no destination-rank gain.", "registered toy fibre multiplicity"),
        ],
        conclusion: "Rational sheet multiplicity is real in the cofactor-eight ambient group, but its additional coefficient directions lie in the pushforward kernel.".into(),
        exploration_boundary: "Complete for four multiplication maps, four frozen bases, and all signed arity-three rows on the registered toy curve.".into(),
    };
    let p256_path = TransferPath {
        name: "P-256 cover or Jacobian pushforward".into(),
        source: "unsupplied cover/Jacobian rational subgroup".into(),
        destination: format!("{CURVE_SLUG}, prime-order E(F_p)"),
        transformation_type: "proposed divisor/Jacobian pushforward group homomorphism".into(),
        field_of_definition: "must be supplied by a candidate".into(),
        source_dimension: None,
        destination_dimension: Some(1),
        relevant_subgroup: "destination subgroup has prime order n and cofactor one".into(),
        obligations: vec![
            obligation("explicit construction", "unknown", "No P-256 cover equation, Jacobian, divisor representation, or executable maps were supplied.", "construction-specific"),
            obligation("rational multiplication-map fibres", "supported", "For gcd(k,n)=1, [k] is a bijection on the prime-order rational group despite algebraic degree k^2.", "P-256 multiplication maps"),
            obligation("kernel-sheet row gain", "refuted", "Rows differing only inside a pushforward kernel have identical destination coefficient rows by the exact quotient theorem.", "all linearized relation rows under a group pushforward"),
            obligation("inverse recovery and exceptional points", "unknown", "Requires explicit formulas and fields for a proposed nontrivial cover.", "construction-specific"),
            obligation("destination decomposition advantage", "unknown", "The experiment does not exclude a cover that makes distinct collapsed decompositions cheaper to find.", "construction-specific"),
            obligation("paired end-to-end cost", "unknown", "No P-256 relation oracle or recovery implementation exists for this family.", "construction-specific"),
        ],
        conclusion: "Algebraic fibre degree cannot be booked as independent P-256 relation rank. A future cover must improve the search for distinct collapsed equations and price the complete maps and recovery.".into(),
        exploration_boundary: "Exact quotient theorem plus P-256 prime-order multiplication-map fact; not a universal search over covers or Jacobians.".into(),
    };
    let mut assessment = TransferAssessment {
        schema: "p256-cover-fiber-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 300,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        methodology_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/references/methodology.md (unavailable)".into(),
        template_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/assets/assessment-template.json (unavailable)".into(),
        typed_correspondence_graph: vec![
            "lifted coefficient rows V --[sheet collapse C]--> destination coefficient rows W".into(),
            "ker(C) --[pushforward]--> 0: deck/fibre differences carry zero destination equation".into(),
            "E(F_1151) --[[k], degree k^2]--> E(F_1151): rational kernel size k for k=1,2,4,8".into(),
            format!("{CURVE_SLUG}/E(F_p) --[[k], gcd(k,n)=1]--> E(F_p): rational bijection"),
            "future cover/Jacobian --[distinct collapsed decompositions plus recovery]--> possible advantage: open".into(),
        ],
        paths: vec![toy_path, p256_path],
        controls: vec![
            "Identity map [1] is the no-lift control.".into(),
            "Maps [2], [4], and [8] provide 2, 4, and 8 rational cofactor sheets on the complete toy ambient group.".into(),
            "Scalar-defined high-multiplicity bases and coordinate-defined bases are both included.".into(),
            "Every source sum is pushed forward with exact curve arithmetic and checked against its destination row.".into(),
        ],
        cost_accounting: vec![
            "Measured: fibre scans, every source-row sum, every row pushforward, modular-rank operations, logical coefficient bytes, process telemetry, and artifact bytes.".into(),
            "Exact analytic: P-256 prime-order rational bijectivity and the coefficient-space quotient identity.".into(),
            "Unset: construction and decomposition costs for any nontrivial P-256 cover, inverse recovery, structured degree, relation yield, sparse linear algebra, and full DLP cost.".into(),
        ],
        weakest_open_obligation: "Supply an explicit P-256 cover or Jacobian whose algorithm finds distinct collapsed destination decompositions more cheaply, then measure formulas, exceptional sets, inverse recovery, degree, yield, and complete one-target cost.".into(),
        narrowest_supported_finding: "Within the complete registered fibre family, every extra sheet-generated row direction lies in the pushforward kernel; quotient relation rank exactly equals collapsed destination rank.".into(),
        semantic_evidence_sha256: String::new(),
        assessment_json_bytes: 0,
    };
    let semantic = json!({
        "curve": assessment.curve,
        "typed_correspondence_graph": assessment.typed_correspondence_graph,
        "paths": assessment.paths,
        "controls": assessment.controls,
        "cost_accounting": assessment.cost_accounting,
        "weakest_open_obligation": assessment.weakest_open_obligation,
        "narrowest_supported_finding": assessment.narrowest_supported_finding,
    });
    assessment.semantic_evidence_sha256 =
        sha256_hex(&serde_json::to_vec(&semantic).expect("assessment semantic JSON"));
    assessment
}

fn boundary_table() -> Vec<BoundaryRow> {
    vec![
        BoundaryRow {
            variant: "Pollard rho".into(),
            effective_log_classes: None,
            cover_sheet_credit: 1,
            optimistic_ratio_to_rho: 1.0,
            scope: "matched generic reference".into(),
            result: "reference".into(),
        },
        BoundaryRow {
            variant: format!("Round 295 {FACTOR_BASE_ID}"),
            effective_log_classes: Some(164),
            cover_sheet_credit: 1,
            optimistic_ratio_to_rho: ROUND295_RATIO_TO_RHO,
            scope: "free-oracle K-th collision lower boundary; all later costs omitted".into(),
            result: "existing independent-log boundary".into(),
        },
        BoundaryRow {
            variant: "Round 300 fibre-lifted width-164 base".into(),
            effective_log_classes: Some(164),
            cover_sheet_credit: 1,
            optimistic_ratio_to_rho: ROUND295_RATIO_TO_RHO,
            scope: "extra sheets quotient to the same destination rows; no cover cost included"
                .into(),
            result: "no boundary movement".into(),
        },
        BoundaryRow {
            variant: format!("registered 17-term {REGISTERED_S17_FACTOR_BASE_ID}"),
            effective_log_classes: Some(131_458),
            cover_sheet_credit: 1,
            optimistic_ratio_to_rho: REGISTERED_S17_RATIO_TO_RHO,
            scope: "free-oracle independent-log comparison".into(),
            result: "unchanged".into(),
        },
    ]
}

fn build_result(cli: &Cli) -> Result<(ResultReceipt, TransferAssessment), String> {
    let (round299_dependency, round299) = checked_dependency(&cli.round299, ROUND299_SHA256)?;
    let (factor_base_dependency, factor_base) =
        checked_dependency(&cli.factor_base, ROUND295_FB_SHA256)?;
    if round299.pointer("/schema").and_then(Value::as_str) != Some("p256-transfer-anchor-screen/v1")
        || round299.pointer("/screening_round").and_then(Value::as_u64) != Some(299)
        || round299.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
    {
        return Err("Round 299 dependency schema or gates mismatch".into());
    }
    if factor_base.pointer("/schema").and_then(Value::as_str)
        != Some("ecbench.factor_base_dump/v1-wide")
        || factor_base
            .pointer("/factor_base/fb_id")
            .and_then(Value::as_str)
            != Some(FACTOR_BASE_ID)
        || factor_base
            .pointer("/factor_base/columns")
            .and_then(Value::as_u64)
            != Some(164)
    {
        return Err("Round 295 factor-base identity mismatch".into());
    }
    let p256 = CurveParams::p256();
    if p256.h.to_string() != "1"
        || p256.n.to_string()
            != "115792089210356248762697446949407573529996955224135760342422259061068512044369"
    {
        return Err("P-256 order or cofactor mismatch".into());
    }
    let group = ToyGroup::build()?;
    let bases = base_specs(&group);
    let mut cells = Vec::new();
    for base in &bases {
        for scalar in MAP_SCALARS {
            cells.push(evaluate_cell(base, scalar, &group)?);
        }
    }
    let total_destination_rows = cells
        .iter()
        .map(|cell| cell.destination_raw_rows)
        .sum::<u64>();
    let total_lifted_rows = cells
        .iter()
        .map(|cell| cell.enumerated_lifted_rows)
        .sum::<u64>();
    let total_fibre_scan_doublings = cells
        .iter()
        .map(|cell| cell.fibre_scan_pushforward_doublings)
        .sum::<u64>();
    let total_destination_additions = cells
        .iter()
        .map(|cell| cell.destination_sum_group_additions)
        .sum::<u64>();
    let total_source_additions = cells
        .iter()
        .map(|cell| cell.source_sum_group_additions)
        .sum::<u64>();
    let total_pushforward_doublings = cells
        .iter()
        .map(|cell| cell.row_pushforward_doublings)
        .sum::<u64>();
    let mut total_rank_operations = RankOperations::default();
    for cell in &cells {
        total_rank_operations.add_assign(cell.rank_field_operations_mod_149);
    }
    let all_cells = cells.len() == bases.len() * MAP_SCALARS.len();
    let all_replayed = cells.iter().all(|cell| cell.replay_failures == 0);
    let zero_errors = cells
        .iter()
        .all(|cell| cell.false_positives == 0 && cell.false_negatives == 0);
    let zero_fibre_errors = cells
        .iter()
        .all(|cell| cell.missing_fibres == 0 && cell.duplicate_lift_rows == 0);
    let quotient_exact = cells
        .iter()
        .all(|cell| cell.buckets_with_quotient_identity_failure == 0);
    let cover_created_rank = cells
        .iter()
        .any(|cell| cell.buckets_with_cover_created_destination_rank != 0);
    if !all_cells || !all_replayed || !zero_errors || !zero_fibre_errors || !quotient_exact {
        return Err("Round 300 correctness gate failed".into());
    }
    let assessment = build_assessment(&cells);
    let gates = Gates {
        dependency_hashes_checked: true,
        toy_curve_and_prime_subgroup_exact: group.points.len() == TOY_GROUP_ORDER
            && group.subgroup_points.len() == SUBGROUP_ORDER as usize,
        all_sixteen_cells_complete: all_cells,
        every_lifted_row_group_replayed: all_replayed,
        zero_false_positives_and_false_negatives: zero_errors,
        zero_missing_fibres_and_duplicates: zero_fibre_errors,
        quotient_identity_holds_everywhere: quotient_exact,
        cover_creates_two_destination_equations_from_one_destination_row: cover_created_rank,
        explicit_p256_cover_and_recovery_implemented: false,
        structured_residual_degree_at_most_5: false,
        parity_at_or_below_rho: false,
        per_usable_relation_below_2_103: false,
        complete_collection_below_2_120: false,
        projected_storage_below_2_50: false,
        promoted: false,
    };
    let maximum_raw_lift_multiplier = cells
        .iter()
        .map(|cell| cell.lift_multiplier_per_destination_row)
        .max()
        .unwrap_or(0);
    let maximum_raw_rank_gain = cells
        .iter()
        .map(|cell| {
            cell.maximum_raw_lifted_coefficient_rank
                .saturating_sub(cell.maximum_collapsed_destination_coefficient_rank)
        })
        .max()
        .unwrap_or(0);
    let maximum_quotient_rank_gain = cells
        .iter()
        .map(|cell| {
            cell.maximum_quotient_coefficient_rank
                .saturating_sub(cell.maximum_collapsed_destination_coefficient_rank)
        })
        .max()
        .unwrap_or(0);
    let classification = "cover-fibre-multiplicity/quotient-rank-conservation/promotion-negative";
    let obstruction = "The largest registered fibre supplies 512 lifted rows per destination row and increases raw lifted-column rank, but every added direction lies in the sheet-collapse kernel. Exact quotient rank equals collapsed destination rank in every complete bucket. On cofactor-one P-256, multiplication maps coprime to n are already bijections on rational points. No explicit cover supplies cheaper distinct collapsed decompositions, degree<=5, or complete cost below rho.";
    let decision = "Do not credit algebraic degree or raw sheet count as relation yield. Close multiplication-map and pure deck-kernel fibre multiplicity as a parity route. Keep open only a cover/Jacobian construction that demonstrably lowers the cost of finding distinct collapsed destination equations and includes executable inverse recovery and end-to-end accounting. Do not attempt an unplanted P-256 relation.";
    let semantic = json!({
        "curve": CURVE_SLUG,
        "factor_base": FACTOR_BASE_ID,
        "registered_s17_factor_base": REGISTERED_S17_FACTOR_BASE_ID,
        "dependencies": [&round299_dependency, &factor_base_dependency],
        "cells": &cells,
        "boundary_table": boundary_table(),
        "assessment": assessment.semantic_evidence_sha256,
        "gates": &gates,
        "classification": classification,
        "dominant_obstruction": obstruction,
        "decision": decision,
    });
    let semantic_evidence_sha256 =
        sha256_hex(&serde_json::to_vec(&semantic).map_err(|error| error.to_string())?);
    Ok((
        ResultReceipt {
            schema: "p256-cover-fiber-screen/v1".into(),
            curve: CURVE_SLUG.into(),
            screening_round: 300,
            execution_status: "complete".into(),
            factor_base: FACTOR_BASE_ID.into(),
            registered_s17_factor_base: REGISTERED_S17_FACTOR_BASE_ID.into(),
            dependencies: vec![round299_dependency, factor_base_dependency],
            theorem: "For sheet collapse C and K=ker(C), dim((span(R)+K)/K)=rank(C(R)); rows differing only by fibre choices produce one destination equation.".into(),
            toy_curve: "E/F_1151: y^2=x^3-3x+241".into(),
            toy_group_order: TOY_GROUP_ORDER as u64,
            cofactor: COFACTOR,
            prime_subgroup_order: SUBGROUP_ORDER,
            subgroup_generator: point_label(group.subgroup_generator),
            ambient_points_enumerated: group.points.len() as u64,
            factor_base_controls: bases,
            cells,
            total_destination_rows,
            total_lifted_rows_replayed: total_lifted_rows,
            total_fibre_scan_pushforward_doublings: total_fibre_scan_doublings,
            total_destination_sum_group_additions: total_destination_additions,
            total_source_sum_group_additions: total_source_additions,
            total_row_pushforward_doublings: total_pushforward_doublings,
            total_rank_field_operations_mod_149: total_rank_operations,
            maximum_raw_lift_multiplier,
            maximum_observed_raw_rank_gain: maximum_raw_rank_gain,
            maximum_observed_quotient_rank_gain: maximum_quotient_rank_gain,
            p256_prime_order_rational_fibre_statement: "E(F_p) has prime order n and cofactor one; [k] with gcd(k,n)=1 is a rational-point bijection even when the algebraic map has degree k^2.".into(),
            boundary_table: boundary_table(),
            relations_reported_on_p256: 0,
            full_depth_unplanted_p256_relation_attempted: false,
            transfer_assessment_semantic_sha256: assessment.semantic_evidence_sha256.clone(),
            gates,
            classification: classification.into(),
            dominant_obstruction: obstruction.into(),
            decision: decision.into(),
            semantic_evidence_sha256,
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
        "round 300: {} cells, {} lifted rows, max lift={}x, quotient gain={}, promoted={}",
        result.cells.len(),
        result.total_lifted_rows_replayed,
        result.maximum_raw_lift_multiplier,
        result.maximum_observed_quotient_rank_gain,
        result.gates.promoted
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_cover_fiber_screen: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complete_toy_group_and_fibres() {
        let group = ToyGroup::build().expect("toy group");
        assert_eq!(group.points.len(), 1192);
        assert_eq!(group.subgroup_points.len(), 149);
        for scalar in MAP_SCALARS {
            let kernel = group
                .points
                .iter()
                .filter(|&&point| power_two_map(point, scalar) == ToyPoint::Infinity)
                .count();
            assert_eq!(kernel, scalar as usize);
        }
    }

    #[test]
    fn quotient_rank_identity_on_explicit_rows() {
        let fibre_size = 4;
        let kernel = deck_kernel_basis(fibre_size);
        let mut rows = Vec::new();
        for sheet in 0..fibre_size {
            let mut row = vec![0i16; COLUMNS * fibre_size];
            row[sheet] = 1;
            row[fibre_size + sheet] = -1;
            rows.push(row);
        }
        let mut union = kernel.clone();
        union.extend(rows.clone());
        let quotient = matrix_rank(&union).rank - matrix_rank(&kernel).rank;
        let collapsed = rows
            .iter()
            .map(|row| collapse_row(row, fibre_size))
            .collect::<Vec<_>>();
        assert_eq!(quotient, matrix_rank(&collapsed).rank);
        assert_eq!(quotient, 1);
    }

    #[test]
    fn complete_identity_cell() {
        let group = ToyGroup::build().expect("toy group");
        let base = base_specs(&group).remove(0);
        let cell = evaluate_cell(&base, 1, &group).expect("identity cell");
        assert_eq!(cell.enumerated_lifted_rows, 160);
        assert_eq!(cell.rational_kernel_size, 1);
        assert_eq!(cell.replay_failures, 0);
        assert_eq!(cell.buckets_with_quotient_identity_failure, 0);
    }
}
