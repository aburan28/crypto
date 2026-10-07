//! Round 299: prime-order transfer-anchor and structured multiplicity screen.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use blake3::Hasher;
use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde::Serialize;
use serde_json::{json, Value};

const ROUND298_SHA256: &str = "8718bbe26751ff0191164b0665618226a97556c64c6e2feeef72468a38206965";
const FACTOR_BASE_ID: &str = "FB1hc72514a2a8d3";
const P: u64 = 1151;
const A: u64 = P - 3;
const B: u64 = 241;
const TOY_GROUP_ORDER: usize = 1192;
const SUBGROUP_ORDER: u64 = 149;
const COFACTOR: u64 = 8;
const PRIMARY_COLUMNS: usize = 6;
const PRIMARY_ARITY: usize = 3;
const REFERENCE_COLUMNS: usize = 4;
const REFERENCE_ARITY: usize = 2;
const FORMAL_RANK_MODULUS: u64 = (1u64 << 61) - 1;
const HASH_SCALAR_BASES: u32 = 4096;
const P256_COLUMNS: u64 = 164;
const RHO_S: f64 = 1.3;

#[derive(Parser)]
#[command(about = "Audit prime-order transport and structured relation multiplicity")]
struct Cli {
    #[arg(long)]
    round298: PathBuf,
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

#[derive(Clone, Serialize)]
struct HistogramRow {
    value: u64,
    count: u64,
}

#[derive(Clone, Copy, Default, Serialize)]
struct FieldOperations {
    additions_or_subtractions: u64,
    multiplications: u64,
    inversions: u64,
}

impl FieldOperations {
    fn add_assign(&mut self, other: Self) {
        self.additions_or_subtractions += other.additions_or_subtractions;
        self.multiplications += other.multiplications;
        self.inversions += other.inversions;
    }
}

#[derive(Clone, Copy)]
struct RankMeasurement {
    rank: usize,
    operations: FieldOperations,
}

#[derive(Clone, Serialize)]
struct BaseMetrics {
    raw_rows: u64,
    distinct_normalized_rows: u64,
    nonempty_target_buckets: u64,
    maximum_distinct_occupancy: u64,
    maximum_occupancy_target: u64,
    maximum_known_target_coefficient_rank_mod_149: u64,
    maximum_homogeneous_rank_mod_149: u64,
    maximum_unknown_target_augmented_rank_mod_149: u64,
    maximum_formal_coefficient_rank_mod_2_61_minus_1: u64,
    maximum_formal_homogeneous_rank_mod_2_61_minus_1: u64,
    maximum_formal_augmented_rank_mod_2_61_minus_1: u64,
    nonzero_buckets_with_full_known_target_rank: u64,
    buckets_with_homogeneous_rank_b_minus_1: u64,
    buckets_with_formal_rank_divergence: u64,
    subgroup_kernel_failures: u64,
    subgroup_rank_ceiling_violations: u64,
    replay_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    scalar_label_additions: u64,
    logical_group_additions: u64,
    rank_field_operations: FieldOperations,
    occupancy_histogram_all_75_targets: Vec<HistogramRow>,
    coefficient_rank_histogram_nonempty_targets: Vec<HistogramRow>,
    homogeneous_rank_histogram_nonempty_targets: Vec<HistogramRow>,
    best_bucket_rows: Vec<Vec<i8>>,
}

#[derive(Clone, Serialize)]
struct BaseReport {
    folded_scalar_labels: Vec<u64>,
    canonical_points: Vec<String>,
    metrics: BaseMetrics,
}

#[derive(Serialize)]
struct FamilySummary {
    name: String,
    selection_uses_discrete_log_labels: bool,
    proposals: u64,
    admissible_proposals: u64,
    unique_bases: u64,
    raw_rows_replayed: u64,
    distinct_rows_checked: u64,
    bases_with_full_known_target_rank: u64,
    bases_with_homogeneous_rank_b_minus_1: u64,
    bases_with_subgroup_ceiling_violation: u64,
    bases_with_formal_rank_divergence: u64,
    subgroup_kernel_failures: u64,
    replay_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    scalar_label_additions: u64,
    rank_field_operations: FieldOperations,
    maximum_distinct_occupancy: u64,
    max_occupancy_histogram: Vec<HistogramRow>,
    max_known_rank_histogram: Vec<HistogramRow>,
    best_base: Option<BaseReport>,
}

#[derive(Serialize)]
struct ReferenceSummary {
    columns: u64,
    arity: u64,
    all_factor_bases: u64,
    raw_rows_per_base: u64,
    distinct_rows_per_base: u64,
    raw_rows_replayed: u64,
    distinct_rows_checked: u64,
    maximum_distinct_occupancy: u64,
    maximum_known_target_rank_mod_149: u64,
    maximum_homogeneous_rank_mod_149: u64,
    maximum_formal_known_target_rank_mod_2_61_minus_1: u64,
    maximum_formal_homogeneous_rank_mod_2_61_minus_1: u64,
    bases_with_full_known_target_rank: u64,
    bases_with_formal_rank_divergence: u64,
    subgroup_kernel_failures: u64,
    subgroup_rank_ceiling_violations: u64,
    replay_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    scalar_label_additions: u64,
    logical_group_additions: u64,
    rank_field_operations: FieldOperations,
    enumeration_blake3: String,
}

#[derive(Serialize)]
struct P256BoundaryCorrection {
    columns: u64,
    prior_round_distinct_rows_required: u64,
    known_nonzero_target_distinct_rows_required: u64,
    equivalent_raw_folded_preimages_required: u64,
    homogeneous_same_target_rank_ceiling: u64,
    exact_global_negation_quotient_domain: String,
    canonical_target_space: String,
    mean_distinct_rows_per_canonical_target: f64,
    random_map_null_log2_union_bound_any_bucket_at_164: f64,
    random_map_null_log2_poisson_expected_buckets_at_164: f64,
    random_map_null_is_measurement: bool,
    one_event_ratio_to_rho_before_omitted_costs: f64,
    two_event_ratio_to_rho_before_omitted_costs: f64,
}

#[derive(Serialize)]
struct Gates {
    dependency_hash_checked: bool,
    toy_group_and_prime_subgroup_exact: bool,
    all_registered_factor_base_families_complete: bool,
    exhaustive_b4_reference_complete: bool,
    every_raw_row_group_replayed: bool,
    zero_false_positives_and_false_negatives: bool,
    subgroup_rank_ceilings_hold: bool,
    independent_anchor_eliminated_by_homomorphic_transfer: bool,
    geometry_defined_sparse_full_rank_candidate_observed: bool,
    p256_geometry_construction_and_scaling_measured: bool,
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
    factor_base: String,
    screening_round: u64,
    execution_status: String,
    dependency: Dependency,
    p256_boundary_correction: P256BoundaryCorrection,
    toy_curve: String,
    toy_group_order: u64,
    cofactor: u64,
    prime_subgroup_order: u64,
    subgroup_generator: String,
    subgroup_generator_source_point: String,
    subgroup_points_enumerated: u64,
    folded_columns: u64,
    primary_columns: u64,
    primary_arity: u64,
    raw_rows_per_primary_base: u64,
    distinct_rows_per_primary_base: u64,
    canonical_subgroup_targets: u64,
    fixed_distinct_mean_occupancy: f64,
    family_summaries: Vec<FamilySummary>,
    global_unique_primary_bases: u64,
    global_primary_raw_rows_replayed: u64,
    global_primary_distinct_rows_checked: u64,
    global_primary_scalar_label_additions: u64,
    global_primary_logical_group_additions: u64,
    global_primary_rank_field_operations: FieldOperations,
    transition_table_exact_group_additions: u64,
    transition_table_field_operations: FieldOperations,
    primary_enumeration_blake3: String,
    exhaustive_reference: ReferenceSummary,
    transfer_assessment_semantic_sha256: String,
    relations_reported_on_p256: u64,
    full_depth_unplanted_p256_relation_attempted: bool,
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
struct TransferCase {
    name: String,
    source: String,
    destination: String,
    transformation_type: String,
    field_of_definition: String,
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
    theorem_cases: Vec<TransferCase>,
    controls: Vec<String>,
    cost_accounting: Vec<String>,
    weakest_open_obligation: String,
    narrowest_supported_finding: String,
    semantic_evidence_sha256: String,
    assessment_json_bytes: u64,
}

#[derive(Default)]
struct FamilyBuild {
    scalar_defined: bool,
    proposals: u64,
    admissible: u64,
    bases: BTreeSet<Vec<u16>>,
}

#[derive(Default)]
struct FamilyAggregate {
    full_rank: u64,
    homogeneous_b_minus_1: u64,
    ceiling_violations: u64,
    formal_rank_divergences: u64,
    raw_rows: u64,
    distinct_rows: u64,
    kernel_failures: u64,
    replay_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    scalar_label_additions: u64,
    rank_field_operations: FieldOperations,
    maximum_occupancy: u64,
    occupancy_histogram: BTreeMap<usize, u64>,
    rank_histogram: BTreeMap<usize, u64>,
    best: Option<(Vec<u16>, BaseMetrics)>,
}

fn checked_dependency(path: &Path) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = sha256_hex(&bytes);
    if digest != ROUND298_SHA256 {
        return Err(format!(
            "{} SHA-256 is {digest}, expected {ROUND298_SHA256}",
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

struct PrimeSubgroup {
    source: ToyPoint,
    generator: ToyPoint,
    points: Vec<ToyPoint>,
    transitions: Vec<Vec<usize>>,
    canonical_points: Vec<ToyPoint>,
    exact_transition_additions: u64,
    transition_field_operations: FieldOperations,
}

impl PrimeSubgroup {
    fn build(roots: &[Vec<u64>]) -> Result<Self, String> {
        let curve_points = all_curve_points(roots);
        if curve_points.len() != TOY_GROUP_ORDER {
            return Err(format!(
                "toy curve has {} points, expected {TOY_GROUP_ORDER}",
                curve_points.len()
            ));
        }
        let (source, generator) = curve_points
            .iter()
            .copied()
            .filter(|point| *point != ToyPoint::Infinity)
            .find_map(|point| {
                let candidate = scalar_mul(point, COFACTOR);
                (candidate != ToyPoint::Infinity
                    && scalar_mul(candidate, SUBGROUP_ORDER) == ToyPoint::Infinity)
                    .then_some((point, candidate))
            })
            .ok_or("failed to find order-149 toy subgroup generator")?;
        let mut points = Vec::with_capacity(SUBGROUP_ORDER as usize);
        let mut point = ToyPoint::Infinity;
        for _ in 0..SUBGROUP_ORDER {
            points.push(point);
            point = add_points(point, generator);
        }
        if point != ToyPoint::Infinity
            || points.iter().copied().collect::<BTreeSet<_>>().len() != SUBGROUP_ORDER as usize
        {
            return Err("order-149 subgroup enumeration is not exact".into());
        }
        let index = points
            .iter()
            .copied()
            .enumerate()
            .map(|(label, point)| (point, label))
            .collect::<BTreeMap<_, _>>();
        let mut transitions = Vec::with_capacity(2 * 74);
        for label in 1..=74usize {
            for signed in [label, SUBGROUP_ORDER as usize - label] {
                let addend = points[signed];
                let mut row = Vec::with_capacity(points.len());
                for &state in &points {
                    let sum = add_points(state, addend);
                    row.push(
                        *index
                            .get(&sum)
                            .ok_or("subgroup transition escaped subgroup")?,
                    );
                }
                transitions.push(row);
            }
        }
        let canonical_points = (1..=74usize)
            .map(|label| points[label].min(negate(points[label])))
            .collect::<Vec<_>>();
        let signed_addends = 2 * 74u64;
        let general_additions = signed_addends * (SUBGROUP_ORDER - 3);
        let doublings = signed_addends;
        Ok(Self {
            source,
            generator,
            points,
            transitions,
            canonical_points,
            exact_transition_additions: (2 * 74 * SUBGROUP_ORDER as usize) as u64,
            transition_field_operations: FieldOperations {
                additions_or_subtractions: 6 * general_additions + 5 * doublings,
                multiplications: 21 * general_additions + 24 * doublings,
                inversions: general_additions + doublings,
            },
        })
    }

    fn replay_row(&self, base: &[u16], row: &[i8]) -> Result<usize, String> {
        let mut state = 0usize;
        for (&label, &coefficient) in base.iter().zip(row) {
            let signed = match coefficient {
                0 => continue,
                1 => 2 * (label as usize - 1),
                -1 => 2 * (label as usize - 1) + 1,
                _ => return Err("coefficient outside {-1,0,1}".into()),
            };
            state = self.transitions[signed][state];
        }
        Ok(state)
    }
}

fn fold_scalar(value: u64) -> Option<u16> {
    let value = value % SUBGROUP_ORDER;
    if value == 0 {
        None
    } else {
        Some(value.min(SUBGROUP_ORDER - value) as u16)
    }
}

fn canonical_base(values: impl IntoIterator<Item = u64>, width: usize) -> Option<Vec<u16>> {
    let mut labels = values
        .into_iter()
        .map(fold_scalar)
        .collect::<Option<Vec<_>>>()?;
    labels.sort_unstable();
    labels.dedup();
    (labels.len() == width).then_some(labels)
}

fn propose(build: &mut FamilyBuild, values: impl IntoIterator<Item = u64>, width: usize) {
    build.proposals += 1;
    if let Some(base) = canonical_base(values, width) {
        build.admissible += 1;
        build.bases.insert(base);
    }
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

fn combinations_from_pool(pool: &[u16], width: usize) -> Vec<Vec<u16>> {
    let mut out = Vec::new();
    let mut combination = (0..width).collect::<Vec<_>>();
    loop {
        let mut base = combination
            .iter()
            .map(|&index| pool[index])
            .collect::<Vec<_>>();
        base.sort_unstable();
        out.push(base);
        if !next_combination(&mut combination, pool.len()) {
            break;
        }
    }
    out
}

fn hash_rank(domain: &[u8], seed: u32, value: u16) -> [u8; 32] {
    let mut bytes = Vec::with_capacity(domain.len() + 6);
    bytes.extend(domain);
    bytes.extend(seed.to_be_bytes());
    bytes.extend(value.to_be_bytes());
    sha256(&bytes)
}

fn build_families(subgroup: &PrimeSubgroup) -> BTreeMap<String, FamilyBuild> {
    let mut families = BTreeMap::<String, FamilyBuild>::new();

    let mut affine = FamilyBuild {
        scalar_defined: true,
        ..FamilyBuild::default()
    };
    for a in 1..SUBGROUP_ORDER {
        for d in 1..SUBGROUP_ORDER {
            propose(
                &mut affine,
                (0..PRIMARY_COLUMNS).map(|j| a + j as u64 * d),
                PRIMARY_COLUMNS,
            );
        }
    }
    families.insert("affine-progression".into(), affine);

    let mut geometric = FamilyBuild {
        scalar_defined: true,
        ..FamilyBuild::default()
    };
    for a in 1..SUBGROUP_ORDER {
        for ratio in 1..SUBGROUP_ORDER {
            let mut value = a;
            let mut labels = Vec::with_capacity(PRIMARY_COLUMNS);
            for _ in 0..PRIMARY_COLUMNS {
                labels.push(value);
                value = value * ratio % SUBGROUP_ORDER;
            }
            propose(&mut geometric, labels, PRIMARY_COLUMNS);
        }
    }
    families.insert("geometric-progression".into(), geometric);

    let mut pair_sum = FamilyBuild {
        scalar_defined: true,
        ..FamilyBuild::default()
    };
    for a in 1..=74u64 {
        for b in a + 1..=74u64 {
            for c in b + 1..=74u64 {
                propose(
                    &mut pair_sum,
                    [a, b, c, a + b, a + c, b + c],
                    PRIMARY_COLUMNS,
                );
            }
        }
    }
    families.insert("pair-sum-closure".into(), pair_sum);

    let mut hash_scalar = FamilyBuild {
        scalar_defined: true,
        ..FamilyBuild::default()
    };
    for seed in 0..HASH_SCALAR_BASES {
        let mut ranked = (1..=74u16)
            .map(|label| (hash_rank(b"round299/hash-scalar/v1", seed, label), label))
            .collect::<Vec<_>>();
        ranked.sort_unstable();
        propose(
            &mut hash_scalar,
            ranked
                .iter()
                .take(PRIMARY_COLUMNS)
                .map(|(_, label)| *label as u64),
            PRIMARY_COLUMNS,
        );
    }
    families.insert("hash-scalar".into(), hash_scalar);

    let mut low_geometry_labels = (1..=74u16).collect::<Vec<_>>();
    low_geometry_labels.sort_by_key(|&label| {
        let point = subgroup.canonical_points[label as usize - 1];
        match point {
            ToyPoint::Infinity => (0, 0, label),
            ToyPoint::Affine(point) => (point.x, point.y, label),
        }
    });
    let low_pool = low_geometry_labels[..12].to_vec();
    let low_bases = combinations_from_pool(&low_pool, PRIMARY_COLUMNS);
    let mut low_geometry = FamilyBuild::default();
    for base in low_bases {
        propose(
            &mut low_geometry,
            base.into_iter().map(u64::from),
            PRIMARY_COLUMNS,
        );
    }
    families.insert("low-coordinate-geometry".into(), low_geometry);

    let mut hash_geometry_labels = (1..=74u16)
        .map(|label| {
            let point = subgroup.canonical_points[label as usize - 1];
            let mut bytes = b"round299/geometry/v1".to_vec();
            match point {
                ToyPoint::Infinity => bytes.push(0),
                ToyPoint::Affine(point) => {
                    bytes.push(1);
                    bytes.extend(point.x.to_be_bytes());
                    bytes.extend(point.y.to_be_bytes());
                }
            }
            (sha256(&bytes), label)
        })
        .collect::<Vec<_>>();
    hash_geometry_labels.sort_unstable();
    let hash_pool = hash_geometry_labels
        .iter()
        .take(12)
        .map(|(_, label)| *label)
        .collect::<Vec<_>>();
    let hash_bases = combinations_from_pool(&hash_pool, PRIMARY_COLUMNS);
    let mut hash_geometry = FamilyBuild::default();
    for base in hash_bases {
        propose(
            &mut hash_geometry,
            base.into_iter().map(u64::from),
            PRIMARY_COLUMNS,
        );
    }
    families.insert("hash-coordinate-geometry".into(), hash_geometry);

    families
}

fn mod_value(value: i16, modulus: u64) -> u64 {
    if value >= 0 {
        value as u64 % modulus
    } else {
        let magnitude = (-value) as u64 % modulus;
        if magnitude == 0 {
            0
        } else {
            modulus - magnitude
        }
    }
}

fn mod_mul_counted(left: u64, right: u64, modulus: u64, operations: &mut FieldOperations) -> u64 {
    operations.multiplications += 1;
    ((left as u128 * right as u128) % modulus as u128) as u64
}

fn mod_pow_counted(
    mut base: u64,
    mut exponent: u64,
    modulus: u64,
    operations: &mut FieldOperations,
) -> u64 {
    let mut out = 1u64;
    while exponent != 0 {
        if exponent & 1 == 1 {
            out = mod_mul_counted(out, base, modulus, operations);
        }
        base = mod_mul_counted(base, base, modulus, operations);
        exponent >>= 1;
    }
    out
}

fn matrix_rank(rows: &[Vec<i16>], modulus: u64) -> RankMeasurement {
    if rows.is_empty() {
        return RankMeasurement {
            rank: 0,
            operations: FieldOperations::default(),
        };
    }
    let mut operations = FieldOperations::default();
    let columns = rows[0].len();
    let mut matrix = rows
        .iter()
        .map(|row| {
            row.iter()
                .map(|&value| mod_value(value, modulus))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut rank = 0usize;
    for column in 0..columns {
        let Some(pivot) = (rank..matrix.len()).find(|&row| matrix[row][column] != 0) else {
            continue;
        };
        matrix.swap(rank, pivot);
        operations.inversions += 1;
        let inverse = mod_pow_counted(matrix[rank][column], modulus - 2, modulus, &mut operations);
        for value in &mut matrix[rank][column..] {
            *value = mod_mul_counted(*value, inverse, modulus, &mut operations);
        }
        let pivot_row = matrix[rank].clone();
        for row in rank + 1..matrix.len() {
            let factor = matrix[row][column];
            if factor == 0 {
                continue;
            }
            for next in column..columns {
                let product = mod_mul_counted(factor, pivot_row[next], modulus, &mut operations);
                operations.additions_or_subtractions += 1;
                matrix[row][next] = (matrix[row][next] + modulus - product) % modulus;
            }
        }
        rank += 1;
        if rank == matrix.len() {
            break;
        }
    }
    RankMeasurement { rank, operations }
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

fn binomial_u64(n: usize, k: usize) -> u64 {
    let k = k.min(n - k);
    let mut value = 1u64;
    for i in 0..k {
        value = value * (n - i) as u64 / (i + 1) as u64;
    }
    value
}

fn evaluate_base(
    base: &[u16],
    arity: usize,
    subgroup: &PrimeSubgroup,
) -> Result<BaseMetrics, String> {
    let columns = base.len();
    let expected_raw = (1u64 << arity) * binomial_u64(columns, arity);
    let mut buckets = BTreeMap::<u16, BTreeSet<Vec<i8>>>::new();
    let mut combination = (0..arity).collect::<Vec<_>>();
    let mut raw_rows = 0u64;
    let mut replay_failures = 0u64;
    loop {
        for mask in 0..1u64 << arity {
            let mut row = vec![0i8; columns];
            let mut scalar = 0i64;
            for (position, &column) in combination.iter().enumerate() {
                let negative = mask >> position & 1 == 1;
                row[column] = if negative { -1 } else { 1 };
                let signed = if negative {
                    -(base[column] as i64)
                } else {
                    base[column] as i64
                };
                scalar += signed;
            }
            let mut scalar = scalar.rem_euclid(SUBGROUP_ORDER as i64) as u16;
            if scalar as u64 > SUBGROUP_ORDER / 2 {
                scalar = SUBGROUP_ORDER as u16 - scalar;
                for coefficient in &mut row {
                    *coefficient = -*coefficient;
                }
            } else if scalar == 0 && row.iter().find(|&&coefficient| coefficient != 0) == Some(&-1)
            {
                for coefficient in &mut row {
                    *coefficient = -*coefficient;
                }
            }
            let state = subgroup.replay_row(base, &row)?;
            if state != scalar as usize {
                replay_failures += 1;
            }
            buckets.entry(scalar).or_default().insert(row);
            raw_rows += 1;
        }
        if !next_combination(&mut combination, columns) {
            break;
        }
    }
    let distinct_rows = buckets.values().map(|rows| rows.len() as u64).sum::<u64>();
    if raw_rows != expected_raw || distinct_rows * 2 != raw_rows {
        return Err(format!(
            "base {base:?}, m={arity}: raw={raw_rows}, distinct={distinct_rows}, expected={expected_raw}"
        ));
    }

    let mut occupancy_histogram = BTreeMap::<usize, u64>::new();
    let mut coefficient_histogram = BTreeMap::<usize, u64>::new();
    let mut homogeneous_histogram = BTreeMap::<usize, u64>::new();
    let mut max_occupancy = 0usize;
    let mut max_target = 0u16;
    let mut max_rows = Vec::new();
    let mut max_coefficient_rank = 0usize;
    let mut max_homogeneous_rank = 0usize;
    let mut max_augmented_rank = 0usize;
    let mut max_formal_coefficient_rank = 0usize;
    let mut max_formal_homogeneous_rank = 0usize;
    let mut max_formal_augmented_rank = 0usize;
    let mut full_rank_nonzero = 0u64;
    let mut homogeneous_b_minus_1 = 0u64;
    let mut formal_rank_divergences = 0u64;
    let mut kernel_failures = 0u64;
    let mut ceiling_violations = 0u64;
    let mut rank_field_operations = FieldOperations::default();
    for (&target, rows) in &buckets {
        let rows = rows.iter().cloned().collect::<Vec<_>>();
        let coefficient = rows
            .iter()
            .map(|row| row.iter().map(|&value| value as i16).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let homogeneous = if target == 0 {
            coefficient.clone()
        } else {
            let first = &rows[0];
            rows.iter()
                .skip(1)
                .map(|row| {
                    row.iter()
                        .zip(first)
                        .map(|(&value, &base_value)| value as i16 - base_value as i16)
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        };
        let augmented = rows
            .iter()
            .map(|row| {
                let mut values = row.iter().map(|&value| value as i16).collect::<Vec<_>>();
                values.push(-1);
                values
            })
            .collect::<Vec<_>>();
        for row in &homogeneous {
            let dot = row
                .iter()
                .zip(base)
                .map(|(&coefficient, &label)| mod_value(coefficient, SUBGROUP_ORDER) * label as u64)
                .sum::<u64>()
                % SUBGROUP_ORDER;
            if dot != 0 {
                kernel_failures += 1;
            }
        }
        for row in &augmented {
            let dot = row[..columns]
                .iter()
                .zip(base)
                .map(|(&coefficient, &label)| mod_value(coefficient, SUBGROUP_ORDER) * label as u64)
                .sum::<u64>()
                % SUBGROUP_ORDER;
            if !(dot + SUBGROUP_ORDER - target as u64).is_multiple_of(SUBGROUP_ORDER) {
                kernel_failures += 1;
            }
        }
        let coefficient_measurement = matrix_rank(&coefficient, SUBGROUP_ORDER);
        let homogeneous_measurement = matrix_rank(&homogeneous, SUBGROUP_ORDER);
        let augmented_measurement = matrix_rank(&augmented, SUBGROUP_ORDER);
        let formal_coefficient_measurement = matrix_rank(&coefficient, FORMAL_RANK_MODULUS);
        let formal_homogeneous_measurement = matrix_rank(&homogeneous, FORMAL_RANK_MODULUS);
        let formal_augmented_measurement = matrix_rank(&augmented, FORMAL_RANK_MODULUS);
        for measurement in [
            coefficient_measurement,
            homogeneous_measurement,
            augmented_measurement,
            formal_coefficient_measurement,
            formal_homogeneous_measurement,
            formal_augmented_measurement,
        ] {
            rank_field_operations.add_assign(measurement.operations);
        }
        let coefficient_rank = coefficient_measurement.rank;
        let homogeneous_rank = homogeneous_measurement.rank;
        let augmented_rank = augmented_measurement.rank;
        let formal_coefficient_rank = formal_coefficient_measurement.rank;
        let formal_homogeneous_rank = formal_homogeneous_measurement.rank;
        let formal_augmented_rank = formal_augmented_measurement.rank;
        if coefficient_rank != formal_coefficient_rank
            || homogeneous_rank != formal_homogeneous_rank
            || augmented_rank != formal_augmented_rank
        {
            formal_rank_divergences += 1;
        }
        if homogeneous_rank > columns.saturating_sub(1) || augmented_rank > columns {
            ceiling_violations += 1;
        }
        if target != 0 && coefficient_rank == columns {
            full_rank_nonzero += 1;
        }
        if homogeneous_rank == columns.saturating_sub(1) {
            homogeneous_b_minus_1 += 1;
        }
        *occupancy_histogram.entry(rows.len()).or_default() += 1;
        *coefficient_histogram.entry(coefficient_rank).or_default() += 1;
        *homogeneous_histogram.entry(homogeneous_rank).or_default() += 1;
        if rows.len() > max_occupancy || (rows.len() == max_occupancy && target < max_target) {
            max_occupancy = rows.len();
            max_target = target;
            max_rows = rows.clone();
        }
        max_coefficient_rank = max_coefficient_rank.max(coefficient_rank);
        max_homogeneous_rank = max_homogeneous_rank.max(homogeneous_rank);
        max_augmented_rank = max_augmented_rank.max(augmented_rank);
        max_formal_coefficient_rank = max_formal_coefficient_rank.max(formal_coefficient_rank);
        max_formal_homogeneous_rank = max_formal_homogeneous_rank.max(formal_homogeneous_rank);
        max_formal_augmented_rank = max_formal_augmented_rank.max(formal_augmented_rank);
    }
    let empty_targets = 75usize
        .checked_sub(buckets.len())
        .ok_or("more than 75 canonical subgroup targets")?;
    if empty_targets != 0 {
        occupancy_histogram.insert(0, empty_targets as u64);
    }
    Ok(BaseMetrics {
        raw_rows,
        distinct_normalized_rows: distinct_rows,
        nonempty_target_buckets: buckets.len() as u64,
        maximum_distinct_occupancy: max_occupancy as u64,
        maximum_occupancy_target: max_target as u64,
        maximum_known_target_coefficient_rank_mod_149: max_coefficient_rank as u64,
        maximum_homogeneous_rank_mod_149: max_homogeneous_rank as u64,
        maximum_unknown_target_augmented_rank_mod_149: max_augmented_rank as u64,
        maximum_formal_coefficient_rank_mod_2_61_minus_1: max_formal_coefficient_rank as u64,
        maximum_formal_homogeneous_rank_mod_2_61_minus_1: max_formal_homogeneous_rank as u64,
        maximum_formal_augmented_rank_mod_2_61_minus_1: max_formal_augmented_rank as u64,
        nonzero_buckets_with_full_known_target_rank: full_rank_nonzero,
        buckets_with_homogeneous_rank_b_minus_1: homogeneous_b_minus_1,
        buckets_with_formal_rank_divergence: formal_rank_divergences,
        subgroup_kernel_failures: kernel_failures,
        subgroup_rank_ceiling_violations: ceiling_violations,
        replay_failures,
        false_positives: replay_failures,
        false_negatives: expected_raw.saturating_sub(raw_rows),
        scalar_label_additions: raw_rows * arity as u64,
        logical_group_additions: raw_rows * arity as u64,
        rank_field_operations,
        occupancy_histogram_all_75_targets: histogram_rows(occupancy_histogram),
        coefficient_rank_histogram_nonempty_targets: histogram_rows(coefficient_histogram),
        homogeneous_rank_histogram_nonempty_targets: histogram_rows(homogeneous_histogram),
        best_bucket_rows: max_rows,
    })
}

fn better_base(
    candidate: &[u16],
    metrics: &BaseMetrics,
    current: &Option<(Vec<u16>, BaseMetrics)>,
) -> bool {
    let candidate_key = (
        metrics.maximum_known_target_coefficient_rank_mod_149 == PRIMARY_COLUMNS as u64,
        metrics.maximum_distinct_occupancy,
        metrics.maximum_homogeneous_rank_mod_149,
    );
    match current {
        None => true,
        Some((base, current_metrics)) => {
            let current_key = (
                current_metrics.maximum_known_target_coefficient_rank_mod_149
                    == PRIMARY_COLUMNS as u64,
                current_metrics.maximum_distinct_occupancy,
                current_metrics.maximum_homogeneous_rank_mod_149,
            );
            candidate_key > current_key || (candidate_key == current_key && candidate < base)
        }
    }
}

fn base_report(base: Vec<u16>, metrics: BaseMetrics, subgroup: &PrimeSubgroup) -> BaseReport {
    BaseReport {
        canonical_points: base
            .iter()
            .map(|&label| point_label(subgroup.canonical_points[label as usize - 1]))
            .collect(),
        folded_scalar_labels: base.iter().map(|&label| label as u64).collect(),
        metrics,
    }
}

fn update_hasher(hasher: &mut Hasher, base: &[u16], metrics: &BaseMetrics) {
    for &label in base {
        hasher.update(&label.to_be_bytes());
    }
    for value in [
        metrics.raw_rows,
        metrics.distinct_normalized_rows,
        metrics.nonempty_target_buckets,
        metrics.maximum_distinct_occupancy,
        metrics.maximum_occupancy_target,
        metrics.maximum_known_target_coefficient_rank_mod_149,
        metrics.maximum_homogeneous_rank_mod_149,
        metrics.maximum_unknown_target_augmented_rank_mod_149,
        metrics.buckets_with_formal_rank_divergence,
        metrics.replay_failures,
        metrics.rank_field_operations.additions_or_subtractions,
        metrics.rank_field_operations.multiplications,
        metrics.rank_field_operations.inversions,
    ] {
        hasher.update(&value.to_be_bytes());
    }
}

fn evaluate_primary(
    families: &BTreeMap<String, FamilyBuild>,
    subgroup: &PrimeSubgroup,
) -> Result<
    (
        Vec<FamilySummary>,
        u64,
        u64,
        u64,
        u64,
        FieldOperations,
        String,
    ),
    String,
> {
    let mut global = BTreeMap::<Vec<u16>, BTreeSet<String>>::new();
    for (name, family) in families {
        for base in &family.bases {
            global.entry(base.clone()).or_default().insert(name.clone());
        }
    }
    let mut aggregates = families
        .keys()
        .map(|name| (name.clone(), FamilyAggregate::default()))
        .collect::<BTreeMap<_, _>>();
    let mut total_raw = 0u64;
    let mut total_distinct = 0u64;
    let mut total_group_additions = 0u64;
    let mut total_rank_field_operations = FieldOperations::default();
    let mut hasher = Hasher::new();
    for (base, memberships) in &global {
        let metrics = evaluate_base(base, PRIMARY_ARITY, subgroup)?;
        update_hasher(&mut hasher, base, &metrics);
        total_raw += metrics.raw_rows;
        total_distinct += metrics.distinct_normalized_rows;
        total_group_additions += metrics.logical_group_additions;
        total_rank_field_operations.add_assign(metrics.rank_field_operations);
        for name in memberships {
            let aggregate = aggregates.get_mut(name).expect("registered family");
            aggregate.raw_rows += metrics.raw_rows;
            aggregate.distinct_rows += metrics.distinct_normalized_rows;
            aggregate.kernel_failures += metrics.subgroup_kernel_failures;
            aggregate.replay_failures += metrics.replay_failures;
            aggregate.false_positives += metrics.false_positives;
            aggregate.false_negatives += metrics.false_negatives;
            aggregate.scalar_label_additions += metrics.scalar_label_additions;
            aggregate
                .rank_field_operations
                .add_assign(metrics.rank_field_operations);
            if metrics.nonzero_buckets_with_full_known_target_rank != 0 {
                aggregate.full_rank += 1;
            }
            if metrics.buckets_with_homogeneous_rank_b_minus_1 != 0 {
                aggregate.homogeneous_b_minus_1 += 1;
            }
            if metrics.subgroup_rank_ceiling_violations != 0 {
                aggregate.ceiling_violations += 1;
            }
            if metrics.buckets_with_formal_rank_divergence != 0 {
                aggregate.formal_rank_divergences += 1;
            }
            aggregate.maximum_occupancy = aggregate
                .maximum_occupancy
                .max(metrics.maximum_distinct_occupancy);
            *aggregate
                .occupancy_histogram
                .entry(metrics.maximum_distinct_occupancy as usize)
                .or_default() += 1;
            *aggregate
                .rank_histogram
                .entry(metrics.maximum_known_target_coefficient_rank_mod_149 as usize)
                .or_default() += 1;
            if better_base(base, &metrics, &aggregate.best) {
                aggregate.best = Some((base.clone(), metrics.clone()));
            }
        }
    }
    let mut summaries = Vec::new();
    for (name, family) in families {
        let aggregate = aggregates.remove(name).expect("family aggregate");
        summaries.push(FamilySummary {
            name: name.clone(),
            selection_uses_discrete_log_labels: family.scalar_defined,
            proposals: family.proposals,
            admissible_proposals: family.admissible,
            unique_bases: family.bases.len() as u64,
            raw_rows_replayed: aggregate.raw_rows,
            distinct_rows_checked: aggregate.distinct_rows,
            bases_with_full_known_target_rank: aggregate.full_rank,
            bases_with_homogeneous_rank_b_minus_1: aggregate.homogeneous_b_minus_1,
            bases_with_subgroup_ceiling_violation: aggregate.ceiling_violations,
            bases_with_formal_rank_divergence: aggregate.formal_rank_divergences,
            subgroup_kernel_failures: aggregate.kernel_failures,
            replay_failures: aggregate.replay_failures,
            false_positives: aggregate.false_positives,
            false_negatives: aggregate.false_negatives,
            scalar_label_additions: aggregate.scalar_label_additions,
            rank_field_operations: aggregate.rank_field_operations,
            maximum_distinct_occupancy: aggregate.maximum_occupancy,
            max_occupancy_histogram: histogram_rows(aggregate.occupancy_histogram),
            max_known_rank_histogram: histogram_rows(aggregate.rank_histogram),
            best_base: aggregate
                .best
                .map(|(base, metrics)| base_report(base, metrics, subgroup)),
        });
    }
    Ok((
        summaries,
        global.len() as u64,
        total_raw,
        total_distinct,
        total_group_additions,
        total_rank_field_operations,
        hasher.finalize().to_hex().to_string(),
    ))
}

fn exhaustive_reference(subgroup: &PrimeSubgroup) -> Result<ReferenceSummary, String> {
    let expected_bases = binomial_u64(74, REFERENCE_COLUMNS);
    let mut combination = (0..REFERENCE_COLUMNS).collect::<Vec<_>>();
    let mut bases = 0u64;
    let mut raw_rows = 0u64;
    let mut distinct_rows = 0u64;
    let mut max_occupancy = 0u64;
    let mut max_coefficient_rank = 0u64;
    let mut max_homogeneous_rank = 0u64;
    let mut max_formal_coefficient_rank = 0u64;
    let mut max_formal_homogeneous_rank = 0u64;
    let mut full_rank = 0u64;
    let mut formal_rank_divergences = 0u64;
    let mut kernel_failures = 0u64;
    let mut ceiling_violations = 0u64;
    let mut replay_failures = 0u64;
    let mut false_positives = 0u64;
    let mut false_negatives = 0u64;
    let mut scalar_label_additions = 0u64;
    let mut group_additions = 0u64;
    let mut rank_field_operations = FieldOperations::default();
    let mut hasher = Hasher::new();
    loop {
        let base = combination
            .iter()
            .map(|&index| index as u16 + 1)
            .collect::<Vec<_>>();
        let metrics = evaluate_base(&base, REFERENCE_ARITY, subgroup)?;
        update_hasher(&mut hasher, &base, &metrics);
        bases += 1;
        raw_rows += metrics.raw_rows;
        distinct_rows += metrics.distinct_normalized_rows;
        max_occupancy = max_occupancy.max(metrics.maximum_distinct_occupancy);
        max_coefficient_rank =
            max_coefficient_rank.max(metrics.maximum_known_target_coefficient_rank_mod_149);
        max_homogeneous_rank = max_homogeneous_rank.max(metrics.maximum_homogeneous_rank_mod_149);
        max_formal_coefficient_rank = max_formal_coefficient_rank
            .max(metrics.maximum_formal_coefficient_rank_mod_2_61_minus_1);
        max_formal_homogeneous_rank = max_formal_homogeneous_rank
            .max(metrics.maximum_formal_homogeneous_rank_mod_2_61_minus_1);
        if metrics.nonzero_buckets_with_full_known_target_rank != 0 {
            full_rank += 1;
        }
        if metrics.buckets_with_formal_rank_divergence != 0 {
            formal_rank_divergences += 1;
        }
        kernel_failures += metrics.subgroup_kernel_failures;
        ceiling_violations += metrics.subgroup_rank_ceiling_violations;
        replay_failures += metrics.replay_failures;
        false_positives += metrics.false_positives;
        false_negatives += metrics.false_negatives;
        scalar_label_additions += metrics.scalar_label_additions;
        group_additions += metrics.logical_group_additions;
        rank_field_operations.add_assign(metrics.rank_field_operations);
        if !next_combination(&mut combination, 74) {
            break;
        }
    }
    if bases != expected_bases {
        return Err(format!(
            "B=4 reference enumerated {bases} bases, expected {expected_bases}"
        ));
    }
    Ok(ReferenceSummary {
        columns: REFERENCE_COLUMNS as u64,
        arity: REFERENCE_ARITY as u64,
        all_factor_bases: bases,
        raw_rows_per_base: (1u64 << REFERENCE_ARITY)
            * binomial_u64(REFERENCE_COLUMNS, REFERENCE_ARITY),
        distinct_rows_per_base: ((1u64 << REFERENCE_ARITY)
            * binomial_u64(REFERENCE_COLUMNS, REFERENCE_ARITY))
            / 2,
        raw_rows_replayed: raw_rows,
        distinct_rows_checked: distinct_rows,
        maximum_distinct_occupancy: max_occupancy,
        maximum_known_target_rank_mod_149: max_coefficient_rank,
        maximum_homogeneous_rank_mod_149: max_homogeneous_rank,
        maximum_formal_known_target_rank_mod_2_61_minus_1: max_formal_coefficient_rank,
        maximum_formal_homogeneous_rank_mod_2_61_minus_1: max_formal_homogeneous_rank,
        bases_with_full_known_target_rank: full_rank,
        bases_with_formal_rank_divergence: formal_rank_divergences,
        subgroup_kernel_failures: kernel_failures,
        subgroup_rank_ceiling_violations: ceiling_violations,
        replay_failures,
        false_positives,
        false_negatives,
        scalar_label_additions,
        logical_group_additions: group_additions,
        rank_field_operations,
        enumeration_blake3: hasher.finalize().to_hex().to_string(),
    })
}

fn log2_factorial(value: u64) -> f64 {
    (2..=value).map(|item| (item as f64).log2()).sum()
}

fn gamma_ratio(collisions: u64) -> f64 {
    let mut ratio = std::f64::consts::PI.sqrt() / 2.0;
    for k in 1..collisions {
        ratio *= (k as f64 + 0.5) / k as f64;
    }
    ratio
}

fn ratio_to_rho(events: u64) -> f64 {
    2.0_f64.sqrt() * gamma_ratio(events) / RHO_S
}

fn p256_boundary(round298: &Value) -> Result<P256BoundaryCorrection, String> {
    let quotient = round298
        .pointer("/p256_boundary/exact_global_negation_quotient_domain")
        .and_then(Value::as_str)
        .ok_or("Round 298 lacks normalized P-256 domain")?;
    let targets = round298
        .pointer("/p256_boundary/canonical_target_space")
        .and_then(Value::as_str)
        .ok_or("Round 298 lacks canonical target space")?;
    let domain = BigUint::parse_bytes(quotient.as_bytes(), 10).ok_or("bad normalized domain")?;
    let target_count = BigUint::parse_bytes(targets.as_bytes(), 10).ok_or("bad target space")?;
    let domain_f = domain
        .to_f64()
        .ok_or("normalized domain does not fit f64")?;
    let targets_f = target_count
        .to_f64()
        .ok_or("target space does not fit f64")?;
    let lambda = domain_f / targets_f;
    let occupancy = P256_COLUMNS;
    let log2_targets = targets_f.log2();
    let union_bound =
        log2_targets + occupancy as f64 * (std::f64::consts::E * lambda / occupancy as f64).log2();
    let poisson_expected = log2_targets - lambda / std::f64::consts::LN_2
        + occupancy as f64 * lambda.log2()
        - log2_factorial(occupancy)
        - (1.0 - lambda / (occupancy + 1) as f64).log2();
    Ok(P256BoundaryCorrection {
        columns: P256_COLUMNS,
        prior_round_distinct_rows_required: P256_COLUMNS + 1,
        known_nonzero_target_distinct_rows_required: P256_COLUMNS,
        equivalent_raw_folded_preimages_required: 2 * P256_COLUMNS,
        homogeneous_same_target_rank_ceiling: P256_COLUMNS - 1,
        exact_global_negation_quotient_domain: quotient.into(),
        canonical_target_space: targets.into(),
        mean_distinct_rows_per_canonical_target: lambda,
        random_map_null_log2_union_bound_any_bucket_at_164: union_bound,
        random_map_null_log2_poisson_expected_buckets_at_164: poisson_expected,
        random_map_null_is_measurement: false,
        one_event_ratio_to_rho_before_omitted_costs: ratio_to_rho(1),
        two_event_ratio_to_rho_before_omitted_costs: ratio_to_rho(2),
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

fn build_assessment() -> TransferAssessment {
    let source = format!("<{CURVE_SLUG}, H=<G>, |H|=n prime>");
    let nonzero = TransferCase {
        name: "nonzero homomorphic restriction".into(),
        source: source.clone(),
        destination: "abelian group A containing phi(H)".into(),
        transformation_type: "group homomorphism induced by a typed correspondence".into(),
        field_of_definition: "must be supplied by the proposed correspondence".into(),
        relevant_subgroup: "H and phi(H), both order n".into(),
        obligations: vec![
            obligation("existence and formulas", "unknown", "No new cover/Jacobian formulas were supplied in Round 299.", "construction-specific"),
            obligation("kernel intersection", "supported", "For prime-order H, phi(G)!=0 makes ker(phi) a proper subgroup of H, hence trivial.", "all group-homomorphic restrictions"),
            obligation("subgroup preservation", "supported", "First isomorphism theorem gives |phi(H)|=n and phi([x]G)=[x]phi(G).", "all nonzero restrictions"),
            obligation("logarithm transport", "supported", "log_{phi(G)} phi([x]G)=x; the unknown scalar is preserved rather than quotiented.", "all nonzero restrictions"),
            obligation("measured decomposition advantage", "unknown", "No executable destination decomposition system or paired P-256 cost was supplied.", "construction-specific"),
            obligation("recovery", "unknown", "Injectivity on H proves uniqueness but does not price or implement an inverse map.", "construction-specific"),
        ],
        conclusion: "A nonzero homomorphic lift relocates the same cyclic DLP. It can help only through a separately demonstrated destination decomposition advantage, not through log quotienting alone.".into(),
        exploration_boundary: "Exact algebraic statement for the restriction to the P-256 prime-order subgroup; not a search over cover/Jacobian constructions.".into(),
    };
    let zero = TransferCase {
        name: "zero homomorphic restriction".into(),
        source: source.clone(),
        destination: "abelian group A".into(),
        transformation_type: "group homomorphism with phi(G)=0".into(),
        field_of_definition: "arbitrary".into(),
        relevant_subgroup: "H is contained in the kernel".into(),
        obligations: vec![
            obligation(
                "kernel intersection",
                "supported",
                "phi([x]G)=[x]phi(G)=0 for every x.",
                "all zero restrictions",
            ),
            obligation(
                "recovery",
                "refuted",
                "All n source elements have the same image; recovery on H is impossible.",
                "all zero restrictions",
            ),
        ],
        conclusion: "Killing the subgroup is not useful transport.".into(),
        exploration_boundary: "Complete for group homomorphisms restricted to H.".into(),
    };
    let image_endomorphism = TransferCase {
        name: "endomorphism of transported image".into(),
        source: "phi(H) ~= C_n".into(),
        destination: "phi(H)".into(),
        transformation_type: "group endomorphism".into(),
        field_of_definition: "irrelevant to the abstract cyclic restriction; executable formulas remain construction-specific".into(),
        relevant_subgroup: "phi(H)".into(),
        obligations: vec![
            obligation("action", "supported", "End(C_n)=Z/nZ: the image of phi(G) uniquely determines scalar [a].", "all image-preserving endomorphisms"),
            obligation("non-generic quotient", "refuted", "Known invertible a gives a scalar orbit; a=0 kills the subgroup; neither creates independent-log identifications for a generic geometric base.", "endomorphism action alone"),
        ],
        conclusion: "Post-transport symmetries on the prime cyclic image are scalar actions.".into(),
        exploration_boundary: "Does not classify nonlinear set maps or additional destination relations.".into(),
    };
    let mut assessment = TransferAssessment {
        schema: "p256-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 299,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        methodology_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/references/methodology.md (unavailable)".into(),
        template_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/assets/assessment-template.json (unavailable)".into(),
        typed_correspondence_graph: vec![
            format!("{CURVE_SLUG}/H --[phi|H=0]--> one point: nonrecoverable"),
            format!("{CURVE_SLUG}/H --[phi|H nonzero, injective]--> phi(H)~=C_n: log scalar preserved"),
            "phi(H) --[any image endomorphism]--> phi(H): scalar [a] mod n".into(),
            "destination representation --[separately required decomposition and recovery evidence]--> possible algorithmic claim: open".into(),
        ],
        theorem_cases: vec![zero, nonzero, image_endomorphism],
        controls: vec![
            "Complete order-149 toy subgroup and exact group replay for every registered row.".into(),
            "Scalar-defined high-multiplicity families are explicit generic-transport controls, not promotable geometric bases.".into(),
            "Low-coordinate and hash-coordinate families select points without toy discrete-log labels.".into(),
            "All C(74,4) factor bases at arity two form the exhaustive implementation reference.".into(),
        ],
        cost_accounting: vec![
            "Measured in the companion result: base proposals, unique bases, raw and normalized rows, logical group additions, exact transition construction, wall/user/system time and RSS.".into(),
            "Unset for P-256: cover construction, conversion, unsuccessful attempts, destination membership, relation yield, verification, inverse recovery, memory and sparse linear algebra.".into(),
        ],
        weakest_open_obligation: "No executable non-endomorphism cover/Jacobian correspondence with destination relation system, inverse recovery and complete paired P-256 cost has been supplied.".into(),
        narrowest_supported_finding: "Any nonzero group-homomorphic restriction from the prime-order P-256 subgroup is injective and preserves discrete-log scalars; representation transport alone is not a logarithm quotient.".into(),
        semantic_evidence_sha256: String::new(),
        assessment_json_bytes: 0,
    };
    let semantic = json!({
        "curve": assessment.curve,
        "typed_correspondence_graph": assessment.typed_correspondence_graph,
        "theorem_cases": assessment.theorem_cases,
        "controls": assessment.controls,
        "cost_accounting": assessment.cost_accounting,
        "weakest_open_obligation": assessment.weakest_open_obligation,
        "narrowest_supported_finding": assessment.narrowest_supported_finding,
    });
    assessment.semantic_evidence_sha256 =
        sha256_hex(&serde_json::to_vec(&semantic).expect("assessment semantic JSON"));
    assessment
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

fn build_result(cli: &Cli) -> Result<(ResultReceipt, TransferAssessment), String> {
    let (dependency, round298) = checked_dependency(&cli.round298)?;
    if round298.pointer("/schema").and_then(Value::as_str) != Some("p256-parity-escape-screen/v2")
        || round298.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
    {
        return Err("Round 298 dependency schema or gate mismatch".into());
    }
    let curve = CurveParams::p256();
    if curve.n.to_string()
        != "115792089210356248762697446949407573529996955224135760342422259061068512044369"
    {
        return Err("P-256 subgroup order mismatch".into());
    }
    let boundary = p256_boundary(&round298)?;
    let roots = square_roots();
    let subgroup = PrimeSubgroup::build(&roots)?;
    let families = build_families(&subgroup);
    let (
        family_summaries,
        global_bases,
        primary_raw,
        primary_distinct,
        primary_additions,
        primary_rank_field_operations,
        primary_digest,
    ) = evaluate_primary(&families, &subgroup)?;
    let reference = exhaustive_reference(&subgroup)?;
    let total_primary_expected =
        global_bases * (1u64 << PRIMARY_ARITY) * binomial_u64(PRIMARY_COLUMNS, PRIMARY_ARITY);
    if primary_raw != total_primary_expected || primary_distinct * 2 != primary_raw {
        return Err("primary domain accounting mismatch".into());
    }
    let expected_proposals = BTreeMap::from([
        ("affine-progression", 148u64 * 148),
        ("geometric-progression", 148u64 * 148),
        ("pair-sum-closure", 64_824),
        ("hash-scalar", HASH_SCALAR_BASES as u64),
        ("low-coordinate-geometry", 924),
        ("hash-coordinate-geometry", 924),
    ]);
    let all_families_complete = family_summaries.iter().all(|summary| {
        expected_proposals.get(summary.name.as_str()) == Some(&summary.proposals)
            && summary.unique_bases != 0
            && summary.raw_rows_replayed == summary.unique_bases * 160
            && summary.distinct_rows_checked == summary.unique_bases * 80
    });
    let all_rows_replayed = family_summaries
        .iter()
        .all(|summary| summary.replay_failures == 0)
        && reference.replay_failures == 0;
    let zero_errors = family_summaries
        .iter()
        .all(|summary| summary.false_positives == 0 && summary.false_negatives == 0)
        && reference.false_positives == 0
        && reference.false_negatives == 0;
    let rank_ceilings = family_summaries.iter().all(|summary| {
        summary.bases_with_subgroup_ceiling_violation == 0 && summary.subgroup_kernel_failures == 0
    }) && reference.subgroup_kernel_failures == 0
        && reference.subgroup_rank_ceiling_violations == 0;
    let geometry_candidate = family_summaries.iter().any(|summary| {
        !summary.selection_uses_discrete_log_labels
            && summary.bases_with_full_known_target_rank != 0
    });
    let assessment = build_assessment();
    let gates = Gates {
        dependency_hash_checked: true,
        toy_group_and_prime_subgroup_exact: subgroup.points.len() == SUBGROUP_ORDER as usize,
        all_registered_factor_base_families_complete: all_families_complete,
        exhaustive_b4_reference_complete: reference.all_factor_bases
            == binomial_u64(74, REFERENCE_COLUMNS),
        every_raw_row_group_replayed: all_rows_replayed,
        zero_false_positives_and_false_negatives: zero_errors,
        subgroup_rank_ceilings_hold: rank_ceilings,
        independent_anchor_eliminated_by_homomorphic_transfer: false,
        geometry_defined_sparse_full_rank_candidate_observed: geometry_candidate,
        p256_geometry_construction_and_scaling_measured: false,
        structured_residual_degree_at_most_5: false,
        parity_at_or_below_rho: false,
        per_usable_relation_below_2_103: false,
        complete_collection_below_2_120: false,
        projected_storage_below_2_50: false,
        promoted: false,
    };
    if !all_families_complete || !all_rows_replayed || !zero_errors || !rank_ceilings {
        return Err("Round 299 correctness gate failed".into());
    }
    let classification = if geometry_candidate {
        "toy-geometry-candidate/anchor-conservation-proved/p256-promotion-negative"
    } else {
        "structured-multiplicity-control/anchor-conservation-proved/promotion-negative"
    };
    let dominant_obstruction = "A nonzero homomorphic lift of the prime-order P-256 subgroup is injective and preserves the discrete-log scalar. Same-target homogeneous rows retain a one-dimensional log kernel. Structured multiplicity can supply known-target rank only if the factor-base construction and decomposition search are both available; no P-256 geometry construction, degree<=5 solver, measured yield or complete cost was established.";
    let decision = if geometry_candidate {
        "Retain the finite geometry-defined toy cell as a next-obligation candidate, but grant no P-256 or rho credit. Require cross-curve scaling, an ICV1 P-256 construction independent of discrete-log labels, and complete decomposition/recovery costs. Do not attempt an unplanted P-256 relation."
    } else {
        "Close homomorphic representation transport as a quotient mechanism and treat scalar-defined multiplicity as a generic control. No geometry-defined toy family passes the sparse full-rank trigger. Do not attempt an unplanted P-256 relation."
    };
    let semantic = json!({
        "curve": CURVE_SLUG,
        "factor_base": FACTOR_BASE_ID,
        "dependency": dependency,
        "p256_boundary_correction": boundary,
        "family_summaries": family_summaries,
        "exhaustive_reference": reference,
        "transfer_assessment_semantic_sha256": assessment.semantic_evidence_sha256,
        "gates": gates,
        "classification": classification,
        "dominant_obstruction": dominant_obstruction,
        "decision": decision,
    });
    let semantic_evidence_sha256 =
        sha256_hex(&serde_json::to_vec(&semantic).map_err(|error| error.to_string())?);
    Ok((
        ResultReceipt {
            schema: "p256-transfer-anchor-screen/v1".into(),
            curve: CURVE_SLUG.into(),
            factor_base: FACTOR_BASE_ID.into(),
            screening_round: 299,
            execution_status: "complete".into(),
            dependency,
            p256_boundary_correction: boundary,
            toy_curve: "E/F_1151: y^2=x^3-3x+241".into(),
            toy_group_order: TOY_GROUP_ORDER as u64,
            cofactor: COFACTOR,
            prime_subgroup_order: SUBGROUP_ORDER,
            subgroup_generator: point_label(subgroup.generator),
            subgroup_generator_source_point: point_label(subgroup.source),
            subgroup_points_enumerated: subgroup.points.len() as u64,
            folded_columns: 74,
            primary_columns: PRIMARY_COLUMNS as u64,
            primary_arity: PRIMARY_ARITY as u64,
            raw_rows_per_primary_base: 160,
            distinct_rows_per_primary_base: 80,
            canonical_subgroup_targets: 75,
            fixed_distinct_mean_occupancy: 80.0 / 75.0,
            family_summaries,
            global_unique_primary_bases: global_bases,
            global_primary_raw_rows_replayed: primary_raw,
            global_primary_distinct_rows_checked: primary_distinct,
            global_primary_scalar_label_additions: primary_additions,
            global_primary_logical_group_additions: primary_additions,
            global_primary_rank_field_operations: primary_rank_field_operations,
            transition_table_exact_group_additions: subgroup.exact_transition_additions,
            transition_table_field_operations: subgroup.transition_field_operations,
            primary_enumeration_blake3: primary_digest,
            exhaustive_reference: reference,
            transfer_assessment_semantic_sha256: assessment.semantic_evidence_sha256.clone(),
            relations_reported_on_p256: 0,
            full_depth_unplanted_p256_relation_attempted: false,
            gates,
            classification: classification.into(),
            dominant_obstruction: dominant_obstruction.into(),
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

fn run(cli: Cli) -> Result<(), String> {
    let (mut result, mut assessment) = build_result(&cli)?;
    write_assessment(&cli.assessment, &mut assessment)?;
    write_result(&cli.out, &mut result)?;
    eprintln!(
        "round 299: {} unique B=6 bases, {} primary raw rows, geometry candidate={}, promoted={}",
        result.global_unique_primary_bases,
        result.global_primary_raw_rows_replayed,
        result
            .gates
            .geometry_defined_sparse_full_rank_candidate_observed,
        result.gates.promoted
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_transfer_anchor_screen: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn toy_group_and_prime_subgroup_replay() {
        let roots = square_roots();
        let subgroup = PrimeSubgroup::build(&roots).expect("subgroup");
        assert_eq!(subgroup.points.len(), 149);
        assert_eq!(scalar_mul(subgroup.generator, 149), ToyPoint::Infinity);
        assert_ne!(subgroup.generator, ToyPoint::Infinity);
    }

    #[test]
    fn normalized_primary_domain_halves_raw_domain() {
        let roots = square_roots();
        let subgroup = PrimeSubgroup::build(&roots).expect("subgroup");
        let metrics = evaluate_base(&[1, 2, 3, 4, 5, 6], 3, &subgroup).expect("evaluation");
        assert_eq!(metrics.raw_rows, 160);
        assert_eq!(metrics.distinct_normalized_rows, 80);
        assert_eq!(metrics.replay_failures, 0);
        assert!(metrics.maximum_homogeneous_rank_mod_149 <= 5);
        assert!(metrics.maximum_unknown_target_augmented_rank_mod_149 <= 6);
    }

    #[test]
    fn family_sizes_are_frozen() {
        let roots = square_roots();
        let subgroup = PrimeSubgroup::build(&roots).expect("subgroup");
        let families = build_families(&subgroup);
        assert_eq!(families["hash-scalar"].proposals, 4096);
        assert_eq!(families["low-coordinate-geometry"].proposals, 924);
        assert_eq!(families["hash-coordinate-geometry"].proposals, 924);
        assert_eq!(families["pair-sum-closure"].proposals, 64_824);
    }

    #[test]
    fn rank_distinguishes_anchor_kernel() {
        let rows = vec![vec![1, 1, -1, 0, 0, 0], vec![0, 1, 1, -1, 0, 0]];
        assert_eq!(matrix_rank(&rows, SUBGROUP_ORDER).rank, 2);
        let labels = [1u64, 2, 3, 5, 8, 13];
        assert_eq!(
            rows[0]
                .iter()
                .zip(labels)
                .map(|(&coefficient, label)| mod_value(coefficient, SUBGROUP_ORDER) * label)
                .sum::<u64>()
                % SUBGROUP_ORDER,
            0
        );
    }
}
