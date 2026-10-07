//! Round 298: exact multi-row and typed-transfer parity escape screen for P-256.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, ToPrimitive};
use serde::Serialize;
use serde_json::{Value, json};

const ROUND297_SHA256: &str = "53949881daf4a48769c9ef8f869fad633c3dc69d11b0d341a5b7416de2e43982";
const ROUND296_SHA256: &str = "52cbea224df1e1d4d6e8406f368ba413fd4e508f5309fcde65f67a98bc91c5ba";
const ROUND25_SHA256: &str = "dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf";
const ROUND28_SHA256: &str = "950d4b7f4997e146ee30fc7c61174070595fa5d344b9a30cbc0ca578ae7e60f4";
const P256_COLUMNS: u64 = 164;
const P256_ARITY: u64 = 109;
const P256_DOMAIN: &str =
    "116383775127750444427048016195109203089116454475801022361283226351238392053760";
const RHO_S: f64 = 1.3;
const P: u64 = 1151;
const A: u64 = P - 3;
const B: u64 = 241;
const RANK_MODULUS: u64 = (1u64 << 61) - 1;
const TOY_SIZES: [(usize, usize); 4] = [(3, 5), (5, 8), (7, 11), (9, 14)];

#[derive(Parser)]
#[command(about = "Screen multi-row and transfer escapes from the P-256 rho boundary")]
struct Cli {
    #[arg(long)]
    round297: PathBuf,
    #[arg(long)]
    round296: PathBuf,
    #[arg(long)]
    round25: PathBuf,
    #[arg(long)]
    round28: PathBuf,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Clone, Serialize)]
struct Dependency {
    path: String,
    bytes: u64,
    sha256: String,
}

#[derive(Serialize)]
struct P256Boundary {
    columns: u64,
    arity: u64,
    exact_signed_domain: String,
    subgroup_order: String,
    mean_target_occupancy: f64,
    conditional_mean_occupancy_given_nonempty_random_map: f64,
    parity_collision_events: u64,
    parity_event_ratio_to_rho: f64,
    two_event_ratio_to_rho: f64,
    minimum_independent_rows_per_event: u64,
    minimum_bucket_occupancy_for_collision_rows: u64,
    random_map_null_log2_union_bound_any_bucket_at_required_occupancy: f64,
    random_map_null_log2_poisson_expected_buckets_at_required_occupancy: f64,
    random_map_null_is_measurement: bool,
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

#[derive(Serialize)]
struct HistogramRow {
    value: u64,
    buckets: u64,
}

#[derive(Serialize)]
struct ToyCase {
    arity: usize,
    columns: usize,
    support_x: Vec<u64>,
    exact_domain: u64,
    candidates_enumerated: u64,
    target_buckets: u64,
    nonidentity_target_buckets: u64,
    identity_bucket_occupancy: u64,
    mean_bucket_occupancy: f64,
    maximum_bucket_occupancy: u64,
    maximum_bucket_target: String,
    maximum_augmented_rank: u64,
    maximum_homogeneous_rank: u64,
    maximum_homogeneous_rank_target: String,
    buckets_with_full_column_rank: u64,
    occupancy_histogram: Vec<HistogramRow>,
    augmented_rank_histogram: Vec<HistogramRow>,
    homogeneous_rank_histogram: Vec<HistogramRow>,
    rank_modulus: u64,
    hadamard_log2_bound: f64,
    rank_transfer_to_p256_exact: bool,
    logical_candidate_group_additions: u64,
    precomputed_group_additions: u64,
    replayed_rows: u64,
    replay_failures: u64,
    enumeration_sha256: String,
    summary_sha256: String,
}

#[derive(Serialize)]
struct Obligation {
    name: String,
    status: String,
    evidence: String,
}

#[derive(Serialize)]
struct TransferPath {
    name: String,
    source: String,
    destination: String,
    transformation_type: String,
    field_of_definition: String,
    subgroup_order: String,
    obligations: Vec<Obligation>,
    quotient_dimension: Option<u64>,
    measured_ratio_to_rho: Option<f64>,
    classification: String,
    exploration_boundary: String,
}

#[derive(Serialize)]
struct TransferAssessment {
    skill_profile: String,
    required_skill_resources_available: bool,
    typed_correspondence_graph: Vec<String>,
    paths: Vec<TransferPath>,
    weakest_open_obligation: String,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_checked: bool,
    complete_toy_domains: bool,
    zero_toy_replay_failures: bool,
    exact_rank_transfer: bool,
    p256_rank_164_event_measured: bool,
    exact_non_generic_log_transport: bool,
    parity_at_or_below_rho: bool,
    zero_false_positives_and_false_negatives_on_complete_relation_instances: bool,
    structured_residual_degree_at_most_5: bool,
    per_usable_relation_below_2_103: bool,
    complete_cost_below_2_120: bool,
    projected_materialized_storage_below_2_50: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultReceipt {
    schema: String,
    curve: String,
    screening_round: u64,
    execution_status: String,
    dependencies: Vec<Dependency>,
    p256_boundary: P256Boundary,
    toy_curve: String,
    toy_support_pool: Vec<u64>,
    toy_support_pool_sha256: String,
    toy_cases: Vec<ToyCase>,
    transfer_assessment: TransferAssessment,
    relations_reported: u64,
    full_depth_unplanted_relation_attempted: bool,
    process_telemetry_recorded_in_isolation_receipt: bool,
    gates: Gates,
    classification: String,
    dominant_obstruction: String,
    decision: String,
    semantic_evidence_sha256: String,
    result_json_bytes: u64,
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
    if value == 0 { 0 } else { P - value }
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

fn toy_iterate(mut value: u64, depth: u32) -> u64 {
    for _ in 0..depth {
        value = subm(mulm(value, value), 2);
    }
    value
}

fn toy_support_pool(roots: &[Vec<u64>]) -> Vec<u64> {
    let mut values = (0..P)
        .filter(|&x| matches!(toy_iterate(x, 5), 0 | 369) && !roots[toy_rhs(x) as usize].is_empty())
        .collect::<Vec<_>>();
    values.sort_by(|left, right| {
        let left_hash = sha256(format!("rr-selector-round296/support/{left}").as_bytes());
        let right_hash = sha256(format!("rr-selector-round296/support/{right}").as_bytes());
        left_hash.cmp(&right_hash).then(left.cmp(right))
    });
    values
}

fn add_points(left: ToyPoint, right: ToyPoint) -> ToyPoint {
    let (left, right) = match (left, right) {
        (ToyPoint::Infinity, value) | (value, ToyPoint::Infinity) => return value,
        (ToyPoint::Affine(left), ToyPoint::Affine(right)) => (left, right),
    };
    let slope = if left.x == right.x {
        if addm(left.y, right.y) == 0 || left.y == 0 {
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

fn negate(point: ToyPoint) -> ToyPoint {
    match point {
        ToyPoint::Infinity => ToyPoint::Infinity,
        ToyPoint::Affine(point) => ToyPoint::Affine(Affine {
            x: point.x,
            y: negm(point.y),
        }),
    }
}

fn low_point(x: u64, roots: &[Vec<u64>]) -> Result<ToyPoint, String> {
    let y = roots[toy_rhs(x) as usize]
        .iter()
        .copied()
        .min()
        .ok_or("support x is not liftable")?;
    Ok(ToyPoint::Affine(Affine { x, y }))
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

struct ToyTransitions {
    points: Vec<ToyPoint>,
    index: BTreeMap<ToyPoint, usize>,
    signed_base: Vec<ToyPoint>,
    transitions: Vec<Vec<usize>>,
    precomputed_additions: u64,
}

impl ToyTransitions {
    fn new(support: &[u64], roots: &[Vec<u64>]) -> Result<Self, String> {
        let points = all_curve_points(roots);
        let index = points
            .iter()
            .copied()
            .enumerate()
            .map(|(position, point)| (point, position))
            .collect::<BTreeMap<_, _>>();
        let mut signed_base = Vec::with_capacity(2 * support.len());
        for &x in support {
            let positive = low_point(x, roots)?;
            signed_base.push(positive);
            signed_base.push(negate(positive));
        }
        let mut transitions = Vec::with_capacity(signed_base.len());
        for &base in &signed_base {
            let mut row = Vec::with_capacity(points.len());
            for &point in &points {
                let sum = add_points(point, base);
                row.push(
                    *index
                        .get(&sum)
                        .ok_or("toy addition left the enumerated curve")?,
                );
            }
            transitions.push(row);
        }
        let precomputed_additions = (signed_base.len() * points.len()) as u64;
        Ok(Self {
            points,
            index,
            signed_base,
            transitions,
            precomputed_additions,
        })
    }

    fn sum_row(&self, row: &[i8]) -> Result<ToyPoint, String> {
        let mut state = *self
            .index
            .get(&ToyPoint::Infinity)
            .ok_or("missing toy identity")?;
        for (column, &coefficient) in row.iter().enumerate() {
            let signed = match coefficient {
                0 => continue,
                1 => 2 * column,
                -1 => 2 * column + 1,
                _ => return Err("toy row coefficient is outside {-1,0,1}".into()),
            };
            state = self.transitions[signed][state];
        }
        Ok(self.points[state])
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

fn binomial_u64(n: usize, k: usize) -> u64 {
    let k = k.min(n - k);
    let mut value = 1u64;
    for i in 0..k {
        value = value * (n - i) as u64 / (i + 1) as u64;
    }
    value
}

fn binomial_big(n: u64, k: u64) -> BigUint {
    let k = k.min(n - k);
    let mut value = BigUint::one();
    for i in 0..k {
        value *= BigUint::from(n - i);
        value /= BigUint::from(i + 1);
    }
    value
}

fn rank_add(left: u64, right: u64) -> u64 {
    ((left as u128 + right as u128) % RANK_MODULUS as u128) as u64
}

fn rank_sub(left: u64, right: u64) -> u64 {
    rank_add(left, if right == 0 { 0 } else { RANK_MODULUS - right })
}

fn rank_mul(left: u64, right: u64) -> u64 {
    ((left as u128 * right as u128) % RANK_MODULUS as u128) as u64
}

fn rank_pow(mut base: u64, mut exponent: u64) -> u64 {
    let mut value = 1u64;
    while exponent != 0 {
        if exponent & 1 == 1 {
            value = rank_mul(value, base);
        }
        base = rank_mul(base, base);
        exponent >>= 1;
    }
    value
}

fn rank_value(value: i16) -> u64 {
    if value >= 0 {
        value as u64
    } else {
        RANK_MODULUS - (-value) as u64
    }
}

fn matrix_rank(rows: &[Vec<i16>]) -> usize {
    if rows.is_empty() {
        return 0;
    }
    let columns = rows[0].len();
    let mut matrix = rows
        .iter()
        .map(|row| {
            row.iter()
                .map(|&value| rank_value(value))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut rank = 0usize;
    for column in 0..columns {
        let Some(pivot) = (rank..matrix.len()).find(|&row| matrix[row][column] != 0) else {
            continue;
        };
        matrix.swap(rank, pivot);
        let inverse = rank_pow(matrix[rank][column], RANK_MODULUS - 2);
        for value in &mut matrix[rank][column..] {
            *value = rank_mul(*value, inverse);
        }
        let pivot_row = matrix[rank].clone();
        for row in rank + 1..matrix.len() {
            let factor = matrix[row][column];
            if factor == 0 {
                continue;
            }
            for next in column..columns {
                matrix[row][next] = rank_sub(matrix[row][next], rank_mul(factor, pivot_row[next]));
            }
        }
        rank += 1;
        if rank == matrix.len() {
            break;
        }
    }
    rank
}

fn point_label(point: ToyPoint) -> String {
    match point {
        ToyPoint::Infinity => "infinity".into(),
        ToyPoint::Affine(point) => format!("{},{}", point.x, point.y),
    }
}

fn histogram_rows(histogram: BTreeMap<usize, u64>) -> Vec<HistogramRow> {
    histogram
        .into_iter()
        .map(|(value, buckets)| HistogramRow {
            value: value as u64,
            buckets,
        })
        .collect()
}

fn toy_case(support: &[u64], m: usize, roots: &[Vec<u64>]) -> Result<ToyCase, String> {
    let columns = support.len();
    let transitions = ToyTransitions::new(support, roots)?;
    if transitions.signed_base.len() != 2 * columns {
        return Err("signed toy base width mismatch".into());
    }
    let expected = (1u64 << m) * binomial_u64(columns, m);
    let mut buckets = BTreeMap::<ToyPoint, Vec<Vec<i8>>>::new();
    let mut enumeration_stream = Vec::with_capacity(expected as usize * (m + 8));
    let mut combination = (0..m).collect::<Vec<_>>();
    let mut candidates = 0u64;
    loop {
        for mask in 0..1u64 << m {
            let mut row = vec![0i8; columns];
            let mut state = *transitions
                .index
                .get(&ToyPoint::Infinity)
                .ok_or("missing toy identity")?;
            for (position, &column) in combination.iter().enumerate() {
                let negative = mask >> position & 1 == 1;
                row[column] = if negative { -1 } else { 1 };
                state = transitions.transitions[2 * column + usize::from(negative)][state];
            }
            let sum = transitions.points[state];
            let negative_sum = negate(sum);
            let target = sum.min(negative_sum);
            if target != sum {
                for coefficient in &mut row {
                    *coefficient = -*coefficient;
                }
            }
            enumeration_stream.extend((candidates as u32).to_be_bytes());
            match target {
                ToyPoint::Infinity => enumeration_stream.push(0),
                ToyPoint::Affine(point) => {
                    enumeration_stream.push(1);
                    enumeration_stream.extend(point.x.to_be_bytes());
                    enumeration_stream.extend(point.y.to_be_bytes());
                }
            }
            enumeration_stream.extend(row.iter().map(|value| *value as u8));
            buckets.entry(target).or_default().push(row);
            candidates += 1;
        }
        if !next_combination(&mut combination, columns) {
            break;
        }
    }
    if candidates != expected
        || buckets.values().map(|rows| rows.len() as u64).sum::<u64>() != expected
    {
        return Err("toy domain or occupancy accounting mismatch".into());
    }

    let mut occupancy_histogram = BTreeMap::<usize, u64>::new();
    let mut augmented_histogram = BTreeMap::<usize, u64>::new();
    let mut homogeneous_histogram = BTreeMap::<usize, u64>::new();
    let mut maximum_occupancy = 0usize;
    let mut maximum_occupancy_target = ToyPoint::Infinity;
    let mut maximum_augmented_rank = 0usize;
    let mut maximum_homogeneous_rank = 0usize;
    let mut maximum_homogeneous_target = ToyPoint::Infinity;
    let mut full_rank_buckets = 0u64;
    let mut replayed = 0u64;
    let mut replay_failures = 0u64;
    let mut identity_occupancy = 0u64;
    let mut summary_stream = Vec::new();
    for (&target, rows) in &buckets {
        if target == ToyPoint::Infinity {
            identity_occupancy = rows.len() as u64;
        }
        for row in rows {
            replayed += 1;
            if transitions.sum_row(row)? != target {
                replay_failures += 1;
            }
        }
        let augmented = rows
            .iter()
            .map(|row| {
                let mut values = row.iter().map(|&value| value as i16).collect::<Vec<_>>();
                values.push(if target == ToyPoint::Infinity { 0 } else { -1 });
                values
            })
            .collect::<Vec<_>>();
        let homogeneous = if target == ToyPoint::Infinity {
            rows.iter()
                .map(|row| row.iter().map(|&value| value as i16).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        } else {
            let first = &rows[0];
            rows.iter()
                .skip(1)
                .map(|row| {
                    row.iter()
                        .zip(first)
                        .map(|(&value, &base)| value as i16 - base as i16)
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        };
        let augmented_rank = matrix_rank(&augmented);
        let homogeneous_rank = matrix_rank(&homogeneous);
        *occupancy_histogram.entry(rows.len()).or_default() += 1;
        *augmented_histogram.entry(augmented_rank).or_default() += 1;
        *homogeneous_histogram.entry(homogeneous_rank).or_default() += 1;
        if rows.len() > maximum_occupancy {
            maximum_occupancy = rows.len();
            maximum_occupancy_target = target;
        }
        maximum_augmented_rank = maximum_augmented_rank.max(augmented_rank);
        if homogeneous_rank > maximum_homogeneous_rank {
            maximum_homogeneous_rank = homogeneous_rank;
            maximum_homogeneous_target = target;
        }
        if homogeneous_rank == columns {
            full_rank_buckets += 1;
        }
        summary_stream.extend(point_label(target).as_bytes());
        summary_stream.push(0);
        summary_stream.extend((rows.len() as u64).to_be_bytes());
        summary_stream.extend((augmented_rank as u64).to_be_bytes());
        summary_stream.extend((homogeneous_rank as u64).to_be_bytes());
    }
    if replay_failures != 0 {
        return Err(format!(
            "toy case m={m}, B={columns} has {replay_failures} replay failures"
        ));
    }
    let augmented_bound = (columns + 1) as f64 / 2.0 * ((m + 1) as f64).log2();
    let difference_bound = columns as f64 * (1.0 + 0.5 * (columns as f64).log2());
    let hadamard_bound = augmented_bound.max(difference_bound);
    if hadamard_bound >= (RANK_MODULUS as f64).log2() {
        return Err("Hadamard rank-transfer bound reaches the rank modulus".into());
    }
    Ok(ToyCase {
        arity: m,
        columns,
        support_x: support.to_vec(),
        exact_domain: expected,
        candidates_enumerated: candidates,
        target_buckets: buckets.len() as u64,
        nonidentity_target_buckets: buckets
            .keys()
            .filter(|&&point| point != ToyPoint::Infinity)
            .count() as u64,
        identity_bucket_occupancy: identity_occupancy,
        mean_bucket_occupancy: candidates as f64 / buckets.len() as f64,
        maximum_bucket_occupancy: maximum_occupancy as u64,
        maximum_bucket_target: point_label(maximum_occupancy_target),
        maximum_augmented_rank: maximum_augmented_rank as u64,
        maximum_homogeneous_rank: maximum_homogeneous_rank as u64,
        maximum_homogeneous_rank_target: point_label(maximum_homogeneous_target),
        buckets_with_full_column_rank: full_rank_buckets,
        occupancy_histogram: histogram_rows(occupancy_histogram),
        augmented_rank_histogram: histogram_rows(augmented_histogram),
        homogeneous_rank_histogram: histogram_rows(homogeneous_histogram),
        rank_modulus: RANK_MODULUS,
        hadamard_log2_bound: hadamard_bound,
        rank_transfer_to_p256_exact: true,
        logical_candidate_group_additions: candidates * m as u64,
        precomputed_group_additions: transitions.precomputed_additions,
        replayed_rows: replayed,
        replay_failures,
        enumeration_sha256: sha256_hex(&enumeration_stream),
        summary_sha256: sha256_hex(&summary_stream),
    })
}

fn gamma_ratio(collisions: u64) -> f64 {
    assert!(collisions >= 1);
    let mut ratio = std::f64::consts::PI.sqrt() / 2.0;
    for k in 1..collisions {
        ratio *= (k as f64 + 0.5) / k as f64;
    }
    ratio
}

fn ratio_to_rho(events: u64) -> f64 {
    2.0_f64.sqrt() * gamma_ratio(events) / RHO_S
}

fn log2_factorial(value: u64) -> f64 {
    (2..=value).map(|item| (item as f64).log2()).sum()
}

fn p256_boundary() -> Result<P256Boundary, String> {
    let curve = CurveParams::p256();
    let domain = (BigUint::one() << P256_ARITY as usize) * binomial_big(P256_COLUMNS, P256_ARITY);
    let frozen = BigUint::parse_bytes(P256_DOMAIN.as_bytes(), 10).ok_or("bad frozen domain")?;
    if domain != frozen {
        return Err("P-256 signed domain does not match the frozen exact value".into());
    }
    let domain_f = domain.to_f64().ok_or("P-256 domain does not fit f64")?;
    let order_f = curve.n.to_f64().ok_or("P-256 order does not fit f64")?;
    let lambda = domain_f / order_f;
    let conditional = lambda / (1.0 - (-lambda).exp());
    let parity_events = (1..=P256_COLUMNS)
        .find(|&events| ratio_to_rho(events) <= 1.0)
        .ok_or("no collision-event count reaches rho")?;
    let rows_per_event = P256_COLUMNS.div_ceil(parity_events);
    let occupancy = rows_per_event + 1;
    let log2_order = order_f.log2();
    let union_bound =
        log2_order + occupancy as f64 * (std::f64::consts::E * lambda / occupancy as f64).log2();
    let poisson_exact = log2_order - lambda / std::f64::consts::LN_2
        + occupancy as f64 * lambda.log2()
        - log2_factorial(occupancy)
        - (1.0 - lambda / (occupancy + 1) as f64).log2();
    Ok(P256Boundary {
        columns: P256_COLUMNS,
        arity: P256_ARITY,
        exact_signed_domain: domain.to_string(),
        subgroup_order: curve.n.to_string(),
        mean_target_occupancy: lambda,
        conditional_mean_occupancy_given_nonempty_random_map: conditional,
        parity_collision_events: parity_events,
        parity_event_ratio_to_rho: ratio_to_rho(parity_events),
        two_event_ratio_to_rho: ratio_to_rho(2),
        minimum_independent_rows_per_event: rows_per_event,
        minimum_bucket_occupancy_for_collision_rows: occupancy,
        random_map_null_log2_union_bound_any_bucket_at_required_occupancy: union_bound,
        random_map_null_log2_poisson_expected_buckets_at_required_occupancy: poisson_exact,
        random_map_null_is_measurement: false,
    })
}

fn require_bool(value: &Value, pointer: &str, expected: bool) -> Result<(), String> {
    let actual = value
        .pointer(pointer)
        .and_then(Value::as_bool)
        .ok_or_else(|| format!("missing Boolean dependency field {pointer}"))?;
    if actual != expected {
        return Err(format!(
            "dependency field {pointer} is {actual}, expected {expected}"
        ));
    }
    Ok(())
}

fn require_u64(value: &Value, pointer: &str, expected: u64) -> Result<(), String> {
    let actual = value
        .pointer(pointer)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("missing integer dependency field {pointer}"))?;
    if actual != expected {
        return Err(format!(
            "dependency field {pointer} is {actual}, expected {expected}"
        ));
    }
    Ok(())
}

fn obligation(name: &str, status: &str, evidence: &str) -> Obligation {
    Obligation {
        name: name.into(),
        status: status.into(),
        evidence: evidence.into(),
    }
}

fn transfer_assessment(round25: &Value, round28: &Value) -> Result<TransferAssessment, String> {
    require_bool(round25, "/cm_screen/factorization_product_verified", true)?;
    require_bool(round25, "/cm_screen/fundamental_discriminant", true)?;
    require_bool(
        round25,
        "/cm_screen/low_degree_noninteger_transport_exists",
        false,
    )?;
    require_u64(
        round25,
        "/factor_base/independent_log_quotient_dimension",
        1,
    )?;
    require_u64(round25, "/factor_base/transport_replay_failures", 0)?;
    require_bool(round28, "/certificate_audit/all_checks_passed", true)?;
    require_u64(round28, "/certificate_audit/degree", 11)?;
    require_u64(round28, "/certificate_audit/edges", 4096)?;
    require_u64(round28, "/certificate_audit/models", 4097)?;
    require_bool(
        round28,
        "/full_depth_screen_exhaustive_over_all_models",
        false,
    )?;
    require_bool(
        round28,
        "/discarded_prefix_branches_counted_as_exhaustive",
        false,
    )?;
    require_bool(
        round28,
        "/promotion_gates/exact_non_generic_log_transport",
        false,
    )?;
    let curve = CurveParams::p256();
    if curve.n.gcd(&BigUint::from(11u8)) != BigUint::one() {
        return Err("degree-11 isogeny kernel can intersect the P-256 subgroup".into());
    }
    let scalar_path = TransferPath {
        name: "base-field endomorphism action".into(),
        source: CURVE_SLUG.into(),
        destination: CURVE_SLUG.into(),
        transformation_type: "endomorphism inducing a prime-subgroup homomorphism".into(),
        field_of_definition: "F_p".into(),
        subgroup_order: curve.n.to_string(),
        obligations: vec![
            obligation("endomorphism-ring construction", "supported", "Round 25 verifies the fundamental CM discriminant, conductor one, and End(E)=Z[pi]."),
            obligation("rational-point action", "supported", "Every a+b*pi acts on E(F_p) as the scalar [a+b]."),
            obligation("non-scalar transport", "refuted", "Within base-field endomorphisms, the rational P-256 subgroup action is always scalar; the minimum noninteger degree is about 2^255.975."),
            obligation("exact quotient", "supported", "The 139,592-column scalar orbit has quotient dimension one and 139,592 replayed transport edges."),
            obligation("complete parity advantage", "refuted", "Known-anchor transport is a generic claw/rho encoding; unknown-anchor elimination needs two events and its free-oracle boundary is 1.446505 times rho."),
        ],
        quotient_dimension: Some(1),
        measured_ratio_to_rho: Some(1.4465047156952717),
        classification: "exact transport, but scalar/generic and above rho after anchor elimination".into(),
        exploration_boundary: "Complete for base-field P-256 endomorphism actions on E(F_p); not a claim about arbitrary covers or higher-dimensional correspondences.".into(),
    };
    let isogeny_path = TransferPath {
        name: "degree-11 isogeny walk".into(),
        source: CURVE_SLUG.into(),
        destination: "4,096 certified adjacent models; staged Dickson full-depth subset".into(),
        transformation_type: "separable degree-11 isogeny inducing a subgroup isomorphism".into(),
        field_of_definition: "F_p".into(),
        subgroup_order: curve.n.to_string(),
        obligations: vec![
            obligation("map existence and certificate", "supported", "Round 28 replays 4,096 continuous degree-11 edges over 4,097 distinct generic-j models with zero certificate failures."),
            obligation("subgroup kernel", "supported", "The target subgroup has prime order n and gcd(n,11)=1, so a degree-11 kernel intersects it trivially."),
            obligation("DLP transport", "supported", "The isogeny restricts to an injective map on the order-n subgroup; recovery requires the certified dual/path and its conversion costs."),
            obligation("factor-base quotient reduction", "refuted", "Every screened geometric column remains an independent log class; Round 28 credits no quotient reduction."),
            obligation("degree and complete-cost advantage", "refuted", "The selected model retains universal degree-four S3 support, worsens the boundary to 494.612 times rho, and fails degree, storage and collection gates."),
            obligation("all-model full-depth search", "unknown", "The walk is exhaustive at depth eight, but prefix selection retains only 64 depth-12 and nine full-depth models; discarded branches are not counted as exhaustive."),
        ],
        quotient_dimension: None,
        measured_ratio_to_rho: Some(494.61185169050134),
        classification: "certified subgroup transfer without logarithm quotient or cost advantage".into(),
        exploration_boundary: "Complete certificate replay for 4,096 edges; staged, non-exhaustive full-depth Dickson inventory as explicitly recorded by Round 28.".into(),
    };
    let open_path = TransferPath {
        name: "extension, cover, or Jacobian correspondence".into(),
        source: CURVE_SLUG.into(),
        destination: "unspecified".into(),
        transformation_type: "unknown until explicit formulas and induced subgroup maps are supplied".into(),
        field_of_definition: "unknown".into(),
        subgroup_order: curve.n.to_string(),
        obligations: vec![
            obligation("executable construction", "unknown", "No new lift, cover, Jacobian map, or recovery formula is supplied by the frozen artifacts."),
            obligation("homomorphism and kernel", "unknown", "Source/destination dimensions, fields, exceptional points, kernel intersection and subgroup image are unset."),
            obligation("paired conversion cost", "unknown", "Construction, conversion, failed attempts, verification, recovery, memory and matched-rho costs are unset."),
            obligation("non-generic quotient", "unknown", "No certified action on factor-base logarithms exists, so no quotient or speed credit is assigned."),
        ],
        quotient_dimension: None,
        measured_ratio_to_rho: None,
        classification: "open obligation; no evidence of advantage".into(),
        exploration_boundary: "Not searched in this round because no concrete correspondence is supplied; this is not a nonexistence result.".into(),
    };
    Ok(TransferAssessment {
        skill_profile: "transfer assessment applied manually; required referenced methodology/template resources were unavailable".into(),
        required_skill_resources_available: false,
        typed_correspondence_graph: vec![
            format!("{CURVE_SLUG} --[a+b*pi on E(F_p) = scalar [a+b]]--> {CURVE_SLUG}"),
            format!("{CURVE_SLUG} --[degree-11 isogeny, kernel disjoint from order-n subgroup]--> certified isogenous model"),
            format!("{CURVE_SLUG} --[no supplied executable map]--> extension/cover/Jacobian: unknown"),
        ],
        paths: vec![scalar_path, isogeny_path, open_path],
        weakest_open_obligation: "No explicit non-endomorphism cover/Jacobian correspondence with subgroup, recovery and full-cost certificates has been supplied; without it, non-generic transport remains unknown and receives no parity credit.".into(),
    })
}

fn build_result(cli: &Cli) -> Result<ResultReceipt, String> {
    let (round297_dependency, round297) = checked_dependency(&cli.round297, ROUND297_SHA256)?;
    let (round296_dependency, round296) = checked_dependency(&cli.round296, ROUND296_SHA256)?;
    let (round25_dependency, round25) = checked_dependency(&cli.round25, ROUND25_SHA256)?;
    let (round28_dependency, round28) = checked_dependency(&cli.round28, ROUND28_SHA256)?;
    require_u64(
        &round297,
        "/stabilizer/independent_log_quotient_dimension",
        P256_COLUMNS,
    )?;
    require_u64(&round297, "/stabilizer/signed_stabilizer_order", 2)?;
    require_bool(&round297, "/gates/scalar_quotient_reduces_boundary", false)?;
    require_bool(
        &round296,
        "/gates/structured_residual_degree_at_most_5",
        false,
    )?;

    let p256_boundary = p256_boundary()?;
    let roots = square_roots();
    let pool = toy_support_pool(&roots);
    if pool.len() < TOY_SIZES.last().expect("nonempty toy ladder").1 {
        return Err("toy support pool is too small".into());
    }
    let mut toy_cases = Vec::new();
    for (m, columns) in TOY_SIZES {
        toy_cases.push(toy_case(&pool[..columns], m, &roots)?);
    }
    let complete_toy_domains = toy_cases
        .iter()
        .all(|case| case.candidates_enumerated == case.exact_domain);
    let zero_replay_failures = toy_cases.iter().all(|case| case.replay_failures == 0);
    let exact_rank_transfer = toy_cases
        .iter()
        .all(|case| case.rank_transfer_to_p256_exact);
    let transfer_assessment = transfer_assessment(&round25, &round28)?;
    let gates = Gates {
        dependency_hashes_checked: true,
        complete_toy_domains,
        zero_toy_replay_failures: zero_replay_failures,
        exact_rank_transfer,
        p256_rank_164_event_measured: false,
        exact_non_generic_log_transport: false,
        parity_at_or_below_rho: false,
        zero_false_positives_and_false_negatives_on_complete_relation_instances: true,
        structured_residual_degree_at_most_5: false,
        per_usable_relation_below_2_103: false,
        complete_cost_below_2_120: false,
        projected_materialized_storage_below_2_50: false,
        promoted: false,
    };
    if !complete_toy_domains || !zero_replay_failures || !exact_rank_transfer {
        return Err("multi-row correctness certificate failed".into());
    }
    let dependencies = vec![
        round297_dependency,
        round296_dependency,
        round25_dependency,
        round28_dependency,
    ];
    let classification = "multi-row-negative/measured-transport-families-closed/promotion-negative";
    let dominant_obstruction = "Rho parity permits only one collision event, so this 164-class base needs one target bucket carrying at least 164 independent rows (165 preimages for collision differences). Its exact fixed-arity domain has mean occupancy only 1.005110 per P-256 target. Complete toy enumeration reaches full row rank only in the dense regime, while certified endomorphism and degree-11 isogeny transports do not supply a non-generic quotient.";
    let decision = "Reject raw multi-relation count as a parity lever without rank-164 event evidence. Close base-field endomorphism and the measured degree-11 isogeny families; retain covers/Jacobians as explicit unknown obligations, not credited improvements. Do not attempt an unplanted P-256 relation.";
    let semantic = json!({
        "curve": CURVE_SLUG,
        "dependencies": dependencies,
        "p256_boundary": p256_boundary,
        "toy_cases": toy_cases,
        "transfer_assessment": transfer_assessment,
        "gates": gates,
        "classification": classification,
        "dominant_obstruction": dominant_obstruction,
        "decision": decision,
    });
    let semantic_evidence_sha256 =
        sha256_hex(&serde_json::to_vec(&semantic).map_err(|error| error.to_string())?);
    Ok(ResultReceipt {
        schema: "p256-parity-escape-screen/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 298,
        execution_status: "complete".into(),
        dependencies,
        p256_boundary,
        toy_curve: "E/F_1151: y^2=x^3-3x+241".into(),
        toy_support_pool_sha256: sha256_hex(
            &pool
                .iter()
                .flat_map(|value| value.to_be_bytes())
                .collect::<Vec<_>>(),
        ),
        toy_support_pool: pool,
        toy_cases,
        transfer_assessment,
        relations_reported: 0,
        full_depth_unplanted_relation_attempted: false,
        process_telemetry_recorded_in_isolation_receipt: true,
        gates,
        classification: classification.into(),
        dominant_obstruction: dominant_obstruction.into(),
        decision: decision.into(),
        semantic_evidence_sha256,
        result_json_bytes: 0,
    })
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
    let mut result = build_result(&cli)?;
    write_result(&cli.out, &mut result)?;
    let largest = result.toy_cases.last().expect("nonempty toy cases");
    eprintln!(
        "round 298: parity needs {} rows/event; largest toy occupancy {}, homogeneous rank {}; transport={} -> {}",
        result.p256_boundary.minimum_independent_rows_per_event,
        largest.maximum_bucket_occupancy,
        largest.maximum_homogeneous_rank,
        result.gates.exact_non_generic_log_transport,
        result.classification
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_parity_escape_screen: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_p256_domain_replays() {
        let domain =
            (BigUint::one() << P256_ARITY as usize) * binomial_big(P256_COLUMNS, P256_ARITY);
        assert_eq!(domain.to_string(), P256_DOMAIN);
    }

    #[test]
    fn parity_requires_one_event() {
        let boundary = p256_boundary().expect("boundary");
        assert_eq!(boundary.parity_collision_events, 1);
        assert_eq!(boundary.minimum_independent_rows_per_event, 164);
        assert_eq!(boundary.minimum_bucket_occupancy_for_collision_rows, 165);
        assert!(boundary.parity_event_ratio_to_rho < 1.0);
        assert!(boundary.two_event_ratio_to_rho > 1.4);
    }

    #[test]
    fn modular_rank_handles_dependencies() {
        let matrix = vec![vec![1, 0, -1], vec![0, 1, -1], vec![1, 1, -2]];
        assert_eq!(matrix_rank(&matrix), 2);
    }

    #[test]
    fn smallest_toy_domain_is_complete_and_replayed() {
        let roots = square_roots();
        let pool = toy_support_pool(&roots);
        let case = toy_case(&pool[..5], 3, &roots).expect("toy case");
        assert_eq!(case.candidates_enumerated, 80);
        assert_eq!(case.replay_failures, 0);
        assert!(case.rank_transfer_to_p256_exact);
    }
}
