//! Independent-start signed-pair batch router for P-256 (round 40).

use std::fs;
use std::hint::black_box;
use std::mem::{size_of, swap};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use clap::Parser;
use crypto_lib::ct_bignum::U256;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::Serialize;
use serde_json::Value;

const CURVE_SLUG: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const COMPARISON_FB: &str = "FB1h2f8621cda105";
const ROUND19_SHA256: &str = "3096540621408e4a48cfa18963ad01da9686d3527ee26776c8b6cf6f45a71114";
const ROUND36_SHA256: &str = "f145231c098da83878e890d340f8ae921598519b0f2626bcc7c4cd8f9446d015";
const ROUND38_SHA256: &str = "4297328eee07822331c68bddd65d8b0155ba1c72c6e0dda044dc97d42aba5cce";
const PAIR_UNIVERSE: u64 = 34_562_148_612;
const COLUMNS: u64 = 131_458;
const REQUESTS_PER_START: u64 = 8;
const RECORD_BYTES: u64 = 16;
const ACCUMULATOR_BYTES: u64 = 16;
const RADIX_BITS: u32 = 12;
const RADIX_BUCKETS: usize = 1 << RADIX_BITS;
const RADIX_PASSES: u64 = 3;
const LOGICAL_BYTES_PER_REQUEST: u64 = 160;
const BASE_RATIO: f64 = 0.998_569_034_150_286;
const LOCAL_RATIO: f64 = 0.964_336_477_130_181;
const MEAN_CAPACITY: f64 = 224.361_249_307_510_55;
const ALLOWED_ADDITIONS: f64 = 0.334_410_508_422_987_64;
const GROUP_ORDER: &str =
    "115792089210356248762697446949407573529996955224135760342422259061068512044369";
const REGISTERED_DEPTHS: [u32; 5] = [12, 14, 16, 18, 20];
const TIMING_DEPTH: u32 = 16;
const TIMING_REPETITIONS: u64 = 5;
const BASELINE_ADDITIONS: u64 = 262_144;
const NATIVE_PAIR_CONTROLS: u64 = 4_096;

#[derive(Parser)]
#[command(about = "Audit independent-start P-256 signed-pair batching")]
struct Cli {
    #[arg(long)]
    round19: PathBuf,
    #[arg(long)]
    round36: PathBuf,
    #[arg(long)]
    round38: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Clone, Serialize)]
struct Dependency {
    round: u64,
    path: String,
    sha256: String,
    schema: String,
}

#[repr(C)]
#[derive(Clone, Copy, Default)]
struct Request {
    pair: u64,
    owner_slot: u64,
}

#[repr(C)]
#[derive(Clone, Copy, Default)]
struct Accumulator {
    checksum: u64,
    mask: u8,
    count: u8,
    padding: [u8; 6],
}

#[derive(Clone, Serialize)]
struct Cell {
    start_bits: u32,
    starts: u64,
    requests: u64,
    distinct_pairs: u64,
    repeated_requests: u64,
    observed_distinct_per_start: f64,
    expected_distinct_pairs: f64,
    expected_distinct_per_start: f64,
    observed_over_expected: f64,
    request_survival_ratio: f64,
    radix_passes: u64,
    route_updates: u64,
    logical_bytes: u64,
    logical_bytes_per_start: u64,
    peak_materialized_bytes: u64,
    incomplete_start_masks: u64,
    checksum_failures: u64,
    negative_controls: u64,
    false_positives: u64,
    false_negatives: u64,
    ordered_request_sha256: String,
    routed_checksum_sha256: String,
}

#[derive(Clone, Serialize)]
struct ToyCell {
    pair_universe: u64,
    starts: u64,
    requests: u64,
    radix_distinct: u64,
    reference_distinct: u64,
    route_failures: u64,
    false_positives: u64,
    false_negatives: u64,
}

#[derive(Clone, Default, Serialize)]
struct GroupOperations {
    scalar_multiplications: u64,
    additions: u64,
    doublings: u64,
}

#[derive(Clone, Serialize)]
struct NativeControls {
    requests_checked: u64,
    decode_failures: u64,
    group_replay_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    operations: GroupOperations,
    control_sha256: String,
}

#[derive(Clone, Serialize)]
struct TimingSample {
    repetition: u64,
    order: String,
    baseline_wall_ns: u64,
    candidate_wall_ns: u64,
    baseline_checksum: u64,
    candidate_checksum: u64,
}

#[derive(Clone, Serialize)]
struct TimingEvidence {
    timing_start_bits: u32,
    timing_starts: u64,
    baseline_additions: u64,
    warmup_completed: bool,
    aa_first_wall_ns: u64,
    aa_second_wall_ns: u64,
    aa_ratio: f64,
    repetitions: Vec<TimingSample>,
    baseline_median_ns_per_addition: f64,
    candidate_median_ns_per_start: f64,
    candidate_median_addition_equivalents_per_start: f64,
    candidate_budget_multiple: f64,
    candidate_projected_ratio_to_rho: f64,
    candidate_logical_throughput_bytes_per_second: f64,
    candidate_checksums_identical: bool,
    complete_ram_route_below_allowed_additions: bool,
    scope: String,
}

#[derive(Clone, Serialize)]
struct ProjectionRow {
    label: String,
    starts: u64,
    starts_log2: f64,
    requests: u64,
    expected_distinct_pairs: f64,
    expected_distinct_per_start: f64,
    pair_construction_ratio_to_rho: f64,
    pair_construction_reaches_rho: bool,
    peak_materialized_bytes: u64,
    peak_materialized_log2_bytes: f64,
    storage_below_2_50: bool,
    logical_routing_bytes: String,
    logical_routing_log2_bytes: f64,
    logical_bytes_per_start: u64,
    routing_addition_headroom_per_start: f64,
    required_routing_bytes_per_addition_time: Option<f64>,
    applicable_tier_measured: bool,
}

#[derive(Clone, Serialize)]
struct FullProjection {
    total_independent_starts_log2: f64,
    batches_at_selected_storage_row_log2: f64,
    complete_logical_routing_bytes_log2: f64,
    useful_relation_rows_requested: u64,
    dickson_relation_collection_operations: Option<String>,
    dickson_cost_per_usable_relation: Option<String>,
    sparse_linear_algebra_operations: Option<String>,
    classification: String,
}

#[derive(Clone, Serialize)]
struct DegreeBoundary {
    comparison_factor_base: String,
    comparison_structured_residual_maximum: u64,
    comparison_degree_gate_passed: bool,
    measured_selector_family: String,
    measured_family_structured_residual_degree: Option<u64>,
    measured_family_degree_gate_passed: bool,
    unsplit_s17_degree_of_regularity: Option<u64>,
}

#[derive(Clone, Serialize)]
struct Gates {
    zero_false_positives_and_false_negatives: bool,
    exact_group_replay: bool,
    structured_residual_degree_at_most_5_for_measured_family: bool,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    measured_complete_time_below_rho_in_projected_tier: bool,
    no_discarded_branch_counted_as_exhaustive: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    measured_family: String,
    comparison_factor_base: String,
    dependencies: Vec<Dependency>,
    frozen_model: Value,
    toy_controls: Vec<ToyCell>,
    depth_ladder: Vec<Cell>,
    native_controls: NativeControls,
    timing: TimingEvidence,
    projection: Vec<ProjectionRow>,
    full_projection: FullProjection,
    degree_boundary: DegreeBoundary,
    semantic_evidence_sha256: String,
    gates: Gates,
    full_depth_unplanted_attempted: bool,
    interpretation: String,
    decision: String,
}

fn json_at<'a>(value: &'a Value, path: &[&str]) -> Result<&'a Value, String> {
    let mut cursor = value;
    for key in path {
        cursor = cursor
            .get(*key)
            .ok_or_else(|| format!("missing JSON path {}", path.join(".")))?;
    }
    Ok(cursor)
}

fn json_str<'a>(value: &'a Value, path: &[&str]) -> Result<&'a str, String> {
    json_at(value, path)?
        .as_str()
        .ok_or_else(|| format!("JSON path {} is not a string", path.join(".")))
}

fn json_u64(value: &Value, path: &[&str]) -> Result<u64, String> {
    json_at(value, path)?
        .as_u64()
        .ok_or_else(|| format!("JSON path {} is not a u64", path.join(".")))
}

fn json_f64(value: &Value, path: &[&str]) -> Result<f64, String> {
    json_at(value, path)?
        .as_f64()
        .ok_or_else(|| format!("JSON path {} is not numeric", path.join(".")))
}

fn load_dependency(
    path: &Path,
    round: u64,
    expected_hash: &str,
    expected_schema: &str,
) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != expected_hash {
        return Err(format!(
            "round-{round} dependency hash mismatch: expected {expected_hash}, got {digest}"
        ));
    }
    let value: Value = serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if json_str(&value, &["schema"])? != expected_schema {
        return Err(format!("round-{round} schema mismatch"));
    }
    Ok((
        Dependency {
            round,
            path: path.display().to_string(),
            sha256: digest,
            schema: expected_schema.into(),
        },
        value,
    ))
}

fn pair_id(start: u64, slot: u8, universe: u64) -> u64 {
    let zone = u64::MAX / universe * universe;
    for retry in 0u32.. {
        let mut input = [0u8; 26];
        input[..13].copy_from_slice(b"icv1-r40-pair");
        input[13..21].copy_from_slice(&start.to_be_bytes());
        input[21] = slot;
        input[22..26].copy_from_slice(&retry.to_be_bytes());
        let digest = sha256(&input);
        let mut word = [0u8; 8];
        word.copy_from_slice(&digest[..8]);
        let candidate = u64::from_be_bytes(word);
        if candidate < zone {
            return candidate % universe;
        }
    }
    unreachable!()
}

fn contribution(pair: u64, slot: u8) -> u64 {
    let mut value = pair ^ 0x9e37_79b9_7f4a_7c15;
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    (value ^ (value >> 31)).rotate_left(u32::from(slot) * 7)
}

fn generate_requests(starts: u64, universe: u64) -> Vec<Request> {
    let mut records = Vec::with_capacity((starts * REQUESTS_PER_START) as usize);
    for start in 0..starts {
        for slot in 0..REQUESTS_PER_START as u8 {
            records.push(Request {
                pair: pair_id(start, slot, universe),
                owner_slot: (start << 3) | u64::from(slot),
            });
        }
    }
    records
}

fn radix_sort(mut source: Vec<Request>) -> Vec<Request> {
    let mut target = vec![Request::default(); source.len()];
    for pass in 0..RADIX_PASSES as u32 {
        let shift = pass * RADIX_BITS;
        let mut counts = vec![0usize; RADIX_BUCKETS];
        for record in &source {
            counts[((record.pair >> shift) as usize) & (RADIX_BUCKETS - 1)] += 1;
        }
        let mut offset = 0usize;
        for count in &mut counts {
            let current = *count;
            *count = offset;
            offset += current;
        }
        for record in &source {
            let digit = ((record.pair >> shift) as usize) & (RADIX_BUCKETS - 1);
            target[counts[digit]] = *record;
            counts[digit] += 1;
        }
        swap(&mut source, &mut target);
    }
    source
}

fn bytes_digest_requests(records: &[Request]) -> String {
    let mut bytes = Vec::with_capacity(records.len() * 16);
    for record in records {
        bytes.extend(record.pair.to_be_bytes());
        bytes.extend(record.owner_slot.to_be_bytes());
    }
    hex::encode(sha256(&bytes))
}

fn bytes_digest_accumulators(accumulators: &[Accumulator]) -> String {
    let mut bytes = Vec::with_capacity(accumulators.len() * 10);
    for accumulator in accumulators {
        bytes.extend(accumulator.checksum.to_be_bytes());
        bytes.push(accumulator.mask);
        bytes.push(accumulator.count);
    }
    hex::encode(sha256(&bytes))
}

struct Routed {
    sorted: Vec<Request>,
    accumulators: Vec<Accumulator>,
    distinct: u64,
}

fn route(starts: u64, universe: u64) -> Routed {
    let sorted = radix_sort(generate_requests(starts, universe));
    let mut accumulators = vec![Accumulator::default(); starts as usize];
    let mut distinct = 0u64;
    let mut previous = None;
    for record in &sorted {
        if previous != Some(record.pair) {
            distinct += 1;
            previous = Some(record.pair);
        }
        let start = (record.owner_slot >> 3) as usize;
        let slot = (record.owner_slot & 7) as u8;
        let accumulator = &mut accumulators[start];
        accumulator.checksum = accumulator
            .checksum
            .wrapping_add(contribution(record.pair, slot));
        accumulator.mask |= 1u8 << slot;
        accumulator.count += 1;
    }
    Routed {
        sorted,
        accumulators,
        distinct,
    }
}

fn expected_accumulator(start: u64, universe: u64) -> Accumulator {
    let mut accumulator = Accumulator::default();
    for slot in 0..REQUESTS_PER_START as u8 {
        accumulator.checksum = accumulator
            .checksum
            .wrapping_add(contribution(pair_id(start, slot, universe), slot));
        accumulator.mask |= 1u8 << slot;
        accumulator.count += 1;
    }
    accumulator
}

fn expected_distinct(starts: u64) -> f64 {
    let requests = (starts * REQUESTS_PER_START) as f64;
    let universe = PAIR_UNIVERSE as f64;
    -universe * (requests * (-1.0 / universe).ln_1p()).exp_m1()
}

fn checked_cell(bits: u32) -> Cell {
    let starts = 1u64 << bits;
    let requests = starts * REQUESTS_PER_START;
    let routed = route(starts, PAIR_UNIVERSE);
    let ordered_request_sha256 = bytes_digest_requests(&routed.sorted);
    let routed_checksum_sha256 = bytes_digest_accumulators(&routed.accumulators);
    let mut incomplete = 0u64;
    let mut checksum_failures = 0u64;
    let mut false_positives = 0u64;
    for start in 0..starts {
        let actual = routed.accumulators[start as usize];
        let expected = expected_accumulator(start, PAIR_UNIVERSE);
        if actual.mask != 0xff || actual.count != 8 {
            incomplete += 1;
        }
        if actual.checksum != expected.checksum {
            checksum_failures += 1;
        }
        let corrupted = actual.checksum.wrapping_add(1);
        if corrupted == expected.checksum {
            false_positives += 1;
        }
    }
    let expected = expected_distinct(starts);
    Cell {
        start_bits: bits,
        starts,
        requests,
        distinct_pairs: routed.distinct,
        repeated_requests: requests - routed.distinct,
        observed_distinct_per_start: routed.distinct as f64 / starts as f64,
        expected_distinct_pairs: expected,
        expected_distinct_per_start: expected / starts as f64,
        observed_over_expected: routed.distinct as f64 / expected,
        request_survival_ratio: routed.distinct as f64 / requests as f64,
        radix_passes: RADIX_PASSES,
        route_updates: requests,
        logical_bytes: requests * LOGICAL_BYTES_PER_REQUEST,
        logical_bytes_per_start: REQUESTS_PER_START * LOGICAL_BYTES_PER_REQUEST,
        peak_materialized_bytes: requests * RECORD_BYTES * 2 + starts * ACCUMULATOR_BYTES,
        incomplete_start_masks: incomplete,
        checksum_failures,
        negative_controls: starts,
        false_positives,
        false_negatives: checksum_failures + incomplete,
        ordered_request_sha256,
        routed_checksum_sha256,
    }
}

fn toy_cell(universe: u64, starts: u64) -> ToyCell {
    let routed = route(starts, universe);
    let mut reference = std::collections::BTreeSet::new();
    for start in 0..starts {
        for slot in 0..REQUESTS_PER_START as u8 {
            reference.insert(pair_id(start, slot, universe));
        }
    }
    let mut route_failures = 0u64;
    let mut false_positives = 0u64;
    for start in 0..starts {
        let expected = expected_accumulator(start, universe);
        let actual = routed.accumulators[start as usize];
        if actual.checksum != expected.checksum || actual.mask != 0xff || actual.count != 8 {
            route_failures += 1;
        }
        if actual.checksum.wrapping_add(1) == expected.checksum {
            false_positives += 1;
        }
    }
    ToyCell {
        pair_universe: universe,
        starts,
        requests: starts * REQUESTS_PER_START,
        radix_distinct: routed.distinct,
        reference_distinct: reference.len() as u64,
        route_failures,
        false_positives,
        false_negatives: route_failures,
    }
}

fn scalar_mul(
    point: &Projective,
    scalar: &BigUint,
    operations: &mut GroupOperations,
) -> Projective {
    operations.scalar_multiplications += 1;
    let mut result = Projective::IDENTITY;
    for bit in (0..scalar.bits()).rev() {
        result = result.double();
        operations.doublings += 1;
        if scalar.bit(bit) {
            result = result.add(point);
            operations.additions += 1;
        }
    }
    result
}

fn points_equal(left: &Projective, right: &Projective, operations: &mut GroupOperations) -> bool {
    operations.additions += 1;
    bool::from(left.add(&right.neg()).is_identity())
}

fn endpoint_scalar(index: u64, modulus: &BigUint) -> BigUint {
    let digest = sha256(format!("{CURVE_SLUG}/round40/endpoint/{index}").as_bytes());
    let mut scalar = BigUint::from_bytes_be(&digest) % modulus;
    if scalar.is_zero() {
        scalar = BigUint::one();
    }
    scalar
}

fn pair_prefix(left: u64) -> u64 {
    left * (2 * COLUMNS - left - 1) / 2
}

fn decode_pair(pair: u64) -> (u64, u64, bool, bool) {
    let ordinal = pair / 4;
    let signs = pair % 4;
    let mut low = 0u64;
    let mut high = COLUMNS - 1;
    while low + 1 < high {
        let middle = low + (high - low) / 2;
        if pair_prefix(middle) <= ordinal {
            low = middle;
        } else {
            high = middle;
        }
    }
    let left = low;
    let right = left + 1 + ordinal - pair_prefix(left);
    (left, right, signs & 1 == 0, signs & 2 == 0)
}

fn native_controls() -> NativeControls {
    let curve = CurveParams::p256();
    let generator = Projective::from_affine(
        &U256::from_biguint(&curve.gx),
        &U256::from_biguint(&curve.gy),
    );
    let mut operations = GroupOperations::default();
    let mut decode_failures = 0u64;
    let mut group_replay_failures = 0u64;
    let mut false_positives = 0u64;
    let mut stream = Vec::new();
    for control in 0..NATIVE_PAIR_CONTROLS {
        let pair = pair_id(control, (control & 7) as u8, PAIR_UNIVERSE);
        let (left, right, left_positive, right_positive) = decode_pair(pair);
        if left >= right || right >= COLUMNS {
            decode_failures += 1;
            continue;
        }
        let left_scalar = endpoint_scalar(left, &curve.n);
        let right_scalar = endpoint_scalar(right, &curve.n);
        let left_point = scalar_mul(&generator, &left_scalar, &mut operations);
        let right_point = scalar_mul(&generator, &right_scalar, &mut operations);
        operations.additions += 1;
        let direct = if left_positive {
            if right_positive {
                left_point.add(&right_point)
            } else {
                left_point.add(&right_point.neg())
            }
        } else if right_positive {
            left_point.neg().add(&right_point)
        } else {
            left_point.neg().add(&right_point.neg())
        };
        let mut combined = if left_positive {
            left_scalar
        } else {
            (&curve.n - left_scalar) % &curve.n
        };
        combined = if right_positive {
            (combined + right_scalar) % &curve.n
        } else {
            (combined + &curve.n - right_scalar) % &curve.n
        };
        let aggregate = scalar_mul(&generator, &combined, &mut operations);
        if !points_equal(&direct, &aggregate, &mut operations) {
            group_replay_failures += 1;
        }
        let corrupted = aggregate.add(&generator);
        operations.additions += 1;
        if points_equal(&direct, &corrupted, &mut operations) {
            false_positives += 1;
        }
        stream.extend(pair.to_be_bytes());
        stream.extend(left.to_be_bytes());
        stream.extend(right.to_be_bytes());
        stream.push(u8::from(left_positive));
        stream.push(u8::from(right_positive));
    }
    NativeControls {
        requests_checked: NATIVE_PAIR_CONTROLS,
        decode_failures,
        group_replay_failures,
        false_positives,
        false_negatives: group_replay_failures + decode_failures,
        operations,
        control_sha256: hex::encode(sha256(&stream)),
    }
}

fn timed_baseline(additions: u64) -> (u64, u64) {
    let curve = CurveParams::p256();
    let generator = Projective::from_affine(
        &U256::from_biguint(&curve.gx),
        &U256::from_biguint(&curve.gy),
    );
    let mut accumulator = generator;
    let start = Instant::now();
    for _ in 0..additions {
        accumulator = black_box(accumulator.add(black_box(&generator)));
    }
    let elapsed = start.elapsed().as_nanos() as u64;
    let checksum = accumulator
        .to_affine()
        .map(|(x, _)| x.0[0])
        .unwrap_or_default();
    (elapsed, checksum)
}

fn timed_candidate(starts: u64) -> (u64, u64) {
    let start = Instant::now();
    let routed = black_box(route(starts, PAIR_UNIVERSE));
    let checksum = routed
        .accumulators
        .iter()
        .fold(routed.distinct, |value, row| value ^ row.checksum);
    (start.elapsed().as_nanos() as u64, black_box(checksum))
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_by(f64::total_cmp);
    values[values.len() / 2]
}

fn timing_evidence() -> TimingEvidence {
    let timing_starts = 1u64 << TIMING_DEPTH;
    let _ = timed_baseline(8_192);
    let _ = timed_candidate(1u64 << 10);
    let (aa_first, _) = timed_baseline(BASELINE_ADDITIONS);
    let (aa_second, _) = timed_baseline(BASELINE_ADDITIONS);
    let mut repetitions = Vec::new();
    for repetition in 0..TIMING_REPETITIONS {
        let (baseline_wall_ns, baseline_checksum, candidate_wall_ns, candidate_checksum, order) =
            if repetition & 1 == 0 {
                let (baseline_ns, baseline_sum) = timed_baseline(BASELINE_ADDITIONS);
                let (candidate_ns, candidate_sum) = timed_candidate(timing_starts);
                (
                    baseline_ns,
                    baseline_sum,
                    candidate_ns,
                    candidate_sum,
                    "baseline-candidate",
                )
            } else {
                let (candidate_ns, candidate_sum) = timed_candidate(timing_starts);
                let (baseline_ns, baseline_sum) = timed_baseline(BASELINE_ADDITIONS);
                (
                    baseline_ns,
                    baseline_sum,
                    candidate_ns,
                    candidate_sum,
                    "candidate-baseline",
                )
            };
        repetitions.push(TimingSample {
            repetition,
            order: order.into(),
            baseline_wall_ns,
            candidate_wall_ns,
            baseline_checksum,
            candidate_checksum,
        });
    }
    let baseline_ns_per_addition = median(
        repetitions
            .iter()
            .map(|row| row.baseline_wall_ns as f64 / BASELINE_ADDITIONS as f64)
            .collect(),
    );
    let candidate_ns_per_start = median(
        repetitions
            .iter()
            .map(|row| row.candidate_wall_ns as f64 / timing_starts as f64)
            .collect(),
    );
    let first_checksum = repetitions[0].candidate_checksum;
    let candidate_checksums_identical = repetitions
        .iter()
        .all(|row| row.candidate_checksum == first_checksum);
    let addition_equivalents = candidate_ns_per_start / baseline_ns_per_addition;
    TimingEvidence {
        timing_start_bits: TIMING_DEPTH,
        timing_starts,
        baseline_additions: BASELINE_ADDITIONS,
        warmup_completed: true,
        aa_first_wall_ns: aa_first,
        aa_second_wall_ns: aa_second,
        aa_ratio: aa_first.max(aa_second) as f64 / aa_first.min(aa_second) as f64,
        repetitions,
        baseline_median_ns_per_addition: baseline_ns_per_addition,
        candidate_median_ns_per_start: candidate_ns_per_start,
        candidate_median_addition_equivalents_per_start: addition_equivalents,
        candidate_budget_multiple: addition_equivalents / ALLOWED_ADDITIONS,
        candidate_projected_ratio_to_rho: pair_ratio(addition_equivalents),
        candidate_logical_throughput_bytes_per_second:
            (REQUESTS_PER_START * LOGICAL_BYTES_PER_REQUEST) as f64
                / (candidate_ns_per_start * 1e-9),
        candidate_checksums_identical,
        complete_ram_route_below_allowed_additions: addition_equivalents <= ALLOWED_ADDITIONS,
        scope: "Matched single-core RAM routing diagnostic at 2^16 starts. It includes SHA-256 request generation, allocation, three radix passes and reconstruction; it is not an external-memory or end-to-end selector measurement."
            .into(),
    }
}

fn pair_ratio(additions_per_start: f64) -> f64 {
    BASE_RATIO + LOCAL_RATIO * additions_per_start / (MEAN_CAPACITY + 1.0)
}

fn peak_bytes(starts: u64) -> u64 {
    starts * (2 * REQUESTS_PER_START * RECORD_BYTES + ACCUMULATOR_BYTES)
}

fn projection_row(label: &str, starts: u64) -> ProjectionRow {
    let expected = expected_distinct(starts);
    let per_start = expected / starts as f64;
    let ratio = pair_ratio(per_start);
    let traffic =
        u128::from(starts) * u128::from(REQUESTS_PER_START) * u128::from(LOGICAL_BYTES_PER_REQUEST);
    let headroom = ALLOWED_ADDITIONS - per_start;
    ProjectionRow {
        label: label.into(),
        starts,
        starts_log2: (starts as f64).log2(),
        requests: starts * REQUESTS_PER_START,
        expected_distinct_pairs: expected,
        expected_distinct_per_start: per_start,
        pair_construction_ratio_to_rho: ratio,
        pair_construction_reaches_rho: ratio <= 1.0,
        peak_materialized_bytes: peak_bytes(starts),
        peak_materialized_log2_bytes: (peak_bytes(starts) as f64).log2(),
        storage_below_2_50: peak_bytes(starts) < (1u64 << 50),
        logical_routing_bytes: traffic.to_string(),
        logical_routing_log2_bytes: (traffic as f64).log2(),
        logical_bytes_per_start: REQUESTS_PER_START * LOGICAL_BYTES_PER_REQUEST,
        routing_addition_headroom_per_start: headroom,
        required_routing_bytes_per_addition_time: (headroom > 0.0)
            .then_some((REQUESTS_PER_START * LOGICAL_BYTES_PER_REQUEST) as f64 / headroom),
        applicable_tier_measured: false,
    }
}

fn parity_threshold() -> u64 {
    let mut low = 1u64;
    let mut high = 1u64 << 50;
    while low + 1 < high {
        let middle = low + (high - low) / 2;
        if expected_distinct(middle) / middle as f64 <= ALLOWED_ADDITIONS {
            high = middle;
        } else {
            low = middle;
        }
    }
    high
}

fn result(cli: &Cli) -> Result<ResultFile, String> {
    if size_of::<Request>() != RECORD_BYTES as usize
        || size_of::<Accumulator>() != ACCUMULATOR_BYTES as usize
    {
        return Err("registered record layout changed".into());
    }
    let (round19, value19) = load_dependency(
        &cli.round19,
        19,
        ROUND19_SHA256,
        "p256.s17_multilevel_selector/v1",
    )?;
    let (round36, value36) = load_dependency(
        &cli.round36,
        36,
        ROUND36_SHA256,
        "p256-mont-affine-table-v1",
    )?;
    let (round38, value38) =
        load_dependency(&cli.round38, 38, ROUND38_SHA256, "p256-codebook-support-v1")?;
    if json_str(&value19, &["curve"])? != CURVE_SLUG
        || json_str(&value19, &["factor_base", "fb_id"])? != COMPARISON_FB
        || json_u64(
            &value19,
            &["residual_degree_evidence", "structured_residual_maximum"],
        )? != 4
        || json_str(&value36, &["curve"])? != CURVE_SLUG
        || json_u64(&value36, &["frozen_model", "pair_entries"])? != PAIR_UNIVERSE
        || json_str(&value38, &["curve"])? != CURVE_SLUG
        || json_u64(&value38, &["frozen_model", "columns"])? != COLUMNS
        || (json_f64(
            &value36,
            &[
                "frozen_model",
                "round33_ratio_to_rho_with_target_correction",
            ],
        )? - BASE_RATIO)
            .abs()
            > 1e-15
        || (json_f64(
            &value36,
            &["frozen_model", "allowed_extra_additions_per_segment"],
        )? - ALLOWED_ADDITIONS)
            .abs()
            > 1e-15
    {
        return Err("frozen dependency facts changed".into());
    }

    let toy_controls = [17, 31, 63]
        .into_iter()
        .flat_map(|universe| [4, 8, 16].map(move |starts| toy_cell(universe, starts)))
        .collect::<Vec<_>>();
    let depth_ladder = REGISTERED_DEPTHS
        .into_iter()
        .map(checked_cell)
        .collect::<Vec<_>>();
    let native = native_controls();
    let timing = timing_evidence();
    let threshold = parity_threshold();
    let maximum_storage_starts =
        ((1u64 << 50) - 1) / (2 * REQUESTS_PER_START * RECORD_BYTES + ACCUMULATOR_BYTES);
    let mut projection = vec![
        projection_row("pair-only parity threshold", threshold),
        projection_row(
            "largest batch below 2^50 materialized bytes",
            maximum_storage_starts,
        ),
    ];
    for bits in [34u32, 36, 38, 40] {
        projection.push(projection_row(
            &format!("power-of-two-2^{bits}"),
            1u64 << bits,
        ));
    }
    projection.sort_by_key(|row| row.starts);
    let selected = projection
        .iter()
        .find(|row| row.label == "largest batch below 2^50 materialized bytes")
        .ok_or("selected storage row missing")?;
    let order = BigUint::parse_bytes(GROUP_ORDER.as_bytes(), 10).ok_or("invalid group order")?;
    let total_starts_log2 = 0.5 * (order.bits() as f64 - 1.0)
        + 0.5
        + (BASE_RATIO * 1.3 / 1.035_498_560_753_378_6).log2();
    let full_projection = FullProjection {
        total_independent_starts_log2: total_starts_log2,
        batches_at_selected_storage_row_log2: total_starts_log2 - selected.starts_log2,
        complete_logical_routing_bytes_log2: total_starts_log2
            + ((REQUESTS_PER_START * LOGICAL_BYTES_PER_REQUEST) as f64).log2(),
        useful_relation_rows_requested: 138_031,
        dickson_relation_collection_operations: None,
        dickson_cost_per_usable_relation: None,
        sparse_linear_algebra_operations: None,
        classification: "The total-start count is the imported known-log one-target model. It is not a projection for collecting 138,031 FB1h2f8621cda105 rows; those fields remain unset."
            .into(),
    };
    let degree_boundary = DegreeBoundary {
        comparison_factor_base: COMPARISON_FB.into(),
        comparison_structured_residual_maximum: 4,
        comparison_degree_gate_passed: true,
        measured_selector_family: "round-33 known-log low-delta/scalar selector".into(),
        measured_family_structured_residual_degree: None,
        measured_family_degree_gate_passed: false,
        unsplit_s17_degree_of_regularity: None,
    };
    let toy_failures = toy_controls
        .iter()
        .map(|row| {
            row.route_failures
                + row.false_positives
                + row.false_negatives
                + u64::from(row.radix_distinct != row.reference_distinct)
        })
        .sum::<u64>();
    let ladder_failures = depth_ladder
        .iter()
        .map(|row| {
            row.incomplete_start_masks
                + row.checksum_failures
                + row.false_positives
                + row.false_negatives
        })
        .sum::<u64>();
    let native_failures = native.decode_failures
        + native.group_replay_failures
        + native.false_positives
        + native.false_negatives;
    let correctness = toy_failures == 0 && ladder_failures == 0 && native_failures == 0;
    let storage_gate = projection
        .iter()
        .any(|row| row.pair_construction_reaches_rho && row.storage_below_2_50);
    let semantic = serde_json::json!({
        "curve": CURVE_SLUG,
        "dependencies": [ROUND19_SHA256, ROUND36_SHA256, ROUND38_SHA256],
        "toy_controls": &toy_controls,
        "depth_ladder": &depth_ladder,
        "native_controls": &native,
        "timing": &timing,
        "projection": &projection,
        "full_projection": &full_projection,
        "degree_boundary": &degree_boundary,
    });
    let semantic_evidence_sha256 = hex::encode(sha256(
        &serde_json::to_vec(&semantic).map_err(|error| error.to_string())?,
    ));
    let gates = Gates {
        zero_false_positives_and_false_negatives: correctness,
        exact_group_replay: native.group_replay_failures == 0,
        structured_residual_degree_at_most_5_for_measured_family: false,
        relation_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        materialized_storage_below_2_50: storage_gate,
        measured_complete_time_below_rho_in_projected_tier: false,
        no_discarded_branch_counted_as_exhaustive: true,
        promoted: false,
    };
    Ok(ResultFile {
        schema: "p256-independent-batch-routing-v1".into(),
        curve: CURVE_SLUG.into(),
        measured_family: "round-33 known-log low-delta/scalar selector".into(),
        comparison_factor_base: COMPARISON_FB.into(),
        dependencies: vec![round19, round36, round38],
        frozen_model: serde_json::json!({
            "columns": COLUMNS,
            "pair_universe": PAIR_UNIVERSE,
            "requests_per_start": REQUESTS_PER_START,
            "record_bytes": RECORD_BYTES,
            "accumulator_bytes": ACCUMULATOR_BYTES,
            "radix_bits": RADIX_BITS,
            "radix_passes": RADIX_PASSES,
            "logical_bytes_per_request": LOGICAL_BYTES_PER_REQUEST,
            "base_ratio_to_rho": BASE_RATIO,
            "local_ratio_to_rho": LOCAL_RATIO,
            "mean_capacity": MEAN_CAPACITY,
            "allowed_extra_additions_per_start": ALLOWED_ADDITIONS,
        }),
        toy_controls,
        depth_ladder,
        native_controls: native,
        timing,
        projection,
        full_projection,
        degree_boundary,
        semantic_evidence_sha256,
        gates,
        full_depth_unplanted_attempted: false,
        interpretation: "Independent batching preserves start entropy and can amortize pair construction below the frozen rho arithmetic boundary in expectation. The crossing requires batches far beyond the measured RAM tier; exact routing contributes 1,280 logical bytes per start and no applicable external-memory timing exists."
            .into(),
        decision: "Do not promote. Preserve independent batching as a surviving pair-construction idea, but parity is not established: the measured family has unknown structured degree, the Dickson relation projection is unset, and multi-terabyte routing is unmeasured."
            .into(),
    })
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match result(&cli) {
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
                if let Err(error) = fs::write(&path, bytes) {
                    eprintln!("{}: {error}", path.display());
                    return ExitCode::FAILURE;
                }
            } else {
                print!("{}", String::from_utf8_lossy(&bytes));
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("P-256 independent batch audit failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn radix_router_matches_complete_reference() {
        for universe in [17, 31, 63] {
            let row = toy_cell(universe, 16);
            assert_eq!(row.radix_distinct, row.reference_distinct);
            assert_eq!(row.route_failures, 0);
            assert_eq!(row.false_positives, 0);
            assert_eq!(row.false_negatives, 0);
        }
    }

    #[test]
    fn pair_decoder_covers_full_universe() {
        for pair in [0, 1, PAIR_UNIVERSE / 2, PAIR_UNIVERSE - 1] {
            let (left, right, _, _) = decode_pair(pair);
            assert!(left < right);
            assert!(right < COLUMNS);
        }
        assert_eq!(pair_prefix(COLUMNS - 1) * 4, PAIR_UNIVERSE);
    }

    #[test]
    fn registered_layout_and_traffic_are_fixed() {
        assert_eq!(size_of::<Request>(), 16);
        assert_eq!(size_of::<Accumulator>(), 16);
        assert_eq!(LOGICAL_BYTES_PER_REQUEST, 160);
        assert_eq!(REQUESTS_PER_START * LOGICAL_BYTES_PER_REQUEST, 1_280);
    }

    #[test]
    fn expected_occupancy_is_monotone_per_start() {
        let shallow = expected_distinct(1 << 20) / (1u64 << 20) as f64;
        let deep = expected_distinct(1 << 40) / (1u64 << 40) as f64;
        assert!(deep < shallow);
        assert!(shallow <= 8.0);
    }
}
