//! Exact P-256 radix-width portfolio for rounds 43 through 290.

use std::collections::BTreeSet;
use std::fs;
use std::hint::black_box;
use std::mem::{size_of, swap};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use blake3::Hasher;
use clap::Parser;
use crypto_lib::ct_bignum::U256;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use serde::Serialize;
use serde_json::Value;

const CURVE_SLUG: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const COMPARISON_FB: &str = "FB1h2f8621cda105";
const ROUND42_SHA256: &str = "82adf1c47c60a9004a83ac061d9a75ce908f3bb9c28f537260b43b7d5016aa2a";
const ROUND42_ISOLATION_SHA256: &str =
    "7a6d11512161cee3020f599e8e705f6ea8fa68ba4e14aada73e0716230d4dad8";
const REQUEST_DIGEST: &str = "b4178238067a37febeeb46116dffc633311549e6325461670bd12db7383e0dd6";
const CONTROL_DIGEST: &str = "98c8f56a06257256bb88f6c6d7cffce79391c461170b91287f8fd193fe5711e6";
const PAIR_UNIVERSE: u64 = 34_562_148_612;
const KEY_BITS: u32 = 36;
const REQUESTS_PER_START: u64 = 8;
const RECORD_BYTES: u64 = 16;
const ACCUMULATOR_BYTES: u64 = 16;
const MIN_WIDTH: u32 = 6;
const MAX_WIDTH: u32 = 36;
const START_BITS: [u32; 8] = [12, 14, 16, 18, 20, 24, 32, 40];
const MAX_EXECUTED_WIDTH: u32 = 24;
const MAX_EXECUTED_START_BITS: u32 = 18;
const FIRST_ROUND: u32 = 43;
const LAST_ROUND: u32 = 290;
const PORTFOLIO_CELLS: usize = 248;
const BASELINE_ADDITIONS: u64 = 262_144;
const TIMING_REPETITIONS: u64 = 7;
const ALLOWED_ADDITIONS: f64 = 0.334_410_508_422_987_64;
const BASE_RATIO: f64 = 0.998_569_034_150_286;
const LOCAL_RATIO: f64 = 0.964_336_477_130_181;
const MEAN_CAPACITY: f64 = 224.361_249_307_510_55;

#[derive(Parser)]
#[command(about = "Execute P-256 radix portfolio rounds 43 through 290")]
struct Cli {
    #[arg(long)]
    round42: PathBuf,
    #[arg(long)]
    isolation: PathBuf,
    #[arg(long)]
    round_dir: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
struct Request {
    pair: u64,
    owner_slot: u64,
}

#[repr(C)]
#[derive(Clone, Copy, Default, PartialEq, Eq)]
struct Accumulator {
    checksum: u64,
    mask: u8,
    count: u8,
    padding: [u8; 6],
}

#[derive(Clone, Serialize)]
struct Dependency {
    path: String,
    sha256: String,
    schema: String,
}

#[derive(Clone, Serialize)]
struct ImportedEvidence {
    round42_request_digest: String,
    round42_native_control_digest: String,
    round42_zero_false_positives_and_false_negatives: bool,
    round42_exact_reference_equality: bool,
    round42_full_universe_algorithmic_bytes: u64,
    round42_materialized_complete_addition_equivalents_per_start: f64,
    round42_direct_complete_addition_equivalents_per_start: f64,
    round42_promoted: bool,
    round42_isolated_runs: u64,
    round42_last_run_uncontended: bool,
}

#[derive(Clone, Serialize)]
struct OperationCounts {
    request_generations_algorithmic: u64,
    request_generations_timed: u64,
    histogram_increments: u64,
    scatter_writes: u64,
    histogram_bucket_scans: u64,
    sorted_record_scans: u64,
    accumulator_updates: u64,
    native_field_operations: u64,
    native_group_operations: u64,
}

#[derive(Clone, Serialize)]
struct Traffic {
    record_bytes: String,
    histogram_bytes: String,
    total_logical_bytes: String,
    disk_bytes: u64,
}

#[derive(Clone, Serialize)]
struct CellGates {
    zero_false_positives_and_false_negatives_for_this_cell: bool,
    structured_residual_degree_at_most_5_for_measured_family: bool,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    measured_complete_selector_below_rho_in_applicable_tier: bool,
    no_discarded_branch_counted_as_exhaustive: bool,
    promoted: bool,
}

#[derive(Clone, Serialize)]
struct ShortlistSummary {
    repetitions: u64,
    median_wall_ns: f64,
    minimum_wall_ns: u64,
    median_ns_per_start: f64,
    routing_addition_equivalents_per_start: f64,
    routing_budget_multiple: f64,
    routing_only_ratio_to_rho: f64,
    checksum_identical: bool,
    excludes_request_generation: bool,
}

#[derive(Clone, Serialize)]
struct RoundRecord {
    schema: String,
    curve: String,
    round: u32,
    family: String,
    comparison_factor_base: String,
    radix_width: u32,
    start_bits: u32,
    starts: u64,
    requests: u64,
    pass_count: u32,
    histogram_bits: Vec<u32>,
    histogram_buckets: Vec<u64>,
    largest_histogram_counters: u64,
    complete_checked: bool,
    evidence_class: String,
    applicable_memory_tier_measured: bool,
    distinct_pairs: Option<u64>,
    expected_distinct_pairs: f64,
    sorted_mismatches: Option<u64>,
    accumulator_mismatches: Option<u64>,
    regenerated_failures: Option<u64>,
    negative_controls: Option<u64>,
    false_positives: Option<u64>,
    false_negatives: Option<u64>,
    sorted_blake3: Option<String>,
    accumulator_blake3: Option<String>,
    reference_sorted_blake3: Option<String>,
    reference_accumulator_blake3: Option<String>,
    one_shot_wall_ns: Option<u64>,
    one_shot_routing_addition_equivalents_per_start: Option<f64>,
    operation_counts: OperationCounts,
    traffic: Traffic,
    peak_algorithmic_bytes: u64,
    peak_log2_bytes: f64,
    time_exponent_across_complete_depths: Option<f64>,
    memory_exponent_across_registered_depths: f64,
    shortlist_timing: Option<ShortlistSummary>,
    structured_residual_degree_for_measured_family: Option<u32>,
    gates: CellGates,
    full_depth_unplanted_attempted: bool,
    interpretation: String,
}

#[derive(Clone, Serialize)]
struct ManifestEntry {
    round: u32,
    radix_width: u32,
    start_bits: u32,
    complete_checked: bool,
    path: String,
    bytes: usize,
    sha256: String,
}

#[derive(Clone, Serialize)]
struct RawShortlistMeasurement {
    width: u32,
    wall_ns: u64,
    checksum: u64,
}

#[derive(Clone, Serialize)]
struct RawShortlistRepetition {
    repetition: u64,
    order: Vec<u32>,
    measurements: Vec<RawShortlistMeasurement>,
}

#[derive(Clone, Serialize)]
struct TimingEvidence {
    baseline_additions: u64,
    aa_first_ns: u64,
    aa_second_ns: u64,
    aa_ratio: f64,
    baseline_repetitions: u64,
    baseline_median_ns_per_addition: f64,
    shortlist_depth_bits: u32,
    shortlist_widths: Vec<u32>,
    shortlist_repetitions: u64,
    raw_shortlist: Vec<RawShortlistRepetition>,
    summaries: Vec<(u32, ShortlistSummary)>,
    time_exponents: Vec<(u32, f64)>,
    memory_exponents: Vec<(u32, f64)>,
    scope: String,
}

#[derive(Clone, Serialize)]
struct PortfolioGates {
    all_248_rounds_present_once: bool,
    all_76_registered_cells_complete: bool,
    zero_false_positives_and_false_negatives_on_complete_cells: bool,
    structured_residual_degree_at_most_5_for_measured_family: bool,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    at_least_one_storage_projection_below_2_50: bool,
    measured_complete_selector_below_rho_in_applicable_tier: bool,
    no_projection_counted_as_exhaustive: bool,
    promoted_rounds: Vec<u32>,
}

#[derive(Serialize)]
struct Manifest {
    schema: String,
    curve: String,
    round_range: [u32; 2],
    family: String,
    dependencies: Vec<Dependency>,
    imported: ImportedEvidence,
    frozen_matrix: Value,
    timing: TimingEvidence,
    entries: Vec<ManifestEntry>,
    complete_checked_rounds: Vec<u32>,
    projection_only_rounds: Vec<u32>,
    best_measured_kernel_round: u32,
    best_measured_kernel_width: u32,
    best_formula_traffic_round: u32,
    best_formula_traffic_width: u32,
    gates: PortfolioGates,
    semantic_evidence_sha256: String,
    full_depth_unplanted_attempted: bool,
    interpretation: String,
    decision: String,
}

struct Route {
    sorted: Vec<Request>,
    accumulators: Vec<Accumulator>,
    distinct: u64,
    sorted_digest: String,
    accumulator_digest: String,
}

struct ExactContext {
    frozen: Vec<Request>,
    canonical_sorted: Vec<Request>,
    expected_accumulators: Vec<Accumulator>,
    reference_regenerated_failures: u64,
    reference_distinct: u64,
    reference_sorted_digest: String,
    reference_accumulator_digest: String,
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

fn json_bool(value: &Value, path: &[&str]) -> Result<bool, String> {
    json_at(value, path)?
        .as_bool()
        .ok_or_else(|| format!("JSON path {} is not a bool", path.join(".")))
}

fn json_u64(value: &Value, path: &[&str]) -> Result<u64, String> {
    json_at(value, path)?
        .as_u64()
        .ok_or_else(|| format!("JSON path {} is not a u64", path.join(".")))
}

fn load_round42(path: &Path) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != ROUND42_SHA256 {
        return Err(format!(
            "round-42 hash mismatch: expected {ROUND42_SHA256}, got {digest}"
        ));
    }
    let value: Value = serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if json_str(&value, &["schema"])? != "p256-direct-address-replay-v1"
        || json_str(&value, &["curve"])? != CURVE_SLUG
    {
        return Err("round-42 identity mismatch".into());
    }
    Ok((
        Dependency {
            path: path.display().to_string(),
            sha256: digest,
            schema: "p256-direct-address-replay-v1".into(),
        },
        value,
    ))
}

fn load_isolation(path: &Path) -> Result<(Dependency, Vec<Value>), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != ROUND42_ISOLATION_SHA256 {
        return Err(format!(
            "round-42 isolation hash mismatch: expected {ROUND42_ISOLATION_SHA256}, got {digest}"
        ));
    }
    let text = String::from_utf8(bytes).map_err(|error| error.to_string())?;
    let records = text
        .lines()
        .map(|line| serde_json::from_str(line).map_err(|error| error.to_string()))
        .collect::<Result<Vec<Value>, String>>()?;
    if records.is_empty() {
        return Err("round-42 isolation has no records".into());
    }
    Ok((
        Dependency {
            path: path.display().to_string(),
            sha256: digest,
            schema: "isolated-bench/1-jsonl".into(),
        },
        records,
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

fn generate_requests(starts: u64) -> Vec<Request> {
    let mut records = Vec::with_capacity((starts * REQUESTS_PER_START) as usize);
    for start in 0..starts {
        for slot in 0..REQUESTS_PER_START as u8 {
            records.push(Request {
                pair: pair_id(start, slot, PAIR_UNIVERSE),
                owner_slot: (start << 3) | u64::from(slot),
            });
        }
    }
    records
}

fn histogram_bits(width: u32) -> Vec<u32> {
    let mut remaining = KEY_BITS;
    let mut bits = Vec::new();
    while remaining > 0 {
        let digit = width.min(remaining);
        bits.push(digit);
        remaining -= digit;
    }
    bits
}

fn histogram_buckets(width: u32) -> Vec<u64> {
    histogram_bits(width)
        .into_iter()
        .map(|bits| 1u64 << bits)
        .collect()
}

fn radix_sort(mut source: Vec<Request>, width: u32) -> Vec<Request> {
    let mut target = vec![Request::default(); source.len()];
    let mut shift = 0u32;
    for bits in histogram_bits(width) {
        let buckets = 1usize << bits;
        let mask = buckets - 1;
        let mut counts = vec![0usize; buckets];
        for record in &source {
            counts[((record.pair >> shift) as usize) & mask] += 1;
        }
        let mut offset = 0usize;
        for count in &mut counts {
            let current = *count;
            *count = offset;
            offset += current;
        }
        for record in &source {
            let digit = ((record.pair >> shift) as usize) & mask;
            target[counts[digit]] = *record;
            counts[digit] += 1;
        }
        swap(&mut source, &mut target);
        shift += bits;
    }
    source
}

fn digest_requests(records: &[Request]) -> String {
    let mut hasher = Hasher::new();
    for record in records {
        hasher.update(&record.pair.to_be_bytes());
        hasher.update(&record.owner_slot.to_be_bytes());
    }
    hasher.finalize().to_hex().to_string()
}

fn digest_accumulators(accumulators: &[Accumulator]) -> String {
    let mut hasher = Hasher::new();
    for accumulator in accumulators {
        hasher.update(&accumulator.checksum.to_be_bytes());
        hasher.update(&[accumulator.mask, accumulator.count]);
    }
    hasher.finalize().to_hex().to_string()
}

fn reconstruct(sorted: Vec<Request>, starts: u64) -> Route {
    let mut accumulators = vec![Accumulator::default(); starts as usize];
    let mut distinct = 0u64;
    let mut previous = None;
    for record in &sorted {
        if previous != Some(record.pair) {
            distinct += 1;
            previous = Some(record.pair);
        }
        let owner = (record.owner_slot >> 3) as usize;
        let slot = (record.owner_slot & 7) as u8;
        let accumulator = &mut accumulators[owner];
        accumulator.checksum = accumulator
            .checksum
            .wrapping_add(contribution(record.pair, slot));
        accumulator.mask |= 1 << slot;
        accumulator.count += 1;
    }
    let sorted_digest = digest_requests(&sorted);
    let accumulator_digest = digest_accumulators(&accumulators);
    Route {
        sorted,
        accumulators,
        distinct,
        sorted_digest,
        accumulator_digest,
    }
}

fn expected_accumulator(start: u64) -> Accumulator {
    let mut accumulator = Accumulator::default();
    for slot in 0..REQUESTS_PER_START as u8 {
        let pair = pair_id(start, slot, PAIR_UNIVERSE);
        accumulator.checksum = accumulator.checksum.wrapping_add(contribution(pair, slot));
        accumulator.mask |= 1 << slot;
        accumulator.count += 1;
    }
    accumulator
}

fn exact_context(start_bits: u32) -> ExactContext {
    let starts = 1u64 << start_bits;
    let frozen = generate_requests(starts);
    let mut canonical_sorted = frozen.clone();
    canonical_sorted.sort_unstable();
    let reference = reconstruct(canonical_sorted.clone(), starts);
    let expected_accumulators = (0..starts).map(expected_accumulator).collect::<Vec<_>>();
    let reference_regenerated_failures = reference
        .accumulators
        .iter()
        .zip(&expected_accumulators)
        .filter(|(left, right)| left != right)
        .count() as u64;
    ExactContext {
        frozen,
        canonical_sorted,
        expected_accumulators,
        reference_regenerated_failures,
        reference_distinct: reference.distinct,
        reference_sorted_digest: reference.sorted_digest,
        reference_accumulator_digest: reference.accumulator_digest,
    }
}

fn expected_distinct(starts: u64) -> f64 {
    let requests = (starts * REQUESTS_PER_START) as f64;
    let universe = PAIR_UNIVERSE as f64;
    -universe * (requests * (-1.0 / universe).ln_1p()).exp_m1()
}

fn round_number(depth_index: usize, width: u32) -> u32 {
    FIRST_ROUND + 31 * depth_index as u32 + (width - MIN_WIDTH)
}

fn peak_bytes(width: u32, starts: u64) -> u64 {
    let requests = starts * REQUESTS_PER_START;
    let max_buckets = histogram_buckets(width).into_iter().max().unwrap();
    2 * RECORD_BYTES * requests + ACCUMULATOR_BYTES * starts + 8 * max_buckets
}

fn traffic(width: u32, starts: u64) -> Traffic {
    let requests = u128::from(starts) * u128::from(REQUESTS_PER_START);
    let passes = u128::from(histogram_bits(width).len() as u64);
    let record = requests * (64 + 32 * passes);
    let histogram = histogram_buckets(width)
        .into_iter()
        .map(|buckets| u128::from(buckets) * 24)
        .sum::<u128>();
    Traffic {
        record_bytes: record.to_string(),
        histogram_bytes: histogram.to_string(),
        total_logical_bytes: (record + histogram).to_string(),
        disk_bytes: 0,
    }
}

fn operations(width: u32, starts: u64) -> OperationCounts {
    let requests = starts * REQUESTS_PER_START;
    let passes = histogram_bits(width).len() as u64;
    OperationCounts {
        request_generations_algorithmic: requests,
        request_generations_timed: 0,
        histogram_increments: requests * passes,
        scatter_writes: requests * passes,
        histogram_bucket_scans: histogram_buckets(width).into_iter().sum(),
        sorted_record_scans: requests,
        accumulator_updates: requests,
        native_field_operations: 0,
        native_group_operations: 0,
    }
}

fn pair_ratio(additions: f64) -> f64 {
    BASE_RATIO + LOCAL_RATIO * additions / (MEAN_CAPACITY + 1.0)
}

fn digest_checksum(route: &Route) -> u64 {
    let bytes = hex::decode(&route.accumulator_digest[..16]).expect("digest hex");
    let mut word = [0u8; 8];
    word.copy_from_slice(&bytes);
    u64::from_be_bytes(word) ^ route.distinct
}

fn timed_route(frozen: &[Request], starts: u64, width: u32) -> (u64, u64) {
    let begin = Instant::now();
    let route = black_box(reconstruct(radix_sort(frozen.to_vec(), width), starts));
    (
        begin.elapsed().as_nanos() as u64,
        black_box(digest_checksum(&route)),
    )
}

fn timed_additions(count: u64) -> (u64, u64) {
    let curve = CurveParams::p256();
    let generator = Projective::from_affine(
        &U256::from_biguint(&curve.gx),
        &U256::from_biguint(&curve.gy),
    );
    let mut accumulator = generator;
    let begin = Instant::now();
    for _ in 0..count {
        accumulator = black_box(accumulator.add(black_box(&generator)));
    }
    let elapsed = begin.elapsed().as_nanos() as u64;
    let checksum = accumulator
        .to_affine()
        .map(|(x, _)| x.0[0])
        .unwrap_or_default();
    (elapsed, checksum)
}

fn median_u64(values: &[u64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    sorted[sorted.len() / 2] as f64
}

fn fit_exponent(points: &[(f64, f64)]) -> f64 {
    let mean_x = points.iter().map(|point| point.0).sum::<f64>() / points.len() as f64;
    let mean_y = points.iter().map(|point| point.1).sum::<f64>() / points.len() as f64;
    points
        .iter()
        .map(|point| (point.0 - mean_x) * (point.1 - mean_y))
        .sum::<f64>()
        / points
            .iter()
            .map(|point| (point.0 - mean_x).powi(2))
            .sum::<f64>()
}

fn projection_record(depth_index: usize, width: u32) -> RoundRecord {
    let start_bits = START_BITS[depth_index];
    let starts = 1u64 << start_bits;
    let requests = starts * REQUESTS_PER_START;
    let bits = histogram_bits(width);
    let buckets = histogram_buckets(width);
    let peak = peak_bytes(width, starts);
    let storage = peak < (1u64 << 50);
    RoundRecord {
        schema: "p256-radix-portfolio-round-v1".into(),
        curve: CURVE_SLUG.into(),
        round: round_number(depth_index, width),
        family: "round-33 known-log selector radix routing".into(),
        comparison_factor_base: COMPARISON_FB.into(),
        radix_width: width,
        start_bits,
        starts,
        requests,
        pass_count: bits.len() as u32,
        histogram_bits: bits,
        largest_histogram_counters: buckets.iter().copied().max().unwrap(),
        histogram_buckets: buckets,
        complete_checked: false,
        evidence_class: "projection_only".into(),
        applicable_memory_tier_measured: false,
        distinct_pairs: None,
        expected_distinct_pairs: expected_distinct(starts),
        sorted_mismatches: None,
        accumulator_mismatches: None,
        regenerated_failures: None,
        negative_controls: None,
        false_positives: None,
        false_negatives: None,
        sorted_blake3: None,
        accumulator_blake3: None,
        reference_sorted_blake3: None,
        reference_accumulator_blake3: None,
        one_shot_wall_ns: None,
        one_shot_routing_addition_equivalents_per_start: None,
        operation_counts: operations(width, starts),
        traffic: traffic(width, starts),
        peak_algorithmic_bytes: peak,
        peak_log2_bytes: (peak as f64).log2(),
        time_exponent_across_complete_depths: None,
        memory_exponent_across_registered_depths: 0.0,
        shortlist_timing: None,
        structured_residual_degree_for_measured_family: None,
        gates: CellGates {
            zero_false_positives_and_false_negatives_for_this_cell: false,
            structured_residual_degree_at_most_5_for_measured_family: false,
            relation_collection_below_2_120: false,
            per_usable_relation_below_2_103: false,
            materialized_storage_below_2_50: storage,
            measured_complete_selector_below_rho_in_applicable_tier: false,
            no_discarded_branch_counted_as_exhaustive: true,
            promoted: false,
        },
        full_depth_unplanted_attempted: false,
        interpretation: "Formula-only cell. Exact pass, traffic, and memory counts are registered, but no allocation, routing, exactness result, or applicable timing was performed."
            .into(),
    }
}

fn exact_record(
    depth_index: usize,
    width: u32,
    context: &ExactContext,
    baseline_ns_per_addition: f64,
) -> RoundRecord {
    let start_bits = START_BITS[depth_index];
    let starts = 1u64 << start_bits;
    let (wall_ns, _) = timed_route(&context.frozen, starts, width);
    let route = reconstruct(radix_sort(context.frozen.clone(), width), starts);
    let sorted_mismatches = route
        .sorted
        .iter()
        .zip(&context.canonical_sorted)
        .filter(|(left, right)| left != right)
        .count() as u64;
    let accumulator_mismatches = route
        .accumulators
        .iter()
        .zip(&context.expected_accumulators)
        .filter(|(left, right)| left != right)
        .count() as u64;
    let regenerated_failures = context.reference_regenerated_failures;
    let false_positives = route
        .accumulators
        .iter()
        .zip(&context.expected_accumulators)
        .filter(|(actual, expected)| {
            let mut negative = **actual;
            negative.checksum = negative.checksum.wrapping_add(1);
            negative == **expected
        })
        .count() as u64;
    let false_negatives = sorted_mismatches
        + accumulator_mismatches
        + regenerated_failures
        + u64::from(route.distinct != context.reference_distinct)
        + u64::from(route.sorted_digest != context.reference_sorted_digest)
        + u64::from(route.accumulator_digest != context.reference_accumulator_digest);
    let equivalents = wall_ns as f64 / starts as f64 / baseline_ns_per_addition;
    let bits = histogram_bits(width);
    let buckets = histogram_buckets(width);
    let peak = peak_bytes(width, starts);
    let exact = false_positives == 0 && false_negatives == 0;
    RoundRecord {
        schema: "p256-radix-portfolio-round-v1".into(),
        curve: CURVE_SLUG.into(),
        round: round_number(depth_index, width),
        family: "round-33 known-log selector radix routing".into(),
        comparison_factor_base: COMPARISON_FB.into(),
        radix_width: width,
        start_bits,
        starts,
        requests: starts * REQUESTS_PER_START,
        pass_count: bits.len() as u32,
        histogram_bits: bits,
        largest_histogram_counters: buckets.iter().copied().max().unwrap(),
        histogram_buckets: buckets,
        complete_checked: true,
        evidence_class: "complete_checked".into(),
        applicable_memory_tier_measured: true,
        distinct_pairs: Some(route.distinct),
        expected_distinct_pairs: expected_distinct(starts),
        sorted_mismatches: Some(sorted_mismatches),
        accumulator_mismatches: Some(accumulator_mismatches),
        regenerated_failures: Some(regenerated_failures),
        negative_controls: Some(starts),
        false_positives: Some(false_positives),
        false_negatives: Some(false_negatives),
        sorted_blake3: Some(route.sorted_digest),
        accumulator_blake3: Some(route.accumulator_digest),
        reference_sorted_blake3: Some(context.reference_sorted_digest.clone()),
        reference_accumulator_blake3: Some(context.reference_accumulator_digest.clone()),
        one_shot_wall_ns: Some(wall_ns),
        one_shot_routing_addition_equivalents_per_start: Some(equivalents),
        operation_counts: operations(width, starts),
        traffic: traffic(width, starts),
        peak_algorithmic_bytes: peak,
        peak_log2_bytes: (peak as f64).log2(),
        time_exponent_across_complete_depths: None,
        memory_exponent_across_registered_depths: 0.0,
        shortlist_timing: None,
        structured_residual_degree_for_measured_family: None,
        gates: CellGates {
            zero_false_positives_and_false_negatives_for_this_cell: exact,
            structured_residual_degree_at_most_5_for_measured_family: false,
            relation_collection_below_2_120: false,
            per_usable_relation_below_2_103: false,
            materialized_storage_below_2_50: peak < (1u64 << 50),
            measured_complete_selector_below_rho_in_applicable_tier: false,
            no_discarded_branch_counted_as_exhaustive: true,
            promoted: false,
        },
        full_depth_unplanted_attempted: false,
        interpretation: "Complete checked routing-kernel cell. Timing includes clone, allocation, radix passes, scan, and reconstruction but excludes the shared deterministic request generation."
            .into(),
    }
}

fn shortlist_timing(
    records: &[RoundRecord],
    baseline_ns_per_addition: f64,
) -> (
    Vec<u32>,
    Vec<RawShortlistRepetition>,
    Vec<(u32, ShortlistSummary)>,
) {
    let mut depth16 = records
        .iter()
        .filter(|record| record.start_bits == 16 && record.complete_checked)
        .collect::<Vec<_>>();
    depth16.sort_by(|left, right| {
        left.one_shot_wall_ns
            .cmp(&right.one_shot_wall_ns)
            .then(
                left.peak_algorithmic_bytes
                    .cmp(&right.peak_algorithmic_bytes),
            )
            .then(left.radix_width.cmp(&right.radix_width))
    });
    let widths = depth16
        .into_iter()
        .take(5)
        .map(|record| record.radix_width)
        .collect::<Vec<_>>();
    let starts = 1u64 << 16;
    let frozen = generate_requests(starts);
    let mut raw = Vec::new();
    for repetition in 0..TIMING_REPETITIONS {
        let rotate = repetition as usize % widths.len();
        let mut order = widths.clone();
        order.rotate_left(rotate);
        let measurements = order
            .iter()
            .map(|width| {
                let (wall_ns, checksum) = timed_route(&frozen, starts, *width);
                RawShortlistMeasurement {
                    width: *width,
                    wall_ns,
                    checksum,
                }
            })
            .collect();
        raw.push(RawShortlistRepetition {
            repetition,
            order,
            measurements,
        });
    }
    let summaries = widths
        .iter()
        .map(|width| {
            let measurements = raw
                .iter()
                .flat_map(|row| &row.measurements)
                .filter(|measurement| measurement.width == *width)
                .collect::<Vec<_>>();
            let values = measurements
                .iter()
                .map(|measurement| measurement.wall_ns)
                .collect::<Vec<_>>();
            let median = median_u64(&values);
            let ns_per_start = median / starts as f64;
            let equivalents = ns_per_start / baseline_ns_per_addition;
            (
                *width,
                ShortlistSummary {
                    repetitions: TIMING_REPETITIONS,
                    median_wall_ns: median,
                    minimum_wall_ns: *values.iter().min().expect("shortlist values"),
                    median_ns_per_start: ns_per_start,
                    routing_addition_equivalents_per_start: equivalents,
                    routing_budget_multiple: equivalents / ALLOWED_ADDITIONS,
                    routing_only_ratio_to_rho: pair_ratio(equivalents),
                    checksum_identical: measurements
                        .iter()
                        .all(|measurement| measurement.checksum == measurements[0].checksum),
                    excludes_request_generation: true,
                },
            )
        })
        .collect();
    (widths, raw, summaries)
}

fn apply_exponents_and_shortlist(
    records: &mut [RoundRecord],
    summaries: &[(u32, ShortlistSummary)],
) -> (Vec<(u32, f64)>, Vec<(u32, f64)>) {
    let mut time_exponents = Vec::new();
    let mut memory_exponents = Vec::new();
    for width in MIN_WIDTH..=MAX_WIDTH {
        let time_points = records
            .iter()
            .filter(|record| record.radix_width == width && record.complete_checked)
            .map(|record| {
                (
                    record.start_bits as f64,
                    (record.one_shot_wall_ns.expect("checked wall") as f64).log2(),
                )
            })
            .collect::<Vec<_>>();
        let time_exponent = if time_points.len() >= 2 {
            Some(fit_exponent(&time_points))
        } else {
            None
        };
        if let Some(exponent) = time_exponent {
            time_exponents.push((width, exponent));
        }
        let memory_points = records
            .iter()
            .filter(|record| record.radix_width == width)
            .map(|record| {
                (
                    record.start_bits as f64,
                    (record.peak_algorithmic_bytes as f64).log2(),
                )
            })
            .collect::<Vec<_>>();
        let memory_exponent = fit_exponent(&memory_points);
        memory_exponents.push((width, memory_exponent));
        for record in records
            .iter_mut()
            .filter(|record| record.radix_width == width)
        {
            record.time_exponent_across_complete_depths = time_exponent;
            record.memory_exponent_across_registered_depths = memory_exponent;
            if record.start_bits == 16 {
                record.shortlist_timing = summaries
                    .iter()
                    .find(|(summary_width, _)| *summary_width == width)
                    .map(|(_, summary)| summary.clone());
            }
        }
    }
    (time_exponents, memory_exponents)
}

fn serialize_pretty<T: Serialize>(value: &T) -> Result<Vec<u8>, String> {
    let mut bytes = serde_json::to_vec_pretty(value).map_err(|error| error.to_string())?;
    bytes.push(b'\n');
    Ok(bytes)
}

fn write_rounds(round_dir: &Path, records: &[RoundRecord]) -> Result<Vec<ManifestEntry>, String> {
    fs::create_dir_all(round_dir)
        .map_err(|error| format!("create {}: {error}", round_dir.display()))?;
    let mut entries = Vec::new();
    for record in records {
        let name = format!(
            "round{:03}-w{:02}-s{:02}.json",
            record.round, record.radix_width, record.start_bits
        );
        let path = round_dir.join(name);
        if path.exists() {
            return Err(format!("refusing to overwrite {}", path.display()));
        }
        let bytes = serialize_pretty(record)?;
        fs::write(&path, &bytes).map_err(|error| format!("{}: {error}", path.display()))?;
        entries.push(ManifestEntry {
            round: record.round,
            radix_width: record.radix_width,
            start_bits: record.start_bits,
            complete_checked: record.complete_checked,
            path: path.display().to_string(),
            bytes: bytes.len(),
            sha256: hex::encode(sha256(&bytes)),
        });
    }
    Ok(entries)
}

fn validate_portfolio(records: &[RoundRecord]) -> Result<(), String> {
    if records.len() != PORTFOLIO_CELLS {
        return Err(format!("expected {PORTFOLIO_CELLS} cells"));
    }
    let rounds = records
        .iter()
        .map(|record| record.round)
        .collect::<BTreeSet<_>>();
    let parameters = records
        .iter()
        .map(|record| (record.radix_width, record.start_bits))
        .collect::<BTreeSet<_>>();
    if rounds.len() != PORTFOLIO_CELLS
        || parameters.len() != PORTFOLIO_CELLS
        || rounds.first() != Some(&FIRST_ROUND)
        || rounds.last() != Some(&LAST_ROUND)
    {
        return Err("round map is incomplete or duplicated".into());
    }
    let complete = records
        .iter()
        .filter(|record| record.complete_checked)
        .count();
    if complete != 76 {
        return Err(format!("expected 76 complete cells, got {complete}"));
    }
    Ok(())
}

fn build(cli: &Cli) -> Result<(Manifest, Vec<u8>), String> {
    if size_of::<Request>() != 16 || size_of::<Accumulator>() != 16 {
        return Err("registered record layout changed".into());
    }
    let (round42, value42) = load_round42(&cli.round42)?;
    let (isolation, isolation_records) = load_isolation(&cli.isolation)?;
    let summaries42 = json_at(&value42, &["timing", "summaries"])?
        .as_array()
        .ok_or("round-42 timing summaries are not an array")?;
    let summary_value = |name: &str| -> Result<f64, String> {
        summaries42
            .iter()
            .find(|row| row.get("variant").and_then(Value::as_str) == Some(name))
            .and_then(|row| row.get("addition_equivalents_per_start"))
            .and_then(Value::as_f64)
            .ok_or_else(|| format!("round-42 summary {name} missing"))
    };
    let imported = ImportedEvidence {
        round42_request_digest: json_str(&value42, &["imported", "round41_request_digest"])?.into(),
        round42_native_control_digest: json_str(
            &value42,
            &["imported", "round41_native_control_digest"],
        )?
        .into(),
        round42_zero_false_positives_and_false_negatives: json_bool(
            &value42,
            &["gates", "zero_false_positives_and_false_negatives"],
        )?,
        round42_exact_reference_equality: json_bool(
            &value42,
            &["gates", "exact_equality_with_materialized_reference"],
        )?,
        round42_full_universe_algorithmic_bytes: json_u64(
            &value42,
            &["frozen_model", "full_universe_algorithmic_bytes"],
        )?,
        round42_materialized_complete_addition_equivalents_per_start: summary_value(
            "complete-materialized-3x12",
        )?,
        round42_direct_complete_addition_equivalents_per_start: summary_value(
            "complete-direct-address-replay",
        )?,
        round42_promoted: json_bool(&value42, &["gates", "promoted"])?,
        round42_isolated_runs: isolation_records.len() as u64,
        round42_last_run_uncontended: json_at(
            isolation_records.last().expect("nonempty isolation"),
            &["run", "contended"],
        )?
        .as_bool()
        .is_some_and(|contended| !contended),
    };
    if imported.round42_request_digest != REQUEST_DIGEST
        || imported.round42_native_control_digest != CONTROL_DIGEST
        || !imported.round42_zero_false_positives_and_false_negatives
        || !imported.round42_exact_reference_equality
        || imported.round42_full_universe_algorithmic_bytes != 557_314_646_392
        || imported.round42_promoted
        || !imported.round42_last_run_uncontended
    {
        return Err("round-42 frozen evidence changed".into());
    }

    for width in [6, 12, 18, 24] {
        let sample = generate_requests(64);
        let _ = timed_route(&sample, 64, width);
    }
    let (aa_first, _) = timed_additions(BASELINE_ADDITIONS);
    let (aa_second, _) = timed_additions(BASELINE_ADDITIONS);
    let baseline_values = (0..TIMING_REPETITIONS)
        .map(|_| timed_additions(BASELINE_ADDITIONS).0)
        .collect::<Vec<_>>();
    let baseline_ns_per_addition = median_u64(&baseline_values) / BASELINE_ADDITIONS as f64;

    let mut records = Vec::with_capacity(PORTFOLIO_CELLS);
    for (depth_index, start_bits) in START_BITS.into_iter().enumerate() {
        if start_bits <= MAX_EXECUTED_START_BITS {
            let context = exact_context(start_bits);
            for width in MIN_WIDTH..=MAX_WIDTH {
                if width <= MAX_EXECUTED_WIDTH {
                    records.push(exact_record(
                        depth_index,
                        width,
                        &context,
                        baseline_ns_per_addition,
                    ));
                } else {
                    records.push(projection_record(depth_index, width));
                }
            }
        } else {
            for width in MIN_WIDTH..=MAX_WIDTH {
                records.push(projection_record(depth_index, width));
            }
        }
    }
    records.sort_by_key(|record| record.round);
    validate_portfolio(&records)?;

    let (shortlist_widths, raw_shortlist, shortlist_summaries) =
        shortlist_timing(&records, baseline_ns_per_addition);
    let (time_exponents, memory_exponents) =
        apply_exponents_and_shortlist(&mut records, &shortlist_summaries);
    let entries = write_rounds(&cli.round_dir, &records)?;

    let complete_checked_rounds = records
        .iter()
        .filter(|record| record.complete_checked)
        .map(|record| record.round)
        .collect::<Vec<_>>();
    let projection_only_rounds = records
        .iter()
        .filter(|record| !record.complete_checked)
        .map(|record| record.round)
        .collect::<Vec<_>>();
    let failures = records
        .iter()
        .filter(|record| record.complete_checked)
        .map(|record| record.false_positives.unwrap_or(1) + record.false_negatives.unwrap_or(1))
        .sum::<u64>();
    let best_measured = shortlist_summaries
        .iter()
        .min_by(|left, right| {
            left.1
                .median_wall_ns
                .total_cmp(&right.1.median_wall_ns)
                .then(left.0.cmp(&right.0))
        })
        .ok_or("empty shortlist")?;
    let best_measured_width = best_measured.0;
    let best_measured_round = round_number(2, best_measured_width);
    let best_formula = records
        .iter()
        .filter(|record| record.start_bits == 40)
        .min_by(|left, right| {
            let left_bytes = left
                .traffic
                .total_logical_bytes
                .parse::<u128>()
                .expect("traffic integer");
            let right_bytes = right
                .traffic
                .total_logical_bytes
                .parse::<u128>()
                .expect("traffic integer");
            left_bytes
                .cmp(&right_bytes)
                .then(
                    left.peak_algorithmic_bytes
                        .cmp(&right.peak_algorithmic_bytes),
                )
                .then(left.radix_width.cmp(&right.radix_width))
        })
        .ok_or("missing depth-40 rows")?;
    let best_formula_round = best_formula.round;
    let best_formula_width = best_formula.radix_width;
    let timing = TimingEvidence {
        baseline_additions: BASELINE_ADDITIONS,
        aa_first_ns: aa_first,
        aa_second_ns: aa_second,
        aa_ratio: aa_first.max(aa_second) as f64 / aa_first.min(aa_second) as f64,
        baseline_repetitions: TIMING_REPETITIONS,
        baseline_median_ns_per_addition: baseline_ns_per_addition,
        shortlist_depth_bits: 16,
        shortlist_widths,
        shortlist_repetitions: TIMING_REPETITIONS,
        raw_shortlist,
        summaries: shortlist_summaries,
        time_exponents,
        memory_exponents,
        scope: "One-thread RAM routing kernels. All timings exclude shared request generation; projected cells have no timings and do not inherit exactness."
            .into(),
    };
    let gates = PortfolioGates {
        all_248_rounds_present_once: entries.len() == PORTFOLIO_CELLS,
        all_76_registered_cells_complete: complete_checked_rounds.len() == 76,
        zero_false_positives_and_false_negatives_on_complete_cells: failures == 0,
        structured_residual_degree_at_most_5_for_measured_family: false,
        relation_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        at_least_one_storage_projection_below_2_50: records
            .iter()
            .any(|record| record.gates.materialized_storage_below_2_50),
        measured_complete_selector_below_rho_in_applicable_tier: false,
        no_projection_counted_as_exhaustive: records
            .iter()
            .filter(|record| !record.complete_checked)
            .all(|record| {
                record.false_positives.is_none()
                    && record.false_negatives.is_none()
                    && record.one_shot_wall_ns.is_none()
                    && !record.applicable_memory_tier_measured
            }),
        promoted_rounds: Vec::new(),
    };
    let semantic = serde_json::json!({
        "curve": CURVE_SLUG,
        "round42_sha256": ROUND42_SHA256,
        "round42_isolation_sha256": ROUND42_ISOLATION_SHA256,
        "imported": &imported,
        "timing": &timing,
        "entries": &entries,
        "gates": &gates,
        "best_measured_round": best_measured_round,
        "best_formula_round": best_formula_round,
    });
    let semantic_evidence_sha256 = hex::encode(sha256(
        &serde_json::to_vec(&semantic).map_err(|error| error.to_string())?,
    ));
    let manifest = Manifest {
        schema: "p256-radix-portfolio-round43-290-v1".into(),
        curve: CURVE_SLUG.into(),
        round_range: [FIRST_ROUND, LAST_ROUND],
        family: "round-33 known-log selector radix routing".into(),
        dependencies: vec![round42, isolation],
        imported,
        frozen_matrix: serde_json::json!({
            "radix_widths": [MIN_WIDTH, MAX_WIDTH],
            "start_bits": START_BITS,
            "round_formula": "43 + 31 * depth_index + (radix_width - 6)",
            "pair_key_bits": KEY_BITS,
            "pair_universe": PAIR_UNIVERSE,
            "requests_per_start": REQUESTS_PER_START,
            "record_bytes": RECORD_BYTES,
            "accumulator_bytes": ACCUMULATOR_BYTES,
            "complete_execution": {
                "maximum_start_bits": MAX_EXECUTED_START_BITS,
                "maximum_radix_width": MAX_EXECUTED_WIDTH,
                "cells": 76
            },
            "projection_only_cells": 172,
            "comparison_factor_base": COMPARISON_FB,
            "comparison_structured_residual_maximum": 4,
            "measured_family_structured_residual_degree": null,
        }),
        timing,
        entries,
        complete_checked_rounds,
        projection_only_rounds,
        best_measured_kernel_round: best_measured_round,
        best_measured_kernel_width: best_measured_width,
        best_formula_traffic_round: best_formula_round,
        best_formula_traffic_width: best_formula_width,
        gates,
        semantic_evidence_sha256,
        full_depth_unplanted_attempted: false,
        interpretation: "Rounds 43–290 are one preregistered finite portfolio, not 248 adaptive claims. Complete cells establish exact routing behavior only; formula-only cells preserve exact resource counts without claiming execution."
            .into(),
        decision: "No round is promoted. The portfolio may select a routing-kernel width, but the measured family still lacks structured-degree and relation-collection evidence, and no complete selector is measured below rho in the applicable tier."
            .into(),
    };
    let bytes = serialize_pretty(&manifest)?;
    Ok((manifest, bytes))
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    if let Some(path) = &cli.out {
        if path.exists() {
            eprintln!("refusing to overwrite {}", path.display());
            return ExitCode::FAILURE;
        }
    }
    match build(&cli) {
        Ok((_manifest, bytes)) => {
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
            eprintln!("P-256 radix portfolio failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_map_is_complete_and_exact() {
        assert_eq!(round_number(0, 6), 43);
        assert_eq!(round_number(7, 36), 290);
        let rounds = START_BITS
            .iter()
            .enumerate()
            .flat_map(|(index, _)| {
                (MIN_WIDTH..=MAX_WIDTH).map(move |width| round_number(index, width))
            })
            .collect::<BTreeSet<_>>();
        assert_eq!(rounds.len(), PORTFOLIO_CELLS);
        assert_eq!(rounds.first(), Some(&FIRST_ROUND));
        assert_eq!(rounds.last(), Some(&LAST_ROUND));
    }

    #[test]
    fn pass_layouts_cover_exactly_36_bits() {
        assert_eq!(histogram_bits(6), vec![6, 6, 6, 6, 6, 6]);
        assert_eq!(histogram_bits(12), vec![12, 12, 12]);
        assert_eq!(histogram_bits(18), vec![18, 18]);
        assert_eq!(histogram_bits(20), vec![20, 16]);
        assert_eq!(histogram_bits(36), vec![36]);
        for width in MIN_WIDTH..=MAX_WIDTH {
            assert_eq!(histogram_bits(width).iter().sum::<u32>(), KEY_BITS);
        }
    }

    #[test]
    fn radix_widths_match_canonical_toy() {
        let starts = 128;
        let frozen = generate_requests(starts);
        let mut canonical = frozen.clone();
        canonical.sort_unstable();
        for width in MIN_WIDTH..=12 {
            assert_eq!(radix_sort(frozen.clone(), width), canonical);
        }
    }

    #[test]
    fn registered_layouts_and_cell_counts_hold() {
        assert_eq!(size_of::<Request>(), 16);
        assert_eq!(size_of::<Accumulator>(), 16);
        assert_eq!((MAX_EXECUTED_WIDTH - MIN_WIDTH + 1) * 4, 76);
        assert_eq!((MAX_WIDTH - MIN_WIDTH + 1) as usize * START_BITS.len(), 248);
    }
}
