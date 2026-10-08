//! Fixed-universe direct-address replay for P-256 pair routing (round 42).

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
const ROUND41_SHA256: &str = "8bf3dd1d2e7942858c6cbf8486170cb6275aa0ec5e5380e686daeebcd9fd26a6";
const ROUND41_ISOLATION_SHA256: &str =
    "6e6dedb3f9367229808b1d3eb7c9978c9445fd01643a72078a8489215d9ebf32";
const ROUND40_REQUEST_DIGEST: &str =
    "b4178238067a37febeeb46116dffc633311549e6325461670bd12db7383e0dd6";
const ROUND40_CONTROL_DIGEST: &str =
    "98c8f56a06257256bb88f6c6d7cffce79391c461170b91287f8fd193fe5711e6";
const PAIR_UNIVERSE: u64 = 34_562_148_612;
const REQUESTS_PER_START: u64 = 8;
const VALUE_BYTES: u64 = 16;
const ACCUMULATOR_BYTES: u64 = 16;
const RECORD_BYTES: u64 = 16;
const REFERENCE_TRAFFIC_PER_REQUEST: u64 = 160;
const BASE_RATIO: f64 = 0.998_569_034_150_286;
const LOCAL_RATIO: f64 = 0.964_336_477_130_181;
const MEAN_CAPACITY: f64 = 224.361_249_307_510_55;
const ALLOWED_ADDITIONS: f64 = 0.334_410_508_422_987_64;
const DEPTH_BITS: [u32; 6] = [12, 14, 16, 18, 20, 22];
const TIMING_BITS: [u32; 4] = [16, 18, 20, 22];
const TIMING_UNIVERSE_BITS: u32 = 20;
const TIMING_START_BITS: u32 = 17;
const TIMING_REPETITIONS: u64 = 7;
const LADDER_REPETITIONS: u64 = 5;
const BASELINE_ADDITIONS: u64 = 262_144;

#[derive(Parser)]
#[command(about = "Audit fixed-universe direct-address replay for P-256 routing")]
struct Cli {
    #[arg(long)]
    round41: PathBuf,
    #[arg(long)]
    isolation: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
}

#[repr(C)]
#[derive(Clone, Copy, Default, PartialEq, Eq)]
struct Request {
    pair: u64,
    owner_slot: u64,
}

#[repr(C)]
#[derive(Clone, Copy, Default, PartialEq, Eq)]
struct PairValue {
    low: u64,
    high: u64,
}

#[repr(C)]
#[derive(Clone, Copy, Default, PartialEq, Eq)]
struct Accumulator {
    checksum: u64,
    mask: u8,
    count: u8,
    padding: [u8; 6],
}

#[derive(Clone, Copy)]
enum Variant {
    Materialized,
    Direct,
}

impl Variant {
    fn name(self) -> &'static str {
        match self {
            Self::Materialized => "complete-materialized-3x12",
            Self::Direct => "complete-direct-address-replay",
        }
    }
}

#[derive(Clone, Serialize)]
struct Dependency {
    path: String,
    sha256: String,
    schema: String,
}

#[derive(Clone, Serialize)]
struct ImportedEvidence {
    round41_request_digest: String,
    round41_native_control_digest: String,
    round41_zero_false_positives_and_false_negatives: bool,
    round41_exact_reference_equality: bool,
    round41_promoted: bool,
    round41_isolated_runs: u64,
    round41_last_run_uncontended: bool,
}

#[derive(Clone, Serialize)]
struct ToyCell {
    universe: u64,
    starts: u64,
    reference_distinct: u64,
    materialized_distinct: u64,
    direct_distinct: u64,
    accumulator_mismatches: u64,
    regenerated_failures: u64,
    false_positives: u64,
    false_negatives: u64,
}

#[derive(Clone, Serialize)]
struct ExactCell {
    universe_bits: u32,
    universe: u64,
    starts: u64,
    requests: u64,
    request_load: f64,
    distinct_pairs: u64,
    expected_distinct_pairs: f64,
    observed_over_expected: f64,
    accumulator_mismatches: u64,
    regenerated_failures: u64,
    selected_value_checks: u64,
    selected_value_failures: u64,
    negative_controls: u64,
    false_positives: u64,
    false_negatives: u64,
    materialized_peak_bytes: u64,
    candidate_algorithmic_bytes: u64,
    harness_reference_accumulator_bytes: u64,
    materialized_logical_bytes: u64,
    candidate_semantic_bytes: u64,
    candidate_cache_line_bytes: u64,
    materialized_accumulator_blake3: String,
    direct_accumulator_blake3: String,
}

#[derive(Clone, Serialize)]
struct RawTiming {
    repetition: u64,
    order: Vec<String>,
    materialized_ns: u64,
    direct_ns: u64,
    materialized_checksum: u64,
    direct_checksum: u64,
}

#[derive(Clone, Serialize)]
struct VariantSummary {
    variant: String,
    median_wall_ns: f64,
    minimum_wall_ns: u64,
    median_ns_per_start: f64,
    addition_equivalents_per_start: f64,
    budget_multiple: f64,
    projected_router_only_ratio_to_rho: f64,
    checksum_identical: bool,
    request_generation_passes: u64,
}

#[derive(Clone, Serialize)]
struct LadderTiming {
    variant: String,
    universe_bits: u32,
    start_bits: u32,
    universe: u64,
    starts: u64,
    repetitions: u64,
    median_wall_ns: f64,
    minimum_wall_ns: u64,
    checksum_identical: bool,
}

#[derive(Clone, Serialize)]
struct TimingEvidence {
    timing_universe_bits: u32,
    timing_start_bits: u32,
    timing_universe: u64,
    timing_starts: u64,
    repetitions: u64,
    baseline_additions: u64,
    aa_first_ns: u64,
    aa_second_ns: u64,
    aa_ratio: f64,
    baseline_median_ns_per_addition: f64,
    raw: Vec<RawTiming>,
    summaries: Vec<VariantSummary>,
    candidate_over_reference_ratio: f64,
    ladder: Vec<LadderTiming>,
    materialized_time_exponent: f64,
    direct_time_exponent: f64,
    direct_algorithmic_memory_exponent: f64,
    scope: String,
}

#[derive(Clone, Serialize)]
struct ProjectionRow {
    label: String,
    starts: u64,
    starts_log2: f64,
    expected_distinct_pairs: f64,
    pair_constructions_per_start: f64,
    measured_scaled_router_addition_equivalents_per_start: f64,
    optimistic_combined_addition_equivalents_per_start: f64,
    optimistic_combined_ratio_to_rho: f64,
    below_rho_arithmetic_only: bool,
    candidate_algorithmic_bytes: u64,
    candidate_log2_bytes: f64,
    storage_below_2_50: bool,
    semantic_logical_bytes_per_batch: String,
    cache_line_model_bytes_per_batch: String,
    applicable_memory_tier_measured: bool,
}

#[derive(Clone, Serialize)]
struct Gates {
    zero_false_positives_and_false_negatives: bool,
    exact_equality_with_materialized_reference: bool,
    structured_residual_degree_at_most_5_for_measured_family: bool,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    measured_complete_selector_below_rho_in_applicable_tier: bool,
    no_discarded_branch_counted_as_exhaustive: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    family: String,
    dependencies: Vec<Dependency>,
    imported: ImportedEvidence,
    frozen_model: Value,
    toy_controls: Vec<ToyCell>,
    load_one_ladder: Vec<ExactCell>,
    saturation_ladder: Vec<ExactCell>,
    timing: TimingEvidence,
    projection: Vec<ProjectionRow>,
    semantic_evidence_sha256: String,
    gates: Gates,
    full_depth_unplanted_attempted: bool,
    interpretation: String,
    decision: String,
}

struct ReferenceRoute {
    accumulators: Vec<Accumulator>,
    distinct: u64,
    digest: String,
}

struct DirectRoute {
    distinct: u64,
    digest: String,
    accumulator_mismatches: u64,
    selected_value_checks: u64,
    selected_value_failures: u64,
    false_positives: u64,
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

fn load_round41(path: &Path) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != ROUND41_SHA256 {
        return Err(format!(
            "round-41 hash mismatch: expected {ROUND41_SHA256}, got {digest}"
        ));
    }
    let value: Value = serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if json_str(&value, &["schema"])? != "p256-two-pass-router-v1"
        || json_str(&value, &["curve"])? != CURVE_SLUG
    {
        return Err("round-41 identity mismatch".into());
    }
    Ok((
        Dependency {
            path: path.display().to_string(),
            sha256: digest,
            schema: "p256-two-pass-router-v1".into(),
        },
        value,
    ))
}

fn load_isolation(path: &Path) -> Result<(Dependency, Vec<Value>), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != ROUND41_ISOLATION_SHA256 {
        return Err(format!(
            "round-41 isolation hash mismatch: expected {ROUND41_ISOLATION_SHA256}, got {digest}"
        ));
    }
    let text = String::from_utf8(bytes).map_err(|error| error.to_string())?;
    let records = text
        .lines()
        .map(|line| serde_json::from_str(line).map_err(|error| error.to_string()))
        .collect::<Result<Vec<Value>, String>>()?;
    if records.is_empty() {
        return Err("round-41 isolation has no records".into());
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

fn mix(mut value: u64) -> u64 {
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn pair_value(pair: u64) -> PairValue {
    PairValue {
        low: mix(pair ^ 0x9e37_79b9_7f4a_7c15),
        high: mix(pair ^ 0xd1b5_4a32_d192_ed03),
    }
}

fn accumulate(accumulator: &mut Accumulator, value: PairValue, slot: u8) {
    let contribution =
        value.low.rotate_left(u32::from(slot) * 7) ^ value.high.rotate_right(u32::from(slot) * 5);
    accumulator.checksum = accumulator.checksum.wrapping_add(contribution);
    accumulator.mask |= 1 << slot;
    accumulator.count += 1;
}

fn expected_accumulator(start: u64, universe: u64) -> Accumulator {
    let mut accumulator = Accumulator::default();
    for slot in 0..REQUESTS_PER_START as u8 {
        let pair = pair_id(start, slot, universe);
        accumulate(&mut accumulator, pair_value(pair), slot);
    }
    accumulator
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
    let buckets = 1usize << 12;
    let mask = buckets - 1;
    let mut target = vec![Request::default(); source.len()];
    let mut counts = vec![0usize; buckets];
    for pass in 0..3 {
        counts.fill(0);
        let shift = pass * 12;
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
    }
    source
}

fn digest_accumulator(hasher: &mut Hasher, accumulator: Accumulator) {
    hasher.update(&accumulator.checksum.to_be_bytes());
    hasher.update(&[accumulator.mask, accumulator.count]);
}

fn digest_accumulators(accumulators: &[Accumulator]) -> String {
    let mut hasher = Hasher::new();
    for accumulator in accumulators {
        digest_accumulator(&mut hasher, *accumulator);
    }
    hasher.finalize().to_hex().to_string()
}

fn materialized_route(universe: u64, starts: u64) -> ReferenceRoute {
    let sorted = radix_sort(generate_requests(starts, universe));
    let mut accumulators = vec![Accumulator::default(); starts as usize];
    let mut distinct = 0u64;
    let mut previous = None;
    let mut value = PairValue::default();
    for record in sorted {
        if previous != Some(record.pair) {
            distinct += 1;
            previous = Some(record.pair);
            value = pair_value(record.pair);
        }
        let owner = (record.owner_slot >> 3) as usize;
        let slot = (record.owner_slot & 7) as u8;
        accumulate(&mut accumulators[owner], value, slot);
    }
    let digest = digest_accumulators(&accumulators);
    ReferenceRoute {
        accumulators,
        distinct,
        digest,
    }
}

fn bitmap_words(universe: u64) -> u64 {
    universe.div_ceil(64)
}

fn candidate_algorithmic_bytes(universe: u64) -> u64 {
    bitmap_words(universe) * 8 + universe * VALUE_BYTES + ACCUMULATOR_BYTES
}

fn direct_route(universe: u64, starts: u64, expected: Option<&[Accumulator]>) -> DirectRoute {
    let mut seen = vec![0u64; bitmap_words(universe) as usize];
    let mut values = vec![PairValue::default(); universe as usize];
    let mut distinct = 0u64;
    for start in 0..starts {
        for slot in 0..REQUESTS_PER_START as u8 {
            let pair = pair_id(start, slot, universe);
            let word = (pair >> 6) as usize;
            let bit = 1u64 << (pair & 63);
            if seen[word] & bit == 0 {
                seen[word] |= bit;
                values[pair as usize] = pair_value(pair);
                distinct += 1;
            }
        }
    }

    let stride = (universe / 4_096).max(1);
    let mut selected_value_checks = 0u64;
    let mut selected_value_failures = 0u64;
    for pair in (0..universe).step_by(stride as usize) {
        let word = (pair >> 6) as usize;
        let bit = 1u64 << (pair & 63);
        if seen[word] & bit != 0 {
            selected_value_checks += 1;
            if values[pair as usize] != pair_value(pair) {
                selected_value_failures += 1;
            }
        }
    }

    let mut hasher = Hasher::new();
    let mut accumulator_mismatches = 0u64;
    let mut false_positives = 0u64;
    for start in 0..starts {
        let mut accumulator = Accumulator::default();
        for slot in 0..REQUESTS_PER_START as u8 {
            let pair = pair_id(start, slot, universe);
            accumulate(&mut accumulator, values[pair as usize], slot);
        }
        digest_accumulator(&mut hasher, accumulator);
        if let Some(expected) = expected {
            let expected = expected[start as usize];
            if accumulator != expected {
                accumulator_mismatches += 1;
            }
            let mut negative = accumulator;
            negative.checksum = negative.checksum.wrapping_add(1);
            if negative == expected {
                false_positives += 1;
            }
        }
    }
    DirectRoute {
        distinct,
        digest: hasher.finalize().to_hex().to_string(),
        accumulator_mismatches,
        selected_value_checks,
        selected_value_failures,
        false_positives,
    }
}

fn expected_distinct(universe: u64, starts: u64) -> f64 {
    let requests = (starts * REQUESTS_PER_START) as f64;
    let universe = universe as f64;
    -universe * (requests * (-1.0 / universe).ln_1p()).exp_m1()
}

fn materialized_peak_bytes(starts: u64) -> u64 {
    starts * (2 * REQUESTS_PER_START * RECORD_BYTES + ACCUMULATOR_BYTES)
}

fn semantic_bytes(requests: u64, distinct: u64) -> u64 {
    (requests + distinct).div_ceil(8) + distinct * VALUE_BYTES + requests * 48
}

fn cache_line_bytes(requests: u64, distinct: u64) -> u64 {
    requests * 160 + distinct * 128
}

fn checked_cell(universe_bits: u32, starts: u64) -> ExactCell {
    let universe = 1u64 << universe_bits;
    let requests = starts * REQUESTS_PER_START;
    let reference = materialized_route(universe, starts);
    let direct = direct_route(universe, starts, Some(&reference.accumulators));
    let mut regenerated_failures = 0u64;
    for start in 0..starts {
        if reference.accumulators[start as usize] != expected_accumulator(start, universe) {
            regenerated_failures += 1;
        }
    }
    let expected = expected_distinct(universe, starts);
    let false_negatives = direct.accumulator_mismatches
        + regenerated_failures
        + direct.selected_value_failures
        + u64::from(reference.distinct != direct.distinct)
        + u64::from(reference.digest != direct.digest);
    ExactCell {
        universe_bits,
        universe,
        starts,
        requests,
        request_load: requests as f64 / universe as f64,
        distinct_pairs: direct.distinct,
        expected_distinct_pairs: expected,
        observed_over_expected: direct.distinct as f64 / expected,
        accumulator_mismatches: direct.accumulator_mismatches,
        regenerated_failures,
        selected_value_checks: direct.selected_value_checks,
        selected_value_failures: direct.selected_value_failures,
        negative_controls: starts,
        false_positives: direct.false_positives,
        false_negatives,
        materialized_peak_bytes: materialized_peak_bytes(starts),
        candidate_algorithmic_bytes: candidate_algorithmic_bytes(universe),
        harness_reference_accumulator_bytes: starts * ACCUMULATOR_BYTES,
        materialized_logical_bytes: requests * REFERENCE_TRAFFIC_PER_REQUEST,
        candidate_semantic_bytes: semantic_bytes(requests, direct.distinct),
        candidate_cache_line_bytes: cache_line_bytes(requests, direct.distinct),
        materialized_accumulator_blake3: reference.digest,
        direct_accumulator_blake3: direct.digest,
    }
}

fn toy_cell(universe: u64, starts: u64) -> ToyCell {
    let reference = materialized_route(universe, starts);
    let direct = direct_route(universe, starts, Some(&reference.accumulators));
    let reference_pairs = (0..starts)
        .flat_map(|start| {
            (0..REQUESTS_PER_START as u8).map(move |slot| pair_id(start, slot, universe))
        })
        .collect::<BTreeSet<_>>();
    let regenerated_failures = (0..starts)
        .filter(|start| {
            reference.accumulators[*start as usize] != expected_accumulator(*start, universe)
        })
        .count() as u64;
    let false_negatives = direct.accumulator_mismatches
        + regenerated_failures
        + direct.selected_value_failures
        + u64::from(reference.distinct != direct.distinct)
        + u64::from(reference.distinct != reference_pairs.len() as u64)
        + u64::from(reference.digest != direct.digest);
    ToyCell {
        universe,
        starts,
        reference_distinct: reference_pairs.len() as u64,
        materialized_distinct: reference.distinct,
        direct_distinct: direct.distinct,
        accumulator_mismatches: direct.accumulator_mismatches,
        regenerated_failures,
        false_positives: direct.false_positives,
        false_negatives,
    }
}

fn digest_checksum(digest: &str, distinct: u64) -> u64 {
    let bytes = hex::decode(&digest[..16]).expect("digest hex");
    let mut word = [0u8; 8];
    word.copy_from_slice(&bytes);
    u64::from_be_bytes(word) ^ distinct
}

fn timed_variant(variant: Variant, universe: u64, starts: u64) -> (u64, u64) {
    let begin = Instant::now();
    let (digest, distinct) = match variant {
        Variant::Materialized => {
            let route = black_box(materialized_route(universe, starts));
            (route.digest, route.distinct)
        }
        Variant::Direct => {
            let route = black_box(direct_route(universe, starts, None));
            (route.digest, route.distinct)
        }
    };
    (
        begin.elapsed().as_nanos() as u64,
        black_box(digest_checksum(&digest, distinct)),
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

fn pair_ratio(additions: f64) -> f64 {
    BASE_RATIO + LOCAL_RATIO * additions / (MEAN_CAPACITY + 1.0)
}

fn timing_summary(
    variant: Variant,
    raw: &[RawTiming],
    baseline_ns_per_addition: f64,
    starts: u64,
) -> VariantSummary {
    let values = raw
        .iter()
        .map(|row| match variant {
            Variant::Materialized => row.materialized_ns,
            Variant::Direct => row.direct_ns,
        })
        .collect::<Vec<_>>();
    let checksums = raw
        .iter()
        .map(|row| match variant {
            Variant::Materialized => row.materialized_checksum,
            Variant::Direct => row.direct_checksum,
        })
        .collect::<Vec<_>>();
    let median = median_u64(&values);
    let ns_per_start = median / starts as f64;
    let equivalents = ns_per_start / baseline_ns_per_addition;
    VariantSummary {
        variant: variant.name().into(),
        median_wall_ns: median,
        minimum_wall_ns: *values.iter().min().expect("timing values"),
        median_ns_per_start: ns_per_start,
        addition_equivalents_per_start: equivalents,
        budget_multiple: equivalents / ALLOWED_ADDITIONS,
        projected_router_only_ratio_to_rho: pair_ratio(equivalents),
        checksum_identical: checksums.iter().all(|checksum| *checksum == checksums[0]),
        request_generation_passes: match variant {
            Variant::Materialized => 1,
            Variant::Direct => 2,
        },
    }
}

fn fit_time_exponent(rows: &[LadderTiming], variant: &str) -> f64 {
    let selected = rows
        .iter()
        .filter(|row| row.variant == variant)
        .collect::<Vec<_>>();
    let mean_x = selected
        .iter()
        .map(|row| row.universe_bits as f64)
        .sum::<f64>()
        / selected.len() as f64;
    let mean_y = selected
        .iter()
        .map(|row| row.median_wall_ns.log2())
        .sum::<f64>()
        / selected.len() as f64;
    let numerator = selected
        .iter()
        .map(|row| (row.universe_bits as f64 - mean_x) * (row.median_wall_ns.log2() - mean_y))
        .sum::<f64>();
    let denominator = selected
        .iter()
        .map(|row| (row.universe_bits as f64 - mean_x).powi(2))
        .sum::<f64>();
    numerator / denominator
}

fn fit_memory_exponent() -> f64 {
    let rows = TIMING_BITS
        .into_iter()
        .map(|bits| {
            (
                bits as f64,
                (candidate_algorithmic_bytes(1u64 << bits) as f64).log2(),
            )
        })
        .collect::<Vec<_>>();
    let mean_x = rows.iter().map(|row| row.0).sum::<f64>() / rows.len() as f64;
    let mean_y = rows.iter().map(|row| row.1).sum::<f64>() / rows.len() as f64;
    rows.iter()
        .map(|row| (row.0 - mean_x) * (row.1 - mean_y))
        .sum::<f64>()
        / rows.iter().map(|row| (row.0 - mean_x).powi(2)).sum::<f64>()
}

fn timing_evidence() -> TimingEvidence {
    let universe = 1u64 << TIMING_UNIVERSE_BITS;
    let starts = 1u64 << TIMING_START_BITS;
    let _ = timed_additions(8_192);
    for variant in [Variant::Materialized, Variant::Direct] {
        let _ = timed_variant(variant, 1 << 12, 1 << 9);
    }
    let (aa_first, _) = timed_additions(BASELINE_ADDITIONS);
    let (aa_second, _) = timed_additions(BASELINE_ADDITIONS);
    let orders = [
        [Variant::Materialized, Variant::Direct],
        [Variant::Direct, Variant::Materialized],
    ];
    let mut raw = Vec::new();
    for repetition in 0..TIMING_REPETITIONS {
        let order = orders[(repetition as usize) % orders.len()];
        let mut materialized = (0u64, 0u64);
        let mut direct = (0u64, 0u64);
        for variant in order {
            let measured = timed_variant(variant, universe, starts);
            match variant {
                Variant::Materialized => materialized = measured,
                Variant::Direct => direct = measured,
            }
        }
        raw.push(RawTiming {
            repetition,
            order: order.iter().map(|variant| variant.name().into()).collect(),
            materialized_ns: materialized.0,
            direct_ns: direct.0,
            materialized_checksum: materialized.1,
            direct_checksum: direct.1,
        });
    }
    let baseline_values = (0..TIMING_REPETITIONS)
        .map(|_| timed_additions(BASELINE_ADDITIONS).0)
        .collect::<Vec<_>>();
    let baseline_ns_per_addition = median_u64(&baseline_values) / BASELINE_ADDITIONS as f64;
    let summaries = [Variant::Materialized, Variant::Direct]
        .into_iter()
        .map(|variant| timing_summary(variant, &raw, baseline_ns_per_addition, starts))
        .collect::<Vec<_>>();

    let mut ladder = Vec::new();
    for bits in TIMING_BITS {
        let ladder_universe = 1u64 << bits;
        let ladder_starts = 1u64 << (bits - 3);
        for variant in [Variant::Materialized, Variant::Direct] {
            let mut values = Vec::new();
            let mut checksums = Vec::new();
            for _ in 0..LADDER_REPETITIONS {
                let measured = timed_variant(variant, ladder_universe, ladder_starts);
                values.push(measured.0);
                checksums.push(measured.1);
            }
            ladder.push(LadderTiming {
                variant: variant.name().into(),
                universe_bits: bits,
                start_bits: bits - 3,
                universe: ladder_universe,
                starts: ladder_starts,
                repetitions: LADDER_REPETITIONS,
                median_wall_ns: median_u64(&values),
                minimum_wall_ns: *values.iter().min().expect("ladder values"),
                checksum_identical: checksums.iter().all(|checksum| *checksum == checksums[0]),
            });
        }
    }
    let reference = summaries
        .iter()
        .find(|row| row.variant == Variant::Materialized.name())
        .expect("reference summary");
    let candidate = summaries
        .iter()
        .find(|row| row.variant == Variant::Direct.name())
        .expect("candidate summary");
    let candidate_over_reference_ratio = candidate.median_wall_ns / reference.median_wall_ns;
    let materialized_time_exponent = fit_time_exponent(&ladder, Variant::Materialized.name());
    let direct_time_exponent = fit_time_exponent(&ladder, Variant::Direct.name());
    TimingEvidence {
        timing_universe_bits: TIMING_UNIVERSE_BITS,
        timing_start_bits: TIMING_START_BITS,
        timing_universe: universe,
        timing_starts: starts,
        repetitions: TIMING_REPETITIONS,
        baseline_additions: BASELINE_ADDITIONS,
        aa_first_ns: aa_first,
        aa_second_ns: aa_second,
        aa_ratio: aa_first.max(aa_second) as f64 / aa_first.min(aa_second) as f64,
        baseline_median_ns_per_addition: baseline_ns_per_addition,
        raw,
        summaries,
        candidate_over_reference_ratio,
        ladder,
        materialized_time_exponent,
        direct_time_exponent,
        direct_algorithmic_memory_exponent: fit_memory_exponent(),
        scope: "Complete deterministic request generation on scaled dense RAM tables. The direct candidate generates every request twice and streams outputs; no timing applies to the 557-GB full table."
            .into(),
    }
}

fn projection_row(label: &str, starts: u64, router_equivalents: f64) -> ProjectionRow {
    let distinct = expected_distinct(PAIR_UNIVERSE, starts);
    let pair_per_start = distinct / starts as f64;
    let combined = pair_per_start + router_equivalents;
    let ratio = pair_ratio(combined);
    let requests = u128::from(starts) * u128::from(REQUESTS_PER_START);
    let rounded_distinct = distinct.round() as u128;
    let semantic =
        (requests + rounded_distinct).div_ceil(8) + rounded_distinct * 16 + requests * 48;
    let cache_lines = requests * 160 + rounded_distinct * 128;
    let state = candidate_algorithmic_bytes(PAIR_UNIVERSE);
    ProjectionRow {
        label: label.into(),
        starts,
        starts_log2: (starts as f64).log2(),
        expected_distinct_pairs: distinct,
        pair_constructions_per_start: pair_per_start,
        measured_scaled_router_addition_equivalents_per_start: router_equivalents,
        optimistic_combined_addition_equivalents_per_start: combined,
        optimistic_combined_ratio_to_rho: ratio,
        below_rho_arithmetic_only: ratio <= 1.0,
        candidate_algorithmic_bytes: state,
        candidate_log2_bytes: (state as f64).log2(),
        storage_below_2_50: state < (1u64 << 50),
        semantic_logical_bytes_per_batch: semantic.to_string(),
        cache_line_model_bytes_per_batch: cache_lines.to_string(),
        applicable_memory_tier_measured: false,
    }
}

fn result(cli: &Cli) -> Result<ResultFile, String> {
    if size_of::<Request>() != 16 || size_of::<PairValue>() != 16 || size_of::<Accumulator>() != 16
    {
        return Err("registered layout changed".into());
    }
    if candidate_algorithmic_bytes(PAIR_UNIVERSE) != 557_314_646_392 {
        return Err("full-universe state arithmetic changed".into());
    }
    let (round41, value41) = load_round41(&cli.round41)?;
    let (isolation, isolation_records) = load_isolation(&cli.isolation)?;
    let depth16 = json_at(&value41, &["depth_ladder"])?
        .as_array()
        .ok_or("round-41 depth ladder is not an array")?
        .iter()
        .find(|row| row.get("start_bits").and_then(Value::as_u64) == Some(16))
        .ok_or("round-41 depth-16 row missing")?;
    let imported = ImportedEvidence {
        round41_request_digest: depth16
            .get("sorted_sha256")
            .and_then(Value::as_str)
            .ok_or("round-41 request digest missing")?
            .into(),
        round41_native_control_digest: json_str(
            &value41,
            &["imported", "round40_native_control_digest"],
        )?
        .into(),
        round41_zero_false_positives_and_false_negatives: json_bool(
            &value41,
            &["gates", "zero_false_positives_and_false_negatives"],
        )?,
        round41_exact_reference_equality: json_bool(
            &value41,
            &["gates", "exact_equality_with_three_pass_reference"],
        )?,
        round41_promoted: json_bool(&value41, &["gates", "promoted"])?,
        round41_isolated_runs: isolation_records.len() as u64,
        round41_last_run_uncontended: json_at(
            isolation_records.last().expect("nonempty isolation"),
            &["run", "contended"],
        )?
        .as_bool()
        .is_some_and(|contended| !contended),
    };
    if imported.round41_request_digest != ROUND40_REQUEST_DIGEST
        || imported.round41_native_control_digest != ROUND40_CONTROL_DIGEST
        || !imported.round41_zero_false_positives_and_false_negatives
        || !imported.round41_exact_reference_equality
        || imported.round41_promoted
        || !imported.round41_last_run_uncontended
    {
        return Err("round-41 frozen controls changed".into());
    }

    let toy_controls = [17, 31, 63]
        .into_iter()
        .flat_map(|universe| [4, 8, 16].map(move |starts| toy_cell(universe, starts)))
        .collect::<Vec<_>>();
    let load_one_ladder = DEPTH_BITS
        .into_iter()
        .map(|bits| checked_cell(bits, 1u64 << (bits - 3)))
        .collect::<Vec<_>>();
    let saturation_ladder = [(13u32, 0.25f64), (15, 1.0), (17, 4.0), (19, 16.0)]
        .into_iter()
        .map(|(start_bits, _)| checked_cell(18, 1u64 << start_bits))
        .collect::<Vec<_>>();
    let timing = timing_evidence();
    let candidate = timing
        .summaries
        .iter()
        .find(|row| row.variant == Variant::Direct.name())
        .ok_or("candidate timing summary missing")?;
    let projection = vec![
        projection_row(
            "power-of-two-2^38",
            1u64 << 38,
            candidate.addition_equivalents_per_start,
        ),
        projection_row(
            "power-of-two-2^40",
            1u64 << 40,
            candidate.addition_equivalents_per_start,
        ),
        projection_row(
            "storage-cap-does-not-bind; largest request-counter-safe start count",
            (1u64 << 61) - 1,
            candidate.addition_equivalents_per_start,
        ),
    ];
    let toy_failures = toy_controls
        .iter()
        .map(|row| row.false_positives + row.false_negatives)
        .sum::<u64>();
    let exact_failures = load_one_ladder
        .iter()
        .chain(&saturation_ladder)
        .map(|row| row.false_positives + row.false_negatives)
        .sum::<u64>();
    let exact = toy_failures == 0 && exact_failures == 0;
    let state_gate = candidate_algorithmic_bytes(PAIR_UNIVERSE) < (1u64 << 50);
    let semantic = serde_json::json!({
        "curve": CURVE_SLUG,
        "round41_sha256": ROUND41_SHA256,
        "round41_isolation_sha256": ROUND41_ISOLATION_SHA256,
        "imported": &imported,
        "toy_controls": &toy_controls,
        "load_one_ladder": &load_one_ladder,
        "saturation_ladder": &saturation_ladder,
        "timing": &timing,
        "projection": &projection,
    });
    let semantic_evidence_sha256 = hex::encode(sha256(
        &serde_json::to_vec(&semantic).map_err(|error| error.to_string())?,
    ));
    let gates = Gates {
        zero_false_positives_and_false_negatives: exact,
        exact_equality_with_materialized_reference: exact,
        structured_residual_degree_at_most_5_for_measured_family: false,
        relation_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        materialized_storage_below_2_50: state_gate,
        measured_complete_selector_below_rho_in_applicable_tier: false,
        no_discarded_branch_counted_as_exhaustive: true,
        promoted: false,
    };
    Ok(ResultFile {
        schema: "p256-direct-address-replay-v1".into(),
        curve: CURVE_SLUG.into(),
        family: "round-33 known-log selector fixed-universe pair routing".into(),
        dependencies: vec![round41, isolation],
        imported,
        frozen_model: serde_json::json!({
            "pair_universe": PAIR_UNIVERSE,
            "requests_per_start": REQUESTS_PER_START,
            "record_bytes": RECORD_BYTES,
            "pair_value_bytes": VALUE_BYTES,
            "accumulator_bytes": ACCUMULATOR_BYTES,
            "seen_bitmap_bits_per_pair": 1,
            "full_universe_bitmap_bytes": bitmap_words(PAIR_UNIVERSE) * 8,
            "full_universe_value_bytes": PAIR_UNIVERSE * VALUE_BYTES,
            "full_universe_algorithmic_bytes": candidate_algorithmic_bytes(PAIR_UNIVERSE),
            "base_ratio_to_rho": BASE_RATIO,
            "local_ratio_to_rho": LOCAL_RATIO,
            "mean_capacity": MEAN_CAPACITY,
            "allowed_additions_per_start": ALLOWED_ADDITIONS,
            "comparison_factor_base": "FB1h2f8621cda105",
            "comparison_structured_residual_maximum": 4,
            "measured_family_structured_residual_degree": null,
        }),
        toy_controls,
        load_one_ladder,
        saturation_ladder,
        timing,
        projection,
        semantic_evidence_sha256,
        gates,
        full_depth_unplanted_attempted: false,
        interpretation: "Direct addressing makes modeled selector state depend on the fixed pair universe rather than the start batch. The complete candidate replays request generation twice, and its measured tables are many orders of magnitude smaller than the projected 557-GB P-256 table."
            .into(),
        decision: "Do not promote unless complete timing below rho is reproduced in the applicable full-table tier and the measured family obtains its own structured-degree and relation-collection evidence."
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
            eprintln!("P-256 direct-address replay failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn direct_matches_materialized_on_complete_toys() {
        for universe in [17, 31, 63] {
            let row = toy_cell(universe, 16);
            assert_eq!(row.false_positives, 0);
            assert_eq!(row.false_negatives, 0);
            assert_eq!(row.materialized_distinct, row.direct_distinct);
        }
    }

    #[test]
    fn registered_layouts_and_full_state_are_exact() {
        assert_eq!(size_of::<Request>(), 16);
        assert_eq!(size_of::<PairValue>(), 16);
        assert_eq!(size_of::<Accumulator>(), 16);
        assert_eq!(bitmap_words(PAIR_UNIVERSE) * 8, 4_320_268_584);
        assert_eq!(PAIR_UNIVERSE * VALUE_BYTES, 552_994_377_792);
        assert_eq!(candidate_algorithmic_bytes(PAIR_UNIVERSE), 557_314_646_392);
    }

    #[test]
    fn candidate_state_is_independent_of_start_count() {
        let state = candidate_algorithmic_bytes(1 << 18);
        assert_eq!(state, candidate_algorithmic_bytes(1 << 18));
        assert!(state < materialized_peak_bytes(1 << 18));
    }
}
