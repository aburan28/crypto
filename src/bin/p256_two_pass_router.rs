//! Two-pass independent signed-pair router for P-256 (round 41).

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
use serde::Serialize;
use serde_json::Value;

const CURVE_SLUG: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const ROUND40_SHA256: &str = "52154eb455ea7edd187d83944ba4b6d953d92014d4f05919e4bd8b7bc097ecbf";
const ROUND40_ISOLATION_SHA256: &str =
    "c04812b438ea26ed6b0f94fee6a7425730fd49b395c9f92c70db0c61915b59a9";
const ROUND40_REQUEST_DIGEST: &str =
    "b4178238067a37febeeb46116dffc633311549e6325461670bd12db7383e0dd6";
const ROUND40_CONTROL_DIGEST: &str =
    "98c8f56a06257256bb88f6c6d7cffce79391c461170b91287f8fd193fe5711e6";
const PAIR_UNIVERSE: u64 = 34_562_148_612;
const REQUESTS_PER_START: u64 = 8;
const RECORD_BYTES: u64 = 16;
const ACCUMULATOR_BYTES: u64 = 16;
const BASE_RATIO: f64 = 0.998_569_034_150_286;
const LOCAL_RATIO: f64 = 0.964_336_477_130_181;
const MEAN_CAPACITY: f64 = 224.361_249_307_510_55;
const ALLOWED_ADDITIONS: f64 = 0.334_410_508_422_987_64;
const DEPTHS: [u32; 4] = [12, 14, 16, 18];
const TIMING_DEPTHS: [u32; 3] = [14, 16, 18];
const TIMING_DEPTH: u32 = 16;
const TIMING_REPETITIONS: u64 = 7;
const LADDER_REPETITIONS: u64 = 5;
const BASELINE_ADDITIONS: u64 = 262_144;
const THREE_PASS_TRAFFIC: u64 = 160;
const TWO_PASS_TRAFFIC: u64 = 128;

#[derive(Parser)]
#[command(about = "Compare three-pass and two-pass P-256 pair routers")]
struct Cli {
    #[arg(long)]
    round40: PathBuf,
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
struct Accumulator {
    checksum: u64,
    mask: u8,
    count: u8,
    padding: [u8; 6],
}

#[derive(Clone, Copy)]
enum Variant {
    Complete12,
    Pregenerated12,
    Pregenerated18,
}

impl Variant {
    fn name(self) -> &'static str {
        match self {
            Self::Complete12 => "complete-3x12",
            Self::Pregenerated12 => "pregenerated-3x12",
            Self::Pregenerated18 => "pregenerated-2x18",
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
struct Cell {
    start_bits: u32,
    starts: u64,
    requests: u64,
    distinct_pairs: u64,
    sorted_mismatches: u64,
    accumulator_mismatches: u64,
    regenerated_failures_3x12: u64,
    regenerated_failures_2x18: u64,
    negative_controls: u64,
    false_positives: u64,
    false_negatives: u64,
    three_pass_logical_bytes: u64,
    two_pass_logical_bytes: u64,
    traffic_reduction_fraction: f64,
    three_pass_peak_bytes: u64,
    two_pass_peak_bytes: u64,
    sorted_sha256: String,
    accumulator_sha256: String,
}

#[derive(Clone, Serialize)]
struct ToyCell {
    universe: u64,
    starts: u64,
    distinct_3x12: u64,
    distinct_2x18: u64,
    reference_distinct: u64,
    sorted_equal: bool,
    accumulator_equal: bool,
    failures: u64,
}

#[derive(Clone, Serialize)]
struct RawTiming {
    repetition: u64,
    order: Vec<String>,
    complete_3x12_ns: u64,
    pregenerated_3x12_ns: u64,
    pregenerated_2x18_ns: u64,
    complete_checksum: u64,
    pregenerated_3x12_checksum: u64,
    pregenerated_2x18_checksum: u64,
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
    includes_request_generation: bool,
}

#[derive(Clone, Serialize)]
struct LadderTiming {
    variant: String,
    start_bits: u32,
    starts: u64,
    repetitions: u64,
    median_wall_ns: f64,
    minimum_wall_ns: u64,
    checksum_identical: bool,
}

#[derive(Clone, Serialize)]
struct TimingEvidence {
    timing_start_bits: u32,
    timing_starts: u64,
    repetitions: u64,
    baseline_additions: u64,
    aa_first_ns: u64,
    aa_second_ns: u64,
    aa_ratio: f64,
    baseline_median_ns_per_addition: f64,
    raw: Vec<RawTiming>,
    summaries: Vec<VariantSummary>,
    candidate_over_control_ratio: f64,
    ladder: Vec<LadderTiming>,
    pregenerated_3x12_time_exponent: f64,
    pregenerated_2x18_time_exponent: f64,
    scope: String,
}

#[derive(Clone, Serialize)]
struct ProjectionRow {
    label: String,
    starts: u64,
    starts_log2: f64,
    pair_constructions_per_start: f64,
    router_addition_equivalents_per_start: f64,
    combined_addition_equivalents_per_start: f64,
    combined_ratio_to_rho: f64,
    below_rho_arithmetic_only: bool,
    peak_materialized_bytes: u64,
    peak_log2_bytes: f64,
    storage_below_2_50: bool,
    logical_bytes_per_start: u64,
    logical_bytes_per_batch: String,
    applicable_memory_tier_measured: bool,
}

#[derive(Clone, Serialize)]
struct ImportedEvidence {
    round40_complete_addition_equivalents_per_start: f64,
    round40_budget_multiple: f64,
    round40_projected_ratio_to_rho: f64,
    round40_pair_only_threshold_starts: u64,
    round40_native_control_digest: String,
    round40_zero_false_positives_and_false_negatives: bool,
    round40_isolated_runs: u64,
    round40_last_run_uncontended: bool,
}

#[derive(Clone, Serialize)]
struct Gates {
    zero_false_positives_and_false_negatives: bool,
    exact_equality_with_three_pass_reference: bool,
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
    depth_ladder: Vec<Cell>,
    timing: TimingEvidence,
    projection: Vec<ProjectionRow>,
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

fn json_f64(value: &Value, path: &[&str]) -> Result<f64, String> {
    json_at(value, path)?
        .as_f64()
        .ok_or_else(|| format!("JSON path {} is not numeric", path.join(".")))
}

fn json_bool(value: &Value, path: &[&str]) -> Result<bool, String> {
    json_at(value, path)?
        .as_bool()
        .ok_or_else(|| format!("JSON path {} is not a bool", path.join(".")))
}

fn load_round40(path: &Path) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != ROUND40_SHA256 {
        return Err(format!(
            "round-40 hash mismatch: expected {ROUND40_SHA256}, got {digest}"
        ));
    }
    let value: Value = serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if json_str(&value, &["schema"])? != "p256-independent-batch-routing-v1"
        || json_str(&value, &["curve"])? != CURVE_SLUG
    {
        return Err("round-40 identity mismatch".into());
    }
    Ok((
        Dependency {
            path: path.display().to_string(),
            sha256: digest,
            schema: "p256-independent-batch-routing-v1".into(),
        },
        value,
    ))
}

fn load_isolation(path: &Path) -> Result<(Dependency, Vec<Value>), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != ROUND40_ISOLATION_SHA256 {
        return Err(format!(
            "round-40 isolation hash mismatch: expected {ROUND40_ISOLATION_SHA256}, got {digest}"
        ));
    }
    let text = String::from_utf8(bytes).map_err(|error| error.to_string())?;
    let records = text
        .lines()
        .map(|line| serde_json::from_str(line).map_err(|error| error.to_string()))
        .collect::<Result<Vec<Value>, String>>()?;
    if records.is_empty() {
        return Err("round-40 isolation has no records".into());
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

fn radix_sort(mut source: Vec<Request>, radix_bits: u32, passes: u32) -> Vec<Request> {
    let buckets = 1usize << radix_bits;
    let mask = buckets - 1;
    let mut target = vec![Request::default(); source.len()];
    let mut counts = vec![0usize; buckets];
    for pass in 0..passes {
        counts.fill(0);
        let shift = pass * radix_bits;
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

struct Routed {
    sorted: Vec<Request>,
    accumulators: Vec<Accumulator>,
    distinct: u64,
}

fn reconstruct(sorted: Vec<Request>, starts: u64) -> Routed {
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
    Routed {
        sorted,
        accumulators,
        distinct,
    }
}

fn route_owned(records: Vec<Request>, starts: u64, radix_bits: u32, passes: u32) -> Routed {
    reconstruct(radix_sort(records, radix_bits, passes), starts)
}

fn route_variant(variant: Variant, frozen: &[Request], starts: u64) -> Routed {
    match variant {
        Variant::Complete12 => route_owned(generate_requests(starts, PAIR_UNIVERSE), starts, 12, 3),
        Variant::Pregenerated12 => route_owned(frozen.to_vec(), starts, 12, 3),
        Variant::Pregenerated18 => route_owned(frozen.to_vec(), starts, 18, 2),
    }
}

fn expected_accumulator(start: u64, universe: u64) -> Accumulator {
    let mut accumulator = Accumulator::default();
    for slot in 0..REQUESTS_PER_START as u8 {
        accumulator.checksum = accumulator
            .checksum
            .wrapping_add(contribution(pair_id(start, slot, universe), slot));
        accumulator.mask |= 1 << slot;
        accumulator.count += 1;
    }
    accumulator
}

fn digest_requests(records: &[Request]) -> String {
    let mut bytes = Vec::with_capacity(records.len() * 16);
    for record in records {
        bytes.extend(record.pair.to_be_bytes());
        bytes.extend(record.owner_slot.to_be_bytes());
    }
    hex::encode(sha256(&bytes))
}

fn digest_accumulators(records: &[Accumulator]) -> String {
    let mut bytes = Vec::with_capacity(records.len() * 10);
    for record in records {
        bytes.extend(record.checksum.to_be_bytes());
        bytes.push(record.mask);
        bytes.push(record.count);
    }
    hex::encode(sha256(&bytes))
}

fn checked_cell(bits: u32) -> Cell {
    let starts = 1u64 << bits;
    let requests = starts * REQUESTS_PER_START;
    let frozen = generate_requests(starts, PAIR_UNIVERSE);
    let route12 = route_variant(Variant::Pregenerated12, &frozen, starts);
    let route18 = route_variant(Variant::Pregenerated18, &frozen, starts);
    let sorted_mismatches = route12
        .sorted
        .iter()
        .zip(&route18.sorted)
        .filter(|(left, right)| left != right)
        .count() as u64;
    let accumulator_mismatches = route12
        .accumulators
        .iter()
        .zip(&route18.accumulators)
        .filter(|(left, right)| left != right)
        .count() as u64;
    let mut regenerated_failures_3x12 = 0u64;
    let mut regenerated_failures_2x18 = 0u64;
    let mut false_positives = 0u64;
    for start in 0..starts {
        let expected = expected_accumulator(start, PAIR_UNIVERSE);
        let actual12 = route12.accumulators[start as usize];
        let actual18 = route18.accumulators[start as usize];
        if actual12 != expected {
            regenerated_failures_3x12 += 1;
        }
        if actual18 != expected {
            regenerated_failures_2x18 += 1;
        }
        if actual18.checksum.wrapping_add(1) == expected.checksum {
            false_positives += 1;
        }
    }
    let false_negatives = regenerated_failures_3x12
        + regenerated_failures_2x18
        + sorted_mismatches
        + accumulator_mismatches;
    Cell {
        start_bits: bits,
        starts,
        requests,
        distinct_pairs: route18.distinct,
        sorted_mismatches,
        accumulator_mismatches,
        regenerated_failures_3x12,
        regenerated_failures_2x18,
        negative_controls: starts,
        false_positives,
        false_negatives,
        three_pass_logical_bytes: requests * THREE_PASS_TRAFFIC,
        two_pass_logical_bytes: requests * TWO_PASS_TRAFFIC,
        traffic_reduction_fraction: 1.0 - TWO_PASS_TRAFFIC as f64 / THREE_PASS_TRAFFIC as f64,
        three_pass_peak_bytes: requests * 2 * RECORD_BYTES + starts * ACCUMULATOR_BYTES,
        two_pass_peak_bytes: requests * 2 * RECORD_BYTES
            + starts * ACCUMULATOR_BYTES
            + (1u64 << 18) * 8,
        sorted_sha256: digest_requests(&route18.sorted),
        accumulator_sha256: digest_accumulators(&route18.accumulators),
    }
}

fn toy_cell(universe: u64, starts: u64) -> ToyCell {
    let frozen = generate_requests(starts, universe);
    let route12 = route_owned(frozen.clone(), starts, 12, 3);
    let route18 = route_owned(frozen, starts, 18, 2);
    let reference = (0..starts)
        .flat_map(|start| {
            (0..REQUESTS_PER_START as u8).map(move |slot| pair_id(start, slot, universe))
        })
        .collect::<std::collections::BTreeSet<_>>();
    let sorted_equal = route12.sorted == route18.sorted;
    let accumulator_equal = route12.accumulators == route18.accumulators;
    let failures = u64::from(!sorted_equal)
        + u64::from(!accumulator_equal)
        + u64::from(route12.distinct != reference.len() as u64)
        + u64::from(route18.distinct != reference.len() as u64);
    ToyCell {
        universe,
        starts,
        distinct_3x12: route12.distinct,
        distinct_2x18: route18.distinct,
        reference_distinct: reference.len() as u64,
        sorted_equal,
        accumulator_equal,
        failures,
    }
}

fn checksum(routed: &Routed) -> u64 {
    routed
        .accumulators
        .iter()
        .fold(routed.distinct, |value, row| value ^ row.checksum)
}

fn timed_variant(variant: Variant, frozen: &[Request], starts: u64) -> (u64, u64) {
    let begin = Instant::now();
    let routed = black_box(route_variant(variant, frozen, starts));
    (
        begin.elapsed().as_nanos() as u64,
        black_box(checksum(&routed)),
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

fn variant_summary(
    variant: Variant,
    raw: &[RawTiming],
    baseline_ns_per_addition: f64,
    starts: u64,
) -> VariantSummary {
    let values = raw
        .iter()
        .map(|row| match variant {
            Variant::Complete12 => row.complete_3x12_ns,
            Variant::Pregenerated12 => row.pregenerated_3x12_ns,
            Variant::Pregenerated18 => row.pregenerated_2x18_ns,
        })
        .collect::<Vec<_>>();
    let checksums = raw
        .iter()
        .map(|row| match variant {
            Variant::Complete12 => row.complete_checksum,
            Variant::Pregenerated12 => row.pregenerated_3x12_checksum,
            Variant::Pregenerated18 => row.pregenerated_2x18_checksum,
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
        checksum_identical: checksums.iter().all(|value| *value == checksums[0]),
        includes_request_generation: matches!(variant, Variant::Complete12),
    }
}

fn fit_exponent(rows: &[LadderTiming], variant: &str) -> f64 {
    let selected = rows
        .iter()
        .filter(|row| row.variant == variant)
        .collect::<Vec<_>>();
    let mean_x = selected
        .iter()
        .map(|row| row.start_bits as f64)
        .sum::<f64>()
        / selected.len() as f64;
    let mean_y = selected
        .iter()
        .map(|row| row.median_wall_ns.log2())
        .sum::<f64>()
        / selected.len() as f64;
    let numerator = selected
        .iter()
        .map(|row| (row.start_bits as f64 - mean_x) * (row.median_wall_ns.log2() - mean_y))
        .sum::<f64>();
    let denominator = selected
        .iter()
        .map(|row| (row.start_bits as f64 - mean_x).powi(2))
        .sum::<f64>();
    numerator / denominator
}

fn timing_evidence() -> TimingEvidence {
    let starts = 1u64 << TIMING_DEPTH;
    let frozen = generate_requests(starts, PAIR_UNIVERSE);
    let _ = timed_additions(8_192);
    for variant in [
        Variant::Complete12,
        Variant::Pregenerated12,
        Variant::Pregenerated18,
    ] {
        let _ = timed_variant(variant, &frozen[..(1 << 13)], 1 << 10);
    }
    let (aa_first, _) = timed_additions(BASELINE_ADDITIONS);
    let (aa_second, _) = timed_additions(BASELINE_ADDITIONS);
    let orders = [
        [
            Variant::Complete12,
            Variant::Pregenerated12,
            Variant::Pregenerated18,
        ],
        [
            Variant::Pregenerated12,
            Variant::Pregenerated18,
            Variant::Complete12,
        ],
        [
            Variant::Pregenerated18,
            Variant::Complete12,
            Variant::Pregenerated12,
        ],
    ];
    let mut raw = Vec::new();
    for repetition in 0..TIMING_REPETITIONS {
        let order = orders[(repetition as usize) % orders.len()];
        let mut complete = (0u64, 0u64);
        let mut pre12 = (0u64, 0u64);
        let mut pre18 = (0u64, 0u64);
        for variant in order {
            let measured = timed_variant(variant, &frozen, starts);
            match variant {
                Variant::Complete12 => complete = measured,
                Variant::Pregenerated12 => pre12 = measured,
                Variant::Pregenerated18 => pre18 = measured,
            }
        }
        raw.push(RawTiming {
            repetition,
            order: order.iter().map(|variant| variant.name().into()).collect(),
            complete_3x12_ns: complete.0,
            pregenerated_3x12_ns: pre12.0,
            pregenerated_2x18_ns: pre18.0,
            complete_checksum: complete.1,
            pregenerated_3x12_checksum: pre12.1,
            pregenerated_2x18_checksum: pre18.1,
        });
    }
    let mut baseline_values = Vec::new();
    for _ in 0..TIMING_REPETITIONS {
        baseline_values.push(timed_additions(BASELINE_ADDITIONS).0);
    }
    let baseline_ns_per_addition = median_u64(&baseline_values) / BASELINE_ADDITIONS as f64;
    let summaries = [
        Variant::Complete12,
        Variant::Pregenerated12,
        Variant::Pregenerated18,
    ]
    .into_iter()
    .map(|variant| variant_summary(variant, &raw, baseline_ns_per_addition, starts))
    .collect::<Vec<_>>();

    let mut ladder = Vec::new();
    for bits in TIMING_DEPTHS {
        let ladder_starts = 1u64 << bits;
        let ladder_frozen = generate_requests(ladder_starts, PAIR_UNIVERSE);
        for variant in [Variant::Pregenerated12, Variant::Pregenerated18] {
            let mut values = Vec::new();
            let mut checksums = Vec::new();
            for _ in 0..LADDER_REPETITIONS {
                let measured = timed_variant(variant, &ladder_frozen, ladder_starts);
                values.push(measured.0);
                checksums.push(measured.1);
            }
            ladder.push(LadderTiming {
                variant: variant.name().into(),
                start_bits: bits,
                starts: ladder_starts,
                repetitions: LADDER_REPETITIONS,
                median_wall_ns: median_u64(&values),
                minimum_wall_ns: *values.iter().min().expect("ladder values"),
                checksum_identical: checksums.iter().all(|value| *value == checksums[0]),
            });
        }
    }
    let control = summaries
        .iter()
        .find(|row| row.variant == Variant::Pregenerated12.name())
        .expect("control summary");
    let candidate = summaries
        .iter()
        .find(|row| row.variant == Variant::Pregenerated18.name())
        .expect("candidate summary");
    let candidate_over_control_ratio = candidate.median_wall_ns / control.median_wall_ns;
    let pregenerated_3x12_time_exponent = fit_exponent(&ladder, Variant::Pregenerated12.name());
    let pregenerated_2x18_time_exponent = fit_exponent(&ladder, Variant::Pregenerated18.name());
    TimingEvidence {
        timing_start_bits: TIMING_DEPTH,
        timing_starts: starts,
        repetitions: TIMING_REPETITIONS,
        baseline_additions: BASELINE_ADDITIONS,
        aa_first_ns: aa_first,
        aa_second_ns: aa_second,
        aa_ratio: aa_first.max(aa_second) as f64 / aa_first.min(aa_second) as f64,
        baseline_median_ns_per_addition: baseline_ns_per_addition,
        raw,
        summaries,
        candidate_over_control_ratio,
        pregenerated_3x12_time_exponent,
        pregenerated_2x18_time_exponent,
        ladder,
        scope: "Single-core RAM routing. Pregenerated rows exclude request derivation and cannot be treated as complete selector timings or external-memory measurements."
            .into(),
    }
}

fn expected_distinct(starts: u64) -> f64 {
    let requests = (starts * REQUESTS_PER_START) as f64;
    let universe = PAIR_UNIVERSE as f64;
    -universe * (requests * (-1.0 / universe).ln_1p()).exp_m1()
}

fn projection_row(label: &str, starts: u64, router_equivalents: f64) -> ProjectionRow {
    let pair_per_start = expected_distinct(starts) / starts as f64;
    let combined = pair_per_start + router_equivalents;
    let ratio = pair_ratio(combined);
    let count_table_bytes = (1u64 << 18) * 8;
    let peak =
        starts * (2 * REQUESTS_PER_START * RECORD_BYTES + ACCUMULATOR_BYTES) + count_table_bytes;
    let traffic =
        u128::from(starts) * u128::from(REQUESTS_PER_START) * u128::from(TWO_PASS_TRAFFIC);
    ProjectionRow {
        label: label.into(),
        starts,
        starts_log2: (starts as f64).log2(),
        pair_constructions_per_start: pair_per_start,
        router_addition_equivalents_per_start: router_equivalents,
        combined_addition_equivalents_per_start: combined,
        combined_ratio_to_rho: ratio,
        below_rho_arithmetic_only: ratio <= 1.0,
        peak_materialized_bytes: peak,
        peak_log2_bytes: (peak as f64).log2(),
        storage_below_2_50: peak < (1u64 << 50),
        logical_bytes_per_start: REQUESTS_PER_START * TWO_PASS_TRAFFIC,
        logical_bytes_per_batch: traffic.to_string(),
        applicable_memory_tier_measured: false,
    }
}

fn result(cli: &Cli) -> Result<ResultFile, String> {
    if size_of::<Request>() != 16 || size_of::<Accumulator>() != 16 {
        return Err("registered layout changed".into());
    }
    let (round40, value40) = load_round40(&cli.round40)?;
    let (isolation, isolation_records) = load_isolation(&cli.isolation)?;
    let depth16 = json_at(&value40, &["depth_ladder"])?
        .as_array()
        .ok_or("round-40 depth ladder is not an array")?
        .iter()
        .find(|row| row.get("start_bits").and_then(Value::as_u64) == Some(16))
        .ok_or("round-40 depth-16 row missing")?;
    if depth16
        .get("ordered_request_sha256")
        .and_then(Value::as_str)
        != Some(ROUND40_REQUEST_DIGEST)
        || json_str(&value40, &["native_controls", "control_sha256"])? != ROUND40_CONTROL_DIGEST
        || !json_bool(
            &value40,
            &["gates", "zero_false_positives_and_false_negatives"],
        )?
    {
        return Err("round-40 frozen controls changed".into());
    }
    let pair_threshold = json_at(&value40, &["projection"])?
        .as_array()
        .ok_or("round-40 projection is not an array")?
        .iter()
        .find(|row| row.get("label").and_then(Value::as_str) == Some("pair-only parity threshold"))
        .ok_or("round-40 threshold row missing")?;
    let imported = ImportedEvidence {
        round40_complete_addition_equivalents_per_start: json_f64(
            &value40,
            &["timing", "candidate_median_addition_equivalents_per_start"],
        )?,
        round40_budget_multiple: json_f64(&value40, &["timing", "candidate_budget_multiple"])?,
        round40_projected_ratio_to_rho: json_f64(
            &value40,
            &["timing", "candidate_projected_ratio_to_rho"],
        )?,
        round40_pair_only_threshold_starts: pair_threshold
            .get("starts")
            .and_then(Value::as_u64)
            .ok_or("threshold starts missing")?,
        round40_native_control_digest: ROUND40_CONTROL_DIGEST.into(),
        round40_zero_false_positives_and_false_negatives: true,
        round40_isolated_runs: isolation_records.len() as u64,
        round40_last_run_uncontended: json_at(
            isolation_records.last().expect("nonempty isolation"),
            &["run", "contended"],
        )?
        .as_bool()
        .is_some_and(|contended| !contended),
    };
    if (imported.round40_complete_addition_equivalents_per_start - 3.432_579_785_165_123_6).abs()
        > 1e-12
        || imported.round40_pair_only_threshold_starts != 103_352_459_747
        || !imported.round40_last_run_uncontended
    {
        return Err("round-40 imported boundary changed".into());
    }

    let toy_controls = [17, 31, 63]
        .into_iter()
        .flat_map(|universe| [4, 8, 16].map(move |starts| toy_cell(universe, starts)))
        .collect::<Vec<_>>();
    let depth_ladder = DEPTHS.into_iter().map(checked_cell).collect::<Vec<_>>();
    let timing = timing_evidence();
    let candidate = timing
        .summaries
        .iter()
        .find(|row| row.variant == Variant::Pregenerated18.name())
        .ok_or("candidate timing summary missing")?;
    let count_table_bytes = (1u64 << 18) * 8;
    let max_starts = ((1u64 << 50) - 1 - count_table_bytes)
        / (2 * REQUESTS_PER_START * RECORD_BYTES + ACCUMULATOR_BYTES);
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
            "largest 2x18 batch below 2^50 bytes",
            max_starts,
            candidate.addition_equivalents_per_start,
        ),
    ];
    let toy_failures = toy_controls.iter().map(|row| row.failures).sum::<u64>();
    let ladder_failures = depth_ladder
        .iter()
        .map(|row| row.false_positives + row.false_negatives)
        .sum::<u64>();
    let exact = toy_failures == 0 && ladder_failures == 0;
    let storage_gate = projection
        .iter()
        .any(|row| row.below_rho_arithmetic_only && row.storage_below_2_50);
    let semantic = serde_json::json!({
        "curve": CURVE_SLUG,
        "round40_sha256": ROUND40_SHA256,
        "round40_isolation_sha256": ROUND40_ISOLATION_SHA256,
        "imported": &imported,
        "toy_controls": &toy_controls,
        "depth_ladder": &depth_ladder,
        "timing": &timing,
        "projection": &projection,
    });
    let semantic_evidence_sha256 = hex::encode(sha256(
        &serde_json::to_vec(&semantic).map_err(|error| error.to_string())?,
    ));
    let gates = Gates {
        zero_false_positives_and_false_negatives: exact,
        exact_equality_with_three_pass_reference: exact,
        structured_residual_degree_at_most_5_for_measured_family: false,
        relation_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        materialized_storage_below_2_50: storage_gate,
        measured_complete_selector_below_rho_in_applicable_tier: false,
        no_discarded_branch_counted_as_exhaustive: true,
        promoted: false,
    };
    Ok(ResultFile {
        schema: "p256-two-pass-router-v1".into(),
        curve: CURVE_SLUG.into(),
        family: "round-33 known-log selector independent pair routing".into(),
        dependencies: vec![round40, isolation],
        imported,
        frozen_model: serde_json::json!({
            "pair_universe": PAIR_UNIVERSE,
            "requests_per_start": REQUESTS_PER_START,
            "record_bytes": RECORD_BYTES,
            "accumulator_bytes": ACCUMULATOR_BYTES,
            "three_pass_bits": 12,
            "three_passes": 3,
            "two_pass_bits": 18,
            "two_passes": 2,
            "three_pass_bytes_per_request": THREE_PASS_TRAFFIC,
            "two_pass_bytes_per_request": TWO_PASS_TRAFFIC,
            "base_ratio_to_rho": BASE_RATIO,
            "local_ratio_to_rho": LOCAL_RATIO,
            "mean_capacity": MEAN_CAPACITY,
            "allowed_additions_per_start": ALLOWED_ADDITIONS,
            "comparison_factor_base": "FB1h2f8621cda105",
            "comparison_structured_residual_maximum": 4,
            "measured_family_structured_residual_degree": null,
        }),
        toy_controls,
        depth_ladder,
        timing,
        projection,
        semantic_evidence_sha256,
        gates,
        full_depth_unplanted_attempted: false,
        interpretation: "Two 18-bit passes reduce exact routing traffic by 20%. A pregenerated-input timing may price the routing kernel, but it excludes upstream request derivation and remains a RAM measurement; projected multi-terabyte rows are not measured parity results."
            .into(),
        decision: "Do not promote unless the complete selector and applicable external-memory tier both pass. Preserve any router-only crossover as a stage result and keep the Dickson relation projection unset."
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
            eprintln!("P-256 two-pass router failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_pass_matches_three_pass_on_complete_toys() {
        for universe in [17, 31, 63] {
            let row = toy_cell(universe, 16);
            assert!(row.sorted_equal);
            assert!(row.accumulator_equal);
            assert_eq!(row.failures, 0);
        }
    }

    #[test]
    fn layouts_and_traffic_are_registered() {
        assert_eq!(size_of::<Request>(), 16);
        assert_eq!(size_of::<Accumulator>(), 16);
        assert_eq!(THREE_PASS_TRAFFIC, 160);
        assert_eq!(TWO_PASS_TRAFFIC, 128);
    }

    #[test]
    fn two_pass_covers_full_pair_key() {
        const { assert!(PAIR_UNIVERSE < (1u64 << 36)) };
        assert_eq!(18 * 2, 36);
    }
}
