//! Execute 289 exact P-256 factor-base screening rounds (rounds 2 through 290).

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use blake3::Hasher;
use clap::Parser;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;
use serde_json::Value;

const CURVE: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const COMPARISON_FB: &str = "FB1h2f8621cda105";
const ROUND1_FB_SHA256: &str = "bf4c0f95ca46bf92b236eea7b8f37a3902f744b07a898a41f6b120faab2e35a4";
const ROUND1_DEGREE_SHA256: &str =
    "efe17e77cd23affe03b7b2ee8146490409afbd779f694058606caf798909c3be";
const ROUND292_SHA256: &str = "19d546fd2e53728f029d349080edbf7337ac2ada785446978b91a62eabbe3821";
const ROUND31_R6935_COEFFICIENT_SHA256: &str =
    "980917981827d813e60484abb0655e8bd527b0beb540ea146d2974ff17303a53";
const COLUMNS: u64 = 131_458;
const ARITY: usize = 17;
const CUTOFF: u64 = 219;
const SIGN_BITS: u32 = 17;
const SIGNS: u64 = 1 << SIGN_BITS;
const LOCAL_RATIO: f64 = 0.964_336_477_130_181;
const ROUND292_GATE_FRACTION: f64 = 0.968_617_146_446_504_4;
const ROUND_FIRST: u32 = 2;
const ROUND_LAST: u32 = 290;
const RARE_FIRST: u64 = 6_647;
const RARE_LAST: u64 = 7_223;
const ROUND_COUNT: usize = (ROUND_LAST - ROUND_FIRST + 1) as usize;

#[derive(Parser)]
#[command(about = "Execute exact P-256 factor-base screens for rounds 2 through 290")]
struct Cli {
    #[arg(long)]
    round1_factor_base: PathBuf,
    #[arg(long)]
    round1_degree: PathBuf,
    #[arg(long)]
    round292: PathBuf,
    #[arg(long)]
    round_dir: PathBuf,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Clone, Serialize)]
struct Dependency {
    path: String,
    sha256: String,
    schema: String,
}

#[derive(Clone, Eq, PartialEq)]
struct Interval {
    start: BigUint,
    end: BigUint,
}

#[derive(Clone)]
struct PreparedTuple {
    attempt: u64,
    columns: Vec<u64>,
    slot_capacities: Vec<u64>,
    capacity: u64,
    original_coefficients: Vec<BigUint>,
    endpoint_coefficients: Vec<BigUint>,
}

struct FactorBaseCandidate {
    candidate_id: String,
    common_edges: u64,
    rare_delta: BigUint,
    anchor_attempts: u64,
    anchor: BigUint,
    coefficients: Vec<BigUint>,
    coefficient_sha256: String,
    distinct_relative_coefficients: u64,
    unique_up_to_sign: bool,
    nonidentity: bool,
    edge_replay_failures: u64,
    construction_wall_ns: u64,
}

struct OrbitSupport {
    merged: Vec<Interval>,
    states_emitted: u64,
    exact_distinct_states: u64,
    source_intervals: u64,
    merged_intervals: u64,
    wrap_splits: u64,
    sort_comparisons: u64,
    merge_comparisons: u64,
    generation_wall_ns: u64,
    union_wall_ns: u64,
    support_digest: String,
    start_digest: String,
}

#[derive(Clone, Serialize)]
struct ToyControl {
    prime: u64,
    arity: u64,
    sign_bits: u32,
    states_emitted: u64,
    direct_distinct_states: u64,
    interval_distinct_states: u64,
    false_positives: u64,
    false_negatives: u64,
    exact: bool,
}

#[derive(Clone, Serialize)]
struct NativeControl {
    rare_edges: u64,
    sign_bits: u32,
    states_emitted: u64,
    direct_distinct_states: u64,
    interval_distinct_states: u64,
    false_positives: u64,
    false_negatives: u64,
    exact: bool,
}

#[derive(Clone, Serialize)]
struct RoundTiming {
    factor_base_wall_ns: u64,
    tuple_wall_ns: u64,
    orbit_generation_wall_ns: u64,
    interval_union_wall_ns: u64,
    total_wall_ns: u64,
    process_cpu_ns: u64,
    peak_rss_bytes: Option<u64>,
}

#[derive(Clone, Serialize)]
struct OperationCounts {
    coefficients_constructed: u64,
    factor_base_edges_replayed: u64,
    signed_start_initial_terms: u64,
    gray_boundary_updates: u64,
    path_additions: u64,
    sign_boundary_additions: u64,
    initial_sum_additions: u64,
    boundary_precomputation_additions: u64,
    target_correction_additions: u64,
    charged_addition_equivalents: u64,
    sort_comparisons: u64,
    merge_comparisons: u64,
}

#[derive(Clone, Serialize)]
struct FactorBaseReceipt {
    id: String,
    columns: u64,
    rare_edges: u64,
    common_edges: u64,
    delta_common: String,
    delta_rare: String,
    anchor_attempts: u64,
    anchor_coefficient: String,
    coefficient_sha256: String,
    coefficient_bytes: u64,
    primitive_cycle: bool,
    distinct_relative_coefficients: u64,
    unique_up_to_sign: bool,
    nonidentity: bool,
    edge_replay_failures: u64,
}

#[derive(Clone, Serialize)]
struct TupleReceipt {
    attempt: u64,
    columns: Vec<u64>,
    slot_capacities: Vec<u64>,
    total_capacity: u64,
    tuple_sha256: String,
}

#[derive(Clone, Serialize)]
struct CoverageReceipt {
    sign_bits: u32,
    sign_segments: u64,
    states_emitted: u64,
    per_segment_distinct_sum: u64,
    exact_union_distinct_states: u64,
    within_segment_duplicates: u64,
    cross_segment_duplicates: u64,
    total_duplicates: u64,
    distinct_fraction: f64,
    source_intervals: u64,
    merged_intervals: u64,
    wrap_splits: u64,
    support_digest: String,
    start_digest: String,
}

#[derive(Clone, Serialize)]
struct BoundaryReceipt {
    rare_fraction: f64,
    hidden_collision_probability: f64,
    hidden_ratio_to_rho: f64,
    support_upper_log2: f64,
    support_upper_fraction_of_group: f64,
    target_retry_lower_bound: f64,
    factor_base_complete_ratio_lower_bound: f64,
    ideal_local_ratio_per_distinct: f64,
    factor_base_corrected_ratio_per_distinct: f64,
}

#[derive(Clone, Serialize)]
struct RoundReceipt {
    schema: String,
    curve: String,
    comparison_factor_base: String,
    screening_round: u32,
    execution_status: String,
    factor_base: FactorBaseReceipt,
    tuple: TupleReceipt,
    coverage: CoverageReceipt,
    boundary: BoundaryReceipt,
    operations: OperationCounts,
    timing: RoundTiming,
    logical_materialized_bytes: u64,
    exact_factor_base: bool,
    exact_coverage: bool,
    selector_gate_passed: bool,
    structured_residual_degree: Option<u32>,
    relations_reported: u64,
    full_depth_unplanted_relation_attempted: bool,
    interpretation: String,
}

#[derive(Clone, Serialize)]
struct ManifestEntry {
    round: u32,
    rare_edges: u64,
    factor_base_id: String,
    path: String,
    bytes: u64,
    sha256: String,
    execution_status: String,
    exact_distinct_states: u64,
    distinct_fraction: f64,
    ideal_local_ratio_per_distinct: f64,
    factor_base_corrected_ratio_per_distinct: f64,
    selector_gate_passed: bool,
}

#[derive(Clone, Serialize)]
struct ExecutionGates {
    exactly_289_rounds_2_through_290_once: bool,
    every_round_complete: bool,
    every_factor_base_exact_and_primitive: bool,
    zero_false_positives_and_false_negatives: bool,
    manifest_hashes_and_byte_counts_replayed: bool,
    no_projection_counted_as_execution: bool,
    peak_materialized_storage_below_2_50: bool,
    execution_gate_passed: bool,
    any_selector_gate_passed: bool,
    same_family_structured_degree_at_most_5: bool,
    per_usable_relation_below_2_103: bool,
    relation_collection_below_2_120: bool,
    complete_non_generic_dlp_below_rho: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct Manifest {
    schema: String,
    curve: String,
    family: String,
    comparison_factor_base: String,
    dependencies: Vec<Dependency>,
    round_range: [u32; 2],
    registered_rounds: u64,
    executed_rounds: u64,
    projection_only_rounds: u64,
    frozen: Value,
    toy_controls: Vec<ToyControl>,
    native_controls: Vec<NativeControl>,
    entries: Vec<ManifestEntry>,
    best_ideal_local_round: u32,
    best_factor_base_corrected_round: u32,
    total_wall_ns: u64,
    total_process_cpu_ns: u64,
    maximum_peak_rss_bytes: Option<u64>,
    total_artifact_bytes: u64,
    gates: ExecutionGates,
    semantic_evidence_sha256: String,
    relations_reported: u64,
    full_depth_unplanted_relation_attempted: bool,
    decision: String,
}

fn json_str<'a>(value: &'a Value, key: &str) -> Result<&'a str, String> {
    value
        .get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("JSON key {key} is missing or not a string"))
}

fn load_dependency(
    path: &Path,
    expected_sha256: &str,
    expected_schema: &str,
) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != expected_sha256 {
        return Err(format!(
            "{} hash mismatch: expected {expected_sha256}, got {digest}",
            path.display()
        ));
    }
    let value: Value = serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if json_str(&value, "schema")? != expected_schema || json_str(&value, "curve")? != CURVE {
        return Err(format!("{} identity mismatch", path.display()));
    }
    Ok((
        Dependency {
            path: path.display().to_string(),
            sha256: digest,
            schema: expected_schema.into(),
        },
        value,
    ))
}

fn fixed_32(value: &BigUint) -> Result<[u8; 32], String> {
    let bytes = value.to_bytes_be();
    if bytes.len() > 32 {
        return Err("value exceeds 32 bytes".into());
    }
    let mut fixed = [0u8; 32];
    fixed[32 - bytes.len()..].copy_from_slice(&bytes);
    Ok(fixed)
}

fn log2_big(value: &BigUint) -> Result<f64, String> {
    if value.is_zero() {
        return Err("log2(0)".into());
    }
    let bits = value.bits();
    let kept = bits.min(53);
    let shift = bits - kept;
    let top = (value >> shift)
        .to_u64()
        .ok_or("top bits did not fit u64")?;
    Ok(shift as f64 + (top as f64).log2())
}

fn process_cpu_ns() -> u64 {
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    if unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage) } != 0 {
        return 0;
    }
    let user = usage.ru_utime.tv_sec as u64 * 1_000_000_000 + usage.ru_utime.tv_usec as u64 * 1_000;
    let system =
        usage.ru_stime.tv_sec as u64 * 1_000_000_000 + usage.ru_stime.tv_usec as u64 * 1_000;
    user + system
}

fn peak_rss_bytes() -> Option<u64> {
    let status = fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|line| line.starts_with("VmHWM:"))?;
    let kib = line.split_whitespace().nth(1)?.parse::<u64>().ok()?;
    Some(kib * 1024)
}

fn screening_rare(round: u32) -> Result<u64, String> {
    if !(ROUND_FIRST..=ROUND_LAST).contains(&round) {
        return Err(format!("screening round {round} is out of range"));
    }
    Ok(RARE_FIRST + 2 * u64::from(round - ROUND_FIRST))
}

fn mechanical_rare(index: u64, rare_edges: u64) -> bool {
    (index + 1) * rare_edges / COLUMNS > index * rare_edges / COLUMNS
}

fn derive_anchor(rare_edges: u64, attempt: u64, modulus: &BigUint) -> BigUint {
    let digest =
        sha256(format!("{CURVE}/low-delta-round31/anchor/{rare_edges}/{attempt}").as_bytes());
    BigUint::from_bytes_be(&digest) % modulus
}

fn construct_factor_base(
    rare_edges: u64,
    modulus: &BigUint,
) -> Result<FactorBaseCandidate, String> {
    let started = Instant::now();
    if rare_edges.gcd(&COLUMNS) != 1 {
        return Err(format!("rare edge count {rare_edges} is not primitive"));
    }
    let common_edges = COLUMNS - rare_edges;
    let inverse = BigUint::from(rare_edges).modpow(&(modulus - 2u8), modulus);
    let rare_delta = (modulus - (BigUint::from(common_edges) * inverse % modulus)) % modulus;
    let mut relative = Vec::with_capacity(COLUMNS as usize);
    let mut current = BigUint::zero();
    let mut counted_rare = 0u64;
    for index in 0..COLUMNS {
        relative.push(current.clone());
        if mechanical_rare(index, rare_edges) {
            current = (current + &rare_delta) % modulus;
            counted_rare += 1;
        } else {
            current += 1u8;
            if current >= *modulus {
                current -= modulus;
            }
        }
    }
    if !current.is_zero() || counted_rare != rare_edges {
        return Err(format!(
            "rare edge count {rare_edges} did not close its cycle"
        ));
    }
    let distinct_relative_coefficients = relative.iter().collect::<BTreeSet<_>>().len() as u64;
    if distinct_relative_coefficients != COLUMNS {
        return Err(format!(
            "rare edge count {rare_edges} is not a single cycle"
        ));
    }

    let mut anchor_attempts = 0u64;
    let mut accepted = None;
    for attempt in 0..4_096u64 {
        anchor_attempts += 1;
        let anchor = derive_anchor(rare_edges, attempt, modulus);
        let coefficients = relative
            .iter()
            .map(|value| (value + &anchor) % modulus)
            .collect::<Vec<_>>();
        if coefficients.iter().any(BigUint::is_zero) {
            continue;
        }
        let signed = coefficients
            .iter()
            .map(|value| value.clone().min(modulus - value))
            .collect::<BTreeSet<_>>();
        if signed.len() == COLUMNS as usize {
            accepted = Some((anchor, coefficients));
            break;
        }
    }
    let (anchor, coefficients) =
        accepted.ok_or_else(|| format!("rare edge count {rare_edges} has no valid anchor"))?;
    let mut edge_replay_failures = 0u64;
    for index in 0..COLUMNS as usize {
        let delta = if mechanical_rare(index as u64, rare_edges) {
            &rare_delta
        } else {
            &BigUint::one()
        };
        if (&coefficients[index] + delta) % modulus != coefficients[(index + 1) % COLUMNS as usize]
        {
            edge_replay_failures += 1;
        }
    }
    let mut coefficient_bytes = Vec::with_capacity(COLUMNS as usize * 32);
    for coefficient in &coefficients {
        coefficient_bytes.extend(fixed_32(coefficient)?);
    }
    let coefficient_sha256 = hex::encode(sha256(&coefficient_bytes));
    let candidate_id = format!("LD2R{rare_edges}h{}", &coefficient_sha256[..12]);
    Ok(FactorBaseCandidate {
        candidate_id,
        common_edges,
        rare_delta,
        anchor_attempts,
        anchor,
        coefficients,
        coefficient_sha256,
        distinct_relative_coefficients,
        unique_up_to_sign: true,
        nonidentity: true,
        edge_replay_failures,
        construction_wall_ns: started.elapsed().as_nanos() as u64,
    })
}

fn residual_capacities(rare_edges: u64) -> Vec<u64> {
    let mut capacities = vec![0u64; COLUMNS as usize];
    let mut distance = 0u64;
    for doubled in (0..2 * COLUMNS).rev() {
        let index = doubled % COLUMNS;
        if mechanical_rare(index, rare_edges) {
            distance = 0;
        } else {
            distance += 1;
        }
        if doubled < COLUMNS {
            capacities[index as usize] = distance;
        }
    }
    capacities
}

fn draw_columns(attempt: u64) -> Vec<u64> {
    let mut columns = BTreeSet::new();
    let mut block = 0u64;
    while columns.len() < ARITY {
        let digest =
            sha256(format!("{CURVE}/biased-restart-round33/columns/{attempt}/{block}").as_bytes());
        for chunk in digest.chunks_exact(4) {
            let value = u32::from_be_bytes(chunk.try_into().expect("four bytes"));
            columns.insert(u64::from(value) % COLUMNS);
            if columns.len() == ARITY {
                break;
            }
        }
        block += 1;
    }
    columns.into_iter().collect()
}

fn add_mod(sum: &mut BigUint, delta: &BigUint, positive: bool, modulus: &BigUint) {
    if positive {
        *sum += delta;
        if *sum >= *modulus {
            *sum -= modulus;
        }
    } else if *sum >= *delta {
        *sum -= delta;
    } else {
        *sum += modulus;
        *sum -= delta;
    }
}

fn select_tuple(
    rare_edges: u64,
    coefficients: &[BigUint],
    modulus: &BigUint,
) -> Result<PreparedTuple, String> {
    let capacities = residual_capacities(rare_edges);
    for attempt in 0..1_000_000u64 {
        let columns = draw_columns(attempt);
        let slot_capacities = columns
            .iter()
            .map(|column| capacities[*column as usize])
            .collect::<Vec<_>>();
        let capacity = slot_capacities.iter().sum::<u64>();
        if capacity < CUTOFF {
            continue;
        }
        let mut original_coefficients = Vec::with_capacity(ARITY);
        let mut endpoint_coefficients = Vec::with_capacity(ARITY);
        for (column, slot_capacity) in columns.iter().zip(&slot_capacities) {
            let original = coefficients[*column as usize].clone();
            let endpoint_column = (column + slot_capacity) % COLUMNS;
            let endpoint = coefficients[endpoint_column as usize].clone();
            let mut expected = original.clone();
            add_mod(&mut expected, &BigUint::from(*slot_capacity), true, modulus);
            if endpoint != expected {
                return Err(format!(
                    "tuple endpoint mismatch for rare count {rare_edges}"
                ));
            }
            original_coefficients.push(original);
            endpoint_coefficients.push(endpoint);
        }
        return Ok(PreparedTuple {
            attempt,
            columns,
            slot_capacities,
            capacity,
            original_coefficients,
            endpoint_coefficients,
        });
    }
    Err(format!(
        "rare edge count {rare_edges} found no tuple at cutoff {CUTOFF}"
    ))
}

fn offset_mod(start: &BigUint, offset: i64, modulus: &BigUint) -> BigUint {
    if offset >= 0 {
        (start + BigUint::from(offset as u64)) % modulus
    } else {
        let amount = BigUint::from(offset.unsigned_abs());
        if start >= &amount {
            start - amount
        } else {
            modulus - (amount - start)
        }
    }
}

fn split_interval(
    start: &BigUint,
    minimum: i64,
    maximum: i64,
    modulus: &BigUint,
) -> (Vec<Interval>, bool) {
    let low = offset_mod(start, minimum, modulus);
    let width = (maximum - minimum) as u64;
    let high_unwrapped = &low + BigUint::from(width);
    if &high_unwrapped < modulus {
        (
            vec![Interval {
                start: low,
                end: high_unwrapped,
            }],
            false,
        )
    } else {
        (
            vec![
                Interval {
                    start: BigUint::zero(),
                    end: &high_unwrapped - modulus,
                },
                Interval {
                    start: low,
                    end: modulus - BigUint::one(),
                },
            ],
            true,
        )
    }
}

fn interval_width(interval: &Interval) -> Result<u64, String> {
    (&interval.end - &interval.start + BigUint::one())
        .to_u64()
        .ok_or_else(|| "interval width exceeds u64".into())
}

fn union_intervals(mut intervals: Vec<Interval>) -> Result<OrbitSupport, String> {
    let started = Instant::now();
    let source_intervals = intervals.len() as u64;
    let mut sort_comparisons = 0u64;
    intervals.sort_unstable_by(|left, right| {
        sort_comparisons += 1;
        left.start
            .cmp(&right.start)
            .then_with(|| left.end.cmp(&right.end))
    });
    let mut merged: Vec<Interval> = Vec::new();
    let mut merge_comparisons = 0u64;
    for interval in intervals {
        if let Some(previous) = merged.last_mut() {
            merge_comparisons += 1;
            if interval.start <= &previous.end + BigUint::one() {
                if interval.end > previous.end {
                    previous.end = interval.end;
                }
                continue;
            }
        }
        merged.push(interval);
    }
    let exact_distinct_states = merged
        .iter()
        .map(interval_width)
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .sum();
    let mut digest = Hasher::new();
    for interval in &merged {
        digest.update(&fixed_32(&interval.start)?);
        digest.update(&fixed_32(&interval.end)?);
    }
    Ok(OrbitSupport {
        merged_intervals: merged.len() as u64,
        merged,
        states_emitted: 0,
        exact_distinct_states,
        source_intervals,
        wrap_splits: 0,
        sort_comparisons,
        merge_comparisons,
        generation_wall_ns: 0,
        union_wall_ns: started.elapsed().as_nanos() as u64,
        support_digest: digest.finalize().to_hex().to_string(),
        start_digest: String::new(),
    })
}

fn initial_corners(tuple: &PreparedTuple, modulus: &BigUint) -> (BigUint, BigUint) {
    let mut forward = BigUint::zero();
    let mut reverse = BigUint::zero();
    for slot in 0..tuple.original_coefficients.len() {
        add_mod(
            &mut forward,
            &tuple.original_coefficients[slot],
            true,
            modulus,
        );
        add_mod(
            &mut reverse,
            &tuple.endpoint_coefficients[slot],
            true,
            modulus,
        );
    }
    (forward, reverse)
}

fn update_corners(
    forward: &mut BigUint,
    reverse: &mut BigUint,
    tuple: &PreparedTuple,
    changed: usize,
    old_positive: bool,
    modulus: &BigUint,
) {
    let boundary =
        (&tuple.original_coefficients[changed] + &tuple.endpoint_coefficients[changed]) % modulus;
    add_mod(forward, &boundary, !old_positive, modulus);
    add_mod(reverse, &boundary, !old_positive, modulus);
}

fn monotone_orbit(
    tuple: &PreparedTuple,
    sign_bits: u32,
    modulus: &BigUint,
) -> Result<OrbitSupport, String> {
    let started = Instant::now();
    let signs = 1u64 << sign_bits;
    let mut intervals = Vec::with_capacity(signs as usize + 2);
    let mut wrap_splits = 0u64;
    let mut start_digest = Hasher::new();
    let (mut forward_corner, mut reverse_corner) = initial_corners(tuple, modulus);
    for ordinal in 0..signs {
        let gray = ordinal ^ (ordinal >> 1);
        let forward = ordinal & 1 == 0;
        let start = if forward {
            &forward_corner
        } else {
            &reverse_corner
        };
        start_digest.update(&ordinal.to_be_bytes());
        start_digest.update(&gray.to_be_bytes());
        start_digest.update(&fixed_32(start)?);
        let (minimum, maximum) = if forward {
            (0, tuple.capacity as i64)
        } else {
            (-(tuple.capacity as i64), 0)
        };
        let (mut pieces, wrapped) = split_interval(start, minimum, maximum, modulus);
        wrap_splits += u64::from(wrapped);
        intervals.append(&mut pieces);
        if ordinal + 1 != signs {
            let next_gray = (ordinal + 1) ^ ((ordinal + 1) >> 1);
            let changed = (gray ^ next_gray).trailing_zeros() as usize;
            let old_positive = gray & (1 << changed) == 0;
            update_corners(
                &mut forward_corner,
                &mut reverse_corner,
                tuple,
                changed,
                old_positive,
                modulus,
            );
        }
    }
    let generation_wall_ns = started.elapsed().as_nanos() as u64;
    let mut support = union_intervals(intervals)?;
    support.states_emitted = signs * (tuple.capacity + 1);
    support.wrap_splits = wrap_splits;
    support.generation_wall_ns = generation_wall_ns;
    support.start_digest = start_digest.finalize().to_hex().to_string();
    Ok(support)
}

fn interval_contains(intervals: &[Interval], value: &BigUint) -> bool {
    let index = intervals.partition_point(|interval| interval.end < *value);
    index < intervals.len() && intervals[index].start <= *value
}

fn direct_monotone_support(
    tuple: &PreparedTuple,
    sign_bits: u32,
    modulus: &BigUint,
) -> BTreeSet<BigUint> {
    let signs = 1u64 << sign_bits;
    let (mut forward_corner, mut reverse_corner) = initial_corners(tuple, modulus);
    let one = BigUint::one();
    let mut support = BTreeSet::new();
    for ordinal in 0..signs {
        let gray = ordinal ^ (ordinal >> 1);
        let forward = ordinal & 1 == 0;
        let mut state = if forward {
            forward_corner.clone()
        } else {
            reverse_corner.clone()
        };
        support.insert(state.clone());
        for _ in 0..tuple.capacity {
            add_mod(&mut state, &one, forward, modulus);
            support.insert(state.clone());
        }
        if ordinal + 1 != signs {
            let next_gray = (ordinal + 1) ^ ((ordinal + 1) >> 1);
            let changed = (gray ^ next_gray).trailing_zeros() as usize;
            let old_positive = gray & (1 << changed) == 0;
            update_corners(
                &mut forward_corner,
                &mut reverse_corner,
                tuple,
                changed,
                old_positive,
                modulus,
            );
        }
    }
    support
}

fn compare_direct(support: &OrbitSupport, direct: &BTreeSet<BigUint>) -> (u64, u64) {
    let intersection = direct
        .iter()
        .filter(|value| interval_contains(&support.merged, value))
        .count() as u64;
    (
        support.exact_distinct_states - intersection,
        direct.len() as u64 - intersection,
    )
}

fn toy_tuple(coefficients: &[u64], capacities: &[u64], prime: u64) -> PreparedTuple {
    PreparedTuple {
        attempt: 0,
        columns: (0..coefficients.len() as u64).collect(),
        slot_capacities: capacities.to_vec(),
        capacity: capacities.iter().sum(),
        original_coefficients: coefficients.iter().copied().map(BigUint::from).collect(),
        endpoint_coefficients: coefficients
            .iter()
            .zip(capacities)
            .map(|(coefficient, capacity)| BigUint::from((coefficient + capacity) % prime))
            .collect(),
    }
}

fn toy_control(prime: u64, coefficients: &[u64], capacities: &[u64]) -> Result<ToyControl, String> {
    let modulus = BigUint::from(prime);
    let tuple = toy_tuple(coefficients, capacities, prime);
    let sign_bits = coefficients.len() as u32;
    let interval = monotone_orbit(&tuple, sign_bits, &modulus)?;
    let direct = direct_monotone_support(&tuple, sign_bits, &modulus);
    let (false_positives, false_negatives) = compare_direct(&interval, &direct);
    Ok(ToyControl {
        prime,
        arity: coefficients.len() as u64,
        sign_bits,
        states_emitted: interval.states_emitted,
        direct_distinct_states: direct.len() as u64,
        interval_distinct_states: interval.exact_distinct_states,
        false_positives,
        false_negatives,
        exact: false_positives == 0 && false_negatives == 0,
    })
}

fn tuple_digest(tuple: &PreparedTuple) -> String {
    let mut bytes = Vec::with_capacity(ARITY * 16);
    for (column, capacity) in tuple.columns.iter().zip(&tuple.slot_capacities) {
        bytes.extend(column.to_be_bytes());
        bytes.extend(capacity.to_be_bytes());
    }
    hex::encode(sha256(&bytes))
}

fn charged_additions(capacity: u64) -> u64 {
    SIGNS * capacity + (SIGNS - 1) + (ARITY - 1) as u64 + ARITY as u64 + SIGNS
}

fn boundary_receipt(
    rare_edges: u64,
    charged: u64,
    distinct: u64,
    modulus: &BigUint,
) -> Result<BoundaryReceipt, String> {
    let rare_fraction = rare_edges as f64 / COLUMNS as f64;
    let common_fraction = 1.0 - rare_fraction;
    let hidden_collision_probability =
        common_fraction * common_fraction + rare_fraction * rare_fraction;
    let hidden_ratio_to_rho = LOCAL_RATIO / hidden_collision_probability.sqrt();
    let support_unclamped = ARITY as f64 * (2.0 * rare_edges as f64).log2()
        + ((2 * ARITY as u64 * COLUMNS + 1) as f64).log2();
    let log_group_order = log2_big(modulus)?;
    let support_upper_log2 = support_unclamped.min(log_group_order);
    let support_upper_fraction_of_group = 2f64.powf((support_unclamped - log_group_order).min(0.0));
    let target_retry_lower_bound = 1.0 / support_upper_fraction_of_group;
    let factor_base_complete_ratio_lower_bound = hidden_ratio_to_rho * target_retry_lower_bound;
    Ok(BoundaryReceipt {
        rare_fraction,
        hidden_collision_probability,
        hidden_ratio_to_rho,
        support_upper_log2,
        support_upper_fraction_of_group,
        target_retry_lower_bound,
        factor_base_complete_ratio_lower_bound,
        ideal_local_ratio_per_distinct: LOCAL_RATIO * charged as f64 / distinct as f64,
        factor_base_corrected_ratio_per_distinct: factor_base_complete_ratio_lower_bound
            * charged as f64
            / distinct as f64,
    })
}

fn round_receipt(round: u32, modulus: &BigUint) -> Result<RoundReceipt, String> {
    let total_started = Instant::now();
    let cpu_started = process_cpu_ns();
    let rare_edges = screening_rare(round)?;
    let factor_base = construct_factor_base(rare_edges, modulus)?;
    let tuple_started = Instant::now();
    let tuple = select_tuple(rare_edges, &factor_base.coefficients, modulus)?;
    let tuple_wall_ns = tuple_started.elapsed().as_nanos() as u64;
    let support = monotone_orbit(&tuple, SIGN_BITS, modulus)?;
    let charged = charged_additions(tuple.capacity);
    let boundary = boundary_receipt(rare_edges, charged, support.exact_distinct_states, modulus)?;
    let distinct_fraction = support.exact_distinct_states as f64 / support.states_emitted as f64;
    let selector_gate_passed = distinct_fraction >= ROUND292_GATE_FRACTION
        && boundary.ideal_local_ratio_per_distinct < 1.0
        && boundary.factor_base_corrected_ratio_per_distinct < 1.0;
    let exact_factor_base = factor_base.distinct_relative_coefficients == COLUMNS
        && factor_base.unique_up_to_sign
        && factor_base.nonidentity
        && factor_base.edge_replay_failures == 0;
    let exact_coverage = support.states_emitted == SIGNS * (tuple.capacity + 1)
        && support.exact_distinct_states <= support.states_emitted;
    let logical_materialized_bytes =
        COLUMNS * 32 + support.source_intervals * 64 + support.merged_intervals * 64 + COLUMNS * 8;
    let timing = RoundTiming {
        factor_base_wall_ns: factor_base.construction_wall_ns,
        tuple_wall_ns,
        orbit_generation_wall_ns: support.generation_wall_ns,
        interval_union_wall_ns: support.union_wall_ns,
        total_wall_ns: total_started.elapsed().as_nanos() as u64,
        process_cpu_ns: process_cpu_ns().saturating_sub(cpu_started),
        peak_rss_bytes: peak_rss_bytes(),
    };
    let operations = OperationCounts {
        coefficients_constructed: COLUMNS,
        factor_base_edges_replayed: COLUMNS,
        signed_start_initial_terms: 2 * ARITY as u64,
        gray_boundary_updates: SIGNS - 1,
        path_additions: SIGNS * tuple.capacity,
        sign_boundary_additions: SIGNS - 1,
        initial_sum_additions: (ARITY - 1) as u64,
        boundary_precomputation_additions: ARITY as u64,
        target_correction_additions: SIGNS,
        charged_addition_equivalents: charged,
        sort_comparisons: support.sort_comparisons,
        merge_comparisons: support.merge_comparisons,
    };
    Ok(RoundReceipt {
        schema: "p256-executed-factor-base-screen-round/v1".into(),
        curve: CURVE.into(),
        comparison_factor_base: COMPARISON_FB.into(),
        screening_round: round,
        execution_status: "complete".into(),
        factor_base: FactorBaseReceipt {
            id: factor_base.candidate_id,
            columns: COLUMNS,
            rare_edges,
            common_edges: factor_base.common_edges,
            delta_common: "1".into(),
            delta_rare: factor_base.rare_delta.to_string(),
            anchor_attempts: factor_base.anchor_attempts,
            anchor_coefficient: factor_base.anchor.to_string(),
            coefficient_sha256: factor_base.coefficient_sha256,
            coefficient_bytes: COLUMNS * 32,
            primitive_cycle: rare_edges.gcd(&COLUMNS) == 1,
            distinct_relative_coefficients: factor_base.distinct_relative_coefficients,
            unique_up_to_sign: factor_base.unique_up_to_sign,
            nonidentity: factor_base.nonidentity,
            edge_replay_failures: factor_base.edge_replay_failures,
        },
        tuple: TupleReceipt {
            attempt: tuple.attempt,
            columns: tuple.columns.clone(),
            slot_capacities: tuple.slot_capacities.clone(),
            total_capacity: tuple.capacity,
            tuple_sha256: tuple_digest(&tuple),
        },
        coverage: CoverageReceipt {
            sign_bits: SIGN_BITS,
            sign_segments: SIGNS,
            states_emitted: support.states_emitted,
            per_segment_distinct_sum: support.states_emitted,
            exact_union_distinct_states: support.exact_distinct_states,
            within_segment_duplicates: 0,
            cross_segment_duplicates: support.states_emitted - support.exact_distinct_states,
            total_duplicates: support.states_emitted - support.exact_distinct_states,
            distinct_fraction,
            source_intervals: support.source_intervals,
            merged_intervals: support.merged_intervals,
            wrap_splits: support.wrap_splits,
            support_digest: support.support_digest,
            start_digest: support.start_digest,
        },
        boundary,
        operations,
        timing,
        logical_materialized_bytes,
        exact_factor_base,
        exact_coverage,
        selector_gate_passed,
        structured_residual_degree: None,
        relations_reported: 0,
        full_depth_unplanted_relation_attempted: false,
        interpretation: "This is a complete full-sign-depth factor-base coverage screen. It is not a summation-polynomial solve, relation, row collection, sparse linear-algebra run, or complete DLP.".into(),
    })
}

fn write_round(receipt: &RoundReceipt, directory: &Path) -> Result<ManifestEntry, String> {
    let filename = format!(
        "round{:03}-r{:05}.json",
        receipt.screening_round, receipt.factor_base.rare_edges
    );
    let path = directory.join(filename);
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    let bytes = serde_json::to_vec_pretty(receipt).map_err(|error| error.to_string())?;
    fs::write(&path, &bytes).map_err(|error| format!("{}: {error}", path.display()))?;
    Ok(ManifestEntry {
        round: receipt.screening_round,
        rare_edges: receipt.factor_base.rare_edges,
        factor_base_id: receipt.factor_base.id.clone(),
        path: path.display().to_string(),
        bytes: bytes.len() as u64,
        sha256: hex::encode(sha256(&bytes)),
        execution_status: receipt.execution_status.clone(),
        exact_distinct_states: receipt.coverage.exact_union_distinct_states,
        distinct_fraction: receipt.coverage.distinct_fraction,
        ideal_local_ratio_per_distinct: receipt.boundary.ideal_local_ratio_per_distinct,
        factor_base_corrected_ratio_per_distinct: receipt
            .boundary
            .factor_base_corrected_ratio_per_distinct,
        selector_gate_passed: receipt.selector_gate_passed,
    })
}

fn replay_manifest_entries(entries: &[ManifestEntry]) -> Result<bool, String> {
    for entry in entries {
        let bytes = fs::read(&entry.path).map_err(|error| format!("{}: {error}", entry.path))?;
        if bytes.len() as u64 != entry.bytes || hex::encode(sha256(&bytes)) != entry.sha256 {
            return Ok(false);
        }
        let receipt: Value = serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
        if receipt.get("screening_round").and_then(Value::as_u64) != Some(u64::from(entry.round))
            || receipt.get("execution_status").and_then(Value::as_str) != Some("complete")
        {
            return Ok(false);
        }
    }
    Ok(true)
}

fn native_controls(modulus: &BigUint) -> Result<Vec<NativeControl>, String> {
    let factor_base = construct_factor_base(6_935, modulus)?;
    if factor_base.coefficient_sha256 != ROUND31_R6935_COEFFICIENT_SHA256 {
        return Err("R=6935 coefficient digest does not match Round 31".into());
    }
    let tuple = select_tuple(6_935, &factor_base.coefficients, modulus)?;
    if tuple.attempt != 240 || tuple.capacity != 227 {
        return Err("R=6935 tuple does not match the frozen Round 291 primary".into());
    }
    let mut controls = Vec::new();
    for sign_bits in [8, 10, 12] {
        let interval = monotone_orbit(&tuple, sign_bits, modulus)?;
        let direct = direct_monotone_support(&tuple, sign_bits, modulus);
        let (false_positives, false_negatives) = compare_direct(&interval, &direct);
        controls.push(NativeControl {
            rare_edges: 6_935,
            sign_bits,
            states_emitted: interval.states_emitted,
            direct_distinct_states: direct.len() as u64,
            interval_distinct_states: interval.exact_distinct_states,
            false_positives,
            false_negatives,
            exact: false_positives == 0 && false_negatives == 0,
        });
    }
    Ok(controls)
}

fn run(cli: Cli) -> Result<(), String> {
    if cli.out.exists() {
        return Err(format!("refusing to overwrite {}", cli.out.display()));
    }
    fs::create_dir_all(&cli.round_dir)
        .map_err(|error| format!("{}: {error}", cli.round_dir.display()))?;
    if fs::read_dir(&cli.round_dir)
        .map_err(|error| error.to_string())?
        .next()
        .is_some()
    {
        return Err(format!(
            "round directory {} is not empty",
            cli.round_dir.display()
        ));
    }

    let (dep_round1_fb, round1_fb) = load_dependency(
        &cli.round1_factor_base,
        ROUND1_FB_SHA256,
        "p256.affine_bitbox_result/v1",
    )?;
    let (dep_round1_degree, _) = load_dependency(
        &cli.round1_degree,
        ROUND1_DEGREE_SHA256,
        "p256.factor_base_degree/v1",
    )?;
    let (dep_round292, round292) = load_dependency(
        &cli.round292,
        ROUND292_SHA256,
        "p256-serpentine-orbit-coverage-round292-v1",
    )?;
    if round1_fb.get("verified").and_then(Value::as_bool) != Some(true)
        || round292
            .get("gates")
            .and_then(|value| value.get("promoted"))
            .and_then(Value::as_bool)
            != Some(false)
    {
        return Err("dependency gate mismatch".into());
    }

    let toy_controls = vec![
        toy_control(257, &[3, 29, 101], &[3, 4, 2])?,
        toy_control(65_537, &[7, 101, 9_997, 31_337, 50_003], &[2, 5, 3, 4, 6])?,
    ];
    if toy_controls.iter().any(|control| !control.exact) {
        return Err("toy monotone support control failed".into());
    }
    let curve = CurveParams::p256();
    let native_controls = native_controls(&curve.n)?;
    if native_controls.iter().any(|control| !control.exact) {
        return Err("native monotone support control failed".into());
    }

    let sweep_started = Instant::now();
    let sweep_cpu_started = process_cpu_ns();
    let mut entries = Vec::with_capacity(ROUND_COUNT);
    let mut every_factor_base_exact = true;
    let mut all_complete = true;
    let mut storage_below = true;
    let mut maximum_peak_rss_bytes = None::<u64>;
    for round in ROUND_FIRST..=ROUND_LAST {
        let receipt = round_receipt(round, &curve.n)?;
        every_factor_base_exact &= receipt.exact_factor_base && receipt.exact_coverage;
        all_complete &= receipt.execution_status == "complete";
        storage_below &= receipt.logical_materialized_bytes < (1u64 << 50);
        maximum_peak_rss_bytes = match (maximum_peak_rss_bytes, receipt.timing.peak_rss_bytes) {
            (Some(left), Some(right)) => Some(left.max(right)),
            (None, right) => right,
            (left, None) => left,
        };
        entries.push(write_round(&receipt, &cli.round_dir)?);
    }
    let total_wall_ns = sweep_started.elapsed().as_nanos() as u64;
    let total_process_cpu_ns = process_cpu_ns().saturating_sub(sweep_cpu_started);
    let hashes_replayed = replay_manifest_entries(&entries)?;
    let exact_range = entries.len() == ROUND_COUNT
        && entries
            .iter()
            .enumerate()
            .all(|(index, entry)| entry.round == ROUND_FIRST + index as u32);
    let controls_exact = toy_controls.iter().all(|control| {
        control.false_positives == 0 && control.false_negatives == 0 && control.exact
    }) && native_controls.iter().all(|control| {
        control.false_positives == 0 && control.false_negatives == 0 && control.exact
    });
    let any_selector_gate_passed = entries.iter().any(|entry| entry.selector_gate_passed);
    let execution_gate_passed = exact_range
        && all_complete
        && every_factor_base_exact
        && controls_exact
        && hashes_replayed
        && storage_below;
    let gates = ExecutionGates {
        exactly_289_rounds_2_through_290_once: exact_range,
        every_round_complete: all_complete,
        every_factor_base_exact_and_primitive: every_factor_base_exact,
        zero_false_positives_and_false_negatives: controls_exact,
        manifest_hashes_and_byte_counts_replayed: hashes_replayed,
        no_projection_counted_as_execution: entries.iter().all(|entry| {
            entry.execution_status == "complete" && !entry.path.contains("projection")
        }),
        peak_materialized_storage_below_2_50: storage_below,
        execution_gate_passed,
        any_selector_gate_passed,
        same_family_structured_degree_at_most_5: false,
        per_usable_relation_below_2_103: false,
        relation_collection_below_2_120: false,
        complete_non_generic_dlp_below_rho: false,
        promoted: false,
    };
    let best_ideal = entries
        .iter()
        .min_by(|left, right| {
            left.ideal_local_ratio_per_distinct
                .total_cmp(&right.ideal_local_ratio_per_distinct)
                .then_with(|| left.round.cmp(&right.round))
        })
        .ok_or("empty executed sweep")?;
    let best_corrected = entries
        .iter()
        .min_by(|left, right| {
            left.factor_base_corrected_ratio_per_distinct
                .total_cmp(&right.factor_base_corrected_ratio_per_distinct)
                .then_with(|| left.round.cmp(&right.round))
        })
        .ok_or("empty executed sweep")?;
    let best_ideal_round = best_ideal.round;
    let best_corrected_round = best_corrected.round;
    let total_artifact_bytes = entries.iter().map(|entry| entry.bytes).sum::<u64>();
    let semantic = serde_json::to_vec(&serde_json::json!({
        "toy_controls": &toy_controls,
        "native_controls": &native_controls,
        "entries": &entries,
        "gates": &gates,
        "best_ideal": best_ideal_round,
        "best_corrected": best_corrected_round,
    }))
    .map_err(|error| error.to_string())?;
    let decision = if any_selector_gate_passed {
        "At least one factor base passed the registered selector gate. Retain only those candidates for same-family degree and relation testing; do not claim end-to-end parity."
    } else {
        "All 289 requested rounds executed, but no factor base passed the selector gate. Publish the best measured candidate as a negative result; do not attempt an unplanted full-depth relation or claim rho parity."
    };
    let manifest = Manifest {
        schema: "p256-executed-factor-base-screen-round2-290/v1".into(),
        curve: CURVE.into(),
        family: "289 exact known-log two-delta monotone-serpentine factor-base screens".into(),
        comparison_factor_base: COMPARISON_FB.into(),
        dependencies: vec![dep_round1_fb, dep_round1_degree, dep_round292],
        round_range: [ROUND_FIRST, ROUND_LAST],
        registered_rounds: ROUND_COUNT as u64,
        executed_rounds: entries.len() as u64,
        projection_only_rounds: 0,
        frozen: serde_json::json!({
            "columns": COLUMNS,
            "arity": ARITY,
            "cutoff": CUTOFF,
            "sign_bits": SIGN_BITS,
            "rare_first": RARE_FIRST,
            "rare_last": RARE_LAST,
            "rare_step": 2,
            "round_formula": "R(r)=6647+2*(r-2)",
            "local_ideal_oracle_ratio": LOCAL_RATIO,
            "distinct_fraction_gate": ROUND292_GATE_FRACTION,
        }),
        toy_controls,
        native_controls,
        entries,
        best_ideal_local_round: best_ideal_round,
        best_factor_base_corrected_round: best_corrected_round,
        total_wall_ns,
        total_process_cpu_ns,
        maximum_peak_rss_bytes,
        total_artifact_bytes,
        gates,
        semantic_evidence_sha256: hex::encode(sha256(&semantic)),
        relations_reported: 0,
        full_depth_unplanted_relation_attempted: false,
        decision: decision.into(),
    };
    let bytes = serde_json::to_vec_pretty(&manifest).map_err(|error| error.to_string())?;
    fs::write(&cli.out, &bytes).map_err(|error| format!("{}: {error}", cli.out.display()))?;
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("P-256 executed screen failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_map_has_289_distinct_odd_candidates() {
        let rare = (ROUND_FIRST..=ROUND_LAST)
            .map(|round| screening_rare(round).unwrap())
            .collect::<Vec<_>>();
        assert_eq!(rare.len(), 289);
        assert_eq!(rare.first(), Some(&RARE_FIRST));
        assert_eq!(rare.last(), Some(&RARE_LAST));
        assert!(rare.windows(2).all(|pair| pair[1] - pair[0] == 2));
        assert!(rare.iter().all(|value| value.gcd(&COLUMNS) == 1));
    }

    #[test]
    fn monotone_toy_support_is_exact() {
        let first = toy_control(257, &[3, 29, 101], &[3, 4, 2]).unwrap();
        let second =
            toy_control(65_537, &[7, 101, 9_997, 31_337, 50_003], &[2, 5, 3, 4, 6]).unwrap();
        assert!(first.exact);
        assert!(second.exact);
        assert_eq!(first.false_positives + first.false_negatives, 0);
        assert_eq!(second.false_positives + second.false_negatives, 0);
    }

    #[test]
    fn charged_formula_contains_every_registered_addition() {
        let capacity = 227;
        let expected = SIGNS * capacity + (SIGNS - 1) + 16 + 17 + SIGNS;
        assert_eq!(charged_additions(capacity), expected);
    }
}
