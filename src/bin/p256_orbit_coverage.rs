//! Exact support accounting for P-256 serpentine sign orbits (round 292).

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
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;
use serde_json::Value;

const CURVE: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const COMPARISON_FB: &str = "FB1h2f8621cda105";
const ROUND291_SHA256: &str = "068be64f508ea2fb2996c265426dd8213088d1cad905859470792fcbc80972da";
const ROUND291_SEMANTIC_SHA256: &str =
    "453d14bc5a9e68ab5c6ba9d73ffa9040bc3c8557be7d73b0eb42d5395d939710";
const ROUND31_SHA256: &str = "9b9ee16c5a24e868d1b9304f08aa80ca39ac23b1b4586ef730534e739172c435";
const COEFFICIENT_SHA256: &str = "980917981827d813e60484abb0655e8bd527b0beb540ea146d2974ff17303a53";
const COLUMNS: u64 = 131_458;
const RARE_EDGES: u64 = 6_935;
const ARITY: usize = 17;
const CUTOFF: u64 = 219;
const FULL_SIGN_BITS: u32 = 17;
const LOCAL_RATIO: f64 = 0.964_336_477_130_181;
const ROUND291_STAGE_RATIO: f64 = 0.968_617_146_446_504_4;
const PRIMARY_DEPTHS: [u32; 6] = [8, 10, 12, 14, 16, 17];
const EXPECTED_ATTEMPTS: [u64; 4] = [240, 1648, 1766, 1920];
const DIRECT_CONTROL_MAX_DEPTH: u32 = 12;

#[derive(Parser)]
#[command(about = "Measure exact P-256 serpentine-orbit support (round 292)")]
struct Cli {
    #[arg(long)]
    round291: PathBuf,
    #[arg(long)]
    round31: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
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
    capacity: u64,
    order: Vec<usize>,
    original_coefficients: Vec<BigUint>,
    endpoint_coefficients: Vec<BigUint>,
}

struct OrbitIntervals {
    attempt: u64,
    capacity: u64,
    sign_bits: u32,
    intervals: Vec<Interval>,
    emitted_states: u64,
    interval_width_sum: u64,
    maximum_interval_width: u64,
    wrap_splits: u64,
    prefix_steps: u64,
    signed_sum_terms: u64,
    generation_wall_ns: u64,
    boundary_digest: String,
}

struct UnionSummary {
    merged: Vec<Interval>,
    width: u64,
    sort_comparisons: u64,
    merge_comparisons: u64,
    wall_ns: u64,
    digest: String,
}

#[derive(Clone, Serialize)]
struct ToyCell {
    prime: u64,
    arity: u64,
    sign_segments: u64,
    states_emitted: u64,
    direct_distinct_states: u64,
    interval_distinct_states: u64,
    direct_duplicate_states: u64,
    interval_duplicate_states: u64,
    false_positives: u64,
    false_negatives: u64,
    exact: bool,
}

#[derive(Clone, Serialize)]
struct CoverageCell {
    class: String,
    tuple_count: u64,
    tuple_attempts: Vec<u64>,
    capacities: Vec<u64>,
    sign_bits: u32,
    sign_segments: u64,
    states_emitted: u64,
    per_segment_distinct_sum: u64,
    exact_union_distinct_states: u64,
    within_segment_duplicates: u64,
    cross_segment_duplicates: u64,
    total_duplicates: u64,
    distinct_fraction: f64,
    maximum_interval_width: u64,
    source_intervals: u64,
    merged_intervals: u64,
    wrap_splits: u64,
    charged_addition_equivalents: u64,
    uncorrected_actual_stage_ratio: f64,
    corrected_actual_stage_ratio_per_distinct: f64,
    corrected_frozen_round291_ratio_per_distinct: f64,
    prefix_steps: u64,
    signed_sum_terms: u64,
    sort_comparisons: u64,
    merge_comparisons: u64,
    logical_interval_endpoint_bytes: u64,
    generation_wall_ns: u64,
    sort_merge_wall_ns: u64,
    peak_rss_bytes: Option<u64>,
    support_digest: String,
    boundary_digests: Vec<String>,
    round291_boundary_digests_match: Option<bool>,
    direct_control_complete: bool,
    direct_distinct_states: Option<u64>,
    false_positives: Option<u64>,
    false_negatives: Option<u64>,
}

#[derive(Clone, Serialize)]
struct Fits {
    primary_depths: Vec<u32>,
    emitted_state_exponent: f64,
    distinct_state_exponent: f64,
    interval_memory_exponent: f64,
    interpretation: String,
}

#[derive(Clone, Serialize)]
struct Gates {
    zero_false_positives_and_false_negatives: bool,
    exact_toy_support_and_duplicates: bool,
    native_boundary_digests_match_round291: bool,
    full_depth_primary_distinct_fraction_at_least_round291_stage_ratio: bool,
    full_depth_primary_corrected_stage_ratio_at_most_one: bool,
    four_tuple_corrected_stage_ratio_at_most_one: bool,
    materialized_storage_below_2_50: bool,
    coverage_hypothesis_passed: bool,
    same_family_structured_degree_at_most_5: bool,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    complete_non_generic_dlp_below_rho: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultRecord {
    schema: String,
    curve: String,
    round: u32,
    family: String,
    comparison_factor_base: String,
    dependencies: Vec<Dependency>,
    frozen: Value,
    exact_interval_argument: String,
    toy_cells: Vec<ToyCell>,
    primary_cells: Vec<CoverageCell>,
    portfolio_cells: Vec<CoverageCell>,
    fits: Fits,
    gates: Gates,
    semantic_evidence_sha256: String,
    relations_reported: u64,
    full_depth_unplanted_relation_attempted: bool,
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

fn json_bool(value: &Value, path: &[&str]) -> Result<bool, String> {
    json_at(value, path)?
        .as_bool()
        .ok_or_else(|| format!("JSON path {} is not a bool", path.join(".")))
}

fn load_json(path: &Path, expected_sha: &str, schema: &str) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != expected_sha {
        return Err(format!(
            "{} hash mismatch: expected {expected_sha}, got {digest}",
            path.display()
        ));
    }
    let value: Value = serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if json_str(&value, &["schema"])? != schema || json_str(&value, &["curve"])? != CURVE {
        return Err(format!("{} identity mismatch", path.display()));
    }
    Ok((
        Dependency {
            path: path.display().to_string(),
            sha256: digest,
            schema: schema.into(),
        },
        value,
    ))
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    BigUint::parse_bytes(value.as_bytes(), 10)
        .ok_or_else(|| format!("invalid frozen integer: {value}"))
}

fn fixed_32(value: &BigUint) -> Result<[u8; 32], String> {
    let raw = value.to_bytes_be();
    if raw.len() > 32 {
        return Err("coefficient exceeded 32 bytes".into());
    }
    let mut out = [0u8; 32];
    out[32 - raw.len()..].copy_from_slice(&raw);
    Ok(out)
}

fn mechanical_rare(index: u64, columns: u64, rare_edges: u64) -> bool {
    (index + 1) * rare_edges / columns > index * rare_edges / columns
}

fn residual_capacities(columns: u64, rare_edges: u64) -> Vec<u64> {
    (0..columns)
        .map(|index| {
            let mut distance = 0u64;
            while !mechanical_rare((index + distance) % columns, columns, rare_edges) {
                distance += 1;
                assert!(distance <= columns);
            }
            distance
        })
        .collect()
}

fn forward_order(columns: &[u64], capacities: &[u64]) -> Result<Vec<usize>, String> {
    let mut positions = columns.to_vec();
    let mut remaining = columns
        .iter()
        .map(|column| capacities[*column as usize])
        .collect::<Vec<_>>();
    let expected = remaining.iter().sum::<u64>();
    let mut order = Vec::with_capacity(expected as usize);
    while let Some(slot) = (0..columns.len())
        .filter(|slot| remaining[*slot] != 0)
        .min_by_key(|slot| positions[*slot])
    {
        order.push(slot);
        positions[slot] = (positions[slot] + 1) % capacities.len() as u64;
        remaining[slot] -= 1;
    }
    if order.len() as u64 != expected {
        return Err("forward order length mismatch".into());
    }
    Ok(order)
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

fn native_candidate(round31: &Value) -> Result<&Value, String> {
    json_at(round31, &["native_candidates"])?
        .as_array()
        .ok_or("round-31 native_candidates is not an array")?
        .iter()
        .find(|row| row.get("rare_edges").and_then(Value::as_u64) == Some(RARE_EDGES))
        .ok_or_else(|| "round-31 R=6935 candidate missing".into())
}

fn native_coefficients(round31: &Value, modulus: &BigUint) -> Result<Vec<BigUint>, String> {
    let candidate = native_candidate(round31)?;
    let common = parse_big(json_str(candidate, &["delta_common"])?)?;
    let rare = parse_big(json_str(candidate, &["delta_rare"])?)?;
    let anchor = parse_big(json_str(candidate, &["anchor_coefficient"])?)?;
    if common != BigUint::one() {
        return Err("registered common delta is not one".into());
    }
    let mut coefficients = Vec::with_capacity(COLUMNS as usize);
    let mut relative = BigUint::zero();
    for index in 0..COLUMNS {
        coefficients.push((&relative + &anchor) % modulus);
        relative = (&relative
            + if mechanical_rare(index, COLUMNS, RARE_EDGES) {
                &rare
            } else {
                &common
            })
            % modulus;
    }
    let mut bytes = Vec::with_capacity(COLUMNS as usize * 32);
    for coefficient in &coefficients {
        bytes.extend(fixed_32(coefficient)?);
    }
    let digest = hex::encode(sha256(&bytes));
    if digest != COEFFICIENT_SHA256 {
        return Err(format!("coefficient digest mismatch: {digest}"));
    }
    Ok(coefficients)
}

fn verify_round291(round291: &Value) -> Result<(), String> {
    if json_u64(round291, &["round"])? != 291
        || json_str(round291, &["comparison_factor_base"])? != COMPARISON_FB
        || json_str(round291, &["semantic_evidence_sha256"])? != ROUND291_SEMANTIC_SHA256
        || json_u64(round291, &["frozen", "columns"])? != COLUMNS
        || json_u64(round291, &["frozen", "arity"])? != ARITY as u64
        || json_u64(round291, &["frozen", "rare_edges"])? != RARE_EDGES
        || json_u64(round291, &["frozen", "cutoff"])? != CUTOFF
        || json_str(round291, &["frozen", "coefficient_sha256"])? != COEFFICIENT_SHA256
        || !json_bool(round291, &["gates", "counted_stage_below_rho"])?
        || json_bool(round291, &["gates", "promoted"])?
    {
        return Err("round-291 frozen identity mismatch".into());
    }
    let ratio = json_at(round291, &["accounting", "conservative_stage_ratio_to_rho"])?
        .as_f64()
        .ok_or("round-291 stage ratio is not f64")?;
    if ratio.to_bits() != ROUND291_STAGE_RATIO.to_bits() {
        return Err("round-291 stage ratio mismatch".into());
    }
    let dependencies = json_at(round291, &["dependencies"])?
        .as_array()
        .ok_or("round-291 dependencies is not an array")?;
    if !dependencies
        .iter()
        .any(|dependency| dependency.get("sha256").and_then(Value::as_str) == Some(ROUND31_SHA256))
    {
        return Err("round-291 does not bind the frozen round-31 source".into());
    }
    Ok(())
}

fn parse_tuples(
    round291: &Value,
    coefficients: &[BigUint],
    capacities: &[u64],
    modulus: &BigUint,
) -> Result<Vec<PreparedTuple>, String> {
    let cells = json_at(round291, &["native_cells"])?
        .as_array()
        .ok_or("round-291 native_cells is not an array")?;
    let mut tuples = Vec::new();
    for expected_attempt in EXPECTED_ATTEMPTS {
        let cell = cells
            .iter()
            .find(|cell| {
                cell.get("tuple_attempt").and_then(Value::as_u64) == Some(expected_attempt)
            })
            .ok_or_else(|| format!("round-291 tuple attempt {expected_attempt} missing"))?;
        let columns = json_at(cell, &["columns"])?
            .as_array()
            .ok_or("round-291 columns is not an array")?
            .iter()
            .map(|value| value.as_u64().ok_or("round-291 column is not u64"))
            .collect::<Result<Vec<_>, _>>()?;
        if columns != draw_columns(expected_attempt) {
            return Err(format!(
                "tuple attempt {expected_attempt} column draw mismatch"
            ));
        }
        let capacity = columns
            .iter()
            .map(|column| capacities[*column as usize])
            .sum::<u64>();
        if capacity != json_u64(cell, &["capacity"])? || capacity < CUTOFF {
            return Err(format!(
                "tuple attempt {expected_attempt} capacity mismatch"
            ));
        }
        let order = forward_order(&columns, capacities)?;
        let mut original_coefficients = Vec::with_capacity(ARITY);
        let mut endpoint_coefficients = Vec::with_capacity(ARITY);
        for column in &columns {
            let original = coefficients[*column as usize].clone();
            let endpoint_column = (column + capacities[*column as usize]) % COLUMNS;
            let endpoint = coefficients[endpoint_column as usize].clone();
            let mut expected = original.clone();
            add_mod(
                &mut expected,
                &BigUint::from(capacities[*column as usize]),
                true,
                modulus,
            );
            if endpoint != expected {
                return Err(format!(
                    "tuple attempt {expected_attempt} endpoint mismatch"
                ));
            }
            original_coefficients.push(original);
            endpoint_coefficients.push(endpoint);
        }
        tuples.push(PreparedTuple {
            attempt: expected_attempt,
            capacity,
            order,
            original_coefficients,
            endpoint_coefficients,
        });
    }
    Ok(tuples)
}

fn signed_start(coefficients: &[BigUint], gray: u64, modulus: &BigUint) -> BigUint {
    let mut sum = BigUint::zero();
    for (slot, coefficient) in coefficients.iter().enumerate() {
        add_mod(&mut sum, coefficient, gray & (1 << slot) == 0, modulus);
    }
    sum
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

fn split_modular_interval(
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

fn path_extrema(tuple: &PreparedTuple, gray: u64, forward: bool) -> (i64, i64) {
    let mut prefix = 0i64;
    let mut minimum = 0i64;
    let mut maximum = 0i64;
    if forward {
        for slot in &tuple.order {
            prefix += if gray & (1 << slot) == 0 { 1 } else { -1 };
            minimum = minimum.min(prefix);
            maximum = maximum.max(prefix);
        }
    } else {
        for slot in tuple.order.iter().rev() {
            prefix += if gray & (1 << slot) == 0 { -1 } else { 1 };
            minimum = minimum.min(prefix);
            maximum = maximum.max(prefix);
        }
    }
    (minimum, maximum)
}

fn generate_orbit(
    tuple: &PreparedTuple,
    sign_bits: u32,
    modulus: &BigUint,
) -> Result<OrbitIntervals, String> {
    let started = Instant::now();
    let signs = 1u64 << sign_bits;
    let mut intervals = Vec::with_capacity(signs as usize + 2);
    let mut interval_width_sum = 0u64;
    let mut maximum_interval_width = 0u64;
    let mut wrap_splits = 0u64;
    let mut digest = Hasher::new();
    for ordinal in 0..signs {
        let gray = ordinal ^ (ordinal >> 1);
        let forward = ordinal & 1 == 0;
        let coefficients = if forward {
            &tuple.original_coefficients
        } else {
            &tuple.endpoint_coefficients
        };
        let start = signed_start(coefficients, gray, modulus);
        digest.update(&ordinal.to_be_bytes());
        digest.update(&gray.to_be_bytes());
        digest.update(&fixed_32(&start)?);
        let (minimum, maximum) = path_extrema(tuple, gray, forward);
        let interval_width = (maximum - minimum + 1) as u64;
        interval_width_sum += interval_width;
        maximum_interval_width = maximum_interval_width.max(interval_width);
        let (mut pieces, wrapped) = split_modular_interval(&start, minimum, maximum, modulus);
        wrap_splits += u64::from(wrapped);
        intervals.append(&mut pieces);
    }
    Ok(OrbitIntervals {
        attempt: tuple.attempt,
        capacity: tuple.capacity,
        sign_bits,
        intervals,
        emitted_states: signs * (tuple.capacity + 1),
        interval_width_sum,
        maximum_interval_width,
        wrap_splits,
        prefix_steps: signs * tuple.capacity,
        signed_sum_terms: signs * ARITY as u64,
        generation_wall_ns: started.elapsed().as_nanos() as u64,
        boundary_digest: digest.finalize().to_hex().to_string(),
    })
}

fn interval_width(interval: &Interval) -> Result<u64, String> {
    (&interval.end - &interval.start + BigUint::one())
        .to_u64()
        .ok_or_else(|| "interval width does not fit u64".into())
}

fn union_intervals(mut intervals: Vec<Interval>) -> Result<UnionSummary, String> {
    let started = Instant::now();
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
    let width = merged
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
    Ok(UnionSummary {
        merged,
        width,
        sort_comparisons,
        merge_comparisons,
        wall_ns: started.elapsed().as_nanos() as u64,
        digest: digest.finalize().to_hex().to_string(),
    })
}

fn interval_contains(intervals: &[Interval], value: &BigUint) -> bool {
    let index = intervals.partition_point(|interval| interval.end < *value);
    index < intervals.len() && intervals[index].start <= *value
}

fn direct_support(tuple: &PreparedTuple, sign_bits: u32, modulus: &BigUint) -> BTreeSet<BigUint> {
    let signs = 1u64 << sign_bits;
    let mut support = BTreeSet::new();
    let one = BigUint::one();
    for ordinal in 0..signs {
        let gray = ordinal ^ (ordinal >> 1);
        let forward = ordinal & 1 == 0;
        let coefficients = if forward {
            &tuple.original_coefficients
        } else {
            &tuple.endpoint_coefficients
        };
        let mut coefficient = signed_start(coefficients, gray, modulus);
        support.insert(coefficient.clone());
        if forward {
            for slot in &tuple.order {
                add_mod(&mut coefficient, &one, gray & (1 << slot) == 0, modulus);
                support.insert(coefficient.clone());
            }
        } else {
            for slot in tuple.order.iter().rev() {
                add_mod(&mut coefficient, &one, gray & (1 << slot) != 0, modulus);
                support.insert(coefficient.clone());
            }
        }
    }
    support
}

fn expected_boundary_digest(round291: &Value, attempt: u64, sign_bits: u32) -> Option<&str> {
    json_at(round291, &["native_cells"])
        .ok()?
        .as_array()?
        .iter()
        .find(|cell| {
            cell.get("tuple_attempt").and_then(Value::as_u64) == Some(attempt)
                && cell.get("sign_bits").and_then(Value::as_u64) == Some(u64::from(sign_bits))
        })?
        .get("boundary_digest")?
        .as_str()
}

fn charged_additions(capacity: u64, sign_bits: u32) -> u64 {
    let signs = 1u64 << sign_bits;
    signs * capacity + (signs - 1) + (ARITY - 1) as u64 + 2 * ARITY as u64 + signs
}

fn peak_rss_bytes() -> Option<u64> {
    let status = fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|line| line.starts_with("VmHWM:"))?;
    let kib = line.split_whitespace().nth(1)?.parse::<u64>().ok()?;
    Some(kib * 1024)
}

fn coverage_cell(
    class: &str,
    orbits: &[&OrbitIntervals],
    round291: &Value,
    direct: Option<&BTreeSet<BigUint>>,
) -> Result<CoverageCell, String> {
    let mut intervals = Vec::new();
    let mut states_emitted = 0u64;
    let mut interval_width_sum = 0u64;
    let mut maximum_interval_width = 0u64;
    let mut wrap_splits = 0u64;
    let mut prefix_steps = 0u64;
    let mut signed_sum_terms = 0u64;
    let mut generation_wall_ns = 0u64;
    let mut charged = 0u64;
    for orbit in orbits {
        intervals.extend(orbit.intervals.iter().cloned());
        states_emitted += orbit.emitted_states;
        interval_width_sum += orbit.interval_width_sum;
        maximum_interval_width = maximum_interval_width.max(orbit.maximum_interval_width);
        wrap_splits += orbit.wrap_splits;
        prefix_steps += orbit.prefix_steps;
        signed_sum_terms += orbit.signed_sum_terms;
        generation_wall_ns += orbit.generation_wall_ns;
        charged += charged_additions(orbit.capacity, orbit.sign_bits);
    }
    let source_intervals = intervals.len() as u64;
    let union = union_intervals(intervals)?;
    let distinct_fraction = union.width as f64 / states_emitted as f64;
    let direct_distinct_states = direct.map(|support| support.len() as u64);
    let (false_positives, false_negatives) = if let Some(support) = direct {
        let intersection = support
            .iter()
            .filter(|value| interval_contains(&union.merged, value))
            .count() as u64;
        (
            Some(union.width - intersection),
            Some(support.len() as u64 - intersection),
        )
    } else {
        (None, None)
    };
    let boundary_matches = orbits
        .iter()
        .map(|orbit| {
            expected_boundary_digest(round291, orbit.attempt, orbit.sign_bits)
                .map(|expected| expected == orbit.boundary_digest)
        })
        .collect::<Vec<_>>();
    let round291_boundary_digests_match = if boundary_matches.iter().all(Option::is_some) {
        Some(boundary_matches.iter().all(|value| value == &Some(true)))
    } else {
        None
    };
    let uncorrected_actual_stage_ratio = LOCAL_RATIO * charged as f64 / states_emitted as f64;
    Ok(CoverageCell {
        class: class.into(),
        tuple_count: orbits.len() as u64,
        tuple_attempts: orbits.iter().map(|orbit| orbit.attempt).collect(),
        capacities: orbits.iter().map(|orbit| orbit.capacity).collect(),
        sign_bits: orbits[0].sign_bits,
        sign_segments: orbits.iter().map(|orbit| 1u64 << orbit.sign_bits).sum(),
        states_emitted,
        per_segment_distinct_sum: interval_width_sum,
        exact_union_distinct_states: union.width,
        within_segment_duplicates: states_emitted - interval_width_sum,
        cross_segment_duplicates: interval_width_sum - union.width,
        total_duplicates: states_emitted - union.width,
        distinct_fraction,
        maximum_interval_width,
        source_intervals,
        merged_intervals: union.merged.len() as u64,
        wrap_splits,
        charged_addition_equivalents: charged,
        uncorrected_actual_stage_ratio,
        corrected_actual_stage_ratio_per_distinct: LOCAL_RATIO * charged as f64
            / union.width as f64,
        corrected_frozen_round291_ratio_per_distinct: ROUND291_STAGE_RATIO / distinct_fraction,
        prefix_steps,
        signed_sum_terms,
        sort_comparisons: union.sort_comparisons,
        merge_comparisons: union.merge_comparisons,
        logical_interval_endpoint_bytes: source_intervals * 64,
        generation_wall_ns,
        sort_merge_wall_ns: union.wall_ns,
        peak_rss_bytes: peak_rss_bytes(),
        support_digest: union.digest,
        boundary_digests: orbits
            .iter()
            .map(|orbit| orbit.boundary_digest.clone())
            .collect(),
        round291_boundary_digests_match,
        direct_control_complete: direct.is_some(),
        direct_distinct_states,
        false_positives,
        false_negatives,
    })
}

fn toy_order(capacities: &[u64]) -> Vec<usize> {
    let mut remaining = capacities.to_vec();
    let mut order = Vec::new();
    loop {
        let mut changed = false;
        for (slot, left) in remaining.iter_mut().enumerate() {
            if *left != 0 {
                order.push(slot);
                *left -= 1;
                changed = true;
            }
        }
        if !changed {
            return order;
        }
    }
}

fn toy_cell(prime: u64, coefficients: &[u64], capacities: &[u64]) -> Result<ToyCell, String> {
    let modulus = BigUint::from(prime);
    let originals = coefficients
        .iter()
        .copied()
        .map(BigUint::from)
        .collect::<Vec<_>>();
    let endpoints = coefficients
        .iter()
        .zip(capacities)
        .map(|(coefficient, capacity)| BigUint::from((coefficient + capacity) % prime))
        .collect::<Vec<_>>();
    let tuple = PreparedTuple {
        attempt: 0,
        capacity: capacities.iter().sum(),
        order: toy_order(capacities),
        original_coefficients: originals,
        endpoint_coefficients: endpoints,
    };
    let sign_bits = coefficients.len() as u32;
    let orbit = generate_orbit(&tuple, sign_bits, &modulus)?;
    let union = union_intervals(orbit.intervals)?;
    let direct = direct_support(&tuple, sign_bits, &modulus);
    let intersection = direct
        .iter()
        .filter(|value| interval_contains(&union.merged, value))
        .count() as u64;
    let false_positives = union.width - intersection;
    let false_negatives = direct.len() as u64 - intersection;
    let direct_duplicates = orbit.emitted_states - direct.len() as u64;
    let interval_duplicates = orbit.emitted_states - union.width;
    Ok(ToyCell {
        prime,
        arity: coefficients.len() as u64,
        sign_segments: 1u64 << sign_bits,
        states_emitted: orbit.emitted_states,
        direct_distinct_states: direct.len() as u64,
        interval_distinct_states: union.width,
        direct_duplicate_states: direct_duplicates,
        interval_duplicate_states: interval_duplicates,
        false_positives,
        false_negatives,
        exact: false_positives == 0
            && false_negatives == 0
            && direct_duplicates == interval_duplicates,
    })
}

fn fit_exponent(xs: &[f64], ys: &[f64]) -> f64 {
    let mean_x = xs.iter().sum::<f64>() / xs.len() as f64;
    let mean_y = ys.iter().sum::<f64>() / ys.len() as f64;
    let numerator = xs
        .iter()
        .zip(ys)
        .map(|(x, y)| (x - mean_x) * (y - mean_y))
        .sum::<f64>();
    let denominator = xs.iter().map(|x| (x - mean_x).powi(2)).sum::<f64>();
    numerator / denominator
}

fn run(cli: Cli) -> Result<(), String> {
    let (dep291, round291) = load_json(
        &cli.round291,
        ROUND291_SHA256,
        "p256-serpentine-sign-orbit-round291-v1",
    )?;
    let (dep31, round31) = load_json(
        &cli.round31,
        ROUND31_SHA256,
        "p256.low_delta_factor_base/v1",
    )?;
    verify_round291(&round291)?;

    let toy_cells = vec![
        toy_cell(257, &[3, 29, 101], &[3, 4, 2])?,
        toy_cell(65_537, &[7, 101, 9_997, 31_337, 50_003], &[2, 5, 3, 4, 6])?,
    ];
    if toy_cells.iter().any(|cell| !cell.exact) {
        return Err("toy interval support control failed".into());
    }

    let curve = CurveParams::p256();
    let coefficients = native_coefficients(&round31, &curve.n)?;
    let capacities = residual_capacities(COLUMNS, RARE_EDGES);
    let tuples = parse_tuples(&round291, &coefficients, &capacities, &curve.n)?;

    let mut primary_cells = Vec::new();
    for depth in PRIMARY_DEPTHS {
        let orbit = generate_orbit(&tuples[0], depth, &curve.n)?;
        let direct = if depth <= DIRECT_CONTROL_MAX_DEPTH {
            Some(direct_support(&tuples[0], depth, &curve.n))
        } else {
            None
        };
        primary_cells.push(coverage_cell(
            "primary-depth",
            &[&orbit],
            &round291,
            direct.as_ref(),
        )?);
    }

    let full_orbits = tuples
        .iter()
        .map(|tuple| generate_orbit(tuple, FULL_SIGN_BITS, &curve.n))
        .collect::<Result<Vec<_>, _>>()?;
    let mut portfolio_cells = Vec::new();
    for tuple_count in 1..=full_orbits.len() {
        let selected = full_orbits[..tuple_count].iter().collect::<Vec<_>>();
        portfolio_cells.push(coverage_cell(
            "cumulative-full-depth-portfolio",
            &selected,
            &round291,
            None,
        )?);
    }

    let native_controls_exact = primary_cells.iter().all(|cell| {
        cell.false_positives.unwrap_or(0) == 0
            && cell.false_negatives.unwrap_or(0) == 0
            && cell
                .direct_distinct_states
                .is_none_or(|direct| direct == cell.exact_union_distinct_states)
    });
    let boundary_digests_match = primary_cells
        .iter()
        .all(|cell| cell.round291_boundary_digests_match == Some(true));
    if !native_controls_exact || !boundary_digests_match {
        return Err("native exact support or Round 291 boundary control failed".into());
    }

    let xs = primary_cells
        .iter()
        .map(|cell| f64::from(cell.sign_bits))
        .collect::<Vec<_>>();
    let emitted_ys = primary_cells
        .iter()
        .map(|cell| (cell.states_emitted as f64).log2())
        .collect::<Vec<_>>();
    let distinct_ys = primary_cells
        .iter()
        .map(|cell| (cell.exact_union_distinct_states as f64).log2())
        .collect::<Vec<_>>();
    let memory_ys = primary_cells
        .iter()
        .map(|cell| (cell.logical_interval_endpoint_bytes as f64).log2())
        .collect::<Vec<_>>();
    let fits = Fits {
        primary_depths: PRIMARY_DEPTHS.to_vec(),
        emitted_state_exponent: fit_exponent(&xs, &emitted_ys),
        distinct_state_exponent: fit_exponent(&xs, &distinct_ys),
        interval_memory_exponent: fit_exponent(&xs, &memory_ys),
        interpretation: "Fits are over sign depth for one frozen tuple. They describe this exact orbit only and are not extrapolated to global tuple coverage or relation probability.".into(),
    };

    let primary_full = primary_cells.last().expect("primary cells nonempty");
    let portfolio_full = portfolio_cells.last().expect("portfolio cells nonempty");
    let storage_below = primary_cells
        .iter()
        .chain(&portfolio_cells)
        .all(|cell| cell.logical_interval_endpoint_bytes < (1u64 << 50));
    let primary_fraction_gate = primary_full.distinct_fraction >= ROUND291_STAGE_RATIO;
    let primary_ratio_gate = primary_full.corrected_frozen_round291_ratio_per_distinct <= 1.0;
    let portfolio_ratio_gate = portfolio_full.corrected_frozen_round291_ratio_per_distinct <= 1.0;
    let coverage_hypothesis_passed = toy_cells.iter().all(|cell| cell.exact)
        && native_controls_exact
        && boundary_digests_match
        && primary_fraction_gate
        && primary_ratio_gate
        && portfolio_ratio_gate
        && storage_below;
    let gates = Gates {
        zero_false_positives_and_false_negatives: native_controls_exact,
        exact_toy_support_and_duplicates: toy_cells.iter().all(|cell| cell.exact),
        native_boundary_digests_match_round291: boundary_digests_match,
        full_depth_primary_distinct_fraction_at_least_round291_stage_ratio: primary_fraction_gate,
        full_depth_primary_corrected_stage_ratio_at_most_one: primary_ratio_gate,
        four_tuple_corrected_stage_ratio_at_most_one: portfolio_ratio_gate,
        materialized_storage_below_2_50: storage_below,
        coverage_hypothesis_passed,
        same_family_structured_degree_at_most_5: false,
        relation_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        complete_non_generic_dlp_below_rho: false,
        promoted: false,
    };
    let semantic = serde_json::to_vec(&serde_json::json!({
        "toy": &toy_cells,
        "primary": &primary_cells,
        "portfolio": &portfolio_cells,
        "fits": &fits,
        "gates": &gates,
    }))
    .map_err(|error| error.to_string())?;
    let decision = if coverage_hypothesis_passed {
        "Coverage gates passed for the frozen orbits. Retain as a stage candidate only and next integrate exact collision detection and relation predicates; do not claim end-to-end parity."
    } else {
        "Coverage gates failed. Round 291 remains a local emitted-state construction result; use the measured duplicate factor as the obstruction and do not claim rho parity or attempt an unplanted full-depth relation."
    };
    let result = ResultRecord {
        schema: "p256-serpentine-orbit-coverage-round292-v1".into(),
        curve: CURVE.into(),
        round: 292,
        family: "round-33/291 known-log two-delta cyclic coefficient base".into(),
        comparison_factor_base: COMPARISON_FB.into(),
        dependencies: vec![dep291, dep31],
        frozen: serde_json::json!({
            "columns": COLUMNS,
            "arity": ARITY,
            "rare_edges": RARE_EDGES,
            "cutoff": CUTOFF,
            "primary_depths": PRIMARY_DEPTHS,
            "portfolio_tuple_attempts": EXPECTED_ATTEMPTS,
            "full_sign_bits": FULL_SIGN_BITS,
            "round291_stage_ratio": ROUND291_STAGE_RATIO,
            "local_ratio": LOCAL_RATIO,
            "coefficient_sha256": COEFFICIENT_SHA256,
            "direct_control_max_depth": DIRECT_CONTROL_MAX_DEPTH,
        }),
        exact_interval_argument: "Each emitted coefficient path changes by exactly +1 or -1. Its visited support is therefore every integer between its minimum and maximum prefix sum. Translating that inclusive interval by the exact signed start coefficient modulo the prime group order, splitting only at zero, and exactly unioning all intervals gives precisely the covered P-256 generator multiples.".into(),
        toy_cells,
        primary_cells,
        portfolio_cells,
        fits,
        gates,
        semantic_evidence_sha256: hex::encode(sha256(&semantic)),
        relations_reported: 0,
        full_depth_unplanted_relation_attempted: false,
        interpretation: "This result measures exact support of selected correlated sign-orbit streams. It does not measure a summation-polynomial root, usable relation probability, structured degree, independent row collection, sparse linear algebra, or complete DLP recovery.".into(),
        decision: decision.into(),
    };
    let bytes = serde_json::to_vec_pretty(&result).map_err(|error| error.to_string())?;
    if let Some(path) = cli.out {
        fs::write(&path, &bytes).map_err(|error| format!("{}: {error}", path.display()))?;
    } else {
        println!(
            "{}",
            String::from_utf8(bytes).map_err(|error| error.to_string())?
        );
    }
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("P-256 orbit coverage audit failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modular_interval_splits_at_zero() {
        let modulus = BigUint::from(17u8);
        let (pieces, wrapped) = split_modular_interval(&BigUint::from(16u8), -1, 2, &modulus);
        assert!(wrapped);
        let union = union_intervals(pieces).unwrap();
        assert_eq!(union.width, 4);
        assert!(interval_contains(&union.merged, &BigUint::from(0u8)));
        assert!(interval_contains(&union.merged, &BigUint::from(15u8)));
        assert!(!interval_contains(&union.merged, &BigUint::from(2u8)));
    }

    #[test]
    fn adjacent_intervals_merge_inclusively() {
        let intervals = vec![
            Interval {
                start: BigUint::from(8u8),
                end: BigUint::from(10u8),
            },
            Interval {
                start: BigUint::from(4u8),
                end: BigUint::from(7u8),
            },
        ];
        let union = union_intervals(intervals).unwrap();
        assert_eq!(union.width, 7);
        assert_eq!(union.merged.len(), 1);
    }

    #[test]
    fn toy_support_is_exact() {
        let first = toy_cell(257, &[3, 29, 101], &[3, 4, 2]).unwrap();
        let second = toy_cell(65_537, &[7, 101, 9_997, 31_337, 50_003], &[2, 5, 3, 4, 6]).unwrap();
        assert!(first.exact);
        assert!(second.exact);
        assert_eq!(first.false_positives + first.false_negatives, 0);
        assert_eq!(second.false_positives + second.false_negatives, 0);
    }
}
