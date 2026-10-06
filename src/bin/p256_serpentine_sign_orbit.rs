//! Table-free serpentine sign-orbit streaming for P-256 (round 291).

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use blake3::Hasher;
use clap::Parser;
use crypto_lib::ct_bignum::U256;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::Zero;
use serde::Serialize;
use serde_json::Value;

const CURVE: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const COMPARISON_FB: &str = "FB1h2f8621cda105";
const ROUND33_SHA256: &str = "9931bccbd9e2821f4498f465ce65f83efb400e7fb3740a239bdb40d3275ffbd8";
const ROUND31_SHA256: &str = "9b9ee16c5a24e868d1b9304f08aa80ca39ac23b1b4586ef730534e739172c435";
const ROUND290_SHA256: &str = "41cdec6ad5265464cc512b2f026f6fa77b333aa209f163e7e77ab78ba665b876";
const COEFFICIENT_SHA256: &str = "980917981827d813e60484abb0655e8bd527b0beb540ea146d2974ff17303a53";
const COLUMNS: u64 = 131_458;
const RARE_EDGES: u64 = 6_935;
const ARITY: usize = 17;
const CUTOFF: u64 = 219;
const FULL_SIGN_BITS: u32 = 17;
const LOCAL_RATIO: f64 = 0.964_336_477_130_181;
const MEAN_CAPACITY: f64 = 224.361_249_307_510_55;
const SIGN_DEPTHS: [u32; 6] = [8, 10, 12, 14, 16, 17];
const HOLDOUT_DEPTHS: [u32; 3] = [8, 10, 12];
const TIMING_DEPTH: u32 = 12;
const TIMING_REPETITIONS: u64 = 7;
const ADDITION_CONTROL_COUNT: u64 = 262_144;

#[derive(Parser)]
#[command(about = "Audit P-256 serpentine sign-orbit streaming (round 291)")]
struct Cli {
    #[arg(long)]
    round33: PathBuf,
    #[arg(long)]
    round31: PathBuf,
    #[arg(long)]
    round290: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Clone, Serialize)]
struct Dependency {
    path: String,
    sha256: String,
    schema: String,
}

#[derive(Clone, Serialize)]
struct ToyCell {
    prime: u64,
    columns: u64,
    rare_edges: u64,
    arity: u64,
    unsigned_tuples: u64,
    signed_segments: u64,
    states_per_reference: u64,
    direct_distinct_states: u64,
    serpentine_distinct_states: u64,
    multiplicity_mismatches: u64,
    false_positives: u64,
    false_negatives: u64,
    transition_failures: u64,
    endpoint_failures: u64,
    direct_digest: String,
    serpentine_digest: String,
}

#[derive(Clone, Serialize)]
struct NativeCell {
    class: String,
    tuple_attempt: u64,
    columns: Vec<u64>,
    capacity: u64,
    sign_bits: u32,
    sign_segments: u64,
    states_emitted: u64,
    path_additions: u64,
    sign_boundary_additions: u64,
    initial_sum_additions: u64,
    endpoint_doublings: u64,
    target_correction_additions: u64,
    charged_addition_equivalents: u64,
    additions_per_state: f64,
    stage_ratio_to_rho: f64,
    direct_boundary_checks: u64,
    transition_failures: u64,
    endpoint_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    wall_ns_with_validation: u64,
    coefficient_checksum: u64,
    boundary_digest: String,
    complete_checked: bool,
}

#[derive(Clone, Serialize)]
struct RawTiming {
    repetition: u64,
    order: Vec<String>,
    direct_wall_ns: u64,
    serpentine_wall_ns: u64,
    direct_checksum: u64,
    serpentine_checksum: u64,
}

#[derive(Clone, Serialize)]
struct TimingSummary {
    depth_bits: u32,
    repetitions: u64,
    aa_first_ns: u64,
    aa_second_ns: u64,
    aa_ratio: f64,
    native_addition_median_ns: f64,
    direct_median_ns: f64,
    serpentine_median_ns: f64,
    direct_minimum_ns: u64,
    serpentine_minimum_ns: u64,
    measured_speedup_direct_over_serpentine: f64,
    checksum_consistent_per_variant: bool,
    includes_target_correction: bool,
    validation_excluded: bool,
    raw: Vec<RawTiming>,
}

#[derive(Clone, Serialize)]
struct Accounting {
    mean_capacity: f64,
    signs_per_tuple: u64,
    asymptotic_additions_per_state: f64,
    conservative_full_orbit_additions_per_state: f64,
    conservative_stage_ratio_to_rho: f64,
    pair_table_entries: u64,
    pair_table_bytes_removed: u64,
    radix_records_removed_per_start: u64,
    peak_candidate_bytes: u64,
    peak_below_2_50: bool,
}

#[derive(Clone, Serialize)]
struct Fits {
    counted_operation_exponent: f64,
    emitted_state_exponent: f64,
    peak_memory_exponent: f64,
    measured_time_exponent: f64,
    depths: Vec<u32>,
    interpretation: String,
}

#[derive(Clone, Serialize)]
struct Gates {
    zero_false_positives_and_false_negatives: bool,
    exact_direct_and_native_replay: bool,
    counted_stage_below_rho: bool,
    measured_complete_selector_below_rho: bool,
    proved_p256_coverage_or_relation_probability: bool,
    structured_residual_degree_at_most_5_same_family: bool,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    storage_below_2_50: bool,
    demonstrably_non_generic_end_to_end: bool,
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
    toy_cells: Vec<ToyCell>,
    native_cells: Vec<NativeCell>,
    timing: TimingSummary,
    accounting: Accounting,
    fits: Fits,
    gates: Gates,
    semantic_evidence_sha256: String,
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

fn enumerate_combinations(
    columns: usize,
    arity: usize,
    start: usize,
    selected: &mut Vec<usize>,
    visit: &mut impl FnMut(&[usize]),
) {
    if selected.len() == arity {
        visit(selected);
        return;
    }
    let needed = arity - selected.len();
    for index in start..=columns - needed {
        selected.push(index);
        enumerate_combinations(columns, arity, index + 1, selected, visit);
        selected.pop();
    }
}

fn toy_coefficients(prime: u64, columns: u64, rare: u64) -> Result<Vec<u64>, String> {
    fn pow_mod(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
        let mut result = 1u64;
        while exponent != 0 {
            if exponent & 1 == 1 {
                result = (u128::from(result) * u128::from(base) % u128::from(modulus)) as u64;
            }
            base = (u128::from(base) * u128::from(base) % u128::from(modulus)) as u64;
            exponent >>= 1;
        }
        result
    }
    let common = columns - rare;
    let inverse = pow_mod(rare, prime - 2, prime);
    let rare_delta =
        (prime - (u128::from(common) * u128::from(inverse) % u128::from(prime)) as u64) % prime;
    let mut relative = Vec::with_capacity(columns as usize);
    let mut current = 0u64;
    for index in 0..columns {
        relative.push(current);
        current = (current
            + if mechanical_rare(index, columns, rare) {
                rare_delta
            } else {
                1
            })
            % prime;
    }
    if current != 0 {
        return Err("toy coefficient cycle did not close".into());
    }
    for anchor in 1..prime {
        let coefficients = relative
            .iter()
            .map(|value| (value + anchor) % prime)
            .collect::<Vec<_>>();
        if coefficients.contains(&0) {
            continue;
        }
        let signed = coefficients
            .iter()
            .map(|value| (*value).min(prime - value))
            .collect::<BTreeSet<_>>();
        if signed.len() == columns as usize {
            return Ok(coefficients);
        }
    }
    Err("toy coefficient cycle has no signed-unique anchor".into())
}

fn toy_path(
    coefficients: &[u64],
    columns: &[usize],
    capacities: &[u64],
    signs: u64,
    forward: bool,
    prime: u64,
) -> Result<Vec<u64>, String> {
    let cols = columns
        .iter()
        .map(|column| *column as u64)
        .collect::<Vec<_>>();
    let order = forward_order(&cols, capacities)?;
    let positions = if forward {
        cols.clone()
    } else {
        cols.iter()
            .map(|column| (column + capacities[*column as usize]) % capacities.len() as u64)
            .collect::<Vec<_>>()
    };
    let mut sum = 0u64;
    for (slot, position) in positions.iter().enumerate() {
        let coefficient = coefficients[*position as usize];
        sum = if signs & (1 << slot) == 0 {
            (sum + coefficient) % prime
        } else {
            (sum + prime - coefficient) % prime
        };
    }
    let mut states = vec![sum];
    let steps: Box<dyn Iterator<Item = &usize>> = if forward {
        Box::new(order.iter())
    } else {
        Box::new(order.iter().rev())
    };
    for slot in steps {
        let positive = signs & (1 << slot) == 0;
        let add_one = positive == forward;
        sum = if add_one {
            (sum + 1) % prime
        } else {
            (sum + prime - 1) % prime
        };
        states.push(sum);
    }
    Ok(states)
}

fn multiset_digest(states: &BTreeMap<u64, u64>) -> String {
    let mut hasher = Hasher::new();
    for (state, multiplicity) in states {
        hasher.update(&state.to_be_bytes());
        hasher.update(&multiplicity.to_be_bytes());
    }
    hasher.finalize().to_hex().to_string()
}

fn toy_cell(prime: u64, columns: u64, rare: u64, arity: usize) -> Result<ToyCell, String> {
    let coefficients = toy_coefficients(prime, columns, rare)?;
    let capacities = residual_capacities(columns, rare);
    let mut direct = BTreeMap::<u64, u64>::new();
    let mut serpentine = BTreeMap::<u64, u64>::new();
    let mut unsigned = 0u64;
    let mut segments = 0u64;
    let mut states = 0u64;
    let mut transition_failures = 0u64;
    let mut endpoint_failures = 0u64;
    enumerate_combinations(
        columns as usize,
        arity,
        0,
        &mut Vec::new(),
        &mut |selected| {
            unsigned += 1;
            for mask in 0..1u64 << arity {
                let path = toy_path(&coefficients, selected, &capacities, mask, true, prime)
                    .expect("valid direct toy path");
                segments += 1;
                states += path.len() as u64;
                for state in path {
                    *direct.entry(state).or_default() += 1;
                }
            }
            let mut previous_end: Option<(u64, u64)> = None;
            for ordinal in 0..1u64 << arity {
                let gray = ordinal ^ (ordinal >> 1);
                let forward = ordinal & 1 == 0;
                let path = toy_path(&coefficients, selected, &capacities, gray, forward, prime)
                    .expect("valid serpentine toy path");
                if let Some((old_gray, old_end)) = previous_end {
                    let changed = (old_gray ^ gray).trailing_zeros() as usize;
                    let position = if forward {
                        selected[changed] as u64
                    } else {
                        (selected[changed] as u64 + capacities[selected[changed]]) % columns
                    };
                    let coefficient = coefficients[position as usize];
                    let expected = if old_gray & (1 << changed) == 0 {
                        (old_end + prime - (2 * coefficient) % prime) % prime
                    } else {
                        (old_end + (2 * coefficient) % prime) % prime
                    };
                    if path.first().copied() != Some(expected) {
                        endpoint_failures += 1;
                    }
                }
                if path.len()
                    != capacities
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| selected.contains(i))
                        .map(|(_, c)| *c)
                        .sum::<u64>() as usize
                        + 1
                {
                    transition_failures += 1;
                }
                previous_end = path.last().copied().map(|end| (gray, end));
                for state in path {
                    *serpentine.entry(state).or_default() += 1;
                }
            }
        },
    );
    let false_positives = serpentine
        .keys()
        .filter(|state| !direct.contains_key(state))
        .count() as u64;
    let false_negatives = direct
        .keys()
        .filter(|state| !serpentine.contains_key(state))
        .count() as u64;
    let multiplicity_mismatches = direct
        .iter()
        .filter(|(state, count)| serpentine.get(state) != Some(count))
        .count() as u64
        + serpentine
            .keys()
            .filter(|state| !direct.contains_key(state))
            .count() as u64;
    Ok(ToyCell {
        prime,
        columns,
        rare_edges: rare,
        arity: arity as u64,
        unsigned_tuples: unsigned,
        signed_segments: segments,
        states_per_reference: states,
        direct_distinct_states: direct.len() as u64,
        serpentine_distinct_states: serpentine.len() as u64,
        multiplicity_mismatches,
        false_positives,
        false_negatives,
        transition_failures,
        endpoint_failures,
        direct_digest: multiset_digest(&direct),
        serpentine_digest: multiset_digest(&serpentine),
    })
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
    if common != BigUint::from(1u8) {
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

fn accepted_tuples(capacities: &[u64], count: usize) -> Vec<(u64, Vec<u64>)> {
    let mut accepted = Vec::new();
    let mut attempt = 0u64;
    while accepted.len() < count {
        let columns = draw_columns(attempt);
        let capacity = columns
            .iter()
            .map(|column| capacities[*column as usize])
            .sum::<u64>();
        if capacity >= CUTOFF {
            accepted.push((attempt, columns));
        }
        attempt += 1;
    }
    accepted
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

fn equivalent(left: &Projective, right: &Projective) -> bool {
    if bool::from(left.is_identity()) || bool::from(right.is_identity()) {
        return bool::from(left.is_identity()) && bool::from(right.is_identity());
    }
    left.x.mul(&right.z) == right.x.mul(&left.z) && left.y.mul(&right.z) == right.y.mul(&left.z)
}

struct PreparedTuple {
    attempt: u64,
    columns: Vec<u64>,
    capacity: u64,
    order: Vec<usize>,
    original_points: Vec<Projective>,
    endpoint_points: Vec<Projective>,
    original_coefficients: Vec<BigUint>,
    endpoint_coefficients: Vec<BigUint>,
}

fn prepare_tuple(
    attempt: u64,
    columns: Vec<u64>,
    coefficients: &[BigUint],
    capacities: &[u64],
    generator: &Projective,
    modulus: &BigUint,
) -> Result<PreparedTuple, String> {
    let order = forward_order(&columns, capacities)?;
    let capacity = order.len() as u64;
    let mut original_points = Vec::with_capacity(ARITY);
    let mut endpoint_points = Vec::with_capacity(ARITY);
    let mut original_coefficients = Vec::with_capacity(ARITY);
    let mut endpoint_coefficients = Vec::with_capacity(ARITY);
    for column in &columns {
        let original_coefficient = coefficients[*column as usize].clone();
        let endpoint_column = (column + capacities[*column as usize]) % COLUMNS;
        let endpoint_coefficient = coefficients[endpoint_column as usize].clone();
        let original_point = generator.scalar_mul_ct(&original_coefficient, 256);
        let endpoint_point = generator.scalar_mul_ct(&endpoint_coefficient, 256);
        original_points.push(original_point);
        endpoint_points.push(endpoint_point);
        original_coefficients.push(original_coefficient);
        endpoint_coefficients.push(endpoint_coefficient);
    }
    for slot in 0..ARITY {
        let mut expected = original_coefficients[slot].clone();
        add_mod(
            &mut expected,
            &BigUint::from(capacities[columns[slot] as usize]),
            true,
            modulus,
        );
        if expected != endpoint_coefficients[slot] {
            return Err("prepared endpoint coefficient mismatch".into());
        }
    }
    Ok(PreparedTuple {
        attempt,
        columns,
        capacity,
        order,
        original_points,
        endpoint_points,
        original_coefficients,
        endpoint_coefficients,
    })
}

fn direct_state(
    tuple: &PreparedTuple,
    gray: u64,
    at_original: bool,
    modulus: &BigUint,
) -> (Projective, BigUint) {
    let points = if at_original {
        &tuple.original_points
    } else {
        &tuple.endpoint_points
    };
    let coefficients = if at_original {
        &tuple.original_coefficients
    } else {
        &tuple.endpoint_coefficients
    };
    let mut point = Projective::IDENTITY;
    let mut coefficient = BigUint::zero();
    for slot in 0..ARITY {
        let positive = gray & (1 << slot) == 0;
        point = if positive {
            point.add(&points[slot])
        } else {
            let negative = points[slot].neg();
            point.add(&negative)
        };
        add_mod(&mut coefficient, &coefficients[slot], positive, modulus);
    }
    (point, coefficient)
}

fn low_word(value: &BigUint) -> u64 {
    value.to_u64_digits().first().copied().unwrap_or(0)
}

fn native_cell(
    class: &str,
    tuple: &PreparedTuple,
    sign_bits: u32,
    generator: &Projective,
    modulus: &BigUint,
) -> NativeCell {
    let signs = 1u64 << sign_bits;
    let negative_generator = generator.neg();
    let (mut point, mut coefficient) = direct_state(tuple, 0, true, modulus);
    let mut transition_failures = 0u64;
    let mut endpoint_failures = 0u64;
    let false_positives = 0u64;
    let mut false_negatives = 0u64;
    let mut checks = 0u64;
    let mut checksum = 0u64;
    let mut digest = Hasher::new();
    let begin = Instant::now();
    for ordinal in 0..signs {
        let gray = ordinal ^ (ordinal >> 1);
        let forward = ordinal & 1 == 0;
        let (reference_point, reference_coefficient) = direct_state(tuple, gray, forward, modulus);
        checks += 1;
        if !equivalent(&point, &reference_point) || coefficient != reference_coefficient {
            endpoint_failures += 1;
        }
        digest.update(&ordinal.to_be_bytes());
        digest.update(&gray.to_be_bytes());
        digest.update(&fixed_32(&coefficient).expect("coefficient fits"));
        let steps: Box<dyn Iterator<Item = &usize>> = if forward {
            Box::new(tuple.order.iter())
        } else {
            Box::new(tuple.order.iter().rev())
        };
        for slot in steps {
            let positive_sign = gray & (1 << slot) == 0;
            let add_positive = positive_sign == forward;
            point = point.add(if add_positive {
                generator
            } else {
                &negative_generator
            });
            add_mod(&mut coefficient, &BigUint::from(1u8), add_positive, modulus);
            checksum = checksum.rotate_left(7) ^ low_word(&coefficient);
        }
        let (end_point, end_coefficient) = direct_state(tuple, gray, !forward, modulus);
        checks += 1;
        if !equivalent(&point, &end_point) || coefficient != end_coefficient {
            transition_failures += 1;
        }
        let mut mutated = coefficient.clone();
        add_mod(&mut mutated, &BigUint::from(1u8), true, modulus);
        if mutated == end_coefficient {
            false_negatives += 1;
        }
        if ordinal + 1 != signs {
            let next_gray = (ordinal + 1) ^ ((ordinal + 1) >> 1);
            let changed = (gray ^ next_gray).trailing_zeros() as usize;
            let current_points = if forward {
                &tuple.endpoint_points
            } else {
                &tuple.original_points
            };
            let current_coefficients = if forward {
                &tuple.endpoint_coefficients
            } else {
                &tuple.original_coefficients
            };
            let old_positive = gray & (1 << changed) == 0;
            let doubled = current_points[changed].double();
            point = if old_positive {
                let negative = doubled.neg();
                point.add(&negative)
            } else {
                point.add(&doubled)
            };
            let twice = (&current_coefficients[changed] << 1usize) % modulus;
            add_mod(&mut coefficient, &twice, !old_positive, modulus);
        }
    }
    let wall_ns = begin.elapsed().as_nanos() as u64;
    let path_additions = signs * tuple.capacity;
    let boundary_additions = signs - 1;
    let initial_additions = (ARITY - 1) as u64;
    let endpoint_doublings = ARITY as u64;
    let target_additions = signs;
    let charged = path_additions
        + boundary_additions
        + initial_additions
        + 2 * endpoint_doublings
        + target_additions;
    let states = signs * (tuple.capacity + 1);
    let additions_per_state = charged as f64 / states as f64;
    NativeCell {
        class: class.into(),
        tuple_attempt: tuple.attempt,
        columns: tuple.columns.clone(),
        capacity: tuple.capacity,
        sign_bits,
        sign_segments: signs,
        states_emitted: states,
        path_additions,
        sign_boundary_additions: boundary_additions,
        initial_sum_additions: initial_additions,
        endpoint_doublings,
        target_correction_additions: target_additions,
        charged_addition_equivalents: charged,
        additions_per_state,
        stage_ratio_to_rho: LOCAL_RATIO * additions_per_state,
        direct_boundary_checks: checks,
        transition_failures,
        endpoint_failures,
        false_positives,
        false_negatives,
        wall_ns_with_validation: wall_ns,
        coefficient_checksum: checksum,
        boundary_digest: digest.finalize().to_hex().to_string(),
        complete_checked: true,
    }
}

fn kernel_direct(
    tuple: &PreparedTuple,
    sign_bits: u32,
    generator: &Projective,
    modulus: &BigUint,
) -> u64 {
    let signs = 1u64 << sign_bits;
    let negative_generator = generator.neg();
    let mut checksum = 0u64;
    for ordinal in 0..signs {
        let gray = ordinal ^ (ordinal >> 1);
        let (mut point, _) = direct_state(tuple, gray, true, modulus);
        for slot in &tuple.order {
            point = point.add(if gray & (1 << slot) == 0 {
                generator
            } else {
                &negative_generator
            });
        }
        let corrected = black_box(point.add(generator));
        checksum ^= corrected.x.to_canonical().0[0].rotate_left((ordinal & 63) as u32);
    }
    black_box(checksum)
}

fn kernel_serpentine(
    tuple: &PreparedTuple,
    sign_bits: u32,
    generator: &Projective,
    modulus: &BigUint,
) -> u64 {
    let signs = 1u64 << sign_bits;
    let negative_generator = generator.neg();
    let (mut point, _) = direct_state(tuple, 0, true, modulus);
    let mut checksum = 0u64;
    for ordinal in 0..signs {
        let gray = ordinal ^ (ordinal >> 1);
        let forward = ordinal & 1 == 0;
        let steps: Box<dyn Iterator<Item = &usize>> = if forward {
            Box::new(tuple.order.iter())
        } else {
            Box::new(tuple.order.iter().rev())
        };
        for slot in steps {
            let positive_sign = gray & (1 << slot) == 0;
            point = point.add(if positive_sign == forward {
                generator
            } else {
                &negative_generator
            });
        }
        let corrected = black_box(point.add(generator));
        checksum ^= corrected.x.to_canonical().0[0].rotate_left((ordinal & 63) as u32);
        if ordinal + 1 != signs {
            let next_gray = (ordinal + 1) ^ ((ordinal + 1) >> 1);
            let changed = (gray ^ next_gray).trailing_zeros() as usize;
            let current_points = if forward {
                &tuple.endpoint_points
            } else {
                &tuple.original_points
            };
            let old_positive = gray & (1 << changed) == 0;
            let doubled = current_points[changed].double();
            point = if old_positive {
                let negative = doubled.neg();
                point.add(&negative)
            } else {
                point.add(&doubled)
            };
        }
    }
    black_box(checksum)
}

fn timed_additions(count: u64, generator: &Projective) -> (u64, u64) {
    let mut point = *generator;
    let begin = Instant::now();
    for _ in 0..count {
        point = black_box(point.add(black_box(generator)));
    }
    (
        begin.elapsed().as_nanos() as u64,
        point.x.to_canonical().0[0],
    )
}

fn median(values: &[u64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    sorted[sorted.len() / 2] as f64
}

fn timing(tuple: &PreparedTuple, generator: &Projective, modulus: &BigUint) -> TimingSummary {
    let _ = black_box(kernel_direct(tuple, 8, generator, modulus));
    let _ = black_box(kernel_serpentine(tuple, 8, generator, modulus));
    let (aa_first_ns, aa_first_checksum) = timed_additions(ADDITION_CONTROL_COUNT, generator);
    let (aa_second_ns, aa_second_checksum) = timed_additions(ADDITION_CONTROL_COUNT, generator);
    assert_eq!(aa_first_checksum, aa_second_checksum);
    let mut raw = Vec::new();
    for repetition in 0..TIMING_REPETITIONS {
        let direct_first = repetition & 1 == 0;
        let (direct_wall_ns, direct_checksum, serpentine_wall_ns, serpentine_checksum) =
            if direct_first {
                let begin = Instant::now();
                let direct_checksum = kernel_direct(tuple, TIMING_DEPTH, generator, modulus);
                let direct_ns = begin.elapsed().as_nanos() as u64;
                let begin = Instant::now();
                let serpentine_checksum =
                    kernel_serpentine(tuple, TIMING_DEPTH, generator, modulus);
                let serpentine_ns = begin.elapsed().as_nanos() as u64;
                (
                    direct_ns,
                    direct_checksum,
                    serpentine_ns,
                    serpentine_checksum,
                )
            } else {
                let begin = Instant::now();
                let serpentine_checksum =
                    kernel_serpentine(tuple, TIMING_DEPTH, generator, modulus);
                let serpentine_ns = begin.elapsed().as_nanos() as u64;
                let begin = Instant::now();
                let direct_checksum = kernel_direct(tuple, TIMING_DEPTH, generator, modulus);
                let direct_ns = begin.elapsed().as_nanos() as u64;
                (
                    direct_ns,
                    direct_checksum,
                    serpentine_ns,
                    serpentine_checksum,
                )
            };
        raw.push(RawTiming {
            repetition,
            order: if direct_first {
                vec!["direct".into(), "serpentine".into()]
            } else {
                vec!["serpentine".into(), "direct".into()]
            },
            direct_wall_ns,
            serpentine_wall_ns,
            direct_checksum,
            serpentine_checksum,
        });
    }
    let direct = raw.iter().map(|row| row.direct_wall_ns).collect::<Vec<_>>();
    let serpentine = raw
        .iter()
        .map(|row| row.serpentine_wall_ns)
        .collect::<Vec<_>>();
    let direct_median = median(&direct);
    let serpentine_median = median(&serpentine);
    let direct_checksum = raw[0].direct_checksum;
    let serpentine_checksum = raw[0].serpentine_checksum;
    TimingSummary {
        depth_bits: TIMING_DEPTH,
        repetitions: TIMING_REPETITIONS,
        aa_first_ns,
        aa_second_ns,
        aa_ratio: aa_first_ns.max(aa_second_ns) as f64 / aa_first_ns.min(aa_second_ns) as f64,
        native_addition_median_ns: median(&[aa_first_ns, aa_second_ns])
            / ADDITION_CONTROL_COUNT as f64,
        direct_median_ns: direct_median,
        serpentine_median_ns: serpentine_median,
        direct_minimum_ns: *direct.iter().min().unwrap(),
        serpentine_minimum_ns: *serpentine.iter().min().unwrap(),
        measured_speedup_direct_over_serpentine: direct_median / serpentine_median,
        checksum_consistent_per_variant: raw.iter().all(|row| {
            row.direct_checksum == direct_checksum && row.serpentine_checksum == serpentine_checksum
        }),
        includes_target_correction: true,
        validation_excluded: true,
        raw,
    }
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

fn verify_dependencies(round33: &Value, round31: &Value, round290: &Value) -> Result<(), String> {
    if json_u64(round33, &["cutoff_sweep", "selected", "cutoff_capacity"])? != CUTOFF
        || json_u64(round33, &["relation_model", "arity"])? != ARITY as u64
        || json_str(round33, &["native_replay", "coefficient_sha256"])? != COEFFICIENT_SHA256
        || json_bool(round33, &["promotion_gates", "promoted"])?
    {
        return Err("round-33 frozen fields mismatch".into());
    }
    let mean = json_at(round33, &["cutoff_sweep", "selected", "mean_capacity"])?
        .as_f64()
        .ok_or("round-33 mean capacity is not f64")?;
    if mean.to_bits() != MEAN_CAPACITY.to_bits()
        || json_str(round33, &["dependency", "sha256"])? != ROUND31_SHA256
        || json_str(round31, &["family"])? != "known-log two-delta cyclic coefficient bases"
        || json_u64(round290, &["best_measured_kernel_round"])? != 108
        || json_u64(round290, &["best_formula_traffic_round"])? != 290
        || !json_at(round290, &["gates", "promoted_rounds"])?
            .as_array()
            .ok_or("round-290 promoted_rounds is not array")?
            .is_empty()
    {
        return Err("dependency boundary mismatch".into());
    }
    Ok(())
}

fn run(cli: Cli) -> Result<(), String> {
    let (dep33, round33) = load_json(
        &cli.round33,
        ROUND33_SHA256,
        "p256.biased_restart_selector/v1",
    )?;
    let (dep31, round31) = load_json(
        &cli.round31,
        ROUND31_SHA256,
        "p256.low_delta_factor_base/v1",
    )?;
    let (dep290, round290) = load_json(
        &cli.round290,
        ROUND290_SHA256,
        "p256-radix-portfolio-round43-290-v1",
    )?;
    verify_dependencies(&round33, &round31, &round290)?;

    let toy_cells = vec![toy_cell(257, 10, 3, 3)?, toy_cell(65_537, 20, 3, 4)?];
    if toy_cells.iter().any(|cell| {
        cell.multiplicity_mismatches != 0
            || cell.false_positives != 0
            || cell.false_negatives != 0
            || cell.transition_failures != 0
            || cell.endpoint_failures != 0
            || cell.direct_digest != cell.serpentine_digest
    }) {
        return Err("toy serpentine reference failed".into());
    }

    let curve = CurveParams::p256();
    let generator = Projective::from_affine(
        &U256::from_biguint(&curve.gx),
        &U256::from_biguint(&curve.gy),
    );
    let coefficients = native_coefficients(&round31, &curve.n)?;
    let capacities = residual_capacities(COLUMNS, RARE_EDGES);
    let accepted = accepted_tuples(&capacities, 4);
    let tuples = accepted
        .into_iter()
        .map(|(attempt, columns)| {
            prepare_tuple(
                attempt,
                columns,
                &coefficients,
                &capacities,
                &generator,
                &curve.n,
            )
        })
        .collect::<Result<Vec<_>, String>>()?;

    let mut native_cells = Vec::new();
    for depth in SIGN_DEPTHS {
        native_cells.push(native_cell(
            "primary", &tuples[0], depth, &generator, &curve.n,
        ));
    }
    for (tuple, depth) in tuples[1..].iter().zip(HOLDOUT_DEPTHS.iter()) {
        native_cells.push(native_cell("holdout", tuple, *depth, &generator, &curve.n));
    }
    let failures = native_cells.iter().any(|cell| {
        cell.transition_failures != 0
            || cell.endpoint_failures != 0
            || cell.false_positives != 0
            || cell.false_negatives != 0
    });
    if failures {
        return Err("native serpentine replay failed".into());
    }

    let timing = timing(&tuples[0], &generator, &curve.n);
    let signs = 1u64 << FULL_SIGN_BITS;
    let full_charged = signs as f64 * MEAN_CAPACITY
        + (signs - 1) as f64
        + signs as f64
        + (ARITY - 1) as f64
        + 2.0 * ARITY as f64;
    let full_states = signs as f64 * (MEAN_CAPACITY + 1.0);
    let additions_per_state = full_charged / full_states;
    let peak_bytes =
        (2 * ARITY * std::mem::size_of::<Projective>() + 4 * ARITY * 32 + 2 * ARITY * 8 + 4096)
            as u64;
    let accounting = Accounting {
        mean_capacity: MEAN_CAPACITY,
        signs_per_tuple: signs,
        asymptotic_additions_per_state: (MEAN_CAPACITY + 2.0) / (MEAN_CAPACITY + 1.0),
        conservative_full_orbit_additions_per_state: additions_per_state,
        conservative_stage_ratio_to_rho: LOCAL_RATIO * additions_per_state,
        pair_table_entries: 34_562_148_612,
        pair_table_bytes_removed: 1_140_550_904_196,
        radix_records_removed_per_start: 8,
        peak_candidate_bytes: peak_bytes,
        peak_below_2_50: peak_bytes < (1u64 << 50),
    };
    let primary = native_cells
        .iter()
        .filter(|cell| cell.class == "primary")
        .collect::<Vec<_>>();
    let xs = primary
        .iter()
        .map(|cell| (cell.sign_segments as f64).log2())
        .collect::<Vec<_>>();
    let operation_ys = primary
        .iter()
        .map(|cell| (cell.charged_addition_equivalents as f64).log2())
        .collect::<Vec<_>>();
    let state_ys = primary
        .iter()
        .map(|cell| (cell.states_emitted as f64).log2())
        .collect::<Vec<_>>();
    let time_ys = primary
        .iter()
        .map(|cell| (cell.wall_ns_with_validation as f64).log2())
        .collect::<Vec<_>>();
    let fits = Fits {
        counted_operation_exponent: fit_exponent(&xs, &operation_ys),
        emitted_state_exponent: fit_exponent(&xs, &state_ys),
        peak_memory_exponent: 0.0,
        measured_time_exponent: fit_exponent(&xs, &time_ys),
        depths: SIGN_DEPTHS.to_vec(),
        interpretation: "Fits are over sign-orbit depth for one tuple. Validation is included in cell time; the matched kernel timing is reported separately. No fit is transferred to tuple coverage or a complete DLP.".into(),
    };
    let exact = !failures
        && toy_cells.iter().all(|cell| {
            cell.multiplicity_mismatches == 0
                && cell.false_positives == 0
                && cell.false_negatives == 0
        });
    let gates = Gates {
        zero_false_positives_and_false_negatives: exact,
        exact_direct_and_native_replay: exact,
        counted_stage_below_rho: accounting.conservative_stage_ratio_to_rho <= 1.0,
        measured_complete_selector_below_rho: false,
        proved_p256_coverage_or_relation_probability: false,
        structured_residual_degree_at_most_5_same_family: false,
        relation_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        storage_below_2_50: accounting.peak_below_2_50,
        demonstrably_non_generic_end_to_end: false,
        promoted: false,
    };
    let semantic = serde_json::to_vec(&serde_json::json!({
        "toy": &toy_cells,
        "native": &native_cells,
        "timing": &timing,
        "accounting": &accounting,
        "fits": &fits,
        "gates": &gates,
    }))
    .map_err(|error| error.to_string())?;
    let result = ResultRecord {
        schema: "p256-serpentine-sign-orbit-round291-v1".into(),
        curve: CURVE.into(),
        round: 291,
        family: "round-33 known-log low-delta/scalar selector".into(),
        comparison_factor_base: COMPARISON_FB.into(),
        dependencies: vec![dep33, dep31, dep290],
        frozen: serde_json::json!({
            "columns": COLUMNS,
            "arity": ARITY,
            "rare_edges": RARE_EDGES,
            "cutoff": CUTOFF,
            "sign_depths": SIGN_DEPTHS,
            "holdout_depths": HOLDOUT_DEPTHS,
            "local_oracle_ratio_to_rho": LOCAL_RATIO,
            "mean_capacity": MEAN_CAPACITY,
            "coefficient_sha256": COEFFICIENT_SHA256,
            "doubling_charge_in_addition_equivalents": 2,
        }),
        toy_cells,
        native_cells,
        timing,
        accounting,
        fits,
        gates,
        semantic_evidence_sha256: hex::encode(sha256(&semantic)),
        full_depth_unplanted_relation_attempted: false,
        interpretation: "Serpentine sign ordering is an exact table-free stage construction. A sub-rho counted stage does not establish achieved P-256 coverage, the same-family structured degree, relation collection, sparse linear algebra, or a non-generic end-to-end algorithm.".into(),
        decision: "Retain Round 291 as a stage candidate only. Do not claim rho parity and do not attempt an unplanted full-depth relation until correlated-stream coverage and every end-to-end gate are established.".into(),
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
            eprintln!("P-256 serpentine sign-orbit audit failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gray_code_flips_one_bit() {
        let mut previous = 0u64;
        for ordinal in 1..1u64 << FULL_SIGN_BITS {
            let gray = ordinal ^ (ordinal >> 1);
            assert_eq!((gray ^ previous).count_ones(), 1);
            previous = gray;
        }
    }

    #[test]
    fn toy_multisets_match() {
        let cell = toy_cell(257, 10, 3, 3).expect("toy cell");
        assert_eq!(cell.direct_digest, cell.serpentine_digest);
        assert_eq!(cell.multiplicity_mismatches, 0);
        assert_eq!(cell.false_positives, 0);
        assert_eq!(cell.false_negatives, 0);
        assert_eq!(cell.endpoint_failures, 0);
    }

    #[test]
    fn registered_accounting_crosses_stage_boundary() {
        let signs = (1u64 << FULL_SIGN_BITS) as f64;
        let charged = signs * MEAN_CAPACITY + signs - 1.0 + signs + 16.0 + 34.0;
        let states = signs * (MEAN_CAPACITY + 1.0);
        assert!(LOCAL_RATIO * charged / states < 1.0);
    }

    #[test]
    fn residual_capacity_shape_matches_round33() {
        let capacities = residual_capacities(COLUMNS, RARE_EDGES);
        for capacity in 0..18 {
            assert_eq!(
                capacities
                    .iter()
                    .filter(|value| **value == capacity)
                    .count(),
                RARE_EDGES as usize
            );
        }
        assert_eq!(
            capacities.iter().filter(|value| **value == 18).count(),
            6_628
        );
    }
}
