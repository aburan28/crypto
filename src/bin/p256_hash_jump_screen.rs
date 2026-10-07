//! Round 293: exact P-256 two-delta closure and sparse hash-jump screen.

use std::collections::{BTreeMap, BTreeSet};
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
const COLUMNS: u64 = 131_458;
const ARITY: usize = 17;
const SIGN_BITS: u32 = 17;
const SIGNS: u64 = 1 << SIGN_BITS;
const LOCAL_RATIO: f64 = 0.964_336_477_130_181;
const COVERAGE_GATE: f64 = 0.968_617_146_446_504_4;
const ROUND189_RARE: u64 = 7_021;
const ROUND189_CAPACITY: u64 = 226;
const ROUND189_CHARGED: u64 = 29_884_448;
const REGISTERED_RARE: [u64; 3] = [1_021, 2_047, 4_093];
const ROUND189_SHA256: &str = "8630dd6ba930c5df7b3b1dd1307bf2fd2289be1d872f0562b5a408c9c4305d0d";
const MANIFEST_SHA256: &str = "8723757f965eec5cfb1310de1091c8668815212c31823fd30b976de48ebadf45";
const ROUND29_SHA256: &str = "9be6e5fc38644c34ad746dec9a6542da7e35ef2d8e11eb39a83de7af0f54f443";
const ROUND30_SHA256: &str = "8927791b3b149c60a0d4cfb43b77aca21ede9ff21c9c9451a534247ec1f56ad1";
const ROUND31_SHA256: &str = "9b9ee16c5a24e868d1b9304f08aa80ca39ac23b1b4586ef730534e739172c435";

#[derive(Parser)]
#[command(about = "Run the exact P-256 hash-jump factor-base screen")]
struct Cli {
    #[arg(long)]
    round189: PathBuf,
    #[arg(long)]
    executed_manifest: PathBuf,
    #[arg(long)]
    round29: PathBuf,
    #[arg(long)]
    round30: PathBuf,
    #[arg(long)]
    round31: PathBuf,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Clone, Serialize)]
struct Dependency {
    path: String,
    bytes: u64,
    sha256: String,
}

#[derive(Clone, Eq, PartialEq)]
struct Interval {
    start: BigUint,
    end: BigUint,
}

#[derive(Clone)]
struct PreparedTuple {
    attempt: u64,
    rejected: u64,
    cutoff: u64,
    columns: Vec<u64>,
    slot_capacities: Vec<u64>,
    capacity: u64,
    original_coefficients: Vec<BigUint>,
    endpoint_coefficients: Vec<BigUint>,
}

struct HashJumpBase {
    id: String,
    rare_edges: u64,
    jump_attempts: u64,
    rejected_delta_attempts: u64,
    rejected_cycle_attempts: u64,
    rejected_anchor_attempts: u64,
    anchor_attempts: u64,
    anchor: BigUint,
    coefficients: Vec<BigUint>,
    coefficient_sha256: String,
    rare_delta_sha256: String,
    rare_deltas: Vec<BigUint>,
    exact_delta_classes: u64,
    sign_folded_delta_classes: u64,
    edge_replay_failures: u64,
    endpoint_replay_failures: u64,
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

#[derive(Serialize)]
struct LayerCeiling {
    hamming_weight: u32,
    sign_patterns: u64,
    emitted_states: u64,
    band_multiplier: u64,
    quotient_span: u64,
    support_ceiling: u64,
}

#[derive(Serialize)]
struct TwoDeltaCertificate {
    rare_edges: u64,
    columns_replayed: u64,
    coefficient_sha256: String,
    coefficient_digest_matches: bool,
    cycle_closes: bool,
    scaling_failures: u64,
    endpoint_failures: u64,
    scaled_weight_min: u64,
    scaled_weight_max: u64,
    registered_band_width: u64,
    observed_band_width: u64,
    layers: Vec<LayerCeiling>,
    emitted_states: u64,
    actual_distinct_states: u64,
    tuple_independent_support_ceiling: u64,
    tuple_independent_fraction_ceiling: f64,
    ideal_ratio_floor: f64,
    corrected_ratio_floor: f64,
    tuple_search_can_reach_parity: bool,
    certificate_digest: String,
}

#[derive(Clone, Serialize)]
struct DirectControl {
    kind: String,
    candidate_id: String,
    prime_or_order: String,
    sign_bits: u32,
    capacity: u64,
    states_emitted: u64,
    direct_distinct_states: u64,
    interval_distinct_states: u64,
    false_positives: u64,
    false_negatives: u64,
    exact: bool,
}

#[derive(Serialize)]
struct CandidateReceipt {
    id: String,
    rare_edges: u64,
    common_edges: u64,
    jump_attempts: u64,
    rejected_delta_attempts: u64,
    rejected_cycle_attempts: u64,
    rejected_anchor_attempts: u64,
    anchor_attempts: u64,
    anchor_coefficient: String,
    coefficient_sha256: String,
    coefficient_bytes: u64,
    rare_delta_sha256: String,
    rare_delta_bytes: u64,
    exact_delta_classes: u64,
    sign_folded_delta_classes: u64,
    delta_multiplicities: BTreeMap<String, u64>,
    distinct_coefficients: u64,
    unique_up_to_sign: bool,
    nonidentity: bool,
    edge_replay_failures: u64,
    endpoint_replay_failures: u64,
    tuple_attempt: u64,
    tuple_rejected: u64,
    tuple_cutoff: u64,
    tuple_columns: Vec<u64>,
    tuple_slot_capacities: Vec<u64>,
    tuple_capacity: u64,
    tuple_sha256: String,
    sign_segments: u64,
    states_emitted: u64,
    exact_distinct_states: u64,
    duplicates: u64,
    distinct_fraction: f64,
    source_intervals: u64,
    merged_intervals: u64,
    wrap_splits: u64,
    support_digest: String,
    start_digest: String,
    delta_collision_probability: f64,
    hidden_ratio_to_rho: f64,
    charged_addition_equivalents: u64,
    ideal_local_ratio_per_distinct: f64,
    delta_corrected_ratio_per_distinct: f64,
    selector_gate_passed: bool,
    membership_x_coordinates: u64,
    membership_polynomial_degree_lower_bound: u64,
    structured_residual_degree: Option<u32>,
    non_generic_complete_path: bool,
    promoted: bool,
    coefficients_constructed: u64,
    edges_replayed: u64,
    tuple_draws: u64,
    sort_comparisons: u64,
    merge_comparisons: u64,
    construction_wall_ns: u64,
    tuple_wall_ns: u64,
    orbit_generation_wall_ns: u64,
    interval_union_wall_ns: u64,
    total_wall_ns: u64,
    peak_rss_bytes: Option<u64>,
    logical_materialized_bytes: u64,
    interpretation: String,
}

#[derive(Serialize)]
struct Gates {
    dependencies_hash_checked: bool,
    two_delta_certificate_exact: bool,
    two_delta_tuple_search_closed_above_rho: bool,
    all_registered_candidates_complete: bool,
    zero_false_positives_and_false_negatives: bool,
    peak_materialized_storage_below_2_50: bool,
    any_selector_gate_passed: bool,
    same_family_structured_degree_at_most_5: bool,
    per_usable_relation_below_2_103: bool,
    relation_collection_below_2_120: bool,
    complete_non_generic_dlp_below_rho: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultReceipt {
    schema: String,
    curve: String,
    comparison_factor_base: String,
    screening_round: u32,
    execution_status: String,
    dependencies: Vec<Dependency>,
    frozen: Value,
    two_delta_certificate: TwoDeltaCertificate,
    toy_controls: Vec<DirectControl>,
    native_controls: Vec<DirectControl>,
    candidates: Vec<CandidateReceipt>,
    best_selector_candidate: String,
    total_process_cpu_ns: u64,
    peak_rss_bytes: Option<u64>,
    gates: Gates,
    relations_reported: u64,
    full_depth_unplanted_relation_attempted: bool,
    decision: String,
    semantic_evidence_sha256: String,
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

fn load_dependency(path: &Path, expected: &str) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != expected {
        return Err(format!(
            "{} hash mismatch: expected {expected}, got {digest}",
            path.display()
        ));
    }
    let value = serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    Ok((
        Dependency {
            path: path.display().to_string(),
            bytes: bytes.len() as u64,
            sha256: digest,
        },
        value,
    ))
}

fn json_u64(value: &Value, path: &[&str]) -> Result<u64, String> {
    let mut current = value;
    for key in path {
        current = current
            .get(key)
            .ok_or_else(|| format!("missing JSON path {}", path.join(".")))?;
    }
    current
        .as_u64()
        .ok_or_else(|| format!("JSON path {} is not u64", path.join(".")))
}

fn json_str<'a>(value: &'a Value, path: &[&str]) -> Result<&'a str, String> {
    let mut current = value;
    for key in path {
        current = current
            .get(key)
            .ok_or_else(|| format!("missing JSON path {}", path.join(".")))?;
    }
    current
        .as_str()
        .ok_or_else(|| format!("JSON path {} is not string", path.join(".")))
}

fn mechanical_rare(index: u64, rare_edges: u64) -> bool {
    (index + 1) * rare_edges / COLUMNS > index * rare_edges / COLUMNS
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

fn coefficients_digest(coefficients: &[BigUint]) -> Result<String, String> {
    let mut bytes = Vec::with_capacity(coefficients.len() * 32);
    for coefficient in coefficients {
        bytes.extend(fixed_32(coefficient)?);
    }
    Ok(hex::encode(sha256(&bytes)))
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

fn binomial_17(k: u32) -> u64 {
    let k = k.min(17 - k);
    let mut value = 1u64;
    for index in 0..k {
        value = value * u64::from(17 - index) / u64::from(index + 1);
    }
    value
}

fn two_delta_certificate(
    round189: &Value,
    modulus: &BigUint,
) -> Result<TwoDeltaCertificate, String> {
    if json_u64(round189, &["factor_base", "rare_edges"])? != ROUND189_RARE
        || json_u64(round189, &["tuple", "total_capacity"])? != ROUND189_CAPACITY
        || json_u64(round189, &["operations", "charged_addition_equivalents"])? != ROUND189_CHARGED
    {
        return Err("Round 189 frozen values do not match the protocol".into());
    }
    let rare_delta = BigUint::parse_bytes(
        json_str(round189, &["factor_base", "delta_rare"])?.as_bytes(),
        10,
    )
    .ok_or("invalid Round 189 rare delta")?;
    let anchor = BigUint::parse_bytes(
        json_str(round189, &["factor_base", "anchor_coefficient"])?.as_bytes(),
        10,
    )
    .ok_or("invalid Round 189 anchor")?;
    let expected_digest = json_str(round189, &["factor_base", "coefficient_sha256"])?;
    let capacities = residual_capacities(ROUND189_RARE);
    let mut coefficients = Vec::with_capacity(COLUMNS as usize);
    let mut current = anchor.clone();
    let mut scaling_failures = 0u64;
    let mut endpoint_failures = 0u64;
    let mut weight_min = u64::MAX;
    let mut weight_max = 0u64;
    for index in 0..COLUMNS {
        coefficients.push(current.clone());
        let relative = if current >= anchor {
            &current - &anchor
        } else {
            modulus - (&anchor - &current)
        };
        let scaled = (&relative * ROUND189_RARE) % modulus;
        let rare_before = index * ROUND189_RARE / COLUMNS;
        let expected_scaled = index * ROUND189_RARE - rare_before * COLUMNS;
        if scaled != BigUint::from(expected_scaled) {
            scaling_failures += 1;
        }
        let capacity = capacities[index as usize];
        let scaled_weight = 2 * expected_scaled + ROUND189_RARE * capacity;
        weight_min = weight_min.min(scaled_weight);
        weight_max = weight_max.max(scaled_weight);
        let delta = if mechanical_rare(index, ROUND189_RARE) {
            &rare_delta
        } else {
            &BigUint::one()
        };
        add_mod(&mut current, delta, true, modulus);
    }
    let cycle_closes = current == anchor;
    for index in 0..COLUMNS {
        let capacity = capacities[index as usize];
        let endpoint_index = (index + capacity) % COLUMNS;
        let expected = (&coefficients[index as usize] + BigUint::from(capacity)) % modulus;
        if expected != coefficients[endpoint_index as usize]
            || (0..capacity).any(|edge| mechanical_rare((index + edge) % COLUMNS, ROUND189_RARE))
        {
            endpoint_failures += 1;
        }
    }
    let digest = coefficients_digest(&coefficients)?;
    let registered_band_width = COLUMNS + ROUND189_RARE;
    let observed_band_width = weight_max - weight_min;
    if observed_band_width > registered_band_width {
        return Err("observed two-delta band exceeds registered bound".into());
    }

    let mut layers = Vec::new();
    let mut support_ceiling = 0u64;
    let mut emitted_states = 0u64;
    let mut certificate_bytes = Vec::new();
    for k in 0..=17u32 {
        let patterns = binomial_17(k);
        let emitted = patterns * (ROUND189_CAPACITY + 1);
        let band_multiplier = u64::from(k.min(17 - k));
        let numerator = band_multiplier * registered_band_width;
        let quotient_span = numerator.div_ceil(ROUND189_RARE);
        let arithmetic_ceiling = ROUND189_RARE * (ROUND189_CAPACITY + quotient_span + 2);
        let ceiling = emitted.min(arithmetic_ceiling);
        support_ceiling += ceiling;
        emitted_states += emitted;
        certificate_bytes.extend(k.to_be_bytes());
        certificate_bytes.extend(patterns.to_be_bytes());
        certificate_bytes.extend(ceiling.to_be_bytes());
        layers.push(LayerCeiling {
            hamming_weight: k,
            sign_patterns: patterns,
            emitted_states: emitted,
            band_multiplier,
            quotient_span,
            support_ceiling: ceiling,
        });
    }
    let actual_distinct = json_u64(round189, &["coverage", "exact_union_distinct_states"])?;
    let hidden_ratio = round189
        .get("boundary")
        .and_then(|value| value.get("factor_base_complete_ratio_lower_bound"))
        .and_then(Value::as_f64)
        .ok_or("missing Round 189 hidden ratio")?;
    let ideal_ratio_floor = LOCAL_RATIO * ROUND189_CHARGED as f64 / support_ceiling as f64;
    let corrected_ratio_floor = hidden_ratio * ROUND189_CHARGED as f64 / support_ceiling as f64;
    Ok(TwoDeltaCertificate {
        rare_edges: ROUND189_RARE,
        columns_replayed: COLUMNS,
        coefficient_sha256: digest.clone(),
        coefficient_digest_matches: digest == expected_digest,
        cycle_closes,
        scaling_failures,
        endpoint_failures,
        scaled_weight_min: weight_min,
        scaled_weight_max: weight_max,
        registered_band_width,
        observed_band_width,
        layers,
        emitted_states,
        actual_distinct_states: actual_distinct,
        tuple_independent_support_ceiling: support_ceiling,
        tuple_independent_fraction_ceiling: support_ceiling as f64 / emitted_states as f64,
        ideal_ratio_floor,
        corrected_ratio_floor,
        tuple_search_can_reach_parity: ideal_ratio_floor < 1.0 && corrected_ratio_floor < 1.0,
        certificate_digest: hex::encode(sha256(&certificate_bytes)),
    })
}

fn derive_jump(rare_edges: u64, attempt: u64, ordinal: u64, modulus: &BigUint) -> BigUint {
    let domain = format!("{CURVE}/hash-jump-round293/{rare_edges}/{attempt}/{ordinal}");
    BigUint::from_bytes_be(&sha256(domain.as_bytes())) % modulus
}

fn derive_anchor(
    rare_edges: u64,
    jump_attempt: u64,
    anchor_attempt: u64,
    modulus: &BigUint,
) -> BigUint {
    let domain =
        format!("{CURVE}/hash-jump-round293/anchor/{rare_edges}/{jump_attempt}/{anchor_attempt}");
    BigUint::from_bytes_be(&sha256(domain.as_bytes())) % modulus
}

fn construct_hash_jump_base(rare_edges: u64, modulus: &BigUint) -> Result<HashJumpBase, String> {
    let started = Instant::now();
    let mut rejected_delta_attempts = 0u64;
    let mut rejected_cycle_attempts = 0u64;
    let mut rejected_anchor_attempts = 0u64;
    for jump_attempt in 0..4_096u64 {
        let mut rare_deltas = Vec::with_capacity(rare_edges as usize);
        let mut exact_deltas = BTreeSet::new();
        let mut running = BigUint::from(COLUMNS - rare_edges) % modulus;
        let mut delta_invalid = false;
        for ordinal in 0..rare_edges - 1 {
            let delta = derive_jump(rare_edges, jump_attempt, ordinal, modulus);
            if delta.is_zero() || delta.is_one() || !exact_deltas.insert(delta.clone()) {
                delta_invalid = true;
                break;
            }
            running += &delta;
            running %= modulus;
            rare_deltas.push(delta);
        }
        if delta_invalid {
            rejected_delta_attempts += 1;
            continue;
        }
        let closing = if running.is_zero() {
            BigUint::zero()
        } else {
            modulus - &running
        };
        if closing.is_zero() || closing.is_one() || !exact_deltas.insert(closing.clone()) {
            rejected_delta_attempts += 1;
            continue;
        }
        rare_deltas.push(closing);

        let mut relative = Vec::with_capacity(COLUMNS as usize);
        let mut current = BigUint::zero();
        let mut rare_ordinal = 0usize;
        for index in 0..COLUMNS {
            relative.push(current.clone());
            let delta = if mechanical_rare(index, rare_edges) {
                let delta = &rare_deltas[rare_ordinal];
                rare_ordinal += 1;
                delta
            } else {
                &BigUint::one()
            };
            add_mod(&mut current, delta, true, modulus);
        }
        if !current.is_zero()
            || rare_ordinal != rare_edges as usize
            || relative.iter().collect::<BTreeSet<_>>().len() != COLUMNS as usize
        {
            rejected_cycle_attempts += 1;
            continue;
        }

        let mut accepted = None;
        let mut anchor_attempts = 0u64;
        for anchor_attempt in 0..4_096u64 {
            anchor_attempts += 1;
            let anchor = derive_anchor(rare_edges, jump_attempt, anchor_attempt, modulus);
            let coefficients = relative
                .iter()
                .map(|value| (value + &anchor) % modulus)
                .collect::<Vec<_>>();
            if coefficients.iter().any(BigUint::is_zero) {
                rejected_anchor_attempts += 1;
                continue;
            }
            let signed = coefficients
                .iter()
                .map(|value| value.clone().min(modulus - value))
                .collect::<BTreeSet<_>>();
            if signed.len() != COLUMNS as usize {
                rejected_anchor_attempts += 1;
                continue;
            }
            accepted = Some((anchor, coefficients, anchor_attempts));
            break;
        }
        let Some((anchor, coefficients, anchor_attempts)) = accepted else {
            rejected_cycle_attempts += 1;
            continue;
        };

        let capacities = residual_capacities(rare_edges);
        let mut edge_replay_failures = 0u64;
        let mut endpoint_replay_failures = 0u64;
        let mut rare_ordinal = 0usize;
        for index in 0..COLUMNS as usize {
            let delta = if mechanical_rare(index as u64, rare_edges) {
                let delta = &rare_deltas[rare_ordinal];
                rare_ordinal += 1;
                delta
            } else {
                &BigUint::one()
            };
            if (&coefficients[index] + delta) % modulus
                != coefficients[(index + 1) % COLUMNS as usize]
            {
                edge_replay_failures += 1;
            }
            let capacity = capacities[index];
            let endpoint_index = (index as u64 + capacity) % COLUMNS;
            let expected = (&coefficients[index] + BigUint::from(capacity)) % modulus;
            if expected != coefficients[endpoint_index as usize] {
                endpoint_replay_failures += 1;
            }
        }

        let coefficient_sha256 = coefficients_digest(&coefficients)?;
        let mut delta_bytes = Vec::with_capacity(rare_deltas.len() * 32);
        for delta in &rare_deltas {
            delta_bytes.extend(fixed_32(delta)?);
        }
        let rare_delta_sha256 = hex::encode(sha256(&delta_bytes));
        let sign_folded_delta_classes = rare_deltas
            .iter()
            .map(|delta| delta.clone().min(modulus - delta))
            .chain(std::iter::once(BigUint::one()))
            .collect::<BTreeSet<_>>()
            .len() as u64;
        let id = format!("LDHJR{rare_edges}h{}", &coefficient_sha256[..12]);
        return Ok(HashJumpBase {
            id,
            rare_edges,
            jump_attempts: jump_attempt + 1,
            rejected_delta_attempts,
            rejected_cycle_attempts,
            rejected_anchor_attempts,
            anchor_attempts,
            anchor,
            coefficients,
            coefficient_sha256,
            rare_delta_sha256,
            rare_deltas,
            exact_delta_classes: rare_edges + 1,
            sign_folded_delta_classes,
            edge_replay_failures,
            endpoint_replay_failures,
            construction_wall_ns: started.elapsed().as_nanos() as u64,
        });
    }
    Err(format!(
        "no accepted hash-jump base for rare edge count {rare_edges}"
    ))
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

fn select_tuple(base: &HashJumpBase, modulus: &BigUint) -> Result<PreparedTuple, String> {
    let capacities = residual_capacities(base.rare_edges);
    let cutoff = (219 * 6_935u64).div_ceil(base.rare_edges);
    for attempt in 0..1_000_000u64 {
        let columns = draw_columns(attempt);
        let slot_capacities = columns
            .iter()
            .map(|column| capacities[*column as usize])
            .collect::<Vec<_>>();
        let capacity = slot_capacities.iter().sum::<u64>();
        if capacity < cutoff {
            continue;
        }
        let original_coefficients = columns
            .iter()
            .map(|column| base.coefficients[*column as usize].clone())
            .collect::<Vec<_>>();
        let mut endpoint_coefficients = Vec::with_capacity(ARITY);
        for ((column, slot_capacity), original) in columns
            .iter()
            .zip(&slot_capacities)
            .zip(&original_coefficients)
        {
            let endpoint_column = (column + slot_capacity) % COLUMNS;
            let endpoint = base.coefficients[endpoint_column as usize].clone();
            if (original + BigUint::from(*slot_capacity)) % modulus != endpoint {
                return Err(format!("tuple endpoint mismatch for {}", base.id));
            }
            endpoint_coefficients.push(endpoint);
        }
        return Ok(PreparedTuple {
            attempt,
            rejected: attempt,
            cutoff,
            columns,
            slot_capacities,
            capacity,
            original_coefficients,
            endpoint_coefficients,
        });
    }
    Err(format!("no tuple reached cutoff {cutoff} for {}", base.id))
}

fn tuple_digest(tuple: &PreparedTuple) -> String {
    let mut bytes = Vec::with_capacity(ARITY * 16);
    for (column, capacity) in tuple.columns.iter().zip(&tuple.slot_capacities) {
        bytes.extend(column.to_be_bytes());
        bytes.extend(capacity.to_be_bytes());
    }
    hex::encode(sha256(&bytes))
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
    let high = &low + BigUint::from(width);
    if &high < modulus {
        (
            vec![Interval {
                start: low,
                end: high,
            }],
            false,
        )
    } else {
        (
            vec![
                Interval {
                    start: BigUint::zero(),
                    end: &high - modulus,
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
        intervals.append(&mut pieces);
        wrap_splits += u64::from(wrapped);
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

fn direct_support(tuple: &PreparedTuple, sign_bits: u32, modulus: &BigUint) -> BTreeSet<BigUint> {
    let signs = 1u64 << sign_bits;
    let (mut forward_corner, mut reverse_corner) = initial_corners(tuple, modulus);
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
            add_mod(&mut state, &BigUint::one(), forward, modulus);
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

fn interval_contains(intervals: &[Interval], value: &BigUint) -> bool {
    let index = intervals.partition_point(|interval| interval.end < *value);
    index < intervals.len() && intervals[index].start <= *value
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
        rejected: 0,
        cutoff: 0,
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

fn direct_control(
    kind: &str,
    candidate_id: &str,
    tuple: &PreparedTuple,
    sign_bits: u32,
    modulus: &BigUint,
) -> Result<DirectControl, String> {
    let interval = monotone_orbit(tuple, sign_bits, modulus)?;
    let direct = direct_support(tuple, sign_bits, modulus);
    let (false_positives, false_negatives) = compare_direct(&interval, &direct);
    Ok(DirectControl {
        kind: kind.into(),
        candidate_id: candidate_id.into(),
        prime_or_order: modulus.to_string(),
        sign_bits,
        capacity: tuple.capacity,
        states_emitted: interval.states_emitted,
        direct_distinct_states: direct.len() as u64,
        interval_distinct_states: interval.exact_distinct_states,
        false_positives,
        false_negatives,
        exact: false_positives == 0 && false_negatives == 0,
    })
}

fn truncated_tuple(tuple: &PreparedTuple, modulus: &BigUint) -> PreparedTuple {
    let slot_capacities = tuple
        .slot_capacities
        .iter()
        .map(|capacity| (*capacity).min(2))
        .collect::<Vec<_>>();
    let endpoint_coefficients = tuple
        .original_coefficients
        .iter()
        .zip(&slot_capacities)
        .map(|(coefficient, capacity)| (coefficient + BigUint::from(*capacity)) % modulus)
        .collect::<Vec<_>>();
    PreparedTuple {
        attempt: tuple.attempt,
        rejected: tuple.rejected,
        cutoff: tuple.cutoff,
        columns: tuple.columns.clone(),
        capacity: slot_capacities.iter().sum(),
        slot_capacities,
        original_coefficients: tuple.original_coefficients.clone(),
        endpoint_coefficients,
    }
}

fn charged_additions(capacity: u64) -> u64 {
    SIGNS * capacity + (SIGNS - 1) + (ARITY - 1) as u64 + ARITY as u64 + SIGNS
}

fn candidate_receipt(
    base: HashJumpBase,
    tuple: &PreparedTuple,
    support: OrbitSupport,
) -> CandidateReceipt {
    let common_edges = COLUMNS - base.rare_edges;
    let numerator = (common_edges as u128) * (common_edges as u128) + u128::from(base.rare_edges);
    let denominator = (COLUMNS as u128) * (COLUMNS as u128);
    let delta_collision_probability = numerator as f64 / denominator as f64;
    let hidden_ratio_to_rho = LOCAL_RATIO / delta_collision_probability.sqrt();
    let charged = charged_additions(tuple.capacity);
    let distinct_fraction = support.exact_distinct_states as f64 / support.states_emitted as f64;
    let ideal_local_ratio_per_distinct =
        LOCAL_RATIO * charged as f64 / support.exact_distinct_states as f64;
    let delta_corrected_ratio_per_distinct =
        hidden_ratio_to_rho * charged as f64 / support.exact_distinct_states as f64;
    let selector_gate_passed = distinct_fraction >= COVERAGE_GATE
        && ideal_local_ratio_per_distinct < 1.0
        && delta_corrected_ratio_per_distinct < 1.0;
    let mut multiplicities = BTreeMap::new();
    multiplicities.insert("1".into(), common_edges);
    for delta in &base.rare_deltas {
        multiplicities.insert(delta.to_string(), 1);
    }
    let logical_materialized_bytes = COLUMNS * 32
        + base.rare_edges * 32
        + support.source_intervals * 64
        + support.merged_intervals * 64
        + COLUMNS * 8;
    CandidateReceipt {
        id: base.id,
        rare_edges: base.rare_edges,
        common_edges,
        jump_attempts: base.jump_attempts,
        rejected_delta_attempts: base.rejected_delta_attempts,
        rejected_cycle_attempts: base.rejected_cycle_attempts,
        rejected_anchor_attempts: base.rejected_anchor_attempts,
        anchor_attempts: base.anchor_attempts,
        anchor_coefficient: base.anchor.to_string(),
        coefficient_sha256: base.coefficient_sha256,
        coefficient_bytes: COLUMNS * 32,
        rare_delta_sha256: base.rare_delta_sha256,
        rare_delta_bytes: base.rare_edges * 32,
        exact_delta_classes: base.exact_delta_classes,
        sign_folded_delta_classes: base.sign_folded_delta_classes,
        delta_multiplicities: multiplicities,
        distinct_coefficients: COLUMNS,
        unique_up_to_sign: true,
        nonidentity: true,
        edge_replay_failures: base.edge_replay_failures,
        endpoint_replay_failures: base.endpoint_replay_failures,
        tuple_attempt: tuple.attempt,
        tuple_rejected: tuple.rejected,
        tuple_cutoff: tuple.cutoff,
        tuple_columns: tuple.columns.clone(),
        tuple_slot_capacities: tuple.slot_capacities.clone(),
        tuple_capacity: tuple.capacity,
        tuple_sha256: tuple_digest(tuple),
        sign_segments: SIGNS,
        states_emitted: support.states_emitted,
        exact_distinct_states: support.exact_distinct_states,
        duplicates: support.states_emitted - support.exact_distinct_states,
        distinct_fraction,
        source_intervals: support.source_intervals,
        merged_intervals: support.merged_intervals,
        wrap_splits: support.wrap_splits,
        support_digest: support.support_digest,
        start_digest: support.start_digest,
        delta_collision_probability,
        hidden_ratio_to_rho,
        charged_addition_equivalents: charged,
        ideal_local_ratio_per_distinct,
        delta_corrected_ratio_per_distinct,
        selector_gate_passed,
        membership_x_coordinates: COLUMNS,
        membership_polynomial_degree_lower_bound: COLUMNS,
        structured_residual_degree: None,
        non_generic_complete_path: false,
        promoted: false,
        coefficients_constructed: COLUMNS,
        edges_replayed: COLUMNS,
        tuple_draws: tuple.attempt + 1,
        sort_comparisons: support.sort_comparisons,
        merge_comparisons: support.merge_comparisons,
        construction_wall_ns: base.construction_wall_ns,
        tuple_wall_ns: 0,
        orbit_generation_wall_ns: support.generation_wall_ns,
        interval_union_wall_ns: support.union_wall_ns,
        total_wall_ns: 0,
        peak_rss_bytes: peak_rss_bytes(),
        logical_materialized_bytes,
        interpretation: "Complete selector-stage coverage measurement. The explicit random x-set has no verified low-degree membership structure, and the constant +G path is a group translation; this is not a non-generic index-calculus DLP.".into(),
    }
}

fn semantic_digest(
    certificate: &TwoDeltaCertificate,
    candidates: &[CandidateReceipt],
    controls: &[DirectControl],
) -> String {
    let mut hasher = Hasher::new();
    hasher.update(&certificate.tuple_independent_support_ceiling.to_be_bytes());
    hasher.update(&certificate.corrected_ratio_floor.to_bits().to_be_bytes());
    for candidate in candidates {
        hasher.update(candidate.id.as_bytes());
        hasher.update(&candidate.exact_distinct_states.to_be_bytes());
        hasher.update(
            &candidate
                .delta_corrected_ratio_per_distinct
                .to_bits()
                .to_be_bytes(),
        );
        hasher.update(candidate.support_digest.as_bytes());
    }
    for control in controls {
        hasher.update(control.candidate_id.as_bytes());
        hasher.update(&control.sign_bits.to_be_bytes());
        hasher.update(&control.false_positives.to_be_bytes());
        hasher.update(&control.false_negatives.to_be_bytes());
    }
    hasher.finalize().to_hex().to_string()
}

fn run(cli: Cli) -> Result<(), String> {
    if cli.out.exists() {
        return Err(format!("refusing to overwrite {}", cli.out.display()));
    }
    let (dep_round189, round189) = load_dependency(&cli.round189, ROUND189_SHA256)?;
    let (dep_manifest, manifest) = load_dependency(&cli.executed_manifest, MANIFEST_SHA256)?;
    let (dep_round29, round29) = load_dependency(&cli.round29, ROUND29_SHA256)?;
    let (dep_round30, round30) = load_dependency(&cli.round30, ROUND30_SHA256)?;
    let (dep_round31, round31) = load_dependency(&cli.round31, ROUND31_SHA256)?;
    if json_str(&round189, &["curve"])? != CURVE
        || json_str(&manifest, &["curve"])? != CURVE
        || json_str(&round29, &["curve"])? != CURVE
        || json_str(&round30, &["curve"])? != CURVE
        || json_str(&round31, &["curve"])? != CURVE
    {
        return Err("dependency curve identity mismatch".into());
    }
    let dependencies = vec![
        dep_round189,
        dep_manifest,
        dep_round29,
        dep_round30,
        dep_round31,
    ];
    let curve = CurveParams::p256();
    let process_started = process_cpu_ns();
    let certificate = two_delta_certificate(&round189, &curve.n)?;
    let certificate_exact = certificate.coefficient_digest_matches
        && certificate.cycle_closes
        && certificate.scaling_failures == 0
        && certificate.endpoint_failures == 0
        && certificate.actual_distinct_states <= certificate.tuple_independent_support_ceiling;
    if !certificate_exact {
        return Err("two-delta certificate failed exact replay".into());
    }

    let toy_specs = [
        (257u64, vec![3, 29, 101], vec![3, 4, 2]),
        (
            65_537u64,
            vec![7, 101, 9_997, 31_337, 50_003],
            vec![2, 5, 3, 4, 6],
        ),
    ];
    let mut toy_controls = Vec::new();
    for (prime, coefficients, capacities) in toy_specs {
        let tuple = toy_tuple(&coefficients, &capacities, prime);
        toy_controls.push(direct_control(
            "toy-complete",
            &format!("toy-p{prime}"),
            &tuple,
            coefficients.len() as u32,
            &BigUint::from(prime),
        )?);
    }

    let mut candidates = Vec::new();
    let mut native_controls = Vec::new();
    for rare_edges in REGISTERED_RARE {
        let candidate_started = Instant::now();
        let base = construct_hash_jump_base(rare_edges, &curve.n)?;
        if base.edge_replay_failures != 0 || base.endpoint_replay_failures != 0 {
            return Err(format!("{} exact replay failed", base.id));
        }
        let tuple_started = Instant::now();
        let tuple = select_tuple(&base, &curve.n)?;
        let tuple_wall_ns = tuple_started.elapsed().as_nanos() as u64;
        let control_tuple = truncated_tuple(&tuple, &curve.n);
        for sign_bits in [8, 10, 12] {
            native_controls.push(direct_control(
                "native-truncated-capacity",
                &base.id,
                &control_tuple,
                sign_bits,
                &curve.n,
            )?);
        }
        let support = monotone_orbit(&tuple, SIGN_BITS, &curve.n)?;
        let mut receipt = candidate_receipt(base, &tuple, support);
        receipt.tuple_wall_ns = tuple_wall_ns;
        receipt.total_wall_ns = candidate_started.elapsed().as_nanos() as u64;
        candidates.push(receipt);
    }
    let controls_exact = toy_controls.iter().chain(&native_controls).all(|control| {
        control.exact && control.false_positives == 0 && control.false_negatives == 0
    });
    if !controls_exact {
        return Err("a direct-enumeration control failed".into());
    }
    let best = candidates
        .iter()
        .min_by(|left, right| {
            left.delta_corrected_ratio_per_distinct
                .total_cmp(&right.delta_corrected_ratio_per_distinct)
                .then_with(|| left.rare_edges.cmp(&right.rare_edges))
        })
        .ok_or("no candidates")?;
    let best_selector_candidate = best.id.clone();
    let any_selector_gate_passed = candidates
        .iter()
        .any(|candidate| candidate.selector_gate_passed);
    let storage_gate = candidates
        .iter()
        .all(|candidate| candidate.logical_materialized_bytes < (1u64 << 50));
    let all_complete = candidates.len() == REGISTERED_RARE.len()
        && candidates.iter().all(|candidate| {
            candidate.edge_replay_failures == 0
                && candidate.endpoint_replay_failures == 0
                && candidate.states_emitted == SIGNS * (candidate.tuple_capacity + 1)
        });
    let gates = Gates {
        dependencies_hash_checked: true,
        two_delta_certificate_exact: certificate_exact,
        two_delta_tuple_search_closed_above_rho: !certificate.tuple_search_can_reach_parity
            && certificate.corrected_ratio_floor > 1.0,
        all_registered_candidates_complete: all_complete,
        zero_false_positives_and_false_negatives: controls_exact,
        peak_materialized_storage_below_2_50: storage_gate,
        any_selector_gate_passed,
        same_family_structured_degree_at_most_5: false,
        per_usable_relation_below_2_103: false,
        relation_collection_below_2_120: false,
        complete_non_generic_dlp_below_rho: false,
        promoted: false,
    };
    let all_controls = toy_controls
        .iter()
        .chain(&native_controls)
        .cloned()
        .collect::<Vec<_>>();
    let semantic_evidence_sha256 = semantic_digest(&certificate, &candidates, &all_controls);
    let receipt = ResultReceipt {
        schema: "p256-hash-jump-round293/v1".into(),
        curve: CURVE.into(),
        comparison_factor_base: COMPARISON_FB.into(),
        screening_round: 293,
        execution_status: "complete".into(),
        dependencies,
        frozen: serde_json::json!({
            "columns": COLUMNS,
            "arity": ARITY,
            "sign_bits": SIGN_BITS,
            "registered_rare_edges": REGISTERED_RARE,
            "local_ratio": LOCAL_RATIO,
            "coverage_gate": COVERAGE_GATE,
            "round189_capacity": ROUND189_CAPACITY,
            "round189_charged_additions": ROUND189_CHARGED,
            "tuple_cutoff_formula": "ceil(219*6935/R)",
            "native_controls_truncate_each_slot_capacity_to": 2,
        }),
        two_delta_certificate: certificate,
        toy_controls,
        native_controls,
        candidates,
        best_selector_candidate,
        total_process_cpu_ns: process_cpu_ns().saturating_sub(process_started),
        peak_rss_bytes: peak_rss_bytes(),
        gates,
        relations_reported: 0,
        full_depth_unplanted_relation_attempted: false,
        decision: if any_selector_gate_passed {
            "A decorrelated candidate crosses the selector-stage boundary, but promotion fails: no same-family degree<=5 evidence, relation/linear-algebra cost, or non-generic complete path exists. Do not claim rho parity or attempt an unplanted full-depth relation."
        } else {
            "No candidate crosses the selector-stage boundary. Close both tuple-only repair and the registered sparse hash-jump screen; do not attempt an unplanted full-depth relation."
        }
        .into(),
        semantic_evidence_sha256,
    };
    let bytes = serde_json::to_vec_pretty(&receipt).map_err(|error| error.to_string())?;
    if let Some(parent) = cli.out.parent() {
        fs::create_dir_all(parent).map_err(|error| error.to_string())?;
    }
    fs::write(&cli.out, bytes).map_err(|error| format!("{}: {error}", cli.out.display()))?;
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binomial_layers_sum_to_all_signs() {
        assert_eq!((0..=17).map(binomial_17).sum::<u64>(), SIGNS);
    }

    #[test]
    fn toy_interval_union_matches_direct_enumeration() {
        for (prime, coefficients, capacities) in [
            (257u64, vec![3, 29, 101], vec![3, 4, 2]),
            (65_537u64, vec![7, 101, 9_997, 31_337], vec![2, 5, 3, 4]),
        ] {
            let tuple = toy_tuple(&coefficients, &capacities, prime);
            let control = direct_control(
                "test",
                "toy",
                &tuple,
                coefficients.len() as u32,
                &BigUint::from(prime),
            )
            .unwrap();
            assert!(control.exact);
            assert_eq!(control.false_positives, 0);
            assert_eq!(control.false_negatives, 0);
        }
    }

    #[test]
    fn mechanical_word_has_exact_rare_count() {
        for rare in REGISTERED_RARE.into_iter().chain([ROUND189_RARE]) {
            assert_eq!(
                (0..COLUMNS)
                    .filter(|index| mechanical_rare(*index, rare))
                    .count() as u64,
                rare
            );
        }
    }
}
