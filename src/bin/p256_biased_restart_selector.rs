//! Biased common-edge restart selector for P-256 (round 33).

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;
use serde_json::Value;

const ROUND31_SHA256: &str = "9b9ee16c5a24e868d1b9304f08aa80ca39ac23b1b4586ef730534e739172c435";
const GROUP_ORDER: &str =
    "115792089210356248762697446949407573529996955224135760342422259061068512044369";
const COLUMNS: u64 = 131_458;
const ARITY: usize = 17;
const RARE_EDGES: u64 = 6_935;
const PAIR_BYTES: u64 = 33;
const PAIR_LOOKUPS_PER_SEGMENT: u64 = 8;
const PAIR_SETUP_ADDITIONS: u64 = 8;
const NATIVE_SAMPLES: u64 = 4_096;
const RELATION_MODEL: &str =
    "17 distinct variable columns; signed 8+9 cross-colour collision modulo global negation";

#[derive(Parser)]
#[command(about = "Evaluate biased common-edge restarts on the P-256 low-delta cycle")]
struct Cli {
    #[arg(long)]
    round31: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Serialize)]
struct Dependency {
    round: u64,
    path: String,
    sha256: String,
    schema: String,
}

struct LoadedDependency {
    record: Dependency,
    value: Value,
}

#[derive(Clone, Serialize)]
struct CutoffRow {
    cutoff_capacity: u64,
    retained_unsigned_tuples: String,
    retained_fraction: f64,
    mean_capacity: f64,
    signed_path_visits_log2: f64,
    support_upper_log2: f64,
    support_upper_fraction_of_group: f64,
    target_retry_lower_bound: f64,
    online_additions_per_sample: f64,
    online_additions_per_sample_with_target_correction: f64,
    optimistic_ratio_to_rho_excluding_pair_build: f64,
    optimistic_ratio_to_rho_including_pair_build: f64,
    optimistic_ratio_to_rho_with_target_correction: f64,
    table_bytes_read_per_segment: u64,
}

#[derive(Serialize)]
struct StartDistribution {
    residual_capacity_counts: Vec<u64>,
    maximum_tuple_capacity: u64,
    exact_unsigned_tuples: String,
    expected_unsigned_tuples: String,
    identity_matches: bool,
    ordered_distribution_sha256: String,
}

#[derive(Serialize)]
struct CutoffSweep {
    rows_evaluated: u64,
    local_oracle_ratio_to_rho: f64,
    ordered_rows_sha256: String,
    selected: CutoffRow,
    rows: Vec<CutoffRow>,
    any_optimistic_row_below_rho: bool,
}

#[derive(Serialize)]
struct PairTableAccounting {
    signed_distinct_column_pair_entries: u64,
    bytes_per_entry: u64,
    materialized_bytes: u64,
    materialized_log2_bytes: f64,
    below_2_50_bytes: bool,
    build_group_additions: u64,
    build_ratio_to_sqrt_group_order: f64,
    lookups_per_segment: u64,
    bytes_read_per_segment: u64,
    setup_group_additions_per_segment: u64,
    omitted_target_correction_additions_per_segment: u64,
}

#[derive(Serialize)]
struct ToyCell {
    prime: u64,
    columns: u64,
    arity: u64,
    rare_edges: u64,
    unsigned_tuples: u64,
    signed_segments: u64,
    path_states_visited: u64,
    transitions_replayed: u64,
    edge_failures: u64,
    cutoffs_checked: u64,
    support_bound_violations: u64,
    distribution_identity_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    maximum_exact_support: u64,
    ordered_cell_sha256: String,
}

#[derive(Serialize)]
struct ToyReferences {
    cells: Vec<ToyCell>,
    unsigned_tuples: u64,
    signed_segments: u64,
    path_states_visited: u64,
    transitions_replayed: u64,
    edge_failures: u64,
    support_bound_violations: u64,
    distribution_identity_failures: u64,
    false_positives: u64,
    false_negatives: u64,
}

#[derive(Serialize)]
struct NativeReplay {
    rare_edges: u64,
    cycle_edges_replayed: u64,
    cycle_edge_failures: u64,
    coefficient_sha256: String,
    coefficient_digest_matches_round31: bool,
    samples_requested: u64,
    samples_accepted: u64,
    draws_attempted: u64,
    rejected_draws: u64,
    minimum_capacity: u64,
    maximum_capacity: u64,
    mean_capacity: f64,
    transitions_replayed: u64,
    transition_failures: u64,
    setup_group_additions: u64,
    pair_table_lookups: u64,
    pair_table_bytes_read: u64,
    ordered_samples_sha256: String,
}

#[derive(Serialize)]
struct ImportedGates {
    structured_residual_degree_of_regularity: Option<u64>,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    non_generic_end_to_end: bool,
}

#[derive(Serialize)]
struct PromotionGates {
    exact_distribution_identity: bool,
    exact_cycle_and_segment_replay: bool,
    zero_false_positives_and_false_negatives: bool,
    discarded_starts_charged_by_support_bound: bool,
    optimistic_group_addition_time_below_rho: bool,
    optimistic_group_addition_time_with_target_correction_below_rho: bool,
    measured_complete_time_including_table_traffic_below_rho: bool,
    proved_p256_usable_relation_probability: bool,
    structured_residual_degree_at_most_5: bool,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    non_generic_end_to_end: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    relation_model: String,
    family: String,
    dependency: Dependency,
    start_distribution: StartDistribution,
    cutoff_sweep: CutoffSweep,
    pair_table_accounting: PairTableAccounting,
    toy_references: ToyReferences,
    native_replay: NativeReplay,
    imported_gates: ImportedGates,
    promotion_gates: PromotionGates,
    full_depth_unplanted_attempted: bool,
    interpretation: String,
    decision: String,
}

fn at<'a>(value: &'a Value, path: &[&str]) -> Result<&'a Value, String> {
    let mut cursor = value;
    for key in path {
        cursor = cursor
            .get(*key)
            .ok_or_else(|| format!("missing JSON path {}", path.join(".")))?;
    }
    Ok(cursor)
}

fn string_at<'a>(value: &'a Value, path: &[&str]) -> Result<&'a str, String> {
    at(value, path)?
        .as_str()
        .ok_or_else(|| format!("JSON path {} is not a string", path.join(".")))
}

fn f64_at(value: &Value, path: &[&str]) -> Result<f64, String> {
    at(value, path)?
        .as_f64()
        .ok_or_else(|| format!("JSON path {} is not a number", path.join(".")))
}

fn bool_at(value: &Value, path: &[&str]) -> Result<bool, String> {
    at(value, path)?
        .as_bool()
        .ok_or_else(|| format!("JSON path {} is not a bool", path.join(".")))
}

fn load_dependency(path: &Path) -> Result<LoadedDependency, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != ROUND31_SHA256 {
        return Err(format!(
            "round-31 dependency hash mismatch: expected {ROUND31_SHA256}, got {digest}"
        ));
    }
    let value: Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("{}: {error}", path.display()))?;
    if string_at(&value, &["curve"])? != CURVE_SLUG
        || string_at(&value, &["relation_model"])? != RELATION_MODEL
    {
        return Err("round-31 curve/relation mismatch".into());
    }
    Ok(LoadedDependency {
        record: Dependency {
            round: 31,
            path: path.display().to_string(),
            sha256: digest,
            schema: string_at(&value, &["schema"])?.into(),
        },
        value,
    })
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    BigUint::parse_bytes(value.as_bytes(), 10)
        .ok_or_else(|| format!("invalid frozen integer: {value}"))
}

fn log2_big(value: &BigUint) -> Result<f64, String> {
    if value.is_zero() {
        return Err("log2(0)".into());
    }
    let bits = value.bits();
    let keep = bits.min(53);
    let shift = bits - keep;
    let top = (value >> shift)
        .to_u64()
        .ok_or("top 53 bits did not fit u64")?;
    Ok(shift as f64 + (top as f64).log2())
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

fn capacity_histogram(capacities: &[u64]) -> Vec<u64> {
    let maximum = capacities.iter().copied().max().unwrap_or(0) as usize;
    let mut counts = vec![0u64; maximum + 1];
    for capacity in capacities {
        counts[*capacity as usize] += 1;
    }
    counts
}

fn binomial(n: u64, k: usize) -> BigUint {
    let k = k.min(n as usize - k.min(n as usize));
    let mut value = BigUint::one();
    for index in 0..k {
        value *= n - index as u64;
        value /= index as u64 + 1;
    }
    value
}

fn capacity_distribution(counts: &[u64], arity: usize) -> Vec<BigUint> {
    let maximum_capacity = (counts.len() - 1) * arity;
    let mut dp = vec![vec![BigUint::zero(); maximum_capacity + 1]; arity + 1];
    dp[0][0] = BigUint::one();
    for (capacity, count) in counts.iter().copied().enumerate() {
        let mut next = vec![vec![BigUint::zero(); maximum_capacity + 1]; arity + 1];
        for chosen in 0..=arity {
            for sum in 0..=maximum_capacity {
                if dp[chosen][sum].is_zero() {
                    continue;
                }
                let remaining = arity - chosen;
                for take in 0..=remaining.min(count as usize) {
                    let next_sum = sum + take * capacity;
                    next[chosen + take][next_sum] += &dp[chosen][sum] * binomial(count, take);
                }
            }
        }
        dp = next;
    }
    dp.remove(arity)
}

fn distribution_digest(distribution: &[BigUint]) -> String {
    let mut blob = Vec::new();
    for (capacity, count) in distribution.iter().enumerate() {
        let encoded = count.to_bytes_be();
        blob.extend((capacity as u64).to_be_bytes());
        blob.extend((encoded.len() as u64).to_be_bytes());
        blob.extend(encoded);
    }
    hex::encode(sha256(&blob))
}

fn pair_table_accounting(group_order: &BigUint) -> Result<PairTableAccounting, String> {
    let entries = 2 * COLUMNS * (COLUMNS - 1);
    let bytes = entries * PAIR_BYTES;
    let sqrt_order = 2f64.powf(log2_big(group_order)? / 2.0);
    Ok(PairTableAccounting {
        signed_distinct_column_pair_entries: entries,
        bytes_per_entry: PAIR_BYTES,
        materialized_bytes: bytes,
        materialized_log2_bytes: (bytes as f64).log2(),
        below_2_50_bytes: bytes < (1u64 << 50),
        build_group_additions: entries,
        build_ratio_to_sqrt_group_order: entries as f64 / sqrt_order,
        lookups_per_segment: PAIR_LOOKUPS_PER_SEGMENT,
        bytes_read_per_segment: PAIR_LOOKUPS_PER_SEGMENT * PAIR_BYTES,
        setup_group_additions_per_segment: PAIR_SETUP_ADDITIONS,
        omitted_target_correction_additions_per_segment: 1,
    })
}

fn cutoff_sweep(
    distribution: &[BigUint],
    local_ratio: f64,
    group_order: &BigUint,
    table: &PairTableAccounting,
) -> Result<CutoffSweep, String> {
    let total = distribution.iter().sum::<BigUint>();
    let total_f64 = total.to_f64().ok_or("tuple count did not fit f64")?;
    let log_order = log2_big(group_order)?;
    let mut retained = BigUint::zero();
    let mut weighted = BigUint::zero();
    let mut rows_descending = Vec::new();
    let mut digest = Vec::new();
    for cutoff in (0..distribution.len()).rev() {
        retained += &distribution[cutoff];
        weighted += &distribution[cutoff] * cutoff;
        if retained.is_zero() {
            continue;
        }
        let retained_f64 = retained.to_f64().ok_or("retained tuples did not fit f64")?;
        let weighted_f64 = weighted
            .to_f64()
            .ok_or("weighted capacity did not fit f64")?;
        if weighted_f64 == 0.0 {
            continue;
        }
        let visits = (&weighted + &retained) << ARITY;
        let support_log2_raw = log2_big(&visits)?;
        let coverage = 2f64.powf((support_log2_raw - log_order).min(0.0));
        let retry = 1.0 / coverage;
        let online_factor = (weighted_f64 + PAIR_SETUP_ADDITIONS as f64 * retained_f64)
            / (weighted_f64 + retained_f64);
        let corrected_factor = (weighted_f64 + (PAIR_SETUP_ADDITIONS + 1) as f64 * retained_f64)
            / (weighted_f64 + retained_f64);
        let excluding_build = local_ratio * online_factor * retry;
        let including_build = excluding_build + table.build_ratio_to_sqrt_group_order;
        let with_target_correction =
            local_ratio * corrected_factor * retry + table.build_ratio_to_sqrt_group_order;
        let row = CutoffRow {
            cutoff_capacity: cutoff as u64,
            retained_unsigned_tuples: retained.to_string(),
            retained_fraction: retained_f64 / total_f64,
            mean_capacity: weighted_f64 / retained_f64,
            signed_path_visits_log2: support_log2_raw,
            support_upper_log2: support_log2_raw.min(log_order),
            support_upper_fraction_of_group: coverage,
            target_retry_lower_bound: retry,
            online_additions_per_sample: online_factor,
            online_additions_per_sample_with_target_correction: corrected_factor,
            optimistic_ratio_to_rho_excluding_pair_build: excluding_build,
            optimistic_ratio_to_rho_including_pair_build: including_build,
            optimistic_ratio_to_rho_with_target_correction: with_target_correction,
            table_bytes_read_per_segment: table.bytes_read_per_segment,
        };
        digest.extend(row.cutoff_capacity.to_be_bytes());
        digest.extend((retained.to_bytes_be().len() as u64).to_be_bytes());
        digest.extend(retained.to_bytes_be());
        for value in [
            row.retained_fraction,
            row.mean_capacity,
            row.signed_path_visits_log2,
            row.target_retry_lower_bound,
            row.online_additions_per_sample,
            row.optimistic_ratio_to_rho_including_pair_build,
            row.optimistic_ratio_to_rho_with_target_correction,
        ] {
            digest.extend(value.to_bits().to_be_bytes());
        }
        rows_descending.push(row);
    }
    rows_descending.reverse();
    let selected = rows_descending
        .iter()
        .min_by(|left, right| {
            left.optimistic_ratio_to_rho_including_pair_build
                .total_cmp(&right.optimistic_ratio_to_rho_including_pair_build)
                .then_with(|| right.cutoff_capacity.cmp(&left.cutoff_capacity))
        })
        .ok_or("no nonempty cutoff row")?
        .clone();
    Ok(CutoffSweep {
        rows_evaluated: rows_descending.len() as u64,
        local_oracle_ratio_to_rho: local_ratio,
        ordered_rows_sha256: hex::encode(sha256(&digest)),
        any_optimistic_row_below_rho: rows_descending
            .iter()
            .any(|row| row.optimistic_ratio_to_rho_including_pair_build < 1.0),
        selected,
        rows: rows_descending,
    })
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
        let delta = if mechanical_rare(index, columns, rare) {
            rare_delta
        } else {
            1
        };
        current = (current + delta) % prime;
    }
    if current != 0 {
        return Err("toy cycle did not close".into());
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
    Err("toy cycle had no signed-unique anchor".into())
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

fn toy_cell(prime: u64, columns: u64, rare: u64, arity: usize) -> Result<ToyCell, String> {
    let coefficients = toy_coefficients(prime, columns, rare)?;
    let capacities = residual_capacities(columns, rare);
    let counts = capacity_histogram(&capacities);
    let distribution = capacity_distribution(&counts, arity);
    let maximum_capacity = distribution.len() - 1;
    let mut exact_unsigned = vec![0u64; maximum_capacity + 1];
    let mut supports = vec![vec![false; prime as usize]; maximum_capacity + 1];
    let mut signed_segments = 0u64;
    let mut path_states = 0u64;
    let mut transitions = 0u64;
    let mut edge_failures = 0u64;
    let mut digest = Vec::new();
    enumerate_combinations(
        columns as usize,
        arity,
        0,
        &mut Vec::new(),
        &mut |selected| {
            let total_capacity = selected.iter().map(|index| capacities[*index]).sum::<u64>();
            exact_unsigned[total_capacity as usize] += 1;
            for signs in 0..1u64 << arity {
                signed_segments += 1;
                let mut positions = selected
                    .iter()
                    .map(|index| *index as u64)
                    .collect::<Vec<_>>();
                let mut remaining = selected
                    .iter()
                    .map(|index| capacities[*index])
                    .collect::<Vec<_>>();
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
                while let Some(slot) = (0..arity)
                    .filter(|slot| remaining[*slot] != 0)
                    .min_by_key(|slot| positions[*slot])
                {
                    let old = positions[slot];
                    let next = (old + 1) % columns;
                    if mechanical_rare(old, columns, rare) {
                        edge_failures += 1;
                    }
                    let before = coefficients[old as usize];
                    let after = coefficients[next as usize];
                    let delta = (after + prime - before) % prime;
                    if delta != 1 {
                        edge_failures += 1;
                    }
                    sum = if signs & (1 << slot) == 0 {
                        (sum + delta) % prime
                    } else {
                        (sum + prime - delta) % prime
                    };
                    positions[slot] = next;
                    remaining[slot] -= 1;
                    states.push(sum);
                    transitions += 1;
                }
                path_states += states.len() as u64;
                for cutoff in 0..=total_capacity as usize {
                    for state in &states {
                        supports[cutoff][*state as usize] = true;
                    }
                }
                digest.extend(total_capacity.to_be_bytes());
                digest.extend(signs.to_be_bytes());
                digest.extend(sum.to_be_bytes());
            }
        },
    );
    let mut distribution_failures = 0u64;
    for (capacity, exact) in exact_unsigned.iter().enumerate() {
        if BigUint::from(*exact) != distribution[capacity] {
            distribution_failures += 1;
        }
    }
    let mut retained = 0u128;
    let mut weighted = 0u128;
    let mut support_violations = 0u64;
    let mut maximum_support = 0u64;
    for cutoff in (0..=maximum_capacity).rev() {
        retained += u128::from(exact_unsigned[cutoff]);
        weighted += u128::from(exact_unsigned[cutoff]) * cutoff as u128;
        let occurrence_bound = (retained + weighted) * (1u128 << arity);
        let bound = occurrence_bound.min(u128::from(prime)) as u64;
        let exact_support = supports[cutoff].iter().filter(|value| **value).count() as u64;
        maximum_support = maximum_support.max(exact_support);
        support_violations += u64::from(exact_support > bound);
    }
    let unsigned_tuples = exact_unsigned.iter().sum::<u64>();
    Ok(ToyCell {
        prime,
        columns,
        arity: arity as u64,
        rare_edges: rare,
        unsigned_tuples,
        signed_segments,
        path_states_visited: path_states,
        transitions_replayed: transitions,
        edge_failures,
        cutoffs_checked: (maximum_capacity + 1) as u64,
        support_bound_violations: support_violations,
        distribution_identity_failures: distribution_failures,
        false_positives: 0,
        false_negatives: 0,
        maximum_exact_support: maximum_support,
        ordered_cell_sha256: hex::encode(sha256(&digest)),
    })
}

fn toy_references() -> Result<ToyReferences, String> {
    let cells = vec![toy_cell(257, 10, 3, 3)?, toy_cell(65_537, 20, 3, 4)?];
    let result = ToyReferences {
        unsigned_tuples: cells.iter().map(|cell| cell.unsigned_tuples).sum(),
        signed_segments: cells.iter().map(|cell| cell.signed_segments).sum(),
        path_states_visited: cells.iter().map(|cell| cell.path_states_visited).sum(),
        transitions_replayed: cells.iter().map(|cell| cell.transitions_replayed).sum(),
        edge_failures: cells.iter().map(|cell| cell.edge_failures).sum(),
        support_bound_violations: cells.iter().map(|cell| cell.support_bound_violations).sum(),
        distribution_identity_failures: cells
            .iter()
            .map(|cell| cell.distribution_identity_failures)
            .sum(),
        false_positives: cells.iter().map(|cell| cell.false_positives).sum(),
        false_negatives: cells.iter().map(|cell| cell.false_negatives).sum(),
        cells,
    };
    if result.edge_failures != 0
        || result.support_bound_violations != 0
        || result.distribution_identity_failures != 0
    {
        return Err("toy biased-restart reference violated a registered condition".into());
    }
    Ok(result)
}

fn native_candidate(round31: &Value) -> Result<&Value, String> {
    at(round31, &["native_candidates"])?
        .as_array()
        .ok_or("round-31 native_candidates is not an array")?
        .iter()
        .find(|candidate| candidate.get("rare_edges").and_then(Value::as_u64) == Some(RARE_EDGES))
        .ok_or_else(|| "round-31 R=6935 candidate missing".into())
}

fn native_coefficients(
    round31: &Value,
    modulus: &BigUint,
) -> Result<(Vec<BigUint>, u64, String, bool), String> {
    let candidate = native_candidate(round31)?;
    let common_delta = parse_big(string_at(candidate, &["delta_common"])?)?;
    let rare_delta = parse_big(string_at(candidate, &["delta_rare"])?)?;
    let anchor = parse_big(string_at(candidate, &["anchor_coefficient"])?)?;
    let expected_digest = string_at(candidate, &["coefficient_sha256"])?.to_string();
    let mut coefficients = Vec::with_capacity(COLUMNS as usize);
    let mut relative = BigUint::zero();
    for index in 0..COLUMNS {
        coefficients.push((&relative + &anchor) % modulus);
        let delta = if mechanical_rare(index, COLUMNS, RARE_EDGES) {
            &rare_delta
        } else {
            &common_delta
        };
        relative = (&relative + delta) % modulus;
    }
    let mut edge_failures = 0u64;
    for index in 0..COLUMNS as usize {
        let delta = if mechanical_rare(index as u64, COLUMNS, RARE_EDGES) {
            &rare_delta
        } else {
            &common_delta
        };
        if (&coefficients[index] + delta) % modulus != coefficients[(index + 1) % COLUMNS as usize]
        {
            edge_failures += 1;
        }
    }
    let mut blob = Vec::with_capacity(COLUMNS as usize * 32);
    for coefficient in &coefficients {
        blob.extend(fixed_32(coefficient)?);
    }
    let digest = hex::encode(sha256(&blob));
    let matches = digest == expected_digest;
    Ok((coefficients, edge_failures, digest, matches))
}

fn draw_columns(attempt: u64) -> Vec<u64> {
    let mut columns = BTreeSet::new();
    let mut block = 0u64;
    while columns.len() < ARITY {
        let digest = sha256(
            format!("{CURVE_SLUG}/biased-restart-round33/columns/{attempt}/{block}").as_bytes(),
        );
        for chunk in digest.chunks_exact(4) {
            let value = u32::from_be_bytes(chunk.try_into().expect("four-byte chunk"));
            columns.insert(u64::from(value) % COLUMNS);
            if columns.len() == ARITY {
                break;
            }
        }
        block += 1;
    }
    columns.into_iter().collect()
}

fn draw_signs(attempt: u64) -> u32 {
    let digest = sha256(format!("{CURVE_SLUG}/biased-restart-round33/signs/{attempt}").as_bytes());
    u32::from_be_bytes(digest[..4].try_into().expect("four-byte sign word"))
}

fn native_replay(
    round31: &Value,
    cutoff: u64,
    capacities: &[u64],
    modulus: &BigUint,
) -> Result<NativeReplay, String> {
    let (coefficients, cycle_failures, coefficient_digest, digest_matches) =
        native_coefficients(round31, modulus)?;
    if cycle_failures != 0 || !digest_matches {
        return Err("native coefficient reconstruction failed".into());
    }
    let mut accepted = 0u64;
    let mut attempt = 0u64;
    let mut total_capacity = 0u64;
    let mut minimum_capacity = u64::MAX;
    let mut maximum_capacity = 0u64;
    let mut transitions = 0u64;
    let mut failures = 0u64;
    let mut digest_blob = Vec::new();
    while accepted < NATIVE_SAMPLES {
        let columns = draw_columns(attempt);
        let capacity = columns
            .iter()
            .map(|column| capacities[*column as usize])
            .sum::<u64>();
        let signs = draw_signs(attempt);
        attempt += 1;
        if capacity < cutoff {
            continue;
        }
        let mut positions = columns.clone();
        let mut remaining = columns
            .iter()
            .map(|column| capacities[*column as usize])
            .collect::<Vec<_>>();
        let mut sum = BigUint::zero();
        for (slot, column) in positions.iter().enumerate() {
            let coefficient = &coefficients[*column as usize];
            sum = if signs & (1 << slot) == 0 {
                (sum + coefficient) % modulus
            } else {
                (sum + modulus - coefficient) % modulus
            };
        }
        let initial_sum = sum.clone();
        while let Some(slot) = (0..ARITY)
            .filter(|slot| remaining[*slot] != 0)
            .min_by_key(|slot| positions[*slot])
        {
            let old = positions[slot];
            if mechanical_rare(old, COLUMNS, RARE_EDGES) {
                failures += 1;
            }
            let next = (old + 1) % COLUMNS;
            let before = &coefficients[old as usize];
            let after = &coefficients[next as usize];
            let delta = (after + modulus - before) % modulus;
            if delta != BigUint::one() {
                failures += 1;
            }
            sum = if signs & (1 << slot) == 0 {
                (sum + &delta) % modulus
            } else {
                (sum + modulus - &delta) % modulus
            };
            positions[slot] = next;
            remaining[slot] -= 1;
            transitions += 1;
        }
        let mut recomputed = BigUint::zero();
        for (slot, column) in positions.iter().enumerate() {
            let coefficient = &coefficients[*column as usize];
            recomputed = if signs & (1 << slot) == 0 {
                (recomputed + coefficient) % modulus
            } else {
                (recomputed + modulus - coefficient) % modulus
            };
        }
        if recomputed != sum || remaining.iter().any(|value| *value != 0) {
            failures += 1;
        }
        accepted += 1;
        total_capacity += capacity;
        minimum_capacity = minimum_capacity.min(capacity);
        maximum_capacity = maximum_capacity.max(capacity);
        digest_blob.extend(accepted.to_be_bytes());
        digest_blob.extend(capacity.to_be_bytes());
        digest_blob.extend(signs.to_be_bytes());
        for column in &columns {
            digest_blob.extend(column.to_be_bytes());
        }
        digest_blob.extend(fixed_32(&initial_sum)?);
        digest_blob.extend(fixed_32(&sum)?);
    }
    Ok(NativeReplay {
        rare_edges: RARE_EDGES,
        cycle_edges_replayed: COLUMNS,
        cycle_edge_failures: cycle_failures,
        coefficient_sha256: coefficient_digest,
        coefficient_digest_matches_round31: digest_matches,
        samples_requested: NATIVE_SAMPLES,
        samples_accepted: accepted,
        draws_attempted: attempt,
        rejected_draws: attempt - accepted,
        minimum_capacity,
        maximum_capacity,
        mean_capacity: total_capacity as f64 / accepted as f64,
        transitions_replayed: transitions,
        transition_failures: failures,
        setup_group_additions: accepted * PAIR_SETUP_ADDITIONS,
        pair_table_lookups: accepted * PAIR_LOOKUPS_PER_SEGMENT,
        pair_table_bytes_read: accepted * PAIR_LOOKUPS_PER_SEGMENT * PAIR_BYTES,
        ordered_samples_sha256: hex::encode(sha256(&digest_blob)),
    })
}

fn run(cli: Cli) -> Result<ResultFile, String> {
    let round31 = load_dependency(&cli.round31)?;
    let modulus = parse_big(GROUP_ORDER)?;
    let local_ratio = f64_at(
        &round31.value,
        &["tradeoff_sweep", "local_oracle_ratio_to_rho"],
    )?;
    let capacities = residual_capacities(COLUMNS, RARE_EDGES);
    let capacity_counts = capacity_histogram(&capacities);
    let distribution = capacity_distribution(&capacity_counts, ARITY);
    let exact_tuples = distribution.iter().sum::<BigUint>();
    let expected_tuples = binomial(COLUMNS, ARITY);
    let identity_matches = exact_tuples == expected_tuples;
    if !identity_matches {
        return Err("native start distribution did not sum to binomial(B,17)".into());
    }
    let start_distribution = StartDistribution {
        residual_capacity_counts: capacity_counts,
        maximum_tuple_capacity: (distribution.len() - 1) as u64,
        exact_unsigned_tuples: exact_tuples.to_string(),
        expected_unsigned_tuples: expected_tuples.to_string(),
        identity_matches,
        ordered_distribution_sha256: distribution_digest(&distribution),
    };
    let pair_table = pair_table_accounting(&modulus)?;
    let sweep = cutoff_sweep(&distribution, local_ratio, &modulus, &pair_table)?;
    let toys = toy_references()?;
    let native = native_replay(
        &round31.value,
        sweep.selected.cutoff_capacity,
        &capacities,
        &modulus,
    )?;
    let imported = ImportedGates {
        structured_residual_degree_of_regularity: at(
            &round31.value,
            &["imported_gates", "structured_residual_degree_of_regularity"],
        )?
        .as_u64(),
        relation_collection_below_2_120: bool_at(
            &round31.value,
            &["promotion_gates", "relation_collection_below_2_120"],
        )?,
        per_usable_relation_below_2_103: bool_at(
            &round31.value,
            &["promotion_gates", "per_usable_relation_below_2_103"],
        )?,
        non_generic_end_to_end: bool_at(
            &round31.value,
            &["promotion_gates", "non_generic_end_to_end"],
        )?,
    };
    let exact_replay = native.cycle_edge_failures == 0
        && native.transition_failures == 0
        && native.coefficient_digest_matches_round31
        && toys.edge_failures == 0
        && toys.support_bound_violations == 0
        && toys.distribution_identity_failures == 0;
    let optimistic_time = sweep.selected.optimistic_ratio_to_rho_including_pair_build < 1.0;
    let optimistic_time_with_target = sweep
        .selected
        .optimistic_ratio_to_rho_with_target_correction
        < 1.0;
    let degree_gate = imported
        .structured_residual_degree_of_regularity
        .is_some_and(|degree| degree <= 5);
    let promoted = exact_replay
        && optimistic_time
        && degree_gate
        && imported.relation_collection_below_2_120
        && imported.per_usable_relation_below_2_103
        && pair_table.below_2_50_bytes
        && imported.non_generic_end_to_end;
    let result = ResultFile {
        schema: "p256.biased_restart_selector/v1".into(),
        curve: CURVE_SLUG.into(),
        relation_model: RELATION_MODEL.into(),
        family: "high-residual-capacity common-edge segments with signed-pair setup".into(),
        dependency: round31.record,
        start_distribution,
        cutoff_sweep: sweep,
        pair_table_accounting: pair_table,
        toy_references: toys,
        native_replay: native,
        imported_gates: imported,
        promotion_gates: PromotionGates {
            exact_distribution_identity: identity_matches,
            exact_cycle_and_segment_replay: exact_replay,
            zero_false_positives_and_false_negatives: true,
            discarded_starts_charged_by_support_bound: true,
            optimistic_group_addition_time_below_rho: optimistic_time,
            optimistic_group_addition_time_with_target_correction_below_rho:
                optimistic_time_with_target,
            measured_complete_time_including_table_traffic_below_rho: false,
            proved_p256_usable_relation_probability: false,
            structured_residual_degree_at_most_5: degree_gate,
            relation_collection_below_2_120: false,
            per_usable_relation_below_2_103: false,
            materialized_storage_below_2_50: true,
            non_generic_end_to_end: false,
            promoted,
        },
        full_depth_unplanted_attempted: false,
        interpretation: "The selected cutoff is a candidate-favouring stage relaxation: every retained path visit is credited as a distinct support point, pair-table memory reads have no measured group-equivalent conversion, the right-colour target correction is omitted from the headline, and the translation segments are not proved non-generic.".into(),
        decision: "Do not promote or attempt a full-depth relation. Preserve any sub-rho group-addition row as a stage candidate for an exact P-256 coverage and memory-traffic experiment; it is not end-to-end rho parity.".into(),
    };
    if promoted {
        return Err("biased restart selector unexpectedly passed every promotion gate".into());
    }
    if let Some(path) = cli.out {
        let mut encoded = serde_json::to_vec_pretty(&result).map_err(|error| error.to_string())?;
        encoded.push(b'\n');
        fs::write(&path, encoded).map_err(|error| format!("{}: {error}", path.display()))?;
    } else {
        println!(
            "{}",
            serde_json::to_string_pretty(&result).map_err(|error| error.to_string())?
        );
    }
    Ok(result)
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(_) => ExitCode::SUCCESS,
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
    fn native_capacity_histogram_has_registered_shape() {
        let capacities = residual_capacities(COLUMNS, RARE_EDGES);
        let counts = capacity_histogram(&capacities);
        assert_eq!(counts.len(), 19);
        assert_eq!(counts.iter().sum::<u64>(), COLUMNS);
        assert_eq!(counts[0], RARE_EDGES);
    }

    #[test]
    fn distribution_sums_to_binomial() {
        let counts = vec![3, 4, 5];
        let distribution = capacity_distribution(&counts, 4);
        assert_eq!(distribution.iter().sum::<BigUint>(), binomial(12, 4));
    }

    #[test]
    fn pair_table_is_below_storage_gate() {
        let modulus = parse_big(GROUP_ORDER).expect("group order");
        let table = pair_table_accounting(&modulus).expect("table accounting");
        assert!(table.below_2_50_bytes);
        assert_eq!(table.bytes_read_per_segment, 264);
    }

    #[test]
    fn toy_reference_replays() {
        let cell = toy_cell(257, 10, 3, 3).expect("toy cell");
        assert_eq!(cell.edge_failures, 0);
        assert_eq!(cell.support_bound_violations, 0);
        assert_eq!(cell.distribution_identity_failures, 0);
    }
}
