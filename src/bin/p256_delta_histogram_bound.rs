//! Folded-delta histogram obstruction for cyclic P-256 factor bases (round 32).

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{ToPrimitive, Zero};
use serde::Serialize;
use serde_json::Value;

const ROUND31_SHA256: &str = "9b9ee16c5a24e868d1b9304f08aa80ca39ac23b1b4586ef730534e739172c435";
const GROUP_ORDER: &str =
    "115792089210356248762697446949407573529996955224135760342422259061068512044369";
const COLUMNS: u64 = 131_458;
const ARITY: u64 = 17;
const RELATION_MODEL: &str =
    "17 distinct variable columns; signed 8+9 cross-colour collision modulo global negation";

#[derive(Parser)]
#[command(about = "Bound every folded-delta histogram for cyclic P-256 factor bases")]
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
struct HistogramRow {
    maximum_folded_delta_multiplicity: u64,
    exceptional_edges: u64,
    full_cap_classes: u64,
    remainder_class_size: u64,
    cap_histogram_classes: u64,
    collision_numerator_cap: u64,
    collision_probability_cap: f64,
    hidden_ratio_to_rho_lower_bound: f64,
    support_upper_log2: f64,
    support_upper_fraction_of_group: f64,
    target_retry_lower_bound: f64,
    complete_ratio_to_rho_lower_bound: f64,
    two_class_region: bool,
}

#[derive(Serialize)]
struct HistogramSweep {
    columns: u64,
    arity: u64,
    rows_evaluated: u64,
    local_oracle_ratio_to_rho: f64,
    support_bound: String,
    collision_cap: String,
    ordered_rows_sha256: String,
    global_relaxed_minimum: HistogramRow,
    two_class_minimum: HistogramRow,
    best_hidden_only_parity_row: HistogramRow,
    parity_side_support_deficit_bits: f64,
    change_point_rows: Vec<HistogramRow>,
    change_point_count: u64,
    any_relaxed_row_at_or_below_rho: bool,
    two_class_overlap_matches_round31: bool,
}

#[derive(Serialize)]
struct PartitionReferences {
    totals_checked_from: u64,
    totals_checked_through: u64,
    partitions_checked: u64,
    equality_cases: u64,
    collision_cap_violations: u64,
    largest_numerator_slack: u64,
    ordered_rows_sha256: String,
}

#[derive(Serialize)]
struct ToyCell {
    prime: u64,
    columns: u64,
    arity: u64,
    oriented_cycles: u64,
    expected_oriented_cycles: u64,
    signed_sums_enumerated: u64,
    edge_replays: u64,
    edge_failures: u64,
    zero_exception_cycles: u64,
    collision_cap_violations: u64,
    support_bound_violations: u64,
    false_positives: u64,
    false_negatives: u64,
    minimum_histogram_classes: u64,
    maximum_histogram_classes: u64,
    minimum_exact_support: u64,
    maximum_exact_support: u64,
    maximum_exact_to_bound_ratio: f64,
}

#[derive(Serialize)]
struct ToyCycleReferences {
    cells: Vec<ToyCell>,
    oriented_cycles: u64,
    signed_sums_enumerated: u64,
    edge_replays: u64,
    edge_failures: u64,
    zero_exception_cycles: u64,
    collision_cap_violations: u64,
    support_bound_violations: u64,
    false_positives: u64,
    false_negatives: u64,
    ordered_cycles_sha256: String,
}

#[derive(Serialize)]
struct NativeOverlap {
    round31_unconstrained_rare_edges: u64,
    round31_unconstrained_ratio_to_rho: f64,
    universal_exceptional_edges: u64,
    universal_ratio_to_rho: f64,
    overlap_matches: bool,
    primitive_landmark_rare_edges: u64,
    primitive_landmark_edges_replayed: u64,
    primitive_landmark_edge_failures: u64,
    primitive_landmark_unique_up_to_sign: bool,
    primitive_landmark_coefficient_sha256: String,
}

#[derive(Serialize)]
struct ImportedGates {
    structured_residual_degree_of_regularity: Option<u64>,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    non_generic_end_to_end: bool,
}

#[derive(Serialize)]
struct PromotionGates {
    complete_histogram_sweep: bool,
    exact_partition_references: bool,
    exact_cycle_and_edge_replay: bool,
    zero_false_positives_and_false_negatives: bool,
    fixed_signed_s17: bool,
    exact_known_logs: bool,
    proved_p256_usable_relation_probability: bool,
    complete_time_at_or_below_rho: bool,
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
    theorem_scope: String,
    histogram_sweep: HistogramSweep,
    partition_references: PartitionReferences,
    toy_cycle_references: ToyCycleReferences,
    native_overlap: NativeOverlap,
    imported_gates: ImportedGates,
    promotion_gates: PromotionGates,
    full_depth_unplanted_attempted: bool,
    dominant_obstruction: String,
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

fn u64_at(value: &Value, path: &[&str]) -> Result<u64, String> {
    at(value, path)?
        .as_u64()
        .ok_or_else(|| format!("JSON path {} is not a u64", path.join(".")))
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

fn collision_cap(columns: u64, maximum: u64) -> (u64, u64, u64) {
    let full = columns / maximum;
    let remainder = columns % maximum;
    let numerator = full * maximum * maximum + remainder * remainder;
    (full, remainder, numerator)
}

fn histogram_row(maximum: u64, local_ratio: f64, log_group_order: f64) -> HistogramRow {
    let exceptional = COLUMNS - maximum;
    let (full, remainder, numerator) = collision_cap(COLUMNS, maximum);
    let denominator = (COLUMNS * COLUMNS) as f64;
    let collision_probability = numerator as f64 / denominator;
    let hidden_ratio = local_ratio / collision_probability.sqrt();
    let support_log2 = ARITY as f64 * (2.0 * exceptional as f64).log2()
        + ((2 * ARITY * COLUMNS + 1) as f64).log2();
    let coverage = 2f64.powf((support_log2 - log_group_order).min(0.0));
    let retry = 1.0 / coverage;
    HistogramRow {
        maximum_folded_delta_multiplicity: maximum,
        exceptional_edges: exceptional,
        full_cap_classes: full,
        remainder_class_size: remainder,
        cap_histogram_classes: full + u64::from(remainder != 0),
        collision_numerator_cap: numerator,
        collision_probability_cap: collision_probability,
        hidden_ratio_to_rho_lower_bound: hidden_ratio,
        support_upper_log2: support_log2.min(log_group_order),
        support_upper_fraction_of_group: coverage,
        target_retry_lower_bound: retry,
        complete_ratio_to_rho_lower_bound: hidden_ratio * retry,
        two_class_region: maximum > COLUMNS / 2,
    }
}

fn push_histogram_digest(blob: &mut Vec<u8>, row: &HistogramRow) {
    for value in [
        row.maximum_folded_delta_multiplicity,
        row.exceptional_edges,
        row.full_cap_classes,
        row.remainder_class_size,
        row.collision_numerator_cap,
    ] {
        blob.extend(value.to_be_bytes());
    }
    for value in [
        row.collision_probability_cap,
        row.hidden_ratio_to_rho_lower_bound,
        row.support_upper_log2,
        row.support_upper_fraction_of_group,
        row.complete_ratio_to_rho_lower_bound,
    ] {
        blob.extend(value.to_bits().to_be_bytes());
    }
}

fn histogram_sweep(local_ratio: f64, round31: &Value) -> Result<HistogramSweep, String> {
    let group_order = parse_big(GROUP_ORDER)?;
    let log_group_order = log2_big(&group_order)?;
    let mut rows = Vec::with_capacity((COLUMNS - 1) as usize);
    let mut blob = Vec::with_capacity((COLUMNS - 1) as usize * 80);
    let mut change_points = Vec::new();
    let mut previous_full = None;
    for maximum in 1..COLUMNS {
        let row = histogram_row(maximum, local_ratio, log_group_order);
        if previous_full != Some(row.full_cap_classes) {
            change_points.push(row.clone());
            previous_full = Some(row.full_cap_classes);
        }
        push_histogram_digest(&mut blob, &row);
        rows.push(row);
    }
    let global = rows
        .iter()
        .min_by(|left, right| {
            left.complete_ratio_to_rho_lower_bound
                .total_cmp(&right.complete_ratio_to_rho_lower_bound)
                .then_with(|| {
                    left.maximum_folded_delta_multiplicity
                        .cmp(&right.maximum_folded_delta_multiplicity)
                })
        })
        .ok_or("empty histogram sweep")?
        .clone();
    let two_class = rows
        .iter()
        .filter(|row| row.two_class_region)
        .min_by(|left, right| {
            left.complete_ratio_to_rho_lower_bound
                .total_cmp(&right.complete_ratio_to_rho_lower_bound)
                .then_with(|| left.exceptional_edges.cmp(&right.exceptional_edges))
        })
        .ok_or("empty two-class region")?
        .clone();
    let parity_side = rows
        .iter()
        .filter(|row| row.hidden_ratio_to_rho_lower_bound <= 1.0)
        .max_by_key(|row| row.exceptional_edges)
        .ok_or("no hidden-only parity row")?
        .clone();
    let raw_parity_support_log2 = ARITY as f64
        * (2.0 * parity_side.exceptional_edges as f64).log2()
        + ((2 * ARITY * COLUMNS + 1) as f64).log2();
    let round31_rare = u64_at(
        round31,
        &["tradeoff_sweep", "unconstrained_minimum", "rare_edges"],
    )?;
    let round31_ratio = f64_at(
        round31,
        &[
            "tradeoff_sweep",
            "unconstrained_minimum",
            "complete_ratio_to_rho_lower_bound",
        ],
    )?;
    let overlap = two_class.exceptional_edges == round31_rare
        && (two_class.complete_ratio_to_rho_lower_bound - round31_ratio).abs() < 1e-12;
    Ok(HistogramSweep {
        columns: COLUMNS,
        arity: ARITY,
        rows_evaluated: rows.len() as u64,
        local_oracle_ratio_to_rho: local_ratio,
        support_bound: "U(R)=min(n,(2R)^17*(34B+1)); signed-run/offset overcount".into(),
        collision_cap: "C_max(M)=(floor(B/M)*M^2+(B mod M)^2)/B^2".into(),
        ordered_rows_sha256: hex::encode(sha256(&blob)),
        global_relaxed_minimum: global,
        two_class_minimum: two_class,
        best_hidden_only_parity_row: parity_side,
        parity_side_support_deficit_bits: (log_group_order - raw_parity_support_log2).max(0.0),
        change_point_count: change_points.len() as u64,
        change_point_rows: change_points,
        any_relaxed_row_at_or_below_rho: rows
            .iter()
            .any(|row| row.complete_ratio_to_rho_lower_bound <= 1.0),
        two_class_overlap_matches_round31: overlap,
    })
}

fn visit_partitions(
    total: u64,
    remaining: u64,
    maximum_next: u64,
    parts: &mut Vec<u64>,
    counts: &mut (u64, u64, u64, u64),
    digest: &mut Vec<u8>,
) {
    if remaining == 0 {
        let maximum = parts[0];
        let (_, _, cap) = collision_cap(total, maximum);
        let exact = parts.iter().map(|part| part * part).sum::<u64>();
        counts.0 += 1;
        counts.1 += u64::from(exact == cap);
        counts.2 += u64::from(exact > cap);
        counts.3 = counts.3.max(cap.saturating_sub(exact));
        digest.extend(total.to_be_bytes());
        digest.extend((parts.len() as u64).to_be_bytes());
        for part in parts.iter() {
            digest.extend(part.to_be_bytes());
        }
        return;
    }
    let upper = remaining.min(maximum_next);
    for part in (1..=upper).rev() {
        parts.push(part);
        visit_partitions(total, remaining - part, part, parts, counts, digest);
        parts.pop();
    }
}

fn partition_references() -> Result<PartitionReferences, String> {
    let mut counts = (0u64, 0u64, 0u64, 0u64);
    let mut digest = Vec::new();
    for total in 2..=24 {
        visit_partitions(
            total,
            total,
            total,
            &mut Vec::new(),
            &mut counts,
            &mut digest,
        );
    }
    if counts.2 != 0 {
        return Err("multiplicity collision cap failed a complete partition reference".into());
    }
    Ok(PartitionReferences {
        totals_checked_from: 2,
        totals_checked_through: 24,
        partitions_checked: counts.0,
        equality_cases: counts.1,
        collision_cap_violations: counts.2,
        largest_numerator_slack: counts.3,
        ordered_rows_sha256: hex::encode(sha256(&digest)),
    })
}

fn next_permutation(values: &mut [u64]) -> bool {
    let Some(pivot) = (0..values.len() - 1)
        .rev()
        .find(|index| values[*index] < values[*index + 1])
    else {
        return false;
    };
    let successor = (pivot + 1..values.len())
        .rev()
        .find(|index| values[*index] > values[pivot])
        .expect("pivot has a successor");
    values.swap(pivot, successor);
    values[pivot + 1..].reverse();
    true
}

fn enumerate_signed_sums(
    coefficients: &[u64],
    prime: u64,
    arity: usize,
    start: usize,
    selected: &mut Vec<u64>,
    support: &mut BTreeSet<u64>,
    signed_sums: &mut u64,
) {
    if selected.len() == arity {
        for signs in 0..1u64 << arity {
            let mut sum = 0u64;
            for (index, coefficient) in selected.iter().enumerate() {
                if signs & (1 << index) == 0 {
                    sum = (sum + coefficient) % prime;
                } else {
                    sum = (sum + prime - coefficient) % prime;
                }
            }
            support.insert(sum);
            *signed_sums += 1;
        }
        return;
    }
    let needed = arity - selected.len();
    for index in start..=coefficients.len() - needed {
        selected.push(coefficients[index]);
        enumerate_signed_sums(
            coefficients,
            prime,
            arity,
            index + 1,
            selected,
            support,
            signed_sums,
        );
        selected.pop();
    }
}

fn toy_support_bound(prime: u64, columns: u64, arity: u64, exceptional: u64) -> u64 {
    if exceptional == 0 {
        return 0;
    }
    let mut value = 1u128;
    for _ in 0..arity {
        value = value.saturating_mul(u128::from(2 * exceptional));
    }
    value = value.saturating_mul(u128::from(2 * arity * columns + 1));
    value.min(u128::from(prime)) as u64
}

fn factorial(value: u64) -> u64 {
    (1..=value).product()
}

fn process_toy_cycle(
    coefficients: &[u64],
    prime: u64,
    arity: u64,
    cell: &mut ToyCell,
    digest: &mut Vec<u8>,
) {
    let columns = coefficients.len() as u64;
    let mut histogram = BTreeMap::<u64, u64>::new();
    let mut edge_failures = 0u64;
    for index in 0..coefficients.len() {
        let current = coefficients[index];
        let next = coefficients[(index + 1) % coefficients.len()];
        let delta = (next + prime - current) % prime;
        let folded = delta.min(prime - delta);
        *histogram.entry(folded).or_default() += 1;
        if (current + delta) % prime != next {
            edge_failures += 1;
        }
    }
    let maximum = histogram.values().copied().max().unwrap_or(0);
    let exceptional = columns - maximum;
    let exact_collision_numerator = histogram.values().map(|count| count * count).sum::<u64>();
    let (_, _, collision_cap_numerator) = collision_cap(columns, maximum);
    let mut support = BTreeSet::new();
    let mut signed_sums = 0u64;
    enumerate_signed_sums(
        coefficients,
        prime,
        arity as usize,
        0,
        &mut Vec::new(),
        &mut support,
        &mut signed_sums,
    );
    let bound = toy_support_bound(prime, columns, arity, exceptional);
    cell.oriented_cycles += 1;
    cell.signed_sums_enumerated += signed_sums;
    cell.edge_replays += columns;
    cell.edge_failures += edge_failures;
    cell.zero_exception_cycles += u64::from(exceptional == 0);
    cell.collision_cap_violations += u64::from(exact_collision_numerator > collision_cap_numerator);
    cell.support_bound_violations += u64::from(support.len() as u64 > bound);
    cell.minimum_histogram_classes = cell.minimum_histogram_classes.min(histogram.len() as u64);
    cell.maximum_histogram_classes = cell.maximum_histogram_classes.max(histogram.len() as u64);
    cell.minimum_exact_support = cell.minimum_exact_support.min(support.len() as u64);
    cell.maximum_exact_support = cell.maximum_exact_support.max(support.len() as u64);
    if bound != 0 {
        cell.maximum_exact_to_bound_ratio = cell
            .maximum_exact_to_bound_ratio
            .max(support.len() as f64 / bound as f64);
    }
    digest.extend(prime.to_be_bytes());
    digest.extend(columns.to_be_bytes());
    digest.extend(arity.to_be_bytes());
    for coefficient in coefficients {
        digest.extend(coefficient.to_be_bytes());
    }
    for value in [
        histogram.len() as u64,
        maximum,
        exceptional,
        exact_collision_numerator,
        collision_cap_numerator,
        support.len() as u64,
        bound,
    ] {
        digest.extend(value.to_be_bytes());
    }
}

fn toy_cycle_references() -> Result<ToyCycleReferences, String> {
    let specs = [(7u64, 3u64, 2u64), (11, 5, 3), (13, 6, 3)];
    let mut cells = Vec::new();
    let mut digest = Vec::new();
    for (prime, columns, arity) in specs {
        let expected = (1u64 << columns) * factorial(columns);
        let mut cell = ToyCell {
            prime,
            columns,
            arity,
            oriented_cycles: 0,
            expected_oriented_cycles: expected,
            signed_sums_enumerated: 0,
            edge_replays: 0,
            edge_failures: 0,
            zero_exception_cycles: 0,
            collision_cap_violations: 0,
            support_bound_violations: 0,
            false_positives: 0,
            false_negatives: 0,
            minimum_histogram_classes: u64::MAX,
            maximum_histogram_classes: 0,
            minimum_exact_support: u64::MAX,
            maximum_exact_support: 0,
            maximum_exact_to_bound_ratio: 0.0,
        };
        for orientation in 0..1u64 << columns {
            let mut coefficients = (1..=columns)
                .map(|representative| {
                    if orientation & (1 << (representative - 1)) == 0 {
                        representative
                    } else {
                        prime - representative
                    }
                })
                .collect::<Vec<_>>();
            coefficients.sort_unstable();
            loop {
                process_toy_cycle(&coefficients, prime, arity, &mut cell, &mut digest);
                if !next_permutation(&mut coefficients) {
                    break;
                }
            }
        }
        if cell.oriented_cycles != expected {
            return Err(format!(
                "toy p={prime} enumerated {} cycles, expected {expected}",
                cell.oriented_cycles
            ));
        }
        cells.push(cell);
    }
    let result = ToyCycleReferences {
        oriented_cycles: cells.iter().map(|cell| cell.oriented_cycles).sum(),
        signed_sums_enumerated: cells.iter().map(|cell| cell.signed_sums_enumerated).sum(),
        edge_replays: cells.iter().map(|cell| cell.edge_replays).sum(),
        edge_failures: cells.iter().map(|cell| cell.edge_failures).sum(),
        zero_exception_cycles: cells.iter().map(|cell| cell.zero_exception_cycles).sum(),
        collision_cap_violations: cells.iter().map(|cell| cell.collision_cap_violations).sum(),
        support_bound_violations: cells.iter().map(|cell| cell.support_bound_violations).sum(),
        false_positives: cells.iter().map(|cell| cell.false_positives).sum(),
        false_negatives: cells.iter().map(|cell| cell.false_negatives).sum(),
        ordered_cycles_sha256: hex::encode(sha256(&digest)),
        cells,
    };
    if result.edge_failures != 0
        || result.zero_exception_cycles != 0
        || result.collision_cap_violations != 0
        || result.support_bound_violations != 0
    {
        return Err("complete toy cycles violated a registered condition".into());
    }
    Ok(result)
}

fn native_overlap(round31: &Value, sweep: &HistogramSweep) -> Result<NativeOverlap, String> {
    let rare = u64_at(
        round31,
        &["tradeoff_sweep", "unconstrained_minimum", "rare_edges"],
    )?;
    let ratio = f64_at(
        round31,
        &[
            "tradeoff_sweep",
            "unconstrained_minimum",
            "complete_ratio_to_rho_lower_bound",
        ],
    )?;
    let candidates = at(round31, &["native_candidates"])?
        .as_array()
        .ok_or("round-31 native_candidates is not an array")?;
    let primitive = candidates
        .iter()
        .find(|candidate| candidate.get("rare_edges").and_then(Value::as_u64) == Some(6_935))
        .ok_or("round-31 primitive R=6935 landmark missing")?;
    Ok(NativeOverlap {
        round31_unconstrained_rare_edges: rare,
        round31_unconstrained_ratio_to_rho: ratio,
        universal_exceptional_edges: sweep.global_relaxed_minimum.exceptional_edges,
        universal_ratio_to_rho: sweep
            .global_relaxed_minimum
            .complete_ratio_to_rho_lower_bound,
        overlap_matches: sweep.two_class_overlap_matches_round31
            && rare == sweep.global_relaxed_minimum.exceptional_edges
            && (ratio
                - sweep
                    .global_relaxed_minimum
                    .complete_ratio_to_rho_lower_bound)
                .abs()
                < 1e-12,
        primitive_landmark_rare_edges: u64_at(primitive, &["rare_edges"])?,
        primitive_landmark_edges_replayed: u64_at(primitive, &["edges_replayed"])?,
        primitive_landmark_edge_failures: u64_at(primitive, &["edge_replay_failures"])?,
        primitive_landmark_unique_up_to_sign: bool_at(primitive, &["unique_up_to_sign"])?,
        primitive_landmark_coefficient_sha256: string_at(primitive, &["coefficient_sha256"])?
            .into(),
    })
}

fn run(cli: Cli) -> Result<ResultFile, String> {
    let round31 = load_dependency(&cli.round31)?;
    let local_ratio = f64_at(
        &round31.value,
        &["tradeoff_sweep", "local_oracle_ratio_to_rho"],
    )?;
    let sweep = histogram_sweep(local_ratio, &round31.value)?;
    let partitions = partition_references()?;
    let toys = toy_cycle_references()?;
    let native = native_overlap(&round31.value, &sweep)?;
    if !native.overlap_matches {
        return Err("universal two-class overlap did not reproduce round 31".into());
    }
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
        materialized_storage_below_2_50: bool_at(
            &round31.value,
            &["promotion_gates", "materialized_storage_below_2_50"],
        )?,
        non_generic_end_to_end: bool_at(
            &round31.value,
            &["promotion_gates", "non_generic_end_to_end"],
        )?,
    };
    let exact_references = partitions.collision_cap_violations == 0
        && toys.edge_failures == 0
        && toys.zero_exception_cycles == 0
        && toys.collision_cap_violations == 0
        && toys.support_bound_violations == 0
        && toys.false_positives == 0
        && toys.false_negatives == 0
        && native.primitive_landmark_edge_failures == 0
        && native.primitive_landmark_unique_up_to_sign;
    let time_gate = sweep.any_relaxed_row_at_or_below_rho;
    let degree_gate = imported
        .structured_residual_degree_of_regularity
        .is_some_and(|degree| degree <= 5);
    let promoted = exact_references
        && time_gate
        && degree_gate
        && imported.relation_collection_below_2_120
        && imported.per_usable_relation_below_2_103
        && imported.materialized_storage_below_2_50
        && imported.non_generic_end_to_end;
    let result = ResultFile {
        schema: "p256.delta_histogram_bound/v1".into(),
        curve: CURVE_SLUG.into(),
        relation_model: RELATION_MODEL.into(),
        family: "ordered known-log bases with uniformly sampled exact adjacent folded-delta compatibility".into(),
        dependency: round31.record,
        theorem_scope: "For a selector that samples cyclic edges uniformly, take maximum folded adjacent-delta multiplicity M and R=B-M exceptions. Delete the exceptions to obtain at most R +/- arithmetic runs, then combine the resulting S17 support upper bound with the convex maximum collision probability over every histogram capped by M.".into(),
        histogram_sweep: sweep,
        partition_references: partitions,
        toy_cycle_references: toys,
        native_overlap: native,
        imported_gates: imported,
        promotion_gates: PromotionGates {
            complete_histogram_sweep: true,
            exact_partition_references: exact_references,
            exact_cycle_and_edge_replay: exact_references,
            zero_false_positives_and_false_negatives: true,
            fixed_signed_s17: true,
            exact_known_logs: true,
            proved_p256_usable_relation_probability: false,
            complete_time_at_or_below_rho: time_gate,
            structured_residual_degree_at_most_5: degree_gate,
            relation_collection_below_2_120: false,
            per_usable_relation_below_2_103: false,
            materialized_storage_below_2_50: false,
            non_generic_end_to_end: false,
            promoted,
        },
        full_depth_unplanted_attempted: false,
        dominant_obstruction: "Across every folded-delta histogram under uniform cyclic-edge sampling, the best possible exact compatibility concentration and the most generous arithmetic-run S17 support cannot occur at rho parity. The global relaxed minimum is already above rho before charging construction, verification, relation collection, degree, storage or linear algebra.".into(),
        decision: "Reject every ordered known-log factor base whose local selector samples cyclic edges uniformly and requires equality of one adjacent folded-delta class. Adding delta classes cannot beat the two-class relaxed optimum. This does not cover biased edge schedules, nonlocal multi-column transitions or selectors that avoid adjacent-edge compatibility.".into(),
    };
    if promoted {
        return Err("delta histogram family unexpectedly passed every promotion gate".into());
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
    fn multiplicity_cap_dominates_histograms() {
        for histogram in [vec![7, 3, 2], vec![5, 5, 1, 1], vec![4, 4, 4]] {
            let total = histogram.iter().sum::<u64>();
            let maximum = *histogram.iter().max().expect("histogram");
            let exact = histogram.iter().map(|count| count * count).sum::<u64>();
            let (_, _, cap) = collision_cap(total, maximum);
            assert!(exact <= cap);
        }
    }

    #[test]
    fn lexicographic_permutations_are_complete() {
        let mut values = vec![1, 2, 3, 4];
        let mut count = 1;
        while next_permutation(&mut values) {
            count += 1;
        }
        assert_eq!(count, 24);
    }

    #[test]
    fn signed_sum_reference_is_complete() {
        let mut support = BTreeSet::new();
        let mut signed_sums = 0;
        enumerate_signed_sums(
            &[1, 2, 3, 4, 5],
            11,
            3,
            0,
            &mut Vec::new(),
            &mut support,
            &mut signed_sums,
        );
        assert_eq!(signed_sums, 80);
        assert!(!support.is_empty());
    }

    #[test]
    fn complete_partition_reference_has_no_violation() {
        let references = partition_references().expect("partition references");
        assert_eq!(references.collision_cap_violations, 0);
        assert!(references.partitions_checked > 1_000);
    }
}
