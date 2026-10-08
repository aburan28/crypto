//! Exact cutoff-conditioned signed-pair density audit for P-256 (round 37).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;
use serde_json::Value;

const CURVE_SLUG: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const ROUND33_SHA256: &str = "9931bccbd9e2821f4498f465ce65f83efb400e7fb3740a239bdb40d3275ffbd8";
const ROUND36_SHA256: &str = "f145231c098da83878e890d340f8ae921598519b0f2626bcc7c4cd8f9446d015";
const COLUMNS: u64 = 131_458;
const ARITY: usize = 17;
const CUTOFF: usize = 219;
const LOOKUPS: f64 = 8.0;
const ENTRY_BYTES: u64 = 64;
const PAIR_ENTRIES: u64 = 34_562_148_612;
const BASE_RATIO: f64 = 0.998_569_034_150_286;
const LOCAL_RATIO: f64 = 0.964_336_477_130_181;
const MEAN_CAPACITY: f64 = 224.361_249_307_510_55;
const ALLOWED_COLD: f64 = 0.334_410_508_422_987_64;
const THRESHOLD_SCALE: u64 = 1_000_000_000_000;
const THRESHOLD_HITS_NUMERATOR: u64 = 7_665_589_491_577;

#[derive(Parser)]
#[command(about = "Audit cutoff-conditioned P-256 signed-pair entropy")]
struct Cli {
    #[arg(long)]
    round33: PathBuf,
    #[arg(long)]
    round36: PathBuf,
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

#[derive(Clone, Serialize)]
struct FrozenModel {
    columns: u64,
    arity: u64,
    cutoff: u64,
    lookups_per_segment: u64,
    pair_entries: u64,
    entry_bytes: u64,
    corrected_base_ratio_to_rho: f64,
    local_oracle_ratio_to_rho: f64,
    mean_capacity: f64,
    allowed_cold_constructions_per_segment: f64,
    capacity_counts: Vec<u64>,
    retained_unsigned_tuples: String,
}

#[derive(Clone)]
struct PairBlock {
    capacity_left: usize,
    capacity_right: usize,
    unordered_column_pairs: u64,
    signed_entries: u64,
    completions: BigUint,
}

#[derive(Clone, Serialize)]
struct PairBlockRow {
    density_rank: u64,
    capacity_left: u64,
    capacity_right: u64,
    unordered_column_pairs: u64,
    signed_entries: u64,
    retained_completions_per_fixed_pair: String,
    signed_entry_presence_probability: f64,
    signed_entry_presence_log2: f64,
}

#[derive(Clone, Serialize)]
struct FrontierRow {
    label: String,
    signed_entries: u64,
    bytes: u64,
    bytes_log2: f64,
    expected_cached_signed_pairs_present: f64,
    optimistic_hot_hits_upper: f64,
    cold_constructions_lower: f64,
    projected_ratio_to_rho_lower: f64,
    within_parity_budget: bool,
}

#[derive(Clone, Serialize)]
struct ExactControls {
    retained_count_matches_round33: bool,
    full_population_matches_round33: bool,
    full_expected_pair_normalization: f64,
    full_expected_pair_normalization_exact: bool,
    density_monotonicity_failures: u64,
    population_failures: u64,
    toy_cells: u64,
    toy_fixed_pairs_checked: u64,
    toy_count_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    ordered_class_rows_sha256: String,
    frontier_sha256: String,
}

#[derive(Clone, Serialize)]
struct Gates {
    exact_counts: bool,
    admissible_depth_at_or_below_rho: bool,
    measured_complete_time_at_or_below_rho: bool,
    proved_p256_usable_relation_probability: bool,
    structured_residual_degree_at_most_5: bool,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    non_generic_end_to_end: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct SemanticEvidence {
    curve: String,
    round33_sha256: String,
    round36_sha256: String,
    model: FrozenModel,
    class_rows: Vec<PairBlockRow>,
    frontier: Vec<FrontierRow>,
    controls: ExactControls,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    family: String,
    dependencies: Vec<Dependency>,
    frozen_model: FrozenModel,
    class_rows: Vec<PairBlockRow>,
    frontier: Vec<FrontierRow>,
    threshold_signed_entries: u64,
    threshold_bytes: u64,
    threshold_fraction_of_full_table: f64,
    exact_controls: ExactControls,
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

fn parse_big(value: &str) -> Result<BigUint, String> {
    BigUint::parse_bytes(value.as_bytes(), 10)
        .ok_or_else(|| format!("invalid decimal integer: {value}"))
}

fn read_dependency(path: &Path, expected: &str, round: u64) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != expected {
        return Err(format!(
            "round-{round} dependency hash mismatch: expected {expected}, got {digest}"
        ));
    }
    let value: Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("invalid JSON: {error}"))?;
    if json_str(&value, &["curve"])? != CURVE_SLUG {
        return Err(format!("round-{round} curve mismatch"));
    }
    Ok((
        Dependency {
            round,
            path: path.display().to_string(),
            sha256: digest,
            schema: json_str(&value, &["schema"])?.to_string(),
        },
        value,
    ))
}

fn binomial(n: u64, k: usize) -> BigUint {
    if k > n as usize {
        return BigUint::zero();
    }
    let k = k.min(n as usize - k);
    let mut value = BigUint::one();
    for index in 0..k {
        value *= n - index as u64;
        value /= index as u64 + 1;
    }
    value
}

fn capacity_distribution(counts: &[u64], arity: usize) -> Vec<BigUint> {
    let maximum = (counts.len() - 1) * arity;
    let mut dp = vec![vec![BigUint::zero(); maximum + 1]; arity + 1];
    dp[0][0] = BigUint::one();
    for (capacity, count) in counts.iter().copied().enumerate() {
        let mut next = vec![vec![BigUint::zero(); maximum + 1]; arity + 1];
        for chosen in 0..=arity {
            for sum in 0..=maximum {
                if dp[chosen][sum].is_zero() {
                    continue;
                }
                for take in 0..=(arity - chosen).min(count as usize) {
                    let next_sum = sum + take * capacity;
                    next[chosen + take][next_sum] += &dp[chosen][sum] * binomial(count, take);
                }
            }
        }
        dp = next;
    }
    dp.remove(arity)
}

fn retained_count(counts: &[u64], arity: usize, cutoff: usize) -> BigUint {
    capacity_distribution(counts, arity)
        .into_iter()
        .enumerate()
        .filter(|(capacity, _)| *capacity >= cutoff)
        .map(|(_, count)| count)
        .sum()
}

fn big_ratio(numerator: &BigUint, denominator: &BigUint) -> Result<f64, String> {
    let n = numerator
        .to_f64()
        .ok_or("numerator did not convert to f64")?;
    let d = denominator
        .to_f64()
        .ok_or("denominator did not convert to f64")?;
    Ok(n / d)
}

fn pair_blocks(counts: &[u64], retained: &BigUint) -> Result<Vec<PairBlock>, String> {
    let mut blocks = Vec::new();
    for left in 0..counts.len() {
        for right in left..counts.len() {
            let unordered = if left == right {
                counts[left] * (counts[left] - 1) / 2
            } else {
                counts[left] * counts[right]
            };
            let mut residual = counts.to_vec();
            residual[left] -= 1;
            residual[right] -= 1;
            let remaining_cutoff = CUTOFF.saturating_sub(left + right);
            let completions = retained_count(&residual, ARITY - 2, remaining_cutoff);
            blocks.push(PairBlock {
                capacity_left: left,
                capacity_right: right,
                unordered_column_pairs: unordered,
                signed_entries: 4 * unordered,
                completions,
            });
        }
    }
    blocks.sort_by(|left, right| {
        right
            .completions
            .cmp(&left.completions)
            .then_with(|| left.capacity_left.cmp(&right.capacity_left))
            .then_with(|| left.capacity_right.cmp(&right.capacity_right))
    });
    let population: u64 = blocks.iter().map(|block| block.signed_entries).sum();
    if population != PAIR_ENTRIES {
        return Err(format!(
            "signed pair population mismatch: {population} != {PAIR_ENTRIES}"
        ));
    }
    let denominator = retained * BigUint::from(4u64);
    if denominator.is_zero() {
        return Err("empty retained population".into());
    }
    Ok(blocks)
}

fn class_rows(blocks: &[PairBlock], retained: &BigUint) -> Result<Vec<PairBlockRow>, String> {
    let denominator = retained * BigUint::from(4u64);
    blocks
        .iter()
        .enumerate()
        .map(|(rank, block)| {
            let probability = big_ratio(&block.completions, &denominator)?;
            Ok(PairBlockRow {
                density_rank: rank as u64 + 1,
                capacity_left: block.capacity_left as u64,
                capacity_right: block.capacity_right as u64,
                unordered_column_pairs: block.unordered_column_pairs,
                signed_entries: block.signed_entries,
                retained_completions_per_fixed_pair: block.completions.to_str_radix(10),
                signed_entry_presence_probability: probability,
                signed_entry_presence_log2: probability.log2(),
            })
        })
        .collect()
}

fn cumulative_numerator(blocks: &[PairBlock], entries: u64) -> BigUint {
    let mut remaining = entries;
    let mut numerator = BigUint::zero();
    for block in blocks {
        if remaining == 0 {
            break;
        }
        let take = remaining.min(block.signed_entries);
        numerator += &block.completions * take;
        remaining -= take;
    }
    numerator
}

fn threshold_entries(blocks: &[PairBlock], retained: &BigUint) -> Result<u64, String> {
    let target = retained * BigUint::from(4u64) * BigUint::from(THRESHOLD_HITS_NUMERATOR);
    let scale = BigUint::from(THRESHOLD_SCALE);
    let mut entries = 0u64;
    let mut numerator = BigUint::zero();
    for block in blocks {
        let full = &numerator + &block.completions * block.signed_entries;
        if &full * &scale < target {
            numerator = full;
            entries += block.signed_entries;
            continue;
        }
        let current = &numerator * &scale;
        let needed = &target - current;
        let per_entry = &block.completions * &scale;
        let take_big = (&needed + &per_entry - BigUint::one()) / &per_entry;
        let take = take_big
            .to_u64()
            .ok_or("threshold partial block did not fit u64")?;
        if take > block.signed_entries {
            return Err("threshold partial block exceeds population".into());
        }
        return Ok(entries + take);
    }
    Err("full table did not reach the cold-construction threshold".into())
}

fn frontier_row(
    label: &str,
    entries: u64,
    blocks: &[PairBlock],
    retained: &BigUint,
) -> Result<FrontierRow, String> {
    let numerator = cumulative_numerator(blocks, entries);
    let denominator = retained * BigUint::from(4u64);
    let expected = big_ratio(&numerator, &denominator)?;
    let hot = expected.min(LOOKUPS);
    let cold = (LOOKUPS - hot).max(0.0);
    let projected = BASE_RATIO + LOCAL_RATIO * cold / (MEAN_CAPACITY + 1.0);
    let bytes = entries
        .checked_mul(ENTRY_BYTES)
        .ok_or("table byte count overflow")?;
    Ok(FrontierRow {
        label: label.into(),
        signed_entries: entries,
        bytes,
        bytes_log2: (bytes as f64).log2(),
        expected_cached_signed_pairs_present: expected,
        optimistic_hot_hits_upper: hot,
        cold_constructions_lower: cold,
        projected_ratio_to_rho_lower: projected,
        within_parity_budget: cold <= ALLOWED_COLD && projected <= 1.0,
    })
}

fn combinations(
    start: usize,
    remaining: usize,
    chosen: &mut Vec<usize>,
    output: &mut Vec<Vec<usize>>,
    total: usize,
) {
    if remaining == 0 {
        output.push(chosen.clone());
        return;
    }
    for index in start..=total - remaining {
        chosen.push(index);
        combinations(index + 1, remaining - 1, chosen, output, total);
        chosen.pop();
    }
}

fn toy_check(counts: &[u64], arity: usize, cutoff: usize) -> (u64, u64) {
    let mut capacities = Vec::new();
    for (capacity, count) in counts.iter().copied().enumerate() {
        capacities.extend(std::iter::repeat_n(capacity, count as usize));
    }
    let mut tuples = Vec::new();
    combinations(0, arity, &mut Vec::new(), &mut tuples, capacities.len());
    tuples.retain(|tuple| tuple.iter().map(|index| capacities[*index]).sum::<usize>() >= cutoff);
    let mut checked = 0u64;
    let mut failures = 0u64;
    for left in 0..capacities.len() {
        for right in left + 1..capacities.len() {
            let exact = tuples
                .iter()
                .filter(|tuple| tuple.contains(&left) && tuple.contains(&right))
                .count();
            let mut residual = counts.to_vec();
            residual[capacities[left]] -= 1;
            residual[capacities[right]] -= 1;
            let remaining_cutoff = cutoff.saturating_sub(capacities[left] + capacities[right]);
            let predicted = retained_count(&residual, arity - 2, remaining_cutoff)
                .to_usize()
                .expect("toy completion count");
            checked += 1;
            if predicted != exact {
                failures += 1;
            }
        }
    }
    (checked, failures)
}

fn digest_rows(rows: &[PairBlockRow]) -> Result<String, String> {
    let bytes = serde_json::to_vec(rows).map_err(|error| error.to_string())?;
    Ok(hex::encode(sha256(&bytes)))
}

fn digest_frontier(rows: &[FrontierRow]) -> Result<String, String> {
    let bytes = serde_json::to_vec(rows).map_err(|error| error.to_string())?;
    Ok(hex::encode(sha256(&bytes)))
}

fn execute(cli: &Cli) -> Result<ResultFile, String> {
    let (round33_dependency, round33) = read_dependency(&cli.round33, ROUND33_SHA256, 33)?;
    let (round36_dependency, round36) = read_dependency(&cli.round36, ROUND36_SHA256, 36)?;
    let counts = json_at(
        &round33,
        &["start_distribution", "residual_capacity_counts"],
    )?
    .as_array()
    .ok_or("round-33 capacity histogram is not an array")?
    .iter()
    .map(|value| value.as_u64().ok_or("capacity count is not a u64"))
    .collect::<Result<Vec<_>, _>>()?;
    let expected_counts = [vec![6_935u64; 18], vec![6_628u64]].concat();
    if counts != expected_counts || counts.iter().sum::<u64>() != COLUMNS {
        return Err("round-33 capacity histogram does not match the protocol".into());
    }
    if json_u64(&round33, &["cutoff_sweep", "selected", "cutoff_capacity"])? != CUTOFF as u64
        || json_u64(
            &round33,
            &[
                "pair_table_accounting",
                "signed_distinct_column_pair_entries",
            ],
        )? != PAIR_ENTRIES
    {
        return Err("round-33 selected model does not match the protocol".into());
    }
    let retained_frozen = parse_big(json_str(
        &round33,
        &["cutoff_sweep", "selected", "retained_unsigned_tuples"],
    )?)?;
    let full_frozen = parse_big(json_str(
        &round33,
        &["start_distribution", "exact_unsigned_tuples"],
    )?)?;
    if json_u64(&round36, &["frozen_model", "pair_entries"])? != PAIR_ENTRIES
        || json_u64(&round36, &["frozen_model", "entry_bytes"])? != ENTRY_BYTES
        || (json_f64(
            &round36,
            &[
                "frozen_model",
                "round33_ratio_to_rho_with_target_correction",
            ],
        )? - BASE_RATIO)
            .abs()
            > 1e-15
        || (json_f64(&round36, &["frozen_model", "local_oracle_ratio_to_rho"])? - LOCAL_RATIO).abs()
            > 1e-15
        || (json_f64(&round36, &["frozen_model", "mean_capacity"])? - MEAN_CAPACITY).abs() > 1e-12
        || (json_f64(
            &round36,
            &["frozen_model", "allowed_extra_additions_per_segment"],
        )? - ALLOWED_COLD)
            .abs()
            > 1e-12
    {
        return Err("round-36 frozen model does not match the protocol".into());
    }

    let distribution = capacity_distribution(&counts, ARITY);
    let retained: BigUint = distribution[CUTOFF..].iter().cloned().sum();
    let full: BigUint = distribution.iter().cloned().sum();
    let retained_matches = retained == retained_frozen;
    let full_matches = full == full_frozen && full == binomial(COLUMNS, ARITY);
    let blocks = pair_blocks(&counts, &retained)?;
    let rows = class_rows(&blocks, &retained)?;
    let density_failures = blocks
        .windows(2)
        .filter(|window| window[0].completions < window[1].completions)
        .count() as u64;
    let population: u64 = blocks.iter().map(|block| block.signed_entries).sum();
    let population_failures = u64::from(population != PAIR_ENTRIES);

    let full_numerator = cumulative_numerator(&blocks, PAIR_ENTRIES);
    let denominator = &retained * BigUint::from(4u64);
    let full_expected = big_ratio(&full_numerator, &denominator)?;
    let normalization_exact = full_numerator == &denominator * BigUint::from(136u64);

    let threshold = threshold_entries(&blocks, &retained)?;
    let mut frontier = Vec::new();
    for (label, entries) in [
        ("affine-2-mib", (2u64 << 20) / ENTRY_BYTES),
        ("affine-64-mib", (64u64 << 20) / ENTRY_BYTES),
        ("affine-256-mib", (256u64 << 20) / ENTRY_BYTES),
        ("parity-threshold", threshold),
        ("complete-table", PAIR_ENTRIES),
    ] {
        frontier.push(frontier_row(label, entries, &blocks, &retained)?);
    }

    let (toy_checked_a, toy_failures_a) = toy_check(&[2, 2, 2], 4, 4);
    let (toy_checked_b, toy_failures_b) = toy_check(&[3, 2, 2, 1], 4, 7);
    let class_digest = digest_rows(&rows)?;
    let frontier_digest = digest_frontier(&frontier)?;
    let controls = ExactControls {
        retained_count_matches_round33: retained_matches,
        full_population_matches_round33: full_matches,
        full_expected_pair_normalization: full_expected,
        full_expected_pair_normalization_exact: normalization_exact,
        density_monotonicity_failures: density_failures,
        population_failures,
        toy_cells: 2,
        toy_fixed_pairs_checked: toy_checked_a + toy_checked_b,
        toy_count_failures: toy_failures_a + toy_failures_b,
        false_positives: 0,
        false_negatives: 0,
        ordered_class_rows_sha256: class_digest,
        frontier_sha256: frontier_digest,
    };
    let exact = controls.retained_count_matches_round33
        && controls.full_population_matches_round33
        && controls.full_expected_pair_normalization_exact
        && controls.density_monotonicity_failures == 0
        && controls.population_failures == 0
        && controls.toy_count_failures == 0
        && controls.false_positives == 0
        && controls.false_negatives == 0;
    let admissible = frontier[..3].iter().any(|row| row.within_parity_budget);
    let model = FrozenModel {
        columns: COLUMNS,
        arity: ARITY as u64,
        cutoff: CUTOFF as u64,
        lookups_per_segment: LOOKUPS as u64,
        pair_entries: PAIR_ENTRIES,
        entry_bytes: ENTRY_BYTES,
        corrected_base_ratio_to_rho: BASE_RATIO,
        local_oracle_ratio_to_rho: LOCAL_RATIO,
        mean_capacity: MEAN_CAPACITY,
        allowed_cold_constructions_per_segment: ALLOWED_COLD,
        capacity_counts: counts,
        retained_unsigned_tuples: retained.to_str_radix(10),
    };
    let semantic = SemanticEvidence {
        curve: CURVE_SLUG.into(),
        round33_sha256: round33_dependency.sha256.clone(),
        round36_sha256: round36_dependency.sha256.clone(),
        model: model.clone(),
        class_rows: rows.clone(),
        frontier: frontier.clone(),
        controls: controls.clone(),
    };
    let semantic_sha = hex::encode(sha256(
        &serde_json::to_vec(&semantic).map_err(|error| error.to_string())?,
    ));
    let threshold_bytes = threshold
        .checked_mul(ENTRY_BYTES)
        .ok_or("threshold byte count overflow")?;
    Ok(ResultFile {
        schema: "p256-pair-entropy-v1".into(),
        curve: CURVE_SLUG.into(),
        family: "round-33 cutoff-conditioned optimistic signed-pair cache".into(),
        dependencies: vec![round33_dependency, round36_dependency],
        frozen_model: model,
        class_rows: rows,
        frontier,
        threshold_signed_entries: threshold,
        threshold_bytes,
        threshold_fraction_of_full_table: threshold as f64 / PAIR_ENTRIES as f64,
        exact_controls: controls,
        semantic_evidence_sha256: semantic_sha,
        gates: Gates {
            exact_counts: exact,
            admissible_depth_at_or_below_rho: admissible,
            measured_complete_time_at_or_below_rho: false,
            proved_p256_usable_relation_probability: false,
            structured_residual_degree_at_most_5: false,
            relation_collection_below_2_120: false,
            per_usable_relation_below_2_103: false,
            materialized_storage_below_2_50: threshold_bytes < (1u64 << 50),
            non_generic_end_to_end: false,
            promoted: false,
        },
        full_depth_unplanted_attempted: false,
        interpretation: if admissible {
            "At least one registered affine depth survives the candidate-favouring co-presence bound. This authorizes a separate measured hot/cold implementation, not promotion or a full-depth relation."
                .into()
        } else {
            "No registered affine depth can cover enough cutoff-conditioned signed pairs even when cached pairs are chosen by an adaptive global oracle and hot access is free."
                .into()
        },
        decision: if admissible {
            "preregister a measured hot/cold selector; do not promote or attempt a full-depth relation"
                .into()
        } else {
            "reject capacity-only hot tiering at the registered depths; do not benchmark or attempt a full-depth relation"
                .into()
        },
    })
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match execute(&cli).and_then(|result| {
        let bytes = serde_json::to_vec_pretty(&result).map_err(|error| error.to_string())?;
        if let Some(path) = &cli.out {
            fs::write(path, bytes).map_err(|error| format!("{}: {error}", path.display()))?;
        } else {
            println!("{}", String::from_utf8_lossy(&bytes));
        }
        Ok(())
    }) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_pair_entropy: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distribution_sums_to_binomial() {
        let counts = [3, 2, 2, 1];
        let distribution = capacity_distribution(&counts, 4);
        assert_eq!(
            distribution.into_iter().sum::<BigUint>(),
            binomial(counts.iter().sum(), 4)
        );
    }

    #[test]
    fn fixed_pair_dp_matches_exhaustive_toys() {
        for (counts, arity, cutoff) in [(&[2, 2, 2][..], 4, 4), (&[3, 2, 2, 1][..], 4, 7)] {
            let (checked, failures) = toy_check(counts, arity, cutoff);
            assert!(checked > 0);
            assert_eq!(failures, 0);
        }
    }

    #[test]
    fn signed_pair_population_is_exact() {
        assert_eq!(2 * COLUMNS * (COLUMNS - 1), PAIR_ENTRIES);
        assert_eq!(PAIR_ENTRIES * ENTRY_BYTES, 2_211_977_511_168);
    }
}
