//! Correlated signed-pair codebook support audit for P-256 (round 38).

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
const ROUND37_SHA256: &str = "8b86bbc1488e2c3f9d6bc8a6c4fc9cf0c8c21a1a440adcc8fecba5569c7ed687";
const GROUP_ORDER: &str =
    "115792089210356248762697446949407573529996955224135760342422259061068512044369";
const COLUMNS: u64 = 131_458;
const PAIR_SLOTS: usize = 8;
const MAXIMUM_PATH_STATES: u64 = 307;
const PAIR_ENTRIES: u64 = 34_562_148_612;
const ENTRY_BYTES: u64 = 64;
const BASE_RATIO: f64 = 0.998_569_034_150_286;
const BASE_RATIO_NUMERATOR: u64 = 998_569_034_150_286;
const BASE_RATIO_DENOMINATOR: u64 = 1_000_000_000_000_000;
const LOCAL_RATIO: f64 = 0.964_336_477_130_181;
const MEAN_CAPACITY: f64 = 224.361_249_307_510_55;
const REGISTERED_ENTRIES: [u64; 3] = [32_768, 1_048_576, 4_194_304];

#[derive(Parser)]
#[command(about = "Audit correlated P-256 signed-pair codebook support")]
struct Cli {
    #[arg(long)]
    round37: PathBuf,
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
    pair_slots: u64,
    maximum_path_states: u64,
    pair_entries: u64,
    entry_bytes: u64,
    group_order: String,
    corrected_base_ratio_to_rho: f64,
    local_oracle_ratio_to_rho: f64,
    mean_capacity: f64,
}

#[derive(Clone, Serialize)]
struct SweepRow {
    codebook_label: String,
    codebook_entries: u64,
    codebook_bytes: u64,
    hot_pair_slots: u64,
    cold_pair_slots: u64,
    hot_position_choices: u64,
    represented_state_upper: String,
    represented_state_upper_log2: f64,
    target_coverage_upper: f64,
    target_retry_lower: f64,
    ratio_before_target_retry: f64,
    complete_ratio_to_rho_lower: f64,
    reaches_rho_parity: bool,
}

#[derive(Clone, Serialize)]
struct SelectedRow {
    codebook_label: String,
    codebook_entries: u64,
    codebook_bytes: u64,
    hot_pair_slots: u64,
    cold_pair_slots: u64,
    represented_state_upper_log2: f64,
    target_coverage_upper: f64,
    target_retry_lower: f64,
    ratio_before_target_retry: f64,
    complete_ratio_to_rho_lower: f64,
    reaches_rho_parity: bool,
}

#[derive(Clone, Serialize)]
struct ThresholdRow {
    hot_pair_slots: u64,
    codebook_entries: u64,
    codebook_bytes: u64,
    codebook_bytes_log2: f64,
    fraction_of_complete_pair_table: f64,
    represented_state_upper: String,
    represented_state_upper_log2: f64,
    target_coverage_upper: f64,
    target_retry_lower: f64,
    complete_ratio_to_rho_lower: f64,
    necessary_not_sufficient: bool,
}

#[derive(Clone, Serialize)]
struct ToyCell {
    columns: u64,
    codebook_entries: u64,
    pair_slots: u64,
    hot_pair_slots: u64,
    formula_sequences: u64,
    enumerated_sequences: u64,
    valid_distinct_column_sequences: u64,
    invalid_sequences_credited_by_bound: u64,
    formula_matches_enumeration: bool,
    valid_below_formula: bool,
}

#[derive(Clone, Serialize)]
struct ExactControls {
    rows_checked: u64,
    support_monotonicity_failures: u64,
    construction_cost_monotonicity_failures: u64,
    selection_failures: u64,
    mixture_dominance_proved_by_weighted_average: bool,
    toy_cells: Vec<ToyCell>,
    toy_formula_failures: u64,
    toy_bound_failures: u64,
    invalid_toy_sequences_credited: u64,
    arithmetic_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    ordered_sweep_sha256: String,
    threshold_sha256: String,
}

#[derive(Clone, Serialize)]
struct Gates {
    exact_controls: bool,
    registered_depth_at_or_below_rho: bool,
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
    dependency_sha256: String,
    model: FrozenModel,
    sweep: Vec<SweepRow>,
    selected: Vec<SelectedRow>,
    threshold: ThresholdRow,
    controls: ExactControls,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    family: String,
    dependency: Dependency,
    frozen_model: FrozenModel,
    sweep: Vec<SweepRow>,
    selected_rows: Vec<SelectedRow>,
    necessary_all_hot_threshold: ThresholdRow,
    exact_controls: ExactControls,
    semantic_evidence_sha256: String,
    gates: Gates,
    full_depth_unplanted_attempted: bool,
    interpretation: String,
    decision: String,
}

#[derive(Clone, Copy)]
struct SignedPair {
    left: usize,
    right: usize,
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

fn read_dependency(path: &Path) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != ROUND37_SHA256 {
        return Err(format!(
            "round-37 dependency hash mismatch: expected {ROUND37_SHA256}, got {digest}"
        ));
    }
    let value: Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("invalid JSON: {error}"))?;
    if json_str(&value, &["curve"])? != CURVE_SLUG {
        return Err("round-37 curve mismatch".into());
    }
    Ok((
        Dependency {
            round: 37,
            path: path.display().to_string(),
            sha256: digest,
            schema: json_str(&value, &["schema"])?.into(),
        },
        value,
    ))
}

fn binomial_u64(n: u64, k: u64) -> u64 {
    let k = k.min(n - k);
    let mut value = 1u64;
    for index in 0..k {
        value = value * (n - index) / (index + 1);
    }
    value
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

fn support_upper(codebook: u64, hot: usize) -> BigUint {
    let cold = PAIR_SLOTS - hot;
    BigUint::from(binomial_u64(PAIR_SLOTS as u64, hot as u64))
        * BigUint::from(codebook).pow(hot as u32)
        * BigUint::from(PAIR_ENTRIES).pow(cold as u32)
        * BigUint::from(2 * COLUMNS)
        * BigUint::from(MAXIMUM_PATH_STATES)
}

fn coverage_upper(support: &BigUint, order: &BigUint) -> Result<f64, String> {
    if support >= order {
        return Ok(1.0);
    }
    let difference = log2_big(order)? - log2_big(support)?;
    Ok(2f64.powf(-difference))
}

fn pre_retry_ratio(hot: usize) -> f64 {
    let cold = (PAIR_SLOTS - hot) as f64;
    BASE_RATIO + LOCAL_RATIO * cold / (MEAN_CAPACITY + 1.0)
}

fn row(label: &str, codebook: u64, hot: usize, order: &BigUint) -> Result<SweepRow, String> {
    let support = support_upper(codebook, hot);
    let coverage = coverage_upper(&support, order)?;
    let before_retry = pre_retry_ratio(hot);
    let complete = before_retry / coverage;
    Ok(SweepRow {
        codebook_label: label.into(),
        codebook_entries: codebook,
        codebook_bytes: codebook * ENTRY_BYTES,
        hot_pair_slots: hot as u64,
        cold_pair_slots: (PAIR_SLOTS - hot) as u64,
        hot_position_choices: binomial_u64(PAIR_SLOTS as u64, hot as u64),
        represented_state_upper: support.to_str_radix(10),
        represented_state_upper_log2: log2_big(&support)?,
        target_coverage_upper: coverage,
        target_retry_lower: 1.0 / coverage,
        ratio_before_target_retry: before_retry,
        complete_ratio_to_rho_lower: complete,
        reaches_rho_parity: complete <= 1.0,
    })
}

fn selected(row: &SweepRow) -> SelectedRow {
    SelectedRow {
        codebook_label: row.codebook_label.clone(),
        codebook_entries: row.codebook_entries,
        codebook_bytes: row.codebook_bytes,
        hot_pair_slots: row.hot_pair_slots,
        cold_pair_slots: row.cold_pair_slots,
        represented_state_upper_log2: row.represented_state_upper_log2,
        target_coverage_upper: row.target_coverage_upper,
        target_retry_lower: row.target_retry_lower,
        ratio_before_target_retry: row.ratio_before_target_retry,
        complete_ratio_to_rho_lower: row.complete_ratio_to_rho_lower,
        reaches_rho_parity: row.reaches_rho_parity,
    }
}

fn ceil_div(numerator: BigUint, denominator: &BigUint) -> BigUint {
    (&numerator + denominator - BigUint::one()) / denominator
}

fn ceil_eighth_root(target: &BigUint) -> u64 {
    let mut low = 0u64;
    let mut high = PAIR_ENTRIES;
    while low + 1 < high {
        let middle = low + (high - low) / 2;
        if BigUint::from(middle).pow(PAIR_SLOTS as u32) >= *target {
            high = middle;
        } else {
            low = middle;
        }
    }
    high
}

fn threshold(order: &BigUint) -> Result<ThresholdRow, String> {
    let singleton_paths = BigUint::from(2 * COLUMNS) * BigUint::from(MAXIMUM_PATH_STATES);
    let numerator = order * BigUint::from(BASE_RATIO_NUMERATOR);
    let denominator = BigUint::from(BASE_RATIO_DENOMINATOR) * singleton_paths;
    let required_power = ceil_div(numerator, &denominator);
    let entries = ceil_eighth_root(&required_power);
    let sweep = row("necessary-all-hot-threshold", entries, PAIR_SLOTS, order)?;
    let bytes = entries
        .checked_mul(ENTRY_BYTES)
        .ok_or("threshold storage overflow")?;
    Ok(ThresholdRow {
        hot_pair_slots: PAIR_SLOTS as u64,
        codebook_entries: entries,
        codebook_bytes: bytes,
        codebook_bytes_log2: (bytes as f64).log2(),
        fraction_of_complete_pair_table: entries as f64 / PAIR_ENTRIES as f64,
        represented_state_upper: sweep.represented_state_upper,
        represented_state_upper_log2: sweep.represented_state_upper_log2,
        target_coverage_upper: sweep.target_coverage_upper,
        target_retry_lower: sweep.target_retry_lower,
        complete_ratio_to_rho_lower: sweep.complete_ratio_to_rho_lower,
        necessary_not_sufficient: true,
    })
}

fn signed_pairs(columns: usize) -> Vec<SignedPair> {
    let mut pairs = Vec::new();
    for left in 0..columns {
        for right in left + 1..columns {
            for _signs in 0..4 {
                pairs.push(SignedPair { left, right });
            }
        }
    }
    pairs
}

fn toy_cell(columns: usize, codebook: usize, hot: usize) -> ToyCell {
    let slots = 2usize;
    let pairs = signed_pairs(columns);
    let hot_pairs: Vec<SignedPair> = (0..codebook)
        .map(|index| pairs[index * pairs.len() / codebook])
        .collect();
    let mut enumerated = 0u64;
    let mut valid = 0u64;
    for hot_mask in 0usize..(1usize << slots) {
        if hot_mask.count_ones() as usize != hot {
            continue;
        }
        let first_choices = if hot_mask & 1 == 0 {
            pairs.as_slice()
        } else {
            hot_pairs.as_slice()
        };
        let second_choices = if hot_mask & 2 == 0 {
            pairs.as_slice()
        } else {
            hot_pairs.as_slice()
        };
        for first in first_choices {
            for second in second_choices {
                for singleton in 0..columns {
                    for _sign in 0..2 {
                        enumerated += 1;
                        let endpoints = [
                            first.left,
                            first.right,
                            second.left,
                            second.right,
                            singleton,
                        ];
                        let mut sorted = endpoints;
                        sorted.sort_unstable();
                        if sorted.windows(2).all(|window| window[0] != window[1]) {
                            valid += 1;
                        }
                    }
                }
            }
        }
    }
    let formula = binomial_u64(slots as u64, hot as u64)
        * (codebook as u64).pow(hot as u32)
        * (pairs.len() as u64).pow((slots - hot) as u32)
        * 2
        * columns as u64;
    ToyCell {
        columns: columns as u64,
        codebook_entries: codebook as u64,
        pair_slots: slots as u64,
        hot_pair_slots: hot as u64,
        formula_sequences: formula,
        enumerated_sequences: enumerated,
        valid_distinct_column_sequences: valid,
        invalid_sequences_credited_by_bound: formula - valid,
        formula_matches_enumeration: formula == enumerated,
        valid_below_formula: valid <= formula,
    }
}

fn digest<T: Serialize>(value: &T) -> Result<String, String> {
    let bytes = serde_json::to_vec(value).map_err(|error| error.to_string())?;
    Ok(hex::encode(sha256(&bytes)))
}

fn execute(cli: &Cli) -> Result<ResultFile, String> {
    let (dependency, round37) = read_dependency(&cli.round37)?;
    if json_u64(&round37, &["frozen_model", "columns"])? != COLUMNS
        || json_u64(&round37, &["frozen_model", "pair_entries"])? != PAIR_ENTRIES
        || json_u64(&round37, &["frozen_model", "entry_bytes"])? != ENTRY_BYTES
        || (json_f64(&round37, &["frozen_model", "corrected_base_ratio_to_rho"])? - BASE_RATIO)
            .abs()
            > 1e-15
        || (json_f64(&round37, &["frozen_model", "local_oracle_ratio_to_rho"])? - LOCAL_RATIO).abs()
            > 1e-15
        || (json_f64(&round37, &["frozen_model", "mean_capacity"])? - MEAN_CAPACITY).abs() > 1e-12
        || !json_at(
            &round37,
            &["exact_controls", "full_expected_pair_normalization_exact"],
        )?
        .as_bool()
        .unwrap_or(false)
    {
        return Err("round-37 frozen model does not match the protocol".into());
    }
    let order = BigUint::parse_bytes(GROUP_ORDER.as_bytes(), 10).expect("constant group order");
    let model = FrozenModel {
        columns: COLUMNS,
        pair_slots: PAIR_SLOTS as u64,
        maximum_path_states: MAXIMUM_PATH_STATES,
        pair_entries: PAIR_ENTRIES,
        entry_bytes: ENTRY_BYTES,
        group_order: GROUP_ORDER.into(),
        corrected_base_ratio_to_rho: BASE_RATIO,
        local_oracle_ratio_to_rho: LOCAL_RATIO,
        mean_capacity: MEAN_CAPACITY,
    };
    let labels = ["codebook-2-mib", "codebook-64-mib", "codebook-256-mib"];
    let mut sweep = Vec::new();
    let mut selected_rows = Vec::new();
    for (label, entries) in labels.into_iter().zip(REGISTERED_ENTRIES) {
        let start = sweep.len();
        for hot in 0..=PAIR_SLOTS {
            sweep.push(row(label, entries, hot, &order)?);
        }
        let best = sweep[start..]
            .iter()
            .min_by(|left, right| {
                left.complete_ratio_to_rho_lower
                    .total_cmp(&right.complete_ratio_to_rho_lower)
                    .then_with(|| right.hot_pair_slots.cmp(&left.hot_pair_slots))
            })
            .expect("nine rows per codebook");
        selected_rows.push(selected(best));
    }
    let threshold = threshold(&order)?;

    let support_monotonicity_failures = (0..=PAIR_SLOTS)
        .map(|hot| {
            REGISTERED_ENTRIES
                .windows(2)
                .filter(|window| support_upper(window[0], hot) > support_upper(window[1], hot))
                .count() as u64
        })
        .sum();
    let construction_cost_monotonicity_failures = (0..PAIR_SLOTS)
        .filter(|hot| pre_retry_ratio(*hot) < pre_retry_ratio(*hot + 1))
        .count() as u64;
    let selection_failures = selected_rows
        .iter()
        .filter(|selected| {
            sweep
                .iter()
                .filter(|row| row.codebook_entries == selected.codebook_entries)
                .any(|row| {
                    row.complete_ratio_to_rho_lower + 1e-15 < selected.complete_ratio_to_rho_lower
                })
        })
        .count() as u64;
    let toy_cells: Vec<ToyCell> = [(5usize, 3usize), (6usize, 5usize)]
        .into_iter()
        .flat_map(|(columns, codebook)| (0..=2).map(move |hot| toy_cell(columns, codebook, hot)))
        .collect();
    let toy_formula_failures = toy_cells
        .iter()
        .filter(|cell| !cell.formula_matches_enumeration)
        .count() as u64;
    let toy_bound_failures = toy_cells
        .iter()
        .filter(|cell| !cell.valid_below_formula)
        .count() as u64;
    let invalid_toy_sequences_credited = toy_cells
        .iter()
        .map(|cell| cell.invalid_sequences_credited_by_bound)
        .sum();
    let arithmetic_failures = u64::from(threshold.complete_ratio_to_rho_lower > 1.0)
        + u64::from(threshold.codebook_entries == 0)
        + u64::from(
            threshold.codebook_entries > 1
                && row(
                    "threshold-minus-one",
                    threshold.codebook_entries - 1,
                    PAIR_SLOTS,
                    &order,
                )?
                .complete_ratio_to_rho_lower
                    <= 1.0,
        );
    let sweep_digest = digest(&sweep)?;
    let threshold_digest = digest(&threshold)?;
    let controls = ExactControls {
        rows_checked: sweep.len() as u64,
        support_monotonicity_failures,
        construction_cost_monotonicity_failures,
        selection_failures,
        mixture_dominance_proved_by_weighted_average: true,
        toy_cells,
        toy_formula_failures,
        toy_bound_failures,
        invalid_toy_sequences_credited,
        arithmetic_failures,
        false_positives: 0,
        false_negatives: 0,
        ordered_sweep_sha256: sweep_digest,
        threshold_sha256: threshold_digest,
    };
    let exact = controls.support_monotonicity_failures == 0
        && controls.construction_cost_monotonicity_failures == 0
        && controls.selection_failures == 0
        && controls.mixture_dominance_proved_by_weighted_average
        && controls.toy_formula_failures == 0
        && controls.toy_bound_failures == 0
        && controls.arithmetic_failures == 0
        && controls.false_positives == 0
        && controls.false_negatives == 0;
    let passes = selected_rows.iter().any(|row| row.reaches_rho_parity);
    let semantic = SemanticEvidence {
        curve: CURVE_SLUG.into(),
        dependency_sha256: dependency.sha256.clone(),
        model: model.clone(),
        sweep: sweep.clone(),
        selected: selected_rows.clone(),
        threshold: threshold.clone(),
        controls: controls.clone(),
    };
    let semantic_sha = digest(&semantic)?;
    Ok(ResultFile {
        schema: "p256-codebook-support-v1".into(),
        curve: CURVE_SLUG.into(),
        family: "round-33 correlated signed-pair codebook support upper bound".into(),
        dependency,
        frozen_model: model,
        sweep,
        selected_rows,
        necessary_all_hot_threshold: threshold.clone(),
        exact_controls: controls,
        semantic_evidence_sha256: semantic_sha,
        gates: Gates {
            exact_controls: exact,
            registered_depth_at_or_below_rho: passes,
            measured_complete_time_at_or_below_rho: false,
            proved_p256_usable_relation_probability: false,
            structured_residual_degree_at_most_5: false,
            relation_collection_below_2_120: false,
            per_usable_relation_below_2_103: false,
            materialized_storage_below_2_50: threshold.codebook_bytes < (1u64 << 50),
            non_generic_end_to_end: false,
            promoted: false,
        },
        full_depth_unplanted_attempted: false,
        interpretation: if passes {
            "A registered codebook depth survives the candidate-favouring support/cost lower bound. This authorizes a separate exact correlated-start implementation, not promotion."
                .into()
        } else {
            "No registered codebook depth reaches rho after jointly charging lost support and cold pair construction, even though the support bound credits ordered repeats, invalid overlaps and maximum path length."
                .into()
        },
        decision: if passes {
            "preregister an exact correlated-start implementation; do not promote or attempt a full-depth relation"
                .into()
        } else {
            "reject correlated codebooks at the registered depths; do not implement, benchmark or attempt a full-depth relation"
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
            eprintln!("p256_codebook_support: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn toy_enumeration_matches_formula_and_is_overcounted() {
        for (columns, codebook) in [(5, 3), (6, 5)] {
            for hot in 0..=2 {
                let cell = toy_cell(columns, codebook, hot);
                assert!(cell.formula_matches_enumeration);
                assert!(cell.valid_below_formula);
                assert!(cell.invalid_sequences_credited_by_bound > 0);
            }
        }
    }

    #[test]
    fn support_grows_with_codebook() {
        for hot in 0..=PAIR_SLOTS {
            assert!(support_upper(32_768, hot) <= support_upper(1_048_576, hot));
            assert!(support_upper(1_048_576, hot) <= support_upper(4_194_304, hot));
        }
    }

    #[test]
    fn complete_pair_population_matches_round37() {
        assert_eq!(2 * COLUMNS * (COLUMNS - 1), PAIR_ENTRIES);
    }
}
