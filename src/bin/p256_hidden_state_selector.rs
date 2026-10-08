//! Non-memoryless hidden-state audit for the P-256 scalar-orbit selector (round 30).

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::Serialize;
use serde_json::Value;

const ROUND25_SHA256: &str = "dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf";
const ROUND26_SHA256: &str = "463c100dd951e2a7a7a419118f84bde90fc419358d8279a430a5a0e6191db029";
const ROUND29_SHA256: &str = "9be6e5fc38644c34ad746dec9a6542da7e35ef2d8e11eb39a83de7af0f54f443";
const ROUND26_DELTA_SHA256: &str =
    "9bef47854a134951a3970dc45b31af2153ffea5b37dca9562577310fde01f9a4";
const GROUP_ORDER: &str =
    "115792089210356248762697446949407573529996955224135760342422259061068512044369";
const SCALAR: &str = "379483656589291393959605088584129393463128953089898693989044430844151895048";
const COLUMNS: u64 = 139_592;
const SIGNED_ATOMS: u64 = 279_184;
const LOCAL_ORACLE_RATIO: f64 = 0.964_336_477_130_181;
const RELATION_MODEL: &str =
    "17 distinct variable columns; signed 8+9 cross-colour collision modulo global negation";
const HASH_DOMAIN: &str = concat!(
    "icv1-fp256-t89188191154553853111372247798585809583-f188c491/",
    "hidden-state-round30/class/v1"
);

#[derive(Parser)]
#[command(about = "Audit compressed non-memoryless state for the P-256 S17 selector")]
struct Cli {
    #[arg(long)]
    round25: PathBuf,
    #[arg(long)]
    round26: PathBuf,
    #[arg(long)]
    round29: PathBuf,
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

#[derive(Serialize)]
struct Reconstruction {
    group_order: String,
    scalar: String,
    columns: u64,
    signed_atoms: u64,
    modular_multiplications: u64,
    ordered_delta_sha256: String,
    round26_digest_reproduced: bool,
    distinct_signed_deltas: u64,
    duplicate_signed_deltas: u64,
    scalar_cycle_closed: bool,
    half_cycle_is_negation: bool,
}

#[derive(Serialize)]
struct ExactStateQuotient {
    compatibility: String,
    exact_classes: u64,
    exact_class_information_bits: f64,
    minimum_fixed_key_bits: u32,
    atoms_per_class_minimum: u64,
    atoms_per_class_maximum: u64,
    compatible_pairs_replayed: u64,
    equal_pair_relations: u64,
    opposite_pair_relations: u64,
    replay_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    ordered_class_sha256: String,
    exact_rank_packed_bytes: u64,
    exact_coefficients_bytes: u64,
}

#[derive(Clone, Serialize)]
struct CompressionRow {
    prefix_bits: u32,
    bucket_capacity: u64,
    occupied_buckets: u64,
    empty_buckets: u64,
    maximum_bucket_width: u64,
    mean_occupied_bucket_width: f64,
    sum_squared_bucket_widths: u64,
    false_compatible_class_pairs_before_replay: u64,
    exact_replay_survival_probability: f64,
    exact_replay_rejection_probability: f64,
    balanced_oracle_sum_squared_widths: u64,
    balanced_oracle_survival_probability: f64,
    measured_ratio_to_rho_lower_bound: f64,
    balanced_ratio_to_rho_lower_bound: f64,
    materialized_prefix_key_bytes: u64,
    streamed_hash_input_bytes: u64,
    streamed_prefix_output_bytes: u64,
    complete_for_fixed_delta_set: bool,
}

#[derive(Serialize)]
struct CompressionLadder {
    hash_domain: String,
    rows: Vec<CompressionRow>,
    first_collision_free_prefix_bits: Option<u32>,
    best_measured_prefix_bits: u32,
    best_measured_ratio_to_rho: f64,
    best_balanced_prefix_bits: u32,
    best_balanced_ratio_to_rho: f64,
    exact_state_ratio_to_rho: f64,
    exact_state_overhead_log2_bits: f64,
    no_probabilistic_branch_discarded: bool,
}

#[derive(Serialize)]
struct ImportedBoundary {
    local_update_group_operations: f64,
    local_oracle_ratio_to_rho: f64,
    local_materialized_log2_bytes: f64,
    local_coalesces_without_hidden_state: bool,
    exact_log_quotient_dimension: u64,
    structured_residual_degree_of_regularity: Option<u64>,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    non_generic_end_to_end: bool,
}

#[derive(Serialize)]
struct PromotionGates {
    zero_false_positives_after_replay: bool,
    zero_false_negatives: bool,
    complete_state_key: bool,
    preserves_fixed_signed_s17: bool,
    exact_one_class_log_transport: bool,
    usable_collision_time_at_or_below_rho: bool,
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
    factor_base: String,
    dependencies: Vec<Dependency>,
    reconstruction: Reconstruction,
    exact_state_quotient: ExactStateQuotient,
    compression_ladder: CompressionLadder,
    imported_boundary: ImportedBoundary,
    promotion_gates: PromotionGates,
    full_depth_unplanted_attempted: bool,
    dominant_obstruction: String,
    decision: String,
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    BigUint::parse_bytes(value.as_bytes(), 10)
        .ok_or_else(|| format!("invalid frozen decimal integer: {value}"))
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

fn load_dependency(
    path: &Path,
    round: u64,
    expected_sha256: &str,
) -> Result<LoadedDependency, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != expected_sha256 {
        return Err(format!(
            "round-{round} dependency hash mismatch: expected {expected_sha256}, got {digest}"
        ));
    }
    let value: Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("{}: {error}", path.display()))?;
    let curve = value
        .get("curve")
        .or_else(|| value.get("root_curve"))
        .and_then(Value::as_str)
        .ok_or_else(|| format!("round-{round} dependency has no curve string"))?;
    if curve != CURVE_SLUG {
        return Err(format!("round-{round} curve mismatch: {curve}"));
    }
    if string_at(&value, &["relation_model"])? != RELATION_MODEL {
        return Err(format!("round-{round} relation-model mismatch"));
    }
    Ok(LoadedDependency {
        record: Dependency {
            round,
            path: path.display().to_string(),
            sha256: digest,
            schema: string_at(&value, &["schema"])?.into(),
        },
        value,
    })
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

fn reconstruct_deltas() -> Result<(Vec<BigUint>, Reconstruction), String> {
    let modulus = parse_big(GROUP_ORDER)?;
    let scalar = parse_big(SCALAR)?;
    let scalar_minus_one = &scalar - 1u8;
    let mut power = BigUint::one();
    let mut deltas = Vec::with_capacity(SIGNED_ATOMS as usize);
    let mut blob = Vec::with_capacity(SIGNED_ATOMS as usize * 32);
    for _ in 0..SIGNED_ATOMS {
        let delta = &scalar_minus_one * &power % &modulus;
        blob.extend(fixed_32(&delta)?);
        deltas.push(delta);
        power = power * &scalar % &modulus;
    }
    let cycle_closed = power.is_one();
    let half_cycle = scalar.modpow(&BigUint::from(COLUMNS), &modulus);
    let half_cycle_is_negation = half_cycle == &modulus - 1u8;
    let distinct = deltas.iter().collect::<BTreeSet<_>>().len() as u64;
    let digest = hex::encode(sha256(&blob));
    if digest != ROUND26_DELTA_SHA256 {
        return Err(format!(
            "round-26 delta digest mismatch: expected {ROUND26_DELTA_SHA256}, got {digest}"
        ));
    }
    if !cycle_closed || !half_cycle_is_negation || distinct != SIGNED_ATOMS {
        return Err("scalar-orbit delta reconstruction did not close exactly".into());
    }
    Ok((
        deltas,
        Reconstruction {
            group_order: GROUP_ORDER.into(),
            scalar: SCALAR.into(),
            columns: COLUMNS,
            signed_atoms: SIGNED_ATOMS,
            modular_multiplications: SIGNED_ATOMS * 2,
            ordered_delta_sha256: digest,
            round26_digest_reproduced: true,
            distinct_signed_deltas: distinct,
            duplicate_signed_deltas: SIGNED_ATOMS - distinct,
            scalar_cycle_closed: cycle_closed,
            half_cycle_is_negation,
        },
    ))
}

fn exact_quotient(deltas: &[BigUint]) -> Result<(Vec<BigUint>, ExactStateQuotient), String> {
    let modulus = parse_big(GROUP_ORDER)?;
    let mut classes = BTreeMap::<BigUint, Vec<BigUint>>::new();
    for delta in deltas {
        let negative = &modulus - delta;
        classes
            .entry(delta.clone().min(negative))
            .or_default()
            .push(delta.clone());
    }
    let mut class_blob = Vec::with_capacity(classes.len() * 32);
    let mut min_width = u64::MAX;
    let mut max_width = 0u64;
    let mut compatible_pairs = 0u64;
    let mut equal_pairs = 0u64;
    let mut opposite_pairs = 0u64;
    let mut failures = 0u64;
    for (class, members) in &classes {
        class_blob.extend(fixed_32(class)?);
        let width = members.len() as u64;
        min_width = min_width.min(width);
        max_width = max_width.max(width);
        for left in 0..members.len() {
            for right in left + 1..members.len() {
                compatible_pairs += 1;
                if members[left] == members[right] {
                    equal_pairs += 1;
                } else if (&members[left] + &members[right]) % &modulus == BigUint::zero() {
                    opposite_pairs += 1;
                } else {
                    failures += 1;
                }
            }
        }
    }
    if classes.len() as u64 != COLUMNS || min_width != 2 || max_width != 2 || failures != 0 {
        return Err("exact hidden-state quotient did not have the frozen two-atom classes".into());
    }
    let information_bits = (classes.len() as f64).log2();
    let minimum_bits = u64::BITS - (COLUMNS - 1).leading_zeros();
    let packed_bytes = (COLUMNS * u64::from(minimum_bits)).div_ceil(8);
    Ok((
        classes.keys().cloned().collect(),
        ExactStateQuotient {
            compatibility:
                "K_j=min(D_j,n-D_j); compatible transitions have equal or opposite deltas".into(),
            exact_classes: classes.len() as u64,
            exact_class_information_bits: information_bits,
            minimum_fixed_key_bits: minimum_bits,
            atoms_per_class_minimum: min_width,
            atoms_per_class_maximum: max_width,
            compatible_pairs_replayed: compatible_pairs,
            equal_pair_relations: equal_pairs,
            opposite_pair_relations: opposite_pairs,
            replay_failures: failures,
            false_positives: 0,
            false_negatives: 0,
            ordered_class_sha256: hex::encode(sha256(&class_blob)),
            exact_rank_packed_bytes: packed_bytes,
            exact_coefficients_bytes: COLUMNS * 32,
        },
    ))
}

fn prefix48(class: &BigUint) -> Result<u64, String> {
    let mut input = Vec::with_capacity(HASH_DOMAIN.len() + 32);
    input.extend(HASH_DOMAIN.as_bytes());
    input.extend(fixed_32(class)?);
    let digest = sha256(&input);
    let mut first = [0u8; 8];
    first.copy_from_slice(&digest[..8]);
    Ok(u64::from_be_bytes(first) >> 16)
}

fn bucket_for(prefix: u64, bits: u32) -> u64 {
    if bits == 0 {
        0
    } else {
        prefix >> (48 - bits)
    }
}

fn balanced_sum_squares(items: u64, buckets: u64) -> u64 {
    let used = buckets.min(items).max(1);
    let quotient = items / used;
    let remainder = items % used;
    remainder * (quotient + 1).pow(2) + (used - remainder) * quotient.pow(2)
}

fn compression_ladder(classes: &[BigUint]) -> Result<CompressionLadder, String> {
    let prefixes = classes
        .iter()
        .map(prefix48)
        .collect::<Result<Vec<_>, _>>()?;
    let mut rows = Vec::new();
    let mut first_collision_free = None;
    for bits in 0u32..=48 {
        let capacity = 1u64 << bits;
        let mut buckets = BTreeMap::<u64, u64>::new();
        for prefix in &prefixes {
            *buckets.entry(bucket_for(*prefix, bits)).or_default() += 1;
        }
        let occupied = buckets.len() as u64;
        let maximum = buckets.values().copied().max().unwrap_or(0);
        let sum_squares = buckets.values().map(|width| width * width).sum::<u64>();
        let false_pairs = buckets
            .values()
            .map(|width| width * (width - 1) / 2)
            .sum::<u64>();
        let survival = COLUMNS as f64 / sum_squares as f64;
        let balanced_squares = balanced_sum_squares(COLUMNS, capacity);
        let balanced_survival = COLUMNS as f64 / balanced_squares as f64;
        if false_pairs == 0 && first_collision_free.is_none() {
            first_collision_free = Some(bits);
        }
        let key_bytes = COLUMNS * u64::from(bits).div_ceil(8);
        rows.push(CompressionRow {
            prefix_bits: bits,
            bucket_capacity: capacity,
            occupied_buckets: occupied,
            empty_buckets: capacity - occupied,
            maximum_bucket_width: maximum,
            mean_occupied_bucket_width: COLUMNS as f64 / occupied as f64,
            sum_squared_bucket_widths: sum_squares,
            false_compatible_class_pairs_before_replay: false_pairs,
            exact_replay_survival_probability: survival,
            exact_replay_rejection_probability: 1.0 - survival,
            balanced_oracle_sum_squared_widths: balanced_squares,
            balanced_oracle_survival_probability: balanced_survival,
            measured_ratio_to_rho_lower_bound: LOCAL_ORACLE_RATIO * (sum_squares as f64).sqrt(),
            balanced_ratio_to_rho_lower_bound: LOCAL_ORACLE_RATIO
                * (balanced_squares as f64).sqrt(),
            materialized_prefix_key_bytes: key_bytes,
            streamed_hash_input_bytes: COLUMNS * 32,
            streamed_prefix_output_bytes: key_bytes,
            complete_for_fixed_delta_set: true,
        });
    }
    let best_measured = rows
        .iter()
        .min_by(|left, right| {
            left.measured_ratio_to_rho_lower_bound
                .total_cmp(&right.measured_ratio_to_rho_lower_bound)
                .then_with(|| left.prefix_bits.cmp(&right.prefix_bits))
        })
        .ok_or("empty measured compression ladder")?;
    let best_balanced = rows
        .iter()
        .min_by(|left, right| {
            left.balanced_ratio_to_rho_lower_bound
                .total_cmp(&right.balanced_ratio_to_rho_lower_bound)
                .then_with(|| left.prefix_bits.cmp(&right.prefix_bits))
        })
        .ok_or("empty balanced compression ladder")?;
    Ok(CompressionLadder {
        hash_domain: HASH_DOMAIN.into(),
        rows: rows.clone(),
        first_collision_free_prefix_bits: first_collision_free,
        best_measured_prefix_bits: best_measured.prefix_bits,
        best_measured_ratio_to_rho: best_measured.measured_ratio_to_rho_lower_bound,
        best_balanced_prefix_bits: best_balanced.prefix_bits,
        best_balanced_ratio_to_rho: best_balanced.balanced_ratio_to_rho_lower_bound,
        exact_state_ratio_to_rho: LOCAL_ORACLE_RATIO * (COLUMNS as f64).sqrt(),
        exact_state_overhead_log2_bits: 0.5 * (COLUMNS as f64).log2(),
        no_probabilistic_branch_discarded: true,
    })
}

fn terminal_named<'a>(round29: &'a Value, name: &str) -> Result<&'a Value, String> {
    at(round29, &["terminal_boundaries"])?
        .as_array()
        .ok_or("round-29 terminal_boundaries is not an array")?
        .iter()
        .find(|row| row.get("name").and_then(Value::as_str) == Some(name))
        .ok_or_else(|| format!("round-29 terminal boundary not found: {name}"))
}

fn imported_boundary(
    round25: &Value,
    round26: &Value,
    round29: &Value,
) -> Result<ImportedBoundary, String> {
    let local_ratio = f64_at(round26, &["parity_projection", "local_oracle_ratio_to_rho"])?;
    if (local_ratio - LOCAL_ORACLE_RATIO).abs() > 1e-15 {
        return Err("round-26 local oracle constant changed".into());
    }
    let known_log = terminal_named(round29, "known-log K=1 collision oracle")?;
    Ok(ImportedBoundary {
        local_update_group_operations: f64_at(
            round26,
            &["parity_projection", "local_update_group_operations"],
        )?,
        local_oracle_ratio_to_rho: local_ratio,
        local_materialized_log2_bytes: f64_at(known_log, &["projected_log2_bytes"])?,
        local_coalesces_without_hidden_state: bool_at(
            round26,
            &["parity_projection", "local_coalesces"],
        )?,
        exact_log_quotient_dimension: u64_at(
            round25,
            &["factor_base", "independent_log_quotient_dimension"],
        )?,
        structured_residual_degree_of_regularity: at(
            round29,
            &[
                "degree_boundary",
                "structured_residual_degree_of_regularity",
            ],
        )?
        .as_u64(),
        relation_collection_below_2_120: bool_at(
            round29,
            &["promotion_gates", "relation_collection_below_2_120"],
        )?,
        per_usable_relation_below_2_103: bool_at(
            round29,
            &["promotion_gates", "per_usable_relation_below_2_103"],
        )?,
        materialized_storage_below_2_50: bool_at(
            round29,
            &["promotion_gates", "materialized_storage_below_2_50"],
        )?,
        non_generic_end_to_end: bool_at(round29, &["promotion_gates", "non_generic"])?,
    })
}

fn run(cli: Cli) -> Result<ResultFile, String> {
    let round25 = load_dependency(&cli.round25, 25, ROUND25_SHA256)?;
    let round26 = load_dependency(&cli.round26, 26, ROUND26_SHA256)?;
    let round29 = load_dependency(&cli.round29, 29, ROUND29_SHA256)?;
    if string_at(&round25.value, &["factor_base", "fb_id"])? != "FB1he6b6b6e25de6"
        || u64_at(&round25.value, &["factor_base", "columns"])? != COLUMNS
        || u64_at(&round26.value, &["delta_census", "signed_atoms"])? != SIGNED_ATOMS
        || string_at(
            &round26.value,
            &["delta_census", "delta_coefficients_sha256"],
        )? != ROUND26_DELTA_SHA256
    {
        return Err("frozen scalar-orbit dependency facts changed".into());
    }

    let (deltas, reconstruction) = reconstruct_deltas()?;
    let (classes, exact_state) = exact_quotient(&deltas)?;
    let ladder = compression_ladder(&classes)?;
    let imported = imported_boundary(&round25.value, &round26.value, &round29.value)?;
    let usable_collision_time_at_or_below_rho = ladder.best_balanced_ratio_to_rho <= 1.0;
    let degree_gate = imported
        .structured_residual_degree_of_regularity
        .is_some_and(|degree| degree <= 5);
    let promoted = exact_state.false_positives == 0
        && exact_state.false_negatives == 0
        && usable_collision_time_at_or_below_rho
        && degree_gate
        && imported.relation_collection_below_2_120
        && imported.per_usable_relation_below_2_103
        && imported.materialized_storage_below_2_50
        && imported.non_generic_end_to_end;
    let result = ResultFile {
        schema: "p256.hidden_state_selector/v1".into(),
        curve: CURVE_SLUG.into(),
        relation_model: RELATION_MODEL.into(),
        factor_base: "FB1he6b6b6e25de6".into(),
        dependencies: vec![round25.record, round26.record, round29.record],
        reconstruction,
        exact_state_quotient: exact_state,
        compression_ladder: ladder,
        imported_boundary: imported,
        promotion_gates: PromotionGates {
            zero_false_positives_after_replay: true,
            zero_false_negatives: true,
            complete_state_key: true,
            preserves_fixed_signed_s17: true,
            exact_one_class_log_transport: true,
            usable_collision_time_at_or_below_rho,
            structured_residual_degree_at_most_5: degree_gate,
            relation_collection_below_2_120: false,
            per_usable_relation_below_2_103: false,
            materialized_storage_below_2_50: false,
            non_generic_end_to_end: false,
            promoted,
        },
        full_depth_unplanted_attempted: false,
        dominant_obstruction: "The 279,184 distinct oriented deltas leave 139,592 exact sign-folded compatibility classes. An exact lifted collision therefore pays at least sqrt(139592)=373.620 hidden-state overhead, or 360.296x rho after the optimistic local oracle constant. Compressing classes increases replay retries and cannot lower this bound.".into(),
        decision: "Reject compressed non-memoryless state for the scalar-orbit S17 selector. The exact state quotient is complete and replayable but hundreds of times over rho; shorter keys create incompatible-delta candidates whose exact-replay rejection cancels the apparent collision gain.".into(),
    };
    if promoted {
        return Err("hidden-state selector unexpectedly passed every promotion gate".into());
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
    fn balanced_partition_formula_is_exact() {
        assert_eq!(balanced_sum_squares(10, 1), 100);
        assert_eq!(balanced_sum_squares(10, 3), 34);
        assert_eq!(balanced_sum_squares(10, 10), 10);
        assert_eq!(balanced_sum_squares(10, 100), 10);
    }

    #[test]
    fn prefix_extraction_uses_high_bits() {
        let prefix = 0xaced_1234_5678u64;
        assert_eq!(bucket_for(prefix, 0), 0);
        assert_eq!(bucket_for(prefix, 4), 0xa);
        assert_eq!(bucket_for(prefix, 8), 0xac);
        assert_eq!(bucket_for(prefix, 48), prefix);
    }

    #[test]
    fn minimum_exact_key_is_eighteen_bits() {
        let bits = u64::BITS - (COLUMNS - 1).leading_zeros();
        assert_eq!(bits, 18);
        const {
            assert!(1u64 << 17 < COLUMNS);
            assert!(1u64 << 18 >= COLUMNS);
        }
    }

    #[test]
    fn exact_state_overhead_is_far_above_parity() {
        let ratio = LOCAL_ORACLE_RATIO * (COLUMNS as f64).sqrt();
        assert!(ratio > 360.0);
        assert!(ratio < 361.0);
    }

    #[test]
    fn compression_cannot_beat_exact_classes() {
        for bits in 0..=18 {
            let buckets = 1u64 << bits;
            assert!(balanced_sum_squares(COLUMNS, buckets) >= COLUMNS);
        }
    }
}
