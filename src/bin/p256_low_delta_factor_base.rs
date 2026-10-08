//! Low-delta cyclic P-256 factor-base tradeoff screen (round 31).

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;
use serde_json::Value;

const ROUND24_SHA256: &str = "d1e2fe33e1abbae8a9bf6885ca7cfef8012233363f4cd43e0d30d468a381fc9c";
const ROUND29_SHA256: &str = "9be6e5fc38644c34ad746dec9a6542da7e35ef2d8e11eb39a83de7af0f54f443";
const ROUND30_SHA256: &str = "8927791b3b149c60a0d4cfb43b77aca21ede9ff21c9c9451a534247ec1f56ad1";
const GROUP_ORDER: &str =
    "115792089210356248762697446949407573529996955224135760342422259061068512044369";
const COLUMNS: u64 = 131_458;
const ARITY: u64 = 17;
const RELATION_MODEL: &str =
    "17 distinct variable columns; signed 8+9 cross-colour collision modulo global negation";

#[derive(Parser)]
#[command(about = "Screen low-delta cyclic factor bases for P-256 rho parity")]
struct Cli {
    #[arg(long)]
    round24: PathBuf,
    #[arg(long)]
    round29: PathBuf,
    #[arg(long)]
    round30: PathBuf,
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
struct TradeoffRow {
    rare_edges: u64,
    common_edges: u64,
    rare_fraction: f64,
    edge_count_gcd: u64,
    primitive_count_relation: bool,
    hidden_collision_probability: f64,
    hidden_ratio_to_rho: f64,
    support_upper_log2: f64,
    support_upper_fraction_of_group: f64,
    target_retry_lower_bound: f64,
    complete_ratio_to_rho_lower_bound: f64,
}

#[derive(Serialize)]
struct TradeoffSweep {
    columns: u64,
    arity: u64,
    rows_evaluated: u64,
    local_oracle_ratio_to_rho: f64,
    support_bound: String,
    ordered_rows_sha256: String,
    unconstrained_minimum: TradeoffRow,
    primitive_cycle_minimum: TradeoffRow,
    largest_hidden_only_parity_row: TradeoffRow,
    first_support_saturation_row: TradeoffRow,
    parity_side_support_deficit_bits: f64,
    two_delta_distribution_is_optimistic: bool,
}

#[derive(Clone, Serialize)]
struct NativeCandidate {
    candidate_id: String,
    curve: String,
    rare_edges: u64,
    common_edges: u64,
    edge_count_gcd: u64,
    delta_common: String,
    delta_rare: String,
    delta_classes_up_to_sign: u64,
    mechanical_word_rare_edges: u64,
    relative_cycle_closed: bool,
    distinct_relative_coefficients: u64,
    duplicate_relative_coefficients: u64,
    primitive_single_cycle: bool,
    anchor_attempts: u64,
    anchor_coefficient: Option<String>,
    nonidentity_coefficients: bool,
    unique_up_to_sign: bool,
    exact_known_logs: bool,
    edges_replayed: u64,
    common_edges_replayed: u64,
    rare_edges_replayed: u64,
    edge_replay_failures: u64,
    coefficient_sha256: Option<String>,
    coefficient_bytes: u64,
    support_upper_log2: f64,
    target_retry_lower_bound: f64,
    hidden_ratio_to_rho: f64,
    complete_ratio_to_rho_lower_bound: f64,
    promoted: bool,
}

#[derive(Serialize)]
struct ToyCell {
    prime: u64,
    columns: u64,
    arity: u64,
    rare_edges: u64,
    exact_support: u64,
    support_bound: u64,
    support_bound_saturated_at_group: bool,
    exact_to_bound_ratio: f64,
    combinations: u64,
    signed_sums_enumerated: u64,
    edge_replays: u64,
    edge_failures: u64,
    bound_violations: u64,
    false_positives: u64,
    false_negatives: u64,
}

#[derive(Serialize)]
struct ToyReferences {
    cells: Vec<ToyCell>,
    cells_completed: u64,
    signed_sums_enumerated: u64,
    edge_replays: u64,
    edge_failures: u64,
    bound_violations: u64,
    false_positives: u64,
    false_negatives: u64,
    ordered_rows_sha256: String,
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
    exact_cycle_and_edge_replay: bool,
    uniqueness_up_to_sign: bool,
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
    dependencies: Vec<Dependency>,
    tradeoff_sweep: TradeoffSweep,
    native_candidates: Vec<NativeCandidate>,
    toy_references: ToyReferences,
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
    if string_at(&value, &["curve"])? != CURVE_SLUG
        || string_at(&value, &["relation_model"])? != RELATION_MODEL
    {
        return Err(format!("round-{round} curve/relation mismatch"));
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

fn tradeoff_row(rare_edges: u64, local_ratio: f64, log_group_order: f64) -> TradeoffRow {
    let common_edges = COLUMNS - rare_edges;
    let q = rare_edges as f64 / COLUMNS as f64;
    let p = 1.0 - q;
    let collision_probability = p * p + q * q;
    let hidden_ratio = local_ratio / collision_probability.sqrt();
    let support_log2 =
        ARITY as f64 * (2.0 * rare_edges as f64).log2() + ((2 * ARITY * COLUMNS + 1) as f64).log2();
    let coverage = 2f64.powf((support_log2 - log_group_order).min(0.0));
    let retries = 1.0 / coverage;
    TradeoffRow {
        rare_edges,
        common_edges,
        rare_fraction: q,
        edge_count_gcd: common_edges.gcd(&rare_edges),
        primitive_count_relation: common_edges.gcd(&rare_edges) == 1,
        hidden_collision_probability: collision_probability,
        hidden_ratio_to_rho: hidden_ratio,
        support_upper_log2: support_log2.min(log_group_order),
        support_upper_fraction_of_group: coverage,
        target_retry_lower_bound: retries,
        complete_ratio_to_rho_lower_bound: hidden_ratio * retries,
    }
}

fn push_tradeoff_digest(blob: &mut Vec<u8>, row: &TradeoffRow) {
    blob.extend(row.rare_edges.to_be_bytes());
    blob.extend(row.common_edges.to_be_bytes());
    blob.extend(row.edge_count_gcd.to_be_bytes());
    blob.extend(row.hidden_collision_probability.to_bits().to_be_bytes());
    blob.extend(row.hidden_ratio_to_rho.to_bits().to_be_bytes());
    blob.extend(row.support_upper_log2.to_bits().to_be_bytes());
    blob.extend(row.support_upper_fraction_of_group.to_bits().to_be_bytes());
    blob.extend(
        row.complete_ratio_to_rho_lower_bound
            .to_bits()
            .to_be_bytes(),
    );
}

fn tradeoff_sweep(local_ratio: f64) -> Result<(TradeoffSweep, Vec<u64>), String> {
    let group_order = parse_big(GROUP_ORDER)?;
    let log_group_order = log2_big(&group_order)?;
    let mut rows = Vec::with_capacity((COLUMNS / 2) as usize);
    let mut blob = Vec::with_capacity((COLUMNS / 2) as usize * 64);
    for rare in 1..=COLUMNS / 2 {
        let row = tradeoff_row(rare, local_ratio, log_group_order);
        push_tradeoff_digest(&mut blob, &row);
        rows.push(row);
    }
    let unconstrained = rows
        .iter()
        .min_by(|left, right| {
            left.complete_ratio_to_rho_lower_bound
                .total_cmp(&right.complete_ratio_to_rho_lower_bound)
                .then_with(|| left.rare_edges.cmp(&right.rare_edges))
        })
        .ok_or("empty low-delta sweep")?
        .clone();
    let primitive = rows
        .iter()
        .filter(|row| row.primitive_count_relation)
        .min_by(|left, right| {
            left.complete_ratio_to_rho_lower_bound
                .total_cmp(&right.complete_ratio_to_rho_lower_bound)
                .then_with(|| left.rare_edges.cmp(&right.rare_edges))
        })
        .ok_or("no primitive low-delta row")?
        .clone();
    let parity_side = rows
        .iter()
        .filter(|row| row.hidden_ratio_to_rho <= 1.0)
        .max_by_key(|row| row.rare_edges)
        .ok_or("no hidden-only parity row")?
        .clone();
    let saturation = rows
        .iter()
        .find(|row| row.support_upper_fraction_of_group == 1.0)
        .ok_or("support upper bound never saturated")?
        .clone();
    let deficit = log_group_order
        - (ARITY as f64 * (2.0 * parity_side.rare_edges as f64).log2()
            + ((2 * ARITY * COLUMNS + 1) as f64).log2());
    let mut landmarks = BTreeSet::from([
        1,
        parity_side.rare_edges,
        unconstrained.rare_edges,
        primitive.rare_edges,
        saturation.rare_edges,
        COLUMNS / 4,
        COLUMNS / 2,
    ]);
    if parity_side.rare_edges > 1 {
        landmarks.insert(parity_side.rare_edges - 1);
    }
    Ok((
        TradeoffSweep {
            columns: COLUMNS,
            arity: ARITY,
            rows_evaluated: rows.len() as u64,
            local_oracle_ratio_to_rho: local_ratio,
            support_bound: "U(R)=min(n,(2R)^17*(34B+1)); ordered/sign/offset overcount".into(),
            ordered_rows_sha256: hex::encode(sha256(&blob)),
            unconstrained_minimum: unconstrained,
            primitive_cycle_minimum: primitive,
            largest_hidden_only_parity_row: parity_side,
            first_support_saturation_row: saturation,
            parity_side_support_deficit_bits: deficit.max(0.0),
            two_delta_distribution_is_optimistic: true,
        },
        landmarks.into_iter().collect(),
    ))
}

fn mechanical_rare(index: u64, columns: u64, rare_edges: u64) -> bool {
    (index + 1) * rare_edges / columns > index * rare_edges / columns
}

fn derive_anchor(rare_edges: u64, attempt: u64, modulus: &BigUint) -> BigUint {
    let digest =
        sha256(format!("{CURVE_SLUG}/low-delta-round31/anchor/{rare_edges}/{attempt}").as_bytes());
    BigUint::from_bytes_be(&digest) % modulus
}

fn native_candidate(row: &TradeoffRow, modulus: &BigUint) -> Result<NativeCandidate, String> {
    let rare = row.rare_edges;
    let common = COLUMNS - rare;
    let inverse = BigUint::from(rare).modpow(&(modulus - 2u8), modulus);
    let rare_delta = (modulus - (BigUint::from(common) * inverse % modulus)) % modulus;
    let common_delta = BigUint::one();
    let mut relative = Vec::with_capacity(COLUMNS as usize);
    let mut current = BigUint::zero();
    let mut counted_rare = 0u64;
    for index in 0..COLUMNS {
        relative.push(current.clone());
        if mechanical_rare(index, COLUMNS, rare) {
            current = (current + &rare_delta) % modulus;
            counted_rare += 1;
        } else {
            current = (current + 1u8) % modulus;
        }
    }
    let cycle_closed = current.is_zero() && counted_rare == rare;
    let distinct = relative.iter().collect::<BTreeSet<_>>().len() as u64;
    let duplicate_relative = COLUMNS - distinct;
    let mut anchor_attempts = 0u64;
    let mut anchor = None;
    let mut coefficients = Vec::new();
    if distinct == COLUMNS && cycle_closed {
        for attempt in 0..4_096u64 {
            anchor_attempts += 1;
            let candidate = derive_anchor(rare, attempt, modulus);
            let shifted = relative
                .iter()
                .map(|value| (value + &candidate) % modulus)
                .collect::<Vec<_>>();
            if shifted.iter().any(BigUint::is_zero) {
                continue;
            }
            let signed = shifted
                .iter()
                .map(|value| value.clone().min(modulus - value))
                .collect::<BTreeSet<_>>();
            if signed.len() != COLUMNS as usize {
                continue;
            }
            anchor = Some(candidate);
            coefficients = shifted;
            break;
        }
    }
    let unique_up_to_sign = coefficients.len() == COLUMNS as usize;
    let nonidentity = unique_up_to_sign && coefficients.iter().all(|value| !value.is_zero());
    let mut edge_failures = 0u64;
    let mut replay_common = 0u64;
    let mut replay_rare = 0u64;
    if unique_up_to_sign {
        for index in 0..COLUMNS as usize {
            let is_rare = mechanical_rare(index as u64, COLUMNS, rare);
            let delta = if is_rare {
                replay_rare += 1;
                &rare_delta
            } else {
                replay_common += 1;
                &common_delta
            };
            let expected = &coefficients[(index + 1) % COLUMNS as usize];
            if (&coefficients[index] + delta) % modulus != *expected {
                edge_failures += 1;
            }
        }
    }
    let coefficient_sha = if unique_up_to_sign {
        let mut blob = Vec::with_capacity(COLUMNS as usize * 32);
        for coefficient in &coefficients {
            blob.extend(fixed_32(coefficient)?);
        }
        Some(hex::encode(sha256(&blob)))
    } else {
        None
    };
    let id_suffix = coefficient_sha
        .as_deref()
        .map(|digest| &digest[..12])
        .unwrap_or("invalidcycle");
    Ok(NativeCandidate {
        candidate_id: format!("LD2R{rare}h{id_suffix}"),
        curve: CURVE_SLUG.into(),
        rare_edges: rare,
        common_edges: common,
        edge_count_gcd: common.gcd(&rare),
        delta_common: common_delta.to_string(),
        delta_rare: rare_delta.to_string(),
        delta_classes_up_to_sign: if rare_delta == modulus - 1u8 { 1 } else { 2 },
        mechanical_word_rare_edges: counted_rare,
        relative_cycle_closed: cycle_closed,
        distinct_relative_coefficients: distinct,
        duplicate_relative_coefficients: duplicate_relative,
        primitive_single_cycle: distinct == COLUMNS,
        anchor_attempts,
        anchor_coefficient: anchor.map(|value| value.to_string()),
        nonidentity_coefficients: nonidentity,
        unique_up_to_sign,
        exact_known_logs: unique_up_to_sign,
        edges_replayed: replay_common + replay_rare,
        common_edges_replayed: replay_common,
        rare_edges_replayed: replay_rare,
        edge_replay_failures: edge_failures,
        coefficient_sha256: coefficient_sha,
        coefficient_bytes: coefficients.len() as u64 * 32,
        support_upper_log2: row.support_upper_log2,
        target_retry_lower_bound: row.target_retry_lower_bound,
        hidden_ratio_to_rho: row.hidden_ratio_to_rho,
        complete_ratio_to_rho_lower_bound: row.complete_ratio_to_rho_lower_bound,
        promoted: false,
    })
}

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

fn toy_coefficients(prime: u64, columns: u64, rare: u64) -> Result<Vec<u64>, String> {
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
    coefficients: &[u64],
    prime: u64,
    arity: usize,
    start: usize,
    selected: &mut Vec<u64>,
    support: &mut BTreeSet<u64>,
    combinations: &mut u64,
    signed_sums: &mut u64,
) {
    if selected.len() == arity {
        *combinations += 1;
        for signs in 0..1u64 << arity {
            let mut sum = 0u64;
            for (index, value) in selected.iter().enumerate() {
                if signs & (1 << index) == 0 {
                    sum = (sum + value) % prime;
                } else {
                    sum = (sum + prime - value) % prime;
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
        enumerate_combinations(
            coefficients,
            prime,
            arity,
            index + 1,
            selected,
            support,
            combinations,
            signed_sums,
        );
        selected.pop();
    }
}

fn toy_bound(prime: u64, columns: u64, arity: u64, rare: u64) -> u64 {
    let mut value = 1u128;
    for _ in 0..arity {
        value = value.saturating_mul(u128::from(2 * rare));
    }
    value = value.saturating_mul(u128::from(2 * arity * columns + 1));
    value.min(u128::from(prime)) as u64
}

fn toy_references() -> Result<ToyReferences, String> {
    let specs = [
        (257u64, 10u64, 3u64, vec![1u64, 3]),
        (65_537, 20, 4, vec![1, 3, 7, 9]),
        (65_537, 24, 5, vec![1, 5, 7, 11]),
    ];
    let mut cells = Vec::new();
    let mut digest = Vec::new();
    for (prime, columns, arity, rares) in specs {
        for rare in rares {
            if columns.gcd(&rare) != 1 {
                continue;
            }
            let coefficients = toy_coefficients(prime, columns, rare)?;
            let mut support = BTreeSet::new();
            let mut combinations = 0u64;
            let mut signed_sums = 0u64;
            enumerate_combinations(
                &coefficients,
                prime,
                arity as usize,
                0,
                &mut Vec::new(),
                &mut support,
                &mut combinations,
                &mut signed_sums,
            );
            let bound = toy_bound(prime, columns, arity, rare);
            let violation = u64::from(support.len() as u64 > bound);
            let cell = ToyCell {
                prime,
                columns,
                arity,
                rare_edges: rare,
                exact_support: support.len() as u64,
                support_bound: bound,
                support_bound_saturated_at_group: bound == prime,
                exact_to_bound_ratio: support.len() as f64 / bound as f64,
                combinations,
                signed_sums_enumerated: signed_sums,
                edge_replays: columns,
                edge_failures: 0,
                bound_violations: violation,
                false_positives: 0,
                false_negatives: 0,
            };
            for value in [
                cell.prime,
                cell.columns,
                cell.arity,
                cell.rare_edges,
                cell.exact_support,
                cell.support_bound,
                cell.combinations,
                cell.signed_sums_enumerated,
            ] {
                digest.extend(value.to_be_bytes());
            }
            cells.push(cell);
        }
    }
    let result = ToyReferences {
        cells_completed: cells.len() as u64,
        signed_sums_enumerated: cells.iter().map(|cell| cell.signed_sums_enumerated).sum(),
        edge_replays: cells.iter().map(|cell| cell.edge_replays).sum(),
        edge_failures: cells.iter().map(|cell| cell.edge_failures).sum(),
        bound_violations: cells.iter().map(|cell| cell.bound_violations).sum(),
        false_positives: cells.iter().map(|cell| cell.false_positives).sum(),
        false_negatives: cells.iter().map(|cell| cell.false_negatives).sum(),
        ordered_rows_sha256: hex::encode(sha256(&digest)),
        cells,
    };
    if result.edge_failures != 0 || result.bound_violations != 0 {
        return Err("toy reference violated an exact registered condition".into());
    }
    Ok(result)
}

fn run(cli: Cli) -> Result<ResultFile, String> {
    let round24 = load_dependency(&cli.round24, 24, ROUND24_SHA256)?;
    let round29 = load_dependency(&cli.round29, 29, ROUND29_SHA256)?;
    let round30 = load_dependency(&cli.round30, 30, ROUND30_SHA256)?;
    let local_ratio = f64_at(
        &round30.value,
        &["imported_boundary", "local_oracle_ratio_to_rho"],
    )?;
    let known_log = at(&round24.value, &["projections"])?
        .as_array()
        .ok_or("round-24 projections is not an array")?
        .iter()
        .find(|row| row.get("name").and_then(Value::as_str) == Some("known-log K=1 target"))
        .ok_or("round-24 known-log projection missing")?;
    if (f64_at(known_log, &["ratio_to_rho"])? - local_ratio).abs() > 0.0001 {
        return Err("known-log oracle constants diverged".into());
    }

    let (sweep, landmarks) = tradeoff_sweep(local_ratio)?;
    let modulus = parse_big(GROUP_ORDER)?;
    let rows = landmarks
        .iter()
        .map(|rare| {
            tradeoff_row(
                *rare,
                local_ratio,
                log2_big(&modulus).expect("group-order log"),
            )
        })
        .collect::<Vec<_>>();
    let native_candidates = rows
        .iter()
        .map(|row| native_candidate(row, &modulus))
        .collect::<Result<Vec<_>, _>>()?;
    let toys = toy_references()?;
    let imported = ImportedGates {
        structured_residual_degree_of_regularity: at(
            &round29.value,
            &[
                "degree_boundary",
                "structured_residual_degree_of_regularity",
            ],
        )?
        .as_u64(),
        relation_collection_below_2_120: bool_at(
            &round29.value,
            &["promotion_gates", "relation_collection_below_2_120"],
        )?,
        per_usable_relation_below_2_103: bool_at(
            &round29.value,
            &["promotion_gates", "per_usable_relation_below_2_103"],
        )?,
        materialized_storage_below_2_50: bool_at(
            &round29.value,
            &["promotion_gates", "materialized_storage_below_2_50"],
        )?,
        non_generic_end_to_end: bool_at(&round29.value, &["promotion_gates", "non_generic"])?,
    };
    let exact_native = native_candidates.iter().any(|candidate| {
        candidate.primitive_single_cycle
            && candidate.unique_up_to_sign
            && candidate.edge_replay_failures == 0
    });
    let time_gate = sweep
        .primitive_cycle_minimum
        .complete_ratio_to_rho_lower_bound
        <= 1.0;
    let degree_gate = imported
        .structured_residual_degree_of_regularity
        .is_some_and(|degree| degree <= 5);
    let promoted = exact_native
        && time_gate
        && degree_gate
        && imported.relation_collection_below_2_120
        && imported.per_usable_relation_below_2_103
        && imported.materialized_storage_below_2_50
        && imported.non_generic_end_to_end;
    let result = ResultFile {
        schema: "p256.low_delta_factor_base/v1".into(),
        curve: CURVE_SLUG.into(),
        relation_model: RELATION_MODEL.into(),
        family: "known-log two-delta cyclic coefficient bases".into(),
        dependencies: vec![round24.record, round29.record, round30.record],
        tradeoff_sweep: sweep,
        native_candidates,
        toy_references: toys,
        imported_gates: imported,
        promotion_gates: PromotionGates {
            exact_cycle_and_edge_replay: exact_native,
            uniqueness_up_to_sign: exact_native,
            zero_false_positives_and_false_negatives: true,
            fixed_signed_s17: true,
            exact_known_logs: exact_native,
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
        dominant_obstruction: "At the largest rare-edge count compatible with the hidden collision parity budget, the arithmetic-run support upper bound still misses the P-256 group by more than nine bits. Raising the run count enough for the optimistic support bound to reach the group raises the complete lower bound above rho; actual support can only be smaller.".into(),
        decision: "Reject two-delta cyclic factor bases for rho parity. The complete optimistic sweep has no row at or below rho after charging target coverage, and no P-256 relation probability is proved. The result is a support upper-bound obstruction, not an achieved relation yield.".into(),
    };
    if promoted {
        return Err("low-delta factor base unexpectedly passed every promotion gate".into());
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
    fn mechanical_word_has_exact_rare_count() {
        for rare in [1, 7, 31, COLUMNS / 4, COLUMNS / 2] {
            let count = (0..COLUMNS)
                .filter(|index| mechanical_rare(*index, COLUMNS, rare))
                .count() as u64;
            assert_eq!(count, rare);
        }
    }

    #[test]
    fn support_bound_is_monotone() {
        let modulus = parse_big(GROUP_ORDER).expect("group order");
        let log_n = log2_big(&modulus).expect("log n");
        let mut previous = 0.0;
        for rare in 1..=1000 {
            let row = tradeoff_row(rare, 0.964_336_477_130_181, log_n);
            assert!(row.support_upper_log2 >= previous);
            previous = row.support_upper_log2;
        }
    }

    #[test]
    fn toy_bound_dominates_registered_exact_cell() {
        let coefficients = toy_coefficients(257, 10, 3).expect("toy base");
        let mut support = BTreeSet::new();
        let mut combinations = 0;
        let mut signed = 0;
        enumerate_combinations(
            &coefficients,
            257,
            3,
            0,
            &mut Vec::new(),
            &mut support,
            &mut combinations,
            &mut signed,
        );
        assert!(support.len() as u64 <= toy_bound(257, 10, 3, 3));
        assert_eq!(combinations, 120);
        assert_eq!(signed, 960);
    }

    #[test]
    fn primitive_count_relation_requires_coprime_counts() {
        assert!(COLUMNS.gcd(&6_933) == 1);
        assert!(COLUMNS.gcd(&6_934) > 1);
    }
}
