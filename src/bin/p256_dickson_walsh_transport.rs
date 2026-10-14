//! Complete global Walsh-character transport screen for the registered
//! depth-18 P-256 Dickson factor base (research round 313).

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{build_indexed, CURVE_SLUG};
use crypto_lib::ecc::p256_field::P256FieldElement as Fe;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::{CurveParams, Point};
use crypto_lib::hash::sha256;
use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, Signed, Zero};
use serde::Serialize;
use serde_json::{json, Value};

const SPEC: &str = "dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129";
const FB_ID: &str = "FB1h2f8621cda105";
const FB_SHA256: &str = "2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42";
const POINTS_SHA256: &str = "70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1";
const ROUND2_SHA256: &str = "dcaa1f66f5f89757a3670a4a8e158ca6d5c360b7ca119d384d2c4a0b2794c2c8";
const ROUND19_SHA256: &str = "3096540621408e4a48cfa18963ad01da9686d3527ee26776c8b6cf6f45a71114";
const ROUND312_SHA256: &str = "0cdec91f0097b498d1642c80c358aadac279c8a0bae951ed433381e1bf53db8e";
const PREREGISTERED_BUILDER_SHA256: &str =
    "a8740201f5f5c7dc56489827eb6f6a98d946338a512eb943741979f43eaf58c5";
const DEPTH: u32 = 18;
const FIBRE_SIZE: usize = 1 << DEPTH;
const COLUMNS: usize = 131_458;
const REQUESTED_ROWS: u64 = 138_031;
const RANK_MODULUS: u64 = (1u64 << 61) - 1;
const RANK_ALLOWANCE_NUMERATOR: u64 = 105;
const RANK_ALLOWANCE_DENOMINATOR: u64 = 100;
const RHO_S: f64 = 1.3;
const MULTIPLIERS: [u32; 6] = [1, 2, 3, 4, 8, 16];
const MAX_RANK_ROWS: usize = 256;
const MAX_RANK_ENTRIES: usize = 16_000_000;

#[derive(Parser)]
#[command(about = "Screen complete P-256 Dickson Walsh-character transport")]
struct Cli {
    #[arg(long)]
    round2: PathBuf,
    #[arg(long)]
    round19: PathBuf,
    #[arg(long)]
    round312: PathBuf,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Clone, Debug, Default, Serialize)]
struct WorkCounts {
    walsh_butterflies: u64,
    walsh_group_additions: u64,
    multiple_group_additions: u64,
    multiple_group_doublings: u64,
    normalization_batches: u64,
    batch_inversions: u64,
    normalization_multiplications: u64,
    canonical_key_comparisons: u64,
    output_collision_candidates: u64,
    output_multiple_candidates: u64,
    exact_rows_constructed: u64,
    zero_rows_discarded: u64,
    duplicate_rows_discarded: u64,
    replay_rows: u64,
    replay_additions: u64,
    replay_doublings: u64,
    replay_failures: u64,
    rank_modular_multiplications: u64,
    rank_modular_inversions: u64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct CanonicalPoint {
    x: [u8; 32],
    orientation: i8,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct OutputRecord {
    x: [u8; 32],
    orientation: i8,
    character: u32,
}

#[derive(Clone, Copy, Debug)]
struct MultipleRecord {
    column: u32,
    multiplier: u32,
    orientation: i8,
}

#[derive(Clone, Debug, Serialize)]
struct RelationReceipt {
    source: String,
    character_a: u32,
    character_b: Option<u32>,
    column: Option<u32>,
    multiplier: Option<u32>,
    support: usize,
    maximum_absolute_coefficient: i32,
    row_sha256: String,
    replayed: bool,
}

#[derive(Clone, Debug, Serialize)]
struct DependencyReceipt {
    path: String,
    sha256: String,
    schema: String,
}

#[derive(Clone, Debug, Serialize)]
struct RankReceipt {
    exact_rows: usize,
    active_columns: usize,
    rank_mod_2_61_minus_1: Option<usize>,
    reverse_rank_mod_2_61_minus_1: Option<usize>,
    rank_mod_p256_order: Option<usize>,
    reverse_rank_mod_p256_order: Option<usize>,
    promotion_credited_rank: usize,
    materialization_limit_rows: usize,
    materialization_limit_entries: usize,
    unset_due_to_limit: bool,
}

#[derive(Clone, Debug, Serialize)]
struct CyclicControl {
    modulus: u64,
    fibre_size: usize,
    active_slots: usize,
    optimized_relations: usize,
    reference_relations: usize,
    false_positives: usize,
    false_negatives: usize,
    replay_failures: usize,
    transform_mismatches: usize,
    relation_sha256: String,
}

#[derive(Clone, Debug, Serialize)]
struct P256Control {
    fibre_size: usize,
    active_slots: usize,
    detected_relations: usize,
    identity_character_zero_found: bool,
    character_one_double_match_found: bool,
    replayed_relations: usize,
    replay_failures: usize,
    relation_sha256: String,
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("invalid integer `{value}`"))
}

fn sha256_hex(bytes: &[u8]) -> String {
    hex::encode(sha256(bytes))
}

fn load_dependency(path: &Path, expected: &str) -> Result<(DependencyReceipt, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = sha256_hex(&bytes);
    if digest != expected {
        return Err(format!(
            "dependency hash mismatch for {}: expected {expected}, got {digest}",
            path.display()
        ));
    }
    let value: Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("{}: {error}", path.display()))?;
    let schema = value
        .get("schema")
        .and_then(Value::as_str)
        .ok_or_else(|| format!("{} has no schema", path.display()))?;
    Ok((
        DependencyReceipt {
            path: path.display().to_string(),
            sha256: digest,
            schema: schema.into(),
        },
        value,
    ))
}

fn low_point(row_x: &str, row_y: &str, curve: &CurveParams) -> Result<Projective, String> {
    let point = Point::Affine {
        x: curve.fe(parse_big(row_x)?),
        y: curve.fe(parse_big(row_y)?),
    };
    if !curve.is_on_curve(&point) {
        return Err("factor-base row is off P-256".into());
    }
    Ok(Projective::from_textbook(&point))
}

fn batch_normalize(points: &[Projective], counts: &mut WorkCounts) -> Vec<Option<CanonicalPoint>> {
    counts.normalization_batches += 1;
    let active = points
        .iter()
        .enumerate()
        .filter_map(|(index, point)| (!bool::from(point.is_identity())).then_some(index))
        .collect::<Vec<_>>();
    let mut inverses = vec![Fe::ZERO; points.len()];
    if !active.is_empty() {
        let mut prefixes = Vec::with_capacity(active.len());
        let mut product = Fe::ONE;
        for &index in &active {
            prefixes.push(product);
            product = product.mul(&points[index].z);
            counts.normalization_multiplications += 1;
        }
        let mut inverse = product.inv();
        counts.batch_inversions += 1;
        for position in (0..active.len()).rev() {
            let index = active[position];
            inverses[index] = inverse.mul(&prefixes[position]);
            counts.normalization_multiplications += 1;
            if position != 0 {
                inverse = inverse.mul(&points[index].z);
                counts.normalization_multiplications += 1;
            }
        }
    }
    points
        .iter()
        .enumerate()
        .map(|(index, point)| {
            if bool::from(point.is_identity()) {
                return None;
            }
            let x = point.x.mul(&inverses[index]).to_bytes_be();
            let y = point.y.mul(&inverses[index]).to_bytes_be();
            counts.normalization_multiplications += 2;
            let neg_y = point.y.mul(&inverses[index]).neg().to_bytes_be();
            counts.normalization_multiplications += 1;
            Some(CanonicalPoint {
                x,
                orientation: if y <= neg_y { 1 } else { -1 },
            })
        })
        .collect()
}

fn affine_span_dimension(active: &[u32]) -> usize {
    let Some(&origin) = active.first() else {
        return 0;
    };
    let mut basis = [0u32; DEPTH as usize];
    let mut rank = 0usize;
    for &index in &active[1..] {
        let mut value = index ^ origin;
        while value != 0 {
            let bit = 31 - value.leading_zeros();
            let slot = bit as usize;
            if basis[slot] == 0 {
                basis[slot] = value;
                rank += 1;
                break;
            }
            value ^= basis[slot];
        }
    }
    rank
}

fn walsh_transform(points: &mut [Projective], counts: &mut WorkCounts) {
    assert!(points.len().is_power_of_two());
    let mut half = 1usize;
    while half < points.len() {
        let block = half * 2;
        for start in (0..points.len()).step_by(block) {
            for offset in 0..half {
                let left = points[start + offset];
                let right = points[start + half + offset];
                points[start + offset] = left.add(&right);
                points[start + half + offset] = left.add(&right.neg());
                counts.walsh_butterflies += 1;
                counts.walsh_group_additions += 2;
            }
        }
        half = block;
    }
}

fn walsh_sign(character: u32, index: u32) -> i32 {
    if (character & index).count_ones().is_multiple_of(2) {
        1
    } else {
        -1
    }
}

fn normalize_row(mut row: Vec<i32>) -> Option<Vec<i32>> {
    let mut gcd = 0i32;
    for value in &row {
        gcd = gcd_i32(gcd, value.abs());
    }
    if gcd == 0 {
        return None;
    }
    if gcd > 1 {
        for value in &mut row {
            *value /= gcd;
        }
    }
    if row
        .iter()
        .find(|value| **value != 0)
        .is_some_and(|v| *v < 0)
    {
        for value in &mut row {
            *value = -*value;
        }
    }
    Some(row)
}

fn gcd_i32(mut left: i32, mut right: i32) -> i32 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left.abs()
}

fn row_digest(row: &[i32]) -> String {
    let mut bytes = Vec::new();
    for (column, coefficient) in row.iter().enumerate().filter(|(_, value)| **value != 0) {
        bytes.extend_from_slice(&(column as u32).to_be_bytes());
        bytes.extend_from_slice(&coefficient.to_be_bytes());
    }
    sha256_hex(&bytes)
}

fn collision_row(
    character_a: u32,
    orientation_a: i8,
    character_b: u32,
    orientation_b: i8,
    enumeration: &[Option<u32>],
) -> Vec<i32> {
    let mut row = vec![0i32; COLUMNS];
    for (index, column) in enumeration.iter().enumerate() {
        let Some(column) = column else { continue };
        row[*column as usize] = i32::from(orientation_b) * walsh_sign(character_a, index as u32)
            - i32::from(orientation_a) * walsh_sign(character_b, index as u32);
    }
    row
}

fn identity_row(character: u32, enumeration: &[Option<u32>]) -> Vec<i32> {
    let mut row = vec![0i32; COLUMNS];
    for (index, column) in enumeration.iter().enumerate() {
        if let Some(column) = column {
            row[*column as usize] = walsh_sign(character, index as u32);
        }
    }
    row
}

fn multiple_row(
    character: u32,
    output_orientation: i8,
    multiple: MultipleRecord,
    enumeration: &[Option<u32>],
) -> Vec<i32> {
    let mut row = identity_row(character, enumeration);
    for value in &mut row {
        *value *= i32::from(multiple.orientation);
    }
    row[multiple.column as usize] -= i32::from(output_orientation) * multiple.multiplier as i32;
    row
}

fn small_scalar_mul(point: Projective, mut scalar: u32, counts: &mut WorkCounts) -> Projective {
    let mut result = Projective::IDENTITY;
    let mut base = point;
    while scalar != 0 {
        if scalar & 1 == 1 {
            result = result.add(&base);
            counts.replay_additions += 1;
        }
        scalar >>= 1;
        if scalar != 0 {
            base = base.double();
            counts.replay_doublings += 1;
        }
    }
    result
}

fn replay_row(row: &[i32], columns: &[Projective], counts: &mut WorkCounts) -> bool {
    let mut sum = Projective::IDENTITY;
    for (column, coefficient) in row.iter().enumerate() {
        if *coefficient == 0 {
            continue;
        }
        let mut term = small_scalar_mul(columns[column], coefficient.unsigned_abs(), counts);
        if *coefficient < 0 {
            term = term.neg();
        }
        sum = sum.add(&term);
        counts.replay_additions += 1;
    }
    counts.replay_rows += 1;
    let valid = bool::from(sum.is_identity());
    if !valid {
        counts.replay_failures += 1;
    }
    valid
}

fn insert_relation(
    raw: Vec<i32>,
    mut receipt: RelationReceipt,
    rows: &mut BTreeMap<Vec<i32>, RelationReceipt>,
    columns: &[Projective],
    counts: &mut WorkCounts,
) -> Result<(), String> {
    counts.exact_rows_constructed += 1;
    let Some(row) = normalize_row(raw) else {
        counts.zero_rows_discarded += 1;
        return Ok(());
    };
    if rows.contains_key(&row) {
        counts.duplicate_rows_discarded += 1;
        return Ok(());
    }
    let replayed = replay_row(&row, columns, counts);
    if !replayed {
        return Err(format!(
            "Walsh relation replay failed for {}",
            receipt.source
        ));
    }
    receipt.support = row.iter().filter(|value| **value != 0).count();
    receipt.maximum_absolute_coefficient = row.iter().map(|v| v.abs()).max().unwrap_or(0);
    receipt.row_sha256 = row_digest(&row);
    receipt.replayed = true;
    rows.insert(row, receipt);
    Ok(())
}

fn build_multiples(
    columns: &[Projective],
    counts: &mut WorkCounts,
) -> (Vec<Projective>, Vec<(u32, u32)>) {
    let mut points = Vec::with_capacity(columns.len() * MULTIPLIERS.len());
    let mut metadata = Vec::with_capacity(points.capacity());
    for (column, point) in columns.iter().copied().enumerate() {
        let p2 = point.double();
        let p3 = p2.add(&point);
        let p4 = p2.double();
        let p8 = p4.double();
        let p16 = p8.double();
        counts.multiple_group_additions += 1;
        counts.multiple_group_doublings += 4;
        for (multiplier, value) in [(1, point), (2, p2), (3, p3), (4, p4), (8, p8), (16, p16)] {
            points.push(value);
            metadata.push((column as u32, multiplier));
        }
    }
    (points, metadata)
}

fn coefficient_mod(value: i32, modulus: &BigUint) -> BigUint {
    if value >= 0 {
        BigUint::from(value as u32) % modulus
    } else {
        let magnitude = BigUint::from(value.unsigned_abs()) % modulus;
        if magnitude.is_zero() {
            magnitude
        } else {
            modulus - magnitude
        }
    }
}

fn sparse_rank(
    rows: &[Vec<i32>],
    modulus: &BigUint,
    reverse: bool,
    counts: &mut WorkCounts,
) -> usize {
    let mut pivots = BTreeMap::<u32, BTreeMap<u32, BigUint>>::new();
    let order: Box<dyn Iterator<Item = &Vec<i32>>> = if reverse {
        Box::new(rows.iter().rev())
    } else {
        Box::new(rows.iter())
    };
    for dense in order {
        let mut row = dense
            .iter()
            .enumerate()
            .filter(|(_, value)| **value != 0)
            .map(|(column, value)| (column as u32, coefficient_mod(*value, modulus)))
            .collect::<BTreeMap<_, _>>();
        while let Some((&pivot, coefficient)) = row.first_key_value() {
            let coefficient = coefficient.clone();
            if let Some(base) = pivots.get(&pivot).cloned() {
                for (column, base_coefficient) in base {
                    let product = (&coefficient * base_coefficient) % modulus;
                    counts.rank_modular_multiplications += 1;
                    let current = row.remove(&column).unwrap_or_else(BigUint::zero);
                    let reduced = if current >= product {
                        current - product
                    } else {
                        current + modulus - product
                    };
                    if !reduced.is_zero() {
                        row.insert(column, reduced);
                    }
                }
            } else {
                let inverse = coefficient.modpow(&(modulus - BigUint::from(2u8)), modulus);
                counts.rank_modular_inversions += 1;
                for value in row.values_mut() {
                    *value = (&*value * &inverse) % modulus;
                    counts.rank_modular_multiplications += 1;
                }
                pivots.insert(pivot, row);
                break;
            }
        }
    }
    pivots.len()
}

fn rank_receipt(rows: &[Vec<i32>], counts: &mut WorkCounts) -> RankReceipt {
    let active_columns = rows
        .iter()
        .flat_map(|row| {
            row.iter()
                .enumerate()
                .filter_map(|(column, value)| (*value != 0).then_some(column))
        })
        .collect::<BTreeSet<_>>()
        .len();
    let entries = rows.len().saturating_mul(active_columns);
    let limited = rows.len() > MAX_RANK_ROWS || entries > MAX_RANK_ENTRIES;
    let (rank61, reverse61, rankn, reversen) = if limited {
        (None, None, None, None)
    } else {
        let modulus61 = BigUint::from(RANK_MODULUS);
        let modulusn = CurveParams::p256().n;
        (
            Some(sparse_rank(rows, &modulus61, false, counts)),
            Some(sparse_rank(rows, &modulus61, true, counts)),
            Some(sparse_rank(rows, &modulusn, false, counts)),
            Some(sparse_rank(rows, &modulusn, true, counts)),
        )
    };
    let credit = match (rank61, reverse61, rankn, reversen) {
        (Some(a), Some(b), Some(c), Some(d)) if a == b && c == d => a.min(c),
        _ => 0,
    };
    RankReceipt {
        exact_rows: rows.len(),
        active_columns,
        rank_mod_2_61_minus_1: rank61,
        reverse_rank_mod_2_61_minus_1: reverse61,
        rank_mod_p256_order: rankn,
        reverse_rank_mod_p256_order: reversen,
        promotion_credited_rank: credit,
        materialization_limit_rows: MAX_RANK_ROWS,
        materialization_limit_entries: MAX_RANK_ENTRIES,
        unset_due_to_limit: limited,
    }
}

fn normalize_small_row(mut row: Vec<i64>) -> Option<Vec<i64>> {
    let mut gcd = 0i64;
    for value in &row {
        gcd = gcd_i64(gcd, value.abs());
    }
    if gcd == 0 {
        return None;
    }
    for value in &mut row {
        *value /= gcd;
    }
    if row
        .iter()
        .find(|value| **value != 0)
        .is_some_and(|v| *v < 0)
    {
        for value in &mut row {
            *value = -*value;
        }
    }
    Some(row)
}

fn gcd_i64(mut left: i64, mut right: i64) -> i64 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left.abs()
}

fn scalar_walsh(labels: &[u64], modulus: u64) -> Vec<u64> {
    let mut values = labels.to_vec();
    let mut half = 1usize;
    while half < values.len() {
        let block = half * 2;
        for start in (0..values.len()).step_by(block) {
            for offset in 0..half {
                let left = values[start + offset];
                let right = values[start + half + offset];
                values[start + offset] = ((left as u128 + right as u128) % modulus as u128) as u64;
                values[start + half + offset] =
                    ((left as u128 + modulus as u128 - right as u128) % modulus as u128) as u64;
            }
        }
        half = block;
    }
    values
}

fn direct_scalar_walsh(labels: &[u64], modulus: u64) -> Vec<u64> {
    (0..labels.len() as u32)
        .map(|character| {
            labels.iter().enumerate().fold(0u64, |sum, (index, label)| {
                if walsh_sign(character, index as u32) > 0 {
                    ((sum as u128 + *label as u128) % modulus as u128) as u64
                } else {
                    ((sum as u128 + modulus as u128 - *label as u128) % modulus as u128) as u64
                }
            })
        })
        .collect()
}

fn scalar_relation_set(outputs: &[u64], labels: &[u64], modulus: u64) -> BTreeSet<Vec<i64>> {
    let n = outputs.len();
    let canonical = outputs
        .iter()
        .map(|value| {
            if *value == 0 {
                (0u64, 0i64)
            } else if *value <= modulus - *value {
                (*value, 1)
            } else {
                (modulus - *value, -1)
            }
        })
        .collect::<Vec<_>>();
    let mut rows = BTreeSet::new();
    for character in 0..n {
        if outputs[character] == 0 {
            let row = (0..n)
                .map(|index| i64::from(walsh_sign(character as u32, index as u32)))
                .collect();
            if let Some(row) = normalize_small_row(row) {
                rows.insert(row);
            }
        }
        for other in character + 1..n {
            if canonical[character].0 != canonical[other].0 {
                continue;
            }
            let row = (0..n)
                .map(|index| {
                    canonical[other].1 * i64::from(walsh_sign(character as u32, index as u32))
                        - canonical[character].1 * i64::from(walsh_sign(other as u32, index as u32))
                })
                .collect();
            if let Some(row) = normalize_small_row(row) {
                rows.insert(row);
            }
        }
        for (column, label) in labels.iter().enumerate() {
            for multiplier in [1i64, 2, 3, 4] {
                let value = ((*label as u128 * multiplier as u128) % modulus as u128) as u64;
                let (key, orientation) = if value == 0 {
                    (0, 0)
                } else if value <= modulus - value {
                    (value, 1)
                } else {
                    (modulus - value, -1)
                };
                if key == 0 || key != canonical[character].0 {
                    continue;
                }
                let mut row = (0..n)
                    .map(|index| {
                        orientation * i64::from(walsh_sign(character as u32, index as u32))
                    })
                    .collect::<Vec<_>>();
                row[column] -= canonical[character].1 * multiplier;
                if let Some(row) = normalize_small_row(row) {
                    rows.insert(row);
                }
            }
        }
    }
    rows
}

/// Deliberately direct reference: compare the uncanonicalized scalar outputs
/// and construct the coefficient equations from their defining equalities.
/// This does not reuse the optimized detector's key/orientation algebra.
fn scalar_reference_relation_set(labels: &[u64], modulus: u64) -> BTreeSet<Vec<i64>> {
    let outputs = direct_scalar_walsh(labels, modulus);
    let n = outputs.len();
    let mut rows = BTreeSet::new();
    for character in 0..n {
        if outputs[character] == 0 {
            let row = (0..n)
                .map(|index| i64::from(walsh_sign(character as u32, index as u32)))
                .collect();
            if let Some(row) = normalize_small_row(row) {
                rows.insert(row);
            }
        }
        for other in character + 1..n {
            // Identity outputs are already represented by their defining
            // character rows.  The optimized collision table contains only
            // nonidentity canonical keys, so compare the same disjoint class.
            if outputs[character] == 0 || outputs[other] == 0 {
                continue;
            }
            let sign = if outputs[character] == outputs[other] {
                Some(-1i64)
            } else if outputs[character] == (modulus - outputs[other]) % modulus {
                Some(1i64)
            } else {
                None
            };
            let Some(sign) = sign else { continue };
            let row = (0..n)
                .map(|index| {
                    i64::from(walsh_sign(character as u32, index as u32))
                        + sign * i64::from(walsh_sign(other as u32, index as u32))
                })
                .collect();
            if let Some(row) = normalize_small_row(row) {
                rows.insert(row);
            }
        }
        for (column, label) in labels.iter().enumerate() {
            for multiplier in [1i64, 2, 3, 4] {
                let multiple = ((*label as u128 * multiplier as u128) % modulus as u128) as u64;
                let coefficient = if outputs[character] == multiple {
                    Some(-multiplier)
                } else if outputs[character] == (modulus - multiple) % modulus {
                    Some(multiplier)
                } else {
                    None
                };
                let Some(coefficient) = coefficient else {
                    continue;
                };
                let mut row = (0..n)
                    .map(|index| i64::from(walsh_sign(character as u32, index as u32)))
                    .collect::<Vec<_>>();
                row[column] += coefficient;
                if let Some(row) = normalize_small_row(row) {
                    rows.insert(row);
                }
            }
        }
    }
    rows
}

fn cyclic_control() -> CyclicControl {
    let modulus = 101u64;
    let labels = (0..16u64)
        .map(|index| (7 * index * index + 11 * index + 3) % modulus)
        .collect::<Vec<_>>();
    let optimized_outputs = scalar_walsh(&labels, modulus);
    let direct_outputs = direct_scalar_walsh(&labels, modulus);
    let transform_mismatches = optimized_outputs
        .iter()
        .zip(&direct_outputs)
        .filter(|(left, right)| left != right)
        .count();
    let optimized = scalar_relation_set(&optimized_outputs, &labels, modulus);
    let reference = scalar_reference_relation_set(&labels, modulus);
    let false_positives = optimized.difference(&reference).count();
    let false_negatives = reference.difference(&optimized).count();
    let replay_failures = optimized
        .iter()
        .filter(|row| {
            let sum = row
                .iter()
                .zip(&labels)
                .fold(0i128, |acc, (coefficient, label)| {
                    acc + i128::from(*coefficient) * i128::from(*label)
                });
            sum.rem_euclid(i128::from(modulus)) != 0
        })
        .count();
    let mut bytes = Vec::new();
    for row in &optimized {
        for value in row {
            bytes.extend_from_slice(&value.to_be_bytes());
        }
    }
    CyclicControl {
        modulus,
        fibre_size: labels.len(),
        active_slots: labels.len(),
        optimized_relations: optimized.len(),
        reference_relations: reference.len(),
        false_positives,
        false_negatives,
        replay_failures,
        transform_mismatches,
        relation_sha256: sha256_hex(&bytes),
    }
}

fn mod_bigint(value: BigInt, modulus: &BigUint) -> BigUint {
    let modulus_i = BigInt::from_biguint(Sign::Plus, modulus.clone());
    let reduced = ((value % &modulus_i) + &modulus_i) % &modulus_i;
    reduced.to_biguint().expect("nonnegative modular value")
}

fn p256_control() -> Result<P256Control, String> {
    let curve = CurveParams::p256();
    let generator = Projective::from_textbook(&curve.generator());
    let mut labels = (1u32..=14).map(BigUint::from).collect::<Vec<BigUint>>();
    let sum = labels
        .iter()
        .fold(BigUint::zero(), |acc, value| acc + value);
    let alternating = labels
        .iter()
        .enumerate()
        .fold(BigInt::zero(), |acc, (index, value)| {
            let term = BigInt::from_biguint(Sign::Plus, value.clone());
            if index.is_multiple_of(2) {
                acc + term
            } else {
                acc - term
            }
        });
    let inverse_two = (&curve.n + BigUint::one()) >> 1usize;
    let twice_l0 = BigInt::from(2u8);
    let rhs_sum = -BigInt::from_biguint(Sign::Plus, sum);
    let rhs_alt = twice_l0 - alternating;
    let l14 = mod_bigint(
        (&rhs_sum + &rhs_alt) * BigInt::from_biguint(Sign::Plus, inverse_two.clone()),
        &curve.n,
    );
    let l15 = mod_bigint(
        (rhs_sum - rhs_alt) * BigInt::from_biguint(Sign::Plus, inverse_two),
        &curve.n,
    );
    labels.push(l14);
    labels.push(l15);
    if labels.iter().any(BigUint::is_zero) {
        return Err("P-256 Walsh control produced a zero label".into());
    }
    let mut points = labels
        .iter()
        .map(|label| generator.scalar_mul_ct(label, 256))
        .collect::<Vec<_>>();
    walsh_transform(&mut points, &mut WorkCounts::default());
    let normalized = batch_normalize(&points, &mut WorkCounts::default());
    let mut rows = BTreeMap::<Vec<i32>, RelationReceipt>::new();
    let mut counts = WorkCounts::default();
    let enumeration = (0..16).map(|index| Some(index as u32)).collect::<Vec<_>>();
    let base_points = labels
        .iter()
        .map(|label| generator.scalar_mul_ct(label, 256))
        .collect::<Vec<_>>();
    let mut identity_found = false;
    let mut double_found = false;
    for (character, canonical) in normalized.iter().enumerate() {
        if canonical.is_none() {
            identity_found |= character == 0;
            insert_relation(
                identity_row(character as u32, &enumeration),
                RelationReceipt {
                    source: "p256-control-identity".into(),
                    character_a: character as u32,
                    character_b: None,
                    column: None,
                    multiplier: None,
                    support: 0,
                    maximum_absolute_coefficient: 0,
                    row_sha256: String::new(),
                    replayed: false,
                },
                &mut rows,
                &base_points,
                &mut counts,
            )?;
        }
    }
    let two_p0 = base_points[0].double();
    let two_key = batch_normalize(&[two_p0], &mut counts)[0]
        .ok_or("P-256 control double became the identity")?;
    if let Some(output) = normalized[1] {
        if output.x == two_key.x {
            double_found = true;
            insert_relation(
                multiple_row(
                    1,
                    output.orientation,
                    MultipleRecord {
                        column: 0,
                        multiplier: 2,
                        orientation: two_key.orientation,
                    },
                    &enumeration,
                ),
                RelationReceipt {
                    source: "p256-control-double".into(),
                    character_a: 1,
                    character_b: None,
                    column: Some(0),
                    multiplier: Some(2),
                    support: 0,
                    maximum_absolute_coefficient: 0,
                    row_sha256: String::new(),
                    replayed: false,
                },
                &mut rows,
                &base_points,
                &mut counts,
            )?;
        }
    }
    let mut bytes = Vec::new();
    for row in rows.keys() {
        bytes.extend_from_slice(row_digest(row).as_bytes());
    }
    Ok(P256Control {
        fibre_size: 16,
        active_slots: 16,
        detected_relations: rows.len(),
        identity_character_zero_found: identity_found,
        character_one_double_match_found: double_found,
        replayed_relations: counts.replay_rows as usize,
        replay_failures: counts.replay_failures as usize,
        relation_sha256: sha256_hex(&bytes),
    })
}

fn gamma_ratio(collisions: u64) -> f64 {
    assert!(collisions >= 1);
    let mut ratio = std::f64::consts::PI.sqrt() / 2.0;
    for k in 1..collisions {
        ratio *= (k as f64 + 0.5) / k as f64;
    }
    ratio
}

fn run(cli: Cli) -> Result<(), String> {
    let (round2_receipt, round2) = load_dependency(&cli.round2, ROUND2_SHA256)?;
    let (round19_receipt, round19) = load_dependency(&cli.round19, ROUND19_SHA256)?;
    let (round312_receipt, round312) = load_dependency(&cli.round312, ROUND312_SHA256)?;
    for value in [&round2, &round19, &round312] {
        if value.get("curve").and_then(Value::as_str) != Some(CURVE_SLUG) {
            return Err("dependency curve identity changed".into());
        }
    }
    if round312
        .pointer("/factor_base/fb_id")
        .and_then(Value::as_str)
        != Some(FB_ID)
    {
        return Err("Round-312 factor-base identity changed".into());
    }

    let indexed = build_indexed(SPEC)?;
    let dump = &indexed.built.dump;
    if dump.curve.slug != CURVE_SLUG
        || dump.factor_base.fb_id != FB_ID
        || dump.factor_base.fb_sha256 != FB_SHA256
        || dump.factor_base.points_sha256 != POINTS_SHA256
        || dump.factor_base.columns as usize != COLUMNS
        || dump.points.len() != COLUMNS * 2
        || indexed.enumeration_to_column.len() != FIBRE_SIZE
    {
        return Err("native indexed factor-base rebuild changed".into());
    }
    let active_indices = indexed
        .enumeration_to_column
        .iter()
        .enumerate()
        .filter_map(|(index, column)| column.map(|_| index as u32))
        .collect::<Vec<_>>();
    let span_dimension = affine_span_dimension(&active_indices);
    if span_dimension != DEPTH as usize {
        return Err(format!(
            "active torus indices have affine dimension {span_dimension}"
        ));
    }

    let curve = CurveParams::p256();
    let columns = dump
        .points
        .chunks_exact(2)
        .enumerate()
        .map(|(column, rows)| {
            if rows[0].col as usize != column || rows[0].coef != "1" {
                return Err(format!(
                    "factor-base low-row convention changed at {column}"
                ));
            }
            low_point(&rows[0].x, &rows[0].y, &curve)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut values = indexed
        .enumeration_to_column
        .iter()
        .map(|column| column.map_or(Projective::IDENTITY, |column| columns[column as usize]))
        .collect::<Vec<_>>();
    let mut counts = WorkCounts::default();
    walsh_transform(&mut values, &mut counts);
    if counts.walsh_group_additions != DEPTH as u64 * FIBRE_SIZE as u64 {
        return Err("Walsh transform operation count changed".into());
    }
    let normalized = batch_normalize(&values, &mut counts);
    let mut outputs = Vec::with_capacity(FIBRE_SIZE);
    let mut identity_characters = Vec::new();
    for (character, canonical) in normalized.iter().enumerate() {
        match canonical {
            Some(canonical) => outputs.push(OutputRecord {
                x: canonical.x,
                orientation: canonical.orientation,
                character: character as u32,
            }),
            None => identity_characters.push(character as u32),
        }
    }
    outputs.sort_unstable_by(|left, right| {
        counts.canonical_key_comparisons += 1;
        (left.x, left.character).cmp(&(right.x, right.character))
    });

    let mut rows = BTreeMap::<Vec<i32>, RelationReceipt>::new();
    for character in &identity_characters {
        insert_relation(
            identity_row(*character, &indexed.enumeration_to_column),
            RelationReceipt {
                source: "walsh-identity".into(),
                character_a: *character,
                character_b: None,
                column: None,
                multiplier: None,
                support: 0,
                maximum_absolute_coefficient: 0,
                row_sha256: String::new(),
                replayed: false,
            },
            &mut rows,
            &columns,
            &mut counts,
        )?;
    }
    let mut start = 0usize;
    while start < outputs.len() {
        let key = outputs[start].x;
        let end = start
            + outputs[start..].partition_point(|record| {
                counts.canonical_key_comparisons += 1;
                record.x == key
            });
        for left in start..end {
            for right in left + 1..end {
                counts.output_collision_candidates += 1;
                insert_relation(
                    collision_row(
                        outputs[left].character,
                        outputs[left].orientation,
                        outputs[right].character,
                        outputs[right].orientation,
                        &indexed.enumeration_to_column,
                    ),
                    RelationReceipt {
                        source: "walsh-collision".into(),
                        character_a: outputs[left].character,
                        character_b: Some(outputs[right].character),
                        column: None,
                        multiplier: None,
                        support: 0,
                        maximum_absolute_coefficient: 0,
                        row_sha256: String::new(),
                        replayed: false,
                    },
                    &mut rows,
                    &columns,
                    &mut counts,
                )?;
            }
        }
        start = end;
    }

    let (multiple_points, multiple_meta) = build_multiples(&columns, &mut counts);
    let multiple_keys = batch_normalize(&multiple_points, &mut counts);
    for (canonical, (column, multiplier)) in multiple_keys.iter().zip(multiple_meta) {
        let Some(canonical) = canonical else {
            return Err(format!(
                "[{multiplier}] column {column} became the identity"
            ));
        };
        let begin = outputs.partition_point(|record| {
            counts.canonical_key_comparisons += 1;
            record.x < canonical.x
        });
        let end = outputs.partition_point(|record| {
            counts.canonical_key_comparisons += 1;
            record.x <= canonical.x
        });
        for output in &outputs[begin..end] {
            counts.output_multiple_candidates += 1;
            insert_relation(
                multiple_row(
                    output.character,
                    output.orientation,
                    MultipleRecord {
                        column,
                        multiplier,
                        orientation: canonical.orientation,
                    },
                    &indexed.enumeration_to_column,
                ),
                RelationReceipt {
                    source: "walsh-small-multiple".into(),
                    character_a: output.character,
                    character_b: None,
                    column: Some(column),
                    multiplier: Some(multiplier),
                    support: 0,
                    maximum_absolute_coefficient: 0,
                    row_sha256: String::new(),
                    replayed: false,
                },
                &mut rows,
                &columns,
                &mut counts,
            )?;
        }
    }

    let row_vectors = rows.keys().cloned().collect::<Vec<_>>();
    let rank = rank_receipt(&row_vectors, &mut counts);
    let cyclic = cyclic_control();
    let p256 = p256_control()?;
    if cyclic.false_positives != 0
        || cyclic.false_negatives != 0
        || cyclic.replay_failures != 0
        || cyclic.transform_mismatches != 0
        || !p256.identity_character_zero_found
        || !p256.character_one_double_match_found
        || p256.replay_failures != 0
        || counts.replay_failures != 0
    {
        return Err("Walsh correctness control failed".into());
    }

    let effective = COLUMNS
        .checked_sub(rank.promotion_credited_rank)
        .ok_or("rank exceeds factor-base columns")?;
    if effective == 0 {
        return Err("homogeneous Walsh rows cannot remove the final log scale".into());
    }
    let adjusted_rows =
        (RANK_ALLOWANCE_NUMERATOR * effective as u64).div_ceil(RANK_ALLOWANCE_DENOMINATOR);
    let free_s = 2f64.sqrt() * gamma_ratio(effective as u64);
    let per_relation = round19
        .pointer("/full_projection/optimistic_per_usable_relation_fme_log2")
        .and_then(Value::as_f64)
        .ok_or("Round-19 per-relation projection missing")?;
    let collection = per_relation + (adjusted_rows as f64).log2();
    let storage = round19
        .pointer("/full_projection/structured_frontier_materialized_bytes_log2")
        .and_then(Value::as_f64)
        .ok_or("Round-19 storage projection missing")?;
    let degree = round19
        .pointer("/residual_degree_evidence/structured_residual_maximum")
        .and_then(Value::as_u64)
        .ok_or("Round-19 degree evidence missing")?;
    let ecbench_eligible = rank.promotion_credited_rank != 0 && free_s <= RHO_S;
    let promoted = cyclic.false_positives == 0
        && cyclic.false_negatives == 0
        && counts.replay_failures == 0
        && degree <= 5
        && per_relation < 103.0
        && collection < 120.0
        && storage < 50.0
        && ecbench_eligible;

    let mut relation_bytes = Vec::new();
    for (row, receipt) in &rows {
        relation_bytes.extend_from_slice(receipt.row_sha256.as_bytes());
        relation_bytes.extend_from_slice(&(receipt.support as u64).to_be_bytes());
        relation_bytes.extend_from_slice(row_digest(row).as_bytes());
    }
    let dump_serialized_bytes = serde_json::to_vec(dump)
        .map_err(|error| format!("cannot size factor-base dump: {error}"))?
        .len() as u64;
    let row_materialized_bytes = rows
        .iter()
        .map(|(row, receipt)| {
            row.capacity() as u64 * std::mem::size_of::<i32>() as u64
                + std::mem::size_of_val(receipt) as u64
                + receipt.source.capacity() as u64
                + receipt.row_sha256.capacity() as u64
                + 256
        })
        .sum::<u64>();
    let logical_bytes = dump_serialized_bytes
        + values.capacity() as u64 * std::mem::size_of::<Projective>() as u64
        + normalized.capacity() as u64 * std::mem::size_of::<Option<CanonicalPoint>>() as u64
        + outputs.capacity() as u64 * std::mem::size_of::<OutputRecord>() as u64
        + identity_characters.capacity() as u64 * std::mem::size_of::<u32>() as u64
        + columns.capacity() as u64 * std::mem::size_of::<Projective>() as u64
        + multiple_points.capacity() as u64 * std::mem::size_of::<Projective>() as u64
        + multiple_keys.capacity() as u64 * std::mem::size_of::<Option<CanonicalPoint>>() as u64
        + COLUMNS as u64 * MULTIPLIERS.len() as u64 * std::mem::size_of::<(u32, u32)>() as u64
        + indexed.enumeration_to_column.capacity() as u64
            * std::mem::size_of::<Option<u32>>() as u64
        + active_indices.capacity() as u64 * std::mem::size_of::<u32>() as u64
        + row_materialized_bytes;
    let result = json!({
        "schema": "p256.dickson_walsh_transport/v1",
        "curve": CURVE_SLUG,
        "factor_base": {
            "spec": SPEC,
            "fb_id": dump.factor_base.fb_id,
            "fb_sha256": dump.factor_base.fb_sha256,
            "points_sha256": dump.factor_base.points_sha256,
            "columns": COLUMNS,
            "signed_points": dump.factor_base.signed_points,
            "complete_rebuild_verified": true,
            "preregistered_builder_source_sha256": PREREGISTERED_BUILDER_SHA256,
            "legacy_enumeration_sha256": indexed.built.legacy_enumeration_sha256,
            "enumeration_to_column_sha256": indexed.enumeration_to_column_sha256,
        },
        "dependencies": [round2_receipt, round19_receipt, round312_receipt],
        "walsh_family": {
            "depth": DEPTH,
            "fibre_slots": FIBRE_SIZE,
            "active_slots": active_indices.len(),
            "nonlifting_identity_slots": FIBRE_SIZE - active_indices.len(),
            "characters_exhausted": FIBRE_SIZE,
            "active_index_affine_span_dimension": span_dimension,
            "distinct_character_rows_certified": span_dimension == DEPTH as usize,
            "identity_outputs": identity_characters,
            "nonidentity_outputs": outputs.len(),
            "output_collision_candidates": counts.output_collision_candidates,
            "output_multiple_candidates": counts.output_multiple_candidates,
            "complete_scope": true,
            "discarded_characters": 0,
        },
        "relations": {
            "exact_unique_rows": rows.len(),
            "ordered_sha256": sha256_hex(&relation_bytes),
            "receipts": rows.values().collect::<Vec<_>>(),
            "every_reported_row_replayed": counts.replay_failures == 0,
        },
        "rank": rank,
        "controls": {
            "complete_cyclic": cyclic,
            "p256_scalar_labelled": p256,
        },
        "work": counts,
        "storage": {
            "conservative_logical_materialized_bytes": logical_bytes,
            "conservative_logical_materialized_log2_bytes": (logical_bytes as f64).log2(),
            "factor_base_serialized_bytes_component": dump_serialized_bytes,
            "exact_row_materialized_bytes_component": row_materialized_bytes,
            "disk_traffic_bytes_excluding_result_and_isolation": 0,
        },
        "boundary": {
            "unit": "P-256 group-addition equivalents / sqrt(n)",
            "rho_s": RHO_S,
            "original_columns": COLUMNS,
            "certified_static_rank": rank.promotion_credited_rank,
            "effective_log_classes": effective,
            "original_requested_rows": REQUESTED_ROWS,
            "adjusted_rows_with_five_percent_allowance": adjusted_rows,
            "free_oracle_s": free_s,
            "free_oracle_ratio_to_rho": free_s / RHO_S,
            "imported_per_usable_relation_fme_log2": per_relation,
            "projected_collection_fme_log2": collection,
            "imported_attack_storage_log2_bytes": storage,
            "complete_cold_s": Value::Null,
            "projection_is_measurement": false,
        },
        "ecbench": {
            "candidate_eligible": ecbench_eligible,
            "session_run": false,
            "reason": if ecbench_eligible {
                "nonzero rank reaches the preregistered free-oracle parity gate; a scaled native candidate is required next"
            } else {
                "the complete Walsh family does not reach the free-oracle parity gate, so no logarithm-recovering candidate exists for an ecbench session"
            },
            "p256_prime_path_supported": false,
            "p256_limit": "ecbench prime execution requires p < 2^62; projections are not session records",
        },
        "promotion_gates": {
            "zero_false_positives_and_false_negatives": cyclic.false_positives == 0 && cyclic.false_negatives == 0,
            "exact_p256_replay": counts.replay_failures == 0 && p256.replay_failures == 0,
            "structured_residual_degree_at_most_5": degree <= 5,
            "per_usable_relation_below_2_103": per_relation < 103.0,
            "collection_below_2_120": collection < 120.0,
            "storage_below_2_50": storage < 50.0,
            "complete_cold_s_at_most_rho": false,
            "ecbench_eligible": ecbench_eligible,
            "promoted": promoted,
        },
        "full_depth_unplanted_attempted": false,
        "decision": if promoted {
            "Walsh transport passed every preregistered gate; proceed to the separately preregistered ecbench and full-depth stages"
        } else {
            "Do not attempt a full-depth unplanted relation; preserve the complete Walsh result and leave parity unclaimed"
        },
    });
    let bytes = serde_json::to_vec_pretty(&result).map_err(|error| error.to_string())?;
    if let Some(parent) = cli.out.parent() {
        fs::create_dir_all(parent).map_err(|error| format!("{}: {error}", parent.display()))?;
    }
    fs::write(&cli.out, &bytes).map_err(|error| format!("{}: {error}", cli.out.display()))?;
    println!(
        "{}",
        serde_json::to_string_pretty(&result).map_err(|e| e.to_string())?
    );
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
    fn scalar_walsh_matches_direct_reference() {
        let labels = (0..16u64)
            .map(|i| (9 * i * i + 5) % 101)
            .collect::<Vec<_>>();
        assert_eq!(
            scalar_walsh(&labels, 101),
            direct_scalar_walsh(&labels, 101)
        );
    }

    #[test]
    fn complete_cyclic_control_has_no_discrepancy() {
        let control = cyclic_control();
        assert_eq!(control.transform_mismatches, 0);
        assert_eq!(control.false_positives, 0);
        assert_eq!(control.false_negatives, 0);
        assert_eq!(control.replay_failures, 0);
        assert!(control.optimized_relations > 0);
    }

    #[test]
    fn p256_planted_control_replays() {
        let control = p256_control().unwrap();
        assert!(control.identity_character_zero_found);
        assert!(control.character_one_double_match_found);
        assert_eq!(control.replay_failures, 0);
        assert!(control.detected_relations >= 2);
    }

    #[test]
    fn row_normalization_removes_scale_and_sign() {
        assert_eq!(normalize_row(vec![0, -4, 8, 0]), Some(vec![0, 1, -2, 0]));
        assert_eq!(normalize_row(vec![0, 0]), None);
    }

    #[test]
    fn affine_span_control_is_full() {
        let points = (0..1u32 << 5).collect::<Vec<_>>();
        assert_eq!(affine_span_dimension(&points), 5);
    }
}
