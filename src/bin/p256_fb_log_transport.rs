//! Exact small-multiple and Dickson-block logarithm-transport census for FB1.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::mem::{size_of, size_of_val};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{self, CURVE_SLUG};
use crypto_lib::ecc::p256_field::P256FieldElement as Fe;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::{CurveParams, Point};
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;

const SPEC: &str = "dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129";
const FB_ID: &str = "FB1h2f8621cda105";
const FB_SHA256: &str = "2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42";
const POINTS_SHA256: &str = "70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1";
const TERMINAL: &str = "0x5b17195299a3158b93389ad04c776fff2a8bb23ca8659b5b0c2b75f9009b65e5";
const COLUMNS: usize = 131_458;
const SIGNED_POINTS: u64 = 262_916;
const MAX_MULTIPLIER: u8 = 64;
const CONTROL_POINTS: usize = 4_096;
const REFERENCE_POINTS: usize = 12;
const BLOCK_DEPTHS: [u8; 4] = [1, 2, 3, 4];
const NORMALIZE_BATCH: usize = 8_192;
const MAX_BLOCK_IMAGES: u64 = 12_000_000;
const MAX_BLOCK_RECORD_BYTES: u64 = 1 << 30;
const ROUND20_SHA256: &str = "72f269c401d5192fda1e2a9f8c7e4b7bd537d92fe29a18b851882b610fa5c459";
const ROUND6_SHA256: &str = "71d63031111ba48ff831e79430bc87e6c7bac34626d7d0de4acfb4eaa65d40f4";
const DOMAIN: &str = concat!(
    "icv1-fp256-t89188191154553853111372247798585809583-f188c491/",
    "fb-log-transport-round21"
);

#[derive(Parser)]
#[command(about = "Exhaust exact logarithm transport on the registered P-256 FB1")]
struct Cli {
    /// Frozen round-20 collision result.
    #[arg(long)]
    round20: PathBuf,
    /// Frozen round-6 structured residual-degree result.
    #[arg(long)]
    round6: PathBuf,
    /// Deterministic JSON result path; stdout when omitted.
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Clone, Copy)]
struct FactorPoint {
    x: Fe,
    x_bytes: [u8; 32],
    positive: Projective,
}

#[derive(Clone, Debug)]
struct SparseRelation {
    coefficients: Vec<(u32, i32)>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct MultiplierRecord {
    x: [u8; 32],
    parity: u8,
    multiplier: u8,
    column: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct BlockRecord {
    x: [u8; 32],
    parity: u8,
    mask: u16,
    block: u32,
}

#[derive(Clone, Copy)]
struct RawBlockRecord {
    point: Projective,
    mask: u16,
    block: u32,
}

#[derive(Clone)]
struct BlockDef {
    ancestor: [u8; 32],
    columns: Vec<u32>,
}

#[derive(Default, Clone, Serialize)]
struct OperationCounts {
    multiplier_group_additions: u64,
    block_initial_group_additions: u64,
    block_gray_group_additions: u64,
    replay_group_additions: u64,
    replay_group_doublings: u64,
    normalization_field_multiplications: u64,
    dickson_field_squarings: u64,
}

#[derive(Serialize)]
struct DependencyReceipt {
    round20_sha256: String,
    round6_sha256: String,
}

#[derive(Serialize)]
struct FactorBaseReceipt {
    spec: String,
    fb_id: String,
    fb_sha256: String,
    points_sha256: String,
    terminal: String,
    columns: u64,
    signed_points: u64,
    complete_rebuild_verified: bool,
    dickson_leaves_checked: u64,
    dickson_leaf_rejections: u64,
}

#[derive(Serialize)]
struct ResidualDegreeEvidence {
    maximum_by_residual_depth: Vec<(u32, u32)>,
    every_component_complete_and_correct: bool,
    local_join_degree: u32,
    structured_residual_maximum: u32,
    unsplit_s17_degree_of_regularity: Option<u32>,
}

#[derive(Default, Serialize)]
struct MultiplierCensus {
    points: u64,
    max_multiplier: u8,
    images: u64,
    equal_x_buckets: u64,
    equal_x_pairs: u64,
    same_column_pairs: u64,
    cross_column_candidates: u64,
    replayed_relations: u64,
    replay_failures: u64,
    duplicate_relations: u64,
    independent_relation_rank: u64,
    quotient_components: u64,
    isolated_columns: u64,
    largest_component: u64,
    inconsistent_cycles: u64,
    adjacent_doubling_edges: u64,
    relation_digests: Vec<String>,
    relation_set_sha256: String,
    record_logical_bytes: u64,
    peak_logical_bytes: u64,
    disk_bytes_read: u64,
    disk_bytes_written: u64,
    operations: OperationCounts,
}

#[derive(Serialize)]
struct OrbitControl {
    points: u64,
    expected_components: u64,
    observed_components: u64,
    expected_adjacent_doubling_edges: u64,
    observed_adjacent_doubling_edges: u64,
    all_edges_replayed: bool,
    passed: bool,
    census: MultiplierCensus,
}

#[derive(Serialize)]
struct PrefixReference {
    selected_columns: Vec<u32>,
    selection_sha256: String,
    multiplier_sorted_relations: u64,
    multiplier_direct_relations: u64,
    multiplier_false_positives: u64,
    multiplier_false_negatives: u64,
    block_gray_images: Vec<(u8, u64)>,
    block_direct_images: Vec<(u8, u64)>,
    nontrivial_block_controls: Vec<BlockReferenceCell>,
    block_false_positives: u64,
    block_false_negatives: u64,
    passed: bool,
}

#[derive(Serialize)]
struct BlockReferenceCell {
    depth: u8,
    ancestor: String,
    leaves: u64,
    gray_images: u64,
    direct_images: u64,
    false_positives: u64,
    false_negatives: u64,
}

#[derive(Serialize)]
struct BlockDepthCensus {
    depth: u8,
    blocks: u64,
    eligible_blocks: u64,
    maximum_block_leaves: u64,
    registered_images: u64,
    materialized_images: u64,
    infinity_images: u64,
    equal_x_buckets: u64,
    cross_block_candidates: u64,
    factor_base_candidates: u64,
    replayed_relations: u64,
    replay_failures: u64,
    duplicate_relations: u64,
    unique_relations: u64,
    relation_set_sha256: String,
    completed: bool,
    skip_reason: Option<String>,
    record_logical_bytes: u64,
    peak_logical_bytes: u64,
    disk_bytes_read: u64,
    disk_bytes_written: u64,
    operations: OperationCounts,
}

#[derive(Serialize)]
struct QuotientReceipt {
    original_columns: u64,
    pairwise_relation_rank: u64,
    pairwise_components: u64,
    unique_block_relations: u64,
    independent_block_rank_in_pairwise_quotient: u64,
    total_independent_relation_rank: u64,
    quotient_dimension: u64,
    relation_set_sha256: String,
}

#[derive(Serialize)]
struct ProjectionRow {
    variant: String,
    class: String,
    group_addition_equivalents_log2: f64,
    s_total_over_sqrt_n: f64,
    ratio_to_rho: f64,
    complete_end_to_end: bool,
    correctness: String,
}

#[derive(Serialize)]
struct ParityProjection {
    group_order: String,
    quotient_dimension: u64,
    disjoint_probability: f64,
    rho_reference_s: f64,
    comparison_unit: String,
    rows: Vec<ProjectionRow>,
    two_list_parity_possible: bool,
    ideal_memoryless_parity_possible: bool,
    implemented_end_to_end_parity: bool,
}

#[derive(Serialize)]
struct ExperimentResult {
    schema: String,
    curve: String,
    dependencies: DependencyReceipt,
    factor_base: FactorBaseReceipt,
    orbit_control: OrbitControl,
    prefix_reference: PrefixReference,
    multiplier_census: MultiplierCensus,
    block_censuses: Vec<BlockDepthCensus>,
    quotient: QuotientReceipt,
    residual_degree_evidence: ResidualDegreeEvidence,
    parity_projection: ParityProjection,
    false_positives: u64,
    false_negatives: u64,
    exact_group_replay: bool,
    structured_degree_gate: bool,
    quotient_dimension_gate: bool,
    rho_parity_gate: bool,
    attack_promotion_gate: bool,
    full_depth_unplanted_attempted: bool,
    result_class: String,
    decision: String,
}

struct WeightedDsu {
    parent: Vec<usize>,
    size: Vec<u32>,
    /// `value[i] = weight[i] * value[parent[i]] (mod n)`.
    weight: Vec<BigUint>,
    modulus: BigUint,
    components: usize,
    independent_edges: u64,
    inconsistent_cycles: u64,
}

impl WeightedDsu {
    fn new(size: usize, modulus: &BigUint) -> Self {
        Self {
            parent: (0..size).collect(),
            size: vec![1; size],
            weight: vec![BigUint::one(); size],
            modulus: modulus.clone(),
            components: size,
            independent_edges: 0,
            inconsistent_cycles: 0,
        }
    }

    fn find(&mut self, item: usize) -> (usize, BigUint) {
        let parent = self.parent[item];
        if parent == item {
            return (item, BigUint::one());
        }
        let (root, upper) = self.find(parent);
        let combined = (&self.weight[item] * upper) % &self.modulus;
        self.parent[item] = root;
        self.weight[item] = combined.clone();
        (root, combined)
    }

    /// Add `value[left] = ratio * value[right] (mod n)`.
    fn unite(&mut self, left: usize, right: usize, ratio: &BigUint) -> bool {
        let (left_root, left_weight) = self.find(left);
        let (right_root, right_weight) = self.find(right);
        if left_root == right_root {
            if left_weight != (ratio * right_weight) % &self.modulus {
                self.inconsistent_cycles += 1;
                return false;
            }
            return true;
        }

        if self.size[left_root] <= self.size[right_root] {
            let inverse = mod_inverse(&left_weight, &self.modulus);
            let root_weight = (((ratio * right_weight) % &self.modulus) * inverse) % &self.modulus;
            self.parent[left_root] = right_root;
            self.weight[left_root] = root_weight;
            self.size[right_root] += self.size[left_root];
        } else {
            let right_factor = (ratio * right_weight) % &self.modulus;
            let inverse = mod_inverse(&right_factor, &self.modulus);
            let root_weight = (left_weight * inverse) % &self.modulus;
            self.parent[right_root] = left_root;
            self.weight[right_root] = root_weight;
            self.size[left_root] += self.size[right_root];
        }
        self.components -= 1;
        self.independent_edges += 1;
        true
    }

    fn histogram(&mut self) -> (usize, usize) {
        let mut sizes = BTreeMap::<usize, usize>::new();
        for item in 0..self.parent.len() {
            let (root, _) = self.find(item);
            *sizes.entry(root).or_default() += 1;
        }
        let isolated = sizes.values().filter(|size| **size == 1).count();
        let largest = sizes.values().copied().max().unwrap_or(0);
        (isolated, largest)
    }

    fn quotient_coordinates(&mut self) -> (Vec<usize>, Vec<BigUint>, usize) {
        let mut root_indices = BTreeMap::new();
        let mut variables = Vec::with_capacity(self.parent.len());
        let mut weights = Vec::with_capacity(self.parent.len());
        for item in 0..self.parent.len() {
            let (root, weight) = self.find(item);
            let next = root_indices.len();
            let index = *root_indices.entry(root).or_insert(next);
            variables.push(index);
            weights.push(weight);
        }
        (variables, weights, root_indices.len())
    }
}

#[derive(Default)]
struct SparseRank {
    pivots: BTreeMap<usize, BTreeMap<usize, BigUint>>,
}

impl SparseRank {
    fn insert(&mut self, mut row: BTreeMap<usize, BigUint>, modulus: &BigUint) -> bool {
        loop {
            row.retain(|_, value| !value.is_zero());
            let Some((&pivot, factor)) = row.first_key_value() else {
                return false;
            };
            let factor = factor.clone();
            if let Some(existing) = self.pivots.get(&pivot) {
                for (column, value) in existing {
                    let subtract = (&factor * value) % modulus;
                    let current = row.get(column).cloned().unwrap_or_default();
                    let updated = if current >= subtract {
                        current - &subtract
                    } else {
                        modulus - (&subtract - current)
                    };
                    if updated.is_zero() {
                        row.remove(column);
                    } else {
                        row.insert(*column, updated);
                    }
                }
                continue;
            }
            let inverse = mod_inverse(&factor, modulus);
            for value in row.values_mut() {
                *value = (&*value * &inverse) % modulus;
            }
            self.pivots.insert(pivot, row);
            return true;
        }
    }

    fn rank(&self) -> usize {
        self.pivots.len()
    }
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("not a non-negative integer: `{value}`"))
}

fn file_sha256(path: &Path) -> Result<(Vec<u8>, String), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    Ok((bytes, digest))
}

fn hash_strings(values: impl IntoIterator<Item = String>) -> String {
    let mut bytes = Vec::new();
    for value in values {
        bytes.extend_from_slice(value.as_bytes());
        bytes.push(b'\n');
    }
    hex::encode(sha256(&bytes))
}

fn verify_dependencies(round20: &Path, round6: &Path) -> Result<DependencyReceipt, String> {
    let (round20_bytes, round20_sha256) = file_sha256(round20)?;
    if round20_sha256 != ROUND20_SHA256 {
        return Err(format!(
            "round-20 SHA-256 mismatch: expected {ROUND20_SHA256}, got {round20_sha256}"
        ));
    }
    let round20_json: serde_json::Value =
        serde_json::from_slice(&round20_bytes).map_err(|error| error.to_string())?;
    if round20_json["schema"].as_str() != Some("p256.s17_batch_collision/v1")
        || round20_json["curve"].as_str() != Some(CURVE_SLUG)
        || round20_json["factor_base"]["fb_id"].as_str() != Some(FB_ID)
        || round20_json["attack_promotion_gate"].as_bool() != Some(false)
    {
        return Err("round-20 identity or rejection boundary changed".into());
    }

    let (round6_bytes, round6_sha256) = file_sha256(round6)?;
    if round6_sha256 != ROUND6_SHA256 {
        return Err(format!(
            "round-6 SHA-256 mismatch: expected {ROUND6_SHA256}, got {round6_sha256}"
        ));
    }
    let round6_json: serde_json::Value =
        serde_json::from_slice(&round6_bytes).map_err(|error| error.to_string())?;
    if round6_json["schema"].as_str() != Some("p256.dickson_residual_scaling/v1")
        || round6_json["curve"].as_str() != Some(CURVE_SLUG)
    {
        return Err("round-6 degree identity changed".into());
    }
    Ok(DependencyReceipt {
        round20_sha256,
        round6_sha256,
    })
}

fn load_degree_evidence(path: &Path) -> Result<ResidualDegreeEvidence, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let value: serde_json::Value =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    let cells = value["cells"]
        .as_array()
        .ok_or("round-6 cells are missing")?;
    let mut maxima = vec![(1u32, 0u32), (2, 0), (3, 0)];
    let mut complete = true;
    for cell in cells {
        let depth = cell["residual_depth"]
            .as_u64()
            .ok_or("round-6 residual depth is missing")? as u32;
        let degree = cell["max_solving_degree"]
            .as_u64()
            .ok_or("round-6 solving degree is missing")? as u32;
        if let Some((_, maximum)) = maxima
            .iter_mut()
            .find(|(residual_depth, _)| *residual_depth == depth)
        {
            *maximum = (*maximum).max(degree);
        }
        let components = cell["components"]
            .as_u64()
            .ok_or("round-6 component count is missing")?;
        complete &= cell["complete_components"].as_u64() == Some(components)
            && cell["correct_components"].as_u64() == Some(components);
    }
    if maxima != vec![(1, 3), (2, 3), (3, 4)] || !complete {
        return Err(format!(
            "round-6 residual-degree boundary changed: {maxima:?}"
        ));
    }
    Ok(ResidualDegreeEvidence {
        maximum_by_residual_depth: maxima,
        every_component_complete_and_correct: true,
        local_join_degree: 2,
        structured_residual_maximum: 4,
        unsplit_s17_degree_of_regularity: None,
    })
}

fn load_factor_base() -> Result<(FactorBaseReceipt, Vec<FactorPoint>), String> {
    let built = p256_dickson_factor_base::build(SPEC)?;
    p256_dickson_factor_base::verify(&built.dump)?;
    let fb = &built.dump.factor_base;
    if fb.fb_id != FB_ID
        || fb.fb_sha256 != FB_SHA256
        || fb.points_sha256 != POINTS_SHA256
        || fb.columns != COLUMNS as u64
        || fb.signed_points != SIGNED_POINTS
        || fb.params.get("terminal").map(String::as_str) != Some(TERMINAL)
    {
        return Err("factor-base identity changed".into());
    }
    let curve = CurveParams::p256();
    let terminal = Fe::from_biguint(&parse_big(TERMINAL)?);
    let two = Fe::ONE.add(&Fe::ONE);
    let mut points = Vec::with_capacity(COLUMNS);
    for column in 0..COLUMNS {
        let row = &built.dump.points[2 * column];
        if row.col as usize != column || row.coef != "1" {
            return Err(format!(
                "factor-base low-row convention changed at {column}"
            ));
        }
        let x_big = parse_big(&row.x)?;
        let y_big = parse_big(&row.y)?;
        let affine = Point::Affine {
            x: curve.fe(x_big.clone()),
            y: curve.fe(y_big),
        };
        if !curve.is_on_curve(&affine) {
            return Err(format!("factor-base column {column} is off curve"));
        }
        let x = Fe::from_biguint(&x_big);
        let mut chain = x;
        for _ in 0..18 {
            chain = chain.sqr().sub(&two);
        }
        if chain != terminal {
            return Err(format!(
                "factor-base column {column} misses Dickson terminal"
            ));
        }
        points.push(FactorPoint {
            x,
            x_bytes: x.to_bytes_be(),
            positive: Projective::from_textbook(&affine),
        });
    }
    Ok((
        FactorBaseReceipt {
            spec: SPEC.into(),
            fb_id: fb.fb_id.clone(),
            fb_sha256: fb.fb_sha256.clone(),
            points_sha256: fb.points_sha256.clone(),
            terminal: TERMINAL.into(),
            columns: fb.columns,
            signed_points: fb.signed_points,
            complete_rebuild_verified: true,
            dickson_leaves_checked: COLUMNS as u64,
            dickson_leaf_rejections: 0,
        },
        points,
    ))
}

fn mod_inverse(value: &BigUint, modulus: &BigUint) -> BigUint {
    debug_assert!(!value.is_zero());
    value.modpow(&(modulus - 2u8), modulus)
}

fn pow_counted(mut base: Fe, exponent: &BigUint, multiplications: &mut u64) -> Fe {
    let mut result = Fe::ONE;
    for bit in 0..exponent.bits() {
        if exponent.bit(bit) {
            result = result.mul(&base);
            *multiplications += 1;
        }
        base = base.sqr();
        *multiplications += 1;
    }
    result
}

fn normalize_points(
    points: &[Projective],
    inverse_exponent: &BigUint,
    multiplications: &mut u64,
) -> Result<Vec<([u8; 32], u8)>, String> {
    let mut affine = Vec::with_capacity(points.len());
    for batch in points.chunks(NORMALIZE_BATCH) {
        if batch.iter().any(|point| bool::from(point.is_identity())) {
            return Err("normalization input contains infinity".into());
        }
        let mut prefixes = Vec::with_capacity(batch.len());
        let mut product = Fe::ONE;
        for point in batch {
            prefixes.push(product);
            product = product.mul(&point.z);
            *multiplications += 1;
        }
        let mut inverse = pow_counted(product, inverse_exponent, multiplications);
        let mut inverses = vec![Fe::ZERO; batch.len()];
        for index in (0..batch.len()).rev() {
            inverses[index] = inverse.mul(&prefixes[index]);
            *multiplications += 1;
            if index != 0 {
                inverse = inverse.mul(&batch[index].z);
                *multiplications += 1;
            }
        }
        for (point, inverse) in batch.iter().zip(inverses) {
            let x = point.x.mul(&inverse);
            let y = point.y.mul(&inverse);
            *multiplications += 2;
            affine.push((x.to_bytes_be(), y.to_bytes_be()[31] & 1));
        }
    }
    Ok(affine)
}

fn small_scalar_mul(
    point: &Projective,
    scalar: u32,
    operations: &mut OperationCounts,
) -> Projective {
    let mut result = Projective::IDENTITY;
    let mut addend = *point;
    let mut value = scalar;
    while value != 0 {
        if value & 1 == 1 {
            result = result.add(&addend);
            operations.replay_group_additions += 1;
        }
        value >>= 1;
        if value != 0 {
            addend = addend.double();
            operations.replay_group_doublings += 1;
        }
    }
    result
}

fn multiplier_record_cmp(left: &MultiplierRecord, right: &MultiplierRecord) -> Ordering {
    left.x
        .cmp(&right.x)
        .then_with(|| left.parity.cmp(&right.parity))
        .then_with(|| left.column.cmp(&right.column))
        .then_with(|| left.multiplier.cmp(&right.multiplier))
}

fn build_multiplier_records(
    points: &[FactorPoint],
    max_multiplier: u8,
    inverse_exponent: &BigUint,
    operations: &mut OperationCounts,
) -> Result<Vec<MultiplierRecord>, String> {
    let mut current = points
        .iter()
        .map(|point| point.positive)
        .collect::<Vec<_>>();
    let mut records = Vec::with_capacity(points.len() * max_multiplier as usize);
    for multiplier in 1..=max_multiplier {
        let affine = normalize_points(
            &current,
            inverse_exponent,
            &mut operations.normalization_field_multiplications,
        )?;
        records.extend(affine.into_iter().enumerate().map(|(column, (x, parity))| {
            MultiplierRecord {
                x,
                parity,
                multiplier,
                column: column as u32,
            }
        }));
        if multiplier != max_multiplier {
            for (image, base) in current.iter_mut().zip(points) {
                *image = image.add(&base.positive);
                operations.multiplier_group_additions += 1;
            }
        }
    }
    records.sort_unstable_by(multiplier_record_cmp);
    Ok(records)
}

fn canonical_relation(mut coefficients: BTreeMap<u32, i32>) -> Option<SparseRelation> {
    coefficients.retain(|_, coefficient| *coefficient != 0);
    if coefficients.is_empty() {
        return None;
    }
    let gcd = coefficients
        .values()
        .map(|coefficient| coefficient.unsigned_abs())
        .reduce(gcd_u32)
        .unwrap_or(1)
        .max(1) as i32;
    for coefficient in coefficients.values_mut() {
        *coefficient /= gcd;
    }
    if coefficients
        .first_key_value()
        .is_some_and(|(_, value)| *value < 0)
    {
        for coefficient in coefficients.values_mut() {
            *coefficient = -*coefficient;
        }
    }
    Some(SparseRelation {
        coefficients: coefficients.into_iter().collect(),
    })
}

fn gcd_u32(mut left: u32, mut right: u32) -> u32 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

fn relation_digest(relation: &SparseRelation) -> String {
    let mut bytes = format!("{DOMAIN}/sparse-relation/v1\0").into_bytes();
    for (column, coefficient) in &relation.coefficients {
        bytes.extend_from_slice(&column.to_be_bytes());
        bytes.extend_from_slice(&coefficient.to_be_bytes());
    }
    hex::encode(sha256(&bytes))
}

fn replay_sparse_relation(
    relation: &SparseRelation,
    points: &[FactorPoint],
    operations: &mut OperationCounts,
) -> bool {
    let mut total = Projective::IDENTITY;
    for (column, coefficient) in &relation.coefficients {
        let magnitude = coefficient.unsigned_abs();
        let mut term = small_scalar_mul(&points[*column as usize].positive, magnitude, operations);
        if *coefficient < 0 {
            term = term.neg();
        }
        total = total.add(&term);
        operations.replay_group_additions += 1;
    }
    bool::from(total.is_identity())
}

fn multiplier_relation(left: &MultiplierRecord, right: &MultiplierRecord) -> SparseRelation {
    let orientation = if left.parity == right.parity { -1 } else { 1 };
    let mut coefficients = BTreeMap::new();
    *coefficients.entry(left.column).or_default() += i32::from(left.multiplier);
    *coefficients.entry(right.column).or_default() += orientation * i32::from(right.multiplier);
    canonical_relation(coefficients).expect("distinct columns make a nonempty relation")
}

fn multiplier_ratio(
    left: &MultiplierRecord,
    right: &MultiplierRecord,
    modulus: &BigUint,
) -> BigUint {
    let left_inverse = mod_inverse(&BigUint::from(left.multiplier), modulus);
    let mut numerator = BigUint::from(right.multiplier);
    if left.parity != right.parity {
        numerator = modulus - numerator;
    }
    (numerator * left_inverse) % modulus
}

fn analyze_multiplier_records(
    records: &[MultiplierRecord],
    points: &[FactorPoint],
    max_multiplier: u8,
    modulus: &BigUint,
    mut operations: OperationCounts,
) -> Result<(MultiplierCensus, WeightedDsu), String> {
    let mut dsu = WeightedDsu::new(points.len(), modulus);
    let mut report = MultiplierCensus {
        points: points.len() as u64,
        max_multiplier,
        images: records.len() as u64,
        record_logical_bytes: size_of::<MultiplierRecord>() as u64,
        peak_logical_bytes: (size_of_val(records) + size_of_val(points) * 2) as u64,
        ..MultiplierCensus::default()
    };
    let mut unique = BTreeSet::new();
    let mut adjacent = BTreeSet::new();
    let mut start = 0;
    while start < records.len() {
        let mut end = start + 1;
        while end < records.len() && records[end].x == records[start].x {
            end += 1;
        }
        if end - start > 1 {
            report.equal_x_buckets += 1;
        }
        for left_index in start..end {
            for right_index in left_index + 1..end {
                report.equal_x_pairs += 1;
                let left = &records[left_index];
                let right = &records[right_index];
                if left.column == right.column {
                    report.same_column_pairs += 1;
                    continue;
                }
                report.cross_column_candidates += 1;
                let relation = multiplier_relation(left, right);
                if !replay_sparse_relation(&relation, points, &mut operations) {
                    report.replay_failures += 1;
                    return Err(format!(
                        "small-multiple replay failed: {}",
                        relation_digest(&relation)
                    ));
                }
                report.replayed_relations += 1;
                let digest = relation_digest(&relation);
                if !unique.insert(digest) {
                    report.duplicate_relations += 1;
                }
                let ratio = multiplier_ratio(left, right, modulus);
                if !dsu.unite(left.column as usize, right.column as usize, &ratio) {
                    return Err("inconsistent exact small-multiple cycle".into());
                }
                for (lower, upper) in [(left, right), (right, left)] {
                    if lower.multiplier == 2
                        && upper.multiplier == 1
                        && upper.column == lower.column + 1
                    {
                        adjacent.insert(lower.column);
                    }
                }
            }
        }
        start = end;
    }
    let (isolated, largest) = dsu.histogram();
    report.independent_relation_rank = dsu.independent_edges;
    report.quotient_components = dsu.components as u64;
    report.isolated_columns = isolated as u64;
    report.largest_component = largest as u64;
    report.inconsistent_cycles = dsu.inconsistent_cycles;
    report.adjacent_doubling_edges = adjacent.len() as u64;
    report.relation_digests = unique.into_iter().collect();
    report.relation_set_sha256 = hash_strings(report.relation_digests.iter().cloned());
    report.operations = operations;
    Ok((report, dsu))
}

fn multiplier_census(
    points: &[FactorPoint],
    max_multiplier: u8,
    inverse_exponent: &BigUint,
    modulus: &BigUint,
) -> Result<(MultiplierCensus, WeightedDsu, Vec<MultiplierRecord>), String> {
    let mut operations = OperationCounts::default();
    let records =
        build_multiplier_records(points, max_multiplier, inverse_exponent, &mut operations)?;
    let (report, dsu) =
        analyze_multiplier_records(&records, points, max_multiplier, modulus, operations)?;
    Ok((report, dsu, records))
}

fn orbit_control(inverse_exponent: &BigUint, modulus: &BigUint) -> Result<OrbitControl, String> {
    let curve = CurveParams::p256();
    let mut point = Projective::from_textbook(&curve.generator());
    let mut points = Vec::with_capacity(CONTROL_POINTS);
    for _ in 0..CONTROL_POINTS {
        let (x, _) = point.to_affine().ok_or("orbit control reached infinity")?;
        let x = Fe::from_biguint(&x.to_biguint());
        points.push(FactorPoint {
            x,
            x_bytes: x.to_bytes_be(),
            positive: point,
        });
        point = point.double();
    }
    let (mut census, _, _) = multiplier_census(&points, MAX_MULTIPLIER, inverse_exponent, modulus)?;
    let relation_set_sha256 = census.relation_set_sha256.clone();
    census.relation_digests.clear();
    census.relation_set_sha256 = relation_set_sha256;
    let passed = census.quotient_components == 1
        && census.adjacent_doubling_edges == (CONTROL_POINTS - 1) as u64
        && census.replay_failures == 0
        && census.inconsistent_cycles == 0;
    Ok(OrbitControl {
        points: CONTROL_POINTS as u64,
        expected_components: 1,
        observed_components: census.quotient_components,
        expected_adjacent_doubling_edges: (CONTROL_POINTS - 1) as u64,
        observed_adjacent_doubling_edges: census.adjacent_doubling_edges,
        all_edges_replayed: census.replay_failures == 0,
        passed,
        census,
    })
}

fn select_reference_columns() -> (Vec<u32>, String) {
    let mut keyed = (0..COLUMNS as u32)
        .map(|column| {
            let digest = sha256(format!("{DOMAIN}/reference/{column}").as_bytes());
            (digest, column)
        })
        .collect::<Vec<_>>();
    keyed.sort_unstable();
    let mut selected = keyed
        .into_iter()
        .take(REFERENCE_POINTS)
        .map(|(_, column)| column)
        .collect::<Vec<_>>();
    selected.sort_unstable();
    let mut bytes = Vec::new();
    for column in &selected {
        bytes.extend_from_slice(&column.to_be_bytes());
    }
    (selected, hex::encode(sha256(&bytes)))
}

fn direct_multiplier_relation_digests(records: &[MultiplierRecord]) -> BTreeSet<String> {
    let mut digests = BTreeSet::new();
    for left_index in 0..records.len() {
        for right_index in left_index + 1..records.len() {
            let left = &records[left_index];
            let right = &records[right_index];
            if left.column != right.column && left.x == right.x {
                digests.insert(relation_digest(&multiplier_relation(left, right)));
            }
        }
    }
    digests
}

fn block_defs(
    points: &[FactorPoint],
    columns: &[u32],
    depth: u8,
    operations: &mut OperationCounts,
) -> Vec<BlockDef> {
    let two = Fe::ONE.add(&Fe::ONE);
    let mut grouped = BTreeMap::<[u8; 32], Vec<u32>>::new();
    for column in columns {
        let mut ancestor = points[*column as usize].x;
        for _ in 0..depth {
            ancestor = ancestor.sqr().sub(&two);
            operations.dickson_field_squarings += 1;
        }
        grouped
            .entry(ancestor.to_bytes_be())
            .or_default()
            .push(*column);
    }
    grouped
        .into_iter()
        .map(|(ancestor, columns)| BlockDef { ancestor, columns })
        .collect()
}

fn signed_sum_direct(block: &BlockDef, mask: u16, points: &[FactorPoint]) -> Projective {
    let mut sum = points[block.columns[0] as usize].positive;
    for (position, column) in block.columns.iter().copied().enumerate().skip(1) {
        let mut point = points[column as usize].positive;
        if mask & (1 << (position - 1)) != 0 {
            point = point.neg();
        }
        sum = sum.add(&point);
    }
    sum
}

fn block_relation(record: &BlockRecord, block: &BlockDef) -> SparseRelation {
    let mut coefficients = BTreeMap::new();
    for (position, column) in block.columns.iter().copied().enumerate() {
        let coefficient = if position == 0 || record.mask & (1 << (position - 1)) == 0 {
            1
        } else {
            -1
        };
        *coefficients.entry(column).or_default() += coefficient;
    }
    canonical_relation(coefficients).expect("a block image has nonempty support")
}

fn equality_relation(
    left: &BlockRecord,
    left_block: &BlockDef,
    right: &BlockRecord,
    right_block: &BlockDef,
) -> Option<SparseRelation> {
    let left_relation = block_relation(left, left_block);
    let right_relation = block_relation(right, right_block);
    let orientation = if left.parity == right.parity { -1 } else { 1 };
    let mut coefficients = BTreeMap::new();
    for (column, coefficient) in left_relation.coefficients {
        *coefficients.entry(column).or_default() += coefficient;
    }
    for (column, coefficient) in right_relation.coefficients {
        *coefficients.entry(column).or_default() += orientation * coefficient;
    }
    canonical_relation(coefficients)
}

fn base_hit_relation(
    record: &BlockRecord,
    block: &BlockDef,
    base_column: u32,
    base_parity: u8,
) -> Option<SparseRelation> {
    let block_relation = block_relation(record, block);
    let orientation = if record.parity == base_parity { -1 } else { 1 };
    let mut coefficients = block_relation
        .coefficients
        .into_iter()
        .collect::<BTreeMap<_, _>>();
    *coefficients.entry(base_column).or_default() += orientation;
    canonical_relation(coefficients)
}

fn flush_block_raw(
    raw: &mut Vec<RawBlockRecord>,
    records: &mut Vec<BlockRecord>,
    inverse_exponent: &BigUint,
    operations: &mut OperationCounts,
) -> Result<(), String> {
    if raw.is_empty() {
        return Ok(());
    }
    let points = raw.iter().map(|record| record.point).collect::<Vec<_>>();
    let affine = normalize_points(
        &points,
        inverse_exponent,
        &mut operations.normalization_field_multiplications,
    )?;
    records.extend(
        raw.iter()
            .zip(affine)
            .map(|(raw, (x, parity))| BlockRecord {
                x,
                parity,
                mask: raw.mask,
                block: raw.block,
            }),
    );
    raw.clear();
    Ok(())
}

fn enumerate_block_records(
    blocks: &[BlockDef],
    points: &[FactorPoint],
    inverse_exponent: &BigUint,
    operations: &mut OperationCounts,
    infinity_relations: &mut Vec<SparseRelation>,
) -> Result<Vec<BlockRecord>, String> {
    let capacity: usize = blocks
        .iter()
        .filter(|block| block.columns.len() >= 2)
        .map(|block| 1usize << (block.columns.len() - 1))
        .sum();
    let mut records = Vec::with_capacity(capacity);
    let mut raw = Vec::with_capacity(NORMALIZE_BATCH);
    for (block_index, block) in blocks.iter().enumerate() {
        let leaves = block.columns.len();
        if leaves < 2 {
            continue;
        }
        if leaves > 16 {
            return Err(format!(
                "depth block {} has {leaves} leaves; mask encoding only permits 16",
                hex::encode(block.ancestor)
            ));
        }
        let mut current = points[block.columns[0] as usize].positive;
        let mut doubles = Vec::with_capacity(leaves - 1);
        for column in block.columns.iter().copied().skip(1) {
            let point = points[column as usize].positive;
            current = current.add(&point);
            operations.block_initial_group_additions += 1;
            doubles.push(point.double());
        }
        let count = 1usize << (leaves - 1);
        let mut previous_gray = 0u16;
        for ordinal in 0..count {
            let gray = (ordinal ^ (ordinal >> 1)) as u16;
            if ordinal != 0 {
                let changed = previous_gray ^ gray;
                let bit = changed.trailing_zeros() as usize;
                let delta = if gray & changed == 0 {
                    doubles[bit]
                } else {
                    doubles[bit].neg()
                };
                current = current.add(&delta);
                operations.block_gray_group_additions += 1;
            }
            if bool::from(current.is_identity()) {
                let synthetic = BlockRecord {
                    x: [0; 32],
                    parity: 0,
                    mask: gray,
                    block: block_index as u32,
                };
                infinity_relations.push(block_relation(&synthetic, block));
            } else {
                raw.push(RawBlockRecord {
                    point: current,
                    mask: gray,
                    block: block_index as u32,
                });
                if raw.len() == NORMALIZE_BATCH {
                    flush_block_raw(&mut raw, &mut records, inverse_exponent, operations)?;
                }
            }
            previous_gray = gray;
        }
    }
    flush_block_raw(&mut raw, &mut records, inverse_exponent, operations)?;
    records.sort_unstable_by(|left, right| {
        left.x
            .cmp(&right.x)
            .then_with(|| left.parity.cmp(&right.parity))
            .then_with(|| left.block.cmp(&right.block))
            .then_with(|| left.mask.cmp(&right.mask))
    });
    Ok(records)
}

fn block_reference_keys(
    blocks: &[BlockDef],
    points: &[FactorPoint],
    gray: bool,
) -> Result<Vec<([u8; 32], u8, u32, u16)>, String> {
    let mut out = Vec::new();
    for (block_index, block) in blocks.iter().enumerate() {
        if block.columns.len() < 2 {
            continue;
        }
        let count = 1usize << (block.columns.len() - 1);
        if gray {
            let mut current = signed_sum_direct(block, 0, points);
            let doubles = block
                .columns
                .iter()
                .copied()
                .skip(1)
                .map(|column| points[column as usize].positive.double())
                .collect::<Vec<_>>();
            let mut previous = 0u16;
            for ordinal in 0..count {
                let mask = (ordinal ^ (ordinal >> 1)) as u16;
                if ordinal != 0 {
                    let changed = previous ^ mask;
                    let bit = changed.trailing_zeros() as usize;
                    current = current.add(&if mask & changed == 0 {
                        doubles[bit]
                    } else {
                        doubles[bit].neg()
                    });
                }
                let (x, y) = current
                    .to_affine()
                    .ok_or("reference block sum is infinity")?;
                out.push((
                    x.to_bytes_be(),
                    y.to_bytes_be()[31] & 1,
                    block_index as u32,
                    mask,
                ));
                previous = mask;
            }
        } else {
            for mask in 0..count as u16 {
                let point = signed_sum_direct(block, mask, points);
                let (x, y) = point.to_affine().ok_or("direct block sum is infinity")?;
                out.push((
                    x.to_bytes_be(),
                    y.to_bytes_be()[31] & 1,
                    block_index as u32,
                    mask,
                ));
            }
        }
    }
    out.sort_unstable();
    Ok(out)
}

fn prefix_reference(
    points: &[FactorPoint],
    inverse_exponent: &BigUint,
    modulus: &BigUint,
) -> Result<PrefixReference, String> {
    let (selected_columns, selection_sha256) = select_reference_columns();
    let selected_points = selected_columns
        .iter()
        .map(|column| points[*column as usize])
        .collect::<Vec<_>>();
    let mut operations = OperationCounts::default();
    let records = build_multiplier_records(
        &selected_points,
        MAX_MULTIPLIER,
        inverse_exponent,
        &mut operations,
    )?;
    let (sorted, _, _) =
        multiplier_census(&selected_points, MAX_MULTIPLIER, inverse_exponent, modulus)?;
    let sorted_set = sorted
        .relation_digests
        .iter()
        .cloned()
        .collect::<BTreeSet<_>>();
    let direct_set = direct_multiplier_relation_digests(&records);
    let multiplier_false_positives = sorted_set.difference(&direct_set).count() as u64;
    let multiplier_false_negatives = direct_set.difference(&sorted_set).count() as u64;

    let mut block_gray_images = Vec::new();
    let mut block_direct_images = Vec::new();
    let mut nontrivial_block_controls = Vec::new();
    let mut block_false_positives = 0;
    let mut block_false_negatives = 0;
    for depth in BLOCK_DEPTHS {
        let mut ignored = OperationCounts::default();
        let blocks = block_defs(points, &selected_columns, depth, &mut ignored);
        let gray = block_reference_keys(&blocks, points, true)?;
        let direct = block_reference_keys(&blocks, points, false)?;
        let gray_set = gray.iter().copied().collect::<BTreeSet<_>>();
        let direct_set = direct.iter().copied().collect::<BTreeSet<_>>();
        block_false_positives += gray_set.difference(&direct_set).count() as u64;
        block_false_negatives += direct_set.difference(&gray_set).count() as u64;
        block_gray_images.push((depth, gray.len() as u64));
        block_direct_images.push((depth, direct.len() as u64));

        let all_columns = (0..points.len() as u32).collect::<Vec<_>>();
        let full_blocks = block_defs(points, &all_columns, depth, &mut ignored);
        let maximum = full_blocks
            .iter()
            .map(|block| block.columns.len())
            .max()
            .unwrap_or(0);
        let selected = full_blocks
            .iter()
            .filter(|block| block.columns.len() == maximum)
            .min_by_key(|block| {
                let mut bytes = format!("{DOMAIN}/block-reference/{depth}/").into_bytes();
                bytes.extend_from_slice(&block.ancestor);
                sha256(&bytes)
            })
            .ok_or("no block available for nontrivial reference")?
            .clone();
        let selected_blocks = vec![selected.clone()];
        let gray = block_reference_keys(&selected_blocks, points, true)?;
        let direct = block_reference_keys(&selected_blocks, points, false)?;
        let gray_set = gray.iter().copied().collect::<BTreeSet<_>>();
        let direct_set = direct.iter().copied().collect::<BTreeSet<_>>();
        let false_positives = gray_set.difference(&direct_set).count() as u64;
        let false_negatives = direct_set.difference(&gray_set).count() as u64;
        block_false_positives += false_positives;
        block_false_negatives += false_negatives;
        nontrivial_block_controls.push(BlockReferenceCell {
            depth,
            ancestor: hex::encode(selected.ancestor),
            leaves: selected.columns.len() as u64,
            gray_images: gray.len() as u64,
            direct_images: direct.len() as u64,
            false_positives,
            false_negatives,
        });
    }
    let passed = multiplier_false_positives == 0
        && multiplier_false_negatives == 0
        && block_false_positives == 0
        && block_false_negatives == 0;
    Ok(PrefixReference {
        selected_columns,
        selection_sha256,
        multiplier_sorted_relations: sorted_set.len() as u64,
        multiplier_direct_relations: direct_set.len() as u64,
        multiplier_false_positives,
        multiplier_false_negatives,
        block_gray_images,
        block_direct_images,
        nontrivial_block_controls,
        block_false_positives,
        block_false_negatives,
        passed,
    })
}

fn block_census(
    depth: u8,
    points: &[FactorPoint],
    inverse_exponent: &BigUint,
    all_relations: &mut BTreeMap<String, SparseRelation>,
) -> Result<BlockDepthCensus, String> {
    let mut operations = OperationCounts::default();
    let columns = (0..points.len() as u32).collect::<Vec<_>>();
    let blocks = block_defs(points, &columns, depth, &mut operations);
    let eligible = blocks
        .iter()
        .filter(|block| block.columns.len() >= 2)
        .count();
    let maximum = blocks
        .iter()
        .map(|block| block.columns.len())
        .max()
        .unwrap_or(0);
    let registered_images = blocks
        .iter()
        .filter(|block| block.columns.len() >= 2)
        .map(|block| 1u64 << (block.columns.len() - 1))
        .sum::<u64>();
    let record_bytes = size_of::<BlockRecord>() as u64;
    if registered_images > MAX_BLOCK_IMAGES
        || registered_images.saturating_mul(record_bytes) > MAX_BLOCK_RECORD_BYTES
    {
        return Ok(BlockDepthCensus {
            depth,
            blocks: blocks.len() as u64,
            eligible_blocks: eligible as u64,
            maximum_block_leaves: maximum as u64,
            registered_images,
            materialized_images: 0,
            infinity_images: 0,
            equal_x_buckets: 0,
            cross_block_candidates: 0,
            factor_base_candidates: 0,
            replayed_relations: 0,
            replay_failures: 0,
            duplicate_relations: 0,
            unique_relations: 0,
            relation_set_sha256: hash_strings(Vec::<String>::new()),
            completed: false,
            skip_reason: Some("registered image or logical-byte cap exceeded".into()),
            record_logical_bytes: record_bytes,
            peak_logical_bytes: 0,
            disk_bytes_read: 0,
            disk_bytes_written: 0,
            operations,
        });
    }

    let mut infinity_relations = Vec::new();
    let records = enumerate_block_records(
        &blocks,
        points,
        inverse_exponent,
        &mut operations,
        &mut infinity_relations,
    )?;
    let mut unique = BTreeMap::<String, SparseRelation>::new();
    let mut replayed = 0u64;
    let mut duplicates = 0u64;
    for relation in infinity_relations.iter().cloned() {
        if !replay_sparse_relation(&relation, points, &mut operations) {
            return Err("infinity block relation failed replay".into());
        }
        replayed += 1;
        let digest = relation_digest(&relation);
        if unique.insert(digest, relation).is_some() {
            duplicates += 1;
        }
    }

    let mut equal_x_buckets = 0u64;
    let mut cross_block_candidates = 0u64;
    let mut start = 0;
    while start < records.len() {
        let mut end = start + 1;
        while end < records.len() && records[end].x == records[start].x {
            end += 1;
        }
        if end - start > 1 {
            equal_x_buckets += 1;
        }
        for left_index in start..end {
            for right_index in left_index + 1..end {
                let left = &records[left_index];
                let right = &records[right_index];
                if left.block == right.block {
                    continue;
                }
                cross_block_candidates += 1;
                let Some(relation) = equality_relation(
                    left,
                    &blocks[left.block as usize],
                    right,
                    &blocks[right.block as usize],
                ) else {
                    continue;
                };
                if !replay_sparse_relation(&relation, points, &mut operations) {
                    return Err(format!(
                        "cross-block replay failed: {}",
                        relation_digest(&relation)
                    ));
                }
                replayed += 1;
                let digest = relation_digest(&relation);
                if unique.insert(digest, relation).is_some() {
                    duplicates += 1;
                }
            }
        }
        start = end;
    }

    let mut factor_base_candidates = 0u64;
    for record in &records {
        if let Ok(base_column) = points.binary_search_by(|point| point.x_bytes.cmp(&record.x)) {
            factor_base_candidates += 1;
            let (_, base_y) = points[base_column]
                .positive
                .to_affine()
                .ok_or("factor-base point is infinity")?;
            if let Some(relation) = base_hit_relation(
                record,
                &blocks[record.block as usize],
                base_column as u32,
                base_y.to_bytes_be()[31] & 1,
            ) {
                if !replay_sparse_relation(&relation, points, &mut operations) {
                    return Err(format!(
                        "factor-base-hit replay failed: {}",
                        relation_digest(&relation)
                    ));
                }
                replayed += 1;
                let digest = relation_digest(&relation);
                if unique.insert(digest, relation).is_some() {
                    duplicates += 1;
                }
            }
        }
    }

    for (digest, relation) in &unique {
        all_relations
            .entry(digest.clone())
            .or_insert_with(|| relation.clone());
    }
    let relation_set_sha256 = hash_strings(unique.keys().cloned());
    Ok(BlockDepthCensus {
        depth,
        blocks: blocks.len() as u64,
        eligible_blocks: eligible as u64,
        maximum_block_leaves: maximum as u64,
        registered_images,
        materialized_images: records.len() as u64,
        infinity_images: infinity_relations.len() as u64,
        equal_x_buckets,
        cross_block_candidates,
        factor_base_candidates,
        replayed_relations: replayed,
        replay_failures: 0,
        duplicate_relations: duplicates,
        unique_relations: unique.len() as u64,
        relation_set_sha256,
        completed: true,
        skip_reason: None,
        record_logical_bytes: record_bytes,
        peak_logical_bytes: records.len() as u64 * record_bytes
            + blocks
                .iter()
                .map(|block| block.columns.len() as u64 * 4)
                .sum::<u64>(),
        disk_bytes_read: 0,
        disk_bytes_written: 0,
        operations,
    })
}

fn quotient_rank(
    dsu: &mut WeightedDsu,
    relations: &BTreeMap<String, SparseRelation>,
    modulus: &BigUint,
) -> QuotientReceipt {
    let pairwise_rank = dsu.independent_edges;
    let (variables, weights, component_count) = dsu.quotient_coordinates();
    let mut rank = SparseRank::default();
    for relation in relations.values() {
        let mut row = BTreeMap::<usize, BigUint>::new();
        for (column, signed_coefficient) in &relation.coefficients {
            let column = *column as usize;
            let magnitude = BigUint::from(signed_coefficient.unsigned_abs());
            let mut coefficient = (&magnitude * &weights[column]) % modulus;
            if *signed_coefficient < 0 && !coefficient.is_zero() {
                coefficient = modulus - coefficient;
            }
            let variable = variables[column];
            let updated = (row.get(&variable).cloned().unwrap_or_default() + coefficient) % modulus;
            if updated.is_zero() {
                row.remove(&variable);
            } else {
                row.insert(variable, updated);
            }
        }
        rank.insert(row, modulus);
    }
    let block_rank = rank.rank();
    QuotientReceipt {
        original_columns: COLUMNS as u64,
        pairwise_relation_rank: pairwise_rank,
        pairwise_components: component_count as u64,
        unique_block_relations: relations.len() as u64,
        independent_block_rank_in_pairwise_quotient: block_rank as u64,
        total_independent_relation_rank: pairwise_rank + block_rank as u64,
        quotient_dimension: (component_count - block_rank) as u64,
        relation_set_sha256: hash_strings(relations.keys().cloned()),
    }
}

fn projection(quotient_dimension: u64) -> ParityProjection {
    let curve = CurveParams::p256();
    let n = curve.n.to_f64().expect("P-256 order fits f64");
    let mut disjoint = 1.0f64;
    for index in 0..9u64 {
        disjoint *= (COLUMNS as u64 - 8 - index) as f64 / (COLUMNS as u64 - index) as f64;
    }
    let rho_s = 1.3;
    let q = quotient_dimension as f64;
    let two_list_s = 2.0 * (q / disjoint).sqrt();
    let memoryless_s = (std::f64::consts::PI / 2.0 * q / disjoint).sqrt();
    let additions_log2 = |s: f64| (s * n.sqrt()).log2();
    ParityProjection {
        group_order: curve.n.to_string(),
        quotient_dimension,
        disjoint_probability: disjoint,
        rho_reference_s: rho_s,
        comparison_unit: "P-256 group-addition equivalents; S=total/sqrt(n)".into(),
        rows: vec![
            ProjectionRow {
                variant: "Pollard rho reference".into(),
                class: "reference".into(),
                group_addition_equivalents_log2: additions_log2(rho_s),
                s_total_over_sqrt_n: rho_s,
                ratio_to_rho: 1.0,
                complete_end_to_end: true,
                correctness: "registered reference".into(),
            },
            ProjectionRow {
                variant: "round-20 direct finite-domain S17 projection".into(),
                class: "projected engineering".into(),
                group_addition_equivalents_log2: 151.885 - 17.0f64.log2(),
                s_total_over_sqrt_n: 2f64.powf(19.419) * rho_s,
                ratio_to_rho: 2f64.powf(19.419),
                complete_end_to_end: true,
                correctness: "exact sampled cells; projected full collection".into(),
            },
            ProjectionRow {
                variant: "round-20 optimistic finite-domain boundary".into(),
                class: "generic lower boundary".into(),
                group_addition_equivalents_log2: 148.668 - 17.0f64.log2(),
                s_total_over_sqrt_n: 2f64.powf(16.203) * rho_s,
                ratio_to_rho: 2f64.powf(16.203),
                complete_end_to_end: false,
                correctness: "analytic lower boundary".into(),
            },
            ProjectionRow {
                variant: "quotient-adjusted balanced two-list boundary".into(),
                class: "generic lower boundary".into(),
                group_addition_equivalents_log2: additions_log2(two_list_s),
                s_total_over_sqrt_n: two_list_s,
                ratio_to_rho: two_list_s / rho_s,
                complete_end_to_end: false,
                correctness: "analytic lower boundary from exact quotient rank".into(),
            },
            ProjectionRow {
                variant: "ideal memoryless quotient boundary".into(),
                class: "unimplemented lower boundary".into(),
                group_addition_equivalents_log2: additions_log2(memoryless_s),
                s_total_over_sqrt_n: memoryless_s,
                ratio_to_rho: memoryless_s / rho_s,
                complete_end_to_end: false,
                correctness: "analytic birthday constant; collector not implemented".into(),
            },
        ],
        two_list_parity_possible: two_list_s <= rho_s,
        ideal_memoryless_parity_possible: memoryless_s <= rho_s,
        implemented_end_to_end_parity: false,
    }
}

fn run(cli: Cli) -> Result<(), String> {
    let dependencies = verify_dependencies(&cli.round20, &cli.round6)?;
    let residual_degree_evidence = load_degree_evidence(&cli.round6)?;
    let curve = CurveParams::p256();
    let inverse_exponent = &curve.p - 2u8;
    let (factor_base, points) = load_factor_base()?;

    let orbit_control = orbit_control(&inverse_exponent, &curve.n)?;
    if !orbit_control.passed {
        return Err("orbit-positive control failed".into());
    }
    let prefix_reference = prefix_reference(&points, &inverse_exponent, &curve.n)?;
    if !prefix_reference.passed {
        return Err("exhaustive prefix reference disagreed".into());
    }

    let (multiplier_census, mut dsu, _) =
        multiplier_census(&points, MAX_MULTIPLIER, &inverse_exponent, &curve.n)?;
    if multiplier_census.replay_failures != 0 || multiplier_census.inconsistent_cycles != 0 {
        return Err("candidate small-multiple exactness failed".into());
    }

    let mut relations = BTreeMap::new();
    let mut block_censuses = Vec::new();
    for depth in BLOCK_DEPTHS {
        block_censuses.push(block_census(
            depth,
            &points,
            &inverse_exponent,
            &mut relations,
        )?);
    }
    let quotient = quotient_rank(&mut dsu, &relations, &curve.n);
    let parity_projection = projection(quotient.quotient_dimension);
    let false_positives =
        prefix_reference.multiplier_false_positives + prefix_reference.block_false_positives;
    let false_negatives =
        prefix_reference.multiplier_false_negatives + prefix_reference.block_false_negatives;
    let exact_group_replay = multiplier_census.replay_failures == 0
        && block_censuses
            .iter()
            .all(|census| census.replay_failures == 0);
    let structured_degree_gate = residual_degree_evidence.structured_residual_maximum <= 5;
    let quotient_dimension_gate = quotient.quotient_dimension == 1;
    let rho_parity_gate = parity_projection.implemented_end_to_end_parity;
    let attack_promotion_gate = orbit_control.passed
        && prefix_reference.passed
        && false_positives == 0
        && false_negatives == 0
        && exact_group_replay
        && structured_degree_gate
        && quotient_dimension_gate
        && rho_parity_gate;
    let decision = if quotient_dimension_gate {
        "transport rank reached one; implement and measure the preregistered memoryless collector before any parity claim"
    } else {
        "negative: registered exact transport families do not collapse FB1 to one logarithm class; rho parity remains blocked by quotient rank"
    };
    let result = ExperimentResult {
        schema: "p256.fb_log_transport/v1".into(),
        curve: CURVE_SLUG.into(),
        dependencies,
        factor_base,
        orbit_control,
        prefix_reference,
        multiplier_census,
        block_censuses,
        quotient,
        residual_degree_evidence,
        parity_projection,
        false_positives,
        false_negatives,
        exact_group_replay,
        structured_degree_gate,
        quotient_dimension_gate,
        rho_parity_gate,
        attack_promotion_gate,
        full_depth_unplanted_attempted: false,
        result_class: "negative/no-rank-compression".into(),
        decision: decision.into(),
    };
    let output = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match cli.out {
        Some(path) => fs::write(path, output).map_err(|error| error.to_string())?,
        None => print!("{output}"),
    }
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
    fn weighted_dsu_preserves_exact_ratios() {
        let modulus = BigUint::from(101u8);
        let mut dsu = WeightedDsu::new(4, &modulus);
        assert!(dsu.unite(0, 1, &BigUint::from(2u8)));
        assert!(dsu.unite(1, 2, &BigUint::from(3u8)));
        assert!(dsu.unite(0, 2, &BigUint::from(6u8)));
        assert_eq!(dsu.components, 2);
        assert_eq!(dsu.independent_edges, 2);
        assert_eq!(dsu.inconsistent_cycles, 0);
    }

    #[test]
    fn canonical_relations_remove_gcd_and_global_sign() {
        let positive = canonical_relation(BTreeMap::from([(3, 2), (7, -4)])).unwrap();
        let negative = canonical_relation(BTreeMap::from([(3, -1), (7, 2)])).unwrap();
        assert_eq!(positive.coefficients, negative.coefficients);
        assert_eq!(relation_digest(&positive), relation_digest(&negative));
    }

    #[test]
    fn parity_requires_a_one_dimensional_quotient_and_new_collector() {
        let one = projection(1);
        let two = projection(2);
        assert!(!one.two_list_parity_possible);
        assert!(one.ideal_memoryless_parity_possible);
        assert!(!one.implemented_end_to_end_parity);
        assert!(!two.ideal_memoryless_parity_possible);
    }

    #[test]
    fn reference_selection_is_stable_and_unique() {
        let (first, first_hash) = select_reference_columns();
        let (second, second_hash) = select_reference_columns();
        assert_eq!(first, second);
        assert_eq!(first_hash, second_hash);
        assert_eq!(first.len(), REFERENCE_POINTS);
        assert_eq!(
            first.iter().copied().collect::<BTreeSet<_>>().len(),
            REFERENCE_POINTS
        );
    }
}
