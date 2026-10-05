//! Batched two-sided birthday census for 17-term P-256 relations.

use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::fs;
use std::mem::size_of;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

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
const KEY_BYTES: usize = 33;
const WIDTHS: [u32; 5] = [20, 22, 24, 26, 28];
const PUBLIC_TARGETS: usize = 256;
const PLANTED_TARGET: u16 = u16::MAX;
const NORMALIZE_BATCH: usize = 8_192;
const ROUND19_SHA256: &str = "3096540621408e4a48cfa18963ad01da9686d3527ee26776c8b6cf6f45a71114";
const ROUND6_SHA256: &str = "71d63031111ba48ff831e79430bc87e6c7bac34626d7d0de4acfb4eaa65d40f4";
const DOMAIN: &str = concat!(
    "icv1-fp256-t89188191154553853111372247798585809583-f188c491/",
    "s17-batch-collision-round20"
);

#[derive(Parser)]
#[command(about = "Measure batched two-sided P-256 S17 collision selection")]
struct Cli {
    /// Frozen round-19 selector result.
    #[arg(long)]
    round19: PathBuf,
    /// Frozen round-6 structured residual-degree result.
    #[arg(long)]
    round6: PathBuf,
    /// Deterministic JSON result path; stdout when omitted.
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Clone, Copy)]
struct FactorPoint {
    positive: Projective,
}

#[derive(Clone)]
struct Target {
    scalar: BigUint,
    point: Projective,
}

#[derive(Clone, Copy)]
struct Representation {
    point: Projective,
    columns: [u32; 9],
    negative: u16,
    len: u8,
}

#[derive(Clone, Copy)]
struct RawRecord {
    point: Projective,
    columns: [u32; 9],
    negative: u16,
    len: u8,
    target: u16,
}

#[derive(Clone, Debug)]
struct SampleRecord {
    projected: u32,
    key: [u8; KEY_BYTES],
    columns: [u32; 9],
    negative: u16,
    len: u8,
    target: u16,
}

impl SampleRecord {
    fn representation_cmp(&self, other: &Self) -> Ordering {
        self.columns
            .cmp(&other.columns)
            .then_with(|| self.negative.cmp(&other.negative))
            .then_with(|| self.target.cmp(&other.target))
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize)]
struct RelationWitness {
    target_kind: String,
    target_index: Option<u16>,
    columns: Vec<u32>,
    negative_columns: Vec<u32>,
}

#[derive(Default, Serialize)]
struct OperationCounts {
    left_sum_group_additions: u64,
    right_sum_group_additions: u64,
    planted_target_group_additions: u64,
    target_generation_group_additions: u64,
    target_generation_group_doublings: u64,
    normalization_field_multiplications: u64,
    replay_group_additions: u64,
    replay_equality_field_multiplications: u64,
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
    dickson_survival_ratio: f64,
}

#[derive(Serialize)]
struct DependencyReceipt {
    round19_sha256: String,
    round6_sha256: String,
}

#[derive(Serialize)]
struct ResidualDegreeEvidence {
    maximum_by_residual_depth: Vec<(u32, u32)>,
    every_component_complete_and_correct: bool,
    local_join_degree: u32,
    structured_residual_maximum: u32,
    unsplit_s17_degree_of_regularity: Option<u32>,
}

#[derive(Serialize)]
struct CellResult {
    projected_key_bits: u32,
    samples_per_side: u64,
    sampled_leaf_visits: u64,
    projected_matches: u64,
    forced_planted_matches: u64,
    expected_random_matches: f64,
    observed_over_expected_random: f64,
    full_key_matches: u64,
    overlapping_full_key_matches: u64,
    disjoint_full_key_matches: u64,
    truncated_false_candidates: u64,
    exact_replay_failures: u64,
    planted_relations: u64,
    unplanted_relations: u64,
    duplicate_relations: u64,
    false_positives: u64,
    false_negatives: u64,
    relation_sha256: String,
    full_key_reference_equal: bool,
    direct_pair_reference_complete: bool,
    direct_pair_comparisons: u64,
    direct_pair_reference_equal: Option<bool>,
    record_logical_bytes: u64,
    logical_resident_bytes: u64,
    disk_bytes_read: u64,
    disk_bytes_written: u64,
    peak_rss_bytes: Option<u64>,
    operations: OperationCounts,
    field_multiplication_equivalents: u64,
    wall_milliseconds_unranked: u64,
    exact_for_sampled_lists: bool,
}

#[derive(Serialize)]
struct GrowthFit {
    widths: Vec<u32>,
    samples_per_side: Vec<u64>,
    field_multiplication_equivalents: Vec<u64>,
    logical_resident_bytes: Vec<u64>,
    work_log2_slope_per_key_bit: f64,
    memory_log2_slope_per_key_bit: f64,
    frozen_prediction: f64,
}

#[derive(Serialize)]
struct SparseLinearAlgebraProjection {
    rows: u64,
    columns: u64,
    row_weight: u64,
    nonzeros: u64,
    working_bytes: u64,
    wiedemann_nonzero_additions: u64,
    berlekamp_massey_scalar_operations: u64,
    optimistic_minimum_fme: u64,
}

#[derive(Serialize)]
struct VariantProjection {
    name: String,
    images_log2: f64,
    group_additions_log2: f64,
    collection_fme_log2: f64,
    total_with_sparse_la_fme_log2: f64,
    per_usable_row_fme_log2: f64,
    total_over_rho_log2: f64,
    class: String,
}

#[derive(Serialize)]
struct FullProjection {
    relation_rows: u64,
    factor_base_columns: u64,
    group_order: String,
    disjoint_probability: f64,
    signed_left_domain: String,
    signed_left_domain_log2: f64,
    signed_right_domain: String,
    signed_right_domain_log2: f64,
    withdrawn_symmetric_samples_per_side_log2: f64,
    target_allowance: f64,
    per_target_relation_mean: f64,
    expected_relation_population: f64,
    collision_events_for_distinct_rows: f64,
    duplicate_collision_allowance: f64,
    required_streamed_right_images: f64,
    required_streamed_right_images_log2: f64,
    measured_normalization_fme_per_image: f64,
    exact_record_bytes: u64,
    direct_left_materialized_bytes_log2: f64,
    rho_reference_fme_log2: f64,
    variants: Vec<VariantProjection>,
    structured_degree_gate: bool,
    exactness_gate: bool,
    per_relation_below_2_pow_103: bool,
    collection_below_2_pow_120: bool,
    storage_below_2_pow_50: bool,
    rho_parity_gate: bool,
    complete_projection: bool,
}

#[derive(Serialize)]
struct ExperimentResult {
    schema: String,
    curve: String,
    dependencies: DependencyReceipt,
    factor_base: FactorBaseReceipt,
    target_count: u64,
    target_scalars: Vec<String>,
    target_set_sha256: String,
    target_generation_group_additions: u64,
    target_generation_group_doublings: u64,
    projected_widths: Vec<u32>,
    cells: Vec<CellResult>,
    growth_fit: GrowthFit,
    residual_degree_evidence: ResidualDegreeEvidence,
    sparse_linear_algebra_projection: SparseLinearAlgebraProjection,
    full_projection: FullProjection,
    exact_for_sampled_lists: bool,
    attack_promotion_gate: bool,
    full_depth_unplanted_attempted: bool,
    decision: String,
}

struct HashStream {
    seed: Vec<u8>,
    block: [u8; 32],
    offset: usize,
    counter: u64,
}

impl HashStream {
    fn new(seed: String) -> Self {
        Self {
            seed: seed.into_bytes(),
            block: [0; 32],
            offset: 32,
            counter: 0,
        }
    }

    fn refill(&mut self) {
        let mut input = self.seed.clone();
        input.extend_from_slice(&self.counter.to_be_bytes());
        self.block = sha256(&input);
        self.counter += 1;
        self.offset = 0;
    }

    fn byte(&mut self) -> u8 {
        if self.offset == self.block.len() {
            self.refill();
        }
        let value = self.block[self.offset];
        self.offset += 1;
        value
    }

    fn u32(&mut self) -> u32 {
        u32::from_be_bytes([self.byte(), self.byte(), self.byte(), self.byte()])
    }

    fn index(&mut self) -> u32 {
        let modulus = COLUMNS as u64;
        let limit = (1u64 << 32) / modulus * modulus;
        loop {
            let value = u64::from(self.u32());
            if value < limit {
                return (value % modulus) as u32;
            }
        }
    }
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("not a non-negative integer: `{value}`"))
}

fn lower_hex(value: &BigUint) -> String {
    format!("0x{}", value.to_str_radix(16))
}

fn file_sha256(path: &Path) -> Result<(Vec<u8>, String), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    Ok((bytes, digest))
}

fn verify_dependencies(round19: &Path, round6: &Path) -> Result<DependencyReceipt, String> {
    let (round19_bytes, round19_sha256) = file_sha256(round19)?;
    if round19_sha256 != ROUND19_SHA256 {
        return Err(format!(
            "round-19 SHA-256 mismatch: expected {ROUND19_SHA256}, got {round19_sha256}"
        ));
    }
    let round19_json: serde_json::Value =
        serde_json::from_slice(&round19_bytes).map_err(|error| error.to_string())?;
    if round19_json["schema"].as_str() != Some("p256.s17_multilevel_selector/v1")
        || round19_json["curve"].as_str() != Some(CURVE_SLUG)
        || round19_json["factor_base"]["fb_id"].as_str() != Some(FB_ID)
        || round19_json["attack_promotion_gate"].as_bool() != Some(false)
    {
        return Err("round-19 identity or rejection boundary changed".into());
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
        round19_sha256,
        round6_sha256,
    })
}

fn load_degree_evidence(path: &Path) -> Result<ResidualDegreeEvidence, String> {
    let (bytes, _) = file_sha256(path)?;
    let dependency: serde_json::Value =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    let cells = dependency["cells"]
        .as_array()
        .ok_or("round-6 cells are missing")?;
    let mut maxima = vec![(1u32, 0u32), (2, 0), (3, 0)];
    let mut complete = true;
    for cell in cells {
        let residual = cell["residual_depth"]
            .as_u64()
            .ok_or("round-6 residual depth is missing")? as u32;
        let degree = cell["max_solving_degree"]
            .as_u64()
            .ok_or("round-6 solving degree is missing")? as u32;
        if let Some((_, maximum)) = maxima.iter_mut().find(|(depth, _)| *depth == residual) {
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
        let x = parse_big(&row.x)?;
        let y = parse_big(&row.y)?;
        let affine = Point::Affine {
            x: curve.fe(x.clone()),
            y: curve.fe(y),
        };
        if !curve.is_on_curve(&affine) {
            return Err(format!("factor-base column {column} is off curve"));
        }
        let mut chain = Fe::from_biguint(&x);
        for _ in 0..18 {
            chain = chain.sqr().sub(&two);
        }
        if chain != terminal {
            return Err(format!(
                "factor-base column {column} misses Dickson terminal"
            ));
        }
        points.push(FactorPoint {
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
            dickson_survival_ratio: 1.0,
        },
        points,
    ))
}

fn scalar_mul_counted(
    base: &Projective,
    scalar: &BigUint,
    operations: &mut OperationCounts,
) -> Projective {
    let mut result = Projective::IDENTITY;
    let mut started = false;
    for bit in (0..scalar.bits()).rev() {
        if started {
            result = result.add(&result);
            operations.target_generation_group_doublings += 1;
        }
        if scalar.bit(bit) {
            result = result.add(base);
            operations.target_generation_group_additions += 1;
            started = true;
        }
    }
    result
}

fn make_targets(operations: &mut OperationCounts) -> Result<(Vec<Target>, String), String> {
    let curve = CurveParams::p256();
    let generator = Projective::from_textbook(&curve.generator());
    let mut targets = Vec::with_capacity(PUBLIC_TARGETS);
    let mut receipt = Vec::with_capacity(PUBLIC_TARGETS * 65);
    for index in 0..PUBLIC_TARGETS {
        let mut scalar = BigUint::from_bytes_be(&sha256(
            format!("{DOMAIN}/unplanted-target/{index}").as_bytes(),
        )) % &curve.n;
        if scalar.is_zero() {
            scalar = BigUint::one();
        }
        let point = scalar_mul_counted(&generator, &scalar, operations);
        if bool::from(point.is_identity()) {
            return Err(format!("unplanted target {index} is infinity"));
        }
        receipt.extend_from_slice(&scalar.to_bytes_be());
        receipt.extend_from_slice(&point_key(&point)?);
        targets.push(Target { scalar, point });
    }
    Ok((targets, hex::encode(sha256(&receipt))))
}

fn sample_representation(
    points: &[FactorPoint],
    width: u32,
    side: &str,
    ordinal: u64,
    len: usize,
    excluded: &[u32],
) -> Result<Representation, String> {
    let mut stream = HashStream::new(format!("{DOMAIN}/sample/{width}/{side}/{ordinal}"));
    let mut selected = Vec::with_capacity(len);
    while selected.len() < len {
        let column = stream.index();
        if !selected.contains(&column) && !excluded.contains(&column) {
            selected.push(column);
        }
    }
    let mut signed: Vec<(u32, bool)> = selected
        .into_iter()
        .map(|column| (column, stream.byte() & 1 == 1))
        .collect();
    signed.sort_unstable_by_key(|(column, _)| *column);
    let mut columns = [0u32; 9];
    let mut negative = 0u16;
    let mut point = Projective::IDENTITY;
    for (position, (column, is_negative)) in signed.iter().copied().enumerate() {
        columns[position] = column;
        if is_negative {
            negative |= 1u16 << position;
        }
        let selected_point = if is_negative {
            points[column as usize].positive.neg()
        } else {
            points[column as usize].positive
        };
        if position == 0 {
            point = selected_point;
        } else {
            point = point.add(&selected_point);
        }
    }
    if bool::from(point.is_identity()) {
        return Err(format!("sample {width}/{side}/{ordinal} is infinity"));
    }
    Ok(Representation {
        point,
        columns,
        negative,
        len: len as u8,
    })
}

fn target_index(width: u32, ordinal: u64) -> u16 {
    let digest = sha256(format!("{DOMAIN}/sample/{width}/right/{ordinal}/target").as_bytes());
    u16::from_be_bytes([digest[30], digest[31]]) % PUBLIC_TARGETS as u16
}

fn projected_key(key: &[u8; KEY_BYTES], width: u32) -> u32 {
    let mut input = format!("{DOMAIN}/state-key/v1").into_bytes();
    input.extend_from_slice(key);
    let digest = sha256(&input);
    let word = u32::from_be_bytes(digest[28..32].try_into().expect("four bytes"));
    word & ((1u32 << width) - 1)
}

fn point_key(point: &Projective) -> Result<[u8; KEY_BYTES], String> {
    let (x, y) = point.to_affine().ok_or("point is infinity")?;
    let mut key = [0u8; KEY_BYTES];
    key[0] = 2 | (y.to_bytes_be()[31] & 1);
    key[1..].copy_from_slice(&x.to_bytes_be());
    Ok(key)
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

fn normalize_records(
    raw: &[RawRecord],
    width: u32,
    inverse_exponent: &BigUint,
    multiplications: &mut u64,
) -> Result<Vec<SampleRecord>, String> {
    let mut records = Vec::with_capacity(raw.len());
    for batch in raw.chunks(NORMALIZE_BATCH) {
        if batch
            .iter()
            .any(|record| bool::from(record.point.is_identity()))
        {
            return Err("sample record became infinity".into());
        }
        let mut prefixes = Vec::with_capacity(batch.len());
        let mut product = Fe::ONE;
        for record in batch {
            prefixes.push(product);
            product = product.mul(&record.point.z);
            *multiplications += 1;
        }
        let mut inverse = pow_counted(product, inverse_exponent, multiplications);
        let mut inverses = vec![Fe::ZERO; batch.len()];
        for index in (0..batch.len()).rev() {
            inverses[index] = inverse.mul(&prefixes[index]);
            *multiplications += 1;
            if index != 0 {
                inverse = inverse.mul(&batch[index].point.z);
                *multiplications += 1;
            }
        }
        for (record, inverse) in batch.iter().zip(inverses) {
            let x = record.point.x.mul(&inverse);
            let y = record.point.y.mul(&inverse);
            *multiplications += 2;
            let mut key = [0u8; KEY_BYTES];
            key[0] = 2 | (y.to_bytes_be()[31] & 1);
            key[1..].copy_from_slice(&x.to_bytes_be());
            records.push(SampleRecord {
                projected: projected_key(&key, width),
                key,
                columns: record.columns,
                negative: record.negative,
                len: record.len,
                target: record.target,
            });
        }
    }
    Ok(records)
}

fn record_sort(left: &SampleRecord, right: &SampleRecord) -> Ordering {
    left.projected
        .cmp(&right.projected)
        .then_with(|| left.key.cmp(&right.key))
        .then_with(|| left.representation_cmp(right))
}

fn pair_witness(left: &SampleRecord, right: &SampleRecord) -> Option<RelationWitness> {
    let mut signed = Vec::with_capacity(17);
    for position in 0..left.len as usize {
        signed.push((
            left.columns[position],
            left.negative & (1u16 << position) != 0,
        ));
    }
    for position in 0..right.len as usize {
        signed.push((
            right.columns[position],
            right.negative & (1u16 << position) != 0,
        ));
    }
    signed.sort_unstable_by_key(|(column, _)| *column);
    if signed.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return None;
    }
    Some(RelationWitness {
        target_kind: if right.target == PLANTED_TARGET {
            "planted-control".into()
        } else {
            "hash-unplanted".into()
        },
        target_index: (right.target != PLANTED_TARGET).then_some(right.target),
        columns: signed.iter().map(|(column, _)| *column).collect(),
        negative_columns: signed
            .iter()
            .filter_map(|(column, negative)| negative.then_some(*column))
            .collect(),
    })
}

fn replay_witness(
    witness: &RelationWitness,
    points: &[FactorPoint],
    target: &Projective,
    operations: &mut OperationCounts,
) -> bool {
    let negatives: BTreeSet<u32> = witness.negative_columns.iter().copied().collect();
    let mut sum = Projective::IDENTITY;
    for (position, column) in witness.columns.iter().copied().enumerate() {
        let point = if negatives.contains(&column) {
            points[column as usize].positive.neg()
        } else {
            points[column as usize].positive
        };
        if position == 0 {
            sum = point;
        } else {
            sum = sum.add(&point);
            operations.replay_group_additions += 1;
        }
    }
    if bool::from(sum.is_identity()) || bool::from(target.is_identity()) {
        return bool::from(sum.is_identity()) == bool::from(target.is_identity());
    }
    operations.replay_equality_field_multiplications += 4;
    sum.x.mul(&target.z) == target.x.mul(&sum.z) && sum.y.mul(&target.z) == target.y.mul(&sum.z)
}

fn relation_digest(relations: &BTreeSet<RelationWitness>) -> String {
    let bytes = serde_json::to_vec(relations).expect("relation serialization");
    hex::encode(sha256(&bytes))
}

fn target_for<'a>(
    right: &SampleRecord,
    planted: &'a Projective,
    targets: &'a [Target],
) -> &'a Projective {
    if right.target == PLANTED_TARGET {
        planted
    } else {
        &targets[right.target as usize].point
    }
}

fn exact_relations_from_pair(
    left: &SampleRecord,
    right: &SampleRecord,
    points: &[FactorPoint],
    planted: &Projective,
    targets: &[Target],
    operations: &mut OperationCounts,
    replay_failures: &mut u64,
) -> Option<RelationWitness> {
    let witness = pair_witness(left, right)?;
    if replay_witness(
        &witness,
        points,
        target_for(right, planted, targets),
        operations,
    ) {
        Some(witness)
    } else {
        *replay_failures += 1;
        None
    }
}

fn full_key_reference(
    left: &[SampleRecord],
    right: &[SampleRecord],
    points: &[FactorPoint],
    planted: &Projective,
    targets: &[Target],
) -> Result<BTreeSet<RelationWitness>, String> {
    let mut left_order: Vec<usize> = (0..left.len()).collect();
    let mut right_order: Vec<usize> = (0..right.len()).collect();
    left_order.sort_unstable_by(|a, b| left[*a].key.cmp(&left[*b].key));
    right_order.sort_unstable_by(|a, b| right[*a].key.cmp(&right[*b].key));
    let mut relations = BTreeSet::new();
    let mut operations = OperationCounts::default();
    let mut failures = 0u64;
    let (mut i, mut j) = (0usize, 0usize);
    while i < left_order.len() && j < right_order.len() {
        let left_key = left[left_order[i]].key;
        let right_key = right[right_order[j]].key;
        match left_key.cmp(&right_key) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                let i_end =
                    i + left_order[i..].partition_point(|index| left[*index].key == left_key);
                let j_end =
                    j + right_order[j..].partition_point(|index| right[*index].key == right_key);
                for &left_index in &left_order[i..i_end] {
                    for &right_index in &right_order[j..j_end] {
                        if let Some(witness) = exact_relations_from_pair(
                            &left[left_index],
                            &right[right_index],
                            points,
                            planted,
                            targets,
                            &mut operations,
                            &mut failures,
                        ) {
                            relations.insert(witness);
                        }
                    }
                }
                i = i_end;
                j = j_end;
            }
        }
    }
    if failures != 0 {
        return Err(format!("full-key reference had {failures} replay failures"));
    }
    Ok(relations)
}

fn direct_pair_reference(
    left: &[SampleRecord],
    right: &[SampleRecord],
) -> (BTreeSet<RelationWitness>, u64) {
    let mut relations = BTreeSet::new();
    let mut comparisons = 0u64;
    for left_record in left {
        for right_record in right {
            comparisons += 1;
            if left_record.key == right_record.key {
                if let Some(witness) = pair_witness(left_record, right_record) {
                    relations.insert(witness);
                }
            }
        }
    }
    (relations, comparisons)
}

fn peak_rss_bytes() -> Option<u64> {
    let status = fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|line| line.starts_with("VmHWM:"))?;
    let kib = line.split_whitespace().nth(1)?.parse::<u64>().ok()?;
    Some(kib * 1024)
}

fn run_cell(width: u32, points: &[FactorPoint], targets: &[Target]) -> Result<CellResult, String> {
    let started = Instant::now();
    let samples = 1usize << (width / 2 + 3);
    let left_plant = sample_representation(points, width, "left", 0, 8, &[])?;
    let excluded = &left_plant.columns[..left_plant.len as usize];
    let right_plant = sample_representation(points, width, "right", 0, 9, excluded)?;
    let planted = left_plant.point.add(&right_plant.point);
    if bool::from(planted.is_identity()) {
        return Err(format!("planted target is infinity at width {width}"));
    }
    let mut operations = OperationCounts {
        planted_target_group_additions: 1,
        ..OperationCounts::default()
    };

    let mut raw_left = Vec::with_capacity(samples);
    let mut raw_right = Vec::with_capacity(samples);
    raw_left.push(RawRecord {
        point: left_plant.point,
        columns: left_plant.columns,
        negative: left_plant.negative,
        len: left_plant.len,
        target: PLANTED_TARGET,
    });
    raw_right.push(RawRecord {
        point: planted.add(&right_plant.point.neg()),
        columns: right_plant.columns,
        negative: right_plant.negative,
        len: right_plant.len,
        target: PLANTED_TARGET,
    });
    operations.left_sum_group_additions += 7;
    operations.right_sum_group_additions += 9;

    for ordinal in 1..samples as u64 {
        let left = sample_representation(points, width, "left", ordinal, 8, &[])?;
        operations.left_sum_group_additions += 7;
        raw_left.push(RawRecord {
            point: left.point,
            columns: left.columns,
            negative: left.negative,
            len: left.len,
            target: PLANTED_TARGET,
        });

        let right = sample_representation(points, width, "right", ordinal, 9, &[])?;
        let target_index = target_index(width, ordinal);
        let adjusted = targets[target_index as usize].point.add(&right.point.neg());
        operations.right_sum_group_additions += 9;
        raw_right.push(RawRecord {
            point: adjusted,
            columns: right.columns,
            negative: right.negative,
            len: right.len,
            target: target_index,
        });
    }

    let inverse_exponent = &CurveParams::p256().p - BigUint::from(2u8);
    let mut normalization = 0u64;
    let mut left = normalize_records(&raw_left, width, &inverse_exponent, &mut normalization)?;
    let mut right = normalize_records(&raw_right, width, &inverse_exponent, &mut normalization)?;
    operations.normalization_field_multiplications = normalization;
    drop(raw_left);
    drop(raw_right);

    let reference = full_key_reference(&left, &right, points, &planted, targets)?;
    let direct_reference = if width == WIDTHS[0] {
        Some(direct_pair_reference(&left, &right))
    } else {
        None
    };

    left.sort_unstable_by(record_sort);
    right.sort_unstable_by(record_sort);
    let mut relations = BTreeSet::new();
    let mut projected_matches = 0u64;
    let mut full_key_matches = 0u64;
    let mut overlapping = 0u64;
    let mut disjoint = 0u64;
    let mut replay_failures = 0u64;
    let mut duplicates = 0u64;
    let (mut i, mut j) = (0usize, 0usize);
    while i < left.len() && j < right.len() {
        match left[i].projected.cmp(&right[j].projected) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                let projected = left[i].projected;
                let i_end = i + left[i..].partition_point(|record| record.projected == projected);
                let j_end = j + right[j..].partition_point(|record| record.projected == projected);
                projected_matches += ((i_end - i) * (j_end - j)) as u64;
                for left_record in &left[i..i_end] {
                    for right_record in &right[j..j_end] {
                        if left_record.key != right_record.key {
                            continue;
                        }
                        full_key_matches += 1;
                        let Some(witness) = pair_witness(left_record, right_record) else {
                            overlapping += 1;
                            continue;
                        };
                        disjoint += 1;
                        if replay_witness(
                            &witness,
                            points,
                            target_for(right_record, &planted, targets),
                            &mut operations,
                        ) {
                            if !relations.insert(witness) {
                                duplicates += 1;
                            }
                        } else {
                            replay_failures += 1;
                        }
                    }
                }
                i = i_end;
                j = j_end;
            }
        }
    }

    let full_reference_equal = relations == reference;
    if !full_reference_equal {
        return Err(format!("full-key reference disagreement at width {width}"));
    }
    let direct_equal = direct_reference
        .as_ref()
        .map(|(direct_relations, _)| direct_relations == &relations);
    if direct_equal == Some(false) {
        return Err("direct pair reference disagreement".into());
    }
    let planted_relations = relations
        .iter()
        .filter(|relation| relation.target_kind == "planted-control")
        .count() as u64;
    let unplanted_relations = relations.len() as u64 - planted_relations;
    if planted_relations != 1 {
        return Err(format!(
            "expected one planted relation at width {width}, got {planted_relations}"
        ));
    }
    let expected_random = (samples as f64).powi(2) / 2f64.powi(width as i32);
    let observed_random = projected_matches.saturating_sub(1) as f64;
    let record_bytes = size_of::<SampleRecord>() as u64;
    let logical_bytes = 2 * samples as u64 * record_bytes;
    let fme = (operations.left_sum_group_additions
        + operations.right_sum_group_additions
        + operations.planted_target_group_additions
        + operations.target_generation_group_additions
        + operations.target_generation_group_doublings
        + operations.replay_group_additions)
        .saturating_mul(17)
        .saturating_add(operations.normalization_field_multiplications)
        .saturating_add(operations.replay_equality_field_multiplications);
    let exact = replay_failures == 0 && full_reference_equal && direct_equal != Some(false);
    Ok(CellResult {
        projected_key_bits: width,
        samples_per_side: samples as u64,
        sampled_leaf_visits: 17 * samples as u64,
        projected_matches,
        forced_planted_matches: 1,
        expected_random_matches: expected_random,
        observed_over_expected_random: observed_random / expected_random,
        full_key_matches,
        overlapping_full_key_matches: overlapping,
        disjoint_full_key_matches: disjoint,
        truncated_false_candidates: projected_matches - full_key_matches,
        exact_replay_failures: replay_failures,
        planted_relations,
        unplanted_relations,
        duplicate_relations: duplicates,
        false_positives: replay_failures,
        false_negatives: 0,
        relation_sha256: relation_digest(&relations),
        full_key_reference_equal: full_reference_equal,
        direct_pair_reference_complete: direct_reference.is_some(),
        direct_pair_comparisons: direct_reference.as_ref().map_or(0, |(_, count)| *count),
        direct_pair_reference_equal: direct_equal,
        record_logical_bytes: record_bytes,
        logical_resident_bytes: logical_bytes,
        disk_bytes_read: 0,
        disk_bytes_written: 0,
        peak_rss_bytes: peak_rss_bytes(),
        operations,
        field_multiplication_equivalents: fme,
        wall_milliseconds_unranked: started.elapsed().as_millis() as u64,
        exact_for_sampled_lists: exact,
    })
}

fn slope(xs: &[f64], ys: &[f64]) -> f64 {
    let count = xs.len() as f64;
    let mean_x = xs.iter().sum::<f64>() / count;
    let mean_y = ys.iter().sum::<f64>() / count;
    let numerator = xs
        .iter()
        .zip(ys)
        .map(|(x, y)| (x - mean_x) * (y - mean_y))
        .sum::<f64>();
    let denominator = xs.iter().map(|x| (x - mean_x).powi(2)).sum::<f64>();
    numerator / denominator
}

fn choose_big(n: u64, k: u32) -> BigUint {
    let k = u64::from(k).min(n - u64::from(k));
    let mut value = BigUint::one();
    for index in 0..k {
        value *= BigUint::from(n - index);
        value /= BigUint::from(index + 1);
    }
    value
}

fn projection(
    record_bytes: u64,
    exactness_gate: bool,
    degree_gate: bool,
    normalization_fme_per_image: f64,
    measured_target_additions: u64,
    measured_target_doublings: u64,
) -> (SparseLinearAlgebraProjection, FullProjection) {
    let rows = 138_031u64;
    let row_weight = 17u64;
    let nonzeros = rows * row_weight;
    let csr_bytes = nonzeros * 5 + (rows + 1) * 8;
    let field_vector_bytes = 3 * COLUMNS as u64 * 32;
    let wiedemann = 2 * COLUMNS as u64 * nonzeros;
    let berlekamp = COLUMNS as u64 * COLUMNS as u64;
    let sparse_minimum = wiedemann + berlekamp;
    let sparse = SparseLinearAlgebraProjection {
        rows,
        columns: COLUMNS as u64,
        row_weight,
        nonzeros,
        working_bytes: csr_bytes + field_vector_bytes,
        wiedemann_nonzero_additions: wiedemann,
        berlekamp_massey_scalar_operations: berlekamp,
        optimistic_minimum_fme: sparse_minimum,
    };

    let mut disjoint = 1.0f64;
    for index in 0..9u64 {
        disjoint *= (COLUMNS as u64 - 8 - index) as f64 / (COLUMNS as u64 - index) as f64;
    }
    let curve = CurveParams::p256();
    let order = curve.n.to_f64().expect("P-256 order fits f64");
    let signed_left = choose_big(COLUMNS as u64, 8) << 8usize;
    let signed_right = choose_big(COLUMNS as u64, 9) << 9usize;
    let signed_domain = choose_big(COLUMNS as u64, 17) << 17usize;
    let left_images = signed_left.to_f64().expect("left domain fits f64");
    let right_domain = signed_right.to_f64().expect("right domain fits f64");
    let relation_mean = signed_domain.to_f64().expect("relation domain fits f64") / order;
    let whole_success = 0.964_000_182_074_505_8f64;
    let target_allowance = rows as f64 / whole_success;
    let relation_population = relation_mean * target_allowance;
    let collision_events = -relation_population * (1.0 - rows as f64 / relation_population).ln();
    let duplicate_allowance = collision_events / rows as f64;
    let right_images = collision_events * order / (left_images * disjoint);
    assert!(right_images / target_allowance < right_domain);
    let images = left_images + right_images;
    let withdrawn_symmetric = ((rows as f64 * order) / disjoint).sqrt();
    let rho_fme = 1.3 * order.sqrt() * 17.0;
    let direct_additions = 7.0 * left_images + 9.0 * right_images;
    let optimistic_additions = images;
    let replay_fme = rows as f64 * (16.0 * 17.0 + 4.0);
    let target_ops_per_target =
        (measured_target_additions + measured_target_doublings) as f64 / PUBLIC_TARGETS as f64;
    let target_setup_fme = target_allowance * target_ops_per_target * 17.0;
    let variants = [
        (
            "direct random signed 8+9 sums",
            direct_additions,
            images * normalization_fme_per_image + replay_fme + target_setup_fme,
            "projected algorithmic advance",
        ),
        (
            "optimistic one-addition-per-image generic boundary",
            optimistic_additions,
            0.0,
            "generic lower boundary",
        ),
    ]
    .into_iter()
    .map(|(name, additions, other_collection_fme, class)| {
        let collection_fme = additions * 17.0 + other_collection_fme;
        let total = collection_fme + sparse_minimum as f64;
        VariantProjection {
            name: name.into(),
            images_log2: images.log2(),
            group_additions_log2: additions.log2(),
            collection_fme_log2: collection_fme.log2(),
            total_with_sparse_la_fme_log2: total.log2(),
            per_usable_row_fme_log2: (total / rows as f64).log2(),
            total_over_rho_log2: (total / rho_fme).log2(),
            class: class.into(),
        }
    })
    .collect::<Vec<_>>();
    let direct = &variants[0];
    let storage_log2 = (left_images * record_bytes as f64).log2();
    let per_relation_gate = direct.per_usable_row_fme_log2 < 103.0;
    let collection_gate = direct.total_with_sparse_la_fme_log2 < 120.0;
    let storage_gate = storage_log2 < 50.0;
    let rho_gate = direct.total_over_rho_log2 < 0.0;
    (
        sparse,
        FullProjection {
            relation_rows: rows,
            factor_base_columns: COLUMNS as u64,
            group_order: curve.n.to_string(),
            disjoint_probability: disjoint,
            signed_left_domain: signed_left.to_string(),
            signed_left_domain_log2: left_images.log2(),
            signed_right_domain: signed_right.to_string(),
            signed_right_domain_log2: right_domain.log2(),
            withdrawn_symmetric_samples_per_side_log2: withdrawn_symmetric.log2(),
            target_allowance,
            per_target_relation_mean: relation_mean,
            expected_relation_population: relation_population,
            collision_events_for_distinct_rows: collision_events,
            duplicate_collision_allowance: duplicate_allowance,
            required_streamed_right_images: right_images,
            required_streamed_right_images_log2: right_images.log2(),
            measured_normalization_fme_per_image: normalization_fme_per_image,
            exact_record_bytes: record_bytes,
            direct_left_materialized_bytes_log2: storage_log2,
            rho_reference_fme_log2: rho_fme.log2(),
            variants,
            structured_degree_gate: degree_gate,
            exactness_gate,
            per_relation_below_2_pow_103: per_relation_gate,
            collection_below_2_pow_120: collection_gate,
            storage_below_2_pow_50: storage_gate,
            rho_parity_gate: rho_gate,
            complete_projection: true,
        },
    )
}

fn run(cli: Cli) -> Result<(), String> {
    let dependencies = verify_dependencies(&cli.round19, &cli.round6)?;
    let residual_degree_evidence = load_degree_evidence(&cli.round6)?;
    let (factor_base, points) = load_factor_base()?;
    let mut target_operations = OperationCounts::default();
    let (targets, target_set_sha256) = make_targets(&mut target_operations)?;
    let target_setup = (
        target_operations.target_generation_group_additions,
        target_operations.target_generation_group_doublings,
    );

    let mut cells = Vec::new();
    for width in WIDTHS {
        cells.push(run_cell(width, &points, &targets)?);
    }
    let exact = cells.iter().all(|cell| cell.exact_for_sampled_lists);
    let xs = cells
        .iter()
        .map(|cell| cell.projected_key_bits as f64)
        .collect::<Vec<_>>();
    let work = cells
        .iter()
        .map(|cell| (cell.field_multiplication_equivalents as f64).log2())
        .collect::<Vec<_>>();
    let memory = cells
        .iter()
        .map(|cell| (cell.logical_resident_bytes as f64).log2())
        .collect::<Vec<_>>();
    let growth_fit = GrowthFit {
        widths: cells.iter().map(|cell| cell.projected_key_bits).collect(),
        samples_per_side: cells.iter().map(|cell| cell.samples_per_side).collect(),
        field_multiplication_equivalents: cells
            .iter()
            .map(|cell| cell.field_multiplication_equivalents)
            .collect(),
        logical_resident_bytes: cells
            .iter()
            .map(|cell| cell.logical_resident_bytes)
            .collect(),
        work_log2_slope_per_key_bit: slope(&xs, &work),
        memory_log2_slope_per_key_bit: slope(&xs, &memory),
        frozen_prediction: 0.5,
    };
    let record_bytes = cells[0].record_logical_bytes;
    let degree_gate = residual_degree_evidence.structured_residual_maximum <= 5;
    let largest_cell = cells.last().expect("nonempty width ladder");
    let normalization_fme_per_image = largest_cell.operations.normalization_field_multiplications
        as f64
        / (2 * largest_cell.samples_per_side) as f64;
    let (sparse_linear_algebra_projection, full_projection) = projection(
        record_bytes,
        exact,
        degree_gate,
        normalization_fme_per_image,
        target_setup.0,
        target_setup.1,
    );
    let attack_promotion_gate = exact
        && degree_gate
        && full_projection.per_relation_below_2_pow_103
        && full_projection.collection_below_2_pow_120
        && full_projection.storage_below_2_pow_50
        && full_projection.rho_parity_gate;
    let result = ExperimentResult {
        schema: "p256.s17_batch_collision/v1".into(),
        curve: CURVE_SLUG.into(),
        dependencies,
        factor_base,
        target_count: targets.len() as u64,
        target_scalars: targets
            .iter()
            .map(|target| lower_hex(&target.scalar))
            .collect(),
        target_set_sha256,
        target_generation_group_additions: target_setup.0,
        target_generation_group_doublings: target_setup.1,
        projected_widths: WIDTHS.to_vec(),
        cells,
        growth_fit,
        residual_degree_evidence,
        sparse_linear_algebra_projection,
        full_projection,
        exact_for_sampled_lists: exact,
        attack_promotion_gate,
        full_depth_unplanted_attempted: false,
        decision: if attack_promotion_gate {
            "all gates pass; an exact full-depth unplanted collision run is authorized".into()
        } else {
            "negative: batching removes per-target frontier repetition, but even the generic image floor remains above rho".into()
        },
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
    fn sampler_indices_are_unique_and_in_range() {
        let mut stream = HashStream::new(format!("{DOMAIN}/test/sampler"));
        let mut values = Vec::new();
        while values.len() < 17 {
            let value = stream.index();
            if !values.contains(&value) {
                values.push(value);
            }
        }
        assert_eq!(values.iter().copied().collect::<BTreeSet<_>>().len(), 17);
        assert!(values.iter().all(|value| *value < COLUMNS as u32));
    }

    #[test]
    fn projected_keys_are_nested() {
        let key = [0x5au8; KEY_BYTES];
        let low20 = projected_key(&key, 20);
        let low28 = projected_key(&key, 28);
        assert_eq!(low20, low28 & ((1u32 << 20) - 1));
    }

    #[test]
    fn generic_batch_floor_remains_above_rho() {
        let (_, full) = projection(64, true, true, 5.0, 32_891, 65_064);
        assert!(full.variants[1].total_over_rho_log2 > 15.0);
        assert!(full.required_streamed_right_images_log2 > 144.0);
        assert!(!full.rho_parity_gate);
        assert!(!full.collection_below_2_pow_120);
    }
}
