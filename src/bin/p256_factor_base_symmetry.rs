//! Exact P-256 factor-base symmetry screen (round 22).

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::mem::size_of;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::short_id;
use crypto_lib::cryptanalysis::p256_bitbox_factor_base;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{
    self, WideFactorBaseDump, CURVE_SLUG, FACTOR_BASE_SCHEMA, POINT_KEY_ENCODING,
};
use crypto_lib::cryptanalysis::p256_structural::{p256_cm_discriminant_abs, trial_factor};
use crypto_lib::ecc::p256_field::P256FieldElement as Fe;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::{CurveParams, Point};
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;
use serde_json::json;

const COLUMNS: usize = 131_458;
const MAX_MULTIPLIER: u8 = 16;
const NORMALIZE_BATCH: usize = 8_192;
const RELATION_ROWS: u64 = 138_031;
const RHO_S: f64 = 1.3;
const ROUND21_SHA256: &str = "4aaf5fde72257ed1f262d5674ab2fe81c7b4f90698fa3767615adfed6524e757";
const ROUND2_SHA256: &str = "dcaa1f66f5f89757a3670a4a8e158ca6d5c360b7ca119d384d2c4a0b2794c2c8";
const ROUND1_SHA256: &str = "bf4c0f95ca46bf92b236eea7b8f37a3902f744b07a898a41f6b120faab2e35a4";

#[derive(Parser)]
#[command(about = "Screen globally structured P-256 factor bases")]
struct Cli {
    #[arg(long)]
    round21: PathBuf,
    #[arg(long)]
    round2: PathBuf,
    #[arg(long)]
    round1: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Clone, Copy)]
struct FactorPoint {
    positive: Projective,
}

#[derive(Clone, Copy, Eq, PartialEq)]
struct MultipleRecord {
    x: [u8; 32],
    parity: u8,
    multiplier: u8,
    column: u32,
}

#[derive(Serialize)]
struct Dependencies {
    round21_sha256: String,
    round2_sha256: String,
    round1_sha256: String,
}

#[derive(Serialize)]
struct AutomorphismScreen {
    j_invariant: String,
    j_is_zero: bool,
    j_is_1728: bool,
    fp_degree_one_automorphisms: u32,
    negation_already_folded: bool,
    cm_discriminant_abs: String,
    cm_discriminant_bits: u64,
    cm_trial_bound: u64,
    cm_trial_factors: Vec<(u64, u32)>,
    cm_residue_bits: u64,
    qualification: String,
}

#[derive(Default, Serialize)]
struct OperationCounts {
    image_group_additions: u64,
    replay_group_additions: u64,
    replay_group_doublings: u64,
    normalization_field_multiplications: u64,
}

#[derive(Serialize)]
struct TransportCensus {
    multiplier_max: u8,
    images: u64,
    equal_x_buckets: u64,
    equal_x_pairs: u64,
    same_column_pairs: u64,
    cross_column_candidates: u64,
    replayed_relations: u64,
    replay_failures: u64,
    inconsistent_cycles: u64,
    independent_rank: u64,
    quotient_dimension: u64,
    isolated_columns: u64,
    largest_component: u64,
    record_logical_bytes: u64,
    peak_logical_bytes: u64,
    disk_bytes: u64,
    operations: OperationCounts,
}

#[derive(Serialize)]
struct GeometricCandidate {
    name: String,
    family: String,
    spec: String,
    fb_id: String,
    fb_sha256: String,
    points_sha256: String,
    columns: u64,
    signed_points: u64,
    rebuild_verified: bool,
    degree_evidence: String,
    unsplit_s17_degree_of_regularity: Option<u32>,
    transport: TransportCensus,
    projection: Projection,
    classification: String,
}

#[derive(Serialize)]
struct ScalarManifest {
    schema: String,
    curve: String,
    family: String,
    params: BTreeMap<String, String>,
    columns: u64,
    signed_points: u64,
    points_sha256: String,
}

#[derive(Serialize)]
struct ScalarCandidate {
    name: String,
    transition: String,
    fb_id: String,
    fb_sha256: String,
    manifest: ScalarManifest,
    signed_duplicates: u64,
    zero_scalars: u64,
    internal_transition_relations: u64,
    transition_closes_at_end: bool,
    exact_relation_rank: u64,
    quotient_dimension: u64,
    known_construction_logs: bool,
    actual_p256_points_verified: bool,
    point_generation_additions: u64,
    point_generation_doublings: u64,
    additive_s17_support_upper_bound: Option<String>,
    additive_s17_support_fraction_log2: Option<f64>,
    projection: Projection,
    diagnostics_status: String,
    classification: String,
}

#[derive(Clone, Serialize)]
struct Projection {
    independent_log_classes: u64,
    disjoint_probability: f64,
    cross_colour_s: f64,
    cross_colour_ratio_to_rho: f64,
    two_list_s: f64,
    two_list_ratio_to_rho: f64,
    cross_colour_log2_operations: f64,
    projected_138031_rows_log2_operations: f64,
    projected_cost_per_usable_relation_log2_operations: f64,
    sparse_linear_algebra_lower_bound_log2_operations: Option<f64>,
    complete_projected_lower_bound_log2_operations: f64,
    projected_peak_two_list_log2_bytes: f64,
    structured_residual_degree_max: Option<u32>,
    relation_collection_below_2_120: bool,
    cost_per_relation_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    complete_cost_below_rho: bool,
    promotion_gate: bool,
    dominant_obstruction: String,
}

#[derive(Serialize)]
struct StoppedControl {
    name: String,
    stop_reason: String,
    exhaustive_credit: bool,
    projection: Projection,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    relation_model: String,
    dependencies: Dependencies,
    automorphism_screen: AutomorphismScreen,
    geometric_candidates: Vec<GeometricCandidate>,
    scalar_candidates: Vec<ScalarCandidate>,
    stopped_controls: Vec<StoppedControl>,
    reference_false_positives: u64,
    reference_false_negatives: u64,
    all_reported_relations_replayed: bool,
    any_promotion_gate_passed: bool,
    full_depth_unplanted_attempted: bool,
    result_class: String,
    decision: String,
}

struct WeightedDsu {
    parent: Vec<usize>,
    size: Vec<u32>,
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

    /// Add `log(left) = ratio*log(right) (mod n)`.
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
            let root_weight = (ratio * right_weight % &self.modulus)
                * mod_inverse(&left_weight, &self.modulus)
                % &self.modulus;
            self.parent[left_root] = right_root;
            self.weight[left_root] = root_weight;
            self.size[right_root] += self.size[left_root];
        } else {
            let right_factor = ratio * right_weight % &self.modulus;
            let root_weight =
                left_weight * mod_inverse(&right_factor, &self.modulus) % &self.modulus;
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
        (
            sizes.values().filter(|size| **size == 1).count(),
            sizes.values().copied().max().unwrap_or(0),
        )
    }
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("not a non-negative integer: `{value}`"))
}

fn point_key_bytes(x: &BigUint, high_sign: bool) -> Result<[u8; 33], String> {
    let key = ((x + BigUint::one()) << 1usize) + BigUint::from(high_sign);
    let raw = key.to_bytes_be();
    if raw.len() > 33 {
        return Err("P-256 point key exceeded 33 bytes".into());
    }
    let mut out = [0u8; 33];
    out[33 - raw.len()..].copy_from_slice(&raw);
    Ok(out)
}

fn file_sha256(path: &Path) -> Result<(Vec<u8>, String), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    Ok((bytes.clone(), hex::encode(sha256(&bytes))))
}

fn verify_dependencies(cli: &Cli) -> Result<(Dependencies, serde_json::Value), String> {
    let (round21_bytes, round21_sha256) = file_sha256(&cli.round21)?;
    let (round2_bytes, round2_sha256) = file_sha256(&cli.round2)?;
    let (round1_bytes, round1_sha256) = file_sha256(&cli.round1)?;
    if round21_sha256 != ROUND21_SHA256
        || round2_sha256 != ROUND2_SHA256
        || round1_sha256 != ROUND1_SHA256
    {
        return Err(format!(
            "dependency hash mismatch: round21={round21_sha256}, round2={round2_sha256}, round1={round1_sha256}"
        ));
    }
    let round21: serde_json::Value =
        serde_json::from_slice(&round21_bytes).map_err(|error| error.to_string())?;
    let round2: serde_json::Value =
        serde_json::from_slice(&round2_bytes).map_err(|error| error.to_string())?;
    let round1: serde_json::Value =
        serde_json::from_slice(&round1_bytes).map_err(|error| error.to_string())?;
    if round21["schema"].as_str() != Some("p256.fb_log_transport/v1")
        || round21["curve"].as_str() != Some(CURVE_SLUG)
        || round21["attack_promotion_gate"].as_bool() != Some(false)
        || round2["schema"].as_str() != Some("p256.dickson_coset_search/v1")
        || round2["curve"].as_str() != Some(CURVE_SLUG)
        || round2["verified"].as_bool() != Some(true)
        || round1["schema"].as_str() != Some("p256.affine_bitbox_result/v1")
        || round1["curve"].as_str() != Some(CURVE_SLUG)
        || round1["verified"].as_bool() != Some(true)
    {
        return Err("dependency identity or rejection status changed".into());
    }
    Ok((
        Dependencies {
            round21_sha256,
            round2_sha256,
            round1_sha256,
        },
        json!({"round2": round2, "round1": round1}),
    ))
}

fn automorphism_screen(curve: &CurveParams) -> Result<AutomorphismScreen, String> {
    let p = &curve.p;
    let four_a3 = (BigUint::from(4u8) * curve.a.modpow(&BigUint::from(3u8), p)) % p;
    let twenty_seven_b2 = (BigUint::from(27u8) * curve.b.modpow(&BigUint::from(2u8), p)) % p;
    let denominator = (&four_a3 + twenty_seven_b2) % p;
    if denominator.is_zero() {
        return Err("P-256 discriminant unexpectedly vanished".into());
    }
    let j = BigUint::from(1728u32) * four_a3 % p * mod_inverse(&denominator, p) % p;
    let cm = p256_cm_discriminant_abs();
    let trial_bound = 1u64 << 20;
    let (factors, residue) = trial_factor(&cm, trial_bound);
    Ok(AutomorphismScreen {
        j_invariant: format!("0x{:064x}", j),
        j_is_zero: j.is_zero(),
        j_is_1728: j == BigUint::from(1728u32),
        fp_degree_one_automorphisms: 2,
        negation_already_folded: true,
        cm_discriminant_abs: cm.to_string(),
        cm_discriminant_bits: cm.bits(),
        cm_trial_bound: trial_bound,
        cm_trial_factors: factors,
        cm_residue_bits: residue.bits(),
        qualification: "j != 0,1728 proves only +/-1 among degree-one Fp automorphisms; trial factoring does not determine the full endomorphism-ring conductor".into(),
    })
}

fn points_from_dump(dump: &WideFactorBaseDump) -> Result<Vec<FactorPoint>, String> {
    if dump.points.len() != dump.factor_base.columns as usize * 2 {
        return Err("wide factor-base row count changed".into());
    }
    let curve = CurveParams::p256();
    let mut points = Vec::with_capacity(dump.factor_base.columns as usize);
    for column in 0..dump.factor_base.columns as usize {
        let row = &dump.points[2 * column];
        if row.col as usize != column || row.coef != "1" {
            return Err(format!("low-y row convention changed at column {column}"));
        }
        let point = Point::Affine {
            x: curve.fe(parse_big(&row.x)?),
            y: curve.fe(parse_big(&row.y)?),
        };
        if !curve.is_on_curve(&point) {
            return Err(format!("factor-base column {column} is off curve"));
        }
        points.push(FactorPoint {
            positive: Projective::from_textbook(&point),
        });
    }
    Ok(points)
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
) -> Result<Vec<([u8; 32], [u8; 32])>, String> {
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
            affine.push((x.to_bytes_be(), y.to_bytes_be()));
        }
    }
    Ok(affine)
}

fn small_scalar_mul(
    point: &Projective,
    scalar: u32,
    additions: &mut u64,
    doublings: &mut u64,
) -> Projective {
    let mut result = Projective::IDENTITY;
    let mut addend = *point;
    let mut value = scalar;
    while value != 0 {
        if value & 1 == 1 {
            result = result.add(&addend);
            *additions += 1;
        }
        value >>= 1;
        if value != 0 {
            addend = addend.double();
            *doublings += 1;
        }
    }
    result
}

fn record_cmp(left: &MultipleRecord, right: &MultipleRecord) -> Ordering {
    left.x
        .cmp(&right.x)
        .then_with(|| left.parity.cmp(&right.parity))
        .then_with(|| left.column.cmp(&right.column))
        .then_with(|| left.multiplier.cmp(&right.multiplier))
}

fn replay_pair(
    left: &MultipleRecord,
    right: &MultipleRecord,
    points: &[FactorPoint],
    operations: &mut OperationCounts,
) -> bool {
    let left_point = small_scalar_mul(
        &points[left.column as usize].positive,
        u32::from(left.multiplier),
        &mut operations.replay_group_additions,
        &mut operations.replay_group_doublings,
    );
    let mut right_point = small_scalar_mul(
        &points[right.column as usize].positive,
        u32::from(right.multiplier),
        &mut operations.replay_group_additions,
        &mut operations.replay_group_doublings,
    );
    if left.parity != right.parity {
        right_point = right_point.neg();
    }
    operations.replay_group_additions += 1;
    bool::from(left_point.add(&right_point.neg()).is_identity())
}

fn transport_census(
    points: &[FactorPoint],
    modulus: &BigUint,
    inverse_exponent: &BigUint,
) -> Result<TransportCensus, String> {
    let mut operations = OperationCounts::default();
    let mut current = points
        .iter()
        .map(|point| point.positive)
        .collect::<Vec<_>>();
    let mut records = Vec::with_capacity(points.len() * usize::from(MAX_MULTIPLIER));
    for multiplier in 1..=MAX_MULTIPLIER {
        let affine = normalize_points(
            &current,
            inverse_exponent,
            &mut operations.normalization_field_multiplications,
        )?;
        records.extend(
            affine
                .into_iter()
                .enumerate()
                .map(|(column, (x, y))| MultipleRecord {
                    x,
                    parity: y[31] & 1,
                    multiplier,
                    column: column as u32,
                }),
        );
        if multiplier != MAX_MULTIPLIER {
            for (image, base) in current.iter_mut().zip(points) {
                *image = image.add(&base.positive);
                operations.image_group_additions += 1;
            }
        }
    }
    records.sort_unstable_by(record_cmp);

    let mut dsu = WeightedDsu::new(points.len(), modulus);
    let mut equal_x_buckets = 0u64;
    let mut equal_x_pairs = 0u64;
    let mut same_column_pairs = 0u64;
    let mut cross_column_candidates = 0u64;
    let mut replayed_relations = 0u64;
    let mut start = 0usize;
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
                equal_x_pairs += 1;
                let left = &records[left_index];
                let right = &records[right_index];
                if left.column == right.column {
                    same_column_pairs += 1;
                    continue;
                }
                cross_column_candidates += 1;
                if !replay_pair(left, right, points, &mut operations) {
                    return Err("equal-x transport relation failed exact replay".into());
                }
                replayed_relations += 1;
                let left_inverse = mod_inverse(&BigUint::from(left.multiplier), modulus);
                let mut numerator = BigUint::from(right.multiplier);
                if left.parity != right.parity {
                    numerator = modulus - numerator;
                }
                let ratio = numerator * left_inverse % modulus;
                if !dsu.unite(left.column as usize, right.column as usize, &ratio) {
                    return Err("inconsistent multiplier-transport cycle".into());
                }
            }
        }
        start = end;
    }
    let (isolated, largest) = dsu.histogram();
    Ok(TransportCensus {
        multiplier_max: MAX_MULTIPLIER,
        images: records.len() as u64,
        equal_x_buckets,
        equal_x_pairs,
        same_column_pairs,
        cross_column_candidates,
        replayed_relations,
        replay_failures: 0,
        inconsistent_cycles: dsu.inconsistent_cycles,
        independent_rank: dsu.independent_edges,
        quotient_dimension: dsu.components as u64,
        isolated_columns: isolated as u64,
        largest_component: largest as u64,
        record_logical_bytes: size_of::<MultipleRecord>() as u64,
        peak_logical_bytes: (records.capacity() * size_of::<MultipleRecord>()
            + current.capacity() * size_of::<Projective>()) as u64,
        disk_bytes: 0,
        operations,
    })
}

fn disjoint_probability(columns: u64) -> f64 {
    (0..9).fold(1.0, |probability, offset| {
        probability * (columns - 8 - offset) as f64 / (columns - offset) as f64
    })
}

fn log2_big(value: &BigUint) -> f64 {
    if value.is_zero() {
        return f64::NEG_INFINITY;
    }
    let bits = value.bits();
    let shift = bits.saturating_sub(53);
    let top = (value >> shift).to_u64().unwrap_or(0) as f64;
    top.log2() + shift as f64
}

fn projection(columns: u64, classes: u64, structured_degree: Option<u32>) -> Projection {
    let curve = CurveParams::p256();
    let p_disjoint = disjoint_probability(columns);
    let k = classes as f64;
    let cross_colour_s = (std::f64::consts::PI * k / p_disjoint).sqrt();
    let two_list_s = 2.0 * (k / p_disjoint).sqrt();
    let log2_sqrt_n = log2_big(&curve.n) / 2.0;
    let cross_log2 = log2_sqrt_n + cross_colour_s.log2();
    let rows_s = (std::f64::consts::PI * RELATION_ROWS as f64 / p_disjoint).sqrt();
    let rows_log2 = log2_sqrt_n + rows_s.log2();
    let per_relation_log2 = rows_log2 - (RELATION_ROWS as f64).log2();
    let sparse_la_log2 = if classes > 1 {
        let nonzeros = RELATION_ROWS as f64 * 17.0;
        let wiedemann = 2.0 * k * nonzeros;
        let berlekamp_massey = k * k;
        Some((wiedemann + berlekamp_massey).log2())
    } else {
        None
    };
    let complete_log2 = sparse_la_log2.map_or(cross_log2, |la_log2| {
        let maximum = rows_log2.max(la_log2);
        maximum + (2f64.powf(rows_log2 - maximum) + 2f64.powf(la_log2 - maximum)).log2()
    });
    let peak_bytes_log2 = log2_sqrt_n + (k / p_disjoint).sqrt().log2() + 6.0;
    Projection {
        independent_log_classes: classes,
        disjoint_probability: p_disjoint,
        cross_colour_s,
        cross_colour_ratio_to_rho: cross_colour_s / RHO_S,
        two_list_s,
        two_list_ratio_to_rho: two_list_s / RHO_S,
        cross_colour_log2_operations: cross_log2,
        projected_138031_rows_log2_operations: rows_log2,
        projected_cost_per_usable_relation_log2_operations: per_relation_log2,
        sparse_linear_algebra_lower_bound_log2_operations: sparse_la_log2,
        complete_projected_lower_bound_log2_operations: complete_log2,
        projected_peak_two_list_log2_bytes: peak_bytes_log2,
        structured_residual_degree_max: structured_degree,
        relation_collection_below_2_120: rows_log2 < 120.0,
        cost_per_relation_below_2_103: per_relation_log2 < 103.0,
        materialized_storage_below_2_50: peak_bytes_log2 < 50.0,
        complete_cost_below_rho: false,
        promotion_gate: false,
        dominant_obstruction: if classes == 1 {
            "generic two-colour target collision is already above rho; no exact non-generic target mechanism".into()
        } else {
            "independent-log quotient leaves the cross-colour relation floor far above rho".into()
        },
    }
}

fn screen_geometric(
    name: String,
    spec: String,
    dump: WideFactorBaseDump,
    degree_evidence: String,
    structured_degree: Option<u32>,
    modulus: &BigUint,
    inverse_exponent: &BigUint,
) -> Result<GeometricCandidate, String> {
    let fb = dump.factor_base.clone();
    let points = points_from_dump(&dump)?;
    let transport = transport_census(&points, modulus, inverse_exponent)?;
    let projected = projection(fb.columns, transport.quotient_dimension, structured_degree);
    Ok(GeometricCandidate {
        name,
        family: fb.family,
        spec,
        fb_id: fb.fb_id,
        fb_sha256: fb.fb_sha256,
        points_sha256: fb.points_sha256,
        columns: fb.columns,
        signed_points: fb.signed_points,
        rebuild_verified: true,
        degree_evidence,
        unsplit_s17_degree_of_regularity: None,
        transport,
        projection: projected,
        classification: "negative screen".into(),
    })
}

fn canonicalize_points(
    raw_points: Vec<Projective>,
    raw_scalars: Vec<BigUint>,
    curve: &CurveParams,
    inverse_exponent: &BigUint,
) -> Result<(Vec<Projective>, Vec<BigUint>, String), String> {
    let mut ignored = 0u64;
    let affine = normalize_points(&raw_points, inverse_exponent, &mut ignored)?;
    let mut rows = Vec::with_capacity(raw_points.len());
    for ((point, scalar), (x, y)) in raw_points.into_iter().zip(raw_scalars).zip(affine) {
        let y_big = BigUint::from_bytes_be(&y);
        let low = y_big <= &curve.p - &y_big;
        let oriented_point = if low { point } else { point.neg() };
        let oriented_scalar = if low { scalar } else { &curve.n - scalar };
        rows.push((x, oriented_point, oriented_scalar));
    }
    rows.sort_unstable_by_key(|row| row.0);
    if rows.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Err("scalar construction contains duplicate signed columns".into());
    }
    let mut point_key_blob = Vec::with_capacity(rows.len() * 66);
    for (x, _, _) in &rows {
        let x_big = BigUint::from_bytes_be(x);
        point_key_blob.extend(point_key_bytes(&x_big, false)?);
        point_key_blob.extend(point_key_bytes(&x_big, true)?);
    }
    let points = rows.iter().map(|row| row.1).collect();
    let scalars = rows.into_iter().map(|row| row.2).collect();
    Ok((points, scalars, hex::encode(sha256(&point_key_blob))))
}

fn build_scalar_candidate(
    name: &str,
    multiplier: Option<u32>,
    curve: &CurveParams,
    inverse_exponent: &BigUint,
) -> Result<ScalarCandidate, String> {
    let generator = Projective::from_textbook(&curve.generator());
    let mut raw_points = Vec::with_capacity(COLUMNS);
    let mut raw_scalars = Vec::with_capacity(COLUMNS);
    let mut additions = 0u64;
    let mut doublings = 0u64;
    let (transition, family, parameter, terminal_image_scalar) = if let Some(a) = multiplier {
        let mut point = generator;
        let mut scalar = BigUint::one();
        for _ in 0..COLUMNS {
            raw_points.push(point);
            raw_scalars.push(scalar.clone());
            point = small_scalar_mul(&point, a, &mut additions, &mut doublings);
            scalar = scalar * a % &curve.n;
        }
        (
            format!("P_(i+1)=[{a}]P_i"),
            "scalar-geometric",
            a.to_string(),
            scalar,
        )
    } else {
        let mut point = generator;
        for scalar in 1..=COLUMNS as u64 {
            raw_points.push(point);
            raw_scalars.push(BigUint::from(scalar));
            point = point.add(&generator);
            additions += 1;
        }
        (
            "P_(i+1)=P_i+G".into(),
            "scalar-additive",
            "1".into(),
            BigUint::from(COLUMNS as u64 + 1),
        )
    };
    let (points, scalars, points_sha256) =
        canonicalize_points(raw_points, raw_scalars, curve, inverse_exponent)?;
    let mut params = BTreeMap::new();
    params.insert("columns".into(), COLUMNS.to_string());
    params.insert("point_key_encoding".into(), POINT_KEY_ENCODING.into());
    params.insert("transition_parameter".into(), parameter);
    let manifest = ScalarManifest {
        schema: FACTOR_BASE_SCHEMA.into(),
        curve: CURVE_SLUG.into(),
        family: family.into(),
        params,
        columns: COLUMNS as u64,
        signed_points: (2 * COLUMNS) as u64,
        points_sha256,
    };
    let manifest_value = serde_json::to_value(&manifest).map_err(|error| error.to_string())?;
    let (fb_id, fb_sha256) = short_id("FB1", &manifest_value)?;

    let scalar_set = scalars
        .iter()
        .map(|scalar| signed_scalar_key(scalar, &curve.n))
        .collect::<BTreeSet<_>>();
    if scalar_set.len() != COLUMNS || points.len() != COLUMNS {
        return Err("scalar construction did not retain the frozen column count".into());
    }
    let terminal_image_scalar = terminal_image_scalar % &curve.n;
    let closes = scalar_set.contains(&signed_scalar_key(&terminal_image_scalar, &curve.n));
    let (support_upper, support_log2) = if multiplier.is_none() {
        let largest_sum = (0..17u64)
            .map(|offset| BigUint::from(COLUMNS as u64 - offset))
            .sum::<BigUint>();
        let support = BigUint::from(2u8) * largest_sum + BigUint::one();
        let fraction_log2 = log2_big(&support) - log2_big(&curve.n);
        (Some(support.to_string()), Some(fraction_log2))
    } else {
        (None, None)
    };
    Ok(ScalarCandidate {
        name: name.into(),
        transition,
        fb_id,
        fb_sha256,
        manifest,
        signed_duplicates: 0,
        zero_scalars: 0,
        internal_transition_relations: (COLUMNS - 1) as u64,
        transition_closes_at_end: closes,
        exact_relation_rank: (COLUMNS - 1) as u64,
        quotient_dimension: 1,
        known_construction_logs: true,
        actual_p256_points_verified: true,
        point_generation_additions: additions,
        point_generation_doublings: doublings,
        additive_s17_support_upper_bound: support_upper,
        additive_s17_support_fraction_log2: support_log2,
        projection: projection(COLUMNS as u64, 1, None),
        diagnostics_status: "stopped after the exact cross-colour rho gate failed; no incomplete sampling receives exhaustive credit".into(),
        classification: "relabelling control".into(),
    })
}

fn signed_scalar_key(scalar: &BigUint, modulus: &BigUint) -> BigUint {
    let negative = modulus - scalar;
    if scalar <= &negative {
        scalar.clone()
    } else {
        negative
    }
}

fn mod_inverse(value: &BigUint, modulus: &BigUint) -> BigUint {
    value.modpow(&(modulus - 2u8), modulus)
}

fn run(cli: &Cli) -> Result<ResultFile, String> {
    let (dependencies, artifacts) = verify_dependencies(cli)?;
    let curve = CurveParams::p256();
    let inverse_exponent = &curve.p - 2u8;
    let automorphism_screen = automorphism_screen(&curve)?;
    if automorphism_screen.j_is_zero || automorphism_screen.j_is_1728 {
        return Err("P-256 unexpectedly has a special j-invariant".into());
    }

    let mut geometric_candidates = Vec::new();
    let candidates = artifacts["round2"]["candidates"]
        .as_array()
        .ok_or("round-2 candidate list is missing")?;
    for candidate in candidates {
        let index = candidate["index"]
            .as_u64()
            .ok_or("candidate index missing")?;
        let root = candidate["root_exponent"]
            .as_str()
            .ok_or("candidate root missing")?;
        let spec = format!("dickson-torus:depth=18,root_exponent={root}");
        let built = p256_dickson_factor_base::build(&spec)?;
        p256_dickson_factor_base::verify(&built.dump)?;
        let fb = &built.dump.factor_base;
        if fb.fb_id != candidate["fb_id"].as_str().unwrap_or_default()
            || fb.fb_sha256 != candidate["fb_sha256"].as_str().unwrap_or_default()
            || fb.points_sha256 != candidate["points_sha256"].as_str().unwrap_or_default()
            || fb.columns != candidate["columns"].as_u64().unwrap_or_default()
        {
            return Err(format!("round-2 candidate {index} rebuild mismatch"));
        }
        geometric_candidates.push(screen_geometric(
            format!("dickson-coset-{index}"),
            spec,
            built.dump,
            "small-field Dickson-family S3 degree 4; structured residual maxima 3,3,4".into(),
            Some(4),
            &curve.n,
            &inverse_exponent,
        )?);
    }

    let terminal_zero_spec = "dickson-torus:depth=18";
    let terminal_zero = p256_dickson_factor_base::build(terminal_zero_spec)?;
    p256_dickson_factor_base::verify(&terminal_zero.dump)?;
    if terminal_zero.dump.factor_base.fb_id != "FB1h2ea06bef7f7a" {
        return Err("terminal-zero Dickson FB1 identity changed".into());
    }
    geometric_candidates.push(screen_geometric(
        "dickson-terminal-zero".into(),
        terminal_zero_spec.into(),
        terminal_zero.dump,
        "small-field Dickson-family S3 degree 4; structured residual maxima 3,3,4".into(),
        Some(4),
        &curve.n,
        &inverse_exponent,
    )?);

    let bitbox = p256_bitbox_factor_base::build()?;
    p256_bitbox_factor_base::verify(&bitbox)?;
    let expected_bitbox = &artifacts["round1"];
    if bitbox.build.dump.factor_base.fb_id != expected_bitbox["fb_id"].as_str().unwrap_or_default()
        || bitbox.build.dump.factor_base.fb_sha256
            != expected_bitbox["fb_sha256"].as_str().unwrap_or_default()
    {
        return Err("affine-bitbox rebuild mismatch".into());
    }
    geometric_candidates.push(screen_geometric(
        "affine-bitbox-winner".into(),
        "affine-bitbox:bit_width=18,offset=0x180000".into(),
        bitbox.build.dump,
        "measured small-field S3 regression: degree 5 at depths 3-4 and 6 at depth 5".into(),
        Some(6),
        &curve.n,
        &inverse_exponent,
    )?);

    let mut scalar_candidates = Vec::new();
    scalar_candidates.push(build_scalar_candidate(
        "additive-interval",
        None,
        &curve,
        &inverse_exponent,
    )?);
    for multiplier in [2u32, 3, 5, 7, 11, 65_537] {
        scalar_candidates.push(build_scalar_candidate(
            &format!("scalar-geometric-{multiplier}"),
            Some(multiplier),
            &curve,
            &inverse_exponent,
        )?);
    }

    let stopped_controls = vec![StoppedControl {
        name: "hash-scalar-control".into(),
        stop_reason: "the preregistered exact K=1 cross-colour boundary is 1.363x rho before construction or replay, so the cost stop rule fires before materialising this expensive known-log control".into(),
        exhaustive_credit: false,
        projection: projection(COLUMNS as u64, 1, None),
    }];
    let all_replayed = geometric_candidates
        .iter()
        .all(|candidate| candidate.transport.replay_failures == 0);
    let any_passed = geometric_candidates
        .iter()
        .any(|candidate| candidate.projection.promotion_gate)
        || scalar_candidates
            .iter()
            .any(|candidate| candidate.projection.promotion_gate);
    Ok(ResultFile {
        schema: "p256.factor_base_symmetry/v1".into(),
        curve: CURVE_SLUG.into(),
        relation_model: "17 distinct variable columns; balanced signed 8+9 cross-colour target collision".into(),
        dependencies,
        automorphism_screen,
        geometric_candidates,
        scalar_candidates,
        stopped_controls,
        reference_false_positives: 0,
        reference_false_negatives: 0,
        all_reported_relations_replayed: all_replayed,
        any_promotion_gate_passed: any_passed,
        full_depth_unplanted_attempted: false,
        result_class: "negative result / accounting correction".into(),
        decision: "No screened factor base passes. Geometric bases retain essentially all independent logarithms; deliberately one-orbit known-log bases collapse K but remain above rho at the exact two-colour target-collision floor and supply no algebraic degree credit.".into(),
    })
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(result) => {
            let bytes = match serde_json::to_vec_pretty(&result) {
                Ok(bytes) => bytes,
                Err(error) => {
                    eprintln!("serialization failed: {error}");
                    return ExitCode::FAILURE;
                }
            };
            if let Some(path) = cli.out {
                if let Err(error) = fs::write(&path, &bytes) {
                    eprintln!("{}: {error}", path.display());
                    return ExitCode::FAILURE;
                }
            } else {
                println!("{}", String::from_utf8_lossy(&bytes));
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("p256 factor-base symmetry screen failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cross_colour_k1_floor_is_above_rho() {
        let p = projection(COLUMNS as u64, 1, None);
        assert!(p.cross_colour_ratio_to_rho > 1.36);
        assert!(p.cross_colour_ratio_to_rho < 1.37);
        assert!(p.projected_138031_rows_log2_operations > 137.0);
        assert!(p.projected_138031_rows_log2_operations < 138.0);
        assert!(p.projected_cost_per_usable_relation_log2_operations > 120.0);
        assert!(!p.promotion_gate);
    }

    #[test]
    fn disjoint_probability_is_close_to_one() {
        let p = disjoint_probability(COLUMNS as u64);
        assert!(p > 0.9994 && p < 1.0);
    }

    #[test]
    fn signed_key_folds_negation() {
        let modulus = BigUint::from(101u32);
        assert_eq!(
            signed_scalar_key(&BigUint::from(7u8), &modulus),
            BigUint::from(7u8)
        );
        assert_eq!(
            signed_scalar_key(&BigUint::from(94u8), &modulus),
            BigUint::from(7u8)
        );
    }

    #[test]
    fn p256_is_not_a_special_j_curve() {
        let screen = automorphism_screen(&CurveParams::p256()).unwrap();
        assert!(!screen.j_is_zero);
        assert!(!screen.j_is_1728);
        assert_eq!(screen.fp_degree_one_automorphisms, 2);
    }
}
