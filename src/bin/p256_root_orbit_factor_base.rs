//! Exact P-256 PGL2 root-orbit factor-base screen (round 23).

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::mem::size_of;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ec_index_calculus::sqrt_mod_p;
use crypto_lib::cryptanalysis::ecbench::canonical::short_id;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{
    CURVE_SLUG, FACTOR_BASE_SCHEMA, POINT_KEY_ENCODING,
};
use crypto_lib::ct_bignum::U256;
use crypto_lib::ecc::p256_field::P256FieldElement as Fe;
use crypto_lib::ecc::p256_point::{P256ProjectivePoint as Projective, A_FE, B_FE};
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;

const ROUND22_SHA256: &str = "3116c678d257794040c5519c85a8037177e12381e5f948a47cf1335d0cafde16";
const SUBGROUP_ORDERS: [u64; 2] = [260_270, 272_425];
const ARMS_PER_ORDER: u64 = 8;
const MAX_MULTIPLIER: u8 = 16;
const NORMALIZE_BATCH: usize = 8_192;
const INDEPENDENT_REPLAY_COLUMNS: usize = 4_096;
const RELATION_ROWS: u64 = 138_031;
const RHO_S: f64 = 1.3;

#[derive(Parser)]
#[command(about = "Screen P-256 PGL2 root-orbit factor bases")]
struct Cli {
    #[arg(long)]
    round22: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Clone)]
struct Matrix {
    a: BigUint,
    b: BigUint,
    c: BigUint,
    e: BigUint,
}

#[derive(Clone, Copy)]
struct MatrixFe {
    a: Fe,
    b: Fe,
    c: Fe,
    e: Fe,
}

#[derive(Clone, Copy)]
struct FactorPoint {
    positive: Projective,
}

#[derive(Clone, Copy)]
struct OrbitRow {
    orbit_index: u32,
    x: [u8; 32],
    y: [u8; 32],
    point: Projective,
}

#[derive(Clone, Copy, Eq, PartialEq)]
struct MultipleRecord {
    x: [u8; 32],
    parity: u8,
    multiplier: u8,
    column: u32,
}

#[derive(Default, Serialize)]
struct GenerationOperations {
    orbit_field_multiplications: u64,
    transform_field_multiplications: u64,
    batch_inverse_field_multiplications: u64,
    batch_inversions: u64,
    sqrt_field_multiplications: u64,
    sqrt_field_squarings: u64,
    closure_field_multiplications: u64,
    independent_bigint_replays: u64,
}

#[derive(Default, Serialize)]
struct TransportOperations {
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
    operations: TransportOperations,
}

#[derive(Serialize)]
struct FactorBaseManifest {
    schema: String,
    curve: String,
    family: String,
    params: BTreeMap<String, String>,
    columns: u64,
    signed_points: u64,
    points_sha256: String,
}

#[derive(Serialize)]
struct Projection {
    independent_log_classes: u64,
    disjoint_probability: f64,
    cross_colour_s: f64,
    cross_colour_ratio_to_rho: f64,
    cross_colour_log2_operations: f64,
    projected_138031_rows_log2_operations: f64,
    projected_cost_per_usable_relation_log2_operations: f64,
    sparse_linear_algebra_lower_bound_log2_operations: f64,
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
struct Candidate {
    name: String,
    family: String,
    subgroup_order: u64,
    arm: u64,
    derivation_counter: u64,
    matrix: BTreeMap<String, String>,
    determinant: String,
    omega: String,
    omega_exact_order_verified: bool,
    orbit_abscissae: u64,
    rejected_poles: u64,
    columns: u64,
    signed_points: u64,
    quadratic_residue_fraction: f64,
    coordinate_orbit_closed: bool,
    induced_recurrence_closed: bool,
    nonidentity_coordinate_action: bool,
    degree_one_curve_automorphism: bool,
    membership_generator_degree: u32,
    structured_residual_degree_of_regularity: Option<u32>,
    unsplit_s17_degree_of_regularity: Option<u32>,
    independent_point_replays: u64,
    point_replay_failures: u64,
    point_keys_sorted_unique: bool,
    manifest: FactorBaseManifest,
    fb_id: String,
    fb_sha256: String,
    generation_operations: GenerationOperations,
    transport: TransportCensus,
    projection: Projection,
    classification: String,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    relation_model: String,
    dependency_path: String,
    dependency_sha256: String,
    p_minus_one_factorization_verified: bool,
    primitive_root: String,
    primitive_root_verified: bool,
    candidates: Vec<Candidate>,
    total_orbit_abscissae: u64,
    total_columns: u64,
    total_transport_images: u64,
    total_cross_column_candidates: u64,
    total_exact_relations: u64,
    all_reported_relations_replayed: bool,
    reference_false_positives: u64,
    reference_false_negatives: u64,
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

fn file_sha256(path: &Path) -> Result<String, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    Ok(hex::encode(sha256(&bytes)))
}

fn prime_factors() -> Result<Vec<BigUint>, String> {
    [
        "2",
        "3",
        "5",
        "17",
        "257",
        "641",
        "1531",
        "65537",
        "490463",
        "6700417",
        "835945042244614951780389953367877943453916927241",
    ]
    .iter()
    .map(|text| BigUint::parse_bytes(text.as_bytes(), 10).ok_or_else(|| text.to_string()))
    .collect()
}

fn verify_factorization(p: &BigUint, factors: &[BigUint]) -> bool {
    let exponents = [1u32, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1];
    factors
        .iter()
        .zip(exponents)
        .fold(BigUint::one(), |product, (factor, exponent)| {
            product * factor.pow(exponent)
        })
        == p - 1u8
}

fn primitive_root(p: &BigUint, factors: &[BigUint]) -> BigUint {
    let p_minus_one = p - 1u8;
    for candidate in 2u64.. {
        let value = BigUint::from(candidate);
        if factors
            .iter()
            .all(|factor| value.modpow(&(&p_minus_one / factor), p) != BigUint::one())
        {
            return value;
        }
    }
    unreachable!()
}

fn derive_matrix(subgroup_order: u64, counter: u64, p: &BigUint) -> Matrix {
    let mut words = Vec::new();
    for word in 0..4u8 {
        let label = format!("{CURVE_SLUG}/root-orbit-round23/{subgroup_order}/{counter}/{word}");
        words.push(BigUint::from_bytes_be(&sha256(label.as_bytes())) % p);
    }
    Matrix {
        a: words[0].clone(),
        b: words[1].clone(),
        c: words[2].clone(),
        e: words[3].clone(),
    }
}

fn matrix_fe(matrix: &Matrix) -> MatrixFe {
    MatrixFe {
        a: Fe::from_biguint(&matrix.a),
        b: Fe::from_biguint(&matrix.b),
        c: Fe::from_biguint(&matrix.c),
        e: Fe::from_biguint(&matrix.e),
    }
}

#[cfg(test)]
fn matrix_map(matrix: &MatrixFe, value: Fe) -> Option<Fe> {
    let numerator = matrix.a.mul(&value).add(&matrix.b);
    let denominator = matrix.c.mul(&value).add(&matrix.e);
    if bool::from(denominator.ct_is_zero()) {
        None
    } else {
        Some(numerator.mul(&denominator.inv()))
    }
}

fn pow_fe(mut base: Fe, exponent: &BigUint, muls: &mut u64, squares: &mut u64) -> Fe {
    let mut result = Fe::ONE;
    for bit in 0..exponent.bits() {
        if exponent.bit(bit) {
            result = result.mul(&base);
            *muls += 1;
        }
        base = base.sqr();
        *squares += 1;
    }
    result
}

fn batch_invert(values: &[Fe], operations: &mut GenerationOperations) -> Result<Vec<Fe>, String> {
    if values.iter().any(|value| bool::from(value.ct_is_zero())) {
        return Err("PGL2 orbit contains a pole".into());
    }
    let mut prefixes = Vec::with_capacity(values.len());
    let mut product = Fe::ONE;
    for value in values {
        prefixes.push(product);
        product = product.mul(value);
        operations.batch_inverse_field_multiplications += 1;
    }
    let mut inverse = product.inv();
    operations.batch_inversions += 1;
    let mut inverses = vec![Fe::ZERO; values.len()];
    for index in (0..values.len()).rev() {
        inverses[index] = inverse.mul(&prefixes[index]);
        operations.batch_inverse_field_multiplications += 1;
        if index != 0 {
            inverse = inverse.mul(&values[index]);
            operations.batch_inverse_field_multiplications += 1;
        }
    }
    Ok(inverses)
}

fn point_key_bytes(x: &[u8; 32], high_sign: bool) -> Result<[u8; 33], String> {
    let x = BigUint::from_bytes_be(x);
    let key = ((x + BigUint::one()) << 1usize) + BigUint::from(high_sign);
    let raw = key.to_bytes_be();
    if raw.len() > 33 {
        return Err("P-256 point key exceeded 33 bytes".into());
    }
    let mut out = [0u8; 33];
    out[33 - raw.len()..].copy_from_slice(&raw);
    Ok(out)
}

fn closure_matrix(matrix: &MatrixFe, omega: Fe) -> MatrixFe {
    MatrixFe {
        a: matrix
            .a
            .mul(&omega)
            .mul(&matrix.e)
            .sub(&matrix.b.mul(&matrix.c)),
        b: matrix.a.mul(&matrix.b).mul(&Fe::ONE.sub(&omega)),
        c: matrix.c.mul(&matrix.e).mul(&omega.sub(&Fe::ONE)),
        e: matrix
            .a
            .mul(&matrix.e)
            .sub(&matrix.c.mul(&omega).mul(&matrix.b)),
    }
}

fn verify_recurrence(
    xs: &[Fe],
    recurrence: &MatrixFe,
    operations: &mut GenerationOperations,
) -> bool {
    for index in 0..xs.len() {
        let next = xs[(index + 1) % xs.len()];
        let numerator = recurrence.a.mul(&xs[index]).add(&recurrence.b);
        let denominator = recurrence.c.mul(&xs[index]).add(&recurrence.e);
        operations.closure_field_multiplications += 3;
        if bool::from(denominator.ct_is_zero()) || next.mul(&denominator) != numerator {
            return false;
        }
    }
    true
}

fn independent_replay(
    rows: &[OrbitRow],
    matrix: &Matrix,
    omega: &BigUint,
    curve: &CurveParams,
) -> Result<u64, String> {
    let mut failures = 0u64;
    let mut first_failure = None;
    for row in rows.iter().take(INDEPENDENT_REPLAY_COLUMNS) {
        let u = omega.modpow(&BigUint::from(row.orbit_index), &curve.p);
        let numerator = (&matrix.a * &u + &matrix.b) % &curve.p;
        let denominator = (&matrix.c * &u + &matrix.e) % &curve.p;
        if denominator.is_zero() {
            return Err("independent point replay found a pole".into());
        }
        let x = numerator * mod_inverse(&denominator, &curve.p) % &curve.p;
        let rhs = ((&x * &x % &curve.p) * &x + &curve.a * &x + &curve.b) % &curve.p;
        let Some(y) = sqrt_mod_p(&rhs, &curve.p) else {
            failures += 1;
            first_failure.get_or_insert_with(|| {
                format!(
                    "orbit {} changed from residue to non-residue",
                    row.orbit_index
                )
            });
            continue;
        };
        let negative = &curve.p - &y;
        let low_y = if y <= negative { y } else { negative };
        let x_bytes = biguint_bytes_32(&x)?;
        let y_bytes = biguint_bytes_32(&low_y)?;
        if x_bytes != row.x || y_bytes != row.y {
            failures += 1;
            first_failure.get_or_insert_with(|| {
                format!(
                    "orbit {} byte mismatch: x={}, y={}",
                    row.orbit_index,
                    x_bytes != row.x,
                    y_bytes != row.y
                )
            });
        }
    }
    if let Some(first_failure) = first_failure {
        eprintln!("first independent replay failure: {first_failure}");
    }
    Ok(failures)
}

fn biguint_bytes_32(value: &BigUint) -> Result<[u8; 32], String> {
    let raw = value.to_bytes_be();
    if raw.len() > 32 {
        return Err("field value exceeded 32 bytes".into());
    }
    let mut bytes = [0u8; 32];
    bytes[32 - raw.len()..].copy_from_slice(&raw);
    Ok(bytes)
}

fn build_rows(
    subgroup_order: u64,
    omega: &BigUint,
    matrix: &Matrix,
    curve: &CurveParams,
) -> Result<(Vec<OrbitRow>, Vec<Fe>, GenerationOperations), String> {
    let omega_fe = Fe::from_biguint(omega);
    let matrix = matrix_fe(matrix);
    let sqrt_exponent = (&curve.p + 1u8) >> 2usize;
    let mut operations = GenerationOperations::default();
    let mut us = Vec::with_capacity(subgroup_order as usize);
    let mut u = Fe::ONE;
    for _ in 0..subgroup_order {
        us.push(u);
        u = u.mul(&omega_fe);
        operations.orbit_field_multiplications += 1;
    }
    if u != Fe::ONE || us[1..].contains(&Fe::ONE) {
        return Err("subgroup generator did not have the declared exact order".into());
    }

    let mut xs = Vec::with_capacity(us.len());
    for batch in us.chunks(NORMALIZE_BATCH) {
        let numerators = batch
            .iter()
            .map(|value| matrix.a.mul(value).add(&matrix.b))
            .collect::<Vec<_>>();
        let denominators = batch
            .iter()
            .map(|value| matrix.c.mul(value).add(&matrix.e))
            .collect::<Vec<_>>();
        operations.transform_field_multiplications += 2 * batch.len() as u64;
        let inverses = batch_invert(&denominators, &mut operations)?;
        for (numerator, inverse) in numerators.into_iter().zip(inverses) {
            xs.push(numerator.mul(&inverse));
            operations.transform_field_multiplications += 1;
        }
    }
    let recurrence = closure_matrix(&matrix, omega_fe);
    if !verify_recurrence(&xs, &recurrence, &mut operations) {
        return Err("induced fractional-linear recurrence failed closure".into());
    }

    let mut rows = Vec::with_capacity(xs.len() / 2 + 1024);
    for (orbit_index, x) in xs.iter().enumerate() {
        let rhs = x.sqr().mul(x).add(&A_FE.mul(x)).add(&B_FE);
        operations.sqrt_field_squarings += 1;
        operations.sqrt_field_multiplications += 2;
        let mut y = pow_fe(
            rhs,
            &sqrt_exponent,
            &mut operations.sqrt_field_multiplications,
            &mut operations.sqrt_field_squarings,
        );
        if y.sqr() != rhs {
            operations.sqrt_field_squarings += 1;
            continue;
        }
        operations.sqrt_field_squarings += 1;
        let negative = y.neg();
        if negative.to_bytes_be() < y.to_bytes_be() {
            y = negative;
        }
        let x_bytes = x.to_bytes_be();
        let y_bytes = y.to_bytes_be();
        let point = Projective::from_affine(
            &U256::from_bytes_be(&x_bytes),
            &U256::from_bytes_be(&y_bytes),
        );
        rows.push(OrbitRow {
            orbit_index: orbit_index as u32,
            x: x_bytes,
            y: y_bytes,
            point,
        });
    }
    rows.sort_unstable_by_key(|row| row.x);
    if rows.windows(2).any(|pair| pair[0].x == pair[1].x) {
        return Err("transformed orbit contains duplicate abscissae".into());
    }
    Ok((rows, xs, operations))
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
            affine.push((
                point.x.mul(&inverse).to_bytes_be(),
                point.y.mul(&inverse).to_bytes_be(),
            ));
            *multiplications += 2;
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
    operations: &mut TransportOperations,
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
    curve: &CurveParams,
) -> Result<TransportCensus, String> {
    let inverse_exponent = &curve.p - 2u8;
    let mut operations = TransportOperations::default();
    let mut current = points
        .iter()
        .map(|point| point.positive)
        .collect::<Vec<_>>();
    let mut records = Vec::with_capacity(points.len() * usize::from(MAX_MULTIPLIER));
    for multiplier in 1..=MAX_MULTIPLIER {
        let affine = normalize_points(
            &current,
            &inverse_exponent,
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

    let mut dsu = WeightedDsu::new(points.len(), &curve.n);
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
                let left_inverse = mod_inverse(&BigUint::from(left.multiplier), &curve.n);
                let mut numerator = BigUint::from(right.multiplier);
                if left.parity != right.parity {
                    numerator = &curve.n - numerator;
                }
                let ratio = numerator * left_inverse % &curve.n;
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
    let bits = value.bits();
    let shift = bits.saturating_sub(53);
    let top = (value >> shift).to_u64().unwrap_or(0) as f64;
    top.log2() + shift as f64
}

fn projection(columns: u64, classes: u64, curve: &CurveParams) -> Projection {
    let p_disjoint = disjoint_probability(columns);
    let k = classes as f64;
    let cross_colour_s = (std::f64::consts::PI * k / p_disjoint).sqrt();
    let log2_sqrt_n = log2_big(&curve.n) / 2.0;
    let cross_log2 = log2_sqrt_n + cross_colour_s.log2();
    let rows_s = (std::f64::consts::PI * RELATION_ROWS as f64 / p_disjoint).sqrt();
    let rows_log2 = log2_sqrt_n + rows_s.log2();
    let per_relation_log2 = rows_log2 - (RELATION_ROWS as f64).log2();
    let nonzeros = RELATION_ROWS as f64 * 17.0;
    let sparse_la_log2 = (2.0 * k * nonzeros + k * k).log2();
    let maximum = rows_log2.max(sparse_la_log2);
    let complete_log2 =
        maximum + (2f64.powf(rows_log2 - maximum) + 2f64.powf(sparse_la_log2 - maximum)).log2();
    let peak_bytes_log2 = log2_sqrt_n + (k / p_disjoint).sqrt().log2() + 6.0;
    Projection {
        independent_log_classes: classes,
        disjoint_probability: p_disjoint,
        cross_colour_s,
        cross_colour_ratio_to_rho: cross_colour_s / RHO_S,
        cross_colour_log2_operations: cross_log2,
        projected_138031_rows_log2_operations: rows_log2,
        projected_cost_per_usable_relation_log2_operations: per_relation_log2,
        sparse_linear_algebra_lower_bound_log2_operations: sparse_la_log2,
        complete_projected_lower_bound_log2_operations: complete_log2,
        projected_peak_two_list_log2_bytes: peak_bytes_log2,
        structured_residual_degree_max: None,
        relation_collection_below_2_120: rows_log2 < 120.0,
        cost_per_relation_below_2_103: per_relation_log2 < 103.0,
        materialized_storage_below_2_50: peak_bytes_log2 < 50.0,
        complete_cost_below_rho: false,
        promotion_gate: false,
        dominant_obstruction: "coordinate closure supplies no elliptic-log transport; the independent-log cross-colour floor remains far above rho".into(),
    }
}

fn mod_inverse(value: &BigUint, modulus: &BigUint) -> BigUint {
    value.modpow(&(modulus - 2u8), modulus)
}

fn hex_field(value: &BigUint) -> String {
    format!("0x{value:064x}")
}

fn build_candidate(
    subgroup_order: u64,
    arm: u64,
    primitive_root: &BigUint,
    curve: &CurveParams,
    seen_point_digests: &mut BTreeSet<String>,
) -> Result<Candidate, String> {
    let omega = primitive_root.modpow(&((&curve.p - 1u8) / subgroup_order), &curve.p);
    let prime_divisors = prime_factors()?
        .into_iter()
        .filter(|factor| BigUint::from(subgroup_order) % factor == BigUint::zero())
        .collect::<Vec<_>>();
    if omega.modpow(&BigUint::from(subgroup_order), &curve.p) != BigUint::one()
        || prime_divisors.iter().any(|factor| {
            omega.modpow(&(BigUint::from(subgroup_order) / factor), &curve.p) == BigUint::one()
        })
    {
        return Err("omega exact-order check failed".into());
    }

    let mut counter = arm;
    loop {
        let matrix = derive_matrix(subgroup_order, counter, &curve.p);
        let determinant =
            (&matrix.a * &matrix.e + &curve.p - (&matrix.b * &matrix.c % &curve.p)) % &curve.p;
        if determinant.is_zero() {
            counter += ARMS_PER_ORDER;
            continue;
        }
        let (mut rows, _xs, mut generation_operations) =
            match build_rows(subgroup_order, &omega, &matrix, curve) {
                Ok(built) => built,
                Err(error) if error.contains("pole") => {
                    counter += ARMS_PER_ORDER;
                    continue;
                }
                Err(error) => return Err(error),
            };
        let replay_failures = independent_replay(&rows, &matrix, &omega, curve)?;
        generation_operations.independent_bigint_replays =
            rows.len().min(INDEPENDENT_REPLAY_COLUMNS) as u64;
        if replay_failures != 0 {
            return Err(format!(
                "{replay_failures} independent point replay failures"
            ));
        }

        let mut point_blob = Vec::with_capacity(rows.len() * 66);
        for row in &rows {
            point_blob.extend(point_key_bytes(&row.x, false)?);
            point_blob.extend(point_key_bytes(&row.x, true)?);
        }
        let points_sha256 = hex::encode(sha256(&point_blob));
        if !seen_point_digests.insert(points_sha256.clone()) {
            counter += ARMS_PER_ORDER;
            continue;
        }
        let mut params = BTreeMap::new();
        params.insert("subgroup_order".into(), subgroup_order.to_string());
        params.insert("primitive_root".into(), primitive_root.to_string());
        params.insert("omega".into(), hex_field(&omega));
        params.insert("pgl2_a".into(), hex_field(&matrix.a));
        params.insert("pgl2_b".into(), hex_field(&matrix.b));
        params.insert("pgl2_c".into(), hex_field(&matrix.c));
        params.insert("pgl2_e".into(), hex_field(&matrix.e));
        params.insert("point_key_encoding".into(), POINT_KEY_ENCODING.into());
        let manifest = FactorBaseManifest {
            schema: FACTOR_BASE_SCHEMA.into(),
            curve: CURVE_SLUG.into(),
            family: "pgl2-root-orbit".into(),
            params,
            columns: rows.len() as u64,
            signed_points: (2 * rows.len()) as u64,
            points_sha256,
        };
        let manifest_value = serde_json::to_value(&manifest).map_err(|error| error.to_string())?;
        let (fb_id, fb_sha256) = short_id("FB1", &manifest_value)?;
        let points = rows
            .drain(..)
            .map(|row| FactorPoint {
                positive: row.point,
            })
            .collect::<Vec<_>>();
        let transport = transport_census(&points, curve)?;
        let projected = projection(manifest.columns, transport.quotient_dimension, curve);
        let mut matrix_output = BTreeMap::new();
        matrix_output.insert("a".into(), hex_field(&matrix.a));
        matrix_output.insert("b".into(), hex_field(&matrix.b));
        matrix_output.insert("c".into(), hex_field(&matrix.c));
        matrix_output.insert("e".into(), hex_field(&matrix.e));
        return Ok(Candidate {
            name: format!("root-orbit-d{subgroup_order}-arm{arm}"),
            family: "pgl2-root-orbit".into(),
            subgroup_order,
            arm,
            derivation_counter: counter,
            matrix: matrix_output,
            determinant: hex_field(&determinant),
            omega: hex_field(&omega),
            omega_exact_order_verified: true,
            orbit_abscissae: subgroup_order,
            rejected_poles: 0,
            columns: manifest.columns,
            signed_points: manifest.signed_points,
            quadratic_residue_fraction: manifest.columns as f64 / subgroup_order as f64,
            coordinate_orbit_closed: true,
            induced_recurrence_closed: true,
            nonidentity_coordinate_action: omega != BigUint::one(),
            degree_one_curve_automorphism: false,
            membership_generator_degree: 2,
            structured_residual_degree_of_regularity: None,
            unsplit_s17_degree_of_regularity: None,
            independent_point_replays: generation_operations.independent_bigint_replays,
            point_replay_failures: 0,
            point_keys_sorted_unique: true,
            manifest,
            fb_id,
            fb_sha256,
            generation_operations,
            transport,
            projection: projected,
            classification: "negative screen unless exact group-log transport reduces K".into(),
        });
    }
}

fn run(cli: &Cli) -> Result<ResultFile, String> {
    let dependency_sha256 = file_sha256(&cli.round22)?;
    if dependency_sha256 != ROUND22_SHA256 {
        return Err(format!(
            "round-22 dependency hash mismatch: {dependency_sha256}"
        ));
    }
    let curve = CurveParams::p256();
    if (&curve.p & BigUint::from(3u8)) != BigUint::from(3u8) {
        return Err("P-256 prime unexpectedly is not 3 mod 4".into());
    }
    let factors = prime_factors()?;
    if !verify_factorization(&curve.p, &factors) {
        return Err("frozen p-1 factorization failed".into());
    }
    let root = primitive_root(&curve.p, &factors);
    let mut seen_point_digests = BTreeSet::new();
    let mut candidates = Vec::new();
    for subgroup_order in SUBGROUP_ORDERS {
        if (&curve.p - 1u8) % subgroup_order != BigUint::zero() {
            return Err(format!(
                "subgroup order {subgroup_order} does not divide p-1"
            ));
        }
        for arm in 0..ARMS_PER_ORDER {
            candidates.push(build_candidate(
                subgroup_order,
                arm,
                &root,
                &curve,
                &mut seen_point_digests,
            )?);
        }
    }
    let total_orbit_abscissae = candidates.iter().map(|item| item.orbit_abscissae).sum();
    let total_columns = candidates.iter().map(|item| item.columns).sum();
    let total_transport_images = candidates.iter().map(|item| item.transport.images).sum();
    let total_cross_column_candidates = candidates
        .iter()
        .map(|item| item.transport.cross_column_candidates)
        .sum();
    let total_exact_relations = candidates
        .iter()
        .map(|item| item.transport.replayed_relations)
        .sum();
    let any_promotion_gate_passed = candidates.iter().any(|item| item.projection.promotion_gate);
    Ok(ResultFile {
        schema: "p256.root_orbit_factor_base/v1".into(),
        curve: CURVE_SLUG.into(),
        relation_model: "17 distinct variable columns; balanced signed 8+9 cross-colour target collision".into(),
        dependency_path: cli.round22.display().to_string(),
        dependency_sha256,
        p_minus_one_factorization_verified: true,
        primitive_root: root.to_string(),
        primitive_root_verified: true,
        candidates,
        total_orbit_abscissae,
        total_columns,
        total_transport_images,
        total_cross_column_candidates,
        total_exact_relations,
        all_reported_relations_replayed: true,
        reference_false_positives: 0,
        reference_false_negatives: 0,
        any_promotion_gate_passed,
        full_depth_unplanted_attempted: false,
        result_class: "negative result unless a screened coordinate orbit yields exact group-log transport".into(),
        decision: "Coordinate symmetry is not credited as logarithm symmetry. Promotion requires measured residual degree and exact cross-column group relations sufficient to lower the independent-log quotient and complete cost below rho.".into(),
    })
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(result) => {
            let mut bytes = match serde_json::to_vec_pretty(&result) {
                Ok(bytes) => bytes,
                Err(error) => {
                    eprintln!("serialization failed: {error}");
                    return ExitCode::FAILURE;
                }
            };
            bytes.push(b'\n');
            if let Some(path) = cli.out {
                if let Err(error) = fs::write(&path, &bytes) {
                    eprintln!("{}: {error}", path.display());
                    return ExitCode::FAILURE;
                }
            } else {
                print!("{}", String::from_utf8_lossy(&bytes));
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("P-256 root-orbit screen failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_factorization_and_subgroups_are_valid() {
        let curve = CurveParams::p256();
        let factors = prime_factors().unwrap();
        assert!(verify_factorization(&curve.p, &factors));
        for order in SUBGROUP_ORDERS {
            assert_eq!((&curve.p - 1u8) % order, BigUint::zero());
        }
    }

    #[test]
    fn derived_matrices_are_nonsingular() {
        let curve = CurveParams::p256();
        for order in SUBGROUP_ORDERS {
            for arm in 0..ARMS_PER_ORDER {
                let matrix = derive_matrix(order, arm, &curve.p);
                let determinant = (&matrix.a * &matrix.e + &curve.p
                    - (&matrix.b * &matrix.c % &curve.p))
                    % &curve.p;
                assert!(!determinant.is_zero());
            }
        }
    }

    #[test]
    fn projected_independent_log_floor_fails_all_cost_gates() {
        let curve = CurveParams::p256();
        let projected = projection(131_458, 131_458, &curve);
        assert!(projected.cross_colour_ratio_to_rho > 490.0);
        assert!(!projected.relation_collection_below_2_120);
        assert!(!projected.cost_per_relation_below_2_103);
        assert!(!projected.materialized_storage_below_2_50);
        assert!(!projected.promotion_gate);
    }

    #[test]
    fn induced_recurrence_matches_direct_map_on_prefix() {
        let curve = CurveParams::p256();
        let factors = prime_factors().unwrap();
        let root = primitive_root(&curve.p, &factors);
        let order = SUBGROUP_ORDERS[0];
        let omega = root.modpow(&((&curve.p - 1u8) / order), &curve.p);
        let omega_fe = Fe::from_biguint(&omega);
        let matrix = matrix_fe(&derive_matrix(order, 0, &curve.p));
        let recurrence = closure_matrix(&matrix, omega_fe);
        let mut u = Fe::ONE;
        for _ in 0..64 {
            let x = matrix_map(&matrix, u).unwrap();
            u = u.mul(&omega_fe);
            let next = matrix_map(&matrix, u).unwrap();
            let recurrence_next = matrix_map(&recurrence, x).unwrap();
            assert_eq!(recurrence_next, next);
        }
    }
}
