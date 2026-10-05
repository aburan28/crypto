//! Transfer the indexed final-S3 image gate to the committed P-256 factor base.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::Write as _;
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{self, CURVE_SLUG};
use crypto_lib::ecc::p256_field::P256FieldElement as Fe;
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
const COLUMNS: u64 = 131_458;
const SIGNED_POINTS: u64 = 262_916;
const ROUND13_SHA256: &str = "2008dcb659d3120157480b6096a4873d1f9c23ee30123f4dd3a32d450f53ba1d";
const ROUND14_SHA256: &str = "6eefb6b27768023af5850cd785d75ef3729484ac85e5defef9d1d596834ebccc";
const ROUND15_SHA256: &str = "4220dafa066613630338886a916218602cd61ca914bd66ac7994cde163a53ad5";
const ATOM_COLUMNS: usize = 16;
const COLUMN_INDEX_BITS: usize = 18;
const PACKED_ATOM_BYTES: usize = ATOM_COLUMNS * COLUMN_INDEX_BITS / 8;
const TARGET_PREIMAGE: &str = concat!(
    "icv1-fp256-t89188191154553853111372247798585809583-f188c491",
    "/s3-image-transfer-round8/target/0"
);

#[derive(Parser)]
#[command(about = "Validate the indexed final-S3 gate on the selected P-256 factor base")]
struct Cli {
    /// Compact deterministic JSON output; stdout when omitted.
    #[arg(long)]
    out: Option<PathBuf>,
    /// Run the round-14 P-256 sixteen-leaf width transfer after hash-checking round 13.
    #[arg(long)]
    width_round13: Option<PathBuf>,
    /// Run the round-15 target-indexed compression after hash-checking round 14.
    #[arg(long)]
    compress_round14: Option<PathBuf>,
    /// Run the round-16 packed-atom compression after hash-checking round 15.
    #[arg(long)]
    pack_round15: Option<PathBuf>,
}

#[derive(Clone)]
struct Column {
    x: Fe,
    low: Point,
}

#[derive(Clone)]
struct SignedRow {
    point: Point,
    col: u64,
}

#[derive(Clone, Copy, Default, Serialize)]
struct MultiplicationCounts {
    coefficients: u64,
    discriminants: u64,
    square_roots: u64,
    inversions: u64,
    root_construction: u64,
    total: u64,
}

impl MultiplicationCounts {
    fn finish(&mut self) {
        self.total = self.coefficients
            + self.discriminants
            + self.square_roots
            + self.inversions
            + self.root_construction;
    }

    fn combined(self, other: Self) -> Self {
        let mut combined = Self {
            coefficients: self.coefficients + other.coefficients,
            discriminants: self.discriminants + other.discriminants,
            square_roots: self.square_roots + other.square_roots,
            inversions: self.inversions + other.inversions,
            root_construction: self.root_construction + other.root_construction,
            total: 0,
        };
        combined.finish();
        combined
    }
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
}

#[derive(Serialize)]
struct TargetResult {
    kind: String,
    scalar: Option<String>,
    target: [String; 2],
    rejected_ordered_column_pairs: u64,
    quadratic_solves: u64,
    quadratic_roots_returned: u64,
    algebra_index_lookups: u64,
    linear_degeneracies: u64,
    universal_degeneracies: u64,
    algebra_pairs: usize,
    reference_pairs: usize,
    false_negatives: usize,
    false_positives: usize,
    multiplication_counts: MultiplicationCounts,
    reference_group_subtractions: u64,
    reference_index_lookups: u64,
    verification_group_additions: u64,
    all_emitted_pairs_verified: bool,
    exact: bool,
    pair_sha256: String,
    pairs: Vec<[u64; 2]>,
}

#[derive(Serialize)]
struct ExperimentResult {
    schema: String,
    curve: String,
    field_prime: String,
    curve_a: String,
    curve_b: String,
    factor_base: FactorBaseReceipt,
    target_preimage: String,
    targets: Vec<TargetResult>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct P256ImageState {
    affine: BTreeMap<[u8; 32], Fe>,
    identity: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ReferenceImageState {
    affine: BTreeSet<[u8; 32]>,
    identity: bool,
}

#[derive(Serialize)]
struct WidthSampleResult {
    kind: String,
    columns: Vec<u64>,
    x_coordinates: Vec<String>,
    two_leaf_widths: Vec<usize>,
    four_leaf_widths: Vec<usize>,
    eight_leaf_widths: Vec<usize>,
    sixteen_leaf_width: usize,
    identity_images_by_level: [u64; 4],
    quadratic_solves: u64,
    quadratic_roots_returned: u64,
    linear_degeneracies: u64,
    universal_degeneracies: u64,
    multiplication_counts: MultiplicationCounts,
    reference_group_additions: u64,
    exact: bool,
    final_image_sha256: String,
}

#[derive(Serialize)]
struct WidthStorageModel {
    generic_sixteen_leaf_x_entries: u64,
    bytes_per_x: u64,
    generic_sixteen_leaf_raw_bytes: u64,
    factor_base_columns: u64,
    unordered_factor_base_pairs: u64,
    two_root_pair_image_entry_upper_bound: u64,
    two_root_pair_image_raw_byte_upper_bound: u64,
    identity_aware_pair_image_affine_entries: u64,
    identity_aware_pair_image_raw_bytes: u64,
}

#[derive(Serialize)]
struct WidthExperimentResult {
    schema: String,
    curve: String,
    field_prime: String,
    curve_a: String,
    curve_b: String,
    round13_sha256: String,
    factor_base: FactorBaseReceipt,
    local_maximum_degree: u32,
    generic_widths: [u64; 4],
    samples: Vec<WidthSampleResult>,
    storage_model: WidthStorageModel,
}

#[derive(Serialize)]
struct CompressedTargetResult {
    kind: String,
    scalar: Option<String>,
    target: [String; 2],
    reference_positive: bool,
    candidate_positive: bool,
    algebra_hits: u64,
    quadratic_solves_warm: u64,
    quadratic_solves_cold: u64,
    quadratic_roots_returned: u64,
    algebra_index_lookups: u64,
    linear_degeneracies: u64,
    universal_degeneracies: u64,
    multiplication_counts_warm: MultiplicationCounts,
    multiplication_counts_cold: MultiplicationCounts,
    false_negative: bool,
    false_positive: bool,
    exact: bool,
    hit_sha256: String,
}

#[derive(Serialize)]
struct CompressionSampleResult {
    kind: String,
    columns: Vec<u64>,
    x_coordinates: Vec<String>,
    left_eight_leaf_entries: usize,
    right_eight_leaf_entries: usize,
    retained_entries: usize,
    materialized_sixteen_leaf_entries: usize,
    retained_raw_bytes: usize,
    materialized_raw_bytes: usize,
    retained_to_materialized_ratio: f64,
    build_quadratic_solves: u64,
    full_materialization_quadratic_solves: u64,
    intermediate_reference_group_additions: u64,
    full_reference_group_additions: u64,
    build_linear_degeneracies: u64,
    build_universal_degeneracies: u64,
    build_multiplication_counts: MultiplicationCounts,
    left_image_sha256: String,
    right_image_sha256: String,
    targets: Vec<CompressedTargetResult>,
}

#[derive(Serialize)]
struct CompressionExperimentResult {
    schema: String,
    curve: String,
    field_prime: String,
    curve_a: String,
    curve_b: String,
    round14_sha256: String,
    factor_base: FactorBaseReceipt,
    local_maximum_degree: u32,
    samples: Vec<CompressionSampleResult>,
}

#[derive(Serialize)]
struct PackedAtomSampleResult {
    kind: String,
    packed_atom_hex: String,
    packed_atom_sha256: String,
    bits_per_column: usize,
    packet_bytes: usize,
    round15_retained_raw_bytes: usize,
    materialized_sixteen_leaf_raw_bytes: usize,
    additional_persistent_reduction: f64,
    full_materialized_reduction: f64,
    columns: Vec<u64>,
    x_coordinates: Vec<String>,
    decoded_exact: bool,
    reencoded_exact: bool,
    dependency_sample_exact: bool,
    left_image_sha256: String,
    right_image_sha256: String,
    transient_reconstruction_raw_bytes: usize,
    build_quadratic_solves: u64,
    one_cold_query_quadratic_solves: u64,
    two_query_one_rebuild_quadratic_solves: u64,
    intermediate_reference_group_additions: u64,
    full_reference_group_additions: u64,
    build_linear_degeneracies: u64,
    build_universal_degeneracies: u64,
    build_multiplication_counts: MultiplicationCounts,
    targets: Vec<CompressedTargetResult>,
}

#[derive(Serialize)]
struct PackedAtomExperimentResult {
    schema: String,
    curve: String,
    field_prime: String,
    curve_a: String,
    curve_b: String,
    round15_sha256: String,
    factor_base: FactorBaseReceipt,
    local_maximum_degree: u32,
    factor_base_column_bits: usize,
    packed_atom_bytes: usize,
    samples: Vec<PackedAtomSampleResult>,
}

#[derive(Clone, Copy)]
struct QuadraticRoots {
    roots: [Fe; 2],
    len: usize,
    linear: bool,
    universal: bool,
}

impl QuadraticRoots {
    fn empty() -> Self {
        Self {
            roots: [Fe::ZERO; 2],
            len: 0,
            linear: false,
            universal: false,
        }
    }
}

#[derive(Default)]
struct CountedField {
    multiplications: u64,
}

impl CountedField {
    fn mul(&mut self, left: Fe, right: Fe) -> Fe {
        self.multiplications += 1;
        left.mul(&right)
    }

    fn square(&mut self, value: Fe) -> Fe {
        self.multiplications += 1;
        value.sqr()
    }

    /// Right-to-left binary exponentiation, with every executed field
    /// multiplication charged.  The final unused base square is deliberately
    /// charged, matching the native loop that was run.
    fn pow(&mut self, mut base: Fe, exponent: &BigUint) -> Fe {
        let mut result = Fe::ONE;
        for bit in 0..exponent.bits() {
            if exponent.bit(bit) {
                result = self.mul(result, base);
            }
            base = self.square(base);
        }
        result
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

fn point(x: &BigUint, y: &BigUint, curve: &CurveParams) -> Point {
    Point::Affine {
        x: curve.fe(x.clone()),
        y: curve.fe(y.clone()),
    }
}

fn affine_coordinates(value: &Point) -> Option<(&BigUint, &BigUint)> {
    match value {
        Point::Infinity => None,
        Point::Affine { x, y } => Some((&x.value, &y.value)),
    }
}

fn affine_key(value: &Point) -> Option<([u8; 32], [u8; 32])> {
    let (x, y) = affine_coordinates(value)?;
    Some((big_key(x), big_key(y)))
}

fn big_key(value: &BigUint) -> [u8; 32] {
    let bytes = value.to_bytes_be();
    assert!(bytes.len() <= 32, "P-256 coordinate exceeded 32 bytes");
    let mut raw = [0; 32];
    raw[32 - bytes.len()..].copy_from_slice(&bytes);
    raw
}

fn fe_key(value: Fe) -> [u8; 32] {
    value.to_bytes_be()
}

fn negate(value: &Point, curve: &CurveParams) -> Point {
    match value {
        Point::Infinity => Point::Infinity,
        Point::Affine { x, y } => Point::Affine {
            x: x.clone(),
            y: curve.fe(if y.value.is_zero() {
                BigUint::zero()
            } else {
                &curve.p - &y.value
            }),
        },
    }
}

fn coefficients(u: Fe, t: Fe, a: Fe, b: Fe, count: &mut CountedField) -> (Fe, Fe, Fe) {
    let two_b = b.add(&b);
    let t2 = count.square(t);
    let u2 = count.square(u);
    let a2 = count.square(a);
    let au = count.mul(a, u);
    let qa = count.square(t.sub(&u));

    let mut qb_inner = count.mul(t2, u);
    qb_inner = qb_inner.add(&count.mul(t, u2));
    qb_inner = qb_inner.add(&count.mul(t, a));
    qb_inner = qb_inner.add(&au);
    qb_inner = qb_inner.add(&two_b);
    let qb = qb_inner.add(&qb_inner).neg();

    let mut qc = count.mul(t2, u2);
    let t_inner = count.mul(t, au.add(&two_b));
    qc = qc.sub(&t_inner.add(&t_inner));
    qc = qc.add(&a2);
    let bu = count.mul(b, u);
    qc = qc.sub(&bu.add(&bu).add(&bu.add(&bu)));
    (qa, qb, qc)
}

fn solve_quadratic(
    qa: Fe,
    qb: Fe,
    qc: Fe,
    sqrt_exponent: &BigUint,
    inverse_exponent: &BigUint,
    counts: &mut MultiplicationCounts,
) -> QuadraticRoots {
    if qa == Fe::ZERO {
        if qb == Fe::ZERO {
            return QuadraticRoots {
                universal: qc == Fe::ZERO,
                ..QuadraticRoots::empty()
            };
        }
        let mut inverse_field = CountedField::default();
        let inverse = inverse_field.pow(qb, inverse_exponent);
        counts.inversions += inverse_field.multiplications;
        let mut root_field = CountedField::default();
        let root = root_field.mul(qc.neg(), inverse);
        counts.root_construction += root_field.multiplications;
        return QuadraticRoots {
            roots: [root, Fe::ZERO],
            len: 1,
            linear: true,
            universal: false,
        };
    }

    let mut discriminant_field = CountedField::default();
    let qb2 = discriminant_field.square(qb);
    let ac = discriminant_field.mul(qa, qc);
    let four_ac = ac.add(&ac).add(&ac.add(&ac));
    let discriminant = qb2.sub(&four_ac);
    counts.discriminants += discriminant_field.multiplications;

    let mut sqrt_field = CountedField::default();
    let sqrt = sqrt_field.pow(discriminant, sqrt_exponent);
    let sqrt_check = sqrt_field.square(sqrt);
    counts.square_roots += sqrt_field.multiplications;
    if sqrt_check != discriminant {
        return QuadraticRoots::empty();
    }

    let mut inverse_field = CountedField::default();
    let inverse = inverse_field.pow(qa.add(&qa), inverse_exponent);
    counts.inversions += inverse_field.multiplications;

    let mut root_field = CountedField::default();
    let neg_qb = qb.neg();
    let first = root_field.mul(neg_qb.add(&sqrt), inverse);
    let second = root_field.mul(neg_qb.sub(&sqrt), inverse);
    counts.root_construction += root_field.multiplications;
    if first == second {
        QuadraticRoots {
            roots: [first, Fe::ZERO],
            len: 1,
            linear: false,
            universal: false,
        }
    } else {
        QuadraticRoots {
            roots: [first, second],
            len: 2,
            linear: false,
            universal: false,
        }
    }
}

fn algebra_pairs(
    target: &Point,
    columns: &[Column],
    x_index: &HashMap<[u8; 32], u64>,
    curve: &CurveParams,
) -> Result<(BTreeSet<(u64, u64)>, u64, u64, u64, MultiplicationCounts), String> {
    let (target_x, _) = affine_coordinates(target).ok_or("target is infinity")?;
    let t = Fe::from_biguint(target_x);
    let a = Fe::from_biguint(&curve.a);
    let b = Fe::from_biguint(&curve.b);
    let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
    let inverse_exponent = &curve.p - BigUint::from(2u8);
    let mut pairs = BTreeSet::new();
    let mut roots_returned = 0u64;
    let mut linear = 0u64;
    let mut universal = 0u64;
    let mut counts = MultiplicationCounts::default();
    for (left_col, column) in columns.iter().enumerate() {
        let mut coefficient_field = CountedField::default();
        let (qa, qb, qc) = coefficients(column.x, t, a, b, &mut coefficient_field);
        counts.coefficients += coefficient_field.multiplications;
        let roots = solve_quadratic(qa, qb, qc, &sqrt_exponent, &inverse_exponent, &mut counts);
        linear += u64::from(roots.linear);
        universal += u64::from(roots.universal);
        if roots.universal {
            return Err(format!(
                "universal quadratic at left column {left_col} has no bounded image"
            ));
        }
        roots_returned += roots.len as u64;
        for root in roots.roots[..roots.len].iter().copied() {
            if let Some(&right_col) = x_index.get(&fe_key(root)) {
                pairs.insert((left_col as u64, right_col));
            }
        }
    }
    counts.finish();
    Ok((pairs, roots_returned, linear, universal, counts))
}

fn reference_pairs(
    target: &Point,
    signed_rows: &[SignedRow],
    signed_index: &HashMap<([u8; 32], [u8; 32]), u64>,
    curve: &CurveParams,
) -> BTreeSet<(u64, u64)> {
    let a = curve.a_fe();
    let mut pairs = BTreeSet::new();
    for row in signed_rows {
        let right = target.add_vartime(&negate(&row.point, curve), &a);
        if let Some(key) = affine_key(&right) {
            if let Some(&right_col) = signed_index.get(&key) {
                pairs.insert((row.col, right_col));
            }
        }
    }
    pairs
}

fn verify_pairs(
    pairs: &BTreeSet<(u64, u64)>,
    target: &Point,
    columns: &[Column],
    curve: &CurveParams,
) -> Result<u64, String> {
    let a = curve.a_fe();
    let neg_target = negate(target, curve);
    let mut additions = 0u64;
    for &(left, right) in pairs {
        let left_point = &columns[left as usize].low;
        let right_point = &columns[right as usize].low;
        let mut verified = false;
        for left_sign in [false, true] {
            for right_sign in [false, true] {
                let lp = if left_sign {
                    negate(left_point, curve)
                } else {
                    left_point.clone()
                };
                let rp = if right_sign {
                    negate(right_point, curve)
                } else {
                    right_point.clone()
                };
                let sum = lp.add_vartime(&rp, &a);
                additions += 1;
                verified |= sum == *target || sum == neg_target;
            }
        }
        if !verified {
            return Err(format!(
                "emitted column pair ({left}, {right}) has no signed sum +/-Q"
            ));
        }
    }
    Ok(additions)
}

fn pair_digest(pairs: &BTreeSet<(u64, u64)>) -> String {
    let mut text = String::new();
    for (left, right) in pairs {
        writeln!(text, "{left},{right}").expect("writing to String cannot fail");
    }
    hex::encode(sha256(text.as_bytes()))
}

fn run_target(
    kind: &str,
    scalar: Option<&BigUint>,
    target: &Point,
    columns: &[Column],
    signed_rows: &[SignedRow],
    x_index: &HashMap<[u8; 32], u64>,
    signed_index: &HashMap<([u8; 32], [u8; 32]), u64>,
    curve: &CurveParams,
) -> Result<TargetResult, String> {
    if matches!(target, Point::Infinity) || !curve.is_on_curve(target) {
        return Err(format!("{kind} target is infinity or off curve"));
    }
    let (algebra, roots_returned, linear, universal, multiplication_counts) =
        algebra_pairs(target, columns, x_index, curve)?;
    let reference = reference_pairs(target, signed_rows, signed_index, curve);
    let false_negatives = reference.difference(&algebra).count();
    let false_positives = algebra.difference(&reference).count();
    let verification_group_additions = verify_pairs(&algebra, target, columns, curve)?;
    let (x, y) = affine_coordinates(target).expect("target was checked affine");
    Ok(TargetResult {
        kind: kind.into(),
        scalar: scalar.map(lower_hex),
        target: [lower_hex(x), lower_hex(y)],
        rejected_ordered_column_pairs: COLUMNS * COLUMNS,
        quadratic_solves: COLUMNS,
        quadratic_roots_returned: roots_returned,
        algebra_index_lookups: roots_returned,
        linear_degeneracies: linear,
        universal_degeneracies: universal,
        algebra_pairs: algebra.len(),
        reference_pairs: reference.len(),
        false_negatives,
        false_positives,
        multiplication_counts,
        reference_group_subtractions: SIGNED_POINTS,
        reference_index_lookups: SIGNED_POINTS,
        verification_group_additions,
        all_emitted_pairs_verified: true,
        exact: false_negatives == 0 && false_positives == 0,
        pair_sha256: pair_digest(&algebra),
        pairs: algebra.into_iter().map(|(l, r)| [l, r]).collect(),
    })
}

fn build_indexes(
    curve: &CurveParams,
) -> Result<
    (
        FactorBaseReceipt,
        Vec<Column>,
        Vec<SignedRow>,
        HashMap<[u8; 32], u64>,
        HashMap<([u8; 32], [u8; 32]), u64>,
    ),
    String,
> {
    let built = p256_dickson_factor_base::build(SPEC)?;
    let dump = &built.dump;
    if dump.factor_base.fb_id != FB_ID
        || dump.factor_base.fb_sha256 != FB_SHA256
        || dump.factor_base.points_sha256 != POINTS_SHA256
        || dump.factor_base.columns != COLUMNS
        || dump.factor_base.signed_points != SIGNED_POINTS
        || dump.factor_base.params.get("terminal").map(String::as_str) != Some(TERMINAL)
    {
        return Err("rebuilt factor base does not match the frozen FB1 identity".into());
    }
    p256_dickson_factor_base::verify(dump)?;

    let mut columns: Vec<Option<Column>> = vec![None; COLUMNS as usize];
    let mut signed_rows = Vec::with_capacity(SIGNED_POINTS as usize);
    let mut x_index = HashMap::with_capacity(COLUMNS as usize);
    let mut signed_index = HashMap::with_capacity(SIGNED_POINTS as usize);
    for row in &dump.points {
        let x = parse_big(&row.x)?;
        let y = parse_big(&row.y)?;
        let affine = point(&x, &y, curve);
        if !curve.is_on_curve(&affine) {
            return Err(format!(
                "factor-base row in column {} is off curve",
                row.col
            ));
        }
        let key = affine_key(&affine).expect("factor-base rows are affine");
        if signed_index.insert(key, row.col).is_some() {
            return Err("duplicate signed factor-base point".into());
        }
        signed_rows.push(SignedRow {
            point: affine.clone(),
            col: row.col,
        });
        if row.coef == "1" {
            let x_fe = Fe::from_biguint(&x);
            if x_index.insert(fe_key(x_fe), row.col).is_some() {
                return Err("duplicate factor-base abscissa".into());
            }
            if columns[row.col as usize]
                .replace(Column {
                    x: x_fe,
                    low: affine,
                })
                .is_some()
            {
                return Err(format!("duplicate low-y row for column {}", row.col));
            }
        }
    }
    let columns: Vec<Column> = columns
        .into_iter()
        .enumerate()
        .map(|(index, column)| column.ok_or_else(|| format!("missing column {index}")))
        .collect::<Result<_, _>>()?;
    if columns.len() as u64 != COLUMNS
        || signed_rows.len() as u64 != SIGNED_POINTS
        || x_index.len() as u64 != COLUMNS
        || signed_index.len() as u64 != SIGNED_POINTS
    {
        return Err("factor-base index cardinality mismatch".into());
    }
    Ok((
        FactorBaseReceipt {
            spec: SPEC.into(),
            fb_id: FB_ID.into(),
            fb_sha256: FB_SHA256.into(),
            points_sha256: POINTS_SHA256.into(),
            terminal: TERMINAL.into(),
            columns: COLUMNS,
            signed_points: SIGNED_POINTS,
            complete_rebuild_verified: true,
        },
        columns,
        signed_rows,
        x_index,
        signed_index,
    ))
}

fn full_point_key(value: &Point) -> [u8; 65] {
    let mut key = [0u8; 65];
    if let Some((x, y)) = affine_coordinates(value) {
        key[0] = 1;
        key[1..33].copy_from_slice(&big_key(x));
        key[33..].copy_from_slice(&big_key(y));
    }
    key
}

fn reference_image(
    points: &[Point],
    curve: &CurveParams,
) -> Result<(ReferenceImageState, u64), String> {
    let a = curve.a_fe();
    let mut sums = BTreeMap::from([(full_point_key(&Point::Infinity), Point::Infinity)]);
    let mut additions = 0u64;
    for point in points {
        let negated = negate(point, curve);
        if negated == *point {
            return Err("factor-base leaf has only one signed lift".into());
        }
        let mut next = BTreeMap::new();
        for sum in sums.values() {
            for signed in [point, &negated] {
                let value = sum.add_vartime(signed, &a);
                additions += 1;
                next.insert(full_point_key(&value), value);
            }
        }
        sums = next;
    }
    let mut affine = BTreeSet::new();
    let mut identity = false;
    for point in sums.values() {
        match affine_coordinates(point) {
            Some((x, _)) => {
                affine.insert(big_key(x));
            }
            None => identity = true,
        }
    }
    Ok((ReferenceImageState { affine, identity }, additions))
}

fn image_matches_reference(candidate: &P256ImageState, reference: &ReferenceImageState) -> bool {
    candidate.identity == reference.identity
        && candidate.affine.len() == reference.affine.len()
        && candidate
            .affine
            .keys()
            .zip(reference.affine.iter())
            .all(|(left, right)| left == right)
}

#[allow(clippy::too_many_arguments)]
fn compose_p256_images(
    left: &P256ImageState,
    right: &P256ImageState,
    a: Fe,
    b: Fe,
    sqrt_exponent: &BigUint,
    inverse_exponent: &BigUint,
    solves: &mut u64,
    roots_returned: &mut u64,
    linear: &mut u64,
    universal: &mut u64,
    counts: &mut MultiplicationCounts,
) -> Result<P256ImageState, String> {
    let mut affine = BTreeMap::new();
    if left.identity {
        affine.extend(right.affine.iter().map(|(key, value)| (*key, *value)));
    }
    if right.identity {
        affine.extend(left.affine.iter().map(|(key, value)| (*key, *value)));
    }
    for &u in left.affine.values() {
        for &v in right.affine.values() {
            *solves += 1;
            let mut coefficient_field = CountedField::default();
            let (qa, qb, qc) = coefficients(u, v, a, b, &mut coefficient_field);
            counts.coefficients += coefficient_field.multiplications;
            let roots = solve_quadratic(qa, qb, qc, sqrt_exponent, inverse_exponent, counts);
            *linear += u64::from(roots.linear);
            *universal += u64::from(roots.universal);
            if roots.universal {
                return Err("universal quadratic has no bounded image".into());
            }
            *roots_returned += roots.len as u64;
            for root in roots.roots[..roots.len].iter().copied() {
                affine.insert(fe_key(root), root);
            }
        }
    }
    let shared_affine = left.affine.keys().any(|key| right.affine.contains_key(key));
    Ok(P256ImageState {
        affine,
        identity: (left.identity && right.identity) || shared_affine,
    })
}

fn image_digest(image: &P256ImageState) -> String {
    let mut bytes = Vec::with_capacity(1 + 32 * image.affine.len());
    bytes.push(u8::from(image.identity));
    for key in image.affine.keys() {
        bytes.extend_from_slice(key);
    }
    hex::encode(sha256(&bytes))
}

fn verify_width_node(
    kind: &str,
    level: usize,
    node: usize,
    candidate: &P256ImageState,
    points: &[Point],
    curve: &CurveParams,
    reference_group_additions: &mut u64,
) -> Result<(), String> {
    let (reference, additions) = reference_image(points, curve)?;
    *reference_group_additions += additions;
    if !image_matches_reference(candidate, &reference) {
        return Err(format!(
            "{kind} {level}-leaf node {node} differs from exact signed group image"
        ));
    }
    Ok(())
}

fn run_width_sample(
    kind: &str,
    column_indices: Vec<u64>,
    columns: &[Column],
    curve: &CurveParams,
) -> Result<WidthSampleResult, String> {
    if column_indices.len() != 16
        || column_indices
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .len()
            != 16
    {
        return Err(format!("{kind} does not contain 16 distinct columns"));
    }
    let selected: Vec<&Column> = column_indices
        .iter()
        .map(|&index| {
            columns
                .get(index as usize)
                .ok_or_else(|| format!("{kind} column {index} is out of range"))
        })
        .collect::<Result<_, _>>()?;
    let points: Vec<Point> = selected.iter().map(|column| column.low.clone()).collect();
    let x_coordinates: Vec<String> = selected
        .iter()
        .map(|column| lower_hex(&BigUint::from_bytes_be(&column.x.to_bytes_be())))
        .collect();
    let mut images: Vec<P256ImageState> = selected
        .iter()
        .map(|column| P256ImageState {
            affine: BTreeMap::from([(fe_key(column.x), column.x)]),
            identity: false,
        })
        .collect();
    let a = Fe::from_biguint(&curve.a);
    let b = Fe::from_biguint(&curve.b);
    let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
    let inverse_exponent = &curve.p - BigUint::from(2u8);
    let mut solves = 0u64;
    let mut roots_returned = 0u64;
    let mut linear = 0u64;
    let mut universal = 0u64;
    let mut counts = MultiplicationCounts::default();
    let mut reference_group_additions = 0u64;
    let mut widths_by_level: Vec<Vec<usize>> = Vec::new();
    let mut identity_images_by_level = [0u64; 4];

    for (level_index, block_size) in [2usize, 4, 8, 16].into_iter().enumerate() {
        let mut next = Vec::with_capacity(images.len() / 2);
        for (node, pair) in images.chunks_exact(2).enumerate() {
            let image = compose_p256_images(
                &pair[0],
                &pair[1],
                a,
                b,
                &sqrt_exponent,
                &inverse_exponent,
                &mut solves,
                &mut roots_returned,
                &mut linear,
                &mut universal,
                &mut counts,
            )?;
            let start = node * block_size;
            verify_width_node(
                kind,
                block_size,
                node,
                &image,
                &points[start..start + block_size],
                curve,
                &mut reference_group_additions,
            )?;
            identity_images_by_level[level_index] += u64::from(image.identity);
            next.push(image);
        }
        widths_by_level.push(next.iter().map(|image| image.affine.len()).collect());
        images = next;
    }
    if images.len() != 1 {
        return Err(format!("{kind} did not reduce to one sixteen-leaf image"));
    }
    counts.finish();
    let final_image = &images[0];
    Ok(WidthSampleResult {
        kind: kind.into(),
        columns: column_indices,
        x_coordinates,
        two_leaf_widths: widths_by_level[0].clone(),
        four_leaf_widths: widths_by_level[1].clone(),
        eight_leaf_widths: widths_by_level[2].clone(),
        sixteen_leaf_width: widths_by_level[3][0],
        identity_images_by_level,
        quadratic_solves: solves,
        quadratic_roots_returned: roots_returned,
        linear_degeneracies: linear,
        universal_degeneracies: universal,
        multiplication_counts: counts,
        reference_group_additions,
        exact: true,
        final_image_sha256: image_digest(final_image),
    })
}

fn hash_sample_indices(sample: usize) -> Vec<u64> {
    let mut indices = Vec::with_capacity(16);
    for leaf in 0..16 {
        let mut counter = 0u64;
        loop {
            let preimage = format!(
                "{CURVE_SLUG}/s17-image-width-round14/sample/{sample}/leaf/{leaf}/counter/{counter}"
            );
            let digest = sha256(preimage.as_bytes());
            let index = (BigUint::from_bytes_be(&digest) % BigUint::from(COLUMNS))
                .to_u64()
                .expect("reduced column index fits u64");
            if !indices.contains(&index) {
                indices.push(index);
                break;
            }
            counter += 1;
        }
    }
    indices
}

fn round14_sample_specs() -> Vec<(String, Vec<u64>)> {
    vec![
        ("prefix".to_string(), (0..16).collect::<Vec<u64>>()),
        ("hash-0".to_string(), hash_sample_indices(0)),
        ("hash-1".to_string(), hash_sample_indices(1)),
        ("hash-2".to_string(), hash_sample_indices(2)),
    ]
}

fn run_width(round13_path: &PathBuf, out: Option<PathBuf>) -> Result<(), String> {
    let bytes = std::fs::read(round13_path).map_err(|error| error.to_string())?;
    let round13_sha256 = hex::encode(sha256(&bytes));
    if round13_sha256 != ROUND13_SHA256 {
        return Err(format!(
            "round-13 SHA-256 mismatch: expected {ROUND13_SHA256}, got {round13_sha256}"
        ));
    }
    let curve = CurveParams::p256();
    let (factor_base, columns, _, _, _) = build_indexes(&curve)?;
    let sample_specs = round14_sample_specs();
    let mut samples = Vec::new();
    for (kind, indices) in sample_specs {
        samples.push(run_width_sample(&kind, indices, &columns, &curve)?);
    }
    let unordered_factor_base_pairs = COLUMNS * (COLUMNS + 1) / 2;
    let two_root_pair_image_entry_upper_bound = 2 * unordered_factor_base_pairs;
    let identity_aware_pair_image_affine_entries = COLUMNS * COLUMNS;
    let result = WidthExperimentResult {
        schema: "p256.s17_image_width_transfer/v1".into(),
        curve: CURVE_SLUG.into(),
        field_prime: lower_hex(&curve.p),
        curve_a: lower_hex(&curve.a),
        curve_b: lower_hex(&curve.b),
        round13_sha256,
        factor_base,
        local_maximum_degree: 2,
        generic_widths: [2, 8, 128, 32_768],
        samples,
        storage_model: WidthStorageModel {
            generic_sixteen_leaf_x_entries: 32_768,
            bytes_per_x: 32,
            generic_sixteen_leaf_raw_bytes: 32_768 * 32,
            factor_base_columns: COLUMNS,
            unordered_factor_base_pairs,
            two_root_pair_image_entry_upper_bound,
            two_root_pair_image_raw_byte_upper_bound: two_root_pair_image_entry_upper_bound * 32,
            identity_aware_pair_image_affine_entries,
            identity_aware_pair_image_raw_bytes: identity_aware_pair_image_affine_entries * 32,
        },
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn build_verified_eight_image(
    label: &str,
    xs: &[Fe],
    points: &[Point],
    curve: &CurveParams,
    a: Fe,
    b: Fe,
    sqrt_exponent: &BigUint,
    inverse_exponent: &BigUint,
    solves: &mut u64,
    roots_returned: &mut u64,
    linear: &mut u64,
    universal: &mut u64,
    counts: &mut MultiplicationCounts,
    reference_group_additions: &mut u64,
) -> Result<P256ImageState, String> {
    if xs.len() != 8 || points.len() != 8 {
        return Err(format!("{label} requires exactly eight leaves"));
    }
    let mut images: Vec<P256ImageState> = xs
        .iter()
        .copied()
        .map(|x| P256ImageState {
            affine: BTreeMap::from([(fe_key(x), x)]),
            identity: false,
        })
        .collect();
    for block_size in [2usize, 4, 8] {
        let mut next = Vec::with_capacity(images.len() / 2);
        for (node, pair) in images.chunks_exact(2).enumerate() {
            let image = compose_p256_images(
                &pair[0],
                &pair[1],
                a,
                b,
                sqrt_exponent,
                inverse_exponent,
                solves,
                roots_returned,
                linear,
                universal,
                counts,
            )?;
            let start = node * block_size;
            verify_width_node(
                label,
                block_size,
                node,
                &image,
                &points[start..start + block_size],
                curve,
                reference_group_additions,
            )?;
            next.push(image);
        }
        images = next;
    }
    if images.len() != 1 {
        return Err(format!("{label} did not reduce to one eight-leaf image"));
    }
    Ok(images.pop().expect("one image was checked"))
}

#[allow(clippy::too_many_arguments)]
fn query_compressed_target(
    kind: &str,
    scalar: Option<&BigUint>,
    target: &Point,
    reference_positive: bool,
    left: &P256ImageState,
    right: &P256ImageState,
    curve: &CurveParams,
    build_solves: u64,
    build_counts: MultiplicationCounts,
) -> Result<CompressedTargetResult, String> {
    if matches!(target, Point::Infinity) || !curve.is_on_curve(target) {
        return Err(format!("{kind} target is infinity or off curve"));
    }
    let (target_x, target_y) = affine_coordinates(target).expect("target was checked affine");
    let t = Fe::from_biguint(target_x);
    let a = Fe::from_biguint(&curve.a);
    let b = Fe::from_biguint(&curve.b);
    let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
    let inverse_exponent = &curve.p - BigUint::from(2u8);
    let mut solves = 0u64;
    let mut roots_returned = 0u64;
    let mut linear = 0u64;
    let mut universal = 0u64;
    let mut counts = MultiplicationCounts::default();
    let mut hits = BTreeSet::new();
    for (&left_key, &left_x) in &left.affine {
        solves += 1;
        let mut coefficient_field = CountedField::default();
        let (qa, qb, qc) = coefficients(left_x, t, a, b, &mut coefficient_field);
        counts.coefficients += coefficient_field.multiplications;
        let roots = solve_quadratic(qa, qb, qc, &sqrt_exponent, &inverse_exponent, &mut counts);
        linear += u64::from(roots.linear);
        universal += u64::from(roots.universal);
        if roots.universal {
            return Err(format!("{kind} query produced a universal quadratic"));
        }
        roots_returned += roots.len as u64;
        for root in roots.roots[..roots.len].iter().copied() {
            let right_key = fe_key(root);
            if right.affine.contains_key(&right_key) {
                let mut record = Vec::with_capacity(65);
                record.push(0);
                record.extend_from_slice(&left_key);
                record.extend_from_slice(&right_key);
                hits.insert(record);
            }
        }
    }
    let target_key = big_key(target_x);
    if left.identity && right.affine.contains_key(&target_key) {
        let mut record = Vec::with_capacity(33);
        record.push(1);
        record.extend_from_slice(&target_key);
        hits.insert(record);
    }
    if right.identity && left.affine.contains_key(&target_key) {
        let mut record = Vec::with_capacity(33);
        record.push(2);
        record.extend_from_slice(&target_key);
        hits.insert(record);
    }
    counts.finish();
    let candidate_positive = !hits.is_empty();
    let false_negative = reference_positive && !candidate_positive;
    let false_positive = !reference_positive && candidate_positive;
    let mut hit_bytes = Vec::new();
    for hit in &hits {
        hit_bytes.extend_from_slice(hit);
    }
    let cold_counts = build_counts.combined(counts);
    Ok(CompressedTargetResult {
        kind: kind.into(),
        scalar: scalar.map(lower_hex),
        target: [lower_hex(target_x), lower_hex(target_y)],
        reference_positive,
        candidate_positive,
        algebra_hits: hits.len() as u64,
        quadratic_solves_warm: solves,
        quadratic_solves_cold: build_solves + solves,
        quadratic_roots_returned: roots_returned,
        algebra_index_lookups: roots_returned,
        linear_degeneracies: linear,
        universal_degeneracies: universal,
        multiplication_counts_warm: counts,
        multiplication_counts_cold: cold_counts,
        false_negative,
        false_positive,
        exact: !false_negative && !false_positive,
        hit_sha256: hex::encode(sha256(&hit_bytes)),
    })
}

fn verify_round14_samples(
    dependency: &serde_json::Value,
    sample_specs: &[(String, Vec<u64>)],
    columns: &[Column],
) -> Result<(), String> {
    let samples = dependency
        .get("samples")
        .and_then(serde_json::Value::as_array)
        .ok_or("round-14 dependency has no samples array")?;
    if samples.len() != sample_specs.len() {
        return Err("round-14 dependency sample count changed".into());
    }
    for (sample, (expected_kind, expected_columns)) in samples.iter().zip(sample_specs) {
        if sample.get("kind").and_then(serde_json::Value::as_str) != Some(expected_kind) {
            return Err(format!("round-14 sample name changed for {expected_kind}"));
        }
        let stored_columns: Vec<u64> = sample
            .get("columns")
            .and_then(serde_json::Value::as_array)
            .ok_or_else(|| format!("round-14 {expected_kind} columns are missing"))?
            .iter()
            .map(|value| {
                value
                    .as_u64()
                    .ok_or_else(|| format!("round-14 {expected_kind} has a non-u64 column"))
            })
            .collect::<Result<_, _>>()?;
        if &stored_columns != expected_columns {
            return Err(format!("round-14 {expected_kind} column selection changed"));
        }
        let stored_xs: Vec<&str> = sample
            .get("x_coordinates")
            .and_then(serde_json::Value::as_array)
            .ok_or_else(|| format!("round-14 {expected_kind} x-coordinates are missing"))?
            .iter()
            .map(|value| {
                value
                    .as_str()
                    .ok_or_else(|| format!("round-14 {expected_kind} has a non-string x"))
            })
            .collect::<Result<_, _>>()?;
        let rebuilt_xs: Vec<String> = expected_columns
            .iter()
            .map(|&index| {
                lower_hex(&BigUint::from_bytes_be(
                    &columns[index as usize].x.to_bytes_be(),
                ))
            })
            .collect();
        if stored_xs
            .iter()
            .zip(&rebuilt_xs)
            .any(|(stored, rebuilt)| *stored != rebuilt)
        {
            return Err(format!("round-14 {expected_kind} x-coordinates changed"));
        }
    }
    Ok(())
}

fn run_compression_sample(
    kind: &str,
    column_indices: Vec<u64>,
    columns: &[Column],
    curve: &CurveParams,
) -> Result<CompressionSampleResult, String> {
    let selected: Vec<&Column> = column_indices
        .iter()
        .map(|&index| {
            columns
                .get(index as usize)
                .ok_or_else(|| format!("{kind} column {index} is out of range"))
        })
        .collect::<Result<_, _>>()?;
    if selected.len() != 16 {
        return Err(format!("{kind} requires exactly sixteen columns"));
    }
    let xs: Vec<Fe> = selected.iter().map(|column| column.x).collect();
    let points: Vec<Point> = selected.iter().map(|column| column.low.clone()).collect();
    let x_coordinates: Vec<String> = xs
        .iter()
        .map(|x| lower_hex(&BigUint::from_bytes_be(&x.to_bytes_be())))
        .collect();
    let a = Fe::from_biguint(&curve.a);
    let b = Fe::from_biguint(&curve.b);
    let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
    let inverse_exponent = &curve.p - BigUint::from(2u8);
    let mut build_solves = 0u64;
    let mut build_roots = 0u64;
    let mut build_linear = 0u64;
    let mut build_universal = 0u64;
    let mut build_counts = MultiplicationCounts::default();
    let mut intermediate_reference_group_additions = 0u64;
    let left = build_verified_eight_image(
        &format!("{kind}/left"),
        &xs[..8],
        &points[..8],
        curve,
        a,
        b,
        &sqrt_exponent,
        &inverse_exponent,
        &mut build_solves,
        &mut build_roots,
        &mut build_linear,
        &mut build_universal,
        &mut build_counts,
        &mut intermediate_reference_group_additions,
    )?;
    let right = build_verified_eight_image(
        &format!("{kind}/right"),
        &xs[8..],
        &points[8..],
        curve,
        a,
        b,
        &sqrt_exponent,
        &inverse_exponent,
        &mut build_solves,
        &mut build_roots,
        &mut build_linear,
        &mut build_universal,
        &mut build_counts,
        &mut intermediate_reference_group_additions,
    )?;
    build_counts.finish();
    if left.affine.len() != 128
        || right.affine.len() != 128
        || left.identity
        || right.identity
        || build_solves != 152
        || build_roots != 304
    {
        return Err(format!("{kind} half-image shape or count changed"));
    }

    let (full_reference, full_reference_group_additions) = reference_image(&points, curve)?;
    let curve_a = curve.a_fe();
    let planted = points.iter().fold(Point::Infinity, |sum, point| {
        sum.add_vartime(point, &curve_a)
    });
    if matches!(planted, Point::Infinity) || !curve.is_on_curve(&planted) {
        return Err(format!("{kind} planted target is infinity or off curve"));
    }
    let public_preimage = format!("{CURVE_SLUG}/s17-target-join-round15/sample/{kind}/public");
    let mut public_scalar = BigUint::from_bytes_be(&sha256(public_preimage.as_bytes())) % &curve.n;
    if public_scalar.is_zero() {
        public_scalar = BigUint::one();
    }
    let public = curve
        .generator()
        .scalar_mul_vartime(&public_scalar, &curve_a);

    let reference_membership = |target: &Point| -> Result<bool, String> {
        let (x, _) = affine_coordinates(target).ok_or("compressed target is infinity")?;
        Ok(full_reference.affine.contains(&big_key(x)))
    };
    let targets = vec![
        query_compressed_target(
            "planted-positive",
            None,
            &planted,
            reference_membership(&planted)?,
            &left,
            &right,
            curve,
            build_solves,
            build_counts,
        )?,
        query_compressed_target(
            "hash-public",
            Some(&public_scalar),
            &public,
            reference_membership(&public)?,
            &left,
            &right,
            curve,
            build_solves,
            build_counts,
        )?,
    ];
    if targets.iter().any(|target| !target.exact) {
        return Err(format!("{kind} compressed target classification mismatch"));
    }
    let retained_entries = left.affine.len() + right.affine.len();
    let materialized_sixteen_leaf_entries = 32_768usize;
    Ok(CompressionSampleResult {
        kind: kind.into(),
        columns: column_indices,
        x_coordinates,
        left_eight_leaf_entries: left.affine.len(),
        right_eight_leaf_entries: right.affine.len(),
        retained_entries,
        materialized_sixteen_leaf_entries,
        retained_raw_bytes: retained_entries * 32,
        materialized_raw_bytes: materialized_sixteen_leaf_entries * 32,
        retained_to_materialized_ratio: retained_entries as f64
            / materialized_sixteen_leaf_entries as f64,
        build_quadratic_solves: build_solves,
        full_materialization_quadratic_solves: 16_536,
        intermediate_reference_group_additions,
        full_reference_group_additions,
        build_linear_degeneracies: build_linear,
        build_universal_degeneracies: build_universal,
        build_multiplication_counts: build_counts,
        left_image_sha256: image_digest(&left),
        right_image_sha256: image_digest(&right),
        targets,
    })
}

fn run_compression(round14_path: &PathBuf, out: Option<PathBuf>) -> Result<(), String> {
    let bytes = std::fs::read(round14_path).map_err(|error| error.to_string())?;
    let round14_sha256 = hex::encode(sha256(&bytes));
    if round14_sha256 != ROUND14_SHA256 {
        return Err(format!(
            "round-14 SHA-256 mismatch: expected {ROUND14_SHA256}, got {round14_sha256}"
        ));
    }
    let dependency: serde_json::Value =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    let curve = CurveParams::p256();
    let (factor_base, columns, _, _, _) = build_indexes(&curve)?;
    let sample_specs = round14_sample_specs();
    verify_round14_samples(&dependency, &sample_specs, &columns)?;
    let mut samples = Vec::new();
    for (kind, indices) in sample_specs {
        samples.push(run_compression_sample(&kind, indices, &columns, &curve)?);
    }
    let result = CompressionExperimentResult {
        schema: "p256.s17_target_indexed_compression/v1".into(),
        curve: CURVE_SLUG.into(),
        field_prime: lower_hex(&curve.p),
        curve_a: lower_hex(&curve.a),
        curve_b: lower_hex(&curve.b),
        round14_sha256,
        factor_base,
        local_maximum_degree: 2,
        samples,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn validate_atom_columns(column_indices: &[u64]) -> Result<(), String> {
    if column_indices.len() != ATOM_COLUMNS {
        return Err(format!(
            "packed atom requires {ATOM_COLUMNS} columns, got {}",
            column_indices.len()
        ));
    }
    let mut seen = BTreeSet::new();
    for &column in column_indices {
        if column >= COLUMNS {
            return Err(format!(
                "packed atom column {column} is outside the {COLUMNS}-column dictionary"
            ));
        }
        if !seen.insert(column) {
            return Err(format!("packed atom repeats column {column}"));
        }
    }
    Ok(())
}

fn pack_atom_columns(column_indices: &[u64]) -> Result<[u8; PACKED_ATOM_BYTES], String> {
    validate_atom_columns(column_indices)?;
    let mut packet = [0u8; PACKED_ATOM_BYTES];
    let mut cursor = 0usize;
    for &column in column_indices {
        for shift in (0..COLUMN_INDEX_BITS).rev() {
            let bit = ((column >> shift) & 1) as u8;
            packet[cursor / 8] |= bit << (7 - cursor % 8);
            cursor += 1;
        }
    }
    if cursor != PACKED_ATOM_BYTES * 8 {
        return Err("packed atom bit count changed".into());
    }
    Ok(packet)
}

fn unpack_atom_columns(packet: &[u8]) -> Result<Vec<u64>, String> {
    if packet.len() != PACKED_ATOM_BYTES {
        return Err(format!(
            "packed atom requires {PACKED_ATOM_BYTES} bytes, got {}",
            packet.len()
        ));
    }
    let mut columns = Vec::with_capacity(ATOM_COLUMNS);
    let mut cursor = 0usize;
    for _ in 0..ATOM_COLUMNS {
        let mut column = 0u64;
        for _ in 0..COLUMN_INDEX_BITS {
            let bit = (packet[cursor / 8] >> (7 - cursor % 8)) & 1;
            column = (column << 1) | u64::from(bit);
            cursor += 1;
        }
        columns.push(column);
    }
    if cursor != packet.len() * 8 {
        return Err("packed atom decoder did not consume the packet".into());
    }
    validate_atom_columns(&columns)?;
    Ok(columns)
}

fn round15_sample<'a>(
    dependency: &'a serde_json::Value,
    kind: &str,
) -> Result<&'a serde_json::Value, String> {
    dependency["samples"]
        .as_array()
        .ok_or("round-15 samples are missing")?
        .iter()
        .find(|sample| sample["kind"].as_str() == Some(kind))
        .ok_or_else(|| format!("round-15 sample {kind} is missing"))
}

fn run_packed_atom_sample(
    kind: &str,
    expected_columns: Vec<u64>,
    dependency: &serde_json::Value,
    columns: &[Column],
    curve: &CurveParams,
) -> Result<PackedAtomSampleResult, String> {
    let stored_sample = round15_sample(dependency, kind)?;
    let stored_columns: Vec<u64> = serde_json::from_value(stored_sample["columns"].clone())
        .map_err(|error| format!("round-15 {kind} columns: {error}"))?;
    if stored_columns != expected_columns {
        return Err(format!("round-15 {kind} column selection changed"));
    }

    let packet = pack_atom_columns(&stored_columns)?;
    let decoded = unpack_atom_columns(&packet)?;
    let decoded_exact = decoded == stored_columns;
    let reencoded_exact = pack_atom_columns(&decoded)? == packet;
    if !decoded_exact || !reencoded_exact {
        return Err(format!("{kind} packed atom failed canonical round trip"));
    }

    let recomputed = run_compression_sample(kind, decoded, columns, curve)?;
    let recomputed_value = serde_json::to_value(&recomputed).map_err(|error| error.to_string())?;
    if &recomputed_value != stored_sample {
        return Err(format!(
            "{kind} reconstructed round-15 sample differs from its frozen receipt"
        ));
    }
    if recomputed.targets.len() != 2 || recomputed.targets.iter().any(|target| !target.exact) {
        return Err(format!("{kind} reconstructed targets are not both exact"));
    }
    let one_cold_query_quadratic_solves = recomputed.targets[0].quadratic_solves_cold;
    if recomputed
        .targets
        .iter()
        .any(|target| target.quadratic_solves_cold != one_cold_query_quadratic_solves)
    {
        return Err(format!("{kind} cold-query solve counts differ"));
    }
    let two_query_one_rebuild_quadratic_solves = recomputed.build_quadratic_solves
        + recomputed
            .targets
            .iter()
            .map(|target| target.quadratic_solves_warm)
            .sum::<u64>();
    if packet.len() != 36
        || recomputed.retained_raw_bytes != 8_192
        || recomputed.materialized_raw_bytes != 1_048_576
        || one_cold_query_quadratic_solves != 280
        || two_query_one_rebuild_quadratic_solves != 408
        || recomputed.build_linear_degeneracies != 0
        || recomputed.build_universal_degeneracies != 0
    {
        return Err(format!("{kind} packed-atom boundary changed"));
    }
    let additional_persistent_reduction =
        recomputed.retained_raw_bytes as f64 / packet.len() as f64;
    if additional_persistent_reduction < 30.0 {
        return Err(format!(
            "{kind} compression gate failed: {additional_persistent_reduction:.6}x"
        ));
    }

    Ok(PackedAtomSampleResult {
        kind: recomputed.kind,
        packed_atom_hex: hex::encode(packet),
        packed_atom_sha256: hex::encode(sha256(&packet)),
        bits_per_column: COLUMN_INDEX_BITS,
        packet_bytes: packet.len(),
        round15_retained_raw_bytes: recomputed.retained_raw_bytes,
        materialized_sixteen_leaf_raw_bytes: recomputed.materialized_raw_bytes,
        additional_persistent_reduction,
        full_materialized_reduction: recomputed.materialized_raw_bytes as f64 / packet.len() as f64,
        columns: recomputed.columns,
        x_coordinates: recomputed.x_coordinates,
        decoded_exact,
        reencoded_exact,
        dependency_sample_exact: true,
        left_image_sha256: recomputed.left_image_sha256,
        right_image_sha256: recomputed.right_image_sha256,
        transient_reconstruction_raw_bytes: recomputed.retained_raw_bytes,
        build_quadratic_solves: recomputed.build_quadratic_solves,
        one_cold_query_quadratic_solves,
        two_query_one_rebuild_quadratic_solves,
        intermediate_reference_group_additions: recomputed.intermediate_reference_group_additions,
        full_reference_group_additions: recomputed.full_reference_group_additions,
        build_linear_degeneracies: recomputed.build_linear_degeneracies,
        build_universal_degeneracies: recomputed.build_universal_degeneracies,
        build_multiplication_counts: recomputed.build_multiplication_counts,
        targets: recomputed.targets,
    })
}

fn run_packed_compression(round15_path: &PathBuf, out: Option<PathBuf>) -> Result<(), String> {
    let bytes = std::fs::read(round15_path).map_err(|error| error.to_string())?;
    let round15_sha256 = hex::encode(sha256(&bytes));
    if round15_sha256 != ROUND15_SHA256 {
        return Err(format!(
            "round-15 SHA-256 mismatch: expected {ROUND15_SHA256}, got {round15_sha256}"
        ));
    }
    let dependency: serde_json::Value =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if dependency["schema"].as_str() != Some("p256.s17_target_indexed_compression/v1")
        || dependency["curve"].as_str() != Some(CURVE_SLUG)
        || dependency["factor_base"]["fb_id"].as_str() != Some(FB_ID)
        || dependency["factor_base"]["fb_sha256"].as_str() != Some(FB_SHA256)
        || dependency["factor_base"]["points_sha256"].as_str() != Some(POINTS_SHA256)
    {
        return Err("round-15 identity receipt changed".into());
    }

    let curve = CurveParams::p256();
    let (factor_base, columns, _, _, _) = build_indexes(&curve)?;
    let mut samples = Vec::new();
    for (kind, expected_columns) in round14_sample_specs() {
        samples.push(run_packed_atom_sample(
            &kind,
            expected_columns,
            &dependency,
            &columns,
            &curve,
        )?);
    }
    let result = PackedAtomExperimentResult {
        schema: "p256.s17_packed_atom/v1".into(),
        curve: CURVE_SLUG.into(),
        field_prime: lower_hex(&curve.p),
        curve_a: lower_hex(&curve.a),
        curve_b: lower_hex(&curve.b),
        round15_sha256,
        factor_base,
        local_maximum_degree: 2,
        factor_base_column_bits: COLUMN_INDEX_BITS,
        packed_atom_bytes: PACKED_ATOM_BYTES,
        samples,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn run_transfer(out: Option<PathBuf>) -> Result<(), String> {
    let curve = CurveParams::p256();
    let (factor_base, columns, signed_rows, x_index, signed_index) = build_indexes(&curve)?;
    let a = curve.a_fe();
    let planted = columns[0].low.add_vartime(&columns[1].low, &a);
    if matches!(planted, Point::Infinity) || !curve.is_on_curve(&planted) {
        return Err("planted target is infinity or off curve".into());
    }

    let mut hash_scalar = BigUint::from_bytes_be(&sha256(TARGET_PREIMAGE.as_bytes())) % &curve.n;
    if hash_scalar.is_zero() {
        hash_scalar = BigUint::one();
    }
    let hash_target = curve
        .generator()
        .scalar_mul_vartime(&hash_scalar, &curve.a_fe());

    let targets = vec![
        run_target(
            "planted-positive",
            None,
            &planted,
            &columns,
            &signed_rows,
            &x_index,
            &signed_index,
            &curve,
        )?,
        run_target(
            "hash-public",
            Some(&hash_scalar),
            &hash_target,
            &columns,
            &signed_rows,
            &x_index,
            &signed_index,
            &curve,
        )?,
    ];
    if targets.iter().any(|target| !target.exact) {
        return Err("algebraic and signed-point pair sets differ".into());
    }
    let result = ExperimentResult {
        schema: "p256.s3_image_transfer/v1".into(),
        curve: CURVE_SLUG.into(),
        field_prime: lower_hex(&curve.p),
        curve_a: lower_hex(&curve.a),
        curve_b: lower_hex(&curve.b),
        factor_base,
        target_preimage: TARGET_PREIMAGE.into(),
        targets,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn run(cli: Cli) -> Result<(), String> {
    let Cli {
        out,
        width_round13,
        compress_round14,
        pack_round15,
    } = cli;
    let continuation_modes = usize::from(width_round13.is_some())
        + usize::from(compress_round14.is_some())
        + usize::from(pack_round15.is_some());
    if continuation_modes > 1 {
        return Err("choose only one continuation mode".into());
    }
    if let Some(path) = width_round13 {
        return run_width(&path, out);
    }
    if let Some(path) = compress_round14 {
        return run_compression(&path, out);
    }
    if let Some(path) = pack_round15 {
        return run_packed_compression(&path, out);
    }
    run_transfer(out)
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

    fn s3(u: Fe, v: Fe, t: Fe, a: Fe, b: Fe) -> Fe {
        let q = u.mul(&v);
        let sum = u.add(&v);
        let difference = u.sub(&v);
        let first = t.sqr().mul(&difference.sqr());
        let inner = sum.mul(&q.add(&a)).add(&b.add(&b));
        let second = t.add(&t).mul(&inner).neg();
        let third = q.sub(&a).sqr().sub(&b.add(&b).add(&b.add(&b)).mul(&sum));
        first.add(&second).add(&third)
    }

    #[test]
    fn coefficients_equal_direct_s3() {
        let curve = CurveParams::p256();
        let a = Fe::from_biguint(&curve.a);
        let b = Fe::from_biguint(&curve.b);
        for (u, t) in [(1u64, 2u64), (17, 301), (500, 500), (1_000_003, 42)] {
            let u = Fe::from_biguint(&BigUint::from(u));
            let t = Fe::from_biguint(&BigUint::from(t));
            let mut count = CountedField::default();
            let (qa, qb, qc) = coefficients(u, t, a, b, &mut count);
            for v in [0u64, 1, 2, 99, 1024, 1_000_033] {
                let v = Fe::from_biguint(&BigUint::from(v));
                let polynomial = qa.mul(&v.sqr()).add(&qb.mul(&v)).add(&qc);
                assert_eq!(polynomial, s3(u, v, t, a, b));
            }
        }
    }

    #[test]
    fn quadratic_solver_recovers_constructed_roots_and_linear_case() {
        let curve = CurveParams::p256();
        let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
        let inverse_exponent = &curve.p - BigUint::from(2u8);
        let r1 = Fe::from_biguint(&BigUint::from(123u16));
        let r2 = Fe::from_biguint(&BigUint::from(456u16));
        let qa = Fe::ONE;
        let qb = r1.add(&r2).neg();
        let qc = r1.mul(&r2);
        let mut counts = MultiplicationCounts::default();
        let roots = solve_quadratic(qa, qb, qc, &sqrt_exponent, &inverse_exponent, &mut counts);
        let actual: BTreeSet<[u8; 32]> = roots.roots[..roots.len]
            .iter()
            .copied()
            .map(fe_key)
            .collect();
        let expected = BTreeSet::from([fe_key(r1), fe_key(r2)]);
        assert_eq!(actual, expected);

        let mut counts = MultiplicationCounts::default();
        let linear = solve_quadratic(
            Fe::ZERO,
            Fe::from_biguint(&BigUint::from(7u8)),
            Fe::from_biguint(&BigUint::from(35u8)),
            &sqrt_exponent,
            &inverse_exponent,
            &mut counts,
        );
        assert!(linear.linear);
        assert_eq!(linear.len, 1);
        let seven = Fe::from_biguint(&BigUint::from(7u8));
        let thirty_five = Fe::from_biguint(&BigUint::from(35u8));
        assert_eq!(seven.mul(&linear.roots[0]).add(&thirty_five), Fe::ZERO);
    }

    #[test]
    fn packed_atom_is_canonical_and_round_trips() {
        let columns: Vec<u64> = (0..ATOM_COLUMNS as u64)
            .map(|index| index * 8_001)
            .collect();
        let packet = pack_atom_columns(&columns).unwrap();
        assert_eq!(packet.len(), 36);
        assert_eq!(unpack_atom_columns(&packet).unwrap(), columns);

        let mut aggregate = BigUint::zero();
        for column in &columns {
            aggregate = (aggregate << COLUMN_INDEX_BITS) + BigUint::from(*column);
        }
        let encoded = aggregate.to_bytes_be();
        let mut reference = [0u8; PACKED_ATOM_BYTES];
        reference[PACKED_ATOM_BYTES - encoded.len()..].copy_from_slice(&encoded);
        assert_eq!(packet, reference);
    }

    #[test]
    fn packed_atom_rejects_noncanonical_columns() {
        let mut duplicate: Vec<u64> = (0..ATOM_COLUMNS as u64).collect();
        duplicate[15] = duplicate[0];
        assert!(pack_atom_columns(&duplicate).is_err());

        let mut out_of_range: Vec<u64> = (0..ATOM_COLUMNS as u64).collect();
        out_of_range[15] = COLUMNS;
        assert!(pack_atom_columns(&out_of_range).is_err());
        assert!(unpack_atom_columns(&[0u8; PACKED_ATOM_BYTES - 1]).is_err());
    }
}
