//! Transfer the indexed final-S3 image gate to the committed P-256 factor base.

use std::collections::{BTreeSet, HashMap};
use std::fmt::Write as _;
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{self, CURVE_SLUG};
use crypto_lib::ecc::p256_field::P256FieldElement as Fe;
use crypto_lib::ecc::{CurveParams, Point};
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::Serialize;

const SPEC: &str = "dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129";
const FB_ID: &str = "FB1h2f8621cda105";
const FB_SHA256: &str = "2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42";
const POINTS_SHA256: &str = "70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1";
const TERMINAL: &str = "0x5b17195299a3158b93389ad04c776fff2a8bb23ca8659b5b0c2b75f9009b65e5";
const COLUMNS: u64 = 131_458;
const SIGNED_POINTS: u64 = 262_916;
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

fn run(cli: Cli) -> Result<(), String> {
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
    match cli.out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
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
}
