//! Exact SMT encoding for two-point decomposition over a prime field.
//!
//! This module is deliberately separate from [`super::semaev_sat`]. A
//! literal prime field has no extension coordinates for a Weil descent, so
//! the native object is modular arithmetic rather than a Boolean ANF. The
//! exporter offers three exact backends: overflow-safe `QF_BV`, native
//! `QF_FF` coordinates, and native `QF_FF` with the Semaev `S3`
//! polynomial. The latter constrains abscissae to the exact liftable factor
//! base and recovers signs with the native checker.
//!
//! The encoded question is whether a finite target `R` on
//! `y^2 = x^3 + a*x + b` is `P1 + P2`, with both abscissae below a fixed
//! factor-base bound. Both ordinary addition and doubling are present.
//! Every SAT model can be checked independently with [`verify_witness`].

use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::collections::{HashMap, HashSet};
use std::fmt::{self, Write as _};
use std::fs::File;
use std::hash::Hash;
use std::io::{self, Read};
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

/// JSON schema emitted and accepted by this module.
pub const INSTANCE_SCHEMA: &str = "prime-field-smt.instance/v1";

/// A non-negative integer serialized as an exact decimal string.
///
/// Input additionally accepts JSON integers and `0x`/`0b` strings. Using a
/// string on output avoids the 53-bit precision limit of generic JSON tools.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ExactUint(pub BigUint);

impl ExactUint {
    pub fn parse(value: &str) -> Result<Self, String> {
        let value = value.trim();
        let (digits, radix) = value
            .strip_prefix("0x")
            .or_else(|| value.strip_prefix("0X"))
            .map_or_else(
                || {
                    value
                        .strip_prefix("0b")
                        .or_else(|| value.strip_prefix("0B"))
                        .map_or((value, 10), |digits| (digits, 2))
                },
                |digits| (digits, 16),
            );
        if digits.is_empty() {
            return Err("empty integer".to_string());
        }
        BigUint::parse_bytes(digits.as_bytes(), radix)
            .map(Self)
            .ok_or_else(|| format!("invalid base-{radix} non-negative integer: {value}"))
    }
}

impl From<u64> for ExactUint {
    fn from(value: u64) -> Self {
        Self(BigUint::from(value))
    }
}

impl fmt::Display for ExactUint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl Serialize for ExactUint {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&self.0.to_str_radix(10))
    }
}

struct ExactUintVisitor;

impl<'de> Visitor<'de> for ExactUintVisitor {
    type Value = ExactUint;

    fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("a non-negative JSON integer or decimal/0x/0b string")
    }

    fn visit_u64<E>(self, value: u64) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(value.into())
    }

    fn visit_u128<E>(self, value: u128) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(ExactUint(BigUint::from(value)))
    }

    fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        ExactUint::parse(value).map_err(E::custom)
    }

    fn visit_string<E>(self, value: String) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        self.visit_str(&value)
    }
}

impl<'de> Deserialize<'de> for ExactUint {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_any(ExactUintVisitor)
    }
}

/// A short-Weierstrass curve over an odd prime field.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PrimeFieldCurve {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    pub p: ExactUint,
    pub a: ExactUint,
    pub b: ExactUint,
}

/// A finite affine point.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct AffinePoint {
    pub x: ExactUint,
    pub y: ExactUint,
}

/// One exact two-point decomposition query.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PrimeFieldSmtInstance {
    pub schema: String,
    pub curve: PrimeFieldCurve,
    pub target: AffinePoint,
    pub x_bound: ExactUint,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<serde_json::Value>,
}

/// A complete satisfying assignment returned by an SMT solver.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DecompositionWitness {
    pub x1: ExactUint,
    pub y1: ExactUint,
    pub x2: ExactUint,
    pub y2: ExactUint,
    pub lambda: ExactUint,
}

impl DecompositionWitness {
    pub fn p1(&self) -> AffinePoint {
        AffinePoint {
            x: self.x1.clone(),
            y: self.y1.clone(),
        }
    }

    pub fn p2(&self) -> AffinePoint {
        AffinePoint {
            x: self.x2.clone(),
            y: self.y2.clone(),
        }
    }
}

/// Static checks performed before an instance is encoded or evaluated.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct InstanceValidation {
    pub schema: String,
    pub field_bits: u64,
    pub multiplication_bits: u64,
    pub primality_assurance: String,
    pub nonsingular: bool,
    pub target_on_curve: bool,
    pub valid: bool,
}

fn add_mod(lhs: &BigUint, rhs: &BigUint, p: &BigUint) -> BigUint {
    (lhs + rhs) % p
}

fn sub_mod(lhs: &BigUint, rhs: &BigUint, p: &BigUint) -> BigUint {
    if lhs >= rhs {
        (lhs - rhs) % p
    } else {
        (lhs + p - rhs) % p
    }
}

fn mul_mod(lhs: &BigUint, rhs: &BigUint, p: &BigUint) -> BigUint {
    (lhs * rhs) % p
}

fn inverse_mod(value: &BigUint, p: &BigUint) -> Option<BigUint> {
    if value.is_zero() {
        None
    } else {
        Some(value.modpow(&(p - BigUint::from(2u8)), p))
    }
}

fn sqrt_mod_prime(value: &BigUint, p: &BigUint) -> Option<BigUint> {
    let value = value % p;
    if value.is_zero() {
        return Some(BigUint::zero());
    }
    let one = BigUint::one();
    let p_minus_one = p - &one;
    if value.modpow(&(&p_minus_one >> 1usize), p) != one {
        return None;
    }
    if (p & BigUint::from(3u8)) == BigUint::from(3u8) {
        return Some(value.modpow(&((p + &one) >> 2usize), p));
    }

    let mut q = p_minus_one.clone();
    let mut s = 0u64;
    while (&q & &one).is_zero() {
        q >>= 1usize;
        s += 1;
    }
    let mut z = BigUint::from(2u8);
    while z.modpow(&(&p_minus_one >> 1usize), p) != p_minus_one {
        z += &one;
    }
    let mut m = s;
    let mut c = z.modpow(&q, p);
    let mut t = value.modpow(&q, p);
    let mut r = value.modpow(&((&q + &one) >> 1usize), p);
    while t != one {
        let mut i = 0u64;
        let mut power = t.clone();
        while power != one {
            power = (&power * &power) % p;
            i += 1;
            if i == m {
                return None;
            }
        }
        let exponent = BigUint::one() << (m - i - 1) as usize;
        let b = c.modpow(&exponent, p);
        let b2 = (&b * &b) % p;
        r = (r * &b) % p;
        t = (t * &b2) % p;
        c = b2;
        m = i;
    }
    Some(r)
}

fn miller_rabin_round(n: &BigUint, d: &BigUint, s: u32, base: u64) -> bool {
    let n_minus_one = n - BigUint::one();
    let a = BigUint::from(base) % n;
    if a.is_zero() {
        return true;
    }
    let mut x = a.modpow(d, n);
    if x.is_one() || x == n_minus_one {
        return true;
    }
    for _ in 1..s {
        x = (&x * &x) % n;
        if x == n_minus_one {
            return true;
        }
    }
    false
}

fn probable_prime(value: &BigUint) -> bool {
    if let Some(value) = value.to_u64() {
        return super::residual_walk::is_prime_u64(value);
    }
    const BASES: [u64; 32] = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97, 101, 103, 107, 109, 113, 127, 131,
    ];
    if value < &BigUint::from(2u8) || (value & BigUint::one()).is_zero() {
        return false;
    }
    for prime in BASES {
        let prime = BigUint::from(prime);
        if value == &prime {
            return true;
        }
        if (value % &prime).is_zero() {
            return false;
        }
    }
    let mut d = value - BigUint::one();
    let mut s = 0u32;
    while (&d & BigUint::one()).is_zero() {
        d >>= 1usize;
        s += 1;
    }
    BASES
        .iter()
        .all(|&base| miller_rabin_round(value, &d, s, base))
}

fn point_on_curve(curve: &PrimeFieldCurve, point: &AffinePoint) -> bool {
    let p = &curve.p.0;
    if point.x.0 >= *p || point.y.0 >= *p {
        return false;
    }
    let y2 = mul_mod(&point.y.0, &point.y.0, p);
    let x2 = mul_mod(&point.x.0, &point.x.0, p);
    let x3 = mul_mod(&x2, &point.x.0, p);
    let ax = mul_mod(&curve.a.0, &point.x.0, p);
    y2 == add_mod(&add_mod(&x3, &ax, p), &curve.b.0, p)
}

/// Validate the field, curve, target, and factor-base bound.
pub fn validate_instance(instance: &PrimeFieldSmtInstance) -> Result<InstanceValidation, String> {
    if instance.schema != INSTANCE_SCHEMA {
        return Err(format!(
            "unsupported schema {:?}; expected {INSTANCE_SCHEMA}",
            instance.schema
        ));
    }
    let p = &instance.curve.p.0;
    if p <= &BigUint::from(3u8) || (p & BigUint::one()).is_zero() {
        return Err("p must be an odd prime greater than 3".to_string());
    }
    if p.bits() > 4096 {
        return Err(format!(
            "{}-bit fields exceed the 4096-bit exporter guard",
            p.bits()
        ));
    }
    if !probable_prime(p) {
        return Err("p failed the prime check".to_string());
    }
    if instance.curve.a.0 >= *p || instance.curve.b.0 >= *p {
        return Err("curve coefficients must be canonical field elements".to_string());
    }
    if instance.x_bound.0.is_zero() || instance.x_bound.0 > *p {
        return Err("x_bound must satisfy 1 <= x_bound <= p".to_string());
    }
    if instance.target.x.0 >= *p || instance.target.y.0 >= *p {
        return Err("target coordinates must be canonical field elements".to_string());
    }

    let a2 = mul_mod(&instance.curve.a.0, &instance.curve.a.0, p);
    let a3 = mul_mod(&a2, &instance.curve.a.0, p);
    let b2 = mul_mod(&instance.curve.b.0, &instance.curve.b.0, p);
    let discriminant_term = add_mod(
        &mul_mod(&BigUint::from(4u8), &a3, p),
        &mul_mod(&BigUint::from(27u8), &b2, p),
        p,
    );
    let nonsingular = !discriminant_term.is_zero();
    if !nonsingular {
        return Err("curve is singular: 4*a^3 + 27*b^2 = 0 mod p".to_string());
    }
    let target_on_curve = point_on_curve(&instance.curve, &instance.target);
    if !target_on_curve {
        return Err("target is not on the curve".to_string());
    }
    Ok(InstanceValidation {
        schema: "prime-field-smt.validation/v1".to_string(),
        field_bits: p.bits(),
        multiplication_bits: p.bits() * 2,
        primality_assurance: if p.bits() <= 64 {
            "deterministic-u64-miller-rabin".to_string()
        } else {
            "probable-prime-32-fixed-miller-rabin-bases".to_string()
        },
        nonsingular,
        target_on_curve,
        valid: true,
    })
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum GroupPoint {
    Infinity,
    Affine(AffinePoint),
}

fn add_points(
    curve: &PrimeFieldCurve,
    left: &AffinePoint,
    right: &AffinePoint,
) -> Result<(GroupPoint, Option<BigUint>), String> {
    let p = &curve.p.0;
    if !point_on_curve(curve, left) || !point_on_curve(curve, right) {
        return Err("group addition received a point off the curve".to_string());
    }
    if left.x == right.x {
        if add_mod(&left.y.0, &right.y.0, p).is_zero() {
            return Ok((GroupPoint::Infinity, None));
        }
        if left.y != right.y {
            return Err("equal x-coordinates have neither equal nor inverse y".to_string());
        }
        let denominator = add_mod(&left.y.0, &left.y.0, p);
        let numerator = add_mod(
            &mul_mod(&BigUint::from(3u8), &mul_mod(&left.x.0, &left.x.0, p), p),
            &curve.a.0,
            p,
        );
        let inverse = inverse_mod(&denominator, p)
            .ok_or_else(|| "doubling denominator is zero".to_string())?;
        let lambda = mul_mod(&numerator, &inverse, p);
        let x3 = sub_mod(
            &mul_mod(&lambda, &lambda, p),
            &add_mod(&left.x.0, &left.x.0, p),
            p,
        );
        let y3 = sub_mod(
            &mul_mod(&lambda, &sub_mod(&left.x.0, &x3, p), p),
            &left.y.0,
            p,
        );
        return Ok((
            GroupPoint::Affine(AffinePoint {
                x: ExactUint(x3),
                y: ExactUint(y3),
            }),
            Some(lambda),
        ));
    }

    let denominator = sub_mod(&right.x.0, &left.x.0, p);
    let numerator = sub_mod(&right.y.0, &left.y.0, p);
    let inverse = inverse_mod(&denominator, p)
        .ok_or_else(|| "ordinary-addition denominator is zero".to_string())?;
    let lambda = mul_mod(&numerator, &inverse, p);
    let x3 = sub_mod(
        &sub_mod(&mul_mod(&lambda, &lambda, p), &left.x.0, p),
        &right.x.0,
        p,
    );
    let y3 = sub_mod(
        &mul_mod(&lambda, &sub_mod(&left.x.0, &x3, p), p),
        &left.y.0,
        p,
    );
    Ok((
        GroupPoint::Affine(AffinePoint {
            x: ExactUint(x3),
            y: ExactUint(y3),
        }),
        Some(lambda),
    ))
}

/// Independent checks applied to one parsed solver model.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct VerificationReceipt {
    pub schema: String,
    pub branch: Option<String>,
    pub canonical_field_values: bool,
    pub factor_base_bounds: bool,
    pub symmetry: bool,
    pub points_on_curve: bool,
    pub encoded_branch_equations: bool,
    pub independent_sum_finite: bool,
    pub independent_sum_matches_target: bool,
    pub lambda_matches_independent_addition: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub recomputed_sum: Option<AffinePoint>,
    pub verified: bool,
}

/// Recompute the original curve relation and group law for a SAT model.
pub fn verify_witness(
    instance: &PrimeFieldSmtInstance,
    witness: &DecompositionWitness,
) -> Result<VerificationReceipt, String> {
    validate_instance(instance)?;
    let p = &instance.curve.p.0;
    let values = [
        &witness.x1.0,
        &witness.y1.0,
        &witness.x2.0,
        &witness.y2.0,
        &witness.lambda.0,
    ];
    let canonical_field_values = values.iter().all(|value| *value < p);
    let factor_base_bounds = witness.x1.0 < instance.x_bound.0 && witness.x2.0 < instance.x_bound.0;
    let symmetry = witness.x1.0 <= witness.x2.0;
    let p1 = witness.p1();
    let p2 = witness.p2();
    let points_on_curve =
        point_on_curve(&instance.curve, &p1) && point_on_curve(&instance.curve, &p2);

    let lambda2 = mul_mod(&witness.lambda.0, &witness.lambda.0, p);
    let (branch, encoded_branch_equations) = if witness.x1.0 < witness.x2.0 {
        let slope = mul_mod(
            &witness.lambda.0,
            &sub_mod(&witness.x2.0, &witness.x1.0, p),
            p,
        ) == sub_mod(&witness.y2.0, &witness.y1.0, p);
        let x_equation =
            instance.target.x.0 == sub_mod(&sub_mod(&lambda2, &witness.x1.0, p), &witness.x2.0, p);
        let y_equation = instance.target.y.0
            == sub_mod(
                &mul_mod(
                    &witness.lambda.0,
                    &sub_mod(&witness.x1.0, &instance.target.x.0, p),
                    p,
                ),
                &witness.y1.0,
                p,
            );
        (
            Some("ordinary".to_string()),
            slope && x_equation && y_equation,
        )
    } else if witness.x1 == witness.x2 && witness.y1 == witness.y2 && !witness.y1.0.is_zero() {
        let slope = mul_mod(
            &add_mod(&witness.y1.0, &witness.y1.0, p),
            &witness.lambda.0,
            p,
        ) == add_mod(
            &mul_mod(
                &BigUint::from(3u8),
                &mul_mod(&witness.x1.0, &witness.x1.0, p),
                p,
            ),
            &instance.curve.a.0,
            p,
        );
        let x_equation =
            instance.target.x.0 == sub_mod(&lambda2, &add_mod(&witness.x1.0, &witness.x1.0, p), p);
        let y_equation = instance.target.y.0
            == sub_mod(
                &mul_mod(
                    &witness.lambda.0,
                    &sub_mod(&witness.x1.0, &instance.target.x.0, p),
                    p,
                ),
                &witness.y1.0,
                p,
            );
        (
            Some("doubling".to_string()),
            slope && x_equation && y_equation,
        )
    } else {
        (None, false)
    };

    let (recomputed_sum, independent_lambda) = if points_on_curve {
        match add_points(&instance.curve, &p1, &p2)? {
            (GroupPoint::Affine(point), lambda) => (Some(point), lambda),
            (GroupPoint::Infinity, lambda) => (None, lambda),
        }
    } else {
        (None, None)
    };
    let independent_sum_finite = recomputed_sum.is_some();
    let independent_sum_matches_target = recomputed_sum.as_ref() == Some(&instance.target);
    let lambda_matches_independent_addition =
        independent_lambda.as_ref() == Some(&witness.lambda.0);
    let verified = canonical_field_values
        && factor_base_bounds
        && symmetry
        && points_on_curve
        && encoded_branch_equations
        && independent_sum_finite
        && independent_sum_matches_target
        && lambda_matches_independent_addition;
    Ok(VerificationReceipt {
        schema: "prime-field-smt.verification/v1".to_string(),
        branch,
        canonical_field_values,
        factor_base_bounds,
        symmetry,
        points_on_curve,
        encoded_branch_equations,
        independent_sum_finite,
        independent_sum_matches_target,
        lambda_matches_independent_addition,
        recomputed_sum,
        verified,
    })
}

/// Size and integrity metadata for an emitted SMT-LIB query.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SmtEncodingReceipt {
    pub schema: String,
    pub logic: String,
    pub relation: String,
    pub field_bits: u64,
    pub multiplication_bits: u64,
    pub declared_variables: u32,
    pub field_helpers: u32,
    pub addition_branches: u32,
    pub assertion_count: u32,
    pub canonical_range_constraints: bool,
    pub widened_arithmetic: bool,
    pub symmetry_breaking: String,
    pub canonical_instance_blake3: String,
    pub query_bytes: u64,
    pub query_blake3: String,
}

/// SMT-LIB text and its receipt.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmtQuery {
    pub text: String,
    pub receipt: SmtEncodingReceipt,
}

fn bv(value: &BigUint, width: u64) -> String {
    format!("(_ bv{} {width})", value.to_str_radix(10))
}

fn blake3_hex(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}

/// Emit one portable, overflow-safe `QF_BV` query.
pub fn emit_smt2(instance: &PrimeFieldSmtInstance) -> Result<SmtQuery, String> {
    let validation = validate_instance(instance)?;
    let width = validation.field_bits;
    let wide = width + 1;
    let product = width * 2;
    let p = &instance.curve.p.0;
    let p_narrow = bv(p, width);
    let p_wide = bv(p, wide);
    let p_product = bv(p, product);
    let a = bv(&instance.curve.a.0, width);
    let b = bv(&instance.curve.b.0, width);
    let target_x = bv(&instance.target.x.0, width);
    let target_y = bv(&instance.target.y.0, width);
    let bound = bv(&instance.x_bound.0, width);
    let zero = bv(&BigUint::zero(), width);
    let three = bv(&BigUint::from(3u8), width);
    let mut out = String::new();
    writeln!(out, "; generated by crypto::cryptanalysis::prime_field_smt").unwrap();
    writeln!(out, "; exact finite two-point decomposition over F_p").unwrap();
    writeln!(out, "(set-logic QF_BV)").unwrap();
    writeln!(out, "(set-option :produce-models true)").unwrap();
    writeln!(
        out,
        "(define-fun fp.add ((lhs (_ BitVec {width})) (rhs (_ BitVec {width}))) (_ BitVec {width})"
    )
    .unwrap();
    writeln!(
        out,
        "  ((_ extract {} 0) (bvurem (bvadd ((_ zero_extend 1) lhs) ((_ zero_extend 1) rhs)) {p_wide})))",
        width - 1
    )
    .unwrap();
    writeln!(
        out,
        "(define-fun fp.sub ((lhs (_ BitVec {width})) (rhs (_ BitVec {width}))) (_ BitVec {width})"
    )
    .unwrap();
    writeln!(
        out,
        "  ((_ extract {} 0) (bvurem (bvsub (bvadd ((_ zero_extend 1) lhs) {p_wide}) ((_ zero_extend 1) rhs)) {p_wide})))",
        width - 1
    )
    .unwrap();
    writeln!(
        out,
        "(define-fun fp.mul ((lhs (_ BitVec {width})) (rhs (_ BitVec {width}))) (_ BitVec {width})"
    )
    .unwrap();
    writeln!(
        out,
        "  ((_ extract {} 0) (bvurem (bvmul ((_ zero_extend {width}) lhs) ((_ zero_extend {width}) rhs)) {p_product})))",
        width - 1
    )
    .unwrap();
    for name in ["x1", "y1", "x2", "y2", "lambda"] {
        writeln!(out, "(declare-fun {name} () (_ BitVec {width}))").unwrap();
        writeln!(out, "(assert (bvult {name} {p_narrow}))").unwrap();
    }
    writeln!(out, "(assert (bvult x1 {bound}))").unwrap();
    writeln!(out, "(assert (bvult x2 {bound}))").unwrap();
    writeln!(
        out,
        "(assert (= (fp.mul y1 y1) (fp.add (fp.add (fp.mul (fp.mul x1 x1) x1) (fp.mul {a} x1)) {b})))"
    )
    .unwrap();
    writeln!(
        out,
        "(assert (= (fp.mul y2 y2) (fp.add (fp.add (fp.mul (fp.mul x2 x2) x2) (fp.mul {a} x2)) {b})))"
    )
    .unwrap();
    writeln!(out, "(assert (or").unwrap();
    writeln!(out, "  (and (bvult x1 x2)").unwrap();
    writeln!(
        out,
        "       (= (fp.mul lambda (fp.sub x2 x1)) (fp.sub y2 y1))"
    )
    .unwrap();
    writeln!(
        out,
        "       (= {target_x} (fp.sub (fp.sub (fp.mul lambda lambda) x1) x2))"
    )
    .unwrap();
    writeln!(
        out,
        "       (= {target_y} (fp.sub (fp.mul lambda (fp.sub x1 {target_x})) y1)))"
    )
    .unwrap();
    writeln!(out, "  (and (= x1 x2) (= y1 y2) (not (= y1 {zero}))").unwrap();
    writeln!(
        out,
        "       (= (fp.mul (fp.add y1 y1) lambda) (fp.add (fp.mul {three} (fp.mul x1 x1)) {a}))"
    )
    .unwrap();
    writeln!(
        out,
        "       (= {target_x} (fp.sub (fp.mul lambda lambda) (fp.add x1 x1)))"
    )
    .unwrap();
    writeln!(
        out,
        "       (= {target_y} (fp.sub (fp.mul lambda (fp.sub x1 {target_x})) y1)))))"
    )
    .unwrap();
    writeln!(out, "(check-sat)").unwrap();
    writeln!(out, "(get-value (x1 y1 x2 y2 lambda))").unwrap();

    let canonical = serde_json::to_vec(instance).map_err(|error| error.to_string())?;
    let receipt = SmtEncodingReceipt {
        schema: "prime-field-smt.encoding/v1".to_string(),
        logic: "QF_BV".to_string(),
        relation: "finite P1 + P2 = target; x1,x2 < x_bound".to_string(),
        field_bits: width,
        multiplication_bits: product,
        declared_variables: 5,
        field_helpers: 3,
        addition_branches: 2,
        assertion_count: 10,
        canonical_range_constraints: true,
        widened_arithmetic: true,
        symmetry_breaking: "x1 < x2 for ordinary addition; x1 = x2 for doubling".to_string(),
        canonical_instance_blake3: blake3_hex(&canonical),
        query_bytes: out.len() as u64,
        query_blake3: blake3_hex(out.as_bytes()),
    };
    Ok(SmtQuery { text: out, receipt })
}

/// Integrity and shape metadata for the native finite-field query.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FiniteFieldEncodingReceipt {
    pub schema: String,
    pub logic: String,
    pub relation: String,
    pub field_modulus: ExactUint,
    pub field_bits: u64,
    pub abscissa_bits: u64,
    pub declared_field_variables: u64,
    pub boolean_field_variables: u64,
    pub addition_branches: u32,
    pub assertion_count: u64,
    pub exact_field_arithmetic: bool,
    pub injective_bit_decomposition: bool,
    pub symmetry_breaking: String,
    pub recommended_solver_option: String,
    pub canonical_instance_blake3: String,
    pub query_bytes: u64,
    pub query_blake3: String,
}

/// Native finite-field SMT-LIB text and its receipt.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FiniteFieldQuery {
    pub text: String,
    pub receipt: FiniteFieldEncodingReceipt,
}

fn ff(value: &BigUint, p: &BigUint) -> String {
    format!("(as ff{} F)", (value % p).to_str_radix(10))
}

fn bit_is(name: &str, one: bool, zero: &str, one_value: &str) -> String {
    format!("(= {name} {})", if one { one_value } else { zero })
}

fn strict_bit_order(left: &[String], right: &[String], zero: &str, one: &str) -> String {
    let mut cases = Vec::new();
    for at in (0..left.len()).rev() {
        let mut terms = Vec::new();
        for higher in ((at + 1)..left.len()).rev() {
            terms.push(format!("(= {} {})", left[higher], right[higher]));
        }
        terms.push(bit_is(&left[at], false, zero, one));
        terms.push(bit_is(&right[at], true, zero, one));
        cases.push(format!("(and {})", terms.join(" ")));
    }
    format!("(or {})", cases.join(" "))
}

fn bits_below_bound(bits: &[String], bound: &BigUint, zero: &str, one: &str) -> String {
    if bound.is_one() {
        return bit_is(&bits[0], false, zero, one);
    }
    let bound_minus_one = bound - BigUint::one();
    if (bound & &bound_minus_one).is_zero() {
        return "true".to_string();
    }
    let mut cases = Vec::new();
    for at in (0..bits.len()).rev() {
        if !bound.bit(at as u64) {
            continue;
        }
        let mut terms = Vec::new();
        for higher in ((at + 1)..bits.len()).rev() {
            terms.push(bit_is(&bits[higher], bound.bit(higher as u64), zero, one));
        }
        terms.push(bit_is(&bits[at], false, zero, one));
        cases.push(format!("(and {})", terms.join(" ")));
    }
    format!("(or {})", cases.join(" "))
}

fn bitsum(bits: &[String]) -> String {
    if bits.len() == 1 {
        bits[0].clone()
    } else {
        format!("(ff.bitsum {})", bits.join(" "))
    }
}

/// Emit an exact cvc5 `QF_FF` query.
///
/// Abscissae are represented by Boolean field variables. Their weighted
/// sums are injective because the bit predicates constrain the represented
/// integers below `x_bound <= p`.
pub fn emit_finite_field_smt2(
    instance: &PrimeFieldSmtInstance,
) -> Result<FiniteFieldQuery, String> {
    let validation = validate_instance(instance)?;
    let p = &instance.curve.p.0;
    let bound_minus_one = &instance.x_bound.0 - BigUint::one();
    let abscissa_bits = bound_minus_one.bits().max(1);
    let zero = ff(&BigUint::zero(), p);
    let one = ff(&BigUint::one(), p);
    let a = ff(&instance.curve.a.0, p);
    let b = ff(&instance.curve.b.0, p);
    let target_x = ff(&instance.target.x.0, p);
    let target_y = ff(&instance.target.y.0, p);
    let three = ff(&BigUint::from(3u8), p);
    let x1_bits: Vec<String> = (0..abscissa_bits)
        .map(|index| format!("x1b{index}"))
        .collect();
    let x2_bits: Vec<String> = (0..abscissa_bits)
        .map(|index| format!("x2b{index}"))
        .collect();
    let mut out = String::new();
    writeln!(out, "; generated by crypto::cryptanalysis::prime_field_smt").unwrap();
    writeln!(out, "; exact native finite-field two-point decomposition").unwrap();
    writeln!(out, "(set-logic QF_FF)").unwrap();
    writeln!(out, "(set-option :produce-models true)").unwrap();
    writeln!(out, "(define-sort F () (_ FiniteField {}))", p).unwrap();
    writeln!(
        out,
        "(define-fun fp.sub ((lhs F) (rhs F)) F (ff.add lhs (ff.neg rhs)))"
    )
    .unwrap();
    for name in ["x1", "y1", "x2", "y2", "lambda"] {
        writeln!(out, "(declare-const {name} F)").unwrap();
    }
    for name in x1_bits.iter().chain(&x2_bits) {
        writeln!(out, "(declare-const {name} F)").unwrap();
        writeln!(out, "(assert (or (= {name} {zero}) (= {name} {one})))").unwrap();
    }
    writeln!(out, "(assert (= x1 {}))", bitsum(&x1_bits)).unwrap();
    writeln!(out, "(assert (= x2 {}))", bitsum(&x2_bits)).unwrap();
    writeln!(
        out,
        "(assert {})",
        bits_below_bound(&x1_bits, &instance.x_bound.0, &zero, &one)
    )
    .unwrap();
    writeln!(
        out,
        "(assert {})",
        bits_below_bound(&x2_bits, &instance.x_bound.0, &zero, &one)
    )
    .unwrap();
    writeln!(
        out,
        "(assert (= (ff.mul y1 y1) (ff.add (ff.mul x1 x1 x1) (ff.mul {a} x1) {b})))"
    )
    .unwrap();
    writeln!(
        out,
        "(assert (= (ff.mul y2 y2) (ff.add (ff.mul x2 x2 x2) (ff.mul {a} x2) {b})))"
    )
    .unwrap();
    let ordered = strict_bit_order(&x1_bits, &x2_bits, &zero, &one);
    writeln!(out, "(assert (or").unwrap();
    writeln!(out, "  (and {ordered}").unwrap();
    writeln!(
        out,
        "       (= (ff.mul lambda (fp.sub x2 x1)) (fp.sub y2 y1))"
    )
    .unwrap();
    writeln!(
        out,
        "       (= {target_x} (fp.sub (fp.sub (ff.mul lambda lambda) x1) x2))"
    )
    .unwrap();
    writeln!(
        out,
        "       (= {target_y} (fp.sub (ff.mul lambda (fp.sub x1 {target_x})) y1)))"
    )
    .unwrap();
    writeln!(out, "  (and (= x1 x2) (= y1 y2) (not (= y1 {zero}))").unwrap();
    writeln!(
        out,
        "       (= (ff.mul (ff.add y1 y1) lambda) (ff.add (ff.mul {three} x1 x1) {a}))"
    )
    .unwrap();
    writeln!(
        out,
        "       (= {target_x} (fp.sub (ff.mul lambda lambda) (ff.add x1 x1)))"
    )
    .unwrap();
    writeln!(
        out,
        "       (= {target_y} (fp.sub (ff.mul lambda (fp.sub x1 {target_x})) y1)))))"
    )
    .unwrap();
    writeln!(out, "(check-sat)").unwrap();
    writeln!(out, "(get-value (x1 y1 x2 y2 lambda))").unwrap();

    let canonical = serde_json::to_vec(instance).map_err(|error| error.to_string())?;
    let boolean_field_variables = abscissa_bits * 2;
    let receipt = FiniteFieldEncodingReceipt {
        schema: "prime-field-smt.ff-encoding/v1".to_string(),
        logic: "QF_FF".to_string(),
        relation: "finite P1 + P2 = target; x1,x2 < x_bound".to_string(),
        field_modulus: instance.curve.p.clone(),
        field_bits: validation.field_bits,
        abscissa_bits,
        declared_field_variables: 5 + boolean_field_variables,
        boolean_field_variables,
        addition_branches: 2,
        assertion_count: boolean_field_variables + 7,
        exact_field_arithmetic: true,
        injective_bit_decomposition: true,
        symmetry_breaking: "exact bitwise x1 < x2 for ordinary addition".to_string(),
        recommended_solver_option: "--ff-solver=split".to_string(),
        canonical_instance_blake3: blake3_hex(&canonical),
        query_bytes: out.len() as u64,
        query_blake3: blake3_hex(out.as_bytes()),
    };
    Ok(FiniteFieldQuery { text: out, receipt })
}

/// Metadata for the native finite-field Semaev-S3 query.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SemaevS3EncodingReceipt {
    pub schema: String,
    pub logic: String,
    pub relation: String,
    pub field_modulus: ExactUint,
    pub field_bits: u64,
    pub abscissa_bits: u64,
    pub x_values_scanned: u64,
    pub liftable_abscissae: u64,
    pub liftable_abscissae_blake3: String,
    pub membership_trie_nodes: u64,
    pub declared_field_variables: u64,
    pub boolean_field_variables: u64,
    pub exact_field_arithmetic: bool,
    pub complete_native_sign_lift_required: bool,
    pub symmetry_breaking: String,
    pub recommended_solver_option: String,
    pub canonical_instance_blake3: String,
    pub query_bytes: u64,
    pub query_blake3: String,
}

/// Native finite-field S3 SMT-LIB text and its receipt.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SemaevS3Query {
    pub text: String,
    pub receipt: SemaevS3EncodingReceipt,
}

fn liftable_abscissae(instance: &PrimeFieldSmtInstance) -> Result<Vec<BigUint>, String> {
    let bound = instance
        .x_bound
        .0
        .to_u64()
        .ok_or_else(|| "S3 factor-base enumeration requires x_bound <= u64::MAX".to_string())?;
    if bound > 1_000_000 {
        return Err(format!(
            "S3 factor-base enumeration guard exceeded: {bound} > 1000000"
        ));
    }
    let p = &instance.curve.p.0;
    let mut allowed = Vec::new();
    for x in 0..bound {
        let x = BigUint::from(x);
        let x2 = mul_mod(&x, &x, p);
        let x3 = mul_mod(&x2, &x, p);
        let rhs = add_mod(
            &add_mod(&x3, &mul_mod(&instance.curve.a.0, &x, p), p),
            &instance.curve.b.0,
            p,
        );
        if sqrt_mod_prime(&rhs, p).is_some() {
            allowed.push(x);
        }
    }
    Ok(allowed)
}

fn membership_trie(
    bit_names: &[String],
    values: &[u64],
    at: isize,
    zero: &str,
    one: &str,
    nodes: &mut u64,
) -> String {
    *nodes += 1;
    if values.is_empty() {
        return "false".to_string();
    }
    if at < 0 {
        return "true".to_string();
    }
    let capacity = 1usize << (at as usize + 1);
    if values.len() == capacity {
        return "true".to_string();
    }
    let mask = 1u64 << at as usize;
    let split = values.partition_point(|value| value & mask == 0);
    let low = membership_trie(bit_names, &values[..split], at - 1, zero, one, nodes);
    let high = membership_trie(bit_names, &values[split..], at - 1, zero, one, nodes);
    match (low.as_str(), high.as_str()) {
        ("false", "false") => "false".to_string(),
        (_, "false") => format!("(and (= {} {zero}) {low})", bit_names[at as usize]),
        ("false", _) => format!("(and (= {} {one}) {high})", bit_names[at as usize]),
        ("true", "true") => "true".to_string(),
        _ => format!(
            "(or (and (= {} {zero}) {low}) (and (= {} {one}) {high}))",
            bit_names[at as usize], bit_names[at as usize]
        ),
    }
}

fn digest_abscissae(values: &[BigUint]) -> String {
    let mut hasher = blake3::Hasher::new();
    for value in values {
        hasher.update(value.to_str_radix(10).as_bytes());
        hasher.update(b"\n");
    }
    hasher.finalize().to_hex().to_string()
}

/// Emit an exact finite-field Semaev-S3 query with a liftable-x trie.
pub fn emit_finite_field_s3_smt2(
    instance: &PrimeFieldSmtInstance,
) -> Result<SemaevS3Query, String> {
    let validation = validate_instance(instance)?;
    let allowed = liftable_abscissae(instance)?;
    let allowed_u64: Vec<u64> = allowed
        .iter()
        .map(|value| value.to_u64().unwrap())
        .collect();
    let p = &instance.curve.p.0;
    let bound_minus_one = &instance.x_bound.0 - BigUint::one();
    let abscissa_bits = bound_minus_one.bits().max(1);
    let zero = ff(&BigUint::zero(), p);
    let one = ff(&BigUint::one(), p);
    let two = ff(&BigUint::from(2u8), p);
    let four = ff(&BigUint::from(4u8), p);
    let a = ff(&instance.curve.a.0, p);
    let b = ff(&instance.curve.b.0, p);
    let target_x = ff(&instance.target.x.0, p);
    let x1_bits: Vec<String> = (0..abscissa_bits)
        .map(|index| format!("x1b{index}"))
        .collect();
    let x2_bits: Vec<String> = (0..abscissa_bits)
        .map(|index| format!("x2b{index}"))
        .collect();
    let formal_bits: Vec<String> = (0..abscissa_bits)
        .map(|index| format!("b{index}"))
        .collect();
    let mut trie_nodes = 0u64;
    let trie = membership_trie(
        &formal_bits,
        &allowed_u64,
        abscissa_bits as isize - 1,
        &zero,
        &one,
        &mut trie_nodes,
    );

    let diff = "(fp.sub x1 x2)";
    let sum = "(ff.add x1 x2)";
    let product = "(ff.mul x1 x2)";
    let term1 = format!("(ff.mul {diff} {diff} {target_x} {target_x})");
    let inner = format!("(ff.add (ff.mul {sum} (ff.add {product} {a})) (ff.mul {two} {b}))");
    let term2 = format!("(ff.mul {two} {inner} {target_x})");
    let product_minus_a = format!("(fp.sub {product} {a})");
    let term3 = format!("(ff.mul {product_minus_a} {product_minus_a})");
    let term4 = format!("(ff.mul {four} {b} {sum})");

    let mut out = String::new();
    writeln!(out, "; generated by crypto::cryptanalysis::prime_field_smt").unwrap();
    writeln!(
        out,
        "; native finite-field Semaev S3 with exact liftable-x trie"
    )
    .unwrap();
    writeln!(out, "(set-logic QF_FF)").unwrap();
    writeln!(out, "(set-option :produce-models true)").unwrap();
    writeln!(out, "(define-sort F () (_ FiniteField {}))", p).unwrap();
    writeln!(
        out,
        "(define-fun fp.sub ((lhs F) (rhs F)) F (ff.add lhs (ff.neg rhs)))"
    )
    .unwrap();
    let params = formal_bits
        .iter()
        .map(|name| format!("({name} F)"))
        .collect::<Vec<_>>()
        .join(" ");
    writeln!(out, "(define-fun factor.x ({params}) Bool {trie})").unwrap();
    for name in ["x1", "x2"] {
        writeln!(out, "(declare-const {name} F)").unwrap();
    }
    for name in x1_bits.iter().chain(&x2_bits) {
        writeln!(out, "(declare-const {name} F)").unwrap();
        writeln!(out, "(assert (or (= {name} {zero}) (= {name} {one})))").unwrap();
    }
    writeln!(out, "(assert (= x1 {}))", bitsum(&x1_bits)).unwrap();
    writeln!(out, "(assert (= x2 {}))", bitsum(&x2_bits)).unwrap();
    writeln!(out, "(assert (factor.x {}))", x1_bits.join(" ")).unwrap();
    writeln!(out, "(assert (factor.x {}))", x2_bits.join(" ")).unwrap();
    writeln!(
        out,
        "(assert (or (= x1 x2) {}))",
        strict_bit_order(&x1_bits, &x2_bits, &zero, &one)
    )
    .unwrap();
    writeln!(
        out,
        "(assert (= (ff.add {term1} (ff.neg {term2}) {term3} (ff.neg {term4})) {zero}))"
    )
    .unwrap();
    writeln!(out, "(check-sat)").unwrap();
    writeln!(out, "(get-value (x1 x2))").unwrap();

    let canonical = serde_json::to_vec(instance).map_err(|error| error.to_string())?;
    let boolean_field_variables = abscissa_bits * 2;
    let receipt = SemaevS3EncodingReceipt {
        schema: "prime-field-smt.s3-encoding/v1".to_string(),
        logic: "QF_FF".to_string(),
        relation: "S3(x1,x2,x_target)=0 with exact liftable factor-base x".to_string(),
        field_modulus: instance.curve.p.clone(),
        field_bits: validation.field_bits,
        abscissa_bits,
        x_values_scanned: instance.x_bound.0.to_u64().unwrap(),
        liftable_abscissae: allowed.len() as u64,
        liftable_abscissae_blake3: digest_abscissae(&allowed),
        membership_trie_nodes: trie_nodes,
        declared_field_variables: 2 + boolean_field_variables,
        boolean_field_variables,
        exact_field_arithmetic: true,
        complete_native_sign_lift_required: true,
        symmetry_breaking: "exact bitwise x1 <= x2".to_string(),
        recommended_solver_option: "--ff-solver=split".to_string(),
        canonical_instance_blake3: blake3_hex(&canonical),
        query_bytes: out.len() as u64,
        query_blake3: blake3_hex(out.as_bytes()),
    };
    Ok(SemaevS3Query { text: out, receipt })
}

/// Lift an S3 pair to all signed `F_p` points and return the first full
/// decomposition that independently sums to the original target.
pub fn lift_s3_abscissae(
    instance: &PrimeFieldSmtInstance,
    x1: &ExactUint,
    x2: &ExactUint,
) -> Result<Option<DecompositionWitness>, String> {
    validate_instance(instance)?;
    if x1.0 >= instance.x_bound.0
        || x2.0 >= instance.x_bound.0
        || x1.0 > x2.0
        || x1.0 >= instance.curve.p.0
        || x2.0 >= instance.curve.p.0
    {
        return Ok(None);
    }
    let p = &instance.curve.p.0;
    let roots = |x: &ExactUint| -> Vec<ExactUint> {
        let x_squared = mul_mod(&x.0, &x.0, p);
        let rhs = add_mod(
            &add_mod(
                &mul_mod(&x_squared, &x.0, p),
                &mul_mod(&instance.curve.a.0, &x.0, p),
                p,
            ),
            &instance.curve.b.0,
            p,
        );
        let Some(root) = sqrt_mod_prime(&rhs, p) else {
            return Vec::new();
        };
        if root.is_zero() {
            vec![ExactUint(root)]
        } else {
            vec![ExactUint(root.clone()), ExactUint(p - root)]
        }
    };
    for y1 in roots(x1) {
        for y2 in roots(x2) {
            let p1 = AffinePoint {
                x: x1.clone(),
                y: y1.clone(),
            };
            let p2 = AffinePoint {
                x: x2.clone(),
                y: y2,
            };
            let (sum, lambda) = add_points(&instance.curve, &p1, &p2)?;
            if sum != GroupPoint::Affine(instance.target.clone()) {
                continue;
            }
            let Some(lambda) = lambda else {
                continue;
            };
            let witness = DecompositionWitness {
                x1: p1.x,
                y1: p1.y,
                x2: p2.x,
                y2: p2.y,
                lambda: ExactUint(lambda),
            };
            if verify_witness(instance, &witness)?.verified {
                return Ok(Some(witness));
            }
        }
    }
    Ok(None)
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum SExpr {
    Atom(String),
    List(Vec<SExpr>),
}

fn tokens(input: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut atom = String::new();
    let mut comment = false;
    for ch in input.chars() {
        if comment {
            if ch == '\n' {
                comment = false;
            }
            continue;
        }
        if ch == ';' {
            if !atom.is_empty() {
                tokens.push(std::mem::take(&mut atom));
            }
            comment = true;
        } else if ch == '(' || ch == ')' {
            if !atom.is_empty() {
                tokens.push(std::mem::take(&mut atom));
            }
            tokens.push(ch.to_string());
        } else if ch.is_whitespace() {
            if !atom.is_empty() {
                tokens.push(std::mem::take(&mut atom));
            }
        } else {
            atom.push(ch);
        }
    }
    if !atom.is_empty() {
        tokens.push(atom);
    }
    tokens
}

fn parse_expr(tokens: &[String], at: &mut usize) -> Result<SExpr, String> {
    let token = tokens
        .get(*at)
        .ok_or_else(|| "unexpected end of SMT output".to_string())?;
    *at += 1;
    if token == "(" {
        let mut values = Vec::new();
        while tokens.get(*at).is_some_and(|token| token != ")") {
            values.push(parse_expr(tokens, at)?);
        }
        if tokens.get(*at).is_none() {
            return Err("unterminated list in SMT output".to_string());
        }
        *at += 1;
        Ok(SExpr::List(values))
    } else if token == ")" {
        Err("unexpected ')' in SMT output".to_string())
    } else {
        Ok(SExpr::Atom(token.clone()))
    }
}

fn parse_all(input: &str) -> Result<Vec<SExpr>, String> {
    let tokens = tokens(input);
    let mut at = 0usize;
    let mut expressions = Vec::new();
    while at < tokens.len() {
        expressions.push(parse_expr(&tokens, &mut at)?);
    }
    Ok(expressions)
}

fn finite_field_sort_modulus(expression: &SExpr) -> Option<BigUint> {
    let SExpr::List(values) = expression else {
        return None;
    };
    if values.len() != 3 {
        return None;
    }
    match (&values[0], &values[1], &values[2]) {
        (SExpr::Atom(head), SExpr::Atom(sort), SExpr::Atom(modulus))
            if head == "_" && sort == "FiniteField" =>
        {
            BigUint::parse_bytes(modulus.as_bytes(), 10)
        }
        _ => None,
    }
}

fn finite_field_value(atom: &str, modulus: Option<&BigUint>) -> Option<BigUint> {
    let raw = atom.strip_prefix("ff")?;
    let (negative, magnitude) = raw
        .strip_prefix('-')
        .map_or((false, raw), |magnitude| (true, magnitude));
    let magnitude = BigUint::parse_bytes(magnitude.as_bytes(), 10)?;
    match modulus {
        Some(modulus) => {
            let reduced = magnitude % modulus;
            if negative && !reduced.is_zero() {
                Some(modulus - reduced)
            } else {
                Some(reduced)
            }
        }
        None if !negative => Some(magnitude),
        None => None,
    }
}

fn finite_field_literal(atom: &str, expected_modulus: Option<&BigUint>) -> Option<BigUint> {
    let body = atom.strip_prefix("#f")?;
    let (value, modulus) = body.rsplit_once('m')?;
    let modulus = BigUint::parse_bytes(modulus.as_bytes(), 10)?;
    if expected_modulus.is_some_and(|expected| expected != &modulus) {
        return None;
    }
    finite_field_value(&format!("ff{value}"), Some(&modulus))
}

fn model_value(expression: &SExpr, modulus: Option<&BigUint>) -> Option<BigUint> {
    match expression {
        SExpr::Atom(atom) if atom.starts_with("#f") => finite_field_literal(atom, modulus),
        SExpr::Atom(atom) if atom.starts_with("#b") => {
            BigUint::parse_bytes(&atom.as_bytes()[2..], 2)
        }
        SExpr::Atom(atom) if atom.starts_with("#x") => {
            BigUint::parse_bytes(&atom.as_bytes()[2..], 16)
        }
        SExpr::Atom(atom) => BigUint::parse_bytes(atom.as_bytes(), 10),
        SExpr::List(values) if values.len() == 3 => match (&values[0], &values[1], &values[2]) {
            (SExpr::Atom(head), SExpr::Atom(value), SExpr::Atom(width))
                if head == "_" && value.starts_with("bv") =>
            {
                let value = BigUint::parse_bytes(&value.as_bytes()[2..], 10)?;
                let width: u64 = width.parse().ok()?;
                (value.bits() <= width).then_some(value)
            }
            (SExpr::Atom(head), SExpr::Atom(value), sort) if head == "as" => {
                let embedded = finite_field_sort_modulus(sort);
                if embedded
                    .as_ref()
                    .zip(modulus)
                    .is_some_and(|(embedded, expected)| embedded != expected)
                {
                    return None;
                }
                finite_field_value(value, embedded.as_ref().or(modulus))
            }
            _ => None,
        },
        _ => None,
    }
}

fn assignments(
    expression: &SExpr,
    modulus: Option<&BigUint>,
    out: &mut HashMap<String, BigUint>,
) -> Result<(), String> {
    if let SExpr::List(values) = expression {
        if values.len() == 2 {
            if let SExpr::Atom(name) = &values[0] {
                if ["x1", "y1", "x2", "y2", "lambda"].contains(&name.as_str()) {
                    let value = model_value(&values[1], modulus)
                        .ok_or_else(|| format!("unsupported value for {name}"))?;
                    if let Some(previous) = out.insert(name.clone(), value.clone()) {
                        if previous != value {
                            return Err(format!("conflicting assignments for {name}"));
                        }
                    }
                    return Ok(());
                }
            }
        }
        for value in values {
            assignments(value, modulus, out)?;
        }
    }
    Ok(())
}

/// Parsed `check-sat` result. UNSAT is a solver assertion unless separately
/// checked with [`exhaustive_reference`]; SAT carries a replayable witness.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum ParsedSolverOutput {
    Sat { witness: DecompositionWitness },
    Unsat,
    Unknown,
}

fn parse_solver_output_inner(
    output: &str,
    modulus: Option<&BigUint>,
) -> Result<ParsedSolverOutput, String> {
    let expressions = parse_all(output)?;
    let status = expressions.iter().find_map(|expression| match expression {
        SExpr::Atom(atom) if matches!(atom.as_str(), "sat" | "unsat" | "unknown") => {
            Some(atom.as_str())
        }
        _ => None,
    });
    match status {
        Some("unsat") => return Ok(ParsedSolverOutput::Unsat),
        Some("unknown") => return Ok(ParsedSolverOutput::Unknown),
        Some("sat") => {}
        _ => return Err("SMT output has no sat/unsat/unknown status".to_string()),
    }
    let mut values = HashMap::new();
    for expression in &expressions {
        assignments(expression, modulus, &mut values)?;
    }
    let take = |name: &str| {
        values
            .get(name)
            .cloned()
            .map(ExactUint)
            .ok_or_else(|| format!("SAT output omitted {name}"))
    };
    Ok(ParsedSolverOutput::Sat {
        witness: DecompositionWitness {
            x1: take("x1")?,
            y1: take("y1")?,
            x2: take("x2")?,
            y2: take("y2")?,
            lambda: take("lambda")?,
        },
    })
}

/// Parse cvc5/Z3-style bit-vector `get-value` output without trusting order.
pub fn parse_solver_output(output: &str) -> Result<ParsedSolverOutput, String> {
    parse_solver_output_inner(output, None)
}

/// Parse a solver model with the instance modulus available for signed cvc5
/// finite-field values such as `(as ff-1 F)`.
pub fn parse_solver_output_for_instance(
    output: &str,
    instance: &PrimeFieldSmtInstance,
) -> Result<ParsedSolverOutput, String> {
    parse_solver_output_inner(output, Some(&instance.curve.p.0))
}

/// Parsed result of an S3 query, before native sign lifting.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum ParsedAbscissaOutput {
    Sat { x1: ExactUint, x2: ExactUint },
    Unsat,
    Unknown,
}

/// Parse the two abscissae returned by the S3 backend.
pub fn parse_abscissa_output_for_instance(
    output: &str,
    instance: &PrimeFieldSmtInstance,
) -> Result<ParsedAbscissaOutput, String> {
    let expressions = parse_all(output)?;
    let status = expressions.iter().find_map(|expression| match expression {
        SExpr::Atom(atom) if matches!(atom.as_str(), "sat" | "unsat" | "unknown") => {
            Some(atom.as_str())
        }
        _ => None,
    });
    match status {
        Some("unsat") => return Ok(ParsedAbscissaOutput::Unsat),
        Some("unknown") => return Ok(ParsedAbscissaOutput::Unknown),
        Some("sat") => {}
        _ => return Err("SMT output has no sat/unsat/unknown status".to_string()),
    }
    let mut values = HashMap::new();
    for expression in &expressions {
        assignments(expression, Some(&instance.curve.p.0), &mut values)?;
    }
    Ok(ParsedAbscissaOutput::Sat {
        x1: ExactUint(
            values
                .remove("x1")
                .ok_or_else(|| "SAT output omitted x1".to_string())?,
        ),
        x2: ExactUint(
            values
                .remove("x2")
                .ok_or_else(|| "SAT output omitted x2".to_string())?,
        ),
    })
}

/// Raw result of invoking the pinned cvc5 command-line interface.
#[derive(Clone, Debug)]
pub struct SolverExecution {
    pub arguments: Vec<String>,
    pub elapsed: Duration,
    pub timed_out: bool,
    pub success: bool,
    pub exit_code: Option<i32>,
    pub stdout: Vec<u8>,
    pub stderr: Vec<u8>,
}

/// Run cvc5 with both its own per-query limit and a five-second harness
/// grace period. The caller retains stdout/stderr even on timeout.
pub fn run_cvc5(
    executable: &Path,
    query: &Path,
    timeout: Duration,
) -> Result<SolverExecution, String> {
    run_cvc5_with_options(executable, query, timeout, &[])
}

/// Run cvc5 with explicit, receipt-visible solver options.
pub fn run_cvc5_with_options(
    executable: &Path,
    query: &Path,
    timeout: Duration,
    solver_options: &[String],
) -> Result<SolverExecution, String> {
    let timeout_ms = timeout.as_millis();
    let mut arguments = vec![
        "--lang=smt2".to_string(),
        format!("--tlimit-per={timeout_ms}"),
    ];
    arguments.extend(solver_options.iter().cloned());
    arguments.push(query.display().to_string());
    let started = Instant::now();
    let mut child = Command::new(executable)
        .args(&arguments)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|error| format!("spawn {}: {error}", executable.display()))?;
    let deadline = started + timeout + Duration::from_secs(5);
    let mut timed_out = false;
    loop {
        match child.try_wait() {
            Ok(Some(_)) => break,
            Ok(None) if Instant::now() >= deadline => {
                timed_out = true;
                let _ = child.kill();
                break;
            }
            Ok(None) => std::thread::sleep(Duration::from_millis(5)),
            Err(error) => return Err(format!("wait for cvc5: {error}")),
        }
    }
    let output = child
        .wait_with_output()
        .map_err(|error| format!("collect cvc5 output: {error}"))?;
    Ok(SolverExecution {
        arguments,
        elapsed: started.elapsed(),
        timed_out,
        success: output.status.success() && !timed_out,
        exit_code: output.status.code(),
        stdout: output.stdout,
        stderr: output.stderr,
    })
}

/// Stream a file through BLAKE3 without loading a solver binary into memory.
pub fn blake3_file(path: &Path) -> io::Result<String> {
    let mut file = File::open(path)?;
    let mut hasher = blake3::Hasher::new();
    let mut buffer = [0u8; 64 * 1024];
    loop {
        let count = file.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        hasher.update(&buffer[..count]);
    }
    Ok(hasher.finalize().to_hex().to_string())
}

/// Verdict from the native factor-set reference oracle.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ExhaustiveStatus {
    Sat,
    Unsat,
}

/// Receipt for the exact `R-P` factor-set membership reference.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExhaustiveReceipt {
    pub schema: String,
    pub field_bits: u64,
    pub x_values_enumerated: u64,
    pub factor_points: u64,
    pub candidate_probes: u64,
    pub status: ExhaustiveStatus,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub witness: Option<DecompositionWitness>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub verification: Option<VerificationReceipt>,
}

fn ordered_pair(mut left: AffinePoint, mut right: AffinePoint) -> (AffinePoint, AffinePoint) {
    if (&left.x, &left.y) > (&right.x, &right.y) {
        std::mem::swap(&mut left, &mut right);
    }
    (left, right)
}

/// Exact small-instance reference: enumerate every signed curve point with
/// `x < x_bound`, then test `target - P` membership for every `P`.
///
/// This oracle is intentionally limited to a `u64` modulus. It is a control
/// for the SMT implementation, not the proposed scaling path.
pub fn exhaustive_reference(
    instance: &PrimeFieldSmtInstance,
    max_factor_points: usize,
) -> Result<ExhaustiveReceipt, String> {
    let validation = validate_instance(instance)?;
    let p = instance
        .curve
        .p
        .0
        .to_u64()
        .ok_or_else(|| "exhaustive reference requires p <= u64::MAX".to_string())?;
    let bound = instance
        .x_bound
        .0
        .to_u64()
        .ok_or_else(|| "exhaustive reference requires x_bound <= u64::MAX".to_string())?;
    let a = instance.curve.a.0.to_u64().unwrap();
    let b = instance.curve.b.0.to_u64().unwrap();
    let modulus = p as u128;
    let mul = |left: u64, right: u64| -> u64 { ((left as u128 * right as u128) % modulus) as u64 };
    let mut points = Vec::new();
    for x in 0..bound {
        let x2 = mul(x, x);
        let x3 = mul(x2, x);
        let rhs = ((x3 as u128 + mul(a, x) as u128 + b as u128) % modulus) as u64;
        if let Some(y) = super::residual_walk::sqrt_mod(rhs, p) {
            points.push(AffinePoint {
                x: x.into(),
                y: y.into(),
            });
            if y != 0 {
                points.push(AffinePoint {
                    x: x.into(),
                    y: (p - y).into(),
                });
            }
        }
        if points.len() > max_factor_points {
            return Err(format!(
                "factor-point guard exceeded: {} > {max_factor_points}",
                points.len()
            ));
        }
    }
    points.sort_by(|left, right| (&left.x, &left.y).cmp(&(&right.x, &right.y)));
    let set: HashSet<AffinePoint> = points.iter().cloned().collect();
    let target = &instance.target;
    let mut probes = 0u64;
    for point in &points {
        probes += 1;
        let negated = AffinePoint {
            x: point.x.clone(),
            y: ExactUint(if point.y.0.is_zero() {
                BigUint::zero()
            } else {
                &instance.curve.p.0 - &point.y.0
            }),
        };
        let (needed, _) = add_points(&instance.curve, target, &negated)?;
        let GroupPoint::Affine(needed) = needed else {
            continue;
        };
        if needed.x.0 >= instance.x_bound.0 || !set.contains(&needed) {
            continue;
        }
        let (p1, p2) = ordered_pair(point.clone(), needed);
        let (_, lambda) = add_points(&instance.curve, &p1, &p2)?;
        let Some(lambda) = lambda else {
            continue;
        };
        let witness = DecompositionWitness {
            x1: p1.x,
            y1: p1.y,
            x2: p2.x,
            y2: p2.y,
            lambda: ExactUint(lambda),
        };
        let verification = verify_witness(instance, &witness)?;
        if !verification.verified {
            return Err("native reference constructed an invalid witness".to_string());
        }
        return Ok(ExhaustiveReceipt {
            schema: "prime-field-smt.exhaustive/v1".to_string(),
            field_bits: validation.field_bits,
            x_values_enumerated: bound,
            factor_points: points.len() as u64,
            candidate_probes: probes,
            status: ExhaustiveStatus::Sat,
            witness: Some(witness),
            verification: Some(verification),
        });
    }
    Ok(ExhaustiveReceipt {
        schema: "prime-field-smt.exhaustive/v1".to_string(),
        field_bits: validation.field_bits,
        x_values_enumerated: bound,
        factor_points: points.len() as u64,
        candidate_probes: probes,
        status: ExhaustiveStatus::Unsat,
        witness: None,
        verification: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(x: u64, y: u64) -> AffinePoint {
        AffinePoint {
            x: x.into(),
            y: y.into(),
        }
    }

    fn toy(target: AffinePoint, bound: u64) -> PrimeFieldSmtInstance {
        PrimeFieldSmtInstance {
            schema: INSTANCE_SCHEMA.to_string(),
            curve: PrimeFieldCurve {
                name: None,
                p: 17.into(),
                a: 2.into(),
                b: 2.into(),
            },
            target,
            x_bound: bound.into(),
            source: None,
        }
    }

    fn toy_witness() -> DecompositionWitness {
        DecompositionWitness {
            x1: 0.into(),
            y1: 6.into(),
            x2: 3.into(),
            y2: 1.into(),
            lambda: 4.into(),
        }
    }

    #[test]
    fn exact_uint_accepts_all_supported_forms() {
        for encoded in [r#""17""#, r#""0x11""#, r#""0b10001""#, "17"] {
            let value: ExactUint = serde_json::from_str(encoded).unwrap();
            assert_eq!(value, ExactUint::from(17));
            assert_eq!(serde_json::to_string(&value).unwrap(), r#""17""#);
        }
    }

    #[test]
    fn planted_witness_verifies_and_corruption_does_not() {
        let instance = toy(point(13, 10), 4);
        let receipt = verify_witness(&instance, &toy_witness()).unwrap();
        assert!(receipt.verified);
        assert_eq!(receipt.branch.as_deref(), Some("ordinary"));

        let mut corrupt = toy_witness();
        corrupt.lambda = 5.into();
        assert!(!verify_witness(&instance, &corrupt).unwrap().verified);
    }

    #[test]
    fn exporter_uses_widened_bitvectors_and_both_branches() {
        let query = emit_smt2(&toy(point(13, 10), 4)).unwrap();
        assert_eq!(query.receipt.field_bits, 5);
        assert_eq!(query.receipt.multiplication_bits, 10);
        assert!(query.text.contains("((_ zero_extend 5) lhs)"));
        assert!(query.text.contains("(_ bv17 10)"));
        assert!(query.text.contains("(bvult x1 x2)"));
        assert!(query.text.contains("(= x1 x2) (= y1 y2)"));
        assert_eq!(query.receipt.query_bytes, query.text.len() as u64);
        assert_eq!(
            query.receipt.query_blake3,
            blake3_hex(query.text.as_bytes())
        );
    }

    #[test]
    fn parser_accepts_binary_hex_and_indexed_decimal_values() {
        let output =
            "sat\n((x1 #b00000) (y1 #x06) (x2 (_ bv3 5)) (y2 #b00001) (lambda (_ bv4 5)))\n";
        assert_eq!(
            parse_solver_output(output).unwrap(),
            ParsedSolverOutput::Sat {
                witness: toy_witness()
            }
        );
        assert_eq!(
            parse_solver_output("unsat\n").unwrap(),
            ParsedSolverOutput::Unsat
        );
        assert_eq!(
            parse_solver_output("unknown\n").unwrap(),
            ParsedSolverOutput::Unknown
        );
    }

    #[test]
    fn finite_field_export_and_signed_model_are_exact() {
        let instance = toy(point(13, 10), 4);
        let query = emit_finite_field_smt2(&instance).unwrap();
        assert_eq!(query.receipt.logic, "QF_FF");
        assert_eq!(query.receipt.abscissa_bits, 2);
        assert_eq!(query.receipt.boolean_field_variables, 4);
        assert!(query.text.contains("(ff.bitsum x1b0 x1b1)"));
        assert!(query.text.contains("(_ FiniteField 17)"));
        assert!(query.text.contains("(or (and"));

        let output = "sat\n((x1 (as ff0 F)) (y1 (as ff6 F)) (x2 (as ff3 F)) (y2 (as ff-16 F)) (lambda (as ff4 F)))\n";
        assert_eq!(
            parse_solver_output_for_instance(output, &instance).unwrap(),
            ParsedSolverOutput::Sat {
                witness: toy_witness()
            }
        );
        let literal_output =
            "sat\n((x1 #f0m17) (y1 #f6m17) (x2 #f3m17) (y2 #f1m17) (lambda #f4m17))\n";
        assert_eq!(
            parse_solver_output_for_instance(literal_output, &instance).unwrap(),
            ParsedSolverOutput::Sat {
                witness: toy_witness()
            }
        );
        let wrong_modulus =
            "sat\n((x1 #f0m19) (y1 #f6m19) (x2 #f3m19) (y2 #f1m19) (lambda #f4m19))\n";
        assert!(parse_solver_output_for_instance(wrong_modulus, &instance).is_err());

        let one_bit = emit_finite_field_smt2(&toy(point(3, 1), 1)).unwrap();
        assert!(one_bit.text.contains("(assert (= x1 x1b0))"));
        assert!(!one_bit.text.contains("(ff.bitsum x1b0)"));
    }

    #[test]
    fn semaev_s3_backend_restricts_and_lifts_abscissae() {
        let instance = toy(point(13, 10), 4);
        let query = emit_finite_field_s3_smt2(&instance).unwrap();
        assert_eq!(query.receipt.liftable_abscissae, 2);
        assert_eq!(query.receipt.x_values_scanned, 4);
        assert!(query.text.contains("(define-fun factor.x"));
        assert!(query.text.contains("Semaev S3"));
        assert!(!query.text.contains("(declare-const y1"));

        let output = "sat\n((x1 #f0m17) (x2 #f3m17))\n";
        let ParsedAbscissaOutput::Sat { x1, x2 } =
            parse_abscissa_output_for_instance(output, &instance).unwrap()
        else {
            panic!("expected SAT abscissae")
        };
        let witness = lift_s3_abscissae(&instance, &x1, &x2)
            .unwrap()
            .expect("S3 pair lifts");
        assert!(verify_witness(&instance, &witness).unwrap().verified);
    }

    #[test]
    fn exact_reference_proves_sat_and_unsat_controls() {
        let sat = exhaustive_reference(&toy(point(13, 10), 4), 100).unwrap();
        assert_eq!(sat.status, ExhaustiveStatus::Sat);
        assert!(sat.verification.unwrap().verified);

        let unsat = exhaustive_reference(&toy(point(3, 1), 1), 100).unwrap();
        assert_eq!(unsat.status, ExhaustiveStatus::Unsat);
        assert_eq!(unsat.factor_points, 2);
        assert_eq!(unsat.candidate_probes, 2);
    }

    #[test]
    fn validation_rejects_composite_field_and_off_curve_target() {
        let mut composite = toy(point(13, 10), 4);
        composite.curve.p = 15.into();
        assert!(validate_instance(&composite).is_err());

        let off_curve = toy(point(13, 11), 4);
        assert!(validate_instance(&off_curve).is_err());
    }
}
