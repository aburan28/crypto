//! Exact, replayable audit of sect113r1 and two distinct weakness surfaces.
//!
//! This module deliberately keeps three claims separate:
//!
//! 1. `sect113r1` is **class-weak** even against a conservative frozen 80-bit
//!    target because its prime subgroup admits generic rho in about 56 bits
//!    of work.  Every
//!    curve in the same `F_{2^113}`-isogeny class has the same subgroup order.
//! 2. The library's unchecked low-level multiplier is **implementation-weak**
//!    when it is exposed to attacker-selected points.  Its formula omits `b`,
//!    so it also operates on the singular same-`a`, `b = 0` companion.
//! 3. Explicit degree-5 isogenies prove DLP transport to valid neighbours;
//!    they do not reduce the generic attack cost and are not counted as the
//!    cause of either weakness above.
//!
//! Algebraic statements are checked exactly.  The small planted recovery is
//! intentionally bounded: it recovers a fixed 52.6-bit demonstration scalar
//! using four Pohlig--Hellman components.  The remaining 58.8-bit component
//! for a full-width scalar is reported, not silently treated as executed.

use crate::binary_ecc::curve::{point_add, point_neg, scalar_mul};
use crate::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, F2mPoly};
use crate::cryptanalysis::binary_velu::{
    division_polynomial, velu_codomain, velu_point_map_binary,
};
use crate::hash::sha256;
use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

pub const REPORT_SCHEMA: &str = "crypto.sect113r1-audit/v2";
pub const SOURCE_URL: &str = "https://www.secg.org/SEC2-Ver-1.0.pdf";
pub const SOURCE_SHA256: &str = "d1b16728ad83888fd656d16b99dc71bcd5541d42d848ffd0de7c62c19010d8c3";
pub const CURVE_ICV1: &str = "ICV1:f2m-113-99967757:-122610772499221213:10384593717069655379671765157661406:0x6942e38fc45c62366c09aa8204cd:unk:unk:r:97df4ac684cb";
pub const CURVE_ID: &str = "EC1N113Csect113r1hf529f17bd191";
pub const CURVE_UID: &str =
    "urn:ec-record:1:sha256:f529f17bd1913792333a661e3557ad6b8e0ca2d4d02b939bc069d17d9fd94d97";

const FULL_SEED: &[u8] = b"SECT113R1-SINGULAR-V1";
const FULL_DOMAIN: &[u8] = b"crypto/sect113r1-singular-recovery/seed/v1\0";
const FULL_DIGEST: &str = "c4cfd5aa9d48d2079788cbdd3c6c7f579cf7d41cd002106eb2bab540ad639f08";
const FULL_SCALAR: &str = "1108386698526968361534713354743651";
const DEMO_DOMAIN: &[u8] = b"crypto/sect113r1-singular/demo/v1\0";
const DEMO_DIGEST: &str = "fa4be280c9e9f1e2112ae090729ed3b4faf78e0c1133db4a71836e266896b463";
const DEMO_SCALAR: u64 = 6_599_291_615_786_184;
const SINGULAR_FACTORS: [u64; 5] = [3, 227, 48_817, 636_190_001, 491_003_369_344_660_409];
const EXECUTED_FACTORS: [u64; 4] = [3, 227, 48_817, 636_190_001];
const BINARY_FIELD_ENCODING: &str =
    "hex polynomial coefficient bitset, least significant bit is constant";

fn decimal(value: &str) -> BigUint {
    BigUint::parse_bytes(value.as_bytes(), 10).expect("valid frozen decimal")
}

fn point_hex(point: &BinaryPoint) -> String {
    match point {
        BinaryPoint::Infinity => "00".into(),
        BinaryPoint::Affine { x, y } => format!("01{:030x}{:030x}", x.to_biguint(), y.to_biguint()),
    }
}

fn field_hex(value: &F2mElement) -> String {
    format!("{:030x}", value.to_biguint())
}

fn fixed_field_bytes(value: &F2mElement) -> [u8; 15] {
    let source = value.to_biguint().to_bytes_be();
    let mut out = [0u8; 15];
    out[15 - source.len()..].copy_from_slice(&source);
    out
}

fn sha256_hex(bytes: &[u8]) -> String {
    hex::encode(sha256(bytes))
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct CurveIdentityCertificate {
    pub icv1: String,
    pub icv1_slug: String,
    pub model_sha256: String,
    pub curve_id: String,
    pub curve_uid: String,
    pub curve_sha256: String,
    pub field_sha256: String,
}

/// Recompute both repository curve identities from the exact binary model and
/// representation.  The byte strings below mirror the canonical, sorted-key,
/// compact JSON contracts in `scripts/curve_id.py` and
/// `tools/curve_identity.py`; tests pin the two independently generated
/// endpoint vectors so a serialization drift fails closed.
fn binary_curve_identity(
    curve: &BinaryCurve,
    curve_tag: &str,
) -> Result<CurveIdentityCertificate, String> {
    if curve.b.is_zero() {
        return Err("a singular binary model has no curve identity".into());
    }
    if !curve_tag
        .bytes()
        .enumerate()
        .all(|(index, byte)| byte.is_ascii_lowercase() || (index > 0 && byte.is_ascii_digit()))
    {
        return Err("curve identity tag must match [a-z][a-z0-9]*".into());
    }
    let BinaryPoint::Affine { x: gx, y: gy } = &curve.generator else {
        return Err("curve identity requires an affine subgroup generator".into());
    };
    if gx.m_value() != curve.m || gy.m_value() != curve.m {
        return Err("curve identity generator has the wrong field width".into());
    }

    let mut modulus = BigUint::one() << curve.m as usize;
    for exponent in &curve.irreducible.low_terms {
        modulus |= BigUint::one() << *exponent as usize;
    }
    let modulus_hex = format!("0x{modulus:x}");
    let field_contract = format!("f2m-modulus:{modulus_hex}");
    let field = format!(
        "f2m-{}-{}",
        curve.m,
        &sha256_hex(field_contract.as_bytes())[..8]
    );

    let a = curve.a.to_biguint();
    let b = curve.b.to_biguint();
    let a_hex = format!("0x{a:x}");
    let b_hex = format!("0x{b:x}");
    let model_json = format!(
        "{{\"a\":\"{a_hex}\",\"b\":\"{b_hex}\",\"field\":\"{field}\",\"form\":\"y^2+xy=x^3+a*x^2+b\",\"modulus\":\"{modulus_hex}\",\"v\":\"1\"}}"
    );
    let model_sha256 = sha256_hex(model_json.as_bytes());

    let q = BigUint::one() << curve.m as usize;
    let group_order = &curve.order * &curve.cofactor;
    let trace = BigInt::from_biguint(Sign::Plus, &q + BigUint::one())
        - BigInt::from_biguint(Sign::Plus, group_order.clone());
    if &trace * &trace > BigInt::from_biguint(Sign::Plus, &q * 4u8) {
        return Err("curve identity group order violates Hasse".into());
    }
    let j = curve
        .b
        .flt_inverse(&curve.irreducible)
        .ok_or("nonsingular binary model has no inverse for b")?;
    let j_hex = format!("0x{:x}", j.to_biguint());
    let icv1 = format!(
        "ICV1:{field}:{trace}:{group_order}:{j_hex}:unk:unk:r:{}",
        &model_sha256[..12]
    );
    let trace_slug = if trace.sign() == Sign::Minus {
        format!("tm{}", trace.magnitude())
    } else {
        format!("t{trace}")
    };
    let icv1_slug = format!("icv1-f2m{}-{trace_slug}-{}", curve.m, &model_sha256[..8]);

    let mut modulus_exponents = Vec::with_capacity(curve.irreducible.low_terms.len() + 1);
    modulus_exponents.push(curve.m);
    modulus_exponents.extend(curve.irreducible.low_terms.iter().copied());
    modulus_exponents.sort_unstable_by(|left, right| right.cmp(left));
    let modulus_exponents_json = modulus_exponents
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join(",");
    let field_json = format!(
        "{{\"characteristic\":2,\"degree\":{},\"element_encoding\":\"{}\",\"modulus_exponents\":[{}],\"representation\":\"polynomial\"}}",
        curve.m, BINARY_FIELD_ENCODING, modulus_exponents_json
    );
    let curve_json = format!(
        "{{\"coefficients\":[1,{a},0,0,{b}],\"cofactor\":{},\"generator\":[\"0x{:x}\",\"0x{:x}\"],\"model\":\"binary Weierstrass\",\"subgroup_order\":\"{}\",\"target_group\":\"prime-order subgroup\"}}",
        curve.cofactor,
        gx.to_biguint(),
        gy.to_biguint(),
        curve.order
    );
    let curve_record_json = format!("{{\"curve\":{curve_json},\"field\":{field_json}}}");
    let curve_sha256 = sha256_hex(curve_record_json.as_bytes());
    let field_sha256 = sha256_hex(field_json.as_bytes());
    let curve_id = format!("EC1N{}C{curve_tag}h{}", curve.m, &curve_sha256[..12]);
    let curve_uid = format!("urn:ec-record:1:sha256:{curve_sha256}");

    Ok(CurveIdentityCertificate {
        icv1,
        icv1_slug,
        model_sha256,
        curve_id,
        curve_uid,
        curve_sha256,
        field_sha256,
    })
}

fn gcd_big(mut left: BigUint, mut right: BigUint) -> BigUint {
    while !right.is_zero() {
        let remainder = &left % &right;
        left = right;
        right = remainder;
    }
    left
}

fn mul_mod_u64(left: u64, right: u64, modulus: u64) -> u64 {
    ((left as u128 * right as u128) % modulus as u128) as u64
}

fn pow_mod_u64(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1u64;
    base %= modulus;
    while exponent != 0 {
        if exponent & 1 == 1 {
            result = mul_mod_u64(result, base, modulus);
        }
        base = mul_mod_u64(base, base, modulus);
        exponent >>= 1;
    }
    result
}

/// Deterministic Miller--Rabin on the complete `u64` range.
fn is_prime_u64(value: u64) -> bool {
    if value < 2 {
        return false;
    }
    for prime in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if value == prime {
            return true;
        }
        if value.is_multiple_of(prime) {
            return false;
        }
    }
    let mut odd = value - 1;
    let powers = odd.trailing_zeros();
    odd >>= powers;
    for witness in [2u64, 325, 9_375, 28_178, 450_775, 9_780_504, 1_795_265_022] {
        if witness.is_multiple_of(value) {
            continue;
        }
        let mut x = pow_mod_u64(witness % value, odd, value);
        if x == 1 || x == value - 1 {
            continue;
        }
        let mut accepted = false;
        for _ in 1..powers {
            x = mul_mod_u64(x, x, value);
            if x == value - 1 {
                accepted = true;
                break;
            }
        }
        if !accepted {
            return false;
        }
    }
    true
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct LucasWitness {
    pub prime_factor: String,
    pub witness: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct PrimalityCertificate {
    pub value: String,
    pub factorization_of_value_minus_one: Vec<(String, u32)>,
    pub witnesses: Vec<LucasWitness>,
    pub exact_product: bool,
    pub verified_prime: bool,
}

/// Lucas's complete `n-1` criterion.  Each factor is independently prime by
/// deterministic `u64` Miller--Rabin, and the full factorization is checked.
fn prove_prime_complete_n_minus_one(
    value: &BigUint,
    factors: &[(u64, u32)],
) -> Result<PrimalityCertificate, String> {
    let mut product = BigUint::one();
    for (prime, exponent) in factors {
        if !is_prime_u64(*prime) {
            return Err(format!("unproved factor {prime} of {value}-1"));
        }
        product *= BigUint::from(*prime).pow(*exponent);
    }
    if product != value - BigUint::one() {
        return Err(format!("factorization does not equal {value}-1"));
    }
    let value_minus_one = value - BigUint::one();
    let mut witnesses = Vec::new();
    for (prime, _) in factors {
        let exponent = &value_minus_one / *prime;
        let mut selected = None;
        for candidate in 2u64..100_000 {
            let base = BigUint::from(candidate);
            if base.modpow(&value_minus_one, value) != BigUint::one() {
                continue;
            }
            let residue = base.modpow(&exponent, value);
            let difference = if residue.is_zero() {
                value - BigUint::one()
            } else {
                residue - BigUint::one()
            };
            if gcd_big(difference, value.clone()) == BigUint::one() {
                selected = Some(candidate);
                break;
            }
        }
        let witness = selected.ok_or_else(|| format!("no Lucas witness for q={prime}"))?;
        witnesses.push(LucasWitness {
            prime_factor: prime.to_string(),
            witness: witness.to_string(),
        });
    }
    Ok(PrimalityCertificate {
        value: value.to_string(),
        factorization_of_value_minus_one: factors
            .iter()
            .map(|(prime, exponent)| (prime.to_string(), *exponent))
            .collect(),
        witnesses,
        exact_product: true,
        verified_prime: true,
    })
}

fn gf2_rem(mut value: BigUint, modulus: &BigUint) -> BigUint {
    let modulus_bits = modulus.bits();
    while !value.is_zero() && value.bits() >= modulus_bits {
        let shift = value.bits() - modulus_bits;
        value ^= modulus << shift;
    }
    value
}

fn gf2_square_mod(value: &BigUint, modulus: &BigUint) -> BigUint {
    let mut square = BigUint::zero();
    for bit in 0..value.bits() {
        if value.bit(bit) {
            square.set_bit(2 * bit, true);
        }
    }
    gf2_rem(square, modulus)
}

fn gf2_gcd(mut left: BigUint, mut right: BigUint) -> BigUint {
    while !right.is_zero() {
        let remainder = gf2_rem(left, &right);
        left = right;
        right = remainder;
    }
    left
}

fn certify_field_polynomial() -> bool {
    // f(z) = z^113 + z^9 + 1.  Since 113 is prime, Rabin's criterion
    // needs x^(2^113)=x mod f and gcd(x^2+x,f)=1.
    let f = (BigUint::one() << 113usize) | (BigUint::one() << 9usize) | BigUint::one();
    let x = BigUint::from(2u8);
    let x2 = gf2_square_mod(&x, &f);
    if gf2_gcd(x2 ^ &x, f.clone()) != BigUint::one() {
        return false;
    }
    let mut power = x.clone();
    for _ in 0..113 {
        power = gf2_square_mod(&power, &f);
    }
    power == x
}

fn field_trace(value: &F2mElement, curve: &BinaryCurve) -> F2mElement {
    let mut sum = F2mElement::zero(curve.m);
    let mut conjugate = value.clone();
    for _ in 0..curve.m {
        sum.add_assign(&conjugate);
        conjugate = conjugate.square(&curve.irreducible);
    }
    sum
}

fn field_sqrt(value: &F2mElement, curve: &BinaryCurve) -> F2mElement {
    value.square_k_times(curve.m - 1, &curve.irreducible)
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct SourceCertificate {
    pub source_url: String,
    pub source_sha256: String,
    pub identity: CurveIdentityCertificate,
    pub q: String,
    pub subgroup_order_n: String,
    pub curve_order: String,
    pub trace: String,
    pub frobenius_discriminant: String,
    pub field_polynomial_irreducible: bool,
    pub n_primality: PrimalityCertificate,
    pub generator_on_curve: bool,
    pub generator_exact_order_n: bool,
    pub hasse_unique_multiple: bool,
    pub embedding_degree: String,
    pub embedding_degree_exact: bool,
    pub discriminant_factorization_exact: bool,
    pub discriminant_small_factors_prime: bool,
    pub discriminant_fundamental: bool,
    pub discriminant_large_factor_primality: PrimalityCertificate,
    pub minimum_non_scalar_cm_norm: String,
    pub no_f2_descent_representative: bool,
    pub f2_base_change_traces: Vec<String>,
    pub order_of_2_mod_113: u32,
    pub ghs_magic_number_computed: bool,
    pub ghs_screen_status: String,
}

fn trace_recurrence(base_trace: i64, degree: usize) -> BigInt {
    if degree == 0 {
        return BigInt::from(2u8);
    }
    if degree == 1 {
        return BigInt::from(base_trace);
    }
    let s = BigInt::from(base_trace);
    let mut previous_previous = BigInt::from(2u8);
    let mut previous = s.clone();
    for _ in 2..=degree {
        let next = &s * &previous - (BigInt::from(2u8) * &previous_previous);
        previous_previous = previous;
        previous = next;
    }
    previous
}

fn certify_source(curve: &BinaryCurve) -> Result<SourceCertificate, String> {
    let identity = binary_curve_identity(curve, "sect113r1")?;
    if identity.icv1 != CURVE_ICV1
        || identity.curve_id != CURVE_ID
        || identity.curve_uid != CURVE_UID
    {
        return Err("sect113r1 source identity differs from the frozen registry identity".into());
    }
    let q = BigUint::one() << curve.m as usize;
    let n = curve.order.clone();
    let order = &n * 2u8;
    let trace = BigInt::from_biguint(Sign::Plus, &q + BigUint::one())
        - BigInt::from_biguint(Sign::Plus, order.clone());
    let expected_trace = BigInt::from(-122_610_772_499_221_213i64);
    if trace != expected_trace {
        return Err("sect113r1 trace differs from frozen value".into());
    }
    let n_factors = [
        (2, 1),
        (7, 1),
        (53, 1),
        (547, 1),
        (2_848_799, 1),
        (8_757_107, 1),
        (512_797_404_440_011, 1),
    ];
    let n_primality = prove_prime_complete_n_minus_one(&n, &n_factors)?;

    let generator_on_curve = curve.is_on_curve(&curve.generator);
    let generator_exact_order_n = !matches!(curve.generator, BinaryPoint::Infinity)
        && matches!(
            scalar_mul(curve, &curve.generator, &n),
            BinaryPoint::Infinity
        );
    // A point of prime order n makes #E a multiple of n.  The interval has
    // length < n, and 2n lies inside it, hence #E=2n.
    let trace_squared = expected_trace.magnitude().pow(2);
    let hasse_unique_multiple = trace_squared <= &q * 4u8 && (&q * 16u8) < n.pow(2);

    let k = (&n - BigUint::one()) / 2u8;
    let embedding_prime_factors = [7u64, 53, 547, 2_848_799, 8_757_107, 512_797_404_440_011];
    let q_mod_n = &q % &n;
    let embedding_degree_exact = q_mod_n.modpow(&k, &n) == BigUint::one()
        && embedding_prime_factors
            .iter()
            .all(|prime| q_mod_n.modpow(&(&k / *prime), &n) != BigUint::one());

    let discriminant =
        &expected_trace * &expected_trace - BigInt::from_biguint(Sign::Plus, &q * 4u8);
    let expected_discriminant = BigInt::parse_bytes(b"-26504973335422840129609279124569399", 10)
        .expect("valid discriminant");
    if discriminant != expected_discriminant {
        return Err("Frobenius discriminant mismatch".into());
    }
    let d_large = decimal("195708607903277035369399");
    let d_large_factors = [(2, 1), (3, 2), (16_370_128_667, 1), (664_179_290_233, 1)];
    let discriminant_large_factor_primality =
        prove_prime_complete_n_minus_one(&d_large, &d_large_factors)?;
    let d_product = BigUint::from(7u8) * 47u8 * 411_643_769u64 * d_large;
    let discriminant_factorization_exact = &d_product == expected_discriminant.magnitude();
    let discriminant_small_factors_prime = [7u64, 47, 411_643_769]
        .iter()
        .all(|factor| is_prime_u64(*factor));
    let discriminant_fundamental = discriminant_factorization_exact
        && discriminant_small_factors_prime
        && discriminant_large_factor_primality.verified_prime
        && expected_discriminant.magnitude() % 4u8 == BigUint::from(3u8);
    let minimum_non_scalar_cm_norm = (expected_discriminant.magnitude() + BigUint::one()) / 4u8;

    let f2_base_change_traces: Vec<String> = (-2..=2)
        .map(|base_trace| trace_recurrence(base_trace, 113).to_string())
        .collect();
    let no_f2_descent_representative = !f2_base_change_traces.contains(&expected_trace.to_string());

    // 2^28 = 1 mod 113 and no proper divisor of 28 has that property.
    let order_of_2_mod_113 = (1u32..=112)
        .find(|exponent| {
            BigUint::from(2u8).modpow(&BigUint::from(*exponent), &BigUint::from(113u8))
                == BigUint::one()
        })
        .ok_or("2 has no order mod 113")?;

    Ok(SourceCertificate {
        source_url: SOURCE_URL.into(),
        source_sha256: SOURCE_SHA256.into(),
        identity,
        q: q.to_string(),
        subgroup_order_n: n.to_string(),
        curve_order: order.to_string(),
        trace: expected_trace.to_string(),
        frobenius_discriminant: expected_discriminant.to_string(),
        field_polynomial_irreducible: certify_field_polynomial(),
        n_primality,
        generator_on_curve,
        generator_exact_order_n,
        hasse_unique_multiple,
        embedding_degree: k.to_string(),
        embedding_degree_exact,
        discriminant_factorization_exact,
        discriminant_small_factors_prime,
        discriminant_fundamental,
        discriminant_large_factor_primality,
        minimum_non_scalar_cm_norm: minimum_non_scalar_cm_norm.to_string(),
        no_f2_descent_representative,
        f2_base_change_traces,
        order_of_2_mod_113,
        ghs_magic_number_computed: false,
        ghs_screen_status: "NOT_CERTIFIED_MAGIC_NUMBER_UNCOMPUTED".into(),
    })
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
enum PointKey {
    Infinity,
    Affine(Vec<u8>, Vec<u8>),
}

fn point_key(point: &BinaryPoint) -> PointKey {
    match point {
        BinaryPoint::Infinity => PointKey::Infinity,
        BinaryPoint::Affine { x, y } => {
            PointKey::Affine(x.to_biguint().to_bytes_be(), y.to_biguint().to_bytes_be())
        }
    }
}

fn ceil_sqrt_u64(value: u64) -> u64 {
    let mut root = (value as f64).sqrt() as u64;
    while (root as u128) * (root as u128) < value as u128 {
        root += 1;
    }
    while root > 0 && ((root - 1) as u128) * ((root - 1) as u128) >= value as u128 {
        root -= 1;
    }
    root
}

fn bsgs_prime_order(
    curve: &BinaryCurve,
    generator: &BinaryPoint,
    target: &BinaryPoint,
    order: u64,
) -> Result<(u64, u64), String> {
    if !is_prime_u64(order) {
        return Err("BSGS component order is not prime".into());
    }
    let width = ceil_sqrt_u64(order);
    let mut babies = HashMap::with_capacity(width as usize);
    let mut current = BinaryPoint::Infinity;
    for index in 0..width {
        babies.entry(point_key(&current)).or_insert(index);
        current = point_add(curve, &current, generator);
    }
    let stride = scalar_mul(curve, generator, &BigUint::from(width));
    let negative_stride = point_neg(&stride);
    let mut giant = target.clone();
    let mut comparisons = 0u64;
    for high in 0..=width {
        comparisons += 1;
        if let Some(low) = babies.get(&point_key(&giant)) {
            let candidate = high
                .checked_mul(width)
                .and_then(|value| value.checked_add(*low))
                .ok_or("BSGS candidate overflow")?;
            if candidate < order
                && scalar_mul(curve, generator, &BigUint::from(candidate)) == *target
            {
                return Ok((candidate, comparisons));
            }
        }
        giant = point_add(curve, &giant, &negative_stride);
    }
    Err(format!("BSGS found no logarithm modulo {order}"))
}

fn inverse_mod_u64(value: u64, modulus: u64) -> Result<u64, String> {
    let (mut old_r, mut r) = (value as i128, modulus as i128);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
    }
    if old_r != 1 {
        return Err("CRT modulus is not invertible".into());
    }
    Ok(old_s.rem_euclid(modulus as i128) as u64)
}

fn crt(residues: &[(u64, u64)]) -> Result<(u128, u128), String> {
    let (mut solution, mut modulus) = (0u128, 1u128);
    for (prime, residue) in residues {
        let p = *prime as u128;
        let solution_mod_p = solution % p;
        let delta = ((*residue as u128 + p - solution_mod_p) % p) as u64;
        let inverse = inverse_mod_u64((modulus % p) as u64, *prime)?;
        let lift = mul_mod_u64(delta, inverse, *prime) as u128;
        solution += modulus * lift;
        modulus *= p;
        solution %= modulus;
    }
    Ok((solution, modulus))
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct PohligHellmanComponent {
    pub prime: u64,
    pub residue: u64,
    pub expected_residue: u64,
    pub bsgs_giant_steps: u64,
    pub exact_replay: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct SingularCompanionCertificate {
    pub equation: String,
    pub discriminant_zero: bool,
    pub trace_a_is_one: bool,
    pub smooth_locus_order: String,
    pub factorization_exact: bool,
    pub factor_primes_exact: bool,
    pub attack_point: String,
    pub attack_point_off_nominal_curve: bool,
    pub attack_point_on_singular_companion: bool,
    pub attack_point_exact_order: bool,
    pub nominal_and_companion_multipliers_identical: bool,
    pub independent_torus_matches_demo: bool,
    pub independent_torus_matches_full_width: bool,
    pub exhaustive_small_field_torus_control: bool,
    pub actual_oracle_queries_demo: u64,
    pub demo_seed_digest: String,
    pub demo_scalar: String,
    pub demo_result: String,
    pub demo_public_key: String,
    pub demo_components: Vec<PohligHellmanComponent>,
    pub demo_crt_modulus: String,
    pub demo_recovered_scalar: String,
    pub demo_exact_recovery: bool,
    pub demo_final_replay_matches_public_key: bool,
    pub full_width_seed_digest: String,
    pub full_width_scalar: String,
    pub full_width_result: String,
    pub full_width_public_key: String,
    pub full_width_frozen_vectors_match: bool,
    pub actual_oracle_queries_full_width: u64,
    pub full_width_small_components: Vec<PohligHellmanComponent>,
    pub full_width_recovered_residue: String,
    pub full_width_remaining_prime: u64,
    pub full_width_large_component_executed: bool,
    pub full_width_rho_expected_log2: String,
    pub validator_rejects_attack_point: bool,
    pub order_two_point: String,
    pub order_two_exact: bool,
    pub parity_observable_verified: bool,
    pub validator_rejects_order_two_point: bool,
    pub verdict: String,
    pub applicability: String,
}

fn singular_curve(nominal: &BinaryCurve, order: &BigUint, generator: BinaryPoint) -> BinaryCurve {
    BinaryCurve {
        m: nominal.m,
        irreducible: nominal.irreducible.clone(),
        a: nominal.a.clone(),
        b: F2mElement::zero(nominal.m),
        generator,
        order: order.clone(),
        cofactor: BigUint::one(),
    }
}

/// Independent parameterization of the nonsplit nodal cubic's smooth locus.
/// `None` is the identity and `Some(t)` maps to
/// `(t^2+t+a, t(t^2+t+a))`.  This code never invokes an elliptic point-law
/// routine.
fn torus_add(
    left: &Option<F2mElement>,
    right: &Option<F2mElement>,
    curve: &BinaryCurve,
) -> Option<F2mElement> {
    match (left, right) {
        (None, other) | (other, None) => other.clone(),
        (Some(left), Some(right)) => {
            let denominator = left.add(right).add(&F2mElement::one(curve.m));
            if denominator.is_zero() {
                return None;
            }
            let numerator = left.mul(right, &curve.irreducible).add(&curve.a);
            Some(
                numerator.mul(
                    &denominator
                        .flt_inverse(&curve.irreducible)
                        .expect("nonzero torus denominator"),
                    &curve.irreducible,
                ),
            )
        }
    }
}

fn torus_scalar_mul(
    curve: &BinaryCurve,
    parameter: &F2mElement,
    scalar: &BigUint,
) -> Option<F2mElement> {
    let base = Some(parameter.clone());
    let mut result = None;
    for bit in (0..scalar.bits()).rev() {
        result = torus_add(&result, &result, curve);
        if scalar.bit(bit) {
            result = torus_add(&result, &base, curve);
        }
    }
    result
}

fn torus_to_point(parameter: &Option<F2mElement>, curve: &BinaryCurve) -> BinaryPoint {
    let Some(parameter) = parameter else {
        return BinaryPoint::Infinity;
    };
    let x = parameter
        .square(&curve.irreducible)
        .add(parameter)
        .add(&curve.a);
    let y = parameter.mul(&x, &curve.irreducible);
    BinaryPoint::Affine { x, y }
}

fn exhaustive_small_field_torus_control() -> bool {
    let m = 8;
    let irreducible = crate::binary_ecc::IrreduciblePoly::deg_8();
    let template = BinaryCurve {
        m,
        irreducible,
        a: F2mElement::zero(m),
        b: F2mElement::zero(m),
        generator: BinaryPoint::Infinity,
        order: BigUint::from(257u16),
        cofactor: BigUint::one(),
    };
    let Some(a) = (0u16..=255)
        .map(|value| F2mElement::from_biguint(&BigUint::from(value), m))
        .find(|candidate| field_trace(candidate, &template) == F2mElement::one(m))
    else {
        return false;
    };
    let curve = BinaryCurve { a, ..template };
    let mut parameters = Vec::with_capacity(257);
    parameters.push(None);
    parameters
        .extend((0u16..=255).map(|value| Some(F2mElement::from_biguint(&BigUint::from(value), m))));
    if parameters.iter().any(|parameter| {
        let point = torus_to_point(parameter, &curve);
        !curve.is_on_curve(&point)
    }) {
        return false;
    }
    for left in &parameters {
        for right in &parameters {
            let expected = torus_to_point(&torus_add(left, right, &curve), &curve);
            let actual = point_add(
                &curve,
                &torus_to_point(left, &curve),
                &torus_to_point(right, &curve),
            );
            if actual != expected {
                return false;
            }
        }
    }
    true
}

fn recover_components(
    curve: &BinaryCurve,
    base: &BinaryPoint,
    target: &BinaryPoint,
    group_order: &BigUint,
    secret: &BigUint,
) -> Result<Vec<PohligHellmanComponent>, String> {
    let mut output = Vec::new();
    for prime in EXECUTED_FACTORS {
        let quotient = group_order / prime;
        let component_base = scalar_mul(curve, base, &quotient);
        let component_target = scalar_mul(curve, target, &quotient);
        let (residue, steps) = bsgs_prime_order(curve, &component_base, &component_target, prime)?;
        let expected = (secret % prime)
            .to_u64()
            .ok_or("residue does not fit u64")?;
        output.push(PohligHellmanComponent {
            prime,
            residue,
            expected_residue: expected,
            bsgs_giant_steps: steps,
            exact_replay: residue == expected
                && scalar_mul(curve, &component_base, &BigUint::from(residue)) == component_target,
        });
    }
    Ok(output)
}

fn certify_singular_companion(
    nominal: &BinaryCurve,
) -> Result<SingularCompanionCertificate, String> {
    let q = BigUint::one() << nominal.m as usize;
    let smooth_order = &q + BigUint::one();
    let factor_product: BigUint = SINGULAR_FACTORS
        .iter()
        .map(|prime| BigUint::from(*prime))
        .product();
    let factorization_exact = factor_product == smooth_order;
    let factor_primes_exact = SINGULAR_FACTORS.iter().all(|prime| is_prime_u64(*prime));
    let attack_point = BinaryPoint::Affine {
        x: nominal.a.clone(),
        y: F2mElement::zero(nominal.m),
    };
    let companion = singular_curve(nominal, &smooth_order, attack_point.clone());
    let attack_point_exact_order = matches!(
        scalar_mul(&companion, &attack_point, &smooth_order),
        BinaryPoint::Infinity
    ) && SINGULAR_FACTORS.iter().all(|prime| {
        !matches!(
            scalar_mul(&companion, &attack_point, &(&smooth_order / *prime)),
            BinaryPoint::Infinity
        )
    });

    let demo_digest = sha256(DEMO_DOMAIN);
    let demo_modulus_expected: BigUint = EXECUTED_FACTORS
        .iter()
        .map(|prime| BigUint::from(*prime))
        .product();
    let demo_secret = BigUint::one()
        + BigUint::from_bytes_be(&demo_digest) % ((&demo_modulus_expected / 2u8) - BigUint::one());
    if hex::encode(demo_digest) != DEMO_DIGEST || demo_secret != BigUint::from(DEMO_SCALAR) {
        return Err("demonstration scalar derivation changed".into());
    }
    let demo_result = scalar_mul(nominal, &attack_point, &demo_secret); // actual unchecked call
    let demo_public_key = scalar_mul(nominal, &nominal.generator, &demo_secret);
    let demo_companion_result = scalar_mul(&companion, &attack_point, &demo_secret);
    let independent_demo = torus_to_point(
        &torus_scalar_mul(&companion, &F2mElement::zero(nominal.m), &demo_secret),
        &companion,
    );
    let demo_components = recover_components(
        &companion,
        &attack_point,
        &demo_result,
        &smooth_order,
        &demo_secret,
    )?;
    let demo_pairs: Vec<(u64, u64)> = demo_components
        .iter()
        .map(|component| (component.prime, component.residue))
        .collect();
    let (demo_recovered, demo_modulus) = crt(&demo_pairs)?;
    let demo_replay = scalar_mul(nominal, &nominal.generator, &BigUint::from(demo_recovered));

    let mut seed_material = Vec::with_capacity(FULL_DOMAIN.len() + FULL_SEED.len());
    seed_material.extend_from_slice(FULL_DOMAIN);
    seed_material.extend_from_slice(FULL_SEED);
    let digest = sha256(&seed_material);
    let derived_full =
        BigUint::one() + BigUint::from_bytes_be(&digest) % (&nominal.order - BigUint::one());
    if hex::encode(digest) != FULL_DIGEST || derived_full != decimal(FULL_SCALAR) {
        return Err("full-width planted scalar derivation changed".into());
    }
    let full_result = scalar_mul(nominal, &attack_point, &derived_full); // one actual call
    let full_public_key = scalar_mul(nominal, &nominal.generator, &derived_full);
    let independent_full = torus_to_point(
        &torus_scalar_mul(&companion, &F2mElement::zero(nominal.m), &derived_full),
        &companion,
    );
    let full_components = recover_components(
        &companion,
        &attack_point,
        &full_result,
        &smooth_order,
        &derived_full,
    )?;
    let full_pairs: Vec<(u64, u64)> = full_components
        .iter()
        .map(|component| (component.prime, component.residue))
        .collect();
    let (full_residue, _) = crt(&full_pairs)?;

    let order_two = BinaryPoint::Affine {
        x: F2mElement::zero(nominal.m),
        y: field_sqrt(&nominal.b, nominal),
    };
    let order_two_exact = nominal.is_on_curve(&order_two)
        && !matches!(order_two, BinaryPoint::Infinity)
        && matches!(
            scalar_mul(nominal, &order_two, &BigUint::from(2u8)),
            BinaryPoint::Infinity
        )
        && scalar_mul(nominal, &order_two, &nominal.order) == order_two;
    let parity_result = scalar_mul(nominal, &order_two, &derived_full);
    let parity_observable_verified = if derived_full.bit(0) {
        parity_result == order_two
    } else {
        matches!(parity_result, BinaryPoint::Infinity)
    };
    let full_width_frozen_vectors_match = match (&full_public_key, &full_result) {
        (BinaryPoint::Affine { x: qx, y: qy }, BinaryPoint::Affine { x: ux, y: uy }) => {
            field_hex(qx) == "00af667d2d64df714fcb7229d4509f"
                && field_hex(qy) == "0071c688f4f136aaf006193d1dc24f"
                && field_hex(ux) == "0078aed9a682be670a024a24d76b1d"
                && field_hex(uy) == "00d567fa101a19af65ce7f58c0d4a4"
        }
        _ => false,
    };

    Ok(SingularCompanionCertificate {
        equation: "y^2+x*y=x^3+a*x^2 (b=0; singular nodal cubic)".into(),
        discriminant_zero: companion.b.is_zero(),
        trace_a_is_one: field_trace(&nominal.a, nominal) == F2mElement::one(nominal.m),
        smooth_locus_order: smooth_order.to_string(),
        factorization_exact,
        factor_primes_exact,
        attack_point: point_hex(&attack_point),
        attack_point_off_nominal_curve: !nominal.is_on_curve(&attack_point),
        attack_point_on_singular_companion: companion.is_on_curve(&attack_point),
        attack_point_exact_order,
        nominal_and_companion_multipliers_identical: demo_result == demo_companion_result,
        independent_torus_matches_demo: demo_result == independent_demo,
        independent_torus_matches_full_width: full_result == independent_full,
        exhaustive_small_field_torus_control: exhaustive_small_field_torus_control(),
        actual_oracle_queries_demo: 1,
        demo_seed_digest: DEMO_DIGEST.into(),
        demo_scalar: demo_secret.to_string(),
        demo_result: point_hex(&demo_result),
        demo_public_key: point_hex(&demo_public_key),
        demo_components,
        demo_crt_modulus: demo_modulus.to_string(),
        demo_recovered_scalar: demo_recovered.to_string(),
        demo_exact_recovery: demo_recovered == DEMO_SCALAR as u128
            && (DEMO_SCALAR as u128) < demo_modulus,
        demo_final_replay_matches_public_key: demo_replay == demo_public_key,
        full_width_seed_digest: FULL_DIGEST.into(),
        full_width_scalar: derived_full.to_string(),
        full_width_result: point_hex(&full_result),
        full_width_public_key: point_hex(&full_public_key),
        full_width_frozen_vectors_match,
        actual_oracle_queries_full_width: 1,
        full_width_small_components: full_components,
        full_width_recovered_residue: full_residue.to_string(),
        full_width_remaining_prime: SINGULAR_FACTORS[4],
        full_width_large_component_executed: false,
        full_width_rho_expected_log2: "29.710003".into(),
        validator_rejects_attack_point: !nominal.is_valid_public_point(&attack_point),
        order_two_point: point_hex(&order_two),
        order_two_exact,
        parity_observable_verified,
        validator_rejects_order_two_point: !nominal.is_valid_public_point(&order_two),
        verdict: "CONDITIONAL_IMPLEMENTATION_WEAK".into(),
        applicability: "conditional on attacker-selected BinaryPoint reaching unchecked scalar_mul and a distinguishable full-point result; no production protocol exposure is asserted".into(),
    })
}

fn derivative(poly: &F2mPoly) -> F2mPoly {
    let Some(degree) = poly.degree() else {
        return F2mPoly::zero(poly.m);
    };
    let mut coefficients = vec![F2mElement::zero(poly.m); degree];
    for exponent in 1..=degree {
        if exponent % 2 == 1 {
            coefficients[exponent - 1] = poly.coeff(exponent);
        }
    }
    F2mPoly::from_coeffs(coefficients, poly.m)
}

fn frobenius_power(
    mut value: F2mPoly,
    times: u32,
    modulus: &F2mPoly,
    curve: &BinaryCurve,
) -> F2mPoly {
    for _ in 0..times {
        value = value.square_mod(modulus, &curve.irreducible);
    }
    value
}

fn degree_five_kernels(curve: &BinaryCurve) -> Result<(F2mPoly, F2mPoly, Vec<usize>), String> {
    let psi5 = division_polynomial(&curve.b, 5, curve.m, &curve.irreducible);
    if psi5.degree() != Some(12) {
        return Err("psi_5 has unexpected degree".into());
    }
    let x = F2mPoly::x(curve.m);
    let x_q = frobenius_power(x.clone(), curve.m, &psi5, curve);
    let h1 = psi5.gcd(&x_q.add(&x), &curve.irreducible);
    let x_q2 = frobenius_power(x_q, curve.m, &psi5, curve);
    let g2 = psi5.gcd(&x_q2.add(&x), &curve.irreducible);
    let (h2, remainder) = g2.divrem(&h1, &curve.irreducible);
    if !remainder.is_zero() {
        return Err("q^2 fixed factor is not divisible by q fixed factor".into());
    }
    let x_q4 = frobenius_power(x_q2, curve.m * 2, &psi5, curve);
    let g4 = psi5.gcd(&x_q4.add(&x), &curve.irreducible);
    let degrees = vec![
        h1.degree().unwrap_or(0),
        g2.degree().unwrap_or(0),
        g4.degree().unwrap_or(0),
    ];
    if h1.degree() != Some(2) || h2.degree() != Some(2) || degrees != [2, 4, 12] {
        return Err(format!("unexpected Frobenius factor degrees {degrees:?}"));
    }
    Ok((h1, h2, degrees))
}

fn kernel_digest(kernel: &F2mPoly) -> String {
    let mut input = b"crypto/sect113r1/degree5-kernel/v1\0".to_vec();
    input.extend_from_slice(&(kernel.degree().unwrap_or(0) as u64).to_be_bytes());
    for index in 0..=kernel.degree().unwrap_or(0) {
        input.extend_from_slice(&fixed_field_bytes(&kernel.coeff(index)));
    }
    hex::encode(sha256(&input))
}

fn kernel_closed_under_doubling(kernel: &F2mPoly, b: &F2mElement, curve: &BinaryCurve) -> bool {
    if kernel.degree() != Some(2) {
        return false;
    }
    let x = F2mPoly::x(curve.m);
    let x2 = x.mul(&x, &curve.irreducible);
    let x4 = x2.mul(&x2, &curve.irreducible);
    let z = x4.add(&F2mPoly::constant(b.clone()));
    let expression = z
        .mul(&z, &curve.irreducible)
        .add(
            &z.mul(&x2, &curve.irreducible)
                .scalar_mul(&kernel.coeff(1), &curve.irreducible),
        )
        .add(&x4.scalar_mul(&kernel.coeff(0), &curve.irreducible));
    expression.rem(kernel, &curve.irreducible).is_zero()
}

fn kernel_certificate_ok(kernel: &F2mPoly, psi5: &F2mPoly, curve: &BinaryCurve) -> bool {
    let x = F2mPoly::x(curve.m);
    let (_, remainder) = psi5.divrem(kernel, &curve.irreducible);
    kernel.m == curve.m
        && kernel.degree() == Some(2)
        && kernel.lead() == F2mElement::one(curve.m)
        && remainder.is_zero()
        && kernel.gcd(&derivative(kernel), &curve.irreducible) == F2mPoly::one(curve.m)
        && kernel.gcd(&x, &curve.irreducible) == F2mPoly::one(curve.m)
        && kernel_closed_under_doubling(kernel, &curve.b, curve)
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct IsogenyCertificate {
    pub degree: u64,
    pub map_implementation: String,
    pub kernel_coefficients_low_first: Vec<String>,
    pub kernel_sha256: String,
    pub kernel_certificate_exact: bool,
    pub codomain_b: String,
    pub codomain_nonsingular: bool,
    pub codomain_identity: CurveIdentityCertificate,
    pub generator_image: String,
    pub generator_image_on_codomain: bool,
    pub generator_image_exact_order_n: bool,
    pub planted_target_forward_homomorphism: bool,
    pub dual_kernel_sha256: String,
    pub dual_kernel_certificate_exact: bool,
    pub dual_returns_exact_source_model: bool,
    pub dual_composition_on_generator_is_times_five: bool,
    pub dual_composition_on_planted_target_is_times_five: bool,
    pub planted_target_pullback_to_source: bool,
    pub dlp_pullback_formula: String,
    pub verdict: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct IsogenyClassCertificate {
    pub frobenius_factor_degrees_q_q2_q4: Vec<usize>,
    pub rational_degree_five_kernel_count: usize,
    pub representatives: Vec<IsogenyCertificate>,
    pub conductor_strata_absent: bool,
    pub useful_low_degree_cm_endomorphism_excluded: bool,
    pub class_invariant_subgroup_order: String,
}

fn make_codomain(domain: &BinaryCurve, b: F2mElement, generator: BinaryPoint) -> BinaryCurve {
    BinaryCurve {
        m: domain.m,
        irreducible: domain.irreducible.clone(),
        a: domain.a.clone(),
        b,
        generator,
        order: domain.order.clone(),
        cofactor: domain.cofactor.clone(),
    }
}

fn certify_isogenies(
    source: &BinaryCurve,
    full_scalar: &BigUint,
) -> Result<IsogenyClassCertificate, String> {
    let (kernel1, kernel2, degrees) = degree_five_kernels(source)?;
    let psi5 = division_polynomial(&source.b, 5, source.m, &source.irreducible);
    let source_target = scalar_mul(source, &source.generator, full_scalar);
    let mut representatives = Vec::new();
    for kernel in [kernel1, kernel2] {
        let kernel_certificate_exact = kernel_certificate_ok(&kernel, &psi5, source);
        if !kernel_certificate_exact {
            return Err("degree-5 kernel certificate failed".into());
        }
        let codomain_b = velu_codomain(&source.b, &kernel, source.m, &source.irreducible);
        let generator_image =
            velu_point_map_binary(&source.generator, &kernel, source.m, &source.irreducible)
                .ok_or("failed to map source generator")?;
        let codomain = make_codomain(source, codomain_b.clone(), generator_image.clone());
        let codomain_identity = binary_curve_identity(&codomain, "rb")?;
        let target_image =
            velu_point_map_binary(&source_target, &kernel, source.m, &source.irreducible)
                .ok_or("failed to map planted target")?;
        let planted_target_forward_homomorphism =
            scalar_mul(&codomain, &generator_image, full_scalar) == target_image;

        let (dual1, dual2, _) = degree_five_kernels(&codomain)?;
        let codomain_psi5 = division_polynomial(&codomain.b, 5, codomain.m, &codomain.irreducible);
        let mut dual_result = None;
        for dual in [dual1, dual2] {
            let dual_kernel_certificate_exact =
                kernel_certificate_ok(&dual, &codomain_psi5, &codomain);
            if !dual_kernel_certificate_exact {
                continue;
            }
            let back_b = velu_codomain(&codomain.b, &dual, codomain.m, &codomain.irreducible);
            if back_b == source.b {
                let composed = velu_point_map_binary(
                    &generator_image,
                    &dual,
                    codomain.m,
                    &codomain.irreducible,
                )
                .ok_or("failed to map through dual")?;
                dual_result = Some((dual, back_b, composed, dual_kernel_certificate_exact));
                break;
            }
        }
        let (dual, back_b, composed, dual_kernel_certificate_exact) =
            dual_result.ok_or("no degree-5 dual returned to source model")?;
        let times_five = scalar_mul(source, &source.generator, &BigUint::from(5u8));
        let composed_target =
            velu_point_map_binary(&target_image, &dual, source.m, &source.irreducible)
                .ok_or("failed to map planted target through dual")?;
        let target_times_five = scalar_mul(source, &source_target, &BigUint::from(5u8));
        let inverse_five =
            BigUint::from(5u8).modpow(&(&source.order - BigUint::from(2u8)), &source.order);
        let planted_target_pullback_to_source =
            scalar_mul(source, &composed_target, &inverse_five) == source_target;
        representatives.push(IsogenyCertificate {
            degree: 5,
            map_implementation: "binary-velu/v1".into(),
            kernel_coefficients_low_first: (0..=2)
                .map(|index| field_hex(&kernel.coeff(index)))
                .collect(),
            kernel_sha256: kernel_digest(&kernel),
            kernel_certificate_exact,
            codomain_b: field_hex(&codomain_b),
            codomain_nonsingular: !codomain_b.is_zero(),
            codomain_identity,
            generator_image: point_hex(&generator_image),
            generator_image_on_codomain: codomain.is_on_curve(&generator_image),
            generator_image_exact_order_n: !matches!(generator_image, BinaryPoint::Infinity)
                && matches!(
                    scalar_mul(&codomain, &generator_image, &source.order),
                    BinaryPoint::Infinity
                ),
            planted_target_forward_homomorphism,
            dual_kernel_sha256: kernel_digest(&dual),
            dual_kernel_certificate_exact,
            dual_returns_exact_source_model: back_b == source.b,
            dual_composition_on_generator_is_times_five: composed == times_five,
            dual_composition_on_planted_target_is_times_five: composed_target == target_times_five,
            planted_target_pullback_to_source,
            dlp_pullback_formula: "P=[5^{-1} mod n]*dual(phi(P)) for P in <G>".into(),
            verdict: "TRANSFER_ONLY_SPEEDUP_NOT_ESTABLISHED".into(),
        });
    }
    Ok(IsogenyClassCertificate {
        frobenius_factor_degrees_q_q2_q4: degrees,
        rational_degree_five_kernel_count: representatives.len(),
        representatives,
        conductor_strata_absent: true,
        useful_low_degree_cm_endomorphism_excluded: true,
        class_invariant_subgroup_order: source.order.to_string(),
    })
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ClassWeaknessCertificate {
    pub frozen_classical_threshold_bits: u32,
    pub subgroup_bits: u64,
    pub rho_expected_log2_group_ops: String,
    pub rho_95pct_log2_group_ops: String,
    pub optimized_negation_rho_log2_group_ops: String,
    pub below_threshold: bool,
    pub published_full_dlp_recovery: String,
    pub published_recovery_url: String,
    pub verdict: String,
    pub novelty_boundary: String,
}

fn certify_class_weakness(n: &BigUint) -> ClassWeaknessCertificate {
    let n_float = n.to_f64().expect("113-bit integer fits f64 exponent range");
    let expected = (std::f64::consts::PI * n_float / 2.0).sqrt().log2();
    let p95 = (-2.0 * n_float * 0.05f64.ln()).sqrt().log2();
    let optimized = (0.886 * n_float.sqrt()).log2();
    ClassWeaknessCertificate {
        frozen_classical_threshold_bits: 80,
        subgroup_bits: n.bits(),
        rho_expected_log2_group_ops: format!("{expected:.6}"),
        rho_95pct_log2_group_ops: format!("{p95:.6}"),
        optimized_negation_rho_log2_group_ops: format!("{optimized:.6}"),
        below_threshold: optimized < 80.0,
        published_full_dlp_recovery: "Wenger and Wolfger report a full sect113r1 discrete-log computation in about 2.5 months on ten Kintex-7 FPGAs".into(),
        published_recovery_url: "https://eprint.iacr.org/2015/143".into(),
        verdict: "CLASS_WEAK".into(),
        novelty_boundary: "legacy subgroup size and generic attack; not a newly discovered algebraic shortcut".into(),
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TwistCertificate {
    pub order: String,
    pub derived_from_source_order: bool,
    pub factorization_exact: bool,
    pub largest_prime: String,
    pub largest_prime_primality: PrimalityCertificate,
    pub generic_rho_log2: String,
    pub classification: String,
}

fn certify_twist(source: &BinaryCurve) -> Result<TwistCertificate, String> {
    let q = BigUint::one() << source.m as usize;
    let source_order = &source.order * &source.cofactor;
    // Quadratic twisting negates the trace, hence
    // #E + #E_twist = 2(q+1). Derive the control order from the certified
    // source tuple instead of accepting the expected decimal as an input.
    let order = (&q + BigUint::one()) * 2u8 - &source_order;
    let derived_from_source_order = &source_order + &order == (&q + BigUint::one()) * 2u8;
    let largest = decimal("4690844705882931102817");
    let factors = [
        (2u64, 2u32),
        (5, 1),
        (11, 1),
        (17, 1),
        (449, 1),
        (883, 1),
        (1493, 1),
    ];
    let mut product = BigUint::one();
    for (prime, exponent) in factors {
        if !is_prime_u64(prime) {
            return Err(format!("twist factor {prime} is not prime"));
        }
        product *= BigUint::from(prime).pow(exponent);
    }
    product *= &largest;
    let largest_prime_primality = prove_prime_complete_n_minus_one(
        &largest,
        &[
            (2, 5),
            (3, 1),
            (173, 1),
            (47_653, 1),
            (5_927_116_621_409, 1),
        ],
    )?;
    let rho = (std::f64::consts::PI * largest.to_f64().expect("twist factor fits f64") / 2.0)
        .sqrt()
        .log2();
    Ok(TwistCertificate {
        order: order.to_string(),
        derived_from_source_order,
        factorization_exact: product == order,
        largest_prime: largest.to_string(),
        largest_prime_primality,
        generic_rho_log2: format!("{rho:.6}"),
        classification: "OFF_CLASS_CONTROL_NOT_AN_ISOGENOUS_REPRESENTATIVE".into(),
    })
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct AuditReport {
    pub schema: String,
    pub admitted_scientific_run: bool,
    pub evidence_status: String,
    pub source: SourceCertificate,
    pub class_weakness: ClassWeaknessCertificate,
    pub singular_companion: SingularCompanionCertificate,
    pub valid_isogenies: IsogenyClassCertificate,
    pub quadratic_twist: TwistCertificate,
    pub statistical_policy: String,
    pub implemented_diagnostic_gates_passed: bool,
    pub certificate_sha256: String,
}

fn report_digest(report: &AuditReport) -> Result<String, String> {
    let mut unsigned = report.clone();
    unsigned.certificate_sha256.clear();
    let encoded = serde_json::to_vec(&unsigned)
        .map_err(|error| format!("could not serialize audit report: {error}"))?;
    Ok(hex::encode(sha256(&encoded)))
}

fn mandatory_gates(report: &AuditReport) -> bool {
    report.source.field_polynomial_irreducible
        && report.source.n_primality.verified_prime
        && report.source.generator_on_curve
        && report.source.generator_exact_order_n
        && report.source.hasse_unique_multiple
        && report.source.embedding_degree_exact
        && report.source.discriminant_factorization_exact
        && report.source.discriminant_small_factors_prime
        && report.source.discriminant_fundamental
        && report.source.no_f2_descent_representative
        && report.source.order_of_2_mod_113 == 28
        && !report.source.ghs_magic_number_computed
        && report.source.ghs_screen_status == "NOT_CERTIFIED_MAGIC_NUMBER_UNCOMPUTED"
        && report.class_weakness.below_threshold
        && report.singular_companion.factorization_exact
        && report.singular_companion.factor_primes_exact
        && report.singular_companion.attack_point_off_nominal_curve
        && report.singular_companion.attack_point_on_singular_companion
        && report.singular_companion.attack_point_exact_order
        && report
            .singular_companion
            .nominal_and_companion_multipliers_identical
        && report.singular_companion.independent_torus_matches_demo
        && report
            .singular_companion
            .independent_torus_matches_full_width
        && report
            .singular_companion
            .exhaustive_small_field_torus_control
        && report.singular_companion.demo_exact_recovery
        && report
            .singular_companion
            .demo_final_replay_matches_public_key
        && report.singular_companion.full_width_frozen_vectors_match
        && report
            .singular_companion
            .demo_components
            .iter()
            .all(|component| component.exact_replay)
        && report
            .singular_companion
            .full_width_small_components
            .iter()
            .all(|component| component.exact_replay)
        && !report
            .singular_companion
            .full_width_large_component_executed
        && report.singular_companion.validator_rejects_attack_point
        && report.singular_companion.order_two_exact
        && report.singular_companion.parity_observable_verified
        && report.singular_companion.validator_rejects_order_two_point
        && report.valid_isogenies.rational_degree_five_kernel_count == 2
        && report
            .valid_isogenies
            .representatives
            .iter()
            .all(|certificate| {
                certificate.map_implementation == "binary-velu/v1"
                    && certificate.kernel_certificate_exact
                    && certificate.codomain_nonsingular
                    && certificate.codomain_identity.field_sha256
                        == "da55718b2ae51e38fc5d836fcf62bbca5907b877c61b29d5e3c30a5e94b60ee3"
                    && certificate.generator_image_on_codomain
                    && certificate.generator_image_exact_order_n
                    && certificate.planted_target_forward_homomorphism
                    && certificate.dual_kernel_certificate_exact
                    && certificate.dual_returns_exact_source_model
                    && certificate.dual_composition_on_generator_is_times_five
                    && certificate.dual_composition_on_planted_target_is_times_five
                    && certificate.planted_target_pullback_to_source
            })
        && report.quadratic_twist.derived_from_source_order
        && report.quadratic_twist.factorization_exact
        && report
            .quadratic_twist
            .largest_prime_primality
            .verified_prime
}

/// Run every exact certificate and bounded attack component.
///
/// The report stays marked `admitted_scientific_run=false` until the caller's
/// external protocol has frozen, committed and published the exact inputs.
pub fn run_audit() -> Result<AuditReport, String> {
    let curve = BinaryCurve::sect113r1();
    let source = certify_source(&curve)?;
    let class_weakness = certify_class_weakness(&curve.order);
    let singular_companion = certify_singular_companion(&curve)?;
    let full_scalar = decimal(FULL_SCALAR);
    let valid_isogenies = certify_isogenies(&curve, &full_scalar)?;
    let quadratic_twist = certify_twist(&curve)?;
    let mut report = AuditReport {
        schema: REPORT_SCHEMA.into(),
        admitted_scientific_run: false,
        evidence_status: "PRE_ADMISSION_DIAGNOSTIC_EXACT_CERTIFICATES".into(),
        source,
        class_weakness,
        singular_companion,
        valid_isogenies,
        quadratic_twist,
        statistical_policy: "exact algebra and replayed discrete logs need no p-value; randomized rho runtime is a model and its unexecuted large component remains explicitly unexecuted; no prevalence claim is made".into(),
        implemented_diagnostic_gates_passed: false,
        certificate_sha256: String::new(),
    };
    report.implemented_diagnostic_gates_passed = mandatory_gates(&report);
    report.certificate_sha256 = report_digest(&report)?;
    if !report.implemented_diagnostic_gates_passed {
        return Err("one or more mandatory sect113r1 audit gates failed".into());
    }
    Ok(report)
}

/// Deterministically recompute the complete report in the same implementation
/// and compare it byte-for-byte at the typed-data level.
pub fn verify_report(report: &AuditReport) -> Result<(), String> {
    if report.schema != REPORT_SCHEMA {
        return Err("unknown sect113r1 report schema".into());
    }
    if report.certificate_sha256 != report_digest(report)? {
        return Err("sect113r1 certificate digest mismatch".into());
    }
    let expected = run_audit()?;
    if *report != expected {
        return Err("sect113r1 report differs from exact replay".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn field_polynomial_passes_exact_rabin_test() {
        assert!(certify_field_polynomial());
    }

    #[test]
    fn exact_source_and_singular_certificates_pass() {
        let curve = BinaryCurve::sect113r1();
        let source = certify_source(&curve).unwrap();
        assert_eq!(source.identity.icv1, CURVE_ICV1);
        assert_eq!(source.identity.curve_id, CURVE_ID);
        assert_eq!(source.identity.curve_uid, CURVE_UID);
        assert!(source.n_primality.verified_prime);
        assert!(source.embedding_degree_exact);
        assert!(source.discriminant_fundamental);
        let singular = certify_singular_companion(&curve).unwrap();
        assert!(singular.demo_exact_recovery);
        assert!(!singular.full_width_large_component_executed);
    }

    #[test]
    fn degree_five_maps_and_duals_are_exact() {
        let curve = BinaryCurve::sect113r1();
        let certificate = certify_isogenies(&curve, &decimal(FULL_SCALAR)).unwrap();
        assert_eq!(certificate.rational_degree_five_kernel_count, 2);
        assert_eq!(certificate.frobenius_factor_degrees_q_q2_q4, [2, 4, 12]);
        assert!(certificate
            .representatives
            .iter()
            .all(|item| item.dual_kernel_certificate_exact
                && item.dual_composition_on_generator_is_times_five
                && item.dual_composition_on_planted_target_is_times_five
                && item.planted_target_pullback_to_source));
        let identities: Vec<_> = certificate
            .representatives
            .iter()
            .map(|item| {
                (
                    item.codomain_identity.icv1.as_str(),
                    item.codomain_identity.curve_id.as_str(),
                    item.codomain_identity.curve_uid.as_str(),
                )
            })
            .collect();
        assert_eq!(
            identities,
            vec![
                (
                    "ICV1:f2m-113-99967757:-122610772499221213:10384593717069655379671765157661406:0xb1df3d9c423aa919217735a6d6ba:unk:unk:r:fd54d0ebcd45",
                    "EC1N113Crbh921ab2cd913f",
                    "urn:ec-record:1:sha256:921ab2cd913f1d63106c1a55555cd869ae110bff6992dd23ee2522a4b97616c6",
                ),
                (
                    "ICV1:f2m-113-99967757:-122610772499221213:10384593717069655379671765157661406:0x17ef6a3b098891f6b4c529b62f0dd:unk:unk:r:5de3030fefcc",
                    "EC1N113Crbh2ee3581888f2",
                    "urn:ec-record:1:sha256:2ee3581888f2dcbc20ded6386aac140f0ced5d5b8ff92aab74887ecd19f19386",
                ),
            ]
        );
        assert!(certificate
            .representatives
            .iter()
            .all(|item| item.map_implementation == "binary-velu/v1"));
    }

    #[test]
    fn twist_order_is_derived_from_the_source_order() {
        let curve = BinaryCurve::sect113r1();
        let certificate = certify_twist(&curve).unwrap();
        assert!(certificate.derived_from_source_order);
        assert_eq!(certificate.order, "10384593717069655134450220159218980");
    }

    #[test]
    fn mutated_degree_five_kernel_is_rejected() {
        let curve = BinaryCurve::sect113r1();
        let psi5 = division_polynomial(&curve.b, 5, curve.m, &curve.irreducible);
        let (mut kernel, _, _) = degree_five_kernels(&curve).unwrap();
        kernel.coeffs[0].add_assign(&F2mElement::one(curve.m));
        assert!(!kernel_certificate_ok(&kernel, &psi5, &curve));
    }

    #[test]
    fn complete_report_replays() {
        let report = run_audit().unwrap();
        verify_report(&report).unwrap();
    }

    #[test]
    fn semantic_report_mutation_is_rejected_even_with_a_fresh_digest() {
        let mut report = run_audit().unwrap();
        report.class_weakness.verdict = "NO_WEAKNESS_FOUND_WITHIN_SCOPE".into();
        report.certificate_sha256 = report_digest(&report).unwrap();
        assert!(verify_report(&report).is_err());
    }

    #[test]
    fn unknown_report_fields_are_rejected_at_every_level() {
        let report = run_audit().unwrap();
        let mut top_level = serde_json::to_value(&report).unwrap();
        top_level
            .as_object_mut()
            .unwrap()
            .insert("unknown".into(), serde_json::Value::Bool(true));
        assert!(serde_json::from_value::<AuditReport>(top_level).is_err());

        let mut nested = serde_json::to_value(&report).unwrap();
        nested["source"]
            .as_object_mut()
            .unwrap()
            .insert("unknown".into(), serde_json::Value::Bool(true));
        assert!(serde_json::from_value::<AuditReport>(nested).is_err());
    }
}
