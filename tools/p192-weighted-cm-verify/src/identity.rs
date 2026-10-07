//! Independent replay of the frozen P-192 identity and maximal-order gates.
//!
//! This module deliberately implements its own affine group law and its own
//! Pocklington verifier.  It does not call the repository's curve or
//! primality implementations.

use std::{fs, path::Path};

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};

use crate::{
    encoding::{canonical_json, parse_canonical_json, sha256},
    params::{CURVE_UID, D_DEC, N_DEC, P_DEC, T_DEC},
    Result,
};

pub const EXPERIMENT_ID: &str = "EXP-SCURVE-1a8daf";
pub const PROTOCOL_VERSION: u64 = 2;
pub const IDENTITY_SCHEMA: &str = "p192-wcm-identity-v1";
pub const MAXIMAL_ORDER_SCHEMA: &str = "p192-wcm-maximal-order-v1";

pub const CURVE_NAME: &str = "NIST P-192 / secp192r1";
pub const CURVE_MODEL: &str = "short Weierstrass y^2 = x^3 + a*x + b over F_p";
pub const A_DEC: &str = "6277101735386680763835789423207666416083908700390324961276";
pub const B_DEC: &str = "2455155546008943817740293915197451784769108058161191238065";
pub const GX_HEX: &str = "188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012";
pub const GY_HEX: &str = "07192B95FFC8DA78631011ED6B24CDD573F977A11E794811";
pub const ICV1: &str = "icv1-fp192-t31607402316713927207482677199-52e4af59";
pub const EC1: &str = "EC1P192Cp192h5531c4a08bdb";
pub const C_DEC: &str = "14140398275856956083603613626459809774163796198179979683";
pub const FACTORIZATION: &str = "-5 * 11 * 31 * C";

const SEC2_STANDARD_ID: &str = "SEC2";
const SEC2_EDITION: &str = "2.0";
const SEC2_SECTION: &str = "2.2.2";
const SEC2_URL: &str = "https://www.secg.org/sec2-v2.pdf";
const SEC2_SHA256: &str = "87b8f3703364ed5b21ba8582e411cc0cbf477bcaa3f4f45e0d6580d1c00d9952";
const SEC2_BYTE_LENGTH: u64 = 306_784;
const SEC2_PARAMETER_NAME: &str = "secp192r1";

const IDENTITY_CHECK_IDS: [&str; 6] = [
    "serialized_tuple_matches_protocol",
    "identities_recomputed",
    "nonsingular",
    "generator_nonidentity",
    "generator_order_n",
    "exact_group_order_n",
];

const MAXIMAL_ORDER_CHECK_IDS: [&str; 17] = [
    "p_prime",
    "n_prime",
    "C_prime",
    "p_not_divide_t",
    "trace_absolute_value_below_p",
    "ordinary",
    "discriminant_identity",
    "group_order_identity",
    "discriminant_factorization",
    "D_squarefree",
    "D_congruent_1_mod_4",
    "D_fundamental",
    "frobenius_conductor_one",
    "endomorphism_order_maximal",
    "discriminants_equal",
    "unit_group_plus_minus_one",
    "integral_basis_multiplication",
];

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Verdict {
    pub check_id: String,
    pub status: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CurveTuple {
    pub name: String,
    pub model: String,
    pub p: String,
    pub a: String,
    pub b: String,
    pub generator_x_hex: String,
    pub generator_y_hex: String,
    pub n: String,
    pub cofactor: u64,
    pub trace_t: String,
    pub icv1: String,
    pub ec1: String,
    pub curve_uid: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IdentityDerived {
    pub group_order: String,
    pub frobenius_discriminant: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IdentityManifest {
    pub schema: String,
    pub experiment_id: String,
    pub protocol_version: u64,
    pub protocol_commit: String,
    pub source_commit: String,
    pub curve: CurveTuple,
    pub derived: IdentityDerived,
    pub checks: Vec<Verdict>,
    pub overall_status: String,
}

#[derive(Clone, Debug)]
pub struct VerifiedIdentity {
    pub manifest: IdentityManifest,
    pub canonical_bytes: Vec<u8>,
    pub curve_tuple_sha256: [u8; 32],
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OrderTuple {
    pub p: String,
    pub n: String,
    pub t: String,
    #[serde(rename = "D_pi")]
    pub d_pi: String,
    #[serde(rename = "D_K")]
    pub d_k: String,
    #[serde(rename = "D_End")]
    pub d_end: String,
    #[serde(rename = "D")]
    pub discriminant: String,
    #[serde(rename = "C")]
    pub remaining_prime: String,
    pub factorization: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PocklingtonFactor {
    pub prime: String,
    pub exponent: u64,
    pub proof: Box<PrimalityProof>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PocklingtonWitness {
    pub prime: String,
    pub base: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "method", deny_unknown_fields)]
pub enum PrimalityProof {
    #[serde(rename = "authoritative_standard")]
    AuthoritativeStandard {
        value: String,
        standard_id: String,
        edition: String,
        section: String,
        document_url: String,
        document_sha256: String,
        document_byte_length: u64,
        parameter_name: String,
    },
    #[serde(rename = "trial_division_u32")]
    TrialDivisionU32 { value: String },
    #[serde(rename = "pocklington_v1")]
    PocklingtonV1 {
        value: String,
        cofactor: String,
        factors: Vec<PocklingtonFactor>,
        witnesses: Vec<PocklingtonWitness>,
    },
}

impl PrimalityProof {
    pub fn value(&self) -> &str {
        match self {
            Self::AuthoritativeStandard { value, .. }
            | Self::TrialDivisionU32 { value }
            | Self::PocklingtonV1 { value, .. } => value,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimalityRecord {
    pub subject: String,
    pub value: String,
    pub proof: PrimalityProof,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaximalOrderManifest {
    pub schema: String,
    pub experiment_id: String,
    pub protocol_version: u64,
    pub protocol_commit: String,
    pub source_commit: String,
    pub curve_uid: String,
    pub order: OrderTuple,
    pub primality: Vec<PrimalityRecord>,
    pub checks: Vec<Verdict>,
    pub overall_status: String,
}

#[derive(Clone, Debug)]
pub struct VerifiedMaximalOrder {
    pub manifest: MaximalOrderManifest,
    pub canonical_bytes: Vec<u8>,
    pub cm_order_tuple_sha256: [u8; 32],
}

fn unsigned_decimal(text: &str, name: &str) -> Result<BigUint> {
    if text.is_empty()
        || text.starts_with('+')
        || text.starts_with('-')
        || (text.starts_with('0') && text != "0")
        || !text.bytes().all(|byte| byte.is_ascii_digit())
    {
        return Err(format!("{name} is not a canonical unsigned decimal"));
    }
    BigUint::parse_bytes(text.as_bytes(), 10).ok_or_else(|| format!("invalid {name}"))
}

fn positive_decimal(text: &str, name: &str) -> Result<BigUint> {
    let value = unsigned_decimal(text, name)?;
    if value.is_zero() {
        return Err(format!("{name} must be positive"));
    }
    Ok(value)
}

fn signed_decimal(text: &str, name: &str) -> Result<BigInt> {
    if text.is_empty()
        || text.starts_with('+')
        || text == "-0"
        || text.starts_with("-0")
        || (text.starts_with('0') && text != "0")
    {
        return Err(format!("{name} is not a canonical signed decimal"));
    }
    let unsigned = text.strip_prefix('-').unwrap_or(text);
    if unsigned.is_empty() || !unsigned.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(format!("invalid {name}"));
    }
    BigInt::parse_bytes(text.as_bytes(), 10).ok_or_else(|| format!("invalid {name}"))
}

fn is_lower_hex(text: &str, length: usize) -> bool {
    text.len() == length
        && text
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn validate_envelope(
    schema: &str,
    expected_schema: &str,
    experiment_id: &str,
    protocol_version: u64,
    protocol_commit: &str,
    source_commit: &str,
) -> Result<()> {
    if schema != expected_schema
        || experiment_id != EXPERIMENT_ID
        || protocol_version != PROTOCOL_VERSION
    {
        return Err(format!(
            "{expected_schema} envelope differs from the frozen protocol"
        ));
    }
    if !is_lower_hex(protocol_commit, 40) || !is_lower_hex(source_commit, 40) {
        return Err("protocol/source commit is not lowercase 40-hex".to_owned());
    }
    Ok(())
}

fn require_passes(checks: &[Verdict], expected: &[&str], overall: &str) -> Result<()> {
    if checks.len() != expected.len() {
        return Err("check count differs from the frozen schema".to_owned());
    }
    for (check, expected_id) in checks.iter().zip(expected) {
        if check.check_id != *expected_id || check.status != "PASS" {
            return Err(format!(
                "check {expected_id} is missing, reordered, or not PASS"
            ));
        }
    }
    if overall != "PASS" {
        return Err("overall_status is not PASS".to_owned());
    }
    Ok(())
}

fn mod_inverse(value: &BigUint, modulus: &BigUint) -> Result<BigUint> {
    let value = BigInt::from(value % modulus);
    let modulus_i = BigInt::from(modulus.clone());
    let extended = value.extended_gcd(&modulus_i);
    if !extended.gcd.is_one() {
        return Err("noninvertible field denominator".to_owned());
    }
    extended
        .x
        .mod_floor(&modulus_i)
        .to_biguint()
        .ok_or_else(|| "inverse residue is negative".to_owned())
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum AffinePoint {
    Infinity,
    Finite { x: BigUint, y: BigUint },
}

fn point_add(
    left: &AffinePoint,
    right: &AffinePoint,
    a: &BigUint,
    p: &BigUint,
) -> Result<AffinePoint> {
    let (x1, y1, x2, y2) = match (left, right) {
        (AffinePoint::Infinity, point) | (point, AffinePoint::Infinity) => {
            return Ok(point.clone());
        }
        (AffinePoint::Finite { x: x1, y: y1 }, AffinePoint::Finite { x: x2, y: y2 }) => {
            (x1, y1, x2, y2)
        }
    };

    if x1 == x2 && (y1 + y2) % p == BigUint::zero() {
        return Ok(AffinePoint::Infinity);
    }
    let slope = if x1 == x2 {
        if y1.is_zero() {
            return Ok(AffinePoint::Infinity);
        }
        let numerator = (BigUint::from(3u8) * x1 * x1 + a) % p;
        let denominator = (BigUint::from(2u8) * y1) % p;
        numerator * mod_inverse(&denominator, p)? % p
    } else {
        let numerator = (y2 + p - y1) % p;
        let denominator = (x2 + p - x1) % p;
        numerator * mod_inverse(&denominator, p)? % p
    };
    let x3 = (&slope * &slope + p + p - x1 - x2) % p;
    let y3 = (&slope * ((x1 + p - &x3) % p) + p - y1) % p;
    Ok(AffinePoint::Finite { x: x3, y: y3 })
}

fn scalar_mul(
    point: &AffinePoint,
    scalar: &BigUint,
    a: &BigUint,
    p: &BigUint,
) -> Result<AffinePoint> {
    let mut result = AffinePoint::Infinity;
    let mut addend = point.clone();
    for bit in 0..scalar.bits() {
        if scalar.bit(bit) {
            result = point_add(&result, &addend, a, p)?;
        }
        addend = point_add(&addend, &addend, a, p)?;
    }
    Ok(result)
}

fn curve_identity_values(curve: &CurveTuple) -> Result<(String, String, String)> {
    let p = positive_decimal(&curve.p, "curve.p")?;
    let a = positive_decimal(&curve.a, "curve.a")?;
    let b = positive_decimal(&curve.b, "curve.b")?;
    let n = positive_decimal(&curve.n, "curve.n")?;
    let trace = signed_decimal(&curve.trace_t, "curve.trace_t")?;
    let gx = BigUint::parse_bytes(curve.generator_x_hex.as_bytes(), 16)
        .ok_or_else(|| "invalid generator_x_hex".to_owned())?;
    let gy = BigUint::parse_bytes(curve.generator_y_hex.as_bytes(), 16)
        .ok_or_else(|| "invalid generator_y_hex".to_owned())?;

    let model_json = format!(
        "{{\"a\":\"{a}\",\"b\":\"{b}\",\"field\":\"fp-{p}\",\"form\":\"y^2=x^3+a*x+b\",\"p\":\"{p}\",\"v\":\"1\"}}"
    );
    let model_hash = hex::encode(sha256(model_json.as_bytes()));
    let icv1 = format!("icv1-fp{}-t{trace}-{}", p.bits(), &model_hash[..8]);

    // This is the exact canonical EC1 preimage.  Its very large
    // characteristic is emitted as an integer by the historical identity
    // format, so construct the bytes directly instead of passing it through
    // serde_json's finite numeric representation.
    let ec1_preimage = format!(
        concat!(
            "{{\"curve\":{{\"coefficients\":[0,0,0,{a},{b}],",
            "\"cofactor\":1,\"generator\":[\"0x{gx:x}\",\"0x{gy:x}\"],",
            "\"model\":\"short Weierstrass\",\"subgroup_order\":\"{n}\",",
            "\"target_group\":\"prime-order subgroup\"}},",
            "\"field\":{{\"characteristic\":{p},\"degree\":1,",
            "\"element_encoding\":\"hex integer modulo the characteristic\",",
            "\"representation\":\"prime\"}}}}"
        ),
        a = a,
        b = b,
        gx = gx,
        gy = gy,
        n = n,
        p = p,
    );
    let ec1_hash = hex::encode(sha256(ec1_preimage.as_bytes()));
    let ec1 = format!("EC1P{}Cp192h{}", p.bits(), &ec1_hash[..12]);
    let uid = format!("urn:ec-record:1:sha256:{ec1_hash}");
    Ok((icv1, ec1, uid))
}

fn verify_curve_arithmetic(curve: &CurveTuple) -> Result<()> {
    let p = positive_decimal(&curve.p, "curve.p")?;
    let a = positive_decimal(&curve.a, "curve.a")?;
    let b = positive_decimal(&curve.b, "curve.b")?;
    let n = positive_decimal(&curve.n, "curve.n")?;
    let gx = BigUint::parse_bytes(curve.generator_x_hex.as_bytes(), 16)
        .ok_or_else(|| "invalid generator_x_hex".to_owned())?;
    let gy = BigUint::parse_bytes(curve.generator_y_hex.as_bytes(), 16)
        .ok_or_else(|| "invalid generator_y_hex".to_owned())?;
    if gx >= p || gy >= p {
        return Err("generator coordinate is outside F_p".to_owned());
    }
    let discriminant = (BigUint::from(4u8) * a.modpow(&BigUint::from(3u8), &p)
        + BigUint::from(27u8) * &b * &b)
        % &p;
    if discriminant.is_zero() {
        return Err("P-192 serialization is singular".to_owned());
    }
    if gy.modpow(&BigUint::from(2u8), &p)
        != (gx.modpow(&BigUint::from(3u8), &p) + &a * &gx + &b) % &p
    {
        return Err("generator is not on the serialized curve".to_owned());
    }
    let generator = AffinePoint::Finite { x: gx, y: gy };
    if scalar_mul(&generator, &n, &a, &p)? != AffinePoint::Infinity {
        return Err("[n]G is not the identity".to_owned());
    }
    Ok(())
}

pub fn verify_identity_bytes(bytes: &[u8]) -> Result<VerifiedIdentity> {
    let value = parse_canonical_json(bytes)?;
    let manifest: IdentityManifest =
        serde_json::from_slice(bytes).map_err(|error| format!("identity schema: {error}"))?;
    validate_envelope(
        &manifest.schema,
        IDENTITY_SCHEMA,
        &manifest.experiment_id,
        manifest.protocol_version,
        &manifest.protocol_commit,
        &manifest.source_commit,
    )?;
    if manifest.curve
        != (CurveTuple {
            name: CURVE_NAME.to_owned(),
            model: CURVE_MODEL.to_owned(),
            p: P_DEC.to_owned(),
            a: A_DEC.to_owned(),
            b: B_DEC.to_owned(),
            generator_x_hex: GX_HEX.to_owned(),
            generator_y_hex: GY_HEX.to_owned(),
            n: N_DEC.to_owned(),
            cofactor: 1,
            trace_t: T_DEC.to_owned(),
            icv1: ICV1.to_owned(),
            ec1: EC1.to_owned(),
            curve_uid: CURVE_UID.to_owned(),
        })
    {
        return Err("identity curve tuple differs from the frozen P-192 tuple".to_owned());
    }
    if manifest.derived.group_order != N_DEC || manifest.derived.frobenius_discriminant != D_DEC {
        return Err("identity derived values differ from their recomputation".to_owned());
    }
    let p = positive_decimal(&manifest.curve.p, "curve.p")?;
    let n = positive_decimal(&manifest.curve.n, "curve.n")?;
    let t = signed_decimal(&manifest.curve.trace_t, "curve.trace_t")?;
    if BigInt::from(p.clone()) + 1 - &t != BigInt::from(n)
        || &t * &t - BigInt::from(p) * 4 != signed_decimal(D_DEC, "D")?
    {
        return Err("identity order/discriminant recomputation failed".to_owned());
    }
    let identities = curve_identity_values(&manifest.curve)?;
    if identities.0 != manifest.curve.icv1
        || identities.1 != manifest.curve.ec1
        || identities.2 != manifest.curve.curve_uid
    {
        return Err("curve identity recomputation disagrees with identity.json".to_owned());
    }
    verify_curve_arithmetic(&manifest.curve)?;
    require_passes(
        &manifest.checks,
        &IDENTITY_CHECK_IDS,
        &manifest.overall_status,
    )?;
    let curve_value = value
        .get("curve")
        .ok_or_else(|| "identity curve value is absent".to_owned())?;
    let mut projection = b"P192-WCM-CURVE-TUPLE-v1\0".to_vec();
    projection.extend_from_slice(&canonical_json(curve_value)?);
    Ok(VerifiedIdentity {
        manifest,
        canonical_bytes: bytes.to_vec(),
        curve_tuple_sha256: sha256(&projection),
    })
}

pub fn verify_identity(source_dir: &Path) -> Result<VerifiedIdentity> {
    let path = source_dir.join("identity.json");
    let bytes = fs::read(&path).map_err(|error| format!("read {}: {error}", path.display()))?;
    verify_identity_bytes(&bytes)
}

fn replay_trial_division(value_text: &str) -> Result<()> {
    let value = positive_decimal(value_text, "trial-division value")?;
    let value = value
        .to_u32()
        .ok_or_else(|| "trial_division_u32 value exceeds u32".to_owned())?;
    if value < 2 {
        return Err("trial-division value is not prime".to_owned());
    }
    let mut divisor = 2u32;
    while divisor <= value / divisor {
        if value % divisor == 0 {
            return Err(format!("trial-division value has divisor {divisor}"));
        }
        divisor += if divisor == 2 { 1 } else { 2 };
    }
    Ok(())
}

fn replay_authoritative(proof: &PrimalityProof, subject: &str) -> Result<()> {
    let PrimalityProof::AuthoritativeStandard {
        value,
        standard_id,
        edition,
        section,
        document_url,
        document_sha256,
        document_byte_length,
        parameter_name,
    } = proof
    else {
        return Err("internal authoritative-proof dispatch error".to_owned());
    };
    if !matches!(subject, "p" | "n") {
        return Err("authoritative_standard is permitted only for p and n".to_owned());
    }
    if standard_id != SEC2_STANDARD_ID
        || edition != SEC2_EDITION
        || section != SEC2_SECTION
        || document_url != SEC2_URL
        || document_sha256 != SEC2_SHA256
        || *document_byte_length != SEC2_BYTE_LENGTH
        || parameter_name != SEC2_PARAMETER_NAME
    {
        return Err("authoritative standard metadata differs from pinned SEC 2 v2.0".to_owned());
    }
    positive_decimal(value, "authoritative value")?;
    // The role-1 interface freezes this exact content-addressed SEC 2
    // provenance tuple.  Section 2.2.2 names secp192r1 and binds p and n to
    // the independently compiled constants below.  Runtime URLs or
    // producer-supplied bytes cannot change this projection.
    let table_value = match subject {
        "p" => P_DEC,
        "n" => N_DEC,
        _ => return Err("authoritative_standard subject is not p or n".to_owned()),
    };
    if value != table_value {
        return Err("authoritative standard named value mismatch".to_owned());
    }
    Ok(())
}

fn replay_proof_inner(proof: &PrimalityProof, subject: &str, depth: usize) -> Result<()> {
    if depth > 128 {
        return Err("primality-proof recursion exceeds 128 levels".to_owned());
    }
    match proof {
        PrimalityProof::AuthoritativeStandard { .. } => replay_authoritative(proof, subject),
        PrimalityProof::TrialDivisionU32 { value } => replay_trial_division(value),
        PrimalityProof::PocklingtonV1 {
            value,
            cofactor,
            factors,
            witnesses,
        } => {
            let candidate = positive_decimal(value, "Pocklington value")?;
            let cofactor = positive_decimal(cofactor, "Pocklington cofactor")?;
            if candidate < BigUint::from(3u8) || candidate.is_even() {
                return Err(
                    "Pocklington candidate must be an odd integer at least three".to_owned(),
                );
            }
            if factors.is_empty() || factors.len() != witnesses.len() {
                return Err(
                    "Pocklington factors/witnesses are empty or differ in length".to_owned(),
                );
            }
            let mut factored_part = BigUint::one();
            let mut previous = BigUint::zero();
            for (factor, witness) in factors.iter().zip(witnesses) {
                let prime = positive_decimal(&factor.prime, "Pocklington factor prime")?;
                if prime <= previous {
                    return Err("Pocklington factor primes are not strictly increasing".to_owned());
                }
                if factor.exponent == 0 {
                    return Err("Pocklington factor exponent is zero".to_owned());
                }
                if factor.proof.value() != factor.prime {
                    return Err("recursive proof value differs from factor prime".to_owned());
                }
                if matches!(*factor.proof, PrimalityProof::AuthoritativeStandard { .. }) {
                    return Err(
                        "authoritative_standard is forbidden for recursive factors".to_owned()
                    );
                }
                replay_proof_inner(&factor.proof, "recursive_factor", depth + 1)?;
                if witness.prime != factor.prime {
                    return Err("Pocklington witness order/prime mismatch".to_owned());
                }
                let exponent = u32::try_from(factor.exponent)
                    .map_err(|_| "Pocklington exponent exceeds u32".to_owned())?;
                factored_part *= prime.pow(exponent);

                let base = unsigned_decimal(&witness.base, "Pocklington witness base")?;
                if base <= BigUint::one() || base >= candidate {
                    return Err("Pocklington witness base is outside 1 < a < value".to_owned());
                }
                let candidate_minus_one = &candidate - 1u8;
                if base.modpow(&candidate_minus_one, &candidate) != BigUint::one() {
                    return Err("Pocklington Fermat congruence failed".to_owned());
                }
                let residue = base.modpow(&(&candidate_minus_one / &prime), &candidate);
                let difference = if residue.is_zero() {
                    &candidate - 1u8
                } else {
                    residue - 1u8
                };
                if difference.gcd(&candidate) != BigUint::one() {
                    return Err("Pocklington gcd condition failed".to_owned());
                }
                previous = prime;
            }
            if &cofactor * &factored_part != &candidate - 1u8 {
                return Err("Pocklington value-1 factorization is not exact".to_owned());
            }
            if &factored_part * &factored_part <= candidate {
                return Err("Pocklington factored part does not exceed sqrt(value)".to_owned());
            }
            Ok(())
        }
    }
}

pub fn replay_primality_proof(
    proof: &PrimalityProof,
    expected_value: &str,
    subject: &str,
) -> Result<()> {
    positive_decimal(expected_value, "expected prime")?;
    if proof.value() != expected_value {
        return Err(format!("{subject} proof value differs from its subject"));
    }
    replay_proof_inner(proof, subject, 0)
}

fn proof_method_allowed(subject: &str, proof: &PrimalityProof) -> bool {
    matches!(
        (subject, proof),
        ("p" | "n", PrimalityProof::AuthoritativeStandard { .. })
            | ("C", PrimalityProof::PocklingtonV1 { .. })
    )
}

fn verify_order_arithmetic(order: &OrderTuple) -> Result<()> {
    let p = positive_decimal(&order.p, "order.p")?;
    let n = positive_decimal(&order.n, "order.n")?;
    let t = positive_decimal(&order.t, "order.t")?;
    let d_pi = signed_decimal(&order.d_pi, "order.D_pi")?;
    let d_k = signed_decimal(&order.d_k, "order.D_K")?;
    let d_end = signed_decimal(&order.d_end, "order.D_End")?;
    let discriminant = signed_decimal(&order.discriminant, "order.D")?;
    let remaining_prime = positive_decimal(&order.remaining_prime, "order.C")?;

    if BigInt::from(&t * &t) - BigInt::from(&p * 4u8) != discriminant {
        return Err("D != t^2-4p".to_owned());
    }
    if &p + 1u8 - &t != n {
        return Err("n != p+1-t".to_owned());
    }
    if (&t % &p).is_zero() || t.is_zero() || t >= p {
        return Err("ordinary trace gate failed".to_owned());
    }
    if discriminant != -BigInt::from(BigUint::from(5u8) * 11u8 * 31u8 * &remaining_prime) {
        return Err("D != -5*11*31*C".to_owned());
    }
    if (&remaining_prime % 5u8).is_zero()
        || (&remaining_prime % 11u8).is_zero()
        || (&remaining_prime % 31u8).is_zero()
    {
        return Err("D factorization is not squarefree".to_owned());
    }
    if discriminant.mod_floor(&BigInt::from(4u8)) != BigInt::one() {
        return Err("D is not congruent to one modulo four".to_owned());
    }
    if d_pi != discriminant
        || d_k != discriminant
        || d_end != discriminant
        || d_pi != d_k
        || d_k != d_end
    {
        return Err("CM discriminants are not all equal".to_owned());
    }
    if discriminant >= BigInt::from(-4) {
        return Err("D does not imply unit group {+1,-1}".to_owned());
    }

    // Independently replay pi*pi = -p + t*pi in the integral basis (1,pi).
    let pi_squared = (-BigInt::from(p.clone()), BigInt::from(t.clone()));
    if pi_squared.0 + BigInt::from(p) != BigInt::zero() || pi_squared.1 != BigInt::from(t) {
        return Err("integral-basis multiplication replay failed".to_owned());
    }
    Ok(())
}

pub fn verify_maximal_order_bytes(bytes: &[u8]) -> Result<VerifiedMaximalOrder> {
    let value = parse_canonical_json(bytes)?;
    let manifest: MaximalOrderManifest =
        serde_json::from_slice(bytes).map_err(|error| format!("maximal-order schema: {error}"))?;
    validate_envelope(
        &manifest.schema,
        MAXIMAL_ORDER_SCHEMA,
        &manifest.experiment_id,
        manifest.protocol_version,
        &manifest.protocol_commit,
        &manifest.source_commit,
    )?;
    if manifest.curve_uid != CURVE_UID {
        return Err("maximal-order curve_uid differs from P-192".to_owned());
    }
    let expected_order = OrderTuple {
        p: P_DEC.to_owned(),
        n: N_DEC.to_owned(),
        t: T_DEC.to_owned(),
        d_pi: D_DEC.to_owned(),
        d_k: D_DEC.to_owned(),
        d_end: D_DEC.to_owned(),
        discriminant: D_DEC.to_owned(),
        remaining_prime: C_DEC.to_owned(),
        factorization: FACTORIZATION.to_owned(),
    };
    if manifest.order != expected_order {
        return Err("maximal-order tuple differs from the frozen CM order".to_owned());
    }
    verify_order_arithmetic(&manifest.order)?;
    if manifest.primality.len() != 3 {
        return Err("maximal-order primality array must contain p,n,C".to_owned());
    }
    for ((record, expected_subject), expected_value) in manifest
        .primality
        .iter()
        .zip(["p", "n", "C"])
        .zip([P_DEC, N_DEC, C_DEC])
    {
        if record.subject != expected_subject || record.value != expected_value {
            return Err("primality subjects/values are missing or reordered".to_owned());
        }
        if !proof_method_allowed(expected_subject, &record.proof) {
            return Err(format!(
                "{expected_subject} uses a primality-proof method forbidden by the frozen schema"
            ));
        }
        replay_primality_proof(&record.proof, expected_value, expected_subject)?;
    }
    require_passes(
        &manifest.checks,
        &MAXIMAL_ORDER_CHECK_IDS,
        &manifest.overall_status,
    )?;
    let order_value = value
        .get("order")
        .ok_or_else(|| "maximal-order order value is absent".to_owned())?;
    let mut projection = b"P192-WCM-CM-ORDER-TUPLE-v1\0".to_vec();
    projection.extend_from_slice(&canonical_json(order_value)?);
    Ok(VerifiedMaximalOrder {
        manifest,
        canonical_bytes: bytes.to_vec(),
        cm_order_tuple_sha256: sha256(&projection),
    })
}

pub fn verify_maximal_order(source_dir: &Path) -> Result<VerifiedMaximalOrder> {
    let path = source_dir.join("maximal-order.json");
    let bytes = fs::read(&path).map_err(|error| format!("read {}: {error}", path.display()))?;
    verify_maximal_order_bytes(&bytes)
}

pub fn verify_identity_and_maximal_order(
    source_dir: &Path,
) -> Result<(VerifiedIdentity, VerifiedMaximalOrder)> {
    let identity_path = source_dir.join("identity.json");
    let maximal_path = source_dir.join("maximal-order.json");
    let identity_bytes = fs::read(&identity_path)
        .map_err(|error| format!("read {}: {error}", identity_path.display()))?;
    let maximal_bytes = fs::read(&maximal_path)
        .map_err(|error| format!("read {}: {error}", maximal_path.display()))?;
    verify_identity_and_maximal_order_bytes(&identity_bytes, &maximal_bytes)
}

pub fn verify_identity_and_maximal_order_bytes(
    identity_bytes: &[u8],
    maximal_order_bytes: &[u8],
) -> Result<(VerifiedIdentity, VerifiedMaximalOrder)> {
    let identity = verify_identity_bytes(identity_bytes)?;
    let maximal = verify_maximal_order_bytes(maximal_order_bytes)?;
    if identity.manifest.protocol_commit != maximal.manifest.protocol_commit
        || identity.manifest.source_commit != maximal.manifest.source_commit
        || identity.manifest.curve.curve_uid != maximal.manifest.curve_uid
        || identity.manifest.curve.p != maximal.manifest.order.p
        || identity.manifest.curve.n != maximal.manifest.order.n
        || identity.manifest.curve.trace_t != maximal.manifest.order.t
        || identity.manifest.derived.frobenius_discriminant != maximal.manifest.order.discriminant
    {
        return Err("identity/maximal-order cross-binding failed".to_owned());
    }

    // The replayed proof establishes n prime.  Since G is nonidentity and
    // [n]G=O, its order is n.  Hasse then makes n the exact group order:
    // no other positive multiple of n lies in [p+1-2sqrt(p),p+1+2sqrt(p)].
    let p = positive_decimal(P_DEC, "p")?;
    let n = positive_decimal(N_DEC, "n")?;
    let center = &p + 1u8;
    let lower_delta = &center - &n;
    let twice_n = &n * 2u8;
    if &lower_delta * &lower_delta > &p * 4u8
        || twice_n <= center
        || (&twice_n - &center) * (&twice_n - &center) <= &p * 4u8
        || &center * &center <= &p * 4u8
    {
        return Err("Hasse exact-order implication failed".to_owned());
    }
    Ok((identity, maximal))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn trial(value: &str) -> PrimalityProof {
        PrimalityProof::TrialDivisionU32 {
            value: value.to_owned(),
        }
    }

    fn frozen_curve() -> CurveTuple {
        CurveTuple {
            name: CURVE_NAME.to_owned(),
            model: CURVE_MODEL.to_owned(),
            p: P_DEC.to_owned(),
            a: A_DEC.to_owned(),
            b: B_DEC.to_owned(),
            generator_x_hex: GX_HEX.to_owned(),
            generator_y_hex: GY_HEX.to_owned(),
            n: N_DEC.to_owned(),
            cofactor: 1,
            trace_t: T_DEC.to_owned(),
            icv1: ICV1.to_owned(),
            ec1: EC1.to_owned(),
            curve_uid: CURVE_UID.to_owned(),
        }
    }

    #[test]
    fn replays_trial_division_and_rejects_composite() {
        replay_primality_proof(&trial("65521"), "65521", "recursive_factor").unwrap();
        assert!(replay_primality_proof(&trial("65535"), "65535", "x").is_err());
        assert!(replay_primality_proof(&trial("1"), "1", "x").is_err());
    }

    #[test]
    fn subject_proof_methods_are_exact() {
        let authoritative = PrimalityProof::AuthoritativeStandard {
            value: P_DEC.to_owned(),
            standard_id: SEC2_STANDARD_ID.to_owned(),
            edition: SEC2_EDITION.to_owned(),
            section: SEC2_SECTION.to_owned(),
            document_url: SEC2_URL.to_owned(),
            document_sha256: SEC2_SHA256.to_owned(),
            document_byte_length: SEC2_BYTE_LENGTH,
            parameter_name: SEC2_PARAMETER_NAME.to_owned(),
        };
        let pocklington = PrimalityProof::PocklingtonV1 {
            value: "13".to_owned(),
            cofactor: "1".to_owned(),
            factors: Vec::new(),
            witnesses: Vec::new(),
        };
        assert!(proof_method_allowed("p", &authoritative));
        assert!(proof_method_allowed("n", &authoritative));
        assert!(!proof_method_allowed("C", &authoritative));
        assert!(proof_method_allowed("C", &pocklington));
        assert!(!proof_method_allowed("p", &pocklington));
        assert!(!proof_method_allowed("C", &trial("13")));
    }

    #[test]
    fn replays_recursive_pocklington_and_rejects_mutation() {
        let proof = PrimalityProof::PocklingtonV1 {
            value: "13".to_owned(),
            cofactor: "1".to_owned(),
            factors: vec![
                PocklingtonFactor {
                    prime: "2".to_owned(),
                    exponent: 2,
                    proof: Box::new(trial("2")),
                },
                PocklingtonFactor {
                    prime: "3".to_owned(),
                    exponent: 1,
                    proof: Box::new(trial("3")),
                },
            ],
            witnesses: vec![
                PocklingtonWitness {
                    prime: "2".to_owned(),
                    base: "2".to_owned(),
                },
                PocklingtonWitness {
                    prime: "3".to_owned(),
                    base: "2".to_owned(),
                },
            ],
        };
        replay_primality_proof(&proof, "13", "C").unwrap();
        let mut mutated = proof;
        let PrimalityProof::PocklingtonV1 { witnesses, .. } = &mut mutated else {
            unreachable!();
        };
        witnesses[1].base = "1".to_owned();
        assert!(replay_primality_proof(&mutated, "13", "C").is_err());
    }

    #[test]
    fn recomputes_all_three_curve_identities() {
        let curve = frozen_curve();
        assert_eq!(
            curve_identity_values(&curve).unwrap(),
            (ICV1.to_owned(), EC1.to_owned(), CURVE_UID.to_owned())
        );
        verify_curve_arithmetic(&curve).unwrap();
    }

    #[test]
    fn identity_parser_is_canonical_and_schema_strict() {
        let manifest = IdentityManifest {
            schema: IDENTITY_SCHEMA.to_owned(),
            experiment_id: EXPERIMENT_ID.to_owned(),
            protocol_version: PROTOCOL_VERSION,
            protocol_commit: "a".repeat(40),
            source_commit: "b".repeat(40),
            curve: frozen_curve(),
            derived: IdentityDerived {
                group_order: N_DEC.to_owned(),
                frobenius_discriminant: D_DEC.to_owned(),
            },
            checks: IDENTITY_CHECK_IDS
                .iter()
                .map(|check_id| Verdict {
                    check_id: (*check_id).to_owned(),
                    status: "PASS".to_owned(),
                })
                .collect(),
            overall_status: "PASS".to_owned(),
        };
        let value = serde_json::to_value(manifest).unwrap();
        let bytes = canonical_json(&value).unwrap();
        verify_identity_bytes(&bytes).unwrap();

        let mut with_unknown = value;
        with_unknown
            .as_object_mut()
            .unwrap()
            .insert("unexpected".to_owned(), serde_json::Value::Bool(true));
        let unknown_bytes = canonical_json(&with_unknown).unwrap();
        assert!(verify_identity_bytes(&unknown_bytes).is_err());

        let mut newline = bytes;
        newline.push(b'\n');
        assert!(verify_identity_bytes(&newline).is_err());
    }

    #[test]
    fn replays_frozen_cm_order_arithmetic() {
        verify_order_arithmetic(&OrderTuple {
            p: P_DEC.to_owned(),
            n: N_DEC.to_owned(),
            t: T_DEC.to_owned(),
            d_pi: D_DEC.to_owned(),
            d_k: D_DEC.to_owned(),
            d_end: D_DEC.to_owned(),
            discriminant: D_DEC.to_owned(),
            remaining_prime: C_DEC.to_owned(),
            factorization: FACTORIZATION.to_owned(),
        })
        .unwrap();
    }

    #[test]
    fn replays_exact_authoritative_metadata_tuple() {
        let proof = PrimalityProof::AuthoritativeStandard {
            value: P_DEC.to_owned(),
            standard_id: "SEC2".to_owned(),
            edition: "2.0".to_owned(),
            section: "2.2.2".to_owned(),
            document_url: "https://www.secg.org/sec2-v2.pdf".to_owned(),
            document_sha256: "87b8f3703364ed5b21ba8582e411cc0cbf477bcaa3f4f45e0d6580d1c00d9952"
                .to_owned(),
            document_byte_length: 306_784,
            parameter_name: "secp192r1".to_owned(),
        };
        replay_primality_proof(&proof, P_DEC, "p").unwrap();
    }
}
