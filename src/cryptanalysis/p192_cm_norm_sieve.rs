//! Exact certificate core for a CM norm-lattice relation search.
//!
//! The existing P-192 experiment enumerates exponent vectors under a bound on
//! the *multiplicative* isogeny degree.  That boundary is useful, but it misses
//! a different regime: a principal ideal can have enormous norm while its
//! factorisation gives a short, additively cheap chain of small-degree maps.
//!
//! This module starts from a principal element
//!
//! ```text
//! alpha = u + v*pi,
//! N(alpha) = u^2 + t*u*v + p*v^2,
//! ```
//!
//! factors its norm over an oriented split/ramified prime-ideal base, and
//! certifies the resulting class relation with binary quadratic forms.  It is
//! deliberately only an algebra/candidate layer: a relation is not an ECDLP
//! weakness until explicit maps, subgroup transport, evaluation costs, and a
//! fully charged payoff all pass separately.

use crate::isogeny::class_group::BinaryQuadraticForm;

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const SCHEMA_FACTOR_BASE: &str = "cm.norm_sieve.factor_base/v1";
pub const SCHEMA_CERTIFICATE: &str = "cm.norm_sieve.certificate/v1";

fn mod_u64(value: &BigInt, modulus: u64) -> u64 {
    value
        .mod_floor(&BigInt::from(modulus))
        .to_u64_digits()
        .1
        .first()
        .copied()
        .unwrap_or(0)
}

fn mul_mod(left: u64, right: u64, modulus: u64) -> u64 {
    ((left as u128 * right as u128) % modulus as u128) as u64
}

fn pow_mod(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut accumulator = 1 % modulus;
    while exponent != 0 {
        if exponent & 1 == 1 {
            accumulator = mul_mod(accumulator, base, modulus);
        }
        exponent >>= 1;
        if exponent != 0 {
            base = mul_mod(base, base, modulus);
        }
    }
    accumulator
}

/// Tonelli--Shanks for an odd prime.  Returns the smaller square root.
fn sqrt_mod_prime(value: u64, prime: u64) -> Option<u64> {
    debug_assert!(prime > 2 && prime % 2 == 1);
    let value = value % prime;
    if value == 0 {
        return Some(0);
    }
    if pow_mod(value, (prime - 1) / 2, prime) != 1 {
        return None;
    }
    if prime % 4 == 3 {
        let root = pow_mod(value, (prime + 1) / 4, prime);
        return Some(root.min(prime - root));
    }

    let mut odd = prime - 1;
    let mut power = 0u32;
    while odd % 2 == 0 {
        odd /= 2;
        power += 1;
    }
    let non_residue =
        (2..prime).find(|candidate| pow_mod(*candidate, (prime - 1) / 2, prime) == prime - 1)?;
    let mut c = pow_mod(non_residue, odd, prime);
    let mut x = pow_mod(value, (odd + 1) / 2, prime);
    let mut t = pow_mod(value, odd, prime);
    let mut m = power;
    while t != 1 {
        let mut i = 1u32;
        let mut probe = mul_mod(t, t, prime);
        while probe != 1 {
            probe = mul_mod(probe, probe, prime);
            i += 1;
            if i == m {
                return None;
            }
        }
        let b = pow_mod(c, 1u64 << (m - i - 1), prime);
        x = mul_mod(x, b, prime);
        let b2 = mul_mod(b, b, prime);
        t = mul_mod(t, b2, prime);
        c = b2;
        m = i;
    }
    Some(x.min(prime - x))
}

fn primes_through(limit: u64) -> Result<Vec<u64>, String> {
    let length = usize::try_from(limit)
        .map_err(|_| "factor-base limit does not fit usize".to_owned())?
        .checked_add(1)
        .ok_or_else(|| "factor-base limit overflow".to_owned())?;
    let mut prime = vec![true; length];
    if !prime.is_empty() {
        prime[0] = false;
    }
    if prime.len() > 1 {
        prime[1] = false;
    }
    let mut q = 2usize;
    while q <= (length.saturating_sub(1)) / q {
        if prime[q] {
            let mut multiple = q * q;
            while multiple < length {
                prime[multiple] = false;
                multiple += q;
            }
        }
        q += 1;
    }
    Ok(prime
        .into_iter()
        .enumerate()
        .filter_map(|(value, is_prime)| is_prime.then_some(value as u64))
        .collect())
}

fn characteristic_roots(ell: u64, t: &BigInt, p: &BigInt, d: &BigInt) -> Vec<u64> {
    if ell == 2 {
        return (0..2)
            .filter(|root| {
                let r = *root as u128;
                let tm = mod_u64(t, 2) as u128;
                let pm = mod_u64(p, 2) as u128;
                (r * r + 2 - (tm * r) % 2 + pm) % 2 == 0
            })
            .collect();
    }
    let dm = mod_u64(d, ell);
    let Some(sqrt_d) = sqrt_mod_prime(dm, ell) else {
        return Vec::new();
    };
    let tm = mod_u64(t, ell);
    let inv_two = ell.div_ceil(2);
    let mut roots = vec![
        mul_mod((tm + sqrt_d) % ell, inv_two, ell),
        mul_mod((tm + ell - sqrt_d) % ell, inv_two, ell),
    ];
    roots.sort_unstable();
    roots.dedup();
    roots.retain(|root| {
        let r = *root as u128;
        let modulus = ell as u128;
        let tm = tm as u128;
        let pm = mod_u64(p, ell) as u128;
        (r * r + modulus - (tm * r) % modulus + pm) % modulus == 0
    });
    roots
}

fn form_for_root(
    disc: &BigInt,
    trace: &BigInt,
    field_norm: &BigInt,
    ell: u64,
    root: u64,
) -> Result<BinaryQuadraticForm, String> {
    let a = BigInt::from(ell);
    let r = BigInt::from(root);
    if (&r * &r - trace * &r + field_norm).mod_floor(&a) != BigInt::zero() {
        return Err(format!(
            "root {root} is not a characteristic root modulo {ell}"
        ));
    }
    // I=(ell,pi-r), sqrt(D)=2*pi-t, and
    // (-b+sqrt(D))/2 == pi-r (mod I), so b=2r-t (mod 2ell).
    let two_a = &a * 2;
    let mut b = (&r * BigInt::from(2u8) - trace).mod_floor(&two_a);
    if b > a {
        b -= &two_a;
    }
    let numerator = &b * &b - disc;
    let denominator = &a * 4;
    if !numerator.mod_floor(&denominator).is_zero() {
        return Err(format!(
            "root/form parity mismatch at ell={ell}, root={root}"
        ));
    }
    let c = numerator / denominator;
    BinaryQuadraticForm::new(a, b, c, disc)
        .map(|form| form.reduce())
        .ok_or_else(|| format!("invalid prime-ideal form at ell={ell}, root={root}"))
}

fn form_strings(form: &BinaryQuadraticForm) -> [String; 3] {
    [form.a.to_string(), form.b.to_string(), form.c.to_string()]
}

fn strings_form(values: &[String; 3], disc: &BigInt) -> Result<BinaryQuadraticForm, String> {
    BinaryQuadraticForm::new(
        BigInt::parse_bytes(values[0].as_bytes(), 10)
            .ok_or_else(|| "invalid form coefficient a".to_owned())?,
        BigInt::parse_bytes(values[1].as_bytes(), 10)
            .ok_or_else(|| "invalid form coefficient b".to_owned())?,
        BigInt::parse_bytes(values[2].as_bytes(), 10)
            .ok_or_else(|| "invalid form coefficient c".to_owned())?,
        disc,
    )
    .map(|form| form.reduce())
    .ok_or_else(|| "invalid factor-base form".to_owned())
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct NormSievePrime {
    pub ell: u64,
    pub kind: String,
    pub positive_root: u64,
    pub negative_root: u64,
    pub positive_form: [String; 3],
    pub negative_form: [String; 3],
    /// Exact integer weight used only to order candidates.  The label on the
    /// enclosing factor base says whether this is a proxy or a measurement.
    pub edge_weight: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct NormSieveFactorBase {
    pub schema: String,
    pub discriminant: String,
    pub trace: String,
    pub field_norm: String,
    pub max_ell: u64,
    pub weight_model: String,
    pub primes: Vec<NormSievePrime>,
}

/// Construct all split and ramified prime ideals through `max_ell`.
///
/// `weights` is intentionally explicit.  Missing values are rejected rather
/// than silently replacing unavailable measured map costs with degree or
/// wall-time proxies.
pub fn build_factor_base(
    disc: &BigInt,
    trace: &BigInt,
    field_norm: &BigInt,
    max_ell: u64,
    weight_model: &str,
    weights: &BTreeMap<u64, u64>,
) -> Result<NormSieveFactorBase, String> {
    if trace * trace - field_norm * 4 != *disc || !disc.is_negative() {
        return Err("D must equal t^2-4p and be negative".to_owned());
    }
    let mut records = Vec::new();
    for ell in primes_through(max_ell)? {
        let roots = characteristic_roots(ell, trace, field_norm, disc);
        if roots.is_empty() {
            continue;
        }
        if roots.len() > 2 {
            return Err(format!("more than two characteristic roots at ell={ell}"));
        }
        let edge_weight = *weights
            .get(&ell)
            .ok_or_else(|| format!("missing explicit edge weight for ell={ell}"))?;
        if edge_weight == 0 {
            return Err(format!("zero edge weight is invalid at ell={ell}"));
        }
        let positive_root = roots[0];
        let negative_root = *roots.last().unwrap();
        let positive = form_for_root(disc, trace, field_norm, ell, positive_root)?;
        let negative = form_for_root(disc, trace, field_norm, ell, negative_root)?;
        let kind = if roots.len() == 1 {
            if positive != positive.inverse() {
                return Err(format!("ramified class is not self-inverse at ell={ell}"));
            }
            "ramified"
        } else {
            if positive.inverse() != negative {
                return Err(format!("split orientations are not inverse at ell={ell}"));
            }
            "split"
        };
        records.push(NormSievePrime {
            ell,
            kind: kind.to_owned(),
            positive_root,
            negative_root,
            positive_form: form_strings(&positive),
            negative_form: form_strings(&negative),
            edge_weight,
        });
    }
    Ok(NormSieveFactorBase {
        schema: SCHEMA_FACTOR_BASE.to_owned(),
        discriminant: disc.to_string(),
        trace: trace.to_string(),
        field_norm: field_norm.to_string(),
        max_ell,
        weight_model: weight_model.to_owned(),
        primes: records,
    })
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PrimeIdealFactor {
    pub ell: u64,
    pub exponent: u32,
    pub orientation: String,
    pub root: u64,
    pub edge_weight: u64,
    pub charged_weight: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct NormRelationCertificate {
    pub schema: String,
    pub u: String,
    pub v: String,
    pub primitive: bool,
    pub alpha_norm: String,
    pub factors: Vec<PrimeIdealFactor>,
    pub exponent_vector: Vec<i32>,
    pub residual_norm_factor: String,
    pub complete_factorization: bool,
    pub additive_chain_weight: u64,
    pub reduced_known_ideal_product: [String; 3],
    pub principal_relation_verified: bool,
    /// On an F_p-rational subgroup Frobenius acts as one, so u+v*pi acts as
    /// the scalar u+v.  This field is only that algebraic candidate; it does
    /// not certify a map or subgroup transport.
    pub lambda_mod_subgroup_order: String,
    pub map_status: String,
    pub payoff_status: String,
}

fn parse_decimal(value: &str, name: &str) -> Result<BigInt, String> {
    BigInt::parse_bytes(value.as_bytes(), 10).ok_or_else(|| format!("invalid decimal {name}"))
}

/// Factor and certify one primitive non-scalar CM element over the supplied
/// oriented factor base.
pub fn certify_norm_candidate(
    u: &BigInt,
    v: &BigInt,
    subgroup_order: &BigInt,
    factor_base: &NormSieveFactorBase,
) -> Result<NormRelationCertificate, String> {
    if factor_base.schema != SCHEMA_FACTOR_BASE {
        return Err("unsupported factor-base schema".to_owned());
    }
    let disc = parse_decimal(&factor_base.discriminant, "discriminant")?;
    let trace = parse_decimal(&factor_base.trace, "trace")?;
    let field_norm = parse_decimal(&factor_base.field_norm, "field norm")?;
    if &trace * &trace - &field_norm * 4 != disc {
        return Err("factor-base D/t/p identity failed".to_owned());
    }
    if subgroup_order <= &BigInt::one() {
        return Err("subgroup order must exceed one".to_owned());
    }
    if v.is_zero() {
        return Err("v=0 is a scalar, not a non-scalar relation candidate".to_owned());
    }
    let primitive = u.gcd(v).abs().is_one();
    if !primitive {
        return Err("candidate u,v must be primitive; remove scalar content first".to_owned());
    }
    let norm = u * u + &trace * u * v + &field_norm * v * v;
    if norm <= BigInt::zero() {
        return Err("imaginary-quadratic norm must be positive".to_owned());
    }
    let mut residual = norm
        .to_biguint()
        .ok_or_else(|| "positive norm conversion failed".to_owned())?;
    let mut factors = Vec::new();
    let mut vector = Vec::with_capacity(factor_base.primes.len());
    let mut additive_weight = 0u64;
    let mut product = BinaryQuadraticForm::principal(&disc)
        .ok_or_else(|| "principal form does not exist".to_owned())?;

    let mut previous = 0u64;
    for record in &factor_base.primes {
        if record.ell <= previous {
            return Err("factor-base primes are not strictly increasing".to_owned());
        }
        previous = record.ell;
        let mut exponent = 0u32;
        while (&residual % record.ell).is_zero() {
            residual /= record.ell;
            exponent = exponent
                .checked_add(1)
                .ok_or_else(|| "prime exponent overflow".to_owned())?;
        }
        if exponent == 0 {
            vector.push(0);
            continue;
        }

        let positive_zero = (mod_u64(u, record.ell)
            + mul_mod(mod_u64(v, record.ell), record.positive_root, record.ell))
            % record.ell
            == 0;
        let negative_zero = (mod_u64(u, record.ell)
            + mul_mod(mod_u64(v, record.ell), record.negative_root, record.ell))
            % record.ell
            == 0;
        let (signed_exponent, orientation, root, selected_form) = match record.kind.as_str() {
            "ramified" if positive_zero && negative_zero => (
                i32::try_from(exponent).map_err(|_| "exponent does not fit i32")?,
                "ramified",
                record.positive_root,
                strings_form(&record.positive_form, &disc)?,
            ),
            "split" if positive_zero && !negative_zero => (
                i32::try_from(exponent).map_err(|_| "exponent does not fit i32")?,
                "positive",
                record.positive_root,
                strings_form(&record.positive_form, &disc)?,
            ),
            "split" if negative_zero && !positive_zero => (
                -i32::try_from(exponent).map_err(|_| "exponent does not fit i32")?,
                "negative",
                record.negative_root,
                strings_form(&record.negative_form, &disc)?,
            ),
            "split" if positive_zero && negative_zero => {
                return Err(format!(
                    "primitive candidate is divisible by both ideals above ell={}",
                    record.ell
                ));
            }
            _ => {
                return Err(format!(
                    "norm is divisible by ell={} but no oriented ideal divides alpha",
                    record.ell
                ));
            }
        };
        let charged_weight = record
            .edge_weight
            .checked_mul(exponent as u64)
            .ok_or_else(|| "per-prime weight overflow".to_owned())?;
        additive_weight = additive_weight
            .checked_add(charged_weight)
            .ok_or_else(|| "total weight overflow".to_owned())?;
        product = product.compose(&selected_form.pow(exponent as u64));
        vector.push(signed_exponent);
        factors.push(PrimeIdealFactor {
            ell: record.ell,
            exponent,
            orientation: orientation.to_owned(),
            root,
            edge_weight: record.edge_weight,
            charged_weight,
        });
    }

    let principal = BinaryQuadraticForm::principal(&disc)
        .ok_or_else(|| "principal form does not exist".to_owned())?;
    let complete = residual.is_one();
    let relation_verified = complete && product == principal;
    if complete && !relation_verified {
        return Err(
            "complete norm factorization did not compose to the principal class".to_owned(),
        );
    }
    let lambda = (u + v).mod_floor(subgroup_order);
    Ok(NormRelationCertificate {
        schema: SCHEMA_CERTIFICATE.to_owned(),
        u: u.to_string(),
        v: v.to_string(),
        primitive,
        alpha_norm: norm.to_string(),
        factors,
        exponent_vector: vector,
        residual_norm_factor: residual.to_string(),
        complete_factorization: complete,
        additive_chain_weight: additive_weight,
        reduced_known_ideal_product: form_strings(&product),
        principal_relation_verified: relation_verified,
        lambda_mod_subgroup_order: lambda.to_string(),
        map_status: "not_constructed".to_owned(),
        payoff_status: "not_evaluated".to_owned(),
    })
}

pub fn verify_norm_certificate(
    certificate: &NormRelationCertificate,
    subgroup_order: &BigInt,
    factor_base: &NormSieveFactorBase,
) -> Result<(), String> {
    if certificate.schema != SCHEMA_CERTIFICATE {
        return Err("unsupported norm-certificate schema".to_owned());
    }
    let u = parse_decimal(&certificate.u, "u")?;
    let v = parse_decimal(&certificate.v, "v")?;
    let replay = certify_norm_candidate(&u, &v, subgroup_order, factor_base)?;
    if replay != *certificate {
        return Err("norm-relation certificate differs from exact replay".to_owned());
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CenteredScanConfig {
    pub first_v: u64,
    pub last_v: u64,
    pub x_radius: u64,
    pub max_lattice_points: u64,
    /// Retain complete relations and partials whose residual has at most this
    /// many bits.  This is a retention boundary, not a weakness threshold.
    pub max_residual_bits: u64,
}

/// Exhaust a small, explicit box centred at the minimum of the norm form.
///
/// This is a deterministic correctness/reference enumerator.  A large run
/// should replace trial division with a segmented quadratic sieve while
/// preserving this function as the paired control.
pub fn scan_centered_box(
    config: &CenteredScanConfig,
    subgroup_order: &BigInt,
    factor_base: &NormSieveFactorBase,
) -> Result<Vec<NormRelationCertificate>, String> {
    if config.first_v == 0 || config.first_v > config.last_v {
        return Err("scan requires 1 <= first_v <= last_v".to_owned());
    }
    let width = config
        .x_radius
        .checked_mul(2)
        .and_then(|value| value.checked_add(1))
        .ok_or_else(|| "scan width overflow".to_owned())?;
    let rows = config.last_v - config.first_v + 1;
    let points = width
        .checked_mul(rows)
        .ok_or_else(|| "scan point count overflow".to_owned())?;
    if points > config.max_lattice_points {
        return Err(format!(
            "scan has {points} lattice points, above cap {}",
            config.max_lattice_points
        ));
    }
    let trace = parse_decimal(&factor_base.trace, "trace")?;
    let mut retained = Vec::new();
    for v_u64 in config.first_v..=config.last_v {
        let v = BigInt::from(v_u64);
        let center = (-&trace * &v).div_floor(&BigInt::from(2u8));
        let radius = i128::from(config.x_radius);
        for offset in -radius..=radius {
            let u = &center + BigInt::from(offset);
            if !u.gcd(&v).abs().is_one() {
                continue;
            }
            let certificate = certify_norm_candidate(&u, &v, subgroup_order, factor_base)?;
            let residual = BigUint::parse_bytes(certificate.residual_norm_factor.as_bytes(), 10)
                .ok_or_else(|| "certificate residual is not an unsigned integer".to_owned())?;
            if certificate.complete_factorization || residual.bits() <= config.max_residual_bits {
                retained.push(certificate);
            }
        }
    }
    retained.sort_by(|left, right| {
        let left_residual = BigUint::parse_bytes(left.residual_norm_factor.as_bytes(), 10).unwrap();
        let right_residual =
            BigUint::parse_bytes(right.residual_norm_factor.as_bytes(), 10).unwrap();
        left.additive_chain_weight
            .cmp(&right.additive_chain_weight)
            .then_with(|| left_residual.bits().cmp(&right_residual.bits()))
            .then_with(|| left.v.cmp(&right.v))
            .then_with(|| left.u.cmp(&right.u))
    });
    Ok(retained)
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: &str = "6277101735386680763835789423207666416083908700390324961279";
    const N: &str = "6277101735386680763835789423176059013767194773182842284081";
    const T: &str = "31607402316713927207482677199";
    const D: &str = "-24109379060336110122544161233113975664949272517896865359515";
    const C: &str = "14140398275856956083603613626459809774163796198179979683";

    fn integer(text: &str) -> BigInt {
        BigInt::parse_bytes(text.as_bytes(), 10).unwrap()
    }

    fn proxy_weights(max: u64) -> BTreeMap<u64, u64> {
        primes_through(max)
            .unwrap()
            .into_iter()
            .map(|ell| (ell, ell))
            .collect()
    }

    #[test]
    fn p192_factor_base_matches_frozen_split_and_ramified_primes() {
        let base = build_factor_base(
            &integer(D),
            &integer(T),
            &integer(P),
            113,
            "degree-linear-proxy/test-only",
            &proxy_weights(113),
        )
        .unwrap();
        let got: Vec<(u64, &str, u64, u64)> = base
            .primes
            .iter()
            .map(|entry| {
                (
                    entry.ell,
                    entry.kind.as_str(),
                    entry.positive_root,
                    entry.negative_root,
                )
            })
            .collect();
        assert_eq!(
            got,
            vec![
                (5, "ramified", 2, 2),
                (11, "ramified", 3, 3),
                (13, "split", 2, 5),
                (23, "split", 21, 22),
                (31, "ramified", 7, 7),
                (37, "split", 12, 35),
                (43, "split", 8, 26),
                (73, "split", 60, 67),
                (89, "split", 6, 83),
                (101, "split", 17, 70),
                (103, "split", 5, 36),
                (107, "split", 56, 68),
                (113, "split", 26, 42),
            ]
        );
    }

    #[test]
    fn symbolic_p192_sqrt_discriminant_is_a_partial_not_a_false_hit() {
        let base = build_factor_base(
            &integer(D),
            &integer(T),
            &integer(P),
            113,
            "degree-linear-proxy/test-only",
            &proxy_weights(113),
        )
        .unwrap();
        // sqrt(D)=2*pi-t.  Its missing prime C is deliberately outside the
        // u64 factor base and must remain a residual, not a claimed relation.
        let certificate =
            certify_norm_candidate(&-integer(T), &BigInt::from(2u8), &integer(N), &base).unwrap();
        assert_eq!(certificate.alpha_norm, integer(D).abs().to_string());
        assert_eq!(certificate.residual_norm_factor, C);
        assert!(!certificate.complete_factorization);
        assert!(!certificate.principal_relation_verified);
        assert_eq!(
            certificate
                .factors
                .iter()
                .map(|factor| factor.ell)
                .collect::<Vec<_>>(),
            vec![5, 11, 31]
        );
        assert_eq!(certificate.additive_chain_weight, 47);
        verify_norm_certificate(&certificate, &integer(N), &base).unwrap();
    }

    #[test]
    fn d_minus_23_positive_control_recovers_g2_cubed() {
        // pi^2-pi+6=0 has discriminant -23.  alpha=1+pi has norm 8 and
        // gives the standard non-scalar class relation g_2^3=1.
        let weights = BTreeMap::from([(2u64, 7u64)]);
        let base = build_factor_base(
            &BigInt::from(-23),
            &BigInt::from(1),
            &BigInt::from(6),
            2,
            "synthetic-exact-weight/test-only",
            &weights,
        )
        .unwrap();
        let certificate =
            certify_norm_candidate(&BigInt::one(), &BigInt::one(), &BigInt::from(101), &base)
                .unwrap();
        assert!(certificate.complete_factorization);
        assert!(certificate.principal_relation_verified);
        assert_eq!(certificate.alpha_norm, "8");
        assert_eq!(certificate.exponent_vector, vec![-3]);
        assert_eq!(certificate.additive_chain_weight, 21);
        assert_eq!(certificate.lambda_mod_subgroup_order, "2");
        verify_norm_certificate(&certificate, &BigInt::from(101), &base).unwrap();
    }

    #[test]
    fn centered_scan_is_capped_and_replays_retained_relations() {
        let weights = BTreeMap::from([(2u64, 1u64)]);
        let base = build_factor_base(
            &BigInt::from(-23),
            &BigInt::from(1),
            &BigInt::from(6),
            2,
            "synthetic-exact-weight/test-only",
            &weights,
        )
        .unwrap();
        let config = CenteredScanConfig {
            first_v: 1,
            last_v: 1,
            x_radius: 2,
            max_lattice_points: 5,
            max_residual_bits: 0,
        };
        let retained = scan_centered_box(&config, &BigInt::from(101), &base).unwrap();
        assert!(retained
            .iter()
            .any(|certificate| certificate.u == "1" && certificate.alpha_norm == "8"));
        for certificate in &retained {
            verify_norm_certificate(certificate, &BigInt::from(101), &base).unwrap();
        }
        let mut over_cap = config;
        over_cap.max_lattice_points = 4;
        assert!(scan_centered_box(&over_cap, &BigInt::from(101), &base).is_err());
    }

    #[test]
    fn certificate_mutations_and_nonprimitive_inputs_are_rejected() {
        let weights = BTreeMap::from([(2u64, 7u64)]);
        let base = build_factor_base(
            &BigInt::from(-23),
            &BigInt::from(1),
            &BigInt::from(6),
            2,
            "synthetic-exact-weight/test-only",
            &weights,
        )
        .unwrap();
        let mut certificate =
            certify_norm_candidate(&BigInt::one(), &BigInt::one(), &BigInt::from(101), &base)
                .unwrap();
        certificate.additive_chain_weight += 1;
        assert!(verify_norm_certificate(&certificate, &BigInt::from(101), &base).is_err());
        assert!(certify_norm_candidate(
            &BigInt::from(2),
            &BigInt::from(2),
            &BigInt::from(101),
            &base
        )
        .is_err());
    }
}
