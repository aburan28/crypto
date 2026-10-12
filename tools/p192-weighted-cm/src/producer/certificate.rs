//! Canonical REF-0 candidate and factorization certificates.

use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, Signed, ToPrimitive, Zero};

use super::arithmetic::{centered_u, multiply_hnf, norm, pow_hnf, primitive, principal_hnf, Hnf};
use super::digest::sha256;
use super::factor_base::{characteristic_roots, mod_u64, mul_mod, Entry};
use super::primality;
use super::ref0::ReferenceRecord;
use super::schema::CURVE_UID;
use super::Result;

const CANDIDATE_DOMAIN: &[u8] = b"P192-WCM-CAND-v1\0";
const FACTORIZATION_DOMAIN: &[u8] = b"P192-WCM-FACT-v1\0";

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RationalFactor {
    pub prime: u64,
    pub exponent: u32,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmallEntry {
    pub ell: u64,
    pub kind: u8,
    pub exponent: BigInt,
    /// Selected characteristic root, used by exact HNF replay but omitted
    /// from the CAND small-entry encoding because the manifest fixes it.
    pub root: u64,
    pub positive_root: u64,
    pub negative_root: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LargePrimeEntry {
    pub prime: u64,
    pub smaller_root: u64,
    pub larger_root: u64,
    pub sign: u8,
    pub multiplicity: u32,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RetainedCertificate {
    pub candidate_bytes: Vec<u8>,
    pub candidate_sha256: [u8; 32],
    pub factorization_bytes: Vec<u8>,
    pub factorization_sha256: [u8; 32],
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Disposition {
    NonprimitiveDuplicate,
    Complete(RetainedCertificate),
    OneLargePrime(RetainedCertificate),
    TwoLargePrime(RetainedCertificate),
    Rejected,
    Invalid,
    Unresolved,
}

impl Disposition {
    pub fn status(&self) -> u8 {
        match self {
            Self::NonprimitiveDuplicate => 0,
            Self::Complete(_) => 1,
            Self::OneLargePrime(_) => 2,
            Self::TwoLargePrime(_) => 3,
            Self::Rejected => 4,
            Self::Invalid => 5,
            Self::Unresolved => 6,
        }
    }

    pub fn retained(&self) -> Option<&RetainedCertificate> {
        match self {
            Self::Complete(certificate)
            | Self::OneLargePrime(certificate)
            | Self::TwoLargePrime(certificate) => Some(certificate),
            _ => None,
        }
    }
}

fn append_u32(output: &mut Vec<u8>, value: usize, label: &str) -> Result<()> {
    let value = u32::try_from(value).map_err(|_| format!("{label} length exceeds u32"))?;
    output.extend_from_slice(&value.to_be_bytes());
    Ok(())
}

fn append_nat(output: &mut Vec<u8>, value: &BigUint) -> Result<()> {
    let magnitude = if value.is_zero() {
        Vec::new()
    } else {
        value.to_bytes_be()
    };
    append_u32(output, magnitude.len(), "nat")?;
    output.extend_from_slice(&magnitude);
    Ok(())
}

fn append_sint(output: &mut Vec<u8>, value: &BigInt) -> Result<()> {
    let (sign, magnitude) = value.to_bytes_be();
    output.push(match sign {
        Sign::NoSign => 0,
        Sign::Plus => 1,
        Sign::Minus => 2,
    });
    append_u32(output, magnitude.len(), "sint")?;
    output.extend_from_slice(&magnitude);
    Ok(())
}

pub(crate) fn encoded_sint(value: &BigInt) -> Result<Vec<u8>> {
    let mut output = Vec::new();
    append_sint(&mut output, value)?;
    Ok(output)
}

pub(crate) fn residual_prime_admissible(
    prime: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> bool {
    prime > algebraic_bound
        && characteristic_roots(prime, trace, field_norm, discriminant).len() == 2
}

fn append_rational_factors(output: &mut Vec<u8>, factors: &[RationalFactor]) -> Result<()> {
    append_u32(output, factors.len(), "rational factor count")?;
    for factor in factors {
        append_nat(output, &BigUint::from(factor.prime))?;
        append_nat(output, &BigUint::from(factor.exponent))?;
    }
    Ok(())
}

fn append_small_entries(output: &mut Vec<u8>, entries: &[SmallEntry]) -> Result<()> {
    append_u32(output, entries.len(), "small entry count")?;
    for entry in entries {
        let ell = u32::try_from(entry.ell)
            .map_err(|_| format!("small factor {} exceeds u32", entry.ell))?;
        output.extend_from_slice(&ell.to_be_bytes());
        output.push(entry.kind);
        append_sint(output, &entry.exponent)?;
    }
    Ok(())
}

fn append_lp_entries(output: &mut Vec<u8>, entries: &[LargePrimeEntry]) -> Result<()> {
    append_u32(output, entries.len(), "large-prime entry count")?;
    for entry in entries {
        output.extend_from_slice(
            &u32::try_from(entry.prime)
                .map_err(|_| "large prime exceeds u32".to_owned())?
                .to_be_bytes(),
        );
        output.extend_from_slice(
            &u32::try_from(entry.smaller_root)
                .map_err(|_| "large-prime root exceeds u32".to_owned())?
                .to_be_bytes(),
        );
        output.extend_from_slice(
            &u32::try_from(entry.larger_root)
                .map_err(|_| "large-prime root exceeds u32".to_owned())?
                .to_be_bytes(),
        );
        output.push(entry.sign);
        append_nat(output, &BigUint::from(entry.multiplicity))?;
    }
    Ok(())
}

fn orientation(u: &BigInt, v: &BigInt, ell: u64, roots: &[u64]) -> Result<(u8, BigInt)> {
    if roots.len() == 1 {
        let divides =
            (mod_u64(u, ell) + mul_mod(mod_u64(v, ell), roots[0], ell)).is_multiple_of(ell);
        if !divides {
            return Err(format!("ramified ideal orientation mismatch at {ell}"));
        }
        return Ok((1, BigInt::one()));
    }
    if roots.len() != 2 {
        return Err(format!("split orientation needs two roots at {ell}"));
    }
    let first = (mod_u64(u, ell) + mul_mod(mod_u64(v, ell), roots[0], ell)).is_multiple_of(ell);
    let second = (mod_u64(u, ell) + mul_mod(mod_u64(v, ell), roots[1], ell)).is_multiple_of(ell);
    match (first, second) {
        (true, false) => Ok((1, BigInt::one())),
        (false, true) => Ok((2, -BigInt::one())),
        _ => Err(format!("non-unique split orientation at {ell}")),
    }
}

fn build_retained(
    record: &ReferenceRecord,
    u: &BigInt,
    v: &BigInt,
    alpha_norm: &BigUint,
    residual: u64,
    type_code: u8,
    factor_base_sha256: &[u8; 32],
    rational_factors: &[RationalFactor],
    small_entries: &[SmallEntry],
    lp_entries: &[LargePrimeEntry],
    field_norm: &BigInt,
    trace: &BigInt,
) -> Result<RetainedCertificate> {
    if type_code > 2 {
        return Err("candidate type is outside 0..=2".to_owned());
    }
    let mut previous = 0u64;
    let mut rational_product = BigUint::one();
    for factor in rational_factors {
        if factor.prime < 2 || factor.prime <= previous || factor.exponent == 0 {
            return Err("noncanonical rational-factor list".to_owned());
        }
        previous = factor.prime;
        rational_product *= BigUint::from(factor.prime).pow(factor.exponent);
    }
    if rational_factors.len() != small_entries.len() {
        return Err("rational/small factor counts differ".to_owned());
    }
    previous = 0;
    for (factor, entry) in rational_factors.iter().zip(small_entries) {
        if entry.ell <= previous
            || entry.ell != factor.prime
            || !matches!(entry.kind, 0 | 1)
            || entry.exponent.is_zero()
            || entry.exponent.abs() != BigInt::from(factor.exponent)
            || (entry.kind == 1 && entry.exponent != BigInt::one())
            || entry.root >= entry.ell
            || entry.positive_root > entry.negative_root
            || entry.negative_root >= entry.ell
            || (entry.kind == 0
                && ((entry.exponent.is_positive() && entry.root != entry.positive_root)
                    || (entry.exponent.is_negative() && entry.root != entry.negative_root)
                    || entry.positive_root == entry.negative_root))
            || (entry.kind == 1
                && (entry.root != entry.positive_root
                    || entry.positive_root != entry.negative_root))
        {
            return Err("noncanonical small-entry list".to_owned());
        }
        previous = entry.ell;
    }
    previous = 0;
    let mut lp_product = BigUint::one();
    let mut lp_multiplicity = 0u32;
    for entry in lp_entries {
        if entry.prime <= previous
            || entry.multiplicity == 0
            || entry.smaller_root >= entry.larger_root
            || entry.larger_root >= entry.prime
            || !matches!(entry.sign, 1 | 2)
        {
            return Err("noncanonical large-prime entry list".to_owned());
        }
        previous = entry.prime;
        lp_multiplicity = lp_multiplicity
            .checked_add(entry.multiplicity)
            .ok_or_else(|| "large-prime multiplicity overflow".to_owned())?;
        lp_product *= BigUint::from(entry.prime).pow(entry.multiplicity);
    }
    let residual_big = BigUint::from(residual);
    let type_consistent = match type_code {
        0 => residual == 1 && lp_entries.is_empty() && lp_multiplicity == 0,
        1 => residual > 1 && lp_multiplicity == 1,
        2 => residual > 1 && lp_multiplicity == 2,
        _ => false,
    };
    if !type_consistent
        || lp_product != residual_big
        || rational_product * &residual_big != *alpha_norm
    {
        return Err("candidate type/factor/product invariant failed".to_owned());
    }

    let mut ideal = Hnf::identity();
    for entry in small_entries {
        let magnitude = entry
            .exponent
            .magnitude()
            .to_u32()
            .ok_or_else(|| "small ideal exponent does not fit u32".to_owned())?;
        let power = pow_hnf(
            &Hnf::prime(entry.ell, entry.root)?,
            magnitude,
            field_norm,
            trace,
        )?;
        ideal = multiply_hnf(&ideal, &power, field_norm, trace)?;
    }
    for entry in lp_entries {
        let root = if entry.sign == 1 {
            entry.smaller_root
        } else if entry.sign == 2 {
            entry.larger_root
        } else {
            return Err("invalid large-prime orientation sign".to_owned());
        };
        let power = pow_hnf(
            &Hnf::prime(entry.prime, root)?,
            entry.multiplicity,
            field_norm,
            trace,
        )?;
        ideal = multiply_hnf(&ideal, &power, field_norm, trace)?;
    }
    let principal = principal_hnf(u, v, field_norm, trace)?;
    if ideal != principal
        || ideal.determinant() != BigInt::from_biguint(Sign::Plus, alpha_norm.clone())
    {
        return Err("oriented prime-power HNF does not equal the principal ideal".to_owned());
    }

    let mut candidate = Vec::new();
    candidate.extend_from_slice(CANDIDATE_DOMAIN);
    candidate.push(type_code);
    append_u32(&mut candidate, CURVE_UID.len(), "curve UID")?;
    candidate.extend_from_slice(CURVE_UID.as_bytes());
    candidate.extend_from_slice(factor_base_sha256);
    append_u32(
        &mut candidate,
        super::ref0::REFERENCE_SHELL.len(),
        "shell ID",
    )?;
    candidate.extend_from_slice(super::ref0::REFERENCE_SHELL.as_bytes());
    candidate.extend_from_slice(&record.index.to_be_bytes());
    candidate.extend_from_slice(&record.v.to_be_bytes());
    candidate.extend_from_slice(&record.x.to_be_bytes());
    append_sint(&mut candidate, u)?;
    append_sint(&mut candidate, v)?;
    append_nat(&mut candidate, alpha_norm)?;
    append_nat(&mut candidate, &BigUint::one())?;
    append_rational_factors(&mut candidate, rational_factors)?;
    append_small_entries(&mut candidate, small_entries)?;
    append_lp_entries(&mut candidate, lp_entries)?;
    let candidate_sha256 = sha256(&candidate);

    let mut factorization = Vec::new();
    factorization.extend_from_slice(FACTORIZATION_DOMAIN);
    factorization.extend_from_slice(&candidate_sha256);
    append_nat(&mut factorization, alpha_norm)?;
    append_rational_factors(&mut factorization, rational_factors)?;
    append_nat(&mut factorization, &BigUint::from(residual))?;
    factorization.push(type_code);
    append_lp_entries(&mut factorization, lp_entries)?;
    let factorization_sha256 = sha256(&factorization);
    Ok(RetainedCertificate {
        candidate_bytes: candidate,
        candidate_sha256,
        factorization_bytes: factorization,
        factorization_sha256,
    })
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PreparedCandidate {
    pub(crate) u: BigInt,
    pub(crate) v: BigInt,
    pub(crate) alpha_norm: BigUint,
    pub(crate) residual: BigUint,
    pub(crate) rational_factors: Vec<RationalFactor>,
    pub(crate) small_entries: Vec<SmallEntry>,
    pub(crate) invalid: bool,
}

pub(crate) fn prepare_candidate(
    record: &ReferenceRecord,
    trace: &BigInt,
    field_norm: &BigInt,
) -> Result<Option<PreparedCandidate>> {
    let v = BigInt::from(record.v);
    let u = centered_u(record.x, record.v, trace)?;
    if !primitive(&u, &v) {
        return Ok(None);
    }
    let alpha_norm_signed = norm(&u, &v, trace, field_norm);
    let alpha_norm = alpha_norm_signed
        .to_biguint()
        .ok_or_else(|| "REF-0 norm must be positive".to_owned())?;
    Ok(Some(PreparedCandidate {
        u,
        v,
        residual: alpha_norm.clone(),
        alpha_norm,
        rational_factors: Vec::new(),
        small_entries: Vec::new(),
        invalid: false,
    }))
}

pub(crate) fn divide_prepared_by_entry(
    prepared: &mut PreparedCandidate,
    entry: &Entry,
) -> Result<()> {
    let mut exponent = 0u32;
    while (&prepared.residual % entry.ell).is_zero() {
        prepared.residual /= entry.ell;
        exponent = exponent
            .checked_add(1)
            .ok_or_else(|| "small-prime exponent overflow".to_owned())?;
    }
    if exponent == 0 {
        return Ok(());
    }
    let (sign, orientation_unit) =
        match orientation(&prepared.u, &prepared.v, entry.ell, &entry.roots) {
            Ok(value) => value,
            Err(_) => {
                prepared.invalid = true;
                return Ok(());
            }
        };
    if entry.kind == 1 && exponent != 1 {
        prepared.invalid = true;
        return Ok(());
    }
    let signed_exponent = if entry.kind == 1 {
        BigInt::one()
    } else {
        orientation_unit * BigInt::from(exponent)
    };
    prepared.rational_factors.push(RationalFactor {
        prime: entry.ell,
        exponent,
    });
    prepared.small_entries.push(SmallEntry {
        ell: entry.ell,
        kind: entry.kind,
        exponent: signed_exponent,
        root: if entry.kind == 1 || sign == 1 {
            entry.roots[0]
        } else {
            entry.roots[1]
        },
        positive_root: entry.roots[0],
        negative_root: *entry.roots.last().expect("factor-base entry has roots"),
    });
    Ok(())
}

pub(crate) fn finish_prepared(
    record: &ReferenceRecord,
    prepared: PreparedCandidate,
    factor_base_sha256: &[u8; 32],
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
    algebraic_bound: u64,
    large_prime_bound: u64,
) -> Result<Disposition> {
    if prepared.invalid {
        return Ok(Disposition::Invalid);
    }
    if prepared.residual.is_one() {
        return Ok(
            match build_retained(
                record,
                &prepared.u,
                &prepared.v,
                &prepared.alpha_norm,
                1,
                0,
                factor_base_sha256,
                &prepared.rational_factors,
                &prepared.small_entries,
                &[],
                field_norm,
                trace,
            ) {
                Ok(certificate) => Disposition::Complete(certificate),
                Err(_) => Disposition::Invalid,
            },
        );
    }
    let bound_squared = BigUint::from(large_prime_bound) * large_prime_bound;
    if prepared.residual > bound_squared {
        return Ok(Disposition::Rejected);
    }
    let residual_u64 = prepared
        .residual
        .to_u64()
        .ok_or_else(|| "bounded REF-0 residual does not fit u64".to_owned())?;
    let factors = match primality::factor(residual_u64) {
        Ok(factors) => factors,
        Err(_) => return Ok(Disposition::Unresolved),
    };
    let multiplicity_count = factors
        .iter()
        .try_fold(0u32, |sum, (_, exponent)| sum.checked_add(*exponent))
        .ok_or_else(|| "large-prime multiplicity overflow".to_owned())?;
    if multiplicity_count > 2 || factors.iter().any(|(prime, _)| *prime > large_prime_bound) {
        return Ok(Disposition::Rejected);
    }
    let mut lp_entries = Vec::new();
    for (prime, multiplicity) in &factors {
        // A remaining algebraic-base divisor or an inert norm factor is an
        // invariant violation, not a normal policy rejection.
        let roots = characteristic_roots(*prime, trace, field_norm, discriminant);
        if !residual_prime_admissible(*prime, algebraic_bound, trace, field_norm, discriminant) {
            return Ok(Disposition::Invalid);
        }
        let (sign, _) = match orientation(&prepared.u, &prepared.v, *prime, &roots) {
            Ok(value) => value,
            Err(_) => return Ok(Disposition::Invalid),
        };
        lp_entries.push(LargePrimeEntry {
            prime: *prime,
            smaller_root: roots[0],
            larger_root: roots[1],
            sign,
            multiplicity: *multiplicity,
        });
    }
    let type_code = match multiplicity_count {
        1 => 1,
        2 => 2,
        _ => return Ok(Disposition::Invalid),
    };
    let retained = match build_retained(
        record,
        &prepared.u,
        &prepared.v,
        &prepared.alpha_norm,
        residual_u64,
        type_code,
        factor_base_sha256,
        &prepared.rational_factors,
        &prepared.small_entries,
        &lp_entries,
        field_norm,
        trace,
    ) {
        Ok(certificate) => certificate,
        Err(_) => return Ok(Disposition::Invalid),
    };
    Ok(if type_code == 1 {
        Disposition::OneLargePrime(retained)
    } else {
        Disposition::TwoLargePrime(retained)
    })
}

pub fn classify(
    record: &ReferenceRecord,
    factor_base: &[Entry],
    factor_base_sha256: &[u8; 32],
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
    algebraic_bound: u64,
    large_prime_bound: u64,
) -> Result<Disposition> {
    let Some(mut prepared) = prepare_candidate(record, trace, field_norm)? else {
        return Ok(Disposition::NonprimitiveDuplicate);
    };
    for entry in factor_base {
        divide_prepared_by_entry(&mut prepared, entry)?;
    }
    finish_prepared(
        record,
        prepared,
        factor_base_sha256,
        trace,
        field_norm,
        discriminant,
        algebraic_bound,
        large_prime_bound,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn d23_certificate(small: SmallEntry, u: BigInt, norm: u64) -> Result<RetainedCertificate> {
        build_retained(
            &ReferenceRecord {
                index: 0,
                v: 1,
                x: 3,
            },
            &u,
            &BigInt::one(),
            &BigUint::from(norm),
            1,
            0,
            &[0u8; 32],
            &[RationalFactor {
                prime: 2,
                exponent: 3,
            }],
            &[small],
            &[],
            &BigInt::from(6u8),
            &BigInt::one(),
        )
    }

    fn d23_small() -> SmallEntry {
        SmallEntry {
            ell: 2,
            kind: 0,
            exponent: BigInt::from(-3),
            root: 1,
            positive_root: 0,
            negative_root: 1,
        }
    }

    #[test]
    fn nat_and_sint_are_minimal() {
        let mut bytes = Vec::new();
        append_nat(&mut bytes, &BigUint::zero()).unwrap();
        assert_eq!(bytes, [0, 0, 0, 0]);
        bytes.clear();
        append_sint(&mut bytes, &BigInt::from(-256)).unwrap();
        assert_eq!(bytes, [2, 0, 0, 0, 2, 1, 0]);
    }

    #[test]
    fn split_orientation_uses_the_smaller_root_as_positive() {
        let u = BigInt::from(-2);
        let v = BigInt::one();
        assert_eq!(orientation(&u, &v, 13, &[2, 5]).unwrap().0, 1);
        let u = BigInt::from(-5);
        assert_eq!(orientation(&u, &v, 13, &[2, 5]).unwrap().0, 2);
    }

    #[test]
    fn d23_cand_and_fact_bytes_are_golden() {
        let certificate = d23_certificate(d23_small(), BigInt::one(), 8).unwrap();
        assert_eq!(
            super::super::ref0::hex(&certificate.candidate_bytes),
            "503139322d57434d2d43414e442d763100000000005775726e3a65632d7265636f72643a313a7368613235363a353533316334613038626462363462366538366136653330653961303861613537656465663761663135616335663664346432613833613533626632663634360000000000000000000000000000000000000000000000000000000000000000000000055245462d3000000000000000000000000000000001000000000000000301000000010101000000010100000001080000000101000000010000000102000000010300000001000000020002000000010300000000"
        );
        assert_eq!(
            super::super::digest::sha256_hex(&certificate.candidate_bytes),
            "e5768dc26538c12671abbe21e492b5b4df780defd5403f777199c33bfea6429d"
        );
        assert_eq!(
            super::super::ref0::hex(&certificate.factorization_bytes),
            "503139322d57434d2d464143542d763100e5768dc26538c12671abbe21e492b5b4df780defd5403f777199c33bfea6429d0000000108000000010000000102000000010300000001010000000000"
        );
        assert_eq!(
            super::super::digest::sha256_hex(&certificate.factorization_bytes),
            "495b7c078e288b9fa04428272b1f7f324d5093a02b3c72ff1632c8b8a0ef0c1d"
        );
        assert_eq!(
            &certificate.factorization_bytes
                [b"P192-WCM-FACT-v1\0".len()..b"P192-WCM-FACT-v1\0".len() + 32],
            &certificate.candidate_sha256
        );
    }

    #[test]
    fn d23_structural_mutations_fail_closed() {
        let mut root = d23_small();
        root.root = 0;
        assert!(d23_certificate(root, BigInt::one(), 8).is_err());

        let mut sign = d23_small();
        sign.exponent = BigInt::from(3);
        assert!(d23_certificate(sign, BigInt::one(), 8).is_err());

        let mut exponent = d23_small();
        exponent.exponent = BigInt::from(-2);
        assert!(d23_certificate(exponent, BigInt::one(), 8).is_err());

        assert!(d23_certificate(d23_small(), BigInt::from(2u8), 8).is_err());
        assert!(d23_certificate(d23_small(), BigInt::one(), 4).is_err());
    }
}
