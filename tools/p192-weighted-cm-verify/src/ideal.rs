use std::collections::BTreeMap;

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::{
    encoding::{
        encode_candidate_certificate, encode_factorization_certificate, sha256,
        CandidateCertificate, FactorizationCertificate, LargePrimeEntry, RationalFactor,
        SmallEntry,
    },
    factor_base::{kronecker_at_prime, roots_for_prime, VerifiedFactorBase},
    params::{self, LARGE_PRIME_BOUND, REF0_SHELL_ID},
    Result,
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RetainedCertificates {
    pub candidate: Vec<u8>,
    pub factorization: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CandidateEvaluation {
    pub status: u8,
    pub certificates: Option<RetainedCertificates>,
}

fn mul_mod(left: u64, right: u64, modulus: u64) -> u64 {
    ((u128::from(left) * u128::from(right)) % u128::from(modulus)) as u64
}

fn add_mod(left: u64, right: u64, modulus: u64) -> u64 {
    ((u128::from(left) + u128::from(right)) % u128::from(modulus)) as u64
}

fn pow_mod(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut accumulator = 1 % modulus;
    while exponent != 0 {
        if exponent & 1 == 1 {
            accumulator = mul_mod(accumulator, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exponent >>= 1;
    }
    accumulator
}

pub fn is_prime_u64(value: u64) -> bool {
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
    let trailing = (value - 1).trailing_zeros();
    let odd = (value - 1) >> trailing;
    // Deterministic for the full u64 domain.
    for base in [2u64, 325, 9_375, 28_178, 450_775, 9_780_504, 1_795_265_022] {
        let base = base % value;
        if base == 0 {
            continue;
        }
        let mut witness = pow_mod(base, odd, value);
        if witness == 1 || witness == value - 1 {
            continue;
        }
        let mut composite = true;
        for _ in 1..trailing {
            witness = mul_mod(witness, witness, value);
            if witness == value - 1 {
                composite = false;
                break;
            }
        }
        if composite {
            return false;
        }
    }
    true
}

fn gcd_u64(mut left: u64, mut right: u64) -> u64 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

fn pollard_rho(value: u64) -> Option<u64> {
    if value.is_multiple_of(2) {
        return Some(2);
    }
    if value.is_multiple_of(3) {
        return Some(3);
    }
    // Deterministic schedule.  The iteration ceiling makes inability explicit
    // rather than silently accepting a probable factorization.
    for constant in 1..=96u64 {
        for seed in [2u64, 3, 5, 7, 11] {
            let mut tortoise = seed % value;
            let mut hare = tortoise;
            for _ in 0..2_000_000u64 {
                tortoise = add_mod(mul_mod(tortoise, tortoise, value), constant, value);
                hare = add_mod(mul_mod(hare, hare, value), constant, value);
                hare = add_mod(mul_mod(hare, hare, value), constant, value);
                let difference = tortoise.abs_diff(hare);
                let divisor = gcd_u64(difference, value);
                if divisor == 1 {
                    continue;
                }
                if divisor != value {
                    return Some(divisor);
                }
                break;
            }
        }
    }
    None
}

fn factor_recursive(value: u64, factors: &mut Vec<u64>) -> Result<()> {
    if value == 1 {
        return Ok(());
    }
    if is_prime_u64(value) {
        factors.push(value);
        return Ok(());
    }
    let divisor = pollard_rho(value)
        .ok_or_else(|| format!("deterministic u64 factorization did not resolve {value}"))?;
    factor_recursive(divisor, factors)?;
    factor_recursive(value / divisor, factors)
}

pub fn factor_u64(value: u64) -> Result<Vec<u64>> {
    let mut factors = Vec::new();
    factor_recursive(value, &mut factors)?;
    factors.sort_unstable();
    let product = factors
        .iter()
        .try_fold(1u128, |product, factor| {
            product.checked_mul(u128::from(*factor))
        })
        .ok_or_else(|| "u64 factor product overflow".to_owned())?;
    if product != u128::from(value) {
        return Err("u64 factorization product mismatch".to_owned());
    }
    Ok(factors)
}

fn mod_u64(value: &BigInt, modulus: u64) -> u64 {
    value
        .mod_floor(&BigInt::from(modulus))
        .to_u64()
        .expect("residue fits u64")
}

fn orientation(u: &BigInt, v: &BigInt, prime: u64, roots: &[u64]) -> Result<usize> {
    let hits: Vec<_> = roots
        .iter()
        .enumerate()
        .filter_map(|(index, root)| {
            let residue = (u + v * BigInt::from(*root)).mod_floor(&BigInt::from(prime));
            residue.is_zero().then_some(index)
        })
        .collect();
    if hits.len() != 1 {
        return Err(format!(
            "primitive orientation at prime {prime} has {} roots",
            hits.len()
        ));
    }
    Ok(hits[0])
}

fn extended_gcd_inverse(value: &BigInt, modulus: &BigInt) -> Result<BigInt> {
    let result = value.extended_gcd(modulus);
    if result.gcd != BigInt::one() {
        return Err("coefficient is not invertible modulo the norm".to_owned());
    }
    Ok(result.x.mod_floor(modulus))
}

/// Reconstruct the HNF root of the principal ideal independently from the
/// certificate orientation list.  In basis (1, pi), (alpha) has HNF columns
/// (N, 0), (-r, 1), with r = -u/v modulo N for primitive non-scalars.
fn verify_principal_hnf_root(u: &BigInt, v: &BigInt, norm: &BigUint) -> Result<BigInt> {
    let modulus = BigInt::from_biguint(Sign::Plus, norm.clone());
    let inverse = extended_gcd_inverse(v, &modulus)?;
    let root = (-u * inverse).mod_floor(&modulus);
    if !(u + v * &root).mod_floor(&modulus).is_zero() {
        return Err("principal HNF root does not annihilate alpha".to_owned());
    }
    let characteristic = &root * &root - params::t() * &root + params::p();
    if !characteristic.mod_floor(&modulus).is_zero() {
        return Err("principal HNF root fails the Frobenius polynomial".to_owned());
    }
    Ok(root)
}

fn factor_map(factors: &[u64]) -> BTreeMap<u64, u64> {
    let mut result = BTreeMap::new();
    for factor in factors {
        *result.entry(*factor).or_default() += 1;
    }
    result
}

pub fn evaluate_candidate(
    candidate_index: u64,
    v_box: u64,
    x_box: u64,
    factor_base: &VerifiedFactorBase,
) -> Result<CandidateEvaluation> {
    evaluate_candidate_inner(candidate_index, v_box, x_box, factor_base, None)
}

pub fn evaluate_candidate_segmented(
    candidate_index: u64,
    v_box: u64,
    x_box: u64,
    factor_base: &VerifiedFactorBase,
    marked_entry_indices: &[usize],
) -> Result<CandidateEvaluation> {
    if marked_entry_indices
        .windows(2)
        .any(|pair| pair[0] >= pair[1])
        || marked_entry_indices
            .last()
            .is_some_and(|index| *index >= factor_base.manifest.entries.len())
    {
        return Err("segmented factor-base marks are unordered or out of range".to_owned());
    }
    evaluate_candidate_inner(
        candidate_index,
        v_box,
        x_box,
        factor_base,
        Some(marked_entry_indices),
    )
}

fn evaluate_candidate_inner(
    candidate_index: u64,
    v_box: u64,
    x_box: u64,
    factor_base: &VerifiedFactorBase,
    marked_entry_indices: Option<&[usize]>,
) -> Result<CandidateEvaluation> {
    let t = params::t();
    let numerator = BigInt::from(x_box) - &t * BigInt::from(v_box);
    if numerator.is_odd() {
        return Err("centered candidate produced nonintegral u".to_owned());
    }
    let u: BigInt = numerator / BigInt::from(2u8);
    let v = BigInt::from(v_box);
    if u.abs().gcd(&v) > BigInt::one() {
        return Ok(CandidateEvaluation {
            status: 0,
            certificates: None,
        });
    }

    let norm_signed = &u * &u + &t * &u * &v + params::p() * &v * &v;
    let norm = norm_signed
        .to_biguint()
        .ok_or_else(|| "candidate norm is not positive".to_owned())?;
    let centered = BigInt::from(x_box).pow(2u32) + (-params::d()) * &v * &v;
    if centered != BigInt::from(4u8) * &norm_signed {
        return Err("centered norm identity failed".to_owned());
    }

    let principal_root = verify_principal_hnf_root(&u, &v, &norm)?;
    let mut residual = norm.clone();
    let mut rational_factors = Vec::new();
    let mut small_entries = Vec::new();
    for (entry_index, entry) in factor_base.manifest.entries.iter().enumerate() {
        if marked_entry_indices.is_some_and(|indices| indices.binary_search(&entry_index).is_err())
        {
            continue;
        }
        let prime = BigUint::from(entry.ell);
        let mut exponent = 0u64;
        while (&residual % &prime).is_zero() {
            residual /= &prime;
            exponent += 1;
        }
        if exponent == 0 {
            if marked_entry_indices.is_some() {
                return Err(format!(
                    "segmented sieve falsely marked factor-base entry {}",
                    entry.ell
                ));
            }
            continue;
        }
        rational_factors.push(RationalFactor {
            prime: prime.clone(),
            exponent: BigUint::from(exponent),
        });
        let root_index = orientation(&u, &v, entry.ell, &entry.roots)?;
        if mod_u64(&principal_root, entry.ell) != entry.roots[root_index] {
            return Err(format!("HNF/orientation mismatch at prime {}", entry.ell));
        }
        let signed_exponent = if entry.kind == 1 {
            if exponent != 1 || root_index != 0 {
                return Err(format!("noncanonical ramified valuation at {}", entry.ell));
            }
            BigInt::one()
        } else if root_index == 0 {
            BigInt::from(exponent)
        } else {
            -BigInt::from(exponent)
        };
        small_entries.push(SmallEntry {
            ell: u32::try_from(entry.ell).map_err(|_| "base prime does not fit u32".to_owned())?,
            kind: entry.kind,
            exponent: signed_exponent,
        });
    }
    if let Some(indices) = marked_entry_indices {
        for (entry_index, entry) in factor_base.manifest.entries.iter().enumerate() {
            if indices.binary_search(&entry_index).is_err()
                && (&residual % BigUint::from(entry.ell)).is_zero()
            {
                return Err(format!(
                    "segmented sieve missed factor-base entry {}",
                    entry.ell
                ));
            }
        }
    }

    let limit_squared = BigUint::from(LARGE_PRIME_BOUND).pow(2u32);
    if residual > limit_squared {
        return Ok(CandidateEvaluation {
            status: 4,
            certificates: None,
        });
    }
    let residual_u64 = residual
        .to_u64()
        .ok_or_else(|| "bounded residual does not fit u64".to_owned())?;
    let residual_factors = factor_u64(residual_u64)?;
    let residual_map = factor_map(&residual_factors);

    if residual_factors.len() > 2
        || residual_factors
            .iter()
            .any(|prime| *prime > LARGE_PRIME_BOUND)
    {
        return Ok(CandidateEvaluation {
            status: 4,
            certificates: None,
        });
    }

    let mut lp_entries = Vec::new();
    for (prime, multiplicity) in residual_map {
        if prime <= params::ALGEBRAIC_BOUND {
            return Err(format!("missed algebraic-base divisor {prime}"));
        }
        if kronecker_at_prime(&params::d(), prime) != 1 {
            return Err(format!("non-split retained residual prime {prime}"));
        }
        let roots = roots_for_prime(prime, &params::t(), &params::p(), &params::d())?;
        if roots.len() != 2 {
            return Err(format!("retained split prime {prime} lacks two roots"));
        }
        let root_index = orientation(&u, &v, prime, &roots)?;
        if mod_u64(&principal_root, prime) != roots[root_index] {
            return Err(format!("LP HNF/orientation mismatch at prime {prime}"));
        }
        lp_entries.push(LargePrimeEntry {
            q: u32::try_from(prime).map_err(|_| "LP does not fit u32".to_owned())?,
            smaller_root: u32::try_from(roots[0])
                .map_err(|_| "root does not fit u32".to_owned())?,
            larger_root: u32::try_from(roots[1]).map_err(|_| "root does not fit u32".to_owned())?,
            sign: if root_index == 0 { 1 } else { 2 },
            multiplicity: BigUint::from(multiplicity),
        });
    }

    let candidate_type = match residual_factors.len() {
        0 => 0,
        1 => 1,
        2 => 2,
        _ => unreachable!("length was bounded above"),
    };
    let candidate = CandidateCertificate {
        candidate_type,
        curve_uid: params::CURVE_UID.to_owned(),
        factor_base_sha256: factor_base.sha256,
        shell_id: REF0_SHELL_ID.to_owned(),
        candidate_index,
        v_box,
        x_box,
        u,
        v,
        norm: norm.clone(),
        scalar_content: BigUint::one(),
        rational_factors: rational_factors.clone(),
        small_entries,
        lp_entries: lp_entries.clone(),
    };
    let candidate_bytes = encode_candidate_certificate(&candidate)?;
    let factorization = FactorizationCertificate {
        candidate_sha256: sha256(&candidate_bytes),
        norm,
        rational_factors,
        residual,
        residual_type: candidate_type,
        lp_entries,
    };
    let factorization_bytes = encode_factorization_certificate(&factorization)?;
    Ok(CandidateEvaluation {
        status: candidate_type + 1,
        certificates: Some(RetainedCertificates {
            candidate: candidate_bytes,
            factorization: factorization_bytes,
        }),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primality_and_factorization_controls() {
        assert!(is_prime_u64(2_147_483_647));
        assert!(!is_prime_u64(2_147_483_647u64.pow(2)));
        assert_eq!(factor_u64(65_537 * 65_539).unwrap(), vec![65_537, 65_539]);
        assert_eq!(
            factor_u64(2_147_483_647u64.pow(2)).unwrap(),
            vec![2_147_483_647; 2]
        );
    }
}
