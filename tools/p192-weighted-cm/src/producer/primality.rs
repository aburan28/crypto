//! Deterministic primality and factorization for the REF-0 u64 residual path.

use std::collections::BTreeMap;

use super::factor_base::{mul_mod, pow_mod};
use super::Result;

fn gcd(mut left: u64, mut right: u64) -> u64 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

pub fn is_prime(value: u64) -> bool {
    if value < 2 {
        return false;
    }
    for prime in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if value.is_multiple_of(prime) {
            return value == prime;
        }
    }
    let mut odd = value - 1;
    let shifts = odd.trailing_zeros();
    odd >>= shifts;
    'witness: for base in [2u64, 325, 9_375, 28_178, 450_775, 9_780_504, 1_795_265_022] {
        if base % value == 0 {
            continue;
        }
        let mut probe = pow_mod(base % value, odd, value);
        if probe == 1 || probe == value - 1 {
            continue;
        }
        for _ in 1..shifts {
            probe = mul_mod(probe, probe, value);
            if probe == value - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

fn rho(value: u64) -> Result<u64> {
    if value.is_multiple_of(2) {
        return Ok(2);
    }
    if value.is_multiple_of(3) {
        return Ok(3);
    }
    // Deterministic restarts make the byte output independent of scheduling.
    for constant in 1u64..=256 {
        let mut slow = 2u64;
        let mut fast = 2u64;
        for _ in 0..2_000_000u64 {
            slow = (mul_mod(slow, slow, value) + constant) % value;
            fast = (mul_mod(fast, fast, value) + constant) % value;
            fast = (mul_mod(fast, fast, value) + constant) % value;
            let divisor = gcd(slow.abs_diff(fast), value);
            if divisor == 1 {
                continue;
            }
            if divisor != value {
                return Ok(divisor);
            }
            break;
        }
    }
    Err(format!(
        "deterministic u64 factorization did not converge for {value}"
    ))
}

fn collect(value: u64, output: &mut Vec<u64>) -> Result<()> {
    if value == 1 {
        return Ok(());
    }
    if is_prime(value) {
        output.push(value);
        return Ok(());
    }
    let divisor = rho(value)?;
    collect(divisor, output)?;
    collect(value / divisor, output)
}

pub fn factor(value: u64) -> Result<Vec<(u64, u32)>> {
    if value == 0 {
        return Err("cannot factor zero".to_owned());
    }
    let mut flat = Vec::new();
    collect(value, &mut flat)?;
    let mut counts = BTreeMap::new();
    for prime in flat {
        let exponent = counts.entry(prime).or_insert(0u32);
        *exponent = exponent
            .checked_add(1)
            .ok_or_else(|| "u64 factor exponent overflow".to_owned())?;
    }
    let factors = counts.into_iter().collect::<Vec<_>>();
    let product = factors
        .iter()
        .try_fold(1u128, |product, (prime, exponent)| {
            product.checked_mul((*prime as u128).pow(*exponent))
        });
    if product != Some(value as u128) || factors.iter().any(|(prime, _)| !is_prime(*prime)) {
        return Err("deterministic u64 factorization reconstruction failed".to_owned());
    }
    Ok(factors)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_u64_boundary_factorization() {
        let large = 2_147_483_647u64;
        assert!(is_prime(large));
        assert_eq!(factor(large).unwrap(), vec![(large, 1)]);
        assert_eq!(factor(large * large).unwrap(), vec![(large, 2)]);
        assert!(!is_prime(large * large));
    }

    #[test]
    fn known_pseudoprimes_are_rejected() {
        assert!(!is_prime(341_550_071_728_321));
        assert!(is_prime(18_446_744_073_709_551_557));
    }
}
