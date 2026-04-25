//! Shared utilities: randomness, encoding, and modular arithmetic helpers.

pub mod encoding;
pub mod random;

pub use encoding::*;
pub use random::*;

use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, Zero};

/// Modular exponentiation: `base^exp mod modulus` using square-and-multiply.
/// Works for any size integers; does NOT require `modulus` to be prime.
pub fn mod_pow(base: &BigUint, exp: &BigUint, modulus: &BigUint) -> BigUint {
    let mut result = BigUint::one();
    let mut b = base % modulus;
    let mut e = exp.clone();
    let zero = BigUint::zero();
    let one = BigUint::one();

    while e > zero {
        if &e & &one == one {
            result = (&result * &b) % modulus;
        }
        b = (&b * &b) % modulus;
        e >>= 1;
    }
    result
}

/// Modular inverse via the extended Euclidean algorithm.
/// Returns `Some(x)` such that `a * x ≡ 1 (mod m)`, or `None` if no inverse exists.
/// Works whether or not `m` is prime (unlike Fermat's little theorem).
pub fn mod_inverse(a: &BigUint, m: &BigUint) -> Option<BigUint> {
    if a.is_zero() {
        return None;
    }
    let a_i = BigInt::from_biguint(Sign::Plus, a.clone());
    let m_i = BigInt::from_biguint(Sign::Plus, m.clone());

    let (mut old_r, mut r) = (a_i.clone(), m_i.clone());
    let (mut old_s, mut s) = (BigInt::one(), BigInt::zero());

    while !r.is_zero() {
        let q = &old_r / &r;
        let tmp = old_r - &q * &r;
        old_r = r;
        r = tmp;
        let tmp = old_s - q * &s;
        old_s = s;
        s = tmp;
    }

    // gcd must be 1 for the inverse to exist
    if old_r != BigInt::one() && old_r != BigInt::from(-1i32) {
        return None;
    }

    let result = ((old_s % &m_i) + &m_i) % &m_i;
    result.to_biguint()
}
