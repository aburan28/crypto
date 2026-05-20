//! Shared helpers for Cryptopals Set 8.
//!
//! - `parse_big` — parse the big-decimal constants the challenges
//!   give us as Rust `BigUint`s.
//! - `crt_pair` — combine two congruences `x ≡ a₁ (mod m₁)` and
//!   `x ≡ a₂ (mod m₂)` into one with the Chinese Remainder Theorem.
//! - `crt_combine` — iterate `crt_pair` over a list.
//! - `small_factors` — trial-divide up to a bound and return the
//!   *unique* small prime factors (challenge 57 wants no repeats).
//! - `hmac_sha256` — the MAC the challenges authenticate with.

use num_bigint::{BigInt, BigUint, ToBigInt};
use num_integer::Integer;
use num_traits::{One, Zero};

/// Parse a base-10 string as a `BigUint`.  Panics on garbage input —
/// the only callers pass embedded string literals.
pub fn parse_big(s: &str) -> BigUint {
    s.parse().expect("set8_util::parse_big: invalid decimal literal")
}

/// Compute modular inverse via the extended Euclidean algorithm.
/// Same as `crate::utils::mod_inverse` but kept inline here so the
/// module is self-contained.
fn mod_inv(a: &BigUint, m: &BigUint) -> Option<BigUint> {
    let a_i = a.to_bigint().unwrap();
    let m_i = m.to_bigint().unwrap();
    let g = a_i.extended_gcd(&m_i);
    if g.gcd != BigInt::one() {
        return None;
    }
    let x = ((g.x % &m_i) + &m_i) % &m_i;
    x.to_biguint()
}

/// Combine two congruences with the Chinese Remainder Theorem.
///
/// Given `x ≡ a₁ (mod m₁)` and `x ≡ a₂ (mod m₂)` with
/// `gcd(m₁, m₂) = 1`, return `(a, m)` such that `x ≡ a (mod m)` and
/// `m = m₁ · m₂`.
pub fn crt_pair(
    (a1, m1): (&BigUint, &BigUint),
    (a2, m2): (&BigUint, &BigUint),
) -> (BigUint, BigUint) {
    // x = a1 + m1 · t where t ≡ (a2 - a1) · m1⁻¹ (mod m2).
    let m = m1 * m2;
    let m1_inv = mod_inv(m1, m2).expect("crt_pair: moduli not coprime");
    // (a2 - a1) mod m2, using a wide BigInt to avoid underflow.
    let diff = {
        let a1_i = a1.to_bigint().unwrap();
        let a2_i = a2.to_bigint().unwrap();
        let m2_i = m2.to_bigint().unwrap();
        let d = ((a2_i - a1_i) % &m2_i + &m2_i) % &m2_i;
        d.to_biguint().unwrap()
    };
    let t = (&diff * &m1_inv) % m2;
    ((a1 + m1 * t) % &m, m)
}

/// CRT over a vector of congruences `[(aᵢ, mᵢ)]`.  Iterates
/// `crt_pair` from left to right.  Panics on empty input.
pub fn crt_combine(parts: &[(BigUint, BigUint)]) -> (BigUint, BigUint) {
    assert!(!parts.is_empty(), "crt_combine: empty");
    let mut acc = parts[0].clone();
    for (a, m) in &parts[1..] {
        let (na, nm) = crt_pair((&acc.0, &acc.1), (a, m));
        acc = (na, nm);
    }
    acc
}

/// Trial-divide `n` and return its **distinct** small prime factors
/// up to `bound`.  Repeats are skipped.  Used by cryptopals 57 to
/// pick small subgroup orders.
pub fn small_factors(n: &BigUint, bound: u64) -> Vec<BigUint> {
    let mut out = Vec::new();
    let mut m = n.clone();
    let mut p: u64 = 2;
    while p < bound {
        let bp = BigUint::from(p);
        if (&m % &bp).is_zero() {
            out.push(bp.clone());
            while (&m % &bp).is_zero() {
                m /= &bp;
            }
        }
        p = if p == 2 { 3 } else { p + 2 };
    }
    out
}

/// HMAC-SHA256.  Cryptopals uses HMAC for the MAC that the attacker
/// brute-forces against the small-subgroup point.
pub fn hmac_sha256(key: &[u8], msg: &[u8]) -> [u8; 32] {
    use crate::hash::sha256::sha256;
    const BLOCK: usize = 64;
    let mut k0 = [0u8; BLOCK];
    if key.len() > BLOCK {
        let h = sha256(key);
        k0[..32].copy_from_slice(&h);
    } else {
        k0[..key.len()].copy_from_slice(key);
    }
    let mut ipad = [0u8; BLOCK];
    let mut opad = [0u8; BLOCK];
    for i in 0..BLOCK {
        ipad[i] = k0[i] ^ 0x36;
        opad[i] = k0[i] ^ 0x5c;
    }
    let mut inner = Vec::with_capacity(BLOCK + msg.len());
    inner.extend_from_slice(&ipad);
    inner.extend_from_slice(msg);
    let h1 = sha256(&inner);
    let mut outer = Vec::with_capacity(BLOCK + 32);
    outer.extend_from_slice(&opad);
    outer.extend_from_slice(&h1);
    sha256(&outer)
}

/// Convert a `BigUint` to a fixed-width big-endian byte vector.
/// Pads on the left with zeros to reach `width`.  Used as a
/// deterministic key-derivation input for the MAC.
pub fn biguint_to_bytes_be(x: &BigUint, width: usize) -> Vec<u8> {
    let raw = x.to_bytes_be();
    if raw.len() >= width {
        raw[raw.len() - width..].to_vec()
    } else {
        let mut out = vec![0u8; width - raw.len()];
        out.extend_from_slice(&raw);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crt_two_pairs() {
        let a1 = BigUint::from(2u32);
        let m1 = BigUint::from(3u32);
        let a2 = BigUint::from(3u32);
        let m2 = BigUint::from(5u32);
        let (a, m) = crt_pair((&a1, &m1), (&a2, &m2));
        assert_eq!(m, BigUint::from(15u32));
        assert_eq!(a, BigUint::from(8u32));
    }

    #[test]
    fn crt_chain() {
        let parts = vec![
            (BigUint::from(2u32), BigUint::from(3u32)),
            (BigUint::from(3u32), BigUint::from(5u32)),
            (BigUint::from(2u32), BigUint::from(7u32)),
        ];
        let (a, m) = crt_combine(&parts);
        assert_eq!(m, BigUint::from(105u32));
        assert_eq!(a, BigUint::from(23u32));
    }

    #[test]
    fn small_factors_basic() {
        let n = BigUint::from(2u32 * 3 * 7 * 11 * 13);
        let f = small_factors(&n, 100);
        let expected: Vec<BigUint> =
            [2u32, 3, 7, 11, 13].iter().map(|x| BigUint::from(*x)).collect();
        assert_eq!(f, expected);
    }

    #[test]
    fn hmac_test_vector() {
        // RFC 4231 test case 1.
        let key = vec![0x0bu8; 20];
        let msg = b"Hi There";
        let mac = hmac_sha256(&key, msg);
        let expected = hex::decode(
            "b0344c61d8db38535ca8afceaf0bf12b881dc200c9833da726e9376c2e32cff7",
        )
        .unwrap();
        assert_eq!(mac.as_slice(), expected.as_slice());
    }
}
