//! # Challenge 44 — DSA: recover key from repeated nonce
//!
//! Two signatures `(r₁, s₁, m₁)` and `(r₂, s₂, m₂)` with identical
//! `k` (and therefore identical `r`).  Solve:
//!
//! ```text
//!     k = (m₁ − m₂) / (s₁ − s₂)  mod q
//! ```
//!
//! then plug `k` back into the standard recovery formula.

use crate::cryptopals::challenge43::{dsa_params, dsa_x_from_known_k};
use crate::cryptopals::challenge28::sha1;
use crate::cryptopals::Report;
use num_bigint::{BigInt, BigUint, ToBigInt};
use num_integer::Integer;

fn mod_inv(a: &BigUint, n: &BigUint) -> BigUint {
    let a_i = a.to_bigint().unwrap();
    let n_i = n.to_bigint().unwrap();
    let g = a_i.extended_gcd(&n_i);
    let r = ((g.x % &n_i) + &n_i) % &n_i;
    r.to_biguint().unwrap()
}

pub fn k_from_repeated(
    m1: &[u8],
    m2: &[u8],
    s1: &BigUint,
    s2: &BigUint,
    q: &BigUint,
) -> BigUint {
    let h1 = BigUint::from_bytes_be(&sha1(m1)) % q;
    let h2 = BigUint::from_bytes_be(&sha1(m2)) % q;
    let num_i = h1.to_bigint().unwrap() - h2.to_bigint().unwrap();
    let den_i = s1.to_bigint().unwrap() - s2.to_bigint().unwrap();
    let q_i = q.to_bigint().unwrap();
    let num = ((num_i % &q_i) + &q_i) % &q_i;
    let den = ((den_i % &q_i) + &q_i) % &q_i;
    let num_b = num.to_biguint().unwrap();
    let den_b = den.to_biguint().unwrap();
    let den_inv = mod_inv(&den_b, q);
    (num_b * den_inv) % q
}

pub fn run() -> Report {
    let mut r = Report::new(44, "DSA: recover key from repeated nonce");
    let (_p, q, _g) = dsa_params();
    // Construct two synthetic signatures sharing k.
    let k = BigUint::from(0xABCDu32);
    let r_val = BigUint::from(0x1357u32); // pretend r = g^k mod p mod q
    let x = BigUint::from(0xDEADBEu32);
    let k_inv = mod_inv(&k, &q);
    let m1 = b"alpha";
    let m2 = b"beta";
    let h1 = BigUint::from_bytes_be(&sha1(m1)) % &q;
    let h2 = BigUint::from_bytes_be(&sha1(m2)) % &q;
    let s1 = ((&h1 + &x * &r_val) * &k_inv) % &q;
    let s2 = ((&h2 + &x * &r_val) * &k_inv) % &q;
    let k_rec = k_from_repeated(m1, m2, &s1, &s2, &q);
    assert_eq!(k_rec, k);
    let x_rec = dsa_x_from_known_k(m1, &k_rec, &r_val, &s1, &q);
    assert_eq!(x_rec, x);
    r.line("Recovered k and x from two sigs sharing a nonce.");
    let _ = BigInt::from(0);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn k_from_repeated_works() {
        assert!(super::run().success);
    }
}
