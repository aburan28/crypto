//! # Challenge 35 — DH with malicious group parameters
//!
//! In a negotiated-group DH, Mallory replaces the agreed `g` with
//! a degenerate value:
//!
//! - `g = 1` → both publics = 1, both shared secrets = 1.
//! - `g = p` → both publics = 0, both shared secrets = 0.
//! - `g = p - 1` → publics ∈ {1, p-1}, shared secrets ∈ {1, p-1}.
//!
//! In every case, Mallory knows the shared secret a priori.

use crate::cryptopals::challenge33::nist_p_g;
use crate::cryptopals::Report;
use num_bigint::{BigUint, RandBigInt};
use num_traits::{One, Zero};
use rand::SeedableRng;
use rand::rngs::StdRng;

pub fn run() -> Report {
    let mut r = Report::new(35, "DH with malicious g");
    let (p, _g) = nist_p_g();
    let mut rng = StdRng::seed_from_u64(35);
    let a = rng.gen_biguint_range(&BigUint::from(2u32), &p);
    let b = rng.gen_biguint_range(&BigUint::from(2u32), &p);

    // Case 1: g = 1
    let g = BigUint::one();
    let s_a = g.modpow(&a, &p).modpow(&b, &p);
    r.line(format!("g=1 :  shared = {} (always 1)", s_a));
    assert_eq!(s_a, BigUint::one());

    // Case 2: g = p
    let g = p.clone();
    let s_a = g.modpow(&a, &p).modpow(&b, &p);
    r.line(format!("g=p :  shared = {} (always 0)", s_a));
    assert!(s_a.is_zero());

    // Case 3: g = p-1.  publics = ±1, shared = ±1.
    let g = &p - BigUint::one();
    let s_a = g.modpow(&a, &p).modpow(&b, &p);
    r.line(format!("g=p-1: shared in {{1, p-1}} → {}", s_a));
    assert!(s_a == BigUint::one() || s_a == &p - BigUint::one());

    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn malicious_g_works() {
        assert!(super::run().success);
    }
}
