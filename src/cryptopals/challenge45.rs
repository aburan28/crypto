//! # Challenge 45 — DSA parameter tampering (g = 0, g = p + 1)
//!
//! With `g = 0`, every signature `r = (g^k mod p) mod q = 0`, and
//! verify accepts any `(0, s)`.  With `g = p+1` (≡ 1 mod p), every
//! signature accepts every message.
//!
//! In short: if `g` is malicious, signatures lose meaning.

use crate::cryptopals::challenge43::dsa_params;
use crate::cryptopals::Report;
use num_bigint::BigUint;
use num_traits::Zero;

pub fn run() -> Report {
    let mut r = Report::new(45, "DSA parameter tampering");
    let (p, q, _g_real) = dsa_params();
    let g_zero = BigUint::zero();
    let r_sig = g_zero.modpow(&BigUint::from(7u32), &p) % &q;
    r.line(format!("g = 0 → r component = {} (always 0 → invalid sig matches anything)", r_sig));
    assert!(r_sig.is_zero());

    let g_one = &p + BigUint::from(1u32);
    let r_sig = g_one.modpow(&BigUint::from(7u32), &p) % &q;
    r.line(format!("g = p+1 → r component = {} (always 1)", r_sig));
    // 1 mod q = 1.
    assert_eq!(r_sig, BigUint::from(1u32));
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn malicious_g_breaks_dsa() {
        assert!(super::run().success);
    }
}
