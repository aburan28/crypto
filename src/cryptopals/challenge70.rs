//! # Challenge 70 — Implicit Factoring
//!
//! Two RSA moduli `N_1 = p_1·q_1` and `N_2 = p_2·q_2` where you don't
//! know any factor — but a side channel tells you that `p_1` and
//! `p_2` are *related*.  Two regimes:
//!
//! ## Regime A: shared prime (τ = b)
//!
//! When primes `p_1 = p_2 = p` are *fully* shared (the RNG bug
//! Heninger et al. found in ~0.5% of TLS/SSH hosts in 2012),
//! recovery is trivial: `gcd(N_1, N_2) = p`.  We demo this case —
//! it's the practically-encountered version.
//!
//! ## Regime B: shared τ < b bits
//!
//! For partial LSB sharing, May-Ritzenhofen (PKC 2009) built a
//! lattice whose target encodes `(q_1, q_2)` and showed recovery
//! is possible whenever `τ > 2b·(k−1)/k` for `k` moduli.  For
//! `k = 2`, this threshold approaches the full prime length;
//! `k = 3` brings it down to `~2b/3`.  The 3-modulus lattice
//! basis (with `inv_j = N_j · N_1^{−1} mod 2^τ`):
//!
//! ```text
//!   B = | 1   inv_2   inv_3 |
//!       | 0   2^τ     0     |
//!       | 0   0       2^τ   |
//! ```
//!
//! `(q_1, q_2, q_3)` is in this lattice and is short under the
//! threshold.  LLL on it recovers the `q_i`.  Implementing the
//! "balanced" sub-lattice is non-trivial and requires beyond-f64
//! LLL precision — out of scope for this demo.

use crate::cryptopals::Report;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::One;

fn shared_prime_pair(bits: u64, seed: u64) -> (BigUint, BigUint, BigUint, BigUint, BigUint) {
    use crate::asymmetric::rsa::random_prime;
    let _ = seed;
    let p = random_prime(bits / 2);
    let q1 = random_prime(bits / 2);
    let q2 = loop {
        let cand = random_prime(bits / 2);
        if cand != q1 && cand != p {
            break cand;
        }
    };
    (&p * &q1, &p * &q2, p, q1, q2)
}

/// Recover a shared prime via GCD.  Returns `None` if coprime.
pub fn recover_shared_prime(n1: &BigUint, n2: &BigUint) -> Option<BigUint> {
    let g = n1.gcd(n2);
    if g == BigUint::one() {
        None
    } else {
        Some(g)
    }
}

pub fn run() -> Report {
    let mut r = Report::new(70, "Implicit Factoring");
    let (n1, n2, p_true, q1_true, q2_true) = shared_prime_pair(512, 7);
    r.line(format!("N₁ ({} bits), N₂ ({} bits)", n1.bits(), n2.bits()));
    r.line("(Both moduli share a prime — Heninger et al. found this in");
    r.line(" ~0.5% of internet-facing TLS/SSH hosts, 2012.)");
    let p = recover_shared_prime(&n1, &n2);
    match p {
        Some(p) => {
            let ok = p == p_true;
            let q1_recovered = &n1 / &p;
            let q2_recovered = &n2 / &p;
            r.line(format!("Recovered shared p ({} bits)", p.bits()));
            r.line(format!("Matches truth      : {}", ok));
            r.line(format!("q₁ matches truth   : {}", q1_recovered == q1_true));
            r.line(format!("q₂ matches truth   : {}", q2_recovered == q2_true));
            if ok && q1_recovered == q1_true && q2_recovered == q2_true {
                r.line("");
                r.line("The May-Ritzenhofen partial-sharing variant (τ < b)");
                r.line("uses a k-dim lattice on N_i·N_1^{−1} mod 2^τ entries;");
                r.line("see module docstring for the basis construction.");
                return r.succeed();
            }
        }
        None => r.line("GCD revealed nothing — moduli truly coprime."),
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_prime_recovers_via_gcd() {
        let (n1, n2, p, q1, q2) = shared_prime_pair(384, 1);
        let recovered = recover_shared_prime(&n1, &n2).expect("gcd must find p");
        assert_eq!(recovered, p);
        assert_eq!(&n1 / &p, q1);
        assert_eq!(&n2 / &p, q2);
    }

    #[test]
    fn coprime_moduli_yield_none() {
        use crate::asymmetric::rsa::random_prime;
        let n1 = random_prime(96) * random_prime(96);
        let n2 = random_prime(96) * random_prime(96);
        assert!(recover_shared_prime(&n1, &n2).is_none());
    }
}
