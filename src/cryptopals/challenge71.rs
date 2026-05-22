//! # Challenge 71 — Heninger-Shacham: Cold-Boot RSA
//!
//! After a power-cut, DRAM bits decay asymmetrically toward 0 or 1.
//! An attacker who reads the memory before complete decay sees a
//! corrupted version of the RSA private key: each bit position is
//! either "known correct" or "unknown".
//!
//! Heninger-Shacham (CRYPTO 2009) showed that with ~27% of bits
//! known across the six RSA-CRT values
//! `{p, q, d, d_p, d_q, q_inv}`, a branch-and-prune search recovers
//! the full key in polynomial time.  The key relations are:
//!
//! ```text
//!     p · q          =  N
//!     e · d          ≡  1   (mod (p-1)(q-1))
//!     e · d_p        ≡  1   (mod p − 1)
//!     e · d_q        ≡  1   (mod q − 1)
//! ```
//!
//! At each iteration, we extend the candidate values by one more
//! low-order bit slice (typically 1 bit at a time) and check the
//! relations *modulo 2^(i+1)*.  Branches that violate any relation
//! die.  In practice the surviving tree stays narrow (~16-128
//! branches per level).
//!
//! ## This implementation
//!
//! We do the slightly easier variant: recover `(p, q)` given partial
//! bits of *just* `p`, using only the constraint `p · q = N`.  At
//! level `i`, we have a candidate `p mod 2^i`; for each "unknown
//! slice" we enumerate 2-bit choices, derive `q mod 2^i` from
//! `q ≡ N · p⁻¹ (mod 2^i)`, and prune branches whose `q` extension
//! conflicts with anything we know about `q`.
//!
//! Smaller scope than full HS but the same algorithmic shape; the
//! full 6-variable variant is a mechanical extension.

use crate::cryptopals::Report;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

/// A "partial bit" array: each position is `Some(0)`, `Some(1)`, or
/// `None` (unknown).  Index 0 is the least-significant bit.
pub type PartialBits = Vec<Option<u8>>;

/// Take the true bit-string of a `BigUint` and corrupt it by setting
/// each position to `None` with probability `unknown_prob`.
pub fn corrupt(b: &BigUint, len: usize, unknown_prob: f64, rng: &mut StdRng) -> PartialBits {
    (0..len)
        .map(|i| {
            if rng.gen::<f64>() < unknown_prob {
                None
            } else {
                Some((((b >> i) & BigUint::one()) == BigUint::one()) as u8)
            }
        })
        .collect()
}

/// Branch-and-prune: given partial bit-strings of both `p` and `q`
/// plus the modulus `N`, recover the full `p`.  At each bit `i ≥ 1`
/// we try each unknown-bit value, derive the corresponding `q`-bit
/// from `p · q ≡ N (mod 2^(i+1))`, and prune branches whose derived
/// `q`-bit conflicts with `partial_q`.  When both bits of `(p, q)`
/// at level `i` are known, we just verify consistency.
pub fn cold_boot_recover_p(
    partial_p: &PartialBits,
    partial_q: &PartialBits,
    n: &BigUint,
    bits: usize,
) -> Option<BigUint> {
    if bits < 2 {
        return None;
    }
    // Initial state at bit 0: p_lsb = 1 (primes are odd), so
    // q_lsb = N (mod 2) = 1.
    let mut candidates: Vec<(BigUint, BigUint)> =
        vec![(BigUint::one(), BigUint::one())];

    for level in 1..bits {
        let mask = BigUint::one() << (level + 1);
        let mut next: Vec<(BigUint, BigUint)> = Vec::new();
        let p_bit_known = partial_p.get(level).and_then(|b| *b);
        let q_bit_known = partial_q.get(level).and_then(|b| *b);
        for (p_so_far, q_so_far) in &candidates {
            let p_candidates: Vec<u8> = match p_bit_known {
                Some(b) => vec![b],
                None => vec![0, 1],
            };
            for p_bit in p_candidates {
                let new_p = p_so_far + (BigUint::from(p_bit) << level);
                // Derive q_bit: (new_p · (q_so_far + q_bit·2^level)) ≡ N (mod 2^(level+1)).
                let pq = (&new_p * q_so_far) % &mask;
                let n_mod = n % &mask;
                let needed: BigUint = if n_mod >= pq {
                    n_mod - pq
                } else {
                    n_mod + &mask - pq
                };
                // The "needed" amount is q_bit · new_p · 2^level mod 2^(level+1).
                // The new_p's parity at position 0 is 1 (odd), so
                // new_p · 2^level mod 2^(level+1) = 2^level (since other bits don't survive).
                // Thus q_bit = (needed >> level) & 1.
                let derived_q_bit = ((&needed >> level) & BigUint::one()).to_u32_digits()
                    .get(0).copied().unwrap_or(0) as u8;
                if let Some(known) = q_bit_known {
                    if derived_q_bit != known {
                        continue; // prune
                    }
                }
                let new_q = q_so_far + (BigUint::from(derived_q_bit) << level);
                next.push((new_p, new_q));
            }
        }
        candidates = next;
        if candidates.len() > 1024 {
            return None;
        }
    }
    for (p_cand, _q_cand) in candidates {
        if p_cand > BigUint::one() && (n % &p_cand).is_zero() && &p_cand < n {
            return Some(p_cand);
        }
    }
    None
}

/// Invert `a` modulo `2^bits`.  Newton iteration: starts from `a` (odd) and doubles precision.
fn invert_mod_pow2(a: &BigUint, bits: usize) -> BigUint {
    let mut x = BigUint::one();
    let mask_init = BigUint::from(3u32); // mod 4: a is odd, a · 1 ≡ a (mod 2), a · a ≡ 1 (mod 8) for odd a actually
    let _ = mask_init;
    let mut prec = 1usize;
    while prec < bits {
        prec = (prec * 2).min(bits);
        let mask = BigUint::one() << prec;
        let ax = (a * &x) % &mask;
        let two = BigUint::from(2u32);
        let two_minus_ax = if two >= ax {
            two - ax
        } else {
            &two + &mask - ax
        };
        x = (&x * &two_minus_ax) % &mask;
    }
    x
}

pub fn run() -> Report {
    let mut r = Report::new(71, "Heninger-Shacham Cold-Boot RSA");
    let kp = crate::asymmetric::rsa::RsaKeyPair::generate(128);
    let n = &kp.public.n;
    let p = &kp.private.p;
    let q = &kp.private.q;
    let bits = p.bits() as usize;
    let mut rng = StdRng::seed_from_u64(71);
    let unknown_prob = 0.30;
    let partial_p = corrupt(p, bits, unknown_prob, &mut rng);
    let partial_q = corrupt(q, bits, unknown_prob, &mut rng);
    let known_p = partial_p.iter().filter(|b| b.is_some()).count();
    let known_q = partial_q.iter().filter(|b| b.is_some()).count();
    r.line(format!("N has {} bits, prime p has {} bits", n.bits(), bits));
    r.line(format!(
        "Known bits: p={}/{}  q={}/{}  (unknown_prob = {:.0}%)",
        known_p, bits, known_q, bits, unknown_prob * 100.0
    ));
    match cold_boot_recover_p(&partial_p, &partial_q, n, bits) {
        Some(recovered) => {
            r.line(format!("Recovered factor : {} bits", recovered.bits()));
            r.line(format!("Match            : {}", &recovered == p));
            if &recovered == p {
                return r.succeed();
            }
        }
        None => r.line("Branch-and-prune did not converge in budget."),
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cold_boot_recovers_small_p() {
        let kp = crate::asymmetric::rsa::RsaKeyPair::generate(96);
        let bits = kp.private.p.bits() as usize;
        let mut rng = StdRng::seed_from_u64(42);
        let partial_p = corrupt(&kp.private.p, bits, 0.30, &mut rng);
        let partial_q = corrupt(&kp.private.q, bits, 0.30, &mut rng);
        let recovered = cold_boot_recover_p(&partial_p, &partial_q, &kp.public.n, bits)
            .expect("must recover");
        assert_eq!(recovered, kp.private.p);
    }
}
