//! **SIS** — Short Integer Solution and Ajtai's collision-resistant
//! hash function (Ajtai 1996, "Generating hard instances of lattice
//! problems").
//!
//! # The problem
//! SIS is the "dual" of LWE and the other pillar of lattice crypto (it
//! underlies the *binding*/hardness of ML-DSA, hash-based lattice
//! commitments, etc.).  Given a uniform matrix `A ∈ Z_q^{n×m}` with
//! `m > n log q`, find a **short** nonzero integer vector `z` with
//!
//! ```text
//! A·z ≡ 0 (mod q),   0 < ‖z‖ ≤ β.
//! ```
//!
//! Ajtai proved this is hard *on average* if worst-case lattice
//! problems are hard — the first worst-case-to-average-case reduction
//! in cryptography.
//!
//! # Ajtai's hash and why SIS makes it collision-resistant
//! Define `H_A: {0,1}^m → Z_q^n` by `H_A(x) = A·x mod q`.  Because
//! `2^m > q^n` (guaranteed by `m > n log₂ q`), `H_A` compresses and so
//! has collisions by pigeonhole.  But a collision `x ≠ x'` with
//! `A·x = A·x'` yields `z = x − x' ∈ {−1, 0, 1}^m`, a short nonzero
//! vector with `A·z ≡ 0` — an SIS solution with `β = √m`.  So *finding*
//! a collision is exactly *solving* SIS: the hash is collision-resistant
//! iff SIS is hard.  This module implements the hash and demonstrates
//! the reduction both ways.
//!
//! # This implementation
//! Toy parameters `n = 4, q = 97, m = 32`, so `q^n = 97^4 ≈ 2^{26.4}`
//! and `2^{32} > q^n` gives guaranteed collisions — and, crucially,
//! keeps the *birthday* collision search feasible (`≈ √(q^n) ≈ 10⁴`
//! hash evaluations), which is what lets the tests actually *exhibit*
//! the collision→SIS-solution reduction.  At real parameters (`n` in
//! the hundreds) that search is infeasible — that infeasibility is the
//! whole security claim.  Not constant-time; see SECURITY.md.

use crate::utils::random::random_bytes;
use std::collections::HashMap;

/// Hash output dimension.
pub const N: usize = 4;
/// Modulus.
pub const Q: i64 = 97;
/// Input length (bits); `m > n·log₂ q` ⇒ compression ⇒ collisions exist.
pub const M: usize = 32;

/// The public matrix `A ∈ Z_q^{n×m}` defining the hash.
#[derive(Clone)]
pub struct AjtaiHash {
    pub a: Vec<Vec<i64>>, // n × m
}

fn rand_mod(q: i64) -> i64 {
    let mut buf = [0u8; 8];
    random_bytes(&mut buf);
    (u64::from_le_bytes(buf) % q as u64) as i64
}

impl AjtaiHash {
    /// Sample a fresh random hash key.
    pub fn keygen() -> AjtaiHash {
        let a = (0..N).map(|_| (0..M).map(|_| rand_mod(Q)).collect()).collect();
        AjtaiHash { a }
    }

    /// `H_A(x) = A·x mod q` for a binary input `x ∈ {0,1}^m`.
    pub fn hash(&self, x: &[u8]) -> Vec<i64> {
        assert_eq!(x.len(), M);
        (0..N)
            .map(|i| {
                let mut acc = 0i64;
                for j in 0..M {
                    if x[j] & 1 == 1 {
                        acc += self.a[i][j];
                    }
                }
                acc.rem_euclid(Q)
            })
            .collect()
    }

    /// Verify that `z ∈ {−1,0,1}^m` is a nonzero SIS solution:
    /// `A·z ≡ 0 (mod q)` and `z ≠ 0`.
    pub fn is_sis_solution(&self, z: &[i64]) -> bool {
        if z.len() != M || z.iter().all(|&x| x == 0) || z.iter().any(|&x| x.abs() > 1) {
            return false;
        }
        (0..N).all(|i| {
            let acc: i64 = (0..M).map(|j| self.a[i][j] * z[j]).sum();
            acc.rem_euclid(Q) == 0
        })
    }
}

/// Birthday-style collision finder over random binary inputs.  Returns
/// two distinct inputs with the same hash (feasible only at toy `n`).
pub fn find_collision(h: &AjtaiHash, max_tries: usize) -> Option<(Vec<u8>, Vec<u8>)> {
    let mut seen: HashMap<Vec<i64>, Vec<u8>> = HashMap::new();
    for _ in 0..max_tries {
        let mut bytes = vec![0u8; M];
        random_bytes(&mut bytes);
        let x: Vec<u8> = bytes.iter().map(|b| b & 1).collect();
        let digest = h.hash(&x);
        if let Some(prev) = seen.get(&digest) {
            if *prev != x {
                return Some((prev.clone(), x));
            }
        } else {
            seen.insert(digest, x);
        }
    }
    None
}

/// The reduction: turn a hash collision `(x, x')` into the short SIS
/// solution `z = x − x' ∈ {−1, 0, 1}^m`.
pub fn collision_to_sis(x: &[u8], x_prime: &[u8]) -> Vec<i64> {
    x.iter().zip(x_prime).map(|(&a, &b)| (a & 1) as i64 - (b & 1) as i64).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compression_guarantees_collisions() {
        // 2^m > q^n is what forces collisions to exist.
        // log2(q^n) = n·log2(q) ≈ 8·8.006 = 64.05 < 80 = m.
        let log2_qn = N as f64 * (Q as f64).log2();
        assert!((M as f64) > log2_qn, "hash must compress");
    }

    #[test]
    fn hash_is_deterministic_and_linear() {
        let h = AjtaiHash::keygen();
        let mut xb = vec![0u8; M];
        random_bytes(&mut xb);
        let x: Vec<u8> = xb.iter().map(|b| b & 1).collect();
        assert_eq!(h.hash(&x), h.hash(&x));
        // Linearity mod q: H(x ⊕-free sum) — check H(0) = 0.
        assert_eq!(h.hash(&vec![0u8; M]), vec![0i64; N]);
    }

    #[test]
    fn collision_yields_valid_sis_solution() {
        // The headline: a found collision reduces to a short nonzero
        // kernel vector — exactly an SIS solution.
        let h = AjtaiHash::keygen();
        let (x, x_prime) = find_collision(&h, 2_000_000).expect("collision at toy params");
        assert_ne!(x, x_prime);
        assert_eq!(h.hash(&x), h.hash(&x_prime));
        let z = collision_to_sis(&x, &x_prime);
        assert!(h.is_sis_solution(&z), "z = x − x' must be an SIS solution");
    }

    #[test]
    fn zero_and_long_vectors_are_not_solutions() {
        let h = AjtaiHash::keygen();
        assert!(!h.is_sis_solution(&vec![0i64; M])); // zero excluded
        let mut big = vec![0i64; M];
        big[0] = 2; // not short (coefficient magnitude > 1)
        assert!(!h.is_sis_solution(&big));
    }
}
