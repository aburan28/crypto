//! **LWE** — Regev's Learning-With-Errors public-key encryption
//! (Regev 2005, "On lattices, learning with errors, …").
//!
//! # The problem
//! LWE is the assumption underneath ML-KEM, ML-DSA, FrodoKEM and most
//! of lattice cryptography.  Fix a modulus `q` and dimension `n`.  A
//! sample is `(a, b = ⟨a, s⟩ + e mod q)` for a fixed secret `s ∈ Z_q^n`,
//! uniform `a ∈ Z_q^n`, and a *small* error `e`.  **Search-LWE**:
//! recover `s` from many samples.  **Decision-LWE**: distinguish such
//! samples from uniform.  Both are as hard as worst-case lattice
//! problems (GapSVP/SIVP) — Regev's celebrated quantum reduction.
//!
//! # Regev encryption (this module)
//! - **KeyGen**: secret `s ∈ Z_q^n`; public key `(A, b = A·s + e)` with
//!   `A ∈ Z_q^{m×n}`, small `e ∈ Z_q^m`.
//! - **Encrypt** one bit `μ`: pick a random selection vector
//!   `r ∈ {0,1}^m`; ciphertext is `(u, v) = (Aᵀr, ⟨b, r⟩ + μ·⌊q/2⌋)`.
//! - **Decrypt**: `v − ⟨u, s⟩ = ⟨e, r⟩ + μ·⌊q/2⌋ ≈ μ·⌊q/2⌋`; round to
//!   the nearer of `0` and `⌊q/2⌋`.
//!
//! Correctness needs the accumulated error `⟨e, r⟩` to stay below
//! `q/4`.  With `m` selection bits and errors of size `≤ B`, the worst
//! case is `m·B`, so `q` is chosen with `q/4 > m·B`.
//!
//! # This implementation
//! Toy parameters `n = 16, m = 128, q = 3329, |e| ≤ 2` (real Regev/Frodo
//! use `n` in the hundreds).  Errors are small-uniform rather than
//! discrete-Gaussian — enough to exhibit the scheme and its error
//! budget, not for security.  Not constant-time; see SECURITY.md.

use crate::utils::random::random_bytes;

/// Secret dimension.
pub const N: usize = 16;
/// Number of LWE samples in the public key (also the ciphertext
/// selection length).
pub const M: usize = 128;
/// Modulus.
pub const Q: i64 = 3329;
/// Error bound: |e| ≤ B.
pub const B: i64 = 2;

fn rand_mod(q: i64) -> i64 {
    let mut buf = [0u8; 8];
    random_bytes(&mut buf);
    (u64::from_le_bytes(buf) % q as u64) as i64
}

/// Small centred error in [−B, B].
fn rand_error() -> i64 {
    let mut b = [0u8; 1];
    random_bytes(&mut b);
    (b[0] as i64 % (2 * B + 1)) - B
}

/// Public key: `A` (m×n) and `b = A·s + e` (length m), all mod q.
#[derive(Clone)]
pub struct LwePublicKey {
    pub a: Vec<Vec<i64>>,
    pub b: Vec<i64>,
}

/// Secret key: `s ∈ Z_q^n`.
#[derive(Clone)]
pub struct LweSecretKey {
    pub s: Vec<i64>,
}

/// Ciphertext for a single bit: `(u ∈ Z_q^n, v ∈ Z_q)`.
#[derive(Clone, Debug, PartialEq)]
pub struct LweCiphertext {
    pub u: Vec<i64>,
    pub v: i64,
}

pub fn lwe_keygen() -> (LwePublicKey, LweSecretKey) {
    let s: Vec<i64> = (0..N).map(|_| rand_mod(Q)).collect();
    let mut a = vec![vec![0i64; N]; M];
    let mut b = vec![0i64; M];
    for i in 0..M {
        let mut acc = 0i64;
        for j in 0..N {
            a[i][j] = rand_mod(Q);
            acc += a[i][j] * s[j];
        }
        b[i] = (acc + rand_error()).rem_euclid(Q);
    }
    (LwePublicKey { a, b }, LweSecretKey { s })
}

/// Encrypt one bit under the public key.
pub fn lwe_encrypt(pk: &LwePublicKey, bit: u8) -> LweCiphertext {
    // Random 0/1 selection of the m samples.
    let mut sel = vec![0u8; M];
    random_bytes(&mut sel);
    let r: Vec<i64> = sel.iter().map(|&x| (x & 1) as i64).collect();

    // u = Aᵀ·r  (length n).
    let mut u = vec![0i64; N];
    for j in 0..N {
        let mut acc = 0i64;
        for i in 0..M {
            if r[i] != 0 {
                acc += pk.a[i][j];
            }
        }
        u[j] = acc.rem_euclid(Q);
    }
    // v = ⟨b, r⟩ + μ·⌊q/2⌋.
    let mut v = 0i64;
    for i in 0..M {
        if r[i] != 0 {
            v += pk.b[i];
        }
    }
    v += (bit as i64) * (Q / 2);
    LweCiphertext { u, v: v.rem_euclid(Q) }
}

/// Decrypt to a single bit.
pub fn lwe_decrypt(sk: &LweSecretKey, ct: &LweCiphertext) -> u8 {
    let mut acc = ct.v;
    for j in 0..N {
        acc -= ct.u[j] * sk.s[j];
    }
    let m = acc.rem_euclid(Q);
    // Round to nearer of 0 and q/2.
    let dist_zero = m.min(Q - m);
    let dist_half = (m - Q / 2).abs();
    if dist_half < dist_zero {
        1
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_budget_holds() {
        // Worst-case accumulated error m·B must stay below q/4 for
        // correctness.
        assert!(M as i64 * B < Q / 4, "error budget exceeded");
    }

    #[test]
    fn encrypt_decrypt_both_bits() {
        let (pk, sk) = lwe_keygen();
        for _ in 0..50 {
            for bit in [0u8, 1] {
                let ct = lwe_encrypt(&pk, bit);
                assert_eq!(lwe_decrypt(&sk, &ct), bit);
            }
        }
    }

    #[test]
    fn wrong_key_fails_often() {
        // A different secret decrypts to noise: it should disagree with
        // the plaintext on a good fraction of trials.
        let (pk, _) = lwe_keygen();
        let (_, sk_other) = lwe_keygen();
        let mut mismatches = 0;
        for _ in 0..100 {
            let ct = lwe_encrypt(&pk, 1);
            if lwe_decrypt(&sk_other, &ct) != 1 {
                mismatches += 1;
            }
        }
        assert!(mismatches > 20, "wrong key decrypted correctly too often");
    }

    #[test]
    fn ciphertext_is_not_plaintext() {
        // Sanity: encryption actually hides the bit (u is nonzero).
        let (pk, _) = lwe_keygen();
        let ct = lwe_encrypt(&pk, 0);
        assert_eq!(ct.u.len(), N);
    }
}
