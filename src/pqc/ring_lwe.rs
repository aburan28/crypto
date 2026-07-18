//! **Ring-LWE** — Lyubashevsky–Peikert–Regev public-key encryption
//! (LPR 2010, "On ideal lattices and learning with errors over rings").
//!
//! # Why rings
//! Plain LWE (see `pqc::lwe`) needs an `m×n` matrix — quadratic key
//! size.  Ring-LWE replaces the matrix with a *single ring element*:
//! work in `R_q = Z_q[x]/(xⁿ + 1)`, where multiplication by a fixed
//! `a ∈ R_q` already acts like an `n×n` structured (negacyclic) matrix.
//! One ring element carries `n` scalars, so keys and ciphertexts shrink
//! by a factor of `n`.  This is the structure ML-KEM/Kyber and NewHope
//! are built on; Ring-LWE is the direct ancestor.
//!
//! # LPR encryption (this module)
//! - **KeyGen**: public `a ∈ R_q` (uniform); secret `s`, error `e` with
//!   small coefficients; public key `b = a·s + e`.
//! - **Encrypt** an n-bit message `m ∈ {0,1}ⁿ` (encoded as `⌊q/2⌋·m`):
//!   sample small `r, e₁, e₂`; ciphertext `(u, v) =
//!   (a·r + e₁, b·r + e₂ + ⌊q/2⌋·m)`.
//! - **Decrypt**: `v − u·s = ⌊q/2⌋·m + (e·r + e₂ − e₁·s)`; round each
//!   coefficient to `0` or `⌊q/2⌋`.  Correct while the noise term stays
//!   below `q/4` per coefficient.
//!
//! # This implementation
//! Toy parameters `n = 64, q = 12289, small coeffs in {−1,0,1}`
//! (real Ring-LWE/Kyber use `n = 256`).  Schoolbook negacyclic
//! multiplication for clarity (no NTT — that optimisation lives in
//! `pqc::ml_kem`).  Not constant-time; see SECURITY.md.

use crate::utils::random::random_bytes;

/// Ring degree: R_q = Z_q[x]/(xⁿ + 1).
pub const N: usize = 64;
/// Modulus (NTT-friendly prime, though we use schoolbook mult here).
pub const Q: i64 = 12289;

type Poly = Vec<i64>;

fn rand_mod(q: i64) -> i64 {
    let mut buf = [0u8; 8];
    random_bytes(&mut buf);
    (u64::from_le_bytes(buf) % q as u64) as i64
}

/// Uniform ring element.
fn rand_poly() -> Poly {
    (0..N).map(|_| rand_mod(Q)).collect()
}

/// Small ternary ring element with coefficients in {−1, 0, 1}.
fn small_poly() -> Poly {
    let mut b = vec![0u8; N];
    random_bytes(&mut b);
    b.iter().map(|&x| (x % 3) as i64 - 1).collect()
}

/// Negacyclic product in R_q: `xⁿ = −1`.
fn poly_mul(a: &[i64], b: &[i64]) -> Poly {
    let mut out = vec![0i64; N];
    for i in 0..N {
        for j in 0..N {
            let k = i + j;
            let t = a[i] * b[j];
            if k < N {
                out[k] = (out[k] + t).rem_euclid(Q);
            } else {
                out[k - N] = (out[k - N] - t).rem_euclid(Q);
            }
        }
    }
    out
}

fn poly_add(a: &[i64], b: &[i64]) -> Poly {
    (0..N).map(|i| (a[i] + b[i]).rem_euclid(Q)).collect()
}

fn poly_sub(a: &[i64], b: &[i64]) -> Poly {
    (0..N).map(|i| (a[i] - b[i]).rem_euclid(Q)).collect()
}

/// Public key: `(a, b = a·s + e)`.
#[derive(Clone)]
pub struct RingLwePublicKey {
    pub a: Poly,
    pub b: Poly,
}

#[derive(Clone)]
pub struct RingLweSecretKey {
    pub s: Poly,
}

/// Ciphertext `(u, v)`, each a ring element.
#[derive(Clone, Debug, PartialEq)]
pub struct RingLweCiphertext {
    pub u: Poly,
    pub v: Poly,
}

pub fn ring_lwe_keygen() -> (RingLwePublicKey, RingLweSecretKey) {
    let a = rand_poly();
    let s = small_poly();
    let e = small_poly();
    let b = poly_add(&poly_mul(&a, &s), &e);
    (RingLwePublicKey { a, b }, RingLweSecretKey { s })
}

/// Encode an n-bit message into a ring element scaled by ⌊q/2⌋.
fn encode(msg: &[u8]) -> Poly {
    (0..N).map(|i| (msg[i] as i64 & 1) * (Q / 2)).collect()
}

/// Encrypt an n-bit message (`msg[i] ∈ {0,1}`).
pub fn ring_lwe_encrypt(pk: &RingLwePublicKey, msg: &[u8]) -> RingLweCiphertext {
    assert_eq!(msg.len(), N);
    let r = small_poly();
    let e1 = small_poly();
    let e2 = small_poly();
    let u = poly_add(&poly_mul(&pk.a, &r), &e1);
    let v = poly_add(&poly_add(&poly_mul(&pk.b, &r), &e2), &encode(msg));
    RingLweCiphertext { u, v }
}

/// Decrypt to an n-bit message.
pub fn ring_lwe_decrypt(sk: &RingLweSecretKey, ct: &RingLweCiphertext) -> Vec<u8> {
    let m = poly_sub(&ct.v, &poly_mul(&ct.u, &sk.s));
    m.iter()
        .map(|&c| {
            let c = c.rem_euclid(Q);
            let dist_zero = c.min(Q - c);
            let dist_half = (c - Q / 2).abs();
            if dist_half < dist_zero {
                1u8
            } else {
                0u8
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_message() -> Vec<u8> {
        let mut b = vec![0u8; N];
        random_bytes(&mut b);
        b.iter().map(|x| x & 1).collect()
    }

    #[test]
    fn negacyclic_wraparound() {
        // x^{N-1} · x = x^N = −1.
        let mut a = vec![0i64; N];
        a[N - 1] = 1;
        let mut b = vec![0i64; N];
        b[1] = 1;
        let p = poly_mul(&a, &b);
        assert_eq!(p[0], (Q - 1));
        assert!(p[1..].iter().all(|&c| c == 0));
    }

    #[test]
    fn encrypt_decrypt_roundtrip() {
        let (pk, sk) = ring_lwe_keygen();
        for _ in 0..30 {
            let msg = random_message();
            let ct = ring_lwe_encrypt(&pk, &msg);
            assert_eq!(ring_lwe_decrypt(&sk, &ct), msg);
        }
    }

    #[test]
    fn all_zero_and_all_one_messages() {
        let (pk, sk) = ring_lwe_keygen();
        for bit in [0u8, 1] {
            let msg = vec![bit; N];
            let ct = ring_lwe_encrypt(&pk, &msg);
            assert_eq!(ring_lwe_decrypt(&sk, &ct), msg);
        }
    }

    #[test]
    fn wrong_secret_key_fails() {
        let (pk, _) = ring_lwe_keygen();
        let (_, sk_other) = ring_lwe_keygen();
        let msg = vec![1u8; N];
        let ct = ring_lwe_encrypt(&pk, &msg);
        assert_ne!(ring_lwe_decrypt(&sk_other, &ct), msg);
    }
}
