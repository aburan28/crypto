//! **Saber** — Module Learning-With-Rounding key encapsulation
//! (D'Anvers–Karmakar–Sinha Roy–Vercauteren 2018; NIST round-3
//! finalist).
//!
//! # The distinctive idea: rounding instead of noise
//! ML-KEM/Kyber add an explicit error `e` sampled from a distribution.
//! Saber adds *no* error — it uses **Learning With Rounding (LWR)**:
//! compute `A·s` over a large modulus `q`, then throw away the low bits
//! by rounding to a smaller modulus `p`.  The discarded bits *are* the
//! error, deterministically.  This means:
//!
//! - **Power-of-two moduli** (`q = 2¹³, p = 2¹⁰, T = 2⁴`): reduction is
//!   just bit-shifting — no modular reduction, no NTT-friendly prime.
//!   (Multiplication uses Toom–Cook/Karatsuba in real Saber; we use
//!   schoolbook.)
//! - **No noise sampling** on the "rounding" side — smaller, simpler.
//! - Rounding is a public deterministic function, so the LWR assumption
//!   (rather than LWE) carries the security.
//!
//! # Flow (module rank `l`, ring `R = Z[x]/(xⁿ+1)`)
//! - **KeyGen**: `b = Round_p(Aᵀ·s)`, secret `s` (small, binomial).
//! - **Encaps**: `μ` a 256-bit message; `b' = Round_p(A·s')`;
//!   `v' = bᵀ·s'` (mod p); `c_m = Round_T(v' − ⌊p/2⌋·μ)`; ciphertext
//!   `(b', c_m)`; shared secret `SHA3-256(μ)`.
//! - **Decaps**: `v = b'ᵀ·s` (mod p);
//!   `μ' = Round_2(v − 2^{ep−eT}·c_m)`; secret `SHA3-256(μ')`.
//!
//! # This implementation
//! LightSaber-like toy params `n = 256, l = 2, q = 2¹³, p = 2¹⁰,
//! T = 2⁴`, binomial secret `ψ_5`.  Faithful to the LWR structure and
//! the rounding constants `h1, h2` that centre the rounding error.
//! Schoolbook multiplication; not constant-time; see SECURITY.md.

use crate::hash::sha3::{sha3_256, shake128};
use crate::utils::random::random_bytes;

/// Ring degree.
pub const N: usize = 256;
/// Module rank (LightSaber = 2).
pub const L: usize = 2;
/// log2 of the large modulus q.
pub const EQ: u32 = 13;
/// log2 of the rounded modulus p.
pub const EP: u32 = 10;
/// log2 of the ciphertext modulus T.
pub const ET: u32 = 4;
/// Binomial parameter for the secret.
pub const MU: usize = 5;

const Q: i64 = 1 << EQ;
const P: i64 = 1 << EP;
const T: i64 = 1 << ET;

type Poly = Vec<i64>;
type PolyVec = Vec<Poly>;

/// Negacyclic multiply mod 2^EQ (bit-mask reduction — no prime needed).
fn poly_mul(a: &[i64], b: &[i64]) -> Poly {
    let mask = Q - 1;
    let mut out = vec![0i64; N];
    for i in 0..N {
        if a[i] == 0 {
            continue;
        }
        for j in 0..N {
            let k = i + j;
            let t = a[i] * b[j];
            if k < N {
                out[k] = (out[k] + t) & mask;
            } else {
                out[k - N] = (out[k - N] - t) & mask;
            }
        }
    }
    out.iter().map(|&x| x.rem_euclid(Q)).collect()
}

fn poly_add_mod(a: &[i64], b: &[i64], m: i64) -> Poly {
    (0..a.len()).map(|i| (a[i] + b[i]).rem_euclid(m)).collect()
}

/// Expand seed → l×l matrix of uniform R_q polys.
fn gen_matrix(seed: &[u8]) -> Vec<PolyVec> {
    let bytes = shake128(seed, L * L * N * 2 + 64);
    let mut pos = 0;
    let mut next = || {
        let v = u16::from_le_bytes([bytes[pos], bytes[pos + 1]]) as i64 & (Q - 1);
        pos += 2;
        v
    };
    (0..L).map(|_| (0..L).map(|_| (0..N).map(|_| next()).collect()).collect()).collect()
}

/// Centred-binomial secret polynomial ψ_MU.
fn secret_poly() -> Poly {
    let mut bytes = vec![0u8; N * 2 * MU / 8 + 1];
    random_bytes(&mut bytes);
    let bit = |i: usize| ((bytes[i / 8] >> (i % 8)) & 1) as i64;
    (0..N)
        .map(|i| {
            let base = i * 2 * MU;
            let mut acc = 0i64;
            for j in 0..MU {
                acc += bit(base + j) - bit(base + MU + j);
            }
            acc
        })
        .collect()
}

fn secret_vec() -> PolyVec {
    (0..L).map(|_| secret_poly()).collect()
}

/// Round from modulus 2^from down to 2^to: add the half-ULP centring
/// constant, then shift.
fn round_shift(x: i64, from: u32, to: u32) -> i64 {
    let shift = from - to;
    let h = 1i64 << (shift - 1);
    ((x + h) >> shift) & ((1 << to) - 1)
}

/// Round every coefficient of a poly from 2^from to 2^to.
fn round_poly(p: &[i64], from: u32, to: u32) -> Poly {
    p.iter().map(|&x| round_shift(x, from, to)).collect()
}

/// Matrixᵀ · vec over R_q (used with rounding afterwards).
fn mat_transpose_vec(a: &[PolyVec], s: &PolyVec) -> PolyVec {
    (0..L)
        .map(|i| {
            let mut acc = vec![0i64; N];
            for j in 0..L {
                acc = poly_add_mod(&acc, &poly_mul(&a[j][i], &s[j]), Q);
            }
            acc
        })
        .collect()
}

fn mat_vec(a: &[PolyVec], s: &PolyVec) -> PolyVec {
    (0..L)
        .map(|i| {
            let mut acc = vec![0i64; N];
            for j in 0..L {
                acc = poly_add_mod(&acc, &poly_mul(&a[i][j], &s[j]), Q);
            }
            acc
        })
        .collect()
}

/// Inner product bᵀ·s over R_p.
fn inner_p(b: &[Poly], s: &PolyVec) -> Poly {
    let mut acc = vec![0i64; N];
    for i in 0..L {
        acc = poly_add_mod(&acc, &poly_mul(&b[i], &s[i]), P);
    }
    acc
}

#[derive(Clone)]
pub struct SaberPublicKey {
    pub seed: Vec<u8>,
    pub b: Vec<Poly>, // l polys in R_p
}

#[derive(Clone)]
pub struct SaberSecretKey {
    pub s: PolyVec,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SaberCiphertext {
    pub b_prime: Vec<Poly>, // l polys in R_p
    pub cm: Poly,           // 1 poly in R_T
}

pub fn saber_keygen() -> (SaberPublicKey, SaberSecretKey) {
    let mut seed = vec![0u8; 32];
    random_bytes(&mut seed);
    let a = gen_matrix(&seed);
    let s = secret_vec();
    // b = Round_p(Aᵀ·s).
    let as_ = mat_transpose_vec(&a, &s);
    let b: Vec<Poly> = as_.iter().map(|p| round_poly(p, EQ, EP)).collect();
    (SaberPublicKey { seed, b }, SaberSecretKey { s })
}

fn encode_msg(msg: &[u8; 32]) -> Poly {
    // Each of 256 bits → one coefficient, scaled by p/2.
    (0..N).map(|i| ((msg[i / 8] >> (i % 8)) & 1) as i64 * (P / 2)).collect()
}

pub fn saber_encaps(pk: &SaberPublicKey) -> (SaberCiphertext, [u8; 32]) {
    let mut mu = [0u8; 32];
    random_bytes(&mut mu);
    let a = gen_matrix(&pk.seed);
    let s2 = secret_vec();
    // b' = Round_p(A·s').
    let as2 = mat_vec(&a, &s2);
    let b_prime: Vec<Poly> = as2.iter().map(|p| round_poly(p, EQ, EP)).collect();
    // v' = bᵀ·s' over R_p; c_m = Round_T(v' − (p/2)·μ + h1).
    let vp = inner_p(&pk.b, &s2);
    let msg = encode_msg(&mu);
    let diff: Poly = (0..N).map(|i| (vp[i] - msg[i]).rem_euclid(P)).collect();
    let cm = round_poly(&diff, EP, ET);
    (SaberCiphertext { b_prime, cm }, sha3_256(&mu))
}

pub fn saber_decaps(sk: &SaberSecretKey, ct: &SaberCiphertext) -> [u8; 32] {
    // v = b'ᵀ·s over R_p.
    let v = inner_p(&ct.b_prime, &sk.s);
    // Lift c_m back to R_p and subtract: m ≈ v − 2^{ep−et}·c_m.
    let mut msg = [0u8; 32];
    for i in 0..N {
        let cm_lifted = ct.cm[i] << (EP - ET);
        let d = (v[i] - cm_lifted).rem_euclid(P);
        // Round to bit: 1 if d is closer to p/2 than to 0/p.
        let dist_zero = d.min(P - d);
        let dist_half = (d - P / 2).abs();
        if dist_half < dist_zero {
            msg[i / 8] |= 1 << (i % 8);
        }
    }
    sha3_256(&msg)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rounding_is_bit_shift() {
        // Round_p(x) drops the low EQ−EP bits (with centring).
        assert_eq!(round_shift(0, EQ, EP), 0);
        assert_eq!(round_shift((1 << (EQ - EP)) / 2, EQ, EP), 1); // half-ULP rounds up
        assert_eq!(round_shift(Q - 1, EQ, EP) & (P - 1), 0); // wraps to 0 mod p
    }

    #[test]
    fn no_explicit_error_in_keygen() {
        // Saber's premise: keygen samples only the secret, never a noise
        // poly — the rounding *is* the error. b lives in R_p.
        let (pk, _) = saber_keygen();
        assert_eq!(pk.b.len(), L);
        assert!(pk.b.iter().all(|p| p.iter().all(|&c| c >= 0 && c < P)));
    }

    #[test]
    fn kem_roundtrip() {
        for _ in 0..20 {
            let (pk, sk) = saber_keygen();
            let (ct, ss1) = saber_encaps(&pk);
            let ss2 = saber_decaps(&sk, &ct);
            assert_eq!(ss1, ss2);
        }
    }

    #[test]
    fn ciphertext_moduli_are_respected() {
        let (pk, _) = saber_keygen();
        let (ct, _) = saber_encaps(&pk);
        assert!(ct.b_prime.iter().all(|p| p.iter().all(|&c| c < P)));
        assert!(ct.cm.iter().all(|&c| c < T));
    }

    #[test]
    fn wrong_secret_key_gives_different_secret() {
        let (pk, _) = saber_keygen();
        let (_, sk_other) = saber_keygen();
        let (ct, ss1) = saber_encaps(&pk);
        assert_ne!(saber_decaps(&sk_other, &ct), ss1);
    }
}
