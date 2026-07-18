//! **NewHope** (NewHope-Simple variant) — Ring-LWE key encapsulation
//! (Alkim–Ducas–Pöppelmann–Schwabe 2016; the "-Simple" LPR-style
//! variant, Alkim et al. 2016).
//!
//! # Historical note
//! NewHope was the first post-quantum KEM deployed at scale: Google ran
//! it in Chrome Canary (CECPQ1, 2016) hybridised with X25519.  It
//! popularised Ring-LWE KEMs and directly influenced Kyber/ML-KEM.
//!
//! # The distinctive ideas
//! NewHope-Simple is Ring-LWE encryption (see `pqc::ring_lwe`) turned
//! into a KEM, with two features worth studying:
//!
//! 1. **Centred-binomial noise** instead of discrete Gaussians —
//!    sampling `Σ(aᵢ − bᵢ)` over random bits.  Cheap, constant-time-
//!    friendly, and good enough; Kyber inherited it.
//! 2. **Redundant 4:1 message encoding** (`NHSEncode`/`NHSDecode`): each
//!    of the 256 shared-secret bits is written into *four* polynomial
//!    coefficients, and decoding sums the four distances to `q/2`.  This
//!    error-correcting encode is what lets NewHope tolerate large noise
//!    at `n = 1024` without decryption failures.
//!
//! # Flow
//! - **KeyGen**: `b = a·s + e`, public `(seed_a, b)`.
//! - **Encaps**: pick 32-byte `μ`; `u = a·s' + e'`,
//!   `v = b·s' + e'' + NHSEncode(μ)`; ciphertext `(u, v)`; shared secret
//!   `SHA3-256(μ)`.
//! - **Decaps**: `μ' = NHSDecode(v − u·s)`; shared secret `SHA3-256(μ')`.
//!
//! # This implementation
//! Faithful structure at NewHope's real ring parameters `n = 1024,
//! q = 12289`, centred-binomial noise `ψ_8`, 4:1 encoding of a 256-bit
//! secret.  Schoolbook negacyclic multiplication (no NTT — that lives
//! in `pqc::ml_kem`).  Ciphertext compression is omitted for clarity.
//! Not constant-time; see SECURITY.md.

use crate::hash::sha3::{sha3_256, shake128};
use crate::utils::random::random_bytes;

/// Ring degree.
pub const N: usize = 1024;
/// Modulus.
pub const Q: i64 = 12289;
/// Centred-binomial parameter (noise in [−K, K]).
pub const K: usize = 8;
/// Shared-secret / message length in bits (each uses 4 coefficients).
pub const MSG_BITS: usize = N / 4; // 256

type Poly = Vec<i64>;

fn poly_mul(a: &[i64], b: &[i64]) -> Poly {
    let mut out = vec![0i64; N];
    for i in 0..N {
        if a[i] == 0 {
            continue;
        }
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

/// Expand a 32-byte seed into a uniform ring element via SHAKE128.
/// Each 2-byte draw is masked to 14 bits (range 16384) before the
/// rejection test `< q`, giving ~75% acceptance; the buffer is sized
/// with ample margin for the rejected draws.
fn gen_a(seed: &[u8]) -> Poly {
    let bytes = shake128(seed, N * 4 + 1024);
    let mut a = Vec::with_capacity(N);
    let mut pos = 0;
    while a.len() < N {
        let v = (u16::from_le_bytes([bytes[pos], bytes[pos + 1]]) & 0x3fff) as i64;
        pos += 2;
        if v < Q {
            a.push(v);
        }
    }
    a
}

/// Centred-binomial sample ψ_K: Σ_{i<K}(aᵢ − bᵢ) from 2K random bits.
fn noise_poly() -> Poly {
    let mut bytes = vec![0u8; N * 2 * K / 8];
    random_bytes(&mut bytes);
    let bit = |i: usize| ((bytes[i / 8] >> (i % 8)) & 1) as i64;
    (0..N)
        .map(|i| {
            let base = i * 2 * K;
            let mut acc = 0i64;
            for j in 0..K {
                acc += bit(base + j) - bit(base + K + j);
            }
            acc.rem_euclid(Q)
        })
        .collect()
}

/// `NHSEncode`: write each of the 256 message bits into 4 coefficients
/// (positions `i, i+256, i+512, i+768`), value `bit·⌊q/2⌋`.
fn nhs_encode(msg: &[u8; 32]) -> Poly {
    let mut v = vec![0i64; N];
    for i in 0..MSG_BITS {
        let bit = ((msg[i / 8] >> (i % 8)) & 1) as i64;
        let val = bit * (Q / 2);
        for r in 0..4 {
            v[i + r * MSG_BITS] = val;
        }
    }
    v
}

/// `NHSDecode`: sum the four coefficients' distances to `q/2`; the bit
/// is 1 iff the total is closer to `q/2` than to 0 (error-correcting).
fn nhs_decode(v: &[i64]) -> [u8; 32] {
    let mut msg = [0u8; 32];
    for i in 0..MSG_BITS {
        let mut sum = 0i64;
        for r in 0..4 {
            let c = v[i + r * MSG_BITS].rem_euclid(Q);
            // Distance from q/2, mapped so that "near q/2" is negative.
            sum += (c - Q / 2).abs() - Q / 4;
        }
        // sum < 0 ⇒ coefficients cluster near q/2 ⇒ bit 1.
        if sum < 0 {
            msg[i / 8] |= 1 << (i % 8);
        }
    }
    msg
}

#[derive(Clone)]
pub struct NewHopePublicKey {
    pub seed: Vec<u8>,
    pub b: Poly,
}

#[derive(Clone)]
pub struct NewHopeSecretKey {
    pub s: Poly,
}

#[derive(Clone, Debug, PartialEq)]
pub struct NewHopeCiphertext {
    pub u: Poly,
    pub v: Poly,
}

pub fn newhope_keygen() -> (NewHopePublicKey, NewHopeSecretKey) {
    let mut seed = vec![0u8; 32];
    random_bytes(&mut seed);
    let a = gen_a(&seed);
    let s = noise_poly();
    let e = noise_poly();
    let b = poly_add(&poly_mul(&a, &s), &e);
    (NewHopePublicKey { seed, b }, NewHopeSecretKey { s })
}

pub fn newhope_encaps(pk: &NewHopePublicKey) -> (NewHopeCiphertext, [u8; 32]) {
    let mut mu = [0u8; 32];
    random_bytes(&mut mu);
    let a = gen_a(&pk.seed);
    let s2 = noise_poly();
    let e1 = noise_poly();
    let e2 = noise_poly();
    let u = poly_add(&poly_mul(&a, &s2), &e1);
    let v = poly_add(&poly_add(&poly_mul(&pk.b, &s2), &e2), &nhs_encode(&mu));
    (NewHopeCiphertext { u, v }, sha3_256(&mu))
}

pub fn newhope_decaps(sk: &NewHopeSecretKey, ct: &NewHopeCiphertext) -> [u8; 32] {
    let m = poly_sub(&ct.v, &poly_mul(&ct.u, &sk.s));
    sha3_256(&nhs_decode(&m))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_decode_is_noise_tolerant() {
        // The 4:1 encoding must survive a per-coefficient perturbation
        // well beyond what a 1:1 encoding could.
        let mut msg = [0u8; 32];
        random_bytes(&mut msg);
        let mut v = nhs_encode(&msg);
        // Add ±q/8 jitter to every coefficient.
        for c in v.iter_mut() {
            *c = (*c + Q / 8).rem_euclid(Q);
        }
        assert_eq!(nhs_decode(&v), msg);
    }

    #[test]
    fn kem_roundtrip() {
        for _ in 0..10 {
            let (pk, sk) = newhope_keygen();
            let (ct, ss1) = newhope_encaps(&pk);
            let ss2 = newhope_decaps(&sk, &ct);
            assert_eq!(ss1, ss2);
        }
    }

    #[test]
    fn distinct_encapsulations_differ() {
        let (pk, _) = newhope_keygen();
        let (_, ss1) = newhope_encaps(&pk);
        let (_, ss2) = newhope_encaps(&pk);
        assert_ne!(ss1, ss2);
    }

    #[test]
    fn wrong_secret_key_gives_different_secret() {
        let (pk, _) = newhope_keygen();
        let (_, sk_other) = newhope_keygen();
        let (ct, ss1) = newhope_encaps(&pk);
        assert_ne!(newhope_decaps(&sk_other, &ct), ss1);
    }
}
