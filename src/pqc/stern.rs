//! **Stern identification / signature** — the original code-based
//! zero-knowledge protocol (Stern 1993).
//!
//! # Historical role
//! Stern's protocol is the ancestor of every code-based ZK signature,
//! including the MPC-in-the-head schemes in this library (`pqc::sdith`
//! and friends).  It proves knowledge of a **low-weight solution to a
//! syndrome-decoding instance** without revealing it, using a 3-move
//! commitment/challenge/response structure with a soundness error of
//! `2/3` per round (so it is repeated `≈ 1.7·λ` times, or Fiat–Shamir-
//! compiled into a signature).
//!
//! # The statement
//! Public: a parity-check matrix `H ∈ F₂^{(n−k)×n}` and a syndrome
//! `y = H·sᵀ` for a secret `s` of Hamming weight exactly `w`.  The
//! prover convinces the verifier it knows such an `s`.
//!
//! # One round (soundness error 2/3)
//! The prover picks a random permutation `π` of `[0,n)` and a random
//! masking vector `u ∈ F₂ⁿ`, and commits to three hashes:
//!
//! ```text
//! c₁ = H(π, H·uᵀ)          c₂ = H(π(u))          c₃ = H(π(u ⊕ s)).
//! ```
//!
//! The verifier sends a challenge `b ∈ {0,1,2}`:
//! - `b = 0`: reveal `(π, u)` — check `c₁, c₂`.
//! - `b = 1`: reveal `(π, u ⊕ s)` — check `c₁, c₃` (uses `H·sᵀ = y`).
//! - `b = 2`: reveal `(π(u), π(s))` — check `c₂, c₃` **and** that
//!   `wt(π(s)) = w`.
//!
//! Any single cheating prover can answer at most two of the three
//! challenges, hence soundness `2/3`.  Zero-knowledge: each response
//! reveals only permuted/masked data.
//!
//! # This implementation
//! Toy parameters `n = 64, k = 32, w = 12`, `τ = 40` Fiat–Shamir
//! rounds (soundness `≈ (2/3)^40 ≈ 2⁻²³`).  Commitments and challenges
//! use SHA3/SHAKE.  Not constant-time; see SECURITY.md.

use crate::hash::sha3::{sha3_256, shake256};
use crate::utils::random::random_bytes;

/// Code length.
pub const N: usize = 64;
/// Number of parity checks (`n − k`).
pub const R: usize = 32;
/// Secret Hamming weight.
pub const W: usize = 12;
/// Fiat–Shamir rounds.
pub const TAU: usize = 40;

/// Public key: parity-check matrix and target syndrome.
#[derive(Clone)]
pub struct SternPublicKey {
    pub h: Vec<Vec<u8>>, // R × N over F₂
    pub y: Vec<u8>,      // length R
}

/// Secret key: the weight-`w` preimage `s`.
#[derive(Clone)]
pub struct SternSecretKey {
    pub s: Vec<u8>,
}

fn matvec(h: &[Vec<u8>], x: &[u8]) -> Vec<u8> {
    h.iter().map(|row| row.iter().zip(x).fold(0u8, |a, (&hij, &xj)| a ^ (hij & xj))).collect()
}

fn xor(a: &[u8], b: &[u8]) -> Vec<u8> {
    a.iter().zip(b).map(|(&x, &y)| x ^ y).collect()
}

fn apply_perm(perm: &[usize], x: &[u8]) -> Vec<u8> {
    let mut out = vec![0u8; x.len()];
    for i in 0..x.len() {
        out[i] = x[perm[i]];
    }
    out
}

pub fn stern_keygen() -> (SternPublicKey, SternSecretKey) {
    // Random parity check.
    let mut h = vec![vec![0u8; N]; R];
    for row in h.iter_mut() {
        let mut b = vec![0u8; N];
        random_bytes(&mut b);
        for (c, byte) in row.iter_mut().zip(&b) {
            *c = byte & 1;
        }
    }
    // Weight-w secret.
    let mut s = vec![0u8; N];
    let mut placed = 0;
    while placed < W {
        let mut b = [0u8; 1];
        random_bytes(&mut b);
        let i = b[0] as usize % N;
        if s[i] == 0 {
            s[i] = 1;
            placed += 1;
        }
    }
    let y = matvec(&h, &s);
    (SternPublicKey { h, y }, SternSecretKey { s })
}

/// Deterministic random permutation of `[0, N)` from a seed.
fn perm_from_seed(seed: &[u8]) -> Vec<usize> {
    let stream = shake256(seed, N * 4);
    let mut perm: Vec<usize> = (0..N).collect();
    // Fisher–Yates with seed-derived randomness.
    for i in (1..N).rev() {
        let r = u32::from_le_bytes([
            stream[4 * i],
            stream[4 * i + 1],
            stream[4 * i + 2],
            stream[4 * i + 3],
        ]) as usize;
        perm.swap(i, r % (i + 1));
    }
    perm
}

fn commit(tag: u8, parts: &[&[u8]]) -> [u8; 32] {
    let mut buf = vec![tag];
    for p in parts {
        buf.extend_from_slice(p);
    }
    sha3_256(&buf)
}

fn perm_to_bytes(perm: &[usize]) -> Vec<u8> {
    perm.iter().flat_map(|&x| (x as u32).to_le_bytes()).collect()
}

/// Per-round data the prover keeps to answer the challenge.
struct RoundState {
    perm_seed: Vec<u8>,
    perm: Vec<usize>,
    u: Vec<u8>,
    c1: [u8; 32],
    c2: [u8; 32],
    c3: [u8; 32],
}

/// One round's opened response (union tagged by the challenge).
#[derive(Clone, Debug, PartialEq)]
pub struct RoundResponse {
    pub challenge: u8,
    /// For b=0/1: the permutation seed and a masked vector; for b=2: the
    /// two permuted vectors.
    pub vec_a: Vec<u8>,
    pub vec_b: Vec<u8>,
    pub perm_seed: Vec<u8>,
    /// The unopened commitment (the verifier recomputes the other two).
    pub unopened: [u8; 32],
}

#[derive(Clone, Debug, PartialEq)]
pub struct SternSignature {
    pub rounds: Vec<RoundResponse>,
}

fn build_round(pk: &SternPublicKey, sk: &SternSecretKey) -> RoundState {
    let mut perm_seed = vec![0u8; 32];
    random_bytes(&mut perm_seed);
    let perm = perm_from_seed(&perm_seed);
    let mut ub = vec![0u8; N];
    random_bytes(&mut ub);
    let u: Vec<u8> = ub.iter().map(|b| b & 1).collect();

    let hu = matvec(&pk.h, &u);
    let pu = apply_perm(&perm, &u);
    let pus = apply_perm(&perm, &xor(&u, &sk.s));

    let c1 = commit(1, &[&perm_to_bytes(&perm), &hu]);
    let c2 = commit(2, &[&pu]);
    let c3 = commit(3, &[&pus]);
    RoundState { perm_seed, perm, u, c1, c2, c3 }
}

/// Sign via the Fiat–Shamir transform of `τ` Stern rounds.
pub fn stern_sign(pk: &SternPublicKey, sk: &SternSecretKey, msg: &[u8]) -> SternSignature {
    let states: Vec<RoundState> = (0..TAU).map(|_| build_round(pk, sk)).collect();

    // Derive challenges from a hash of the message and all commitments.
    let mut ch_input = msg.to_vec();
    for st in &states {
        ch_input.extend_from_slice(&st.c1);
        ch_input.extend_from_slice(&st.c2);
        ch_input.extend_from_slice(&st.c3);
    }
    let ch_stream = shake256(&ch_input, TAU);

    let rounds = states
        .iter()
        .enumerate()
        .map(|(i, st)| {
            let b = ch_stream[i] % 3;
            match b {
                0 => RoundResponse {
                    challenge: 0,
                    vec_a: st.u.clone(),
                    vec_b: Vec::new(),
                    perm_seed: st.perm_seed.clone(),
                    unopened: st.c3,
                },
                1 => RoundResponse {
                    challenge: 1,
                    vec_a: xor(&st.u, &sk.s),
                    vec_b: Vec::new(),
                    perm_seed: st.perm_seed.clone(),
                    unopened: st.c2,
                },
                _ => RoundResponse {
                    challenge: 2,
                    vec_a: apply_perm(&st.perm, &st.u),
                    vec_b: apply_perm(&st.perm, &sk.s),
                    perm_seed: Vec::new(),
                    unopened: st.c1,
                },
            }
        })
        .collect();

    SternSignature { rounds }
}

/// Recompute the three commitments for a round from its response, or
/// `None` if the response is malformed (e.g. wrong weight at `b=2`).
fn recompute_commitments(
    pk: &SternPublicKey,
    resp: &RoundResponse,
) -> Option<([u8; 32], [u8; 32], [u8; 32])> {
    match resp.challenge {
        0 => {
            let perm = perm_from_seed(&resp.perm_seed);
            let u = &resp.vec_a;
            if u.len() != N {
                return None;
            }
            let c1 = commit(1, &[&perm_to_bytes(&perm), &matvec(&pk.h, u)]);
            let c2 = commit(2, &[&apply_perm(&perm, u)]);
            Some((c1, c2, resp.unopened))
        }
        1 => {
            let perm = perm_from_seed(&resp.perm_seed);
            let us = &resp.vec_a; // u ⊕ s
            if us.len() != N {
                return None;
            }
            // H·(u⊕s)ᵀ = H·uᵀ ⊕ y, so c1 is recoverable from us and y.
            let hus = xor(&matvec(&pk.h, us), &pk.y);
            let c1 = commit(1, &[&perm_to_bytes(&perm), &hus]);
            let c3 = commit(3, &[&apply_perm(&perm, us)]);
            Some((c1, resp.unopened, c3))
        }
        2 => {
            let pu = &resp.vec_a; // π(u)
            let ps = &resp.vec_b; // π(s)
            if pu.len() != N || ps.len() != N {
                return None;
            }
            // The soundness-carrying weight check.
            if ps.iter().filter(|&&x| x == 1).count() != W {
                return None;
            }
            let c2 = commit(2, &[pu]);
            let c3 = commit(3, &[&xor(pu, ps)]); // π(u) ⊕ π(s) = π(u ⊕ s)
            Some((resp.unopened, c2, c3))
        }
        _ => None,
    }
}

pub fn stern_verify(pk: &SternPublicKey, msg: &[u8], sig: &SternSignature) -> bool {
    if sig.rounds.len() != TAU {
        return false;
    }
    // Recompute every round's (c1, c2, c3), then re-derive the challenge
    // hash and check it selects exactly the openings the prover gave.
    let mut ch_input = msg.to_vec();
    let mut commitments = Vec::with_capacity(TAU);
    for resp in &sig.rounds {
        let (c1, c2, c3) = match recompute_commitments(pk, resp) {
            Some(c) => c,
            None => return false,
        };
        ch_input.extend_from_slice(&c1);
        ch_input.extend_from_slice(&c2);
        ch_input.extend_from_slice(&c3);
        commitments.push((c1, c2, c3));
    }
    let ch_stream = shake256(&ch_input, TAU);
    sig.rounds.iter().enumerate().all(|(i, resp)| ch_stream[i] % 3 == resp.challenge)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keygen_syndrome_is_consistent() {
        let (pk, sk) = stern_keygen();
        assert_eq!(sk.s.iter().filter(|&&x| x == 1).count(), W);
        assert_eq!(matvec(&pk.h, &sk.s), pk.y);
    }

    #[test]
    fn permutation_is_a_bijection() {
        let mut seed = vec![0u8; 32];
        random_bytes(&mut seed);
        let perm = perm_from_seed(&seed);
        let mut sorted = perm.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..N).collect::<Vec<_>>());
    }

    #[test]
    fn sign_verify_roundtrip() {
        let (pk, sk) = stern_keygen();
        let sig = stern_sign(&pk, &sk, b"prove knowledge of a low-weight word");
        assert!(stern_verify(&pk, b"prove knowledge of a low-weight word", &sig));
    }

    #[test]
    fn wrong_message_rejected() {
        let (pk, sk) = stern_keygen();
        let sig = stern_sign(&pk, &sk, b"one");
        assert!(!stern_verify(&pk, b"two", &sig));
    }

    #[test]
    fn tampered_response_rejected() {
        let (pk, sk) = stern_keygen();
        let msg = b"tamper";
        let mut sig = stern_sign(&pk, &sk, msg);
        // Flip a bit in a b=2 round's permuted secret (breaks the weight
        // check or the commitment).
        if let Some(r) = sig.rounds.iter_mut().find(|r| r.challenge == 2) {
            r.vec_b[0] ^= 1;
        }
        assert!(!stern_verify(&pk, msg, &sig));
    }

    #[test]
    fn forged_challenge_selection_rejected() {
        // Flipping a stored challenge without the matching openings must
        // fail the Fiat–Shamir consistency check.
        let (pk, sk) = stern_keygen();
        let msg = b"m";
        let mut sig = stern_sign(&pk, &sk, msg);
        sig.rounds[0].challenge = (sig.rounds[0].challenge + 1) % 3;
        assert!(!stern_verify(&pk, msg, &sig));
    }
}
