//! **Picnic** — signatures from a block cipher in the head
//! (Chase et al. 2017; NIST round-3 *alternate* candidate), built on
//! the **LowMC** cipher and this library's MPC-in-the-head engine.
//!
//! # The idea (cf. FAEST)
//! Picnic proves knowledge of a secret key `k` with `LowMC_k(p) = c`
//! (public plaintext/ciphertext) without revealing `k`, via
//! MPC-in-the-head (`pqc::mpcith`) — exactly like `pqc::faest`, but with
//! **LowMC** instead of AES.  LowMC (Albrecht et al. 2015) was designed
//! for MPC/FHE/ZK: it has a **minimal number of AND gates**, because in
//! these "in the head" proofs the cost is dominated by the nonlinear
//! (multiplication) gates.  It achieves this with a *partial* S-box
//! layer (only the top `3m` bits are S-boxed each round) and heavy,
//! random-looking `F₂`-linear layers.
//!
//! # LowMC and its AND gates
//! Each round applies `m` copies of the 3-bit S-box
//! ```text
//! S(a,b,c) = (a ⊕ b·c,  a ⊕ b ⊕ a·c,  a ⊕ b ⊕ c ⊕ a·b),
//! ```
//! then a linear layer `L_r`, a round constant `C_r`, and a round key
//! `K_r·k`.  Its only nonlinear gates are the three ANDs per S-box.
//! The **extended witness** is `k` together with every S-box *output*;
//! given those, the whole cipher is `F₂`-linear, so each MPC party
//! recomputes it from its share.  The AND relations reduce to
//! `b·c = d⊕a`, `a·c = e⊕a⊕b`, `a·b = f⊕a⊕b⊕c` (with `(d,e,f)` the
//! S-box output), and all are batched — under verifier randomness `ε` —
//! into the single sacrificed-triple dot-product check that the
//! `mpcith` engine provides.
//!
//! # This implementation
//! Toy LowMC `n = 16` (block = key), `m = 2` S-boxes/round, `R = 5`
//! rounds; the public matrices/constants are expanded from a fixed seed
//! with SHAKE (as in the LowMC spec).  Faithful to the partial-S-box
//! structure and the AND-gate proof; plus the engine simplifications in
//! `pqc::mpcith`.  Toy scale, not constant-time; see SECURITY.md.

use super::mpcith::{gf_mul, mpcith_prove, mpcith_verify, MpcRelation, MpcithProof, PartyView};
use crate::hash::sha3::{shake256};
use crate::utils::random::random_bytes;

/// Block / key size in bits.
pub const N: usize = 16;
/// S-boxes per round (partial layer covers the top `3·m` bits).
pub const M: usize = 2;
/// Rounds.
pub const R: usize = 5;
/// Bits touched by the S-box layer.
const SB_BITS: usize = 3 * M; // 6
/// AND gates total = 3 per S-box.
const GATES: usize = 3 * M * R; // 30
/// Extended witness = key ‖ all S-box outputs.
const WITNESS: usize = N + SB_BITS * R; // 16 + 30 = 46

// ── Public LowMC parameters (expanded from a fixed seed) ─────────────────────

/// Linear layers `L_r`, key matrices `K_r`, and round constants `C_r`.
struct LowMcParams {
    lin: Vec<Vec<Vec<u8>>>, // R matrices, each N×N over F₂
    key: Vec<Vec<Vec<u8>>>, // R+1 matrices (whitening + per round), N×N
    rc: Vec<Vec<u8>>,       // R constants, each N bits
}

fn expand_params() -> LowMcParams {
    // Deterministic public parameters (every signer/verifier agrees).
    let stream = shake256(b"Picnic-LowMC-toy-params", (R + (R + 1)) * N * N / 8 + R * N + 64);
    let mut pos = 0usize;
    let mut bit = || {
        let b = (stream[pos / 8] >> (pos % 8)) & 1;
        pos += 1;
        b
    };
    let mut matrix = |bitref: &mut dyn FnMut() -> u8| {
        (0..N).map(|_| (0..N).map(|_| bitref()).collect()).collect::<Vec<Vec<u8>>>()
    };
    let lin = (0..R).map(|_| matrix(&mut bit)).collect();
    let key = (0..R + 1).map(|_| matrix(&mut bit)).collect();
    let rc = (0..R).map(|_| (0..N).map(|_| bit()).collect()).collect();
    LowMcParams { lin, key, rc }
}

fn mat_vec(m: &[Vec<u8>], v: &[u8]) -> Vec<u8> {
    m.iter().map(|row| row.iter().zip(v).fold(0u8, |a, (&r, &x)| a ^ (r & x)) & 1).collect()
}

fn xor(a: &[u8], b: &[u8]) -> Vec<u8> {
    a.iter().zip(b).map(|(&x, &y)| x ^ y).collect()
}

/// The 3-bit LowMC S-box.
fn sbox(a: u8, b: u8, c: u8) -> (u8, u8, u8) {
    (a ^ (b & c), a ^ b ^ (a & c), a ^ b ^ c ^ (a & b))
}

/// Run LowMC, returning the ciphertext and the extended witness's S-box
/// outputs (flattened, `SB_BITS·R` bits).
fn lowmc_run(p: &LowMcParams, key: &[u8], plaintext: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut state = xor(plaintext, &mat_vec(&p.key[0], key)); // whitening
    let mut sbox_outs = Vec::with_capacity(SB_BITS * R);
    for r in 0..R {
        // Partial S-box layer on the top 3m bits.
        let mut new_state = state.clone();
        for s in 0..M {
            let (a, b, c) = (state[3 * s], state[3 * s + 1], state[3 * s + 2]);
            let (d, e, f) = sbox(a, b, c);
            new_state[3 * s] = d;
            new_state[3 * s + 1] = e;
            new_state[3 * s + 2] = f;
            sbox_outs.extend_from_slice(&[d, e, f]);
        }
        // Linear layer, round constant, round key.
        state = xor(&mat_vec(&p.lin[r], &new_state), &p.rc[r]);
        state = xor(&state, &mat_vec(&p.key[r + 1], key));
    }
    (state, sbox_outs)
}

// ── Keys ──────────────────────────────────────────────────────────────────────

#[derive(Clone)]
pub struct PicnicPublicKey {
    pub plaintext: Vec<u8>,
    pub ciphertext: Vec<u8>,
}

#[derive(Clone)]
pub struct PicnicSecretKey {
    key: Vec<u8>,
}

pub type PicnicSignature = MpcithProof;

pub fn picnic_keygen() -> (PicnicPublicKey, PicnicSecretKey) {
    let params = expand_params();
    let mut kb = vec![0u8; N];
    random_bytes(&mut kb);
    let key: Vec<u8> = kb.iter().map(|b| b & 1).collect();
    let mut pb = vec![0u8; N];
    random_bytes(&mut pb);
    let plaintext: Vec<u8> = pb.iter().map(|b| b & 1).collect();
    let (ciphertext, _) = lowmc_run(&params, &key, &plaintext);
    (PicnicPublicKey { plaintext, ciphertext }, PicnicSecretKey { key })
}

/// Build the extended witness `key ‖ sbox_outputs` for the secret key.
fn extended_witness(params: &LowMcParams, sk: &PicnicSecretKey, pk: &PicnicPublicKey) -> Vec<u8> {
    let (_, sbox_outs) = lowmc_run(params, &sk.key, &pk.plaintext);
    let mut w = sk.key.clone();
    w.extend_from_slice(&sbox_outs);
    w
}

// ── The MPC relation: LowMC_k(p) = c with the AND gates checked ───────────────

struct PicnicRelation {
    params: LowMcParams,
    plaintext: Vec<u8>,
    ciphertext: Vec<u8>,
}

impl MpcRelation for PicnicRelation {
    fn witness_len(&self) -> usize {
        WITNESS
    }
    fn dot_len(&self) -> usize {
        GATES
    }
    fn eps_len(&self) -> usize {
        GATES
    }
    fn lin_len(&self) -> usize {
        N
    }
    fn party_compute(&self, wshare: &[u8], leader: bool, eps: &[u8]) -> PartyView {
        let key = &wshare[..N];
        let sbox_outs = &wshare[N..];
        // Simulate LowMC linearly on this party's share; whenever an
        // S-box is reached, take its output from the witness and record
        // the three AND-gate relations.
        let mut state = xor(
            &if leader { self.plaintext.clone() } else { vec![0u8; N] },
            &mat_vec(&self.params.key[0], key),
        );
        let mut u = vec![0u8; GATES]; // ε ∘ (left AND operand)
        let mut v = vec![0u8; GATES]; // right AND operand
        let mut t = 0u8; // Σ ε_g · target_g
        let mut gate = 0usize;
        for r in 0..R {
            let mut new_state = state.clone();
            for s in 0..M {
                let (a, b, c) = (state[3 * s], state[3 * s + 1], state[3 * s + 2]);
                let base = (r * M + s) * 3;
                let (d, e, f) = (sbox_outs[base], sbox_outs[base + 1], sbox_outs[base + 2]);
                // AND gates and their (linear) targets:
                //   b·c = d⊕a,  a·c = e⊕a⊕b,  a·b = f⊕a⊕b⊕c.
                let gates: [(u8, u8, u8); 3] =
                    [(b, c, d ^ a), (a, c, e ^ a ^ b), (a, b, f ^ a ^ b ^ c)];
                for (left, right, target) in gates {
                    u[gate] = gf_mul(eps[gate], left);
                    v[gate] = right;
                    t ^= gf_mul(eps[gate], target);
                    gate += 1;
                }
                new_state[3 * s] = d;
                new_state[3 * s + 1] = e;
                new_state[3 * s + 2] = f;
            }
            state = xor(&mat_vec(&self.params.lin[r], &new_state), &mat_vec(&self.params.key[r + 1], key));
            if leader {
                state = xor(&state, &self.params.rc[r]);
            }
        }
        // Linear check: the final state must equal the ciphertext.
        let lin = if leader { xor(&state, &self.ciphertext) } else { state };
        PartyView { u, v, t, lin }
    }
}

fn statement(pk: &PicnicPublicKey, msg: &[u8]) -> Vec<u8> {
    let mut input = b"Picnic-toy".to_vec();
    input.extend_from_slice(&pk.plaintext);
    input.extend_from_slice(&pk.ciphertext);
    input.extend_from_slice(msg);
    shake256(&input, 32)
}

pub fn picnic_sign(pk: &PicnicPublicKey, sk: &PicnicSecretKey, msg: &[u8]) -> PicnicSignature {
    let params = expand_params();
    let witness = extended_witness(&params, sk, pk);
    let rel =
        PicnicRelation { params, plaintext: pk.plaintext.clone(), ciphertext: pk.ciphertext.clone() };
    mpcith_prove(&rel, &witness, &statement(pk, msg))
}

pub fn picnic_verify(pk: &PicnicPublicKey, msg: &[u8], sig: &PicnicSignature) -> bool {
    let rel = PicnicRelation {
        params: expand_params(),
        plaintext: pk.plaintext.clone(),
        ciphertext: pk.ciphertext.clone(),
    };
    mpcith_verify(&rel, &statement(pk, msg), sig)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sbox_and_relations_hold() {
        // The three AND relations the proof checks must be identities of
        // the S-box.
        for v in 0u8..8 {
            let (a, b, c) = (v & 1, (v >> 1) & 1, (v >> 2) & 1);
            let (d, e, f) = sbox(a, b, c);
            assert_eq!(b & c, d ^ a);
            assert_eq!(a & c, e ^ a ^ b);
            assert_eq!(a & b, f ^ a ^ b ^ c);
        }
    }

    #[test]
    fn lowmc_is_deterministic() {
        let p = expand_params();
        let key = vec![1u8, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1];
        let pt = vec![0u8, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0];
        assert_eq!(lowmc_run(&p, &key, &pt).0, lowmc_run(&p, &key, &pt).0);
    }

    #[test]
    fn keygen_encrypts_correctly() {
        let (pk, sk) = picnic_keygen();
        let params = expand_params();
        assert_eq!(lowmc_run(&params, &sk.key, &pk.plaintext).0, pk.ciphertext);
    }

    #[test]
    fn sign_verify_roundtrip() {
        let (pk, sk) = picnic_keygen();
        let msg = b"LowMC in the head";
        let sig = picnic_sign(&pk, &sk, msg);
        assert!(picnic_verify(&pk, msg, &sig));
    }

    #[test]
    fn wrong_message_rejected() {
        let (pk, sk) = picnic_keygen();
        let sig = picnic_sign(&pk, &sk, b"one");
        assert!(!picnic_verify(&pk, b"two", &sig));
    }

    #[test]
    fn wrong_key_rejected() {
        let (pk_a, sk_a) = picnic_keygen();
        let (pk_b, _) = picnic_keygen();
        let sig = picnic_sign(&pk_a, &sk_a, b"m");
        assert!(!picnic_verify(&pk_b, b"m", &sig));
    }

    #[test]
    fn forged_witness_rejected() {
        // A wrong key gives a witness whose final state ≠ ciphertext, so
        // the linear check fails.
        let (pk, sk) = picnic_keygen();
        let params = expand_params();
        let mut bad = sk.clone();
        bad.key[0] ^= 1;
        let rel = PicnicRelation {
            params: expand_params(),
            plaintext: pk.plaintext.clone(),
            ciphertext: pk.ciphertext.clone(),
        };
        let (_, sbox_outs) = lowmc_run(&params, &bad.key, &pk.plaintext);
        let mut w = bad.key.clone();
        w.extend_from_slice(&sbox_outs);
        let msg = statement(&pk, b"m");
        let sig = mpcith_prove(&rel, &w, &msg);
        assert!(!mpcith_verify(&rel, &msg, &sig));
    }

    #[test]
    fn tampered_proof_rejected() {
        let (pk, sk) = picnic_keygen();
        let msg = b"tamper";
        let mut sig = picnic_sign(&pk, &sk, msg);
        sig.reps[0].hidden_lin[0] ^= 1;
        assert!(!picnic_verify(&pk, msg, &sig));
    }
}
