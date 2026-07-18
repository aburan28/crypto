//! **LMS / HSS** — Leighton–Micali (hierarchical) hash-based signatures
//! (RFC 8554; NIST SP 800-208).
//!
//! # Why stateful hash signatures
//! LMS and XMSS are the *stateful* hash-based signatures — the most
//! conservative post-quantum signatures that exist, relying only on the
//! second-preimage resistance of a hash function.  Unlike SLH-DSA
//! (stateless SPHINCS+), they require the signer to **never reuse a
//! one-time key**, tracked as mutable state.  In exchange they are far
//! smaller and faster than SLH-DSA.  Standardised by NIST in SP 800-208
//! for firmware/boot signing, where the signer is a controlled
//! environment that can maintain state.
//!
//! # Construction
//! Two layers:
//!
//! 1. **LM-OTS** (a Winternitz one-time signature): to sign a hash `h`,
//!    split it (plus a checksum) into `w`-bit digits `d₁…d_p`.  The
//!    secret key is `p` random strings `x₁…x_p`; the public key hashes
//!    `yᵢ = H^{2^w−1}(xᵢ)`.  The signature reveals `H^{dᵢ}(xᵢ)`, and the
//!    verifier finishes the chain `H^{2^w−1−dᵢ}(sigᵢ) = yᵢ`.  The
//!    checksum makes it infeasible to forge a *different* digest (any
//!    increase in one digit forces a decrease in the checksum, which the
//!    one-way chains forbid).
//! 2. **LMS** (a Merkle tree): `2^h` LM-OTS key pairs form the leaves;
//!    each leaf is `H(OTS public key)`; internal nodes hash their
//!    children; the root is the LMS public key.  A signature is
//!    `(index q, LM-OTS signature, authentication path)`.
//!
//! **HSS** stacks LMS trees (a top tree signs the roots of lower trees)
//! to get a huge key space without a huge single tree; this module
//! includes a 2-level HSS on top of LMS.
//!
//! # This implementation
//! Toy params: SHA-256, Winternitz `w = 4` (so `2^w = 16`, digits are
//! nibbles), tree height `h = 4` (16 one-time keys per tree).  Faithful
//! to RFC 8554's chaining + checksum + Merkle structure; domain-
//! separation tags are simplified.  **Stateful**: `LmsPrivateKey`
//! tracks the next leaf and refuses to reuse one.  Not constant-time;
//! see SECURITY.md.

use crate::hash::sha256::sha256;
use crate::utils::random::random_bytes;

/// Hash output / chain element length.
pub const HASH: usize = 32;
/// Winternitz parameter: digits are `w` bits.
pub const W: usize = 4;
/// Chain length `2^w`.
const CHAIN: usize = 1 << W;
/// Number of message digits (256-bit hash in 4-bit digits).
const MSG_DIGITS: usize = (HASH * 8) / W; // 64
/// Checksum digits (covers the max checksum value in base 2^w).
const CKSUM_DIGITS: usize = 3;
/// Total LM-OTS chains.
pub const P: usize = MSG_DIGITS + CKSUM_DIGITS; // 67
/// Merkle tree height.
pub const H: usize = 4;
/// Leaves per LMS tree.
pub const LEAVES: usize = 1 << H;

// ── LM-OTS (Winternitz one-time signature) ───────────────────────────────────

/// Apply the hash chain `H^n` starting from `x`, with domain separation
/// by the chain index and iteration (a simplified LM-OTS tag).
fn chain(x: &[u8; HASH], chain_idx: usize, start: usize, steps: usize) -> [u8; HASH] {
    let mut cur = *x;
    for i in start..start + steps {
        let mut buf = Vec::with_capacity(HASH + 8);
        buf.extend_from_slice(b"LMOTS");
        buf.extend_from_slice(&(chain_idx as u32).to_be_bytes());
        buf.extend_from_slice(&(i as u32).to_be_bytes());
        buf.extend_from_slice(&cur);
        cur = sha256(&buf);
    }
    cur
}

/// Split a 32-byte digest (plus checksum) into `P` base-`2^w` digits.
fn digits(msg_hash: &[u8; HASH]) -> Vec<usize> {
    let mut d = Vec::with_capacity(P);
    // Message digits: high nibble then low nibble of each byte (w = 4).
    for &byte in msg_hash.iter() {
        d.push((byte >> 4) as usize);
        d.push((byte & 0x0f) as usize);
    }
    // Checksum = Σ (2^w − 1 − dᵢ), spread over CKSUM_DIGITS base-2^w digits.
    let mut cksum: usize = d.iter().map(|&di| (CHAIN - 1) - di).sum();
    let mut cd = vec![0usize; CKSUM_DIGITS];
    for slot in cd.iter_mut().rev() {
        *slot = cksum & (CHAIN - 1);
        cksum >>= W;
    }
    d.extend_from_slice(&cd);
    d
}

/// One LM-OTS key pair, derived deterministically from a per-leaf seed.
fn lmots_secret(seed: &[u8], leaf: usize) -> Vec<[u8; HASH]> {
    (0..P)
        .map(|i| {
            let mut buf = Vec::new();
            buf.extend_from_slice(seed);
            buf.extend_from_slice(&(leaf as u32).to_be_bytes());
            buf.extend_from_slice(&(i as u32).to_be_bytes());
            sha256(&buf)
        })
        .collect()
}

/// LM-OTS public key = H(all chain tops concatenated).
fn lmots_public(sk: &[[u8; HASH]]) -> [u8; HASH] {
    let mut buf = Vec::with_capacity(P * HASH);
    for (i, x) in sk.iter().enumerate() {
        buf.extend_from_slice(&chain(x, i, 0, CHAIN - 1));
    }
    sha256(&buf)
}

fn lmots_sign(sk: &[[u8; HASH]], msg_hash: &[u8; HASH]) -> Vec<[u8; HASH]> {
    let d = digits(msg_hash);
    (0..P).map(|i| chain(&sk[i], i, 0, d[i])).collect()
}

fn lmots_public_from_sig(sig: &[[u8; HASH]], msg_hash: &[u8; HASH]) -> [u8; HASH] {
    let d = digits(msg_hash);
    let mut buf = Vec::with_capacity(P * HASH);
    for i in 0..P {
        // Finish the chain: from d[i] up to 2^w − 1.
        buf.extend_from_slice(&chain(&sig[i], i, d[i], (CHAIN - 1) - d[i]));
    }
    sha256(&buf)
}

// ── LMS (Merkle tree over LM-OTS leaves) ─────────────────────────────────────

fn leaf_hash(leaf: usize, ots_pub: &[u8; HASH]) -> [u8; HASH] {
    let mut buf = Vec::new();
    buf.extend_from_slice(b"LMSLEAF");
    buf.extend_from_slice(&(leaf as u32).to_be_bytes());
    buf.extend_from_slice(ots_pub);
    sha256(&buf)
}

fn node_hash(left: &[u8; HASH], right: &[u8; HASH]) -> [u8; HASH] {
    let mut buf = Vec::new();
    buf.extend_from_slice(b"LMSNODE");
    buf.extend_from_slice(left);
    buf.extend_from_slice(right);
    sha256(&buf)
}

/// Build the full Merkle tree; returns all leaf hashes and the root.
fn build_tree(seed: &[u8]) -> (Vec<[u8; HASH]>, [u8; HASH]) {
    let leaves: Vec<[u8; HASH]> = (0..LEAVES)
        .map(|l| leaf_hash(l, &lmots_public(&lmots_secret(seed, l))))
        .collect();
    let mut level = leaves.clone();
    while level.len() > 1 {
        level = level.chunks(2).map(|pair| node_hash(&pair[0], &pair[1])).collect();
    }
    (leaves, level[0])
}

/// Authentication path (sibling hashes) for `leaf`.
fn auth_path(leaves: &[[u8; HASH]], leaf: usize) -> Vec<[u8; HASH]> {
    let mut path = Vec::with_capacity(H);
    let mut level = leaves.to_vec();
    let mut idx = leaf;
    while level.len() > 1 {
        let sibling = if idx % 2 == 0 { idx + 1 } else { idx - 1 };
        path.push(level[sibling]);
        level = level.chunks(2).map(|pair| node_hash(&pair[0], &pair[1])).collect();
        idx /= 2;
    }
    path
}

fn root_from_path(leaf: usize, leaf_h: &[u8; HASH], path: &[[u8; HASH]]) -> [u8; HASH] {
    let mut cur = *leaf_h;
    let mut idx = leaf;
    for sib in path {
        cur = if idx % 2 == 0 { node_hash(&cur, sib) } else { node_hash(sib, &cur) };
        idx /= 2;
    }
    cur
}

/// The LMS public key: the Merkle root.
#[derive(Clone, Debug, PartialEq)]
pub struct LmsPublicKey {
    pub root: [u8; HASH],
}

/// The LMS private key: the seed, the precomputed leaves, and the
/// **mutable next-leaf counter** — the state that must never rewind.
#[derive(Clone)]
pub struct LmsPrivateKey {
    seed: Vec<u8>,
    leaves: Vec<[u8; HASH]>,
    pub next_leaf: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct LmsSignature {
    pub leaf: usize,
    pub ots_sig: Vec<[u8; HASH]>,
    pub auth: Vec<[u8; HASH]>,
}

pub fn lms_keygen() -> (LmsPublicKey, LmsPrivateKey) {
    let mut seed = vec![0u8; 32];
    random_bytes(&mut seed);
    let (leaves, root) = build_tree(&seed);
    (LmsPublicKey { root }, LmsPrivateKey { seed, leaves, next_leaf: 0 })
}

/// Sign, consuming the next one-time leaf.  Returns `None` when the tree
/// is exhausted (all `2^h` one-time keys used) — the signer must never
/// reuse a leaf.
pub fn lms_sign(sk: &mut LmsPrivateKey, msg: &[u8]) -> Option<LmsSignature> {
    if sk.next_leaf >= LEAVES {
        return None;
    }
    let leaf = sk.next_leaf;
    sk.next_leaf += 1; // advance state before returning: one-time use
    let msg_hash = sha256(msg);
    let ots_sk = lmots_secret(&sk.seed, leaf);
    let ots_sig = lmots_sign(&ots_sk, &msg_hash);
    let auth = auth_path(&sk.leaves, leaf);
    Some(LmsSignature { leaf, ots_sig, auth })
}

pub fn lms_verify(pk: &LmsPublicKey, msg: &[u8], sig: &LmsSignature) -> bool {
    if sig.leaf >= LEAVES || sig.ots_sig.len() != P || sig.auth.len() != H {
        return false;
    }
    let msg_hash = sha256(msg);
    let ots_pub = lmots_public_from_sig(&sig.ots_sig, &msg_hash);
    let leaf_h = leaf_hash(sig.leaf, &ots_pub);
    root_from_path(sig.leaf, &leaf_h, &sig.auth) == pk.root
}

// ── HSS: a 2-level hierarchy (top tree signs the bottom tree's root) ──────────

#[derive(Clone)]
pub struct HssPrivateKey {
    /// The top tree.  In this two-level toy the bottom certificate is
    /// precomputed at keygen, so `top` is retained only as part of the
    /// key material (a deeper HSS would re-sign new bottom roots with it).
    #[allow(dead_code)]
    top: LmsPrivateKey,
    bottom: LmsPrivateKey,
    bottom_pub: LmsPublicKey,
    /// The top-tree signature over the bottom tree's public key, made
    /// once and reused for every message signed under the bottom tree.
    bottom_cert: LmsSignature,
}

#[derive(Clone, Debug, PartialEq)]
pub struct HssSignature {
    pub bottom_pub: LmsPublicKey,
    pub bottom_cert: LmsSignature,
    pub msg_sig: LmsSignature,
}

pub fn hss_keygen() -> (LmsPublicKey, HssPrivateKey) {
    let (top_pub, mut top) = lms_keygen();
    let (bottom_pub, bottom) = lms_keygen();
    // Top tree certifies the bottom tree's root (consumes top leaf 0).
    let bottom_cert = lms_sign(&mut top, &bottom_pub.root).unwrap();
    (top_pub, HssPrivateKey { top, bottom, bottom_pub, bottom_cert })
}

pub fn hss_sign(sk: &mut HssPrivateKey, msg: &[u8]) -> Option<HssSignature> {
    let msg_sig = lms_sign(&mut sk.bottom, msg)?;
    Some(HssSignature {
        bottom_pub: sk.bottom_pub.clone(),
        bottom_cert: sk.bottom_cert.clone(),
        msg_sig,
    })
}

pub fn hss_verify(top_pub: &LmsPublicKey, msg: &[u8], sig: &HssSignature) -> bool {
    // 1. Top tree must have certified this bottom public key.
    if !lms_verify(top_pub, &sig.bottom_pub.root, &sig.bottom_cert) {
        return false;
    }
    // 2. Bottom tree must have signed the message.
    lms_verify(&sig.bottom_pub, msg, &sig.msg_sig)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ots_sign_verify() {
        let seed = b"test-seed-0000000000000000000000";
        let sk = lmots_secret(seed, 0);
        let pk = lmots_public(&sk);
        let mh = sha256(b"one-time message");
        let sig = lmots_sign(&sk, &mh);
        assert_eq!(lmots_public_from_sig(&sig, &mh), pk);
    }

    #[test]
    fn ots_forgery_on_different_message_fails() {
        // The checksum defeats forging a different digest from one sig.
        let seed = b"test-seed-1111111111111111111111";
        let sk = lmots_secret(seed, 0);
        let pk = lmots_public(&sk);
        let sig = lmots_sign(&sk, &sha256(b"real"));
        let forged_pub = lmots_public_from_sig(&sig, &sha256(b"forged"));
        assert_ne!(forged_pub, pk);
    }

    #[test]
    fn lms_sign_verify() {
        let (pk, mut sk) = lms_keygen();
        let sig = lms_sign(&mut sk, b"hello LMS").unwrap();
        assert!(lms_verify(&pk, b"hello LMS", &sig));
        assert!(!lms_verify(&pk, b"tampered", &sig));
    }

    #[test]
    fn state_advances_and_exhausts() {
        let (pk, mut sk) = lms_keygen();
        let mut seen = std::collections::HashSet::new();
        for i in 0..LEAVES {
            let sig = lms_sign(&mut sk, format!("msg {i}").as_bytes()).unwrap();
            assert!(seen.insert(sig.leaf), "leaf {} reused!", sig.leaf);
            assert!(lms_verify(&pk, format!("msg {i}").as_bytes(), &sig));
        }
        // Tree exhausted: no more one-time keys.
        assert!(lms_sign(&mut sk, b"overflow").is_none());
    }

    #[test]
    fn different_leaves_authenticate_to_same_root() {
        let (pk, mut sk) = lms_keygen();
        let s0 = lms_sign(&mut sk, b"a").unwrap();
        let s1 = lms_sign(&mut sk, b"b").unwrap();
        assert_ne!(s0.leaf, s1.leaf);
        assert!(lms_verify(&pk, b"a", &s0));
        assert!(lms_verify(&pk, b"b", &s1));
    }

    #[test]
    fn hss_two_level_sign_verify() {
        let (top_pub, mut sk) = hss_keygen();
        let sig = hss_sign(&mut sk, b"hierarchical").unwrap();
        assert!(hss_verify(&top_pub, b"hierarchical", &sig));
        assert!(!hss_verify(&top_pub, b"wrong", &sig));
    }

    #[test]
    fn hss_rejects_uncertified_bottom_key() {
        let (top_pub, mut sk) = hss_keygen();
        let mut sig = hss_sign(&mut sk, b"m").unwrap();
        // Swap in a foreign bottom tree the top never certified.
        let (other_pub, mut other) = lms_keygen();
        sig.bottom_pub = other_pub;
        sig.msg_sig = lms_sign(&mut other, b"m").unwrap();
        assert!(!hss_verify(&top_pub, b"m", &sig));
    }
}
