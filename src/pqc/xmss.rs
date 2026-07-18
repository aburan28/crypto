//! **XMSS** — the eXtended Merkle Signature Scheme
//! (Buchmann–Dahmen–Hülsing 2011; RFC 8391; NIST SP 800-208).
//!
//! # Relationship to LMS
//! XMSS is the other NIST-standardised stateful hash signature (see
//! `pqc::lms`).  Same skeleton — a Winternitz OTS at each leaf of a
//! Merkle tree — but with **randomised (masked) hashing** throughout,
//! which lets its security proof rely only on *second-preimage
//! resistance* (and PRF security) rather than collision resistance.
//! That is the technical reason XMSS predates and is often preferred to
//! LMS in the academic literature.
//!
//! # WOTS+ (the distinctive part)
//! Where LM-OTS iterates a plain hash, **WOTS+** interleaves a public
//! per-position **bitmask** and **key** into each chain step:
//!
//! ```text
//! c_{i}(x) = H( key_i ⊕-masked ( c_{i-1}(x) XOR mask_{i} ) ).
//! ```
//!
//! The masks/keys come from a public seed via a PRF, so every hash call
//! is on freshly randomised input — defeating multi-target attacks that
//! a fixed hash chain is vulnerable to.  The tree nodes are likewise
//! combined with per-level bitmasks (an "L-tree"/hash-tree with
//! randomised hashing).
//!
//! # This implementation
//! Toy params: SHA-256, Winternitz `w = 16` (4-bit digits), tree height
//! `h = 4` (16 one-time keys).  Faithful to WOTS+ masked chaining and
//! the bitmasked tree; the `PRF`/`H`/`H_msg` domain separation follows
//! RFC 8391 in spirit with simplified address encoding.  **Stateful**:
//! `XmssPrivateKey` tracks and consumes the leaf index.  Not
//! constant-time; see SECURITY.md.

use crate::hash::sha256::sha256;
use crate::utils::random::random_bytes;

pub const HASH: usize = 32;
/// Winternitz digit width (bits).
pub const W: usize = 4;
const CHAIN: usize = 1 << W;
const MSG_DIGITS: usize = (HASH * 8) / W; // 64
const CKSUM_DIGITS: usize = 3;
/// WOTS+ chains per leaf.
pub const LEN: usize = MSG_DIGITS + CKSUM_DIGITS; // 67
pub const H: usize = 4;
pub const LEAVES: usize = 1 << H;

// ── Keyed hash primitives (RFC 8391 style, simplified addressing) ────────────

/// PRF: derive a pseudorandom HASH-byte value from a key and address.
fn prf(key: &[u8], addr: &[u8]) -> [u8; HASH] {
    let mut buf = Vec::with_capacity(1 + key.len() + addr.len());
    buf.push(0x03); // PRF domain tag
    buf.extend_from_slice(key);
    buf.extend_from_slice(addr);
    sha256(&buf)
}

/// Keyed, bitmasked hash `H(key, m) = SHA256(tag ‖ key ‖ (m XOR mask))`
/// where `key` and `mask` are supplied by the caller.
fn hash_masked(key: &[u8; HASH], mask: &[u8; HASH], m: &[u8; HASH]) -> [u8; HASH] {
    let mut buf = Vec::with_capacity(1 + 2 * HASH);
    buf.push(0x00); // F/H domain tag
    buf.extend_from_slice(key);
    for i in 0..HASH {
        buf.push(m[i] ^ mask[i]);
    }
    sha256(&buf)
}

fn addr(kind: u8, a: usize, b: usize, c: usize) -> [u8; 16] {
    let mut out = [0u8; 16];
    out[0] = kind;
    out[1..5].copy_from_slice(&(a as u32).to_be_bytes());
    out[5..9].copy_from_slice(&(b as u32).to_be_bytes());
    out[9..13].copy_from_slice(&(c as u32).to_be_bytes());
    out
}

// ── WOTS+ ────────────────────────────────────────────────────────────────────

/// The masked WOTS+ chain: apply `steps` chaining steps to `x` starting
/// at position `start`, in chain `chain_idx` of leaf `leaf`.  Each step
/// draws a fresh key and bitmask from `pub_seed` via the PRF.
fn wots_chain(
    x: &[u8; HASH],
    pub_seed: &[u8],
    leaf: usize,
    chain_idx: usize,
    start: usize,
    steps: usize,
) -> [u8; HASH] {
    let mut cur = *x;
    for i in start..start + steps {
        let key = prf(pub_seed, &addr(0, leaf, chain_idx, 2 * i));
        let mask = prf(pub_seed, &addr(0, leaf, chain_idx, 2 * i + 1));
        cur = hash_masked(&key, &mask, &cur);
    }
    cur
}

fn digits(msg_hash: &[u8; HASH]) -> Vec<usize> {
    let mut d = Vec::with_capacity(LEN);
    for &byte in msg_hash.iter() {
        d.push((byte >> 4) as usize);
        d.push((byte & 0x0f) as usize);
    }
    let mut cksum: usize = d.iter().map(|&di| (CHAIN - 1) - di).sum();
    let mut cd = vec![0usize; CKSUM_DIGITS];
    for slot in cd.iter_mut().rev() {
        *slot = cksum & (CHAIN - 1);
        cksum >>= W;
    }
    d.extend_from_slice(&cd);
    d
}

/// WOTS+ secret chains for a leaf, derived from the secret seed.
fn wots_secret(sk_seed: &[u8], leaf: usize) -> Vec<[u8; HASH]> {
    (0..LEN).map(|i| prf(sk_seed, &addr(1, leaf, i, 0))).collect()
}

/// WOTS+ public key = the LEN chain tops.
fn wots_public(sk: &[[u8; HASH]], pub_seed: &[u8], leaf: usize) -> Vec<[u8; HASH]> {
    (0..LEN).map(|i| wots_chain(&sk[i], pub_seed, leaf, i, 0, CHAIN - 1)).collect()
}

fn wots_sign(
    sk: &[[u8; HASH]],
    pub_seed: &[u8],
    leaf: usize,
    msg_hash: &[u8; HASH],
) -> Vec<[u8; HASH]> {
    let d = digits(msg_hash);
    (0..LEN).map(|i| wots_chain(&sk[i], pub_seed, leaf, i, 0, d[i])).collect()
}

fn wots_public_from_sig(
    sig: &[[u8; HASH]],
    pub_seed: &[u8],
    leaf: usize,
    msg_hash: &[u8; HASH],
) -> Vec<[u8; HASH]> {
    let d = digits(msg_hash);
    (0..LEN)
        .map(|i| wots_chain(&sig[i], pub_seed, leaf, i, d[i], (CHAIN - 1) - d[i]))
        .collect()
}

// ── Bitmasked Merkle tree ────────────────────────────────────────────────────

/// Compress a WOTS+ public key (LEN hashes) into a single leaf via a
/// masked hash of the concatenation.
fn wots_pk_to_leaf(pk: &[[u8; HASH]], pub_seed: &[u8], leaf: usize) -> [u8; HASH] {
    let mut acc = [0u8; HASH];
    for (i, node) in pk.iter().enumerate() {
        let key = prf(pub_seed, &addr(2, leaf, i, 0));
        let mask = prf(pub_seed, &addr(2, leaf, i, 1));
        // Fold each chain top into the accumulator with a masked hash.
        let mut combined = [0u8; HASH];
        for j in 0..HASH {
            combined[j] = acc[j] ^ node[j];
        }
        acc = hash_masked(&key, &mask, &combined);
    }
    acc
}

/// Masked internal node: `H(key, (left‖right) masked)` — here folded as
/// two masked hashes so the node depends on both children with fresh
/// randomisation.
fn tree_node(left: &[u8; HASH], right: &[u8; HASH], pub_seed: &[u8], height: usize, index: usize) -> [u8; HASH] {
    let key = prf(pub_seed, &addr(3, height, index, 0));
    let mask_l = prf(pub_seed, &addr(3, height, index, 1));
    let mask_r = prf(pub_seed, &addr(3, height, index, 2));
    let mut buf = Vec::with_capacity(1 + HASH * 3);
    buf.push(0x01); // tree-hash tag
    buf.extend_from_slice(&key);
    for i in 0..HASH {
        buf.push(left[i] ^ mask_l[i]);
    }
    for i in 0..HASH {
        buf.push(right[i] ^ mask_r[i]);
    }
    sha256(&buf)
}

fn build_tree(sk_seed: &[u8], pub_seed: &[u8]) -> (Vec<[u8; HASH]>, [u8; HASH]) {
    let leaves: Vec<[u8; HASH]> = (0..LEAVES)
        .map(|l| {
            let pk = wots_public(&wots_secret(sk_seed, l), pub_seed, l);
            wots_pk_to_leaf(&pk, pub_seed, l)
        })
        .collect();
    let mut level = leaves.clone();
    let mut height = 0;
    while level.len() > 1 {
        level = level
            .chunks(2)
            .enumerate()
            .map(|(i, pair)| tree_node(&pair[0], &pair[1], pub_seed, height, i))
            .collect();
        height += 1;
    }
    (leaves, level[0])
}

fn auth_path(leaves: &[[u8; HASH]], pub_seed: &[u8], leaf: usize) -> Vec<[u8; HASH]> {
    let mut path = Vec::with_capacity(H);
    let mut level = leaves.to_vec();
    let mut idx = leaf;
    let mut height = 0;
    while level.len() > 1 {
        let sibling = if idx % 2 == 0 { idx + 1 } else { idx - 1 };
        path.push(level[sibling]);
        level = level
            .chunks(2)
            .enumerate()
            .map(|(i, pair)| tree_node(&pair[0], &pair[1], pub_seed, height, i))
            .collect();
        idx /= 2;
        height += 1;
    }
    path
}

fn root_from_path(
    leaf: usize,
    leaf_h: &[u8; HASH],
    path: &[[u8; HASH]],
    pub_seed: &[u8],
) -> [u8; HASH] {
    let mut cur = *leaf_h;
    let mut idx = leaf;
    for (height, sib) in path.iter().enumerate() {
        let parent_index = idx / 2;
        cur = if idx % 2 == 0 {
            tree_node(&cur, sib, pub_seed, height, parent_index)
        } else {
            tree_node(sib, &cur, pub_seed, height, parent_index)
        };
        idx /= 2;
    }
    cur
}

// ── Public API ────────────────────────────────────────────────────────────────

#[derive(Clone, Debug, PartialEq)]
pub struct XmssPublicKey {
    pub root: [u8; HASH],
    pub pub_seed: Vec<u8>,
}

#[derive(Clone)]
pub struct XmssPrivateKey {
    sk_seed: Vec<u8>,
    pub_seed: Vec<u8>,
    leaves: Vec<[u8; HASH]>,
    pub next_leaf: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct XmssSignature {
    pub leaf: usize,
    pub wots_sig: Vec<[u8; HASH]>,
    pub auth: Vec<[u8; HASH]>,
}

pub fn xmss_keygen() -> (XmssPublicKey, XmssPrivateKey) {
    let mut sk_seed = vec![0u8; 32];
    let mut pub_seed = vec![0u8; 32];
    random_bytes(&mut sk_seed);
    random_bytes(&mut pub_seed);
    let (leaves, root) = build_tree(&sk_seed, &pub_seed);
    (
        XmssPublicKey { root, pub_seed: pub_seed.clone() },
        XmssPrivateKey { sk_seed, pub_seed, leaves, next_leaf: 0 },
    )
}

pub fn xmss_sign(sk: &mut XmssPrivateKey, msg: &[u8]) -> Option<XmssSignature> {
    if sk.next_leaf >= LEAVES {
        return None;
    }
    let leaf = sk.next_leaf;
    sk.next_leaf += 1; // one-time: advance before returning
    let msg_hash = sha256(msg);
    let wsk = wots_secret(&sk.sk_seed, leaf);
    let wots_sig = wots_sign(&wsk, &sk.pub_seed, leaf, &msg_hash);
    let auth = auth_path(&sk.leaves, &sk.pub_seed, leaf);
    Some(XmssSignature { leaf, wots_sig, auth })
}

pub fn xmss_verify(pk: &XmssPublicKey, msg: &[u8], sig: &XmssSignature) -> bool {
    if sig.leaf >= LEAVES || sig.wots_sig.len() != LEN || sig.auth.len() != H {
        return false;
    }
    let msg_hash = sha256(msg);
    let wots_pk = wots_public_from_sig(&sig.wots_sig, &pk.pub_seed, sig.leaf, &msg_hash);
    let leaf_h = wots_pk_to_leaf(&wots_pk, &pk.pub_seed, sig.leaf);
    root_from_path(sig.leaf, &leaf_h, &sig.auth, &pk.pub_seed) == pk.root
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wots_chain_uses_fresh_masks() {
        // Two positions in a chain must use different masks, so the same
        // input hashes differently at different steps.
        let x = [7u8; HASH];
        let seed = b"pubseed-000000000000000000000000";
        let a = wots_chain(&x, seed, 0, 0, 0, 1);
        let b = wots_chain(&x, seed, 0, 0, 1, 1);
        assert_ne!(a, b);
    }

    #[test]
    fn wots_sign_verify() {
        let sk_seed = b"skseed-0000000000000000000000000";
        let pub_seed = b"pubseed-000000000000000000000000";
        let sk = wots_secret(sk_seed, 0);
        let pk = wots_public(&sk, pub_seed, 0);
        let mh = sha256(b"wots message");
        let sig = wots_sign(&sk, pub_seed, 0, &mh);
        assert_eq!(wots_public_from_sig(&sig, pub_seed, 0, &mh), pk);
    }

    #[test]
    fn xmss_sign_verify() {
        let (pk, mut sk) = xmss_keygen();
        let sig = xmss_sign(&mut sk, b"hello XMSS").unwrap();
        assert!(xmss_verify(&pk, b"hello XMSS", &sig));
        assert!(!xmss_verify(&pk, b"tampered", &sig));
    }

    #[test]
    fn state_advances_and_exhausts() {
        let (pk, mut sk) = xmss_keygen();
        let mut seen = std::collections::HashSet::new();
        for i in 0..LEAVES {
            let sig = xmss_sign(&mut sk, format!("m{i}").as_bytes()).unwrap();
            assert!(seen.insert(sig.leaf));
            assert!(xmss_verify(&pk, format!("m{i}").as_bytes(), &sig));
        }
        assert!(xmss_sign(&mut sk, b"overflow").is_none());
    }

    #[test]
    fn forged_message_rejected() {
        let (pk, mut sk) = xmss_keygen();
        let sig = xmss_sign(&mut sk, b"real").unwrap();
        assert!(!xmss_verify(&pk, b"forged", &sig));
    }

    #[test]
    fn wrong_pub_seed_rejected() {
        // The bitmasks depend on pub_seed; a different one breaks verify.
        let (mut pk, mut sk) = xmss_keygen();
        let sig = xmss_sign(&mut sk, b"m").unwrap();
        pk.pub_seed[0] ^= 1;
        assert!(!xmss_verify(&pk, b"m", &sig));
    }
}
