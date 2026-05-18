//! # Challenge 52 — Iterated Hash Function Multicollisions
//!
//! Antoine Joux (CRYPTO 2004) observed that Merkle–Damgård hashes
//! have a structural weakness: finding *one* collision in the
//! compression function lets you bootstrap to `2^n` collisions in
//! the full hash with only `n · 2^(b/2)` work — not the `2^(b/2·n)`
//! one would naïvely expect.
//!
//! The trick: a collision between two message blocks at every
//! Merkle–Damgård level acts as an independent fork point.  After
//! `n` levels you have `2^n` message *strings* that all hash to the
//! same final value.
//!
//! ## Cascading is not security
//!
//! Naïve advice: "use h(x) = f(x) || g(x)" — combine two
//! independent hash functions to add their security widths.  Joux
//! refutes this for iterated hashes: generate a `2^(b₂/2)`-way
//! multicollision in `f` (cheap), then do a birthday search on `g`
//! within that set.  Total cost ≈ `(b₂/2) · 2^(b_f/2) + 2^(b_g/2)`
//! — dominated by the larger hash's birthday bound, not the sum.
//!
//! ## This module's weak hash
//!
//! Use AES-128 as the compression-function block cipher.  State is
//! truncated to `b` bits.  The compression function is
//! `C(state, block) = truncate(AES_K(block) ^ block, b/8 bytes)`
//! where the *block* is the input message-block and the *key* is
//! the running state, padded to 16 bytes.  `b` is configurable —
//! the demo uses 16-bit `f` and 24-bit `g` so birthday work fits
//! in seconds.

use crate::cryptopals::Report;
use crate::symmetric::aes::{encrypt_block, AesKey};
use std::collections::HashMap;

/// One Davies–Meyer–style step over a 16-byte message block.  The
/// state is padded to 16 bytes by repeating, used as the AES key.
/// The output is truncated to `b_bytes` bytes.
pub fn compress(state: &[u8], block: &[u8; 16], b_bytes: usize) -> Vec<u8> {
    let mut key_bytes = [0u8; 16];
    for (i, slot) in key_bytes.iter_mut().enumerate() {
        *slot = state[i % state.len()];
    }
    let key = AesKey::new(&key_bytes).unwrap();
    let mut out = encrypt_block(block, &key);
    // XOR Davies–Meyer feedforward.
    for i in 0..16 {
        out[i] ^= block[i];
    }
    out[..b_bytes].to_vec()
}

/// Hash a multi-block message under our toy compression function.
/// `b_bytes` controls digest width (and therefore birthday cost).
/// IV is hard-coded to `0`; the message must be a multiple of 16
/// bytes (callers prepad as needed — multicollision construction
/// stays inside block boundaries so this is fine).
pub fn weak_hash(message: &[u8], iv: &[u8], b_bytes: usize) -> Vec<u8> {
    assert!(message.len() % 16 == 0);
    let mut state = iv.to_vec();
    for chunk in message.chunks_exact(16) {
        let mut block = [0u8; 16];
        block.copy_from_slice(chunk);
        state = compress(&state, &block, b_bytes);
    }
    state
}

/// Find a single-block collision under `compress(state, ·, b_bytes)`.
/// Returns `(blockA, blockB, next_state)`.  Uses a hash-map birthday
/// search; expected ~`2^(b_bytes·4)` queries.
pub fn find_block_collision(
    state: &[u8],
    b_bytes: usize,
) -> ([u8; 16], [u8; 16], Vec<u8>) {
    let mut seen: HashMap<Vec<u8>, [u8; 16]> = HashMap::new();
    let mut block = [0u8; 16];
    let mut counter: u64 = 0;
    loop {
        block[..8].copy_from_slice(&counter.to_le_bytes());
        let out = compress(state, &block, b_bytes);
        if let Some(prev) = seen.get(&out) {
            if *prev != block {
                return (*prev, block, out);
            }
        } else {
            seen.insert(out, block);
        }
        counter += 1;
        // Safety valve so a misbehaving cipher can't trap us forever.
        if counter > (1u64 << 30) {
            panic!("birthday search exceeded 2^30 queries");
        }
    }
}

/// Joux multicollision: produce `2^n` messages, each `n` blocks
/// long, all hashing to the same `b_bytes`-byte digest.  Returns
/// `(per-level [a,b] pairs, final state)`.  Caller expands the
/// pairs into the `2^n` messages.
pub fn joux_multicollision(
    n: usize,
    iv: &[u8],
    b_bytes: usize,
) -> (Vec<([u8; 16], [u8; 16])>, Vec<u8>) {
    let mut state = iv.to_vec();
    let mut pairs = Vec::with_capacity(n);
    for _ in 0..n {
        let (a, b, next) = find_block_collision(&state, b_bytes);
        pairs.push((a, b));
        state = next;
    }
    (pairs, state)
}

/// Expand a multicollision pair-list into all `2^n` messages.
pub fn expand_multicollision(pairs: &[([u8; 16], [u8; 16])]) -> Vec<Vec<u8>> {
    let mut messages = vec![Vec::new()];
    for (a, b) in pairs {
        let mut next = Vec::with_capacity(2 * messages.len());
        for m in &messages {
            let mut ma = m.clone();
            ma.extend_from_slice(a);
            let mut mb = m.clone();
            mb.extend_from_slice(b);
            next.push(ma);
            next.push(mb);
        }
        messages = next;
    }
    messages
}

/// Concrete cascade attack: find two messages that collide under
/// `h(x) = weak_hash_f(x) || weak_hash_g(x)`.  Strategy:
///
/// 1. Generate a `2^(b_g · 4)` multicollision in the `f`-hash.
/// 2. Probe each multicollision member under `g`; expect a
///    collision after ~`2^(b_g · 4)` probes (birthday in `g`).
pub fn cascade_collision(
    f_bytes: usize,
    g_bytes: usize,
) -> Option<(Vec<u8>, Vec<u8>)> {
    // Joux multicollision of size 2^n where n = b_g · 8 / 2.
    let n = g_bytes * 4;
    let iv_f = vec![0u8; f_bytes];
    let iv_g = vec![1u8; g_bytes];
    let (pairs, _) = joux_multicollision(n, &iv_f, f_bytes);
    let candidates = expand_multicollision(&pairs);
    let mut seen: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
    for m in candidates {
        let g = weak_hash(&m, &iv_g, g_bytes);
        if let Some(prev) = seen.get(&g) {
            if *prev != m {
                return Some((prev.clone(), m));
            }
        } else {
            seen.insert(g, m);
        }
    }
    None
}

pub fn run() -> Report {
    let mut r = Report::new(52, "Iterated Hash Function Multicollisions");
    // 16-bit toy hash so the birthday search runs in milliseconds.
    let b_bytes = 2;
    let iv = vec![0u8; b_bytes];

    // ── Part 1: build a 2^4 = 16-way multicollision in `f`. ──
    let n = 4;
    let (pairs, final_state) = joux_multicollision(n, &iv, b_bytes);
    let messages = expand_multicollision(&pairs);
    r.line(format!(
        "Built {}-way Joux multicollision (n={}, b={} bits)",
        messages.len(),
        n,
        b_bytes * 8
    ));
    r.line(format!("Final state    : {}", hex::encode(&final_state)));
    let h0 = weak_hash(&messages[0], &iv, b_bytes);
    let all_same = messages.iter().all(|m| weak_hash(m, &iv, b_bytes) == h0);
    r.line(format!("All messages share hash {}: {}", hex::encode(&h0), all_same));
    assert!(all_same);

    // ── Part 2: cascade attack.  f = 16-bit, g = 24-bit. ──
    let f_bytes = 2;
    let g_bytes = 3;
    r.line("");
    r.line(format!(
        "Cascading h(x) = f(x) || g(x) with f={} bits, g={} bits …",
        f_bytes * 8,
        g_bytes * 8
    ));
    let (m1, m2) = cascade_collision(f_bytes, g_bytes).expect("cascade attack must succeed");
    let iv_f = vec![0u8; f_bytes];
    let iv_g = vec![1u8; g_bytes];
    let cascaded = |m: &[u8]| {
        let mut h = weak_hash(m, &iv_f, f_bytes);
        h.extend_from_slice(&weak_hash(m, &iv_g, g_bytes));
        h
    };
    let h1 = cascaded(&m1);
    let h2 = cascaded(&m2);
    r.line(format!("m1 length      : {} bytes", m1.len()));
    r.line(format!("m2 length      : {} bytes", m2.len()));
    r.line(format!("h(m1) = {}", hex::encode(&h1)));
    r.line(format!("h(m2) = {}", hex::encode(&h2)));
    r.line(format!("colliding      : {}", h1 == h2 && m1 != m2));
    assert_eq!(h1, h2);
    assert_ne!(m1, m2);
    r.succeed()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn block_collision_works_in_16_bits() {
        let (a, b, _) = find_block_collision(&[0u8; 2], 2);
        assert_ne!(a, b);
    }

    #[test]
    fn joux_pairs_collide() {
        let iv = vec![0u8; 2];
        let (pairs, _) = joux_multicollision(3, &iv, 2);
        let msgs = expand_multicollision(&pairs);
        let h0 = weak_hash(&msgs[0], &iv, 2);
        for m in &msgs[1..] {
            assert_eq!(weak_hash(m, &iv, 2), h0);
        }
    }

    #[test]
    fn cascade_attack_finds_collision() {
        let (m1, m2) = cascade_collision(2, 3).unwrap();
        assert_ne!(m1, m2);
        let iv_f = vec![0u8; 2];
        let iv_g = vec![1u8; 3];
        assert_eq!(weak_hash(&m1, &iv_f, 2), weak_hash(&m2, &iv_f, 2));
        assert_eq!(weak_hash(&m1, &iv_g, 3), weak_hash(&m2, &iv_g, 3));
    }
}
