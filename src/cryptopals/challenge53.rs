//! # Challenge 53 — Kelsey & Schneier "Expandable Messages"
//!
//! Long-message second-pre-image attack on Merkle-Damgård.  For a
//! `b`-bit hash you'd hope second pre-image costs `2^b`, but for a
//! message containing many blocks the attacker only needs to
//! collide with *any* intermediate state — not the final tag —
//! reducing the search to roughly `2^(b/2) + 2^(b - k)` for a
//! target with `2^k` blocks.
//!
//! ## Expandable message construction
//!
//! An `(k, k + 2^k - 1)`-expandable message is a set of length-pairs
//! `{(1, 2^(k-1)+1), (1, 2^(k-2)+1), …, (1, 2)}` of single-block /
//! many-block collisions chained together.  Each level lets the
//! attacker spend either 1 block or `2^i+1` blocks while landing on
//! the same intermediate state.  Composing all levels gives a
//! message of *any* length in `[k, k + 2^k - 1]` blocks that hashes
//! to the same state.
//!
//! ## Second pre-image
//!
//! Given a long target message `M` with `2^k` blocks:
//!
//! 1. Hash `M`; record every intermediate state and its block
//!    index in a map `state → index`.
//! 2. Build the expandable message (final state `s`).
//! 3. Find a single "bridge block" `B` such that
//!    `compress(s, B) = state-at-index-j` for some `j > k`.
//! 4. Use the expandable message to produce a prefix of length
//!    `j - 1` blocks, then `B`, then `M[j..]`.  Total length = `M`'s
//!    length, so padding matches, so the full hash matches.

use crate::cryptopals::Report;
use crate::cryptopals::challenge52::{compress, weak_hash};
use std::collections::HashMap;

/// One level of the expandable-message construction.  Given a
/// starting state `s`, find a single-block message and a
/// `(2^i + 1)`-block message that produce the same next state.
///
/// Algorithm:
/// 1. Run `2^i` arbitrary "dummy" blocks from `s` to produce some
///    state `s'`.
/// 2. Birthday: generate random single blocks from `s` (call those
///    states `Aₘ`) and random single blocks from `s'` (states
///    `Bₙ`).  Look for `Aₘ == Bₙ` — the collision.
///
/// Returns `(short_block, long_blocks, next_state)`.
pub fn expandable_level(
    state: &[u8],
    i: usize,
    b_bytes: usize,
) -> ([u8; 16], Vec<[u8; 16]>, Vec<u8>) {
    // Step 1: choose 2^i dummy blocks (we use a deterministic ramp
    // so the result is reproducible) and run them.
    let n_dummies = 1usize << i;
    let mut dummies: Vec<[u8; 16]> = Vec::with_capacity(n_dummies);
    let mut s_prime = state.to_vec();
    for d in 0..n_dummies {
        let mut block = [0u8; 16];
        block[..8].copy_from_slice(&(d as u64).to_le_bytes());
        block[15] = 0xAA;
        s_prime = compress(&s_prime, &block, b_bytes);
        dummies.push(block);
    }
    // Step 2: birthday between (single-block from state) and
    // (single-block from s_prime).
    let mut short_map: HashMap<Vec<u8>, [u8; 16]> = HashMap::new();
    let mut counter: u64 = 0;
    loop {
        let mut sb = [0u8; 16];
        sb[..8].copy_from_slice(&counter.to_le_bytes());
        sb[15] = 0x55;
        let s_short = compress(state, &sb, b_bytes);
        short_map.insert(s_short, sb);
        let mut lb = [0u8; 16];
        lb[..8].copy_from_slice(&counter.to_le_bytes());
        lb[15] = 0x77;
        let s_long = compress(&s_prime, &lb, b_bytes);
        if let Some(short_block) = short_map.get(&s_long).cloned() {
            // The "long path" = dummies || lb (total = 2^i + 1 blocks).
            let mut long_blocks = dummies;
            long_blocks.push(lb);
            return (short_block, long_blocks, s_long);
        }
        counter += 1;
        if counter > (1u64 << 24) {
            panic!("expandable_level: birthday search blew its budget");
        }
    }
}

/// Build a `(k, k + 2^k - 1)`-expandable message from initial state
/// `iv`.  Returns `(levels, final_state)`.  Each `level[i]` is
/// `(short_block, long_blocks_for_level)` with
/// `long_blocks.len() = 2^(k-1-i) + 1`.
pub fn build_expandable(
    k: usize,
    iv: &[u8],
    b_bytes: usize,
) -> (Vec<([u8; 16], Vec<[u8; 16]>)>, Vec<u8>) {
    let mut state = iv.to_vec();
    let mut levels = Vec::with_capacity(k);
    for i in (0..k).rev() {
        let (sb, lb, next) = expandable_level(&state, i, b_bytes);
        state = next;
        levels.push((sb, lb));
    }
    (levels, state)
}

/// Given an expandable message and a desired total length `n`
/// (in blocks, with `k ≤ n ≤ k + 2^k - 1`), produce the actual
/// message.  Greedily picks each level: if remaining target length
/// minus remaining short-only cost is at least the long length,
/// choose long; otherwise short.
pub fn assemble_expandable(
    levels: &[([u8; 16], Vec<[u8; 16]>)],
    target_blocks: usize,
) -> Vec<u8> {
    let k = levels.len();
    assert!(target_blocks >= k);
    let max_extra = (1usize << k) - 1;
    let mut extra = target_blocks - k;
    assert!(extra <= max_extra);
    let mut out = Vec::new();
    for (i, (short, long)) in levels.iter().enumerate() {
        // Level i lets us add either 1 block (short) or 2^(k-1-i)+1
        // blocks (long).  Long contributes (long.len() - 1) extra
        // blocks beyond the short choice.
        let level_extra = long.len() - 1;
        // Greedy: take long when there's room; this matches the
        // standard k-bit binary representation of `extra` from MSB.
        let bit = 1usize << (k - 1 - i);
        if extra >= bit {
            for b in long {
                out.extend_from_slice(b);
            }
            extra -= level_extra;
        } else {
            out.extend_from_slice(short);
        }
    }
    assert_eq!(extra, 0);
    out
}

/// Catalogue intermediate states of `message` under `weak_hash`.
/// Returns `state → block_index_after_consuming` map for indices
/// `>= min_index`.
pub fn catalogue_states(
    message: &[u8],
    iv: &[u8],
    b_bytes: usize,
    min_index: usize,
) -> HashMap<Vec<u8>, usize> {
    let mut state = iv.to_vec();
    let mut out = HashMap::new();
    for (i, chunk) in message.chunks_exact(16).enumerate() {
        let mut block = [0u8; 16];
        block.copy_from_slice(chunk);
        state = compress(&state, &block, b_bytes);
        if i + 1 >= min_index {
            // Map state -> block index *after* this block.
            out.insert(state.clone(), i + 1);
        }
    }
    out
}

/// Find a "bridge" single block from `state` whose compress output
/// lies in `catalogue`.  Returns `(bridge_block, matched_index)`.
pub fn find_bridge(
    state: &[u8],
    b_bytes: usize,
    catalogue: &HashMap<Vec<u8>, usize>,
) -> ([u8; 16], usize) {
    let mut counter: u64 = 0;
    loop {
        let mut block = [0u8; 16];
        block[..8].copy_from_slice(&counter.to_le_bytes());
        block[15] = 0xCC;
        let next = compress(state, &block, b_bytes);
        if let Some(&idx) = catalogue.get(&next) {
            return (block, idx);
        }
        counter += 1;
        if counter > (1u64 << 24) {
            panic!("find_bridge: search blew its budget");
        }
    }
}

/// End-to-end second-pre-image: given a long target message,
/// produce a different message of the same block length whose
/// hash matches.
pub fn second_preimage(
    target: &[u8],
    iv: &[u8],
    b_bytes: usize,
) -> Vec<u8> {
    let n_blocks = target.len() / 16;
    // Pick the smallest k with 2^k >= n_blocks.
    let k = (n_blocks as f64).log2().ceil() as usize;
    let (levels, expand_state) = build_expandable(k, iv, b_bytes);
    // Catalogue states from block index k+1 onwards (the
    // expandable can produce that many blocks at minimum).
    let catalogue = catalogue_states(target, iv, b_bytes, k + 1);
    let (bridge, j) = find_bridge(&expand_state, b_bytes, &catalogue);
    // Assemble the expandable message to length j-1 blocks
    // (because bridge contributes 1 and we want total length =
    // n_blocks · 16 bytes).
    let prefix = assemble_expandable(&levels, j - 1);
    let mut forged = prefix;
    forged.extend_from_slice(&bridge);
    forged.extend_from_slice(&target[j * 16..]);
    debug_assert_eq!(forged.len(), target.len());
    forged
}

pub fn run() -> Report {
    let mut r = Report::new(53, "Kelsey-Schneier Expandable Messages");
    // 16-bit hash, target message of 2^k = 16 blocks.
    let b_bytes = 2;
    let iv = vec![0u8; b_bytes];
    let k = 4;
    let n_blocks = 1usize << k;
    // A "long" target message.  Content doesn't matter — fill with
    // a counter so blocks are distinct.
    let mut target = Vec::with_capacity(n_blocks * 16);
    for i in 0..n_blocks {
        target.extend_from_slice(&(i as u128).to_le_bytes());
    }
    let target_hash = weak_hash(&target, &iv, b_bytes);
    r.line(format!(
        "Target message: {} blocks ({} bytes), hash = {}",
        n_blocks,
        target.len(),
        hex::encode(&target_hash)
    ));

    let forged = second_preimage(&target, &iv, b_bytes);
    let forged_hash = weak_hash(&forged, &iv, b_bytes);
    r.line(format!(
        "Forged message: {} blocks, hash = {}",
        forged.len() / 16,
        hex::encode(&forged_hash)
    ));
    r.line(format!("Same length   : {}", forged.len() == target.len()));
    r.line(format!("Same hash     : {}", forged_hash == target_hash));
    r.line(format!("Different msg : {}", forged != target));
    assert_eq!(forged_hash, target_hash);
    assert_eq!(forged.len(), target.len());
    assert_ne!(forged, target);
    r.succeed()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expandable_message_hashes_consistently() {
        let iv = vec![0u8; 2];
        let (levels, end_state) = build_expandable(4, &iv, 2);
        // Every assembled length lands on `end_state`.
        for len in 4..=(4 + 15) {
            let m = assemble_expandable(&levels, len);
            assert_eq!(m.len() / 16, len);
            assert_eq!(weak_hash(&m, &iv, 2), end_state);
        }
    }

    #[test]
    fn second_preimage_attack_lands() {
        let b_bytes = 2;
        let iv = vec![0u8; 2];
        let k = 4;
        let n = 1usize << k;
        let mut m = Vec::new();
        for i in 0..n {
            m.extend_from_slice(&(i as u128).to_le_bytes());
        }
        let h = weak_hash(&m, &iv, b_bytes);
        let forged = second_preimage(&m, &iv, b_bytes);
        assert_ne!(forged, m);
        assert_eq!(forged.len(), m.len());
        assert_eq!(weak_hash(&forged, &iv, b_bytes), h);
    }
}
