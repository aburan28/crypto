//! # Challenge 54 — Kelsey & Kohno "Nostradamus" / herding attack
//!
//! The setup is *commit-then-reveal*: a fortune-teller publishes
//! `H(prediction)` before the football season, then at the end
//! produces a message that, supposedly, made that very prediction.
//! With a real hash you can't choose the prediction after the fact.
//! With Merkle–Damgård you can — at a cost.
//!
//! ## The "diamond structure"
//!
//! Pre-commit phase: build a funnel of `2^k` leaf states down to a
//! single root state.  At each level pair up adjacent states and
//! find a single-block collision merging each pair.  This is the
//! "diamond" — a balanced binary tree of `k` collision-finding
//! levels.
//!
//! After the diamond exists, hash a padding block from the root
//! and *publish* the result as your commitment.
//!
//! ## Forgery phase
//!
//! Later, when the actual event has happened:
//!
//! 1. Write whatever message `P` you want as your "prediction".
//! 2. Append a single "glue" block `g` such that
//!    `compress(state_after_P, g)` lands on one of the `2^k` leaves.
//!    This is a `2^(b-k)` work birthday search.
//! 3. Walk the path from that leaf down to the root, concatenating
//!    each level's stored collision block.
//! 4. Append the padding block.  Result hashes to the published
//!    commitment.
//!
//! Total work: pre-commit `2^(k + b/2 + 1)`-ish, forgery `2^(b-k)`.

use crate::cryptopals::Report;
use crate::cryptopals::challenge52::{compress, weak_hash};
use std::collections::HashMap;

/// Find a single block from each of `a` and `b` whose outputs
/// collide.  Returns `(block_a, block_b, common_next)`.
pub fn merge_states(
    a: &[u8],
    b: &[u8],
    b_bytes: usize,
) -> ([u8; 16], [u8; 16], Vec<u8>) {
    let mut from_a: HashMap<Vec<u8>, [u8; 16]> = HashMap::new();
    let mut counter: u64 = 0;
    loop {
        let mut blk_a = [0u8; 16];
        blk_a[..8].copy_from_slice(&counter.to_le_bytes());
        blk_a[15] = 0x11;
        let next_a = compress(a, &blk_a, b_bytes);
        from_a.insert(next_a, blk_a);

        let mut blk_b = [0u8; 16];
        blk_b[..8].copy_from_slice(&counter.to_le_bytes());
        blk_b[15] = 0x22;
        let next_b = compress(b, &blk_b, b_bytes);
        if let Some(blk_a_match) = from_a.get(&next_b).cloned() {
            return (blk_a_match, blk_b, next_b);
        }
        counter += 1;
        if counter > (1u64 << 24) {
            panic!("merge_states: ran out of birthday budget");
        }
    }
}

/// Single layer of the diamond.  Each leaf-state maps to either
/// (left-child, right-child, single block to play, common parent).
#[derive(Debug, Clone)]
pub struct DiamondNode {
    pub state: Vec<u8>,
    pub block: [u8; 16], // block to play *from this node* to reach parent
}

/// Build a diamond of `2^k` leaves down to a single root.  Returns:
///
/// - `leaves`: the initial `2^k` states (each chosen as
///   `compress(iv, distinct-block)` so they're real reachable states).
/// - `tree`: layer-by-layer.  `tree[layer][i]` is the DiamondNode
///   for the `i`th state at that layer (block = move from this node
///   to its parent in the next layer).  `tree.last()` has exactly
///   one entry — that's the root.
/// - `leaf_seed_blocks`: the blocks we used to spawn each leaf from
///   `iv`, so the attacker can include them in the prefix later.
pub struct Diamond {
    pub leaves: Vec<Vec<u8>>,
    pub leaf_seed_blocks: Vec<[u8; 16]>,
    pub tree: Vec<Vec<DiamondNode>>,
    pub root_state: Vec<u8>,
}

pub fn build_diamond(k: usize, iv: &[u8], b_bytes: usize) -> Diamond {
    let n_leaves = 1usize << k;
    // Generate distinct leaf states by hashing `iv` over distinct seed blocks.
    let mut leaves = Vec::with_capacity(n_leaves);
    let mut seed_blocks = Vec::with_capacity(n_leaves);
    for i in 0..n_leaves {
        let mut block = [0u8; 16];
        block[..8].copy_from_slice(&(i as u64).to_le_bytes());
        block[15] = 0xDD;
        seed_blocks.push(block);
        leaves.push(compress(iv, &block, b_bytes));
    }
    let mut tree: Vec<Vec<DiamondNode>> = Vec::new();
    // First layer: the leaves themselves; the "block from this leaf
    // to its parent" gets filled by the next loop iteration.
    let mut current: Vec<Vec<u8>> = leaves.clone();
    let mut current_layer: Vec<DiamondNode> = current
        .iter()
        .map(|s| DiamondNode {
            state: s.clone(),
            block: [0u8; 16],
        })
        .collect();
    while current.len() > 1 {
        let mut next: Vec<Vec<u8>> = Vec::with_capacity(current.len() / 2);
        let pairs: Vec<(usize, usize)> =
            (0..current.len()).step_by(2).map(|i| (i, i + 1)).collect();
        for (i, j) in &pairs {
            let (blk_i, blk_j, parent) =
                merge_states(&current[*i], &current[*j], b_bytes);
            // Backfill the block for both children.
            current_layer[*i].block = blk_i;
            current_layer[*j].block = blk_j;
            next.push(parent);
        }
        tree.push(current_layer);
        current = next;
        current_layer = current
            .iter()
            .map(|s| DiamondNode {
                state: s.clone(),
                block: [0u8; 16],
            })
            .collect();
    }
    // current_layer at this point is the single root (block = unused).
    let root_state = current[0].clone();
    tree.push(current_layer);
    Diamond {
        leaves,
        leaf_seed_blocks: seed_blocks,
        tree,
        root_state,
    }
}

/// Walk a path from leaf index `leaf_idx` to the root, concatenating
/// blocks.  Returns the byte string `leaf-to-root-blocks`.
pub fn path_to_root(diamond: &Diamond, leaf_idx: usize) -> Vec<u8> {
    let mut out = Vec::new();
    let mut idx = leaf_idx;
    for layer in &diamond.tree[..diamond.tree.len() - 1] {
        out.extend_from_slice(&layer[idx].block);
        idx /= 2;
    }
    out
}

/// Commitment = hash `iv → leaf → root → padding-block`.  We pick
/// any fixed "commitment block" (e.g., all-`0x99`).  Once that
/// block is fixed, the forger has to land on `root_state`.
pub fn commitment(diamond: &Diamond, b_bytes: usize) -> Vec<u8> {
    let pad = [0x99u8; 16];
    compress(&diamond.root_state, &pad, b_bytes)
}

/// Birthday-search for a glue block landing on any leaf.  Returns
/// `(glue, leaf_idx)`.
pub fn find_glue(
    state_after_prediction: &[u8],
    diamond: &Diamond,
    b_bytes: usize,
) -> ([u8; 16], usize) {
    let leaf_set: HashMap<Vec<u8>, usize> = diamond
        .leaves
        .iter()
        .enumerate()
        .map(|(i, s)| (s.clone(), i))
        .collect();
    let mut counter: u64 = 0;
    loop {
        let mut block = [0u8; 16];
        block[..8].copy_from_slice(&counter.to_le_bytes());
        block[15] = 0xEE;
        let next = compress(state_after_prediction, &block, b_bytes);
        if let Some(&idx) = leaf_set.get(&next) {
            return (block, idx);
        }
        counter += 1;
        if counter > (1u64 << 24) {
            panic!("find_glue: ran out of budget");
        }
    }
}

/// Forge a full message that starts with `prediction`, has the
/// expected diamond-derived tail, and hashes to the commitment.
pub fn forge(
    diamond: &Diamond,
    prediction: &[u8],
    iv: &[u8],
    b_bytes: usize,
) -> Vec<u8> {
    // Hash the prediction up to its final state.  The prediction
    // should already be 16-byte aligned (caller's responsibility).
    assert!(prediction.len() % 16 == 0);
    let state_after_pred = weak_hash(prediction, iv, b_bytes);
    let (glue, leaf_idx) = find_glue(&state_after_pred, diamond, b_bytes);
    let mut out = prediction.to_vec();
    out.extend_from_slice(&glue);
    out.extend_from_slice(&path_to_root(diamond, leaf_idx));
    out.extend_from_slice(&[0x99u8; 16]); // commitment padding block
    out
}

pub fn run() -> Report {
    let mut r = Report::new(54, "Kelsey-Kohno Nostradamus Attack");
    let b_bytes = 2;
    let iv = vec![0u8; b_bytes];
    let k = 4;
    r.line(format!(
        "Building diamond: k = {} (2^k = {} leaves), {}-bit hash",
        k,
        1 << k,
        b_bytes * 8
    ));
    let diamond = build_diamond(k, &iv, b_bytes);
    let commit = commitment(&diamond, b_bytes);
    r.line(format!("Public commitment: {}", hex::encode(&commit)));
    r.line(format!("Root state       : {}", hex::encode(&diamond.root_state)));

    // After the event has happened, write whatever "prediction" we like.
    let prediction = b"The Wu-Tang Clan will rule the AFC East in week 17\0\0\0\0\0\0\0\0\0\0\0\0\0\0";
    // Align to 16 bytes — trim/pad as needed.
    let mut padded = prediction.to_vec();
    while padded.len() % 16 != 0 {
        padded.push(0);
    }
    let forged = forge(&diamond, &padded, &iv, b_bytes);
    let h = weak_hash(&forged, &iv, b_bytes);
    r.line("");
    r.line("Forgery:");
    r.line(format!(
        "  prediction prefix = {:?}",
        String::from_utf8_lossy(&padded).trim_end_matches('\0')
    ));
    r.line(format!("  total length      = {} bytes", forged.len()));
    r.line(format!("  forged hash       = {}", hex::encode(&h)));
    r.line(format!(
        "  matches commit    = {}",
        h == commit
    ));
    assert_eq!(h, commit);
    r.succeed()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diamond_root_matches() {
        let iv = vec![0u8; 2];
        let d = build_diamond(3, &iv, 2);
        // Each leaf, then walked up via the per-node blocks, should
        // land on the root.
        for (i, leaf) in d.leaves.iter().enumerate() {
            let path = path_to_root(&d, i);
            let mut state = leaf.clone();
            for chunk in path.chunks_exact(16) {
                let mut block = [0u8; 16];
                block.copy_from_slice(chunk);
                state = compress(&state, &block, 2);
            }
            assert_eq!(state, d.root_state, "leaf {} should reach root", i);
        }
    }

    #[test]
    fn forgery_matches_commitment() {
        let iv = vec![0u8; 2];
        let d = build_diamond(3, &iv, 2);
        let commit = commitment(&d, 2);
        let pred = b"prediction-block";
        let forged = forge(&d, pred, &iv, 2);
        assert_eq!(weak_hash(&forged, &iv, 2), commit);
    }
}
