//! **Binary Merkle tree** over SHA-256.
//!
//! A Merkle tree is a **vector commitment**: a single 32-byte root
//! commits to an ordered list of `N` leaves such that any single
//! leaf can be revealed and verified against the root using only
//! `⌈log₂ N⌉` "sibling" hashes (an inclusion proof).
//!
//! Used everywhere:
//!
//! - **STARK proofs** commit to evaluation tables this way.
//! - **Bitcoin / Ethereum** commit to transactions per block.
//! - **Certificate Transparency** commits to issued TLS certs.
//! - **Git** uses Merkle trees over content-addressed object DAGs.
//!
//! # Domain separation
//!
//! Naive Merkle trees are vulnerable to **second-preimage attacks**
//! (Sze 2008): an attacker can prove that an **internal node** is a
//! "leaf" of the same root.  Fix: prepend a different domain byte
//! for leaves (`0x00`) vs internal nodes (`0x01`):
//!
//! ```text
//! leaf_hash(d)       = SHA256(0x00 ‖ d)
//! node_hash(L, R)    = SHA256(0x01 ‖ L ‖ R)
//! ```
//!
//! This is RFC 6962 (Certificate Transparency).  We adopt it
//! verbatim — the same hashing rules used by Google's CT logs.
//!
//! # Unbalanced trees
//!
//! When `N` is not a power of 2, the standard convention is to
//! **duplicate the last node** at each level until pairs form.
//! Bitcoin and CT both use this convention.  Note this is
//! second-preimage-resistant only because of the leaf/node tag
//! distinction.

use crate::hash::sha256::sha256;

/// A 32-byte Merkle hash.
pub type Hash32 = [u8; 32];

const LEAF_TAG: u8 = 0x00;
const NODE_TAG: u8 = 0x01;

/// Hash a single leaf with domain separation.
fn hash_leaf(data: &[u8]) -> Hash32 {
    let mut input = Vec::with_capacity(1 + data.len());
    input.push(LEAF_TAG);
    input.extend_from_slice(data);
    sha256(&input)
}

/// Hash an internal node from its two children.
fn hash_node(left: &Hash32, right: &Hash32) -> Hash32 {
    let mut input = Vec::with_capacity(1 + 32 + 32);
    input.push(NODE_TAG);
    input.extend_from_slice(left);
    input.extend_from_slice(right);
    sha256(&input)
}

/// A Merkle inclusion proof.
///
/// `path` lists the sibling hashes from the leaf level up to the
/// root.  `directions[i]` is `true` if the sibling at level `i` is
/// the **left** child (i.e., our node is on the right).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MerkleProof {
    pub leaf_index: usize,
    pub path: Vec<Hash32>,
    pub directions: Vec<bool>,
}

/// A Merkle tree precomputed from a list of leaves.  Stores the
/// full tree so multiple inclusion proofs can be extracted cheaply.
#[derive(Clone, Debug)]
pub struct MerkleTree {
    /// Level 0 = leaves, level `levels-1` = root.
    /// `nodes[level][i]` is the `i`th hash at that level.
    /// For the empty tree, `nodes` is empty and `root()` returns
    /// `SHA-256("")` (RFC 6962 convention).
    pub nodes: Vec<Vec<Hash32>>,
}

impl MerkleTree {
    /// Build a Merkle tree from raw leaf data.  Empty list →
    /// `nodes` is empty; `root()` returns `SHA-256("")` (RFC 6962).
    pub fn from_leaves(leaves: &[Vec<u8>]) -> Self {
        if leaves.is_empty() {
            return Self { nodes: vec![] };
        }
        let leaf_hashes: Vec<Hash32> = leaves.iter().map(|d| hash_leaf(d)).collect();
        Self::from_leaf_hashes(&leaf_hashes)
    }

    /// Build a tree from already-hashed leaves (e.g., when the
    /// caller wants a different leaf-hash convention).
    pub fn from_leaf_hashes(leaf_hashes: &[Hash32]) -> Self {
        if leaf_hashes.is_empty() {
            return Self { nodes: vec![] };
        }
        let mut nodes: Vec<Vec<Hash32>> = vec![leaf_hashes.to_vec()];
        while nodes.last().unwrap().len() > 1 {
            let prev = nodes.last().unwrap();
            let mut next: Vec<Hash32> = Vec::with_capacity(prev.len().div_ceil(2));
            let mut i = 0;
            while i < prev.len() {
                let left = &prev[i];
                // Duplicate last node when level has odd length.
                let right = if i + 1 < prev.len() {
                    &prev[i + 1]
                } else {
                    &prev[i]
                };
                next.push(hash_node(left, right));
                i += 2;
            }
            nodes.push(next);
        }
        Self { nodes }
    }

    /// Number of leaves.
    pub fn num_leaves(&self) -> usize {
        self.nodes.first().map(|v| v.len()).unwrap_or(0)
    }

    /// Tree root.  For an empty tree, returns `SHA-256("")` per
    /// RFC 6962.
    pub fn root(&self) -> Hash32 {
        if self.nodes.is_empty() {
            return sha256(&[]);
        }
        *self.nodes.last().unwrap().first().unwrap()
    }

    /// Extract an inclusion proof for leaf `leaf_index`.
    pub fn proof(&self, leaf_index: usize) -> Option<MerkleProof> {
        if leaf_index >= self.num_leaves() {
            return None;
        }
        let mut path = Vec::new();
        let mut directions = Vec::new();
        let mut idx = leaf_index;
        for level in 0..(self.nodes.len() - 1) {
            let level_nodes = &self.nodes[level];
            let sibling_idx = if idx % 2 == 0 {
                // We're the LEFT child; sibling is right (idx+1 or
                // duplicate of self if at boundary).
                if idx + 1 < level_nodes.len() {
                    idx + 1
                } else {
                    idx
                }
            } else {
                // We're the RIGHT child; sibling is left.
                idx - 1
            };
            path.push(level_nodes[sibling_idx]);
            // direction: did sibling come from the LEFT?
            directions.push(sibling_idx < idx);
            idx /= 2;
        }
        Some(MerkleProof {
            leaf_index,
            path,
            directions,
        })
    }
}

/// Convenience: compute Merkle root of a leaf list.
pub fn merkle_root(leaves: &[Vec<u8>]) -> Hash32 {
    MerkleTree::from_leaves(leaves).root()
}

/// Convenience: produce an inclusion proof for `leaf_index`.
pub fn merkle_proof(leaves: &[Vec<u8>], leaf_index: usize) -> Option<MerkleProof> {
    MerkleTree::from_leaves(leaves).proof(leaf_index)
}

/// **Verify** an inclusion proof.  Given the claimed leaf data, the
/// proof, and the expected root, recompute the root by hashing
/// upward and compare.
pub fn merkle_verify(leaf_data: &[u8], proof: &MerkleProof, expected_root: &Hash32) -> bool {
    if proof.path.len() != proof.directions.len() {
        return false;
    }
    if proof.path.is_empty() && proof.leaf_index != 0 {
        return false;
    }
    for (level, sibling_is_left) in proof.directions.iter().enumerate() {
        let index_bit = ((proof.leaf_index >> level) & 1) == 1;
        if index_bit != *sibling_is_left {
            return false;
        }
    }
    let mut current = hash_leaf(leaf_data);
    for (sibling, sibling_is_left) in proof.path.iter().zip(&proof.directions) {
        current = if *sibling_is_left {
            hash_node(sibling, &current)
        } else {
            hash_node(&current, sibling)
        };
    }
    &current == expected_root
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Single-leaf tree: root is just the (tagged) hash of the leaf.
    #[test]
    fn single_leaf_tree() {
        let leaves = vec![b"hello".to_vec()];
        let tree = MerkleTree::from_leaves(&leaves);
        assert_eq!(tree.num_leaves(), 1);
        assert_eq!(tree.root(), hash_leaf(b"hello"));
        let proof = tree.proof(0).unwrap();
        assert_eq!(proof.path.len(), 0);
        assert!(merkle_verify(b"hello", &proof, &tree.root()));
    }

    /// Two-leaf tree: root = node(leaf(a), leaf(b)).
    #[test]
    fn two_leaf_tree() {
        let leaves = vec![b"a".to_vec(), b"b".to_vec()];
        let tree = MerkleTree::from_leaves(&leaves);
        let expected_root = hash_node(&hash_leaf(b"a"), &hash_leaf(b"b"));
        assert_eq!(tree.root(), expected_root);

        let proof_a = tree.proof(0).unwrap();
        assert_eq!(proof_a.path.len(), 1);
        assert_eq!(proof_a.path[0], hash_leaf(b"b"));
        assert!(
            !proof_a.directions[0],
            "sibling is on the right for index 0"
        );
        assert!(merkle_verify(b"a", &proof_a, &tree.root()));

        let proof_b = tree.proof(1).unwrap();
        assert_eq!(proof_b.path[0], hash_leaf(b"a"));
        assert!(proof_b.directions[0], "sibling is on the left for index 1");
        assert!(merkle_verify(b"b", &proof_b, &tree.root()));
    }

    /// Three-leaf tree: odd-leaf duplication kicks in.
    #[test]
    fn three_leaf_tree_with_duplication() {
        let leaves: Vec<Vec<u8>> = (0..3).map(|i| vec![i as u8]).collect();
        let tree = MerkleTree::from_leaves(&leaves);
        // Should still be verifiable for every leaf.
        for i in 0..3 {
            let proof = tree.proof(i).unwrap();
            assert!(
                merkle_verify(&leaves[i], &proof, &tree.root()),
                "proof for index {} should verify",
                i
            );
        }
        // Out-of-bounds index → None.
        assert!(tree.proof(3).is_none());
    }

    /// Eight-leaf tree (perfect binary, log₂(8) = 3 levels).
    #[test]
    fn eight_leaf_tree_perfect_binary() {
        let leaves: Vec<Vec<u8>> = (0..8u8).map(|i| vec![i, i * 2]).collect();
        let tree = MerkleTree::from_leaves(&leaves);
        assert_eq!(tree.nodes.len(), 4); // leaves + 3 internal levels = 4 levels
        for i in 0..8 {
            let proof = tree.proof(i).unwrap();
            assert_eq!(proof.path.len(), 3, "proof depth should be 3 for 8 leaves");
            assert!(merkle_verify(&leaves[i], &proof, &tree.root()));
        }
    }

    /// **Soundness**: tampering with a path hash breaks the proof.
    #[test]
    fn tampered_path_fails_verification() {
        let leaves: Vec<Vec<u8>> = (0..4u8).map(|i| vec![i]).collect();
        let tree = MerkleTree::from_leaves(&leaves);
        let mut proof = tree.proof(2).unwrap();
        // Flip a bit in the first sibling.
        proof.path[0][0] ^= 0x01;
        assert!(!merkle_verify(&leaves[2], &proof, &tree.root()));
    }

    /// **Soundness**: tampering with the leaf breaks the proof.
    #[test]
    fn tampered_leaf_fails_verification() {
        let leaves: Vec<Vec<u8>> = (0..4u8).map(|i| vec![i]).collect();
        let tree = MerkleTree::from_leaves(&leaves);
        let proof = tree.proof(1).unwrap();
        let wrong = vec![99u8];
        assert!(!merkle_verify(&wrong, &proof, &tree.root()));
    }

    /// **Domain-separation soundness**: an internal node cannot be
    /// passed off as a leaf with the same bytes.
    ///
    /// Sze (2008) — without leaf/node tags, the bytes of an
    /// intermediate hash could collide with a hypothetical leaf.
    /// With the `0x00`/`0x01` tag prefix, leaves and internal nodes
    /// occupy different parts of the hash function's input space.
    #[test]
    fn domain_separation_prevents_node_as_leaf() {
        let leaves: Vec<Vec<u8>> = vec![b"foo".to_vec(), b"bar".to_vec()];
        let tree = MerkleTree::from_leaves(&leaves);
        // Pretend the internal node L = node(leaf(foo), leaf(bar))
        // is itself a "leaf" of the same root.  Construct a "proof"
        // that the L's preimage hashes to the root with empty path.
        let internal = hash_node(&hash_leaf(b"foo"), &hash_leaf(b"bar"));
        let fake_proof = MerkleProof {
            leaf_index: 0,
            path: vec![],
            directions: vec![],
        };
        // Verifier hashes `internal` as a LEAF (tag 0x00), not as the
        // bare bytes, and gets a different value than the root.
        assert!(!merkle_verify(&internal, &fake_proof, &tree.root()));
    }

    /// Large tree: 100 leaves, all proofs verify.
    #[test]
    fn large_tree_all_proofs_verify() {
        let leaves: Vec<Vec<u8>> = (0..100u32).map(|i| i.to_be_bytes().to_vec()).collect();
        let tree = MerkleTree::from_leaves(&leaves);
        for i in 0..100 {
            let proof = tree.proof(i).unwrap();
            assert!(
                merkle_verify(&leaves[i], &proof, &tree.root()),
                "leaf {} should verify in a 100-leaf tree",
                i
            );
        }
    }

    #[test]
    fn proof_leaf_index_is_bound_to_directions() {
        let leaves: Vec<Vec<u8>> = (0..4u8).map(|i| vec![i]).collect();
        let tree = MerkleTree::from_leaves(&leaves);
        let mut proof = tree.proof(2).unwrap();
        assert!(merkle_verify(&leaves[2], &proof, &tree.root()));
        proof.leaf_index = 0;
        assert!(!merkle_verify(&leaves[2], &proof, &tree.root()));
    }

    /// Empty list: degenerate root, no leaves.
    #[test]
    fn empty_tree() {
        let tree = MerkleTree::from_leaves(&[]);
        assert_eq!(tree.num_leaves(), 0);
        assert!(tree.proof(0).is_none());
    }
}
