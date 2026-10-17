//! **STARK** (Ben-Sasson, Bentov, Horesh, Riabzev 2018) — scalable
//! transparent argument of knowledge.
//!
//! No trusted setup.  No pairings.  Post-quantum (collision-
//! resistant-hash assumption only).  Proof size `O(log² N)` for an
//! `N`-step computation.
//!
//! ## Components
//!
//! 1. **AIR** (Algebraic Intermediate Representation): a
//!    computation as a sequence of trace rows `t_0, t_1, …` over a
//!    finite field `F`, where consecutive rows are related by a
//!    transition polynomial constraint `C(t_i, t_{i+1}) = 0`.
//! 2. **Trace LDE**: extend the trace from the original domain
//!    `H` (size `n`, multiplicative) to a coset of a larger domain
//!    `L` (size `n · β`, the blow-up factor) via low-degree
//!    extension.
//! 3. **Composition polynomial** `C(X)` evaluated on `L`, equal
//!    to zero on `H` iff the trace is valid.
//! 4. **FRI** (Fast Reed-Solomon IOP of Proximity, Ben-Sasson et
//!    al. 2018): prove that `C(X) / Z_H(X)` is a polynomial of
//!    low degree.  This is the soundness-bearing step.
//! 5. **Merkle commitment** to all evaluation tables; verifier
//!    queries random positions and the prover opens with Merkle
//!    paths.
//!
//! ## What this module ships
//!
//! - [`Air`] — the algebraic intermediate representation: trace
//!   width, transition constraint, public boundary constraints.
//! - [`fri_commit`] / [`fri_verify`] — minimal FRI prover/verifier
//!   for the proximity-test step.
//! - [`prove`] / [`verify`] — the spot-check argument described under
//!   "Honest scope" below.
//!
//! ## Worked example
//!
//! A 1-column AIR `trace[i+1] = transition(trace[i])` with public
//! boundaries `trace[0]` and `trace[n−1]` (the tests use doubling).
//!
//! ## Honest scope: this is a spot-check argument, not a STARK
//!
//! There is no low-degree extension and no composition polynomial:
//! FRI runs on the trace itself with no blow-up, and every table of
//! `n` values on an `n`-point domain has degree `< n`, so FRI here
//! proves **nothing** about the trace.  What the verifier does check:
//!
//! - both boundary rows, opened against the trace commitment;
//! - the transition between row `i` and row `i + 1` at every queried
//!   row `i`, with the query positions re-derived from the
//!   Fiat-Shamir transcript (the prover cannot choose them) and every
//!   Merkle opening bound to its position.
//!
//! A trace with `k` broken transitions survives `q` queries with
//! probability about `(1 − k/n)^q`, and a prover can grind its
//! commitment to re-roll the positions, so a single bad row is caught
//! only with probability about `q/n` per attempt.  That is spot-check
//! soundness, not succinct soundness.  Use [`super::stark_v2`] (coset
//! LDE, composition polynomial, degree-bounded FRI) for anything that
//! needs a proof; neither module is audited, and neither may back a
//! research claim.

use super::merkle::{merkle_verify, MerkleProof, MerkleTree};
use super::polynomial::{fr_add, fr_mul, fr_sub};
use crate::bls12_381::fq::scalar_modulus;
use crate::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::One;

const FS_TAG: &str = "ZK-STARK/v1";

// ── Algebraic Intermediate Representation ─────────────────────────

/// An AIR for a 1-column trace with a polynomial transition:
/// `trace[i+1] = transition(trace[i])` for a chosen scalar-to-scalar
/// `transition` closure (specified by the prover and verifier via
/// the AIR's constraint polynomial).
///
/// For Fibonacci with two columns we expose a different API; this
/// module ships the 1-column case as a pedagogical anchor and
/// the verification logic generalizes naturally.
#[derive(Clone)]
pub struct Air {
    /// Trace length (number of rows).
    pub n: usize,
    /// Boundary constraint: trace[0] must equal this value.
    pub boundary_start: BigUint,
    /// Boundary constraint: trace[n−1] must equal this value
    /// (the "public output").
    pub boundary_end: BigUint,
    /// Transition: given trace[i], the expected trace[i+1].
    /// In a real STARK this is a polynomial constraint; for the
    /// 1-column toy AIR we use a closure.
    pub transition: std::sync::Arc<dyn Fn(&BigUint) -> BigUint + Send + Sync>,
}

impl Air {
    /// Build a trace satisfying the boundary and transition
    /// constraints (called by the prover).
    pub fn execute(&self) -> Vec<BigUint> {
        let mut trace = Vec::with_capacity(self.n);
        trace.push(self.boundary_start.clone());
        for i in 1..self.n {
            trace.push((self.transition)(&trace[i - 1]));
        }
        trace
    }

    /// Verify a candidate trace satisfies the AIR's constraints.
    pub fn check_trace(&self, trace: &[BigUint]) -> bool {
        if trace.len() != self.n {
            return false;
        }
        if trace[0] != self.boundary_start {
            return false;
        }
        if trace[self.n - 1] != self.boundary_end {
            return false;
        }
        for i in 0..(self.n - 1) {
            let expected = (self.transition)(&trace[i]);
            if expected != trace[i + 1] {
                return false;
            }
        }
        true
    }
}

// ── FRI (Fast Reed-Solomon IOP of Proximity) ──────────────────────

/// A FRI commitment to a polynomial: a series of Merkle roots
/// (one per folding layer) plus the final constant.
#[derive(Clone, Debug)]
pub struct FriCommitment {
    /// Merkle root of each layer's evaluation table.
    pub layer_roots: Vec<[u8; 32]>,
    /// The final-layer constant (degree-0 polynomial value).
    pub final_constant: BigUint,
    /// Per-layer evaluation tables (kept by the prover for opening).
    pub layers: Vec<Vec<BigUint>>,
    /// Per-layer Merkle trees.
    pub layer_trees: Vec<MerkleTree>,
}

/// A FRI query opening: evaluations at the queried position and
/// its sibling at each layer, plus Merkle inclusion proofs.
#[derive(Clone, Debug)]
pub struct FriQueryOpening {
    /// `(layer_eval, sibling_eval, merkle_proof_for_eval, merkle_proof_for_sibling)`.
    pub openings: Vec<(BigUint, BigUint, MerkleProof, MerkleProof)>,
}

fn scalar_to_bytes(s: &BigUint) -> Vec<u8> {
    let mut buf = [0u8; 32];
    let bytes = s.to_bytes_be();
    buf[32 - bytes.len()..].copy_from_slice(&bytes);
    buf.to_vec()
}

/// **FRI commit**: take a polynomial as evaluations on a power-of-2
/// domain and produce the folding-layer commitments.  Each layer
/// halves the evaluation count by combining pairs:
/// `f(x) → α · f_e(x²) + f_o(x²)` for challenge `α`, where
/// `f_e, f_o` are the even/odd halves.
pub fn fri_commit(
    initial_evals: &[BigUint],
    initial_domain: &[BigUint],
    transcript_state: &mut Vec<u8>,
) -> FriCommitment {
    fri_commit_to(initial_evals, initial_domain, transcript_state, 1)
}

/// **FRI commit with a degree bound**: fold until a layer has
/// `final_size` entries.  For evaluations of a polynomial of degree
/// `< initial_len / final_size` the final layer is constant, which is
/// what [`fri_verify`] checks at every query; folding to `final_size =
/// 1` checks no degree bound at all.
pub fn fri_commit_to(
    initial_evals: &[BigUint],
    initial_domain: &[BigUint],
    transcript_state: &mut Vec<u8>,
    final_size: usize,
) -> FriCommitment {
    assert!(initial_evals.len().is_power_of_two());
    assert!(final_size.is_power_of_two() && final_size <= initial_evals.len());
    let mut layers = vec![initial_evals.to_vec()];
    let mut domains = vec![initial_domain.to_vec()];
    let mut layer_trees = Vec::new();
    let mut layer_roots = Vec::new();

    // Layer 0 commitment.
    let leaves: Vec<Vec<u8>> = initial_evals.iter().map(scalar_to_bytes).collect();
    let tree = MerkleTree::from_leaves(&leaves);
    layer_roots.push(tree.root());
    layer_trees.push(tree);

    while layers.last().unwrap().len() > final_size {
        // Fiat-Shamir challenge α for this folding step.
        transcript_state.extend_from_slice(&layer_roots.last().unwrap()[..]);
        let alpha_bytes = sha256(transcript_state);
        transcript_state.extend_from_slice(&alpha_bytes);
        let alpha = BigUint::from_bytes_be(&alpha_bytes) % scalar_modulus();

        let cur_evals = layers.last().unwrap();
        let cur_domain = domains.last().unwrap();
        let half = cur_evals.len() / 2;
        let mut next_evals = Vec::with_capacity(half);
        let mut next_domain = Vec::with_capacity(half);
        for i in 0..half {
            // f(x_i), f(-x_i) — but in multiplicative domain, the
            // pair is index i and i+half (where x_{i+half} = -x_i).
            let f_x = &cur_evals[i];
            let f_neg_x = &cur_evals[i + half];
            // f(x) = (f_e(x²) + x · f_o(x²)) → combine with α:
            //   new(x²) = (f(x) + f(-x))/2 + α · (f(x) - f(-x))/(2x).
            // For simplicity (and matching most STARK implementations
            // that use multiplicative subgroups), we use:
            //   new(x²) = ((f(x) + f(-x)) + α · (f(x) − f(-x)) · x⁻¹) / 2.
            let sum = fr_add(f_x, f_neg_x);
            let diff = fr_sub(f_x, f_neg_x);
            let x_inv = super::polynomial::fr_inv(&cur_domain[i]).expect("non-zero domain element");
            let alpha_term = fr_mul(&alpha, &fr_mul(&diff, &x_inv));
            let two_inv =
                super::polynomial::fr_inv(&BigUint::from(2u32)).expect("2 is invertible in F_r");
            let folded = fr_mul(&fr_add(&sum, &alpha_term), &two_inv);
            next_evals.push(folded);
            next_domain.push(fr_mul(&cur_domain[i], &cur_domain[i]));
        }
        let leaves: Vec<Vec<u8>> = next_evals.iter().map(scalar_to_bytes).collect();
        let tree = MerkleTree::from_leaves(&leaves);
        layer_roots.push(tree.root());
        layer_trees.push(tree);
        layers.push(next_evals);
        domains.push(next_domain);
    }
    let final_constant = layers.last().unwrap()[0].clone();
    FriCommitment {
        layer_roots,
        final_constant,
        layers,
        layer_trees,
    }
}

/// FRI query: open the queried position and its sibling at each
/// layer.
pub fn fri_query(commitment: &FriCommitment, mut position: usize) -> FriQueryOpening {
    let mut openings = Vec::new();
    for layer in 0..(commitment.layers.len() - 1) {
        let evals = &commitment.layers[layer];
        let half = evals.len() / 2;
        let pos_main = position % evals.len();
        let pos_sibling = if pos_main < half {
            pos_main + half
        } else {
            pos_main - half
        };
        let proof_main = commitment.layer_trees[layer]
            .proof(pos_main)
            .expect("position in range");
        let proof_sibling = commitment.layer_trees[layer]
            .proof(pos_sibling)
            .expect("sibling in range");
        openings.push((
            evals[pos_main].clone(),
            evals[pos_sibling].clone(),
            proof_main,
            proof_sibling,
        ));
        position = pos_main.min(pos_sibling); // next layer index
    }
    FriQueryOpening { openings }
}

/// The transcript state after [`fri_commit_to`] has absorbed `commitment`
/// (each folded layer's root and its challenge), so a verifier can
/// derive the query positions exactly as the prover did.
pub fn fri_transcript_after(
    transcript_state_initial: &[u8],
    commitment: &FriCommitment,
) -> Vec<u8> {
    let mut transcript = transcript_state_initial.to_vec();
    for root in &commitment.layer_roots[..commitment.layer_roots.len().saturating_sub(1)] {
        transcript.extend_from_slice(&root[..]);
        let alpha_bytes = sha256(&transcript);
        transcript.extend_from_slice(&alpha_bytes);
    }
    transcript
}

/// Fiat-Shamir query positions in `[0, domain_size)`, appended to
/// `transcript` (shared by provers and verifiers).
pub fn query_positions(transcript: &mut Vec<u8>, count: usize, domain_size: usize) -> Vec<usize> {
    let mut positions = Vec::with_capacity(count);
    for q in 0..count {
        transcript.extend_from_slice(&(q as u64).to_be_bytes());
        let h = sha256(transcript);
        transcript.extend_from_slice(&h);
        let pos = (BigUint::from_bytes_be(&h) % BigUint::from(domain_size as u64))
            .iter_u64_digits()
            .next()
            .unwrap_or(0) as usize;
        positions.push(pos);
    }
    positions
}

/// FRI verify: replay the Fiat-Shamir challenges, check every Merkle
/// path at the position the folding dictates, and verify the folding
/// is consistent down to the final layer, which must be constant.
/// The caller decides how many queries it requires and at which
/// positions; this checks only the openings it is given.
pub fn fri_verify(
    commitment: &FriCommitment,
    initial_domain: &[BigUint],
    queries: &[(usize, FriQueryOpening)],
    transcript_state_initial: &[u8],
) -> bool {
    if commitment.layer_roots.is_empty() {
        return false;
    }
    let folds = commitment.layer_roots.len() - 1;
    // Recompute the α challenges.
    let mut transcript = transcript_state_initial.to_vec();
    let mut alphas = Vec::new();
    for root in &commitment.layer_roots[..commitment.layer_roots.len() - 1] {
        transcript.extend_from_slice(&root[..]);
        let alpha_bytes = sha256(&transcript);
        transcript.extend_from_slice(&alpha_bytes);
        alphas.push(BigUint::from_bytes_be(&alpha_bytes) % scalar_modulus());
    }

    for (init_position, opening) in queries {
        // One opening per folding step: no layer may be skipped.
        if opening.openings.len() != folds || *init_position >= initial_domain.len() {
            return false;
        }
        // Track the current domain size and position.
        let mut layer_size = initial_domain.len();
        let mut cur_x = initial_domain[init_position % initial_domain.len()].clone();
        let mut cur_pos = init_position % initial_domain.len();
        for (layer_idx, (eval_main, eval_sibling, proof_main, proof_sibling)) in
            opening.openings.iter().enumerate()
        {
            // Check Merkle paths for both eval and sibling against
            // the layer's committed root, at the positions the folding
            // dictates (an opening of some other leaf proves nothing).
            let half_here = layer_size / 2;
            let sibling_pos = if cur_pos < half_here {
                cur_pos + half_here
            } else {
                cur_pos - half_here
            };
            if proof_main.leaf_index != cur_pos || proof_sibling.leaf_index != sibling_pos {
                return false;
            }
            let root = &commitment.layer_roots[layer_idx];
            if !merkle_verify(&scalar_to_bytes(eval_main), proof_main, root) {
                return false;
            }
            if !merkle_verify(&scalar_to_bytes(eval_sibling), proof_sibling, root) {
                return false;
            }
            // Verify the folding step produces the next layer's
            // value at the corresponding position.
            let alpha = &alphas[layer_idx];
            // x⁻¹ where x is the queried domain value at this layer.
            let x_inv = super::polynomial::fr_inv(&cur_x).expect("non-zero domain element");
            let half = layer_size / 2;
            // eval_main is always f(cur_x); eval_sibling is always
            // f(-cur_x) — regardless of whether cur_pos is in the
            // low or high half.
            let f_x = eval_main;
            let f_neg_x = eval_sibling;
            let sum = fr_add(f_x, f_neg_x);
            let diff = fr_sub(f_x, f_neg_x);
            let alpha_term = fr_mul(alpha, &fr_mul(&diff, &x_inv));
            let two_inv = super::polynomial::fr_inv(&BigUint::from(2u32)).unwrap();
            let folded_expected = fr_mul(&fr_add(&sum, &alpha_term), &two_inv);
            // The next layer's value at position min(cur_pos, cur_pos+half)
            // should equal folded_expected.
            if layer_idx + 1 < opening.openings.len() {
                // Next layer's claimed eval_main.
                let next_eval = &opening.openings[layer_idx + 1].0;
                if next_eval != &folded_expected {
                    return false;
                }
            } else {
                // Last layer: check folded equals the FRI commitment's
                // final_constant.
                if folded_expected != commitment.final_constant {
                    return false;
                }
            }
            // Advance: domain squares, layer size halves.
            cur_x = fr_mul(&cur_x, &cur_x);
            cur_pos = cur_pos.min(if cur_pos >= half {
                cur_pos - half
            } else {
                cur_pos
            });
            layer_size = half;
        }
    }
    true
}

// ── STARK end-to-end (toy: 1-column trace) ─────────────────────────

/// A STARK proof: commitments + queries.
#[derive(Clone, Debug)]
pub struct StarkProof {
    pub trace_commitment_root: [u8; 32],
    pub fri: FriCommitment,
    /// Trace row at each Fiat-Shamir query position.
    pub trace_queries: Vec<(usize, BigUint, MerkleProof)>,
    /// The successor row `(position + 1) mod n` of each query.
    pub trace_next: Vec<(usize, BigUint, MerkleProof)>,
    /// Rows `0` and `n − 1`, checked against the public boundaries.
    pub boundary_openings: Vec<(usize, BigUint, MerkleProof)>,
    pub fri_queries: Vec<(usize, FriQueryOpening)>,
}

/// **Prove** that `air.execute()` produces a trace satisfying the
/// AIR's constraints, as the spot-check argument of the module docs:
/// the trace is committed via a Merkle tree, the boundary rows are
/// opened, and the verifier checks the transition at `n_queries`
/// Fiat-Shamir positions.  FRI on the trace adds no soundness here.
pub fn prove(air: &Air, n_queries: usize) -> Option<StarkProof> {
    let trace = air.execute();
    let n = trace.len();
    if !n.is_power_of_two() {
        return None;
    }
    let omega = super::plonk::primitive_nth_root(n);
    let mut domain = Vec::with_capacity(n);
    let mut acc = BigUint::one();
    for _ in 0..n {
        domain.push(acc.clone());
        acc = fr_mul(&acc, &omega);
    }

    // Commit to the trace via Merkle.
    let leaves: Vec<Vec<u8>> = trace.iter().map(scalar_to_bytes).collect();
    let trace_tree = MerkleTree::from_leaves(&leaves);
    let trace_root = trace_tree.root();

    // Start the Fiat-Shamir transcript with the trace commitment.
    let mut transcript_state = FS_TAG.as_bytes().to_vec();
    transcript_state.extend_from_slice(&trace_root[..]);

    // FRI commit on the trace polynomial evaluations.
    let fri = fri_commit(&trace, &domain, &mut transcript_state);

    // Query phase: pick n_queries positions via Fiat-Shamir.
    let query_positions = query_positions(&mut transcript_state, n_queries, n);

    // Open the trace at each queried position, its successor, and the
    // two boundary rows.
    let open = |pos: usize| -> (usize, BigUint, MerkleProof) {
        (
            pos,
            trace[pos].clone(),
            trace_tree.proof(pos).expect("position in range"),
        )
    };
    let trace_queries: Vec<(usize, BigUint, MerkleProof)> =
        query_positions.iter().map(|&pos| open(pos)).collect();
    let trace_next: Vec<(usize, BigUint, MerkleProof)> = query_positions
        .iter()
        .map(|&pos| open((pos + 1) % n))
        .collect();
    let boundary_openings = vec![open(0), open(n - 1)];

    // Open FRI at the same positions.
    let fri_queries: Vec<(usize, FriQueryOpening)> = query_positions
        .iter()
        .map(|&pos| (pos, fri_query(&fri, pos)))
        .collect();

    Some(StarkProof {
        trace_commitment_root: trace_root,
        fri,
        trace_queries,
        trace_next,
        boundary_openings,
        fri_queries,
    })
}

/// **Verify** the spot-check argument (module docs, "Honest scope").
///
/// Rejects unless: there is at least one query; every query position
/// is the one the transcript dictates; every opening is a Merkle path
/// at its own position; both boundary rows hold the public values;
/// the transition holds between each queried row and its successor;
/// and the FRI openings (on the committed trace) fold consistently.
pub fn verify(air: &Air, proof: &StarkProof) -> bool {
    let n = air.n;
    if !n.is_power_of_two() || n < 2 {
        return false;
    }
    let q = proof.trace_queries.len();
    if q == 0
        || proof.trace_next.len() != q
        || proof.fri_queries.len() != q
        || proof.boundary_openings.len() != 2
    {
        return false;
    }
    // FRI runs on the committed trace, so its first layer is the trace.
    if proof.fri.layer_roots.first() != Some(&proof.trace_commitment_root) {
        return false;
    }
    let opened = |(pos, val, path): &(usize, BigUint, MerkleProof), at: usize| -> bool {
        *pos == at
            && path.leaf_index == at
            && merkle_verify(&scalar_to_bytes(val), path, &proof.trace_commitment_root)
    };

    // Replay Fiat-Shamir and derive the positions the prover had to use.
    let mut transcript_state = FS_TAG.as_bytes().to_vec();
    transcript_state.extend_from_slice(&proof.trace_commitment_root[..]);
    let fri_transcript_initial = transcript_state.clone();
    let mut after_fri = fri_transcript_after(&transcript_state, &proof.fri);
    let positions = query_positions(&mut after_fri, q, n);

    // Boundary rows.
    let (first, last) = (&proof.boundary_openings[0], &proof.boundary_openings[1]);
    if !opened(first, 0)
        || first.1 != air.boundary_start
        || !opened(last, n - 1)
        || last.1 != air.boundary_end
    {
        return false;
    }

    // Transition at every queried row.
    for (i, &pos) in positions.iter().enumerate() {
        let (cur, next) = (&proof.trace_queries[i], &proof.trace_next[i]);
        if !opened(cur, pos) || !opened(next, (pos + 1) % n) || proof.fri_queries[i].0 != pos {
            return false;
        }
        if pos + 1 < n && (air.transition)(&cur.1) != next.1 {
            return false;
        }
    }

    // Domain (for the FRI folding check).
    let omega = super::plonk::primitive_nth_root(n);
    let mut domain = Vec::with_capacity(n);
    let mut acc = BigUint::one();
    for _ in 0..n {
        domain.push(acc.clone());
        acc = fr_mul(&acc, &omega);
    }
    fri_verify(
        &proof.fri,
        &domain,
        &proof.fri_queries,
        &fri_transcript_initial,
    )
}

#[cfg(test)]
mod tests {
    use super::super::polynomial::Poly;
    use super::*;
    use std::sync::Arc;

    /// Simple AIR: trace[i+1] = 2 · trace[i] (geometric sequence).
    fn doubling_air(n: usize, start: u64) -> Air {
        let trace_target = start * (1u64 << (n - 1));
        let transition: Arc<dyn Fn(&BigUint) -> BigUint + Send + Sync> =
            Arc::new(|x| fr_mul(x, &BigUint::from(2u32)));
        Air {
            n,
            boundary_start: BigUint::from(start),
            boundary_end: BigUint::from(trace_target),
            transition,
        }
    }

    /// AIR execution produces the expected sequence.
    #[test]
    fn air_execution_doubling() {
        let air = doubling_air(8, 1);
        let trace = air.execute();
        let expected: Vec<BigUint> = (0..8).map(|i| BigUint::from(1u64 << i)).collect();
        assert_eq!(trace, expected);
        assert!(air.check_trace(&trace));
    }

    /// AIR check rejects a tampered trace.
    #[test]
    fn air_check_rejects_bad_trace() {
        let air = doubling_air(8, 1);
        let mut trace = air.execute();
        trace[3] = BigUint::from(99u32);
        assert!(!air.check_trace(&trace));
    }

    /// FRI commit + verify on a simple polynomial evaluation.
    #[test]
    fn fri_commit_and_verify_polynomial() {
        let n = 16;
        let omega = super::super::plonk::primitive_nth_root(n);
        let mut domain = Vec::with_capacity(n);
        let mut acc = BigUint::one();
        for _ in 0..n {
            domain.push(acc.clone());
            acc = fr_mul(&acc, &omega);
        }
        // p(X) = 1 + 2X + 3X²  (degree 2 < n)
        let p = Poly::from_coeffs(vec![
            BigUint::from(1u32),
            BigUint::from(2u32),
            BigUint::from(3u32),
        ]);
        let evals: Vec<BigUint> = domain.iter().map(|x| p.evaluate(x)).collect();

        let mut transcript_state = b"FRI-TEST".to_vec();
        let initial_transcript = transcript_state.clone();
        let commitment = fri_commit(&evals, &domain, &mut transcript_state);

        // Query at a few positions.
        let positions = [0usize, 3, 7];
        let queries: Vec<(usize, FriQueryOpening)> = positions
            .iter()
            .map(|&p| (p, fri_query(&commitment, p)))
            .collect();

        assert!(fri_verify(
            &commitment,
            &domain,
            &queries,
            &initial_transcript
        ));
    }

    /// **End-to-end STARK** on the doubling AIR.
    #[test]
    fn stark_prove_verify_doubling() {
        let air = doubling_air(8, 1);
        let proof = prove(&air, 3).unwrap();
        assert!(verify(&air, &proof));
    }

    // Negative controls: each builds a proof of a FALSE statement (or
    // a mangled proof) and requires the verifier to reject it.  These
    // are the mutants a 2026-10-16 sheriff audit found accepted
    // (crypto-autoresearcher research/sheriff/SHF-20261016-3caa6f).

    /// 1 doubled seven times is 128, not 12345: a prover that breaks
    /// one transition to hit the false output is caught when a query
    /// lands on that row (64 queries over 8 rows).
    #[test]
    fn rejects_false_output_with_broken_transition() {
        let false_air = doubling_air(8, 1);
        let false_air = Air {
            boundary_end: BigUint::from(12345u32),
            ..false_air
        };
        assert!(!false_air.check_trace(&false_air.execute()));
        let cheat = Air {
            transition: Arc::new(|x: &BigUint| {
                if *x == BigUint::from(64u32) {
                    BigUint::from(12345u32)
                } else {
                    fr_mul(x, &BigUint::from(2u32))
                }
            }),
            ..false_air.clone()
        };
        let forged = prove(&cheat, 64).unwrap();
        assert!(!verify(&false_air, &forged));
    }

    /// A proof with no queries proves nothing.
    #[test]
    fn rejects_empty_proof() {
        let air = doubling_air(8, 1);
        let mut proof = prove(&air, 8).unwrap();
        proof.trace_queries.clear();
        proof.trace_next.clear();
        proof.fri_queries.clear();
        assert!(!verify(&air, &proof));
    }

    /// An honest proof (output 128) presented for another output.
    #[test]
    fn rejects_proof_for_another_output() {
        let proof = prove(&doubling_air(8, 1), 8).unwrap();
        let other = Air {
            boundary_end: BigUint::from(256u32),
            ..doubling_air(8, 1)
        };
        assert!(!verify(&other, &proof));
    }

    /// The prover may not choose its query positions.
    #[test]
    fn rejects_moved_query_position() {
        let air = doubling_air(8, 1);
        let mut proof = prove(&air, 8).unwrap();
        let (pos, _, _) = proof.trace_queries[0].clone();
        let moved = (pos + 1) % 8;
        proof.trace_queries[0] = (
            moved,
            air.execute()[moved].clone(),
            proof.trace_queries[0].2.clone(),
        );
        assert!(!verify(&air, &proof));
    }
}
