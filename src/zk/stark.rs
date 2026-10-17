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
//! **Educational and unaudited.  Not a STARK, not zero-knowledge (the
//! opened trace rows are revealed), no low-degree extension, no
//! out-of-domain (DEEP) sample.  It must not back a research claim.**
//!
//! There is no low-degree extension and no composition polynomial:
//! FRI runs on the trace itself with no blow-up, and every table of
//! `n` values on an `n`-point domain has degree `< n`, so FRI here
//! proves **nothing** about the trace.  What the verifier does check:
//!
//! - both boundary rows, opened against the trace commitment;
//! - the transition between row `i` and row `i + 1` at every queried
//!   row `i`, with the query positions re-derived from the
//!   Fiat-Shamir transcript and every Merkle opening bound to its
//!   position;
//! - that every opened value is a canonical field element (`< r`);
//! - the FRI layer count (`log₂ n` folds) and the final layer.
//!
//! **Transcript.**  Before the trace commitment the transcript absorbs
//! the whole statement: a domain-separation tag, an AIR identifier
//! (the transition evaluated at fixed public points, see
//! [`air_fingerprint`]), `n`, both boundary values, the query count `q`
//! and the FRI fold count.  It then absorbs the trace root, every FRI
//! layer root and the final constant before the query positions are
//! drawn.  The verifier takes `q` and the fold count from its own
//! arguments, never from the proof, and rejects a proof whose query
//! count, layer count or final layer differs.
//!
//! **Spot-check bound.**  Query positions are uniform on `[0, n)`; a
//! query on row `n − 1` checks no transition.  For a committed trace
//! with `k` broken transitions, one proof with a verifier-fixed `q`
//! escapes with probability `(1 − k/n)^q`, so a single bad row is
//! caught with probability `1 − (1 − 1/n)^q` per proof.  A prover can
//! grind: re-committing a different cheating trace re-rolls the
//! positions, so `G` attempts succeed with probability up to
//! `1 − (1 − (1 − k/n)^q)^G ≤ G · (1 − k/n)^q`.  Size `q` with
//! [`spot_check_queries`]: `q = 20` at `n = 8` is grindable in a few
//! hundred attempts (`tests::spot_check_at_q20_is_grindable`).  That
//! is spot-check soundness, not succinct soundness.  Use
//! [`super::stark_v2`] (coset LDE, composition polynomial,
//! degree-bounded FRI) for anything that needs a proof; neither
//! module is audited, and neither may back a research claim.

use super::merkle::{merkle_verify, MerkleProof, MerkleTree};
use super::polynomial::{fr_add, fr_mul, fr_sub};
use crate::bls12_381::fq::scalar_modulus;
use crate::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::One;
use std::sync::OnceLock;

/// Transcript domain-separation tag.  The `statement-bound` suffix marks
/// the transcript that absorbs the statement before any commitment.
const FS_TAG: &str = "ZK-STARK/v1/statement-bound";

/// Number of fixed public points at which [`air_fingerprint`] evaluates
/// a transition.
const AIR_PROBES: u64 = 4;

fn modulus() -> &'static BigUint {
    static R: OnceLock<BigUint> = OnceLock::new();
    R.get_or_init(scalar_modulus)
}

/// The 32-byte big-endian encoding of a **canonical** field element
/// (`s < r`), or `None` for any other integer (`s ≥ r`, including
/// values wider than 32 bytes).  Verifiers use this so that no field
/// element has two encodings and no malformed value can panic.
pub fn canonical_bytes(s: &BigUint) -> Option<[u8; 32]> {
    if s >= modulus() {
        return None;
    }
    let bytes = s.to_bytes_be();
    let mut buf = [0u8; 32];
    buf[32 - bytes.len()..].copy_from_slice(&bytes);
    Some(buf)
}

/// Reduce `s` mod `r` and encode it (prover side: never panics).
pub(crate) fn reduced_bytes(s: &BigUint) -> Vec<u8> {
    let reduced = s % modulus();
    canonical_bytes(&reduced)
        .expect("reduced value is canonical")
        .to_vec()
}

/// Absorb a canonical field element; `false` if `s` is not canonical.
pub(crate) fn absorb_scalar(transcript: &mut Vec<u8>, s: &BigUint) -> bool {
    match canonical_bytes(s) {
        Some(b) => {
            transcript.extend_from_slice(&b);
            true
        }
        None => false,
    }
}

/// Absorb a labelled `u64` (sizes, counts, positions).
pub(crate) fn absorb_u64(transcript: &mut Vec<u8>, label: &[u8], x: u64) {
    transcript.extend_from_slice(label);
    transcript.extend_from_slice(&x.to_be_bytes());
}

/// A fixed public probe point `H(tag ‖ j) mod r`, used to fingerprint a
/// transition closure.  The points do not depend on any proof.
pub(crate) fn probe_point(tag: &str, j: u64) -> BigUint {
    let mut m = tag.as_bytes().to_vec();
    m.extend_from_slice(&j.to_be_bytes());
    BigUint::from_bytes_be(&sha256(&m)) % modulus()
}

/// Draw a Fiat-Shamir challenge: `h = SHA-256(transcript)`, absorb `h`,
/// return `h mod r`.
pub(crate) fn draw_challenge(transcript: &mut Vec<u8>) -> BigUint {
    let h = sha256(transcript);
    transcript.extend_from_slice(&h);
    BigUint::from_bytes_be(&h) % modulus()
}

/// Is a Merkle opening `(pos, value, path)` for leaf `at` of a tree of
/// `2^depth` leaves valid against `root`?  Checks the position, the path
/// length, that the value is canonical, and the path itself.
pub(crate) fn opening_ok(
    (pos, val, path): &(usize, BigUint, MerkleProof),
    at: usize,
    depth: usize,
    root: &[u8; 32],
) -> bool {
    *pos == at && leaf_ok(val, path, at, depth, root)
}

/// [`opening_ok`] without the claimed position: `val` at leaf `at`.
fn leaf_ok(val: &BigUint, path: &MerkleProof, at: usize, depth: usize, root: &[u8; 32]) -> bool {
    if path.leaf_index != at || path.path.len() != depth || path.directions.len() != depth {
        return false;
    }
    match canonical_bytes(val) {
        Some(b) => merkle_verify(&b, path, root),
        None => false,
    }
}

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

/// **AIR identifier**: the transition evaluated (mod `r`) at `0`, `1`
/// and [`AIR_PROBES`] fixed public points `H("ZK-STARK/v1/air-probe" ‖ j)`.
/// A transition of degree `≤ 1` is determined exactly by its values at
/// `0` and `1`; for higher degree this is a fingerprint, and two
/// transitions with equal fingerprints agree at all six points.
pub fn air_fingerprint(air: &Air) -> Vec<u8> {
    let mut points = vec![BigUint::from(0u32), BigUint::one()];
    points.extend((0..AIR_PROBES).map(|j| probe_point("ZK-STARK/v1/air-probe", j)));
    let mut out = b"air:v1-unary-transition".to_vec();
    for x in &points {
        out.extend_from_slice(&reduced_bytes(&(air.transition)(x)));
    }
    out
}

/// Number of FRI folds the v1 verifier requires: `log₂ n` (down to one
/// value).  Fixed by the verifier, absorbed into the transcript.
fn v1_folds(n: usize) -> usize {
    n.trailing_zeros() as usize
}

/// The Fiat-Shamir transcript state before any commitment: tag, AIR
/// identifier, `n`, both boundary values, `q` and the FRI fold count.
/// `None` if `n` is not a power of two `≥ 2`, `q = 0`, or a boundary
/// value is not canonical.
pub fn statement_transcript(air: &Air, query_count: usize) -> Option<Vec<u8>> {
    let n = air.n;
    if !n.is_power_of_two() || n < 2 || query_count == 0 {
        return None;
    }
    let mut t = FS_TAG.as_bytes().to_vec();
    t.extend_from_slice(&air_fingerprint(air));
    absorb_u64(&mut t, b"n", n as u64);
    t.extend_from_slice(b"boundary_start");
    if !absorb_scalar(&mut t, &air.boundary_start) {
        return None;
    }
    t.extend_from_slice(b"boundary_end");
    if !absorb_scalar(&mut t, &air.boundary_end) {
        return None;
    }
    absorb_u64(&mut t, b"q", query_count as u64);
    absorb_u64(&mut t, b"fri_folds", v1_folds(n) as u64);
    t.extend_from_slice(b"trace_commitment");
    Some(t)
}

/// Queries needed so that one proof with a single broken transition
/// escapes with probability at most `2^−bits`: the least `q` with
/// `(1 − 1/n)^q ≤ 2^−bits`.  A prover making `2^g` attempts then
/// succeeds with probability at most `2^(g − bits)` (module docs).
pub fn spot_check_queries(n: usize, bits: u32) -> usize {
    assert!(n >= 2, "n must be at least 2");
    let per_query = -(1.0 - 1.0 / n as f64).log2();
    (bits as f64 / per_query).ceil() as usize
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

/// Absorb the last FRI layer: its root and the final constant (both
/// bound before any query position is drawn).
fn absorb_fri_final(transcript: &mut Vec<u8>, last_root: &[u8; 32], final_constant: &[u8]) {
    transcript.extend_from_slice(b"fri_final");
    transcript.extend_from_slice(last_root);
    transcript.extend_from_slice(final_constant);
}

/// The Merkle root of a final layer of `size` copies of `constant`.
fn constant_layer_root(constant: &[u8; 32], size: usize) -> [u8; 32] {
    MerkleTree::from_leaves(&vec![constant.to_vec(); size]).root()
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
///
/// Input values are reduced mod `r` first.  Every layer root is
/// absorbed (each folded one before its challenge), then the last root
/// and the final constant, so the query positions drawn afterwards
/// depend on the whole commitment.
pub fn fri_commit_to(
    initial_evals: &[BigUint],
    initial_domain: &[BigUint],
    transcript_state: &mut Vec<u8>,
    final_size: usize,
) -> FriCommitment {
    assert!(initial_evals.len().is_power_of_two());
    assert!(final_size.is_power_of_two() && final_size <= initial_evals.len());
    let initial: Vec<BigUint> = initial_evals.iter().map(|v| v % modulus()).collect();
    let mut layers = vec![initial];
    let mut domains = vec![initial_domain.to_vec()];
    let mut layer_trees = Vec::new();
    let mut layer_roots = Vec::new();

    // Layer 0 commitment.
    let leaves: Vec<Vec<u8>> = layers[0].iter().map(reduced_bytes).collect();
    let tree = MerkleTree::from_leaves(&leaves);
    layer_roots.push(tree.root());
    layer_trees.push(tree);

    while layers.last().unwrap().len() > final_size {
        // Fiat-Shamir challenge α for this folding step.
        transcript_state.extend_from_slice(&layer_roots.last().unwrap()[..]);
        let alpha = draw_challenge(transcript_state);

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
        let leaves: Vec<Vec<u8>> = next_evals.iter().map(reduced_bytes).collect();
        let tree = MerkleTree::from_leaves(&leaves);
        layer_roots.push(tree.root());
        layer_trees.push(tree);
        layers.push(next_evals);
        domains.push(next_domain);
    }
    let final_constant = layers.last().unwrap()[0].clone();
    absorb_fri_final(
        transcript_state,
        layer_roots.last().unwrap(),
        &reduced_bytes(&final_constant),
    );
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
/// (each folded layer's root and its challenge, then the last root and
/// the final constant), so a verifier can derive the query positions
/// exactly as the prover did.  `None` if the commitment has no layer or
/// a non-canonical final constant.
pub fn fri_transcript_after(
    transcript_state_initial: &[u8],
    commitment: &FriCommitment,
) -> Option<Vec<u8>> {
    let (last, folded) = commitment.layer_roots.split_last()?;
    let final_constant = canonical_bytes(&commitment.final_constant)?;
    let mut transcript = transcript_state_initial.to_vec();
    for root in folded {
        transcript.extend_from_slice(&root[..]);
        draw_challenge(&mut transcript);
    }
    absorb_fri_final(&mut transcript, last, &final_constant);
    Some(transcript)
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
/// is consistent down to the final layer, which must be constant: the
/// last layer root must commit to `initial_len / 2^folds` copies of
/// the final constant.  Every opened value must be canonical.  The
/// caller fixes the number of folds (by checking
/// `commitment.layer_roots.len()`), how many queries it requires and
/// at which positions; this checks only the openings it is given.
/// Returns `false`, never panics, on malformed input.
pub fn fri_verify(
    commitment: &FriCommitment,
    initial_domain: &[BigUint],
    queries: &[(usize, FriQueryOpening)],
    transcript_state_initial: &[u8],
) -> bool {
    let domain_len = initial_domain.len();
    if commitment.layer_roots.is_empty() || !domain_len.is_power_of_two() {
        return false;
    }
    let folds = commitment.layer_roots.len() - 1;
    let log_len = domain_len.trailing_zeros() as usize;
    if folds > log_len {
        return false;
    }
    // The final layer: `domain_len >> folds` copies of the constant.
    let Some(final_bytes) = canonical_bytes(&commitment.final_constant) else {
        return false;
    };
    if commitment.layer_roots[folds] != constant_layer_root(&final_bytes, domain_len >> folds) {
        return false;
    }
    // Recompute the α challenges.
    let mut transcript = transcript_state_initial.to_vec();
    let mut alphas = Vec::new();
    for root in &commitment.layer_roots[..folds] {
        transcript.extend_from_slice(&root[..]);
        alphas.push(draw_challenge(&mut transcript));
    }
    let two_inv = super::polynomial::fr_inv(&BigUint::from(2u32)).expect("2 is invertible");

    for (init_position, opening) in queries {
        // One opening per folding step: no layer may be skipped.
        if opening.openings.len() != folds || *init_position >= domain_len {
            return false;
        }
        // Track the current domain size and position.
        let mut layer_size = domain_len;
        let mut cur_x = initial_domain[*init_position].clone();
        let mut cur_pos = *init_position;
        for (layer_idx, (eval_main, eval_sibling, proof_main, proof_sibling)) in
            opening.openings.iter().enumerate()
        {
            // Check Merkle paths for both eval and sibling against
            // the layer's committed root, at the positions the folding
            // dictates (an opening of some other leaf proves nothing).
            let half = layer_size / 2;
            let depth = layer_size.trailing_zeros() as usize;
            let sibling_pos = if cur_pos < half {
                cur_pos + half
            } else {
                cur_pos - half
            };
            let root = &commitment.layer_roots[layer_idx];
            if !leaf_ok(eval_main, proof_main, cur_pos, depth, root)
                || !leaf_ok(eval_sibling, proof_sibling, sibling_pos, depth, root)
            {
                return false;
            }
            // Verify the folding step produces the next layer's
            // value at the corresponding position.
            let alpha = &alphas[layer_idx];
            // x⁻¹ where x is the queried domain value at this layer.
            let Some(x_inv) = super::polynomial::fr_inv(&cur_x) else {
                return false;
            };
            // eval_main is always f(cur_x); eval_sibling is always
            // f(-cur_x) — regardless of whether cur_pos is in the
            // low or high half.
            let sum = fr_add(eval_main, eval_sibling);
            let diff = fr_sub(eval_main, eval_sibling);
            let alpha_term = fr_mul(alpha, &fr_mul(&diff, &x_inv));
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
            if cur_pos >= half {
                cur_pos -= half;
            }
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

fn subgroup(n: usize) -> Vec<BigUint> {
    let omega = super::plonk::primitive_nth_root(n);
    let mut domain = Vec::with_capacity(n);
    let mut acc = BigUint::one();
    for _ in 0..n {
        domain.push(acc.clone());
        acc = fr_mul(&acc, &omega);
    }
    domain
}

/// **Prove** that `air.execute()` produces a trace satisfying the
/// AIR's constraints, as the spot-check argument of the module docs:
/// the trace is committed via a Merkle tree, the boundary rows are
/// opened, and the verifier checks the transition at `n_queries`
/// Fiat-Shamir positions.  FRI on the trace adds no soundness here.
/// `None` if `n` is not a power of two `≥ 2`, `n_queries = 0`, or a
/// boundary value is not canonical.
pub fn prove(air: &Air, n_queries: usize) -> Option<StarkProof> {
    prove_trace(air, &air.execute(), n_queries, 1)
}

/// The prover on an explicit trace (reduced mod `r` before it is
/// committed), folding FRI down to `fri_final_size` values.  [`prove`]
/// uses the executed trace and `fri_final_size = 1`; tests use other
/// traces and depths as cheating provers.
fn prove_trace(
    air: &Air,
    trace: &[BigUint],
    n_queries: usize,
    fri_final_size: usize,
) -> Option<StarkProof> {
    let n = air.n;
    if trace.len() != n || !fri_final_size.is_power_of_two() || fri_final_size > n {
        return None;
    }
    let mut transcript_state = statement_transcript(air, n_queries)?;
    let trace: Vec<BigUint> = trace.iter().map(|v| v % modulus()).collect();
    let domain = subgroup(n);

    // Commit to the trace via Merkle.
    let leaves: Vec<Vec<u8>> = trace.iter().map(reduced_bytes).collect();
    let trace_tree = MerkleTree::from_leaves(&leaves);
    let trace_root = trace_tree.root();
    transcript_state.extend_from_slice(&trace_root[..]);

    // FRI commit on the trace polynomial evaluations.
    let fri = fri_commit_to(&trace, &domain, &mut transcript_state, fri_final_size);

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

/// **Verify** the spot-check argument (module docs, "Honest scope")
/// with the verifier's own query count `query_count` (never the
/// proof's; size it with [`spot_check_queries`]).
///
/// Rejects unless: the proof has exactly `query_count` queries and
/// `log₂ n + 1` FRI layers; every query position is the one the
/// statement-bound transcript dictates; every opening is a Merkle path
/// at its own position with a canonical value; both boundary rows hold
/// the public values; the transition holds between each queried row
/// and its successor; and the FRI openings (on the committed trace)
/// fold consistently to the committed final layer.  Malformed input
/// returns `false`; it never panics.
pub fn verify(air: &Air, query_count: usize, proof: &StarkProof) -> bool {
    let n = air.n;
    let q = query_count;
    // Checks n, q and that the statement's values are canonical.
    let Some(mut transcript_state) = statement_transcript(air, q) else {
        return false;
    };
    if proof.trace_queries.len() != q
        || proof.trace_next.len() != q
        || proof.fri_queries.len() != q
        || proof.boundary_openings.len() != 2
    {
        return false;
    }
    // FRI depth is the verifier's: log₂ n folds down to one value.
    if proof.fri.layer_roots.len() != v1_folds(n) + 1 {
        return false;
    }
    // FRI runs on the committed trace, so its first layer is the trace.
    if proof.fri.layer_roots[0] != proof.trace_commitment_root {
        return false;
    }
    let depth = v1_folds(n);
    let root = proof.trace_commitment_root;
    let opened = |o: &(usize, BigUint, MerkleProof), at: usize| opening_ok(o, at, depth, &root);

    // Replay Fiat-Shamir and derive the positions the prover had to use.
    transcript_state.extend_from_slice(&root[..]);
    let fri_transcript_initial = transcript_state.clone();
    let Some(mut after_fri) = fri_transcript_after(&transcript_state, &proof.fri) else {
        return false;
    };
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
        if pos + 1 < n && (air.transition)(&cur.1) % modulus() != next.1 {
            return false;
        }
    }

    fri_verify(
        &proof.fri,
        &subgroup(n),
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

    /// **End-to-end STARK** on the doubling AIR, at a few query counts
    /// including one sized for 64 bits per proof.
    #[test]
    fn stark_prove_verify_doubling() {
        let air = doubling_air(8, 1);
        for q in [3, 20, spot_check_queries(8, 64)] {
            let proof = prove(&air, q).unwrap();
            assert!(verify(&air, q, &proof), "honest proof rejected at q = {q}");
        }
        let air16 = doubling_air(16, 3);
        let proof = prove(&air16, 20).unwrap();
        assert!(verify(&air16, 20, &proof));
    }

    #[test]
    fn spot_check_query_counts() {
        // (7/8)^q <= 2^-64 first at q = 333.
        assert_eq!(spot_check_queries(8, 64), 333);
        assert!((7.0f64 / 8.0).powi(333) <= 2f64.powi(-64));
        assert!((7.0f64 / 8.0).powi(332) > 2f64.powi(-64));
    }

    // Negative controls: each builds a proof of a FALSE statement (or
    // a mangled proof) and requires the verifier to reject it.  The
    // first four are the mutants a 2026-10-16 sheriff audit found
    // accepted (crypto-autoresearcher research/sheriff/SHF-20261016-3caa6f);
    // the rest port the re-audit's mutants (SHF-20261017-e4d7f6, reader
    // "stark", new-mutants/src/main.rs, rows `v1.*`).

    fn false_air() -> Air {
        Air {
            boundary_end: BigUint::from(12345u32),
            ..doubling_air(8, 1)
        }
    }

    /// The re-audit's cheating transition: 64 -> 12345, else doubling.
    fn cheat_air() -> Air {
        Air {
            transition: Arc::new(|x: &BigUint| {
                if *x == BigUint::from(64u32) {
                    BigUint::from(12345u32)
                } else {
                    fr_mul(x, &BigUint::from(2u32))
                }
            }),
            ..false_air()
        }
    }

    /// 1 doubled seven times is 128, not 12345: a prover that breaks
    /// one transition to hit the false output is caught when a query
    /// lands on that row (64 queries over 8 rows).
    #[test]
    fn rejects_false_output_with_broken_transition() {
        let false_air = false_air();
        assert!(!false_air.check_trace(&false_air.execute()));
        let forged = prove(&cheat_air(), 64).unwrap();
        assert!(!verify(&false_air, 64, &forged));
    }

    /// N7a: the same false statement proved with ONE query (the forger
    /// picks q).  The verifier uses its own q, so the proof is rejected.
    #[test]
    fn rejects_q1_forgery_verifier_fixes_query_count() {
        let f1 = prove(&cheat_air(), 1).unwrap();
        assert!(!verify(&false_air(), 20, &f1));
        assert!(!verify(&false_air(), 64, &f1));
    }

    /// N7b: a 20-query proof cut to a prefix.  Positions are bound to
    /// q through the transcript, so a prefix is not a proof for any q.
    #[test]
    fn rejects_truncated_proof_prefix() {
        let truncate = |mut p: StarkProof, k: usize| {
            p.trace_queries.truncate(k);
            p.trace_next.truncate(k);
            p.fri_queries.truncate(k);
            p
        };
        let forged = truncate(prove(&cheat_air(), 20).unwrap(), 3);
        assert!(!verify(&false_air(), 20, &forged));
        assert!(!verify(&false_air(), 3, &forged));
        // Even an honest proof's prefix is rejected at the prefix's q.
        let air = doubling_air(8, 1);
        let honest = truncate(prove(&air, 20).unwrap(), 3);
        assert!(!verify(&air, 20, &honest));
        assert!(!verify(&air, 3, &honest));
    }

    /// A proof with more or fewer queries than the verifier's q.
    #[test]
    fn rejects_query_count_mismatch() {
        let air = doubling_air(8, 1);
        let proof = prove(&air, 20).unwrap();
        assert!(verify(&air, 20, &proof));
        assert!(!verify(&air, 10, &proof));
        assert!(!verify(&air, 40, &proof));
        assert!(!verify(&air, 0, &proof));
    }

    /// A proof with no queries proves nothing.
    #[test]
    fn rejects_empty_proof() {
        let air = doubling_air(8, 1);
        let mut proof = prove(&air, 8).unwrap();
        proof.trace_queries.clear();
        proof.trace_next.clear();
        proof.fri_queries.clear();
        assert!(!verify(&air, 8, &proof));
        assert!(!verify(&air, 0, &proof));
    }

    /// An honest proof (output 128) presented for another output, or
    /// (N2b) another input.
    #[test]
    fn rejects_proof_for_another_output() {
        let proof = prove(&doubling_air(8, 1), 8).unwrap();
        let other = Air {
            boundary_end: BigUint::from(256u32),
            ..doubling_air(8, 1)
        };
        assert!(!verify(&other, 8, &proof));
        let other_start = Air {
            boundary_start: BigUint::from(3u32),
            ..doubling_air(8, 1)
        };
        assert!(!verify(&other_start, 8, &proof));
        // N4: a valid proof of 2 -> 256 presented for 1 -> 128.
        let p_other = prove(&doubling_air(8, 2), 8).unwrap();
        assert!(!verify(&doubling_air(8, 1), 8, &p_other));
    }

    /// N3a / N3b: an honest doubling proof under a tripling AIR (a
    /// false statement: 3^7 != 128) or a 16-row AIR.
    #[test]
    fn rejects_wrong_air() {
        let air = doubling_air(8, 1);
        let proof = prove(&air, 20).unwrap();
        let tripling = Air {
            transition: Arc::new(|x: &BigUint| fr_mul(x, &BigUint::from(3u32))),
            ..air.clone()
        };
        // 1 tripled seven times is 2187, not 128: the statement is false.
        assert!(!tripling.check_trace(&tripling.execute()));
        assert!(!verify(&tripling, 20, &proof));
        let air16 = Air {
            n: 16,
            boundary_end: BigUint::from(32768u32),
            ..air.clone()
        };
        assert!(!verify(&air16, 20, &proof));
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
        assert!(!verify(&air, 8, &proof));
    }

    /// First index whose position differs from query 0's.
    fn distinct_from_first(p: &StarkProof) -> usize {
        (1..p.trace_queries.len())
            .find(|&j| p.trace_queries[j].0 != p.trace_queries[0].0)
            .expect("two distinct positions")
    }

    /// N6a / N6b: reordered and duplicated queries.
    #[test]
    fn rejects_reordered_or_duplicated_query() {
        let air = doubling_air(8, 1);
        let proof = prove(&air, 20).unwrap();
        let j = distinct_from_first(&proof);
        let mut p = proof.clone();
        p.trace_queries.swap(0, j);
        p.trace_next.swap(0, j);
        p.fri_queries.swap(0, j);
        assert!(!verify(&air, 20, &p));
        let mut p = proof.clone();
        p.trace_queries[j] = p.trace_queries[0].clone();
        p.trace_next[j] = p.trace_next[0].clone();
        p.fri_queries[j] = p.fri_queries[0].clone();
        assert!(!verify(&air, 20, &p));
    }

    /// N5a / N5b / N5c: FRI with a dropped root, dropped openings, or
    /// zero folds with the transcript re-derived (the prover may not
    /// choose the FRI depth).
    #[test]
    fn rejects_truncated_or_zero_fold_fri() {
        let air = doubling_air(8, 1);
        let proof = prove(&air, 20).unwrap();
        let mut p = proof.clone();
        p.fri.layer_roots.pop();
        assert!(!verify(&air, 20, &p));
        let mut p = proof.clone();
        for fq in p.fri_queries.iter_mut() {
            fq.1.openings.pop();
        }
        assert!(!verify(&air, 20, &p));
        // Zero folds and one fold fewer, from the real prover.
        for final_size in [8, 2] {
            let p = prove_trace(&air, &air.execute(), 20, final_size).unwrap();
            assert!(!verify(&air, 20, &p), "final size {final_size} accepted");
        }
    }

    /// The last FRI root and the final constant are bound: zeroing the
    /// root (re-audit v2.M1, never checked before) or bumping the
    /// constant is rejected.
    #[test]
    fn rejects_tampered_fri_final_layer() {
        let air = doubling_air(8, 1);
        let proof = prove(&air, 20).unwrap();
        let mut p = proof.clone();
        let last = p.fri.layer_roots.len() - 1;
        p.fri.layer_roots[last] = [0u8; 32];
        assert!(!verify(&air, 20, &p));
        let mut p = proof.clone();
        p.fri.final_constant = fr_add(&p.fri.final_constant, &BigUint::one());
        assert!(!verify(&air, 20, &p));
    }

    /// Non-canonical encodings (v + r) are rejected, as are statements
    /// with non-canonical values.
    #[test]
    fn rejects_non_canonical_values() {
        let r = scalar_modulus();
        assert!(canonical_bytes(&r).is_none());
        assert!(canonical_bytes(&(&r - 1u32)).is_some());
        let air = doubling_air(8, 1);
        let proof = prove(&air, 20).unwrap();
        let mut p = proof.clone();
        p.trace_queries[0].1 += &r;
        assert!(!verify(&air, 20, &p));
        let mut p = proof.clone();
        p.boundary_openings[0].1 += &r;
        assert!(!verify(&air, 20, &p));
        let mut p = proof.clone();
        p.fri.final_constant += &r;
        assert!(!verify(&air, 20, &p));
        let shifted = Air {
            boundary_start: BigUint::one() + &r,
            ..air.clone()
        };
        assert!(!verify(&shifted, 20, &proof));
        assert!(prove(&shifted, 20).is_none());
    }

    /// An opened value of 2^300 (wider than 32 bytes) returns false
    /// instead of panicking, wherever it appears.
    #[test]
    fn oversized_values_return_false_without_panic() {
        let air = doubling_air(8, 1);
        let proof = prove(&air, 20).unwrap();
        let huge: BigUint = BigUint::one() << 300;
        let mut p = proof.clone();
        p.trace_queries[0].1 = huge.clone();
        assert!(!verify(&air, 20, &p));
        let mut p = proof.clone();
        p.trace_next[0].1 = huge.clone();
        assert!(!verify(&air, 20, &p));
        let mut p = proof.clone();
        p.boundary_openings[1].1 = huge.clone();
        assert!(!verify(&air, 20, &p));
        let mut p = proof.clone();
        p.fri_queries[0].1.openings[0].1 = huge.clone();
        assert!(!verify(&air, 20, &p));
        let mut p = proof.clone();
        p.fri.final_constant = huge.clone();
        assert!(!verify(&air, 20, &p));
        let big_statement = Air {
            boundary_end: huge,
            ..air.clone()
        };
        assert!(!verify(&big_statement, 20, &proof));
    }

    /// The re-audit's two-broken-row grinding family (rows 2 -> 3 and
    /// 6 -> 7 broken, free value v; v1.N8b).
    fn two_break_trace(i: u64) -> Vec<BigUint> {
        let v = probe_point("sheriff-SHF-20261016-3caa6f/reader-stark/v1-grind-2row/", i);
        let two = BigUint::from(2u32);
        let mut t = vec![BigUint::one(), two.clone(), BigUint::from(4u32), v.clone()];
        for _ in 0..3 {
            let next = fr_mul(t.last().unwrap(), &two);
            t.push(next);
        }
        t.push(BigUint::from(12345u32));
        t
    }

    /// N8b with the verifier's q sized by [`spot_check_queries`] for 64
    /// bits: 300 re-committed forgeries (per-attempt escape
    /// (6/8)^333 < 2^-138) are all rejected.
    #[test]
    fn grinding_fails_at_sized_query_count() {
        let q = spot_check_queries(8, 64);
        let fa = false_air();
        for i in 0..300 {
            let p = prove_trace(&fa, &two_break_trace(i), q, 1).unwrap();
            assert!(!verify(&fa, q, &p), "forgery accepted at attempt {i}");
        }
    }

    /// The spot-check bound is real, not a defect of the transcript: at
    /// q = 20 and n = 8 a two-row forgery escapes each attempt with
    /// probability (6/8)^20 ≈ 0.0032, so grinding finds an accepted
    /// proof of the FALSE statement 1 -> 12345 within a few thousand
    /// attempts.  This pins the module docs: q must be sized.
    #[test]
    fn spot_check_at_q20_is_grindable() {
        let fa = false_air();
        let found = (0..20_000u64).find(|&i| {
            let p = prove_trace(&fa, &two_break_trace(i), 20, 1).unwrap();
            verify(&fa, 20, &p)
        });
        assert!(found.is_some(), "no q = 20 forgery in 20000 attempts");
        // The same forgery is rejected by a verifier with a sized q.
        let i = found.unwrap();
        let p = prove_trace(&fa, &two_break_trace(i), 20, 1).unwrap();
        assert!(!verify(&fa, spot_check_queries(8, 64), &p));
    }
}
