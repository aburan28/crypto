//! **STARK v2**: a scalable transparent argument of knowledge with
//! a coset low-degree extension, a composition polynomial that
//! encodes the AIR, and degree-bounded FRI.  Educational and not
//! audited: it must not back a research claim.
//!
//! 1. **Coset LDE.**  The trace polynomial `f(X)`, interpolated on
//!    `H_trace = ⟨ω⟩` of size `n`, is evaluated on the coset
//!    `g · H_lde` of size `L = n · β` (`β` = blow-up, a power of two
//!    `≥ 2`; `g = 7` generates `F_r*`, so the coset misses `H_trace`).
//!    Because `ω = ω_L^β`, `f(ω·x)` sits `β` positions after `f(x)`.
//!
//! 2. **Composition polynomial**, evaluated pointwise on the coset:
//!
//!    ```text
//!    Comp(x) = α₁ · C(f(x), f(ω·x)) / Z(x)  +  Σ_k α₂^k · (f(x) − v_k) / (x − ω^{p_k})
//!    Z(x)    = (xⁿ − 1) / (x − ω^{n−1})
//!    ```
//!
//!    For a valid trace and a transition `C` of degree 1 in the trace
//!    values, `Comp` is a polynomial of degree `< n`; if a transition
//!    or a boundary fails, the division leaves a non-polynomial and
//!    `Comp` is far from every degree-`< n` polynomial.
//!
//! 3. **Degree-bounded FRI** folds `Comp`'s evaluations `log₂ n` times
//!    to a final layer of `β` values, which must be constant.
//!
//! The verifier re-derives `α₁, α₂`, the FRI challenges and the query
//! positions from the transcript; at each of `query_count` positions
//! it opens `f(x)`, `f(ω·x)` and `Comp(x)` (each Merkle path bound to
//! its position), recomputes `Comp(x)` from the AIR, checks that FRI's
//! first layer holds that value, and checks the folding down to the
//! constant final layer.
//!
//! **Soundness.**  Under the Reed–Solomon proximity conjecture each
//! query catches a cheating prover with probability about
//! `1 − 1/β`, so `q` queries give about `q · log₂ β` bits (240 for
//! the default `β = 8, q = 80`; 40 for [`StarkConfig::fast`]).  Proven
//! bounds are weaker (about half that), there is no grinding, no
//! out-of-domain (DEEP) sample and no zero-knowledge masking, the
//! transition must be of degree 1, and SHA-256 and the 255-bit field
//! cap everything.
//!
//! ## API
//!
//! - [`AirV2`]: boundary constraints and a closure-based transition
//!   on adjacent trace rows.
//! - [`StarkConfig`]: `(blow_up, query_count)`.
//! - [`prove`] / [`verify`].

use super::merkle::{merkle_verify, MerkleProof, MerkleTree};
use super::polynomial::{fr_add, fr_inv, fr_mul, fr_sub, Poly};
use super::stark::{
    fri_commit_to, fri_query, fri_transcript_after, fri_verify, query_positions, FriCommitment,
    FriQueryOpening,
};
use crate::bls12_381::fq::scalar_modulus;
use crate::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::One;
use std::sync::Arc;

const FS_TAG: &str = "ZK-STARK-V2/v1";

/// Configuration: blow-up factor and number of FRI queries.
#[derive(Clone, Debug)]
pub struct StarkConfig {
    pub blow_up: usize,
    pub query_count: usize,
}

impl StarkConfig {
    /// Default: `blow_up = 8`, `query_count = 80`: about 240 bits
    /// under the proximity conjecture (module docs, "Soundness").
    #[allow(clippy::should_implement_trait)]
    pub fn default() -> Self {
        Self {
            blow_up: 8,
            query_count: 80,
        }
    }
    /// Faster: `blow_up = 4`, `query_count = 20`, about 40 bits
    /// conjectured.  Good for tests.
    pub fn fast() -> Self {
        Self {
            blow_up: 4,
            query_count: 20,
        }
    }
}

/// An AIR with boundary constraints and a transition operating on
/// adjacent trace rows (`f(X), f(ω·X)`).
#[derive(Clone)]
pub struct AirV2 {
    /// Number of trace rows.  Must be a power of 2.
    pub n: usize,
    /// Boundary constraints: `(position, expected_value)` pairs.
    /// Each `position` is in `[0, n)`.
    pub boundary: Vec<(usize, BigUint)>,
    /// Transition: given two consecutive trace values, return the
    /// constraint value (should be zero on valid traces).
    pub transition: Arc<dyn Fn(&BigUint, &BigUint) -> BigUint + Send + Sync>,
    /// Generator of the trace (the prover's actual computation).
    pub generate_trace: Arc<dyn Fn() -> Vec<BigUint> + Send + Sync>,
}

impl AirV2 {
    /// **Doubling AIR**: trace[i+1] = 2·trace[i], starting from `start`.
    /// Boundary: `trace[0] = start`.
    pub fn doubling(n: usize, start: u64) -> Self {
        let s = start;
        let s_for_gen = s;
        let gen = move || -> Vec<BigUint> {
            let mut t = Vec::with_capacity(n);
            t.push(BigUint::from(s_for_gen));
            for i in 1..n {
                t.push(fr_mul(&t[i - 1], &BigUint::from(2u32)));
            }
            t
        };
        let transition = Arc::new(|prev: &BigUint, next: &BigUint| -> BigUint {
            // Constraint: next − 2·prev = 0.
            fr_sub(next, &fr_mul(prev, &BigUint::from(2u32)))
        });
        let generate_trace = Arc::new(gen);
        let boundary = vec![(0, BigUint::from(s))];
        Self {
            n,
            boundary,
            transition,
            generate_trace,
        }
    }
}

// ── Domain helpers ────────────────────────────────────────────────

/// Build the multiplicative subgroup of size `n` (`n` must be a
/// power of 2 dividing 2³²).
fn subgroup(n: usize) -> Vec<BigUint> {
    let omega = super::plonk::primitive_nth_root(n);
    let mut out = Vec::with_capacity(n);
    let mut acc = BigUint::one();
    for _ in 0..n {
        out.push(acc.clone());
        acc = fr_mul(&acc, &omega);
    }
    out
}

/// `g · H` where `g` is a non-trivial coset shift.  We use
/// `g = 7` (multiplicative generator of F_r*) to ensure the coset
/// is disjoint from the trace domain.
fn coset(n: usize, shift: &BigUint) -> Vec<BigUint> {
    subgroup(n).into_iter().map(|x| fr_mul(&x, shift)).collect()
}

fn coset_shift() -> BigUint {
    BigUint::from(7u32)
}

// ── Polynomial evaluation on subgroups ─────────────────────────────

/// Naive O(n²) DFT: evaluate polynomial `p` on the domain `dom`.
fn eval_on_domain(p: &Poly, dom: &[BigUint]) -> Vec<BigUint> {
    dom.iter().map(|x| p.evaluate(x)).collect()
}

/// Naive interpolation: given evaluations of a polynomial of degree
/// `< n` on a multiplicative subgroup of size `n`, recover the
/// polynomial.  Uses Lagrange via [`Poly::interpolate`].
fn interpolate_from_subgroup(dom: &[BigUint], evals: &[BigUint]) -> Poly {
    let pts: Vec<(BigUint, BigUint)> = dom
        .iter()
        .zip(evals)
        .map(|(d, e)| (d.clone(), e.clone()))
        .collect();
    Poly::interpolate(&pts)
}

// ── Composition polynomial ─────────────────────────────────────────

/// `Comp(x)` from the trace values `f(x)`, `f(ω·x)` at one coset point
/// (module docs, item 2).  Prover and verifier call the same function.
fn composition_at(
    air: &AirV2,
    h_trace: &[BigUint],
    x: &BigUint,
    fx: &BigUint,
    f_wx: &BigUint,
    alpha1: &BigUint,
    alpha2: &BigUint,
) -> BigUint {
    let n = air.n;
    // xⁿ by squaring (n is a power of two).
    let mut xn = x.clone();
    let mut k = 1;
    while k < n {
        xn = fr_mul(&xn, &xn);
        k *= 2;
    }
    // Z(x) = (xⁿ − 1) / (x − ω^{n−1}); the coset never meets H_trace.
    let z = fr_mul(
        &fr_sub(&xn, &BigUint::one()),
        &fr_inv(&fr_sub(x, &h_trace[n - 1])).expect("coset point off the trace domain"),
    );
    let transition = fr_mul(
        &(air.transition)(fx, f_wx),
        &fr_inv(&z).expect("coset point off the trace domain"),
    );
    let mut acc = fr_mul(alpha1, &transition);
    let mut alpha2_pow = BigUint::one();
    for (pos, val) in &air.boundary {
        let term = fr_mul(
            &fr_sub(fx, val),
            &fr_inv(&fr_sub(x, &h_trace[*pos])).expect("coset point off the trace domain"),
        );
        acc = fr_add(&acc, &fr_mul(&alpha2_pow, &term));
        alpha2_pow = fr_mul(&alpha2_pow, alpha2);
    }
    acc
}

/// Fiat-Shamir state up to FRI: absorbs the trace commitment, draws
/// `α₁, α₂`, absorbs the composition commitment.
fn transcript_to_fri(
    trace_commit: &[u8; 32],
    composition_commit: &[u8; 32],
) -> (Vec<u8>, BigUint, BigUint) {
    let mut transcript = FS_TAG.as_bytes().to_vec();
    transcript.extend_from_slice(&trace_commit[..]);
    let alpha1_bytes = sha256(&transcript);
    transcript.extend_from_slice(&alpha1_bytes);
    let alpha1 = BigUint::from_bytes_be(&alpha1_bytes) % scalar_modulus();
    let alpha2_bytes = sha256(&transcript);
    transcript.extend_from_slice(&alpha2_bytes);
    let alpha2 = BigUint::from_bytes_be(&alpha2_bytes) % scalar_modulus();
    transcript.extend_from_slice(&composition_commit[..]);
    (transcript, alpha1, alpha2)
}

// ── Prover / Verifier ──────────────────────────────────────────────

#[derive(Clone, Debug)]
pub struct StarkProofV2 {
    pub trace_commit: [u8; 32],
    pub composition_commit: [u8; 32],
    pub fri: FriCommitment,
    /// Per-query openings on the trace and composition LDE tables:
    /// `f(x)`, `f(ω·x)` (position `+ β`) and `Comp(x)`.
    pub trace_openings: Vec<(usize, BigUint, MerkleProof)>,
    pub trace_next_openings: Vec<(usize, BigUint, MerkleProof)>,
    pub composition_openings: Vec<(usize, BigUint, MerkleProof)>,
    pub fri_openings: Vec<(usize, FriQueryOpening)>,
}

fn scalar_to_bytes(s: &BigUint) -> Vec<u8> {
    let mut buf = [0u8; 32];
    let bytes = s.to_bytes_be();
    buf[32 - bytes.len()..].copy_from_slice(&bytes);
    buf.to_vec()
}

/// **Prove** the AIR is satisfied by an executed trace.  The prover
/// does not check the trace: an invalid one yields a proof that
/// [`verify`] rejects.
pub fn prove(air: &AirV2, config: &StarkConfig) -> Option<StarkProofV2> {
    let n = air.n;
    let beta = config.blow_up;
    if !n.is_power_of_two() || n < 2 || !beta.is_power_of_two() || beta < 2 {
        return None;
    }
    let trace = (air.generate_trace)();
    if trace.len() != n {
        return None;
    }
    let h_trace = subgroup(n);

    // Trace polynomial via Lagrange, evaluated on the coset g · H_lde.
    let f_poly = interpolate_from_subgroup(&h_trace, &trace);
    let lde_size = n * beta;
    let h_lde = coset(lde_size, &coset_shift());
    let f_lde = eval_on_domain(&f_poly, &h_lde);

    let trace_leaves: Vec<Vec<u8>> = f_lde.iter().map(scalar_to_bytes).collect();
    let trace_tree = MerkleTree::from_leaves(&trace_leaves);
    let trace_commit = trace_tree.root();

    // α₁, α₂ depend only on the trace commitment.
    let (_, alpha1, alpha2) = transcript_to_fri(&trace_commit, &[0u8; 32]);
    let composition_lde: Vec<BigUint> = (0..lde_size)
        .map(|i| {
            composition_at(
                air,
                &h_trace,
                &h_lde[i],
                &f_lde[i],
                &f_lde[(i + beta) % lde_size],
                &alpha1,
                &alpha2,
            )
        })
        .collect();
    let comp_leaves: Vec<Vec<u8>> = composition_lde.iter().map(scalar_to_bytes).collect();
    let comp_tree = MerkleTree::from_leaves(&comp_leaves);
    let composition_commit = comp_tree.root();
    let (mut transcript, _, _) = transcript_to_fri(&trace_commit, &composition_commit);

    // FRI to a constant final layer of β values: degree < n.
    let fri = fri_commit_to(&composition_lde, &h_lde, &mut transcript, beta);
    let positions = query_positions(&mut transcript, config.query_count, lde_size);

    let open = |tree: &MerkleTree, table: &[BigUint], p: usize| {
        (p, table[p].clone(), tree.proof(p).unwrap())
    };
    Some(StarkProofV2 {
        trace_commit,
        composition_commit,
        trace_openings: positions
            .iter()
            .map(|&p| open(&trace_tree, &f_lde, p))
            .collect(),
        trace_next_openings: positions
            .iter()
            .map(|&p| open(&trace_tree, &f_lde, (p + beta) % lde_size))
            .collect(),
        composition_openings: positions
            .iter()
            .map(|&p| open(&comp_tree, &composition_lde, p))
            .collect(),
        fri_openings: positions.iter().map(|&p| (p, fri_query(&fri, p))).collect(),
        fri,
    })
}

/// **Verify** a STARK v2 proof (module docs, last paragraph).
pub fn verify(air: &AirV2, config: &StarkConfig, proof: &StarkProofV2) -> bool {
    let n = air.n;
    let beta = config.blow_up;
    let q = config.query_count;
    if !n.is_power_of_two() || n < 2 || !beta.is_power_of_two() || beta < 2 || q == 0 {
        return false;
    }
    if air.boundary.iter().any(|(pos, _)| *pos >= n) {
        return false;
    }
    if proof.trace_openings.len() != q
        || proof.trace_next_openings.len() != q
        || proof.composition_openings.len() != q
        || proof.fri_openings.len() != q
    {
        return false;
    }
    // log₂ n folds: the final layer of β values must be constant.
    if proof.fri.layer_roots.len() != n.trailing_zeros() as usize + 1 {
        return false;
    }
    let lde_size = n * beta;
    let h_lde = coset(lde_size, &coset_shift());
    let h_trace = subgroup(n);

    let (fri_transcript_initial, alpha1, alpha2) =
        transcript_to_fri(&proof.trace_commit, &proof.composition_commit);
    let mut after_fri = fri_transcript_after(&fri_transcript_initial, &proof.fri);
    let positions = query_positions(&mut after_fri, q, lde_size);

    let opened =
        |(pos, val, path): &(usize, BigUint, MerkleProof), at: usize, root: &[u8; 32]| -> bool {
            *pos == at && path.leaf_index == at && merkle_verify(&scalar_to_bytes(val), path, root)
        };
    for (i, &p) in positions.iter().enumerate() {
        let (fx, f_wx, comp) = (
            &proof.trace_openings[i],
            &proof.trace_next_openings[i],
            &proof.composition_openings[i],
        );
        if !opened(fx, p, &proof.trace_commit)
            || !opened(f_wx, (p + beta) % lde_size, &proof.trace_commit)
            || !opened(comp, p, &proof.composition_commit)
        {
            return false;
        }
        // The committed composition is the AIR's composition here ...
        if comp.1 != composition_at(air, &h_trace, &h_lde[p], &fx.1, &f_wx.1, &alpha1, &alpha2) {
            return false;
        }
        // ... and FRI's first layer is the committed composition.
        let (fri_pos, fri_open) = &proof.fri_openings[i];
        if *fri_pos != p || fri_open.openings.first().map(|o| &o.0) != Some(&comp.1) {
            return false;
        }
    }
    if proof.fri.layer_roots[0] != proof.composition_commit {
        return false;
    }
    fri_verify(
        &proof.fri,
        &h_lde,
        &proof.fri_openings,
        &fri_transcript_initial,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Zero;

    /// AIR doubling generates the expected trace.
    #[test]
    fn air_doubling_generates_powers_of_two() {
        let air = AirV2::doubling(8, 3);
        let trace = (air.generate_trace)();
        let expected: Vec<BigUint> = (0..8).map(|i| BigUint::from(3 * (1u64 << i))).collect();
        assert_eq!(trace, expected);
    }

    /// AIR transition constraint zeroes on valid adjacent pairs.
    #[test]
    fn air_transition_zero_on_valid_pair() {
        let air = AirV2::doubling(8, 3);
        // Valid pair: (3, 6).  3·2 = 6.
        let c = (air.transition)(&BigUint::from(3u32), &BigUint::from(6u32));
        assert!(c.is_zero());
    }

    /// AIR transition is non-zero on invalid pair.
    #[test]
    fn air_transition_nonzero_on_invalid_pair() {
        let air = AirV2::doubling(8, 3);
        // Invalid pair: (3, 7).
        let c = (air.transition)(&BigUint::from(3u32), &BigUint::from(7u32));
        assert!(!c.is_zero());
    }

    /// **End-to-end STARK v2** with fast config (blow_up=4, queries=20).
    #[test]
    fn stark_v2_prove_verify_doubling_fast() {
        let air = AirV2::doubling(8, 1);
        let cfg = StarkConfig::fast();
        let proof = prove(&air, &cfg).expect("proving");
        assert!(
            verify(&air, &cfg, &proof),
            "v2 verify must accept honest proof"
        );
    }

    /// Default config (blow_up=8, queries=80) on a longer trace.
    #[test]
    fn stark_v2_prove_verify_doubling_default() {
        let air = AirV2::doubling(16, 3);
        let cfg = StarkConfig::default();
        let proof = prove(&air, &cfg).expect("proving");
        assert!(verify(&air, &cfg, &proof));
    }

    // Negative controls: each is a proof of a FALSE statement or a
    // mangled proof, and the verifier must reject it.  The first two
    // are the mutants a 2026-10-16 sheriff audit found accepted
    // (crypto-autoresearcher research/sheriff/SHF-20261016-3caa6f).

    /// An honest proof (trace starts at 1) for another boundary (5).
    #[test]
    fn rejects_proof_reused_for_another_boundary() {
        let cfg = StarkConfig::fast();
        let proof = prove(&AirV2::doubling(8, 1), &cfg).unwrap();
        assert!(!verify(&AirV2::doubling(8, 5), &cfg, &proof));
    }

    /// A proof with no queries proves nothing.
    #[test]
    fn rejects_empty_proof() {
        let air = AirV2::doubling(8, 1);
        let cfg = StarkConfig::fast();
        let mut proof = prove(&air, &cfg).unwrap();
        proof.trace_openings.clear();
        proof.trace_next_openings.clear();
        proof.composition_openings.clear();
        proof.fri_openings.clear();
        assert!(!verify(&air, &cfg, &proof));
    }

    /// A trace that breaks one transition: the composition is not a
    /// polynomial of degree < n, and degree-bounded FRI rejects it.
    #[test]
    fn rejects_trace_with_broken_transition() {
        let cfg = StarkConfig::fast();
        let mut cheat = AirV2::doubling(8, 1);
        cheat.generate_trace = Arc::new(|| {
            let mut t: Vec<BigUint> = (0..8).map(|i| BigUint::from(1u64 << i)).collect();
            t[5] = BigUint::from(33u32);
            t
        });
        let forged = prove(&cheat, &cfg).unwrap();
        assert!(!verify(&AirV2::doubling(8, 1), &cfg, &forged));
    }

    /// Dropping a FRI layer drops the degree bound.
    #[test]
    fn rejects_truncated_fri() {
        let air = AirV2::doubling(8, 1);
        let cfg = StarkConfig::fast();
        let mut proof = prove(&air, &cfg).unwrap();
        proof.fri.layer_roots.pop();
        assert!(!verify(&air, &cfg, &proof));
    }

    /// Queries must come in the order the transcript dictates.
    #[test]
    fn rejects_swapped_queries() {
        let air = AirV2::doubling(8, 1);
        let cfg = StarkConfig::fast();
        let mut proof = prove(&air, &cfg).unwrap();
        let (a, b) = (
            0,
            (1..cfg.query_count)
                .find(|&j| proof.trace_openings[j].0 != proof.trace_openings[0].0)
                .expect("two distinct positions"),
        );
        proof.trace_openings.swap(a, b);
        proof.trace_next_openings.swap(a, b);
        proof.composition_openings.swap(a, b);
        proof.fri_openings.swap(a, b);
        assert!(!verify(&air, &cfg, &proof));
    }
}
