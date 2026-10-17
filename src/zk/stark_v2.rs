//! **STARK v2**: a scalable transparent argument of knowledge with
//! a coset low-degree extension, a composition polynomial that
//! encodes the AIR, and degree-bounded FRI.
//!
//! **Educational and unaudited.  Not zero-knowledge (no masking: the
//! opened trace-LDE values leak information about the witness), no
//! out-of-domain (DEEP) sample, no grinding, no proven soundness
//! bound.  It must not back a research claim.**
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
//!    `Comp` is far from every degree-`< n` polynomial — unless the
//!    poles cancel, which is why `α₁, α₂` must be drawn only after the
//!    statement (the `v_k`) is fixed.
//!
//! 3. **Degree-bounded FRI** folds `Comp`'s evaluations `log₂ n` times
//!    to a final layer of `β` values, which must be constant.
//!
//! **Transcript (Fiat-Shamir).**  It absorbs, in order: a
//! domain-separation tag; the whole statement — an AIR identifier
//! ([`air_fingerprint`]), `n`, `β`, `q`, the FRI fold count `log₂ n`
//! and final-layer size `β`, and every boundary `(p_k, v_k)`
//! ([`statement_transcript`]); the trace commitment; then draws `α₁`,
//! `α₂`; absorbs the composition commitment; each FRI layer root
//! before its folding challenge; the last layer root and the final
//! constant; and only then draws the `q` query positions.  A prover
//! therefore cannot pick a public input after seeing `α₁, α₂` (the
//! "frozen heart" forgery of the 2026-10-17 sheriff re-audit).
//!
//! The verifier takes `β`, `q` and the FRI depth from its own
//! [`StarkConfig`], never from the proof, and rejects a proof whose
//! query count, layer count or final layer differs.  At each position
//! it opens `f(x)`, `f(ω·x)` and `Comp(x)` (each Merkle path bound to
//! its position, each value a canonical field element `< r`),
//! recomputes `Comp(x)` from the AIR, checks that FRI's first layer
//! holds that value, and checks the folding down to the committed
//! constant final layer.  Malformed input returns `false`; it never
//! panics.
//!
//! **Soundness: no security level is claimed for this code.**  The
//! usual *conjectured* analysis of a FRI-based STARK assumes the
//! Reed–Solomon proximity conjecture and Fiat–Shamir in the random
//! oracle model, and is normally stated for DEEP-ALI constructions,
//! which this is not.  It bounds the soundness error of one proof by
//! roughly
//!
//! ```text
//! ε  ≲  ρ^q · 2^(−g)  +  ε_field  +  ε_hash,      ρ = 1/β (code rate)
//! ```
//!
//! with `q` queries and `g` grinding bits.  Here `g = 0`, so the
//! conjectured query term is `β^(−q)` (`8^(−80)` for
//! [`StarkConfig::default`], `4^(−20)` for [`StarkConfig::fast`]).
//! `ε_field` covers unlucky draws of `α₁`, `α₂` and the `log₂ n`
//! folding challenges from `F_r` (BLS12-381 scalar field, `r ≈
//! 2^254.86`); it depends on `L` and `r` and is not analysed here.
//! `ε_hash`: every commitment is a SHA-256 Merkle tree, and a SHA-256
//! collision (generic cost about `2^128`) lets a prover open one leaf
//! two ways, so nothing above 128 bits can hold whatever `q` and `β`
//! are; a prover making `T` hash queries also gains a factor of about
//! `T` through Fiat–Shamir.  Proven (non-conjectural) FRI bounds are
//! weaker than the conjectured query term; this module does not state
//! one.  The transition must be of degree 1 in the trace values.
//!
//! ## API
//!
//! - [`AirV2`]: boundary constraints and a closure-based transition
//!   on adjacent trace rows.
//! - [`StarkConfig`]: `(blow_up, query_count)`.
//! - [`prove`] / [`verify`].

use super::merkle::{MerkleProof, MerkleTree};
use super::polynomial::{fr_add, fr_inv, fr_mul, fr_sub, Poly};
use super::stark::{
    absorb_scalar, absorb_u64, draw_challenge, fri_commit_to, fri_query, fri_transcript_after,
    fri_verify, opening_ok, probe_point, query_positions, reduced_bytes, FriCommitment,
    FriQueryOpening,
};
use crate::bls12_381::fq::scalar_modulus;
use num_bigint::BigUint;
use num_traits::One;
use std::sync::Arc;

/// Transcript domain-separation tag.  The `statement-bound` suffix marks
/// the transcript that absorbs the statement before any challenge.
const FS_TAG: &str = "ZK-STARK-V2/v2/statement-bound";

/// Number of fixed pseudo-random points (beyond `(0,0)`, `(1,0)`,
/// `(0,1)`) at which [`air_fingerprint`] evaluates a transition.
const AIR_PROBES: u64 = 2;

/// Configuration: blow-up factor and number of FRI queries.  The
/// verifier's copy is authoritative: proofs made with another
/// configuration are rejected.
#[derive(Clone, Debug)]
pub struct StarkConfig {
    pub blow_up: usize,
    pub query_count: usize,
}

impl StarkConfig {
    /// Default: `blow_up = 8`, `query_count = 80`.  Conjectured query
    /// term `8^(−80)` (module docs, "Soundness"); that is not a
    /// security level for this code, which claims none.
    #[allow(clippy::should_implement_trait)]
    pub fn default() -> Self {
        Self {
            blow_up: 8,
            query_count: 80,
        }
    }
    /// Faster: `blow_up = 4`, `query_count = 20`, conjectured query
    /// term `4^(−20)`.  Good for tests.
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

/// **AIR identifier**: the transition evaluated (mod `r`) at `(0,0)`,
/// `(1,0)`, `(0,1)` and [`AIR_PROBES`] fixed public points
/// `H("ZK-STARK-V2/air-probe" ‖ j)`.  A transition of degree `≤ 1`
/// (the only kind this module supports) is determined exactly by its
/// first three values; for higher degree this is a fingerprint.
pub fn air_fingerprint(air: &AirV2) -> Vec<u8> {
    let (zero, one) = (BigUint::from(0u32), BigUint::one());
    let mut points = vec![
        (zero.clone(), zero.clone()),
        (one.clone(), zero.clone()),
        (zero, one),
    ];
    points.extend((0..AIR_PROBES).map(|j| {
        (
            probe_point("ZK-STARK-V2/air-probe", 2 * j),
            probe_point("ZK-STARK-V2/air-probe", 2 * j + 1),
        )
    }));
    let mut out = b"air:v2-adjacent-row-transition".to_vec();
    for (a, b) in &points {
        out.extend_from_slice(&reduced_bytes(&(air.transition)(a, b)));
    }
    out
}

/// The Fiat-Shamir transcript state before the trace commitment: tag,
/// AIR identifier, `n`, `β`, `q`, FRI fold count and final-layer size,
/// and every boundary `(position, value)`.  `None` if the configuration
/// is invalid, a boundary position is `≥ n`, or a boundary value is not
/// canonical (`≥ r`).
pub fn statement_transcript(air: &AirV2, config: &StarkConfig) -> Option<Vec<u8>> {
    let (n, beta, q) = (air.n, config.blow_up, config.query_count);
    if !n.is_power_of_two() || n < 2 || !beta.is_power_of_two() || beta < 2 || q == 0 {
        return None;
    }
    let mut t = FS_TAG.as_bytes().to_vec();
    t.extend_from_slice(&air_fingerprint(air));
    absorb_u64(&mut t, b"n", n as u64);
    absorb_u64(&mut t, b"blow_up", beta as u64);
    absorb_u64(&mut t, b"q", q as u64);
    absorb_u64(&mut t, b"fri_folds", n.trailing_zeros() as u64);
    absorb_u64(&mut t, b"fri_final_size", beta as u64);
    absorb_u64(&mut t, b"boundaries", air.boundary.len() as u64);
    for (pos, val) in &air.boundary {
        if *pos >= n {
            return None;
        }
        absorb_u64(&mut t, b"p", *pos as u64);
        if !absorb_scalar(&mut t, val) {
            return None;
        }
    }
    t.extend_from_slice(b"trace_commit");
    Some(t)
}

/// Absorb the trace commitment and draw `(α₁, α₂)`.
fn absorb_trace_draw_alphas(t: &mut Vec<u8>, trace_commit: &[u8; 32]) -> (BigUint, BigUint) {
    t.extend_from_slice(trace_commit);
    let alpha1 = draw_challenge(t);
    let alpha2 = draw_challenge(t);
    t.extend_from_slice(b"composition_commit");
    (alpha1, alpha2)
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

/// Deviations from the honest prover, used by the tests as cheating
/// provers.  [`prove`] uses the default (no deviation).
#[derive(Clone, Debug, Default)]
struct ProverKnobs {
    /// Fold FRI to this many values instead of `β`.
    fri_final_size: Option<usize>,
    /// Build the composition with `α₁ = 0` (drops the transition).
    alpha1_zero: bool,
    /// Commit the zero polynomial as the composition.
    composition_zero: bool,
}

/// **Prove** the AIR is satisfied by an executed trace.  The prover
/// does not check the trace: an invalid one yields a proof that
/// [`verify`] rejects.  `None` if the configuration or the statement is
/// invalid (see [`statement_transcript`]) or the trace has the wrong
/// length; trace values are reduced mod `r`.
pub fn prove(air: &AirV2, config: &StarkConfig) -> Option<StarkProofV2> {
    prove_with(air, config, &ProverKnobs::default())
}

fn prove_with(air: &AirV2, config: &StarkConfig, knobs: &ProverKnobs) -> Option<StarkProofV2> {
    let mut transcript = statement_transcript(air, config)?;
    let n = air.n;
    let beta = config.blow_up;
    let r = scalar_modulus();
    let trace: Vec<BigUint> = (air.generate_trace)().iter().map(|v| v % &r).collect();
    if trace.len() != n {
        return None;
    }
    let h_trace = subgroup(n);

    // Trace polynomial via Lagrange, evaluated on the coset g · H_lde.
    let f_poly = interpolate_from_subgroup(&h_trace, &trace);
    let lde_size = n * beta;
    let h_lde = coset(lde_size, &coset_shift());
    let f_lde: Vec<BigUint> = eval_on_domain(&f_poly, &h_lde)
        .iter()
        .map(|v| v % &r)
        .collect();

    let trace_leaves: Vec<Vec<u8>> = f_lde.iter().map(reduced_bytes).collect();
    let trace_tree = MerkleTree::from_leaves(&trace_leaves);
    let trace_commit = trace_tree.root();

    // α₁, α₂ depend on the statement and the trace commitment.
    let (alpha1, alpha2) = absorb_trace_draw_alphas(&mut transcript, &trace_commit);
    let alpha1 = if knobs.alpha1_zero {
        BigUint::from(0u32)
    } else {
        alpha1
    };
    let composition_lde: Vec<BigUint> = (0..lde_size)
        .map(|i| {
            if knobs.composition_zero {
                return BigUint::from(0u32);
            }
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
    let comp_leaves: Vec<Vec<u8>> = composition_lde.iter().map(reduced_bytes).collect();
    let comp_tree = MerkleTree::from_leaves(&comp_leaves);
    let composition_commit = comp_tree.root();
    transcript.extend_from_slice(&composition_commit);

    // FRI to a constant final layer of β values: degree < n.
    let final_size = knobs.fri_final_size.unwrap_or(beta);
    if !final_size.is_power_of_two() || final_size > lde_size {
        return None;
    }
    let fri = fri_commit_to(&composition_lde, &h_lde, &mut transcript, final_size);
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

/// **Verify** a STARK v2 proof (module docs, "Transcript" and the
/// paragraph after it) against the verifier's own `config`.
pub fn verify(air: &AirV2, config: &StarkConfig, proof: &StarkProofV2) -> bool {
    // Checks n, β, q, boundary positions and canonical boundary values.
    let Some(mut transcript) = statement_transcript(air, config) else {
        return false;
    };
    let n = air.n;
    let beta = config.blow_up;
    let q = config.query_count;
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
    if proof.fri.layer_roots[0] != proof.composition_commit {
        return false;
    }
    let lde_size = n * beta;
    let depth = lde_size.trailing_zeros() as usize;
    let h_lde = coset(lde_size, &coset_shift());
    let h_trace = subgroup(n);

    let (alpha1, alpha2) = absorb_trace_draw_alphas(&mut transcript, &proof.trace_commit);
    transcript.extend_from_slice(&proof.composition_commit);
    let fri_transcript_initial = transcript.clone();
    let Some(mut after_fri) = fri_transcript_after(&fri_transcript_initial, &proof.fri) else {
        return false;
    };
    let positions = query_positions(&mut after_fri, q, lde_size);

    for (i, &p) in positions.iter().enumerate() {
        let (fx, f_wx, comp) = (
            &proof.trace_openings[i],
            &proof.trace_next_openings[i],
            &proof.composition_openings[i],
        );
        if !opening_ok(fx, p, depth, &proof.trace_commit)
            || !opening_ok(f_wx, (p + beta) % lde_size, depth, &proof.trace_commit)
            || !opening_ok(comp, p, depth, &proof.composition_commit)
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

    type T2 = Arc<dyn Fn(&BigUint, &BigUint) -> BigUint + Send + Sync>;

    fn big(x: u64) -> BigUint {
        BigUint::from(x)
    }
    fn dbl() -> T2 {
        Arc::new(|p: &BigUint, nx: &BigUint| fr_sub(nx, &fr_mul(p, &big(2))))
    }
    fn trip() -> T2 {
        Arc::new(|p: &BigUint, nx: &BigUint| fr_sub(nx, &fr_mul(p, &big(3))))
    }
    fn air2(n: usize, boundary: Vec<(usize, BigUint)>, t: T2, trace: Vec<BigUint>) -> AirV2 {
        AirV2 {
            n,
            boundary,
            transition: t,
            generate_trace: Arc::new(move || trace.clone()),
        }
    }
    fn doubling_trace(n: usize, start: u64) -> Vec<BigUint> {
        (0..n).map(|i| big(start << i)).collect()
    }
    /// Doubling AIR with boundaries (0, 1) and (7, 128), honest trace.
    fn two_boundary_air() -> AirV2 {
        air2(
            8,
            vec![(0, big(1)), (7, big(128))],
            dbl(),
            doubling_trace(8, 1),
        )
    }
    /// The false statement "1 doubled seven times is 12345", proved from
    /// the trace [1, 2, .., 64, 12345] (one broken transition).
    fn false_output_air() -> AirV2 {
        let mut t = doubling_trace(8, 1);
        t[7] = big(12345);
        air2(8, vec![(0, big(1)), (7, big(12345))], dbl(), t)
    }

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
        let air2b = two_boundary_air();
        let proof = prove(&air2b, &cfg).unwrap();
        assert!(verify(&air2b, &cfg, &proof));
    }

    /// Default config (blow_up=8, queries=80), n = 8 and a longer trace.
    #[test]
    fn stark_v2_prove_verify_doubling_default() {
        let cfg = StarkConfig::default();
        for air in [AirV2::doubling(8, 1), AirV2::doubling(16, 3)] {
            let proof = prove(&air, &cfg).expect("proving");
            assert!(verify(&air, &cfg, &proof));
        }
    }

    // Negative controls: each is a proof of a FALSE statement or a
    // mangled proof, and the verifier must reject it.  The first two
    // are the mutants a 2026-10-16 sheriff audit found accepted
    // (crypto-autoresearcher research/sheriff/SHF-20261016-3caa6f); the
    // rest port the re-audit's mutants (SHF-20261017-e4d7f6, reader
    // "stark", new-mutants/src/main.rs, rows `v2.*`).

    /// An honest proof (trace starts at 1) for another boundary (5);
    /// N4: a valid proof of doubling(8, 3) for doubling(8, 1).
    #[test]
    fn rejects_proof_reused_for_another_boundary() {
        let cfg = StarkConfig::fast();
        let proof = prove(&AirV2::doubling(8, 1), &cfg).unwrap();
        assert!(!verify(&AirV2::doubling(8, 5), &cfg, &proof));
        let p3 = prove(&AirV2::doubling(8, 3), &cfg).unwrap();
        assert!(!verify(&AirV2::doubling(8, 1), &cfg, &p3));
    }

    /// N2a / N2b: wrong public output or input.
    #[test]
    fn rejects_wrong_public_values() {
        let cfg = StarkConfig::fast();
        let proof = prove(&two_boundary_air(), &cfg).unwrap();
        let t = doubling_trace(8, 1);
        let out = air2(8, vec![(0, big(1)), (7, big(129))], dbl(), t.clone());
        assert!(!verify(&out, &cfg, &proof));
        let inp = air2(8, vec![(0, big(3)), (7, big(128))], dbl(), t.clone());
        assert!(!verify(&inp, &cfg, &proof));
        // Dropping a boundary changes the statement too.
        let fewer = air2(8, vec![(0, big(1))], dbl(), t);
        assert!(!verify(&fewer, &cfg, &proof));
    }

    /// N3a / N3b: wrong AIR (tripling: 3^7 != 128, a false statement;
    /// and a 16-row AIR).
    #[test]
    fn rejects_wrong_air() {
        let cfg = StarkConfig::fast();
        let proof = prove(&two_boundary_air(), &cfg).unwrap();
        let tripling = air2(
            8,
            vec![(0, big(1)), (7, big(128))],
            trip(),
            doubling_trace(8, 1),
        );
        assert_ne!(
            air_fingerprint(&tripling),
            air_fingerprint(&two_boundary_air())
        );
        assert!(!verify(&tripling, &cfg, &proof));
        let p8 = prove(&AirV2::doubling(8, 1), &cfg).unwrap();
        assert!(!verify(&AirV2::doubling(16, 1), &cfg, &p8));
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
    /// Also N1a (tampered row, true statement) and N1b (false output).
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
        let mut t = doubling_trace(8, 1);
        t[3] = big(99);
        let tampered = prove(&air2(8, vec![(0, big(1))], dbl(), t), &cfg).unwrap();
        assert!(!verify(&AirV2::doubling(8, 1), &cfg, &tampered));
        let fa = false_output_air();
        assert!(!verify(&fa, &cfg, &prove(&fa, &cfg).unwrap()));
    }

    /// Dropping a FRI layer drops the degree bound (N5a); so does
    /// dropping each query's last opening (N5b) or folding once fewer
    /// with the transcript re-derived (N5c, true and false statements).
    #[test]
    fn rejects_truncated_fri() {
        let air = AirV2::doubling(8, 1);
        let cfg = StarkConfig::fast();
        let mut proof = prove(&air, &cfg).unwrap();
        proof.fri.layer_roots.pop();
        assert!(!verify(&air, &cfg, &proof));
        let mut proof = prove(&air, &cfg).unwrap();
        for fq in proof.fri_openings.iter_mut() {
            fq.1.openings.pop();
        }
        assert!(!verify(&air, &cfg, &proof));
        let fewer = ProverKnobs {
            fri_final_size: Some(2 * cfg.blow_up),
            ..ProverKnobs::default()
        };
        for a in [air, false_output_air()] {
            let p = prove_with(&a, &cfg, &fewer).unwrap();
            assert!(!verify(&a, &cfg, &p));
        }
    }

    /// N5d / M1: the final constant and the last FRI root are bound
    /// (zeroing the last root was accepted before: it was never
    /// absorbed or checked).
    #[test]
    fn rejects_tampered_fri_final_layer() {
        let air = two_boundary_air();
        let cfg = StarkConfig::fast();
        let proof = prove(&air, &cfg).unwrap();
        let mut p = proof.clone();
        p.fri.final_constant = fr_add(&p.fri.final_constant, &BigUint::one());
        assert!(!verify(&air, &cfg, &p));
        let mut p = proof.clone();
        let last = p.fri.layer_roots.len() - 1;
        p.fri.layer_roots[last] = [0u8; 32];
        assert!(!verify(&air, &cfg, &p));
    }

    /// First index whose position differs from query 0's.
    fn distinct_from_first(p: &StarkProofV2) -> usize {
        (1..p.trace_openings.len())
            .find(|&j| p.trace_openings[j].0 != p.trace_openings[0].0)
            .expect("two distinct positions")
    }

    /// Queries must come in the order the transcript dictates.
    #[test]
    fn rejects_swapped_queries() {
        let air = AirV2::doubling(8, 1);
        let cfg = StarkConfig::fast();
        let mut proof = prove(&air, &cfg).unwrap();
        let (a, b) = (0, distinct_from_first(&proof));
        proof.trace_openings.swap(a, b);
        proof.trace_next_openings.swap(a, b);
        proof.composition_openings.swap(a, b);
        proof.fri_openings.swap(a, b);
        assert!(!verify(&air, &cfg, &proof));
    }

    /// N6b: a query replaced by a copy of another.
    #[test]
    fn rejects_duplicated_query() {
        let air = two_boundary_air();
        let cfg = StarkConfig::fast();
        let mut p = prove(&air, &cfg).unwrap();
        let j = distinct_from_first(&p);
        p.trace_openings[j] = p.trace_openings[0].clone();
        p.trace_next_openings[j] = p.trace_next_openings[0].clone();
        p.composition_openings[j] = p.composition_openings[0].clone();
        p.fri_openings[j] = p.fri_openings[0].clone();
        assert!(!verify(&air, &cfg, &p));
    }

    /// N7a / N7b / N7c: the verifier's q and β are authoritative.
    #[test]
    fn rejects_query_count_or_blowup_mismatch() {
        let air = two_boundary_air();
        let cfg = StarkConfig::fast();
        let q10 = StarkConfig {
            query_count: 10,
            ..cfg.clone()
        };
        let p10 = prove(&air, &q10).unwrap();
        assert!(verify(&air, &q10, &p10));
        assert!(!verify(&air, &cfg, &p10));
        let mut p = prove(&air, &cfg).unwrap();
        p.trace_openings.truncate(10);
        p.trace_next_openings.truncate(10);
        p.composition_openings.truncate(10);
        p.fri_openings.truncate(10);
        assert!(!verify(&air, &cfg, &p));
        assert!(!verify(&air, &q10, &p), "a prefix is not a q = 10 proof");
        let b8 = StarkConfig {
            blow_up: 8,
            query_count: 20,
        };
        assert!(!verify(&air, &b8, &prove(&air, &cfg).unwrap()));
    }

    /// N9 / N10: a committed composition of zero, or one built with
    /// α₁ = 0, for the false statement 1 -> 12345.
    #[test]
    fn rejects_composition_not_bound_to_trace() {
        let cfg = StarkConfig::fast();
        let fa = false_output_air();
        for knobs in [
            ProverKnobs {
                composition_zero: true,
                ..ProverKnobs::default()
            },
            ProverKnobs {
                alpha1_zero: true,
                ..ProverKnobs::default()
            },
        ] {
            let p = prove_with(&fa, &cfg, &knobs).unwrap();
            assert!(!verify(&fa, &cfg, &p), "{knobs:?}");
        }
    }

    /// Is the composition, built with the given α₁, α₂, of degree < n on
    /// the coset (FRI folds to a constant final layer of β values)?
    fn composition_is_low_degree(air: &AirV2, beta: usize, a1: &BigUint, a2: &BigUint) -> bool {
        let n = air.n;
        let h = subgroup(n);
        let f = interpolate_from_subgroup(&h, &(air.generate_trace)());
        let l = n * beta;
        let hl = coset(l, &coset_shift());
        let fl = eval_on_domain(&f, &hl);
        let comp: Vec<BigUint> = (0..l)
            .map(|i| composition_at(air, &h, &hl[i], &fl[i], &fl[(i + beta) % l], a1, a2))
            .collect();
        let fri = fri_commit_to(&comp, &hl, &mut Vec::new(), beta);
        fri.layers
            .last()
            .unwrap()
            .iter()
            .all(|v| *v == fri.final_constant)
    }

    /// The frozen-heart forgery (re-audit N8a/N8b/N8c).  A cheating
    /// prover commits to a trace breaking the transition at row `p`,
    /// derives α₁, α₂ from everything it can absorb before choosing the
    /// public value at row `p` (the statement with a placeholder 0
    /// there), and picks the value `v` that cancels the composition's
    /// pole at ω^p under those α's.  Returns (false AIR, its proof from
    /// the real prover, whether the composition is low-degree under the
    /// attacker's α's).  Before statement binding the attacker's α's
    /// were the real ones and the proof was accepted.
    fn frozen_heart(
        cheat: &[BigUint],
        p: usize,
        other: &[(usize, BigUint)],
        k_index: usize,
        cfg: &StarkConfig,
        nudge: u64,
    ) -> (AirV2, StarkProofV2, bool) {
        let n = cheat.len();
        let with_value = |v: BigUint| {
            let mut boundary = other.to_vec();
            boundary.insert(k_index, (p, v));
            air2(n, boundary, dbl(), cheat.to_vec())
        };
        let placeholder = with_value(BigUint::zero());
        let tc = prove(&placeholder, cfg).unwrap().trace_commit;
        let mut t = statement_transcript(&placeholder, cfg).unwrap();
        let (a1, a2) = absorb_trace_draw_alphas(&mut t, &tc);
        let h = subgroup(n);
        let mut zp = BigUint::one();
        for (i, hi) in h.iter().enumerate().take(n - 1) {
            if i != p {
                zp = fr_mul(&zp, &fr_sub(&h[p], hi));
            }
        }
        let c_p = (dbl())(&cheat[p], &cheat[p + 1]);
        let mut a2k = BigUint::one();
        for _ in 0..k_index {
            a2k = fr_mul(&a2k, &a2);
        }
        let shift = fr_mul(&fr_mul(&a1, &c_p), &fr_inv(&fr_mul(&a2k, &zp)).unwrap());
        let v = fr_add(&fr_add(&cheat[p], &shift), &big(nudge));
        let air_false = with_value(v);
        let low = composition_is_low_degree(&air_false, cfg.blow_up, &a1, &a2);
        let proof = prove(&air_false, cfg).unwrap();
        (air_false, proof, low)
    }

    /// N8a (fast and default configs): "doubling with t0 = v, t7 = 128"
    /// for v != 1, from a trace [5, 2, 4, .., 128] that breaks row 0.
    #[test]
    fn rejects_frozen_heart_forged_input() {
        let mut cheat = doubling_trace(8, 1);
        cheat[0] = big(5);
        for cfg in [StarkConfig::fast(), StarkConfig::default()] {
            let (af, pf, low) = frozen_heart(&cheat, 0, &[(7, big(128))], 0, &cfg, 0);
            // Attack control: under the α's the attacker can compute the
            // pole cancels, so only statement binding stops this proof.
            assert!(low, "forgery machinery broken");
            assert_ne!(af.boundary[0].1, big(1), "statement must be false");
            assert!(!verify(&af, &cfg, &pf), "frozen-heart forgery accepted");
        }
    }

    /// N8b: "doubling with t0 = 1, t6 = v" for v != 64 (row 6 broken).
    #[test]
    fn rejects_frozen_heart_forged_row6() {
        let mut cheat = doubling_trace(8, 1);
        cheat[7] = big(999);
        let cfg = StarkConfig::fast();
        let (af, pf, low) = frozen_heart(&cheat, 6, &[(0, big(1))], 1, &cfg, 0);
        assert!(low, "forgery machinery broken");
        assert_ne!(af.boundary[1].1, big(64));
        assert!(!verify(&af, &cfg, &pf));
    }

    /// N8c: the attack control.  With v + 1 the pole does not cancel even
    /// under the attacker's α's, so the composition is not low degree.
    #[test]
    fn frozen_heart_control_v_plus_one_is_not_low_degree() {
        let mut cheat = doubling_trace(8, 1);
        cheat[0] = big(5);
        let cfg = StarkConfig::fast();
        let (af, pf, low) = frozen_heart(&cheat, 0, &[(7, big(128))], 0, &cfg, 1);
        assert!(!low);
        assert!(!verify(&af, &cfg, &pf));
    }

    /// M2: non-canonical encodings (v + r) of opened values are
    /// rejected; so are non-canonical statement values.
    #[test]
    fn rejects_non_canonical_values() {
        let r = scalar_modulus();
        let air = two_boundary_air();
        let cfg = StarkConfig::fast();
        let proof = prove(&air, &cfg).unwrap();
        let mut p = proof.clone();
        p.trace_openings[0].1 += &r;
        assert!(!verify(&air, &cfg, &p));
        let mut p = proof.clone();
        p.trace_next_openings[0].1 += &r;
        assert!(!verify(&air, &cfg, &p));
        let mut p = proof.clone();
        p.composition_openings[0].1 += &r;
        p.fri_openings[0].1.openings[0].0 += &r;
        assert!(!verify(&air, &cfg, &p));
        let mut p = proof.clone();
        p.fri_openings[0].1.openings[0].1 += &r;
        assert!(!verify(&air, &cfg, &p));
        let shifted = air2(
            8,
            vec![(0, BigUint::one() + &r), (7, big(128))],
            dbl(),
            doubling_trace(8, 1),
        );
        assert!(!verify(&shifted, &cfg, &proof));
        assert!(prove(&shifted, &cfg).is_none());
    }

    /// M3: an opened value of 2^300 returns false instead of panicking,
    /// wherever it appears; so does an out-of-range boundary position.
    #[test]
    fn oversized_values_return_false_without_panic() {
        let air = two_boundary_air();
        let cfg = StarkConfig::fast();
        let proof = prove(&air, &cfg).unwrap();
        let huge: BigUint = BigUint::one() << 300;
        let mut p = proof.clone();
        p.trace_openings[0].1 = huge.clone();
        assert!(!verify(&air, &cfg, &p));
        let mut p = proof.clone();
        p.composition_openings[0].1 = huge.clone();
        assert!(!verify(&air, &cfg, &p));
        let mut p = proof.clone();
        p.fri_openings[0].1.openings[1].0 = huge.clone();
        assert!(!verify(&air, &cfg, &p));
        let mut p = proof.clone();
        p.fri.final_constant = huge.clone();
        assert!(!verify(&air, &cfg, &p));
        let big_value = air2(8, vec![(0, huge)], dbl(), doubling_trace(8, 1));
        assert!(!verify(&big_value, &cfg, &proof));
        let bad_pos = air2(8, vec![(8, big(1))], dbl(), doubling_trace(8, 1));
        assert!(!verify(&bad_pos, &cfg, &proof));
        assert!(prove(&bad_pos, &cfg).is_none());
    }
}
