//! # Symmetrised decomposition systems for Koblitz curves.
//!
//! The follow-through on `RESEARCH_EXOTIC_COORDINATES.md`: the coordinate
//! search found that on a Koblitz curve the rational 2-torsion point
//! `T = (0, 1)` acts on the `x`-line as `x ↦ 1/x`, that the frame
//! `u = 1/(x + 1)` (defined over `F₂`, so the Frobenius survives) turns
//! that into `u ↦ u + 1`, and that the summation polynomial written in the
//! invariants `w = u² + u`, `s = Σu` has half the degree.  That was a
//! statement about polynomials.  This module turns it into decomposition
//! systems the existing oracles can solve, so the claim "smaller system"
//! can be tested against the claim that matters, "faster solve".
//!
//! ## The system
//!
//! Factor base `F_u = { P : u(P) ∈ V }` for a Frobenius-stable subspace
//! `V ∋ 1` (the divisor of `xⁿ − 1` must contain `x − 1`).  Choose a basis
//! `1, b₂, …, b_ℓ` of `V`.  Every `u ∈ V` is `ε + Σ_t c_t b_t`, and the
//! Artin–Schreier map `AS(v) = v² + v` is `F₂`-linear with kernel `{0, 1}`,
//! so `w = AS(u) = Σ_t c_t AS(b_t)` — the bit `ε` drops out.  For `m`
//! summands the unknowns are the `c_{i,t}` (`m(ℓ − 1)` bits) and one bit
//! `ε = Σ_i ε_i`; `w_i` is linear in them, and so is `s = Σ_i u_i + u_R`.
//! The symmetrised polynomial, interpolated and verified in
//! [`crate::cryptanalysis::coordinate_search`], is then Weil-restricted
//! exactly as the `x`-system is: `n` Boolean equations.
//!
//! ```text
//!   m = 2:  w₁ w₂ w_R + w₁ + w₂ + w_R + s = 0                (bilinear, s linear)
//!   m = 3:  the 18-term S₄ of degree 2 in each w, s linear  (Boolean degree 4)
//! ```
//!
//! A root gives each `w_i`, hence `u_i` up to `u_i ↦ u_i + 1` with the
//! parity fixed by `ε`, hence `x_i = 1/u_i + 1` up to sign.  The lifted
//! points sum to `R` or to `R + T`; both are relations over `F_u` since
//! `T ∈ F_u` (`u(T) = 1 ∈ V`), and the oracle reports which.
//!
//! ## What is compared
//!
//! [`paired_bench`] runs, on the same curve and the same subspace `V`:
//!
//! - **x-chained**: the production system,
//!   [`crate::cryptanalysis::koblitz_groebner::build_decomposition_system`]
//!   on `F_x = {x(P) ∈ V}` — `S₃` for `m = 2`, chained `S₃` links with
//!   `(m − 2)·n` free intermediate unknowns for `m = 3`;
//! - **x-direct** (`m = 3`): the 24-term `S₄` in `x`, no chaining, Boolean
//!   degree 7 — separates "no chaining" from "symmetry";
//! - **sym**: the symmetrised system on `F_u`.
//!
//! Every verdict is gated against exhaustive enumeration on its own base,
//! every returned relation is re-summed in the group, and timings are
//! reported separately for found and refuted targets, as the scaling
//! harness does.  The two bases differ as sets — same `V`, different
//! coordinate — so their found/refuted mixes differ too; the counts are
//! printed so nobody compares a median across regimes.
//!
//! ## Honest scope
//!
//! Toy `n ≤ 24` (the curve constructor's limit), `m ≤ 3`.  The
//! Kosters–Yeo trace row the SAT path adds to the `x`-system has no
//! linear analogue in `u`, so the SAT arms run without it on both sides.

use std::collections::HashMap;
use std::time::Instant;

use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::binary_ecc::curve::point_neg;
use crate::binary_ecc::{BinaryPoint, F2mElement};
use crate::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, first_fall_degree, solve_boolean_system_filtered, FieldStructure,
    SolveOptions, SolveStats, SolverEngine, SymElement, MAX_VARS,
};
use crate::cryptanalysis::koblitz_index_calculus::{
    all_factors_of_x_n_minus_1, build_frobenius_factor_base_from_divisor, enumerate_decompose,
    groebner_decompose, point_key, points_with_x, span_f2, subspace_basis_for_divisor,
    KoblitzCurve,
};
use crate::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use crate::cryptanalysis::sat::SolveResult;
use crate::cryptanalysis::semaev_sat::{encode_boolean_system_with, XorEncoding};

fn bits_of(e: &F2mElement) -> u64 {
    e.raw_bits().first().copied().unwrap_or(0)
}

fn from_bits(v: u64, n: u32) -> F2mElement {
    F2mElement::from_biguint(&BigUint::from(v), n)
}

// ── The frame ──────────────────────────────────────────────────────

/// `u = 1/(x + 1)`; `None` at `x = 1`, where `u = ∞`.
pub fn u_of_x(x: &F2mElement, kc: &KoblitzCurve) -> Option<F2mElement> {
    let d = x.add(&F2mElement::one(kc.n));
    d.flt_inverse(&kc.curve.irreducible)
}

/// `x = 1/u + 1`; `None` at `u = 0`, where `x = ∞`.
pub fn x_of_u(u: &F2mElement, kc: &KoblitzCurve) -> Option<F2mElement> {
    Some(
        u.flt_inverse(&kc.curve.irreducible)?
            .add(&F2mElement::one(kc.n)),
    )
}

/// `w = u² + u`, the Artin–Schreier invariant of `u ↦ u + 1`.
pub fn artin_schreier(u: &F2mElement, kc: &KoblitzCurve) -> F2mElement {
    u.square(&kc.curve.irreducible).add(u)
}

// ── Factor base in the u-frame ─────────────────────────────────────

/// `F_u = { P : u(P) ∈ V }` with its bases.
#[derive(Clone, Debug)]
pub struct SymmetrisedFactorBase {
    /// Basis of `V`; `v_basis[0] = 1`.
    pub v_basis: Vec<F2mElement>,
    /// `AS(v_basis[t])` for `t ≥ 1`: a basis of `AS(V)`, dimension `ℓ − 1`.
    pub w_basis: Vec<F2mElement>,
    /// `ℓ = dim V`.
    pub ell: usize,
    pub points: Vec<BinaryPoint>,
    pub index_of: HashMap<(BigUint, BigUint), usize>,
    /// `T = (0, 1)`.
    pub two_torsion: BinaryPoint,
    pub two_torsion_index: usize,
}

/// Indices into [`all_factors_of_x_n_minus_1`] of a divisor containing
/// `x − 1` whose root space has dimension as close as possible to
/// `target_dim` (ties to the smaller).  `None` if `x − 1` is not a listed
/// factor.
pub fn divisor_for_dimension(n: u32, target_dim: u32) -> Option<Vec<usize>> {
    let factors = all_factors_of_x_n_minus_1(n);
    let x1 = factors.iter().position(|&f| f == 0b11)?;
    let others: Vec<(usize, u32)> = factors
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != x1)
        .map(|(i, &f)| (i, 63 - f.leading_zeros()))
        .collect();
    let mut best: Option<(u32, Vec<usize>)> = None;
    for mask in 0..(1u32 << others.len()) {
        let dim = 1 + others
            .iter()
            .enumerate()
            .filter(|(k, _)| (mask >> k) & 1 == 1)
            .map(|(_, &(_, d))| d)
            .sum::<u32>();
        if dim >= n {
            continue;
        }
        let dist = (dim as i64 - target_dim as i64).abs();
        let better = match &best {
            None => true,
            Some((bd, _)) => {
                let bdist = (*bd as i64 - target_dim as i64).abs();
                dist < bdist || (dist == bdist && dim < *bd)
            }
        };
        if better {
            let mut idx = vec![x1];
            idx.extend(
                others
                    .iter()
                    .enumerate()
                    .filter(|(k, _)| (mask >> k) & 1 == 1)
                    .map(|(_, &(i, _))| i),
            );
            best = Some((dim, idx));
        }
    }
    best.map(|(_, idx)| idx)
}

/// Build `F_u` from a divisor of `xⁿ − 1` (indices as in
/// [`all_factors_of_x_n_minus_1`]).  Returns `None` if `1 ∉ V`.
pub fn build_symmetrised_factor_base(
    kc: &KoblitzCurve,
    divisor_indices: &[usize],
) -> Option<SymmetrisedFactorBase> {
    let n = kc.n;
    let raw = subspace_basis_for_divisor(n, divisor_indices, &kc.curve.irreducible)?;
    let span = span_f2(&raw, n);
    let span_bits: std::collections::HashSet<u64> = span.iter().map(bits_of).collect();
    if !span_bits.contains(&1) {
        return None;
    }
    // Re-basis with 1 first: greedy independence over the F_2-span.
    let mut v_basis: Vec<F2mElement> = vec![F2mElement::one(n)];
    let mut current: std::collections::HashSet<u64> = [0u64, 1].into_iter().collect();
    for b in &raw {
        let bb = bits_of(b);
        if current.contains(&bb) {
            continue;
        }
        let added: Vec<u64> = current.iter().map(|&c| c ^ bb).collect();
        current.extend(added);
        v_basis.push(b.clone());
    }
    let ell = v_basis.len();
    if ell != raw.len() {
        return None;
    }
    let w_basis: Vec<F2mElement> = v_basis[1..].iter().map(|b| artin_schreier(b, kc)).collect();
    if span_f2(&w_basis, n).len() != 1 << (ell - 1) {
        return None;
    }
    let mut points = Vec::new();
    for u in &span {
        if u.is_zero() {
            continue;
        }
        let x = x_of_u(u, kc)?;
        points.extend(points_with_x(&kc.curve, &x));
    }
    let index_of: HashMap<(BigUint, BigUint), usize> = points
        .iter()
        .enumerate()
        .map(|(i, p)| (point_key(p), i))
        .collect();
    let two_torsion = points_with_x(&kc.curve, &F2mElement::zero(n))
        .into_iter()
        .next()?;
    let two_torsion_index = *index_of.get(&point_key(&two_torsion))?;
    Some(SymmetrisedFactorBase {
        v_basis,
        w_basis,
        ell,
        points,
        index_of,
        two_torsion,
        two_torsion_index,
    })
}

impl SymmetrisedFactorBase {
    /// Sum of the indexed points.
    pub fn sum(&self, kc: &KoblitzCurve, idxs: &[usize]) -> BinaryPoint {
        idxs.iter().fold(BinaryPoint::Infinity, |acc, &i| {
            kc.add(&acc, &self.points[i])
        })
    }
}

/// Exhaustive decomposition over `F_u`: `Some((indices, used_t))` with
/// `Σ = R` (`used_t = false`) or `Σ = R + T` (`used_t = true`, and `T`'s
/// index appended so the returned list really sums to `R`).
pub fn enumerate_symmetrised(
    kc: &KoblitzCurve,
    fb: &SymmetrisedFactorBase,
    target: &BinaryPoint,
    m: usize,
) -> Option<(Vec<usize>, bool)> {
    fn rec(
        kc: &KoblitzCurve,
        fb: &SymmetrisedFactorBase,
        target: &BinaryPoint,
        m: usize,
        start: usize,
    ) -> Option<Vec<usize>> {
        if m == 1 {
            let idx = *fb.index_of.get(&point_key(target))?;
            return (idx >= start).then(|| vec![idx]);
        }
        for i in start..fb.points.len() {
            let rest = kc.add(target, &point_neg(&fb.points[i]));
            if let Some(mut tail) = rec(kc, fb, &rest, m - 1, i) {
                let mut out = vec![i];
                out.append(&mut tail);
                return Some(out);
            }
        }
        None
    }
    if let Some(v) = rec(kc, fb, target, m, 0) {
        return Some((v, false));
    }
    let shifted = kc.add(target, &fb.two_torsion);
    rec(kc, fb, &shifted, m, 0).map(|mut v| {
        v.push(fb.two_torsion_index);
        (v, true)
    })
}

// ── The polynomials ────────────────────────────────────────────────

/// Exponent vectors of the symmetrised summation polynomial over
/// `(w₁, …, w_m, w_R, s)`, all coefficients `1`.  Interpolated by
/// `coordinate_search` on a Koblitz curve; `F₂` coefficients, so valid
/// for every `n` and both `a`.  A test pins them to the interpolation.
pub fn symmetrised_terms(m: usize) -> Option<Vec<Vec<u32>>> {
    Some(match m {
        2 => vec![
            vec![1, 1, 1, 0],
            vec![1, 0, 0, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![0, 0, 0, 1],
        ],
        3 => vec![
            vec![2, 2, 2, 0, 0],
            vec![2, 2, 0, 2, 0],
            vec![2, 0, 2, 2, 0],
            vec![0, 2, 2, 2, 0],
            vec![2, 1, 1, 1, 0],
            vec![1, 2, 1, 1, 0],
            vec![1, 1, 2, 1, 0],
            vec![1, 1, 1, 2, 0],
            vec![1, 1, 1, 1, 1],
            vec![2, 0, 0, 0, 0],
            vec![0, 2, 0, 0, 0],
            vec![0, 0, 2, 0, 0],
            vec![0, 0, 0, 2, 0],
            vec![1, 0, 0, 0, 0],
            vec![0, 1, 0, 0, 0],
            vec![0, 0, 1, 0, 0],
            vec![0, 0, 0, 1, 0],
            vec![0, 0, 0, 0, 1],
        ],
        _ => return None,
    })
}

/// Exponent vectors of Semaev's `S_{m+1}` in `x` over `(x₁, …, x_m, x_R)`
/// for `b = 1`, all coefficients `1` (the constant term is `b`).
pub fn plain_terms(m: usize) -> Option<Vec<Vec<u32>>> {
    Some(match m {
        2 => vec![
            vec![2, 2, 0],
            vec![2, 0, 2],
            vec![0, 2, 2],
            vec![1, 1, 1],
            vec![0, 0, 0],
        ],
        3 => vec![
            vec![4, 4, 4, 0],
            vec![4, 4, 0, 4],
            vec![4, 0, 4, 4],
            vec![0, 4, 4, 4],
            vec![4, 2, 2, 2],
            vec![3, 3, 3, 1],
            vec![3, 3, 1, 3],
            vec![3, 1, 3, 3],
            vec![2, 4, 2, 2],
            vec![2, 2, 4, 2],
            vec![2, 2, 2, 4],
            vec![1, 3, 3, 3],
            vec![3, 1, 1, 1],
            vec![2, 2, 2, 0],
            vec![2, 2, 0, 2],
            vec![2, 0, 2, 2],
            vec![1, 3, 1, 1],
            vec![1, 1, 3, 1],
            vec![1, 1, 1, 3],
            vec![0, 2, 2, 2],
            vec![4, 0, 0, 0],
            vec![0, 4, 0, 0],
            vec![0, 0, 4, 0],
            vec![0, 0, 0, 4],
        ],
        _ => return None,
    })
}

/// `v^e` for a symbolic field element: `e = 2^k` is `k` linear squarings,
/// odd exponents multiply once more.
fn sym_pow(v: &SymElement, e: u32, st: &FieldStructure) -> SymElement {
    let mut result: Option<SymElement> = None;
    let mut base = v.clone();
    let mut e = e;
    while e > 0 {
        if e & 1 == 1 {
            result = Some(match result {
                None => base.clone(),
                Some(r) => r.mul(&base, st),
            });
        }
        e >>= 1;
        if e > 0 {
            base = base.square(st);
        }
    }
    result.expect("positive exponent")
}

/// Weil restriction of `Σ_terms Π_i vars[i]^{e_i} = 0` to `n` Boolean
/// equations.  Squares are linear and cost nothing; the Boolean degree of
/// a term is `Σ popcount(e_i)` over the non-constant variables.
pub fn sym_polynomial(
    vars: &[SymElement],
    terms: &[Vec<u32>],
    st: &FieldStructure,
) -> Vec<F2BoolPoly> {
    let n_vars = vars[0].coords[0].n_vars;
    let mut acc = SymElement::zero(st.n, n_vars);
    for term in terms {
        let mut prod: Option<SymElement> = None;
        for (i, &e) in term.iter().enumerate() {
            if e == 0 {
                continue;
            }
            let p = sym_pow(&vars[i], e, st);
            prod = Some(match prod {
                None => p,
                Some(q) => q.mul(&p, st),
            });
        }
        let t = prod.unwrap_or_else(|| SymElement::constant(&F2mElement::one(st.n), st.n, n_vars));
        acc = acc.add(&t);
    }
    acc.coords
}

// ── Systems ────────────────────────────────────────────────────────

/// A decomposition system plus what its roots mean.
#[derive(Clone, Debug)]
pub struct SymmetrisedSystem {
    pub equations: Vec<F2BoolPoly>,
    pub n_vars: usize,
    /// Bits per summand: `ℓ − 1`.
    pub ell_w: usize,
    pub m: usize,
    /// `u(R)`.
    pub u_r: F2mElement,
}

impl SymmetrisedSystem {
    /// The `u`-coordinate of summand `i` that a root gives, before the
    /// `ε_i` choice: `Σ_t c_{i,t} b_{t+1}`.
    pub fn summand_u0(
        &self,
        fb: &SymmetrisedFactorBase,
        root: u64,
        i: usize,
        n: u32,
    ) -> F2mElement {
        let mut acc = F2mElement::zero(n);
        for t in 0..self.ell_w {
            if (root >> (i * self.ell_w + t)) & 1 == 1 {
                acc = acc.add(&fb.v_basis[t + 1]);
            }
        }
        acc
    }
    /// The parity bit `ε = Σ ε_i`.
    pub fn parity(&self, root: u64) -> bool {
        (root >> (self.m * self.ell_w)) & 1 == 1
    }
}

/// Build the symmetrised system for target `R` and `m` summands.  `None`
/// if `m ∉ {2, 3}`, `u(R)` is `0` or `∞`, or the layout exceeds
/// [`MAX_VARS`].
pub fn build_symmetrised_system(
    kc: &KoblitzCurve,
    fb: &SymmetrisedFactorBase,
    target: &BinaryPoint,
    m: usize,
    st: &FieldStructure,
) -> Option<SymmetrisedSystem> {
    let terms = symmetrised_terms(m)?;
    let x_r = match target {
        BinaryPoint::Affine { x, .. } => x,
        BinaryPoint::Infinity => return None,
    };
    let u_r = u_of_x(x_r, kc)?;
    if u_r.is_zero() {
        return None;
    }
    let n = kc.n;
    let ell_w = fb.ell - 1;
    let n_vars = m * ell_w + 1;
    if n_vars > MAX_VARS {
        return None;
    }
    let eps = (m * ell_w) as u32;
    let mut vars: Vec<SymElement> = (0..m)
        .map(|i| SymElement::from_subspace_vars(&fb.w_basis, i * ell_w, n, n_vars))
        .collect();
    vars.push(SymElement::constant(&artin_schreier(&u_r, kc), n, n_vars));
    // s = ε + Σ_i Σ_t c_{i,t} b_{t+1} + u_R
    let mut s = SymElement::constant(&u_r, n, n_vars);
    for i in 0..m {
        s = s.add(&SymElement::from_subspace_vars(
            &fb.v_basis[1..],
            i * ell_w,
            n,
            n_vars,
        ));
    }
    let mut eps_el = SymElement::zero(n, n_vars);
    eps_el.coords[0] = F2BoolPoly::from_monos(vec![F2BoolMono::var(eps)], n_vars);
    s = s.add(&eps_el);
    vars.push(s);
    let equations = sym_polynomial(&vars, &terms, st);
    Some(SymmetrisedSystem {
        equations,
        n_vars,
        ell_w,
        m,
        u_r,
    })
}

/// The direct (unchained) `x`-system: Semaev's `S_{m+1}` in `x` with every
/// summand in the span of `basis`.  Same variable layout as the chained
/// system's first `m·ℓ` bits, so
/// [`crate::cryptanalysis::koblitz_groebner::DecompositionSystem::summand_x`]
/// logic applies.
pub fn build_direct_x_system(
    basis: &[F2mElement],
    x_r: &F2mElement,
    m: usize,
    st: &FieldStructure,
) -> Option<(Vec<F2BoolPoly>, usize)> {
    let terms = plain_terms(m)?;
    let ell = basis.len();
    let n_vars = m * ell;
    if n_vars > MAX_VARS {
        return None;
    }
    let mut vars: Vec<SymElement> = (0..m)
        .map(|i| SymElement::from_subspace_vars(basis, i * ell, st.n, n_vars))
        .collect();
    vars.push(SymElement::constant(x_r, st.n, n_vars));
    Some((sym_polynomial(&vars, &terms, st), n_vars))
}

fn subspace_x(basis: &[F2mElement], root: u64, i: usize, n: u32) -> F2mElement {
    let ell = basis.len();
    let mut acc = F2mElement::zero(n);
    for t in 0..ell {
        if (root >> (i * ell + t)) & 1 == 1 {
            acc = acc.add(&basis[t]);
        }
    }
    acc
}

// ── Lifting ────────────────────────────────────────────────────────

/// Walk sign choices over `points_with_x` for each abscissa; succeed when
/// the sum hits `target`.
fn lift_signs(
    kc: &KoblitzCurve,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    points: &[BinaryPoint],
    xs: &[F2mElement],
    target: &BinaryPoint,
) -> Option<Vec<usize>> {
    fn walk(
        kc: &KoblitzCurve,
        index_of: &HashMap<(BigUint, BigUint), usize>,
        points: &[BinaryPoint],
        xs: &[F2mElement],
        depth: usize,
        acc: &BinaryPoint,
        chosen: &mut Vec<usize>,
        target: &BinaryPoint,
    ) -> bool {
        if depth == xs.len() {
            return acc == target;
        }
        for p in points_with_x(&kc.curve, &xs[depth]) {
            let Some(&idx) = index_of.get(&point_key(&p)) else {
                continue;
            };
            chosen.push(idx);
            let next = kc.add(acc, &points[idx]);
            if walk(kc, index_of, points, xs, depth + 1, &next, chosen, target) {
                return true;
            }
            chosen.pop();
        }
        false
    }
    let mut chosen = Vec::new();
    walk(
        kc,
        index_of,
        points,
        xs,
        0,
        &BinaryPoint::Infinity,
        &mut chosen,
        target,
    )
    .then_some(chosen)
}

/// Lift a root of the symmetrised system to factor-base points summing to
/// `R` (`used_t = false`) or to `R + T` (`used_t = true`, `T` appended).
pub fn lift_symmetrised_root(
    kc: &KoblitzCurve,
    fb: &SymmetrisedFactorBase,
    sys: &SymmetrisedSystem,
    root: u64,
    target: &BinaryPoint,
) -> Option<(Vec<usize>, bool)> {
    let m = sys.m;
    let u0: Vec<F2mElement> = (0..m).map(|i| sys.summand_u0(fb, root, i, kc.n)).collect();
    let parity = sys.parity(root);
    let shifted = kc.add(target, &fb.two_torsion);
    for mask in 0..(1u32 << m) {
        if (mask.count_ones() % 2 == 1) != parity {
            continue;
        }
        let mut xs = Vec::with_capacity(m);
        let mut ok = true;
        for i in 0..m {
            let mut u = u0[i].clone();
            if (mask >> i) & 1 == 1 {
                u = u.add(&F2mElement::one(kc.n));
            }
            match x_of_u(&u, kc) {
                Some(x) => xs.push(x),
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if !ok {
            continue;
        }
        if let Some(v) = lift_signs(kc, &fb.index_of, &fb.points, &xs, target) {
            return Some((v, false));
        }
        if let Some(mut v) = lift_signs(kc, &fb.index_of, &fb.points, &xs, &shifted) {
            v.push(fb.two_torsion_index);
            return Some((v, true));
        }
    }
    None
}

// ── Oracles ────────────────────────────────────────────────────────

/// Outcome of one oracle call.
#[derive(Clone, Debug)]
pub struct OracleOutcome {
    /// Factor-base indices summing to the target (with `T` appended when
    /// the root lifted to `R + T`).
    pub relation: Option<Vec<usize>>,
    /// The lifted relation went through `R + T`.
    pub used_t: bool,
    /// The solve completed: a `None` relation is a proof of absence.
    pub complete: bool,
    pub n_vars: usize,
    pub n_equations: usize,
    /// Boolean degree of the system.
    pub degree: u32,
    /// Engine-specific effort: F4 splits, or SAT conflicts.
    pub effort: u64,
}

fn system_degree(eqs: &[F2BoolPoly]) -> u32 {
    eqs.iter()
        .flat_map(|p| p.terms.iter())
        .map(|t| t.mask.count_ones())
        .max()
        .unwrap_or(0)
}

/// Symmetrised decomposition by matrix-F4 with splitting.
pub fn symmetrised_groebner_decompose(
    kc: &KoblitzCurve,
    fb: &SymmetrisedFactorBase,
    st: &FieldStructure,
    target: &BinaryPoint,
    m: usize,
    engine: SolverEngine,
    node_budget: usize,
) -> Option<OracleOutcome> {
    let sys = build_symmetrised_system(kc, fb, target, m, st)?;
    let opts = SolveOptions {
        engine,
        max_solutions: usize::MAX,
        node_budget,
    };
    let mut found: Option<(Vec<usize>, bool)> = None;
    let (_, stats): (Vec<u64>, SolveStats) =
        solve_boolean_system_filtered(&sys.equations, sys.n_vars, &opts, |root| {
            match lift_symmetrised_root(kc, fb, &sys, root, target) {
                Some(r) => {
                    found = Some(r);
                    true
                }
                None => false,
            }
        });
    let used_t = found.as_ref().map(|r| r.1).unwrap_or(false);
    Some(OracleOutcome {
        relation: found.map(|r| r.0),
        used_t,
        complete: !stats.exhausted,
        n_vars: sys.n_vars,
        n_equations: sys.equations.len(),
        degree: system_degree(&sys.equations),
        effort: stats.splits as u64,
    })
}

/// Direct (unchained) `x`-system decomposition by matrix-F4 with
/// splitting, on the `x`-base `fb_x` (its subspace basis and points).
pub fn direct_x_groebner_decompose(
    kc: &KoblitzCurve,
    basis: &[F2mElement],
    points: &[BinaryPoint],
    index_of: &HashMap<(BigUint, BigUint), usize>,
    st: &FieldStructure,
    target: &BinaryPoint,
    m: usize,
    engine: SolverEngine,
    node_budget: usize,
) -> Option<OracleOutcome> {
    let x_r = match target {
        BinaryPoint::Affine { x, .. } => x.clone(),
        BinaryPoint::Infinity => return None,
    };
    let (equations, n_vars) = build_direct_x_system(basis, &x_r, m, st)?;
    let opts = SolveOptions {
        engine,
        max_solutions: usize::MAX,
        node_budget,
    };
    let mut found: Option<Vec<usize>> = None;
    let (_, stats) = solve_boolean_system_filtered(&equations, n_vars, &opts, |root| {
        let xs: Vec<F2mElement> = (0..m).map(|i| subspace_x(basis, root, i, kc.n)).collect();
        match lift_signs(kc, index_of, points, &xs, target) {
            Some(v) => {
                found = Some(v);
                true
            }
            None => false,
        }
    });
    Some(OracleOutcome {
        relation: found,
        used_t: false,
        complete: !stats.exhausted,
        n_vars,
        n_equations: equations.len(),
        degree: system_degree(&equations),
        effort: stats.splits as u64,
    })
}

/// Generic CDCL loop over a Boolean system: enumerate models, lift each
/// with `lift`, block the ones that do not lift.
fn sat_loop(
    equations: &[F2BoolPoly],
    n_vars: usize,
    conflict_budget: u64,
    max_models: usize,
    mut lift: impl FnMut(u64) -> Option<(Vec<usize>, bool)>,
) -> OracleOutcome {
    let mut enc = encode_boolean_system_with(n_vars, equations, &[], XorEncoding::Native);
    enc.solver.conflict_budget = conflict_budget;
    let degree = system_degree(equations);
    let mut models = 0usize;
    let outcome = |relation: Option<(Vec<usize>, bool)>, complete: bool, effort: u64| {
        let used_t = relation.as_ref().map(|r| r.1).unwrap_or(false);
        OracleOutcome {
            relation: relation.map(|r| r.0),
            used_t,
            complete,
            n_vars,
            n_equations: equations.len(),
            degree,
            effort,
        }
    };
    loop {
        if enc.trivially_unsat {
            return outcome(None, true, 0);
        }
        match enc.solver.solve() {
            SolveResult::Unsat => return outcome(None, true, enc.solver.conflicts()),
            SolveResult::Unknown => return outcome(None, false, enc.solver.conflicts()),
            SolveResult::Sat => {
                let root = enc.model_assignment();
                if !equations.iter().all(|e| e.eval(root) == 0) {
                    // A bad model is an encoding bug; fail closed.
                    return outcome(None, false, enc.solver.conflicts());
                }
                models += 1;
                if let Some(r) = lift(root) {
                    return outcome(Some(r), true, enc.solver.conflicts());
                }
                if models >= max_models {
                    return outcome(None, false, enc.solver.conflicts());
                }
                let clause: Vec<i32> = (0..n_vars)
                    .map(|i| {
                        let lit = (i + 1) as i32;
                        if (root >> i) & 1 == 1 {
                            -lit
                        } else {
                            lit
                        }
                    })
                    .collect();
                enc.solver.add_clause(clause);
            }
        }
    }
}

/// Symmetrised decomposition by CDCL SAT with native XOR rows.
pub fn symmetrised_sat_decompose(
    kc: &KoblitzCurve,
    fb: &SymmetrisedFactorBase,
    st: &FieldStructure,
    target: &BinaryPoint,
    m: usize,
    conflict_budget: u64,
    max_models: usize,
) -> Option<OracleOutcome> {
    let sys = build_symmetrised_system(kc, fb, target, m, st)?;
    Some(sat_loop(
        &sys.equations,
        sys.n_vars,
        conflict_budget,
        max_models,
        |root| lift_symmetrised_root(kc, fb, &sys, root, target),
    ))
}

/// The production `x`-system (chained for `m = 3`) by CDCL SAT, without
/// the trace row, so it is the same encoding as the symmetrised arm.
pub fn chained_x_sat_decompose(
    kc: &KoblitzCurve,
    basis: &[F2mElement],
    points: &[BinaryPoint],
    index_of: &HashMap<(BigUint, BigUint), usize>,
    st: &FieldStructure,
    target: &BinaryPoint,
    m: usize,
    conflict_budget: u64,
    max_models: usize,
) -> Option<OracleOutcome> {
    let x_r = match target {
        BinaryPoint::Affine { x, .. } => x.clone(),
        BinaryPoint::Infinity => return None,
    };
    let sys = build_decomposition_system(basis, &x_r, &kc.curve.b, m, st)?;
    Some(sat_loop(
        &sys.equations,
        sys.n_vars,
        conflict_budget,
        max_models,
        |root| {
            let xs: Vec<F2mElement> = (0..m)
                .map(|i| sys.summand_x(basis, root, i, kc.n))
                .collect();
            lift_signs(kc, index_of, points, &xs, target).map(|v| (v, false))
        },
    ))
}

// ── Paired benchmark ───────────────────────────────────────────────

/// One arm's aggregate over the targets.
#[derive(Clone, Debug)]
pub struct ArmRun {
    pub arm: &'static str,
    pub n_vars: usize,
    pub n_equations: usize,
    pub degree: u32,
    pub found: usize,
    pub refuted: usize,
    pub inconclusive: usize,
    /// Relations that lifted through `R + T`.
    pub via_t: usize,
    pub median_found_ms: f64,
    pub median_refuted_ms: f64,
    pub total_ms: f64,
    pub median_effort: f64,
    /// Verdicts contradicting exhaustive enumeration on the arm's own
    /// base, plus relations that did not re-sum to the target.
    pub gate_failures: usize,
    /// First fall degree of the system at the first target, if computed.
    pub first_fall_degree: Option<u32>,
}

/// The paired comparison on one instance.
#[derive(Clone, Debug)]
pub struct PairedBench {
    pub a: u8,
    pub n: u32,
    pub ell: usize,
    pub m: usize,
    pub divisor: Vec<usize>,
    pub fb_x_size: usize,
    pub fb_u_size: usize,
    pub targets: usize,
    pub arms: Vec<ArmRun>,
}

/// Knobs for [`paired_bench`].
#[derive(Clone, Debug)]
pub struct PairedOptions {
    pub targets: usize,
    pub seed: u64,
    pub node_budget: usize,
    pub sat: bool,
    pub conflict_budget: u64,
    pub max_models: usize,
    /// Include the direct (unchained) `x`-`S₄` arm for `m = 3`.
    pub direct_x: bool,
    /// Compute first fall degrees up to this Macaulay degree (0 = skip).
    pub ffd_max_degree: u32,
    /// Draw targets from the prime-order subgroup (as index calculus
    /// does) rather than from the whole group.
    pub targets_in_subgroup: bool,
}

impl Default for PairedOptions {
    fn default() -> Self {
        Self {
            targets: 8,
            seed: 0x5EED,
            node_budget: 200_000,
            sat: true,
            conflict_budget: 2_000_000,
            max_models: 64,
            direct_x: true,
            ffd_max_degree: 4,
            targets_in_subgroup: true,
        }
    }
}

fn median(mut v: Vec<f64>) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

struct Acc {
    arm: &'static str,
    n_vars: usize,
    n_eqs: usize,
    degree: u32,
    found: usize,
    refuted: usize,
    inconclusive: usize,
    via_t: usize,
    found_ms: Vec<f64>,
    refuted_ms: Vec<f64>,
    total_ms: f64,
    effort: Vec<f64>,
    gate_failures: usize,
    ffd: Option<u32>,
}

impl Acc {
    fn new(arm: &'static str) -> Self {
        Acc {
            arm,
            n_vars: 0,
            n_eqs: 0,
            degree: 0,
            found: 0,
            refuted: 0,
            inconclusive: 0,
            via_t: 0,
            found_ms: Vec::new(),
            refuted_ms: Vec::new(),
            total_ms: 0.0,
            effort: Vec::new(),
            gate_failures: 0,
            ffd: None,
        }
    }
    /// Record one outcome against the truth `(exists, target_sum_check)`.
    fn record(&mut self, o: &OracleOutcome, ms: f64, truth_exists: bool, sum_ok: Option<bool>) {
        self.n_vars = o.n_vars;
        self.n_eqs = o.n_equations;
        self.degree = o.degree;
        self.total_ms += ms;
        self.effort.push(o.effort as f64);
        match (&o.relation, o.complete) {
            (Some(_), _) => {
                self.found += 1;
                self.found_ms.push(ms);
                if o.used_t {
                    self.via_t += 1;
                }
                if !truth_exists || sum_ok == Some(false) {
                    self.gate_failures += 1;
                }
            }
            (None, true) => {
                self.refuted += 1;
                self.refuted_ms.push(ms);
                if truth_exists {
                    self.gate_failures += 1;
                }
            }
            (None, false) => self.inconclusive += 1,
        }
    }
    fn finish(self) -> ArmRun {
        ArmRun {
            arm: self.arm,
            n_vars: self.n_vars,
            n_equations: self.n_eqs,
            degree: self.degree,
            found: self.found,
            refuted: self.refuted,
            inconclusive: self.inconclusive,
            via_t: self.via_t,
            median_found_ms: median(self.found_ms),
            median_refuted_ms: median(self.refuted_ms),
            total_ms: self.total_ms,
            median_effort: median(self.effort),
            gate_failures: self.gate_failures,
            first_fall_degree: self.ffd,
        }
    }
}

/// Run every arm on the same targets and gate every verdict.
pub fn paired_bench(a: u8, n: u32, m: usize, opts: &PairedOptions) -> Option<PairedBench> {
    let kc = KoblitzCurve::new(a, n)?;
    let target_dim = (n + 1).div_ceil(m as u32);
    let divisor = divisor_for_dimension(n, target_dim)?;
    let fb_u = build_symmetrised_factor_base(&kc, &divisor)?;
    let fb_x = build_frobenius_factor_base_from_divisor(&kc, &divisor)?;
    let index_x = fb_x.index_map();
    let st = FieldStructure::new(kc.n, &kc.curve.irreducible);

    let mut rng = StdRng::seed_from_u64(opts.seed);
    let g = kc.generator().clone();
    let r_u64 = kc
        .subgroup_order
        .to_u64_digits()
        .first()
        .copied()
        .unwrap_or(2)
        .max(2);
    let all_points: Vec<BinaryPoint> = if opts.targets_in_subgroup {
        Vec::new()
    } else {
        (0..(1u64 << n))
            .flat_map(|x| points_with_x(&kc.curve, &from_bits(x, n)))
            .collect()
    };
    let mut targets: Vec<BinaryPoint> = Vec::new();
    while targets.len() < opts.targets {
        let p = if opts.targets_in_subgroup {
            kc.mul(&g, &BigUint::from(rng.gen_range(1..r_u64)))
        } else {
            all_points[rng.gen_range(0..all_points.len())].clone()
        };
        // Both systems need u(R) finite and non-zero.
        if let BinaryPoint::Affine { x, .. } = &p {
            if *x == F2mElement::one(n) {
                continue;
            }
        } else {
            continue;
        }
        targets.push(p);
    }

    let mut x_chain = Acc::new("x-chained F4");
    let mut x_direct = Acc::new("x-direct F4");
    let mut sym = Acc::new("sym F4");
    let mut x_chain_sat = Acc::new("x-chained SAT");
    let mut sym_sat = Acc::new("sym SAT");
    let engine = SolverEngine::default();

    for (ti, target) in targets.iter().enumerate() {
        // Truth on each base.
        let truth_x = enumerate_decompose(&kc, &fb_x, &index_x, target, m).is_some();
        let truth_u = enumerate_symmetrised(&kc, &fb_u, target, m).is_some();
        let t_plus = kc.add(target, &fb_u.two_torsion);
        let sum_check_u = |o: &OracleOutcome| -> Option<bool> {
            o.relation.as_ref().map(|idx| {
                let s = fb_u.sum(&kc, idx);
                if o.used_t {
                    // T appended: the list itself sums to R.
                    s == *target || s == t_plus
                } else {
                    s == *target
                }
            })
        };
        let sum_check_x = |idx: &Option<Vec<usize>>| -> Option<bool> {
            idx.as_ref().map(|v| {
                v.iter().fold(BinaryPoint::Infinity, |acc, &i| {
                    kc.add(&acc, &fb_x.points[i])
                }) == *target
            })
        };

        // x-chained (production path).
        let t0 = Instant::now();
        let (rel, stats) = groebner_decompose(
            &kc,
            &fb_x,
            &index_x,
            &st,
            target,
            m,
            engine,
            opts.node_budget,
        );
        let ms = t0.elapsed().as_secs_f64() * 1e3;
        let x_r = match target {
            BinaryPoint::Affine { x, .. } => x.clone(),
            BinaryPoint::Infinity => unreachable!(),
        };
        let chained_sys =
            build_decomposition_system(&fb_x.subspace_basis, &x_r, &kc.curve.b, m, &st)?;
        let o = OracleOutcome {
            used_t: false,
            complete: !stats.exhausted,
            n_vars: chained_sys.n_vars,
            n_equations: chained_sys.equations.len(),
            degree: system_degree(&chained_sys.equations),
            effort: stats.splits as u64,
            relation: rel.clone(),
        };
        x_chain.record(&o, ms, truth_x, sum_check_x(&rel));
        if ti == 0 && opts.ffd_max_degree >= 2 {
            x_chain.ffd = first_fall_degree(
                &chained_sys.equations,
                chained_sys.n_vars,
                opts.ffd_max_degree,
            )
            .0;
        }

        // x-direct S₄.
        if m == 3 && opts.direct_x {
            let t0 = Instant::now();
            if let Some(o) = direct_x_groebner_decompose(
                &kc,
                &fb_x.subspace_basis,
                &fb_x.points,
                &index_x,
                &st,
                target,
                m,
                engine,
                opts.node_budget,
            ) {
                let ms = t0.elapsed().as_secs_f64() * 1e3;
                let ok = sum_check_x(&o.relation);
                x_direct.record(&o, ms, truth_x, ok);
                if ti == 0 && opts.ffd_max_degree >= 2 {
                    if let Some((eqs, nv)) =
                        build_direct_x_system(&fb_x.subspace_basis, &x_r, m, &st)
                    {
                        x_direct.ffd = first_fall_degree(&eqs, nv, opts.ffd_max_degree).0;
                    }
                }
            }
        }

        // sym F4.
        let t0 = Instant::now();
        if let Some(o) =
            symmetrised_groebner_decompose(&kc, &fb_u, &st, target, m, engine, opts.node_budget)
        {
            let ms = t0.elapsed().as_secs_f64() * 1e3;
            let ok = sum_check_u(&o);
            sym.record(&o, ms, truth_u, ok);
            if ti == 0 && opts.ffd_max_degree >= 2 {
                if let Some(s) = build_symmetrised_system(&kc, &fb_u, target, m, &st) {
                    sym.ffd = first_fall_degree(&s.equations, s.n_vars, opts.ffd_max_degree).0;
                }
            }
        }

        if opts.sat {
            let t0 = Instant::now();
            if let Some(o) = chained_x_sat_decompose(
                &kc,
                &fb_x.subspace_basis,
                &fb_x.points,
                &index_x,
                &st,
                target,
                m,
                opts.conflict_budget,
                opts.max_models,
            ) {
                let ms = t0.elapsed().as_secs_f64() * 1e3;
                let ok = sum_check_x(&o.relation);
                x_chain_sat.record(&o, ms, truth_x, ok);
            }
            let t0 = Instant::now();
            if let Some(o) = symmetrised_sat_decompose(
                &kc,
                &fb_u,
                &st,
                target,
                m,
                opts.conflict_budget,
                opts.max_models,
            ) {
                let ms = t0.elapsed().as_secs_f64() * 1e3;
                let ok = sum_check_u(&o);
                sym_sat.record(&o, ms, truth_u, ok);
            }
        }
    }

    let mut arms = vec![x_chain.finish()];
    if m == 3 && opts.direct_x {
        arms.push(x_direct.finish());
    }
    arms.push(sym.finish());
    if opts.sat {
        arms.push(x_chain_sat.finish());
        arms.push(sym_sat.finish());
    }
    Some(PairedBench {
        a,
        n,
        ell: fb_u.ell,
        m,
        divisor,
        fb_x_size: fb_x.points.len(),
        fb_u_size: fb_u.points.len(),
        targets: targets.len(),
        arms,
    })
}

/// Text table of a paired benchmark.
pub fn format_paired(b: &PairedBench) -> String {
    use std::fmt::Write as _;
    let mut s = String::new();
    let _ = writeln!(
        s,
        "== K_{}/F_2^{}  m = {}  dim V = {}  |F_x| = {}  |F_u| = {}  targets = {}",
        b.a, b.n, b.m, b.ell, b.fb_x_size, b.fb_u_size, b.targets
    );
    let _ = writeln!(
        s,
        "   {:<14} {:>4} {:>4} {:>3} {:>5} {:>7} {:>6} {:>5} {:>12} {:>12} {:>9} {:>4} {:>5}",
        "arm",
        "vars",
        "eqs",
        "deg",
        "found",
        "refuted",
        "inconc",
        "via T",
        "found ms",
        "refuted ms",
        "effort",
        "ffd",
        "gate"
    );
    let fmt = |v: f64| {
        if v.is_nan() {
            "-".to_string()
        } else if v >= 100.0 {
            format!("{v:.0}")
        } else {
            format!("{v:.2}")
        }
    };
    for a in &b.arms {
        let _ = writeln!(
            s,
            "   {:<14} {:>4} {:>4} {:>3} {:>5} {:>7} {:>6} {:>5} {:>12} {:>12} {:>9} {:>4} {:>5}",
            a.arm,
            a.n_vars,
            a.n_equations,
            a.degree,
            a.found,
            a.refuted,
            a.inconclusive,
            a.via_t,
            fmt(a.median_found_ms),
            fmt(a.median_refuted_ms),
            fmt(a.median_effort),
            a.first_fall_degree
                .map(|d| d.to_string())
                .unwrap_or_else(|| "-".into()),
            if a.gate_failures == 0 {
                "ok".to_string()
            } else {
                format!("FAIL {}", a.gate_failures)
            }
        );
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::coordinate_search::{
        detect_symmetries, interpolate_relation, linearising_frame, CoordinateSystem, Curve, Frame,
        Rng64, SymmetryKind,
    };

    fn interpolated_terms(m: usize, symmetrised: bool) -> Vec<Vec<u32>> {
        let k = Curve::koblitz(1, 9);
        let pts = k.affine_points();
        let mut rng = Rng64::new(2);
        let cs = if symmetrised {
            let g = detect_symmetries(&k, &pts, &mut rng)
                .into_iter()
                .find_map(|s| match s.kind {
                    SymmetryKind::TwoTorsion(_) => s.mobius,
                    _ => None,
                })
                .unwrap();
            let frame = linearising_frame(&k.f, &g);
            CoordinateSystem::with_involution(&k.f, &frame, &g, false)
        } else {
            CoordinateSystem::plain(&Frame::weierstrass(), false)
        };
        let r = interpolate_relation(&k, &pts, &cs, m, &mut rng).unwrap();
        let mut terms: Vec<Vec<u32>> = r
            .poly
            .terms
            .iter()
            .map(|(e, c)| {
                assert_eq!(*c, 1, "coefficients must be in F_2");
                e.clone()
            })
            .collect();
        terms.sort();
        terms
    }

    #[test]
    fn hardcoded_polynomials_match_the_interpolation() {
        for m in [2usize, 3] {
            let mut sym = symmetrised_terms(m).unwrap();
            sym.sort();
            assert_eq!(sym, interpolated_terms(m, true), "symmetrised m = {m}");
            let mut plain = plain_terms(m).unwrap();
            plain.sort();
            assert_eq!(plain, interpolated_terms(m, false), "plain m = {m}");
        }
    }

    #[test]
    fn u_factor_base_is_closed_under_t_and_frobenius() {
        for (a, n) in [(0u8, 9u32), (1, 9)] {
            let kc = KoblitzCurve::new(a, n).unwrap();
            let div = divisor_for_dimension(n, 5).unwrap();
            let fb = build_symmetrised_factor_base(&kc, &div).unwrap();
            assert_eq!(fb.v_basis[0], F2mElement::one(n));
            assert!(fb.index_of.contains_key(&point_key(&fb.two_torsion)));
            for p in &fb.points {
                let q = kc.add(p, &fb.two_torsion);
                if q != BinaryPoint::Infinity {
                    assert!(fb.index_of.contains_key(&point_key(&q)), "not T-closed");
                }
                assert!(
                    fb.index_of.contains_key(&point_key(&kc.frobenius(p))),
                    "not π-closed"
                );
                assert!(
                    fb.index_of.contains_key(&point_key(&point_neg(p))),
                    "not −-closed"
                );
            }
        }
    }

    #[test]
    fn divisor_without_x_minus_1_is_rejected() {
        let kc = KoblitzCurve::new(1, 9).unwrap();
        let factors = all_factors_of_x_n_minus_1(9);
        let other = (0..factors.len()).find(|&i| factors[i] != 0b11).unwrap();
        assert!(build_symmetrised_factor_base(&kc, &[other]).is_none());
    }

    #[test]
    fn direct_x_s3_matches_the_production_s3_system() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let div = divisor_for_dimension(9, 5).unwrap();
        let fb_x = build_frobenius_factor_base_from_divisor(&kc, &div).unwrap();
        let st = FieldStructure::new(9, &kc.curve.irreducible);
        let x_r = from_bits(0x1a3, 9);
        let (direct, nv) = build_direct_x_system(&fb_x.subspace_basis, &x_r, 2, &st).unwrap();
        let chained =
            build_decomposition_system(&fb_x.subspace_basis, &x_r, &kc.curve.b, 2, &st).unwrap();
        assert_eq!(nv, chained.n_vars);
        for root in 0..(1u64 << nv) {
            for (p, q) in direct.iter().zip(&chained.equations) {
                assert_eq!(p.eval(root), q.eval(root));
            }
        }
    }

    fn agreement(a: u8, n: u32, m: usize, sat: bool) {
        let opts = PairedOptions {
            targets: 6,
            seed: 7,
            node_budget: 50_000,
            sat,
            conflict_budget: 500_000,
            max_models: 32,
            direct_x: true,
            ffd_max_degree: 0,
            targets_in_subgroup: false,
        };
        let b = paired_bench(a, n, m, &opts).unwrap();
        for arm in &b.arms {
            assert_eq!(arm.gate_failures, 0, "{}\n{}", arm.arm, format_paired(&b));
            assert_eq!(arm.inconclusive, 0, "{}\n{}", arm.arm, format_paired(&b));
        }
    }

    #[test]
    fn symmetrised_oracles_agree_with_enumeration_m2() {
        agreement(0, 9, 2, true);
        agreement(1, 9, 2, true);
    }

    #[test]
    fn symmetrised_oracles_agree_with_enumeration_m3() {
        agreement(0, 9, 3, true);
    }

    #[test]
    fn symmetrised_system_has_the_predicted_shape() {
        let kc = KoblitzCurve::new(1, 15).unwrap();
        let st = FieldStructure::new(15, &kc.curve.irreducible);
        let div = divisor_for_dimension(15, 8).unwrap();
        let fb = build_symmetrised_factor_base(&kc, &div).unwrap();
        let target = kc.generator().clone();
        let s2 = build_symmetrised_system(&kc, &fb, &target, 2, &st).unwrap();
        assert_eq!(s2.n_vars, 2 * (fb.ell - 1) + 1);
        assert_eq!(s2.equations.len(), 15);
        assert_eq!(system_degree(&s2.equations), 2);
        let s3 = build_symmetrised_system(&kc, &fb, &target, 3, &st).unwrap();
        assert_eq!(s3.n_vars, 3 * (fb.ell - 1) + 1);
        assert_eq!(system_degree(&s3.equations), 4);
    }
}
