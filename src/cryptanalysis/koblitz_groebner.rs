//! # Semaev polynomial systems for Koblitz-curve index calculus, solved by Gröbner basis.
//!
//! The decomposition oracle behind
//! [`crate::cryptanalysis::koblitz_index_calculus`].  Its question is
//!
//! > given `R ∈ E(F_{2^n})`, is `R = P_1 + … + P_m` with every `P_i` in
//! > the Frobenius-invariant factor base `F`?
//!
//! Galbraith–Granger–Merz–Petit answer it the way every index calculus
//! on elliptic curves since Gaudry and Diem does: write the condition
//! as a **Semaev summation polynomial**, restrict the unknowns to the
//! factor base (here, to an `F_2`-subspace `V ⊆ F_{2^n}`), take the Weil
//! restriction of the single `F_{2^n}`-equation down to `n` equations
//! over `F_2`, and hand the resulting Boolean system to a Gröbner-basis
//! engine.  This module is that pipeline; the caller keeps its
//! exhaustive-search oracle as a cross-check.
//!
//! ## The system
//!
//! For `m = 2` and `E : y² + xy = x³ + a x² + b`, the third summation
//! polynomial is (see [`crate::cryptanalysis::binary_semaev`])
//!
//! ```text
//!     S₃(x₁, x₂, x₃) = (x₁+x₂)² x₃² + x₁x₂ x₃ + (x₁x₂)² + b,
//! ```
//!
//! which vanishes exactly when some choice of `y`-coordinates makes
//! `±P₁ ± P₂ ± P₃ = O`.  Setting `x₃ := x(R)` and requiring
//! `x₁, x₂ ∈ V` gives one equation in two `V`-valued unknowns.  For
//! `m ≥ 3` the polynomial is **chained** rather than resolved — `S_{k}`
//! has degree `2^{k−2}` per variable and is impractical past `k = 3` —
//! so `m` summands become `m − 1` links of `S₃` joined by `m − 2`
//! intermediate unknowns ranging over the whole field:
//!
//! ```text
//!     S₃(x₁, x₂, e₁) = S₃(e₁, x₃, e₂) = … = S₃(e_{m−2}, x_m, x(R)) = 0.
//! ```
//!
//! ## Why the system stays low-degree
//!
//! Writing `x_i = Σ_t u_{i,t} · b_t` over an `F_2`-basis `b_t` of `V`
//! makes each `x_i` **linear** in the Boolean unknowns.  Squaring is
//! `F_2`-linear in characteristic 2 — and in the Boolean quotient
//! `v² = v` a polynomial satisfies `p² = p`, so a symbolic square costs
//! a permutation of coordinates and nothing else — while `x₁x₂` is
//! bilinear.
//!
//! For `m = 2` the third argument of `S₃` is the *known* `x(R)`, so
//! every term has degree ≤ 2: the Weil restriction is a **quadratic**
//! system of `n` equations in `2ℓ` unknowns — the shape Faugère et al.
//! analyse.  For `m ≥ 3` the chain's intermediate points are unknowns
//! too, and the term `x₁x₂·e` makes the system **cubic** in
//! `m·ℓ + (m−2)·n` unknowns.  Degree 3 is why the `m ≥ 3` instances
//! cost so much more here, in both engines.
//!
//! Restricting the unknowns to `V` is what keeps the degree this low at
//! all, and it is the property the linearised-polynomial factor base
//! has and a weight-bounded one does not.
//!
//! ## What the Gröbner step buys
//!
//! The exhaustive oracle costs `|F|^{m−1}` group operations per target
//! whether or not a decomposition exists.  Here, an **unsolvable**
//! target is usually rejected by the Gröbner computation alone: the
//! basis reduces to `{1}`, no search happens.  Solvable targets are
//! read off the basis by splitting (below).  The cost is governed by
//! the degree of regularity of the system rather than by `|F|`, which
//! is what makes the real algorithm sub-exponential.
//!
//! ## Solving the Boolean system
//!
//! [`groebner_basis_f2`] returns a DegRevLex basis, not a triangular
//! one, and this module deliberately does **not** call the brute-force
//! `solve_system_f2` extractor that would defeat the point.  Instead
//! [`solve_boolean_system`] runs the standard *Gröbner-with-splitting*
//! (`DPLL(GB)`) loop: compute a basis; if it contains the constant `1`
//! the branch is infeasible; propagate every basis element that has
//! collapsed to `v_i` or `v_i + 1`; otherwise split on the lowest free
//! variable and recurse.  Every branch is bounded by a node budget, and
//! solutions are verified against the original system before they are
//! returned.
//!
//! ## Honest scope
//!
//! - The engine underneath is Buchberger with Gebauer–Möller pruning
//!   ([`crate::cryptanalysis::pq_groebner_f2`]), not F4/F5 with sparse
//!   linear algebra.  It is the right *algorithm* and the wrong
//!   *constant*; the toy parameters (`m·ℓ + (m−2)·n ≤ 64` unknowns, and
//!   in practice far fewer) are chosen accordingly.
//! - `S₃` is a *necessary* condition on `x`-coordinates only.  A root
//!   fixes the summands up to sign, so each candidate is lifted to
//!   actual points and the group identity `P_1 + … + P_m = R` is
//!   re-checked before a relation is emitted.  No false relation can
//!   escape.
//!
//! ## References
//!
//! - S. D. Galbraith, R. Granger, S.-P. Merz, C. Petit, *On index
//!   calculus algorithms for subfield curves*, SAC 2020.
//! - I. Semaev, *Summation polynomials and the discrete logarithm
//!   problem on elliptic curves*, ePrint 2004/031.
//! - C. Diem, *On the discrete logarithm problem in elliptic curves*,
//!   Compositio Math. 147 (2011).
//! - J.-C. Faugère, L. Perret, C. Petit, G. Renault, *Improving the
//!   complexity of index calculus algorithms in elliptic curves over
//!   binary fields*, EUROCRYPT 2012 — the Weil-restriction-to-Boolean
//!   -system analysis this module implements.
//! - G. Bard, *Algebraic Cryptanalysis*, Springer 2009, ch. 13 — the
//!   Boolean-ring representation.

use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use crate::cryptanalysis::pq_groebner_f2::{cmp_mono, groebner_basis_f2, F2BoolMono, F2BoolPoly};

/// Hard cap: Boolean monomials are `u64` bitmasks in
/// [`crate::cryptanalysis::pq_groebner_f2`].
pub const MAX_VARS: usize = 64;

// ── Symbolic F_{2^n} arithmetic ────────────────────────────────────

/// Multiplication table of the polynomial basis `1, z, …, z^{n−1}`:
/// `reduced[i][j]` is the bit pattern of `z^i · z^j mod f(z)`.
///
/// This is all a symbolic multiplication needs: the coordinates of a
/// product are `F_2`-bilinear in the coordinates of the operands, with
/// these structure constants.
#[derive(Clone, Debug)]
pub struct FieldStructure {
    /// Extension degree.
    pub n: u32,
    /// `reduced[i][j] = z^{i+j} mod f`, as a bitmask.
    pub reduced: Vec<Vec<u64>>,
    /// `squares[k] = z^{2k} mod f`, as a bitmask — the coordinate
    /// permutation-with-reduction that implements symbolic squaring.
    pub squares: Vec<u64>,
}

impl FieldStructure {
    /// Tabulate the structure constants of `F_2[z]/(f)`.
    pub fn new(n: u32, irr: &IrreduciblePoly) -> Self {
        let bits = |e: &F2mElement| e.raw_bits().first().copied().unwrap_or(0);
        let mono = |k: u32| F2mElement::from_bit_positions(&[k], n);
        let mut reduced = vec![vec![0u64; n as usize]; n as usize];
        for i in 0..n as usize {
            for j in 0..n as usize {
                reduced[i][j] = bits(&mono(i as u32).mul(&mono(j as u32), irr));
            }
        }
        let squares = (0..n).map(|k| reduced[k as usize][k as usize]).collect();
        Self {
            n,
            reduced,
            squares,
        }
    }
}

/// An element of `F_{2^n}` whose coordinates are Boolean polynomials —
/// i.e. a symbolic field element in the Weil restriction.
#[derive(Clone, Debug)]
pub struct SymElement {
    /// `coords[k]` multiplies `z^k`.  Length `n`.
    pub coords: Vec<F2BoolPoly>,
}

/// Product of two Boolean polynomials (the engine ships `mul_mono`
/// only; `p·q` is the sum of `p·m` over the monomials `m` of `q`).
pub fn poly_mul(p: &F2BoolPoly, q: &F2BoolPoly) -> F2BoolPoly {
    if p.is_zero() || q.is_zero() {
        return F2BoolPoly::zero(p.n_vars);
    }
    let (p, q) = if p.terms.len() <= q.terms.len() {
        (p, q)
    } else {
        (q, p)
    };
    let mut acc = F2BoolPoly::zero(p.n_vars);
    for m in &q.terms {
        acc = acc.add(&p.mul_mono(*m));
    }
    acc
}

impl SymElement {
    /// The zero element.
    pub fn zero(n: u32, n_vars: usize) -> Self {
        Self {
            coords: (0..n).map(|_| F2BoolPoly::zero(n_vars)).collect(),
        }
    }

    /// A known field element, lifted to a constant symbolic one.
    pub fn constant(v: &F2mElement, n: u32, n_vars: usize) -> Self {
        let bits = v.raw_bits().first().copied().unwrap_or(0);
        let coords = (0..n)
            .map(|k| {
                if (bits >> k) & 1 == 1 {
                    F2BoolPoly::one(n_vars)
                } else {
                    F2BoolPoly::zero(n_vars)
                }
            })
            .collect();
        Self { coords }
    }

    /// An unknown ranging over the `F_2`-span of `basis`, using the
    /// Boolean variables `offset … offset + basis.len() − 1`:
    /// `x = Σ_t v_{offset+t} · basis[t]`.
    ///
    /// This is where the factor base enters the algebra — restricting
    /// the unknowns to the Frobenius-invariant subspace `V` is exactly
    /// what makes the system solvable at all.
    pub fn from_subspace_vars(basis: &[F2mElement], offset: usize, n: u32, n_vars: usize) -> Self {
        let mut coords = vec![Vec::new(); n as usize];
        for (t, b) in basis.iter().enumerate() {
            let bits = b.raw_bits().first().copied().unwrap_or(0);
            for k in 0..n as usize {
                if (bits >> k) & 1 == 1 {
                    coords[k].push(F2BoolMono::var((offset + t) as u32));
                }
            }
        }
        Self {
            coords: coords
                .into_iter()
                .map(|ms| F2BoolPoly::from_monos(ms, n_vars))
                .collect(),
        }
    }

    /// A wholly unknown field element on variables `offset … offset+n−1`
    /// — used for the intermediate points of a chained `S₃`.
    pub fn from_free_vars(offset: usize, n: u32, n_vars: usize) -> Self {
        Self {
            coords: (0..n as usize)
                .map(|k| F2BoolPoly::from_monos(vec![F2BoolMono::var((offset + k) as u32)], n_vars))
                .collect(),
        }
    }

    /// Coordinate-wise sum (char 2: also the difference).
    pub fn add(&self, other: &Self) -> Self {
        Self {
            coords: self
                .coords
                .iter()
                .zip(&other.coords)
                .map(|(a, b)| a.add(b))
                .collect(),
        }
    }

    /// Field multiplication via the structure constants.
    pub fn mul(&self, other: &Self, st: &FieldStructure) -> Self {
        let n = st.n as usize;
        let n_vars = self.coords[0].n_vars;
        // Cross products a_i · b_j, computed once.
        let mut prod = vec![vec![F2BoolPoly::zero(n_vars); n]; n];
        for i in 0..n {
            if self.coords[i].is_zero() {
                continue;
            }
            for j in 0..n {
                if other.coords[j].is_zero() {
                    continue;
                }
                prod[i][j] = poly_mul(&self.coords[i], &other.coords[j]);
            }
        }
        let mut coords = vec![F2BoolPoly::zero(n_vars); n];
        for i in 0..n {
            for j in 0..n {
                if prod[i][j].is_zero() {
                    continue;
                }
                let red = st.reduced[i][j];
                for (k, coord) in coords.iter_mut().enumerate() {
                    if (red >> k) & 1 == 1 {
                        *coord = coord.add(&prod[i][j]);
                    }
                }
            }
        }
        Self { coords }
    }

    /// Field squaring — an `F_2`-linear map, so it never raises the
    /// degree of the Boolean system.
    ///
    /// In characteristic 2, `(Σ_k a_k z^k)² = Σ_k a_k² z^{2k}`, and in
    /// the Boolean quotient `p² = p` for every polynomial `p` (the
    /// cross terms cancel in pairs and `v² = v`).  So squaring just
    /// re-indexes coordinates through `z^{2k} mod f`.
    pub fn square(&self, st: &FieldStructure) -> Self {
        let n = st.n as usize;
        let n_vars = self.coords[0].n_vars;
        let mut coords = vec![F2BoolPoly::zero(n_vars); n];
        for (k, c) in self.coords.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let sq = st.squares[k];
            for (t, coord) in coords.iter_mut().enumerate() {
                if (sq >> t) & 1 == 1 {
                    *coord = coord.add(c);
                }
            }
        }
        Self { coords }
    }

    /// Evaluate at a Boolean assignment, giving a concrete field element.
    pub fn eval(&self, point: u64, n: u32) -> F2mElement {
        let mut positions = Vec::new();
        for (k, c) in self.coords.iter().enumerate() {
            if c.eval(point) == 1 {
                positions.push(k as u32);
            }
        }
        F2mElement::from_bit_positions(&positions, n)
    }
}

// ── The Semaev system ──────────────────────────────────────────────

/// Symbolic `S₃(x₁, x₂, x₃) = (x₁+x₂)² x₃² + x₁x₂ x₃ + (x₁x₂)² + b`.
///
/// Returns the `n` Boolean coordinate equations of the single
/// `F_{2^n}`-equation `S₃ = 0`.
pub fn sym_semaev_s3(
    x1: &SymElement,
    x2: &SymElement,
    x3: &SymElement,
    b: &F2mElement,
    st: &FieldStructure,
) -> Vec<F2BoolPoly> {
    let n_vars = x1.coords[0].n_vars;
    let sum_sq = x1.add(x2).square(st);
    let x3_sq = x3.square(st);
    let prod = x1.mul(x2, st);

    let mut acc = sum_sq.mul(&x3_sq, st);
    acc = acc.add(&prod.mul(x3, st));
    acc = acc.add(&prod.square(st));
    acc = acc.add(&SymElement::constant(b, st.n, n_vars));
    acc.coords
}

/// The Boolean system whose roots are the `m`-point decompositions of a
/// target with abscissa `x_r` over the subspace spanned by `basis`.
#[derive(Clone, Debug)]
pub struct DecompositionSystem {
    /// The equations, `n` per `S₃` link.
    pub equations: Vec<F2BoolPoly>,
    /// Total Boolean unknowns.
    pub n_vars: usize,
    /// Unknowns per summand (`ℓ = dim V`).
    pub ell: usize,
    /// Number of summands.
    pub m: usize,
}

impl DecompositionSystem {
    /// The `x`-coordinate that Boolean assignment `point` gives to
    /// summand `i`, as an element of the subspace.
    pub fn summand_x(&self, basis: &[F2mElement], point: u64, i: usize, n: u32) -> F2mElement {
        let mut positions = Vec::new();
        for t in 0..self.ell {
            if (point >> (i * self.ell + t)) & 1 == 1 {
                let bits = basis[t].raw_bits().first().copied().unwrap_or(0);
                for k in 0..n {
                    if (bits >> k) & 1 == 1 {
                        positions.push(k);
                    }
                }
            }
        }
        // XOR-fold: a coordinate hit twice cancels.
        let mut acc = F2mElement::zero(n);
        for k in positions {
            acc = acc.add(&F2mElement::from_bit_positions(&[k], n));
        }
        acc
    }
}

/// **Build the decomposition system** for `R = P_1 + … + P_m`, all
/// `x(P_i)` constrained to the `F_2`-span of `basis`.
///
/// Variable layout: `m · ℓ` subspace coordinates first (summand `i`
/// owns `[i·ℓ, (i+1)·ℓ)`), then `(m − 2)·n` coordinates for the
/// intermediate points of the chain.
///
/// Returns `None` if the layout would exceed [`MAX_VARS`] or if
/// `m < 2`.
pub fn build_decomposition_system(
    basis: &[F2mElement],
    x_r: &F2mElement,
    b: &F2mElement,
    m: usize,
    st: &FieldStructure,
) -> Option<DecompositionSystem> {
    if m < 2 {
        return None;
    }
    let ell = basis.len();
    let n = st.n;
    let n_vars = m * ell + m.saturating_sub(2) * n as usize;
    if n_vars > MAX_VARS {
        return None;
    }

    let summands: Vec<SymElement> = (0..m)
        .map(|i| SymElement::from_subspace_vars(basis, i * ell, n, n_vars))
        .collect();
    let target = SymElement::constant(x_r, n, n_vars);

    let mut equations = Vec::new();
    if m == 2 {
        equations.extend(sym_semaev_s3(&summands[0], &summands[1], &target, b, st));
    } else {
        // Chain: S₃(x₁, x₂, e₁), S₃(e_i, x_{i+2}, e_{i+1}), …, last
        // link closes on x(R).
        let inter: Vec<SymElement> = (0..m - 2)
            .map(|i| SymElement::from_free_vars(m * ell + i * n as usize, n, n_vars))
            .collect();
        equations.extend(sym_semaev_s3(&summands[0], &summands[1], &inter[0], b, st));
        for i in 0..m - 3 {
            equations.extend(sym_semaev_s3(
                &inter[i],
                &summands[i + 2],
                &inter[i + 1],
                b,
                st,
            ));
        }
        equations.extend(sym_semaev_s3(
            &inter[m - 3],
            &summands[m - 1],
            &target,
            b,
            st,
        ));
    }

    Some(DecompositionSystem {
        equations,
        n_vars,
        ell,
        m,
    })
}

// ── Matrix-F4 over the Boolean ring ────────────────────────────────

/// All monomials of degree `≤ deg` over `n_vars` Boolean variables, as
/// masks.  (In the Boolean ring a monomial is a subset of variables.)
fn monomials_up_to(n_vars: usize, deg: u32) -> Vec<u64> {
    let mut out = vec![0u64];
    let mut level = vec![0u64];
    for _ in 0..deg {
        let mut next = Vec::new();
        for &m in &level {
            let top = 64 - m.leading_zeros();
            for v in top..n_vars as u32 {
                next.push(m | (1u64 << v));
            }
        }
        out.extend(next.iter().copied());
        level = next;
        if level.is_empty() {
            break;
        }
    }
    out
}

/// Limits on one Macaulay matrix, so a too-large system degrades to
/// splitting instead of exhausting memory.
const MAX_F4_ROWS: usize = 20_000;
const MAX_F4_COLS: usize = 40_000;

/// Macaulay size caps, overridable for experiments through the
/// `F4_F2_MAX_ROWS` / `F4_F2_MAX_COLS` environment variables.
fn max_f4_rows() -> usize {
    std::env::var("F4_F2_MAX_ROWS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(MAX_F4_ROWS)
}
fn max_f4_cols() -> usize {
    std::env::var("F4_F2_MAX_COLS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(MAX_F4_COLS)
}

/// **Matrix-F4 step over `F_2`**: multiply every input polynomial by
/// every monomial that keeps the degree `≤ degree`, row-reduce the
/// resulting Macaulay matrix, and return the reduced rows as
/// polynomials.
///
/// This is F4's engine — Buchberger's S-pair reductions replaced by one
/// linear-algebra pass over all products at a bounded degree — without
/// the symbolic preprocessing that chooses a minimal row set.  Every
/// returned polynomial lies in the ideal generated by `polys` (it is an
/// `F_2`-combination of multiples of them), so a returned constant `1`
/// is a genuine certificate of infeasibility.
///
/// Returns `None` if the matrix would exceed the size limits.
pub fn matrix_f4_f2(polys: &[F2BoolPoly], n_vars: usize, degree: u32) -> Option<Vec<F2BoolPoly>> {
    if polys.is_empty() {
        return Some(Vec::new());
    }
    let (cols, mut matrix) = build_macaulay(polys, n_vars, degree)?;
    if matrix.is_empty() {
        return Some(Vec::new());
    }
    let rank = rref_f2(&mut matrix, cols.len());

    let n_vars_out = polys[0].n_vars;
    let words = cols.len().div_ceil(64);
    let _ = words;
    let mut out = Vec::with_capacity(rank);
    for row in matrix.iter().take(rank) {
        let monos: Vec<F2BoolMono> = (0..cols.len())
            .filter(|c| row[c / 64] & (1u64 << (c % 64)) != 0)
            .map(|c| F2BoolMono::from_mask(cols[c]))
            .collect();
        if !monos.is_empty() {
            out.push(F2BoolPoly::from_monos(monos, n_vars_out));
        }
    }
    Some(out)
}

/// Build the Macaulay matrix: every product `p · m` with
/// `deg(p·m) ≤ degree`, as bit-rows over the monomials that occur.
///
/// Returns the column monomials (DegRevLex descending) and the rows.
/// `None` if the matrix would exceed the size limits.
fn build_macaulay(
    polys: &[F2BoolPoly],
    n_vars: usize,
    degree: u32,
) -> Option<(Vec<u64>, Vec<Vec<u64>>)> {
    let mut rows_monos: Vec<Vec<u64>> = Vec::new();
    for p in polys {
        let pdeg = p
            .terms
            .iter()
            .map(|t| t.mask.count_ones())
            .max()
            .unwrap_or(0);
        if pdeg > degree {
            continue;
        }
        for mult in monomials_up_to(n_vars, degree - pdeg) {
            // Multiplying by a monomial is a union of masks, so two
            // distinct terms of `p` can collide — and collide means
            // cancel, in characteristic 2.  Keep the odd multiplicities.
            let mut all: Vec<u64> = p.terms.iter().map(|t| t.mask | mult).collect();
            all.sort_unstable();
            let mut row: Vec<u64> = Vec::with_capacity(all.len());
            let mut i = 0;
            while i < all.len() {
                let mut j = i;
                while j < all.len() && all[j] == all[i] {
                    j += 1;
                }
                if (j - i) % 2 == 1 {
                    row.push(all[i]);
                }
                i = j;
            }
            if !row.is_empty() {
                rows_monos.push(row);
            }
            if rows_monos.len() > max_f4_rows() {
                return None;
            }
        }
    }
    if rows_monos.is_empty() {
        return Some((Vec::new(), Vec::new()));
    }

    let mut cols: Vec<u64> = rows_monos.iter().flatten().copied().collect();
    cols.sort_unstable();
    cols.dedup();
    if cols.len() > max_f4_cols() {
        return None;
    }
    cols.sort_by(|a, b| cmp_mono(F2BoolMono::from_mask(*a), F2BoolMono::from_mask(*b)).reverse());
    let index: std::collections::HashMap<u64, usize> =
        cols.iter().enumerate().map(|(i, m)| (*m, i)).collect();

    let words = cols.len().div_ceil(64);
    let matrix: Vec<Vec<u64>> = rows_monos
        .iter()
        .map(|monos| {
            let mut row = vec![0u64; words];
            for m in monos {
                let c = index[m];
                row[c / 64] |= 1 << (c % 64);
            }
            row
        })
        .collect();
    Some((cols, matrix))
}

/// Reduced row echelon form over `F_2`; returns the rank, with the
/// pivot rows moved to the front of `matrix`.
fn rref_f2(matrix: &mut [Vec<u64>], n_cols: usize) -> usize {
    let words = n_cols.div_ceil(64);
    let mut pivot_row = 0usize;
    for c in 0..n_cols {
        let (w, bit) = (c / 64, 1u64 << (c % 64));
        let piv = (pivot_row..matrix.len()).find(|&r| matrix[r][w] & bit != 0);
        let piv = match piv {
            Some(p) => p,
            None => continue,
        };
        matrix.swap(pivot_row, piv);
        for r in 0..matrix.len() {
            if r != pivot_row && matrix[r][w] & bit != 0 {
                for k in 0..words {
                    matrix[r][k] ^= matrix[pivot_row][k];
                }
            }
        }
        pivot_row += 1;
        if pivot_row == matrix.len() {
            break;
        }
    }
    pivot_row
}

// ── Macaulay profile / first fall degree ───────────────────────────

/// Rank measurement of one Macaulay matrix.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MacaulayProfile {
    /// Degree the matrix was built at.
    pub degree: u32,
    /// Rows actually constructed.
    pub rows: usize,
    /// Distinct monomials occurring, i.e. columns.
    pub cols: usize,
    /// Rank over `F_2`.
    pub rank: usize,
}

impl MacaulayProfile {
    /// `rows − rank`: the number of independent syzygies among the
    /// Macaulay-shifted equations at this degree.
    pub fn syzygies(&self) -> usize {
        self.rows.saturating_sub(self.rank)
    }
}

/// Rank profile of the Macaulay matrix of `polys` at one degree.
pub fn macaulay_profile(
    polys: &[F2BoolPoly],
    n_vars: usize,
    degree: u32,
) -> Option<MacaulayProfile> {
    let (cols, mut matrix) = build_macaulay(polys, n_vars, degree)?;
    let rows = matrix.len();
    let rank = if matrix.is_empty() {
        0
    } else {
        rref_f2(&mut matrix, cols.len())
    };
    Some(MacaulayProfile {
        degree,
        rows,
        cols: cols.len(),
        rank,
    })
}

/// **First fall degree** of `polys`: the smallest `D ≥ 2` whose Macaulay
/// matrix has `rank < rows` *and* `rank < cols` — a non-trivial syzygy
/// appears and the system has not saturated.
///
/// This is the same operational definition
/// [`crate::cryptanalysis::ffd_harness`] uses for the full-field
/// Weil descent of `S₃`, so the numbers are directly comparable: that
/// harness measures the `2n`-variable system, this one the system
/// restricted to a Frobenius-invariant subspace, which is the version
/// the subfield-curve attack actually solves.
///
/// Returns the fall degree (if any up to `d_max`) and the per-degree
/// profiles.  A profile is omitted for degrees whose matrix exceeded
/// the size limits.
pub fn first_fall_degree(
    polys: &[F2BoolPoly],
    n_vars: usize,
    d_max: u32,
) -> (Option<u32>, Vec<MacaulayProfile>) {
    let mut fall = None;
    let mut profiles = Vec::new();
    for d in 2..=d_max {
        let prof = match macaulay_profile(polys, n_vars, d) {
            Some(p) => p,
            None => break,
        };
        if fall.is_none() && prof.rank < prof.rows && prof.rank < prof.cols {
            fall = Some(d);
        }
        profiles.push(prof);
    }
    (fall, profiles)
}

// ── Gröbner solve with splitting ───────────────────────────────────

/// Specialise `p` by setting variable `var` to `value`.
fn substitute(p: &F2BoolPoly, var: u32, value: bool) -> F2BoolPoly {
    let bit = 1u64 << var;
    let mut monos = Vec::with_capacity(p.terms.len());
    for t in &p.terms {
        if t.mask & bit == 0 {
            monos.push(*t);
        } else if value {
            monos.push(F2BoolMono::from_mask(t.mask & !bit));
        }
        // v = 0 kills every term containing v.
    }
    F2BoolPoly::from_monos(monos, p.n_vars)
}

/// A basis element that has collapsed to `v_i` or `v_i + 1` forces its
/// variable; return `(var, value)` if `p` has that shape.
fn forced_assignment(p: &F2BoolPoly) -> Option<(u32, bool)> {
    let (mut var, mut has_const) = (None, false);
    for t in &p.terms {
        if t.mask == 0 {
            has_const = true;
        } else if t.mask.count_ones() == 1 && var.is_none() {
            var = Some(t.mask.trailing_zeros());
        } else {
            return None;
        }
    }
    // p = v  ⇒ v = 0;  p = v + 1 ⇒ v = 1.
    var.map(|v| (v, has_const))
}

/// Which algebraic engine reduces the system at each node.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolverEngine {
    /// Boolean matrix-F4 ([`matrix_f4_f2`]) at degrees `2 ..= max_degree`.
    /// Fast, and the engine the real algorithm uses.
    MatrixF4 {
        /// Highest Macaulay degree to build before splitting.
        max_degree: u32,
    },
    /// Full Buchberger Gröbner basis
    /// ([`crate::cryptanalysis::pq_groebner_f2::groebner_basis_f2`]) at
    /// every node.  The reference engine: same answers, far slower on
    /// dense systems, kept so the fast path can be tested against a
    /// textbook one.
    Buchberger,
}

impl Default for SolverEngine {
    fn default() -> Self {
        SolverEngine::MatrixF4 { max_degree: 3 }
    }
}

/// Knobs for [`solve_boolean_system`].
#[derive(Clone, Copy, Debug)]
pub struct SolveOptions {
    /// Reduction engine.
    pub engine: SolverEngine,
    /// Stop after this many solutions.
    pub max_solutions: usize,
    /// Cap on algebraic reductions (one per splitting node, plus one
    /// per propagation round).
    pub node_budget: usize,
    /// How the splitter picks its variable when the algebra stalls.
    pub split_rule: SplitRule,
}

/// Which free variable the splitter branches on.
///
/// Every rule picks *some* unassigned variable, so the search stays
/// exhaustive whichever is chosen; they differ only in how quickly the
/// branch closes.  `LowestFree` is the historical behaviour and stays
/// the default so existing measurements keep their meaning.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SplitRule {
    /// Lowest-indexed unassigned variable.  On the descent systems the
    /// variable index runs point by point, so this fixes one summand's
    /// bits before touching the next.
    #[default]
    LowestFree,
    /// The free variable occurring in the most monomials of the current
    /// system.  Substituting it removes the most terms.
    MostFrequent,
    /// Occurrences weighted by `2^{1-d}` for a degree-`d` monomial, so a
    /// variable sitting in short monomials outranks one buried in long
    /// ones (the classic MOM rule).  A short monomial is closer to
    /// forcing an assignment, so this favours propagation over bulk
    /// term removal.
    MinTermWeight,
}

/// The splitting rule the decomposition entry points use when their
/// caller does not build [`SolveOptions`] itself, overridable through
/// the `SOLVER_SPLIT_RULE` environment variable (`lowest`, `frequent`,
/// `mom`).  Unset means [`SplitRule::LowestFree`], so the production
/// paths and every earlier measurement are unchanged by default.
pub fn split_rule_default() -> SplitRule {
    match std::env::var("SOLVER_SPLIT_RULE").ok().as_deref() {
        Some("frequent") => SplitRule::MostFrequent,
        Some("mom") => SplitRule::MinTermWeight,
        _ => SplitRule::LowestFree,
    }
}

/// Pick the variable to split on among those still unassigned.
///
/// Returns `None` only when every variable is assigned.  A free
/// variable that no longer occurs in the system is a don't-care: the
/// occurrence-based rules score it zero, and the lowest such variable
/// is taken when nothing scores higher.
fn choose_split(
    system: &[F2BoolPoly],
    assignment: &[Option<bool>],
    rule: SplitRule,
) -> Option<usize> {
    let lowest = assignment.iter().position(|a| a.is_none())?;
    if rule == SplitRule::LowestFree {
        return Some(lowest);
    }
    let mut score = vec![0f64; assignment.len()];
    for p in system {
        for t in &p.terms {
            let d = t.mask.count_ones();
            if d == 0 {
                continue;
            }
            let w = match rule {
                SplitRule::MinTermWeight => (2.0f64).powi(1 - d as i32),
                _ => 1.0,
            };
            let mut m = t.mask;
            while m != 0 {
                let v = m.trailing_zeros() as usize;
                m &= m - 1;
                if v < score.len() && assignment[v].is_none() {
                    score[v] += w;
                }
            }
        }
    }
    // Highest score wins; ties go to the lowest index, so the rule is
    // deterministic and degrades to LowestFree on an empty system.
    let best = score
        .iter()
        .enumerate()
        .filter(|(v, _)| assignment[*v].is_none())
        .fold(
            (lowest, 0f64),
            |(bv, bs), (v, &sc)| {
                if sc > bs {
                    (v, sc)
                } else {
                    (bv, bs)
                }
            },
        );
    Some(best.0)
}

impl Default for SolveOptions {
    fn default() -> Self {
        Self {
            engine: SolverEngine::default(),
            max_solutions: 32,
            node_budget: 4096,
            split_rule: SplitRule::default(),
        }
    }
}

/// Statistics from a solve, so callers can report what the algebra
/// actually cost.
#[derive(Clone, Debug, Default)]
pub struct SolveStats {
    /// Algebraic reductions performed (F4 passes or Gröbner bases).
    pub reductions: usize,
    /// Branches closed by the reduction producing the constant `1` —
    /// the infeasibility certificates that replace exhaustive search.
    pub infeasible_branches: usize,
    /// Variables fixed by propagation rather than by splitting.
    pub propagations: usize,
    /// Splitting decisions made.
    pub splits: usize,
    /// True if the node budget ran out, so results may be incomplete.
    pub exhausted: bool,
    /// Highest Macaulay degree whose matrix was actually built.
    pub max_degree_built: u32,
    /// Reductions at which the next Macaulay matrix exceeded the size
    /// caps (the engine then split on what it had).
    pub oversize: usize,
}

/// Reduce `system`, returning polynomials in the same ideal — either a
/// Gröbner basis or the reduced Macaulay rows.  `None` means the F4
/// matrix would have been too large.
fn reduce_system(
    system: &[F2BoolPoly],
    n_vars: usize,
    engine: SolverEngine,
    stats: &mut SolveStats,
) -> Option<Vec<F2BoolPoly>> {
    stats.reductions += 1;
    match engine {
        SolverEngine::Buchberger => Some(groebner_basis_f2(system.to_vec(), n_vars)),
        SolverEngine::MatrixF4 { max_degree } => {
            let base = system
                .iter()
                .flat_map(|p| p.terms.iter())
                .map(|t| t.mask.count_ones())
                .max()
                .unwrap_or(0)
                .max(2);
            let mut best: Option<Vec<F2BoolPoly>> = None;
            for d in base..=max_degree.max(base) {
                match matrix_f4_f2(system, n_vars, d) {
                    Some(rows) => {
                        stats.max_degree_built = stats.max_degree_built.max(d);
                        let decisive = rows.iter().any(|p| is_constant_one(p))
                            || rows.iter().any(|p| forced_assignment(p).is_some());
                        best = Some(rows);
                        if decisive {
                            break;
                        }
                    }
                    None => {
                        // matrix too large: use what we have
                        stats.oversize += 1;
                        break;
                    }
                }
            }
            best
        }
    }
}

/// Is `p` the constant `1` — the infeasibility certificate?
fn is_constant_one(p: &F2BoolPoly) -> bool {
    p.terms.len() == 1 && p.terms[0].mask == 0
}

/// **Solve a Boolean system** by algebraic reduction plus splitting.
///
/// At each node the system is reduced ([`matrix_f4_f2`] by default,
/// Buchberger on request).  A reduction that yields the constant `1`
/// closes the branch — no search happens at all, which is the whole
/// point of doing algebra instead of enumeration.  Rows that have
/// collapsed to `v_i` or `v_i + 1` are propagated; when the algebra
/// stalls, the solver splits on the lowest free variable and recurses.
///
/// Every returned assignment is verified against the original
/// equations, so no spurious root can escape.
pub fn solve_boolean_system(
    equations: &[F2BoolPoly],
    n_vars: usize,
    opts: &SolveOptions,
) -> (Vec<u64>, SolveStats) {
    solve_boolean_system_filtered(equations, n_vars, opts, |_| false)
}

/// As [`solve_boolean_system`], but `accept` is called on each verified
/// solution as it is found; returning `true` stops the search
/// immediately (used to stop at the first root that lifts to a real
/// point decomposition).
pub fn solve_boolean_system_filtered(
    equations: &[F2BoolPoly],
    n_vars: usize,
    opts: &SolveOptions,
    mut accept: impl FnMut(u64) -> bool,
) -> (Vec<u64>, SolveStats) {
    let mut stats = SolveStats::default();
    let mut out = Vec::new();
    let mut stop = false;
    solve_rec(
        equations.to_vec(),
        equations,
        vec![None; n_vars],
        n_vars,
        opts,
        &mut stats,
        &mut out,
        &mut accept,
        &mut stop,
    );
    (out, stats)
}

#[allow(clippy::too_many_arguments)]
fn solve_rec(
    mut system: Vec<F2BoolPoly>,
    original: &[F2BoolPoly],
    mut assignment: Vec<Option<bool>>,
    n_vars: usize,
    opts: &SolveOptions,
    stats: &mut SolveStats,
    out: &mut Vec<u64>,
    accept: &mut impl FnMut(u64) -> bool,
    stop: &mut bool,
) {
    if *stop || out.len() >= opts.max_solutions {
        return;
    }
    if stats.reductions >= opts.node_budget {
        stats.exhausted = true;
        return;
    }

    // Reduce, propagate, repeat until the algebra stops learning.
    loop {
        system.retain(|p| !p.is_zero());
        if system.iter().any(is_constant_one) {
            stats.infeasible_branches += 1;
            return;
        }
        let reduced = match reduce_system(&system, n_vars, opts.engine, stats) {
            Some(r) => r,
            None => break, // no reduction available; split instead
        };
        if reduced.iter().any(is_constant_one) {
            stats.infeasible_branches += 1;
            return;
        }
        let forced: Vec<(u32, bool)> = reduced.iter().filter_map(forced_assignment).collect();
        if forced.is_empty() {
            break;
        }
        for (v, val) in forced {
            match assignment[v as usize] {
                Some(existing) if existing != val => {
                    stats.infeasible_branches += 1;
                    return;
                }
                Some(_) => {}
                None => {
                    stats.propagations += 1;
                    assignment[v as usize] = Some(val);
                    system = system.iter().map(|p| substitute(p, v, val)).collect();
                }
            }
        }
        if stats.reductions >= opts.node_budget {
            stats.exhausted = true;
            return;
        }
    }

    match choose_split(&system, &assignment, opts.split_rule) {
        None => {
            let mut pt = 0u64;
            for (i, a) in assignment.iter().enumerate() {
                if *a == Some(true) {
                    pt |= 1 << i;
                }
            }
            // Verify against the untouched system before accepting.
            if original.iter().all(|e| e.eval(pt) == 0) {
                out.push(pt);
                if accept(pt) {
                    *stop = true;
                }
            }
        }
        Some(free) => {
            stats.splits += 1;
            for value in [false, true] {
                let mut branch = assignment.clone();
                branch[free] = Some(value);
                let specialised: Vec<F2BoolPoly> = system
                    .iter()
                    .map(|p| substitute(p, free as u32, value))
                    .collect();
                solve_rec(
                    specialised,
                    original,
                    branch,
                    n_vars,
                    opts,
                    stats,
                    out,
                    accept,
                    stop,
                );
                if *stop || out.len() >= opts.max_solutions {
                    return;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_ecc::BinaryCurve;
    use crate::cryptanalysis::binary_semaev::binary_semaev_s3;
    use crate::cryptanalysis::koblitz_index_calculus::{find_irreducible, KoblitzCurve};

    fn fe(v: u64, n: u32) -> F2mElement {
        F2mElement::from_biguint(&num_bigint::BigUint::from(v), n)
    }

    /// Symbolic multiplication must agree with the field's own.
    #[test]
    fn symbolic_mul_matches_field_mul() {
        let n = 9;
        let irr = find_irreducible(n).unwrap();
        let st = FieldStructure::new(n, &irr);
        for (u, v) in [(1u64, 1u64), (5, 7), (300, 47), (511, 256), (0, 123)] {
            let a = fe(u, n);
            let b = fe(v, n);
            let sym = SymElement::constant(&a, n, 4).mul(&SymElement::constant(&b, n, 4), &st);
            assert_eq!(sym.eval(0, n), a.mul(&b, &irr), "{u} · {v}");
        }
    }

    /// Symbolic squaring is the linear shortcut, and equals `x · x`.
    #[test]
    fn symbolic_square_matches_mul_by_self() {
        let n = 9;
        let irr = find_irreducible(n).unwrap();
        let st = FieldStructure::new(n, &irr);
        let basis: Vec<F2mElement> = (0..n)
            .map(|k| F2mElement::from_bit_positions(&[k], n))
            .collect();
        let x = SymElement::from_subspace_vars(&basis, 0, n, n as usize);
        let sq = x.square(&st);
        let by_mul = x.mul(&x, &st);
        for point in [0u64, 1, 2, 5, 170, 341, 511] {
            assert_eq!(sq.eval(point, n), by_mul.eval(point, n), "point {point}");
            // …and against the field's squaring of the same value.
            assert_eq!(sq.eval(point, n), x.eval(point, n).square(&irr));
        }
    }

    /// The symbolic S₃ coordinates evaluate to the scalar S₃.
    #[test]
    fn symbolic_s3_matches_scalar_s3() {
        let n = 9;
        let irr = find_irreducible(n).unwrap();
        let st = FieldStructure::new(n, &irr);
        let b = F2mElement::one(n);
        let basis: Vec<F2mElement> = (0..n)
            .map(|k| F2mElement::from_bit_positions(&[k], n))
            .collect();

        let x1 = SymElement::from_subspace_vars(&basis, 0, n, 2 * n as usize);
        let x2 = SymElement::from_subspace_vars(&basis, n as usize, n, 2 * n as usize);
        let x3 = SymElement::constant(&fe(37, n), n, 2 * n as usize);
        let eqs = sym_semaev_s3(&x1, &x2, &x3, &b, &st);

        for point in [0u64, 1, 0b101, 0b1010_1010_1, 12345, 0o777] {
            let v1 = x1.eval(point, n);
            let v2 = x2.eval(point, n);
            let scalar = binary_semaev_s3(&v1, &v2, &fe(37, n), &b, &irr);
            let mut got = Vec::new();
            for (k, e) in eqs.iter().enumerate() {
                if e.eval(point) == 1 {
                    got.push(k as u32);
                }
            }
            assert_eq!(F2mElement::from_bit_positions(&got, n), scalar);
        }
    }

    /// A real decomposition is a root of the built system.
    #[test]
    fn known_decomposition_is_a_root_of_the_system() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = crate::cryptanalysis::koblitz_index_calculus::build_frobenius_factor_base(&kc, 0)
            .unwrap();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let basis = &fb.subspace_basis;

        // Pick two factor-base points and decompose their sum.
        let p1 = fb.points[0].clone();
        let p2 = fb.points[3].clone();
        let r = kc.add(&p1, &p2);
        let x_r = match &r {
            crate::binary_ecc::BinaryPoint::Affine { x, .. } => x.clone(),
            _ => panic!("P1 + P2 = O"),
        };
        let sys = build_decomposition_system(basis, &x_r, &kc.curve.b, 2, &st).unwrap();

        // Coordinates of x(P_i) in the subspace basis, as a Boolean point.
        let coords = |x: &F2mElement| -> u64 {
            for mask in 0..(1u64 << basis.len()) {
                let mut acc = F2mElement::zero(kc.n);
                for (t, bv) in basis.iter().enumerate() {
                    if (mask >> t) & 1 == 1 {
                        acc = acc.add(bv);
                    }
                }
                if &acc == x {
                    return mask;
                }
            }
            panic!("x not in the subspace");
        };
        let (x1, x2) = match (&p1, &p2) {
            (
                crate::binary_ecc::BinaryPoint::Affine { x: a, .. },
                crate::binary_ecc::BinaryPoint::Affine { x: b, .. },
            ) => (a.clone(), b.clone()),
            _ => panic!("factor base holds affine points"),
        };
        let point = coords(&x1) | (coords(&x2) << basis.len());
        assert!(
            sys.equations.iter().all(|e| e.eval(point) == 0),
            "the true decomposition must satisfy every equation"
        );
        // …and the layout helper reads the same x-coordinates back.
        assert_eq!(sys.summand_x(basis, point, 0, kc.n), x1);
        assert_eq!(sys.summand_x(basis, point, 1, kc.n), x2);
    }

    /// Degree: quadratic for `m = 2`, cubic once the chain introduces
    /// intermediate unknowns.
    #[test]
    fn system_degree_is_two_for_m2_and_three_when_chained() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = crate::cryptanalysis::koblitz_index_calculus::build_frobenius_factor_base(&kc, 0)
            .unwrap();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let sys =
            build_decomposition_system(&fb.subspace_basis, &fe(11, kc.n), &kc.curve.b, 2, &st)
                .unwrap();
        let deg = sys
            .equations
            .iter()
            .flat_map(|e| e.terms.iter())
            .map(|t| t.mask.count_ones())
            .max()
            .unwrap();
        assert_eq!(deg, 2);
        assert_eq!(sys.n_vars, 2 * fb.subspace_basis.len());
        assert_eq!(sys.equations.len(), kc.n as usize);

        // m = 3: one chaining unknown per link, and the term x₁x₂·e
        // pushes the degree to 3.
        let chained =
            build_decomposition_system(&fb.subspace_basis, &fe(11, kc.n), &kc.curve.b, 3, &st)
                .unwrap();
        let deg3 = chained
            .equations
            .iter()
            .flat_map(|e| e.terms.iter())
            .map(|t| t.mask.count_ones())
            .max()
            .unwrap();
        assert_eq!(deg3, 3);
        assert_eq!(chained.n_vars, 3 * fb.subspace_basis.len() + kc.n as usize);
        assert_eq!(chained.equations.len(), 2 * kc.n as usize);
    }

    /// The splitting solver finds exactly the roots, and closes
    /// infeasible branches by Gröbner rather than by enumeration.
    #[test]
    fn splitting_solver_agrees_with_exhaustive_evaluation() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = crate::cryptanalysis::koblitz_index_calculus::build_frobenius_factor_base(&kc, 0)
            .unwrap();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let basis = &fb.subspace_basis;

        for raw in [1u64, 9, 23, 64, 300] {
            let sys =
                build_decomposition_system(basis, &fe(raw, kc.n), &kc.curve.b, 2, &st).unwrap();
            let opts = SolveOptions {
                max_solutions: 64,
                ..SolveOptions::default()
            };
            let (sols, stats) = solve_boolean_system(&sys.equations, sys.n_vars, &opts);
            assert!(!stats.exhausted, "budget should suffice at this size");

            let mut brute = Vec::new();
            for pt in 0..(1u64 << sys.n_vars) {
                if sys.equations.iter().all(|e| e.eval(pt) == 0) {
                    brute.push(pt);
                }
            }
            let mut got = sols.clone();
            got.sort_unstable();
            brute.sort_unstable();
            assert_eq!(got, brute, "x_R = {raw}");
        }
    }

    /// An unsatisfiable system is rejected by the basis, not by search.
    #[test]
    fn infeasible_system_is_closed_by_the_basis() {
        let n_vars = 4;
        // v0 = 0, v0 = 1 — contradictory.
        let p = F2BoolPoly::from_monos(vec![F2BoolMono::var(0)], n_vars);
        let q = F2BoolPoly::from_monos(vec![F2BoolMono::var(0), F2BoolMono::one()], n_vars);
        let (sols, stats) = solve_boolean_system(&[p, q], n_vars, &SolveOptions::default());
        assert!(sols.is_empty());
        assert_eq!(stats.infeasible_branches, 1);
        // The algebra settled it; no splitting happened.
        assert_eq!(stats.splits, 0);
    }

    /// `FieldStructure` really is the multiplication table of the field.
    #[test]
    fn structure_constants_match_the_curve_field() {
        let curve = BinaryCurve::test_curve_f256();
        let st = FieldStructure::new(curve.m, &curve.irreducible);
        for i in 0..curve.m {
            for j in 0..curve.m {
                let zi = F2mElement::from_bit_positions(&[i], curve.m);
                let zj = F2mElement::from_bit_positions(&[j], curve.m);
                let want = zi.mul(&zj, &curve.irreducible);
                assert_eq!(st.reduced[i as usize][j as usize], want.raw_bits()[0]);
            }
        }
    }

    /// A splitting rule may change how fast the search closes, never
    /// what it finds: on random Boolean systems all three rules return
    /// the same solution set, and it is the true one.
    #[test]
    fn every_split_rule_finds_exactly_the_solutions() {
        use crate::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
        let mut seed = 0x9E3779B97F4A7C15u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for case in 0..12 {
            let n_vars = 6 + (case % 3);
            let n_eqs = 3 + (case % 4);
            let equations: Vec<F2BoolPoly> = (0..n_eqs)
                .map(|_| {
                    let terms: Vec<F2BoolMono> = (0..4)
                        .map(|_| {
                            let r = next();
                            // A degree-≤2 monomial over the first n_vars.
                            let a = (r % n_vars as u64) as u32;
                            let b = ((r >> 8) % n_vars as u64) as u32;
                            if (r >> 16) & 1 == 0 {
                                F2BoolMono::var(a)
                            } else {
                                F2BoolMono {
                                    mask: (1u64 << a) | (1u64 << b),
                                }
                            }
                        })
                        .collect();
                    F2BoolPoly::from_monos(terms, n_vars)
                })
                .collect();

            // Truth by enumeration.
            let mut truth: Vec<u64> = Vec::new();
            for pt in 0u64..(1u64 << n_vars) {
                if equations.iter().all(|e| e.eval(pt) == 0) {
                    truth.push(pt);
                }
            }

            for rule in [
                SplitRule::LowestFree,
                SplitRule::MostFrequent,
                SplitRule::MinTermWeight,
            ] {
                let opts = SolveOptions {
                    engine: SolverEngine::MatrixF4 { max_degree: 3 },
                    max_solutions: usize::MAX,
                    node_budget: 100_000,
                    split_rule: rule,
                };
                let (mut got, stats) = solve_boolean_system(&equations, n_vars, &opts);
                got.sort_unstable();
                got.dedup();
                assert!(
                    !stats.exhausted,
                    "case {case} rule {rule:?} ran out of budget"
                );
                assert_eq!(
                    got, truth,
                    "case {case} rule {rule:?} disagrees with enumeration"
                );
            }
        }
    }
}
