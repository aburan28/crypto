//! # Matrix-F5 over the Boolean ring: the F5 criterion for the Koblitz Macaulay matrices.
//!
//! [`super::koblitz_groebner::matrix_f4_f2`] builds, at a degree `d`, every
//! product `t·f_i` with `deg t ≤ d − deg f_i` and row-reduces the lot.  Most
//! of what makes that expensive is rows that reduce to **zero**: products
//! that are already linear combinations of the others.  Faugère's F5
//! criterion names a large family of them in advance — the ones a *trivial*
//! syzygy predicts — so they are never built, never packed and never
//! eliminated.  This module is that criterion, in the matrix form Bardet,
//! Faugère and Salvy use to analyse semi-regular systems, adapted to the
//! Boolean ring the Koblitz oracle works in.
//!
//! ## The criterion, and why the Boolean ring needs its own
//!
//! Order the generators `f_1, …, f_m` and write `d_i = deg f_i`.  Over a
//! polynomial ring the F5 criterion says: the row `t·f_i` is redundant at
//! degree `d` whenever `t` is the leading monomial of some element of the
//! ideal `⟨f_1, …, f_{i−1}⟩` in degree `d − d_i`, because with
//! `g = t + (lower terms)` in that ideal,
//!
//! ```text
//!     t·f_i = (t + g)·f_i + g·f_i,
//! ```
//!
//! where `(t+g)·f_i` is a sum of rows of `f_i` with smaller multipliers
//! and `g·f_i` lies in the span of the rows of `f_1, …, f_{i−1}`.  Induction
//! over the signature order `(i, t)` makes every pruned row a combination
//! of the rows kept.  Over `F_2` with the field equations `x² = x` there is
//! a second trivial syzygy, the Frobenius one: `f_i² = f_i`, i.e.
//! `(f_i + 1)·f_i = 0`, and with it `s·(f_i + 1)·f_i = 0` for every
//! monomial `s`.
//!
//! The textbook adaptation — "prune `t·f_i` when `t ∈ LM⟨f_1, …, f_i⟩`",
//! with `f_i` itself included to stand for the Frobenius syzygy — is **not
//! sound** in the Boolean ring, because multiplying by a monomial is not
//! order-compatible there: `x₁·(x₁x₂ + x₁ + x₂) = x₁`, a *smaller* leading
//! monomial.  The induction over signatures then runs the wrong way and a
//! pruned row can leave the row space (the test
//! `naive_f2_criterion_is_unsound_in_the_boolean_ring` keeps the
//! counterexample).  What is sound is to use the Frobenius syzygies as the
//! polynomials they are.  Let
//!
//! ```text
//!     W_i(d) = V_{i−1}(d − d_i)  +  span{ s·(f_i + 1) : deg s ≤ d − 2·d_i },
//!     V_j(e) = span{ s·f_k : k ≤ j, deg s ≤ e − d_k }.
//! ```
//!
//! Then for `t ∈ LM(W_i(d))` with witness `g = g' + σ` (`g' ∈ V_{i−1}`,
//! `σ` in the Frobenius span, `LM(g) = t`),
//!
//! ```text
//!     t·f_i = (t + g)·f_i + g'·f_i + σ·f_i = (t + g)·f_i + g'·f_i,
//! ```
//!
//! since `σ·f_i = 0` identically.  Every monomial of `t + g` is below `t`
//! in a degree-compatible order, so has degree `≤ d − d_i` and names a row
//! of `f_i` with smaller signature, and `g'·f_i` is a combination of rows of
//! `f_1, …, f_{i−1}` at degree `d`, because `(s·f_i)·f_k` expands into
//! multiples of `f_k` of degree `≤ d`.  That is the criterion implemented
//! here: **prune `(t, i)` iff `t ∈ LM(W_i(d))`.**  Pruning never changes the
//! row space, so a solver that reads consequences off the reduced rows gets
//! the same consequences from fewer rows.
//!
//! ## Cost accounting
//!
//! The leading-monomial sets come from echelonising the lower-degree
//! Macaulay matrices `V_j(e)` incrementally in generator order, one pass
//! per distinct `e = d − d_i`.  That work is real and is charged in the
//! same unit as the reduction it saves — 64-bit word XORs — through
//! [`F5Criterion::word_ops`], so a caller adding the two sees the
//! method's cost, not the phase's.  On a quadratic system at `d = 3` the
//! lower levels are the linear equations (level 1) and the equations
//! themselves plus linear equations times variables (level 2): a few
//! hundred rows over a few hundred columns against the degree-3 matrix's
//! thousands.
//!
//! ## What this is and is not
//!
//! This is the F5 *criterion* as a row filter over the Macaulay matrix,
//! with the same elimination kernels and the same output as the F4 step —
//! the ingredient of F5 that removes reductions to zero, without F5's
//! incremental basis construction or signature-compatible reduction, which
//! a fixed-degree matrix that only needs the row space does not require.
//! The zero reductions it removes are the trivial ones.  Semaev systems
//! have many non-trivial syzygies at low degree (their first fall degree
//! is where the attack lives), and those are not predicted here.
//! [`F5Report`] separates the two so a reader sees which share was
//! removable in principle.
//!
//! ## References
//!
//! - J.-C. Faugère, *A new efficient algorithm for computing Gröbner bases
//!   without reduction to zero (F5)*, ISSAC 2002.
//! - M. Bardet, J.-C. Faugère, B. Salvy, *On the complexity of Gröbner
//!   basis computation of semi-regular overdetermined algebraic
//!   equations*, ICPSS 2004 — the matrix-F5 formulation.
//! - M. Bardet, J.-C. Faugère, B. Salvy, B.-Y. Yang, *Asymptotic behaviour
//!   of the degree of regularity of semi-regular polynomial systems*,
//!   MEGA 2005 — the `F_2` variant with field equations.
//! - C. Eder, J.-C. Faugère, *A survey on signature-based algorithms for
//!   computing Gröbner bases*, J. Symbolic Comput. 80 (2017).

use crate::cryptanalysis::koblitz_groebner::{
    all_variable_mask, macaulay_columns, monomials_up_to_mask, occurring_vars,
};
use crate::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use std::collections::{HashMap, HashSet};

/// Degree of a Boolean polynomial (`0` for a constant or zero).
fn poly_degree(p: &F2BoolPoly) -> u32 {
    p.terms
        .iter()
        .map(|t| t.mask.count_ones())
        .max()
        .unwrap_or(0)
}

/// Product `p · s` as a set of monomial masks with even multiplicities
/// cancelled — the same row construction as the F4 builder, kept here so
/// the criterion's rows and the matrix's rows cannot disagree.
fn product_monos(p: &F2BoolPoly, s: u64) -> Vec<u64> {
    let mut all: Vec<u64> = p.terms.iter().map(|t| t.mask | s).collect();
    all.sort_unstable();
    let mut out = Vec::with_capacity(all.len());
    let mut i = 0;
    while i < all.len() {
        let mut j = i;
        while j < all.len() && all[j] == all[i] {
            j += 1;
        }
        if (j - i) % 2 == 1 {
            out.push(all[i]);
        }
        i = j;
    }
    out
}

/// A row echelon basis grown one row at a time, in generator order.
///
/// Columns are monomials in descending order, so a packed row's leading
/// monomial is its lowest set column.  Inserting a row reduces its leading
/// column against existing pivots until it either becomes a new pivot or
/// vanishes.  The pivot columns present after the rows of `f_1, …, f_j`
/// have been inserted are exactly `LM(V_j(e))`, and `created_at` records
/// which prefix each pivot first appeared in so that every prefix's
/// leading-monomial set can be read off one structure.
struct PrefixEchelon {
    columns: Vec<u64>,
    index: HashMap<u64, usize>,
    words: usize,
    /// `pivot[c]` is the index into `rows` of the row leading at column `c`.
    pivot: Vec<Option<usize>>,
    rows: Vec<Vec<u64>>,
    /// Prefix (generator index, 1-based) whose insertion created the pivot.
    created_at: Vec<usize>,
    word_ops: u64,
    inserted: u64,
    zero_reductions: u64,
}

fn leading_column(row: &[u64]) -> Option<usize> {
    row.iter()
        .enumerate()
        .find(|(_, &w)| w != 0)
        .map(|(i, &w)| i * 64 + w.trailing_zeros() as usize)
}

impl PrefixEchelon {
    fn new(columns: Vec<u64>) -> Self {
        let index = columns.iter().enumerate().map(|(i, &m)| (m, i)).collect();
        let words = columns.len().div_ceil(64).max(1);
        Self {
            pivot: vec![None; columns.len()],
            columns,
            index,
            words,
            rows: Vec::new(),
            created_at: Vec::new(),
            word_ops: 0,
            inserted: 0,
            zero_reductions: 0,
        }
    }

    fn pack(&self, monos: &[u64]) -> Option<Vec<u64>> {
        let mut row = vec![0u64; self.words];
        for m in monos {
            let &c = self.index.get(m)?;
            row[c / 64] |= 1u64 << (c % 64);
        }
        Some(row)
    }

    /// Reduce `row`'s leading column against `self` (and `extra`, a
    /// temporary echelon over the same columns) until it is a fresh
    /// leading column or zero.  Returns that column.
    fn reduce_leading(&mut self, row: &mut [u64], extra: Option<&TempEchelon>) -> Option<usize> {
        loop {
            let c = leading_column(row)?;
            let w = c / 64;
            if let Some(r) = self.pivot[c] {
                let pivot = &self.rows[r];
                for (dst, &src) in row[w..].iter_mut().zip(&pivot[w..]) {
                    *dst ^= src;
                }
                self.word_ops += (self.words - w) as u64;
                continue;
            }
            if let Some(extra) = extra {
                if let Some(r) = extra.pivot.get(&c) {
                    let pivot = &extra.rows[*r];
                    for (dst, &src) in row[w..].iter_mut().zip(&pivot[w..]) {
                        *dst ^= src;
                    }
                    self.word_ops += (self.words - w) as u64;
                    continue;
                }
            }
            return Some(c);
        }
    }

    /// Insert one row of generator `prefix` (1-based).
    fn insert(&mut self, mut row: Vec<u64>, prefix: usize) {
        self.inserted += 1;
        match self.reduce_leading(&mut row, None) {
            Some(c) => {
                self.pivot[c] = Some(self.rows.len());
                self.rows.push(row);
                self.created_at.push(prefix);
            }
            None => self.zero_reductions += 1,
        }
    }

    /// Leading monomials of `V_j(e)` for `j = prefix`.
    fn leading_monomials(&self, prefix: usize) -> impl Iterator<Item = u64> + '_ {
        self.rows
            .iter()
            .zip(&self.created_at)
            .filter(move |(_, &created)| created <= prefix)
            .filter_map(|(row, _)| leading_column(row))
            .map(|c| self.columns[c])
    }
}

/// Pivots added on top of a [`PrefixEchelon`] for one generator's
/// Frobenius span, then discarded.
#[derive(Default)]
struct TempEchelon {
    pivot: HashMap<usize, usize>,
    rows: Vec<Vec<u64>>,
}

/// The F5 criterion evaluated for one system at one degree: which
/// multipliers of which generator are pruned, and what it cost to know.
#[derive(Clone, Debug)]
pub struct F5Criterion {
    /// `pruned[i]` is the set of multiplier monomials `t` with
    /// `t ∈ LM(W_i(d))`; `t·f_i` is redundant at this degree.
    pruned: Vec<HashSet<u64>>,
    /// Word XORs spent echelonising the lower-degree matrices.
    word_ops: u64,
    /// Rows inserted into the lower-degree echelons.
    lower_rows: u64,
    /// Of those, rows that reduced to zero.
    lower_zero_reductions: u64,
    /// Pruned multipliers predicted by the Koszul part `V_{i−1}(d − d_i)`.
    koszul_pruned: u64,
    /// Pruned multipliers predicted only once the Frobenius span is added.
    frobenius_pruned: u64,
}

impl F5Criterion {
    /// Evaluate the criterion for `polys` (in the given order) at `degree`,
    /// with multipliers drawn from `multiplier_mask`.
    ///
    /// The lower-level matrices are built over the variables occurring in
    /// the system together with the multiplier mask, which is the subring
    /// every row of the degree-`degree` matrix lives in.
    pub fn new(polys: &[F2BoolPoly], n_vars: usize, degree: u32, multiplier_mask: u64) -> Self {
        let m = polys.len();
        let mut out = Self {
            pruned: vec![HashSet::new(); m],
            word_ops: 0,
            lower_rows: 0,
            lower_zero_reductions: 0,
            koszul_pruned: 0,
            frobenius_pruned: 0,
        };
        if polys.is_empty() {
            return out;
        }
        let degrees: Vec<u32> = polys.iter().map(poly_degree).collect();
        // A constant generator (`1`, or zero) is outside the criterion's
        // hypotheses; the solver never passes one, so prune nothing.
        if degrees.contains(&0) {
            return out;
        }
        let column_mask = occurring_vars(polys) | (multiplier_mask & all_variable_mask(n_vars));

        // Distinct lower levels e = d − d_i, each echelonised once.
        let mut levels: Vec<u32> = degrees
            .iter()
            .filter(|&&di| di <= degree)
            .map(|&di| degree - di)
            .collect();
        levels.sort_unstable();
        levels.dedup();

        for e in levels {
            let columns = match macaulay_columns(&[monomials_up_to_mask(column_mask, e)]) {
                Some(c) => c,
                None => continue,
            };
            let mut echelon = PrefixEchelon::new(columns);
            for (i, (p, &di)) in polys.iter().zip(&degrees).enumerate() {
                let prefix = i + 1;
                // The criterion for f_i at this level is read before f_i's
                // own rows enter: W_i(d) = V_{i−1}(e) + Frobenius span of f_i.
                if di <= degree && degree - di == e {
                    let mut pruned: HashSet<u64> = echelon.leading_monomials(i).collect();
                    out.koszul_pruned += pruned.len() as u64;
                    if degree >= 2 * di {
                        let gap = degree - 2 * di;
                        let mut temp = TempEchelon::default();
                        let f_plus_one = p.add(&F2BoolPoly::one(p.n_vars));
                        for s in monomials_up_to_mask(multiplier_mask, gap) {
                            let monos = product_monos(&f_plus_one, s);
                            if monos.is_empty() {
                                continue;
                            }
                            let Some(mut row) = echelon.pack(&monos) else {
                                continue;
                            };
                            if let Some(c) = echelon.reduce_leading(&mut row, Some(&temp)) {
                                temp.pivot.insert(c, temp.rows.len());
                                temp.rows.push(row);
                            }
                        }
                        for &c in temp.pivot.keys() {
                            if pruned.insert(echelon.columns[c]) {
                                out.frobenius_pruned += 1;
                            }
                        }
                    }
                    out.pruned[i] = pruned;
                }
                if di > e {
                    continue;
                }
                for s in monomials_up_to_mask(multiplier_mask, e - di) {
                    let monos = product_monos(p, s);
                    if monos.is_empty() {
                        continue;
                    }
                    if let Some(row) = echelon.pack(&monos) {
                        echelon.insert(row, prefix);
                    }
                }
            }
            out.word_ops += echelon.word_ops;
            out.lower_rows += echelon.inserted;
            out.lower_zero_reductions += echelon.zero_reductions;
        }
        out
    }

    /// Is the row `t·f_i` (generator index `i`, 0-based) pruned?
    pub fn prunes(&self, i: usize, t: u64) -> bool {
        self.pruned.get(i).is_some_and(|set| set.contains(&t))
    }

    /// Total multipliers pruned across all generators.
    pub fn pruned_count(&self) -> u64 {
        self.pruned.iter().map(|s| s.len() as u64).sum()
    }

    /// Word XORs spent evaluating the criterion.
    pub fn word_ops(&self) -> u64 {
        self.word_ops
    }

    /// Rows inserted into the lower-degree echelons and how many of them
    /// reduced to zero.
    pub fn lower_level_rows(&self) -> (u64, u64) {
        (self.lower_rows, self.lower_zero_reductions)
    }

    /// How many pruned multipliers each part of the criterion predicted:
    /// `(koszul, frobenius)`.
    pub fn pruned_by_part(&self) -> (u64, u64) {
        (self.koszul_pruned, self.frobenius_pruned)
    }
}

/// What one matrix-F5 step did, beside its rows.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct F5Report {
    /// Degree the matrix was built at.
    pub degree: u32,
    /// Rows the plain F4 step would have built.
    pub rows_f4: u64,
    /// Rows the criterion removed.
    pub rows_pruned: u64,
    /// Rows actually built and reduced.
    pub rows_built: u64,
    /// Columns of the matrix.
    pub cols: u64,
    /// Rank of the reduced matrix (equal to the F4 rank: same row space).
    pub rank: u64,
    /// Rows built that still reduced to zero — the non-trivial syzygies
    /// the criterion does not see.
    pub zero_reductions: u64,
    /// Word XORs in the degree-`degree` elimination.
    pub reduce_word_ops: u64,
    /// Word XORs spent evaluating the criterion.
    pub criterion_word_ops: u64,
    /// Rows of the lower-degree echelons the criterion built.
    pub criterion_rows: u64,
}

impl F5Report {
    /// Total word XORs charged to the step: elimination plus criterion.
    pub fn word_ops(&self) -> u64 {
        self.reduce_word_ops + self.criterion_word_ops
    }
}

/// **Matrix-F5 step over `F_2`**: the Macaulay matrix of `polys` at
/// `degree` with the F5 criterion applied, reduced to row echelon form and
/// returned as polynomials — the same row space as
/// [`super::koblitz_groebner::matrix_f4_f2`] from fewer rows.
///
/// Returns `None` if the matrix would exceed the size limits the F4 step
/// uses.
pub fn matrix_f5_f2(
    polys: &[F2BoolPoly],
    n_vars: usize,
    degree: u32,
) -> Option<(Vec<F2BoolPoly>, F5Report)> {
    matrix_f5_f2_timed(polys, n_vars, degree).map(|(rows, report, _)| (rows, report))
}

/// Wall time of each phase of one [`matrix_f5_f2`] step, in nanoseconds.
/// Kept out of [`F5Report`], which is compared and stored as a record of
/// what a step did.
#[derive(Clone, Copy, Debug, Default, serde::Serialize)]
pub struct F5Timings {
    /// Evaluating the criterion: the lower-degree echelons.
    pub criterion_ns: u64,
    /// Listing the surviving rows, the column map and the packing.
    pub build_ns: u64,
    /// The elimination.
    pub reduce_ns: u64,
    /// Turning the pivot rows back into polynomials.
    pub unpack_ns: u64,
    /// Whether the full-column direct packed-row builder was used.
    pub direct_pack_used: bool,
    /// Whether the scalar preallocated direct-write unpack path was used.
    pub direct_unpack_used: bool,
    /// Whether AVX2 was selected for the shared GF(2) table builder.
    pub avx2_table_build_used: bool,
}

/// Row form requested from a matrix-F5 step. Both forms span the same
/// space; only `Reduced` has a unique list of returned polynomials.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum F5OutputForm {
    Reduced,
    Echelon,
    /// Echelon only for degree-4 systems with at least 20 variables,
    /// where the saved reduction outweighed unpacking in the first study.
    SelectiveEchelon,
}

/// Expand a packed pivot row into its canonical, descending monomial list.
fn unpack_row_scalar(row: &[u64], cols: &[u64], n_vars: usize) -> F2BoolPoly {
    let terms = row.iter().map(|w| w.count_ones() as usize).sum();
    let mut monos = Vec::with_capacity(terms);
    for (w, &word) in row.iter().enumerate() {
        let mut bits = word;
        while bits != 0 {
            let c = w * 64 + bits.trailing_zeros() as usize;
            bits &= bits - 1;
            monos.push(F2BoolMono::from_mask(cols[c]));
        }
    }
    F2BoolPoly {
        terms: monos,
        n_vars,
    }
}

/// Expand a packed row into the exact preallocated number of terms.
fn unpack_row_scalar_direct(row: &[u64], cols: &[u64], n_vars: usize) -> F2BoolPoly {
    let terms = row.iter().map(|w| w.count_ones() as usize).sum();
    let mut monos = Vec::<F2BoolMono>::with_capacity(terms);
    let out = monos.as_mut_ptr();
    let mut written = 0;
    for (w, &word) in row.iter().enumerate() {
        let mut bits = word;
        while bits != 0 {
            let c = w * 64 + bits.trailing_zeros() as usize;
            bits &= bits - 1;
            let mono = F2BoolMono::from_mask(cols[c]);
            // SAFETY: `written` increases once per set bit and `terms` is
            // exactly the number of set bits in the whole row.
            unsafe { out.add(written).write(mono) };
            written += 1;
        }
    }
    debug_assert_eq!(written, terms);
    // SAFETY: all `terms` slots have been written above.
    unsafe { monos.set_len(terms) };
    F2BoolPoly {
        terms: monos,
        n_vars,
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn unpack_row_avx512(row: &[u64], cols: &[u64], n_vars: usize) -> F2BoolPoly {
    use std::arch::x86_64::{_mm512_loadu_si512, _mm512_mask_compressstoreu_epi64};

    let terms: usize = row.iter().map(|w| w.count_ones() as usize).sum();
    let mut monos: Vec<F2BoolMono> = Vec::with_capacity(terms);
    let out = monos.as_mut_ptr().cast::<u64>();
    let mut written = 0;
    for (w, &word) in row.iter().enumerate() {
        for byte in 0..8 {
            let bits = (word >> (byte * 8)) as u8;
            if bits == 0 {
                continue;
            }
            let base = w * 64 + byte * 8;
            if base + 8 <= cols.len() {
                // SAFETY: the eight columns exist; `out` has capacity for
                // every set bit and F2BoolMono is transparent over u64.
                let values = unsafe { _mm512_loadu_si512(cols.as_ptr().add(base).cast()) };
                unsafe { _mm512_mask_compressstoreu_epi64(out.add(written).cast(), bits, values) };
                written += bits.count_ones() as usize;
            } else {
                let mut tail = bits;
                while tail != 0 {
                    let c = base + tail.trailing_zeros() as usize;
                    debug_assert!(c < cols.len());
                    unsafe { out.add(written).write(cols[c]) };
                    written += 1;
                    tail &= tail - 1;
                }
            }
        }
    }
    debug_assert_eq!(written, terms);
    // SAFETY: each of the `terms` slots was initialized exactly once.
    unsafe { monos.set_len(terms) };
    F2BoolPoly {
        terms: monos,
        n_vars,
    }
}

fn unpack_row(row: &[u64], cols: &[u64], n_vars: usize, avx512: bool, direct: bool) -> F2BoolPoly {
    #[cfg(target_arch = "x86_64")]
    if avx512 {
        // SAFETY: the caller selects this path only after a runtime check.
        return unsafe { unpack_row_avx512(row, cols, n_vars) };
    }
    let _ = avx512;
    if direct {
        return unpack_row_scalar_direct(row, cols, n_vars);
    }
    unpack_row_scalar(row, cols, n_vars)
}

/// [`matrix_f5_f2`] with the wall time of each phase.
pub fn matrix_f5_f2_timed(
    polys: &[F2BoolPoly],
    n_vars: usize,
    degree: u32,
) -> Option<(Vec<F2BoolPoly>, F5Report, F5Timings)> {
    matrix_f5_f2_with_form_timed(polys, n_vars, degree, F5OutputForm::Reduced)
}

/// A timed matrix-F5 step with an explicit output form. `Echelon` skips
/// above-pivot reduction; callers that need a canonical basis can request
/// `Reduced` or reduce the returned rows themselves.
pub fn matrix_f5_f2_with_form_timed(
    polys: &[F2BoolPoly],
    n_vars: usize,
    degree: u32,
    form: F5OutputForm,
) -> Option<(Vec<F2BoolPoly>, F5Report, F5Timings)> {
    use std::time::Instant;
    let mut timings = F5Timings::default();
    use crate::cryptanalysis::koblitz_groebner::{
        echelon_f2_counted, f5_rows_monos_with_f4_count, f5_rows_monos_with_mask,
        f5_rows_packed_full_columns, macaulay_row_count, pack_rows, rref_f2_counted,
    };
    let mut report = F5Report {
        degree,
        ..Default::default()
    };
    if polys.is_empty() {
        return Some((Vec::new(), report, timings));
    }
    let t = Instant::now();
    let mask = all_variable_mask(n_vars);
    let criterion = F5Criterion::new(polys, n_vars, degree, mask);
    timings.criterion_ns = t.elapsed().as_nanos() as u64;
    let t = Instant::now();
    report.criterion_word_ops = criterion.word_ops();
    report.criterion_rows = criterion.lower_level_rows().0;
    static FUSED_BUILD: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    let fused_build =
        *FUSED_BUILD.get_or_init(|| std::env::var("KIC_F5_FUSED_BUILD").as_deref() == Ok("1"));
    static DIRECT_PACK: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    let direct_pack =
        *DIRECT_PACK.get_or_init(|| std::env::var("KIC_F5_DIRECT_PACK").as_deref() == Ok("1"));
    let (cols, mut matrix) = if direct_pack {
        if let Some((rows_f4, cols, matrix)) =
            f5_rows_packed_full_columns(polys, n_vars, degree, &criterion)
        {
            timings.direct_pack_used = true;
            report.rows_f4 = rows_f4 as u64;
            (cols, matrix)
        } else {
            let (rows_f4, rows) =
                f5_rows_monos_with_f4_count(polys, n_vars, degree, mask, &criterion)?;
            report.rows_f4 = rows_f4 as u64;
            let cols = macaulay_columns(&rows)?;
            let matrix = pack_rows(&rows, &cols);
            (cols, matrix)
        }
    } else {
        let rows_monos = if fused_build {
            let (rows_f4, rows) =
                f5_rows_monos_with_f4_count(polys, n_vars, degree, mask, &criterion)?;
            report.rows_f4 = rows_f4 as u64;
            rows
        } else {
            report.rows_f4 = macaulay_row_count(polys, n_vars, degree)? as u64;
            f5_rows_monos_with_mask(polys, n_vars, degree, mask, &criterion)?
        };
        if rows_monos.is_empty() {
            return Some((Vec::new(), report, timings));
        }
        let cols = macaulay_columns(&rows_monos)?;
        let matrix = pack_rows(&rows_monos, &cols);
        (cols, matrix)
    };
    report.rows_built = matrix.len() as u64;
    report.rows_pruned = report.rows_f4 - report.rows_built;
    if matrix.is_empty() {
        return Some((Vec::new(), report, timings));
    }
    timings.build_ns = t.elapsed().as_nanos() as u64;
    let t = Instant::now();
    let mut word_ops = 0u64;
    let echelon = matches!(form, F5OutputForm::Echelon)
        || (matches!(form, F5OutputForm::SelectiveEchelon) && degree == 4 && n_vars >= 20);
    timings.avx2_table_build_used = echelon
        && matrix.len() >= 128
        && cols.len() >= 256
        && crate::cryptanalysis::gf2_elim::avx2_table_build_enabled();
    let rank = if !echelon {
        rref_f2_counted(&mut matrix, cols.len(), &mut word_ops)
    } else if matrix.len() >= 128 && cols.len() >= 256 {
        crate::cryptanalysis::gf2_elim::echelon_counted(&mut matrix, cols.len(), &mut word_ops)
    } else {
        echelon_f2_counted(&mut matrix, cols.len(), &mut word_ops)
    };
    timings.reduce_ns = t.elapsed().as_nanos() as u64;
    let t = Instant::now();
    report.cols = cols.len() as u64;
    report.rank = rank as u64;
    report.zero_reductions = report.rows_built - rank as u64;
    report.reduce_word_ops = word_ops;
    let n_vars_out = polys[0].n_vars;
    // One polynomial per pivot row. Dense rows make this a material part
    // of the complete call even when back-substitution is skipped.
    use rayon::prelude::*;
    static AVX512_UNPACK: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    let avx512_unpack = *AVX512_UNPACK.get_or_init(|| {
        #[cfg(target_arch = "x86_64")]
        {
            std::env::var("KIC_F5_AVX512_UNPACK").as_deref() == Ok("1")
                && std::arch::is_x86_feature_detected!("avx512f")
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            false
        }
    });
    static DIRECT_UNPACK: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    let direct_unpack =
        *DIRECT_UNPACK.get_or_init(|| std::env::var("KIC_F5_UNPACK_DIRECT").as_deref() == Ok("1"));
    timings.direct_unpack_used = direct_unpack && !avx512_unpack;
    let out = matrix[..rank]
        .par_iter()
        .map(|row| {
            let p = unpack_row(row, &cols, n_vars_out, avx512_unpack, direct_unpack);
            debug_assert!(p.is_canonical());
            p
        })
        .filter(|p| !p.is_zero())
        .collect();
    timings.unpack_ns = t.elapsed().as_nanos() as u64;
    Some((out, report, timings))
}

/// Fingerprint the unique RREF of a list of Boolean rows. This is a
/// diagnostic for comparing reduced and echelon output on frozen systems;
/// its work is outside the timed F5 step.
pub fn canonical_row_space_fingerprint(rows: &[F2BoolPoly]) -> Option<u64> {
    use crate::cryptanalysis::koblitz_groebner::{pack_rows, rref_f2_counted};
    use std::hash::{Hash, Hasher};
    let rows_monos: Vec<Vec<u64>> = rows
        .iter()
        .map(|p| p.terms.iter().map(|t| t.mask).collect())
        .collect();
    let cols = macaulay_columns(&rows_monos)?;
    let mut matrix = pack_rows(&rows_monos, &cols);
    let mut ops = 0;
    let rank = rref_f2_counted(&mut matrix, cols.len(), &mut ops);
    let mut hash = std::collections::hash_map::DefaultHasher::new();
    for row in &matrix[..rank] {
        for (w, &word) in row.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                let c = w * 64 + bits.trailing_zeros() as usize;
                bits &= bits - 1;
                cols[c].hash(&mut hash);
            }
        }
        u64::MAX.hash(&mut hash);
    }
    Some(hash.finish())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::koblitz_groebner::{
        f5_rows_monos_with_f4_count, f5_rows_packed_full_columns, matrix_f4_f2, pack_rows,
    };

    #[test]
    fn direct_scalar_unpack_matches_push_path() {
        for cols_len in [1usize, 7, 64, 65, 127, 129, 4097] {
            let cols: Vec<u64> = (0..cols_len)
                .map(|i| (i as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15))
                .collect();
            for salt in [0u64, 1, 3, 7, u64::MAX] {
                let row: Vec<u64> = (0..cols_len.div_ceil(64))
                    .map(|w| {
                        let bits = (w as u64).wrapping_mul(0xd6e8_feb8_6659_fd93) ^ salt;
                        let valid = (cols_len - w * 64).min(64);
                        bits & (u64::MAX >> (64 - valid))
                    })
                    .collect();
                assert_eq!(
                    unpack_row_scalar_direct(&row, &cols, 24),
                    unpack_row_scalar(&row, &cols, 24),
                );
            }
        }
    }

    #[test]
    fn direct_packed_f5_build_matches_sorted_builder_and_falls_back_on_sparse_columns() {
        let n_vars = 4;
        let degree = 3;
        let mask = all_variable_mask(n_vars);
        let dense = F2BoolPoly::from_monos(
            monomials_up_to_mask(mask, 2)
                .into_iter()
                .map(F2BoolMono::from_mask)
                .collect(),
            n_vars,
        );
        let criterion = F5Criterion::new(std::slice::from_ref(&dense), n_vars, degree, mask);
        let (direct_count, direct_cols, direct_rows) =
            f5_rows_packed_full_columns(std::slice::from_ref(&dense), n_vars, degree, &criterion)
                .expect("dense polynomial spans the complete column universe");
        let (normal_count, rows_monos) = f5_rows_monos_with_f4_count(
            std::slice::from_ref(&dense),
            n_vars,
            degree,
            mask,
            &criterion,
        )
        .unwrap();
        let normal_cols = macaulay_columns(&rows_monos).unwrap();
        let normal_rows = pack_rows(&rows_monos, &normal_cols);
        assert_eq!(direct_count, normal_count);
        assert_eq!(direct_cols, normal_cols);
        assert_eq!(direct_rows, normal_rows);

        let sparse = poly(n_vars, &[&[0], &[1]]);
        let sparse_criterion =
            F5Criterion::new(std::slice::from_ref(&sparse), n_vars, degree, mask);
        assert!(f5_rows_packed_full_columns(
            std::slice::from_ref(&sparse),
            n_vars,
            degree,
            &sparse_criterion
        )
        .is_none());
    }

    #[test]
    fn avx512_unpack_matches_scalar_for_sparse_dense_and_partial_words() {
        #[cfg(target_arch = "x86_64")]
        if std::arch::is_x86_feature_detected!("avx512f") {
            for cols_len in [1usize, 7, 8, 9, 63, 64, 65, 127, 128, 129, 4097] {
                let cols: Vec<u64> = (0..cols_len)
                    .map(|i| (i as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15))
                    .collect();
                for salt in [0u64, 1, 3, 7, u64::MAX] {
                    let row: Vec<u64> = (0..cols_len.div_ceil(64))
                        .map(|w| {
                            let bits = (w as u64).wrapping_mul(0xd6e8_feb8_6659_fd93) ^ salt;
                            let valid = (cols_len - w * 64).min(64);
                            bits & (u64::MAX >> (64 - valid))
                        })
                        .collect();
                    let scalar = unpack_row_scalar(&row, &cols, 24);
                    let simd = unsafe { unpack_row_avx512(&row, &cols, 24) };
                    assert_eq!(simd, scalar, "cols_len={cols_len}, salt={salt}");
                }
            }
        }
    }

    fn poly(n_vars: usize, monos: &[&[u32]]) -> F2BoolPoly {
        F2BoolPoly::from_monos(
            monos
                .iter()
                .map(|vars| {
                    F2BoolMono::from_mask(vars.iter().fold(0u64, |acc, &v| acc | (1u64 << v)))
                })
                .collect(),
            n_vars,
        )
    }

    /// Deterministic random Boolean polynomial of degree ≤ `deg`.
    fn random_poly(n_vars: usize, deg: u32, terms: usize, seed: &mut u64) -> F2BoolPoly {
        let mut next = || {
            *seed ^= *seed << 13;
            *seed ^= *seed >> 7;
            *seed ^= *seed << 17;
            *seed
        };
        let mut monos = Vec::new();
        for _ in 0..terms {
            let d = (next() % (deg as u64 + 1)) as u32;
            let mut mask = 0u64;
            for _ in 0..d {
                mask |= 1u64 << (next() % n_vars as u64);
            }
            monos.push(F2BoolMono::from_mask(mask));
        }
        F2BoolPoly::from_monos(monos, n_vars)
    }

    /// Reduced row echelon basis of a polynomial list (no multiplication),
    /// as a canonical set.
    fn row_space_canonical(polys: &[F2BoolPoly], n_vars: usize) -> Vec<F2BoolPoly> {
        use crate::cryptanalysis::koblitz_groebner::{pack_rows, rref_f2_counted};
        let rows_monos: Vec<Vec<u64>> = polys
            .iter()
            .filter(|p| !p.is_zero())
            .map(|p| p.terms.iter().map(|t| t.mask).collect())
            .collect();
        if rows_monos.is_empty() {
            return Vec::new();
        }
        let cols = macaulay_columns(&rows_monos).unwrap();
        let mut matrix = pack_rows(&rows_monos, &cols);
        let mut ops = 0;
        let rank = rref_f2_counted(&mut matrix, cols.len(), &mut ops);
        let mut rows: Vec<F2BoolPoly> = matrix
            .iter()
            .take(rank)
            .map(|row| {
                let monos: Vec<F2BoolMono> = (0..cols.len())
                    .filter(|c| row[c / 64] & (1u64 << (c % 64)) != 0)
                    .map(|c| F2BoolMono::from_mask(cols[c]))
                    .collect();
                F2BoolPoly::from_monos(monos, n_vars)
            })
            .collect();
        rows.sort_by(|a, b| format!("{a:?}").cmp(&format!("{b:?}")));
        rows
    }

    #[test]
    fn f5_rows_span_the_f4_row_space_on_random_systems() {
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        for trial in 0..40 {
            let n_vars = 5 + (trial % 4);
            let m = 3 + (trial % 5);
            let polys: Vec<F2BoolPoly> = (0..m)
                .map(|k| {
                    let deg = if k % 3 == 0 { 1 } else { 2 };
                    random_poly(n_vars, deg, 4 + k, &mut seed)
                })
                .filter(|p| poly_degree(p) >= 1)
                .collect();
            if polys.is_empty() {
                continue;
            }
            for degree in 2..=4u32 {
                let f4 = matrix_f4_f2(&polys, n_vars, degree).unwrap();
                let (f5, report) = matrix_f5_f2(&polys, n_vars, degree).unwrap();
                assert_eq!(
                    row_space_canonical(&f4, n_vars),
                    row_space_canonical(&f5, n_vars),
                    "trial {trial} degree {degree}: row spaces differ"
                );
                assert_eq!(f4.len() as u64, report.rank, "rank must agree");
                assert_eq!(report.rows_built + report.rows_pruned, report.rows_f4);
            }
        }
    }

    #[test]
    fn fused_row_build_matches_two_pass_rows_and_full_count() {
        use crate::cryptanalysis::koblitz_groebner::{
            f5_rows_monos_with_f4_count, f5_rows_monos_with_mask, macaulay_row_count,
        };
        let mut seed = 0x7a11_f5c0_1d5e_2028u64;
        for trial in 0..32 {
            let n_vars = 5 + trial % 4;
            let mut polys: Vec<F2BoolPoly> = (0..(3 + trial % 5))
                .map(|k| random_poly(n_vars, 2, 4 + k, &mut seed))
                .filter(|p| poly_degree(p) >= 1)
                .collect();
            // Distinct terms can map to the same mask after multiplication
            // by x_1, so the F4 count must exclude a cancelled zero row.
            polys.push(poly(n_vars, &[&[0], &[0, 1]]));
            polys.push(poly(n_vars, &[&[0], &[0, 1], &[2]]));
            let mask = all_variable_mask(n_vars);
            for degree in 2..=4 {
                let criterion = F5Criterion::new(&polys, n_vars, degree, mask);
                let old_count = macaulay_row_count(&polys, n_vars, degree).unwrap();
                let old_rows =
                    f5_rows_monos_with_mask(&polys, n_vars, degree, mask, &criterion).unwrap();
                let (new_count, new_rows) =
                    f5_rows_monos_with_f4_count(&polys, n_vars, degree, mask, &criterion).unwrap();
                assert_eq!(new_count, old_count, "trial {trial} degree {degree}");
                assert_eq!(new_rows, old_rows, "trial {trial} degree {degree}");
            }
        }
    }

    #[test]
    fn echelon_option_preserves_the_f5_row_space() {
        let mut seed = 0x8102_34ab_cdef_9876u64;
        for trial in 0..16 {
            let n_vars = 5 + trial % 4;
            let polys: Vec<F2BoolPoly> = (0..(4 + trial % 3))
                .map(|k| random_poly(n_vars, 2, 5 + k, &mut seed))
                .filter(|p| poly_degree(p) >= 1)
                .collect();
            for degree in [3, 4] {
                let (reduced, r_report, _) =
                    matrix_f5_f2_with_form_timed(&polys, n_vars, degree, F5OutputForm::Reduced)
                        .unwrap();
                let (echelon, e_report, _) =
                    matrix_f5_f2_with_form_timed(&polys, n_vars, degree, F5OutputForm::Echelon)
                        .unwrap();
                assert_eq!(r_report.rank, e_report.rank);
                assert_eq!(r_report.rows_built, e_report.rows_built);
                assert_eq!(r_report.rows_pruned, e_report.rows_pruned);
                assert_eq!(
                    row_space_canonical(&reduced, n_vars),
                    row_space_canonical(&echelon, n_vars)
                );
                assert_eq!(
                    canonical_row_space_fingerprint(&reduced),
                    canonical_row_space_fingerprint(&echelon)
                );
            }
        }
    }

    #[test]
    fn frobenius_criterion_prunes_the_leading_monomial_multiple() {
        // One quadratic generator at degree 4: the only trivial syzygy is
        // Frobenius, so exactly one multiplier — LM(f) — is pruned, and the
        // pruned row is still in the span.
        let n_vars = 4;
        let f = poly(n_vars, &[&[0, 1], &[2], &[]]);
        let (_, report) = matrix_f5_f2(std::slice::from_ref(&f), n_vars, 4).unwrap();
        assert_eq!(report.rows_pruned, 1);
        let c = F5Criterion::new(
            std::slice::from_ref(&f),
            n_vars,
            4,
            all_variable_mask(n_vars),
        );
        assert!(c.prunes(0, f.lt().unwrap().mask));
        assert_eq!(c.pruned_by_part(), (0, 1));
    }

    #[test]
    fn koszul_criterion_prunes_leading_monomials_of_earlier_generators() {
        // Two linear generators, degree 2: t·f_2 is pruned for t = LM(f_1)
        // (Koszul) and t = LM(f_2) (Frobenius, since 2 ≥ 2·1).
        let n_vars = 4;
        let f1 = poly(n_vars, &[&[0], &[1]]);
        let f2 = poly(n_vars, &[&[2], &[3], &[]]);
        let c = F5Criterion::new(
            &[f1.clone(), f2.clone()],
            n_vars,
            2,
            all_variable_mask(n_vars),
        );
        assert!(c.prunes(1, f1.lt().unwrap().mask));
        assert!(c.prunes(1, f2.lt().unwrap().mask));
        assert!(c.prunes(0, f1.lt().unwrap().mask));
        assert!(!c.prunes(0, f2.lt().unwrap().mask));
        let f4 = matrix_f4_f2(&[f1.clone(), f2.clone()], n_vars, 2).unwrap();
        let (f5, report) = matrix_f5_f2(&[f1, f2], n_vars, 2).unwrap();
        assert_eq!(
            row_space_canonical(&f4, n_vars),
            row_space_canonical(&f5, n_vars)
        );
        assert_eq!(report.rows_pruned, 3);
    }

    /// The `F_2` criterion quoted from the semi-regularity literature —
    /// "prune `t·f_i` for `t ∈ LM⟨f_1..f_i⟩`", with `f_i` included to stand
    /// for the Frobenius syzygy — is unsound in the Boolean ring, because
    /// `x₁·(x₁x₂ + x₁ + x₂) = x₁` has a smaller leading monomial than its
    /// multiplier.  Pinned so nobody "simplifies" the criterion back to it.
    #[test]
    fn naive_f2_criterion_is_unsound_in_the_boolean_ring() {
        let n_vars = 3;
        let f = poly(n_vars, &[&[0, 1], &[0], &[1]]);
        let degree = 5;
        // Naive pruned set: LM of V_1(3) = span{s·f : deg s ≤ 1}.
        let level: Vec<F2BoolPoly> = monomials_up_to_mask(all_variable_mask(n_vars), 1)
            .into_iter()
            .map(|s| f.mul_mono(F2BoolMono::from_mask(s)))
            .collect();
        let naive: HashSet<u64> = row_space_canonical(&level, n_vars)
            .iter()
            .filter_map(|p| p.lt().map(|m| m.mask))
            .collect();
        assert!(
            naive.contains(&0b01) && naive.contains(&0b10),
            "x₀ and x₁ are naive LMs"
        );
        // Rows kept by the naive rule, and their span:
        let kept: Vec<F2BoolPoly> = monomials_up_to_mask(all_variable_mask(n_vars), degree - 2)
            .into_iter()
            .filter(|t| !naive.contains(t))
            .map(|t| f.mul_mono(F2BoolMono::from_mask(t)))
            .collect();
        let full = matrix_f4_f2(std::slice::from_ref(&f), n_vars, degree).unwrap();
        let naive_space = row_space_canonical(&kept, n_vars);
        assert!(
            naive_space.len() < full.len(),
            "the naive criterion must lose rank here, or the counterexample is gone"
        );
        // The sound criterion keeps the full row space.
        let (f5, report) = matrix_f5_f2(&[f], n_vars, degree).unwrap();
        assert_eq!(f5.len(), full.len());
        assert_eq!(
            row_space_canonical(&full, n_vars),
            row_space_canonical(&f5, n_vars)
        );
        assert!(report.rows_pruned > 0);
    }

    #[test]
    fn criterion_prunes_nothing_for_quadratics_at_degree_three() {
        // W_i(3) = V_{i−1}(1): with no linear generators there is nothing
        // to prune, so the F5 step is exactly the F4 step.
        let n_vars = 6;
        let mut seed = 42;
        let polys: Vec<F2BoolPoly> = (0..5)
            .map(|_| random_poly(n_vars, 2, 6, &mut seed))
            .filter(|p| poly_degree(p) == 2)
            .collect();
        let (_, report) = matrix_f5_f2(&polys, n_vars, 3).unwrap();
        assert_eq!(report.rows_pruned, 0);
    }
}
