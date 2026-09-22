//! # Inherited matrix-F4: the reduced Macaulay basis specialised down the splitting tree.
//!
//! The Koblitz decomposition oracle
//! ([`super::koblitz_groebner::solve_boolean_system`]) is a DPLL loop
//! around a Macaulay reduction: at every node it takes the current system,
//! builds the degree-`D` Macaulay matrix from scratch, reduces it, reads
//! off the refutation or the forced variables, and recurses on
//! `v := 0` and `v := 1`.  Every child rebuilds and re-reduces a matrix
//! that is a *specialisation* of the one its parent has just reduced.
//!
//! F5's organising idea is to never redo a reduction the algorithm can
//! predict.  Here the prediction is exact.  Substitution `φ_c : v ↦ c` is a
//! ring homomorphism of the Boolean ring, and the Macaulay rows of the
//! specialised system are the images of the parent's rows:
//!
//! ```text
//!     V(S|_{v=c})  =  φ_c( V(S) )  =  span{ φ_c(ρ) : ρ a reduced row of the parent }
//! ```
//!
//! (multipliers range over every variable; a generator whose degree
//! *drops* under `φ_c` gains multipliers — see [`ReducedBasis::specialise`]
//! for the completion rows).  So the child's row space is spanned by the
//! `rank` reduced rows of its parent, specialised — and most of them are
//! already in echelon form:
//!
//! - under `v := 0` every monomial containing `v` is deleted; a row whose
//!   leading monomial does not contain `v` keeps it;
//! - under `v := 1` every monomial `m ∋ v` folds into `m ∖ v`, which is
//!   *smaller* in any degree-compatible order, so again a leading monomial
//!   without `v` survives and cannot be cancelled.
//!
//! Only the rows whose **pivot contains `v`** lose their place and are
//! re-reduced against the rest — about `D/n` of them.  No Macaulay matrix
//! is built below the root, and no reduction is repeated.  Because the
//! row space is the same one the from-scratch step would compute, the
//! linear tail — hence every refutation, forced variable, split and
//! verdict — is identical; only the work per node changes.
//!
//! ## Accounting
//!
//! The stage's unit is 64-bit word XORs in the elimination.  Everything
//! this module does to a packed row is charged in that unit: the root
//! reduction through the shared kernel, every XOR of the re-reduction, and
//! the specialisation itself at one word operation per word read and per
//! word written — a charge the from-scratch path never pays, since its
//! matrix build is not counted.  Wall time is reported beside it.
//!
//! ## What this does not change
//!
//! Nothing about the algebra: same degree, same row space, same splitting
//! rule, same node budget.  Nothing about the decomposition oracle's cost
//! at `n = 131`, which is bounded by how many candidate tuples must be
//! ruled out and not by how one node's matrix is reduced.

use crate::cryptanalysis::koblitz_groebner::{
    all_variable_mask, echelon_f2_counted, macaulay_columns, macaulay_rows_monos,
    monomials_up_to_mask, pack_rows, rref_f2_counted,
};
use crate::cryptanalysis::pq_groebner_f2::{cmp_mono, F2BoolMono, F2BoolPoly};
use std::collections::HashMap;

/// Degree of a Boolean polynomial (`0` for a constant or zero).
fn poly_degree(p: &F2BoolPoly) -> u32 {
    p.terms.iter().map(|t| t.mask.count_ones()).max().unwrap_or(0)
}

/// Leading (lowest-index) set column of a packed row.
fn leading_column(row: &[u64]) -> Option<usize> {
    row.iter()
        .enumerate()
        .find(|(_, &w)| w != 0)
        .map(|(i, &w)| i * 64 + w.trailing_zeros() as usize)
}

/// Work counters for one basis operation, all in 64-bit word operations.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct InheritCost {
    /// XORs in the root reduction or in re-reducing displaced rows.
    pub reduce_word_ops: u64,
    /// Word reads and writes performed by specialisation (deleting or
    /// folding the columns of the assigned variable).
    pub specialise_word_ops: u64,
    /// Rows whose pivot contained the assigned variable and were re-reduced.
    pub displaced_rows: u64,
    /// Completion rows added because a generator's degree dropped.
    pub completion_rows: u64,
}

impl InheritCost {
    /// Everything charged, in word operations.
    pub fn word_ops(&self) -> u64 {
        self.reduce_word_ops + self.specialise_word_ops
    }
    fn add(&mut self, other: InheritCost) {
        self.reduce_word_ops += other.reduce_word_ops;
        self.specialise_word_ops += other.specialise_word_ops;
        self.displaced_rows += other.displaced_rows;
        self.completion_rows += other.completion_rows;
    }
}

/// The reduced degree-`degree` Macaulay row space of a system, kept as an
/// echelon basis with distinct leading monomials, together with the
/// system it belongs to.
#[derive(Clone, Debug)]
pub struct ReducedBasis {
    /// Macaulay degree.
    pub degree: u32,
    /// Number of Boolean variables of the ambient ring.
    pub n_vars: usize,
    /// The system whose Macaulay row space this is.
    pub system: Vec<F2BoolPoly>,
    /// Degree of each generator when the basis was last completed.
    generator_degrees: Vec<u32>,
    /// Variables specialised away on the path from the root.  Multipliers
    /// range over the complement: a multiplier containing an assigned
    /// variable is not the image of any parent row, and contributes
    /// nothing to the linear tail (see [`ReducedBasis::specialise`]).
    assigned: u64,
    /// Column monomials in descending order.
    columns: Vec<u64>,
    column_index: HashMap<u64, usize>,
    words: usize,
    /// Echelon rows; every row is nonzero and their leading columns are
    /// pairwise distinct.
    rows: Vec<Vec<u64>>,
    /// `pivot_of[c]` is the row leading at column `c`.
    pivot_of: Vec<Option<u32>>,
}

impl ReducedBasis {
    /// Build and reduce the Macaulay matrix of `system` at `degree` from
    /// scratch, with multipliers over every variable — the root of a
    /// splitting tree, or any node whose parent had no basis.  `None` if
    /// the matrix exceeds the F4 size caps.
    pub fn from_system(
        system: &[F2BoolPoly],
        n_vars: usize,
        degree: u32,
    ) -> Option<(Self, InheritCost)> {
        let system: Vec<F2BoolPoly> = system.iter().filter(|p| !p.is_zero()).cloned().collect();
        let generator_degrees: Vec<u32> = system.iter().map(poly_degree).collect();
        let rows_monos = macaulay_rows_monos(&system, n_vars, degree)?;
        let mut cost = InheritCost::default();
        if rows_monos.is_empty() {
            return Some((
                Self {
                    degree,
                    n_vars,
                    system,
                    generator_degrees,
                    assigned: 0,
                    columns: Vec::new(),
                    column_index: HashMap::new(),
                    words: 1,
                    rows: Vec::new(),
                    pivot_of: Vec::new(),
                },
                cost,
            ));
        }
        let columns = macaulay_columns(&rows_monos)?;
        let mut matrix = pack_rows(&rows_monos, &columns);
        // A basis needs distinct leading columns and nothing more, so an
        // echelon form without back-substitution would do.  The fully
        // reduced form costs more here but makes every later displaced row
        // reduce in exactly as many XORs as it has pivot bits, which pays
        // back over a deep tree.  The default mirrors the from-scratch
        // step's own kernel choice on the same shape — full RREF below 24
        // variables, echelon-only at and above it — so the root costs what
        // it always cost and every descendant is the saving.
        // `KIC_F4_INHERIT_ROOT=rref|ref` pins either for controls.
        let full = match std::env::var("KIC_F4_INHERIT_ROOT").as_deref() {
            Ok("rref") => true,
            Ok("ref") => false,
            _ => n_vars < 24,
        };
        let rank = if full {
            rref_f2_counted(&mut matrix, columns.len(), &mut cost.reduce_word_ops)
        } else {
            echelon_f2_counted(&mut matrix, columns.len(), &mut cost.reduce_word_ops)
        };
        matrix.truncate(rank);
        let words = columns.len().div_ceil(64).max(1);
        let column_index = columns.iter().enumerate().map(|(i, &m)| (m, i)).collect();
        let mut pivot_of = vec![None; columns.len()];
        for (r, row) in matrix.iter().enumerate() {
            let c = leading_column(row).expect("rank rows are nonzero");
            debug_assert!(pivot_of[c].is_none());
            pivot_of[c] = Some(r as u32);
        }
        Some((
            Self {
                degree,
                n_vars,
                system,
                generator_degrees,
                assigned: 0,
                columns,
                column_index,
                words,
                rows: matrix,
                pivot_of,
            },
            cost,
        ))
    }

    /// Variables specialised away since the root, as a mask.
    pub fn assigned(&self) -> u64 {
        self.assigned
    }

    /// Rank of the row space.
    pub fn rank(&self) -> usize {
        self.rows.len()
    }

    /// Number of columns.
    pub fn columns(&self) -> usize {
        self.columns.len()
    }

    /// Index of the first column of degree `≤ 1` (the linear tail).
    fn low_start(&self) -> usize {
        self.columns
            .iter()
            .position(|m| m.count_ones() <= 1)
            .unwrap_or(self.columns.len())
    }

    /// Insert a packed row, reducing its leading term against the existing
    /// pivots until it becomes a new pivot or vanishes.
    fn insert(&mut self, mut row: Vec<u64>, cost: &mut InheritCost) {
        loop {
            let Some(c) = leading_column(&row) else {
                return;
            };
            match self.pivot_of[c] {
                Some(r) => {
                    let w = c / 64;
                    let pivot = &self.rows[r as usize];
                    for (dst, &src) in row[w..].iter_mut().zip(&pivot[w..]) {
                        *dst ^= src;
                    }
                    cost.reduce_word_ops += (self.words - w) as u64;
                }
                None => {
                    self.pivot_of[c] = Some(self.rows.len() as u32);
                    self.rows.push(row);
                    return;
                }
            }
        }
    }

    /// The basis of `system|_{var = value}`: the specialised rows, with
    /// the displaced ones re-reduced, plus completion rows for any
    /// generator whose degree dropped.
    ///
    /// **Completion.**  The Macaulay matrix multiplies `f_i` by every
    /// monomial of degree `≤ D − deg f_i`.  If `deg f_i|_{v=c} < deg f_i`,
    /// the specialised system has multipliers of higher degree that no
    /// parent row maps to; those products are built and inserted so the
    /// row space is exactly the one a from-scratch step on the specialised
    /// system would reduce with multipliers over the unassigned variables.
    ///
    /// **Assigned variables as multipliers.**  The from-scratch step
    /// multiplies by every variable, assigned ones included.  Such a row
    /// `x·g`, with `x` occurring in no generator, lies in the linear tail
    /// only if `g ∈ {0, 1}`, and `g = 1` is already a refutation — so the
    /// two row spaces have the same tail and the solver behaves
    /// identically.
    pub fn specialise(&self, var: u32, value: bool) -> (Self, InheritCost) {
        let mut cost = InheritCost::default();
        let bit = 1u64 << var;

        // New system.
        let system: Vec<F2BoolPoly> = self
            .system
            .iter()
            .map(|p| substitute(p, var, value))
            .filter(|p| !p.is_zero())
            .collect();
        let new_degrees: Vec<u32> = system.iter().map(poly_degree).collect();

        // New column layout: images of the old columns.
        let mut new_columns: Vec<u64> = self
            .columns
            .iter()
            .filter(|&&m| value || m & bit == 0)
            .map(|&m| m & !bit)
            .collect();
        new_columns.sort_unstable();
        new_columns.dedup();
        new_columns.sort_by(|a, b| cmp_mono(F2BoolMono::from_mask(*a), F2BoolMono::from_mask(*b)).reverse());
        let column_index: HashMap<u64, usize> =
            new_columns.iter().enumerate().map(|(i, &m)| (m, i)).collect();
        let new_words = new_columns.len().div_ceil(64).max(1);
        // Old column -> new column (None: deleted).
        let remap: Vec<Option<usize>> = self
            .columns
            .iter()
            .map(|&m| {
                if !value && m & bit != 0 {
                    None
                } else {
                    Some(column_index[&(m & !bit)])
                }
            })
            .collect();

        let mut out = Self {
            degree: self.degree,
            n_vars: self.n_vars,
            system,
            generator_degrees: new_degrees.clone(),
            assigned: self.assigned | bit,
            columns: new_columns,
            column_index,
            words: new_words,
            rows: Vec::with_capacity(self.rows.len()),
            pivot_of: Vec::new(),
        };
        out.pivot_of = vec![None; out.columns.len()];

        // Specialise every row; keep the ones whose pivot survives, queue
        // the displaced ones.
        let mut displaced: Vec<Vec<u64>> = Vec::new();
        for (r, row) in self.rows.iter().enumerate() {
            let mut image = vec![0u64; new_words];
            for (wi, &word) in row.iter().enumerate() {
                let mut bits = word;
                while bits != 0 {
                    let b = bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    if let Some(c) = remap[wi * 64 + b] {
                        image[c / 64] ^= 1u64 << (c % 64);
                    }
                }
            }
            cost.specialise_word_ops += (self.words + new_words) as u64;
            let old_pivot = self.columns[leading_column(row).expect("basis rows are nonzero")];
            if old_pivot & bit == 0 {
                // Leading monomial survives and stays distinct.
                let c = remap[self.column_index[&old_pivot]].expect("kept column");
                debug_assert_eq!(leading_column(&image), Some(c));
                debug_assert!(out.pivot_of[c].is_none(), "row {r} collides");
                out.pivot_of[c] = Some(out.rows.len() as u32);
                out.rows.push(image);
            } else if image.iter().any(|&w| w != 0) {
                displaced.push(image);
            }
        }
        cost.displaced_rows = displaced.len() as u64;
        for row in displaced {
            out.insert(row, &mut cost);
        }

        // Completion rows for generators whose degree dropped, with
        // multipliers over the unassigned variables.
        let multiplier_mask = all_variable_mask(out.n_vars) & !out.assigned;
        let mut completion: Vec<Vec<u64>> = Vec::new();
        {
            // Pair each surviving generator with its pre-specialisation degree.
            let mut old = self
                .system
                .iter()
                .zip(&self.generator_degrees)
                .map(|(p, &d)| (substitute(p, var, value), d))
                .filter(|(p, _)| !p.is_zero());
            for (p, &new_degree) in out.system.iter().zip(&new_degrees) {
                let (_, old_degree) = old.next().expect("systems align");
                if new_degree == 0 || new_degree >= old_degree || old_degree > self.degree {
                    continue;
                }
                let old_gap = self.degree - old_degree;
                let new_gap = self.degree - new_degree;
                for t in monomials_up_to_mask(multiplier_mask, new_gap) {
                    if t.count_ones() <= old_gap {
                        continue;
                    }
                    let mut monos: Vec<u64> = p.terms.iter().map(|m| m.mask | t).collect();
                    monos.sort_unstable();
                    let mut row_monos = Vec::with_capacity(monos.len());
                    let mut i = 0;
                    while i < monos.len() {
                        let mut j = i;
                        while j < monos.len() && monos[j] == monos[i] {
                            j += 1;
                        }
                        if (j - i) % 2 == 1 {
                            row_monos.push(monos[i]);
                        }
                        i = j;
                    }
                    if !row_monos.is_empty() {
                        completion.push(row_monos);
                    }
                }
            }
        }
        if !completion.is_empty() {
            out.extend_columns(&completion, &mut cost);
            cost.completion_rows = completion.len() as u64;
            for monos in completion {
                let mut row = vec![0u64; out.words];
                for m in monos {
                    let c = out.column_index[&m];
                    row[c / 64] |= 1u64 << (c % 64);
                }
                out.insert(row, &mut cost);
            }
        }
        (out, cost)
    }

    /// Add any monomials of `rows_monos` missing from the layout, keeping
    /// the descending column order and re-packing the existing rows.
    fn extend_columns(&mut self, rows_monos: &[Vec<u64>], cost: &mut InheritCost) {
        let mut missing: Vec<u64> = rows_monos
            .iter()
            .flatten()
            .copied()
            .filter(|m| !self.column_index.contains_key(m))
            .collect();
        if missing.is_empty() {
            return;
        }
        missing.sort_unstable();
        missing.dedup();
        let mut columns = self.columns.clone();
        columns.extend(missing);
        columns.sort_by(|a, b| cmp_mono(F2BoolMono::from_mask(*a), F2BoolMono::from_mask(*b)).reverse());
        let column_index: HashMap<u64, usize> =
            columns.iter().enumerate().map(|(i, &m)| (m, i)).collect();
        let words = columns.len().div_ceil(64).max(1);
        let remap: Vec<usize> = self.columns.iter().map(|m| column_index[m]).collect();
        let mut rows = Vec::with_capacity(self.rows.len());
        let mut pivot_of = vec![None; columns.len()];
        for row in &self.rows {
            let mut image = vec![0u64; words];
            for (wi, &word) in row.iter().enumerate() {
                let mut bits = word;
                while bits != 0 {
                    let b = bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    let c = remap[wi * 64 + b];
                    image[c / 64] |= 1u64 << (c % 64);
                }
            }
            cost.specialise_word_ops += (self.words + words) as u64;
            let c = leading_column(&image).expect("nonzero");
            pivot_of[c] = Some(rows.len() as u32);
            rows.push(image);
        }
        self.columns = columns;
        self.column_index = column_index;
        self.words = words;
        self.rows = rows;
        self.pivot_of = pivot_of;
    }

    /// The linear tail — the intersection of the row space with the
    /// span of the degree-`≤ 1` monomials — in reduced row echelon form,
    /// as polynomials.  Exactly what the from-scratch step's readback
    /// consumes.
    pub fn linear_tail(&self, cost: &mut InheritCost) -> Vec<F2BoolPoly> {
        let low_start = self.low_start();
        let low_width = self.columns.len() - low_start;
        if low_width == 0 {
            return Vec::new();
        }
        let low_words = low_width.div_ceil(64).max(1);
        let mut low: Vec<Vec<u64>> = Vec::new();
        for row in &self.rows {
            let c = leading_column(row).expect("nonzero");
            if c < low_start {
                continue;
            }
            let mut packed = vec![0u64; low_words];
            for column in c..self.columns.len() {
                if row[column / 64] & (1u64 << (column % 64)) != 0 {
                    let k = column - low_start;
                    packed[k / 64] |= 1u64 << (k % 64);
                }
            }
            low.push(packed);
        }
        if low.is_empty() {
            return Vec::new();
        }
        let rank = rref_f2_counted(&mut low, low_width, &mut cost.reduce_word_ops);
        let n_vars = self.n_vars;
        low.iter()
            .take(rank)
            .map(|row| {
                let monos: Vec<F2BoolMono> = (0..low_width)
                    .filter(|k| row[k / 64] & (1u64 << (k % 64)) != 0)
                    .map(|k| F2BoolMono::from_mask(self.columns[low_start + k]))
                    .collect();
                F2BoolPoly::from_monos(monos, n_vars)
            })
            .collect()
    }

    /// The rows the splitting solver acts on: the constant `1`, and every
    /// `v` or `v + 1` of the linear tail.
    pub fn decisive_rows(&self, cost: &mut InheritCost) -> Vec<F2BoolPoly> {
        self.linear_tail(cost)
            .into_iter()
            .filter(|p| {
                (p.terms.len() == 1 && p.terms[0].mask == 0)
                    || (p.terms.len() <= 2
                        && p.terms.iter().filter(|t| t.mask.count_ones() == 1).count() == 1
                        && p.terms.iter().all(|t| t.mask.count_ones() <= 1))
            })
            .collect()
    }

    /// Specialise by several assignments in sequence.
    pub fn specialise_all(&self, assignments: &[(u32, bool)]) -> (Self, InheritCost) {
        let mut cost = InheritCost::default();
        let mut basis = self.clone();
        for &(v, val) in assignments {
            let (next, c) = basis.specialise(v, val);
            cost.add(c);
            basis = next;
        }
        (basis, cost)
    }
}

/// Specialise `p` by setting variable `var` to `value`.
pub fn substitute(p: &F2BoolPoly, var: u32, value: bool) -> F2BoolPoly {
    let bit = 1u64 << var;
    let mut monos = Vec::with_capacity(p.terms.len());
    for t in &p.terms {
        if t.mask & bit == 0 {
            monos.push(*t);
        } else if value {
            monos.push(F2BoolMono::from_mask(t.mask & !bit));
        }
    }
    F2BoolPoly::from_monos(monos, p.n_vars)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::koblitz_groebner::matrix_f4_f2;

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

    fn canonical(mut rows: Vec<F2BoolPoly>) -> Vec<F2BoolPoly> {
        rows.retain(|p| !p.is_zero());
        rows.sort_by(|a, b| format!("{a:?}").cmp(&format!("{b:?}")));
        rows
    }

    /// RREF of a list of polynomials, as a canonical set.
    fn rref_polys(polys: &[F2BoolPoly], n_vars: usize) -> Vec<F2BoolPoly> {
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
        canonical(
            matrix
                .iter()
                .take(rank)
                .map(|row| {
                    F2BoolPoly::from_monos(
                        (0..cols.len())
                            .filter(|c| row[c / 64] & (1u64 << (c % 64)) != 0)
                            .map(|c| F2BoolMono::from_mask(cols[c]))
                            .collect(),
                        n_vars,
                    )
                })
                .collect(),
        )
    }

    /// From-scratch Macaulay row space with multipliers over `mask`.
    fn masked_space(system: &[F2BoolPoly], n_vars: usize, degree: u32, mask: u64) -> Vec<F2BoolPoly> {
        use crate::cryptanalysis::koblitz_groebner::macaulay_rows_monos_with_mask;
        let rows = macaulay_rows_monos_with_mask(system, n_vars, degree, mask, None).unwrap();
        let polys: Vec<F2BoolPoly> = rows
            .iter()
            .map(|monos| {
                F2BoolPoly::from_monos(monos.iter().map(|&m| F2BoolMono::from_mask(m)).collect(), n_vars)
            })
            .collect();
        rref_polys(&polys, n_vars)
    }

    fn tail_of(polys: &[F2BoolPoly]) -> Vec<F2BoolPoly> {
        polys
            .iter()
            .filter(|p| p.terms.iter().all(|t| t.mask.count_ones() <= 1))
            .cloned()
            .collect()
    }

    /// The solver's view of a tail: a refutation if `1` is present,
    /// otherwise the tail itself.  The legacy step multiplies by assigned
    /// variables too, and once `1` is in the row space those junk rows put
    /// `x·1 = x` into its tail; the solver never reads past the `1`.
    fn solver_view(tail: &[F2BoolPoly], n_vars: usize) -> Vec<F2BoolPoly> {
        let tail = rref_polys(tail, n_vars);
        if tail.iter().any(|p| p.terms.len() == 1 && p.terms[0].mask == 0) {
            vec![F2BoolPoly::one(n_vars)]
        } else {
            tail
        }
    }

    /// The solver refutes a system containing the constant `1` before any
    /// reduction, so the exact row-space invariant is only claimed for
    /// systems it would actually reduce.
    fn has_constant(system: &[F2BoolPoly]) -> bool {
        system.iter().any(|p| p.terms.len() == 1 && p.terms[0].mask == 0)
    }

    fn basis_polys(b: &ReducedBasis) -> Vec<F2BoolPoly> {
        b.rows
            .iter()
            .map(|row| {
                F2BoolPoly::from_monos(
                    (0..b.columns.len())
                        .filter(|c| row[c / 64] & (1u64 << (c % 64)) != 0)
                        .map(|c| F2BoolMono::from_mask(b.columns[c]))
                        .collect(),
                    b.n_vars,
                )
            })
            .collect()
    }

    #[test]
    fn specialised_basis_spans_the_from_scratch_row_space() {
        let mut seed = 0x1234_5678_9abc_def1u64;
        for trial in 0..60 {
            let n_vars = 5 + trial % 4;
            let m = 3 + trial % 4;
            let system: Vec<F2BoolPoly> = (0..m)
                .map(|k| random_poly(n_vars, if k % 4 == 3 { 1 } else { 2 }, 5 + k, &mut seed))
                .filter(|p| poly_degree(p) >= 1)
                .collect();
            if system.is_empty() {
                continue;
            }
            for degree in 2..=3u32 {
                let (root, _) = ReducedBasis::from_system(&system, n_vars, degree).unwrap();
                // Root spans the F4 row space.
                let f4 = matrix_f4_f2(&system, n_vars, degree).unwrap();
                assert_eq!(rref_polys(&f4, n_vars), rref_polys(&basis_polys(&root), n_vars));
                // Every one-variable specialisation, both values, then a second one.
                let occurring = system.iter().flat_map(|p| p.terms.iter()).fold(0, |a, t| a | t.mask);
                let all = all_variable_mask(n_vars);
                for v in 0..n_vars as u32 {
                    if occurring & (1 << v) == 0 {
                        continue;
                    }
                    for value in [false, true] {
                        let (child, _) = root.specialise(v, value);
                        let child_system: Vec<F2BoolPoly> = system
                            .iter()
                            .map(|p| substitute(p, v, value))
                            .filter(|p| !p.is_zero())
                            .collect();
                        assert_eq!(child.system, child_system);
                        // Exact row space: multipliers over the unassigned variables.
                        if !has_constant(&child_system) {
                            assert_eq!(
                                masked_space(&child_system, n_vars, degree, all & !(1 << v)),
                                rref_polys(&basis_polys(&child), n_vars),
                                "trial {trial} degree {degree} v{v}={value}: row space differs"
                            );
                        }
                        // Tail agreement with the legacy all-variable step's readback.
                        let legacy = matrix_f4_f2(&child_system, n_vars, degree).unwrap_or_default();
                        let mut cost = InheritCost::default();
                        let tail = child.linear_tail(&mut cost);
                        assert_eq!(
                            solver_view(&tail_of(&legacy), n_vars),
                            solver_view(&tail, n_vars),
                            "trial {trial} degree {degree} v{v}={value}: tail differs"
                        );
                        // One more level.
                        let w = (v + 1) % n_vars as u32;
                        let (grandchild, _) = child.specialise(w, !value);
                        let gc_system: Vec<F2BoolPoly> = child_system
                            .iter()
                            .map(|p| substitute(p, w, !value))
                            .filter(|p| !p.is_zero())
                            .collect();
                        if !has_constant(&gc_system) {
                            assert_eq!(
                                masked_space(&gc_system, n_vars, degree, all & !(1 << v) & !(1 << w)),
                                rref_polys(&basis_polys(&grandchild), n_vars),
                                "trial {trial} degree {degree} v{v}={value} w{w}: grandchild differs"
                            );
                        }
                        let legacy = matrix_f4_f2(&gc_system, n_vars, degree).unwrap_or_default();
                        let tail = grandchild.linear_tail(&mut cost);
                        assert_eq!(
                            solver_view(&tail_of(&legacy), n_vars),
                            solver_view(&tail, n_vars),
                            "trial {trial} degree {degree} v{v}={value} w{w}: grandchild tail differs"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn degree_drop_completion_rows_are_added() {
        // f = x0·x1 + x2: setting x0 = 0 makes it linear, so degree-2
        // multipliers appear in the from-scratch matrix at degree 3.
        let n_vars = 4;
        let f = F2BoolPoly::from_monos(
            vec![F2BoolMono::from_mask(0b0011), F2BoolMono::from_mask(0b0100)],
            n_vars,
        );
        let g = F2BoolPoly::from_monos(
            vec![F2BoolMono::from_mask(0b1100), F2BoolMono::from_mask(0b0010), F2BoolMono::one()],
            n_vars,
        );
        let (root, _) = ReducedBasis::from_system(&[f, g], n_vars, 3).unwrap();
        let (child, cost) = root.specialise(0, false);
        assert!(cost.completion_rows > 0, "f dropped to degree 1; completion expected");
        assert_eq!(
            masked_space(&child.system, n_vars, 3, all_variable_mask(n_vars) & !1),
            rref_polys(&basis_polys(&child), n_vars)
        );
    }
}
