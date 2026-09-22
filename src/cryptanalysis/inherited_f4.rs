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
//! (up to multipliers that cannot reach the linear tail, and up to a
//! generator whose degree *drops* under `φ_c` gaining multipliers — see
//! [`ReducedBasis::specialise`] for both).  So the child's row space is
//! spanned by the `rank` reduced rows of its parent, specialised — and most
//! of them are already in echelon form:
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
//! ## What decides the cost: which variable the solver splits on
//!
//! Re-reducing a displaced row costs one XOR per pivot it crosses, so the
//! method's cost is the number of displaced rows times the density of the
//! basis.  Both depend on the **split variable**.  A reduced row's pivot
//! is its largest monomial under DegRevLex, so pivots are biased toward
//! the large variables; splitting on the lowest-indexed free variable (the
//! historical `LowestFree` rule) picks exactly those, and deep in the tree
//! half the basis is displaced per level.  Splitting on the smallest free
//! variable ([`super::koblitz_groebner::SplitRule::HighestFree`]) displaces
//! the fewest rows: the engine's work halves on every quadratic rung and
//! its loss on the chained cubic `m = 3` systems — whose 40%-rank-deficient
//! matrices make the reduced rows dense — turns into a `1.9×` win on the
//! deep `K_0/2^15` cell (`RESEARCH_INHERITED_F4.md` §3.5).  Restoring
//! reduced form after each level, cost-triggered fallbacks and
//! drift-triggered rebuilds were all tried and rejected; each cost more
//! than it saved.
//!
//! ## What this does not change
//!
//! Nothing about the algebra: same degree, same row space, same splitting
//! rule, same node budget.  Nothing about the decomposition oracle's cost
//! at `n = 131`, which is bounded by how many candidate tuples must be
//! ruled out and not by how one node's matrix is reduced.

use crate::cryptanalysis::koblitz_groebner::{
    all_variable_mask, echelon_f2_counted, macaulay_columns, macaulay_rows_monos_with_mask,
    monomials_up_to_mask, pack_rows, rref_f2_counted,
};
use crate::cryptanalysis::pq_groebner_f2::{cmp_mono, F2BoolMono, F2BoolPoly};
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

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
    /// XORs in the root reduction, in re-reducing displaced rows, and in
    /// restoring reduced form ([`ReducedBasis::reduce_fully`]).
    pub reduce_word_ops: u64,
    /// The part of `reduce_word_ops` spent restoring reduced form.
    pub rref_word_ops: u64,
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
        self.rref_word_ops += other.rref_word_ops;
        self.specialise_word_ops += other.specialise_word_ops;
        self.displaced_rows += other.displaced_rows;
        self.completion_rows += other.completion_rows;
    }
}

/// One change of column layout: how the columns of one epoch map to the
/// columns of the next.  `None` is a deleted column (`v := 0`); two
/// columns may map to the same target (`v := 1` folds `m ∋ v` onto
/// `m ∖ v`), and a materialisation XORs them together.
#[derive(Debug)]
struct LayoutStep {
    /// Words per row in the epoch this step maps *from*.
    from_words: usize,
    map: Rc<[Option<u32>]>,
    /// Composed maps from earlier epochs to the epoch this step maps *to*,
    /// built on demand.  Shared with every basis below this step, so a
    /// child composes one step on top of what its parent already built.
    to_here: RefCell<HashMap<u32, Rc<[Option<u32>]>>>,
}

/// A basis row: its content, valid in the layout of epoch `version`.
///
/// Rows are shared by reference between a basis and the children
/// specialised from it; a row is rewritten into a newer layout — and
/// charged for it — only when something needs its content there.
#[derive(Clone, Debug)]
struct LazyRow {
    version: u32,
    data: Rc<Vec<u64>>,
}

/// The reduced degree-`degree` Macaulay row space of a system, kept as an
/// echelon basis with distinct leading monomials, together with the
/// system it belongs to.
///
/// **Lazy materialisation.**  Only a fraction of the rows a level keeps
/// are ever touched at that level — displaced, hit as a pivot while a
/// displaced row reduces, or read for the tail; on `K_1/2^23` three
/// quarters are not.  So a kept row is not rewritten into the child's
/// layout: the child shares the parent's row by reference and records
/// the layout step, and a row is materialised only when needed, through
/// the composition of every step since it was last written, in one pass
/// charged once.  Its pivot column is tracked through the steps meanwhile,
/// which is all `specialise` needs to decide whether it stays a pivot.
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
    /// Variables specialised away on the path from the root.
    assigned: u64,
    /// Column monomials of the current layout, in descending order.
    columns: Vec<u64>,
    /// Monomial → column, built only when a row has to be packed from
    /// monomials (completion), which most levels never do.
    column_index: Option<HashMap<u64, usize>>,
    words: usize,
    /// Layout steps since the root: `history[e]` maps epoch `e` to `e + 1`.
    /// The current epoch is `history.len()`.
    history: Vec<Rc<LayoutStep>>,
    /// Echelon rows; every row is nonzero and their leading columns are
    /// pairwise distinct.
    rows: Vec<LazyRow>,
    /// `pivot_col[r]` is row `r`'s leading column in the current layout.
    pivot_col: Vec<u32>,
    /// `pivot_of[c]` is the row leading at column `c`.
    pivot_of: Vec<Option<u32>>,
}

impl ReducedBasis {
    /// Build and reduce the Macaulay matrix of `system` at `degree` from
    /// scratch — the root of a splitting tree, or any node whose parent had
    /// no basis.  Multipliers range over the variables occurring in the
    /// system, the from-scratch step's own active-multiplier policy: a
    /// multiplier containing a variable that occurs nowhere adds rows that
    /// cannot reach the linear tail (see [`ReducedBasis::specialise`]).
    /// `None` if the matrix exceeds the F4 size caps.
    pub fn from_system(
        system: &[F2BoolPoly],
        n_vars: usize,
        degree: u32,
    ) -> Option<(Self, InheritCost)> {
        let system: Vec<F2BoolPoly> = system.iter().filter(|p| !p.is_zero()).cloned().collect();
        let generator_degrees: Vec<u32> = system.iter().map(poly_degree).collect();
        let occurring = system
            .iter()
            .flat_map(|p| p.terms.iter())
            .fold(0u64, |acc, t| acc | t.mask)
            & all_variable_mask(n_vars);
        let rows_monos =
            macaulay_rows_monos_with_mask(&system, n_vars, degree, occurring, None)?;
        let mut cost = InheritCost::default();
        let empty = |columns: Vec<u64>| Self {
            degree,
            n_vars,
            system: system.clone(),
            generator_degrees: generator_degrees.clone(),
            assigned: 0,
            column_index: None,
            words: columns.len().div_ceil(64).max(1),
            pivot_of: vec![None; columns.len()],
            columns,
            history: Vec::new(),
            rows: Vec::new(),
            pivot_col: Vec::new(),
        };
        if rows_monos.is_empty() {
            return Some((empty(Vec::new()), cost));
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
        let mut out = empty(columns);
        for row in matrix {
            let c = leading_column(&row).expect("rank rows are nonzero");
            debug_assert!(out.pivot_of[c].is_none());
            out.pivot_of[c] = Some(out.rows.len() as u32);
            out.pivot_col.push(c as u32);
            out.rows.push(LazyRow {
                version: 0,
                data: Rc::new(row),
            });
        }
        Some((out, cost))
    }

    /// Monomial → current column, built on first use after a layout change.
    fn column_index(&mut self) -> &HashMap<u64, usize> {
        if self.column_index.is_none() {
            self.column_index =
                Some(self.columns.iter().enumerate().map(|(i, &m)| (m, i)).collect());
        }
        self.column_index.as_ref().expect("just built")
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

    /// Depth of this basis below its root: the number of variables
    /// specialised away.
    pub fn depth(&self) -> u32 {
        self.assigned.count_ones()
    }

    /// The current epoch: the number of layout steps since the root.
    fn epoch(&self) -> u32 {
        self.history.len() as u32
    }

    /// Index of the first column of degree `≤ 1` (the linear tail).
    fn low_start(&self) -> usize {
        self.columns
            .iter()
            .position(|m| m.count_ones() <= 1)
            .unwrap_or(self.columns.len())
    }

    /// The composed map from the columns of epoch `from` to those of epoch
    /// `to`.  Memoised on the step that maps into `to`, which every basis
    /// below that step shares, so a node composes at most one step on top
    /// of what its ancestors already built.  Index bookkeeping, not
    /// charged, as the from-scratch path's column indexing never was.
    fn composite_to(&self, to: u32, from: u32) -> Rc<[Option<u32>]> {
        debug_assert!(from < to);
        let step = &self.history[(to - 1) as usize];
        if from + 1 == to {
            return step.map.clone();
        }
        if let Some(c) = step.to_here.borrow().get(&from) {
            return c.clone();
        }
        let prev = self.composite_to(to - 1, from);
        let composed: Rc<[Option<u32>]> = prev
            .iter()
            .map(|c| c.and_then(|c| step.map[c as usize]))
            .collect();
        step.to_here.borrow_mut().insert(from, composed.clone());
        composed
    }

    /// The composed map from the columns of epoch `from` to the current
    /// layout.
    fn composite(&self, from: u32) -> Rc<[Option<u32>]> {
        self.composite_to(self.epoch(), from)
    }

    /// Rewrite `row` from the layout of epoch `version` into the current
    /// one.  Charged one word operation per word read and per word written.
    fn materialised(&self, row: &LazyRow, cost: &mut InheritCost) -> Rc<Vec<u64>> {
        debug_assert!(row.version < self.epoch());
        let map = self.composite(row.version);
        let from_words = self.history[row.version as usize].from_words;
        let mut out = vec![0u64; self.words];
        for (wi, &word) in row.data.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                let b = bits.trailing_zeros() as usize;
                bits &= bits - 1;
                if let Some(c) = map[wi * 64 + b] {
                    out[c as usize / 64] ^= 1u64 << (c % 64);
                }
            }
        }
        cost.specialise_word_ops += (from_words + self.words) as u64;
        Rc::new(out)
    }

    /// Bring row `r` into the current layout if it is not there yet.
    fn ensure_current(&mut self, r: usize, cost: &mut InheritCost) {
        if self.rows[r].version != self.epoch() {
            let data = self.materialised(&self.rows[r], cost);
            self.rows[r] = LazyRow {
                version: self.epoch(),
                data,
            };
        }
    }

    /// Bring every row into the current layout.
    pub fn materialise_all(&mut self, cost: &mut InheritCost) {
        for r in 0..self.rows.len() {
            self.ensure_current(r, cost);
        }
    }

    /// Row `r`'s content in the current layout, materialising it if needed.
    fn current_row(&mut self, r: usize, cost: &mut InheritCost) -> Rc<Vec<u64>> {
        self.ensure_current(r, cost);
        self.rows[r].data.clone()
    }

    /// Restore reduced row echelon form: clear every pivot column from
    /// every row but its own.
    ///
    /// One pass in **decreasing** pivot-column order suffices.  A row's
    /// bits all sit at or after its leading column, so when pivot column
    /// `c` is processed every pivot column beyond `c` has already been
    /// cleared from its row, and XORing it into another row introduces no
    /// pivot-column bit.  Only XORs are charged, as in the from-scratch
    /// kernels; the pivot-column tests are not.
    ///
    /// Why it is worth paying: after a specialisation the basis is in
    /// echelon form but not reduced — a fold `m ∋ v ↦ m ∖ v` can land on
    /// another row's pivot column, and a freshly inserted pivot has never
    /// been cleared from the rows above it.  A displaced row reduced
    /// against such a basis picks up those bits and cascades; against a
    /// reduced basis it clears one pivot column per XOR and stops.
    pub fn reduce_fully(&mut self, cost: &mut InheritCost) {
        self.materialise_all(cost);
        let words = self.words;
        let epoch = self.epoch();
        for c in (0..self.pivot_of.len()).rev() {
            let Some(r) = self.pivot_of[c] else {
                continue;
            };
            let r = r as usize;
            let (w, bit) = (c / 64, 1u64 << (c % 64));
            let pivot = self.rows[r].data.clone();
            for (i, row) in self.rows.iter_mut().enumerate() {
                if i != r && row.data[w] & bit != 0 {
                    let mut data = row.data.to_vec();
                    for (dst, &src) in data[w..].iter_mut().zip(&pivot[w..]) {
                        *dst ^= src;
                    }
                    *row = LazyRow {
                        version: epoch,
                        data: Rc::new(data),
                    };
                    cost.reduce_word_ops += (words - w) as u64;
                    cost.rref_word_ops += (words - w) as u64;
                }
            }
        }
    }

    /// Should a basis at `depth` be restored to reduced form after its
    /// specialisation?  `KIC_F4_INHERIT_RREF` names the policy: `0` never
    /// (the default), `k ≥ 1` every `k` levels.
    ///
    /// A retained negative control.  Restoring reduced form does cut the
    /// cascade — on `K_1/2^23` a displaced row then hits `9.8` pivots
    /// instead of `40.5` — but the pass itself costs `72 M` word XORs
    /// against the `49 M` it saves, because a reduced basis is *denser*
    /// than an echelon one and every specialisation refills some twenty
    /// pivot columns per row (`RESEARCH_INHERITED_F4.md` §3.5).
    fn rref_every() -> u32 {
        std::env::var("KIC_F4_INHERIT_RREF")
            .ok()
            .and_then(|v| v.parse::<u32>().ok())
            .unwrap_or(0)
    }

    /// Insert a row given in the current layout, reducing its leading term
    /// against the existing pivots until it becomes a new pivot or vanishes.
    /// Every pivot it hits is materialised first.
    fn insert(&mut self, mut row: Vec<u64>, cost: &mut InheritCost) {
        loop {
            let Some(c) = leading_column(&row) else {
                return;
            };
            match self.pivot_of[c] {
                Some(r) => {
                    let w = c / 64;
                    let pivot = self.current_row(r as usize, cost);
                    for (dst, &src) in row[w..].iter_mut().zip(&pivot[w..]) {
                        *dst ^= src;
                    }
                    cost.reduce_word_ops += (self.words - w) as u64;
                }
                None => {
                    self.pivot_of[c] = Some(self.rows.len() as u32);
                    self.pivot_col.push(c as u32);
                    self.rows.push(LazyRow {
                        version: self.epoch(),
                        data: Rc::new(row),
                    });
                    return;
                }
            }
        }
    }

    /// Adopt a new column layout: record the step, remap every pivot
    /// column and rebuild the pivot index.  Rows are left where they are.
    fn adopt_layout(&mut self, columns: Vec<u64>, map: Vec<Option<u32>>) {
        let step = LayoutStep {
            from_words: self.words,
            map: map.into(),
            to_here: RefCell::new(HashMap::new()),
        };
        self.column_index = None;
        self.words = columns.len().div_ceil(64).max(1);
        self.pivot_of = vec![None; columns.len()];
        for (r, pc) in self.pivot_col.iter_mut().enumerate() {
            let c = step.map[*pc as usize].expect("a kept pivot column survives the step");
            debug_assert!(self.pivot_of[c as usize].is_none(), "row {r} collides");
            self.pivot_of[c as usize] = Some(r as u32);
            *pc = c;
        }
        self.columns = columns;
        self.history.push(Rc::new(step));
    }

    /// The basis of `system|_{var = value}`: the kept rows shared by
    /// reference, the displaced ones materialised and re-reduced, plus
    /// completion rows for any generator whose degree dropped.
    ///
    /// **Completion.**  The Macaulay matrix multiplies `f_i` by every
    /// monomial of degree `≤ D − deg f_i`.  If `deg f_i|_{v=c} < deg f_i`,
    /// the specialised system has multipliers of higher degree that no
    /// parent row maps to; those products are built and inserted, with
    /// multipliers over the variables occurring in the specialised system.
    ///
    /// **Which row space, exactly.**  Write `V_M(S)` for the Macaulay row
    /// space of `S` with multipliers over the variable set `M`.  The
    /// inherited rows span `V_{unassigned}` for the generators whose degree
    /// did not drop and the completion adds `V_{occurring}` for the rest,
    /// so the basis lies between `V_{occurring}(S')` and `V_{all}(S')`, the
    /// latter being what the from-scratch step reduces.  Every space in
    /// that sandwich has the **same linear tail**: a row `x·g` with `x`
    /// occurring in no generator lies in the tail only if `g ∈ {0, 1}`,
    /// and `g = 1` is already a refutation.  The solver reads nothing but
    /// the tail, so it behaves identically.
    pub fn specialise(&self, var: u32, value: bool) -> (Self, InheritCost) {
        let mut cost = InheritCost::default();
        let bit = 1u64 << var;

        // New system, with each survivor's degree before and after.
        let mut system = Vec::with_capacity(self.system.len());
        let mut new_degrees = Vec::with_capacity(self.system.len());
        let mut dropped: Vec<(usize, u32, u32)> = Vec::new();
        for (p, &old_degree) in self.system.iter().zip(&self.generator_degrees) {
            let q = substitute(p, var, value);
            if q.is_zero() {
                continue;
            }
            let new_degree = poly_degree(&q);
            if new_degree != 0 && new_degree < old_degree && old_degree <= self.degree {
                dropped.push((system.len(), old_degree, new_degree));
            }
            system.push(q);
            new_degrees.push(new_degree);
        }

        // New column layout: the images of the old columns.  Deleting the
        // columns that contain `v` keeps the rest in order; folding
        // `m ∋ v ↦ m ∖ v` keeps the folded ones in order among themselves
        // (every degree drops by one, and the symmetric difference of two
        // of them is unchanged), so the new layout is a merge of two sorted
        // lists and the old → new map falls out of the merge.
        let n_old = self.columns.len();
        let mut map: Vec<Option<u32>> = vec![None; n_old];
        let mut new_columns: Vec<u64> = Vec::with_capacity(n_old);
        if value {
            let (mut i, mut j) = (0usize, 0usize);
            let next_keep = |i: &mut usize| -> Option<usize> {
                while *i < n_old && self.columns[*i] & bit != 0 {
                    *i += 1;
                }
                (*i < n_old).then_some(*i)
            };
            let next_fold = |j: &mut usize| -> Option<usize> {
                while *j < n_old && self.columns[*j] & bit == 0 {
                    *j += 1;
                }
                (*j < n_old).then_some(*j)
            };
            loop {
                let k = next_keep(&mut i);
                let f = next_fold(&mut j);
                let target = new_columns.len() as u32;
                match (k, f) {
                    (None, None) => break,
                    (Some(k), None) => {
                        new_columns.push(self.columns[k]);
                        map[k] = Some(target);
                        i += 1;
                    }
                    (None, Some(f)) => {
                        new_columns.push(self.columns[f] & !bit);
                        map[f] = Some(target);
                        j += 1;
                    }
                    (Some(k), Some(f)) => {
                        let (mk, mf) = (self.columns[k], self.columns[f] & !bit);
                        match cmp_mono(F2BoolMono::from_mask(mk), F2BoolMono::from_mask(mf)) {
                            std::cmp::Ordering::Greater => {
                                new_columns.push(mk);
                                map[k] = Some(target);
                                i += 1;
                            }
                            std::cmp::Ordering::Less => {
                                new_columns.push(mf);
                                map[f] = Some(target);
                                j += 1;
                            }
                            std::cmp::Ordering::Equal => {
                                new_columns.push(mk);
                                map[k] = Some(target);
                                map[f] = Some(target);
                                i += 1;
                                j += 1;
                            }
                        }
                    }
                }
            }
        } else {
            for (k, &m) in self.columns.iter().enumerate() {
                if m & bit == 0 {
                    map[k] = Some(new_columns.len() as u32);
                    new_columns.push(m);
                }
            }
        }
        debug_assert!(new_columns.windows(2).all(|w| {
            cmp_mono(F2BoolMono::from_mask(w[0]), F2BoolMono::from_mask(w[1]))
                == std::cmp::Ordering::Greater
        }));

        // The child starts as a copy of the parent's bookkeeping — rows by
        // reference — minus the displaced rows, then adopts the new layout.
        let mut out = Self {
            degree: self.degree,
            n_vars: self.n_vars,
            system,
            generator_degrees: new_degrees,
            assigned: self.assigned | bit,
            columns: self.columns.clone(),
            column_index: None,
            words: self.words,
            history: self.history.clone(),
            rows: Vec::with_capacity(self.rows.len()),
            pivot_col: Vec::with_capacity(self.rows.len()),
            pivot_of: Vec::new(),
        };
        let mut displaced: Vec<LazyRow> = Vec::new();
        for (row, &pc) in self.rows.iter().zip(&self.pivot_col) {
            if self.columns[pc as usize] & bit == 0 {
                // Leading monomial survives and stays distinct; the row's
                // content is not needed here.
                out.rows.push(row.clone());
                out.pivot_col.push(pc);
            } else {
                displaced.push(row.clone());
            }
        }
        out.adopt_layout(new_columns, map);
        cost.displaced_rows = displaced.len() as u64;
        for row in displaced {
            let image = out.materialised(&row, &mut cost);
            if image.iter().any(|&w| w != 0) {
                out.insert(Rc::try_unwrap(image).unwrap_or_else(|rc| (*rc).clone()), &mut cost);
            }
        }

        // Completion rows for generators whose degree dropped, with
        // multipliers over the variables occurring in the specialised
        // system (the from-scratch step's own active-multiplier policy):
        // a multiplier containing a variable that occurs nowhere adds only
        // rows `x·g` that cannot reach the linear tail.
        if !dropped.is_empty() {
            let multiplier_mask = out
                .system
                .iter()
                .flat_map(|p| p.terms.iter())
                .fold(0u64, |acc, t| acc | t.mask)
                & all_variable_mask(out.n_vars);
            let mut completion: Vec<Vec<u64>> = Vec::new();
            for &(index, old_degree, new_degree) in &dropped {
                let p = &out.system[index];
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
            if !completion.is_empty() {
                out.extend_columns(&completion);
                cost.completion_rows = completion.len() as u64;
                for monos in completion {
                    let words = out.words;
                    let index = out.column_index();
                    let mut row = vec![0u64; words];
                    for m in monos {
                        let c = index[&m];
                        row[c / 64] |= 1u64 << (c % 64);
                    }
                    out.insert(row, &mut cost);
                }
            }
        }
        let every = Self::rref_every();
        if every > 0 && out.depth() % every == 0 {
            out.reduce_fully(&mut cost);
        }
        (out, cost)
    }

    /// Add any monomials of `rows_monos` missing from the layout, keeping
    /// the descending column order.  A layout step like any other: rows
    /// are not touched.
    fn extend_columns(&mut self, rows_monos: &[Vec<u64>]) {
        let index = self.column_index();
        let mut missing: Vec<u64> = rows_monos
            .iter()
            .flatten()
            .copied()
            .filter(|m| !index.contains_key(m))
            .collect();
        if missing.is_empty() {
            return;
        }
        missing.sort_unstable();
        missing.dedup();
        let mut columns = self.columns.clone();
        columns.extend(missing);
        columns.sort_by(|a, b| cmp_mono(F2BoolMono::from_mask(*a), F2BoolMono::from_mask(*b)).reverse());
        let index: HashMap<u64, usize> = columns.iter().enumerate().map(|(i, &m)| (m, i)).collect();
        let map: Vec<Option<u32>> = self.columns.iter().map(|m| Some(index[m] as u32)).collect();
        self.adopt_layout(columns, map);
    }

    /// The linear tail — the intersection of the row space with the
    /// span of the degree-`≤ 1` monomials — in reduced row echelon form,
    /// as polynomials.  Exactly what the from-scratch step's readback
    /// consumes.  Only the rows leading in the tail are materialised.
    pub fn linear_tail(&mut self, cost: &mut InheritCost) -> Vec<F2BoolPoly> {
        let low_start = self.low_start();
        let low_width = self.columns.len() - low_start;
        if low_width == 0 {
            return Vec::new();
        }
        let low_words = low_width.div_ceil(64).max(1);
        let tail_rows: Vec<usize> = (0..self.rows.len())
            .filter(|&r| self.pivot_col[r] as usize >= low_start)
            .collect();
        let mut low: Vec<Vec<u64>> = Vec::with_capacity(tail_rows.len());
        for r in tail_rows {
            let row = self.current_row(r, cost);
            let c = self.pivot_col[r] as usize;
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
    pub fn decisive_rows(&mut self, cost: &mut InheritCost) -> Vec<F2BoolPoly> {
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

    fn occurring(system: &[F2BoolPoly]) -> u64 {
        system.iter().flat_map(|p| p.terms.iter()).fold(0, |a, t| a | t.mask)
    }

    /// Does every polynomial of `inner` reduce to zero against the RREF
    /// `outer`?  (`inner ⊆ span(outer)`.)
    fn contained(inner: &[F2BoolPoly], outer: &[F2BoolPoly], n_vars: usize) -> bool {
        let joint = rref_polys(&[outer, inner].concat(), n_vars);
        joint.len() == rref_polys(outer, n_vars).len()
    }

    /// The sandwich `V_occurring(S) ⊆ basis ⊆ V_unassigned(S)` the module
    /// documents, for a node reached by assigning `assigned`.
    fn assert_sandwich(basis: &ReducedBasis, system: &[F2BoolPoly], n_vars: usize, assigned: u64, what: &str) {
        let rows = rref_polys(&basis_polys(basis), n_vars);
        let lower = masked_space(system, n_vars, basis.degree, occurring(system) & !assigned);
        let upper = masked_space(system, n_vars, basis.degree, all_variable_mask(n_vars) & !assigned);
        assert!(contained(&lower, &rows, n_vars), "{what}: V_occurring not contained in the basis");
        assert!(contained(&rows, &upper, n_vars), "{what}: basis not contained in V_unassigned");
    }

    /// Every row, materialised into the current layout, as polynomials.
    fn basis_polys(b: &ReducedBasis) -> Vec<F2BoolPoly> {
        let mut b = b.clone();
        let mut cost = InheritCost::default();
        b.materialise_all(&mut cost);
        b.rows
            .iter()
            .map(|row| {
                F2BoolPoly::from_monos(
                    (0..b.columns.len())
                        .filter(|c| row.data[c / 64] & (1u64 << (c % 64)) != 0)
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
                let (mut root, _) = ReducedBasis::from_system(&system, n_vars, degree).unwrap();
                // Root: the sandwich with nothing assigned, and the legacy tail.
                assert_sandwich(&root, &system, n_vars, 0, &format!("trial {trial} degree {degree} root"));
                let f4 = matrix_f4_f2(&system, n_vars, degree).unwrap();
                let mut root_cost = InheritCost::default();
                assert_eq!(
                    solver_view(&tail_of(&f4), n_vars),
                    solver_view(&root.linear_tail(&mut root_cost), n_vars),
                    "trial {trial} degree {degree}: root tail differs"
                );
                // Every one-variable specialisation, both values, then a second one.
                let occurring = system.iter().flat_map(|p| p.terms.iter()).fold(0, |a, t| a | t.mask);
                let all = all_variable_mask(n_vars);
                for v in 0..n_vars as u32 {
                    if occurring & (1 << v) == 0 {
                        continue;
                    }
                    for value in [false, true] {
                        let (mut child, _) = root.specialise(v, value);
                        let child_system: Vec<F2BoolPoly> = system
                            .iter()
                            .map(|p| substitute(p, v, value))
                            .filter(|p| !p.is_zero())
                            .collect();
                        assert_eq!(child.system, child_system);
                        if !has_constant(&child_system) {
                            assert_sandwich(
                                &child,
                                &child_system,
                                n_vars,
                                1 << v,
                                &format!("trial {trial} degree {degree} v{v}={value}"),
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
                        let (mut grandchild, _) = child.specialise(w, !value);
                        let gc_system: Vec<F2BoolPoly> = child_system
                            .iter()
                            .map(|p| substitute(p, w, !value))
                            .filter(|p| !p.is_zero())
                            .collect();
                        if !has_constant(&gc_system) {
                            assert_sandwich(
                                &grandchild,
                                &gc_system,
                                n_vars,
                                (1 << v) | (1 << w),
                                &format!("trial {trial} degree {degree} v{v}={value} w{w}: grandchild"),
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
        let system = child.system.clone();
        assert_sandwich(&child, &system, n_vars, 1, "degree drop");
        // With every unassigned variable occurring, the sandwich is an equality.
        assert_eq!(occurring(&system), all_variable_mask(n_vars) & !1);
        assert_eq!(
            masked_space(&system, n_vars, 3, all_variable_mask(n_vars) & !1),
            rref_polys(&basis_polys(&child), n_vars)
        );
    }
}
