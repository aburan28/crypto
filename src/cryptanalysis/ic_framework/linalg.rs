//! # Relation matrices, as plug-ins.
//!
//! The relation loop takes any [`RelationSolver`]; this module ships a
//! second one so the linear-algebra stage is a real lever rather than a
//! fixed cost.
//!
//! | name | method | storage | pivot | work unit |
//! |:--|:--|:--|:--|:--|
//! | `incremental-gauss` | dense reduced row echelon, maintained as rows arrive | `rank × |F|` entries | leftmost non-zero | `row_ops` |
//! | `incremental-gauss-full-rank` | same dense matrix, continue after target pin until every column is independent | `rank × |F|` entries | leftmost non-zero | `row_ops` |
//! | `incremental-gauss-full-rank-checked` | same full-rank matrix, then point-check one representative per factor-base column | `rank × |F|` entries | leftmost non-zero | `row_ops` |
//! | `structured-gauss` | sparse reduced row echelon, maintained as rows arrive | the non-zeros | lightest column (Markowitz) | `row_ops` |
//! | `wiedemann` | rows kept as they arrive; singleton-filtered core solved in one batch by Wiedemann, checked exactly | the non-zeros | none | `row_ops` |
//!
//! Both count a multiply-subtract on a **non-zero** entry as one
//! `row_op`, so the two `work` figures are directly comparable and
//! differ only through the pivot choice: pivoting on the lightest
//! column limits fill-in, which on the `m`-non-zero rows an
//! index-calculus relation phase produces is worth a few per cent, not
//! a factor (538 against 582 `row_ops` on 63 two-summand rows over 121
//! columns, same rank, same answer).  What the dense matrix pays and
//! does not count is the scan of every column of every pivot row, and
//! `rank × |F|` words of storage against the sparse one's non-zeros;
//! that is the cost a structured elimination exists to avoid, and it
//! shows in `wall_ns` and `peak` before it shows in `row_ops`.
//!
//! ## The contract both obey
//!
//! `add_row` reduces the new row against the pivots held so far, and
//! either adopts it as a new pivot (`Independent`), finds it reduces to
//! `0 = 0` (`Dependent`) or to `0 = c ≠ 0` (`Inconsistent`).  `pinned`
//! answers whether the matrix determines one column on its own.  The
//! relation loop calls them in that order after every relation, so an
//! implementation must answer incrementally; a batch method would have
//! to re-solve on every call and is not what this trait is for.

use crate::cryptanalysis::ic_boundary::{IncrementalGauss, RelationSolver, RowStatus};
use crate::cryptanalysis::koblitz_strong_rho::RhoScalar;

#[cfg(test)]
fn addmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 + b as u128) % m as u128) as u64
}
fn submod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 + m as u128 - b as u128) % m as u128) as u64
}
fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}
fn invmod(a: u64, m: u64) -> Option<u64> {
    // Extended Euclid over i128; `m` is prime, so any non-zero `a` inverts.
    let (mut r0, mut r1) = (m as i128, (a % m) as i128);
    let (mut t0, mut t1) = (0i128, 1i128);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (t0, t1) = (t1, t0 - q * t1);
    }
    if r0 != 1 {
        return None;
    }
    Some(((t0 % m as i128 + m as i128) % m as i128) as u64)
}

/// A sparse row: `(column, coefficient)` sorted by column, no zeros.
type Sparse = Vec<(usize, u64)>;

/// `a − f·b`, both sparse and sorted, counting one `row_op` per
/// non-zero touched.
fn sparse_sub_scaled(a: &Sparse, f: u64, b: &Sparse, m: u64, ops: &mut u64) -> Sparse {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() || j < b.len() {
        match (a.get(i), b.get(j)) {
            (Some(&(ca, va)), Some(&(cb, vb))) if ca == cb => {
                let v = submod(va, mulmod(f, vb, m), m);
                *ops += 1;
                if v != 0 {
                    out.push((ca, v));
                }
                i += 1;
                j += 1;
            }
            (Some(&(ca, va)), Some(&(cb, _))) if ca < cb => {
                out.push((ca, va));
                i += 1;
            }
            (Some(&(ca, va)), None) => {
                out.push((ca, va));
                i += 1;
            }
            (_, Some(&(cb, vb))) => {
                let v = submod(0, mulmod(f, vb, m), m);
                *ops += 1;
                if v != 0 {
                    out.push((cb, v));
                }
                j += 1;
            }
            (None, None) => break,
        }
    }
    out
}

/// Sparse incremental elimination, pivoting on the lightest column.
///
/// Rows are kept sparse and the pivot of a new row is its non-zero
/// column of least global weight — the Markowitz-style choice that
/// keeps fill-in down.  Every multiply-subtract on a stored non-zero is
/// one `row_op`; nothing is charged for the zeros the dense method
/// walks over.  On relation rows with `m` non-zeros in `|F|` columns
/// that is the difference between `O(|F|)` and `O(m)` per reduction
/// step, which is where a real index-calculus linear algebra spends its
/// time and its memory.
pub struct StructuredGauss {
    modulus: u64,
    cols: usize,
    rows: Vec<Sparse>,
    rhs: Vec<u64>,
    pivot_of_row: Vec<usize>,
    row_of_pivot: Vec<Option<usize>>,
    /// Non-zeros per column across stored rows, for the pivot choice.
    weight: Vec<u32>,
    pub row_ops: u64,
    pub dependent: u64,
    pub inconsistent: u64,
    /// Total non-zeros stored: the method's footprint, which the dense
    /// method cannot report because it has none to speak of.
    pub nonzeros: u64,
    pub peak_nonzeros: u64,
}

impl StructuredGauss {
    pub fn new(cols: usize, modulus: u64) -> Self {
        Self {
            modulus,
            cols,
            rows: Vec::new(),
            rhs: Vec::new(),
            pivot_of_row: Vec::new(),
            row_of_pivot: vec![None; cols],
            weight: vec![0; cols],
            row_ops: 0,
            dependent: 0,
            inconsistent: 0,
            nonzeros: 0,
            peak_nonzeros: 0,
        }
    }

    fn note_weight(&mut self, row: &Sparse, delta: i32) {
        for &(c, _) in row {
            let w = &mut self.weight[c];
            *w = (*w as i32 + delta).max(0) as u32;
        }
    }
}

impl RelationSolver for StructuredGauss {
    fn add_row(&mut self, dense: Vec<u64>, mut rhs: u64) -> RowStatus {
        let m = self.modulus;
        debug_assert_eq!(dense.len(), self.cols);
        let mut row: Sparse = dense
            .iter()
            .enumerate()
            .filter(|(_, &v)| v != 0)
            .map(|(c, &v)| (c, v % m))
            .collect();
        // Reduce against every pivot the new row touches.  A sparse row
        // touches few, and the reduction adds only the pivot row's
        // non-zeros, so this is the step the structure pays for.
        loop {
            let hit = row
                .iter()
                .find_map(|&(c, f)| self.row_of_pivot[c].map(|i| (i, f)));
            let Some((i, f)) = hit else { break };
            let pivot_row = self.rows[i].clone();
            row = sparse_sub_scaled(&row, f, &pivot_row, m, &mut self.row_ops);
            rhs = submod(rhs, mulmod(f, self.rhs[i], m), m);
            self.row_ops += 1;
        }
        if row.is_empty() {
            if rhs == 0 {
                self.dependent += 1;
                return RowStatus::Dependent;
            }
            self.inconsistent += 1;
            return RowStatus::Inconsistent;
        }
        // Pivot on the lightest column the row touches.
        let (pc, pv) = row
            .iter()
            .copied()
            .min_by_key(|&(c, _)| (self.weight[c], c))
            .expect("non-empty row");
        let inv = invmod(pv, m).expect("prime modulus");
        for (_, v) in row.iter_mut() {
            *v = mulmod(*v, inv, m);
            self.row_ops += 1;
        }
        rhs = mulmod(rhs, inv, m);
        self.row_ops += 1;
        // Eliminate the new pivot column from every stored row that
        // holds it; only those rows are touched.
        for i in 0..self.rows.len() {
            let f = match self.rows[i].binary_search_by_key(&pc, |&(c, _)| c) {
                Ok(k) => self.rows[i][k].1,
                Err(_) => continue,
            };
            let before = self.rows[i].clone();
            let after = sparse_sub_scaled(&before, f, &row, m, &mut self.row_ops);
            self.rhs[i] = submod(self.rhs[i], mulmod(f, rhs, m), m);
            self.row_ops += 1;
            self.nonzeros = self.nonzeros - before.len() as u64 + after.len() as u64;
            self.note_weight(&before, -1);
            self.note_weight(&after, 1);
            self.rows[i] = after;
        }
        self.row_of_pivot[pc] = Some(self.rows.len());
        self.pivot_of_row.push(pc);
        self.nonzeros += row.len() as u64;
        self.note_weight(&row, 1);
        self.peak_nonzeros = self.peak_nonzeros.max(self.nonzeros);
        self.rows.push(row);
        self.rhs.push(rhs);
        RowStatus::Independent
    }

    fn pinned(&self, col: usize) -> Option<u64> {
        let i = self.row_of_pivot[col]?;
        let row = &self.rows[i];
        (row.len() == 1 && row[0].0 == col).then_some(self.rhs[i])
    }

    fn rank(&self) -> usize {
        self.rows.len()
    }

    fn dependent(&self) -> u64 {
        self.dependent
    }

    fn work(&self) -> (u64, &'static str) {
        (self.row_ops, "row_ops")
    }
}

/// Sparse batch solve by Wiedemann's method, for when the relation
/// matrix is too large to keep in echelon form.
///
/// Both eliminations above keep a reduced echelon form as rows arrive.
/// On index-calculus rows that form fills in: the dense matrix walks
/// `rank × |F|` entries per pivot, and the `k = 4` sessions in
/// characteristic 2 measured its `row_ops` growing as `|F|^{3.0}`
/// (`research/gaudry_k4_char2_crossover_20261015`).  This matrix keeps
/// the rows as they came and solves only when they can determine the
/// target:
///
/// 1. **Filter.** Columns no row touches are dropped. A column one row
///    touches is dropped with that row, repeatedly, since such a row
///    can only determine its own singleton. The target column is never
///    dropped.  What is left is the core: `k` rows on `n` columns.
/// 2. **Wait** until `k ≥ n`; before that, the core cannot have full
///    column rank and the solve is not attempted.
/// 3. **Solve** `A z = b` with `A = Rᵀ D R` (`D` a random non-zero
///    diagonal) and `b = Rᵀ D a`. `A` is symmetric, so the Krylov
///    sequence `s_i = bᵀ A^i b` takes `n` products for `2n` terms
///    (`s_{2j} = w_j·w_j`, `s_{2j+1} = w_j·w_{j+1}`, `w_j = A^j b`).
///    Berlekamp–Massey gives its minimal polynomial `f`; if `f(0) ≠ 0`,
///    `z = −f(0)^{-1} Σ_{j≥1} f_j A^{j−1} b`, `deg f − 1` more products.
/// 4. **Check** `R z = a` on every core row, exactly. A failed check or
///    `f(0) = 0` (a singular core) is a failed attempt; the next one
///    waits for `max(4, n/64)` more rows.
///
/// If `R` has full column rank, `A` is non-singular with probability
/// about `1 − n/r`, and then `z` is the unique solution, so the target
/// is pinned. `f(0) ≠ 0` with a singular `A` needs the random start to
/// miss the generalised kernel, probability about `1/r`; the relation
/// loop's group check (`[d]G = Q`) catches that case.
///
/// Work is counted in the dense matrices' unit: one `row_op` per
/// multiplication modulo `r` (products, dot products, Berlekamp–Massey,
/// the back-substitution and the check), failed attempts included.
/// Storing a row costs nothing. The method needs the core at full
/// column rank, where the echelon matrices pin the target as soon as it
/// lies in the row space, so it can take a few more relations; the
/// sessions report both.
pub struct Wiedemann<S: RhoScalar = u64> {
    modulus: S,
    cols: usize,
    rows: Vec<Vec<(usize, S)>>,
    rhs: Vec<S>,
    /// Non-zeros per column over all stored rows.
    weight: Vec<u32>,
    active: usize,
    state: std::cell::RefCell<WiedemannState<S>>,
}

struct WiedemannState<S> {
    row_ops: u64,
    attempts: u64,
    failed_attempts: u64,
    next_attempt_rows: usize,
    solution: Option<Vec<(usize, S)>>,
    core_rows: usize,
    core_cols: usize,
    matvecs: u64,
}

impl<S> Default for WiedemannState<S> {
    fn default() -> Self {
        Self {
            row_ops: 0,
            attempts: 0,
            failed_attempts: 0,
            next_attempt_rows: 0,
            solution: None,
            core_rows: 0,
            core_cols: 0,
            matvecs: 0,
        }
    }
}

impl<S: RhoScalar> Wiedemann<S> {
    pub fn new(cols: usize, modulus: S) -> Self {
        Self {
            modulus,
            cols,
            rows: Vec::new(),
            rhs: Vec::new(),
            weight: vec![0; cols],
            active: 0,
            state: std::cell::RefCell::new(WiedemannState::default()),
        }
    }

    /// `(attempts, failed attempts, core rows, core columns, products)`
    /// of the solve so far.
    pub fn stats(&self) -> (u64, u64, usize, usize, u64) {
        let s = self.state.borrow();
        (
            s.attempts,
            s.failed_attempts,
            s.core_rows,
            s.core_cols,
            s.matvecs,
        )
    }

    /// Store a sparse row (columns increasing, entries reduced and
    /// non-zero).
    pub fn add_sparse(&mut self, row: Vec<(usize, S)>, rhs: S) {
        for &(c, _) in &row {
            if self.weight[c] == 0 {
                self.active += 1;
            }
            self.weight[c] += 1;
        }
        self.rows.push(row);
        self.rhs.push(rhs);
    }

    /// The value of `col`, once the core determines it (see the type's
    /// documentation for when a solve is attempted).
    pub fn pinned_value(&self, col: usize) -> Option<S> {
        let mut st = self.state.borrow_mut();
        if let Some(sol) = &st.solution {
            return sol.iter().find(|&&(c, _)| c == col).map(|&(_, v)| v);
        }
        if self.rows.len() < self.active.max(st.next_attempt_rows) {
            return None;
        }
        let sol = self.attempt(col, &mut st)?;
        let v = sol.iter().find(|&&(c, _)| c == col).map(|&(_, v)| v);
        st.solution = Some(sol);
        v
    }

    pub fn rows_stored(&self) -> usize {
        self.rows.len()
    }

    pub fn solved_rank(&self) -> usize {
        let st = self.state.borrow();
        if st.solution.is_some() {
            st.core_cols
        } else {
            self.rows.len()
        }
    }

    pub fn row_ops(&self) -> u64 {
        self.state.borrow().row_ops
    }

    /// The core after singleton filtering: kept row indices and kept
    /// columns, in increasing order.
    fn core(&self, keep: usize) -> (Vec<usize>, Vec<usize>) {
        let mut row_alive = vec![true; self.rows.len()];
        let mut weight = self.weight.clone();
        let mut rows_of: Vec<Vec<usize>> = vec![Vec::new(); self.cols];
        for (i, row) in self.rows.iter().enumerate() {
            for &(c, _) in row {
                rows_of[c].push(i);
            }
        }
        let mut queue: Vec<usize> = (0..self.cols)
            .filter(|&c| c != keep && weight[c] == 1)
            .collect();
        while let Some(c) = queue.pop() {
            if weight[c] != 1 {
                continue;
            }
            let Some(&i) = rows_of[c].iter().find(|&&i| row_alive[i]) else {
                continue;
            };
            row_alive[i] = false;
            for &(c2, _) in &self.rows[i] {
                weight[c2] -= 1;
                if c2 != keep && weight[c2] == 1 {
                    queue.push(c2);
                }
            }
        }
        let rows = (0..self.rows.len()).filter(|&i| row_alive[i]).collect();
        let cols = (0..self.cols).filter(|&c| weight[c] > 0).collect();
        (rows, cols)
    }

    fn attempt(&self, target: usize, st: &mut WiedemannState<S>) -> Option<Vec<(usize, S)>> {
        let m = self.modulus;
        let ctx = S::mul_ctx(m);
        let mul = |a: S, b: S| S::mul_in(&ctx, a, b);
        let add = |a: S, b: S| S::add_mod(a, b, m);
        let (rows, cols) = self.core(target);
        st.core_rows = rows.len();
        st.core_cols = cols.len();
        let n = cols.len();
        if rows.len() < n || !cols.contains(&target) {
            st.next_attempt_rows = self.rows.len() + (n - rows.len().min(n)).max(1);
            return None;
        }
        st.attempts += 1;
        let mut ops = 0u64;
        let mut index = vec![usize::MAX; self.cols];
        for (j, &c) in cols.iter().enumerate() {
            index[c] = j;
        }
        let r: Vec<Vec<(usize, S)>> = rows
            .iter()
            .map(|&i| self.rows[i].iter().map(|&(c, v)| (index[c], v)).collect())
            .collect();
        let a: Vec<S> = rows.iter().map(|&i| self.rhs[i]).collect();
        // Deterministic randomness: replays must reproduce every count.
        // One 64-bit draw per entry, reduced below the modulus and made
        // non-zero; a modulus past 2^64 takes two words.
        let mut seed = 0x9e37_79b9_7f4a_7c15u64 ^ st.attempts.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let m128 = m.to_u128();
        let d: Vec<S> = (0..r.len())
            .map(|_| {
                let w = if m128 <= u64::MAX as u128 {
                    u128::from(next())
                } else {
                    u128::from(next()) | (u128::from(next()) << 64)
                };
                S::from_biguint(&num_bigint::BigUint::from(1 + w % (m128 - 1))).expect("below the modulus")
            })
            .collect();
        // y = Rᵀ D (R x)
        let apply = |x: &[S], ops: &mut u64| -> Vec<S> {
            let mut y = vec![S::zero(); n];
            for (row, &di) in r.iter().zip(&d) {
                let mut t = S::zero();
                for &(j, v) in row {
                    t = add(t, mul(v, x[j]));
                }
                *ops += row.len() as u64;
                if t == S::zero() {
                    continue;
                }
                let t = mul(t, di);
                *ops += 1;
                for &(j, v) in row {
                    y[j] = add(y[j], mul(v, t));
                }
                *ops += row.len() as u64;
            }
            y
        };
        let dot = |x: &[S], y: &[S], ops: &mut u64| -> S {
            *ops += n as u64;
            x.iter().zip(y).fold(S::zero(), |acc, (&p, &q)| add(acc, mul(p, q)))
        };
        let mut b = vec![S::zero(); n];
        for ((row, &di), &ai) in r.iter().zip(&d).zip(&a) {
            let t = mul(di, ai);
            for &(j, v) in row {
                b[j] = add(b[j], mul(v, t));
            }
            ops += 1 + row.len() as u64;
        }
        // s_0 .. s_{2n-1} from w_0 .. w_n.
        let mut seq = Vec::with_capacity(2 * n);
        let mut w = b.clone();
        for _ in 0..n {
            let next = apply(&w, &mut ops);
            st.matvecs += 1;
            seq.push(dot(&w, &w, &mut ops));
            seq.push(dot(&w, &next, &mut ops));
            w = next;
        }
        let c = berlekamp_massey(&seq, m, &mut ops);
        let l = c.len() - 1;
        let solution = (|| {
            if l == 0 || c[l] == S::zero() {
                return None;
            }
            // f(λ) = Σ_j f_j λ^j with f_j = c_{l−j}; f_0 = c_l.
            let f0_inv = S::inverse_mod(c[l], m)?;
            let mut z = vec![S::zero(); n];
            let mut w = b.clone();
            for j in 1..=l {
                let fj = c[l - j];
                if fj != S::zero() {
                    for (zi, &wi) in z.iter_mut().zip(&w) {
                        *zi = add(*zi, mul(fj, wi));
                    }
                    ops += n as u64;
                }
                if j < l {
                    w = apply(&w, &mut ops);
                    st.matvecs += 1;
                }
            }
            let scale = S::sub_mod(S::zero(), f0_inv, m);
            for zi in z.iter_mut() {
                *zi = mul(*zi, scale);
            }
            ops += n as u64;
            // The exact check, on every core row.
            for (row, &ai) in r.iter().zip(&a) {
                let mut t = S::zero();
                for &(j, v) in row {
                    t = add(t, mul(v, z[j]));
                }
                ops += row.len() as u64;
                if t != ai {
                    return None;
                }
            }
            Some(cols.iter().zip(&z).map(|(&c, &v)| (c, v)).collect())
        })();
        st.row_ops += ops;
        if solution.is_none() {
            st.failed_attempts += 1;
            st.next_attempt_rows = self.rows.len() + (n / 64).max(4);
        }
        solution
    }
}

/// The shortest linear recurrence of `s` modulo `m`, as its connection
/// polynomial `C` (`C[0] = 1`, `Σ_i C[i] s[k−i] = 0`), counting one op
/// per multiplication.
fn berlekamp_massey<S: RhoScalar>(s: &[S], m: S, ops: &mut u64) -> Vec<S> {
    let ctx = S::mul_ctx(m);
    let mul = |a: S, b: S| S::mul_in(&ctx, a, b);
    let mut c = vec![S::one()];
    let mut bpoly = vec![S::one()];
    let mut l = 0usize;
    let mut shift = 1usize;
    let mut bd = S::one();
    for k in 0..s.len() {
        let mut d = s[k];
        for i in 1..=l.min(c.len() - 1) {
            d = S::add_mod(d, mul(c[i], s[k - i]), m);
        }
        *ops += l as u64;
        if d == S::zero() {
            shift += 1;
            continue;
        }
        let coef = mul(d, S::inverse_mod(bd, m).expect("prime modulus"));
        *ops += 2;
        let old = c.clone();
        if c.len() < bpoly.len() + shift {
            c.resize(bpoly.len() + shift, S::zero());
        }
        for (i, &bv) in bpoly.iter().enumerate() {
            c[i + shift] = S::sub_mod(c[i + shift], mul(coef, bv), m);
        }
        *ops += bpoly.len() as u64;
        if 2 * l <= k {
            l = k + 1 - l;
            bpoly = old;
            bd = d;
            shift = 1;
        } else {
            shift += 1;
        }
    }
    c.truncate(l + 1);
    c.resize(l + 1, S::zero());
    c
}

impl RelationSolver for Wiedemann<u64> {
    fn add_row(&mut self, dense: Vec<u64>, rhs: u64) -> RowStatus {
        let m = self.modulus;
        debug_assert_eq!(dense.len(), self.cols);
        let row: Vec<(usize, u64)> = dense
            .iter()
            .enumerate()
            .filter(|(_, &v)| v % m != 0)
            .map(|(c, &v)| (c, v % m))
            .collect();
        self.add_sparse(row, rhs % m);
        RowStatus::Independent
    }

    fn pinned(&self, col: usize) -> Option<u64> {
        self.pinned_value(col)
    }

    fn rank(&self) -> usize {
        // The core's dimension once solved; before that the method has
        // computed no rank, and the rows stored are reported.
        self.solved_rank()
    }

    fn dependent(&self) -> u64 {
        0
    }

    fn work(&self) -> (u64, &'static str) {
        (self.row_ops(), "row_ops")
    }
}

/// The matrices the framework knows, by name.
pub const MATRIX_NAMES: &[(&str, &str)] = &[
    (
        "incremental-gauss",
        "dense reduced row echelon maintained as rows arrive, pivot on the leftmost non-zero; one row_op per non-zero multiply-subtract",
    ),
    (
        "incremental-gauss-full-rank",
        "same dense incremental matrix, but collect until every factor-base column and the target column are independent",
    ),
    (
        "incremental-gauss-full-rank-checked",
        "same full-rank matrix, then charge a full-point logarithm check for one representative of every factor-base column",
    ),
    (
        "structured-gauss",
        "sparse reduced row echelon maintained as rows arrive, pivot on the lightest column; one row_op per non-zero multiply-subtract, so the two are comparable",
    ),
    (
        "wiedemann",
        "rows kept sparse, singleton-filtered core solved by Wiedemann (A = R^T D R, Berlekamp-Massey) once it has as many rows as columns, checked exactly; one row_op per multiplication mod r",
    ),
];

/// Build one by name, for a matrix with `cols` columns over `Z/rZ`.
pub fn matrix_by_name(
    name: &str,
    cols: usize,
    modulus: u64,
) -> Result<Box<dyn RelationSolver>, String> {
    match name {
        "incremental-gauss"
        | "incremental-gauss-full-rank"
        | "incremental-gauss-full-rank-checked" => {
            Ok(Box::new(IncrementalGauss::new(cols, modulus)))
        }
        "structured-gauss" => Ok(Box::new(StructuredGauss::new(cols, modulus))),
        "wiedemann" => Ok(Box::new(Wiedemann::<u64>::new(cols, modulus))),
        other => Err(format!(
            "unknown relation matrix `{other}`; known: {}",
            MATRIX_NAMES
                .iter()
                .map(|(n, _)| *n)
                .collect::<Vec<_>>()
                .join(", ")
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    /// **Both matrices must pin the same unknown to the same value on
    /// the same rows, at the same moment.**  A second implementation
    /// that agreed on the answer but pinned it a row later would change
    /// every trial count downstream while looking correct.
    #[test]
    fn structured_and_dense_agree_row_by_row() {
        let modulus = 1_000_003u64;
        let cols = 12usize;
        let target = cols - 1;
        let mut rng = StdRng::seed_from_u64(20260922);
        // A planted solution, so every row is consistent and the value
        // pinned can be checked against something.
        let x: Vec<u64> = (0..cols).map(|_| rng.gen_range(0..modulus)).collect();
        let mut dense = IncrementalGauss::new(cols, modulus);
        let mut sparse = StructuredGauss::new(cols, modulus);
        let mut pinned_at = None;
        for row_no in 0..40 {
            // Sparse relation rows: two or three non-zeros plus the target.
            let mut row = vec![0u64; cols];
            for _ in 0..rng.gen_range(2..=3) {
                row[rng.gen_range(0..cols)] = rng.gen_range(1..modulus);
            }
            row[target] = rng.gen_range(1..modulus);
            let rhs = row.iter().zip(&x).fold(0u64, |acc, (&a, &b)| {
                addmod(acc, mulmod(a, b, modulus), modulus)
            });
            let a = dense.add_row(row.clone(), rhs);
            let b = sparse.add_row(row, rhs);
            assert_eq!(
                std::mem::discriminant(&a),
                std::mem::discriminant(&b),
                "row {row_no}: the two matrices classified it differently"
            );
            assert_eq!(
                dense.rank(),
                RelationSolver::rank(&sparse),
                "row {row_no}: rank"
            );
            let pd = dense.pinned(target);
            let ps = sparse.pinned(target);
            assert_eq!(pd, ps, "row {row_no}: pinned value");
            if let Some(v) = ps {
                assert_eq!(v, x[target], "pinned the wrong value");
                pinned_at.get_or_insert(row_no);
            }
        }
        assert!(
            pinned_at.is_some(),
            "forty rows on twelve columns must pin the target"
        );
        // And the structure must have bought something: fewer
        // multiply-subtracts than the dense walk over every column.
        assert!(
            sparse.row_ops < dense.row_ops,
            "structured {} vs dense {} row_ops",
            sparse.row_ops,
            dense.row_ops
        );
    }

    /// An inconsistent row must be reported, not absorbed.
    #[test]
    fn an_inconsistent_row_is_reported() {
        let modulus = 101u64;
        let mut m = StructuredGauss::new(3, modulus);
        assert!(matches!(
            m.add_row(vec![1, 0, 0], 5),
            RowStatus::Independent
        ));
        assert!(matches!(m.add_row(vec![1, 0, 0], 5), RowStatus::Dependent));
        assert!(matches!(
            m.add_row(vec![2, 0, 0], 11),
            RowStatus::Inconsistent
        ));
        assert_eq!(m.pinned(0), Some(5));
    }

    #[test]
    fn the_registry_resolves_every_matrix_it_advertises() {
        for (name, _) in MATRIX_NAMES {
            assert!(matrix_by_name(name, 4, 101).is_ok(), "{name}");
        }
        assert!(matrix_by_name("lanczos", 4, 101).is_err());
    }

    /// Wiedemann pins the planted value once the core is square, and
    /// never pins a wrong one.
    #[test]
    fn wiedemann_pins_the_planted_value() {
        let modulus = 1_000_000_007u64;
        for seed in 0..20u64 {
            let cols = 60usize;
            let target = cols - 1;
            let mut rng = StdRng::seed_from_u64(20261015 + seed);
            let x: Vec<u64> = (0..cols).map(|_| rng.gen_range(0..modulus)).collect();
            let mut w = Wiedemann::<u64>::new(cols, modulus);
            let mut dense = IncrementalGauss::new(cols, modulus);
            let mut pinned = None;
            let mut dense_pin = None;
            for row_no in 0..200 {
                let mut row = vec![0u64; cols];
                for _ in 0..4 {
                    row[rng.gen_range(0..cols - 1)] = rng.gen_range(1..modulus);
                }
                row[target] = rng.gen_range(1..modulus);
                let rhs = row.iter().zip(&x).fold(0u64, |acc, (&a, &b)| {
                    addmod(acc, mulmod(a, b, modulus), modulus)
                });
                dense.add_row(row.clone(), rhs);
                w.add_row(row, rhs);
                if dense_pin.is_none() && dense.pinned(target).is_some() {
                    dense_pin = Some(row_no);
                }
                if let Some(v) = w.pinned(target) {
                    assert_eq!(v, x[target], "seed {seed}: wrong value");
                    pinned = Some(row_no);
                    break;
                }
            }
            let (p, d) = (pinned.expect("pins"), dense_pin.expect("dense pins"));
            assert!(p >= d, "seed {seed}: Wiedemann cannot pin before the row space allows");
            assert!(p <= d + 16, "seed {seed}: pinned at {p}, dense at {d}");
        }
    }

    #[test]
    fn berlekamp_massey_finds_a_known_recurrence() {
        // s_k = 3 s_{k-1} + 5 s_{k-2} mod 101
        let m = 101u64;
        let mut s = vec![1u64, 2];
        for k in 2..12 {
            s.push((3 * s[k - 1] + 5 * s[k - 2]) % m);
        }
        let mut ops = 0;
        let c = berlekamp_massey(&s, m, &mut ops);
        assert_eq!(c, vec![1, m - 3, m - 5]);
    }
}
