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
//! | `wiedemann` | Wiedemann on the transpose for the target column only, attempted once rows outnumber touched columns ([`Wiedemann`]) | the non-zeros | none | `row_ops` (every multiplication mod `r`) |
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

/// Wiedemann on the relation matrix's transpose, solved for the target
/// column only.
///
/// The rows arrive as `Σ c_j·x_j + t·d = rhs` with `d` the target's
/// logarithm in the **last** column.  Write `C` for the `m × k` block of
/// the factor-base columns some row touches.  Once `m > k` the left
/// kernel of `C` is non-zero, and any `λ` in it gives
/// `d·(λ·t) = λ·rhs`, so `d = (λ·rhs)/(λ·t)` whenever `λ·t ≠ 0`.  `λ` is
/// found by Wiedemann's method on the square `m × m` matrix `B = [Cᵀ; 0]`:
/// the sequence `uᵀBⁱz` for `i < 2m`, Berlekamp–Massey for its minimal
/// polynomial `Xᵉ·g(X)`, `v = g(B)z`, and the last non-zero `Bʲv`.
/// Every product `B·w` costs the stored non-zeros, so a solve costs about
/// `3m` products, `2m` dot products and the `O(m²)` of Berlekamp–Massey:
/// `O(m·nnz + m²)`, against the `O(m²·k)` of an elimination.  This is the
/// sparse linear algebra the `q²` term of the index-calculus analysis
/// assumes.
///
/// **Counting.**  Every multiplication modulo `r` — in a product, a dot
/// product, Berlekamp–Massey, the Horner evaluation, the final division —
/// is one `row_op`, the unit the eliminations count their multiply-
/// subtracts in, and nothing else is.  Additions are not counted, as in
/// the eliminations.  Failed attempts are counted like successful ones.
///
/// **When it solves.**  An attempt runs when a row arrives and `m ≥ k + 1`,
/// and after an attempt that pins nothing the next waits for
/// `max(1, m/16)` more rows, so the attempts' total is a constant factor
/// on the last one.  An attempt that finds `λ·t = 0` or whose
/// polynomial was not the minimal one retries once with fresh vectors.
/// The randomness is a fixed function of the matrix's shape, so a replay
/// reproduces every count.
///
/// **Shadow elimination.**  The same rows also go, until it pins the
/// target, to an [`IncrementalGauss`] that is reported (`shadow_gauss_*`)
/// and not priced: the dense elimination's count on the same relation
/// stream, for comparison.  `rank` is `min(m, k)`, an upper bound: the
/// method does not compute the rank.
pub struct Wiedemann {
    modulus: u64,
    cols: usize,
    /// Factor-base part of each row (target column excluded).
    rows: Vec<Sparse>,
    tcoef: Vec<u64>,
    rhs: Vec<u64>,
    /// Compressed index of each touched factor-base column.
    index_of: Vec<Option<usize>>,
    touched: usize,
    next_attempt: usize,
    pinned_value: Option<u64>,
    attempt_no: u64,
    pub row_ops: u64,
    pub attempts: u64,
    pub retries: u64,
    pub matvecs: u64,
    pub last_dimension: u64,
    pub last_row_ops: u64,
    shadow: Option<IncrementalGauss>,
    shadow_pin: Option<(u64, u64, u64)>,
    agrees: Option<bool>,
}

impl Wiedemann {
    pub fn new(cols: usize, modulus: u64) -> Self {
        assert!(cols >= 2, "a target column and at least one base column");
        Self {
            modulus,
            cols,
            rows: Vec::new(),
            tcoef: Vec::new(),
            rhs: Vec::new(),
            index_of: vec![None; cols - 1],
            touched: 0,
            next_attempt: 0,
            pinned_value: None,
            attempt_no: 0,
            row_ops: 0,
            attempts: 0,
            retries: 0,
            matvecs: 0,
            last_dimension: 0,
            last_row_ops: 0,
            shadow: Some(IncrementalGauss::new(cols, modulus)),
            shadow_pin: None,
            agrees: None,
        }
    }

    fn mul(&mut self, a: u64, b: u64) -> u64 {
        self.row_ops += 1;
        mulmod(a, b, self.modulus)
    }

    /// `out = B·w`, `B = [Cᵀ; 0]` (`m × m`): `out[idx(c)] += C[i][c]·w[i]`.
    fn apply(&mut self, w: &[u64]) -> Vec<u64> {
        let m = self.modulus;
        let mut out = vec![0u64; w.len()];
        let mut ops = 0u64;
        for (i, row) in self.rows.iter().enumerate() {
            let wi = w[i];
            if wi == 0 {
                continue;
            }
            for &(c, v) in row {
                let j = self.index_of[c].expect("touched column");
                out[j] = addmod_lin(out[j], mulmod(v, wi, m), m);
                ops += 1;
            }
        }
        self.row_ops += ops;
        self.matvecs += 1;
        out
    }

    fn dot(&mut self, a: &[u64], b: &[u64]) -> u64 {
        let m = self.modulus;
        let mut s = 0u64;
        let mut ops = 0u64;
        for (&x, &y) in a.iter().zip(b) {
            if x != 0 && y != 0 {
                s = addmod_lin(s, mulmod(x, y, m), m);
                ops += 1;
            }
        }
        self.row_ops += ops;
        s
    }

    /// Berlekamp–Massey over `F_r`: the connection polynomial
    /// `1 + c₁X + … + c_LX^L` of `s`, as `(coefficients, L)`.
    fn berlekamp_massey(&mut self, s: &[u64]) -> (Vec<u64>, usize) {
        let m = self.modulus;
        let mut c = vec![1u64];
        let mut b = vec![1u64];
        let (mut l, mut shift, mut bd) = (0usize, 1usize, 1u64);
        for n in 0..s.len() {
            let mut d = s[n];
            for i in 1..=l.min(c.len() - 1) {
                let t = self.mul(c[i], s[n - i]);
                d = addmod_lin(d, t, m);
            }
            if d == 0 {
                shift += 1;
                continue;
            }
            let coef = self.mul(d, invmod(bd, m).expect("prime modulus"));
            let old = c.clone();
            if c.len() < b.len() + shift {
                c.resize(b.len() + shift, 0);
            }
            for (i, &bi) in b.iter().enumerate() {
                if bi != 0 {
                    let t = self.mul(coef, bi);
                    c[i + shift] = submod(c[i + shift], t, m);
                }
            }
            if 2 * l <= n {
                l = n + 1 - l;
                b = old;
                bd = d;
                shift = 1;
            } else {
                shift += 1;
            }
        }
        c.resize(l + 1, 0);
        (c, l)
    }

    /// One attempt with the `try_no`-th pair of random vectors: the
    /// target's value, or `None`.
    fn try_solve(&mut self, try_no: u64) -> Option<u64> {
        use rand::{Rng, SeedableRng};
        let m = self.modulus;
        let n = self.rows.len();
        let mut rng = rand::rngs::StdRng::seed_from_u64(
            0x5749_4544 ^ (self.cols as u64) << 20 ^ self.attempt_no << 8 ^ try_no ^ m,
        );
        let u: Vec<u64> = (0..n).map(|_| rng.gen_range(1..m)).collect();
        let z: Vec<u64> = (0..n).map(|_| rng.gen_range(0..m)).collect();
        // s_i = uᵀ Bⁱ z, i < 2n.
        let mut seq = Vec::with_capacity(2 * n);
        let mut w = z.clone();
        for i in 0..2 * n {
            seq.push(self.dot(&u, &w));
            if i + 1 < 2 * n {
                w = self.apply(&w);
            }
        }
        let (conn, l) = self.berlekamp_massey(&seq);
        // f(X) = X^L·conn(1/X): f_k = conn[L−k].  f = Xᵉ·g with g(0) ≠ 0,
        // e the number of trailing zeros of conn read from the top.
        let f: Vec<u64> = (0..=l).map(|k| conn[l - k]).collect();
        let e = f.iter().take_while(|&&x| x == 0).count();
        if e == 0 || e > l {
            return None;
        }
        let g = &f[e..];
        // v = g(B)·z by Horner from the leading coefficient.
        let mut v = vec![0u64; n];
        for k in (0..g.len()).rev() {
            if k + 1 < g.len() {
                v = self.apply(&v);
            }
            let gk = g[k];
            if gk != 0 {
                for i in 0..n {
                    if z[i] != 0 {
                        let t = self.mul(gk, z[i]);
                        v[i] = addmod_lin(v[i], t, m);
                    }
                }
            }
        }
        if v.iter().all(|&x| x == 0) {
            return None;
        }
        let mut lambda = v;
        let mut found = false;
        for _ in 0..=e {
            let next = self.apply(&lambda);
            if next.iter().all(|&x| x == 0) {
                found = true;
                break;
            }
            lambda = next;
        }
        if !found {
            return None;
        }
        let tc = self.tcoef.clone();
        let rh = self.rhs.clone();
        let den = self.dot(&lambda, &tc);
        if den == 0 {
            return None;
        }
        let num = self.dot(&lambda, &rh);
        Some(self.mul(num, invmod(den, m).expect("prime modulus")))
    }

    fn attempt(&mut self) {
        let before = self.row_ops;
        self.attempts += 1;
        self.last_dimension = self.rows.len() as u64;
        for try_no in 0..2 {
            if try_no > 0 {
                self.retries += 1;
            }
            if let Some(d) = self.try_solve(try_no) {
                self.pinned_value = Some(d);
                break;
            }
        }
        self.attempt_no += 1;
        self.last_row_ops = self.row_ops - before;
    }
}

fn addmod_lin(a: u64, b: u64, m: u64) -> u64 {
    let s = a as u128 + b as u128;
    (if s >= m as u128 { s - m as u128 } else { s }) as u64
}

impl RelationSolver for Wiedemann {
    fn add_row(&mut self, dense: Vec<u64>, rhs: u64) -> RowStatus {
        let m = self.modulus;
        debug_assert_eq!(dense.len(), self.cols);
        let last = self.cols - 1;
        if let Some(shadow) = self.shadow.as_mut() {
            shadow.add_row(dense.clone(), rhs);
            if let Some(v) = shadow.pinned(last) {
                self.shadow_pin = Some((shadow.rank() as u64, shadow.row_ops, v));
                self.shadow = None;
            }
        }
        if self.pinned_value.is_some() {
            return RowStatus::Independent;
        }
        let row: Sparse = dense[..last]
            .iter()
            .enumerate()
            .filter(|(_, &v)| v % m != 0)
            .map(|(c, &v)| (c, v % m))
            .collect();
        for &(c, _) in &row {
            if self.index_of[c].is_none() {
                self.index_of[c] = Some(self.touched);
                self.touched += 1;
            }
        }
        self.rows.push(row);
        self.tcoef.push(dense[last] % m);
        self.rhs.push(rhs % m);
        let n = self.rows.len();
        if n > self.touched && n >= self.next_attempt {
            self.attempt();
            if self.pinned_value.is_none() {
                self.next_attempt = n + (n / 16).max(1);
            }
        }
        if let (Some(d), Some((_, _, v))) = (self.pinned_value, self.shadow_pin) {
            self.agrees.get_or_insert(d == v);
        }
        RowStatus::Independent
    }

    fn pinned(&self, col: usize) -> Option<u64> {
        (col == self.cols - 1)
            .then_some(self.pinned_value)
            .flatten()
    }

    fn rank(&self) -> usize {
        self.rows.len().min(self.touched)
    }

    fn dependent(&self) -> u64 {
        0
    }

    fn work(&self) -> (u64, &'static str) {
        (self.row_ops, "row_ops")
    }

    fn native(&self) -> Vec<(&'static str, u64)> {
        let mut out = vec![
            ("wiedemann_attempts", self.attempts),
            ("wiedemann_retries", self.retries),
            ("wiedemann_matvecs", self.matvecs),
            ("wiedemann_dimension", self.last_dimension),
            ("wiedemann_touched_columns", self.touched as u64),
            ("wiedemann_last_attempt_row_ops", self.last_row_ops),
            (
                "wiedemann_nonzeros",
                self.rows.iter().map(|r| r.len() as u64).sum(),
            ),
        ];
        if let Some((rank, ops, _)) = self.shadow_pin {
            out.push(("shadow_gauss_pin_rank", rank));
            out.push(("shadow_gauss_row_ops", ops));
        }
        if let Some(a) = self.agrees {
            out.push(("shadow_gauss_agrees", u64::from(a)));
        }
        out
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
        "Wiedemann on the transpose for the target column alone, attempted once the rows outnumber the touched columns; one row_op per multiplication mod r (products, dot products, Berlekamp–Massey), O(m·nnz + m²); the dense elimination on the same rows is reported as shadow_gauss_* and not priced",
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
        "wiedemann" => Ok(Box::new(Wiedemann::new(cols, modulus))),
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
        assert!(matrix_by_name("block-lanczos", 4, 101).is_err());
    }

    /// Wiedemann pins the planted target value on sparse planted rows,
    /// the same value the dense elimination pins, and replays its count.
    #[test]
    fn wiedemann_pins_the_planted_target() {
        let modulus = 1_000_003u64;
        for (cols, seed) in [(12usize, 1u64), (40, 2), (90, 3)] {
            let mut rng = StdRng::seed_from_u64(seed);
            let x: Vec<u64> = (0..cols).map(|_| rng.gen_range(0..modulus)).collect();
            let run = || {
                let mut rng = StdRng::seed_from_u64(seed ^ 0xabc);
                let mut w = Wiedemann::new(cols, modulus);
                let mut rows = 0;
                while w.pinned(cols - 1).is_none() {
                    let mut row = vec![0u64; cols];
                    for _ in 0..4 {
                        let c = rng.gen_range(0..cols - 1);
                        row[c] =
                            (row[c] + if rng.gen_bool(0.5) { 1 } else { modulus - 1 }) % modulus;
                    }
                    row[cols - 1] = rng.gen_range(1..modulus);
                    let rhs = row.iter().zip(&x).fold(0u64, |s, (&a, &b)| {
                        addmod(s, mulmod(a, b, modulus), modulus)
                    });
                    w.add_row(row, rhs);
                    rows += 1;
                    assert!(rows < 10 * cols, "never pinned");
                }
                w
            };
            let w = run();
            assert_eq!(w.pinned(cols - 1), Some(x[cols - 1]));
            let native: std::collections::BTreeMap<_, _> = w.native().into_iter().collect();
            assert_eq!(native.get("shadow_gauss_agrees"), Some(&1));
            assert_eq!(run().row_ops, w.row_ops, "deterministic count");
            assert!(w.row_ops > 0);
        }
    }
}
