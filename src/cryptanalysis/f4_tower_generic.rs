//! # F4 over a tower quotient ring, generic over the prime field
//!
//! The algorithm of [`f4_fp_tower`](super::f4_fp_tower) — the same
//! monomials, pairs, criteria, symbolic preprocessing, Faugère–Lachartre
//! elimination, full-rank exit and trace — over any
//! [`PrimeField`](super::fp_backend::PrimeField) backend, so that a PKM
//! tower system can be solved at a prime above `2³²`: Kummer towers at
//! `c·2^k + 1` with `k ≥ 40` ([`Fp64`](super::fp_backend::Fp64)), at NIST
//! P-224 ([`Fp224`](super::fp_backend::Fp224)), and Dickson torus towers
//! at P-192 / P-256 ([`FpWide`](super::fp_backend::FpWide)).
//!
//! The control flow is the `u32` engine's, step for step: the same pairs
//! are selected in the same order, the same rows are formed, and the
//! eliminations visit rows and pivots in the same order. Only the
//! arithmetic changes:
//!
//! - coefficients are `F::Elem` (`u32` for [`Fp32`](super::fp_backend::Fp32),
//!   Montgomery words above);
//! - the dense kernels (`B' = A⁻¹B`, `D − C·B'`, the residues' echelon) run
//!   on `F::Wide` accumulators through `row_load` / `row_sub` /
//!   `row_reduce` / `row_reduce_scaled`, one lazy product per pivot row
//!   and one reduction per entry when the row is read;
//! - the sparse accumulator of the input echelon subtracts a pivot row with
//!   `row_scatter_axpy` and sets its bitset in a separate pass.
//!
//! Because every value the control flow reads (a coefficient's being zero,
//! a pivot's position, a row's length) is a field fact, and the
//! elimination's block and chunk sizes are the `u32` engine's whatever the
//! accumulator's width, the engine at one
//! prime gives the same report — solving degree, widths, steps,
//! multiply-adds, basis — under every backend that can hold that prime,
//! and at `Fp32` the same report as `f4_fp_tower`. The tests check both.
//!
//! **Arithmetic cost.** At `Fp32` with `p < 2³¹` the dense kernel is the
//! `u32` engine's own lazy loop (AVX2 / AVX-512 dispatch). At
//! `2³¹ ≤ p < 2³²` it is Shoup's product with a carry fold, where the `u32`
//! engine reduces every step. At `Fp64` it is a `u128` accumulator kept
//! below `p·2⁶⁴`; at `FpWide<N>` a `2N + 1`-word sum reduced once.
//!
//! **Lazy capacity.** Every accumulator row counts the products added to it
//! since it was loaded, and is reduced and reloaded in place when the count
//! reaches the backend's
//! [`lazy_capacity`](super::fp_backend::PrimeField::lazy_capacity). For
//! `Fp32`, `Fp64` and `FpWide` (capacity at least `2⁶⁴ − 2`) that never
//! happens; the IFMA backend [`Fp52`](super::fp_backend_ifma::Fp52) has a
//! capacity of 4095 products. A reload changes no value, so the report is
//! the same either way (tested with a capacity of one).
//!
//! Limits as in `f4_fp_tower`: at most 64 tower variables and 4 free
//! variables with exponents below 256.

use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use std::sync::atomic::{AtomicBool, Ordering as AtomicOrdering};
use std::sync::Arc;
use std::time::Instant;

use rayon::prelude::*;

use super::f4_fp_tower::{
    add_free, any_between, coprime, degree_of, divides, free_of, lcm, mask_of, mono, quotient,
    staircase_at_most, Mono, MonoHasher, MonoMap, RPoly, SquareRule, StepTrace, TowerF4Options,
    TowerF4Report, TowerRing, MAX_FREE, NONE, ONE,
};
use super::fp_backend::PrimeField;

// ── The ring ─────────────────────────────────────────────────────────

/// How one tower variable's square rewrites:
/// `y_j² = a·y_j·y_next + b·y_next + c·y_j + d`. On a tower's top level
/// `next` is `None` and `a = b = 0`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GSquareRule<E> {
    pub next: Option<u8>,
    pub a: E,
    pub b: E,
    pub c: E,
    pub d: E,
}

/// `F_p[y_0, …, y_{n−1}, u_0, …] / (y_j² − rule_j)` over the backend `F`.
#[derive(Clone, Debug)]
pub struct GTowerRing<F: PrimeField> {
    pub field: F,
    pub rules: Vec<GSquareRule<F::Elem>>,
    pub n_free: usize,
}

/// An element of `R` in normal form: monomials strictly descending,
/// coefficients non-zero.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GPoly<E> {
    pub monos: Vec<Mono>,
    pub coefs: Vec<E>,
}

impl<E> Default for GPoly<E> {
    fn default() -> Self {
        GPoly {
            monos: Vec::new(),
            coefs: Vec::new(),
        }
    }
}

impl<E: Copy> GPoly<E> {
    pub fn len(&self) -> usize {
        self.monos.len()
    }
    pub fn is_empty(&self) -> bool {
        self.monos.is_empty()
    }
    pub fn lm(&self) -> Option<Mono> {
        self.monos.first().copied()
    }
    pub fn degree(&self) -> u32 {
        self.monos.iter().map(|&m| degree_of(m)).max().unwrap_or(0)
    }
    pub fn is_unit(&self) -> bool {
        self.monos.len() == 1 && self.monos[0] == ONE
    }

    /// Sort, merge repeats and drop zeros.
    pub fn from_terms<F: PrimeField<Elem = E>>(f: &F, mut terms: Vec<(Mono, E)>) -> Self {
        terms.sort_unstable_by_key(|t| std::cmp::Reverse(t.0));
        let mut monos = Vec::with_capacity(terms.len());
        let mut coefs = Vec::with_capacity(terms.len());
        let mut i = 0;
        while i < terms.len() {
            let m = terms[i].0;
            let mut c = f.zero();
            while i < terms.len() && terms[i].0 == m {
                c = f.add(c, terms[i].1);
                i += 1;
            }
            if !f.is_zero(c) {
                monos.push(m);
                coefs.push(c);
            }
        }
        GPoly { monos, coefs }
    }

    /// The polynomial of a `u32`-engine [`RPoly`], coefficient by
    /// coefficient (`x mod p`).
    pub fn from_rpoly<F: PrimeField<Elem = E>>(f: &F, g: &RPoly) -> Self {
        GPoly {
            monos: g.monos.clone(),
            coefs: g.coefs.iter().map(|&c| f.from_u64(u64::from(c))).collect(),
        }
    }

    /// The `u32`-engine polynomial, when every coefficient is below `2³²`.
    pub fn to_rpoly<F: PrimeField<Elem = E>>(&self, f: &F) -> Option<RPoly> {
        let coefs = self
            .coefs
            .iter()
            .map(|&c| f.to_u64(c).and_then(|x| u32::try_from(x).ok()))
            .collect::<Option<Vec<u32>>>()?;
        Some(RPoly {
            monos: self.monos.clone(),
            coefs,
        })
    }
}

/// Combine repeated masks and drop zero coefficients, in place.
fn merge_masks<F: PrimeField>(f: &F, terms: &mut Vec<(u64, F::Elem)>) {
    if terms.len() < 2 {
        return;
    }
    terms.sort_unstable_by_key(|t| t.0);
    let mut out = 0usize;
    let mut i = 0usize;
    while i < terms.len() {
        let m = terms[i].0;
        let mut c = f.zero();
        while i < terms.len() && terms[i].0 == m {
            c = f.add(c, terms[i].1);
            i += 1;
        }
        if !f.is_zero(c) {
            terms[out] = (m, c);
            out += 1;
        }
    }
    terms.truncate(out);
}

impl<F: PrimeField> GTowerRing<F> {
    /// The ring of a `u32`-engine [`TowerRing`] over the backend `field`,
    /// which must hold the same prime.
    pub fn from_tower_ring(field: F, ring: &TowerRing) -> Self {
        assert_eq!(
            field.modulus(),
            &num_bigint::BigUint::from(ring.p),
            "f4_tower_generic: the backend's prime is not the ring's"
        );
        let rules = ring
            .rules
            .iter()
            .map(|r: &SquareRule| GSquareRule {
                next: r.next,
                a: field.from_u64(r.a),
                b: field.from_u64(r.b),
                c: field.from_u64(r.c),
                d: field.from_u64(r.d),
            })
            .collect();
        GTowerRing {
            field,
            rules,
            n_free: ring.n_free,
        }
    }

    pub fn n_tower(&self) -> usize {
        self.rules.len()
    }

    pub fn n_vars(&self) -> usize {
        self.rules.len() + self.n_free
    }

    /// Panics unless the rules describe a tower this engine handles.
    pub fn check(&self) {
        let f = &self.field;
        assert!(
            self.rules.len() <= 64,
            "f4_tower_generic: at most 64 tower variables"
        );
        assert!(
            self.n_free <= MAX_FREE,
            "f4_tower_generic: at most {MAX_FREE} free variables"
        );
        for (j, r) in self.rules.iter().enumerate() {
            match r.next {
                Some(nx) => assert!(
                    (nx as usize) > j && (nx as usize) < self.rules.len(),
                    "f4_tower_generic: rule {j} must point at a later tower variable"
                ),
                None => assert!(
                    f.is_zero(r.a) && f.is_zero(r.b),
                    "f4_tower_generic: top rule {j} cannot use a next level"
                ),
            }
        }
    }

    /// `coef · y^m · y_j` in normal form (tower part only), pushed onto
    /// `out` as masks with coefficients (repeats allowed).
    fn mul_var_into(&self, m: u64, j: usize, coef: F::Elem, out: &mut Vec<(u64, F::Elem)>) {
        let f = &self.field;
        let bit = 1u64 << j;
        if m & bit == 0 {
            out.push((m | bit, coef));
            return;
        }
        let r = self.rules[j];
        let m0 = m & !bit;
        if !f.is_zero(r.d) {
            out.push((m0, f.mul(coef, r.d)));
        }
        if !f.is_zero(r.c) {
            out.push((m, f.mul(coef, r.c)));
        }
        if let Some(nx) = r.next {
            let nx = nx as usize;
            if !f.is_zero(r.b) {
                self.mul_var_into(m0, nx, f.mul(coef, r.b), out);
            }
            if !f.is_zero(r.a) {
                self.mul_var_into(m, nx, f.mul(coef, r.a), out);
            }
        }
    }

    /// `coef · y^ma · y^mb` in normal form (tower parts only).
    fn mask_product_into(&self, ma: u64, mb: u64, coef: F::Elem, out: &mut Vec<(u64, F::Elem)>) {
        let common = ma & mb;
        if common == 0 {
            out.push((ma | mb, coef));
            return;
        }
        let mut cur = vec![(ma | mb, coef)];
        let mut next = Vec::new();
        let mut bits = common;
        while bits != 0 {
            let j = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            next.clear();
            for &(m, c) in &cur {
                self.mul_var_into(m, j, c, &mut next);
            }
            merge_masks(&self.field, &mut next);
            std::mem::swap(&mut cur, &mut next);
        }
        out.extend(cur);
    }

    /// `u · g` in normal form.
    pub fn mul_mono(&self, g: &GPoly<F::Elem>, u: Mono) -> GPoly<F::Elem> {
        let terms = g.monos.iter().copied().zip(g.coefs.iter().copied());
        self.mul_terms(terms, g.len(), u)
    }

    fn mul_elem(&self, g: &Elem<F::Elem>, u: Mono) -> GPoly<F::Elem> {
        self.mul_terms(g.terms(self.field.zero()), g.len(), u)
    }

    /// `u · g` for `g` given by its terms, descending, non-zero.
    fn mul_terms(
        &self,
        g: impl Iterator<Item = (Mono, F::Elem)>,
        len: usize,
        u: Mono,
    ) -> GPoly<F::Elem> {
        let (mu, fu) = (mask_of(u), free_of(u));
        let mut terms: Vec<(Mono, F::Elem)> = Vec::with_capacity(len);
        let mut carried = false;
        let mut tmp: Vec<(u64, F::Elem)> = Vec::new();
        for (m, c) in g {
            let mm = mask_of(m);
            let fr = add_free(free_of(m), fu);
            if mm & mu == 0 {
                terms.push((mono(mm | mu, fr), c));
            } else {
                carried = true;
                tmp.clear();
                self.mask_product_into(mm, mu, c, &mut tmp);
                for &(x, cc) in &tmp {
                    terms.push((mono(x, fr), cc));
                }
            }
        }
        if carried {
            GPoly::from_terms(&self.field, terms)
        } else {
            GPoly {
                monos: terms.iter().map(|t| t.0).collect(),
                coefs: terms.iter().map(|t| t.1).collect(),
            }
        }
    }

    /// The normal form of a polynomial given by exponent vectors (length
    /// [`n_vars`](Self::n_vars)) and coefficients.
    pub fn from_raw(&self, terms: &[(Vec<u32>, F::Elem)]) -> GPoly<F::Elem> {
        let f = &self.field;
        let n_t = self.n_tower();
        let mut out: Vec<(Mono, F::Elem)> = Vec::new();
        let mut cur: Vec<(u64, F::Elem)> = Vec::new();
        let mut next: Vec<(u64, F::Elem)> = Vec::new();
        for (e, c) in terms {
            assert_eq!(e.len(), self.n_vars(), "f4_tower_generic: exponent length");
            if f.is_zero(*c) {
                continue;
            }
            let mut fr = [0u8; MAX_FREE];
            for k in 0..self.n_free {
                fr[k] =
                    u8::try_from(e[n_t + k]).expect("f4_tower_generic: free exponent above 255");
            }
            cur.clear();
            cur.push((0, *c));
            for (j, &ej) in e.iter().take(n_t).enumerate() {
                for _ in 0..ej {
                    next.clear();
                    for &(m, cc) in &cur {
                        self.mul_var_into(m, j, cc, &mut next);
                    }
                    merge_masks(f, &mut next);
                    std::mem::swap(&mut cur, &mut next);
                }
            }
            for &(m, cc) in &cur {
                out.push((mono(m, fr), cc));
            }
        }
        GPoly::from_terms(f, out)
    }

    /// `g` at a point: one value per variable, tower variables first.
    pub fn eval(&self, g: &GPoly<F::Elem>, x: &[F::Elem]) -> F::Elem {
        let f = &self.field;
        let n_t = self.n_tower();
        let mut acc = f.zero();
        for (&m, &c) in g.monos.iter().zip(&g.coefs) {
            let mut t = c;
            let mut bits = mask_of(m);
            while bits != 0 {
                let j = bits.trailing_zeros() as usize;
                bits &= bits - 1;
                t = f.mul(t, x[j]);
            }
            for (k, &e) in free_of(m).iter().enumerate().take(self.n_free) {
                if e > 0 {
                    t = f.mul(t, f.pow_u64(x[n_t + k], u64::from(e)));
                }
            }
            acc = f.add(acc, t);
        }
        acc
    }
}

/// The number of standard monomials of `R / (lms)` (tower variables
/// square-free), or `None` above `bound` or when infinite: with a radical
/// zero-dimensional ideal, its number of points.
pub fn staircase(lms: &[Mono], n_tower: usize, n_free: usize, bound: usize) -> Option<usize> {
    staircase_at_most(lms, n_tower, n_free, bound)
}

// ── Basis elements ───────────────────────────────────────────────────

/// A basis element as the elimination leaves it: the coefficients of
/// `cols[start..]`, zeros included (see `f4_fp_tower::Elem`).
#[derive(Clone, Debug)]
struct Elem<E> {
    cols: Arc<[Mono]>,
    start: u32,
    vals: Vec<E>,
    len: u32,
}

impl<E: Copy + PartialEq> Elem<E> {
    fn new(zero: E, cols: Arc<[Mono]>, start: u32, vals: Vec<E>) -> Self {
        debug_assert!(vals.first().is_some_and(|&v| v != zero));
        debug_assert!(start as usize + vals.len() <= cols.len());
        let len = vals.iter().filter(|&&v| v != zero).count() as u32;
        Elem {
            cols,
            start,
            vals,
            len,
        }
    }

    fn from_poly(p: GPoly<E>) -> Self {
        let len = p.len() as u32;
        Elem {
            cols: p.monos.into(),
            start: 0,
            vals: p.coefs,
            len,
        }
    }

    fn lm(&self) -> Mono {
        self.cols[self.start as usize]
    }

    fn len(&self) -> usize {
        self.len as usize
    }

    fn terms(&self, zero: E) -> impl Iterator<Item = (Mono, E)> + '_ {
        self.cols[self.start as usize..]
            .iter()
            .zip(&self.vals)
            .filter(move |&(_, &v)| v != zero)
            .map(|(&m, &v)| (m, v))
    }

    fn to_poly(&self, zero: E) -> GPoly<E> {
        let (monos, coefs) = self.terms(zero).unzip();
        GPoly { monos, coefs }
    }
}

/// The report: the `u32` engine's counters and trace (its `basis` left
/// empty) and the basis over `F`.
#[derive(Clone, Debug)]
pub struct GTowerF4Report<E> {
    pub stats: TowerF4Report,
    /// The active basis when the run ended (`[1]` when inconsistent).
    pub basis: Vec<GPoly<E>>,
}

// ── Sparse rows and the eliminations ─────────────────────────────────

/// A sparse row: ascending column indices and non-zero coefficients.
#[derive(Clone, Debug)]
struct Row<E> {
    cols: Vec<u32>,
    vals: Vec<E>,
}

/// A dense lazy accumulator with a bitset of the columns it may be
/// non-zero in (see `f4_fp_tower::Acc`).
struct Acc<W> {
    vals: Vec<W>,
    bits: Vec<u64>,
    hi: usize,
}

impl<W: Copy> Acc<W> {
    fn new(n_cols: usize, zero: W) -> Self {
        Acc {
            vals: vec![zero; n_cols],
            bits: vec![0; n_cols.div_ceil(64).max(1)],
            hi: 0,
        }
    }
}

fn acc_load<F: PrimeField>(f: &F, acc: &mut Acc<F::Wide>, row: &Row<F::Elem>) {
    for (&c, &v) in row.cols.iter().zip(&row.vals) {
        acc.vals[c as usize] = f.wide_from(v);
        acc.bits[c as usize / 64] |= 1u64 << (c % 64);
    }
    acc.hi = row.cols.last().map_or(0, |&c| c as usize / 64);
}

/// Subtract `v · row` (without its lead): one scattered lazy product per
/// entry, then the bitset in a separate pass.
fn acc_sub_tail<F: PrimeField>(
    f: &F,
    acc: &mut Acc<F::Wide>,
    v: F::Elem,
    row: &Row<F::Elem>,
    cur_word: usize,
    word: &mut u64,
) -> u64 {
    f.row_scatter_axpy(&mut acc.vals, f.neg(v), &row.cols[1..], &row.vals[1..]);
    for &c in &row.cols[1..] {
        let ci = c as usize;
        let (w, b) = (ci / 64, 1u64 << (ci % 64));
        if w == cur_word {
            *word |= b;
        } else {
            acc.bits[w] |= b;
        }
    }
    if let Some(&c) = row.cols.last() {
        acc.hi = acc.hi.max(c as usize / 64);
    }
    row.cols.len() as u64 - 1
}

/// Move every non-zero entry from word `w` on into `out`, scaled.
fn acc_drain_from<F: PrimeField>(
    f: &F,
    acc: &mut Acc<F::Wide>,
    w: usize,
    scale: F::Elem,
    out: &mut Row<F::Elem>,
) {
    let mut w = w;
    let wz = f.wide_zero();
    while w <= acc.hi {
        let mut word = acc.bits[w];
        acc.bits[w] = 0;
        while word != 0 {
            let b = word.trailing_zeros() as usize;
            word &= word - 1;
            let c = w * 64 + b;
            let v = f.wide_reduce(&acc.vals[c]);
            acc.vals[c] = wz;
            if !f.is_zero(v) {
                out.cols.push(c as u32);
                out.vals.push(f.mul(v, scale));
            }
        }
        w += 1;
    }
}

/// Echelon form of `rows` (the input echelon), as in `f4_fp_tower`.
fn echelon<F: PrimeField>(
    f: &F,
    mut rows: Vec<Row<F::Elem>>,
    n_cols: usize,
    muladds: &mut u64,
    deadline: Option<Instant>,
) -> Option<Vec<Row<F::Elem>>> {
    rows.retain(|r| !r.cols.is_empty());
    rows.sort_by_key(|r| (r.cols[0], r.cols.len()));
    let mut pivot_of = vec![NONE; n_cols];
    let mut pivots: Vec<Row<F::Elem>> = Vec::new();
    let mut acc = Acc::new(n_cols, f.wide_zero());
    let wz = f.wide_zero();
    let cap = f.lazy_capacity();
    for (k, row) in rows.into_iter().enumerate() {
        if k % 256 == 0 && deadline.is_some_and(|d| Instant::now() >= d) {
            return None;
        }
        if pivot_of[row.cols[0] as usize] == NONE {
            let inv = f.inv(row.vals[0]).expect("a row's lead is non-zero");
            let vals = row.vals.iter().map(|&v| f.mul(v, inv)).collect();
            pivot_of[row.cols[0] as usize] = pivots.len() as u32;
            pivots.push(Row {
                cols: row.cols,
                vals,
            });
            continue;
        }
        acc_load(f, &mut acc, &row);
        let mut used = 0u64;
        let mut w = row.cols[0] as usize / 64;
        let mut found: Option<Row<F::Elem>> = None;
        'scan: while w <= acc.hi {
            let mut word = acc.bits[w];
            acc.bits[w] = 0;
            while word != 0 {
                let b = word.trailing_zeros() as usize;
                word &= word - 1;
                let c = w * 64 + b;
                let v = f.wide_reduce(&acc.vals[c]);
                acc.vals[c] = wz;
                if f.is_zero(v) {
                    continue;
                }
                let p = pivot_of[c];
                if p != NONE {
                    if used >= cap {
                        // Every entry outside the bitset is zero.
                        renormalise(f, &mut acc.vals);
                        used = 0;
                    }
                    used += 1;
                    *muladds += acc_sub_tail(f, &mut acc, v, &pivots[p as usize], w, &mut word);
                    continue;
                }
                let inv = f.inv(v).expect("non-zero");
                let mut out = Row {
                    cols: vec![c as u32],
                    vals: vec![f.one()],
                };
                acc.bits[w] = word;
                acc_drain_from(f, &mut acc, w, inv, &mut out);
                found = Some(out);
                break 'scan;
            }
            w += 1;
        }
        if let Some(out) = found {
            pivot_of[out.cols[0] as usize] = pivots.len() as u32;
            pivots.push(out);
        }
    }
    Some(pivots)
}

enum Stop {
    Deadline,
    Oversize,
}

#[derive(Default)]
struct ElimStats {
    residual_rows: usize,
    full_rank_after: Option<usize>,
    skipped_rows: usize,
    q_cols: usize,
    dense_entries: u64,
    muladds: u64,
    ms_reducers: f64,
    ms_residues: f64,
}

/// Forms row `k` of a set of rows on demand.
type FormRow<'a, E> = &'a (dyn Fn(usize) -> Row<E> + Sync);

/// The pivots of an elimination (see `f4_fp_tower::Echelon`).
struct Echelon<E> {
    q_col: Vec<u32>,
    leads: Vec<u32>,
    rows: Vec<Vec<E>>,
}

/// `row[off..] −= v·b` (lazy), reducing and reloading `row` first when it
/// has taken `cap` products since its load (`used` counts them); the
/// multiply-adds it counts.
#[inline]
fn sub_row<F: PrimeField>(
    f: &F,
    row: &mut [F::Wide],
    off: usize,
    used: &mut u64,
    cap: u64,
    v: F::Elem,
    b: &[F::Elem],
) -> u64 {
    if *used >= cap {
        renormalise(f, row);
        *used = 0;
    }
    *used += 1;
    f.row_sub(&mut row[off..], v, b);
    b.len() as u64
}

/// Reduce a row of accumulators and reload it, so that its product count
/// starts again; returns whether any entry is non-zero.
#[inline]
fn renormalise<F: PrimeField>(f: &F, row: &mut [F::Wide]) -> bool {
    let mut any = false;
    for x in row.iter_mut() {
        let v = f.wide_reduce(x);
        any |= !f.is_zero(v);
        *x = f.wide_from(v);
    }
    any
}

/// One step's matrix, eliminated in the manner of Faugère–Lachartre; the
/// `u32` engine's `eliminate` over `F` (its documentation applies).
#[allow(clippy::too_many_arguments)]
fn eliminate<F: PrimeField>(
    f: &F,
    n_cols: usize,
    piv_leads: &[u32],
    piv_row: FormRow<F::Elem>,
    pivot_of: &[u32],
    n_s: usize,
    s_row: FormRow<F::Elem>,
    deadline: Option<Instant>,
    max_dense: Option<u64>,
    full_rank_exit: bool,
    st: &mut ElimStats,
) -> Result<Echelon<F::Elem>, Stop> {
    let wz = f.wide_zero();
    let cap = f.lazy_capacity();
    let mut q_of = vec![NONE; n_cols];
    let mut q_start = vec![0u32; n_cols];
    let mut q_col: Vec<u32> = Vec::new();
    for c in 0..n_cols {
        q_start[c] = q_col.len() as u32;
        if pivot_of[c] == NONE {
            q_of[c] = q_col.len() as u32;
            q_col.push(c as u32);
        }
    }
    let nq = q_col.len();
    st.q_cols = nq;
    let tail_len = |c: u32| nq - q_start[c as usize] as usize;
    st.dense_entries = piv_leads.iter().map(|&c| tail_len(c) as u64).sum();
    if max_dense.is_some_and(|cap| st.dense_entries > cap) {
        return Err(Stop::Oversize);
    }
    let expired = AtomicBool::new(false);
    let past = || {
        if expired.load(AtomicOrdering::Relaxed) {
            return true;
        }
        if deadline.is_some_and(|d| Instant::now() >= d) {
            expired.store(true, AtomicOrdering::Relaxed);
            return true;
        }
        false
    };

    // 1. B' = A⁻¹B, the rightmost lead first.
    let started = Instant::now();
    let n_piv = piv_leads.len();
    let mut order: Vec<u32> = (0..n_piv as u32).collect();
    order.sort_unstable_by_key(|&k| std::cmp::Reverse(piv_leads[k as usize]));
    let mut pos = vec![0u32; n_piv];
    for (i, &k) in order.iter().enumerate() {
        pos[k as usize] = i as u32;
    }
    let mut bred: Vec<Vec<F::Elem>> = vec![Vec::new(); n_piv];
    const BLOCK: usize = 64;
    for (b, block) in order.chunks(BLOCK).enumerate() {
        let first = (b * BLOCK) as u32;
        type Partial<W, E> = (Vec<W>, Vec<(u32, E)>, u64, u64);
        let partial: Vec<Partial<F::Wide, F::Elem>> = block
            .par_iter()
            .map(|&k| {
                let s0 = q_start[piv_leads[k as usize] as usize] as usize;
                let mut acc = vec![wz; nq - s0];
                let mut later = Vec::new();
                let mut work = 0u64;
                let mut used = 0u64;
                if past() {
                    return (acc, later, work, used);
                }
                let row = piv_row(k as usize);
                debug_assert_eq!(row.cols[0], piv_leads[k as usize]);
                debug_assert!(f.is_one(row.vals[0]), "a pivot row must be monic");
                for (&c, &v) in row.cols[1..].iter().zip(&row.vals[1..]) {
                    if q_of[c as usize] != NONE {
                        acc[q_of[c as usize] as usize - s0] = f.wide_from(v);
                    }
                }
                for (&c, &v) in row.cols[1..].iter().zip(&row.vals[1..]) {
                    let cu = c as usize;
                    if q_of[cu] != NONE {
                        continue;
                    }
                    let j = pivot_of[cu];
                    if pos[j as usize] < first {
                        let off = q_start[cu] as usize - s0;
                        work += sub_row(f, &mut acc, off, &mut used, cap, v, &bred[j as usize]);
                    } else {
                        later.push((j, v));
                    }
                }
                (acc, later, work, used)
            })
            .collect();
        if expired.load(AtomicOrdering::Relaxed) {
            return Err(Stop::Deadline);
        }
        for (&k, (mut acc, later, work, mut used)) in block.iter().zip(partial) {
            st.muladds += work;
            let s0 = q_start[piv_leads[k as usize] as usize] as usize;
            for (j, v) in later {
                let off = q_start[piv_leads[j as usize] as usize] as usize - s0;
                st.muladds += sub_row(f, &mut acc, off, &mut used, cap, v, &bred[j as usize]);
            }
            let mut out = vec![f.zero(); acc.len()];
            f.row_reduce(&acc, &mut out);
            bred[k as usize] = out;
        }
    }
    let residues_started = Instant::now();
    st.ms_reducers = (residues_started - started).as_secs_f64() * 1e3;

    // 2 and 3. The S-rows, in chunks.
    let mut ech: Vec<Vec<F::Elem>> = vec![Vec::new(); nq];
    let mut leads: Vec<u32> = Vec::new();
    // The `u32` engine's block and chunk sizes, whatever the accumulator's
    // width: they decide which pivots reduce which rows first, so they set
    // the multiply-add count, and a report must not depend on the backend.
    let per_task = ((1usize << 18) / (8 * nq.max(1))).clamp(4, 64);
    let chunk = per_task * 2 * rayon::current_num_threads().max(1);
    type Reduced<W> = (Vec<(Option<Vec<W>>, bool)>, u64);
    let all: Vec<usize> = (0..n_s).collect();
    let mut done = 0usize;
    if nq == 0 {
        st.full_rank_after = Some(0);
    }
    for rows in all.chunks(chunk) {
        if full_rank_exit && leads.len() == nq {
            st.skipped_rows += n_s - done;
            break;
        }
        let reduced: Vec<Reduced<F::Wide>> = rows
            .par_chunks(per_task)
            .map(|ks| {
                if past() {
                    return (Vec::new(), 0);
                }
                let blk: Vec<Row<F::Elem>> = ks.iter().map(|&k| s_row(k)).collect();
                let mut acc = vec![wz; blk.len() * nq];
                let mut used = vec![0u64; blk.len()];
                let mut work = 0u64;
                let mut ents: Vec<(u32, u32, F::Elem)> = Vec::new();
                for (r, row) in blk.iter().enumerate() {
                    for (&c, &v) in row.cols.iter().zip(&row.vals) {
                        let q = q_of[c as usize];
                        if q != NONE {
                            acc[r * nq + q as usize] = f.wide_from(v);
                        } else {
                            ents.push((c, r as u32, v));
                        }
                    }
                }
                // (column, row) is unique, so this is the `u32` engine's order.
                ents.sort_unstable_by_key(|e| (e.0, e.1));
                for &(c, r, v) in &ents {
                    let r = r as usize;
                    work += sub_row(
                        f,
                        &mut acc[r * nq..(r + 1) * nq],
                        q_start[c as usize] as usize,
                        &mut used[r],
                        cap,
                        v,
                        &bred[pivot_of[c as usize] as usize],
                    );
                }
                let residual: Vec<bool> = (0..blk.len())
                    .map(|r| {
                        acc[r * nq..(r + 1) * nq]
                            .iter()
                            .any(|x| !f.is_zero(f.wide_reduce(x)))
                    })
                    .collect();
                let mut live: Vec<usize> = (0..blk.len()).filter(|&r| residual[r]).collect();
                for &l in &leads {
                    let l = l as usize;
                    for &r in &live {
                        let a = f.wide_reduce(&acc[r * nq + l]);
                        if !f.is_zero(a) {
                            let row = &mut acc[r * nq..(r + 1) * nq];
                            work += sub_row(f, row, l, &mut used[r], cap, a, &ech[l]);
                        }
                    }
                }
                let mut out: Vec<(Option<Vec<F::Wide>>, bool)> =
                    residual.iter().map(|&x| (None, x)).collect();
                live.retain(|&r| renormalise(f, &mut acc[r * nq..(r + 1) * nq]));
                for r in live {
                    out[r].0 = Some(acc[r * nq..(r + 1) * nq].to_vec());
                }
                (out, work)
            })
            .collect();
        if expired.load(AtomicOrdering::Relaxed) {
            return Err(Stop::Deadline);
        }
        let mut new_leads: Vec<u32> = Vec::new();
        let mut residues: Vec<(Option<Vec<F::Wide>>, bool)> = Vec::with_capacity(rows.len());
        for (out, work) in reduced {
            st.muladds += work;
            residues.extend(out);
        }
        for (i, (acc, residual)) in residues.into_iter().enumerate() {
            if full_rank_exit && leads.len() + new_leads.len() == nq {
                st.skipped_rows += rows.len() - i;
                break;
            }
            st.residual_rows += usize::from(residual);
            let Some(mut acc) = acc else { continue };
            // Renormalised when it left its block: no products yet.
            let mut used = 0u64;
            for &l in &new_leads {
                let l = l as usize;
                let a = f.wide_reduce(&acc[l]);
                if !f.is_zero(a) {
                    st.muladds += sub_row(f, &mut acc, l, &mut used, cap, a, &ech[l]);
                }
            }
            let Some(l) = acc.iter().position(|x| !f.is_zero(f.wide_reduce(x))) else {
                continue;
            };
            let inv = f.inv(f.wide_reduce(&acc[l])).expect("non-zero");
            let mut row = vec![f.zero(); nq - l];
            f.row_reduce_scaled(&acc[l..], inv, &mut row);
            ech[l] = row;
            let at = new_leads.partition_point(|&m| m < l as u32);
            new_leads.insert(at, l as u32);
            if leads.len() + new_leads.len() == nq {
                st.full_rank_after = Some(done + i + 1);
            }
        }
        done += rows.len();
        leads.extend(new_leads);
        leads.sort_unstable();
    }
    st.ms_residues = residues_started.elapsed().as_secs_f64() * 1e3;

    Ok(Echelon {
        q_col,
        leads,
        rows: ech,
    })
}

// ── Pairs and the basis ──────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
enum PairKind {
    Critical(u32, u32),
    Tower(u32, u8),
}

#[derive(Clone, Copy, Debug)]
struct Pair {
    kind: PairKind,
    lcm: Mono,
    deg: u32,
}

struct State<E> {
    polys: Vec<Elem<E>>,
    lm: Vec<Mono>,
    active: Vec<bool>,
    lm_index: MonoMap<u32>,
    pairs: Vec<Pair>,
}

impl<E: Copy + PartialEq> State<E> {
    /// Becker–Weispfenning `UPDATE` plus the tower pairs (as in
    /// `f4_fp_tower`).
    fn insert(&mut self, h_poly: Elem<E>, rep: &mut TowerF4Report) {
        let h = self.polys.len();
        let lh = h_poly.lm();
        self.polys.push(h_poly);
        self.lm.push(lh);
        self.active.push(false);

        let dh = degree_of(lh);
        let mut bits = mask_of(lh);
        while bits != 0 {
            let v = bits.trailing_zeros();
            bits &= bits - 1;
            self.pairs.push(Pair {
                kind: PairKind::Tower(h as u32, v as u8),
                lcm: lh,
                deg: dh + 1,
            });
        }
        let mut groups: MonoMap<(u32, bool)> = MonoMap::default();
        for g in 0..h {
            if !self.active[g] {
                continue;
            }
            let l = lcm(lh, self.lm[g]);
            let e = groups.entry(l).or_insert((g as u32, false));
            if coprime(lh, self.lm[g]) {
                e.1 = true;
            }
        }
        let mut kept: Vec<(u32, Mono)> = Vec::new();
        for (&l, &(g, has_coprime)) in &groups {
            if any_between(lh, l, |c| groups.contains_key(&c)) {
                rep.pairs_chain_skipped += 1;
                continue;
            }
            if has_coprime {
                rep.pairs_product_skipped += 1;
                continue;
            }
            kept.push((g, l));
        }
        let lm = &self.lm;
        let mut dropped = 0u64;
        self.pairs.retain(|p| match p.kind {
            PairKind::Tower(..) => true,
            PairKind::Critical(i, j) => {
                let keep = !divides(lh, p.lcm)
                    || lcm(lm[i as usize], lh) == p.lcm
                    || lcm(lm[j as usize], lh) == p.lcm;
                if !keep {
                    dropped += 1;
                }
                keep
            }
        });
        rep.pairs_chain_skipped += dropped;
        kept.sort_unstable_by_key(|&(g, _)| g);
        for (g, l) in kept {
            self.pairs.push(Pair {
                kind: PairKind::Critical(g, h as u32),
                lcm: l,
                deg: degree_of(l),
            });
        }
        for g in 0..h {
            if self.active[g] && divides(lh, self.lm[g]) {
                self.active[g] = false;
                if self.lm_index.get(&self.lm[g]) == Some(&(g as u32)) {
                    self.lm_index.remove(&self.lm[g]);
                }
            }
        }
        self.active[h] = true;
        self.lm_index.insert(lh, h as u32);
    }

    /// An active element whose leading monomial divides `m`, the shortest.
    fn reducer_for(&self, m: Mono) -> Option<u32> {
        let mask = mask_of(m);
        let fm = free_of(m);
        let mut best: Option<u32> = None;
        let mut s = mask;
        loop {
            let mut fv = [0u8; MAX_FREE];
            loop {
                if let Some(&g) = self.lm_index.get(&mono(s, fv)) {
                    if best
                        .is_none_or(|b| self.polys[g as usize].len() < self.polys[b as usize].len())
                    {
                        best = Some(g);
                    }
                }
                let mut k = 0;
                loop {
                    if k == MAX_FREE {
                        break;
                    }
                    if fv[k] < fm[k] {
                        fv[k] += 1;
                        break;
                    }
                    fv[k] = 0;
                    k += 1;
                }
                if k == MAX_FREE {
                    break;
                }
            }
            if s == 0 {
                break;
            }
            s = (s - 1) & mask;
        }
        best
    }

    fn active_lms(&self) -> Vec<Mono> {
        (0..self.polys.len())
            .filter(|&g| self.active[g])
            .map(|g| self.lm[g])
            .collect()
    }
}

#[derive(Clone, Copy, Debug)]
struct RowDesc {
    mult: Mono,
    g: u32,
    lead: Mono,
    len: u32,
}

#[allow(clippy::too_many_arguments)]
fn scan_products<F: PrimeField>(
    ring: &GTowerRing<F>,
    polys: &[Elem<F::Elem>],
    todo: &[(Mono, u32)],
    skip: usize,
    out: &mut Vec<RowDesc>,
    examined: &mut MonoMap<()>,
    frontier: &mut Vec<Mono>,
) -> u64 {
    let mut terms = 0u64;
    for batch in todo.chunks(1024) {
        let prods: Vec<GPoly<F::Elem>> = batch
            .par_iter()
            .map(|&(mult, g)| ring.mul_elem(&polys[g as usize], mult))
            .collect();
        for (&(mult, g), p) in batch.iter().zip(&prods) {
            let Some(lead) = p.lm() else { continue };
            terms += p.len() as u64;
            out.push(RowDesc {
                mult,
                g,
                lead,
                len: p.len() as u32,
            });
            for &m in p.monos.iter().skip(skip) {
                if examined.insert(m, ()).is_none() {
                    frontier.push(m);
                }
            }
        }
    }
    terms
}

// ── The algorithm ────────────────────────────────────────────────────

/// **F4 on `R` over `F`** for the ideal generated by `input` (tower normal
/// forms): `f4_fp_tower::f4_tower` with the arithmetic of the backend.
pub fn f4_tower_generic<F: PrimeField>(
    input: &[GPoly<F::Elem>],
    ring: &GTowerRing<F>,
    opts: &TowerF4Options,
) -> GTowerF4Report<F::Elem> {
    ring.check();
    let f = &ring.field;
    assert!(
        f.lazy_capacity() >= 1,
        "f4_tower_generic: the backend's lazy capacity is zero"
    );
    let zero = f.zero();
    let started = Instant::now();
    let mut rep = TowerF4Report::default();
    let finish = |mut rep: TowerF4Report, basis: Vec<GPoly<F::Elem>>| {
        rep.ms = started.elapsed().as_secs_f64() * 1e3;
        GTowerF4Report { stats: rep, basis }
    };
    let unit = || {
        vec![GPoly {
            monos: vec![ONE],
            coefs: vec![f.one()],
        }]
    };

    let inputs: Vec<GPoly<F::Elem>> = input.iter().filter(|p| !p.is_empty()).cloned().collect();
    if inputs.iter().any(|p| p.is_unit()) {
        rep.inconsistent = true;
        return finish(rep, unit());
    }

    let start = {
        let mut all: Vec<Mono> = inputs
            .iter()
            .flat_map(|p| p.monos.iter().copied())
            .collect();
        all.sort_unstable_by(|a, b| b.cmp(a));
        all.dedup();
        let col_of: MonoMap<u32> = all
            .iter()
            .enumerate()
            .map(|(i, &m)| (m, i as u32))
            .collect();
        let rows: Vec<Row<F::Elem>> = inputs
            .iter()
            .map(|p| Row {
                cols: p.monos.iter().map(|m| col_of[m]).collect(),
                vals: p.coefs.clone(),
            })
            .collect();
        match echelon(f, rows, all.len(), &mut rep.muladds, opts.deadline) {
            Some(piv) => piv
                .into_iter()
                .map(|r| GPoly {
                    monos: r.cols.iter().map(|&c| all[c as usize]).collect(),
                    coefs: r.vals,
                })
                .collect::<Vec<_>>(),
            None => {
                rep.timed_out = true;
                return finish(rep, inputs);
            }
        }
    };
    if start.iter().any(|p| p.is_unit()) {
        rep.inconsistent = true;
        return finish(rep, unit());
    }
    let mut start = start;
    start.sort_by_key(|p| std::cmp::Reverse(p.lm().unwrap()));

    let mut s: State<F::Elem> = State {
        polys: Vec::new(),
        lm: Vec::new(),
        active: Vec::new(),
        lm_index: MonoMap::default(),
        pairs: Vec::new(),
    };
    for p in start {
        s.insert(Elem::from_poly(p), &mut rep);
    }
    let mut learned = (0u32, 0usize, 0usize);
    let mut max_cols_so_far = 0usize;

    let active_basis = |s: &State<F::Elem>| -> Vec<GPoly<F::Elem>> {
        (0..s.polys.len())
            .filter(|&g| s.active[g])
            .map(|g| s.polys[g].to_poly(zero))
            .collect()
    };
    let kept = |s: &State<F::Elem>, tr: &mut StepTrace| {
        tr.basis_active = s.active.iter().filter(|&&a| a).count();
        tr.basis_kept = s.polys.len();
        tr.basis_entries = s.polys.iter().map(|p| p.vals.len() as u64).sum();
        tr.pairs_left = s.pairs.len();
    };
    let record = |rep: &mut TowerF4Report, learned: (u32, usize, usize)| {
        rep.solving_degree_max = learned.0;
        rep.max_cols_to_solution = learned.1;
        rep.steps_to_solution = learned.2;
    };

    while !s.pairs.is_empty() {
        if opts.deadline.is_some_and(|d| Instant::now() >= d) {
            rep.timed_out = true;
            break;
        }
        let d = s.pairs.iter().map(|p| p.deg).min().unwrap();
        if d > opts.max_degree {
            rep.pairs_above_bound = s.pairs.len();
            break;
        }
        let step_started = Instant::now();
        let (selected, rest): (Vec<Pair>, Vec<Pair>) = s.pairs.drain(..).partition(|p| p.deg == d);
        s.pairs = rest;
        rep.steps += 1;
        rep.degree_reached = rep.degree_reached.max(d);
        let mut tr = StepTrace {
            degree: d,
            ..Default::default()
        };

        let mut seen: HashMap<(Mono, u32), (), BuildHasherDefault<MonoHasher>> = HashMap::default();
        let mut jobs: Vec<(Mono, u32)> = Vec::new();
        let mut lcm_cols: MonoMap<()> = MonoMap::default();
        for p in &selected {
            match p.kind {
                PairKind::Critical(i, j) => {
                    tr.critical_pairs += 1;
                    lcm_cols.insert(p.lcm, ());
                    for g in [i, j] {
                        let mult = quotient(p.lcm, s.lm[g as usize]);
                        if seen.insert((mult, g), ()).is_none() {
                            jobs.push((mult, g));
                        }
                    }
                }
                PairKind::Tower(g, v) => {
                    tr.tower_pairs += 1;
                    let mult = mono(1u64 << v, [0; MAX_FREE]);
                    if seen.insert((mult, g), ()).is_none() {
                        jobs.push((mult, g));
                    }
                }
            }
        }
        rep.critical_pairs_reduced += tr.critical_pairs as u64;
        rep.tower_pairs_reduced += tr.tower_pairs as u64;

        let mut examined: MonoMap<()> = lcm_cols.clone();
        let mut frontier: Vec<Mono> = Vec::new();
        let mut s_desc: Vec<RowDesc> = Vec::with_capacity(jobs.len());
        let mut nnz = scan_products(
            ring,
            &s.polys,
            &jobs,
            0,
            &mut s_desc,
            &mut examined,
            &mut frontier,
        );
        let mut r_desc: Vec<RowDesc> = Vec::new();
        let mut aborted = false;
        while !frontier.is_empty() {
            if opts.deadline.is_some_and(|dl| Instant::now() >= dl) {
                aborted = true;
                break;
            }
            let found: Vec<(Mono, Option<u32>)> = frontier
                .par_iter()
                .map(|&m| (m, s.reducer_for(m)))
                .collect();
            frontier.clear();
            let mut todo: Vec<(Mono, u32)> = Vec::new();
            for &(m, g) in &found {
                if let Some(g) = g {
                    todo.push((quotient(m, s.lm[g as usize]), g));
                }
            }
            nnz += scan_products(
                ring,
                &s.polys,
                &todo,
                1,
                &mut r_desc,
                &mut examined,
                &mut frontier,
            );
            if opts.max_nnz.is_some_and(|cap| nnz > cap) {
                rep.oversize = true;
                aborted = true;
                break;
            }
        }
        if aborted {
            if !rep.oversize {
                rep.timed_out = true;
            }
            s.pairs.extend(selected);
            tr.nnz = nnz;
            kept(&s, &mut tr);
            tr.ms = step_started.elapsed().as_secs_f64() * 1e3;
            if let Some(cb) = opts.on_stop {
                cb(rep.steps, rep.oversize, &tr);
            }
            rep.stopped_step = Some(tr);
            break;
        }
        rep.reducer_rows += r_desc.len() as u64;

        let mut cols: Vec<Mono> = examined.keys().copied().collect();
        drop(examined);
        cols.par_sort_unstable_by(|a, b| b.cmp(a));
        let col_of: MonoMap<u32> = cols
            .iter()
            .enumerate()
            .map(|(i, &m)| (m, i as u32))
            .collect();
        let n_cols = cols.len();
        let mut pivot_of = vec![NONE; n_cols];
        let mut piv: Vec<RowDesc> = r_desc;
        for (k, d) in piv.iter().enumerate() {
            pivot_of[col_of[&d.lead] as usize] = k as u32;
        }
        let n_reducers = piv.len();
        s_desc.sort_by_key(|d| (std::cmp::Reverse(d.lead), d.len));
        let mut s_list: Vec<RowDesc> = Vec::with_capacity(s_desc.len());
        for d in s_desc {
            let c = col_of[&d.lead] as usize;
            if pivot_of[c] == NONE && lcm_cols.contains_key(&d.lead) {
                pivot_of[c] = piv.len() as u32;
                piv.push(d);
            } else {
                s_list.push(d);
            }
        }
        tr.promoted_rows = piv.len() - n_reducers;
        tr.s_rows = s_list.len();
        tr.reducer_rows = n_reducers;
        tr.cols = n_cols;
        tr.nnz = nnz;
        rep.max_rows = rep.max_rows.max(s_list.len() + piv.len());
        rep.max_cols = rep.max_cols.max(n_cols);
        rep.max_nnz = rep.max_nnz.max(nnz);
        max_cols_so_far = max_cols_so_far.max(n_cols);

        tr.ms_rows = step_started.elapsed().as_secs_f64() * 1e3;

        let polys = &s.polys;
        let form = |d: &RowDesc| -> Row<F::Elem> {
            let p = ring.mul_elem(&polys[d.g as usize], d.mult);
            debug_assert_eq!(p.len(), d.len as usize);
            Row {
                cols: p.monos.iter().map(|m| col_of[m]).collect(),
                vals: p.coefs,
            }
        };
        let piv_row = |k: usize| -> Row<F::Elem> {
            let mut r = form(&piv[k]);
            if !f.is_one(r.vals[0]) {
                debug_assert!(k >= n_reducers);
                let inv = f.inv(r.vals[0]).expect("non-zero lead");
                for v in &mut r.vals {
                    *v = f.mul(*v, inv);
                }
            }
            r
        };
        let s_row = |k: usize| form(&s_list[k]);
        let piv_leads: Vec<u32> = piv.iter().map(|d| col_of[&d.lead]).collect();
        let mut st = ElimStats::default();
        let result = eliminate(
            f,
            n_cols,
            &piv_leads,
            &piv_row,
            &pivot_of,
            s_list.len(),
            &s_row,
            opts.deadline,
            opts.max_dense.or(opts.max_nnz),
            opts.full_rank_exit,
            &mut st,
        );
        tr.muladds = st.muladds;
        rep.muladds += st.muladds;
        tr.residual_rows = st.residual_rows;
        tr.full_rank_after = st.full_rank_after;
        tr.skipped_rows = st.skipped_rows;
        rep.rows_skipped += st.skipped_rows as u64;
        tr.residual_cols = st.q_cols;
        tr.dense_entries = st.dense_entries;
        tr.ms_reduce = st.ms_reducers;
        tr.ms_echelon = st.ms_residues;
        rep.max_residual_rows = rep.max_residual_rows.max(st.residual_rows);
        rep.max_dense_entries = rep.max_dense_entries.max(st.dense_entries);
        let ech = match result {
            Ok(ech) => ech,
            Err(stop) => {
                match stop {
                    Stop::Deadline => rep.timed_out = true,
                    Stop::Oversize => rep.oversize = true,
                }
                s.pairs.extend(selected);
                kept(&s, &mut tr);
                tr.ms = step_started.elapsed().as_secs_f64() * 1e3;
                if let Some(cb) = opts.on_stop {
                    cb(rep.steps, rep.oversize, &tr);
                }
                rep.stopped_step = Some(tr);
                break;
            }
        };
        let update_started = Instant::now();

        let Echelon {
            q_col,
            leads,
            rows: mut ech_rows,
        } = ech;
        let q_monos: Arc<[Mono]> = q_col.iter().map(|&c| cols[c as usize]).collect();
        let mut fresh: Vec<Elem<F::Elem>> = leads
            .iter()
            .map(|&l| {
                Elem::new(
                    zero,
                    q_monos.clone(),
                    l,
                    std::mem::take(&mut ech_rows[l as usize]),
                )
            })
            .collect();
        drop(ech_rows);
        tr.fresh = fresh.len();
        tr.fresh_min_degree = fresh.iter().map(|p| degree_of(p.lm())).min().unwrap_or(0);
        if !fresh.is_empty() {
            learned = (learned.0.max(d), max_cols_so_far, rep.steps);
            rep.last_productive_degree = d;
        }
        if fresh.iter().any(|p| p.lm() == ONE) {
            rep.inconsistent = true;
            record(&mut rep, learned);
            kept(&s, &mut tr);
            tr.ms = step_started.elapsed().as_secs_f64() * 1e3;
            if let Some(cb) = opts.on_step {
                cb(rep.steps, &tr);
            }
            rep.trace.push(tr);
            return finish(rep, unit());
        }
        fresh.sort_by_key(|p| std::cmp::Reverse(p.lm()));
        for p in fresh {
            debug_assert!(f.is_one(p.vals[0]), "an echelon row is monic");
            s.insert(p, &mut rep);
        }
        kept(&s, &mut tr);
        tr.ms_update = update_started.elapsed().as_secs_f64() * 1e3;
        tr.ms = step_started.elapsed().as_secs_f64() * 1e3;
        if let Some(cb) = opts.on_step {
            cb(rep.steps, &tr);
        }
        rep.trace.push(tr);

        if let Some(bound) = opts.stop_staircase {
            if !s.pairs.is_empty() {
                if let Some(count) =
                    staircase_at_most(&s.active_lms(), ring.n_tower(), ring.n_free, bound)
                {
                    rep.staircase_at_stop = Some(count);
                    break;
                }
            }
        }
    }
    record(&mut rep, learned);
    let basis = active_basis(&s);
    finish(rep, basis)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::f4_fp_tower::f4_tower;
    use crate::cryptanalysis::fp_backend::{Fp224, Fp32, Fp64, FpWide, SqrtPre};
    use crate::cryptanalysis::fp_backend_ifma::Fp52;
    use num_bigint::BigUint;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    /// The report with every wall-clock field zeroed and the basis left
    /// out, as text: two runs of one computation must agree on it.
    fn fingerprint(rep: &TowerF4Report) -> String {
        let mut r = rep.clone();
        r.ms = 0.0;
        r.basis.clear();
        for t in r.trace.iter_mut().chain(r.stopped_step.iter_mut()) {
            t.ms = 0.0;
            t.ms_rows = 0.0;
            t.ms_reduce = 0.0;
            t.ms_echelon = 0.0;
            t.ms_update = 0.0;
        }
        format!("{r:?}")
    }

    /// The basis with canonical integer coefficients.
    fn canonical<F: PrimeField>(f: &F, basis: &[GPoly<F::Elem>]) -> Vec<(Vec<Mono>, Vec<BigUint>)> {
        basis
            .iter()
            .map(|g| {
                (
                    g.monos.clone(),
                    g.coefs.iter().map(|&c| f.to_biguint(c)).collect(),
                )
            })
            .collect()
    }

    fn is_prime_u64(n: u64) -> bool {
        if n < 2 {
            return false;
        }
        for q in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
            if n.is_multiple_of(q) {
                return n == q;
            }
        }
        let mulm = |a: u64, b: u64| ((u128::from(a) * u128::from(b)) % u128::from(n)) as u64;
        let powm = |mut a: u64, mut e: u64| {
            let mut r = 1u64;
            while e > 0 {
                if e & 1 == 1 {
                    r = mulm(r, a);
                }
                a = mulm(a, a);
                e >>= 1;
            }
            r
        };
        let (mut d, mut s) = (n - 1, 0);
        while d % 2 == 0 {
            d /= 2;
            s += 1;
        }
        'w: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
            let mut x = powm(a, d);
            if x == 1 || x == n - 1 {
                continue;
            }
            for _ in 1..s {
                x = mulm(x, x);
                if x == n - 1 {
                    continue 'w;
                }
            }
            return false;
        }
        true
    }

    /// The smallest prime `c·2^k + 1` with `c` odd.
    fn proth_prime(k: u32) -> u64 {
        (1u64..)
            .step_by(2)
            .map(|c| (c << k) + 1)
            .find(|&p| is_prime_u64(p))
            .unwrap()
    }

    fn p224() -> BigUint {
        (BigUint::from(1u32) << 224u32) - (BigUint::from(1u32) << 96u32) + 1u32
    }

    /// A Kummer block of `t` levels from variable `base` over `f` (needs
    /// `2^t | p − 1`), and its `2^t` points.
    fn kummer_block<F: PrimeField>(
        f: &F,
        t: usize,
        base: usize,
        rng: &mut StdRng,
    ) -> (Vec<GSquareRule<F::Elem>>, Vec<Vec<F::Elem>>) {
        let pm1 = f.modulus() - 1u32;
        let e = &pm1 >> t;
        let zeta = loop {
            let h = f.from_u64(rng.gen_range(2..u64::MAX));
            let z = f.pow_biguint(h, &e);
            if !f.is_one(f.pow_u64(z, 1 << (t - 1))) {
                break z;
            }
        };
        let rules = (0..t)
            .map(|j| GSquareRule {
                next: (j + 1 < t).then_some((base + j + 1) as u8),
                a: f.zero(),
                b: if j + 1 < t { f.one() } else { f.zero() },
                c: f.zero(),
                d: if j + 1 < t { f.zero() } else { f.one() },
            })
            .collect();
        let points = (0..1u64 << t)
            .map(|k| {
                let mut y = vec![f.pow_u64(zeta, k)];
                for j in 1..t {
                    y.push(f.sqr(y[j - 1]));
                }
                y
            })
            .collect();
        (rules, points)
    }

    /// A random bidegree-(2, 2) polynomial in `y_{0,0}`, `y_{1,0}` (the
    /// shape of `S₃`), vanishing at `zero_at` when given.
    fn biquadratic<F: PrimeField>(
        ring: &GTowerRing<F>,
        t: usize,
        rng: &mut StdRng,
        zero_at: Option<&[F::Elem]>,
    ) -> GPoly<F::Elem> {
        let f = &ring.field;
        let n = ring.n_vars();
        let mut raw = Vec::new();
        for a in 0..3u32 {
            for b in 0..3u32 {
                let mut e = vec![0u32; n];
                e[0] = a;
                e[t] = b;
                raw.push((e, f.from_u64(rng.gen_range(1..u64::MAX))));
            }
        }
        let mut g = ring.from_raw(&raw);
        if let Some(pt) = zero_at {
            let v = ring.eval(&g, pt);
            raw.push((vec![0; n], f.neg(v)));
            g = ring.from_raw(&raw);
        }
        g
    }

    /// A random `u32`-engine ring: general rules or Kummer-like ones.
    fn random_ring(
        rng: &mut StdRng,
        p: u64,
        blocks: usize,
        t: usize,
        n_free: usize,
        general: bool,
    ) -> TowerRing {
        let mut rules = Vec::new();
        for b in 0..blocks {
            for j in 0..t {
                let top = j + 1 == t;
                rules.push(if general {
                    SquareRule {
                        next: (!top).then_some((b * t + j + 1) as u8),
                        a: if top { 0 } else { rng.gen_range(0..p) },
                        b: if top { 0 } else { rng.gen_range(0..p) },
                        c: rng.gen_range(0..p),
                        d: rng.gen_range(0..p),
                    }
                } else {
                    SquareRule {
                        next: (!top).then_some((b * t + j + 1) as u8),
                        a: 0,
                        b: u64::from(!top),
                        c: 0,
                        d: if top { rng.gen_range(1..p) } else { 0 },
                    }
                });
            }
        }
        TowerRing { p, rules, n_free }
    }

    /// Random sparse raw polynomials of degree at most 3.
    fn random_raw(rng: &mut StdRng, n: usize, p: u64, terms: usize) -> Vec<(Vec<u32>, u64)> {
        (0..terms)
            .map(|_| {
                let mut e = vec![0u32; n];
                let mut budget = 3u32;
                for x in e.iter_mut() {
                    if budget > 0 && rng.gen_range(0..10) < 3 {
                        *x = rng.gen_range(1..=budget.min(2));
                        budget -= *x;
                    }
                }
                (e, rng.gen_range(1..p))
            })
            .collect()
    }

    /// The generic engine over `f` against the `u32` engine on one system:
    /// the same report, field for field, and the same basis.
    fn same_as_u32<F: PrimeField>(f: F, ring: &TowerRing, input: &[RPoly], opts: &TowerF4Options) {
        let want = f4_tower(input, ring, opts);
        let gring = GTowerRing::from_tower_ring(f.clone(), ring);
        let ginput: Vec<GPoly<F::Elem>> = input.iter().map(|g| GPoly::from_rpoly(&f, g)).collect();
        let got = f4_tower_generic(&ginput, &gring, opts);
        assert_eq!(
            fingerprint(&got.stats),
            fingerprint(&want),
            "{} p = {}",
            f.name(),
            ring.p
        );
        let got_basis: Vec<RPoly> = got.basis.iter().map(|g| g.to_rpoly(&f).unwrap()).collect();
        assert_eq!(got_basis, want.basis, "{} p = {}", f.name(), ring.p);
    }

    /// Every backend that can hold the prime gives the `u32` engine's
    /// report and basis, on random rings and systems, below and above
    /// `2³¹`, with and without the full-rank exit.
    #[test]
    fn every_backend_reproduces_the_u32_engine() {
        let mut rng = StdRng::seed_from_u64(0x4634);
        for &p in &[97u64, 786_433, 2_013_265_921, 3_221_225_473, 4_294_967_291] {
            for trial in 0..8 {
                let general = trial % 2 == 1;
                let n_free = trial % 3 / 2;
                let t = 2 + trial % 3;
                let ring = random_ring(&mut rng, p, 2, t, n_free, general);
                let n = ring.n_vars();
                let input: Vec<RPoly> = (0..2 + trial % 2)
                    .map(|_| ring.from_raw(&random_raw(&mut rng, n, p, 6)))
                    .collect();
                let mut opts = TowerF4Options::new(n as u32 + 8);
                opts.full_rank_exit = trial % 4 != 3;
                same_as_u32(Fp32::new(p), &ring, &input, &opts);
                same_as_u32(Fp32::new_scalar(p), &ring, &input, &opts);
                same_as_u32(Fp64::new(p), &ring, &input, &opts);
                same_as_u32(FpWide::<1>::new(&BigUint::from(p)), &ring, &input, &opts);
                same_as_u32(FpWide::<4>::new(&BigUint::from(p)), &ring, &input, &opts);
                same_as_u32(Fp224::new(&BigUint::from(p)), &ring, &input, &opts);
                same_as_u32(Fp52::new(p), &ring, &input, &opts);
                same_as_u32(Fp52::new_scalar(p), &ring, &input, &opts);
                same_as_u32(
                    SmallCap {
                        inner: Fp32::new(p),
                        cap: 1,
                    },
                    &ring,
                    &input,
                    &opts,
                );
                same_as_u32(
                    SmallCap {
                        inner: Fp52::new(p),
                        cap: 3,
                    },
                    &ring,
                    &input,
                    &opts,
                );
            }
        }
    }

    /// A backend whose lazy capacity is `cap`: the engine must reduce and
    /// reload its rows that often and still give the same report.
    #[derive(Clone, Debug)]
    struct SmallCap<F> {
        inner: F,
        cap: u64,
    }

    impl<F: PrimeField> PrimeField for SmallCap<F> {
        type Elem = F::Elem;
        type Wide = F::Wide;
        fn name(&self) -> &'static str {
            "small-capacity"
        }
        fn modulus(&self) -> &BigUint {
            self.inner.modulus()
        }
        fn zero(&self) -> F::Elem {
            self.inner.zero()
        }
        fn one(&self) -> F::Elem {
            self.inner.one()
        }
        fn from_u64(&self, x: u64) -> F::Elem {
            self.inner.from_u64(x)
        }
        fn from_biguint(&self, x: &BigUint) -> F::Elem {
            self.inner.from_biguint(x)
        }
        fn to_biguint(&self, a: F::Elem) -> BigUint {
            self.inner.to_biguint(a)
        }
        fn add(&self, a: F::Elem, b: F::Elem) -> F::Elem {
            self.inner.add(a, b)
        }
        fn sub(&self, a: F::Elem, b: F::Elem) -> F::Elem {
            self.inner.sub(a, b)
        }
        fn neg(&self, a: F::Elem) -> F::Elem {
            self.inner.neg(a)
        }
        fn mul(&self, a: F::Elem, b: F::Elem) -> F::Elem {
            self.inner.mul(a, b)
        }
        fn inv(&self, a: F::Elem) -> Option<F::Elem> {
            self.inner.inv(a)
        }
        fn sqrt_pre(&self) -> &SqrtPre<F::Elem> {
            self.inner.sqrt_pre()
        }
        fn lazy_capacity(&self) -> u64 {
            self.cap
        }
        fn wide_zero(&self) -> F::Wide {
            self.inner.wide_zero()
        }
        fn wide_from(&self, a: F::Elem) -> F::Wide {
            self.inner.wide_from(a)
        }
        fn wide_mul_add(&self, w: &mut F::Wide, a: F::Elem, b: F::Elem) {
            self.inner.wide_mul_add(w, a, b)
        }
        fn wide_reduce(&self, w: &F::Wide) -> F::Elem {
            self.inner.wide_reduce(w)
        }
        fn row_axpy(&self, acc: &mut [F::Wide], c: F::Elem, row: &[F::Elem]) {
            self.inner.row_axpy(acc, c, row)
        }
        fn row_scatter_axpy(
            &self,
            acc: &mut [F::Wide],
            c: F::Elem,
            cols: &[u32],
            vals: &[F::Elem],
        ) {
            self.inner.row_scatter_axpy(acc, c, cols, vals)
        }
    }

    /// On Kummer systems large enough to use the blocked, chunked and
    /// parallel elimination paths, the same.
    #[test]
    fn a_larger_kummer_system_matches_the_u32_engine() {
        let mut rng = StdRng::seed_from_u64(77);
        for &p in &[2_013_265_921u64, 3_221_225_473] {
            let t = 6usize;
            let ring = random_ring(&mut rng, p, 2, t, 0, false);
            let n = ring.n_vars();
            let mut raw = Vec::new();
            for a in 0..3u32 {
                for b in 0..3u32 {
                    let mut e = vec![0u32; n];
                    e[0] = a;
                    e[t] = b;
                    raw.push((e, rng.gen_range(1..p)));
                }
            }
            let input = vec![ring.from_raw(&raw)];
            let opts = TowerF4Options::new(n as u32 + 8);
            same_as_u32(Fp32::new(p), &ring, &input, &opts);
            same_as_u32(Fp64::new(p), &ring, &input, &opts);
            same_as_u32(Fp224::new(&BigUint::from(p)), &ring, &input, &opts);
            same_as_u32(Fp52::new(p), &ring, &input, &opts);
        }
    }

    /// Refutes exactly the systems without a point on `V × V` and keeps
    /// every point of the rest, at primes no `u32` engine holds.
    fn refutes_exactly_at<F: PrimeField>(f: F, seed: u64) {
        let mut rng = StdRng::seed_from_u64(seed);
        for t in 1..=3usize {
            let (r0, pts) = kummer_block(&f, t, 0, &mut rng);
            let (r1, _) = kummer_block(&f, t, t, &mut rng);
            let ring = GTowerRing {
                field: f.clone(),
                rules: [r0, r1].concat(),
                n_free: 0,
            };
            for trial in 0..6 {
                let plant: Option<Vec<F::Elem>> = (trial % 2 == 0).then(|| {
                    let a = &pts[rng.gen_range(0..pts.len())];
                    let b = &pts[rng.gen_range(0..pts.len())];
                    [a.clone(), b.clone()].concat()
                });
                let g = biquadratic(&ring, t, &mut rng, plant.as_deref());
                let mut zeros = Vec::new();
                for a in &pts {
                    for b in &pts {
                        let x = [a.clone(), b.clone()].concat();
                        if f.is_zero(ring.eval(&g, &x)) {
                            zeros.push(x);
                        }
                    }
                }
                let rep =
                    f4_tower_generic(std::slice::from_ref(&g), &ring, &TowerF4Options::new(64));
                assert!(!rep.stats.timed_out && rep.stats.pairs_above_bound == 0);
                assert_eq!(
                    rep.stats.inconsistent,
                    zeros.is_empty(),
                    "{} t {t}",
                    f.name()
                );
                if !rep.stats.inconsistent {
                    for z in &zeros {
                        for b in &rep.basis {
                            assert!(f.is_zero(ring.eval(b, z)), "a solution was lost");
                        }
                    }
                    let lms: Vec<Mono> = rep.basis.iter().map(|b| b.lm().unwrap()).collect();
                    assert_eq!(
                        staircase_at_most(&lms, ring.n_tower(), 0, 1 << 12),
                        Some(zeros.len())
                    );
                }
            }
        }
    }

    /// One system over another backend (rules and input moved through
    /// canonical integers), fingerprinted with its basis.
    fn rerun<F: PrimeField, G: PrimeField>(
        ring: &GTowerRing<F>,
        g_field: &G,
        input: &[GPoly<F::Elem>],
    ) -> String {
        let f = &ring.field;
        let mv = |x: F::Elem| g_field.from_biguint(&f.to_biguint(x));
        let gring = GTowerRing {
            field: g_field.clone(),
            rules: ring
                .rules
                .iter()
                .map(|r| GSquareRule {
                    next: r.next,
                    a: mv(r.a),
                    b: mv(r.b),
                    c: mv(r.c),
                    d: mv(r.d),
                })
                .collect(),
            n_free: ring.n_free,
        };
        let ginput: Vec<GPoly<G::Elem>> = input
            .iter()
            .map(|g| GPoly {
                monos: g.monos.clone(),
                coefs: g.coefs.iter().map(|&c| mv(c)).collect(),
            })
            .collect();
        let rep = f4_tower_generic(&ginput, &gring, &TowerF4Options::new(64));
        format!(
            "{}|{:?}",
            fingerprint(&rep.stats),
            canonical(g_field, &rep.basis)
        )
    }

    #[test]
    fn refutes_exactly_at_a_40_bit_proth_prime() {
        let p = proth_prime(40);
        assert!(p > 1 << 40 && p < 1 << 63);
        refutes_exactly_at(Fp64::new(p), 40);
        refutes_exactly_at(Fp52::new(p), 40);
    }

    #[test]
    fn refutes_exactly_at_p224() {
        refutes_exactly_at(Fp224::new(&p224()), 224);
    }

    /// One system at a 56-bit prime and one at P-224, run under every
    /// backend that holds the prime: identical reports and bases.
    #[test]
    fn backends_agree_above_2_32() {
        let mut rng = StdRng::seed_from_u64(5656);
        let p = proth_prime(56);
        let f = Fp64::new(p);
        let t = 4;
        let (r0, _) = kummer_block(&f, t, 0, &mut rng);
        let (r1, _) = kummer_block(&f, t, t, &mut rng);
        let ring = GTowerRing {
            field: f.clone(),
            rules: [r0, r1].concat(),
            n_free: 0,
        };
        let g = biquadratic(&ring, t, &mut rng, None);
        let want = rerun(&ring, &f, std::slice::from_ref(&g));
        assert_eq!(
            rerun(
                &ring,
                &FpWide::<1>::new(&BigUint::from(p)),
                std::slice::from_ref(&g)
            ),
            want
        );
        assert_eq!(
            rerun(
                &ring,
                &FpWide::<4>::new(&BigUint::from(p)),
                std::slice::from_ref(&g)
            ),
            want
        );
        assert_eq!(
            rerun(
                &ring,
                &Fp224::new(&BigUint::from(p)),
                std::slice::from_ref(&g)
            ),
            want
        );
        for k in [40u32, 48] {
            let q = proth_prime(k);
            let f = Fp64::new(q);
            let (r0, _) = kummer_block(&f, t, 0, &mut rng);
            let (r1, _) = kummer_block(&f, t, t, &mut rng);
            let ring = GTowerRing {
                field: f.clone(),
                rules: [r0, r1].concat(),
                n_free: 0,
            };
            let g = biquadratic(&ring, t, &mut rng, None);
            let want = rerun(&ring, &f, std::slice::from_ref(&g));
            assert_eq!(rerun(&ring, &Fp52::new(q), std::slice::from_ref(&g)), want);
            assert_eq!(
                rerun(
                    &ring,
                    &SmallCap {
                        inner: Fp52::new(q),
                        cap: 2
                    },
                    std::slice::from_ref(&g)
                ),
                want
            );
            assert_eq!(
                rerun(
                    &ring,
                    &Fp224::new(&BigUint::from(q)),
                    std::slice::from_ref(&g)
                ),
                want
            );
        }

        let big = p224();
        let f = Fp224::new(&big);
        let (r0, _) = kummer_block(&f, t, 0, &mut rng);
        let (r1, _) = kummer_block(&f, t, t, &mut rng);
        let ring = GTowerRing {
            field: f.clone(),
            rules: [r0, r1].concat(),
            n_free: 0,
        };
        let g = biquadratic(&ring, t, &mut rng, None);
        let want = rerun(&ring, &f, std::slice::from_ref(&g));
        assert_eq!(
            rerun(&ring, &FpWide::<4>::new(&big), std::slice::from_ref(&g)),
            want
        );
    }

    /// `from_raw`, `mul_mono` and `eval` agree with the `u32` ring's.
    #[test]
    fn ring_operations_match_the_u32_ring() {
        let mut rng = StdRng::seed_from_u64(9);
        let p = 1_000_003u64;
        let f = Fp64::new(p);
        for trial in 0..40 {
            let ring = random_ring(&mut rng, p, 2, 3, 1, trial % 2 == 0);
            let gring = GTowerRing::from_tower_ring(f.clone(), &ring);
            let n = ring.n_vars();
            let raw = random_raw(&mut rng, n, p, 6);
            let want = ring.from_raw(&raw);
            let graw: Vec<(Vec<u32>, u64)> = raw.clone();
            let got = gring.from_raw(
                &graw
                    .iter()
                    .map(|(e, c)| (e.clone(), f.from_u64(*c)))
                    .collect::<Vec<_>>(),
            );
            assert_eq!(got.to_rpoly(&f).unwrap(), want);
            let u = mono(rng.gen_range(0..64), [rng.gen_range(0..2), 0, 0, 0]);
            assert_eq!(
                gring.mul_mono(&got, u).to_rpoly(&f).unwrap(),
                ring.mul_mono(&want, u)
            );
            let x: Vec<u64> = (0..n).map(|_| rng.gen_range(0..p)).collect();
            let gx: Vec<u64> = x.iter().map(|&v| f.from_u64(v)).collect();
            assert_eq!(
                f.to_u64(gring.eval(&got, &gx)).unwrap(),
                ring.eval(&want, &x)
            );
        }
    }
}
