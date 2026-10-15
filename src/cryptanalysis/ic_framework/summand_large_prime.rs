//! # Single and double large primes over any summand oracle.
//!
//! Gaudry's index calculus on `E(F_{qⁿ})` reaches `Õ(q^{2−2/n})` only with
//! the double large-prime variation (Gaudry 2009, after Gaudry–Thériault–
//! Thomé–Diem for hyperelliptic curves); without it the linear algebra on
//! a `q`-column matrix costs `q²`, the exponent of rho.  The device needs
//! nothing from the oracle: a decomposition into `n` points of the full
//! base `V` is kept whole when every summand lies in a smaller **small
//! base** `F ⊂ V`, kept as a **partial** when one (or, with `large_primes
//! = 2`, two) of its columns lie in `V \ F`, and partials are combined
//! until their large columns cancel.  The matrix then has `|F|` columns.
//!
//! This wrapper does that for any [`DecompositionOracle`] that returns
//! [`Relation::Summands`] (the `gaudry-d4` oracle on `E(F_{p⁴})` is the
//! one it was written for).  The oracle is unchanged and is asked the
//! same question on the same targets; only what happens to its answers
//! differs.
//!
//! ## Combining partials
//!
//! A relation is the row `Σ c_col·x_col = h·a + h·b·d` over `Z/rZ`; any
//! linear combination of rows is again a relation, with the combined
//! `(a, b)`, so a merged row holds as a group equation by construction
//! (the run's final check `[d]G = Q` catches any slip).  Partials are
//! reduced against stored **pivots on the large columns**: a pivot is a
//! partial normalised to coefficient 1 on its pivot column, with at most
//! one other large column (its *parent*).  A new partial is reduced while
//! one of its large columns has a pivot, which replaces that column by its
//! parent, so its large support never grows.  What is left is
//!
//! * no large column: a **full relation**, handed to the matrix;
//! * one or two large columns without pivots: a new pivot (two large
//!   columns join two trees of the large-prime graph; one large column
//!   ties its tree to the small base).
//!
//! These are the cycles of the classical large-prime graph, found by
//! union-by-elimination; with coefficients `±1` a cycle that does not
//! pass through the small base can leave one large column with a non-zero
//! coefficient, which is then exactly a single-large partial on that
//! column, as the elimination treats it.
//!
//! ## Counting
//!
//! Each multiplication modulo `r` in a reduction is one `row_op`, the
//! unit the matrix counts in; the total is reported in the relation phase
//! both as `lp_merge_row_ops` (descriptive) and as `row_ops`, which
//! `price_phase` prices.  So merging is priced in the relation phase, and
//! the matrix phase holds only the matrix.

use std::collections::{BTreeMap, HashMap};

use super::stages::{DecompositionOracle, InstanceCtx, Params, SolverTotals, SystemShape};
use crate::cryptanalysis::ic_boundary::{
    CountedGroup, FactorBase, GroupOps, OracleCounters, Relation,
};

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}
fn submod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 + m as u128 - b as u128) % m as u128) as u64
}
fn invmod(a: u64, m: u64) -> Option<u64> {
    let (mut r0, mut r1) = (m as i128, (a % m) as i128);
    let (mut t0, mut t1) = (0i128, 1i128);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (t0, t1) = (t1, t0 - q * t1);
    }
    (r0 == 1).then(|| ((t0 % m as i128 + m as i128) % m as i128) as u64)
}

/// A relation row in progress: sparse columns (sorted), and `(a, b)`.
#[derive(Clone, Debug)]
struct Row {
    cols: Vec<(usize, u64)>,
    a: u64,
    b: u64,
}

#[derive(Clone, Copy, Debug, Default)]
struct Counters {
    decompositions: u64,
    full_direct: u64,
    partials_one: u64,
    partials_two: u64,
    merged_full: u64,
    discarded: u64,
    trivial: u64,
    merge_row_ops: u64,
    reductions: u64,
    max_depth: u64,
}

/// Which columns of the base are small, and how many large primes a
/// partial may carry.
#[derive(Clone, Debug)]
pub struct SmallBase {
    pub small: Vec<bool>,
    pub large_primes: u32,
}

impl SmallBase {
    /// `large_primes` (1 or 2), and the small columns: the zero-log
    /// columns listed in `always_small`, plus the first
    /// `⌈scale·N^exponent⌉` of the `N` columns not listed there (capped
    /// at `N`, at least 1).
    pub fn first_columns(
        columns: usize,
        always_small: &[usize],
        large_primes: u32,
        exponent: f64,
        scale: f64,
    ) -> Result<Self, String> {
        if !(1..=2).contains(&large_primes) {
            return Err("large_primes is 1 or 2".into());
        }
        if !(exponent > 0.0 && exponent <= 1.0) || !(scale > 0.0) {
            return Err("lp_small_exp is in (0, 1] and lp_small_scale is positive".into());
        }
        let mut small = vec![false; columns];
        for &c in always_small {
            small[c] = true;
        }
        let n = columns - always_small.len();
        let k = ((scale * (n as f64).powf(exponent)).ceil() as usize).clamp(1, n.max(1));
        let mut taken = 0;
        for c in 0..columns {
            if taken == k {
                break;
            }
            if !small[c] {
                small[c] = true;
                taken += 1;
            }
        }
        Ok(Self {
            small,
            large_primes,
        })
    }

    pub fn small_count(&self) -> usize {
        self.small.iter().filter(|&&s| s).count()
    }
}

/// The elimination on the large columns, independent of the group: rows
/// in, full relations out.
pub struct LargePrimeForest {
    base: SmallBase,
    r: u64,
    /// Large column → its pivot row (coefficient 1 on that column).
    pivots: HashMap<usize, Row>,
    counters: Counters,
}

impl LargePrimeForest {
    pub fn new(base: SmallBase, r: u64) -> Self {
        Self {
            base,
            r,
            pivots: HashMap::new(),
            counters: Counters::default(),
        }
    }

    /// One decomposition's row (`Σ coef·x_col = h·a + h·b·d`, columns
    /// sorted, no zeros): a full relation now, or nothing (stored as a
    /// partial, discarded for too many large columns, or trivial).
    pub fn ingest(&mut self, cols: Vec<(usize, u64)>, a: u64, b: u64) -> Option<Relation> {
        self.counters.decompositions += 1;
        let m = self.r;
        let row = Row {
            cols,
            a: a % m,
            b: b % m,
        };
        match self.large_cols(&row).len() {
            0 => {}
            1 => self.counters.partials_one += 1,
            2 if self.base.large_primes == 2 => self.counters.partials_two += 1,
            _ => {
                self.counters.discarded += 1;
                return None;
            }
        }
        self.absorb(row)
    }

    fn counter_list(&self) -> Vec<(&'static str, u64)> {
        let c = &self.counters;
        vec![
            ("lp_decompositions", c.decompositions),
            ("lp_full_direct", c.full_direct),
            ("lp_partials_one_large", c.partials_one),
            ("lp_partials_two_large", c.partials_two),
            ("lp_merged_full", c.merged_full),
            ("lp_discarded", c.discarded),
            ("lp_trivial", c.trivial),
            ("lp_pivots", self.pivots.len() as u64),
            ("lp_reductions", c.reductions),
            ("lp_max_merge_depth", c.max_depth),
            ("lp_small_columns", self.base.small_count() as u64),
            ("lp_merge_row_ops", c.merge_row_ops),
            ("row_ops", c.merge_row_ops),
        ]
    }

    fn large_cols(&self, row: &Row) -> Vec<usize> {
        row.cols
            .iter()
            .filter(|(c, _)| !self.base.small[*c])
            .map(|&(c, _)| c)
            .collect()
    }

    /// `row − f·pivot`, counting one `row_op` per multiplication.
    fn sub_scaled(&mut self, row: &Row, f: u64, pivot: &Row) -> Row {
        let m = self.r;
        let mut out = Vec::with_capacity(row.cols.len() + pivot.cols.len());
        let (mut i, mut j) = (0, 0);
        let mut ops = 0u64;
        while i < row.cols.len() || j < pivot.cols.len() {
            match (row.cols.get(i), pivot.cols.get(j)) {
                (Some(&(ca, va)), Some(&(cb, vb))) if ca == cb => {
                    let v = submod(va, mulmod(f, vb, m), m);
                    ops += 1;
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
                    ops += 1;
                    if v != 0 {
                        out.push((cb, v));
                    }
                    j += 1;
                }
                (None, None) => break,
            }
        }
        let a = submod(row.a, mulmod(f, pivot.a, m), m);
        let b = submod(row.b, mulmod(f, pivot.b, m), m);
        ops += 2;
        self.counters.merge_row_ops += ops;
        Row { cols: out, a, b }
    }

    /// Reduce, then emit a full relation or store a pivot.
    fn absorb(&mut self, mut row: Row) -> Option<Relation> {
        let mut depth = 0u64;
        loop {
            let hit = row
                .cols
                .iter()
                .find(|(c, _)| !self.base.small[*c] && self.pivots.contains_key(c))
                .copied();
            let Some((c, f)) = hit else { break };
            let pivot = self.pivots[&c].clone();
            row = self.sub_scaled(&row, f, &pivot);
            depth += 1;
        }
        self.counters.reductions += depth;
        self.counters.max_depth = self.counters.max_depth.max(depth);
        let large = self.large_cols(&row);
        if large.is_empty() {
            if row.cols.is_empty() && row.b == 0 {
                // 0 = h·a: nothing learnt (a ≠ 0 cannot happen for a
                // true combination of true relations).
                self.counters.trivial += 1;
                return None;
            }
            if depth > 0 {
                self.counters.merged_full += 1;
            } else {
                self.counters.full_direct += 1;
            }
            return Some(Relation::Row {
                columns: row.cols,
                a: row.a,
                b: row.b,
            });
        }
        // One or two large columns, none pivoted: a new pivot on the first.
        let pc = large[0];
        let pv = row
            .cols
            .iter()
            .find(|(c, _)| *c == pc)
            .map(|&(_, v)| v)
            .expect("pivot column in row");
        let inv = invmod(pv, self.r).expect("prime modulus");
        let m = self.r;
        for (_, v) in row.cols.iter_mut() {
            *v = mulmod(*v, inv, m);
        }
        row.a = mulmod(row.a, inv, m);
        row.b = mulmod(row.b, inv, m);
        self.counters.merge_row_ops += row.cols.len() as u64 + 2;
        self.pivots.insert(pc, row);
        None
    }
}

/// The wrapper: the inner oracle and the forest.
pub struct SummandLargePrime<'o, G: CountedGroup> {
    inner: &'o mut dyn DecompositionOracle<G>,
    forest: LargePrimeForest,
}

impl<'o, G: CountedGroup> SummandLargePrime<'o, G> {
    pub fn new(inner: &'o mut dyn DecompositionOracle<G>, base: SmallBase, r: u64) -> Self {
        Self {
            inner,
            forest: LargePrimeForest::new(base, r),
        }
    }
}

impl<G: CountedGroup> DecompositionOracle<G> for SummandLargePrime<'_, G> {
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn summands(&self) -> u32 {
        self.inner.summands()
    }
    fn describe(&self, params: &Params) -> String {
        format!(
            "{}; {} large prime{} over a small base of {} of {} columns",
            self.inner.describe(params),
            self.forest.base.large_primes,
            if self.forest.base.large_primes == 1 {
                ""
            } else {
                "s"
            },
            self.forest.base.small_count(),
            self.forest.base.small.len()
        )
    }
    fn prepare(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        if self.forest.base.small.len() != fb.columns {
            return Err(format!(
                "small base covers {} columns, the base has {}",
                self.forest.base.small.len(),
                fb.columns
            ));
        }
        self.inner.prepare(ctx, fb, params, ops)
    }
    fn setup_native(&self) -> BTreeMap<String, u64> {
        self.inner.setup_native()
    }
    fn decompose(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: G::Elt,
    ) -> Option<Vec<usize>> {
        self.inner.decompose(ctx, fb, ops, counters, point)
    }
    fn decompose_relations(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: G::Elt,
        a: u64,
        b: u64,
    ) -> Vec<Relation> {
        let Some(summands) = self.inner.decompose(ctx, fb, ops, counters, point) else {
            return Vec::new();
        };
        let m = self.forest.r;
        let mut acc: BTreeMap<usize, u64> = BTreeMap::new();
        for &i in &summands {
            let e = acc.entry(fb.col_of[i]).or_insert(0);
            *e = ((*e as u128 + (fb.coef_of[i] % m) as u128) % m as u128) as u64;
        }
        let cols = acc.into_iter().filter(|&(_, v)| v != 0).collect();
        self.forest.ingest(cols, a, b).into_iter().collect()
    }
    fn relation_counters(&self) -> BTreeMap<String, u64> {
        let mut out = self.inner.relation_counters();
        for (k, v) in self.forest.counter_list() {
            *out.entry(k.into()).or_insert(0) += v;
        }
        out
    }
    fn last_system(&self) -> Option<SystemShape> {
        self.inner.last_system()
    }
    fn solver_totals(&self) -> Option<SolverTotals> {
        self.inner.solver_totals()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    /// Every relation the forest emits has no large column and holds for
    /// the planted logarithms; partials do merge into relations.
    #[test]
    fn merged_relations_hold_and_carry_no_large_column() {
        let r = 1_000_003u64;
        for large_primes in [1u32, 2] {
            let cols = 40usize;
            let base =
                SmallBase::first_columns(cols, &[cols - 1], large_primes, 0.75, 1.0).unwrap();
            assert_eq!(base.small_count(), 1 + (39f64.powf(0.75).ceil() as usize));
            let small = base.small.clone();
            let mut forest = LargePrimeForest::new(base, r);
            let mut rng = StdRng::seed_from_u64(7 + large_primes as u64);
            let mut x: Vec<u64> = (0..cols).map(|_| rng.gen_range(0..r)).collect();
            x[cols - 1] = 0;
            let d = rng.gen_range(0..r);
            let (mut emitted, mut merged) = (0, 0);
            for _ in 0..2000 {
                let mut acc: BTreeMap<usize, u64> = BTreeMap::new();
                for _ in 0..4 {
                    let c = rng.gen_range(0..cols);
                    let s = if rng.gen_bool(0.5) { 1 } else { r - 1 };
                    let e = acc.entry(c).or_insert(0);
                    *e = (*e + s) % r;
                }
                let row: Vec<(usize, u64)> = acc.into_iter().filter(|&(_, v)| v != 0).collect();
                let lhs = row
                    .iter()
                    .fold(0u64, |s, &(c, v)| (s + mulmod(v, x[c], r)) % r);
                let b = rng.gen_range(1..r);
                let a = submod(lhs, mulmod(b, d, r), r);
                let before = forest.counters.merged_full;
                if let Some(Relation::Row { columns, a, b }) = forest.ingest(row, a, b) {
                    emitted += 1;
                    merged += (forest.counters.merged_full > before) as u32;
                    assert!(columns.iter().all(|&(c, _)| small[c]));
                    let lhs = columns
                        .iter()
                        .fold(0u64, |s, &(c, v)| (s + mulmod(v, x[c], r)) % r);
                    assert_eq!(lhs, (a + mulmod(b, d, r)) % r);
                }
            }
            assert!(
                emitted > 20 && merged > 5,
                "{large_primes}: {emitted} {merged}"
            );
        }
    }
}
