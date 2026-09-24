//! **Decomposition systems past 64 unknowns** — plan target C of
//! `docs/ic/perf/OPTIMIZATION_PLAN.md`.
//!
//! Every `F_2` solver in [`crate::cryptanalysis::koblitz_groebner`] keeps a
//! monomial in a `u64`, so the chained `S₃` system of an `m`-point
//! decomposition — `m·ℓ + (m − 2)·n` unknowns — cannot be written down
//! once that exceeds 64: already at `m = 4`, `n = 31`, `ℓ = 6` (86
//! unknowns), and at every `m ≥ 4` rung the index-calculus family cares
//! about.  This module is the same system in `u128` monomials (up to 128
//! unknowns) and a solver for it, kept separate so the 64-bit engines,
//! their counters and their frozen benchmarks are untouched.
//!
//! ## The solver
//!
//! A root of the system is wanted only for its summand abscissae: once
//! those `m·ℓ` bits are fixed, [`lift_candidate`] settles in the group
//! whether signs exist that make the summands add to the target.  So the
//! search branches on summand bits only, lowest first — which fixes one
//! summand at a time — and never on the `(m − 2)·n` intermediate
//! unknowns.  At every node the equations, with the branch's assignment
//! substituted, are **linearised**: their Macaulay matrix at their own
//! degree, columns in descending degree, reduced to row echelon form by
//! [`crate::cryptanalysis::gf2_elim`].  A reduced row equal to `1`
//! closes the branch; a row of degree `≤ 1` is a linear consequence,
//! eliminated by substitution, and the step repeats until nothing new
//! falls out.  That is the degree-`d` Macaulay step of the 64-bit
//! engines' default `max_degree = 3` on cubic systems, without the
//! splitting on intermediate unknowns that their measured failure at
//! `n = 31, m = 3` traced to.
//!
//! What this does not claim: that it is fast.  Its purpose is to make
//! the `m ≥ 4` cells of `examples/pdp_bench.rs` *measurable*, so the
//! Gröbner route's frontier at `m = 4, 5, 6` is a table rather than a
//! "not buildable".

use std::collections::HashMap;

use num_bigint::BigUint;

use crate::binary_ecc::{BinaryPoint, F2mElement};
use crate::cryptanalysis::fx_hash::FxMap;
use crate::cryptanalysis::gf2_elim;
use crate::cryptanalysis::koblitz_groebner::FieldStructure;
use crate::cryptanalysis::koblitz_index_calculus::{
    lift_candidate, point_key, points_with_x, FrobeniusFactorBase, KoblitzCurve,
};

/// Most unknowns a [`WPoly`] can carry.
pub const MAX_WIDE_VARS: usize = 128;

/// A monomial: bit `v` set means variable `v` divides it.
pub type Mono = u128;

/// A polynomial over `F_2` in Boolean variables (`v² = v`): its
/// monomials, ascending as integers, no repeats.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct WPoly {
    pub terms: Vec<Mono>,
}

impl WPoly {
    pub fn zero() -> Self {
        Self { terms: Vec::new() }
    }

    pub fn one() -> Self {
        Self { terms: vec![0] }
    }

    pub fn var(v: usize) -> Self {
        Self {
            terms: vec![1u128 << v],
        }
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn is_one(&self) -> bool {
        self.terms == [0]
    }

    /// Highest monomial degree (0 for constants and zero).
    pub fn degree(&self) -> u32 {
        self.terms.iter().map(|m| m.count_ones()).max().unwrap_or(0)
    }

    /// From any list of monomials: sorted, and pairs cancelled.
    pub fn from_monos(mut monos: Vec<Mono>) -> Self {
        monos.sort_unstable();
        let mut terms = Vec::with_capacity(monos.len());
        let mut i = 0;
        while i < monos.len() {
            let mut j = i;
            while j < monos.len() && monos[j] == monos[i] {
                j += 1;
            }
            if (j - i) % 2 == 1 {
                terms.push(monos[i]);
            }
            i = j;
        }
        Self { terms }
    }

    /// Sum: the symmetric difference of the term sets.
    pub fn add(&self, other: &Self) -> Self {
        let (a, b) = (&self.terms, &other.terms);
        let mut out = Vec::with_capacity(a.len() + b.len());
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            match a[i].cmp(&b[j]) {
                std::cmp::Ordering::Less => {
                    out.push(a[i]);
                    i += 1;
                }
                std::cmp::Ordering::Greater => {
                    out.push(b[j]);
                    j += 1;
                }
                std::cmp::Ordering::Equal => {
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend_from_slice(&a[i..]);
        out.extend_from_slice(&b[j..]);
        Self { terms: out }
    }

    /// Product in the Boolean quotient: monomials multiply by union.
    pub fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let mut monos = Vec::with_capacity(self.terms.len() * other.terms.len());
        for &a in &self.terms {
            for &b in &other.terms {
                monos.push(a | b);
            }
        }
        Self::from_monos(monos)
    }

    /// Substitute the affine form `form` for variable `v`.
    pub fn substitute(&self, v: usize, form: &WPoly) -> Self {
        let bit = 1u128 << v;
        if self.terms.iter().all(|m| m & bit == 0) {
            return self.clone();
        }
        let (with, without): (Vec<Mono>, Vec<Mono>) =
            self.terms.iter().partition(|&&m| m & bit != 0);
        let cofactor = Self::from_monos(with.into_iter().map(|m| m & !bit).collect());
        Self { terms: without }.add(&cofactor.mul(form))
    }

    /// Evaluate at a full assignment (`bits` bit `v` = value of variable `v`).
    pub fn eval(&self, bits: u128) -> bool {
        self.terms.iter().filter(|&&m| m & !bits == 0).count() % 2 == 1
    }
}

/// A field element whose `n` coordinates are [`WPoly`]s.
#[derive(Clone, Debug)]
pub struct WSym {
    pub coords: Vec<WPoly>,
}

impl WSym {
    fn constant(v: &F2mElement, n: u32) -> Self {
        let bits = v.raw_bits().first().copied().unwrap_or(0);
        Self {
            coords: (0..n)
                .map(|k| {
                    if bits >> k & 1 == 1 {
                        WPoly::one()
                    } else {
                        WPoly::zero()
                    }
                })
                .collect(),
        }
    }

    fn from_subspace_vars(basis: &[F2mElement], offset: usize, n: u32) -> Self {
        let mut coords = vec![Vec::new(); n as usize];
        for (t, b) in basis.iter().enumerate() {
            let bits = b.raw_bits().first().copied().unwrap_or(0);
            for (k, c) in coords.iter_mut().enumerate() {
                if bits >> k & 1 == 1 {
                    c.push(1u128 << (offset + t));
                }
            }
        }
        Self {
            coords: coords.into_iter().map(WPoly::from_monos).collect(),
        }
    }

    fn from_free_vars(offset: usize, n: u32) -> Self {
        Self {
            coords: (0..n as usize).map(|k| WPoly::var(offset + k)).collect(),
        }
    }

    fn add(&self, other: &Self) -> Self {
        Self {
            coords: self
                .coords
                .iter()
                .zip(&other.coords)
                .map(|(a, b)| a.add(b))
                .collect(),
        }
    }

    fn mul(&self, other: &Self, st: &FieldStructure) -> Self {
        let n = st.n as usize;
        let mut acc: Vec<Vec<Mono>> = vec![Vec::new(); n];
        for (i, a) in self.coords.iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for (j, b) in other.coords.iter().enumerate() {
                if b.is_zero() {
                    continue;
                }
                let prod = a.mul(b);
                let red = st.reduced[i][j];
                for (k, slot) in acc.iter_mut().enumerate() {
                    if red >> k & 1 == 1 {
                        slot.extend_from_slice(&prod.terms);
                    }
                }
            }
        }
        Self {
            coords: acc.into_iter().map(WPoly::from_monos).collect(),
        }
    }

    fn square(&self, st: &FieldStructure) -> Self {
        let n = st.n as usize;
        let mut acc: Vec<Vec<Mono>> = vec![Vec::new(); n];
        for (k, c) in self.coords.iter().enumerate() {
            let sq = st.squares[k];
            for (t, slot) in acc.iter_mut().enumerate() {
                if sq >> t & 1 == 1 {
                    slot.extend_from_slice(&c.terms);
                }
            }
        }
        Self {
            coords: acc.into_iter().map(WPoly::from_monos).collect(),
        }
    }
}

/// `S₃(x₁, x₂, x₃) = (x₁+x₂)² x₃² + x₁x₂x₃ + (x₁x₂)² + b`, coordinatewise.
fn s3(x1: &WSym, x2: &WSym, x3: &WSym, b: &F2mElement, st: &FieldStructure) -> Vec<WPoly> {
    let prod = x1.mul(x2, st);
    x1.add(x2)
        .square(st)
        .mul(&x3.square(st), st)
        .add(&prod.mul(x3, st))
        .add(&prod.square(st))
        .add(&WSym::constant(b, st.n))
        .coords
}

/// The chained `S₃` system of an `m`-point decomposition of a target with
/// abscissa `x_r`, summands in the span of `basis` — the layout of
/// `koblitz_groebner::build_decomposition_system`: summand `i` owns
/// variables `[i·ℓ, (i+1)·ℓ)`, then `(m − 2)·n` intermediate unknowns.
#[derive(Clone, Debug)]
pub struct WideSystem {
    pub equations: Vec<WPoly>,
    pub n_vars: usize,
    pub ell: usize,
    pub m: usize,
}

impl WideSystem {
    pub fn build(
        basis: &[F2mElement],
        x_r: &F2mElement,
        b: &F2mElement,
        m: usize,
        st: &FieldStructure,
    ) -> Option<Self> {
        if m < 2 || st.n > 64 {
            return None;
        }
        let n = st.n;
        let ell = basis.len();
        let n_vars = m * ell + (m - 2) * n as usize;
        if n_vars > MAX_WIDE_VARS {
            return None;
        }
        let xs: Vec<WSym> = (0..m)
            .map(|i| WSym::from_subspace_vars(basis, i * ell, n))
            .collect();
        let inter: Vec<WSym> = (0..m - 2)
            .map(|i| WSym::from_free_vars(m * ell + i * n as usize, n))
            .collect();
        let target = WSym::constant(x_r, n);
        let mut equations = Vec::new();
        if m == 2 {
            equations.extend(s3(&xs[0], &xs[1], &target, b, st));
        } else {
            equations.extend(s3(&xs[0], &xs[1], &inter[0], b, st));
            for i in 0..m - 3 {
                equations.extend(s3(&inter[i], &xs[i + 2], &inter[i + 1], b, st));
            }
            equations.extend(s3(&inter[m - 3], &xs[m - 1], &target, b, st));
        }
        equations.retain(|p| !p.is_zero());
        Some(Self {
            equations,
            n_vars,
            ell,
            m,
        })
    }

    /// The last `k` links of an `m`-point chain, `m > k`, starting from a
    /// **free** abscissa: `w_0` stands for `x(P_1 + … + P_{m−k})`, and the
    /// links are `S₃(w_0, y_1, w_1), …, S₃(w_{k−1}, y_k, x_R)`.  Summand
    /// `y_i` owns variables `[(i−1)·ℓ, i·ℓ)`, then the `k` intermediates
    /// `w_0 … w_{k−1}`: `k·(ℓ + n)` unknowns.  Every decomposition of the
    /// target satisfies it, whatever the first `m − k` summands are, so
    /// it prunes soundly; it cannot finish a decomposition on its own.
    pub fn build_suffix(
        basis: &[F2mElement],
        x_r: &F2mElement,
        b: &F2mElement,
        k: usize,
        st: &FieldStructure,
    ) -> Option<Self> {
        if k == 0 || st.n > 64 {
            return None;
        }
        let n = st.n;
        let ell = basis.len();
        let n_vars = k * (ell + n as usize);
        if n_vars > MAX_WIDE_VARS {
            return None;
        }
        let ys: Vec<WSym> = (0..k)
            .map(|i| WSym::from_subspace_vars(basis, i * ell, n))
            .collect();
        let ws: Vec<WSym> = (0..k)
            .map(|i| WSym::from_free_vars(k * ell + i * n as usize, n))
            .collect();
        let target = WSym::constant(x_r, n);
        let mut equations = Vec::new();
        for i in 0..k {
            let next = if i + 1 < k { &ws[i + 1] } else { &target };
            equations.extend(s3(&ws[i], &ys[i], next, b, st));
        }
        equations.retain(|p| !p.is_zero());
        Some(Self {
            equations,
            n_vars,
            ell,
            m: k,
        })
    }
}

/// What a wide solve did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct WideStats {
    /// Search nodes visited.
    pub nodes: usize,
    /// Linearisation passes (one Macaulay reduction each).
    pub reductions: usize,
    /// Branches closed by a reduced `1`.
    pub refuted: usize,
    /// Leaves whose summands were lifted in the group.
    pub leaves: usize,
    /// Largest linearisation matrix, rows × columns.
    pub max_rows: usize,
    pub max_cols: usize,
    /// The node budget ran out: a `None` then says nothing.
    pub exhausted: bool,
}

/// Columns above which a node skips linearisation and just splits: a
/// dense matrix that wide would cost more than the branch it might close.
const MAX_LINEARISATION_COLS: usize = 1 << 20;

/// Linearise `eqs` to a fixpoint: reduce their Macaulay matrix at their
/// own degree, substitute every linear consequence, repeat.  Returns
/// `false` when a reduction yields `1` (the branch has no root), and
/// otherwise the substitutions made (variable, affine form) in order.
fn linearise(eqs: &mut Vec<WPoly>, stats: &mut WideStats, subs: &mut Vec<(usize, WPoly)>) -> bool {
    loop {
        eqs.retain(|p| !p.is_zero());
        if eqs.iter().any(WPoly::is_one) {
            return false;
        }
        if eqs.is_empty() {
            return true;
        }
        // Columns: every monomial, descending degree, then descending
        // integer (any fixed order within a degree serves), so a reduced
        // row whose leading column has degree ≤ 1 is linear.
        //
        // Sorted as plain integers to deduplicate, then laid out by degree
        // with a counting pass: a comparator that recounts bits on every
        // compare was 40 % of a search.
        let mut all: Vec<Mono> = eqs.iter().flat_map(|p| p.terms.iter().copied()).collect();
        all.sort_unstable();
        all.dedup();
        if all.len() > MAX_LINEARISATION_COLS {
            return true;
        }
        let max_deg = all.iter().map(|m| m.count_ones()).max().unwrap_or(0) as usize;
        let mut start = vec![0usize; max_deg + 2];
        for m in &all {
            // Descending degree: degree d goes after every higher one.
            start[max_deg - m.count_ones() as usize + 1] += 1;
        }
        for d in 1..start.len() {
            start[d] += start[d - 1];
        }
        let mut cols = vec![0 as Mono; all.len()];
        let mut index: FxMap<Mono, usize> = FxMap::default();
        index.reserve(all.len());
        for &m in all.iter().rev() {
            let slot = &mut start[max_deg - m.count_ones() as usize];
            cols[*slot] = m;
            index.insert(m, *slot);
            *slot += 1;
        }
        let words = cols.len().div_ceil(64);
        let mut matrix: Vec<Vec<u64>> = eqs
            .iter()
            .map(|p| {
                let mut row = vec![0u64; words];
                for m in &p.terms {
                    let c = index[m];
                    row[c / 64] |= 1 << (c % 64);
                }
                row
            })
            .collect();
        stats.reductions += 1;
        stats.max_rows = stats.max_rows.max(matrix.len());
        stats.max_cols = stats.max_cols.max(cols.len());
        let mut ops = 0;
        let rank = gf2_elim::eliminate(
            &mut matrix,
            cols.len(),
            true,
            gf2_elim::Config::from_env(),
            &mut ops,
        );
        let to_poly = |row: &Vec<u64>| {
            let mut terms: Vec<Mono> = Vec::new();
            for (w, &word) in row.iter().enumerate() {
                let mut bits = word;
                while bits != 0 {
                    let b = bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    terms.push(cols[w * 64 + b]);
                }
            }
            WPoly::from_monos(terms)
        };
        let reduced: Vec<WPoly> = matrix[..rank].iter().map(to_poly).collect();
        if reduced.iter().any(WPoly::is_one) {
            return false;
        }
        // Linear rows: leading (highest) monomial of degree 1.  In the
        // reduced form each has a distinct leading variable, so they
        // can be eliminated one after another.
        let mut linear: Vec<(usize, WPoly)> = Vec::new();
        for p in &reduced {
            if p.degree() == 1 {
                let lead = p
                    .terms
                    .iter()
                    .filter(|m| m.count_ones() == 1)
                    .max()
                    .copied()
                    .expect("degree one");
                let v = lead.trailing_zeros() as usize;
                let form = p.add(&WPoly::var(v)); // v = rest
                linear.push((v, form));
            }
        }
        if linear.is_empty() {
            *eqs = reduced;
            return true;
        }
        *eqs = reduced;
        for (v, form) in linear {
            // Earlier substitutions may already mention `v` in `form`'s
            // variables' place; apply them to the form first.
            let mut f = form;
            for (u, g) in subs.iter() {
                f = f.substitute(*u, g);
            }
            if f.terms.iter().any(|m| m >> v & 1 == 1) {
                // `v = … + v …` cannot happen for a linear form with `v`
                // removed; guard anyway.
                continue;
            }
            for p in eqs.iter_mut() {
                *p = p.substitute(v, &f);
            }
            for (_, g) in subs.iter_mut() {
                *g = g.substitute(v, &f);
            }
            subs.push((v, f));
        }
    }
}

/// Decompose `target` into `m` factor-base points through the wide
/// chained system; `None` with `stats.exhausted` unset is a complete
/// search that found nothing, with it set says nothing.
///
/// When the whole chain needs more than [`MAX_WIDE_VARS`] unknowns, the
/// search takes the largest `k` trailing summands whose chain suffix
/// ([`WideSystem::build_suffix`]) fits, branches on those with the suffix
/// as the pruning system, and at each of its leaves decomposes the
/// remainder `target − ΣP_i` into the other `m − k` summands, by the
/// same function.
#[allow(clippy::too_many_arguments)]
pub fn wide_groebner_decompose(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    st: &FieldStructure,
    target: &BinaryPoint,
    m: usize,
    node_budget: usize,
) -> (Option<Vec<usize>>, WideStats) {
    let mut stats = WideStats::default();
    let opts = SearchOptions::from_env();
    let found = decompose_rec(
        kc,
        fb,
        index_of,
        st,
        target,
        m,
        node_budget,
        &opts,
        &mut stats,
    );
    (found, stats)
}

/// The search's switches, each a same-binary control read once per call.
#[derive(Clone, Copy)]
struct SearchOptions {
    /// Only summands in non-decreasing order (`KIC_WIDE_ORDER=0` off).
    order: bool,
    /// Each summand from its highest coordinate (`KIC_WIDE_BRANCH=low` off).
    high_first: bool,
    /// Look the last summand up instead of branching (`KIC_WIDE_FINISH=0` off).
    finish: bool,
    /// Widest system built in one piece; above it the suffix recursion
    /// takes over.  [`MAX_WIDE_VARS`] except in tests, which lower it to
    /// reach the recursion on small curves.
    max_vars: usize,
}

impl SearchOptions {
    fn from_env() -> Self {
        Self {
            order: std::env::var("KIC_WIDE_ORDER").as_deref() != Ok("0"),
            high_first: std::env::var("KIC_WIDE_BRANCH").as_deref() != Ok("low"),
            finish: std::env::var("KIC_WIDE_FINISH").as_deref() != Ok("0"),
            max_vars: MAX_WIDE_VARS,
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn decompose_rec(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    st: &FieldStructure,
    target: &BinaryPoint,
    m: usize,
    node_budget: usize,
    opts: &SearchOptions,
    stats: &mut WideStats,
) -> Option<Vec<usize>> {
    if m == 0 {
        return (*target == BinaryPoint::Infinity).then(Vec::new);
    }
    if m == 1 {
        return index_of.get(&point_key(target)).map(|&i| vec![i]);
    }
    let x_r = match target {
        BinaryPoint::Affine { x, .. } => x.clone(),
        BinaryPoint::Infinity => return None,
    };
    let ell = fb.subspace_basis.len();
    let full_vars = m * ell + (m - 2) * st.n as usize;
    let whole = if full_vars <= opts.max_vars {
        WideSystem::build(&fb.subspace_basis, &x_r, &kc.curve.b, m, st)
    } else {
        None
    };
    if let Some(sys) = whole {
        // The whole chain fits: its leaves finish the decomposition.
        let ell = sys.ell;
        let finish_head = opts.finish.then_some(m - 1);
        return search(
            &sys,
            m,
            finish_head,
            node_budget,
            opts,
            stats,
            &mut |values, is_head, _stats| {
                let count = if is_head { m - 1 } else { m };
                let xs: Vec<F2mElement> = (0..count)
                    .map(|i| summand_x(&fb.subspace_basis, ell, values, i, kc.n))
                    .collect();
                if is_head {
                    finish_last(kc, fb, index_of, &xs, target)
                } else {
                    lift_candidate(kc, fb, index_of, &xs, target)
                }
            },
        );
    }
    // Too wide: prune with the longest chain suffix that fits.
    let k = (1..m - 1)
        .rev()
        .find(|&k| k * (ell + st.n as usize) <= opts.max_vars)?;
    let sys = WideSystem::build_suffix(&fb.subspace_basis, &x_r, &kc.curve.b, k, st)?;
    search(
        &sys,
        k,
        None,
        node_budget,
        opts,
        stats,
        &mut |values, _is_head, stats| {
            // The k trailing summands are known up to sign: for each
            // choice of base points, decompose what is left.
            let choices: Vec<Vec<usize>> = (0..k)
                .map(|i| {
                    let x = summand_x(&fb.subspace_basis, ell, values, i, kc.n);
                    points_with_x(&kc.curve, &x)
                        .iter()
                        .filter_map(|p| index_of.get(&point_key(p)).copied())
                        .collect()
                })
                .collect();
            if choices.iter().any(Vec::is_empty) {
                return None;
            }
            let mut pick = vec![0usize; k];
            loop {
                let chosen: Vec<usize> = (0..k).map(|i| choices[i][pick[i]]).collect();
                let rest = chosen.iter().fold(target.clone(), |acc, &i| {
                    kc.add(&acc, &negate(&fb.points[i]))
                });
                if let Some(mut head) =
                    decompose_rec(kc, fb, index_of, st, &rest, m - k, node_budget, opts, stats)
                {
                    head.extend(chosen);
                    return Some(head);
                }
                if stats.exhausted {
                    return None;
                }
                // Next sign pattern.
                let mut i = 0;
                loop {
                    if i == k {
                        return None;
                    }
                    pick[i] += 1;
                    if pick[i] < choices[i].len() {
                        break;
                    }
                    pick[i] = 0;
                    i += 1;
                }
            }
        },
    )
}

/// `−P` on a binary curve: `(x, x + y)`.
fn negate(p: &BinaryPoint) -> BinaryPoint {
    match p {
        BinaryPoint::Affine { x, y } => BinaryPoint::Affine {
            x: x.clone(),
            y: x.add(y),
        },
        BinaryPoint::Infinity => BinaryPoint::Infinity,
    }
}

/// Summand `i`'s abscissa under an assignment of the summand bits.
fn summand_x(basis: &[F2mElement], ell: usize, bits: u128, i: usize, n: u32) -> F2mElement {
    let mut acc = 0u64;
    for t in 0..ell {
        if bits >> (i * ell + t) & 1 == 1 {
            acc ^= basis[t].raw_bits().first().copied().unwrap_or(0);
        }
    }
    let positions: Vec<u32> = (0..n).filter(|&k| acc >> k & 1 == 1).collect();
    F2mElement::from_bit_positions(&positions, n)
}

/// Depth-first over the summand bits of `sys` (`count` summands), each
/// node linearised.  A node whose first `head` summands are all fixed is
/// handed to `leaf` with `is_head = true` (the last-summand shortcut); a
/// node with every summand fixed, with `false`.  `leaf` returning `Some`
/// ends the search.
#[allow(clippy::type_complexity)]
fn search(
    sys: &WideSystem,
    count: usize,
    head: Option<usize>,
    node_budget: usize,
    opts: &SearchOptions,
    stats: &mut WideStats,
    leaf: &mut dyn FnMut(u128, bool, &mut WideStats) -> Option<Vec<usize>>,
) -> Option<Vec<usize>> {
    let ell = sys.ell;
    let summand_vars = count * ell;
    struct Frame {
        eqs: Vec<WPoly>,
        subs: Vec<(usize, WPoly)>,
        fixed: u128,
        values: u128,
    }
    let mut stack = vec![Frame {
        eqs: sys.equations.clone(),
        subs: Vec::new(),
        fixed: 0,
        values: 0,
    }];
    let head_mask = head.map(|h| {
        let bits = h * ell;
        if bits >= 128 {
            u128::MAX
        } else {
            (1u128 << bits) - 1
        }
    });
    while let Some(mut f) = stack.pop() {
        if stats.nodes >= node_budget {
            stats.exhausted = true;
            return None;
        }
        stats.nodes += 1;
        if !linearise(&mut f.eqs, stats, &mut f.subs) {
            stats.refuted += 1;
            continue;
        }
        // Summand bits fixed by substitution to constants count as fixed.
        let mut fixed = f.fixed;
        let mut values = f.values;
        for (v, form) in &f.subs {
            if *v < summand_vars && form.terms.iter().all(|&m| m == 0) {
                fixed |= 1u128 << v;
                if form.is_one() {
                    values |= 1u128 << v;
                }
            }
        }
        if opts.order && !summands_ordered(fixed, values, count, ell) {
            stats.refuted += 1;
            continue;
        }
        if let Some(mask) = head_mask {
            if fixed & mask == mask {
                stats.leaves += 1;
                if let Some(found) = leaf(values, true, stats) {
                    return Some(found);
                }
                if stats.exhausted {
                    return None;
                }
                continue;
            }
        }
        // Each summand from its highest coordinate down, so the order
        // test above decides as early as it can.
        let next = if opts.high_first {
            (0..count)
                .flat_map(|i| (0..ell).rev().map(move |t| i * ell + t))
                .find(|&v| fixed >> v & 1 == 0)
        } else {
            (0..summand_vars).find(|&v| fixed >> v & 1 == 0)
        };
        let Some(v) = next else {
            stats.leaves += 1;
            if let Some(found) = leaf(values, false, stats) {
                return Some(found);
            }
            if stats.exhausted {
                return None;
            }
            continue;
        };
        for value in [true, false] {
            let form = if value { WPoly::one() } else { WPoly::zero() };
            let mut eqs: Vec<WPoly> = f.eqs.iter().map(|p| p.substitute(v, &form)).collect();
            eqs.retain(|p| !p.is_zero());
            let mut subs = f.subs.clone();
            for (_, g) in subs.iter_mut() {
                *g = g.substitute(v, &form);
            }
            subs.push((v, form));
            stack.push(Frame {
                eqs,
                subs,
                fixed: fixed | 1u128 << v,
                values: if value { values | 1u128 << v } else { values },
            });
        }
    }
    None
}

/// With the abscissae of all but the last summand known, find the last:
/// for every choice of factor-base points with those abscissae (each is
/// one of `±P`), the remainder `target − ΣP_i` must itself be a base
/// point.  Returns the indices, last summand included.
fn finish_last(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    xs: &[F2mElement],
    target: &BinaryPoint,
) -> Option<Vec<usize>> {
    let choices: Vec<Vec<usize>> = xs
        .iter()
        .map(|x| {
            points_with_x(&kc.curve, x)
                .iter()
                .filter_map(|p| index_of.get(&point_key(p)).copied())
                .collect()
        })
        .collect();
    if choices.iter().any(Vec::is_empty) {
        return None;
    }
    fn walk(
        kc: &KoblitzCurve,
        fb: &FrobeniusFactorBase,
        index_of: &HashMap<(BigUint, BigUint), usize>,
        choices: &[Vec<usize>],
        depth: usize,
        rest: &BinaryPoint,
        chosen: &mut Vec<usize>,
    ) -> bool {
        if depth == choices.len() {
            if let Some(&last) = index_of.get(&point_key(rest)) {
                chosen.push(last);
                return true;
            }
            return false;
        }
        for &i in &choices[depth] {
            let neg = negate(&fb.points[i]);
            chosen.push(i);
            if walk(
                kc,
                fb,
                index_of,
                choices,
                depth + 1,
                &kc.add(rest, &neg),
                chosen,
            ) {
                return true;
            }
            chosen.pop();
        }
        false
    }
    let mut chosen = Vec::with_capacity(xs.len() + 1);
    walk(kc, fb, index_of, &choices, 0, target, &mut chosen).then_some(chosen)
}

/// Whether the fixed bits already contradict `x_0 ≤ x_1 ≤ … ≤ x_{m−1}`,
/// comparing each adjacent pair from its highest coordinate down:
/// `false` only on a definite violation — the first coordinate where both
/// are fixed and differ has the earlier summand's bit set and the later
/// one's clear, with every higher coordinate fixed and equal in both.
fn summands_ordered(fixed: u128, values: u128, m: usize, ell: usize) -> bool {
    for i in 1..m {
        let (lo, hi) = ((i - 1) * ell, i * ell);
        for t in (0..ell).rev() {
            let (a, b) = (lo + t, hi + t);
            if fixed >> a & 1 == 0 || fixed >> b & 1 == 0 {
                break;
            }
            let (va, vb) = (values >> a & 1, values >> b & 1);
            if va != vb {
                if va > vb {
                    return false;
                }
                break;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::koblitz_groebner::build_decomposition_system;
    use crate::cryptanalysis::koblitz_index_calculus::build_frobenius_factor_base;
    use rand::{rngs::StdRng, Rng, SeedableRng};

    #[test]
    fn wide_system_matches_the_64_bit_builder() {
        // Where both fit, the equations are the same polynomials.
        for (a, n, m) in [
            (0u8, 9u32, 2usize),
            (0, 9, 3),
            (1, 17, 2),
            (0, 13, 3),
            (0, 9, 4),
        ] {
            let kc = KoblitzCurve::new(a, n).unwrap();
            let fb = build_frobenius_factor_base(&kc, 0).unwrap();
            let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
            let x = match kc.mul(kc.generator(), &BigUint::from(4242u32)) {
                BinaryPoint::Affine { x, .. } => x,
                _ => unreachable!(),
            };
            let Some(narrow) =
                build_decomposition_system(&fb.subspace_basis, &x, &kc.curve.b, m, &st)
            else {
                continue;
            };
            let wide = WideSystem::build(&fb.subspace_basis, &x, &kc.curve.b, m, &st).unwrap();
            assert_eq!(wide.n_vars, narrow.n_vars);
            let mut want: Vec<Vec<Mono>> = narrow
                .equations
                .iter()
                .filter(|p| !p.terms.is_empty())
                .map(|p| {
                    let mut t: Vec<Mono> = p.terms.iter().map(|m| u128::from(m.mask)).collect();
                    t.sort_unstable();
                    t
                })
                .collect();
            let mut got: Vec<Vec<Mono>> = wide.equations.iter().map(|p| p.terms.clone()).collect();
            want.sort();
            got.sort();
            assert_eq!(got, want, "K_{a}/2^{n} m={m}");
        }
    }

    #[test]
    fn substitution_and_arithmetic_agree_with_evaluation() {
        let mut rng = StdRng::seed_from_u64(5);
        for _ in 0..200 {
            let rand_poly = |rng: &mut StdRng| {
                WPoly::from_monos(
                    (0..rng.gen_range(0..8))
                        .map(|_| {
                            (0..3).fold(0u128, |m, _| {
                                if rng.gen_bool(0.7) {
                                    m | 1u128 << rng.gen_range(0..100)
                                } else {
                                    m
                                }
                            })
                        })
                        .collect(),
                )
            };
            let (p, q) = (rand_poly(&mut rng), rand_poly(&mut rng));
            let v = rng.gen_range(0..100);
            let form = WPoly::from_monos(vec![0, 1u128 << rng.gen_range(0..100)]);
            for _ in 0..20 {
                let bits: u128 = rng.gen();
                assert_eq!(p.add(&q).eval(bits), p.eval(bits) ^ q.eval(bits));
                assert_eq!(p.mul(&q).eval(bits), p.eval(bits) & q.eval(bits));
                // Substituting then evaluating is evaluating with v set to form(bits).
                let fv = form.eval(bits);
                let with = if fv {
                    bits | 1u128 << v
                } else {
                    bits & !(1u128 << v)
                };
                assert_eq!(p.substitute(v, &form).eval(bits), p.eval(with));
            }
        }
    }

    #[test]
    fn the_order_test_refuses_only_definite_violations() {
        let ell = 3;
        // Summand 0 = 0b101, summand 1 = 0b011: 5 > 3, a violation.
        let all = 0b111_111u128;
        assert!(!summands_ordered(all, 0b011_101, 2, ell));
        assert!(summands_ordered(all, 0b101_011, 2, ell));
        assert!(summands_ordered(all, 0b101_101, 2, ell)); // equal is allowed
                                                           // Top bits differ in the violating direction but a lower bit is
                                                           // unfixed: still a violation, the top bit decides.
        assert!(!summands_ordered(0b100_100, 0b000_100, 2, ell));
        // Top bit of the later summand unfixed: undecided, not refused.
        assert!(summands_ordered(0b000_111, 0b000_101, 2, ell));
    }

    #[test]
    fn the_suffix_recursion_finds_planted_decompositions() {
        // Force the recursion on a small curve by lowering the one-piece
        // cap below the whole chain: at n = 9, m = 4 the chain has
        // 4l + 18 unknowns, the suffix k(l + 9).
        for (a, n, m) in [(0u8, 9u32, 4usize), (0, 9, 5), (1, 17, 3)] {
            let kc = KoblitzCurve::new(a, n).unwrap();
            let fb = build_frobenius_factor_base(&kc, 0).unwrap();
            let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
            let index_of = fb.index_map();
            let ell = fb.subspace_basis.len();
            let opts = SearchOptions {
                max_vars: m * ell + (m - 2) * n as usize - 1,
                ..SearchOptions::from_env()
            };
            let mut rng = StdRng::seed_from_u64(91 + n as u64 + m as u64);
            for _ in 0..4 {
                let t = (0..m).fold(BinaryPoint::Infinity, |acc, _| {
                    kc.add(&acc, &fb.points[rng.gen_range(0..fb.points.len())])
                });
                if t == BinaryPoint::Infinity {
                    continue;
                }
                let mut stats = WideStats::default();
                let idxs =
                    decompose_rec(&kc, &fb, &index_of, &st, &t, m, 1 << 22, &opts, &mut stats)
                        .unwrap_or_else(|| panic!("K_{a}/2^{n} m={m}: missed, {stats:?}"));
                assert_eq!(idxs.len(), m);
                let sum = idxs
                    .iter()
                    .fold(BinaryPoint::Infinity, |s, &i| kc.add(&s, &fb.points[i]));
                assert_eq!(sum, t, "K_{a}/2^{n} m={m}");
            }
        }
    }

    #[test]
    fn finds_planted_decompositions_on_small_curves() {
        for (a, n, m) in [(0u8, 9u32, 2usize), (0, 9, 3), (1, 17, 2)] {
            let kc = KoblitzCurve::new(a, n).unwrap();
            let fb = build_frobenius_factor_base(&kc, 0).unwrap();
            let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
            let index_of = fb.index_map();
            let mut rng = StdRng::seed_from_u64(77 + n as u64);
            for _ in 0..6 {
                let t = (0..m).fold(BinaryPoint::Infinity, |acc, _| {
                    kc.add(&acc, &fb.points[rng.gen_range(0..fb.points.len())])
                });
                if t == BinaryPoint::Infinity {
                    continue;
                }
                let (found, stats) =
                    wide_groebner_decompose(&kc, &fb, &index_of, &st, &t, m, 1 << 20);
                let idxs = found.unwrap_or_else(|| panic!("K_{a}/2^{n} m={m}: missed, {stats:?}"));
                let sum = idxs
                    .iter()
                    .fold(BinaryPoint::Infinity, |s, &i| kc.add(&s, &fb.points[i]));
                assert_eq!(sum, t);
            }
        }
    }
}
