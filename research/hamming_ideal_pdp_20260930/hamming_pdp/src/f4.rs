//! A degree-truncated Boolean F4: the paper's `Groebner(H, d)` inside a
//! budget, and the tame/wild verdict of `GroebnerSafe`.
//!
//! Pairs are pruned with the Gebauer–Möller lcm criteria only (the product
//! criterion is false in the Boolean quotient); every field pair
//! `(x_v + 1)·tail(f)` is reduced. Elements whose leading monomial a newer
//! one divides are retired after their pair with it is queued. Linear
//! basis elements are propagated into the whole basis and the computation
//! restarted, which is the interreduction any F4 performs. The verdict is
//! sound whatever the truncation does: `1` is only reported when it is in
//! the ideal, and a linear verdict is checked by the caller against the
//! original equations.

use crate::boolpoly::{active_limbs, cmp_mono, limbs_covering, LimbsGuard, Mono, Poly, W};
use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};

/// Multiply-add mixer. SipHash of a wide monomial dominates the reducer
/// index; the inputs are our own monomials, not attacker-chosen keys.
#[derive(Clone, Default)]
struct FxHasher(u64);

impl Hasher for FxHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.write_u64(b as u64);
        }
    }
    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.0 = self.0.rotate_left(5) ^ i.wrapping_mul(0x517cc1b727220a95);
    }
    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.write_u64(i as u64);
    }
    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.write_u64(i as u64);
    }
}

type FxBuild = BuildHasherDefault<FxHasher>;
type MonoMap<V> = HashMap<Mono, V, FxBuild>;
type MonoSet = HashSet<Mono, FxBuild>;

fn mono_map<V>() -> MonoMap<V> {
    HashMap::with_hasher(FxBuild::default())
}
fn mono_set() -> MonoSet {
    HashSet::with_hasher(FxBuild::default())
}

#[derive(Clone, Debug, Default)]
pub struct F4Stats {
    pub steps: u32,
    pub rows_max: usize,
    pub cols_max: usize,
    pub xor_words: u64,
    pub max_degree_built: u32,
    pub basis_len: usize,
    pub pairs_reduced: u64,
    pub restarts: u32,
    pub matrix_cap_hits: u64,
    pub prep_ns: u64,
    pub elim_ns: u64,
    pub add_ns: u64,
    pub rref_ns: u64,
    pub subst_ns: u64,
    pub post_ns: u64,
}

pub enum Outcome {
    /// `1` is in the ideal: no `F_2`-solution.
    Inconsistent,
    /// Every basis element reduces to a linear polynomial: the RREF rows.
    Linear(Vec<Poly>),
    /// Truncated basis is not linear.
    Wild,
    /// The XOR budget was exhausted (the paper's timeout).
    Budget,
}

enum Pair {
    S(usize, usize, Mono),
    Row(Poly),
}

fn lowest_col(bits: &[u64]) -> Option<usize> {
    for (k, w) in bits.iter().enumerate() {
        if *w != 0 {
            return Some(k * 64 + w.trailing_zeros() as usize);
        }
    }
    None
}

struct Core {
    basis: Vec<Poly>,
    alive: Vec<bool>,
    /// Leading monomial of each live element -> its index.
    lm_index: MonoMap<usize>,
    queue: BTreeMap<u32, Vec<Pair>>,
    d_max: u32,
    /// Row-words allowed in one matrix (memory cap; exceeding it is a
    /// budget failure like the XOR budget).
    matrix_cap_words: u64,
    /// Whether a new linear element ends the run for propagation.
    propagate: bool,
    /// Scratch for Gebauer–Möller pair selection, reused so a basis of
    /// thousands of elements does not allocate a map per element.
    spair_cand: Vec<(usize, Mono)>,
    spair_by_deg: Vec<Vec<usize>>,
    spair_keep: Vec<bool>,
    spair_minimal: Vec<Mono>,
    spair_first: MonoMap<usize>,
}

/// Largest polynomial a linear substitution may expand to before the
/// substitution is left to the matrix instead.
const MAX_SUBST_TERMS: usize = 1024;

/// Set-bit positions of `m`, low variable first.
/// `None` when the degree exceeds 64 (the stack buffer).
fn support64(m: &Mono, idx: &mut [u16; 64]) -> Option<usize> {
    let mut k = 0usize;
    let limbs = active_limbs();
    for limb in 0..limbs {
        let mut w = m.0[limb];
        while w != 0 {
            if k == 64 {
                return None;
            }
            let b = w.trailing_zeros();
            idx[k] = (limb as u16) * 64 + b as u16;
            k += 1;
            w &= w - 1;
        }
    }
    Some(k)
}

fn support_all(m: &Mono) -> Vec<u16> {
    let mut idx = Vec::new();
    let limbs = active_limbs();
    for limb in 0..limbs {
        let mut w = m.0[limb];
        while w != 0 {
            let b = w.trailing_zeros();
            idx.push((limb as u16) * 64 + b as u16);
            w &= w - 1;
        }
    }
    idx
}

/// `true` when the subset `a` precedes `b` in increasing mask order on `idx`.
fn earlier_subset(a: &Mono, b: &Mono, idx: &[u16]) -> bool {
    for &vi in idx {
        let v = vi as usize;
        let ba = (a.0[v / 64] >> (v % 64)) & 1;
        let bb = (b.0[v / 64] >> (v % 64)) & 1;
        if ba != bb {
            return ba < bb;
        }
    }
    false
}

/// Visit the submasks (divisors) of a square-free monomial in increasing
/// subset-index order; stop early when the visitor returns `true`.
/// Criterion M only calls this on lcms of degree `≤ d_max` (6).
fn for_each_submask(m: &Mono, mut f: impl FnMut(&Mono) -> bool) {
    let mut idx = [0u16; 64];
    let k = support64(m, &mut idx).expect("submask enumeration past degree 64");
    let mut cur = [0u64; W];
    let total = 1u64 << k;
    let mut a = 0u64;
    loop {
        if f(&Mono(cur)) {
            return;
        }
        a += 1;
        if a == total {
            return;
        }
        let tz = a.trailing_zeros() as usize;
        for i in 0..tz {
            let v = idx[i] as usize;
            cur[v / 64] &= !(1u64 << (v % 64));
        }
        let v = idx[tz] as usize;
        cur[v / 64] |= 1u64 << (v % 64);
    }
}

enum CoreEnd {
    Done,
    Inconsistent,
    Budget,
    /// A linear element appeared: propagate it before spending more.
    Propagate,
}

impl Core {
    fn new(d_max: u32, matrix_cap_words: u64) -> Core {
        Core {
            basis: vec![],
            alive: vec![],
            lm_index: mono_map(),
            queue: BTreeMap::new(),
            d_max,
            matrix_cap_words,
            propagate: true,
            spair_cand: Vec::new(),
            spair_by_deg: vec![Vec::new(); d_max as usize + 1],
            spair_keep: Vec::new(),
            spair_minimal: Vec::new(),
            spair_first: mono_map(),
        }
    }

    /// Queued singleton rows (field pairs, equal-leading-monomial
    /// differences): relations not yet in the basis, carried across a
    /// propagation restart so that nothing queued is lost.
    fn drain_rows(&mut self) -> Vec<Poly> {
        let mut out = Vec::new();
        for (_, v) in self.queue.iter_mut() {
            for pr in v.drain(..) {
                if let Pair::Row(p) = pr {
                    out.push(p);
                }
            }
        }
        out
    }

    /// The live reducer of `m` with the fewest terms, if any.
    ///
    /// Ties keep the divisor visited first in subset-index order (the order
    /// `for_each_submask` uses). A one-term reducer is optimal, so the scan
    /// stops at the first one. When the leading-monomial index is smaller
    /// than the divisor lattice, the same `(length, subset index)` order is
    /// applied by walking the index instead of the lattice.
    fn find_reducer(&self, m: &Mono) -> Option<usize> {
        let mut idx = [0u16; 64];
        let k_opt = support64(m, &mut idx);
        let enumerate = match k_opt {
            Some(k) if k < 20 => (1usize << k) <= self.lm_index.len().saturating_mul(2),
            _ => false,
        };
        if enumerate {
            let k = k_opt.unwrap();
            let mut cur = [0u64; W];
            let total = 1u64 << k;
            let mut a = 0u64;
            let mut best: Option<usize> = None;
            let mut best_len = usize::MAX;
            loop {
                if let Some(&gi) = self.lm_index.get(&Mono(cur)) {
                    let l = self.basis[gi].terms.len();
                    if l < best_len {
                        best_len = l;
                        best = Some(gi);
                        if best_len == 1 {
                            return best;
                        }
                    }
                }
                a += 1;
                if a == total {
                    break;
                }
                let tz = a.trailing_zeros() as usize;
                for i in 0..tz {
                    let v = idx[i] as usize;
                    cur[v / 64] &= !(1u64 << (v % 64));
                }
                let v = idx[tz] as usize;
                cur[v / 64] |= 1u64 << (v % 64);
            }
            return best;
        }
        let owned;
        let idx_ref: &[u16] = match k_opt {
            Some(k) => &idx[..k],
            None => {
                owned = support_all(m);
                &owned
            }
        };
        self.find_reducer_by_index(m, idx_ref)
    }

    fn find_reducer_by_index(&self, m: &Mono, idx: &[u16]) -> Option<usize> {
        let mut best: Option<usize> = None;
        let mut best_lm = Mono::ONE;
        let mut best_len = usize::MAX;
        for (lm, &gi) in &self.lm_index {
            if !lm.divides(m) {
                continue;
            }
            let l = self.basis[gi].terms.len();
            if l < best_len || (l == best_len && earlier_subset(lm, &best_lm, idx)) {
                best_len = l;
                best_lm = *lm;
                best = Some(gi);
            }
        }
        best
    }

    fn add_element(&mut self, p: Poly, queue_spairs: bool) -> bool {
        // returns false if p == 1
        if p.is_one() {
            return false;
        }
        let lm = *p.lm().unwrap();
        let reducible = self.find_reducer(&lm).is_some();
        // A leading monomial some live element already divides: the
        // element is redundant as a basis element, but its content is not;
        // queue it as a row to be reduced.
        if reducible {
            let d = p.degree();
            if d <= self.d_max {
                self.queue.entry(d).or_default().push(Pair::Row(p));
            }
            return true;
        }
        let k = self.basis.len();
        // A propagation restart discards every S-pair (`drain_rows` keeps
        // only rows) and keeps retired elements anyway. Gebauer–Möller and
        // the new S-pair queue are then dead work. Field-pair rows, the
        // redundant-element rows above, and retirement still decide which
        // polynomials land in the restarted system and in what order, so
        // those stay. `queue_spairs` is false only on that path.
        if queue_spairs {
            // Gebauer–Möller criterion B on old pairs: drop (i, j) when lm_k
            // divides lcm(i, j) strictly finer than both lcm(i, k), lcm(j, k).
            // Only queue degrees at or above deg(lm_k) can be affected.
            let lm_deg = lm.degree();
            for (d, v) in self.queue.iter_mut() {
                if *d < lm_deg {
                    continue;
                }
                v.retain(|pr| match pr {
                    Pair::S(i, j, l) => {
                        if !lm.divides(l) {
                            return true;
                        }
                        let lik = self.basis[*i].lm().unwrap().mul(&lm);
                        let ljk = self.basis[*j].lm().unwrap().mul(&lm);
                        !(lik != *l && ljk != *l)
                    }
                    Pair::Row(_) => true,
                });
            }
            // New pairs (i, k) with criteria M and F. F keeps the earliest
            // partner with a given lcm. M drops an lcm that a different
            // candidate lcm properly divides. A proper divisor has lower
            // degree, so the division-minimal lcms of lower degree are a
            // complete test; same-degree lcms are an antichain and only F
            // applies. Once that minimal set is larger than a submask scan,
            // fall back to probing divisors in `first_with`.
            self.spair_cand.clear();
            self.spair_first.clear();
            self.spair_minimal.clear();
            for b in &mut self.spair_by_deg {
                b.clear();
            }
            for i in 0..k {
                if !self.alive[i] {
                    continue;
                }
                let l = self.basis[i].lm().unwrap().mul(&lm);
                let d = l.degree();
                if d > self.d_max {
                    continue;
                }
                let a = self.spair_cand.len();
                if let Entry::Vacant(e) = self.spair_first.entry(l) {
                    e.insert(a);
                    self.spair_by_deg[d as usize].push(a);
                }
                self.spair_cand.push((i, l));
            }
            self.spair_keep.clear();
            self.spair_keep.resize(self.spair_cand.len(), false);
            for di in 0..self.spair_by_deg.len() {
                let mut born: Vec<Mono> = Vec::new();
                for aj in 0..self.spair_by_deg[di].len() {
                    let a = self.spair_by_deg[di][aj];
                    let l = self.spair_cand[a].1;
                    let dropped = if self.spair_minimal.len() <= 40 {
                        self.spair_minimal.iter().any(|s| s.divides(&l))
                    } else {
                        let mut yes = false;
                        let first = &self.spair_first;
                        for_each_submask(&l, |d| {
                            if d != &l && first.contains_key(d) {
                                yes = true;
                                return true;
                            }
                            false
                        });
                        yes
                    };
                    if !dropped {
                        self.spair_keep[a] = true;
                        born.push(l);
                    }
                }
                self.spair_minimal.extend(born);
            }
            for (a, (i, l)) in self.spair_cand.iter().enumerate() {
                if self.spair_keep[a] {
                    self.queue
                        .entry(l.degree())
                        .or_default()
                        .push(Pair::S(*i, k, *l));
                }
            }
        }
        // Field pairs: (x_v + 1) · tail(p) for each variable of lm(p).
        if p.terms.len() > 1 {
            let tail = Poly {
                terms: p.terms[1..].to_vec(),
            };
            for v in lm.vars() {
                let row = tail.mul_mono(&Mono::var(v)).add(&tail);
                if !row.is_zero() {
                    let d = row.degree();
                    if d <= self.d_max {
                        self.queue.entry(d).or_default().push(Pair::Row(row));
                    }
                }
            }
        }
        // Retire elements whose leading monomial the new one divides; their
        // pair with the new element is queued above (criterion M keeps the
        // pair with lcm = lm_i unless an equal-lcm pair precedes it).
        for i in 0..k {
            if self.alive[i] && lm.divides(self.basis[i].lm().unwrap()) {
                self.alive[i] = false;
                self.lm_index.remove(self.basis[i].lm().unwrap());
            }
        }
        self.lm_index.insert(lm, k);
        self.basis.push(p);
        self.alive.push(true);
        true
    }

    fn run(&mut self, stats: &mut F4Stats, budget_xor: u64) -> CoreEnd {
        loop {
            let deg = match self.queue.keys().next() {
                Some(&d) if d <= self.d_max => d,
                _ => return CoreEnd::Done,
            };
            let pairs = self.queue.remove(&deg).unwrap();
            let t_step = std::time::Instant::now();
            stats.steps += 1;
            stats.max_degree_built = stats.max_degree_built.max(deg);

            let mut rows: Vec<Poly> = Vec::new();
            let mut multiple: Vec<bool> = Vec::new();
            for pr in pairs {
                match pr {
                    Pair::S(i, j, l) => {
                        if !self.alive[i] && !self.alive[j] {
                            continue;
                        }
                        let li = *self.basis[i].lm().unwrap();
                        let lj = *self.basis[j].lm().unwrap();
                        rows.push(self.basis[i].mul_mono(&l.div(&li)));
                        rows.push(self.basis[j].mul_mono(&l.div(&lj)));
                        multiple.push(true);
                        multiple.push(true);
                        stats.pairs_reduced += 1;
                    }
                    Pair::Row(p) => {
                        rows.push(p);
                        multiple.push(false);
                        stats.pairs_reduced += 1;
                    }
                }
            }
            {
                let keep: Vec<bool> = rows.iter().map(|r| !r.is_zero()).collect();
                let mut k = 0;
                rows.retain(|_| {
                    k += 1;
                    keep[k - 1]
                });
                k = 0;
                multiple.retain(|_| {
                    k += 1;
                    keep[k - 1]
                });
            }
            if rows.is_empty() {
                continue;
            }

            // Symbolic preprocessing: every monomial that is not the leading
            // monomial of a multiple row gets a reducer if one exists
            // (the live element with the fewest terms).
            let mut done = mono_set();
            done.extend(
                rows.iter()
                    .zip(&multiple)
                    .filter(|(_, m)| **m)
                    .map(|(r, _)| *r.lm().unwrap()),
            );
            let mut seen = mono_set();
            let mut todo: Vec<Mono> = Vec::new();
            for r in &rows {
                for m in &r.terms {
                    if seen.insert(*m) {
                        todo.push(*m);
                    }
                }
            }
            while let Some(m) = todo.pop() {
                if done.contains(&m) {
                    continue;
                }
                // The matrix cap, checked as the matrix grows: rows and
                // columns only increase during preprocessing, so an
                // early exceedance is a final one.
                if (rows.len() as u64) * (seen.len().div_ceil(64) as u64) > self.matrix_cap_words {
                    stats.matrix_cap_hits += 1;
                    stats.prep_ns += t_step.elapsed().as_nanos() as u64;
                    return CoreEnd::Budget;
                }
                if let Some(gi) = self.find_reducer(&m) {
                    let row = self.basis[gi].mul_mono(&m.div(self.basis[gi].lm().unwrap()));
                    done.insert(m);
                    for t in &row.terms {
                        if seen.insert(*t) {
                            todo.push(*t);
                        }
                    }
                    rows.push(row);
                    multiple.push(true);
                }
            }
            let mut input_lms = mono_set();
            input_lms.extend(
                rows.iter()
                    .zip(&multiple)
                    .filter(|(_, m)| **m)
                    .map(|(r, _)| *r.lm().unwrap()),
            );

            let mut cols: Vec<Mono> = seen.into_iter().collect();
            cols.sort_by(|a, b| cmp_mono(b, a));
            let mut col_index: MonoMap<usize> = mono_map();
            for (i, m) in cols.iter().enumerate() {
                col_index.insert(*m, i);
            }
            let ncols = cols.len();
            let words = ncols.div_ceil(64);
            stats.rows_max = stats.rows_max.max(rows.len());
            stats.cols_max = stats.cols_max.max(ncols);
            if (rows.len() as u64) * (words as u64) > self.matrix_cap_words {
                stats.matrix_cap_hits += 1;
                return CoreEnd::Budget;
            }
            stats.prep_ns += t_step.elapsed().as_nanos() as u64;
            let t_elim = std::time::Instant::now();

            let mut mat: Vec<Vec<u64>> = Vec::with_capacity(rows.len());
            for r in &rows {
                let mut bits = vec![0u64; words];
                for m in &r.terms {
                    let c = col_index[m];
                    bits[c / 64] |= 1 << (c % 64);
                }
                mat.push(bits);
            }
            let mut order: Vec<usize> = (0..mat.len()).collect();
            order.sort_by_key(|&i| lowest_col(&mat[i]).unwrap_or(usize::MAX));
            let mut pivots: Vec<Option<usize>> = vec![None; ncols];
            let mut pivot_rows: Vec<usize> = Vec::new();
            for i in order {
                loop {
                    let c = match lowest_col(&mat[i]) {
                        Some(c) => c,
                        None => break,
                    };
                    match pivots[c] {
                        Some(p) => {
                            let (a, b) = if p < i {
                                let (l, r) = mat.split_at_mut(i);
                                (&mut r[0], &l[p])
                            } else {
                                let (l, r) = mat.split_at_mut(p);
                                (&mut l[i], &r[0])
                            };
                            let start = c / 64;
                            for k in start..words {
                                a[k] ^= b[k];
                            }
                            stats.xor_words += (words - start) as u64;
                            if stats.xor_words > budget_xor {
                                return CoreEnd::Budget;
                            }
                        }
                        None => {
                            pivots[c] = Some(i);
                            pivot_rows.push(i);
                            break;
                        }
                    }
                }
            }

            stats.elim_ns += t_elim.elapsed().as_nanos() as u64;
            let t_post = std::time::Instant::now();
            let mut new_polys: Vec<Poly> = Vec::new();
            for &i in &pivot_rows {
                let c = lowest_col(&mat[i]).unwrap();
                if input_lms.contains(&cols[c]) {
                    continue;
                }
                let mut terms = Vec::new();
                for (k, w) in mat[i].iter().enumerate() {
                    let mut w = *w;
                    while w != 0 {
                        let b = w.trailing_zeros() as usize;
                        terms.push(cols[k * 64 + b]);
                        w &= w - 1;
                    }
                }
                new_polys.push(Poly::from_terms(terms));
            }
            if std::env::var("F4_TRACE").is_ok() {
                eprintln!(
                    "  f4 step {} deg {} rows {} cols {} new {} basis {} queue {:?}",
                    stats.steps,
                    deg,
                    rows.len(),
                    ncols,
                    new_polys.len(),
                    self.alive.iter().filter(|a| **a).count(),
                    self.queue
                        .iter()
                        .map(|(d, v)| (*d, v.len()))
                        .collect::<Vec<_>>()
                );
            }
            // Smallest first so that retirement happens early.
            new_polys.sort_by(|a, b| cmp_mono(a.lm().unwrap(), b.lm().unwrap()));
            let mut linear_found = false;
            for p in new_polys {
                if p.degree() <= 1 {
                    linear_found = true;
                }
                if !self.add_element(p, true) {
                    return CoreEnd::Inconsistent;
                }
            }
            stats.post_ns += t_post.elapsed().as_nanos() as u64;
            if linear_found && self.propagate {
                return CoreEnd::Propagate;
            }
        }
    }
}

pub fn truncated_groebner(
    input: &[Poly],
    d_max: u32,
    budget_xor: u64,
    matrix_cap_words: u64,
) -> (Outcome, F4Stats) {
    // Every monomial this call builds is an OR of monomials already in
    // `input`, so the limb width is fixed for the whole substitution loop.
    let _limbs = LimbsGuard::set(limbs_covering(input));
    let mut stats = F4Stats::default();
    let mut system: Vec<Poly> = input.iter().filter(|p| !p.is_zero()).cloned().collect();
    // Accumulated linear relations, kept in RREF over pivot variables.
    let mut linear_acc: Vec<Poly> = Vec::new();
    loop {
        let mut core = Core::new(d_max, matrix_cap_words);
        let t_add = std::time::Instant::now();
        // Known before insertion: a linear input returns Propagate without
        // `core.run`, and that path drops every S-pair.
        let has_linear_input = system.iter().any(|p| p.degree() <= 1);
        for p in system.drain(..) {
            if !core.add_element(p, !has_linear_input) {
                return (Outcome::Inconsistent, stats);
            }
        }
        stats.add_ns += t_add.elapsed().as_nanos() as u64;
        // Propagate any linear input before spending anything.
        let end = if has_linear_input {
            CoreEnd::Propagate
        } else {
            core.run(&mut stats, budget_xor)
        };
        let complete = match end {
            CoreEnd::Inconsistent => return (Outcome::Inconsistent, stats),
            CoreEnd::Budget => return (Outcome::Budget, stats),
            CoreEnd::Done => true,
            CoreEnd::Propagate => false,
        };
        // After a complete run the retired elements are reducible by the
        // live ones; after an early exit their pairs are still pending, so
        // they are carried along with the queued rows.
        let alive = std::mem::take(&mut core.alive);
        let basis = std::mem::take(&mut core.basis);
        let mut live: Vec<Poly> = basis
            .into_iter()
            .zip(alive)
            .filter(|(_, a)| *a || !complete)
            .map(|(p, _)| p)
            .collect();
        live.extend(core.drain_rows());
        stats.basis_len = live.len();
        let mut lin: Vec<Poly> = Vec::new();
        let mut nonlin: Vec<Poly> = Vec::new();
        for p in live {
            if p.degree() <= 1 {
                lin.push(p);
            } else {
                nonlin.push(p);
            }
        }
        if lin.iter().any(|p| p.is_one()) {
            return (Outcome::Inconsistent, stats);
        }
        if lin.is_empty() {
            if std::env::var("F4_TRACE").is_ok() {
                for p in &nonlin {
                    eprintln!("  basis: {p}");
                }
            }
            return (
                if nonlin.is_empty() {
                    Outcome::Linear(rref_linear(&mut linear_acc).unwrap())
                } else {
                    Outcome::Wild
                },
                stats,
            );
        }
        // Propagate: fold the new linear relations into the accumulated
        // RREF, substitute into the nonlinear part, restart.
        let mut all = std::mem::take(&mut linear_acc);
        all.append(&mut lin);
        let t_rref = std::time::Instant::now();
        let rref = match rref_linear(&mut all) {
            Some(r) => r,
            None => return (Outcome::Inconsistent, stats),
        };
        stats.rref_ns += t_rref.elapsed().as_nanos() as u64;
        let t_subst = std::time::Instant::now();
        let mut next: Vec<Poly> = Vec::new();
        let mut too_big = false;
        let map = subst_map(&rref);
        let mut prepared = PreparedSubst::from_map(&map);
        for p in &nonlin {
            match prepared.apply(p, MAX_SUBST_TERMS) {
                Some(q) => {
                    if q.is_one() {
                        return (Outcome::Inconsistent, stats);
                    }
                    if !q.is_zero() {
                        next.push(q);
                    }
                }
                None => {
                    too_big = true;
                    break;
                }
            }
        }
        stats.subst_ns += t_subst.elapsed().as_nanos() as u64;
        if too_big {
            // Keep the linear relations as basis elements and let the
            // matrix do the reduction: no further propagation restarts.
            let mut sys = nonlin;
            sys.extend(rref.iter().cloned());
            let mut core = Core::new(d_max, matrix_cap_words);
            for p in sys.drain(..) {
                if !core.add_element(p, true) {
                    return (Outcome::Inconsistent, stats);
                }
            }
            core.propagate = false;
            match core.run(&mut stats, budget_xor) {
                CoreEnd::Inconsistent => return (Outcome::Inconsistent, stats),
                CoreEnd::Budget => return (Outcome::Budget, stats),
                CoreEnd::Done | CoreEnd::Propagate => {}
            }
            let live: Vec<Poly> = core
                .basis
                .iter()
                .zip(&core.alive)
                .filter(|(_, a)| **a)
                .map(|(p, _)| p.clone())
                .collect();
            if live.iter().any(|p| p.is_one()) {
                return (Outcome::Inconsistent, stats);
            }
            let mut lin2: Vec<Poly> = live.iter().filter(|p| p.degree() <= 1).cloned().collect();
            let rref2 = match rref_linear(&mut lin2) {
                Some(r) => r,
                None => return (Outcome::Inconsistent, stats),
            };
            let map2 = subst_map(&rref2);
            for p in live.iter().filter(|p| p.degree() > 1) {
                match substitute_rref(p, &map2, MAX_SUBST_TERMS) {
                    Some(q) if q.is_zero() => {}
                    Some(q) if q.is_one() => return (Outcome::Inconsistent, stats),
                    _ => return (Outcome::Wild, stats),
                }
            }
            return (Outcome::Linear(rref2), stats);
        }
        linear_acc = rref;
        if next.is_empty() {
            if std::env::var("F4_TRACE").is_ok() {
                for p in &linear_acc {
                    eprintln!("  basis: {p}");
                }
            }
            return (Outcome::Linear(linear_acc), stats);
        }
        // Substituted polynomials may now be linear themselves; loop.
        system = next;
        stats.restarts += 1;
    }
}

/// Substitute every pivot variable of an RREF system at once; `None` if an
/// intermediate product exceeds `max_terms` (the substitution is then left
/// to the matrix, where a linear reducer costs one row per monomial).
pub fn subst_map(rref: &[Poly]) -> HashMap<usize, Poly> {
    rref.iter()
        .map(|l| {
            (
                l.lm().unwrap().vars()[0],
                Poly {
                    terms: l.terms[1..].to_vec(),
                },
            )
        })
        .collect()
}

/// Pivot tails indexed by variable, so a term test is a bit-mask instead of
/// a `vars()` allocation and a hash lookup per variable. Multiplication
/// order is ascending variable index, matching `m.vars()`.
struct PreparedSubst {
    pivot: [u64; W],
    tails: Vec<Poly>,
    /// Merge buffer for the substitution accumulator. `Poly::add` would
    /// allocate a fresh vector per term; this one is swapped back in.
    scratch: Vec<Mono>,
}

impl PreparedSubst {
    fn from_map(map: &HashMap<usize, Poly>) -> Self {
        let mut max_v = 0usize;
        for &v in map.keys() {
            max_v = max_v.max(v);
        }
        let mut tails = vec![Poly::zero(); max_v + 1];
        let mut pivot = [0u64; W];
        for (&v, tail) in map {
            pivot[v / 64] |= 1u64 << (v % 64);
            tails[v] = tail.clone();
        }
        Self {
            pivot,
            tails,
            scratch: Vec::new(),
        }
    }

    fn xor_acc(&mut self, acc: &mut Poly, extra: &Poly) {
        self.scratch.clear();
        self.scratch.reserve(acc.terms.len() + extra.terms.len());
        let (mut i, mut j) = (0, 0);
        while i < acc.terms.len() && j < extra.terms.len() {
            match cmp_mono(&acc.terms[i], &extra.terms[j]) {
                Ordering::Greater => {
                    self.scratch.push(acc.terms[i]);
                    i += 1;
                }
                Ordering::Less => {
                    self.scratch.push(extra.terms[j]);
                    j += 1;
                }
                Ordering::Equal => {
                    i += 1;
                    j += 1;
                }
            }
        }
        self.scratch.extend_from_slice(&acc.terms[i..]);
        self.scratch.extend_from_slice(&extra.terms[j..]);
        std::mem::swap(&mut acc.terms, &mut self.scratch);
    }

    fn apply(&mut self, p: &Poly, max_terms: usize) -> Option<Poly> {
        let mut acc = Poly::zero();
        let mut plain: Vec<Mono> = Vec::new();
        let limbs = active_limbs();
        for m in &p.terms {
            let mut hit = [0u64; W];
            let mut any = false;
            for k in 0..limbs {
                hit[k] = m.0[k] & self.pivot[k];
                any |= hit[k] != 0;
            }
            if !any {
                plain.push(*m);
                continue;
            }
            let mut base = m.0;
            for k in 0..limbs {
                base[k] &= !self.pivot[k];
            }
            let mut prod = Poly {
                terms: vec![Mono(base)],
            };
            let mut vanished = false;
            for k in 0..limbs {
                let mut w = hit[k];
                while w != 0 {
                    let b = w.trailing_zeros() as usize;
                    let v = k * 64 + b;
                    prod = prod.mul(&self.tails[v]);
                    if prod.terms.len() > max_terms {
                        return None;
                    }
                    if prod.is_zero() {
                        vanished = true;
                        break;
                    }
                    w &= w - 1;
                }
                if vanished {
                    break;
                }
            }
            if !vanished {
                self.xor_acc(&mut acc, &prod);
                if acc.terms.len() > max_terms {
                    return None;
                }
            }
        }
        let plain = Poly::from_terms(plain);
        self.xor_acc(&mut acc, &plain);
        Some(acc)
    }
}

pub fn substitute_rref(p: &Poly, map: &HashMap<usize, Poly>, max_terms: usize) -> Option<Poly> {
    PreparedSubst::from_map(map).apply(p, max_terms)
}

/// Substitute the polynomial `rest` for variable `v` (i.e. `x_v := rest`).
pub fn substitute_poly(p: &Poly, v: usize, rest: &Poly) -> Poly {
    let mut with = Vec::new();
    let mut without = Vec::new();
    for m in &p.terms {
        if m.has(v) {
            with.push(m.without(v));
        } else {
            without.push(*m);
        }
    }
    if with.is_empty() {
        return p.clone();
    }
    let a = Poly::from_terms(with).mul(rest);
    a.add(&Poly { terms: without })
}

/// Reduced row echelon form of linear polynomials (each `x_{v} + …`), with
/// the pivot variable as leading term; `None` if `1 = 0` appears.
pub fn rref_linear(lin: &mut Vec<Poly>) -> Option<Vec<Poly>> {
    let mut out: Vec<Poly> = Vec::new();
    for p in lin.drain(..) {
        let mut q = p;
        loop {
            let mut changed = false;
            for r in &out {
                let pv = r.lm().unwrap();
                if q.terms.iter().any(|m| m == pv) {
                    q = q.add(r);
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        if q.is_zero() {
            continue;
        }
        if q.is_one() {
            return None;
        }
        let pv = *q.lm().unwrap();
        for r in out.iter_mut() {
            if r.terms.iter().any(|m| *m == pv) {
                *r = r.add(&q);
            }
        }
        out.push(q);
    }
    Some(out)
}

/// Enumerate the affine solution set of an RREF linear system over the
/// variables `0..n_vars` if it has at most `2^max_free` points.
pub fn linear_solutions(rref: &[Poly], n_vars: usize, max_free: u32) -> Option<Vec<[u64; W]>> {
    let mut is_pivot = vec![false; n_vars];
    for r in rref {
        is_pivot[r.lm().unwrap().vars()[0]] = true;
    }
    let free: Vec<usize> = (0..n_vars).filter(|v| !is_pivot[*v]).collect();
    if free.len() as u32 > max_free {
        return None;
    }
    let mut sols = Vec::new();
    for a in 0u64..(1u64 << free.len()) {
        let mut pt = [0u64; W];
        for (i, &v) in free.iter().enumerate() {
            if (a >> i) & 1 == 1 {
                pt[v / 64] |= 1 << (v % 64);
            }
        }
        for r in rref {
            let v = r.lm().unwrap().vars()[0];
            let rest = Poly {
                terms: r.terms[1..].to_vec(),
            };
            if rest.eval(&pt) {
                pt[v / 64] |= 1 << (v % 64);
            }
        }
        sols.push(pt);
    }
    Some(sols)
}
