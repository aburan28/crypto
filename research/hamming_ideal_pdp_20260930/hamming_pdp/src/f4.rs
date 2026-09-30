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

use crate::boolpoly::{cmp_mono, Mono, Poly, W};
use std::collections::{BTreeMap, HashMap, HashSet};

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
    lm_index: HashMap<Mono, usize>,
    queue: BTreeMap<u32, Vec<Pair>>,
    d_max: u32,
    /// Row-words allowed in one matrix (memory cap; exceeding it is a
    /// budget failure like the XOR budget).
    matrix_cap_words: u64,
    /// Whether a new linear element ends the run for propagation.
    propagate: bool,
}

/// Largest polynomial a linear substitution may expand to before the
/// substitution is left to the matrix instead.
const MAX_SUBST_TERMS: usize = 1024;

/// Visit the submasks (divisors) of a square-free monomial; stop early
/// when the visitor returns `true`.
fn for_each_submask(m: &Mono, mut f: impl FnMut(&Mono) -> bool) {
    let vars = m.vars();
    let k = vars.len();
    for a in 0..(1u64 << k) {
        let mut d = Mono::ONE;
        for (i, v) in vars.iter().enumerate() {
            if (a >> i) & 1 == 1 {
                d = d.mul(&Mono::var(*v));
            }
        }
        if f(&d) {
            return;
        }
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
        Core { basis: vec![], alive: vec![], lm_index: HashMap::new(), queue: BTreeMap::new(), d_max, matrix_cap_words, propagate: true }
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
    fn find_reducer(&self, m: &Mono) -> Option<usize> {
        let mut best: Option<usize> = None;
        let mut best_len = usize::MAX;
        for_each_submask(m, |d| {
            if let Some(&gi) = self.lm_index.get(d) {
                let l = self.basis[gi].terms.len();
                if l < best_len {
                    best_len = l;
                    best = Some(gi);
                }
            }
            false
        });
        best
    }

    fn add_element(&mut self, p: Poly) -> bool {
        // returns false if p == 1
        if p.is_one() {
            return false;
        }
        let lm = *p.lm().unwrap();
        // A leading monomial some live element already divides: the
        // element is redundant as a basis element, but its content is not;
        // queue it as a row to be reduced.
        if self.find_reducer(&lm).is_some() {
            let d = p.degree();
            if d <= self.d_max {
                self.queue.entry(d).or_default().push(Pair::Row(p));
            }
            return true;
        }
        let k = self.basis.len();
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
        // New pairs (i, k) with criteria M and F, through a hash of the
        // candidate lcms: a candidate is dropped when a proper divisor of
        // its lcm is another candidate's lcm (M), or when an earlier
        // candidate has the same lcm (F).
        let mut cand: Vec<(usize, Mono)> = Vec::new();
        let mut first_with: HashMap<Mono, usize> = HashMap::new();
        for i in 0..k {
            if !self.alive[i] {
                continue;
            }
            let l = self.basis[i].lm().unwrap().mul(&lm);
            if l.degree() > self.d_max {
                continue;
            }
            first_with.entry(l).or_insert(cand.len());
            cand.push((i, l));
        }
        for (a, (i, l)) in cand.iter().enumerate() {
            if first_with[l] != a {
                continue; // F
            }
            let mut dropped = false;
            for_each_submask(l, |d| {
                if d != l && first_with.contains_key(d) {
                    dropped = true; // M
                    return true;
                }
                false
            });
            if dropped {
                continue;
            }
            self.queue.entry(l.degree()).or_default().push(Pair::S(*i, k, *l));
        }
        // Field pairs: (x_v + 1) · tail(p) for each variable of lm(p).
        if p.terms.len() > 1 {
            let tail = Poly { terms: p.terms[1..].to_vec() };
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
            let mut done: HashSet<Mono> = rows
                .iter()
                .zip(&multiple)
                .filter(|(_, m)| **m)
                .map(|(r, _)| *r.lm().unwrap())
                .collect();
            let mut seen: HashSet<Mono> = HashSet::new();
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
            let input_lms: HashSet<Mono> = rows
                .iter()
                .zip(&multiple)
                .filter(|(_, m)| **m)
                .map(|(r, _)| *r.lm().unwrap())
                .collect();

            let mut cols: Vec<Mono> = seen.into_iter().collect();
            cols.sort_by(|a, b| cmp_mono(b, a));
            let col_index: HashMap<Mono, usize> = cols.iter().enumerate().map(|(i, m)| (*m, i)).collect();
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
                    self.queue.iter().map(|(d, v)| (*d, v.len())).collect::<Vec<_>>()
                );
            }
            // Smallest first so that retirement happens early.
            new_polys.sort_by(|a, b| cmp_mono(a.lm().unwrap(), b.lm().unwrap()));
            let mut linear_found = false;
            for p in new_polys {
                if p.degree() <= 1 {
                    linear_found = true;
                }
                if !self.add_element(p) {
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

pub fn truncated_groebner(input: &[Poly], d_max: u32, budget_xor: u64, matrix_cap_words: u64) -> (Outcome, F4Stats) {
    let mut stats = F4Stats::default();
    let mut system: Vec<Poly> = input.iter().filter(|p| !p.is_zero()).cloned().collect();
    // Accumulated linear relations, kept in RREF over pivot variables.
    let mut linear_acc: Vec<Poly> = Vec::new();
    loop {
        let mut core = Core::new(d_max, matrix_cap_words);
        let t_add = std::time::Instant::now();
        for p in &system {
            if !core.add_element(p.clone()) {
                return (Outcome::Inconsistent, stats);
            }
        }
        stats.add_ns += t_add.elapsed().as_nanos() as u64;
        // Propagate any linear input before spending anything.
        let has_linear_input = system.iter().any(|p| p.degree() <= 1);
        let end = if has_linear_input { CoreEnd::Propagate } else { core.run(&mut stats, budget_xor) };
        let complete = match end {
            CoreEnd::Inconsistent => return (Outcome::Inconsistent, stats),
            CoreEnd::Budget => return (Outcome::Budget, stats),
            CoreEnd::Done => true,
            CoreEnd::Propagate => false,
        };
        // After a complete run the retired elements are reducible by the
        // live ones; after an early exit their pairs are still pending, so
        // they are carried along with the queued rows.
        let mut live: Vec<Poly> = core
            .basis
            .iter()
            .zip(&core.alive)
            .filter(|(_, a)| **a || !complete)
            .map(|(p, _)| p.clone())
            .collect();
        live.extend(core.drain_rows());
        stats.basis_len = live.len();
        let mut lin: Vec<Poly> = live.iter().filter(|p| p.degree() <= 1).cloned().collect();
        let nonlin: Vec<Poly> = live.iter().filter(|p| p.degree() > 1).cloned().collect();
        if lin.iter().any(|p| p.is_one()) {
            return (Outcome::Inconsistent, stats);
        }
        if lin.is_empty() {
            if std::env::var("F4_TRACE").is_ok() {
                for p in &live {
                    eprintln!("  basis: {p}");
                }
            }
            return (if nonlin.is_empty() { Outcome::Linear(rref_linear(&mut linear_acc).unwrap()) } else { Outcome::Wild }, stats);
        }
        // Propagate: fold the new linear relations into the accumulated
        // RREF, substitute into the nonlinear part, restart.
        let mut all = linear_acc.clone();
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
        for p in &nonlin {
            match substitute_rref(p, &map, MAX_SUBST_TERMS) {
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
            for p in &sys {
                if !core.add_element(p.clone()) {
                    return (Outcome::Inconsistent, stats);
                }
            }
            core.propagate = false;
            match core.run(&mut stats, budget_xor) {
                CoreEnd::Inconsistent => return (Outcome::Inconsistent, stats),
                CoreEnd::Budget => return (Outcome::Budget, stats),
                CoreEnd::Done | CoreEnd::Propagate => {}
            }
            let live: Vec<Poly> = core.basis.iter().zip(&core.alive).filter(|(_, a)| **a).map(|(p, _)| p.clone()).collect();
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
        .map(|l| (l.lm().unwrap().vars()[0], Poly { terms: l.terms[1..].to_vec() }))
        .collect()
}

pub fn substitute_rref(p: &Poly, map: &HashMap<usize, Poly>, max_terms: usize) -> Option<Poly> {
    let mut acc = Poly::zero();
    let mut plain: Vec<Mono> = Vec::new();
    for m in &p.terms {
        let vs: Vec<usize> = m.vars().into_iter().filter(|v| map.contains_key(v)).collect();
        if vs.is_empty() {
            plain.push(*m);
            continue;
        }
        let mut base = *m;
        for v in &vs {
            base = base.without(*v);
        }
        let mut prod = Poly { terms: vec![base] };
        for v in &vs {
            prod = prod.mul(&map[v]);
            if prod.terms.len() > max_terms {
                return None;
            }
            if prod.is_zero() {
                break;
            }
        }
        acc.add_assign(&prod);
        if acc.terms.len() > max_terms {
            return None;
        }
    }
    Some(acc.add(&Poly::from_terms(plain)))
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
            let rest = Poly { terms: r.terms[1..].to_vec() };
            if rest.eval(&pt) {
                pt[v / 64] |= 1 << (v % 64);
            }
        }
        sols.push(pt);
    }
    Some(sols)
}
