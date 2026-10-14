//! Prime-track root-finders for the PDP on tower factor bases
//! (`y² = x³ + ax + b` over `F_p`, `p < 2^62`, `V = {g·ζ^e}` for a Kummer
//! tower or an explicit list), by the group law on the Montgomery one-word
//! field, never by a polynomial system.
//!
//! * `subtract2` (`m = 2`): `R − P_k ∈ F`?  `|F|` batched additions and
//!   lookups in the base's abscissa set.
//! * `pairs3` (`m = 3`, memoryless): for each `P₁` the two row points
//!   `U± = R ± P₁`, then for each `P₂ ≥ P₁` the candidates `x(U± ± P₂)`,
//!   tested against the abscissa set of `V` in `O(1)` each.
//! * `mitm3` (`m = 3`, a 2-sum table): `x(R ± P_k)` looked up in the
//!   stored abscissae of `P_i ± P_j`.
//! * `mitm4` (`m = 4`): `x(R − (P_i ± P_j))` for every pair, looked up in
//!   the same 2-sum table.
//!
//! Two independent references, sharing no code with the oracles' loops:
//! `brute_rational` walks every triple of base points with every sign
//! pattern through the general group law, and `algebraic_count` counts the
//! ordered tuples of `V^m` that solve the summation system through the
//! resultant form the PKM checker uses (`pkm_tower_check`), so that a cell
//! reconstructed from a committed F4 row is checked against the same
//! quantity the row reports.

use crate::fpset::{contains_batch, FpSet};
use crypto_lib::cryptanalysis::prime_fast::{FastCurve, Fp64, Pt};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Ops {
    pub mul: u64,
    pub sqr: u64,
    pub inv: u64,
    pub lookups: u64,
    pub candidates: u64,
}

impl Ops {
    pub fn add(&mut self, o: Ops) {
        self.mul += o.mul;
        self.sqr += o.sqr;
        self.inv += o.inv;
        self.lookups += o.lookups;
        self.candidates += o.candidates;
    }
}

#[inline(always)]
fn hash_key(x: u64) -> u64 {
    x.wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(32)
}

/// Open-addressing set of non-zero `u64` keys (Montgomery abscissae).
pub struct KeySet {
    slots: Vec<u64>,
    mask: usize,
}

impl KeySet {
    pub fn new(keys: &[u64]) -> Self {
        let cap = (keys.len() * 4).max(16).next_power_of_two();
        let mut s = KeySet {
            slots: vec![0u64; cap],
            mask: cap - 1,
        };
        for &k in keys {
            assert!(k != 0, "0 is the empty slot");
            let mut i = (hash_key(k) as usize) & s.mask;
            loop {
                if s.slots[i] == 0 {
                    s.slots[i] = k;
                    break;
                }
                if s.slots[i] == k {
                    break;
                }
                i = (i + 1) & s.mask;
            }
        }
        s
    }
    #[inline(always)]
    pub fn contains(&self, k: u64) -> bool {
        let mut i = (hash_key(k) as usize) & self.mask;
        loop {
            let v = self.slots[i];
            if v == k {
                return true;
            }
            if v == 0 {
                return false;
            }
            i = (i + 1) & self.mask;
        }
    }
}

/// A tower cell: the curve, `V` (Montgomery abscissae) and the base points.
pub struct Cell {
    pub p: u64,
    pub a: u64,
    pub b: u64,
    pub curve: FastCurve,
    /// Abscissae of `V`, Montgomery form.
    pub v: Vec<u64>,
    pub vset: KeySet,
    /// Base points (one of the two lifts), Montgomery form, `x` ascending.
    pub base: Vec<Pt>,
    pub base_x: KeySet,
}

impl Cell {
    pub fn from_v(p: u64, a: u64, b: u64, v_canonical: &[u64]) -> Self {
        let curve = FastCurve::new("tower-cell", p, a, b, p);
        let f = curve.f;
        let mut v: Vec<u64> = v_canonical.iter().map(|&x| f.to_mont(x)).collect();
        v.sort_unstable();
        v.dedup();
        let mut base: Vec<Pt> = v_canonical
            .iter()
            .filter_map(|&x| curve.lift_x(x))
            .collect();
        base.sort_by_key(|pt| pt.x);
        base.dedup_by_key(|pt| pt.x);
        let vset = KeySet::new(&v);
        let base_x = KeySet::new(&base.iter().map(|pt| pt.x).collect::<Vec<_>>());
        Cell {
            p,
            a,
            b,
            curve,
            v,
            vset,
            base,
            base_x,
        }
    }

    /// A Kummer cell from the committed row fields: `V = {g·ζ^e : e < 2^t}`.
    pub fn kummer(p: u64, a: u64, b: u64, t: u32, g: u64, zeta: u64) -> Self {
        let f = Fp64::new(p);
        let gm = f.to_mont(g);
        let zm = f.to_mont(zeta);
        let mut v = Vec::with_capacity(1 << t);
        let mut acc = gm;
        for _ in 0..(1u64 << t) {
            v.push(f.from_mont(acc));
            acc = f.mul(acc, zm);
        }
        Self::from_v(p, a, b, &v)
    }

    pub fn f(&self) -> Fp64 {
        self.curve.f
    }

    pub fn lift(&self, x_canonical: u64) -> Option<Pt> {
        self.curve.lift_x(x_canonical)
    }
}

/// `x(P + Q)` and `x(P − Q)` for `x_P ≠ x_Q`, from `inv = (x_Q − x_P)⁻¹`.
#[inline(always)]
fn two_x(f: &Fp64, p: Pt, q: Pt, inv: u64) -> (u64, u64) {
    let s = f.add(p.x, q.x);
    let lam = f.mul(f.sub(q.y, p.y), inv);
    let xs = f.sub(f.sqr(lam), s);
    let lam2 = f.mul(f.neg(f.add(q.y, p.y)), inv);
    let xd = f.sub(f.sqr(lam2), s);
    (xs, xd)
}

/// The 2-sum table: abscissae of `P_i ± P_j` (`i ≤ j`), sorted, with a
/// presence filter; witnesses `(key, i, j, sign)` when asked for.
pub struct PairSumTable {
    pub keys: Vec<u64>,
    pub witnesses: Option<Vec<(u64, u32, u32, bool)>>,
    pub fp: FpSet,
    pub build_ops: Ops,
}

impl PairSumTable {
    pub fn build(cell: &Cell, with_witnesses: bool) -> Self {
        let f = cell.f();
        let base = &cell.base;
        let m = base.len();
        let mut keys = Vec::with_capacity(m * (m + 1));
        let mut wit = Vec::new();
        let mut ops = Ops::default();
        let mut den = Vec::with_capacity(m);
        let mut scratch = Vec::with_capacity(m);
        for i in 0..m {
            let pi = base[i];
            if let Some(d) = cell.curve.double(pi) {
                ops.inv += 1;
                ops.mul += 3;
                ops.sqr += 2;
                keys.push(d.x);
                if with_witnesses {
                    wit.push((d.x, i as u32, i as u32, true));
                }
            }
            den.clear();
            for pj in &base[i + 1..] {
                den.push(f.sub(pj.x, pi.x));
            }
            if den.is_empty() {
                continue;
            }
            f.batch_inv(&mut den, &mut scratch);
            ops.inv += 1;
            ops.mul += 3 * (den.len() as u64 - 1);
            for (k, pj) in base[i + 1..].iter().enumerate() {
                let (xs, xd) = two_x(&f, pi, *pj, den[k]);
                ops.mul += 2;
                ops.sqr += 2;
                keys.push(xs);
                keys.push(xd);
                if with_witnesses {
                    let j = (i + 1 + k) as u32;
                    wit.push((xs, i as u32, j, true));
                    wit.push((xd, i as u32, j, false));
                }
            }
        }
        keys.sort_unstable();
        let mut fp = FpSet::new(keys.len());
        for &k in &keys {
            fp.insert(k);
        }
        let witnesses = if with_witnesses {
            wit.sort_unstable();
            Some(wit)
        } else {
            None
        };
        PairSumTable {
            keys,
            witnesses,
            fp,
            build_ops: ops,
        }
    }

    #[inline(always)]
    pub fn contains(&self, x: u64) -> bool {
        self.fp.maybe(x) && self.keys.binary_search(&x).is_ok()
    }

    /// Every candidate at once, fingerprint lines prefetched ahead.
    pub fn contains_batch(&self, cands: &[u64], hits: &mut Vec<bool>) {
        contains_batch(&self.fp, &self.keys, cands, hits)
    }

    pub fn bytes(&self) -> usize {
        self.keys.len() * 8 + self.fp.bytes()
    }
}

/// A found decomposition as sorted canonical abscissae of the summands.
pub type Triple = Vec<u64>;

fn sorted_canon(f: &Fp64, xs: &[u64]) -> Triple {
    let mut t: Vec<u64> = xs.iter().map(|&x| f.from_mont(x)).collect();
    t.sort_unstable();
    t
}

/// `m = 2`: `R − P_k` or `R + P_k` in the base.
pub fn subtract2(cell: &Cell, r: Pt, all: bool, out: &mut Vec<Triple>) -> (bool, Ops) {
    let f = cell.f();
    let base = &cell.base;
    let m = base.len();
    let mut ops = Ops::default();
    let mut den: Vec<u64> = base.iter().map(|pk| f.sub(pk.x, r.x)).collect();
    if den.contains(&0) {
        return subtract2_slow(cell, r, all, out, ops);
    }
    let mut scratch = Vec::with_capacity(m);
    f.batch_inv(&mut den, &mut scratch);
    ops.inv += 1;
    ops.mul += 3 * (m as u64 - 1);
    let mut found = false;
    for k in 0..m {
        let (xs, xd) = two_x(&f, r, base[k], den[k]);
        ops.mul += 2;
        ops.sqr += 2;
        ops.lookups += 2;
        ops.candidates += 2;
        for x in [xs, xd] {
            if cell.base_x.contains(x) {
                found = true;
                out.push(sorted_canon(&f, &[base[k].x, x]));
                if !all {
                    return (true, ops);
                }
            }
        }
    }
    (found, ops)
}

#[cold]
fn subtract2_slow(
    cell: &Cell,
    r: Pt,
    all: bool,
    out: &mut Vec<Triple>,
    mut ops: Ops,
) -> (bool, Ops) {
    let f = cell.f();
    let mut found = false;
    for pk in &cell.base {
        for q in [
            cell.curve.add(Some(r), Some(cell.curve.neg(*pk))),
            cell.curve.add(Some(r), Some(*pk)),
        ] {
            ops.inv += 1;
            ops.candidates += 1;
            if let Some(q) = q {
                if cell.base_x.contains(q.x) {
                    found = true;
                    out.push(sorted_canon(&f, &[pk.x, q.x]));
                    if !all {
                        return (true, ops);
                    }
                }
            }
        }
    }
    (found, ops)
}

/// `m = 3`, memoryless pairs-and-solve by the group law.
pub fn pairs3(cell: &Cell, r: Pt, all: bool, out: &mut Vec<Triple>) -> (bool, Ops) {
    let f = cell.f();
    let base = &cell.base;
    let m = base.len();
    let mut ops = Ops::default();
    // Row points U± = R ± P₁ for every P₁, one batched inversion.
    let mut den: Vec<u64> = base.iter().map(|p1| f.sub(p1.x, r.x)).collect();
    let mut scratch = Vec::with_capacity(2 * m);
    let degenerate_rows = den.contains(&0);
    let mut rows: Vec<Option<(Pt, Pt)>> = Vec::with_capacity(m);
    if degenerate_rows {
        for p1 in base {
            let up = cell.curve.add(Some(r), Some(*p1));
            let um = cell.curve.add(Some(r), Some(cell.curve.neg(*p1)));
            ops.inv += 2;
            rows.push(match (up, um) {
                (Some(a), Some(b)) => Some((a, b)),
                _ => None,
            });
        }
    } else {
        f.batch_inv(&mut den, &mut scratch);
        ops.inv += 1;
        ops.mul += 3 * (m as u64 - 1);
        for (k, p1) in base.iter().enumerate() {
            let inv = den[k];
            let up = cell.curve.add_with_inv(r, *p1, inv);
            let um = cell.curve.add_with_inv(r, cell.curve.neg(*p1), inv);
            ops.mul += 4;
            ops.sqr += 2;
            rows.push(Some((up, um)));
        }
    }
    let mut found = false;
    let mut den2: Vec<u64> = Vec::with_capacity(2 * m);
    for i1 in 0..m {
        let Some((up, um)) = rows[i1] else {
            // R = ∓P₁: a two-summand relation; the three-summand search
            // for this row is done exactly through the general group law.
            found |= row_slow(cell, r, i1, all, out, &mut ops);
            if found && !all {
                return (true, ops);
            }
            continue;
        };
        let cnt = m - i1;
        den2.clear();
        for p2 in &base[i1..] {
            den2.push(f.sub(p2.x, up.x));
        }
        for p2 in &base[i1..] {
            den2.push(f.sub(p2.x, um.x));
        }
        if den2.contains(&0) {
            found |= row_slow(cell, r, i1, all, out, &mut ops);
            if found && !all {
                return (true, ops);
            }
            continue;
        }
        f.batch_inv(&mut den2, &mut scratch);
        ops.inv += 1;
        ops.mul += 3 * (2 * cnt as u64 - 1);
        for k in 0..cnt {
            let p2 = base[i1 + k];
            let (t1, t2) = two_x(&f, up, p2, den2[k]);
            let (t3, t4) = two_x(&f, um, p2, den2[cnt + k]);
            ops.mul += 4;
            ops.sqr += 4;
            ops.candidates += 4;
            let h1 = cell.vset.contains(t1);
            let h2 = cell.vset.contains(t2);
            let h3 = cell.vset.contains(t3);
            let h4 = cell.vset.contains(t4);
            if h1 | h2 | h3 | h4 {
                for (hit, t) in [(h1, t1), (h2, t2), (h3, t3), (h4, t4)] {
                    if hit {
                        found = true;
                        out.push(sorted_canon(&f, &[base[i1].x, p2.x, t]));
                        if !all {
                            return (true, ops);
                        }
                    }
                }
            }
        }
    }
    (found, ops)
}

#[cold]
fn row_slow(
    cell: &Cell,
    r: Pt,
    i1: usize,
    all: bool,
    out: &mut Vec<Triple>,
    ops: &mut Ops,
) -> bool {
    let f = cell.f();
    let base = &cell.base;
    let p1 = base[i1];
    let mut found = false;
    for p2 in &base[i1..] {
        for s1 in [1i8, -1] {
            for s2 in [1i8, -1] {
                let a = if s1 > 0 { p1 } else { cell.curve.neg(p1) };
                let b = if s2 > 0 { *p2 } else { cell.curve.neg(*p2) };
                let q = cell.curve.add(
                    cell.curve.add(Some(r), Some(cell.curve.neg(a))),
                    Some(cell.curve.neg(b)),
                );
                ops.inv += 2;
                ops.candidates += 1;
                if let Some(q) = q {
                    if cell.vset.contains(q.x) {
                        found = true;
                        out.push(sorted_canon(&f, &[p1.x, p2.x, q.x]));
                        if !all {
                            return true;
                        }
                    }
                }
            }
        }
    }
    found
}

/// `m = 3` through the 2-sum table.
pub fn mitm3(
    cell: &Cell,
    table: &PairSumTable,
    r: Pt,
    all: bool,
    out: &mut Vec<Triple>,
) -> (bool, Ops) {
    let f = cell.f();
    let base = &cell.base;
    let m = base.len();
    let mut ops = Ops::default();
    let mut den: Vec<u64> = base.iter().map(|pk| f.sub(pk.x, r.x)).collect();
    if den.contains(&0) {
        return mitm3_slow(cell, table, r, all, out, ops);
    }
    let mut scratch = Vec::with_capacity(m);
    f.batch_inv(&mut den, &mut scratch);
    ops.inv += 1;
    ops.mul += 3 * (m as u64 - 1);
    // All 2|F| candidates first, then the probes in one prefetched pass.
    let mut cands: Vec<u64> = Vec::with_capacity(2 * m);
    for k in 0..m {
        let (xs, xd) = two_x(&f, r, base[k], den[k]);
        cands.push(xs);
        cands.push(xd);
    }
    ops.mul += 2 * m as u64;
    ops.sqr += 2 * m as u64;
    ops.lookups += 2 * m as u64;
    ops.candidates += 2 * m as u64;
    let mut hits: Vec<bool> = Vec::with_capacity(2 * m);
    table.contains_batch(&cands, &mut hits);
    let mut found = false;
    for (i, &h) in hits.iter().enumerate() {
        if h {
            found = true;
            if let Some(w) = &table.witnesses {
                push_witnesses(cell, w, cands[i], base[i / 2].x, out);
            }
            if !all {
                return (true, ops);
            }
        }
    }
    (found, ops)
}

#[cold]
fn mitm3_slow(
    cell: &Cell,
    table: &PairSumTable,
    r: Pt,
    all: bool,
    out: &mut Vec<Triple>,
    mut ops: Ops,
) -> (bool, Ops) {
    let mut found = false;
    for pk in &cell.base {
        for q in [
            cell.curve.add(Some(r), Some(cell.curve.neg(*pk))),
            cell.curve.add(Some(r), Some(*pk)),
        ] {
            ops.inv += 1;
            ops.lookups += 1;
            if let Some(q) = q {
                if table.contains(q.x) {
                    found = true;
                    if let Some(w) = &table.witnesses {
                        push_witnesses(cell, w, q.x, pk.x, out);
                    }
                    if !all {
                        return (true, ops);
                    }
                }
            }
        }
    }
    (found, ops)
}

fn push_witnesses(
    cell: &Cell,
    w: &[(u64, u32, u32, bool)],
    key: u64,
    xk: u64,
    out: &mut Vec<Triple>,
) {
    let f = cell.f();
    let lo = w.partition_point(|e| e.0 < key);
    let mut idx = lo;
    while idx < w.len() && w[idx].0 == key {
        let (_, i, j, _) = w[idx];
        out.push(sorted_canon(
            &f,
            &[cell.base[i as usize].x, cell.base[j as usize].x, xk],
        ));
        idx += 1;
    }
}

/// `m = 4` through the 2-sum table: `R − (P_i ± P_j)` for every pair.
pub fn mitm4(
    cell: &Cell,
    table: &PairSumTable,
    r: Pt,
    all: bool,
    out: &mut Vec<Triple>,
) -> (bool, Ops) {
    let f = cell.f();
    let base = &cell.base;
    let m = base.len();
    let mut ops = Ops::default();
    let mut found = false;
    let mut den = Vec::with_capacity(m);
    let mut scratch = Vec::with_capacity(m);
    for i in 0..m {
        let pi = base[i];
        // Pair sums S = P_i ± P_j for j ≥ i (j = i: 2P_i), then R − S and R + S.
        let mut sums: Vec<(Pt, usize, bool)> = Vec::with_capacity(2 * (m - i));
        if let Some(d) = cell.curve.double(pi) {
            ops.inv += 1;
            ops.mul += 3;
            ops.sqr += 2;
            sums.push((d, i, true));
        }
        den.clear();
        for pj in &base[i + 1..] {
            den.push(f.sub(pj.x, pi.x));
        }
        if !den.is_empty() {
            f.batch_inv(&mut den, &mut scratch);
            ops.inv += 1;
            ops.mul += 3 * (den.len() as u64 - 1);
            for (k, pj) in base[i + 1..].iter().enumerate() {
                let inv = den[k];
                let s = cell.curve.add_with_inv(pi, *pj, inv);
                let d = cell.curve.add_with_inv(pi, cell.curve.neg(*pj), inv);
                ops.mul += 4;
                ops.sqr += 2;
                sums.push((s, i + 1 + k, true));
                sums.push((d, i + 1 + k, false));
            }
        }
        // R ∓ S with one batched inversion over the sums.
        den.clear();
        for (s, _, _) in &sums {
            den.push(f.sub(s.x, r.x));
        }
        if den.contains(&0) {
            for (s, j, sign) in &sums {
                for q in [
                    cell.curve.add(Some(r), Some(cell.curve.neg(*s))),
                    cell.curve.add(Some(r), Some(*s)),
                ] {
                    ops.inv += 1;
                    ops.lookups += 1;
                    ops.candidates += 1;
                    if let Some(q) = q {
                        if table.contains(q.x) {
                            found = true;
                            if let Some(w) = &table.witnesses {
                                push_witnesses4(cell, w, q.x, base[i].x, base[*j].x, *sign, out);
                            }
                            if !all {
                                return (true, ops);
                            }
                        }
                    }
                }
            }
            continue;
        }
        f.batch_inv(&mut den, &mut scratch);
        ops.inv += 1;
        ops.mul += 3 * (den.len() as u64).saturating_sub(1);
        let mut cands: Vec<u64> = Vec::with_capacity(2 * sums.len());
        for (k, (s, _, _)) in sums.iter().enumerate() {
            let (xs, xd) = two_x(&f, r, *s, den[k]);
            cands.push(xs);
            cands.push(xd);
        }
        let ns = sums.len() as u64;
        ops.mul += 2 * ns;
        ops.sqr += 2 * ns;
        ops.lookups += 2 * ns;
        ops.candidates += 2 * ns;
        let mut hits: Vec<bool> = Vec::with_capacity(2 * sums.len());
        table.contains_batch(&cands, &mut hits);
        for (c, &h) in hits.iter().enumerate() {
            if h {
                found = true;
                let (_, j, sign) = sums[c / 2];
                if let Some(w) = &table.witnesses {
                    push_witnesses4(cell, w, cands[c], base[i].x, base[j].x, sign, out);
                }
                if !all {
                    return (true, ops);
                }
            }
        }
    }
    (found, ops)
}

fn push_witnesses4(
    cell: &Cell,
    w: &[(u64, u32, u32, bool)],
    key: u64,
    xi: u64,
    xj: u64,
    _sign: bool,
    out: &mut Vec<Triple>,
) {
    let f = cell.f();
    let lo = w.partition_point(|e| e.0 < key);
    let mut idx = lo;
    while idx < w.len() && w[idx].0 == key {
        let (_, a, b, _) = w[idx];
        out.push(sorted_canon(
            &f,
            &[xi, xj, cell.base[a as usize].x, cell.base[b as usize].x],
        ));
        idx += 1;
    }
}

/// Independent reference: every `m`-multiset of base points with every sign
/// pattern through the general group law; sorted canonical abscissae.
pub fn brute_rational(cell: &Cell, r: Pt, m: usize) -> Vec<Triple> {
    let f = cell.f();
    let base = &cell.base;
    let n = base.len();
    let neg_r = cell.curve.neg(r);
    let mut out = Vec::new();
    let mut idx = vec![0usize; m];
    // Enumerate non-decreasing index tuples.
    fn next(idx: &mut [usize], n: usize) -> bool {
        let m = idx.len();
        let mut k = m;
        while k > 0 {
            k -= 1;
            if idx[k] + 1 < n {
                idx[k] += 1;
                for j in k + 1..m {
                    idx[j] = idx[k];
                }
                return true;
            }
        }
        false
    }
    loop {
        for signs in 0..(1u32 << (m - 1)) {
            let mut acc = Some(base[idx[0]]);
            for (s, &i) in idx[1..].iter().enumerate() {
                let p = if (signs >> s) & 1 == 1 {
                    cell.curve.neg(base[i])
                } else {
                    base[i]
                };
                acc = cell.curve.add(acc, Some(p));
            }
            if let Some(q) = acc {
                if q == r || q == neg_r {
                    let xs: Vec<u64> = idx.iter().map(|&i| base[i].x).collect();
                    out.push(sorted_canon(&f, &xs));
                }
            }
        }
        if !next(&mut idx, n) {
            break;
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

// ── the PKM checker's algebraic count (resultant form), canonical arithmetic ──

#[derive(Clone, Copy)]
struct Fc(u64);
impl Fc {
    fn add(self, a: u64, b: u64) -> u64 {
        (a + b) % self.0
    }
    fn sub(self, a: u64, b: u64) -> u64 {
        (a + self.0 - b % self.0) % self.0
    }
    fn mul(self, a: u64, b: u64) -> u64 {
        (u128::from(a) * u128::from(b) % u128::from(self.0)) as u64
    }
    fn neg(self, a: u64) -> u64 {
        self.sub(0, a)
    }
}

fn s3_coeffs(f: Fc, x1: u64, x2: u64, a: u64, b: u64) -> (u64, u64, u64) {
    let d = f.sub(x1, x2);
    let s = f.add(x1, x2);
    let pr = f.mul(x1, x2);
    let aa = f.mul(d, d);
    let bb = f.neg(f.mul(2, f.add(f.mul(s, f.add(pr, a)), f.mul(2, b))));
    let q = f.sub(pr, a);
    let cc = f.sub(f.mul(q, q), f.mul(f.mul(4, b), s));
    (aa, bb, cc)
}

fn s3(f: Fc, x1: u64, x2: u64, x3: u64, a: u64, b: u64) -> u64 {
    let (aa, bb, cc) = s3_coeffs(f, x1, x2, a, b);
    f.add(f.add(f.mul(f.mul(aa, x3), x3), f.mul(bb, x3)), cc)
}

fn res2(f: Fc, g: (u64, u64, u64), h: (u64, u64, u64)) -> u64 {
    let (a, b, c) = g;
    let (a2, b2, c2) = h;
    let u = f.sub(f.mul(a, c2), f.mul(a2, c));
    let v = f.sub(f.mul(a, b2), f.mul(a2, b));
    let w = f.sub(f.mul(b, c2), f.mul(b2, c));
    f.sub(f.mul(u, u), f.mul(v, w))
}

/// Ordered tuples of `V^m` (canonical abscissae) solving the summation
/// system of the target `x_r`, as `pkm_tower_check` counts them (`m ≤ 3`).
pub fn algebraic_count(cell: &Cell, x_r: u64, m: usize) -> Option<u64> {
    let f = Fc(cell.p);
    let fm = cell.f();
    let v: Vec<u64> = cell.v.iter().map(|&x| fm.from_mont(x)).collect();
    let (a, b) = (cell.a, cell.b);
    match m {
        2 => Some(
            v.iter()
                .flat_map(|&x1| v.iter().map(move |&x2| (x1, x2)))
                .filter(|&(x1, x2)| s3(f, x1, x2, x_r, a, b) == 0)
                .count() as u64,
        ),
        3 => {
            let tails: Vec<(u64, u64, u64)> =
                v.iter().map(|&x3| s3_coeffs(f, x3, x_r, a, b)).collect();
            let mut n = 0u64;
            for &x1 in &v {
                for &x2 in &v {
                    let g = s3_coeffs(f, x1, x2, a, b);
                    n += tails.iter().filter(|&&h| res2(f, g, h) == 0).count() as u64;
                }
            }
            Some(n)
        }
        _ => None,
    }
}
