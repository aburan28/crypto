//! Binary-track root-finders for the PDP (`m = 3`) on the PS-0 rungs:
//! the pair-sum table (`rf-mitm`) and a tuned pairs-and-solve (`rf-pairs`).
//!
//! Both answer the PS-0 question exactly — every sorted triple
//! `(x₁ ≤ x₂ ≤ x₃)` of base abscissae with `±P₁ ± P₂ ± P₃ = ±R` — by the
//! group law on the one-word field, never by a polynomial system.
//!
//! * `rf-mitm`: `prepare` stores `x(P_i + P_j)` and `x(P_i − P_j)` for every
//!   `i ≤ j` (the 2-sums), then `decide(R)` asks for each base point `P_k`
//!   whether `x(R − P_k)` or `x(R + P_k)` is a stored 2-sum: `|F|` batched
//!   group additions and `2|F|` lookups a target, against the `|F|²/2`
//!   pairs of pairs-and-solve.  A presence filter (bitset over a hash of
//!   the key) answers the misses, which is almost all lookups, without
//!   touching the sorted key array.
//! * `rf-pairs`: the PS-2 row kernel with the denominators of both row
//!   points shared (one batched inversion of `x_R + x_k` serves `R ± P_k`
//!   when the row point is formed) and the inner loop kept branch-free
//!   until a candidate lands in `V`.
//!
//! Counted units (`Ops`): field multiplications, squarings, inversions
//! (the one inversion inside each batched inversion; its `3(k − 1)`
//! multiplications are counted as multiplications), table lookups and
//! candidate abscissae tested.

use crate::fpset::{contains_batch, FpSet};
use crate::group::*;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;

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

/// The 2-sum table of a base: the abscissae of `P_i ± P_j`, `i ≤ j`.
pub struct PairSumTable {
    /// Sorted abscissae of every stored 2-sum (duplicates kept).
    pub keys: Vec<u64>,
    /// `(key, i, j, sign)` sorted by key, only when witnesses are wanted.
    pub witnesses: Option<Vec<(u64, u32, u32, bool)>>,
    pub fp: FpSet,
    pub build_ops: Ops,
}

impl PairSumTable {
    /// Build from the base points (`x` ascending), with one batched
    /// inversion per row `i` over `x_i + x_j`, `j > i`; the doubling `2P_i`
    /// is one inversion of its own.
    pub fn build(gf: &Gf2, base: &[(u64, u64)], with_witnesses: bool) -> Self {
        let m = base.len();
        let cap = m * (m + 1);
        let mut keys: Vec<u64> = Vec::with_capacity(cap);
        let mut wit: Vec<(u64, u32, u32, bool)> = if with_witnesses {
            Vec::with_capacity(cap)
        } else {
            Vec::new()
        };
        let mut ops = Ops::default();
        let mut den: Vec<u64> = Vec::with_capacity(m);
        let mut scratch: Vec<u64> = Vec::with_capacity(m);
        for i in 0..m {
            let (xi, yi) = base[i];
            // 2P_i (x_i ≠ 0 on the base: x = 0 is the 2-torsion point).
            if let Pt::Aff(x2, _) = add(gf, Pt::Aff(xi, yi), Pt::Aff(xi, yi)) {
                ops.inv += 1;
                ops.mul += 2;
                ops.sqr += 2;
                keys.push(x2);
                if with_witnesses {
                    wit.push((x2, i as u32, i as u32, true));
                }
            }
            den.clear();
            for &(xj, _) in &base[i + 1..] {
                den.push(xi ^ xj);
            }
            if den.is_empty() {
                continue;
            }
            gf.batch_inv(&mut den, &mut scratch);
            ops.inv += 1;
            ops.mul += 3 * (den.len() as u64 - 1);
            for (k, &(xj, yj)) in base[i + 1..].iter().enumerate() {
                let j = i + 1 + k;
                let inv = den[k];
                let d = xi ^ xj;
                // P_i + P_j: λ = (y_i + y_j)/(x_i + x_j); x = λ² + λ + d + 1.
                let lam = gf.mul(yi ^ yj, inv);
                let xs = gf.sqr(lam) ^ lam ^ d ^ 1;
                // P_i − P_j: λ' = λ + x_j/(x_i + x_j).
                let lam2 = lam ^ gf.mul(xj, inv);
                let xd = gf.sqr(lam2) ^ lam2 ^ d ^ 1;
                ops.mul += 2;
                ops.sqr += 2;
                keys.push(xs);
                keys.push(xd);
                if with_witnesses {
                    wit.push((xs, i as u32, j as u32, true));
                    wit.push((xd, i as u32, j as u32, false));
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

    /// Every candidate at once, fingerprint lines prefetched ahead.
    pub fn contains_batch(&self, cands: &[u64], hits: &mut Vec<bool>) {
        contains_batch(&self.fp, &self.keys, cands, hits)
    }

    pub fn bytes(&self) -> usize {
        self.keys.len() * 8 + self.fp.bytes()
    }
}

/// The `rf-mitm` arm.
pub struct Mitm<'a> {
    pub gf: &'a Gf2,
    pub base: Vec<(u64, u64)>,
    pub table: PairSumTable,
}

impl<'a> Mitm<'a> {
    pub fn new(gf: &'a Gf2, solver: &ArtinSchreier, l: u32, with_witnesses: bool) -> Self {
        let base = base_points(gf, solver, l);
        let table = PairSumTable::build(gf, &base, with_witnesses);
        Mitm { gf, base, table }
    }

    /// Decide `R`: the first hit stops unless `all`; every hit's triples are
    /// pushed when witnesses are stored.  Returns (found, ops).
    pub fn run(&self, r: (u64, u64), all: bool, out: &mut Vec<[u64; 3]>) -> (bool, Ops) {
        #[cfg(target_arch = "x86_64")]
        if self.gf.folds() && std::arch::is_x86_feature_detected!("pclmulqdq") {
            // SAFETY: pclmulqdq was just detected.
            return unsafe { self.run_pclmul(r, all, out) };
        }
        self.run_body(r, all, out)
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "pclmulqdq")]
    unsafe fn run_pclmul(&self, r: (u64, u64), all: bool, out: &mut Vec<[u64; 3]>) -> (bool, Ops) {
        self.run_body(r, all, out)
    }

    #[inline(always)]
    fn run_body(&self, r: (u64, u64), all: bool, out: &mut Vec<[u64; 3]>) -> (bool, Ops) {
        let gf = self.gf;
        let (xr, yr) = r;
        let base = &self.base;
        let m = base.len();
        let mut ops = Ops::default();
        let mut den: Vec<u64> = Vec::with_capacity(m);
        let mut scratch: Vec<u64> = Vec::with_capacity(m);
        for &(xk, _) in base {
            den.push(xr ^ xk);
        }
        // x_R ∉ V for the targets routed here, so every denominator is non-zero.
        gf.batch_inv(&mut den, &mut scratch);
        ops.inv += 1;
        ops.mul += 3 * (m as u64 - 1);
        // All 2|F| candidates first, then the probes in one prefetched pass.
        let mut cands: Vec<u64> = Vec::with_capacity(2 * m);
        for k in 0..m {
            let (xk, yk) = base[k];
            let inv = den[k];
            let d = xr ^ xk;
            // R + P_k and R − P_k share the inversion.
            let lam = gf.mul(yr ^ yk, inv);
            let xs = gf.sqr(lam) ^ lam ^ d ^ 1;
            let lam2 = lam ^ gf.mul(xk, inv);
            let xd = gf.sqr(lam2) ^ lam2 ^ d ^ 1;
            cands.push(xs);
            cands.push(xd);
        }
        ops.mul += 2 * m as u64;
        ops.sqr += 2 * m as u64;
        ops.lookups += 2 * m as u64;
        ops.candidates += 2 * m as u64;
        let mut hits: Vec<bool> = Vec::with_capacity(2 * m);
        self.table.contains_batch(&cands, &mut hits);
        let mut found = false;
        for (i, &h) in hits.iter().enumerate() {
            if h {
                found = true;
                if let Some(w) = &self.table.witnesses {
                    self.push_witnesses(w, cands[i], base[i / 2].0, out);
                }
                if !all {
                    return (true, ops);
                }
            }
        }
        (found, ops)
    }

    fn push_witnesses(
        &self,
        w: &[(u64, u32, u32, bool)],
        key: u64,
        xk: u64,
        out: &mut Vec<[u64; 3]>,
    ) {
        let lo = w.partition_point(|e| e.0 < key);
        let mut idx = lo;
        while idx < w.len() && w[idx].0 == key {
            let (_, i, j, _) = w[idx];
            let (xi, xj) = (self.base[i as usize].0, self.base[j as usize].0);
            let mut t = [xi, xj, xk];
            t.sort_unstable();
            out.push(t);
            idx += 1;
        }
    }
}

/// The `rf-pairs` arm: pairs-and-solve with shared row inversions.
pub struct Pairs<'a> {
    pub gf: &'a Gf2,
    pub l: u32,
    pub base: Vec<(u64, u64)>,
}

impl<'a> Pairs<'a> {
    pub fn new(gf: &'a Gf2, solver: &ArtinSchreier, l: u32) -> Self {
        Pairs {
            gf,
            l,
            base: base_points(gf, solver, l),
        }
    }

    pub fn run(&self, r: (u64, u64), all: bool, out: &mut Vec<[u64; 3]>) -> (bool, Ops) {
        #[cfg(target_arch = "x86_64")]
        if self.gf.folds() && std::arch::is_x86_feature_detected!("pclmulqdq") {
            // SAFETY: pclmulqdq was just detected.
            return unsafe { self.run_pclmul(r, all, out) };
        }
        self.run_body(r, all, out)
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "pclmulqdq")]
    unsafe fn run_pclmul(&self, r: (u64, u64), all: bool, out: &mut Vec<[u64; 3]>) -> (bool, Ops) {
        self.run_body(r, all, out)
    }

    #[inline(always)]
    fn run_body(&self, r: (u64, u64), all: bool, out: &mut Vec<[u64; 3]>) -> (bool, Ops) {
        let gf = self.gf;
        let span = 1u64 << self.l;
        let base = &self.base;
        let m = base.len();
        let (xr, yr) = r;
        let mut ops = Ops::default();
        // Row points U± = R ± P_k for every k, from one batched inversion.
        let mut den: Vec<u64> = base.iter().map(|&(xk, _)| xr ^ xk).collect();
        let mut scratch: Vec<u64> = Vec::with_capacity(2 * m);
        gf.batch_inv(&mut den, &mut scratch);
        ops.inv += 1;
        ops.mul += 3 * (m as u64 - 1);
        let mut rows: Vec<(u64, u64, u64, u64)> = Vec::with_capacity(m);
        for k in 0..m {
            let (xk, yk) = base[k];
            let inv = den[k];
            let d = xr ^ xk;
            let lam = gf.mul(yr ^ yk, inv);
            let xu = gf.sqr(lam) ^ lam ^ d ^ 1;
            let yu = gf.mul(lam, xr ^ xu) ^ xu ^ yr;
            let lam2 = lam ^ gf.mul(xk, inv);
            let xv = gf.sqr(lam2) ^ lam2 ^ d ^ 1;
            let yv = gf.mul(lam2, xr ^ xv) ^ xv ^ yr;
            ops.mul += 4;
            ops.sqr += 2;
            rows.push((xu, yu, xv, yv));
        }
        let mut found = false;
        let mut den2: Vec<u64> = Vec::with_capacity(2 * m);
        for i1 in 0..m {
            let (xu, yu, xv, yv) = rows[i1];
            let cnt = m - i1;
            den2.clear();
            for &(x2, _) in &base[i1..] {
                den2.push(xu ^ x2);
            }
            for &(x2, _) in &base[i1..] {
                den2.push(xv ^ x2);
            }
            if den2.contains(&0) {
                // U = ±P₂ for some P₂ (x_U ∈ V): the exact slow path for the row.
                for i2 in i1..m {
                    ops.candidates += 4;
                    for t in thirds(gf, base[i1], base[i2], r).into_iter().flatten() {
                        if t < span {
                            found = true;
                            out.push(sorted3(base[i1].0, base[i2].0, t));
                            if !all {
                                return (true, ops);
                            }
                        }
                    }
                }
                continue;
            }
            gf.batch_inv(&mut den2, &mut scratch);
            ops.inv += 1;
            ops.mul += 3 * (2 * cnt as u64 - 1);
            for k in 0..cnt {
                let (x2, y2) = base[i1 + k];
                let (j, jj) = (den2[k], den2[cnt + k]);
                let d = xu ^ x2;
                let lam = gf.mul(yu ^ y2, j);
                let t1 = gf.sqr(lam) ^ lam ^ d ^ 1;
                let lam2 = lam ^ gf.mul(x2, j);
                let t2 = gf.sqr(lam2) ^ lam2 ^ d ^ 1;
                let dd = xv ^ x2;
                let mu = gf.mul(yv ^ y2, jj);
                let t3 = gf.sqr(mu) ^ mu ^ dd ^ 1;
                let mu2 = mu ^ gf.mul(x2, jj);
                let t4 = gf.sqr(mu2) ^ mu2 ^ dd ^ 1;
                ops.mul += 4;
                ops.sqr += 4;
                ops.candidates += 4;
                if (t1 < span) | (t2 < span) | (t3 < span) | (t4 < span) {
                    for t in [t1, t2, t3, t4] {
                        if t < span {
                            found = true;
                            out.push(sorted3(base[i1].0, x2, t));
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
}

pub fn sorted3(a: u64, b: u64, c: u64) -> [u64; 3] {
    let mut t = [a, b, c];
    t.sort_unstable();
    t
}
