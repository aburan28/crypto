//! The pair-table oracle: ground truth for the verdicts and the solution
//! sets, organised differently from the kernel.
//!
//! The panel's exact oracle decides `R = p + a + b` over base points by
//! looking `x(R − p)` up in a table of the abscissae of every pair sum
//! `a + b`.  Its code runs on the library's slow inversion and is out of
//! reach for the larger bases (see the README), so this module computes the same
//! predicate with the same meet-in-the-middle structure but one pass over the
//! pair sums for all targets at once: the keys are `x(R ± p)` for every
//! target and base point, and a pair sum `a ± b` (no row hoisting of the target)
//! that lands on a key is a decomposition.  It also reproduces the panel's
//! `ref_count` (the raw sum of the table multiplicities), which is a
//! fingerprint of the whole pair table the panel built.

use crate::base::Base;
use crate::curve::{add, neg, sorted3, Pt};
use crate::field::Onb;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

#[derive(Default)]
pub struct Fx(u64);

impl Hasher for Fx {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(5) ^ b as u64).wrapping_mul(0x517c_c1b7_2722_0a95);
        }
    }
    fn write_u64(&mut self, x: u64) {
        self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(0x517c_c1b7_2722_0a95);
        self.0 ^= self.0 >> 29;
    }
}

pub type FxMap<V> = HashMap<u64, V, BuildHasherDefault<Fx>>;
pub type FxSet = std::collections::HashSet<u64, BuildHasherDefault<Fx>>;

pub struct O3 {
    /// For each target: every sorted triple of abscissae with `±P_a ± P_b ± P_c = R`.
    pub sols: Vec<Vec<[u64; 3]>>,
    /// For each target: the panel's `ref_count` (raw).
    pub ref_count: Vec<u64>,
    pub keys: usize,
    pub pair_sums: u64,
}

struct Acc<'a> {
    keys: &'a FxMap<Vec<(u16, u32)>>,
    b: &'a Base,
    sols: Vec<Vec<[u64; 3]>>,
    ref_count: Vec<u64>,
}

impl Acc<'_> {
    #[inline]
    fn visit(&mut self, x: u64, i: usize, j: usize, weight: u64) {
        if let Some(list) = self.keys.get(&x) {
            for &(t, p) in list {
                self.sols[t as usize].push(sorted3(
                    self.b.xs[p as usize],
                    self.b.xs[i],
                    self.b.xs[j],
                ));
                self.ref_count[t as usize] += weight;
            }
        }
    }
}

/// `targets` are `(x_R, y_R)`.
pub fn o3(f: &Onb, b: &Base, targets: &[(u64, u64)]) -> O3 {
    let m = b.xs.len();
    // keys: x(R + s·p) -> (target, p index).  The point x = 0 is its own negative.
    let mut keys: FxMap<Vec<(u16, u32)>> = FxMap::default();
    let mut nkeys = 0usize;
    for (t, &(xr, yr)) in targets.iter().enumerate() {
        let r = Pt::Aff(xr, yr);
        for pi in 0..m {
            let p = Pt::Aff(b.xs[pi], b.ys[pi]);
            let list: Vec<Pt> = if b.xs[pi] == 0 {
                vec![p]
            } else {
                vec![p, neg(p)]
            };
            for q in list {
                if let Pt::Aff(qx, _) = add(f, r, q) {
                    keys.entry(qx).or_default().push((t as u16, pi as u32));
                    nkeys += 1;
                }
            }
        }
    }
    let mut acc = Acc {
        keys: &keys,
        b,
        sols: vec![Vec::new(); targets.len()],
        ref_count: vec![0u64; targets.len()],
    };
    let mut den: Vec<u64> = Vec::with_capacity(m);
    let mut scratch: Vec<u64> = Vec::with_capacity(m);
    let mut pair_sums = 0u64;
    for i in 0..m {
        // a_i + a_i: a candidate for the solutions, not a pair of distinct points.
        if b.xs[i] != 0 {
            let lam = b.xs[i] ^ f.mul_raw(b.ys[i], f.inv_c(b.xs[i]));
            let x2 = f.sqr_raw(lam) ^ lam;
            pair_sums += 1;
            acc.visit(x2, i, i, 0);
        }
        den.clear();
        den.extend((i + 1..m).map(|j| b.xs[i] ^ b.xs[j]));
        f.batch_inv(&mut den, &mut scratch);
        for (k, &invd) in den.iter().enumerate() {
            let j = i + 1 + k;
            let d = b.xs[i] ^ b.xs[j];
            let lam = f.mul_raw(b.ys[i] ^ b.ys[j], invd);
            let xa = f.sqr_raw(lam) ^ lam ^ d;
            let lam2 = lam ^ f.mul_raw(b.xs[j], invd);
            let xb = f.sqr_raw(lam2) ^ lam2 ^ d;
            // Distinct base points: a_i, a_j with both signs; x = 0 is one point.
            let wgt = if b.xs[i] == 0 || b.xs[j] == 0 { 1 } else { 2 };
            pair_sums += 2;
            acc.visit(xa, i, j, wgt);
            acc.visit(xb, i, j, wgt);
        }
    }
    let Acc {
        mut sols,
        ref_count,
        ..
    } = acc;
    for s in sols.iter_mut() {
        s.sort_unstable();
        s.dedup();
    }
    O3 {
        sols,
        ref_count,
        keys: nkeys,
        pair_sums,
    }
}
