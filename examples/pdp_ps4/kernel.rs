//! The scalar group-law pair kernel over the normal-basis field.
//!
//! For a target `R` and a first summand `P₁` let `U₊ = P₁ + R`, `U₋ = P₁ − R`
//! (two exact additions per row).  For every second summand `P₂` the four
//! candidate thirds are `x(U₊ ± P₂)` and `x(U₋ ± P₂)`; with `a = 0` and
//! `λ = (y_U + y₂)/(x_U + x₂)`, `d = x_U + x₂`
//!
//! ```text
//!   x(U + P₂) = λ² + λ + d,        x(U − P₂) = λ̃² + λ̃ + d,   λ̃ = λ + x₂/d
//! ```
//!
//! The denominators of a row (`2` per pair) are inverted together by
//! Montgomery's trick.  The leaf is a closure: `wt(t) ≤ w` (one popcount) for the
//! Hamming-weight base.  Squaring is a free coordinate permutation here.  The
//! row is the unit of the loop: a loop is a set of first summands and, for
//! each, either every base point as second summand or those from the first
//! one on (the PS-0/PS-2 convention, `m(m+1)/2` pairs).

use crate::base::Base;
use crate::curve::{add, neg, sorted3, thirds, Pt};
use crate::field::{Counts, Onb};
use std::time::Instant;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Second {
    /// Every base point (`m` pairs per row).
    All,
    /// From the first summand on (`m − i₁` pairs per row).
    FromFirst,
}

pub struct LoopResult {
    /// Pairs `(P₁, P₂)` visited.
    pub pairs: u64,
    /// Candidate thirds tested by the leaf.
    pub thirds: u64,
    /// Pairs that went through the exact (slow) path.
    pub exact_pairs: u64,
    pub rows: u64,
    /// Sorted triples of abscissae, with repeats (one per candidate hit).
    pub hits: Vec<[u64; 3]>,
    pub found: bool,
    pub counts: Counts,
    pub nanos: u128,
}

#[allow(clippy::too_many_arguments)]
fn exact_pair<L: Fn(u64) -> bool>(
    f: &Onb,
    b: &Base,
    i1: usize,
    i2: usize,
    r: (u64, u64),
    leaf: &L,
    res: &mut LoopResult,
) {
    res.exact_pairs += 1;
    let p1 = (b.xs[i1], b.ys[i1]);
    let p2 = (b.xs[i2], b.ys[i2]);
    for t in thirds(f, p1, p2, r).into_iter().flatten() {
        res.thirds += 1;
        if leaf(t) {
            res.found = true;
            res.hits.push(sorted3(p1.0, p2.0, t));
        }
    }
}

/// Run the loop.  `all` collects every hit; otherwise it stops after the pair
/// that produced the first one (the pair's four leaf tests are all charged).
pub fn run<L: Fn(u64) -> bool>(
    f: &Onb,
    b: &Base,
    r: (u64, u64),
    firsts: &[usize],
    second: Second,
    leaf: &L,
    all: bool,
) -> LoopResult {
    let c0 = f.snapshot();
    let t0 = Instant::now();
    let m = b.xs.len();
    let rp = Pt::Aff(r.0, r.1);
    let nr = neg(rp);
    let mut den: Vec<u64> = Vec::with_capacity(2 * m);
    let mut scratch: Vec<u64> = Vec::with_capacity(2 * m);
    let mut res = LoopResult {
        pairs: 0,
        thirds: 0,
        exact_pairs: 0,
        rows: 0,
        hits: Vec::new(),
        found: false,
        counts: Counts::default(),
        nanos: 0,
    };
    'rows: for &i1 in firsts {
        res.rows += 1;
        let p1 = Pt::Aff(b.xs[i1], b.ys[i1]);
        let start = match second {
            Second::All => 0,
            Second::FromFirst => i1,
        };
        let cnt = m - start;
        let (pu, pv) = (add(f, p1, rp), add(f, p1, nr));
        let (Pt::Aff(xu, yu), Pt::Aff(xv, yv)) = (pu, pv) else {
            // P₁ = ±R: only when x_R is in the base.  Exact for the row.
            for i2 in start..m {
                res.pairs += 1;
                exact_pair(f, b, i1, i2, r, leaf, &mut res);
                if res.found && !all {
                    break 'rows;
                }
            }
            continue;
        };
        den.clear();
        den.resize(2 * cnt, 0);
        for k in 0..cnt {
            let x2 = b.xs[start + k];
            den[k] = f.add_c(xu, x2);
            den[cnt + k] = f.add_c(xv, x2);
        }
        f.batch_inv(&mut den[..2 * cnt], &mut scratch);
        for k in 0..cnt {
            let i2 = start + k;
            let (x2, y2) = (b.xs[i2], b.ys[i2]);
            res.pairs += 1;
            let (j, jj) = (den[k], den[cnt + k]);
            if j == 0 || jj == 0 {
                exact_pair(f, b, i1, i2, r, leaf, &mut res);
                if res.found && !all {
                    break 'rows;
                }
                continue;
            }
            let d = f.add_c(xu, x2);
            let lam = f.mul_c(f.add_c(yu, y2), j);
            let t1 = f.add_c(f.add_c(f.sqr_c(lam), lam), d);
            let lam2 = f.add_c(lam, f.mul_c(x2, j));
            let t2 = f.add_c(f.add_c(f.sqr_c(lam2), lam2), d);
            let dd = f.add_c(xv, x2);
            let mu = f.mul_c(f.add_c(yv, y2), jj);
            let t3 = f.add_c(f.add_c(f.sqr_c(mu), mu), dd);
            let mu2 = f.add_c(mu, f.mul_c(x2, jj));
            let t4 = f.add_c(f.add_c(f.sqr_c(mu2), mu2), dd);
            res.thirds += 4;
            let h = [leaf(t1), leaf(t2), leaf(t3), leaf(t4)];
            if h[0] | h[1] | h[2] | h[3] {
                res.found = true;
                for (t, &hit) in [t1, t2, t3, t4].iter().zip(&h) {
                    if hit {
                        res.hits.push(sorted3(b.xs[i1], x2, *t));
                    }
                }
                if !all {
                    break 'rows;
                }
            }
        }
    }
    res.counts = f.snapshot().sub(&c0);
    res.nanos = t0.elapsed().as_nanos();
    res
}

pub fn dedup(mut v: Vec<[u64; 3]>) -> Vec<[u64; 3]> {
    v.sort_unstable();
    v.dedup();
    v
}
