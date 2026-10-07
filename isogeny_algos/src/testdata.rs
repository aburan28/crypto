//! Deterministic workload generators shared by tests and benchmarks.
use crate::curve::*;
use crate::field::*;

/// Random curve over F_p with a rational point of exact prime order `ell`.
pub fn curve_with_point(fp: &Zp, ell: u64, rng: &mut Rng) -> (Curve<u64>, Pt<u64>) {
    loop {
        let c = Curve::new(fp.random(rng), fp.random(rng));
        if !is_smooth(fp, &c) {
            continue;
        }
        let n = order(fp, &c, rng);
        if n % ell != 0 {
            continue;
        }
        let r = c.random_point(fp, rng);
        let q = pmul(fp, &c, &r, (n / ell) as u128);
        if q != Pt::Inf && pmul(fp, &c, &q, ell as u128) == Pt::Inf {
            return (c, q);
        }
    }
}

use crate::path::graph::*;

/// Non-backtracking random walk on the j-line over the given primes.
pub fn random_walk<F: Field>(f: &F, cache: &PhiCache, ells: &[usize], j: F::E, steps: usize, rng: &mut Rng) -> Path<F::E> {
    let mut path = Path { js: vec![j], ells: vec![] };
    let mut prev = None;
    let mut cur = j;
    for _ in 0..steps {
        let ns = neighbors(f, cache, ells, cur, rng);
        if ns.is_empty() {
            break;
        }
        let nb: Vec<_> = ns.iter().copied().filter(|&(_, n)| Some(n) != prev).collect();
        let ns = if nb.is_empty() { ns } else { nb };
        let (l, nx) = ns[rng.below(ns.len() as u64) as usize];
        prev = Some(cur);
        cur = nx;
        path.js.push(cur);
        path.ells.push(l);
    }
    path
}

/// Random ordinary curve with the Frobenius trace; `pred(t)` filters on the trace.
pub fn curve_with_trace(fp: &Zp, rng: &mut Rng, pred: impl Fn(i64) -> bool) -> (Curve<u64>, i64) {
    loop {
        let c = Curve::new(fp.random(rng), fp.random(rng));
        if !is_smooth(fp, &c) {
            continue;
        }
        let j = jinv(fp, &c);
        if j == 0 || j == 1728 {
            continue;
        }
        let t = fp.p as i64 + 1 - order(fp, &c, rng) as i64;
        if pred(t) {
            return (c, t);
        }
    }
}
