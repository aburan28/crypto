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
pub fn random_walk<F: Field>(
    f: &F,
    cache: &PhiCache<F>,
    ells: &[usize],
    j: F::E,
    steps: usize,
    rng: &mut Rng,
) -> Path<F::E> {
    let mut path = Path {
        js: vec![j],
        ells: vec![],
    };
    let mut prev = None;
    let mut cur = j;
    for _ in 0..steps {
        let ns = neighbors(f, cache, ells, cur, rng);
        if ns.is_empty() {
            break;
        }
        let nb: Vec<_> = ns
            .iter()
            .copied()
            .filter(|&(_, n)| Some(n) != prev)
            .collect();
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

/// A supersingular curve over a large prime field with rational points of several prime orders.
pub struct SsWorkload<const N: usize> {
    pub f: crate::fpm::FpM<N>,
    /// y^2 = x^3 + a x + b with #E(F_p) = p + 1 (j != 0, 1728)
    pub e: Curve<[u64; N]>,
    /// (l, point of exact order l) for each requested l
    pub pts: Vec<(u64, Pt<[u64; N]>)>,
}

/// p = 4 * prod(ells) * c - 1 prime with exactly `bits` bits (p = 3 mod 4), so y^2 = x^3 + x is
/// supersingular with p + 1 points and every l in `ells` divides #E(F_p). The curve is then moved
/// off j = 1728 by a few rational 3-isogenies (3 must be in `ells`).
pub fn supersingular_workload<const N: usize>(
    ells: &[u64],
    bits: usize,
    rng: &mut Rng,
) -> SsWorkload<N> {
    use crate::bigint::{is_probable_prime, Big};
    use crate::fpm::FpM;
    assert!(ells.contains(&3));
    let mut prod = Big::from_u64(4);
    for &l in ells {
        prod = prod.mul_small(l);
    }
    let cbits = bits - prod.bits();
    let p = loop {
        let mut limbs: Vec<u64> = (0..(cbits + 63) / 64).map(|_| rng.next()).collect();
        let top = (cbits - 1) % 64;
        let last = limbs.len() - 1;
        limbs[last] &= (1u64 << top << 1).wrapping_sub(1);
        limbs[last] |= 1u64 << top;
        let p = prod.mul(&Big::from_limbs(&limbs)).sub_small(1);
        if p.bits() == bits && is_probable_prime(&p) {
            break p;
        }
    };
    let f = FpM::<N>::new(&p);
    let p1 = p.add_small(1);
    let point_of_order = |e: &Curve<[u64; N]>, ell: u64, rng: &mut Rng| loop {
        let r = random_point_f(&f, e, rng);
        let q = pmul_big(&f, e, &r, &p1.divrem_small(ell).0);
        if q != Pt::Inf {
            assert!(pmul(&f, e, &q, ell as u128) == Pt::Inf);
            return q;
        }
    };
    let mut e = Curve::new(f.one(), f.zero());
    for _ in 0..6 {
        let k = point_of_order(&e, 3, rng);
        let Pt::Aff(x0, _) = k else { unreachable!() };
        e = crate::kernel::xonly::velu_xonly_fast(&f, &e, x0, 3).unwrap().cod;
    }
    let j = jinv(&f, &e);
    assert!(j != f.zero() && j != f.from_u64(1728));
    let pts = ells.iter().map(|&l| (l, point_of_order(&e, l, rng))).collect();
    SsWorkload { f, e, pts }
}
