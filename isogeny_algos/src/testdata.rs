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

// ------------------------------------------------------------------ Kani instances

/// A point of exact order m (m | p + 1) on a supersingular curve over F_{p^2} whose group is
/// (Z/(p+1))^2.
pub fn ss_point_of_order(f: &Zp2, e: &Curve<(u64, u64)>, p: u64, m: u64, rng: &mut Rng) -> Pt<(u64, u64)> {
    let fac = factor_u64(m);
    loop {
        let r = random_point_f(f, e, rng);
        let q = pmul(f, e, &r, ((p + 1) / m) as u128);
        if fac.iter().all(|&(l, _)| pmul(f, e, &q, (m / l) as u128) != Pt::Inf) {
            return q;
        }
    }
}

/// Basis of E[2^m] on such a curve.
pub fn ss_basis2(f: &Zp2, e: &Curve<(u64, u64)>, p: u64, m: u32, rng: &mut Rng) -> (Pt<(u64, u64)>, Pt<(u64, u64)>) {
    let n = 1u64 << m;
    let a = ss_point_of_order(f, e, p, n, rng);
    let a2 = pmul(f, e, &a, (n / 2) as u128);
    loop {
        let b = ss_point_of_order(f, e, p, n, rng);
        if pmul(f, e, &b, (n / 2) as u128) != a2 {
            return (a, b);
        }
    }
}

/// Cyclic isogeny with kernel <r> of order m as prime-degree Velu steps: codomain and images.
pub fn cyclic_isogeny_chain<F: Field>(f: &F, e: &Curve<F::E>, r: &Pt<F::E>, m: u64, pts: &[Pt<F::E>]) -> (Curve<F::E>, Vec<Pt<F::E>>) {
    use crate::kernel::velu::velu_cyclic;
    let mut cur = *e;
    let mut ker = *r;
    let mut left = m;
    let mut pts = pts.to_vec();
    for (l, k) in factor_u64(m) {
        for _ in 0..k {
            let q = pmul(f, &cur, &ker, (left / l) as u128);
            let iso = velu_cyclic(f, &cur, &q, l);
            ker = iso.eval(f, &ker);
            for x in pts.iter_mut() {
                *x = iso.eval(f, x);
            }
            cur = iso.cod;
            left /= l;
        }
    }
    (cur, pts)
}

/// Kani diamond over F_{p^2}, p = 2^(a+2) 3^b d2 c - 1: phi: E0 -> E of degree 3^b,
/// gamma: E0 -> C of degree d2 = 2^a - 3^b, X = E0 / (ker phi + ker gamma); kernel generators
/// (gamma P, phi P), (gamma Q, phi Q) of order 2^(a+2) on C x E.
pub struct KaniInstance {
    pub p: u64,
    pub f: Zp2,
    pub a: u32,
    pub e0: Curve<(u64, u64)>,
    pub c: Curve<(u64, u64)>,
    pub e: Curve<(u64, u64)>,
    pub x: Curve<(u64, u64)>,
    pub k: [(Pt<(u64, u64)>, Pt<(u64, u64)>); 2],
    /// the second generator twisted by M = [[1, 2], [0, 1]]: isotropic, not from a diamond
    pub k_twisted: [(Pt<(u64, u64)>, Pt<(u64, u64)>); 2],
}

pub fn kani_instance(a: u32, b: u32, seed: u64) -> KaniInstance {
    let d1 = 3u64.pow(b);
    let d2 = (1u64 << a) - d1;
    let base = (1u64 << (a + 2)) * d1 * d2;
    let p = (1..).map(|c| base * c - 1).find(|&p| is_prime(p)).unwrap();
    let f = Zp2::new(p);
    let e0 = Curve::new(f.one(), f.zero());
    let mut rng = Rng::new(seed);
    let (pp, qq) = ss_basis2(&f, &e0, p, a + 2, &mut rng);
    let kphi = ss_point_of_order(&f, &e0, p, d1, &mut rng);
    let kgam = ss_point_of_order(&f, &e0, p, d2, &mut rng);
    let (e, im_phi) = cyclic_isogeny_chain(&f, &e0, &kphi, d1, &[pp, qq, kgam]);
    let (c, im_gam) = cyclic_isogeny_chain(&f, &e0, &kgam, d2, &[pp, qq]);
    let (x, _) = cyclic_isogeny_chain(&f, &e, &im_phi[2], d2, &[]);
    let q2 = padd(&f, &e, &im_phi[1], &pmul(&f, &e, &im_phi[0], 2));
    KaniInstance {
        p,
        f,
        a,
        e0,
        c,
        e,
        x,
        k: [(im_gam[0], im_phi[0]), (im_gam[1], im_phi[1])],
        k_twisted: [(im_gam[0], im_phi[0]), (im_gam[1], q2)],
    }
}
