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

/// Kani diamond with an endomorphism as the auxiliary isogeny, usable at cryptographic size:
/// E0: y^2 = x^3 + x over F_{p^2} (p = 2^(a+2) 3^b c - 1), phi: E0 -> E of degree 3^b (Velu chain),
/// gamma = u + v i in End(E0) of degree d2 = u^2 + v^2 = 2^a - 3^b (i: (x, y) -> (-x, iota y)).
/// The (2^a, 2^a)-isogeny of E0 x E with kernel {(gamma P, phi P)} has codomain E0 x X.
pub struct KaniEndoInstance<F: Field> {
    pub e0: Curve<F::E>,
    pub e: Curve<F::E>,
    pub a: u32,
    pub k: [(Pt<F::E>, Pt<F::E>); 2],
    pub k_twisted: [(Pt<F::E>, Pt<F::E>); 2],
}

/// `p1` = p + 1; `iota` a square root of -1 in the field.
#[allow(clippy::too_many_arguments)]
pub fn kani_endomorphism_instance<F: Field>(
    f: &F,
    p1: &crate::bigint::Big,
    a: u32,
    b: u32,
    u: &crate::bigint::Big,
    v: &crate::bigint::Big,
    iota: F::E,
    rng: &mut Rng,
) -> KaniEndoInstance<F> {
    use crate::bigint::Big;
    use crate::kernel::chain::{ell_power_isogeny, Strategy};
    let e0 = Curve::new(f.one(), f.zero());
    let two_e = Big::from_u64(1).shl((a + 2) as usize);
    let three_b = (0..b).fold(Big::from_u64(1), |acc, _| acc.mul(&Big::from_u64(3)));
    let (cof2, _) = p1.divrem(&two_e);
    let (cof3, _) = p1.divrem(&three_b);
    let half = Big::from_u64(1).shl((a + 1) as usize);
    let third = (0..b - 1).fold(Big::from_u64(1), |acc, _| acc.mul(&Big::from_u64(3)));
    let point2 = |rng: &mut Rng| loop {
        let q = pmul_big(f, &e0, &random_point_f(f, &e0, rng), &cof2);
        if pmul_big(f, &e0, &q, &half) != Pt::Inf {
            return q;
        }
    };
    let pp = point2(rng);
    let pp_h = pmul_big(f, &e0, &pp, &half);
    let qq = loop {
        let q = point2(rng);
        if pmul_big(f, &e0, &q, &half) != pp_h {
            break q;
        }
    };
    let k3 = loop {
        let q = pmul_big(f, &e0, &random_point_f(f, &e0, rng), &cof3);
        if pmul_big(f, &e0, &q, &third) != Pt::Inf {
            break q;
        }
    };
    let (ch, pushed, _) = ell_power_isogeny(f, &e0, &k3, 3, b as usize, &[pp, qq], Strategy::Balanced);
    let e = ch.cod;
    let iota_map = |p: &Pt<F::E>| match *p {
        Pt::Inf => Pt::Inf,
        Pt::Aff(x, y) => Pt::Aff(f.neg(x), f.mul(iota, y)),
    };
    let gamma = |p: &Pt<F::E>| padd(f, &e0, &pmul_big(f, &e0, p, u), &pmul_big(f, &e0, &iota_map(p), v));
    let (gp, gq) = (gamma(&pp), gamma(&qq));
    let q2 = padd(f, &e, &pushed[1], &pmul(f, &e, &pushed[0], 2));
    KaniEndoInstance { e0, e, a, k: [(gp, pushed[0]), (gq, pushed[1])], k_twisted: [(gp, pushed[0]), (gq, q2)] }
}

/// c = a1^2 + a2^2 + a3^2 + a4^2 (randomised: the remainder after two random squares is split by
/// Cornacchia when it is a prime = 1 mod 4).
pub fn four_squares(c: &crate::int::Int, rng: &mut Rng) -> [crate::int::Int; 4] {
    use crate::int::Int;
    use crate::quat::klpt::cornacchia_sum_two_squares;
    let s = c.isqrt();
    let bits = s.mag().bits().max(1);
    loop {
        let rand_below = |rng: &mut Rng| {
            let mut v = Int::zero();
            for _ in 0..(bits + 63) / 64 {
                v = &(&v * &Int::from(1i64 << 32)) * &Int::from(1i64 << 32);
                v = &v + &Int::from((rng.next() >> 1) as i64);
            }
            v.modulo(&(&s + &Int::one()))
        };
        let a1 = rand_below(rng);
        let r1 = c - &(&a1 * &a1);
        if r1.is_neg() {
            continue;
        }
        let r1s = r1.isqrt();
        let a2 = if r1s.is_zero() { Int::zero() } else { rand_below(rng).modulo(&(&r1s + &Int::one())) };
        let r = &r1 - &(&a2 * &a2);
        if r.is_neg() {
            continue;
        }
        if (a1.is_zero() && a2.is_zero()) || r.is_zero() {
            continue;
        }
        if let Some((a3, a4)) = cornacchia_sum_two_squares(&r) {
            return [a1, a2, a3, a4];
        }
    }
}

/// [a] P + [b] i(P) on y^2 = x^3 + x (i: (x, y) -> (-x, iota y)), a, b signed.
pub fn gaussian_action<F: Field>(f: &F, e0: &Curve<F::E>, a: &crate::int::Int, b: &crate::int::Int, iota: F::E, p: &Pt<F::E>) -> Pt<F::E> {
    let smul = |k: &crate::int::Int, q: &Pt<F::E>| {
        let r = pmul_big(f, e0, q, k.mag());
        if k.is_neg() {
            neg(f, &r)
        } else {
            r
        }
    };
    let ip = match *p {
        Pt::Inf => Pt::Inf,
        Pt::Aff(x, y) => Pt::Aff(f.neg(x), f.mul(iota, y)),
    };
    padd(f, e0, &smul(a, p), &smul(b, &ip))
}

/// Dimension-4 Kani embedding (Robert 2022) of phi: E0 -> E of degree 3^b with
/// c = 2^n - 3^b = |u|^2 + |w|^2, u, w in Z[i] acting on E0: alpha = [[u, -conj(w)], [w, conj(u)]]
/// on E0^2, kernel {(alpha(P), Phi(P)) : P in E0^2[2^n]} in E0 x E0 x E x E (generators of order
/// 2^(n+2)). The (2^n, ..., 2^n)-isogeny splits (codomain E0^2 x a surface).
pub struct Kani4Instance<F: Field> {
    pub curves: Vec<Curve<F::E>>,
    pub k: Vec<Vec<Pt<F::E>>>,
    /// Phi on the second basis point twisted by M = [[1, 2], [0, 1]] (still isotropic)
    pub k_twisted: Vec<Vec<Pt<F::E>>>,
}

/// Generators of the dimension-4 Kani kernel from phi(B1), phi(B2) (B1, B2 a basis of
/// E0[2^(n+2)]) and c = |u|^2 + |w|^2 given as [a1, a2, a3, a4] (u = a1 + a2 i, w = a3 + a4 i).
#[allow(clippy::too_many_arguments)]
pub fn kani4_kernel<F: Field>(
    f: &F,
    e0: &Curve<F::E>,
    b1: &Pt<F::E>,
    b2: &Pt<F::E>,
    pb1: &Pt<F::E>,
    pb2: &Pt<F::E>,
    sq: &[crate::int::Int; 4],
    iota: F::E,
) -> Vec<Vec<Pt<F::E>>> {
    let [a1, a2, a3, a4] = sq;
    let (na2, na3) = (-a2, -a3);
    let u = |p: &Pt<F::E>| gaussian_action(f, e0, a1, a2, iota, p);
    let ub = |p: &Pt<F::E>| gaussian_action(f, e0, a1, &na2, iota, p);
    let w = |p: &Pt<F::E>| gaussian_action(f, e0, a3, a4, iota, p);
    let mwb = |p: &Pt<F::E>| gaussian_action(f, e0, &na3, a4, iota, p); // -conj(w) = -a3 + a4 i
    vec![
        vec![u(b1), w(b1), *pb1, Pt::Inf],
        vec![u(b2), w(b2), *pb2, Pt::Inf],
        vec![mwb(b1), ub(b1), Pt::Inf, *pb1],
        vec![mwb(b2), ub(b2), Pt::Inf, *pb2],
    ]
}

#[allow(clippy::too_many_arguments)]
pub fn kani4_instance<F: Field>(f: &F, p1: &crate::bigint::Big, n: u32, b: u32, iota: F::E, rng: &mut Rng) -> Kani4Instance<F> {
    use crate::bigint::Big;
    use crate::int::Int;
    use crate::kernel::chain::{ell_power_isogeny, Strategy};
    let e0 = Curve::new(f.one(), f.zero());
    let two_e = Big::from_u64(1).shl((n + 2) as usize);
    let three_b = (0..b).fold(Big::from_u64(1), |acc, _| acc.mul(&Big::from_u64(3)));
    let (cof2, _) = p1.divrem(&two_e);
    let (cof3, _) = p1.divrem(&three_b);
    let half = Big::from_u64(1).shl((n + 1) as usize);
    let third = (0..b.saturating_sub(1)).fold(Big::from_u64(1), |acc, _| acc.mul(&Big::from_u64(3)));
    let point2 = |rng: &mut Rng| loop {
        let q = pmul_big(f, &e0, &random_point_f(f, &e0, rng), &cof2);
        if pmul_big(f, &e0, &q, &half) != Pt::Inf {
            return q;
        }
    };
    let b1 = point2(rng);
    let h1 = pmul_big(f, &e0, &b1, &half);
    let b2 = loop {
        let q = point2(rng);
        if pmul_big(f, &e0, &q, &half) != h1 {
            break q;
        }
    };
    let k3 = loop {
        let q = pmul_big(f, &e0, &random_point_f(f, &e0, rng), &cof3);
        if pmul_big(f, &e0, &q, &third) != Pt::Inf {
            break q;
        }
    };
    let (ch, pushed, _) = ell_power_isogeny(f, &e0, &k3, 3, b as usize, &[b1, b2], Strategy::Balanced);
    let e = ch.cod;
    let c = &Int::from_big(&Big::from_u64(1).shl(n as usize)) - &Int::from_big(&three_b);
    let sq = four_squares(&c, rng);
    let k = kani4_kernel(f, &e0, &b1, &b2, &pushed[0], &pushed[1], &sq, iota);
    let tw = padd(f, &e, &pushed[1], &pmul(f, &e, &pushed[0], 2));
    let k_twisted = kani4_kernel(f, &e0, &b1, &b2, &pushed[0], &tw, &sq, iota);
    Kani4Instance { curves: vec![e0, e0, e, e], k, k_twisted }
}
