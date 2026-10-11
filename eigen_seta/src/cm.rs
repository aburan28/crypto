//! Variable-time helpers for public data and key generation: point sampling,
//! Jacobian-coordinate scalar multiplication with public scalars, torsion
//! bases, the CM(-7) start curve E_D and its degree-2 endomorphism omega.

use crate::field::Field;
use crate::weier::{Isog2, WCurve, WIso, WPoint};
use num_bigint::BigUint;
use rand_core::RngCore;

/// Jacobian point (X : Y : Z), x = X/Z^2, y = Y/Z^3.  Variable time.
#[derive(Clone, Copy, Debug)]
struct Jac<F: Field> {
    x: F,
    y: F,
    z: F,
}

fn jac_dbl<F: Field>(e: &WCurve<F>, p: &Jac<F>) -> Jac<F> {
    let z2 = p.z.square();
    let y2 = p.y.square();
    let m = p.x.square().mul_small(3) + (e.a2 * p.x * z2).double() + e.a4 * z2.square();
    let z3 = (p.y * p.z).double();
    let x3 = m.square() - e.a2 * z3.square() - (p.x * y2).mul_small(8);
    let y3 = m * ((p.x * y2).mul_small(4) - x3) - y2.square().mul_small(8);
    Jac { x: x3, y: y3, z: z3 }
}

fn jac_add<F: Field>(e: &WCurve<F>, p: &Jac<F>, q: &Jac<F>) -> Jac<F> {
    if bool::from(p.z.is_zero()) {
        return *q;
    }
    if bool::from(q.z.is_zero()) {
        return *p;
    }
    let z1s = p.z.square();
    let z2s = q.z.square();
    let u1 = p.x * z2s;
    let u2 = q.x * z1s;
    let s1 = p.y * z2s * q.z;
    let s2 = q.y * z1s * p.z;
    let h = u2 - u1;
    let r = s2 - s1;
    if bool::from(h.is_zero()) {
        if bool::from(r.is_zero()) {
            return jac_dbl(e, p);
        }
        return Jac { x: F::one(), y: F::one(), z: F::zero() };
    }
    let z3 = h * p.z * q.z;
    let h2 = h.square();
    let x3 = r.square() - e.a2 * z3.square() - (u1 + u2) * h2;
    let y3 = r * (u1 * h2 - x3) - s1 * h2 * h;
    Jac { x: x3, y: y3, z: z3 }
}

fn to_jac<F: Field>(p: &WPoint<F>) -> Jac<F> {
    if bool::from(p.inf) {
        Jac { x: F::one(), y: F::one(), z: F::zero() }
    } else {
        Jac { x: p.x, y: p.y, z: F::one() }
    }
}

fn from_jac<F: Field>(p: &Jac<F>) -> WPoint<F> {
    if bool::from(p.z.is_zero()) {
        return WPoint::infinity();
    }
    let zi = p.z.invert();
    let zi2 = zi.square();
    WPoint::new(p.x * zi2, p.y * zi2 * zi)
}

/// [k]P for a public scalar.  Variable time.
pub fn mul_vartime<F: Field>(e: &WCurve<F>, p: &WPoint<F>, k: &BigUint) -> WPoint<F> {
    let pj = to_jac(p);
    let mut r = Jac { x: F::one(), y: F::one(), z: F::zero() };
    for i in (0..k.bits()).rev() {
        r = jac_dbl(e, &r);
        if k.bit(i) {
            r = jac_add(e, &r, &pj);
        }
    }
    from_jac(&r)
}

pub fn add_vartime<F: Field>(e: &WCurve<F>, p: &WPoint<F>, q: &WPoint<F>) -> WPoint<F> {
    from_jac(&jac_add(e, &to_jac(p), &to_jac(q)))
}

pub fn point_from_x<F: Field>(e: &WCurve<F>, x: F) -> Option<WPoint<F>> {
    let (y, ok) = e.rhs(&x).sqrt();
    if bool::from(ok) {
        Some(WPoint::new(x, y))
    } else {
        None
    }
}

pub fn random_point<F: Field, R: RngCore>(e: &WCurve<F>, rng: &mut R) -> WPoint<F> {
    loop {
        let x = F::from_counter(rng.next_u64()) * F::from_counter(rng.next_u64());
        if let Some(p) = point_from_x(e, x) {
            return p;
        }
    }
}

pub fn p_plus_1<F: Field>() -> BigUint {
    let l = F::char_limbs();
    let mut b = BigUint::default();
    for x in l.iter().rev() {
        b = (b << 64u32) + BigUint::from(*x);
    }
    b + 1u32
}

/// Point of exact order l^k (l prime, l^k | p + 1) from random sampling.
pub fn point_of_order<F: Field, R: RngCore>(e: &WCurve<F>, l: u32, k: usize, rng: &mut R) -> WPoint<F> {
    let lk = BigUint::from(l).pow(k as u32);
    let cof = p_plus_1::<F>() / &lk;
    let lk1 = BigUint::from(l).pow(k as u32 - 1);
    loop {
        let r = mul_vartime(e, &random_point(e, rng), &cof);
        if !bool::from(mul_vartime(e, &r, &lk1).inf) {
            return r;
        }
    }
}

/// Basis (P, Q) of E[l^k].
pub fn torsion_basis<F: Field, R: RngCore>(e: &WCurve<F>, l: u32, k: usize, rng: &mut R) -> (WPoint<F>, WPoint<F>) {
    let lk1 = BigUint::from(l).pow(k as u32 - 1);
    let p = point_of_order(e, l, k, rng);
    let pl = mul_vartime(e, &p, &lk1);
    loop {
        let q = point_of_order(e, l, k, rng);
        let ql = mul_vartime(e, &q, &lk1);
        // independent iff [l^(k-1)]Q is not a multiple of [l^(k-1)]P
        let mut dep = false;
        let mut m = WPoint::infinity();
        for _ in 0..l {
            m = add_vartime(e, &m, &pl);
            if bool::from(m.ct_eq(&ql)) {
                dep = true;
            }
        }
        if !dep {
            return (p, q);
        }
    }
}

/// The degree-2 endomorphism omega of E_D (omega^2 - t omega + 2 = 0).
#[derive(Clone, Copy, Debug)]
pub struct Omega<F: Field> {
    phi: Isog2<F>,
    iso: WIso<F>,
    pub trace: i64,
}

impl<F: Field> Omega<F> {
    pub fn eval(&self, p: &WPoint<F>) -> WPoint<F> {
        self.iso.eval(&self.phi.eval(p))
    }

    pub fn find<R: RngCore>(ed: &WCurve<F>, rng: &mut R) -> Self {
        let half = p_plus_1::<F>() >> 1;
        let mut xs: Vec<F> = vec![];
        while xs.len() < 3 {
            let t = mul_vartime(ed, &random_point(ed, rng), &half);
            if !bool::from(t.inf) && !xs.iter().any(|x| *x == t.x) {
                xs.push(t.x);
            }
        }
        let jd = ed.j_invariant();
        for x0 in xs {
            let phi = Isog2::new(ed, &x0);
            if phi.codomain.j_invariant() != jd {
                continue;
            }
            let base = WIso::find(&phi.codomain, ed);
            for iso in [base, base.negate()] {
                let om = Omega { phi, iso, trace: 0 };
                let r = random_point(ed, rng);
                let w = om.eval(&r);
                let ww = om.eval(&w);
                let two_r = ed.add(&r, &r);
                for t in [1i64, -1] {
                    // omega^2 - t omega + 2 = 0
                    let tw = if t == 1 { w } else { w.neg() };
                    let lhs = ed.add(&ed.sub(&ww, &tw), &two_r);
                    if bool::from(lhs.inf) {
                        return Omega { trace: t, ..om };
                    }
                }
            }
        }
        panic!("no degree-2 endomorphism found (is p inert in Q(sqrt(-7))?)");
    }
}

/// [k]P for a public signed small integer.
pub fn mul_i64<F: Field>(e: &WCurve<F>, p: &WPoint<F>, k: i64) -> WPoint<F> {
    let r = mul_vartime(e, p, &BigUint::from(k.unsigned_abs()));
    if k < 0 {
        r.neg()
    } else {
        r
    }
}

/// Eigenvector of an endomorphism `w` (given as a closure) for eigenvalue
/// `lam` with `other` the second eigenvalue: P = (w - other)(R), retried until
/// P has exact order 2^k.  Variable time (public / keygen data).
pub fn eigenvector_2k<F: Field, R: RngCore>(
    e: &WCurve<F>,
    w: &dyn Fn(&WPoint<F>) -> WPoint<F>,
    other: &BigUint,
    k: usize,
    rng: &mut R,
) -> WPoint<F> {
    let half = BigUint::from(1u32) << (k - 1);
    loop {
        let r = point_of_order(e, 2, k, rng);
        let p = add_vartime(e, &w(&r), &mul_vartime(e, &r, other).neg());
        if !bool::from(mul_vartime(e, &p, &half).inf) {
            return p;
        }
    }
}

/// Odd-degree (l = 3) Vélu isogeny on full points, variable time
/// (tests / key generation only).  Kernel {O, Q, -Q} with Q of order 3.
#[derive(Clone, Copy, Debug)]
pub struct Isog3W<F: Field> {
    xq: F,
    vq: F,
    uq: F,
    pub codomain: WCurve<F>,
}

impl<F: Field> Isog3W<F> {
    pub fn new(e: &WCurve<F>, q: &WPoint<F>) -> Self {
        let gx = (q.x.mul_small(3) + e.a2.double()) * q.x + e.a4;
        let vq = gx.double();
        let uq = q.y.square().mul_small(4);
        let w = uq + q.x * vq;
        let a4 = e.a4 - vq.mul_small(5);
        let a6 = e.a6 - e.a2.mul_small(4) * vq - w.mul_small(7);
        Isog3W { xq: q.x, vq, uq, codomain: WCurve::new(e.a2, a4, a6) }
    }
    pub fn eval(&self, p: &WPoint<F>) -> WPoint<F> {
        if bool::from(p.inf) {
            return *p;
        }
        let dx = p.x - self.xq;
        if bool::from(dx.is_zero()) {
            return WPoint::infinity();
        }
        let i = dx.invert();
        let i2 = i.square();
        let x = p.x + self.vq * i + self.uq * i2;
        let dxp = F::one() - self.vq * i2 - (self.uq * i2 * i).double();
        WPoint::new(x, p.y * dxp)
    }
}

/// Chain of 3-isogenies with kernel <K> (K of order 3^k), pushing `pts`.
pub fn chain3_vartime<F: Field>(e: &WCurve<F>, k: &WPoint<F>, n: usize, pts: &mut [WPoint<F>]) -> WCurve<F> {
    let mut cur = *e;
    let mut ker = *k;
    for i in 0..n {
        let q = mul_vartime(&cur, &ker, &BigUint::from(3u32).pow((n - 1 - i) as u32));
        let phi = Isog3W::new(&cur, &q);
        ker = phi.eval(&ker);
        for p in pts.iter_mut() {
            *p = phi.eval(p);
        }
        cur = phi.codomain;
    }
    cur
}
