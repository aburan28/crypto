//! Full points on y^2 = x^3 + a2 x^2 + a4 x + a6, affine with an explicit
//! infinity flag.  Group law, 2-isogenies (Vélu) and isomorphisms are
//! computed branch-free: exceptional cases are computed alongside the
//! generic case and selected with `subtle`.

use crate::ct::{Choice, ConditionallySelectable};
use crate::field::Field;

#[derive(Clone, Copy, Debug)]
pub struct WCurve<F: Field> {
    pub a2: F,
    pub a4: F,
    pub a6: F,
}

#[derive(Clone, Copy, Debug)]
pub struct WPoint<F: Field> {
    pub x: F,
    pub y: F,
    pub inf: Choice,
}

impl<F: Field> ConditionallySelectable for WPoint<F> {
    fn conditional_select(a: &Self, b: &Self, c: Choice) -> Self {
        WPoint {
            x: F::conditional_select(&a.x, &b.x, c),
            y: F::conditional_select(&a.y, &b.y, c),
            inf: Choice::conditional_select(&a.inf, &b.inf, c),
        }
    }
}

impl<F: Field> WPoint<F> {
    pub fn new(x: F, y: F) -> Self {
        WPoint { x, y, inf: Choice::from(0) }
    }
    pub fn infinity() -> Self {
        WPoint { x: F::zero(), y: F::one(), inf: Choice::from(1) }
    }
    pub fn neg(&self) -> Self {
        WPoint { x: self.x, y: -self.y, inf: self.inf }
    }
    /// Constant-time equality (both infinity, or equal affine coordinates).
    pub fn ct_eq(&self, o: &Self) -> Choice {
        let both_inf = self.inf & o.inf;
        let both_fin = !self.inf & !o.inf;
        both_inf | (both_fin & self.x.ct_eq(&o.x) & self.y.ct_eq(&o.y))
    }
}

impl<F: Field> WCurve<F> {
    pub fn new(a2: F, a4: F, a6: F) -> Self {
        WCurve { a2, a4, a6 }
    }
    /// Montgomery y^2 = x^3 + A x^2 + x.
    pub fn from_montgomery(a: F) -> Self {
        WCurve { a2: a, a4: F::one(), a6: F::zero() }
    }
    pub fn rhs(&self, x: &F) -> F {
        ((*x + self.a2) * *x + self.a4) * *x + self.a6
    }
    pub fn on_curve(&self, p: &WPoint<F>) -> Choice {
        p.inf | p.y.square().ct_eq(&self.rhs(&p.x))
    }
    /// Short Weierstrass coefficients (A, B) of the curve shifted by a2/3.
    pub fn short_ab(&self) -> (F, F) {
        let inv3 = F::from_i64(3).invert();
        let inv27 = F::from_i64(27).invert();
        let a2sq = self.a2.square();
        let a = self.a4 - a2sq * inv3;
        let b = self.a6 - self.a2 * self.a4 * inv3 + (a2sq * self.a2).double() * inv27;
        (a, b)
    }
    pub fn j_invariant(&self) -> F {
        let (a, b) = self.short_ab();
        let a3 = a.square() * a * F::from_i64(4);
        F::from_i64(1728) * a3 * (a3 + F::from_i64(27) * b.square()).invert()
    }

    /// P + Q, all cases, constant time (one inversion).
    pub fn add(&self, p: &WPoint<F>, q: &WPoint<F>) -> WPoint<F> {
        let same_x = p.x.ct_eq(&q.x);
        // one inversion: denominator x2 - x1 (addition) or 2 y1 (doubling)
        let den = F::conditional_select(&(q.x - p.x), &p.y.double(), same_x);
        let num_add = q.y - p.y;
        let num_dbl = (p.x.mul_small(3) + self.a2.double()) * p.x + self.a4;
        let num = F::conditional_select(&num_add, &num_dbl, same_x);
        let lam = num * den.invert();
        let x3 = lam.square() - self.a2 - p.x - q.x;
        let y3 = lam * (p.x - x3) - p.y;
        let mut r = WPoint::new(x3, y3);
        // P = -Q (including doubling a 2-torsion point) -> infinity
        let opposite = same_x & (p.y + q.y).is_zero();
        r = WPoint::conditional_select(&r, &WPoint::infinity(), opposite);
        r = WPoint::conditional_select(&r, q, p.inf);
        r = WPoint::conditional_select(&r, p, q.inf & !p.inf);
        r
    }

    pub fn sub(&self, p: &WPoint<F>, q: &WPoint<F>) -> WPoint<F> {
        self.add(p, &q.neg())
    }

    /// [2^n]P via Jacobian coordinates (no inversion until the end).
    pub fn dbl_iter(&self, p: &WPoint<F>, n: usize) -> WPoint<F> {
        let (mut x, mut y, mut z) = (p.x, p.y, F::conditional_select(&F::one(), &F::zero(), p.inf));
        for _ in 0..n {
            let z2 = z.square();
            let y2 = y.square();
            let m = x.square().mul_small(3) + (self.a2 * x * z2).double() + self.a4 * z2.square();
            let z3 = (y * z).double();
            let x3 = m.square() - self.a2 * z3.square() - (x * y2).mul_small(8);
            let y3 = m * ((x * y2).mul_small(4) - x3) - y2.square().mul_small(8);
            x = x3;
            y = y3;
            z = z3;
        }
        let zi = z.invert();
        let zi2 = zi.square();
        let inf = z.is_zero();
        let r = WPoint { x: x * zi2, y: y * zi2 * zi, inf: Choice::from(0) };
        WPoint::conditional_select(&r, &WPoint::infinity(), inf)
    }

    /// [k]P for a small public k (keygen / tests only).
    pub fn mul_small_vartime(&self, p: &WPoint<F>, k: u64) -> WPoint<F> {
        let mut r = WPoint::infinity();
        let mut b = *p;
        let mut k = k;
        while k > 0 {
            if k & 1 == 1 {
                r = self.add(&r, &b);
            }
            b = self.add(&b, &b);
            k >>= 1;
        }
        r
    }
}

/// 2-isogeny with kernel {O, (x0, 0)} (Vélu), constant time.
#[derive(Clone, Copy, Debug)]
pub struct Isog2<F: Field> {
    pub x0: F,
    pub t: F,
    pub codomain: WCurve<F>,
}

impl<F: Field> Isog2<F> {
    pub fn new(e: &WCurve<F>, x0: &F) -> Self {
        let t = (x0.mul_small(3) + e.a2.double()) * *x0 + e.a4;
        let a4 = e.a4 - t.mul_small(5);
        let a6 = e.a6 - (e.a2.mul_small(4) * t) - (*x0 * t).mul_small(7);
        Isog2 { x0: *x0, t, codomain: WCurve { a2: e.a2, a4, a6 } }
    }
    pub fn eval(&self, p: &WPoint<F>) -> WPoint<F> {
        let dx = p.x - self.x0;
        let in_kernel = dx.is_zero() & !p.inf;
        let i = dx.invert();
        let ti = self.t * i;
        let x = p.x + ti;
        let y = p.y * (F::one() - ti * i);
        let r = WPoint { x, y, inf: p.inf };
        WPoint::conditional_select(&r, &WPoint::infinity(), in_kernel)
    }
}

/// Isomorphism (x, y) -> (u^2 (x + s1) - s2, u^3 y) between two curves of
/// this shape, s1 = a2/3, s2 = a2'/3.
#[derive(Clone, Copy, Debug)]
pub struct WIso<F: Field> {
    pub u2: F,
    pub u3: F,
    pub s1: F,
    pub s2: F,
}

impl<F: Field> WIso<F> {
    /// Some isomorphism from `e1` to `e2` (the other one is its negative).
    /// Constant time; garbage if the curves are not isomorphic or have
    /// j = 0 / 1728.
    pub fn find(e1: &WCurve<F>, e2: &WCurve<F>) -> Self {
        let (a1, b1) = e1.short_ab();
        let (a2, b2) = e2.short_ab();
        // u^4 = a2/a1, u^6 = b2/b1  ->  u^2 = (b2 a1)/(b1 a2)
        let u2 = (b2 * a1) * (b1 * a2).invert();
        let (u, _) = u2.sqrt();
        let inv3 = F::from_i64(3).invert();
        WIso { u2, u3: u2 * u, s1: e1.a2 * inv3, s2: e2.a2 * inv3 }
    }
    pub fn eval(&self, p: &WPoint<F>) -> WPoint<F> {
        WPoint { x: self.u2 * (p.x + self.s1) - self.s2, y: self.u3 * p.y, inf: p.inf }
    }
    pub fn negate(&self) -> Self {
        WIso { u2: self.u2, u3: -self.u3, s1: self.s1, s2: self.s2 }
    }
}
