//! x-only arithmetic on Montgomery curves y^2 = x^3 + A x^2 + x (B = 1),
//! projective (X : Z).  All routines are branch-free in the point values;
//! loops run over a public number of scalar bits.

use crate::ct::{choice_from_bit, Choice, ConditionallySelectable};
use crate::field::Field;
use crate::scalar::Scalar;

/// Projective x-coordinate (X : Z); Z = 0 is the point at infinity.
#[derive(Clone, Copy, Debug)]
pub struct XPoint<F: Field> {
    pub x: F,
    pub z: F,
}

impl<F: Field> ConditionallySelectable for XPoint<F> {
    fn conditional_select(a: &Self, b: &Self, c: Choice) -> Self {
        XPoint { x: F::conditional_select(&a.x, &b.x, c), z: F::conditional_select(&a.z, &b.z, c) }
    }
}

impl<F: Field> XPoint<F> {
    pub fn from_affine(x: F) -> Self {
        XPoint { x, z: F::one() }
    }
    pub fn infinity() -> Self {
        XPoint { x: F::one(), z: F::zero() }
    }
    pub fn is_infinity(&self) -> Choice {
        self.z.is_zero()
    }
    /// Affine x (0 for infinity).  One inversion.
    pub fn affine(&self) -> F {
        self.x * self.z.invert()
    }
    pub fn cswap(a: &mut Self, b: &mut Self, c: Choice) {
        let ta = *a;
        *a = Self::conditional_select(a, b, c);
        *b = Self::conditional_select(b, &ta, c);
    }
    /// Projective equality of x-coordinates.
    pub fn x_eq(&self, o: &Self) -> Choice {
        (self.x * o.z).ct_eq(&(o.x * self.z))
    }
}

/// Montgomery curve with projective constants: A24p = A + 2C, A24m = A - 2C, C24 = 4C.
#[derive(Clone, Copy, Debug)]
pub struct MontCurve<F: Field> {
    pub a: F, // affine A (C = 1)
}

#[derive(Clone, Copy, Debug)]
pub struct ProjCurve<F: Field> {
    pub a24p: F,
    pub a24m: F,
}

impl<F: Field> ProjCurve<F> {
    pub fn from_affine(a: F) -> Self {
        let two = F::from_i64(2);
        ProjCurve { a24p: a + two, a24m: a - two }
    }
    /// C24 = A24p - A24m (= 4C).
    pub fn c24(&self) -> F {
        self.a24p - self.a24m
    }
    /// Affine A = 2 (A24p + A24m) / (A24p - A24m).  One inversion.
    pub fn affine_a(&self) -> F {
        (self.a24p + self.a24m).double() * (self.a24p - self.a24m).invert()
    }
}

impl<F: Field> MontCurve<F> {
    pub fn new(a: F) -> Self {
        MontCurve { a }
    }
    pub fn proj(&self) -> ProjCurve<F> {
        ProjCurve::from_affine(self.a)
    }
    pub fn j_invariant(&self) -> F {
        let a2 = self.a.square();
        let t = a2 - F::from_i64(3);
        F::from_i64(256) * t.square() * t * (a2 - F::from_i64(4)).invert()
    }
    /// Right-hand side x^3 + A x^2 + x.
    pub fn rhs(&self, x: F) -> F {
        ((x + self.a) * x + F::one()) * x
    }
}

/// [2]P.
#[inline]
pub fn xdbl<F: Field>(p: &XPoint<F>, c: &ProjCurve<F>) -> XPoint<F> {
    let c24 = c.c24();
    let t0 = (p.x - p.z).square();
    let t1 = (p.x + p.z).square();
    let z2 = c24 * t0;
    let x2 = z2 * t1;
    let t1 = t1 - t0;
    let t0 = c.a24p * t1;
    let z2 = (z2 + t0) * t1;
    XPoint { x: x2, z: z2 }
}

/// x(P + Q) from x(P), x(Q), x(P - Q) (all projective).
#[inline]
pub fn xadd<F: Field>(p: &XPoint<F>, q: &XPoint<F>, d: &XPoint<F>) -> XPoint<F> {
    let a = p.x + p.z;
    let b = p.x - p.z;
    let cc = q.x + q.z;
    let dd = q.x - q.z;
    let da = dd * a;
    let cb = cc * b;
    XPoint { x: d.z * (da + cb).square(), z: d.x * (da - cb).square() }
}

/// ([2]P, P + Q) given x(P - Q).
#[inline]
pub fn xdbladd<F: Field>(p: &XPoint<F>, q: &XPoint<F>, d: &XPoint<F>, c: &ProjCurve<F>) -> (XPoint<F>, XPoint<F>) {
    (xdbl(p, c), xadd(p, q, d))
}

/// ([k]P, [k+1]P) for a scalar of `nbits` public bits.  Constant time.
pub fn ladder<F: Field>(p: &XPoint<F>, k: &Scalar, nbits: usize, c: &ProjCurve<F>) -> (XPoint<F>, XPoint<F>) {
    let mut r0 = XPoint::infinity();
    let mut r1 = *p;
    let mut prev = 0u64;
    for i in (0..nbits).rev() {
        let b = k.bit(i);
        XPoint::cswap(&mut r0, &mut r1, choice_from_bit(b ^ prev));
        prev = b;
        let (d, a) = xdbladd(&r0, &r1, p, c);
        r0 = d;
        r1 = a;
    }
    XPoint::cswap(&mut r0, &mut r1, choice_from_bit(prev));
    (r0, r1)
}

/// x(P + [m]Q) from x(P), x(Q), x(P - Q).  Constant time over `nbits`.
pub fn ladder3pt<F: Field>(
    xp: &XPoint<F>,
    xq: &XPoint<F>,
    xpq: &XPoint<F>,
    m: &Scalar,
    nbits: usize,
    c: &ProjCurve<F>,
) -> XPoint<F> {
    // A = [2^i]Q, B = P + [m mod 2^i]Q, D = B - A
    let mut a = *xq;
    let mut b = *xp;
    let mut d = *xpq;
    for i in 0..nbits {
        let bit = choice_from_bit(m.bit(i));
        // bit = 1: B <- B + A (diff D);  bit = 0: D <- D - A = xADD(D, A; diff B)
        XPoint::cswap(&mut b, &mut d, !bit);
        let nb = xadd(&b, &a, &d);
        b = nb;
        XPoint::cswap(&mut b, &mut d, !bit);
        a = xdbl(&a, c);
    }
    b
}

/// [3]P (Costello–Hisil / SIKE formulas).
#[inline]
pub fn xtpl<F: Field>(p: &XPoint<F>, c: &ProjCurve<F>) -> XPoint<F> {
    let t0 = p.x - p.z;
    let t2 = t0.square();
    let t1 = p.x + p.z;
    let t3 = t1.square();
    let t4 = t1 + t0;
    let t0 = t1 - t0;
    let t1 = t4.square() - t3 - t2;
    let t5 = t3 * c.a24p;
    let t3 = t3 * t5;
    let t6 = t2 * c.a24m;
    let t2 = t2 * t6;
    let t3 = t2 - t3;
    let t2 = t5 - t6;
    let t1 = t2 * t1;
    let t2 = (t3 + t1).square();
    let x3 = t2 * t4;
    let t1 = (t3 - t1).square();
    let z3 = t1 * t0;
    XPoint { x: x3, z: z3 }
}

pub fn xtpl_iter<F: Field>(p: &XPoint<F>, n: usize, c: &ProjCurve<F>) -> XPoint<F> {
    let mut r = *p;
    for _ in 0..n {
        r = xtpl(&r, c);
    }
    r
}

pub fn xdbl_iter<F: Field>(p: &XPoint<F>, n: usize, c: &ProjCurve<F>) -> XPoint<F> {
    let mut r = *p;
    for _ in 0..n {
        r = xdbl(&r, c);
    }
    r
}

/// 3-isogeny with kernel <P> (P of order 3): codomain constants and the
/// evaluation coefficients.
pub fn get_3_isog<F: Field>(p: &XPoint<F>) -> (ProjCurve<F>, [F; 2]) {
    let k0 = p.x - p.z;
    let t0 = k0.square();
    let k1 = p.x + p.z;
    let t1 = k1.square();
    let t2 = t0 + t1;
    let t3 = (k0 + k1).square() - t2;
    let t2 = t1 + t3;
    let t3 = t3 + t0;
    let t4 = (t3 + t0).double() + t1;
    let a24m = t2 * t4;
    let t4 = (t1 + t2).double() + t0;
    let a24p = t3 * t4;
    (ProjCurve { a24p, a24m }, [k0, k1])
}

pub fn eval_3_isog<F: Field>(q: &XPoint<F>, k: &[F; 2]) -> XPoint<F> {
    let t0 = k[0] * (q.x + q.z);
    let t1 = k[1] * (q.x - q.z);
    let t2 = (t0 + t1).square();
    let t0 = (t1 - t0).square();
    XPoint { x: q.x * t2, z: q.z * t0 }
}

/// Okeya–Sakurai y-recovery on y^2 = x^3 + A x^2 + x: given P = (xp, yp)
/// affine, Q = x([k]P), R = x([k+1]P), return projective (X : Y : Z) of [k]P.
pub fn recover_y<F: Field>(a: &F, xp: &F, yp: &F, q: &XPoint<F>, r: &XPoint<F>) -> (F, F, F) {
    let v1 = *xp * q.z;
    let v2 = q.x + v1;
    let v3 = (q.x - v1).square() * r.x;
    let v1 = a.double() * q.z;
    let v2 = v2 + v1;
    let v4 = *xp * q.x + q.z;
    let v2 = v2 * v4;
    let v1 = v1 * q.z;
    let v2 = (v2 - v1) * r.z;
    let y = v2 - v3;
    let v1 = yp.double() * q.z * r.z;
    let x = v1 * q.x;
    let z = v1 * q.z;
    (x, y, z)
}
