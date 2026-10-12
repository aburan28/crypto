//! x-only arithmetic on Montgomery curves `B y² = x³ + A x² + x`.
//!
//! Points are projective `(X : Z)` on the Kummer line; `(1 : 0)` is the
//! point at infinity.  The x-only formulas do not depend on `B`, so a point
//! of the quadratic twist is handled exactly like a point of the curve.
//!
//! [`xdbl`], [`xadd`] and [`ladder`] are generic over [`FieldArith`], so the
//! same code can run over an extension field; this module implements
//! [`FieldArith`] only for [`PrimeField`].  Costs (in `M`): `xDBL` 5 (with
//! `a24 = (A + 2)/4` given), `xADD` 6, ladder `11` per scalar bit.

use super::counted_fp::{OpCounter, PrimeField};
use super::IsogenyEvalError;

/// The ring operations the x-only formulas need.
pub trait FieldArith {
    /// A field element.
    type Elem: Clone + PartialEq + std::fmt::Debug;
    /// `a·b`.
    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;
    /// `a + b`.
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;
    /// `a − b`.
    fn sub(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;
    /// `0`.
    fn zero(&self) -> Self::Elem;
    /// `1`.
    fn one(&self) -> Self::Elem;
    /// Whether `a = 0`.
    fn is_zero(&self, a: &Self::Elem) -> bool;
}

impl<C: OpCounter> FieldArith for PrimeField<C> {
    type Elem = u64;
    #[inline]
    fn mul(&self, a: &u64, b: &u64) -> u64 {
        PrimeField::mul(self, *a, *b)
    }
    #[inline]
    fn add(&self, a: &u64, b: &u64) -> u64 {
        PrimeField::add(self, *a, *b)
    }
    #[inline]
    fn sub(&self, a: &u64, b: &u64) -> u64 {
        PrimeField::sub(self, *a, *b)
    }
    #[inline]
    fn zero(&self) -> u64 {
        0
    }
    #[inline]
    fn one(&self) -> u64 {
        1
    }
    #[inline]
    fn is_zero(&self, a: &u64) -> bool {
        *a == 0
    }
}

/// A projective x-only point `(X : Z)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct XZ<E> {
    /// `X`.
    pub x: E,
    /// `Z`; zero exactly for the point at infinity.
    pub z: E,
}

impl XZ<u64> {
    /// The affine point `(x : 1)`.
    pub fn affine(x: u64) -> Self {
        XZ { x, z: 1 }
    }

    /// The point at infinity `(1 : 0)`.
    pub fn infinity() -> Self {
        XZ { x: 1, z: 0 }
    }

    /// Whether this is the point at infinity.
    pub fn is_infinity(&self) -> bool {
        self.z == 0
    }

    /// Projective equality `X₁Z₂ = X₂Z₁` (two `M`).
    pub fn same_as<C: OpCounter>(&self, f: &PrimeField<C>, o: &XZ<u64>) -> bool {
        match (self.z == 0, o.z == 0) {
            (true, true) => true,
            (false, false) => f.mul(self.x, o.z) == f.mul(o.x, self.z),
            _ => false,
        }
    }

    /// The affine `x = X/Z` (one `I` and one `M`), `None` at infinity.
    pub fn to_affine<C: OpCounter>(&self, f: &PrimeField<C>) -> Option<u64> {
        Some(f.mul(self.x, f.inv(self.z)?))
    }
}

/// A Montgomery curve `B y² = x³ + A x² + x` over `F_p`, `A² ≠ 4`.  `B` is
/// not stored: x-only arithmetic does not need it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MontgomeryCurve {
    /// `A`.
    pub a: u64,
    /// `(A + 2)/4`.
    pub a24: u64,
}

impl MontgomeryCurve {
    /// The curve with coefficient `a` (reduced mod `p`).  Errors with
    /// [`IsogenyEvalError::SingularCurve`] if `A² = 4`.
    pub fn new<C: OpCounter>(f: &PrimeField<C>, a: u64) -> Result<Self, IsogenyEvalError> {
        let a = f.elem(a);
        if f.sqr(a) == 4 % f.p() {
            return Err(IsogenyEvalError::SingularCurve);
        }
        let inv4 = f.inv(4).ok_or(IsogenyEvalError::SingularCurve)?;
        let a24 = f.mul(f.add(a, 2), inv4);
        Ok(MontgomeryCurve { a, a24 })
    }

    /// `j = 256 (A² − 3)³ / (A² − 4)`.
    pub fn j_invariant<C: OpCounter>(&self, f: &PrimeField<C>) -> u64 {
        let a2 = f.sqr(self.a);
        let t = f.sub(a2, 3);
        let num = f.mul_small(256, f.mul(t, f.sqr(t)));
        // A² ≠ 4 by construction, so the inverse exists.
        f.div(num, f.sub(a2, 4)).unwrap_or(0)
    }

    /// The short Weierstrass model `y² = x³ + a x + b` of the `B = 1`
    /// curve, through `x_W = x + A/3`: `a = 1 − A²/3`,
    /// `b = (2A³ − 9A)/27`.  Returns `(a, b, A/3)`.
    pub fn short_weierstrass<C: OpCounter>(&self, f: &PrimeField<C>) -> (u64, u64, u64) {
        // p ≥ 5, so 3 and 27 are invertible.
        let inv3 = f.inv(3).unwrap_or(0);
        let shift = f.mul(self.a, inv3);
        let a2 = f.sqr(self.a);
        let aw = f.sub(1, f.mul(a2, inv3));
        let num = f.sub(f.mul_small(2, f.mul(a2, self.a)), f.mul_small(9, self.a));
        let bw = f.mul(num, f.inv(27 % f.p()).unwrap_or(0));
        (aw, bw, shift)
    }

    /// Whether `x³ + A x² + x` is a square, i.e. whether `x` comes from a
    /// point of the `B = 1` curve rather than of its twist.
    pub fn x_on_b1_curve<C: OpCounter>(&self, f: &PrimeField<C>, x: u64) -> bool {
        let rhs = f.mul(x, f.add(f.mul(x, f.add(x, self.a)), 1));
        f.legendre(rhs) >= 0
    }
}

/// `x([2]P)`: `X = (X+Z)²(X−Z)²`, `Z = 4XZ·((X−Z)² + a24·4XZ)`.
pub fn xdbl<F: FieldArith>(f: &F, p: &XZ<F::Elem>, a24: &F::Elem) -> XZ<F::Elem> {
    let t0 = f.add(&p.x, &p.z);
    let t1 = f.sub(&p.x, &p.z);
    let t0 = f.mul(&t0, &t0);
    let t1 = f.mul(&t1, &t1);
    let x = f.mul(&t0, &t1);
    let t = f.sub(&t0, &t1);
    let z = f.mul(&t, &f.add(&t1, &f.mul(a24, &t)));
    XZ { x, z }
}

/// `x(P + Q)` from `x(P)`, `x(Q)` and `x(P − Q)` (differential addition).
/// Requires `P − Q ≠ O`, i.e. `P ≠ Q` (use [`xdbl`] for doubling).
pub fn xadd<F: FieldArith>(
    f: &F,
    p: &XZ<F::Elem>,
    q: &XZ<F::Elem>,
    diff: &XZ<F::Elem>,
) -> XZ<F::Elem> {
    let u = f.mul(&f.sub(&p.x, &p.z), &f.add(&q.x, &q.z));
    let v = f.mul(&f.add(&p.x, &p.z), &f.sub(&q.x, &q.z));
    let s = f.add(&u, &v);
    let t = f.sub(&u, &v);
    let x = f.mul(&diff.z, &f.mul(&s, &s));
    let z = f.mul(&diff.x, &f.mul(&t, &t));
    XZ { x, z }
}

/// `x([n]P)` by the Montgomery ladder.  `n = 0` gives the point at
/// infinity; `P` may be projective.  If `P` is the point at infinity the
/// result is too.
pub fn ladder<F: FieldArith>(f: &F, p: &XZ<F::Elem>, n: u128, a24: &F::Elem) -> XZ<F::Elem> {
    if n == 0 || f.is_zero(&p.z) {
        return XZ {
            x: f.one(),
            z: f.zero(),
        };
    }
    let mut r0 = p.clone();
    let mut r1 = xdbl(f, p, a24);
    let bits = 128 - n.leading_zeros();
    for i in (0..bits - 1).rev() {
        if (n >> i) & 1 == 1 {
            r0 = xadd(f, &r1, &r0, p);
            r1 = xdbl(f, &r1, a24);
        } else {
            r1 = xadd(f, &r1, &r0, p);
            r0 = xdbl(f, &r0, a24);
        }
    }
    r0
}

/// Affine `x` of every point with one inversion (Montgomery's trick:
/// `3(n−1)` `M` + 1 `I` + `n` `M`).  `None` if any point is at infinity.
pub fn batch_to_affine<C: OpCounter>(f: &PrimeField<C>, pts: &[XZ<u64>]) -> Option<Vec<u64>> {
    let n = pts.len();
    if n == 0 {
        return Some(Vec::new());
    }
    if pts.iter().any(|q| q.z == 0) {
        return None;
    }
    let mut pre = vec![1u64; n];
    let mut acc = 1u64;
    for (i, q) in pts.iter().enumerate() {
        pre[i] = acc;
        acc = f.mul(acc, q.z);
    }
    let mut inv = f.inv(acc)?;
    let mut out = vec![0u64; n];
    for i in (0..n).rev() {
        out[i] = f.mul(pts[i].x, f.mul(inv, pre[i]));
        inv = f.mul(inv, pts[i].z);
    }
    Some(out)
}

/// Whether `x` is the x-coordinate of a point of exact order `ell` (an odd
/// prime) on the curve or its twist: `P ≠ O`, `[ell]P = O`.  Costs one
/// ladder (`≈ 11 log₂ ell` `M`).  Returns `false` if `ell` is not an odd
/// prime.
pub fn is_kernel_generator<C: OpCounter>(
    f: &PrimeField<C>,
    curve: &MontgomeryCurve,
    ell: u64,
    x: u64,
) -> bool {
    if ell < 3 || !super::counted_fp::is_prime_u64(ell) || x >= f.p() {
        return false;
    }
    ladder(f, &XZ::affine(x), ell as u128, &curve.a24).is_infinity()
}
