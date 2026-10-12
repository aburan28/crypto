//! F_{p^2} = F_p[i]/(i^2 + 1) for p ≡ 3 (mod 4).  Constant-time.

use crate::ct::{Choice, ConditionallySelectable, ConstantTimeEq};
use crate::fp::{fp_ops, Fp, FpConfig};
use core::fmt;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use zeroize::Zeroize;

pub struct Fp2<C: FpConfig<N>, const N: usize> {
    pub c0: Fp<C, N>,
    pub c1: Fp<C, N>,
}

impl<C: FpConfig<N>, const N: usize> Clone for Fp2<C, N> {
    #[inline(always)]
    fn clone(&self) -> Self {
        *self
    }
}
impl<C: FpConfig<N>, const N: usize> Copy for Fp2<C, N> {}

impl<C: FpConfig<N>, const N: usize> Zeroize for Fp2<C, N> {
    fn zeroize(&mut self) {
        self.c0.zeroize();
        self.c1.zeroize();
    }
}

impl<C: FpConfig<N>, const N: usize> fmt::Debug for Fp2<C, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:?} + i*{:?})", self.c0, self.c1)
    }
}

impl<C: FpConfig<N>, const N: usize> Fp2<C, N> {
    pub const ZERO: Self = Fp2 {
        c0: Fp::ZERO,
        c1: Fp::ZERO,
    };
    pub const ONE: Self = Fp2 {
        c0: Fp::ONE,
        c1: Fp::ZERO,
    };

    #[inline(always)]
    pub fn new(c0: Fp<C, N>, c1: Fp<C, N>) -> Self {
        Fp2 { c0, c1 }
    }
    #[inline(always)]
    pub fn zero() -> Self {
        Self::ZERO
    }
    #[inline(always)]
    pub fn one() -> Self {
        Self::ONE
    }
    pub fn i() -> Self {
        Fp2 {
            c0: Fp::ZERO,
            c1: Fp::ONE,
        }
    }
    pub fn from_fp(a: Fp<C, N>) -> Self {
        Fp2 {
            c0: a,
            c1: Fp::ZERO,
        }
    }
    pub fn from_i64(v: i64) -> Self {
        Self::from_fp(Fp::from_i64(v))
    }
    pub fn from_u64(v: u64) -> Self {
        Self::from_fp(Fp::from_u64(v))
    }

    #[inline]
    pub fn add(self, o: &Self) -> Self {
        Fp2 {
            c0: self.c0.add(&o.c0),
            c1: self.c1.add(&o.c1),
        }
    }
    #[inline]
    pub fn sub(self, o: &Self) -> Self {
        Fp2 {
            c0: self.c0.sub(&o.c0),
            c1: self.c1.sub(&o.c1),
        }
    }
    #[inline]
    pub fn neg(self) -> Self {
        Fp2 {
            c0: self.c0.neg(),
            c1: self.c1.neg(),
        }
    }
    #[inline]
    pub fn double(&self) -> Self {
        self.add(self)
    }
    /// Karatsuba: 3 base multiplications.
    #[inline]
    pub fn mul(self, o: &Self) -> Self {
        let v0 = self.c0.mul(&o.c0);
        let v1 = self.c1.mul(&o.c1);
        let s = self.c0.add(&self.c1).mul(&o.c0.add(&o.c1));
        Fp2 {
            c0: v0.sub(&v1),
            c1: s.sub(&v0).sub(&v1),
        }
    }
    /// (a+bi)^2 = (a+b)(a-b) + 2ab i: 2 base multiplications.
    #[inline]
    pub fn square(&self) -> Self {
        let t = self.c0.add(&self.c1).mul(&self.c0.sub(&self.c1));
        let u = self.c0.mul(&self.c1);
        Fp2 {
            c0: t,
            c1: u.double(),
        }
    }
    pub fn mul_fp(&self, s: &Fp<C, N>) -> Self {
        Fp2 {
            c0: self.c0.mul(s),
            c1: self.c1.mul(s),
        }
    }
    pub fn mul_small(&self, k: i64) -> Self {
        self.mul_fp(&Fp::from_i64(k))
    }
    pub fn conj(&self) -> Self {
        Fp2 {
            c0: self.c0,
            c1: self.c1.neg(),
        }
    }
    /// Frobenius x -> x^p (= conjugation since p ≡ 3 mod 4).
    pub fn frob(&self) -> Self {
        self.conj()
    }
    pub fn norm(&self) -> Fp<C, N> {
        self.c0.square().add(&self.c1.square())
    }
    /// 1/x; maps 0 to 0.  Constant time (one F_p Fermat inversion).
    pub fn invert(&self) -> Self {
        let n = self.norm().invert();
        Fp2 {
            c0: self.c0.mul(&n),
            c1: self.c1.neg().mul(&n),
        }
    }
    pub fn is_zero(&self) -> Choice {
        self.c0.is_zero() & self.c1.is_zero()
    }
    /// Squares of F_{p^2}: exactly the elements whose norm is a square in F_p.
    pub fn is_square(&self) -> Choice {
        self.norm().is_square()
    }

    /// x^e for a public exponent (little-endian limbs).
    pub fn pow_public(&self, e: &[u64]) -> Self {
        let mut r = Self::ONE;
        let mut started = false;
        for i in (0..e.len()).rev() {
            for bit in (0..64).rev() {
                if started {
                    r = r.square();
                }
                if (e[i] >> bit) & 1 == 1 {
                    r = if started { r.mul(self) } else { *self };
                    started = true;
                }
            }
        }
        r
    }

    /// Square root (Adj–Rodríguez-Henríquez, Algorithm 9, q = p^2,
    /// p ≡ 3 mod 4), constant time: both branches are computed and the
    /// result selected.  Returns (root, is_square).  sqrt(0) = 0.
    pub fn sqrt(&self) -> (Self, Choice) {
        let a1 = self.pow_public(&C::EXP_P34); // a^((p-3)/4)
        let alpha = a1.square().mul(self); // a^((p-1)/2)
        let a0 = alpha.frob().mul(&alpha); // a^((p^2-1)/2) up to the norm map
        let minus_one = Self::ONE.neg();
        let nonsquare = a0.ct_eq(&minus_one);
        let x0 = a1.mul(self);
        // branch 1: alpha == -1  ->  x = i * x0
        let x_neg = Self::i().mul(&x0);
        // branch 2: x = (1 + alpha)^((p-1)/2) * x0
        let b = Self::ONE.add(&alpha).pow_public(&C::EXP_LEGENDRE);
        let x_pos = b.mul(&x0);
        let x = Self::conditional_select(&x_pos, &x_neg, alpha.ct_eq(&minus_one));
        let ok = !nonsquare & x.square().ct_eq(self);
        (x, ok)
    }

    pub fn to_bytes_be(&self) -> Vec<u8> {
        let mut v = self.c0.to_bytes_be();
        v.extend(self.c1.to_bytes_be());
        v
    }
    pub fn from_bytes_be(b: &[u8]) -> (Self, Choice) {
        if b.len() != 2 * C::BYTES {
            return (Self::ZERO, Choice::from(0));
        }
        let (x, y) = b.split_at(C::BYTES);
        let (c0, ok0) = Fp::from_bytes_be(x);
        let (c1, ok1) = Fp::from_bytes_be(y);
        (Fp2 { c0, c1 }, ok0 & ok1)
    }
    pub fn byte_len() -> usize {
        2 * C::BYTES
    }
}

/// Batch inversion (Montgomery's trick); zero inputs give zero outputs and do
/// not disturb the others.  Constant time in the values.
pub fn batch_invert<C: FpConfig<N>, const N: usize>(xs: &[Fp2<C, N>]) -> Vec<Fp2<C, N>> {
    let n = xs.len();
    let mut pre = Vec::with_capacity(n);
    let mut acc = Fp2::ONE;
    let mut zs = Vec::with_capacity(n);
    for x in xs {
        let z = x.is_zero();
        zs.push(z);
        let xx = Fp2::conditional_select(x, &Fp2::ONE, z);
        pre.push(acc);
        acc = acc.mul(&xx);
    }
    let mut inv = acc.invert();
    let mut out = vec![Fp2::ZERO; n];
    for k in (0..n).rev() {
        let xx = Fp2::conditional_select(&xs[k], &Fp2::ONE, zs[k]);
        let v = inv.mul(&pre[k]);
        out[k] = Fp2::conditional_select(&v, &Fp2::ZERO, zs[k]);
        inv = inv.mul(&xx);
    }
    out
}

impl<C: FpConfig<N>, const N: usize> ConditionallySelectable for Fp2<C, N> {
    #[inline(always)]
    fn conditional_select(a: &Self, b: &Self, c: Choice) -> Self {
        Fp2 {
            c0: Fp::conditional_select(&a.c0, &b.c0, c),
            c1: Fp::conditional_select(&a.c1, &b.c1, c),
        }
    }
}
impl<C: FpConfig<N>, const N: usize> ConstantTimeEq for Fp2<C, N> {
    #[inline(always)]
    fn ct_eq(&self, o: &Self) -> Choice {
        self.c0.ct_eq(&o.c0) & self.c1.ct_eq(&o.c1)
    }
}
impl<C: FpConfig<N>, const N: usize> PartialEq for Fp2<C, N> {
    fn eq(&self, o: &Self) -> bool {
        bool::from(self.ct_eq(o))
    }
}
impl<C: FpConfig<N>, const N: usize> Eq for Fp2<C, N> {}

fp_ops!(Fp2);
