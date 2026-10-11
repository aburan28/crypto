//! The field interface used by all curve code (implemented by `Fp2`).

use crate::ct::{Choice, ConditionallySelectable, ConstantTimeEq};
use crate::fp::FpConfig;
use crate::fp2::Fp2;
use core::fmt::Debug;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

pub trait Field:
    Copy
    + Debug
    + PartialEq
    + ConditionallySelectable
    + ConstantTimeEq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + zeroize::Zeroize
{
    fn zero() -> Self;
    fn one() -> Self;
    fn from_i64(v: i64) -> Self;
    fn square(&self) -> Self;
    /// 1/x, 0 -> 0, constant time.
    fn invert(&self) -> Self;
    /// (root, is_square), constant time.
    fn sqrt(&self) -> (Self, Choice);
    fn is_zero(&self) -> Choice;
    fn is_square(&self) -> Choice;
    fn pow_public(&self, e: &[u64]) -> Self;
    fn frob(&self) -> Self;
    fn to_bytes(&self) -> Vec<u8>;
    fn from_bytes(b: &[u8]) -> (Self, Choice);
    fn byte_len() -> usize;
    /// Characteristic p as little-endian limbs (public).
    fn char_limbs() -> Vec<u64>;
    /// Deterministic element from a counter (public data only).
    fn from_counter(k: u64) -> Self;

    #[inline(always)]
    fn double(&self) -> Self {
        *self + *self
    }
    #[inline(always)]
    fn mul_small(&self, k: i64) -> Self {
        *self * Self::from_i64(k)
    }
}

impl<C: FpConfig<N>, const N: usize> Field for Fp2<C, N> {
    fn zero() -> Self {
        Fp2::ZERO
    }
    fn one() -> Self {
        Fp2::ONE
    }
    fn from_i64(v: i64) -> Self {
        Fp2::from_i64(v)
    }
    fn square(&self) -> Self {
        Fp2::square(self)
    }
    fn invert(&self) -> Self {
        Fp2::invert(self)
    }
    fn sqrt(&self) -> (Self, Choice) {
        Fp2::sqrt(self)
    }
    fn is_zero(&self) -> Choice {
        Fp2::is_zero(self)
    }
    fn is_square(&self) -> Choice {
        Fp2::is_square(self)
    }
    fn pow_public(&self, e: &[u64]) -> Self {
        Fp2::pow_public(self, e)
    }
    fn frob(&self) -> Self {
        Fp2::frob(self)
    }
    fn to_bytes(&self) -> Vec<u8> {
        self.to_bytes_be()
    }
    fn from_bytes(b: &[u8]) -> (Self, Choice) {
        Fp2::from_bytes_be(b)
    }
    fn byte_len() -> usize {
        Fp2::<C, N>::byte_len()
    }
    fn char_limbs() -> Vec<u64> {
        C::MODULUS.to_vec()
    }
    fn from_counter(k: u64) -> Self {
        // i + k: never in F_p, so lifts hit both twists' x-coordinates evenly
        Fp2::new(crate::fp::Fp::from_u64(k), crate::fp::Fp::ONE)
    }
}

/// Select helper for any ConditionallySelectable.
#[inline(always)]
pub fn csel<T: ConditionallySelectable>(a: &T, b: &T, c: Choice) -> T {
    T::conditional_select(a, b, c)
}
