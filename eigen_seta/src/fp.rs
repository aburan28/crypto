//! Constant-time arithmetic in F_p, p odd, in Montgomery form with N 64-bit limbs.
//!
//! All derived constants (R mod p, R^2 mod p, -p^{-1} mod 2^64 and the public
//! exponents) are computed at compile time from the modulus by `const fn`s.
//!
//! Constant-time contract: `add`, `sub`, `neg`, `mul`, `square`, `ct_eq`,
//! `is_zero`, `conditional_select`, and every `pow_public` call execute the
//! same instruction sequence and memory accesses for all field elements.
//! `pow_public` branches on bits of its *exponent*, which must be public
//! (it is only called with p-2, (p+1)/4, (p-1)/2, (p-3)/4 or other
//! parameter-derived constants).

use crate::ct::{choice_from_bit, mask, Choice, ConditionallySelectable, ConstantTimeEq};
use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use zeroize::Zeroize;

// ---------------------------------------------------------------------------
// compile-time helpers on [u64; N]
// ---------------------------------------------------------------------------

/// -p^{-1} mod 2^64 for odd p0.
pub const fn mont_pinv(p0: u64) -> u64 {
    let mut inv: u64 = 1;
    let mut i = 0;
    while i < 7 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(p0.wrapping_mul(inv)));
        i += 1;
    }
    inv.wrapping_neg()
}

const fn geq<const N: usize>(a: &[u64; N], b: &[u64; N]) -> bool {
    let mut i = N;
    while i > 0 {
        i -= 1;
        if a[i] > b[i] {
            return true;
        }
        if a[i] < b[i] {
            return false;
        }
    }
    true
}

const fn sub_const<const N: usize>(a: &[u64; N], b: &[u64; N]) -> [u64; N] {
    let mut r = [0u64; N];
    let mut borrow = 0u64;
    let mut i = 0;
    while i < N {
        let (d1, b1) = a[i].overflowing_sub(b[i]);
        let (d2, b2) = d1.overflowing_sub(borrow);
        r[i] = d2;
        borrow = (b1 | b2) as u64;
        i += 1;
    }
    r
}

const fn double_mod<const N: usize>(a: &[u64; N], p: &[u64; N]) -> [u64; N] {
    let mut r = [0u64; N];
    let mut carry = 0u64;
    let mut i = 0;
    while i < N {
        r[i] = (a[i] << 1) | carry;
        carry = a[i] >> 63;
        i += 1;
    }
    if carry == 1 || geq(&r, p) {
        r = sub_const(&r, p);
    }
    r
}

/// 2^(64 N k) mod p.
pub const fn r_pow_mod_p<const N: usize>(p: &[u64; N], k: usize) -> [u64; N] {
    let mut r = [0u64; N];
    r[0] = 1;
    let mut i = 0;
    while i < 64 * N * k {
        r = double_mod(&r, p);
        i += 1;
    }
    r
}

/// a - small (a > small).
pub const fn sub_small<const N: usize>(a: &[u64; N], s: u64) -> [u64; N] {
    let mut b = [0u64; N];
    b[0] = s;
    sub_const(a, &b)
}

/// a + small.
pub const fn add_small<const N: usize>(a: &[u64; N], s: u64) -> [u64; N] {
    let mut r = *a;
    let mut carry = s;
    let mut i = 0;
    while i < N {
        let (v, c) = r[i].overflowing_add(carry);
        r[i] = v;
        carry = c as u64;
        i += 1;
    }
    r
}

/// a >> k for k < 64.
pub const fn shr_const<const N: usize>(a: &[u64; N], k: u32) -> [u64; N] {
    let mut r = [0u64; N];
    let mut i = 0;
    while i < N {
        r[i] = a[i] >> k;
        if k > 0 && i + 1 < N {
            r[i] |= a[i + 1] << (64 - k);
        }
        i += 1;
    }
    r
}

pub const fn bit_length<const N: usize>(a: &[u64; N]) -> usize {
    let mut i = N;
    while i > 0 {
        i -= 1;
        if a[i] != 0 {
            return 64 * i + (64 - a[i].leading_zeros() as usize);
        }
    }
    0
}

// ---------------------------------------------------------------------------
// field configuration
// ---------------------------------------------------------------------------

/// A prime field p ≡ 3 (mod 4).  Implementors give `MODULUS`; everything
/// else has compile-time defaults.
pub trait FpConfig<const N: usize>: 'static + Copy + Clone + Send + Sync + fmt::Debug {
    const MODULUS: [u64; N];
    const NAME: &'static str;
    const PINV: u64 = mont_pinv(Self::MODULUS[0]);
    const R: [u64; N] = r_pow_mod_p(&Self::MODULUS, 1);
    const R2: [u64; N] = r_pow_mod_p(&Self::MODULUS, 2);
    /// p - 2 (Fermat inversion).
    const EXP_INV: [u64; N] = sub_small(&Self::MODULUS, 2);
    /// (p - 1) / 2 (Euler criterion).
    const EXP_LEGENDRE: [u64; N] = shr_const(&sub_small(&Self::MODULUS, 1), 1);
    /// (p + 1) / 4 (square root, p ≡ 3 mod 4).
    const EXP_SQRT: [u64; N] = shr_const(&add_small(&Self::MODULUS, 1), 2);
    /// (p - 3) / 4 (F_{p^2} square root).
    const EXP_P34: [u64; N] = shr_const(&sub_small(&Self::MODULUS, 3), 2);
    const BITS: usize = bit_length(&Self::MODULUS);
    /// R^3 mod p (to return safegcd output to Montgomery form).
    const R3: [u64; N] = r_pow_mod_p(&Self::MODULUS, 3);
    const BYTES: usize = (bit_length(&Self::MODULUS) + 7) / 8;
    /// Number of low 64-bit limbs of p + 1 that are zero ("Montgomery-friendly"
    /// primes c*2^a*3^b - 1).  When positive, p = -1 mod 2^64 and the
    /// reduction skips those limbs.
    const P1_ZERO_LIMBS: usize = p1_zero_limbs(&Self::MODULUS);
}

const fn p1_zero_limbs<const N: usize>(p: &[u64; N]) -> usize {
    let mut z = 0;
    while z < N && p[z] == u64::MAX {
        z += 1;
    }
    z
}

/// Largest limb count with a stack product buffer (L5 uses 42).
pub const MAX_LIMBS: usize = 48;

/// Element of F_p in Montgomery form.
pub struct Fp<C: FpConfig<N>, const N: usize> {
    pub(crate) l: [u64; N],
    _c: PhantomData<C>,
}

impl<C: FpConfig<N>, const N: usize> Clone for Fp<C, N> {
    #[inline(always)]
    fn clone(&self) -> Self {
        *self
    }
}
impl<C: FpConfig<N>, const N: usize> Copy for Fp<C, N> {}

impl<C: FpConfig<N>, const N: usize> Zeroize for Fp<C, N> {
    fn zeroize(&mut self) {
        self.l.zeroize();
    }
}

impl<C: FpConfig<N>, const N: usize> fmt::Debug for Fp<C, N> {
    /// Debug output reveals the value; never format secret elements.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let b = self.to_bytes_be();
        write!(f, "Fp(0x")?;
        for x in b {
            write!(f, "{:02x}", x)?;
        }
        write!(f, ")")
    }
}

#[inline(always)]
fn mac(a: u64, b: u64, c: u64, carry: u64) -> (u64, u64) {
    // a + b*c + carry, never overflows u128
    let t = (a as u128) + (b as u128) * (c as u128) + (carry as u128);
    (t as u64, (t >> 64) as u64)
}

#[inline(always)]
fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let t = (a as u128) + (b as u128) + (carry as u128);
    (t as u64, (t >> 64) as u64)
}

#[inline(always)]
fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let t = (a as u128).wrapping_sub((b as u128) + (borrow as u128));
    (t as u64, ((t >> 64) as u64) & 1)
}

impl<C: FpConfig<N>, const N: usize> Fp<C, N> {
    #[inline(always)]
    pub(crate) const fn from_raw(l: [u64; N]) -> Self {
        Fp { l, _c: PhantomData }
    }

    pub const ZERO: Self = Self::from_raw([0u64; N]);
    pub const ONE: Self = Self::from_raw(C::R);

    #[inline(always)]
    pub fn zero() -> Self {
        Self::ZERO
    }
    #[inline(always)]
    pub fn one() -> Self {
        Self::ONE
    }

    /// r = t - p if t >= p (or `hi` set), else t; branch-free.
    #[inline(always)]
    fn reduce_once(t: [u64; N], hi: u64) -> [u64; N] {
        let mut d = [0u64; N];
        let mut borrow = 0u64;
        for i in 0..N {
            let (v, b) = sbb(t[i], C::MODULUS[i], borrow);
            d[i] = v;
            borrow = b;
        }
        // keep d iff hi == 1 or no borrow
        let (_, b2) = sbb(hi, 0, borrow);
        let use_t = choice_from_bit(b2); // borrow out of the (hi, t) - p subtraction
        let m = mask(use_t);
        let mut r = [0u64; N];
        for i in 0..N {
            r[i] = (t[i] & m) | (d[i] & !m);
        }
        r
    }

    #[inline]
    pub fn add(self, o: &Self) -> Self {
        let mut t = [0u64; N];
        let mut carry = 0u64;
        for i in 0..N {
            let (v, c) = adc(self.l[i], o.l[i], carry);
            t[i] = v;
            carry = c;
        }
        Self::from_raw(Self::reduce_once(t, carry))
    }

    #[inline]
    pub fn sub(self, o: &Self) -> Self {
        let mut t = [0u64; N];
        let mut borrow = 0u64;
        for i in 0..N {
            let (v, b) = sbb(self.l[i], o.l[i], borrow);
            t[i] = v;
            borrow = b;
        }
        // add p back iff borrow
        let m = 0u64.wrapping_sub(borrow);
        let m = core::hint::black_box(m);
        let mut carry = 0u64;
        for i in 0..N {
            let (v, c) = adc(t[i], C::MODULUS[i] & m, carry);
            t[i] = v;
            carry = c;
        }
        Self::from_raw(t)
    }

    #[inline]
    pub fn neg(self) -> Self {
        Self::ZERO.sub(&self)
    }

    #[inline]
    pub fn double(&self) -> Self {
        self.add(self)
    }

    /// Montgomery product (CIOS), branch-free.
    #[inline]
    pub fn mul(self, o: &Self) -> Self {
        if C::P1_ZERO_LIMBS > 0 && N <= MAX_LIMBS {
            self.mul_friendly(o)
        } else {
            self.mul_cios(o)
        }
    }

    /// Operand-scanning product then Montgomery reduction specialised to
    /// p + 1 = 0 mod 2^(64 Z): with p = -1 mod 2^64 the quotient digit is
    /// m = t[i] and t + m p 2^(64i) = t - m 2^(64i) + m (p+1) 2^(64i), where
    /// the first term clears t[i] exactly and the second touches only limbs
    /// i+Z .. i+N-1.  Cost N^2 + N (N - Z) word products instead of 2 N^2.
    /// All loop bounds are compile-time constants: constant time.
    #[inline]
    fn mul_friendly(self, o: &Self) -> Self {
        let a = &self.l;
        let b = &o.l;
        let p = &C::MODULUS;
        let z = C::P1_ZERO_LIMBS;
        let mut t = [0u64; 2 * MAX_LIMBS + 1];
        for i in 0..N {
            let mut c = 0u64;
            for j in 0..N {
                let (v, cc) = mac(t[i + j], a[j], b[i], c);
                t[i + j] = v;
                c = cc;
            }
            t[i + N] = c;
        }
        // p + 1 limbs: zero below z, p[z] + 1 at z (p[z] != MAX), p[j] above.
        let p1z = p[z % N].wrapping_add(1);
        let mut pending = 0u64; // carry into t[i + N] from the previous row
        for i in 0..N {
            let m = t[i];
            t[i] = 0;
            let mut c = 0u64;
            for j in z..N {
                let pj = if j == z { p1z } else { p[j] };
                let (v, cc) = mac(t[i + j], m, pj, c);
                t[i + j] = v;
                c = cc;
            }
            let (v, c1) = adc(t[i + N], c, pending);
            t[i + N] = v;
            pending = c1;
        }
        let mut r = [0u64; N];
        r.copy_from_slice(&t[N..2 * N]);
        Self::from_raw(Self::reduce_once(r, pending))
    }

    fn mul_cios(self, o: &Self) -> Self {
        let a = &self.l;
        let b = &o.l;
        let p = &C::MODULUS;
        let mut t = [0u64; N];
        let mut tn: u64 = 0;
        for i in 0..N {
            let mut c = 0u64;
            for j in 0..N {
                let (v, cc) = mac(t[j], a[j], b[i], c);
                t[j] = v;
                c = cc;
            }
            let (v, tn1) = adc(tn, c, 0);
            tn = v;
            let m = t[0].wrapping_mul(C::PINV);
            let (_, mut c) = mac(t[0], m, p[0], 0);
            for j in 1..N {
                let (v, cc) = mac(t[j], m, p[j], c);
                t[j - 1] = v;
                c = cc;
            }
            let (v, cc) = adc(tn, c, 0);
            t[N - 1] = v;
            tn = tn1 + cc;
        }
        Self::from_raw(Self::reduce_once(t, tn))
    }

    #[inline]
    pub fn square(&self) -> Self {
        self.mul(self)
    }

    /// x^e for a *public* exponent given as little-endian limbs.
    /// Branches only on exponent bits.
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

    /// 1/x via Fermat; maps 0 to 0.  Constant time.
    pub fn invert_fermat(&self) -> Self {
        self.pow_public(&C::EXP_INV)
    }

    /// 1/x via constant-time safegcd; maps 0 to 0.  The Montgomery value
    /// x R is inverted as an integer, then multiplied by R^3 (Montgomery
    /// product -> x^{-1} R).
    pub fn invert(&self) -> Self {
        let inv = crate::safegcd::modinv(&self.l, &C::MODULUS, C::PINV, C::BITS);
        let mut l = [0u64; N];
        l.copy_from_slice(&inv);
        Self::from_raw(l).mul(&Self::from_raw(C::R3))
    }

    /// Square root candidate x^((p+1)/4) and whether it is a root.
    pub fn sqrt(&self) -> (Self, Choice) {
        let r = self.pow_public(&C::EXP_SQRT);
        let ok = r.square().ct_eq(self);
        (r, ok)
    }

    /// Euler criterion: true for squares (including 0).
    pub fn is_square(&self) -> Choice {
        let l = self.pow_public(&C::EXP_LEGENDRE);
        l.ct_eq(&Self::ONE) | self.is_zero()
    }

    pub fn is_zero(&self) -> Choice {
        let mut acc = 0u64;
        for i in 0..N {
            acc |= self.l[i];
        }
        crate::ct::is_zero_u64(acc)
    }

    /// Small integer (may be negative) into the field.
    pub fn from_i64(v: i64) -> Self {
        let r = Self::from_u64(v.unsigned_abs());
        let neg = choice_from_bit((v < 0) as u64);
        Self::conditional_select(&r, &r.neg(), neg)
    }

    pub fn from_u64(v: u64) -> Self {
        let mut l = [0u64; N];
        l[0] = v;
        // v < p always holds for our parameter sets (p > 2^64)
        Self::from_raw(l).mul(&Self::from_raw(C::R2))
    }

    /// From canonical little-endian limbs (value must be < p).
    pub fn from_canonical_limbs(l: [u64; N]) -> Self {
        Self::from_raw(l).mul(&Self::from_raw(C::R2))
    }

    /// Canonical little-endian limbs (out of Montgomery form).
    pub fn to_canonical_limbs(&self) -> [u64; N] {
        let mut one = [0u64; N];
        one[0] = 1;
        self.mul(&Self::from_raw(one)).l
    }

    /// Big-endian encoding of the canonical value on `C::BYTES` bytes.
    pub fn to_bytes_be(&self) -> Vec<u8> {
        let l = self.to_canonical_limbs();
        let mut out = vec![0u8; C::BYTES];
        for (k, byte) in out.iter_mut().rev().enumerate() {
            *byte = (l[k / 8] >> (8 * (k % 8))) as u8;
        }
        out
    }

    /// Decode big-endian bytes; the Choice is false if the value is >= p or
    /// the length is wrong (then the element is zero).  Constant time in the
    /// value.
    pub fn from_bytes_be(b: &[u8]) -> (Self, Choice) {
        if b.len() != C::BYTES {
            return (Self::ZERO, Choice::from(0));
        }
        let mut l = [0u64; N];
        for (k, byte) in b.iter().rev().enumerate() {
            l[k / 8] |= (*byte as u64) << (8 * (k % 8));
        }
        // canonical iff l < p
        let mut borrow = 0u64;
        for i in 0..N {
            let (_, bb) = sbb(l[i], C::MODULUS[i], borrow);
            borrow = bb;
        }
        let ok = choice_from_bit(borrow);
        let v = Self::from_canonical_limbs(l);
        (Self::conditional_select(&Self::ZERO, &v, ok), ok)
    }

    /// Reduce 2*BYTES uniformly random bytes to a field element (bias
    /// < 2^-(8*BYTES)); used for sampling.
    pub fn from_wide_bytes(b: &[u8]) -> Self {
        assert_eq!(b.len(), 2 * C::BYTES);
        let (hi, lo) = b.split_at(C::BYTES);
        let mut acc = Self::ZERO;
        let base = Self::from_u64(256);
        for &x in hi.iter().chain(lo.iter()) {
            acc = acc.mul(&base).add(&Self::from_u64(x as u64));
        }
        acc
    }
}

impl<C: FpConfig<N>, const N: usize> ConditionallySelectable for Fp<C, N> {
    #[inline(always)]
    fn conditional_select(a: &Self, b: &Self, c: Choice) -> Self {
        let m = mask(c);
        let m = core::hint::black_box(m);
        let mut l = [0u64; N];
        for i in 0..N {
            l[i] = (a.l[i] & !m) | (b.l[i] & m);
        }
        Self::from_raw(l)
    }
}

impl<C: FpConfig<N>, const N: usize> ConstantTimeEq for Fp<C, N> {
    #[inline(always)]
    fn ct_eq(&self, o: &Self) -> Choice {
        let mut acc = 0u64;
        for i in 0..N {
            acc |= self.l[i] ^ o.l[i];
        }
        crate::ct::is_zero_u64(acc)
    }
}

impl<C: FpConfig<N>, const N: usize> PartialEq for Fp<C, N> {
    /// Constant-time comparison exposed as `==` for test convenience.
    fn eq(&self, o: &Self) -> bool {
        bool::from(self.ct_eq(o))
    }
}
impl<C: FpConfig<N>, const N: usize> Eq for Fp<C, N> {}

macro_rules! fp_ops {
    ($T:ident) => {
        impl<C: FpConfig<N>, const N: usize> Add for $T<C, N> {
            type Output = Self;
            #[inline(always)]
            fn add(self, o: Self) -> Self {
                $T::add(self, &o)
            }
        }
        impl<C: FpConfig<N>, const N: usize> Sub for $T<C, N> {
            type Output = Self;
            #[inline(always)]
            fn sub(self, o: Self) -> Self {
                $T::sub(self, &o)
            }
        }
        impl<C: FpConfig<N>, const N: usize> Mul for $T<C, N> {
            type Output = Self;
            #[inline(always)]
            fn mul(self, o: Self) -> Self {
                $T::mul(self, &o)
            }
        }
        impl<C: FpConfig<N>, const N: usize> Neg for $T<C, N> {
            type Output = Self;
            #[inline(always)]
            fn neg(self) -> Self {
                $T::neg(self)
            }
        }
        impl<C: FpConfig<N>, const N: usize> AddAssign for $T<C, N> {
            #[inline(always)]
            fn add_assign(&mut self, o: Self) {
                *self = $T::add(*self, &o);
            }
        }
        impl<C: FpConfig<N>, const N: usize> SubAssign for $T<C, N> {
            #[inline(always)]
            fn sub_assign(&mut self, o: Self) {
                *self = $T::sub(*self, &o);
            }
        }
        impl<C: FpConfig<N>, const N: usize> MulAssign for $T<C, N> {
            #[inline(always)]
            fn mul_assign(&mut self, o: Self) {
                *self = $T::mul(*self, &o);
            }
        }
    };
}
pub(crate) use fp_ops;
fp_ops!(Fp);

/// Parse a big-endian hex string (no prefix) into little-endian limbs at
/// compile time.
pub const fn hex_to_limbs<const N: usize>(s: &str) -> [u64; N] {
    let b = s.as_bytes();
    let mut l = [0u64; N];
    let mut k = 0; // nibble index from the right
    let mut i = b.len();
    while i > 0 {
        i -= 1;
        let c = b[i];
        let v = if c >= b'0' && c <= b'9' {
            c - b'0'
        } else if c >= b'a' && c <= b'f' {
            c - b'a' + 10
        } else if c >= b'A' && c <= b'F' {
            c - b'A' + 10
        } else {
            panic!("bad hex digit")
        } as u64;
        assert!(k / 16 < N, "hex constant does not fit");
        l[k / 16] |= v << (4 * (k % 16));
        k += 1;
    }
    l
}

#[macro_export]
macro_rules! hexlimbs {
    ($s:expr, $n:expr) => {
        $crate::fp::hex_to_limbs::<$n>($s)
    };
}
