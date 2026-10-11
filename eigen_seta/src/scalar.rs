//! Fixed-width integers for secret scalars, with constant-time arithmetic
//! modulo 2^k and modulo a public odd modulus (3^f).  Widths and moduli are
//! public; values may be secret.

use crate::ct::{choice_from_bit, mask, Choice, ConditionallySelectable, ConstantTimeEq};
use num_bigint::BigUint;
use zeroize::{Zeroize, ZeroizeOnDrop};

/// Little-endian limbs, fixed length.
#[derive(Clone, Debug, Zeroize, ZeroizeOnDrop)]
pub struct Scalar {
    pub l: Vec<u64>,
}

impl Scalar {
    pub fn zero(limbs: usize) -> Self {
        Scalar { l: vec![0u64; limbs] }
    }
    pub fn from_u64(v: u64, limbs: usize) -> Self {
        let mut s = Self::zero(limbs);
        s.l[0] = v;
        s
    }
    pub fn limbs(&self) -> usize {
        self.l.len()
    }
    /// Bit i (public index).
    #[inline(always)]
    pub fn bit(&self, i: usize) -> u64 {
        if i / 64 >= self.l.len() {
            return 0;
        }
        (self.l[i / 64] >> (i % 64)) & 1
    }
    pub fn from_biguint(b: &BigUint, limbs: usize) -> Self {
        let mut s = Self::zero(limbs);
        for (i, d) in b.to_u64_digits().iter().enumerate() {
            assert!(i < limbs, "value does not fit");
            s.l[i] = *d;
        }
        s
    }
    /// Variable-time conversion; use only for public values or in keygen.
    pub fn to_biguint(&self) -> BigUint {
        let mut b = BigUint::default();
        for i in (0..self.l.len()).rev() {
            b = (b << 64u32) + BigUint::from(self.l[i]);
        }
        b
    }
    /// Little-endian bytes into a scalar of `limbs` limbs (input must fit).
    pub fn from_le_bytes(b: &[u8], limbs: usize) -> Self {
        let mut s = Self::zero(limbs);
        for (k, byte) in b.iter().enumerate() {
            if k / 8 < limbs {
                s.l[k / 8] |= (*byte as u64) << (8 * (k % 8));
            }
        }
        s
    }
    pub fn to_le_bytes(&self, nbytes: usize) -> Vec<u8> {
        (0..nbytes).map(|k| (self.l[k / 8] >> (8 * (k % 8))) as u8).collect()
    }

    /// Reduce modulo 2^k (k public).
    pub fn mod_2k(&self, k: usize) -> Self {
        let mut r = self.clone();
        for i in 0..r.l.len() {
            let lo = 64 * i;
            if lo >= k {
                r.l[i] = 0;
            } else if lo + 64 > k {
                r.l[i] &= (1u64 << (k - lo)) - 1;
            }
        }
        r
    }

    pub fn add_mod_2k(&self, o: &Self, k: usize) -> Self {
        let mut r = Self::zero(self.l.len());
        let mut c = 0u128;
        for i in 0..self.l.len() {
            let t = self.l[i] as u128 + o.l[i] as u128 + c;
            r.l[i] = t as u64;
            c = t >> 64;
        }
        r.mod_2k(k)
    }

    pub fn neg_mod_2k(&self, k: usize) -> Self {
        // 2^k - x = (~x + 1) mod 2^k
        let mut r = Self::zero(self.l.len());
        let mut c = 1u128;
        for i in 0..self.l.len() {
            let t = (!self.l[i]) as u128 + c;
            r.l[i] = t as u64;
            c = t >> 64;
        }
        r.mod_2k(k)
    }

    pub fn sub_mod_2k(&self, o: &Self, k: usize) -> Self {
        self.add_mod_2k(&o.neg_mod_2k(k), k)
    }

    /// Truncated product modulo 2^k; constant time (schoolbook over all limbs).
    pub fn mul_mod_2k(&self, o: &Self, k: usize) -> Self {
        let n = self.l.len();
        let mut r = vec![0u64; n];
        for i in 0..n {
            let mut c = 0u128;
            for j in 0..(n - i) {
                let t = r[i + j] as u128 + (self.l[i] as u128) * (o.l[j] as u128) + c;
                r[i + j] = t as u64;
                c = t >> 64;
            }
        }
        Scalar { l: r }.mod_2k(k)
    }

    /// Inverse of an odd x modulo 2^k by Newton iteration (fixed count).
    pub fn inv_mod_2k(&self, k: usize) -> Self {
        let n = self.l.len();
        let two = Self::from_u64(2, n);
        let mut y = Self::from_u64(1, n); // correct mod 2
        let mut prec = 1usize;
        while prec < k {
            // y <- y (2 - x y)
            let xy = self.mul_mod_2k(&y, k);
            let t = two.sub_mod_2k(&xy, k);
            y = y.mul_mod_2k(&t, k);
            prec *= 2;
        }
        y
    }

    pub fn is_odd(&self) -> Choice {
        choice_from_bit(self.l[0] & 1)
    }

    pub fn ct_is_zero(&self) -> Choice {
        let mut acc = 0u64;
        for x in &self.l {
            acc |= *x;
        }
        crate::ct::is_zero_u64(acc)
    }
}

impl ConstantTimeEq for Scalar {
    fn ct_eq(&self, o: &Self) -> Choice {
        assert_eq!(self.l.len(), o.l.len());
        let mut acc = 0u64;
        for i in 0..self.l.len() {
            acc |= self.l[i] ^ o.l[i];
        }
        crate::ct::is_zero_u64(acc)
    }
}

impl Scalar {
    pub fn select(a: &Self, b: &Self, c: Choice) -> Self {
        assert_eq!(a.l.len(), b.l.len());
        let m = core::hint::black_box(mask(c));
        Scalar { l: a.l.iter().zip(&b.l).map(|(x, y)| (x & !m) | (y & m)).collect() }
    }
}

/// Constant-time arithmetic modulo a public odd modulus (runtime width).
#[derive(Clone, Debug)]
pub struct OddModulus {
    pub m: Vec<u64>,
    minv: u64,
    r2: Vec<u64>,
    pub bits: usize,
}

impl OddModulus {
    pub fn new(m: &BigUint) -> Self {
        let n = m.to_u64_digits().len();
        let mut ml = vec![0u64; n];
        for (i, d) in m.to_u64_digits().iter().enumerate() {
            ml[i] = *d;
        }
        let minv = crate::fp::mont_pinv(ml[0]);
        let r2b = (BigUint::from(1u32) << (128 * n)) % m;
        let mut r2 = vec![0u64; n];
        for (i, d) in r2b.to_u64_digits().iter().enumerate() {
            r2[i] = *d;
        }
        OddModulus { m: ml, minv, r2, bits: m.bits() as usize }
    }
    pub fn limbs(&self) -> usize {
        self.m.len()
    }

    fn mont_mul(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        let n = self.m.len();
        let mut t = vec![0u64; n];
        let mut tn = 0u64;
        for i in 0..n {
            let mut c = 0u128;
            for j in 0..n {
                let v = t[j] as u128 + (a[j] as u128) * (b[i] as u128) + c;
                t[j] = v as u64;
                c = v >> 64;
            }
            let v = tn as u128 + c;
            tn = v as u64;
            let tn1 = (v >> 64) as u64;
            let mm = t[0].wrapping_mul(self.minv);
            let v = t[0] as u128 + (mm as u128) * (self.m[0] as u128);
            let mut c = v >> 64;
            for j in 1..n {
                let v = t[j] as u128 + (mm as u128) * (self.m[j] as u128) + c;
                t[j - 1] = v as u64;
                c = v >> 64;
            }
            let v = tn as u128 + c;
            t[n - 1] = v as u64;
            tn = tn1 + (v >> 64) as u64;
        }
        // conditional subtraction
        let mut d = vec![0u64; n];
        let mut borrow = 0u64;
        for i in 0..n {
            let v = (t[i] as u128).wrapping_sub(self.m[i] as u128 + borrow as u128);
            d[i] = v as u64;
            borrow = ((v >> 64) as u64) & 1;
        }
        let v = (tn as u128).wrapping_sub(borrow as u128);
        let use_t = choice_from_bit(((v >> 64) as u64) & 1);
        let msk = core::hint::black_box(mask(use_t));
        (0..n).map(|i| (t[i] & msk) | (d[i] & !msk)).collect()
    }

    /// Reduce a scalar (any width, value < 2^(64*limbs)) to [0, m) via
    /// Montgomery: x * R^2 * R^-1 ... applied limb-block-wise.  Constant time.
    pub fn reduce(&self, x: &Scalar) -> Scalar {
        let n = self.m.len();
        // process x in n-limb blocks from the top: acc = acc * 2^(64n) + block
        let mut acc = vec![0u64; n];
        let blocks = x.l.len().div_ceil(n);
        // R mod m in Montgomery form is R^2 mont-reduced once... compute shift = R mod m
        for bi in (0..blocks).rev() {
            let mut blk = vec![0u64; n];
            for j in 0..n {
                let idx = bi * n + j;
                if idx < x.l.len() {
                    blk[j] = x.l[idx];
                }
            }
            // acc <- acc * R + blk  (mod m), working in normal form:
            // acc*R mod m = mont_mul(acc, R^2)
            let shifted = self.mont_mul(&acc, &self.r2);
            // blk may be >= m: reduce blk = mont_mul(mont_mul(blk, R^2), 1)
            let mut one = vec![0u64; n];
            one[0] = 1;
            let blk_r = self.mont_mul(&self.mont_mul(&blk, &self.r2), &one);
            acc = self.add(&shifted, &blk_r);
        }
        Scalar { l: acc }
    }

    pub fn add(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        let n = self.m.len();
        let mut t = vec![0u64; n];
        let mut c = 0u128;
        for i in 0..n {
            let v = a[i] as u128 + b[i] as u128 + c;
            t[i] = v as u64;
            c = v >> 64;
        }
        let mut d = vec![0u64; n];
        let mut borrow = 0u128;
        for i in 0..n {
            let v = (t[i] as u128).wrapping_sub(self.m[i] as u128 + borrow);
            d[i] = v as u64;
            borrow = (v >> 64) & 1;
        }
        let v = c.wrapping_sub(borrow);
        let use_t = choice_from_bit(((v >> 64) as u64) & 1);
        let msk = core::hint::black_box(mask(use_t));
        (0..n).map(|i| (t[i] & msk) | (d[i] & !msk)).collect()
    }

    /// a * b mod m for reduced a, b.
    pub fn mul(&self, a: &Scalar, b: &Scalar) -> Scalar {
        let ar = self.mont_mul(&a.l, &self.r2);
        Scalar { l: self.mont_mul(&ar, &b.l) }
    }

    /// a^e mod m, e public.
    pub fn pow_public(&self, a: &Scalar, e: &BigUint) -> Scalar {
        let n = self.m.len();
        let mut one = vec![0u64; n];
        one[0] = 1;
        let mut r = self.mont_mul(&one, &self.r2); // R mod m
        let am = self.mont_mul(&a.l, &self.r2);
        for i in (0..e.bits()).rev() {
            r = self.mont_mul(&r, &r);
            if e.bit(i) {
                r = self.mont_mul(&r, &am);
            }
        }
        Scalar { l: self.mont_mul(&r, &one) }
    }

    pub fn to_scalar_width(&self, x: &Scalar, limbs: usize) -> Scalar {
        let mut s = Scalar::zero(limbs);
        for i in 0..x.l.len().min(limbs) {
            s.l[i] = x.l[i];
        }
        s
    }
}

impl ConditionallySelectable for crate::scalar::Wrap {
    fn conditional_select(a: &Self, b: &Self, c: Choice) -> Self {
        Wrap(u64::conditional_select(&a.0, &b.0, c))
    }
}
/// Copyable word wrapper used for CT table selection of small digits.
#[derive(Clone, Copy, Debug, Default)]
pub struct Wrap(pub u64);
