//! F_p and F_{p^2} = F_p[i]/(i^2+1) arithmetic.
//! Requires p ≡ 3 (mod 4) so that -1 is a non-residue and i² + 1 is irreducible.

use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Zero};

// ---- F_p -----------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct Fp {
    pub p: BigInt,
}

impl Fp {
    pub fn new(p: BigInt) -> Self {
        assert!(&p % BigInt::from(4) == BigInt::from(3), "need p ≡ 3 (mod 4)");
        Fp { p }
    }
    pub fn r(&self, x: &BigInt) -> BigInt {
        let r = x % &self.p;
        if r.sign() == Sign::Minus { r + &self.p } else { r }
    }
    pub fn add(&self, a: &BigInt, b: &BigInt) -> BigInt { self.r(&(a + b)) }
    pub fn sub(&self, a: &BigInt, b: &BigInt) -> BigInt { self.r(&(a - b)) }
    pub fn mul(&self, a: &BigInt, b: &BigInt) -> BigInt { self.r(&(a * b)) }
    pub fn neg(&self, a: &BigInt) -> BigInt { self.r(&-a) }
    pub fn inv(&self, a: &BigInt) -> BigInt {
        let eg = self.r(a).extended_gcd(&self.p);
        assert!(eg.gcd.is_one(), "non-invertible element");
        self.r(&eg.x)
    }
    pub fn is_square(&self, a: &BigInt) -> bool {
        if a.is_zero() { return true; }
        a.modpow(&((&self.p - 1) / 2), &self.p).is_one()
    }
    pub fn sqrt(&self, a: &BigInt) -> Option<BigInt> {
        if a.is_zero() { return Some(BigInt::zero()); }
        if !self.is_square(a) { return None; }
        Some(a.modpow(&((&self.p + 1) / 4), &self.p))
    }
}

// ---- F_{p^2} = F_p[i]/(i^2+1) --------------------------------------------

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct F2 {
    pub a: BigInt,  // real part
    pub b: BigInt,  // imag part (coefficient of i)
}

#[derive(Clone)]
pub struct Fp2 {
    pub fp: Fp,
}

impl Fp2 {
    pub fn new(fp: Fp) -> Self { Fp2 { fp } }
    pub fn zero(&self) -> F2 { F2 { a: BigInt::zero(), b: BigInt::zero() } }
    pub fn one(&self) -> F2 { F2 { a: BigInt::one(), b: BigInt::zero() } }
    pub fn from_int(&self, n: i64) -> F2 {
        F2 { a: self.fp.r(&BigInt::from(n)), b: BigInt::zero() }
    }
    pub fn from_fp(&self, x: &BigInt) -> F2 {
        F2 { a: self.fp.r(x), b: BigInt::zero() }
    }
    pub fn i(&self) -> F2 { F2 { a: BigInt::zero(), b: BigInt::one() } }

    pub fn is_zero(&self, x: &F2) -> bool { x.a.is_zero() && x.b.is_zero() }

    pub fn add(&self, x: &F2, y: &F2) -> F2 {
        F2 { a: self.fp.add(&x.a, &y.a), b: self.fp.add(&x.b, &y.b) }
    }
    pub fn sub(&self, x: &F2, y: &F2) -> F2 {
        F2 { a: self.fp.sub(&x.a, &y.a), b: self.fp.sub(&x.b, &y.b) }
    }
    pub fn neg(&self, x: &F2) -> F2 {
        F2 { a: self.fp.neg(&x.a), b: self.fp.neg(&x.b) }
    }
    pub fn mul(&self, x: &F2, y: &F2) -> F2 {
        // (a+bi)(c+di) = (ac-bd) + (ad+bc)i
        let ac = self.fp.mul(&x.a, &y.a);
        let bd = self.fp.mul(&x.b, &y.b);
        let ad = self.fp.mul(&x.a, &y.b);
        let bc = self.fp.mul(&x.b, &y.a);
        F2 { a: self.fp.sub(&ac, &bd), b: self.fp.add(&ad, &bc) }
    }
    pub fn sq(&self, x: &F2) -> F2 { self.mul(x, x) }
    pub fn inv(&self, x: &F2) -> F2 {
        // (a+bi)^{-1} = (a - bi)/(a² + b²)
        let n = self.fp.add(&self.fp.mul(&x.a, &x.a), &self.fp.mul(&x.b, &x.b));
        let ni = self.fp.inv(&n);
        F2 {
            a: self.fp.mul(&x.a, &ni),
            b: self.fp.neg(&self.fp.mul(&x.b, &ni)),
        }
    }
    pub fn div(&self, x: &F2, y: &F2) -> F2 { self.mul(x, &self.inv(y)) }

    /// Multiply F_p² element by an integer scalar.
    pub fn mul_scalar(&self, x: &F2, n: i64) -> F2 {
        let s = self.fp.r(&BigInt::from(n));
        F2 { a: self.fp.mul(&x.a, &s), b: self.fp.mul(&x.b, &s) }
    }

    pub fn pow(&self, x: &F2, e: &BigInt) -> F2 {
        let mut r = self.one();
        let mut b = x.clone();
        let mut k = e.clone();
        while k.sign() == Sign::Plus {
            if k.is_odd() { r = self.mul(&r, &b); }
            b = self.sq(&b);
            k >>= 1;
        }
        r
    }

    /// Is x a square in F_{p²}?
    pub fn is_square(&self, x: &F2) -> bool {
        if self.is_zero(x) { return true; }
        let exp = (&self.fp.p * &self.fp.p - BigInt::one()) / 2;
        self.pow(x, &exp) == self.one()
    }

    /// Square root in F_{p²} for p ≡ 3 (mod 4). Decompose √(a + bi) = c + di
    /// via the norm: solve c² − d² = a, 2cd = b using F_p arithmetic.
    pub fn sqrt(&self, x: &F2) -> Option<F2> {
        if self.is_zero(x) { return Some(self.zero()); }
        if x.b.is_zero() {
            // Pure F_p element
            if let Some(rt) = self.fp.sqrt(&x.a) {
                return Some(F2 { a: rt, b: BigInt::zero() });
            }
            // x.a is a non-square in F_p → its square root is i · √(−x.a).
            let neg = self.fp.sub(&BigInt::zero(), &x.a);
            let rt = self.fp.sqrt(&neg)?;
            return Some(F2 { a: BigInt::zero(), b: rt });
        }
        let two_inv = self.fp.inv(&BigInt::from(2));
        let norm = self.fp.add(
            &self.fp.mul(&x.a, &x.a),
            &self.fp.mul(&x.b, &x.b),
        );
        let s = self.fp.sqrt(&norm)?;
        for sign in [s.clone(), self.fp.sub(&BigInt::zero(), &s)] {
            let half = self.fp.mul(&self.fp.add(&x.a, &sign), &two_inv);
            if let Some(c) = self.fp.sqrt(&half) {
                if c.is_zero() { continue; }
                let inv = self.fp.inv(&self.fp.mul(&BigInt::from(2), &c));
                let d = self.fp.mul(&x.b, &inv);
                let cand = F2 { a: c, b: d };
                if self.sq(&cand) == *x { return Some(cand); }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn ctx() -> Fp2 {
        Fp2::new(Fp::new(BigInt::from(431u64)))
    }

    #[test]
    fn fp_basic() {
        let fp = Fp::new(BigInt::from(431u64));
        let a = BigInt::from(200);
        let b = BigInt::from(300);
        assert_eq!(fp.add(&a, &b), BigInt::from(69));   // 500 mod 431
        assert_eq!(fp.mul(&a, &fp.inv(&a)), BigInt::one());
        let s = fp.sqrt(&BigInt::from(2)).unwrap();
        assert_eq!(fp.mul(&s, &s), BigInt::from(2));
    }

    #[test]
    fn fp2_basic() {
        let fp2 = ctx();
        let i = fp2.i();
        // i^2 = -1
        let i2 = fp2.sq(&i);
        let neg1 = fp2.neg(&fp2.one());
        assert_eq!(i2, neg1);
        // (3 + 4i)(3 - 4i) = 9 + 16 = 25
        let z = F2 { a: BigInt::from(3), b: BigInt::from(4) };
        let zbar = F2 { a: BigInt::from(3), b: fp2.fp.neg(&BigInt::from(4)) };
        assert_eq!(fp2.mul(&z, &zbar), fp2.from_int(25));
        // inverse
        assert_eq!(fp2.mul(&z, &fp2.inv(&z)), fp2.one());
    }
}
