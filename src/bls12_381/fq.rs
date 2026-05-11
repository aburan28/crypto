//! Base field `F_p` for BLS12-381.
//!
//! `p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab`
//! (381 bits).

use num_bigint::BigUint;
use num_traits::{One, Zero};

/// BLS12-381 base-field prime `p`.
pub fn modulus() -> BigUint {
    BigUint::parse_bytes(
        b"1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
        16,
    )
    .unwrap()
}

/// BLS12-381 scalar-field prime `r` (group order on G1, G2).
pub fn scalar_modulus() -> BigUint {
    BigUint::parse_bytes(
        b"73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
        16,
    )
    .unwrap()
}

/// An element of `F_p`.  Stored normalised in `[0, p)`.
#[derive(Clone, Debug, PartialEq)]
pub struct Fq {
    pub value: BigUint,
}

impl Fq {
    pub fn new(v: BigUint) -> Self {
        let p = modulus();
        Self { value: v % &p }
    }
    pub fn zero() -> Self {
        Self { value: BigUint::zero() }
    }
    pub fn one() -> Self {
        Self { value: BigUint::one() }
    }
    pub fn is_zero(&self) -> bool {
        self.value.is_zero()
    }
    pub fn add(&self, other: &Self) -> Self {
        Self::new(&self.value + &other.value)
    }
    pub fn sub(&self, other: &Self) -> Self {
        let p = modulus();
        Self::new(((&self.value + &p) - (&other.value % &p)) % &p)
    }
    pub fn neg(&self) -> Self {
        let p = modulus();
        if self.value.is_zero() {
            Self::zero()
        } else {
            Self { value: &p - &self.value }
        }
    }
    pub fn mul(&self, other: &Self) -> Self {
        Self::new(&self.value * &other.value)
    }
    pub fn square(&self) -> Self {
        self.mul(self)
    }
    pub fn inverse(&self) -> Option<Self> {
        let p = modulus();
        crate::utils::mod_inverse(&self.value, &p).map(|v| Self { value: v })
    }
    /// `self ^ exp` via square-and-multiply.
    pub fn pow(&self, exp: &BigUint) -> Self {
        Self {
            value: self.value.modpow(exp, &modulus()),
        }
    }
    /// Test whether `self` is a quadratic residue mod p via Euler's
    /// criterion `self^((p-1)/2) == 1`.
    pub fn is_qr(&self) -> bool {
        if self.is_zero() {
            return true;
        }
        let p = modulus();
        let exp = (&p - BigUint::one()) / BigUint::from(2u32);
        self.value.modpow(&exp, &p) == BigUint::one()
    }
    /// Modular square root using `p ≡ 3 (mod 4)`:
    /// `sqrt(a) = a^{(p+1)/4} (mod p)` when one exists.
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::zero());
        }
        if !self.is_qr() {
            return None;
        }
        let p = modulus();
        // BLS12-381 prime satisfies p ≡ 3 (mod 4): verified by p mod 4 = 3.
        let exp = (&p + BigUint::one()) / BigUint::from(4u32);
        Some(Self { value: self.value.modpow(&exp, &p) })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modulus_has_381_bits() {
        let p = modulus();
        assert_eq!(p.bits(), 381);
    }

    #[test]
    fn scalar_modulus_has_255_bits() {
        let r = scalar_modulus();
        assert_eq!(r.bits(), 255);
    }

    #[test]
    fn modulus_is_3_mod_4() {
        // Required for our sqrt() formula.
        let p = modulus();
        let four = BigUint::from(4u32);
        assert_eq!(&p % &four, BigUint::from(3u32));
    }

    #[test]
    fn add_sub_neg_basic() {
        let a = Fq::new(BigUint::from(10u32));
        let b = Fq::new(BigUint::from(3u32));
        let c = a.add(&b);
        assert_eq!(c.value, BigUint::from(13u32));
        let d = a.sub(&b);
        assert_eq!(d.value, BigUint::from(7u32));
        let e = b.sub(&a);
        let expected = modulus() - BigUint::from(7u32);
        assert_eq!(e.value, expected);
    }

    #[test]
    fn mul_square_consistency() {
        let a = Fq::new(BigUint::from(123456u64));
        let sq = a.square();
        let mm = a.mul(&a);
        assert_eq!(sq, mm);
    }

    #[test]
    fn inverse_roundtrip() {
        let a = Fq::new(BigUint::from(2u32).pow(200));
        let ai = a.inverse().unwrap();
        let prod = a.mul(&ai);
        assert_eq!(prod, Fq::one());
    }

    #[test]
    fn sqrt_roundtrip() {
        let a = Fq::new(BigUint::from(2u32).pow(100));
        let sq = a.square();
        let root = sq.sqrt().unwrap();
        // sqrt returns ±a; either way root² = sq.
        assert_eq!(root.square(), sq);
    }

    #[test]
    fn non_qr_returns_none() {
        // Find a non-QR by trial; about half are non-QR.
        let mut x = Fq::new(BigUint::from(2u32));
        for _ in 0..20 {
            if !x.is_qr() {
                assert!(x.sqrt().is_none());
                return;
            }
            x = x.add(&Fq::one());
        }
        panic!("could not find a non-QR in 20 trials (statistically improbable)");
    }
}
