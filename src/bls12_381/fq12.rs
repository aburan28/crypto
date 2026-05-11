//! Quadratic extension `F_{p¹²} = F_{p⁶}[w]/(w² − v)`.
//!
//! Element: `c₀ + c₁·w` with `c_i ∈ F_{p⁶}`.  This is the
//! **target group** of the BLS12-381 pairing.

use super::fq6::Fq6;
use num_bigint::BigUint;
use num_traits::{One, Zero};

#[derive(Clone, Debug, PartialEq)]
pub struct Fq12 {
    pub c0: Fq6,
    pub c1: Fq6,
}

impl Fq12 {
    pub fn new(c0: Fq6, c1: Fq6) -> Self {
        Self { c0, c1 }
    }
    pub fn zero() -> Self {
        Self { c0: Fq6::zero(), c1: Fq6::zero() }
    }
    pub fn one() -> Self {
        Self { c0: Fq6::one(), c1: Fq6::zero() }
    }
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }
    pub fn add(&self, other: &Self) -> Self {
        Self { c0: self.c0.add(&other.c0), c1: self.c1.add(&other.c1) }
    }
    pub fn sub(&self, other: &Self) -> Self {
        Self { c0: self.c0.sub(&other.c0), c1: self.c1.sub(&other.c1) }
    }
    pub fn neg(&self) -> Self {
        Self { c0: self.c0.neg(), c1: self.c1.neg() }
    }
    /// Multiplication using `w² = v`:
    /// `(a + bw)(c + dw) = (ac + bd·v) + (ad + bc)·w`.
    pub fn mul(&self, other: &Self) -> Self {
        let ac = self.c0.mul(&other.c0);
        let bd = self.c1.mul(&other.c1);
        let ad = self.c0.mul(&other.c1);
        let bc = self.c1.mul(&other.c0);
        Self {
            c0: ac.add(&bd.mul_by_v()),
            c1: ad.add(&bc),
        }
    }
    pub fn square(&self) -> Self {
        self.mul(self)
    }
    /// Conjugate: `a + bw ↦ a − bw`.
    pub fn conjugate(&self) -> Self {
        Self { c0: self.c0.clone(), c1: self.c1.neg() }
    }
    pub fn inverse(&self) -> Option<Self> {
        // For (a + bw)(a − bw) = a² − b²·v in F_{p⁶}.
        let t = self.c0.square().sub(&self.c1.square().mul_by_v());
        let t_inv = t.inverse()?;
        Some(Self {
            c0: self.c0.mul(&t_inv),
            c1: self.c1.neg().mul(&t_inv),
        })
    }
    pub fn pow(&self, exp: &BigUint) -> Self {
        let mut result = Self::one();
        let mut base = self.clone();
        let mut e = exp.clone();
        while !e.is_zero() {
            if (&e & BigUint::one()) == BigUint::one() {
                result = result.mul(&base);
            }
            base = base.square();
            e >>= 1;
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::fq::Fq;
    use super::super::fq2::Fq2;
    use num_bigint::BigUint;

    fn sample() -> Fq12 {
        let a = Fq6::new(
            Fq2::new(Fq::new(BigUint::from(1u32)), Fq::new(BigUint::from(2u32))),
            Fq2::new(Fq::new(BigUint::from(3u32)), Fq::new(BigUint::from(4u32))),
            Fq2::new(Fq::new(BigUint::from(5u32)), Fq::new(BigUint::from(6u32))),
        );
        let b = Fq6::new(
            Fq2::new(Fq::new(BigUint::from(7u32)), Fq::new(BigUint::from(8u32))),
            Fq2::new(Fq::new(BigUint::from(9u32)), Fq::new(BigUint::from(10u32))),
            Fq2::new(Fq::new(BigUint::from(11u32)), Fq::new(BigUint::from(12u32))),
        );
        Fq12::new(a, b)
    }

    #[test]
    fn mul_one_is_identity() {
        let a = sample();
        assert_eq!(a.mul(&Fq12::one()), a);
    }

    #[test]
    fn square_matches_mul() {
        let a = sample();
        assert_eq!(a.square(), a.mul(&a));
    }

    #[test]
    fn inverse_roundtrip() {
        let a = sample();
        let ai = a.inverse().unwrap();
        assert_eq!(a.mul(&ai), Fq12::one());
    }
}
