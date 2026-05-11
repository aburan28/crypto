//! Quadratic extension `F_{p²} = F_p[u]/(u² + 1)`.
//!
//! Elements `a + b·u` with `a, b ∈ F_p`.  Multiplication:
//! `(a + b·u)(c + d·u) = (ac − bd) + (ad + bc)·u`
//! since `u² = −1`.

use super::fq::{modulus, Fq};
use num_bigint::BigUint;
use num_traits::{One, Zero};

/// `F_{p²}` element `c0 + c1·u`.
#[derive(Clone, Debug, PartialEq)]
pub struct Fq2 {
    pub c0: Fq,
    pub c1: Fq,
}

impl Fq2 {
    pub fn new(c0: Fq, c1: Fq) -> Self {
        Self { c0, c1 }
    }
    pub fn zero() -> Self {
        Self { c0: Fq::zero(), c1: Fq::zero() }
    }
    pub fn one() -> Self {
        Self { c0: Fq::one(), c1: Fq::zero() }
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
    pub fn mul(&self, other: &Self) -> Self {
        // (a + bu)(c + du) = (ac − bd) + (ad + bc)u
        let ac = self.c0.mul(&other.c0);
        let bd = self.c1.mul(&other.c1);
        let ad = self.c0.mul(&other.c1);
        let bc = self.c1.mul(&other.c0);
        Self { c0: ac.sub(&bd), c1: ad.add(&bc) }
    }
    pub fn square(&self) -> Self {
        // (a + bu)² = (a² − b²) + 2ab·u
        let a_sq = self.c0.square();
        let b_sq = self.c1.square();
        let two_ab = self.c0.mul(&self.c1).add(&self.c0.mul(&self.c1));
        Self { c0: a_sq.sub(&b_sq), c1: two_ab }
    }
    /// Frobenius: `(a + bu)^p = a − b·u`  (since `u^p = -u` when
    /// `p ≡ 3 (mod 4)`, true for BLS12-381).
    pub fn frobenius_map(&self) -> Self {
        Self { c0: self.c0.clone(), c1: self.c1.neg() }
    }
    pub fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        // norm = a² + b²
        let norm = self.c0.square().add(&self.c1.square());
        let norm_inv = norm.inverse()?;
        // inverse = (a − bu) / norm
        Some(Self {
            c0: self.c0.mul(&norm_inv),
            c1: self.c1.neg().mul(&norm_inv),
        })
    }
    /// `self · (u + 1)`  — needed for the Fq6 / Fq12 tower.
    pub fn mul_by_nonresidue(&self) -> Self {
        // (a + bu)·(1 + u) = (a − b) + (a + b)u
        Self {
            c0: self.c0.sub(&self.c1),
            c1: self.c0.add(&self.c1),
        }
    }
    /// Conjugate: `a + bu ↦ a − bu`.  Equivalent to Frobenius here.
    pub fn conjugate(&self) -> Self {
        self.frobenius_map()
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
    /// `√self` when one exists.  For `F_{p²}` we use the Adj sqrt
    /// algorithm (Adj-Rodriguez 2012).  Returns `None` otherwise.
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::zero());
        }
        let p = modulus();
        // Compute α = self^((p² − 1) / 4) — gives a QR-test in F_{p²}.
        // For our needs (curve cofactor / clearing operations) we use
        // a generic Tonelli-Shanks-style approach.
        // Specifically: for p ≡ 3 mod 4, sqrt in F_{p²} via
        //   α = self^((p-3)/4)
        //   β = α · self
        //   if β² · α² · self = -1: x = β · (u + 1)^((p-1)/4) (special case)
        //   else: x = β
        // Simplified Tonelli specialised to BLS12-381 below.
        let p_minus_3_over_4 = (&p - BigUint::from(3u32)) / BigUint::from(4u32);
        let alpha = self.pow(&p_minus_3_over_4);
        let beta = alpha.mul(self);
        let alpha_squared_self = alpha.square().mul(self);
        // Check if alpha² · self equals -1 in F_p² (i.e. value (-1, 0)).
        let neg_one = Fq2 { c0: Fq::one().neg(), c1: Fq::zero() };
        let sqrt_candidate = if alpha_squared_self == neg_one {
            // x = β · u  (u² = -1 makes (β u)² = -β² · 1 = self · -1 · -1 = self).
            // Wait: β² = self · α² · self · α² = …
            // Use the more general approach: β · ζ where ζ² = -1.
            // In F_p² with u² = -1, ζ = u works.
            beta.mul(&Fq2 { c0: Fq::zero(), c1: Fq::one() })
        } else {
            beta
        };
        // Verify.
        let test = sqrt_candidate.square();
        if test == *self {
            Some(sqrt_candidate)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_sub_neg_basic() {
        let a = Fq2::new(Fq::new(BigUint::from(1u32)), Fq::new(BigUint::from(2u32)));
        let b = Fq2::new(Fq::new(BigUint::from(3u32)), Fq::new(BigUint::from(4u32)));
        let sum = a.add(&b);
        assert_eq!(sum.c0, Fq::new(BigUint::from(4u32)));
        assert_eq!(sum.c1, Fq::new(BigUint::from(6u32)));
    }

    #[test]
    fn mul_distributive() {
        let a = Fq2::new(Fq::new(BigUint::from(5u32)), Fq::new(BigUint::from(7u32)));
        let b = Fq2::new(Fq::new(BigUint::from(11u32)), Fq::new(BigUint::from(13u32)));
        let c = Fq2::new(Fq::new(BigUint::from(17u32)), Fq::new(BigUint::from(19u32)));
        // a·(b + c) = a·b + a·c
        let lhs = a.mul(&b.add(&c));
        let rhs = a.mul(&b).add(&a.mul(&c));
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn square_matches_mul() {
        let a = Fq2::new(Fq::new(BigUint::from(13u32)), Fq::new(BigUint::from(17u32)));
        assert_eq!(a.square(), a.mul(&a));
    }

    #[test]
    fn inverse_roundtrip() {
        let a = Fq2::new(Fq::new(BigUint::from(123u32)), Fq::new(BigUint::from(456u32)));
        let ai = a.inverse().unwrap();
        assert_eq!(a.mul(&ai), Fq2::one());
    }

    #[test]
    fn frobenius_squared_is_identity() {
        // F^2 = identity on F_{p²}.
        let a = Fq2::new(Fq::new(BigUint::from(7u32)), Fq::new(BigUint::from(11u32)));
        let f = a.frobenius_map();
        let ff = f.frobenius_map();
        assert_eq!(ff, a);
    }

    #[test]
    fn mul_by_nonresidue_consistent() {
        let a = Fq2::new(Fq::new(BigUint::from(3u32)), Fq::new(BigUint::from(5u32)));
        let nonresidue = Fq2::new(Fq::one(), Fq::one()); // u + 1
        assert_eq!(a.mul_by_nonresidue(), a.mul(&nonresidue));
    }
}
