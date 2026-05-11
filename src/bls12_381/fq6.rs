//! Cubic extension `F_{p⁶} = F_{p²}[v]/(v³ − (u + 1))`.
//!
//! Element: `c₀ + c₁·v + c₂·v²` with `c_i ∈ F_{p²}`.
//! Multiplication uses `v³ = u + 1`.

use super::fq2::Fq2;

#[derive(Clone, Debug, PartialEq)]
pub struct Fq6 {
    pub c0: Fq2,
    pub c1: Fq2,
    pub c2: Fq2,
}

impl Fq6 {
    pub fn new(c0: Fq2, c1: Fq2, c2: Fq2) -> Self {
        Self { c0, c1, c2 }
    }
    pub fn zero() -> Self {
        Self { c0: Fq2::zero(), c1: Fq2::zero(), c2: Fq2::zero() }
    }
    pub fn one() -> Self {
        Self { c0: Fq2::one(), c1: Fq2::zero(), c2: Fq2::zero() }
    }
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }
    pub fn add(&self, other: &Self) -> Self {
        Self {
            c0: self.c0.add(&other.c0),
            c1: self.c1.add(&other.c1),
            c2: self.c2.add(&other.c2),
        }
    }
    pub fn sub(&self, other: &Self) -> Self {
        Self {
            c0: self.c0.sub(&other.c0),
            c1: self.c1.sub(&other.c1),
            c2: self.c2.sub(&other.c2),
        }
    }
    pub fn neg(&self) -> Self {
        Self { c0: self.c0.neg(), c1: self.c1.neg(), c2: self.c2.neg() }
    }
    /// Schoolbook multiplication using `v³ = ξ = u + 1`.
    /// (a + b·v + c·v²)(d + e·v + f·v²) =
    ///   (ad + ξ·(bf + ce)) + (ae + bd + ξ·cf)·v + (af + be + cd)·v²
    pub fn mul(&self, other: &Self) -> Self {
        let ad = self.c0.mul(&other.c0);
        let ae = self.c0.mul(&other.c1);
        let af = self.c0.mul(&other.c2);
        let bd = self.c1.mul(&other.c0);
        let be = self.c1.mul(&other.c1);
        let bf = self.c1.mul(&other.c2);
        let cd = self.c2.mul(&other.c0);
        let ce = self.c2.mul(&other.c1);
        let cf = self.c2.mul(&other.c2);

        let c0 = ad.add(&bf.add(&ce).mul_by_nonresidue());
        let c1 = ae.add(&bd).add(&cf.mul_by_nonresidue());
        let c2 = af.add(&be).add(&cd);
        Self { c0, c1, c2 }
    }
    pub fn square(&self) -> Self {
        self.mul(self)
    }
    /// Multiply by `v` (the cubic non-residue's variable).
    pub fn mul_by_v(&self) -> Self {
        // (a + b·v + c·v²) · v = a·v + b·v² + c·v³ = c·ξ + a·v + b·v²
        Self {
            c0: self.c2.mul_by_nonresidue(),
            c1: self.c0.clone(),
            c2: self.c1.clone(),
        }
    }
    pub fn inverse(&self) -> Option<Self> {
        // Use the formula for inversion in cubic extensions:
        // t0 = c0² − ξ·c1·c2
        // t1 = ξ·c2² − c0·c1
        // t2 = c1² − c0·c2
        // norm = c0·t0 + ξ·c2·t1 + ξ·c1·t2
        // inv = (t0, t1, t2) / norm
        let t0 = self.c0.square().sub(&self.c1.mul(&self.c2).mul_by_nonresidue());
        let t1 = self.c2.square().mul_by_nonresidue().sub(&self.c0.mul(&self.c1));
        let t2 = self.c1.square().sub(&self.c0.mul(&self.c2));
        let norm = self.c0.mul(&t0)
            .add(&self.c2.mul(&t1).mul_by_nonresidue())
            .add(&self.c1.mul(&t2).mul_by_nonresidue());
        let norm_inv = norm.inverse()?;
        Some(Self {
            c0: t0.mul(&norm_inv),
            c1: t1.mul(&norm_inv),
            c2: t2.mul(&norm_inv),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::fq::Fq;
    use num_bigint::BigUint;

    fn sample() -> Fq6 {
        Fq6::new(
            Fq2::new(Fq::new(BigUint::from(1u32)), Fq::new(BigUint::from(2u32))),
            Fq2::new(Fq::new(BigUint::from(3u32)), Fq::new(BigUint::from(4u32))),
            Fq2::new(Fq::new(BigUint::from(5u32)), Fq::new(BigUint::from(6u32))),
        )
    }

    #[test]
    fn add_neg_inverse() {
        let a = sample();
        let neg_a = a.neg();
        assert_eq!(a.add(&neg_a), Fq6::zero());
    }

    #[test]
    fn mul_one_is_identity() {
        let a = sample();
        assert_eq!(a.mul(&Fq6::one()), a);
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
        assert_eq!(a.mul(&ai), Fq6::one());
    }

    #[test]
    fn mul_distributive() {
        let a = sample();
        let b = Fq6::new(
            Fq2::new(Fq::new(BigUint::from(7u32)), Fq::new(BigUint::from(8u32))),
            Fq2::new(Fq::new(BigUint::from(9u32)), Fq::new(BigUint::from(10u32))),
            Fq2::new(Fq::new(BigUint::from(11u32)), Fq::new(BigUint::from(12u32))),
        );
        let c = Fq6::new(
            Fq2::new(Fq::new(BigUint::from(13u32)), Fq::new(BigUint::from(14u32))),
            Fq2::new(Fq::new(BigUint::from(15u32)), Fq::new(BigUint::from(16u32))),
            Fq2::new(Fq::new(BigUint::from(17u32)), Fq::new(BigUint::from(18u32))),
        );
        let lhs = a.mul(&b.add(&c));
        let rhs = a.mul(&b).add(&a.mul(&c));
        assert_eq!(lhs, rhs);
    }
}
