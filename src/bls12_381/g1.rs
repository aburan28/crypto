//! BLS12-381 G1: curve `E(F_p): y² = x³ + 4`.
//!
//! Affine representation with explicit identity (point at infinity).

use super::fq::{scalar_modulus, Fq};
use num_bigint::BigUint;
use num_traits::{One, Zero};

#[derive(Clone, Debug, PartialEq)]
pub enum G1Point {
    Infinity,
    Affine { x: Fq, y: Fq },
}

impl G1Point {
    /// Standard generator (from Zcash protocol spec §5.4.9.1).
    pub fn generator() -> Self {
        let x = Fq::new(BigUint::parse_bytes(
            b"17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb",
            16,
        ).unwrap());
        let y = Fq::new(BigUint::parse_bytes(
            b"08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1",
            16,
        ).unwrap());
        Self::Affine { x, y }
    }

    pub fn is_infinity(&self) -> bool {
        matches!(self, Self::Infinity)
    }

    /// Verify `y² = x³ + 4` (and infinity is trivially on curve).
    pub fn is_on_curve(&self) -> bool {
        match self {
            Self::Infinity => true,
            Self::Affine { x, y } => {
                let lhs = y.square();
                let rhs = x.square().mul(x).add(&Fq::new(BigUint::from(4u32)));
                lhs == rhs
            }
        }
    }

    /// Negation: `(x, y) ↦ (x, −y)`.
    pub fn neg(&self) -> Self {
        match self {
            Self::Infinity => Self::Infinity,
            Self::Affine { x, y } => Self::Affine { x: x.clone(), y: y.neg() },
        }
    }

    /// Point addition.  Affine; handles all special cases.
    pub fn add(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Infinity, p) | (p, Self::Infinity) => p.clone(),
            (Self::Affine { x: x1, y: y1 }, Self::Affine { x: x2, y: y2 }) => {
                if x1 == x2 {
                    if y1 == &y2.neg() || y1.is_zero() {
                        return Self::Infinity;
                    }
                    return self.double();
                }
                let lambda = y2.sub(y1).mul(&x2.sub(x1).inverse().unwrap());
                let x3 = lambda.square().sub(x1).sub(x2);
                let y3 = lambda.mul(&x1.sub(&x3)).sub(y1);
                Self::Affine { x: x3, y: y3 }
            }
        }
    }

    pub fn double(&self) -> Self {
        match self {
            Self::Infinity => Self::Infinity,
            Self::Affine { x, y } => {
                if y.is_zero() {
                    return Self::Infinity;
                }
                let three = Fq::new(BigUint::from(3u32));
                let two = Fq::new(BigUint::from(2u32));
                let lambda = three.mul(&x.square()).mul(&two.mul(y).inverse().unwrap());
                let x3 = lambda.square().sub(x).sub(x);
                let y3 = lambda.mul(&x.sub(&x3)).sub(y);
                Self::Affine { x: x3, y: y3 }
            }
        }
    }

    /// `[k] · self` via double-and-add (variable-time).
    pub fn scalar_mul(&self, k: &BigUint) -> Self {
        let mut result = Self::Infinity;
        let mut addend = self.clone();
        let mut e = k.clone();
        while !e.is_zero() {
            if (&e & BigUint::one()) == BigUint::one() {
                result = result.add(&addend);
            }
            addend = addend.double();
            e >>= 1;
        }
        result
    }

    /// Subgroup check: `[r] · self == O` where `r` is the group order.
    /// For G1 of BLS12-381, the cofactor is 1 in the primary subgroup
    /// — but the curve has full group of order `r · h₁` with
    /// `h₁ ≠ 1`, so this matters in some embeddings.  The well-known
    /// G1 generator above has cofactor cleared.
    pub fn is_in_correct_subgroup(&self) -> bool {
        let r = scalar_modulus();
        self.scalar_mul(&r).is_infinity()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generator_is_on_curve() {
        let g = G1Point::generator();
        assert!(g.is_on_curve());
    }

    #[test]
    fn double_matches_self_add() {
        let g = G1Point::generator();
        let dbl = g.double();
        let added = g.add(&g);
        assert_eq!(dbl, added);
    }

    #[test]
    fn add_neg_is_infinity() {
        let g = G1Point::generator();
        let neg_g = g.neg();
        assert_eq!(g.add(&neg_g), G1Point::Infinity);
    }

    #[test]
    fn scalar_mul_1_is_self() {
        let g = G1Point::generator();
        assert_eq!(g.scalar_mul(&BigUint::one()), g);
    }

    #[test]
    fn scalar_mul_associates() {
        let g = G1Point::generator();
        let two_g = g.scalar_mul(&BigUint::from(2u32));
        let three_g_a = g.scalar_mul(&BigUint::from(3u32));
        let three_g_b = two_g.add(&g);
        assert_eq!(three_g_a, three_g_b);
    }

    #[test]
    fn subgroup_check_on_generator() {
        let g = G1Point::generator();
        assert!(g.is_in_correct_subgroup(), "[r]·G == O for the BLS12-381 generator");
    }
}
