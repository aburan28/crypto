//! BLS12-381 G2: curve `E'(F_{p²}): y² = x³ + 4(u + 1)` over `F_{p²}`.
//!
//! This is the **sextic twist** of G1.  The pairing `e(P, Q)` takes
//! `P ∈ G1`, `Q ∈ G2`.

use super::fq::{scalar_modulus, Fq};
use super::fq2::Fq2;
use num_bigint::BigUint;
use num_traits::{One, Zero};

#[derive(Clone, Debug, PartialEq)]
pub enum G2Point {
    Infinity,
    Affine { x: Fq2, y: Fq2 },
}

impl G2Point {
    /// Standard G2 generator (Zcash protocol spec §5.4.9.2).
    pub fn generator() -> Self {
        // x = (x0 + x1·u), y = (y0 + y1·u)
        let x0 = Fq::new(BigUint::parse_bytes(
            b"024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8",
            16,
        ).unwrap());
        let x1 = Fq::new(BigUint::parse_bytes(
            b"13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e",
            16,
        ).unwrap());
        let y0 = Fq::new(BigUint::parse_bytes(
            b"0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801",
            16,
        ).unwrap());
        let y1 = Fq::new(BigUint::parse_bytes(
            b"0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be",
            16,
        ).unwrap());
        Self::Affine {
            x: Fq2::new(x0, x1),
            y: Fq2::new(y0, y1),
        }
    }

    pub fn is_infinity(&self) -> bool {
        matches!(self, Self::Infinity)
    }

    /// Verify `y² = x³ + 4(u + 1)`.
    pub fn is_on_curve(&self) -> bool {
        match self {
            Self::Infinity => true,
            Self::Affine { x, y } => {
                let lhs = y.square();
                let four_u_plus_1 = Fq2::new(Fq::new(BigUint::from(4u32)), Fq::new(BigUint::from(4u32)));
                let rhs = x.square().mul(x).add(&four_u_plus_1);
                lhs == rhs
            }
        }
    }

    pub fn neg(&self) -> Self {
        match self {
            Self::Infinity => Self::Infinity,
            Self::Affine { x, y } => Self::Affine { x: x.clone(), y: y.neg() },
        }
    }

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
                let three = Fq2::new(Fq::new(BigUint::from(3u32)), Fq::zero());
                let two = Fq2::new(Fq::new(BigUint::from(2u32)), Fq::zero());
                let lambda = three.mul(&x.square()).mul(&two.mul(y).inverse().unwrap());
                let x3 = lambda.square().sub(x).sub(x);
                let y3 = lambda.mul(&x.sub(&x3)).sub(y);
                Self::Affine { x: x3, y: y3 }
            }
        }
    }

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

    /// `[r] · self == O` test for primary-subgroup membership.
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
        let g = G2Point::generator();
        assert!(g.is_on_curve());
    }

    #[test]
    fn double_matches_self_add() {
        let g = G2Point::generator();
        assert_eq!(g.double(), g.add(&g));
    }

    #[test]
    fn add_neg_is_infinity() {
        let g = G2Point::generator();
        assert_eq!(g.add(&g.neg()), G2Point::Infinity);
    }

    #[test]
    fn scalar_mul_associates() {
        let g = G2Point::generator();
        let two_g = g.scalar_mul(&BigUint::from(2u32));
        let three_a = g.scalar_mul(&BigUint::from(3u32));
        let three_b = two_g.add(&g);
        assert_eq!(three_a, three_b);
    }

    #[test]
    #[ignore] // Slow; ~ minutes due to BigUint multiplication on Fq2.
    fn subgroup_check_on_generator() {
        let g = G2Point::generator();
        assert!(g.is_in_correct_subgroup());
    }
}
