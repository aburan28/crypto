//! Self-contained prime-field P-192 group arithmetic for the finite census.
//! Variable time: all points and scalars here are public research inputs.

use num_bigint::BigUint;
use num_traits::Zero;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FieldElement {
    pub value: BigUint,
    pub modulus: BigUint,
}

impl FieldElement {
    pub fn new(value: BigUint, modulus: BigUint) -> Self {
        Self {
            value: value % &modulus,
            modulus,
        }
    }
    pub fn zero(modulus: BigUint) -> Self {
        Self::new(BigUint::zero(), modulus)
    }
    pub fn is_zero(&self) -> bool {
        self.value.is_zero()
    }
    pub fn add(&self, rhs: &Self) -> Self {
        assert_eq!(self.modulus, rhs.modulus);
        Self::new(&self.value + &rhs.value, self.modulus.clone())
    }
    pub fn sub(&self, rhs: &Self) -> Self {
        assert_eq!(self.modulus, rhs.modulus);
        Self::new(
            (&self.value + &self.modulus) - &rhs.value,
            self.modulus.clone(),
        )
    }
    pub fn mul(&self, rhs: &Self) -> Self {
        assert_eq!(self.modulus, rhs.modulus);
        Self::new(&self.value * &rhs.value, self.modulus.clone())
    }
    pub fn pow(&self, exp: &BigUint) -> Self {
        Self::new(self.value.modpow(exp, &self.modulus), self.modulus.clone())
    }
    fn inv(&self) -> Self {
        assert!(!self.is_zero());
        self.pow(&(&self.modulus - BigUint::from(2u32)))
    }
    fn neg(&self) -> Self {
        Self::zero(self.modulus.clone()).sub(self)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Point {
    Infinity,
    Affine { x: FieldElement, y: FieldElement },
}

impl Point {
    pub fn neg(&self) -> Self {
        match self {
            Self::Infinity => Self::Infinity,
            Self::Affine { x, y } => Self::Affine {
                x: x.clone(),
                y: y.neg(),
            },
        }
    }
    pub fn add_vartime(&self, other: &Self, a: &FieldElement) -> Self {
        match (self, other) {
            (Self::Infinity, _) => other.clone(),
            (_, Self::Infinity) => self.clone(),
            (Self::Affine { x: x1, y: y1 }, Self::Affine { x: x2, y: y2 }) => {
                if x1 == x2 {
                    return if y1 == y2 {
                        self.double_vartime(a)
                    } else {
                        Self::Infinity
                    };
                }
                let slope = y2.sub(y1).mul(&x2.sub(x1).inv());
                let x3 = slope.mul(&slope).sub(x1).sub(x2);
                let y3 = slope.mul(&x1.sub(&x3)).sub(y1);
                Self::Affine { x: x3, y: y3 }
            }
        }
    }
    pub fn double_vartime(&self, a: &FieldElement) -> Self {
        match self {
            Self::Infinity => Self::Infinity,
            Self::Affine { x, y } => {
                if y.is_zero() {
                    return Self::Infinity;
                }
                let p = x.modulus.clone();
                let two = FieldElement::new(BigUint::from(2u32), p.clone());
                let three = FieldElement::new(BigUint::from(3u32), p);
                let slope = three.mul(&x.mul(x)).add(a).mul(&two.mul(y).inv());
                let x3 = slope.mul(&slope).sub(&two.mul(x));
                let y3 = slope.mul(&x.sub(&x3)).sub(y);
                Self::Affine { x: x3, y: y3 }
            }
        }
    }
    pub fn scalar_mul(&self, scalar: &BigUint, a: &FieldElement) -> Self {
        let mut acc = Self::Infinity;
        for i in (0..scalar.bits()).rev() {
            acc = acc.double_vartime(a);
            if scalar.bit(i) {
                acc = acc.add_vartime(self, a);
            }
        }
        acc
    }
}

pub struct CurveParams {
    pub p: BigUint,
    pub a: BigUint,
    pub b: BigUint,
    pub n: BigUint,
}

impl CurveParams {
    pub fn p192() -> Self {
        fn hex(s: &str) -> BigUint {
            BigUint::parse_bytes(s.as_bytes(), 16).unwrap()
        }
        Self {
            p: hex("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF"),
            a: hex("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFC"),
            b: hex("64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1"),
            n: hex("FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831"),
        }
    }
    pub fn fe(&self, value: BigUint) -> FieldElement {
        FieldElement::new(value, self.p.clone())
    }
    pub fn generator(&self) -> Point {
        let gx =
            BigUint::parse_bytes(b"188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012", 16).unwrap();
        let gy =
            BigUint::parse_bytes(b"07192B95FFC8DA78631011ED6B24CDD573F977A11E794811", 16).unwrap();
        Point::Affine {
            x: self.fe(gx),
            y: self.fe(gy),
        }
    }
    pub fn a_fe(&self) -> FieldElement {
        self.fe(self.a.clone())
    }
    pub fn is_on_curve(&self, point: &Point) -> bool {
        match point {
            Point::Infinity => true,
            Point::Affine { x, y } => {
                y.mul(y)
                    == x.mul(&x.mul(x))
                        .add(&self.a_fe().mul(x))
                        .add(&self.fe(self.b.clone()))
            }
        }
    }
}

pub mod curve_zoo {
    pub use super::CurveParams;
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn p192_generator_and_order() {
        let curve = CurveParams::p192();
        let g = curve.generator();
        assert!(curve.is_on_curve(&g));
        assert_eq!(g.scalar_mul(&curve.n, &curve.a_fe()), Point::Infinity);
        assert_eq!(g.add_vartime(&g.neg(), &curve.a_fe()), Point::Infinity);
    }
}
