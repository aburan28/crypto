//! **Polynomial arithmetic over the BLS12-381 scalar field `F_r`.**
//!
//! Supports addition, multiplication, division-with-remainder,
//! evaluation, and Lagrange interpolation.  Used as the substrate for
//! KZG commitments and Groth16's QAP construction.

use crate::bls12_381::fq::scalar_modulus;
use num_bigint::BigUint;
use num_traits::{One, Zero};

/// Reduce a scalar mod `r`.
pub fn fr_reduce(v: &BigUint) -> BigUint {
    v % scalar_modulus()
}

pub fn fr_add(a: &BigUint, b: &BigUint) -> BigUint {
    (a + b) % scalar_modulus()
}
pub fn fr_sub(a: &BigUint, b: &BigUint) -> BigUint {
    let r = scalar_modulus();
    ((a + &r) - (b % &r)) % &r
}
pub fn fr_mul(a: &BigUint, b: &BigUint) -> BigUint {
    (a * b) % scalar_modulus()
}
pub fn fr_neg(a: &BigUint) -> BigUint {
    let r = scalar_modulus();
    (&r - (a % &r)) % &r
}
pub fn fr_inv(a: &BigUint) -> Option<BigUint> {
    crate::utils::mod_inverse(a, &scalar_modulus())
}
pub fn fr_pow(base: &BigUint, exp: &BigUint) -> BigUint {
    base.modpow(exp, &scalar_modulus())
}

/// A polynomial over `F_r`, represented in coefficient form (low to
/// high degree).
#[derive(Clone, Debug, PartialEq)]
pub struct Poly {
    pub coeffs: Vec<BigUint>,
}

impl Poly {
    /// Zero polynomial.
    pub fn zero() -> Self {
        Self { coeffs: vec![] }
    }

    /// Polynomial equal to the constant `c`.
    pub fn constant(c: BigUint) -> Self {
        if c.is_zero() {
            Self::zero()
        } else {
            Self { coeffs: vec![fr_reduce(&c)] }
        }
    }

    pub fn from_coeffs(coeffs: Vec<BigUint>) -> Self {
        let mut p = Self { coeffs: coeffs.into_iter().map(|c| fr_reduce(&c)).collect() };
        p.trim();
        p
    }

    /// Strip trailing zero coefficients.
    fn trim(&mut self) {
        while let Some(last) = self.coeffs.last() {
            if last.is_zero() {
                self.coeffs.pop();
            } else {
                break;
            }
        }
    }

    /// Polynomial degree (`-1` for zero polynomial; reported as
    /// `None`).
    pub fn degree(&self) -> Option<usize> {
        if self.coeffs.is_empty() {
            None
        } else {
            Some(self.coeffs.len() - 1)
        }
    }

    /// Evaluate at `x ∈ F_r` via Horner's rule.
    pub fn evaluate(&self, x: &BigUint) -> BigUint {
        let mut acc = BigUint::zero();
        for c in self.coeffs.iter().rev() {
            acc = fr_add(&fr_mul(&acc, x), c);
        }
        acc
    }

    pub fn add(&self, other: &Self) -> Self {
        let n = self.coeffs.len().max(other.coeffs.len());
        let mut out = vec![BigUint::zero(); n];
        for (i, c) in self.coeffs.iter().enumerate() {
            out[i] = fr_add(&out[i], c);
        }
        for (i, c) in other.coeffs.iter().enumerate() {
            out[i] = fr_add(&out[i], c);
        }
        Self::from_coeffs(out)
    }

    pub fn sub(&self, other: &Self) -> Self {
        let n = self.coeffs.len().max(other.coeffs.len());
        let mut out = vec![BigUint::zero(); n];
        for (i, c) in self.coeffs.iter().enumerate() {
            out[i] = fr_add(&out[i], c);
        }
        for (i, c) in other.coeffs.iter().enumerate() {
            out[i] = fr_sub(&out[i], c);
        }
        Self::from_coeffs(out)
    }

    pub fn mul(&self, other: &Self) -> Self {
        if self.coeffs.is_empty() || other.coeffs.is_empty() {
            return Self::zero();
        }
        let mut out = vec![BigUint::zero(); self.coeffs.len() + other.coeffs.len() - 1];
        for (i, a) in self.coeffs.iter().enumerate() {
            if a.is_zero() { continue; }
            for (j, b) in other.coeffs.iter().enumerate() {
                out[i + j] = fr_add(&out[i + j], &fr_mul(a, b));
            }
        }
        Self::from_coeffs(out)
    }

    pub fn scale(&self, c: &BigUint) -> Self {
        Self::from_coeffs(self.coeffs.iter().map(|x| fr_mul(x, c)).collect())
    }

    /// Polynomial division: returns `(q, rem)` with `self = q · divisor + rem`.
    pub fn divmod(&self, divisor: &Self) -> (Self, Self) {
        if divisor.coeffs.is_empty() {
            panic!("division by zero polynomial");
        }
        let d_deg = divisor.degree().unwrap();
        let lead_inv = fr_inv(&divisor.coeffs[d_deg]).expect("leading coef must be invertible");

        let mut rem = self.clone();
        let mut q_coeffs = vec![BigUint::zero(); self.coeffs.len().saturating_sub(d_deg).max(1)];
        while let Some(r_deg) = rem.degree() {
            if r_deg < d_deg {
                break;
            }
            let shift = r_deg - d_deg;
            let q_coef = fr_mul(&rem.coeffs[r_deg], &lead_inv);
            if shift >= q_coeffs.len() {
                q_coeffs.resize(shift + 1, BigUint::zero());
            }
            q_coeffs[shift] = fr_add(&q_coeffs[shift], &q_coef);
            // rem -= q_coef · shift · divisor
            for (i, dc) in divisor.coeffs.iter().enumerate() {
                let term = fr_mul(&q_coef, dc);
                let idx = i + shift;
                while rem.coeffs.len() <= idx {
                    rem.coeffs.push(BigUint::zero());
                }
                rem.coeffs[idx] = fr_sub(&rem.coeffs[idx], &term);
            }
            rem.trim();
        }
        (Self::from_coeffs(q_coeffs), rem)
    }

    /// Lagrange interpolation through `(x_i, y_i)` points.
    pub fn interpolate(points: &[(BigUint, BigUint)]) -> Self {
        let mut result = Self::zero();
        for (i, (xi, yi)) in points.iter().enumerate() {
            // Build basis polynomial L_i(X) = ∏_{j ≠ i} (X − x_j) / (x_i − x_j)
            let mut num = Self::constant(BigUint::one());
            let mut denom = BigUint::one();
            for (j, (xj, _)) in points.iter().enumerate() {
                if i == j { continue; }
                // (X − x_j)
                let factor = Self::from_coeffs(vec![fr_neg(xj), BigUint::one()]);
                num = num.mul(&factor);
                denom = fr_mul(&denom, &fr_sub(xi, xj));
            }
            let denom_inv = fr_inv(&denom).expect("distinct x_i required");
            let term = num.scale(&fr_mul(yi, &denom_inv));
            result = result.add(&term);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evaluation_simple() {
        // p(x) = 1 + 2x + 3x²; p(1) = 6; p(2) = 17.
        let p = Poly::from_coeffs(vec![
            BigUint::from(1u32), BigUint::from(2u32), BigUint::from(3u32),
        ]);
        assert_eq!(p.evaluate(&BigUint::from(1u32)), BigUint::from(6u32));
        assert_eq!(p.evaluate(&BigUint::from(2u32)), BigUint::from(17u32));
    }

    #[test]
    fn add_mul_distributive() {
        let p = Poly::from_coeffs(vec![BigUint::from(1u32), BigUint::from(2u32)]);
        let q = Poly::from_coeffs(vec![BigUint::from(3u32), BigUint::from(4u32)]);
        let r = Poly::from_coeffs(vec![BigUint::from(5u32), BigUint::from(6u32)]);
        let lhs = p.mul(&q.add(&r));
        let rhs = p.mul(&q).add(&p.mul(&r));
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn divmod_roundtrip() {
        // p = x³ + 1; divisor = x − 1; expect q · (x − 1) + r = p.
        let p = Poly::from_coeffs(vec![BigUint::from(1u32), BigUint::zero(), BigUint::zero(), BigUint::from(1u32)]);
        let d = Poly::from_coeffs(vec![fr_neg(&BigUint::one()), BigUint::one()]); // -1 + x
        let (q, rem) = p.divmod(&d);
        assert_eq!(q.mul(&d).add(&rem), p);
    }

    #[test]
    fn interpolate_passes_through_points() {
        // Pick 4 distinct (x, y); interpolation evaluates to y at x.
        let points: Vec<(BigUint, BigUint)> = vec![
            (BigUint::from(1u32), BigUint::from(2u32)),
            (BigUint::from(2u32), BigUint::from(5u32)),
            (BigUint::from(3u32), BigUint::from(10u32)),
            (BigUint::from(4u32), BigUint::from(17u32)),
        ];
        let p = Poly::interpolate(&points);
        for (x, y) in &points {
            assert_eq!(p.evaluate(x), *y);
        }
    }
}
