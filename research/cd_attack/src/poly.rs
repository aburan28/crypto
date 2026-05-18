//! Polynomial arithmetic over F_{p^2}.
//!
//! Coefficients are stored low-to-high (coeffs[i] = coefficient of x^i).
//! Normalized form has no trailing zeros; the zero polynomial is an empty vec.

use crate::field::{F2, Fp2};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly {
    pub coeffs: Vec<F2>,
}

impl Poly {
    /// Build from raw coefficient vector and normalize (strip trailing zeros).
    pub fn new(mut coeffs: Vec<F2>, fp2: &Fp2) -> Self {
        while coeffs.last().map_or(false, |c| fp2.is_zero(c)) {
            coeffs.pop();
        }
        Poly { coeffs }
    }
    pub fn zero() -> Self { Poly { coeffs: vec![] } }
    pub fn one(fp2: &Fp2) -> Self { Poly { coeffs: vec![fp2.one()] } }
    pub fn x(fp2: &Fp2) -> Self {
        Poly { coeffs: vec![fp2.zero(), fp2.one()] }
    }
    pub fn constant(c: F2, fp2: &Fp2) -> Self {
        if fp2.is_zero(&c) { Self::zero() } else { Poly { coeffs: vec![c] } }
    }
    /// Build monomial c · x^n.
    pub fn monomial(c: F2, n: usize, fp2: &Fp2) -> Self {
        if fp2.is_zero(&c) { return Self::zero(); }
        let mut v = vec![fp2.zero(); n];
        v.push(c);
        Poly { coeffs: v }
    }

    pub fn is_zero(&self) -> bool { self.coeffs.is_empty() }
    /// Returns None for the zero polynomial.
    pub fn degree(&self) -> Option<usize> {
        if self.is_zero() { None } else { Some(self.coeffs.len() - 1) }
    }
    pub fn leading(&self, fp2: &Fp2) -> F2 {
        self.coeffs.last().cloned().unwrap_or_else(|| fp2.zero())
    }

    pub fn eval(&self, x: &F2, fp2: &Fp2) -> F2 {
        // Horner
        let mut acc = fp2.zero();
        for c in self.coeffs.iter().rev() {
            acc = fp2.add(&fp2.mul(&acc, x), c);
        }
        acc
    }

    pub fn neg(&self, fp2: &Fp2) -> Self {
        Poly::new(self.coeffs.iter().map(|c| fp2.neg(c)).collect(), fp2)
    }

    pub fn add(&self, other: &Self, fp2: &Fp2) -> Self {
        let n = self.coeffs.len().max(other.coeffs.len());
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let a = self.coeffs.get(i).cloned().unwrap_or_else(|| fp2.zero());
            let b = other.coeffs.get(i).cloned().unwrap_or_else(|| fp2.zero());
            out.push(fp2.add(&a, &b));
        }
        Poly::new(out, fp2)
    }

    pub fn sub(&self, other: &Self, fp2: &Fp2) -> Self {
        let n = self.coeffs.len().max(other.coeffs.len());
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let a = self.coeffs.get(i).cloned().unwrap_or_else(|| fp2.zero());
            let b = other.coeffs.get(i).cloned().unwrap_or_else(|| fp2.zero());
            out.push(fp2.sub(&a, &b));
        }
        Poly::new(out, fp2)
    }

    pub fn mul(&self, other: &Self, fp2: &Fp2) -> Self {
        if self.is_zero() || other.is_zero() { return Self::zero(); }
        let n = self.coeffs.len() + other.coeffs.len() - 1;
        let mut out = vec![fp2.zero(); n];
        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in other.coeffs.iter().enumerate() {
                out[i + j] = fp2.add(&out[i + j], &fp2.mul(a, b));
            }
        }
        Poly::new(out, fp2)
    }

    pub fn scalar_mul(&self, c: &F2, fp2: &Fp2) -> Self {
        if fp2.is_zero(c) { return Self::zero(); }
        Poly::new(self.coeffs.iter().map(|x| fp2.mul(x, c)).collect(), fp2)
    }

    /// Polynomial long division. Returns (quotient, remainder) with
    /// self = q·other + r,  deg(r) < deg(other).  Panics if `other` is zero.
    pub fn div_rem(&self, other: &Self, fp2: &Fp2) -> (Self, Self) {
        assert!(!other.is_zero(), "polynomial division by zero");
        if self.degree().map_or(true, |d| d < other.degree().unwrap()) {
            return (Self::zero(), self.clone());
        }
        let mut r = self.coeffs.clone();
        let d_o = other.coeffs.len() - 1;
        let lc_inv = fp2.inv(&other.leading(fp2));
        let mut q = vec![fp2.zero(); self.coeffs.len() - d_o];
        while r.len() > d_o {
            let lc_r = r.last().unwrap().clone();
            if fp2.is_zero(&lc_r) {
                r.pop();
                continue;
            }
            let scale = fp2.mul(&lc_r, &lc_inv);
            let pos = r.len() - d_o - 1;
            q[pos] = scale.clone();
            // r -= scale · x^pos · other
            for (i, b) in other.coeffs.iter().enumerate() {
                let idx = pos + i;
                r[idx] = fp2.sub(&r[idx], &fp2.mul(&scale, b));
            }
            // Top coefficient should now be zero — pop it.
            r.pop();
        }
        (Poly::new(q, fp2), Poly::new(r, fp2))
    }

    pub fn rem(&self, other: &Self, fp2: &Fp2) -> Self {
        self.div_rem(other, fp2).1
    }

    pub fn make_monic(&self, fp2: &Fp2) -> Self {
        if self.is_zero() { return Self::zero(); }
        let inv = fp2.inv(&self.leading(fp2));
        self.scalar_mul(&inv, fp2)
    }

    /// Formal derivative: d/dx of the polynomial.
    pub fn derivative(&self, fp2: &Fp2) -> Self {
        let coeffs: Vec<F2> = self
            .coeffs
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| fp2.mul_scalar(c, i as i64))
            .collect();
        Poly::new(coeffs, fp2)
    }
}

// Greatest common divisor (Euclidean). Returned as a monic polynomial.
pub fn gcd(a: &Poly, b: &Poly, fp2: &Fp2) -> Poly {
    let mut x = a.clone();
    let mut y = b.clone();
    while !y.is_zero() {
        let r = x.rem(&y, fp2);
        x = y;
        y = r;
    }
    x.make_monic(fp2)
}

/// Extended Euclidean: returns (g, u, v) with u·a + v·b = g where g is monic.
pub fn ext_gcd(a: &Poly, b: &Poly, fp2: &Fp2) -> (Poly, Poly, Poly) {
    let (mut r0, mut r1) = (a.clone(), b.clone());
    let (mut s0, mut s1) = (Poly::one(fp2), Poly::zero());
    let (mut t0, mut t1) = (Poly::zero(), Poly::one(fp2));
    while !r1.is_zero() {
        let (q, r) = r0.div_rem(&r1, fp2);
        let s_new = s0.sub(&q.mul(&s1, fp2), fp2);
        let t_new = t0.sub(&q.mul(&t1, fp2), fp2);
        r0 = r1; r1 = r;
        s0 = s1; s1 = s_new;
        t0 = t1; t1 = t_new;
    }
    // Normalize to monic g; scale s, t accordingly.
    if r0.is_zero() { return (Poly::zero(), s0, t0); }
    let lc_inv = fp2.inv(&r0.leading(fp2));
    let g = r0.scalar_mul(&lc_inv, fp2);
    let u = s0.scalar_mul(&lc_inv, fp2);
    let v = t0.scalar_mul(&lc_inv, fp2);
    (g, u, v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{Fp, Fp2};
    use num_bigint::BigInt;

    fn ctx() -> Fp2 { Fp2::new(Fp::new(BigInt::from(431u64))) }

    #[test]
    fn poly_basics() {
        let fp2 = ctx();
        let x = Poly::x(&fp2);
        let one = Poly::one(&fp2);
        // (x + 1)^2 = x^2 + 2x + 1
        let xp1 = x.add(&one, &fp2);
        let sq = xp1.mul(&xp1, &fp2);
        let expected = Poly::new(
            vec![fp2.from_int(1), fp2.from_int(2), fp2.from_int(1)], &fp2);
        assert_eq!(sq, expected);
    }

    #[test]
    fn poly_div() {
        let fp2 = ctx();
        // (x^2 + 2x + 1) / (x + 1) = (x + 1) remainder 0
        let num = Poly::new(
            vec![fp2.from_int(1), fp2.from_int(2), fp2.from_int(1)], &fp2);
        let den = Poly::new(vec![fp2.from_int(1), fp2.from_int(1)], &fp2);
        let (q, r) = num.div_rem(&den, &fp2);
        assert_eq!(q, den);
        assert!(r.is_zero());
    }

    #[test]
    fn poly_div_rem_random() {
        let fp2 = ctx();
        // num = (x + 2)(x + 3) + (5) = x^2 + 5x + 11
        let den = Poly::new(vec![fp2.from_int(2), fp2.from_int(1)], &fp2);
        let q_expected = Poly::new(vec![fp2.from_int(3), fp2.from_int(1)], &fp2);
        let r_expected = Poly::constant(fp2.from_int(5), &fp2);
        let num = den.mul(&q_expected, &fp2).add(&r_expected, &fp2);
        let (q, r) = num.div_rem(&den, &fp2);
        assert_eq!(q, q_expected);
        assert_eq!(r, r_expected);
    }

    #[test]
    fn poly_ext_gcd() {
        let fp2 = ctx();
        // a = (x+1)(x+2)(x+3),  b = (x+2)(x+5).  gcd should be x+2 (monic).
        let f = |c: i64| Poly::new(vec![fp2.from_int(c), fp2.one()], &fp2);
        let a = f(1).mul(&f(2), &fp2).mul(&f(3), &fp2);
        let b = f(2).mul(&f(5), &fp2);
        let (g, u, v) = ext_gcd(&a, &b, &fp2);
        assert_eq!(g, f(2));
        // Bezout: u·a + v·b == g
        let lhs = u.mul(&a, &fp2).add(&v.mul(&b, &fp2), &fp2);
        assert_eq!(lhs, g);
    }
}
