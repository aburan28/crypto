//! Genus-2 hyperelliptic curve Jacobian arithmetic via Cantor's algorithm.
//!
//! Curve model:  C: y² = f(x)  with deg f = 5  (imaginary hyperelliptic).
//! Mumford representation of a divisor class:  D = ⟨u(x), v(x)⟩  where
//!   - u monic, deg u ≤ 2
//!   - deg v < deg u
//!   - v² ≡ f  (mod u)     (this is what makes the divisor "principal-modulo")
//! Identity: ⟨1, 0⟩.    Negation: ⟨u, v⟩ → ⟨u, −v mod u⟩.
//!
//! References:
//!   Cantor, "Computing in the Jacobian of a hyperelliptic curve" (1987).
//!   Galbraith, "Mathematics of Public Key Cryptography" Ch. 10.

use crate::field::{F2, Fp2};
use crate::poly::{ext_gcd, gcd, Poly};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::Zero;

/// Hyperelliptic curve  y² = f(x),  deg f = 5.
#[derive(Clone)]
pub struct Curve {
    pub f: Poly,
}

impl Curve {
    /// Generic constructor — accepts any non-zero f. For Cantor arithmetic on
    /// this Curve, the caller should pass monic f (the reduction loop assumes
    /// it). Non-monic f is accepted to allow representing Richelot codomains
    /// which may emerge with non-square leading coefficient, where monicizing
    /// would produce a twist instead of the genuine isogenous curve.
    pub fn new(f: Poly, _fp2: &Fp2) -> Self {
        assert!(f.degree().is_some(), "f cannot be zero");
        Curve { f }
    }
    /// Strict constructor for curves used in Cantor arithmetic: insists on
    /// deg f = 5 and monic.
    pub fn new_monic_deg5(f: Poly, fp2: &Fp2) -> Self {
        assert_eq!(f.degree(), Some(5), "expected deg f = 5");
        assert_eq!(f.leading(fp2), fp2.one(), "f must be monic");
        Curve { f }
    }
}

/// Mumford representation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Div {
    pub u: Poly,
    pub v: Poly,
}

impl Div {
    pub fn identity(fp2: &Fp2) -> Self {
        Div { u: Poly::one(fp2), v: Poly::zero() }
    }
    pub fn is_identity(&self, fp2: &Fp2) -> bool {
        self.u == Poly::one(fp2) && self.v.is_zero()
    }
    pub fn neg(&self, fp2: &Fp2) -> Self {
        let neg_v = self.v.neg(fp2);
        Div { u: self.u.clone(), v: if self.u.is_zero() { neg_v } else { neg_v.rem(&self.u, fp2) } }
    }
    /// Sanity-check the Mumford invariants.
    pub fn is_valid(&self, c: &Curve, fp2: &Fp2) -> bool {
        if self.u.is_zero() { return false; }
        if self.u.leading(fp2) != fp2.one() { return false; }
        let du = self.u.degree().unwrap();
        if du > 2 { return false; }
        if !self.v.is_zero() && self.v.degree().unwrap() >= du { return false; }
        // v² ≡ f (mod u)
        let v2 = self.v.mul(&self.v, fp2);
        let lhs = v2.rem(&self.u, fp2);
        let rhs = c.f.rem(&self.u, fp2);
        lhs == rhs
    }
}

/// Cantor's algorithm: composition then reduction. Returns D_1 + D_2 in Jac(C).
pub fn add(c: &Curve, d1: &Div, d2: &Div, fp2: &Fp2) -> Div {
    if d1.is_identity(fp2) { return d2.clone(); }
    if d2.is_identity(fp2) { return d1.clone(); }

    let (u1, v1) = (&d1.u, &d1.v);
    let (u2, v2) = (&d2.u, &d2.v);
    let v_sum = v1.add(v2, fp2);

    // Step 1: d1 = gcd(u1, u2) = e1·u1 + e2·u2
    let (d1g, e1, e2) = ext_gcd(u1, u2, fp2);

    // Step 2: d = gcd(d1g, v_sum) = c1·d1g + c2·v_sum
    let (d, c1, c2) = ext_gcd(&d1g, &v_sum, fp2);

    // Step 3: s1 = c1·e1,  s2 = c1·e2,  s3 = c2
    let s1 = c1.mul(&e1, fp2);
    let s2 = c1.mul(&e2, fp2);
    let s3 = c2;

    // Step 4: u = u1·u2 / d²
    let u_num = u1.mul(u2, fp2);
    let d2_sq = d.mul(&d, fp2);
    let (mut u, rem) = u_num.div_rem(&d2_sq, fp2);
    debug_assert!(rem.is_zero(), "Cantor composition: d² should divide u1·u2");
    u = u.make_monic(fp2);

    // Step 5: v = (s1·u1·v2 + s2·u2·v1 + s3·(v1·v2 + f)) / d  (mod u)
    let term1 = s1.mul(u1, fp2).mul(v2, fp2);
    let term2 = s2.mul(u2, fp2).mul(v1, fp2);
    let v1v2_plus_f = v1.mul(v2, fp2).add(&c.f, fp2);
    let term3 = s3.mul(&v1v2_plus_f, fp2);
    let v_num = term1.add(&term2, fp2).add(&term3, fp2);
    let (v_div_d, vr) = v_num.div_rem(&d, fp2);
    debug_assert!(vr.is_zero(), "Cantor composition: d should divide v_num");
    let v = if u.is_zero() { v_div_d } else { v_div_d.rem(&u, fp2) };

    reduce(c, &u, &v, fp2)
}

/// Reduction loop: bring (u, v) into Mumford form with deg u ≤ 2.
pub fn reduce(c: &Curve, u: &Poly, v: &Poly, fp2: &Fp2) -> Div {
    let mut u = u.clone();
    let mut v = v.clone();
    while u.degree().map_or(false, |d| d > 2) {
        // u' = (f − v²) / u
        let v2 = v.mul(&v, fp2);
        let num = c.f.sub(&v2, fp2);
        let (u_new, rem) = num.div_rem(&u, fp2);
        debug_assert!(rem.is_zero(), "reduce: u must divide (f − v²)");
        // v' = −v mod u'
        let v_new = v.neg(fp2).rem(&u_new, fp2);
        u = u_new.make_monic(fp2);
        v = v_new;
    }
    // Always normalize u to monic for canonical Mumford form. v is taken mod u
    // so this final reduction keeps the divisor unchanged.
    if !u.is_zero() {
        u = u.make_monic(fp2);
        v = v.rem(&u, fp2);
    }
    Div { u, v }
}

/// Scalar multiplication via left-to-right double-and-add.
pub fn scalar_mul(c: &Curve, k: &BigInt, d: &Div, fp2: &Fp2) -> Div {
    if k.is_zero() || d.is_identity(fp2) { return Div::identity(fp2); }
    let one = BigInt::from(1);
    let mut r = Div::identity(fp2);
    let bits = k.bits();
    for i in (0..bits).rev() {
        r = add(c, &r, &r, fp2);
        if ((k >> i) & &one) == one {
            r = add(c, &r, d, fp2);
        }
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{Fp, Fp2};

    fn ctx() -> Fp2 { Fp2::new(Fp::new(BigInt::from(431u64))) }

    /// A specific genus-2 curve over F_p ⊂ F_{p²}, f(x) = x⁵ + x + 1.
    fn test_curve(fp2: &Fp2) -> Curve {
        let f = Poly::new(vec![
            fp2.from_int(1), // 1
            fp2.from_int(1), // x
            fp2.from_int(0), // x²
            fp2.from_int(0), // x³
            fp2.from_int(0), // x⁴
            fp2.from_int(1), // x⁵
        ], fp2);
        Curve::new(f, fp2)
    }

    /// Build a divisor from a single affine point (x₀, y₀) on C.
    fn divisor_from_point(x0: &F2, y0: &F2, fp2: &Fp2) -> Div {
        // u = (x − x₀),  v = y₀
        let u = Poly::new(vec![fp2.neg(x0), fp2.one()], fp2);
        let v = Poly::constant(y0.clone(), fp2);
        Div { u, v }
    }

    /// Sample a random point on the curve over F_p (real part only) for tests.
    fn sample_point(c: &Curve, x_try: i64, fp2: &Fp2) -> Option<(F2, F2)> {
        let x = fp2.from_int(x_try);
        let rhs = c.f.eval(&x, fp2);
        if rhs.b != BigInt::zero() { return None; }
        let y_real = fp2.fp.sqrt(&rhs.a)?;
        Some((x, F2 { a: y_real, b: BigInt::zero() }))
    }

    #[test]
    fn identity_neg() {
        let fp2 = ctx();
        let c = test_curve(&fp2);
        let zero = Div::identity(&fp2);
        // Find a divisor: combine two points.
        let mut pts = vec![];
        for k in 0..50i64 {
            if let Some(p) = sample_point(&c, k, &fp2) {
                pts.push(p);
                if pts.len() >= 2 { break; }
            }
        }
        assert!(pts.len() >= 2, "couldn't find 2 points on test curve");
        let d1 = divisor_from_point(&pts[0].0, &pts[0].1, &fp2);
        assert!(d1.is_valid(&c, &fp2));
        // D + 0 = D
        assert_eq!(add(&c, &d1, &zero, &fp2), d1);
        // 0 + D = D
        assert_eq!(add(&c, &zero, &d1, &fp2), d1);
        // D + (−D) = 0
        let d1_neg = d1.neg(&fp2);
        let sum = add(&c, &d1, &d1_neg, &fp2);
        assert!(sum.is_identity(&fp2),
                "D + (−D) should be 0, got u={:?} v={:?}", sum.u, sum.v);
    }

    #[test]
    fn commutative_associative() {
        let fp2 = ctx();
        let c = test_curve(&fp2);
        let mut pts = vec![];
        for k in 0..200i64 {
            if let Some(p) = sample_point(&c, k, &fp2) {
                pts.push(p);
                if pts.len() >= 6 { break; }
            }
        }
        assert!(pts.len() >= 6, "need 6 points on test curve");
        let d1 = divisor_from_point(&pts[0].0, &pts[0].1, &fp2);
        let d2 = divisor_from_point(&pts[1].0, &pts[1].1, &fp2);
        let d3 = divisor_from_point(&pts[2].0, &pts[2].1, &fp2);
        // Commutativity
        assert_eq!(add(&c, &d1, &d2, &fp2), add(&c, &d2, &d1, &fp2));
        // Associativity: (D1 + D2) + D3 = D1 + (D2 + D3)
        let lhs = add(&c, &add(&c, &d1, &d2, &fp2), &d3, &fp2);
        let rhs = add(&c, &d1, &add(&c, &d2, &d3, &fp2), &fp2);
        assert_eq!(lhs, rhs);
        // Both should pass Mumford validity.
        assert!(lhs.is_valid(&c, &fp2));
    }

    #[test]
    fn scalar_mult_matches_repeated_add() {
        let fp2 = ctx();
        let c = test_curve(&fp2);
        let mut pts = vec![];
        for k in 0..50i64 {
            if let Some(p) = sample_point(&c, k, &fp2) {
                pts.push(p);
                if pts.len() >= 2 { break; }
            }
        }
        let d1 = divisor_from_point(&pts[0].0, &pts[0].1, &fp2);
        let d2 = divisor_from_point(&pts[1].0, &pts[1].1, &fp2);
        let big = add(&c, &d1, &d2, &fp2);
        // Compute [n]D by repeated addition and by scalar_mul; compare.
        let mut acc = Div::identity(&fp2);
        for n in 1..=10u64 {
            acc = add(&c, &acc, &big, &fp2);
            let smul = scalar_mul(&c, &BigInt::from(n), &big, &fp2);
            assert_eq!(acc, smul, "[{n}]D mismatch");
        }
    }

    /// #C(F_p) by brute enumeration. f acts on F_p (real-valued) since our
    /// test curve has F_p coefficients.
    fn count_points_fp(c: &Curve, fp2: &Fp2) -> u64 {
        let p_u: u64 = fp2.fp.p.to_string().parse().unwrap();
        let mut n: u64 = 1; // ∞
        for x in 0..p_u {
            let xe = fp2.from_int(x as i64);
            let rhs = c.f.eval(&xe, fp2);
            if rhs.b != BigInt::zero() { unreachable!("test curve is F_p-rational"); }
            if rhs.a.is_zero() { n += 1; }
            else if fp2.fp.is_square(&rhs.a) { n += 2; }
        }
        n
    }

    /// #C(F_{p²}) by brute enumeration.
    fn count_points_fp2(c: &Curve, fp2: &Fp2) -> u64 {
        let p_u: u64 = fp2.fp.p.to_string().parse().unwrap();
        let mut n: u64 = 1;
        for a in 0..p_u {
            for b in 0..p_u {
                let x = F2 { a: BigInt::from(a), b: BigInt::from(b) };
                let rhs = c.f.eval(&x, fp2);
                if fp2.is_zero(&rhs) { n += 1; }
                else {
                    // is_square in F_{p²}: rhs^((p²−1)/2) == 1
                    let exp = (&fp2.fp.p * &fp2.fp.p - BigInt::from(1)) / 2;
                    if fp2.pow(&rhs, &exp) == fp2.one() { n += 2; }
                }
            }
        }
        n
    }

    #[test]
    fn order_annihilates() {
        let fp2 = ctx();
        let c = test_curve(&fp2);
        // #J(F_p) = p² + 1 + (p+1)·a₁ + a₂
        //   with a₁ = #C(F_p) − (p+1),  a₂ = (#C(F_{p²}) − p² − 1 − a₁²)/2
        let n1 = count_points_fp(&c, &fp2) as i128;
        let n2 = count_points_fp2(&c, &fp2) as i128;
        let p = 431i128;
        let a1 = n1 - (p + 1);
        // c₂ = (s₁² − s₂)/2 with s₁ = −a₁, s₂ = p² + 1 − N₂
        //   ⇒ c₂ = (a₁² + N₂ − p² − 1)/2
        let a2 = (a1 * a1 + n2 - p * p - 1) / 2;
        let n_j = p * p + 1 + (p + 1) * a1 + a2;
        println!("#C(F_p)={n1}  #C(F_p²)={n2}  a1={a1} a2={a2}  #J(F_p)={n_j}");
        let n_j_big = BigInt::from(n_j as u64);

        // Build a non-trivial divisor and verify [#J]·D = 0.
        let mut pts = vec![];
        for k in 0..200i64 {
            if let Some(p) = sample_point(&c, k, &fp2) {
                pts.push(p);
                if pts.len() >= 2 { break; }
            }
        }
        let d = add(
            &c,
            &divisor_from_point(&pts[0].0, &pts[0].1, &fp2),
            &divisor_from_point(&pts[1].0, &pts[1].1, &fp2),
            &fp2,
        );
        let annihilated = scalar_mul(&c, &n_j_big, &d, &fp2);
        assert!(
            annihilated.is_identity(&fp2),
            "[#J]·D should be 0, got u={:?} v={:?}", annihilated.u, annihilated.v
        );
    }

    #[test]
    fn distributivity() {
        let fp2 = ctx();
        let c = test_curve(&fp2);
        let mut pts = vec![];
        for k in 0..200i64 {
            if let Some(p) = sample_point(&c, k, &fp2) {
                pts.push(p);
                if pts.len() >= 4 { break; }
            }
        }
        let d = add(
            &c,
            &divisor_from_point(&pts[0].0, &pts[0].1, &fp2),
            &divisor_from_point(&pts[1].0, &pts[1].1, &fp2),
            &fp2,
        );
        for (a, b) in [(3u64, 5), (7, 11), (12, 13)] {
            let lhs = scalar_mul(&c, &BigInt::from(a + b), &d, &fp2);
            let rhs = add(
                &c,
                &scalar_mul(&c, &BigInt::from(a), &d, &fp2),
                &scalar_mul(&c, &BigInt::from(b), &d, &fp2),
                &fp2,
            );
            assert_eq!(lhs, rhs, "[a+b]D ≠ [a]D + [b]D for a={a} b={b}");
            // [a][b]D = [a·b]D
            let inner = scalar_mul(&c, &BigInt::from(b), &d, &fp2);
            let nested = scalar_mul(&c, &BigInt::from(a), &inner, &fp2);
            let direct = scalar_mul(&c, &BigInt::from(a * b), &d, &fp2);
            assert_eq!(nested, direct, "[a]([b]D) ≠ [a·b]D for a={a} b={b}");
        }
    }
}
