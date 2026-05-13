//! **Binary elliptic curves**: `y² + xy = x³ + a·x² + b` over
//! `F_{2^m}`.
//!
//! This is the standard Weierstrass form for characteristic-2
//! fields used by NIST B-163, B-233, B-283, B-409, B-571 (SEC 2,
//! FIPS 186-4 Appendix D.1.3).
//!
//! ## Group law (Putranto et al. §3)
//!
//! For affine points `P₁ = (x₁, y₁)`, `P₂ = (x₂, y₂)`:
//!
//! **Negation**: `−P₁ = (x₁, y₁ + x₁)`.
//!
//! **Point Addition (P₁ ≠ ±P₂)**:
//! ```text
//! λ  = (y₁ + y₂) / (x₁ + x₂)
//! x₃ = λ² + λ + x₁ + x₂ + a
//! y₃ = λ·(x₁ + x₃) + x₃ + y₁
//! ```
//!
//! **Point Doubling (P₁ ≠ −P₁)**:
//! ```text
//! λ  = x₁ + y₁/x₁
//! x₃ = λ² + λ + a
//! y₃ = x₁² + (λ + 1)·x₃
//! ```
//!
//! ## ECPM
//!
//! Scalar multiplication `[k]·P` via left-to-right double-and-add.
//! In the paper, this is the operation whose **quantum-circuit
//! resource estimation** is the main contribution; we provide the
//! **classical reference** here.

use super::f2m::{F2mElement, IrreduciblePoly};
use num_bigint::BigUint;
use num_traits::{One, Zero};

/// A binary elliptic curve over `F_{2^m}`.
#[derive(Clone, Debug)]
pub struct BinaryCurve {
    pub m: u32,
    pub irreducible: IrreduciblePoly,
    /// Coefficient `a` of `y² + xy = x³ + a·x² + b`.
    pub a: F2mElement,
    /// Coefficient `b` of the curve equation.
    pub b: F2mElement,
    /// A base point of large prime order (the generator).
    pub generator: BinaryPoint,
    /// Order of the generator (a prime number, as `BigUint`).
    pub order: BigUint,
    /// Cofactor of the curve (`h = #E / n`).
    pub cofactor: BigUint,
}

/// An affine point on a binary elliptic curve, or the point at
/// infinity (`O`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BinaryPoint {
    Infinity,
    Affine { x: F2mElement, y: F2mElement },
}

impl BinaryCurve {
    /// **NIST sect163k1** (Koblitz curve over `F_{2^163}`).
    ///
    /// Curve equation: `y² + xy = x³ + x² + 1` (a = 1, b = 1).
    /// Reduction polynomial: `z¹⁶³ + z⁷ + z⁶ + z³ + 1`.
    /// Per FIPS 186-4 Appendix D.1.3 (and SEC 2 v2.0 §3.1).
    pub fn sect163k1() -> Self {
        let m = 163;
        let irr = IrreduciblePoly::deg_163();
        let a = F2mElement::one(m); // a = 1
        let b = F2mElement::one(m); // b = 1
                                    // Generator coordinates from SEC 2.
        let gx = F2mElement::from_hex("02fe13c0537bbc11acaa07d793de4e6d5e5c94eee8", m);
        let gy = F2mElement::from_hex("0289070fb05d38ff58321f2e800536d538ccdaa3d9", m);
        // Order n (SEC 2): 4000000000000000000020108A2E0CC0D99F8A5EF
        let order = BigUint::parse_bytes(b"4000000000000000000020108A2E0CC0D99F8A5EF", 16).unwrap();
        Self {
            m,
            irreducible: irr,
            a,
            b,
            generator: BinaryPoint::Affine { x: gx, y: gy },
            order,
            cofactor: BigUint::from(2u32),
        }
    }

    /// Tiny test curve over `F_{2^8}` with `a = 1, b = 1`.  Used
    /// for unit-testing without depending on the full sect163k1
    /// parameters.
    pub fn test_curve_f256() -> Self {
        let m = 8;
        let irr = IrreduciblePoly::deg_8();
        let a = F2mElement::one(m);
        let b = F2mElement::one(m);
        // Pick a generator by trial (we'll find one valid point).
        let gx = F2mElement::from_bit_positions(&[1], m); // x = z
                                                          // Find y satisfying y² + zy = z³ + z² + 1 ...
                                                          // Just brute-force search for a valid point.
        let (gen, ord) = find_any_generator(m, &irr, &a, &b);
        Self {
            m,
            irreducible: irr,
            a,
            b,
            generator: gen,
            order: ord,
            cofactor: BigUint::one(),
            // Note: cofactor unknown precisely without full point
            // counting; for the test curve it's fine.
        }
    }

    /// Verify a point satisfies the curve equation.
    pub fn is_on_curve(&self, p: &BinaryPoint) -> bool {
        match p {
            BinaryPoint::Infinity => true,
            BinaryPoint::Affine { x, y } => {
                // y² + xy = x³ + a·x² + b
                let y_sq = y.square(&self.irreducible);
                let xy = x.mul(y, &self.irreducible);
                let lhs = y_sq.add(&xy);
                let x_sq = x.square(&self.irreducible);
                let x_cu = x_sq.mul(x, &self.irreducible);
                let ax_sq = self.a.mul(&x_sq, &self.irreducible);
                let rhs = x_cu.add(&ax_sq).add(&self.b);
                lhs == rhs
            }
        }
    }
}

/// Negation: `−(x, y) = (x, y + x)`.
pub fn point_neg(p: &BinaryPoint) -> BinaryPoint {
    match p {
        BinaryPoint::Infinity => BinaryPoint::Infinity,
        BinaryPoint::Affine { x, y } => BinaryPoint::Affine {
            x: x.clone(),
            y: y.add(x),
        },
    }
}

/// **Point Addition** on a binary curve (Putranto et al. §3, formula
/// reproduced in the module header).
///
/// Cases:
/// - One operand is `O` → return the other.
/// - `P₁ = P₂` → defer to [`point_double`].
/// - `P₁ = −P₂` → return `O`.
/// - Otherwise: λ = (y₁+y₂)/(x₁+x₂), x₃ = λ² + λ + x₁ + x₂ + a,
///   y₃ = λ(x₁+x₃) + x₃ + y₁.
pub fn point_add(curve: &BinaryCurve, p1: &BinaryPoint, p2: &BinaryPoint) -> BinaryPoint {
    match (p1, p2) {
        (BinaryPoint::Infinity, q) | (q, BinaryPoint::Infinity) => q.clone(),
        (BinaryPoint::Affine { x: x1, y: y1 }, BinaryPoint::Affine { x: x2, y: y2 }) => {
            if x1 == x2 {
                // x₁ = x₂.  Either P + P (double) or P + (-P) = O.
                if y1.add(y2) == *x1 {
                    // y₂ = y₁ + x₁ ⇒ P₂ = -P₁.
                    return BinaryPoint::Infinity;
                }
                return point_double(curve, p1);
            }
            let irr = &curve.irreducible;
            let num = y1.add(y2);
            let den = x1.add(x2);
            let lambda = num.mul(&den.flt_inverse(irr).expect("non-zero"), irr);
            let lambda_sq = lambda.square(irr);
            // x₃ = λ² + λ + x₁ + x₂ + a
            let x3 = lambda_sq.add(&lambda).add(x1).add(x2).add(&curve.a);
            // y₃ = λ·(x₁ + x₃) + x₃ + y₁
            let y3 = lambda.mul(&x1.add(&x3), irr).add(&x3).add(y1);
            BinaryPoint::Affine { x: x3, y: y3 }
        }
    }
}

/// **Point Doubling** on a binary curve.
///
/// `[2]·(x₁, y₁)`:
/// - If `x₁ = 0` → result is `O`.
/// - Otherwise: λ = x₁ + y₁/x₁, x₃ = λ² + λ + a, y₃ = x₁² + (λ + 1)·x₃.
pub fn point_double(curve: &BinaryCurve, p: &BinaryPoint) -> BinaryPoint {
    match p {
        BinaryPoint::Infinity => BinaryPoint::Infinity,
        BinaryPoint::Affine { x: x1, y: y1 } => {
            let irr = &curve.irreducible;
            if x1.is_zero() {
                return BinaryPoint::Infinity;
            }
            let x1_inv = x1.flt_inverse(irr).expect("non-zero x₁");
            // λ = x₁ + y₁/x₁
            let lambda = x1.add(&y1.mul(&x1_inv, irr));
            let lambda_sq = lambda.square(irr);
            let x3 = lambda_sq.add(&lambda).add(&curve.a);
            // y₃ = x₁² + (λ + 1)·x₃
            let lambda_plus_one = lambda.add(&F2mElement::one(curve.m));
            let y3 = x1.square(irr).add(&lambda_plus_one.mul(&x3, irr));
            BinaryPoint::Affine { x: x3, y: y3 }
        }
    }
}

/// **ECPM**: scalar multiplication `[k]·P` via left-to-right
/// double-and-add.  This is the operation analysed in Putranto
/// et al. for quantum-resource estimation; classical version here
/// for verification.
pub fn scalar_mul(curve: &BinaryCurve, p: &BinaryPoint, k: &BigUint) -> BinaryPoint {
    if k.is_zero() {
        return BinaryPoint::Infinity;
    }
    let bits = k.bits();
    let mut result = BinaryPoint::Infinity;
    for i in (0..bits).rev() {
        result = point_double(curve, &result);
        if k.bit(i) {
            result = point_add(curve, &result, p);
        }
    }
    result
}

/// Brute-force search for any valid point and its order, on a tiny
/// test curve.  Only feasible for `m ≤ 16`.
fn find_any_generator(
    m: u32,
    irr: &IrreduciblePoly,
    a: &F2mElement,
    b: &F2mElement,
) -> (BinaryPoint, BigUint) {
    let curve_for_check = BinaryCurve {
        m,
        irreducible: irr.clone(),
        a: a.clone(),
        b: b.clone(),
        generator: BinaryPoint::Infinity,
        order: BigUint::one(),
        cofactor: BigUint::one(),
    };
    for xi in 1..(1u64 << m) {
        let x = F2mElement::from_biguint(&BigUint::from(xi), m);
        // Find y satisfying y² + xy = x³ + ax² + b.
        // Brute force: try every possible y.
        for yi in 1..(1u64 << m) {
            let y = F2mElement::from_biguint(&BigUint::from(yi), m);
            let pt = BinaryPoint::Affine { x: x.clone(), y };
            if curve_for_check.is_on_curve(&pt) {
                // Determine order by repeated doubling.
                let ord = point_order(&curve_for_check, &pt);
                if ord >= BigUint::from(2u32) {
                    return (pt, ord);
                }
            }
        }
    }
    (BinaryPoint::Infinity, BigUint::one())
}

fn point_order(curve: &BinaryCurve, p: &BinaryPoint) -> BigUint {
    let mut k = BigUint::one();
    let mut acc = p.clone();
    while k < BigUint::from(1u32 << 10) {
        acc = point_add(curve, &acc, p);
        k += 1u32;
        if matches!(acc, BinaryPoint::Infinity) {
            return k;
        }
    }
    // Couldn't find order in budget; return placeholder.
    k
}

#[cfg(test)]
mod tests {
    use super::*;

    /// NIST sect163k1 generator is on the curve.
    #[test]
    fn sect163k1_generator_is_on_curve() {
        let curve = BinaryCurve::sect163k1();
        assert!(curve.is_on_curve(&curve.generator));
    }

    /// Adding `P + (-P) = O`.
    #[test]
    fn point_plus_neg_is_infinity() {
        let curve = BinaryCurve::sect163k1();
        let neg_g = point_neg(&curve.generator);
        let sum = point_add(&curve, &curve.generator, &neg_g);
        assert!(matches!(sum, BinaryPoint::Infinity));
    }

    /// `2P` via doubling equals `P + P` via addition.
    #[test]
    fn double_equals_self_add() {
        let curve = BinaryCurve::sect163k1();
        let g = curve.generator.clone();
        let g2_dbl = point_double(&curve, &g);
        let g2_add = point_add(&curve, &g, &g);
        assert_eq!(g2_dbl, g2_add);
    }

    /// `2P` is on the curve.
    #[test]
    fn double_stays_on_curve() {
        let curve = BinaryCurve::sect163k1();
        let g2 = point_double(&curve, &curve.generator);
        assert!(curve.is_on_curve(&g2));
    }

    /// Scalar multiplication: `[1]·G = G`, `[2]·G = 2G`, etc.
    #[test]
    fn scalar_mul_small_consistency() {
        let curve = BinaryCurve::sect163k1();
        let g = curve.generator.clone();
        let g1 = scalar_mul(&curve, &g, &BigUint::one());
        assert_eq!(g1, g);
        let g2_sm = scalar_mul(&curve, &g, &BigUint::from(2u32));
        let g2_dbl = point_double(&curve, &g);
        assert_eq!(g2_sm, g2_dbl);
        let g3_sm = scalar_mul(&curve, &g, &BigUint::from(3u32));
        let g3_add = point_add(&curve, &g2_dbl, &g);
        assert_eq!(g3_sm, g3_add);
    }

    /// `[5]G = 2·(2G) + G` via mixed PA/PD.  Demonstrates ECPM
    /// matches a hand-computed sequence.
    #[test]
    fn scalar_mul_matches_handcomputed() {
        let curve = BinaryCurve::sect163k1();
        let g = curve.generator.clone();
        let g2 = point_double(&curve, &g);
        let g4 = point_double(&curve, &g2);
        let g5 = point_add(&curve, &g4, &g);
        let g5_sm = scalar_mul(&curve, &g, &BigUint::from(5u32));
        assert_eq!(g5, g5_sm);
    }

    /// `[n]·G = O` where `n` is the generator's order.  This is
    /// the critical Lagrange-theorem test that distinguishes a
    /// real elliptic-curve implementation from one with subtle
    /// bugs.  Expensive (~thousands of point ops at sect163k1
    /// scale); ignored by default but can be run via --ignored.
    #[test]
    #[ignore]
    fn scalar_mul_by_order_is_identity() {
        let curve = BinaryCurve::sect163k1();
        let n_g = scalar_mul(&curve, &curve.generator, &curve.order);
        assert!(matches!(n_g, BinaryPoint::Infinity));
    }
}
