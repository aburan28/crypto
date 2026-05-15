//! # Semaev summation polynomials for binary elliptic curves.
//!
//! Companion to [`crate::cryptanalysis::ec_index_calculus`] (which targets
//! prime-field curves `y² = x³ + ax + b`).  This module handles the
//! characteristic-2 case `E: y² + xy = x³ + a x² + b` over `F_{2^m}`,
//! which is what shows up in subfield-defined w-curves like the X9.62
//! `c2pnb*w1` family and in any binary-curve Weil-descent attack target.
//!
//! ## The polynomial
//!
//! Three points `P_i = (x_i, y_i) ∈ E(F_{2^m})` sum to `O` iff they are
//! collinear, i.e. lie on a single line `y = λ·x + μ`.  Substituting
//! into the curve equation and reducing in char 2:
//!
//! ```text
//!     x³ + (λ² + λ + a) x² + μ x + (μ² + b) = 0
//! ```
//!
//! The three `x_i` are exactly the roots of this cubic, so Vieta gives
//!
//! ```text
//!     s₁ := x₁ + x₂ + x₃        = λ² + λ + a
//!     s₂ := x₁x₂ + x₁x₃ + x₂x₃  = μ
//!     s₃ := x₁ · x₂ · x₃         = μ² + b
//! ```
//!
//! Eliminating `μ`:  `s₃ = s₂² + b`.  Therefore the algebraic vanishing
//! condition for `∃ y_i: P₁ + P₂ + P₃ = O` is
//!
//! ```text
//!     S₃(x₁, x₂, x₃) := s₃ + s₂² + b
//!                    = x₁x₂x₃ + (x₁x₂ + x₁x₃ + x₂x₃)² + b
//!                    = (x₁+x₂)² · x₃² + x₁x₂ · x₃ + (x₁x₂)² + b.
//! ```
//!
//! Note that `S₃` does **not** depend on the curve parameter `a` — the
//! `a`-contribution lives entirely in the *side condition*
//! `Tr_{F_{2^m}/F_2}(s₁ + a) = 0`, which is the Artin-Schreier
//! solvability of `t² + t = s₁ + a` for the slope `λ`.  Callers using
//! `S₃ = 0` as a sieve test must verify the trace condition separately
//! before back-solving `y_i`.
//!
//! ## What we ship
//!
//! - [`binary_semaev_s3`] — evaluate `S₃(x₁, x₂, x₃)`.
//! - [`binary_semaev_s3_in_x3`] — return `(A, B, C)` with `S₃ = A·X₃² +
//!   B·X₃ + C` as a quadratic in `X₃`, ready for an index-calculus sweep
//!   that fixes `(x₁, x₂)` and asks "is there an `x₃` in the factor base
//!   with `S₃ = 0`?"
//! - [`solve_artin_schreier`] — solve `t² + t = c` over `F_{2^m}` via
//!   the half-trace formula (odd `m`) with a brute-force fallback for
//!   even `m`.
//! - [`solve_quadratic_f2m`] — solve `A·X² + B·X + C = 0` over `F_{2^m}`
//!   by the standard `X = (B/A)·z` substitution to an AS equation.
//!
//! ## References
//!
//! - I. Semaev, *Summation polynomials and the discrete logarithm
//!   problem on elliptic curves*, IACR ePrint 2004/031.
//! - C. Diem, *On the discrete logarithm problem in elliptic curves
//!   over non-prime finite fields*, LMS J. Comput. Math. 14 (2011).
//! - J.-C. Faugère, P. Gaudry, L. Huot, G. Renault, *Sub-cubic change
//!   of ordering for Gröbner basis*, ISSAC 2013 — the binary-curve
//!   S₃ is in §6.

use crate::binary_ecc::{F2mElement, IrreduciblePoly};

// ── Semaev S₃ for E: y² + xy = x³ + a x² + b ───────────────────────

/// Evaluate `S₃(x₁, x₂, x₃) = (x₁+x₂)² x₃² + x₁x₂ x₃ + (x₁x₂)² + b`.
///
/// Returns `0` iff there exist `y_i ∈ F_{2^m}` (possibly in a quadratic
/// extension) with `(x_i, y_i) ∈ E` and `±P₁ ± P₂ ± P₃ = O`.  In the
/// binary-curve form `y² + xy = …` the "negation" of `P = (x, y)` is
/// `(x, x + y)`, so the sign ambiguity collapses into the choice of
/// `y_i` from the two roots of `T² + x_i T + (x_i³ + a x_i² + b) = 0`.
pub fn binary_semaev_s3(
    x1: &F2mElement,
    x2: &F2mElement,
    x3: &F2mElement,
    b: &F2mElement,
    irr: &IrreduciblePoly,
) -> F2mElement {
    let sum12 = x1.add(x2);
    let sum12_sq = sum12.square(irr);
    let prod12 = x1.mul(x2, irr);
    let prod12_sq = prod12.square(irr);
    let x3_sq = x3.square(irr);
    let t1 = sum12_sq.mul(&x3_sq, irr);
    let t2 = prod12.mul(x3, irr);
    let t3 = prod12_sq.add(b);
    t1.add(&t2).add(&t3)
}

/// Return `(A, B, C)` such that `S₃(x₁, x₂, X₃) = A · X₃² + B · X₃ + C`
/// when `(x₁, x₂)` are fixed.  Plug into [`solve_quadratic_f2m`] to find
/// candidate `X₃` values.
pub fn binary_semaev_s3_in_x3(
    x1: &F2mElement,
    x2: &F2mElement,
    b: &F2mElement,
    irr: &IrreduciblePoly,
) -> (F2mElement, F2mElement, F2mElement) {
    let sum12 = x1.add(x2);
    let sum12_sq = sum12.square(irr);
    let prod12 = x1.mul(x2, irr);
    let prod12_sq = prod12.square(irr);
    let c = prod12_sq.add(b);
    (sum12_sq, prod12, c)
}

// ── Artin-Schreier and quadratic solvers in F_{2^m} ────────────────

/// Solve `t² + t = c` over `F_{2^m}`.
///
/// In characteristic 2 this Artin-Schreier equation has a solution iff
/// `Tr_{F_{2^m}/F_2}(c) = 0`.  When it does, the two solutions differ
/// by `1`; we return one (the caller adds `1` for the other).
///
/// For odd `m` we use the half-trace formula
/// `t = c + c^{2²} + c^{2⁴} + … + c^{2^{m-1}}`, which gives a valid
/// solution whenever `Tr(c) = 0`.  For even `m` the half-trace recipe
/// requires a fixed non-AS element; we fall back to brute force,
/// which is fine for the toy fields this module targets (`m ≤ 16`).
pub fn solve_artin_schreier(
    c: &F2mElement,
    m: u32,
    irr: &IrreduciblePoly,
) -> Option<F2mElement> {
    // Compute Tr(c) = c + c² + c⁴ + … + c^{2^{m-1}}.
    let mut acc = c.clone();
    let mut tr = c.clone();
    for _ in 1..m {
        acc = acc.square(irr);
        tr = tr.add(&acc);
    }
    let one = F2mElement::one(m);
    if tr == one {
        return None; // Tr = 1, no F_{2^m} solution.
    }
    debug_assert!(tr.is_zero(), "trace should lie in F_2");

    if m % 2 == 1 {
        // Half-trace:  H(c) = c + c^{2²} + c^{2⁴} + … + c^{2^{m-1}}.
        // This satisfies H(c)² + H(c) = c whenever Tr(c) = 0.
        let mut t = c.clone();
        let mut acc = c.clone();
        let half = (m - 1) / 2;
        for _ in 0..half {
            acc = acc.square(irr).square(irr);
            t = t.add(&acc);
        }
        return Some(t);
    }

    // Even m: brute force (sufficient for the toy field sizes here).
    if m > 20 {
        return None;
    }
    for v in 0u64..(1u64 << m) {
        let bits: Vec<u32> = (0..m).filter(|i| (v >> i) & 1 == 1).collect();
        let cand = F2mElement::from_bit_positions(&bits, m);
        let lhs = cand.square(irr).add(&cand);
        if &lhs == c {
            return Some(cand);
        }
    }
    None
}

/// Solve `A · X² + B · X + C = 0` over `F_{2^m}`.
///
/// - `A = 0, B ≠ 0`: linear, single root `X = C / B`.
/// - `A = 0, B = 0`: degenerate (the polynomial is the constant `C`);
///    no usable roots reported.
/// - `A ≠ 0, B = 0`: pure square, `X = √(C/A)`.  In char 2 the square
///    root is the absolute Frobenius `m-1` times.
/// - `A ≠ 0, B ≠ 0`: substitute `X = (B/A) · z` to get
///    `z² + z = A · C / B²`, an Artin-Schreier equation in `z`.
///
/// Returns up to two solutions (no guaranteed order).
pub fn solve_quadratic_f2m(
    a_coef: &F2mElement,
    b_coef: &F2mElement,
    c_coef: &F2mElement,
    m: u32,
    irr: &IrreduciblePoly,
) -> Vec<F2mElement> {
    if a_coef.is_zero() {
        if b_coef.is_zero() {
            return vec![];
        }
        let inv = b_coef.flt_inverse(irr).expect("B ≠ 0");
        return vec![c_coef.mul(&inv, irr)];
    }
    if b_coef.is_zero() {
        // A X² = C  ⇒  X² = C/A  ⇒  X = √(C/A).
        let a_inv = a_coef.flt_inverse(irr).expect("A ≠ 0");
        let val = c_coef.mul(&a_inv, irr);
        let mut r = val;
        for _ in 0..(m - 1) {
            r = r.square(irr);
        }
        return vec![r];
    }
    // Generic case.
    let b_sq = b_coef.square(irr);
    let b_sq_inv = b_sq.flt_inverse(irr).expect("B² ≠ 0");
    let arg = a_coef.mul(c_coef, irr).mul(&b_sq_inv, irr);
    match solve_artin_schreier(&arg, m, irr) {
        None => vec![],
        Some(z) => {
            let one = F2mElement::one(m);
            let z_other = z.add(&one);
            let scale = b_coef.mul(&a_coef.flt_inverse(irr).expect("A ≠ 0"), irr);
            vec![scale.mul(&z, irr), scale.mul(&z_other, irr)]
        }
    }
}

// ── Tests ──────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_ecc::IrreduciblePoly;
    use crate::cryptanalysis::ghs_descent::{ECurve, Pt};
    use num_bigint::BigUint;

    fn f64_irr() -> IrreduciblePoly {
        IrreduciblePoly { degree: 6, low_terms: vec![0, 1] }
    }

    /// Symmetry: `S₃(x₁, x₂, x₃)` equals every permutation.
    #[test]
    fn s3_symmetric() {
        let irr = f64_irr();
        let m = 6;
        let b = F2mElement::from_bit_positions(&[0, 2], m);
        let x1 = F2mElement::from_bit_positions(&[1], m);
        let x2 = F2mElement::from_bit_positions(&[0, 3], m);
        let x3 = F2mElement::from_bit_positions(&[2, 4], m);
        let v123 = binary_semaev_s3(&x1, &x2, &x3, &b, &irr);
        let v132 = binary_semaev_s3(&x1, &x3, &x2, &b, &irr);
        let v213 = binary_semaev_s3(&x2, &x1, &x3, &b, &irr);
        let v231 = binary_semaev_s3(&x2, &x3, &x1, &b, &irr);
        let v312 = binary_semaev_s3(&x3, &x1, &x2, &b, &irr);
        let v321 = binary_semaev_s3(&x3, &x2, &x1, &b, &irr);
        assert_eq!(v123, v132);
        assert_eq!(v123, v213);
        assert_eq!(v123, v231);
        assert_eq!(v123, v312);
        assert_eq!(v123, v321);
    }

    /// `binary_semaev_s3_in_x3` agrees with `binary_semaev_s3` at any `x₃`.
    #[test]
    fn s3_assembled_form_matches() {
        let irr = f64_irr();
        let m = 6;
        let b = F2mElement::from_bit_positions(&[1, 4], m);
        let x1 = F2mElement::from_bit_positions(&[0, 3], m);
        let x2 = F2mElement::from_bit_positions(&[2], m);
        let (a_coef, b_coef, c_coef) = binary_semaev_s3_in_x3(&x1, &x2, &b, &irr);
        for v in 0u64..(1u64 << m) {
            let bits: Vec<u32> = (0..m).filter(|i| (v >> i) & 1 == 1).collect();
            let x3 = F2mElement::from_bit_positions(&bits, m);
            let direct = binary_semaev_s3(&x1, &x2, &x3, &b, &irr);
            let x3_sq = x3.square(&irr);
            let assembled = a_coef
                .mul(&x3_sq, &irr)
                .add(&b_coef.mul(&x3, &irr))
                .add(&c_coef);
            assert_eq!(direct, assembled, "mismatch at x3 = {:?}", x3.to_biguint());
        }
    }

    /// Artin-Schreier round-trip: pick t, compute c = t² + t, solve, recover t (or t+1).
    #[test]
    fn artin_schreier_round_trip() {
        let irr = f64_irr(); // F_{2^6}, m = 6, even — exercises brute-force path.
        let m = 6;
        let one = F2mElement::one(m);
        let mut tested = 0;
        for v in 1u64..(1u64 << m) {
            let bits: Vec<u32> = (0..m).filter(|i| (v >> i) & 1 == 1).collect();
            let t = F2mElement::from_bit_positions(&bits, m);
            let c = t.square(&irr).add(&t);
            let solved = solve_artin_schreier(&c, m, &irr)
                .expect("Tr(t² + t) = 0 always, so a solution exists");
            assert!(solved == t || solved == t.add(&one), "round-trip failed");
            tested += 1;
            if tested > 20 {
                break;
            }
        }
    }

    /// Quadratic solver round-trip: build a known-root quadratic, recover the roots.
    #[test]
    fn quadratic_round_trip() {
        let irr = f64_irr();
        let m = 6;
        let r1 = F2mElement::from_bit_positions(&[0, 2, 5], m);
        let r2 = F2mElement::from_bit_positions(&[1, 3], m);
        // (X + r1)(X + r2) = X² + (r1+r2) X + r1·r2.
        let a_coef = F2mElement::one(m);
        let b_coef = r1.add(&r2);
        let c_coef = r1.mul(&r2, &irr);
        let mut roots = solve_quadratic_f2m(&a_coef, &b_coef, &c_coef, m, &irr);
        roots.sort_by_key(|r| r.to_biguint());
        let mut expected = vec![r1, r2];
        expected.sort_by_key(|r| r.to_biguint());
        assert_eq!(roots, expected);
    }

    /// `S₃` vanishes on a collinear triple of E-points.
    ///
    /// Pick `P₁` and `P₂` on a small binary curve `E`, compute
    /// `P₃ = -(P₁ + P₂)` via the group law, then check that
    /// `S₃(x₁, x₂, x₃) = 0`.
    #[test]
    fn s3_vanishes_on_collinear_triple() {
        let irr = f64_irr();
        let m = 6;
        // Pick a curve.
        let a = F2mElement::one(m);
        let b = F2mElement::from_bit_positions(&[0, 1, 4], m);
        let curve = ECurve::new(m, irr.clone(), a, b.clone());
        // Collect a handful of affine points; we'll pick a non-degenerate
        // pair (one whose sum isn't O and which is distinct from one
        // another's negative).
        let mut points = Vec::new();
        for xi in 1u64..(1u64 << m) {
            let x = F2mElement::from_biguint(&BigUint::from(xi), m);
            for yi in 0u64..(1u64 << m) {
                let y = F2mElement::from_biguint(&BigUint::from(yi), m);
                let pt = Pt::Aff { x: x.clone(), y };
                if curve.is_on_curve(&pt) {
                    points.push(pt);
                    if points.len() >= 8 {
                        break;
                    }
                }
            }
            if points.len() >= 8 {
                break;
            }
        }
        assert!(points.len() >= 2, "need at least two points on E");
        // Find a non-degenerate (P₁, P₂) pair.
        let mut tested = false;
        'outer: for i in 0..points.len() {
            for j in (i + 1)..points.len() {
                let p1 = points[i].clone();
                let p2 = points[j].clone();
                let sum = curve.add(&p1, &p2);
                if matches!(sum, Pt::Inf) {
                    continue; // P₂ = -P₁, skip.
                }
                let p3 = curve.neg(&sum);
                let (x1, x2, x3) = match (&p1, &p2, &p3) {
                    (Pt::Aff { x: a, .. }, Pt::Aff { x: b, .. }, Pt::Aff { x: c, .. }) =>
                        (a.clone(), b.clone(), c.clone()),
                    _ => continue,
                };
                // Also skip degenerate doubling (x₁ = x₂).
                if x1 == x2 {
                    continue;
                }
                let s3 = binary_semaev_s3(&x1, &x2, &x3, &b, &irr);
                assert!(
                    s3.is_zero(),
                    "S₃ should vanish on a collinear E-triple; got {:?} at \
                     (x1={:?}, x2={:?}, x3={:?})",
                    s3.to_biguint(),
                    x1.to_biguint(),
                    x2.to_biguint(),
                    x3.to_biguint()
                );
                tested = true;
                break 'outer;
            }
        }
        assert!(tested, "no non-degenerate (P₁, P₂) pair found");
    }
}
