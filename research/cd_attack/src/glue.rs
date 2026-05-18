//! Castryck-Decru "gluing" step: build a genus-2 hyperelliptic curve from
//! the (2,2)-torsion data of two elliptic curves.
//!
//! Setup: two elliptic curves E_α and E_β, each with three non-identity
//! 2-torsion points giving x-coordinates (α_1, α_2, α_3) and (β_1, β_2, β_3).
//! The construction produces a genus-2 curve h(x) whose Jacobian J(h) is
//! (2,2)-isogenous to E_α × E_β.  The matching pairs the 2-torsion: roughly,
//!   "(α_i, _) glues with (β_i, _)" for i = 1, 2, 3.
//!
//! This is the start of a Castryck-Decru (2^a, 2^a)-isogeny chain on
//! abelian surfaces. After the gluing step, subsequent steps apply Richelot
//! (Jac → Jac) isogenies until the final step which should split back into
//! a product if the input matches a "Kani-glueable" structure.
//!
//! Reference: Castryck-Decru ePrint 2022/975, Magma code `richelot_aux.m`
//! function `FromProdToJac` (lines 1-103).

use crate::field::{F2, Fp2};
use crate::jacobian::Curve;
use crate::poly::Poly;

/// Compute the deg-6 polynomial h(x) for the genus-2 curve gluing E_α and E_β
/// via Kani-style (2,2)-isogeny. The α's and β's must be the x-coords of
/// non-identity 2-torsion points of two short-Weierstrass elliptic curves
/// (three points each, all distinct from each other and across curves
/// — specifically α_i ≠ α_j and β_i ≠ β_j for i ≠ j).
///
/// Returns the codomain Curve { f = h(x) }, with deg h = 6. Same as
/// `from_prod_to_jac_with_partition` but discards the natural partition.
pub fn from_prod_to_jac(alpha: &[F2; 3], beta: &[F2; 3], fp2: &Fp2) -> Curve {
    from_prod_to_jac_with_partition(alpha, beta, fp2).0
}

/// Same as `from_prod_to_jac` but also returns the natural three-quadratic
/// partition (f₁, f₂, f₃) with f₁·f₂·f₃ = h. The "natural Richelot" using
/// this partition splits J(h) back into the elliptic product (Kani inverse).
pub fn from_prod_to_jac_with_partition(
    alpha: &[F2; 3],
    beta: &[F2; 3],
    fp2: &Fp2,
) -> (Curve, [Poly; 3]) {
    // Auxiliary differences (matching the Magma variable naming).
    // d_alp[i][j] = alpha[i] - alpha[j]
    let alp = |i: usize, j: usize| -> F2 { fp2.sub(&alpha[i], &alpha[j]) };
    let bet = |i: usize, j: usize| -> F2 { fp2.sub(&beta[i], &beta[j]) };

    // a1 = (α₃-α₂)²/(β₃-β₂) + (α₂-α₁)²/(β₂-β₁) + (α₁-α₃)²/(β₁-β₃)
    let a1 = fp2.add(
        &fp2.add(
            &fp2.div(&fp2.sq(&alp(2, 1)), &bet(2, 1)),
            &fp2.div(&fp2.sq(&alp(1, 0)), &bet(1, 0)),
        ),
        &fp2.div(&fp2.sq(&alp(0, 2)), &bet(0, 2)),
    );
    // b1: symmetric, swap α↔β
    let b1 = fp2.add(
        &fp2.add(
            &fp2.div(&fp2.sq(&bet(2, 1)), &alp(2, 1)),
            &fp2.div(&fp2.sq(&bet(1, 0)), &alp(1, 0)),
        ),
        &fp2.div(&fp2.sq(&bet(0, 2)), &alp(0, 2)),
    );
    // a2 = α₁(β₃-β₂) + α₂(β₁-β₃) + α₃(β₂-β₁)
    let a2 = fp2.add(
        &fp2.add(
            &fp2.mul(&alpha[0], &bet(2, 1)),
            &fp2.mul(&alpha[1], &bet(0, 2)),
        ),
        &fp2.mul(&alpha[2], &bet(1, 0)),
    );
    // b2 = β₁(α₃-α₂) + β₂(α₁-α₃) + β₃(α₂-α₁)
    let b2 = fp2.add(
        &fp2.add(
            &fp2.mul(&beta[0], &alp(2, 1)),
            &fp2.mul(&beta[1], &alp(0, 2)),
        ),
        &fp2.mul(&beta[2], &alp(1, 0)),
    );

    // Δ_α = (α₁-α₂)² (α₁-α₃)² (α₂-α₃)²
    let delta_alp = fp2.mul(
        &fp2.mul(&fp2.sq(&alp(0, 1)), &fp2.sq(&alp(0, 2))),
        &fp2.sq(&alp(1, 2)),
    );
    let delta_bet = fp2.mul(
        &fp2.mul(&fp2.sq(&bet(0, 1)), &fp2.sq(&bet(0, 2))),
        &fp2.sq(&bet(1, 2)),
    );

    // A = Δ_β · a1 / a2,   B = Δ_α · b1 / b2
    let big_a = fp2.div(&fp2.mul(&delta_bet, &a1), &a2);
    let big_b = fp2.div(&fp2.mul(&delta_alp, &b1), &b2);

    // h(x) = − (A(α₂-α₁)(α₁-α₃)x² + B(β₂-β₁)(β₁-β₃))
    //         · (A(α₃-α₂)(α₂-α₁)x² + B(β₃-β₂)(β₂-β₁))
    //         · (A(α₁-α₃)(α₃-α₂)x² + B(β₁-β₃)(β₃-β₂))
    let factor = |aa: &F2, bb: &F2| -> Poly {
        // x² · aa + bb  →  Poly with coeffs [bb, 0, aa]
        Poly::new(vec![bb.clone(), fp2.zero(), aa.clone()], fp2)
    };

    let f1 = factor(
        &fp2.mul(&big_a, &fp2.mul(&alp(1, 0), &alp(0, 2))),
        &fp2.mul(&big_b, &fp2.mul(&bet(1, 0), &bet(0, 2))),
    );
    let f2 = factor(
        &fp2.mul(&big_a, &fp2.mul(&alp(2, 1), &alp(1, 0))),
        &fp2.mul(&big_b, &fp2.mul(&bet(2, 1), &bet(1, 0))),
    );
    let f3 = factor(
        &fp2.mul(&big_a, &fp2.mul(&alp(0, 2), &alp(2, 1))),
        &fp2.mul(&big_b, &fp2.mul(&bet(0, 2), &bet(2, 1))),
    );

    let prod = f1.mul(&f2, fp2).mul(&f3, fp2);
    // Negate (the leading "−" in the Magma formula).
    let h = prod.neg(fp2);
    // Natural partition (f₁_neg, f₂, f₃) with product = -f₁f₂f₃ = h.
    let f1_neg = f1.neg(fp2);
    (Curve::new(h, fp2), [f1_neg, f2, f3])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{Fp, Fp2};
    use num_bigint::BigInt;
    use num_integer::Integer;
    use num_traits::{One, Zero};

    fn ctx() -> Fp2 { Fp2::new(Fp::new(BigInt::from(431u64))) }

    /// 2-torsion x-coords of a Montgomery curve y² = x³ + Ax² + x. These are
    /// the roots of x(x² + Ax + 1) (excluding the point at infinity).
    /// Returns the three points (x = 0, and the two roots of x² + Ax + 1).
    fn mont_2_torsion_xs(big_a: &F2, fp2: &Fp2) -> Option<[F2; 3]> {
        // Solve x² + A x + 1 = 0 ⇒ x = (-A ± √(A² - 4)) / 2
        let two = fp2.from_int(2);
        let four = fp2.from_int(4);
        let disc = fp2.sub(&fp2.sq(big_a), &four);
        let sqrt_disc = fp2.sqrt(&disc)?;
        let two_inv = fp2.inv(&two);
        let r1 = fp2.mul(&fp2.add(&fp2.neg(big_a), &sqrt_disc), &two_inv);
        let r2 = fp2.mul(&fp2.sub(&fp2.neg(big_a), &sqrt_disc), &two_inv);
        Some([fp2.zero(), r1, r2])
    }

    #[test]
    fn glue_produces_deg6_curve() {
        let fp2 = ctx();
        // E_α: Montgomery A=0  →  2-torsion at 0, ±i  (i² = -1 in F_{p²}).
        let alpha = mont_2_torsion_xs(&fp2.zero(), &fp2).expect("E_α 2-torsion");
        // E_β: Montgomery A=6  →  2-torsion at 0, 52, 373  (computed by hand).
        let beta = mont_2_torsion_xs(&fp2.from_int(6), &fp2).expect("E_β 2-torsion");
        let c = from_prod_to_jac(&alpha, &beta, &fp2);
        assert_eq!(c.f.degree(), Some(6), "glued curve must be degree 6");
    }

    #[test]
    fn glue_curve_has_distinct_roots() {
        let fp2 = ctx();
        let alpha = mont_2_torsion_xs(&fp2.zero(), &fp2).unwrap();
        let beta = mont_2_torsion_xs(&fp2.from_int(6), &fp2).unwrap();
        let c = from_prod_to_jac(&alpha, &beta, &fp2);
        // Check disc(h) ≠ 0  ⇔  gcd(h, h') = constant.
        let h_prime = c.f.derivative(&fp2);
        let g = crate::poly::gcd(&c.f, &h_prime, &fp2);
        assert_eq!(g.degree(), Some(0),
                   "glued curve has repeated root: gcd(h, h') has degree {:?}",
                   g.degree());
    }

    /// Point count on a Montgomery curve y² = x(x² + Ax + 1) over F_p, when A ∈ F_p.
    fn count_mont_fp(big_a: &BigInt, fp: &Fp) -> u64 {
        let p_u: u64 = fp.p.to_string().parse().unwrap();
        let mut n: u64 = 1; // ∞
        for x_int in 0..p_u {
            let x = BigInt::from(x_int);
            let x_sq = fp.mul(&x, &x);
            let rhs = fp.mul(&x, &fp.add(&fp.add(&x_sq, &fp.mul(big_a, &x)), &BigInt::one()));
            if rhs.is_zero() { n += 1; }
            else if fp.is_square(&rhs) { n += 2; }
        }
        n
    }

    /// Point count of #C(F_p) for our curve h ∈ F_p[x] (treating h with F_p coefs).
    fn count_jac_curve_fp(c: &Curve, fp2: &Fp2) -> u64 {
        let p_u: u64 = fp2.fp.p.to_string().parse().unwrap();
        let mut n: u64 = 0;
        // h has degree 6, so 2 points at infinity if leading is a square in F_p, else 0.
        let lc = c.f.leading(fp2);
        if lc.b.is_zero() && fp2.fp.is_square(&lc.a) { n += 2; }
        // else: 0 (or 1 if non-square but we'd get F_{p²}-points which we're not counting)
        for x_int in 0..p_u {
            let x = fp2.from_int(x_int as i64);
            let rhs = c.f.eval(&x, fp2);
            assert!(rhs.b.is_zero(), "test expects h ∈ F_p[x]");
            if rhs.a.is_zero() { n += 1; }
            else if fp2.fp.is_square(&rhs.a) { n += 2; }
        }
        n
    }

    /// L-polynomial Jacobian order from #C(F_p), #C(F_{p²}).
    /// #J(F_p) = p² + 1 + (p+1)a₁ + a₂, a₁ = N₁ - (p+1), a₂ = (a₁² + N₂ - p² - 1)/2.
    fn jac_order_for_fp_curve(c: &Curve, fp2: &Fp2) -> i128 {
        let n1 = count_jac_curve_fp(c, fp2) as i128;
        // For #C(F_{p²}), iterate all x ∈ F_{p²}. Slow but okay for p = 431.
        let p_u: u64 = fp2.fp.p.to_string().parse().unwrap();
        let lc = c.f.leading(fp2);
        // Deg-6: leading is non-zero; F_{p²}-square check
        let exp = (&fp2.fp.p * &fp2.fp.p - BigInt::one()) / 2;
        let mut n2: u64 = 0;
        // Infinity points: 2 if leading is a square in F_{p²} (always true for non-zero).
        if !fp2.is_zero(&lc) {
            // Every non-zero element of F_{p²} is a square in F_{p²} when |F_{p²}*| is even,
            // which holds whenever p > 2. Specifically, x is a square iff x^((p²-1)/2) = 1.
            if fp2.pow(&lc, &exp) == fp2.one() { n2 += 2; }
        }
        for a_u in 0..p_u {
            for b_u in 0..p_u {
                let x = F2 { a: BigInt::from(a_u), b: BigInt::from(b_u) };
                let rhs = c.f.eval(&x, fp2);
                if fp2.is_zero(&rhs) { n2 += 1; }
                else if fp2.pow(&rhs, &exp) == fp2.one() { n2 += 2; }
            }
        }
        let p = 431i128;
        let a1 = n1 - (p + 1);
        let a2 = (a1 * a1 + n2 as i128 - p * p - 1) / 2;
        p * p + 1 + (p + 1) * a1 + a2
    }

    /// Build a divisor on J(h) from two distinct affine x-coords on the curve.
    fn random_deg2_div_on(c: &Curve, fp2: &Fp2, x1_try: i64, x2_try: i64) -> Option<crate::jacobian::Div> {
        let lift = |k: i64| -> Option<(F2, F2)> {
            let x = fp2.from_int(k);
            let rhs = c.f.eval(&x, fp2);
            if !rhs.b.is_zero() { return None; }
            if rhs.a.is_zero() { return None; }
            let y_re = fp2.fp.sqrt(&rhs.a)?;
            Some((x, F2 { a: y_re, b: BigInt::zero() }))
        };
        let (x1, y1) = lift(x1_try)?;
        let (x2, y2) = lift(x2_try)?;
        if x1 == x2 { return None; }
        // u(x) = (x - x1)(x - x2)
        let neg_x1 = fp2.neg(&x1);
        let neg_x2 = fp2.neg(&x2);
        let u = Poly::new(vec![
            fp2.mul(&x1, &x2),
            fp2.add(&neg_x1, &neg_x2),
            fp2.one(),
        ], fp2);
        // v(x) linear with v(x1) = y1, v(x2) = y2
        // v(x) = (y2 - y1)/(x2 - x1) * (x - x1) + y1
        let slope = fp2.div(&fp2.sub(&y2, &y1), &fp2.sub(&x2, &x1));
        let intercept = fp2.sub(&y1, &fp2.mul(&slope, &x1));
        let v = Poly::new(vec![intercept, slope], fp2);
        Some(crate::jacobian::Div { u, v })
    }

    #[test]
    fn cantor_works_on_glued_deg6_curve() {
        // Build the glued curve and verify Cantor arithmetic (add, double,
        // scalar mult) works on it, with the resulting divisors satisfying
        // the Mumford invariant.
        let fp2 = ctx();
        let big_a_alp = fp2.from_int(6);
        let big_a_bet = fp2.from_int(10);
        let alpha = mont_2_torsion_xs(&big_a_alp, &fp2).expect("E_α");
        let beta = mont_2_torsion_xs(&big_a_bet, &fp2).expect("E_β");
        let c = from_prod_to_jac(&alpha, &beta, &fp2);
        assert_eq!(c.f.degree(), Some(6));

        // Build a random deg-2 divisor by sampling x-coords until we find
        // F_p-rational points on J(h). h ∈ F_p[x] (coefs in F_p) so this works.
        let d1 = (3..300i64)
            .filter_map(|k| random_deg2_div_on(&c, &fp2, k, k + 7))
            .next()
            .expect("couldn't find a deg-2 divisor on glued curve");
        assert!(d1.is_valid(&c, &fp2), "constructed divisor must be valid Mumford");

        // [n] D1 must be valid for various small n.
        for n in [2u64, 3, 5, 7, 13, 17] {
            let nd = crate::jacobian::scalar_mul(&c, &BigInt::from(n), &d1, &fp2);
            assert!(
                nd.is_valid(&c, &fp2),
                "[{n}]D must be valid on J(h) (deg-6 curve); got u={:?} v={:?}",
                nd.u, nd.v
            );
        }

        // [2]D = D + D
        let two_d = crate::jacobian::add(&c, &d1, &d1, &fp2);
        let smul2 = crate::jacobian::scalar_mul(&c, &BigInt::from(2u64), &d1, &fp2);
        assert_eq!(two_d, smul2, "[2]D = D+D on glued curve");
    }

    #[test]
    fn natural_partition_splits_glued_jacobian() {
        // The Kani-glued J(h) IS (2,2)-isogenous to E_α × E_β. The natural
        // partition (f₁, f₂, f₃) returned by `from_prod_to_jac_with_partition`
        // should give the inverse isogeny, which means δ = 0 (the partition
        // is split → codomain decomposes as elliptic product).
        let fp2 = ctx();
        let alpha = mont_2_torsion_xs(&fp2.from_int(6), &fp2).expect("E_α");
        let beta = mont_2_torsion_xs(&fp2.from_int(10), &fp2).expect("E_β");
        let (c, partition) = from_prod_to_jac_with_partition(&alpha, &beta, &fp2);

        // Verify partition product = h
        let prod = partition[0].mul(&partition[1], &fp2).mul(&partition[2], &fp2);
        assert_eq!(prod, c.f, "partition factors must multiply to h");

        // Each factor should be quadratic
        for f in &partition {
            assert_eq!(f.degree(), Some(2),
                       "natural partition factors must be quadratic; got {:?}", f.degree());
        }

        // δ = 0 ⇔ Kani inverse split
        let splitting = crate::richelot::Splitting { g: partition };
        let d = crate::richelot::delta(&splitting, &fp2);
        assert!(
            fp2.is_zero(&d),
            "Kani natural partition must give δ = 0 (split back to product); got δ = {:?}", d
        );
    }

    #[test]
    fn glue_preserves_jacobian_order() {
        // The (2,2)-isogeny J(h) → E_α × E_β preserves point count:
        //   #J(h)(F_p) = #E_α(F_p) · #E_β(F_p)
        let fp2 = ctx();
        // Choose both curves with F_p-rational 2-torsion so coefficients of h
        // stay in F_p and point counting works on F_p only.
        // E_α: A=6 → 2-torsion {0, 52, 373} ⊂ F_p
        // E_β: A=10 → discriminant 96, 96 mod 431 = 96. is_square(96, 431)?
        //   96 = 2^5 · 3, (96/431) = (2/431)^5 (3/431) = 1·1 = 1.  Square.
        // So E_β with A=10 has F_p-rational 2-torsion too.
        let big_a_alp = fp2.from_int(6);
        let big_a_bet = fp2.from_int(10);
        let alpha = mont_2_torsion_xs(&big_a_alp, &fp2).expect("E_α");
        let beta = mont_2_torsion_xs(&big_a_bet, &fp2).expect("E_β");
        // Ensure each α_i, β_i is F_p-rational
        for v in alpha.iter().chain(beta.iter()) {
            assert!(v.b.is_zero(), "torsion x-coord must be in F_p for this test");
        }
        let c = from_prod_to_jac(&alpha, &beta, &fp2);

        let e_alp_order = count_mont_fp(&BigInt::from(6), &fp2.fp) as i128;
        let e_bet_order = count_mont_fp(&BigInt::from(10), &fp2.fp) as i128;
        let expected_j = e_alp_order * e_bet_order;
        let actual_j = jac_order_for_fp_curve(&c, &fp2);
        println!(
            "#E_α(F_p) = {}, #E_β(F_p) = {}, product = {}, #J(h) = {}",
            e_alp_order, e_bet_order, expected_j, actual_j
        );
        assert_eq!(
            actual_j, expected_j,
            "Kani gluing: #J(h) should equal #E_α · #E_β"
        );
    }
}
