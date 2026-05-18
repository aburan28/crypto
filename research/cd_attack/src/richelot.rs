//! Richelot (2,2)-isogenies on genus-2 Jacobians.
//!
//! Setup: f(x) factors as G_1(x) · G_2(x) · G_3(x) where each G_i is at most
//! quadratic. This partition of the 6 (or 5+∞) Weierstrass points into three
//! pairs is exactly an isotropic 2-dimensional subgroup of J[2] — a Richelot
//! kernel.
//!
//! Construction (Bost-Mestre 1988):
//!   δ      = det(coefficient matrix of G_1, G_2, G_3)
//!   H_i(X) = G_j(X) G_k'(X) − G_j'(X) G_k(X)   (j, k complementary to i)
//!   f̃(X)  = δ⁻¹ · H_1(X) · H_2(X) · H_3(X)        (codomain polynomial)
//!
//! Splitting condition: δ = 0 ⇔ Richelot codomain decomposes as elliptic ×
//! elliptic (this is the splitting test used by the Castryck-Decru attack
//! as its verification signal).

use crate::field::{F2, Fp2};
use crate::jacobian::{self, Curve, Div};
use crate::poly::{ext_gcd, Poly};
use num_bigint::BigInt;

/// A factorization f = G_0 · G_1 · G_2 (we use 0-indexed throughout).
pub struct Splitting {
    pub g: [Poly; 3],
}

impl Splitting {
    /// Verify the product equals the curve's f.
    pub fn matches(&self, c: &Curve, fp2: &Fp2) -> bool {
        let prod = self.g[0].mul(&self.g[1], fp2).mul(&self.g[2], fp2);
        prod == c.f
    }
}

/// Extract (c₂, c₁, c₀) from a polynomial of degree ≤ 2.
fn coeffs3(p: &Poly, fp2: &Fp2) -> (F2, F2, F2) {
    (
        p.coeffs.get(2).cloned().unwrap_or_else(|| fp2.zero()),
        p.coeffs.get(1).cloned().unwrap_or_else(|| fp2.zero()),
        p.coeffs.get(0).cloned().unwrap_or_else(|| fp2.zero()),
    )
}

/// δ = det [[a_i, b_i, c_i]] where G_i = a_i·x² + b_i·x + c_i.
/// Vanishing δ means the codomain factors as a product of elliptic curves.
pub fn delta(s: &Splitting, fp2: &Fp2) -> F2 {
    let (a0, b0, c0) = coeffs3(&s.g[0], fp2);
    let (a1, b1, c1) = coeffs3(&s.g[1], fp2);
    let (a2, b2, c2) = coeffs3(&s.g[2], fp2);
    let m00 = fp2.sub(&fp2.mul(&b1, &c2), &fp2.mul(&b2, &c1));
    let m01 = fp2.sub(&fp2.mul(&a1, &c2), &fp2.mul(&a2, &c1));
    let m02 = fp2.sub(&fp2.mul(&a1, &b2), &fp2.mul(&a2, &b1));
    let t1 = fp2.mul(&a0, &m00);
    let t2 = fp2.mul(&b0, &m01);
    let t3 = fp2.mul(&c0, &m02);
    fp2.add(&fp2.sub(&t1, &t2), &t3)
}

pub fn is_split(s: &Splitting, fp2: &Fp2) -> bool {
    fp2.is_zero(&delta(s, fp2))
}

/// The three "co-Wronskian" polynomials. H[i] omits G[i] and uses the
/// Wronskian of the other two:  H[i] = G[j]·G[k]' − G[j]'·G[k].
/// Returned in same indexing convention (H[0], H[1], H[2]).
pub fn co_wronskians(s: &Splitting, fp2: &Fp2) -> [Poly; 3] {
    let gd: [Poly; 3] = [
        s.g[0].derivative(fp2),
        s.g[1].derivative(fp2),
        s.g[2].derivative(fp2),
    ];
    let w = |j: usize, k: usize| -> Poly {
        s.g[j].mul(&gd[k], fp2).sub(&gd[j].mul(&s.g[k], fp2), fp2)
    };
    [
        w(1, 2), // omit index 0
        w(2, 0), // omit index 1
        w(0, 1), // omit index 2
    ]
}

/// Codomain curve C̃: y² = f̃(x).  Panics if the isogeny is split (δ = 0) —
/// callers should `is_split` first if they want to handle that case.
pub fn codomain(s: &Splitting, fp2: &Fp2) -> Curve {
    let d = delta(s, fp2);
    assert!(!fp2.is_zero(&d), "Richelot isogeny is split (δ = 0)");
    let d_inv = fp2.inv(&d);
    let h = co_wronskians(s, fp2);
    let prod = h[0].mul(&h[1], fp2).mul(&h[2], fp2);
    let f_tilde = prod.scalar_mul(&d_inv, fp2);
    // Do NOT monicize: if the leading coefficient is a non-square in F_p²,
    // monicizing twists the curve and changes #J. Keep f̃ as-is — its
    // Jacobian is the genuine isogenous one.
    Curve::new(f_tilde, fp2)
}

/// Push a divisor D ∈ Jac(C) through the Richelot isogeny to Jac(C̃).
///
/// Correspondence (Castryck-Decru Magma `FromJacToJac`):
///   x-part:  G_1(x) H_1(X) + G_2(x) H_2(X) = 0       (degree 2 in X)
///   y-part:  y₁ · y₂ = G_1(x) H_1(X) · (x − X)
///
/// For D = ⟨u(x), v(x)⟩ with u monic deg 2:
///   1. Reduce Φ(x,X) mod u(x) to linear in x:  α(X) x + β(X)
///   2. U_new(X) = u₀ α² − u₁ α β + β²                (degree ≤ 4)
///   3. V_new(X) = (−G_1*(X) · H_1(X) · (β + αX)) / (α² · v*(X))   mod U_new
///      where G_1* = β² − g₁₁ αβ + g₁₀ α² and v* = v₀α − v₁β
///   4. Cantor-reduce ⟨U_new, V_new⟩ on C̃ to standard Mumford form.
pub fn pushforward(s: &Splitting, d: &Div, fp2: &Fp2) -> (Curve, Div) {
    let c_tilde = codomain(s, fp2);

    // Identity → identity
    if d.is_identity(fp2) {
        return (c_tilde, Div::identity(fp2));
    }

    let du = d.u.degree().unwrap();
    if du != 2 {
        // Degree-1 u (single affine point) would need its own branch.
        // For C-D the kernels we propagate are always deg-2 generically.
        unimplemented!("pushforward only implemented for deg-2 u (got deg {du})");
    }

    let h = co_wronskians(s, fp2);
    let h1 = &h[0];
    let h2 = &h[1];
    let (g1_2, g1_1, g1_0) = coeffs3(&s.g[0], fp2);
    let (g2_2, g2_1, g2_0) = coeffs3(&s.g[1], fp2);

    // Φ(x, X) = G_1(x) H_1(X) + G_2(x) H_2(X)
    //        = A2(X) x² + A1(X) x + A0(X)
    let big_a2 = h1.scalar_mul(&g1_2, fp2).add(&h2.scalar_mul(&g2_2, fp2), fp2);
    let big_a1 = h1.scalar_mul(&g1_1, fp2).add(&h2.scalar_mul(&g2_1, fp2), fp2);
    let big_a0 = h1.scalar_mul(&g1_0, fp2).add(&h2.scalar_mul(&g2_0, fp2), fp2);

    // u(x) = x² + u₁ x + u₀  (monic)
    let u0 = d.u.coeffs[0].clone();
    let u1 = d.u.coeffs[1].clone();

    // Reduce Φ mod u(x): use x² = −u₁x − u₀.
    //   α(X) = A1 − A2·u₁,  β(X) = A0 − A2·u₀
    let alpha = big_a1.sub(&big_a2.scalar_mul(&u1, fp2), fp2);
    let beta = big_a0.sub(&big_a2.scalar_mul(&u0, fp2), fp2);

    // U_new(X) = Res_x(Φ(x,X), u(x)) = u₀ α² − u₁ αβ + β²
    let alpha_sq = alpha.mul(&alpha, fp2);
    let alpha_beta = alpha.mul(&beta, fp2);
    let beta_sq = beta.mul(&beta, fp2);
    let u_new = alpha_sq.scalar_mul(&u0, fp2)
        .sub(&alpha_beta.scalar_mul(&u1, fp2), fp2)
        .add(&beta_sq, fp2);

    // v(x) = v₁ x + v₀
    let v0 = d.v.coeffs.get(0).cloned().unwrap_or_else(|| fp2.zero());
    let v1 = d.v.coeffs.get(1).cloned().unwrap_or_else(|| fp2.zero());

    // G_1*(X) := G_1(x) with x ← −β/α  (cleared denominator α²)
    //          = β² − g₁₁ αβ + g₁₀ α²    (since G_1 = x² + g₁₁x + g₁₀)
    // (Assuming G_1 monic; g1_2 = 1.)
    let g1_star = beta_sq.clone()
        .sub(&alpha_beta.scalar_mul(&g1_1, fp2), fp2)
        .add(&alpha_sq.scalar_mul(&g1_0, fp2), fp2);

    // v*(X) := v(x) with x ← −β/α, cleared α: v* = v₀ α − v₁ β
    let v_star = alpha.scalar_mul(&v0, fp2).sub(&beta.scalar_mul(&v1, fp2), fp2);

    // β + α·X  (linear-in-X factor for the (x − X) term, after substitution)
    let alpha_times_x = alpha.mul(&Poly::x(fp2), fp2);
    let beta_plus_ax = beta.add(&alpha_times_x, fp2);

    // N(X) = −G_1*(X) · H_1(X) · (β + αX)
    let n_poly = g1_star.mul(h1, fp2).mul(&beta_plus_ax, fp2).neg(fp2);
    // D(X) = α² · v*
    let d_poly = alpha_sq.mul(&v_star, fp2);

    // V_new(X) ≡ N · D⁻¹ (mod U_new).
    // The denominator D and the unreduced U_new can share a common factor that
    // comes from the algebraic identity Σ Gᵢ Hᵢ = c·(x−X)²  — geometrically,
    // "ideal" image points that don't lie on a generic fibre. Cancel the gcd
    // from N, D, and U_new before inverting; what remains is the genuine
    // unreduced image divisor.
    // The denominator D and U_new can share a common factor that comes from
    // the algebraic identity Σ Gᵢ Hᵢ ∝ (x−X)² — geometrically, "ideal" image
    // points that don't lie on a generic fibre. When the source divisor has
    // a 2-torsion component, v* vanishes at those image X-values; cancelling
    // the gcd from U_new, N, D removes those ideal points before inversion.
    let g = crate::poly::gcd(&d_poly, &u_new, fp2);
    let (u_new, n_poly, d_poly) = if g.degree() == Some(0) {
        (u_new, n_poly, d_poly)
    } else {
        let (un, r1) = u_new.div_rem(&g, fp2);
        let (np, _r2) = n_poly.div_rem(&g, fp2);
        let (dp, r3) = d_poly.div_rem(&g, fp2);
        assert!(r1.is_zero() && r3.is_zero(),
                "gcd should divide D and U_new exactly");
        (un, np, dp)
    };

    let (gcd2, inv_d, _) = ext_gcd(&d_poly, &u_new, fp2);
    assert!(
        gcd2.degree() == Some(0),
        "D(X) still not invertible mod U_new after gcd cancellation (gcd deg {:?})",
        gcd2.degree()
    );
    let v_unreduced = n_poly.mul(&inv_d, fp2).rem(&u_new, fp2);

    // Cantor reduce on the codomain curve.
    let pushed = jacobian::reduce(&c_tilde, &u_new, &v_unreduced, fp2);
    (c_tilde, pushed)
}

/// Castryck-Decru chain-splitting check (analog of `Does22ChainSplit` from
/// the published Magma reference).
///
/// Given a genus-2 curve `h` (produced by `glue::from_prod_to_jac`) and two
/// divisors `d1, d2` of order 2^(a-1) representing the Kani-mapped torsion
/// data, runs `a-2` Richelot steps and checks whether the final
/// partition of the resulting curve is split (δ = 0).
///
/// Returns:
///   true   — chain runs to completion and final step gives δ = 0
///   false  — premature split or degenerate state somewhere in the chain
///
/// In the C-D attack, a correct guess of Alice's secret bits gives `true`;
/// wrong guesses give `false` because the chain splits prematurely.
pub fn does_22_chain_split(
    h: &Curve,
    d1: &Div,
    d2: &Div,
    a: u32,
    fp2: &Fp2,
) -> bool {
    assert!(a >= 2, "chain length must be at least 2");
    let mut cur_curve = h.clone();
    let mut cur_d1 = d1.clone();
    let mut cur_d2 = d2.clone();

    // Run a-2 intermediate Richelot steps. At iteration i (1-indexed),
    // the chain parameter is a-i, and we double d1, d2 by 2^(a-i-1) to
    // extract their 2-torsion piece.
    for i in 1..=(a as i64 - 2) {
        let chain_param = (a as i64 - i) as u32;
        let pow_exp = chain_param - 1;
        let pow = if pow_exp == 0 {
            BigInt::from(1)
        } else {
            BigInt::from(1u64 << pow_exp)
        };

        let t1 = jacobian::scalar_mul(&cur_curve, &pow, &cur_d1, fp2);
        let t2 = jacobian::scalar_mul(&cur_curve, &pow, &cur_d2, fp2);

        if t1.u.degree() != Some(2) || t2.u.degree() != Some(2) {
            // Degenerate: cannot form deg-2 partition factors.
            return false;
        }
        let g1 = t1.u.clone();
        let g2 = t2.u.clone();
        let g1g2 = g1.mul(&g2, fp2);
        let (g3, rem) = cur_curve.f.div_rem(&g1g2, fp2);
        if !rem.is_zero() {
            return false;
        }
        let s = Splitting { g: [g1, g2, g3] };
        if is_split(&s, fp2) {
            return false; // premature split → wrong guess in C-D
        }
        let (new_curve, new_d1) = pushforward(&s, &cur_d1, fp2);
        let (_, new_d2) = pushforward(&s, &cur_d2, fp2);
        cur_curve = new_curve;
        cur_d1 = new_d1;
        cur_d2 = new_d2;
    }

    // Final step: divisors should now have order exactly 2 (their u-polys
    // ARE the partition factors directly).
    if cur_d1.u.degree() != Some(2) || cur_d2.u.degree() != Some(2) {
        return false;
    }
    let g1 = cur_d1.u.clone();
    let g2 = cur_d2.u.clone();
    let g1g2 = g1.mul(&g2, fp2);
    let (g3, rem) = cur_curve.f.div_rem(&g1g2, fp2);
    if !rem.is_zero() {
        return false;
    }
    let s_final = Splitting { g: [g1, g2, g3] };
    is_split(&s_final, fp2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{Fp, Fp2};
    use num_bigint::BigInt;
    use num_integer::Integer;
    use num_traits::Zero;

    fn ctx() -> Fp2 { Fp2::new(Fp::new(BigInt::from(431u64))) }

    /// Build a deg-2 polynomial from int coefficients [c, b, a] meaning a·x² + b·x + c.
    fn quad(c: i64, b: i64, a: i64, fp2: &Fp2) -> Poly {
        Poly::new(vec![fp2.from_int(c), fp2.from_int(b), fp2.from_int(a)], fp2)
    }
    fn linr(c: i64, b: i64, fp2: &Fp2) -> Poly {
        Poly::new(vec![fp2.from_int(c), fp2.from_int(b)], fp2)
    }

    /// A deg-5 curve with explicit factorization (x²+1)(x²+3)(x−2).
    /// The two quadratics share their x-coefficient (both zero) — this
    /// forces H[2] to drop from generic deg 2 to deg 1, keeping the
    /// codomain f̃ at deg 5 so it fits our current Curve invariant.
    fn split_curve(fp2: &Fp2) -> (Curve, Splitting) {
        let g0 = quad(1, 0, 1, fp2);      // x² + 1
        let g1 = quad(3, 0, 1, fp2);      // x² + 3
        let g2 = linr(-2, 1, fp2);        // x − 2
        let f = g0.mul(&g1, fp2).mul(&g2, fp2);
        let c = Curve::new_monic_deg5(f, fp2);
        let s = Splitting { g: [g0, g1, g2] };
        assert!(s.matches(&c, fp2));
        (c, s)
    }

    #[test]
    fn splitting_matches_product() {
        let fp2 = ctx();
        let (_c, s) = split_curve(&fp2);
        // matches() was used inside the helper; explicit re-check.
        let prod = s.g[0].mul(&s.g[1], &fp2).mul(&s.g[2], &fp2);
        let expected = split_curve(&fp2).0.f;
        assert_eq!(prod, expected);
    }

    #[test]
    fn delta_nonzero_general_case() {
        let fp2 = ctx();
        let (_c, s) = split_curve(&fp2);
        let d = delta(&s, &fp2);
        assert!(!fp2.is_zero(&d), "δ should be nonzero for this splitting");
        assert!(!is_split(&s, &fp2));
    }

    #[test]
    fn delta_zero_when_factors_dependent() {
        let fp2 = ctx();
        // G_2 = α·G_0 + β·G_1 (linear combo) ⇒ δ = 0.
        // Try G_0 = x²+1,  G_1 = x²+2,  G_2 = 2·x² + 3 = G_0 + G_1.
        let g0 = quad(1, 0, 1, &fp2);
        let g1 = quad(2, 0, 1, &fp2);
        // For G_2 to be a literal sum (coefficient-wise): leading 2.
        let g2 = quad(3, 0, 2, &fp2);
        let s = Splitting { g: [g0, g1, g2] };
        let d = delta(&s, &fp2);
        assert!(fp2.is_zero(&d), "linearly dependent rows ⇒ δ = 0");
        assert!(is_split(&s, &fp2));
    }

    #[test]
    fn co_wronskians_match_manual() {
        let fp2 = ctx();
        let (_c, s) = split_curve(&fp2);
        let h = co_wronskians(&s, &fp2);
        // Manual: H[2] = G[0] · G[1]' − G[0]' · G[1] with G[0]=x²+1, G[1]=x²+3.
        //   G[0]' = 2x,  G[1]' = 2x
        //   G[0] G[1]' = (x²+1)(2x) = 2x³ + 2x
        //   G[0]' G[1] = 2x(x²+3)   = 2x³ + 6x
        //   H[2] = (2x³+2x) − (2x³+6x) = −4x
        let expected_h2 = linr(0, -4, &fp2);
        assert_eq!(h[2], expected_h2);
    }

    #[test]
    fn codomain_is_degree_5_for_mixed_split() {
        let fp2 = ctx();
        let (_c, s) = split_curve(&fp2);
        let c_tilde = codomain(&s, &fp2);
        // (quad × quad × linear) → H₀ deg 2, H₁ deg 2, H₂ has the cancellation
        // that drops it to linear since G_0, G_1 are both monic quadratics with
        // different x-coefficient ⇒ H[2] is genuinely deg 2 here. So H[0]H[1]H[2]
        // ends up deg 5 (with the linear G_2 making one of the H's linear).
        //
        // For this particular split G_2 = x − 2 → H[0] uses G_2, H[1] uses G_2,
        // both ending up linear/quadratic. Total degree depends on cancellation;
        // empirically for this case it's 5.
        assert_eq!(c_tilde.f.degree(), Some(5));
    }

    /// Count #C(F_p) for a deg-5 curve y² = f(x) with f having F_p coefficients.
    fn count_points_fp(c: &Curve, fp2: &Fp2) -> u64 {
        let p_u: u64 = fp2.fp.p.to_string().parse().unwrap();
        let mut n: u64 = 1; // ∞
        for x in 0..p_u {
            let xe = fp2.from_int(x as i64);
            let rhs = c.f.eval(&xe, fp2);
            assert!(rhs.b.is_zero());
            if rhs.a.is_zero() { n += 1; }
            else if fp2.fp.is_square(&rhs.a) { n += 2; }
        }
        n
    }

    /// Count #C(F_{p²}) — used to derive #J(F_p) via the L-polynomial.
    fn count_points_fp2(c: &Curve, fp2: &Fp2) -> u64 {
        let p_u: u64 = fp2.fp.p.to_string().parse().unwrap();
        let mut n: u64 = 1;
        let exp = (&fp2.fp.p * &fp2.fp.p - BigInt::from(1)) / 2;
        for a in 0..p_u {
            for b in 0..p_u {
                let x = F2 { a: BigInt::from(a), b: BigInt::from(b) };
                let rhs = c.f.eval(&x, fp2);
                if fp2.is_zero(&rhs) { n += 1; }
                else if fp2.pow(&rhs, &exp) == fp2.one() { n += 2; }
            }
        }
        n
    }

    fn jac_order(c: &Curve, fp2: &Fp2) -> i128 {
        let n1 = count_points_fp(c, fp2) as i128;
        let n2 = count_points_fp2(c, fp2) as i128;
        let p = 431i128;
        let a1 = n1 - (p + 1);
        let a2 = (a1 * a1 + n2 - p * p - 1) / 2;
        p * p + 1 + (p + 1) * a1 + a2
    }

    #[test]
    fn jacobian_order_preserved_under_isogeny() {
        let fp2 = ctx();
        let (c, s) = split_curve(&fp2);
        let c_tilde = codomain(&s, &fp2);
        let j_c = jac_order(&c, &fp2);
        let j_ct = jac_order(&c_tilde, &fp2);
        println!("|J(C)| = {j_c},  |J(C̃)| = {j_ct}");
        assert_eq!(j_c, j_ct, "(2,2)-isogeny must preserve Jacobian order");
    }

    /// Build a divisor on C from a single affine point (x₀, y₀) on C.
    fn pt_to_div(x0: &F2, y0: &F2, fp2: &Fp2) -> Div {
        let u = Poly::new(vec![fp2.neg(x0), fp2.one()], fp2);
        let v = Poly::constant(y0.clone(), fp2);
        Div { u, v }
    }

    fn sample_pt(c: &Curve, x_try: i64, fp2: &Fp2) -> Option<(F2, F2)> {
        let x = fp2.from_int(x_try);
        let rhs = c.f.eval(&x, fp2);
        if !rhs.b.is_zero() { return None; }
        // Skip 2-torsion points (y = 0) — they cause v* to vanish identically
        // in the pushforward, exercising a degenerate code path.
        if rhs.a.is_zero() { return None; }
        let y_re = fp2.fp.sqrt(&rhs.a)?;
        Some((x, F2 { a: y_re, b: num_bigint::BigInt::from(0) }))
    }

    #[test]
    fn pushforward_yields_valid_mumford_form() {
        let fp2 = ctx();
        let (c, s) = split_curve(&fp2);
        // Build a deg-2 divisor by combining two affine points on C.
        let mut pts = vec![];
        for k in 0..200i64 {
            if let Some(p) = sample_pt(&c, k, &fp2) {
                pts.push(p);
                if pts.len() == 2 { break; }
            }
        }
        let d1 = pt_to_div(&pts[0].0, &pts[0].1, &fp2);
        let d2 = pt_to_div(&pts[1].0, &pts[1].1, &fp2);
        let d = jacobian::add(&c, &d1, &d2, &fp2);
        assert_eq!(d.u.degree(), Some(2), "test divisor should have deg-2 u");
        assert!(d.is_valid(&c, &fp2), "test divisor must be valid Mumford on C");

        let (c_tilde, d_image) = pushforward(&s, &d, &fp2);
        // The image must satisfy Mumford invariants on the codomain.
        assert!(
            d_image.is_valid(&c_tilde, &fp2),
            "image divisor must satisfy v² ≡ f̃ (mod u) on codomain.\n\
             u_image = {:?}\n  v_image = {:?}",
            d_image.u, d_image.v
        );
    }

    #[test]
    fn pushforward_respects_doubling() {
        // φ([2]D) = [2] φ(D). This is the key linearity test.
        let fp2 = ctx();
        let (c, s) = split_curve(&fp2);
        let mut pts = vec![];
        for k in 0..200i64 {
            if let Some(p) = sample_pt(&c, k, &fp2) {
                pts.push(p);
                if pts.len() == 3 { break; }
            }
        }
        let d1 = pt_to_div(&pts[0].0, &pts[0].1, &fp2);
        let d2 = pt_to_div(&pts[1].0, &pts[1].1, &fp2);
        let d = jacobian::add(&c, &d1, &d2, &fp2);

        let (c_tilde, phi_d) = pushforward(&s, &d, &fp2);
        let two_d = jacobian::add(&c, &d, &d, &fp2);
        let (c_tilde2, phi_two_d) = pushforward(&s, &two_d, &fp2);
        let two_phi_d = jacobian::add(&c_tilde, &phi_d, &phi_d, &fp2);

        assert_eq!(c_tilde.f, c_tilde2.f, "codomain must be deterministic");
        assert_eq!(
            phi_two_d, two_phi_d,
            "φ should be a group homomorphism: φ([2]D) = [2]φ(D)"
        );
    }

    #[test]
    fn chain_a2_split_partition_detected() {
        // With a=2, the chain executor's intermediate loop runs 0 times and
        // it just performs the final splitting check. Construct two 2-torsion
        // divisors whose u-polys are linearly dependent with the third
        // factor's coefficients — i.e., a split partition (δ=0).
        let fp2 = ctx();
        let (c, _s) = split_curve(&fp2);
        // We need 2-torsion divisors d1, d2 on J(c) whose u-polys G_1, G_2
        // give δ=0 with G_3 = c.f/(G_1·G_2). Use a known split partition:
        //   G_1 = x²+1, G_2 = x²+3, G_3 = (linear)*(linear)
        //   Their coefficient matrix has middle column [0, 0, ...] making δ=0.
        // For our test curve f = (x²+1)(x²+3)(x-2):
        //   G_1 = x²+1, G_2 = x²+3 (both have x-coef 0)
        //   G_3 = x-2 (only 1 quadratic among the 3 factors though — won't have δ=0)
        // Better: use the curve f = (x²+1)(x²+2)(2-x²-1) ... too complex.
        //
        // Use a simpler construction: build a fresh curve where the natural
        // splitting IS δ=0 by construction (linearly dependent rows).
        let g0 = Poly::new(vec![fp2.from_int(1), fp2.from_int(0), fp2.from_int(1)], &fp2); // x²+1
        let g1 = Poly::new(vec![fp2.from_int(2), fp2.from_int(0), fp2.from_int(1)], &fp2); // x²+2
        // g2 = g0+g1 = 2x²+3 (gives δ=0 in 3x3 cofactor matrix with first 2 rows
        // having the same x-coef pattern). Need product = monic deg 6.
        let g2_unnorm = Poly::new(vec![fp2.from_int(3), fp2.from_int(0), fp2.from_int(2)], &fp2);
        let f_split = g0.mul(&g1, &fp2).mul(&g2_unnorm, &fp2);
        let curve_with_split = jacobian::Curve::new(f_split, &fp2);

        // 2-torsion divisor on curve_with_split with u = g0:
        // u = x²+1, v = 0. v² = 0, f mod u = ?  (f is divisible by g0, so f mod g0 = 0)
        let d_a = Div { u: g0.clone(), v: Poly::zero() };
        let d_b = Div { u: g1.clone(), v: Poly::zero() };
        assert!(d_a.is_valid(&curve_with_split, &fp2));
        assert!(d_b.is_valid(&curve_with_split, &fp2));

        // a=2 chain: just final splitting check on (G_1, G_2, G_3).
        let split = does_22_chain_split(&curve_with_split, &d_a, &d_b, 2, &fp2);
        assert!(split, "constructed split partition should yield does_22_chain_split = true");
    }

    #[test]
    fn chain_a2_non_split_partition_detected() {
        let fp2 = ctx();
        let (c, _s) = split_curve(&fp2);
        // Use the test split's natural G_1, G_2 — these don't give δ=0.
        // G_1 = x²+1, G_2 = x²+3 (the curve's natural factors). δ ≠ 0 here.
        let g0 = Poly::new(vec![fp2.from_int(1), fp2.from_int(0), fp2.from_int(1)], &fp2);
        let g1 = Poly::new(vec![fp2.from_int(3), fp2.from_int(0), fp2.from_int(1)], &fp2);
        let d_a = Div { u: g0.clone(), v: Poly::zero() };
        let d_b = Div { u: g1.clone(), v: Poly::zero() };
        assert!(d_a.is_valid(&c, &fp2));
        assert!(d_b.is_valid(&c, &fp2));
        let split = does_22_chain_split(&c, &d_a, &d_b, 2, &fp2);
        assert!(!split, "non-split partition should yield does_22_chain_split = false");
    }

    #[test]
    fn pushforward_kills_kernel_element() {
        // A 2-torsion divisor whose u-poly equals G_1 (one of the kernel
        // generators) should be killed by the Richelot isogeny.
        let fp2 = ctx();
        let (_c, s) = split_curve(&fp2);
        // G_1 = x² + 1.  D = ⟨G_1, 0⟩ is a 2-torsion divisor in ker(φ).
        let g1 = s.g[0].clone();
        let kernel_div = Div { u: g1, v: Poly::zero() };
        // u has degree 2 — meets our pushforward precondition.
        let (_, phi_d) = pushforward(&s, &kernel_div, &fp2);
        assert!(
            phi_d.is_identity(&fp2),
            "kernel element should map to identity, got u={:?} v={:?}",
            phi_d.u, phi_d.v
        );
    }
}
