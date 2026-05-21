// Semaev-style instance generators for index-calculus ECDLP.
//
// Current status:
//   - cyclic-n: working (Phase 0 benchmark for the RL pipeline)
//   - semaev_s3 / weil_descent: scaffolded, not implemented
//
// The math is summarized at the bottom of this file so a future implementor
// can fill in `semaev_s3_system` without re-reading three papers.

use crate::field::Fp;
use crate::monomial::Monomial;
use crate::poly::{Poly, Term};

// Cyclic-n: classical Gröbner benchmark family. Use this as the Phase-0 instance
// distribution to validate the RL pipeline against published baselines (sugar etc.)
// before moving to actual Semaev systems.
//
//   f_i = sum_{j=0}^{n-1} prod_{k=0}^{i-1} x_{(j+k) mod n}     for i = 1..n
//   f_n = x_0 * x_1 * ... * x_{n-1} - 1
pub fn cyclic(n: usize) -> Vec<Poly> {
    let mut polys = Vec::with_capacity(n);
    for i in 1..n {
        let mut terms: Vec<Term> = Vec::with_capacity(n);
        for j in 0..n {
            let mut exp = vec![0u32; n];
            for k in 0..i {
                exp[(j + k) % n] += 1;
            }
            terms.push(Term { coef: Fp::one(), mono: Monomial { exp } });
        }
        polys.push(Poly::from_terms(terms, n));
    }
    let prod = Term { coef: Fp::one(), mono: Monomial { exp: vec![1u32; n] } };
    let minus_one = Term { coef: -Fp::one(), mono: Monomial::one(n) };
    polys.push(Poly::from_terms(vec![prod, minus_one], n));
    polys
}

// === Real Semaev systems — not implemented ============================

// Curve parameters for short Weierstrass y^2 = x^3 + a*x + b over the base field.
// Once we move beyond F_p with small p, switch Fp -> a generic field trait
// (or ark-ff backed).
pub struct Curve {
    pub a: Fp,
    pub b: Fp,
}

// === Semaev summation polynomials ====================================
//
// S_m(X_1, ..., X_m) vanishes iff there exist points P_i = (X_i, Y_i) on the
// short-Weierstrass curve y^2 = x^3 + a*x + b with P_1 + ... + P_m = O.
//
// Closed forms (Semaev 2004, characteristic ≠ 2, 3):
//   S_2(X1, X2) = X1 - X2
//   S_3(X1, X2, X3) derived from P1+P2 = -P3 by eliminating y1, y2 via
//     y_i^2 = x_i^3 + a x_i + b. With A := P1(X1) + P2(X2) - (X1-X2)^2 (X1+X2+X3)
//     where P_i(X) = X^3 + a X + b, we get
//         S_3 = A^2 - 4 * P1(X1) * P2(X2)
//     (Symmetric in X1, X2, X3 after expansion.)
//   S_{m+1}(X1, ..., X_{m+1}) = Res_Y( S_{m-k}(X1, ..., X_{m-k-1}, Y),
//                                       S_{k+2}(X_{m-k}, ..., X_{m+1}, Y) )
//   The recursion uses a resultant in Y; pick m-k = ceil(m/2) for best balance.

// Build a 3-variable monomial X_i^p.
fn xi_pow(i: usize, p: u32) -> Monomial {
    let mut e = vec![0u32; 3];
    e[i] = p;
    Monomial { exp: e }
}

fn p_of_x(curve: &Curve, i: usize) -> Poly {
    // P(X_i) = X_i^3 + a * X_i + b   as a 3-variable polynomial.
    let one = Monomial::one(3);
    Poly::from_terms(vec![
        Term { coef: Fp::one(),  mono: xi_pow(i, 3) },
        Term { coef: curve.a,    mono: xi_pow(i, 1) },
        Term { coef: curve.b,    mono: one          },
    ], 3)
}

pub fn semaev_s2(_curve: &Curve) -> Poly {
    // 3-variable representation for uniformity with S_3 (X3 unused).
    Poly::from_terms(vec![
        Term { coef:  Fp::one(),  mono: xi_pow(0, 1) },
        Term { coef: -Fp::one(),  mono: xi_pow(1, 1) },
    ], 3)
}

pub fn semaev_s3(curve: &Curve) -> Poly {
    let one = Fp::one();
    let p1 = p_of_x(curve, 0);
    let p2 = p_of_x(curve, 1);

    // (X1 - X2)
    let diff = Poly::from_terms(vec![
        Term { coef:  one, mono: xi_pow(0, 1) },
        Term { coef: -one, mono: xi_pow(1, 1) },
    ], 3);
    // (X1 + X2 + X3)
    let sum3 = Poly::from_terms(vec![
        Term { coef: one, mono: xi_pow(0, 1) },
        Term { coef: one, mono: xi_pow(1, 1) },
        Term { coef: one, mono: xi_pow(2, 1) },
    ], 3);

    // A = P1 + P2 - (X1 - X2)^2 * (X1 + X2 + X3)
    let diff_sq = diff.mul(&diff);
    let term3 = diff_sq.mul(&sum3);
    let a_poly = p1.add(&p2).sub(&term3);

    // S_3 = A^2 - 4 * P1 * P2
    let a_sq = a_poly.mul(&a_poly);
    let four_p1p2 = p1.mul(&p2).scale(Fp::new(4));
    a_sq.sub(&four_p1p2)
}

// S_4(X1, X2, X3, X4) via resultant elimination of a shared dummy Y:
//   S_4 = Res_Y( S_3(X1, X2, Y),  S_3(X3, X4, Y) )
// Returned in 5-variable ambient [X1, X2, X3, X4, Y] where Y always has
// exponent 0 in every output term — the Y axis is just a workspace.
pub fn semaev_s4(curve: &Curve) -> Poly {
    let s3 = semaev_s3(curve);  // 3-var: [X1, X2, X3]
    // Embed both copies into a 5-var ambient [X1, X2, X3, X4, Y].
    let s3_left  = s3.embed(5, &[0, 1, 4]);  // (X1, X2, Y)
    let s3_right = s3.embed(5, &[2, 3, 4]);  // (X3, X4, Y)
    crate::resultant::resultant(&s3_left, &s3_right, 4)
}

// Generic dispatch.
pub fn semaev_s(curve: &Curve, m: usize) -> Poly {
    match m {
        2 => semaev_s2(curve),
        3 => semaev_s3(curve),
        4 => semaev_s4(curve),
        _ => unimplemented!("S_m for m ≥ 5: chain the resultant recursion further"),
    }
}

// Build the polynomial system whose common zeros are 2-decompositions of a
// target R into the factor base F, i.e. R = P_1 + P_2 with x(P_i) a root of
// `factor_base_poly`.
//
// Constructed equations (over F_p, 2 variables [x_1, x_2]):
//   1. S_3(x_R, x_1, x_2) = 0      // -P_1 - P_2 has x-coord x_R
//   2. F(x_1) = 0                  // P_1 is in the factor base
//   3. F(x_2) = 0                  // P_2 is in the factor base
//
// `factor_base_poly` is a univariate polynomial (nvars = 1) whose roots define
// the factor-base x-coordinates V ⊂ F_p. Typical choice: F(x) = prod (x - v)
// for v ∈ V.
//
// This is the simplest m=2 case over the base field — no Weil descent needed.
// Generalizing to m≥3 needs higher Semaev polynomials (S_4 via resultant);
// generalizing to F_{q^n} needs Weil restriction (`weil_restrict`).
pub fn decomposition_system_2(curve: &Curve, x_r: Fp, factor_base_poly: &Poly) -> Vec<Poly> {
    assert_eq!(factor_base_poly.nvars, 1, "factor base poly must be univariate");

    // S_3 starts as 3-var (x_1, x_2, x_3); fix x_3 := x_R and drop that var.
    // After elimination, vars = [x_1, x_2] (2 vars).
    let s3 = semaev_s3(curve);
    let s3_specialized = s3.eliminate(2, x_r);

    // Lift F from 1-var to 2-var as F(x_1) and F(x_2).
    let f_x1 = factor_base_poly.embed(2, &[0]);
    let f_x2 = factor_base_poly.embed(2, &[1]);

    vec![s3_specialized, f_x1, f_x2]
}

// 3-decomposition system: R = P_1 + P_2 + P_3 with all P_i in factor base.
//
// S_4 lives in 5-var ambient [X1, X2, X3, X4, Y] with Y always at exp 0.
// We specialize X4 := x_R and drop X4 and Y to get a polynomial in (X1,X2,X3).
// Then we add three factor-base constraints F(X_i) = 0.
pub fn decomposition_system_3(curve: &Curve, x_r: Fp, factor_base_poly: &Poly) -> Vec<Poly> {
    assert_eq!(factor_base_poly.nvars, 1, "factor base poly must be univariate");

    // S_4 is 5-var: [X1, X2, X3, X4, Y]. Y has exp 0 — eliminate it (constant 0).
    let s4 = semaev_s4(curve);          // nvars=5
    let s4_no_y = s4.eliminate(4, Fp::zero());   // nvars=4
    // Now [X1, X2, X3, X4]. Specialize X4 := x_R.
    let s4_spec = s4_no_y.eliminate(3, x_r);     // nvars=3

    // Lift F into 3-var ambient: F(X1), F(X2), F(X3).
    let f_x1 = factor_base_poly.embed(3, &[0]);
    let f_x2 = factor_base_poly.embed(3, &[1]);
    let f_x3 = factor_base_poly.embed(3, &[2]);

    vec![s4_spec, f_x1, f_x2, f_x3]
}

// 3-decomposition instance with small fixed factor base (suitable for testing).
pub fn semaev_decomposition_3_instance(fb_size: usize) -> Vec<Poly> {
    let curve = Curve { a: Fp::new(1), b: Fp::new(1) };
    let fb_xs = first_n_curve_xs(&curve, fb_size, 0);
    let fb_poly = univariate_product(&fb_xs);
    let mut x_r = Fp::new(24990);
    while fb_xs.iter().any(|v| *v == x_r) {
        x_r = x_r + Fp::new(1);
    }
    decomposition_system_3(&curve, x_r, &fb_poly)
}

// === Default Semaev training-instance distribution ====================
//
// 2-decomposition systems where the factor base size determines the LM
// structure (and hence Buchberger step count). Varying x_R alone changes
// only constants, not the basis structure — observed in practice that all
// such instances yield identical step counts. Varying |F| does change the
// LM structure (degree of F polys is |F|), giving real training diversity.
//
// Layout:
//   - Curve: y² = x³ + x + 1 over F_p (same as decomp_demo / ecdlp_demo).
//   - Factor base for `semaev2-{n}`: first n x-values starting at 0 that
//     yield curve points (curve has affine points there).
//   - x_R: a fixed value (24990), shifted up if it collides with the FB.

pub fn semaev_decomposition_instance(fb_size: usize) -> Vec<Poly> {
    let curve = Curve { a: Fp::new(1), b: Fp::new(1) };
    let fb_xs = first_n_curve_xs(&curve, fb_size, 0);
    let fb_poly = univariate_product(&fb_xs);
    let mut x_r = Fp::new(24990);
    while fb_xs.iter().any(|v| *v == x_r) {
        x_r = x_r + Fp::new(1);
    }
    decomposition_system_2(&curve, x_r, &fb_poly)
}

// Scan x = start, start+1, ... and return the first n that yield curve points.
fn first_n_curve_xs(curve: &Curve, n: usize, start: u64) -> Vec<Fp> {
    let mut out: Vec<Fp> = Vec::new();
    let mut xi: u64 = start;
    let limit = start.saturating_add(50_000);
    while out.len() < n && xi < limit {
        let x = Fp(xi % crate::field::P);
        let rhs = x * x * x + curve.a * x + curve.b;
        if rhs.is_zero() || is_square_fp(rhs) {
            out.push(x);
        }
        xi += 1;
    }
    assert_eq!(out.len(), n, "couldn't find {} affine x's starting at {}", n, start);
    out
}

// === Diverse Semaev training-instance distribution ====================
//
// `semaev2d-{seed}` builds a decomposition system whose curve, factor base,
// and target are all deterministically derived from `seed`. Used to test
// whether RL policies overfit to a fixed factor-base structure or actually
// generalize across the Semaev problem family.
//
// Mapping seed → instance:
//   a, b   : seed-hashed Fp values; skip degenerate (4a³+27b² == 0) curves
//   |F|    : 6 + (seed mod 6), so {6, 7, 8, 9, 10, 11}
//   FB     : first |F| valid x-coords starting from offset = (seed*7) % 200
//   x_R    : seed-hashed Fp, skipped to ensure x_R ∉ FB
pub fn semaev_decomposition_diverse(seed: u64) -> Vec<Poly> {
    let p_i = crate::field::P as i64;
    // Try seed, seed+1, ... until we get a non-singular curve.
    let mut s = seed;
    loop {
        let a_val = (s.wrapping_mul(31).wrapping_add(7)) % crate::field::P;
        let b_val = (s.wrapping_mul(53).wrapping_add(11)) % crate::field::P;
        let a = Fp(a_val);
        let b = Fp(b_val);
        // discriminant 4a³ + 27b² ≠ 0 ⇔ non-singular over Fp (char ≠ 2, 3).
        let disc = Fp::new(4) * a * a * a + Fp::new(27) * b * b;
        if disc.is_zero() {
            s = s.wrapping_add(1);
            continue;
        }
        let curve = Curve { a, b };
        let fb_size = 6 + ((s % 6) as usize);
        let offset = (s.wrapping_mul(7)) % 200;
        let fb_xs = first_n_curve_xs(&curve, fb_size, offset);
        let fb_poly = univariate_product(&fb_xs);
        // Target x_R, skip until disjoint from FB.
        let mut t = (s.wrapping_mul(97).wrapping_add(13)) % crate::field::P;
        loop {
            let x_r = Fp(t);
            if !fb_xs.iter().any(|v| *v == x_r) {
                let _ = p_i;
                return decomposition_system_2(&curve, x_r, &fb_poly);
            }
            t = (t + 1) % crate::field::P;
        }
    }
}

// Re-export of field::is_qr under the old name so the existing callsites
// in this module keep compiling.
fn is_square_fp(a: Fp) -> bool {
    crate::field::is_qr(a)
}

// Univariate product F(x) = prod_i (x - v_i).
fn univariate_product(vs: &[Fp]) -> Poly {
    let mut acc = Poly::from_terms(
        vec![Term { coef: Fp::one(), mono: Monomial::one(1) }],
        1,
    );
    for v in vs {
        let factor = Poly::from_terms(
            vec![
                Term { coef:  Fp::one(), mono: Monomial { exp: vec![1] } },
                Term { coef: -(*v),      mono: Monomial { exp: vec![0] } },
            ],
            1,
        );
        acc = acc.mul(&factor);
    }
    acc
}

// === Weil descent — scaffold =======================================
//
// The Weil restriction maps a polynomial system over F_{q^n} into n× as many
// equations over F_q, by expanding each F_{q^n}-variable in an F_q-basis of a
// chosen subspace V ⊂ F_{q^n}. This is the central trick in the Diem /
// Faugère-Petit-Renault / Gaudry-Huot index-calculus attacks: a small system
// over the big field (m variables) becomes a much bigger but lower-degree
// system over the small field (m*n variables), where F4 / Buchberger have a
// real chance of solving.
//
// What this scaffold doesn't have yet:
//   - F_{q^n} arithmetic: our Fp is fixed to one prime. Need either ark-ff
//     extension fields or a hand-rolled struct ExtFp<P, N> over an irreducible
//     defining polynomial of degree n.
//   - A representation for the subspace V: a basis {b_1, ..., b_n} ⊂ F_{q^n}
//     given as columns of an n×n matrix over F_q.
//   - Polynomial substitution and re-collection: substitute
//        x_i = sum_k v_{i,k} * b_k       (v_{i,k} are new variables in F_q)
//     into each input polynomial, then split each resulting F_{q^n}-coefficient
//     into its n F_q-components against the basis B. Each input equation
//     becomes n equations.
//
// Implementation sketch (~250–400 lines once the field side is done):
//
//   1. Build a multiplication table for the F_{q^n} basis B:
//        b_i * b_j = sum_k T_{i,j,k} * b_k        (T_{i,j,k} ∈ F_q)
//      This is what lets us re-collect into F_q components after products.
//
//   2. For each input polynomial f(x_1,...,x_m) over F_{q^n}:
//      a. Substitute x_i = sum_k v_{i,k} * b_k. The output is a polynomial in
//         {v_{i,k}} whose coefficients live in F_{q^n}.
//      b. Reduce each coefficient via T into its n F_q components, generating
//         one F_q-polynomial per basis element.
//      c. Append all n F_q-polynomials to the output list.
//
//   3. Pass the resulting list to the existing Buchberger / F4 engine. The
//      RL strategy hook applies unchanged.
//
// References to follow:
//   - Faugère, Gaudry, Huot, Renault (2014), "Using symmetries in the index
//     calculus..." — the variable-elimination subtleties that turn Weil
//     descent from "theoretically possible" into "practically fast".
//   - Diem (2011), "On the discrete logarithm problem in elliptic curves".
//     The complexity analysis that motivates picking n carefully.
//   - Petit, Quisquater (2012) — concrete F_{2^n} attack pipeline.

/// Placeholder type for "a basis of V ⊂ F_{q^n} expressed as Fp coordinates".
/// Real impl will replace with a proper ExtField wrapper.
pub type SubfieldBasis = Vec<Vec<Fp>>;  // rows = basis elements, cols = F_q coords

pub fn weil_restrict(
    _system_over_extension: Vec<Poly>,
    _basis: SubfieldBasis,
) -> Vec<Poly> {
    unimplemented!("Weil restriction — see scaffold notes above for the 3-step recipe")
}

// References to read before filling these in:
//   - Semaev, "Summation polynomials and the discrete logarithm problem on
//     elliptic curves" (2004). Section 2 has the explicit S_3 formula and the
//     resultant recursion for S_m.
//   - Diem, "On the discrete logarithm problem in elliptic curves" (2011).
//     Asymptotic analysis of when index calculus is subexponential.
//   - Faugère, Gaudry, Huot, Renault, "Using symmetries in the index calculus
//     for elliptic curves discrete logarithm" (2014). Symmetry-aware S_m
//     simplifications and Weil descent variable orderings — biggest practical
//     speedups currently known.
//   - Petit, Quisquater, "On polynomial systems arising from a Weil descent"
//     (2012). The concrete pipeline this module is meant to reproduce.
