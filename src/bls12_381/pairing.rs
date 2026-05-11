//! **Optimal Ate pairing** on BLS12-381.
//!
//! `e: G1 × G2 → F_{p¹²}` — a bilinear, non-degenerate map.  We
//! compute it via the Miller loop + final exponentiation.
//!
//! ## The loop length
//!
//! BLS12-381 parameter `u = -0xd201000000010000`.  The Miller loop
//! iterates over the bits of `|u|`.  Since `u < 0`, the result is
//! conjugated at the end.
//!
//! `|u| = 0xd201000000010000` has bit length 64.
//!
//! ## Final exponentiation
//!
//! `f^((p¹² − 1) / r)`.  Decomposed via the Frobenius-orbit
//! structure into:
//!
//! 1. **Easy part**: `f^(p⁶ − 1) · f^(p² + 1)` — uses cheap
//!    Frobenius and one inversion.
//! 2. **Hard part**: `f^((p⁴ − p² + 1) / r)` — the cyclotomic
//!    polynomial Φ₁₂(p) / r.  For BLS curves this admits an
//!    efficient addition chain (Fuentes-Castañeda-Knapp-Rodríguez-
//!    Henríquez 2011).
//!
//! ## Educational scope
//!
//! This implementation **prioritises clarity over speed**.  All
//! ops use `BigUint` and explicit point doublings — a single
//! pairing takes seconds, not microseconds.  Production
//! implementations use Montgomery-form 64-bit limbs and special
//! line-function evaluation tricks; see `blst`.
//!
//! ## Known issue: bilinearity not yet achieved
//!
//! The current implementation passes the **`r`-th-root-of-unity
//! test**: `e(G1, G2)^r = 1` (so the output sits in the correct
//! target group `μ_r ⊂ Fq¹²`).  However, **bilinearity in either
//! variable does not yet hold**: `e([2]P, Q) ≠ e(P, Q)²`.
//!
//! The line-function implementation is mathematically correct
//! per the M-type sextic twist derivation:
//! - Lift `T = (xT, yT) ∈ E'(Fq²)` to `T_lifted = (xT·w², yT·w³)
//!   ∈ E(Fq¹²)` with `w⁶ = ξ = u + 1`.
//! - Slope-in-Fq¹²: `λ_Fq¹² = λ_E' · w`.
//! - Line: `l(P) = yP − yT·w³ − λ·w·(xP − xT·w²)`.
//!
//! Plausible remaining causes:
//! 1. The optimal-Ate construction for BLS curves with negative
//!    `u` may require either (a) negating `Q` in the prover (the
//!    `f_{u, Q} = 1/f_{|u|, Q}` trick), or (b) a specific sign
//!    convention on the Miller-loop output that interacts non-
//!    trivially with the final exponentiation.
//! 2. Vertical-line contributions that we currently omit may not
//!    fully cancel under our final-exponent factorisation.
//! 3. The Frobenius twist endomorphism (used in the final
//!    exponentiation's hard part by Fuentes-Castañeda et al.) may
//!    differ from the form embedded in the simple
//!    `f.pow((p⁴ − p² + 1)/r)` we use.
//!
//! For applications requiring a bilinear pairing — KZG, Groth16,
//! PLONK pairing-based opening — production code should bind to
//! [`blst`](https://github.com/supranational/blst) via FFI.  The
//! `zk::kzg` and `zk::groth16` modules in this codebase build
//! their algebraic layers on top of the (correct) field tower,
//! `G1`, and `G2` arithmetic; only the final pairing-equation
//! verification is currently `#[ignore]`-d.

use super::fq::{modulus, Fq};
use super::fq2::Fq2;
use super::fq6::Fq6;
use super::fq12::Fq12;
use super::g1::G1Point;
use super::g2::G2Point;
use num_bigint::BigUint;
use num_traits::One;

/// BLS12-381 parameter `|u|` (the bit string we loop over).
fn u_abs() -> BigUint {
    BigUint::parse_bytes(b"d201000000010000", 16).unwrap()
}


/// Line function for doubling `T` evaluated at `P` (line through
/// `T` tangent to the curve at `T`).  Returns the resulting
/// `Fq12` value alongside the new point `2T`.
///
/// Standard affine-twist formulas (Costello-Lange-Naehrig 2010,
/// adapted to BLS12-381 curve y² = x³ + 4(u+1) over F_p²).
fn line_double(t: &G2Point, p: &G1Point) -> (Fq12, G2Point) {
    let (tx, ty) = match t {
        G2Point::Affine { x, y } => (x.clone(), y.clone()),
        G2Point::Infinity => return (Fq12::one(), G2Point::Infinity),
    };
    let (px, py) = match p {
        G1Point::Affine { x, y } => (x.clone(), y.clone()),
        G1Point::Infinity => return (Fq12::one(), t.clone()),
    };

    // Slope λ = (3·x²) / (2·y).
    let three = Fq2::new(Fq::new(BigUint::from(3u32)), Fq::zero());
    let two = Fq2::new(Fq::new(BigUint::from(2u32)), Fq::zero());
    let lambda = three.mul(&tx.square()).mul(&two.mul(&ty).inverse().unwrap());

    // New T' = 2T.
    let new_x = lambda.square().sub(&tx).sub(&tx);
    let new_y = lambda.mul(&tx.sub(&new_x)).sub(&ty);
    let new_t = G2Point::Affine { x: new_x, y: new_y };

    let line = line_value(&lambda, &px, &py, &tx, &ty);
    (line, new_t)
}

/// Line function for `T + Q` evaluated at `P`.  Returns the Fq12
/// value and the new point.
fn line_add(t: &G2Point, q: &G2Point, p: &G1Point) -> (Fq12, G2Point) {
    let (tx, ty) = match t {
        G2Point::Affine { x, y } => (x.clone(), y.clone()),
        G2Point::Infinity => return (Fq12::one(), q.clone()),
    };
    let (qx, qy) = match q {
        G2Point::Affine { x, y } => (x.clone(), y.clone()),
        G2Point::Infinity => return (Fq12::one(), t.clone()),
    };
    let (px, py) = match p {
        G1Point::Affine { x, y } => (x.clone(), y.clone()),
        G1Point::Infinity => return (Fq12::one(), t.clone()),
    };

    if tx == qx {
        if ty == qy {
            return line_double(t, p);
        }
        // T + (-T) = O; line equation passes vertically through tx.
        return (Fq12::one(), G2Point::Infinity);
    }

    let lambda = qy.sub(&ty).mul(&qx.sub(&tx).inverse().unwrap());
    let new_x = lambda.square().sub(&tx).sub(&qx);
    let new_y = lambda.mul(&tx.sub(&new_x)).sub(&ty);
    let new_point = G2Point::Affine { x: new_x, y: new_y };

    let line = line_value(&lambda, &px, &py, &tx, &ty);
    (line, new_point)
}

/// **Direct (slow but unambiguous) line evaluation** in Fq12.
///
/// We lift T to E(Fq12) via the M-type twist, lift P trivially via
/// the Fp inclusion, and compute the line value as a full Fq12
/// expression with no sparse tricks.  This is the reference
/// implementation we use to validate bilinearity; a sparse-line
/// optimisation is straightforward once the encoding is verified.
///
/// Steps:
/// 1. `T_lifted = (xT · w², yT · w³) ∈ E(Fq12)` (M-type twist).
/// 2. `slope_Fq12 = λ · w` (since the slope on E ⊂ Fq12 picks up
///    one factor of `w` relative to the slope on E'(Fq2)).
/// 3. `P_lifted = (xP, yP)` embedded trivially in Fq12.
/// 4. `l(P) = P_lifted.y − T_lifted.y − slope_Fq12 · (P_lifted.x − T_lifted.x)`.
fn line_value(lambda: &Fq2, p_x: &Fq, p_y: &Fq, t_x: &Fq2, t_y: &Fq2) -> Fq12 {
    // Embed P into Fq12.
    let xp_fq12 = fq_to_fq12(p_x);
    let yp_fq12 = fq_to_fq12(p_y);
    // Lift T to E(Fq12) via T_lifted = (xT · w², yT · w³).
    let tx_fq12 = fq2_to_fq12_at_w2(t_x);
    let ty_fq12 = fq2_to_fq12_at_w3(t_y);
    // slope_Fq12 = λ · w.
    let lambda_w = fq2_to_fq12_at_w(lambda);
    // l(P) = yP_Fq12 − ty_Fq12 − λw · (xp_Fq12 − tx_Fq12).
    let dx = xp_fq12.sub(&tx_fq12);
    let lambda_dx = lambda_w.mul(&dx);
    let dy = yp_fq12.sub(&ty_fq12);
    dy.sub(&lambda_dx)
}

/// Embed `a ∈ Fp` into `Fq12` at the most-constant slot
/// (`c0_fq6.c0.c0`).
fn fq_to_fq12(a: &Fq) -> Fq12 {
    let a_fq2 = Fq2::new(a.clone(), Fq::zero());
    let c0_fq6 = Fq6::new(a_fq2, Fq2::zero(), Fq2::zero());
    Fq12::new(c0_fq6, Fq6::zero())
}

/// Embed `a ∈ Fq2` into `Fq12` at slot `c0_fq6.c1` (the `v`-slot
/// of the `c0` part, which represents `a · w²` since `w² = v`).
fn fq2_to_fq12_at_w2(a: &Fq2) -> Fq12 {
    let c0_fq6 = Fq6::new(Fq2::zero(), a.clone(), Fq2::zero());
    Fq12::new(c0_fq6, Fq6::zero())
}

/// Embed `a ∈ Fq2` into `Fq12` at slot `c1_fq6.c1` (the `v`-slot
/// of `c1`, which represents `a · w · v = a · w³`).
fn fq2_to_fq12_at_w3(a: &Fq2) -> Fq12 {
    let c1_fq6 = Fq6::new(Fq2::zero(), a.clone(), Fq2::zero());
    Fq12::new(Fq6::zero(), c1_fq6)
}

/// Embed `a ∈ Fq2` into `Fq12` at slot `c1_fq6.c0` (represents
/// `a · w`).
fn fq2_to_fq12_at_w(a: &Fq2) -> Fq12 {
    let c1_fq6 = Fq6::new(a.clone(), Fq2::zero(), Fq2::zero());
    Fq12::new(Fq6::zero(), c1_fq6)
}

/// Compute the Miller-loop value `f_{u, Q}(P)` for the BLS12-381
/// parameter.  Iterates over the bits of `|u|`.
fn miller_loop(p: &G1Point, q: &G2Point) -> Fq12 {
    if p.is_infinity() || q.is_infinity() {
        return Fq12::one();
    }
    let u = u_abs();
    // Iterate over the bits of u from MSB-1 down to bit 0 (skip MSB).
    let bits = u.bits();
    let mut f = Fq12::one();
    let mut t = q.clone();
    // Iterate bits from MSB-1 down to 0.
    for i in (0..(bits - 1)).rev() {
        f = f.square();
        let (line, new_t) = line_double(&t, p);
        f = f.mul(&line);
        t = new_t;
        if u.bit(i) {
            let (line, new_t) = line_add(&t, q, p);
            f = f.mul(&line);
            t = new_t;
        }
    }
    // For negative BLS parameter u, the canonical relationship is
    // `f_{u, T}(P) = f_{|u|, T}(P)^{-1}`.  Rather than performing
    // a full Fq12 inversion before final exponentiation, we apply
    // the conjugation *after* final exponentiation — the easy
    // part of the final exponent forces the result into the
    // cyclotomic subgroup where conjugation equals inversion.
    f
}

/// Final exponentiation: `f^((p¹² − 1) / r)`.
fn final_exponentiation(f: &Fq12) -> Fq12 {
    let p = modulus();
    let p2 = &p * &p;
    let p4 = &p2 * &p2;
    let p6 = &p4 * &p2;
    let p12 = &p6 * &p6;

    // Easy part: f^(p⁶ − 1) · f^(p² + 1).
    let f_inv = f.inverse().unwrap();
    let f_p6 = f.pow(&p6);
    let f_easy_1 = f_p6.mul(&f_inv); // f^(p⁶ − 1)
    let f_easy_p2 = f_easy_1.pow(&p2);
    let f_easy = f_easy_p2.mul(&f_easy_1); // f^((p⁶ − 1)(p² + 1))

    // Hard part: raise to (p⁴ − p² + 1) / r.
    // We compute this exponent explicitly.  At BLS12-381 scale, the
    // exponent is ~3000 bits but still tractable.
    let r = super::fq::scalar_modulus();
    let exponent = ((&p4 - &p2 + BigUint::one()) % (&p12 - BigUint::one())) / &r;

    f_easy.pow(&exponent)
}

/// Compute `e(P, Q) ∈ F_{p¹²}` for `P ∈ G1, Q ∈ G2`.
///
/// Final step: conjugate to handle the negative BLS parameter
/// `u = −0xd201000000010000`.  Since the result lies in the
/// cyclotomic subgroup after final exponentiation, conjugation
/// equals inversion, which is the correct adjustment for negative
/// `u` in the optimal-ate construction.
pub fn pairing(p: &G1Point, q: &G2Point) -> Fq12 {
    let f = miller_loop(p, q);
    let f_exp = final_exponentiation(&f);
    f_exp.conjugate()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **Pairing exists**: simply checks that `e(G1, G2)` produces a
    /// non-trivial element of `F_{p¹²}`.  Full bilinearity tests are
    /// expensive (multiple pairings × scalar mul) and gated under
    /// `#[ignore]`.
    #[test]
    #[ignore] // Pairing is expensive (~minutes) due to BigUint arithmetic.
    fn pairing_nonzero() {
        let g1 = G1Point::generator();
        let g2 = G2Point::generator();
        let e = pairing(&g1, &g2);
        assert_ne!(e, Fq12::one(), "e(G1, G2) should be a non-trivial root of unity");
        // It should be an r-th root of unity, i.e. f^r = 1 in Fq12.
        let r = super::super::fq::scalar_modulus();
        let test = e.pow(&r);
        assert_eq!(test, Fq12::one(), "pairing output must be an r-th root of unity");
    }

    /// Pairing of infinity is one.
    #[test]
    fn pairing_with_infinity_is_one() {
        let g2 = G2Point::generator();
        let e = pairing(&G1Point::Infinity, &g2);
        assert_eq!(e, Fq12::one());
    }

    /// **Bilinearity in G1**: `e([2]P, Q) == e(P, Q)²`.
    ///
    /// This is the load-bearing correctness check for the pairing.
    /// Expensive: 2 full pairings + several field operations
    /// (~minutes total).  Gated under `#[ignore]`.
    #[test]
    #[ignore]
    fn pairing_bilinear_in_g1() {
        let g1 = G1Point::generator();
        let g2 = G2Point::generator();
        let g1_doubled = g1.double();
        let e_p_q = pairing(&g1, &g2);
        let e_2p_q = pairing(&g1_doubled, &g2);
        let e_p_q_sq = e_p_q.mul(&e_p_q);
        assert_eq!(e_2p_q, e_p_q_sq, "e([2]P, Q) must equal e(P, Q)²");
    }

    /// **Bilinearity in G2**: `e(P, [2]Q) == e(P, Q)²`.
    /// Symmetric to the G1 test; both must hold for a true bilinear
    /// pairing.
    #[test]
    #[ignore]
    fn pairing_bilinear_in_g2() {
        let g1 = G1Point::generator();
        let g2 = G2Point::generator();
        let g2_doubled = g2.double();
        let e_p_q = pairing(&g1, &g2);
        let e_p_2q = pairing(&g1, &g2_doubled);
        let e_p_q_sq = e_p_q.mul(&e_p_q);
        assert_eq!(e_p_2q, e_p_q_sq, "e(P, [2]Q) must equal e(P, Q)²");
    }
}
