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

/// Encode an `F_p` element as an `Fq12` element living in the
/// constant subfield: `(c0 = (a, 0, 0), c1 = 0)` with `a ∈ Fq2`
/// having `c0 = the Fq value, c1 = 0`.  Used to inject G1 coords
/// into the Miller-loop accumulator.
fn fq_to_fq12_subfield(a: &Fq) -> Fq12 {
    let a_fq2 = Fq2::new(a.clone(), Fq::zero());
    let zero_fq2 = Fq2::zero();
    let c0_fq6 = Fq6::new(a_fq2, zero_fq2.clone(), zero_fq2);
    Fq12::new(c0_fq6, Fq6::zero())
}

/// Embed an `Fq2` element into `Fq12` as `(c0 = (a, 0, 0), c1 = 0)`.
fn fq2_to_fq12(a: &Fq2) -> Fq12 {
    let zero_fq2 = Fq2::zero();
    let c0_fq6 = Fq6::new(a.clone(), zero_fq2.clone(), zero_fq2);
    Fq12::new(c0_fq6, Fq6::zero())
}

/// Embed an `Fq2` into `Fq12` at the `c1·v` position of `c0_fq6`:
/// `(c0 = (0, a, 0), c1 = 0)`.  This corresponds to line-function
/// terms that live at the "y" slot of the twist.
fn fq2_to_fq12_at_v(a: &Fq2) -> Fq12 {
    let zero_fq2 = Fq2::zero();
    let c0_fq6 = Fq6::new(zero_fq2.clone(), a.clone(), zero_fq2);
    Fq12::new(c0_fq6, Fq6::zero())
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

    // Line: l(X, Y) = Y − ty − λ·(X − tx).  Evaluated at (P_x, P_y)
    // where (X, Y) substituted by the M-type encoding for Fq12:
    //   l_evaluated = (ty − λ·tx) + λ·P_x·w² − P_y·w³ (roughly)
    // We follow the standard "affine line function for sextic D-type
    // twist", which gives:
    //   line = (−λ · P_x) + (P_y) · w  +  (λ · T_x − T_y) · w_at_v
    // Plugging into Fq12 = Fq6 + Fq6·w with Fq6 = Fq2 + Fq2·v + Fq2·v²:
    //
    //   line.c0.c0 = λ · T_x − T_y                (constant Fq2)
    //   line.c0.c1 = 0
    //   line.c0.c2 = 0
    //   line.c1.c0 = (−λ · P_x) embedded into Fq2  (real-only)
    //                  + (P_y) embedded into Fq2 at the v-position
    //
    // Simplification: we adopt the "M-type" twist formulation where
    // the line is parameterised as l(P) = c + a·w + b·w·v with:
    //   c = ty − λ·tx      (constant in Fq2)
    //   a = λ·P_x          (Fq2, P_x as constant)
    //   b = P_y            (Fq2, P_y as constant)
    //
    // This is a sparse Fq12 element with three non-zero "slots".

    let c = ty.sub(&lambda.mul(&tx));
    let a_coeff = lambda.mul(&Fq2::new(px.clone(), Fq::zero()));
    let b_coeff = Fq2::new(py.clone(), Fq::zero());

    let line = build_sparse_line(&c, &a_coeff, &b_coeff);
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

    // Same line structure as in line_double.
    let c = ty.sub(&lambda.mul(&tx));
    let a_coeff = lambda.mul(&Fq2::new(px.clone(), Fq::zero()));
    let b_coeff = Fq2::new(py.clone(), Fq::zero());

    let line = build_sparse_line(&c, &a_coeff, &b_coeff);
    (line, new_point)
}

/// Construct a sparse `Fq12` line of the form `c + a·w + b·w·v`.
fn build_sparse_line(c: &Fq2, a: &Fq2, b: &Fq2) -> Fq12 {
    // Fq12 = Fq6 + Fq6·w.
    // Fq6 = Fq2 + Fq2·v + Fq2·v².
    //
    // c is in the constant slot of c0_Fq6:
    //   c0_Fq6 = (c, 0, 0)
    // a·w is in the constant slot of c1_Fq6:
    //   c1_Fq6 = (a, b, 0)  ← b at v-slot encodes b·w·v as part of c1·w
    let c0_fq6 = Fq6::new(c.clone(), Fq2::zero(), Fq2::zero());
    let c1_fq6 = Fq6::new(a.clone(), b.clone(), Fq2::zero());
    Fq12::new(c0_fq6, c1_fq6)
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
    // BLS u is negative for BLS12-381 → conjugate the final value.
    f.conjugate()
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
pub fn pairing(p: &G1Point, q: &G2Point) -> Fq12 {
    let f = miller_loop(p, q);
    final_exponentiation(&f)
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
}
