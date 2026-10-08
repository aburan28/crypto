//! Radical isogenies (Castryck–Decru–Vercauteren, ASIACRYPT 2020): chains of N-isogenies where
//! each step costs one N-th root and a few field operations, with no torsion-point search.
//!
//!  N = 3: E: y^2 + a1 xy + a3 y = x^3 (P = (0,0) a flex point), rho = -a3, alpha = rho^(1/3):
//!         a1' = -6 alpha + a1,  a3' = 3 a1 alpha^2 - a1^2 alpha + 9 a3.
//!  N = 5: Tate normal form y^2 + (1-b) xy - b y = x^3 - b x^2, rho = b, alpha = rho^(1/5):
//!         b' = alpha (alpha^4 + 3 alpha^3 + 4 alpha^2 + 2 alpha + 1)
//!                   / (alpha^4 - 2 alpha^3 + 4 alpha^2 - 3 alpha + 1).
//! Each step continues with the P-distinguished point (non-backtracking); over F_q with
//! gcd(N, q - 1) = 1 the N-th root is unique, alpha = rho^(N^-1 mod q-1), and the chain stays
//! F_q-rational. On a supersingular curve over F_p (CSIDH) the distinguished point is in E(F_p),
//! so the chain is the action of the same ideal class l_N^k.
use crate::bigint::Big;
use crate::curve::{Curve, Pt};
use crate::field::Field;
use crate::int::Int;

/// j-invariant of y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6.
pub fn weierstrass_j<F: Field>(f: &F, a: [F::E; 5]) -> F::E {
    let [a1, a2, a3, a4, a6] = a;
    let c = |v: u64| f.from_u64(v);
    let b2 = f.add(f.mul(a1, a1), f.mul(c(4), a2));
    let b4 = f.add(f.mul(c(2), a4), f.mul(a1, a3));
    let b6 = f.add(f.mul(a3, a3), f.mul(c(4), a6));
    let c4 = f.sub(f.mul(b2, b2), f.mul(c(24), b4));
    let c6 = f.add(f.sub(f.mul(c(36), f.mul(b2, b4)), f.mul(b2, f.mul(b2, b2))), f.neg(f.mul(c(216), b6)));
    let c43 = f.mul(c4, f.mul(c4, c4));
    let disc1728 = f.sub(c43, f.mul(c6, c6)); // 1728 Delta
    f.div(f.mul(c(1728), c43), disc1728)
}

/// Short Weierstrass model of a general Weierstrass curve: y^2 = x^3 - 27 c4 x - 54 c6.
pub fn to_short<F: Field>(f: &F, a: [F::E; 5]) -> Curve<F::E> {
    let [a1, a2, a3, a4, a6] = a;
    let c = |v: u64| f.from_u64(v);
    let b2 = f.add(f.mul(a1, a1), f.mul(c(4), a2));
    let b4 = f.add(f.mul(c(2), a4), f.mul(a1, a3));
    let b6 = f.add(f.mul(a3, a3), f.mul(c(4), a6));
    let c4 = f.sub(f.mul(b2, b2), f.mul(c(24), b4));
    let c6 = f.add(f.sub(f.mul(c(36), f.mul(b2, b4)), f.mul(b2, f.mul(b2, b2))), f.neg(f.mul(c(216), b6)));
    Curve::new(f.neg(f.mul(c(27), c4)), f.neg(f.mul(c(54), c6)))
}

/// From y^2 = x^3 + A x + B and a point P of order N >= 3: move P to (0,0) with horizontal
/// tangent, giving y^2 + a1 xy + a3 y = x^3 + a2 x^2 (a2 = 0 exactly when N = 3).
pub fn tangent_form<F: Field>(f: &F, e: &Curve<F::E>, p: &Pt<F::E>) -> Option<(F::E, F::E, F::E)> {
    let Pt::Aff(x0, y0) = *p else { return None };
    if f.is_zero(y0) {
        return None;
    }
    let c = |v: u64| f.from_u64(v);
    let lam = f.div(f.add(f.mul(c(3), f.mul(x0, x0)), e.a), f.add(y0, y0));
    let a1 = f.add(lam, lam);
    let a3 = f.add(y0, y0);
    let a2 = f.sub(f.mul(c(3), x0), f.mul(lam, lam));
    Some((a1, a2, a3))
}

/// Tate normal form parameters (b, c) of y^2 + a1 xy + a3 y = x^3 + a2 x^2 (a2, a3 != 0):
/// b = -a2^3/a3^2, c = 1 - a1 a2/a3.
pub fn tate_bc<F: Field>(f: &F, a1: F::E, a2: F::E, a3: F::E) -> (F::E, F::E) {
    let b = f.neg(f.div(f.mul(a2, f.mul(a2, a2)), f.mul(a3, a3)));
    let cc = f.sub(f.one(), f.div(f.mul(a1, a2), a3));
    (b, cc)
}

/// The exponent N^-1 mod (q - 1) (None if gcd(N, q - 1) != 1): x -> x^e is the unique N-th root.
pub fn root_exponent<F: Field>(f: &F, n: u64) -> Option<Big> {
    let qm1 = Int::from_big(&f.q().sub_small(1));
    let e = Int::from(n).inv_mod(&qm1)?;
    Some(e.mag().clone())
}

/// One radical 3-isogeny step on (a1, a3).
pub fn step3<F: Field>(f: &F, a1: F::E, a3: F::E, e3: &Big) -> (F::E, F::E) {
    let alpha = f.pow_big(f.neg(a3), e3);
    let c = |v: u64| f.from_u64(v);
    let a1n = f.sub(a1, f.mul(c(6), alpha));
    let al2 = f.mul(alpha, alpha);
    let a3n = f.add(f.sub(f.mul(c(3), f.mul(a1, al2)), f.mul(f.mul(a1, a1), alpha)), f.mul(c(9), a3));
    (a1n, a3n)
}

/// One radical 5-isogeny step on b (Tate normal form with c = b).
pub fn step5<F: Field>(f: &F, b: F::E, e5: &Big) -> F::E {
    let a = f.pow_big(b, e5);
    let c = |v: u64| f.from_u64(v);
    let a2 = f.mul(a, a);
    let a3 = f.mul(a2, a);
    let a4 = f.mul(a3, a);
    let num = f.add(f.add(f.add(f.add(a4, f.mul(c(3), a3)), f.mul(c(4), a2)), f.mul(c(2), a)), f.one());
    let den = f.add(f.sub(f.add(f.sub(a4, f.mul(c(2), a3)), f.mul(c(4), a2)), f.mul(c(3), a)), f.one());
    f.div(f.mul(a, num), den)
}

/// Evaluate an integer-coefficient polynomial (low->high degree) at t over F.
fn ipoly<F: Field>(f: &F, coeffs: &[i64], t: F::E) -> F::E {
    let mut acc = f.zero();
    for &c in coeffs.iter().rev() {
        let cc = if c >= 0 { f.from_u64(c as u64) } else { f.neg(f.from_u64((-c) as u64)) };
        acc = f.add(f.mul(acc, t), cc);
    }
    acc
}

/// X_1(7) family in Tate normal form: parameter t gives b = t^3 - t^2, c = t^2 - t
/// (the order-7 locus b^2 - b c - c^3 = 0, parametrised through the node). P = (0,0) has order 7.
pub fn family7<F: Field>(f: &F, t: F::E) -> (F::E, F::E) {
    let t2 = f.mul(t, t);
    let t3 = f.mul(t2, t);
    (f.sub(t3, t2), f.sub(t2, t)) // (b, c)
}

/// One radical 7-isogeny step on the X_1(7) parameter t. The radicand is rho = t (t-1)^2 and
/// alpha = rho^(1/7) (unique when gcd(7, q-1) = 1, via e7 = 7^-1 mod q-1); then
///   t' = ( P0(t) + P1 alpha + ... + P6 alpha^6 ) / D(t),
/// with the integer polynomials below (derived in this work: radicand from the discriminant
/// ramification and the F_7 split-locus, update map by rational reconstruction, cross-checked
/// against the degree-7 modular correspondence and the Velu 7-isogeny orbit). The chain follows
/// the P-distinguished direction and is non-backtracking.
pub fn step7<F: Field>(f: &F, t: F::E, e7: &Big) -> F::E {
    // rho = t (t-1)^2
    let tm1 = f.sub(t, f.one());
    let rho = f.mul(t, f.mul(tm1, tm1));
    let a = f.pow_big(rho, e7); // alpha
    // numerator coefficients P_k(t) (cleared by *7); denominator D(t) likewise.
    const D: [i64; 5] = [1, 4, -13, 9, -1];
    const P: [&[i64]; 7] = [
        &[0, 12, -28, 19, -3],
        &[1, 0, -6, 5],
        &[2, -14, 16, -4],
        &[0, -7, 7],
        &[-5, 9, -3],
        &[-4, 10, -1],
        &[7],
    ];
    let mut num = f.zero();
    let mut ak = f.one(); // alpha^k
    for pk in P.iter() {
        num = f.add(num, f.mul(ipoly(f, pk, t), ak));
        ak = f.mul(ak, a);
    }
    f.div(num, ipoly(f, &D, t))
}

/// Tate-normal-form Weierstrass coefficients of the X_1(7) family member with parameter t.
pub fn model7<F: Field>(f: &F, t: F::E) -> [F::E; 5] {
    let (b, c) = family7(f, t);
    tate_model(f, b, c)
}

/// Weierstrass coefficients of the N = 3 model and of the Tate normal form.
pub fn model3<F: Field>(f: &F, a1: F::E, a3: F::E) -> [F::E; 5] {
    [a1, f.zero(), a3, f.zero(), f.zero()]
}
pub fn tate_model<F: Field>(f: &F, b: F::E, c: F::E) -> [F::E; 5] {
    [f.sub(f.one(), c), f.neg(b), f.neg(b), f.zero(), f.zero()]
}
