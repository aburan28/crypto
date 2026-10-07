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

/// Weierstrass coefficients of the N = 3 model and of the Tate normal form.
pub fn model3<F: Field>(f: &F, a1: F::E, a3: F::E) -> [F::E; 5] {
    [a1, f.zero(), a3, f.zero(), f.zero()]
}
pub fn tate_model<F: Field>(f: &F, b: F::E, c: F::E) -> [F::E; 5] {
    [f.sub(f.one(), c), f.neg(b), f.neg(b), f.zero(), f.zero()]
}
