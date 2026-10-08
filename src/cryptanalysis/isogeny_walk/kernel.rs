//! An `ℓ`-isogeny's kernel polynomial from a root of `Φ_ℓ(j, Y)` (Elkies),
//! and the verifier that certifies an edge from the kernel polynomial
//! alone.
//!
//! # Construction
//!
//! With `E4 = −a/3`, `E6 = −b/2` and `J = θj = −jE6/E4` (`θ = q d/dq`),
//! differentiating `Φ_ℓ(j(τ), j(ℓτ)) = 0` once gives `J̃ = θj(ℓτ)/ℓ =
//! −ΦX J / (ℓ ΦY)`, hence `E4' = J̃² / (j'(j' − 1728))`, `E6' = −J̃ E4'/j'`
//! and the normalised codomain `a' = −3ℓ⁴E4'`, `b' = −2ℓ⁶E6'`.
//! Differentiating twice and using Ramanujan's `θE4 = (E2E4 − E6)/3`,
//! `θE6 = (E2E6 − E4²)/2` isolates the weight-two form
//!
//! ```text
//! E2 − ℓE2' = 6 / (ℓ ΦY J̃) · [ ΦXX J² + 2ℓ ΦXY J J̃ + ℓ² ΦYY J̃²
//!             + ΦX (⅔ j E6²/E4² + ½ j E4) + ℓ² ΦY (⅔ j' E6'²/E4'² + ½ j' E4') ],
//! ```
//!
//! and the `q`-expansion of `℘` at the `ℓ`-torsion gives
//! `Σ_{Q ∈ K∖0} x(Q) = −ℓ(E2 − ℓE2')`.  The remaining power sums follow
//! from Vélu's `℘_{E'}(z) = ℘_E(z) + Σ_Q [℘_E(z + Q) − ℘_E(Q)]`: the
//! `z^{2k}` coefficients give `Σ_Q ℘^{(2k)}(Q)/(2k)! = c'_k − c_k`, and
//! `℘^{(2k)}` is a polynomial of degree `k + 1` in `℘`.
//!
//! # Verification
//!
//! [`verify_kernel`] never looks at `Φ_ℓ` or at how `h` was found.  It
//! accepts `h` only if `h` is squarefree of degree `(ℓ − 1)/2`, divides the
//! `ℓ`-division polynomial, and its roots are closed under `P ↦ [g]P` for a
//! generator `g` of `(Z/ℓ)^*/±1` — which makes them exactly the `x`-
//! coordinates of one cyclic subgroup of order `ℓ`, stable under Galois
//! because `h ∈ F_p[x]`.  Vélu's formulas then give the codomain, which
//! must match the recorded target up to an `F_p`-isomorphism.

use super::curve::{self, Model};
use super::field::{Fe, Field};
use super::modpoly::ModPoly;
use super::poly::{self, Poly};

/// An `ℓ`-isogeny `E → E'` with its kernel polynomial.
#[derive(Clone, Debug)]
pub struct Isogeny {
    pub ell: u64,
    /// Monic, degree `(ℓ − 1)/2`, in the source model.
    pub kernel: Poly,
    /// Vélu's (normalised) codomain.
    pub codomain: Model,
}

/// Why a construction or a certificate failed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EdgeError {
    SpecialJ,
    MultipleRoot,
    NotSquarefree,
    WrongDegree,
    NotTorsion,
    NotSubgroup,
    CodomainMismatch,
    NotIsomorphic,
}

/// `ψ`'s odd-index factors: `f_n = ψ_n` for odd `n` and `ψ_n / 2y` for even
/// `n`, as polynomials in `x` reduced mod `m`, `n = 0..=top`.
pub fn division_polys_mod(f: &Field, e: &Model, top: usize, m: &Poly) -> Vec<Poly> {
    let c = |v: i64| f.from_i64(v);
    let (a, b) = (e.a, e.b);
    let a2 = f.sqr(&a);
    let a3 = f.mul(&a2, &a);
    let b2 = f.sqr(&b);
    let ab = f.mul(&a, &b);
    // F = 4(x³ + ax + b) = (2y)².
    let big_f = poly::rem(
        f,
        &poly::trim(f, vec![f.mul(&c(4), &b), f.mul(&c(4), &a), f.zero(), c(4)]),
        m,
    );
    let f3 = poly::trim(
        f,
        vec![
            f.neg(&a2),
            f.mul(&c(12), &b),
            f.mul(&c(6), &a),
            f.zero(),
            c(3),
        ],
    );
    let f4 = poly::trim(
        f,
        vec![
            f.mul(&c(2), &f.sub(&f.neg(&a3), &f.mul(&c(8), &b2))),
            f.mul(&c(-8), &ab),
            f.mul(&c(-10), &a2),
            f.mul(&c(40), &b),
            f.mul(&c(10), &a),
            f.zero(),
            c(2),
        ],
    );
    let mut v: Vec<Poly> = vec![
        Vec::new(),
        poly::rem(f, &vec![f.one()], m),
        poly::rem(f, &vec![f.one()], m),
        poly::rem(f, &f3, m),
        poly::rem(f, &f4, m),
    ];
    let mm = |x: &Poly, y: &Poly| poly::mulmod(f, x, y, m);
    let cube = |x: &Poly| mm(&mm(x, x), x);
    let f2 = mm(&big_f, &big_f);
    for n in 5..=top {
        let k = n / 2;
        let next = if n % 2 == 1 {
            // n = 2k + 1.
            let (t1, t2) = (mm(&v[k + 2], &cube(&v[k])), mm(&v[k - 1], &cube(&v[k + 1])));
            if k % 2 == 0 {
                poly::sub(f, &mm(&f2, &t1), &t2)
            } else {
                poly::sub(f, &t1, &mm(&f2, &t2))
            }
        } else {
            // n = 2k.
            let t1 = mm(&v[k + 2], &mm(&v[k - 1], &v[k - 1]));
            let t2 = mm(&v[k - 2], &mm(&v[k + 1], &v[k + 1]));
            mm(&v[k], &poly::sub(f, &t1, &t2))
        };
        v.push(next);
    }
    v
}

/// The least `g ≥ 2` generating `(Z/ℓ)^*/±1`.
pub fn half_generator(ell: u64) -> u64 {
    let d = (ell - 1) / 2;
    (2..ell)
        .find(|&g| {
            let mut x = 1u64;
            for k in 1..=d {
                x = x * g % ell;
                if (x == 1 || x == ell - 1) && k < d {
                    return false;
                }
            }
            true
        })
        .unwrap_or(1)
}

/// Vélu's codomain from a kernel polynomial of odd degree `ℓ`.
pub fn velu_codomain(f: &Field, e: &Model, h: &Poly) -> Model {
    let d = poly::degree(h).unwrap_or(0);
    let ps = poly::power_sums(f, h, 3);
    let (p1, p2, p3) = (ps[0], ps[1], ps[2]);
    let dd = f.from_u64(d as u64);
    // v = Σ 2(3x² + a), w = Σ [4(x³ + ax + b) + 2x(3x² + a)].
    let v = f.add(
        &f.mul(&f.from_u64(6), &p2),
        &f.mul(&f.from_u64(2), &f.mul(&e.a, &dd)),
    );
    let w = f.add(
        &f.add(
            &f.mul(&f.from_u64(10), &p3),
            &f.mul(&f.from_u64(6), &f.mul(&e.a, &p1)),
        ),
        &f.mul(&f.from_u64(4), &f.mul(&e.b, &dd)),
    );
    Model {
        a: f.sub(&e.a, &f.mul(&f.from_u64(5), &v)),
        b: f.sub(&e.b, &f.mul(&f.from_u64(7), &w)),
    }
}

/// Certify that `h` is the kernel polynomial of an `F_p`-rational
/// `ℓ`-isogeny from `e`, and return Vélu's codomain.  Independent of `Φ_ℓ`.
pub fn verify_kernel(f: &Field, e: &Model, ell: u64, h: &Poly) -> Result<Model, EdgeError> {
    let d = ((ell - 1) / 2) as usize;
    if ell < 3 || ell.is_multiple_of(2) || poly::degree(h) != Some(d) || h[d] != f.one() {
        return Err(EdgeError::WrongDegree);
    }
    if poly::degree(&poly::gcd(f, h, &poly::derivative(f, h))) != Some(0) {
        return Err(EdgeError::NotSquarefree);
    }
    let g = half_generator(ell) as usize;
    let top = (ell as usize).max(g + 1);
    let fs = division_polys_mod(f, e, top, h);
    if !fs[ell as usize].is_empty() {
        return Err(EdgeError::NotTorsion);
    }
    if d > 1 {
        // x([g]P) = x − F f_{g−1} f_{g+1} / f_g²  (g odd)
        //         = x − f_{g−1} f_{g+1} / (F f_g²) (g even), F = 4(x³+ax+b).
        let big_f = poly::trim(
            f,
            vec![
                f.mul(&f.from_u64(4), &e.b),
                f.mul(&f.from_u64(4), &e.a),
                f.zero(),
                f.from_u64(4),
            ],
        );
        let mm = |x: &Poly, y: &Poly| poly::mulmod(f, x, y, h);
        let mut num = mm(&fs[g - 1], &fs[g + 1]);
        let mut den = mm(&fs[g], &fs[g]);
        if g % 2 == 1 {
            num = mm(&num, &big_f);
        } else {
            den = mm(&den, &big_f);
        }
        let den_inv = poly::invmod(f, &den, h).ok_or(EdgeError::NotSubgroup)?;
        let xg = poly::sub(f, &poly::x(f), &mm(&num, &den_inv));
        let image = poly::eval_mod(f, h, &poly::rem(f, &xg, h), h);
        if !image.is_empty() {
            return Err(EdgeError::NotSubgroup);
        }
    }
    Ok(velu_codomain(f, e, h))
}

/// `℘ = z^{-2} + Σ_{k≥1} c_k z^{2k}` for `y² = x³ + ax + b`, `k = 1..=n`.
fn wp_coefficients(f: &Field, m: &Model, n: usize) -> Vec<Fe> {
    let mut c = vec![f.zero(); n + 1];
    if n >= 1 {
        c[1] = f.div(&f.neg(&m.a), &f.from_u64(5)).unwrap();
    }
    if n >= 2 {
        c[2] = f.div(&f.neg(&m.b), &f.from_u64(7)).unwrap();
    }
    for k in 3..=n {
        let mut acc = f.zero();
        for j in 1..=k - 2 {
            acc = f.add(&acc, &f.mul(&c[j], &c[k - 1 - j]));
        }
        let den = f.from_u64(((k - 2) * (2 * k + 3)) as u64);
        c[k] = f.div(&f.mul(&f.from_u64(3), &acc), &den).unwrap();
    }
    c
}

/// The Elkies isogeny from `e` (with `j = j(e)`) to the curve with
/// `j`-invariant `j2`, a simple root of `Φ_ℓ(j, Y)`.  The result is
/// re-certified by [`verify_kernel`] before it is returned.
pub fn elkies_isogeny(f: &Field, e: &Model, phi: &ModPoly, j2: &Fe) -> Result<Isogeny, EdgeError> {
    let ell = phi.ell;
    let j = e.j(f).ok_or(EdgeError::SpecialJ)?;
    let k1728 = f.from_u64(1728);
    if f.is_zero(&j) || j == k1728 || f.is_zero(j2) || *j2 == k1728 {
        return Err(EdgeError::SpecialJ);
    }
    let l = f.from_u64(ell);
    let (px, py) = (phi.partial(f, &j, j2, 1, 0), phi.partial(f, &j, j2, 0, 1));
    if f.is_zero(&py) || f.is_zero(&px) {
        return Err(EdgeError::MultipleRoot);
    }
    let pxx = phi.partial(f, &j, j2, 2, 0);
    let pxy = phi.partial(f, &j, j2, 1, 1);
    let pyy = phi.partial(f, &j, j2, 0, 2);
    let div = |a: &Fe, b: &Fe| f.div(a, b).ok_or(EdgeError::SpecialJ);
    let e4 = div(&f.neg(&e.a), &f.from_u64(3))?;
    let e6 = div(&f.neg(&e.b), &f.from_u64(2))?;
    let jj = f.neg(&div(&f.mul(&j, &e6), &e4)?); // J = θj
    let jt = f.neg(&div(&f.mul(&px, &jj), &f.mul(&l, &py))?); // J̃
    let e4t = div(&f.sqr(&jt), &f.mul(j2, &f.sub(j2, &k1728)))?;
    let e6t = f.neg(&div(&f.mul(&jt, &e4t), j2)?);
    let l2 = f.sqr(&l);
    let l4 = f.sqr(&l2);
    let l6 = f.mul(&l4, &l2);
    let codomain = Model {
        a: f.mul(&f.from_i64(-3), &f.mul(&l4, &e4t)),
        b: f.mul(&f.from_i64(-2), &f.mul(&l6, &e6t)),
    };
    let two_thirds = div(&f.from_u64(2), &f.from_u64(3))?;
    let half = div(&f.one(), &f.from_u64(2))?;
    let term = |jv: &Fe, e4v: &Fe, e6v: &Fe| -> Result<Fe, EdgeError> {
        let r = div(&f.sqr(e6v), &f.sqr(e4v))?;
        Ok(f.add(
            &f.mul(&two_thirds, &f.mul(jv, &r)),
            &f.mul(&half, &f.mul(jv, e4v)),
        ))
    };
    let mut bracket = f.mul(&pxx, &f.sqr(&jj));
    bracket = f.add(
        &bracket,
        &f.mul(&f.from_u64(2), &f.mul(&l, &f.mul(&pxy, &f.mul(&jj, &jt)))),
    );
    bracket = f.add(&bracket, &f.mul(&l2, &f.mul(&pyy, &f.sqr(&jt))));
    bracket = f.add(&bracket, &f.mul(&px, &term(&j, &e4, &e6)?));
    bracket = f.add(&bracket, &f.mul(&l2, &f.mul(&py, &term(j2, &e4t, &e6t)?)));
    let g = div(
        &f.mul(&f.from_u64(6), &bracket),
        &f.mul(&l, &f.mul(&py, &jt)),
    )?; // E2 − ℓE2'
    let s1 = f.neg(&f.mul(&l, &g)); // Σ_{Q≠0} x(Q)
    let d = ((ell - 1) / 2) as usize;
    let kernel = kernel_from_sums(f, e, &codomain, s1, d)?;
    let certified = verify_kernel(f, e, ell, &kernel)?;
    if certified != codomain {
        return Err(EdgeError::CodomainMismatch);
    }
    Ok(Isogeny {
        ell,
        kernel,
        codomain,
    })
}

/// The kernel polynomial from the codomain and `S₁ = Σ_{Q≠0} x(Q)`.
fn kernel_from_sums(f: &Field, e: &Model, e2: &Model, s1: Fe, d: usize) -> Result<Poly, EdgeError> {
    if d == 0 {
        return Err(EdgeError::WrongDegree);
    }
    let c = wp_coefficients(f, e, d);
    let ct = wp_coefficients(f, e2, d);
    // F_0 = X; F_{k+1} = F_k''·(4X³ + 4aX + 4b) + F_k'·(6X² + 2a).
    let cubic = poly::trim(
        f,
        vec![
            f.mul(&f.from_u64(4), &e.b),
            f.mul(&f.from_u64(4), &e.a),
            f.zero(),
            f.from_u64(4),
        ],
    );
    let quad = poly::trim(
        f,
        vec![f.mul(&f.from_u64(2), &e.a), f.zero(), f.from_u64(6)],
    );
    // S_i = Σ_{Q≠0} x(Q)^i over all ℓ − 1 points.
    let mut s = vec![f.from_u64((2 * d) as u64), s1];
    let mut fk = poly::x(f);
    let mut fact = f.one(); // (2k)!
    for k in 1..d {
        let d1 = poly::derivative(f, &fk);
        let d2 = poly::derivative(f, &d1);
        fk = poly::add(f, &poly::mul(f, &d2, &cubic), &poly::mul(f, &d1, &quad));
        fact = f.mul(&fact, &f.from_u64(((2 * k - 1) * (2 * k)) as u64));
        // Σ_Q F_k(x_Q) = (2k)! (c'_k − c_k).
        let mut rhs = f.mul(&fact, &f.sub(&ct[k], &c[k]));
        for i in 0..=k {
            if let Some(coef) = fk.get(i) {
                rhs = f.sub(&rhs, &f.mul(coef, &s[i]));
            }
        }
        let lead = fk.get(k + 1).copied().ok_or(EdgeError::WrongDegree)?;
        s.push(f.div(&rhs, &lead).ok_or(EdgeError::WrongDegree)?);
    }
    let half = f.inv(&f.from_u64(2)).unwrap();
    let halves: Vec<Fe> = s[1..=d].iter().map(|v| f.mul(v, &half)).collect();
    Ok(poly::from_power_sums(f, &halves, d))
}

/// The `F_p`-isomorphism from Vélu's codomain to the recorded target, as
/// `u²`, or [`EdgeError::NotIsomorphic`].
pub fn link_to_target(f: &Field, codomain: &Model, target: &Model) -> Result<Fe, EdgeError> {
    curve::isomorphism(f, codomain, target).ok_or(EdgeError::NotIsomorphic)
}
