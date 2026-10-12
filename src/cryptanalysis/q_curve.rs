//! # Q-curves of degree 2 and 3 over `F_{p²}` (E14)
//!
//! A **Q-curve of degree `d`** is an elliptic curve `E/F_{p²}` with a
//! `d`-isogeny `φ: E → E^σ` to its Frobenius conjugate (`σ` the
//! `p`-power map on coefficients).  Composing with the `p`-power
//! Frobenius `π: E^σ → E` gives an endomorphism `ψ = π ∘ ι ∘ φ` of `E`
//! defined over `F_{p²}`, of degree `dp` (Smith, "Families of fast
//! elliptic curves from Q-curves", ASIACRYPT 2013; `ι: φ(E) ≅ E^σ` the
//! isomorphism that matches Vélu's codomain to the conjugate).  On
//! `E(F_{p²})` the `p²`-Frobenius is trivial, so `ψ²` is an endomorphism
//! of degree `d²`, `[±d]` for these curves, and `ψ` acts on a
//! prime-order subgroup `⟨G⟩` as an eigenvalue `λ` with `λ² ≡ ±d`.
//!
//! The question E14 asks of them is the type-C one of E4: does `ψ` keep
//! any point of a factor base in the base?  A GLS endomorphism (`d = 1`)
//! fixes a line `x ∈ u·s·F_p` and folds it; a degree-`d` isogeny's
//! `x`-map is a rational function of degree `d`, so no line is stable
//! and the prediction is the chance overlap.
//!
//! **Construction.**  The `j`-invariants of degree-`d` Q-curves are the
//! `j ∈ F_{p²} ∖ F_p` with `Φ_d(j, j^p) = 0`, `Φ_d` the classical modular
//! polynomial ([`super::modular_polynomial`]).  `Φ_d` is symmetric with
//! integer coefficients, so `Φ_d(j, j^p)` lies in `F_p`: one equation in
//! the two coordinates of `j`, with about `(d + 1)·p` solutions.  For
//! each, `E_j: y² = x³ + 3kx + 2k` (`k = j/(1728 − j)`) and its quadratic
//! twist are counted, the kernel of the `d`-isogeny to `E^σ` is found
//! among the `F_{p²}`-roots of the `d`-division polynomial and **checked**
//! by Vélu's codomain (`j(φ(E)) = j^p`, independent of `Φ_d`), and `ψ`'s
//! eigenvalue is read off `ψ(G)` by baby-step giant-step and verified on
//! random points like every other map in this thread.  `p ≤ 2^{12}`:
//! roots and orders are found by scanning `F_{p²}`.
//!
//! Companion to `research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`
//! §8.9.

use crate::cryptanalysis::ext_curve::{jacobi, ExtField, ExtPoint, Fp2, Fp2El};
use crate::cryptanalysis::gls_fp2::{Fp2Curve, Fp2Point};
use crate::cryptanalysis::glv_invariant_base::{factor_u64, mulm, Endomorphism};
use crate::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use crate::cryptanalysis::modular_polynomial::{phi_2, phi_3, ModularPolynomial};
use num_bigint::BigInt;
use num_traits::ToPrimitive;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;
use std::collections::HashMap;

fn is_prime(n: u64) -> bool {
    n >= 2 && factor_u64(n) == vec![(n, 1)]
}

/// `Φ_d` with its coefficients reduced modulo `p`, as `(i, j, c)`.
fn reduced_phi(d: u64, p: u64) -> Result<Vec<(usize, usize, u64)>, String> {
    let phi: ModularPolynomial = match d {
        2 => phi_2(),
        3 => phi_3(),
        _ => return Err(format!("degree {d}: only Φ₂ and Φ₃ are tabulated")),
    };
    let m = BigInt::from(p);
    let top = d as u32 + 1;
    let mut out = Vec::new();
    for i in 0..=top {
        for j in 0..=top {
            let c = ((phi.coeff(i, j) % &m) + &m) % &m;
            let c = c.to_u64().expect("reduced below p");
            if c != 0 {
                out.push((i as usize, j as usize, c));
            }
        }
    }
    Ok(out)
}

/// `Φ_d(x, y)` at `(x, y) ∈ F_{p²}²`.
fn eval_phi(f: &Fp2, phi: &[(usize, usize, u64)], x: Fp2El, y: Fp2El) -> Fp2El {
    let mut xp = [f.one(); 5];
    let mut yp = [f.one(); 5];
    for k in 1..5 {
        xp[k] = f.mul(xp[k - 1], x);
        yp[k] = f.mul(yp[k - 1], y);
    }
    phi.iter().fold(f.zero(), |acc, &(i, j, c)| {
        f.add(acc, f.scale(f.mul(xp[i], yp[j]), c))
    })
}

/// `j(E)` for `y² = x³ + ax + b`; `None` on a singular curve.
pub fn j_invariant(f: &Fp2, a: Fp2El, b: Fp2El) -> Option<Fp2El> {
    let a3 = f.scale(f.mul(f.sqr(a), a), 4);
    let den = f.add(a3, f.scale(f.sqr(b), 27));
    Some(f.scale(f.mul(a3, f.inv(den)?), 1728))
}

/// Every root in `F_{p²}` of a small polynomial (coefficients low degree
/// first), by scanning the field.
fn roots_by_scan(f: &Fp2, coeffs: &[Fp2El]) -> Vec<Fp2El> {
    let p = f.p;
    let mut out = Vec::new();
    for a0 in 0..p {
        for a1 in 0..p {
            let x = [a0, a1];
            let v = coeffs
                .iter()
                .rev()
                .fold(f.zero(), |acc, &c| f.add(f.mul(acc, x), c));
            if f.is_zero(v) {
                out.push(x);
            }
        }
    }
    out
}

/// `#E(F_{p²})`: the quadratic character of `F_{p²}` is the Legendre
/// symbol of the norm.
fn count_points(curve: &Fp2Curve) -> u64 {
    let f = &curve.f;
    let p = f.p;
    let mut n = 1u64;
    for a0 in 0..p {
        for a1 in 0..p {
            let rhs = curve.rhs([a0, a1]);
            n += (1 + jacobi(f.norm(rhs), p)) as u64;
        }
    }
    n
}

/// One point of the kernel's half `S` (Vélu): `x_Q`, `t_Q` and `u_Q`.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct VeluTerm {
    pub x: Fp2El,
    pub t: Fp2El,
    pub u: Fp2El,
}

/// `ψ = π ∘ ι ∘ φ` on a degree-`d` Q-curve.
#[derive(Clone, Debug, Serialize)]
pub struct QCurveEndomorphism {
    pub d: u64,
    pub p: u64,
    pub kernel: Vec<VeluTerm>,
    /// `ι: (X, Y) ↦ (μ²X, μ³Y)`.
    pub mu2: Fp2El,
    pub mu3: Fp2El,
    pub eigenvalue: u64,
    /// `λ² mod r`, and the sign `s` with `λ² ≡ s·d`.
    pub eigenvalue_squared: u64,
    pub square_is_plus_d: Option<bool>,
}

impl Endomorphism<Fp2Curve> for QCurveEndomorphism {
    fn name(&self) -> String {
        format!("q-curve-psi[d={},lambda={}]", self.d, self.eigenvalue)
    }
    fn degree(&self) -> u64 {
        self.d * self.p
    }
    fn eigenvalue(&self) -> u64 {
        self.eigenvalue
    }
    fn apply(&self, g: &Fp2Curve, pt: Fp2Point) -> Fp2Point {
        let f = &g.f;
        if pt.infinity {
            return pt;
        }
        // Vélu: X = x + Σ (t/(x − x_Q) + u/(x − x_Q)²),
        //       Y = y·(1 − Σ (t/(x − x_Q)² + 2u/(x − x_Q)³)).
        let (mut x_out, mut dfac) = (pt.x, f.one());
        for q in &self.kernel {
            let Some(di) = f.inv(f.sub(pt.x, q.x)) else {
                return ExtPoint::infinity(f.zero());
            };
            let di2 = f.sqr(di);
            let di3 = f.mul(di2, di);
            x_out = f.add(x_out, f.add(f.mul(q.t, di), f.mul(q.u, di2)));
            dfac = f.sub(dfac, f.add(f.mul(q.t, di2), f.scale(f.mul(q.u, di3), 2)));
        }
        let y_out = f.mul(pt.y, dfac);
        ExtPoint::affine(
            f.frob(f.mul(self.mu2, x_out)),
            f.frob(f.mul(self.mu3, y_out)),
        )
    }
}

/// Why a kernel that maps to `E^σ` gave no map, or how the search went.
#[derive(Clone, Debug, Default, Serialize)]
pub struct QCurveSearch {
    pub j_scanned: u64,
    pub q_curve_j_found: u64,
    pub curves_counted: u64,
    pub rejected_order: u64,
    /// Kernels whose codomain is `E^σ` only over `F_{p⁴}` (`μ² ∉ (F_{p²})²`):
    /// `ψ` would not be defined over `F_{p²}`.
    pub kernels_twist_only: u64,
    /// Roots of the division polynomial whose codomain is not `E^σ`.
    pub kernels_other_codomain: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct QCurveInstance {
    pub name: String,
    pub p: u64,
    pub d: u64,
    pub j: Fp2El,
    pub twisted: bool,
    pub curve: Fp2Curve,
    pub group_order: u64,
    pub r: u64,
    pub cofactor: u64,
    pub generator: Fp2Point,
    pub maps: Vec<QCurveEndomorphism>,
    pub search: QCurveSearch,
}

fn random_point(curve: &Fp2Curve, rng: &mut StdRng) -> Fp2Point {
    let p = curve.f.p;
    loop {
        let x = [rng.gen_range(0..p), rng.gen_range(0..p)];
        if let Some(&pt) = curve.lift_x(x).first() {
            return pt;
        }
    }
}

/// `log_G(T)` in `⟨G⟩` of prime order `r`, by baby-step giant-step.
pub fn dlog_bsgs(curve: &Fp2Curve, g: Fp2Point, t: Fp2Point, r: u64) -> Option<u64> {
    let mut ops = GroupOps::default();
    let m = (r as f64).sqrt().ceil() as u64 + 1;
    let mut baby: HashMap<u64, u64> = HashMap::with_capacity(m as usize);
    let mut cur = curve.identity();
    for k in 0..m {
        baby.entry(curve.key(&cur)).or_insert(k);
        cur = curve.add(&mut ops, cur, g);
    }
    let step = curve.neg(curve.mul(&mut ops, g, m));
    let mut gamma = t;
    for i in 0..=m {
        if let Some(&k) = baby.get(&curve.key(&gamma)) {
            let cand = (i * m + k) % r;
            if curve.mul(&mut ops, g, cand) == t {
                return Some(cand);
            }
        }
        gamma = curve.add(&mut ops, gamma, step);
    }
    None
}

/// The kernel half `S` of every `d`-isogeny from `curve` whose codomain is
/// `E^σ` over `F_{p²}`, with `ι`'s scalings.  `d = 2`: the roots of
/// `x³ + ax + b`; `d = 3`: the roots of `ψ₃ = 3x⁴ + 6ax² + 12bx − a²`, one
/// per `{Q, −Q}`.
fn conjugate_isogenies(
    curve: &Fp2Curve,
    d: u64,
    search: &mut QCurveSearch,
) -> Vec<(Vec<VeluTerm>, Fp2El, Fp2El)> {
    let f = curve.f;
    let (a, b) = (curve.a, curve.b);
    let poly: Vec<Fp2El> = match d {
        2 => vec![b, a, f.zero(), f.one()],
        _ => vec![
            f.neg(f.sqr(a)),
            f.scale(b, 12),
            f.scale(a, 6),
            f.zero(),
            f.from_fp(3),
        ],
    };
    let target_j = j_invariant(&f, f.frob(a), f.frob(b));
    let mut out = Vec::new();
    for x0 in roots_by_scan(&f, &poly) {
        let gx = f.add(f.scale(f.sqr(x0), 3), a);
        let term = if d == 2 {
            VeluTerm {
                x: x0,
                t: gx,
                u: f.zero(),
            }
        } else {
            VeluTerm {
                x: x0,
                t: f.scale(gx, 2),
                u: f.scale(curve.rhs(x0), 4),
            }
        };
        let t = term.t;
        let w = f.add(term.u, f.mul(term.x, term.t));
        let a2 = f.sub(a, f.scale(t, 5));
        let b2 = f.sub(b, f.scale(w, 7));
        if j_invariant(&f, a2, b2) != target_j || target_j.is_none() {
            search.kernels_other_codomain += 1;
            continue;
        }
        // ι: μ⁴ a' = a^σ, μ⁶ b' = b^σ, so μ² = b^σ a' / (a^σ b').
        let (asg, bsg) = (f.frob(a), f.frob(b));
        let (Some(ia), Some(ib)) = (f.inv(f.mul(asg, b2)), f.inv(a2)) else {
            continue;
        };
        let mu2 = f.mul(f.mul(bsg, a2), ia);
        let _ = ib;
        if f.mul(f.sqr(mu2), a2) != asg {
            // μ⁴a' ≠ a^σ: the codomain shares E^σ's j but not its model.
            search.kernels_other_codomain += 1;
            continue;
        }
        let Some(mu) = f.sqrt(mu2) else {
            search.kernels_twist_only += 1;
            continue;
        };
        out.push((vec![term], mu2, f.mul(mu2, mu)));
    }
    out
}

/// A degree-`d` Q-curve over `F_{p²}`, `p` of `p_bits` bits, with a
/// prime-order subgroup of cofactor at most `max_cofactor` and every
/// `ψ` the curve carries over `F_{p²}`, each verified.  Deterministic in
/// `seed`.
pub fn generate_q_curve_instance(
    d: u64,
    p_bits: u32,
    seed: u64,
    max_cofactor: u64,
) -> Result<QCurveInstance, String> {
    if !(6..=12).contains(&p_bits) {
        return Err(format!("p_bits = {p_bits} outside 6..=12"));
    }
    let mut rng = StdRng::seed_from_u64(seed ^ 0x5143_5552_5645 ^ (d << 40) ^ (p_bits as u64));
    let mut search = QCurveSearch::default();
    for _ in 0..64 {
        let p = loop {
            let c = rng.gen_range((1u64 << (p_bits - 1))..(1u64 << p_bits)) | 1;
            if c > 3 && is_prime(c) {
                break c;
            }
        };
        let f = Fp2::new(p)?;
        let phi = reduced_phi(d, p)?;
        // Scan j = j₀ + j₁s from a random j₀.
        let start = rng.gen_range(0..p);
        for da in 0..p {
            let j0 = (start + da) % p;
            for j1 in 1..p {
                search.j_scanned += 1;
                let j = [j0, j1];
                let v = eval_phi(&f, &phi, j, f.frob(j));
                debug_assert_eq!(v[1], 0, "Φ_d(j, j^p) lies in F_p");
                if !f.is_zero(v) {
                    continue;
                }
                search.q_curve_j_found += 1;
                let Some(k) = f.inv(f.sub(f.from_fp(1728 % p), j)).map(|i| f.mul(j, i)) else {
                    continue;
                };
                let base = Fp2Curve {
                    f,
                    a: f.scale(k, 3),
                    b: f.scale(k, 2),
                };
                let delta = f.nonresidue();
                let twist = Fp2Curve {
                    f,
                    a: f.mul(base.a, f.sqr(delta)),
                    b: f.mul(base.b, f.mul(f.sqr(delta), delta)),
                };
                for (twisted, curve) in [(false, base), (true, twist)] {
                    search.curves_counted += 1;
                    let n = count_points(&curve);
                    let factors = factor_u64(n);
                    let Some(&(r, 1)) = factors.last() else {
                        search.rejected_order += 1;
                        continue;
                    };
                    let h = n / r;
                    if h > max_cofactor || r <= p {
                        search.rejected_order += 1;
                        continue;
                    }
                    let isos = conjugate_isogenies(&curve, d, &mut search);
                    if isos.is_empty() {
                        continue;
                    }
                    let mut ops = GroupOps::default();
                    let generator = loop {
                        let g = curve.mul(&mut ops, random_point(&curve, &mut rng), h);
                        if !g.infinity {
                            break g;
                        }
                    };
                    if !curve.mul(&mut ops, generator, r).infinity {
                        return Err(format!("[r]G ≠ O at p = {p}: the count {n} is wrong"));
                    }
                    let mut maps = Vec::new();
                    for (kernel, mu2, mu3) in isos {
                        let mut psi = QCurveEndomorphism {
                            d,
                            p,
                            kernel,
                            mu2,
                            mu3,
                            eigenvalue: 0,
                            eigenvalue_squared: 0,
                            square_is_plus_d: None,
                        };
                        let image = psi.apply(&curve, generator);
                        let Some(lambda) = dlog_bsgs(&curve, generator, image, r) else {
                            return Err(format!(
                                "ψ(G) is not in ⟨G⟩ at p = {p}: the map is not an endomorphism"
                            ));
                        };
                        psi.eigenvalue = lambda;
                        psi.eigenvalue_squared = mulm(lambda, lambda, r);
                        psi.square_is_plus_d = if psi.eigenvalue_squared == d % r {
                            Some(true)
                        } else if psi.eigenvalue_squared == r - d % r {
                            Some(false)
                        } else {
                            None
                        };
                        maps.push(psi);
                    }
                    let tag = j0 ^ (j1 << 16) ^ (twisted as u64) << 40;
                    return Ok(QCurveInstance {
                        name: format!("qcurve-d{d}-p{p_bits}bit-p{p}-{tag:x}"),
                        p,
                        d,
                        j,
                        twisted,
                        curve,
                        group_order: n,
                        r,
                        cofactor: h,
                        generator,
                        maps,
                        search,
                    });
                }
            }
        }
    }
    Err(format!(
        "no degree-{d} Q-curve with cofactor ≤ {max_cofactor} at {p_bits} bits"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::glv_invariant_base::verify_endomorphism;

    #[test]
    fn modular_polynomials_vanish_on_known_isogenous_pairs() {
        // 1728 and 287496 = 66³ are 2-isogenous; 0 is 3-isogenous to
        // itself (Φ₃(0, 0) = 0) and 2-isogenous to 54000.
        let p = 1_000_003u64;
        let f = Fp2::new(p).unwrap();
        let phi2 = reduced_phi(2, p).unwrap();
        let phi3 = reduced_phi(3, p).unwrap();
        let v = |phi: &[(usize, usize, u64)], x: u64, y: u64| {
            eval_phi(&f, phi, f.from_fp(x % p), f.from_fp(y % p))
        };
        assert!(f.is_zero(v(&phi2, 1728, 287496)));
        assert!(f.is_zero(v(&phi2, 0, 54000)));
        assert!(f.is_zero(v(&phi3, 0, 0)));
        assert!(!f.is_zero(v(&phi2, 1728, 1729)));
    }

    #[test]
    fn q_curves_of_degree_2_and_3_carry_a_verified_endomorphism() {
        for d in [2u64, 3] {
            let inst = generate_q_curve_instance(d, 8, 1, 16).unwrap();
            let f = inst.curve.f;
            assert!(inst.j[1] != 0, "j ∉ F_p");
            assert!(!inst.maps.is_empty());
            for psi in &inst.maps {
                // The codomain check is independent of Φ_d.
                verify_endomorphism(&inst.curve, inst.generator, inst.r, psi, 20, 7)
                    .unwrap_or_else(|e| panic!("d = {d}: {e}"));
                assert!(psi.square_is_plus_d.is_some(), "λ² ≡ ±{d} (mod r)");
                let mut ops = GroupOps::default();
                let twice = psi.apply(&inst.curve, psi.apply(&inst.curve, inst.generator));
                let s = if psi.square_is_plus_d == Some(true) {
                    d
                } else {
                    inst.r - d
                };
                assert_eq!(twice, inst.curve.mul(&mut ops, inst.generator, s));
            }
            let _ = f;
        }
    }
}
