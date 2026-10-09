//! # The 2-torsion symmetry of the decomposition system on a Frobenius line (E13).
//!
//! Faugère, Gaudry, Huot and Renault (ePrint 2012/199) use a rational
//! 2-torsion point `T` to make the summation polynomials smaller: on a
//! coordinate `Y` with `Y(P + T) = −Y(P)` and `Y(−P) = Y(P)`, adding `T`
//! to an even number of summands leaves their sum alone, so `S_{m+1}`
//! is invariant under the group `D_m = (Z/2)^{m−1} ⋊ S_m` of even sign
//! changes and permutations of `(Y₁, …, Y_m)`, not only under `S_m`.
//! Written in the invariants of `D_m` instead of the elementary
//! symmetric functions, the system has `2^{m−1}` times fewer solutions.
//!
//! This module puts that symmetry on the Frobenius line of
//! `research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`
//! (§6.5, §8.2, §8.5), so the base fold and the system symmetry can be
//! measured together (E13 of §8.1):
//!
//! * **The coordinate.**  On a subfield curve `y² = x³ + ax + b` over
//!   `F_p` with a rational root `x₀` of the cubic and `f'(x₀) = 3x₀² + a
//!   = c²` a square in `F_p`, `x(P + T) − x₀ = c² / (x − x₀)` for
//!   `T = (x₀, 0)`, so `w = (x − x₀)/c` satisfies `w(P + T) = 1/w` and
//!   `Y = (w − 1)/(w + 1) = (x − x₀ − c)/(x − x₀ + c)` satisfies
//!   `Y(P + T) = −Y(P)`.  Its coefficients are in `F_p`, so it commutes
//!   with the Frobenius.
//! * **The base.**  `{P : Y(P) ∈ L·F_p^*}` for the line `L ∈ {s, s²}` of
//!   `F_{p³}` whose points survive the cofactor.  It is stable under
//!   negation, under `π` (`π(L t) = ω L t`) and under `τ_T: P ↦ P + T`
//!   (`Y ↦ −Y`).  `P` and `P + T` have the same `⟨G⟩`-component (`T` is
//!   in the cofactor), so `τ_T` folds with coefficient 1: the base folds
//!   `12` to a column under `⟨−1, π, τ_T⟩` against `6` under `⟨−1, π⟩`.
//! * **The system.**  `S₄(x₁, x₂, x₃, x_R)` from `gaudry_cubic` with
//!   `x_i = (A Y_i + B)/(C Y_i + D)` substituted (denominators cleared to
//!   degree 4 in each `Y_i`), then symmetrised two ways: in
//!   `e_k(Y)` (the `S₃` presentation, 64 solutions, solved by the
//!   existing Macaulay solver — [`YLineS4Oracle`]) and in the `D₃`
//!   invariants `p₁ = ΣY_i²`, `p₂ = ΣY_i²Y_j²`, `p₃ = Y₁Y₂Y₃` (16
//!   solutions — [`FghrOracle`]).  On the line `Y_i = L t_i` both scale
//!   to `F_p` unknowns: `e_k = L^k σ_k(t)` and `p₁ = L² q₁`,
//!   `p₂ = L⁴ q₂`, `p₃ = L³ q₃`.
//! * **The `D₃` solver.**  Each descended equation has weighted degree
//!   `≤ 4` for weights `(2, 2, 1)`: as a polynomial in `(q₁, q₂)` it is a
//!   conic whose coefficients are polynomials in `q₃`.  Three conics meet
//!   exactly where their resultant `R(q₃)` vanishes: `R` is computed by
//!   Macaulay's formula at enough values of `q₃` and interpolated, its
//!   roots found by Cantor–Zassenhaus, and the conics solved at each
//!   root.  Then `t_i²` are the roots of `Z³ − q₁Z² + q₂Z − q₃²` and the
//!   signs of `t_i` are fixed by `t₁t₂t₃ = q₃`, up to the even sign
//!   changes that the symmetry quotients out.
//!
//! Every decomposition either oracle returns is lifted to base points by
//! group arithmetic, so nothing is returned that does not sum to the
//! target.

use std::collections::HashMap;
use std::time::Instant;

use rand::rngs::StdRng;
use rand::SeedableRng;
use serde::Serialize;

use crate::cryptanalysis::ext_curve::{jacobi, ExtField, ExtPoint, Fp3El};
use crate::cryptanalysis::gaudry_cubic::{
    s4_terms, solve_s4_subspace_with, Curve3, Fp3 as SolverField, Instance3, Pt3, SolveMode,
    SolveStats, SymmetrisedS4, E3,
};
use crate::cryptanalysis::glv_invariant_base::{
    addm, fold_by_endomorphisms, mulm, poly_roots_fp_counted, subm, Closure, Endomorphism,
    FoldReport, Negation,
};
use crate::cryptanalysis::ic_boundary::{
    lift_abscissae, CountedGroup, FactorBase, GroupOps, OracleCounters,
};
use crate::cryptanalysis::ic_framework::stages::{
    DecompositionOracle, InstanceCtx, Params, SolverCost, SolverTotals, SystemShape,
};
use crate::cryptanalysis::residual_walk::{inv_mod, sqrt_mod};
use crate::cryptanalysis::subfield_fp3::{
    generate_subfield_instance, Fp3Curve, Fp3Point, SubfieldInstance,
};

// ── The instance and the coordinate ────────────────────────────────

/// A subfield instance with a rational 2-torsion point whose `f'` is a
/// square, the `Y` coordinate, and the line it is taken on.
#[derive(Clone, Debug, Serialize)]
pub struct FghrInstance {
    pub sub: SubfieldInstance,
    pub x0: u64,
    /// `f'(x₀) = 3x₀² + a = c²`.
    pub fprime: u64,
    pub c: u64,
    pub t: Fp3Point,
    /// `x = (A Y + B)/(C Y + D)`: `A = c − x₀`, `B = x₀ + c`, `C = −1`,
    /// `D = 1`.
    pub mobius: [u64; 4],
    /// The direction `L` of the `Y`-line.
    pub yline: Fp3El,
}

impl FghrInstance {
    /// `x` from `Y` on `E(F_{p³})`, on the field's counted arithmetic.
    pub fn x_of_y(&self, y: Fp3El) -> Fp3El {
        let f = &self.sub.curve.f;
        let [a, b, c, d] = self.mobius;
        let num = f.add(f.scale(y, a), f.from_fp(b));
        let den = f.add(f.scale(y, c), f.from_fp(d));
        f.mul(num, f.inv(den).expect("Y = 1 is not on the line"))
    }

    /// `Y` from `x`.
    pub fn y_of_x(&self, x: Fp3El) -> Fp3El {
        let f = &self.sub.curve.f;
        let p = f.p;
        let num = f.sub(x, f.from_fp(addm(self.x0, self.c, p)));
        let den = f.add(x, f.from_fp(subm(self.c, self.x0, p)));
        f.mul(num, f.inv(den).expect("x = x₀ − c has no Y"))
    }
}

/// A `Y`-line: the map from the line's `F_p^*` parameter `t` to the
/// abscissa of the base points above it.  The oracles below see an
/// instance only through it, so one oracle serves the subfield line of
/// E13 and the full-group line of E17 (`fghr_full`).
pub trait YLine {
    /// `x` above the line's parameter `t`, on the field's counted
    /// arithmetic; `None` where `t` is the Möbius map's pole (`Y = 1`,
    /// the point at infinity), which only a line through `1` meets.
    fn x_at(&self, t: u64) -> Option<Fp3El>;
}

impl YLine for FghrInstance {
    fn x_at(&self, t: u64) -> Option<Fp3El> {
        Some(self.x_of_y(self.sub.curve.f.scale(self.yline, t)))
    }
}

/// A subfield instance (`generate_subfield_instance`) whose base curve
/// has a root `x₀ ∈ F_p` of `x³ + ax + b` with `3x₀² + a` a nonzero
/// square; tries successive seeds.
pub fn generate_fghr_instance(
    p_bits: u32,
    seed: u64,
    max_cofactor_ratio: u64,
) -> Result<FghrInstance, String> {
    for k in 0..400u64 {
        let Ok(sub) = generate_subfield_instance(
            p_bits,
            seed.wrapping_mul(1009).wrapping_add(k),
            false,
            max_cofactor_ratio,
        ) else {
            continue;
        };
        let p = sub.p;
        let (a, b) = (sub.base_a, sub.base_b);
        let mut muls = 0u64;
        let roots = poly_roots_fp_counted(&[b, a, 0, 1], p, seed ^ k, &mut muls);
        for x0 in roots {
            let fprime = addm(mulm(3, mulm(x0, x0, p), p), a, p);
            if fprime == 0 || jacobi(fprime, p) != 1 {
                continue;
            }
            let Some(c) = sqrt_mod(fprime, p) else {
                continue;
            };
            let f = sub.curve.f;
            let t = ExtPoint::affine(f.from_fp(x0), f.zero());
            let mobius = [subm(c, x0, p), addm(x0, c, p), p - 1, 1];
            let mut inst = FghrInstance {
                sub,
                x0,
                fprime,
                c,
                t,
                mobius,
                yline: f.s(),
            };
            // The line whose points survive the cofactor.
            let survivors = |inst: &FghrInstance, line: Fp3El| -> usize {
                let mut n = 0usize;
                let mut o = GroupOps::default();
                for tt in 1..inst.sub.p.min(64) {
                    let x = inst.x_of_y(f.scale(line, tt));
                    for pt in inst.sub.curve.lift_x(x) {
                        if !inst.sub.curve.mul(&mut o, pt, inst.sub.cofactor).infinity {
                            n += 1;
                        }
                    }
                }
                n
            };
            let s1 = f.s();
            let s2 = f.mul(s1, s1);
            if survivors(&inst, s2) > survivors(&inst, s1) {
                inst.yline = s2;
            }
            // The translation acts as Y ↦ −Y: checked on points, not assumed.
            let mut rng = StdRng::seed_from_u64(seed ^ 0xF6);
            let mut ops = GroupOps::default();
            for _ in 0..8 {
                let pt = crate::cryptanalysis::ext_curve::random_point(&inst.sub.curve, &mut rng);
                let q = inst.sub.curve.add(&mut ops, pt, inst.t);
                if q.infinity || f.is_zero(f.add(q.x, f.from_fp(subm(inst.c, x0, p)))) {
                    continue;
                }
                let (y1, y2) = (inst.y_of_x(pt.x), inst.y_of_x(q.x));
                if f.add(y1, y2) != f.zero() {
                    return Err("Y(P + T) ≠ −Y(P): the coordinate is wrong".into());
                }
            }
            return Ok(inst);
        }
    }
    Err(format!(
        "no subfield instance with a usable rational 2-torsion point at p_bits = {p_bits}"
    ))
}

/// `τ_T: P ↦ P + T`, eigenvalue 1 on `⟨G⟩` (`T` is in the cofactor).
pub struct Translation {
    pub t: Fp3Point,
}

impl Endomorphism<Fp3Curve> for Translation {
    fn name(&self) -> String {
        "tau_T".into()
    }
    fn degree(&self) -> u64 {
        1
    }
    fn eigenvalue(&self) -> u64 {
        1
    }
    fn apply(&self, g: &Fp3Curve, p: Fp3Point) -> Fp3Point {
        let mut ops = GroupOps::default();
        g.add(&mut ops, p, self.t)
    }
}

/// Which group the `Y`-line base is folded by.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum FghrFold {
    Negation,
    Frobenius,
    FrobeniusTranslation,
}

impl FghrFold {
    pub fn name(self) -> &'static str {
        match self {
            Self::Negation => "negation",
            Self::Frobenius => "frobenius",
            Self::FrobeniusTranslation => "frobenius+translation",
        }
    }
}

/// The base `{P : Y(P) ∈ L·F_p^*}`, folded by `fold`; the set is checked
/// to be invariant under the group, not assumed.  The same points in the
/// same order for every fold.
pub fn fghr_line_base(
    inst: &FghrInstance,
    fold: FghrFold,
) -> Result<(FactorBase<Fp3Point>, FoldReport), String> {
    let start = Instant::now();
    let sub = &inst.sub;
    let curve = &sub.curve;
    let f = &curve.f;
    let mut seed = Vec::with_capacity(2 * sub.p as usize);
    for t in 1..sub.p {
        seed.extend(curve.lift_x(inst.x_of_y(f.scale(inst.yline, t))));
    }
    let neg = Negation { r: sub.r };
    let tau = Translation { t: inst.t };
    let mut gens: Vec<&dyn Endomorphism<Fp3Curve>> = vec![&neg];
    match fold {
        FghrFold::Negation => {}
        FghrFold::Frobenius => gens.push(&sub.frobenius),
        FghrFold::FrobeniusTranslation => {
            gens.push(&sub.frobenius);
            gens.push(&tau);
        }
    }
    let (mut fb, report) = fold_by_endomorphisms(
        curve,
        sub.r,
        sub.cofactor,
        seed,
        &gens,
        Closure::Strict,
        32,
        |p| curve.key(p),
        |p| f.pack(p.x),
        format!(
            "the Y-line Y ∈ L·F_p ({} values), Y = (x − x0 − c)/(x − x0 + c), folded by {}",
            sub.p - 1,
            fold.name()
        ),
    )?;
    fb.cost.wall_ns = start.elapsed().as_nanos() as u64;
    fb.cost.count("sqrt_solves", sub.p - 1);
    Ok((fb, report))
}

// ── The summation polynomial in Y, symmetrised two ways ────────────

type Poly4 = HashMap<[u8; 4], E3>;
type Poly3 = HashMap<[u8; 3], E3>;

fn add_term4(f: &SolverField, p: &mut Poly4, e: [u8; 4], c: E3) {
    let entry = p.entry(e).or_insert(SolverField::ZERO);
    *entry = f.add(entry, &c);
    if *entry == SolverField::ZERO {
        p.remove(&e);
    }
}

fn mul3(f: &SolverField, a: &Poly3, b: &Poly3) -> Poly3 {
    let mut out: Poly3 = HashMap::new();
    for (ea, ca) in a {
        for (eb, cb) in b {
            let e = [ea[0] + eb[0], ea[1] + eb[1], ea[2] + eb[2]];
            let entry = out.entry(e).or_insert(SolverField::ZERO);
            *entry = f.add(entry, &f.mul(ca, cb));
        }
    }
    out.retain(|_, c| *c != SolverField::ZERO);
    out
}

fn monomial3(e: [u8; 3]) -> Poly3 {
    HashMap::from([(e, SolverField::ONE)])
}

fn sum3(f: &SolverField, parts: &[Poly3]) -> Poly3 {
    let mut out: Poly3 = HashMap::new();
    for part in parts {
        for (e, c) in part {
            let entry = out.entry(*e).or_insert(SolverField::ZERO);
            *entry = f.add(entry, c);
        }
    }
    out.retain(|_, c| *c != SolverField::ZERO);
    out
}

/// `S₄` with `x_i = (A Y_i + B)/(C Y_i + D)` for `i = 1, 2, 3`,
/// multiplied through by `Π (C Y_i + D)⁴`: degree ≤ 4 in each `Y_i`;
/// the fourth variable stays `x_R`.
pub fn s4_in_y(f: &SolverField, s4x: &Poly4, mobius: [u64; 4]) -> Poly4 {
    let p = f.p;
    let [a, b, c, d] = mobius;
    // u_k(Y) = (A Y + B)^k (C Y + D)^{4 − k}, coefficients in F_p.
    let mul1 = |u: &[u64], lin: [u64; 2]| -> Vec<u64> {
        let mut out = vec![0u64; u.len() + 1];
        for (i, &x) in u.iter().enumerate() {
            out[i] = addm(out[i], mulm(x, lin[0], p), p);
            out[i + 1] = addm(out[i + 1], mulm(x, lin[1], p), p);
        }
        out
    };
    let u: Vec<Vec<u64>> = (0..=4usize)
        .map(|k| {
            let mut v = vec![1u64];
            for _ in 0..k {
                v = mul1(&v, [b, a]);
            }
            for _ in k..4 {
                v = mul1(&v, [d, c]);
            }
            v
        })
        .collect();
    let mut out: Poly4 = HashMap::new();
    for (e, coef) in s4x {
        let (u1, u2, u3) = (&u[e[0] as usize], &u[e[1] as usize], &u[e[2] as usize]);
        for (i, &c1) in u1.iter().enumerate() {
            if c1 == 0 {
                continue;
            }
            for (j, &c2) in u2.iter().enumerate() {
                if c2 == 0 {
                    continue;
                }
                let c12 = mulm(c1, c2, p);
                for (k, &c3) in u3.iter().enumerate() {
                    if c3 == 0 {
                        continue;
                    }
                    let scalar = mulm(c12, c3, p);
                    add_term4(
                        f,
                        &mut out,
                        [i as u8, j as u8, k as u8, e[3]],
                        f.scale(coef, scalar),
                    );
                }
            }
        }
    }
    out
}

/// [`s4_in_y`] with a Möbius map whose coefficients are in `F_{p³}`
/// (the full-group line of `fghr_full`, where `x₀ ∉ F_p`); the same
/// substitution, on the solver field's arithmetic.
pub fn s4_in_y_ext(f: &SolverField, s4x: &Poly4, mobius: [E3; 4]) -> Poly4 {
    let [a, b, c, d] = mobius;
    let mul1 = |u: &[E3], lin: [E3; 2]| -> Vec<E3> {
        let mut out = vec![SolverField::ZERO; u.len() + 1];
        for (i, x) in u.iter().enumerate() {
            out[i] = f.add(&out[i], &f.mul(x, &lin[0]));
            out[i + 1] = f.add(&out[i + 1], &f.mul(x, &lin[1]));
        }
        out
    };
    let u: Vec<Vec<E3>> = (0..=4usize)
        .map(|k| {
            let mut v = vec![SolverField::ONE];
            for _ in 0..k {
                v = mul1(&v, [b, a]);
            }
            for _ in k..4 {
                v = mul1(&v, [d, c]);
            }
            v
        })
        .collect();
    let mut out: Poly4 = HashMap::new();
    for (e, coef) in s4x {
        let (u1, u2, u3) = (&u[e[0] as usize], &u[e[1] as usize], &u[e[2] as usize]);
        for (i, c1) in u1.iter().enumerate() {
            if *c1 == SolverField::ZERO {
                continue;
            }
            for (j, c2) in u2.iter().enumerate() {
                if *c2 == SolverField::ZERO {
                    continue;
                }
                let c12 = f.mul(c1, c2);
                for (k, c3) in u3.iter().enumerate() {
                    if *c3 == SolverField::ZERO {
                        continue;
                    }
                    let scalar = f.mul(&c12, c3);
                    add_term4(
                        f,
                        &mut out,
                        [i as u8, j as u8, k as u8, e[3]],
                        f.mul(coef, &scalar),
                    );
                }
            }
        }
    }
    out
}

/// Rewrite a polynomial symmetric in its first three variables in a set
/// of invariant generators by repeatedly cancelling the lex-leading
/// term: `lead_to_key` maps a leading exponent to the generator exponent
/// whose product has that leading term (or refuses it), `expand` gives
/// that product as a polynomial in `Y₁, Y₂, Y₃`.
fn rewrite_in_invariants(
    f: &SolverField,
    mut g: Poly4,
    lead_to_key: impl Fn([u8; 3]) -> Result<[u8; 3], String>,
    expand: impl Fn([u8; 3]) -> Poly3,
) -> Result<Poly4, String> {
    let mut cache: HashMap<[u8; 3], Poly3> = HashMap::new();
    let mut out: Poly4 = HashMap::new();
    while let Some((&lead, &coef)) = g.iter().max_by_key(|(e, _)| **e) {
        let key = lead_to_key([lead[0], lead[1], lead[2]])?;
        let expansion = cache.entry(key).or_insert_with(|| expand(key)).clone();
        add_term4(f, &mut out, [key[0], key[1], key[2], lead[3]], coef);
        for (e, c) in &expansion {
            add_term4(
                f,
                &mut g,
                [e[0], e[1], e[2], lead[3]],
                f.neg(&f.mul(c, &coef)),
            );
        }
        if g.contains_key(&lead) {
            return Err(format!("the leading term {lead:?} did not cancel"));
        }
    }
    Ok(out)
}

/// In `(e₁, e₂, e₃, x_R)`: the `S₃` presentation.
pub fn symmetrise_s3(f: &SolverField, g: Poly4) -> Result<Poly4, String> {
    let e1 = sum3(
        f,
        &[
            monomial3([1, 0, 0]),
            monomial3([0, 1, 0]),
            monomial3([0, 0, 1]),
        ],
    );
    let e2 = sum3(
        f,
        &[
            monomial3([1, 1, 0]),
            monomial3([1, 0, 1]),
            monomial3([0, 1, 1]),
        ],
    );
    let e3 = monomial3([1, 1, 1]);
    rewrite_in_invariants(
        f,
        g,
        |[a, b, c]| {
            if a >= b && b >= c {
                Ok([a - b, b - c, c])
            } else {
                Err(format!("not symmetric: leading exponent ({a}, {b}, {c})"))
            }
        },
        |[i, j, k]| {
            let mut acc = monomial3([0, 0, 0]);
            for _ in 0..i {
                acc = mul3(f, &acc, &e1);
            }
            for _ in 0..j {
                acc = mul3(f, &acc, &e2);
            }
            for _ in 0..k {
                acc = mul3(f, &acc, &e3);
            }
            acc
        },
    )
}

/// In `(p₁, p₂, p₃, x_R)` with `p₁ = ΣY_i²`, `p₂ = ΣY_i²Y_j²`,
/// `p₃ = Y₁Y₂Y₃`: the `D₃` presentation.  Refuses a polynomial that is
/// not invariant under even sign changes.
pub fn symmetrise_d3(f: &SolverField, g: Poly4) -> Result<Poly4, String> {
    let p1 = sum3(
        f,
        &[
            monomial3([2, 0, 0]),
            monomial3([0, 2, 0]),
            monomial3([0, 0, 2]),
        ],
    );
    let p2 = sum3(
        f,
        &[
            monomial3([2, 2, 0]),
            monomial3([2, 0, 2]),
            monomial3([0, 2, 2]),
        ],
    );
    let p3 = monomial3([1, 1, 1]);
    rewrite_in_invariants(
        f,
        g,
        |[a, b, c]| {
            if !(a >= b && b >= c) {
                Err(format!("not symmetric: leading exponent ({a}, {b}, {c})"))
            } else if (a - b) % 2 != 0 || (b - c) % 2 != 0 {
                Err(format!(
                    "not invariant under even sign changes: leading exponent ({a}, {b}, {c})"
                ))
            } else {
                Ok([(a - b) / 2, (b - c) / 2, c])
            }
        },
        |[i, j, k]| {
            let mut acc = monomial3([0, 0, 0]);
            for _ in 0..i {
                acc = mul3(f, &acc, &p1);
            }
            for _ in 0..j {
                acc = mul3(f, &acc, &p2);
            }
            for _ in 0..k {
                acc = mul3(f, &acc, &p3);
            }
            acc
        },
    )
}

/// The polynomials an instance's oracles need, built once.
pub struct FghrPolynomials {
    pub field: SolverField,
    pub curve3: Curve3,
    /// `S₃` presentation on the line (`e_k ↦ L^k σ_k`), for the Macaulay solver.
    pub s3_on_line: SymmetrisedS4,
    /// `D₃` presentation on the line: `(i, j, k, d) ↦ coefficient` of
    /// `q₁^i q₂^j q₃^k x_R^d`.
    pub d3_on_line: Poly4,
    pub y_terms: usize,
    pub s3_terms: usize,
    pub d3_terms: usize,
    /// `F_p` multiplications the construction cost (the solver field's
    /// counter): a once-per-curve set-up.
    pub setup_muls: u64,
}

pub fn fghr_polynomials(inst: &FghrInstance) -> Result<FghrPolynomials, String> {
    let sub = &inst.sub;
    fghr_polynomials_on(
        &sub.curve,
        sub.r,
        sub.generator,
        inst.yline,
        Mobius::Fp(inst.mobius),
    )
}

/// The coefficients of `x = (A Y + B)/(C Y + D)`: in `F_p` on a subfield
/// curve (E13), in `F_{p³}` on the full group (E17).
#[derive(Clone, Copy, Debug)]
pub enum Mobius {
    Fp([u64; 4]),
    Fp3([Fp3El; 4]),
}

/// [`FghrPolynomials`] for any curve over `F_{p³}` with a `Y`-coordinate
/// `x = Mobius(Y)` on which a rational `2`-torsion translation acts as
/// `Y ↦ −Y`, taken on the line `Y ∈ L·F_p`.
pub fn fghr_polynomials_on(
    curve: &Fp3Curve,
    r: u64,
    generator: Fp3Point,
    yline: Fp3El,
    mobius: Mobius,
) -> Result<FghrPolynomials, String> {
    let field = SolverField::with_cube_nonresidue(curve.f.p, curve.f.nu)?;
    let g = Pt3::affine(E3(generator.x), E3(generator.y));
    let curve3 = Curve3::new(field.clone(), E3(curve.a), E3(curve.b), r, g);
    let f = &curve3.field;
    f.reset_muls();
    let s4x = s4_terms(&curve3);
    let y = match mobius {
        Mobius::Fp(m) => s4_in_y(f, &s4x, m),
        Mobius::Fp3(m) => s4_in_y_ext(f, &s4x, m.map(E3)),
    };
    let y_terms = y.len();
    let s3 = symmetrise_s3(f, y.clone())?;
    let d3 = symmetrise_d3(f, y)?;
    let l = E3(yline);
    let mut pow = vec![SolverField::ONE];
    for _ in 0..24 {
        let last = *pow.last().unwrap();
        pow.push(f.mul(&last, &l));
    }
    let s3_terms = s3.len();
    let s3_on_line = SymmetrisedS4::from_terms(s3).on_line(f, &l);
    let d3_on_line: Poly4 = d3
        .into_iter()
        .map(|(e, c)| {
            let k = 2 * e[0] as usize + 4 * e[1] as usize + 3 * e[2] as usize;
            (e, f.mul(&c, &pow[k]))
        })
        .collect();
    for e in d3_on_line.keys() {
        if 2 * e[0] + 2 * e[1] + e[2] > 4 {
            return Err(format!(
                "the D₃ presentation has a term of weighted degree > 4: {e:?}"
            ));
        }
    }
    let setup_muls = f.muls();
    Ok(FghrPolynomials {
        field: field.clone(),
        d3_terms: d3_on_line.len(),
        curve3,
        s3_on_line,
        d3_on_line,
        y_terms,
        s3_terms,
        setup_muls,
    })
}

// ── The D₃ solver: three conics with coefficients in F_p[q₃] ───────

/// `F_p` arithmetic with a multiplication counter.
struct Fp {
    p: u64,
    muls: u64,
}

impl Fp {
    fn mul(&mut self, a: u64, b: u64) -> u64 {
        self.muls += 1;
        mulm(a, b, self.p)
    }
    fn add(&self, a: u64, b: u64) -> u64 {
        addm(a, b, self.p)
    }
    fn sub(&self, a: u64, b: u64) -> u64 {
        subm(a, b, self.p)
    }
    fn inv(&mut self, a: u64) -> u64 {
        self.muls += 1;
        inv_mod(a, self.p)
    }
    fn eval(&mut self, poly: &[u64], z: u64) -> u64 {
        let mut acc = 0u64;
        for &c in poly.iter().rev() {
            let t = self.mul(acc, z);
            acc = self.add(t, c);
        }
        acc
    }
    /// Determinant by Gaussian elimination.
    fn det(&mut self, mut m: Vec<Vec<u64>>) -> u64 {
        let n = m.len();
        let mut det = 1u64;
        for col in 0..n {
            let Some(piv) = (col..n).find(|&r| m[r][col] != 0) else {
                return 0;
            };
            if piv != col {
                m.swap(piv, col);
                det = self.sub(0, det);
            }
            det = self.mul(det, m[col][col]);
            let inv = self.inv(m[col][col]);
            for r in col + 1..n {
                if m[r][col] == 0 {
                    continue;
                }
                let factor = self.mul(m[r][col], inv);
                for k in col..n {
                    let t = self.mul(factor, m[col][k]);
                    m[r][k] = self.sub(m[r][k], t);
                }
            }
        }
        det
    }
}

/// The six coefficients of a conic in `(X₀, X₁, X₂) = (1, q₁, q₂)`,
/// order `[X₀², X₀X₁, X₀X₂, X₁², X₁X₂, X₂²]`.
type Conic = [u64; 6];

/// One equation `Σ c_{ijk} q₁^i q₂^j q₃^k` as a conic in `(q₁, q₂)` with
/// coefficients in `F_p[q₃]`: `[slot] → polynomial in q₃` (low to high).
struct ConicPoly {
    slots: [Vec<u64>; 6],
}

fn slot_of(i: u8, j: u8) -> usize {
    match (i, j) {
        (0, 0) => 0,
        (1, 0) => 1,
        (0, 1) => 2,
        (2, 0) => 3,
        (1, 1) => 4,
        (0, 2) => 5,
        _ => unreachable!("weighted degree ≤ 4 keeps i + j ≤ 2"),
    }
}

impl ConicPoly {
    fn from_terms(terms: &HashMap<[u8; 3], u64>) -> Self {
        let mut slots: [Vec<u64>; 6] = Default::default();
        for (e, &c) in terms {
            let s = slot_of(e[0], e[1]);
            let k = e[2] as usize;
            if slots[s].len() <= k {
                slots[s].resize(k + 1, 0);
            }
            slots[s][k] = c;
        }
        Self { slots }
    }
    fn at(&self, fp: &mut Fp, z: u64) -> Conic {
        let mut out = [0u64; 6];
        for (s, poly) in self.slots.iter().enumerate() {
            out[s] = fp.eval(poly, z);
        }
        out
    }
    fn degree_in_q3(&self) -> usize {
        self.slots
            .iter()
            .map(|s| s.len().saturating_sub(1))
            .max()
            .unwrap_or(0)
    }
}

/// The resultant of three ternary quadratic forms by Macaulay's formula:
/// the `15 × 15` matrix of `(m / X_v²)·f_v` over the degree-4 monomials
/// `m` (`v` the first variable with `X_v² | m`), divided by the `3 × 3`
/// minor on the monomials divisible by two squares.  `None` when that
/// minor vanishes at this point.
fn conic_resultant(fp: &mut Fp, f: &[Conic; 3]) -> Option<u64> {
    // Quadratic monomials in slot order.
    const Q: [[u8; 3]; 6] = [
        [2, 0, 0],
        [1, 1, 0],
        [1, 0, 1],
        [0, 2, 0],
        [0, 1, 1],
        [0, 0, 2],
    ];
    let mut quartic: Vec<[u8; 3]> = Vec::with_capacity(15);
    for a in (0..=4u8).rev() {
        for b in (0..=(4 - a)).rev() {
            quartic.push([a, b, 4 - a - b]);
        }
    }
    let index = |m: [u8; 3]| quartic.iter().position(|&q| q == m).expect("degree 4");
    let mut mat = vec![vec![0u64; 15]; 15];
    let mut reduced = [true; 15];
    for (row, &m) in quartic.iter().enumerate() {
        let squares = (0..3).filter(|&v| m[v] >= 2).count();
        reduced[row] = squares == 1;
        let v = (0..3)
            .find(|&v| m[v] >= 2)
            .expect("a degree-4 monomial has a square");
        let mut shift = m;
        shift[v] -= 2;
        for (s, q) in Q.iter().enumerate() {
            let col = index([shift[0] + q[0], shift[1] + q[1], shift[2] + q[2]]);
            mat[row][col] = f[v][s];
        }
    }
    let non: Vec<usize> = (0..15).filter(|&i| !reduced[i]).collect();
    let minor: Vec<Vec<u64>> = non
        .iter()
        .map(|&r| non.iter().map(|&c| mat[r][c]).collect())
        .collect();
    let dm = fp.det(minor);
    if dm == 0 {
        return None;
    }
    let d = fp.det(mat);
    let inv = fp.inv(dm);
    Some(fp.mul(d, inv))
}

/// Newton interpolation through `(x_i, y_i)`: coefficients low to high.
fn interpolate(fp: &mut Fp, xs: &[u64], ys: &[u64]) -> Vec<u64> {
    let n = xs.len();
    let mut coef = ys.to_vec();
    for j in 1..n {
        for i in (j..n).rev() {
            let num = fp.sub(coef[i], coef[i - 1]);
            let den = fp.sub(xs[i], xs[i - j]);
            let inv = fp.inv(den);
            coef[i] = fp.mul(num, inv);
        }
    }
    // Newton form to monomial form.
    let mut out = vec![0u64; n];
    for i in (0..n).rev() {
        // out = out · (x − x_i) + coef[i]
        let mut next = vec![0u64; n];
        for k in 0..n {
            if out[k] == 0 {
                continue;
            }
            if k + 1 < n {
                next[k + 1] = fp.add(next[k + 1], out[k]);
            }
            let t = fp.mul(out[k], xs[i]);
            next[k] = fp.sub(next[k], t);
        }
        next[0] = fp.add(next[0], coef[i]);
        out = next;
    }
    while out.last() == Some(&0) {
        out.pop();
    }
    out
}

/// A univariate polynomial in `q₂` whose coefficients are polynomials in
/// `q₁`: `A₂ q₂² + A₁(q₁) q₂ + A₀(q₁)` for a conic at fixed `q₃`.
fn conic_in_q2(c: &Conic) -> [Vec<u64>; 3] {
    [vec![c[0], c[1], c[3]], vec![c[2], c[4]], vec![c[5]]]
}

fn padd(fp: &Fp, a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut out = vec![0u64; a.len().max(b.len())];
    for (i, &x) in a.iter().enumerate() {
        out[i] = fp.add(out[i], x);
    }
    for (i, &x) in b.iter().enumerate() {
        out[i] = fp.add(out[i], x);
    }
    out
}

fn psub(fp: &Fp, a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut out = vec![0u64; a.len().max(b.len())];
    for (i, &x) in a.iter().enumerate() {
        out[i] = fp.add(out[i], x);
    }
    for (i, &x) in b.iter().enumerate() {
        out[i] = fp.sub(out[i], x);
    }
    out
}

fn pmul(fp: &mut Fp, a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            let t = fp.mul(x, y);
            out[i + j] = fp.add(out[i + j], t);
        }
    }
    out
}

fn pscale(fp: &mut Fp, a: &[u64], k: u64) -> Vec<u64> {
    a.iter().map(|&x| fp.mul(x, k)).collect()
}

fn is_zero_poly(a: &[u64]) -> bool {
    a.iter().all(|&x| x == 0)
}

fn conic_value(fp: &mut Fp, c: &Conic, q1: u64, q2: u64) -> u64 {
    let (q11, q12, q22) = (fp.mul(q1, q1), fp.mul(q1, q2), fp.mul(q2, q2));
    let terms = [
        c[0],
        fp.mul(c[1], q1),
        fp.mul(c[2], q2),
        fp.mul(c[3], q11),
        fp.mul(c[4], q12),
        fp.mul(c[5], q22),
    ];
    terms.iter().fold(0, |acc, &t| fp.add(acc, t))
}

/// Common zeros `(q₁, q₂) ∈ F_p²` of three conics.
fn solve_three_conics(fp: &mut Fp, cs: &[Conic; 3], seed: u64) -> Vec<(u64, u64)> {
    let p = fp.p;
    let mut out = Vec::new();
    for (ia, ib) in [(0usize, 1usize), (0, 2), (1, 2)] {
        let (ga, gb) = (conic_in_q2(&cs[ia]), conic_in_q2(&cs[ib]));
        let (a2, b2) = (ga[2][0], gb[2][0]);
        // Eliminate q₂²: lin = b₂·ga − a₂·gb is linear in q₂.
        let (lin, quad) = if a2 == 0 && b2 == 0 {
            (ga.clone(), gb.clone())
        } else {
            let (a0, b0) = (pscale(fp, &ga[0], b2), pscale(fp, &gb[0], a2));
            let (a1, b1) = (pscale(fp, &ga[1], b2), pscale(fp, &gb[1], a2));
            let l: [Vec<u64>; 3] = [psub(fp, &a0, &b0), psub(fp, &a1, &b1), vec![0]];
            (l, if a2 != 0 { ga.clone() } else { gb.clone() })
        };
        let (h0, h1) = (&lin[0], &lin[1]);
        // P(q₁) = A₂ H₀² − A₁ H₀ H₁ + A₀ H₁²: zero at every common solution.
        let h0h0 = pmul(fp, h0, h0);
        let h0h1 = pmul(fp, h0, h1);
        let h1h1 = pmul(fp, h1, h1);
        let t2 = pscale(fp, &h0h0, quad[2][0]);
        let t1 = pmul(fp, &quad[1], &h0h1);
        let t0 = pmul(fp, &quad[0], &h1h1);
        let poly = padd(fp, &psub(fp, &t2, &t1), &t0);
        if is_zero_poly(&poly) {
            continue;
        }
        let q1s = poly_roots_fp_counted(&poly, p, seed, &mut fp.muls);
        for q1 in q1s {
            // q₂ from the first conic at this q₁, checked on all three.
            let coeffs_a = conic_in_q2(&cs[ia]);
            let at: Vec<u64> = (0..3).map(|k| fp.eval(&coeffs_a[k], q1)).collect();
            let mut cands = if is_zero_poly(&at) {
                let coeffs = conic_in_q2(&cs[ib]);
                let v: Vec<u64> = (0..3).map(|k| fp.eval(&coeffs[k], q1)).collect();
                poly_roots_fp_counted(&v, p, seed ^ 1, &mut fp.muls)
            } else {
                poly_roots_fp_counted(&at, p, seed ^ 2, &mut fp.muls)
            };
            cands.sort();
            cands.dedup();
            for q2 in cands {
                if cs.iter().all(|c| conic_value(fp, c, q1, q2) == 0) && !out.contains(&(q1, q2)) {
                    out.push((q1, q2));
                }
            }
        }
        return out;
    }
    out
}

/// Statistics of the `D₃` solver.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct FghrStats {
    pub calls: u64,
    pub fp_muls: u64,
    /// Weil restriction (`F_{p³}` multiplications, at 11 `F_p` each).
    pub restrict_muls: u64,
    pub resultant_evaluations: u64,
    /// Calls whose degree-16 interpolation failed the check point and
    /// fell back to the coarse bound.
    pub degree_check_failures: u64,
    pub resultant_degree_sum: u64,
    pub resultant_degree_max: u64,
    /// Calls whose resultant vanished identically or whose minor
    /// vanished everywhere tried: skipped and reported, never guessed.
    pub unsolved: u64,
    pub q3_roots: u64,
    /// `(q₁, q₂, q₃)` solutions in `F_p³`.
    pub solutions: u64,
    /// Of those, the ones whose `t_i` are all in `F_p^*`.
    pub split_solutions: u64,
    pub unliftable: u64,
}

/// Solve the `D₃` system at one target abscissa: the `t`-triples
/// `(t₁, t₂, t₃) ∈ (F_p^*)³` up to even sign changes.
pub fn solve_d3(
    polys: &FghrPolynomials,
    x_r: Fp3El,
    seed: u64,
    stats: &mut FghrStats,
) -> Option<Vec<[u64; 3]>> {
    let f = &polys.curve3.field;
    let p = f.p;
    let mut fp = Fp { p, muls: 0 };
    stats.calls += 1;
    // Weil restriction at x_R.
    let xr = E3(x_r);
    let mut pow = [SolverField::ONE; 5];
    for i in 1..5 {
        pow[i] = f.mul(&pow[i - 1], &xr);
    }
    let mut comps: [HashMap<[u8; 3], u64>; 3] = Default::default();
    for (e, c) in &polys.d3_on_line {
        let v = f.mul(c, &pow[e[3] as usize]);
        for k in 0..3 {
            let entry = comps[k].entry([e[0], e[1], e[2]]).or_insert(0);
            *entry = addm(*entry, v.0[k], p);
        }
    }
    let restrict = (polys.d3_on_line.len() as u64 + 4) * 11;
    stats.restrict_muls += restrict;
    fp.muls += restrict;
    let conics = [
        ConicPoly::from_terms(&comps[0]),
        ConicPoly::from_terms(&comps[1]),
        ConicPoly::from_terms(&comps[2]),
    ];
    // R(q₃) has degree 16 = 64 / 4 (the weighted Bézout number of three
    // conics whose q₃-weights are (2, 2, 1)): interpolate through 17
    // points where Macaulay's minor does not vanish and confirm on an
    // 18th.  A failed check falls back to the coarse bound
    // 4·Σ deg_{q₃} (48), so the shortcut can cost time but never a root.
    let coarse = 4 * conics.iter().map(|c| c.degree_in_q3()).sum::<usize>();
    let mut z = 1u64;
    let sample = |fp: &mut Fp,
                  want: usize,
                  xs: &mut Vec<u64>,
                  ys: &mut Vec<u64>,
                  z: &mut u64,
                  stats: &mut FghrStats| {
        while xs.len() < want && *z < p {
            let cs = [
                conics[0].at(fp, *z),
                conics[1].at(fp, *z),
                conics[2].at(fp, *z),
            ];
            stats.resultant_evaluations += 1;
            if let Some(v) = conic_resultant(fp, &cs) {
                xs.push(*z);
                ys.push(v);
            }
            *z += 1;
        }
    };
    let (mut xs, mut ys) = (
        Vec::with_capacity(coarse + 1),
        Vec::with_capacity(coarse + 1),
    );
    sample(&mut fp, 18, &mut xs, &mut ys, &mut z, stats);
    let mut r = if xs.len() == 18 {
        let r17 = interpolate(&mut fp, &xs[..17], &ys[..17]);
        if fp.eval(&r17, xs[17]) == ys[17] {
            Some(r17)
        } else {
            None
        }
    } else {
        None
    };
    if r.is_none() {
        stats.degree_check_failures += 1;
        sample(&mut fp, coarse + 1, &mut xs, &mut ys, &mut z, stats);
        if xs.len() <= coarse {
            stats.unsolved += 1;
            stats.fp_muls += fp.muls;
            return None;
        }
        r = Some(interpolate(&mut fp, &xs, &ys));
    }
    let r = r.expect("set above");
    if r.is_empty() {
        stats.unsolved += 1;
        stats.fp_muls += fp.muls;
        return None;
    }
    let deg = (r.len() - 1) as u64;
    stats.resultant_degree_sum += deg;
    stats.resultant_degree_max = stats.resultant_degree_max.max(deg);
    let roots = poly_roots_fp_counted(&r, p, seed, &mut fp.muls);
    let mut out = Vec::new();
    for q3 in roots {
        if q3 == 0 {
            continue;
        }
        stats.q3_roots += 1;
        let cs = [
            conics[0].at(&mut fp, q3),
            conics[1].at(&mut fp, q3),
            conics[2].at(&mut fp, q3),
        ];
        for (q1, q2) in solve_three_conics(&mut fp, &cs, seed ^ q3) {
            stats.solutions += 1;
            // t_i² are the roots of Z³ − q₁Z² + q₂Z − q₃².
            let q3sq = fp.mul(q3, q3);
            let cubic = vec![fp.sub(0, q3sq), q2, fp.sub(0, q1), 1];
            let distinct = poly_roots_fp_counted(&cubic, p, seed ^ 3, &mut fp.muls);
            let mut zs = Vec::new();
            let mut rest = cubic.clone();
            for &root in &distinct {
                while rest.len() > 1 && fp.eval(&rest, root) == 0 {
                    // rest /= (Z − root), synthetic division.
                    let n = rest.len() - 1;
                    let mut q = vec![0u64; n];
                    q[n - 1] = rest[n];
                    for k in (1..n).rev() {
                        let t = fp.mul(q[k], root);
                        q[k - 1] = fp.add(rest[k], t);
                    }
                    rest = q;
                    zs.push(root);
                }
            }
            if zs.len() != 3 || zs.iter().any(|&z| z == 0 || jacobi(z, p) != 1) {
                continue;
            }
            let mut t: Vec<u64> = zs
                .iter()
                .map(|&z| sqrt_mod(z, p).expect("a square"))
                .collect();
            fp.muls += 3 * (64 - p.leading_zeros() as u64);
            let t01 = fp.mul(t[0], t[1]);
            let prod = fp.mul(t01, t[2]);
            if prod == q3 {
            } else if prod == fp.sub(0, q3) {
                t[2] = fp.sub(0, t[2]);
            } else {
                continue;
            }
            stats.split_solutions += 1;
            out.push([t[0], t[1], t[2]]);
        }
    }
    stats.fp_muls += fp.muls;
    Some(out)
}

// ── The oracles ────────────────────────────────────────────────────

fn lift_triples(
    inst: &dyn YLine,
    ctx: &InstanceCtx<'_, Fp3Curve>,
    fb: &FactorBase<Fp3Point>,
    ops: &mut GroupOps,
    counters: &mut OracleCounters,
    target: Fp3Point,
    triples: &[[u64; 3]],
) -> Option<Vec<usize>> {
    let f = &ctx.group.f;
    for t in triples {
        // A solution at the pole has a summand at infinity: not a
        // decomposition into three base points.
        let (Some(x1), Some(x2), Some(x3)) = (inst.x_at(t[0]), inst.x_at(t[1]), inst.x_at(t[2]))
        else {
            counters.lift_failures += 1;
            continue;
        };
        let xs = [f.pack(x1), f.pack(x2), f.pack(x3)];
        if let Some(idx) = lift_abscissae(ctx.group, fb, ops, &xs, target) {
            return Some(idx);
        }
        counters.lift_failures += 1;
    }
    None
}

/// The `D₃`-symmetrised oracle: 16 solutions a system.
pub struct FghrOracle<'i> {
    inst: &'i dyn YLine,
    polys: &'i FghrPolynomials,
    pub stats: FghrStats,
    totals: SolverTotals,
    seed: u64,
}

impl<'i> FghrOracle<'i> {
    pub fn new(inst: &'i dyn YLine, polys: &'i FghrPolynomials, seed: u64) -> Self {
        Self {
            inst,
            polys,
            stats: FghrStats::default(),
            totals: SolverTotals {
                solver: "d3-conic-resultant".into(),
                ..Default::default()
            },
            seed,
        }
    }
}

impl DecompositionOracle<Fp3Curve> for FghrOracle<'_> {
    fn name(&self) -> &str {
        "fghr-d3"
    }
    fn summands(&self) -> u32 {
        3
    }
    fn describe(&self, _params: &Params) -> String {
        "S₄ in Y = (x − x0 − c)/(x − x0 + c), symmetrised by the even sign changes and permutations of (Y₁, Y₂, Y₃), Weil-descended on the Y-line to three conics in (q₁, q₂) over F_p[q₃], solved by Macaulay's resultant in q₃".into()
    }
    fn decompose(
        &mut self,
        ctx: &InstanceCtx<'_, Fp3Curve>,
        fb: &FactorBase<Fp3Point>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: Fp3Point,
    ) -> Option<Vec<usize>> {
        if point.infinity {
            return None;
        }
        let started = Instant::now();
        let before = self.stats.fp_muls;
        let seed = self.seed ^ ctx.group.key(&point);
        let triples = solve_d3(self.polys, point.x, seed, &mut self.stats);
        let cost = SolverCost {
            ops: self.stats.fp_muls - before,
            op_unit: "F_p multiplications".into(),
            wall_ns: started.elapsed().as_nanos() as u64,
            peak_bytes: 0,
            degree_reached: None,
            solving_degree: None,
            timed_out: false,
            extra: Default::default(),
        };
        let shape = SystemShape {
            n_vars: 3,
            n_equations: 3,
            degrees: vec![4; 3],
            semi_regular_degree: None,
        };
        self.totals.absorb(&shape, Some(&cost), triples.is_none());
        let triples = triples?;
        if triples.is_empty() {
            return None;
        }
        // The solver returns one member of each class of even sign
        // changes, always the same one; an arm that keeps P and P + T in
        // different columns needs all four, so pick one per target from
        // the target's key (all four are decompositions of the target).
        let p = ctx.group.f.p;
        let flip = |t: u64| if t == 0 { 0 } else { p - t };
        let variant = (seed.rotate_right(17) ^ seed) % 4;
        let triples: Vec<[u64; 3]> = triples
            .into_iter()
            .map(|[a, b, c]| match variant {
                1 => [flip(a), flip(b), c],
                2 => [flip(a), b, flip(c)],
                3 => [a, flip(b), flip(c)],
                _ => [a, b, c],
            })
            .collect();
        let found = lift_triples(self.inst, ctx, fb, ops, counters, point, &triples);
        if found.is_none() {
            self.stats.unliftable += 1;
            counters.unliftable_systems += 1;
        }
        found
    }
    fn solver_totals(&self) -> Option<SolverTotals> {
        Some(self.totals.clone())
    }
}

/// The `S₃`-symmetrised oracle on the same `Y`-line: 64 solutions a
/// system, the existing Macaulay solver — the E11 oracle in the `Y`
/// coordinate, so the two presentations meet on one base.
pub struct YLineS4Oracle<'i> {
    inst: &'i dyn YLine,
    solver_inst: Instance3,
    polys: &'i FghrPolynomials,
    rng: StdRng,
    pub stats: SolveStats,
    pub unliftable: u64,
    totals: SolverTotals,
}

impl<'i> YLineS4Oracle<'i> {
    pub fn new(inst: &'i dyn YLine, polys: &'i FghrPolynomials, seed: u64) -> Self {
        Self {
            inst,
            solver_inst: Instance3 {
                curve: polys.curve3.clone(),
                q: Pt3::INFINITY,
                d: 0,
            },
            polys,
            rng: StdRng::seed_from_u64(seed ^ 0x5933),
            stats: SolveStats::default(),
            unliftable: 0,
            totals: SolverTotals {
                solver: "y-line-s4-macaulay".into(),
                ..Default::default()
            },
        }
    }
}

impl DecompositionOracle<Fp3Curve> for YLineS4Oracle<'_> {
    fn name(&self) -> &str {
        "y-line-s3-macaulay"
    }
    fn summands(&self) -> u32 {
        3
    }
    fn describe(&self, _params: &Params) -> String {
        "S₄ in Y symmetrised by the permutations of (Y₁, Y₂, Y₃) only, Weil-descended on the Y-line, solved by the Macaulay matrix and the eigenvalues of e₁".into()
    }
    fn decompose(
        &mut self,
        ctx: &InstanceCtx<'_, Fp3Curve>,
        fb: &FactorBase<Fp3Point>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: Fp3Point,
    ) -> Option<Vec<usize>> {
        if point.infinity {
            return None;
        }
        let started = Instant::now();
        let before = self.stats.fp_muls;
        let triples = solve_s4_subspace_with(
            &self.solver_inst,
            &self.polys.s3_on_line,
            &E3(point.x),
            &mut self.rng,
            &mut self.stats,
            SolveMode::default(),
        );
        let cost = SolverCost {
            ops: self.stats.fp_muls - before,
            op_unit: "F_p multiplications".into(),
            wall_ns: started.elapsed().as_nanos() as u64,
            peak_bytes: 0,
            degree_reached: Some(13),
            solving_degree: Some(10),
            timed_out: false,
            extra: Default::default(),
        };
        let shape = SystemShape {
            n_vars: 3,
            n_equations: 3,
            degrees: vec![4; 3],
            semi_regular_degree: None,
        };
        self.totals.absorb(&shape, Some(&cost), triples.is_none());
        let triples = triples?;
        if triples.is_empty() {
            return None;
        }
        let found = lift_triples(self.inst, ctx, fb, ops, counters, point, &triples);
        if found.is_none() {
            self.unliftable += 1;
            counters.unliftable_systems += 1;
        }
        found
    }
    fn solver_totals(&self) -> Option<SolverTotals> {
        Some(self.totals.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_framework::plugins::MitmOracle;

    fn ctx_of(inst: &FghrInstance) -> InstanceCtx<'_, Fp3Curve> {
        InstanceCtx {
            group: &inst.sub.curve,
            generator: inst.sub.generator,
            target: inst.sub.generator,
            r: inst.sub.r,
            cofactor: inst.sub.cofactor,
            group_order: inst.sub.group_order,
            name: inst.sub.name.clone(),
            field_degree: Some(3),
        }
    }

    #[test]
    fn the_y_line_base_folds_twelve_to_a_column() {
        let inst = generate_fghr_instance(7, 1, 8).unwrap();
        let (f12, r12) = fghr_line_base(&inst, FghrFold::FrobeniusTranslation).unwrap();
        let (f6, r6) = fghr_line_base(&inst, FghrFold::Frobenius).unwrap();
        let (f2, _) = fghr_line_base(&inst, FghrFold::Negation).unwrap();
        assert_eq!(f12.points.len(), f2.points.len());
        assert!((r12.points_per_orbit - 12.0).abs() < 1e-9, "{r12:?}");
        assert!((r6.points_per_orbit - 6.0).abs() < 1e-9, "{r6:?}");
        assert_eq!(f6.columns, 2 * f12.columns);
    }

    #[test]
    fn the_s4_in_y_is_invariant_under_even_sign_changes() {
        let inst = generate_fghr_instance(7, 2, 8).unwrap();
        let polys = fghr_polynomials(&inst).unwrap();
        // The D₃ rewrite would have refused a non-invariant polynomial;
        // check the sizes too: the D₃ form is the smallest.
        assert!(
            polys.d3_terms < polys.s3_terms,
            "{} vs {}",
            polys.d3_terms,
            polys.s3_terms
        );
        assert!(polys.s3_terms < polys.y_terms);
    }

    /// Both algebraic oracles agree with the pair table on whether a
    /// target decomposes, and every returned triple sums to its target.
    #[test]
    fn the_d3_and_s3_oracles_agree_with_the_pair_table() {
        let inst = generate_fghr_instance(7, 3, 8).unwrap();
        let (fb, _) = fghr_line_base(&inst, FghrFold::FrobeniusTranslation).unwrap();
        let polys = fghr_polynomials(&inst).unwrap();
        let ctx = ctx_of(&inst);
        let mut ops = GroupOps::default();
        let mut d3 = FghrOracle::new(&inst, &polys, 1);
        let mut s3 = YLineS4Oracle::new(&inst, &polys, 1);
        let mut mitm = MitmOracle::new(3);
        let mut params = Params::default();
        params.set("negation_folded", "1");
        mitm.prepare(&ctx, &fb, &params, &mut ops).unwrap();
        let mut c = [
            OracleCounters::default(),
            OracleCounters::default(),
            OracleCounters::default(),
        ];
        let (mut hits, mut d3_dis, mut s3_dis) = (0, 0, 0);
        for k in 2..120u64 {
            let pt = inst.sub.curve.mul(&mut ops, inst.sub.generator, k);
            let a = d3.decompose(&ctx, &fb, &mut ops, &mut c[0], pt);
            let b = s3.decompose(&ctx, &fb, &mut ops, &mut c[1], pt);
            let m = mitm.decompose(&ctx, &fb, &mut ops, &mut c[2], pt);
            for idx in [&a, &b].into_iter().flatten() {
                let sum = idx.iter().fold(inst.sub.curve.identity(), |acc, &i| {
                    inst.sub.curve.add(&mut ops, acc, fb.points[i])
                });
                assert_eq!(sum, pt, "k = {k}");
            }
            if m.is_some() {
                hits += 1;
            }
            if a.is_some() != m.is_some() {
                d3_dis += 1;
            }
            if b.is_some() != m.is_some() {
                s3_dis += 1;
            }
        }
        assert!(hits > 5, "{hits}");
        assert!(
            d3_dis <= 2,
            "D₃ disagrees on {d3_dis} targets; {:?}",
            d3.stats
        );
        assert!(s3_dis <= 4, "S₃ disagrees on {s3_dis} targets");
        assert!(
            d3.stats.resultant_degree_max <= 24,
            "resultant degree {}",
            d3.stats.resultant_degree_max
        );
    }
}

#[cfg(test)]
mod rank_tests {
    use super::*;
    use crate::cryptanalysis::glv_invariant_experiments::full_rank_stream;

    /// The `D₃` solver returns one member of each class of even sign
    /// changes; picked per target, every member is used, so an arm that
    /// keeps `P` and `P + T` apart (`⟨−1, π⟩`) still reaches full rank.
    /// With a fixed member it stopped at 17 of 21 on this instance.
    #[test]
    fn the_d3_oracle_lets_the_six_fold_arm_reach_full_rank() {
        let inst = generate_fghr_instance(7, 1, 8).unwrap();
        let sub = &inst.sub;
        let (f6, _) = fghr_line_base(&inst, FghrFold::Frobenius).unwrap();
        let (c2, _) = fghr_line_base(&inst, FghrFold::Negation).unwrap();
        let polys = fghr_polynomials(&inst).unwrap();
        let mut ops = GroupOps::default();
        let target = sub.curve.mul(&mut ops, sub.generator, 1234);
        let ctx = InstanceCtx {
            group: &sub.curve,
            generator: sub.generator,
            target,
            r: sub.r,
            cofactor: sub.cofactor,
            group_order: sub.group_order,
            name: sub.name.clone(),
            field_degree: Some(3),
        };
        let mut d3 = FghrOracle::new(&inst, &polys, 1);
        let rep = full_rank_stream(
            &sub.curve,
            sub.generator,
            target,
            sub.r,
            sub.cofactor,
            1234,
            &f6,
            &c2,
            1,
            10_000_000,
            3,
            None,
            |o, c, pt| d3.decompose(&ctx, &f6, o, c, pt),
        )
        .unwrap();
        assert_eq!(rep.folded.rank_final as usize, rep.folded.columns + 1);
        assert!(rep.folded.verified && rep.control.verified);
    }
}
