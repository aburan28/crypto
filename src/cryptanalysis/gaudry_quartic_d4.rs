//! **Gaudry's `k = 4` solve with FGHR's 2-torsion symmetry: the `D₄`
//! invariants of `S₅` (candidate C8 of
//! `research/pdp_root_finder_requirements_20261012` §8).**
//!
//! [`super::gaudry_quartic`] measures `C₄`, one decomposition of a point of
//! `E(F_{p⁴})` into four points with abscissae in `F_p`: `S₅` symmetrised by
//! `S₄` in the elementary symmetric functions `e₁, …, e₄` of the abscissae,
//! Weil-restricted to four equations of degree 8 over `F_p`, Bézout number
//! `8⁴ = 4096`, Macaulay degree 29 (`1.227·10¹²` `F_p` multiplications).
//!
//! Here the curve has a rational 2-torsion point `T = (x₀, 0)` with
//! `f'(x₀) = c²`, `c ∈ F_p^*`, and the base is FGHR's line in the
//! coordinate `Y = (x − x₀ − c)/(x − x₀ + c)` (as in [`super::fghr_full`]
//! at `k = 3`): `Y(P + T) = −Y(P)`, `Y(−P) = Y(P)`, base `{P : Y(P) ∈ F_p}`.
//! Since `2T = 0`, translating an even number of the four summands by `T`
//! keeps their sum, so `H(Y₁, …, Y₄; Y_R) = S₅(x(Y₁), …, x(Y_R)) ·
//! Π(1 − Y_i)⁸ (1 − Y_R)⁸` is invariant under permutations and even sign
//! changes of `(Y₁, …, Y₄)`: the Weyl group `D₄ = (Z/2)³ ⋊ S₄`, of order
//! 192 against `S₄`'s 24.  Its invariant ring is polynomial in
//! `σ_k = e_k(Y₁², …, Y₄²)` (`k = 1, 2, 3`) and `π = Y₁Y₂Y₃Y₄`; with
//! `σ₄ = π²` adjoined as a fifth unknown, every monomial of `H` (degree
//! `≤ 8` in each `Y_i`) has total degree `≤ 4` in `(σ₁, σ₂, σ₃, σ₄, π)`.
//! The system is the four Weil components of degree 4 and `π² − σ₄`:
//! Bézout `4⁴ · 2 = 512`, Macaulay degree `4·3 + 1 + 1 = 14`.
//!
//! `H`'s coefficients on the 105 invariant monomials are found by
//! interpolation, as `gaudry_quartic` finds the `S₄` form, and checked on
//! fresh points: a failure of the invariance would show there.  Every solve
//! can be checked against an independent meet-in-the-middle enumeration of
//! the decompositions over the same base ([`mitm_y_decompositions`]).

use std::collections::{BTreeSet, HashMap};

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use super::gaudry_cubic::rref_mod_p;
use super::gaudry_quartic::{
    solve_system_w_inner, split_quartic, Curve4, Fp4, Pt4, QuarticSolve, ReplayOutcome, ShapeTrace,
    E4,
};
use super::ic_boundary::{CountedGroup, FactorBase, GroupOps, OracleCounters};
use super::ic_framework::stages::{
    DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, SolverTotals, SystemShape,
};
use super::residual_walk::{is_prime_u64, pow_mod, sqrt_mod};

/// The interpolation nodes for `Y_R` (`Y = 1` is the pole of `x(Y)`).
const NODES: [u64; 9] = [2, 3, 4, 5, 6, 7, 8, 9, 10];

/// A curve over `F_{p⁴}` with a rational 2-torsion point `T = (x₀, 0)`,
/// `f'(x₀) = c²`.
pub struct D4Instance {
    pub curve: Curve4,
    pub x0: E4,
    pub c: u64,
}

/// `x₀ ∈ F_{p⁴} \ F_p` and `c ∈ F_p^*` at random, `a = c² − 3x₀²`,
/// `b = −x₀³ − a x₀`; coefficients outside `F_{p²}` as in
/// [`Curve4::random`].  Deterministic in `seed`.
pub fn generate_d4_instance(p: u64, seed: u64) -> D4Instance {
    assert!(p % 4 == 1 && p < 1 << 15, "p ≡ 1 (mod 4), below 2^15");
    let mut rng = StdRng::seed_from_u64(seed ^ 0xd4d4);
    let f = Fp4::new(p);
    loop {
        let x0 = f.random(&mut rng);
        if x0.in_fp() {
            continue;
        }
        let c = rng.gen_range(1..p);
        let a = f.sub(&f.from_fp(c * c % p), &f.scale(&f.sq(&x0), 3));
        let b = f.neg(&f.add(&f.mul(&f.sq(&x0), &x0), &f.mul(&a, &x0)));
        if a.0[1] == 0 && a.0[3] == 0 || b.0[1] == 0 && b.0[3] == 0 {
            continue;
        }
        let a3 = f.mul(&f.sq(&a), &a);
        if f.add(&f.scale(&a3, 4), &f.scale(&f.sq(&b), 27)).is_zero() {
            continue;
        }
        let curve = Curve4::with_params(p, a, b);
        assert!(curve.rhs(&x0).is_zero(), "x₀ is a root of x³ + ax + b");
        return D4Instance { curve, x0, c };
    }
}

impl D4Instance {
    /// `x(Y) = x₀ + c(1 + Y)/(1 − Y)`.
    pub fn x_of_y(&self, y: &E4) -> E4 {
        let f = &self.curve.f;
        let c = f.from_fp(self.c);
        let num = f.mul(&c, &f.add(&E4::ONE, y));
        let den = f.sub(&E4::ONE, y);
        f.add(&self.x0, &f.mul(&num, &f.inv(&den)))
    }

    /// `Y(x) = (x − x₀ − c)/(x − x₀ + c)`.
    pub fn y_of_x(&self, x: &E4) -> E4 {
        let f = &self.curve.f;
        let c = f.from_fp(self.c);
        let u = f.sub(x, &self.x0);
        f.mul(&f.sub(&u, &c), &f.inv(&f.add(&u, &c)))
    }

    /// `H(Y₁, …, Y₄; Y_R)`, denominators cleared.
    fn h(&self, ys: &[E4; 4], yr: &E4) -> E4 {
        let f = &self.curve.f;
        let x = [
            self.x_of_y(&ys[0]),
            self.x_of_y(&ys[1]),
            self.x_of_y(&ys[2]),
            self.x_of_y(&ys[3]),
            self.x_of_y(yr),
        ];
        let mut v = self.curve.s5(&x);
        for y in ys.iter().chain(std::iter::once(yr)) {
            let d = f.sub(&E4::ONE, y);
            let d2 = f.sq(&d);
            let d4 = f.sq(&d2);
            v = f.mul(&v, &f.sq(&d4));
        }
        v
    }
}

/// The invariant monomials `σ₁^a σ₂^b σ₃^c π^d` of degree `≤ 8` in each
/// `Y_i`: `2(a + b + c) + d ≤ 8`, 105 of them.
pub fn d4_monomials() -> Vec<[u8; 4]> {
    let mut v = Vec::new();
    for d in 0..=8u8 {
        let k = (8 - d) / 2;
        for a in 0..=k {
            for b in 0..=(k - a) {
                for c in 0..=(k - a - b) {
                    v.push([a, b, c, d]);
                }
            }
        }
    }
    v
}

/// `(σ₁, σ₂, σ₃, π)` of `(Y₁, …, Y₄)`.
fn invariants(f: &Fp4, y: &[E4; 4]) -> [E4; 4] {
    let s: Vec<E4> = y.iter().map(|v| f.sq(v)).collect();
    let mut e1 = E4::ZERO;
    let mut e2 = E4::ZERO;
    let mut e3 = E4::ZERO;
    for i in 0..4 {
        e1 = f.add(&e1, &s[i]);
        for j in (i + 1)..4 {
            let sij = f.mul(&s[i], &s[j]);
            e2 = f.add(&e2, &sij);
            for k in (j + 1)..4 {
                e3 = f.add(&e3, &f.mul(&sij, &s[k]));
            }
        }
    }
    let pi = f.mul(&f.mul(&y[0], &y[1]), &f.mul(&y[2], &y[3]));
    [e1, e2, e3, pi]
}

fn mono_values(f: &Fp4, inv: &[E4; 4], monos: &[[u8; 4]]) -> Vec<E4> {
    let mut pows = [[E4::ONE; 9]; 4];
    for i in 0..4 {
        for k in 1..9 {
            pows[i][k] = f.mul(&pows[i][k - 1], &inv[i]);
        }
    }
    monos
        .iter()
        .map(|m| {
            f.mul(
                &f.mul(&pows[0][m[0] as usize], &pows[1][m[1] as usize]),
                &f.mul(&pows[2][m[2] as usize], &pows[3][m[3] as usize]),
            )
        })
        .collect()
}

/// `H = Σ_m coef[m](Y_R) · m(σ₁, σ₂, σ₃, π)`, nine `F_{p⁴}` coefficients in
/// `Y_R` per invariant monomial.
pub struct SymmetrisedD4 {
    pub monos: Vec<[u8; 4]>,
    pub coef: Vec<[E4; 9]>,
    /// Fresh points on which the interpolated form reproduced `H`.
    pub checked_points: usize,
    pub precompute_muls: u64,
}

impl SymmetrisedD4 {
    fn coefficients_in_last(inst: &D4Instance, y: &[E4; 4], vinv: &[[u64; 9]; 9]) -> [E4; 9] {
        let f = &inst.curve.f;
        let vals: Vec<E4> = NODES.iter().map(|&k| inst.h(y, &f.from_fp(k))).collect();
        core::array::from_fn(|j| {
            vals.iter()
                .enumerate()
                .fold(E4::ZERO, |acc, (k, v)| f.add(&acc, &f.scale(v, vinv[j][k])))
        })
    }

    fn random_y(f: &Fp4, rng: &mut StdRng) -> [E4; 4] {
        core::array::from_fn(|_| loop {
            let v = f.random(rng);
            if v != E4::ONE {
                break v;
            }
        })
    }

    /// Interpolate on 105 random points of `F_{p⁴}⁴` and check on
    /// `checks` fresh ones.  `Err` if the check fails, i.e. `H` is not in
    /// the span of the `D₄` invariants.
    pub fn precompute(
        inst: &D4Instance,
        rng: &mut StdRng,
        checks: usize,
    ) -> Result<SymmetrisedD4, String> {
        let f = &inst.curve.f;
        let p = f.p;
        let before = f.muls();
        let vinv: [[u64; 9]; 9] = {
            let mut m: Vec<Vec<u64>> = NODES
                .iter()
                .enumerate()
                .map(|(i, &k)| {
                    let mut row: Vec<u64> = (0..9).map(|j| pow_mod(k, j, p)).collect();
                    row.extend((0..9).map(|c| u64::from(c == i)));
                    row
                })
                .collect();
            let mut scratch = 0u64;
            rref_mod_p(&mut m, p, &mut scratch);
            core::array::from_fn(|j| core::array::from_fn(|k| m[j][9 + k]))
        };
        let monos = d4_monomials();
        let n = monos.len();
        loop {
            let mut sys: Vec<Vec<E4>> = Vec::with_capacity(n);
            for _ in 0..n {
                let y = Self::random_y(f, rng);
                let mut row = mono_values(f, &invariants(f, &y), &monos);
                row.extend(Self::coefficients_in_last(inst, &y, &vinv));
                sys.push(row);
            }
            let mut ok = true;
            for c in 0..n {
                let Some(r) = (c..n).find(|&r| !sys[r][c].is_zero()) else {
                    ok = false;
                    break;
                };
                sys.swap(r, c);
                let inv = f.inv(&sys[c][c]);
                let pivot: Vec<E4> = sys[c].iter().map(|v| f.mul(v, &inv)).collect();
                sys[c] = pivot.clone();
                for r in 0..n {
                    if r == c || sys[r][c].is_zero() {
                        continue;
                    }
                    let k = sys[r][c];
                    for j in c..pivot.len() {
                        if !pivot[j].is_zero() {
                            let v = f.mul(&k, &pivot[j]);
                            sys[r][j] = f.sub(&sys[r][j], &v);
                        }
                    }
                }
            }
            if !ok {
                continue;
            }
            let coef: Vec<[E4; 9]> = (0..n)
                .map(|r| core::array::from_fn(|j| sys[r][n + j]))
                .collect();
            for t in 0..checks {
                let y = Self::random_y(f, rng);
                let want = Self::coefficients_in_last(inst, &y, &vinv);
                let mv = mono_values(f, &invariants(f, &y), &monos);
                for j in 0..9 {
                    let got = mv
                        .iter()
                        .zip(&coef)
                        .fold(E4::ZERO, |acc, (m, c)| f.add(&acc, &f.mul(m, &c[j])));
                    if got != want[j] {
                        return Err(format!(
                            "H is not in the span of the D₄ invariants (check point {t}, coefficient Y_R^{j})"
                        ));
                    }
                }
            }
            return Ok(SymmetrisedD4 {
                monos,
                coef,
                checked_points: checks,
                precompute_muls: f.muls() - before,
            });
        }
    }

    /// The five equations over `F_p` in `(σ₁, σ₂, σ₃, σ₄, π)`: the four
    /// components of `H(·; Y_R)` with `π^d = σ₄^{⌊d/2⌋} π^{d mod 2}`, and
    /// `π² − σ₄`.
    pub fn system(&self, f: &Fp4, y_r: &E4) -> Vec<HashMap<[u8; 5], u64>> {
        let mut pows = [E4::ONE; 9];
        for k in 1..9 {
            pows[k] = f.mul(&pows[k - 1], y_r);
        }
        let mut out: Vec<HashMap<[u8; 5], u64>> = vec![HashMap::new(); 4];
        for (m, c) in self.monos.iter().zip(&self.coef) {
            let v = (0..9).fold(E4::ZERO, |acc, j| f.add(&acc, &f.mul(&c[j], &pows[j])));
            let key = [m[0], m[1], m[2], m[3] / 2, m[3] % 2];
            for (i, comp) in out.iter_mut().enumerate() {
                if v.0[i] != 0 {
                    comp.insert(key, v.0[i]);
                }
            }
        }
        let mut rel = HashMap::new();
        rel.insert([0, 0, 0, 0, 2], 1);
        rel.insert([0, 0, 0, 1, 0], f.p - 1);
        out.push(rel);
        out
    }
}

impl SymmetrisedD4 {
    /// The four Weil components of `H(·; Y_R)` in `(σ₁, σ₂, σ₃, π)`, with
    /// `σ₄ = π²` substituted: `σ₁^a σ₂^b σ₃^c π^d`, `2(a + b + c) + d ≤ 8`,
    /// weighted homogeneous of degree 8 at most for the weights
    /// [`WEIGHTS`] = `(2, 2, 2, 1)` (the degree in each `Y_i`).  The
    /// same Weil restriction as [`Self::system`], one equation fewer.
    pub fn system_weighted(&self, f: &Fp4, y_r: &E4) -> Vec<HashMap<[u8; 4], u64>> {
        let mut pows = [E4::ONE; 9];
        for k in 1..9 {
            pows[k] = f.mul(&pows[k - 1], y_r);
        }
        let mut out: Vec<HashMap<[u8; 4], u64>> = vec![HashMap::new(); 4];
        for (m, c) in self.monos.iter().zip(&self.coef) {
            let v = (0..9).fold(E4::ZERO, |acc, j| f.add(&acc, &f.mul(&c[j], &pows[j])));
            for (i, comp) in out.iter_mut().enumerate() {
                if v.0[i] != 0 {
                    comp.insert(*m, v.0[i]);
                }
            }
        }
        out
    }
}

/// The weights of `(σ₁, σ₂, σ₃, π)`: each `σ_k` has degree 2 in every
/// `Y_i` it involves, `π` degree 1.  With them the substituted system's
/// Hilbert series (four equations of weighted degree 8, if regular) is
/// `(1 + t² + t⁴ + t⁶)³ (1 + t + … + t⁷)`: degree 25, value 512 at
/// `t = 1`, so weighted Macaulay degree 26, 4,200 columns against the
/// unweighted five-unknown matrix's 11,628.
pub const WEIGHTS: [u8; 4] = [2, 2, 2, 1];

/// How a `D₄` solve runs.  All off is the committed sessions' oracle,
/// operation for operation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct D4Options {
    /// Four unknowns `(σ₁, σ₂, σ₃, π)` with `σ₄ = π²` substituted, in the
    /// weights [`WEIGHTS`] (weighted Macaulay degree 26), instead of five
    /// unknowns and `π² − σ₄`, unweighted (degree 14).
    pub weighted: bool,
    /// The eigen step along the chains of the multiplication unknown
    /// (`π` when weighted, `σ₁` otherwise): only the border normal forms,
    /// the characteristic polynomial as `det A(λ)` of a matrix with one
    /// row per chain, eigenvectors as its kernels.
    pub chains: bool,
    /// The elimination's pivot in each column is the sparsest candidate.
    pub sparsest: bool,
}

/// What one `D₄` solve cost and found.
#[derive(Clone, Debug, Default, Serialize)]
pub struct D4Solve {
    pub solve: QuarticSolve,
    /// `F_p`-rational `(σ₁, σ₂, σ₃, σ₄, π)` of the system.
    pub rational_solutions: usize,
    /// Of those, the ones whose `Y_i` lie in `F_p` and give base points
    /// summing (with signs) to `R`: the decompositions, as sorted
    /// `Y`-quadruples.
    pub decompositions: Vec<[u64; 4]>,
    /// The first decomposition found, as four signed points summing to `R`.
    #[serde(skip)]
    pub witness: Option<[Pt4; 4]>,
}

/// The base points `{P : Y(P) ∈ F_p \ {±1}}`, one per `Y` (with `−P`
/// implied), and their `Y`.
pub fn y_base(inst: &D4Instance, rng: &mut StdRng) -> Vec<(u64, Pt4)> {
    let f = &inst.curve.f;
    let p = f.p;
    (0..p)
        .filter(|&y| y != 1 && y != p - 1)
        .filter_map(|y| {
            let x = inst.x_of_y(&f.from_fp(y));
            inst.curve.lift_x(&x, rng).map(|pt| (y, pt))
        })
        .collect()
}

/// The signs that make the four points sum to `r`, applied, if any do.
fn signed_sum(curve: &Curve4, pts: &[Pt4; 4], r: &Pt4) -> Option<[Pt4; 4]> {
    for s in 0..16u32 {
        let signed: [Pt4; 4] = core::array::from_fn(|i| {
            if s >> i & 1 == 1 {
                curve.neg(&pts[i])
            } else {
                pts[i]
            }
        });
        let mut acc = Pt4::INF;
        for q in &signed {
            acc = curve.add(&acc, q);
        }
        if !acc.inf && acc.x == r.x && acc.y == r.y {
            return Some(signed);
        }
    }
    None
}

/// One decomposition by the `D₄` system: Weil restriction, the Macaulay
/// solve at degree `macaulay` (14 for a generic system), the split of each
/// rational solution (`T⁴ − σ₁T³ + σ₂T² − σ₃T + σ₄` into `Y_i²`, square
/// roots with `ΠY_i = π`) and the check that the points sum to `R`.
/// `C₄` is `out.solve.total_muls` (Weil restriction, solve, split).
pub fn solve_d4(
    inst: &D4Instance,
    pre: &SymmetrisedD4,
    r: &Pt4,
    macaulay: u8,
    rng: &mut StdRng,
    debug: bool,
) -> Option<D4Solve> {
    solve_d4_or_cost(inst, pre, r, macaulay, rng, debug).ok()
}

/// [`solve_d4`], with a system the degree-`macaulay` matrix does not solve
/// (a degenerate one) returned as `Err` carrying what the attempt cost,
/// phases summed into `total_muls`: work done whether or not it solved.
#[allow(clippy::result_large_err)]
pub fn solve_d4_or_cost(
    inst: &D4Instance,
    pre: &SymmetrisedD4,
    r: &Pt4,
    macaulay: u8,
    rng: &mut StdRng,
    debug: bool,
) -> Result<D4Solve, QuarticSolve> {
    solve_d4_form(inst, pre, r, D4Options::default(), macaulay, rng, debug)
}

/// The shape of a first solved system, replayed on the later ones
/// (`gaudry_quartic`'s [`ShapeTrace`], F4Remake's idea for a Macaulay
/// matrix): the rows that reduced to zero are not built.  A replay whose
/// leads, row count, rank or dimension differ runs the full solve as well,
/// both charged; one that matches spans the same row space (a subset of
/// the rows of equal rank), so it is exact.
#[derive(Clone, Debug, Default)]
pub struct D4Trace {
    pub shape: Option<ShapeTrace>,
    pub replays: u64,
    pub fallbacks: u64,
    /// Why replays fell back, by the check that failed.
    pub fallback_reasons: std::collections::BTreeMap<String, u64>,
}

#[allow(clippy::too_many_arguments)]
fn solve_traced<const N: usize>(
    comps: &[HashMap<[u8; N], u64>],
    w: &[u8; N],
    opts: D4Options,
    macaulay: u8,
    p: u64,
    rng: &mut StdRng,
    st: &mut QuarticSolve,
    debug: bool,
    trace: Option<&mut D4Trace>,
) -> Option<Vec<[u64; N]>> {
    let run = |rng: &mut StdRng,
               st: &mut QuarticSolve,
               replay: Option<&ShapeTrace>,
               record: Option<&mut ShapeTrace>| {
        solve_system_w_inner::<N>(
            comps,
            w,
            opts.chains,
            opts.sparsest,
            true,
            0,
            macaulay,
            p,
            rng,
            st,
            debug,
            replay,
            record,
        )
    };
    let Some(tr) = trace else {
        return run(rng, st, None, None).ok().flatten();
    };
    if let Some(shape) = &tr.shape {
        match run(rng, st, Some(shape), None) {
            Ok(sols) => {
                tr.replays += 1;
                st.replay = ReplayOutcome::Replayed;
                return sols;
            }
            Err(()) => {
                tr.fallbacks += 1;
                *tr.fallback_reasons.entry(st.outcome.clone()).or_insert(0) += 1;
            }
        }
        let sols = run(rng, st, None, None).ok().flatten();
        st.replay = ReplayOutcome::FellBack;
        return sols;
    }
    let mut shape = ShapeTrace::default();
    let sols = run(rng, st, None, Some(&mut shape)).ok().flatten();
    if sols.is_some() {
        tr.shape = Some(shape);
    }
    sols
}

/// [`solve_d4_or_cost`] with [`D4Options`]: all off is that function,
/// operation for operation; `weighted` substitutes `σ₄ = π²` and solves the
/// four equations in the weighted grading [`WEIGHTS`] at weighted Macaulay
/// degree `macaulay` (26 for a regular system), the multiplication matrix
/// being that of `π`.
#[allow(clippy::result_large_err)]
pub fn solve_d4_form(
    inst: &D4Instance,
    pre: &SymmetrisedD4,
    r: &Pt4,
    opts: D4Options,
    macaulay: u8,
    rng: &mut StdRng,
    debug: bool,
) -> Result<D4Solve, QuarticSolve> {
    solve_d4_traced(inst, pre, r, opts, macaulay, rng, debug, None)
}

/// [`solve_d4_form`] with a [`D4Trace`] carried from solve to solve
/// (`None`: no replay, the same arithmetic as [`solve_d4_form`]).
#[allow(clippy::too_many_arguments)]
#[allow(clippy::result_large_err)]
pub fn solve_d4_traced(
    inst: &D4Instance,
    pre: &SymmetrisedD4,
    r: &Pt4,
    opts: D4Options,
    macaulay: u8,
    rng: &mut StdRng,
    debug: bool,
    trace: Option<&mut D4Trace>,
) -> Result<D4Solve, QuarticSolve> {
    let f = &inst.curve.f;
    let p = f.p;
    let mut out = D4Solve::default();
    let before = f.muls();
    let y_r = inst.y_of_x(&r.x);
    let sols: Option<Vec<[u64; 5]>> = if opts.weighted {
        let comps = pre.system_weighted(f, &y_r);
        out.solve.weil_muls = f.muls() - before;
        solve_traced::<4>(
            &comps,
            &WEIGHTS,
            opts,
            macaulay,
            p,
            rng,
            &mut out.solve,
            debug,
            trace,
        )
        .map(|v| {
            v.into_iter()
                .map(|s| {
                    // σ₄ = π², one multiplication.
                    out.solve.split_muls += 1;
                    [s[0], s[1], s[2], s[3] * s[3] % p, s[3]]
                })
                .collect()
        })
    } else {
        let comps = pre.system(f, &y_r);
        out.solve.weil_muls = f.muls() - before;
        solve_traced::<5>(
            &comps,
            &[1; 5],
            opts,
            macaulay,
            p,
            rng,
            &mut out.solve,
            debug,
            trace,
        )
    };
    let Some(sols) = sols else {
        if debug {
            eprintln!(
                "  unsolved: {} (dim {}, pivots {}, zero rows {})",
                out.solve.outcome, out.solve.dim, out.solve.pivots, out.solve.zero_rows
            );
        }
        out.solve.close_public();
        return Err(out.solve);
    };
    out.solve.equation_degree = 4;
    out.rational_solutions = sols.len();
    let mut found = BTreeSet::new();
    for s in sols {
        let mut sm = 0u64;
        let Some(sq) = split_quartic(&[s[0], s[1], s[2], s[3]], p, rng, &mut sm) else {
            out.solve.split_muls += sm;
            continue;
        };
        out.solve.split_muls += sm;
        let roots: Option<Vec<u64>> = sq.iter().map(|&v| sqrt_mod(v, p)).collect();
        let Some(roots) = roots else {
            continue;
        };
        for signs in 0..16u32 {
            let ys: [u64; 4] = core::array::from_fn(|i| {
                if signs >> i & 1 == 1 {
                    (p - roots[i]) % p
                } else {
                    roots[i]
                }
            });
            let prod = ys.iter().fold(1u64, |a, &v| a * v % p);
            if prod != s[4] {
                continue;
            }
            if ys.iter().any(|&v| v == 1 || v == p - 1) {
                continue;
            }
            let mut sorted = ys;
            sorted.sort_unstable();
            if sorted.windows(2).any(|w| w[0] == w[1]) {
                continue;
            }
            let pts: Option<Vec<Pt4>> = ys
                .iter()
                .map(|&v| inst.curve.lift_x(&inst.x_of_y(&f.from_fp(v)), rng))
                .collect();
            let Some(pts) = pts else {
                continue;
            };
            let pts: [Pt4; 4] = [pts[0], pts[1], pts[2], pts[3]];
            if let Some(signed) = signed_sum(&inst.curve, &pts, r) {
                found.insert(sorted);
                if out.witness.is_none() {
                    out.witness = Some(signed);
                }
            }
        }
    }
    out.decompositions = found.into_iter().collect();
    out.solve.close_public();
    Ok(out)
}

/// Every set of four distinct base points `{P_i, P_j, P_k, P_l}` with
/// `±P_i ± P_j ± P_k ± P_l = R`, as sorted `Y`-quadruples, by meet in the
/// middle over pairs: the independent check of [`solve_d4`].
pub fn mitm_y_decompositions(curve: &Curve4, base: &[(u64, Pt4)], r: &Pt4) -> BTreeSet<[u64; 4]> {
    let n = base.len();
    let mut table: HashMap<E4, Vec<(usize, usize)>> = HashMap::new();
    for i in 0..n {
        for j in (i + 1)..n {
            for q in [
                curve.add(&base[i].1, &base[j].1),
                curve.sub(&base[i].1, &base[j].1),
            ] {
                if !q.inf {
                    table.entry(q.x).or_default().push((i, j));
                }
            }
        }
    }
    let mut out = BTreeSet::new();
    for k in 0..n {
        for l in (k + 1)..n {
            for (sk, sl) in [(false, false), (false, true), (true, false), (true, true)] {
                let pk = if sk { curve.neg(&base[k].1) } else { base[k].1 };
                let pl = if sl { curve.neg(&base[l].1) } else { base[l].1 };
                let v = curve.sub(&curve.sub(r, &pk), &pl);
                if v.inf {
                    continue;
                }
                if let Some(pairs) = table.get(&v.x) {
                    for &(i, j) in pairs {
                        if i == k || i == l || j == k || j == l {
                            continue;
                        }
                        let mut ys = [base[i].0, base[j].0, base[k].0, base[l].0];
                        ys.sort_unstable();
                        out.insert(ys);
                    }
                }
            }
        }
    }
    out
}

// ── In ecbench: the curve, its counted group, the base and the oracle ─────

/// `x ↦ x₀ + p x₁ + p² x₂ + p³ x₃ < p⁴`, an element as one word.
pub fn pack_e4(p: u64, e: &E4) -> u64 {
    e.0[0] + p * (e.0[1] + p * (e.0[2] + p * e.0[3]))
}

/// The inverse of [`pack_e4`].
pub fn unpack_e4(p: u64, mut v: u64) -> E4 {
    let mut out = [0u64; 4];
    for o in &mut out {
        *o = v % p;
        v /= p;
    }
    E4(out)
}

/// `E(F_{p⁴})` as a counted group: every addition and doubling one unit
/// on the ledger.  A point's key is its packed abscissa and one bit for
/// the ordinate (`y` or `−y`, the larger packed value), so it identifies
/// the point exactly for `p < 2^15`.
impl CountedGroup for Curve4 {
    type Elt = Pt4;
    fn identity(&self) -> Pt4 {
        Pt4::INF
    }
    fn is_identity(&self, p: &Pt4) -> bool {
        p.inf
    }
    fn add(&self, ops: &mut GroupOps, p: Pt4, q: Pt4) -> Pt4 {
        if p.inf || q.inf || p.x != q.x {
            ops.adds += 1;
        } else {
            ops.doubles += 1;
        }
        Curve4::add(self, &p, &q)
    }
    fn double(&self, ops: &mut GroupOps, p: Pt4) -> Pt4 {
        ops.doubles += 1;
        Curve4::add(self, &p, &p)
    }
    fn neg(&self, p: Pt4) -> Pt4 {
        Curve4::neg(self, &p)
    }
    fn key(&self, p: &Pt4) -> u64 {
        if p.inf {
            return u64::MAX;
        }
        let q = self.f.p;
        let y = pack_e4(q, &p.y);
        let ny = pack_e4(q, &self.f.neg(&p.y));
        pack_e4(q, &p.x) * 2 + u64::from(y > ny)
    }
}

/// A curve of [`generate_d4_instance`]'s family with its group order
/// `#E(F_{p⁴}) = h·r`, `r` prime, and a generator of the order-`r`
/// subgroup.
pub struct D4Curve {
    pub inst: D4Instance,
    /// The instance seed [`generate_d4_instance`] was called with.
    pub instance_seed: u64,
    pub group_order: u64,
    pub r: u64,
    pub cofactor: u64,
    pub generator: Pt4,
}

impl D4Curve {
    pub fn p(&self) -> u64 {
        self.inst.curve.f.p
    }
    /// `T = (x₀, 0)`.
    pub fn torsion(&self) -> Pt4 {
        Pt4 {
            x: self.inst.x0,
            y: E4::ZERO,
            inf: false,
        }
    }
}

/// `#E(F_{p⁴})`: the one multiple in the Hasse interval
/// `[p⁴ + 1 − 2p², p⁴ + 1 + 2p²]` that kills every one of a few random
/// points, by a walk through the interval.  `None` if the candidates do
/// not narrow to one.
pub fn group_order_d4(curve: &Curve4, rng: &mut StdRng) -> Option<u64> {
    let p = curve.f.p;
    let q = p.pow(4);
    let (lo, hi) = (q + 1 - 2 * p * p, q + 1 + 2 * p * p);
    let mut cands: Option<Vec<u64>> = None;
    for _ in 0..24 {
        let pt = loop {
            let x = curve.f.random(rng);
            if let Some(pt) = curve.lift_x(&x, rng) {
                break pt;
            }
        };
        let mut scratch = GroupOps::default();
        let mut acc = CountedGroup::mul(curve, &mut scratch, pt, lo);
        let mut here = Vec::new();
        for m in lo..=hi {
            if acc.inf {
                here.push(m);
            }
            acc = Curve4::add(curve, &acc, &pt);
        }
        cands = Some(match cands {
            None => here,
            Some(c) => c.into_iter().filter(|m| here.contains(m)).collect(),
        });
        if cands.as_ref().is_some_and(|c| c.len() == 1) {
            return cands.map(|c| c[0]);
        }
    }
    None
}

/// The first instance seed from `seed` whose curve has `#E = h·r`, `r`
/// prime, `h ≤ max_cofactor` (even: `T` is rational), and a generator of
/// the order-`r` subgroup.  Deterministic in `(p, seed, max_cofactor)`.
pub fn find_d4_curve(p: u64, seed: u64, max_cofactor: u64) -> Result<D4Curve, String> {
    for s in seed..seed + 2_000 {
        let inst = generate_d4_instance(p, s);
        let mut rng = StdRng::seed_from_u64(s ^ 0x0c0f);
        let Some(order) = group_order_d4(&inst.curve, &mut rng) else {
            continue;
        };
        let Some(h) = (2..=max_cofactor).find(|&h| order % h == 0 && is_prime_u64(order / h))
        else {
            continue;
        };
        let r = order / h;
        let mut scratch = GroupOps::default();
        let generator = loop {
            let x = inst.curve.f.random(&mut rng);
            let Some(pt) = inst.curve.lift_x(&x, &mut rng) else {
                continue;
            };
            let g = CountedGroup::mul(&inst.curve, &mut scratch, pt, h);
            if !g.inf {
                break g;
            }
        };
        if !CountedGroup::mul(&inst.curve, &mut scratch, generator, r).inf {
            continue;
        }
        return Ok(D4Curve {
            inst,
            instance_seed: s,
            group_order: order,
            r,
            cofactor: h,
            generator,
        });
    }
    Err(format!(
        "no curve with rational 2-torsion and cofactor ≤ {max_cofactor} over F_{{{p}^4}} from seed {seed}"
    ))
}

/// F_p multiplications one affine addition of two distinct points costs
/// on this curve, by [`Fp4`]'s own count: the ratio the `D₄` oracle's
/// work is priced at.
pub fn fp_muls_per_add(curve: &Curve4) -> u64 {
    let mut rng = StdRng::seed_from_u64(0xadd);
    let pts: Vec<Pt4> = (0..2)
        .map(|_| loop {
            let x = curve.f.random(&mut rng);
            if let Some(pt) = curve.lift_x(&x, &mut rng) {
                break pt;
            }
        })
        .collect();
    let before = curve.f.muls();
    let _ = Curve4::add(curve, &pts[0], &pts[1]);
    curve.f.muls() - before
}

/// **`fghr-y-line`**: `{P : Y(P) ∈ F_p \ {0, ±1}}`, a column per pair
/// `{Y, −Y}`: the four points `±P, ±(P + T)`, with `P + T` computed (one
/// addition, charged) and `[h](P + T) = [h]P` because `h` is even.
pub struct YLineBase<'a> {
    pub curve: &'a D4Curve,
}

impl FactorBaseBuilder<Curve4> for YLineBase<'_> {
    fn name(&self) -> &str {
        "fghr-y-line"
    }
    fn describe(&self, params: &Params) -> String {
        if params.u64_or("zero", 0).unwrap_or(0) == 1 {
            "FGHR's Y-line: Y(P) ∈ F_p, columns {±P, ±(P + T)}, and {±P₀} for Y(P₀) = 0".into()
        } else {
            "FGHR's Y-line: Y(P) ∈ F_p, columns {±P, ±(P + T)}".into()
        }
    }
    fn parameters(&self) -> &[(&str, &str)] {
        &[(
            "zero",
            "1: also the point with Y = 0 (order 4, P₀ + T = −P₀), its own column; default 0, the base without it",
        )]
    }
    fn build(
        &self,
        _ctx: &InstanceCtx<Curve4>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<FactorBase<Pt4>, String> {
        let zero = match params.u64_or("zero", 0)? {
            0 => false,
            1 => true,
            v => return Err(format!("fghr-y-line: zero is 0 or 1, not {v}")),
        };
        let d = self.curve;
        let inst = &d.inst;
        let c = &inst.curve;
        let p = d.p();
        let t = d.torsion();
        let mut rng = StdRng::seed_from_u64(0x9b);
        let (mut points, mut col_of, mut coef_of) = (Vec::new(), Vec::new(), Vec::new());
        let (mut columns, mut lifts) = (0usize, 0u64);
        for y in 2..=(p - 1) / 2 {
            let x = inst.x_of_y(&c.f.from_fp(y));
            lifts += 1;
            let Some(pt) = c.lift_x(&x, &mut rng) else {
                continue;
            };
            // The lift's sign is the sampler's: fix the one with the
            // larger key, so the base does not depend on it.
            let pt = if c.key(&pt) & 1 == 1 {
                pt
            } else {
                Curve4::neg(c, &pt)
            };
            let pt_t = CountedGroup::add(c, ops, pt, t);
            for (q, coef) in [
                (pt, 1),
                (Curve4::neg(c, &pt), d.r - 1),
                (pt_t, 1),
                (Curve4::neg(c, &pt_t), d.r - 1),
            ] {
                points.push(q);
                col_of.push(columns);
                coef_of.push(coef);
            }
            columns += 1;
        }
        // Y = 0: P₀ + T = −P₀ (Y(P + T) = −Y(P)), so 2P₀ = T, P₀ has order
        // 4, 4 | h and [h]P₀ = O.  Its column's logarithm is 0; a relation
        // through it is exact with the summand in the row or out of it.
        if zero {
            let x = inst.x_of_y(&c.f.from_fp(0));
            lifts += 1;
            if let Some(pt) = c.lift_x(&x, &mut rng) {
                let pt = if c.key(&pt) & 1 == 1 {
                    pt
                } else {
                    Curve4::neg(c, &pt)
                };
                let pt_t = CountedGroup::add(c, ops, pt, t);
                if pt_t != Curve4::neg(c, &pt) {
                    return Err("the point with Y = 0 is not of order 4".into());
                }
                for (q, coef) in [(pt, 1), (Curve4::neg(c, &pt), d.r - 1)] {
                    points.push(q);
                    col_of.push(columns);
                    coef_of.push(coef);
                }
                columns += 1;
            }
        }
        if columns == 0 {
            return Err("no Y in F_p carries a point".into());
        }
        let mut fb = FactorBase::from_column_map(
            if zero {
                format!("FGHR Y-line over F_{p}^4 with Y = 0, {columns} columns {{±P, ±(P + T)}}")
            } else {
                format!("FGHR Y-line over F_{p}^4, {columns} columns {{±P, ±(P + T)}}")
            },
            points,
            col_of,
            coef_of,
            columns,
            |q| c.key(q),
            |q| c.key(&Curve4::neg(c, q)),
            |q| pack_e4(p, &q.x),
        )?;
        fb.cost.count("sqrt_solves", lifts);
        Ok(fb)
    }
}

/// **`gaudry-d4`**: the four-summand decomposition by [`solve_d4`] at
/// Macaulay degree `macaulay`.  Its `F_p` multiplications are reported as
/// solver work in the unit `F_p multiplications`, with the curve's own
/// `fp_muls_per_add` in `extra` so the runner prices them by count;
/// the additions that check a decomposition are charged as additions.
pub struct D4Oracle<'a> {
    pub curve: &'a D4Curve,
    pub macaulay: u8,
    pub seed: u64,
    /// Charge the work of a system the matrix did not solve
    /// (`charge_unsolved=1`).  Off by default: the first sessions
    /// (`research/gaudry_quartic_d4_ecbench_20261013`) were recorded
    /// without it, and their replays must keep their counts; their `S`
    /// is understated by those attempts, which the study reports.
    pub charge_unsolved: bool,
    /// `weighted=1`, `chains=1`, `sparsest=1` ([`D4Options`]); all off
    /// by default, the committed sessions' oracle.
    pub opts: D4Options,
    /// `replay=1`: the first solved system's shape replayed on the rest
    /// ([`D4Trace`]); off by default.
    pub replay: bool,
    trace: D4Trace,
    pre: Option<SymmetrisedD4>,
    totals: SolverTotals,
    per_add: u64,
}

impl<'a> D4Oracle<'a> {
    pub fn new(curve: &'a D4Curve, macaulay: u8, seed: u64) -> Self {
        D4Oracle {
            curve,
            macaulay,
            seed,
            charge_unsolved: false,
            opts: D4Options::default(),
            replay: false,
            trace: D4Trace::default(),
            pre: None,
            totals: SolverTotals {
                solver: "d4-macaulay".into(),
                op_unit: "F_p multiplications".into(),
                ..Default::default()
            },
            per_add: fp_muls_per_add(&curve.inst.curve),
        }
    }

    fn shape(&self) -> SystemShape {
        if !self.opts.weighted {
            SystemShape {
                n_vars: 5,
                n_equations: 5,
                degrees: vec![4, 4, 4, 4, 2],
                semi_regular_degree: Some(14),
            }
        } else {
            // Weighted degrees, weights (2, 2, 2, 1).
            SystemShape {
                n_vars: 4,
                n_equations: 4,
                degrees: vec![8, 8, 8, 8],
                semi_regular_degree: Some(26),
            }
        }
    }
}

impl DecompositionOracle<Curve4> for D4Oracle<'_> {
    fn name(&self) -> &str {
        "gaudry-d4"
    }
    fn summands(&self) -> u32 {
        4
    }
    fn describe(&self, _params: &Params) -> String {
        let mut d = if self.opts.weighted {
            format!(
                "S₅ in FGHR's D₄ invariants (σ₁, σ₂, σ₃, π), σ₄ = π², weights (2, 2, 2, 1): 512 solutions, weighted Macaulay degree {}",
                self.macaulay
            )
        } else {
            format!(
                "S₅ in FGHR's D₄ invariants (σ₁, σ₂, σ₃, σ₄, π): Bézout 512, Macaulay degree {}",
                self.macaulay
            )
        };
        if self.opts.chains {
            d.push_str(", eigen step along the chains");
        }
        if self.opts.sparsest {
            d.push_str(", sparsest pivots");
        }
        d
    }
    fn prepare(
        &mut self,
        _ctx: &InstanceCtx<Curve4>,
        _fb: &FactorBase<Pt4>,
        _params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<(), String> {
        let mut rng = StdRng::seed_from_u64(self.seed ^ 0x9e);
        let pre = SymmetrisedD4::precompute(&self.curve.inst, &mut rng, 8)?;
        // The interpolation is solver work done once: priced with it.
        self.totals.ops += pre.precompute_muls;
        self.totals
            .extra
            .insert("precompute_fp_muls".into(), pre.precompute_muls);
        self.totals
            .extra
            .insert("fp_muls_per_add".into(), self.per_add);
        self.pre = Some(pre);
        Ok(())
    }
    fn decompose(
        &mut self,
        _ctx: &InstanceCtx<Curve4>,
        fb: &FactorBase<Pt4>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: Pt4,
    ) -> Option<Vec<usize>> {
        let pre = self.pre.as_ref()?;
        let c = &self.curve.inst.curve;
        let mut rng =
            StdRng::seed_from_u64(self.seed ^ self.totals.calls.wrapping_mul(0x9e37_79b9));
        let adds_before = c.ops();
        let t = std::time::Instant::now();
        let sol = solve_d4_traced(
            &self.curve.inst,
            pre,
            &point,
            self.opts,
            self.macaulay,
            &mut rng,
            false,
            self.replay.then_some(&mut self.trace),
        );
        if self.replay {
            self.totals
                .extra
                .insert("shape_replays".into(), self.trace.replays);
            self.totals
                .extra
                .insert("shape_fallbacks".into(), self.trace.fallbacks);
        }
        self.totals.wall_ns += t.elapsed().as_nanos() as u64;
        // The additions that test a candidate's signs.
        ops.adds += c.ops() - adds_before;
        self.totals.calls += 1;
        self.totals.shape = Some(self.shape());
        let sol = match sol {
            Ok(sol) => sol,
            Err(cost) => {
                *self
                    .totals
                    .extra
                    .entry("degenerate_systems".into())
                    .or_insert(0) += 1;
                if self.charge_unsolved {
                    self.totals.ops += cost.total_muls;
                    *self
                        .totals
                        .extra
                        .entry("degenerate_fp_muls".into())
                        .or_insert(0) += cost.total_muls;
                }
                return None;
            }
        };
        self.totals.ops += sol.solve.total_muls;
        for (k, v) in [
            ("echelon_fp_muls", sol.solve.echelon_muls),
            ("normal_form_fp_muls", sol.solve.nf_muls),
            ("charpoly_fp_muls", sol.solve.charpoly_muls),
            ("rational_solutions", sol.rational_solutions as u64),
        ] {
            *self.totals.extra.entry(k.into()).or_insert(0) += v;
        }
        self.totals.solving_degree_sum += u64::from(self.macaulay);
        self.totals.solving_degree_max =
            self.totals.solving_degree_max.max(u32::from(self.macaulay));
        self.totals.calls_with_a_degree += 1;
        let w = sol.witness?;
        let idx: Option<Vec<usize>> = w.iter().map(|q| fb.index_of_key(c.key(q))).collect();
        if idx.is_none() {
            counters.lift_failures += 1;
        }
        idx
    }
    fn last_system(&self) -> Option<SystemShape> {
        Some(self.shape())
    }
    fn solver_totals(&self) -> Option<SolverTotals> {
        Some(self.totals.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn y_coordinate_negates_under_the_two_torsion_translation() {
        let inst = generate_d4_instance(13, 1);
        let f = &inst.curve.f;
        let mut rng = StdRng::seed_from_u64(3);
        let t = Pt4 {
            x: inst.x0,
            y: E4::ZERO,
            inf: false,
        };
        let mut checked = 0;
        while checked < 20 {
            let x = f.random(&mut rng);
            let Some(pt) = inst.curve.lift_x(&x, &mut rng) else {
                continue;
            };
            let q = inst.curve.add(&pt, &t);
            if q.inf {
                continue;
            }
            assert_eq!(inst.y_of_x(&q.x), f.neg(&inst.y_of_x(&pt.x)));
            assert_eq!(inst.x_of_y(&inst.y_of_x(&pt.x)), pt.x);
            checked += 1;
        }
    }

    #[test]
    fn the_ecbench_curve_and_its_identity_match_the_registry() {
        // `scripts/curve_id.py`'s `extension_id` on the same model gives
        // this slug and ICV1 (the registry's entry for SPEC-k4-p29).
        let c = find_d4_curve(29, 1, 16).unwrap();
        assert_eq!((c.group_order, c.cofactor, c.r), (708_208, 16, 44_263));
        let f = &c.inst.curve.f;
        let id = crate::cryptanalysis::curve_id::extension(
            29,
            &[(29 - f.c % 29) % 29, 0, 0, 0],
            &c.inst.curve.a.0,
            &c.inst.curve.b.0,
            &num_bigint::BigUint::from(c.group_order),
        )
        .unwrap();
        assert_eq!(id.slug, "icv1-fp5k4-tm926-7bc5b359");
        assert_eq!(
            id.icv1,
            "ICV1:fpk-29-4-9dbb234a:-926:708208:5,20,2,25:unk:unk:r:7bc5b3595444"
        );
    }

    #[test]
    fn the_y_zero_column_is_order_four_and_optional() {
        let d = find_d4_curve(29, 1, 16).unwrap();
        let c = &d.inst.curve;
        let ctx = InstanceCtx {
            group: c,
            generator: d.generator,
            target: d.generator,
            r: d.r,
            cofactor: d.cofactor,
            group_order: d.group_order,
            name: "test".into(),
            field_degree: None,
        };
        let base = YLineBase { curve: &d };
        let build = |spec: &str| {
            let (_, params) = Params::parse_spec(spec).unwrap();
            base.build(&ctx, &params, &mut GroupOps::default()).unwrap()
        };
        let plain = build("fghr-y-line");
        let zero = build("fghr-y-line:zero=1");
        // The spec the committed session used keeps its base.
        assert_eq!(plain.points, build("fghr-y-line:zero=0").points);
        assert_eq!(zero.columns, plain.columns + 1);
        assert_eq!(&zero.points[..plain.points.len()], &plain.points[..]);
        let p0 = zero.points[plain.points.len()];
        assert!(inst_y_is_zero(&d, &p0));
        let mut ops = GroupOps::default();
        let four = CountedGroup::mul(c, &mut ops, p0, 4);
        assert!(four.inf && !CountedGroup::mul(c, &mut ops, p0, 2).inf);
        assert!(CountedGroup::mul(c, &mut ops, p0, d.cofactor).inf);
    }

    fn inst_y_is_zero(d: &D4Curve, q: &Pt4) -> bool {
        d.inst.y_of_x(&q.x).is_zero()
    }

    /// A target that decomposes: four base points with signs, summed.
    fn planted(inst: &D4Instance, base: &[(u64, Pt4)], rng: &mut StdRng) -> Pt4 {
        loop {
            let mut idx: Vec<usize> = Vec::new();
            while idx.len() < 4 {
                let i = rng.gen_range(0..base.len());
                if !idx.contains(&i) {
                    idx.push(i);
                }
            }
            let mut acc = Pt4::INF;
            for &i in &idx {
                let q = if rng.gen_bool(0.5) {
                    inst.curve.neg(&base[i].1)
                } else {
                    base[i].1
                };
                acc = inst.curve.add(&acc, &q);
            }
            if !acc.inf {
                return acc;
            }
        }
    }

    #[test]
    fn every_option_finds_the_meet_in_the_middle_decompositions() {
        // The ecbench curve (p = 29): F_p small, so `det A` is scanned on
        // all of F_p, and π-collisions between solutions are likely.
        let d = find_d4_curve(29, 1, 16).unwrap();
        let inst = &d.inst;
        let mut rng = StdRng::seed_from_u64(7);
        let pre = SymmetrisedD4::precompute(inst, &mut rng, 8).unwrap();
        let base = y_base(inst, &mut rng);
        let all = D4Options {
            weighted: true,
            chains: true,
            sparsest: true,
        };
        let mut trace = D4Trace::default();
        let mut agree = 0;
        for k in 0..3 {
            let r = planted(inst, &base, &mut rng);
            let mitm = mitm_y_decompositions(&inst.curve, &base, &r);
            let mut s1 = StdRng::seed_from_u64(k);
            let plain = solve_d4_form(inst, &pre, &r, D4Options::default(), 14, &mut s1, false);
            let mut s2 = StdRng::seed_from_u64(k);
            let new = solve_d4_traced(inst, &pre, &r, all, 26, &mut s2, false, Some(&mut trace));
            if let (Ok(a), Ok(b)) = (plain, new) {
                let a: BTreeSet<[u64; 4]> = a.decompositions.into_iter().collect();
                let b: BTreeSet<[u64; 4]> = b.decompositions.into_iter().collect();
                assert_eq!(a, mitm);
                assert_eq!(b, mitm);
                assert!(!b.is_empty());
                agree += 1;
            }
        }
        assert!(agree >= 2, "{agree} of 3 planted targets solved by both");
        assert!(trace.replays + trace.fallbacks >= 1);
    }

    #[test]
    fn the_chain_step_matches_the_dense_eigen_step_beyond_the_scan() {
        // p = 613 > dim + 1: `det A` is interpolated and its roots found.
        let inst = generate_d4_instance(613, 1);
        let mut rng = StdRng::seed_from_u64(11);
        let pre = SymmetrisedD4::precompute(&inst, &mut rng, 8).unwrap();
        let base = y_base(&inst, &mut rng);
        for k in 0..2 {
            let r = planted(&inst, &base, &mut rng);
            let run = |chains: bool| {
                let opts = D4Options {
                    weighted: true,
                    chains,
                    sparsest: true,
                };
                let mut s = StdRng::seed_from_u64(k);
                solve_d4_form(&inst, &pre, &r, opts, 26, &mut s, false).unwrap()
            };
            let (dense, chain) = (run(false), run(true));
            assert_eq!(dense.solve.dim, 512);
            assert_eq!(dense.rational_solutions, chain.rational_solutions);
            assert_eq!(dense.decompositions, chain.decompositions);
            assert_eq!(
                dense.solve.rational_eigenvalues,
                chain.solve.rational_eigenvalues
            );
            assert!(chain.solve.total_muls * 3 < dense.solve.total_muls);
            let mitm = mitm_y_decompositions(&inst.curve, &base, &r);
            assert_eq!(
                chain.decompositions.into_iter().collect::<BTreeSet<_>>(),
                mitm
            );
        }
    }

    #[test]
    fn s5_in_y_lies_in_the_span_of_the_d4_invariants() {
        assert_eq!(d4_monomials().len(), 105);
        let inst = generate_d4_instance(13, 2);
        let mut rng = StdRng::seed_from_u64(5);
        let pre = SymmetrisedD4::precompute(&inst, &mut rng, 8).expect("D₄-invariant");
        assert_eq!(pre.checked_points, 8);
    }
}
