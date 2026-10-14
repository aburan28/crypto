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
    solve_system_n_inner, split_quartic, Curve4, Fp4, Pt4, QuarticSolve, E4,
};
use super::residual_walk::{pow_mod, sqrt_mod};

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

fn signed_sum_is(curve: &Curve4, pts: &[Pt4; 4], r: &Pt4) -> bool {
    for s in 0..16u32 {
        let mut acc = Pt4 {
            x: E4::ZERO,
            y: E4::ZERO,
            inf: true,
        };
        for (i, q) in pts.iter().enumerate() {
            let q = if s >> i & 1 == 1 { curve.neg(q) } else { *q };
            acc = curve.add(&acc, &q);
        }
        if !acc.inf && acc.x == r.x && acc.y == r.y {
            return true;
        }
    }
    false
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
    let f = &inst.curve.f;
    let p = f.p;
    let mut out = D4Solve::default();
    let before = f.muls();
    let y_r = inst.y_of_x(&r.x);
    let comps = pre.system(f, &y_r);
    out.solve.weil_muls = f.muls() - before;
    let sols = solve_system_n_inner::<5>(
        &comps,
        0,
        macaulay,
        p,
        rng,
        &mut out.solve,
        debug,
        None,
        None,
    )
    .ok()?;
    let Some(sols) = sols else {
        if debug {
            eprintln!(
                "  unsolved: {} (dim {}, pivots {}, zero rows {})",
                out.solve.outcome, out.solve.dim, out.solve.pivots, out.solve.zero_rows
            );
        }
        return None;
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
            if signed_sum_is(&inst.curve, &pts, r) {
                found.insert(sorted);
            }
        }
    }
    out.decompositions = found.into_iter().collect();
    out.solve.close_public();
    Some(out)
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
    fn s5_in_y_lies_in_the_span_of_the_d4_invariants() {
        assert_eq!(d4_monomials().len(), 105);
        let inst = generate_d4_instance(13, 2);
        let mut rng = StdRng::seed_from_u64(5);
        let pre = SymmetrisedD4::precompute(&inst, &mut rng, 8).expect("D₄-invariant");
        assert_eq!(pre.checked_points, 8);
    }
}
