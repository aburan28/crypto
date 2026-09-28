//! **Joux–Vitse `(k − 1)`-point decompositions at `k = 4`: the one route
//! whose ratio to rho is a constant, measured.**
//!
//! Companion to `research/notes/index-calculus/RESEARCH_RHO_PARITY_PROGRAMME.md`.
//! On `E(F_{p⁴})` with the subspace base `F = {P : x(P) ∈ F_p}`, a residual
//! is a sum of *three* base points about `1/(6p)` of the time, so a
//! relation phase built on three-point decompositions tries `≈ 3p²`
//! residuals — `∝ n^{1/2}`, rho's own exponent — and the linear algebra over
//! `|F| ≈ p/2` unknowns is `∝ p² = n^{1/2}` too.  The method's `S / rho`
//! therefore tends to a constant, `C′ / (S_rho · c_add) · 6 + r∞` with `C′`
//! the cost of one three-point test, and that constant is the distance to
//! parity at every size.  `RESEARCH_RESIDUAL_WALKS.md` §11.16 derived it as
//! `≈ 36×` on the assumption `C′ ≈ 1,513` (the two-unknown pair test of
//! `k = 3`); this module measures `C′` and the whole method.
//!
//! The test: `S₄(x₁, x₂, x₃, x_R) = 0` symmetrised in `(x₁, x₂, x₃)` is a
//! polynomial `H(e₁, e₂, e₃; x_R)` of total degree `≤ 4` in the elementary
//! symmetric functions, found by interpolation from the resultant
//! ([`SymmetrisedS4Q`]).  Its four `F_p`-components are four equations of
//! degree `≤ 4` in three unknowns — overdetermined — solved by the
//! repository's F4 ([`super::f4_fp::solve`]); an `e`-solution whose cubic
//! `T³ − e₁T² + e₂T − e₃` splits over the base's abscissae is a
//! decomposition, and the signs are read off with three group operations.
//! Every test in the cost experiment, and a sample in the end-to-end run,
//! is checked against an independent meet-in-the-middle oracle over the
//! pair table ([`mitm3_signed`]); every relation is verified in the group
//! before it is stored.

use std::collections::{BTreeSet, HashMap};
use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::Serialize;

use super::f4_fp::{self, F4Options, Ordering as F4Ordering, Verdict};
use super::gaudry_cubic::{square_core, wiedemann_u64, PolyRing, SparseRel, UPoly};
use super::gaudry_quartic::{
    factor_base, generate_instance4, pair_table, rho4, Curve4, Fp4, PairTable, Pt4, RhoRun4, E4,
};

// ── The symmetrised S₄ over F_{p⁴}, by interpolation ──────────────────────

/// Monomials `e₁^a e₂^b e₃^c` of total degree `≤ d`.
fn exponents3(d: u8) -> Vec<[u8; 3]> {
    let mut v = Vec::new();
    for a in 0..=d {
        for b in 0..=(d - a) {
            for c in 0..=(d - a - b) {
                v.push([a, b, c]);
            }
        }
    }
    v
}

fn elementary3(f: &Fp4, x: &[E4; 3]) -> [E4; 3] {
    let e1 = f.add(&f.add(&x[0], &x[1]), &x[2]);
    let e2 = f.add(
        &f.add(&f.mul(&x[0], &x[1]), &f.mul(&x[0], &x[2])),
        &f.mul(&x[1], &x[2]),
    );
    let e3 = f.mul(&f.mul(&x[0], &x[1]), &x[2]);
    [e1, e2, e3]
}

fn mono_values3(f: &Fp4, e: &[E4; 3], monos: &[[u8; 3]]) -> Vec<E4> {
    let mut pows = [[E4::ONE; 5]; 3];
    for i in 0..3 {
        for k in 1..5 {
            pows[i][k] = f.mul(&pows[i][k - 1], &e[i]);
        }
    }
    monos
        .iter()
        .map(|m| {
            f.mul(
                &f.mul(&pows[0][m[0] as usize], &pows[1][m[1] as usize]),
                &pows[2][m[2] as usize],
            )
        })
        .collect()
}

/// `S₄(x₁, x₂, x₃, X) = H(e₁, e₂, e₃, X)`, `H` of total degree `≤ 4` in the
/// `e`s and degree `≤ 4` in `X`: `35 · 5` coefficients in `F_{p⁴}`, found by
/// interpolation of the resultant `Res_Z(S₃(x₁, x₂, Z), S₃(x₃, X, Z))` at
/// `35` random points and checked against fresh evaluations.
pub struct SymmetrisedS4Q {
    pub monos: Vec<[u8; 3]>,
    pub coef: Vec<[E4; 5]>,
    /// `F_p` multiplications spent building it (once per curve).
    pub precompute_muls: u64,
}

impl SymmetrisedS4Q {
    fn coefficients_in_last(curve: &Curve4, x: &[E4; 3]) -> [E4; 5] {
        let mut v = curve.s4_in_last(&x[0], &x[1], &x[2]);
        v.resize(5, E4::ZERO);
        core::array::from_fn(|j| v[j])
    }

    pub fn precompute(curve: &Curve4, rng: &mut StdRng) -> SymmetrisedS4Q {
        let f = &curve.f;
        let before = f.muls();
        let monos = exponents3(4);
        let n = monos.len();
        loop {
            let mut sys: Vec<Vec<E4>> = Vec::with_capacity(n);
            for _ in 0..n {
                let x: [E4; 3] = core::array::from_fn(|_| f.random(rng));
                let e = elementary3(f, &x);
                let mut row = mono_values3(f, &e, &monos);
                row.extend(Self::coefficients_in_last(curve, &x));
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
            let coef: Vec<[E4; 5]> = (0..n)
                .map(|r| core::array::from_fn(|j| sys[r][n + j]))
                .collect();
            let pre = SymmetrisedS4Q {
                monos: monos.clone(),
                coef,
                precompute_muls: 0,
            };
            let good = (0..8).all(|_| {
                let x: [E4; 3] = core::array::from_fn(|_| f.random(rng));
                let want = Self::coefficients_in_last(curve, &x);
                let e = elementary3(f, &x);
                let mv = mono_values3(f, &e, &pre.monos);
                (0..5).all(|j| {
                    let got = mv
                        .iter()
                        .zip(&pre.coef)
                        .fold(E4::ZERO, |acc, (m, c)| f.add(&acc, &f.mul(m, &c[j])));
                    got == want[j]
                })
            });
            assert!(good, "the symmetrised S₄ does not reproduce the resultant");
            return SymmetrisedS4Q {
                precompute_muls: f.muls() - before,
                ..pre
            };
        }
    }

    /// `H(e; x_R)` at an `F_{p⁴}` point (tests and cross-checks).
    pub fn eval(&self, f: &Fp4, e: &[E4; 3], x_r: &E4) -> E4 {
        let mut pows = [E4::ONE; 5];
        for k in 1..5 {
            pows[k] = f.mul(&pows[k - 1], x_r);
        }
        let mv = mono_values3(f, e, &self.monos);
        mv.iter().zip(&self.coef).fold(E4::ZERO, |acc, (m, c)| {
            let v = (0..5).fold(E4::ZERO, |a, j| f.add(&a, &f.mul(&c[j], &pows[j])));
            f.add(&acc, &f.mul(m, &v))
        })
    }

    /// The four `F_p`-components of `H(e₁, e₂, e₃, x_R)`: four polynomials
    /// of total degree `≤ 4` in three `F_p` unknowns.
    pub fn weil_restrict(&self, f: &Fp4, x_r: &E4) -> [HashMap<[u8; 3], u64>; 4] {
        let mut pows = [E4::ONE; 5];
        for k in 1..5 {
            pows[k] = f.mul(&pows[k - 1], x_r);
        }
        let mut out: [HashMap<[u8; 3], u64>; 4] = Default::default();
        for (m, c) in self.monos.iter().zip(&self.coef) {
            let v = (0..5).fold(E4::ZERO, |acc, j| f.add(&acc, &f.mul(&c[j], &pows[j])));
            for (i, comp) in out.iter_mut().enumerate() {
                if v.0[i] != 0 {
                    comp.insert(*m, v.0[i]);
                }
            }
        }
        out
    }
}

// ── The three-point test ──────────────────────────────────────────────────

/// One decomposition `R = Σ s_t P_{i_t}` over three distinct base points.
pub type Triple = [(usize, i64); 3];

/// What one three-point test cost and did.
#[derive(Clone, Debug, Default, Serialize)]
pub struct Jv3Cost {
    /// `F_p` multiplications of the Weil restriction.
    pub weil_muls: u64,
    /// `F_p` multiplications of F4's row reductions (this report's own
    /// difference of the process counter: exact when tests run one at a
    /// time, interleaved when they run concurrently — the batch runners
    /// charge the batch difference instead).
    pub f4_muls: u64,
    /// `F_p` multiplications of the cubic splits.
    pub roots_muls: u64,
    /// Group operations of the sign resolution.
    pub group_ops: u64,
    pub solving_degree: u32,
    pub max_rows: usize,
    pub max_cols: usize,
    pub f4_ms: f64,
    pub inconsistent: bool,
    pub undetermined: bool,
    /// `e`-solutions F4 returned (before the split).
    pub e_solutions: usize,
    pub decompositions: usize,
}

fn polys_of(comps: &[HashMap<[u8; 3], u64>; 4], p: u64) -> Vec<f4_fp::Poly> {
    comps
        .iter()
        .map(|c| {
            let terms: Vec<(Vec<u32>, u64)> = c
                .iter()
                .map(|(e, v)| (e.iter().map(|&x| x as u32).collect(), *v))
                .collect();
            f4_fp::normalise(&terms, p, F4Ordering::Grevlex)
        })
        .filter(|f| !f.is_empty())
        .collect()
}

/// The F4 options every test uses: grevlex, pairs bounded at degree `12`
/// (the certificate of an inconsistent generic quartic system in three
/// unknowns lives at degree `7`), no deadline.
pub fn jv3_options() -> F4Options {
    F4Options::new(F4Ordering::Grevlex, 12)
}

/// The signs of a candidate triple by group arithmetic: `R − s₂P₂ − s₃P₃ =
/// ±P₁` for one of four `(s₂, s₃)`, `s₁` read from the `y`-coordinate.
fn resolve_signs(curve: &Curve4, base: &[Pt4], r: &Pt4, idx: [usize; 3]) -> Option<Triple> {
    let [i1, i2, i3] = idx;
    for s2 in [1i64, -1] {
        for s3 in [1i64, -1] {
            let p2 = if s2 == 1 {
                base[i2]
            } else {
                curve.neg(&base[i2])
            };
            let p3 = if s3 == 1 {
                base[i3]
            } else {
                curve.neg(&base[i3])
            };
            let t = curve.sub(&curve.sub(r, &p2), &p3);
            if t.inf || t.x != base[i1].x {
                continue;
            }
            let s1 = if t.y == base[i1].y { 1 } else { -1 };
            let mut terms = [(i1, s1), (i2, s2), (i3, s3)];
            terms.sort_unstable();
            return Some(terms);
        }
    }
    None
}

/// Every three-point decomposition of `r` over the base, by the symmetrised
/// `S₄` and F4, with the cost split.  `curve.f` and `curve` carry the
/// multiplication and group-operation counters the cost is read from.
pub fn jv3_decompose(
    curve: &Curve4,
    base: &[Pt4],
    by_x: &HashMap<u64, usize>,
    pre: &SymmetrisedS4Q,
    r: &Pt4,
    opts: &F4Options,
    rng: &mut StdRng,
) -> (Vec<Triple>, Jv3Cost) {
    let f = &curve.f;
    let p = f.p;
    let mut cost = Jv3Cost::default();
    let mut out: BTreeSet<Triple> = BTreeSet::new();
    if r.inf {
        return (Vec::new(), cost);
    }
    let m0 = f.muls();
    let comps = pre.weil_restrict(f, &r.x);
    cost.weil_muls = f.muls() - m0;
    let polys = polys_of(&comps, p);
    let rep = f4_fp::solve(&polys, 3, p, opts);
    cost.f4_muls = rep.field_ops;
    cost.solving_degree = rep.solving_degree;
    cost.max_rows = rep.max_rows;
    cost.max_cols = rep.max_cols;
    cost.f4_ms = rep.ms;
    let es = match rep.verdict {
        Verdict::Inconsistent => {
            cost.inconsistent = true;
            Vec::new()
        }
        Verdict::Undetermined => {
            cost.undetermined = true;
            Vec::new()
        }
        Verdict::Solutions(s) => s,
    };
    cost.e_solutions = es.len();
    let ring = PolyRing::new(p);
    let g0 = curve.ops();
    for e in &es {
        // T³ − e₁T² + e₂T − e₃, low to high.
        let cubic = UPoly(vec![(p - e[2]) % p, e[1], (p - e[0]) % p, 1]);
        let roots = ring.roots(&cubic, rng);
        if roots.len() != 3 {
            continue;
        }
        let Some(idx) = roots
            .iter()
            .map(|x| by_x.get(x).copied())
            .collect::<Option<Vec<usize>>>()
        else {
            continue;
        };
        if idx[0] == idx[1] || idx[1] == idx[2] || idx[0] == idx[2] {
            continue;
        }
        if let Some(t) = resolve_signs(curve, base, r, [idx[0], idx[1], idx[2]]) {
            out.insert(t);
        }
    }
    cost.roots_muls = ring.muls.get();
    cost.group_ops = curve.ops() - g0;
    cost.decompositions = out.len();
    (out.into_iter().collect(), cost)
}

/// The independent oracle: every `R = s_k P_k + t (P_i + s P_j)` by meet
/// in the middle against the pair table, `2|F|` group operations.
pub fn mitm3_signed(curve: &Curve4, base: &[Pt4], table: &PairTable, r: &Pt4) -> Vec<Triple> {
    let mut out: BTreeSet<Triple> = BTreeSet::new();
    for k in 0..base.len() {
        for sk in [1i64, -1] {
            let pk = if sk == 1 {
                base[k]
            } else {
                curve.neg(&base[k])
            };
            let v = curve.sub(r, &pk);
            if v.inf {
                continue;
            }
            let Some(cands) = table.map.get(&v.x) else {
                continue;
            };
            for &(i, j, s, y) in cands {
                let (i, j) = (i as usize, j as usize);
                if i == k || j == k {
                    continue;
                }
                let t: i64 = if v.y == y { 1 } else { -1 };
                let mut terms = [(k, sk), (i, t), (j, t * s as i64)];
                terms.sort_unstable();
                out.insert(terms);
            }
        }
    }
    out.into_iter().collect()
}

fn verify_triple(curve: &Curve4, base: &[Pt4], r: &Pt4, t: &Triple) -> bool {
    let mut acc = *r;
    for &(i, s) in t {
        acc = if s == 1 {
            curve.sub(&acc, &base[i])
        } else {
            curve.add(&acc, &base[i])
        };
    }
    acc.inf
}

// ── Experiment A: C′, the cost of one three-point test ────────────────────

#[derive(Clone, Debug, Default, Serialize)]
pub struct CostStats {
    pub mean: f64,
    pub min: u64,
    pub max: u64,
}

fn stats(v: &[u64]) -> CostStats {
    if v.is_empty() {
        return CostStats::default();
    }
    CostStats {
        mean: v.iter().sum::<u64>() as f64 / v.len() as f64,
        min: *v.iter().min().unwrap(),
        max: *v.iter().max().unwrap(),
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct CprimeReport {
    pub p: u64,
    pub seed: u64,
    pub n: u64,
    pub bits: f64,
    pub base: usize,
    pub fp_muls_per_add: f64,
    pub precompute_muls: u64,
    pub random_residuals: usize,
    pub constructed_residuals: usize,
    /// Random residuals with at least one three-point decomposition, and
    /// the rate against `1/(6p)`.
    pub random_decomposable: usize,
    pub decomposition_rate: f64,
    pub expected_rate: f64,
    /// Constructed residuals whose planted triple the test returned.
    pub planted_found: usize,
    /// Residuals (random and constructed) where the F4 test and the
    /// meet-in-the-middle oracle disagree.
    pub mismatches: usize,
    pub undetermined: usize,
    /// Per-test cost over the random residuals, `F_p` multiplications:
    /// total `C′` and its split.
    pub c_prime: CostStats,
    pub weil_muls: CostStats,
    pub f4_muls: CostStats,
    pub roots_muls: CostStats,
    pub sign_group_ops: CostStats,
    pub f4_ms: CostStats,
    pub solving_degree: CostStats,
    pub max_rows: CostStats,
    pub max_cols: CostStats,
    /// The same over the constructed (decomposable) residuals.
    pub c_prime_decomposable: CostStats,
    pub solving_degree_decomposable: CostStats,
    /// The oracle's cost per residual, group operations, for scale.
    pub mitm3_group_ops: f64,
    pub wall_ms: f64,
}

/// `C′` measured: `random` residuals `aG + bQ` and `constructed` sums of
/// three base points, every one cross-checked against [`mitm3_signed`].
pub fn run_jv4_cprime(p: u64, seed: u64, random: usize, constructed: usize) -> CprimeReport {
    let start = Instant::now();
    let inst = generate_instance4(p, seed);
    let curve = &inst.curve;
    let f = &curve.f;
    let n = inst.n;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x3C01);
    let base = factor_base(curve, &mut rng);
    let by_x: HashMap<u64, usize> = base
        .iter()
        .enumerate()
        .map(|(i, q)| (q.x.0[0], i))
        .collect();
    let fp_per_add = {
        f.reset_muls();
        let mut acc = base[0];
        for _ in 0..64 {
            acc = curve.add(&acc, &base[1]);
        }
        f.muls() as f64 / 64.0
    };
    let pre = SymmetrisedS4Q::precompute(curve, &mut rng);
    let table = pair_table(curve, &base);
    let opts = jv3_options();
    let mut rep = CprimeReport {
        p,
        seed,
        n,
        bits: (n as f64).log2(),
        base: base.len(),
        fp_muls_per_add: fp_per_add,
        precompute_muls: pre.precompute_muls,
        random_residuals: random,
        constructed_residuals: constructed,
        expected_rate: 1.0 / (6.0 * p as f64),
        ..Default::default()
    };
    let mut costs: Vec<Jv3Cost> = Vec::new();
    let mut costs_dec: Vec<Jv3Cost> = Vec::new();
    let mut mitm_ops = 0u64;
    let mut check_rng = StdRng::seed_from_u64(seed ^ 0x5EED);
    let mut check = |r: &Pt4, want_planted: Option<Triple>, rep: &mut CprimeReport| -> Jv3Cost {
        let (found, cost) = jv3_decompose(curve, &base, &by_x, &pre, r, &opts, &mut check_rng);
        let g0 = curve.ops();
        let oracle = mitm3_signed(curve, &base, &table, r);
        mitm_ops += curve.ops() - g0;
        if found != oracle {
            rep.mismatches += 1;
        }
        if cost.undetermined {
            rep.undetermined += 1;
        }
        if let Some(t) = want_planted {
            if found.contains(&t) {
                rep.planted_found += 1;
            }
        }
        cost
    };
    for _ in 0..random {
        let a = rng.gen_range(0..n);
        let b = rng.gen_range(1..n);
        let r = curve.add(&curve.mul(&inst.g, a), &curve.mul(&inst.q, b));
        let c = check(&r, None, &mut rep);
        if c.decompositions > 0 {
            rep.random_decomposable += 1;
        }
        costs.push(c);
    }
    for _ in 0..constructed {
        let idx: [usize; 3] = loop {
            let idx: [usize; 3] = core::array::from_fn(|_| rng.gen_range(0..base.len()));
            if idx[0] != idx[1] && idx[1] != idx[2] && idx[0] != idx[2] {
                break idx;
            }
        };
        let signs: [i64; 3] = core::array::from_fn(|_| if rng.gen_bool(0.5) { 1 } else { -1 });
        let mut r = Pt4::INF;
        for t in 0..3 {
            let q = if signs[t] == 1 {
                base[idx[t]]
            } else {
                curve.neg(&base[idx[t]])
            };
            r = curve.add(&r, &q);
        }
        let mut terms = [(idx[0], signs[0]), (idx[1], signs[1]), (idx[2], signs[2])];
        terms.sort_unstable();
        let c = check(&r, Some(terms), &mut rep);
        costs_dec.push(c);
    }
    let total = |c: &Jv3Cost| {
        c.weil_muls + c.f4_muls + c.roots_muls + (c.group_ops as f64 * fp_per_add) as u64
    };
    let col = |v: &[Jv3Cost], g: &dyn Fn(&Jv3Cost) -> u64| -> CostStats {
        stats(&v.iter().map(g).collect::<Vec<u64>>())
    };
    rep.decomposition_rate = rep.random_decomposable as f64 / random.max(1) as f64;
    rep.c_prime = col(&costs, &total);
    rep.weil_muls = col(&costs, &|c| c.weil_muls);
    rep.f4_muls = col(&costs, &|c| c.f4_muls);
    rep.roots_muls = col(&costs, &|c| c.roots_muls);
    rep.sign_group_ops = col(&costs, &|c| c.group_ops);
    rep.f4_ms = col(&costs, &|c| (c.f4_ms * 1000.0) as u64);
    rep.f4_ms.mean /= 1000.0;
    rep.solving_degree = col(&costs, &|c| c.solving_degree as u64);
    rep.max_rows = col(&costs, &|c| c.max_rows as u64);
    rep.max_cols = col(&costs, &|c| c.max_cols as u64);
    rep.c_prime_decomposable = col(&costs_dec, &total);
    rep.solving_degree_decomposable = col(&costs_dec, &|c| c.solving_degree as u64);
    rep.mitm3_group_ops = mitm_ops as f64 / (random + constructed).max(1) as f64;
    rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
    rep
}

// ── Experiment B: the method end to end ───────────────────────────────────

#[derive(Clone, Debug, Default, Serialize)]
pub struct Jv4DlpReport {
    pub p: u64,
    pub seed: u64,
    pub n: u64,
    pub bits: f64,
    pub base: usize,
    pub fp_muls_per_add: f64,
    /// Residuals the walk produced and tested.
    pub residuals: u64,
    pub decompositions: u64,
    pub decomposition_rate: f64,
    pub expected_rate: f64,
    /// Residuals with a decomposition, against the floor `|F| + 1` a square
    /// system needs.
    pub relations: usize,
    pub floor_residuals: f64,
    pub residual_ratio: f64,
    pub undetermined: u64,
    /// Residuals also run through the oracle, and disagreements.
    pub cross_checked: u64,
    pub mismatches: u64,
    /// Every phase, in `F_p` multiplications.
    pub setup_muls: u64,
    pub precompute_muls: u64,
    pub walk_group_ops: u64,
    pub walk_muls: u64,
    pub weil_muls: u64,
    pub f4_muls: u64,
    pub roots_muls: u64,
    pub sign_group_ops: u64,
    pub sign_muls: u64,
    pub verify_group_ops: u64,
    pub verify_muls: u64,
    pub oracle_muls: u64,
    /// `C′` as the run paid it: oracle multiplications per residual.
    pub c_prime: f64,
    pub unknowns: usize,
    pub filtered_out: usize,
    pub phi: f64,
    pub row_weight: f64,
    pub la_attempts: u64,
    pub la_ops: u64,
    pub la_muls: u64,
    pub wiedemann_constant: f64,
    pub total_muls: u64,
    pub s: f64,
    pub solved: bool,
    pub correct: bool,
    pub rho: Vec<RhoRun4>,
    pub rho_s_mean: f64,
    pub rho_s_sd: f64,
    /// `S / rho`, and its two parts: the relation phase (walk and oracle)
    /// and the linear algebra (`r`, against §11.19's `0.518`).
    pub s_over_rho: f64,
    pub relation_over_rho: f64,
    pub r: f64,
    /// The constant §11.16's formula predicts from this run's own `C′`:
    /// `6·C′ / (S_rho · c_add) + r`.
    pub predicted_s_over_rho: f64,
    pub wall_ms: f64,
}

/// The method end to end: the residual stream `R_0 + i·M` over `⟨G, Q⟩`
/// (one group operation per residual), the three-point test on every residual,
/// verified relations, filtering and Wiedemann as `run_k4_la` does, the
/// logarithm checked, and `rho_runs` rho walks on the same group.  Every
/// `check_every`-th residual is also run through the oracle.
pub fn run_jv4_dlp(p: u64, seed: u64, rho_runs: usize, check_every: u64) -> Jv4DlpReport {
    let start = Instant::now();
    let inst = generate_instance4(p, seed);
    let curve = &inst.curve;
    let f = &curve.f;
    let n = inst.n;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x3D1F);
    let base = factor_base(curve, &mut rng);
    let by_x: HashMap<u64, usize> = base
        .iter()
        .enumerate()
        .map(|(i, q)| (q.x.0[0], i))
        .collect();
    let fp_per_add = {
        f.reset_muls();
        let mut acc = base[0];
        for _ in 0..64 {
            acc = curve.add(&acc, &base[1]);
        }
        f.muls() as f64 / 64.0
    };
    let table = pair_table(curve, &base);
    curve.reset_ops();
    f.reset_muls();
    let pre = SymmetrisedS4Q::precompute(curve, &mut rng);
    let spec = inst.spec();
    let small = base.len();
    let unknowns = small + 1;
    let opts = jv3_options();
    let mut rep = Jv4DlpReport {
        p,
        seed,
        n,
        bits: (n as f64).log2(),
        base: small,
        fp_muls_per_add: fp_per_add,
        precompute_muls: pre.precompute_muls,
        expected_rate: 1.0 / (6.0 * p as f64),
        floor_residuals: unknowns as f64 * 6.0 * p as f64,
        ..Default::default()
    };
    // The stream's start and step, charged as setup.  Residuals are the
    // arithmetic progression `R_i = R_0 + i·M`, one group operation each:
    // unlike an r-adding walk it cannot revisit a point before `n` steps
    // (a walk on a group of order `2^32` cycles after `≈ √n ≈ 10⁵` steps,
    // fewer than the `3p²` residuals the method needs, and a repeated
    // residual with a second `(a, b)` is rho's collision, not a relation).
    f.reset_muls();
    let g0 = curve.ops();
    let (al, be) = (rng.gen_range(1..n), rng.gen_range(1..n));
    let step = curve.add(&curve.mul(&inst.g, al), &curve.mul(&inst.q, be));
    let mut a = rng.gen_range(0..n);
    let mut b = rng.gen_range(1..n);
    let mut r = curve.add(&curve.mul(&inst.g, a), &curve.mul(&inst.q, b));
    let setup_ops = curve.ops() - g0;
    rep.setup_muls = (setup_ops as f64 * fp_per_add) as u64 + f.muls();
    f.reset_muls();
    let walk_start = curve.ops();
    let mut full_rels: Vec<SparseRel> = Vec::new();
    let mut next_attempt = unknowns;
    let mut la_ops = 0u64;
    let cap = 200 * p * p;
    'collect: while rep.residuals < cap {
        // A batch of residuals from the stream, one operation each.
        let mut batch: Vec<(u64, u64, Pt4, u64)> = Vec::with_capacity(64);
        for _ in 0..64 {
            let idx = rep.residuals + batch.len() as u64;
            batch.push((a, b, r, idx));
            r = curve.add(&r, &step);
            a = (a + al) % n;
            b = (b + be) % n;
        }
        let f4_before = f4_fp::field_ops_total();
        let found: Vec<(
            u64,
            u64,
            Pt4,
            Vec<Triple>,
            Jv3Cost,
            u64,
            Option<Vec<Triple>>,
        )> = batch
            .par_iter()
            .map(|&(a, b, r, idx)| {
                let c = Curve4::with_params(spec.p, spec.a, spec.b);
                let mut lrng =
                    StdRng::seed_from_u64(seed ^ idx.wrapping_mul(0x9E37_79B9_7F4A_7C15));
                let (decs, cost) = jv3_decompose(&c, &base, &by_x, &pre, &r, &opts, &mut lrng);
                let oracle = (check_every > 0 && idx % check_every == 0)
                    .then(|| mitm3_signed(&c, &base, &table, &r));
                (a, b, r, decs, cost, c.f.muls(), oracle)
            })
            .collect();
        rep.f4_muls += f4_fp::field_ops_total() - f4_before;
        for (a, b, r, decs, cost, _muls, oracle) in found {
            rep.residuals += 1;
            rep.weil_muls += cost.weil_muls;
            rep.roots_muls += cost.roots_muls;
            rep.sign_group_ops += cost.group_ops;
            if cost.undetermined {
                rep.undetermined += 1;
            }
            if let Some(o) = oracle {
                rep.cross_checked += 1;
                if o != decs {
                    rep.mismatches += 1;
                }
            }
            let mut any = false;
            for t in decs {
                let v0 = curve.ops();
                let ok = verify_triple(curve, &base, &r, &t);
                rep.verify_group_ops += curve.ops() - v0;
                if !ok {
                    continue;
                }
                any = true;
                rep.decompositions += 1;
                let mut cols: Vec<(usize, u64)> = t
                    .iter()
                    .map(|&(i, s)| (i, if s == 1 { 1 } else { n - 1 }))
                    .collect();
                cols.push((small, (n - b) % n));
                cols.sort_unstable();
                full_rels.push(SparseRel { cols, rhs: a });
            }
            if any {
                rep.relations += 1;
            }
            if full_rels.len() >= next_attempt {
                if let Some(core) = square_core(&full_rels, small) {
                    rep.la_attempts += 1;
                    let mut map = vec![usize::MAX; unknowns];
                    let mut k = 0;
                    for c in 0..unknowns {
                        if core.columns.get(c).copied().unwrap_or(false) {
                            map[c] = k;
                            k += 1;
                        }
                    }
                    let sel: Vec<SparseRel> = core
                        .rows
                        .iter()
                        .map(|&i| SparseRel {
                            cols: full_rels[i]
                                .cols
                                .iter()
                                .map(|&(c, v)| (map[c], v))
                                .collect(),
                            rhs: full_rels[i].rhs,
                        })
                        .collect();
                    let x = wiedemann_u64(&sel, sel.len(), n, &mut rng, &mut la_ops);
                    let dd = x.map(|x| x[map[small]]);
                    if dd.is_some_and(|dd| curve.mul(&inst.g, dd) == inst.q) {
                        rep.solved = true;
                        rep.correct = dd == Some(inst.d);
                        rep.unknowns = sel.len();
                        rep.filtered_out = full_rels.len() - sel.len();
                        break 'collect;
                    }
                    next_attempt = full_rels.len() + (sel.len() / 20).max(1);
                } else {
                    next_attempt = full_rels.len() + (unknowns / 20).max(1);
                }
            }
        }
    }
    // The walk's own operations: everything on the main curve since the
    // walk started, less the verification and the final check.
    let main_ops = curve.ops() - walk_start;
    rep.walk_group_ops = main_ops
        .saturating_sub(rep.verify_group_ops)
        .saturating_sub(64);
    rep.walk_muls = (rep.walk_group_ops as f64 * fp_per_add) as u64;
    rep.sign_muls = (rep.sign_group_ops as f64 * fp_per_add) as u64;
    rep.verify_muls = (rep.verify_group_ops as f64 * fp_per_add) as u64;
    rep.oracle_muls = rep.weil_muls + rep.f4_muls + rep.roots_muls + rep.sign_muls;
    rep.c_prime = rep.oracle_muls as f64 / rep.residuals.max(1) as f64;
    rep.la_ops = la_ops;
    rep.la_muls = la_ops * 16;
    rep.decomposition_rate = rep.relations as f64 / rep.residuals.max(1) as f64;
    rep.residual_ratio = rep.residuals as f64 / rep.floor_residuals;
    rep.phi = rep.unknowns as f64 / small as f64;
    rep.row_weight = full_rels.iter().map(|r| r.cols.len()).sum::<usize>() as f64
        / full_rels.len().max(1) as f64;
    rep.wiedemann_constant = la_ops as f64 / (rep.unknowns.max(1) as f64).powi(2);
    rep.total_muls = rep.setup_muls
        + rep.precompute_muls
        + rep.walk_muls
        + rep.oracle_muls
        + rep.verify_muls
        + rep.la_muls;
    let sqrt_n = (n as f64).sqrt();
    rep.s = rep.total_muls as f64 / (fp_per_add * sqrt_n);
    rep.rho = (0..rho_runs as u64)
        .into_par_iter()
        .map(|k| rho4(&spec, seed * 1000 + k))
        .collect();
    let ss: Vec<f64> = rep.rho.iter().map(|r| r.s).collect();
    rep.rho_s_mean = ss.iter().sum::<f64>() / ss.len().max(1) as f64;
    rep.rho_s_sd = (ss.iter().map(|s| (s - rep.rho_s_mean).powi(2)).sum::<f64>()
        / (ss.len().max(2) - 1) as f64)
        .sqrt();
    let rho_muls = rep.rho_s_mean * sqrt_n * fp_per_add;
    rep.s_over_rho = rep.total_muls as f64 / rho_muls;
    rep.relation_over_rho =
        (rep.setup_muls + rep.precompute_muls + rep.walk_muls + rep.oracle_muls + rep.verify_muls)
            as f64
            / rho_muls;
    rep.r = rep.la_muls as f64 / rho_muls;
    rep.predicted_s_over_rho = 6.0 * rep.c_prime / (rep.rho_s_mean * fp_per_add) + rep.r;
    rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
    rep
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small_setup() -> (
        super::super::gaudry_quartic::Instance4,
        Vec<Pt4>,
        HashMap<u64, usize>,
        SymmetrisedS4Q,
        StdRng,
    ) {
        let inst = generate_instance4(29, 5);
        let mut rng = StdRng::seed_from_u64(11);
        let base = factor_base(&inst.curve, &mut rng);
        let by_x: HashMap<u64, usize> = base
            .iter()
            .enumerate()
            .map(|(i, q)| (q.x.0[0], i))
            .collect();
        let pre = SymmetrisedS4Q::precompute(&inst.curve, &mut rng);
        (inst, base, by_x, pre, rng)
    }

    #[test]
    fn symmetrised_s4_vanishes_on_three_point_sums() {
        let (inst, base, _, pre, mut rng) = small_setup();
        let curve = &inst.curve;
        let f = &curve.f;
        assert_eq!(pre.monos.len(), 35);
        for _ in 0..8 {
            let i: [usize; 3] = core::array::from_fn(|_| rng.gen_range(0..base.len()));
            if i[0] == i[1] || i[1] == i[2] || i[0] == i[2] {
                continue;
            }
            let sum = curve.add(&curve.add(&base[i[0]], &base[i[1]]), &base[i[2]]);
            let e = elementary3(f, &[base[i[0]].x, base[i[1]].x, base[i[2]].x]);
            assert!(pre.eval(f, &e, &sum.x).is_zero());
            // and not at a random point
            let x: [E4; 3] = core::array::from_fn(|_| f.random(&mut rng));
            let e = elementary3(f, &x);
            assert!(!pre.eval(f, &e, &sum.x).is_zero());
        }
    }

    #[test]
    fn the_test_finds_planted_triples_and_agrees_with_the_oracle() {
        let (inst, base, by_x, pre, mut rng) = small_setup();
        let curve = &inst.curve;
        let table = pair_table(curve, &base);
        let opts = jv3_options();
        let mut planted = 0;
        for _ in 0..12 {
            let i: [usize; 3] = core::array::from_fn(|_| rng.gen_range(0..base.len()));
            if i[0] == i[1] || i[1] == i[2] || i[0] == i[2] {
                continue;
            }
            let s: [i64; 3] = core::array::from_fn(|_| if rng.gen_bool(0.5) { 1 } else { -1 });
            let mut r = Pt4::INF;
            for t in 0..3 {
                let q = if s[t] == 1 {
                    base[i[t]]
                } else {
                    curve.neg(&base[i[t]])
                };
                r = curve.add(&r, &q);
            }
            let mut want = [(i[0], s[0]), (i[1], s[1]), (i[2], s[2])];
            want.sort_unstable();
            let (found, cost) = jv3_decompose(curve, &base, &by_x, &pre, &r, &opts, &mut rng);
            assert!(!cost.undetermined, "{cost:?}");
            assert!(found.contains(&want), "planted {want:?} not in {found:?}");
            assert_eq!(found, mitm3_signed(curve, &base, &table, &r));
            assert!(found.iter().all(|t| verify_triple(curve, &base, &r, t)));
            planted += 1;
        }
        assert!(planted >= 8);
        // Random residuals: the oracle and the test agree (mostly "no").
        let mut inconsistent = 0;
        for _ in 0..40 {
            let a = rng.gen_range(0..inst.n);
            let b = rng.gen_range(1..inst.n);
            let r = curve.add(&curve.mul(&inst.g, a), &curve.mul(&inst.q, b));
            let (found, cost) = jv3_decompose(curve, &base, &by_x, &pre, &r, &opts, &mut rng);
            assert_eq!(found, mitm3_signed(curve, &base, &table, &r));
            if cost.inconsistent {
                inconsistent += 1;
            }
            assert!(cost.f4_muls > 0 && cost.weil_muls > 0);
        }
        assert!(inconsistent >= 20, "{inconsistent}");
    }

    #[test]
    fn cprime_report_is_consistent() {
        let rep = run_jv4_cprime(29, 5, 30, 6);
        assert_eq!(rep.mismatches, 0, "{rep:?}");
        assert_eq!(rep.planted_found, 6, "{rep:?}");
        assert_eq!(rep.undetermined, 0);
        assert!(rep.c_prime.mean > 0.0 && rep.f4_muls.mean > 0.0);
        assert_eq!(rep.fp_muls_per_add, 97.0);
    }

    #[test]
    fn the_method_recovers_the_planted_logarithm() {
        let rep = run_jv4_dlp(29, 3, 4, 8);
        assert!(rep.solved && rep.correct, "{rep:?}");
        assert_eq!(rep.mismatches, 0, "{rep:?}");
        assert!(rep.cross_checked > 0);
        assert!(rep.rho.iter().all(|r| r.correct));
        assert!(rep.residuals > 0 && rep.oracle_muls > 0 && rep.la_ops > 0);
        assert!(
            rep.phi > 0.5,
            "the core must carry most of the base: {}",
            rep.phi
        );
        assert!(
            (rep.row_weight - 4.0).abs() < 1e-9,
            "three points and d: {}",
            rep.row_weight
        );
        assert!(rep.s > 0.0 && rep.s_over_rho > 0.0);
        assert_eq!(
            rep.total_muls,
            rep.setup_muls
                + rep.precompute_muls
                + rep.walk_muls
                + rep.oracle_muls
                + rep.verify_muls
                + rep.la_muls
        );
    }
}
