//! # Folds of the base on the `k = 4` linear algebra (E18).
//!
//! `RESEARCH_RESIDUAL_WALKS.md` §11.19 measured the asymptotic ratio of a
//! plain `k = 4` index calculus on `E(F_{p⁴})` to rho: `r∞ = 0.518`, the
//! sparse linear algebra over `≈ p/2` unknowns against an unfolded rho
//! walk.  Past the handover the method beats rho by `1/r∞` and no more.
//! This module measures what a fold of the base does to that ratio
//! (`research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`
//! §9.1, registered before the code was written):
//!
//! * **`ψ` (`j = 0`).**  `y² = x³ + b` over `F_{p⁴}` with `b ∉ F_{p²}`,
//!   `p ≡ 1 (mod 12)`, prime order `n`.  `ψ(x, y) = (ωx, y)` with
//!   `ω ∈ F_p` of order 3 keeps the base `x ∈ F_p^*`, so its columns are
//!   the `⟨ψ⟩`-orbits of abscissae and a point's logarithm is `±λ^k` times
//!   its orbit representative's (`ψ = [λ]` on `⟨G⟩`).  The matched rho
//!   walks the classes `{±ψ^k P}`.
//! * **`τ_T` (2-torsion).**  `fghr_full`'s construction over `F_{p⁴}`:
//!   `x₀ ∉ F_{p²}`, `c ∈ F_p^*`, `a = c² − 3x₀²`, `b = −x₀³ − a x₀`, so
//!   `T = (x₀, 0)` and `x(P + T) − x₀ = c²/(x − x₀)`; `#E = h·n`.  The base
//!   `x ∈ x₀ + F_p` (offset `u ∉ {0, ±c}`) pairs `u` with `c²/u`; `P` and
//!   `P + T` share a column with coefficient `±1`.  Rho cannot use `τ_T`,
//!   so the matched rho is the negation walk.
//!
//! One relation stream — §11.19's meet-in-the-middle oracle, every
//! decomposition re-added in the group — feeds a control arm (one column
//! per abscissa, §11.19's matrix) and the folded arm.  Each arm is filtered
//! to a square core and solved by §11.7's Wiedemann, and its logarithm is
//! accepted only if `[d]G = Q`.  Rho runs unfolded, negation-folded and (on
//! `j = 0`) `ψ`-folded on one code path, every field multiplication
//! counted, so the three walks are priced alike.

use std::collections::HashMap;
use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::Serialize;

use super::gaudry_cubic::{square_core, wiedemann_u64, SparseRel};
use super::gaudry_quartic::{mitm_signed, pair_table, Curve4, Pt4, E4};
use super::glv_invariant_base::factor_u64;
use super::residual_walk::{inv_mod, is_prime_u64, mix64, pow_mod, sqrt_mod};

fn mulm(a: u64, b: u64, n: u64) -> u64 {
    ((a as u128 * b as u128) % n as u128) as u64
}

/// Which fold an instance is built for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum Family4 {
    /// `j = 0`, prime order, `ψ` of order 3.
    Psi,
    /// A rational 2-torsion point `T` whose translation keeps the base.
    Translation,
}

impl Family4 {
    pub fn name(self) -> &'static str {
        match self {
            Self::Psi => "j0-psi",
            Self::Translation => "2torsion-tau",
        }
    }
}

/// A known-answer instance over `F_{p⁴}`.
#[derive(Clone, Copy, Debug)]
pub struct FoldInstance4 {
    pub family: Family4,
    pub p: u64,
    pub a: E4,
    pub b: E4,
    pub group_order: u64,
    pub n: u64,
    pub cofactor: u64,
    pub g: Pt4,
    pub d: u64,
    pub q: Pt4,
    /// `ω ∈ F_p` of order 3 and `λ` with `ψ = [λ]` on `⟨G⟩` (`Psi`).
    pub omega: u64,
    pub lambda: u64,
    /// `T = (x₀, 0)`, `f'(x₀) = c²` (`Translation`).
    pub x0: E4,
    pub c: u64,
}

impl FoldInstance4 {
    /// A fresh curve with its own counters.
    pub fn curve(&self) -> Curve4 {
        Curve4::with_params(self.p, self.a, self.b)
    }

    pub fn psi(&self, curve: &Curve4, pt: &Pt4, k: u32) -> Pt4 {
        if pt.inf || k.is_multiple_of(3) {
            return *pt;
        }
        let w = pow_mod(self.omega, k as u64 % 3, self.p);
        Pt4 {
            x: curve.f.scale(&pt.x, w),
            y: pt.y,
            inf: false,
        }
    }
}

fn random_point(curve: &Curve4, rng: &mut StdRng) -> Pt4 {
    loop {
        let x = curve.f.random(rng);
        if let Some(pt) = curve.lift_x(&x, rng) {
            if !pt.y.is_zero() {
                return pt;
            }
        }
    }
}

/// The multiple of a random point's order in the Hasse interval
/// `p⁴ + 1 ± 2p²`, when it is unique there (it is `#E` once the order
/// exceeds the interval's width).
fn hasse_multiple(curve: &Curve4, rng: &mut StdRng) -> Option<u64> {
    let p = curve.f.p;
    let q = p.pow(4);
    let lo = q + 1 - 2 * p * p;
    let width = 4 * p * p;
    let pt = random_point(curve, rng);
    let steps = (width as f64).sqrt() as u64 + 1;
    let mut table: HashMap<(E4, E4), u64> = HashMap::new();
    let mut jp = Pt4::INF;
    for j in 0..steps {
        if !jp.inf {
            table.entry((jp.x, jp.y)).or_insert(j);
        }
        jp = curve.add(&jp, &pt);
    }
    let giant = curve.mul(&pt, steps);
    let mut t = curve.mul(&pt, lo);
    let mut found = None;
    let mut i = 0u64;
    while i * steps <= width + steps {
        let m = if t.inf {
            Some(lo + i * steps)
        } else {
            let neg = curve.neg(&t);
            table.get(&(neg.x, neg.y)).map(|&j| lo + i * steps + j)
        };
        if let Some(m) = m {
            if m <= lo + width {
                if found.is_some_and(|f| f != m) {
                    return None;
                }
                found.get_or_insert(m);
            }
        }
        t = curve.add(&t, &giant);
        i += 1;
    }
    found
}

fn not_in_fp2(a: &E4) -> bool {
    a.0[1] != 0 || a.0[3] != 0
}

/// A `j = 0` instance of prime order with `ψ = [λ]` verified on `G`.
pub fn generate_psi_instance(p: u64, seed: u64) -> Result<FoldInstance4, String> {
    if p % 12 != 1 || !is_prime_u64(p) {
        return Err(format!("p = {p} must be a prime ≡ 1 (mod 12)"));
    }
    let mut rng = StdRng::seed_from_u64(seed ^ 0x9510_4F0D);
    let omega = (2..p)
        .map(|g| pow_mod(g, (p - 1) / 3, p))
        .find(|&w| w != 1)
        .ok_or("no cube root of unity")?;
    for _ in 0..10_000 {
        let probe = Curve4::with_params(p, E4::ZERO, E4::ZERO);
        let b = probe.f.random(&mut rng);
        if b.is_zero() || !not_in_fp2(&b) {
            continue;
        }
        let curve = Curve4::with_params(p, E4::ZERO, b);
        let Some(n) = hasse_multiple(&curve, &mut rng) else {
            continue;
        };
        if !is_prime_u64(n) || n % 3 != 1 {
            continue;
        }
        let g = random_point(&curve, &mut rng);
        if !curve.mul(&g, n).inf {
            continue;
        }
        let s = sqrt_mod(n - 3, n).ok_or("−3 is not a square modulo n")?;
        let half = inv_mod(2, n);
        let roots = [
            mulm((n - 1 + s) % n, half, n),
            mulm((2 * n - 1 - s) % n, half, n),
        ];
        let mut inst = FoldInstance4 {
            family: Family4::Psi,
            p,
            a: E4::ZERO,
            b,
            group_order: n,
            n,
            cofactor: 1,
            g,
            d: 0,
            q: Pt4::INF,
            omega,
            lambda: 0,
            x0: E4::ZERO,
            c: 0,
        };
        let psi_g = inst.psi(&curve, &g, 1);
        if !curve.on_curve(&psi_g) {
            return Err("ψ(G) is not on the curve".into());
        }
        let Some(&lambda) = roots.iter().find(|&&l| curve.mul(&g, l) == psi_g) else {
            return Err("neither root of λ² + λ + 1 is ψ's eigenvalue".into());
        };
        inst.lambda = lambda;
        inst.d = rng.gen_range(1..n);
        inst.q = curve.mul(&g, inst.d);
        return Ok(inst);
    }
    Err(format!("no prime-order j = 0 curve found at p = {p}"))
}

/// A 2-torsion instance with `#E = h·n`, `h ≤ max_cofactor`, the smallest
/// `h` the search meets within its budget preferred.
pub fn generate_translation_instance(
    p: u64,
    seed: u64,
    max_cofactor: u64,
) -> Result<FoldInstance4, String> {
    if p % 4 != 1 || !is_prime_u64(p) {
        return Err(format!("p = {p} must be a prime ≡ 1 (mod 4)"));
    }
    let mut rng = StdRng::seed_from_u64(seed ^ 0x7A07_2B05);
    let mut fallback: Option<FoldInstance4> = None;
    for attempt in 0..4_000u32 {
        if attempt == 2_000 && fallback.is_some() {
            break;
        }
        let probe = Curve4::with_params(p, E4::ZERO, E4::ZERO);
        let f = &probe.f;
        let x0 = f.random(&mut rng);
        if !not_in_fp2(&x0) {
            continue;
        }
        let c = rng.gen_range(1..p);
        let x0sq = f.sq(&x0);
        let a = f.sub(&f.from_fp(c * c % p), &f.scale(&x0sq, 3));
        let b = f.neg(&f.add(&f.mul(&x0sq, &x0), &f.mul(&a, &x0)));
        let a3x4 = f.scale(&f.mul(&f.sq(&a), &a), 4);
        let disc = f.add(&a3x4, &f.scale(&f.sq(&b), 27));
        if disc.is_zero() {
            continue;
        }
        let j = f.mul(&a3x4, &f.inv(&disc));
        if !not_in_fp2(&j) {
            continue;
        }
        let curve = Curve4::with_params(p, a, b);
        let Some(order) = hasse_multiple(&curve, &mut rng) else {
            continue;
        };
        if order % 2 != 0 {
            return Err(format!(
                "#E = {order} is odd although T = (x₀, 0) is rational"
            ));
        }
        let Some(&(n, _)) = factor_u64(order).last() else {
            continue;
        };
        let h = order / n;
        if h > max_cofactor || order % (n * n) == 0 || n <= p * p * p {
            continue;
        }
        let g = loop {
            let g = curve.mul(&random_point(&curve, &mut rng), h);
            if !g.inf {
                break g;
            }
        };
        if !curve.mul(&g, n).inf {
            return Err("n·G ≠ O: the order is wrong".into());
        }
        let t = Pt4 {
            x: x0,
            y: E4::ZERO,
            inf: false,
        };
        if !curve.on_curve(&t) {
            return Err("T is not on the curve".into());
        }
        let d = rng.gen_range(1..n);
        let inst = FoldInstance4 {
            family: Family4::Translation,
            p,
            a,
            b,
            group_order: order,
            n,
            cofactor: h,
            g,
            d,
            q: curve.mul(&g, d),
            omega: 0,
            lambda: 0,
            x0,
            c,
        };
        if h == 2 {
            return Ok(inst);
        }
        if fallback.as_ref().is_none_or(|f| h < f.cofactor) {
            fallback = Some(inst);
        }
    }
    fallback.ok_or_else(|| format!("no 2-torsion instance with h ≤ {max_cofactor} at p = {p}"))
}

/// The base, the same points for both arms, and the folded arm's column
/// map: `log P_i = coef_i · X_{col_i}`.
#[derive(Clone, Debug)]
pub struct FoldedBase4 {
    pub points: Vec<Pt4>,
    pub col: Vec<usize>,
    pub coef: Vec<u64>,
    pub columns: usize,
    /// Points per column, mean.
    pub points_per_column: f64,
}

/// Build the base and check every column coefficient on points.
pub fn folded_base(inst: &FoldInstance4, rng: &mut StdRng) -> Result<FoldedBase4, String> {
    let curve = inst.curve();
    let f = &curve.f;
    let p = inst.p;
    let n = inst.n;
    let mut points = Vec::new();
    let mut key: Vec<u64> = Vec::new(); // the F_p parameter of each point
    match inst.family {
        Family4::Psi => {
            for x in 1..p {
                if let Some(pt) = curve.lift_x(&f.from_fp(x), rng) {
                    points.push(pt);
                    key.push(x);
                }
            }
        }
        Family4::Translation => {
            for u in 1..p {
                if u == inst.c || u == p - inst.c {
                    continue;
                }
                let x = f.add(&inst.x0, &f.from_fp(u));
                if let Some(pt) = curve.lift_x(&x, rng) {
                    points.push(pt);
                    key.push(u);
                }
            }
        }
    }
    let index: HashMap<u64, usize> = key.iter().enumerate().map(|(i, &k)| (k, i)).collect();
    let mut col = vec![usize::MAX; points.len()];
    let mut coef = vec![0u64; points.len()];
    let mut columns = 0usize;
    for i in 0..points.len() {
        if col[i] != usize::MAX {
            continue;
        }
        // i is a fresh representative: walk its orbit.
        let rep = points[i];
        col[i] = columns;
        coef[i] = 1;
        let images: Vec<(u64, Pt4, u64)> = match inst.family {
            Family4::Psi => (1..3u32)
                .map(|k| {
                    let w = pow_mod(inst.omega, k as u64, p);
                    (
                        mulm(key[i], w, p),
                        inst.psi(&curve, &rep, k),
                        pow_mod(inst.lambda, k as u64, n),
                    )
                })
                .collect(),
            Family4::Translation => {
                let t = Pt4 {
                    x: inst.x0,
                    y: E4::ZERO,
                    inf: false,
                };
                let u2 = mulm(mulm(inst.c, inst.c, p), inv_mod(key[i], p), p);
                vec![(u2, curve.add(&rep, &t), 1)]
            }
        };
        for (k2, img, mult) in images {
            let Some(&j) = index.get(&k2) else {
                return Err(format!("the orbit of {} leaves the base", key[i]));
            };
            if j == i {
                continue;
            }
            if points[j].x != img.x {
                return Err("an image's abscissa is not the base point's".into());
            }
            let s = if points[j].y == img.y {
                mult
            } else if points[j].y == f.neg(&img.y) {
                (n - mult) % n
            } else {
                return Err("an image is neither P nor −P".into());
            };
            if col[j] != usize::MAX && (col[j] != columns || coef[j] != s) {
                return Err("an orbit closes inconsistently".into());
            }
            col[j] = columns;
            coef[j] = s;
        }
        columns += 1;
    }
    let ppc = points.len() as f64 / columns.max(1) as f64;
    Ok(FoldedBase4 {
        points,
        col,
        coef,
        columns,
        points_per_column: ppc,
    })
}

/// One arm's matrix and its solve.
#[derive(Clone, Debug, Default, Serialize)]
pub struct ArmReport {
    pub arm: String,
    pub columns: usize,
    pub relations: usize,
    pub residuals_at_solve: u64,
    pub unknowns: usize,
    pub filtered_out: usize,
    pub phi: f64,
    pub row_weight: f64,
    pub la_attempts: u64,
    /// Multiplications mod `n`, every attempt.
    pub la_ops: u64,
    pub la_ops_last: u64,
    pub wiedemann_constant: f64,
    pub solved: bool,
    pub correct: bool,
    /// Rows identical to one already kept, skipped.  On the `ψ` arm two
    /// decompositions of one target that differ by `P + ψP + ψ²P = O`
    /// fold to the same row; the control keeps them in distinct columns.
    pub duplicate_rows: u64,
}

struct Arm {
    rep: ArmReport,
    rels: Vec<SparseRel>,
    seen: std::collections::HashSet<(Vec<(usize, u64)>, u64)>,
    small: usize,
    next_attempt: usize,
}

impl Arm {
    fn new(name: &str, columns: usize) -> Self {
        Arm {
            rep: ArmReport {
                arm: name.into(),
                columns,
                ..Default::default()
            },
            rels: Vec::new(),
            seen: Default::default(),
            small: columns,
            next_attempt: columns + 1,
        }
    }

    /// Keep a row unless an identical one is already in.
    fn push(&mut self, row: SparseRel) {
        if self.seen.insert((row.cols.clone(), row.rhs)) {
            self.rels.push(row);
        } else {
            self.rep.duplicate_rows += 1;
        }
    }

    /// Try the solve once enough rows are in.
    fn attempt(&mut self, inst: &FoldInstance4, curve: &Curve4, rng: &mut StdRng, residuals: u64) {
        if self.rep.solved || self.rels.len() < self.next_attempt {
            return;
        }
        let unknowns = self.small + 1;
        let n = inst.n;
        let Some(core) = square_core(&self.rels, self.small) else {
            self.next_attempt = self.rels.len() + (unknowns / 20).max(1);
            return;
        };
        self.rep.la_attempts += 1;
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
                cols: self.rels[i]
                    .cols
                    .iter()
                    .map(|&(c, v)| (map[c], v))
                    .collect(),
                rhs: self.rels[i].rhs,
            })
            .collect();
        let before = self.rep.la_ops;
        let x = wiedemann_u64(&sel, sel.len(), n, rng, &mut self.rep.la_ops);
        let dd = x.map(|x| x[map[self.small]]);
        if dd.is_some_and(|dd| curve.mul(&inst.g, dd) == inst.q) {
            self.rep.solved = true;
            self.rep.correct = dd == Some(inst.d);
            self.rep.unknowns = sel.len();
            self.rep.filtered_out = self.rels.len() - sel.len();
            self.rep.la_ops_last = self.rep.la_ops - before;
            self.rep.relations = self.rels.len();
            self.rep.residuals_at_solve = residuals;
            self.rep.phi = sel.len() as f64 / self.small as f64;
            self.rep.row_weight = self.rels.iter().map(|r| r.cols.len()).sum::<usize>() as f64
                / self.rels.len().max(1) as f64;
            self.rep.wiedemann_constant =
                self.rep.la_ops as f64 / (self.rep.unknowns.max(1) as f64).powi(2);
        } else {
            self.next_attempt = self.rels.len() + (sel.len() / 20).max(1);
        }
    }
}

/// Which classes a rho walk folds by.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum RhoFold {
    None,
    Negation,
    /// `{±ψ^k P}` (`j = 0` instances only).
    Psi,
}

#[derive(Clone, Debug, Serialize)]
pub struct RhoRun4F {
    pub fold: RhoFold,
    pub seed: u64,
    pub steps: u64,
    pub group_ops: u64,
    /// Every `F_p` multiplication the walk spent, set-up and the class
    /// representative included (an `F_{p⁴}` inversion at `40`).
    pub fp_muls: u64,
    /// Fruitless cycles met and escaped by doubling.
    pub cycle_escapes: u64,
    pub correct: bool,
}

/// An r-adding walk with exhaustive storage on class representatives.
/// A step whose successor would land in its own partition takes the next
/// partition instead (the look-ahead that keeps a folded walk out of
/// fruitless 2-cycles); a self-collision — a longer fruitless cycle — is
/// escaped by doubling the current point.  Every group operation and field
/// multiplication is charged, set-up included.
pub fn rho4_folded(inst: &FoldInstance4, fold: RhoFold, seed: u64) -> RhoRun4F {
    assert!(fold != RhoFold::Psi || inst.family == Family4::Psi);
    let curve = inst.curve();
    let f = &curve.f;
    let n = inst.n;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x8D04_F01D);
    let r = 32usize;
    let lam = [1, inst.lambda % n, mulm(inst.lambda, inst.lambda, n)];
    let canon = |pt: Pt4| -> (Pt4, u64) {
        match fold {
            RhoFold::None => (pt, 1),
            RhoFold::Negation => {
                let ny = f.neg(&pt.y);
                if ny < pt.y {
                    (Pt4 { y: ny, ..pt }, n - 1)
                } else {
                    (pt, 1)
                }
            }
            RhoFold::Psi => {
                let mut best = (pt, 1u64);
                for (k, &l) in lam.iter().enumerate().skip(1) {
                    let img = inst.psi(&curve, &pt, k as u32);
                    if img.x < best.0.x {
                        best = (img, l);
                    }
                }
                let ny = f.neg(&best.0.y);
                if ny < best.0.y {
                    (Pt4 { y: ny, ..best.0 }, (n - best.1) % n)
                } else {
                    best
                }
            }
        }
    };
    let part = |pt: &Pt4| -> usize {
        let h = pt.x.0[0] ^ pt.x.0[1].wrapping_mul(0x9E37_79B9) ^ pt.x.0[2].rotate_left(17);
        (mix64(h ^ pt.x.0[3].rotate_left(29)) % r as u64) as usize
    };
    let mults: Vec<(u64, u64, Pt4)> = (0..r)
        .map(|_| {
            let al = rng.gen_range(0..n);
            let be = rng.gen_range(1..n);
            (
                al,
                be,
                curve.add(&curve.mul(&inst.g, al), &curve.mul(&inst.q, be)),
            )
        })
        .collect();
    let mut table: HashMap<Pt4, (u64, u64)> = HashMap::new();
    let (mut steps, mut cycle_escapes) = (0u64, 0u64);
    let mut found = None;
    let cap = 64 * (n as f64).sqrt() as u64 + 1_000_000;
    let mut a = rng.gen_range(0..n);
    let mut b = rng.gen_range(1..n);
    let start = curve.add(&curve.mul(&inst.g, a), &curve.mul(&inst.q, b));
    let (mut l, m) = canon(start);
    a = mulm(a, m, n);
    b = mulm(b, m, n);
    while steps <= cap {
        steps += 1;
        if let Some(&(a2, b2)) = table.get(&l) {
            let db = (b + n - b2) % n;
            if db != 0 {
                let da = (a2 + n - a) % n;
                found = Some(mulm(da, inv_mod(db, n), n));
                break;
            }
            // A fruitless cycle: escape by doubling, which costs one group
            // operation where a restart costs two scalar multiplications.
            cycle_escapes += 1;
            let (dbl, m) = canon(curve.add(&l, &l));
            l = dbl;
            a = mulm(mulm(2, a, n), m, n);
            b = mulm(mulm(2, b, n), m, n);
            continue;
        }
        table.insert(l, (a, b));
        let mut j = part(&l);
        let (mut next, mut m) = canon(curve.add(&l, &mults[j].2));
        if fold != RhoFold::None && part(&next) == j {
            j = (j + 1) % r;
            (next, m) = canon(curve.add(&l, &mults[j].2));
        }
        a = mulm((a + mults[j].0) % n, m, n);
        b = mulm((b + mults[j].1) % n, m, n);
        l = next;
    }
    RhoRun4F {
        fold,
        seed,
        steps,
        group_ops: curve.ops(),
        fp_muls: f.muls(),
        cycle_escapes,
        correct: found == Some(inst.d),
    }
}

/// One E18 run: the instance, both arms on one relation stream, and the
/// matched walks.
#[derive(Clone, Debug, Serialize)]
pub struct K4FoldReport {
    pub family: Family4,
    pub p: u64,
    pub seed: u64,
    pub n: u64,
    pub bits: f64,
    pub group_order: u64,
    pub cofactor: u64,
    pub base_points: usize,
    pub points_per_column: f64,
    pub residuals: u64,
    /// Every verified decomposition the oracle returned.
    pub decompositions: u64,
    /// Targets with at least one; each contributes one decomposition, drawn
    /// uniformly from its verified ones, to both arms.
    pub targets_decomposed: u64,
    /// How many verified decompositions targets had (`[k]` = targets with
    /// `k`, the last bucket `≥ 9`).  On a curve with rational `T` one
    /// decomposition comes with the `7` others that add `T` to an even
    /// number of summands, and the eight rows have rank `5`.
    pub decompositions_per_target: [u64; 10],
    pub decomposition_rate: f64,
    pub fp_muls_per_add: f64,
    pub control: ArmReport,
    pub folded: ArmReport,
    pub rho: Vec<RhoRun4F>,
    pub wall_ms: f64,
}

pub fn run_k4_folds(
    family: Family4,
    p: u64,
    seed: u64,
    rho_runs: usize,
) -> Result<K4FoldReport, String> {
    let start = Instant::now();
    let inst = match family {
        Family4::Psi => generate_psi_instance(p, seed)?,
        Family4::Translation => generate_translation_instance(p, seed, 4)?,
    };
    let curve = inst.curve();
    let n = inst.n;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x1A4F_01D5);
    let fb = folded_base(&inst, &mut rng)?;
    let base = &fb.points;
    let fp_per_add = {
        curve.f.reset_muls();
        let mut acc = base[0];
        for _ in 0..64 {
            acc = curve.add(&acc, &base[1]);
        }
        curve.f.muls() as f64 / 64.0
    };
    let table = pair_table(&curve, base);
    let mut control = Arm::new("control", base.len());
    let mut folded = Arm::new(family.name(), fb.columns);
    let (mut residuals, mut decompositions, mut targets_decomposed) = (0u64, 0u64, 0u64);
    let mut per_target = [0u64; 10];
    while residuals < 50_000_000 && !(control.rep.solved && folded.rep.solved) {
        let batch: Vec<(u64, u64)> = (0..64)
            .map(|_| (rng.gen_range(0..n), rng.gen_range(1..n)))
            .collect();
        let found: Vec<(u64, u64, Pt4, Vec<Vec<(usize, i64)>>)> = batch
            .par_iter()
            .map(|&(a, b)| {
                let c = inst.curve();
                let r = c.add(&c.mul(&inst.g, a), &c.mul(&inst.q, b));
                let decs = mitm_signed(&c, base, &table, &r);
                (a, b, r, decs)
            })
            .collect();
        for (a, b, r, decs) in found {
            residuals += 1;
            let mut good: Vec<Vec<(usize, i64)>> = Vec::new();
            for terms in decs {
                // Verify aG + bQ = Σ s P before using it.
                let mut acc = r;
                for &(i, s) in &terms {
                    acc = if s == 1 {
                        curve.sub(&acc, &base[i])
                    } else {
                        curve.add(&acc, &base[i])
                    };
                }
                if !acc.inf {
                    continue;
                }
                decompositions += 1;
                good.push(terms);
            }
            per_target[good.len().min(9)] += 1;
            if !good.is_empty() {
                targets_decomposed += 1;
                // One row per target, drawn uniformly from its verified
                // decompositions: the first in the oracle's sorted order would
                // prefer the lower index of every `{P, P + T}` pair and so half
                // fold the control.
                let terms = good.swap_remove(rng.gen_range(0..good.len()));
                let qcol = |small: usize| (small, (n - b) % n);
                if !control.rep.solved {
                    let mut cols: Vec<(usize, u64)> = terms
                        .iter()
                        .map(|&(i, s)| (i, if s == 1 { 1 } else { n - 1 }))
                        .collect();
                    cols.push(qcol(control.small));
                    cols.sort_unstable();
                    control.push(SparseRel { cols, rhs: a });
                }
                if !folded.rep.solved {
                    let mut merged: HashMap<usize, u64> = HashMap::new();
                    for &(i, s) in &terms {
                        let v = if s == 1 {
                            fb.coef[i]
                        } else {
                            (n - fb.coef[i]) % n
                        };
                        let e = merged.entry(fb.col[i]).or_insert(0);
                        *e = (*e + v) % n;
                    }
                    let mut cols: Vec<(usize, u64)> =
                        merged.into_iter().filter(|&(_, v)| v != 0).collect();
                    cols.push(qcol(folded.small));
                    cols.sort_unstable();
                    folded.push(SparseRel { cols, rhs: a });
                }
            }
            control.attempt(&inst, &curve, &mut rng, residuals);
            folded.attempt(&inst, &curve, &mut rng, residuals);
        }
    }
    let mut folds = vec![RhoFold::None, RhoFold::Negation];
    if family == Family4::Psi {
        folds.push(RhoFold::Psi);
    }
    let jobs: Vec<(RhoFold, u64)> = folds
        .iter()
        .flat_map(|&f| (0..rho_runs as u64).map(move |k| (f, seed * 1000 + k)))
        .collect();
    let rho: Vec<RhoRun4F> = jobs
        .par_iter()
        .map(|&(f, s)| rho4_folded(&inst, f, s))
        .collect();
    Ok(K4FoldReport {
        family,
        p,
        seed,
        n,
        bits: (n as f64).log2(),
        group_order: inst.group_order,
        cofactor: inst.cofactor,
        base_points: base.len(),
        points_per_column: fb.points_per_column,
        residuals,
        decompositions,
        targets_decomposed,
        decompositions_per_target: per_target,
        decomposition_rate: decompositions as f64 / residuals.max(1) as f64,
        fp_muls_per_add: fp_per_add,
        control: control.rep,
        folded: folded.rep,
        rho,
        wall_ms: start.elapsed().as_secs_f64() * 1e3,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn psi_columns_are_orbits_of_three_with_verified_coefficients() {
        let inst = generate_psi_instance(37, 1).unwrap();
        let mut rng = StdRng::seed_from_u64(3);
        let fb = folded_base(&inst, &mut rng).unwrap();
        assert_eq!(fb.points_per_column, 3.0);
        let curve = inst.curve();
        assert_eq!(
            curve.mul(&inst.g, inst.lambda),
            inst.psi(&curve, &inst.g, 1)
        );
    }

    #[test]
    fn translation_columns_pair_p_and_p_plus_t() {
        let inst = generate_translation_instance(37, 2, 4).unwrap();
        let mut rng = StdRng::seed_from_u64(4);
        let fb = folded_base(&inst, &mut rng).unwrap();
        assert!(fb.points_per_column > 1.9 && fb.points_per_column <= 2.0);
        assert_eq!(inst.group_order, inst.cofactor * inst.n);
    }

    #[test]
    fn both_arms_recover_the_logarithm_and_every_walk_is_correct() {
        for family in [Family4::Psi, Family4::Translation] {
            let rep = run_k4_folds(family, 37, 1, 2).unwrap();
            assert!(
                rep.control.solved && rep.control.correct,
                "{family:?} control"
            );
            assert!(rep.folded.solved && rep.folded.correct, "{family:?} folded");
            assert!(rep.rho.iter().all(|w| w.correct), "{family:?} rho");
            assert!(rep.folded.columns * 2 <= rep.control.columns + 2);
        }
    }
}

#[cfg(test)]
mod log_tests {
    use super::*;

    fn dlog(curve: &Curve4, g: &Pt4, h: &Pt4, n: u64) -> u64 {
        let m = (n as f64).sqrt() as u64 + 1;
        let mut baby: HashMap<Pt4, u64> = HashMap::new();
        let mut cur = Pt4::INF;
        for j in 0..m {
            baby.entry(cur).or_insert(j);
            cur = curve.add(&cur, g);
        }
        let giant = curve.neg(&curve.mul(g, m));
        let mut t = *h;
        for i in 0..=m {
            if let Some(&j) = baby.get(&t) {
                return (i * m + j) % n;
            }
            t = curve.add(&t, &giant);
        }
        panic!("no logarithm");
    }

    /// Every base point's logarithm (of its `⟨G⟩`-component, `h⁻¹·log(hP)`)
    /// is its column coefficient times its representative's.
    #[test]
    fn fold_coefficients_match_true_logarithms() {
        for inst in [
            generate_psi_instance(37, 1).unwrap(),
            generate_translation_instance(37, 2, 4).unwrap(),
        ] {
            let curve = inst.curve();
            let n = inst.n;
            let mut rng = StdRng::seed_from_u64(9);
            let fb = folded_base(&inst, &mut rng).unwrap();
            let hinv = inv_mod(inst.cofactor % n, n);
            let logs: Vec<u64> = fb
                .points
                .iter()
                .map(|p| {
                    let hp = curve.mul(p, inst.cofactor);
                    mulm(dlog(&curve, &inst.g, &hp, n), hinv, n)
                })
                .collect();
            let mut rep = vec![usize::MAX; fb.columns];
            for i in 0..fb.points.len() {
                if fb.coef[i] == 1 && rep[fb.col[i]] == usize::MAX {
                    rep[fb.col[i]] = i;
                }
            }
            for i in 0..fb.points.len() {
                let r = rep[fb.col[i]];
                assert_eq!(
                    mulm(fb.coef[i], logs[r], n),
                    logs[i],
                    "{:?} point {i}",
                    inst.family
                );
            }
        }
    }
}

#[cfg(test)]
mod rho_tests {
    use super::*;

    /// Each walk takes the steps its class size predicts: `√(π n / 2w)`,
    /// `w = 1, 2, 6`, so the unfolded walk takes `√2` and `√6` times the
    /// steps of the negation and `ψ` walks (bands for `80` walks each).
    #[test]
    fn folded_walks_take_the_steps_their_classes_predict() {
        let inst = generate_psi_instance(97, 1).unwrap();
        let mean_steps = |fold: RhoFold| -> f64 {
            let runs: Vec<RhoRun4F> = (0..80u64)
                .into_par_iter()
                .map(|s| rho4_folded(&inst, fold, s))
                .collect();
            assert!(runs.iter().all(|r| r.correct), "{fold:?}");
            runs.iter().map(|r| r.steps as f64).sum::<f64>() / runs.len() as f64
        };
        let none = mean_steps(RhoFold::None);
        let neg = none / mean_steps(RhoFold::Negation);
        let psi = none / mean_steps(RhoFold::Psi);
        assert!((1.15..1.75).contains(&neg), "unfolded / negation = {neg}");
        assert!((1.9..3.1).contains(&psi), "unfolded / psi = {psi}");
    }
}
