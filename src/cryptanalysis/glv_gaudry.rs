//! **GLV / `C₃` experiments on the `E(F_{p³})` index-calculus harness.**
//!
//! Companion to `RESEARCH_GLV_INDEX_CALCULUS.md`.  Four experiments on
//! `j = 0` curves `y² = x³ + b` over `F_{p³}`, `p ≡ 1 (mod 3)`, where the
//! order-3 automorphism `ψ(x, y) = (ωx, y)` (`ω ∈ F_p` a primitive cube
//! root of unity) acts on the prime-order group as multiplication by
//! `λ`, `λ² + λ + 1 ≡ 0 (mod n)`.  Because `ω ∈ F_p`, `ψ` preserves the
//! subspace factor base `F = {P : x(P) ∈ F_p}` of `gaudry_cubic`, so
//! the base splits into exact `ψ`-orbits of size three and every base
//! point is `ψ^k` of its orbit representative, i.e. `λ^k` times it.
//!
//! 1. [`run_quotient_experiment`] — the **factor-base quotient**: one
//!    column per `⟨ψ⟩`-orbit, orientation encoded as `λ^k`, measured
//!    against the ordinary `⟨−1⟩`-folded base on the *same residual
//!    stream* (relations, columns, non-zeros, Wiedemann time, memory).
//! 2. [`run_canonical_experiment`] — **canonical relation generation**:
//!    residuals and decompositions are reduced to a canonical
//!    representative of their `⟨−1, ψ⟩`-orbit before the solver and the
//!    relation store see them; duplicates and solver calls saved are
//!    counted on a random stream and on a pair-enumeration stream.
//! 3. [`run_graded_experiment`] — the **`Z/3`-graded Macaulay matrix**:
//!    monomials partitioned by their `ψ`-weight `Σ wᵢeᵢ (mod 3)`,
//!    block-aware elimination against the plain matrix, on the ordinary
//!    per-residual system and on the `ψ`-homogeneous orbit system.
//! 4. [`run_invariant_experiment`] — the **`C₃`-invariant formulation**:
//!    invariant monomials and generators enumerated, the PDP rewritten
//!    in invariants of the diagonal action, and compared with the
//!    ordinary symmetrised Semaev `S₄` solve and with a prime-field
//!    function-first (Nagao `L(4O)` coefficient) system solved by F4.
//!
//! Costs are in the harness's unit (`F_p` multiplications, converted to
//! affine additions at the measured rate) so rows are comparable with
//! `RESEARCH_RESIDUAL_WALKS.md` §11 and with rho on the same group.

use std::collections::{HashMap, HashSet};
use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use super::f4_fp::{self, F4Options, Ordering as F4Ordering, Verdict};
use super::gaudry_cubic::{
    run_rho3, s4_terms, solve_s4_subspace, square_core, subspace_triple_oracle_groebner,
    unique_hasse_multiple3, wiedemann_u64, Curve3, Fp3, Instance3, Pt3, RhoReport3, SolveStats,
    SparseRel, SubspaceBase, SymmetrisedS4, E3,
};
use super::residual_walk::{inv_mod, is_prime_u64, pow_mod, sqrt_mod};

// ── F_p helpers ─────────────────────────────────────────────────────────

#[inline]
fn mm(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}
#[inline]
fn am(a: u64, b: u64, p: u64) -> u64 {
    let s = a + b;
    if s >= p {
        s - p
    } else {
        s
    }
}
#[inline]
fn sm(a: u64, b: u64, p: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a + p - b
    }
}

/// Peak resident set size of this process, in bytes (`getrusage`); 0 on
/// targets without `getrusage`.
pub fn peak_rss_bytes() -> u64 {
    #[cfg(unix)]
    // SAFETY: getrusage writes into the zeroed struct we hand it.
    unsafe {
        let mut ru: libc::rusage = std::mem::zeroed();
        libc::getrusage(libc::RUSAGE_SELF, &mut ru);
        let unit = if cfg!(target_os = "macos") { 1 } else { 1024 };
        (ru.ru_maxrss as u64) * unit
    }
    #[cfg(not(unix))]
    0
}

// ── j = 0 instances over F_{p³} ─────────────────────────────────────────

/// A known-answer instance on `y² = x³ + b` over `F_{p³}` with its
/// order-3 automorphism `ψ(x, y) = (ωx, y) = [λ](x, y)`.
#[derive(Clone, Debug, Serialize)]
pub struct GlvInstance3 {
    pub inst: Instance3,
    /// Primitive cube root of unity in `F_p`.
    pub omega: u64,
    /// `ψ = [λ]` on the group, `λ² + λ + 1 ≡ 0 (mod n)`.
    pub lambda: u64,
}

impl GlvInstance3 {
    pub fn omega_pow(&self, k: u8) -> u64 {
        pow_mod(self.omega, k as u64 % 3, self.inst.curve.field.p)
    }
    pub fn lambda_pow(&self, k: u8) -> u64 {
        pow_mod(self.lambda, k as u64 % 3, self.inst.curve.n)
    }
    /// `ψ(P) = (ωx, y)`.
    pub fn psi(&self, pt: &Pt3) -> Pt3 {
        self.psi_k(pt, 1)
    }
    /// `ψ^k(P) = (ω^k x, y)`.
    pub fn psi_k(&self, pt: &Pt3, k: u8) -> Pt3 {
        if pt.inf {
            return *pt;
        }
        let f = &self.inst.curve.field;
        Pt3::affine(f.scale(&pt.x, self.omega_pow(k)), pt.y)
    }
    /// The unit `u = ±λ^k` with `P = u · C` for the canonical
    /// representative `C` of the `⟨−1, ψ⟩`-orbit of `P`, and `C`.
    pub fn canonical(&self, pt: &Pt3) -> (Pt3, u64) {
        let f = &self.inst.curve.field;
        let n = self.inst.curve.n;
        if pt.inf {
            return (*pt, 1);
        }
        // Smallest x among {x, ωx, ω²x}, then the smaller y.
        let mut best_k = 0u8;
        let mut best_x = pt.x;
        for k in 1..3u8 {
            let xk = f.scale(&pt.x, self.omega_pow(k));
            if xk.0 < best_x.0 {
                best_x = xk;
                best_k = k;
            }
        }
        // C = ψ^{best_k}(σP) with σ = ±1 chosen so that y is small.
        let (y, sigma) = if f.is_small(&pt.y) {
            (pt.y, 1u64)
        } else {
            (f.neg(&pt.y), n - 1)
        };
        // P = σ · ψ^{−best_k}(C) = σ · λ^{(3 − best_k) mod 3} · C.
        let inv_k = (3 - best_k) % 3;
        let u = mm(sigma, self.lambda_pow(inv_k), n);
        (Pt3::affine(best_x, y), u)
    }
    /// Six-word key of the canonical representative.
    pub fn canonical_key(&self, pt: &Pt3) -> [u64; 6] {
        let (c, _) = self.canonical(pt);
        [c.x.0[0], c.x.0[1], c.x.0[2], c.y.0[0], c.y.0[1], c.y.0[2]]
    }
}

/// Random prime-order `j = 0` curve `y² = x³ + b` over `F_{p³}` for a
/// prime `p ≡ 1 (mod 3)`, `p < 2^20`, with `ω`, `λ` and a random
/// known-answer target.  Prime order forces `n ≡ 1 (mod 3)` (the
/// automorphism acts as a scalar of order 3) and `b` a non-square in
/// `F_{p³}` (no point with `x = 0`, which `ψ` would fix), so every
/// subspace-base orbit has exactly three points.
///
/// Unlike a generic curve, a `j = 0` curve over `F_{p³}` has only the
/// six sextic-twist orders to choose from, fixed by `p`; for many `p`
/// none is prime (`p = 67`: none; `p = 73, 271, 541, 1051, 2113`: one
/// each), so the search over `b` is bounded and panics, naming the
/// problem, when `p` admits no prime-order twist.
pub fn generate_j0_instance3(p: u64, seed: u64) -> GlvInstance3 {
    assert!(p % 3 == 1 && is_prime_u64(p) && p < (1 << 20));
    let mut rng = StdRng::seed_from_u64(seed ^ 0x6A0D_0000_4B17);
    let mut attempts = 0u32;
    let omega = loop {
        let g = rng.gen_range(2..p);
        let w = pow_mod(g, (p - 1) / 3, p);
        if w != 1 {
            break w;
        }
    };
    loop {
        attempts += 1;
        assert!(
            attempts <= 240,
            "no prime-order j = 0 twist over F_{{{p}³}} found in 240 draws of b: \
             the six sextic-twist orders are fixed by p, pick another p ≡ 1 (mod 3)"
        );
        let field = Fp3::new(p, &mut rng);
        let rnd = |rng: &mut StdRng| {
            E3([
                rng.gen_range(0..p),
                rng.gen_range(0..p),
                rng.gen_range(0..p),
            ])
        };
        let b = rnd(&mut rng);
        if b == Fp3::ZERO {
            continue;
        }
        let curve = Curve3::new(field.clone(), Fp3::ZERO, b, 0, Pt3::INFINITY);
        let g = loop {
            let x = rnd(&mut rng);
            if let Some(pt) = curve.lift(&x, &mut rng) {
                break pt;
            }
        };
        let Some(order) = unique_hasse_multiple3(&curve, &g) else {
            continue;
        };
        if !is_prime_u64(order) || order % 3 != 1 {
            continue;
        }
        let curve = Curve3::new(field, Fp3::ZERO, b, order, g);
        assert!(curve.mul(&g, order).inf);
        let inst = Instance3 {
            curve,
            q: Pt3::INFINITY,
            d: 0,
        };
        let mut glv = GlvInstance3 {
            inst,
            omega,
            lambda: 1,
        };
        // λ = (−1 ± √−3) / 2 mod n; the root that matches ψ on g.
        let root = sqrt_mod(order - 3, order).expect("−3 is a square mod n when n ≡ 1 mod 3");
        let half = inv_mod(2, order);
        let psi_g = glv.psi(&g);
        let lambda = [root, order - root]
            .into_iter()
            .map(|r| mm((r + order - 1) % order, half, order))
            .find(|&l| glv.inst.curve.mul(&g, l) == psi_g)
            .expect("one root of λ² + λ + 1 matches ψ");
        glv.lambda = lambda;
        let d = rng.gen_range(1..order);
        let q = glv.inst.curve.mul(&g, d);
        glv.inst.q = q;
        glv.inst.d = d;
        glv.inst.curve.reset();
        return glv;
    }
}

// ── Orbit quotient of the subspace base ─────────────────────────────────

/// The subspace base modulo `⟨ψ⟩`: one representative per orbit (the
/// smallest `x`), and for every base point its orbit and orientation
/// `k` with `x = ω^k · x_rep`, so that `P = ψ^k(P_rep) = [λ^k] P_rep`.
#[derive(Clone, Debug)]
pub struct OrbitBase {
    pub reps: Vec<usize>,
    pub member: Vec<(usize, u8)>,
}

impl OrbitBase {
    pub fn build(glv: &GlvInstance3, base: &SubspaceBase) -> OrbitBase {
        let p = glv.inst.curve.field.p;
        let by_x: HashMap<u64, usize> = base
            .points
            .iter()
            .enumerate()
            .map(|(i, pt)| {
                assert!(glv.inst.curve.field.is_base(&pt.x));
                (pt.x.0[0], i)
            })
            .collect();
        let mut orbit_of_rep_x: HashMap<u64, usize> = HashMap::new();
        let mut reps = Vec::new();
        let mut member = vec![(usize::MAX, 0u8); base.len()];
        for (i, pt) in base.points.iter().enumerate() {
            let x = pt.x.0[0];
            let xs = [
                x,
                mm(x, glv.omega, p),
                mm(x, mm(glv.omega, glv.omega, p), p),
            ];
            let rep_x = *xs.iter().min().unwrap();
            let orbit = *orbit_of_rep_x.entry(rep_x).or_insert_with(|| {
                let rep = *by_x
                    .get(&rep_x)
                    .expect("ψ-orbit is not closed inside the subspace base");
                reps.push(rep);
                reps.len() - 1
            });
            // x = ω^k · rep_x.
            let k = (0..3u8)
                .find(|&k| mm(rep_x, glv.omega_pow(k), p) == x)
                .expect("x is in the orbit of its representative");
            member[i] = (orbit, k);
            debug_assert_eq!(glv.psi_k(&base.points[reps[orbit]], k), *pt);
        }
        OrbitBase { reps, member }
    }
    pub fn len(&self) -> usize {
        self.reps.len()
    }
    pub fn is_empty(&self) -> bool {
        self.reps.is_empty()
    }
    /// Fold decomposition terms `(base index, ±c)` into orbit columns
    /// with coefficients `c · λ^k`, merging terms that land in the same
    /// orbit; returns the sorted columns and the number of merges.
    pub fn fold_terms(
        &self,
        glv: &GlvInstance3,
        terms: &[(usize, i64)],
    ) -> (Vec<(usize, u64)>, u64) {
        let n = glv.inst.curve.n;
        let mut acc: HashMap<usize, u64> = HashMap::new();
        let mut merged = 0u64;
        for &(i, c) in terms {
            let (orbit, k) = self.member[i];
            let v = signed_mod(c, n);
            let v = mm(v, glv.lambda_pow(k), n);
            match acc.entry(orbit) {
                std::collections::hash_map::Entry::Occupied(mut e) => {
                    merged += 1;
                    *e.get_mut() = am(*e.get(), v, n);
                }
                std::collections::hash_map::Entry::Vacant(e) => {
                    e.insert(v);
                }
            }
        }
        let mut cols: Vec<(usize, u64)> = acc.into_iter().filter(|&(_, v)| v != 0).collect();
        cols.sort_unstable();
        (cols, merged)
    }
}

fn signed_mod(c: i64, n: u64) -> u64 {
    let v = c.unsigned_abs() % n;
    if c >= 0 {
        v
    } else {
        (n - v) % n
    }
}

// ── Relation stores with filtering and Wiedemann ────────────────────────

/// What one relation store (control or quotient) reports.
#[derive(Clone, Debug, Serialize, Default)]
pub struct StoreReport {
    pub label: String,
    /// Unknowns, including the `d` column.
    pub columns: usize,
    pub relations: u64,
    /// Decomposition terms merged because two points of one orbit
    /// appeared in the same decomposition.
    pub merged_entries: u64,
    /// Rows offered to the store that folded to `0 = 0` (a
    /// decomposition made of `ψ`-conjugates of the residual's own
    /// summands) and rows that were unit multiples of a stored row;
    /// neither is inserted.
    pub rows_zero: u64,
    pub rows_duplicate: u64,
    /// Non-zeros over every stored row; mean row weight; bytes of the
    /// store at 16 bytes per non-zero plus 32 per row.
    pub nnz: u64,
    pub avg_weight: f64,
    pub bytes: u64,
    pub residuals_at_solve: u64,
    pub solver_calls_at_solve: u64,
    pub group_ops_at_solve: u64,
    pub oracle_fp_muls_at_solve: u64,
    pub la_attempts: u64,
    /// Multiplications modulo `n` inside filtering and Wiedemann.
    pub la_ops: u64,
    pub core_rows: usize,
    pub core_nnz: u64,
    pub core_bytes: u64,
    /// Wall time of every linear-algebra attempt, summed.
    pub wiedemann_ms: f64,
    pub solved: bool,
    pub correct: Option<bool>,
    /// `group_ops + (oracle + precompute) / fp_per_add + la_ops · 9 / fp_per_add`.
    pub total_ops: f64,
    pub s: f64,
    /// `columns / decomposition rate`: residuals a square system needs.
    pub floor_residuals: f64,
    pub residual_ratio: f64,
}

struct Store {
    report: StoreReport,
    rels: Vec<SparseRel>,
    d_col: usize,
    next_attempt: usize,
    la_rng: StdRng,
    n: u64,
    row_keys: HashSet<Vec<(usize, u64)>>,
}

impl Store {
    fn new(label: &str, columns_without_d: usize, seed: u64, n: u64) -> Store {
        Store {
            report: StoreReport {
                label: label.to_string(),
                columns: columns_without_d + 1,
                ..StoreReport::default()
            },
            rels: Vec::new(),
            d_col: columns_without_d,
            next_attempt: columns_without_d + 1,
            la_rng: StdRng::seed_from_u64(seed ^ 0x1A_5EED),
            n,
            row_keys: HashSet::new(),
        }
    }
    /// Insert unless the row is zero or a unit multiple of a stored row.
    fn insert(&mut self, rel: SparseRel, merged: u64) -> bool {
        self.report.merged_entries += merged;
        if rel.cols.is_empty() {
            self.report.rows_zero += 1;
            return false;
        }
        if !self.row_keys.insert(row_key(&rel, self.n)) {
            self.report.rows_duplicate += 1;
            return false;
        }
        self.report.relations += 1;
        self.report.nnz += rel.cols.len() as u64;
        self.rels.push(rel);
        true
    }
    fn due(&self) -> bool {
        !self.report.solved && self.rels.len() >= self.next_attempt
    }
    /// Filter to a square core and run Wiedemann; `true` when the
    /// logarithm came out and the group confirms it.
    fn try_solve(&mut self, glv: &GlvInstance3) -> bool {
        let n = glv.inst.curve.n;
        let t0 = Instant::now();
        let mut ops = 0u64;
        let outcome = match square_core(&self.rels, self.d_col) {
            Some(core) => {
                self.report.la_attempts += 1;
                let unknowns = self.report.columns;
                let mut map = vec![usize::MAX; unknowns];
                let mut k = 0;
                for c in 0..unknowns {
                    if core.columns[c] {
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
                self.report.core_rows = sel.len();
                self.report.core_nnz = sel.iter().map(|r| r.cols.len() as u64).sum();
                self.report.core_bytes = 16 * self.report.core_nnz + 32 * sel.len() as u64;
                let x = wiedemann_u64(&sel, sel.len(), n, &mut self.la_rng, &mut ops);
                let dd = x.map(|x| x[map[self.d_col]]);
                if dd.is_some_and(|dd| glv.inst.curve.mul(&glv.inst.curve.g, dd) == glv.inst.q) {
                    self.report.solved = true;
                    self.report.correct = Some(dd == Some(glv.inst.d));
                    true
                } else {
                    self.next_attempt = self.rels.len() + (sel.len() / 20).max(1);
                    false
                }
            }
            None => {
                self.next_attempt = self.rels.len() + (self.report.columns / 20).max(1);
                false
            }
        };
        self.report.la_ops += ops;
        self.report.wiedemann_ms += t0.elapsed().as_secs_f64() * 1e3;
        outcome
    }
    fn finish(&mut self, fp_per_add: f64, precompute: u64, n: u64, rate: f64) {
        let r = &mut self.report;
        r.avg_weight = if r.relations == 0 {
            0.0
        } else {
            r.nnz as f64 / r.relations as f64
        };
        r.bytes = 16 * r.nnz + 32 * r.relations;
        r.total_ops = r.group_ops_at_solve as f64
            + (r.oracle_fp_muls_at_solve + precompute) as f64 / fp_per_add
            + r.la_ops as f64 * 9.0 / fp_per_add;
        r.s = r.total_ops / (n as f64).sqrt();
        r.floor_residuals = if rate > 0.0 {
            r.columns as f64 / rate
        } else {
            f64::INFINITY
        };
        r.residual_ratio = r.residuals_at_solve as f64 / r.floor_residuals;
    }
}

// ── Experiment 1: the factor-base quotient ──────────────────────────────

#[derive(Clone, Debug, Serialize)]
pub struct QuotientReport {
    pub p: u64,
    pub n: u64,
    pub bits: f64,
    pub omega: u64,
    pub lambda: u64,
    pub base: usize,
    pub orbits: usize,
    pub residuals: u64,
    pub solver_calls: u64,
    pub decompositions: u64,
    pub relations_verified: u64,
    pub decomposition_rate: f64,
    pub fp_muls_per_add: f64,
    pub precompute_fp_muls: u64,
    pub solve_stats: SolveStats,
    pub control: StoreReport,
    pub quotient: StoreReport,
    pub wall_ms: f64,
    pub peak_rss_bytes: u64,
    pub rho: RhoReport3,
    /// Plain rho on the same instance with `rho_runs` walk seeds: every
    /// run's `S`, their mean, minimum and maximum, so the reference
    /// carries its own spread instead of one draw.
    pub rho_runs: Vec<RhoReport3>,
    pub rho_s_mean: f64,
    pub rho_s_min: f64,
    pub rho_s_max: f64,
}

/// One residual stream, two relation stores: the ordinary `⟨−1⟩`-folded
/// base (control) and the `⟨−1, ψ⟩` orbit quotient.  Each store is
/// solved as soon as its own square core determines `d`, so the two
/// residual counts are what two stand-alone runs on the same stream
/// would need; the solver cost per residual is identical by
/// construction.
pub fn run_quotient_experiment(
    glv: &GlvInstance3,
    seed: u64,
    max_residuals: u64,
) -> QuotientReport {
    run_quotient_experiment_with_rho(glv, seed, max_residuals, 1)
}

/// As [`run_quotient_experiment`], with `rho_runs` independent rho
/// walks on the instance for the reference column.
pub fn run_quotient_experiment_with_rho(
    glv: &GlvInstance3,
    seed: u64,
    max_residuals: u64,
    rho_runs: u32,
) -> QuotientReport {
    let inst = &glv.inst;
    let curve = &inst.curve;
    let n = curve.n;
    let start = Instant::now();
    let mut rng = StdRng::seed_from_u64(seed ^ 0x9A0D);
    let base = SubspaceBase::build(inst, &mut rng);
    let orbit = OrbitBase::build(glv, &base);
    curve.reset();
    let fp_per_add = {
        let f = &curve.field;
        f.reset_muls();
        let mut acc = curve.g;
        for _ in 0..64 {
            acc = curve.add(&acc, &inst.q);
        }
        let m = f.muls() as f64 / 64.0;
        curve.reset();
        m
    };
    let pre = SymmetrisedS4::precompute(curve);
    let precompute = curve.field.muls();
    curve.reset();
    let mut stats = SolveStats::default();
    let mut control = Store::new("control ⟨−1⟩", base.len(), seed, n);
    let mut quotient = Store::new("quotient ⟨−1, ψ⟩", orbit.len(), seed, n);
    let mut stream = StdRng::seed_from_u64(seed ^ 0x57_AEA1);
    let mut residuals = 0u64;
    let mut decompositions = 0u64;
    let mut verified = 0u64;
    while residuals < max_residuals && !(control.report.solved && quotient.report.solved) {
        residuals += 1;
        let a = stream.gen_range(0..n);
        let b = stream.gen_range(1..n);
        let r = curve.add(&curve.mul(&curve.g, a), &curve.mul(&inst.q, b));
        let decs = subspace_triple_oracle_groebner(inst, &base, &pre, &r, &mut rng, &mut stats);
        for terms in decs {
            decompositions += 1;
            let mut acc = r;
            for &(i, c) in &terms {
                acc = curve.sub(&acc, &curve.mul_signed(&base.points[i], c));
            }
            if !acc.inf {
                continue;
            }
            verified += 1;
            if !control.report.solved {
                let mut cols: Vec<(usize, u64)> =
                    terms.iter().map(|&(i, c)| (i, signed_mod(c, n))).collect();
                cols.push((base.len(), (n - b) % n));
                cols.sort_unstable();
                control.insert(SparseRel { cols, rhs: a }, 0);
            }
            if !quotient.report.solved {
                let (mut cols, merged) = orbit.fold_terms(glv, &terms);
                cols.push((orbit.len(), (n - b) % n));
                cols.sort_unstable();
                quotient.insert(SparseRel { cols, rhs: a }, merged);
            }
        }
        for store in [&mut control, &mut quotient] {
            if store.due() && store.try_solve(glv) {
                store.report.residuals_at_solve = residuals;
                store.report.solver_calls_at_solve = residuals;
                store.report.group_ops_at_solve = curve.ops();
                store.report.oracle_fp_muls_at_solve = stats.fp_muls;
            }
        }
    }
    for store in [&mut control, &mut quotient] {
        if !store.report.solved {
            store.report.residuals_at_solve = residuals;
            store.report.solver_calls_at_solve = residuals;
            store.report.group_ops_at_solve = curve.ops();
            store.report.oracle_fp_muls_at_solve = stats.fp_muls;
        }
    }
    let rate = decompositions as f64 / residuals.max(1) as f64;
    control.finish(fp_per_add, precompute, n, rate);
    quotient.finish(fp_per_add, precompute, n, rate);
    let wall_ms = start.elapsed().as_secs_f64() * 1e3;
    let rho = run_rho3(inst, seed);
    let runs: Vec<RhoReport3> = (0..rho_runs.max(1))
        .map(|k| {
            run_rho3(
                inst,
                seed.wrapping_mul(1_000_003).wrapping_add(k as u64 + 1),
            )
        })
        .collect();
    let rho_s_mean = runs.iter().map(|r| r.s).sum::<f64>() / runs.len() as f64;
    let rho_s_min = runs.iter().map(|r| r.s).fold(f64::INFINITY, f64::min);
    let rho_s_max = runs.iter().map(|r| r.s).fold(0.0, f64::max);
    QuotientReport {
        p: curve.field.p,
        n,
        bits: (n as f64).log2(),
        omega: glv.omega,
        lambda: glv.lambda,
        base: base.len(),
        orbits: orbit.len(),
        residuals,
        solver_calls: residuals,
        decompositions,
        relations_verified: verified,
        decomposition_rate: rate,
        fp_muls_per_add: fp_per_add,
        precompute_fp_muls: precompute,
        solve_stats: stats,
        control: control.report,
        quotient: quotient.report,
        wall_ms,
        peak_rss_bytes: peak_rss_bytes(),
        rho,
        rho_runs: runs,
        rho_s_mean,
        rho_s_min,
        rho_s_max,
    }
}

// ── Experiment 2: canonical relation generation ─────────────────────────

#[derive(Clone, Debug, Serialize, Default)]
pub struct StreamReport {
    pub label: String,
    pub residuals: u64,
    /// Distinct `⟨−1, ψ⟩`-orbits among the residuals.
    pub distinct_orbits: u64,
    /// Residuals whose orbit had been seen: solver calls the canonical
    /// generator skips.
    pub duplicates: u64,
    pub solver_calls_naive: u64,
    pub solver_calls_canonical: u64,
    pub solver_calls_saved: u64,
    /// `C(T, 2) · 6 / n`: orbit collisions a uniform stream of `T`
    /// residuals is expected to contain.
    pub expected_duplicates_uniform: f64,
    /// Duplicates that were solved anyway and whose canonical rows
    /// were checked to be unit multiples of the stored ones.
    pub verified_duplicates: u64,
    pub verification_mismatches: u64,
    /// Rows produced by the canonical generator; how many folded to
    /// `0 = 0` (decompositions into `ψ`-conjugates of the residual's own
    /// summands) and how many were unit multiples of a row already in
    /// the store; neither is inserted.
    pub rows_produced: u64,
    pub rows_zero: u64,
    pub rows_duplicate: u64,
    pub merged_entries: u64,
    pub group_ops: u64,
    pub oracle_fp_muls: u64,
    pub wall_ms: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct CanonicalReport {
    pub p: u64,
    pub n: u64,
    pub bits: f64,
    pub base: usize,
    pub orbits: usize,
    pub random: StreamReport,
    pub pairs: StreamReport,
}

/// `r2 == u · r1` for some unit `u ∈ {±λ^k}`.
fn unit_multiple(glv: &GlvInstance3, r1: &SparseRel, r2: &SparseRel) -> bool {
    let n = glv.inst.curve.n;
    if r1.cols.len() != r2.cols.len() {
        return false;
    }
    for k in 0..3u8 {
        for sign in [1u64, n - 1] {
            let u = mm(sign, glv.lambda_pow(k), n);
            let same_cols = r1
                .cols
                .iter()
                .zip(&r2.cols)
                .all(|(&(c1, v1), &(c2, v2))| c1 == c2 && mm(u, v1, n) == v2);
            if same_cols && mm(u, r1.rhs, n) == r2.rhs {
                return true;
            }
        }
    }
    false
}

/// Key of a row up to scaling: the canonical row's columns plus its
/// right-hand side.
fn row_key(r: &SparseRel, n: u64) -> Vec<(usize, u64)> {
    let ck = canonical_row(r, n);
    let mut key = ck.cols;
    key.push((usize::MAX, ck.rhs));
    key
}

/// Canonical form of a quotient row: scaled so that its first non-zero
/// coefficient is `1` (unit multiples of one relation then coincide).
fn canonical_row(r: &SparseRel, n: u64) -> SparseRel {
    let Some(&(_, lead)) = r.cols.first() else {
        return r.clone();
    };
    let inv = inv_mod(lead, n);
    SparseRel {
        cols: r.cols.iter().map(|&(c, v)| (c, mm(v, inv, n))).collect(),
        rhs: mm(r.rhs, inv, n),
    }
}

struct CanonicalStream<'a> {
    glv: &'a GlvInstance3,
    base: &'a SubspaceBase,
    orbit: &'a OrbitBase,
    pre: &'a SymmetrisedS4,
    seen: HashMap<[u64; 6], Vec<SparseRel>>,
    row_keys: HashSet<Vec<(usize, u64)>>,
    report: StreamReport,
    verify_cap: u64,
}

impl<'a> CanonicalStream<'a> {
    /// Feed one residual `r` with its `(a, b)` tag (`a = b = 0` for a
    /// base-only relation); returns whether the solver ran.
    fn feed(
        &mut self,
        r: &Pt3,
        a: u64,
        b: u64,
        extra_terms: &[(usize, i64)],
        rng: &mut StdRng,
        stats: &mut SolveStats,
    ) -> bool {
        let glv = self.glv;
        let curve = &glv.inst.curve;
        let n = curve.n;
        self.report.residuals += 1;
        let key = glv.canonical_key(r);
        let (canon, u) = glv.canonical(r);
        let is_dup = self.seen.contains_key(&key);
        if is_dup {
            self.report.duplicates += 1;
            if self.report.verified_duplicates + self.report.verification_mismatches
                >= self.verify_cap
            {
                return false;
            }
        }
        // Solve the canonical representative; its decompositions are
        // those of `r` scaled by `u`.
        let decs =
            subspace_triple_oracle_groebner(&glv.inst, self.base, self.pre, &canon, rng, stats);
        let mut rows = Vec::new();
        for terms in decs {
            let mut acc = canon;
            for &(i, c) in &terms {
                acc = curve.sub(&acc, &curve.mul_signed(&self.base.points[i], c));
            }
            if !acc.inf {
                continue;
            }
            // canon = Σ c_i P_i  ⇒  r = u · canon = Σ u c_i P_i, and
            // r = aG + bQ + Σ extra.
            let (mut cols, merged) = self.orbit.fold_terms(glv, &terms);
            for e in cols.iter_mut() {
                e.1 = mm(e.1, u, n);
            }
            let (extra, merged_extra) = self.orbit.fold_terms(glv, extra_terms);
            let mut acc: HashMap<usize, u64> = cols.into_iter().collect();
            for (c, v) in extra {
                let e = acc.entry(c).or_insert(0);
                *e = sm(*e, v, n);
            }
            let mut cols: Vec<(usize, u64)> = acc.into_iter().filter(|&(_, v)| v != 0).collect();
            if b != 0 {
                cols.push((self.orbit.len(), (n - b) % n));
            }
            cols.sort_unstable();
            let rel = SparseRel { cols, rhs: a };
            self.report.merged_entries += merged + merged_extra;
            rows.push(rel);
        }
        if is_dup {
            let stored = &self.seen[&key];
            let ok = rows.len() == stored.len()
                && rows
                    .iter()
                    .all(|r2| stored.iter().any(|r1| unit_multiple(glv, r1, r2)));
            if ok {
                self.report.verified_duplicates += 1;
            } else {
                self.report.verification_mismatches += 1;
            }
            return true;
        }
        for rel in &rows {
            self.report.rows_produced += 1;
            if rel.cols.is_empty() {
                self.report.rows_zero += 1;
            } else if !self.row_keys.insert(row_key(rel, n)) {
                self.report.rows_duplicate += 1;
            }
        }
        self.seen.insert(key, rows);
        true
    }
    fn finish(mut self, calls: u64, group_ops: u64, muls: u64, t0: Instant) -> StreamReport {
        let n = self.glv.inst.curve.n as f64;
        let t = self.report.residuals as f64;
        self.report.distinct_orbits = self.seen.len() as u64;
        self.report.solver_calls_naive = self.report.residuals;
        self.report.solver_calls_canonical = self.report.distinct_orbits;
        self.report.solver_calls_saved = self.report.residuals - self.report.distinct_orbits;
        self.report.expected_duplicates_uniform = t * (t - 1.0) / 2.0 * 6.0 / n;
        let _ = calls;
        self.report.group_ops = group_ops;
        self.report.oracle_fp_muls = muls;
        self.report.wall_ms = t0.elapsed().as_secs_f64() * 1e3;
        self.report
    }
}

/// Canonical relation generation on two residual streams: `T` uniform
/// residuals `aG + bQ`, and the first `T'` base pairs `P_i + P_j`
/// (`i < j`, index order) as a sieve-style generator.  Every residual
/// is reduced to its `⟨−1, ψ⟩`-canonical form before the solver; a
/// residual whose canonical form was already solved is a saved solver
/// call, and up to `verify_cap` of them are solved anyway to check that
/// their rows are unit multiples of the stored ones.
pub fn run_canonical_experiment(
    glv: &GlvInstance3,
    seed: u64,
    random_residuals: u64,
    pair_residuals: u64,
    verify_cap: u64,
) -> CanonicalReport {
    let inst = &glv.inst;
    let curve = &inst.curve;
    let n = curve.n;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x9A0D);
    let base = SubspaceBase::build(inst, &mut rng);
    let orbit = OrbitBase::build(glv, &base);
    let pre = SymmetrisedS4::precompute(curve);
    let new_stream = |label: &str| CanonicalStream {
        glv,
        base: &base,
        orbit: &orbit,
        pre: &pre,
        seen: HashMap::new(),
        row_keys: HashSet::new(),
        report: StreamReport {
            label: label.to_string(),
            ..StreamReport::default()
        },
        verify_cap,
    };
    // Random stream.
    curve.reset();
    let mut stats = SolveStats::default();
    let t0 = Instant::now();
    let mut stream = StdRng::seed_from_u64(seed ^ 0x57_AEA1);
    let mut random = new_stream("uniform aG + bQ");
    let mut calls = 0u64;
    for _ in 0..random_residuals {
        let a = stream.gen_range(0..n);
        let b = stream.gen_range(1..n);
        let r = curve.add(&curve.mul(&curve.g, a), &curve.mul(&inst.q, b));
        if random.feed(&r, a, b, &[], &mut rng, &mut stats) {
            calls += 1;
        }
    }
    let random = random.finish(calls, curve.ops(), stats.fp_muls, t0);
    // Pair-enumeration stream.
    curve.reset();
    let mut stats = SolveStats::default();
    let t0 = Instant::now();
    let mut pairs = new_stream("base pairs P_i + P_j");
    let mut calls = 0u64;
    let mut count = 0u64;
    'outer: for i in 0..base.len() {
        for j in (i + 1)..base.len() {
            if count >= pair_residuals {
                break 'outer;
            }
            count += 1;
            let r = curve.add(&base.points[i], &base.points[j]);
            if pairs.feed(&r, 0, 0, &[(i, 1), (j, 1)], &mut rng, &mut stats) {
                calls += 1;
            }
        }
    }
    let pairs = pairs.finish(calls, curve.ops(), stats.fp_muls, t0);
    CanonicalReport {
        p: curve.field.p,
        n,
        bits: (n as f64).log2(),
        base: base.len(),
        orbits: orbit.len(),
        random,
        pairs,
    }
}

// ── Weighted polynomial systems and the Macaulay instrument ─────────────

/// A polynomial system over `F_p` with a `Z/3`-weight on every
/// variable: monomial `x^e` has weight `Σ wᵢ eᵢ (mod 3)`.
#[derive(Clone, Debug, Serialize)]
pub struct WeightedSystem {
    pub name: String,
    pub n_vars: usize,
    pub weights: Vec<u8>,
    /// Each polynomial as `(exponent, coefficient)` terms.
    pub polys: Vec<Vec<(Vec<u8>, u64)>>,
    /// Every polynomial has all its monomials in one weight class.
    pub homogeneous: bool,
    pub poly_weights: Vec<Option<u8>>,
    pub max_degree: u8,
    pub terms: usize,
}

fn mono_weight(e: &[u8], w: &[u8]) -> u8 {
    (e.iter()
        .zip(w)
        .map(|(&a, &b)| a as u32 * b as u32)
        .sum::<u32>()
        % 3) as u8
}

fn mono_degree(e: &[u8]) -> u8 {
    e.iter().sum()
}

impl WeightedSystem {
    pub fn new(name: &str, weights: Vec<u8>, polys: Vec<HashMap<Vec<u8>, u64>>) -> WeightedSystem {
        let n_vars = weights.len();
        let mut out_polys = Vec::new();
        let mut poly_weights = Vec::new();
        let mut homogeneous = true;
        let mut max_degree = 0u8;
        let mut terms = 0usize;
        for poly in polys {
            let mut v: Vec<(Vec<u8>, u64)> = poly.into_iter().filter(|(_, c)| *c != 0).collect();
            v.sort_by(|a, b| grevlex_cmp(&b.0, &a.0));
            let ws: HashSet<u8> = v.iter().map(|(e, _)| mono_weight(e, &weights)).collect();
            let pw = if ws.len() == 1 {
                ws.into_iter().next()
            } else {
                homogeneous = false;
                None
            };
            for (e, _) in &v {
                max_degree = max_degree.max(mono_degree(e));
            }
            terms += v.len();
            poly_weights.push(pw);
            out_polys.push(v);
        }
        WeightedSystem {
            name: name.to_string(),
            n_vars,
            weights,
            polys: out_polys,
            homogeneous,
            poly_weights,
            max_degree,
            terms,
        }
    }
    /// The system in the F4 engine's representation.
    pub fn f4_polys(&self, p: u64) -> Vec<f4_fp::Poly> {
        self.polys
            .iter()
            .map(|f| {
                let terms: Vec<(Vec<u32>, u64)> = f
                    .iter()
                    .map(|(e, c)| (e.iter().map(|&x| x as u32).collect(), *c))
                    .collect();
                f4_fp::normalise(&terms, p, F4Ordering::Grevlex)
            })
            .filter(|f| !f.is_empty())
            .collect()
    }
    /// Evaluate polynomial `i` at a point.
    pub fn eval(&self, i: usize, x: &[u64], p: u64) -> u64 {
        let mut acc = 0u64;
        for (e, c) in &self.polys[i] {
            let mut t = *c;
            for (k, &v) in x.iter().enumerate() {
                for _ in 0..e[k] {
                    t = mm(t, v, p);
                }
            }
            acc = am(acc, t, p);
        }
        acc
    }
    pub fn vanishes_at(&self, x: &[u64], p: u64) -> bool {
        (0..self.polys.len()).all(|i| self.eval(i, x, p) == 0)
    }
}

/// Graded reverse lexicographic comparison (`Greater` = larger monomial).
fn grevlex_cmp(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
    let (da, db) = (mono_degree(a), mono_degree(b));
    if da != db {
        return da.cmp(&db);
    }
    for i in (0..a.len()).rev() {
        if a[i] != b[i] {
            // The monomial with the *smaller* last differing exponent is larger.
            return b[i].cmp(&a[i]);
        }
    }
    std::cmp::Ordering::Equal
}

/// All monomials of total degree `≤ d` in `v` variables, largest first
/// in grevlex.
pub fn monomials_upto(v: usize, d: u8) -> Vec<Vec<u8>> {
    fn rec(v: usize, left: u8, cur: &mut Vec<u8>, out: &mut Vec<Vec<u8>>) {
        if cur.len() == v {
            out.push(cur.clone());
            return;
        }
        for e in 0..=left {
            cur.push(e);
            rec(v, left - e, cur, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    rec(v, d, &mut Vec::new(), &mut out);
    out.sort_by(|a, b| grevlex_cmp(b, a));
    out
}

/// Forward elimination to echelon form over `F_p`, sparse-aware: only
/// rows with a non-zero in the pivot column are touched and only the
/// non-zeros of the pivot row are multiplied; `muls` counts exactly
/// those.  Returns the pivot columns.
fn echelon(m: &mut [Vec<u64>], p: u64, muls: &mut u64) -> Vec<usize> {
    let rows = m.len();
    let cols = if rows > 0 { m[0].len() } else { 0 };
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..cols {
        if r >= rows {
            break;
        }
        let Some(pr) = (r..rows).find(|&i| m[i][c] != 0) else {
            continue;
        };
        m.swap(r, pr);
        let inv = inv_mod(m[r][c], p);
        let mut nz: Vec<usize> = Vec::new();
        for j in c..cols {
            if m[r][j] != 0 {
                m[r][j] = mm(m[r][j], inv, p);
                nz.push(j);
            }
        }
        *muls += nz.len() as u64;
        let pivot_row = m[r].clone();
        for i in (r + 1)..rows {
            if m[i][c] != 0 {
                let fct = m[i][c];
                for &j in &nz {
                    m[i][j] = sm(m[i][j], mm(fct, pivot_row[j], p), p);
                }
                *muls += nz.len() as u64;
            }
        }
        pivots.push(c);
        r += 1;
    }
    pivots
}

#[derive(Clone, Debug, Serialize, Default)]
pub struct BlockProfile {
    pub weight: u8,
    pub rows: usize,
    pub cols: usize,
    pub rank: usize,
    pub field_muls: u64,
}

/// The Macaulay matrix of a system at one degree, reduced.
#[derive(Clone, Debug, Serialize)]
pub struct MacaulayProfile {
    pub degree: u8,
    pub rows: usize,
    pub cols: usize,
    pub rank: usize,
    /// Non-pivot columns of degree `< degree`: the quotient dimension
    /// the truncation exhibits.
    pub standard: usize,
    /// Every product of a standard monomial by a variable is a pivot
    /// or standard: the multiplication matrices exist, the system is
    /// solved at this degree.
    pub closed: bool,
    pub field_muls: u64,
    /// Dense cells allocated (`rows × cols`, summed over blocks).
    pub dense_cells: u64,
    pub ms: f64,
    /// Rows whose support meets more than one weight class.
    pub coupled_rows: usize,
    pub block_aware: bool,
    pub blocks: Vec<BlockProfile>,
}

/// Build and reduce the Macaulay matrix of `sys` at `degree`.  With
/// `block_aware`, the columns are partitioned by weight and each block
/// is reduced on its own — only possible when the system is
/// weight-homogeneous (returns `None` otherwise, and `None` when the
/// dense matrix would exceed `cell_cap` cells).
pub fn macaulay_profile(
    sys: &WeightedSystem,
    p: u64,
    degree: u8,
    block_aware: bool,
    cell_cap: u64,
) -> Option<MacaulayProfile> {
    let t0 = Instant::now();
    let cols = monomials_upto(sys.n_vars, degree);
    let col_index: HashMap<Vec<u8>, usize> = cols
        .iter()
        .enumerate()
        .map(|(i, m)| (m.clone(), i))
        .collect();
    // Rows: every shift m · f_i with deg m ≤ degree − deg f_i.
    let mut row_terms: Vec<Vec<(usize, u64)>> = Vec::new();
    let mut row_weight: Vec<u8> = Vec::new();
    let mut coupled = 0usize;
    for (i, f) in sys.polys.iter().enumerate() {
        let df = f.iter().map(|(e, _)| mono_degree(e)).max().unwrap_or(0);
        if df > degree {
            continue;
        }
        for sh in monomials_upto(sys.n_vars, degree - df) {
            let mut terms: Vec<(usize, u64)> = Vec::with_capacity(f.len());
            let mut ws: HashSet<u8> = HashSet::new();
            for (e, c) in f {
                let m: Vec<u8> = e.iter().zip(&sh).map(|(a, b)| a + b).collect();
                ws.insert(mono_weight(&m, &sys.weights));
                terms.push((col_index[&m], *c));
            }
            if ws.len() > 1 {
                coupled += 1;
            }
            row_weight.push(
                sys.poly_weights[i]
                    .map(|w| (w as u32 + mono_weight(&sh, &sys.weights) as u32) as u8 % 3)
                    .unwrap_or(0),
            );
            row_terms.push(terms);
        }
    }
    let n_rows = row_terms.len();
    if block_aware && coupled > 0 {
        return None;
    }
    let mut field_muls = 0u64;
    let mut pivot_set: HashSet<usize> = HashSet::new();
    let mut blocks = Vec::new();
    let mut dense_cells = 0u64;
    if block_aware {
        for w in 0..3u8 {
            let bcols: Vec<usize> = (0..cols.len())
                .filter(|&c| mono_weight(&cols[c], &sys.weights) == w)
                .collect();
            let local: HashMap<usize, usize> =
                bcols.iter().enumerate().map(|(i, &c)| (c, i)).collect();
            let brows: Vec<usize> = (0..n_rows).filter(|&r| row_weight[r] == w).collect();
            let cells = brows.len() as u64 * bcols.len() as u64;
            if cells > cell_cap {
                return None;
            }
            dense_cells += cells;
            let mut m: Vec<Vec<u64>> = brows
                .iter()
                .map(|&r| {
                    let mut row = vec![0u64; bcols.len()];
                    for &(c, v) in &row_terms[r] {
                        let lc = local[&c];
                        row[lc] = am(row[lc], v, p);
                    }
                    row
                })
                .collect();
            let mut muls = 0u64;
            let piv = echelon(&mut m, p, &mut muls);
            field_muls += muls;
            blocks.push(BlockProfile {
                weight: w,
                rows: brows.len(),
                cols: bcols.len(),
                rank: piv.len(),
                field_muls: muls,
            });
            for c in piv {
                pivot_set.insert(bcols[c]);
            }
        }
    } else {
        let cells = n_rows as u64 * cols.len() as u64;
        if cells > cell_cap {
            return None;
        }
        dense_cells = cells;
        let mut m: Vec<Vec<u64>> = row_terms
            .iter()
            .map(|terms| {
                let mut row = vec![0u64; cols.len()];
                for &(c, v) in terms {
                    row[c] = am(row[c], v, p);
                }
                row
            })
            .collect();
        let piv = echelon(&mut m, p, &mut field_muls);
        pivot_set.extend(piv);
    }
    let standard: Vec<usize> = (0..cols.len())
        .filter(|&c| !pivot_set.contains(&c) && mono_degree(&cols[c]) < degree)
        .collect();
    let standard_set: HashSet<usize> = standard.iter().copied().collect();
    let closed = !standard.is_empty()
        && standard.iter().all(|&bc| {
            (0..sys.n_vars).all(|v| {
                let mut m = cols[bc].clone();
                m[v] += 1;
                match col_index.get(&m) {
                    Some(&c) => pivot_set.contains(&c) || standard_set.contains(&c),
                    None => false,
                }
            })
        });
    Some(MacaulayProfile {
        degree,
        rows: n_rows,
        cols: cols.len(),
        rank: pivot_set.len(),
        standard: standard.len(),
        closed,
        field_muls,
        dense_cells,
        ms: t0.elapsed().as_secs_f64() * 1e3,
        coupled_rows: coupled,
        block_aware,
        blocks,
    })
}

/// Macaulay profiles from `d_min` up, stopping at the first degree at
/// which the truncated quotient closes.
#[derive(Clone, Debug, Serialize)]
pub struct DegreeSweep {
    pub system: String,
    pub n_vars: usize,
    pub weights: Vec<u8>,
    pub homogeneous: bool,
    pub block_aware: bool,
    pub profiles: Vec<MacaulayProfile>,
    pub solve_degree: Option<u8>,
    pub quotient_dim: Option<usize>,
    pub peak_rows: usize,
    pub peak_cols: usize,
    pub peak_cells: u64,
    /// Field multiplications at the solving degree alone, and summed
    /// over every degree tried.
    pub solve_muls: u64,
    pub total_muls: u64,
    pub total_ms: f64,
    /// `Some(reason)` when the sweep stopped without closing.
    pub stopped: Option<String>,
}

pub fn degree_sweep(
    sys: &WeightedSystem,
    p: u64,
    d_min: u8,
    d_max: u8,
    block_aware: bool,
    cell_cap: u64,
) -> DegreeSweep {
    let mut sweep = DegreeSweep {
        system: sys.name.clone(),
        n_vars: sys.n_vars,
        weights: sys.weights.clone(),
        homogeneous: sys.homogeneous,
        block_aware,
        profiles: Vec::new(),
        solve_degree: None,
        quotient_dim: None,
        peak_rows: 0,
        peak_cols: 0,
        peak_cells: 0,
        solve_muls: 0,
        total_muls: 0,
        total_ms: 0.0,
        stopped: None,
    };
    if block_aware && !sys.homogeneous {
        sweep.stopped = Some("not weight-homogeneous: no block structure".into());
        return sweep;
    }
    for d in d_min..=d_max {
        let Some(prof) = macaulay_profile(sys, p, d, block_aware, cell_cap) else {
            sweep.stopped = Some(format!("degree {d}: matrix over the cell cap"));
            break;
        };
        sweep.peak_rows = sweep.peak_rows.max(prof.rows);
        sweep.peak_cols = sweep.peak_cols.max(prof.cols);
        sweep.peak_cells = sweep.peak_cells.max(prof.dense_cells);
        sweep.total_muls += prof.field_muls;
        sweep.total_ms += prof.ms;
        let closed = prof.closed;
        let dim = prof.standard;
        let muls = prof.field_muls;
        sweep.profiles.push(prof);
        if closed {
            sweep.solve_degree = Some(d);
            sweep.quotient_dim = Some(dim);
            sweep.solve_muls = muls;
            return sweep;
        }
    }
    if sweep.stopped.is_none() {
        sweep.stopped = Some(format!("no closure up to degree {d_max}"));
    }
    sweep
}

// ── System builders on a residual ───────────────────────────────────────

/// Split an `F_{p³}`-coefficient polynomial into its three `F_p`
/// components (Weil restriction of the equation).
fn weil_split(poly: &HashMap<Vec<u8>, E3>, p: u64) -> Vec<HashMap<Vec<u8>, u64>> {
    let mut out = vec![HashMap::new(), HashMap::new(), HashMap::new()];
    for (e, c) in poly {
        for k in 0..3 {
            if c.0[k] != 0 {
                let entry = out[k].entry(e.clone()).or_insert(0);
                *entry = am(*entry, c.0[k], p);
            }
        }
    }
    for o in out.iter_mut() {
        o.retain(|_, v| *v != 0);
    }
    out
}

/// The ordinary per-residual system: the three Weil components of the
/// symmetrised `S₄(e₁, e₂, e₃, x_R)`, weights `(1, 2, 0)`.
pub fn ordinary_system(pre: &SymmetrisedS4, f: &Fp3, x_r: &E3) -> WeightedSystem {
    let comps = pre.weil_restrict(f, x_r);
    let polys = comps
        .iter()
        .map(|c| c.iter().map(|(e, v)| (e.to_vec(), *v)).collect())
        .collect();
    WeightedSystem::new("ordinary", vec![1, 2, 0], polys)
}

/// Terms of the symmetrised `S₄` by `ψ`-weight `a + 2b + d (mod 3)` of
/// `e₁^a e₂^b e₃^c x₄^d`; the polynomial is `ψ`-invariant exactly when
/// every term has weight `0`.
pub fn s4_weight_histogram(pre: &SymmetrisedS4) -> [usize; 3] {
    let mut h = [0usize; 3];
    for e in pre.terms().keys() {
        h[(e[0] as usize + 2 * e[1] as usize + e[3] as usize) % 3] += 1;
    }
    h
}

/// The orbit system: `S₄(e₁, e₂, e₃, z·x_R)` with `z³ = 1`, `z^d`
/// reduced modulo `z³ − 1`, Weil-restricted, plus `z³ − 1`; unknowns
/// `(e₁, e₂, e₃, z)` with weights `(1, 2, 0, 1)`.  Its solutions are the
/// decompositions of `R`, `ψR` and `ψ²R` together (`z = 1, ω, ω²`).
pub fn orbit_system(pre: &SymmetrisedS4, f: &Fp3, x_r: &E3) -> WeightedSystem {
    let mut pow = [Fp3::ONE; 5];
    for i in 1..5 {
        pow[i] = f.mul(&pow[i - 1], x_r);
    }
    let mut h: HashMap<Vec<u8>, E3> = HashMap::new();
    for (e, c) in pre.terms() {
        let v = f.mul(c, &pow[e[3] as usize]);
        let key = vec![e[0], e[1], e[2], e[3] % 3];
        let entry = h.entry(key).or_insert(Fp3::ZERO);
        *entry = f.add(entry, &v);
    }
    h.retain(|_, v| *v != Fp3::ZERO);
    let mut polys = weil_split(&h, f.p);
    let mut z3 = HashMap::new();
    z3.insert(vec![0, 0, 0, 3], 1u64);
    z3.insert(vec![0, 0, 0, 0], f.p - 1);
    polys.push(z3);
    WeightedSystem::new("orbit (z-homogenised)", vec![1, 2, 0, 1], polys)
}

/// The `C₃`-invariant formulation.  With `X` the residual abscissa,
/// `u₁ = e₁/X`, `u₂ = e₂/X²`, `u₃ = e₃` are invariant under the
/// diagonal action (`X ↦ ωX` moves `R` along its orbit), and because
/// every term of `S₄` has weight `0`, `S₄(u₁X, u₂X², u₃, X)` is a
/// polynomial in `u` and the orbit invariant `c_R = X³` only.  The
/// `F_p`-rationality of `e₁, e₂` then pins `u₁ ∈ X⁻¹F_p`,
/// `u₂ ∈ X⁻²F_p`; substituting `u₁ = ẽ₁/X`, `u₂ = ẽ₂/X²` and
/// Weil-restricting gives the returned system in `(ẽ₁, ẽ₂, ẽ₃)`.
/// `orbit_invariant_only` reports whether the intermediate polynomial
/// used `x_R` through `c_R` alone.
pub fn invariant_system(pre: &SymmetrisedS4, f: &Fp3, x_r: &E3) -> (WeightedSystem, bool) {
    let c_r = f.mul(&f.sq(x_r), x_r);
    let inv_x = f.inv(x_r);
    let mut only_orbit_invariant = true;
    // H̃(u; c_R): coefficient of u₁^a u₂^b u₃^c is coef · c_R^{(a+2b+d)/3}.
    let mut h: HashMap<Vec<u8>, E3> = HashMap::new();
    for (e, c) in pre.terms() {
        let w = e[0] as u32 + 2 * e[1] as u32 + e[3] as u32;
        if !w.is_multiple_of(3) {
            only_orbit_invariant = false;
        }
        let v = f.mul(c, &f.pow(&c_r, (w / 3) as u128));
        let key = vec![e[0], e[1], e[2]];
        let entry = h.entry(key).or_insert(Fp3::ZERO);
        *entry = f.add(entry, &v);
    }
    // Rationality: u₁ = ẽ₁ X⁻¹, u₂ = ẽ₂ X⁻².
    let mut g: HashMap<Vec<u8>, E3> = HashMap::new();
    for (e, c) in h {
        let scale = f.pow(&inv_x, (e[0] as u32 + 2 * e[1] as u32) as u128);
        let v = f.mul(&c, &scale);
        if v != Fp3::ZERO {
            g.insert(e, v);
        }
    }
    let polys = weil_split(&g, f.p);
    (
        WeightedSystem::new("invariant (u-coordinates)", vec![1, 2, 0], polys),
        only_orbit_invariant,
    )
}

/// Minimal generators of the monoid of weight-`0` monomials up to
/// `max_degree` (a weight-`0` monomial that is not a product of two
/// non-constant weight-`0` monomials), and the number of weight-`0`
/// monomials in each degree.
pub fn invariant_generators(weights: &[u8], max_degree: u8) -> (Vec<Vec<u8>>, Vec<usize>) {
    let all = monomials_upto(weights.len(), max_degree);
    let inv: Vec<Vec<u8>> = all
        .into_iter()
        .filter(|m| mono_degree(m) > 0 && mono_weight(m, weights) == 0)
        .collect();
    let inv_set: HashSet<Vec<u8>> = inv.iter().cloned().collect();
    let mut hilbert = vec![0usize; max_degree as usize + 1];
    for m in &inv {
        hilbert[mono_degree(m) as usize] += 1;
    }
    let mut gens: Vec<Vec<u8>> = inv
        .iter()
        .filter(|m| {
            // Any proper non-constant divisor of weight 0 makes m a product.
            let mut divisor_found = false;
            let v = m.len();
            let mut d = vec![0u8; v];
            'enumerate: loop {
                // advance d as a mixed-radix counter bounded by m
                let mut i = 0;
                loop {
                    if i == v {
                        break 'enumerate;
                    }
                    if d[i] < m[i] {
                        d[i] += 1;
                        break;
                    }
                    d[i] = 0;
                    i += 1;
                }
                if d.as_slice() != m.as_slice() && mono_degree(&d) > 0 && inv_set.contains(&d) {
                    let rest: Vec<u8> = m.iter().zip(&d).map(|(a, b)| a - b).collect();
                    if inv_set.contains(&rest) {
                        divisor_found = true;
                        break;
                    }
                }
            }
            !divisor_found
        })
        .cloned()
        .collect();
    gens.sort();
    (gens, hilbert)
}

/// Render a monomial with variable names.
pub fn render_monomial(e: &[u8], names: &[&str]) -> String {
    let parts: Vec<String> = e
        .iter()
        .zip(names)
        .filter(|(&k, _)| k > 0)
        .map(|(&k, n)| {
            if k == 1 {
                n.to_string()
            } else {
                format!("{n}^{k}")
            }
        })
        .collect();
    if parts.is_empty() {
        "1".into()
    } else {
        parts.join("·")
    }
}

/// Multivariate polynomial over `F_p` unknowns with `F_{p³}`
/// coefficients (for the function-first construction).
#[derive(Clone, Debug)]
struct EPoly {
    terms: HashMap<Vec<u8>, E3>,
    n_vars: usize,
}

impl EPoly {
    fn zero(n_vars: usize) -> EPoly {
        EPoly {
            terms: HashMap::new(),
            n_vars,
        }
    }
    fn constant(n_vars: usize, c: E3) -> EPoly {
        let mut p = EPoly::zero(n_vars);
        if c != Fp3::ZERO {
            p.terms.insert(vec![0; n_vars], c);
        }
        p
    }
    fn var(n_vars: usize, i: usize, coef: E3) -> EPoly {
        let mut p = EPoly::zero(n_vars);
        let mut e = vec![0u8; n_vars];
        e[i] = 1;
        p.terms.insert(e, coef);
        p
    }
    fn add_term(&mut self, f: &Fp3, e: Vec<u8>, c: E3) {
        let entry = self.terms.entry(e.clone()).or_insert(Fp3::ZERO);
        *entry = f.add(entry, &c);
        if *entry == Fp3::ZERO {
            self.terms.remove(&e);
        }
    }
    fn add(&self, f: &Fp3, o: &EPoly) -> EPoly {
        let mut out = self.clone();
        for (e, c) in &o.terms {
            out.add_term(f, e.clone(), *c);
        }
        out
    }
    fn sub(&self, f: &Fp3, o: &EPoly) -> EPoly {
        let mut out = self.clone();
        for (e, c) in &o.terms {
            out.add_term(f, e.clone(), f.neg(c));
        }
        out
    }
    fn mul(&self, f: &Fp3, o: &EPoly) -> EPoly {
        let mut out = EPoly::zero(self.n_vars);
        for (e1, c1) in &self.terms {
            for (e2, c2) in &o.terms {
                let e: Vec<u8> = e1.iter().zip(e2).map(|(a, b)| a + b).collect();
                out.add_term(f, e, f.mul(c1, c2));
            }
        }
        out
    }
    fn scale(&self, f: &Fp3, k: &E3) -> EPoly {
        let mut out = EPoly::zero(self.n_vars);
        for (e, c) in &self.terms {
            out.add_term(f, e.clone(), f.mul(c, k));
        }
        out
    }
}

/// The prime-field **function-first** (Nagao `L(4O)`) system for the
/// residual `R = (x_R, y_R)`: a function `f = y + c₂x² + c₁x + c₀`
/// through `R` has three further zeros `P₁, P₂, P₃` with
/// `P₁ + P₂ + P₃ + R = O`; their abscissae are the roots of the cubic
/// `N(f)/(x − x_R)`, `N(f) = (c₂x² + c₁x + c₀)² − x³ − b`, so the
/// decomposition lies in the subspace base iff that cubic is
/// `c₂²·(x³ + t₁x² + t₂x + t₃)` with `t ∈ F_p`.  Unknowns: the six
/// `F_p`-coordinates of `c₁, c₂` and `t₁, t₂, t₃`; nine equations of
/// degree `≤ 3` (three Weil components of three coefficient
/// identities).  `c₀` is eliminated by `f(R) = 0`, which also makes the
/// division exact — the remainder is checked to vanish identically.
/// The roots' elementary symmetric functions are `e₁ = −t₁, e₂ = t₂,
/// e₃ = −t₃`.
pub fn function_first_system(curve: &Curve3, r: &Pt3) -> WeightedSystem {
    let f = &curve.field;
    let nv = 9;
    let t = E3([0, 1, 0]);
    let t2 = E3([0, 0, 1]);
    let one = Fp3::ONE;
    let c1 = EPoly::var(nv, 0, one)
        .add(f, &EPoly::var(nv, 1, t))
        .add(f, &EPoly::var(nv, 2, t2));
    let c2 = EPoly::var(nv, 3, one)
        .add(f, &EPoly::var(nv, 4, t))
        .add(f, &EPoly::var(nv, 5, t2));
    let t1 = EPoly::var(nv, 6, one);
    let t2p = EPoly::var(nv, 7, one);
    let t3 = EPoly::var(nv, 8, one);
    let xr = r.x;
    let xr2 = f.sq(&xr);
    // c₀ = −(y_R + c₂ x_R² + c₁ x_R).
    let c0 = EPoly::constant(nv, r.y)
        .add(f, &c2.scale(f, &xr2))
        .add(f, &c1.scale(f, &xr))
        .scale(f, &f.neg(&one));
    let two = f.from_base(2);
    // N = n₄x⁴ + n₃x³ + n₂x² + n₁x + n₀.
    let n4 = c2.mul(f, &c2);
    let n3 = c1
        .mul(f, &c2)
        .scale(f, &two)
        .sub(f, &EPoly::constant(nv, one));
    let n2 = c1.mul(f, &c1).add(f, &c0.mul(f, &c2).scale(f, &two));
    let n1 = c0.mul(f, &c1).scale(f, &two);
    let n0 = c0.mul(f, &c0).sub(f, &EPoly::constant(nv, curve.b));
    // Synthetic division by (x − x_R).
    let q3 = n4.clone();
    let q2 = n3.add(f, &q3.scale(f, &xr));
    let q1 = n2.add(f, &q2.scale(f, &xr));
    let q0 = n1.add(f, &q1.scale(f, &xr));
    let rem = n0.add(f, &q0.scale(f, &xr));
    assert!(
        rem.terms.is_empty(),
        "N(f) is divisible by x − x_R once f(R) = 0"
    );
    let eqs = [
        q2.sub(f, &q3.mul(f, &t1)),
        q1.sub(f, &q3.mul(f, &t2p)),
        q0.sub(f, &q3.mul(f, &t3)),
    ];
    let mut polys = Vec::new();
    for e in &eqs {
        polys.extend(weil_split(&e.terms, f.p));
    }
    // ψ-weights: c₁ ↦ ωc₁ (weight 1), c₂ ↦ ω²c₂ (weight 2), t₁ = −e₁
    // weight 1, t₂ = e₂ weight 2, t₃ weight 0.
    WeightedSystem::new(
        "function-first (Nagao L(4O))",
        vec![1, 1, 1, 2, 2, 2, 1, 2, 0],
        polys,
    )
}

// ── Experiment 3: the graded Macaulay matrix ────────────────────────────

#[derive(Clone, Debug, Serialize)]
pub struct GradedRow {
    pub residual_index: u64,
    pub x_r: [u64; 3],
    /// The harness's own solve on this residual: `F_p`-rational triples
    /// found, its quotient dimension, its `F_p` multiplications.
    pub harness_triples: Option<usize>,
    pub harness_quotient_dim: u64,
    pub harness_fp_muls: u64,
    pub harness_macaulay_rows: u64,
    /// Fraction of ordinary Macaulay rows (at its solving degree, or the
    /// last degree tried) meeting more than one weight class.
    pub ordinary_coupled_fraction: f64,
    pub ordinary: DegreeSweep,
    pub orbit_vanilla: DegreeSweep,
    pub orbit_block: DegreeSweep,
    /// Ratios at the solving degrees: orbit-block / ordinary.
    pub muls_ratio_block_vs_ordinary: Option<f64>,
    pub muls_ratio_block_vs_orbit_vanilla: Option<f64>,
    pub cells_ratio_block_vs_orbit_vanilla: Option<f64>,
    /// Degree-by-degree F4 (with substitution solving) on the ordinary
    /// system, on the orbit system with plain reduction, and on the
    /// orbit system with the matrix of every step split by weight.
    pub ordinary_f4: F4Summary,
    pub orbit_f4_vanilla: F4Summary,
    pub orbit_f4_block: F4Summary,
    /// F4 field multiplications, orbit-block over ordinary and over
    /// orbit-vanilla; wall time orbit-block over orbit-vanilla.
    pub f4_ops_ratio_block_vs_ordinary: Option<f64>,
    pub f4_ops_ratio_block_vs_vanilla: Option<f64>,
    pub f4_ms_ratio_block_vs_vanilla: Option<f64>,
    /// F4 found three times the ordinary solutions on the orbit system,
    /// with and without blocks.
    pub orbit_f4_triples_ordinary: Option<bool>,
}

#[derive(Clone, Debug, Serialize)]
pub struct GradedReport {
    pub p: u64,
    pub n: u64,
    pub bits: f64,
    pub s4_terms: usize,
    pub s4_weight_histogram: [usize; 3],
    pub s4_weight_homogeneous: bool,
    /// Residuals drawn and skipped because the harness found no
    /// `F_p`-rational decomposition (rows are on decomposable residuals).
    pub residuals_skipped: u64,
    pub f4_max_degree: u32,
    pub f4_budget_secs: f64,
    pub rows: Vec<GradedRow>,
}

/// Macaulay profiles of the ordinary and orbit systems on `residuals`
/// random residuals: vanilla elimination on both and block-aware
/// elimination on the orbit system.
#[allow(clippy::too_many_arguments)]
pub fn run_graded_experiment(
    glv: &GlvInstance3,
    seed: u64,
    residuals: u64,
    d_max_ordinary: u8,
    d_max_orbit: u8,
    cell_cap: u64,
    f4_max_degree: u32,
    f4_budget_secs: f64,
) -> GradedReport {
    let inst = &glv.inst;
    let curve = &inst.curve;
    let f = &curve.field;
    let n = curve.n;
    let p = f.p;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x9A0D);
    let pre = SymmetrisedS4::precompute(curve);
    let hist = s4_weight_histogram(&pre);
    let mut stream = StdRng::seed_from_u64(seed ^ 0x57_AEA1);
    let mut rows = Vec::new();
    let mut skipped = 0u64;
    let mut idx = 0u64;
    while (rows.len() as u64) < residuals {
        let a = stream.gen_range(0..n);
        let b = stream.gen_range(1..n);
        let r = curve.add(&curve.mul(&curve.g, a), &curve.mul(&inst.q, b));
        let mut stats = SolveStats::default();
        let triples = solve_s4_subspace(inst, &pre, &r.x, &mut rng, &mut stats);
        if triples.as_ref().is_none_or(|t| t.is_empty()) {
            skipped += 1;
            continue;
        }
        idx += 1;
        let ordinary = ordinary_system(&pre, f, &r.x);
        let orbit = orbit_system(&pre, f, &r.x);
        let (ord_f4, ord_sols) =
            f4_summary(&ordinary, p, f4_max_degree, f4_budget_secs, true, false);
        let (orb_v_f4, _) = f4_summary(&orbit, p, f4_max_degree, f4_budget_secs, true, false);
        let (orb_b_f4, _) = f4_summary(&orbit, p, f4_max_degree, f4_budget_secs, true, true);
        let _ = ord_sols;
        let f4_ratio = |a: &F4Summary, b: &F4Summary| -> Option<f64> {
            if !a.timed_out && !b.timed_out && b.field_ops > 0 {
                Some(a.field_ops as f64 / b.field_ops as f64)
            } else {
                None
            }
        };
        let f4_ms_ratio = if !orb_b_f4.timed_out && !orb_v_f4.timed_out && orb_v_f4.ms > 0.0 {
            Some(orb_b_f4.ms / orb_v_f4.ms)
        } else {
            None
        };
        let orbit_triples = match (orb_v_f4.solutions, orb_b_f4.solutions, ord_f4.solutions) {
            (Some(v), Some(bl), Some(e)) => Some(v == 3 * e && bl == 3 * e),
            _ => None,
        };
        let ord_sweep = degree_sweep(&ordinary, p, 4, d_max_ordinary, false, cell_cap);
        let orb_v = degree_sweep(&orbit, p, 4, d_max_orbit, false, cell_cap);
        let orb_b = degree_sweep(&orbit, p, 4, d_max_orbit, true, cell_cap);
        let coupled = ord_sweep
            .profiles
            .last()
            .map(|pr| pr.coupled_rows as f64 / pr.rows.max(1) as f64)
            .unwrap_or(0.0);
        let ratio = |num: &DegreeSweep, den: &DegreeSweep| -> Option<f64> {
            if num.solve_degree.is_some() && den.solve_degree.is_some() && den.solve_muls > 0 {
                Some(num.solve_muls as f64 / den.solve_muls as f64)
            } else {
                None
            }
        };
        let cells_ratio = if orb_b.solve_degree.is_some() && orb_v.solve_degree.is_some() {
            let cb = orb_b.profiles.last().map(|x| x.dense_cells).unwrap_or(0);
            let cv = orb_v.profiles.last().map(|x| x.dense_cells).unwrap_or(0);
            if cv > 0 {
                Some(cb as f64 / cv as f64)
            } else {
                None
            }
        } else {
            None
        };
        let row = GradedRow {
            residual_index: idx - 1,
            x_r: r.x.0,
            harness_triples: triples.map(|t| t.len()),
            harness_quotient_dim: stats.quotient_dim_total,
            harness_fp_muls: stats.fp_muls,
            harness_macaulay_rows: stats.macaulay_rows,
            ordinary_coupled_fraction: coupled,
            muls_ratio_block_vs_ordinary: ratio(&orb_b, &ord_sweep),
            muls_ratio_block_vs_orbit_vanilla: ratio(&orb_b, &orb_v),
            cells_ratio_block_vs_orbit_vanilla: cells_ratio,
            f4_ops_ratio_block_vs_ordinary: f4_ratio(&orb_b_f4, &ord_f4),
            f4_ops_ratio_block_vs_vanilla: f4_ratio(&orb_b_f4, &orb_v_f4),
            f4_ms_ratio_block_vs_vanilla: f4_ms_ratio,
            orbit_f4_triples_ordinary: orbit_triples,
            ordinary_f4: ord_f4,
            orbit_f4_vanilla: orb_v_f4,
            orbit_f4_block: orb_b_f4,
            ordinary: ord_sweep,
            orbit_vanilla: orb_v,
            orbit_block: orb_b,
        };
        rows.push(row);
    }
    GradedReport {
        p,
        n,
        bits: (n as f64).log2(),
        s4_terms: pre.terms().len(),
        s4_weight_histogram: hist,
        s4_weight_homogeneous: hist[1] == 0 && hist[2] == 0,
        residuals_skipped: skipped,
        f4_max_degree,
        f4_budget_secs,
        rows,
    }
}

// ── Experiment 4: the invariant formulation against the baselines ───────

#[derive(Clone, Debug, Serialize, Default)]
pub struct F4Summary {
    pub system: String,
    pub n_vars: usize,
    pub equations: usize,
    pub solving_degree: u32,
    pub degree_reached: u32,
    pub max_rows: usize,
    pub max_cols: usize,
    pub f4_runs: usize,
    pub field_ops: u64,
    pub ms: f64,
    pub timed_out: bool,
    pub inconsistent: bool,
    /// `F_p`-rational solutions found (`None` when undetermined).
    pub solutions: Option<usize>,
    /// Size of the degree-bounded basis (basis-only runs).
    pub basis_size: Option<usize>,
    pub pairs_above_bound: usize,
    /// Reduction partitioned by `ψ`-weight, and the number of steps
    /// (top-level run) that were reduced block by block.
    pub block_aware: bool,
    pub blocked_steps: usize,
    /// Basis-only runs: the highest step degree at which the run learned a
    /// new element, and the number of standard monomials of the basis
    /// (`None` when some variable has no pure power in the leading
    /// ideal, i.e. the basis is not zero-dimensional or is truncated).
    pub solving_degree_max: u32,
    pub staircase: Option<usize>,
}

/// Standard monomials of a Gröbner basis: monomials divisible by no
/// leading monomial, counted under the pure-power caps of the leading
/// ideal (`None` when a variable has no pure power).
pub fn staircase_size(basis: &[f4_fp::Poly], n_vars: usize) -> Option<usize> {
    let mut cap = vec![0u32; n_vars];
    for (k, c) in cap.iter_mut().enumerate() {
        *c = basis
            .iter()
            .map(|g| &g[0].0)
            .filter(|m| m.iter().enumerate().all(|(i, &e)| (i == k) == (e > 0)))
            .map(|m| m[k])
            .min()?;
    }
    fn walk(k: usize, mono: &mut Vec<u32>, cap: &[u32], lms: &[&Vec<u32>], count: &mut usize) {
        if lms
            .iter()
            .any(|l| l.iter().zip(mono.iter()).all(|(a, b)| a <= b))
        {
            return;
        }
        if k == cap.len() {
            *count += 1;
            return;
        }
        for e in 0..cap[k] {
            mono[k] = e;
            walk(k + 1, mono, cap, lms, count);
        }
        mono[k] = 0;
    }
    let lms: Vec<&Vec<u32>> = basis.iter().map(|g| &g[0].0).collect();
    let mut count = 0usize;
    let mut mono = vec![0u32; n_vars];
    walk(0, &mut mono, &cap, &lms, &mut count);
    Some(count)
}

/// Run the F4 engine on a system.  With `solve`, the substitution
/// solver extracts the `F_p`-rational solutions (cheap for three or
/// four unknowns); without it only the degree-bounded basis is
/// computed and the matrix profile reported.
fn f4_summary(
    sys: &WeightedSystem,
    p: u64,
    max_degree: u32,
    budget_secs: f64,
    solve: bool,
    block_aware: bool,
) -> (F4Summary, Option<Vec<Vec<u64>>>) {
    let polys = sys.f4_polys(p);
    let mut opts = F4Options::new(F4Ordering::Grevlex, max_degree)
        .with_budget(std::time::Duration::from_secs_f64(budget_secs));
    if block_aware {
        opts = opts.with_weights(sys.weights.clone());
    }
    let mut summary = F4Summary {
        system: if block_aware {
            format!("{} [block-aware]", sys.name)
        } else {
            sys.name.clone()
        },
        n_vars: sys.n_vars,
        equations: polys.len(),
        block_aware,
        ..F4Summary::default()
    };
    if solve {
        let rep = f4_fp::solve(&polys, sys.n_vars, p, &opts);
        let (solutions, sols, inconsistent) = match rep.verdict {
            Verdict::Solutions(s) => (Some(s.len()), Some(s), false),
            Verdict::Inconsistent => (Some(0), Some(Vec::new()), true),
            Verdict::Undetermined => (None, None, false),
        };
        summary.solving_degree = rep.solving_degree;
        summary.degree_reached = rep.degree_reached;
        summary.max_rows = rep.max_rows;
        summary.max_cols = rep.max_cols;
        summary.f4_runs = rep.f4_runs;
        summary.field_ops = rep.field_ops;
        summary.ms = rep.ms;
        summary.timed_out = rep.timed_out;
        summary.inconsistent = inconsistent;
        summary.solutions = solutions;
        summary.blocked_steps = rep.blocked_steps;
        (summary, sols)
    } else {
        let rep = f4_fp::f4(&polys, sys.n_vars, p, &opts);
        summary.solving_degree = rep.solving_degree;
        summary.degree_reached = rep.degree_reached;
        summary.max_rows = rep.max_rows;
        summary.max_cols = rep.max_cols;
        summary.f4_runs = 1;
        summary.field_ops = rep.field_ops;
        summary.ms = rep.ms;
        summary.timed_out = rep.timed_out;
        summary.inconsistent = rep.inconsistent;
        summary.basis_size = Some(rep.basis.len());
        summary.pairs_above_bound = rep.pairs_above_bound;
        summary.blocked_steps = rep.blocked_steps;
        summary.solving_degree_max = rep.solving_degree_max;
        summary.staircase = if rep.timed_out || rep.inconsistent {
            None
        } else {
            staircase_size(&rep.basis, sys.n_vars)
        };
        (summary, None)
    }
}

/// A witness for the function-first system from a harness triple:
/// signs `s_i` and the interpolating `f = y + c₂x² + c₁x + c₀` through
/// `s_iP_i` and `−R`; returns the nine-variable point when `f(R) = 0`
/// and the three coefficient identities hold, i.e. when the harness's
/// decomposition satisfies the function-first system.
pub fn function_first_witness(
    sys: &WeightedSystem,
    curve: &Curve3,
    base_y: impl Fn(u64) -> Option<E3>,
    r: &Pt3,
    triple: [u64; 3],
) -> Option<Vec<u64>> {
    let f = &curve.field;
    let p = f.p;
    let ys: Option<Vec<E3>> = triple.iter().map(|&x| base_y(x)).collect();
    let ys = ys?;
    for signs in 0..8u8 {
        let pts: Vec<Pt3> = (0..3)
            .map(|i| {
                let y = if (signs >> i) & 1 == 0 {
                    ys[i]
                } else {
                    f.neg(&ys[i])
                };
                Pt3::affine(f.from_base(triple[i]), y)
            })
            .collect();
        // Σ pts + R must be O for f to exist through all four.
        let mut acc = *r;
        for pt in &pts {
            acc = curve.add(&acc, pt);
        }
        if !acc.inf {
            continue;
        }
        // Interpolate c₂x² + c₁x + c₀ = −y through three points with
        // distinct x (Vandermonde); a repeated abscissa is skipped.
        if triple[0] == triple[1] || triple[1] == triple[2] || triple[0] == triple[2] {
            continue;
        }
        let x: Vec<E3> = pts.iter().map(|q| q.x).collect();
        let v: Vec<E3> = pts.iter().map(|q| f.neg(&q.y)).collect();
        // Lagrange: coefficients of Σ v_i Π_{j≠i} (X − x_j)/(x_i − x_j).
        let mut c = [Fp3::ZERO; 3];
        for i in 0..3 {
            let (j, k) = match i {
                0 => (1, 2),
                1 => (0, 2),
                _ => (0, 1),
            };
            let denom = f.inv(&f.mul(&f.sub(&x[i], &x[j]), &f.sub(&x[i], &x[k])));
            let s = f.mul(&v[i], &denom);
            // (X − x_j)(X − x_k) = X² − (x_j + x_k) X + x_j x_k
            c[2] = f.add(&c[2], &s);
            c[1] = f.sub(&c[1], &f.mul(&s, &f.add(&x[j], &x[k])));
            c[0] = f.add(&c[0], &f.mul(&s, &f.mul(&x[j], &x[k])));
        }
        // t from the cubic Π (X − x_i) = X³ + t₁X² + t₂X + t₃.
        let e1 = am(am(triple[0], triple[1], p), triple[2], p);
        let e2 = am(
            am(mm(triple[0], triple[1], p), mm(triple[0], triple[2], p), p),
            mm(triple[1], triple[2], p),
            p,
        );
        let e3 = mm(mm(triple[0], triple[1], p), triple[2], p);
        let point = vec![
            c[1].0[0],
            c[1].0[1],
            c[1].0[2],
            c[2].0[0],
            c[2].0[1],
            c[2].0[2],
            (p - e1) % p,
            e2,
            (p - e3) % p,
        ];
        if sys.vanishes_at(&point, p) {
            return Some(point);
        }
    }
    None
}

#[derive(Clone, Debug, Serialize)]
pub struct InvariantRow {
    pub residual_index: u64,
    pub x_r: [u64; 3],
    pub harness_triples: Option<Vec<[u64; 3]>>,
    pub harness_e_solutions: u64,
    pub harness_fp_muls: u64,
    /// The intermediate invariant polynomial depended on `R` through
    /// `c_R = x_R³` only.
    pub invariant_uses_orbit_invariant_only: bool,
    /// After the rationality substitution the invariant system equals
    /// the ordinary one term for term.
    pub invariant_identical_to_ordinary: bool,
    pub ordinary_macaulay: DegreeSweep,
    pub invariant_macaulay: DegreeSweep,
    pub function_first_macaulay: DegreeSweep,
    pub ordinary_f4: F4Summary,
    pub invariant_f4: F4Summary,
    pub orbit_f4: F4Summary,
    pub function_first_f4: F4Summary,
    /// Basis-only F4 (the complete grevlex basis, no substitution runs)
    /// on the same four formulations, with the staircase size: the cost
    /// of one Gröbner basis against another on the same variety.
    pub ordinary_f4_basis: F4Summary,
    pub invariant_f4_basis: F4Summary,
    pub orbit_f4_basis: F4Summary,
    pub function_first_f4_basis: F4Summary,
    /// F4 on the ordinary system found the harness's `e`-solutions.
    pub ordinary_f4_matches_harness: Option<bool>,
    /// The orbit system has three times the ordinary solutions.
    pub orbit_f4_triples_ordinary: Option<bool>,
    /// Every harness triple with distinct abscissae yields a point of
    /// the function-first system (interpolated `f`, then `t` from the
    /// cubic): the system is satisfied by the decompositions it is
    /// meant to find.
    pub function_first_witnessed: Option<bool>,
    pub function_first_witnesses: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct InvariantReport {
    pub p: u64,
    pub n: u64,
    pub bits: f64,
    /// Generators of the invariant monomials for weights `(1, 2, 0, 1)`
    /// on `(e₁, e₂, e₃, X)` and for the diagonal action `(1, 1, 1, 1)`
    /// on `(x₁, x₂, x₃, X)`, with the count of invariant monomials per
    /// degree.
    pub generators_symmetrised: Vec<String>,
    pub hilbert_symmetrised: Vec<usize>,
    pub generators_diagonal: Vec<String>,
    pub hilbert_diagonal: Vec<usize>,
    /// Residuals drawn and skipped because the harness found no
    /// `F_p`-rational decomposition.
    pub residuals_skipped: u64,
    pub f4_max_degree: u32,
    pub f4_budget_secs: f64,
    pub rows: Vec<InvariantRow>,
}

pub fn run_invariant_experiment(
    glv: &GlvInstance3,
    seed: u64,
    residuals: u64,
    f4_max_degree: u32,
    f4_budget_secs: f64,
    ff_d_max: u8,
    cell_cap: u64,
) -> InvariantReport {
    let inst = &glv.inst;
    let curve = &inst.curve;
    let f = &curve.field;
    let n = curve.n;
    let p = f.p;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x9A0D);
    let pre = SymmetrisedS4::precompute(curve);
    let (gs, hs) = invariant_generators(&[1, 2, 0, 1], 4);
    let (gd, hd) = invariant_generators(&[1, 1, 1, 1], 4);
    let sym_names = ["e1", "e2", "e3", "X"];
    let diag_names = ["x1", "x2", "x3", "X"];
    let mut stream = StdRng::seed_from_u64(seed ^ 0x57_AEA1);
    let mut rows = Vec::new();
    let mut skipped = 0u64;
    let mut idx = 0u64;
    while (rows.len() as u64) < residuals {
        let a = stream.gen_range(0..n);
        let b = stream.gen_range(1..n);
        let r = curve.add(&curve.mul(&curve.g, a), &curve.mul(&inst.q, b));
        let mut stats = SolveStats::default();
        let triples = solve_s4_subspace(inst, &pre, &r.x, &mut rng, &mut stats);
        if triples.as_ref().is_none_or(|t| t.is_empty()) {
            skipped += 1;
            continue;
        }
        idx += 1;
        let ordinary = ordinary_system(&pre, f, &r.x);
        let (invariant, orbit_only) = invariant_system(&pre, f, &r.x);
        let orbit = orbit_system(&pre, f, &r.x);
        let ff = function_first_system(curve, &r);
        let identical = ordinary.polys == invariant.polys;
        let ord_m = degree_sweep(&ordinary, p, 4, 13, false, cell_cap);
        let inv_m = degree_sweep(&invariant, p, 4, 13, false, cell_cap);
        let ff_m = degree_sweep(&ff, p, 2, ff_d_max, false, cell_cap);
        let (ord_f4, ord_sols) =
            f4_summary(&ordinary, p, f4_max_degree, f4_budget_secs, true, false);
        let (inv_f4, _) = f4_summary(&invariant, p, f4_max_degree, f4_budget_secs, true, false);
        let (orb_f4, _) = f4_summary(&orbit, p, f4_max_degree, f4_budget_secs, true, false);
        let (ff_f4, _) = f4_summary(&ff, p, f4_max_degree, f4_budget_secs, false, false);
        let (ord_b, _) = f4_summary(&ordinary, p, f4_max_degree, f4_budget_secs, false, false);
        let (inv_b, _) = f4_summary(&invariant, p, f4_max_degree, f4_budget_secs, false, false);
        let (orb_b, _) = f4_summary(&orbit, p, f4_max_degree, f4_budget_secs, false, false);
        let ordinary_matches = ord_sols.as_ref().map(|s| {
            s.iter().all(|x| ordinary.vanishes_at(x, p)) && s.len() as u64 == stats.e_solutions
        });
        let orbit_triples = match (orb_f4.solutions, ord_f4.solutions) {
            (Some(o), Some(e)) => Some(o == 3 * e),
            _ => None,
        };
        let base_y = |x: u64| -> Option<E3> {
            let xe = f.from_base(x);
            let rhs = curve.rhs(&xe);
            let mut r2 = StdRng::seed_from_u64(x ^ 0xB0B);
            f.sqrt(&rhs, &mut r2)
        };
        let (ff_witnessed, ff_witnesses) = match &triples {
            Some(tr) => {
                let distinct: Vec<[u64; 3]> = tr
                    .iter()
                    .copied()
                    .filter(|t| t[0] != t[1] && t[1] != t[2] && t[0] != t[2])
                    .collect();
                let found = distinct
                    .iter()
                    .filter(|&&t| function_first_witness(&ff, curve, base_y, &r, t).is_some())
                    .count();
                (Some(found == distinct.len()), found)
            }
            None => (None, 0),
        };
        rows.push(InvariantRow {
            residual_index: idx - 1,
            x_r: r.x.0,
            harness_triples: triples,
            harness_e_solutions: stats.e_solutions,
            harness_fp_muls: stats.fp_muls,
            invariant_uses_orbit_invariant_only: orbit_only,
            invariant_identical_to_ordinary: identical,
            ordinary_macaulay: ord_m,
            invariant_macaulay: inv_m,
            function_first_macaulay: ff_m,
            ordinary_f4: ord_f4,
            invariant_f4: inv_f4,
            orbit_f4: orb_f4,
            function_first_f4_basis: ff_f4.clone(),
            function_first_f4: ff_f4,
            ordinary_f4_basis: ord_b,
            invariant_f4_basis: inv_b,
            orbit_f4_basis: orb_b,
            ordinary_f4_matches_harness: ordinary_matches,
            orbit_f4_triples_ordinary: orbit_triples,
            function_first_witnessed: ff_witnessed,
            function_first_witnesses: ff_witnesses,
        });
    }
    InvariantReport {
        p,
        n,
        bits: (n as f64).log2(),
        generators_symmetrised: gs.iter().map(|g| render_monomial(g, &sym_names)).collect(),
        hilbert_symmetrised: hs,
        generators_diagonal: gd.iter().map(|g| render_monomial(g, &diag_names)).collect(),
        hilbert_diagonal: hd,
        residuals_skipped: skipped,
        f4_max_degree,
        f4_budget_secs,
        rows,
    }
}

// ── The Veronese (unsymmetrised diagonal-invariant) formulation ─────────

/// The cubic monomials in `(x₁, x₂, x₃, z)`: the twenty generators of
/// the invariant ring of the diagonal `C₃` action with weights
/// `(1, 1, 1, 1)`, in grevlex order (largest first).
pub fn veronese_generators() -> Vec<[u8; 4]> {
    let mut out: Vec<[u8; 4]> = Vec::new();
    for a in 0..=3u8 {
        for b in 0..=3 - a {
            for c in 0..=3 - a - b {
                out.push([a, b, c, 3 - a - b - c]);
            }
        }
    }
    out.sort_by(|x, y| grevlex_cmp(y, x));
    out
}

/// Write a monomial of degree `3k` in `(x₁, x₂, x₃, z)` as a product of
/// `k` cubic generators: peel the grevlex-largest generator dividing what
/// is left, `k` times.  Returns the exponent vector over the twenty
/// generators.
fn veronese_factor(mono: [u8; 4], gens: &[[u8; 4]]) -> Vec<u8> {
    let mut rest = mono;
    let mut e = vec![0u8; gens.len()];
    while rest.iter().map(|&v| v as u32).sum::<u32>() > 0 {
        let (i, g) = gens
            .iter()
            .enumerate()
            .find(|(_, g)| g.iter().zip(&rest).all(|(a, b)| a <= b))
            .expect("a monomial of degree 3k is divisible by a cubic");
        e[i] += 1;
        for (r, &gv) in rest.iter_mut().zip(g) {
            *r -= gv;
        }
    }
    e
}

/// The orbit PDP in the invariants of the diagonal `C₃` action on
/// `(x₁, x₂, x₃, z)`: unknowns `m_α = x^α`, one per cubic monomial
/// (twenty, all `F_p`-valued since `z ∈ F_p` with `z³ = 1`); equations
/// the three Weil components of `S₄(x₁, x₂, x₃, z·x_R)`, each of weight
/// `0` and hence a polynomial in the `m_α` of degree `≤ 4`, together with
/// `m_{z³} = 1` and the toric relations `m_α m_β = m_γ m_δ` of the cubic
/// Veronese (`α + β = γ + δ`).  Its points are the `(x, z)` of the orbit
/// system with `x` ordered and the three points of each `ψ`-orbit
/// identified (every cubic monomial is `ζ`-invariant), so its staircase
/// is `6 × 64 = 384` against the orbit system's `192`.
pub fn veronese_system(curve: &Curve3, x_r: &E3) -> WeightedSystem {
    let f = &curve.field;
    let p = f.p;
    let gens = veronese_generators();
    let n_vars = gens.len();
    // S₄(x₁, x₂, x₃, z·x_R) with z^d reduced modulo z³ − 1.
    let mut pow = [Fp3::ONE; 5];
    for i in 1..5 {
        pow[i] = f.mul(&pow[i - 1], x_r);
    }
    let mut h: HashMap<[u8; 4], E3> = HashMap::new();
    for (e, c) in s4_terms(curve) {
        let v = f.mul(&c, &pow[e[3] as usize]);
        let key = [e[0], e[1], e[2], e[3] % 3];
        let entry = h.entry(key).or_insert(Fp3::ZERO);
        *entry = f.add(entry, &v);
    }
    h.retain(|_, v| *v != Fp3::ZERO);
    let mut in_m: HashMap<Vec<u8>, E3> = HashMap::new();
    for (mono, c) in h {
        assert!(
            mono.iter().map(|&v| v as u32).sum::<u32>() % 3 == 0,
            "S₄ on a j = 0 curve is ψ-invariant"
        );
        let key = veronese_factor(mono, &gens);
        let entry = in_m.entry(key).or_insert(Fp3::ZERO);
        *entry = f.add(entry, &c);
    }
    let mut polys = weil_split(&in_m, p);
    // m_{z³} = 1.
    let z3 = gens
        .iter()
        .position(|g| *g == [0, 0, 0, 3])
        .expect("z³ is a generator");
    let mut eq = HashMap::new();
    let mut e = vec![0u8; n_vars];
    e[z3] = 1;
    eq.insert(e, 1u64);
    eq.insert(vec![0u8; n_vars], p - 1);
    polys.push(eq);
    // Toric relations: for every pair of generators, the product's
    // canonical factorisation must agree with the pair.
    let mut seen: HashSet<Vec<u8>> = HashSet::new();
    for i in 0..n_vars {
        for j in i..n_vars {
            let prod = [
                gens[i][0] + gens[j][0],
                gens[i][1] + gens[j][1],
                gens[i][2] + gens[j][2],
                gens[i][3] + gens[j][3],
            ];
            let canon = veronese_factor(prod, &gens);
            let mut pair = vec![0u8; n_vars];
            pair[i] += 1;
            pair[j] += 1;
            if pair == canon {
                continue;
            }
            let mut rel: Vec<u8> = pair.clone();
            rel.extend(&canon);
            if !seen.insert(rel) {
                continue;
            }
            let mut eq = HashMap::new();
            eq.insert(pair, 1u64);
            eq.insert(canon, p - 1);
            polys.push(eq);
        }
    }
    WeightedSystem::new("veronese (20 cubic invariants)", vec![0; n_vars], polys)
}

#[derive(Clone, Debug, Serialize)]
pub struct VeroneseRow {
    pub residual_index: u64,
    pub x_r: [u64; 3],
    pub harness_triples: Option<Vec<[u64; 3]>>,
    pub harness_fp_muls: u64,
    pub equations: usize,
    pub toric_relations: usize,
    /// Every harness triple, with the signs the harness confirms, gives a
    /// point of the Veronese system (`m_α = x^α` at `z = 1`).
    pub witnessed: Option<bool>,
    pub witnesses: usize,
    pub ordinary_f4_basis: F4Summary,
    pub orbit_f4_basis: F4Summary,
    pub veronese_f4_basis: F4Summary,
    pub mults_vs_ordinary: Option<f64>,
    pub mults_vs_orbit: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct VeroneseReport {
    pub p: u64,
    pub n: u64,
    pub bits: f64,
    pub generators: Vec<String>,
    pub residuals_skipped: u64,
    pub f4_max_degree: u32,
    pub f4_budget_secs: f64,
    pub rows: Vec<VeroneseRow>,
}

/// The Veronese formulation against the ordinary and orbit bases on the
/// same decomposable residuals (basis-only F4 with staircase sizes).
pub fn run_veronese_experiment(
    glv: &GlvInstance3,
    seed: u64,
    residuals: u64,
    f4_max_degree: u32,
    f4_budget_secs: f64,
) -> VeroneseReport {
    let inst = &glv.inst;
    let curve = &inst.curve;
    let f = &curve.field;
    let n = curve.n;
    let p = f.p;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x9A0D);
    let pre = SymmetrisedS4::precompute(curve);
    let gens = veronese_generators();
    let names = ["x1", "x2", "x3", "z"];
    let mut stream = StdRng::seed_from_u64(seed ^ 0x57_AEA1);
    let mut rows = Vec::new();
    let mut skipped = 0u64;
    let mut idx = 0u64;
    while (rows.len() as u64) < residuals {
        let a = stream.gen_range(0..n);
        let b = stream.gen_range(1..n);
        let r = curve.add(&curve.mul(&curve.g, a), &curve.mul(&inst.q, b));
        let mut stats = SolveStats::default();
        let triples = solve_s4_subspace(inst, &pre, &r.x, &mut rng, &mut stats);
        if triples.as_ref().is_none_or(|t| t.is_empty()) {
            skipped += 1;
            continue;
        }
        idx += 1;
        let ordinary = ordinary_system(&pre, f, &r.x);
        let orbit = orbit_system(&pre, f, &r.x);
        let veronese = veronese_system(curve, &r.x);
        let toric = veronese.polys.len() - 4;
        // Witnesses: every harness triple, in each of its six orders, at
        // z = 1, evaluated in the generators.
        let (witnessed, witnesses) = match &triples {
            Some(tr) => {
                let mut ok = true;
                let mut count = 0usize;
                for t in tr {
                    let pt: Vec<u64> = gens
                        .iter()
                        .map(|g| {
                            let mut v = 1u64;
                            for (k, &e) in g.iter().enumerate() {
                                let base = if k < 3 { t[k] } else { 1 };
                                for _ in 0..e {
                                    v = mm(v, base, p);
                                }
                            }
                            v
                        })
                        .collect();
                    if veronese.vanishes_at(&pt, p) {
                        count += 1;
                    } else {
                        ok = false;
                    }
                }
                (Some(ok), count)
            }
            None => (None, 0),
        };
        let (ord_b, _) = f4_summary(&ordinary, p, f4_max_degree, f4_budget_secs, false, false);
        let (orb_b, _) = f4_summary(&orbit, p, f4_max_degree, f4_budget_secs, false, false);
        let (ver_b, _) = f4_summary(&veronese, p, f4_max_degree, f4_budget_secs, false, false);
        let ratio = |a: &F4Summary, b: &F4Summary| -> Option<f64> {
            if !a.timed_out && !b.timed_out && b.field_ops > 0 {
                Some(a.field_ops as f64 / b.field_ops as f64)
            } else {
                None
            }
        };
        rows.push(VeroneseRow {
            residual_index: idx - 1,
            x_r: r.x.0,
            harness_triples: triples,
            harness_fp_muls: stats.fp_muls,
            equations: veronese.polys.len(),
            toric_relations: toric,
            witnessed,
            witnesses,
            mults_vs_ordinary: ratio(&ver_b, &ord_b),
            mults_vs_orbit: ratio(&ver_b, &orb_b),
            ordinary_f4_basis: ord_b,
            orbit_f4_basis: orb_b,
            veronese_f4_basis: ver_b,
        });
    }
    VeroneseReport {
        p,
        n,
        bits: (n as f64).log2(),
        generators: gens.iter().map(|g| render_monomial(g, &names)).collect(),
        residuals_skipped: skipped,
        f4_max_degree,
        f4_budget_secs,
        rows,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny() -> GlvInstance3 {
        generate_j0_instance3(73, 3)
    }

    #[test]
    fn j0_instance_has_order_three_automorphism() {
        let glv = tiny();
        let c = &glv.inst.curve;
        assert!(is_prime_u64(c.n) && c.n % 3 == 1);
        assert_eq!(c.a, Fp3::ZERO);
        assert_eq!(pow_mod(glv.omega, 3, c.field.p), 1);
        assert_ne!(glv.omega, 1);
        let l = glv.lambda;
        assert_eq!((mm(l, l, c.n) + l + 1) % c.n, 0);
        for k in 1..20u64 {
            let pt = c.mul(&c.g, k);
            assert!(c.is_on_curve(&glv.psi(&pt)));
            assert_eq!(glv.psi(&pt), c.mul(&pt, l));
            assert_eq!(glv.psi_k(&pt, 3), pt);
        }
        assert_eq!(c.mul(&c.g, glv.inst.d), glv.inst.q);
    }

    #[test]
    fn canonical_form_is_orbit_invariant_and_unit_correct() {
        let glv = tiny();
        let c = &glv.inst.curve;
        for k in 1..30u64 {
            let pt = c.mul(&c.g, k);
            let (canon, u) = glv.canonical(&pt);
            assert_eq!(c.mul(&canon, u), pt);
            let key = glv.canonical_key(&pt);
            for j in 0..3u8 {
                let s = glv.psi_k(&pt, j);
                assert_eq!(glv.canonical_key(&s), key);
                assert_eq!(glv.canonical_key(&c.neg(&s)), key);
            }
        }
    }

    #[test]
    fn orbit_base_is_a_three_to_one_quotient() {
        let glv = tiny();
        let mut rng = StdRng::seed_from_u64(1);
        let base = SubspaceBase::build(&glv.inst, &mut rng);
        let orbit = OrbitBase::build(&glv, &base);
        assert_eq!(base.len(), 3 * orbit.len());
        let c = &glv.inst.curve;
        for (i, &(o, k)) in orbit.member.iter().enumerate() {
            let rep = base.points[orbit.reps[o]];
            assert_eq!(glv.psi_k(&rep, k), base.points[i]);
            assert_eq!(c.mul(&rep, glv.lambda_pow(k)), base.points[i]);
        }
        // Folded terms reproduce the point: Σ c λ^k Rep = Σ c P.
        let terms = vec![(0usize, 1i64), (1, -1), (2, 1)];
        let (cols, _) = orbit.fold_terms(&glv, &terms);
        let mut lhs = Pt3::INFINITY;
        for &(o, v) in &cols {
            lhs = c.add(&lhs, &c.mul(&base.points[orbit.reps[o]], v));
        }
        let mut rhs = Pt3::INFINITY;
        for &(i, s) in &terms {
            rhs = c.add(&rhs, &c.mul_signed(&base.points[i], s));
        }
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn s4_is_psi_invariant_on_j0() {
        let glv = tiny();
        let pre = SymmetrisedS4::precompute(&glv.inst.curve);
        let h = s4_weight_histogram(&pre);
        assert_eq!(h[1] + h[2], 0, "histogram {h:?}");
        let f = &glv.inst.curve.field;
        let mut rng = StdRng::seed_from_u64(7);
        let x = E3([
            rng.gen_range(0..f.p),
            rng.gen_range(0..f.p),
            rng.gen_range(0..f.p),
        ]);
        let orbit = orbit_system(&pre, f, &x);
        assert!(orbit.homogeneous);
        let (inv, only) = invariant_system(&pre, f, &x);
        assert!(only);
        let ord = ordinary_system(&pre, f, &x);
        assert_eq!(ord.polys, inv.polys);
        assert!(!ord.homogeneous);
    }

    #[test]
    fn block_and_vanilla_elimination_agree_on_the_orbit_system() {
        let glv = tiny();
        let f = &glv.inst.curve.field;
        let pre = SymmetrisedS4::precompute(&glv.inst.curve);
        let x = E3([5, 9, 2]);
        let orbit = orbit_system(&pre, f, &x);
        let v = macaulay_profile(&orbit, f.p, 7, false, 1 << 28).unwrap();
        let b = macaulay_profile(&orbit, f.p, 7, true, 1 << 28).unwrap();
        assert_eq!(v.rank, b.rank);
        assert_eq!(v.standard, b.standard);
        assert_eq!(v.rows, b.rows);
        assert_eq!(v.cols, b.blocks.iter().map(|x| x.cols).sum::<usize>());
        assert!(b.dense_cells < v.dense_cells);
        assert_eq!(v.coupled_rows, 0);
    }

    #[test]
    fn block_aware_f4_solves_the_orbit_system_like_plain_f4() {
        let glv = tiny();
        let f = &glv.inst.curve.field;
        let pre = SymmetrisedS4::precompute(&glv.inst.curve);
        let x = E3([5, 9, 2]);
        let orbit = orbit_system(&pre, f, &x);
        let ordinary = ordinary_system(&pre, f, &x);
        let (o, _) = f4_summary(&ordinary, f.p, 14, 60.0, true, false);
        let (v, _) = f4_summary(&orbit, f.p, 14, 60.0, true, false);
        let (b, _) = f4_summary(&orbit, f.p, 14, 60.0, true, true);
        assert!(!o.timed_out && !v.timed_out && !b.timed_out);
        assert_eq!(v.solutions, b.solutions);
        assert_eq!(v.solutions, o.solutions.map(|k| 3 * k));
        assert!(b.blocked_steps > 0);
        assert_eq!(v.solving_degree, b.solving_degree);
    }

    #[test]
    fn veronese_system_is_witnessed_by_harness_triples() {
        let glv = tiny();
        let c = &glv.inst.curve;
        let f = &c.field;
        let pre = SymmetrisedS4::precompute(c);
        let mut rng = StdRng::seed_from_u64(11);
        let mut stats = SolveStats::default();
        for k in 2..600u64 {
            let r = c.mul(&c.g, k);
            let Some(triples) = solve_s4_subspace(&glv.inst, &pre, &r.x, &mut rng, &mut stats)
            else {
                continue;
            };
            if triples.is_empty() {
                continue;
            }
            let ver = veronese_system(c, &r.x);
            assert_eq!(ver.n_vars, 20);
            let gens = veronese_generators();
            for t in &triples {
                let pt: Vec<u64> = gens
                    .iter()
                    .map(|g| {
                        let mut v = 1u64;
                        for (i, &e) in g.iter().enumerate() {
                            let base = if i < 3 { t[i] } else { 1 };
                            for _ in 0..e {
                                v = mm(v, base, f.p);
                            }
                        }
                        v
                    })
                    .collect();
                assert!(ver.vanishes_at(&pt, f.p));
            }
            // 20 generators, 3 Weil components, z³ = 1, and the 126
            // toric quadrics of the cubic Veronese.
            assert_eq!(ver.polys.len(), 130);
            return;
        }
        panic!("no decomposable multiple of g found");
    }

    #[test]
    fn invariant_generators_match_the_hand_count() {
        let (gs, _) = invariant_generators(&[1, 2, 0, 1], 4);
        let names = ["e1", "e2", "e3", "X"];
        let rendered: Vec<String> = gs.iter().map(|g| render_monomial(g, &names)).collect();
        for want in [
            "e3", "e1·e2", "e2·X", "e1^3", "X^3", "e1^2·X", "e1·X^2", "e2^3",
        ] {
            assert!(rendered.contains(&want.to_string()), "{rendered:?}");
        }
        assert_eq!(gs.len(), 8);
        let (gd, hd) = invariant_generators(&[1, 1, 1, 1], 4);
        assert_eq!(gd.len(), 20);
        assert_eq!(hd[3], 20);
    }

    #[test]
    fn function_first_system_is_witnessed_by_harness_decompositions() {
        let glv = tiny();
        let c = &glv.inst.curve;
        let f = &c.field;
        let pre = SymmetrisedS4::precompute(c);
        let mut rng = StdRng::seed_from_u64(11);
        let mut stats = SolveStats::default();
        let base_y = |x: u64| -> Option<E3> {
            let mut r2 = StdRng::seed_from_u64(x ^ 0xB0B);
            f.sqrt(&c.rhs(&f.from_base(x)), &mut r2)
        };
        let mut witnessed = 0;
        for k in 2..600u64 {
            let r = c.mul(&c.g, k);
            let Some(triples) = solve_s4_subspace(&glv.inst, &pre, &r.x, &mut rng, &mut stats)
            else {
                continue;
            };
            let ff = function_first_system(c, &r);
            assert_eq!(ff.n_vars, 9);
            assert_eq!(ff.polys.len(), 9);
            for t in triples {
                if t[0] != t[1] && t[1] != t[2] && t[0] != t[2] {
                    assert!(
                        function_first_witness(&ff, c, base_y, &r, t).is_some(),
                        "harness triple {t:?} does not satisfy the function-first system"
                    );
                    witnessed += 1;
                }
            }
            if witnessed >= 2 {
                break;
            }
        }
        assert!(witnessed >= 1, "no decomposable multiple of g found");
    }

    #[test]
    fn quotient_pipeline_recovers_the_logarithm() {
        let glv = tiny();
        let rep = run_quotient_experiment(&glv, 1, 20_000);
        assert!(rep.quotient.solved && rep.quotient.correct == Some(true));
        assert!(rep.control.solved && rep.control.correct == Some(true));
        assert_eq!(rep.base, 3 * rep.orbits);
        assert!(rep.quotient.residuals_at_solve <= rep.control.residuals_at_solve);
    }

    #[test]
    fn canonical_generation_saves_pair_solver_calls() {
        let glv = tiny();
        let rep = run_canonical_experiment(&glv, 1, 40, 300, 50);
        assert_eq!(rep.random.verification_mismatches, 0);
        assert_eq!(rep.pairs.verification_mismatches, 0);
        assert!(rep.pairs.duplicates > 0);
        assert_eq!(
            rep.pairs.solver_calls_saved,
            rep.pairs.residuals - rep.pairs.distinct_orbits
        );
    }
}
