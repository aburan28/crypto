//! # Index-calculus boundary ledger: one unit, three regimes, every phase priced.
//!
//! This module answers one question with numbers instead of asymptotics:
//! **how far is every index-calculus pipeline in this repository from the
//! boundary it cannot cross, and which phase puts it there?**  It runs the
//! same discrete logarithm three ways — a generic prime-field curve, a
//! generic binary curve, and a Koblitz curve — and prices every phase of
//! each in one unit against two boundaries, following the reporting rule
//! in `AGENTS.md` (boundary, table, ratio).
//!
//! ## The unit
//!
//! ```text
//!     S = total operations / √r
//! ```
//!
//! where `r` is the prime order of the subgroup the logarithm lives in and
//! an *operation* is one affine group addition on the curve in question
//! (a doubling counts as one).  Everything the method does is inside the
//! number: building the factor base, drawing the targets, running the
//! decomposition oracle, the pair table, the linear algebra, the
//! verification.  Costs that are not group additions — a modular square
//! root, an Artin–Schreier solve, one pair of the pairs-and-solve loop, a
//! hash probe, one multiply-subtract of the elimination — are counted
//! **exactly** in their own native unit and converted with a factor
//! **measured on the host at run time** (nanoseconds per native unit over
//! nanoseconds per addition).  The native counts and the factors are both
//! in the report, so a reader can re-convert; the count is the
//! measurement and the factor is the only thing hardware touches.
//!
//! ## The boundaries
//!
//! - **Floor.**  A generic algorithm on a group of prime order `r` with an
//!   automorphism group of order `A` needs about `√(πr/2A)` operations
//!   (Shoup's bound gives the `√r`; the automorphisms give the `√A`), so
//!   the floor in the unit is `S_floor = √(π/2A)`: `0.886` when only
//!   negation is available, `√(π/4n)` on a Koblitz curve over `F_{2^n}`.
//!   The relation phase has its own, counting, floor: a base of `F`
//!   signed points has at most `C(F+m−1, m)` `m`-sums, so a random target
//!   decomposes with probability at most `C(F+m−1, m)/#E`, and the
//!   `K + 1` relations the linear algebra needs cost at least
//!   `(K+1)·#E/C(F+m−1, m)` targets.  The report carries the measured
//!   trials over that number; a value near one says the base is as good
//!   as a random one and nothing structural is being exploited.
//! - **Reference.**  Pollard rho, run on the same instance in the same
//!   process with the same accounting, counted exactly: an
//!   r-adding walk with distinguished points for the prime and generic
//!   binary regimes (`A = 1` measured, `A = 2` stated as the floor), and
//!   the repository's signed-Frobenius walk (`A = 2n`, measured) on the
//!   Koblitz curves.  Every rho run verifies `[d]G = Q`; a rho that does
//!   not finish is not a reference.
//!
//! ## The regimes and their variants
//!
//! | regime | curve | factor base | oracles |
//! |:--|:--|:--|:--|
//! | `prime` | `y² = x³ + ax + b` over `F_p`, prime order | smallest abscissae | Semaev `S₃` roots (`m = 2`); direct subtraction (`m = 2`); meet in the middle (`m = 3`) |
//! | `char2` | random `y² + xy = x³ + ax² + b` over `F_{2^n}` | low-order subspace `⟨1, z, …, z^{l−1}⟩`, `l = n/3` | Semaev `S₄` pairs-and-solve (`m = 3`); meet in the middle (`m = 3`) |
//! | `koblitz` | `K_a` over `F_{2^n}` | Frobenius-invariant subspace or orbit union | meet in the middle over signed-orbit columns (`m = 3`); the same without the orbit fold; `S₄` pairs-and-solve on the invariant subspace |
//!
//! All three share one relation loop (`R = [a]G + [b]Q`, oracle, lift,
//! one row) and one incremental Gauss–Jordan over `Z/rZ` that stops the
//! moment the target's column is pinned, so the only thing that differs
//! between a `char2` row and a `koblitz` row on the same curve is the
//! column map — which is exactly what the Frobenius is supposed to buy.
//!
//! ## What the numbers are for
//!
//! Exponents fitted over the ladder (`ops ∝ r^α`, per phase and in
//! total) say which phase dominates and where; the `S` column and its
//! ratios say how far each variant is from rho *today*.  Nothing here is
//! a claim about a deployed curve.  The largest instance is a 41-bit
//! Koblitz subgroup, and the point of the exercise is to have a concrete
//! boundary to iterate against rather than an opinion.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::time::Instant;

use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use crate::binary_ecc::{BinaryCurve, BinaryPoint, IrreduciblePoly};
use crate::cryptanalysis::koblitz_fast::{FastCurve, FastPoint};
use crate::cryptanalysis::koblitz_index_calculus::{
    all_factors_of_x_n_minus_1, available_subspace_dimensions,
    build_frobenius_factor_base_from_divisor, build_frobenius_union_factor_base,
    find_irreducible_sparse, koblitz_signed_frobenius_rho_reference,
    saturate_factor_base_two_torsion, FrobeniusFactorBase, KoblitzCurve,
    KoblitzSignedRhoOptions,
};
use crate::cryptanalysis::research_bench::{bench_curves, linear_fit};
use crate::cryptanalysis::residual_walk::is_prime_u64;
use crate::cryptanalysis::semaev_decomp::{Gf2, SubspaceOracle};

// ── Units, ledgers, boundaries ─────────────────────────────────────

/// Exact group-operation ledger.  A doubling is an operation like an
/// addition (both are one inversion and a few multiplications in affine
/// coordinates); they are kept apart so the walk's shape stays visible.
#[derive(Clone, Copy, Debug, Default, Serialize, PartialEq, Eq)]
pub struct GroupOps {
    pub adds: u64,
    pub doubles: u64,
    pub scalar_mults: u64,
}

impl GroupOps {
    /// Group-addition equivalents: additions plus doublings.
    pub fn gae(&self) -> f64 {
        (self.adds + self.doubles) as f64
    }
    fn merge(&mut self, other: GroupOps) {
        self.adds += other.adds;
        self.doubles += other.doubles;
        self.scalar_mults += other.scalar_mults;
    }
}

/// Nanoseconds per native unit, measured on the host for the instance at
/// hand.  Every conversion the report performs goes through one of these
/// and the report carries them, so the counts stay the measurement.
#[derive(Clone, Debug, Default, Serialize)]
pub struct Calibration {
    /// One affine group addition on this curve (random operands).
    pub ns_per_add: f64,
    /// One affine doubling.
    pub ns_per_double: f64,
    /// One modular square root (Tonelli–Shanks), prime regime.
    pub ns_per_sqrt: Option<f64>,
    /// One Artin–Schreier solve `u² + u = c`, binary regimes.
    pub ns_per_as_solve: Option<f64>,
    /// One pair of the `S₄` pairs-and-solve loop, binary regimes.
    pub ns_per_s4_pair: Option<f64>,
    /// One probe into the pair table (hash lookup), measured on the real
    /// table when one is built.
    pub ns_per_lookup: f64,
    /// One multiply-subtract of the modular elimination.
    pub ns_per_row_op: f64,
    /// A 64-bit word XOR of the Boolean Macaulay elimination (algebraic
    /// oracle pricing); measured once per process.
    pub ns_per_word_xor: Option<f64>,
    /// One Legendre symbol (Jacobi algorithm), prime regime.
    pub ns_per_legendre: Option<f64>,
    /// One modular inversion (extended Euclid), prime regime.
    pub ns_per_inversion: Option<f64>,
    /// One Frobenius map of a point (two squarings), binary regimes.
    pub ns_per_frobenius: Option<f64>,
}

impl Calibration {
    /// Group-addition equivalents of `count` native units at
    /// `ns_per_unit`.
    pub fn gae(&self, count: u64, ns_per_unit: f64) -> f64 {
        if self.ns_per_add <= 0.0 {
            return count as f64;
        }
        count as f64 * ns_per_unit / self.ns_per_add
    }
}

/// `√(πr/2A)`: expected operations of a generic collision search on a
/// group of order `r` with `A` usable automorphisms.
pub fn generic_floor_ops(r: f64, automorphisms: f64) -> f64 {
    (std::f64::consts::PI * r / (2.0 * automorphisms)).sqrt()
}

/// The same floor in the unit: `√(π/2A)`.
pub fn generic_floor_s(automorphisms: f64) -> f64 {
    (std::f64::consts::PI / (2.0 * automorphisms)).sqrt()
}

/// `C(n, k)` as a float, for the counting ceiling.
fn binomial_f64(n: f64, k: u32) -> f64 {
    let mut acc = 1.0f64;
    for i in 0..k {
        acc *= (n - i as f64) / (i as f64 + 1.0);
    }
    acc.max(0.0)
}

/// Upper bound on the probability that a uniform target of a space of
/// size `space` is an `m`-sum of `signed_points` base points:
/// `min(1, C(F + m − 1, m) / space)`.
pub fn decomposition_probability_ceiling(signed_points: u64, m: u32, space: f64) -> f64 {
    let sums = binomial_f64(signed_points as f64 + m as f64 - 1.0, m);
    (sums / space).min(1.0)
}

/// Fewest targets that can yield `columns + 1` relations under the
/// counting ceiling.
pub fn trials_floor(columns: u64, signed_points: u64, m: u32, space: f64) -> f64 {
    (columns as f64 + 1.0) / decomposition_probability_ceiling(signed_points, m, space)
}

/// Cost of one phase: what it did (exact), what it took (wall) and what
/// it comes to in the unit (converted).
#[derive(Clone, Debug, Default, Serialize)]
pub struct PhaseCost {
    pub wall_ns: u64,
    pub group_ops: GroupOps,
    /// Native counters by name: `sqrt_solves`, `as_solves`, `s4_pairs`,
    /// `lookups`, `row_ops`, `abscissae_scanned`, `trials`, …
    pub native: std::collections::BTreeMap<String, u64>,
    /// Group-addition equivalents after conversion.
    pub gae: f64,
}

impl PhaseCost {
    fn count(&mut self, name: &str, by: u64) {
        *self.native.entry(name.to_string()).or_insert(0) += by;
    }
    fn get(&self, name: &str) -> u64 {
        self.native.get(name).copied().unwrap_or(0)
    }
}

/// One variant of one regime on one instance, priced end to end.
#[derive(Clone, Debug, Serialize)]
pub struct VariantResult {
    pub name: String,
    pub oracle: String,
    pub summands: u32,
    pub factor_base: PhaseCost,
    /// Base size (signed points), distinct abscissae, relation columns,
    /// subspace dimension when the base is a subspace.
    pub signed_points: u64,
    pub abscissae: u64,
    pub columns: u64,
    pub dimension: Option<u32>,
    pub relations: PhaseCost,
    pub linear_algebra: PhaseCost,
    pub verify: PhaseCost,
    pub trials: u64,
    pub relations_found: u64,
    pub independent: u64,
    pub dependent: u64,
    pub rank: u64,
    /// Fewest targets the counting ceiling allows for `K + 1` relations.
    pub trials_floor: f64,
    /// Measured trials over that floor.  Below one when the elimination
    /// pinned the logarithm before `K + 1` relations (a repeated
    /// decomposition does that at small sizes).
    pub trials_over_floor: f64,
    /// Upper bound on the decomposition probability, `C(F+m−1, m)/#E`.
    pub decomposition_probability_ceiling: f64,
    /// Measured relations per trial over that ceiling: how much of the
    /// counting bound the base actually realises.
    pub yield_over_ceiling: f64,
    pub total_gae: f64,
    pub total_wall_ns: u64,
    /// `total_gae / √r`.
    pub s: f64,
    /// `total_wall / (ns_per_add · √r)`: the same number from wall
    /// time, a practicality note, never the metric.
    pub s_wall: f64,
    pub ratio_to_rho: f64,
    pub ratio_to_floor: f64,
    pub recovered: Option<u64>,
    pub verified: bool,
}

/// A counted rho run.
#[derive(Clone, Debug, Serialize)]
pub struct RhoResult {
    pub method: String,
    pub automorphisms: u32,
    pub group_ops: GroupOps,
    pub steps: u64,
    pub walks: u64,
    pub distinguished_points: u64,
    pub gae: f64,
    /// Everything: setup, walk, verification.
    pub s: f64,
    /// Walk steps alone, the number the `√(πr/2A)` floor is about.
    pub s_walk: f64,
    pub expected_steps: f64,
    pub steps_over_expected: f64,
    pub wall_ns: u64,
    pub recovered: Option<u64>,
    pub verified: bool,
}

/// One instance of one regime: the curve, the boundaries, the reference
/// and every variant, with the repeats folded to means and kept raw.
#[derive(Clone, Debug, Serialize)]
pub struct RegimeInstance {
    pub regime: String,
    pub curve: serde_json::Value,
    pub r: u64,
    pub log2_r: f64,
    pub group_order: f64,
    pub cofactor: u64,
    /// Automorphisms a generic algorithm may use on this curve.
    pub automorphisms_generic: u32,
    pub floor_s: f64,
    pub floor_ops: f64,
    pub calibration: Calibration,
    /// Every rho run (best reference first within a repeat).
    pub rho: Vec<RhoResult>,
    /// Mean `S` of the best-reference rho over the repeats.
    pub rho_s_mean: f64,
    pub rho_verified_all: bool,
    /// Every variant run; several rows per variant when `repeats > 1`.
    pub variants: Vec<VariantResult>,
    pub seeds: Vec<u64>,
    pub targets: Vec<u64>,
}

/// A fitted exponent `ops ∝ r^α` for one phase of one variant, and the
/// same fit against the full group order `#E` (the two differ on
/// Koblitz curves with large cofactors, where the relation phase is
/// priced by `#E` and the unit by `r`).
#[derive(Clone, Debug, Serialize)]
pub struct ExponentFit {
    pub regime: String,
    pub variant: String,
    pub phase: String,
    pub alpha: f64,
    pub r_squared: f64,
    pub alpha_group_order: f64,
    pub r_squared_group_order: f64,
    pub points: usize,
}

/// The whole ledger.
#[derive(Clone, Debug, Serialize)]
pub struct BoundaryLedger {
    pub unit: String,
    pub floor: String,
    pub reference: String,
    pub instances: Vec<RegimeInstance>,
    pub fits: Vec<ExponentFit>,
}

// ── A fast hasher for point keys ───────────────────────────────────

/// Multiplicative mixing of one `u64`; the keys are already packed
/// points, and SipHash would dominate a table probe.
#[derive(Default, Clone, Copy)]
pub struct MixHasher(u64);

impl Hasher for MixHasher {
    fn finish(&self) -> u64 {
        let mut h = self.0;
        h ^= h >> 32;
        h = h.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        h ^= h >> 29;
        h
    }
    fn write(&mut self, bytes: &[u8]) {
        for chunk in bytes.chunks(8) {
            let mut w = [0u8; 8];
            w[..chunk.len()].copy_from_slice(chunk);
            self.write_u64(u64::from_le_bytes(w));
        }
    }
    fn write_u64(&mut self, v: u64) {
        self.0 = (self.0.rotate_left(5) ^ v).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
}

type FastMap<V> = HashMap<u64, V, BuildHasherDefault<MixHasher>>;

fn fast_map<V>(capacity: usize) -> FastMap<V> {
    HashMap::with_capacity_and_hasher(capacity, BuildHasherDefault::default())
}

// ── Modular arithmetic on u64 ──────────────────────────────────────

#[inline]
fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

#[inline]
fn addmod(a: u64, b: u64, m: u64) -> u64 {
    let s = a as u128 + b as u128;
    (s % m as u128) as u64
}

#[inline]
fn submod(a: u64, b: u64, m: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        m - (b - a)
    }
}

fn powmod(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1u64 % m;
    a %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, m);
        }
        a = mulmod(a, a, m);
        e >>= 1;
    }
    r
}

/// Modular inverse by the extended Euclidean algorithm.
pub fn invmod(a: u64, m: u64) -> Option<u64> {
    let (mut old_r, mut r) = (a as i128 % m as i128, m as i128);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        (old_r, r) = (r, old_r - q * r);
        (old_s, s) = (s, old_s - q * s);
    }
    if old_r != 1 {
        return None;
    }
    let mut v = old_s % m as i128;
    if v < 0 {
        v += m as i128;
    }
    Some(v as u64)
}

// ── Incremental Gauss–Jordan over Z/rZ ─────────────────────────────

/// What adding a row did.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RowStatus {
    Independent,
    Dependent,
    Inconsistent,
}

/// Reduced row-echelon form maintained as relations arrive, over a
/// prime modulus that fits a word.  Every multiply-subtract on an entry
/// is counted, so the linear algebra is priced in its own native unit.
pub struct IncrementalGauss {
    modulus: u64,
    cols: usize,
    rows: Vec<Vec<u64>>,
    rhs: Vec<u64>,
    pivot_of_row: Vec<usize>,
    row_of_pivot: Vec<Option<usize>>,
    /// Multiply-subtracts performed on matrix entries.
    pub row_ops: u64,
    pub dependent: u64,
    pub inconsistent: u64,
}

impl IncrementalGauss {
    pub fn new(cols: usize, modulus: u64) -> Self {
        Self {
            modulus,
            cols,
            rows: Vec::new(),
            rhs: Vec::new(),
            pivot_of_row: Vec::new(),
            row_of_pivot: vec![None; cols],
            row_ops: 0,
            dependent: 0,
            inconsistent: 0,
        }
    }

    pub fn rank(&self) -> usize {
        self.rows.len()
    }

    pub fn columns(&self) -> usize {
        self.cols
    }

    /// Add one equation `Σ row[j]·x_j ≡ rhs (mod r)`.
    pub fn add_row(&mut self, mut row: Vec<u64>, mut rhs: u64) -> RowStatus {
        let m = self.modulus;
        debug_assert_eq!(row.len(), self.cols);
        // Reduce against every existing pivot.
        for i in 0..self.rows.len() {
            let pc = self.pivot_of_row[i];
            let f = row[pc];
            if f == 0 {
                continue;
            }
            let pivot = &self.rows[i];
            for j in 0..self.cols {
                let v = pivot[j];
                if v != 0 {
                    row[j] = submod(row[j], mulmod(f, v, m), m);
                    self.row_ops += 1;
                }
            }
            rhs = submod(rhs, mulmod(f, self.rhs[i], m), m);
            self.row_ops += 1;
        }
        let Some(pc) = (0..self.cols).find(|&j| row[j] != 0) else {
            if rhs == 0 {
                self.dependent += 1;
                return RowStatus::Dependent;
            }
            self.inconsistent += 1;
            return RowStatus::Inconsistent;
        };
        // Normalise.
        let inv = invmod(row[pc], m).expect("prime modulus");
        for v in row.iter_mut() {
            if *v != 0 {
                *v = mulmod(*v, inv, m);
                self.row_ops += 1;
            }
        }
        rhs = mulmod(rhs, inv, m);
        self.row_ops += 1;
        // Eliminate the new pivot column from every existing row.
        let nnz: Vec<usize> = (0..self.cols).filter(|&j| row[j] != 0).collect();
        for i in 0..self.rows.len() {
            let f = self.rows[i][pc];
            if f == 0 {
                continue;
            }
            for &j in &nnz {
                let v = mulmod(f, row[j], m);
                self.rows[i][j] = submod(self.rows[i][j], v, m);
                self.row_ops += 1;
            }
            self.rhs[i] = submod(self.rhs[i], mulmod(f, rhs, m), m);
            self.row_ops += 1;
        }
        self.row_of_pivot[pc] = Some(self.rows.len());
        self.pivot_of_row.push(pc);
        self.rows.push(row);
        self.rhs.push(rhs);
        RowStatus::Independent
    }

    /// The value of unknown `col`, if the reduced system determines it
    /// on its own (a pivot row with no other non-zero entry).
    pub fn pinned(&self, col: usize) -> Option<u64> {
        let i = self.row_of_pivot[col]?;
        let row = &self.rows[i];
        if row.iter().enumerate().all(|(j, &v)| j == col || v == 0) {
            Some(self.rhs[i])
        } else {
            None
        }
    }
}

// ── A counted group ────────────────────────────────────────────────

/// The little a collision search and a relation loop need from a group,
/// with every operation charged to a ledger.
pub trait CountedGroup {
    type Elt: Copy + PartialEq + std::fmt::Debug;
    fn identity(&self) -> Self::Elt;
    fn is_identity(&self, p: &Self::Elt) -> bool;
    fn add(&self, ops: &mut GroupOps, p: Self::Elt, q: Self::Elt) -> Self::Elt;
    fn double(&self, ops: &mut GroupOps, p: Self::Elt) -> Self::Elt;
    fn neg(&self, p: Self::Elt) -> Self::Elt;
    /// A key that identifies the element (both walks must agree on it).
    fn key(&self, p: &Self::Elt) -> u64;
    /// `[k]P` by double-and-add, counted exactly.
    fn mul(&self, ops: &mut GroupOps, p: Self::Elt, k: u64) -> Self::Elt {
        ops.scalar_mults += 1;
        if k == 0 || self.is_identity(&p) {
            return self.identity();
        }
        let bits = 64 - k.leading_zeros();
        let mut acc = p;
        for i in (0..bits - 1).rev() {
            acc = self.double(ops, acc);
            if (k >> i) & 1 == 1 {
                acc = self.add(ops, acc, p);
            }
        }
        acc
    }
}

fn mix(v: u64) -> u64 {
    let mut h = v ^ 0x2545_F491_4F6C_DD1D;
    h = h.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    h ^= h >> 31;
    h = h.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    h ^= h >> 29;
    h
}

/// **Pollard rho, counted.**  Van Oorschot–Wiener: independent
/// r-adding walks (32 jumps `[a_j]G + [b_j]Q`) run until a distinguished
/// point, stored with their coefficients; two walks reaching one point
/// with different coefficients give the logarithm, which is verified as
/// `[d]G = Q` before it is returned.  No automorphism is used
/// (`A = 1`), so the expected step count is `√(πr/2)`.
pub fn rho_reference<G: CountedGroup>(
    g: &G,
    generator: G::Elt,
    target: G::Elt,
    r: u64,
    seed: u64,
    max_steps: u64,
) -> RhoResult {
    let start = Instant::now();
    let mut ops = GroupOps::default();
    let mut rng = StdRng::seed_from_u64(seed);
    let bits = 64 - r.leading_zeros();
    // Walks of about √r / 8 steps, so the per-walk start (two scalar
    // multiplications) stays a small fraction of the walk.
    let dp_bits = (bits / 2).saturating_sub(3).min(24);
    let dp_mask = if dp_bits == 0 { 0 } else { (1u64 << dp_bits) - 1 };
    // Sixteen jumps (Teske's r-adding walk needs about that many to be
    // close to a random walk) with 12-bit coefficients: the jumps only
    // need to be group elements with known coefficients, and short
    // scalars keep the table's setup from dominating small instances.
    let jump_count = 16usize;
    let short = r.min(1u64 << 12);
    let jumps: Vec<(G::Elt, u64, u64)> = (0..jump_count)
        .map(|_| {
            let a = rng.gen_range(1..short);
            let b = rng.gen_range(1..short);
            let ag = g.mul(&mut ops, generator, a);
            let bq = g.mul(&mut ops, target, b);
            (g.add(&mut ops, ag, bq), a, b)
        })
        .collect();
    let mut table: FastMap<(u64, u64)> = fast_map(1024);
    let expected = generic_floor_ops(r as f64, 1.0);
    let mut steps = 0u64;
    let mut walks = 0u64;
    let mut dps = 0u64;
    let mut recovered = None;
    let walk_cap = (1u64 << dp_bits) * 20 + 64;

    // A walk steps first and tests for a distinguished point after, so
    // even when every point is distinguished (tiny groups) each walk
    // does work and the step cap is reached.
    'outer: while steps < max_steps && walks < max_steps {
        walks += 1;
        let mut a = rng.gen_range(1..r);
        let mut b = rng.gen_range(1..r);
        let ag = g.mul(&mut ops, generator, a);
        let bq = g.mul(&mut ops, target, b);
        let mut x = g.add(&mut ops, ag, bq);
        let mut len = 0u64;
        loop {
            if g.is_identity(&x) {
                break;
            }
            let h = mix(g.key(&x));
            let j = ((h >> 24) % jump_count as u64) as usize;
            let (m, aj, bj) = jumps[j];
            x = g.add(&mut ops, x, m);
            a = addmod(a, aj, r);
            b = addmod(b, bj, r);
            steps += 1;
            len += 1;
            if g.is_identity(&x) {
                break;
            }
            let h = mix(g.key(&x));
            if h & dp_mask == 0 {
                dps += 1;
                match table.get(&g.key(&x)) {
                    Some(&(a2, b2)) => {
                        if b2 != b {
                            // a + b d = a2 + b2 d  ⇒  d = (a − a2)/(b2 − b)
                            let num = submod(a, a2, r);
                            let den = submod(b2, b, r);
                            if let Some(inv) = invmod(den, r) {
                                let d = mulmod(num, inv, r);
                                let check = g.mul(&mut ops, generator, d);
                                if check == target {
                                    recovered = Some(d);
                                    break 'outer;
                                }
                            }
                        }
                        break;
                    }
                    None => {
                        table.insert(g.key(&x), (a, b));
                        break;
                    }
                }
            }
            if len > walk_cap || steps >= max_steps {
                break;
            }
        }
    }
    let wall_ns = start.elapsed().as_nanos() as u64;
    let gae = ops.gae();
    RhoResult {
        method: "r-adding walk, distinguished points, 16 jumps".into(),
        automorphisms: 1,
        group_ops: ops,
        steps,
        walks,
        distinguished_points: dps,
        gae,
        s: gae / (r as f64).sqrt(),
        s_walk: steps as f64 / (r as f64).sqrt(),
        expected_steps: expected,
        steps_over_expected: steps as f64 / expected,
        wall_ns,
        verified: recovered.is_some(),
        recovered,
    }
}

// ── Prime-field curves on one word ─────────────────────────────────

/// A point on `y² = x³ + ax + b` over `F_p`, `p < 2^63`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PrimePoint {
    pub x: u64,
    pub y: u64,
    pub infinity: bool,
}

impl PrimePoint {
    pub const INFINITY: Self = Self {
        x: 0,
        y: 0,
        infinity: true,
    };
    pub fn affine(x: u64, y: u64) -> Self {
        Self {
            x,
            y,
            infinity: false,
        }
    }
}

/// `y² = x³ + ax + b` over `F_p` in single-word arithmetic.
#[derive(Clone, Debug, Serialize)]
pub struct PrimeCurve {
    pub p: u64,
    pub a: u64,
    pub b: u64,
}

impl PrimeCurve {
    pub fn is_on_curve(&self, pt: PrimePoint) -> bool {
        if pt.infinity {
            return true;
        }
        let p = self.p;
        let lhs = mulmod(pt.y, pt.y, p);
        let rhs = addmod(
            addmod(mulmod(mulmod(pt.x, pt.x, p), pt.x, p), mulmod(self.a, pt.x, p), p),
            self.b,
            p,
        );
        lhs == rhs
    }

    /// `x³ + ax + b`.
    pub fn rhs(&self, x: u64) -> u64 {
        let p = self.p;
        addmod(
            addmod(mulmod(mulmod(x, x, p), x, p), mulmod(self.a, x, p), p),
            self.b,
            p,
        )
    }

    /// Legendre symbol of `a` modulo the prime `p`, by the Jacobi
    /// algorithm.
    pub fn legendre(&self, a: u64) -> i32 {
        jacobi(a % self.p, self.p)
    }

    /// A square root modulo `p` (Tonelli–Shanks), or `None` for a
    /// non-residue.
    pub fn sqrt(&self, n: u64) -> Option<u64> {
        let p = self.p;
        let n = n % p;
        if n == 0 {
            return Some(0);
        }
        if p % 4 == 3 {
            let r = powmod(n, (p + 1) / 4, p);
            return (mulmod(r, r, p) == n).then_some(r);
        }
        if powmod(n, (p - 1) / 2, p) != 1 {
            return None;
        }
        let mut q = p - 1;
        let mut s = 0u32;
        while q % 2 == 0 {
            q /= 2;
            s += 1;
        }
        let mut z = 2u64;
        while powmod(z, (p - 1) / 2, p) != p - 1 {
            z += 1;
        }
        let mut m = s;
        let mut c = powmod(z, q, p);
        let mut t = powmod(n, q, p);
        let mut r = powmod(n, (q + 1) / 2, p);
        loop {
            if t == 1 {
                return Some(r);
            }
            let mut i = 0u32;
            let mut tt = t;
            while tt != 1 {
                tt = mulmod(tt, tt, p);
                i += 1;
                if i == m {
                    return None;
                }
            }
            let b = powmod(c, 1u64 << (m - i - 1), p);
            m = i;
            c = mulmod(b, b, p);
            t = mulmod(t, c, p);
            r = mulmod(r, b, p);
        }
    }

    fn add_raw(&self, a: PrimePoint, b: PrimePoint) -> PrimePoint {
        if a.infinity {
            return b;
        }
        if b.infinity {
            return a;
        }
        let p = self.p;
        if a.x == b.x {
            if addmod(a.y, b.y, p) == 0 {
                return PrimePoint::INFINITY;
            }
            return self.double_raw(a);
        }
        let inv = invmod(submod(b.x, a.x, p), p).expect("prime field");
        let lambda = mulmod(submod(b.y, a.y, p), inv, p);
        let x3 = submod(submod(mulmod(lambda, lambda, p), a.x, p), b.x, p);
        let y3 = submod(mulmod(lambda, submod(a.x, x3, p), p), a.y, p);
        PrimePoint::affine(x3, y3)
    }

    fn double_raw(&self, a: PrimePoint) -> PrimePoint {
        if a.infinity || a.y == 0 {
            return PrimePoint::INFINITY;
        }
        let p = self.p;
        let num = addmod(mulmod(3, mulmod(a.x, a.x, p), p), self.a, p);
        let inv = invmod(mulmod(2, a.y, p), p).expect("prime field");
        let lambda = mulmod(num, inv, p);
        let x3 = submod(mulmod(lambda, lambda, p), mulmod(2, a.x, p), p);
        let y3 = submod(mulmod(lambda, submod(a.x, x3, p), p), a.y, p);
        PrimePoint::affine(x3, y3)
    }

    /// `#E(F_p)` by the Legendre symbol at every abscissa: `O(p)`.
    pub fn point_count(&self) -> u64 {
        let mut count = 1u64;
        for x in 0..self.p {
            count += (1 + self.legendre(self.rhs(x))) as u64;
        }
        count
    }
}

fn jacobi(mut a: u64, mut n: u64) -> i32 {
    let mut result = 1i32;
    a %= n;
    while a != 0 {
        while a % 2 == 0 {
            a /= 2;
            if n % 8 == 3 || n % 8 == 5 {
                result = -result;
            }
        }
        std::mem::swap(&mut a, &mut n);
        if a % 4 == 3 && n % 4 == 3 {
            result = -result;
        }
        a %= n;
    }
    if n == 1 {
        result
    } else {
        0
    }
}

impl CountedGroup for PrimeCurve {
    type Elt = PrimePoint;
    fn identity(&self) -> PrimePoint {
        PrimePoint::INFINITY
    }
    fn is_identity(&self, p: &PrimePoint) -> bool {
        p.infinity
    }
    fn add(&self, ops: &mut GroupOps, p: PrimePoint, q: PrimePoint) -> PrimePoint {
        ops.adds += 1;
        self.add_raw(p, q)
    }
    fn double(&self, ops: &mut GroupOps, p: PrimePoint) -> PrimePoint {
        ops.doubles += 1;
        self.double_raw(p)
    }
    fn neg(&self, p: PrimePoint) -> PrimePoint {
        if p.infinity {
            p
        } else {
            PrimePoint::affine(p.x, if p.y == 0 { 0 } else { self.p - p.y })
        }
    }
    fn key(&self, p: &PrimePoint) -> u64 {
        if p.infinity {
            0
        } else {
            ((p.x + 1) << 1) | u64::from(p.y > self.p - p.y)
        }
    }
}

/// A prime-field instance: curve, prime-order subgroup, generator.
#[derive(Clone, Debug, Serialize)]
pub struct PrimeInstance {
    pub name: String,
    pub curve: PrimeCurve,
    pub group_order: u64,
    pub r: u64,
    pub cofactor: u64,
    pub generator: (u64, u64),
}

impl PrimeInstance {
    pub fn generator_point(&self) -> PrimePoint {
        PrimePoint::affine(self.generator.0, self.generator.1)
    }
}

/// The repository's bench roster (`research_bench::bench_curves`),
/// converted to single words, keyed by the roster's bit size.
pub fn roster_prime_instance(bits: u32) -> Option<PrimeInstance> {
    let (_, c) = bench_curves().into_iter().find(|(b, _)| *b == bits)?;
    let p = c.p.to_u64()?;
    let curve = PrimeCurve {
        p,
        a: c.a.to_u64()?,
        b: c.b.to_u64()?,
    };
    let n = c.n.to_u64()?;
    Some(PrimeInstance {
        name: c.name.to_string(),
        curve,
        group_order: n * c.h as u64,
        r: n,
        cofactor: c.h as u64,
        generator: (c.gx.to_u64()?, c.gy.to_u64()?),
    })
}

/// **Search a prime-order curve** of about `bits` bits: random prime
/// `p`, random `(a, b)`, `O(p)` point count, until `#E` is prime.
/// Deterministic in the seed; the count of curves tried is not returned
/// because it is not a cost of the attack.
pub fn find_prime_order_curve(bits: u32, seed: u64) -> PrimeInstance {
    assert!((8..=32).contains(&bits));
    let mut rng = StdRng::seed_from_u64(seed ^ (bits as u64) << 40);
    loop {
        let p = loop {
            let cand = (rng.gen::<u64>() & ((1u64 << bits) - 1)) | (1u64 << (bits - 1)) | 1;
            if is_prime_u64(cand) {
                break cand;
            }
        };
        for _ in 0..64 {
            let a = rng.gen_range(0..p);
            let b = rng.gen_range(1..p);
            let curve = PrimeCurve { p, a, b };
            // Non-singular?
            let disc = addmod(
                mulmod(4, mulmod(mulmod(a, a, p), a, p), p),
                mulmod(27, mulmod(b, b, p), p),
                p,
            );
            if disc == 0 {
                continue;
            }
            let order = curve.point_count();
            if !is_prime_u64(order) {
                continue;
            }
            // Any point generates.
            let g = loop {
                let x = rng.gen_range(0..p);
                if let Some(y) = curve.sqrt(curve.rhs(x)) {
                    if curve.legendre(curve.rhs(x)) >= 0 {
                        break PrimePoint::affine(x, y);
                    }
                }
            };
            let mut ops = GroupOps::default();
            debug_assert!(curve.mul(&mut ops, g, order).infinity);
            return PrimeInstance {
                name: format!("generated-{bits}bit-{p}"),
                curve,
                group_order: order,
                r: order,
                cofactor: 1,
                generator: (g.x, g.y),
            };
        }
    }
}

/// Semaev's third summation polynomial for `y² = x³ + ax + b` as a
/// quadratic in its last argument: `S₃(x₁, x₂, X) = qa·X² + qb·X + qc`.
pub fn semaev_s3_quadratic(curve: &PrimeCurve, x1: u64, x2: u64) -> (u64, u64, u64) {
    let p = curve.p;
    let d = submod(x1, x2, p);
    let qa = mulmod(d, d, p);
    let s = addmod(x1, x2, p);
    let pr = mulmod(x1, x2, p);
    // qb = −2((x1 + x2)(x1x2 + a) + 2b)
    let inner = addmod(mulmod(s, addmod(pr, curve.a, p), p), mulmod(2, curve.b, p), p);
    let qb = submod(0, mulmod(2, inner, p), p);
    // qc = (x1x2 − a)² − 4b(x1 + x2)
    let t = submod(pr, curve.a, p);
    let qc = submod(mulmod(t, t, p), mulmod(mulmod(4, curve.b, p), s, p), p);
    (qa, qb, qc)
}

// ── Binary curves on one word ──────────────────────────────────────

/// Solver for `u² + u = c` over `F_{2^n}`: the map is `F₂`-linear with
/// kernel `{0, 1}`, so one echelon form of its matrix answers every
/// query in `O(n)` word operations.
#[derive(Clone, Debug)]
pub struct ArtinSchreier {
    n: u32,
    by_lead: Vec<Option<(u64, u64)>>,
    trace_mask: u64,
}

impl ArtinSchreier {
    pub fn new(gf: &Gf2) -> Self {
        let n = gf.n;
        let mut by_lead: Vec<Option<(u64, u64)>> = vec![None; 64];
        for i in 0..n {
            let v = 1u64 << i;
            let mut img = gf.sqr(v) ^ v;
            let mut pre = v;
            while img != 0 {
                let lb = 63 - img.leading_zeros() as usize;
                match by_lead[lb] {
                    Some((pi, pp)) => {
                        img ^= pi;
                        pre ^= pp;
                    }
                    None => {
                        by_lead[lb] = Some((img, pre));
                        break;
                    }
                }
            }
        }
        // Tr(z^j) for every basis element, as a mask.
        let mut trace_mask = 0u64;
        for j in 0..n {
            let mut t = 1u64 << j;
            let mut acc = 0u64;
            for _ in 0..n {
                acc ^= t;
                t = gf.sqr(t);
            }
            debug_assert!(acc == 0 || acc == 1);
            if acc == 1 {
                trace_mask |= 1u64 << j;
            }
        }
        Self {
            n,
            by_lead,
            trace_mask,
        }
    }

    /// Absolute trace of `v`.
    #[inline]
    pub fn trace(&self, v: u64) -> u64 {
        ((v & self.trace_mask).count_ones() & 1) as u64
    }

    /// A solution of `u² + u = c`, or `None` when `Tr(c) = 1`.
    pub fn solve(&self, mut c: u64) -> Option<u64> {
        let mut u = 0u64;
        while c != 0 {
            let lb = 63 - c.leading_zeros() as usize;
            let (pi, pp) = self.by_lead[lb]?;
            c ^= pi;
            u ^= pp;
        }
        Some(u)
    }

    pub fn degree(&self) -> u32 {
        self.n
    }
}

/// A binary-curve instance in single-word arithmetic: the curve, the
/// prime-order subgroup, a generator, and — for a Koblitz curve — the
/// repository's [`KoblitzCurve`] carrying its Frobenius structure.
pub struct BinaryInstance {
    pub name: String,
    pub n: u32,
    pub irreducible: IrreduciblePoly,
    pub fast: FastCurve,
    pub gf: Gf2,
    pub artin_schreier: ArtinSchreier,
    pub a: u64,
    pub b: u64,
    pub group_order: u64,
    pub r: u64,
    pub cofactor: u64,
    pub generator: FastPoint,
    pub koblitz: Option<KoblitzCurve>,
}

impl BinaryInstance {
    /// Both points with abscissa `x`, or none.
    pub fn points_with_x(&self, x: u64) -> Vec<FastPoint> {
        let f = &self.fast.field;
        if x == 0 {
            let y = f.sqr_k(self.b, self.n - 1);
            return vec![FastPoint::affine(0, y)];
        }
        let inv = f.inv(x);
        let c = x ^ self.a ^ f.mul(self.b, f.sqr(inv));
        match self.artin_schreier.solve(c) {
            Some(u) => {
                let y = f.mul(x, u);
                let p = FastPoint::affine(x, y);
                vec![p, self.fast.neg(p)]
            }
            None => Vec::new(),
        }
    }

    pub fn describe(&self) -> serde_json::Value {
        serde_json::json!({
            "name": self.name,
            "field": {"kind": "binary", "degree": self.n, "polynomial_low_terms": self.irreducible.low_terms},
            "a": self.a, "b": format!("0x{:x}", self.b),
            "group_order": self.group_order, "subgroup_order": self.r, "cofactor": self.cofactor,
            "generator": {"x": format!("0x{:x}", self.generator.x), "y": format!("0x{:x}", self.generator.y)},
            "koblitz": self.koblitz.is_some(),
        })
    }
}

/// A counted wrapper around [`FastCurve`].
pub struct BinaryGroup<'a>(pub &'a FastCurve);

impl CountedGroup for BinaryGroup<'_> {
    type Elt = FastPoint;
    fn identity(&self) -> FastPoint {
        FastPoint::INFINITY
    }
    fn is_identity(&self, p: &FastPoint) -> bool {
        p.infinity
    }
    fn add(&self, ops: &mut GroupOps, p: FastPoint, q: FastPoint) -> FastPoint {
        ops.adds += 1;
        self.0.add(p, q)
    }
    fn double(&self, ops: &mut GroupOps, p: FastPoint) -> FastPoint {
        ops.doubles += 1;
        self.0.double(p)
    }
    fn neg(&self, p: FastPoint) -> FastPoint {
        self.0.neg(p)
    }
    fn key(&self, p: &FastPoint) -> u64 {
        p.pack()
    }
}

/// Trial-division factorisation of a word.
fn factorise_u64(mut v: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut d = 2u64;
    while d * d <= v {
        let mut e = 0;
        while v % d == 0 {
            v /= d;
            e += 1;
        }
        if e > 0 {
            out.push((d, e));
        }
        d += if d == 2 { 1 } else { 2 };
    }
    if v > 1 {
        out.push((v, 1));
    }
    out
}

/// `#E` for `y² + xy = x³ + ax² + b` over `F_{2^n}`: one point at
/// infinity, one at `x = 0`, and two above every non-zero `x` with
/// `Tr(x + a + b/x²) = 0`, found with a batched inversion.
pub fn binary_point_count(gf: &Gf2, ash: &ArtinSchreier, a: u64, b: u64) -> u64 {
    let n = gf.n;
    let total = 1u64 << n;
    let chunk = 4096usize;
    let mut count = 2u64;
    let mut xs: Vec<u64> = Vec::with_capacity(chunk);
    let mut scratch = Vec::with_capacity(chunk);
    let mut x = 1u64;
    let a_tr = ash.trace(a);
    while x < total {
        xs.clear();
        while x < total && xs.len() < chunk {
            xs.push(x);
            x += 1;
        }
        let originals = xs.clone();
        gf.batch_inv(&mut xs, &mut scratch);
        for (&xo, &inv) in originals.iter().zip(&xs) {
            let c = xo ^ gf.mul(b, gf.sqr(inv));
            if ash.trace(c) ^ a_tr == 0 {
                count += 2;
            }
        }
    }
    count
}

/// **A random binary curve** over `F_{2^n}` with a large prime-order
/// subgroup (`r² ∤ #E`, cofactor at most `max_cofactor`), found by
/// counting points of random `(a, b)` until one qualifies.  The search
/// is deterministic in the seed.
pub fn random_binary_instance(n: u32, seed: u64, max_cofactor: u64) -> Option<BinaryInstance> {
    let irreducible = find_irreducible_sparse(n)?;
    let gf = Gf2::new(&irreducible);
    let ash = ArtinSchreier::new(&gf);
    let mut rng = StdRng::seed_from_u64(seed ^ ((n as u64) << 48) ^ 0xC2A2);
    for _attempt in 0..4096 {
        let a = rng.gen_range(0..2u64);
        let b = (rng.gen::<u64>() & gf.mask) | 2;
        if b == 0 || b == 1 {
            continue;
        }
        let order = binary_point_count(&gf, &ash, a, b);
        let factors = factorise_u64(order);
        let Some(&(r, e)) = factors.last() else {
            continue;
        };
        if e != 1 || r < 5 {
            continue;
        }
        let h = order / r;
        if h > max_cofactor {
            continue;
        }
        let curve = BinaryCurve {
            m: n,
            irreducible: irreducible.clone(),
            a: gf.to_element(a),
            b: gf.to_element(b),
            generator: BinaryPoint::Infinity,
            order: BigUint::from(r),
            cofactor: BigUint::from(h),
        };
        let fast = FastCurve::new(&curve)?;
        let inst = BinaryInstance {
            name: format!("random-binary-n{n}-b{b:x}"),
            n,
            irreducible: irreducible.clone(),
            fast,
            gf: gf.clone(),
            artin_schreier: ash.clone(),
            a,
            b,
            group_order: order,
            r,
            cofactor: h,
            generator: FastPoint::INFINITY,
            koblitz: None,
        };
        let group = BinaryGroup(&inst.fast);
        let mut ops = GroupOps::default();
        for _ in 0..4096 {
            let x = rng.gen::<u64>() & gf.mask;
            let Some(&p) = inst.points_with_x(x).first() else {
                continue;
            };
            let g = group.mul(&mut ops, p, h);
            if g.infinity {
                continue;
            }
            if !group.mul(&mut ops, g, r).infinity {
                continue;
            }
            return Some(BinaryInstance {
                generator: g,
                ..inst
            });
        }
    }
    None
}

/// A Koblitz instance, `K_a / F_{2^n}`, from the repository's
/// constructor; `None` when its largest prime factor is too small.
pub fn koblitz_instance(a: u8, n: u32) -> Option<BinaryInstance> {
    let kc = KoblitzCurve::new(a, n)?;
    let fast = FastCurve::new(&kc.curve)?;
    let gf = Gf2::new(&kc.curve.irreducible);
    let ash = ArtinSchreier::new(&gf);
    let generator = fast.lift(kc.generator());
    Some(BinaryInstance {
        name: kc.label(),
        n,
        irreducible: kc.curve.irreducible.clone(),
        a: u64::from(a),
        b: 1,
        group_order: kc.group_order.to_u64()?,
        r: kc.subgroup_order.to_u64()?,
        cofactor: kc.cofactor.to_u64()?,
        generator,
        fast,
        gf,
        artin_schreier: ash,
        koblitz: Some(kc),
    })
}

// ── Factor bases ───────────────────────────────────────────────────

/// How factor-base points map to relation columns.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum ColumnFold {
    /// One column per abscissa; `−P` carries coefficient `−1`.
    Abscissa,
    /// One column per signed Frobenius orbit; `(−1)^s π^k(rep)` carries
    /// `(−1)^s λ^k`.  Koblitz only.
    SignedFrobeniusOrbit,
}

/// A materialised factor base with its column map.  Both `P` and `−P`
/// are entries, so an oracle can name a signed point by index and a
/// relation row is a sum of `coef_of` over the summands.
pub struct FactorBase<E> {
    pub points: Vec<E>,
    pub neg_index: Vec<usize>,
    pub col_of: Vec<usize>,
    pub coef_of: Vec<u64>,
    pub columns: usize,
    pub abscissae: usize,
    pub dimension: Option<u32>,
    pub description: String,
    /// Every abscissa of the base, sorted.
    pub abscissa_list: Vec<u64>,
    /// Abscissa → indices of the points above it.
    pub x_index: FastMap<Vec<usize>>,
    /// Packed point → index.
    pub point_index: FastMap<usize>,
    pub cost: PhaseCost,
    /// Subspace basis (single words) when the base is a linear subspace.
    pub subspace_basis: Option<Vec<u64>>,
}

impl<E: Copy> FactorBase<E> {
    fn empty(description: String) -> Self {
        Self {
            points: Vec::new(),
            neg_index: Vec::new(),
            col_of: Vec::new(),
            coef_of: Vec::new(),
            columns: 0,
            abscissae: 0,
            dimension: None,
            description,
            abscissa_list: Vec::new(),
            x_index: fast_map(0),
            point_index: fast_map(0),
            cost: PhaseCost::default(),
            subspace_basis: None,
        }
    }
}

/// The prime-field base: the `size` smallest abscissae that are on the
/// curve, both signs, one column per abscissa.
pub fn prime_factor_base(inst: &PrimeInstance, size: usize) -> FactorBase<PrimePoint> {
    let start = Instant::now();
    let curve = &inst.curve;
    let mut fb = FactorBase::empty(format!("smallest {size} abscissae"));
    let mut x = 1u64;
    let mut scanned = 0u64;
    let mut sqrts = 0u64;
    while fb.abscissae < size && x < curve.p {
        scanned += 1;
        let rhs = curve.rhs(x);
        if curve.legendre(rhs) >= 0 {
            sqrts += 1;
            if let Some(y) = curve.sqrt(rhs) {
                let y = if y == 0 { 0 } else { y.min(curve.p - y) };
                let p = PrimePoint::affine(x, y);
                let q = curve.neg(p);
                let col = fb.abscissae;
                let i = fb.points.len();
                fb.points.push(p);
                fb.col_of.push(col);
                fb.coef_of.push(1);
                if q != p {
                    fb.points.push(q);
                    fb.col_of.push(col);
                    fb.coef_of.push(inst.r - 1);
                    fb.neg_index.push(i + 1);
                    fb.neg_index.push(i);
                    fb.x_index.insert(x, vec![i, i + 1]);
                } else {
                    fb.neg_index.push(i);
                    fb.x_index.insert(x, vec![i]);
                }
                fb.abscissae += 1;
            }
        }
        x += 1;
    }
    fb.columns = fb.abscissae;
    for (i, p) in fb.points.iter().enumerate() {
        fb.point_index.insert(curve.key(p), i);
    }
    fb.abscissa_list = fb.x_index.keys().copied().collect();
    fb.abscissa_list.sort_unstable();
    fb.cost.wall_ns = start.elapsed().as_nanos() as u64;
    fb.cost.count("abscissae_scanned", scanned);
    fb.cost.count("legendre_symbols", scanned);
    fb.cost.count("sqrt_solves", sqrts);
    fb
}

/// A binary base over the span of `basis` (single words), every point
/// above every non-zero element, columns by abscissa.
pub fn binary_subspace_factor_base(inst: &BinaryInstance, basis: &[u64]) -> FactorBase<FastPoint> {
    let start = Instant::now();
    let l = basis.len() as u32;
    let mut fb = FactorBase::empty(format!("subspace of dimension {l}"));
    let mut solves = 0u64;
    for idx in 1..(1u64 << l) {
        let mut x = 0u64;
        for (j, &e) in basis.iter().enumerate() {
            if (idx >> j) & 1 == 1 {
                x ^= e;
            }
        }
        solves += 1;
        let pts = inst.points_with_x(x);
        if pts.is_empty() {
            continue;
        }
        let col = fb.abscissae;
        let i = fb.points.len();
        let mut idxs = Vec::new();
        for (k, &p) in pts.iter().enumerate() {
            fb.points.push(p);
            fb.col_of.push(col);
            fb.coef_of.push(if k == 0 { 1 } else { inst.r - 1 });
            idxs.push(i + k);
        }
        if pts.len() == 2 {
            fb.neg_index.push(i + 1);
            fb.neg_index.push(i);
        } else {
            fb.neg_index.push(i);
        }
        fb.x_index.insert(x, idxs);
        fb.abscissae += 1;
    }
    fb.columns = fb.abscissae;
    fb.dimension = Some(l);
    fb.subspace_basis = Some(basis.to_vec());
    for (i, p) in fb.points.iter().enumerate() {
        fb.point_index.insert(p.pack(), i);
    }
    fb.abscissa_list = fb.x_index.keys().copied().collect();
    fb.abscissa_list.sort_unstable();
    fb.cost.wall_ns = start.elapsed().as_nanos() as u64;
    fb.cost.count("abscissae_scanned", (1u64 << l) - 1);
    fb.cost.count("as_solves", solves);
    fb
}

/// A base from the repository's Frobenius-invariant constructor, with
/// the chosen column fold.  The construction's own cost (one
/// Artin–Schreier solve per abscissa, one Frobenius map per point to
/// walk the orbits) is counted from the base's size, since the
/// constructor does not carry a ledger; its wall time is measured.
pub fn koblitz_factor_base(
    inst: &BinaryInstance,
    fb: &FrobeniusFactorBase,
    fold: ColumnFold,
    description: String,
) -> Option<FactorBase<FastPoint>> {
    let start = Instant::now();
    let kc = inst.koblitz.as_ref()?;
    let lambda = kc.lambda.to_u64()?;
    let r = inst.r;
    let mut out = FactorBase::empty(description);
    let mut abscissa_col: FastMap<usize> = fast_map(fb.points.len());
    let mut key_to_index: FastMap<usize> = fast_map(fb.points.len());
    for (i, p) in fb.points.iter().enumerate() {
        let fp = inst.fast.lift(p);
        key_to_index.insert(fp.pack(), i);
        out.points.push(fp);
        let (col, coef) = match fold {
            ColumnFold::SignedFrobeniusOrbit => {
                let (orbit, k, negated) = fb.signed_orbit_of[i];
                let c = powmod(lambda, k as u64, r);
                (orbit, if negated { submod(0, c, r) } else { c })
            }
            ColumnFold::Abscissa => {
                let next = abscissa_col.len();
                let col = *abscissa_col.entry(fp.x).or_insert(next);
                // The first point seen above an abscissa is the
                // positive representative; the second is its negative.
                let negated = out.x_index.contains_key(&fp.x);
                (col, if negated { r - 1 } else { 1 })
            }
        };
        out.col_of.push(col);
        out.coef_of.push(coef);
        out.x_index.entry(fp.x).or_default().push(i);
    }
    out.neg_index = out
        .points
        .iter()
        .map(|p| *key_to_index.get(&inst.fast.neg(*p).pack()).expect("base closed under negation"))
        .collect();
    out.point_index = key_to_index;
    out.abscissae = out.x_index.len();
    out.abscissa_list = out.x_index.keys().copied().collect();
    out.abscissa_list.sort_unstable();
    out.columns = match fold {
        ColumnFold::SignedFrobeniusOrbit => fb.signed_orbits.len(),
        ColumnFold::Abscissa => abscissa_col.len(),
    };
    if !fb.uses_ambient_basis() {
        out.dimension = Some(fb.ell);
        out.subspace_basis = Some(
            fb.subspace_basis
                .iter()
                .map(|e| inst.gf.from_element(e))
                .collect(),
        );
    }
    out.cost.wall_ns = start.elapsed().as_nanos() as u64;
    out.cost.count("abscissae_scanned", fb.subspace.len() as u64);
    out.cost.count("as_solves", fb.subspace.len() as u64);
    out.cost.count("frobenius_maps", fb.points.len() as u64);
    out.cost.count("derived_counts", 1);
    Some(out)
}

/// **Can `m` base points sum into the order-`r` subgroup at all?**
///
/// A target `R ∈ ⟨G⟩` has `[r]R = O`, so `R = Σ P_i` forces
/// `Σ [r]P_i = O`: every summand contributes its class in the cofactor
/// group and the classes must cancel.  With `h ≤ 8` there are at most
/// eight classes and the `m`-fold sumset is a handful of additions; the
/// Koblitz regime uses the repository's orbit-closed version of the same
/// test.  Necessary, not sufficient, for a relation to exist.
pub fn summands_admissible<G: CountedGroup>(
    g: &G,
    ops: &mut GroupOps,
    points: &[G::Elt],
    r: u64,
    m: u32,
) -> bool {
    let mut classes: Vec<G::Elt> = Vec::new();
    let mut seen: FastMap<()> = fast_map(64);
    for &p in points {
        let c = g.mul(ops, p, r);
        if seen.insert(g.key(&c), ()).is_none() {
            classes.push(c);
        }
        if classes.len() > 4096 {
            // Too many classes to enumerate sumsets of; assume admissible.
            return true;
        }
    }
    let mut layer: Vec<G::Elt> = vec![g.identity()];
    for _ in 0..m {
        let mut next: Vec<G::Elt> = Vec::new();
        let mut next_seen: FastMap<()> = fast_map(64);
        for &a in &layer {
            for &c in &classes {
                let s = g.add(ops, a, c);
                if next_seen.insert(g.key(&s), ()).is_none() {
                    next.push(s);
                }
            }
        }
        layer = next;
    }
    layer.iter().any(|p| g.is_identity(p))
}

/// **How many of `samples` random subgroup targets decompose** over
/// the base with `m` summands, by the meet-in-the-middle oracle on a
/// built table.  A direct check that relations exist, which the
/// cofactor-class test cannot give: a base that is a subfield's whole
/// point group passes the class test and yields nothing.  Not charged
/// to any phase.
pub fn census_hits(
    inst: &BinaryInstance,
    fb: &FactorBase<FastPoint>,
    table: &PairTable,
    m: u32,
    samples: usize,
    seed: u64,
) -> usize {
    let g = BinaryGroup(&inst.fast);
    let mut rng = StdRng::seed_from_u64(seed ^ 0xCE_A5_05);
    let mut ops = GroupOps::default();
    let mut ctr = OracleCounters::default();
    (0..samples)
        .filter(|_| {
            let q = g.mul(&mut ops, inst.generator, rng.gen_range(1..inst.r));
            decompose_mitm(&g, fb, table, m, &mut ops, &mut ctr, q).is_some()
        })
        .count()
}

// ── Oracles ────────────────────────────────────────────────────────

/// A table of all pair sums of a base, keyed by the sum's abscissa so
/// `S` and `−S` share one entry; the probe checks which it found.
pub struct PairTable {
    map: FastMap<(u32, u32)>,
    pub entries: u64,
    pub build_ops: GroupOps,
    pub build_wall_ns: u64,
}

impl PairTable {
    pub fn build<G: CountedGroup>(g: &G, fb: &FactorBase<G::Elt>) -> Self {
        let start = Instant::now();
        let n = fb.points.len();
        let mut ops = GroupOps::default();
        let mut map: FastMap<(u32, u32)> = fast_map(n * (n + 1) / 4 + 1);
        for i in 0..n {
            for j in i..n {
                let s = g.add(&mut ops, fb.points[i], fb.points[j]);
                if g.is_identity(&s) {
                    continue;
                }
                map.entry(g.key(&s) >> 1).or_insert((i as u32, j as u32));
            }
        }
        Self {
            entries: map.len() as u64,
            map,
            build_ops: ops,
            build_wall_ns: start.elapsed().as_nanos() as u64,
        }
    }

    /// Indices `(i, j)` with `P_i + P_j = target`, if the table holds a
    /// pair summing to `±target`; one probe, at most one addition.
    pub fn probe<G: CountedGroup>(
        &self,
        g: &G,
        fb: &FactorBase<G::Elt>,
        ops: &mut GroupOps,
        lookups: &mut u64,
        target: G::Elt,
    ) -> Option<(usize, usize)> {
        *lookups += 1;
        let &(i, j) = self.map.get(&(g.key(&target) >> 1))?;
        let (i, j) = (i as usize, j as usize);
        let s = g.add(ops, fb.points[i], fb.points[j]);
        if s == target {
            Some((i, j))
        } else if s == g.neg(target) {
            Some((fb.neg_index[i], fb.neg_index[j]))
        } else {
            None
        }
    }

    pub fn calibrate_lookup(&self, samples: u64) -> f64 {
        let mut rng = StdRng::seed_from_u64(7);
        let keys: Vec<u64> = (0..1024).map(|_| rng.gen::<u64>() >> 1).collect();
        let start = Instant::now();
        let mut hits = 0u64;
        for i in 0..samples {
            if self.map.contains_key(&keys[(i % 1024) as usize]) {
                hits += 1;
            }
        }
        let ns = start.elapsed().as_nanos() as f64 / samples as f64;
        std::hint::black_box(hits);
        ns
    }
}

/// Which decomposition oracle a variant uses.
pub enum Oracle<'a> {
    /// Semaev `S₃` as a quadratic in the second abscissa: one square
    /// root per base abscissa per target (`m = 2`, prime).
    SemaevS3Roots,
    /// `R − P_i ∈ F` for every signed base point (`m = 2`).
    Subtract,
    /// Meet in the middle: a pair table probed once per target (`m = 2`)
    /// or once per `R − P_i` (`m = 3`).
    Mitm { table: &'a PairTable, m: u32 },
    /// Semaev `S₄` pairs-and-solve over a subspace (`m = 3`, binary).
    SemaevS4 { oracle: &'a SubspaceOracle },
}

impl Oracle<'_> {
    pub fn name(&self) -> String {
        match self {
            Oracle::SemaevS3Roots => "semaev_s3_roots".into(),
            Oracle::Subtract => "direct_subtraction".into(),
            Oracle::Mitm { m, .. } => format!("meet_in_the_middle_m{m}"),
            Oracle::SemaevS4 { .. } => "semaev_s4_pairs_and_solve".into(),
        }
    }
    pub fn summands(&self) -> u32 {
        match self {
            Oracle::SemaevS3Roots | Oracle::Subtract => 2,
            Oracle::Mitm { m, .. } => *m,
            Oracle::SemaevS4 { .. } => 3,
        }
    }
}

/// Native counters of the oracle, per phase.
#[derive(Clone, Copy, Debug, Default)]
pub struct OracleCounters {
    pub lookups: u64,
    pub sqrt_solves: u64,
    pub inversions: u64,
    pub s4_pairs: u64,
    pub s4_hits: u64,
    pub lift_failures: u64,
}

/// Resolve the signs of a witness of abscissae: which points above them
/// sum to the target.  Costs at most `2·2^m` additions.
fn lift_abscissae<G: CountedGroup>(
    g: &G,
    fb: &FactorBase<G::Elt>,
    ops: &mut GroupOps,
    xs: &[u64],
    target: G::Elt,
) -> Option<Vec<usize>> {
    let mut choices: Vec<&Vec<usize>> = Vec::with_capacity(xs.len());
    for x in xs {
        choices.push(fb.x_index.get(x)?);
    }
    let mut current = vec![0usize; xs.len()];
    loop {
        let mut acc = g.identity();
        for (k, &c) in current.iter().enumerate() {
            acc = g.add(ops, acc, fb.points[choices[k][c]]);
        }
        if acc == target {
            return Some(current.iter().enumerate().map(|(k, &c)| choices[k][c]).collect());
        }
        // Next combination.
        let mut k = 0;
        loop {
            if k == current.len() {
                return None;
            }
            current[k] += 1;
            if current[k] < choices[k].len() {
                break;
            }
            current[k] = 0;
            k += 1;
        }
    }
}

/// Run one oracle on one target: indices of a decomposition, or `None`.
fn decompose_prime(
    curve: &PrimeCurve,
    fb: &FactorBase<PrimePoint>,
    oracle: &Oracle,
    ops: &mut GroupOps,
    ctr: &mut OracleCounters,
    target: PrimePoint,
) -> Option<Vec<usize>> {
    match oracle {
        Oracle::SemaevS3Roots => {
            let p = curve.p;
            let xr = target.x;
            for &xi in &fb.abscissa_list {
                let (qa, qb, qc) = semaev_s3_quadratic(curve, xr, xi);
                let mut candidates: [u64; 2] = [u64::MAX; 2];
                if qa == 0 {
                    if qb == 0 {
                        continue;
                    }
                    ctr.inversions += 1;
                    let Some(inv) = invmod(qb, p) else {
                        continue;
                    };
                    candidates[0] = mulmod(submod(0, qc, p), inv, p);
                } else {
                    let disc = submod(
                        mulmod(qb, qb, p),
                        mulmod(4, mulmod(qa, qc, p), p),
                        p,
                    );
                    ctr.sqrt_solves += 1;
                    let Some(s) = curve.sqrt(disc) else {
                        continue;
                    };
                    ctr.inversions += 1;
                    let Some(inv) = invmod(mulmod(2, qa, p), p) else {
                        continue;
                    };
                    candidates[0] = mulmod(addmod(submod(0, qb, p), s, p), inv, p);
                    candidates[1] = mulmod(submod(submod(0, qb, p), s, p), inv, p);
                }
                for xj in candidates {
                    if xj == u64::MAX {
                        continue;
                    }
                    if !fb.x_index.contains_key(&xj) {
                        continue;
                    }
                    ctr.lookups += 1;
                    if let Some(w) = lift_abscissae(curve, fb, ops, &[xi, xj], target) {
                        return Some(w);
                    }
                    ctr.lift_failures += 1;
                }
            }
            None
        }
        Oracle::Subtract => {
            for i in 0..fb.points.len() {
                let s = curve.add(ops, target, curve.neg(fb.points[i]));
                ctr.lookups += 1;
                if let Some(&j) = fb.point_index.get(&curve.key(&s)) {
                    return Some(vec![i, j]);
                }
            }
            None
        }
        Oracle::Mitm { table, m } => decompose_mitm(curve, fb, table, *m, ops, ctr, target),
        Oracle::SemaevS4 { .. } => None,
    }
}

/// Meet in the middle on a built pair table: one probe for `m = 2`,
/// one `R − P_i` and probe per signed base point for `m = 3`.
pub fn decompose_mitm<G: CountedGroup>(
    g: &G,
    fb: &FactorBase<G::Elt>,
    table: &PairTable,
    m: u32,
    ops: &mut GroupOps,
    ctr: &mut OracleCounters,
    target: G::Elt,
) -> Option<Vec<usize>> {
    match m {
        2 => table
            .probe(g, fb, ops, &mut ctr.lookups, target)
            .map(|(i, j)| vec![i, j]),
        3 => {
            for i in 0..fb.points.len() {
                let s = g.add(ops, target, g.neg(fb.points[i]));
                if g.is_identity(&s) {
                    continue;
                }
                if let Some((j, k)) = table.probe(g, fb, ops, &mut ctr.lookups, s) {
                    return Some(vec![i, j, k]);
                }
            }
            None
        }
        _ => None,
    }
}

fn decompose_binary(
    inst: &BinaryInstance,
    fb: &FactorBase<FastPoint>,
    oracle: &Oracle,
    ops: &mut GroupOps,
    ctr: &mut OracleCounters,
    target: FastPoint,
) -> Option<Vec<usize>> {
    let g = BinaryGroup(&inst.fast);
    match oracle {
        Oracle::Mitm { table, m } => decompose_mitm(&g, fb, table, *m, ops, ctr, target),
        Oracle::Subtract => {
            for i in 0..fb.points.len() {
                let s = g.add(ops, target, g.neg(fb.points[i]));
                ctr.lookups += 1;
                if let Some(&j) = fb.point_index.get(&s.pack()) {
                    return Some(vec![i, j]);
                }
            }
            None
        }
        Oracle::SemaevS4 { oracle } => {
            let (witness, pairs) = oracle.decompose(target.x, &inst.gf);
            ctr.s4_pairs += pairs;
            let xs = witness?;
            ctr.s4_hits += 1;
            match lift_abscissae(&g, fb, ops, &xs, target) {
                Some(w) => Some(w),
                None => {
                    ctr.lift_failures += 1;
                    None
                }
            }
        }
        Oracle::SemaevS3Roots => None,
    }
}

// ── The relation loop ──────────────────────────────────────────────

/// What the shared loop returns.
pub struct PipelineOutcome {
    pub relations: PhaseCost,
    pub linear_algebra: PhaseCost,
    pub verify: PhaseCost,
    pub trials: u64,
    pub relations_found: u64,
    pub independent: u64,
    pub dependent: u64,
    pub rank: u64,
    pub recovered: Option<u64>,
    pub verified: bool,
    pub counters: OracleCounters,
}

/// Draw `R = [a]G + [b]Q`, ask the oracle, turn a decomposition into a
/// row, and stop the moment the elimination pins the logarithm.
#[allow(clippy::too_many_arguments)]
fn collect_and_solve<G: CountedGroup>(
    g: &G,
    generator: G::Elt,
    target: G::Elt,
    r: u64,
    h: u64,
    fb: &FactorBase<G::Elt>,
    seed: u64,
    max_trials: u64,
    mut oracle: impl FnMut(&mut GroupOps, &mut OracleCounters, G::Elt) -> Option<Vec<usize>>,
) -> PipelineOutcome {
    let mut rng = StdRng::seed_from_u64(seed ^ 0x5245_4C41_5449_4F4E);
    let mut rel = PhaseCost::default();
    let mut la = PhaseCost::default();
    let mut ver = PhaseCost::default();
    let mut ctr = OracleCounters::default();
    let cols = fb.columns + 1;
    let d_col = fb.columns;
    let mut gauss = IncrementalGauss::new(cols, r);
    let mut trials = 0u64;
    let mut found = 0u64;
    let mut recovered = None;
    let mut verified = false;
    let rel_start = Instant::now();
    let mut la_ns = 0u64;
    let h_mod = h % r;

    while trials < max_trials {
        trials += 1;
        let a = rng.gen_range(1..r);
        let b = rng.gen_range(1..r);
        let ag = g.mul(&mut rel.group_ops, generator, a);
        let bq = g.mul(&mut rel.group_ops, target, b);
        let point = g.add(&mut rel.group_ops, ag, bq);
        if g.is_identity(&point) {
            rel.count("direct_relations_skipped", 1);
            continue;
        }
        let Some(summands) = oracle(&mut rel.group_ops, &mut ctr, point) else {
            continue;
        };
        found += 1;
        let la_start = Instant::now();
        let mut row = vec![0u64; cols];
        for &i in &summands {
            let c = fb.col_of[i];
            row[c] = addmod(row[c], fb.coef_of[i], r);
        }
        // h·a + h·b·d = Σ coef·x  ⇒  Σ coef·x − h·b·d = h·a
        row[d_col] = submod(0, mulmod(h_mod, b, r), r);
        let rhs = mulmod(h_mod, a, r);
        match gauss.add_row(row, rhs) {
            RowStatus::Inconsistent => {
                la_ns += la_start.elapsed().as_nanos() as u64;
                la.count("inconsistent", 1);
                break;
            }
            RowStatus::Dependent => {}
            RowStatus::Independent => {}
        }
        if let Some(d) = gauss.pinned(d_col) {
            la_ns += la_start.elapsed().as_nanos() as u64;
            let v_start = Instant::now();
            let check = g.mul(&mut ver.group_ops, generator, d);
            ver.wall_ns = v_start.elapsed().as_nanos() as u64;
            if check == target {
                recovered = Some(d);
                verified = true;
                break;
            }
            ver.count("verification_failures", 1);
            break;
        }
        la_ns += la_start.elapsed().as_nanos() as u64;
    }
    rel.wall_ns = rel_start.elapsed().as_nanos() as u64 - la_ns - ver.wall_ns;
    rel.count("trials", trials);
    rel.count("relations", found);
    rel.count("lookups", ctr.lookups);
    rel.count("sqrt_solves", ctr.sqrt_solves);
    rel.count("inversions", ctr.inversions);
    rel.count("s4_pairs", ctr.s4_pairs);
    rel.count("s4_hits", ctr.s4_hits);
    rel.count("lift_failures", ctr.lift_failures);
    la.wall_ns = la_ns;
    la.count("row_ops", gauss.row_ops);
    la.count("rows", found);
    la.count("columns", cols as u64);
    la.count("rank", gauss.rank() as u64);
    PipelineOutcome {
        relations: rel,
        linear_algebra: la,
        verify: ver,
        trials,
        relations_found: found,
        independent: gauss.rank() as u64,
        dependent: gauss.dependent,
        rank: gauss.rank() as u64,
        recovered,
        verified,
        counters: ctr,
    }
}

// ── Calibration ────────────────────────────────────────────────────

fn calibrate_group<G: CountedGroup>(g: &G, points: &[G::Elt], calib: &mut Calibration) {
    let mut ops = GroupOps::default();
    let n = points.len();
    let samples = 200_000u64;
    let mut acc = points[0];
    let start = Instant::now();
    for i in 0..samples {
        acc = g.add(&mut ops, acc, points[(i as usize * 7 + 1) % n]);
        if g.is_identity(&acc) {
            acc = points[(i as usize + 3) % n];
        }
    }
    calib.ns_per_add = start.elapsed().as_nanos() as f64 / samples as f64;
    std::hint::black_box(acc);
    let mut acc = points[1];
    let start = Instant::now();
    for i in 0..samples {
        acc = g.double(&mut ops, acc);
        if g.is_identity(&acc) {
            acc = points[(i as usize + 5) % n];
        }
    }
    calib.ns_per_double = start.elapsed().as_nanos() as f64 / samples as f64;
    std::hint::black_box(acc);
}

fn calibrate_row_ops(modulus: u64, calib: &mut Calibration) {
    let mut rng = StdRng::seed_from_u64(11);
    let cols = 2048usize;
    let pivot: Vec<u64> = (0..cols).map(|_| rng.gen_range(0..modulus)).collect();
    let mut row: Vec<u64> = (0..cols).map(|_| rng.gen_range(0..modulus)).collect();
    let reps = 500u64;
    let start = Instant::now();
    for k in 0..reps {
        let f = (k + 3) % modulus;
        for j in 0..cols {
            row[j] = submod(row[j], mulmod(f, pivot[j], modulus), modulus);
        }
    }
    calib.ns_per_row_op = start.elapsed().as_nanos() as f64 / (reps * cols as u64) as f64;
    std::hint::black_box(&row);
}

/// A 64-bit word XOR of one row into another, the unit of the Boolean
/// Macaulay elimination.
pub fn calibrate_word_xor() -> f64 {
    let words = 4096usize;
    let mut a: Vec<u64> = (0..words as u64).map(mix).collect();
    let b: Vec<u64> = (0..words as u64).map(|i| mix(i + 77)).collect();
    let reps = 2000u64;
    let start = Instant::now();
    for k in 0..reps {
        let shift = (k % 7) as usize;
        for j in 0..words {
            a[j] ^= b[(j + shift) % words];
        }
    }
    let ns = start.elapsed().as_nanos() as f64 / (reps * words as u64) as f64;
    std::hint::black_box(&a);
    ns
}

// ── Pricing a run ──────────────────────────────────────────────────

fn price_phase(phase: &mut PhaseCost, calib: &Calibration) {
    let mut gae = phase.group_ops.gae();
    let conv = |name: &str, ns: Option<f64>| -> f64 {
        match ns {
            Some(ns) => calib.gae(phase.get(name), ns),
            None => 0.0,
        }
    };
    gae += conv("sqrt_solves", calib.ns_per_sqrt);
    gae += conv("as_solves", calib.ns_per_as_solve);
    gae += conv("s4_pairs", calib.ns_per_s4_pair);
    gae += conv("lookups", Some(calib.ns_per_lookup));
    gae += conv("row_ops", Some(calib.ns_per_row_op));
    gae += conv("legendre_symbols", calib.ns_per_legendre);
    gae += conv("inversions", calib.ns_per_inversion);
    gae += conv("frobenius_maps", calib.ns_per_frobenius);
    phase.gae = gae;
}

#[allow(clippy::too_many_arguments)]
fn assemble_variant(
    name: &str,
    oracle_name: String,
    summands: u32,
    fb_cost: &PhaseCost,
    signed_points: u64,
    abscissae: u64,
    columns: u64,
    dimension: Option<u32>,
    table: Option<&PairTable>,
    outcome: PipelineOutcome,
    r: u64,
    space: f64,
    calib: &Calibration,
    rho_s: f64,
    floor_s: f64,
) -> VariantResult {
    let mut fb_phase = fb_cost.clone();
    if let Some(t) = table {
        fb_phase.group_ops.merge(t.build_ops);
        fb_phase.wall_ns += t.build_wall_ns;
        fb_phase.count("pair_table_entries", t.entries);
    }
    price_phase(&mut fb_phase, calib);
    let mut rel = outcome.relations;
    price_phase(&mut rel, calib);
    let mut la = outcome.linear_algebra;
    price_phase(&mut la, calib);
    let mut ver = outcome.verify;
    price_phase(&mut ver, calib);
    let total_gae = fb_phase.gae + rel.gae + la.gae + ver.gae;
    let total_wall = fb_phase.wall_ns + rel.wall_ns + la.wall_ns + ver.wall_ns;
    let sqrt_r = (r as f64).sqrt();
    let s = total_gae / sqrt_r;
    let s_wall = if calib.ns_per_add > 0.0 {
        total_wall as f64 / calib.ns_per_add / sqrt_r
    } else {
        f64::NAN
    };
    let tf = trials_floor(columns, signed_points, summands, space);
    let p_ceiling = decomposition_probability_ceiling(signed_points, summands, space);
    let measured_yield = outcome.relations_found as f64 / outcome.trials.max(1) as f64;
    VariantResult {
        name: name.to_string(),
        oracle: oracle_name,
        summands,
        factor_base: fb_phase,
        signed_points,
        abscissae,
        columns,
        dimension,
        relations: rel,
        linear_algebra: la,
        verify: ver,
        trials: outcome.trials,
        relations_found: outcome.relations_found,
        independent: outcome.independent,
        dependent: outcome.dependent,
        rank: outcome.rank,
        trials_floor: tf,
        trials_over_floor: outcome.trials as f64 / tf,
        decomposition_probability_ceiling: p_ceiling,
        yield_over_ceiling: measured_yield / p_ceiling,
        total_gae,
        total_wall_ns: total_wall,
        s,
        s_wall,
        ratio_to_rho: s / rho_s,
        ratio_to_floor: s / floor_s,
        recovered: outcome.recovered,
        verified: outcome.verified,
    }
}

// ── Configuration ──────────────────────────────────────────────────

/// What to run.
#[derive(Clone, Debug, Serialize)]
pub struct BoundaryConfig {
    pub seed: u64,
    /// Repeats per instance (fresh target each; rho and every variant
    /// see the same target within a repeat).
    pub repeats: usize,
    /// Rho walks per repeat on that target (different walk seeds).  The
    /// reference is cheap and its single-run spread is about half its
    /// mean, so it gets more samples than the variants.
    pub rho_repeats: usize,
    pub prime_bits: Vec<u32>,
    pub char2_degrees: Vec<u32>,
    pub koblitz_degrees: Vec<u32>,
    /// Cap on targets drawn per variant run.
    pub max_trials: u64,
    /// A run also stops after this multiple of the counting floor's
    /// trial count (never fewer than 4096 trials), so a base that cannot
    /// yield relations fails in bounded time instead of at the cap.
    pub trial_floor_multiple: f64,
    /// Largest Koblitz degree at which the un-folded control (one column
    /// per abscissa) is run; it needs `2n` times the relations.
    pub koblitz_no_fold_max_degree: u32,
    /// Largest degree at which the `S₄` pairs-and-solve oracle is run.
    pub s4_max_degree: u32,
    /// Largest base the pair table may hold (signed points).
    pub max_table_points: usize,
    /// Rho step cap as a multiple of the expected count.
    pub rho_cap_multiple: f64,
}

impl Default for BoundaryConfig {
    fn default() -> Self {
        Self {
            seed: 0x1C_B0_0B_DA_7A,
            repeats: 3,
            rho_repeats: 8,
            prime_bits: vec![10, 12, 14, 16, 18, 20, 22, 24],
            char2_degrees: vec![15, 18, 21, 24, 27],
            // The degrees whose K_a has a usable prime-order subgroup
            // (largest prime factor above the cofactor, squarefree);
            // 21, 25, 27, 33 and 35 have none.
            koblitz_degrees: vec![11, 13, 15, 17, 19, 23, 29, 31, 37, 39, 41],
            max_trials: 50_000_000,
            trial_floor_multiple: 64.0,
            koblitz_no_fold_max_degree: 31,
            s4_max_degree: 31,
            max_table_points: 6000,
            rho_cap_multiple: 64.0,
        }
    }
}

impl BoundaryConfig {
    /// A configuration small enough for a unit test.
    pub fn quick() -> Self {
        Self {
            repeats: 1,
            rho_repeats: 2,
            prime_bits: vec![10, 12],
            char2_degrees: vec![12],
            koblitz_degrees: vec![11, 15],
            koblitz_no_fold_max_degree: 15,
            s4_max_degree: 15,
            ..Self::default()
        }
    }
}

fn rho_cap(r: u64, multiple: f64) -> u64 {
    (generic_floor_ops(r as f64, 1.0) * multiple) as u64 + 4096
}

/// Trials a variant may draw: the configured multiple of the counting
/// floor, at least 4096, at most the absolute cap.
fn trial_budget<E>(cfg: &BoundaryConfig, fb: &FactorBase<E>, m: u32, space: f64) -> u64 {
    let floor = trials_floor(fb.columns as u64, fb.points.len() as u64, m, space);
    ((cfg.trial_floor_multiple * floor).ceil() as u64).clamp(4096, cfg.max_trials.max(4096))
}

// ── The prime regime ───────────────────────────────────────────────

/// Run every prime-field variant on one instance.
pub fn run_prime_instance(inst: &PrimeInstance, cfg: &BoundaryConfig) -> RegimeInstance {
    let curve = &inst.curve;
    let r = inst.r;
    let g = inst.generator_point();
    let mut calib = Calibration::default();
    let mut rng = StdRng::seed_from_u64(cfg.seed ^ r);

    // Calibration on random points of the subgroup.
    let mut ops = GroupOps::default();
    let pts: Vec<PrimePoint> = (0..64)
        .map(|_| curve.mul(&mut ops, g, rng.gen_range(1..r)))
        .collect();
    calibrate_group(curve, &pts, &mut calib);
    calibrate_row_ops(r, &mut calib);
    {
        let mut residues = Vec::new();
        while residues.len() < 2048 {
            let v = rng.gen_range(1..curve.p);
            if curve.legendre(v) == 1 {
                residues.push(v);
            }
        }
        let start = Instant::now();
        let mut acc = 0u64;
        for _ in 0..8 {
            for &v in &residues {
                acc ^= curve.sqrt(v).unwrap_or(0);
            }
        }
        calib.ns_per_sqrt = Some(start.elapsed().as_nanos() as f64 / (8 * residues.len()) as f64);
        std::hint::black_box(acc);
        let values: Vec<u64> = (0..2048).map(|_| rng.gen_range(1..curve.p)).collect();
        let start = Instant::now();
        let mut acc = 0i32;
        for _ in 0..8 {
            for &v in &values {
                acc += curve.legendre(v);
            }
        }
        calib.ns_per_legendre = Some(start.elapsed().as_nanos() as f64 / (8 * values.len()) as f64);
        std::hint::black_box(acc);
        let start = Instant::now();
        let mut acc = 0u64;
        for _ in 0..8 {
            for &v in &values {
                acc ^= invmod(v, curve.p).unwrap_or(0);
            }
        }
        calib.ns_per_inversion = Some(start.elapsed().as_nanos() as f64 / (8 * values.len()) as f64);
        std::hint::black_box(acc);
    }
    calib.ns_per_word_xor = Some(calibrate_word_xor());

    let bits = 64 - r.leading_zeros();
    let fb_size = 1usize << bits.div_ceil(3);
    let fb = prime_factor_base(inst, fb_size);
    let table = PairTable::build(curve, &fb);
    calib.ns_per_lookup = table.calibrate_lookup(1_000_000);

    let floor_s = generic_floor_s(2.0);
    let mut out = RegimeInstance {
        regime: "prime".into(),
        curve: serde_json::json!({
            "name": inst.name, "p": curve.p, "a": curve.a, "b": curve.b,
            "group_order": inst.group_order, "subgroup_order": r, "cofactor": inst.cofactor,
            "generator": {"x": inst.generator.0, "y": inst.generator.1},
        }),
        r,
        log2_r: (r as f64).log2(),
        group_order: inst.group_order as f64,
        cofactor: inst.cofactor,
        automorphisms_generic: 2,
        floor_s,
        floor_ops: generic_floor_ops(r as f64, 2.0),
        calibration: calib.clone(),
        rho: Vec::new(),
        rho_s_mean: 0.0,
        rho_verified_all: true,
        variants: Vec::new(),
        seeds: Vec::new(),
        targets: Vec::new(),
    };

    let space = inst.group_order as f64;
    for rep in 0..cfg.repeats {
        let seed = cfg.seed.wrapping_add(rep as u64 * 0x9E37);
        let d = rng.gen_range(1..r);
        let target = curve.mul(&mut ops, g, d);
        out.seeds.push(seed);
        out.targets.push(d);
        let mut rho_s = 0.0;
        for k in 0..cfg.rho_repeats.max(1) {
            let rho = rho_reference(curve, g, target, r, seed ^ (k as u64 * 0x5EED), rho_cap(r, cfg.rho_cap_multiple));
            out.rho_verified_all &= rho.verified && rho.recovered == Some(d);
            rho_s += rho.s / cfg.rho_repeats.max(1) as f64;
            out.rho.push(rho);
        }

        let variants: Vec<(&str, Oracle)> = vec![
            ("semaev_s3_roots_m2", Oracle::SemaevS3Roots),
            ("direct_subtraction_m2", Oracle::Subtract),
            (
                "mitm_m2",
                Oracle::Mitm {
                    table: &table,
                    m: 2,
                },
            ),
            (
                "mitm_m3",
                Oracle::Mitm {
                    table: &table,
                    m: 3,
                },
            ),
        ];
        for (name, oracle) in variants {
            let budget = trial_budget(cfg, &fb, oracle.summands(), space);
            let outcome = collect_and_solve(
                curve,
                g,
                target,
                r,
                inst.cofactor,
                &fb,
                seed,
                budget,
                |ops, ctr, point| decompose_prime(curve, &fb, &oracle, ops, ctr, point),
            );
            let uses_table = matches!(oracle, Oracle::Mitm { .. });
            let mut v = assemble_variant(
                name,
                oracle.name(),
                oracle.summands(),
                &fb.cost,
                fb.points.len() as u64,
                fb.abscissae as u64,
                fb.columns as u64,
                None,
                uses_table.then_some(&table),
                outcome,
                r,
                space,
                &calib,
                rho_s,
                floor_s,
            );
            v.verified &= v.recovered == Some(d);
            out.variants.push(v);
        }
    }
    out.rho_s_mean = out.rho.iter().map(|x| x.s).sum::<f64>() / out.rho.len().max(1) as f64;
    out
}

// ── The binary regimes ─────────────────────────────────────────────

/// The binary-regime conversion factors for one instance: addition,
/// doubling, Artin–Schreier solve, Frobenius map, row operation, word
/// XOR.  The lookup and `S₄` factors are added by the callers that
/// build the table and the oracle.
pub fn calibrate_binary_instance(inst: &BinaryInstance) -> Calibration {
    let mut c = Calibration::default();
    calibrate_binary(inst, &mut c);
    c
}

fn calibrate_binary(inst: &BinaryInstance, calib: &mut Calibration) {
    let g = BinaryGroup(&inst.fast);
    let mut rng = StdRng::seed_from_u64(inst.r ^ 0xB1);
    let mut ops = GroupOps::default();
    let pts: Vec<FastPoint> = (0..64)
        .map(|_| g.mul(&mut ops, inst.generator, rng.gen_range(1..inst.r)))
        .collect();
    calibrate_group(&g, &pts, calib);
    calibrate_row_ops(inst.r, calib);
    let samples = 200_000u64;
    let cs: Vec<u64> = (0..1024).map(|_| rng.gen::<u64>() & inst.gf.mask).collect();
    let start = Instant::now();
    let mut acc = 0u64;
    for i in 0..samples {
        acc ^= inst.artin_schreier.solve(cs[(i % 1024) as usize]).unwrap_or(1);
    }
    calib.ns_per_as_solve = Some(start.elapsed().as_nanos() as f64 / samples as f64);
    std::hint::black_box(acc);
    let start = Instant::now();
    let mut acc = pts[0];
    for _ in 0..samples {
        acc = inst.fast.frobenius_k(acc, 1);
    }
    calib.ns_per_frobenius = Some(start.elapsed().as_nanos() as f64 / samples as f64);
    std::hint::black_box(acc);
    calib.ns_per_word_xor = Some(calibrate_word_xor());
}

/// Nanoseconds per pair of the `S₄` pairs-and-solve loop, measured on
/// random targets of this instance.
pub fn calibrate_s4(inst: &BinaryInstance, oracle: &SubspaceOracle, calib: &mut Calibration) {
    let mut rng = StdRng::seed_from_u64(inst.r ^ 0x54);
    let mut pairs = 0u64;
    let start = Instant::now();
    let mut tries = 0;
    while tries < 3 || (pairs < 200_000 && tries < 64) {
        let xr = rng.gen::<u64>() & inst.gf.mask;
        let (_, p) = oracle.decompose(xr, &inst.gf);
        pairs += p;
        tries += 1;
    }
    calib.ns_per_s4_pair = Some(start.elapsed().as_nanos() as f64 / pairs.max(1) as f64);
}

/// The Koblitz base for degree `n`: the Frobenius-invariant subspace
/// whose dimension is nearest `⌈n/3⌉` (at most 12), or, where none
/// exists, the Frobenius closure of a random seed space sized to about
/// `2^{⌈n/3⌉}` abscissae — capped by the pair-table budget.
pub fn choose_koblitz_base(
    inst: &BinaryInstance,
    seed: u64,
    max_points: usize,
) -> Option<(FrobeniusFactorBase, String)> {
    let kc = inst.koblitz.as_ref()?;
    let n = inst.n;
    let target = n.div_ceil(3).max(3);
    let max_dim = 12u32.min((max_points as f64 / 2.0).log2().floor() as u32 + 1);
    let usable_dims: Vec<u32> = available_subspace_dimensions(n)
        .into_iter()
        .filter(|&d| d >= 2 && d < n && d <= max_dim && (d as i64 - target as i64).abs() <= 2)
        .collect();
    if !usable_dims.is_empty() {
        // Every divisor of x^n − 1 whose degree is a usable dimension,
        // nearest the target first.  A Frobenius-invariant subspace can
        // hold almost no curve abscissae (the dimension-5 subspace of
        // K_0/F_2^15 holds three points), so a candidate is kept only
        // when it carries at least a quarter of the points a random set
        // of its size would.
        let factors = all_factors_of_x_n_minus_1(n);
        let degs: Vec<u32> = factors.iter().map(|&f| 63 - f.leading_zeros()).collect();
        let count = factors.len().min(16);
        let mut candidates: Vec<(i64, u32, Vec<usize>)> = Vec::new();
        for mask in 1u32..(1u32 << count) {
            let idx: Vec<usize> = (0..count).filter(|&i| (mask >> i) & 1 == 1).collect();
            let total: u32 = idx.iter().map(|&i| degs[i]).sum();
            if usable_dims.contains(&total) {
                candidates.push(((total as i64 - target as i64).abs(), total, idx));
            }
        }
        candidates.sort();
        let mut best: Option<(FrobeniusFactorBase, String)> = None;
        for (_, dim, idx) in candidates {
            if let Some(fb) = build_frobenius_factor_base_from_divisor(kc, &idx) {
                // A subspace that is a subfield F_{2^d} carries the whole
                // point group E(F_{2^d}), a subgroup that meets the target
                // subgroup only at O: no relation exists at any cost.
                let is_subfield = n % dim == 0
                    && fb.subspace_basis.iter().all(|v| {
                        v.square_k_times(dim, &kc.curve.irreducible) == *v
                    });
                let healthy = !is_subfield
                    && fb.points.len() >= (1usize << dim) / 4
                    && fb.points.len() <= max_points;
                if healthy {
                    let better = best
                        .as_ref()
                        .is_none_or(|(b, _)| b.points.len() < fb.points.len() && b.ell == dim);
                    if better {
                        best = Some((
                            fb,
                            format!("invariant subspace, divisor {idx:?}, dimension {dim}"),
                        ));
                    }
                    if best.as_ref().is_some_and(|(b, _)| b.ell != dim) {
                        break;
                    }
                }
            }
        }
        if let Some(found) = best {
            return Some(found);
        }
    }
    // Orbit union of a random seed space, sized like the subspaces: about
    // `2^{⌈n/3⌉+1}` abscissae (`n · 2^s`), within the table budget.
    let mut rng = StdRng::seed_from_u64(seed ^ (n as u64));
    let log2_n = 63 - (n as u64).leading_zeros();
    let mut s = (target + 1).saturating_sub(log2_n).clamp(1, 12);
    while s > 1 && (n as usize) * (1usize << s) > max_points {
        s -= 1;
    }
    for _ in 0..64 {
        let basis: Vec<_> = (0..s)
            .map(|_| {
                let bits: Vec<u32> = (0..n).filter(|_| rng.gen::<bool>()).collect();
                crate::binary_ecc::F2mElement::from_bit_positions(&bits, n)
            })
            .collect();
        if let Some(fb) = build_frobenius_union_factor_base(kc, &basis) {
            if fb.points.len() <= max_points && fb.points.len() >= 8 {
                return Some((fb, format!("Frobenius orbit union of a seed space of dimension {s}")));
            }
        }
    }
    None
}

fn binary_regime_shell(inst: &BinaryInstance, regime: &str, automorphisms: u32) -> RegimeInstance {
    RegimeInstance {
        regime: regime.into(),
        curve: inst.describe(),
        r: inst.r,
        log2_r: (inst.r as f64).log2(),
        group_order: inst.group_order as f64,
        cofactor: inst.cofactor,
        automorphisms_generic: automorphisms,
        floor_s: generic_floor_s(automorphisms as f64),
        floor_ops: generic_floor_ops(inst.r as f64, automorphisms as f64),
        calibration: Calibration::default(),
        rho: Vec::new(),
        rho_s_mean: 0.0,
        rho_verified_all: true,
        variants: Vec::new(),
        seeds: Vec::new(),
        targets: Vec::new(),
    }
}

/// Run every generic-binary variant on one instance: `S₄`
/// pairs-and-solve and meet in the middle over the low-order subspace
/// of dimension `⌈n/3⌉`, columns by abscissa.  `None` when no summand
/// count in `{2, 3}` can reach the subgroup from this base.
pub fn run_char2_instance(inst: &BinaryInstance, cfg: &BoundaryConfig) -> Option<RegimeInstance> {
    let g = BinaryGroup(&inst.fast);
    let r = inst.r;
    let mut calib = Calibration::default();
    calibrate_binary(inst, &mut calib);
    let l = inst.n.div_ceil(3);
    let basis: Vec<u64> = (0..l).map(|i| 1u64 << i).collect();
    let mut fb = binary_subspace_factor_base(inst, &basis);
    let mut class_ops = GroupOps::default();
    let admissible: Vec<u32> = [3u32, 2]
        .into_iter()
        .filter(|&m| summands_admissible(&g, &mut class_ops, &fb.points, r, m))
        .collect();
    fb.cost.group_ops.merge(class_ops);
    fb.cost.count("admissibility_scalar_mults", class_ops.scalar_mults);
    let table = PairTable::build(&g, &fb);
    calib.ns_per_lookup = table.calibrate_lookup(1_000_000);
    let census: Vec<(u32, usize)> = admissible
        .iter()
        .map(|&m| (m, census_hits(inst, &fb, &table, m, 64, cfg.seed)))
        .collect();
    let m_used = census.iter().find(|(_, hits)| *hits > 0).map(|(m, _)| *m)?;
    let s4 = (m_used == 3).then(|| SubspaceOracle::new(&basis, inst.b, &inst.gf));
    if let Some(o) = &s4 {
        calibrate_s4(inst, o, &mut calib);
    }

    let mut out = binary_regime_shell(inst, "char2", 2);
    out.calibration = calib.clone();
    out.curve["factor_base"] = serde_json::json!({
        "description": fb.description,
        "dimension": l,
        "signed_points": fb.points.len(),
        "abscissae": fb.abscissae,
        "columns": fb.columns,
        "summands_admissible": admissible,
        "census_hits_of_64": census,
        "summands_used": m_used,
    });
    let space = inst.group_order as f64;
    let mut rng = StdRng::seed_from_u64(cfg.seed ^ r ^ 0xC2);
    let mut ops = GroupOps::default();
    for rep in 0..cfg.repeats {
        let seed = cfg.seed.wrapping_add(rep as u64 * 0x9E37);
        let d = rng.gen_range(1..r);
        let target = g.mul(&mut ops, inst.generator, d);
        out.seeds.push(seed);
        out.targets.push(d);
        let mut rho_s = 0.0;
        for k in 0..cfg.rho_repeats.max(1) {
            let rho = rho_reference(&g, inst.generator, target, r, seed ^ (k as u64 * 0x5EED), rho_cap(r, cfg.rho_cap_multiple));
            out.rho_verified_all &= rho.verified && rho.recovered == Some(d);
            rho_s += rho.s / cfg.rho_repeats.max(1) as f64;
            out.rho.push(rho);
        }

        let mut variants: Vec<(String, Oracle)> = vec![(
            format!("mitm_m{m_used}"),
            Oracle::Mitm {
                table: &table,
                m: m_used,
            },
        )];
        if let Some(o) = &s4 {
            if inst.n <= cfg.s4_max_degree {
                variants.push(("semaev_s4_pairs_and_solve_m3".into(), Oracle::SemaevS4 { oracle: o }));
            }
        }
        for (name, oracle) in variants {
            let budget = trial_budget(cfg, &fb, oracle.summands(), space);
            let outcome = collect_and_solve(
                &g,
                inst.generator,
                target,
                r,
                inst.cofactor,
                &fb,
                seed,
                budget,
                |ops, ctr, point| decompose_binary(inst, &fb, &oracle, ops, ctr, point),
            );
            let uses_table = matches!(oracle, Oracle::Mitm { .. });
            let mut v = assemble_variant(
                &name,
                oracle.name(),
                oracle.summands(),
                &fb.cost,
                fb.points.len() as u64,
                fb.abscissae as u64,
                fb.columns as u64,
                fb.dimension,
                uses_table.then_some(&table),
                outcome,
                r,
                space,
                &calib,
                rho_s,
                out.floor_s,
            );
            v.verified &= v.recovered == Some(d);
            out.variants.push(v);
        }
    }
    out.rho_s_mean = out.rho.iter().map(|x| x.s).sum::<f64>() / out.rho.len().max(1) as f64;
    Some(out)
}

/// Run every Koblitz variant on one instance: meet in the middle over
/// signed-orbit columns, the same without the fold (the `char2` pipeline
/// on the Koblitz curve), and `S₄` pairs-and-solve on the invariant
/// subspace when the base is one.  The reference is the signed-Frobenius
/// rho (`A = 2n`); the plain walk is run too, for the comparison with
/// the other regimes.
pub fn run_koblitz_instance(inst: &BinaryInstance, cfg: &BoundaryConfig) -> Option<RegimeInstance> {
    let kc = inst.koblitz.as_ref()?;
    let g = BinaryGroup(&inst.fast);
    let r = inst.r;
    let n = inst.n;
    let mut calib = Calibration::default();
    calibrate_binary(inst, &mut calib);
    let (mut frob, mut description) = choose_koblitz_base(inst, cfg.seed, cfg.max_table_points)?;
    // Which summand counts actually reach the subgroup from this base?
    // The cofactor-class test is the cheap necessary filter; the census
    // on the built pair table is the real one.  When neither count
    // yields, the two-torsion saturation `B ∪ (B + T)` changes the
    // cofactor classes without adding projected points.
    let admissible_of = |fb: &FrobeniusFactorBase| -> Vec<u32> {
        [3u32, 2]
            .into_iter()
            .filter(|&m| fb.m_can_decompose(kc, m as usize))
            .collect()
    };
    let mut admissible = admissible_of(&frob);
    let mut folded = koblitz_factor_base(inst, &frob, ColumnFold::SignedFrobeniusOrbit, description.clone())?;
    let mut table = PairTable::build(&g, &folded);
    let census_of = |fb: &FactorBase<FastPoint>, table: &PairTable, adm: &[u32]| -> Vec<(u32, usize)> {
        adm.iter()
            .map(|&m| (m, census_hits(inst, fb, table, m, 64, cfg.seed)))
            .collect()
    };
    let mut census = census_of(&folded, &table, &admissible);
    if !census.iter().any(|(_, hits)| *hits > 0) {
        let saturated = saturate_factor_base_two_torsion(kc, &frob)?;
        if saturated.points.len() > cfg.max_table_points {
            return None;
        }
        admissible = admissible_of(&saturated);
        frob = saturated;
        description.push_str(" + two-torsion saturation");
        folded = koblitz_factor_base(inst, &frob, ColumnFold::SignedFrobeniusOrbit, description.clone())?;
        table = PairTable::build(&g, &folded);
        census = census_of(&folded, &table, &admissible);
    }
    let m_used = census.iter().find(|(_, hits)| *hits > 0).map(|(m, _)| *m)?;
    let unfolded = koblitz_factor_base(inst, &frob, ColumnFold::Abscissa, description.clone())?;
    calib.ns_per_lookup = table.calibrate_lookup(1_000_000);
    let s4 = folded
        .subspace_basis
        .as_ref()
        .filter(|_| n <= cfg.s4_max_degree && m_used == 3)
        .map(|basis| SubspaceOracle::new(basis, 1, &inst.gf));
    if let Some(o) = &s4 {
        calibrate_s4(inst, o, &mut calib);
    }

    let mut out = binary_regime_shell(inst, "koblitz", 2 * n);
    out.calibration = calib.clone();
    out.curve["factor_base"] = serde_json::json!({
        "description": description,
        "domain": format!("{:?}", frob.domain),
        "signed_points": folded.points.len(),
        "abscissae": folded.abscissae,
        "signed_orbits": folded.columns,
        "dimension": folded.dimension,
        "lambda": kc.lambda.to_string(),
        "summands_admissible": admissible,
        "census_hits_of_64": census,
        "summands_used": m_used,
    });
    let space = inst.group_order as f64;
    let mut rng = StdRng::seed_from_u64(cfg.seed ^ r ^ 0x4B);
    let mut ops = GroupOps::default();
    let bits = (r as f64).log2();
    for rep in 0..cfg.repeats {
        let seed = cfg.seed.wrapping_add(rep as u64 * 0x9E37);
        let d = rng.gen_range(1..r);
        let target = g.mul(&mut ops, inst.generator, d);
        out.seeds.push(seed);
        out.targets.push(d);

        // Signed-Frobenius rho: the reference on this curve.
        let target_big = inst.fast.lower(target);
        let mut rho_s = 0.0;
        for k in 0..cfg.rho_repeats.max(1) {
            let opts = KoblitzSignedRhoOptions {
                seed: seed ^ (k as u64 * 0x5EED),
                ..Default::default()
            };
            let start = Instant::now();
            let report = koblitz_signed_frobenius_rho_reference(kc, &target_big, &opts, &mut |_| {});
            let wall = start.elapsed().as_nanos() as u64;
            let c = &report.charges;
            // Scalar multiplications are charged at 1.5·log2(r) additions
            // each, the double-and-add average; walk additions are exact.
            let scalar_adds = (c.setup_scalar_multiplications
                + c.candidate_verification_scalar_multiplications) as f64
                * 1.5
                * bits;
            let gops = GroupOps {
                adds: c.walk_group_additions + c.setup_group_additions,
                doubles: 0,
                scalar_mults: c.setup_scalar_multiplications
                    + c.candidate_verification_scalar_multiplications,
            };
            let gae = gops.gae() + scalar_adds;
            let expected = generic_floor_ops(r as f64, 2.0 * n as f64);
            let recovered = report.recovered_log.as_ref().and_then(|v| v.to_u64());
            let signed = RhoResult {
                method: "signed-Frobenius r-adding walk, distinguished points (koblitz_signed_frobenius_rho_reference)".into(),
                automorphisms: 2 * n,
                group_ops: gops,
                steps: report.iterations,
                walks: report.parallel_walks as u64 * report.restarts_attempted as u64,
                distinguished_points: c.collisions + c.failed_collisions,
                gae,
                s: gae / (r as f64).sqrt(),
                s_walk: report.iterations as f64 / (r as f64).sqrt(),
                expected_steps: expected,
                steps_over_expected: report.iterations as f64 / expected,
                wall_ns: wall,
                recovered,
                verified: report.verified && recovered == Some(d),
            };
            out.rho_verified_all &= signed.verified;
            rho_s += signed.s / cfg.rho_repeats.max(1) as f64;
            out.rho.push(signed);
            // The plain walk, for the cross-regime comparison.
            let plain = rho_reference(&g, inst.generator, target, r, seed ^ (k as u64 * 0x5EED), rho_cap(r, cfg.rho_cap_multiple));
            out.rho.push(plain);
        }

        let mut variants: Vec<(String, &FactorBase<FastPoint>, Oracle)> = vec![(
            format!("mitm_m{m_used}_signed_orbit_columns"),
            &folded,
            Oracle::Mitm {
                table: &table,
                m: m_used,
            },
        )];
        if n <= cfg.koblitz_no_fold_max_degree {
            variants.push((
                format!("mitm_m{m_used}_abscissa_columns_control"),
                &unfolded,
                Oracle::Mitm {
                    table: &table,
                    m: m_used,
                },
            ));
        }
        if let Some(o) = &s4 {
            variants.push((
                "semaev_s4_pairs_and_solve_m3_signed_orbit_columns".into(),
                &folded,
                Oracle::SemaevS4 { oracle: o },
            ));
        }
        for (name, fb, oracle) in variants {
            let budget = trial_budget(cfg, fb, oracle.summands(), space);
            let outcome = collect_and_solve(
                &g,
                inst.generator,
                target,
                r,
                inst.cofactor,
                fb,
                seed,
                budget,
                |ops, ctr, point| decompose_binary(inst, fb, &oracle, ops, ctr, point),
            );
            let uses_table = matches!(oracle, Oracle::Mitm { .. });
            let mut v = assemble_variant(
                &name,
                oracle.name(),
                oracle.summands(),
                &fb.cost,
                fb.points.len() as u64,
                fb.abscissae as u64,
                fb.columns as u64,
                fb.dimension,
                uses_table.then_some(&table),
                outcome,
                r,
                space,
                &calib,
                rho_s,
                out.floor_s,
            );
            v.verified &= v.recovered == Some(d);
            out.variants.push(v);
        }
    }
    // Mean over the signed reference only.
    let signed: Vec<f64> = out
        .rho
        .iter()
        .filter(|x| x.automorphisms > 1)
        .map(|x| x.s)
        .collect();
    out.rho_s_mean = signed.iter().sum::<f64>() / signed.len().max(1) as f64;
    Some(out)
}

// ── Ladders and fits ───────────────────────────────────────────────

/// A prime instance for `bits`: the roster's curve when it has one,
/// otherwise a generated prime-order curve.
pub fn prime_instance_for(bits: u32, seed: u64) -> PrimeInstance {
    roster_prime_instance(bits).unwrap_or_else(|| find_prime_order_curve(bits, seed))
}

pub fn run_prime_ladder(cfg: &BoundaryConfig, mut progress: impl FnMut(&str)) -> Vec<RegimeInstance> {
    let mut out = Vec::new();
    for &bits in &cfg.prime_bits {
        progress(&format!("prime: {bits}-bit instance"));
        let inst = prime_instance_for(bits, cfg.seed);
        let res = run_prime_instance(&inst, cfg);
        progress(&format!(
            "prime: {} r=2^{:.1} rho S={:.2} best IC S={:.1}",
            inst.name,
            res.log2_r,
            res.rho_s_mean,
            res.variants.iter().map(|v| v.s).fold(f64::INFINITY, f64::min)
        ));
        out.push(res);
    }
    out
}

pub fn run_char2_ladder(cfg: &BoundaryConfig, mut progress: impl FnMut(&str)) -> Vec<RegimeInstance> {
    let mut out = Vec::new();
    for &n in &cfg.char2_degrees {
        progress(&format!("char2: searching a curve over GF(2^{n})"));
        let Some(inst) = random_binary_instance(n, cfg.seed, 8) else {
            progress(&format!("char2: no curve found at n = {n}"));
            continue;
        };
        let Some(res) = run_char2_instance(&inst, cfg) else {
            progress(&format!("char2: no admissible summand count at n = {n}"));
            continue;
        };
        progress(&format!(
            "char2: {} r=2^{:.1} rho S={:.2} best IC S={:.1}",
            inst.name,
            res.log2_r,
            res.rho_s_mean,
            res.variants.iter().map(|v| v.s).fold(f64::INFINITY, f64::min)
        ));
        out.push(res);
    }
    out
}

/// `K_0` or `K_1` over `F_{2^n}`, whichever has the smaller cofactor
/// (ties to `K_0`); `None` when neither has a usable subgroup.
pub fn koblitz_instance_best(n: u32) -> Option<BinaryInstance> {
    let k0 = koblitz_instance(0, n);
    let k1 = koblitz_instance(1, n);
    match (k0, k1) {
        (Some(a), Some(b)) => Some(if b.cofactor < a.cofactor { b } else { a }),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    }
}

pub fn run_koblitz_ladder(cfg: &BoundaryConfig, mut progress: impl FnMut(&str)) -> Vec<RegimeInstance> {
    let mut out = Vec::new();
    for &n in &cfg.koblitz_degrees {
        let inst = match koblitz_instance_best(n) {
            Some(i) => i,
            None => {
                progress(&format!("koblitz: no usable K_a at n = {n}"));
                continue;
            }
        };
        progress(&format!("koblitz: {} r=2^{:.1}", inst.name, (inst.r as f64).log2()));
        match run_koblitz_instance(&inst, cfg) {
            Some(res) => {
                progress(&format!(
                    "koblitz: {} rho(2n) S={:.3} best IC S={:.1}",
                    inst.name,
                    res.rho_s_mean,
                    res.variants.iter().map(|v| v.s).fold(f64::INFINITY, f64::min)
                ));
                out.push(res);
            }
            None => progress(&format!("koblitz: no factor base at n = {n}")),
        }
    }
    out
}

/// Fit `log₂(gae) = α·log₂(r) + c` per regime, variant and phase over
/// every instance carrying that variant, using the mean over repeats.
pub fn fit_exponents(instances: &[RegimeInstance]) -> Vec<ExponentFit> {
    // (log2 r, log2 #E, log2 value)
    let mut cells: HashMap<(String, String, String), Vec<(f64, f64, f64)>> = HashMap::new();
    for inst in instances {
        let log2_e = inst.group_order.max(inst.r as f64).log2();
        let mut by_variant: HashMap<String, Vec<&VariantResult>> = HashMap::new();
        for v in &inst.variants {
            by_variant.entry(v.name.clone()).or_default().push(v);
        }
        for (name, runs) in by_variant {
            let mean = |f: &dyn Fn(&VariantResult) -> f64| {
                runs.iter().map(|v| f(v)).sum::<f64>() / runs.len() as f64
            };
            let phases: Vec<(&str, f64)> = vec![
                ("total", mean(&|v| v.total_gae)),
                ("factor_base", mean(&|v| v.factor_base.gae)),
                ("relations", mean(&|v| v.relations.gae)),
                ("linear_algebra", mean(&|v| v.linear_algebra.gae)),
                ("trials", mean(&|v| v.trials as f64)),
            ];
            for (phase, value) in phases {
                if value > 0.0 {
                    cells
                        .entry((inst.regime.clone(), name.clone(), phase.to_string()))
                        .or_default()
                        .push((inst.log2_r, log2_e, value.log2()));
                }
            }
        }
        let rho_mean = inst.rho_s_mean * (inst.r as f64).sqrt();
        if rho_mean > 0.0 {
            cells
                .entry((inst.regime.clone(), "rho_reference".into(), "total".into()))
                .or_default()
                .push((inst.log2_r, log2_e, rho_mean.log2()));
        }
    }
    let mut fits: Vec<ExponentFit> = cells
        .into_iter()
        .filter_map(|((regime, variant, phase), pts)| {
            if pts.len() < 3 {
                return None;
            }
            let xs: Vec<f64> = pts.iter().map(|p| p.0).collect();
            let es: Vec<f64> = pts.iter().map(|p| p.1).collect();
            let ys: Vec<f64> = pts.iter().map(|p| p.2).collect();
            let (slope, _, r2) = linear_fit(&xs, &ys)?;
            let (slope_e, _, r2_e) = linear_fit(&es, &ys).unwrap_or((f64::NAN, 0.0, f64::NAN));
            Some(ExponentFit {
                regime,
                variant,
                phase,
                alpha: slope,
                r_squared: r2,
                alpha_group_order: slope_e,
                r_squared_group_order: r2_e,
                points: pts.len(),
            })
        })
        .collect();
    fits.sort_by(|a, b| {
        (&a.regime, &a.variant, &a.phase).cmp(&(&b.regime, &b.variant, &b.phase))
    });
    fits
}

/// Run every regime in the configuration.
pub fn run_all(cfg: &BoundaryConfig, mut progress: impl FnMut(&str)) -> BoundaryLedger {
    let mut instances = Vec::new();
    instances.extend(run_prime_ladder(cfg, &mut progress));
    instances.extend(run_char2_ladder(cfg, &mut progress));
    instances.extend(run_koblitz_ladder(cfg, &mut progress));
    let fits = fit_exponents(&instances);
    BoundaryLedger {
        unit: "S = group-addition equivalents / sqrt(r); native counts exact, conversions measured per instance".into(),
        floor: "generic: sqrt(pi/(2A)) in S, A = usable automorphisms; relations: trials >= (K+1)/min(1, C(F+m-1,m)/#E)".into(),
        reference: "Pollard rho on the same instance, same process, exact operation count, verified [d]G = Q".into(),
        instances,
        fits,
    }
}

// ── Rendering ──────────────────────────────────────────────────────

fn fmt_s(v: f64) -> String {
    if !v.is_finite() {
        "—".into()
    } else if v >= 1000.0 {
        format!("{:.0}", v)
    } else if v >= 10.0 {
        format!("{:.1}", v)
    } else {
        format!("{:.2}", v)
    }
}

/// One table, one unit: every variant on every instance as a row, with
/// the reference and both ratios.
pub fn format_markdown(ledger: &BoundaryLedger) -> String {
    let mut out = String::new();
    out.push_str("| regime | instance | log₂ r | variant | m | \\|F\\| | K | trials | yield/ceiling | S | S rho | rho walk | vs rho | vs floor | FB | rel | LA | ok |\n");
    out.push_str("|:--|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|\n");
    for inst in &ledger.instances {
        let name = inst.curve["name"].as_str().unwrap_or("?").to_string();
        let mut seen: Vec<String> = Vec::new();
        for v in &inst.variants {
            if seen.contains(&v.name) {
                continue;
            }
            seen.push(v.name.clone());
            let runs: Vec<&VariantResult> = inst.variants.iter().filter(|w| w.name == v.name).collect();
            let k = runs.len() as f64;
            let mean = |f: &dyn Fn(&VariantResult) -> f64| runs.iter().map(|w| f(w)).sum::<f64>() / k;
            let ok = runs.iter().all(|w| w.verified);
            let reference: Vec<&RhoResult> = inst
                .rho
                .iter()
                .filter(|x| x.automorphisms == inst.automorphisms_generic || inst.automorphisms_generic <= 2)
                .collect();
            let rho_walk = reference.iter().map(|x| x.s_walk).sum::<f64>() / reference.len().max(1) as f64;
            out.push_str(&format!(
                "| {} | {} | {:.1} | {} | {} | {} | {} | {:.0} | {:.2} | {} | {} | {} | {}× | {}× | {} | {} | {} | {} |\n",
                inst.regime,
                name,
                inst.log2_r,
                v.name,
                v.summands,
                v.signed_points,
                v.columns,
                mean(&|w| w.trials as f64),
                mean(&|w| w.yield_over_ceiling),
                fmt_s(mean(&|w| w.s)),
                fmt_s(inst.rho_s_mean),
                fmt_s(rho_walk),
                fmt_s(mean(&|w| w.s) / inst.rho_s_mean),
                fmt_s(mean(&|w| w.s) / inst.floor_s),
                fmt_s(mean(&|w| w.factor_base.gae) / (inst.r as f64).sqrt()),
                fmt_s(mean(&|w| w.relations.gae) / (inst.r as f64).sqrt()),
                fmt_s(mean(&|w| w.linear_algebra.gae) / (inst.r as f64).sqrt()),
                if ok { "✓" } else { "✗" },
            ));
        }
    }
    out.push_str("\n| regime | variant | phase | α (ops ∝ r^α) | R² | α vs #E | points |\n|:--|:--|:--|--:|--:|--:|--:|\n");
    for f in &ledger.fits {
        out.push_str(&format!(
            "| {} | {} | {} | {:.3} | {:.3} | {:.3} | {} |\n",
            f.regime, f.variant, f.phase, f.alpha, f.r_squared, f.alpha_group_order, f.points
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floors_are_the_textbook_numbers() {
        assert!((generic_floor_s(1.0) - 1.2533).abs() < 1e-3);
        assert!((generic_floor_s(2.0) - 0.8862).abs() < 1e-3);
        let r = 1u64 << 40;
        assert!((generic_floor_ops(r as f64, 1.0) / 2f64.powi(20) - 1.2533).abs() < 1e-3);
        // Counting ceiling: 8 signed points, m = 2, 36 multisets over 100.
        assert!((decomposition_probability_ceiling(8, 2, 100.0) - 0.36).abs() < 1e-9);
        assert!((trials_floor(3, 8, 2, 100.0) - 4.0 / 0.36).abs() < 1e-9);
        assert_eq!(decomposition_probability_ceiling(1000, 3, 10.0), 1.0);
    }

    #[test]
    fn incremental_gauss_pins_the_last_unknown() {
        // x + y = 3, x − y = 1 over Z/101 ⇒ x = 2, y = 1.
        let mut g = IncrementalGauss::new(2, 101);
        assert_eq!(g.add_row(vec![1, 1], 3), RowStatus::Independent);
        assert_eq!(g.pinned(1), None);
        assert_eq!(g.add_row(vec![1, 100], 1), RowStatus::Independent);
        assert_eq!(g.pinned(1), Some(1));
        assert_eq!(g.pinned(0), Some(2));
        assert_eq!(g.add_row(vec![2, 2], 6), RowStatus::Dependent);
        assert_eq!(g.add_row(vec![2, 2], 7), RowStatus::Inconsistent);
        assert!(g.row_ops > 0);
    }

    #[test]
    fn prime_arithmetic_is_a_group_of_the_right_order() {
        let inst = roster_prime_instance(12).unwrap();
        let curve = &inst.curve;
        let g = inst.generator_point();
        assert!(curve.is_on_curve(g));
        let mut ops = GroupOps::default();
        assert!(curve.mul(&mut ops, g, inst.r).infinity);
        assert_eq!(ops.scalar_mults, 1);
        assert!(ops.doubles > 0 && ops.adds > 0);
        // Naive count agrees with the roster order (prime order, h = 1).
        assert_eq!(curve.point_count(), inst.group_order);
        // Semaev S₃ roots are the abscissae of P ± Q.
        let p = curve.mul(&mut ops, g, 17);
        let q = curve.mul(&mut ops, g, 91);
        let sum = curve.add(&mut ops, p, q);
        let diff = curve.add(&mut ops, p, curve.neg(q));
        let (qa, qb, qc) = semaev_s3_quadratic(curve, p.x, q.x);
        for x3 in [sum.x, diff.x] {
            let v = addmod(
                addmod(mulmod(qa, mulmod(x3, x3, curve.p), curve.p), mulmod(qb, x3, curve.p), curve.p),
                qc,
                curve.p,
            );
            assert_eq!(v, 0, "S₃ must vanish on x(P ± Q)");
        }
    }

    #[test]
    fn counted_rho_recovers_and_verifies_on_a_prime_curve() {
        let inst = roster_prime_instance(14).unwrap();
        let g = inst.generator_point();
        let mut ops = GroupOps::default();
        let d = 4321u64 % inst.r;
        let q = inst.curve.mul(&mut ops, g, d);
        let res = rho_reference(&inst.curve, g, q, inst.r, 5, 1 << 24);
        assert_eq!(res.recovered, Some(d));
        assert!(res.verified);
        assert!(res.steps > 0 && res.group_ops.adds >= res.steps);
        assert!(res.steps_over_expected < 8.0, "{res:?}");
    }

    #[test]
    fn artin_schreier_solver_and_point_count_agree_with_the_general_code() {
        use crate::cryptanalysis::koblitz_index_calculus::points_with_x;
        let n = 11u32;
        let irr = find_irreducible_sparse(n).unwrap();
        let gf = Gf2::new(&irr);
        let ash = ArtinSchreier::new(&gf);
        // Every trace-zero c has a solution and every trace-one c none.
        let mut solved = 0;
        for c in 0..(1u64 << n) {
            match ash.solve(c) {
                Some(u) => {
                    assert_eq!(gf.sqr(u) ^ u, c);
                    assert_eq!(ash.trace(c), 0);
                    solved += 1;
                }
                None => assert_eq!(ash.trace(c), 1),
            }
        }
        assert_eq!(solved, 1 << (n - 1));
        // The batched count agrees with a brute-force one.
        let (a, b) = (1u64, 0x3A5u64 & gf.mask);
        let curve = BinaryCurve {
            m: n,
            irreducible: irr.clone(),
            a: gf.to_element(a),
            b: gf.to_element(b),
            generator: BinaryPoint::Infinity,
            order: BigUint::from(1u32),
            cofactor: BigUint::from(1u32),
        };
        let brute = 1 + (0..(1u64 << n))
            .map(|x| points_with_x(&curve, &gf.to_element(x)).len() as u64)
            .sum::<u64>();
        assert_eq!(binary_point_count(&gf, &ash, a, b), brute);
    }

    #[test]
    fn random_binary_instance_has_a_subgroup_of_the_claimed_order() {
        let inst = random_binary_instance(13, 3, 8).expect("a curve at n = 13");
        let g = BinaryGroup(&inst.fast);
        let mut ops = GroupOps::default();
        assert!(inst.fast.is_on_curve(inst.generator));
        assert!(g.mul(&mut ops, inst.generator, inst.r).infinity);
        assert!(!g.mul(&mut ops, inst.generator, 1).infinity);
        assert!(is_prime_u64(inst.r));
        assert_eq!(inst.r * inst.cofactor, inst.group_order);
        assert!(inst.cofactor <= 8);
    }

    #[test]
    fn prime_pipeline_recovers_the_planted_logarithm_with_every_oracle() {
        let cfg = BoundaryConfig {
            repeats: 1,
            ..BoundaryConfig::quick()
        };
        let inst = roster_prime_instance(12).unwrap();
        let res = run_prime_instance(&inst, &cfg);
        assert!(res.rho_verified_all, "{:?}", res.rho);
        assert_eq!(res.variants.len(), 4);
        for v in &res.variants {
            assert!(v.verified, "{} did not verify: {v:?}", v.name);
            assert!(v.total_gae > 0.0 && v.s.is_finite());
            assert!(v.relations.get("trials") == v.trials);
            assert!(v.rank >= 1 && v.relations_found >= v.rank);
        }
        let s3 = res.variants.iter().find(|v| v.name == "semaev_s3_roots_m2").unwrap();
        assert!(s3.relations.get("sqrt_solves") > 0);
        let mitm = res.variants.iter().find(|v| v.name == "mitm_m3").unwrap();
        assert!(mitm.factor_base.get("pair_table_entries") > 0);
    }

    #[test]
    fn char2_pipeline_recovers_the_planted_logarithm() {
        let cfg = BoundaryConfig {
            repeats: 1,
            ..BoundaryConfig::quick()
        };
        let inst = random_binary_instance(12, 9, 8).expect("a curve at n = 12");
        let res = run_char2_instance(&inst, &cfg).expect("an admissible summand count");
        assert!(res.rho_verified_all, "{:?}", res.rho);
        assert!(!res.variants.is_empty());
        for v in &res.variants {
            assert!(v.verified, "{} did not verify: {v:?}", v.name);
        }
        if let Some(s4) = res
            .variants
            .iter()
            .find(|v| v.name == "semaev_s4_pairs_and_solve_m3")
        {
            assert!(s4.relations.get("s4_pairs") > 0);
            assert!(res.calibration.ns_per_s4_pair.unwrap() > 0.0);
        }
    }

    #[test]
    #[ignore = "diagnostic: stage timings of the Koblitz pipeline at n = 15"]
    fn koblitz_stage_timing_probe() {
        let t = Instant::now();
        let inst = koblitz_instance_best(15).expect("K_a / 2^15");
        eprintln!("instance {} r={} h={} in {:?}", inst.name, inst.r, inst.cofactor, t.elapsed());
        let t = Instant::now();
        let (frob, desc) = choose_koblitz_base(&inst, 1, 6000).expect("base");
        eprintln!("base {desc}: {} points in {:?}", frob.points.len(), t.elapsed());
        let t = Instant::now();
        let folded = koblitz_factor_base(&inst, &frob, ColumnFold::SignedFrobeniusOrbit, desc.clone()).unwrap();
        eprintln!("folded: {} points, {} columns in {:?}", folded.points.len(), folded.columns, t.elapsed());
        let g = BinaryGroup(&inst.fast);
        let t = Instant::now();
        let table = PairTable::build(&g, &folded);
        eprintln!("table: {} entries in {:?}", table.entries, t.elapsed());
        let mut calib = Calibration::default();
        let t = Instant::now();
        calibrate_binary(&inst, &mut calib);
        eprintln!("calibration in {:?}: {:?}", t.elapsed(), calib);
        let mut ops = GroupOps::default();
        let d = 123 % inst.r;
        let target = g.mul(&mut ops, inst.generator, d);
        let t = Instant::now();
        let plain = rho_reference(&g, inst.generator, target, inst.r, 5, rho_cap(inst.r, 64.0));
        eprintln!("plain rho: {:?} in {:?}", plain.recovered, t.elapsed());
        // Ground truth: how many random subgroup targets are 3-sums of
        // base points, by exhaustive triples, and does the oracle agree?
        let mut rng = StdRng::seed_from_u64(77);
        let (mut truth_hits, mut oracle_hits, mut checked) = (0, 0, 0);
        for _ in 0..40 {
            let q = g.mul(&mut ops, inst.generator, rng.gen_range(1..inst.r));
            let n = folded.points.len();
            let mut truth = false;
            'triples: for i in 0..n {
                for j in i..n {
                    let pij = g.add(&mut ops, folded.points[i], folded.points[j]);
                    for k in j..n {
                        if g.add(&mut ops, pij, folded.points[k]) == q {
                            truth = true;
                            break 'triples;
                        }
                    }
                }
            }
            let mut ctr = OracleCounters::default();
            let found = decompose_mitm(&g, &folded, &table, 3, &mut ops, &mut ctr, q);
            checked += 1;
            truth_hits += truth as u32;
            oracle_hits += found.is_some() as u32;
            if truth != found.is_some() {
                eprintln!("DISAGREE on target {q:?}: truth {truth} oracle {found:?}");
            }
        }
        eprintln!("ground truth: {truth_hits}/{checked} targets decompose; oracle found {oracle_hits}");
        let t = Instant::now();
        let outcome = collect_and_solve(&g, inst.generator, target, inst.r, inst.cofactor, &folded, 5, 1_000_000, |ops, ctr, point| {
            decompose_binary(&inst, &folded, &Oracle::Mitm { table: &table, m: 3 }, ops, ctr, point)
        });
        eprintln!(
            "mitm m3 folded: recovered {:?} verified {} trials {} relations {} independent {} dependent {} inconsistent {} in {:?}",
            outcome.recovered, outcome.verified, outcome.trials, outcome.relations_found, outcome.independent, outcome.dependent,
            outcome.linear_algebra.get("inconsistent"), t.elapsed()
        );
        let kc = inst.koblitz.as_ref().unwrap();
        let t = Instant::now();
        let report = koblitz_signed_frobenius_rho_reference(kc, &inst.fast.lower(target), &KoblitzSignedRhoOptions { seed: 5, ..Default::default() }, &mut |_| {});
        eprintln!("signed rho: {:?} verified {} iterations {} in {:?}", report.recovered_log, report.verified, report.iterations, t.elapsed());
    }

    #[test]
    fn koblitz_pipeline_recovers_the_planted_logarithm_and_folds_columns() {
        let cfg = BoundaryConfig {
            repeats: 1,
            ..BoundaryConfig::quick()
        };
        let inst = koblitz_instance_best(15).expect("K_a / 2^15");
        let res = run_koblitz_instance(&inst, &cfg).expect("a factor base");
        assert!(res.rho_verified_all, "{:?}", res.rho);
        let folded = res
            .variants
            .iter()
            .find(|v| v.name.starts_with("mitm_m") && v.name.ends_with("signed_orbit_columns"))
            .unwrap();
        let control = res
            .variants
            .iter()
            .find(|v| v.name.starts_with("mitm_m") && v.name.ends_with("abscissa_columns_control"))
            .unwrap();
        assert!(folded.verified && control.verified, "{res:?}");
        assert!(folded.columns < control.columns, "the fold must cut the columns");
        assert_eq!(folded.signed_points, control.signed_points);
        assert!(res.automorphisms_generic == 30);
        assert!(res.floor_s < generic_floor_s(2.0));
    }

    #[test]
    fn exponent_fit_reads_a_planted_slope() {
        let mk = |log2_r: f64, gae: f64| RegimeInstance {
            regime: "t".into(),
            curve: serde_json::json!({"name": "x"}),
            r: 2f64.powf(log2_r) as u64,
            log2_r,
            group_order: 0.0,
            cofactor: 1,
            automorphisms_generic: 1,
            floor_s: 1.0,
            floor_ops: 1.0,
            calibration: Calibration::default(),
            rho: Vec::new(),
            rho_s_mean: 0.0,
            rho_verified_all: true,
            variants: vec![VariantResult {
                name: "v".into(),
                oracle: "o".into(),
                summands: 2,
                factor_base: PhaseCost::default(),
                signed_points: 1,
                abscissae: 1,
                columns: 1,
                dimension: None,
                relations: PhaseCost::default(),
                linear_algebra: PhaseCost::default(),
                verify: PhaseCost::default(),
                trials: 1,
                relations_found: 1,
                independent: 1,
                dependent: 0,
                rank: 1,
                trials_floor: 1.0,
                trials_over_floor: 1.0,
                decomposition_probability_ceiling: 1.0,
                yield_over_ceiling: 1.0,
                total_gae: gae,
                total_wall_ns: 0,
                s: 1.0,
                s_wall: 1.0,
                ratio_to_rho: 1.0,
                ratio_to_floor: 1.0,
                recovered: Some(1),
                verified: true,
            }],
            seeds: Vec::new(),
            targets: Vec::new(),
        };
        let insts: Vec<RegimeInstance> = [10.0, 14.0, 18.0, 22.0]
            .iter()
            .map(|&b| mk(b, 3.0 * 2f64.powf(0.75 * b)))
            .collect();
        let fits = fit_exponents(&insts);
        let total = fits.iter().find(|f| f.variant == "v" && f.phase == "total").unwrap();
        assert!((total.alpha - 0.75).abs() < 1e-9, "{total:?}");
        assert!(total.r_squared > 0.999);
    }
}
