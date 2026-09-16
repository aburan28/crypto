//! # Index calculus vs Pollard rho on the same hyperelliptic Jacobian.
//!
//! [`crate::cryptanalysis::hyperelliptic_index_calculus`] implements the
//! Adleman–DeMarrais–Huang / Gaudry attack but claims nothing about its
//! cost.  This module supplies the missing half: the **reference**
//! `AGENTS.md` requires — Pollard rho run on the *same* group, the same
//! `D₁`, `D₂` and the same prime order `N`, with the same operation
//! accounting — and a head-to-head table in one unit.
//!
//! ## Unit
//!
//! ```text
//!   S = total group operations / sqrt(N)
//! ```
//!
//! which holds rho flat at every size: a collision arrives after about
//! `sqrt(π N / 2) ≈ 1.2533·sqrt(N)` distinct points, and Floyd spends
//! three group operations per iteration, so this reference sits near
//! `S ≈ 3` (a distinguished-point rho would sit near `1.25`).  Either
//! way it is a constant, which is what makes "is index calculus better
//! than rho here" one column to read.
//!
//! Every cost each side incurs is inside `S`:
//!
//! - **rho** — the `r` precomputed walk steps `R_j = a_j D₁ + b_j D₂`
//!   plus every tortoise/hare step.
//! - **index calculus** — the relation-stage group operations, *and* the
//!   smoothness oracle's field work, *and* the linear algebra, the last
//!   two converted into group-operation equivalents by a **measured**
//!   factor
//!   ([`calibrate_modmuls_per_group_op`]), not an assumed one.  Leaving
//!   the solve out is exactly the relabelling failure mode: the work
//!   does not disappear, it moves somewhere the headline does not look.
//!
//! ## Boundaries
//!
//! Two, both derived before measuring:
//!
//! - **Floor (index calculus).**  `m + 1` unknowns need `m + 1`
//!   independent relations.  Under heuristic **H1** — a reduced divisor
//!   drawn uniformly from `Jac(C)(F_p)` splits into degree-1 places with
//!   probability `≈ 1/g!` — the expected number of trials is at least
//!   `(m+1)·g!`, and no walk can produce a trial for less than one group
//!   operation.  So
//!
//!   ```text
//!     T = (m + 1)·g!                        (trials)
//!     floor_ops = T  +  (T·g + m²)/c
//!   ```
//!
//!   one group operation per trial, plus the field work no
//!   implementation can skip: the smoothness oracle must read each
//!   candidate's `u`, whose `g` coefficients cost `g` multiplications,
//!   and any dense solve must read its own `m × m` matrix once.  Those
//!   `T·g + m²` mul-mods are divided by the measured
//!   mul-mods-per-group-op factor `c` so every term is in the same unit
//!   as the measurement.  It moves only with `m` and `g`, so it cannot
//!   be tuned away; it is a floor on *this* algorithm, not on the DLP.
//! - **Reference.**  Pollard rho, measured, on the same instances.
//!
//! The ratio `S_ic / S_rho` and the ratio `S_ic` to the floor are the
//! progress measures.  This module *measures*; the numbers it produces
//! at toy `p` say nothing about cryptographic sizes, where the relation
//! stage's `O(p)` root finding alone would dominate.

use std::time::Instant;

use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::rngs::StdRng;
use rand::{RngCore, SeedableRng};

use crate::prime_hyperelliptic::{HyperellipticCurveP, MumfordDivisorP};
use crate::utils::mod_inverse;

use super::hyperelliptic_index_calculus::{
    build_factor_base, hec_index_calculus_dlp, HecIndexCalculusParams,
};

// ── Pollard rho on the Jacobian ────────────────────────────────────────

/// Number of branches in the additive walk.  Teske's analysis puts the
/// cycle length of an `r`-adding walk within a few percent of the
/// random-map ideal once `r ≥ 16`, so a smaller `r` would flatter the
/// index calculus by handicapping its reference.
const RHO_BRANCHES: usize = 16;

/// Bit length of the rho branch coefficients `a_j, b_j`.
///
/// The index-calculus walk builds its steps from short coefficients
/// because precomputation is its largest term once a trial costs one
/// operation.  The same trick is available to rho and costs it nothing,
/// so the reference gets it too: optimising one side's precomputation
/// and not the other's would be a rigged comparison, and the rigging
/// would favour the side under study.
const RHO_STEP_BITS: u32 = 8;

/// Outcome of one rho run.
#[derive(Clone, Debug)]
pub struct RhoResult {
    pub k: Option<BigUint>,
    /// Jacobian group operations: walk steps plus the precomputed
    /// branch divisors.
    pub group_ops: usize,
    /// Of those, the branch precomputation — a fixed cost that does not
    /// grow with `sqrt(N)` and so dominates `S` at toy sizes.  Split
    /// out because a reader comparing against the `sqrt(π/2)` ideal
    /// needs the walk alone, and hiding it inside one number would make
    /// rho look worse than it is.
    pub precompute_ops: usize,
    /// Of those, the tortoise/hare steps.
    pub walk_ops: usize,
    pub restarts: usize,
    pub wall_ms: f64,
}

/// Branch selector: a cheap hash of the Mumford representation.
///
/// Any function of the *class* works, as long as equal classes select
/// equal branches — that is what makes the walk deterministic and the
/// collision meaningful.
fn branch_of(d: &MumfordDivisorP) -> usize {
    let mut acc: u64 = 0;
    for poly in [&d.u, &d.v] {
        for limb in poly.coeffs.iter().flat_map(|c| c.to_u64_digits()) {
            acc = acc.wrapping_mul(0x9e37_79b9_7f4a_7c15).wrapping_add(limb);
        }
        acc = acc.wrapping_mul(0x9e37_79b9_7f4a_7c15).wrapping_add(1);
    }
    (acc >> 32) as usize % RHO_BRANCHES
}

fn rand_below(rng: &mut StdRng, n: &BigUint) -> BigUint {
    let bytes = (n.bits() as usize).div_ceil(8) + 8;
    let mut buf = vec![0u8; bytes];
    rng.fill_bytes(&mut buf);
    BigUint::from_bytes_be(&buf) % n
}

/// **Solve `D₂ = k·D₁`** by Pollard rho with an `r`-adding walk and
/// Floyd cycle finding, over the order-`N` subgroup generated by `D₁`.
///
/// `n` must be prime (the final division inverts mod `n`).  A collision
/// with `b₁ ≡ b₂` carries no information, so the walk restarts from a
/// fresh seed; `max_steps` bounds the total.  The returned `k` is
/// verified before it is returned.
pub fn pollard_rho_jacobian(
    curve: &HyperellipticCurveP,
    d1: &MumfordDivisorP,
    d2: &MumfordDivisorP,
    n: &BigUint,
    seed: u64,
    max_steps: usize,
) -> RhoResult {
    let started = Instant::now();
    let mut group_ops = 0usize;
    let mut precompute_ops = 0usize;
    let mut restarts = 0usize;
    let mut steps_left = max_steps;

    for attempt in 0..64u64 {
        let mut rng = StdRng::seed_from_u64(seed ^ (attempt << 32));

        // Precomputed branches R_j = a_j·D₁ + b_j·D₂.  Charged.
        let mut branch: Vec<(BigUint, BigUint, MumfordDivisorP)> = Vec::with_capacity(RHO_BRANCHES);
        let step_bound = BigUint::one() << RHO_STEP_BITS;
        for _ in 0..RHO_BRANCHES {
            let a = rand_below(&mut rng, &step_bound) + BigUint::one();
            let b = rand_below(&mut rng, &step_bound) + BigUint::one();
            let r = d1
                .scalar_mul(&a, curve)
                .add(&d2.scalar_mul(&b, curve), curve);
            group_ops += 2 * RHO_STEP_BITS as usize;
            precompute_ops += 2 * RHO_STEP_BITS as usize;
            branch.push((a, b, r));
        }

        let step = |a: &BigUint, b: &BigUint, x: &MumfordDivisorP| {
            let (aj, bj, rj) = &branch[branch_of(x)];
            ((a + aj) % n, (b + bj) % n, x.add(rj, curve))
        };

        // Start away from the identity so the first branch is not
        // degenerate.
        let a0 = rand_below(&mut rng, n);
        let b0 = rand_below(&mut rng, n);
        let x0 = d1
            .scalar_mul(&a0, curve)
            .add(&d2.scalar_mul(&b0, curve), curve);
        group_ops += 2 * (n.bits() as usize + 1);
        precompute_ops += 2 * (n.bits() as usize + 1);

        let (mut at, mut bt, mut xt) = (a0.clone(), b0.clone(), x0.clone());
        let (mut ah, mut bh, mut xh) = (a0, b0, x0);

        while steps_left > 0 {
            let (na, nb, nx) = step(&at, &bt, &xt);
            at = na;
            bt = nb;
            xt = nx;
            group_ops += 1;

            let (na, nb, nx) = step(&ah, &bh, &xh);
            let (na, nb, nx) = step(&na, &nb, &nx);
            ah = na;
            bh = nb;
            xh = nx;
            group_ops += 2;
            steps_left = steps_left.saturating_sub(3);

            if xt == xh {
                // a_t + b_t k ≡ a_h + b_h k  ⟹  k ≡ (a_t − a_h)/(b_h − b_t).
                let db = (&bh + n - &bt) % n;
                if db.is_zero() {
                    restarts += 1;
                    break; // useless collision; new walk
                }
                let da = (&at + n - &ah) % n;
                let k = match mod_inverse(&db, n) {
                    Some(inv) => (da * inv) % n,
                    None => {
                        restarts += 1;
                        break;
                    }
                };
                if &d1.scalar_mul(&k, curve) == d2 {
                    return RhoResult {
                        k: Some(k),
                        group_ops,
                        precompute_ops,
                        walk_ops: group_ops - precompute_ops,
                        restarts,
                        wall_ms: started.elapsed().as_secs_f64() * 1e3,
                    };
                }
                restarts += 1;
                break;
            }
        }
        if steps_left == 0 {
            break;
        }
    }

    RhoResult {
        k: None,
        group_ops,
        precompute_ops,
        walk_ops: group_ops - precompute_ops,
        restarts,
        wall_ms: started.elapsed().as_secs_f64() * 1e3,
    }
}

// ── Unit conversion ────────────────────────────────────────────────────

/// Measure how many `BigUint` mul-mod-`n` operations fit in the time of
/// one Jacobian group operation on this curve.
///
/// The linear-algebra stage is counted in mul-mods and the walks in
/// group operations; without this factor the two cannot be added, and
/// adding them anyway is how a cost gets hidden rather than removed.
/// Reported alongside every table so a reader can redo the arithmetic.
pub fn calibrate_modmuls_per_group_op(
    curve: &HyperellipticCurveP,
    d: &MumfordDivisorP,
    n: &BigUint,
    rounds: usize,
) -> f64 {
    let mut acc = d.clone();
    let start = Instant::now();
    for _ in 0..rounds {
        acc = acc.add(d, curve);
    }
    let per_group_op = start.elapsed().as_secs_f64() / rounds as f64;
    std::hint::black_box(&acc);

    let a: BigUint = (n - BigUint::one()) >> 1u32;
    let mut x = a.clone();
    let start = Instant::now();
    for _ in 0..rounds {
        x = (&x * &a) % n;
    }
    let per_modmul = start.elapsed().as_secs_f64() / rounds as f64;
    std::hint::black_box(&x);

    if per_modmul <= 0.0 {
        return f64::NAN;
    }
    per_group_op / per_modmul
}

// ── Head-to-head ───────────────────────────────────────────────────────

/// One row of the comparison table: both algorithms, one instance, one
/// unit.
#[derive(Clone, Debug)]
pub struct HeadToHeadRow {
    pub p: u64,
    pub genus: u32,
    pub jac_order: BigUint,
    /// Prime order of `D₁` — the group both algorithms actually work in.
    pub n: BigUint,
    pub factor_base_size: usize,

    /// Means over `HeadToHeadTrials::ic` runs.
    pub ic_relation_ops: f64,
    /// Of those, walk-step precomputation and restarts.
    pub ic_precompute_ops: f64,
    /// Group operations per trial, precomputation excluded.
    pub ic_ops_per_trial: f64,
    /// Linear algebra, in mul-mods.
    pub ic_la_modmuls: f64,
    /// Smoothness oracle (root finding + decomposition), in mul-mods.
    pub ic_oracle_modmuls: f64,
    /// The same, converted to group-operation equivalents.
    pub ic_oracle_group_equiv: f64,
    /// The same, converted to group-operation equivalents.
    pub ic_la_group_equiv: f64,
    pub ic_total_group_ops: f64,
    pub ic_s: f64,
    pub ic_wall_ms: f64,
    pub ic_correct: bool,
    pub ic_smoothness_rate: f64,

    /// Means over `HeadToHeadTrials::rho` runs.
    pub rho_group_ops: f64,
    pub rho_precompute_ops: f64,
    pub rho_s: f64,
    /// `S` counting only the walk — this is the column that should sit
    /// at the `sqrt(π/2) ≈ 1.25` ideal.
    pub rho_walk_s: f64,
    pub rho_wall_ms: f64,
    pub rho_correct: bool,

    /// `(m+1)·g! + m²/c`, derived — see the module docs.
    pub ic_floor_ops: f64,
    pub ic_floor_s: f64,
    pub modmuls_per_group_op: f64,
}

impl HeadToHeadRow {
    /// `S_ic / S_rho`: below 1 means index calculus won this instance.
    pub fn ratio_to_reference(&self) -> f64 {
        self.ic_s / self.rho_s
    }

    /// `S_ic` against the derived floor.  A ratio near 1 means the
    /// implementation is close to what the method can do; it says
    /// nothing about whether the method is any good.
    pub fn ratio_to_floor(&self) -> f64 {
        self.ic_s / self.ic_floor_s
    }
}

fn factorial(g: u32) -> f64 {
    (1..=g as u64).map(|v| v as f64).product::<f64>().max(1.0)
}

fn sqrt_big(n: &BigUint) -> f64 {
    // Toy sizes; f64 is exact enough for a reporting denominator and
    // the alternative (integer sqrt) hides nothing here.
    let as_f = n
        .to_u64_digits()
        .iter()
        .rev()
        .fold(0f64, |acc, &d| acc * 2f64.powi(64) + d as f64);
    as_f.sqrt()
}

/// How many independent runs of each side a row averages over.
///
/// Both algorithms are randomised and a single run of rho has a very
/// wide distribution — its step count is a Rayleigh-like variable, not
/// a constant — so a one-seed row would report sampling noise as a
/// difference between methods.  These are the fixed sample counts of
/// this comparison; changing them changes the protocol.
#[derive(Clone, Debug)]
pub struct HeadToHeadTrials {
    pub ic: usize,
    pub rho: usize,
}

impl Default for HeadToHeadTrials {
    fn default() -> Self {
        // Rho's step count is a wide random variable and the measured
        // gap is now a factor of ~2, not ~10, so the sample has to be
        // large enough that the gap is not the sampling noise.  25 runs
        // put the standard error of the mean near 20% of one run's
        // spread; 9 did not.
        Self { ic: 5, rho: 25 }
    }
}

/// Run both algorithms on one instance and score them in `S`, averaging
/// each side over `trials` independent seeds.
///
/// `d1` must have prime order `n`, and `d2 = k·D₁` for the `k` being
/// recovered.  Every run is verified independently; `*_correct` is true
/// only if **all** runs of that side returned the right `k`, and a row
/// without it is not a result, per `AGENTS.md`.
#[allow(clippy::too_many_arguments)]
pub fn head_to_head(
    curve: &HyperellipticCurveP,
    d1: &MumfordDivisorP,
    d2: &MumfordDivisorP,
    n: &BigUint,
    ic_params: &HecIndexCalculusParams,
    rho_seed: u64,
    rho_max_steps: usize,
    expected_k: &BigUint,
    trials: &HeadToHeadTrials,
) -> HeadToHeadRow {
    let p_u = curve.p.to_u64_digits().first().copied().unwrap_or(0);
    let conv = calibrate_modmuls_per_group_op(curve, d1, n, 2_000);

    let ic_runs = trials.ic.max(1);
    let mut ic_relation_ops = 0f64;
    let mut ic_la_modmuls = 0f64;
    let mut ic_oracle_modmuls = 0f64;
    let mut ic_precompute_ops = 0f64;
    let mut ic_ops_per_trial = 0f64;
    let mut ic_wall_ms = 0f64;
    let mut ic_smooth = 0f64;
    let mut ic_fb = 0usize;
    let mut ic_correct = true;
    for t in 0..ic_runs {
        let mut params = ic_params.clone();
        params.seed = ic_params.seed.wrapping_add(t as u64);
        let ic_start = Instant::now();
        let (ic_k, rep) = hec_index_calculus_dlp(curve, d1, d2, n, &params);
        ic_wall_ms += ic_start.elapsed().as_secs_f64() * 1e3;
        ic_relation_ops += rep.jacobian_ops as f64;
        ic_la_modmuls += rep.solve_row_ops as f64;
        ic_oracle_modmuls += rep.smoothness_field_ops as f64;
        ic_precompute_ops += rep.precompute_ops as f64;
        ic_ops_per_trial += rep.ops_per_trial();
        ic_smooth += rep.smoothness_rate();
        ic_fb = rep.factor_base_size;
        ic_correct &= ic_k.as_ref() == Some(expected_k);
    }
    let f = ic_runs as f64;
    let (ic_precompute_ops, ic_ops_per_trial) = (ic_precompute_ops / f, ic_ops_per_trial / f);
    let (ic_relation_ops, ic_la_modmuls, ic_oracle_modmuls, ic_wall_ms, ic_smooth) = (
        ic_relation_ops / f,
        ic_la_modmuls / f,
        ic_oracle_modmuls / f,
        ic_wall_ms / f,
        ic_smooth / f,
    );

    let rho_runs = trials.rho.max(1);
    let mut rho_ops = 0f64;
    let mut rho_pre = 0f64;
    let mut rho_walk = 0f64;
    let mut rho_wall_ms = 0f64;
    let mut rho_correct = true;
    for t in 0..rho_runs {
        let r = pollard_rho_jacobian(
            curve,
            d1,
            d2,
            n,
            rho_seed.wrapping_add(t as u64),
            rho_max_steps,
        );
        rho_ops += r.group_ops as f64;
        rho_pre += r.precompute_ops as f64;
        rho_walk += r.walk_ops as f64;
        rho_wall_ms += r.wall_ms;
        rho_correct &= r.k.as_ref() == Some(expected_k);
    }
    let g = rho_runs as f64;
    let (rho_ops, rho_pre, rho_walk, rho_wall_ms) =
        (rho_ops / g, rho_pre / g, rho_walk / g, rho_wall_ms / g);

    let root_n = sqrt_big(n);
    let la_group_equiv = ic_la_modmuls / conv.max(f64::MIN_POSITIVE);
    let oracle_group_equiv = ic_oracle_modmuls / conv.max(f64::MIN_POSITIVE);
    let ic_total = ic_relation_ops + la_group_equiv + oracle_group_equiv;

    let m = ic_fb as f64;
    // Both terms in group-operation equivalents: the solve's `m²`
    // matrix read is `m²` mul-mods, which is `m²/conv` group ops.  A
    // floor stated in a different unit from the measurement is not a
    // floor — it was one of these that first read as "below the
    // floor", which is impossible and meant the unit, not the run.
    // The oracle is not optional either: every trial must be tested for
    // smoothness, and the cheapest test that reads `u` at all touches
    // its `g` coefficients.  Charging zero for it would put the floor
    // below what any implementation can reach.
    let trials_floor = (m + 1.0) * factorial(curve.genus);
    let floor_ops =
        trials_floor + (trials_floor * curve.genus as f64 + m * m) / conv.max(f64::MIN_POSITIVE);

    HeadToHeadRow {
        p: p_u,
        genus: curve.genus,
        jac_order: BigUint::zero(), // filled by the caller when known
        n: n.clone(),
        factor_base_size: ic_fb,

        ic_relation_ops,
        ic_precompute_ops,
        ic_ops_per_trial,
        ic_la_modmuls,
        ic_la_group_equiv: la_group_equiv,
        ic_oracle_modmuls,
        ic_oracle_group_equiv: oracle_group_equiv,
        ic_total_group_ops: ic_total,
        ic_s: ic_total / root_n,
        ic_wall_ms,
        ic_correct,
        ic_smoothness_rate: ic_smooth,

        rho_group_ops: rho_ops,
        rho_precompute_ops: rho_pre,
        rho_s: rho_ops / root_n,
        rho_walk_s: rho_walk / root_n,
        rho_wall_ms,
        rho_correct,

        ic_floor_ops: floor_ops,
        ic_floor_s: floor_ops / root_n,
        modmuls_per_group_op: conv,
    }
}

/// Convenience: the factor-base size for `curve`, which is what the
/// floor and the solve both scale in.
pub fn full_factor_base_size(curve: &HyperellipticCurveP) -> usize {
    build_factor_base(curve, usize::MAX).len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::hyperelliptic_index_calculus::{
        prime_order_of, subgroup_generator, LinearAlgebra, RelationSearch, SmoothnessTest,
    };
    use crate::prime_hyperelliptic::{brute_force_jac_order_via_lpoly, FpPoly};

    fn toy_curve(p: u64) -> HyperellipticCurveP {
        let p = BigUint::from(p);
        let f = FpPoly::from_coeffs(
            vec![
                BigUint::from(1u32),
                BigUint::from(1u32),
                BigUint::from(2u32),
                BigUint::from(3u32),
                BigUint::zero(),
                BigUint::from(1u32),
            ],
            p.clone(),
        );
        HyperellipticCurveP::new(p, f, 2)
    }

    fn instance(
        p: u64,
    ) -> (
        HyperellipticCurveP,
        MumfordDivisorP,
        MumfordDivisorP,
        BigUint,
        BigUint,
    ) {
        let curve = toy_curve(p);
        let jac = brute_force_jac_order_via_lpoly(&curve);
        let fb = build_factor_base(&curve, usize::MAX);
        // Largest prime factor of #Jac.
        let mut l = BigUint::one();
        let mut rest = jac.clone();
        let mut d = BigUint::from(2u32);
        while &d * &d <= rest {
            while (&rest % &d).is_zero() {
                rest /= &d;
                if d > l {
                    l = d.clone();
                }
            }
            d += 1u32;
        }
        if rest > l {
            l = rest;
        }
        let d1 = subgroup_generator(&curve, &fb, &jac, &l).expect("generator");
        assert_eq!(prime_order_of(&curve, &d1, &jac).as_ref(), Some(&l));
        let k = &l / BigUint::from(3u32) + BigUint::from(5u32);
        let d2 = d1.scalar_mul(&k, &curve);
        (curve, d1, d2, l, k)
    }

    #[test]
    fn rho_solves_the_same_instance() {
        let (curve, d1, d2, n, k) = instance(41);
        let r = pollard_rho_jacobian(&curve, &d1, &d2, &n, 1, 1_000_000);
        assert_eq!(r.k, Some(k));
        assert!(r.group_ops > 0);
    }

    #[test]
    fn head_to_head_scores_both_sides() {
        let (curve, d1, d2, n, k) = instance(41);
        let params = HecIndexCalculusParams {
            fb_size: usize::MAX,
            extra_relations: 5,
            max_trials: 50_000,
            seed: 20260916,
            search: RelationSearch::Random,
            linear_algebra: LinearAlgebra::Dense,
            smoothness: SmoothnessTest::Scan,
        };
        let row = head_to_head(
            &curve,
            &d1,
            &d2,
            &n,
            &params,
            1,
            1_000_000,
            &k,
            &HeadToHeadTrials { ic: 2, rho: 3 },
        );
        assert!(
            row.ic_correct,
            "index calculus row without a verified answer"
        );
        assert!(row.rho_correct, "rho row without a verified answer");
        assert!(row.ic_s > 0.0 && row.rho_s > 0.0);
        assert!(row.modmuls_per_group_op.is_finite());
        // The floor is a floor: no run may come in under it.
        assert!(
            row.ic_total_group_ops >= row.ic_floor_ops * 0.999,
            "measured {} below the derived floor {}",
            row.ic_total_group_ops,
            row.ic_floor_ops
        );
    }
}
