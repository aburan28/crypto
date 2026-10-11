//! **Residual matching** (lead L3 of `docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md`):
//! relations from collisions between *partial* decompositions instead of
//! from full decompositions of single targets.
//!
//! A state is `(a, b, row)` with residual `L = [a]G + [b]Q − Σ_j P_j`, the
//! row recording which base columns were subtracted with which
//! coefficients.  Two states with the same residual (up to sign) give
//!
//! ```text
//!   (a − a')G + (b − b')Q = Σ_j P_j − Σ_j P'_j ,
//! ```
//!
//! a relation over the base alone, exactly as the large-prime variation
//! matches two relations on their large prime.  A residual that itself
//! lies in `±F` is a complete decomposition and enters as a classic row.
//!
//! The walk is strategy C1 of `residual_walk.rs` on the framework's own
//! groups and bases: an r-adding walk with sixteen tracked jumps, and at
//! every step one base point chosen by the residual's hash is subtracted
//! (one more addition), so equal residuals have equal successors and the
//! walk is collision-preserving.  Every residual is stored (exact table,
//! no distinguished points), so a collision is found the first time it
//! happens; the table is capped and the walk keeps checking after the
//! cap.  A collision restarts the walk from a fresh `[a]G + [b]Q`.
//!
//! What the row costs is the honest point of the lead: a collision on
//! `r` points is a birthday event, so the expected steps per relation are
//! of the order of one rho solve, whatever the base.  The measurement
//! reports steps and additions per relation beside the baseline's.
use std::collections::HashMap;
use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::cryptanalysis::ic_boundary::{
    CountedGroup, FactorBase, GroupOps, OracleCounters, PhaseCost, PipelineOutcome, RelationSolver,
    RowStatus,
};

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}
fn addmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 + b as u128) % m as u128) as u64
}
fn submod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 + m as u128 - b as u128) % m as u128) as u64
}
fn mix(v: u64) -> u64 {
    let mut z = v.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

struct Stored {
    a: u64,
    b: u64,
    row: Vec<u64>,
}

/// Largest number of residuals kept; beyond it the walk still probes.
pub const TABLE_CAP: usize = 1 << 22;

/// Collect relations by residual collisions until the target column is
/// pinned, then verify.  `max_trials` bounds the walk steps.
#[allow(clippy::too_many_arguments)]
pub fn collect_residual_and_solve<G: CountedGroup, L: RelationSolver + ?Sized>(
    g: &G,
    generator: G::Elt,
    target: G::Elt,
    r: u64,
    h: u64,
    fb: &FactorBase<G::Elt>,
    seed: u64,
    max_trials: u64,
    matrix: &mut L,
) -> PipelineOutcome {
    let mut rng = StdRng::seed_from_u64(seed ^ 0x5245_5349_4455_414C);
    let mut rel = PhaseCost::default();
    let mut la = PhaseCost::default();
    let mut ver = PhaseCost::default();
    let ctr = OracleCounters::default();
    let d_col = fb.columns;
    let cols = d_col + 1;
    let h_mod = h % r;
    let f = fb.points.len();
    let rel_start = Instant::now();
    let mut la_ns = 0u64;
    let (mut trials, mut found, mut direct_hits, mut collisions, mut self_collisions) =
        (0u64, 0u64, 0u64, 0u64, 0u64);
    let (mut restarts, mut table_probes) = (0u64, 0u64);
    let mut recovered = None;
    let mut verified = false;
    if f == 0 {
        return outcome(
            rel, la, ver, trials, found, matrix, recovered, verified, ctr,
        );
    }

    // Sixteen jumps with tracked coefficients.
    let jump_count = 16usize;
    let jumps: Vec<(G::Elt, u64, u64)> = (0..jump_count)
        .map(|_| {
            let c = rng.gen_range(1..r);
            let d = rng.gen_range(1..r);
            let cg = g.mul(&mut rel.group_ops, generator, c);
            let dq = g.mul(&mut rel.group_ops, target, d);
            let p = g.add(&mut rel.group_ops, cg, dq);
            (p, c, d)
        })
        .collect();

    let mut table: HashMap<u64, Stored> = HashMap::new();
    let mut fresh = true;
    let (mut point, mut a, mut b) = (g.identity(), 0u64, 0u64);
    let mut row = vec![0u64; d_col];

    let add_relation = |matrix: &mut L,
                        la: &mut PhaseCost,
                        ver: &mut PhaseCost,
                        rel_ops: &mut GroupOps,
                        row_part: Vec<u64>,
                        da: u64,
                        db: u64,
                        la_ns: &mut u64,
                        recovered: &mut Option<u64>,
                        verified: &mut bool|
     -> bool {
        let _ = rel_ops;
        let la_start = Instant::now();
        let mut full = row_part;
        full.push(submod(0, mulmod(h_mod, db, r), r));
        let rhs = mulmod(h_mod, da, r);
        if let RowStatus::Inconsistent = matrix.add_row(full, rhs) {
            la.count("inconsistent", 1);
            *la_ns += la_start.elapsed().as_nanos() as u64;
            return true;
        }
        if let Some(d) = matrix.pinned(d_col) {
            *la_ns += la_start.elapsed().as_nanos() as u64;
            let v_start = Instant::now();
            let check = g.mul(&mut ver.group_ops, generator, d);
            ver.wall_ns = v_start.elapsed().as_nanos() as u64;
            if check == target {
                *recovered = Some(d);
                *verified = true;
            } else {
                ver.count("verification_failures", 1);
            }
            return true;
        }
        *la_ns += la_start.elapsed().as_nanos() as u64;
        false
    };

    while trials < max_trials {
        if fresh {
            a = rng.gen_range(1..r);
            b = rng.gen_range(1..r);
            let ag = g.mul(&mut rel.group_ops, generator, a);
            let bq = g.mul(&mut rel.group_ops, target, b);
            point = g.add(&mut rel.group_ops, ag, bq);
            row.iter_mut().for_each(|c| *c = 0);
            fresh = false;
            restarts += 1;
        } else {
            let hsh = mix(g.key(&point));
            let j = (hsh as usize >> 8) % jump_count;
            let (jp, aj, bj) = jumps[j];
            point = g.add(&mut rel.group_ops, point, jp);
            a = addmod(a, aj, r);
            b = addmod(b, bj, r);
        }
        if g.is_identity(&point) {
            fresh = true;
            continue;
        }
        // Subtract one base point chosen by the residual itself, so that
        // equal residuals take equal steps.
        let i = (mix(g.key(&point) ^ 0xA5A5) as usize >> 4) % f;
        point = g.add(&mut rel.group_ops, point, fb.points[fb.neg_index[i]]);
        let c = fb.col_of[i];
        row[c] = addmod(row[c], fb.coef_of[i], r);
        trials += 1;
        if g.is_identity(&point) {
            // [a]G + [b]Q = Σ P_j exactly: a classic relation.
            direct_hits += 1;
            found += 1;
            if add_relation(
                matrix,
                &mut la,
                &mut ver,
                &mut rel.group_ops,
                row.clone(),
                a,
                b,
                &mut la_ns,
                &mut recovered,
                &mut verified,
            ) {
                break;
            }
            fresh = true;
            continue;
        }
        // A residual in ±F is a complete decomposition.
        if let Some(k) = fb.index_of_key(g.key(&point)) {
            direct_hits += 1;
            found += 1;
            let mut full = row.clone();
            let kc = fb.col_of[k];
            full[kc] = addmod(full[kc], fb.coef_of[k], r);
            if add_relation(
                matrix,
                &mut la,
                &mut ver,
                &mut rel.group_ops,
                full,
                a,
                b,
                &mut la_ns,
                &mut recovered,
                &mut verified,
            ) {
                break;
            }
            fresh = true;
            continue;
        }
        // Canonical sign: store the representative with the smaller key,
        // with its coefficients negated when that is −L.
        let kp = g.key(&point);
        let kn = g.key(&g.neg(point));
        let (key, sa, sb, srow): (u64, u64, u64, Vec<u64>) = if kp <= kn {
            (kp, a, b, row.clone())
        } else {
            (
                kn,
                submod(0, a, r),
                submod(0, b, r),
                row.iter().map(|&c| submod(0, c, r)).collect(),
            )
        };
        table_probes += 1;
        if let Some(prev) = table.get(&key) {
            // Same residual: (a − a')G + (b − b')Q = Σ row − Σ row'.
            let da = submod(sa, prev.a, r);
            let db = submod(sb, prev.b, r);
            let diff: Vec<u64> = srow
                .iter()
                .zip(&prev.row)
                .map(|(&x, &y)| submod(x, y, r))
                .collect();
            if da == 0 && db == 0 && diff.iter().all(|&c| c == 0) {
                self_collisions += 1;
            } else {
                collisions += 1;
                found += 1;
                if add_relation(
                    matrix,
                    &mut la,
                    &mut ver,
                    &mut rel.group_ops,
                    diff,
                    da,
                    db,
                    &mut la_ns,
                    &mut recovered,
                    &mut verified,
                ) {
                    break;
                }
            }
            fresh = true;
            continue;
        }
        if table.len() < TABLE_CAP {
            table.insert(
                key,
                Stored {
                    a: sa,
                    b: sb,
                    row: srow,
                },
            );
        }
    }
    rel.wall_ns = rel_start.elapsed().as_nanos() as u64 - la_ns - ver.wall_ns;
    rel.count("trials", trials);
    rel.count("relations", found);
    rel.count("residual_collisions", collisions);
    rel.count("residual_self_collisions", self_collisions);
    rel.count("residual_direct_hits", direct_hits);
    rel.count("residual_table_entries", table.len() as u64);
    rel.count(
        "residual_table_bytes_uncharged",
        (table.len() * (16 + 8 * d_col)) as u64,
    );
    rel.count("lookups", table_probes + trials);
    rel.count("walk_restarts", restarts);
    rel.count("walk_jumps", jump_count as u64);
    la.wall_ns = la_ns;
    let (work, work_unit) = matrix.work();
    la.count(work_unit, work);
    la.count("rows", found);
    la.count("columns", cols as u64);
    la.count("rank", matrix.rank() as u64);
    outcome(
        rel, la, ver, trials, found, matrix, recovered, verified, ctr,
    )
}

#[allow(clippy::too_many_arguments)]
fn outcome<L: RelationSolver + ?Sized>(
    relations: PhaseCost,
    linear_algebra: PhaseCost,
    verify: PhaseCost,
    trials: u64,
    found: u64,
    matrix: &mut L,
    recovered: Option<u64>,
    verified: bool,
    counters: OracleCounters,
) -> PipelineOutcome {
    PipelineOutcome {
        relations,
        linear_algebra,
        verify,
        trials,
        relations_found: found,
        independent: matrix.rank() as u64,
        dependent: matrix.dependent(),
        rank: matrix.rank() as u64,
        recovered,
        verified,
        counters,
    }
}

#[cfg(test)]
mod tests {
    use super::super::plugins::{MitmOracle, PrimeAbscissaBase};
    use super::super::stages::{InstanceCtx, Params, Targets};
    use super::super::{run_pipeline, PipelineSpec};
    use crate::cryptanalysis::ic_boundary::{
        roster_prime_instance, Calibration, CountedGroup, GroupOps,
    };

    /// **The residual collector recovers the planted logarithm** on a
    /// prime instance, with the same base and matrix the walk collector
    /// uses, and reports its own counters.
    #[test]
    fn residual_collector_recovers_the_planted_logarithm() {
        let inst = roster_prime_instance(14).unwrap();
        let g = inst.generator_point();
        let mut ops = GroupOps::default();
        let planted = 30_011 % inst.r;
        let q = inst.curve.mul(&mut ops, g, planted);
        let ctx = InstanceCtx {
            group: &inst.curve,
            generator: g,
            target: q,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: None,
        };
        let base = PrimeAbscissaBase { instance: &inst };
        let mut spec = PipelineSpec {
            factor_base: "prime-abscissa".into(),
            oracle: "mitm".into(),
            targets: Targets::Residual,
            max_trials: 4_000_000,
            seed: 11,
            ..Default::default()
        };
        spec.factor_base_params = Params::default();
        spec.factor_base_params.set("size", "24");
        spec.oracle_params.set("negation_folded", "1");
        let mut oracle = MitmOracle::new(2);
        let report = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut oracle,
            planted,
            &Calibration::default(),
            None,
        )
        .expect("the residual configuration should run");
        assert!(report.verified, "did not recover the planted logarithm");
        assert_eq!(report.recovered, Some(planted));
        let native = &report.decomposition.cost.native;
        assert!(native.get("residual_table_entries").copied().unwrap_or(0) > 0);
        assert!(native.get("trials").copied().unwrap_or(0) >= report.decomposition.relations_found);
        assert_eq!(report.spec.targets.name(), "residual");
    }
}
