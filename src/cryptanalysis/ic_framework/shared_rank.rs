//! Target-blind full-rank setup for a folded Koblitz factor base.
//!
//! This is the shared part of a multi-target index-calculus run. Every
//! relation is over `[a]G`, with no public target in the probe or matrix.
//! The output includes the complete accepted-row transcript and charges
//! base construction, folded-table construction, failed probes, elimination,
//! relation checks, and one full-point check per solved base column. A later
//! descent can reuse these verified logs; this module makes no IC/rho claim.

use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use super::check_full_rank_base_columns;
use super::linalg::matrix_by_name;
use super::plugins::{CompactOrbitScanBase, FrobeniusMitmOracle};
use super::stages::{DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params};
use crate::cryptanalysis::ecbench::canonical::sha256_hex;
use crate::cryptanalysis::ic_boundary::{
    price_phase, BinaryGroup, BinaryInstance, Calibration, CountedGroup, FactorBase, GroupOps,
    OracleCounters, PhaseCost, RowStatus,
};
use crate::cryptanalysis::ic_measurement;
use crate::cryptanalysis::koblitz_fast::FastPoint;

/// Frozen target-independent setup policy. The source base and oracle are
/// the same `compact-orbit-scan` / `mitm-frobenius-counted:m=3` components
/// used in the n37 one-target sessions.
#[derive(Clone, Debug, Serialize)]
pub struct SharedRankSpec {
    pub columns: usize,
    pub raw_x_cap: u64,
    pub rank_seed: u64,
    pub max_trials: u64,
}

/// One accepted relation, sufficient for a separate implementation to
/// reconstruct `[a]G`, recheck the witness sum, and replay the matrix.
#[derive(Clone, Debug, Serialize)]
pub struct RankRelation {
    pub trial: u64,
    pub scalar: u64,
    pub point_x: u64,
    pub point_y: u64,
    pub witness: Vec<usize>,
    pub row: Vec<u64>,
    pub rhs: u64,
    pub rank_after: usize,
    pub independent: bool,
}

/// A complete target-blind setup record, including failures. `total_gae`
/// remains a lower bound when field, hash and allocation work is unpriced.
#[derive(Clone, Debug, Serialize)]
pub struct SharedRankReport {
    pub curve: String,
    pub spec: SharedRankSpec,
    pub base_sha256: String,
    pub signed_points: usize,
    pub columns: usize,
    pub trials: u64,
    pub hits: u64,
    pub rank: usize,
    pub relations: Vec<RankRelation>,
    pub column_logs: Vec<u64>,
    pub base: PhaseCost,
    pub table: PhaseCost,
    pub rank_search: PhaseCost,
    pub linear_algebra: PhaseCost,
    pub verification: PhaseCost,
    pub total_gae: f64,
    pub verified: bool,
    pub exhausted: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct SharedTargetSpec {
    pub residual_seed: u64,
    pub max_attempts: u32,
}

#[derive(Clone, Debug, Serialize)]
pub struct SharedTargetAttempt {
    pub number: u32,
    pub residual_scalar: u64,
    pub residual_x: u64,
    pub residual_y: u64,
    pub witness: Option<Vec<usize>>,
    pub row: Option<Vec<u64>>,
}

#[derive(Clone, Debug, Serialize)]
pub struct SharedTargetRow {
    pub index: usize,
    pub point_x: u64,
    pub point_y: u64,
    pub attempts: Vec<SharedTargetAttempt>,
    pub recovered_log: Option<u64>,
    pub verified: bool,
    pub query: PhaseCost,
    pub pdp: PhaseCost,
    pub relation_check: PhaseCost,
    pub descent: PhaseCost,
    pub recovery_check: PhaseCost,
    pub online_wall_ns: u64,
    pub online_gae: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct SharedTargetReport {
    pub curve: String,
    pub spec: SharedTargetSpec,
    pub rank: SharedRankReport,
    pub targets: Vec<SharedTargetRow>,
    pub points_sha256: String,
    pub total_gae: f64,
    pub verified: bool,
}

#[inline]
fn mulmod(a: u64, b: u64, r: u64) -> u64 {
    ((a as u128 * b as u128) % r as u128) as u64
}

/// Hash the exact point order and column map in a byte format independent
/// of JSON rendering. Every integer is little-endian `u64`, prefixed by the
/// domain tag; the next replay must use the same ordering and encoding.
fn base_digest(
    fb: &crate::cryptanalysis::ic_boundary::FactorBase<
        crate::cryptanalysis::koblitz_fast::FastPoint,
    >,
) -> String {
    let mut bytes = b"ic-shared-rank-base-v1\0".to_vec();
    bytes.extend_from_slice(&(fb.columns as u64).to_le_bytes());
    for (p, (&col, &coef)) in fb.points.iter().zip(fb.col_of.iter().zip(&fb.coef_of)) {
        bytes.extend_from_slice(&p.x.to_le_bytes());
        bytes.extend_from_slice(&p.y.to_le_bytes());
        bytes.push(u8::from(p.infinity));
        bytes.extend_from_slice(&(col as u64).to_le_bytes());
        bytes.extend_from_slice(&coef.to_le_bytes());
    }
    sha256_hex(&bytes)
}

/// Build and certify a reusable Koblitz rank database without looking at
/// any discrete-log target. The caller supplies a calibration; unpriced
/// native counters remain visible in the phase records.
pub fn run_shared_rank(
    inst: &BinaryInstance,
    spec: &SharedRankSpec,
    calib: &Calibration,
) -> Result<SharedRankReport, String> {
    prepare_shared_rank(inst, spec, calib).map(|(report, _, _)| report)
}

fn prepare_shared_rank<'i>(
    inst: &'i BinaryInstance,
    spec: &SharedRankSpec,
    calib: &Calibration,
) -> Result<
    (
        SharedRankReport,
        FactorBase<FastPoint>,
        FrobeniusMitmOracle<'i>,
    ),
    String,
> {
    if inst.koblitz.is_none() || inst.r < 3 || spec.columns == 0 || spec.max_trials == 0 {
        return Err("shared rank needs a Koblitz instance, positive columns and trial cap".into());
    }
    let group = BinaryGroup(&inst.fast);
    let ctx = InstanceCtx {
        group: &group,
        generator: inst.generator,
        // The builder and folded oracle are target independent. This
        // placeholder is never used to construct a probe or relation.
        target: group.identity(),
        r: inst.r,
        cofactor: inst.cofactor,
        group_order: inst.group_order,
        name: inst.name.clone(),
        field_degree: Some(inst.n),
    };
    let mut base_params = Params::default();
    base_params.set("columns", spec.columns.to_string());
    base_params.set("raw_x_cap", spec.raw_x_cap.to_string());
    let builder = CompactOrbitScanBase { instance: inst };
    let mut base_ops = GroupOps::default();
    let base_start = Instant::now();
    let fb = builder.build(&ctx, &base_params, &mut base_ops)?;
    if fb.columns != spec.columns {
        return Err(format!(
            "base has {} columns, requested {}",
            fb.columns, spec.columns
        ));
    }
    let mut base = fb.cost.clone();
    base.group_ops.merge(base_ops);
    base.wall_ns = base_start.elapsed().as_nanos() as u64;
    let base_sha256 = base_digest(&fb);

    let mut oracle = FrobeniusMitmOracle::new_counted(3, inst);
    let table_start = Instant::now();
    let mut table_ops = GroupOps::default();
    oracle.prepare(&ctx, &fb, &Params::default(), &mut table_ops)?;
    let mut table = PhaseCost {
        group_ops: table_ops,
        native: oracle.setup_native(),
        wall_ns: table_start.elapsed().as_nanos() as u64,
        ..PhaseCost::default()
    };

    let mut matrix = matrix_by_name("incremental-gauss", fb.columns, inst.r)?;
    let mut rng = StdRng::seed_from_u64(spec.rank_seed ^ 0x5348_4152_4544_524b);
    let mut search = PhaseCost::default();
    let mut la = PhaseCost::default();
    let mut verification = PhaseCost::default();
    let mut counters = OracleCounters::default();
    let mut relations = Vec::with_capacity(fb.columns.saturating_mul(2));
    let mut hits = 0u64;
    let mut trials = 0u64;
    let search_start = Instant::now();

    while trials < spec.max_trials && matrix.rank() < fb.columns {
        let trial = trials;
        trials += 1;
        let scalar = rng.gen_range(1..inst.r);
        let point = group.mul(&mut search.group_ops, inst.generator, scalar);
        let Some(witness) =
            oracle.decompose(&ctx, &fb, &mut search.group_ops, &mut counters, point)
        else {
            continue;
        };
        hits += 1;
        let verify_start = Instant::now();
        let mut sum = group.identity();
        for &i in &witness {
            let Some(&p) = fb.points.get(i) else {
                return Err(format!(
                    "oracle witness index {i} outside base at trial {trial}"
                ));
            };
            sum = group.add(&mut verification.group_ops, sum, p);
        }
        verification.wall_ns += verify_start.elapsed().as_nanos() as u64;
        verification.count("relation_witnesses_checked", 1);
        if sum != point {
            verification.count("relation_witness_failures", 1);
            return Err(format!(
                "oracle witness does not sum to [a]G at trial {trial}"
            ));
        }
        let mut row = vec![0u64; fb.columns];
        for &i in &witness {
            let col = fb.col_of[i];
            row[col] = ((row[col] as u128 + fb.coef_of[i] as u128) % inst.r as u128) as u64;
        }
        let rhs = mulmod(inst.cofactor % inst.r, scalar, inst.r);
        let la_start = Instant::now();
        let status = matrix.add_row(row.clone(), rhs);
        la.wall_ns += la_start.elapsed().as_nanos() as u64;
        if status == RowStatus::Inconsistent {
            la.count("inconsistent_rows", 1);
            return Err(format!("inconsistent relation at trial {trial}"));
        }
        relations.push(RankRelation {
            trial,
            scalar,
            point_x: point.x,
            point_y: point.y,
            witness,
            row,
            rhs,
            rank_after: matrix.rank(),
            independent: status == RowStatus::Independent,
        });
    }

    search.wall_ns = search_start
        .elapsed()
        .as_nanos()
        .saturating_sub(la.wall_ns as u128 + verification.wall_ns as u128)
        as u64;
    search.count("trials", trials);
    search.count("hits", hits);
    search.count("lookups", counters.lookups);
    search.count("canonicalisations", counters.canonicalisations);
    search.count("frobfold_mismatches", counters.frobfold_mismatches);
    search.count("lift_failures", counters.lift_failures);
    la.count("rows", hits);
    la.count("rank", matrix.rank() as u64);
    la.count("dependent_rows", matrix.dependent());
    let (work, unit) = matrix.work();
    la.count(unit, work);

    let full_rank = matrix.rank() == fb.columns;
    let verified =
        full_rank && check_full_rank_base_columns(&ctx, &fb, matrix.as_ref(), &mut verification);
    let column_logs = if verified {
        (0..fb.columns)
            .map(|col| matrix.pinned(col).expect("full rank pins every column"))
            .collect()
    } else {
        Vec::new()
    };
    for phase in [
        &mut base,
        &mut table,
        &mut search,
        &mut la,
        &mut verification,
    ] {
        price_phase(phase, calib);
    }
    let total_gae = base.gae + table.gae + search.gae + la.gae + verification.gae;
    Ok((
        SharedRankReport {
            curve: inst.name.clone(),
            spec: spec.clone(),
            base_sha256,
            signed_points: fb.points.len(),
            columns: fb.columns,
            trials,
            hits,
            rank: matrix.rank(),
            relations,
            column_logs,
            base,
            table,
            rank_search: search,
            linear_algebra: la,
            verification,
            total_gae,
            verified,
            exhausted: !full_rank,
        },
        fb,
        oracle,
    ))
}

/// Recover point-only targets using the target-blind rank setup and its
/// already prepared folded table. The input point digest is supplied by the
/// caller after reading the point-only file; no fixture scalar is read here.
pub fn run_shared_rank_targets(
    inst: &BinaryInstance,
    rank_spec: &SharedRankSpec,
    target_spec: &SharedTargetSpec,
    points: &[FastPoint],
    points_sha256: String,
    calib: &Calibration,
) -> Result<SharedTargetReport, String> {
    if target_spec.max_attempts == 0 || points.is_empty() {
        return Err("target gate needs points and a positive attempt cap".into());
    }
    let profile_online = std::env::var("ECBENCH_CALLGRIND_TARGET").as_deref() == Ok("1");
    if profile_online && points.len() != 1 {
        return Err("target-only Callgrind profile requires exactly one target".into());
    }
    let (rank, fb, mut oracle) = prepare_shared_rank(inst, rank_spec, calib)?;
    if !rank.verified || rank.column_logs.len() != fb.columns {
        return Err("target-blind rank setup did not verify every base column".into());
    }
    let group = BinaryGroup(&inst.fast);
    let inverse_h = modpow(inst.cofactor % inst.r, inst.r - 2, inst.r);
    if mulmod(inverse_h, inst.cofactor % inst.r, inst.r) != 1 {
        return Err("subgroup cofactor is not invertible".into());
    }
    let mut targets = Vec::with_capacity(points.len());
    let mut total_gae = rank.total_gae;
    for (index, &target) in points.iter().enumerate() {
        let ctx = InstanceCtx {
            group: &group,
            generator: inst.generator,
            target,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: Some(inst.n),
        };
        if profile_online {
            ic_measurement::callgrind_dump(b"ecbench_before_online\0");
        }
        let started = Instant::now();
        let mut query = PhaseCost::default();
        let mut pdp = PhaseCost::default();
        let mut relation_check = PhaseCost::default();
        let mut descent = PhaseCost::default();
        let mut recovery_check = PhaseCost::default();
        let mut counters = OracleCounters::default();
        let mut attempts = Vec::new();
        let mut recovered_log = None;
        let mut verified = false;
        let mut rng = StdRng::seed_from_u64(target_spec.residual_seed ^ index as u64);
        let validation_start = Instant::now();
        let valid = !target.infinity
            && inst.fast.is_on_curve(target)
            && group.is_identity(&group.mul(&mut query.group_ops, target, inst.r));
        query.wall_ns += validation_start.elapsed().as_nanos() as u64;
        if !valid {
            return Err(format!("target {index} is not a valid subgroup point"));
        }
        for number in 0..target_spec.max_attempts {
            let query_start = Instant::now();
            let a = if number == 0 {
                0
            } else {
                rng.gen_range(1..inst.r)
            };
            let residual = if a == 0 {
                target
            } else {
                let shift = group.mul(&mut query.group_ops, inst.generator, a);
                group.add(&mut query.group_ops, target, shift)
            };
            query.wall_ns += query_start.elapsed().as_nanos() as u64;
            query.count("attempts", 1);
            let pdp_start = Instant::now();
            let witness = oracle.decompose(&ctx, &fb, &mut pdp.group_ops, &mut counters, residual);
            pdp.wall_ns += pdp_start.elapsed().as_nanos() as u64;
            let mut row = None;
            if let Some(indices) = &witness {
                let check_start = Instant::now();
                let mut sum = group.identity();
                let mut dense = vec![0u64; fb.columns];
                for &i in indices {
                    let point = *fb
                        .points
                        .get(i)
                        .ok_or("target witness index outside base")?;
                    sum = group.add(&mut relation_check.group_ops, sum, point);
                    let column = fb.col_of[i];
                    dense[column] =
                        ((dense[column] as u128 + fb.coef_of[i] as u128) % inst.r as u128) as u64;
                }
                relation_check.wall_ns += check_start.elapsed().as_nanos() as u64;
                relation_check.count("witnesses_checked", 1);
                if sum != residual {
                    return Err(format!("target {index} witness does not sum to residual"));
                }
                let descent_start = Instant::now();
                let mut projected_log = 0u64;
                for (&coefficient, &column_log) in dense.iter().zip(&rank.column_logs) {
                    projected_log = ((projected_log as u128
                        + mulmod(coefficient, column_log, inst.r) as u128)
                        % inst.r as u128) as u64;
                }
                let log = (mulmod(projected_log, inverse_h, inst.r) + inst.r - a) % inst.r;
                descent.wall_ns += descent_start.elapsed().as_nanos() as u64;
                let replay_start = Instant::now();
                verified = group.mul(&mut recovery_check.group_ops, inst.generator, log) == target;
                recovery_check.wall_ns += replay_start.elapsed().as_nanos() as u64;
                recovery_check.count("scalar_replays", 1);
                if !verified {
                    return Err(format!(
                        "target {index} recovered scalar failed full-point check"
                    ));
                }
                recovered_log = Some(log);
                row = Some(dense);
            }
            attempts.push(SharedTargetAttempt {
                number,
                residual_scalar: a,
                residual_x: residual.x,
                residual_y: residual.y,
                witness,
                row,
            });
            if verified {
                break;
            }
        }
        pdp.count("lookups", counters.lookups);
        pdp.count("canonicalisations", counters.canonicalisations);
        pdp.count("frobfold_mismatches", counters.frobfold_mismatches);
        pdp.count("lift_failures", counters.lift_failures);
        let mut online_wall_ns = started.elapsed().as_nanos() as u64;
        if profile_online {
            ic_measurement::callgrind_dump(b"ecbench_online\0");
        }
        let timed_phases = query.wall_ns
            + pdp.wall_ns
            + relation_check.wall_ns
            + descent.wall_ns
            + recovery_check.wall_ns;
        if timed_phases <= online_wall_ns {
            query.wall_ns += online_wall_ns - timed_phases;
        } else {
            online_wall_ns = timed_phases;
        }
        for phase in [
            &mut query,
            &mut pdp,
            &mut relation_check,
            &mut descent,
            &mut recovery_check,
        ] {
            price_phase(phase, calib);
        }
        let online_gae =
            query.gae + pdp.gae + relation_check.gae + descent.gae + recovery_check.gae;
        total_gae += online_gae;
        targets.push(SharedTargetRow {
            index,
            point_x: target.x,
            point_y: target.y,
            attempts,
            recovered_log,
            verified,
            query,
            pdp,
            relation_check,
            descent,
            recovery_check,
            online_wall_ns,
            online_gae,
        });
    }
    let verified = targets.iter().all(|row| row.verified);
    Ok(SharedTargetReport {
        curve: inst.name.clone(),
        spec: target_spec.clone(),
        rank,
        targets,
        points_sha256,
        total_gae,
        verified,
    })
}

fn modpow(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = mulmod(result, base, modulus);
        }
        base = mulmod(base, base, modulus);
        exponent >>= 1;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::koblitz_instance;

    #[test]
    fn small_koblitz_rank_is_target_blind_and_point_checked() {
        let inst = koblitz_instance(1, 17).expect("test curve");
        let spec = SharedRankSpec {
            columns: 4,
            raw_x_cap: 10_000,
            rank_seed: 7,
            max_trials: 10_000,
        };
        let report = run_shared_rank(&inst, &spec, &Calibration::default()).unwrap();
        assert!(
            report.verified,
            "rank {} after {} trials",
            report.rank, report.trials
        );
        assert_eq!(report.rank, spec.columns);
        assert_eq!(report.column_logs.len(), spec.columns);
        assert_eq!(
            report.verification.get("base_log_columns_checked"),
            spec.columns as u64
        );
        assert_eq!(report.verification.get("relation_witness_failures"), 0);
        assert_eq!(report.relations.len() as u64, report.hits);
    }

    #[test]
    fn shared_rank_recovers_point_only_targets_with_full_point_replay() {
        let inst = koblitz_instance(1, 17).expect("test curve");
        let rank_spec = SharedRankSpec {
            columns: 4,
            raw_x_cap: 10_000,
            rank_seed: 7,
            max_trials: 10_000,
        };
        let target_spec = SharedTargetSpec {
            residual_seed: 19,
            max_attempts: 64,
        };
        let scalars = [7u64, 113, 257];
        let points: Vec<_> = scalars
            .iter()
            .map(|&scalar| inst.fast.mul_u64(inst.generator, scalar))
            .collect();
        let report = run_shared_rank_targets(
            &inst,
            &rank_spec,
            &target_spec,
            &points,
            "test-point-only-file".into(),
            &Calibration::default(),
        )
        .unwrap();
        assert!(report.verified);
        for (target, scalar) in report.targets.iter().zip(scalars) {
            assert_eq!(target.recovered_log, Some(scalar % inst.r));
            assert!(!target.attempts.is_empty());
            assert_eq!(target.recovery_check.get("scalar_replays"), 1);
        }
    }
}
