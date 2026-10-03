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
    price_phase, BinaryGroup, BinaryInstance, Calibration, CountedGroup, GroupOps, OracleCounters,
    PhaseCost, RowStatus,
};

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
    let mut table = PhaseCost::default();
    table.group_ops = table_ops;
    table.native = oracle.setup_native();
    table.wall_ns = table_start.elapsed().as_nanos() as u64;

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
    Ok(SharedRankReport {
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
    })
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
}
