//! Negation-map interval BSGS for the legitimate P-192 residual.
//!
//! The production plan is frozen by `EXP-SCURVE-29040c`: one shared
//! `2^24` open-addressed baby table, `m=5_800_019`, and stride
//! `2m-1=11_600_037`.  The two residue orientations are advanced in
//! lockstep and every 64-bit fingerprint hit is checked with a complete
//! native scalar multiplication before it can be returned.

use crate::cryptanalysis::bsgs_fast::BabyTable;
use crate::cryptanalysis::p192_native::{
    affine_equal, batch_to_affine, scalar_mul_u64, x_fingerprint, AffinePoint, JacobianPoint,
    NativeOpCounts, P192Arithmetic,
};
use num_bigint::BigUint;
use num_traits::Zero;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

pub const FROZEN_TABLE_BITS: u32 = 24;
pub const FROZEN_M: u64 = 5_800_019;
pub const FROZEN_STRIDE: u64 = 11_600_037;
pub const FROZEN_FINITE_BABIES: u64 = 5_800_018;
pub const FROZEN_MAX_GIANT_STEPS: u64 = 5_800_018;

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct IntervalPlan {
    pub table_bits: u32,
    pub m: u64,
    pub stride: u64,
    pub workers: usize,
    pub batch_lanes: usize,
}

impl IntervalPlan {
    pub fn frozen(workers: usize, batch_lanes: usize) -> Result<Self, String> {
        let plan = Self {
            table_bits: FROZEN_TABLE_BITS,
            m: FROZEN_M,
            stride: FROZEN_STRIDE,
            workers,
            batch_lanes,
        };
        plan.validate_frozen()?;
        Ok(plan)
    }

    pub fn reduced(table_bits: u32, m: u64, workers: usize, batch_lanes: usize) -> Self {
        Self {
            table_bits,
            m,
            stride: 2 * m - 1,
            workers: workers.max(1),
            batch_lanes: batch_lanes.max(1),
        }
    }

    pub fn validate_common(&self) -> Result<(), String> {
        if self.m < 2 || self.m > u32::MAX as u64 {
            return Err("BSGS m must be in 2..=u32::MAX".into());
        }
        if self.stride != 2 * self.m - 1 {
            return Err("negation-BSGS stride must equal 2m-1".into());
        }
        if self.table_bits < 2 || self.table_bits >= usize::BITS {
            return Err("unsupported BSGS table width".into());
        }
        if (1u128 << self.table_bits) < 2 * self.m as u128 {
            return Err("baby table would exceed one-half load".into());
        }
        if self.workers == 0 || self.workers > 64 {
            return Err("BSGS worker count must be in 1..=64".into());
        }
        if self.batch_lanes == 0 || self.batch_lanes > 65_536 {
            return Err("BSGS batch lanes must be in 1..=65536".into());
        }
        Ok(())
    }

    pub fn validate_frozen(&self) -> Result<(), String> {
        self.validate_common()?;
        if self.table_bits != FROZEN_TABLE_BITS
            || self.m != FROZEN_M
            || self.stride != FROZEN_STRIDE
        {
            return Err("the crypto-scale run must use the frozen BSGS geometry".into());
        }
        if self.workers > 4 {
            return Err("the frozen budget permits at most four BSGS workers".into());
        }
        if self.batch_lanes != 256 {
            return Err("the frozen crypto-scale run requires 256 batch lanes".into());
        }
        Ok(())
    }

    pub fn table_bytes(&self) -> u64 {
        8u64 << self.table_bits
    }
}

#[derive(Clone, Copy, Debug)]
pub struct IntervalTarget {
    pub orientation: i8,
    pub point: Option<AffinePoint>,
    pub max_k: u64,
}

/// Legitimate-curve replay required before a BSGS fingerprint hit may stop
/// the search.  This prevents the interval solver from accepting only the
/// quotient equation `[k]H=T`; it must also reconstruct `d=r+kM`, range-check
/// it, and verify `[d]G=Q`.
#[derive(Clone, Debug)]
pub struct FinalReplayContext {
    pub generator: AffinePoint,
    pub public_key: AffinePoint,
    pub plus_residue: BigUint,
    pub minus_residue: BigUint,
    pub modulus: BigUint,
    pub order: BigUint,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct IntervalSolution {
    pub orientation: i8,
    pub k: u64,
    pub giant_index: u64,
    pub baby_index: u64,
    pub baby_sign: i8,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct IntervalStats {
    pub table_bits: u32,
    pub table_bytes: u64,
    pub baby_entries: u64,
    pub baby_chain_steps: u64,
    /// Shifted giant centers visited in each lane.  Center `i` is
    /// `(m-1)+i(2m-1)`, so one position covers exactly one contiguous
    /// block of `2m-1` candidate scalars.
    pub giant_positions_plus: u64,
    pub giant_positions_minus: u64,
    pub fingerprint_hits: u64,
    pub false_fingerprint_hits: u64,
    pub exact_replays: u64,
    pub baby_operations: NativeOpCounts,
    pub giant_and_replay_operations: NativeOpCounts,
    pub operations: NativeOpCounts,
}

struct BuiltBabyTable {
    table: BabyTable,
    stats: IntervalStats,
}

fn build_baby_table(base: AffinePoint, plan: IntervalPlan) -> Result<BuiltBabyTable, String> {
    plan.validate_common()?;
    let table = BabyTable::new(plan.table_bits);
    let entries = plan.m - 1;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(plan.workers)
        .build()
        .map_err(|error| format!("failed to create BSGS worker pool: {error}"))?;

    let partials: Result<Vec<(u64, NativeOpCounts)>, String> = pool.install(|| {
        (0..plan.workers)
            .into_par_iter()
            .map(|worker| {
                let start = 1 + entries * worker as u64 / plan.workers as u64;
                let end = 1 + entries * (worker as u64 + 1) / plan.workers as u64;
                let mut arithmetic = P192Arithmetic::new();
                let mut current = scalar_mul_u64(base, start, &mut arithmetic);
                let mut inserted = 0u64;
                let mut position = start;
                while position < end {
                    let take = (end - position).min(plan.batch_lanes as u64) as usize;
                    let mut batch = Vec::with_capacity(take);
                    for _ in 0..take {
                        batch.push(current);
                        current = current.add_affine(base, &mut arithmetic);
                    }
                    let affine = batch_to_affine(&batch, &mut arithmetic);
                    for (offset, point) in affine.into_iter().enumerate() {
                        let point = point.ok_or_else(|| {
                            "unexpected Infinity in the P-192 baby interval".to_string()
                        })?;
                        let index = position + offset as u64;
                        if !table.insert(x_fingerprint(point, &arithmetic), index as u32) {
                            return Err("P-192 baby table became full".into());
                        }
                        inserted += 1;
                    }
                    position += take as u64;
                }
                Ok((inserted, arithmetic.counts))
            })
            .collect()
    });

    let mut stats = IntervalStats {
        table_bits: plan.table_bits,
        table_bytes: plan.table_bytes(),
        ..IntervalStats::default()
    };
    for (inserted, counts) in partials? {
        stats.baby_entries += inserted;
        stats.baby_chain_steps += inserted;
        stats.baby_operations.merge(counts);
        stats.operations.merge(counts);
    }
    if stats.baby_entries != entries || table.len() as u64 != entries {
        return Err(format!(
            "baby table cardinality mismatch: inserted {}, table {}, expected {entries}",
            stats.baby_entries,
            table.len()
        ));
    }
    Ok(BuiltBabyTable { table, stats })
}

fn last_giant_index(max_k: u64, plan: IntervalPlan) -> u64 {
    max_k / plan.stride
}

fn exact_candidate(
    base: AffinePoint,
    target: IntervalTarget,
    k: u64,
    replay: &FinalReplayContext,
    arithmetic: &mut P192Arithmetic,
) -> bool {
    if k > target.max_k {
        return false;
    }
    let candidate = scalar_mul_u64(base, k, arithmetic);
    let quotient_matches = match target.point {
        None => candidate.is_infinity(arithmetic),
        Some(point) => affine_equal(
            candidate,
            JacobianPoint::from_affine(point, arithmetic),
            arithmetic,
        ),
    };
    if !quotient_matches {
        return false;
    }
    let residue = if target.orientation == 1 {
        &replay.plus_residue
    } else {
        &replay.minus_residue
    };
    let d = residue + &replay.modulus * BigUint::from(k);
    if d.is_zero() || d >= replay.order {
        return false;
    }
    affine_equal(
        super::p192_native::scalar_mul(replay.generator, &d, arithmetic),
        JacobianPoint::from_affine(replay.public_key, arithmetic),
        arithmetic,
    )
}

fn probe_position(
    built: &BabyTable,
    base: AffinePoint,
    target: IntervalTarget,
    point: Option<AffinePoint>,
    giant_index: u64,
    plan: IntervalPlan,
    replay: &FinalReplayContext,
    arithmetic: &mut P192Arithmetic,
    stats: &mut IntervalStats,
) -> Option<IntervalSolution> {
    // The m-1 shift is correctness-critical.  With j in [0,m-1], this
    // center covers [i*stride,(i+1)*stride-1] without gaps.
    let center = (plan.m - 1).checked_add(giant_index.checked_mul(plan.stride)?)?;
    if point.is_none() {
        stats.exact_replays += 1;
        if exact_candidate(base, target, center, replay, arithmetic) {
            return Some(IntervalSolution {
                orientation: target.orientation,
                k: center,
                giant_index,
                baby_index: 0,
                baby_sign: 0,
            });
        }
    }
    let point = point?;
    let fingerprint = x_fingerprint(point, arithmetic);
    let mut indices = Vec::new();
    built.for_each_match(fingerprint, |index| indices.push(index));
    // Parallel insertion can permute colliding table entries.  Sort public
    // baby indices so verification replays the same candidates and operation
    // counts on every scheduler/host.
    indices.sort_unstable();
    for index in indices {
        stats.fingerprint_hits += 1;
        let baby = index as u64;
        let candidates = [
            (center.checked_add(baby), 1i8),
            (center.checked_sub(baby), -1i8),
        ];
        let mut any_replay = false;
        for (candidate, sign) in candidates {
            let Some(candidate) = candidate else {
                continue;
            };
            if candidate > target.max_k {
                continue;
            }
            any_replay = true;
            stats.exact_replays += 1;
            if exact_candidate(base, target, candidate, replay, arithmetic) {
                return Some(IntervalSolution {
                    orientation: target.orientation,
                    k: candidate,
                    giant_index,
                    baby_index: baby,
                    baby_sign: sign,
                });
            }
        }
        if any_replay {
            stats.false_fingerprint_hits += 1;
        }
    }
    None
}

/// Solve the two x-negation-compatible orientations with one baby table.
/// The returned statistics include table construction and both giant lanes.
pub fn solve_two_orientations(
    base: AffinePoint,
    plus: IntervalTarget,
    minus: IntervalTarget,
    plan: IntervalPlan,
    replay: &FinalReplayContext,
) -> Result<(IntervalSolution, IntervalStats), String> {
    plan.validate_common()?;
    if plus.orientation != 1 || minus.orientation != -1 {
        return Err("interval targets must be ordered as +1 then -1".into());
    }
    let BuiltBabyTable { table, mut stats } = build_baby_table(base, plan)?;

    let max_plus = last_giant_index(plus.max_k, plan);
    let max_minus = last_giant_index(minus.max_k, plan);
    let max_index = max_plus.max(max_minus);
    if max_plus + 1 > FROZEN_MAX_GIANT_STEPS && plan.m == FROZEN_M {
        return Err("plus orientation exceeds the frozen giant-step bound".into());
    }
    if max_minus + 1 > FROZEN_MAX_GIANT_STEPS && plan.m == FROZEN_M {
        return Err("minus orientation exceeds the frozen giant-step bound".into());
    }

    let mut arithmetic = P192Arithmetic::new();
    let stride = scalar_mul_u64(base, plan.stride, &mut arithmetic)
        .to_affine(&mut arithmetic)
        .ok_or_else(|| "frozen BSGS stride unexpectedly maps to Infinity".to_string())?;
    let negative_stride = stride.neg(&mut arithmetic);
    let center_zero = scalar_mul_u64(base, plan.m - 1, &mut arithmetic)
        .to_affine(&mut arithmetic)
        .ok_or_else(|| "initial shifted BSGS center unexpectedly maps to Infinity".to_string())?;
    let negative_center_zero = center_zero.neg(&mut arithmetic);
    let plus_unshifted = plus
        .point
        .map(|point| JacobianPoint::from_affine(point, &arithmetic))
        .unwrap_or_else(|| JacobianPoint::infinity(&arithmetic));
    let minus_unshifted = minus
        .point
        .map(|point| JacobianPoint::from_affine(point, &arithmetic))
        .unwrap_or_else(|| JacobianPoint::infinity(&arithmetic));
    let mut plus_current = plus_unshifted.add_affine(negative_center_zero, &mut arithmetic);
    let mut minus_current = minus_unshifted.add_affine(negative_center_zero, &mut arithmetic);

    let mut giant_index = 0u64;
    while giant_index <= max_index {
        let end = (giant_index + plan.batch_lanes as u64 - 1).min(max_index);
        let len = (end - giant_index + 1) as usize;
        let mut plus_batch = Vec::with_capacity(len);
        let mut minus_batch = Vec::with_capacity(len);
        for offset in 0..len {
            let index = giant_index + offset as u64;
            plus_batch.push(plus_current);
            minus_batch.push(minus_current);
            if index < max_index {
                plus_current = plus_current.add_affine(negative_stride, &mut arithmetic);
                minus_current = minus_current.add_affine(negative_stride, &mut arithmetic);
            }
        }
        let plus_affine = batch_to_affine(&plus_batch, &mut arithmetic);
        let minus_affine = batch_to_affine(&minus_batch, &mut arithmetic);
        for offset in 0..len {
            let index = giant_index + offset as u64;
            if index <= max_plus {
                stats.giant_positions_plus += 1;
                if let Some(solution) = probe_position(
                    &table,
                    base,
                    plus,
                    plus_affine[offset],
                    index,
                    plan,
                    replay,
                    &mut arithmetic,
                    &mut stats,
                ) {
                    stats.giant_and_replay_operations = arithmetic.counts;
                    stats.operations.merge(arithmetic.counts);
                    return Ok((solution, stats));
                }
            }
            if index <= max_minus {
                stats.giant_positions_minus += 1;
                if let Some(solution) = probe_position(
                    &table,
                    base,
                    minus,
                    minus_affine[offset],
                    index,
                    plan,
                    replay,
                    &mut arithmetic,
                    &mut stats,
                ) {
                    stats.giant_and_replay_operations = arithmetic.counts;
                    stats.operations.merge(arithmetic.counts);
                    return Ok((solution, stats));
                }
            }
        }
        giant_index = end + 1;
    }
    stats.giant_and_replay_operations = arithmetic.counts;
    stats.operations.merge(arithmetic.counts);
    Err("interval BSGS exhausted both frozen orientations".into())
}

/// Construct `T=Q-[r]G` and the inclusive residual interval bound.
pub fn make_target(
    orientation: i8,
    q: AffinePoint,
    g: AffinePoint,
    residue: &BigUint,
    modulus: &BigUint,
    order: &BigUint,
    arithmetic: &mut P192Arithmetic,
) -> Result<IntervalTarget, String> {
    if orientation != 1 && orientation != -1 {
        return Err("orientation must be +1 or -1".into());
    }
    if residue >= modulus || residue >= order {
        return Err("residue is outside the legitimate scalar interval".into());
    }
    let r_g = super::p192_native::scalar_mul(g, residue, arithmetic)
        .to_affine(arithmetic)
        .ok_or_else(|| "nonzero residue mapped to Infinity".to_string())?;
    let target =
        JacobianPoint::from_affine(q, arithmetic).add_affine(r_g.neg(arithmetic), arithmetic);
    let target = target.to_affine(arithmetic);
    let max_k_big = (order - BigUint::from(1u8) - residue) / modulus;
    let digits: Vec<_> = max_k_big.iter_u64_digits().collect();
    let max_k = match digits.as_slice() {
        [] => 0,
        [value] => *value,
        _ => return Err("residual interval does not fit u64".into()),
    };
    Ok(IntervalTarget {
        orientation,
        point: target,
        max_k,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::p192_native::{affine_equal, p192_order, scalar_mul};

    fn replay_context(
        generator: AffinePoint,
        public_key: AffinePoint,
        plus_residue: u64,
        minus_residue: u64,
    ) -> FinalReplayContext {
        FinalReplayContext {
            generator,
            public_key,
            plus_residue: BigUint::from(plus_residue),
            minus_residue: BigUint::from(minus_residue),
            modulus: BigUint::from(1u8),
            order: p192_order(),
        }
    }

    #[test]
    fn plan_freezes_requested_geometry() {
        let plan = IntervalPlan::frozen(4, 256).unwrap();
        assert_eq!(plan.stride, 2 * plan.m - 1);
        assert_eq!(plan.m - 1, FROZEN_FINITE_BABIES);
        assert_eq!(plan.table_bytes(), 128 * 1024 * 1024);
        assert!(IntervalPlan::frozen(5, 256).is_err());
        assert!(IntervalPlan::frozen(4, 128).is_err());
    }

    #[test]
    fn reduced_width_shared_table_recovers_either_orientation() {
        let mut arithmetic = P192Arithmetic::new();
        let g = AffinePoint::generator(&arithmetic);
        let k = 12_345u64;
        let target = scalar_mul_u64(g, k, &mut arithmetic)
            .to_affine(&mut arithmetic)
            .unwrap();
        let decoy = scalar_mul_u64(g, 50_001, &mut arithmetic)
            .to_affine(&mut arithmetic)
            .unwrap();
        let plan = IntervalPlan::reduced(8, 64, 2, 16);
        let replay = replay_context(g, target, 0, 0);
        let (solution, stats) = solve_two_orientations(
            g,
            IntervalTarget {
                orientation: 1,
                point: Some(target),
                max_k: 20_000,
            },
            IntervalTarget {
                orientation: -1,
                point: Some(decoy),
                max_k: 20_000,
            },
            plan,
            &replay,
        )
        .unwrap();
        assert_eq!(solution.orientation, 1);
        assert_eq!(solution.k, k);
        assert_eq!(stats.baby_entries, 63);
        let replay = scalar_mul(g, &BigUint::from(solution.k), &mut arithmetic);
        assert!(affine_equal(
            replay,
            JacobianPoint::from_affine(target, &arithmetic),
            &mut arithmetic
        ));
    }

    #[test]
    fn endpoint_and_infinity_special_cases_are_covered() {
        let arithmetic = P192Arithmetic::new();
        let g = AffinePoint::generator(&arithmetic);
        let plan = IntervalPlan::reduced(6, 16, 1, 8);
        let mut work = P192Arithmetic::new();
        let edge_k = 2 * plan.stride + plan.m - 1;
        let edge = scalar_mul_u64(g, edge_k, &mut work)
            .to_affine(&mut work)
            .unwrap();
        let edge_replay = replay_context(g, edge, 0, 0);
        let (edge_solution, _) = solve_two_orientations(
            g,
            IntervalTarget {
                orientation: 1,
                point: Some(edge),
                max_k: edge_k,
            },
            IntervalTarget {
                orientation: -1,
                point: Some(edge),
                max_k: edge_k,
            },
            plan,
            &edge_replay,
        )
        .unwrap();
        assert_eq!(edge_solution.k, edge_k);

        let one_public = g;
        let zero_replay = replay_context(g, one_public, 1, 1);
        let (zero_solution, _) = solve_two_orientations(
            g,
            IntervalTarget {
                orientation: 1,
                point: None,
                max_k: 0,
            },
            IntervalTarget {
                orientation: -1,
                point: Some(edge),
                max_k: 0,
            },
            plan,
            &zero_replay,
        )
        .unwrap();
        assert_eq!(zero_solution.k, 0);

        let center_k = plan.m - 1;
        let center_target = scalar_mul_u64(g, center_k, &mut work)
            .to_affine(&mut work)
            .unwrap();
        let center_public = scalar_mul_u64(g, center_k + 1, &mut work)
            .to_affine(&mut work)
            .unwrap();
        let center_replay = replay_context(g, center_public, 1, 1);
        let (center_solution, _) = solve_two_orientations(
            g,
            IntervalTarget {
                orientation: 1,
                point: Some(center_target),
                max_k: center_k,
            },
            IntervalTarget {
                orientation: -1,
                point: Some(edge),
                max_k: center_k,
            },
            plan,
            &center_replay,
        )
        .unwrap();
        assert_eq!(center_solution.k, center_k);
        assert_eq!(center_solution.baby_index, 0);
    }
}
