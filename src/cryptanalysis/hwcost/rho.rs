//! Parallel Pollard rho with distinguished points (van Oorschot–Wiener) as
//! a [`Workload`], from the kernel's exact operation counts per iteration.
//!
//! - Expected iterations to the first collision: `√(πn/2)`, or `√(πn/4)`
//!   with the negation map (whose fruitless-cycle handling is not priced
//!   here and is listed as a note), plus `W · 2^d` for the `W` concurrent
//!   walks to reach a distinguished point after it.
//! - Distinguished-point bits `d`: the largest with `W · 2^d ≤ 1%` of the
//!   expected iterations, so the tail stays a 1% overhead.  No upper cap:
//!   a deployment tests `d` bits of a hash of `x`, not the low word the
//!   timed kernel tests.
//! - Start points: `W` scalar multiplications at about `1.5 log₂ n` group
//!   operations each, priced at the same per-iteration operations.
//! - Each distinguished point uploads `x` plus a 16-byte walk seed and an
//!   8-byte length, in messages of 1024 points over every machine's own
//!   link, and stays stored for the run.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use super::kernel::{ELEMS_MOVED_PER_ITERATION, ELEMS_RESIDENT_PER_WALK, TABLE};
use super::measure::{field_spec, KernelCounts, FIELDS};
use super::{Count, Phase, Storage, Transfer, Workload, WORKLOAD_SCHEMA};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RhoParams {
    /// A measured field (`fp64`, `w1` … `w4`).
    pub field: String,
    /// `log₂` of the prime group order `n`.
    pub group_log2: f64,
    /// Walks per core, advanced together with one batch inversion.
    pub batch: usize,
    pub machines: f64,
    pub cores_per_machine: u32,
    pub negation_map: bool,
    /// Distinguished points per upload message.
    pub dp_per_message: f64,
}

/// The smallest measured field whose prime has at least `group_log2` bits
/// (a prime-order curve over `F_p` has `n ≈ p`).
pub fn field_for(group_log2: f64) -> Result<&'static str, String> {
    // The general backends only: `fp64bf`, `fp64cm` and `pf64` are
    // comparison arms, not what a group is priced on.
    let mut fs: Vec<_> = FIELDS
        .iter()
        .filter(|f| matches!(f.name, "fp64" | "w1" | "w2" | "w3" | "w4"))
        .collect();
    fs.sort_by_key(|f| f.bits());
    fs.into_iter()
        .find(|f| f.bits() as f64 >= group_log2)
        .map(|f| f.name)
        .ok_or_else(|| format!("no measured field holds a 2^{group_log2} group"))
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RhoShape {
    pub expected_iterations_log2: f64,
    pub walks: f64,
    pub dp_bits: u32,
    pub iterations_log2: f64,
    pub distinguished_points_log2: f64,
    pub dp_record_bytes: f64,
}

pub fn shape(p: &RhoParams) -> Result<RhoShape, String> {
    let spec = field_spec(&p.field)?;
    let n_log2 = p.group_log2;
    // √(πn/c) in log₂.
    let c = if p.negation_map { 4.0 } else { 2.0 };
    let expected = 0.5 * (n_log2 + (PI / c).log2());
    let walks = p.batch as f64 * p.machines * f64::from(p.cores_per_machine);
    let d = (expected + 0.01f64.log2() - walks.log2()).floor().max(0.0) as u32;
    let tail = walks.log2() + d as f64;
    let iters = (expected.exp2() + tail.exp2()).log2();
    Ok(RhoShape {
        expected_iterations_log2: expected,
        walks,
        dp_bits: d,
        iterations_log2: iters,
        distinguished_points_log2: iters - d as f64,
        dp_record_bytes: spec.elem_bytes() as f64 + 24.0,
    })
}

/// The workload, with `counts` the kernel's exact operations per iteration
/// at this batch size ([`super::measure::counts_for`]).
pub fn workload(p: &RhoParams, counts: &KernelCounts) -> Result<(Workload, RhoShape), String> {
    if counts.field != p.field || counts.batch != p.batch {
        return Err("counts are for another field or batch size".into());
    }
    let spec = field_spec(&p.field)?;
    let sh = shape(p)?;
    let e = spec.elem_bytes() as f64;
    let f = &p.field;
    let mut ops = BTreeMap::new();
    ops.insert(format!("{f}.mul"), counts.mul);
    ops.insert(format!("{f}.mul@latency"), counts.mul_dep);
    ops.insert(format!("{f}.sqr"), counts.sqr);
    ops.insert(format!("{f}.add"), counts.add);
    ops.insert(format!("{f}.inv@latency"), counts.inv);
    ops.insert(format!("{f}.rho_residual"), 1.0);
    let walk = Phase {
        name: "walk".into(),
        iterations: Count::Log2 {
            log2: sh.iterations_log2,
        },
        ops_per_iteration: ops.clone(),
        overhead_cycles_per_iteration: 0.0,
        overhead_source: String::new(),
        bytes_per_iteration: ELEMS_MOVED_PER_ITERATION as f64 * e,
        working_set_bytes_per_core: ELEMS_RESIDENT_PER_WALK as f64 * e * p.batch as f64
            + 2.0 * e * TABLE as f64,
        parallel_fraction: 1.0,
        pipes_overlap: false,
    };
    let start = Phase {
        name: "start-points".into(),
        iterations: Count::Plain(sh.walks * 1.5 * p.group_log2),
        ..walk.clone()
    };
    let dps = Count::Log2 {
        log2: sh.distinguished_points_log2,
    };
    let dp_bytes = Count::Log2 {
        log2: sh.distinguished_points_log2 + sh.dp_record_bytes.log2(),
    };
    let wl = Workload {
        schema: WORKLOAD_SCHEMA.into(),
        name: format!(
            "rho {} n=2^{} B={} ({}×{} cores{})",
            spec.prime_label,
            p.group_log2,
            p.batch,
            p.machines,
            p.cores_per_machine,
            if p.negation_map { ", negation map" } else { "" }
        ),
        phases: vec![start, walk],
        transfers: vec![Transfer {
            name: "dp-upload".into(),
            channel: "network".into(),
            bytes: dp_bytes,
            messages: Count::Log2 { log2: dps.log2() - p.dp_per_message.log2() },
            overlaps_compute: true,
            per_machine_links: true,
        }],
        storage: vec![Storage { name: "dp-table".into(), bytes: dp_bytes, months: None }],
        egress_bytes: Count::Plain(0.0),
        notes: vec![
            format!(
                "expected 2^{:.2} iterations + {} walks × 2^{} tail; 2^{:.2} distinguished points of {} B",
                sh.expected_iterations_log2, sh.walks, sh.dp_bits, sh.distinguished_points_log2, sh.dp_record_bytes
            ),
            "server-side collision lookup and the final logarithm are not priced (negligible beside the walk)".into(),
            "start points: about 1.5 log2(n) group operations each, priced as walk iterations (approximate)".into(),
            if p.negation_map {
                "negation map: fruitless-cycle detection and escape are not priced".into()
            } else {
                "no negation map".into()
            },
        ],
    };
    Ok((wl, sh))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shape_matches_sqrt_pi_n_over_2_and_keeps_tail_under_one_percent() {
        let p = RhoParams {
            field: "w2".into(),
            group_log2: 100.0,
            batch: 1024,
            machines: 1.0,
            cores_per_machine: 4,
            negation_map: false,
            dp_per_message: 1024.0,
        };
        let s = shape(&p).unwrap();
        assert!((s.expected_iterations_log2 - (50.0 + 0.5 * (PI / 2.0).log2())).abs() < 1e-12);
        let tail = s.walks * (s.dp_bits as f64).exp2();
        assert!(tail <= 0.01 * s.expected_iterations_log2.exp2());
        assert!(
            2.0 * tail > 0.01 * s.expected_iterations_log2.exp2(),
            "d is the largest such"
        );
    }

    #[test]
    fn large_groups_get_as_many_dp_bits_as_the_one_percent_rule_allows() {
        let p = RhoParams {
            field: "w4".into(),
            group_log2: 256.0,
            batch: 1024,
            machines: 1.0,
            cores_per_machine: 4,
            negation_map: false,
            dp_per_message: 1024.0,
        };
        let s = shape(&p).unwrap();
        assert!(s.dp_bits > 64);
        assert!(s.distinguished_points_log2 < 20.0);
    }

    #[test]
    fn field_for_picks_the_smallest_field_that_holds_the_group() {
        assert_eq!(field_for(60.0).unwrap(), "fp64");
        assert_eq!(field_for(64.0).unwrap(), "w1");
        assert_eq!(field_for(100.0).unwrap(), "w2");
        assert_eq!(field_for(256.0).unwrap(), "w4");
        assert!(field_for(300.0).is_err());
    }
}
