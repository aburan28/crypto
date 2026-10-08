//! Bounded planted-source witness diagnostic; not an ordinary-yield sample.
use std::time::Instant;

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, KoblitzCurve,
};
use crypto_lib::cryptanalysis::wide_groebner::{wide_groebner_decompose, TwoWordFieldStructure};
use num_bigint::BigUint;
use serde_json::json;

fn peak_rss_bytes() -> Option<u64> {
    #[cfg(any(target_os = "macos", target_os = "linux"))]
    {
        let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
        if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } != 0 {
            return None;
        }
        let usage = unsafe { usage.assume_init() };
        #[cfg(target_os = "macos")]
        return Some(usage.ru_maxrss as u64);
        #[cfg(target_os = "linux")]
        return Some((usage.ru_maxrss as u64) * 1024);
    }
    #[cfg(not(any(target_os = "macos", target_os = "linux")))]
    None
}

fn main() {
    const NODE_BUDGET: usize = 128;
    let setup = Instant::now();
    let kc = KoblitzCurve::known_n83_k0().expect("pinned K0 curve");
    let source = build_standard_subspace_factor_base(&kc, 12).expect("source subspace");
    assert_eq!(source.points.len(), 4_057);
    let planted = [0usize, 2, 4];
    let target = planted.into_iter().fold(BinaryPoint::Infinity, |acc, i| {
        kc.add(&acc, &source.points[i])
    });
    assert_ne!(target, BinaryPoint::Infinity);
    let four = BigUint::from(4u32);
    let subgroup_target = kc.mul(&target, &four);
    let table = TwoWordFieldStructure::new(83, &kc.curve.irreducible).expect("two-word table");
    let index_of = source.index_map();
    let setup_ns = setup.elapsed().as_nanos();

    let query = Instant::now();
    let (witness, stats) =
        wide_groebner_decompose(&kc, &source, &index_of, &table, &target, 3, NODE_BUDGET);
    let query_ns = query.elapsed().as_nanos();
    let replayed = witness.as_ref().is_some_and(|indices| {
        let source_sum = indices.iter().fold(BinaryPoint::Infinity, |acc, &i| {
            kc.add(&acc, &source.points[i])
        });
        source_sum == target && kc.mul(&source_sum, &four) == subgroup_target
    });
    assert!(witness.is_none() || replayed);
    let outcome = if witness.is_some() {
        "found"
    } else if stats.unsupported {
        "unsupported"
    } else if stats.exhausted {
        "exhausted"
    } else {
        "refuted"
    };
    println!(
        "{}",
        json!({
            "curve_id": kc.label(),
            "source_points": source.points.len(),
            "summands": 3,
            "variables": 119,
            "planted_indices": planted,
            "node_budget": NODE_BUDGET,
            "setup_ns": setup_ns,
            "query_ns": query_ns,
            "outcome": outcome,
            "witness": witness,
            "source_and_projection_replayed": replayed,
            "nodes": stats.nodes,
            "reductions": stats.reductions,
            "refuted_branches": stats.refuted,
            "leaves": stats.leaves,
            "max_rows": stats.max_rows,
            "max_cols": stats.max_cols,
            "peak_rss_bytes": peak_rss_bytes(),
            "claim_scope": "planted_block_diagnostic_only"
        })
    );
}
