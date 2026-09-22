//! Complete the selected public n=59 logarithm system from a uniform
//! relation prefix plus a forced-column pair-lookup tail.

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, cofactor_project_factor_base, CollectedRelation,
    DecompositionStrategy, KoblitzCurve, KoblitzIcOptions, LinearAlgebra, PairSumTable,
    ProjectedFactorBase, RelationCollector, RelationWorkUnit,
};
use crypto_lib::cryptanalysis::koblitz_sparse_la::SparseSolveOptions;
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn process_resources() -> (f64, u64) {
    unsafe {
        let mut usage: libc::rusage = std::mem::zeroed();
        if libc::getrusage(libc::RUSAGE_SELF, &mut usage) != 0 {
            return (0.0, 0);
        }
        let seconds = |time: libc::timeval| {
            time.tv_sec as f64 + time.tv_usec as f64 / 1_000_000.0
        };
        #[cfg(target_os = "macos")]
        let rss_multiplier = 1u64;
        #[cfg(not(target_os = "macos"))]
        let rss_multiplier = 1024u64;
        (
            seconds(usage.ru_utime) + seconds(usage.ru_stime),
            usage.ru_maxrss.max(0) as u64 * rss_multiplier,
        )
    }
}

fn read_prefix(directory: &Path, units: usize) -> Result<Vec<CollectedRelation>, String> {
    let mut relations = Vec::new();
    for unit in 0..units {
        let path = directory.join(format!("unit-{unit:05}.json"));
        let value: Value = serde_json::from_slice(&fs::read(&path).map_err(|e| e.to_string())?)
            .map_err(|e| e.to_string())?;
        let mut rows: Vec<CollectedRelation> =
            serde_json::from_value(value["relations"].clone()).map_err(|e| e.to_string())?;
        relations.append(&mut rows);
    }
    Ok(relations)
}

fn decimal(value: &str) -> BigUint {
    BigUint::parse_bytes(value.as_bytes(), 10).expect("decimal integer")
}

fn main() -> Result<(), String> {
    let relation_directory = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(
                "/private/tmp/ic-n59-projected-l15-witnessed-cached-full/run/relations",
            )
        });
    let prefix_units = std::env::args()
        .nth(2)
        .and_then(|value| value.parse().ok())
        .unwrap_or(25usize);
    let chunk_trials = std::env::args()
        .nth(3)
        .and_then(|value| value.parse().ok())
        .unwrap_or(100_000u64);
    let max_chunks = std::env::args()
        .nth(4)
        .and_then(|value| value.parse().ok())
        .unwrap_or(10usize);

    let whole_begin = Instant::now();
    let selection_begin = Instant::now();
    let kc = KoblitzCurve::new(1, 59).ok_or("K_1/F_2^59 unavailable")?;
    let parent = build_standard_subspace_factor_base(&kc, 15)?;
    let fb = cofactor_project_factor_base(&kc, &parent)?;
    let projected = ProjectedFactorBase::new(&kc, &fb);
    let selection_seconds = selection_begin.elapsed().as_secs_f64();

    let opts = KoblitzIcOptions {
        m: 3,
        descent_m: Some(2),
        collection_window: Some(1024),
        strategy: DecompositionStrategy::PairTable,
        collapse_negation: true,
        collapse_projected_orbits: true,
        allow_direct_relation: false,
        linear_algebra: LinearAlgebra::Sparse(SparseSolveOptions::default()),
        ..KoblitzIcOptions::default()
    };

    let pair_begin = Instant::now();
    let pair = PairSumTable::build_witnessed_compact_within(&kc, &fb, 8u128 << 30)
        .ok_or("witnessed compact table did not fit")?;
    let pair_seconds = pair_begin.elapsed().as_secs_f64();

    let prefix = read_prefix(&relation_directory, prefix_units)?;
    let prefix_begin = Instant::now();
    let mut solver = projected.log_solver(&opts).ok_or("projected log solver")?;
    solver.push(&prefix);
    let prefix_verification_seconds = prefix_begin.elapsed().as_secs_f64();
    if solver.try_solve().is_some() {
        return Err("uniform prefix unexpectedly solved before the targeted tail".into());
    }
    let initially_uncovered = solver.uncovered_columns();
    if initially_uncovered.len() != 1 {
        return Err(format!(
            "expected one uncovered prefix column, found {:?}",
            initially_uncovered
        ));
    }
    let fixed_column = initially_uncovered[0];
    let fixed_point_index = projected
        .factor_point_for_column(fixed_column)
        .ok_or("uncovered column has no factor-base point")?;
    let collector = RelationCollector::with_pair_table(&kc, &fb, &opts, Some(&pair))
        .ok_or("targeted relation collector")?;

    let targeted_begin = Instant::now();
    let targeted_seed = 0x5441_5247_4554_5901u64;
    let mut targeted_relations = Vec::new();
    let mut targeted_trials = 0usize;
    let mut pair_lookups = 0u64;
    let mut targeted_collection_seconds = 0.0f64;
    let mut solved = None;
    let mut chunks_used = 0usize;
    for chunk in 0..max_chunks {
        let (relations, report) = collector
            .collect_with_forced_point(
                RelationWorkUnit {
                    seed: targeted_seed,
                    start: chunk as u64 * chunk_trials,
                    count: chunk_trials,
                },
                fixed_point_index,
            )
            .ok_or("targeted collection unsupported")?;
        targeted_trials += report.trials;
        pair_lookups += report.pair_lookups;
        targeted_collection_seconds += report.elapsed_seconds;
        targeted_relations.extend(relations.iter().cloned());
        solver.push(&relations);
        chunks_used = chunk + 1;
        if let Some(result) = solver.try_solve() {
            solved = Some(result);
            break;
        }
    }
    let targeted_seconds = targeted_begin.elapsed().as_secs_f64();
    let (table, log_report) = solved.ok_or("targeted tail did not determine every column")?;
    if !table.verify(&kc) {
        return Err("targeted log table failed group certification".into());
    }
    let remaining_uncovered = solver.uncovered_columns();
    drop(solver);

    let target = BinaryPoint::Affine {
        x: F2mElement::from_biguint(&decimal("551007654412192707"), 59),
        y: F2mElement::from_biguint(&decimal("76595654768815112"), 59),
    };
    let descent_begin = Instant::now();
    let individual = projected
        .individual_log_solver(&table, &opts, Some(&pair))
        .ok_or("individual log solver")?;
    let (recovered, descent_report) = individual.solve(&target).ok_or("target descent failed")?;
    let descent_seconds = descent_begin.elapsed().as_secs_f64();
    if kc.mul(kc.generator(), &recovered) != target {
        return Err("recovered scalar failed [d]G = Q".into());
    }
    let (core_seconds, peak_rss_bytes) = process_resources();

    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "schema_version": 1,
            "operation": "koblitz_n59_targeted_relation_tail",
            "factor_base": {
                "points": fb.points.len(),
                "projected_columns": projected.columns(),
                "logs_known_by_construction": false,
                "target_subgroup_enumerated": false,
            },
            "pair_table": {"tier": pair.tier(), "entries": pair.len()},
            "uniform_prefix": {
                "units": prefix_units,
                "relations": prefix.len(),
                "verification_seconds": prefix_verification_seconds,
                "uncovered_columns": initially_uncovered,
            },
            "targeted_tail": {
                "seed": targeted_seed,
                "fixed_column": fixed_column,
                "fixed_point_index": fixed_point_index,
                "chunks_used": chunks_used,
                "trials": targeted_trials,
                "pair_lookups": pair_lookups,
                "relations": targeted_relations,
                "collection_seconds": targeted_collection_seconds,
                "elapsed_seconds": targeted_seconds,
                "remaining_uncovered_columns": remaining_uncovered,
            },
            "linear_algebra": {
                "relations": log_report.relations,
                "attempts": log_report.solve_attempts,
                "seconds": log_report.linear_algebra_seconds,
                "sparse": log_report.sparse_report,
                "verified": log_report.verified,
            },
            "solution": {
                "recovered_scalar": recovered.to_string(),
                "verified": true,
                "descent_trials": descent_report.trials,
                "descent_seconds": descent_seconds,
            },
            "timing": {
                "selection_seconds": selection_seconds,
                "pair_table_seconds": pair_seconds,
                "whole_seconds": whole_begin.elapsed().as_secs_f64(),
            },
            "resources": {
                "total_core_seconds": core_seconds,
                "peak_rss_bytes": peak_rss_bytes,
            },
        }))
        .map_err(|e| e.to_string())?
    );
    Ok(())
}
