//! Fully charged unknown-scalar run over a GGMP linearised-polynomial factor base.
//!
//! The factor base is selected from a factor of `T^n-1`; this program never
//! constructs factor-base logarithms or a subgroup log table.  The scalar is
//! used only by the harness to create a public target and verify the answer.

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_groebner::SolverEngine;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_frobenius_factor_base, factor_x_n_minus_1, koblitz_index_calculus_dlp_with_factor_base,
    DecompositionStrategy, FactorBaseDomain, KoblitzCurve, KoblitzIcOptions,
    SatDecompositionOptions,
};
use crypto_lib::cryptanalysis::pollard_rho::{pollard_rho_dlp, RhoOptions};
use crypto_lib::cryptanalysis::semaev_sat::XorEncoding;
use num_bigint::BigUint;
use serde_json::json;
use std::time::Instant;

fn point_json(point: &BinaryPoint) -> serde_json::Value {
    match point {
        BinaryPoint::Infinity => serde_json::Value::Null,
        BinaryPoint::Affine { x, y } => json!({
            "x":x.to_biguint().to_string(),
            "y":y.to_biguint().to_string()
        }),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(
        args.len(),
        9,
        "usage: koblitz_algebraic_e2e <ic|rho> <n> <a> <m> <secret> <seed> <conflict-budget> <max-trials>"
    );
    let mode = args[1].as_str();
    let n: u32 = args[2].parse().expect("n");
    let a: u8 = args[3].parse().expect("a");
    let m: usize = args[4].parse().expect("m");
    let secret: u64 = args[5].parse().expect("secret");
    let seed: u64 = args[6].parse().expect("seed");
    let conflict_budget: u64 = args[7].parse().expect("conflict budget");
    let max_trials: usize = args[8].parse().expect("max trials");
    assert!(matches!(mode, "ic" | "rho"));

    let total = Instant::now();
    let curve_start = Instant::now();
    let curve = KoblitzCurve::new(a, n).expect("usable prime-order subgroup");
    let curve_ns = curve_start.elapsed().as_nanos();
    assert!(BigUint::from(secret) < curve.subgroup_order && secret > 0);
    let target_start = Instant::now();
    let target = curve.mul(curve.generator(), &BigUint::from(secret));
    let target_ns = target_start.elapsed().as_nanos();

    if mode == "rho" {
        let rho_start = Instant::now();
        let solution = pollard_rho_dlp(
            curve.generator(),
            &target,
            &curve.subgroup_order,
            |left, right| curve.add(left, right),
            |left, right| left == right,
            |point| match point {
                BinaryPoint::Infinity => 0,
                BinaryPoint::Affine { x, .. } => {
                    (x.to_biguint().iter_u64_digits().next().unwrap_or(0) % 3) as u8
                }
            },
            |point, scalar| curve.mul(point, scalar),
            &RhoOptions {
                max_iterations: 1u64 << 28,
                max_restarts: 64,
                seed: Some(seed),
            },
        )
        .expect("generic Pollard rho run");
        let rho_ns = rho_start.elapsed().as_nanos();
        let verified = solution.x == BigUint::from(secret)
            && curve.mul(curve.generator(), &solution.x) == target;
        println!(
            "{}",
            serde_json::to_string_pretty(&json!({
                "kind":"koblitz_unknown_scalar_generic_rho_control",
                "mode":"rho",
                "scope":"same-curve generic Pollard rho control; not automorphism optimized",
                "n":n,"a":a,"subgroup_order":curve.subgroup_order.to_string(),
                "target":point_json(&target),"seed":seed,
                "iterations":solution.iterations,"verified_unknown_scalar_recovery":verified,
                "automorphism_optimized":false,
                "timing_ns":{"curve_construction":curve_ns,"target_construction":target_ns,"rho":rho_ns,"end_to_end":total.elapsed().as_nanos()}
            })).unwrap()
        );
        assert!(verified);
        return;
    }

    let factor_base_start = Instant::now();
    let factor_index = 0usize;
    let factor = factor_x_n_minus_1(n)[factor_index];
    let factor_exponents: Vec<u32> = (0..64).filter(|i| (factor >> i) & 1 == 1).collect();
    let factor_base = build_frobenius_factor_base(&curve, factor_index)
        .expect("GGMP linearised-polynomial factor base");
    let factor_base_ns = factor_base_start.elapsed().as_nanos();
    assert_eq!(factor_base.domain, FactorBaseDomain::LinearSubspace);

    let options = KoblitzIcOptions {
        m,
        factor_index,
        extra_relations: 2,
        max_trials,
        seed,
        strategy: DecompositionStrategy::Sat,
        engine: SolverEngine::default(),
        node_budget: 0,
        max_models: 64,
        sat_macaulay_degree: None,
        sat_options: SatDecompositionOptions {
            encoding: XorEncoding::Native,
            branch_on_summands: true,
            restrict_to_factor_base: true,
            trace_constraint: true,
            conflict_budget,
        },
        collapse_negation: true,
        stop_on_verified_rank: true,
        relation_batch_size: 1,
    };
    let solve_start = Instant::now();
    let report =
        koblitz_index_calculus_dlp_with_factor_base(&curve, &target, &factor_base, &options)
            .expect("index-calculus report");
    let solve_ns = solve_start.elapsed().as_nanos();
    let verified = report.log.as_ref() == Some(&BigUint::from(secret))
        && report
            .log
            .as_ref()
            .is_some_and(|value| curve.mul(curve.generator(), value) == target);
    let output = json!({
        "kind":"koblitz_unknown_scalar_algebraic_factor_base_e2e",
        "mode":"ic",
        "scope":"complete toy index-calculus pipeline over a factor base with no constructed point logs",
        "n":n,"a":a,"m":m,"subgroup_order":curve.subgroup_order.to_string(),
        "target":point_json(&target),"seed":seed,
        "factor_base_predicate":{
            "kind":"ggmp_linearised_kernel",
            "factor_index":factor_index,
            "factor_bitmask":factor,
            "linearised_exponents":factor_exponents,
            "enumerates_target_subgroup":false,
            "uses_discrete_log_labels":false,
            "factor_base_logs_constructed":false
        },
        "factor_base":{"ell":factor_base.ell,"points":factor_base.points.len(),"signed_frobenius_orbits":factor_base.unknowns()},
        "options":{"conflict_budget_per_target":conflict_budget,"max_trials":max_trials,"max_models":64,"parallel_threads":1},
        "report":{
            "relations":report.relations,"trials":report.trials,"relation_batches":report.relation_batches,
            "linear_solve_attempts":report.linear_solve_attempts,"sat_calls":report.sat_calls,
            "sat_models":report.sat_models,"sat_refutations":report.sat_refutations,
            "sat_unknowns":report.sat_unknowns,"sat_invalid_models":report.sat_invalid_models,
            "sat_conflicts":report.sat_conflicts,"direct_relation":report.direct_relation,
            "verified_unknown_scalar_recovery":verified
        },
        "timing_ns":{
            "curve_and_subgroup_construction":curve_ns,"target_construction":target_ns,
            "factor_base_predicate_and_materialisation":factor_base_ns,
            "relation_collection":report.relation_collection_ns,"linear_algebra":report.linear_algebra_ns,
            "driver_solve":solve_ns,"end_to_end":total.elapsed().as_nanos()
        },
        "interpretation":"A failed or resource-capped run is operationally inconclusive; success is accepted only after the recovered scalar reproduces the public target point"
    });
    println!("{}", serde_json::to_string_pretty(&output).unwrap());
    assert!(
        !report.direct_relation,
        "direct relation bypassed the matrix"
    );
    assert_eq!(report.sat_invalid_models, 0);
    assert!(verified, "unknown scalar was not recovered and verified");
}
