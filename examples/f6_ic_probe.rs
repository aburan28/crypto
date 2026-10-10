//! Same-binary timing probe for the PDP3 Gröbner decomposers: inherited F4
//! against F6-IC on a standard-subspace base of one Koblitz curve.
//!
//! ```bash
//! cargo run --release --example f6_ic_probe -- [--a 0] [--n 17] [--dim 6] [--m 3] \
//!     [--targets 8] [--budget 20000] [--engines f4,f6] [--reps 2]
//! ```
//!
//! Targets are `[k]G` for `k = 1 …`.  Each arm decides every target with
//! the same node budget; the probe prints one JSON line per arm and target
//! (wall, verdict, reductions, splits, geometric counters, word XORs) and a
//! per-arm total, and checks that the arms agree on every verdict.  It is a
//! development probe: wall times on a shared host are indicative only.

use crypto_lib::cryptanalysis::koblitz_groebner::{
    f4_word_ops_thread, FieldStructure, SolveStats, SolverEngine,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, groebner_decompose, groebner_decompose_f6_ic, KoblitzCurve,
};
use num_bigint::BigUint;
use std::time::Instant;

fn flag(name: &str) -> Option<String> {
    let args: Vec<String> = std::env::args().collect();
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1).cloned())
}

fn main() {
    let a: u8 = flag("--a").and_then(|v| v.parse().ok()).unwrap_or(0);
    let n: u32 = flag("--n").and_then(|v| v.parse().ok()).unwrap_or(17);
    let dim: u32 = flag("--dim").and_then(|v| v.parse().ok()).unwrap_or(6);
    let m: usize = flag("--m").and_then(|v| v.parse().ok()).unwrap_or(3);
    let targets: u64 = flag("--targets").and_then(|v| v.parse().ok()).unwrap_or(8);
    let budget: usize = flag("--budget")
        .and_then(|v| v.parse().ok())
        .unwrap_or(20_000);
    let reps: usize = flag("--reps").and_then(|v| v.parse().ok()).unwrap_or(1);
    let engines = flag("--engines").unwrap_or_else(|| "f4,f6".to_string());
    let engines: Vec<&str> = engines.split(',').collect();

    let kc = KoblitzCurve::new(a, n).expect("curve has a usable subgroup");
    let fb = build_standard_subspace_factor_base(&kc, dim).expect("factor base");
    let index_of = fb.index_map();
    let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
    let engine = SolverEngine::InheritedF4 { max_degree: 3 };
    eprintln!(
        "K({a}, {n}), standard subspace dim {dim}: {} base points, m = {m}, {targets} targets, budget {budget}",
        fb.points.len()
    );

    let mut verdicts: Vec<Vec<bool>> = Vec::new();
    for arm in &engines {
        let mut arm_verdicts = Vec::new();
        let mut total_ns = 0u128;
        let mut total = SolveStats::default();
        let mut total_xors = 0u64;
        for k in 1..=targets {
            let target = kc.mul(kc.generator(), &BigUint::from(k));
            let mut best_ns = u128::MAX;
            let mut last: Option<(Option<Vec<usize>>, SolveStats, u64)> = None;
            for _ in 0..reps {
                let x0 = f4_word_ops_thread();
                let t0 = Instant::now();
                let (hit, stats) = match *arm {
                    "f4" => {
                        groebner_decompose(&kc, &fb, &index_of, &st, &target, m, engine, budget)
                    }
                    "f6" => groebner_decompose_f6_ic(
                        &kc, &fb, &index_of, &st, &target, m, engine, budget,
                    ),
                    other => panic!("unknown engine {other}; use f4 or f6"),
                };
                let ns = t0.elapsed().as_nanos();
                let xors = f4_word_ops_thread() - x0;
                best_ns = best_ns.min(ns);
                last = Some((hit, stats, xors));
            }
            let (hit, stats, xors) = last.unwrap();
            if let Some(indices) = &hit {
                let sum = indices
                    .iter()
                    .fold(crypto_lib::binary_ecc::BinaryPoint::Infinity, |acc, &i| {
                        kc.add(&acc, &fb.points[i])
                    });
                assert_eq!(sum, target, "{arm}: witness for target {k} does not add up");
            }
            println!(
                "{{\"arm\":\"{arm}\",\"target\":{k},\"ms\":{:.3},\"hit\":{},\"exhausted\":{},\"reductions\":{},\"splits\":{},\"propagations\":{},\"infeasible\":{},\"geo_refute\":{},\"geo_witness\":{},\"partial_refute\":{},\"geo_adds\":{},\"word_xors\":{}}}",
                best_ns as f64 / 1e6,
                hit.is_some(),
                stats.exhausted,
                stats.reductions,
                stats.splits,
                stats.propagations,
                stats.infeasible_branches,
                stats.geometric_refutations,
                stats.geometric_witnesses,
                stats.geometric_partial_support_refutations,
                stats.geometric_group_additions,
                xors
            );
            arm_verdicts.push(hit.is_some());
            total_ns += best_ns;
            total.reductions += stats.reductions;
            total.splits += stats.splits;
            total.geometric_refutations += stats.geometric_refutations;
            total.geometric_witnesses += stats.geometric_witnesses;
            total.geometric_group_additions += stats.geometric_group_additions;
            total_xors += xors;
        }
        println!(
            "{{\"arm\":\"{arm}\",\"total_ms\":{:.3},\"hits\":{},\"reductions\":{},\"splits\":{},\"geo_refute\":{},\"geo_witness\":{},\"geo_adds\":{},\"word_xors\":{}}}",
            total_ns as f64 / 1e6,
            arm_verdicts.iter().filter(|h| **h).count(),
            total.reductions,
            total.splits,
            total.geometric_refutations,
            total.geometric_witnesses,
            total.geometric_group_additions,
            total_xors
        );
        verdicts.push(arm_verdicts);
    }
    for w in verdicts.windows(2) {
        assert_eq!(w[0], w[1], "the arms disagree on a verdict");
    }
}
