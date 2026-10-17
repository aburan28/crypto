//! Same-binary timing probe for the PDP3 Gröbner decomposers: inherited F4
//! against F6-IC on a standard-subspace base of one Koblitz curve.
//!
//! ```bash
//! cargo run --release --example f6_ic_probe -- [--a 0] [--n 17] [--dim 6] [--m 3] \
//!     [--targets 8] [--budget 20000] [--engines f4,f6,f6p] [--reps 2]
//! ```
//!
//! Arms: `f4` inherited F4, `f6` F6-IC, `f6p` F6-IC with the shared pair-sum
//! index.  Targets are `[k]G` for `k = 1 …`.  Each arm decides every target with
//! the same node budget; the probe prints one JSON line per arm and target
//! (wall, verdict, reductions, splits, geometric counters, word XORs) and a
//! per-arm total, and checks that the arms agree on every verdict.  It is a
//! development probe: wall times on a shared host are indicative only.
//!
//! With `KIC_F4_LANES=k > 1` (lit2026 W5) the targets are decided in
//! batches of `k`: each repetition of a batch first builds the batch's
//! roots together (`prebuild_decomposition_roots`, timed with the batch)
//! and then decides its targets, which take those roots.  A target's `ms`
//! is then its batch's best time divided by the batch size.  The per-arm
//! total also reports the process CPU time the arm used (`cpu_ms`, user and
//! system, all repetitions), the Macaulay rows generated, the stage units
//! charged to root eliminations, and the nodes and specialisations the
//! gate-first order (default; `KIC_F6_GATE_FIRST=0` turns it off) skipped.

use crypto_lib::cryptanalysis::inherited_f4::{
    f4_root_elimination_thread, f4_rows_generated_thread,
};
use crypto_lib::cryptanalysis::koblitz_groebner::{
    clear_lane_roots, f4_lanes, f4_profile, f4_profile_reset, f4_word_ops_thread,
    gate_first_skips_thread, lane_roots_pending, FieldStructure, SolveStats, SolverEngine,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, groebner_decompose, groebner_decompose_f6_ic,
    groebner_decompose_f6_ic_pair_index, prebuild_decomposition_roots, KoblitzCurve,
};
use num_bigint::BigUint;
use std::time::Instant;

fn flag(name: &str) -> Option<String> {
    let args: Vec<String> = std::env::args().collect();
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1).cloned())
}

/// User plus system CPU time of this process, in nanoseconds.
fn cpu_ns() -> u128 {
    // SAFETY: `getrusage` writes into the zeroed struct it is handed.
    let usage = unsafe {
        let mut usage: libc::rusage = std::mem::zeroed();
        libc::getrusage(libc::RUSAGE_SELF, &mut usage);
        usage
    };
    let tv = |t: libc::timeval| t.tv_sec as u128 * 1_000_000_000 + t.tv_usec as u128 * 1_000;
    tv(usage.ru_utime) + tv(usage.ru_stime)
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
    let lanes = f4_lanes();

    let kc = KoblitzCurve::new(a, n).expect("curve has a usable subgroup");
    let fb = build_standard_subspace_factor_base(&kc, dim).expect("factor base");
    let index_of = fb.index_map();
    let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
    let engine = SolverEngine::InheritedF4 { max_degree: 3 };
    eprintln!(
        "K({a}, {n}), standard subspace dim {dim}: {} base points, m = {m}, {targets} targets, budget {budget}, lanes {lanes}",
        fb.points.len()
    );
    let points: Vec<_> = (1..=targets)
        .map(|k| kc.mul(kc.generator(), &BigUint::from(k)))
        .collect();

    let mut verdicts: Vec<Vec<bool>> = Vec::new();
    for arm in &engines {
        let decide = |target| match *arm {
            "f4" => groebner_decompose(&kc, &fb, &index_of, &st, target, m, engine, budget),
            "f6" => groebner_decompose_f6_ic(&kc, &fb, &index_of, &st, target, m, engine, budget),
            "f6p" => groebner_decompose_f6_ic_pair_index(
                &kc, &fb, &index_of, &st, target, m, engine, budget,
            ),
            other => panic!("unknown engine {other}; use f4, f6 or f6p"),
        };
        let mut arm_verdicts = Vec::new();
        let mut total_ns = 0u128;
        let mut total = SolveStats::default();
        let mut total_xors = 0u64;
        let mut total_oracle_ns = 0u64;
        let mut lane_built = 0usize;
        f4_profile_reset();
        let cpu0 = cpu_ns();
        let rows0 = f4_rows_generated_thread();
        let root0 = f4_root_elimination_thread();
        let skips0 = gate_first_skips_thread();
        for (b, batch) in points.chunks(lanes).enumerate() {
            let mut best_ns = u128::MAX;
            let mut last: Vec<(Option<Vec<usize>>, SolveStats, u64)> = Vec::new();
            for _ in 0..reps {
                last.clear();
                let t0 = Instant::now();
                if lanes > 1 {
                    clear_lane_roots();
                    lane_built += prebuild_decomposition_roots(&kc, &fb, &st, batch, m, engine).1;
                }
                for target in batch {
                    let x0 = f4_word_ops_thread();
                    let (hit, stats) = decide(target);
                    last.push((hit, stats, f4_word_ops_thread() - x0));
                }
                best_ns = best_ns.min(t0.elapsed().as_nanos());
                if lanes > 1 {
                    assert_eq!(lane_roots_pending(), 0, "a lane-built root was not used");
                }
            }
            for (i, (hit, stats, xors)) in last.into_iter().enumerate() {
                let k = (b * lanes + i + 1) as u64;
                let target = &points[k as usize - 1];
                if let Some(indices) = &hit {
                    let sum = indices
                        .iter()
                        .fold(crypto_lib::binary_ecc::BinaryPoint::Infinity, |acc, &i| {
                            kc.add(&acc, &fb.points[i])
                        });
                    assert_eq!(
                        &sum, target,
                        "{arm}: witness for target {k} does not add up"
                    );
                }
                println!(
                    "{{\"arm\":\"{arm}\",\"target\":{k},\"ms\":{:.3},\"hit\":{},\"indices\":{:?},\"exhausted\":{},\"reductions\":{},\"splits\":{},\"propagations\":{},\"infeasible\":{},\"geo_refute\":{},\"geo_witness\":{},\"partial_refute\":{},\"geo_adds\":{},\"oracle_ms\":{:.3},\"pair_builds\":{},\"pair_lookups\":{},\"word_xors\":{}}}",
                    best_ns as f64 / 1e6 / batch.len() as f64,
                    hit.is_some(),
                    hit.clone().unwrap_or_default(),
                    stats.exhausted,
                    stats.reductions,
                    stats.splits,
                    stats.propagations,
                    stats.infeasible_branches,
                    stats.geometric_refutations,
                    stats.geometric_witnesses,
                    stats.geometric_partial_support_refutations,
                    stats.geometric_group_additions,
                    stats.geometric_oracle_ns as f64 / 1e6,
                    stats.geometric_pair_index_builds,
                    stats.geometric_pair_index_lookups,
                    xors
                );
                arm_verdicts.push(hit.is_some());
                total.reductions += stats.reductions;
                total.splits += stats.splits;
                total.geometric_refutations += stats.geometric_refutations;
                total.geometric_witnesses += stats.geometric_witnesses;
                total.geometric_group_additions += stats.geometric_group_additions;
                total_oracle_ns += stats.geometric_oracle_ns;
                total_xors += xors;
            }
            total_ns += best_ns;
        }
        let cpu = cpu_ns() - cpu0;
        let rows = f4_rows_generated_thread();
        let skips = gate_first_skips_thread();
        println!(
            "{{\"arm\":\"{arm}\",\"total_ms\":{:.3},\"cpu_ms\":{:.3},\"reps\":{reps},\"lanes\":{lanes},\"lane_built\":{lane_built},\"hits\":{},\"reductions\":{},\"splits\":{},\"geo_refute\":{},\"geo_witness\":{},\"geo_adds\":{},\"oracle_ms\":{:.3},\"word_xors\":{},\"root_elim_xors\":{},\"rows_generated\":{},\"row_words_generated\":{},\"gate_first_nodes\":{},\"gate_first_specialisations\":{}}}",
            total_ns as f64 / 1e6,
            cpu as f64 / 1e6,
            arm_verdicts.iter().filter(|h| **h).count(),
            total.reductions,
            total.splits,
            total.geometric_refutations,
            total.geometric_witnesses,
            total.geometric_group_additions,
            total_oracle_ns as f64 / 1e6,
            total_xors,
            f4_root_elimination_thread() - root0,
            rows.0 - rows0.0,
            rows.1 - rows0.1,
            skips.0 - skips0.0,
            skips.1 - skips0.1,
        );
        // The solver's own timers over the arm (all repetitions): the
        // Macaulay build (specialisation included), the reduction and the
        // read-back, and the word reads and writes of specialisation.
        let prof = f4_profile();
        println!(
            "{{\"arm\":\"{arm}\",\"profile\":{{\"calls\":{},\"build_ms\":{:.3},\"reduce_ms\":{:.3},\"readback_ms\":{:.3},\"rows\":{},\"cols\":{},\"word_ops\":{},\"specialise_word_ops\":{}}}}}",
            prof.calls,
            prof.build_ns as f64 / 1e6,
            prof.reduce_ns as f64 / 1e6,
            prof.readback_ns as f64 / 1e6,
            prof.rows,
            prof.cols,
            prof.word_ops,
            prof.specialise_word_ops
        );
        verdicts.push(arm_verdicts);
    }
    for w in verdicts.windows(2) {
        assert_eq!(w[0], w[1], "the arms disagree on a verdict");
    }
}
