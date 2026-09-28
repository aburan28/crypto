//! **GLV / `C₃` experiments on the `E(F_{p³})` index-calculus harness —
//! measurement.**
//!
//! Companion to `crypto_lib::cryptanalysis::glv_gaudry` and
//! `RESEARCH_GLV_INDEX_CALCULUS.md`.  Every run generates a prime-order
//! `j = 0` curve over `F_{p³}` (`p ≡ 1 mod 3`) with its order-3
//! automorphism and runs one of the four experiments:
//!
//! ```text
//! cargo run --release --example glv_gaudry_bench -- --exp quotient  --sizes 271,541,1051 --seeds 2 --json experiments/22_glv_quotient.json
//! cargo run --release --example glv_gaudry_bench -- --exp canonical --sizes 271,541,1051 --seeds 2 --json experiments/22_glv_canonical.json
//! cargo run --release --example glv_gaudry_bench -- --exp graded    --sizes 271,541,1051 --seeds 1 --residuals 6 --json experiments/22_glv_graded.json
//! cargo run --release --example glv_gaudry_bench -- --exp invariant --sizes 271,541,1051 --seeds 1 --residuals 4 --json experiments/22_glv_invariant.json
//! ```

use std::env;
use std::fs;

use crypto_lib::cryptanalysis::glv_gaudry::{
    generate_j0_instance3, run_canonical_experiment, run_graded_experiment,
    run_invariant_experiment, run_quotient_experiment, CanonicalReport, GradedReport,
    InvariantReport, QuotientReport,
};
use serde::Serialize;

#[derive(Serialize)]
#[serde(untagged)]
enum Row {
    Quotient { seed: u64, report: QuotientReport },
    Canonical { seed: u64, report: CanonicalReport },
    Graded { seed: u64, report: GradedReport },
    Invariant { seed: u64, report: InvariantReport },
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut exp = "quotient".to_string();
    let mut sizes: Vec<u64> = vec![271, 541, 1051];
    let mut seeds = 2u64;
    let mut json: Option<String> = None;
    let mut max_residuals = 1_000_000u64;
    let mut residuals = 6u64;
    let mut random_residuals = 0u64;
    let mut pair_residuals = 3000u64;
    let mut verify_cap = 200u64;
    let mut d_max_ordinary = 13u8;
    let mut d_max_orbit = 19u8;
    let mut ff_d_max = 6u8;
    let mut f4_max_degree = 14u32;
    let mut f4_budget = 120.0f64;
    let mut cell_cap = 1u64 << 27;
    let mut i = 0;
    while i < args.len() {
        let mut next = |i: &mut usize| -> String {
            *i += 1;
            args[*i].clone()
        };
        match args[i].as_str() {
            "--exp" => exp = next(&mut i),
            "--sizes" => {
                sizes = next(&mut i)
                    .split(',')
                    .map(|v| v.parse().expect("--sizes"))
                    .collect()
            }
            "--seeds" => seeds = next(&mut i).parse().expect("--seeds"),
            "--json" => json = Some(next(&mut i)),
            "--max-residuals" => max_residuals = next(&mut i).parse().expect("--max-residuals"),
            "--residuals" => residuals = next(&mut i).parse().expect("--residuals"),
            "--random-residuals" => {
                random_residuals = next(&mut i).parse().expect("--random-residuals")
            }
            "--pair-residuals" => pair_residuals = next(&mut i).parse().expect("--pair-residuals"),
            "--verify-cap" => verify_cap = next(&mut i).parse().expect("--verify-cap"),
            "--d-max-ordinary" => d_max_ordinary = next(&mut i).parse().expect("--d-max-ordinary"),
            "--d-max-orbit" => d_max_orbit = next(&mut i).parse().expect("--d-max-orbit"),
            "--ff-d-max" => ff_d_max = next(&mut i).parse().expect("--ff-d-max"),
            "--f4-max-degree" => f4_max_degree = next(&mut i).parse().expect("--f4-max-degree"),
            "--f4-budget" => f4_budget = next(&mut i).parse().expect("--f4-budget"),
            "--cell-cap" => cell_cap = next(&mut i).parse().expect("--cell-cap"),
            other => {
                eprintln!("unknown argument {other}");
                std::process::exit(2);
            }
        }
        i += 1;
    }
    let mut rows: Vec<Row> = Vec::new();
    for &p in &sizes {
        for seed in 1..=seeds {
            let glv = generate_j0_instance3(p, seed);
            let n = glv.inst.curve.n;
            match exp.as_str() {
                "quotient" => {
                    let r = run_quotient_experiment(&glv, seed, max_residuals);
                    for st in [&r.control, &r.quotient] {
                        println!(
                            "p={p:>5} n=2^{:5.1} seed={seed} {:<18} cols={:>5} rels={:>5} (zero {:>4} dup {:>4} merged {:>4}) residuals={:>6} nnz={:>7} avg_w={:>5.2} bytes={:>9} core={:>5}x nnz {:>7} LA attempts={} ops={:>10} wiedemann={:>8.1}ms | group_ops={:>9} oracle_muls={:>13} total_ops={:>12.0} S={:>9.1} floor_res={:>8.1} ratio={:>5.2} ok={:?}",
                            r.bits, st.label, st.columns, st.relations, st.rows_zero, st.rows_duplicate, st.merged_entries, st.residuals_at_solve, st.nnz, st.avg_weight, st.bytes,
                            st.core_rows, st.core_nnz, st.la_attempts, st.la_ops, st.wiedemann_ms,
                            st.group_ops_at_solve, st.oracle_fp_muls_at_solve, st.total_ops, st.s, st.floor_residuals, st.residual_ratio, st.correct
                        );
                    }
                    println!(
                        "         base={} orbits={} rate={:.3} solves={} unsolved={} Fp/add={:.1} peak_rss={} MB wall={:.0} ms | rho S={:.2} ok={:?}",
                        r.base, r.orbits, r.decomposition_rate, r.solve_stats.solves, r.solve_stats.unsolved,
                        r.fp_muls_per_add, r.peak_rss_bytes >> 20, r.wall_ms, r.rho.s, r.rho.correct
                    );
                    rows.push(Row::Quotient { seed, report: r });
                }
                "canonical" => {
                    let t = if random_residuals == 0 {
                        // About what the quotient pipeline needs: |F|/(3ρ) with ρ ≈ 1/6.
                        (p as u64 / 2 / 3 * 6).max(50)
                    } else {
                        random_residuals
                    };
                    let r = run_canonical_experiment(&glv, seed, t, pair_residuals, verify_cap);
                    for st in [&r.random, &r.pairs] {
                        println!(
                            "p={p:>5} n=2^{:5.1} seed={seed} {:<22} residuals={:>6} orbits={:>6} dups={:>5} expected_uniform={:>8.4} saved={:>5} verified={:>4}/{} mismatches={} rows={:>5} zero={:>4} duplicate={:>4} merged={:>3} | {:>8.0} ms",
                            r.bits, st.label, st.residuals, st.distinct_orbits, st.duplicates, st.expected_duplicates_uniform,
                            st.solver_calls_saved, st.verified_duplicates, st.duplicates, st.verification_mismatches,
                            st.rows_produced, st.rows_zero, st.rows_duplicate, st.merged_entries, st.wall_ms
                        );
                    }
                    rows.push(Row::Canonical { seed, report: r });
                }
                "graded" => {
                    let r = run_graded_experiment(
                        &glv,
                        seed,
                        residuals,
                        d_max_ordinary,
                        d_max_orbit,
                        cell_cap,
                        f4_max_degree,
                        f4_budget,
                    );
                    println!(
                        "p={p:>5} n=2^{:5.1} seed={seed} S4 terms={} weight histogram={:?} homogeneous={}",
                        r.bits, r.s4_terms, r.s4_weight_histogram, r.s4_weight_homogeneous
                    );
                    for row in &r.rows {
                        let show = |s: &crypto_lib::cryptanalysis::glv_gaudry::DegreeSweep| {
                            format!(
                                "D={:?} dim={:?} peak={}x{} cells={} muls@D={} muls_total={} ms={:.0}{}",
                                s.solve_degree, s.quotient_dim, s.peak_rows, s.peak_cols, s.peak_cells,
                                s.solve_muls, s.total_muls, s.total_ms,
                                s.stopped.as_ref().map(|x| format!(" [{x}]")).unwrap_or_default()
                            )
                        };
                        println!(
                            "  residual {} harness: triples={:?} dim={} muls={} | ordinary coupled={:.2} {} | orbit vanilla {} | orbit block {} | ratios block/ordinary={:?} block/vanilla muls={:?} cells={:?}",
                            row.residual_index, row.harness_triples, row.harness_quotient_dim, row.harness_fp_muls,
                            row.ordinary_coupled_fraction, show(&row.ordinary), show(&row.orbit_vanilla), show(&row.orbit_block),
                            row.muls_ratio_block_vs_ordinary, row.muls_ratio_block_vs_orbit_vanilla, row.cells_ratio_block_vs_orbit_vanilla
                        );
                        for s in [&row.ordinary_f4, &row.orbit_f4_vanilla, &row.orbit_f4_block] {
                            println!(
                                "    F4 {}: D_solve={} D_reach={} max={}x{} runs={} ops={} sols={:?} blocked_steps={} ms={:.0}{}",
                                s.system, s.solving_degree, s.degree_reached, s.max_rows, s.max_cols, s.f4_runs, s.field_ops,
                                s.solutions, s.blocked_steps, s.ms, if s.timed_out { " TIMED OUT" } else { "" }
                            );
                        }
                        println!(
                            "    F4 ratios: ops block/ordinary={:?} block/vanilla={:?} ms block/vanilla={:?}; orbit = 3 × ordinary {:?}",
                            row.f4_ops_ratio_block_vs_ordinary, row.f4_ops_ratio_block_vs_vanilla, row.f4_ms_ratio_block_vs_vanilla, row.orbit_f4_triples_ordinary
                        );
                    }
                    println!(
                        "  (residuals without an F_p-rational decomposition skipped: {})",
                        r.residuals_skipped
                    );
                    rows.push(Row::Graded { seed, report: r });
                }
                "invariant" => {
                    let r = run_invariant_experiment(
                        &glv,
                        seed,
                        residuals,
                        f4_max_degree,
                        f4_budget,
                        ff_d_max,
                        cell_cap,
                    );
                    println!(
                        "p={p:>5} n=2^{:5.1} seed={seed} generators (1,2,0,1): {:?} hilbert {:?}; diagonal (1,1,1,1): {} generators, hilbert {:?}",
                        r.bits, r.generators_symmetrised, r.hilbert_symmetrised, r.generators_diagonal.len(), r.hilbert_diagonal
                    );
                    for row in &r.rows {
                        let f4 = |s: &crypto_lib::cryptanalysis::glv_gaudry::F4Summary| {
                            format!(
                                "{}: vars={} eqs={} D_solve={} D_reach={} max={}x{} runs={} ops={} sols={:?} basis={:?} ms={:.0}{}",
                                s.system, s.n_vars, s.equations, s.solving_degree, s.degree_reached, s.max_rows, s.max_cols,
                                s.f4_runs, s.field_ops, s.solutions, s.basis_size, s.ms, if s.timed_out { " TIMED OUT" } else { "" }
                            )
                        };
                        println!(
                            "  residual {} harness triples={:?} e_sols={} muls={} | orbit-invariant-only={} identical={} | ordinary Macaulay D={:?} dim={:?} muls={} | invariant Macaulay D={:?} dim={:?} muls={} | function-first Macaulay D={:?} peak={}x{} muls={} {:?}",
                            row.residual_index, row.harness_triples.as_ref().map(|t| t.len()), row.harness_e_solutions, row.harness_fp_muls,
                            row.invariant_uses_orbit_invariant_only, row.invariant_identical_to_ordinary,
                            row.ordinary_macaulay.solve_degree, row.ordinary_macaulay.quotient_dim, row.ordinary_macaulay.solve_muls,
                            row.invariant_macaulay.solve_degree, row.invariant_macaulay.quotient_dim, row.invariant_macaulay.solve_muls,
                            row.function_first_macaulay.solve_degree, row.function_first_macaulay.peak_rows, row.function_first_macaulay.peak_cols,
                            row.function_first_macaulay.total_muls, row.function_first_macaulay.stopped
                        );
                        println!("    F4 {}", f4(&row.ordinary_f4));
                        println!("    F4 {}", f4(&row.invariant_f4));
                        println!("    F4 {}", f4(&row.orbit_f4));
                        println!("    F4 {}", f4(&row.function_first_f4));
                        println!(
                            "    checks: ordinary F4 = harness {:?}; orbit = 3 × ordinary {:?}; function-first witnessed by harness triples {:?} ({} witnesses) [skipped residuals so far: {}]",
                            row.ordinary_f4_matches_harness, row.orbit_f4_triples_ordinary, row.function_first_witnessed, row.function_first_witnesses, r.residuals_skipped
                        );
                    }
                    rows.push(Row::Invariant { seed, report: r });
                }
                other => {
                    eprintln!("unknown experiment {other}");
                    std::process::exit(2);
                }
            }
            let _ = n;
        }
    }
    if let Some(path) = json {
        fs::write(&path, serde_json::to_string_pretty(&rows).unwrap()).expect("write json");
        println!("wrote {path}");
    }
}
