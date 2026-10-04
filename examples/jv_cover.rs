//! **Cover and decomposition on `E(F_{p⁶})`: `C_cov`, the cost of one Nagao test, and the method end to end.**
//!
//! Companion to `research/notes/index-calculus/RESEARCH_COVER_DECOMPOSITION_LEDGER.md`.
//!
//! ```text
//! cargo run --release --example jv_cover -- --exp ccov --sizes 53,101,251,503,1009 --seeds 2 --residuals 200 --constructed 40 --max-degree 24 --budget 300 --json experiments/30_jv_cover_ccov.json
//! cargo run --release --example jv_cover -- --exp dlp --sizes 53,101,251 --seeds 2 --rho-runs 16 --check-every 64 --json experiments/30_jv_cover_dlp.json
//! cargo run --release --example jv_cover -- --exp ccov --stop 64 --sizes 53,61,71 --seeds 2 --residuals 300 --constructed 60 --oracle-base 45 --json experiments/31_jv_cover_stop_ccov_oracle.json
//! cargo run --release --example jv_cover -- --exp dlp --sizes 503 --seeds 1 --rho-runs 0 --rho-ref <pooled S_rho of the smaller sizes> --json experiments/30_jv_cover_dlp_503.json
//! cargo run --release --example jv_cover -- --exp sieve --stop 64 --sizes 251,503 --seeds 2 --rho-runs 4 --json experiments/32_jv_cover_sieve_dlp.json
//! ```

use std::env;
use std::fs;

use crypto_lib::cryptanalysis::jv_cover::{
    run_cover_ccov, run_cover_dlp, CcovReport, CoverDlpReport,
};
use crypto_lib::cryptanalysis::jv_sieve::{run_cover_sieve_dlp, SieveDlpReport};
use serde::Serialize;

#[derive(Serialize)]
#[serde(untagged)]
#[allow(clippy::large_enum_variant)]
enum Row {
    Ccov(CcovReport),
    Dlp(CoverDlpReport),
    Sieve(SieveDlpReport),
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut exp = "ccov".to_string();
    let mut sizes: Vec<u64> = vec![53, 101, 251];
    let mut seeds = 1u64;
    let mut first_seed = 1u64;
    let mut residuals = 100usize;
    let mut constructed = 24usize;
    let mut max_degree = 24u32;
    let mut budget = 300f64;
    let mut rho_runs = 8usize;
    let mut check_every = 0u64;
    let mut oracle_base = 40usize;
    let mut rho_ref = 1.3f64;
    let mut stop = 0usize;
    let mut margin = 1.25f64;
    let mut m_override = 0usize;
    let mut json: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        let next = |i: &mut usize| -> String {
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
            "--first-seed" => first_seed = next(&mut i).parse().expect("--first-seed"),
            "--residuals" => residuals = next(&mut i).parse().expect("--residuals"),
            "--constructed" => constructed = next(&mut i).parse().expect("--constructed"),
            "--max-degree" => max_degree = next(&mut i).parse().expect("--max-degree"),
            "--budget" => budget = next(&mut i).parse().expect("--budget"),
            "--rho-runs" => rho_runs = next(&mut i).parse().expect("--rho-runs"),
            "--check-every" => check_every = next(&mut i).parse().expect("--check-every"),
            "--rho-ref" => rho_ref = next(&mut i).parse().expect("--rho-ref"),
            "--stop" => stop = next(&mut i).parse().expect("--stop"),
            "--margin" => margin = next(&mut i).parse().expect("--margin"),
            "--m" => m_override = next(&mut i).parse().expect("--m"),
            "--oracle-base" => oracle_base = next(&mut i).parse().expect("--oracle-base"),
            "--json" => json = Some(next(&mut i)),
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    let mut rows: Vec<Row> = Vec::new();
    for &p in &sizes {
        for seed in first_seed..first_seed + seeds {
            match exp.as_str() {
                "ccov" => {
                    let r = run_cover_ccov(
                        p,
                        seed,
                        residuals,
                        constructed,
                        max_degree,
                        budget,
                        oracle_base,
                        (stop > 0).then_some(stop),
                    );
                    println!(
                        "p={:>5} seed={} l=2^{:.1} |F|={:>4} c_add E {:.0} J {:.0} stop {:?} ({} stopped) | random {} decomposable {} (1/720 = {:.5}) | planted {}/{} | oracle {} mismatches {} unverified {} incomplete {} timed out {} | C_cov {:.3e} [{:.2e}, {:.2e}] = weil {:.0} + f4 {:.3e} + lin {:.3e} + post {:.0} | delta {:.0} F4 degree {:.1} rows {:.0} cols {:.0} {:.1} ms + lin {:.1} ms = {:.1} ms | decomposable C_cov {:.3e} | {:.0} s",
                        r.p, r.seed, r.bits, r.base, r.c_add_e, r.c_add_j, r.stop_staircase, r.stopped, r.random_residuals, r.random_decomposable,
                        r.expected_rate, r.planted_found, r.constructed_residuals, r.oracle_checked, r.mismatches,
                        r.unverified, r.incomplete, r.timed_out, r.c_cov.mean, r.c_cov.min as f64, r.c_cov.max as f64,
                        r.weil_muls.mean, r.f4_muls.mean, r.lin_muls.mean, r.post_muls.mean, r.delta.mean,
                        r.degree_reached.mean, r.max_rows.mean, r.max_cols.mean, r.f4_ms.mean, r.lin_ms.mean, r.total_ms.mean,
                        r.c_cov_decomposable.mean, r.wall_ms / 1e3
                    );
                    rows.push(Row::Ccov(r));
                }
                "dlp" => {
                    let r = run_cover_dlp(
                        p,
                        seed,
                        rho_runs,
                        check_every,
                        rho_ref,
                        (stop > 0).then_some(stop),
                    );
                    println!(
                        "p={:>5} seed={} l=2^{:.1} |F|={:>4} c_add E {:.0} J {:.0} stop {:?} ({} stopped) | residuals {} decomposable {} ({:.5} vs 1/720 {:.5}) relations {} (floor {:.0}, ratio {:.2}) | incomplete {} timed out {} cross-checked {} mismatches {} | C_cov {:.3e} | unknowns {} (filtered {}) phi {:.2} row {:.1} LA attempts {} ops {} ({:.1}·u²) | solved {} correct {} | S {:.3e} | rho S {:.3} ± {:.3} ({} runs) | S/rho {:.3} = relation {:.3} + la {:.4} (last {:.4}) | predicted {:.3} | {:.0} s",
                        r.p, r.seed, r.bits, r.base, r.c_add_e, r.c_add_j, r.stop_staircase, r.stopped, r.residuals, r.decompositions,
                        r.decomposition_rate, r.expected_rate, r.relations, r.floor_residuals, r.residual_ratio,
                        r.incomplete, r.timed_out, r.cross_checked, r.mismatches, r.c_cov, r.unknowns, r.filtered_out,
                        r.phi, r.row_weight, r.la_attempts, r.la_ops, r.wiedemann_constant, r.solved, r.correct, r.s,
                        r.rho_s_mean, r.rho_s_sd, r.rho.len(), r.s_over_rho, r.relation_over_rho, r.r, r.r_last,
                        r.predicted_s_over_rho, r.wall_ms / 1e3
                    );
                    rows.push(Row::Dlp(r));
                }
                "sieve" => {
                    let r = run_cover_sieve_dlp(
                        p,
                        seed,
                        rho_runs,
                        rho_ref,
                        (stop > 0).then_some(stop),
                        margin,
                        (m_override > 0).then_some(m_override),
                    );
                    println!(
                        "p={:>5} seed={} l=2^{:.1} |F|={:>4} c_add E {:.0} J {:.0} m {} (rule {}, {:.0} available) | B's {} lines {} steps {} (base {}) roots {} hits {} false {} | relations {} ({:.3e}/line vs p/m! {:.3e}: ratio {:.2}) verified {} failed {} duplicates {} per-m {:?} | C_rel {:.3e} = enum {:.3e} + sieve {:.3e} + extract {:.3e} + verify {:.3e} (per relation) | adds {:.3e} lookups {:.3e} | descent: residuals {} successes {} ({:.0}/success) C_cov {:.3e} muls {:.3e} stopped {} incomplete {} timed out {} | unknowns {} (filtered {}) row {:.1} LA attempts {} ops {} | solved {} correct {} exhausted {} | S {:.3e} S+ {:.3e} | rho S {:.3} ± {:.3} ({} runs) | S/rho {:.4} (S+/rho {:.4}) = relation {:.4} + descent {:.4} + la {:.4} | {:.0} s",
                        r.p, r.seed, r.bits, r.base, r.c_add_e, r.c_add_j, r.m_final, r.m_rule, r.relations_available,
                        r.bs, r.lines, r.sieve_steps, r.base_steps, r.roots, r.hits, r.false_hits,
                        r.relations, r.rels_per_line, r.expected_rels_per_line, r.rate_ratio, r.rels_verified, r.rels_failed_verify, r.duplicates, r.per_m,
                        r.c_rel, r.enum_muls as f64 / r.relations.max(1) as f64, r.sieve_muls as f64 / r.relations.max(1) as f64,
                        r.extract_muls as f64 / r.relations.max(1) as f64, r.verify_muls as f64 / r.relations.max(1) as f64,
                        r.sieve_adds as f64, r.sieve_lookups as f64,
                        r.descent_residuals, r.descent_successes, r.descent_tests_per_success, r.descent_c_cov, r.descent_muls as f64,
                        r.descent_stopped, r.descent_incomplete, r.descent_timed_out,
                        r.unknowns, r.filtered_out, r.row_weight, r.la_attempts, r.la_ops, r.solved, r.correct, r.exhausted,
                        r.s, r.s_plus, r.rho_s_mean, r.rho_s_sd, r.rho.len(), r.s_over_rho, r.s_plus_over_rho,
                        r.relation_over_rho, r.descent_over_rho, r.la_over_rho, r.wall_ms / 1e3
                    );
                    rows.push(Row::Sieve(r));
                }
                other => panic!("unknown experiment {other}"),
            }
            if let Some(path) = &json {
                fs::write(path, serde_json::to_string_pretty(&rows).unwrap()).expect("write json");
            }
        }
    }
}
