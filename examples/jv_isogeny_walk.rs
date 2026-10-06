//! **The isogeny walk to a weak curve over `F_{p⁶}`, priced at toy sizes.**
//!
//! Companion to `research/notes/index-calculus/RESEARCH_COVER_DECOMPOSITION_LEDGER.md` §13
//! (the first walk) and §17 (`--v2`: the rebuilt walk).
//!
//! ```text
//! cargo run --release --example jv_isogeny_walk -- --sizes 7,11,13,17,23,31,53 --trials 40 --cap 200000 --samples 20000 --json experiments/34_jv_isogeny_walk.json
//! cargo run --release --example jv_isogeny_walk -- --v2 --sizes 7,11,13,17,23,31,53,101,251,503 --trials 40 --cap-mult 3 --jumps 3,5,7 --samples 20000 --json experiments/39_jv_isogeny_walk_v2.json
//! ```
//!
//! `--v2` prints, besides the log line, the Markdown row of ledger §17.5
//! (rho at `p` is `1.3 · (p³/2) · 331`, as section G prices it).

use std::env;
use std::fs;

use crypto_lib::cryptanalysis::jv_isogeny_walk::{run_walk, run_walk2, Walk2Report, WalkReport};

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut sizes: Vec<u64> = vec![7, 11, 13];
    let mut trials = 20usize;
    let mut cap = 100_000u64;
    let mut samples = 10_000u64;
    let mut seed = 1u64;
    let mut two_only = false;
    let mut v2 = false;
    let mut closure = false;
    let mut cap_mult = 3u64;
    let mut jumps: Vec<u64> = vec![3, 5, 7];
    let mut json: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        let next = |i: &mut usize| -> String {
            *i += 1;
            args[*i].clone()
        };
        match args[i].as_str() {
            "--sizes" => {
                sizes = next(&mut i)
                    .split(',')
                    .map(|v| v.parse().expect("--sizes"))
                    .collect()
            }
            "--trials" => trials = next(&mut i).parse().expect("--trials"),
            "--cap" => cap = next(&mut i).parse().expect("--cap"),
            "--samples" => samples = next(&mut i).parse().expect("--samples"),
            "--seed" => seed = next(&mut i).parse().expect("--seed"),
            "--two-only" => two_only = true,
            "--v2" => v2 = true,
            "--closure" => closure = true,
            "--cap-mult" => cap_mult = next(&mut i).parse().expect("--cap-mult"),
            "--jumps" => {
                jumps = next(&mut i)
                    .split(',')
                    .filter(|v| !v.is_empty())
                    .map(|v| v.parse().expect("--jumps"))
                    .collect()
            }
            "--json" => json = Some(next(&mut i)),
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    if v2 {
        let mut rows: Vec<Walk2Report> = Vec::new();
        println!("| p | q | weak × q (sampled) | jumps | walks found / capped / exhausted (start weak) | success | curves met median / mean (q/3) | first component mean, weak share | components mean | jumps 3/5/7, wasted | muls per curve | muls per jump | muls per found walk | rho at p | walk / rho (found) | walk / rho (all) |");
        println!("|---:|--:|--:|:--|:--|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|");
        for &p in &sizes {
            let r = run_walk2(p, seed, trials, cap_mult * p * p, &jumps, samples, closure);
            eprintln!(
                "p={:>5} q={:>8} | weak×q {:.2} ({}) | found {} capped {} exhausted {} start-weak {} of {} | success {:.2} | curves median {:.0} mean {:.0} (q/3 {:.0}) | first comp {:.1} weak {:.2} | comps {:.1} | jumps {:?} wasted {} | c_curve {:.0} c_jump {:.3e} | muls/found walk {:.3e} all {:.3e} | walk/rho {:.4} all {:.4} | {:.0} s",
                r.p, r.q, r.weak_fraction_times_q, r.sampled, r.found, r.capped, r.exhausted, r.start_weak, r.trials, r.success_fraction, r.median_curves, r.mean_curves, r.q_over_3,
                r.mean_first_component, r.frac_first_component_weak, r.mean_components, r.total_jumps, r.total_wasted_jumps,
                r.c_curve, r.c_jump, r.mean_muls, r.mean_muls_all, r.walk_over_rho, r.walk_over_rho_all, r.wall_ms / 1e3
            );
            println!(
                "| {} | {} | {:.2} ({}) | {}{} | {} / {} / {} ({}) | {:.2} | {:.0} / {:.0} ({:.0}) | {:.1}, {:.2} | {:.1} | {}/{}/{}, {} | {:.0} | {:.2e} | {:.3e} | {:.3e} | {:.4} | {:.4} |",
                r.p, r.q, r.weak_fraction_times_q, r.sampled,
                r.jump_degrees.iter().map(|d| d.to_string()).collect::<Vec<_>>().join(","),
                if r.closure_mode { " (closure)" } else { "" },
                r.found, r.capped, r.exhausted, r.start_weak, r.success_fraction, r.median_curves, r.mean_curves, r.q_over_3,
                r.mean_first_component, r.frac_first_component_weak, r.mean_components,
                r.total_jumps[0], r.total_jumps[1], r.total_jumps[2], r.total_wasted_jumps,
                r.c_curve, r.c_jump, r.mean_muls, r.rho_at_p, r.walk_over_rho, r.walk_over_rho_all
            );
            rows.push(r);
            if let Some(path) = &json {
                fs::write(path, serde_json::to_string_pretty(&rows).unwrap()).expect("write json");
            }
        }
        return;
    }
    let mut rows: Vec<WalkReport> = Vec::new();
    for &p in &sizes {
        let r = run_walk(p, seed, trials, cap, !two_only, samples);
        println!(
            "p={:>4} q={:>6} | weak fraction {:.5} (× q = {:.2}; {} sampled) | walks {} found {} exhausted {} capped {} | steps mean {:.1} median {:.0} (q/3 = {:.0}, cited q = {}) | distinct j mean {:.1} | 2-steps {} 3-steps {} restarts {} | F_p muls per step {:.0}, per walk {:.3e} | {:.0} s",
            r.p, r.q, r.weak_fraction, r.weak_fraction_times_q, r.sampled, r.trials, r.found, r.exhausted, r.capped, r.mean_steps, r.median_steps, r.q_over_3, r.cited_q,
            r.mean_distinct, r.rows.iter().map(|t| t.two_steps).sum::<u64>(), r.rows.iter().map(|t| t.three_steps).sum::<u64>(),
            r.rows.iter().map(|t| t.stuck_restarts).sum::<u64>(), r.mean_muls_per_step, r.mean_muls, r.wall_ms / 1e3
        );
        rows.push(r);
        if let Some(path) = &json {
            fs::write(path, serde_json::to_string_pretty(&rows).unwrap()).expect("write json");
        }
    }
}
