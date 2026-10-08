//! **The isogeny walk to a weak curve over `F_{p⁶}`, priced at toy sizes.**
//!
//! Companion to `research/notes/index-calculus/RESEARCH_COVER_DECOMPOSITION_LEDGER.md` §12.
//!
//! ```text
//! cargo run --release --example jv_isogeny_walk -- --sizes 7,11,13,17,23,31,53 --trials 40 --cap 200000 --samples 20000 --json experiments/34_jv_isogeny_walk.json
//! ```

use std::env;
use std::fs;

use crypto_lib::cryptanalysis::jv_isogeny_walk::{run_walk, WalkReport};

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut sizes: Vec<u64> = vec![7, 11, 13];
    let mut trials = 20usize;
    let mut cap = 100_000u64;
    let mut samples = 10_000u64;
    let mut seed = 1u64;
    let mut two_only = false;
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
            "--json" => json = Some(next(&mut i)),
            other => panic!("unknown argument {other}"),
        }
        i += 1;
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
