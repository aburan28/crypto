//! `pdp_split_dreg` — solving degree and cost of unsplit, split and hybrid
//! point-decomposition systems on Koblitz curves
//! (`research/pdp_split_dreg_20261011/PROTOCOL.md`).
//!
//! ```sh
//! cargo run --release --bin pdp_split_dreg -- \
//!   --cells 3:1:11:4,4:1:11:3 --guesses 0,2,4 --targets 4 \
//!   --seed 20261011 --budget-secs 120 --out runs.jsonl
//! ```
//!
//! A cell is `m:a:n:ℓ`: `m` summands on the Koblitz curve `K_a / GF(2^n)`
//! over the standard subspace of dimension `ℓ`.  Every target is solved
//! under every plan of its `m` (`[3]` and `[2,1]` for `m = 3`; `[2,2]` and
//! `[2,1,1]` for `m = 4`; `--plans` overrides) and every guess count, plans
//! interleaved per target.  One JSON line per cell header and per solve.

use std::fs::OpenOptions;
use std::io::Write;
use std::time::Duration;

use crypto_lib::cryptanalysis::ic_boundary::koblitz_instance;
use crypto_lib::cryptanalysis::pdp_split::{
    build_for_instance, exhaustive_decompositions, factor_base_size, run_hybrid,
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn default_plans(m: usize) -> Vec<Vec<usize>> {
    match m {
        2 => vec![vec![2]],
        3 => vec![vec![3], vec![2, 1]],
        4 => vec![vec![2, 2], vec![2, 1, 1]],
        _ => {
            let mut chain = vec![2];
            chain.extend(std::iter::repeat_n(1, m - 2));
            vec![chain]
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let flag = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let num = |name: &str, default: u64| {
        flag(name)
            .map(|v| {
                v.parse()
                    .unwrap_or_else(|_| panic!("{name} takes a number"))
            })
            .unwrap_or(default)
    };
    let cells: Vec<(usize, u8, u32, usize)> = flag("--cells")
        .expect("--cells m:a:n:ell[,...]")
        .split(',')
        .map(|c| {
            let p: Vec<u64> = c
                .split(':')
                .map(|x| x.parse().expect("cell is m:a:n:ell"))
                .collect();
            assert_eq!(p.len(), 4, "cell is m:a:n:ell");
            (p[0] as usize, p[1] as u8, p[2] as u32, p[3] as usize)
        })
        .collect();
    let guesses: Vec<usize> = flag("--guesses")
        .unwrap_or_else(|| "0".into())
        .split(',')
        .map(|g| g.parse().expect("--guesses g[,g...]"))
        .collect();
    let plan_override: Option<Vec<Vec<usize>>> = flag("--plans").map(|s| {
        s.split('|')
            .map(|p| {
                p.split('+')
                    .map(|x| x.parse().expect("plan is a+b+c"))
                    .collect()
            })
            .collect()
    });
    let targets = num("--targets", 4);
    let seed = num("--seed", 20261011);
    let budget = Duration::from_secs(num("--budget-secs", 120));
    let max_free = num("--max-free", 24) as usize;
    let lambda_targets = num("--lambda-targets", 256);
    let out_path = flag("--out").expect("--out FILE (must not exist)");
    let mut out = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&out_path)
        .unwrap_or_else(|e| panic!("{out_path}: {e} (runs are write-once)"));
    // `--replay-summary` prints, at the end, one JSON object of the
    // host-independent figures (solving degrees, F4 word XORs, verified
    // decompositions) for a `cairn` replay objective; timed-out solves are
    // host-dependent and are listed so a replay can refuse them.
    let replay_summary = args.iter().any(|a| a == "--replay-summary");
    let mut profile: Vec<String> = Vec::new();
    let (mut xors_total, mut verified_total, mut timeouts_total) = (0u64, 0usize, 0usize);
    let mut all_exact = true;
    let mut emit = |v: serde_json::Value| {
        writeln!(out, "{v}").expect("write");
        out.flush().expect("flush");
    };

    emit(serde_json::json!({
        "kind": "config",
        "schema": "pdp_split_dreg/v1",
        "cells": flag("--cells"),
        "guesses": guesses,
        "plans": flag("--plans"),
        "targets": targets,
        "seed": seed,
        "budget_secs": budget.as_secs(),
        "max_free": max_free,
        "lambda_targets": lambda_targets,
        "solver": "pq_f4_f2::groebner_basis_f4",
        "factor_base": "standard subspace span(1, z, ..., z^(l-1))",
    }));

    for &(m, a, n, ell) in &cells {
        let inst = koblitz_instance(a, n).unwrap_or_else(|| panic!("no K_{a} at n = {n}"));
        let words: Vec<u64> = (0..ell).map(|t| 1u64 << t).collect();
        let fb = factor_base_size(&inst, &words);
        let plans = plan_override.clone().unwrap_or_else(|| default_plans(m));
        let curve_id = inst.curve_id();
        // The decomposition yield λ (decompositions of a random subgroup
        // target over V), by exhaustive enumeration on its own target
        // stream: independent of every solver.
        let mut lrng = StdRng::seed_from_u64(
            seed ^ 0x4c41_4d42 ^ ((n as u64) << 40) ^ ((m * 64 + ell) as u64),
        );
        let mut lambda_hits = 0usize;
        let mut lambda_decomps = 0usize;
        for _ in 0..lambda_targets {
            let k = lrng.gen_range(1..inst.r);
            let x = inst.fast.mul_u64(inst.generator, k).x;
            let d = exhaustive_decompositions(&inst, &words, m, x).len();
            lambda_decomps += d;
            lambda_hits += usize::from(d > 0);
        }
        emit(serde_json::json!({
            "kind": "cell",
            "m": m, "a": a, "n": n, "ell": ell,
            "curve": inst.name,
            "curve_icv1": curve_id.icv1,
            "curve_slug": curve_id.slug,
            "group_order": inst.group_order,
            "r": inst.r,
            "cofactor": inst.cofactor,
            "factor_base_liftable_abscissae": fb,
            "surplus_n_minus_m_ell": n as i64 - (m * ell) as i64,
            "lambda_targets": lambda_targets,
            "lambda_decompositions": lambda_decomps,
            "lambda_targets_decomposable": lambda_hits,
            "plans": plans,
        }));
        eprintln!("cell m={m} K_{a}/GF(2^{n}) l={ell}: r={} fb={fb}", inst.r);
        let cell_seed = seed ^ ((m as u64) << 56) ^ ((n as u64) << 40) ^ ((ell as u64) << 32);
        let mut rng = StdRng::seed_from_u64(cell_seed);
        for t in 0..targets {
            let k = rng.gen_range(1..inst.r);
            let r_pt = inst.fast.mul_u64(inst.generator, k);
            let x_r = r_pt.x;
            let truth = exhaustive_decompositions(&inst, &words, m, x_r);
            for &g in &guesses {
                for plan in &plans {
                    if g > m * ell {
                        continue;
                    }
                    let (sys, build_ns) = match build_for_instance(&inst, plan, ell, x_r) {
                        Ok(s) => s,
                        Err(e) => {
                            emit(serde_json::json!({
                                "kind": "refused", "m": m, "n": n, "ell": ell,
                                "plan": plan, "g": g, "target": t, "reason": e,
                            }));
                            continue;
                        }
                    };
                    let run = run_hybrid(
                        &inst,
                        &sys,
                        &words,
                        x_r,
                        g,
                        Some(budget),
                        max_free,
                        &truth,
                        build_ns,
                    );
                    eprintln!(
                        "  t={t} {:<14} g={g}: vars={:>2} sd_max={} wall={:.3}s xors={} \
                         verified={}/{} exact={:?} timeouts={}",
                        run.plan,
                        run.n_vars,
                        run.solving_degree_max,
                        run.f4_wall_ns as f64 * 1e-9,
                        run.word_xors,
                        run.verified.len(),
                        truth.len(),
                        run.complete_and_exact,
                        run.branches_timed_out,
                    );
                    profile.push(format!(
                        "{m}:{n}:{ell}:{t}:{}:g{g}=sd{}/x{}/v{}/e{}",
                        run.plan,
                        run.solving_degree_max,
                        run.word_xors,
                        run.verified.len(),
                        run.complete_and_exact.map_or("-".into(), |e| e.to_string()),
                    ));
                    xors_total += run.word_xors;
                    verified_total += run.verified.len();
                    timeouts_total += run.branches_timed_out;
                    all_exact &= run.complete_and_exact == Some(true);
                    let mut v = serde_json::to_value(&run).expect("serialize");
                    let o = v.as_object_mut().expect("object");
                    o.insert("kind".into(), "run".into());
                    o.insert("m".into(), m.into());
                    o.insert("a".into(), a.into());
                    o.insert("n".into(), n.into());
                    o.insert("ell".into(), ell.into());
                    o.insert("target".into(), t.into());
                    o.insert("x_r".into(), x_r.into());
                    o.insert("truth_count".into(), truth.len().into());
                    emit(v);
                }
            }
        }
    }
    if replay_summary {
        println!(
            "{}",
            serde_json::json!({
                "runs": profile.len(),
                "profile": profile.join(";"),
                "word_xors_total": xors_total,
                "verified_total": verified_total,
                "branches_timed_out": timeouts_total,
                "all_complete_and_exact": all_exact,
            })
        );
    }
}
