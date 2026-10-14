//! **`C₄` with the `D₄` symmetry (candidate C8): one decomposition of a
//! point of `E(F_{p⁴})` over FGHR's `Y`-line base, measured and checked.**
//!
//! Companion to `research/gaudry_quartic_d4_20261013/`.  For each target:
//! the `D₄` solve of [`gaudry_quartic_d4::solve_d4`] (its `F_p`
//! multiplications by phase, the matrix shape, the quotient dimension) and
//! the independent meet-in-the-middle enumeration of the decompositions over
//! the same base; the two sets of decompositions must be equal.
//!
//! ```text
//! cargo run --release --example gaudry_quartic_d4 -- --p 269 --seed 1 --targets 5 --json OUT.jsonl
//! # decomposition rates of the Y-line base and of Gaudry's x-line base on
//! # the same curve, by meet in the middle over N random targets
//! cargo run --release --example gaudry_quartic_d4 -- --p 269 --seed 1 --hit-rate 4000
//! # targets planted as a signed sum of four distinct base points
//! cargo run --release --example gaudry_quartic_d4 -- --p 269 --seed 1 --targets 5 --planted --json OUT.jsonl
//! ```

use std::env;
use std::fs::OpenOptions;
use std::io::Write;
use std::time::Instant;

use crypto_lib::cryptanalysis::gaudry_quartic::{
    factor_base, mitm_decompositions, Pt4, QuarticSolve,
};
use crypto_lib::cryptanalysis::gaudry_quartic_d4::{
    generate_d4_instance, mitm_y_decompositions, solve_d4, y_base, SymmetrisedD4,
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

#[derive(Serialize)]
struct Row {
    p: u64,
    seed: u64,
    target: usize,
    planted: bool,
    macaulay_degree: u8,
    base_points: usize,
    solve: QuarticSolve,
    rational_solutions: usize,
    solver_decompositions: usize,
    mitm_decompositions: usize,
    decompositions_agree: bool,
    mitm_ms: f64,
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let get = |k: &str, d: u64| {
        args.iter()
            .position(|a| a == k)
            .map_or(d, |i| args[i + 1].parse().expect(k))
    };
    let p = get("--p", 269);
    let seed = get("--seed", 1);
    let targets = get("--targets", 3) as usize;
    let macaulay = get("--macaulay", 14) as u8;
    let debug = args.iter().any(|a| a == "--debug");
    let planted = args.iter().any(|a| a == "--planted");
    let out = args
        .iter()
        .position(|a| a == "--json")
        .map(|i| args[i + 1].clone());
    let inst = generate_d4_instance(p, seed);
    let mut rng = StdRng::seed_from_u64(seed ^ 0x55);
    let t = Instant::now();
    let pre = SymmetrisedD4::precompute(&inst, &mut rng, 8).expect("H is D₄-invariant");
    eprintln!(
        "p = {p}: D₄ form of S₅ on {} monomials, checked on {} points, {:.1} s",
        pre.monos.len(),
        pre.checked_points,
        t.elapsed().as_secs_f64()
    );
    let base = y_base(&inst, &mut rng);
    let hit_rate = get("--hit-rate", 0) as usize;
    if hit_rate > 0 {
        // The same targets for both bases.
        let xbase = factor_base(&inst.curve, &mut rng);
        let mut trng = StdRng::seed_from_u64(seed ^ 0x417);
        let (mut hy, mut hx, mut dy, mut dx) = (0usize, 0usize, 0usize, 0usize);
        for _ in 0..hit_rate {
            let r = loop {
                let x = inst.curve.f.random(&mut trng);
                if let Some(pt) = inst.curve.lift_x(&x, &mut trng) {
                    break pt;
                }
            };
            let ny = mitm_y_decompositions(&inst.curve, &base, &r).len();
            let nx = mitm_decompositions(&inst.curve, &xbase, &r).len();
            hy += usize::from(ny > 0);
            hx += usize::from(nx > 0);
            dy += ny;
            dx += nx;
        }
        let line = serde_json::json!({
            "p": p, "seed": seed, "targets": hit_rate,
            "y_base_values": base.len(), "x_base_values": xbase.len(),
            "y_hit_targets": hy, "x_hit_targets": hx,
            "y_decompositions": dy, "x_decompositions": dx,
        });
        println!("{line}");
        if let Some(path) = &out {
            let mut f = OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
                .expect("open --json");
            writeln!(f, "{line}").expect("write");
        }
        return;
    }
    for k in 0..targets {
        let r: Pt4 = if planted {
            loop {
                let mut idx: Vec<usize> = Vec::new();
                while idx.len() < 4 {
                    let i = rng.gen_range(0..base.len());
                    if !idx.contains(&i) {
                        idx.push(i);
                    }
                }
                let mut acc = Pt4 {
                    x: inst.curve.f.from_fp(0),
                    y: inst.curve.f.from_fp(0),
                    inf: true,
                };
                for &i in &idx {
                    let q = if rng.gen_bool(0.5) {
                        inst.curve.neg(&base[i].1)
                    } else {
                        base[i].1
                    };
                    acc = inst.curve.add(&acc, &q);
                }
                if !acc.inf {
                    break acc;
                }
            }
        } else {
            loop {
                let x = inst.curve.f.random(&mut rng);
                if let Some(pt) = inst.curve.lift_x(&x, &mut rng) {
                    break pt;
                }
            }
        };
        let mut srng = StdRng::seed_from_u64(seed + k as u64);
        let Some(sol) = solve_d4(&inst, &pre, &r, macaulay, &mut srng, debug) else {
            eprintln!(
                "target {k}: the Macaulay matrix at degree {macaulay} did not close the quotient"
            );
            continue;
        };
        let tm = Instant::now();
        let mitm = mitm_y_decompositions(&inst.curve, &base, &r);
        let mitm_ms = tm.elapsed().as_secs_f64() * 1e3;
        let mine: std::collections::BTreeSet<[u64; 4]> =
            sol.decompositions.iter().copied().collect();
        let row = Row {
            p,
            seed,
            target: k,
            planted,
            macaulay_degree: macaulay,
            base_points: base.len(),
            rational_solutions: sol.rational_solutions,
            solver_decompositions: mine.len(),
            mitm_decompositions: mitm.len(),
            decompositions_agree: mine == mitm,
            mitm_ms,
            solve: sol.solve,
        };
        let s = &row.solve;
        println!(
            "target {k}: {} × {} matrix, dim {}, {} rational solutions, decompositions {} (mitm {}, agree {}) | muls {:.3e} (echelon {:.2e}, nf {:.2e}, charpoly {:.2e}) | {:.1} s",
            s.rows,
            s.cols,
            s.dim,
            row.rational_solutions,
            row.solver_decompositions,
            row.mitm_decompositions,
            row.decompositions_agree,
            s.total_muls as f64,
            s.echelon_muls as f64,
            s.nf_muls as f64,
            s.charpoly_muls as f64,
            s.total_ms / 1e3
        );
        if let Some(path) = &out {
            let mut f = OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
                .expect("open --json");
            writeln!(f, "{}", serde_json::to_string(&row).expect("row")).expect("write");
        }
    }
}
