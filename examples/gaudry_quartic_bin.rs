//! `C₄` in characteristic 2: the `D₄`-symmetric `k = 4` decomposition on
//! `E/F_{q⁴}`, `q = 2^m`, over the `t`-line base
//! ([`gaudry_quartic_bin`]), every solve checked against meet in the
//! middle.  One JSON row per closed solve with `--json`.
//!
//! ```sh
//! cargo run --release --example gaudry_quartic_bin -- --m 6 --seed 1 --targets 5 --planted --weighted --sparsest --replay
//! ```

use std::env;
use std::fs::OpenOptions;
use std::io::Write;
use std::time::Instant;

use crypto_lib::cryptanalysis::gaudry_quartic_bin::{
    find_bin_curve, generate_instance, mitm_t_decompositions, solve_d4_bin, t_base, BinOptions,
    BinSolve, BinTrace, PtB, SymmetrisedBin,
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

#[derive(Serialize)]
struct Row {
    m: u32,
    seed: u64,
    target: usize,
    planted: bool,
    macaulay_degree: u8,
    form: BinOptions,
    base_points: usize,
    solve: BinSolve,
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
    let flag = |k: &str| args.iter().any(|a| a == k);
    let m = get("--m", 6) as u32;
    let seed = get("--seed", 1);
    let targets = get("--targets", 3) as usize;
    let form = BinOptions {
        weighted: flag("--weighted") || flag("--frobenius"),
        frobenius: flag("--frobenius"),
        sparsest: flag("--sparsest"),
        replay: flag("--replay"),
    };
    let macaulay = get(
        "--macaulay",
        if form.frobenius {
            11
        } else if form.weighted {
            26
        } else {
            14
        },
    ) as u8;
    let planted = flag("--planted");
    let out = args
        .iter()
        .position(|a| a == "--json")
        .map(|i| args[i + 1].clone());
    let find = get("--find", 0);
    let curve = if find > 0 {
        let c = find_bin_curve(m, seed, find).expect("a curve");
        eprintln!(
            "f2m4({m}, {seed}, {find}): instance seed {}, #E = {} = {} · {}",
            c.instance_seed, c.group_order, c.cofactor, c.r
        );
        c.curve
    } else {
        generate_instance(m, seed)
    };
    let mut rng = StdRng::seed_from_u64(seed ^ 0x55);
    let t = Instant::now();
    let pre = SymmetrisedBin::precompute(&curve, &mut rng, 8).expect("H is D₄-invariant");
    eprintln!(
        "m = {m}: D₄ form of S₅ on {} monomials, checked on {} points, {:.3e} F_q muls, {:.1} s",
        pre.monos.len(),
        pre.checked_points,
        pre.precompute_muls as f64,
        t.elapsed().as_secs_f64()
    );
    eprintln!(
        "Frobenius-reducible (even powers of s₂, s₃, π): {}",
        pre.frobenius_reducible()
    );
    let base = t_base(&curve, &mut rng);
    let mut trace = BinTrace::default();
    for k in 0..targets {
        let r: PtB = if planted {
            loop {
                let mut idx: Vec<usize> = Vec::new();
                while idx.len() < 4 {
                    let i = rng.gen_range(0..base.len());
                    if !idx.contains(&i) {
                        idx.push(i);
                    }
                }
                let mut acc = PtB::INF;
                for &i in &idx {
                    let q = if rng.gen_bool(0.5) {
                        curve.neg(&base[i].1)
                    } else {
                        base[i].1
                    };
                    acc = curve.add(&acc, &q);
                }
                if !acc.inf {
                    break acc;
                }
            }
        } else {
            loop {
                let x = curve.f.random(&mut rng);
                if let Some(p) = curve.lift_x(&x, &mut rng) {
                    break p;
                }
            }
        };
        let mut srng = StdRng::seed_from_u64(seed + k as u64);
        let sol = match solve_d4_bin(
            &curve,
            &pre,
            &r,
            form,
            macaulay,
            &mut srng,
            form.replay.then_some(&mut trace),
        ) {
            Ok(s) => s,
            Err(cost) => {
                let t_r = curve.t_of_x(&r.x);
                let mitm = mitm_t_decompositions(&curve, &base, &r);
                let paired = mitm
                    .iter()
                    .filter(|d| (0..4).any(|i| (0..4).any(|j| d[i] ^ 1 == d[j])))
                    .count();
                eprintln!(
                    "target {k}: mitm {} decompositions, {} with a pair {{t, t + 1}}",
                    mitm.len(),
                    paired
                );
                eprintln!(
                    "target {k}: did not close ({}, dim {}, {:.3e} muls; t_R in F_q {}, in F_q² {})",
                    cost.outcome,
                    cost.dim,
                    cost.total_muls as f64,
                    t_r.in_fq(),
                    curve.f.in_fq2(&t_r)
                );
                continue;
            }
        };
        let tm = Instant::now();
        let mitm = if flag("--no-mitm") {
            sol.decompositions.iter().copied().collect()
        } else {
            mitm_t_decompositions(&curve, &base, &r)
        };
        let mitm_ms = tm.elapsed().as_secs_f64() * 1e3;
        let mine: std::collections::BTreeSet<[u16; 4]> =
            sol.decompositions.iter().copied().collect();
        let row = Row {
            m,
            seed,
            target: k,
            planted,
            macaulay_degree: macaulay,
            form,
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
            "target {k}: {} × {} matrix, dim {}, {} chains, {} rational solutions, decompositions {} (mitm {}, agree {}) | muls {:.3e} (echelon {:.2e}, nf {:.2e}, det A {:.2e}, roots {:.2e}) | {} | {:.2} s",
            s.rows,
            s.cols,
            s.dim,
            s.chains,
            row.rational_solutions,
            row.solver_decompositions,
            row.mitm_decompositions,
            row.decompositions_agree,
            s.total_muls as f64,
            s.echelon_muls as f64,
            s.nf_muls as f64,
            s.charpoly_muls as f64,
            s.roots_muls as f64,
            s.replay,
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
    if form.replay {
        eprintln!(
            "shape replays: {} replayed, {} fell back",
            trace.replays, trace.fallbacks
        );
    }
}
