//! **Shape replay for the `k = 4` Macaulay solve: does it cut `C₄`?**
//!
//! Companion to `research/gaudry_quartic_replay_20261013/` and
//! `research/ic_speed_and_rho_scope_20261013/` §3–4.  Every trial of
//! Gaudry's `k = 4` index calculus solves a system of the same monomial
//! support, so the Macaulay matrix has the same shape each time: the rows
//! that reduce to zero on the first system reduce to zero on the next.
//! `solve_system4_traced` records the shape on a first system and, on each
//! later one, does not build those rows; a replay whose pivots or quotient
//! dimension differ falls back to the full solve, charged.
//!
//! For each system this runs the full solve (the baseline, the counter
//! §11.17 measured `C₄` with) and the replay, on the same coefficients, and
//! writes both counts, both wall times and whether the solution sets agree.
//!
//! ```text
//! # generic planted systems of degree d (the §11.17 scaling panel's family)
//! cargo run --release --example gaudry_quartic_replay -- --generic 3,4,5,6 --systems 5 --json OUT.jsonl
//! # the S₅ systems of real residuals (degree 8; about 30 min a solve)
//! cargo run --release --example gaudry_quartic_replay -- --s5 --p 269 --seed 1 --systems 3 --json OUT.jsonl
//! ```

use std::collections::HashMap;
use std::env;
use std::fs::OpenOptions;
use std::io::Write;

use crypto_lib::cryptanalysis::gaudry_quartic::{
    generate_instance4, solve_system4_traced, QuarticSolve, ReplayOutcome, ShapeTrace,
    SymmetrisedS5,
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

#[derive(Serialize)]
struct Row {
    family: String,
    equation_degree: u8,
    system: usize,
    baseline: QuarticSolve,
    replay: QuarticSolve,
    replay_outcome: ReplayOutcome,
    solutions_agree: bool,
    planted_found_baseline: Option<bool>,
    planted_found_replay: Option<bool>,
    baseline_solutions: usize,
    replay_solutions: usize,
    trace_zero_rows: usize,
    ratio_total_muls: f64,
    ratio_wall: f64,
}

/// A generic planted system of degree `d`: every monomial of degree at most
/// `d`, random coefficients, constant terms set so that `root` is a zero.
fn generic_system(d: u8, p: u64, rng: &mut StdRng) -> (Vec<HashMap<[u8; 4], u64>>, [u64; 4]) {
    let root: [u64; 4] = std::array::from_fn(|_| rng.gen_range(0..p));
    let mut monos = Vec::new();
    for a in 0..=d {
        for b in 0..=(d - a) {
            for c in 0..=(d - a - b) {
                for e in 0..=(d - a - b - c) {
                    monos.push([a, b, c, e]);
                }
            }
        }
    }
    let eval = |c: &HashMap<[u8; 4], u64>| -> u64 {
        c.iter().fold(0u64, |acc, (m, &k)| {
            let mut t = k;
            for i in 0..4 {
                for _ in 0..m[i] {
                    t = t * root[i] % p;
                }
            }
            (acc + t) % p
        })
    };
    let comps = (0..4)
        .map(|_| {
            let mut c: HashMap<[u8; 4], u64> =
                monos.iter().map(|m| (*m, rng.gen_range(1..p))).collect();
            let v = eval(&c);
            let e = c.entry([0, 0, 0, 0]).or_insert(0);
            *e = (*e + p - v) % p;
            c
        })
        .collect();
    (comps, root)
}

fn emit(out: &Option<String>, row: &Row) {
    let line = serde_json::to_string(row).expect("row serialises");
    println!(
        "{} d={} system {}: {:?} | muls {:.3e} → {:.3e} ({:.3}×) | wall {:.1} → {:.1} ms ({:.3}×) | zero rows skipped {} | agree {}",
        row.family,
        row.equation_degree,
        row.system,
        row.replay_outcome,
        row.baseline.total_muls as f64,
        row.replay.total_muls as f64,
        row.ratio_total_muls,
        row.baseline.total_ms,
        row.replay.total_ms,
        row.ratio_wall,
        row.trace_zero_rows,
        row.solutions_agree
    );
    if let Some(path) = out {
        let mut f = OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .expect("open --json");
        writeln!(f, "{line}").expect("write --json");
    }
}

#[allow(clippy::too_many_arguments)]
fn run_pair(
    family: &str,
    d: u8,
    system: usize,
    comps: &[HashMap<[u8; 4], u64>],
    macaulay: u8,
    p: u64,
    seed: u64,
    root: Option<[u64; 4]>,
    trace: &mut Option<ShapeTrace>,
    out: &Option<String>,
) {
    // Baseline: the full solve, exactly as `solve_system4`.  The first
    // system's baseline also records the trace the later ones replay.
    let mut base_st = QuarticSolve::default();
    let mut rng = StdRng::seed_from_u64(seed);
    let mut recorded = ShapeTrace::default();
    let first = trace.is_none();
    let base = solve_system4_traced(
        comps,
        d,
        macaulay,
        p,
        &mut rng,
        &mut base_st,
        false,
        None,
        if first { Some(&mut recorded) } else { None },
    );
    if first {
        *trace = Some(recorded);
    }
    let t = trace.as_ref().expect("recorded above");
    let mut rep_st = QuarticSolve::default();
    let mut rng = StdRng::seed_from_u64(seed);
    let rep = solve_system4_traced(
        comps,
        d,
        macaulay,
        p,
        &mut rng,
        &mut rep_st,
        false,
        Some(t),
        None,
    );
    let (bs, rs) = (
        base.clone().unwrap_or_default(),
        rep.clone().unwrap_or_default(),
    );
    let row = Row {
        family: family.into(),
        equation_degree: d,
        system,
        replay_outcome: rep_st.replay,
        solutions_agree: base.is_some() && bs == rs,
        planted_found_baseline: root.map(|r| bs.contains(&r)),
        planted_found_replay: root.map(|r| rs.contains(&r)),
        baseline_solutions: bs.len(),
        replay_solutions: rs.len(),
        trace_zero_rows: t.zero_ids.len(),
        ratio_total_muls: rep_st.total_muls as f64 / base_st.total_muls.max(1) as f64,
        ratio_wall: rep_st.total_ms / base_st.total_ms.max(1e-9),
        baseline: base_st,
        replay: rep_st,
    };
    emit(out, &row);
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut generic: Vec<u8> = Vec::new();
    let mut s5 = false;
    let mut p = 269u64;
    let mut seed = 1u64;
    let mut systems = 5usize;
    let mut out: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--generic" => {
                i += 1;
                generic = args[i]
                    .split(',')
                    .map(|v| v.parse().expect("--generic"))
                    .collect();
            }
            "--s5" => s5 = true,
            "--p" => {
                i += 1;
                p = args[i].parse().expect("--p");
            }
            "--seed" => {
                i += 1;
                seed = args[i].parse().expect("--seed");
            }
            "--systems" => {
                i += 1;
                systems = args[i].parse().expect("--systems");
            }
            "--json" => {
                i += 1;
                out = Some(args[i].clone());
            }
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    for d in generic {
        let mut rng = StdRng::seed_from_u64(seed ^ ((d as u64) << 32));
        let mut trace = None;
        for k in 0..systems {
            let (comps, root) = generic_system(d, 269, &mut rng);
            run_pair(
                "generic",
                d,
                k,
                &comps,
                4 * (d - 1) + 1,
                269,
                seed + k as u64,
                Some(root),
                &mut trace,
                &out,
            );
        }
    }
    if s5 {
        let inst = generate_instance4(p, seed);
        let curve = &inst.curve;
        let mut rng = StdRng::seed_from_u64(seed ^ 0x55);
        let pre = SymmetrisedS5::precompute(curve, &mut rng);
        let mut trace = None;
        for k in 0..systems {
            let r = loop {
                let x = curve.f.random(&mut rng);
                if let Some(pt) = curve.lift_x(&x, &mut rng) {
                    break pt;
                }
            };
            let comps = pre.weil_restrict(&curve.f, &r.x);
            run_pair(
                "s5",
                8,
                k,
                &comps[..],
                29,
                p,
                seed + k as u64,
                None,
                &mut trace,
                &out,
            );
        }
    }
}
