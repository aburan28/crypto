//! **The `D₄`-symmetrised `S₅` in four weighted unknowns against C8's
//! five-unknown form, on the same targets, with ground truth.**
//!
//! Companion to `research/d4_symmetrised_s5_20261015/`.  For each target
//! `R` of a curve over `F_{p⁴}` with a rational 2-torsion point (instance of
//! `gaudry_quartic_d4::generate_d4_instance`):
//!
//! * arm `w`: [`d4_weighted_s5::solve_weighted`] at weighted degree 26
//!   (`(σ₁, σ₂, σ₃, π)`, weights `(2, 2, 2, 1)`), then
//!   [`d4_weighted_s5::lift_back`], every multiplication counted;
//! * arm `c8`: `gaudry_quartic_d4::solve_d4` at degree 14 (C8's form,
//!   `(σ₁, σ₂, σ₃, σ₄, π)` with `π² = σ₄`);
//! * ground truth: exhaustive enumeration of four-subsets (`--brute`) and
//!   meet in the middle over pairs (`--mitm`) on the same base.
//!
//! ```text
//! d4_weighted_s5 --p 269 --seed 1 --targets 10 --mitm --json OUT.jsonl
//! d4_weighted_s5 --p 13 --seed 1 --targets 10 --brute --mitm --json OUT.jsonl
//! d4_weighted_s5 --generic 269 --seeds 3 --json OUT.jsonl     # solver positive control
//! d4_weighted_s5 --invariance-control 13 --seed 1             # negative control
//! ```
//!
//! Targets alternate planted (a signed sum of four distinct base points)
//! and uniform (a random point of `E(F_{p⁴})`).

use std::collections::BTreeSet;
use std::env;
use std::fs::OpenOptions;
use std::io::Write;
use std::time::Instant;

use crypto_lib::cryptanalysis::d4_weighted_s5::{
    brute_force_decompositions, lift_back, planted_generic, solve_weighted, weighted_system,
    LiftBack, MACAULAY_BOUND,
};
use crypto_lib::cryptanalysis::gaudry_quartic::{Pt4, QuarticSolve};
use crypto_lib::cryptanalysis::gaudry_quartic_d4::{
    generate_d4_instance, mitm_y_decompositions, solve_d4, y_base, D4Instance, SymmetrisedD4,
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
    base_points: usize,
    precompute_muls: u64,
    // arm w
    w_solved: bool,
    w_degree_used: u16,
    w_solve: QuarticSolve,
    w_rational_solutions: usize,
    w_lift: LiftBack,
    w_decompositions: usize,
    // arm c8
    c8_run: bool,
    c8_solved: bool,
    c8_solve: Option<QuarticSolve>,
    c8_rational_solutions: usize,
    c8_decompositions: usize,
    // ground truth
    brute: Option<usize>,
    brute_ms: Option<f64>,
    mitm: Option<usize>,
    mitm_ms: Option<f64>,
    w_equals_brute: Option<bool>,
    w_equals_mitm: Option<bool>,
    c8_equals_w: Option<bool>,
    planted_found_w: Option<bool>,
    planted_found_c8: Option<bool>,
}

fn append(path: &Option<String>, line: &str) {
    if let Some(path) = path {
        let mut f = OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .expect("open --json");
        writeln!(f, "{line}").expect("write");
    }
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let get = |k: &str, d: u64| {
        args.iter()
            .position(|a| a == k)
            .map_or(d, |i| args[i + 1].parse().expect(k))
    };
    let has = |k: &str| args.iter().any(|a| a == k);
    let out = args
        .iter()
        .position(|a| a == "--json")
        .map(|i| args[i + 1].clone());
    let seed = get("--seed", 1);
    let debug = has("--debug");
    let no_selection = has("--no-row-selection");

    // Solver positive control: planted generic weighted systems.
    let gp = get("--generic", 0);
    if gp > 0 {
        for s in 0..get("--seeds", 3) {
            let (comps, root) = planted_generic(gp, seed + s);
            let mut rng = StdRng::seed_from_u64(s);
            let mut st = QuarticSolve::default();
            let (sols, deg) = solve_weighted(
                &comps,
                MACAULAY_BOUND,
                gp,
                &mut rng,
                &mut st,
                !no_selection,
                debug,
            );
            let found = sols.as_ref().is_some_and(|v| v.contains(&root));
            let line = serde_json::json!({
                "mode": "generic", "p": gp, "seed": seed + s, "degree_used": deg,
                "dim": st.dim, "rows": st.rows, "cols": st.cols,
                "total_muls": st.total_muls, "echelon_muls": st.echelon_muls,
                "nf_muls": st.nf_muls, "charpoly_muls": st.charpoly_muls,
                "planted_root_found": found, "solutions": sols.map(|v| v.len()),
                "total_ms": st.total_ms,
            });
            println!("{line}");
            append(&out, &line.to_string());
        }
        return;
    }

    // Negative control for the invariance check: a wrong `c` (not a square
    // root of f'(x₀)) gives a coordinate that τ_T does not negate, and H must
    // then fall outside the span of the D₄ invariants.
    let ic = get("--invariance-control", 0);
    if ic > 0 {
        let inst = generate_d4_instance(ic, seed);
        let mut rng = StdRng::seed_from_u64(seed ^ 0x55);
        let good = SymmetrisedD4::precompute(&inst, &mut rng, 8).is_ok();
        let bad_c = (inst.c % (ic - 1)) + 1;
        let bad_c = if bad_c * bad_c % ic == inst.c * inst.c % ic {
            (bad_c % (ic - 1)) + 1
        } else {
            bad_c
        };
        let wrong = D4Instance {
            curve: inst.curve.clone(),
            x0: inst.x0,
            c: bad_c,
        };
        let bad = SymmetrisedD4::precompute(&wrong, &mut rng, 8);
        let line = serde_json::json!({
            "mode": "invariance_control", "p": ic, "seed": seed,
            "true_c_in_span": good, "wrong_c": bad_c,
            "wrong_c_in_span": bad.is_ok(),
            "wrong_c_error": bad.err(),
        });
        println!("{line}");
        append(&out, &line.to_string());
        return;
    }

    let p = get("--p", 269);
    let targets = get("--targets", 4) as usize;
    let degree = get("--degree", MACAULAY_BOUND as u64) as u16;
    let run_c8 = !has("--no-c8");
    let brute = has("--brute");
    let mitm = has("--mitm");
    let inst = generate_d4_instance(p, seed);
    let mut rng = StdRng::seed_from_u64(seed ^ 0x55);
    let t = Instant::now();
    let pre = SymmetrisedD4::precompute(&inst, &mut rng, 8).expect("H is D₄-invariant");
    eprintln!(
        "p = {p}, seed {seed}: D₄ form of S₅ on {} monomials, checked on {} points, {:.1} s",
        pre.monos.len(),
        pre.checked_points,
        t.elapsed().as_secs_f64()
    );
    let base = y_base(&inst, &mut rng);
    let mut trng = StdRng::seed_from_u64(seed ^ 0x7a12);
    for k in 0..targets {
        let planted = k % 2 == 0;
        let mut planted_ys: Option<[u64; 4]> = None;
        let r: Pt4 = if planted {
            loop {
                let mut idx: Vec<usize> = Vec::new();
                while idx.len() < 4 {
                    let i = trng.gen_range(0..base.len());
                    if !idx.contains(&i) {
                        idx.push(i);
                    }
                }
                let mut acc = Pt4::INF;
                for &i in &idx {
                    let q = if trng.gen_bool(0.5) {
                        inst.curve.neg(&base[i].1)
                    } else {
                        base[i].1
                    };
                    acc = inst.curve.add(&acc, &q);
                }
                if !acc.inf {
                    let mut ys: [u64; 4] = core::array::from_fn(|j| base[idx[j]].0);
                    ys.sort_unstable();
                    planted_ys = Some(ys);
                    break acc;
                }
            }
        } else {
            loop {
                let x = inst.curve.f.random(&mut trng);
                if let Some(pt) = inst.curve.lift_x(&x, &mut trng) {
                    break pt;
                }
            }
        };
        // Arm w.
        let f = &inst.curve.f;
        let mut st = QuarticSolve::default();
        let before = f.muls();
        let y_r = inst.y_of_x(&r.x);
        let comps = weighted_system(&pre, f, &y_r);
        st.weil_muls = f.muls() - before;
        let mut srng = StdRng::seed_from_u64(seed * 1000 + k as u64);
        let (sols, deg_used) =
            solve_weighted(&comps, degree, p, &mut srng, &mut st, !no_selection, debug);
        let w_solved = sols.is_some();
        let sols = sols.unwrap_or_default();
        let lb = lift_back(&inst, &sols, &r, &mut srng);
        let w_set: BTreeSet<[u64; 4]> = lb.decompositions.iter().copied().collect();
        // Arm c8.
        let (c8_solved, c8_solve, c8_rat, c8_set) = if run_c8 {
            let mut crng = StdRng::seed_from_u64(seed * 1000 + k as u64);
            match solve_d4(&inst, &pre, &r, 14, &mut crng, debug) {
                Some(s) => (
                    true,
                    Some(s.solve),
                    s.rational_solutions,
                    s.decompositions.iter().copied().collect::<BTreeSet<_>>(),
                ),
                None => (false, None, 0, BTreeSet::new()),
            }
        } else {
            (false, None, 0, BTreeSet::new())
        };
        // Ground truth.
        let (bset, bms) = if brute {
            let tb = Instant::now();
            let s = brute_force_decompositions(&inst.curve, &base, &r);
            (Some(s), Some(tb.elapsed().as_secs_f64() * 1e3))
        } else {
            (None, None)
        };
        let (mset, mms) = if mitm {
            let tm = Instant::now();
            let s = mitm_y_decompositions(&inst.curve, &base, &r);
            (Some(s), Some(tm.elapsed().as_secs_f64() * 1e3))
        } else {
            (None, None)
        };
        let row = Row {
            p,
            seed,
            target: k,
            planted,
            base_points: base.len(),
            precompute_muls: pre.precompute_muls,
            w_solved,
            w_degree_used: deg_used,
            w_rational_solutions: sols.len(),
            w_decompositions: w_set.len(),
            w_lift: lb,
            w_solve: st,
            c8_run: run_c8,
            c8_solved,
            c8_rational_solutions: c8_rat,
            c8_decompositions: c8_set.len(),
            c8_solve,
            brute: bset.as_ref().map(|s| s.len()),
            brute_ms: bms,
            mitm: mset.as_ref().map(|s| s.len()),
            mitm_ms: mms,
            w_equals_brute: bset.as_ref().map(|s| *s == w_set),
            w_equals_mitm: mset.as_ref().map(|s| *s == w_set),
            c8_equals_w: if run_c8 && c8_solved {
                Some(c8_set == w_set)
            } else {
                None
            },
            planted_found_w: planted_ys.map(|ys| w_set.contains(&ys)),
            planted_found_c8: if run_c8 {
                planted_ys.map(|ys| c8_set.contains(&ys))
            } else {
                None
            },
        };
        println!(
            "p {p} target {k} ({}): w deg {} {}×{} dim {} C {:.3e} lift {:.2e} sols {} dec {} | c8 {} C {} dec {} | brute {:?} mitm {:?} | w=brute {:?} w=mitm {:?} c8=w {:?} | {:.1} s",
            if planted { "planted" } else { "uniform" },
            row.w_degree_used,
            row.w_solve.rows,
            row.w_solve.cols,
            row.w_solve.dim,
            row.w_solve.total_muls as f64,
            row.w_lift.muls as f64,
            row.w_rational_solutions,
            row.w_decompositions,
            row.c8_solved,
            row.c8_solve
                .as_ref()
                .map_or("-".to_string(), |s| format!("{:.3e} dim {}", s.total_muls as f64, s.dim)),
            row.c8_decompositions,
            row.brute,
            row.mitm,
            row.w_equals_brute,
            row.w_equals_mitm,
            row.c8_equals_w,
            row.w_solve.total_ms / 1e3
        );
        append(&out, &serde_json::to_string(&row).expect("row"));
    }
}
