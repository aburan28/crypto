//! **Joux–Vitse `(n − 1)`-decomposition at `n = 5` on `E(F_{p⁵})`: the
//! cost of one trial and of one relation, measured and checked.**
//!
//! Companion to `research/joux_vitse_n1_20261015/` and
//! [`crypto_lib::cryptanalysis::joux_vitse`].
//!
//! ```text
//! # one rung of the ladder: F4 trials (each cross-checked by meet in the
//! # middle), then MITM-only trials for the decomposition rate
//! cargo run --release --example joux_vitse_ladder -- trials --p 31 --seed 1 \
//!     --f4-trials 64 --mitm-trials 20000 --json OUT.jsonl
//! # a complete DLP at toy size, paired with Pollard rho on the same target
//! cargo run --release --example joux_vitse_ladder -- dlp --p 11 --seed 1 --json OUT.jsonl
//! ```
//!
//! Every relation is re-added in the group; every F4 trial's set of
//! decompositions into four distinct base points must equal the
//! meet-in-the-middle set.  Costs are `F_p` multiplications.

use std::collections::{BTreeMap, HashMap};
use std::env;
use std::fs::OpenOptions;
use std::io::Write;
use std::time::Instant;

use crypto_lib::cryptanalysis::f4_fp;
use crypto_lib::cryptanalysis::joux_vitse::{
    factor_base, mitm_decompositions, rho, solve_logs, trial, trial_options, verify_relation,
    Curve5, PairTable, Pt5, SymmetrisedS5x5, TrialReport,
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::Serialize;
use serde_json::json;

fn arg<T: std::str::FromStr>(args: &[String], k: &str, d: T) -> T {
    args.iter()
        .position(|a| a == k)
        .map_or(d, |i| args[i + 1].parse().ok().expect(k))
}

fn emit(path: &Option<String>, v: &serde_json::Value) {
    println!("{v}");
    if let Some(p) = path {
        let mut f = OpenOptions::new()
            .create(true)
            .append(true)
            .open(p)
            .expect("open json");
        writeln!(f, "{v}").unwrap();
    }
}

struct Setup {
    curve: Curve5,
    sym: SymmetrisedS5x5,
    base: Vec<Pt5>,
    index: HashMap<u64, usize>,
    pairs: PairTable,
    trace: Option<f4_fp::F4Trace>,
    setup: serde_json::Value,
}

/// Curve, symmetrised `S₅`, factor base, pair table, and the F4 trace
/// recorded on the first inconsistent system (Joux–Vitse's precomputation).
fn setup(p: u64, seed: u64, rng: &mut StdRng, record_trace: bool) -> Setup {
    let t0 = Instant::now();
    let curve = Curve5::random(p, rng);
    let m0 = curve.f.muls();
    let sym = SymmetrisedS5x5::precompute(&curve, rng);
    let sym_ms = t0.elapsed().as_secs_f64() * 1e3;
    let m1 = curve.f.muls();
    let base = factor_base(&curve);
    let base_muls = curve.f.muls() - m1;
    let index: HashMap<u64, usize> = base
        .iter()
        .enumerate()
        .map(|(i, pt)| (pt.x.0[0], i))
        .collect();
    let m2 = curve.f.muls();
    let pairs = PairTable::new(&curve, &base);
    let pair_muls = curve.f.muls() - m2;
    let opts = trial_options();
    let (mut trace, mut trace_muls, mut trace_tries) = (None, 0u64, 0u64);
    let t1 = Instant::now();
    if record_trace {
        while trace.is_none() {
            trace_tries += 1;
            let r = random_point(&curve, rng);
            let sys = sym.weil_restrict(&curve.f, &r.x);
            let (rep, tr) = f4_fp::f4_record(&sys, 4, p, &opts);
            trace_muls += rep.field_ops;
            if rep.inconsistent {
                trace = Some(tr);
            }
        }
    }
    let trace_ms = t1.elapsed().as_secs_f64() * 1e3;
    let setup = json!({
        "kind": "setup", "p": p, "seed": seed,
        "a": curve.a.0, "b": curve.b.0, "c": curve.f.c,
        "base_points": base.len(),
        "sym_precompute_muls": sym.precompute_muls, "sym_ms": sym_ms,
        "curve_muls_before_sym": m0,
        "base_muls": base_muls, "pair_table_muls": pair_muls,
        "trace_record_muls": trace_muls, "trace_record_tries": trace_tries,
        "trace_ms": trace_ms,
        "trace_steps": trace.as_ref().map(|t| t.steps.len()),
    });
    Setup {
        curve,
        sym,
        base,
        index,
        pairs,
        trace,
        setup,
    }
}

fn random_point(curve: &Curve5, rng: &mut StdRng) -> Pt5 {
    loop {
        if let Some(pt) = curve.lift_x(&curve.f.random(rng)) {
            return pt;
        }
    }
}

#[derive(Serialize)]
struct TrialRow<'a> {
    kind: &'static str,
    p: u64,
    seed: u64,
    trial: usize,
    report: &'a TrialReport,
    mitm_decompositions: usize,
    agree: bool,
    relations_verified: bool,
}

fn trials_mode(args: &[String]) {
    let p: u64 = arg(args, "--p", 31);
    let seed: u64 = arg(args, "--seed", 1);
    let n_f4: usize = arg(args, "--f4-trials", 16);
    let n_mitm: usize = arg(args, "--mitm-trials", 0);
    let replay = arg(args, "--replay", 1u8) == 1;
    let json: Option<String> = args
        .iter()
        .position(|a| a == "--json")
        .map(|i| args[i + 1].clone());
    let mut rng = StdRng::seed_from_u64(seed ^ (p << 20));
    let s = setup(p, seed, &mut rng, replay);
    emit(&json, &s.setup);
    let opts = trial_options();
    // Targets: R = [a]G + [b]H for random G, H on the curve.
    let g = random_point(&s.curve, &mut rng);
    let h = random_point(&s.curve, &mut rng);
    let targets: Vec<(u128, u128)> = (0..n_f4 + n_mitm)
        .map(|_| (rng.gen_range(0..1u128 << 64), rng.gen_range(0..1u128 << 64)))
        .collect();
    let (a, b, params) = (s.curve.a, s.curve.b, (p, s.curve.f.c));
    let _ = params;
    // F4 trials, sequential (each F4 runs its own threads), each checked.
    let t0 = Instant::now();
    let mut agg = BTreeMap::<&str, u64>::new();
    for (k, &(ta, tb)) in targets[..n_f4].iter().enumerate() {
        let c = Curve5::with_params(p, a, b);
        let r = c.add(&c.mul(&g, ta), &c.mul(&h, tb));
        let rep = trial(&c, &s.sym, &s.base, &s.index, &r, &opts, s.trace.as_ref());
        let mitm = mitm_decompositions(&c, &s.base, &s.pairs, &r);
        let agree = mitm == rep.decompositions;
        let verified = rep
            .relations
            .iter()
            .all(|rel| verify_relation(&c, &s.base, rel, &r));
        *agg.entry("muls_system").or_default() += rep.muls_system;
        *agg.entry("muls_f4").or_default() += rep.muls_f4;
        *agg.entry("muls_solve").or_default() += rep.muls_solve;
        *agg.entry("consistent").or_default() += rep.consistent as u64;
        *agg.entry("replayed").or_default() += rep.replayed as u64;
        *agg.entry("replay_mismatch").or_default() += rep.replay_mismatch as u64;
        *agg.entry("relations").or_default() += rep.relations.len() as u64;
        *agg.entry("trials_with_relation").or_default() += (!rep.relations.is_empty()) as u64;
        *agg.entry("disagree").or_default() += (!agree) as u64;
        *agg.entry("unverified").or_default() += (!verified) as u64;
        *agg.entry("undetermined").or_default() += rep.undetermined as u64;
        *agg.entry("solutions").or_default() += rep.split.solutions as u64;
        *agg.entry("not_split").or_default() += rep.split.not_split as u64;
        let row = TrialRow {
            kind: "f4_trial",
            p,
            seed,
            trial: k,
            report: &rep,
            mitm_decompositions: mitm.len(),
            agree,
            relations_verified: verified,
        };
        emit(&json, &serde_json::to_value(&row).unwrap());
    }
    let f4_ms = t0.elapsed().as_secs_f64() * 1e3;
    // MITM-only trials, in parallel: the decomposition rate.
    let t1 = Instant::now();
    let hits: Vec<usize> = targets[n_f4..]
        .par_iter()
        .map(|&(ta, tb)| {
            let c = Curve5::with_params(p, a, b);
            let r = c.add(&c.mul(&g, ta), &c.mul(&h, tb));
            mitm_decompositions(&c, &s.base, &s.pairs, &r).len()
        })
        .collect();
    let mitm_ms = t1.elapsed().as_secs_f64() * 1e3;
    let mitm_hit_trials = hits.iter().filter(|&&h| h > 0).count();
    let mitm_decomps: usize = hits.iter().sum();
    emit(
        &json,
        &json!({
            "kind": "summary", "p": p, "seed": seed,
            "base_points": s.base.len(),
            "f4_trials": n_f4, "f4_agg": agg, "f4_wall_ms": f4_ms,
            "mitm_trials": n_mitm, "mitm_hit_trials": mitm_hit_trials,
            "mitm_decompositions": mitm_decomps, "mitm_wall_ms": mitm_ms,
            "threads": rayon::current_num_threads(),
        }),
    );
}

/// Smallest prime factor; the group orders here are below `2^64`.
fn factor(mut n: u128) -> Vec<u128> {
    let mut out = Vec::new();
    let mut d = 2u128;
    while d * d <= n {
        while n % d == 0 {
            out.push(d);
            n /= d;
        }
        d += 1;
    }
    if n > 1 {
        out.push(n);
    }
    out
}

fn dlp_mode(args: &[String]) {
    let p: u64 = arg(args, "--p", 11);
    let seed: u64 = arg(args, "--seed", 1);
    let max_trials: usize = arg(args, "--max-trials", 200_000);
    let min_ell_bits: u32 = arg(args, "--min-ell-bits", 10);
    let json: Option<String> = args
        .iter()
        .position(|a| a == "--json")
        .map(|i| args[i + 1].clone());
    let mut rng = StdRng::seed_from_u64(seed ^ (p << 20) ^ 0xd1);
    // Draw curves until the group has a prime factor ℓ of the wanted size.
    let mut curve_tries = 0;
    let (s, order, ell) = loop {
        curve_tries += 1;
        let s = setup(p, seed, &mut rng, true);
        let order = s.curve.order_by_count();
        let ell = *factor(order).last().unwrap();
        if 128 - ell.leading_zeros() >= min_ell_bits && ell > (s.base.len() as u128) * 4 {
            break (s, order, ell);
        }
    };
    let h = order / ell;
    emit(&json, &s.setup);
    let c = &s.curve;
    // P of order ℓ, Q = kP with a planted k.
    let pp = loop {
        let g = random_point(c, &mut rng);
        let q = c.mul(&g, h);
        if !q.inf {
            break q;
        }
    };
    assert!(c.mul(&pp, ell).inf, "P has order ℓ");
    let k_secret = rng.gen_range(1..ell);
    let qq = c.mul(&pp, k_secret);
    emit(
        &json,
        &json!({"kind": "instance", "p": p, "seed": seed, "order": order.to_string(),
                "ell": ell.to_string(), "cofactor": h.to_string(),
                "curve_tries": curve_tries, "P": pp, "Q": qq,
                "k_secret": k_secret.to_string()}),
    );
    // Rho on the same target, counted in group operations and F_p mults.
    let rc = Curve5::with_params(p, c.a, c.b);
    let mut rrng = StdRng::seed_from_u64(seed ^ 0x5eed);
    let m0 = rc.f.muls();
    let tr = Instant::now();
    let (k_rho, rho_ops) = rho(&rc, &pp, &qq, ell, &mut rrng);
    let rho_muls = rc.f.muls() - m0;
    let rho_ms = tr.elapsed().as_secs_f64() * 1e3;
    let rho_ok = k_rho.map(|k| rc.mul(&pp, k) == qq).unwrap_or(false);
    emit(
        &json,
        &json!({"kind": "rho", "p": p, "seed": seed, "k": k_rho.map(|k| k.to_string()),
                "verified": rho_ok, "group_ops": rho_ops, "fp_muls": rho_muls,
                "fp_muls_per_op": rho_muls as f64 / rho_ops.max(1) as f64, "ms": rho_ms}),
    );
    // Index calculus: trials on R = aP + bQ until k is determined.
    let opts = trial_options();
    let m = s.base.len();
    let mut rows: Vec<(Vec<(usize, i8)>, u128, u128)> = Vec::new();
    let mut agg = BTreeMap::<&str, u64>::new();
    let ti = Instant::now();
    let mut trials = 0usize;
    let mut k_ic = None;
    let mut la_ops_total = 0u64;
    let mut disagree = 0u64;
    let mut gen_muls = 0u64;
    // Trials run in parallel batches.  Under concurrency the F4 counter is
    // read as one difference over the batch (`f4_fp::field_ops_total`);
    // every other count comes from per-thread curve copies.
    let batch = 2 * rayon::current_num_threads();
    let (ca, cb) = (c.a, c.b);
    while trials < max_trials && k_ic.is_none() {
        let jobs: Vec<(u128, u128)> = (0..batch)
            .map(|_| (rng.gen_range(0..ell), rng.gen_range(0..ell)))
            .collect();
        let f4_0 = f4_fp::field_ops_total();
        let done: Vec<(u128, u128, Pt5, TrialReport, bool, u64, u64)> = jobs
            .par_iter()
            .map(|&(a, b)| {
                let cc = Curve5::with_params(p, ca, cb);
                let r = cc.add(&cc.mul(&pp, a), &cc.mul(&qq, b));
                let g = cc.f.muls();
                let rep = trial(&cc, &s.sym, &s.base, &s.index, &r, &opts, s.trace.as_ref());
                let other = cc.f.muls() - g; // system build, splitting, re-adding
                let mitm = mitm_decompositions(&cc, &s.base, &s.pairs, &r);
                let agree = mitm == rep.decompositions;
                (a, b, r, rep, agree, g, other)
            })
            .collect();
        let f4_batch = f4_fp::field_ops_total() - f4_0;
        *agg.entry("muls_f4").or_default() += f4_batch;
        for (a, b, r, rep, agree, g, other) in done {
            trials += 1;
            gen_muls += g;
            disagree += (!agree) as u64;
            *agg.entry("muls_system_split_readd").or_default() += other;
            *agg.entry("consistent").or_default() += rep.consistent as u64;
            *agg.entry("replay_mismatch").or_default() += rep.replay_mismatch as u64;
            *agg.entry("solutions").or_default() += rep.split.solutions as u64;
            *agg.entry("not_split").or_default() += rep.split.not_split as u64;
            *agg.entry("split_repeated").or_default() += rep.split.split_repeated as u64;
            *agg.entry("undetermined").or_default() += rep.undetermined as u64;
            for rel in &rep.relations {
                assert!(
                    verify_relation(c, &s.base, rel, &r),
                    "relation re-adds to R"
                );
                rows.push((rel.terms.clone(), a, b));
                emit(
                    &json,
                    &json!({"kind": "relation", "p": p, "seed": seed, "trial": trials,
                            "a": a.to_string(), "b": b.to_string(), "terms": rel.terms,
                            "verified": true}),
                );
            }
            if !rep.relations.is_empty() && rows.len() > m && k_ic.is_none() {
                let (k, ops) = solve_logs(&rows, m, h, ell);
                la_ops_total += ops;
                k_ic = k;
            }
        }
    }
    let ic_ms = ti.elapsed().as_secs_f64() * 1e3;
    let ic_ok = k_ic.map(|k| c.mul(&pp, k) == qq).unwrap_or(false);
    let trial_muls = agg["muls_system_split_readd"] + agg["muls_f4"];
    let setup_muls = s.setup["sym_precompute_muls"].as_u64().unwrap()
        + s.setup["base_muls"].as_u64().unwrap()
        + s.setup["trace_record_muls"].as_u64().unwrap();
    emit(
        &json,
        &json!({"kind": "ic", "p": p, "seed": seed, "k": k_ic.map(|k| k.to_string()),
                "verified": ic_ok, "trials": trials, "relations": rows.len(),
                "unknowns": m + 1, "agg": agg, "mitm_disagree": disagree,
                "target_gen_muls": gen_muls, "trial_muls": trial_muls,
                "linear_algebra_mulmods": la_ops_total,
                "setup_muls_target_independent": setup_muls,
                "online_muls": trial_muls + gen_muls + la_ops_total,
                "ms": ic_ms,
                "ratio_online_ic_over_rho": (trial_muls + gen_muls + la_ops_total) as f64 / rho_muls.max(1) as f64}),
    );
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("trials") => trials_mode(&args),
        Some("dlp") => dlp_mode(&args),
        _ => eprintln!("usage: joux_vitse_ladder (trials|dlp) [--p P] [--seed S] ..."),
    }
}
