//! **The `k = 5` 2-torsion Edwards route with F4 trace replay: the trial,
//! paired against the committed test, and complete DLPs paired with rho.**
//!
//! Companion to `research/joux_vitse_n1_20261015/` (protocol
//! `PROTOCOL_EDWARDS.md`) and [`crypto_lib::cryptanalysis::joux_vitse_edwards`].
//!
//! ```text
//! # each residual through the committed test and through the replay
//! cargo run --release --example joux_vitse_edwards -- cost --p 271 --seed 1 --residuals 32 --json OUT.jsonl
//! # a complete DLP, paired with rho on the same target
//! cargo run --release --example joux_vitse_edwards -- dlp --p 31 --seed 1 --json OUT.jsonl
//! ```

use std::collections::{BTreeMap, HashMap};
use std::env;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::time::Instant;

use crypto_lib::cryptanalysis::f4_fp;
use crypto_lib::cryptanalysis::joux_vitse::solve_logs;
use crypto_lib::cryptanalysis::joux_vitse_edwards::{
    check, dlp_row, record_trace, rho, trial, trial_options, TrialE,
};
use crypto_lib::cryptanalysis::jv_quintic_edwards::{
    factor_base, generate_instance_edwards, pair_table, y2_key, Ed5, InstanceE, PairTableE, PtE,
    QuadE, SymmetrisedS5Ed,
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde_json::{json, Value};

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
    inst: InstanceE,
    base: Vec<PtE>,
    by_y2: HashMap<u64, usize>,
    pre: SymmetrisedS5Ed,
    table: PairTableE,
    trace: f4_fp::F4Trace,
    rng: StdRng,
    row: serde_json::Value,
}

/// The committed experiment's set-up, in its order (instance, base and the
/// symmetrised `S₅` from `seed ^ 0xED02`), then the trace.
fn setup(p: u64, seed: u64) -> Setup {
    let t0 = Instant::now();
    let inst = generate_instance_edwards(p, seed);
    let curve = &inst.curve;
    let f = &curve.f;
    let mut rng = StdRng::seed_from_u64(seed ^ 0xED02);
    let base = factor_base(curve, &mut rng);
    let by_y2: HashMap<u64, usize> = base
        .iter()
        .enumerate()
        .map(|(i, q)| (y2_key(curve, q), i))
        .collect();
    let fp_per_add = {
        f.reset_muls();
        let mut acc = base[0];
        for _ in 0..64 {
            acc = curve.add(&acc, &base[1]);
        }
        f.muls() as f64 / 64.0
    };
    let pre = SymmetrisedS5Ed::precompute(curve, &mut rng);
    let table = pair_table(curve, &base);
    let opts = trial_options();
    let n = inst.n;
    let mut trng = StdRng::seed_from_u64(seed ^ 0x7ACE);
    let (trace, trace_muls, trace_tries) = record_trace(curve, &pre, &opts, || {
        let a = trng.gen_range(0..n);
        curve.mul(&inst.g, a)
    });
    let row = json!({
        "kind": "setup", "p": p, "seed": seed, "n": n.to_string(),
        "log2_n": (n as f64).log2(), "d": inst.curve.d.0,
        "base_points": base.len(), "fp_muls_per_add": fp_per_add,
        "sym_precompute_muls": pre.precompute_muls,
        "trace_record_muls": trace_muls, "trace_record_tries": trace_tries,
        "trace_steps": trace.steps.len(), "ms": t0.elapsed().as_secs_f64() * 1e3,
    });
    Setup {
        inst,
        base,
        by_y2,
        pre,
        table,
        trace,
        rng,
        row,
    }
}

fn cost_mode(args: &[String]) {
    let p: u64 = arg(args, "--p", 271);
    let seed: u64 = arg(args, "--seed", 1);
    let residuals: usize = arg(args, "--residuals", 32);
    let json: Option<String> = args
        .iter()
        .position(|a| a == "--json")
        .map(|i| args[i + 1].clone());
    let mut s = setup(p, seed);
    emit(&json, &s.row);
    let curve = &s.inst.curve;
    let opts = trial_options();
    let n = s.inst.n;
    let mut rrng = StdRng::seed_from_u64(seed ^ 0x0F1E);
    for k in 0..residuals {
        let a = s.rng.gen_range(0..n);
        let b = s.rng.gen_range(1..n);
        let r = curve.add(&curve.mul(&s.inst.g, a), &curve.mul(&s.inst.q, b));
        // the committed test (no trace), then the replay, on the same r
        let (q0, t0) = trial(curve, &s.base, &s.by_y2, &s.pre, &r, &opts, None, &mut rrng);
        let (q1, t1) = trial(
            curve,
            &s.base,
            &s.by_y2,
            &s.pre,
            &r,
            &opts,
            Some(&s.trace),
            &mut rrng,
        );
        let (agree, verified) = check(curve, &s.base, &s.table, &r, &q1);
        emit(
            &json,
            &json!({"kind": "cost_trial", "p": p, "seed": seed, "trial": k,
                    "full": t0, "replay": t1, "same_decompositions": q0 == q1,
                    "oracle_agree": agree, "verified": verified,
                    "decompositions": q1.len()}),
        );
    }
}

fn dlp_mode(args: &[String]) {
    let p: u64 = arg(args, "--p", 31);
    let seed: u64 = arg(args, "--seed", 1);
    let max_trials: usize = arg(args, "--max-trials", 5_000_000);
    let json: Option<String> = args
        .iter()
        .position(|a| a == "--json")
        .map(|i| args[i + 1].clone());
    // Checkpoint after every batch (`<json>.ckpt`); `--resume` continues from
    // it.  A resumed run draws the same targets in the same order (the target
    // generator is fast-forwarded past the completed batches), so its record
    // is the one an uninterrupted run would write; only wall time differs.
    let ckpt_path = json.as_ref().map(|j| format!("{j}.ckpt"));
    let resume = args.iter().any(|a| a == "--resume");
    let ck: Option<Value> = if resume {
        ckpt_path
            .as_ref()
            .and_then(|c| fs::read_to_string(c).ok())
            .and_then(|t| serde_json::from_str(&t).ok())
    } else {
        None
    };
    if ck.is_none() {
        if let Some(j) = &json {
            fs::File::create(j).expect("truncate the record");
        }
    }
    let s = setup(p, seed);
    let curve = &s.inst.curve;
    let (n, g, q, d) = (s.inst.n, s.inst.g, s.inst.q, curve.d);
    let (rho_muls, per_op) = match &ck {
        Some(c) => (
            c["rho_fp_muls"].as_u64().unwrap(),
            c["rho_fp_muls_per_op"].as_f64().unwrap(),
        ),
        None => {
            emit(&json, &s.row);
            emit(
                &json,
                &json!({"kind": "instance", "p": p, "seed": seed, "n": n.to_string(),
                        "k_secret": s.inst.d.to_string()}),
            );
            // Rho on the same target.
            let rc = Ed5::with_params(p, d);
            let mut rrng = StdRng::seed_from_u64(seed ^ 0x5EED);
            let m0 = rc.f.muls();
            let (k_rho, rho_ops) = rho(&rc, &g, &q, n, &mut rrng);
            let rho_muls = rc.f.muls() - m0;
            let rho_ok = k_rho.map(|k| rc.mul(&g, k) == q).unwrap_or(false);
            let per_op = rho_muls as f64 / rho_ops.max(1) as f64;
            emit(
                &json,
                &json!({"kind": "rho", "p": p, "seed": seed, "k": k_rho.map(|k| k.to_string()),
                        "verified": rho_ok, "group_ops": rho_ops, "fp_muls": rho_muls,
                        "fp_muls_per_op": per_op,
                        "expected_walk_ops": (std::f64::consts::PI * n as f64 / 2.0).sqrt()}),
            );
            (rho_muls, per_op)
        }
    };
    // Index calculus, trials in parallel batches; F4 counted over the batch.
    let opts = trial_options();
    let m = s.base.len();
    let mut tgen = StdRng::seed_from_u64(seed ^ 0xD1D1);
    let mut rows: Vec<(Vec<(usize, i8)>, u128, u128)> = Vec::new();
    let mut agg = BTreeMap::<String, u64>::new();
    let (mut trials, mut la_ops, mut disagree, mut unverified) = (0usize, 0u64, 0u64, 0u64);
    let mut k_ic: Option<u128> = None;
    let mut batch = 4 * rayon::current_num_threads();
    let (mut batches_done, mut wall_before, mut resumes) = (0u64, 0f64, 0u64);
    if let Some(c) = &ck {
        batch = c["batch"].as_u64().unwrap() as usize;
        batches_done = c["batches_done"].as_u64().unwrap();
        trials = c["trials"].as_u64().unwrap() as usize;
        la_ops = c["la_ops"].as_u64().unwrap();
        disagree = c["disagree"].as_u64().unwrap();
        unverified = c["unverified"].as_u64().unwrap();
        wall_before = c["wall_ms"].as_f64().unwrap();
        resumes = c["resumes"].as_u64().unwrap() + 1;
        agg = serde_json::from_value(c["agg"].clone()).unwrap();
        for r in c["rows"].as_array().unwrap() {
            let terms: Vec<(usize, i8)> = serde_json::from_value(r["terms"].clone()).unwrap();
            let a: u128 = r["a"].as_str().unwrap().parse().unwrap();
            let b: u128 = r["b"].as_str().unwrap().parse().unwrap();
            rows.push((terms, a, b));
        }
        // the draws of the completed batches, in the loop's order
        for _ in 0..batches_done * batch as u64 {
            let _: (u64, u64, u64) = (tgen.gen_range(0..n), tgen.gen_range(0..n), tgen.gen());
        }
        emit(
            &json,
            &json!({"kind": "resume", "p": p, "seed": seed, "batches_done": batches_done,
                    "trials": trials, "relations": rows.len()}),
        );
    }
    let bump = |agg: &mut BTreeMap<String, u64>, k: &str, v: u64| {
        *agg.entry(k.to_string()).or_default() += v;
    };
    let t0 = Instant::now();
    while trials < max_trials && k_ic.is_none() {
        let jobs: Vec<(u64, u64, u64)> = (0..batch)
            .map(|_| (tgen.gen_range(0..n), tgen.gen_range(0..n), tgen.gen()))
            .collect();
        let f4_0 = f4_fp::field_ops_total();
        type Done = (u64, u64, Vec<QuadE>, TrialE, bool, bool, u64);
        let done: Vec<Done> = jobs
            .par_iter()
            .map(|&(a, b, rs)| {
                let c = Ed5::with_params(p, d);
                let r = c.add(&c.mul(&g, a), &c.mul(&q, b));
                let gen = c.f.muls();
                let mut rng = StdRng::seed_from_u64(rs);
                let (qs, t) = trial(
                    &c,
                    &s.base,
                    &s.by_y2,
                    &s.pre,
                    &r,
                    &opts,
                    Some(&s.trace),
                    &mut rng,
                );
                let (agree, verified) = check(&c, &s.base, &s.table, &r, &qs);
                (a, b, qs, t, agree, verified, gen)
            })
            .collect();
        bump(&mut agg, "f4_muls", f4_fp::field_ops_total() - f4_0);
        for (a, b, qs, t, agree, verified, gen) in done {
            trials += 1;
            bump(&mut agg, "target_gen_muls", gen);
            bump(&mut agg, "weil_muls", t.weil_muls);
            bump(&mut agg, "roots_muls", t.roots_muls);
            bump(&mut agg, "resolve_group_ops", t.group_ops);
            bump(&mut agg, "consistent", t.consistent as u64);
            bump(&mut agg, "replay_mismatch", t.replay_mismatch as u64);
            bump(&mut agg, "e_solutions", t.e_solutions as u64);
            bump(&mut agg, "not_split", t.not_split as u64);
            bump(&mut agg, "undetermined", t.undetermined as u64);
            disagree += (!agree) as u64;
            unverified += (!verified) as u64;
            for qd in &qs {
                rows.push(dlp_row(qd, a, b));
                emit(
                    &json,
                    &json!({"kind": "relation", "p": p, "seed": seed, "trial": trials,
                            "a": a.to_string(), "b": b.to_string(), "quad": qd,
                            "verified": verified}),
                );
            }
            if !qs.is_empty() && rows.len() > m && k_ic.is_none() {
                let (k, ops) = solve_logs(&rows, m, 4, n as u128);
                la_ops += ops;
                k_ic = k;
            }
        }
        batches_done += 1;
        if let Some(cp) = &ckpt_path {
            let ck = json!({
                "p": p, "seed": seed, "batch": batch, "batches_done": batches_done,
                "trials": trials, "la_ops": la_ops, "disagree": disagree,
                "unverified": unverified, "agg": agg, "resumes": resumes,
                "wall_ms": wall_before + t0.elapsed().as_secs_f64() * 1e3,
                "rho_fp_muls": rho_muls, "rho_fp_muls_per_op": per_op,
                "rows": rows.iter().map(|(t, a, b)| json!({"terms": t, "a": a.to_string(), "b": b.to_string()})).collect::<Vec<_>>(),
            });
            let tmp = format!("{cp}.tmp");
            fs::write(&tmp, ck.to_string()).expect("write checkpoint");
            fs::rename(&tmp, cp).expect("commit checkpoint");
        }
    }
    let ic_ms = wall_before + t0.elapsed().as_secs_f64() * 1e3;
    let k_ic = k_ic.map(|k| k as u64);
    let ic_ok = k_ic.map(|k| curve.mul(&g, k) == q).unwrap_or(false);
    // the resolve's group operations at the measured cost of one
    let resolve_muls = (agg["resolve_group_ops"] as f64 * per_op) as u64;
    let online = agg["f4_muls"]
        + agg["weil_muls"]
        + agg["roots_muls"]
        + agg["target_gen_muls"]
        + resolve_muls
        + la_ops;
    let expected_rho = (std::f64::consts::PI * n as f64 / 2.0).sqrt() * per_op;
    emit(
        &json,
        &json!({"kind": "ic", "p": p, "seed": seed, "k": k_ic.map(|k| k.to_string()),
                "verified": ic_ok, "trials": trials, "relations": rows.len(),
                "unknowns": m + 1, "agg": agg, "resolve_muls": resolve_muls,
                "linear_algebra_mulmods": la_ops, "oracle_disagree": disagree,
                "relations_unverified": unverified, "online_muls": online,
                "ms": ic_ms, "resumes": resumes,
                "ratio_online_ic_over_rho": online as f64 / rho_muls.max(1) as f64,
                "ratio_online_ic_over_expected_rho_walk": online as f64 / expected_rho}),
    );
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("cost") => cost_mode(&args),
        Some("dlp") => dlp_mode(&args),
        _ => eprintln!("usage: joux_vitse_edwards (cost|dlp) [--p P] [--seed S] ..."),
    }
}
