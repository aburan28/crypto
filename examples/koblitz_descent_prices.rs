//! Ledger §22: what a Koblitz descent spends per target, part by part.
//!
//! Ledger §21 left the descent the largest phase at six of nine sizes,
//! and below `2^37` its cost is not probes: at `n = 19` a target takes
//! about seven trials and 1,461 units.  This prices the descent of every
//! target of a workflow parameter file, on one thread, in the thread's
//! unit (one batched affine addition, `add_many` over 1,024 points)
//! measured in the same process:
//!
//! - the whole `IndividualLogSolver::solve_report`, and its trials, which
//!   must equal the pricer's for the same file;
//! - the same call under a measurement session, split into walk set-up
//!   and stepping (`target_query`), keying and lookups (`target_pdp`),
//!   relation assembly (`target_descent`) and the recovery check
//!   (`recovery_check`);
//! - the parts suspected of being fixed per target, each timed alone:
//!   the three scalar multiplications and 63 additions that start the 64
//!   walks, and `[d]G = Q` in the general (big-integer) arithmetic and in
//!   the single-word one.
//!
//! The log table comes from an `ic workflow` run on the same file, so
//! nothing here depends on the collection.
//!
//!     ic workflow --params P --dir D
//!     cargo run --release --example koblitz_descent_prices -- P D/logs.json [tier] [rounds]
//!
//! `tier` is the pair table the pipeline built for P (`folded`, `compact`,
//! `witnessed_compact` or `full`), as the pricer's report names it.
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::ic_measurement::Session;
use crypto_lib::cryptanalysis::koblitz_fast::{BatchScratch, FastCurve, FastPoint};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_subgroup_orbit_factor_base, DecompositionStrategy, FactorBaseLogTable,
    IndividualLogSolver, KoblitzCurve, KoblitzIcOptions, PairSumTable,
};
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::time::Instant;

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

fn hex(s: &str) -> BigUint {
    BigUint::parse_bytes(s.trim_start_matches("0x").as_bytes(), 16).expect("hex coordinate")
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let params_path = args.get(1).expect("a workflow parameter file");
    let logs_path = args.get(2).expect("the run's logs.json");
    let tier = args.get(3).map(String::as_str).unwrap_or("folded");
    let rounds: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(5);
    let p: Value =
        serde_json::from_str(&std::fs::read_to_string(params_path).expect("params")).expect("json");
    let logs: Value =
        serde_json::from_str(&std::fs::read_to_string(logs_path).expect("logs")).expect("json");

    let n = p["curve"]["degree"].as_u64().expect("degree") as u32;
    let a = p["curve"]["curve_a"].as_u64().unwrap_or(0) as u8;
    let spec = &p["factor_base"]["spec"];
    let points = spec["points"].as_u64().expect("points") as usize;
    let base_seed = spec["seed"].as_u64().expect("seed");
    let seed = p["seed"].as_u64().expect("seed");

    let kc = KoblitzCurve::new(a, n).expect("a Koblitz curve");
    let fc = FastCurve::new(&kc.curve).expect("single-word arithmetic");
    let fb = build_subgroup_orbit_factor_base(&kc, base_seed, points).expect("a base");
    let budget = PairSumTable::DEFAULT_BYTE_BUDGET;
    let pair = match tier {
        "folded" => PairSumTable::build_folded_within(&kc, &fb, budget),
        "compact" => PairSumTable::build_compact_within(&kc, &fb, budget),
        "witnessed_compact" => PairSumTable::build_witnessed_compact_within(&kc, &fb, budget),
        "full" => PairSumTable::build_full_within(&kc, &fb, budget),
        other => panic!("unknown tier {other}"),
    }
    .expect("the pair table");
    let table = FactorBaseLogTable {
        columns: logs["columns"]
            .as_array()
            .expect("columns")
            .iter()
            .map(|c| {
                let x = F2mElement::from_biguint(&hex(c["x"].as_str().expect("x")), n);
                let y = F2mElement::from_biguint(&hex(c["y"].as_str().expect("y")), n);
                let log = c["log"]
                    .as_str()
                    .expect("log")
                    .parse::<BigUint>()
                    .expect("log");
                (BinaryPoint::Affine { x, y }, log)
            })
            .collect(),
    };
    let opts = KoblitzIcOptions {
        m: p["summands"].as_u64().unwrap_or(3) as usize,
        descent_m: p["descent_summands"].as_u64().map(|m| m as usize),
        collection_window: p["collection_window"].as_u64().map(|w| w as usize),
        strategy: DecompositionStrategy::PairTable,
        collapse_negation: true,
        collapse_projected_orbits: true,
        crossbred: None,
        allow_direct_relation: false,
        max_trials: p["max_trials"].as_u64().unwrap_or(100_000) as usize,
        seed,
        ..KoblitzIcOptions::default()
    };
    let solver =
        IndividualLogSolver::new(&kc, &fb, &table, &opts, Some(&pair)).expect("a descent solver");

    // The workflow's targets: a known scalar from each random seed.
    let r_u64 = kc.subgroup_order.to_u64_digits()[0];
    let targets: Vec<(BinaryPoint, BigUint)> = p["targets"]
        .as_array()
        .expect("targets")
        .iter()
        .map(|t| {
            let s = t["random_seed"].as_u64().expect("random_seed");
            let mut rng = StdRng::seed_from_u64(s ^ 0x534f_4c56_4552_5447);
            let k = BigUint::from(rng.gen_range(1..r_u64));
            (kc.mul(kc.generator(), &k), k)
        })
        .collect();

    // The unit, as `ic price` measures it.
    let g = fc.lift(kc.generator());
    let batch: Vec<FastPoint> = (0..1024u64)
        .map(|i| fc.mul(g, &BigUint::from(1 + i * 7919)))
        .collect();
    let mut scratch = BatchScratch::default();
    let mut out = Vec::with_capacity(1024);
    let mut unit = || {
        let t = Instant::now();
        for _ in 0..1000 {
            out.clear();
            fc.add_many(g, &batch, &mut out, &mut scratch);
        }
        t.elapsed().as_nanos() as f64 / 1_024_000.0
    };
    unit();

    // What each part costs, per target, in ns; converted at the end.
    let phases = [
        "target_query",
        "target_pdp",
        "target_descent",
        "recovery_check",
    ];
    let mut whole: Vec<Vec<f64>> = vec![Vec::new(); targets.len()];
    let mut by_phase: Vec<Vec<Vec<f64>>> = vec![vec![Vec::new(); phases.len()]; targets.len()];
    let mut starts: Vec<f64> = Vec::new();
    let (mut one_mul, mut sequential, mut batched) = (Vec::new(), Vec::new(), Vec::new());
    let mut batched_out = Vec::with_capacity(64);
    let mut batch_scratch = BatchScratch::default();
    let mut check_general: Vec<f64> = Vec::new();
    let mut check_single: Vec<f64> = Vec::new();
    let mut trials = vec![0usize; targets.len()];
    let mut units_ns = Vec::new();
    let stride = (r_u64 / 64).max(1 << 20);
    for _ in 0..rounds {
        units_ns.push(unit());
        for (i, (q, k)) in targets.iter().enumerate() {
            let t = Instant::now();
            let report = solver.solve_report(q);
            whole[i].push(t.elapsed().as_nanos() as f64);
            assert_eq!(report.log.as_ref(), Some(k), "target {i} was not recovered");
            trials[i] = report.trials;

            let session = Session::begin().expect("a measurement session");
            let again = solver.solve_report(q);
            let snap = session.finish().expect("a snapshot");
            assert_eq!(again.trials, report.trials);
            for (j, name) in phases.iter().enumerate() {
                by_phase[i][j].push(snap.phases_ns[name].unwrap_or(0) as f64);
            }

            // The walk start, as the solver builds it: three scalar
            // multiplications and 63 additions (its own random a0, b).
            let qf = fc.lift(q);
            let t = Instant::now();
            let stride_point = fc.mul_u64(g, stride);
            let mut state = fc.add(
                fc.mul_u64(g, 1 + (i as u64) % (r_u64 - 1)),
                fc.mul_u64(qf, 3),
            );
            for _ in 0..63 {
                state = fc.add(state, stride_point);
            }
            std::hint::black_box(state);
            starts.push(t.elapsed().as_nanos() as f64);

            // Its pieces: one scalar multiplication, the 63 additions one
            // at a time, and the same 63 points as one batched addition of
            // the start to precomputed multiples of the stride.
            let t = Instant::now();
            std::hint::black_box(fc.mul_u64(qf, 3 + (i as u64) % (r_u64 - 3)));
            one_mul.push(t.elapsed().as_nanos() as f64);
            let t = Instant::now();
            let mut s = state;
            for _ in 0..63 {
                s = fc.add(s, stride_point);
            }
            std::hint::black_box(s);
            sequential.push(t.elapsed().as_nanos() as f64);
            let multiples: Vec<FastPoint> = (1..64u64).map(|j| fc.mul_u64(g, j * stride)).collect();
            let t = Instant::now();
            batched_out.clear();
            fc.add_many(state, &multiples, &mut batched_out, &mut batch_scratch);
            std::hint::black_box(&batched_out);
            batched.push(t.elapsed().as_nanos() as f64);

            // [d]G = Q, general and single-word.
            let t = Instant::now();
            std::hint::black_box(kc.mul(kc.generator(), k) == *q);
            check_general.push(t.elapsed().as_nanos() as f64);
            let t = Instant::now();
            std::hint::black_box(fc.mul(g, k) == qf);
            check_single.push(t.elapsed().as_nanos() as f64);
        }
        units_ns.push(unit());
    }
    let u = median(&mut units_ns);
    let per_target = |v: &mut Vec<f64>| median(v) / u;
    let mut totals: Vec<f64> = whole.iter_mut().map(|v| median(v) / u).collect();
    let phase_units: Vec<Value> = phases
        .iter()
        .enumerate()
        .map(|(j, name)| {
            let mut v: Vec<f64> = by_phase.iter_mut().map(|t| median(&mut t[j]) / u).collect();
            let mean = v.iter().sum::<f64>() / v.len() as f64;
            json!({"phase": name, "mean_units_per_target": mean, "median_units_per_target": median(&mut v)})
        })
        .collect();
    let mean_total = totals.iter().sum::<f64>() / totals.len() as f64;
    let out = json!({
        "params": params_path, "logs": logs_path, "tier": tier, "rounds": rounds,
        "n": n, "a": a, "points": fb.points.len(), "columns": table.columns.len(),
        "unit_ns": u,
        "trials_per_target": trials,
        "trials_total": trials.iter().sum::<usize>(),
        "solve_mean_units_per_target": mean_total,
        "solve_median_units_per_target": median(&mut totals),
        "phases": phase_units,
        "walk_start_units": per_target(&mut starts),
        "one_scalar_multiplication_units": per_target(&mut one_mul),
        "sixty_three_additions_sequential_units": per_target(&mut sequential),
        "sixty_three_additions_batched_units": per_target(&mut batched),
        "recovery_check_general_units": per_target(&mut check_general),
        "recovery_check_single_word_units": per_target(&mut check_single),
    });
    println!("{}", serde_json::to_string_pretty(&out).expect("json"));
}
