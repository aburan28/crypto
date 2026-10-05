//! **Pilot timings for the presentation-vs-curve experiment.**
//!
//! One member, one `(m, l)` setting, a handful of probes: how long a
//! decomposition call takes, how much Gröbner work it does, the first fall
//! degree, and the yield the full-group instrument would need to be measured
//! at.  The design is in
//! `research/notes/koblitz-isogeny/presentation-vs-curve-design-20261003.md`.
//!
//! Run: `cargo run --release --example koblitz_ic_pilot <n> <a2> <a6> <m> <l> <probes> [full|exact] [basis_seed|g<seed>]`
//! (`KOBLITZ_PILOT_MAX_TRIALS=k` adds a measured full solve)
//! Prints one JSON line.

use crypto_lib::cryptanalysis::koblitz_isogeny_cost::*;
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let p = |i: usize| -> u64 { a[i].parse().expect("numeric argument") };
    let (n, a2, a6, m, l, probes) = (
        p(0) as u32,
        p(1) as u8,
        p(2),
        p(3) as usize,
        p(4) as u32,
        p(5) as usize,
    );
    let full = a.get(6).map(|s| s == "full").unwrap_or(true);
    // `basis_seed` = `s` for a random V, `g<s>` for a geometric V
    let geometric = a.get(7).is_some_and(|s| s.starts_with('g'));
    let basis_seed: Option<u64> = a
        .get(7)
        .and_then(|s| s.trim_start_matches('g').parse().ok());
    // `KOBLITZ_PILOT_MAX_TRIALS=k` also runs the closure phase (a measured
    // full solve checked against BSGS) instead of projecting it.
    let max_trials: usize = std::env::var("KOBLITZ_PILOT_MAX_TRIALS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);
    let irr = field_for(n).expect("field");
    let known_order = koblitz_family_order(n, a2) as u64;
    let opts = IcCostOptions {
        l,
        m,
        yield_probes: probes,
        max_trials,
        ffd_targets: 4,
        ffd_d_max: 6,
        known_order: Some(known_order),
        full_group_probes: full,
        basis_seed,
        geometric_basis: geometric,
        ..Default::default()
    };
    let t = Instant::now();
    let row = measure_member(n, &irr, a2, a6, &opts).expect("member");
    let secs = t.elapsed().as_secs_f64();
    // Heuristic yield for the full-group instrument: each unordered m-set of
    // factor-base abscissae reaches 2^m signed sums, so about
    // (2|F|)^m / (m!·#E) of the group per probe.
    let fact: f64 = (1..=m).map(|k| k as f64).product();
    let pred_yield =
        (2.0 * row.factor_base_points as f64).powi(m as i32) / (fact * known_order as f64);
    let ns_call = row.groebner_ns as f64 / row.groebner_calls.max(1) as f64;
    println!(
        "{}",
        serde_json::json!({
            "n": n, "a2": a2, "a6": a6, "m": m, "l": l, "mode": if full {"full"} else {"exact"},
            "basis_seed": basis_seed, "geometric": geometric, "probes": row.trials, "seconds": secs,
            "factor_base_points": row.factor_base_points, "unknowns": row.unknowns,
            "us_per_call": ns_call / 1e3,
            "reductions_per_call": row.reductions as f64 / row.groebner_calls.max(1) as f64,
            "relations": row.relations, "unknown_probes": row.unknown_probes,
            "first_fall_hist": row.first_fall_hist, "d_star_hist": row.d_star_hist,
            "predicted_yield_per_probe": pred_yield,
            "probes_for_100_relations": 100.0 / pred_yield.max(1e-300),
            "core_hours_for_100_relations": 100.0 / pred_yield.max(1e-300) * ns_call / 1e9 / 3600.0,
            "max_trials": max_trials, "closure_trials": row.closure_trials,
            "closure_relations": row.closure_relations, "solved_at_probe": row.solved_at_probe,
            "log": row.log, "bsgs_log": row.bsgs_log, "planted": row.planted, "verified": row.verified,
            "groebner_seconds": row.groebner_ns as f64 / 1e9,
            "linear_algebra_seconds": row.linear_algebra_ns as f64 / 1e9,
        })
    );
}
