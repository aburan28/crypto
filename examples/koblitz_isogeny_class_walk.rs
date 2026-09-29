//! **Walk the whole isogeny class by explicit isogenies, carrying one ECDLP
//! instance to every vertex.**
//!
//! For each class: build `E[ℓ]` in `F_{q^m}` (`m = ord_ℓ λ`), enumerate the
//! `ℓ + 1` kernels at every vertex where Frobenius is scalar, apply Vélu, and
//! walk descending edges from the Koblitz curve.  The instance carried is the
//! one the cost sweep solved — `d = class_secret(DEFAULT_SEED, n, r)`.
//!
//! The reached set is then compared with the **exhaustive trace-scan census**
//! of the same class, which shares no code with the walk: one enumerates
//! curves by counting points, the other by computing isogenies.
//!
//! Run: `cargo run --release --example koblitz_isogeny_class_walk [n:a2 ...]`
//! Snapshot → `experiments/koblitz_isogeny_class_walk.json`

use std::collections::BTreeSet;
use std::time::Instant;

use crypto_lib::cryptanalysis::binary_torsion_walk::*;
use crypto_lib::cryptanalysis::isogeny_class_search::{
    class_number_of_conductor, koblitz_isogeny_class,
};
use num_bigint::BigInt;

/// The volcano's prediction: `π` is scalar on `E[ℓ]` exactly at vertices whose
/// endomorphism conductor `f` has `ℓ | c/f`, so the rank-2 count is
/// `Σ h(O_f)` over those `f`.
fn predicted_rank2(factors: &[(BigInt, u32)], ell: u64) -> BigInt {
    // enumerate divisors f of c by exponent vectors
    let mut total = BigInt::from(0);
    let mut idx = vec![0u32; factors.len()];
    loop {
        let f_factors: Vec<(BigInt, u32)> = factors
            .iter()
            .zip(&idx)
            .filter(|(_, &e)| e > 0)
            .map(|((p, _), &e)| (p.clone(), e))
            .collect();
        let ell_b = BigInt::from(ell);
        let v_c = factors
            .iter()
            .find(|(p, _)| *p == ell_b)
            .map(|(_, e)| *e)
            .unwrap_or(0);
        let v_f = f_factors
            .iter()
            .find(|(p, _)| *p == ell_b)
            .map(|(_, e)| *e)
            .unwrap_or(0);
        if v_c > v_f {
            total += class_number_of_conductor(&f_factors);
        }
        // next exponent vector
        let mut i = 0;
        loop {
            if i == idx.len() {
                return total;
            }
            if idx[i] < factors[i].1 {
                idx[i] += 1;
                break;
            }
            idx[i] = 0;
            i += 1;
        }
    }
}
use crypto_lib::cryptanalysis::koblitz_isogeny_cost::{
    enumerate_class_exact, field_for, DEFAULT_SEED,
};

/// `(n, a₂)`: the prime-degree classes first, then the two-prime class that
/// exercises multi-level walking, then the two the cost sweep measured.
const DEFAULT_CASES: &[(u32, u8)] = &[(11, 1), (16, 0), (16, 1), (19, 1), (17, 1), (23, 1)];

/// Scan budget: the exhaustive trace scan is `O(4^n)`; past `2^40` it is
/// not run and the census falls back to the CM class number.
const SCAN_LOG2_BUDGET: u32 = 40;

/// The census the walk is checked against, in decreasing strength:
///
/// 1. the cost sweep's snapshot, when it covers `(n, a₂)` — an exhaustive
///    trace scan already done and committed;
/// 2. a fresh exhaustive scan, when `4^n ≤ 2^40`;
/// 3. otherwise only the **size**, `Σ_{f|c} h(O_f)` from the CM side.  At
///    that point the walk is the only thing in the repository that can name
///    the members at all, and the check is on the count, not the set.
fn census(n: u32, a2: u8) -> (Option<BTreeSet<u64>>, usize, &'static str) {
    if let Ok(text) = std::fs::read_to_string("experiments/koblitz_isogeny_cost_sweep.json") {
        if let Ok(v) = serde_json::from_str::<serde_json::Value>(&text) {
            for s in v["sweeps"].as_array().into_iter().flatten() {
                if s["n"].as_u64() == Some(n as u64) && s["a2"].as_u64() == Some(a2 as u64) {
                    let set: BTreeSet<u64> = s["rows"]
                        .as_array()
                        .into_iter()
                        .flatten()
                        .filter_map(|r| r["a6"].as_u64())
                        .collect();
                    if !set.is_empty() {
                        let len = set.len();
                        return (
                            Some(set),
                            len,
                            "cost-sweep snapshot (exhaustive trace scan)",
                        );
                    }
                }
            }
        }
    }
    if 2 * n <= SCAN_LOG2_BUDGET {
        let irr = field_for(n).expect("field");
        let c = enumerate_class_exact(n, &irr, a2);
        let set: BTreeSet<u64> = c.members.into_iter().collect();
        let len = set.len();
        return (Some(set), len, "fresh exhaustive trace scan");
    }
    let size: usize = koblitz_isogeny_class(n, 5_000_000)
        .class_size
        .to_string()
        .parse()
        .expect("class size fits");
    (
        None,
        size,
        "CM class number only — the 4^n scan is past budget",
    )
}

fn main() {
    let args: Vec<(u32, u8)> = std::env::args()
        .skip(1)
        .filter_map(|s| {
            let (a, b) = s.split_once(':')?;
            Some((a.parse().ok()?, b.parse().ok()?))
        })
        .collect();
    let cases: Vec<(u32, u8)> = if args.is_empty() {
        DEFAULT_CASES.to_vec()
    } else {
        args
    };

    println!("════════════════════════════════════════════════════════════════");
    println!("Walking isogeny classes by explicit isogenies, instance in hand");
    println!("════════════════════════════════════════════════════════════════");
    println!("E[ℓ] is built in F_(q^m), m = ord_ℓ(λ), where π acts as the scalar λ;");
    println!("every kernel, every t and every transported X must land in F_q.");

    let mut json_cases = Vec::new();
    for (n, a2) in cases {
        let irr = field_for(n).expect("field");
        let cls = koblitz_isogeny_class(n, 5_000_000);
        let ells: Vec<u64> = cls
            .nontrivial_isogeny_degrees()
            .iter()
            .filter_map(|d| d.to_string().parse().ok())
            .collect();
        let Some(inst) = koblitz_instance(n, &irr, a2, DEFAULT_SEED) else {
            println!("\n── n={n} a₂={a2}: no instance; skipped");
            continue;
        };

        println!("\n════════════════════════════════════════════════════════════════");
        println!(
            "── n = {n}, a₂ = {a2}: class of {} over conductor {} (ℓ ∈ {:?}) ──",
            cls.class_size, cls.conductor, ells
        );
        println!(
            "   instance on K_a: r = {}, d = {} (the cost sweep's class secret)",
            inst.3, inst.2
        );

        let t0 = Instant::now();
        let mut last = 0usize;
        let report = walk_class(
            n,
            &irr,
            a2,
            &ells,
            inst,
            4,
            usize::MAX,
            DEFAULT_SEED ^ ((n as u64) << 8) ^ a2 as u64,
            std::env::var_os("KOBLITZ_WALK_CHECKPOINT")
                .map(|d| std::path::PathBuf::from(d).join(format!("walk_{n}_{a2}.ckpt")))
                .as_deref(),
            |r| {
                if r.reached.len() >= last + 64 {
                    last = r.reached.len();
                    println!(
                        "     … {} vertices, {} edges, {:.1}s",
                        r.reached.len(),
                        r.edges.len(),
                        r.seconds
                    );
                }
            },
        );
        let walk_s = t0.elapsed().as_secs_f64();

        println!("   splitting fields: {:?}  (ℓ ↦ m)", report.ext_degree);
        for &ell in &ells {
            let pred = predicted_rank2(&cls.conductor_factors, ell);
            let got = report.rank2.get(&ell).copied().unwrap_or(0);
            println!(
                "   ℓ = {ell:>3}: E[ℓ] proved rank 2 at {got} vertices (volcano predicts {pred}) {},  rank 1 at {},  undecided {}",
                if pred == BigInt::from(got) { "✓" } else { "✗" },
                report.rank1.get(&ell).copied().unwrap_or(0),
                report.undecided.get(&ell).copied().unwrap_or(0)
            );
        }
        let depths = report
            .reached
            .values()
            .fold(std::collections::BTreeMap::new(), |mut m, d| {
                *m.entry(*d).or_insert(0usize) += 1;
                m
            });
        println!(
            "   reached {} vertices by {} explicit isogenies in {walk_s:.1}s; depth histogram {:?}",
            report.reached.len(),
            report.edges.len(),
            depths
        );
        println!(
            "   rationality: {} kernels with t ∉ F_q, {} images with X ∉ F_q   (both must be 0)",
            report.irrational_kernels, report.irrational_images
        );
        let ok_edges = report.edges.len() - report.transport_failures();
        println!(
            "   transport: {ok_edges} of {} edges carried (P, Q = [d]P) with #E' = #E, ord φ(P) = r, x([d]φ(P)) = x(φ(Q))",
            report.edges.len()
        );

        let (cen, cen_size, src) = census(n, a2);
        let reached: BTreeSet<u64> = report.reached.keys().copied().collect();
        let agrees = match &cen {
            Some(set) => {
                let ok = &reached == set;
                if ok {
                    println!(
                        "   census ({src}): {cen_size} members → walk REACHED EXACTLY THE CLASS"
                    );
                } else {
                    println!(
                        "   census ({src}): {cen_size} members → walk DIFFERS: {} missing (e.g. {:?}), {} extra (e.g. {:?})",
                        set.difference(&reached).count(),
                        set.difference(&reached).take(8).collect::<Vec<_>>(),
                        reached.difference(set).count(),
                        reached.difference(set).take(8).collect::<Vec<_>>()
                    );
                }
                ok
            }
            None => {
                let ok = reached.len() == cen_size;
                println!(
                    "   census ({src}): {cen_size} members → walk reached {} {}",
                    reached.len(),
                    if ok {
                        "= the class size; the walk is the only census there is at this n"
                    } else {
                        "≠ the class size"
                    }
                );
                ok
            }
        };
        let cen = cen.unwrap_or_default();
        json_cases.push(serde_json::json!({
            "n": n,
            "a2": a2,
            "ells": ells,
            "conductor": cls.conductor.to_string(),
            "class_size_cm": cls.class_size.to_string(),
            "splitting_field_degree": report.ext_degree.iter().map(|(l, m)| (l.to_string(), *m)).collect::<std::collections::BTreeMap<_, _>>(),
            "instance": {"x_p": inst.0, "x_q": inst.1, "d": inst.2, "r": inst.3},
            "reached": report.reached.len(),
            "edges": report.edges.len(),
            "edges_transported": ok_edges,
            "irrational_kernels": report.irrational_kernels,
            "irrational_images": report.irrational_images,
            "rank2_vertices": report.rank2.iter().map(|(l, c)| (l.to_string(), *c)).collect::<std::collections::BTreeMap<_, _>>(),
            "rank1_vertices": report.rank1.iter().map(|(l, c)| (l.to_string(), *c)).collect::<std::collections::BTreeMap<_, _>>(),
            "undecided_vertices": report.undecided.iter().map(|(l, c)| (l.to_string(), *c)).collect::<std::collections::BTreeMap<_, _>>(),
            "predicted_rank2": ells.iter().map(|&l| (l.to_string(), predicted_rank2(&cls.conductor_factors, l).to_string())).collect::<std::collections::BTreeMap<_, _>>(),
            "census_size": cen_size,
            "census_is_a_set": !cen.is_empty(),
            "census_source": src,
            "reached_equals_census": agrees,
            "seconds": walk_s,
            "edge_list": report.edges.iter().map(|e| serde_json::json!({
                "ell": e.ell, "from": e.from, "to": e.to,
                "x_rational": e.x_rational, "order_preserved": e.order_preserved,
                "image_order_ok": e.image_order_ok, "transported": e.transported,
            })).collect::<Vec<_>>(),
        }));
        // Write after every case so a long run leaves its finished cases on disk.
        let doc = serde_json::json!({
            "schema": "koblitz_isogeny_class_walk/v1",
            "seed": DEFAULT_SEED,
            "cases": json_cases,
        });
        let _ = std::fs::write(
            "experiments/koblitz_isogeny_class_walk.json",
            serde_json::to_string_pretty(&doc).unwrap() + "\n",
        );
    }
    println!("\nwrote experiments/koblitz_isogeny_class_walk.json");
}
