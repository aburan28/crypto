//! Prime-field ECDLP speed-up bench: the native-width engine in
//! `cryptanalysis::prime_fast` against the crate's `BigUint` rho and
//! index-calculus paths, on the same `a = −3` bench curves and the same
//! seeded known-answer targets.
//!
//! ```bash
//! cargo run --release --example prime_ecdlp_fast_bench -- [--quick] \
//!     [--baseline-ic-bits 24] [--baseline-s4-bits 16] [--rho-max-bits 56] \
//!     [--ic-max-bits 32] [--seeds 3] [--json out.json] [--sections micro,rho,ic2,ic3]
//! ```
//!
//! Sections:
//! 1. field / point micro-benchmarks (ns per op),
//! 2. rho: `pollard_rho_ecdlp` (BigUint, Floyd) vs `rho_parallel` on the
//!    16–28-bit ladder, then `rho_parallel` alone on generated 32–56-bit
//!    curves (Hasse-interval BSGS point counting),
//! 3. index calculus, 2-decomposition: `ec_index_calculus_dlp_staged`
//!    (Semaev S₃ sweep) vs `ic_solve(summands = 2)` (direct subtraction),
//! 4. index calculus, 3-decomposition: `find_one_relation_s4_counted`
//!    per-relation cost (S₄ quartic roots) vs `ic_solve(summands = 3)`
//!    (pair table), end to end.
//!
//! Public synthetic known-answer only; every recovered log is checked
//! against the planted scalar.  Rho remains the faster solver at every
//! size — the point is the per-operation cost, not a vs_rho claim.

use std::time::Instant;

use crypto_lib::cryptanalysis::ec_index_calculus::{
    build_factor_base, ec_index_calculus_dlp_staged, find_one_relation_s4_counted,
    pollard_rho_ecdlp,
};
use crypto_lib::cryptanalysis::prime_fast::{
    find_a3_curve, ic_solve, rho_parallel, FastCurve, Fp64, Pt, RhoConfig,
};
use crypto_lib::cryptanalysis::research_bench::bench_curves_a_minus_3;
use crypto_lib::ecc::curve::CurveParams;
use crypto_lib::ecc::field::FieldElement;
use crypto_lib::ecc::point::Point;
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;

struct Args {
    quick: bool,
    baseline_ic_bits: u32,
    baseline_s4_bits: u32,
    rho_max_bits: u32,
    ic_max_bits: u32,
    seeds: u64,
    json: Option<String>,
    /// Comma-separated subset of `micro,rho,ic2,ic3` (default: all).
    sections: Vec<String>,
}

fn parse_args() -> Args {
    let argv: Vec<String> = std::env::args().collect();
    let mut a = Args {
        quick: false,
        baseline_ic_bits: 20,
        baseline_s4_bits: 16,
        rho_max_bits: 56,
        ic_max_bits: 32,
        seeds: 3,
        json: None,
        sections: vec!["micro".into(), "rho".into(), "ic2".into(), "ic3".into()],
    };
    let mut i = 1;
    while i < argv.len() {
        let next = |i: &mut usize| -> String {
            *i += 1;
            argv.get(*i).cloned().unwrap_or_default()
        };
        match argv[i].as_str() {
            "--quick" => a.quick = true,
            "--baseline-ic-bits" => a.baseline_ic_bits = next(&mut i).parse().unwrap_or(20),
            "--baseline-s4-bits" => a.baseline_s4_bits = next(&mut i).parse().unwrap_or(16),
            "--rho-max-bits" => a.rho_max_bits = next(&mut i).parse().unwrap_or(56),
            "--ic-max-bits" => a.ic_max_bits = next(&mut i).parse().unwrap_or(32),
            "--seeds" => a.seeds = next(&mut i).parse().unwrap_or(3),
            "--json" => a.json = Some(next(&mut i)),
            "--sections" => a.sections = next(&mut i).split(',').map(|s| s.trim().to_string()).collect(),
            other => eprintln!("ignoring unknown argument {other}"),
        }
        i += 1;
    }
    if a.quick {
        a.baseline_ic_bits = a.baseline_ic_bits.min(16);
        a.rho_max_bits = a.rho_max_bits.min(40);
        a.ic_max_bits = a.ic_max_bits.min(24);
        a.seeds = 1;
    }
    a
}

/// The committed a = −3 P-256-class ladder rung at `bits`, or — when the
/// ladder has no rung there — the deterministic curve `find_a3_curve`
/// produces for the same policy, converted to textbook parameters so the
/// BigUint baselines can run on it too.
fn ladder_curve(bits: u32) -> Option<CurveParams> {
    if let Some((_, c)) = bench_curves_a_minus_3()
        .into_iter()
        .find(|(b, c)| *b == bits && c.name.contains("p256class"))
    {
        return Some(c);
    }
    let (fc, g) = find_a3_curve(bits, 3)?;
    let (gx, gy) = fc.canonical(g);
    Some(CurveParams {
        name: "a3-generated-p256class",
        p: BigUint::from(fc.f.p),
        a: BigUint::from(fc.f.p - 3),
        b: BigUint::from(fc.f.from_mont(fc.b)),
        gx: BigUint::from(gx),
        gy: BigUint::from(gy),
        n: BigUint::from(fc.n),
        h: 1,
    })
}

fn secs<T>(f: impl FnOnce() -> T) -> (T, f64) {
    let t = Instant::now();
    let out = f();
    (out, t.elapsed().as_secs_f64())
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if v.is_empty() {
        0.0
    } else {
        v[v.len() / 2]
    }
}

// ── 1. micro-benchmarks ───────────────────────────────────────────────

fn micro(out: &mut Vec<serde_json::Value>) {
    println!("## 1. Field and point micro-benchmarks (ns per operation)\n");
    println!("| p bits | op | BigUint `FieldElement`/`Point` | `Fp64`/`FastCurve` | speed-up |");
    println!("|---|---|---:|---:|---:|");
    let primes: Vec<(u32, u64)> = vec![(28, 268_435_399), (61, (1u64 << 61) - 1)];
    for (bits, p) in primes {
        let f = Fp64::new(p);
        let pb = BigUint::from(p);
        let mut rng = StdRng::seed_from_u64(1);
        let xs: Vec<u64> = (0..1024).map(|_| rng.gen_range(1..p)).collect();
        let xm: Vec<u64> = xs.iter().map(|&x| f.to_mont(x)).collect();
        let xb: Vec<FieldElement> = xs
            .iter()
            .map(|&x| FieldElement::new(BigUint::from(x), pb.clone()))
            .collect();

        // mul
        let iters = 200_000usize;
        let (_, tb) = secs(|| {
            let mut acc = xb[0].clone();
            for i in 0..iters {
                acc = acc.mul(&xb[i & 1023]);
            }
            std::hint::black_box(acc)
        });
        let (_, tf) = secs(|| {
            let mut acc = xm[0];
            for i in 0..iters {
                acc = f.mul(acc, xm[i & 1023]);
            }
            std::hint::black_box(acc)
        });
        let (nb, nf) = (tb / iters as f64 * 1e9, tf / iters as f64 * 1e9);
        println!("| {bits} | mul | {nb:.1} | {nf:.2} | {:.0}× |", nb / nf);
        out.push(json!({"section":"micro","bits":bits,"op":"mul","biguint_ns":nb,"fast_ns":nf}));

        // inv
        let iters = 20_000usize;
        let (_, tb) = secs(|| {
            let mut acc = BigUint::from(0u32);
            for i in 0..iters {
                acc += &xb[i & 1023].inv().unwrap().value;
            }
            std::hint::black_box(acc)
        });
        let (_, tf) = secs(|| {
            let mut acc = 0u64;
            for i in 0..iters {
                acc ^= f.inv(xm[i & 1023]);
            }
            std::hint::black_box(acc)
        });
        let (nb, nf) = (tb / iters as f64 * 1e9, tf / iters as f64 * 1e9);
        println!("| {bits} | inv | {nb:.0} | {nf:.1} | {:.0}× |", nb / nf);
        out.push(json!({"section":"micro","bits":bits,"op":"inv","biguint_ns":nb,"fast_ns":nf}));

        // batched inversion, W = 1024 (per element)
        let (_, tf) = secs(|| {
            let mut v = xm.clone();
            let mut scratch = Vec::new();
            for _ in 0..200 {
                f.batch_inv(&mut v, &mut scratch);
            }
            std::hint::black_box(v)
        });
        let nf_b = tf / (200.0 * 1024.0) * 1e9;
        println!("| {bits} | inv, batched W=1024 | — | {nf_b:.2} | {:.0}× vs BigUint inv |", nb / nf_b);
        out.push(json!({"section":"micro","bits":bits,"op":"inv_batched_1024","fast_ns":nf_b}));
    }

    // affine add on the 28-bit ladder curve and a 56-bit generated curve
    for (label, curve, g) in [
        {
            let c = ladder_curve(28).unwrap();
            let fc = FastCurve::from_params(&c).unwrap();
            let g = fc.from_textbook(&c.generator()).unwrap();
            ("28-bit ladder", fc, g)
        },
        {
            let (fc, g) = find_a3_curve(56, 3).unwrap();
            ("56-bit generated", fc, g)
        },
    ] {
        let c = CurveParams {
            name: "tmp",
            p: BigUint::from(curve.f.p),
            a: BigUint::from(curve.f.from_mont(curve.a)),
            b: BigUint::from(curve.f.from_mont(curve.b)),
            gx: BigUint::from(curve.canonical(g).0),
            gy: BigUint::from(curve.canonical(g).1),
            n: BigUint::from(curve.n),
            h: 1,
        };
        let a_fe = c.a_fe();
        let tg = c.generator();
        let pts_fast: Vec<Pt> = (1..=1024u64).map(|k| curve.mul(Some(g), k).unwrap()).collect();
        let pts_text: Vec<Point> = pts_fast.iter().map(|p| curve.to_textbook(Some(*p))).collect();
        let iters = 20_000usize;
        let (_, tb) = secs(|| {
            let mut acc = tg.clone();
            for i in 0..iters {
                acc = acc.add(&pts_text[i & 1023], &a_fe);
            }
            std::hint::black_box(acc)
        });
        let (_, tf) = secs(|| {
            let mut acc = Some(g);
            for i in 0..iters {
                acc = curve.add(acc, Some(pts_fast[i & 1023]));
            }
            std::hint::black_box(acc)
        });
        // batched: 1024 independent walkers each adding a table point
        let f = curve.f;
        let (_, tbat) = secs(|| {
            let mut xs: Vec<u64> = pts_fast.iter().map(|p| p.x).collect();
            let mut ys: Vec<u64> = pts_fast.iter().map(|p| p.y).collect();
            let mut d = vec![0u64; 1024];
            let mut scratch = Vec::new();
            for round in 0..200 {
                let t = pts_fast[(round * 7 + 3) & 1023];
                for i in 0..1024 {
                    let dd = f.sub(t.x, xs[i]);
                    d[i] = if dd == 0 { f.one } else { dd };
                }
                f.batch_inv(&mut d, &mut scratch);
                for i in 0..1024 {
                    let s = curve.add_with_inv(Pt { x: xs[i], y: ys[i] }, t, d[i]);
                    xs[i] = s.x;
                    ys[i] = s.y;
                }
            }
            std::hint::black_box((xs, ys))
        });
        let nb = tb / iters as f64 * 1e9;
        let nf = tf / iters as f64 * 1e9;
        let nbat = tbat / (200.0 * 1024.0) * 1e9;
        println!("| {label} | affine add (own inverse) | {nb:.0} | {nf:.0} | {:.0}× |", nb / nf);
        println!("| {label} | affine add, batched W=1024 | — | {nbat:.1} | {:.0}× vs BigUint add |", nb / nbat);
        out.push(json!({"section":"micro","curve":label,"op":"affine_add","biguint_ns":nb,"fast_ns":nf,"fast_batched_ns":nbat}));
    }
    println!();
}

// ── 2. rho ────────────────────────────────────────────────────────────

fn rho_section(args: &Args, out: &mut Vec<serde_json::Value>) {
    println!("## 2. Pollard rho: BigUint Floyd (`pollard_rho_ecdlp`) vs native parallel DP (`rho_parallel`)\n");
    println!("Known-answer targets `Q = [k]G`, `k` from seed; medians over {} seed(s).  Baseline uses the ladder examples' step cap `2^(bits/2+8)`.\n", args.seeds);
    println!("| bits | curve | baseline wall s | fast wall s | speed-up | fast steps/s (total) | threads × walkers | dp bits | fast steps (median) |");
    println!("|---|---|---:|---:|---:|---:|---|---|---:|");
    let mut bits_list: Vec<u32> = vec![16, 20, 24, 28];
    let mut b = 32;
    while b <= args.rho_max_bits {
        bits_list.push(b);
        b += 4;
    }
    for bits in bits_list {
        let (fc, g, textbook): (FastCurve, Pt, Option<CurveParams>) = if bits <= 28 {
            let c = ladder_curve(bits).unwrap();
            let fc = FastCurve::from_params(&c).unwrap();
            let g = fc.from_textbook(&c.generator()).unwrap();
            (fc, g, Some(c))
        } else {
            let (t, gen_secs) = secs(|| find_a3_curve(bits, 3));
            let Some((fc, g)) = t else {
                eprintln!("no curve at {bits} bits");
                continue;
            };
            println!(
                "<!-- generated {}: p={} b={} n={} G=({}, {}) in {:.2}s -->",
                fc.name,
                fc.f.p,
                fc.f.from_mont(fc.b),
                fc.n,
                fc.canonical(g).0,
                fc.canonical(g).1,
                gen_secs
            );
            (fc, g, None)
        };
        let mut base_walls = Vec::new();
        let mut fast_walls = Vec::new();
        let mut rates = Vec::new();
        let mut steps = Vec::new();
        let mut cfg = RhoConfig::for_bits(bits, true);
        for seed in 0..args.seeds {
            let mut rng = StdRng::seed_from_u64(20261007 + seed * 1000 + bits as u64);
            let k = rng.gen_range(1..fc.n);
            let q = fc.mul(Some(g), k).unwrap();
            if let Some(c) = &textbook {
                let a_fe = c.a_fe();
                let tg = c.generator();
                let tq = tg.scalar_mul(&BigUint::from(k), &a_fe);
                assert_eq!(fc.to_textbook(Some(q)), tq, "fast and textbook [k]G disagree");
                let max_steps = (1usize << (bits / 2 + 8)).min(5_000_000);
                let (res, wall) = secs(|| pollard_rho_ecdlp(c, &tg, &tq, max_steps));
                match res {
                    Some(r) => assert_eq!(r, BigUint::from(k), "baseline rho wrong answer"),
                    None => eprintln!("baseline rho timed out at {bits} bits seed {seed}"),
                }
                base_walls.push(wall);
            }
            cfg.seed = 7 + seed;
            let res = rho_parallel(&fc, g, q, &cfg).expect("fast rho solves");
            assert_eq!(res.log, k, "fast rho wrong answer at {bits} bits");
            fast_walls.push(res.wall_secs);
            rates.push(res.steps_per_sec);
            steps.push(res.total_steps as f64);
        }
        let bw = median(&mut base_walls);
        let fw = median(&mut fast_walls);
        let rate = median(&mut rates);
        let st = median(&mut steps);
        let speedup = if textbook.is_some() { format!("{:.0}×", bw / fw) } else { "—".into() };
        let bws = if textbook.is_some() { format!("{bw:.3}") } else { "—".into() };
        println!(
            "| {bits} | {} | {bws} | {fw:.4} | {speedup} | {:.2e} | {} × {} | {} | {:.3e} |",
            fc.name, rate, cfg.threads, cfg.walkers, cfg.dp_bits, st
        );
        out.push(json!({"section":"rho","bits":bits,"curve":fc.name,"p":fc.f.p,"n":fc.n,
            "baseline_wall_s": if textbook.is_some() { Some(bw) } else { None },
            "fast_wall_s":fw,"fast_steps_per_s":rate,"fast_steps_median":st,
            "threads":cfg.threads,"walkers":cfg.walkers,"dp_bits":cfg.dp_bits,"seeds":args.seeds}));
    }
    println!();
}

// ── 3. index calculus, 2-decomposition ────────────────────────────────

fn ic2_section(args: &Args, out: &mut Vec<serde_json::Value>) {
    println!("## 3. Index calculus, 2-decomposition: Semaev S₃ sweep (BigUint) vs direct subtraction (native)\n");
    println!("Factor base policy `m = 120` smallest-x points, `+6` extra relations (the a=−3 ladder policy).  Both drivers verify `[x]G = Q`.\n");
    println!("| bits | baseline total s (fb / rel / LA) | fast total s (fb / rel / LA) | speed-up | fast targets | fast rows (1/2) | rho fast s |");
    println!("|---|---|---|---:|---:|---|---:|");
    let mut bits_list: Vec<u32> = vec![16, 20, 24, 28];
    let mut b = 32;
    while b <= args.ic_max_bits {
        bits_list.push(b);
        b += 4;
    }
    for bits in bits_list {
        let (fc, g, textbook) = if bits <= 28 {
            let c = ladder_curve(bits).unwrap();
            let fc = FastCurve::from_params(&c).unwrap();
            let g = fc.from_textbook(&c.generator()).unwrap();
            (fc, g, Some(c))
        } else {
            let Some((fc, g)) = find_a3_curve(bits, 3) else { continue };
            (fc, g, None)
        };
        let mut rng = StdRng::seed_from_u64(20261007 + bits as u64);
        let k = rng.gen_range(1..fc.n);
        let q = fc.mul(Some(g), k).unwrap();
        let mut base = None;
        if bits <= args.baseline_ic_bits {
            if let Some(c) = &textbook {
                let tg = c.generator();
                let tq = tg.scalar_mul(&BigUint::from(k), &c.a_fe());
                let max_trials = 1usize << (bits + 6);
                let (res, wall) = secs(|| {
                    ec_index_calculus_dlp_staged(c, &tg, &tq, 120, 6, max_trials.min(5_000_000), 8)
                });
                match res {
                    Some((x, rep)) => {
                        assert_eq!(x, BigUint::from(k));
                        base = Some((wall, rep.factor_base_ms, rep.relations_ms, rep.linear_algebra_ms, rep.trials_total));
                    }
                    None => eprintln!("baseline IC failed at {bits} bits"),
                }
            }
        }
        let (res, wall) = secs(|| ic_solve(&fc, g, q, 120, 6, 2, u64::MAX, 3));
        let Some((x, rep)) = res else {
            eprintln!("fast IC failed at {bits} bits");
            continue;
        };
        assert_eq!(x, k, "fast IC wrong answer");
        let mut cfg = RhoConfig::for_bits(bits, true);
        cfg.seed = 3;
        let rho = rho_parallel(&fc, g, q, &cfg).unwrap();
        assert_eq!(rho.log, k);
        let base_str = match base {
            Some((w, fbm, rm, lam, _)) => format!("{w:.2} ({:.3} / {:.2} / {:.3})", fbm / 1e3, rm / 1e3, lam / 1e3),
            None => "—".into(),
        };
        let speedup = match base {
            Some((w, ..)) => format!("{:.0}×", w / wall),
            None => "—".into(),
        };
        println!(
            "| {bits} | {base_str} | {wall:.3} ({:.4} / {:.3} / {:.4}) | {speedup} | {} | {}/{} | {:.4} |",
            rep.factor_base_ms / 1e3,
            rep.relations_ms / 1e3,
            rep.linear_algebra_ms / 1e3,
            rep.targets,
            rep.rows_one,
            rep.rows_two,
            rho.wall_secs
        );
        out.push(json!({"section":"ic2","bits":bits,"curve":fc.name,"fb":rep.factor_base_size,"relations":rep.relations,
            "baseline_total_s": base.map(|b| b.0), "baseline_trials": base.map(|b| b.4),
            "fast_total_s":wall,"fast_fb_ms":rep.factor_base_ms,"fast_relations_ms":rep.relations_ms,
            "fast_la_ms":rep.linear_algebra_ms,"fast_targets":rep.targets,"rows_one":rep.rows_one,"rows_two":rep.rows_two,
            "rho_fast_s":rho.wall_secs}));
    }
    println!();
}

// ── 4. index calculus, 3-decomposition ────────────────────────────────

fn ic3_section(args: &Args, out: &mut Vec<serde_json::Value>) {
    println!("## 4. Index calculus, 3-decomposition: S₄ quartic roots (BigUint) vs pair table (native)\n");
    println!("Baseline cost is measured per relation with `find_one_relation_s4_counted` (fb = 120); the fast driver runs end to end with the same base and `+6` extra relations.\n");
    println!("| bits | baseline s/relation (relations timed) | fast total s (table / rel / LA) | fast s/relation | per-relation speed-up | pair table entries | fast targets | rows (1/2/3) |");
    println!("|---|---:|---|---:|---:|---:|---:|---|");
    let mut bits_list: Vec<u32> = vec![16, 20, 24, 28];
    let mut b = 32;
    while b <= args.ic_max_bits {
        bits_list.push(b);
        b += 4;
    }
    for bits in bits_list {
        let (fc, g, textbook) = if bits <= 28 {
            let c = ladder_curve(bits).unwrap();
            let fc = FastCurve::from_params(&c).unwrap();
            let g = fc.from_textbook(&c.generator()).unwrap();
            (fc, g, Some(c))
        } else {
            let Some((fc, g)) = find_a3_curve(bits, 3) else { continue };
            (fc, g, None)
        };
        let mut rng = StdRng::seed_from_u64(20261007 + bits as u64);
        let k = rng.gen_range(1..fc.n);
        let q = fc.mul(Some(g), k).unwrap();
        let mut base_per_rel = None;
        let mut base_count = 0usize;
        if bits <= args.baseline_s4_bits {
            if let Some(c) = &textbook {
                let tg = c.generator();
                let tq = tg.scalar_mul(&BigUint::from(k), &c.a_fe());
                let fb = build_factor_base(c, 120);
                let timed = if bits <= 16 { 2 } else { 1 };
                let (_, wall) = secs(|| {
                    for _ in 0..timed {
                        let r = find_one_relation_s4_counted(c, &tg, &tq, &fb, 1 << 16);
                        assert!(r.is_some(), "baseline S4 relation search failed");
                    }
                });
                base_per_rel = Some(wall / timed as f64);
                base_count = timed;
            }
        }
        let (res, wall) = secs(|| ic_solve(&fc, g, q, 120, 6, 3, u64::MAX, 5));
        let Some((x, rep)) = res else {
            eprintln!("fast S4 IC failed at {bits} bits");
            continue;
        };
        assert_eq!(x, k, "fast 3-decomp IC wrong answer");
        let fast_per_rel = rep.relations_ms / 1e3 / rep.relations as f64;
        let base_str = match base_per_rel {
            Some(b) => format!("{b:.3} ({base_count})"),
            None => "—".into(),
        };
        let speedup = match base_per_rel {
            Some(b) => format!("{:.0}×", b / fast_per_rel),
            None => "—".into(),
        };
        println!(
            "| {bits} | {base_str} | {wall:.3} ({:.3} / {:.3} / {:.4}) | {fast_per_rel:.6} | {speedup} | {} | {} | {}/{}/{} |",
            rep.pair_table_ms / 1e3,
            rep.relations_ms / 1e3,
            rep.linear_algebra_ms / 1e3,
            rep.pair_table_entries,
            rep.targets,
            rep.rows_one,
            rep.rows_two,
            rep.rows_three
        );
        out.push(json!({"section":"ic3","bits":bits,"curve":fc.name,"fb":rep.factor_base_size,"relations":rep.relations,
            "baseline_s_per_relation":base_per_rel,"baseline_relations_timed":base_count,
            "fast_total_s":wall,"fast_pair_table_ms":rep.pair_table_ms,"fast_relations_ms":rep.relations_ms,
            "fast_la_ms":rep.linear_algebra_ms,"fast_s_per_relation":fast_per_rel,"pair_table_entries":rep.pair_table_entries,
            "fast_targets":rep.targets,"rows_one":rep.rows_one,"rows_two":rep.rows_two,"rows_three":rep.rows_three}));
    }
    println!();
}

fn main() {
    let args = parse_args();
    let mut out: Vec<serde_json::Value> = Vec::new();
    let host = json!({
        "threads": std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1),
        "arch": std::env::consts::ARCH,
        "os": std::env::consts::OS,
    });
    println!("# Prime-field ECDLP speed-up bench\n");
    println!("Host: {} {} with {} hardware threads.\n", host["os"], host["arch"], host["threads"]);
    let on = |name: &str| args.sections.iter().any(|s| s == name);
    if on("micro") {
        micro(&mut out);
    }
    if on("rho") {
        rho_section(&args, &mut out);
    }
    if on("ic2") {
        ic2_section(&args, &mut out);
    }
    if on("ic3") {
        ic3_section(&args, &mut out);
    }
    if let Some(path) = &args.json {
        let doc = json!({"host": host, "results": out});
        std::fs::write(path, serde_json::to_string_pretty(&doc).unwrap()).expect("write json");
        println!("JSON written to {path}");
    }
}
