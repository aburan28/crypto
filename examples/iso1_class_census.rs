//! Exact weak-form representatives and trace-level CM features for ISO-1.
//!
//! This is a diagnostic census, not an isogeny walk or a DLP benchmark.
//! `p` is prime, `q = p^2`, and the curve field is `F_{q^3}`.

use std::collections::{BTreeMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time::Instant;

use crypto_lib::cryptanalysis::jv_cover::{Fld, Fq3, E2, E6};
use crypto_lib::cryptanalysis::jv_isogeny_walk::{curve_order, random_curve, Curve2};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::prelude::*;

fn factors(mut n: u128) -> Vec<(u128, u32)> {
    let mut out = Vec::new();
    let mut d = 2;
    while d * d <= n {
        if n % d == 0 {
            let mut e = 0;
            while n % d == 0 {
                n /= d;
                e += 1;
            }
            out.push((d, e));
        }
        d = if d == 2 { 3 } else { d + 2 };
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

fn cm_features(p: u64, t: i128) -> (i128, u32, i8, u32, bool) {
    let q3 = (p as i128).pow(6);
    let d = (t / 2).pow(2) - q3;
    let fac = factors((-d) as u128);
    let mut squarefree = -1_i128;
    let mut square_part = 1_i128;
    for &(prime, exp) in &fac {
        if exp % 2 == 1 {
            squarefree *= prime as i128;
        }
        square_part *= (prime as i128).pow(exp / 2);
    }
    let fundamental = if squarefree.rem_euclid(4) == 1 {
        squarefree
    } else {
        4 * squarefree
    };
    let f_pi = if fundamental == squarefree {
        2 * square_part
    } else {
        square_part
    };
    let depth = f_pi.trailing_zeros();
    let split2 = if fundamental % 2 == 0 {
        0
    } else if fundamental.rem_euclid(8) == 1 {
        1
    } else {
        -1
    };
    let h_odd = factors((-fundamental) as u128).len() == 1;
    debug_assert_eq!(f_pi * f_pi * fundamental, 4 * d);
    (
        fundamental,
        depth,
        split2,
        (-d as u128).trailing_zeros(),
        h_odd,
    )
}

fn exact_order_by_squares(f: &Fq3, c: &Curve2) -> u128 {
    exact_order_with_squares(f, c, &field_squares(f))
}

fn field_squares(f: &Fq3) -> HashSet<E6> {
    let p = f.f.p;
    let q = p * p;
    let q3 = q * q * q;
    let e2 = |k: u64| E2([k % p, k / p]);
    let e6 = |k: u64| E6([e2(k % q), e2((k / q) % q), e2(k / (q * q))]);
    let mut squares = HashSet::new();
    for k in 0..q3 {
        squares.insert(f.sq(&e6(k)));
    }
    squares
}

fn exact_order_with_squares(f: &Fq3, c: &Curve2, squares: &HashSet<E6>) -> u128 {
    let p = f.f.p;
    let q = p * p;
    let q3 = q * q * q;
    let e2 = |k: u64| E2([k % p, k / p]);
    let e6 = |k: u64| E6([e2(k % q), e2((k / q) % q), e2(k / (q * q))]);
    let mut order = 1_u128;
    for k in 0..q3 {
        let x = e6(k);
        let rhs = f.mul(
            &f.mul(&f.sub(&x, &c.e[0]), &f.sub(&x, &c.e[1])),
            &f.sub(&x, &c.e[2]),
        );
        order += if rhs == E6::ZERO {
            1
        } else if squares.contains(&rhs) {
            2
        } else {
            0
        };
    }
    order
}

fn absolute_orbit_weight(f: &Fq3, lambda: E6, p: u64) -> Option<u64> {
    // The absolute p-Frobenius preserves the point count and the norm-one
    // torus. Together with inversion it gives a twelve-element action.
    let lambda_p = f.pow(&lambda, p as u128);
    let lambda_q = f.sigma(&lambda);
    let lambda_pq = f.sigma(&lambda_p);
    let lambda_q2 = f.sigma(&lambda_q);
    let lambda_pq2 = f.sigma(&lambda_pq);
    let lambda_inv = f.mul(&lambda_q, &lambda_q2);
    let lambda_p_inv = f.mul(&lambda_pq, &lambda_pq2);
    let lambda_inv_q = f.sigma(&lambda_inv);
    let lambda_p_inv_q = f.sigma(&lambda_p_inv);
    let mut orbit = [
        lambda,
        lambda_p,
        lambda_q,
        lambda_pq,
        lambda_q2,
        lambda_pq2,
        lambda_inv,
        lambda_p_inv,
        lambda_inv_q,
        lambda_p_inv_q,
        f.sigma(&lambda_inv_q),
        f.sigma(&lambda_p_inv_q),
    ];
    debug_assert_ne!(lambda, E6::ONE);
    debug_assert_eq!(f.mul(&lambda, &lambda_inv), E6::ONE);
    orbit.sort_unstable();
    if orbit[0] != lambda {
        return None;
    }
    let weight = 1 + orbit.windows(2).filter(|pair| pair[0] != pair[1]).count();
    debug_assert!(matches!(weight, 2 | 6 | 12));
    Some(weight as u64)
}

fn census(
    p: u64,
    output: &str,
    verify_trace: Option<i128>,
    derive_twists: bool,
    orbit_quotient: bool,
    absolute_orbit_quotient: bool,
    exact_squares: bool,
) {
    assert!(!(orbit_quotient && absolute_orbit_quotient));
    assert!(!(orbit_quotient || absolute_orbit_quotient) || derive_twists);
    let start = Instant::now();
    let f = Fq3::new(p);
    assert!(
        !exact_squares || p <= 13,
        "exact square-table census requires p <= 13"
    );
    let squares = exact_squares.then(|| field_squares(&f));
    let q = p * p;
    let q3 = (p as i128).pow(6);
    let e2 = |k: u64| E2([k % p, k / p]);
    // Every nonzero F_p element is a square in F_{p²}.  The second
    // representative must be a nonsquare in F_{p²}, not f.f.w ∈ F_p.
    let w = (1..q).map(e2).find(|x| !f.f.is_square(x)).unwrap();
    assert!(!f.f.is_square(&w));
    let branches = if derive_twists {
        vec![(false, E2::ONE), (true, E2::ONE)]
    } else {
        vec![(false, E2::ONE), (false, w), (true, E2::ONE), (true, w)]
    };
    let jobs: Vec<_> = branches
        .into_iter()
        .flat_map(|(branch, val)| (0..q).map(move |a0| (branch, val, a0)))
        .collect();
    let (mut weak, mut reps, calls, muls, witness) = jobs
        .into_par_iter()
        .map(|(branch, val, a0)| {
            let f = Fq3::new(p);
            let mut weak = BTreeMap::<i128, u64>::new();
            let mut witness = None;
            let mut reps = 0_u64;
            let mut calls = 0_u64;
            let mut rng = StdRng::seed_from_u64(
                0x1501_u64
                    ^ p.wrapping_mul(0x9e3779b1)
                    ^ a0
                    ^ ((branch as u64) << 48)
                    ^ ((val == w) as u64) << 49,
            );
            let count = if branch { 1 } else { q };
            for a2 in 0..count {
                let alpha = if branch {
                    E6([e2(a0), E2::ZERO, val])
                } else {
                    E6([e2(a0), val, e2(a2)])
                };
                let c = Curve2 {
                    e: [E6::ZERO, alpha, f.sigma(&alpha)],
                };
                debug_assert!(c.weak_by_norms(&f));
                // Hilbert 90 identifies these normalized alpha values with
                // T = ker(N_{F_(q^3)/F_q}) minus 1 via lambda = sigma(alpha)/alpha.
                // The trace is unchanged by q-Frobenius or lambda inversion;
                // twist derivation makes the latter valid regardless of the
                // square class of this particular alpha representative.
                let weight = if absolute_orbit_quotient {
                    let lambda = f.mul(&c.e[2], &f.inv(&alpha));
                    let Some(weight) = absolute_orbit_weight(&f, lambda, p) else {
                        continue;
                    };
                    weight
                } else if orbit_quotient {
                    let lambda = f.mul(&c.e[2], &f.inv(&alpha));
                    let lambda_q = f.sigma(&lambda);
                    let lambda_q2 = f.sigma(&lambda_q);
                    let lambda_inv = f.mul(&lambda_q, &lambda_q2);
                    let orbit = [
                        lambda,
                        lambda_q,
                        lambda_q2,
                        lambda_inv,
                        f.sigma(&lambda_inv),
                        f.sigma(&f.sigma(&lambda_inv)),
                    ];
                    debug_assert_ne!(lambda, E6::ONE);
                    debug_assert_eq!(f.mul(&lambda, &lambda_inv), E6::ONE);
                    if orbit.iter().any(|&x| x < lambda) {
                        continue;
                    }
                    if lambda.in_fq() {
                        2
                    } else {
                        6
                    }
                } else {
                    1
                };
                let order = match &squares {
                    Some(table) => exact_order_with_squares(&f, &c, table),
                    None => curve_order(&f, &c, &mut rng),
                };
                let trace = q3 + 1 - order as i128;
                calls += 1;
                if Some(trace) == verify_trace && witness.is_none() {
                    witness = Some(c);
                }
                *weak.entry(trace).or_default() += weight;
                reps += weight;
            }
            (weak, reps, calls, f.muls(), witness)
        })
        .reduce(
            || (BTreeMap::new(), 0, 0, 0, None),
            |(mut a, na, ca, ma, wa), (b, nb, cb, mb, wb)| {
                for (t, c) in b {
                    *a.entry(t).or_default() += c;
                }
                (a, na + nb, ca + cb, ma + mb, wa.or(wb))
            },
        );
    if derive_twists {
        for (t, count) in weak.clone() {
            *weak.entry(-t).or_default() += count;
        }
        reps *= 2;
    }
    let mut out = BufWriter::new(File::create(output).expect("create output"));
    writeln!(out, "p,q,trace,weak_representatives,fundamental_discriminant,frobenius_2_depth,two_split,v2_delta_quarter,class_number_odd,trace_status").unwrap();
    let mut ordinary = 0_u64;
    let lim = 2 * (p as i128).pow(3);
    for t in (-lim..=lim).step_by(4) {
        if t.abs() == lim || t % p as i128 == 0 {
            writeln!(
                out,
                "{p},{q},{t},{},,,,,,nonordinary_or_unrealized",
                weak.get(&t).copied().unwrap_or(0)
            )
            .unwrap();
        } else {
            let (fundamental, depth, split2, v2_d, h_odd) = cm_features(p, t);
            writeln!(
                out,
                "{p},{q},{t},{},{fundamental},{depth},{split2},{v2_d},{},ordinary",
                weak.get(&t).copied().unwrap_or(0),
                u8::from(h_odd)
            )
            .unwrap();
            ordinary += 1;
        }
    }
    out.flush().unwrap();
    if let Some(t) = verify_trace {
        let c = witness.expect("no witness for requested trace");
        let exact = exact_order_by_squares(&f, &c);
        eprintln!(
            "audit p={p} provisional_trace={t} exact_trace={} roots={:?}",
            q3 + 1 - exact as i128,
            c.e
        );
    }
    eprintln!("p={p} q={q} reps={reps} expected={} point_counts={calls} weak_classes={} ordinary_traces={ordinary} all_traces={} muls={} wall_s={:.3} method={} counter={} output={output}",
        2*q*q+2*q, weak.len(), lim / 2 + 1, muls, start.elapsed().as_secs_f64(), if absolute_orbit_quotient { "absolute-orbit" } else if orbit_quotient { "twist-orbit" } else if derive_twists { "twist-derived" } else { "full" }, if exact_squares { "exact-squares" } else { "randomized-two-point" });
    assert_eq!(reps, 2 * q * q + 2 * q);
    assert_eq!(
        calls,
        if absolute_orbit_quotient {
            (q * q + 3 * q + 8) / 12
        } else if orbit_quotient {
            (q * q + q + 4) / 6
        } else if derive_twists {
            q * q + q
        } else {
            2 * q * q + 2 * q
        }
    );
}

#[derive(Clone, Copy)]
struct Row {
    p: u64,
    trace: i128,
    weak: bool,
    depth: u32,
    split: i8,
    h_odd: bool,
}

fn read_rows(path: &str) -> Vec<Row> {
    let file = File::open(path).expect("open census CSV");
    BufReader::new(file)
        .lines()
        .skip(1)
        .map(|line| line.expect("CSV line"))
        .filter_map(|line| {
            let parts: Vec<_> = line.split(',').collect();
            if parts[9] != "ordinary" {
                return None;
            }
            Some(Row {
                p: parts[0].parse().unwrap(),
                trace: parts[2].parse().unwrap(),
                weak: parts[3].parse::<u64>().unwrap() > 0,
                depth: parts[5].parse().unwrap(),
                split: parts[6].parse().unwrap(),
                h_odd: parts[8] == "1",
            })
        })
        .collect()
}

fn key(row: &Row, rule: &str) -> Vec<i128> {
    let delta_quarter = (row.trace / 2).pow(2) - (row.p as i128).pow(6);
    match rule {
        "split_2" => vec![row.split as i128],
        "frobenius_depth" => vec![row.depth as i128],
        "trace_mod_8" => vec![row.trace.rem_euclid(8)],
        "trace_mod_16" => vec![delta_quarter.rem_euclid(16)],
        "class_number_parity" => vec![i128::from(row.h_odd)],
        "combined" => vec![
            row.depth as i128,
            row.split as i128,
            i128::from(row.h_odd),
            delta_quarter.rem_euclid(32),
        ],
        _ => panic!("unknown rule"),
    }
}

fn fit(paths: &[String]) {
    assert_eq!(paths.len(), 3, "fit TRAIN1.csv TRAIN2.csv HOLDOUT.csv");
    let train: Vec<Row> = paths[..2].iter().flat_map(|path| read_rows(path)).collect();
    let holdout = read_rows(&paths[2]);
    for (name, rows) in [("train", &train), ("holdout", &holdout)] {
        let weak = rows.iter().filter(|r| r.weak).count();
        let depth1_weak = rows.iter().filter(|r| r.depth == 1 && r.weak).count();
        let high_missing = rows.iter().filter(|r| r.depth >= 2 && !r.weak).count();
        println!("set={name} n={} weak={weak} depth1_weak={depth1_weak} depth_ge2_missing={high_missing}", rows.len());
    }
    println!("rule,holdout_correct,holdout_total,tp,fp,fn,tn,unseen_feature_rows");
    for rule in [
        "depth_ge_2",
        "split_2",
        "frobenius_depth",
        "trace_mod_8",
        "trace_mod_16",
        "class_number_parity",
        "combined",
    ] {
        let mut counts = BTreeMap::<Vec<i128>, (u64, u64)>::new();
        if rule != "depth_ge_2" {
            for row in &train {
                let entry = counts.entry(key(row, rule)).or_default();
                if row.weak {
                    entry.0 += 1;
                } else {
                    entry.1 += 1;
                }
            }
        }
        let (mut tp, mut fp, mut fn_, mut tn, mut unseen) = (0, 0, 0, 0, 0);
        for row in &holdout {
            let predicted = if rule == "depth_ge_2" {
                row.depth >= 2
            } else if let Some(&(pos, neg)) = counts.get(&key(row, rule)) {
                pos > neg
            } else {
                unseen += 1;
                false
            };
            match (predicted, row.weak) {
                (true, true) => tp += 1,
                (true, false) => fp += 1,
                (false, true) => fn_ += 1,
                (false, false) => tn += 1,
            }
        }
        println!(
            "{rule},{},{},{tp},{fp},{fn_},{tn},{unseen}",
            tp + tn,
            holdout.len()
        );
    }
}

fn sample(p: u64, path: &str, n: u64) {
    let file = File::open(path).expect("open census CSV");
    let weak: HashSet<i128> = BufReader::new(file)
        .lines()
        .skip(1)
        .map(|line| line.unwrap())
        .filter_map(|line| {
            let parts: Vec<_> = line.split(',').collect();
            (parts[3].parse::<u64>().unwrap() > 0).then(|| parts[2].parse().unwrap())
        })
        .collect();
    let f = Fq3::new(p);
    let mut rng = StdRng::seed_from_u64(1 ^ 0xE8AC7);
    let q3 = (p as i128).pow(6);
    let mut hit = 0_u64;
    let mut distinct = HashSet::new();
    for _ in 0..n {
        let c = random_curve(&f, &mut rng);
        let t = q3 + 1 - curve_order(&f, &c, &mut rng) as i128;
        distinct.insert(t);
        if weak.contains(&t) {
            hit += 1;
        }
    }
    let z = 1.96_f64;
    let phat = hit as f64 / n as f64;
    let denom = 1.0 + z * z / n as f64;
    let center = (phat + z * z / (2.0 * n as f64)) / denom;
    let half =
        z * (phat * (1.0 - phat) / n as f64 + z * z / (4.0 * (n as f64).powi(2))).sqrt() / denom;
    println!("p={p} n={n} hit={hit} fraction={phat:.6} wilson95=[{:.6},{:.6}] distinct_traces={} muls={}", center - half, center + half, distinct.len(), f.muls());
}

fn probe_square_lambda(p: u64, n: u64) {
    let f = Fq3::new(p);
    let mut rng = StdRng::seed_from_u64(0x1501_5a);
    let q3 = (p as i128).pow(6);
    let mut depth1 = 0_u64;
    let mut ordinary = 0_u64;
    for _ in 0..n {
        let lambda = loop {
            let x = f.random(&mut rng);
            let lambda = f.sq(&x);
            if lambda != E6::ZERO && lambda != E6::ONE {
                break lambda;
            }
        };
        let c = Curve2 {
            e: [E6::ZERO, E6::ONE, lambda],
        };
        let t = q3 + 1 - curve_order(&f, &c, &mut rng) as i128;
        if t % p as i128 != 0 && t.abs() < 2 * (p as i128).pow(3) {
            ordinary += 1;
            if cm_features(p, t).1 == 1 {
                depth1 += 1;
            }
        }
    }
    println!("p={p} sampled_square_lambda={n} ordinary={ordinary} depth1={depth1}");
}

fn visual(paths: &[String]) {
    assert!(paths.len() >= 2, "visual INPUT.csv ... OUTPUT.svg");
    let output = paths.last().unwrap();
    let height = 280 + 58 * (paths.len() - 1);
    let mut svg = format!("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"920\" height=\"{height}\" viewBox=\"0 0 920 {height}\">\n<rect width=\"920\" height=\"{height}\" fill=\"#fbfcff\"/>\n<g font-family=\"Arial,Helvetica,sans-serif\" fill=\"#14233b\"><text x=\"36\" y=\"42\" font-size=\"23\" font-weight=\"700\">ISO-1: trace-level weak-class census</text><text x=\"36\" y=\"69\" font-size=\"14\">Ordinary full-2-torsion trace candidates over F_(p^6); corrected F_(p^2) square classes</text></g>\n<rect x=\"38\" y=\"91\" width=\"228\" height=\"56\" rx=\"8\" fill=\"#e9eef8\" stroke=\"#8fa2c1\"/><text x=\"51\" y=\"115\" font-size=\"14\" font-family=\"Arial Unicode MS,DejaVu Sans,Arial\" fill=\"#14233b\">trace t, t ≡ 2 (mod 4)</text><text x=\"51\" y=\"134\" font-size=\"12\" font-family=\"Arial Unicode MS,DejaVu Sans,Arial\" fill=\"#43536c\">class invariant: Δ = t² − 4p⁶</text>\n<path d=\"M 266 119 H 312\" stroke=\"#677a99\" stroke-width=\"2\"/><path d=\"M 312 119 L 304 115 M 312 119 L 304 123\" stroke=\"#677a99\" stroke-width=\"2\"/>\n<rect x=\"318\" y=\"91\" width=\"250\" height=\"56\" rx=\"8\" fill=\"#ffe7dd\" stroke=\"#d88470\"/><text x=\"330\" y=\"115\" font-size=\"14\" font-family=\"Arial Unicode MS,DejaVu Sans,Arial\" fill=\"#14233b\">v₂(f_Frob) = 1</text><text x=\"330\" y=\"134\" font-size=\"12\" font-family=\"Arial Unicode MS,DejaVu Sans,Arial\" fill=\"#43536c\">no weak reps at p = 11, 13, 17, 37</text>\n<rect x=\"580\" y=\"91\" width=\"301\" height=\"56\" rx=\"8\" fill=\"#e1f4ec\" stroke=\"#60ad87\"/><text x=\"592\" y=\"115\" font-size=\"14\" font-family=\"Arial Unicode MS,DejaVu Sans,Arial\" fill=\"#14233b\">v₂(f_Frob) ≥ 2</text><text x=\"592\" y=\"134\" font-size=\"12\" font-family=\"Arial Unicode MS,DejaVu Sans,Arial\" fill=\"#43536c\">mostly weak; residual zero rows remain</text>\n<text x=\"36\" y=\"185\" font-size=\"16\" font-family=\"Arial Unicode MS,DejaVu Sans,Arial\" font-weight=\"700\" fill=\"#14233b\">Class labels by stratum</text>\n");
    svg = svg.replace(
        "no weak reps at p = 11, 13, 17, 37",
        "proved absent for every odd p",
    );
    svg = svg.replace("y=\"185\"", "y=\"215\"");
    svg.push_str("<rect x=\"38\" y=\"151\" width=\"843\" height=\"27\" rx=\"6\" fill=\"#e7efff\" stroke=\"#8fa2c1\"/><text x=\"49\" y=\"169\" font-size=\"12\" font-family=\"Arial Unicode MS,DejaVu Sans,Arial\" fill=\"#14233b\">Proof: norm-one λ → fourth power → 2-isogenous full 4-torsion → t ≡ ±(p⁶ + 1) mod 16</text>\n");
    for (i, path) in paths[..paths.len() - 1].iter().enumerate() {
        let rows = read_rows(path);
        let p = rows.first().unwrap().p;
        let n = rows.len() as f64;
        let low = rows.iter().filter(|r| r.depth == 1).count();
        let high_weak = rows.iter().filter(|r| r.depth >= 2 && r.weak).count();
        let high_zero = rows.iter().filter(|r| r.depth >= 2 && !r.weak).count();
        let y = 237 + 58 * i;
        let w_low = 620.0 * low as f64 / n;
        let w_weak = 620.0 * high_weak as f64 / n;
        let w_zero = 620.0 * high_zero as f64 / n;
        svg.push_str(&format!("<text x=\"38\" y=\"{}\" font-size=\"16\" font-family=\"Arial Unicode MS,DejaVu Sans,Arial\" fill=\"#14233b\">p = {p}</text><rect x=\"142\" y=\"{}\" width=\"{w_low:.2}\" height=\"28\" fill=\"#e99b84\"/><rect x=\"{:.2}\" y=\"{}\" width=\"{w_weak:.2}\" height=\"28\" fill=\"#65b991\"/><rect x=\"{:.2}\" y=\"{}\" width=\"{w_zero:.2}\" height=\"28\" fill=\"#e6bf64\"/><text x=\"776\" y=\"{}\" font-size=\"12\" font-family=\"Arial Unicode MS,DejaVu Sans,Arial\" fill=\"#14233b\">{high_weak} / {}</text>\n", y + 20, y, 142.0 + w_low, y, 142.0 + w_low + w_weak, y, y + 20, rows.len()));
    }
    let legend_y = 242 + 58 * (paths.len() - 1);
    svg.push_str(&format!("<g font-family=\"Arial Unicode MS,DejaVu Sans,Arial\" font-size=\"12\" fill=\"#14233b\"><rect x=\"142\" y=\"{legend_y}\" width=\"14\" height=\"14\" fill=\"#e99b84\"/><text x=\"163\" y=\"{}\">depth 1, zero</text><rect x=\"325\" y=\"{legend_y}\" width=\"14\" height=\"14\" fill=\"#65b991\"/><text x=\"346\" y=\"{}\">depth ≥ 2, weak</text><rect x=\"535\" y=\"{legend_y}\" width=\"14\" height=\"14\" fill=\"#e6bf64\"/><text x=\"556\" y=\"{}\">depth ≥ 2, zero</text></g></svg>\n", legend_y + 12, legend_y + 12, legend_y + 12));
    std::fs::write(output, svg).expect("write SVG");
}

fn main() {
    let mut args = env::args().skip(1);
    let first = args.next().expect("usage: iso1_class_census P OUTPUT.csv | fit TRAIN1.csv TRAIN2.csv HOLDOUT.csv | sample P CSV N");
    if first == "fit" {
        fit(&args.collect::<Vec<_>>());
        return;
    }
    if first == "sample" {
        let p = args.next().expect("p").parse().unwrap();
        let path = args.next().expect("CSV");
        let n = args.next().expect("N").parse().unwrap();
        assert!(args.next().is_none());
        sample(p, &path, n);
        return;
    }
    if first == "probe-square-lambda" {
        let p = args.next().expect("p").parse().unwrap();
        let n = args.next().expect("N").parse().unwrap();
        assert!(args.next().is_none());
        probe_square_lambda(p, n);
        return;
    }
    if first == "visual" {
        visual(&args.collect::<Vec<_>>());
        return;
    }
    let p: u64 = first.parse().expect("prime p");
    let output = args.next().expect("output CSV");
    let rest: Vec<_> = args.collect();
    let derive_twists = rest.iter().any(|x| x == "--derive-twists");
    let orbit_quotient = rest.iter().any(|x| x == "--orbit-quotient");
    let absolute_orbit_quotient = rest.iter().any(|x| x == "--absolute-orbit-quotient");
    let exact_squares = rest.iter().any(|x| x == "--exact-squares");
    let positional: Vec<_> = rest
        .iter()
        .filter(|x| {
            x.as_str() != "--derive-twists"
                && x.as_str() != "--orbit-quotient"
                && x.as_str() != "--absolute-orbit-quotient"
                && x.as_str() != "--exact-squares"
        })
        .collect();
    assert!(
        positional.len() <= 1
            && rest.len()
                == positional.len()
                    + usize::from(derive_twists)
                    + usize::from(orbit_quotient)
                    + usize::from(absolute_orbit_quotient)
                    + usize::from(exact_squares),
        "unexpected argument"
    );
    let verify_trace = positional.first().map(|t| t.parse().expect("verify trace"));
    census(
        p,
        &output,
        verify_trace,
        derive_twists,
        orbit_quotient,
        absolute_orbit_quotient,
        exact_squares,
    );
}
