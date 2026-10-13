//! Volcano level and exact 2-sum factor-base yield of every member of one ell-isogeny class of
//! ordinary binary curves y^2 + xy = x^3 + a2 x^2 + a6 over GF(2^n): one JSON line per member.
//!
//! 1. sample a6 (and the twist a2) from `--seed` until the Frobenius conductor f_pi has an odd
//!    prime ell <= `--ell-max` with v_ell(f_pi) >= 2 (smallest such ell of the first hit);
//! 2. enumerate the connected component of j(E) in the ell-graph (roots of Phi_ell mod 2),
//!    at most `--max-vertices` vertices, choosing for each j the twist with the class order;
//! 3. per vertex: Kohel level (non-backtracking walks to the floor) in the volcano of every odd
//!    prime ell' | f_pi, ell' <= ell-max, the factor base V = {x < 2^l}, the exact 2-sum yield,
//!    the order-r subgroup yield and the cofactor-balance count (`binary_volcano.rs`).
//!
//! cargo run --release --example binary_level_yield -- --n 51 --l 12 --ell-max 13 \
//!     --max-samples 2000 --max-vertices 3000 --seed 1 --out results/binary_level_yield_n51.jsonl
use isogeny_algos::binary::BinCurve;
use isogeny_algos::binary_volcano::*;
use isogeny_algos::curve::Pt;
use isogeny_algos::field::{Field, Rng};
use isogeny_algos::gf2n::GF2n;
use isogeny_algos::json::Json;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

struct Args {
    n: u32,
    l: u32,
    ell_max: u64,
    max_samples: usize,
    max_vertices: usize,
    seed: u64,
    out: String,
    threads: usize,
}

fn parse_args() -> Args {
    let mut a = Args {
        n: 51,
        l: 12,
        ell_max: 13,
        max_samples: 2000,
        max_vertices: 3000,
        seed: 1,
        out: String::new(),
        threads: std::thread::available_parallelism().map_or(1, |p| p.get()),
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let key = argv[i].as_str();
        let val = argv.get(i + 1).unwrap_or_else(|| {
            eprintln!("missing value for {key}");
            std::process::exit(2)
        });
        match key {
            "--n" => a.n = val.parse().expect("--n"),
            "--l" => a.l = val.parse().expect("--l"),
            "--ell-max" => a.ell_max = val.parse().expect("--ell-max"),
            "--max-samples" => a.max_samples = val.parse().expect("--max-samples"),
            "--max-vertices" => a.max_vertices = val.parse().expect("--max-vertices"),
            "--seed" => a.seed = val.parse().expect("--seed"),
            "--out" => a.out = val.clone(),
            "--threads" => a.threads = val.parse().expect("--threads"),
            _ => {
                eprintln!("unknown option {key}");
                std::process::exit(2)
            }
        }
        i += 2;
    }
    if a.out.is_empty() {
        a.out = format!("results/binary_level_yield_n{}.jsonl", a.n);
    }
    assert!(a.l <= a.n, "--l must not exceed --n");
    a
}

fn hex(v: u64) -> Json {
    Json::str(format!("{v:#x}"))
}
fn num<T: TryInto<i64> + Copy>(v: T) -> Json {
    v.try_into().map_or(Json::Null, Json::Num)
}

/// Per-ell volcano data of one vertex: (ell, height, level, rational isogenies).
struct Position {
    ell: u64,
    height: u32,
    level: Option<u32>,
    rational: usize,
}

fn main() {
    let args = parse_args();
    let t0 = Instant::now();
    let f = GF2n::new(args.n);
    let q = f.size();
    let tau = trace_one_element(&f);
    let a2_pair = (0u64, tau);
    let modulus = (1u64 << args.n) | f.red;
    let mut rng = Rng::new(args.seed);
    let admissible = admissible_ells(q, args.ell_max);
    eprintln!(
        "n = {}, modulus = {modulus:#x}, twist a2 in {{0, {tau:#x}}}, odd ell <= {} that can divide D: {admissible:?}",
        args.n, args.ell_max
    );

    // 1. sample until some odd ell <= ell_max has v_ell(f_pi) >= 2
    let mut samples = 0usize;
    let mut found: Option<(BinCurve, u128, FrobDisc, u64)> = None;
    while samples < args.max_samples {
        samples += 1;
        let a6 = loop {
            let v = f.random(&mut rng);
            if v != 0 {
                break v;
            }
        };
        let a2 = if rng.next() & 1 == 1 { tau } else { 0 };
        let e = BinCurve::new(a2, 0, a6);
        let order = e.order(&f, &mut rng);
        let disc = frobenius_disc(q, order);
        let hit = disc.odd_heights(args.ell_max).into_iter().find(|&(_, h)| h >= 2);
        if samples % 100 == 0 {
            eprintln!(
                "  sample {samples}: t = {}, f_pi = {} ({:?}) [{:.1}s]",
                disc.trace,
                disc.f_pi,
                disc.f_pi_factors,
                t0.elapsed().as_secs_f64()
            );
        }
        if let Some((ell, _)) = hit {
            found = Some((e, order, disc, ell));
            break;
        }
    }
    let t_sample = t0.elapsed();
    let Some((e0, order, disc, ell)) = found else {
        eprintln!(
            "no class with v_ell(f_pi) >= 2 (odd ell <= {}) in {samples} samples [{:.1}s]; nothing written",
            args.ell_max,
            t_sample.as_secs_f64()
        );
        std::process::exit(3);
    };
    let height = disc.height(ell);
    let j0 = e0.j(&f);
    eprintln!(
        "sample {samples}: a2 = {:#x}, a6 = {:#x}, j = {j0:#x}, #E = {order}, t = {}, D = {}, D_K = {}, f_pi = {} = {:?}, ell = {ell}, height = {height} [{:.1}s]",
        e0.a2, e0.a6, disc.trace, disc.d, disc.d_k, disc.f_pi, disc.f_pi_factors, t_sample.as_secs_f64()
    );
    let (r, cofactor) = largest_prime_factor(order).map_or((None, None), |(r, c)| (Some(r), Some(c)));
    let order_factors: Vec<(u64, u32)> = if order <= u64::MAX as u128 {
        isogeny_algos::field::factor_u64(order as u64)
    } else {
        vec![]
    };

    // 2. Phi_ell mod 2 for the chosen ell and for every other reported ell'
    let t1 = Instant::now();
    let g = EllGraph::new(&f, ell, height);
    let others: Vec<EllGraph> = disc
        .odd_heights(args.ell_max)
        .into_iter()
        .filter(|&(l, _)| l != ell)
        .map(|(l, h)| EllGraph::new(&f, l, h))
        .collect();
    eprintln!("Phi_ell mod 2 built for ell in {:?} [{:.2}s]",
        std::iter::once(ell).chain(others.iter().map(|g| g.ell)).collect::<Vec<_>>(), t1.elapsed().as_secs_f64());

    // 3. the component
    let t2 = Instant::now();
    let comp = enumerate_component(&g, j0, e0.a2, order, a2_pair, args.max_vertices, &mut rng, |done, total| {
        if done % 200 == 0 {
            eprintln!("  component: {done} expanded / {total} discovered [{:.1}s]", t2.elapsed().as_secs_f64());
        }
    });
    let t_component = t2.elapsed();
    let mut level_counts: Vec<usize> = vec![0; height as usize + 1];
    let mut level_unknown = 0;
    for v in &comp.vertices {
        match v.level {
            Some(l) => level_counts[l as usize] += 1,
            None => level_unknown += 1,
        }
    }
    eprintln!(
        "component: {} vertices (truncated: {}), levels 0..={height}: {level_counts:?}, unknown: {level_unknown} [{:.1}s]",
        comp.vertices.len(),
        comp.truncated,
        t_component.as_secs_f64()
    );
    // positions in the other volcanoes (walks only; cheap)
    let mut other_positions: Vec<Vec<Position>> = Vec::with_capacity(comp.vertices.len());
    for v in &comp.vertices {
        let mut ps = vec![];
        for og in &others {
            let ns = og.neighbors(v.j, &mut rng);
            ps.push(Position {
                ell: og.ell,
                height: og.height,
                level: og.level(v.j, &mut rng),
                rational: ns.len(),
            });
        }
        other_positions.push(ps);
    }
    let class_id = format!("n{}-ell{ell}-t{}-j{j0:#x}", args.n, disc.trace);

    // 4. yields, in parallel, one row per vertex as it completes
    let file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&args.out)
        .unwrap_or_else(|e| panic!("cannot open {}: {e}", args.out));
    let out = Mutex::new(std::io::BufWriter::new(file));
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    let t3 = Instant::now();
    let nv = comp.vertices.len();
    let threads = args.threads.clamp(1, nv.max(1));
    std::thread::scope(|s| {
        for tid in 0..threads {
            let (f, comp, other_positions, disc, order_factors, admissible, args) =
                (&f, &comp, &other_positions, &disc, &order_factors, &admissible, &args);
            let (out, next, done, class_id) = (&out, &next, &done, &class_id);
            s.spawn(move || {
                let mut rng = Rng::new(args.seed ^ (0x5eed << 8) ^ tid as u64);
                loop {
                    let i = next.fetch_add(1, Ordering::SeqCst);
                    if i >= nv {
                        break;
                    }
                    let v = &comp.vertices[i];
                    let tv = Instant::now();
                    let a2 = v.a2;
                    let a6 = f.inv(v.j);
                    let mut row: Vec<(&str, Json)> = vec![
                        ("n", num(args.n)),
                        ("modulus", hex(modulus)),
                        ("a2", a2.map_or(Json::Null, hex)),
                        ("a6", hex(a6)),
                        ("j", hex(v.j)),
                        ("order", num(order)),
                        ("trace", num(disc.trace)),
                        ("cofactor", cofactor.map_or(Json::Null, num)),
                        ("r", r.map_or(Json::Null, num)),
                        ("order_factored", Json::Bool(r.is_some())),
                        (
                            "order_factors",
                            Json::Arr(order_factors.iter().map(|&(p, e)| Json::Arr(vec![num(p), num(e)])).collect()),
                        ),
                        ("D", num(disc.d)),
                        ("D_K", num(disc.d_k)),
                        ("f_pi", num(disc.f_pi)),
                        (
                            "f_pi_factored",
                            Json::Arr(disc.f_pi_factors.iter().map(|&(p, e)| Json::Arr(vec![num(p), num(e)])).collect()),
                        ),
                        ("f_pi_unfactored", disc.unfactored.map_or(Json::Null, num)),
                        ("admissible_ells", Json::Arr(admissible.iter().map(|&l| num(l)).collect())),
                        ("ell", num(ell)),
                        ("height_ell", num(height)),
                        ("v_ell_f_E", v.level.map_or(Json::Null, num)),
                        ("rational_ell_isogenies", num(v.rational)),
                        ("kronecker_DK_ell", num(kronecker_odd(disc.d_k, ell))),
                    ];
                    // every reported ell' (the chosen one first)
                    let mut volc: Vec<(String, Json)> = vec![(
                        ell.to_string(),
                        Json::obj(vec![
                            ("v_ell_f_E", v.level.map_or(Json::Null, num)),
                            ("height_ell", num(height)),
                            ("rational_ell_isogenies", num(v.rational)),
                            ("kronecker_DK_ell", num(kronecker_odd(disc.d_k, ell))),
                        ]),
                    )];
                    for p in &other_positions[i] {
                        volc.push((
                            p.ell.to_string(),
                            Json::obj(vec![
                                ("v_ell_f_E", p.level.map_or(Json::Null, num)),
                                ("height_ell", num(p.height)),
                                ("rational_ell_isogenies", num(p.rational)),
                                ("kronecker_DK_ell", num(kronecker_odd(disc.d_k, p.ell))),
                            ]),
                        ));
                    }
                    row.push(("volcano", Json::Obj(volc)));
                    // Kohel consistency: ell + 1 rational isogenies above the floor, 1 on it
                    let level_consistent = v.level.is_some_and(|l| {
                        if l == height {
                            v.rational == 1
                        } else {
                            v.rational == ell as usize + 1
                        }
                    }) && other_positions[i].iter().all(|p| {
                        p.level.is_some_and(|l| if l == p.height { p.rational == 1 } else { p.rational == p.ell as usize + 1 })
                    });
                    row.push(("level_consistent", Json::Bool(level_consistent)));
                    match a2 {
                        None => {
                            row.push(("order_check", Json::Bool(false)));
                            row.push(("error", Json::str("no twist of j has the class order")));
                        }
                        Some(a2) => {
                            let e = BinCurve::new(a2, 0, a6);
                            // the order once more, on two fresh random points
                            let ok = (0..2).all(|_| e.mul(f, &e.random_point(f, &mut rng), order) == Pt::Inf)
                                && e.j(f) == v.j;
                            row.push(("order_check", Json::Bool(ok)));
                            let fb = factor_base(f, &e, args.l);
                            let y = two_sum_yield(f, &e, &fb, r);
                            let qf = q as f64;
                            row.extend([
                                ("fb_dim", num(args.l)),
                                ("fb_points", num(fb.points())),
                                ("fb_abscissae", num(fb.abscissae())),
                                ("fb_pairs", num(y.pairs)),
                                ("yield2_count", num(y.count)),
                                ("yield2_exact", Json::Float(y.count as f64 / qf)),
                                ("yield2_subgroup_count", y.subgroup_count.map_or(Json::Null, num)),
                                (
                                    "yield2_subgroup_exact",
                                    y.subgroup_count.map_or(Json::Null, |c| Json::Float(c as f64 / qf)),
                                ),
                                ("subgroup_combos", y.subgroup_combos.map_or(Json::Null, num)),
                                ("cofactor_balance_C", y.cofactor_balance_c.map_or(Json::Null, num)),
                            ]);
                            if r.is_none() {
                                row.push(("subgroup_note", Json::str("order not factored: no r")));
                            }
                        }
                    }
                    row.extend([
                        ("class_component_id", Json::str(class_id)),
                        ("member_index", num(i)),
                        ("parent_index", v.parent.map_or(Json::Null, num)),
                        ("component_size", num(nv)),
                        ("component_truncated", Json::Bool(comp.truncated)),
                        ("samples_tried", num(samples)),
                        ("seed", num(args.seed)),
                        ("yield_ms", num(tv.elapsed().as_millis() as u64)),
                    ]);
                    let line = Json::obj(row).dump();
                    {
                        let mut w = out.lock().unwrap();
                        writeln!(w, "{line}").unwrap();
                        w.flush().unwrap();
                    }
                    let d = done.fetch_add(1, Ordering::SeqCst) + 1;
                    if d % 50 == 0 || d == nv {
                        eprintln!("  yields: {d}/{nv} rows [{:.1}s]", t3.elapsed().as_secs_f64());
                    }
                }
            });
        }
    });
    let t_yield = t3.elapsed();
    eprintln!(
        "done: {} rows -> {} | sampling {:.1}s ({samples} samples), Phi {:.2}s, component {:.1}s, yields {:.1}s on {threads} threads, total {:.1}s",
        nv,
        args.out,
        t_sample.as_secs_f64(),
        t1.elapsed().as_secs_f64() - t_component.as_secs_f64() - t_yield.as_secs_f64(),
        t_component.as_secs_f64(),
        t_yield.as_secs_f64(),
        t0.elapsed().as_secs_f64()
    );
}
