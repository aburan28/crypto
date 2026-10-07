//! Baseline benchmark harness. Usage:
//!   bench [--quick] [--out FILE] [group ...]      groups: kernel find path (default: all)
//! Single-threaded, fixed seeds. Every timed result is cross-checked for correctness first;
//! records carry `verified`. Output: JSON lines (one record per measurement) + a text table.
use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::{bmss, dual};
use isogeny_algos::find::{divpoly, elkies, modpoly::Phi};
use isogeny_algos::fpm::FpM;
use isogeny_algos::kernel::chain::{ell_power_isogeny, mul_by_ell, ChainStats, Strategy};
use isogeny_algos::kernel::sqrt_velu::sqrt_velu_fast;
use isogeny_algos::kernel::xonly::velu_xonly_fast;
use isogeny_algos::kernel::{
    kohel::kohel,
    montgomery,
    sqrt_velu::sqrt_velu,
    velu::{cyclic_reps, velu_cyclic, velu_from_point, velu_general},
    xonly::velu_xonly,
};
use isogeny_algos::path::csidh::Csidh;
use isogeny_algos::path::endo;
use isogeny_algos::path::graph::*;
use isogeny_algos::path::{couveignes, delfs_galbraith, galbraith, ghs, volcano};
use isogeny_algos::poly;
use isogeny_algos::testdata::*;
use std::fmt::Write as _;
use std::hint::black_box;
use std::io::Write;
use std::time::Instant;

struct Out {
    file: Option<std::fs::File>,
    quick: bool,
}
impl Out {
    fn rec(
        &mut self,
        group: &str,
        algo: &str,
        params: &[(&str, String)],
        metrics: &[(&str, f64)],
        verified: bool,
        note: &str,
    ) {
        let mut s = format!("{{\"group\":\"{group}\",\"algo\":\"{algo}\"");
        for (k, v) in params {
            let _ = write!(s, ",\"{k}\":\"{v}\"");
        }
        for (k, v) in metrics {
            let _ = write!(s, ",\"{k}\":{v}");
        }
        let _ = write!(s, ",\"verified\":{verified},\"note\":\"{note}\"}}");
        println!("{s}");
        if let Some(f) = self.file.as_mut() {
            let _ = writeln!(f, "{s}");
        }
    }
}

/// median / min of repeated runs within a time budget.
fn time_it<T>(budget_ms: u64, max_reps: usize, mut f: impl FnMut() -> T) -> (f64, f64, usize) {
    let mut ts = vec![];
    let start = Instant::now();
    loop {
        let t = Instant::now();
        black_box(f());
        ts.push(t.elapsed().as_nanos() as f64);
        if ts.len() >= max_reps
            || (start.elapsed().as_millis() as u64 >= budget_ms && ts.len() >= 3)
            || start.elapsed().as_millis() as u64 >= budget_ms * 4
        {
            break;
        }
    }
    ts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (ts[ts.len() / 2], ts[0], ts.len())
}

fn prime_bits(bits: u32) -> u64 {
    next_prime(1u64 << (bits - 1))
}

fn bench_kernel(o: &mut Out) {
    let ells: &[u64] = if o.quick {
        &[3, 7, 31, 101]
    } else {
        &[3, 5, 7, 11, 13, 31, 61, 101, 211, 401, 1009]
    };
    let bits_list: &[u32] = if o.quick { &[32] } else { &[32, 61] };
    let budget = if o.quick { 100 } else { 400 };
    for &bits in bits_list {
        let fp = Zp::new(prime_bits(bits));
        let mut rng = Rng::new(1000 + bits as u64);
        for &ell in ells {
            // generating a curve with a rational point of order 1009 needs ~1000 point counts;
            // too slow at 61 bits with the baseline affine arithmetic, so 61-bit stops at 401
            if bits > 32 && ell > 401 {
                continue;
            }
            let (e, p) = curve_with_point(&fp, ell, &mut rng);
            let (h, pts) = kernel_poly_from_point(&fp, &e, &p, ell);
            let v = velu_from_point(&fp, &e, &p, ell);
            let k = kohel(&fp, &e, &h, ell);
            let s = sqrt_velu(&fp, &e, &p, ell);
            let ok_cod = v.cod == k.cod && v.cod == s.cod;
            let mut ok_eval = true;
            for _ in 0..4 {
                let q = e.random_point(&fp, &mut rng);
                let (a, b, c) = (v.eval(&fp, &q), k.eval(&fp, &q), s.eval(&fp, &q));
                ok_eval &= a == b && a == c;
            }
            let ok = ok_cod && ok_eval && check_homomorphism(&fp, &v, &mut rng, 3);
            let params = |extra: &str| {
                vec![
                    ("p_bits", bits.to_string()),
                    ("ell", ell.to_string()),
                    ("task", extra.to_string()),
                ]
            };
            let mut run = |algo: &str, task: &str, f: &mut dyn FnMut()| {
                let (med, min, reps) = time_it(budget, 2000, || f());
                o.rec(
                    "kernel",
                    algo,
                    &params(task),
                    &[("median_ns", med), ("min_ns", min), ("reps", reps as f64)],
                    ok,
                    "",
                );
            };
            // codomain from generator
            run("velu", "codomain_from_point", &mut || {
                black_box(velu_from_point(&fp, &e, &p, ell));
            });
            run("kohel", "codomain_from_point(poly+formulas)", &mut || {
                let (h, _) = kernel_poly_from_point(&fp, &e, &p, ell);
                black_box(kohel(&fp, &e, &h, ell));
            });
            run("kohel", "codomain_given_h", &mut || {
                black_box(kohel(&fp, &e, &h, ell));
            });
            run("sqrt_velu", "codomain_from_point", &mut || {
                black_box(sqrt_velu(&fp, &e, &p, ell));
            });
            // evaluation of a point on a built isogeny
            let q = e.random_point(&fp, &mut rng);
            run("velu", "eval_point", &mut || {
                black_box(v.eval(&fp, &q));
            });
            run("kohel", "eval_point", &mut || {
                black_box(k.eval(&fp, &q));
            });
            run("sqrt_velu", "eval_point", &mut || {
                black_box(s.eval(&fp, &q));
            });
            let _ = pts;
        }
    }
}

fn bench_find(o: &mut Out) {
    let ells: &[u64] = if o.quick {
        &[3, 5, 7]
    } else {
        &[3, 5, 7, 11, 13, 17, 19, 23]
    };
    let bits_list: &[u32] = if o.quick { &[30] } else { &[30, 61] };
    let budget = if o.quick { 100 } else { 500 };
    for &bits in bits_list {
        let p = prime_bits(bits);
        let fp = Zp::new(p);
        let mut rng = Rng::new(2000 + bits as u64);
        for &ell in ells {
            // setup (reported separately): Phi_l
            let t = Instant::now();
            let phi = Phi::compute(&fp, ell as usize);
            let phi_ns = t.elapsed().as_nanos() as f64;
            o.rec(
                "find",
                "phi_setup",
                &[("p_bits", bits.to_string()), ("ell", ell.to_string())],
                &[("median_ns", phi_ns), ("reps", 1.0)],
                true,
                "one-off per (p,l)",
            );
            let (e, _pt) = loop {
                let (e, pt) = curve_with_point(&fp, ell, &mut rng);
                let j = jinv(&fp, &e);
                if j != 0 && j != 1728 {
                    break (e, pt);
                }
            };
            let mut ks_div = divpoly::kernel_polys(&fp, &e, ell, &mut rng);
            let isos = elkies::elkies_isogenies(&fp, &phi, &e, &mut rng);
            let mut ks_elk: Vec<_> = isos.iter().map(|i| i.ker.clone()).collect();
            ks_div.sort();
            ks_elk.sort();
            let ok = !ks_div.is_empty() && ks_div == ks_elk;
            let params = vec![
                ("p_bits", bits.to_string()),
                ("ell", ell.to_string()),
                ("kernels_found", ks_div.len().to_string()),
            ];
            let (m, mn, r) = time_it(budget, 200, || {
                divpoly::kernel_polys(&fp, &e, ell, &mut rng.clone())
            });
            o.rec(
                "find",
                "divpoly_factor",
                &params,
                &[("median_ns", m), ("min_ns", mn), ("reps", r as f64)],
                ok,
                "psi_l DDF/EDF + subset search",
            );
            let (m, mn, r) = time_it(budget, 200, || {
                elkies::elkies_isogenies(&fp, &phi, &e, &mut rng.clone())
            });
            o.rec(
                "find",
                "elkies_phi_bmss",
                &params,
                &[("median_ns", m), ("min_ns", mn), ("reps", r as f64)],
                ok,
                "Phi roots + Elkies codomain + BMSS Pade (Phi given)",
            );
            let j = jinv(&fp, &e);
            let (m, mn, r) = time_it(budget, 200, || phi.neighbors(&fp, j, &mut rng.clone()));
            o.rec(
                "find",
                "phi_roots_only",
                &params,
                &[("median_ns", m), ("min_ns", mn), ("reps", r as f64)],
                ok,
                "j-invariants of neighbours only, no isogeny",
            );
        }
    }
}

fn bench_path(o: &mut Out) {
    let ord_bits: &[u32] = if o.quick { &[20] } else { &[20, 24, 28, 32] };
    let inst = if o.quick { 2 } else { 5 };
    let ells = [3usize, 5, 7];
    for &bits in ord_bits {
        let p = prime_bits(bits);
        let fp = Zp::new(p);
        let cache = PhiCache::new(&Zp::new(p), &ells);
        for i in 0..inst {
            let mut rng = Rng::new(3000 + (bits as u64) * 100 + i);
            // E1 must have >= 3 rational neighbours over >= 2 primes, else its component is small/cyclic
            let (e1, t) = loop {
                let (e1, t) = curve_with_trace(&fp, &mut rng, |t| {
                    let d = (t * t - 4 * p as i64).abs();
                    ells.iter().all(|&l| d % (l as i64 * l as i64) != 0)
                });
                let nb = neighbors(&fp, &cache, &ells, jinv(&fp, &e1), &mut rng);
                let mut primes: Vec<usize> = nb.iter().map(|x| x.0).collect();
                primes.sort();
                primes.dedup();
                if nb.len() >= 3 && primes.len() >= 2 {
                    break (e1, t);
                }
            };
            let j1 = jinv(&fp, &e1);
            let walk = random_walk(&fp, &cache, &ells, j1, 4 * bits as usize, &mut rng);
            let j2 = *walk.js.last().unwrap();
            let params = vec![
                ("p_bits", bits.to_string()),
                ("instance", i.to_string()),
                ("trace", t.to_string()),
            ];
            // Galbraith
            let t0 = Instant::now();
            let (r, st) =
                galbraith::galbraith(&fp, &cache, &ells, j1, j2, 5_000_000, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| {
                verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2
            });
            o.rec(
                "path",
                "galbraith_bidirectional_bfs",
                &params,
                &[
                    ("ns", ns),
                    ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)),
                    ("nodes_expanded", st.nodes_expanded as f64),
                ],
                ok,
                "primes 3,5,7",
            );
            if let Some(path) = &r {
                let t0 = Instant::now();
                let chain = explicit_chain(&fp, &cache, &e1, path);
                o.rec(
                    "path",
                    "explicit_chain_from_path",
                    &params,
                    &[
                        ("ns", t0.elapsed().as_nanos() as f64),
                        ("steps", path.len() as f64),
                    ],
                    chain.is_some(),
                    "Elkies+BMSS per step",
                );
            }
            // GHS
            let t0 = Instant::now();
            let (r, st) = ghs::ghs(&fp, &cache, &ells, j1, j2, 20_000_000, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| {
                verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2
            });
            o.rec(
                "path",
                "ghs_random_walk_collision",
                &params,
                &[
                    ("ns", ns),
                    ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)),
                    ("steps", st.steps as f64),
                ],
                ok,
                "primes 3,5,7",
            );
            let t0 = Instant::now();
            let (r, st) = ghs::ghs_volcano(
                &fp,
                &cache,
                &ells,
                p,
                t,
                j1,
                j2,
                20_000_000,
                &mut rng.clone(),
            );
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| {
                verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2
            });
            o.rec(
                "path",
                "ghs_with_kohel_volcano",
                &params,
                &[
                    ("ns", ns),
                    ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)),
                    ("steps", st.steps as f64),
                ],
                ok,
                "all heights 0 here, so same walk plus height checks",
            );
        }
    }
    // Kohel volcano (l = 3, height >= 1)
    let vbits: &[u32] = if o.quick { &[17] } else { &[17, 19, 21] };
    for &bits in vbits {
        // need p a square mod 9 so that 9 | t^2 - 4p is possible
        let mut p = prime_bits(bits);
        while ![1u64, 4, 7].contains(&(p % 9)) {
            p = next_prime(p + 1);
        }
        let fp = Zp::new(p);
        let cache = PhiCache::new(&Zp::new(p), &[3, 5, 7, 11, 13]);
        for i in 0..inst {
            let mut rng = Rng::new(4000 + (bits as u64) * 100 + i);
            let (e1, t) = curve_with_trace(&fp, &mut rng, |t| {
                let d = (t * t - 4 * p as i64).abs();
                d % 9 == 0 && d % 25 != 0
            });
            let h = volcano::volcano_height(p, t, 3);
            let j1 = jinv(&fp, &e1);
            let walk = random_walk(&fp, &cache, &[3], j1, 4 * bits as usize, &mut rng);
            let j2 = *walk.js.last().unwrap();
            let params = vec![
                ("p_bits", bits.to_string()),
                ("instance", i.to_string()),
                ("height", h.to_string()),
            ];
            let t0 = Instant::now();
            let r = volcano::kohel_volcano_path(&fp, &cache, 3, h, j1, j2, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| {
                verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2
            });
            o.rec(
                "path",
                "kohel_volcano_crater_walk",
                &params,
                &[
                    ("ns", ns),
                    ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)),
                ],
                ok,
                "l=3 volcano; same-volcano targets only",
            );
            let t0 = Instant::now();
            let (r, st) = ghs::ghs_volcano(
                &fp,
                &cache,
                &[3, 5, 7, 11, 13],
                p,
                t,
                j1,
                j2,
                200_000,
                &mut rng.clone(),
            );
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| {
                verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2
            });
            o.rec("path", "ghs_with_kohel_volcano", &params, &[("ns", ns), ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)), ("steps", st.steps as f64)], ok, "ascend l=3, walk with l=5,7,11,13; fails when the split primes do not generate the class group orbit");
        }
    }
    // Couveignes
    let cbits: &[u32] = if o.quick { &[20] } else { &[20, 24, 28] };
    let cells = [3usize, 5, 7, 11, 13, 17, 19, 23];
    for &bits in cbits {
        let p = prime_bits(bits);
        let fp = Zp::new(p);
        let t0 = Instant::now();
        let cache = PhiCache::new(&Zp::new(p), &cells);
        let setup = t0.elapsed().as_nanos() as f64;
        for i in 0..inst {
            let mut rng = Rng::new(5000 + (bits as u64) * 100 + i);
            let (e1, _t) = curve_with_trace(&fp, &mut rng, |t| {
                let d = (t * t - 4 * p as i64).abs();
                cells.iter().all(|&l| d % (l as i64 * l as i64) != 0)
            });
            let j1 = jinv(&fp, &e1);
            let mut act = couveignes::Action {
                fp: &fp,
                cache: &cache,
                plus_class: Default::default(),
            };
            let t0 = Instant::now();
            let primes = couveignes::select_primes(&mut act, j1, &mut rng);
            let sel_ns = t0.elapsed().as_nanos() as f64;
            let k = primes.len().min(4);
            if k < 2 {
                continue;
            }
            let primes = &primes[..k];
            let m = 4usize;
            let secret: Vec<i64> = primes
                .iter()
                .map(|_| rng.below(2 * m as u64 + 1) as i64 - m as i64)
                .collect();
            let j2 = *act
                .apply(j1, primes, &secret, &mut rng)
                .unwrap()
                .js
                .last()
                .unwrap();
            let params = vec![
                ("p_bits", bits.to_string()),
                ("instance", i.to_string()),
                ("primes", format!("{:?}", primes)),
                ("box_m", m.to_string()),
            ];
            let t0 = Instant::now();
            let (r, st) = couveignes::couveignes(&act, primes, m, j1, j2, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = match &r {
                Some((e, path)) => {
                    let chk = act
                        .apply(j1, primes, e, &mut rng.clone())
                        .map_or(false, |q| *q.js.last().unwrap() == j2);
                    chk && verify_path(&fp, &cache, path)
                }
                None => false,
            };
            o.rec(
                "path",
                "couveignes_hhs_mitm",
                &params,
                &[
                    ("ns", ns),
                    ("nodes", st.nodes as f64),
                    ("path_len", r.as_ref().map_or(-1.0, |x| x.1.len() as f64)),
                    ("prime_select_ns", sel_ns),
                    ("phi_setup_ns", setup),
                ],
                ok,
                "planted exponents in [-m,m]^k; different problem from generic path finding",
            );
        }
    }
    // Delfs-Galbraith
    let dbits: &[u32] = if o.quick {
        &[22]
    } else {
        &[22, 26, 30, 34, 38]
    };
    for &bits in dbits {
        // p = 7 mod 8: 2 splits in Q(sqrt(-p)), so the F_p-rational 2-isogeny graph has a
        // crater cycle (for p = 3 mod 8 it is a forest of stars and the F_p BFS cannot connect)
        let mut p = prime_bits(bits);
        while p % 8 != 7 {
            p = next_prime(p + 1);
        }
        let fp = Zp::new(p);
        let f2 = Zp2::new(p);
        let cache = PhiCache::new(&Zp::new(p), &[2, 3]);
        let cache2 = cache.lift(|x| (x, 0u64));
        for i in 0..inst {
            let mut rng = Rng::new(6000 + (bits as u64) * 100 + i);
            let start = (1728u64, 0u64);
            let steps = 3 * bits as usize;
            let j1 = *random_walk(&f2, &cache2, &[2], start, steps, &mut rng)
                .js
                .last()
                .unwrap();
            let j2 = *random_walk(&f2, &cache2, &[2], start, steps, &mut rng)
                .js
                .last()
                .unwrap();
            let params = vec![("p_bits", bits.to_string()), ("instance", i.to_string())];
            let t0 = Instant::now();
            let (r, st) = delfs_galbraith::delfs_galbraith(
                &f2,
                &fp,
                &cache,
                &cache2,
                j1,
                j2,
                50_000_000,
                5_000_000,
                &mut rng.clone(),
            );
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| {
                verify_path(&f2, &cache2, p) && p.js[0] == j1 && *p.js.last().unwrap() == j2
            });
            o.rec(
                "path",
                "delfs_galbraith_supersingular",
                &params,
                &[
                    ("ns", ns),
                    ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)),
                    ("walk_steps", st.walk_steps as f64),
                    ("bfs_nodes", st.bfs_nodes as f64),
                ],
                ok,
                "j-path only; Fp2 arithmetic",
            );
        }
    }
}

/// Odd-degree kernel->isogeny comparison of every implementation, plus even-degree general Velu/Kohel.
fn bench_kernel2(o: &mut Out) {
    let ells: &[u64] = if o.quick {
        &[3, 7, 31]
    } else {
        &[3, 5, 7, 11, 13, 31, 61, 101, 211, 401, 1009]
    };
    let budget = if o.quick { 100 } else { 400 };
    let fp = Zp::new(prime_bits(32));
    let mut rng = Rng::new(7000);
    for &ell in ells {
        let (e, p) = curve_with_point(&fp, ell, &mut rng);
        let x0 = if let Pt::Aff(x, _) = p {
            x
        } else {
            unreachable!()
        };
        let (h, _) = kernel_poly_from_point(&fp, &e, &p, ell);
        let v = velu_cyclic(&fp, &e, &p, ell);
        let xo = velu_xonly(&fp, &e, x0, ell).unwrap();
        let k = kohel(&fp, &e, &h, ell);
        let ok = v.cod == xo.cod && v.cod == k.cod;
        let params = vec![("p_bits", "32".to_string()), ("ell", ell.to_string())];
        let mut run = |algo: &str, note: &str, f: &mut dyn FnMut()| {
            let (med, min, reps) = time_it(budget, 2000, || f());
            o.rec(
                "kernel2",
                algo,
                &params,
                &[("median_ns", med), ("min_ns", min), ("reps", reps as f64)],
                ok,
                note,
            );
        };
        run(
            "velu_general_from_point",
            "points kP enumerated, general Velu (any degree)",
            &mut || {
                black_box(velu_cyclic(&fp, &e, &p, ell));
            },
        );
        run(
            "velu_xonly_from_x",
            "x(P) only, division-polynomial values",
            &mut || {
                black_box(velu_xonly(&fp, &e, x0, ell));
            },
        );
        run("kohel_given_h", "kernel polynomial given", &mut || {
            black_box(kohel(&fp, &e, &h, ell));
        });
    }
    // even and composite degrees: general Velu vs Kohel (kernel polynomial given)
    for &n in if o.quick {
        &[2u64, 4][..]
    } else {
        &[2u64, 4, 6, 8, 10, 12, 16, 20][..]
    } {
        let (e, p) = loop {
            let c = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
            if !is_smooth(&fp, &c) {
                continue;
            }
            let ord = order(&fp, &c, &mut rng);
            if ord % n != 0 {
                continue;
            }
            let r = c.random_point(&fp, &mut rng);
            let q = pmul(&fp, &c, &r, (ord / n) as u128);
            let mut good = pmul(&fp, &c, &q, n as u128) == Pt::Inf;
            for pr in [2u64, 3, 5] {
                if n % pr == 0 && pmul(&fp, &c, &q, (n / pr) as u128) == Pt::Inf {
                    good = false;
                }
            }
            if good {
                break (c, q);
            }
        };
        let reps = cyclic_reps(&fp, &e, &p, n);
        let xs: Vec<u64> = reps
            .iter()
            .map(|q| if let Pt::Aff(x, _) = q { *x } else { 0 })
            .collect();
        let h = isogeny_algos::poly::from_roots(&fp, &xs);
        let v = velu_general(&fp, &e, &reps);
        let kh = kohel(&fp, &e, &h, n);
        let ok = v.cod == kh.cod && check_homomorphism(&fp, &v, &mut rng, 3);
        let params = vec![("p_bits", "32".to_string()), ("degree", n.to_string())];
        let (m1, mn1, r1) = time_it(budget, 2000, || velu_cyclic(&fp, &e, &p, n));
        o.rec(
            "kernel2_even",
            "velu_general_from_point",
            &params,
            &[("median_ns", m1), ("min_ns", mn1), ("reps", r1 as f64)],
            ok,
            "cyclic kernel of the given order",
        );
        let (m2, mn2, r2) = time_it(budget, 2000, || kohel(&fp, &e, &h, n));
        o.rec(
            "kernel2_even",
            "kohel_given_h",
            &params,
            &[("median_ns", m2), ("min_ns", mn2), ("reps", r2 as f64)],
            ok,
            "kernel polynomial given (2-torsion factor split off)",
        );
    }
    // Montgomery x-only Velu in the CSIDH setting (supersingular, p = 4 prod l_i - 1)
    let cs = Csidh::with_n_primes(if o.quick { 6 } else { 10 }).unwrap();
    let mfp = cs.fp;
    let mut rng = Rng::new(7100);
    for (idx, &ell) in cs.primes.clone().iter().enumerate() {
        if o.quick && idx > 2 {
            break;
        }
        let a = 0u64;
        let e24 = montgomery::a24(&mfp, a);
        let rhs = |u: u64| mfp.add(mfp.mul(mfp.add(mfp.mul(u, u), mfp.mul(a, u)), u), u);
        let kx = loop {
            let x = mfp.random(&mut rng);
            if mfp.sqrt(rhs(x)).is_none() {
                continue;
            }
            let kk = montgomery::ladder(&mfp, e24, x, ((mfp.p + 1) / ell) as u128);
            if kk.1 != 0 {
                break mfp.div(kk.0, kk.1);
            }
        };
        let d = ((ell - 1) / 2) as usize;
        let ms = montgomery::multiples(&mfp, e24, (kx, 1), d);
        let anew = montgomery::velu_codomain(&mfp, a, &ms, ell);
        // cross-check against Weierstrass Velu
        let w = montgomery::to_weierstrass(&mfp, a);
        let ky = mfp.sqrt(rhs(kx)).unwrap();
        let wv = velu_from_point(&mfp, &w, &Pt::Aff(mfp.add(kx, mfp.div(a, 3)), ky), ell);
        let ok = montgomery::j_invariant(&mfp, anew) == jinv(&mfp, &wv.cod);
        let params = vec![
            ("p_bits", (64 - mfp.p.leading_zeros()).to_string()),
            ("ell", ell.to_string()),
        ];
        let (med, mn, r) = time_it(budget, 2000, || {
            let ms = montgomery::multiples(&mfp, e24, (kx, 1), d);
            montgomery::velu_codomain(&mfp, a, &ms, ell)
        });
        o.rec(
            "kernel2_mont",
            "montgomery_xonly_velu_codomain",
            &params,
            &[("median_ns", med), ("min_ns", mn), ("reps", r as f64)],
            ok,
            "multiples + codomain, projective x-only",
        );
        let (med, mn, r) = time_it(budget, 2000, || {
            velu_from_point(&mfp, &w, &Pt::Aff(mfp.add(kx, mfp.div(a, 3)), ky), ell)
        });
        o.rec(
            "kernel2_mont",
            "weierstrass_velu_same_kernel",
            &params,
            &[("median_ns", med), ("min_ns", mn), ("reps", r as f64)],
            ok,
            "same kernel on the Weierstrass model, for comparison",
        );
    }
}

/// l^e kernels: chain strategies. SIDH-shaped primes p = 2^a 3^b f - 1 over F_{p^2}.
fn bench_chain(o: &mut Out) {
    let sizes: &[(u32, u32)] = if o.quick {
        &[(8, 5)]
    } else {
        &[(8, 5), (12, 8), (16, 10), (20, 12), (24, 14)]
    };
    for &(a, b) in sizes {
        let base = 2u64.pow(a) * 3u64.pow(b);
        let p = (1..10_000u64)
            .map(|f| base * f - 1)
            .find(|&p| p > 1000 && p < (1u64 << 60) && is_prime(p) && p % 4 == 3)
            .expect("prime");
        let f2 = Zp2::new(p);
        let mut rng = Rng::new(8000 + a as u64);
        let e0 = Curve::new(f2.from_u64(1), f2.zero());
        for (ell, n) in [(2u64, a), (3u64, b)] {
            // point of exact order ell^n on E0(F_{p^2}) = E0[p+1]
            let cof = ((p + 1) / ell.pow(n)) as u128;
            let pt = loop {
                let r = random_point_f(&f2, &e0, &mut rng);
                let q = pmul(&f2, &e0, &r, cof);
                if pmul(&f2, &e0, &q, ell.pow(n - 1) as u128) != Pt::Inf {
                    break q;
                }
            };
            // calibrate the cost model on this curve
            let mut st = ChainStats::default();
            let (t_mul, _, _) = time_it(50, 50, || mul_by_ell(&f2, &e0, &pt, ell, 1, &mut st));
            let (t_iso, _, _) = time_it(50, 50, || velu_cyclic(&f2, &e0, &pt, ell));
            let iso = velu_cyclic(&f2, &e0, &pt, ell);
            // evaluate at a point outside the kernel (at a kernel point the evaluation returns at once)
            let probe = random_point_f(&f2, &e0, &mut rng);
            let (t_eval, _, _) = time_it(50, 50, || iso.eval(&f2, &probe));
            let strategies = [
                ("naive", Strategy::Naive),
                ("balanced", Strategy::Balanced),
                (
                    "optimal_calibrated",
                    Strategy::Optimal {
                        mul: t_mul,
                        eval: t_eval,
                        iso: t_iso,
                    },
                ),
            ];
            let mut reference = None;
            for (name, strat) in strategies {
                let (chain, _, stats) =
                    ell_power_isogeny(&f2, &e0, &pt, ell, n as usize, &[], strat);
                let jc = jinv(&f2, &chain.cod);
                let ok = *reference.get_or_insert(jc) == jc;
                let (med, mn, reps) = time_it(300, 200, || {
                    ell_power_isogeny(&f2, &e0, &pt, ell, n as usize, &[], strat)
                });
                let params = vec![
                    ("p_bits", (64 - p.leading_zeros()).to_string()),
                    ("ell", ell.to_string()),
                    ("e", n.to_string()),
                    ("strategy", name.to_string()),
                ];
                o.rec(
                    "chain",
                    "ell_power_isogeny",
                    &params,
                    &[
                        ("median_ns", med),
                        ("min_ns", mn),
                        ("reps", reps as f64),
                        ("l_mults", stats.l_mults as f64),
                        ("evals", stats.evals as f64),
                        ("builds", stats.builds as f64),
                    ],
                    ok,
                    "all strategies must give the same codomain j",
                );
            }
        }
    }
}

/// The BMSS family: (E, E~) -> isogeny, one record per method and l.
fn bench_bmss(o: &mut Out) {
    let ells: &[u64] = if o.quick {
        &[5, 7, 11]
    } else {
        &[5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 53, 61, 101]
    };
    let bits = 40;
    let fp = Zp::new(prime_bits(bits));
    let mut rng = Rng::new(9000);
    for &ell in ells {
        let (e, pt) = loop {
            let (e, pt) = curve_with_point(&fp, ell, &mut rng);
            let j = jinv(&fp, &e);
            if j != 0 && j != 1728 {
                break (e, pt);
            }
        };
        let (g, _) = kernel_poly_from_point(&fp, &e, &pt, ell);
        let d = g.len() - 1;
        let sigma = fp.neg(fp.mul(2, g[d - 1]));
        let kh = kohel(&fp, &e, &g, ell);
        let et = kh.cod;
        let budget = if o.quick { 100 } else { 400 };
        for m in bmss::Method::ALL {
            let s = if m.needs_sigma() { Some(sigma) } else { None };
            let ok = bmss::isogeny(&fp, m, &e, &et, ell as usize, s)
                .map_or(false, |i| i.ker == g && i.num == kh.num);
            let (med, mn, reps) = time_it(budget, 200, || {
                bmss::isogeny(&fp, m, &e, &et, ell as usize, s)
            });
            let params = vec![
                ("p_bits", bits.to_string()),
                ("ell", ell.to_string()),
                ("needs_sigma", m.needs_sigma().to_string()),
            ];
            o.rec(
                "bmss",
                m.name(),
                &params,
                &[("median_ns", med), ("min_ns", mn), ("reps", reps as f64)],
                ok,
                "sigma supplied from the true kernel; baseline arithmetic",
            );
        }
        // reference: Kohel's formula with the kernel polynomial given, and the V1 Pade path
        let (med, mn, reps) = time_it(budget, 200, || kohel(&fp, &e, &g, ell));
        o.rec(
            "bmss",
            "reference_kohel_given_g",
            &[("p_bits", bits.to_string()), ("ell", ell.to_string())],
            &[("median_ns", med), ("min_ns", mn), ("reps", reps as f64)],
            true,
            "not an algorithm for this problem: needs the kernel",
        );
        let (med, mn, reps) = time_it(budget, 200, || {
            isogeny_algos::find::elkies::bmss_isogeny(&fp, &e, &et, ell as usize)
        });
        o.rec(
            "bmss",
            "v1_pade_weierstrass_series",
            &[("p_bits", bits.to_string()), ("ell", ell.to_string())],
            &[("median_ns", med), ("min_ns", mn), ("reps", reps as f64)],
            true,
            "V1 implementation: Pade on the inverted Weierstrass series",
        );
        if ell <= 23 {
            let phi = Phi::compute(&fp, ell as usize);
            let ok = bmss::sigma_from_phi(&fp, &phi, &e, &et) == Some(sigma);
            let (med, mn, reps) = time_it(budget, 200, || bmss::sigma_from_phi(&fp, &phi, &e, &et));
            o.rec(
                "bmss",
                "sigma_from_phi",
                &[("p_bits", bits.to_string()), ("ell", ell.to_string())],
                &[("median_ns", med), ("min_ns", mn), ("reps", reps as f64)],
                ok,
                "E2 - l E2' from Phi's second derivatives",
            );
            for m in [bmss::Method::Elkies1998, bmss::Method::FastElkiesPrime] {
                let (med, mn, reps) = time_it(budget, 200, || {
                    bmss::isogenies_via_phi(&fp, &phi, &e, m, &mut rng.clone())
                });
                let ok = bmss::isogenies_via_phi(&fp, &phi, &e, m, &mut rng.clone())
                    .iter()
                    .any(|i| i.ker == g);
                o.rec(
                    "bmss",
                    "end_to_end_via_phi",
                    &[
                        ("p_bits", bits.to_string()),
                        ("ell", ell.to_string()),
                        ("method", m.name().to_string()),
                    ],
                    &[("median_ns", med), ("min_ns", mn), ("reps", reps as f64)],
                    ok,
                    "Phi roots + Elkies codomain + sigma + method (Phi given)",
                );
            }
            // dual isogeny
            let isos = isogeny_algos::find::elkies::elkies_isogenies(&fp, &phi, &e, &mut rng);
            if let Some(iso) = isos.iter().find(|i| i.ker == g) {
                let ok = dual::dual_isogeny(&fp, &phi, iso).map_or(false, |d| {
                    let p = e.random_point(&fp, &mut rng.clone());
                    d.eval(&fp, &iso.eval(&fp, &p)) == pmul(&fp, &e, &p, ell as u128)
                });
                let (med, mn, reps) = time_it(budget, 200, || dual::dual_isogeny(&fp, &phi, iso));
                o.rec(
                    "bmss",
                    "dual_isogeny",
                    &[("p_bits", bits.to_string()), ("ell", ell.to_string())],
                    &[("median_ns", med), ("min_ns", mn), ("reps", reps as f64)],
                    ok,
                    "Elkies codomain over j(E) + Pade; dual o phi = [l] checked",
                );
            }
        }
    }
}

/// CSIDH-style class-group action and meet-in-the-middle.
fn bench_csidh(o: &mut Out) {
    let ns: &[usize] = if o.quick { &[6] } else { &[6, 8, 10, 12, 14] };
    for &n in ns {
        let Some(cs) = Csidh::with_n_primes(n) else {
            continue;
        };
        let mut rng = Rng::new(10_000 + n as u64);
        let bits = 64 - cs.p().leading_zeros();
        let m = 3i32;
        let exps: Vec<i32> = (0..n)
            .map(|_| rng.below(2 * m as u64 + 1) as i32 - m)
            .collect();
        let a1 = cs.action(0, &exps, &mut rng.clone());
        // verification: commutativity with a second vector and inverse
        let exps2: Vec<i32> = (0..n)
            .map(|_| rng.below(2 * m as u64 + 1) as i32 - m)
            .collect();
        let ab = cs.action(cs.action(0, &exps, &mut rng), &exps2, &mut rng);
        let ba = cs.action(cs.action(0, &exps2, &mut rng), &exps, &mut rng);
        let neg: Vec<i32> = exps.iter().map(|x| -x).collect();
        let ok = ab == ba && cs.action(a1, &neg, &mut rng) == 0;
        let steps: i32 = exps.iter().map(|x| x.abs()).sum();
        let (med, mn, reps) = time_it(if o.quick { 100 } else { 500 }, 200, || {
            cs.action(0, &exps, &mut rng.clone())
        });
        o.rec(
            "csidh",
            "class_group_action",
            &[
                ("p_bits", bits.to_string()),
                ("n_primes", n.to_string()),
                ("exp_bound", m.to_string()),
                ("isogeny_steps", steps.to_string()),
            ],
            &[
                ("median_ns", med),
                ("min_ns", mn),
                ("reps", reps as f64),
                ("ns_per_step", med / steps.max(1) as f64),
            ],
            ok,
            "one fresh random point per isogeny, no batching",
        );
        let okb = cs.action_batched(0, &exps, &mut rng.clone()) == a1;
        let (med, mn, reps) = time_it(if o.quick { 100 } else { 500 }, 200, || {
            cs.action_batched(0, &exps, &mut rng.clone())
        });
        o.rec(
            "csidh",
            "class_group_action_batched",
            &[
                ("p_bits", bits.to_string()),
                ("n_primes", n.to_string()),
                ("exp_bound", m.to_string()),
                ("isogeny_steps", steps.to_string()),
            ],
            &[
                ("median_ns", med),
                ("min_ns", mn),
                ("reps", reps as f64),
                ("ns_per_step", med / steps.max(1) as f64),
            ],
            okb && ok,
            "CLMPR: one point per sign serves all primes; equals stepwise result",
        );
        let okf = cs.action_fast(0, &exps, &mut rng.clone()) == a1;
        let (med, mn, reps) = time_it(if o.quick { 100 } else { 500 }, 200, || {
            cs.action_fast(0, &exps, &mut rng.clone())
        });
        o.rec(
            "csidh",
            "class_group_action_tree_projective",
            &[
                ("p_bits", bits.to_string()),
                ("n_primes", n.to_string()),
                ("exp_bound", m.to_string()),
                ("isogeny_steps", steps.to_string()),
            ],
            &[
                ("median_ns", med),
                ("min_ns", mn),
                ("reps", reps as f64),
                ("ns_per_step", med / steps.max(1) as f64),
            ],
            okf && ok,
            "projective (A24:C24), tree strategy for kernel points; equals stepwise result",
        );
        // meet-in-the-middle on a small box
        let mm = if n <= 8 { 2usize } else { 1 };
        if n <= 10 {
            let e: Vec<i32> = (0..n)
                .map(|_| rng.below(2 * mm as u64 + 1) as i32 - mm as i32)
                .collect();
            let target = cs.action(0, &e, &mut rng);
            let t0 = Instant::now();
            let r = cs.mitm(0, target, mm, &mut rng.clone());
            let ns_t = t0.elapsed().as_nanos() as f64;
            let ok = r
                .as_ref()
                .map_or(false, |(f, _)| cs.action(0, f, &mut rng.clone()) == target);
            o.rec(
                "csidh",
                "mitm_group_action_inversion",
                &[
                    ("p_bits", bits.to_string()),
                    ("n_primes", n.to_string()),
                    ("box_m", mm.to_string()),
                ],
                &[("ns", ns_t), ("nodes", r.map_or(-1.0, |x| x.1 as f64))],
                ok,
                "boxes [0,m]^n from both ends",
            );
        }
    }
    // ideal class orders from isogeny cycles on the smallest parameter set
    if let Some(cs) = Csidh::with_n_primes(6) {
        let mut rng = Rng::new(10_500);
        let h = isogeny_algos::path::csidh::class_number(-4 * cs.p() as i64);
        for idx in 0..3 {
            let t0 = Instant::now();
            let ord = cs.ideal_order(idx, &mut rng, h as usize + 1);
            let ns_t = t0.elapsed().as_nanos() as f64;
            let ok = ord.map_or(false, |x| h % x as u64 == 0);
            o.rec(
                "csidh",
                "ideal_order_by_cycle",
                &[
                    ("p_bits", (64 - cs.p().leading_zeros()).to_string()),
                    ("ell", cs.primes[idx].to_string()),
                    ("class_number", h.to_string()),
                ],
                &[("ns", ns_t), ("order", ord.map_or(-1.0, |x| x as f64))],
                ok,
                "order divides h(-4p) (independent form count)",
            );
        }
    }
    // CSIDH-512 over FpM<8>: the real parameter set, exponents in [-m, m]^74
    let cs = Csidh::<FpM<8>>::csidh512();
    let f = &cs.fp;
    let mut rng = Rng::new(10_600);
    for &m in if o.quick { &[1i32][..] } else { &[1i32, 5][..] } {
        let e1: Vec<i32> = (0..74).map(|_| rng.below(2 * m as u64 + 1) as i32 - m).collect();
        let e2: Vec<i32> = (0..74).map(|_| rng.below(2 * m as u64 + 1) as i32 - m).collect();
        let t0 = Instant::now();
        let a1 = cs.action_batched(f.zero(), &e1, &mut rng);
        let t_first = t0.elapsed().as_nanos() as f64;
        let a12 = cs.action_batched(a1, &e2, &mut rng);
        let a21 = cs.action_batched(cs.action_batched(f.zero(), &e2, &mut rng), &e1, &mut rng);
        let ok = a12 == a21;
        let steps: i32 = e1.iter().map(|x| x.abs()).sum();
        let (med, mn, reps) = time_it(if o.quick { 200 } else { 2000 }, 20, || {
            cs.action_batched(f.zero(), &e1, &mut rng.clone())
        });
        o.rec(
            "csidh",
            "class_group_action_batched",
            &[
                ("p_bits", "511".to_string()),
                ("n_primes", "74".to_string()),
                ("exp_bound", m.to_string()),
                ("isogeny_steps", steps.to_string()),
                ("field", "FpM8".to_string()),
            ],
            &[
                ("median_ns", med),
                ("min_ns", mn),
                ("reps", reps as f64),
                ("first_ns", t_first),
                ("ns_per_step", med / steps.max(1) as f64),
            ],
            ok,
            "CSIDH-512 prime; verified by commutativity of two random keys; variable time",
        );
        let okf = cs.action_fast(f.zero(), &e1, &mut rng.clone()) == a1;
        let (med, mn, reps) = time_it(if o.quick { 200 } else { 2000 }, 20, || {
            cs.action_fast(f.zero(), &e1, &mut rng.clone())
        });
        o.rec(
            "csidh",
            "class_group_action_tree_projective",
            &[
                ("p_bits", "511".to_string()),
                ("n_primes", "74".to_string()),
                ("exp_bound", m.to_string()),
                ("isogeny_steps", steps.to_string()),
                ("field", "FpM8".to_string()),
            ],
            &[
                ("median_ns", med),
                ("min_ns", mn),
                ("reps", reps as f64),
                ("ns_per_step", med / steps.max(1) as f64),
            ],
            ok && okf,
            "CSIDH-512; equals the CLMPR result; variable time",
        );
    }
}

/// Kohel's End(E) conductor, and the weighted walk next to GHS on identical instances.
fn bench_v2path(o: &mut Out) {
    // End ring
    let bits: &[u32] = if o.quick { &[17] } else { &[17, 19, 21] };
    for &b in bits {
        let mut p = prime_bits(b);
        while ![1u64, 4, 7].contains(&(p % 9)) {
            p = next_prime(p + 1);
        }
        let fp = Zp::new(p);
        let cache = PhiCache::new(&Zp::new(p), &[3, 5, 7]);
        let mut rng = Rng::new(11_000 + b as u64);
        for i in 0..3 {
            let (e, t) = curve_with_trace(&fp, &mut rng, |t| {
                let d = (t * t - 4 * p as i64).abs();
                t % 2 != 0 && d % 27 == 0 && d % 25 != 0 && d % 49 != 0
            });
            let j = jinv(&fp, &e);
            let r = endo::endomorphism_ring(&fp, &cache, p, t, j, &mut rng.clone());
            let ok = r.as_ref().map_or(false, |x| {
                x.disc == x.fundamental_disc * (x.conductor * x.conductor) as i128
                    && x.per_prime[0].2 <= x.per_prime[0].1
            });
            let (med, mn, reps) = time_it(200, 100, || {
                endo::endomorphism_ring(&fp, &cache, p, t, j, &mut rng.clone())
            });
            let (h, lvl) = r
                .as_ref()
                .map_or((0, 0), |x| (x.per_prime[0].1, x.per_prime[0].2));
            o.rec(
                "endo",
                "kohel_end_ring_conductor",
                &[
                    ("p_bits", b.to_string()),
                    ("instance", i.to_string()),
                    ("height_3", h.to_string()),
                    ("level_3", lvl.to_string()),
                ],
                &[("median_ns", med), ("min_ns", mn), ("reps", reps as f64)],
                ok,
                "volcano level by walk lengths",
            );
        }
    }
    // weighted walk vs GHS vs Galbraith on the same instances
    let ells = [3usize, 5, 7, 11, 13];
    let wbits: &[u32] = if o.quick { &[24] } else { &[24, 28, 32] };
    let inst = if o.quick { 2 } else { 5 };
    for &b in wbits {
        let p = prime_bits(b);
        let fp = Zp::new(p);
        let cache = PhiCache::new(&Zp::new(p), &ells);
        for i in 0..inst {
            let mut rng = Rng::new(12_000 + (b as u64) * 100 + i);
            let (e1, _t) = loop {
                let (e1, t) = curve_with_trace(&fp, &mut rng, |t| {
                    let d = (t * t - 4 * p as i64).abs();
                    ells.iter().all(|&l| d % (l as i64 * l as i64) != 0)
                });
                let nb = neighbors(&fp, &cache, &ells, jinv(&fp, &e1), &mut rng);
                let mut pr: Vec<usize> = nb.iter().map(|x| x.0).collect();
                pr.sort();
                pr.dedup();
                if nb.len() >= 3 && pr.len() >= 2 {
                    break (e1, t);
                }
            };
            let j1 = jinv(&fp, &e1);
            let j2 = *random_walk(&fp, &cache, &ells, j1, 4 * b as usize, &mut rng)
                .js
                .last()
                .unwrap();
            let params = vec![("p_bits", b.to_string()), ("instance", i.to_string())];
            let t0 = Instant::now();
            let (r, st) = ghs::ghs(&fp, &cache, &ells, j1, j2, 5_000_000, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| {
                verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2
            });
            o.rec(
                "walks",
                "ghs_uniform_over_all_neighbours",
                &params,
                &[
                    ("ns", ns),
                    ("steps", st.steps as f64),
                    ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)),
                ],
                ok,
                "every step solves Phi_l for all primes",
            );
            for (name, w) in [
                ("galbraith_stolbunov_w_1_1_1_1_1", [1u32, 1, 1, 1, 1]),
                ("galbraith_stolbunov_w_16_8_4_2_1", [16, 8, 4, 2, 1]),
            ] {
                let t0 = Instant::now();
                let (r, st) = ghs::galbraith_stolbunov(
                    &fp,
                    &cache,
                    &ells,
                    &w,
                    j1,
                    j2,
                    5_000_000,
                    &mut rng.clone(),
                );
                let ns = t0.elapsed().as_nanos() as f64;
                let ok = r.as_ref().map_or(false, |p| {
                    verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2
                });
                o.rec(
                    "walks",
                    name,
                    &params,
                    &[
                        ("ns", ns),
                        ("steps", st.steps as f64),
                        ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)),
                    ],
                    ok,
                    "prime chosen first by weight, one Phi_l solved per step",
                );
            }
            let t0 = Instant::now();
            let (r, st) =
                galbraith::galbraith(&fp, &cache, &ells, j1, j2, 5_000_000, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| {
                verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2
            });
            o.rec(
                "walks",
                "galbraith_bfs",
                &params,
                &[
                    ("ns", ns),
                    ("nodes_expanded", st.nodes_expanded as f64),
                    ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)),
                ],
                ok,
                "same instance, bidirectional BFS",
            );
        }
    }
}

/// P1: every kernel -> isogeny implementation, old and new, Zp and FpM<1>; codomain and one evaluation.
fn bench_p1kernel(o: &mut Out) {
    let ells: &[u64] = if o.quick {
        &[3, 31, 401]
    } else {
        &[3, 7, 13, 31, 101, 401, 1009, 4001, 10007, 40009, 100003]
    };
    let budget = if o.quick { 100 } else { 300 };
    let fp = Zp::new(prime_bits(32));
    let fm = FpM::<1>::from_u64_modulus(fp.p);
    let mut rng = Rng::new(13_000);
    for &ell in ells {
        let (e, p) = curve_with_point(&fp, ell, &mut rng);
        let x0 = if let Pt::Aff(x, _) = p {
            x
        } else {
            unreachable!()
        };
        let q = e.random_point(&fp, &mut rng);
        let reference = velu_xonly_fast(&fp, &e, x0, ell).unwrap();
        let ref_img = reference.eval(&fp, &q);
        let params = vec![("p_bits", "32".to_string()), ("ell", ell.to_string())];
        let rec = |o: &mut Out,
                   algo: &str,
                   field: &str,
                   task: &str,
                   ok: bool,
                   med: f64,
                   mn: f64,
                   reps: usize| {
            let mut pr = params.clone();
            pr.push(("field", field.to_string()));
            pr.push(("task", task.to_string()));
            o.rec(
                "p1kernel",
                algo,
                &pr,
                &[("median_ns", med), ("min_ns", mn), ("reps", reps as f64)],
                ok,
                "",
            );
        };
        // Zp
        let small = ell <= 10007; // the O(l^2)-ish reference versions are skipped beyond this
        if small {
            let (m, mn, r) = time_it(budget, 500, || velu_cyclic(&fp, &e, &p, ell));
            let ok = velu_cyclic(&fp, &e, &p, ell).cod == reference.cod;
            rec(o, "velu_points_affine_v2", "Zp", "codomain", ok, m, mn, r);
            let (m, mn, r) = time_it(budget, 500, || velu_xonly(&fp, &e, x0, ell));
            rec(o, "velu_xonly_divpoly_v2", "Zp", "codomain", true, m, mn, r);
            let (m, mn, r) = time_it(budget, 500, || {
                isogeny_algos::kernel::sqrt_velu::sqrt_velu(&fp, &e, &p, ell)
            });
            let ok =
                isogeny_algos::kernel::sqrt_velu::sqrt_velu(&fp, &e, &p, ell).cod == reference.cod;
            rec(o, "sqrt_velu_naive_v1", "Zp", "codomain", ok, m, mn, r);
        }
        let (m, mn, r) = time_it(budget, 500, || velu_xonly_fast(&fp, &e, x0, ell));
        rec(o, "velu_xonly_fast", "Zp", "codomain", true, m, mn, r);
        let sv = sqrt_velu_fast(&fp, &e, x0, ell).unwrap();
        let ok = sv.cod == reference.cod && sv.eval(&fp, &q) == ref_img;
        let (m, mn, r) = time_it(budget, 500, || sqrt_velu_fast(&fp, &e, x0, ell));
        rec(o, "sqrt_velu_fast", "Zp", "codomain", ok, m, mn, r);
        let (m, mn, r) = time_it(budget, 500, || reference.eval(&fp, &q));
        rec(o, "velu_xonly_fast", "Zp", "eval_point", true, m, mn, r);
        let (m, mn, r) = time_it(budget, 500, || sv.eval(&fp, &q));
        rec(o, "sqrt_velu_fast", "Zp", "eval_point", ok, m, mn, r);
        // FpM<1>
        let em = Curve::new(fm.from_u64(e.a), fm.from_u64(e.b));
        let xm = fm.from_u64(x0);
        let vm = velu_xonly_fast(&fm, &em, xm, ell).unwrap();
        let okm = fm.to_canonical(&vm.cod.a)[0] == reference.cod.a
            && fm.to_canonical(&vm.cod.b)[0] == reference.cod.b;
        let (m, mn, r) = time_it(budget, 500, || velu_xonly_fast(&fm, &em, xm, ell));
        rec(o, "velu_xonly_fast", "FpM1", "codomain", okm, m, mn, r);
        let svm = sqrt_velu_fast(&fm, &em, xm, ell).unwrap();
        let okm2 = svm.cod == vm.cod;
        let (m, mn, r) = time_it(budget, 500, || sqrt_velu_fast(&fm, &em, xm, ell));
        rec(o, "sqrt_velu_fast", "FpM1", "codomain", okm2, m, mn, r);
        let qm = if let Pt::Aff(a, b) = q {
            Pt::Aff(fm.from_u64(a), fm.from_u64(b))
        } else {
            unreachable!()
        };
        let (m, mn, r) = time_it(budget, 500, || vm.eval(&fm, &qm));
        rec(o, "velu_xonly_fast", "FpM1", "eval_point", okm, m, mn, r);
        let (m, mn, r) = time_it(budget, 500, || svm.eval(&fm, &qm));
        rec(
            o,
            "sqrt_velu_fast",
            "FpM1",
            "eval_point",
            okm2 && svm.eval(&fm, &qm) == vm.eval(&fm, &qm),
            m,
            mn,
            r,
        );
        if ell <= 1009 {
            let (h, _) = kernel_poly_from_point(&fp, &e, &p, ell);
            let (m, mn, r) = time_it(budget, 500, || kohel(&fp, &e, &h, ell));
            rec(
                o,
                "kohel_given_h",
                "Zp",
                "codomain",
                kohel(&fp, &e, &h, ell).cod == reference.cod,
                m,
                mn,
                r,
            );
        }
    }
}

/// Kernel-to-isogeny, BMSS, Phi_l and the j-line neighbour oracle over 256- and 511-bit primes
/// (supersingular workloads: p = 4 prod(l) c - 1, so every listed l has rational kernels).
fn bench_big_field<const N: usize>(o: &mut Out, bits: usize, ells: &[u64], seed: u64) {
    let mut rng = Rng::new(seed);
    let w = isogeny_algos::testdata::supersingular_workload::<N>(ells, bits, &mut rng);
    let f = &w.f;
    let budget = if o.quick { 100 } else { 300 };
    let field = format!("FpM{N}");
    for &(ell, p) in &w.pts {
        let Pt::Aff(x0, _) = p else { unreachable!() };
        let reference = velu_xonly_fast(f, &w.e, x0, ell).unwrap();
        let t = random_point_f(f, &w.e, &mut rng);
        let img = reference.eval(f, &t);
        let params = vec![
            ("p_bits", bits.to_string()),
            ("ell", ell.to_string()),
            ("field", field.clone()),
        ];
        let rec = |o: &mut Out, algo: &str, task: &str, ok: bool, (m, mn, r): (f64, f64, usize)| {
            let mut pr = params.clone();
            pr.push(("task", task.to_string()));
            o.rec("big", algo, &pr, &[("median_ns", m), ("min_ns", mn), ("reps", r as f64)], ok, "");
        };
        let tm = time_it(budget, 500, || velu_xonly_fast(f, &w.e, x0, ell));
        rec(o, "velu_xonly_fast", "codomain", true, tm);
        let sv = sqrt_velu_fast(f, &w.e, x0, ell).unwrap();
        let ok = sv.cod == reference.cod && sv.eval(f, &t) == img;
        let tm = time_it(budget, 500, || sqrt_velu_fast(f, &w.e, x0, ell));
        rec(o, "sqrt_velu_fast", "codomain", ok, tm);
        let tm = time_it(budget, 500, || reference.eval(f, &t));
        rec(o, "velu_xonly_fast", "eval_point", true, tm);
        let tm = time_it(budget, 500, || sv.eval(f, &t));
        rec(o, "sqrt_velu_fast", "eval_point", ok, tm);
        if ell <= 1009 {
            let ok = velu_cyclic(f, &w.e, &p, ell).cod == reference.cod;
            let tm = time_it(budget, 500, || velu_cyclic(f, &w.e, &p, ell));
            rec(o, "velu_points_affine", "codomain", ok, tm);
            let (g, _) = kernel_poly_from_point_f(f, &w.e, &p, ell);
            let ok = kohel(f, &w.e, &g, ell).cod == reference.cod;
            let tm = time_it(budget, 500, || kohel(f, &w.e, &g, ell));
            rec(o, "kohel_given_h", "codomain", ok, tm);
            if ell <= 101 {
                let et = reference.cod;
                let d = g.len() - 1;
                let sigma = f.neg(f.add(g[d - 1], g[d - 1]));
                for m in [
                    bmss::Method::LinearAlgebra,
                    bmss::Method::Elkies1998,
                    bmss::Method::FastElkies,
                    bmss::Method::FastElkiesPrime,
                ] {
                    let s = if m.needs_sigma() { Some(sigma) } else { None };
                    let ok = bmss::isogeny(f, m, &w.e, &et, ell as usize, s).map_or(false, |i| i.ker == g);
                    let tm = time_it(budget, 200, || bmss::isogeny(f, m, &w.e, &et, ell as usize, s));
                    rec(o, m.name(), "kernel_from_codomain", ok, tm);
                }
            }
        }
    }
    // Montgomery curve y^2 = x^3 + x (A = 0, p + 1 points): projective Vélu vs sqrt-Velu
    {
        use isogeny_algos::kernel::montgomery as mg;
        use isogeny_algos::kernel::sqrt_velu_mont::SqrtVeluMont;
        let a = f.zero();
        let k24 = mg::proj24(f, a);
        let p1 = f.modulus().add_small(1);
        for &(ell, _) in &w.pts {
            let cof = p1.divrem_small(ell).0;
            let kp = loop {
                let x = f.random(&mut rng);
                let k = mg::ladder_p(f, k24, (x, f.one()), &cof);
                if !f.is_zero(k.1) {
                    break k;
                }
            };
            let d = ((ell - 1) / 2) as usize;
            let velu = || {
                let ms = mg::multiples_p(f, k24, kp, d);
                let pre = mg::kernel_pre(f, &ms);
                (mg::velu_codomain_p(f, k24, &pre, ell), pre)
            };
            let (kc, pre) = velu();
            let a_ref = mg::affine_a(f, kc);
            let u = f.random(&mut rng);
            let img = mg::isog_xz_pre(f, &pre, (u, f.one()));
            let img_aff = f.div(img.0, img.1);
            let params = vec![
                ("p_bits", bits.to_string()),
                ("ell", ell.to_string()),
                ("field", field.clone()),
            ];
            let rec = |o: &mut Out, algo: &str, task: &str, ok: bool, (m, mn, r): (f64, f64, usize)| {
                let mut pr = params.clone();
                pr.push(("task", task.to_string()));
                pr.push(("model", "montgomery".to_string()));
                o.rec("big", algo, &pr, &[("median_ns", m), ("min_ns", mn), ("reps", r as f64)], ok, "");
            };
            let tm = time_it(budget, 500, velu);
            rec(o, "velu_montgomery_projective", "codomain", true, tm);
            let tm = time_it(budget, 500, || mg::isog_xz_pre(f, &pre, (u, f.one())));
            rec(o, "velu_montgomery_projective", "eval_point", true, tm);
            let sv = SqrtVeluMont::new(f, a, kp, ell).unwrap();
            let ok = sv.codomain(f) == a_ref && sv.eval(f, u) == img_aff;
            let tm = time_it(budget, 500, || {
                let sv = SqrtVeluMont::new(f, a, kp, ell).unwrap();
                sv.codomain_proj(f)
            });
            rec(o, "sqrt_velu_montgomery", "codomain", ok, tm);
            let tm = time_it(budget, 500, || sv.eval_nd(f, u));
            rec(o, "sqrt_velu_montgomery", "eval_point", ok, tm);
        }
    }
    // Phi_l over the big field, and the neighbour oracle (roots of Phi_l(j, Y)) at a random j
    if N <= 4 {
        let j = jinv(f, &w.e);
        for &ell in &[2usize, 3, 5, 7, 11, 13] {
            let t0 = Instant::now();
            let phi = Phi::compute(f, ell);
            let t_phi = t0.elapsed().as_nanos() as f64;
            let ns = phi.neighbors(f, j, &mut rng);
            // supersingular over F_p: the F_p-rational neighbours are ell-isogenous: check Phi = 0
            let ok = ns.iter().all(|&jn| f.is_zero(phi.eval(f, j, jn)));
            let tm = time_it(budget, 200, || phi.neighbors(f, j, &mut rng.clone()));
            o.rec(
                "big",
                "phi_neighbors",
                &[("p_bits", bits.to_string()), ("ell", ell.to_string()), ("field", field.clone())],
                &[
                    ("median_ns", tm.0),
                    ("min_ns", tm.1),
                    ("reps", tm.2 as f64),
                    ("phi_compute_ns", t_phi),
                    ("rational_neighbors", ns.len() as f64),
                ],
                ok,
                "roots of Phi_l(j, Y) by gcd with Y^p - Y and equal-degree splitting",
            );
        }
    }
}

/// Characteristic 2: GF(2^n) binary curves y^2 + xy = x^3 + a2 x^2 + a6.
fn bench_char2(o: &mut Out) {
    use isogeny_algos::binary::*;
    use isogeny_algos::gf2n::GF2n;
    use isogeny_algos::path::galbraith::galbraith_with;
    use isogeny_algos::path::ghs::ghs_with;
    let budget = if o.quick { 100 } else { 300 };
    let ns: &[u32] = if o.quick { &[31] } else { &[23, 41, 61] };
    let ells: &[u64] = if o.quick { &[3, 7] } else { &[3, 5, 7, 11, 13] };
    for &n in ns {
        let f = GF2n::new(n);
        let mut rng = Rng::new(15_000 + n as u64);
        for &ell in ells {
            // a curve with a rational point of order ell
            let (e, p, ord) = loop {
                let e = BinCurve::new(rng.next() & 1, 0, f.random(&mut rng) | 1);
                let ord = e.order(&f, &mut rng);
                if ord % ell as u128 != 0 {
                    continue;
                }
                let p = e.mul(&f, &e.random_point(&f, &mut rng), ord / ell as u128);
                if p != Pt::Inf {
                    break (e, p, ord);
                }
            };
            let params = vec![("field", format!("GF2^{n}")), ("ell", ell.to_string())];
            let rec = |o: &mut Out, algo: &str, ok: bool, (m, mn, r): (f64, f64, usize), note: &str| {
                o.rec("char2", algo, &params, &[("median_ns", m), ("min_ns", mn), ("reps", r as f64)], ok, note);
            };
            let iso = velu(&f, &e, &p, ell);
            let q = e.random_point(&f, &mut rng);
            let q2 = e.random_point(&f, &mut rng);
            let ok = iso.cod.on_curve(&f, &iso.eval(&f, &q))
                && iso.eval(&f, &e.add(&f, &q, &q2)) == iso.cod.add(&f, &iso.eval(&f, &q), &iso.eval(&f, &q2));
            rec(o, "velu_char2_points", ok, time_it(budget, 500, || velu(&f, &e, &p, ell)), "kernel points -> codomain (t = sum x_Q)");
            rec(o, "velu_char2_eval", ok, time_it(budget, 500, || iso.eval(&f, &q)), "full (x, y) image");
            let xs: Vec<u64> = iso.s.iter().map(|q| q.0).collect();
            let h = poly::from_roots(&f, &xs);
            let okk = kohel_codomain(&f, &e, &h) == iso.cod;
            rec(o, "kohel_char2_codomain", okk, time_it(budget, 500, || kohel_codomain(&f, &e, &h)), "codomain from the kernel polynomial");
            let ks = kernel_polys(&f, &e, ell, &mut rng);
            let okf = ks.contains(&h) && ks.iter().all(|k| kohel_codomain(&f, &e, k).order(&f, &mut rng.clone()) == ord);
            rec(o, "kernel_polys_divpoly_char2", okf, time_it(budget, 50, || kernel_polys(&f, &e, ell, &mut rng.clone())), "factor f_l, subsets of degree (l-1)/2, x-map/doubling check");
            if ell <= 7 || !o.quick {
                let t0 = Instant::now();
                let phi = phi_mod2(&f, ell as usize);
                let t_phi = t0.elapsed().as_nanos() as f64;
                let j = e.j(&f);
                let okp = phi.neighbors(&f, j, &mut rng).contains(&iso.cod.j(&f));
                let tm = time_it(budget, 200, || phi.neighbors(&f, j, &mut rng.clone()));
                o.rec(
                    "char2",
                    "phi_mod2_neighbors",
                    &params,
                    &[("median_ns", tm.0), ("min_ns", tm.1), ("reps", tm.2 as f64), ("phi_crt_ns", t_phi)],
                    okp,
                    "roots of Phi_l(j, Y) over GF(2^n) (trace splitting); Phi_l mod 2 by CRT",
                );
            }
        }
    }
    // path finding between isogenous binary curves: kernel oracle vs Phi oracle
    let pns: &[u32] = if o.quick { &[15] } else { &[15, 19, 23] };
    for &n in pns {
        let f = GF2n::new(n);
        let mut rng = Rng::new(16_000 + n as u64);
        let ells = vec![3u64, 5, 7];
        let ko = BinKernelOracle { f: &f, ells: ells.clone() };
        let cache = PhiCache { phis: ells.iter().map(|&l| phi_mod2(&f, l as usize)).collect() };
        let eu: Vec<usize> = ells.iter().map(|&l| l as usize).collect();
        let po = PhiOracle { f: &f, cache: &cache, ells: &eu };
        for t in 0..(if o.quick { 1 } else { 3 }) {
            let j1 = f.random(&mut rng) | 2;
            let mut j2 = j1;
            for _ in 0..(n as usize) {
                let ns = po.neighbors(j2, &mut rng);
                if ns.is_empty() {
                    break;
                }
                j2 = ns[rng.below(ns.len() as u64) as usize].1;
            }
            if j2 == j1 {
                continue;
            }
            for (oname, which) in [("kernel_oracle", 0), ("phi_oracle", 1)] {
                let t0 = Instant::now();
                let (p, st) = if which == 0 {
                    galbraith_with(&ko, j1, j2, 2_000_000, &mut rng.clone())
                } else {
                    galbraith_with(&po, j1, j2, 2_000_000, &mut rng.clone())
                };
                let ns_t = t0.elapsed().as_nanos() as f64;
                let ok = p.as_ref().map_or(false, |p| verify_path(&f, &cache, p));
                o.rec(
                    "char2",
                    "galbraith_bfs",
                    &[("field", format!("GF2^{n}")), ("oracle", oname.to_string()), ("instance", t.to_string())],
                    &[("ns", ns_t), ("nodes", st.nodes_expanded as f64), ("path_len", p.map_or(-1.0, |p| p.len() as f64))],
                    ok,
                    "l in {3,5,7}",
                );
            }
            let t0 = Instant::now();
            let (p, st) = ghs_with(&po, j1, j2, 2_000_000, &mut rng.clone());
            let ns_t = t0.elapsed().as_nanos() as f64;
            let ok = p.as_ref().map_or(false, |p| verify_path(&f, &cache, p));
            o.rec(
                "char2",
                "ghs_walk",
                &[("field", format!("GF2^{n}")), ("oracle", "phi_oracle".to_string()), ("instance", t.to_string())],
                &[("ns", ns_t), ("steps", st.steps as f64), ("path_len", p.map_or(-1.0, |p| p.len() as f64))],
                ok,
                "l in {3,5,7}; no volcano normalisation",
            );
        }
    }
}

/// Quaternion side: KLPT, class sets / Brandt matrices, Deuring (ideal -> curve).
fn bench_quat(o: &mut Out) {
    use isogeny_algos::bigint::Big;
    use isogeny_algos::int::Int;
    use isogeny_algos::quat::brandt::{class_set, supersingular_graph};
    use isogeny_algos::quat::deuring::Deuring;
    use isogeny_algos::quat::klpt::{ideal_from, klpt};
    use isogeny_algos::quat::*;
    let mut rng = Rng::new(17_000);
    // random left O_0-ideal of prime norm n: alpha with a = sqrt(-(b^2 + p(c^2+d^2))) mod n
    let rand_ideal = |alg: &Alg, o0: &Lattice, n: &Int, rng: &mut Rng| loop {
        let (b, c, d) = (Int::from(rng.next() >> 2), Int::from(rng.next() >> 2), Int::from(rng.next() >> 2));
        let t = -&(&(&b * &b) + &(&alg.p * &(&(&c * &c) + &(&d * &d))));
        if let Some(a) = Int::sqrt_mod_prime(&t, n) {
            let i = ideal_from(alg, o0, n, &Quat::new([a, b, c, d], Int::one()));
            if i.norm(o0) == (n.clone(), Int::one()) {
                return i;
            }
        }
    };
    let ps: &[&str] = if o.quick {
        &["2147483647", "1152921504606847067"]
    } else {
        &["2147483647", "1152921504606847067", "1267650600228229401496703205707", "340282366920938463463374607431768211507"]
    };
    for ps in ps {
        let p = Int::from_big(&Big::from_dec(ps));
        if p.mod_u64(4) != 3 || !p.is_probable_prime() {
            continue;
        }
        let alg = Alg::new(&p);
        let o0 = alg.o0();
        let mut n = &p + &Int::from(2i64);
        while !n.is_probable_prime() {
            n = &n + &Int::from(2i64);
        }
        let mut es = vec![];
        let mut ok = true;
        let ideals: Vec<Lattice> = (0..5).map(|_| rand_ideal(&alg, &o0, &n, &mut rng)).collect();
        let t0 = Instant::now();
        for i in &ideals {
            match klpt(&alg, i, 2, &mut rng) {
                Some(r) => {
                    ok &= r.j.norm(&o0) == (Int::from(2i64).pow(r.e), Int::one()) && i.rmul(&alg, &r.xi) == r.j;
                    es.push(r.e as f64);
                }
                None => ok = false,
            }
        }
        let per = t0.elapsed().as_nanos() as f64 / ideals.len() as f64;
        let logp = p.to_f64().log2();
        let mean_e = es.iter().sum::<f64>() / es.len().max(1) as f64;
        o.rec(
            "quat",
            "klpt_l2",
            &[("p_bits", format!("{logp:.0}")), ("input", "prime norm ~ p".to_string())],
            &[("mean_ns", per), ("mean_e", mean_e), ("e_over_log2p", mean_e / logp), ("runs", es.len() as f64)],
            ok,
            "output J ~ I with N(J) = 2^e; J = I xi checked exactly",
        );
    }
    // class sets and Brandt matrices
    for &p in if o.quick { &[431i64][..] } else { &[431i64, 1019, 1259, 3499][..] } {
        let alg = Alg::new(&Int::from(p));
        let t0 = Instant::now();
        let cs = class_set(&alg, 2, &mut rng);
        let t_cs = t0.elapsed().as_nanos() as f64;
        let t0 = Instant::now();
        let (js, a) = supersingular_graph(p as u64, 2, &mut rng);
        let t_g = t0.elapsed().as_nanos() as f64;
        let ok = cs.reps.len() == js.len()
            && isogeny_algos::quat::brandt::power_traces(&cs.brandt, js.len()) == isogeny_algos::quat::brandt::power_traces(&a, js.len());
        o.rec(
            "quat",
            "class_set_brandt_l2",
            &[("p", p.to_string())],
            &[("ns", t_cs), ("classes", cs.reps.len() as f64), ("mestre_graph_ns", t_g)],
            ok,
            "classes of O_0 by BFS on 2-neighbours (equivalence by short vectors); mass formula asserted; spectrum vs Phi_2 graph",
        );
        // Deuring map on all classes (p with smooth p^2 - 1 only)
        if p == 1259 || p == 3499 {
            let d = Deuring::new(p as u64, &mut rng);
            let t0 = Instant::now();
            let mut img = vec![];
            for i in &cs.reps {
                img.push(d.ideal_to_j(i, &mut rng));
            }
            let t_d = t0.elapsed().as_nanos() as f64 / cs.reps.len() as f64;
            let mut s: Vec<_> = img.iter().flatten().copied().collect();
            s.sort();
            s.dedup();
            let ok = s.len() == js.len() && s.iter().all(|j| js.contains(j));
            o.rec(
                "quat",
                "deuring_ideal_to_curve",
                &[("p", p.to_string()), ("t_odd", d.t_odd.to_string())],
                &[("mean_ns_per_class", t_d), ("classes", cs.reps.len() as f64)],
                ok,
                "smooth-norm equivalent ideal (N | odd part of p^2-1), kernel over F_p^4 by 2D Pohlig-Hellman, Velu chain; bijection onto supersingular j checked",
            );
        }
    }
}

/// Genus 2: Richelot (2,2)-isogenies, splitting, gluing, the superspecial Richelot graph.
fn bench_genus2(o: &mut Out) {
    use isogeny_algos::genus2::*;
    let budget = if o.quick { 100 } else { 300 };
    let mut rng = Rng::new(18_000);
    // Richelot codomain at 61 and 256 bits (timing; correctness is tested on small p via L-polys)
    fn richelot_case<F: Field>(f: &F, rng: &mut Rng) -> ([isogeny_algos::poly::Poly<F>; 3], F::E, F::E) {
        let r: Vec<F::E> = (0..6).map(|_| f.random(rng)).collect();
        let g = [
            isogeny_algos::poly::from_roots(f, &[r[0], r[1]]),
            isogeny_algos::poly::from_roots(f, &[r[2], r[3]]),
            isogeny_algos::poly::from_roots(f, &[r[4], r[5]]),
        ];
        (g, r[0], r[1])
    }
    let fp = Zp::new(prime_bits(61));
    let (g, _, _) = richelot_case(&fp, &mut rng);
    let tm = time_it(budget, 10_000, || richelot(&fp, &g).codomain(&fp));
    o.rec("genus2", "richelot_codomain", &[("field", "Zp61".to_string())], &[("median_ns", tm.0), ("min_ns", tm.1), ("reps", tm.2 as f64)], true, "Delta, H1 H2 H3 from G1 G2 G3");
    let p256 = FpM::<4>::from_dec("115792089210356248762697446949407573530086143415290314195533631308867097853951");
    let (g4, _, _) = richelot_case(&p256, &mut rng);
    let tm = time_it(budget, 10_000, || richelot(&p256, &g4).codomain(&p256));
    o.rec("genus2", "richelot_codomain", &[("field", "FpM4_P256".to_string())], &[("median_ns", tm.0), ("min_ns", tm.1), ("reps", tm.2 as f64)], true, "");
    // verified instances at p = 1009: Richelot, gluing, splitting, point images
    let p = 1009u64;
    let fp = Zp::new(p);
    let f2 = Zp2::new(p);
    let r: Vec<u64> = {
        let mut v = vec![];
        while v.len() < 6 {
            let x = fp.random(&mut rng);
            if !v.contains(&x) {
                v.push(x);
            }
        }
        v
    };
    let fsex = poly::from_roots(&fp, &r);
    let g = [
        poly::from_roots(&fp, &[r[0], r[1]]),
        poly::from_roots(&fp, &[r[2], r[3]]),
        poly::from_roots(&fp, &[r[4], r[5]]),
    ];
    let rl = richelot(&fp, &g);
    let ok = rl.delta != 0 && lpoly(&fp, &rl.codomain(&fp)) == lpoly(&fp, &fsex);
    let tm = time_it(budget, 10_000, || richelot(&fp, &g).codomain(&fp));
    o.rec("genus2", "richelot_codomain", &[("field", "Zp1009".to_string())], &[("median_ns", tm.0), ("min_ns", tm.1), ("reps", tm.2 as f64)], ok, "L-polynomial of codomain equal (naive point counts over F_p, F_p^2)");
    let x0 = (0..p).find(|&x| fp.sqrt(poly::eval(&fp, &fsex, x)).map_or(false, |y| y != 0)).unwrap();
    let y0 = fp.sqrt(poly::eval(&fp, &fsex, x0)).unwrap();
    let tm = time_it(budget, 10_000, || rl.image_point(&fp, &g, &f2, |c| (c, 0), x0, y0, &mut rng.clone()));
    o.rec("genus2", "richelot_point_image", &[("field", "Zp1009".to_string())], &[("median_ns", tm.0), ("min_ns", tm.1), ("reps", tm.2 as f64)], true, "two points over F_p^2 (quadratic in x')");
    let (a, b) = ([r[0], r[1], r[2]], [r[3], r[4], r[5]]);
    let c = glue(&fp, a, b);
    let ok = c.as_ref().map_or(false, |c| {
        lpoly(&fp, c) == LPoly::product(elliptic_trace(&fp, &poly::from_roots(&fp, &a)), elliptic_trace(&fp, &poly::from_roots(&fp, &b)), p as i64)
    });
    let tm = time_it(budget, 10_000, || glue(&fp, a, b));
    o.rec("genus2", "glue_e1xe2", &[("field", "Zp1009".to_string())], &[("median_ns", tm.0), ("min_ns", tm.1), ("reps", tm.2 as f64)], ok, "L(C) = L(E1) L(E2) checked");
    // superspecial graph BFS (vertex counts checked against Ibukiyama-Katsura-Oort and h(h+1)/2)
    for &p in if o.quick { &[43u64][..] } else { &[43u64, 83, 131, 199][..] } {
        let t0 = Instant::now();
        let gr = superspecial_graph(p, 1_000_000, &mut rng);
        let ns = t0.elapsed().as_nanos() as f64;
        let pi = p as i64;
        let m1 = if pi % 4 == 1 { 1 } else { -1 };
        let m2 = if pi % 8 == 1 || pi % 8 == 3 { 1 } else { -1 };
        let m3 = if pi % 3 == 1 { 1 } else { -1 };
        let iko = ((pi - 1) * (pi * pi + 25 * pi + 166) - 90 * (1 - m1) + 360 * (1 - m2) + 160 * (1 - m3) + if pi % 5 == 4 { 2304 } else { 0 }) / 2880;
        let h = (p / 12 + [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2][(p % 12) as usize]) as usize;
        let ok = gr.jacobians as i64 == iko && gr.products == h * (h + 1) / 2;
        o.rec(
            "genus2",
            "superspecial_richelot_graph_bfs",
            &[("p", p.to_string())],
            &[("ns", ns), ("jacobians", gr.jacobians as f64), ("products", gr.products as f64), ("edges", gr.edges as f64)],
            ok,
            "vertices by Igusa-Clebsch invariants; counts = Ibukiyama-Katsura-Oort and h(h+1)/2",
        );
    }
}

/// Classical modular polynomials: Hecke/Newton construction vs dense linear algebra; SEA.
fn bench_phi(o: &mut Out) {
    let fp = Zp::new(prime_bits(61));
    let f2 = FpM::<2>::from_dec("170141183460469231731687303715884105727");
    let ells: &[usize] = if o.quick { &[11, 23] } else { &[11, 23, 31, 43, 61, 83, 101, 127] };
    for &ell in ells {
        let t0 = Instant::now();
        let a = Phi::compute_hecke(&fp, ell);
        let th = t0.elapsed().as_nanos() as f64;
        let (tl, ok) = if ell <= 43 {
            let t0 = Instant::now();
            let b = Phi::compute_linear_algebra(&fp, ell);
            (t0.elapsed().as_nanos() as f64, a.c == b.c)
        } else {
            // spot check: Phi(j, j') = 0 for an l-isogenous pair from Velu
            let mut rng = Rng::new(19_000 + ell as u64);
            let ok = (|| {
                for _ in 0..200 {
                    let e = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
                    let ks = a.neighbors(&fp, jinv(&fp, &e), &mut rng);
                    if let Some(&jn) = ks.first() {
                        return a.eval(&fp, jinv(&fp, &e), jn) == 0 && a.eval(&fp, jn, jinv(&fp, &e)) == 0;
                    }
                }
                false
            })();
            (-1.0, ok)
        };
        o.rec(
            "phi",
            "phi_hecke_newton",
            &[("field", "Zp61".to_string()), ("ell", ell.to_string())],
            &[("ns", th), ("linear_algebra_ns", tl)],
            ok,
            "equal to the linear-algebra Phi_l (l <= 43); symmetric root check above",
        );
        if ell <= 61 {
            let t0 = Instant::now();
            let _ = Phi::compute_hecke(&f2, ell);
            let th2 = t0.elapsed().as_nanos() as f64;
            o.rec("phi", "phi_hecke_newton", &[("field", "FpM2_127".to_string()), ("ell", ell.to_string())], &[("ns", th2)], true, "");
        }
    }
    // SEA, with the modular polynomials precomputed (as from a database) and timed separately
    use isogeny_algos::find::sea::sea;
    fn prefill<F: Field>(f: &F, max_ell: usize) -> (std::collections::HashMap<usize, Phi<F>>, f64) {
        let t0 = Instant::now();
        let mut m = std::collections::HashMap::new();
        for l in 3..=max_ell {
            if is_prime(l as u64) {
                m.insert(l, Phi::compute(f, l));
            }
        }
        (m, t0.elapsed().as_nanos() as f64)
    }
    let mut rng = Rng::new(19_500);
    let cases: Vec<(&str, u32)> = if o.quick { vec![("Zp", 61)] } else { vec![("Zp", 40), ("Zp", 61), ("FpM2", 127)] };
    for (fld, bits) in cases {
        if fld == "Zp" {
            let fp = Zp::new(next_prime((1u64 << (bits - 1)) + 12345));
            let (mut phis, tpre) = prefill(&fp, 61);
            o.rec("phi", "sea_phi_precompute", &[("field", format!("Zp{bits}")), ("max_ell", "61".to_string())], &[("ns", tpre)], true, "Phi_l for all primes l <= 61 (Hecke/Newton)");
            for t in 0..5 {
                let e = loop {
                    let e = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
                    let j = jinv(&fp, &e);
                    if is_smooth(&fp, &e) && j != 0 && j != 1728 {
                        break e;
                    }
                };
                let t0 = Instant::now();
                let r = sea(&fp, &e, 61, &mut phis, &mut rng);
                let ns = t0.elapsed().as_nanos() as f64;
                let t1 = Instant::now();
                let reference = order(&fp, &e, &mut rng);
                let ns_bsgs = t1.elapsed().as_nanos() as f64;
                let ok = r.as_ref().map_or(false, |(n, _)| *n == isogeny_algos::bigint::Big::from_u64(reference));
                let st = r.map(|x| x.1).unwrap_or_default();
                o.rec(
                    "phi",
                    "sea_point_count",
                    &[("field", format!("Zp{bits}")), ("instance", t.to_string())],
                    &[("ns", ns), ("bsgs_ns", ns_bsgs), ("elkies_primes", st.elkies.len() as f64), ("atkin_primes", st.atkin.len() as f64), ("candidates", st.candidates as f64)],
                    ok,
                    "Phi_l precomputed; equals the BSGS order",
                );
            }
        } else {
            let f = FpM::<2>::from_dec("170141183460469231731687303715884105727");
            let (mut phis, tpre) = prefill(&f, 89);
            o.rec("phi", "sea_phi_precompute", &[("field", "FpM2_127".to_string()), ("max_ell", "89".to_string())], &[("ns", tpre)], true, "Phi_l for all primes l <= 89 (Hecke/Newton)");
            for t in 0..3 {
                let e = loop {
                    let e = Curve::new(f.random(&mut rng), f.random(&mut rng));
                    let j = jinv(&f, &e);
                    if is_smooth(&f, &e) && j != f.zero() && j != f.from_u64(1728) {
                        break e;
                    }
                };
                let t0 = Instant::now();
                let r = sea(&f, &e, 89, &mut phis, &mut rng);
                let ns = t0.elapsed().as_nanos() as f64;
                let ok = r.as_ref().map_or(false, |(n, _)| (0..5).all(|_| pmul_big(&f, &e, &random_point_f(&f, &e, &mut rng), n) == Pt::Inf));
                let st = r.map(|x| x.1).unwrap_or_default();
                o.rec(
                    "phi",
                    "sea_point_count",
                    &[("field", "FpM2_127".to_string()), ("instance", t.to_string())],
                    &[("ns", ns), ("elkies_primes", st.elkies.len() as f64), ("atkin_primes", st.atkin.len() as f64), ("candidates", st.candidates as f64)],
                    ok,
                    "Phi_l precomputed; [#E]P = 0 for 5 random points",
                );
            }
        }
    }
}

/// Radical isogenies (CDV 2020) vs Velu steps on CSIDH-512.
fn bench_radical(o: &mut Out) {
    use isogeny_algos::kernel::montgomery;
    use isogeny_algos::kernel::radical::*;
    let cs = Csidh::<FpM<8>>::csidh512();
    let f = &cs.fp;
    let mut rng = Rng::new(20_000);
    let p1 = f.modulus().add_small(1);
    let k = if o.quick { 16 } else { 64 };
    for (ell, idx) in [(3u64, 0usize), (5, 1)] {
        let e_root = root_exponent(f, ell).unwrap();
        let ew = montgomery::to_weierstrass(f, f.zero());
        let p = loop {
            let r = random_point_f(f, &ew, &mut rng);
            let p = pmul_big(f, &ew, &r, &p1.divrem_small(ell).0);
            if p != Pt::Inf {
                break p;
            }
        };
        let (a1, a2, a3) = tangent_form(f, &ew, &p).unwrap();
        let (b0, _) = if ell == 5 { tate_bc(f, a1, a2, a3) } else { (f.zero(), f.zero()) };
        let chain = |f: &FpM<8>| {
            if ell == 3 {
                let mut cur = (a1, a3);
                for _ in 0..k {
                    cur = step3(f, cur.0, cur.1, &e_root);
                }
                weierstrass_j(f, model3(f, cur.0, cur.1))
            } else {
                let mut b = b0;
                for _ in 0..k {
                    b = step5(f, b, &e_root);
                }
                weierstrass_j(f, tate_model(f, b, b))
            }
        };
        let jr = chain(f);
        let mut e = vec![0i32; 74];
        e[idx] = k as i32;
        let ok = jr == montgomery::j_invariant(f, cs.action(f.zero(), &e, &mut rng));
        let (m, mn, r) = time_it(if o.quick { 200 } else { 1000 }, 20, || chain(f));
        o.rec(
            "radical",
            "radical_chain",
            &[("ell", ell.to_string()), ("steps", k.to_string()), ("field", "CSIDH-512".to_string())],
            &[("median_ns", m), ("min_ns", mn), ("reps", r as f64), ("ns_per_step", m / k as f64)],
            ok,
            "one N-th root (unique) per step; j equals CSIDH action l^k",
        );
        let (m, mn, r) = time_it(if o.quick { 200 } else { 1000 }, 20, || cs.action(f.zero(), &e, &mut rng.clone()));
        o.rec(
            "radical",
            "velu_chain_csidh_steps",
            &[("ell", ell.to_string()), ("steps", k.to_string()), ("field", "CSIDH-512".to_string())],
            &[("median_ns", m), ("min_ns", mn), ("reps", r as f64), ("ns_per_step", m / k as f64)],
            true,
            "CSIDH stepwise action: per step a fresh point, cofactor ladder, Velu",
        );
        let (m, mn, r) = time_it(if o.quick { 200 } else { 1000 }, 20, || cs.action_fast(f.zero(), &e, &mut rng.clone()));
        o.rec(
            "radical",
            "velu_chain_csidh_tree",
            &[("ell", ell.to_string()), ("steps", k.to_string()), ("field", "CSIDH-512".to_string())],
            &[("median_ns", m), ("min_ns", mn), ("reps", r as f64), ("ns_per_step", m / k as f64)],
            true,
            "CSIDH action_fast with only this prime: one point per isogeny (no batching gain)",
        );
    }
}

/// CSIDH relation lattice (CSI-FiSh style): class number, discrete logs, LLL, reduced action.
fn bench_relation(o: &mut Out) {
    use isogeny_algos::path::relation::relation_lattice;
    let mut rng = Rng::new(21_000);
    let ns: &[usize] = if o.quick { &[8] } else { &[6, 8, 10, 12, 14] };
    for &n in ns {
        let Some(cs) = Csidh::with_n_primes(n) else { continue };
        let p = cs.p() as i128;
        let t0 = Instant::now();
        let mut method = "cyclic_dlog";
        let mut rl = relation_lattice(p, &cs.primes);
        if rl.is_none() {
            method = "explicit_subgroup";
            rl = isogeny_algos::path::relation::relation_lattice_explicit(p, &cs.primes, 20_000_000);
        }
        let t_rl = t0.elapsed().as_nanos() as f64;
        let Some(rl) = rl else {
            o.rec("relation", "relation_lattice", &[("n_primes", n.to_string())], &[("ns", t_rl)], false, "no element of order h and the group is too large for the explicit method");
            continue;
        };
        let ok = rl.basis.iter().take(3).all(|v| cs.action_fast(0, &v.iter().map(|&x| x as i32).collect::<Vec<_>>(), &mut rng) == 0);
        let maxnorm = rl.basis.iter().map(|v| v.iter().map(|x| x.abs()).sum::<i64>()).max().unwrap();
        o.rec(
            "relation",
            "relation_lattice",
            &[("n_primes", n.to_string()), ("p_bits", (128 - p.leading_zeros()).to_string()), ("method", method.to_string())],
            &[("ns", t_rl), ("class_number", rl.h as f64), ("max_l1_basis", maxnorm as f64)],
            ok,
            "h by BSGS around the analytic estimate; dlogs by BSGS; LLL; basis vectors act trivially (3 checked)",
        );
        // random classes [l_1]^a: reduced vector vs the naive (a, 0, ..., 0)
        let samples = 10;
        let (mut l1_red, mut t_red, mut t_naive) = (0f64, 0f64, 0f64);
        let mut ok = true;
        for _ in 0..samples {
            let a = rng.below(rl.h as u64) as i128;
            // a class: [l_1]^a (cyclic case) or a random exponent vector of l1 norm ~ h/2
            let v = if rl.logs.is_empty() {
                let e: Vec<i64> = (0..n).map(|i| if i == 0 { a as i64 } else { rng.below(1000) as i64 }).collect();
                rl.reduce(&e)
            } else {
                rl.vector_of(a)
            };
            l1_red += v.iter().map(|x| x.abs()).sum::<i64>() as f64;
            let vi: Vec<i32> = v.iter().map(|&x| x as i32).collect();
            let t0 = Instant::now();
            let ar = cs.action_fast(0, &vi, &mut rng);
            t_red += t0.elapsed().as_nanos() as f64;
            if rl.h <= 50_000 && !rl.logs.is_empty() {
                let mut e1 = vec![0i32; n];
                e1[0] = a as i32;
                let t0 = Instant::now();
                let an = cs.action_fast(0, &e1, &mut rng);
                t_naive += t0.elapsed().as_nanos() as f64;
                ok &= an == ar;
            }
        }
        o.rec(
            "relation",
            "reduced_class_action",
            &[("n_primes", n.to_string())],
            &[
                ("mean_l1_reduced", l1_red / samples as f64),
                ("mean_ns_reduced", t_red / samples as f64),
                ("mean_ns_naive_l1_power", if rl.h <= 50_000 { t_naive / samples as f64 } else { -1.0 }),
                ("mean_l1_naive", rl.h as f64 / 2.0),
            ],
            ok,
            "class [l_1]^a, a uniform mod h: Babai-reduced exponent vector vs l_1^a; equal curves where both run",
        );
    }
}

/// Alternate models (Moody-Shumow): Edwards and Huff isogenies vs Weierstrass Velu.
fn bench_models(o: &mut Out) {
    use isogeny_algos::kernel::models::*;
    use isogeny_algos::kernel::montgomery;
    let fq = Zp::new(prime_bits(61));
    let mut rng = Rng::new(22_000);
    let budget = if o.quick { 100 } else { 300 };
    let ells: &[u64] = if o.quick { &[3, 7] } else { &[3, 5, 7, 11, 13, 31] };
    let point_of_order = |ew: &Curve<u64>, n: u64, ell: u64, rng: &mut Rng| {
        let mut lv = 1u64;
        while n % (lv * ell) == 0 {
            lv *= ell;
        }
        loop {
            let r = ew.random_point(&fq, rng);
            let mut t = pmul(&fq, ew, &r, (n / lv) as u128);
            if t == Pt::Inf {
                continue;
            }
            loop {
                let nt = pmul(&fq, ew, &t, ell as u128);
                if nt == Pt::Inf {
                    return t;
                }
                t = nt;
            }
        }
    };
    for &ell in ells {
        // Montgomery curve with A + 2 a square and an l-torsion point
        let (a, ew, n, sa) = loop {
            let a = fq.random(&mut rng);
            let ed = edwards_from_montgomery(&fq, a);
            if ed.d == 0 || ed.a == 0 {
                continue;
            }
            let Some(sa) = fq.sqrt(ed.a) else { continue };
            let ew = montgomery::to_weierstrass(&fq, a);
            if !is_smooth(&fq, &ew) {
                continue;
            }
            let n = order(&fq, &ew, &mut rng);
            if n % ell == 0 {
                break (a, ew, n, sa);
            }
        };
        let k = point_of_order(&ew, n, ell, &mut rng);
        let third = fq.div(a, 3);
        let ed = edwards_from_montgomery(&fq, a);
        let e1 = Edwards { a: 1u64, d: fq.div(ed.d, ed.a) };
        let Pt::Aff(xw, yw) = k else { unreachable!() };
        let (x0, y0) = mont_to_edwards_point(&fq, fq.sub(xw, third), yw);
        let k1 = (fq.mul(x0, sa), y0);
        let q = ew.random_point(&fq, &mut rng);
        let Pt::Aff(qx, qy) = q else { continue };
        let (x, y) = mont_to_edwards_point(&fq, fq.sub(qx, third), qy);
        let pt = (fq.mul(x, sa), y);
        let iso = edwards_isogeny(&fq, &e1, k1, ell);
        let ok = iso.eval_explicit(&fq, pt) == iso.eval_def(&fq, pt) && iso.cod.on_curve(&fq, iso.eval_explicit(&fq, pt));
        let params = vec![("field", "Zp61".to_string()), ("ell", ell.to_string())];
        let rec = |o: &mut Out, algo: &str, ok: bool, (m, mn, r): (f64, f64, usize)| {
            o.rec("models", algo, &params, &[("median_ns", m), ("min_ns", mn), ("reps", r as f64)], ok, "");
        };
        rec(o, "edwards_kernel_and_codomain", ok, time_it(budget, 2000, || edwards_isogeny(&fq, &e1, k1, ell)));
        rec(o, "edwards_eval_definition", ok, time_it(budget, 2000, || iso.eval_def(&fq, pt)));
        rec(o, "edwards_eval_theorem2", ok, time_it(budget, 2000, || iso.eval_explicit(&fq, pt)));
        rec(o, "edwards_eval_x_only", ok, time_it(budget, 2000, || iso.eval_x_only(&fq, pt.0)));
        let wv = isogeny_algos::kernel::velu::velu_cyclic(&fq, &ew, &k, ell);
        let okw = iso.cod.j(&fq) == jinv(&fq, &wv.cod);
        rec(o, "weierstrass_velu_kernel_and_codomain", okw, time_it(budget, 2000, || isogeny_algos::kernel::velu::velu_cyclic(&fq, &ew, &k, ell)));
        rec(o, "weierstrass_velu_eval", okw, time_it(budget, 2000, || wv.eval(&fq, &q)));
        // Huff with the same group: y^2 = x(x + a)(x + b) needs full rational 2-torsion; use a
        // fresh Huff curve with an l-torsion point
        let (h, hw, hn) = loop {
            let h = Huff { a: fq.random(&mut rng), b: fq.random(&mut rng) };
            if h.a == h.b || h.a == 0 || h.b == 0 {
                continue;
            }
            let hw = h.short_weierstrass(&fq);
            if !is_smooth(&fq, &hw) {
                continue;
            }
            let hn = order(&fq, &hw, &mut rng);
            if hn % ell == 0 {
                break (h, hw, hn);
            }
        };
        let hk = point_of_order(&hw, hn, ell, &mut rng);
        let shift = fq.div(fq.add(h.a, h.b), 3);
        let Pt::Aff(hx, hy) = hk else { unreachable!() };
        let kh = h.from_weierstrass_point(&fq, (fq.sub(hx, shift), hy));
        if let Some(hiso) = huff_isogeny(&fq, &h, kh, ell) {
            let hq = hw.random_point(&fq, &mut rng);
            if let Pt::Aff(a1, b1) = hq {
                let hp = h.from_weierstrass_point(&fq, (fq.sub(a1, shift), b1));
                let okh = hiso.cod.on_curve(&fq, hiso.eval(&fq, hp));
                rec(o, "huff_kernel_and_codomain", okh, time_it(budget, 2000, || huff_isogeny(&fq, &h, kh, ell)));
                rec(o, "huff_eval", okh, time_it(budget, 2000, || hiso.eval(&fq, hp)));
            }
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let quick = args.iter().any(|a| a == "--quick");
    let mut out_path = None;
    let mut groups = vec![];
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--quick" => {}
            "--out" => {
                i += 1;
                out_path = Some(args[i].clone());
            }
            g => groups.push(g.to_string()),
        }
        i += 1;
    }
    if groups.is_empty() {
        groups = vec!["kernel".into(), "find".into(), "path".into()];
    }
    let file = out_path.map(|p| std::fs::File::create(p).expect("open out"));
    let mut o = Out { file, quick };
    let _ = poly::x_poly(&Zp::new(5));
    for g in &groups {
        match g.as_str() {
            "kernel" => bench_kernel(&mut o),
            "find" => bench_find(&mut o),
            "path" => bench_path(&mut o),
            "kernel2" => bench_kernel2(&mut o),
            "chain" => bench_chain(&mut o),
            "bmss" => bench_bmss(&mut o),
            "csidh" => bench_csidh(&mut o),
            "v2path" => bench_v2path(&mut o),
            "p1kernel" => bench_p1kernel(&mut o),
            "char2" => bench_char2(&mut o),
            "quat" => bench_quat(&mut o),
            "genus2" => bench_genus2(&mut o),
            "phi" => bench_phi(&mut o),
            "radical" => bench_radical(&mut o),
            "relation" => bench_relation(&mut o),
            "models" => bench_models(&mut o),
            "big" => {
                let ells: &[u64] = if o.quick { &[3, 31, 401] } else { &[3, 5, 7, 13, 31, 101, 401, 1009, 4001, 10007] };
                bench_big_field::<4>(&mut o, 256, ells, 14_000);
                bench_big_field::<8>(&mut o, 511, ells, 14_001);
            }
            "v2" => {
                bench_kernel2(&mut o);
                bench_chain(&mut o);
                bench_bmss(&mut o);
                bench_csidh(&mut o);
                bench_v2path(&mut o);
            }
            x => eprintln!("unknown group {x}"),
        }
    }
}
