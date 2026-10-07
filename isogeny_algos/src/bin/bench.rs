//! Baseline benchmark harness. Usage:
//!   bench [--quick] [--out FILE] [group ...]      groups: kernel find path (default: all)
//! Single-threaded, fixed seeds. Every timed result is cross-checked for correctness first;
//! records carry `verified`. Output: JSON lines (one record per measurement) + a text table.
use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::{bmss, dual};
use isogeny_algos::find::{divpoly, elkies, modpoly::Phi};
use isogeny_algos::kernel::chain::{ell_power_isogeny, mul_by_ell, ChainStats, Strategy};
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
        let cache = PhiCache::new(p, &ells);
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
        let cache = PhiCache::new(p, &[3, 5, 7, 11, 13]);
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
        let cache = PhiCache::new(p, &cells);
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
        let cache = PhiCache::new(p, &[2, 3]);
        for i in 0..inst {
            let mut rng = Rng::new(6000 + (bits as u64) * 100 + i);
            let start = (1728u64, 0u64);
            let steps = 3 * bits as usize;
            let j1 = *random_walk(&f2, &cache, &[2], start, steps, &mut rng)
                .js
                .last()
                .unwrap();
            let j2 = *random_walk(&f2, &cache, &[2], start, steps, &mut rng)
                .js
                .last()
                .unwrap();
            let params = vec![("p_bits", bits.to_string()), ("instance", i.to_string())];
            let t0 = Instant::now();
            let (r, st) = delfs_galbraith::delfs_galbraith(
                &f2,
                &fp,
                &cache,
                j1,
                j2,
                50_000_000,
                5_000_000,
                &mut rng.clone(),
            );
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| {
                verify_path(&f2, &cache, p) && p.js[0] == j1 && *p.js.last().unwrap() == j2
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
        let cache = PhiCache::new(p, &[3, 5, 7]);
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
        let cache = PhiCache::new(p, &ells);
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
