//! Baseline benchmark harness. Usage:
//!   bench [--quick] [--out FILE] [group ...]      groups: kernel find path (default: all)
//! Single-threaded, fixed seeds. Every timed result is cross-checked for correctness first;
//! records carry `verified`. Output: JSON lines (one record per measurement) + a text table.
use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::{divpoly, elkies, modpoly::Phi};
use isogeny_algos::kernel::{kohel::kohel, sqrt_velu::sqrt_velu, velu::velu_from_point};
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
    fn rec(&mut self, group: &str, algo: &str, params: &[(&str, String)], metrics: &[(&str, f64)], verified: bool, note: &str) {
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
        if ts.len() >= max_reps || (start.elapsed().as_millis() as u64 >= budget_ms && ts.len() >= 3) || start.elapsed().as_millis() as u64 >= budget_ms * 4 {
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
    let ells: &[u64] = if o.quick { &[3, 7, 31, 101] } else { &[3, 5, 7, 11, 13, 31, 61, 101, 211, 401, 1009] };
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
            let params = |extra: &str| vec![("p_bits", bits.to_string()), ("ell", ell.to_string()), ("task", extra.to_string())];
            let mut run = |algo: &str, task: &str, f: &mut dyn FnMut()| {
                let (med, min, reps) = time_it(budget, 2000, || f());
                o.rec("kernel", algo, &params(task), &[("median_ns", med), ("min_ns", min), ("reps", reps as f64)], ok, "");
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
    let ells: &[u64] = if o.quick { &[3, 5, 7] } else { &[3, 5, 7, 11, 13, 17, 19, 23] };
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
            o.rec("find", "phi_setup", &[("p_bits", bits.to_string()), ("ell", ell.to_string())], &[("median_ns", phi_ns), ("reps", 1.0)], true, "one-off per (p,l)");
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
            let params = vec![("p_bits", bits.to_string()), ("ell", ell.to_string()), ("kernels_found", ks_div.len().to_string())];
            let (m, mn, r) = time_it(budget, 200, || divpoly::kernel_polys(&fp, &e, ell, &mut rng.clone()));
            o.rec("find", "divpoly_factor", &params, &[("median_ns", m), ("min_ns", mn), ("reps", r as f64)], ok, "psi_l DDF/EDF + subset search");
            let (m, mn, r) = time_it(budget, 200, || elkies::elkies_isogenies(&fp, &phi, &e, &mut rng.clone()));
            o.rec("find", "elkies_phi_bmss", &params, &[("median_ns", m), ("min_ns", mn), ("reps", r as f64)], ok, "Phi roots + Elkies codomain + BMSS Pade (Phi given)");
            let j = jinv(&fp, &e);
            let (m, mn, r) = time_it(budget, 200, || phi.neighbors(&fp, j, &mut rng.clone()));
            o.rec("find", "phi_roots_only", &params, &[("median_ns", m), ("min_ns", mn), ("reps", r as f64)], ok, "j-invariants of neighbours only, no isogeny");
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
            let params = vec![("p_bits", bits.to_string()), ("instance", i.to_string()), ("trace", t.to_string())];
            // Galbraith
            let t0 = Instant::now();
            let (r, st) = galbraith::galbraith(&fp, &cache, &ells, j1, j2, 5_000_000, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2);
            o.rec("path", "galbraith_bidirectional_bfs", &params, &[("ns", ns), ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)), ("nodes_expanded", st.nodes_expanded as f64)], ok, "primes 3,5,7");
            if let Some(path) = &r {
                let t0 = Instant::now();
                let chain = explicit_chain(&fp, &cache, &e1, path);
                o.rec("path", "explicit_chain_from_path", &params, &[("ns", t0.elapsed().as_nanos() as f64), ("steps", path.len() as f64)], chain.is_some(), "Elkies+BMSS per step");
            }
            // GHS
            let t0 = Instant::now();
            let (r, st) = ghs::ghs(&fp, &cache, &ells, j1, j2, 20_000_000, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2);
            o.rec("path", "ghs_random_walk_collision", &params, &[("ns", ns), ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)), ("steps", st.steps as f64)], ok, "primes 3,5,7");
            let t0 = Instant::now();
            let (r, st) = ghs::ghs_volcano(&fp, &cache, &ells, p, t, j1, j2, 20_000_000, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2);
            o.rec("path", "ghs_with_kohel_volcano", &params, &[("ns", ns), ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)), ("steps", st.steps as f64)], ok, "all heights 0 here, so same walk plus height checks");
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
            let params = vec![("p_bits", bits.to_string()), ("instance", i.to_string()), ("height", h.to_string())];
            let t0 = Instant::now();
            let r = volcano::kohel_volcano_path(&fp, &cache, 3, h, j1, j2, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2);
            o.rec("path", "kohel_volcano_crater_walk", &params, &[("ns", ns), ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64))], ok, "l=3 volcano; same-volcano targets only");
            let t0 = Instant::now();
            let (r, st) = ghs::ghs_volcano(&fp, &cache, &[3, 5, 7, 11, 13], p, t, j1, j2, 200_000, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| verify_path(&fp, &cache, p) && *p.js.last().unwrap() == j2);
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
            let mut act = couveignes::Action { fp: &fp, cache: &cache, plus_class: Default::default() };
            let t0 = Instant::now();
            let primes = couveignes::select_primes(&mut act, j1, &mut rng);
            let sel_ns = t0.elapsed().as_nanos() as f64;
            let k = primes.len().min(4);
            if k < 2 {
                continue;
            }
            let primes = &primes[..k];
            let m = 4usize;
            let secret: Vec<i64> = primes.iter().map(|_| rng.below(2 * m as u64 + 1) as i64 - m as i64).collect();
            let j2 = *act.apply(j1, primes, &secret, &mut rng).unwrap().js.last().unwrap();
            let params = vec![("p_bits", bits.to_string()), ("instance", i.to_string()), ("primes", format!("{:?}", primes)), ("box_m", m.to_string())];
            let t0 = Instant::now();
            let (r, st) = couveignes::couveignes(&act, primes, m, j1, j2, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = match &r {
                Some((e, path)) => {
                    let chk = act.apply(j1, primes, e, &mut rng.clone()).map_or(false, |q| *q.js.last().unwrap() == j2);
                    chk && verify_path(&fp, &cache, path)
                }
                None => false,
            };
            o.rec("path", "couveignes_hhs_mitm", &params, &[("ns", ns), ("nodes", st.nodes as f64), ("path_len", r.as_ref().map_or(-1.0, |x| x.1.len() as f64)), ("prime_select_ns", sel_ns), ("phi_setup_ns", setup)], ok, "planted exponents in [-m,m]^k; different problem from generic path finding");
        }
    }
    // Delfs-Galbraith
    let dbits: &[u32] = if o.quick { &[22] } else { &[22, 26, 30, 34, 38] };
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
            let j1 = *random_walk(&f2, &cache, &[2], start, steps, &mut rng).js.last().unwrap();
            let j2 = *random_walk(&f2, &cache, &[2], start, steps, &mut rng).js.last().unwrap();
            let params = vec![("p_bits", bits.to_string()), ("instance", i.to_string())];
            let t0 = Instant::now();
            let (r, st) = delfs_galbraith::delfs_galbraith(&f2, &fp, &cache, j1, j2, 50_000_000, 5_000_000, &mut rng.clone());
            let ns = t0.elapsed().as_nanos() as f64;
            let ok = r.as_ref().map_or(false, |p| verify_path(&f2, &cache, p) && p.js[0] == j1 && *p.js.last().unwrap() == j2);
            o.rec("path", "delfs_galbraith_supersingular", &params, &[("ns", ns), ("path_len", r.as_ref().map_or(-1.0, |p| p.len() as f64)), ("walk_steps", st.walk_steps as f64), ("bfs_nodes", st.bfs_nodes as f64)], ok, "j-path only; Fp2 arithmetic");
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
            x => eprintln!("unknown group {x}"),
        }
    }
}
