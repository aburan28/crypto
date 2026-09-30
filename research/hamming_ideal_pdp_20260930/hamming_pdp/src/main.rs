//! Hamming-ideal point decomposition on `K_0/F_{2^n}`: see PROTOCOL.md.

mod boolpoly;
mod curve;
mod f4;
mod gf2n;
mod hamming;
mod multisolve;
mod selftest;
mod semaev;

use boolpoly::{Poly, W};
use curve::{Curve, Point};
use gf2n::Field;
use hamming::{hamming_ideal, Alloc, Encoding};
use multisolve::{multisolve, Config, Instance};
use semaev::{s3_descent, s3_value, SymElem};
use std::collections::HashSet;
use std::io::Write;
use std::time::Instant;

struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Rng {
        Rng(seed.wrapping_mul(0x9E3779B97F4A7C15) ^ 0xD1B54A32D192ED03 | 1)
    }
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

struct Base {
    name: String,
    /// abscissae
    xs: HashSet<u64>,
    points: Vec<Point>,
    /// weight cap (normal basis) or None for a subspace
    w: Option<u32>,
    /// subspace basis (polynomial-basis bits), if any
    sub: Vec<u64>,
    x_digest: u64,
}

fn fnv(acc: u64, v: u64) -> u64 {
    (acc ^ v).wrapping_mul(0x100000001b3)
}

fn weight_base(f: &Field, c: &Curve, w: u32) -> Base {
    let mut xs = HashSet::new();
    let mut points = Vec::new();
    let mut digest = 0xcbf29ce484222325u64;
    for x in 0..(1u64 << f.n) {
        if f.to_normal(x).count_ones() <= w {
            if let Some(ps) = c.lift(x) {
                xs.insert(x);
                digest = fnv(digest, x);
                points.push(ps[0]);
                if ps[1] != ps[0] {
                    points.push(ps[1]);
                }
            }
        }
    }
    Base { name: format!("WT{w}"), xs, points, w: Some(w), sub: vec![], x_digest: digest }
}

fn subspace_base(f: &Field, c: &Curve, l: usize, rng: &mut Rng) -> Base {
    let mut basis: Vec<u64> = Vec::new();
    while basis.len() < l {
        let cand = rng.next() & f.mask;
        let mut trial = basis.clone();
        trial.push(cand);
        if gf2n::rank_f2(&trial) == trial.len() {
            basis = trial;
        }
    }
    let mut xs = HashSet::new();
    let mut points = Vec::new();
    let mut digest = 0xcbf29ce484222325u64;
    for u in 0..(1u64 << l) {
        let mut x = 0;
        for (i, b) in basis.iter().enumerate() {
            if (u >> i) & 1 == 1 {
                x ^= b;
            }
        }
        if let Some(ps) = c.lift(x) {
            xs.insert(x);
            digest = fnv(digest, x);
            points.push(ps[0]);
            if ps[1] != ps[0] {
                points.push(ps[1]);
            }
        }
    }
    Base { name: format!("SUB{l}"), xs, points, w: None, sub: basis, x_digest: digest }
}

struct Exhaustive {
    ordered_pairs: u64,
    first_hit: Option<u64>,
    scanned: u64,
}

fn exhaustive(c: &Curve, base: &Base, r: Point) -> Exhaustive {
    let mut count = 0;
    let mut first = None;
    let mut scanned = 0;
    for (i, &p) in base.points.iter().enumerate() {
        scanned += 1;
        let q = c.sub(r, p);
        if let Some(xq) = c.x(q) {
            if base.xs.contains(&xq) {
                count += 1;
                if first.is_none() {
                    first = Some(i as u64 + 1);
                }
            }
        }
    }
    Exhaustive { ordered_pairs: count, first_hit: first, scanned }
}

fn random_point(c: &Curve, rng: &mut Rng) -> Point {
    loop {
        let x = rng.next() & c.f.mask;
        if let Some(ps) = c.lift(x) {
            let p = ps[rng.below(2) as usize];
            let q = c.mul(4, p);
            if q != Point::Inf {
                return q;
            }
        }
    }
}

fn json_list<T: std::fmt::Display>(v: &[T]) -> String {
    let s: Vec<String> = v.iter().map(|x| x.to_string()).collect();
    format!("[{}]", s.join(","))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut n: u32 = 7;
    let mut w: u32 = 2;
    let mut l: usize = 0;
    let mut targets: u64 = 24;
    let mut encodings: Vec<Encoding> = vec![Encoding::Sub, Encoding::Mono, Encoding::C, Encoding::Fc, Encoding::Qfc];
    let mut d_max = 6;
    let mut budget: u64 = 1 << 30;
    let mut max_calls: u64 = 1 << 13;
    let mut matrix_cap_log2: u32 = 27;
    let mut out = String::from("/dev/stdout");
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--selftest" => { selftest::run(args[i + 1].parse().unwrap(), 300); return; }
            "--n" => { n = args[i + 1].parse().unwrap(); i += 2; }
            "--w" => { w = args[i + 1].parse().unwrap(); i += 2; }
            "--l" => { l = args[i + 1].parse().unwrap(); i += 2; }
            "--targets" => { targets = args[i + 1].parse().unwrap(); i += 2; }
            "--d" => { d_max = args[i + 1].parse().unwrap(); i += 2; }
            "--budget-log2" => { budget = 1u64 << args[i + 1].parse::<u32>().unwrap(); i += 2; }
            "--matrix-cap-log2" => { matrix_cap_log2 = args[i + 1].parse().unwrap(); i += 2; }
            "--max-calls" => { max_calls = args[i + 1].parse().unwrap(); i += 2; }
            "--encodings" => { encodings = args[i + 1].split(',').map(|s| Encoding::parse(s).expect("encoding")).collect(); i += 2; }
            "--out" => { out = args[i + 1].clone(); i += 2; }
            _ => panic!("unknown arg {}", args[i]),
        }
    }
    let t_setup = Instant::now();
    let f = Field::new(n);
    let c = Curve::new(&f);
    let n_x_wt: u64 = (0..=w).map(|j| binom(n as u64, j as u64)).sum();
    if l == 0 {
        l = (64 - (n_x_wt - 1).leading_zeros()) as usize; // ceil(log2)
    }
    let mut rng_base = Rng::new(0x5eed_0000 + n as u64);
    let wt = weight_base(&f, &c, w);
    let sub = subspace_base(&f, &c, l, &mut rng_base);
    let setup_ms = t_setup.elapsed().as_millis();
    eprintln!(
        "n={n} irr={:#x} alpha={:#x} #E={} odd={} |F_wt|={} (x:{}) |F_sub|={} (x:{}) l={l} setup={}ms",
        f.irr, f.nb[0], c.order, c.odd_order, wt.points.len(), wt.xs.len(), sub.points.len(), sub.xs.len(), setup_ms
    );

    // Targets: even seeds planted from the weight base, odd seeds uniform.
    let mut tg: Vec<(u64, Point, &str)> = Vec::new();
    for s in 0..targets {
        let mut rng = Rng::new((n as u64) << 32 | s);
        if s % 2 == 0 {
            loop {
                let p1 = wt.points[rng.below(wt.points.len() as u64) as usize];
                let p2 = wt.points[rng.below(wt.points.len() as u64) as usize];
                let r = c.add(p1, p2);
                if r != Point::Inf {
                    tg.push((s, r, "planted"));
                    break;
                }
            }
        } else {
            tg.push((s, random_point(&c, &mut rng), "uniform"));
        }
    }

    let mut file = std::fs::File::create(&out).expect("open output");
    let cfg = Config { d_max, budget_xor: budget, max_calls, max_free: 12, matrix_cap_words: 1u64 << matrix_cap_log2 };
    let nn = n as usize;

    // Manifest line.
    writeln!(
        file,
        "{{\"kind\":\"manifest\",\"n\":{n},\"irr\":\"{:#x}\",\"normal_alpha\":\"{:#x}\",\"order\":{},\"odd_order\":{},\"w\":{w},\"l\":{l},\"wt_points\":{},\"wt_x\":{},\"wt_x_digest\":\"{:016x}\",\"sub_points\":{},\"sub_x\":{},\"sub_x_digest\":\"{:016x}\",\"sub_basis\":{},\"d_max\":{d_max},\"budget_xor\":{budget},\"max_calls\":{max_calls},\"matrix_cap_words\":{},\"targets\":{targets},\"setup_ms\":{setup_ms}}}",
        f.irr, f.nb[0], c.order, c.odd_order, wt.points.len(), wt.xs.len(), wt.x_digest, sub.points.len(), sub.xs.len(), sub.x_digest,
        json_list(&sub.sub.iter().map(|b| format!("\"{b:#x}\"")).collect::<Vec<_>>()),
        1u64 << matrix_cap_log2
    ).unwrap();

    for enc in &encodings {
        let base = if *enc == Encoding::Sub { &sub } else { &wt };
        let coords_per = if *enc == Encoding::Sub { l } else { nn };
        for (seed, r, kind) in &tg {
            let xr = c.x(*r).unwrap();
            let t0 = Instant::now();
            let ex = exhaustive(&c, base, *r);
            let ex_ms = t0.elapsed().as_micros();

            // Build the system.
            let t1 = Instant::now();
            let mut alloc = Alloc { next: 2 * coords_per };
            let sv: Vec<Vec<usize>> = (0..2).map(|i| (0..coords_per).map(|j| i * coords_per + j).collect()).collect();
            let (x1, x2) = if *enc == Encoding::Sub {
                (SymElem::from_vars(&f, &sub.sub, &sv[0]), SymElem::from_vars(&f, &sub.sub, &sv[1]))
            } else {
                (SymElem::from_vars(&f, &f.nb, &sv[0]), SymElem::from_vars(&f, &f.nb, &sv[1]))
            };
            let mut eqs = s3_descent(&f, &x1, &x2, xr);
            let s3_eqs = eqs.len();
            let mut aux = 0;
            for i in 0..2 {
                let b = hamming_ideal(*enc, &sv[i], w, &mut alloc);
                aux += b.aux_vars;
                eqs.extend(b.eqs);
            }
            let n_vars = alloc.next;
            assert!(n_vars <= boolpoly::MAX_VARS);
            let max_deg = eqs.iter().map(|e| e.degree()).max().unwrap_or(0);
            let build_ms = t1.elapsed().as_micros();
            let mut branch = Vec::new();
            for i in 0..2 {
                for &v in &sv[i] {
                    branch.push((v, i));
                }
            }
            let inst = Instance {
                eqs,
                n_vars,
                branch,
                summand_coords: sv.clone(),
                weight_cap: base.w,
            };

            // Verification closure: lift a Boolean point to abscissae and points.
            let decode = |pt: &[u64; W], i: usize| -> u64 {
                let mut bits = 0u64;
                for (j, &v) in sv[i].iter().enumerate() {
                    if (pt[v / 64] >> (v % 64)) & 1 == 1 {
                        bits |= 1 << j;
                    }
                }
                if *enc == Encoding::Sub {
                    let mut x = 0;
                    for (j, b) in sub.sub.iter().enumerate() {
                        if (bits >> j) & 1 == 1 {
                            x ^= b;
                        }
                    }
                    x
                } else {
                    f.from_normal(bits)
                }
            };
            let mut witness: Option<(u64, u64)> = None;
            let mut encoding_violations = 0u64;
            let mut verify = |pt: &[u64; W]| -> bool {
                let a = decode(pt, 0);
                let b = decode(pt, 1);
                if let Some(wc) = base.w {
                    if f.to_normal(a).count_ones() > wc || f.to_normal(b).count_ones() > wc {
                        encoding_violations += 1;
                        return false;
                    }
                }
                if s3_value(&f, a, b, xr) != 0 {
                    return false;
                }
                let (pa, pb) = match (c.lift(a), c.lift(b)) {
                    (Some(pa), Some(pb)) => (pa, pb),
                    _ => return false,
                };
                for p in pa {
                    for q in pb {
                        if c.add(p, q) == *r {
                            witness = Some((a, b));
                            return true;
                        }
                    }
                }
                false
            };
            let t2 = Instant::now();
            let st = multisolve(&inst, &cfg, &mut verify);
            let solve_ms = t2.elapsed().as_micros();
            let found = !st.solutions.is_empty();
            let label = ex.ordered_pairs > 0;
            let agree = if st.exhausted { "timeout" } else if found == label { "ok" } else { "MISMATCH" };
            let mean_tame_depth = if st.tame_depths.is_empty() { -1.0 } else { st.tame_depths.iter().sum::<usize>() as f64 / st.tame_depths.len() as f64 };
            eprintln!(
                "n={n} {} {} seed={seed} {kind} label={} pairs={} found={found} {agree} calls={} tame={} wild={} budget={} cap={} depth_mean={:.2} max_depth={} xor={} rows={} cols={} deg={} prep={}ms elim={}ms add={}ms rref={}ms subst={}ms post={}ms branch={}ms {}ms",
                enc.name(), base.name, label, ex.ordered_pairs, st.calls, st.tame, st.wild, st.budget, st.matrix_cap_hits, mean_tame_depth, st.max_depth, st.xor_total, st.rows_max, st.cols_max, st.max_degree, st.prep_ns / 1_000_000, st.elim_ns / 1_000_000, st.add_ns / 1_000_000, st.rref_ns / 1_000_000, st.subst_ns / 1_000_000, st.post_ns / 1_000_000, st.branch_ns / 1_000_000, solve_ms / 1000
            );
            writeln!(
                file,
                "{{\"kind\":\"cell\",\"n\":{n},\"encoding\":\"{}\",\"base\":\"{}\",\"seed\":{seed},\"target_kind\":\"{kind}\",\"x_r\":\"{xr:#x}\",\"label\":{label},\"exhaustive_ordered_pairs\":{},\"exhaustive_first_hit\":{},\"exhaustive_scanned\":{},\"exhaustive_us\":{ex_ms},\"vars\":{n_vars},\"aux_vars\":{aux},\"gens\":{},\"s3_gens\":{s3_eqs},\"max_gen_degree\":{max_deg},\"build_us\":{build_ms},\"found\":{found},\"agree\":\"{agree}\",\"witness\":{},\"calls\":{},\"tame\":{},\"wild\":{},\"budget_calls\":{},\"inconsistent\":{},\"max_depth\":{},\"tame_depths\":{},\"xor_words\":{},\"rows_max\":{},\"cols_max\":{},\"max_f4_degree\":{},\"basis_max\":{},\"pruned\":{},\"forced\":{},\"unresolved_leaves\":{},\"candidates_rejected\":{},\"encoding_violations\":{encoding_violations},\"matrix_cap_hits\":{},\"leaf_direct\":{},\"prep_ns\":{},\"elim_ns\":{},\"exhausted\":{},\"solve_us\":{solve_ms}}}",
                enc.name(), base.name, ex.ordered_pairs,
                ex.first_hit.map(|v| v.to_string()).unwrap_or("null".into()), ex.scanned,
                inst.eqs.len(),
                witness.map(|(a, b)| format!("[\"{a:#x}\",\"{b:#x}\"]")).unwrap_or("null".into()),
                st.calls, st.tame, st.wild, st.budget, st.inconsistent, st.max_depth, json_list(&st.tame_depths),
                st.xor_total, st.rows_max, st.cols_max, st.max_degree, st.basis_max, st.pruned, st.forced, st.unresolved_leaves, st.candidates_rejected, st.matrix_cap_hits, st.leaf_direct, st.prep_ns, st.elim_ns, st.exhausted
            ).unwrap();
            let _ = Poly::zero();
        }
    }
}

fn binom(n: u64, k: u64) -> u64 {
    let mut r = 1u64;
    for i in 0..k {
        r = r * (n - i) / (i + 1);
    }
    r
}
