//! dorv: replay, controls and verdict analysis for the binary DoR campaign
//! (research/factor_base_dor_campaign_20261013). See ../PROTOCOL.md.
mod engine;
mod json;

use engine::*;
use json::{dump, n as jn, obj, parse, st, V};
use std::collections::BTreeMap;
use std::time::Instant;

struct SplitMix(u64);
impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 { self.next() % n }
}

fn hexs(v: &str) -> u64 { u64::from_str_radix(v.trim_start_matches("0x"), 16).unwrap() }
fn read(p: &str) -> V { parse(&std::fs::read_to_string(p).unwrap_or_else(|e| panic!("{}: {}", p, e))) }

fn profile_v(r: &Refute) -> (V, V) {
    let p = V::Obj(r.profile.iter().map(|(d, k)| (d.to_string(), jn(*k as f64))).collect());
    let s = V::Obj(r.shapes.iter().map(|(d, a, b)| (d.to_string(), V::Arr(vec![jn(*a as f64), jn(*b as f64)]))).collect());
    (p, s)
}

fn same_profile(rec: Option<&V>, mine: &V) -> bool {
    match rec {
        None => false,
        Some(r) => dump(r).replace([' ', '\n'], "") == dump(mine).replace([' ', '\n'], ""),
    }
}

// ── replay of one recorded m=2 grid cell ──────────────────────────────────
fn replay_grid(path: &str, seed: u64, out: &str, deadline: f64, t0: &Instant) {
    let rec = read(path);
    let n = rec.get("field_degree").unwrap().num().unwrap() as usize;
    let dim = rec.get("dim").unwrap().num().unwrap() as usize;
    let dmax = rec.get("dmax").unwrap().num().unwrap() as u32;
    let nt = rec.get("targets").unwrap().num().unwrap() as usize;
    let f = Field::new(n);
    let tg = targets(n, nt, seed);
    let mut structs = vec![];
    for r in rec.get("results").unwrap().arr().unwrap() {
        let name = r.get("structure").unwrap().str().unwrap().to_string();
        if r.get("skipped").is_some() { continue; }
        let (b1, b2) = structure_bases(&f, &name, dim, seed).expect("structure");
        let rb1: Vec<u64> = r.get("basis").unwrap().arr().unwrap().iter().map(|x| hexs(x.str().unwrap())).collect();
        let rb2: Vec<u64> = match r.get("basis2") { Some(V::Arr(a)) => a.iter().map(|x| hexs(x.str().unwrap())).collect(), _ => rb1.clone() };
        let basis_match = rb1 == b1 && rb2 == b2;
        let m1 = span(&b1);
        let m2 = span(&b2);
        let unsat_list: Vec<u64> = tg.iter().cloned().filter(|&t| unsat(&f, &m1, &m2, t, 1)).collect();
        let frac = 1.0 - unsat_list.len() as f64 / tg.len() as f64;
        let pts = fb_points(&f, &b1);
        let (r_vv, vv) = product_dim(&f, &b1, &b2);
        let mut sum = b1.clone(); sum.extend_from_slice(&b2);
        let l_pred = rank_u64(&sum);
        let rec_targets = r.get("targets").unwrap().arr().unwrap();
        let rec_list: Vec<u64> = rec_targets.iter().map(|x| hexs(x.get("target").unwrap().str().unwrap())).collect();
        let mut rows = vec![];
        let mut nmatch = 0; let mut ndone = 0; let mut e1_ok = 0; let mut nshape = 0;
        for (k, &t) in unsat_list.iter().enumerate() {
            if t0.elapsed().as_secs_f64() > deadline { break; }
            let eqs = descent_equations(&f, &b1, &b2, t, 1);
            let nvars = b1.len() + b2.len();
            let rf = refutation_degree(&eqs, nvars, dmax);
            let (rank_eq, l, q) = part_ranks(&eqs, nvars);
            let x3in = in_span(&vv, t);
            let q_pred = r_vv - if x3in { 1 } else { 0 };
            let e1 = l == l_pred && q == q_pred && rank_eq <= 1 + l + q;
            if e1 { e1_ok += 1; }
            let (pv, sv) = profile_v(&rf);
            let recd = rec_targets.get(k);
            let m = recd.map(|x| {
                hexs(x.get("target").unwrap().str().unwrap()) == t
                    && x.get("d_star").and_then(|v| v.num()).map(|v| v as u32) == rf.dstar
                    && same_profile(x.get("rank_profile"), &pv)
                    && (x.get("macaulay_shape").is_none() || same_profile(x.get("macaulay_shape"), &sv))
            }).unwrap_or(false);
            let shape_rec = recd.map(|x| x.get("macaulay_shape").is_some()).unwrap_or(false);
            if shape_rec { nshape += 1; }
            if m { nmatch += 1; }
            ndone += 1;
            rows.push(obj(vec![
                ("target", st(&format!("0x{:x}", t))),
                ("d_star", rf.dstar.map(|d| jn(d as f64)).unwrap_or(V::Null)),
                ("rank_profile", pv), ("macaulay_shape", sv),
                ("l", jn(l as f64)), ("q", jn(q as f64)), ("rank_eq", jn(rank_eq as f64)),
                ("x3_in_VV", V::Bool(x3in)), ("e1_ok", V::Bool(e1)), ("replay_match", V::Bool(m)),
            ]));
        }
        let rec_frac = r.get("decomposable_fraction").unwrap().num().unwrap();
        let rec_pts = r.get("factor_base_points").unwrap().num().unwrap() as usize;
        structs.push(obj(vec![
            ("structure", st(&name)), ("dim", jn(dim as f64)),
            ("basis_match", V::Bool(basis_match)),
            ("unsat_list_match", V::Bool(rec_list == unsat_list)),
            ("decomposable_fraction", jn(frac)), ("fraction_match", V::Bool((rec_frac - frac).abs() < 1e-12)),
            ("factor_base_points", jn(pts as f64)), ("points_match", V::Bool(rec_pts == pts)),
            ("dim_VV", jn(r_vv as f64)), ("dim_V1_plus_V2", jn(l_pred as f64)),
            ("unsat_targets", jn(unsat_list.len() as f64)), ("replayed", jn(ndone as f64)),
            ("replay_matches", jn(nmatch as f64)), ("shape_recorded", jn(nshape as f64)), ("e1_ok", jn(e1_ok as f64)),
            ("targets", V::Arr(rows)),
        ]));
    }
    let o = obj(vec![
        ("source", st(path)), ("field_degree", jn(n as f64)), ("modulus", st(&format!("0x{:x}", f.f))),
        ("dim", jn(dim as f64)), ("seed", jn(seed as f64)), ("dmax", jn(dmax as f64)),
        ("results", V::Arr(structs)),
    ]);
    std::fs::write(out, dump(&o) + "\n").unwrap();
}

// ── replay of a delta sweep (standard, m=2, seed 1) with per-target D* ────
fn replay_delta(path: &str, out: &str, deadline: f64, t0: &Instant) {
    let rec = read(path);
    let n = rec.get("field_degree").unwrap().num().unwrap() as usize;
    let dmax = rec.get("dmax").unwrap().num().unwrap() as u32;
    let nt = rec.get("targets").unwrap().num().unwrap() as usize;
    let f = Field::new(n);
    let tg = targets(n, nt, 1);
    let mut rows = vec![];
    for r in rec.get("rows").unwrap().arr().unwrap() {
        if r.get("skipped").is_some() { continue; }
        let delta = r.get("delta").unwrap().num().unwrap();
        let dim = r.get("dim").unwrap().num().unwrap() as usize;
        let (b1, b2) = structure_bases(&f, "standard", dim, 1).unwrap();
        let m1 = span(&b1);
        let un: Vec<u64> = tg.iter().cloned().filter(|&t| unsat(&f, &m1, &m1, t, 1)).collect();
        let frac = 1.0 - un.len() as f64 / tg.len() as f64;
        let mut ds = vec![]; let mut rows_sh = vec![]; let mut complete = true;
        for &t in &un {
            if t0.elapsed().as_secs_f64() > deadline { complete = false; break; }
            let eqs = descent_equations(&f, &b1, &b2, t, 1);
            let rf = refutation_degree(&eqs, 2 * dim, dmax);
            ds.push(rf.dstar.map(|d| jn(d as f64)).unwrap_or(V::Null));
            if let Some(d) = rf.dstar { let s = rf.shapes.iter().find(|x| x.0 == d).unwrap(); rows_sh.push((s.1, s.2)); }
        }
        let solved: Vec<f64> = ds.iter().filter_map(|v| v.num()).collect();
        let mean = if solved.is_empty() { f64::NAN } else { solved.iter().sum::<f64>() / solved.len() as f64 };
        let rm = rows_sh.iter().map(|x| x.0 as f64).sum::<f64>() / rows_sh.len().max(1) as f64;
        let cm = rows_sh.iter().map(|x| x.1 as f64).sum::<f64>() / rows_sh.len().max(1) as f64;
        let rec_mean = r.get("d_star_mean").and_then(|v| v.num()).unwrap_or(f64::NAN);
        let rec_frac = r.get("decomposable_fraction").unwrap().num().unwrap();
        rows.push(obj(vec![
            ("delta", jn(delta)), ("dim", jn(dim as f64)), ("complete", V::Bool(complete)),
            ("decomposable_fraction", jn(frac)), ("fraction_match", V::Bool((frac - rec_frac).abs() < 1e-12)),
            ("d_star_mean", jn(mean)), ("mean_match", V::Bool(complete && (mean - rec_mean).abs() < 1e-9)),
            ("macaulay_rows_mean", jn(rm)), ("macaulay_cols_mean", jn(cm)),
            ("d_star", V::Arr(ds)),
        ]));
    }
    let o = obj(vec![("source", st(path)), ("field_degree", jn(n as f64)), ("rows", V::Arr(rows))]);
    std::fs::write(out, dump(&o) + "\n").unwrap();
}

// ── N2 invariant null: shape-matched random bilinear systems ──────────────
static mut SYM_LIN: bool = false;
fn null_system(rng: &mut SplitMix, neq: usize, d1: usize, d2: usize, l: usize, q: usize) -> Vec<Poly> {
    let nv = d1 + d2;
    let sym = unsafe { SYM_LIN };
    // l independent linear forms over nv vars (exploratory variant: in the a_i + b_i forms)
    let lin: Vec<u64> = loop {
        let v: Vec<u64> = (0..l).map(|_| if sym {
            let m = rng.below(1u64 << d1); m | (m << d1)
        } else { rng.below(1u64 << nv) }).collect();
        if rank_u64(&v) == l { break v; }
    };
    // q independent bilinear forms as d1*d2-bit matrices
    let qb = d1 * d2;
    let quad: Vec<u64> = loop {
        let v: Vec<u64> = (0..q).map(|_| if qb == 64 { rng.next() } else { rng.below(1u64 << qb) }).collect();
        if rank_u64(&v) == q { break v; }
    };
    let mut eqs = vec![];
    while eqs.len() < neq {
        let mut lm = 0u64; let mut qm = 0u64;
        for &x in &lin { if rng.next() & 1 == 1 { lm ^= x; } }
        for &x in &quad { if rng.next() & 1 == 1 { qm ^= x; } }
        let c = rng.next() & 1;
        let mut p: Poly = vec![];
        if c == 1 { p.push(0); }
        for i in 0..nv { if (lm >> i) & 1 == 1 { p.push(1u32 << i); } }
        for i in 0..d1 { for j in 0..d2 { if (qm >> (i * d2 + j)) & 1 == 1 { p.push((1u32 << i) | (1u32 << (d1 + j))); } } }
        p.sort_unstable();
        if !p.is_empty() { eqs.push(p); }
    }
    eqs
}

fn bool_sat(eqs: &[Poly], nv: usize) -> bool {
    'a: for a in 0u32..(1u32 << nv) {
        for e in eqs {
            let mut v = 0;
            for &m in e { if a & m == m { v ^= 1; } }
            if v == 1 { continue 'a; }
        }
        return true;
    }
    false
}

fn null_cell(neq: usize, d1: usize, d2: usize, l: usize, q: usize, dmax: u32, count: usize, seed: u64) -> (Vec<Option<u32>>, usize) {
    let mut rng = SplitMix(seed);
    let mut out = vec![];
    let mut rejected = 0;
    while out.len() < count {
        let eqs = null_system(&mut rng, neq, d1, d2, l, q);
        if bool_sat(&eqs, d1 + d2) { rejected += 1; if rejected > 200 * count { break; } continue; }
        out.push(refutation_degree(&eqs, d1 + d2, dmax).dstar);
    }
    (out, rejected)
}

fn run_null(replay_dir: &str, out: &str, count: usize, deadline: f64, t0: &Instant) {
    let mut cells = vec![];
    let mut files: Vec<_> = std::fs::read_dir(replay_dir).unwrap().filter_map(|e| e.ok()).map(|e| e.path())
        .filter(|p| p.file_name().unwrap().to_str().unwrap().starts_with("m2_")).collect();
    files.sort();
    // pool seeds: key (n, dim, structure) -> modal (l, q) and list of seeds
    let mut groups: BTreeMap<(usize, usize, String), Vec<(usize, usize)>> = BTreeMap::new();
    let mut dmaxs: BTreeMap<usize, u32> = BTreeMap::new();
    for p in &files {
        let v = read(p.to_str().unwrap());
        let n = v.get("field_degree").unwrap().num().unwrap() as usize;
        if n > 15 { continue; } // PROTOCOL: N2 on R1 cells (n <= 15) only
        let dim = v.get("dim").unwrap().num().unwrap() as usize;
        dmaxs.insert(n, v.get("dmax").unwrap().num().unwrap() as u32);
        for r in v.get("results").unwrap().arr().unwrap() {
            let name = r.get("structure").unwrap().str().unwrap().to_string();
            if unsafe { SYM_LIN } && name == "nb_break" { continue; }
            for t in r.get("targets").unwrap().arr().unwrap() {
                let l = t.get("l").unwrap().num().unwrap() as usize;
                let q = t.get("q").unwrap().num().unwrap() as usize;
                groups.entry((n, dim, name.clone())).or_default().push((l, q));
            }
        }
    }
    for ((n, dim, name), lq) in groups {
        if t0.elapsed().as_secs_f64() > deadline { break; }
        let mut cnt: BTreeMap<(usize, usize), usize> = BTreeMap::new();
        for x in &lq { *cnt.entry(*x).or_default() += 1; }
        let (&(l, q), _) = cnt.iter().max_by_key(|x| *x.1).unwrap();
        let (d1, d2) = (dim, dim);
        let seed = 20261015u64 ^ ((n as u64) << 32) ^ ((dim as u64) << 16) ^ (q as u64) << 8 ^ (l as u64);
        let (ds, rej) = null_cell(n, d1, d2, l, q, dmaxs[&n], count, seed);
        cells.push(obj(vec![
            ("field_degree", jn(n as f64)), ("dim", jn(dim as f64)), ("structure", st(&name)),
            ("l", jn(l as f64)), ("q", jn(q as f64)), ("modal_share", jn(*cnt.get(&(l, q)).unwrap() as f64 / lq.len() as f64)),
            ("count", jn(ds.len() as f64)), ("rejected_sat", jn(rej as f64)),
            ("d_star", V::Arr(ds.iter().map(|d| d.map(|x| jn(x as f64)).unwrap_or(V::Null)).collect())),
        ]));
    }
    // P2b: planted q small vs large at n=13, dim=4 (l=4): q=4 vs q=13
    let mut p2 = vec![];
    for &q in &[4usize, 7, 10, 13] {
        let (ds, _) = null_cell(13, 4, 4, 4, q, 6, count, 777 + q as u64);
        p2.push(obj(vec![("n", jn(13.0)), ("dim", jn(4.0)), ("l", jn(4.0)), ("q", jn(q as f64)),
            ("d_star", V::Arr(ds.iter().map(|d| d.map(|x| jn(x as f64)).unwrap_or(V::Null)).collect()))]));
    }
    let o = obj(vec![("schema", st("dorv-null-v1")), ("cells", V::Arr(cells)), ("planted_q_control", V::Arr(p2))]);
    std::fs::write(out, dump(&o) + "\n").unwrap();
}

// ── E2: product-space dimension over all / sampled subspaces ──────────────
fn rref_subspaces(n: usize, d: usize, f: &mut dyn FnMut(&[u64])) {
    // pivots p_0 > p_1 > ... (highest bit leads); free bits: non-pivot positions below each pivot
    fn rec(n: usize, d: usize, start: usize, piv: &mut Vec<usize>, f: &mut dyn FnMut(&[u64])) {
        if piv.len() == d {
            let pset: u64 = piv.iter().fold(0, |a, &p| a | (1 << p));
            let mut frees: Vec<Vec<usize>> = vec![];
            for &p in piv.iter() { frees.push((0..p).filter(|b| (pset >> b) & 1 == 0).collect()); }
            let total: usize = frees.iter().map(|x| x.len()).sum();
            let mut rows = vec![0u64; d];
            for mask in 0u64..(1u64 << total) {
                let mut k = 0;
                for (i, &p) in piv.iter().enumerate() {
                    let mut r = 1u64 << p;
                    for &b in &frees[i] { if (mask >> k) & 1 == 1 { r |= 1 << b; } k += 1; }
                    rows[i] = r;
                }
                f(&rows);
            }
            return;
        }
        for p in start..n { piv.push(p); rec(n, d, p + 1, piv, f); piv.pop(); }
    }
    let mut piv = vec![];
    rec(n, d, 0, &mut piv, f);
}

fn run_e2(out: &str) {
    let mut rows = vec![];
    for &(n, d, exhaustive) in &[(7usize, 2usize, true), (7, 3, true), (9, 2, true), (9, 3, true), (11, 2, true), (11, 3, true),
                                 (13, 3, false), (13, 4, false), (13, 5, false), (13, 6, false)] {
        let f = Field::new(n);
        let mut hist: BTreeMap<usize, u64> = BTreeMap::new();
        let mut count = 0u64;
        let mut tally = |b: &[u64]| {
            let (r, _) = product_dim(&f, b, b);
            *hist.entry(r).or_default() += 1;
            count += 1;
        };
        if exhaustive {
            rref_subspaces(n, d, &mut tally);
        } else {
            let mut rng = SplitMix(4242 + (n * 100 + d) as u64);
            let mut k = 0;
            while k < 200000 {
                let g: Vec<u64> = (0..d).map(|_| 1 + rng.below((1u64 << n) - 1)).collect();
                let b = subspace(&g);
                if b.len() == d { tally(&b); k += 1; }
            }
        }
        let std_b: Vec<u64> = (0..d).map(|i| 1u64 << i).collect();
        let (rstd, _) = product_dim(&f, &std_b, &std_b);
        let min = *hist.keys().next().unwrap();
        rows.push(obj(vec![
            ("n", jn(n as f64)), ("d", jn(d as f64)), ("exhaustive", V::Bool(exhaustive)), ("subspaces", jn(count as f64)),
            ("min_dim_VV", jn(min as f64)), ("bound_min_n_2d_1", jn(std::cmp::min(n, 2 * d - 1) as f64)),
            ("standard_dim_VV", jn(rstd as f64)), ("generic_d_d1_2", jn(std::cmp::min(n, d * (d + 1) / 2) as f64)),
            ("histogram", V::Obj(hist.iter().map(|(k, v)| (k.to_string(), jn(*v as f64))).collect())),
        ]));
    }
    std::fs::write(out, dump(&V::Arr(rows)) + "\n").unwrap();
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let t0 = Instant::now();
    match a[1].as_str() {
        "replay-grid" => replay_grid(&a[2], a[3].parse().unwrap(), &a[4], a[5].parse().unwrap(), &t0),
        "replay-delta" => replay_delta(&a[2], &a[3], a[4].parse().unwrap(), &t0),
        "null-sym" => { unsafe { SYM_LIN = true; } run_null(&a[2], &a[3], a[4].parse().unwrap(), a[5].parse().unwrap(), &t0) }
        "null" => run_null(&a[2], &a[3], a[4].parse().unwrap(), a[5].parse().unwrap(), &t0),
        "e2" => run_e2(&a[2]),
        "analyze" => analyze::run(&a[2], &a[3], &a[4]),
        "mt" => { let mut r = PyRandom::new(a[2].parse().unwrap()); for _ in 0..5 { print!("{} ", r.randrange(1, 2048)); } println!(); }
        _ => panic!("unknown"),
    }
    eprintln!("{} done in {:.1}s", a[1], t0.elapsed().as_secs_f64());
}

mod analyze;
