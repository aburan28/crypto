//! `volcano-rho` — one fixed ECDLP instance carried down every stratum of
//! an `ℓ`-isogeny volcano, and Pollard rho timed on each curve.
//!
//! 1. Search a random prime-field curve `E/F_p` (`p ≡ 3 mod 4`) whose
//!    Frobenius discriminant `D = t² − 4p = f²·D_K` has `v_ℓ(f) ≥ depth`,
//!    with `#E = h·n`, `n` prime and `h` small.  The `ℓ`-volcano of `E`
//!    then has `v_ℓ(f)` levels below its crater.
//! 2. Walk the `ℓ`-isogeny graph breadth first from `E`.  Each edge is an
//!    explicit Vélu isogeny from a rational root of the `ℓ`-division
//!    polynomial (`ℓ = 2`: the cubic; `ℓ = 3`: `ψ₃`), keyed by `j`.
//!    A vertex's level is `depth − dist(v, floor)`, the floor being the
//!    vertices with exactly one rational `ℓ`-kernel.
//! 3. Fix `P` of order `n` and `Q = kP` on `E`, and push `(P, Q)` through
//!    the edge maps to every curve, so each curve holds the same discrete
//!    log `k` (`ℓ ∤ n`, so every map is injective on `⟨P⟩`).  Every edge is
//!    checked: image on the codomain, `φ(P + Q) = φ(P) + φ(Q)`, `n·φ(P) = O`.
//! 4. Run single-walker distinguished-point Pollard rho (32-adding walk,
//!    no negation map) on every curve with the same trial seeds, plus
//!    isomorphic re-scalings of the root as a representation control, and
//!    record group operations and wall time per solve.
//!
//! Isogenous curves share `#E(F_p)`, so the expected rho cost is the same
//! on every curve; this tool measures how the observed costs spread.

use std::collections::{HashMap, HashSet, VecDeque};
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use crypto_lib::cryptanalysis::isogeny_walk::field::Field;
use crypto_lib::cryptanalysis::isogeny_walk::poly;
use num_bigint::BigUint;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde_json::json;

// ---------------------------------------------------------------- F_p, p < 2^62

#[derive(Clone, Copy)]
struct Fp {
    p: u64,
}

impl Fp {
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p { s - self.p } else { s }
    }
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b { a - b } else { a + self.p - b }
    }
    fn neg(&self, a: u64) -> u64 {
        if a == 0 { 0 } else { self.p - a }
    }
    fn mul(&self, a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % self.p as u128) as u64
    }
    fn sm(&self, k: u64, a: u64) -> u64 {
        self.mul(k % self.p, a)
    }
    fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }
    fn inv(&self, a: u64) -> u64 {
        inv_mod(a, self.p)
    }
    /// `p ≡ 3 mod 4`.
    fn sqrt(&self, a: u64) -> Option<u64> {
        let r = self.pow(a, (self.p + 1) / 4);
        (self.mul(r, r) == a).then_some(r)
    }
}

fn inv_mod(a: u64, m: u64) -> u64 {
    let (mut r0, mut r1) = (m as i128, a as i128);
    let (mut s0, mut s1) = (0i128, 1i128);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (s0, s1) = (s1, s0 - q * s1);
    }
    assert_eq!(r0, 1, "not invertible");
    s0.rem_euclid(m as i128) as u64
}

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

fn powmod(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1 % m;
    a %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, m);
        }
        a = mulmod(a, a, m);
        e >>= 1;
    }
    r
}

/// Deterministic Miller–Rabin for `n < 2^64`.
fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for q in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % q == 0 {
            return n == q;
        }
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    'w: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'w;
            }
        }
        return false;
    }
    true
}

fn factor(mut n: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut d = 2u64;
    while d * d <= n {
        if n % d == 0 {
            let mut e = 0;
            while n % d == 0 {
                n /= d;
                e += 1;
            }
            out.push((d, e));
        }
        d += if d == 2 { 1 } else { 2 };
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

// ---------------------------------------------------------------- curves

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Curve {
    a: u64,
    b: u64,
}

type Pt = Option<(u64, u64)>;

impl Curve {
    fn rhs(&self, f: &Fp, x: u64) -> u64 {
        f.add(f.mul(f.add(f.mul(x, x), self.a), x), self.b)
    }
    fn on(&self, f: &Fp, p: Pt) -> bool {
        p.is_none_or(|(x, y)| f.mul(y, y) == self.rhs(f, x))
    }
    fn j(&self, f: &Fp) -> Option<u64> {
        let a3 = f.sm(4, f.mul(f.mul(self.a, self.a), self.a));
        let den = f.add(a3, f.sm(27, f.mul(self.b, self.b)));
        (den != 0).then(|| f.mul(f.sm(1728, a3), f.inv(den)))
    }
    fn add(&self, f: &Fp, p: Pt, q: Pt) -> Pt {
        let (Some((x1, y1)), Some((x2, y2))) = (p, q) else {
            return p.or(q);
        };
        let l = if x1 == x2 {
            if f.add(y1, y2) == 0 {
                return None;
            }
            f.mul(f.add(f.sm(3, f.mul(x1, x1)), self.a), f.inv(f.add(y1, y1)))
        } else {
            f.mul(f.sub(y2, y1), f.inv(f.sub(x2, x1)))
        };
        let x3 = f.sub(f.sub(f.mul(l, l), x1), x2);
        Some((x3, f.sub(f.mul(l, f.sub(x1, x3)), y1)))
    }
    fn mul(&self, f: &Fp, p: Pt, k: u64) -> Pt {
        let (mut r, mut b, mut k) = (None, p, k);
        while k > 0 {
            if k & 1 == 1 {
                r = self.add(f, r, b);
            }
            b = self.add(f, b, b);
            k >>= 1;
        }
        r
    }
    fn neg(f: &Fp, p: Pt) -> Pt {
        p.map(|(x, y)| (x, f.neg(y)))
    }
    fn random_point(&self, f: &Fp, rng: &mut SmallRng) -> (u64, u64) {
        loop {
            let x = rng.gen_range(0..f.p);
            if let Some(y) = f.sqrt(self.rhs(f, x)) {
                if y != 0 {
                    return (x, y);
                }
            }
        }
    }
}

/// The order of `E(F_p)`: every `M` in the Hasse interval with `M·P = O`
/// by baby-step giant-step, intersected over random points until one is
/// left.  `None` if ten points do not pin it down.
fn group_order(f: &Fp, e: &Curve, rng: &mut SmallRng) -> Option<u64> {
    let s = (4.0 * (f.p as f64).sqrt()).ceil() as u64 + 2;
    let lo = f.p + 1 - s / 2;
    let m = (s as f64).sqrt().ceil() as u64 + 1;
    let mut cands: Option<HashSet<u64>> = None;
    for _ in 0..10 {
        let p = Some(e.random_point(f, rng));
        let mut baby: HashMap<Pt, u64> = HashMap::with_capacity(m as usize);
        let mut jp: Pt = None;
        for j in 0..m {
            baby.entry(Curve::neg(f, jp)).or_insert(j);
            jp = e.add(f, jp, p);
        }
        let step = e.mul(f, p, m);
        let mut g = e.mul(f, p, lo);
        let mut found = HashSet::new();
        for i in 0..=m {
            if let Some(&j) = baby.get(&g) {
                // Every j with −jP = g differs by a multiple of ord(P) < m
                // only when ord(P) is tiny; the intersection handles that.
                found.insert(lo + i * m + j);
            }
            g = e.add(f, g, step);
        }
        cands = Some(match cands {
            None => found,
            Some(c) => c.intersection(&found).copied().collect(),
        });
        let c = cands.as_ref().unwrap();
        let hasse: Vec<u64> = c
            .iter()
            .copied()
            .filter(|&v| (v as i128 - (f.p as i128 + 1)).pow(2) <= 4 * f.p as i128)
            .collect();
        if hasse.len() == 1 {
            return Some(hasse[0]);
        }
    }
    None
}

// ---------------------------------------------------------------- isogenies

/// The rational `x`-coordinates generating `ℓ`-kernels (`ℓ ∈ {2, 3}`).
fn kernel_xs(f: &Fp, field: &Field, e: &Curve, ell: u64) -> Vec<u64> {
    let fe = |v: u64| field.from_u64(v);
    let h = if ell == 2 {
        // x³ + ax + b
        vec![fe(e.b), fe(e.a), fe(0), fe(1)]
    } else {
        // ψ₃ = 3x⁴ + 6ax² + 12bx − a²
        vec![
            fe(f.neg(f.mul(e.a, e.a))),
            fe(f.sm(12, e.b)),
            fe(f.sm(6, e.a)),
            fe(0),
            fe(3),
        ]
    };
    let h = poly::trim(field, h);
    poly::roots(field, &h, 0x5eed)
        .iter()
        .map(|r| field.to_big(r).iter_u64_digits().next().unwrap_or(0))
        .collect()
}

/// Vélu's isogeny with kernel `⟨(x0, ·)⟩` of order `ℓ ∈ {2, 3}`.
#[derive(Clone, Copy, Debug)]
struct Velu {
    ell: u64,
    x0: u64,
    tq: u64,
    uq: u64,
    codomain: Curve,
}

fn velu(f: &Fp, e: &Curve, ell: u64, x0: u64) -> Velu {
    let gx = f.add(f.sm(3, f.mul(x0, x0)), e.a);
    let (tq, uq) = if ell == 2 {
        (gx, 0)
    } else {
        (f.add(gx, gx), f.sm(4, e.rhs(f, x0)))
    };
    let w = f.add(uq, f.mul(x0, tq));
    Velu {
        ell,
        x0,
        tq,
        uq,
        codomain: Curve {
            a: f.sub(e.a, f.sm(5, tq)),
            b: f.sub(e.b, f.sm(7, w)),
        },
    }
}

impl Velu {
    fn map(&self, f: &Fp, p: Pt) -> Pt {
        let (x, y) = p?;
        if x == self.x0 {
            return None;
        }
        let i = f.inv(f.sub(x, self.x0));
        let i2 = f.mul(i, i);
        let i3 = f.mul(i2, i);
        let xx = f.add(x, f.add(f.mul(self.tq, i), f.mul(self.uq, i2)));
        let yy = f.sub(
            y,
            f.mul(y, f.add(f.mul(f.add(self.uq, self.uq), i3), f.mul(self.tq, i2))),
        );
        Some((xx, yy))
    }
}

// ---------------------------------------------------------------- rho

struct Instance {
    curve: Curve,
    p: (u64, u64),
    q: (u64, u64),
}

struct Solve {
    k: Option<u64>,
    iters: u64,
    restarts: u32,
    ns: u128,
}

fn rho(f: &Fp, inst: &Instance, n: u64, seed: u64, dp_bits: u32) -> Solve {
    let e = &inst.curve;
    let (p, q) = (Some(inst.p), Some(inst.q));
    let mut rng = SmallRng::seed_from_u64(seed);
    let t0 = Instant::now();
    let steps: Vec<(u64, u64, (u64, u64))> = (0..32)
        .map(|_| loop {
            let (a, b) = (rng.gen_range(0..n), rng.gen_range(0..n));
            if let Some(r) = e.add(f, e.mul(f, p, a), e.mul(f, q, b)) {
                break (a, b, r);
            }
        })
        .collect();
    let mut table: HashMap<u64, (u64, u64, u64)> = HashMap::new();
    let (mut iters, mut restarts) = (0u64, 0u32);
    let dp_mask = (1u64 << dp_bits) - 1;
    let max_run = 40u64 << dp_bits;
    'restart: loop {
        let (mut a, mut b) = (rng.gen_range(0..n), rng.gen_range(0..n));
        let Some(mut x) = e.add(f, e.mul(f, p, a), e.mul(f, q, b)) else {
            restarts += 1;
            continue;
        };
        let mut run = 0u64;
        loop {
            let s = &steps[(x.0.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 59) as usize];
            iters += 1;
            run += 1;
            a = (a + s.0) % n;
            b = (b + s.1) % n;
            match e.add(f, Some(x), Some(s.2)) {
                Some(nx) => x = nx,
                None => {
                    // aP + bQ = O.
                    if b != 0 {
                        let k = mulmod(n - a % n, inv_mod(b, n), n) % n;
                        return Solve { k: Some(k), iters, restarts, ns: t0.elapsed().as_nanos() };
                    }
                    restarts += 1;
                    continue 'restart;
                }
            }
            if x.0.wrapping_mul(0xC2B2_AE3D_27D4_EB4F) >> 40 & dp_mask == 0 {
                if let Some(&(y2, a2, b2)) = table.get(&x.0) {
                    // Same point: (a − a2)P = (b2 − b)Q; opposite: (a + a2)P = −(b + b2)Q.
                    let (num, den) = if y2 == x.1 {
                        ((a + n - a2) % n, (b2 + n - b) % n)
                    } else {
                        ((n - (a + a2) % n) % n, (b + b2) % n)
                    };
                    if den != 0 {
                        let k = mulmod(num, inv_mod(den, n), n);
                        return Solve { k: Some(k), iters, restarts, ns: t0.elapsed().as_nanos() };
                    }
                    restarts += 1;
                    continue 'restart;
                }
                table.insert(x.0, (x.1, a, b));
                run = 0;
            }
            if run > max_run {
                restarts += 1;
                continue 'restart;
            }
        }
    }
}

// ---------------------------------------------------------------- statistics

fn stats(v: &[f64]) -> serde_json::Value {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len() as f64;
    let mean = s.iter().sum::<f64>() / n;
    let sd = (s.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0).max(1.0)).sqrt();
    let q = |p: f64| s[((p * (n - 1.0)).round() as usize).min(s.len() - 1)];
    json!({
        "n": s.len(), "min": s[0], "max": s[s.len() - 1], "mean": mean,
        "median": q(0.5), "p10": q(0.1), "p90": q(0.9), "sd": sd,
        "ci95_half": 1.96 * sd / n.sqrt(),
    })
}

// ---------------------------------------------------------------- main

struct Args {
    bits: u32,
    ell: u64,
    depth: u32,
    seed: u64,
    trials: u64,
    per_level: usize,
    dp_bits: u32,
    max_cofactor: u64,
    controls: usize,
    out: PathBuf,
}

fn parse_args() -> Args {
    let mut a = Args {
        bits: 40,
        ell: 3,
        depth: 3,
        seed: 1,
        trials: 100,
        per_level: 10,
        dp_bits: 8,
        max_cofactor: 64,
        controls: 3,
        out: PathBuf::from("volcano-rho-out"),
    };
    let v: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < v.len() {
        let val = v.get(i + 1).cloned().unwrap_or_default();
        match v[i].as_str() {
            "--bits" => a.bits = val.parse().unwrap(),
            "--ell" => a.ell = val.parse().unwrap(),
            "--depth" => a.depth = val.parse().unwrap(),
            "--seed" => a.seed = val.parse().unwrap(),
            "--trials" => a.trials = val.parse().unwrap(),
            "--per-level" => a.per_level = val.parse().unwrap(),
            "--dp-bits" => a.dp_bits = val.parse().unwrap(),
            "--max-cofactor" => a.max_cofactor = val.parse().unwrap(),
            "--controls" => a.controls = val.parse().unwrap(),
            "--out" => a.out = PathBuf::from(val),
            "-h" | "--help" => {
                println!(
                    "volcano-rho [--bits 40] [--ell 3] [--depth 3] [--seed 1] [--trials 100] \
                     [--per-level 10] [--dp-bits 8] [--max-cofactor 64] [--controls 3] [--out DIR]"
                );
                std::process::exit(0);
            }
            o => panic!("unknown argument {o}"),
        }
        i += 2;
    }
    assert!(a.ell == 2 || a.ell == 3, "--ell must be 2 or 3");
    assert!(a.depth >= 1 && (20..=60).contains(&a.bits));
    a
}

fn vl(mut x: u128, l: u128) -> u32 {
    let mut v = 0;
    while x > 0 && x % l == 0 {
        x /= l;
        v += 1;
    }
    v
}

struct Vertex {
    j: u64,
    curve: Curve,
    dist: u32,
    parent: Option<usize>,
    edge: Option<Velu>,
    kernels: Vec<u64>,
    nbrs: Vec<usize>,
    p: (u64, u64),
    q: (u64, u64),
}


struct Found {
    f: Fp,
    root: Curve,
    order: u64,
    n: u64,
    t: i64,
    d_abs: u64,
    depth: u32,
    tried: u64,
}

/// Random curves over random primes until one has an `ℓ`-volcano of the
/// requested depth and `#E = h·n` with `n` prime, `h ≤ max_cofactor`.
fn search(args: &Args, rng: &mut SmallRng) -> Found {
    let mut tried = 0u64;
    loop {
        let p = loop {
            let c = rng.gen_range((1u64 << (args.bits - 1))..(1u64 << args.bits)) | 3;
            if is_prime(c) {
                break c;
            }
        };
        let f = Fp { p };
        for _ in 0..200 {
            tried += 1;
            let e = Curve { a: rng.gen_range(1..p), b: rng.gen_range(1..p) };
            match e.j(&f) {
                Some(j) if j != 0 && j != 1728 % p => {}
                _ => continue,
            }
            let Some(order) = group_order(&f, &e, rng) else { continue };
            let t = p as i64 + 1 - order as i64;
            let d_abs = (4 * p as i128 - (t as i128) * (t as i128)) as u64;
            let depth = vl(d_abs as u128, args.ell as u128) / 2;
            if depth != args.depth {
                continue;
            }
            let n = factor(order).last().unwrap().0;
            if order / n > args.max_cofactor || n % args.ell == 0 || n < (1 << 20) {
                continue;
            }
            return Found { f, root: e, order, n, t, d_abs, depth, tried };
        }
    }
}

/// `(D_K, f)` with `−d_abs = f²·D_K`, `D_K` fundamental.
fn fundamental(d_abs: u64) -> (i64, u64) {
    let (mut core, mut f) = (1u64, 1u64);
    for (q, e) in factor(d_abs) {
        core *= q.pow(e % 2);
        f *= q.pow(e / 2);
    }
    let d0 = -(core as i64);
    if d0.rem_euclid(4) == 1 {
        (d0, f)
    } else {
        (4 * d0, f / 2)
    }
}

fn crater_kind(dk: i64, ell: u64) -> &'static str {
    let r = if ell == 2 {
        match dk.rem_euclid(8) {
            1 => 1,
            5 => -1,
            _ => 0,
        }
    } else {
        match dk.rem_euclid(3) {
            1 => 1,
            2 => -1,
            _ => 0,
        }
    };
    match r {
        1 => "split (crater is a cycle; two horizontal edges per crater vertex)",
        -1 => "inert (crater is one vertex; no horizontal edges)",
        _ => "ramified (one horizontal edge per crater vertex)",
    }
}

fn main() {
    let args = parse_args();
    let mut rng = SmallRng::seed_from_u64(args.seed);
    let ell = args.ell;
    let t_search = Instant::now();
    let found = search(&args, &mut rng);
    let search_s = t_search.elapsed().as_secs_f64();
    let (f, n, depth) = (found.f, found.n, found.depth);
    let (dk, cond) = fundamental(found.d_abs);
    let field = Field::new(&BigUint::from(f.p)).expect("p prime");
    eprintln!(
        "curve: p={} a={} b={} #E={} n={} h={} t={} D=-{} = {}^2*{}  depth={} ({} curves tried, {:.1}s)",
        f.p, found.root.a, found.root.b, found.order, n, found.order / n, found.t,
        found.d_abs, cond, dk, depth, found.tried, search_s
    );

    // The instance: P of order n, Q = kP.
    let e0 = found.root;
    let cof = found.order / n;
    let p0 = loop {
        if let Some(pt) = e0.mul(&f, Some(e0.random_point(&f, &mut rng)), cof) {
            break pt;
        }
    };
    assert!(e0.mul(&f, Some(p0), n).is_none());
    let k_secret = rng.gen_range(1..n);
    let q0 = e0.mul(&f, Some(p0), k_secret).unwrap();

    // 2. Breadth-first walk with explicit Vélu maps.
    let r_meas = 2 * depth + 2;
    let r_explore = r_meas + depth;
    let mut vs: Vec<Vertex> = vec![Vertex {
        j: e0.j(&f).unwrap(),
        curve: e0,
        dist: 0,
        parent: None,
        edge: None,
        kernels: Vec::new(),
        nbrs: Vec::new(),
        p: p0,
        q: q0,
    }];
    let mut by_j: HashMap<u64, usize> = HashMap::from([(vs[0].j, 0)]);
    let mut queue = VecDeque::from([0usize]);
    let mut edges_checked = 0u64;
    let mut skipped_special = 0u64;
    while let Some(v) = queue.pop_front() {
        let e = vs[v].curve;
        let kernels = kernel_xs(&f, &field, &e, ell);
        vs[v].kernels = kernels.clone();
        for x0 in kernels {
            let iso = velu(&f, &e, ell, x0);
            let c = iso.codomain;
            let Some(j2) = c.j(&f) else { continue };
            if j2 == 0 || j2 == 1728 % f.p {
                skipped_special += 1;
                continue;
            }
            let w = match by_j.get(&j2) {
                Some(&w) => w,
                None => {
                    if vs[v].dist >= r_explore || vs.len() >= 20_000 {
                        continue;
                    }
                    let (pp, qq) = (Some(vs[v].p), Some(vs[v].q));
                    let (ip, iq) = (iso.map(&f, pp), iso.map(&f, qq));
                    let isum = iso.map(&f, e.add(&f, pp, qq));
                    assert!(c.on(&f, ip) && c.on(&f, iq), "image off the codomain");
                    assert_eq!(c.add(&f, ip, iq), isum, "Vélu map is not a homomorphism");
                    assert!(ip.is_some() && c.mul(&f, ip, n).is_none(), "order not preserved");
                    edges_checked += 1;
                    let w = vs.len();
                    vs.push(Vertex {
                        j: j2,
                        curve: c,
                        dist: vs[v].dist + 1,
                        parent: Some(v),
                        edge: Some(iso),
                        kernels: Vec::new(),
                        nbrs: Vec::new(),
                        p: ip.unwrap(),
                        q: iq.unwrap(),
                    });
                    by_j.insert(j2, w);
                    queue.push_back(w);
                    w
                }
            };
            if w != v && !vs[v].nbrs.contains(&w) {
                vs[v].nbrs.push(w);
            }
            if w != v && !vs[w].nbrs.contains(&v) {
                vs[w].nbrs.push(v);
            }
        }
    }

    // Levels: depth − distance to the floor (vertices with one ℓ-kernel).
    let level_of = |v: usize| -> Option<u32> {
        let mut seen = HashMap::from([(v, 0u32)]);
        let mut q = VecDeque::from([v]);
        while let Some(u) = q.pop_front() {
            let du = seen[&u];
            if vs[u].kernels.len() == 1 {
                return Some(depth.saturating_sub(du));
            }
            if du == depth {
                continue;
            }
            for &w in &vs[u].nbrs {
                if let std::collections::hash_map::Entry::Vacant(s) = seen.entry(w) {
                    s.insert(du + 1);
                    q.push_back(w);
                }
            }
        }
        None
    };
    let levels: Vec<Option<u32>> =
        (0..vs.len()).map(|v| (vs[v].dist <= r_meas).then(|| level_of(v)).flatten()).collect();
    for (v, lv) in levels.iter().enumerate() {
        let Some(lv) = lv else { continue };
        for &w in &vs[v].nbrs {
            if let Some(lw) = levels[w] {
                assert!(lv.abs_diff(lw) <= 1, "edge skips a level");
            }
        }
    }
    let mut per_level: Vec<Vec<usize>> = vec![Vec::new(); depth as usize + 1];
    for (v, lv) in levels.iter().enumerate() {
        if let Some(lv) = lv {
            per_level[*lv as usize].push(v);
        }
    }
    eprintln!(
        "walk: {} curves, {} verified edges; per level (0 = crater): {:?}",
        vs.len(),
        edges_checked,
        per_level.iter().map(|l| l.len()).collect::<Vec<_>>()
    );

    // 3. Curves to measure: the root, then up to `per_level` per level in BFS order.
    let mut measured: Vec<usize> = vec![0];
    for l in &per_level {
        for &v in l.iter().filter(|&&v| v != 0).take(args.per_level) {
            measured.push(v);
        }
    }
    measured.sort();
    struct Arm {
        label: String,
        vertex: Option<usize>,
        level: Option<u32>,
        inst: Instance,
    }
    let mut arms: Vec<Arm> = measured
        .iter()
        .map(|&v| Arm {
            label: format!("v{v}"),
            vertex: Some(v),
            level: levels[v],
            inst: Instance { curve: vs[v].curve, p: vs[v].p, q: vs[v].q },
        })
        .collect();
    for c in 0..args.controls {
        // Isomorphic model (u⁴a, u⁶b), (x, y) ↦ (u²x, u³y): same curve, new coordinates.
        let u = rng.gen_range(2..f.p);
        let (u2, u3) = (f.mul(u, u), f.mul(f.mul(u, u), u));
        let cv = Curve { a: f.mul(f.mul(u2, u2), e0.a), b: f.mul(f.mul(u3, u3), e0.b) };
        let mv = |pt: (u64, u64)| (f.mul(u2, pt.0), f.mul(u3, pt.1));
        let (cp, cq) = (mv(p0), mv(q0));
        assert!(cv.on(&f, Some(cp)) && cv.on(&f, Some(cq)));
        arms.push(Arm {
            label: format!("iso{c}"),
            vertex: None,
            level: levels[0],
            inst: Instance { curve: cv, p: cp, q: cq },
        });
    }

    // 4. Pollard rho: every arm × the same trial seeds, interleaved.
    let mut jobs: Vec<(usize, u64)> =
        (0..arms.len()).flat_map(|a| (0..args.trials).map(move |s| (a, s))).collect();
    jobs.shuffle(&mut rng);
    let seed_base = args.seed.wrapping_mul(0x1000_0000_01b3);
    eprintln!("rho: {} arms x {} trials = {} solves", arms.len(), args.trials, jobs.len());
    let t_rho = Instant::now();
    let results: Vec<(usize, u64, Solve)> = jobs
        .par_iter()
        .map(|&(a, s)| (a, s, rho(&f, &arms[a].inst, n, seed_base ^ s, args.dp_bits)))
        .collect();
    let rho_s = t_rho.elapsed().as_secs_f64();
    let wrong = results.iter().filter(|r| r.2.k != Some(k_secret)).count();
    assert_eq!(wrong, 0, "a rho solve returned the wrong logarithm");

    // 5. Records and summaries.
    fs::create_dir_all(&args.out).unwrap();
    let mut csv = String::from("arm,vertex,level,seed,iters,restarts,ns\n");
    let mut by_arm: Vec<Vec<&Solve>> = vec![Vec::new(); arms.len()];
    let mut sorted: Vec<&(usize, u64, Solve)> = results.iter().collect();
    sorted.sort_by_key(|r| (r.0, r.1));
    for (a, s, r) in &sorted {
        by_arm[*a].push(r);
        csv += &format!(
            "{},{},{},{},{},{},{}\n",
            arms[*a].label,
            arms[*a].vertex.map_or(String::new(), |v| v.to_string()),
            arms[*a].level.map_or(String::new(), |l| l.to_string()),
            seed_base ^ s,
            r.iters,
            r.restarts,
            r.ns
        );
    }
    fs::write(args.out.join("trials.csv"), csv).unwrap();

    let summarize = |rs: &[&Solve]| {
        let it: Vec<f64> = rs.iter().map(|r| r.iters as f64).collect();
        let ms: Vec<f64> = rs.iter().map(|r| r.ns as f64 / 1e6).collect();
        let nsop: Vec<f64> = rs.iter().map(|r| r.ns as f64 / r.iters as f64).collect();
        json!({"iters": stats(&it), "ms": stats(&ms), "ns_per_op": stats(&nsop)})
    };
    let arm_json: Vec<serde_json::Value> = arms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            json!({
                "arm": a.label, "vertex": a.vertex, "level": a.level,
                "kind": if a.vertex.is_some() { "isogenous" } else { "isomorphic-control" },
                "curve": {"a": a.inst.curve.a, "b": a.inst.curve.b, "j": a.inst.curve.j(&f)},
                "P": [a.inst.p.0, a.inst.p.1], "Q": [a.inst.q.0, a.inst.q.1],
                "summary": summarize(&by_arm[i]),
            })
        })
        .collect();
    let mut level_json = Vec::new();
    for l in 0..=depth {
        let rs: Vec<&Solve> = arms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.vertex.is_some() && a.level == Some(l))
            .flat_map(|(i, _)| by_arm[i].iter().copied())
            .collect();
        if !rs.is_empty() {
            let curves = arms.iter().filter(|a| a.vertex.is_some() && a.level == Some(l)).count();
            level_json.push(json!({"level": l, "curves": curves, "summary": summarize(&rs)}));
        }
    }
    let all_iso: Vec<&Solve> = arms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.vertex.is_some())
        .flat_map(|(i, _)| by_arm[i].iter().copied())
        .collect();
    let ctrl: Vec<&Solve> = arms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.vertex.is_none())
        .flat_map(|(i, _)| by_arm[i].iter().copied())
        .collect();
    let vertices: Vec<serde_json::Value> = vs
        .iter()
        .enumerate()
        .map(|(i, v)| {
            json!({
                "id": i, "j": v.j, "a": v.curve.a, "b": v.curve.b, "dist": v.dist,
                "level": levels[i], "kernels": v.kernels.len(), "nbrs": v.nbrs,
                "parent": v.parent,
                "edge": v.edge.map(|e| json!({"ell": e.ell, "kernel_x": e.x0})),
            })
        })
        .collect();
    let expected = (std::f64::consts::PI * n as f64 / 2.0).sqrt();
    let host = fs::read_to_string("/proc/cpuinfo")
        .ok()
        .and_then(|c| c.lines().find(|l| l.starts_with("model name")).map(|l| l.split(':').nth(1).unwrap().trim().to_string()));
    let out = json!({
        "schema": "volcano-rho/v1",
        "command": std::env::args().collect::<Vec<_>>().join(" "),
        "params": {"bits": args.bits, "ell": ell, "depth": args.depth, "seed": args.seed,
                   "trials": args.trials, "per_level": args.per_level, "dp_bits": args.dp_bits,
                   "max_cofactor": args.max_cofactor, "controls": args.controls},
        "host": {"cpu": host, "threads": rayon::current_num_threads()},
        "problem": {"p": f.p, "root": {"a": e0.a, "b": e0.b}, "order": found.order, "n": n,
                    "cofactor": found.order / n, "trace": found.t, "disc": -(found.d_abs as i64),
                    "fundamental_disc": dk, "conductor": cond, "volcano_depth": depth,
                    "crater": crater_kind(dk, ell), "P": [p0.0, p0.1], "Q": [q0.0, q0.1], "k": k_secret,
                    "search_curves_tried": found.tried, "search_s": search_s},
        "walk": {"radius_measured": r_meas, "radius_explored": r_explore, "vertices": vs.len(),
                 "edges_verified": edges_checked, "special_j_skipped": skipped_special,
                 "per_level_found": per_level.iter().map(|l| l.len()).collect::<Vec<_>>()},
        "rho": {"walk": "32-adding (Teske), single walker, distinguished points, no negation map",
                "expected_iters_sqrt_pi_n_over_2": expected, "solves": results.len(),
                "all_correct": wrong == 0, "wall_s": rho_s},
        "arms": arm_json,
        "levels": level_json,
        "all_isogenous": summarize(&all_iso),
        "controls": if ctrl.is_empty() { serde_json::Value::Null } else { summarize(&ctrl) },
        "vertices": vertices,
    });
    fs::write(args.out.join("results.json"), serde_json::to_string_pretty(&out).unwrap()).unwrap();

    println!("\nPollard rho on one instance (k = {k_secret}, n = {n}) across the {ell}-volcano");
    println!("expected iterations ~ sqrt(pi n / 2) = {expected:.0} (+ DP tail ~2^{})", args.dp_bits);
    println!("{:<8}{:>6}{:>12}{:>12}{:>12}{:>10}{:>10}{:>10}{:>9}", "arm", "level", "min it", "mean it", "max it", "min ms", "mean ms", "max ms", "ns/op");
    let row = |name: &str, lv: String, s: &serde_json::Value| {
        println!(
            "{:<8}{:>6}{:>12.0}{:>12.0}{:>12.0}{:>10.2}{:>10.2}{:>10.2}{:>9.1}",
            name, lv, s["iters"]["min"].as_f64().unwrap(), s["iters"]["mean"].as_f64().unwrap(),
            s["iters"]["max"].as_f64().unwrap(), s["ms"]["min"].as_f64().unwrap(),
            s["ms"]["mean"].as_f64().unwrap(), s["ms"]["max"].as_f64().unwrap(),
            s["ns_per_op"]["mean"].as_f64().unwrap()
        );
    };
    for a in &arm_json {
        row(a["arm"].as_str().unwrap(), a["level"].to_string(), &a["summary"]);
    }
    println!("-- pooled by level --");
    for l in &level_json {
        row(&format!("L{}", l["level"]), format!("{}c", l["curves"]), &l["summary"]);
    }
    row("ALL", "iso".into(), &out["all_isogenous"]);
    if !ctrl.is_empty() {
        row("CTRL", "iso~".into(), &out["controls"]);
    }
    eprintln!("wrote {}", args.out.display());
}
