//! Two probes on the question "where could a prime-field ECDLP
//! breakthrough come from?", run on the native engine in
//! `cryptanalysis::prime_fast`.
//!
//! **A. Preprocessing rho (Bernstein–Lange 2013) on prime-field curves.**
//! The only known way to get sub-√n *per target* on a generic group:
//! a target-independent walk (table points are multiples of `G` only)
//! precomputes `S` distinguished points; each online target walks until
//! it hits one. With `S ≈ n^{1/3}` and `2^dp ≈ √(n/S)` the online cost is
//! `≈ 2√(n/S) = 2n^{1/3}` steps against `√(πn/4)` for plain rho. This is
//! non-uniform (the precomputation is charged once per curve) and sits
//! exactly on the generic `S·T² ≈ n` line, so it is a *control*, not a
//! breakthrough; the repo never ran it on prime curves.
//!
//! **B. Lattice (Coppersmith / Jochemsz–May) decomposition of `S₃`.**
//! The one non-generic lever over prime fields: recover `(x_i, x_j)`
//! with `x_i, x_j ≤ B = p^δ` from `S₃(x_i, x_j, x_R) ≡ 0 (mod p)` by
//! small-root lattice reduction instead of enumerating the factor base.
//! For index calculus to beat rho this needs `δ ≥ 1/4` (so that the
//! `p/(2B²)` trials per relation fall below `√p`). The basic lattice for
//! the `[0,2]²` monomial rectangle predicts `δ < 1/18, 1/10, 0.119`
//! for `m = 1, 2, 3` levels and `1/6` asymptotically. The probe plants
//! roots and measures the empirical success threshold.
//!
//! ```bash
//! cargo run --release --example prime_ecdlp_breakthrough_probe -- [--probe a|b|c|abc]
//!     [--rho-bits 32,36,40,44] [--lattice-bits 40,56] [--targets 16] [--json out.json]
//! ```
//!
//! Public synthetic known-answer instances only; every log recovered in A is
//! checked against the planted scalar, every root in B against the planted root.

use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Mutex;
use std::time::Instant;

use crypto_lib::cryptanalysis::lattice::lll_reduce;
use crypto_lib::cryptanalysis::prime_fast::{
    find_a3_curve, rho_parallel, FastCurve, Pt, RhoConfig,
};
use num_bigint::{BigInt, BigUint};
use num_traits::{One, ToPrimitive, Zero};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;

// ── shared walk machinery (same hashing as prime_fast::rho_parallel) ─

#[inline(always)]
fn branch_index(x: u64, table_bits: u32) -> usize {
    (x.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - table_bits)) as usize
}

#[inline(always)]
fn is_distinguished(x: u64, dp_bits: u32) -> bool {
    if dp_bits == 0 {
        return true;
    }
    (x.rotate_left(23).wrapping_mul(0xD6E8_FEB8_6659_FD93) >> (64 - dp_bits)) == 0
}

#[inline(always)]
fn add_mod(a: u64, b: u64, n: u64) -> u64 {
    let s = a + b;
    if s >= n {
        s - n
    } else {
        s
    }
}

#[inline(always)]
fn neg_mod(a: u64, n: u64) -> u64 {
    if a == 0 {
        0
    } else {
        n - a
    }
}

fn inv_mod(a: u64, n: u64) -> u64 {
    // n prime; a != 0
    let (mut old_r, mut r) = (a as i128, n as i128);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        let t = old_r - q * r;
        old_r = r;
        r = t;
        let t = old_s - q * s;
        old_s = s;
        s = t;
    }
    old_s.rem_euclid(n as i128) as u64
}

/// Target-independent walk: table entries are `[u_j]G` only.
struct GWalk {
    table: Vec<(Pt, u64)>,
    table_bits: u32,
    dp_bits: u32,
}

impl GWalk {
    fn new(curve: &FastCurve, g: Pt, table_bits: u32, dp_bits: u32, seed: u64) -> Self {
        let mut rng = StdRng::seed_from_u64(seed ^ 0x5151_1234_ABCD_0001);
        let mut table = Vec::with_capacity(1 << table_bits);
        while table.len() < (1 << table_bits) {
            let u = rng.gen_range(1..curve.n);
            if let Some(pt) = curve.mul(Some(g), u) {
                table.push((pt, u));
            }
        }
        GWalk {
            table,
            table_bits,
            dp_bits,
        }
    }
}

// ── Probe A: preprocessing rho ────────────────────────────────────────

struct PreTable {
    /// x (Montgomery) → a with `[a]G` = the canonical point.
    map: HashMap<u64, u64>,
    precompute_steps: u64,
    precompute_wall: f64,
}

/// Precompute `S = 2^table_entries_bits` distinguished points of `G`-only walks.
fn precompute(
    curve: &FastCurve,
    g: Pt,
    walk: &GWalk,
    entries_bits: u32,
    threads: usize,
    walkers: usize,
    seed: u64,
) -> PreTable {
    let t0 = Instant::now();
    let want = 1usize << entries_bits;
    let map: Mutex<HashMap<u64, u64>> = Mutex::new(HashMap::with_capacity(want));
    let done = AtomicBool::new(false);
    let steps = AtomicU64::new(0);
    let f = curve.f;
    let n = curve.n;
    std::thread::scope(|s| {
        for tid in 0..threads {
            let map = &map;
            let done = &done;
            let steps = &steps;
            s.spawn(move || {
                let mut rng =
                    StdRng::seed_from_u64(seed ^ ((tid as u64 + 11) * 0x9E37_79B9_7F4A_7C15));
                let w = walkers;
                let cap: u32 = (24u64 << walk.dp_bits).max(256) as u32;
                let mut xs = vec![0u64; w];
                let mut ys = vec![0u64; w];
                let mut as_ = vec![0u64; w];
                let mut prev = vec![u64::MAX; w];
                let mut cnt = vec![0u32; w];
                let mut d = vec![0u64; w];
                let mut js = vec![0usize; w];
                let mut reset = vec![false; w];
                let mut scratch = Vec::with_capacity(w);
                let restart = |i: usize,
                               rng: &mut StdRng,
                               xs: &mut [u64],
                               ys: &mut [u64],
                               as_: &mut [u64],
                               prev: &mut [u64],
                               cnt: &mut [u32]| loop {
                    let a = rng.gen_range(1..n);
                    if let Some(pt) = curve.mul(Some(g), a) {
                        let (pt, flipped) = curve.canonical_sign(pt);
                        xs[i] = pt.x;
                        ys[i] = pt.y;
                        as_[i] = if flipped { neg_mod(a, n) } else { a };
                        prev[i] = u64::MAX;
                        cnt[i] = 0;
                        return;
                    }
                };
                for i in 0..w {
                    restart(i, &mut rng, &mut xs, &mut ys, &mut as_, &mut prev, &mut cnt);
                }
                let mut local = 0u64;
                'outer: loop {
                    if done.load(Ordering::Relaxed) {
                        break;
                    }
                    for i in 0..w {
                        let j = branch_index(xs[i], walk.table_bits);
                        js[i] = j;
                        let dd = f.sub(walk.table[j].0.x, xs[i]);
                        reset[i] = dd == 0;
                        d[i] = if dd == 0 { f.one } else { dd };
                    }
                    f.batch_inv(&mut d, &mut scratch);
                    for i in 0..w {
                        if reset[i] {
                            restart(i, &mut rng, &mut xs, &mut ys, &mut as_, &mut prev, &mut cnt);
                            continue;
                        }
                        let (t, u) = walk.table[js[i]];
                        let lam = f.mul(f.sub(t.y, ys[i]), d[i]);
                        let x3 = f.sub(f.sub(f.sqr(lam), xs[i]), t.x);
                        let mut y3 = f.sub(f.mul(lam, f.sub(xs[i], x3)), ys[i]);
                        let mut a = add_mod(as_[i], u, n);
                        if y3 > f.half {
                            y3 = f.p - y3;
                            a = neg_mod(a, n);
                        }
                        if x3 == prev[i] {
                            let (cx, cy, ca) = if xs[i] < x3 {
                                (xs[i], ys[i], as_[i])
                            } else {
                                (x3, y3, a)
                            };
                            match curve.double(Pt { x: cx, y: cy }) {
                                Some(pt) => {
                                    let (pt, flipped) = curve.canonical_sign(pt);
                                    let a2 = add_mod(ca, ca, n);
                                    xs[i] = pt.x;
                                    ys[i] = pt.y;
                                    as_[i] = if flipped { neg_mod(a2, n) } else { a2 };
                                    prev[i] = u64::MAX;
                                }
                                None => {
                                    restart(
                                        i, &mut rng, &mut xs, &mut ys, &mut as_, &mut prev,
                                        &mut cnt,
                                    );
                                    continue;
                                }
                            }
                        } else {
                            prev[i] = xs[i];
                            xs[i] = x3;
                            ys[i] = y3;
                            as_[i] = a;
                        }
                        cnt[i] += 1;
                        if is_distinguished(xs[i], walk.dp_bits) {
                            let mut m = map.lock().unwrap();
                            if m.len() >= want {
                                done.store(true, Ordering::Relaxed);
                                drop(m);
                                break 'outer;
                            }
                            m.entry(xs[i]).or_insert(as_[i]);
                            drop(m);
                            restart(i, &mut rng, &mut xs, &mut ys, &mut as_, &mut prev, &mut cnt);
                        } else if cnt[i] > cap {
                            restart(i, &mut rng, &mut xs, &mut ys, &mut as_, &mut prev, &mut cnt);
                        }
                    }
                    local += w as u64;
                    if local >= 1 << 16 {
                        steps.fetch_add(local, Ordering::Relaxed);
                        local = 0;
                    }
                }
                steps.fetch_add(local, Ordering::Relaxed);
            });
        }
    });
    PreTable {
        map: map.into_inner().unwrap(),
        precompute_steps: steps.load(Ordering::Relaxed),
        precompute_wall: t0.elapsed().as_secs_f64(),
    }
}

/// Online phase for one target: walk from `[a₀]G + Q` with full `(a, b)`
/// coefficients (`a·G + b·Q`, both mod n) and *the same* fruitless-cycle
/// escape as the precomputation, so an online trail that merges into a
/// precomputed trail follows it all the way to its distinguished point.
fn online(
    curve: &FastCurve,
    g: Pt,
    q: Pt,
    walk: &GWalk,
    table: &PreTable,
    seed: u64,
    max_steps: u64,
) -> Option<(u64, u64)> {
    let f = curve.f;
    let n = curve.n;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x0F0F_1234_5678_9ABC);
    let cap = (24u64 << walk.dp_bits).max(256);
    let mut total = 0u64;
    while total < max_steps {
        let a0 = rng.gen_range(1..n);
        let Some(start) = curve.add(curve.mul(Some(g), a0), Some(q)) else {
            continue;
        };
        let (mut cur, flipped) = curve.canonical_sign(start);
        let mut a = if flipped { neg_mod(a0, n) } else { a0 };
        let mut b = if flipped { n - 1 } else { 1 };
        let mut prev = u64::MAX;
        let mut cnt = 0u64;
        loop {
            if is_distinguished(cur.x, walk.dp_bits) {
                if let Some(&at) = table.map.get(&cur.x) {
                    // a_t G = a G + b Q  ⟹  k = (a_t − a) / b
                    let num = if at >= a { at - a } else { at + n - a };
                    if b != 0 {
                        let k = ((num as u128 * inv_mod(b, n) as u128) % n as u128) as u64;
                        if curve.mul(Some(g), k) == Some(q) {
                            return Some((k, total));
                        }
                    }
                }
                break; // DP not in table: restart from a fresh point
            }
            if cnt > cap {
                break;
            }
            let j = branch_index(cur.x, walk.table_bits);
            let (t, u) = walk.table[j];
            let dd = f.sub(t.x, cur.x);
            if dd == 0 {
                break;
            }
            let inv = f.inv(dd);
            let nxt = curve.add_with_inv(cur, t, inv);
            let mut na = add_mod(a, u, n);
            let mut nb = b;
            let (mut np, fl) = curve.canonical_sign(nxt);
            if fl {
                na = neg_mod(na, n);
                nb = neg_mod(nb, n);
            }
            if np.x == prev {
                // Same deterministic escape as `precompute`: double the member with smaller x.
                let (cp, ca, cb) = if cur.x < np.x {
                    (cur, a, b)
                } else {
                    (np, na, nb)
                };
                match curve.double(cp) {
                    Some(pt) => {
                        let (pt, fl2) = curve.canonical_sign(pt);
                        let a2 = add_mod(ca, ca, n);
                        let b2 = add_mod(cb, cb, n);
                        np = pt;
                        na = if fl2 { neg_mod(a2, n) } else { a2 };
                        nb = if fl2 { neg_mod(b2, n) } else { b2 };
                        prev = u64::MAX;
                    }
                    None => break,
                }
            } else {
                prev = cur.x;
            }
            cur = np;
            a = na;
            b = nb;
            cnt += 1;
            total += 1;
        }
    }
    None
}

fn probe_a(bits_list: &[u32], targets: usize, out: &mut Vec<serde_json::Value>) {
    println!("## A. Preprocessing rho on prime-field curves (native engine)\n");
    println!("`S = 2^s` table entries, `2^dp ≈ √(n/S)`, online walks single-threaded; plain rho = `rho_parallel` on the same targets.\n");
    println!("| bits | S (2^s) | dp | precompute steps | precompute wall s | online steps mean (targets) | 2√(n/S) | plain rho steps √(πn/4) | per-target speed-up (steps) | online wall ms mean | plain rho wall s | break-even targets |");
    println!("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|");
    let avail = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    for &bits in bits_list {
        let Some((curve, g)) = find_a3_curve(bits, 3) else {
            continue;
        };
        let n = curve.n as f64;
        let s_bits = ((n.ln() / 3.0) / 2f64.ln()).round() as u32; // S ≈ n^{1/3}
        let s = (1u64 << s_bits) as f64;
        let dp_bits = ((n / s).sqrt().log2()).round().max(1.0) as u32;
        let walk = GWalk::new(&curve, g, 10, dp_bits, 77);
        let table = precompute(&curve, g, &walk, s_bits, avail, 512, 99);
        let mut rng = StdRng::seed_from_u64(20_261_008 + bits as u64);
        let mut online_steps = Vec::new();
        let mut online_walls = Vec::new();
        let mut rho_steps = Vec::new();
        let mut rho_walls = Vec::new();
        let mut solved = 0usize;
        for t in 0..targets {
            let k = rng.gen_range(1..curve.n);
            let q = curve.mul(Some(g), k).unwrap();
            let t0 = Instant::now();
            let res = online(&curve, g, q, &walk, &table, 1000 + t as u64, 1u64 << 40);
            let wall = t0.elapsed().as_secs_f64();
            match res {
                Some((kk, st)) => {
                    assert_eq!(kk, k, "preprocessing rho wrong answer");
                    solved += 1;
                    online_steps.push(st as f64);
                    online_walls.push(wall);
                }
                None => eprintln!("online failed at {bits} bits target {t}"),
            }
            if t < 3 {
                let mut cfg = RhoConfig::for_bits(bits, true);
                cfg.seed = 5 + t as u64;
                let r = rho_parallel(&curve, g, q, &cfg).unwrap();
                assert_eq!(r.log, k);
                rho_steps.push(r.total_steps as f64);
                rho_walls.push(r.wall_secs);
            }
        }
        let mean = |v: &[f64]| {
            if v.is_empty() {
                0.0
            } else {
                v.iter().sum::<f64>() / v.len() as f64
            }
        };
        let on_mean = mean(&online_steps);
        let rho_mean = mean(&rho_steps);
        let rho_theory = (std::f64::consts::PI * n / 4.0).sqrt();
        let on_theory = 2.0 * (n / s).sqrt();
        let on_wall = mean(&online_walls);
        let rho_wall = mean(&rho_walls);
        let breakeven = if rho_wall > on_wall {
            table.precompute_wall / (rho_wall - on_wall)
        } else {
            f64::INFINITY
        };
        println!(
            "| {bits} | 2^{s_bits} | {dp_bits} | {:.3e} | {:.1} | {:.3e} ({solved}/{targets}) | {:.3e} | {:.3e} (measured {:.3e}) | {:.0}× | {:.2} | {:.3} | {:.0} |",
            table.precompute_steps as f64,
            table.precompute_wall,
            on_mean,
            on_theory,
            rho_theory,
            rho_mean,
            rho_theory / on_mean.max(1.0),
            on_wall * 1e3,
            rho_wall,
            breakeven
        );
        out.push(json!({"probe":"A","bits":bits,"p":curve.f.p,"n":curve.n,"table_bits":s_bits,"dp_bits":dp_bits,
            "precompute_steps":table.precompute_steps,"precompute_wall_s":table.precompute_wall,"table_entries":table.map.len(),
            "online_steps_mean":on_mean,"online_theory_2sqrt_n_over_s":on_theory,"online_wall_s_mean":on_wall,
            "targets":targets,"solved":solved,"rho_steps_theory":rho_theory,"rho_steps_measured_mean":rho_mean,
            "rho_wall_s_mean":rho_wall,"breakeven_targets":breakeven}));
    }
    println!();
}

// ── Probe B: lattice decomposition of S₃ ──────────────────────────────

/// Coefficients of `S₃(u, v, z)` as a polynomial in `(u, v)` over F_p
/// (canonical integers), indexed `[a][b]` for `u^a v^b`, `a, b ∈ 0..=2`.
fn s3_coeffs(curve: &FastCurve, z_mont: u64) -> [[u64; 3]; 3] {
    let f = curve.f;
    let z = z_mont;
    let two = f.to_mont(2);
    let four = f.to_mont(4);
    let a = curve.a;
    let cb = curve.b;
    let z2 = f.sqr(z);
    let m2z = f.neg(f.mul(two, z));
    let c = |v: u64| f.from_mont(v);
    let mut out = [[0u64; 3]; 3];
    out[2][2] = c(f.one);
    out[2][1] = c(m2z);
    out[1][2] = c(m2z);
    out[2][0] = c(z2);
    out[0][2] = c(z2);
    out[1][1] = c(f.sub(f.neg(f.mul(two, z2)), f.mul(two, a)));
    let lin = f.sub(f.neg(f.mul(f.mul(two, a), z)), f.mul(four, cb));
    out[1][0] = c(lin);
    out[0][1] = c(lin);
    out[0][0] = c(f.add(f.neg(f.mul(f.mul(four, cb), z)), f.sqr(a)));
    out
}

type BiPoly = HashMap<(usize, usize), BigInt>;

fn poly_mul(p: &BiPoly, q: &BiPoly) -> BiPoly {
    let mut out: BiPoly = HashMap::new();
    for ((a1, b1), c1) in p {
        for ((a2, b2), c2) in q {
            *out.entry((a1 + a2, b1 + b2)).or_insert_with(BigInt::zero) += c1 * c2;
        }
    }
    out.retain(|_, c| !c.is_zero());
    out
}

fn poly_eval(p: &BiPoly, u: &BigInt, v: &BigInt) -> BigInt {
    let mut acc = BigInt::zero();
    for ((a, b), c) in p {
        acc += c * u.pow(*a as u32) * v.pow(*b as u32);
    }
    acc
}

/// One lattice trial: plant `x_i, x_j ≤ bound`, build the level-`m`
/// lattice, LLL, count reduced rows that vanish at the root over Z.
fn lattice_trial(
    curve: &FastCurve,
    m: usize,
    bound: u64,
    rng: &mut StdRng,
) -> Option<(usize, bool)> {
    let p = BigInt::from(curve.f.p);
    // plant
    let pick = |rng: &mut StdRng| -> Option<Pt> {
        for _ in 0..10_000 {
            let x = rng.gen_range(1..=bound);
            if let Some(pt) = curve.lift_x(x) {
                return Some(pt);
            }
        }
        None
    };
    let fi = pick(rng)?;
    let fj = pick(rng)?;
    let fj_signed = if rng.gen_bool(0.5) { fj } else { curve.neg(fj) };
    let r = curve.add(Some(fi), Some(fj_signed))?;
    let (xi, _) = curve.canonical(fi);
    let (xj, _) = curve.canonical(fj);
    let coeffs = s3_coeffs(curve, r.x);
    let mut f: BiPoly = HashMap::new();
    for a in 0..3 {
        for b in 0..3 {
            if coeffs[a][b] != 0 {
                f.insert((a, b), BigInt::from(coeffs[a][b]));
            }
        }
    }
    // sanity: f(xi, xj) ≡ 0 mod p
    let root_u = BigInt::from(xi);
    let root_v = BigInt::from(xj);
    debug_assert!((poly_eval(&f, &root_u, &root_v) % &p).is_zero());
    // powers of f
    let mut fpow: Vec<BiPoly> = vec![HashMap::from([((0usize, 0usize), BigInt::one())])];
    for k in 1..=m {
        let next = poly_mul(&fpow[k - 1], &f);
        fpow.push(next);
    }
    let x_bound = BigInt::from(bound);
    let dim = (2 * m + 1) * (2 * m + 1);
    let mono = |a: usize, b: usize| a * (2 * m + 1) + b;
    let mut basis: Vec<Vec<BigInt>> = Vec::with_capacity(dim);
    let mut rows_poly: Vec<BiPoly> = Vec::with_capacity(dim);
    for alpha in 0..=(2 * m) {
        for beta in 0..=(2 * m) {
            let k = (alpha / 2).min(beta / 2).min(m);
            let (a, b) = (alpha - 2 * k, beta - 2 * k);
            let pk = p.pow((m - k) as u32);
            let mut poly: BiPoly = HashMap::new();
            for ((pa, pb), c) in &fpow[k] {
                poly.insert((pa + a, pb + b), c * &pk);
            }
            let mut row = vec![BigInt::zero(); dim];
            for ((pa, pb), c) in &poly {
                row[mono(*pa, *pb)] = c * x_bound.pow(*pa as u32) * x_bound.pow(*pb as u32);
            }
            basis.push(row);
            rows_poly.push(poly);
        }
    }
    if lll_reduce(&mut basis, 0.99).is_err() {
        return None;
    }
    // unscale reduced rows into polynomials and test vanishing at the root over Z
    let mut vanishing = 0usize;
    let mut first_two_hg = true;
    let pm = p.pow(m as u32);
    let dimf = (dim as f64).sqrt();
    for (ri, row) in basis.iter().enumerate() {
        let mut poly: BiPoly = HashMap::new();
        let mut norm2 = BigInt::zero();
        for alpha in 0..=(2 * m) {
            for beta in 0..=(2 * m) {
                let c = &row[mono(alpha, beta)];
                if c.is_zero() {
                    continue;
                }
                norm2 += c * c;
                let scale = x_bound.pow(alpha as u32) * x_bound.pow(beta as u32);
                poly.insert((alpha, beta), c / &scale);
            }
        }
        if poly.is_empty() {
            continue;
        }
        let val = poly_eval(&poly, &root_u, &root_v);
        if val.is_zero() {
            vanishing += 1;
        }
        if ri < 2 {
            // Howgrave-Graham: ||g|| · √dim < p^m guarantees vanishing over Z.
            let norm = norm2.to_f64().unwrap_or(f64::INFINITY).sqrt();
            let ok = norm * dimf < pm.to_f64().unwrap_or(f64::INFINITY);
            first_two_hg &= ok;
        }
    }
    Some((vanishing, first_two_hg))
}

fn probe_b(bits_list: &[u32], out: &mut Vec<serde_json::Value>) {
    println!("## B. Lattice (Coppersmith) decomposition of S₃: empirical small-root reach\n");
    println!("Planted `R = F_i ± F_j` with `x_i, x_j ≤ p^δ`; level-`m` lattice of dimension `(2m+1)²`; success = at least two LLL-reduced polynomials vanish at the root over Z (the Coppersmith resultant condition). 8 trials per cell. Theory (basic lattice): δ < 1/18 (m=1), 1/10 (m=2), 0.119 (m=3), 1/6 (m→∞); IC-2 needs δ ≥ 1/4 to beat rho.\n");
    println!("| bits | m | dim | δ | B | successes/8 (≥2 vanishing) | HG condition met (first two rows) | mean vanishing rows | ms per LLL |");
    println!("|---|---|---:|---:|---:|---:|---:|---:|---:|");
    let deltas = [0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.20, 0.25];
    for &bits in bits_list {
        let Some((curve, _g)) = find_a3_curve(bits, 3) else {
            continue;
        };
        for m in 1..=3usize {
            for &delta in &deltas {
                let bound = ((curve.f.p as f64).powf(delta)).floor().max(2.0) as u64;
                let mut rng = StdRng::seed_from_u64(
                    0xB00B + bits as u64 * 100 + m as u64 * 10 + (delta * 100.0) as u64,
                );
                let mut succ = 0usize;
                let mut hg = 0usize;
                let mut vsum = 0usize;
                let mut done = 0usize;
                let t0 = Instant::now();
                for _ in 0..8 {
                    if let Some((v, h)) = lattice_trial(&curve, m, bound, &mut rng) {
                        done += 1;
                        vsum += v;
                        if v >= 2 {
                            succ += 1;
                        }
                        if h {
                            hg += 1;
                        }
                    }
                }
                let ms = t0.elapsed().as_secs_f64() * 1e3 / done.max(1) as f64;
                println!(
                    "| {bits} | {m} | {} | {delta:.2} | {bound} | {succ}/{done} | {hg}/{done} | {:.1} | {ms:.0} |",
                    (2 * m + 1) * (2 * m + 1),
                    vsum as f64 / done.max(1) as f64
                );
                out.push(json!({"probe":"B","bits":bits,"m":m,"dim":(2*m+1)*(2*m+1),"delta":delta,"bound":bound,
                    "trials":done,"successes":succ,"hg_first_two":hg,"mean_vanishing":vsum as f64 / done.max(1) as f64,"ms_per_lll":ms}));
                if succ == 0 && delta >= 0.16 {
                    break; // no point going higher once it has failed well past theory
                }
            }
        }
    }
    println!();
}

/// Which lattice to build for the `S₃` small-root problem.
#[derive(Clone, Copy, Debug)]
enum LatticeKind {
    /// Jochemsz–May basic/extended: monomials `u^α v^β`, `α ≤ 2m + t`, `β ≤ 2m`,
    /// row `u^{α−2k} v^{β−2k} f^k p^{m−k}` with `k = min(⌊α/2⌋, ⌊β/2⌋, m)`.
    Rect { m: usize, t: usize },
    /// Symmetric rewriting in `(e₁, e₂) = (x₁+x₂, x₁x₂)`, bounds `(2B, B²)`,
    /// triangle `α + β ≤ 2m`, row `e₁^α e₂^{β−2k} f^k p^{m−k}`, `k = min(⌊β/2⌋, m)`.
    Sym { m: usize },
}

/// One generalised lattice trial; returns (vanishing rows, HG-on-first-two, dim).
fn lattice_trial_kind(
    curve: &FastCurve,
    kind: LatticeKind,
    bound: u64,
    rng: &mut StdRng,
) -> Option<(usize, bool, usize)> {
    let p = BigInt::from(curve.f.p);
    let pick = |rng: &mut StdRng| -> Option<Pt> {
        for _ in 0..10_000 {
            let x = rng.gen_range(1..=bound);
            if let Some(pt) = curve.lift_x(x) {
                return Some(pt);
            }
        }
        None
    };
    let fi = pick(rng)?;
    let fj = pick(rng)?;
    let fj_signed = if rng.gen_bool(0.5) { fj } else { curve.neg(fj) };
    let r = curve.add(Some(fi), Some(fj_signed))?;
    let (xi, _) = curve.canonical(fi);
    let (xj, _) = curve.canonical(fj);
    let c = s3_coeffs(curve, r.x);
    let big = |v: u64| BigInt::from(v);
    // polynomial, bounds and monomial set per kind
    let (f, bounds, monos, m): (BiPoly, (BigInt, BigInt), Vec<(usize, usize)>, usize) = match kind {
        LatticeKind::Rect { m, t } => {
            let mut f: BiPoly = HashMap::new();
            for a in 0..3 {
                for b in 0..3 {
                    if c[a][b] != 0 {
                        f.insert((a, b), big(c[a][b]));
                    }
                }
            }
            let mut monos = Vec::new();
            for alpha in 0..=(2 * m + t) {
                for beta in 0..=(2 * m) {
                    monos.push((alpha, beta));
                }
            }
            (f, (big(bound), big(bound)), monos, m)
        }
        LatticeKind::Sym { m } => {
            // f(u,v) = Σ c_ab u^a v^b, symmetric: rewrite in e1 = u+v, e2 = uv.
            // S3 = (e1² − 4e2) z² − 2[e1(e2 + A) + 2B] z + (e2 − A)² − 4B e1.
            // Read the needed scalars back from the (u,v) coefficients:
            //   z² = c[2][0], −2z = c[2][1], −2z² − 2A = c[1][1], lin = c[1][0], const = c[0][0].
            let pm = &p;
            let md = |x: BigInt| ((x % pm) + pm) % pm;
            let z2 = big(c[2][0]);
            let m2z = big(c[2][1]); // −2z
            let c11 = big(c[1][1]); // −2z² − 2A
            let lin = big(c[1][0]); // −2Az − 4B
            let c00 = big(c[0][0]); // −4Bz + A²
                                    // e1²: z²; e2: −4z² + (coefficient of uv from (uv−A)² and −2(u+v)(uv)z ... ) — derive directly:
                                    // (e1² − 4e2) z²  → e1²·z², e2·(−4z²)
                                    // −2[e1 e2 + A e1 + 2B] z → e1e2·(−2z), e1·(−2Az), const·(−4Bz)
                                    // (e2 − A)² → e2²·1, e2·(−2A), const·A²
                                    // −4B e1 → e1·(−4B)
                                    // In terms of available scalars: e1 e2 coefficient = −2z = m2z; e1 coefficient = lin;
                                    // const = c00; e2 coefficient = −4z² − 2A = c11 − 2z² ; e2² = 1; e1² = z2.
            let mut f: BiPoly = HashMap::new();
            f.insert((2, 0), md(z2.clone()));
            f.insert((0, 2), BigInt::one());
            f.insert((1, 1), md(m2z));
            f.insert((1, 0), md(lin));
            f.insert((0, 1), md(c11 - BigInt::from(2u8) * &z2));
            f.insert((0, 0), md(c00));
            let mut monos = Vec::new();
            for alpha in 0..=(2 * m) {
                for beta in 0..=(2 * m) {
                    if alpha + beta <= 2 * m {
                        monos.push((alpha, beta));
                    }
                }
            }
            (f, (big(2 * bound), big(bound) * big(bound)), monos, m)
        }
    };
    // roots in the chosen variables
    let (ru, rv) = match kind {
        LatticeKind::Rect { .. } => (big(xi), big(xj)),
        LatticeKind::Sym { .. } => (big(xi) + big(xj), big(xi) * big(xj)),
    };
    debug_assert!((poly_eval(&f, &ru, &rv) % &p).is_zero());
    let mut fpow: Vec<BiPoly> = vec![HashMap::from([((0usize, 0usize), BigInt::one())])];
    for k in 1..=m {
        let next = poly_mul(&fpow[k - 1], &f);
        fpow.push(next);
    }
    let dim = monos.len();
    let index: HashMap<(usize, usize), usize> =
        monos.iter().enumerate().map(|(i, &mn)| (mn, i)).collect();
    let (xb, yb) = bounds;
    let mut basis: Vec<Vec<BigInt>> = Vec::with_capacity(dim);
    for &(alpha, beta) in &monos {
        let k = match kind {
            LatticeKind::Rect { m, .. } => (alpha / 2).min(beta / 2).min(m),
            LatticeKind::Sym { m } => (beta / 2).min(m),
        };
        let (a, b) = match kind {
            LatticeKind::Rect { .. } => (alpha - 2 * k, beta - 2 * k),
            LatticeKind::Sym { .. } => (alpha, beta - 2 * k),
        };
        let pk = p.pow((m - k) as u32);
        let mut row = vec![BigInt::zero(); dim];
        for ((pa, pb), cf) in &fpow[k] {
            let mono = (pa + a, pb + b);
            let &col = index.get(&mono)?; // shape mismatch guard
            row[col] += cf * &pk * xb.pow(mono.0 as u32) * yb.pow(mono.1 as u32);
        }
        basis.push(row);
    }
    if lll_reduce(&mut basis, 0.99).is_err() {
        return None;
    }
    let mut vanishing = 0usize;
    let mut hg = true;
    let pm = p.pow(m as u32);
    let dimf = (dim as f64).sqrt();
    for (ri, row) in basis.iter().enumerate() {
        let mut poly: BiPoly = HashMap::new();
        let mut norm2 = BigInt::zero();
        for (col, &(alpha, beta)) in monos.iter().enumerate() {
            let cf = &row[col];
            if cf.is_zero() {
                continue;
            }
            norm2 += cf * cf;
            poly.insert(
                (alpha, beta),
                cf / (xb.pow(alpha as u32) * yb.pow(beta as u32)),
            );
        }
        if poly.is_empty() {
            continue;
        }
        if poly_eval(&poly, &ru, &rv).is_zero() {
            vanishing += 1;
        }
        if ri < 2 {
            let norm = norm2.to_f64().unwrap_or(f64::INFINITY).sqrt();
            hg &= norm * dimf < pm.to_f64().unwrap_or(f64::INFINITY);
        }
    }
    Some((vanishing, hg, dim))
}

fn probe_c(bits: u32, out: &mut Vec<serde_json::Value>) {
    println!("## C. Best lattice reach for S₃: extended Jochemsz–May shifts and the symmetric (e₁, e₂) lattice\n");
    println!("Same planted instances and success criteria as probe B. `Rect{{m,t}}` = basic lattice at level m with t extra u-shifts; `Sym{{m}}` = symmetric formulation (bounds 2B, B²). 8 trials per cell on the {bits}-bit curve.\n");
    println!("| lattice | dim | δ | B | ≥2 vanishing / trials | HG first two | mean vanishing | ms per LLL |");
    println!("|---|---:|---:|---:|---:|---:|---:|---:|");
    let Some((curve, _g)) = find_a3_curve(bits, 3) else {
        return;
    };
    let kinds = [
        LatticeKind::Rect { m: 2, t: 1 },
        LatticeKind::Rect { m: 2, t: 2 },
        LatticeKind::Rect { m: 3, t: 1 },
        LatticeKind::Sym { m: 1 },
        LatticeKind::Sym { m: 2 },
        LatticeKind::Sym { m: 3 },
    ];
    let deltas = [0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20];
    for kind in kinds {
        for &delta in &deltas {
            let bound = ((curve.f.p as f64).powf(delta)).floor().max(2.0) as u64;
            let mut rng =
                StdRng::seed_from_u64(0xC0DE + bits as u64 * 100 + (delta * 100.0) as u64);
            let (mut succ, mut hgc, mut vsum, mut done, mut dim) =
                (0usize, 0usize, 0usize, 0usize, 0usize);
            let t0 = Instant::now();
            for _ in 0..8 {
                if let Some((v, h, d)) = lattice_trial_kind(&curve, kind, bound, &mut rng) {
                    done += 1;
                    vsum += v;
                    dim = d;
                    if v >= 2 {
                        succ += 1;
                    }
                    if h {
                        hgc += 1;
                    }
                }
            }
            let ms = t0.elapsed().as_secs_f64() * 1e3 / done.max(1) as f64;
            println!("| {kind:?} | {dim} | {delta:.2} | {bound} | {succ}/{done} | {hgc}/{done} | {:.1} | {ms:.0} |", vsum as f64 / done.max(1) as f64);
            out.push(json!({"probe":"C","bits":bits,"lattice":format!("{kind:?}"),"dim":dim,"delta":delta,"bound":bound,"trials":done,"successes":succ,"hg_first_two":hgc,"mean_vanishing":vsum as f64/done.max(1) as f64,"ms_per_lll":ms}));
            if succ == 0 && delta >= 0.14 {
                break;
            }
        }
    }
    println!();
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let mut probe = "ab".to_string();
    let mut rho_bits = vec![32u32, 36, 40, 44];
    let mut lattice_bits = vec![40u32, 56];
    let mut targets = 16usize;
    let mut json_path: Option<String> = None;
    let mut i = 1;
    while i < argv.len() {
        match argv[i].as_str() {
            "--probe" => {
                i += 1;
                probe = argv.get(i).cloned().unwrap_or_else(|| "ab".into());
            }
            "--rho-bits" => {
                i += 1;
                rho_bits = argv
                    .get(i)
                    .map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect())
                    .unwrap_or(rho_bits);
            }
            "--lattice-bits" => {
                i += 1;
                lattice_bits = argv
                    .get(i)
                    .map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect())
                    .unwrap_or(lattice_bits);
            }
            "--targets" => {
                i += 1;
                targets = argv.get(i).and_then(|s| s.parse().ok()).unwrap_or(16);
            }
            "--json" => {
                i += 1;
                json_path = argv.get(i).cloned();
            }
            other => eprintln!("ignoring {other}"),
        }
        i += 1;
    }
    let mut out = Vec::new();
    println!("# Prime-field ECDLP breakthrough probe\n");
    println!(
        "Host threads: {}.\n",
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    );
    if probe.contains('a') {
        probe_a(&rho_bits, targets, &mut out);
    }
    if probe.contains('b') {
        probe_b(&lattice_bits, &mut out);
    }
    if probe.contains('c') {
        probe_c(lattice_bits[0], &mut out);
    }
    if let Some(path) = json_path {
        std::fs::write(
            &path,
            serde_json::to_string_pretty(&json!({"results": out})).unwrap(),
        )
        .unwrap();
        println!("JSON written to {path}");
    }
    let _ = BigUint::zero();
}
