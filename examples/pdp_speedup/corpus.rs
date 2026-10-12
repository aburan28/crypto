//! The frozen verdict corpus of PS-0: targets, an independent exhaustive
//! root enumerator, digests and the on-disk format.
//!
//! The enumerator here deliberately shares nothing with the baseline's
//! `decompose` beyond the coefficients of the quartic (`quartic_in_x3`) and
//! the field: it normalises each quartic, forms `t^(2^i) mod f` by generic
//! polynomial squaring, takes a generic Euclidean gcd with `L_V`, and falls
//! back to evaluating over `V`.  The baseline's Horner form, Frobenius table
//! and pseudo-remainder gcd are all absent.

use crate::common::*;
use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::binary_semaev_s4::{elementary_symmetric_3, symmetrised_s4_eval};
use crypto_lib::cryptanalysis::semaev_decomp::{quartic_in_x3, subspace_poly, Gf2};
use crypto_lib::hash::sha256::sha256;
use std::sync::atomic::{AtomicU64, Ordering};

pub const SEED_BASE: u64 = 20261011;
pub const N_UNIFORM: usize = 64;
pub const N_PLANTED: usize = 16;

#[derive(Clone, Debug)]
pub struct Target {
    pub planted: bool,
    pub idx: u32,
    pub xr: u64,
    pub yr: u64,
    /// The planted triple (sorted), if any.
    pub triple: Option<[u64; 3]>,
}

#[derive(Clone, Debug)]
pub struct Record {
    pub t: Target,
    /// Every sorted triple `(x₁ ≤ x₂ ≤ x₃)` in `V³` with `S₄(x₁,x₂,x₃,x_R) = 0`,
    /// sorted lexicographically.
    pub roots: Vec<[u64; 3]>,
    /// For each root: whether it re-adds to `±R` in the group.
    pub rational: Vec<bool>,
}

pub fn make_targets(rung: &Rung, n_uniform: usize, n_planted: usize) -> Vec<Target> {
    let mut rng = SplitMix64(SEED_BASE + rung.l as u64);
    let mut out = Vec::new();
    // Uniform: a uniformly random abscissa among those that lift to a rational point.
    while out.len() < n_uniform {
        let x = rng.next() & rung.mask();
        if x == 0 {
            continue;
        }
        let ys = rung.lifts(x);
        if ys.is_empty() {
            continue;
        }
        let y = ys[(rng.next() as usize) % ys.len()];
        out.push(Target {
            planted: false,
            idx: out.len() as u32,
            xr: x,
            yr: y,
            triple: None,
        });
    }
    // Planted: R = P₁ + P₂ + P₃ with every abscissa in V.
    let mut planted = 0u32;
    while (planted as usize) < n_planted {
        let mut pts = Vec::new();
        let mut xs = [0u64; 3];
        let mut ok = true;
        for slot in xs.iter_mut() {
            let x = rng.next() & (rung.span() - 1);
            let ys = if x == 0 { vec![] } else { rung.lifts(x) };
            if ys.is_empty() {
                ok = false;
                break;
            }
            let y = ys[(rng.next() as usize) % ys.len()];
            *slot = x;
            pts.push(rung.point(x, y));
        }
        if !ok {
            continue;
        }
        let r = rung.add(&rung.add(&pts[0], &pts[1]), &pts[2]);
        let Some((xr, yr)) = rung.xy(&r) else {
            continue;
        };
        if xr == 0 {
            continue;
        }
        out.push(Target {
            planted: true,
            idx: planted,
            xr,
            yr,
            triple: Some(sorted3(xs[0], xs[1], xs[2])),
        });
        planted += 1;
    }
    out
}

// ── generic polynomial arithmetic over the one-word field ───────────────

type P5 = [u64; 5];

fn deg5(p: &P5) -> Option<usize> {
    (0..5).rev().find(|&i| p[i] != 0)
}

/// Remainder of `a` by `b` (`b ≠ 0`), with one inversion of `lead(b)`.
fn rem5(mut a: P5, b: &P5, gf: &Gf2) -> P5 {
    let db = deg5(b).expect("divisor is non-zero");
    let inv = gf.inv(b[db]);
    while let Some(da) = deg5(&a) {
        if da < db {
            break;
        }
        let f = gf.mul(a[da], inv);
        for i in 0..=db {
            a[da - db + i] ^= gf.mul(b[i], f);
        }
        a[da] = 0;
    }
    a
}

fn gcd5(mut a: P5, mut b: P5, gf: &Gf2) -> P5 {
    while deg5(&b).is_some() {
        let r = rem5(a, &b, gf);
        a = b;
        b = r;
    }
    a
}

/// Square a residue modulo a monic `f` of degree `d` (`1 ≤ d ≤ 4`).
fn sqr_mod(u: &[u64; 4], f: &P5, d: usize, gf: &Gf2) -> [u64; 4] {
    let mut wide = [0u64; 8];
    for i in 0..d {
        wide[2 * i] = gf.sqr(u[i]);
    }
    for k in (d..2 * d).rev() {
        let c = wide[k];
        if c != 0 {
            for j in 0..d {
                wide[k - d + j] ^= gf.mul(f[j], c);
            }
            wide[k] = 0;
        }
    }
    [wide[0], wide[1], wide[2], wide[3]]
}

/// All sorted triples `(x₁ ≤ x₂ ≤ x₃)` in `V³` with `S₄ = 0` that arise from
/// the pair `(x₁, x₂)`, `x₁ ≤ x₂`; `x₃` ranges over all of `V`, so triples
/// are re-sorted on output and de-duplicated by the caller.
pub fn pair_roots(rung: &Rung, lv: &[u64], xr: u64, x1: u64, x2: u64, out: &mut Vec<[u64; 3]>) {
    let gf = &rung.gf;
    let q = quartic_in_x3(x1, x2, xr, gf);
    let Some(d) = q.deg() else {
        for t in 0..rung.span() {
            out.push(sorted3(x1, x2, t));
        }
        return;
    };
    if d == 0 {
        return;
    }
    let inv = gf.inv(q.c[d]);
    let mut f: P5 = [0; 5];
    for i in 0..=d.min(4) {
        f[i] = gf.mul(q.c[i], inv);
    }
    debug_assert!(d <= 4);
    // pow = t mod f.
    let mut pow = [0u64; 4];
    if d == 1 {
        pow[0] = f[0];
    } else {
        pow[1] = 1;
    }
    let mut acc: P5 = [0; 5];
    for (i, &ai) in lv.iter().enumerate() {
        if ai != 0 {
            for j in 0..d {
                acc[j] ^= gf.mul(pow[j], ai);
            }
        }
        if i + 1 < lv.len() {
            pow = sqr_mod(&pow, &f, d, gf);
        }
    }
    let g = if deg5(&acc).is_none() {
        f
    } else {
        gcd5(f, acc, gf)
    };
    match deg5(&g) {
        None | Some(0) => {}
        Some(1) => {
            let t = gf.mul(g[0], gf.inv(g[1]));
            debug_assert!(t < rung.span());
            out.push(sorted3(x1, x2, t));
        }
        Some(_) => {
            // Several subspace roots: evaluate f over V.
            for t in 0..rung.span() {
                let mut v = 0u64;
                for j in (0..=d).rev() {
                    v = gf.mul(v, t) ^ f[j];
                }
                if v == 0 {
                    out.push(sorted3(x1, x2, t));
                }
            }
        }
    }
}

/// Independent exhaustive enumeration of the root set of one target.
pub fn all_roots(rung: &Rung, xr: u64, threads: usize) -> Vec<[u64; 3]> {
    let lv = subspace_poly(rung.l, &rung.gf);
    let next = AtomicU64::new(0);
    let span = rung.span();
    let mut all: Vec<[u64; 3]> = Vec::new();
    std::thread::scope(|s| {
        let hs: Vec<_> = (0..threads)
            .map(|_| {
                s.spawn(|| {
                    let mut local = Vec::new();
                    loop {
                        let x1 = next.fetch_add(1, Ordering::Relaxed);
                        if x1 >= span {
                            break;
                        }
                        for x2 in x1..span {
                            pair_roots(rung, &lv, xr, x1, x2, &mut local);
                        }
                    }
                    local
                })
            })
            .collect();
        for h in hs {
            all.extend(h.join().unwrap());
        }
    });
    all.sort_unstable();
    all.dedup();
    all
}

/// Brute force over every sorted triple with the crate's general
/// (`Vec`-backed) `S₄` evaluator: the cross-check for small `ℓ`.
pub fn brute_roots(rung: &Rung, xr: u64, threads: usize) -> Vec<[u64; 3]> {
    let span = rung.span();
    let irr = &rung.irr;
    let base: Vec<F2mElement> = (0..span).map(|v| rung.gf.to_element(v)).collect();
    let xre = rung.gf.to_element(xr);
    let next = AtomicU64::new(0);
    let mut all = Vec::new();
    std::thread::scope(|s| {
        let hs: Vec<_> = (0..threads)
            .map(|_| {
                s.spawn(|| {
                    let mut local = Vec::new();
                    loop {
                        let a = next.fetch_add(1, Ordering::Relaxed);
                        if a >= span {
                            break;
                        }
                        for b in a..span {
                            for c in b..span {
                                let (e1, e2, e3) = elementary_symmetric_3(
                                    &base[a as usize],
                                    &base[b as usize],
                                    &base[c as usize],
                                    irr,
                                );
                                if symmetrised_s4_eval(&e1, &e2, &e3, &xre, irr).is_zero() {
                                    local.push([a, b, c]);
                                }
                            }
                        }
                    }
                    local
                })
            })
            .collect();
        for h in hs {
            all.extend(h.join().unwrap());
        }
    });
    all.sort_unstable();
    all
}

pub fn build_records(
    rung: &Rung,
    threads: usize,
    n_uniform: usize,
    n_planted: usize,
) -> Vec<Record> {
    make_targets(rung, n_uniform, n_planted)
        .into_iter()
        .map(|t| {
            let roots = all_roots(rung, t.xr, threads);
            let r = rung.point(t.xr, t.yr);
            assert!(rung.curve.is_on_curve(&r), "target is not on the curve");
            let rational: Vec<bool> = roots
                .iter()
                .map(|&tr| match rung.readd(tr, &r) {
                    Readd::Sums => true,
                    Readd::NotRational => false,
                    Readd::RationalButNoSum => {
                        panic!("root {tr:?} lifts but does not sum to ±R: oracle or lift bug")
                    }
                })
                .collect();
            if let Some(p) = t.triple {
                assert!(
                    roots.binary_search(&p).is_ok(),
                    "planted triple {p:?} missing from the enumerated root set"
                );
            }
            Record { t, roots, rational }
        })
        .collect()
}

// ── digests ─────────────────────────────────────────────────────────────

/// SHA-256 of the verdict stream: per target `[planted u8][x_R u64 LE][found u8]`.
pub fn digest_verdicts(recs: &[Record]) -> String {
    let mut buf = Vec::new();
    for r in recs {
        buf.push(r.t.planted as u8);
        buf.extend_from_slice(&r.t.xr.to_le_bytes());
        buf.push(!r.roots.is_empty() as u8);
    }
    hex(&sha256(&buf))
}

/// SHA-256 of the root-set stream: the verdict record followed by
/// `[count u32 LE]` and the sorted triples (three `u64 LE` each).
pub fn digest_roots(recs: &[Record]) -> String {
    let mut buf = Vec::new();
    for r in recs {
        buf.push(r.t.planted as u8);
        buf.extend_from_slice(&r.t.xr.to_le_bytes());
        buf.push(!r.roots.is_empty() as u8);
        buf.extend_from_slice(&(r.roots.len() as u32).to_le_bytes());
        for t in &r.roots {
            for w in t {
                buf.extend_from_slice(&w.to_le_bytes());
            }
        }
    }
    hex(&sha256(&buf))
}

// ── on-disk format ──────────────────────────────────────────────────────

pub fn write_corpus(rung: &Rung, recs: &[Record]) -> String {
    let mut s = String::new();
    s.push_str("# pdp_speedup corpus v1 (research/pdp_oracle_speedup_20261011)\n");
    s.push_str(&format!("l {}\nn {}\n", rung.l, rung.n));
    s.push_str(&format!("modulus 0x{:x}\n", rung.modulus()));
    s.push_str(&format!("seed {}\n", SEED_BASE + rung.l as u64));
    s.push_str("curve y^2+xy=x^3+x^2+1\n");
    s.push_str(&format!("targets {}\n", recs.len()));
    for r in recs {
        let planted = match r.t.triple {
            Some(t) => format!("{:x},{:x},{:x}", t[0], t[1], t[2]),
            None => "-".into(),
        };
        s.push_str(&format!(
            "T {} {} {:x} {:x} {} {} {} {}\n",
            r.t.idx,
            if r.t.planted { "p" } else { "u" },
            r.t.xr,
            r.t.yr,
            !r.roots.is_empty() as u8,
            r.roots.len(),
            r.rational.iter().filter(|&&b| b).count(),
            planted
        ));
        for (t, &rat) in r.roots.iter().zip(&r.rational) {
            s.push_str(&format!(
                "R {:x} {:x} {:x} {}\n",
                t[0], t[1], t[2], rat as u8
            ));
        }
    }
    s.push_str(&format!("digest_verdicts {}\n", digest_verdicts(recs)));
    s.push_str(&format!("digest_roots {}\n", digest_roots(recs)));
    s
}

pub struct Loaded {
    pub l: u32,
    pub modulus: u64,
    pub recs: Vec<Record>,
    pub file_digest_verdicts: String,
    pub file_digest_roots: String,
}

pub fn read_corpus(text: &str) -> Loaded {
    let mut l = 0u32;
    let mut modulus = 0u64;
    let mut recs: Vec<Record> = Vec::new();
    let mut dv = String::new();
    let mut dr = String::new();
    let h = |s: &str| u64::from_str_radix(s, 16).unwrap();
    for line in text.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        match f.first().copied() {
            Some("l") => l = f[1].parse().unwrap(),
            Some("modulus") => modulus = h(f[1].trim_start_matches("0x")),
            Some("digest_verdicts") => dv = f[1].to_string(),
            Some("digest_roots") => dr = f[1].to_string(),
            Some("T") => {
                let triple = if f[8] == "-" {
                    None
                } else {
                    let p: Vec<u64> = f[8].split(',').map(h).collect();
                    Some([p[0], p[1], p[2]])
                };
                recs.push(Record {
                    t: Target {
                        planted: f[2] == "p",
                        idx: f[1].parse().unwrap(),
                        xr: h(f[3]),
                        yr: h(f[4]),
                        triple,
                    },
                    roots: vec![],
                    rational: vec![],
                });
            }
            Some("R") => {
                let r = recs.last_mut().unwrap();
                r.roots.push([h(f[1]), h(f[2]), h(f[3])]);
                r.rational.push(f[4] == "1");
            }
            _ => {}
        }
    }
    Loaded {
        l,
        modulus,
        recs,
        file_digest_verdicts: dv,
        file_digest_roots: dr,
    }
}
