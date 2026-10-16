//! Invariants finer than the `L`-polynomial for Smith's `full_success`
//! (`research/smith_class_invariant_20261015/PROTOCOL.md`).
//!
//! ```text
//! cargo run --release --example smith_class_invariant -- compute --p 7 --out DIR
//! cargo run --release --example smith_class_invariant -- summarize --dir DIR 7 11
//! ```
//!
//! `compute` reads the frozen per-isomorphism-class rows
//! `research/smith_iso_granularity_20261015/data/iso_p<p>.csv.gz` and, for
//! every tractable class, computes on the representative model `y^2 = t f(x)`:
//!
//! - X1: the Frobenius cycle type on the eight Weierstrass points;
//! - X2: the 2-Sylow subgroup of `J(F_p)`, from the class group of the
//!   function field (places of degree <= 3, relations `div(a)` and
//!   `div(y - v)`, Smith form over `Z/2^E`, certified by its order);
//! - X3 / X3+: halvability of the `F_{p^k}`-rational 2-torsion points by
//!   2-descent over `F_{p^k}` for every `k | 2m` (`m` the order of Frobenius
//!   on the roots), giving `#J(F_{p^k})[4]` and the halvable sets in a
//!   canonical form under relabelings commuting with Frobenius.
//!
//! It writes `inv_p<p>.csv.gz` (one row per isomorphism class, in the frozen
//! order) and `controls_p<p>.json` (K2-K5). `summarize` runs the conditioned
//! and placebo permutation nulls, control K1, and writes `summary_p<p>.json`
//! and `verdict.json`.

use crypto_lib::cryptanalysis::smith_prop24 as sp;
use flate2::{write::GzEncoder, Compression};
use rayon::prelude::*;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::{fs, time::Instant};

const FROZEN: &str = "research/smith_iso_granularity_20261015/data";
const PERMS: usize = 999;
const SEED: u64 = 20261015;
const PLACEBO_SEED: u64 = 20261016;
const K2_SAMPLE: usize = 1000;
const MAX_ROUNDS: usize = 40;

fn arg(args: &[String], k: &str) -> Option<String> {
    args.iter()
        .position(|a| a == k)
        .and_then(|i| args.get(i + 1).cloned())
}

// ---------------------------------------------------------------------------
// Polynomials over F_p (low degree first, trimmed).

type P = Vec<u64>;

fn trim(a: &mut P) {
    while a.last() == Some(&0) {
        a.pop();
    }
}

fn pmul(a: &[u64], b: &[u64], p: u64) -> P {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = (r[i + j] + x * y) % p;
        }
    }
    trim(&mut r);
    r
}

fn psub(a: &[u64], b: &[u64], p: u64) -> P {
    let n = a.len().max(b.len());
    let mut r: P = (0..n)
        .map(|i| {
            let x = a.get(i).copied().unwrap_or(0);
            let y = b.get(i).copied().unwrap_or(0);
            (x + p - y) % p
        })
        .collect();
    trim(&mut r);
    r
}

fn prem(a: &[u64], m: &[u64], p: u64) -> P {
    let mut r = a.to_vec();
    trim(&mut r);
    let dm = m.len() - 1;
    let li = sp::inv_mod(m[dm], p);
    while r.len() > dm {
        let s = r.len() - 1 - dm;
        let c = r[r.len() - 1] * li % p;
        for (i, &mi) in m.iter().enumerate() {
            r[s + i] = (r[s + i] + p * p - c * mi % p) % p;
        }
        trim(&mut r);
    }
    r
}

fn ppowmod(a: &[u64], mut e: u128, m: &[u64], p: u64) -> P {
    let mut r = vec![1u64];
    let mut b = prem(a, m, p);
    while e > 0 {
        if e & 1 == 1 {
            r = prem(&pmul(&r, &b, p), m, p);
        }
        b = prem(&pmul(&b, &b, p), m, p);
        e >>= 1;
    }
    r
}

// ---------------------------------------------------------------------------
// F_{p^m}, m <= 15, as F_p[z]/(h).

type E = [u64; 16];

#[derive(Clone)]
struct Fq {
    p: u64,
    m: usize,
    h: Vec<u64>,
}

impl Fq {
    fn new(p: u64, m: usize, rng: &mut sp::SplitMix64) -> Fq {
        assert!((1..=15).contains(&m));
        if m == 1 {
            return Fq {
                p,
                m,
                h: vec![0, 1],
            };
        }
        loop {
            let mut h: P = (0..m).map(|_| rng.below(p)).collect();
            h.push(1);
            if let Some(fs) = sp::factor_squarefree(&h, p, rng) {
                if fs.len() == 1 {
                    return Fq { p, m, h };
                }
            }
        }
    }
    fn q(&self) -> u128 {
        (self.p as u128).pow(self.m as u32)
    }
    fn c(&self, x: u64) -> E {
        let mut a = [0; 16];
        a[0] = x % self.p;
        a
    }
    fn one(&self) -> E {
        self.c(1)
    }
    fn is_zero(&self, a: &E) -> bool {
        a[..self.m].iter().all(|&x| x == 0)
    }
    fn add(&self, a: &E, b: &E) -> E {
        let mut r = [0; 16];
        for i in 0..self.m {
            r[i] = (a[i] + b[i]) % self.p;
        }
        r
    }
    fn sub(&self, a: &E, b: &E) -> E {
        let mut r = [0; 16];
        for i in 0..self.m {
            r[i] = (a[i] + self.p - b[i]) % self.p;
        }
        r
    }
    fn neg(&self, a: &E) -> E {
        self.sub(&[0; 16], a)
    }
    fn mul(&self, a: &E, b: &E) -> E {
        let (m, p) = (self.m, self.p);
        let mut r = [0u64; 31];
        for i in 0..m {
            if a[i] == 0 {
                continue;
            }
            for j in 0..m {
                r[i + j] += a[i] * b[j];
            }
        }
        for k in (m..2 * m - 1).rev() {
            let c = r[k] % p;
            if c == 0 {
                continue;
            }
            for i in 0..m {
                r[k - m + i] += (p - self.h[i]) * c;
            }
        }
        let mut out = [0u64; 16];
        for i in 0..m {
            out[i] = r[i] % p;
        }
        out
    }
    fn pow(&self, a: &E, mut e: u128) -> E {
        let mut r = self.one();
        let mut b = *a;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(&r, &b);
            }
            b = self.mul(&b, &b);
            e >>= 1;
        }
        r
    }
    fn inv(&self, a: &E) -> E {
        assert!(!self.is_zero(a), "inverse of zero");
        self.pow(a, self.q() - 2)
    }
    fn rand(&self, rng: &mut sp::SplitMix64) -> E {
        let mut a = [0; 16];
        for x in a.iter_mut().take(self.m) {
            *x = rng.below(self.p);
        }
        a
    }
    /// `x^{p^j}`.
    fn frob(&self, a: &E, j: u32) -> E {
        self.pow(a, (self.p as u128).pow(j))
    }
    /// Quadratic character over the subfield `F_{p^g}` (`a` must lie in it).
    fn chi_sub(&self, a: &E, g: u32) -> Option<i8> {
        let qg = (self.p as u128).pow(g);
        if self.pow(a, qg) != *a {
            return None;
        }
        if self.is_zero(a) {
            return Some(0);
        }
        let t = self.pow(a, (qg - 1) / 2);
        if t == self.one() {
            Some(1)
        } else if t == self.neg(&self.one()) {
            Some(-1)
        } else {
            None
        }
    }
    /// Trace from `F_{p^m}` down to `F_{p^g}`.
    fn trace_to(&self, a: &E, g: u32) -> E {
        let mut s = [0; 16];
        let mut x = *a;
        for _ in 0..(self.m as u32 / g) {
            s = self.add(&s, &x);
            x = self.frob(&x, g);
        }
        s
    }
}

// Polynomials over F_q (low degree first).

fn qtrim(k: &Fq, a: &mut Vec<E>) {
    while a.last().is_some_and(|x| k.is_zero(x)) {
        a.pop();
    }
}

fn qdivrem(k: &Fq, a: &[E], b: &[E]) -> (Vec<E>, Vec<E>) {
    let mut r = a.to_vec();
    qtrim(k, &mut r);
    let mut b = b.to_vec();
    qtrim(k, &mut b);
    let db = b.len() - 1;
    if r.len() <= db {
        return (vec![], r);
    }
    let li = k.inv(&b[db]);
    let mut q = vec![[0; 16]; r.len() - db];
    while r.len() > db {
        let s = r.len() - 1 - db;
        let c = k.mul(&r[r.len() - 1], &li);
        q[s] = c;
        for (i, bi) in b.iter().enumerate() {
            r[s + i] = k.sub(&r[s + i], &k.mul(&c, bi));
        }
        qtrim(k, &mut r);
    }
    (q, r)
}

fn qmulmod(k: &Fq, a: &[E], b: &[E], m: &[E]) -> Vec<E> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![[0; 16]; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            r[i + j] = k.add(&r[i + j], &k.mul(x, y));
        }
    }
    qdivrem(k, &r, m).1
}

fn qmonic(k: &Fq, a: &[E]) -> Vec<E> {
    let mut a = a.to_vec();
    qtrim(k, &mut a);
    if let Some(l) = a.last().copied() {
        let li = k.inv(&l);
        for x in a.iter_mut() {
            *x = k.mul(x, &li);
        }
    }
    a
}

fn qgcd(k: &Fq, a: &[E], b: &[E]) -> Vec<E> {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    qtrim(k, &mut a);
    qtrim(k, &mut b);
    while !b.is_empty() {
        let r = qdivrem(k, &a, &b).1;
        a = b;
        b = r;
    }
    qmonic(k, &a)
}

/// One root in `F_q` of a squarefree polynomial over `F_q` that splits in `F_q`.
fn one_root(k: &Fq, f: &[E], rng: &mut sp::SplitMix64) -> E {
    let mut f = qmonic(k, f);
    let e = (k.q() - 1) / 2;
    loop {
        let n = f.len() - 1;
        if n == 1 {
            return k.neg(&f[0]);
        }
        let mut a: Vec<E> = (0..n).map(|_| k.rand(rng)).collect();
        qtrim(k, &mut a);
        if a.len() < 2 {
            continue;
        }
        let mut r = vec![k.one()];
        let mut b = qdivrem(k, &a, &f).1;
        let mut ee = e;
        while ee > 0 {
            if ee & 1 == 1 {
                r = qmulmod(k, &r, &b, &f);
            }
            b = qmulmod(k, &b, &b, &f);
            ee >>= 1;
        }
        if r.is_empty() {
            r.push([0; 16]);
        }
        r[0] = k.sub(&r[0], &k.one());
        let h = qgcd(k, &r, &f);
        let dh = h.len() - 1;
        if dh > 0 && dh < n {
            let q = qmonic(k, &qdivrem(k, &f, &h).0);
            f = if dh < q.len() { h } else { q };
        }
    }
}

/// Square root in the subfield `F_{p^g}` (Tonelli-Shanks); `a` a nonzero square there.
fn sqrt_sub(k: &Fq, a: &E, g: u32, rng: &mut sp::SplitMix64) -> E {
    let qg = (k.p as u128).pow(g);
    let mut z;
    loop {
        z = k.trace_to(&k.rand(rng), g);
        if k.chi_sub(&z, g) == Some(-1) {
            break;
        }
    }
    let (mut qq, mut s) = (qg - 1, 0u32);
    while qq % 2 == 0 {
        qq /= 2;
        s += 1;
    }
    let mut mm = s;
    let mut c = k.pow(&z, qq);
    let mut t = k.pow(a, qq);
    let mut r = k.pow(a, qq.div_ceil(2));
    let one = k.one();
    while t != one {
        let mut i = 0;
        let mut tt = t;
        while tt != one {
            tt = k.mul(&tt, &tt);
            i += 1;
        }
        let mut b = c;
        for _ in 0..(mm - i - 1) {
            b = k.mul(&b, &b);
        }
        mm = i;
        c = k.mul(&b, &b);
        t = k.mul(&t, &c);
        r = k.mul(&r, &b);
    }
    r
}

// ---------------------------------------------------------------------------
// Roots of f in its splitting field, with the Frobenius labelling.

struct Split {
    k: Fq,
    /// Roots `e_j`; factor `t` of degree `d` owns a block `r, r^p, ...`.
    roots: Vec<E>,
    /// Frobenius on root indices.
    sigma: [usize; 8],
    /// Root index -> (cycle id, position).
    cycles: Vec<Vec<usize>>,
    ct: Vec<u32>,
}

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn lcm(a: u64, b: u64) -> u64 {
    a / gcd(a, b) * b
}

fn split(f: &[u64], p: u64, rng: &mut sp::SplitMix64) -> Option<Split> {
    let fs = sp::factor_squarefree(f, p, rng)?;
    let m = fs.iter().fold(1u64, |a, g| lcm(a, (g.len() - 1) as u64)) as usize;
    let k = Fq::new(p, m, rng);
    let mut roots = vec![];
    let mut sigma = [0usize; 8];
    let mut cycles = vec![];
    for g in &fs {
        let ge: Vec<E> = g.iter().map(|&c| k.c(c)).collect();
        let mut r = one_root(&k, &ge, rng);
        let d = g.len() - 1;
        let base = roots.len();
        let mut cyc = vec![];
        for i in 0..d {
            roots.push(r);
            sigma[base + i] = base + (i + 1) % d;
            cyc.push(base + i);
            r = k.frob(&r, 1);
        }
        cycles.push(cyc);
    }
    assert_eq!(roots.len(), 8);
    let mut ct: Vec<u32> = fs.iter().map(|g| (g.len() - 1) as u32).collect();
    ct.sort_unstable();
    Some(Split {
        k,
        roots,
        sigma,
        cycles,
        ct,
    })
}

fn ct_string(ct: &[u32]) -> String {
    ct.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(".")
}

fn perm_pow(s: &[usize; 8], k: u64) -> [usize; 8] {
    let mut r = [0usize; 8];
    for (j, rj) in r.iter_mut().enumerate() {
        let mut x = j;
        for _ in 0..k {
            x = s[x];
        }
        *rj = x;
    }
    r
}

fn apply_mask(perm: &[usize; 8], mask: u32) -> u32 {
    (0..8)
        .filter(|&j| mask >> j & 1 == 1)
        .fold(0, |a, j| a | 1 << perm[j])
}

// ---------------------------------------------------------------------------
// 2-descent over K = F_{p^k}.

#[derive(Default, Clone)]
struct DescentK {
    k: u64,
    /// `#J(K)[2]`.
    n2: u32,
    /// Rational 2-torsion points with trivial descent image (0 included).
    halvable: u32,
    /// Halvable set, bit `mask >> 1` for the representative mask without root 0.
    set: u128,
    /// Field-membership or character failures (must stay 0).
    errors: u32,
    /// Every Frobenius orbit on the roots over K has even length.
    all_even: bool,
    /// Corrected halvable count (Amendment 1): `h_lo == h_hi` when exact.
    h_lo: u32,
    h_hi: u32,
    /// Descent-trivial stable points with `B(T, y_0) = 1` (all of `set` if not all-even).
    stable_true: u128,
}

/// Element `β ∈ F_{q^2} \ F_q` with `β^{q+1} = lc`, `q = p^k`, `2k | m`.
fn find_beta(k: &Fq, kk: u32, lc: u64, rng: &mut sp::SplitMix64) -> E {
    let q = (k.p as u128).pow(kk);
    let lce = k.c(lc);
    loop {
        let w = k.trace_to(&k.rand(rng), 2 * kk);
        if k.is_zero(&w) || k.pow(&w, q) == w {
            continue;
        }
        let n = k.pow(&w, q + 1);
        let ratio = k.mul(&lce, &k.inv(&n));
        if k.chi_sub(&ratio, kk) != Some(1) {
            continue;
        }
        let mu = sqrt_sub(k, &ratio, kk, rng);
        let beta = k.mul(&w, &mu);
        assert!(k.pow(&beta, q + 1) == lce && k.pow(&beta, q) != beta);
        return beta;
    }
}

fn descent(sp_: &Split, lc: u64, kk: u64, rng: &mut sp::SplitMix64) -> DescentK {
    let k = &sp_.k;
    let e = &sp_.roots;
    let m = k.m as u64;
    let sk = perm_pow(&sp_.sigma, kk);
    // Orbits of σ^k: representative and size.
    let mut orb_rep = vec![];
    let mut seen = [false; 8];
    for j in 0..8 {
        if seen[j] {
            continue;
        }
        let mut n = 0u64;
        let mut x = j;
        while !seen[x] {
            seen[x] = true;
            n += 1;
            x = sk[x];
        }
        orb_rep.push((j, n));
    }
    // (mask, stable, per-orbit characters, descent-trivial)
    let mut recs: Vec<(u32, bool, Vec<i8>, bool)> = vec![];
    let lce = k.c(lc);
    let fprime: Vec<E> = (0..8)
        .map(|j| {
            (0..8)
                .filter(|&i| i != j)
                .fold(lce, |a, i| k.mul(&a, &k.sub(&e[j], &e[i])))
        })
        .collect();
    let mut beta: Option<E> = None;
    let mut out = DescentK {
        k: kk,
        ..Default::default()
    };
    for mask in (0u32..256).step_by(2) {
        if mask.count_ones() % 2 == 1 {
            continue;
        }
        let comp = !mask & 0xff;
        let im = apply_mask(&sk, mask);
        let stable = im == mask;
        if !stable && im != comp {
            continue;
        }
        out.n2 += 1;
        if mask == 0 {
            out.halvable += 1;
            out.set |= 1;
            continue;
        }
        let mut chis = vec![];
        for &(j, nj) in &orb_rep {
            let u = if stable {
                let mut u = k.one();
                for i in 0..8 {
                    if i != j && mask >> i & 1 == 1 {
                        u = k.mul(&u, &k.sub(&e[i], &e[j]));
                    }
                }
                if mask >> j & 1 == 1 {
                    u = k.mul(&u, &fprime[j]);
                }
                u
            } else {
                // σ^k S = S^c: rational representative 2 D∞ - (R_1 + ... + R_4)
                // from y - β g_S(x) with β^{q+1} = lc.
                assert!(
                    m.is_multiple_of(2 * kk),
                    "swap case needs F_(p^2k) inside F_(p^m)"
                );
                let b = *beta.get_or_insert_with(|| find_beta(k, kk as u32, lc, rng));
                let b2 = k.mul(&b, &b);
                let gs = |set: u32| {
                    (0..8)
                        .filter(|&i| set >> i & 1 == 1)
                        .fold(k.one(), |a, i| k.mul(&a, &k.sub(&e[j], &e[i])))
                };
                let r = if mask >> j & 1 == 1 {
                    k.mul(&lce, &gs(comp))
                } else {
                    k.neg(&k.mul(&b2, &gs(mask)))
                };
                k.mul(&r, &k.inv(&k.sub(&lce, &b2)))
            };
            let n = kk * nj; // [K(e_j) : F_p]
            let g = gcd(m, n) as u32;
            match k.chi_sub(&u, g) {
                Some(c) if c != 0 => {
                    let c = if (n / g as u64).is_multiple_of(2) {
                        1
                    } else {
                        c
                    };
                    chis.push((c, nj));
                }
                _ => {
                    out.errors += 1;
                    chis.push((0, nj));
                }
            }
        }
        let all_one = chis.iter().all(|&(c, _)| c == 1);
        let all_twisted = chis
            .iter()
            .all(|&(c, nj)| c == if nj % 2 == 0 { 1 } else { -1 });
        if all_one || all_twisted {
            out.halvable += 1;
            out.set |= 1u128 << (mask >> 1);
        }
        recs.push((
            mask,
            stable,
            chis.iter().map(|c| c.0).collect(),
            all_one || all_twisted,
        ));
    }
    // Amendment 1: resolve the descent kernel on all-even factorizations.
    out.all_even = orb_rep.iter().all(|&(_, n)| n % 2 == 0);
    if !out.all_even {
        out.h_lo = out.halvable;
        out.h_hi = out.halvable;
        out.stable_true = out.set;
        return out;
    }
    let Some(y0) = recs.iter().find(|r| !r.1).cloned() else {
        out.errors += 1;
        return out;
    };
    let q3 = (k.p as u128).pow(kk as u32) % 4 == 3;
    let (mut sd, mut st, mut sw) = (1u32, 1u32, 0u32); // T = 0 is stable and halvable
    out.stable_true = 1;
    for (mask, stable, _, triv) in &recs {
        if !triv {
            continue;
        }
        if !stable {
            sw += 1;
            continue;
        }
        sd += 1;
        // B(y0, T) = ∏_{O ⊂ S_T} χ_O(δ_O(y0)); B(T, y0) = e_2(T, y0)^[q≡3] B(y0, T).
        let mut b: i8 = 1;
        for (o, &(j, _)) in orb_rep.iter().enumerate() {
            if mask >> j & 1 == 1 {
                b *= y0.2[o];
            }
        }
        if q3 && (mask & y0.0).count_ones() % 2 == 1 {
            b = -b;
        }
        if b == 1 {
            st += 1;
            out.stable_true |= 1u128 << (mask >> 1);
        } else if b != -1 {
            out.errors += 1;
        }
    }
    debug_assert_eq!(sd + sw, out.halvable);
    if sw == 0 {
        out.h_lo = st;
        out.h_hi = st;
    } else if st < sd {
        out.h_lo = out.halvable / 2;
        out.h_hi = out.halvable / 2;
    } else {
        out.h_lo = out.halvable / 2;
        out.h_hi = out.halvable;
    }
    out
}

/// Every relabeling of the roots commuting with σ (σ in block form).
fn centralizer(cycles: &[Vec<usize>]) -> Vec<[usize; 8]> {
    let mut by_len: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (c, cyc) in cycles.iter().enumerate() {
        by_len.entry(cyc.len()).or_default().push(c);
    }
    let mut perms: Vec<[usize; 8]> = vec![[0, 1, 2, 3, 4, 5, 6, 7]];
    for (l, cs) in by_len {
        let mut next = vec![];
        let mut order: Vec<usize> = (0..cs.len()).collect();
        let mut all_orders = vec![];
        permutations(&mut order, 0, &mut all_orders);
        let nrot = (l as u64).pow(cs.len() as u32);
        for base in &perms {
            for ord in &all_orders {
                for rot in 0..nrot {
                    let mut pm = *base;
                    let mut r = rot;
                    for (a, &src) in cs.iter().enumerate() {
                        let dst = cs[ord[a]];
                        let o = (r % l as u64) as usize;
                        r /= l as u64;
                        for i in 0..l {
                            pm[cycles[src][i]] = cycles[dst][(i + o) % l];
                        }
                    }
                    next.push(pm);
                }
            }
        }
        perms = next;
    }
    perms
}

fn permutations(v: &mut Vec<usize>, i: usize, out: &mut Vec<Vec<usize>>) {
    if i == v.len() {
        out.push(v.clone());
        return;
    }
    for j in i..v.len() {
        v.swap(i, j);
        permutations(v, i + 1, out);
        v.swap(i, j);
    }
}

fn canonical_sets(sets: &[u128], perms: &[[usize; 8]]) -> Vec<u128> {
    let mut best: Option<Vec<u128>> = None;
    for pm in perms {
        let t: Vec<u128> = sets
            .iter()
            .map(|&s| {
                let mut o = 0u128;
                for idx in 0..128u32 {
                    if s >> idx & 1 == 1 {
                        let mut m2 = apply_mask(pm, idx << 1);
                        if m2 & 1 == 1 {
                            m2 = !m2 & 0xff;
                        }
                        o |= 1u128 << (m2 >> 1);
                    }
                }
                o
            })
            .collect();
        if best.as_ref().is_none_or(|b| t < *b) {
            best = Some(t);
        }
    }
    best.unwrap()
}

fn divisors(n: u64) -> Vec<u64> {
    (1..=n).filter(|d| n.is_multiple_of(*d)).collect()
}

struct J4 {
    ct: String,
    /// `dim J(F_p)[2]`.
    r: u32,
    /// `(k, #J(F_{p^k})[2], descent-trivial count)` for `k | 2m`.
    prof: Vec<(u64, u32, u32)>,
    /// Corrected `(h_lo, h_hi)` per `k` (Amendment 1).
    corr: Vec<(u32, u32)>,
    x3: String,
    x3p: String,
    /// Canonical (descent-trivial set, stable halvable subset) per `k`.
    x3cp_sets: String,
    errors: u32,
}

fn j4(f: &[u64], lc: u64, p: u64, seed: u64) -> J4 {
    let mut rng = sp::SplitMix64(seed);
    let s = split(f, p, &mut rng).expect("squarefree");
    // Instrument: every root is a root of f.
    let mut errors = 0u32;
    for r in &s.roots {
        let v = f
            .iter()
            .rev()
            .fold([0u64; 16], |a, &c| s.k.add(&s.k.mul(&a, r), &s.k.c(c)));
        errors += !s.k.is_zero(&v) as u32;
    }
    let m = s.k.m as u64;
    let ds: Vec<DescentK> = divisors(2 * m)
        .into_iter()
        .map(|kk| descent(&s, lc, kk, &mut rng))
        .collect();
    errors += ds.iter().map(|d| d.errors).sum::<u32>();
    let ct = ct_string(&s.ct);
    let prof: Vec<(u64, u32, u32)> = ds.iter().map(|d| (d.k, d.n2, d.halvable)).collect();
    let x3 = format!(
        "{ct}|{}",
        prof.iter()
            .map(|&(k, n2, h)| format!("{k}:{}", n2 * h))
            .collect::<Vec<_>>()
            .join(",")
    );
    let perms = centralizer(&s.cycles);
    let sets: Vec<u128> = ds.iter().map(|d| d.set).collect();
    let can = canonical_sets(&sets, &perms);
    let x3p = format!(
        "{ct}|{}",
        can.iter()
            .map(|x| format!("{x:x}"))
            .collect::<Vec<_>>()
            .join(",")
    );
    let pairs: Vec<u128> = ds.iter().flat_map(|d| [d.set, d.stable_true]).collect();
    let canp = canonical_sets(&pairs, &perms);
    let x3cp_sets = canp
        .iter()
        .map(|x| format!("{x:x}"))
        .collect::<Vec<_>>()
        .join(",");
    J4 {
        ct,
        r: ds[0].n2.trailing_zeros(),
        corr: ds.iter().map(|d| (d.h_lo, d.h_hi)).collect(),
        prof,
        x3,
        x3p,
        x3cp_sets,
        errors,
    }
}

// ---------------------------------------------------------------------------
// 2-Sylow of J(F_p) from the class group of the function field.

struct Places {
    irr: Vec<P>,
    index: HashMap<P, usize>,
}

fn places(p: u64) -> Places {
    let mut irr = vec![];
    let mut rng = sp::SplitMix64(1);
    for d in 1..=3u32 {
        for code in 0..p.pow(d) {
            let mut a: P = (0..d).map(|i| code / p.pow(i) % p).collect();
            a.push(1);
            if d == 1 || sp::factor_squarefree(&a, p, &mut rng).is_some_and(|f| f.len() == 1) {
                irr.push(a);
            }
        }
    }
    let index = irr
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, a)| (a, i))
        .collect();
    Places { irr, index }
}

#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Ram,
    Split,
    Inert,
    Big,
}

/// Valuations of the Smith form over `Z/2^e`; returns (nonzero-pivot valuations, zero count).
fn snf_2adic(mut rows: Vec<Vec<u64>>, ncols: usize, e: u32) -> (Vec<u32>, usize) {
    let mask = if e >= 64 { u64::MAX } else { (1u64 << e) - 1 };
    let inv_odd = |a: u64| {
        let mut x = a; // Newton: x = x (2 - a x)
        for _ in 0..6 {
            x = x.wrapping_mul(2u64.wrapping_sub(a.wrapping_mul(x)));
        }
        x & mask
    };
    let mut col_done = vec![false; ncols];
    let mut alive: Vec<bool> = vec![true; rows.len()];
    let mut vals = vec![];
    // Phase 1: unit pivots.
    for c in 0..ncols {
        let Some(r) = (0..rows.len()).find(|&r| alive[r] && rows[r][c] & 1 == 1) else {
            continue;
        };
        alive[r] = false;
        col_done[c] = true;
        vals.push(0);
        let iv = inv_odd(rows[r][c]);
        let piv: Vec<u64> = rows[r].iter().map(|&x| x.wrapping_mul(iv) & mask).collect();
        for r2 in 0..rows.len() {
            if !alive[r2] || rows[r2][c] == 0 {
                continue;
            }
            let f = rows[r2][c];
            for j in 0..ncols {
                if piv[j] != 0 {
                    rows[r2][j] = rows[r2][j].wrapping_sub(f.wrapping_mul(piv[j])) & mask;
                }
            }
        }
    }
    // Phase 2: minimal-valuation pivots on what is left.
    loop {
        let mut best: Option<(u32, usize, usize)> = None;
        for r in 0..rows.len() {
            if !alive[r] {
                continue;
            }
            for c in 0..ncols {
                if col_done[c] || rows[r][c] == 0 {
                    continue;
                }
                let v = rows[r][c].trailing_zeros();
                if best.is_none_or(|b| v < b.0) {
                    best = Some((v, r, c));
                }
            }
        }
        let Some((v, r, c)) = best else { break };
        alive[r] = false;
        col_done[c] = true;
        vals.push(v);
        let iv = inv_odd(rows[r][c] >> v);
        let piv: Vec<u64> = rows[r].iter().map(|&x| x.wrapping_mul(iv) & mask).collect();
        for r2 in 0..rows.len() {
            if !alive[r2] || rows[r2][c] == 0 {
                continue;
            }
            let f = rows[r2][c] >> v;
            for j in 0..ncols {
                if piv[j] != 0 {
                    rows[r2][j] = rows[r2][j].wrapping_sub(f.wrapping_mul(piv[j])) & mask;
                }
            }
        }
    }
    let zeros = col_done.iter().filter(|&&d| !d).count();
    (vals.into_iter().filter(|&v| v > 0).collect(), zeros)
}

struct Sylow {
    /// Exponents `v_i` of `⊕ Z/2^{v_i}`, descending.
    v: Vec<u32>,
    certified: bool,
    rounds: usize,
    relations: usize,
    errors: u32,
}

fn sylow2(f: &[u64], t: u64, p: u64, e2: u32, pl: &Places, seed: u64) -> Sylow {
    let big_e = e2 + 1;
    let mut rng = sp::SplitMix64(seed);
    let ff: P = f.iter().map(|&c| c * t % p).collect(); // F = t f
    let split_inf = sp::chi(t, p) == 1;
    let n_inf = if split_inf { 2 } else { 1 };
    let mut kind = vec![Kind::Big; pl.irr.len()];
    let mut col = vec![usize::MAX; pl.irr.len()];
    let mut ncols = n_inf;
    for (i, a) in pl.irr.iter().enumerate() {
        let d = (a.len() - 1) as u32;
        let r = prem(&ff, a, p);
        kind[i] = if r.is_empty() {
            Kind::Ram
        } else if ppowmod(&r, ((p as u128).pow(d) - 1) / 2, a, p) == vec![1] {
            Kind::Split
        } else if d == 1 {
            Kind::Inert
        } else {
            Kind::Big
        };
        if matches!(kind[i], Kind::Ram | Kind::Split) {
            col[i] = ncols;
            ncols += 1;
        }
    }
    // D∞ as a vector.
    let add_dinf = |row: &mut Vec<i64>, c: i64| {
        row[0] += c;
        if split_inf {
            row[1] += c;
        }
    };
    let mut fixed: Vec<Vec<i64>> = vec![];
    for (i, a) in pl.irr.iter().enumerate() {
        if kind[i] == Kind::Ram {
            let mut row = vec![0i64; ncols];
            row[col[i]] = 2;
            add_dinf(&mut row, -((a.len() - 1) as i64));
            fixed.push(row);
        }
    }
    let mut root: HashMap<usize, P> = HashMap::new();
    let mut rels: Vec<Vec<i64>> = vec![];
    let mut errors = 0u32;
    let mut target = ncols + 30;
    let mut rounds = 0;
    let mut tries = 0u64;
    loop {
        rounds += 1;
        while rels.len() < target && tries < 50_000_000 {
            tries += 1;
            // Any v of degree <= 4: L(4 D∞) is spanned by 1, x, .., x^4, y.
            let mut v: P = (0..5).map(|_| rng.below(p)).collect();
            // v_4 = ±√t: y - v loses pole order at one point at infinity.
            let typ_b = split_inf && (v[4] == 1 || v[4] == p - 1);
            let eps_plus = v[4] == 1;
            trim(&mut v);
            let n = psub(&pmul(&v, &v, p), &ff, p);
            if n.is_empty() {
                continue;
            }
            let Some(fs) = sp::factor_squarefree(&n, p, &mut rng) else {
                continue;
            };
            if fs.iter().any(|g| g.len() > 4) {
                continue;
            }
            let mut row = vec![0i64; ncols];
            if typ_b {
                let dn = (n.len() - 1) as i64;
                let (ce, cm) = if eps_plus { (0, 1) } else { (1, 0) };
                row[ce] += 4 - dn;
                row[cm] += -4;
            } else {
                add_dinf(&mut row, -4);
            }
            let mut ok = true;
            for g in &fs {
                let i = pl.index[g];
                let d = (g.len() - 1) as i64;
                match kind[i] {
                    Kind::Ram => row[col[i]] += 1,
                    Kind::Split => {
                        let rv = prem(&v, g, p);
                        let r0 = root.entry(i).or_insert_with(|| rv.clone()).clone();
                        if rv == r0 {
                            row[col[i]] += 1;
                        } else if psub(&[], &rv, p) == r0 {
                            row[col[i]] -= 1;
                            add_dinf(&mut row, d);
                        } else {
                            ok = false;
                        }
                    }
                    _ => ok = false,
                }
            }
            if !ok {
                errors += 1;
                continue;
            }
            rels.push(row);
        }
        let modulus = 1i64 << big_e;
        let rows: Vec<Vec<u64>> = fixed
            .iter()
            .chain(rels.iter())
            .map(|r| r.iter().map(|&x| x.rem_euclid(modulus) as u64).collect())
            .collect();
        let (vals, zeros) = snf_2adic(rows, ncols, big_e);
        let sum: u32 = vals.iter().sum::<u32>() + zeros as u32 * big_e;
        if sum < big_e + e2 {
            errors += 1; // impossible: the quotient surjects onto Z/2^E ⊕ J[2^∞]
        }
        if zeros == 1 && vals.iter().sum::<u32>() == e2 {
            let mut v = vals;
            v.sort_unstable_by(|a, b| b.cmp(a));
            return Sylow {
                v,
                certified: true,
                rounds,
                relations: rels.len(),
                errors,
            };
        }
        if rounds >= MAX_ROUNDS || tries >= 50_000_000 {
            return Sylow {
                v: vec![],
                certified: false,
                rounds,
                relations: rels.len(),
                errors,
            };
        }
        target += 100;
    }
}

// ---------------------------------------------------------------------------
// L-polynomial arithmetic.

/// `v2(#J(F_{p^k}))` from `L` (exact mod 2^126).
fn v2_jk(a: [i64; 3], p: u64, k: u64) -> u32 {
    let p = p as i128;
    let c = [
        a[0] as i128,
        a[1] as i128,
        a[2] as i128,
        p * a[1] as i128,
        p * p * a[0] as i128,
        p * p * p,
    ];
    let n = (3 * k) as usize;
    let mut s = vec![0i128; n + 1];
    for i in 1..=n {
        let mut acc: i128 = 0;
        for j in 1..=6.min(i) {
            let sj = if j == i { i as i128 } else { s[i - j] };
            acc = acc.wrapping_add(c[j - 1].wrapping_mul(sj));
        }
        s[i] = acc.wrapping_neg();
    }
    let k = k as usize;
    let (s1, s2, s3) = (s[k], s[2 * k], s[3 * k]);
    let inv3: i128 = {
        let mut x: i128 = 3;
        for _ in 0..7 {
            x = x.wrapping_mul(2i128.wrapping_sub(3i128.wrapping_mul(x)));
        }
        x
    };
    // Newton: e1 = s1, 2 e2 = e1 s1 - s2, 3 e3 = e2 s1 - e1 s2 + s3.
    let e1 = s1;
    let e2 = (e1.wrapping_mul(s1).wrapping_sub(s2)) >> 1; // exact mod 2^127
    let e3 = e2
        .wrapping_mul(s1)
        .wrapping_sub(e1.wrapping_mul(s2))
        .wrapping_add(s3)
        .wrapping_mul(inv3);
    let q = p.wrapping_pow(k as u32);
    let nj = 1i128
        .wrapping_sub(e1)
        .wrapping_add(e2)
        .wrapping_sub(e3)
        .wrapping_add(q.wrapping_mul(e2))
        .wrapping_sub(q.wrapping_mul(q).wrapping_mul(e1))
        .wrapping_add(q.wrapping_mul(q).wrapping_mul(q));
    let nj = nj & ((1i128 << 126) - 1);
    assert!(nj != 0, "#J(F_p^k) vanishes mod 2^126");
    nj.trailing_zeros()
}

// ---------------------------------------------------------------------------
// Orbits (copied from smith_iso_granularity.rs).

fn shift(f: &[u64; 9], a: u64, p: u64) -> [u64; 9] {
    let mut g = *f;
    for i in 0..9 {
        for j in (i..8).rev() {
            g[j] = (g[j] + a * g[j + 1]) % p;
        }
    }
    g
}

fn orbit(code: u64, tw: bool, p: u64) -> Vec<(u64, bool)> {
    let f = decode(code, p);
    let inv8 = sp::inv_mod(8 % p, p);
    let mut out = Vec::new();
    let mut images = vec![f];
    for a in 0..p {
        let g = shift(&f, a, p);
        if g[0] != 0 {
            let mut h = g;
            h.reverse();
            images.push(h);
        }
    }
    for h in images {
        let lc = h[8];
        let tw2 = tw ^ (sp::chi(lc, p) == -1);
        let il = sp::inv_mod(lc, p);
        let mut m = h.map(|c| c * il % p);
        let a = (p - m[7]) * inv8 % p;
        m = shift(&m, a, p);
        for u in 1..p {
            let ui = sp::inv_mod(u, p);
            let mut s = 1u64;
            let mut code2 = 0u64;
            let mut coeffs = [0u64; 7];
            for i in (0..8).rev() {
                s = s * ui % p;
                if i < 7 {
                    coeffs[i] = m[i] * s % p;
                }
            }
            for i in (0..7).rev() {
                code2 = code2 * p + coeffs[i];
            }
            out.push((code2, tw2));
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn decode(code: u64, p: u64) -> [u64; 9] {
    let mut f = [0u64; 9];
    let mut v = code;
    for c in f.iter_mut().take(7) {
        *c = v % p;
        v /= p;
    }
    f[8] = 1;
    f
}

// ---------------------------------------------------------------------------
// compute

#[derive(Clone)]
struct Row {
    code: u64,
    tw: bool,
    key: [i64; 3],
    stratum: u8,
    full: bool,
    size: u64,
}

fn read_rows(path: &Path) -> Vec<Row> {
    let mut text = String::new();
    flate2::read::GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text.lines()
        .skip(1)
        .map(|line| {
            let v: Vec<i64> = line.split(',').map(|x| x.parse().unwrap()).collect();
            Row {
                code: v[0] as u64,
                tw: v[1] == 1,
                key: [v[2], v[3], v[4]],
                stratum: v[5] as u8,
                full: v[6] == 1,
                size: v[7] as u64,
            }
        })
        .collect()
}

struct Out {
    j4: J4,
    syl: Sylow,
    e2: u32,
    k2_j4: Option<bool>,
    k2_syl: Option<bool>,
    k5_fail: u32,
    /// Amendment 1: corrected X3c / X3c+ strings, K4' and K5 for X3c.
    x3c: String,
    x3cp: String,
    k4c_fail: bool,
    k4c_checked: bool,
    k5c_fail: u32,
    ambiguous_k: u32,
}

fn model(code: u64, tw: bool, p: u64) -> (P, u64) {
    let f = decode(code, p).to_vec();
    let t = if tw { sp::least_nonsquare(p) } else { 1 };
    (f, t)
}

fn compute_one(row: &Row, p: u64, pl: &Places, k2_syl: bool) -> Out {
    let (f, t) = model(row.code, row.tw, p);
    let seed = row.code * 2 + row.tw as u64;
    let j = j4(&f, t, p, seed);
    let nj_v2 = {
        let [a1, a2, a3] = row.key;
        let pi = p as i64;
        let nj = 1 + a1 + a2 + a3 + pi * a2 + pi * pi * a1 + pi * pi * pi;
        nj.trailing_zeros()
    };
    let syl = sylow2(&f, t, p, nj_v2, pl, seed);
    // Second orbit member (K2).
    let orb = orbit(row.code, row.tw, p);
    // The last orbit member other than the representative (if any).
    let &(c2, t2) = orb
        .iter()
        .rev()
        .find(|&&x| x != (row.code, row.tw))
        .unwrap_or(&(row.code, row.tw));
    let (f2, tt2) = model(c2, t2, p);
    let seed2 = c2 * 2 + t2 as u64;
    let j2 = j4(&f2, tt2, p, seed2);
    let k2_j4 = Some(
        j2.ct == j.ct
            && j2.x3 == j.x3
            && j2.x3p == j.x3p
            && j2.corr == j.corr
            && j2.x3cp_sets == j.x3cp_sets,
    );
    let k2_s = if k2_syl {
        let s2 = sylow2(&f2, tt2, p, nj_v2, pl, seed2);
        Some(s2.certified && s2.v == syl.v)
    } else {
        None
    };
    // K5.
    let mut k5_fail = 0;
    for &(k, n2, h) in &j.prof {
        let v = v2_jk(row.key, p, k);
        let l4 = (n2 * h).trailing_zeros();
        let l2 = n2.trailing_zeros();
        if l4 > v || (v == l2 && l4 != l2) {
            k5_fail += 1;
        }
    }
    // Amendment 1: X3c / X3c+ with the k = 1 count from the class group.
    let n4 = syl.v.iter().filter(|&&v| v >= 2).count() as u32;
    let truth1 = syl.certified.then_some(1u32 << n4);
    let (lo1, hi1) = j.corr[0];
    let k4c_checked = truth1.is_some();
    let k4c_fail = match truth1 {
        Some(t) => !(t == lo1 || t == hi1) || (lo1 == hi1 && lo1 != t),
        None => false,
    };
    let mut ambiguous_k = 0;
    let mut k5c_fail = 0;
    let mut entries = vec![];
    for (i, &(k, n2, _)) in j.prof.iter().enumerate() {
        let (lo, hi) = if i == 0 {
            truth1.map(|t| (t, t)).unwrap_or((lo1, hi1))
        } else {
            j.corr[i]
        };
        let v = v2_jk(row.key, p, k);
        let l2 = n2.trailing_zeros();
        let ok = |h: u32| {
            let l4 = (n2 * h).trailing_zeros();
            !(l4 > v || (v == l2 && l4 != l2))
        };
        if !(ok(lo) || ok(hi)) {
            k5c_fail += 1;
        }
        if lo == hi {
            entries.push(format!("{k}:{}", n2 * lo));
        } else {
            ambiguous_k += 1;
            entries.push(format!("{k}:{}~{}", n2 * lo, n2 * hi));
        }
    }
    let x3c = format!("{}|{}", j.ct, entries.join(","));
    let x3cp = format!(
        "{}|{}|{}",
        j.ct,
        truth1.map(|t| t.to_string()).unwrap_or("?".into()),
        j.x3cp_sets
    );
    Out {
        j4: j,
        syl,
        e2: nj_v2,
        k2_j4,
        k2_syl: k2_s,
        k5_fail,
        x3c,
        x3cp,
        k4c_fail,
        k4c_checked,
        k5c_fail,
        ambiguous_k,
    }
}

fn compute_cmd(args: &[String]) {
    let p: u64 = arg(args, "--p").expect("--p").parse().unwrap();
    let out = PathBuf::from(arg(args, "--out").expect("--out DIR"));
    let limit: Option<usize> = arg(args, "--limit").map(|s| s.parse().unwrap());
    fs::create_dir_all(&out).unwrap();
    let t0 = Instant::now();
    let rows = read_rows(&Path::new(FROZEN).join(format!("iso_p{p}.csv.gz")));
    let pl = places(p);
    let mut tract: Vec<usize> = (0..rows.len()).filter(|&i| rows[i].stratum != 0).collect();
    if let Some(l) = limit {
        tract.truncate(l);
    }
    // K2 sample for the class group: first K2_SAMPLE of a seeded shuffle.
    let mut order = tract.clone();
    let mut rng = sp::SplitMix64(SEED);
    for a in (1..order.len()).rev() {
        let b = rng.below(a as u64 + 1) as usize;
        order.swap(a, b);
    }
    let sample: BTreeSet<usize> = order.into_iter().take(K2_SAMPLE).collect();
    let done = std::sync::atomic::AtomicUsize::new(0);
    let n_tr = tract.len();
    let outs: Vec<(usize, Out)> = tract
        .par_iter()
        .map(|&i| {
            let o = compute_one(&rows[i], p, &pl, sample.contains(&i));
            let d = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
            if d.is_multiple_of(5000) {
                eprintln!(
                    "p={p}: {d}/{n_tr} classes, {:.0}s",
                    t0.elapsed().as_secs_f64()
                );
            }
            (i, o)
        })
        .collect();
    let secs = t0.elapsed().as_secs_f64();
    let by_idx: HashMap<usize, &Out> = outs.iter().map(|(i, o)| (*i, o)).collect();

    let mut csv = String::from(
        "rep_code,rep_twist,a1,a2,a3,stratum,full,size,ct,r,e2,sylow,certified,rounds,relations,x3,x3p,x3c,x3cp\n",
    );
    let mut c = BTreeMap::<&str, u64>::new();
    for (i, r) in rows.iter().enumerate() {
        let base = format!(
            "{},{},{},{},{},{},{},{}",
            r.code, r.tw as u8, r.key[0], r.key[1], r.key[2], r.stratum, r.full as u8, r.size
        );
        let Some(o) = by_idx.get(&i) else {
            csv.push_str(&format!("{base},,,,,,,,,,,\n"));
            continue;
        };
        let syl = o
            .syl
            .v
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(".");
        csv.push_str(&format!(
            "{base},{},{},{},{},{},{},{},{},{},{},{}\n",
            o.j4.ct,
            o.j4.r,
            o.e2,
            if syl.is_empty() { "0".into() } else { syl },
            o.syl.certified as u8,
            o.syl.rounds,
            o.syl.relations,
            o.j4.x3,
            o.j4.x3p,
            o.x3c,
            o.x3cp
        ));
        *c.entry("tractable").or_default() += 1;
        *c.entry("j4_instrument_errors").or_default() += o.j4.errors as u64;
        *c.entry("sylow_instrument_errors").or_default() += o.syl.errors as u64;
        *c.entry("k2_j4_checked").or_default() += o.k2_j4.is_some() as u64;
        *c.entry("k2_j4_mismatch").or_default() += (o.k2_j4 == Some(false)) as u64;
        *c.entry("k2_sylow_checked").or_default() += o.k2_syl.is_some() as u64;
        *c.entry("k2_sylow_mismatch").or_default() += (o.k2_syl == Some(false)) as u64;
        *c.entry("k3_uncertified").or_default() += (!o.syl.certified) as u64;
        if o.syl.certified {
            let sum: u32 = o.syl.v.iter().sum();
            *c.entry("k3_sum_mismatch").or_default() += (sum != o.e2) as u64;
            *c.entry("k3_rank_mismatch").or_default() += (o.syl.v.len() as u32 != o.j4.r) as u64;
            let n4 = o.syl.v.iter().filter(|&&v| v >= 2).count() as u32;
            let h1 = o.j4.prof[0].2;
            *c.entry("k4_checked").or_default() += 1;
            *c.entry("k4_mismatch").or_default() += (h1 != 1 << n4) as u64;
        }
        *c.entry("k5_failures").or_default() += o.k5_fail as u64;
        *c.entry("k4c_checked").or_default() += o.k4c_checked as u64;
        *c.entry("k4c_mismatch").or_default() += o.k4c_fail as u64;
        *c.entry("k5c_failures").or_default() += o.k5c_fail as u64;
        *c.entry("x3c_classes_with_ambiguous_k").or_default() += (o.ambiguous_k > 0) as u64;
        *c.entry("x3c_ambiguous_entries").or_default() += o.ambiguous_k as u64;
    }
    let mut gz = GzEncoder::new(Vec::new(), Compression::best());
    gz.write_all(csv.as_bytes()).unwrap();
    fs::write(out.join(format!("inv_p{p}.csv.gz")), gz.finish().unwrap()).unwrap();
    let z = |k: &str| c.get(k).copied().unwrap_or(0);
    let git = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default();
    let ctl = json!({
        "p": p,
        "git_head": git,
        "iso_classes": rows.len(),
        "tractable_classes": z("tractable"),
        "limit": limit,
        "counts": c,
        "k2_pass_j4": z("k2_j4_mismatch") == 0 && z("k2_j4_checked") == z("tractable"),
        "k2_pass_sylow": z("k2_sylow_mismatch") == 0,
        "k3_pass": z("k3_uncertified") == 0 && z("k3_sum_mismatch") == 0 && z("k3_rank_mismatch") == 0
            && z("sylow_instrument_errors") == 0,
        "k4_pass": z("k4_mismatch") == 0 && z("k4_checked") == z("tractable") && z("j4_instrument_errors") == 0,
        "k5_pass": z("k5_failures") == 0,
        "k4c_pass": z("k4c_mismatch") == 0 && z("k4c_checked") == z("tractable") && z("j4_instrument_errors") == 0,
        "k5c_pass": z("k5c_failures") == 0,
        "seconds": secs,
        "threads": rayon::current_num_threads(),
    });
    fs::write(
        out.join(format!("controls_p{p}.json")),
        serde_json::to_string_pretty(&ctl).unwrap() + "\n",
    )
    .unwrap();
    println!("{}", serde_json::to_string_pretty(&ctl).unwrap());
}

// ---------------------------------------------------------------------------
// summarize

struct Iso {
    lid: Vec<usize>,
    keys: Vec<[i64; 3]>,
    stratum: Vec<u8>,
    full: Vec<bool>,
    size: Vec<u64>,
    n_l: usize,
    inv: Vec<[String; 6]>,
}

fn load(dir: &Path, p: u64) -> Iso {
    let mut text = String::new();
    flate2::read::GzDecoder::new(fs::File::open(dir.join(format!("inv_p{p}.csv.gz"))).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    let mut ids: HashMap<[i64; 3], usize> = HashMap::new();
    let mut iso = Iso {
        lid: vec![],
        keys: vec![],
        stratum: vec![],
        full: vec![],
        size: vec![],
        n_l: 0,
        inv: vec![],
    };
    for line in text.lines().skip(1) {
        let v: Vec<&str> = line.split(',').collect();
        let n = |i: usize| v[i].parse::<i64>().unwrap();
        let k = [n(2), n(3), n(4)];
        let len = ids.len();
        iso.lid.push(*ids.entry(k).or_insert(len));
        iso.keys.push(k);
        iso.stratum.push(n(5) as u8);
        iso.full.push(n(6) == 1);
        iso.size.push(n(7) as u64);
        iso.inv.push([
            v[8].to_string(),
            v[11].to_string(),
            v[15].to_string(),
            v[16].to_string(),
            v[17].to_string(),
            v[18].to_string(),
        ]);
    }
    iso.n_l = ids.len();
    iso
}

fn class_share(iso: &Iso, full: &[bool], total: u64) -> f64 {
    let mut hit = vec![false; iso.n_l];
    for (i, &f) in full.iter().enumerate() {
        if f {
            hit[iso.lid[i]] = true;
        }
    }
    let s: u64 = (0..full.len())
        .filter(|&i| hit[iso.lid[i]])
        .map(|i| iso.size[i])
        .sum();
    s as f64 / total as f64
}

/// Identical to `smith_iso_granularity::null` (K1 depends on it).
fn null(iso: &Iso, group: &[u64], total: u64, rng: &mut sp::SplitMix64) -> Vec<f64> {
    let mut groups: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    for (i, &g) in group.iter().enumerate() {
        groups.entry(g).or_default().push(i);
    }
    let mut out = Vec::with_capacity(PERMS);
    let mut f = iso.full.clone();
    for _ in 0..PERMS {
        for idx in groups.values() {
            let mut lab: Vec<bool> = idx.iter().map(|&i| iso.full[i]).collect();
            for a in (1..lab.len()).rev() {
                let b = rng.below(a as u64 + 1) as usize;
                lab.swap(a, b);
            }
            for (j, &i) in idx.iter().enumerate() {
                f[i] = lab[j];
            }
        }
        out.push(class_share(iso, &f, total));
    }
    out.sort_by(|a, b| a.partial_cmp(b).unwrap());
    out
}

fn quant(v: &[f64], x: f64) -> f64 {
    let h = x * (v.len() - 1) as f64;
    let (lo, hi) = (h.floor() as usize, h.ceil() as usize);
    v[lo] + (h - lo as f64) * (v[hi] - v[lo])
}

fn null_json(ns: &[f64], c: f64) -> Value {
    json!({
        "perms": PERMS, "seed": SEED,
        "median": quant(ns, 0.5), "lo": quant(ns, 0.025), "hi": quant(ns, 0.975),
        "frac_le_observed": ns.iter().filter(|&&x| x <= c).count() as f64 / ns.len() as f64,
    })
}

const CANDIDATES: [(&str, usize); 6] = [
    ("X1_cycle_type", 0),
    ("X2_sylow2", 1),
    ("X3_j4_profile", 2),
    ("X3plus_j4_sets", 3),
    ("X3c_j4_profile_corrected", 4),
    ("X3cplus_j4_descent_sets_corrected", 5),
];

fn summarize(dir: &Path, p: u64) -> Value {
    let ctl: Value =
        serde_json::from_str(&fs::read_to_string(dir.join(format!("controls_p{p}.json"))).unwrap())
            .unwrap();
    let frozen: Value = serde_json::from_str(
        &fs::read_to_string(Path::new(FROZEN).join(format!("summary_p{p}.json"))).unwrap(),
    )
    .unwrap();
    let iso = load(dir, p);
    let n = iso.size.len();
    let total: u64 = iso.size.iter().sum();
    let c = class_share(&iso, &iso.full, total);
    // K1.
    let g0: Vec<u64> = iso.stratum.iter().map(|&s| s as u64).collect();
    let n0 = null(&iso, &g0, total, &mut sp::SplitMix64(SEED));
    let fz = &frozen["null_primary"];
    let k1 = frozen["C"].as_f64() == Some(c)
        && fz["lo"].as_f64() == Some(quant(&n0, 0.025))
        && fz["median"].as_f64() == Some(quant(&n0, 0.5))
        && fz["hi"].as_f64() == Some(quant(&n0, 0.975));
    let med0 = quant(&n0, 0.5);
    let valid = |cand: usize| -> bool {
        let b = |k: &str| ctl[k].as_bool() == Some(true);
        k1 && match cand {
            0 => b("k2_pass_j4"),
            1 => b("k2_pass_sylow") && b("k3_pass"),
            2 | 3 => b("k2_pass_j4") && b("k4_pass") && b("k5_pass"),
            _ => b("k2_pass_j4") && b("k4c_pass") && b("k5c_pass") && b("k3_pass"),
        }
    };
    let mut res = BTreeMap::new();
    for (name, ci) in CANDIDATES {
        // Cells (stratum, X).
        let mut cell_id: HashMap<(u8, String), u64> = HashMap::new();
        let group: Vec<u64> = (0..n)
            .map(|i| {
                if iso.stratum[i] == 0 {
                    return 0;
                }
                let key = (iso.stratum[i], iso.inv[i][ci].clone());
                let len = cell_id.len() as u64 + 1;
                *cell_id.entry(key).or_insert(len)
            })
            .collect();
        let mut cell_n: HashMap<u64, (u64, u64)> = HashMap::new();
        for i in 0..n {
            if group[i] != 0 {
                let e = cell_n.entry(group[i]).or_default();
                e.0 += 1;
                e.1 += iso.full[i] as u64;
            }
        }
        let tractable = cell_n.values().map(|e| e.0).sum::<u64>();
        let singleton = cell_n.values().filter(|e| e.0 == 1).count() as u64;
        let degenerate = cell_n
            .values()
            .filter(|e| e.1 == 0 || e.1 == e.0)
            .map(|e| e.0)
            .sum::<u64>();
        let nx = null(&iso, &group, total, &mut sp::SplitMix64(SEED));
        // Placebo: L-class buckets with n_s values per stratum.
        let mut ns: BTreeMap<u8, BTreeSet<&str>> = BTreeMap::new();
        for i in 0..n {
            if iso.stratum[i] != 0 {
                ns.entry(iso.stratum[i])
                    .or_default()
                    .insert(&iso.inv[i][ci]);
            }
        }
        let pgroup: Vec<u64> = (0..n)
            .map(|i| {
                let s = iso.stratum[i];
                if s == 0 {
                    return 0;
                }
                let k = iso.keys[i];
                let h = sp::SplitMix64(
                    PLACEBO_SEED
                        ^ (k[0] as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
                        ^ (k[1] as u64).wrapping_mul(0xC2B2_AE3D_27D4_EB4F)
                        ^ (k[2] as u64).wrapping_mul(0x1656_67B1_9E37_79F9)
                        ^ (s as u64) << 56,
                )
                .next_u64();
                let b = h % ns[&s].len() as u64;
                1 + s as u64 * 1_000_000 + b
            })
            .collect();
        let np = null(&iso, &pgroup, total, &mut sp::SplitMix64(SEED));
        let (lo, hi, med) = (quant(&nx, 0.025), quant(&nx, 0.975), quant(&nx, 0.5));
        let gap = (med0 - med) / (med0 - c);
        let verdict = if !valid(ci) {
            "invalid"
        } else if c > hi {
            "over-explains"
        } else if c >= lo {
            if c < quant(&np, 0.025) {
                "explains"
            } else {
                "confounded with L-class"
            }
        } else if gap >= 0.5 {
            "partial"
        } else {
            "does not explain"
        };
        res.insert(
            name.to_string(),
            json!({
                "cells": cell_n.len(),
                "values_per_stratum": ns.iter().map(|(s, v)| (sp::STRATA[*s as usize].to_string(), v.len())).collect::<BTreeMap<_,_>>(),
                "tractable_iso_classes": tractable,
                "singleton_cells": singleton,
                "share_in_singleton_cells": singleton as f64 / tractable as f64,
                "share_in_degenerate_cells": degenerate as f64 / tractable as f64,
                "null": null_json(&nx, c),
                "placebo": null_json(&np, c),
                "gap_share": gap,
                "verdict": verdict,
            }),
        );
    }
    // Descriptive: success by cycle type.
    let mut by_ct: BTreeMap<String, [u64; 2]> = BTreeMap::new();
    for i in 0..n {
        if iso.stratum[i] != 0 {
            let e = by_ct
                .entry(format!(
                    "s{}|{}",
                    sp::STRATA[iso.stratum[i] as usize],
                    iso.inv[i][0]
                ))
                .or_default();
            e[0] += 1;
            e[1] += iso.full[i] as u64;
        }
    }
    let mut by_syl: BTreeMap<String, [u64; 2]> = BTreeMap::new();
    for i in 0..n {
        if iso.stratum[i] != 0 {
            let e = by_syl
                .entry(format!(
                    "s{}|{}",
                    sp::STRATA[iso.stratum[i] as usize],
                    iso.inv[i][1]
                ))
                .or_default();
            e[0] += 1;
            e[1] += iso.full[i] as u64;
        }
    }
    json!({
        "p": p,
        "controls": ctl,
        "k1_reproduces_frozen_null": k1,
        "C": c,
        "null_unconditioned": null_json(&n0, c),
        "candidates": res,
        "descriptive_by_stratum_and_cycle_type": by_ct,
        "descriptive_by_stratum_and_sylow2": by_syl,
    })
}

fn summarize_cmd(args: &[String]) {
    let dir = PathBuf::from(arg(args, "--dir").expect("--dir DIR"));
    let ps: Vec<u64> = args[4..].iter().map(|s| s.parse().unwrap()).collect();
    let mut per: BTreeMap<String, BTreeMap<String, Value>> = BTreeMap::new();
    for p in ps {
        let s = summarize(&dir, p);
        fs::write(
            dir.join(format!("summary_p{p}.json")),
            serde_json::to_string_pretty(&s).unwrap() + "\n",
        )
        .unwrap();
        for (name, _) in CANDIDATES {
            let v = s["candidates"][name]["verdict"].clone();
            eprintln!(
                "p={p} {name}: C={:.4} null [{:.4}, {:.4}] placebo lo {:.4} gap {:.3} -> {}",
                s["C"],
                s["candidates"][name]["null"]["lo"],
                s["candidates"][name]["null"]["hi"],
                s["candidates"][name]["placebo"]["lo"],
                s["candidates"][name]["gap_share"],
                v
            );
            per.entry(name.to_string())
                .or_default()
                .insert(p.to_string(), v);
        }
    }
    let overall: BTreeMap<String, Value> = per
        .iter()
        .map(|(name, m)| {
            let vs: BTreeSet<String> = m.values().map(|v| v.to_string()).collect();
            let o = if vs.len() == 1 {
                m.values().next().unwrap().clone()
            } else {
                json!("prime-dependent")
            };
            (name.clone(), json!({"per_prime": m, "overall": o}))
        })
        .collect();
    fs::write(
        dir.join("verdict.json"),
        serde_json::to_string_pretty(&overall).unwrap() + "\n",
    )
    .unwrap();
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("compute") => compute_cmd(&args),
        Some("summarize") => summarize_cmd(&args),
        _ => eprintln!("usage: compute --p P --out DIR [--limit N] | summarize --dir DIR P..."),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roots_of_an_irreducible_quintic_and_cubic() {
        let p: u64 = 7;
        let mut rng = sp::SplitMix64(3);
        // Find a model whose factor type is 3+5.
        for code in 0..p.pow(7) {
            let f = decode(code, p).to_vec();
            if let Some(fs) = sp::factor_squarefree(&f, p, &mut rng) {
                let mut d: Vec<usize> = fs.iter().map(|g| g.len() - 1).collect();
                d.sort();
                if d == vec![3, 5] {
                    let s = split(&f, p, &mut rng).unwrap();
                    assert_eq!(s.k.m, 15);
                    for r in &s.roots {
                        let v = f
                            .iter()
                            .rev()
                            .fold([0u64; 16], |a, &c| s.k.add(&s.k.mul(&a, r), &s.k.c(c)));
                        assert!(s.k.is_zero(&v));
                    }
                    return;
                }
            }
        }
        panic!("no 3+5 model");
    }

    #[test]
    fn v2_jk_matches_l1_at_k1() {
        let a = [3i64, -2, 5];
        let p = 7i64;
        let nj = 1 + a[0] + a[1] + a[2] + p * a[1] + p * p * a[0] + p * p * p;
        assert_eq!(v2_jk(a, 7, 1), nj.trailing_zeros());
    }

    #[test]
    fn sylow_and_descent_agree_on_the_first_tractable_classes() {
        let p: u64 = 7;
        let rows = read_rows(&Path::new(FROZEN).join("iso_p7.csv.gz"));
        let pl = places(p);
        for r in rows.iter().filter(|r| r.stratum != 0).take(40) {
            let o = compute_one(r, p, &pl, true);
            assert!(o.syl.certified);
            assert_eq!(o.syl.v.iter().sum::<u32>(), o.e2);
            assert_eq!(o.syl.v.len() as u32, o.j4.r);
            // K4' (Amendment 1): the corrected k = 1 count brackets the class group.
            assert!(!o.k4c_fail);
            assert_eq!(o.k2_j4, Some(true));
            assert_eq!(o.k2_syl, Some(true));
            assert_eq!(o.k5c_fail, 0);
            assert_eq!(o.j4.errors + o.syl.errors, 0);
        }
    }
}
