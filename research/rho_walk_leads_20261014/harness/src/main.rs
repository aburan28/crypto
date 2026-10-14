//! Pollard rho walk leads (2026-10-14): one harness for four proposals.
//!
//!   e1   selection walks: step to the best-scoring of k candidates
//!        (P+R vs P−R and variants) against the plain r-adding walk, single
//!        walk with full storage, steps to the first repeated point; plus an
//!        abstract random-function control with the same selection rule.
//!   e2   lazy reduction under the distinguished-point method: index or DP
//!        test computed on a redundant (x or x+p) representation that depends
//!        on how the point was produced; also a time-shared index sequence
//!        (correlated walks) and a precomputed pair table R_i+R_j.
//!   e4   negation map with the class representative rule fixed, chosen per
//!        point, or switched with a global epoch.
//!   bench  batched-affine step cost (Montgomery trick over m walks) for the
//!        plain step, the P±R best-of-two step and the best-of-four step.
//!   summary  read the JSON lines the other commands wrote and print, per
//!        cell, means, paired ratios against the reference arm with
//!        bootstrap intervals, and the predicted value where one exists.
//!
//! Curves: y^2 = x^3 + ax + b over F_p, p = 3 mod 4, prime group order N found
//! by baby-step giant-step in the Hasse interval. Every arm of a trial uses
//! the same curve, secret, table and start (paired by trial seed). Every
//! collision is turned into a discrete log and checked against the secret.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufRead, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

// ---------------------------------------------------------------- utilities

#[inline]
fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        mix(self.0)
    }
    fn below(&mut self, n: u64) -> u64 {
        ((self.next() as u128 * n as u128) >> 64) as u64
    }
}

#[derive(Default)]
struct FxH(u64);
impl Hasher for FxH {
    fn finish(&self) -> u64 {
        mix(self.0)
    }
    fn write(&mut self, b: &[u8]) {
        for c in b.chunks(8) {
            let mut w = [0u8; 8];
            w[..c.len()].copy_from_slice(c);
            self.write_u64(u64::from_le_bytes(w));
        }
    }
    fn write_u64(&mut self, x: u64) {
        self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(0x517C_C1B7_2722_0A95);
    }
}
type FastMap<K, V> = HashMap<K, V, BuildHasherDefault<FxH>>;

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    (a as u128 * b as u128 % m as u128) as u64
}
fn powmod(mut b: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1 % m;
    b %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, b, m);
        }
        b = mulmod(b, b, m);
        e >>= 1;
    }
    r
}
fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for q in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % q == 0 {
            return n == q;
        }
    }
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
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
fn isqrt(n: u64) -> u64 {
    let mut x = (n as f64).sqrt() as u64;
    while x * x > n {
        x -= 1;
    }
    while (x + 1) * (x + 1) <= n {
        x += 1;
    }
    x
}

// ------------------------------------------------------- Montgomery F_p

#[derive(Clone, Copy)]
struct Fp {
    p: u64,
    ninv: u64,
    r2: u64,
    one: u64,
}
impl Fp {
    fn new(p: u64) -> Fp {
        let mut inv = p;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        }
        let r = ((1u128 << 64) % p as u128) as u64;
        Fp { p, ninv: inv.wrapping_neg(), r2: mulmod(r, r, p), one: r }
    }
    #[inline(always)]
    fn mul(&self, a: u64, b: u64) -> u64 {
        let t = a as u128 * b as u128;
        let m = (t as u64).wrapping_mul(self.ninv);
        let u = ((t + m as u128 * self.p as u128) >> 64) as u64;
        if u >= self.p { u - self.p } else { u }
    }
    #[inline(always)]
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p { s - self.p } else { s }
    }
    #[inline(always)]
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b { a - b } else { a + self.p - b }
    }
    #[inline(always)]
    fn neg(&self, a: u64) -> u64 {
        if a == 0 { 0 } else { self.p - a }
    }
    fn to_m(&self, a: u64) -> u64 {
        self.mul(a % self.p, self.r2)
    }
    fn pow(&self, mut b: u64, mut e: u64) -> u64 {
        let mut r = self.one;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, b);
            }
            b = self.mul(b, b);
            e >>= 1;
        }
        r
    }
    fn inv(&self, a: u64) -> u64 {
        self.pow(a, self.p - 2)
    }
}

// ------------------------------------------------------------------ curve

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Pt {
    x: u64,
    y: u64,
    inf: bool,
}
const INF: Pt = Pt { x: 0, y: 0, inf: true };

#[derive(Clone, Copy)]
struct Curve {
    f: Fp,
    a: u64,
    b: u64,
    n: u64,
    g: Pt,
    bits: u32,
}
impl Curve {
    fn rhs(&self, x: u64) -> u64 {
        let f = &self.f;
        f.add(f.mul(f.add(f.mul(x, x), self.a), x), self.b)
    }
    fn neg(&self, p: Pt) -> Pt {
        if p.inf { p } else { Pt { x: p.x, y: self.f.neg(p.y), inf: false } }
    }
    fn dbl(&self, p: Pt) -> Pt {
        let f = &self.f;
        if p.inf || p.y == 0 {
            return INF;
        }
        let x2 = f.mul(p.x, p.x);
        let num = f.add(f.add(f.add(x2, x2), x2), self.a);
        let l = f.mul(num, f.inv(f.add(p.y, p.y)));
        let x3 = f.sub(f.sub(f.mul(l, l), p.x), p.x);
        let y3 = f.sub(f.mul(l, f.sub(p.x, x3)), p.y);
        Pt { x: x3, y: y3, inf: false }
    }
    fn add(&self, p: Pt, q: Pt) -> Pt {
        let f = &self.f;
        if p.inf {
            return q;
        }
        if q.inf {
            return p;
        }
        if p.x == q.x {
            return if p.y == q.y { self.dbl(p) } else { INF };
        }
        let l = f.mul(f.sub(q.y, p.y), f.inv(f.sub(q.x, p.x)));
        let x3 = f.sub(f.sub(f.mul(l, l), p.x), q.x);
        let y3 = f.sub(f.mul(l, f.sub(p.x, x3)), p.y);
        Pt { x: x3, y: y3, inf: false }
    }
    fn mul(&self, mut k: u64, p: Pt) -> Pt {
        let (mut r, mut b) = (INF, p);
        while k > 0 {
            if k & 1 == 1 {
                r = self.add(r, b);
            }
            b = self.dbl(b);
            k >>= 1;
        }
        r
    }
    fn random_point(&self, rng: &mut Rng) -> Pt {
        let f = &self.f;
        loop {
            let x = f.to_m(rng.below(f.p));
            let r = self.rhs(x);
            if r == 0 {
                continue;
            }
            let y = f.pow(r, (f.p + 1) / 4);
            if f.mul(y, y) == r {
                return Pt { x, y, inf: false };
            }
        }
    }
    /// Unique multiple m of P in the Hasse interval with mP = O, if prime.
    fn prime_order_of(&self, p: Pt) -> Option<u64> {
        let q = self.f.p;
        let s = isqrt(q) + 1;
        let lo = q + 1 - 2 * s;
        let w = 4 * s + 1;
        let bstep = isqrt(w) + 1;
        let mut tab: FastMap<u64, (u64, u64)> = FastMap::default();
        let mut jp = p;
        for j in 1..=bstep {
            tab.insert(jp.x, (j, jp.y));
            jp = self.add(jp, p);
        }
        let giant = self.mul(bstep, p);
        let mut qk = self.mul(lo, p);
        let mut cands = Vec::new();
        for k in 0..=(w / bstep + 1) {
            let base = lo + k * bstep;
            if qk.inf {
                cands.push(base);
            } else if let Some(&(j, y)) = tab.get(&qk.x) {
                if y == self.f.neg(qk.y) {
                    cands.push(base + j);
                }
                if y == qk.y && base >= j {
                    cands.push(base - j);
                }
            }
            qk = self.add(qk, giant);
        }
        cands.retain(|&m| m >= lo && m < lo + w);
        cands.sort();
        cands.dedup();
        if cands.len() == 1 && is_prime(cands[0]) { Some(cands[0]) } else { None }
    }
    fn generate(bits: u32, seed: u64) -> Curve {
        let mut rng = Rng(seed ^ 0xC0FFEE);
        loop {
            let p = loop {
                let c = (1u64 << (bits - 1)) | rng.below(1u64 << (bits - 1)) | 3;
                if is_prime(c) {
                    break c;
                }
            };
            let f = Fp::new(p);
            for _ in 0..400 {
                let a = f.to_m(rng.below(p));
                let b = f.to_m(rng.below(p));
                // 4a^3 + 27b^2 != 0
                let d = f.add(
                    f.mul(f.to_m(4), f.mul(a, f.mul(a, a))),
                    f.mul(f.to_m(27), f.mul(b, b)),
                );
                if d == 0 {
                    continue;
                }
                let mut c = Curve { f, a, b, n: 0, g: INF, bits };
                let g = c.random_point(&mut rng);
                if let Some(n) = c.prime_order_of(g) {
                    c.n = n;
                    c.g = g;
                    assert!(c.mul(n, g).inf);
                    return c;
                }
            }
        }
    }
}

// ------------------------------------------------------- scalars mod N

fn sadd(a: u64, b: u64, n: u64) -> u64 {
    ((a as u128 + b as u128) % n as u128) as u64
}
fn sneg(a: u64, n: u64) -> u64 {
    if a == 0 { 0 } else { n - a }
}
fn solve(c: &Curve, a1: u64, b1: u64, a2: u64, b2: u64) -> Option<u64> {
    // a1 + b1 k = a2 + b2 k (mod N)
    let n = c.n;
    let db = sadd(b2, sneg(b1, n), n);
    if db == 0 {
        return None;
    }
    Some(mulmod(sadd(a1, sneg(a2, n), n), powmod(db, n - 2, n), n))
}

// ------------------------------------------------------------- instance

#[derive(Clone, Copy)]
struct W {
    p: Pt,
    a: u64,
    b: u64,
}

struct Inst {
    c: Curve,
    k: u64,
    q: Pt,
    tab: Vec<W>,
    keys: [u64; 8],
    rng: Rng,
}
impl Inst {
    fn new(c: &Curve, seed: u64, r: usize) -> Inst {
        let mut rng = Rng(mix(seed) ^ 0x5EED);
        let k = 1 + rng.below(c.n - 1);
        let q = c.mul(k, c.g);
        let mut keys = [0u64; 8];
        for x in keys.iter_mut() {
            *x = rng.next();
        }
        let mut s = Inst { c: *c, k, q, tab: Vec::new(), keys, rng };
        for _ in 0..r {
            let w = s.random_w();
            s.tab.push(w);
        }
        s
    }
    fn random_w(&mut self) -> W {
        let n = self.c.n;
        let a = self.rng.below(n);
        let b = self.rng.below(n);
        let p = self.c.add(self.c.mul(a, self.c.g), self.c.mul(b, self.q));
        W { p, a, b }
    }
    fn wadd(&self, w: W, t: W, minus: bool) -> W {
        let n = self.c.n;
        if minus {
            W { p: self.c.add(w.p, self.c.neg(t.p)), a: sadd(w.a, sneg(t.a, n), n), b: sadd(w.b, sneg(t.b, n), n) }
        } else {
            W { p: self.c.add(w.p, t.p), a: sadd(w.a, t.a, n), b: sadd(w.b, t.b, n) }
        }
    }
    #[inline]
    fn idx(&self, x: u64, key: usize, r: usize) -> usize {
        ((mix(x ^ self.keys[key]) as u128 * r as u128) >> 64) as usize
    }
    #[inline]
    fn score(&self, p: &Pt) -> u64 {
        if p.inf { 0 } else { mix(p.x ^ self.keys[4]) }
    }
    fn check(&self, a1: u64, b1: u64, a2: u64, b2: u64) -> bool {
        solve(&self.c, a1, b1, a2, b2) == Some(self.k)
    }
}

fn key(p: &Pt) -> (u64, u64) {
    if p.inf { (u64::MAX, u64::MAX) } else { (p.x, p.y) }
}

struct Rec {
    steps: u64,
    adds: u64,
    solved: bool,
    censored: bool,
    escapes: u64,
    restarts: u64,
}

// -------------------------------------------------------------------- e1

fn e1_step(s: &Inst, arm: &str, w: W) -> (W, u64) {
    let r = s.tab.len();
    let i = s.idx(w.p.x, 0, r);
    let best = |cands: &[W]| -> W {
        *cands.iter().max_by_key(|c| s.score(&c.p)).unwrap()
    };
    match arm {
        "base" => (s.wadd(w, s.tab[i], false), 1),
        "pm_coin" => (s.wadd(w, s.tab[i], mix(w.p.x ^ s.keys[3]) & 1 == 1), 1),
        "pm_score" => (best(&[s.wadd(w, s.tab[i], false), s.wadd(w, s.tab[i], true)]), 1),
        "two_score" => {
            let j = s.idx(w.p.x, 1, r);
            (best(&[s.wadd(w, s.tab[i], false), s.wadd(w, s.tab[j], false)]), 2)
        }
        "best4_score" => {
            let j = s.idx(w.p.x, 1, r);
            (
                best(&[
                    s.wadd(w, s.tab[i], false),
                    s.wadd(w, s.tab[i], true),
                    s.wadd(w, s.tab[j], false),
                    s.wadd(w, s.tab[j], true),
                ]),
                2,
            )
        }
        _ => panic!("unknown e1 arm {arm}"),
    }
}

/// Single walk, every point stored; the first repeat with different
/// coefficients ends the trial. A repeat with the same coefficients is a
/// fruitless cycle: with `escape` the walk leaves it by doubling (counted),
/// without it the trial ends unsolved there.
fn e1_trial(c: &Curve, arm: &str, seed: u64, r: usize, cap: u64, mode: &str) -> Rec {
    let n = c.n;
    let escape = mode != "off";
    let mut s = Inst::new(c, seed, r);
    let shift = s.random_w(); // fixed escape offset for mode "shift"
    let mut w = s.random_w();
    let mut seen: FastMap<(u64, u64), (u64, u64)> = FastMap::default();
    let (mut steps, mut adds, mut escapes) = (0u64, 0u64, 0u64);
    loop {
        if let Some(&(a, b)) = seen.get(&key(&w.p)) {
            if (a, b) != (w.a, w.b) || !escape {
                let solved = s.check(a, b, w.a, w.b);
                return Rec { steps, adds, solved, censored: false, escapes, restarts: 0 };
            }
            w = if mode == "shift" {
                s.wadd(w, shift, false)
            } else {
                W { p: c.dbl(w.p), a: sadd(w.a, w.a, n), b: sadd(w.b, w.b, n) }
            };
            escapes += 1;
            adds += 1;
            continue;
        }
        seen.insert(key(&w.p), (w.a, w.b));
        let (nw, k) = e1_step(&s, arm, w);
        w = nw;
        steps += 1;
        adds += k;
        if steps >= cap {
            return Rec { steps, adds, solved: false, censored: true, escapes, restarts: 0 };
        }
    }
}

/// Abstract random-function control: f(x) = best-scoring of k hashed images.
fn rf_trial(n: u64, k: usize, seed: u64, cap: u64) -> Rec {
    let mut rng = Rng(mix(seed) ^ 0xAB57);
    let keys: Vec<u64> = (0..k + 1).map(|_| rng.next()).collect();
    let mut x = rng.below(n);
    let mut seen: FastMap<u64, ()> = FastMap::default();
    let mut steps = 0u64;
    loop {
        if seen.insert(x, ()).is_some() {
            return Rec { steps, adds: steps, solved: true, censored: false, escapes: 0, restarts: 0 };
        }
        x = (0..k)
            .map(|c| ((mix(x ^ keys[c]) as u128 * n as u128) >> 64) as u64)
            .max_by_key(|y| mix(y ^ keys[k]))
            .unwrap();
        steps += 1;
        if steps >= cap {
            return Rec { steps, adds: steps, solved: false, censored: true, escapes: 0, restarts: 0 };
        }
    }
}

// --------------------------------------------- e2 / e4: distinguished points

struct DpCfg {
    walks: usize,
    dpbits: u32,
    r: usize,
    cap: u64,
}

/// Representative choice for the negation map: true means take -P.
fn want_neg(c: &Curve, p: &Pt, flip: bool) -> bool {
    !p.inf && ((p.y > c.f.neg(p.y)) ^ flip)
}

fn dp_trial(c: &Curve, arm: &str, seed: u64, cfg: &DpCfg) -> Rec {
    let mut s = Inst::new(c, seed, cfg.r);
    let n = c.n;
    let neg = arm.starts_with("neg_");
    let sel = arm.strip_prefix("sel_"); // e1 selection walks under the DP method
    let mask = (1u64 << cfg.dpbits) - 1;
    let restart_after = 20u64 << cfg.dpbits;
    let epoch_len = 4u64 << cfg.dpbits;
    let p = c.f.p;
    // pair table R_i + R_j (r^2 entries) for the "pair_table" arm
    let pairs: Vec<W> = if arm == "pair_table" {
        let mut v = Vec::with_capacity(cfg.r * cfg.r);
        for i in 0..cfg.r {
            for j in 0..cfg.r {
                v.push(s.wadd(s.tab[i], s.tab[j], false));
            }
        }
        v
    } else {
        Vec::new()
    };
    let canon = |w: W, flip: bool| -> W {
        if want_neg(c, &w.p, flip) {
            W { p: c.neg(w.p), a: sneg(w.a, n), b: sneg(w.b, n) }
        } else {
            w
        }
    };
    let flip_of = |x: u64, round: u64| -> bool {
        match arm {
            "neg_hashrule" => mix(x ^ 0xF11F) & 1 == 1,
            "neg_epoch" => (round / epoch_len) & 1 == 1,
            _ => false,
        }
    };
    struct St {
        w: W,
        prev: Pt,
        c: u64,
        since: u64,
    }
    let fresh = |s: &mut Inst, round: u64| -> St {
        let mut w = s.random_w();
        if neg {
            w = canon(w, flip_of(w.p.x, round));
        }
        St { w, prev: INF, c: 0, since: 0 }
    };
    let mut walks: Vec<St> = (0..cfg.walks).map(|_| fresh(&mut s, 0)).collect();
    let mut dps: FastMap<(u64, u64), (u64, u64)> = FastMap::default();
    let (mut steps, mut adds, mut escapes, mut restarts) = (0u64, 0u64, 0u64, 0u64);
    let mut round = 0u64;
    loop {
        for wi in 0..walks.len() {
            let st = &mut walks[wi];
            let x = st.w.p.x;
            let rep = x + st.c * p; // lazy representative of x: x or x + p
            let i = match arm {
                "lazy_index" => s.idx(rep, 0, cfg.r),
                "shared_index" => ((mix(round ^ s.keys[5]) as u128 * cfg.r as u128) >> 64) as usize,
                _ => s.idx(x, 0, cfg.r),
            };
            let mut nw = if let Some(sa) = sel {
                let (nw, k) = e1_step(&s, sa, st.w);
                adds += k - 1;
                nw
            } else if arm == "pair_table" {
                let j = s.idx(x, 1, cfg.r);
                s.wadd(st.w, pairs[i * cfg.r + j], false)
            } else {
                s.wadd(st.w, s.tab[i], false)
            };
            adds += 1;
            // the redundant form of the new x depends only on what produced it
            let nc = mix(x ^ (i as u64) ^ s.keys[6]) & 1;
            if neg {
                nw = canon(nw, flip_of(nw.p.x, round));
            }
            if neg || sel.is_some() {
                if !st.prev.inf && nw.p == st.prev {
                    // fruitless 2-cycle {P, f(P)}: leave by doubling the smaller x
                    let m = if st.w.p.x < nw.p.x { st.w } else { nw };
                    nw = W { p: c.dbl(m.p), a: sadd(m.a, m.a, n), b: sadd(m.b, m.b, n) };
                    if neg {
                        nw = canon(nw, flip_of(nw.p.x, round));
                    }
                    adds += 1;
                    escapes += 1;
                }
            }
            st.prev = st.w.p;
            st.w = nw;
            st.c = nc;
            st.since += 1;
            steps += 1;
            let dpv = if arm == "lazy_dp" { st.w.p.x + st.c * p } else { st.w.p.x };
            if mix(dpv ^ s.keys[7]) & mask == 0 {
                st.since = 0;
                let k = key(&st.w.p);
                match dps.get(&k) {
                    Some(&(a, b)) if (a, b) != (st.w.a, st.w.b) => {
                        let solved = s.check(a, b, st.w.a, st.w.b);
                        return Rec { steps, adds, solved, censored: false, escapes, restarts };
                    }
                    Some(_) => {}
                    None => {
                        dps.insert(k, (st.w.a, st.w.b));
                    }
                }
            }
            if st.since > restart_after {
                *st = fresh(&mut s, round);
                restarts += 1;
            }
            if steps >= cfg.cap {
                return Rec { steps, adds, solved: false, censored: true, escapes, restarts };
            }
        }
        round += 1;
    }
}

// ----------------------------------------------------------------- bench

/// Branch-free L-limb Montgomery field for timing only: every operation runs
/// the same instructions whatever the data, so short cycles in a walk cannot
/// make it look cheaper. Values need not be curve points; the cost does not
/// depend on them.
#[derive(Clone, Copy)]
struct MF<const L: usize> {
    p: [u64; L],
    ninv: u64,
    one: [u64; L],
}
impl<const L: usize> MF<L> {
    fn new(p: [u64; L]) -> Self {
        let mut inv = p[0];
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(p[0].wrapping_mul(inv)));
        }
        let mut f = MF { p, ninv: inv.wrapping_neg(), one: [0; L] };
        let mut r = [0u64; L];
        r[0] = 1;
        for _ in 0..64 * L {
            r = f.add(&r, &r);
        }
        f.one = r;
        f
    }
    #[inline(always)]
    fn sub_p_masked(&self, t: &[u64; L], carry: u64) -> [u64; L] {
        // t - p if t + carry*2^(64L) >= p, else t; branch-free
        let mut d = [0u64; L];
        let mut br = 0u64;
        for i in 0..L {
            let (x, b1) = t[i].overflowing_sub(self.p[i]);
            let (y, b2) = x.overflowing_sub(br);
            d[i] = y;
            br = (b1 | b2) as u64;
        }
        let keep_t = (br & !carry & 1).wrapping_neg(); // all ones: keep t
        let mut o = [0u64; L];
        for i in 0..L {
            o[i] = (t[i] & keep_t) | (d[i] & !keep_t);
        }
        o
    }
    #[inline(always)]
    fn add(&self, a: &[u64; L], b: &[u64; L]) -> [u64; L] {
        let mut t = [0u64; L];
        let mut c = 0u64;
        for i in 0..L {
            let s = a[i] as u128 + b[i] as u128 + c as u128;
            t[i] = s as u64;
            c = (s >> 64) as u64;
        }
        self.sub_p_masked(&t, c)
    }
    #[inline(always)]
    fn sub(&self, a: &[u64; L], b: &[u64; L]) -> [u64; L] {
        let mut t = [0u64; L];
        let mut br = 0u64;
        for i in 0..L {
            let (x, b1) = a[i].overflowing_sub(b[i]);
            let (y, b2) = x.overflowing_sub(br);
            t[i] = y;
            br = (b1 | b2) as u64;
        }
        let m = br.wrapping_neg();
        let mut c = 0u64;
        for i in 0..L {
            let s = t[i] as u128 + (self.p[i] & m) as u128 + c as u128;
            t[i] = s as u64;
            c = (s >> 64) as u64;
        }
        t
    }
    #[inline(always)]
    fn neg(&self, a: &[u64; L]) -> [u64; L] {
        self.sub(&[0; L], a)
    }
    /// CIOS Montgomery multiplication.
    #[inline(always)]
    fn mul(&self, a: &[u64; L], b: &[u64; L]) -> [u64; L] {
        let mut t = [0u64; L];
        let mut t_hi = 0u64;
        for i in 0..L {
            let mut c = 0u64;
            for j in 0..L {
                let s = t[j] as u128 + a[j] as u128 * b[i] as u128 + c as u128;
                t[j] = s as u64;
                c = (s >> 64) as u64;
            }
            let s = t_hi as u128 + c as u128;
            let (tl, th) = (s as u64, (s >> 64) as u64);
            let m = t[0].wrapping_mul(self.ninv);
            let mut c = ((t[0] as u128 + m as u128 * self.p[0] as u128) >> 64) as u64;
            for j in 1..L {
                let s = t[j] as u128 + m as u128 * self.p[j] as u128 + c as u128;
                t[j - 1] = s as u64;
                c = (s >> 64) as u64;
            }
            let s = tl as u128 + c as u128;
            t[L - 1] = s as u64;
            t_hi = th + (s >> 64) as u64;
        }
        self.sub_p_masked(&t, t_hi)
    }
    fn inv(&self, a: &[u64; L]) -> [u64; L] {
        // a^(p-2), fixed exponent: square-and-multiply on the bits of p-2
        let mut e = self.p;
        e[0] -= 2; // p is odd and > 2, no borrow for the primes used here
        let mut r = self.one;
        for i in (0..L).rev() {
            for k in (0..64).rev() {
                r = self.mul(&r, &r);
                if (e[i] >> k) & 1 == 1 {
                    r = self.mul(&r, a);
                }
            }
        }
        r
    }
    fn random(&self, rng: &mut Rng) -> [u64; L] {
        let mut v = [0u64; L];
        for x in v.iter_mut() {
            *x = rng.next();
        }
        v[L - 1] %= self.p[L - 1].max(1);
        v
    }
}

#[inline(always)]
fn sel<const L: usize>(c: bool, a: &[u64; L], b: &[u64; L]) -> [u64; L] {
    let m = (c as u64).wrapping_neg();
    let mut o = [0u64; L];
    for i in 0..L {
        o[i] = (a[i] & m) | (b[i] & !m);
    }
    o
}

/// Batched-affine step cost over m walks (Montgomery simultaneous inversion),
/// branch-free arithmetic and selection.
fn bench_arm<const L: usize>(f: &MF<L>, arm: &str, m: usize, steps: usize, seed: u64) -> (f64, u64) {
    let mut rng = Rng(seed);
    let r = 32;
    let tab: Vec<([u64; L], [u64; L])> = (0..r).map(|_| (f.random(&mut rng), f.random(&mut rng))).collect();
    let mut xs: Vec<[u64; L]> = (0..m).map(|_| f.random(&mut rng)).collect();
    let mut ys: Vec<[u64; L]> = (0..m).map(|_| f.random(&mut rng)).collect();
    let k2 = arm == "best4";
    let nd = if k2 { 2 * m } else { m };
    let mut den = vec![[0u64; L]; nd];
    let mut pre = vec![[0u64; L]; nd];
    let mut ii = vec![0usize; nd];
    let ks = 0x5C0E_u64;
    let t0 = Instant::now();
    for _ in 0..steps {
        for w in 0..m {
            let i = ((mix(xs[w][0] ^ 0x1D) as u128 * r as u128) >> 64) as usize;
            ii[w] = i;
            den[w] = f.sub(&tab[i].0, &xs[w]);
            if k2 {
                let j = ((mix(xs[w][0] ^ 0x2D) as u128 * r as u128) >> 64) as usize;
                ii[m + w] = j;
                den[m + w] = f.sub(&tab[j].0, &xs[w]);
            }
        }
        // simultaneous inversion: 3(nd-1) multiplications + one inversion
        let mut acc = f.one;
        for d in 0..nd {
            pre[d] = acc;
            acc = f.mul(&acc, &den[d]);
        }
        let mut inv = f.inv(&acc);
        for d in (0..nd).rev() {
            let id = f.mul(&inv, &pre[d]);
            inv = f.mul(&inv, &den[d]);
            den[d] = id;
        }
        for w in 0..m {
            let (xp, yp) = (xs[w], ys[w]);
            let (tx, ty) = tab[ii[w]];
            let (l, x3) = match arm {
                "base" => {
                    let l = f.mul(&f.sub(&ty, &yp), &den[w]);
                    (l, f.sub(&f.sub(&f.mul(&l, &l), &xp), &tx))
                }
                "pm_score" => {
                    let lp = f.mul(&f.sub(&ty, &yp), &den[w]);
                    let lm = f.mul(&f.neg(&f.add(&ty, &yp)), &den[w]);
                    let xp3 = f.sub(&f.sub(&f.mul(&lp, &lp), &xp), &tx);
                    let xm3 = f.sub(&f.sub(&f.mul(&lm, &lm), &xp), &tx);
                    let c = mix(xp3[0] ^ ks) >= mix(xm3[0] ^ ks);
                    (sel(c, &lp, &lm), sel(c, &xp3, &xm3))
                }
                "best4" => {
                    let (ux, uy) = tab[ii[m + w]];
                    let mut bl = [0u64; L];
                    let mut bx = [0u64; L];
                    let mut bs = 0u64;
                    for (cx, cy, inv) in [(tx, ty, den[w]), (ux, uy, den[m + w])] {
                        for sgn in [false, true] {
                            let yy = sel(sgn, &f.neg(&cy), &cy);
                            let l = f.mul(&f.sub(&yy, &yp), &inv);
                            let x3 = f.sub(&f.sub(&f.mul(&l, &l), &xp), &cx);
                            let sc = mix(x3[0] ^ ks);
                            let c = sc >= bs;
                            bl = sel(c, &l, &bl);
                            bx = sel(c, &x3, &bx);
                            bs = if c { sc } else { bs };
                        }
                    }
                    (bl, bx)
                }
                _ => panic!("unknown bench arm {arm}"),
            };
            ys[w] = f.sub(&f.mul(&l, &f.sub(&xp, &x3)), &yp);
            xs[w] = x3;
        }
    }
    let ns = t0.elapsed().as_nanos() as f64 / (steps * m) as f64;
    let chk = xs.iter().fold(0u64, |a, x| a ^ x[0]);
    (ns, chk)
}

fn bench_run<const L: usize>(f: &MF<L>, label: &str, arms: &[String], m: usize, steps: usize, reps: usize, seed: u64) {
    bench_arm(f, "base", m, steps / 4, seed); // warm up
    for rep in 0..reps {
        // rotate the arm order every repetition
        for k in 0..arms.len() {
            let arm = &arms[(k + rep) % arms.len()];
            let (ns, chk) = bench_arm(f, arm, m, steps, seed + rep as u64);
            println!(
                "{{\"exp\":\"bench\",\"bits\":{},\"field\":\"{label}\",\"arm\":\"{arm}\",\"rep\":{rep},\"walks\":{m},\"steps\":{steps},\"ns_per_step\":{ns:.3},\"chk\":{chk}}}",
                64 * L
            );
        }
    }
}

// ----------------------------------------------------------------- driver

fn arg<T: std::str::FromStr>(args: &[String], name: &str, def: T) -> T {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(def)
}
fn arg_s(args: &[String], name: &str, def: &str) -> String {
    args.iter().position(|a| a == name).and_then(|i| args.get(i + 1)).cloned().unwrap_or_else(|| def.into())
}

fn threads() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
}

fn run_trials(exp: &str, bits: u32, arms: &[String], t0: usize, trials: usize, seed0: u64, f: &(dyn Fn(&str, u64) -> Rec + Sync), extra: &str) {
    let jobs: Vec<(usize, &str)> = (t0..t0 + trials).flat_map(|t| arms.iter().map(move |a| (t, a.as_str()))).collect();
    let next = AtomicUsize::new(0);
    let out = Mutex::new(std::io::stdout());
    std::thread::scope(|sc| {
        for _ in 0..threads() {
            sc.spawn(|| loop {
                let j = next.fetch_add(1, Ordering::Relaxed);
                if j >= jobs.len() {
                    break;
                }
                let (t, arm) = jobs[j];
                let seed = seed0.wrapping_add(t as u64);
                let t0 = Instant::now();
                let r = f(arm, seed);
                let line = format!(
                    "{{\"exp\":\"{exp}\",\"bits\":{bits},\"arm\":\"{arm}\",\"trial\":{t},\"seed\":{seed},\"steps\":{},\"adds\":{},\"solved\":{},\"censored\":{},\"escapes\":{},\"restarts\":{},\"wall_ms\":{:.1}{extra}}}",
                    r.steps, r.adds, r.solved, r.censored, r.escapes, r.restarts,
                    t0.elapsed().as_secs_f64() * 1e3
                );
                let mut o = out.lock().unwrap();
                writeln!(o, "{line}").unwrap();
            });
        }
    });
}

fn curve_line(c: &Curve, seed: u64) {
    println!(
        "{{\"exp\":\"curve\",\"bits\":{},\"seed\":{seed},\"p\":{},\"n\":{},\"a_mont\":{},\"b_mont\":{},\"gx_mont\":{},\"gy_mont\":{},\"threads\":{}}}",
        c.bits, c.f.p, c.n, c.a, c.b, c.g.x, c.g.y, threads()
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cmd = args.get(1).map(|s| s.as_str()).unwrap_or("help");
    let bits: u32 = arg(&args, "--bits", 28);
    let trials: usize = arg(&args, "--trials", 100);
    // first trial index: shards of one run share the curve (--seed) and split trials
    let t0: usize = arg(&args, "--t0", 0);
    let seed: u64 = arg(&args, "--seed", 20261014);
    let arms_s = arg_s(&args, "--arms", "");
    let arms: Vec<String> = arms_s.split(',').filter(|s| !s.is_empty()).map(String::from).collect();
    match cmd {
        "e1" => {
            let c = Curve::generate(bits, seed);
            curve_line(&c, seed);
            let r: usize = arg(&args, "--r", 1024);
            let mode = arg_s(&args, "--escape", "double");
            let exp = match mode.as_str() { "double" => "e1", "off" => "e1noesc", "shift" => "e1shift", _ => panic!("--escape double|shift|off") };
            let cap = (40.0 * (c.n as f64).sqrt()) as u64;
            let arms = if arms.is_empty() {
                ["base", "pm_coin", "pm_score", "two_score", "best4_score"].map(String::from).to_vec()
            } else {
                arms
            };
            run_trials(exp, bits, &arms, t0, trials, seed, &|arm, s| e1_trial(&c, arm, s, r, cap, &mode), &format!(",\"r\":{r},\"escape\":\"{mode}\""));
        }
        "rf" => {
            let c = Curve::generate(bits, seed);
            let n = c.n; // same group size as the curve experiments
            let cap = (40.0 * (n as f64).sqrt()) as u64;
            let arms = if arms.is_empty() { ["k1", "k2", "k4"].map(String::from).to_vec() } else { arms };
            run_trials("rf", bits, &arms, t0, trials, seed, &|arm, s| rf_trial(n, arm[1..].parse().unwrap(), s, cap), "");
        }
        "e2" | "e4" | "e1dp" => {
            let c = Curve::generate(bits, seed);
            curve_line(&c, seed);
            let cfg = DpCfg {
                walks: arg(&args, "--walks", 64),
                dpbits: arg(&args, "--dpbits", 8),
                r: arg(&args, "--r", if cmd == "e2" { 32 } else { 1024 }),
                cap: (arg::<f64>(&args, "--cap", 30.0) * (c.n as f64).sqrt()) as u64,
            };
            let arms = if !arms.is_empty() {
                arms
            } else if cmd == "e1dp" {
                ["sel_base", "sel_pm_coin", "sel_pm_score", "sel_two_score", "sel_best4_score"].map(String::from).to_vec()
            } else if cmd == "e2" {
                ["canon", "lazy_dp", "lazy_index", "shared_index", "pair_table"].map(String::from).to_vec()
            } else {
                ["nonneg", "neg_canon", "neg_hashrule", "neg_epoch"].map(String::from).to_vec()
            };
            let extra = format!(",\"walks\":{},\"dpbits\":{},\"r\":{}", cfg.walks, cfg.dpbits, cfg.r);
            run_trials(cmd, bits, &arms, t0, trials, seed, &|arm, s| dp_trial(&c, arm, s, &cfg), &extra);
        }
        "bench" => {
            let m: usize = arg(&args, "--walks", 512);
            let steps: usize = arg(&args, "--steps", 2000);
            let reps: usize = arg(&args, "--reps", 11);
            let arms = if arms.is_empty() { ["base", "pm_score", "best4"].map(String::from).to_vec() } else { arms };
            // 1 limb: the 62-bit prime 2^62 - 57; 4 limbs: the P-256 prime
            let f1 = MF::<1>::new([(1u64 << 62) - 57]);
            bench_run(&f1, "2^62-57", &arms, m, steps, reps, seed);
            let f4 = MF::<4>::new([0xFFFF_FFFF_FFFF_FFFF, 0x0000_0000_FFFF_FFFF, 0, 0xFFFF_FFFF_0000_0001]);
            bench_run(&f4, "P-256", &arms, m, steps / 4, reps, seed);
        }
        "selftest" => {
            fn t<const L: usize>(f: &MF<L>, name: &str) {
                let mut rng = Rng(7);
                for _ in 0..2000 {
                    let (a, b, c) = (f.random(&mut rng), f.random(&mut rng), f.random(&mut rng));
                    assert!(f.mul(&f.inv(&a), &a) == f.one, "{name} inverse");
                    assert!(f.sub(&f.add(&a, &b), &b) == a, "{name} add/sub");
                    assert!(f.mul(&a, &f.add(&b, &c)) == f.add(&f.mul(&a, &b), &f.mul(&a, &c)), "{name} distributive");
                    assert!(f.add(&a, &f.neg(&a)) == [0; L], "{name} neg");
                }
                println!("{{\"exp\":\"selftest\",\"field\":\"{name}\",\"ok\":true}}");
            }
            assert!(is_prime((1u64 << 62) - 57));
            t(&MF::<1>::new([(1u64 << 62) - 57]), "2^62-57");
            t(&MF::<4>::new([0xFFFF_FFFF_FFFF_FFFF, 0x0000_0000_FFFF_FFFF, 0, 0xFFFF_FFFF_0000_0001]), "P-256");
        }
        "summary" => summary(&args),
        _ => eprintln!("usage: rho_walk_leads e1|rf|e2|e4|bench|summary [--bits B] [--trials T] [--seed S] [--arms a,b]"),
    }
}

// --------------------------------------------------------------- summary

fn jget<'a>(l: &'a str, k: &str) -> Option<&'a str> {
    let pat = format!("\"{k}\":");
    let i = l.find(&pat)? + pat.len();
    let rest = &l[i..];
    let rest = rest.strip_prefix('"').unwrap_or(rest);
    let end = rest.find(|ch| ch == ',' || ch == '}' || ch == '"').unwrap_or(rest.len());
    Some(&rest[..end])
}

/// Predicted steps ratio reference/arm (a speedup when above 1).
fn predicted(exp: &str, arm: &str) -> Option<f64> {
    let (exp, arm) = if exp.starts_with("e1") { ("e1", arm.strip_prefix("sel_").unwrap_or(arm)) } else { (exp, arm) };
    let k2 = |k: f64| (k * k / (2.0 * k - 1.0)).sqrt();
    match (exp, arm) {
        ("e1", "pm_score") | ("e1", "two_score") | ("rf", "k2") => Some(k2(2.0)),
        ("e1", "best4_score") | ("rf", "k4") => Some(k2(4.0)),
        ("e1", "pm_coin") | ("e2", "lazy_dp") | ("e2", "pair_table") => Some(1.0),
        ("e2", "lazy_index") => Some(std::f64::consts::FRAC_1_SQRT_2),
        ("e4", "neg_canon") | ("e4", "neg_hashrule") => Some(std::f64::consts::SQRT_2),
        _ => None,
    }
}

fn summary(args: &[String]) {
    let path = arg_s(args, "--in", "-");
    let rd: Box<dyn BufRead> = if path == "-" {
        Box::new(std::io::BufReader::new(std::io::stdin()))
    } else {
        Box::new(std::io::BufReader::new(std::fs::File::open(&path).expect("open")))
    };
    // (exp, bits, arm) -> trial -> (steps, solved, censored, escapes, restarts, adds)
    let mut cells: std::collections::BTreeMap<(String, u32, String), std::collections::BTreeMap<u64, [f64; 6]>> = Default::default();
    let mut bench: std::collections::BTreeMap<(u32, String), Vec<f64>> = Default::default();
    for l in rd.lines() {
        let l = l.unwrap();
        let Some(exp) = jget(&l, "exp") else { continue };
        let bits: u32 = jget(&l, "bits").and_then(|v| v.parse().ok()).unwrap_or(0);
        if exp == "bench" {
            let arm = jget(&l, "arm").unwrap().to_string();
            bench.entry((bits, arm)).or_default().push(jget(&l, "ns_per_step").unwrap().parse().unwrap());
            continue;
        }
        if exp == "curve" {
            continue;
        }
        let arm = jget(&l, "arm").unwrap().to_string();
        let trial: u64 = jget(&l, "trial").unwrap().parse().unwrap();
        let g = |k: &str| -> f64 {
            match jget(&l, k).unwrap_or("0") {
                "true" => 1.0,
                "false" => 0.0,
                v => v.parse().unwrap_or(0.0),
            }
        };
        cells.entry((exp.into(), bits, arm)).or_default().insert(
            trial,
            [g("steps"), g("solved"), g("censored"), g("escapes"), g("restarts"), g("adds")],
        );
    }
    let reference = |exp: &str| match exp {
        "e1dp" => "sel_base",
        e if e.starts_with("e1") => "base",
        "rf" => "k1",
        "e2" => "canon",
        "e4" => "nonneg",
        _ => "",
    };
    for ((exp, bits, arm), rows) in &cells {
        let n = rows.len() as f64;
        let steps: Vec<f64> = rows.values().map(|r| r[0]).collect();
        let mean = steps.iter().sum::<f64>() / n;
        let sd = (steps.iter().map(|s| (s - mean).powi(2)).sum::<f64>() / (n - 1.0).max(1.0)).sqrt();
        let solved = rows.values().filter(|r| r[1] == 1.0).count();
        let cens = rows.values().filter(|r| r[2] == 1.0).count();
        let esc: f64 = rows.values().map(|r| r[3]).sum();
        let rst: f64 = rows.values().map(|r| r[4]).sum();
        let adds: f64 = rows.values().map(|r| r[5]).sum();
        let tot: f64 = steps.iter().sum();
        let mut line = format!(
            "{{\"exp\":\"summary\",\"of\":\"{exp}\",\"bits\":{bits},\"arm\":\"{arm}\",\"n\":{},\"solved\":{solved},\"censored\":{cens},\"mean_steps\":{mean:.1},\"sd_steps\":{sd:.1},\"se_steps\":{:.1},\"adds_per_step\":{:.4},\"escapes_per_1e6_steps\":{:.2},\"restarts\":{rst}",
            rows.len(),
            sd / n.sqrt(),
            adds / tot,
            esc / tot * 1e6
        );
        let rname = reference(exp);
        if let Some(rrows) = cells.get(&(exp.clone(), *bits, rname.to_string())) {
            if rname != arm {
                // paired by trial: ratio of means ref/arm, bootstrap over trials
                let common: Vec<(f64, f64)> =
                    rows.iter().filter_map(|(t, r)| rrows.get(t).map(|q| (q[0], r[0]))).collect();
                let m = common.len();
                let ratio = |idx: &mut dyn Iterator<Item = usize>| {
                    let (mut a, mut b) = (0.0, 0.0);
                    for i in idx {
                        a += common[i].0;
                        b += common[i].1;
                    }
                    a / b
                };
                let point = ratio(&mut (0..m));
                let mut rng = Rng(0xB007);
                let mut boots: Vec<f64> = (0..4000)
                    .map(|_| {
                        let v: Vec<usize> = (0..m).map(|_| rng.below(m as u64) as usize).collect();
                        ratio(&mut v.into_iter())
                    })
                    .collect();
                boots.sort_by(|a, b| a.partial_cmp(b).unwrap());
                line += &format!(
                    ",\"ref\":\"{rname}\",\"paired_n\":{m},\"speedup_ref_over_arm\":{point:.4},\"ci95\":[{:.4},{:.4}]",
                    boots[(0.025 * boots.len() as f64) as usize],
                    boots[(0.975 * boots.len() as f64) as usize]
                );
                if let Some(p) = predicted(exp, arm) {
                    line += &format!(",\"predicted\":{p:.4}");
                }
            }
        }
        if *exp == "e1" || *exp == "e1shift" || *exp == "rf" {
            let th = (std::f64::consts::PI * 2f64.powi(*bits as i32) / 2.0).sqrt();
            line += &format!(",\"sqrt_pi_n_over_2_at_2^bits\":{th:.1}");
        }
        println!("{line}}}");
    }
    for ((bits, arm), v) in &bench {
        let mut v = v.clone();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med = v[v.len() / 2];
        let mut line = format!(
            "{{\"exp\":\"summary\",\"of\":\"bench\",\"bits\":{bits},\"arm\":\"{arm}\",\"reps\":{},\"median_ns_per_step\":{med:.3},\"min\":{:.3},\"max\":{:.3}",
            v.len(),
            v[0],
            v[v.len() - 1]
        );
        if let Some(b) = bench.get(&(*bits, "base".to_string())) {
            let mut b = b.clone();
            b.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let bm = b[b.len() / 2];
            line += &format!(",\"cost_ratio_vs_base\":{:.4}", med / bm);
            // net work ratio: steps saved per unit cost, using the e1 speedups at the largest size
            let e1arm = match arm.as_str() {
                "pm_score" => Some("pm_score"),
                "best4" => Some("best4_score"),
                _ => None,
            };
            // predicted net: the derived steps speedup over the measured cost ratio
            if let Some(ea) = e1arm {
                let k = if ea == "pm_score" { 2.0f64 } else { 4.0 };
                let pred = (k * k / (2.0 * k - 1.0)).sqrt();
                line += &format!(",\"predicted_net_ratio\":{:.4}", pred / (med / bm));
            }
        }
        println!("{line}}}");
    }
}
