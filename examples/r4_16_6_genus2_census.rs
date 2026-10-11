//! **Lead 16.6: the affine saturation law on a genus-2 Jacobian.**
//!
//! Ports the exact census of `RESEARCH_ISOGENY_FIBERS_ECDLP.md` §15.2 (the
//! j = 0 elliptic census of `round3/automorphism_vs_translation.py`) to the
//! Jacobian of `C: y² = x⁵ + b` over `F_p`, `p ≡ 1 (mod 10)`. There
//! `ζ: (x, y) ↦ (ζ₅x, y)` and the hyperelliptic involution generate a rational
//! automorphism group `A` of order 10. Divisor classes are Mumford pairs
//! `(u, v)`, added by Cantor's algorithm (own `u64` implementation below,
//! cross-checked against `crypto_lib::cryptanalysis::jv_cover::Hyp`).
//!
//! For each base `F` (a set of divisor classes, closed under negation):
//! * `D_m`: the exact number of distinct nonzero classes of the order-`r`
//!   subgroup `G` that are sums of `m` base elements (with repetition);
//! * `R_m`: the number of such multisets (representations);
//! * `U`: columns, the orbits of `[h]F` under `A`;
//! * `c = D_m · m! / |F|^m`, and the mitm cost model `U · (r / D_m) · L_m`
//!   with `L_2 = 1`, `L_3 = |F|`, exactly as in §15.2.
//!
//! ```text
//! cargo run --release --example r4_16_6_genus2_census -- search --lo 4000 --hi 4400
//! cargo run --release --example r4_16_6_genus2_census -- census --p 4111 --b 3 \
//!     --size 80 --reps 30 --seed 1 --json out.json
//! ```

use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;

use crypto_lib::cryptanalysis::jv_cover::{Div as LDiv, Fld, Fq, Hyp, E2};
use crypto_lib::cryptanalysis::residual_walk::is_prime_u64;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

// ── F_p ────────────────────────────────────────────────────────────────────

#[derive(Clone)]
struct Fp {
    p: u64,
    /// Legendre symbol table: 0, 1 or p − 1 (as −1).
    chi: Vec<i8>,
    /// A fixed non-residue (for F_{p²} = F_p[t]/(t² − w)).
    w: u64,
}

impl Fp {
    fn new(p: u64) -> Fp {
        let mut chi = vec![-1i8; p as usize];
        chi[0] = 0;
        for x in 1..p {
            chi[(x * x % p) as usize] = 1;
        }
        let w = (2..p).find(|&w| chi[w as usize] == -1).unwrap();
        Fp { p, chi, w }
    }
    fn add(&self, a: u64, b: u64) -> u64 {
        (a + b) % self.p
    }
    fn sub(&self, a: u64, b: u64) -> u64 {
        (a + self.p - b) % self.p
    }
    fn neg(&self, a: u64) -> u64 {
        (self.p - a) % self.p
    }
    fn mul(&self, a: u64, b: u64) -> u64 {
        a * b % self.p
    }
    fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1;
        a %= self.p;
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
        assert!(!a.is_multiple_of(self.p), "inverse of zero");
        self.pow(a, self.p - 2)
    }
    fn is_qr(&self, a: u64) -> bool {
        self.chi[(a % self.p) as usize] >= 0
    }
    /// Tonelli–Shanks.
    fn sqrt(&self, a: u64) -> Option<u64> {
        let p = self.p;
        let a = a % p;
        if a == 0 {
            return Some(0);
        }
        if self.chi[a as usize] != 1 {
            return None;
        }
        let mut q = p - 1;
        let mut s = 0;
        while q.is_multiple_of(2) {
            q /= 2;
            s += 1;
        }
        let mut m = s;
        let mut c = self.pow(self.w, q);
        let mut t = self.pow(a, q);
        let mut r = self.pow(a, q.div_ceil(2));
        while t != 1 {
            let mut i = 0;
            let mut tt = t;
            while tt != 1 {
                tt = self.mul(tt, tt);
                i += 1;
            }
            let mut b = c;
            for _ in 0..(m - i - 1) {
                b = self.mul(b, b);
            }
            m = i;
            c = self.mul(b, b);
            t = self.mul(t, c);
            r = self.mul(r, b);
        }
        Some(r)
    }
    // F_{p²} = F_p[t]/(t² − w), elements (a0, a1).
    fn m2(&self, a: (u64, u64), b: (u64, u64)) -> (u64, u64) {
        let p = self.p;
        (
            (a.0 * b.0 % p + self.w * (a.1 * b.1 % p)) % p,
            (a.0 * b.1 + a.1 * b.0) % p,
        )
    }
    fn norm2(&self, a: (u64, u64)) -> u64 {
        self.sub(self.mul(a.0, a.0), self.mul(self.w, self.mul(a.1, a.1)))
    }
    /// Square root in F_{p²} by the norm method; verified before return.
    fn sqrt2(&self, a: (u64, u64)) -> Option<(u64, u64)> {
        if a.1 == 0 {
            if let Some(s) = self.sqrt(a.0) {
                return Some((s, 0));
            }
            let s = self.sqrt(self.mul(a.0, self.inv(self.w)))?;
            return Some((0, s));
        }
        let n = self.norm2(a);
        let s = self.sqrt(n)?;
        let half = self.inv(2);
        let mut x0sq = self.mul(self.add(a.0, s), half);
        if !self.is_qr(x0sq) {
            x0sq = self.mul(self.sub(a.0, s), half);
        }
        let x0 = self.sqrt(x0sq)?;
        if x0 == 0 {
            return None;
        }
        let x1 = self.mul(a.1, self.inv(self.mul(2, x0)));
        let r = (x0, x1);
        if self.m2(r, r) == a {
            Some(r)
        } else {
            None
        }
    }
}

// ── Polynomials over F_p (trimmed Vec, low degree first) ─────────────────

type Poly = Vec<u64>;

fn trim(mut a: Poly) -> Poly {
    while a.last() == Some(&0) {
        a.pop();
    }
    a
}
fn padd(k: &Fp, a: &[u64], b: &[u64]) -> Poly {
    let n = a.len().max(b.len());
    trim(
        (0..n)
            .map(|i| k.add(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0)))
            .collect(),
    )
}
fn psub(k: &Fp, a: &[u64], b: &[u64]) -> Poly {
    let n = a.len().max(b.len());
    trim(
        (0..n)
            .map(|i| k.sub(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0)))
            .collect(),
    )
}
fn pmul(k: &Fp, a: &[u64], b: &[u64]) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = (r[i + j] + x * y) % k.p;
        }
    }
    trim(r)
}
fn pscale(k: &Fp, a: &[u64], c: u64) -> Poly {
    trim(a.iter().map(|&x| k.mul(x, c)).collect())
}
fn pdivrem(k: &Fp, a: &[u64], b: &[u64]) -> (Poly, Poly) {
    assert!(!b.is_empty(), "division by zero polynomial");
    let mut r: Poly = trim(a.to_vec());
    if r.len() < b.len() {
        return (Vec::new(), r);
    }
    let db = b.len() - 1;
    let li = k.inv(b[db]);
    let mut q = vec![0u64; r.len() - db];
    while r.len() >= b.len() {
        let d = r.len() - 1;
        let c = k.mul(r[d], li);
        q[d - db] = c;
        for (j, &bj) in b.iter().enumerate() {
            r[d - db + j] = k.sub(r[d - db + j], k.mul(c, bj));
        }
        r = trim(r);
    }
    (trim(q), r)
}
fn pmonic(k: &Fp, a: &[u64]) -> Poly {
    match a.last() {
        None => Vec::new(),
        Some(&lc) => pscale(k, a, k.inv(lc)),
    }
}
/// `(d, s, t)` with `d = s·a + t·b` monic.
fn pxgcd(k: &Fp, a: &[u64], b: &[u64]) -> (Poly, Poly, Poly) {
    let (mut r0, mut r1) = (trim(a.to_vec()), trim(b.to_vec()));
    let (mut s0, mut s1): (Poly, Poly) = (vec![1], Vec::new());
    let (mut t0, mut t1): (Poly, Poly) = (Vec::new(), vec![1]);
    while !r1.is_empty() {
        let (q, r) = pdivrem(k, &r0, &r1);
        let s = psub(k, &s0, &pmul(k, &q, &s1));
        let t = psub(k, &t0, &pmul(k, &q, &t1));
        r0 = std::mem::replace(&mut r1, r);
        s0 = std::mem::replace(&mut s1, s);
        t0 = std::mem::replace(&mut t1, t);
    }
    match r0.last() {
        None => (r0, s0, t0),
        Some(&lc) => {
            let i = k.inv(lc);
            (pscale(k, &r0, i), pscale(k, &s0, i), pscale(k, &t0, i))
        }
    }
}
fn peval(k: &Fp, a: &[u64], x: u64) -> u64 {
    a.iter().rev().fold(0, |acc, &c| k.add(k.mul(acc, x), c))
}

// ── The Jacobian of y² = x⁵ + b ───────────────────────────────────────────

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Dv {
    u: Poly,
    v: Poly,
}

struct Jac {
    k: Fp,
    b: u64,
    f: Poly,
    /// A primitive 5th root of unity in F_p.
    z5: u64,
}

impl Jac {
    fn new(p: u64, b: u64) -> Jac {
        let k = Fp::new(p);
        let f = vec![b % p, 0, 0, 0, 0, 1];
        let z5 = (2..p)
            .map(|g| k.pow(g, (p - 1) / 5))
            .find(|&z| z != 1)
            .unwrap();
        Jac { k, b: b % p, f, z5 }
    }
    fn p(&self) -> u64 {
        self.k.p
    }
    fn zero(&self) -> Dv {
        Dv {
            u: vec![1],
            v: Vec::new(),
        }
    }
    fn is_zero(&self, d: &Dv) -> bool {
        d.u.len() == 1
    }
    fn fx(&self, x: u64) -> u64 {
        peval(&self.k, &self.f, x)
    }
    fn key(&self, d: &Dv) -> u64 {
        let g = |a: &Poly, i: usize| *a.get(i).unwrap_or(&0);
        let deg = (d.u.len() - 1) as u64;
        deg | g(&d.u, 0) << 2 | g(&d.u, 1) << 16 | g(&d.v, 0) << 30 | g(&d.v, 1) << 44
    }
    fn valid(&self, d: &Dv) -> bool {
        let k = &self.k;
        if d.u.last() != Some(&1) || d.u.len() > 3 || d.v.len() >= d.u.len() {
            return false;
        }
        pdivrem(k, &psub(k, &self.f, &pmul(k, &d.v, &d.v)), &d.u)
            .1
            .is_empty()
    }
    fn neg(&self, d: &Dv) -> Dv {
        Dv {
            u: d.u.clone(),
            v: trim(d.v.iter().map(|&c| self.k.neg(c)).collect()),
        }
    }
    /// The degree-1 divisor class `[(x, y) − ∞]`.
    fn pt(&self, x: u64, y: u64) -> Dv {
        let p = self.p();
        Dv {
            u: vec![(p - x % p) % p, 1],
            v: trim(vec![y % p]),
        }
    }
    fn add(&self, a: &Dv, b: &Dv) -> Dv {
        if self.is_zero(a) {
            return b.clone();
        }
        if self.is_zero(b) {
            return a.clone();
        }
        let k = &self.k;
        let (d1, e1, e2) = pxgcd(k, &a.u, &b.u);
        let vs = padd(k, &a.v, &b.v);
        let (d, c1, c2) = pxgcd(k, &d1, &vs);
        let s1 = pmul(k, &c1, &e1);
        let s2 = pmul(k, &c1, &e2);
        let d2 = pmul(k, &d, &d);
        let (u, r) = pdivrem(k, &pmul(k, &a.u, &b.u), &d2);
        debug_assert!(r.is_empty());
        let t1 = pmul(k, &pmul(k, &s1, &a.u), &b.v);
        let t2 = pmul(k, &pmul(k, &s2, &b.u), &a.v);
        let t3 = pmul(k, &c2, &padd(k, &pmul(k, &a.v, &b.v), &self.f));
        let num = padd(k, &padd(k, &t1, &t2), &t3);
        let (vq, vr) = pdivrem(k, &num, &d);
        debug_assert!(vr.is_empty());
        let mut u = pmonic(k, &u);
        let mut v = pdivrem(k, &vq, &u).1;
        while u.len() > 3 {
            let num = psub(k, &self.f, &pmul(k, &v, &v));
            let (u2, r) = pdivrem(k, &num, &u);
            debug_assert!(r.is_empty());
            u = pmonic(k, &u2);
            let nv: Poly = trim(v.iter().map(|&c| k.neg(c)).collect());
            v = pdivrem(k, &nv, &u).1;
        }
        Dv { u, v }
    }
    fn mul(&self, a: &Dv, mut n: u64) -> Dv {
        let mut acc = self.zero();
        let mut base = a.clone();
        while n > 0 {
            if n & 1 == 1 {
                acc = self.add(&acc, &base);
            }
            n >>= 1;
            if n > 0 {
                base = self.add(&base, &base);
            }
        }
        acc
    }
    /// `ζ: (x, y) ↦ (ζ₅x, y)` on Mumford pairs: `u ↦ ζ^d u(x/ζ)`, `v ↦ v(x/ζ)`.
    fn zeta(&self, a: &Dv) -> Dv {
        let k = &self.k;
        let d = a.u.len() - 1;
        let zi = k.inv(self.z5);
        let u =
            a.u.iter()
                .enumerate()
                .map(|(i, &c)| k.mul(c, k.pow(self.z5, (d - i) as u64)))
                .collect();
        let v = trim(
            a.v.iter()
                .enumerate()
                .map(|(i, &c)| k.mul(c, k.pow(zi, i as u64)))
                .collect(),
        );
        Dv { u, v }
    }
    /// All reduced classes with `u = x² + a1·x + a0`, each verified.
    fn deg2(&self, a1: u64, a0: u64, fails: &mut u64) -> Vec<Dv> {
        let k = &self.k;
        let half = k.inv(2);
        let u: Poly = vec![a0, a1, 1];
        let disc = k.sub(k.mul(a1, a1), k.mul(4, a0));
        let mut cand: Vec<Poly> = Vec::new();
        if disc == 0 {
            let e = k.mul(k.neg(a1), half);
            let fe = self.fx(e);
            if fe != 0 {
                if let Some(s) = k.sqrt(fe) {
                    // f'(e) = 5 e^4
                    let fd = k.mul(5, k.pow(e, 4));
                    for s in [s, k.neg(s)] {
                        let v1 = k.mul(fd, k.inv(k.mul(2, s)));
                        let v0 = k.sub(s, k.mul(v1, e));
                        cand.push(trim(vec![v0, v1]));
                    }
                }
            }
        } else if let Some(sq) = k.sqrt(disc) {
            let e1 = k.mul(k.add(k.neg(a1), sq), half);
            let e2 = k.mul(k.sub(k.neg(a1), sq), half);
            if let (Some(s1), Some(s2)) = (k.sqrt(self.fx(e1)), k.sqrt(self.fx(e2))) {
                let di = k.inv(k.sub(e1, e2));
                for y1 in [s1, k.neg(s1)] {
                    for y2 in [s2, k.neg(s2)] {
                        let v1 = k.mul(k.sub(y1, y2), di);
                        let v0 = k.sub(y1, k.mul(v1, e1));
                        cand.push(trim(vec![v0, v1]));
                    }
                }
            }
        } else {
            // θ = (−a1 + κ t)/2 with κ² = disc / w.
            let kap = k.sqrt(k.mul(disc, k.inv(k.w))).expect("disc/w is a square");
            let th = (k.mul(k.neg(a1), half), k.mul(kap, half));
            let th2 = k.m2(th, th);
            let th4 = k.m2(th2, th2);
            let th5 = k.m2(th4, th);
            let fth = (k.add(th5.0, self.b), th5.1);
            let roots: Vec<(u64, u64)> = if fth == (0, 0) {
                vec![(0, 0)]
            } else if k.is_qr(k.norm2(fth)) {
                match k.sqrt2(fth) {
                    Some(s) => vec![s, (k.neg(s.0), k.neg(s.1))],
                    None => {
                        *fails += 1;
                        Vec::new()
                    }
                }
            } else {
                Vec::new()
            };
            for s in roots {
                let v1 = k.mul(s.1, k.inv(th.1));
                let v0 = k.sub(s.0, k.mul(v1, th.0));
                cand.push(trim(vec![v0, v1]));
            }
        }
        let mut out: Vec<Dv> = Vec::new();
        for v in cand {
            let d = Dv { u: u.clone(), v };
            if !self.valid(&d) {
                *fails += 1;
                continue;
            }
            if !out.contains(&d) {
                out.push(d);
            }
        }
        out
    }
    /// Exhaustive count of reduced Mumford pairs (every `u` of degree ≤ 2).
    fn exhaustive_order(&self) -> (u64, u64) {
        let p = self.p();
        let mut fails = 0;
        let mut n = 1u64;
        for a in 0..p {
            let fa = self.fx(a);
            if fa == 0 {
                n += 1;
            } else if self.k.is_qr(fa) {
                n += 2;
            }
        }
        for a1 in 0..p {
            for a0 in 0..p {
                n += self.deg2(a1, a0, &mut fails).len() as u64;
            }
        }
        (n, fails)
    }
    /// `#C(F_p)`, `#C(F_{p²})` (one point at infinity) and `#J = L(1)`.
    fn l_poly_order(&self) -> (u64, u64, i128, i128, u64) {
        let k = &self.k;
        let p = self.p();
        let mut n1 = 1u64;
        for x in 0..p {
            n1 += (1 + k.chi[self.fx(x) as usize] as i64) as u64;
        }
        let mut n2 = 1u64;
        for x0 in 0..p {
            for x1 in 0..p {
                let x = (x0, x1);
                let x2 = k.m2(x, x);
                let x4 = k.m2(x2, x2);
                let x5 = k.m2(x4, x);
                let fx = (k.add(x5.0, self.b), x5.1);
                if fx == (0, 0) {
                    n2 += 1;
                } else {
                    n2 += (1 + k.chi[k.norm2(fx) as usize] as i64) as u64;
                }
            }
        }
        let (pi, n1i, n2i) = (p as i128, n1 as i128, n2 as i128);
        let s1 = pi + 1 - n1i;
        let s2 = pi * pi + 1 - n2i;
        let e2 = (s1 * s1 - s2) / 2;
        let nj = 1 - s1 + e2 - pi * s1 + pi * pi;
        (n1, n2, s1, e2, nj as u64)
    }
    /// A random class: a random `u` of degree 2 with a solution.
    fn random(&self, rng: &mut StdRng) -> Dv {
        let p = self.p();
        let mut fails = 0;
        loop {
            let s = self.deg2(rng.gen_range(0..p), rng.gen_range(0..p), &mut fails);
            if !s.is_empty() {
                return s[rng.gen_range(0..s.len())].clone();
            }
        }
    }
    /// Rational finite Weierstrass points.
    fn roots(&self) -> Vec<u64> {
        (0..self.p()).filter(|&x| self.fx(x) == 0).collect()
    }
    /// Both points over `x` (none if `f(x)` is a non-residue or zero).
    fn lift(&self, x: u64) -> Vec<Dv> {
        let fx = self.fx(x);
        if fx == 0 {
            return Vec::new();
        }
        match self.k.sqrt(fx) {
            Some(s) => vec![self.pt(x, s), self.pt(x, self.k.neg(s))],
            None => Vec::new(),
        }
    }
}

fn factor(mut n: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut d = 2;
    while d * d <= n {
        let mut e = 0;
        while n.is_multiple_of(d) {
            n /= d;
            e += 1;
        }
        if e > 0 {
            out.push((d, e));
        }
        d += 1;
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

// ── Search ────────────────────────────────────────────────────────────────

#[derive(Serialize)]
struct Cand {
    p: u64,
    b: u64,
    order: u64,
    h: u64,
    r: u64,
    rational_2_torsion: usize,
    b_is_square: bool,
}

fn search(lo: u64, hi: u64) -> Vec<Cand> {
    let mut out = Vec::new();
    for p in lo..hi {
        if p % 10 != 1 || !is_prime_u64(p) || p >= 1 << 14 {
            continue;
        }
        let k = Fp::new(p);
        let g = (2..p)
            .find(|&g| {
                factor(p - 1)
                    .iter()
                    .all(|&(q, _)| k.pow(g, (p - 1) / q) != 1)
            })
            .unwrap();
        for e in 0..10 {
            let b = k.pow(g, e);
            let jac = Jac::new(p, b);
            let (_, _, _, _, n) = jac.l_poly_order();
            let fs = factor(n);
            let (r, er) = *fs.last().unwrap();
            if er != 1 {
                continue;
            }
            out.push(Cand {
                p,
                b,
                order: n,
                h: n / r,
                r,
                rational_2_torsion: jac.roots().len(),
                b_is_square: k.is_qr(b),
            });
        }
    }
    out
}

// ── Census ────────────────────────────────────────────────────────────────

#[derive(Serialize, Clone)]
struct MStat {
    distinct: u64,
    /// Distinct hits from sums of pairwise distinct base elements only.
    distinct_strict: u64,
    representations: u64,
    hit_rate: f64,
    c: f64,
    cost_model: f64,
}

#[derive(Serialize, Clone)]
struct Census {
    arm: String,
    rep: String,
    points: usize,
    degree1_points: usize,
    columns: usize,
    classes: usize,
    m2: MStat,
    m3: MStat,
}

struct Ctx<'a> {
    jac: &'a Jac,
    r: u64,
    h: u64,
}

fn census(cx: &Ctx, arm: &str, rep: &str, base: &[Dv]) -> Census {
    let jac = cx.jac;
    let n = base.len();
    let keys: HashSet<u64> = base.iter().map(|d| jac.key(d)).collect();
    assert_eq!(keys.len(), n, "{arm}: duplicate base elements");
    for d in base {
        assert!(jac.valid(d), "{arm}: invalid base element");
        assert!(keys.contains(&jac.key(&jac.neg(d))), "{arm}: not signed");
    }
    let cls: Vec<Dv> = base.iter().map(|d| jac.mul(d, cx.r)).collect();
    let classes: HashSet<u64> = cls.iter().map(|d| jac.key(d)).collect();
    // columns: orbits of [h]D under A = <ζ, −1>
    let mut seen: HashSet<u64> = HashSet::new();
    let mut cols = 0;
    for d in base {
        let q = jac.mul(d, cx.h);
        if jac.is_zero(&q) || seen.contains(&jac.key(&q)) {
            continue;
        }
        cols += 1;
        let mut z = q;
        for _ in 0..5 {
            seen.insert(jac.key(&z));
            seen.insert(jac.key(&jac.neg(&z)));
            z = jac.zeta(&z);
        }
    }
    // pair table keyed by class sum: (j, k, sum), j ≤ k
    let mut pairs: HashMap<u64, Vec<(usize, usize, Dv)>> = HashMap::new();
    let zero_key = jac.key(&jac.zero());
    let mut hits2: HashSet<u64> = HashSet::new();
    let mut hits2s: HashSet<u64> = HashSet::new();
    let mut reps2 = 0u64;
    for i in 0..n {
        for j in i..n {
            let s = jac.add(&base[i], &base[j]);
            let c = jac.key(&jac.add(&cls[i], &cls[j]));
            if c == zero_key && !jac.is_zero(&s) {
                hits2.insert(jac.key(&s));
                if i < j {
                    hits2s.insert(jac.key(&s));
                }
                reps2 += 1;
            }
            pairs.entry(c).or_default().push((i, j, s));
        }
    }
    let mut hits3: HashSet<u64> = HashSet::new();
    let mut hits3s: HashSet<u64> = HashSet::new();
    let mut reps3 = 0u64;
    for i in 0..n {
        let need = jac.key(&jac.neg(&cls[i]));
        if let Some(v) = pairs.get(&need) {
            for (j, k, s) in v {
                if *j < i {
                    continue;
                }
                let t = jac.add(s, &base[i]);
                if !jac.is_zero(&t) {
                    hits3.insert(jac.key(&t));
                    if i < *j && j < k {
                        hits3s.insert(jac.key(&t));
                    }
                    reps3 += 1;
                }
            }
        }
    }
    let mk = |d: usize, ds: usize, reps: u64, fact: f64, m: i32, l: f64| MStat {
        distinct: d as u64,
        distinct_strict: ds as u64,
        representations: reps,
        hit_rate: d as f64 / cx.r as f64,
        c: d as f64 * fact / (n as f64).powi(m),
        cost_model: if d > 0 {
            cols as f64 * (cx.r as f64 / d as f64) * l
        } else {
            f64::NAN
        },
    };
    Census {
        arm: arm.to_string(),
        rep: rep.to_string(),
        points: n,
        degree1_points: base.iter().filter(|d| d.u.len() == 2).count(),
        columns: cols,
        classes: classes.len(),
        m2: mk(hits2.len(), hits2s.len(), reps2, 2.0, 2, 1.0),
        m3: mk(hits3.len(), hits3s.len(), reps3, 6.0, 3, n as f64),
    }
}

#[derive(Serialize)]
struct Validation {
    n1: u64,
    n2: u64,
    l_s1: i128,
    l_e2: i128,
    order_lpoly: u64,
    order_exhaustive: Option<u64>,
    exhaustive_solver_failures: u64,
    random_order_checks: usize,
    assoc_checks: usize,
    library_crosschecks: usize,
    zeta_hom_checks: usize,
    zeta_eigenvalue_on_g: u64,
    rational_2_torsion_points: usize,
    two_sylow_order: u64,
    index_2j: u64,
}

#[derive(Serialize)]
struct Out {
    p: u64,
    b: u64,
    order: u64,
    factorisation: Vec<(u64, u32)>,
    r: u64,
    h: u64,
    zeta5: u64,
    kernels: Vec<(String, usize)>,
    validation: Validation,
    size: usize,
    reps: usize,
    seed: u64,
    summary: Vec<ArmSummary>,
    census: Vec<Census>,
}

fn crosscheck_library(jac: &Jac, rng: &mut StdRng, n: usize) -> usize {
    let fq = Fq::new(jac.p());
    let lift = |a: &Poly| -> Vec<E2> { a.iter().map(|&c| fq.from_fp(c)).collect() };
    let h: Vec<E2> = lift(&jac.f);
    let hyp = Hyp {
        f: &fq,
        h: h.clone(),
        g: 2,
    };
    let mut ok = 0;
    for _ in 0..n {
        let a = jac.random(rng);
        let b = jac.random(rng);
        let s = jac.add(&a, &b);
        let la = LDiv {
            u: lift(&a.u),
            v: lift(&a.v),
        };
        let lb = LDiv {
            u: lift(&b.u),
            v: lift(&b.v),
        };
        let ls = hyp.add(&la, &lb);
        let mine = LDiv {
            u: lift(&s.u),
            v: lift(&s.v),
        };
        let norm = |d: &LDiv<E2>| -> (Vec<E2>, Vec<E2>) {
            let mut u = d.u.clone();
            while u.last() == Some(&fq.zero()) {
                u.pop();
            }
            let mut v = d.v.clone();
            while v.last() == Some(&fq.zero()) {
                v.pop();
            }
            (u, v)
        };
        assert_eq!(
            norm(&ls),
            norm(&mine),
            "Cantor disagrees with jv_cover::Hyp"
        );
        ok += 1;
    }
    ok
}

/// Subgroup generated by `gens` (closure under addition), as keys.
fn closure(jac: &Jac, gens: &[Dv]) -> HashMap<u64, Dv> {
    let mut set: HashMap<u64, Dv> = HashMap::new();
    set.insert(jac.key(&jac.zero()), jac.zero());
    for g in gens {
        if set.contains_key(&jac.key(g)) {
            continue;
        }
        let cur: Vec<Dv> = set.values().cloned().collect();
        let mut mult = g.clone();
        let mut all = cur.clone();
        while !set.contains_key(&jac.key(&mult)) {
            for c in &cur {
                let s = jac.add(c, &mult);
                all.push(s);
            }
            mult = jac.add(&mult, g);
        }
        for d in all {
            set.insert(jac.key(&d), d);
        }
    }
    set
}

fn pick_points<F: Fn(&Dv) -> bool>(
    jac: &Jac,
    want: usize,
    rng: Option<&mut StdRng>,
    used5: &mut HashSet<u64>,
    distinct_orbit: bool,
    keep: F,
) -> Option<Vec<Dv>> {
    let p = jac.p();
    let mut out = Vec::new();
    let mut seen_x: HashSet<u64> = HashSet::new();
    let mut xs: Vec<u64> = (1..p).collect();
    if let Some(rng) = rng {
        for i in (1..xs.len()).rev() {
            let j = rng.gen_range(0..=i);
            xs.swap(i, j);
        }
    }
    for x in xs {
        if out.len() >= want {
            break;
        }
        if seen_x.contains(&x) {
            continue;
        }
        let x5 = jac.k.pow(x, 5);
        if distinct_orbit && used5.contains(&x5) {
            continue;
        }
        let l = jac.lift(x);
        if l.len() != 2 || !keep(&l[0]) {
            continue;
        }
        seen_x.insert(x);
        used5.insert(x5);
        out.extend(l);
    }
    if out.len() < want {
        return None;
    }
    out.truncate(want);
    Some(out)
}

/// A ratio of means over the random bases, with a 2000-resample
/// percentile bootstrap interval over bases, and the value on the
/// deterministic small-x base.
#[derive(Serialize)]
struct Stat {
    value: f64,
    ci95: (f64, f64),
    smallx: f64,
}

#[derive(Serialize)]
struct ArmSummary {
    arm: String,
    /// The law's predictions relative to the plain base: `[U, c₂, c₃, cost₂, cost₃]`.
    predicted: [f64; 5],
    n_reps: usize,
    columns: Stat,
    c2: Stat,
    c3: Stat,
    /// Same as `c2`, `c3`, counting only sums of pairwise distinct base elements.
    c2_strict: Stat,
    c3_strict: Stat,
    cost2: Stat,
    cost3: Stat,
    /// Cost ratios with the plain base's columns set to their nominal
    /// `|F|/2` (negation folding only), removing the plain base's
    /// accidental `ζ`-coincidences.
    cost2_nominal: Stat,
    cost3_nominal: Stat,
    hit_rate3: Stat,
    /// Plain base: columns as a fraction of `|F|/2`.
    plain_columns_frac: Stat,
    /// Plain base: fraction of the m = 2 hits that are doublings `2P`.
    plain_diag_frac2: Stat,
}

/// Predicted ratios to plain: `[U, c₂, c₃, cost₂, cost₃]`.
fn predicted(arm: &str, h: u64, index_2j: u64) -> [f64; 5] {
    let law = |a: f64, mu: f64, i: f64| {
        [
            1.0 / (a * mu),
            i / mu,
            i / (mu * mu),
            1.0 / (a * i),
            mu / (a * i),
        ]
    };
    match arm {
        "A" => law(5.0, 1.0, 1.0),
        "M2" => law(1.0, 2.0, 1.0),
        "M4iso" | "M4non" => law(1.0, 4.0, 1.0),
        "M5" => law(1.0, 5.0, 1.0),
        "AM5" => law(5.0, 5.0, 1.0),
        "2J" => law(1.0, 1.0, index_2j as f64),
        "G" => law(1.0, 1.0, h as f64),
        _ => [1.0; 5],
    }
}

fn summarise(res: &[Census], h: u64, index_2j: u64, seed: u64) -> Vec<ArmSummary> {
    let mut arms: Vec<String> = Vec::new();
    for c in res {
        if c.arm != "plain" && !arms.contains(&c.arm) {
            arms.push(c.arm.clone());
        }
    }
    let plain = |rep: &str| res.iter().find(|c| c.arm == "plain" && c.rep == rep);
    // per base: [U_a, D2_a, D3_a, D2s_a, D3s_a, hit3_a, U_p, D2_p, D3_p, D2s_p, D3s_p, |F|/2]
    let row = |c: &Census, p: &Census| -> [f64; 12] {
        [
            c.columns as f64,
            c.m2.distinct as f64,
            c.m3.distinct as f64,
            c.m2.distinct_strict as f64,
            c.m3.distinct_strict as f64,
            c.m3.hit_rate,
            p.columns as f64,
            p.m2.distinct as f64,
            p.m3.distinct as f64,
            p.m2.distinct_strict as f64,
            p.m3.distinct_strict as f64,
            p.points as f64 / 2.0,
        ]
    };
    // metrics from (aggregated) rows
    let metrics = |m: &[f64; 12]| -> [f64; 12] {
        let u = m[0] / m[6];
        [
            u,
            m[1] / m[7],
            m[2] / m[8],
            m[3] / m[9],
            m[4] / m[10],
            u * m[7] / m[1],
            u * m[8] / m[2],
            (m[0] / m[11]) * m[7] / m[1],
            (m[0] / m[11]) * m[8] / m[2],
            m[5],
            m[6] / m[11],
            1.0 - m[9] / m[7],
        ]
    };
    let mean_rows = |rows: &[[f64; 12]], idx: &[usize]| -> [f64; 12] {
        let mut s = [0.0; 12];
        for &i in idx {
            for k in 0..12 {
                s[k] += rows[i][k];
            }
        }
        s.map(|x| x / idx.len() as f64)
    };
    let mut rng = StdRng::seed_from_u64(seed ^ 0x5eed);
    let mut out = Vec::new();
    for arm in arms {
        let mut rows: Vec<[f64; 12]> = Vec::new();
        let mut sx = [f64::NAN; 12];
        for c in res.iter().filter(|c| c.arm == arm) {
            let Some(p) = plain(&c.rep) else { continue };
            if c.rep == "smallx" {
                sx = metrics(&row(c, p));
            } else {
                rows.push(row(c, p));
            }
        }
        if rows.is_empty() {
            continue;
        }
        let all: Vec<usize> = (0..rows.len()).collect();
        let point = metrics(&mean_rows(&rows, &all));
        let mut boots: Vec<[f64; 12]> = Vec::new();
        for _ in 0..2000 {
            let idx: Vec<usize> = (0..rows.len())
                .map(|_| rng.gen_range(0..rows.len()))
                .collect();
            boots.push(metrics(&mean_rows(&rows, &idx)));
        }
        let col = |k: usize| {
            let mut v: Vec<f64> = boots.iter().map(|b| b[k]).collect();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            Stat {
                value: point[k],
                ci95: (v[50], v[1949]),
                smallx: sx[k],
            }
        };
        out.push(ArmSummary {
            predicted: predicted(&arm, h, index_2j),
            arm,
            n_reps: rows.len(),
            columns: col(0),
            c2: col(1),
            c3: col(2),
            c2_strict: col(3),
            c3_strict: col(4),
            cost2: col(5),
            cost3: col(6),
            cost2_nominal: col(7),
            cost3_nominal: col(8),
            hit_rate3: col(9),
            plain_columns_frac: col(10),
            plain_diag_frac2: col(11),
        });
    }
    out
}

fn arg(args: &[String], name: &str, default: &str) -> String {
    args.iter()
        .position(|a| a == name)
        .map(|i| args[i + 1].clone())
        .unwrap_or_else(|| default.to_string())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mode = args.get(1).cloned().unwrap_or_default();
    if mode == "search" {
        let lo: u64 = arg(&args, "--lo", "1000").parse().unwrap();
        let hi: u64 = arg(&args, "--hi", "1100").parse().unwrap();
        let maxh: u64 = arg(&args, "--maxh", "200").parse().unwrap();
        for c in search(lo, hi) {
            if c.h <= maxh {
                println!("{}", serde_json::to_string(&c).unwrap());
            }
        }
        return;
    }
    assert_eq!(mode, "census", "usage: search | census");
    let p: u64 = arg(&args, "--p", "1051").parse().unwrap();
    let b: u64 = arg(&args, "--b", "1").parse().unwrap();
    let size: usize = arg(&args, "--size", "80").parse().unwrap();
    let reps: usize = arg(&args, "--reps", "20").parse().unwrap();
    let seed: u64 = arg(&args, "--seed", "1").parse().unwrap();
    let exh_max: u64 = arg(&args, "--exhaustive-max-p", "5000").parse().unwrap();
    let arms_arg = arg(&args, "--arms", "plain,A,M2,M4iso,M4non,M5,AM5,2J,G");
    let json = arg(&args, "--json", "");
    assert!(p % 10 == 1 && is_prime_u64(p) && p < 1 << 14);
    let jac = Jac::new(p, b);
    let mut rng = StdRng::seed_from_u64(seed);

    // ── validation ──
    let (n1, n2, s1, e2, order) = jac.l_poly_order();
    let (order_exh, fails) = if p <= exh_max {
        let (n, f) = jac.exhaustive_order();
        (Some(n), f)
    } else {
        (None, 0)
    };
    if let Some(n) = order_exh {
        assert_eq!(n, order, "exhaustive count disagrees with L(1)");
    }
    assert_eq!(fails, 0);
    let fs = factor(order);
    let (r, er) = *fs.last().unwrap();
    assert_eq!(er, 1, "r² divides #J");
    let h = order / r;
    let mut rc = 0;
    for _ in 0..40 {
        let d = jac.random(&mut rng);
        assert!(jac.is_zero(&jac.mul(&d, order)), "[#J]D != 0");
        rc += 1;
    }
    let mut ac = 0;
    for _ in 0..300 {
        let (a, b2, c) = (
            jac.random(&mut rng),
            jac.random(&mut rng),
            jac.random(&mut rng),
        );
        assert_eq!(
            jac.add(&jac.add(&a, &b2), &c),
            jac.add(&a, &jac.add(&b2, &c))
        );
        assert_eq!(jac.add(&a, &b2), jac.add(&b2, &a));
        assert!(jac.is_zero(&jac.add(&a, &jac.neg(&a))));
        ac += 1;
    }
    let lc = crosscheck_library(&jac, &mut rng, 300);
    let mut zc = 0;
    for _ in 0..200 {
        let (a, b2) = (jac.random(&mut rng), jac.random(&mut rng));
        assert_eq!(
            jac.zeta(&jac.add(&a, &b2)),
            jac.add(&jac.zeta(&a), &jac.zeta(&b2))
        );
        let mut z = a.clone();
        for _ in 0..5 {
            z = jac.zeta(&z);
        }
        assert_eq!(z, a);
        zc += 1;
    }
    // generator of G and the eigenvalue of ζ on it
    let g = loop {
        let d = jac.mul(&jac.random(&mut rng), h);
        if !jac.is_zero(&d) {
            break d;
        }
    };
    assert!(jac.is_zero(&jac.mul(&g, r)));
    let zg = jac.zeta(&g);
    let lam = {
        let kk = Fp::new(r);
        let mut found = 0;
        let gen = (2..r)
            .map(|x| kk.pow(x, (r - 1) / 5))
            .find(|&z| z != 1)
            .expect("r ≡ 1 mod 5");
        let mut l = gen;
        for _ in 0..4 {
            if jac.mul(&g, l) == zg {
                found = l;
            }
            l = kk.mul(l, gen);
        }
        assert!(found != 0, "ζ is not a scalar on G");
        found
    };
    let roots = jac.roots();
    // 2-Sylow
    let e2v = fs.iter().find(|f| f.0 == 2).map(|f| f.1).unwrap_or(0);
    let odd = order >> e2v;
    let mut gens = Vec::new();
    let mut syl: HashMap<u64, Dv> = closure(&jac, &[]);
    while (syl.len() as u64) < (1u64 << e2v) {
        gens.push(jac.mul(&jac.random(&mut rng), odd));
        syl = closure(&jac, &gens);
    }
    let two_syl: HashSet<u64> = syl.values().map(|d| jac.key(&jac.add(d, d))).collect();
    let index_2j = syl.len() as u64 / two_syl.len() as u64;
    let in_2j = |d: &Dv| two_syl.contains(&jac.key(&jac.mul(d, odd)));
    let in_g = |d: &Dv| jac.is_zero(&jac.mul(d, r));

    // ── kernels ──
    let w: Vec<Dv> = roots.iter().map(|&e| jac.pt(e, 0)).collect();
    let mut kernels: Vec<(String, Vec<Dv>)> = Vec::new();
    if !w.is_empty() {
        kernels.push(("M2".into(), vec![jac.zero(), w[0].clone()]));
    }
    if w.len() == 5 {
        // Richelot kernel {e1,e2},{e3,e4},{e5,∞}: maximal isotropic.
        let t1 = w[4].clone();
        let t2 = jac.add(&w[0], &w[1]);
        let t12 = jac.add(&t1, &t2);
        assert_eq!(t12, jac.add(&w[2], &w[3]));
        kernels.push(("M4iso".into(), vec![jac.zero(), t1, t2, t12]));
        // {e1 − ∞, e2 − ∞}: Weil pairing −1, not isotropic.
        let s = jac.add(&w[0], &w[1]);
        kernels.push((
            "M4non".into(),
            vec![jac.zero(), w[0].clone(), w[1].clone(), s],
        ));
    }
    if let Some(s) = jac.k.sqrt(jac.b) {
        let t5 = jac.pt(0, s);
        let mut m = vec![jac.zero()];
        for _ in 0..4 {
            let nx = jac.add(m.last().unwrap(), &t5);
            m.push(nx);
        }
        assert!(jac.is_zero(&jac.add(m.last().unwrap(), &t5)));
        assert_eq!(jac.zeta(&t5), t5);
        kernels.push(("M5".into(), m));
    }
    for (name, m) in &kernels {
        let ks: HashSet<u64> = m.iter().map(|d| jac.key(d)).collect();
        for a in m {
            for b2 in m {
                assert!(ks.contains(&jac.key(&jac.add(a, b2))), "{name} not closed");
            }
        }
    }
    let kern = |name: &str| kernels.iter().find(|k| k.0 == name).map(|k| k.1.clone());

    let cx = Ctx { jac: &jac, r, h };
    let arms: Vec<String> = arms_arg.split(',').map(|s| s.to_string()).collect();
    let mut results = Vec::new();
    for rep in 0..=reps {
        // rep 0 is the deterministic small-x base of §15.2
        let rep_name = if rep == 0 {
            "smallx".to_string()
        } else {
            format!("rand{rep}")
        };
        let mut rr = StdRng::seed_from_u64(seed * 1_000_003 + rep as u64);
        for arm in &arms {
            let mut used = HashSet::new();
            let rng_opt = if rep == 0 { None } else { Some(&mut rr) };
            let base: Option<Vec<Dv>> = (|| match arm.as_str() {
                "plain" => pick_points(&jac, size, rng_opt, &mut used, false, |_| true),
                "A" => {
                    if !size.is_multiple_of(10) {
                        None
                    } else {
                        let s = pick_points(&jac, size / 5, rng_opt, &mut used, true, |_| true)?;
                        let mut b = Vec::new();
                        for d in &s {
                            let mut z = d.clone();
                            for _ in 0..5 {
                                b.push(z.clone());
                                z = jac.zeta(&z);
                            }
                        }
                        Some(b)
                    }
                }
                "M2" | "M4iso" | "M4non" | "M5" => match kern(arm) {
                    Some(m) if size.is_multiple_of(2 * m.len()) => {
                        let s =
                            pick_points(&jac, size / m.len(), rng_opt, &mut used, false, |_| true)?;
                        let mut b = Vec::new();
                        for t in &m {
                            for d in &s {
                                b.push(jac.add(d, t));
                            }
                        }
                        Some(b)
                    }
                    _ => None,
                },
                "AM5" => match kern("M5") {
                    Some(m) if size.is_multiple_of(50) => {
                        let s = pick_points(&jac, size / 25, rng_opt, &mut used, true, |_| true)?;
                        let mut b = Vec::new();
                        for d in &s {
                            let mut z = d.clone();
                            for _ in 0..5 {
                                for t in &m {
                                    b.push(jac.add(&z, t));
                                }
                                z = jac.zeta(&z);
                            }
                        }
                        Some(b)
                    }
                    _ => None,
                },
                "2J" => {
                    if index_2j > 1 && index_2j != h {
                        pick_points(&jac, size, rng_opt, &mut used, false, in_2j)
                    } else {
                        None
                    }
                }
                "G" => {
                    if h > 1 {
                        pick_points(&jac, size, rng_opt, &mut used, false, in_g)
                    } else {
                        None
                    }
                }
                other => panic!("unknown arm {other}"),
            })();
            if let Some(base) = base {
                let c = census(&cx, arm, &rep_name, &base);
                eprintln!(
                    "{rep_name:>7} {arm:>6} |F|={} U={} cls={} D2={} D3={} c2={:.4} c3={:.4}",
                    c.points, c.columns, c.classes, c.m2.distinct, c.m3.distinct, c.m2.c, c.m3.c
                );
                results.push(c);
            }
        }
    }
    let out = Out {
        p,
        b: jac.b,
        order,
        factorisation: fs.clone(),
        r,
        h,
        zeta5: jac.z5,
        kernels: kernels.iter().map(|k| (k.0.clone(), k.1.len())).collect(),
        validation: Validation {
            n1,
            n2,
            l_s1: s1,
            l_e2: e2,
            order_lpoly: order,
            order_exhaustive: order_exh,
            exhaustive_solver_failures: fails,
            random_order_checks: rc,
            assoc_checks: ac,
            library_crosschecks: lc,
            zeta_hom_checks: zc,
            zeta_eigenvalue_on_g: lam,
            rational_2_torsion_points: roots.len(),
            two_sylow_order: syl.len() as u64,
            index_2j,
        },
        size,
        reps,
        seed,
        summary: summarise(&results, h, index_2j, seed),
        census: results,
    };
    let s = serde_json::to_string_pretty(&out).unwrap();
    if json.is_empty() {
        println!("{s}");
    } else {
        fs::write(&json, s).unwrap();
    }
}
