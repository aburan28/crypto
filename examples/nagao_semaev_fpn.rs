//! Nagao's Riemann–Roch decomposition against Semaev's symmetrised summation
//! polynomial on `E / F_{p^n}`, odd characteristic, under one F4 engine.
//!
//! Companion to `research/ic_nagao_semaev_fpn_20261016/PROTOCOL.md`.
//!
//! Factor base `F = {P ∈ E(F_{p^n}) : x(P) ∈ F_p}` (Gaudry, Diem).  A target
//! `R` decomposes as `R = P₁ + … + P_m`, `P_i ∈ F`, when some function
//! `f ∈ L((m+1)·O)` vanishes at `P₁, …, P_m, −R`.  Two systems over `F_p`
//! decide that, both solved by `jv_cover::solve_zero_dim` (grevlex F4, then
//! multiplication matrices) with the same options:
//!
//! * `sem`: `S_{m+1}(x₁, …, x_m, x_R)` rewritten in the elementary symmetric
//!   functions `e₁, …, e_m ∈ F_p` of the summand abscissae and Weil-descended:
//!   `m` unknowns, `n` equations, total degree `≤ 2^{m−1}`.  The symmetrised
//!   polynomial is built once per curve (tensor-grid interpolation of nested
//!   resultants, then the fundamental theorem of symmetric polynomials) and
//!   specialised at `x_R` per target.
//! * `rr`: Nagao's form.  `f = φ_lead + Σ c_j φ_j` monic in the basis of
//!   `L((m+1)·O)`, the constant coefficient eliminated by `f(−R) = 0`, the
//!   other `m − 1` coefficients Weil-restricted (`n(m−1)` unknowns over
//!   `F_p`).  The norm `A² − B²·h = (x − x_R)·g(x)` with `g` monic of degree
//!   `m`; the system asks that every coefficient of `g` lie in `F_p`:
//!   `m(n−1)` quadratic equations.
//!
//! Every solution of either system is turned into the polynomial
//! `g = ∏(x − x_i)` and checked by one verifier: `g` splits over `F_p`, each
//! `h(x_i)` is a square in `F_{p^n}`, and some choice of signs gives
//! `Σ ±P_i = R` on the curve.
//!
//! Modes:
//! * `--validate --p P`: small-`p` cross-check of both arms against brute
//!   force over all `m`-multisets of the factor base, every cell.
//! * `--n N --m M --p P --seed S --targets T --kind planted|random
//!   --budget-s B --dmax D --out FILE`: one cell, one JSON line per
//!   (target, arm).

use crypto_lib::cryptanalysis::f4_fp::{self, F4Options, Ordering, Poly, Verdict};
use crypto_lib::cryptanalysis::jv_cover::{solve_zero_dim, ZeroDimStats};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Write;
use std::time::{Duration, Instant};

// ── F_p ─────────────────────────────────────────────────────────────

fn mulp(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

fn powp(mut b: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1 % p;
    b %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mulp(r, b, p);
        }
        b = mulp(b, b, p);
        e >>= 1;
    }
    r
}

fn invp(a: u64, p: u64) -> u64 {
    assert!(a % p != 0, "inverse of zero");
    powp(a, p - 2, p)
}

// ── F_q = F_p[t]/(μ), q = p^n, n ≤ 5 ────────────────────────────────

const NMAX: usize = 5;

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default, PartialOrd, Ord)]
struct Fq([u64; NMAX]);

#[derive(Clone, Debug)]
struct Field {
    p: u64,
    n: usize,
    /// Monic modulus `μ = t^n + Σ mu[i] t^i`.
    mu: Vec<u64>,
    q: u128,
}

impl Field {
    fn new(p: u64, n: usize, rng: &mut StdRng) -> Field {
        loop {
            let mu: Vec<u64> = (0..n).map(|_| rng.gen_range(0..p)).collect();
            if mu[0] == 0 && n > 1 {
                continue;
            }
            if irreducible(&mu, p) {
                return Field {
                    p,
                    n,
                    mu,
                    q: (p as u128).pow(n as u32),
                };
            }
        }
    }
    fn zero(&self) -> Fq {
        Fq::default()
    }
    fn one(&self) -> Fq {
        self.c(1)
    }
    fn c(&self, v: u64) -> Fq {
        let mut r = Fq::default();
        r.0[0] = v % self.p;
        r
    }
    fn theta_pow(&self, s: usize) -> Fq {
        let mut r = Fq::default();
        r.0[s] = 1;
        r
    }
    fn is_zero(&self, a: Fq) -> bool {
        a == Fq::default()
    }
    fn in_fp(&self, a: Fq) -> bool {
        a.0[1..].iter().all(|&v| v == 0)
    }
    fn add(&self, a: Fq, b: Fq) -> Fq {
        let mut r = Fq::default();
        for i in 0..self.n {
            r.0[i] = (a.0[i] + b.0[i]) % self.p;
        }
        r
    }
    fn sub(&self, a: Fq, b: Fq) -> Fq {
        let mut r = Fq::default();
        for i in 0..self.n {
            r.0[i] = (a.0[i] + self.p - b.0[i]) % self.p;
        }
        r
    }
    fn neg(&self, a: Fq) -> Fq {
        self.sub(self.zero(), a)
    }
    fn scale(&self, a: Fq, s: u64) -> Fq {
        let mut r = Fq::default();
        for i in 0..self.n {
            r.0[i] = mulp(a.0[i], s, self.p);
        }
        r
    }
    fn mul(&self, a: Fq, b: Fq) -> Fq {
        let (n, p) = (self.n, self.p);
        let mut acc = [0u128; 2 * NMAX];
        for i in 0..n {
            if a.0[i] == 0 {
                continue;
            }
            for j in 0..n {
                acc[i + j] += a.0[i] as u128 * b.0[j] as u128;
            }
        }
        let mut r: Vec<u64> = acc[..2 * n]
            .iter()
            .map(|&v| (v % p as u128) as u64)
            .collect();
        for d in (n..2 * n - 1).rev() {
            let c = r[d];
            if c == 0 {
                continue;
            }
            r[d] = 0;
            for i in 0..n {
                let k = d - n + i;
                r[k] = (r[k] + p - mulp(c, self.mu[i], p)) % p;
            }
        }
        let mut out = Fq::default();
        out.0[..n].copy_from_slice(&r[..n]);
        out
    }
    fn sqr(&self, a: Fq) -> Fq {
        self.mul(a, a)
    }
    fn pow(&self, mut b: Fq, mut e: u128) -> Fq {
        let mut r = self.one();
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, b);
            }
            b = self.sqr(b);
            e >>= 1;
        }
        r
    }
    fn inv(&self, a: Fq) -> Fq {
        assert!(!self.is_zero(a), "inverse of zero in F_q");
        self.pow(a, self.q - 2)
    }
    fn is_square(&self, a: Fq) -> bool {
        self.is_zero(a) || self.pow(a, (self.q - 1) / 2) == self.one()
    }
    /// Tonelli–Shanks; `None` for a non-square.
    fn sqrt(&self, a: Fq, rng: &mut StdRng) -> Option<Fq> {
        if self.is_zero(a) {
            return Some(a);
        }
        if !self.is_square(a) {
            return None;
        }
        let mut s = 0u32;
        let mut t = self.q - 1;
        while t % 2 == 0 {
            t /= 2;
            s += 1;
        }
        let z = loop {
            let z = self.random(rng);
            if !self.is_square(z) {
                break z;
            }
        };
        let mut m = s;
        let mut c = self.pow(z, t);
        let mut x = self.pow(a, (t + 1) / 2);
        let mut b = self.pow(a, t);
        while b != self.one() {
            let mut i = 0;
            let mut b2 = b;
            while b2 != self.one() {
                b2 = self.sqr(b2);
                i += 1;
            }
            let mut g = c;
            for _ in 0..(m - i - 1) {
                g = self.sqr(g);
            }
            x = self.mul(x, g);
            c = self.sqr(g);
            b = self.mul(b, c);
            m = i;
        }
        Some(x)
    }
    fn random(&self, rng: &mut StdRng) -> Fq {
        let mut r = Fq::default();
        for i in 0..self.n {
            r.0[i] = rng.gen_range(0..self.p);
        }
        r
    }
}

// Univariate polynomials over F_p (low degree first), for the irreducibility test.
fn up_trim(mut a: Vec<u64>) -> Vec<u64> {
    while a.last() == Some(&0) {
        a.pop();
    }
    a
}

fn up_rem(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
    let mut r = up_trim(a.to_vec());
    let b = up_trim(b.to_vec());
    let inv = invp(*b.last().unwrap(), p);
    while r.len() >= b.len() {
        let c = mulp(*r.last().unwrap(), inv, p);
        let sh = r.len() - b.len();
        for (i, &bi) in b.iter().enumerate() {
            r[sh + i] = (r[sh + i] + p - mulp(c, bi, p)) % p;
        }
        r = up_trim(r);
    }
    r
}

fn up_gcd(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
    let (mut a, mut b) = (up_trim(a.to_vec()), up_trim(b.to_vec()));
    while !b.is_empty() {
        let r = up_rem(&a, &b, p);
        a = b;
        b = r;
    }
    a
}

/// Rabin's test for `t^n + Σ mu[i] t^i` over `F_p`.
fn irreducible(mu: &[u64], p: u64) -> bool {
    let n = mu.len();
    if n == 1 {
        return true;
    }
    let f = Field {
        p,
        n,
        mu: mu.to_vec(),
        q: (p as u128).pow(n as u32),
    };
    let mut full: Vec<u64> = mu.to_vec();
    full.push(1);
    let t = f.theta_pow(1);
    let mut tp = t;
    for i in 1..=n {
        tp = f.pow(tp, p as u128);
        let diff = f.sub(tp, t);
        let g = up_gcd(&full, &diff.0[..n], p);
        if i < n && g.len() > 1 {
            return false;
        }
        if i == n && !diff.0[..n].iter().all(|&v| v == 0) {
            return false;
        }
    }
    true
}

// ── The curve ───────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Pt {
    Inf,
    A(Fq, Fq),
}

struct Curve<'a> {
    f: &'a Field,
    a: Fq,
    b: Fq,
}

impl<'a> Curve<'a> {
    fn h(&self, x: Fq) -> Fq {
        let f = self.f;
        f.add(f.mul(f.add(f.sqr(x), self.a), x), self.b)
    }
    fn neg(&self, p: Pt) -> Pt {
        match p {
            Pt::Inf => Pt::Inf,
            Pt::A(x, y) => Pt::A(x, self.f.neg(y)),
        }
    }
    fn add(&self, p: Pt, q: Pt) -> Pt {
        let f = self.f;
        match (p, q) {
            (Pt::Inf, _) => q,
            (_, Pt::Inf) => p,
            (Pt::A(x1, y1), Pt::A(x2, y2)) => {
                let lam = if x1 == x2 {
                    if f.add(y1, y2) == f.zero() {
                        return Pt::Inf;
                    }
                    let num = f.add(f.scale(f.sqr(x1), 3), self.a);
                    f.mul(num, f.inv(f.scale(y1, 2)))
                } else {
                    f.mul(f.sub(y2, y1), f.inv(f.sub(x2, x1)))
                };
                let x3 = f.sub(f.sub(f.sqr(lam), x1), x2);
                let y3 = f.sub(f.mul(lam, f.sub(x1, x3)), y1);
                Pt::A(x3, y3)
            }
        }
    }
}

// ── Multivariate polynomials over F_q in F_p-valued unknowns ────────

type Mp = BTreeMap<Vec<u32>, Fq>;

fn mp_add(f: &Field, a: &Mp, b: &Mp) -> Mp {
    let mut r = a.clone();
    for (e, c) in b {
        let v = r.entry(e.clone()).or_insert(f.zero());
        *v = f.add(*v, *c);
    }
    r.retain(|_, c| !f.is_zero(*c));
    r
}

fn mp_scale(f: &Field, a: &Mp, s: Fq) -> Mp {
    let mut r = Mp::new();
    if f.is_zero(s) {
        return r;
    }
    for (e, c) in a {
        r.insert(e.clone(), f.mul(*c, s));
    }
    r
}

fn mp_mul(f: &Field, a: &Mp, b: &Mp) -> Mp {
    let mut r = Mp::new();
    for (ea, ca) in a {
        for (eb, cb) in b {
            let e: Vec<u32> = ea.iter().zip(eb).map(|(x, y)| x + y).collect();
            let v = r.entry(e).or_insert(f.zero());
            *v = f.add(*v, f.mul(*ca, *cb));
        }
    }
    r.retain(|_, c| !f.is_zero(*c));
    r
}

fn mp_const(f: &Field, nv: usize, c: Fq) -> Mp {
    let mut r = Mp::new();
    if !f.is_zero(c) {
        r.insert(vec![0; nv], c);
    }
    r
}

/// Component `k` of every coefficient, as an `F_p` polynomial (grevlex, monic).
fn component(f: &Field, a: &Mp, k: usize) -> Poly {
    let terms: Vec<(Vec<u32>, u64)> = a
        .iter()
        .filter(|(_, c)| c.0[k] != 0)
        .map(|(e, c)| (e.clone(), c.0[k]))
        .collect();
    f4_fp::normalise(&terms, f.p, Ordering::Grevlex)
}

fn mp_eval_fp(f: &Field, a: &Mp, x: &[u64]) -> Fq {
    let mut acc = f.zero();
    for (e, c) in a {
        let mut m = 1u64;
        for (xi, &ei) in x.iter().zip(e) {
            m = mulp(m, powp(*xi, ei as u64, f.p), f.p);
        }
        acc = f.add(acc, f.scale(*c, m));
    }
    acc
}

// Polynomials in a curve variable `x` with `Mp` coefficients (low degree first).
type XPoly = Vec<Mp>;

fn xp_mul(f: &Field, a: &XPoly, b: &XPoly) -> XPoly {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r: XPoly = vec![Mp::new(); a.len() + b.len() - 1];
    for (i, ai) in a.iter().enumerate() {
        if ai.is_empty() {
            continue;
        }
        for (j, bj) in b.iter().enumerate() {
            if bj.is_empty() {
                continue;
            }
            let prod = mp_mul(f, ai, bj);
            r[i + j] = mp_add(f, &r[i + j], &prod);
        }
    }
    r
}

fn xp_sub(f: &Field, a: &XPoly, b: &XPoly) -> XPoly {
    let len = a.len().max(b.len());
    let minus1 = f.neg(f.one());
    (0..len)
        .map(|i| {
            let ai = a.get(i).cloned().unwrap_or_default();
            let bi = b.get(i).map(|x| mp_scale(f, x, minus1)).unwrap_or_default();
            mp_add(f, &ai, &bi)
        })
        .collect()
}

// ── Verifier shared by both arms ────────────────────────────────────

/// Roots in `F_p` of a monic polynomial (low degree first), with
/// multiplicity, by trial over `F_p` (validation sizes) or by
/// `gcd(g, x^p − x)` splitting (measurement sizes).
fn split_roots(g: &[u64], p: u64) -> Option<Vec<u64>> {
    let m = g.len() - 1;
    let mut roots = Vec::new();
    let mut cur = g.to_vec();
    // Trial division is exact; p ≤ 2^17 keeps it cheap next to F4.
    let mut x = 0u64;
    while x < p && roots.len() < m {
        // Horner value and synthetic quotient.
        let d = cur.len() - 1;
        let mut q = vec![0u64; d];
        let mut acc = cur[d];
        for i in (0..d).rev() {
            q[i] = acc;
            acc = (cur[i] + mulp(acc, x, p)) % p;
        }
        if acc == 0 {
            roots.push(x);
            cur = q;
            continue; // test the same x again for multiplicity
        }
        x += 1;
    }
    if roots.len() == m {
        Some(roots)
    } else {
        None
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Outcome {
    /// `g` does not split over `F_p`.
    NoSplit,
    /// Splits, but some `h(x_i)` is a non-square: the summand is on the twist.
    Twist,
    /// Splits and lifts, but no sign pattern sums to `R`.
    NoSigns,
    /// A verified relation `Σ ±P_i = R`.
    Relation(Vec<u64>),
}

fn verify(c: &Curve, r: Pt, g: &[u64], rng: &mut StdRng) -> Outcome {
    let f = c.f;
    let Some(mut xs) = split_roots(g, f.p) else {
        return Outcome::NoSplit;
    };
    xs.sort();
    let mut pts = Vec::new();
    for &x in &xs {
        let xq = f.c(x);
        match f.sqrt(c.h(xq), rng) {
            Some(y) => pts.push(Pt::A(xq, y)),
            None => return Outcome::Twist,
        }
    }
    let m = pts.len();
    for mask in 0u32..(1 << m) {
        let mut s = Pt::Inf;
        for (i, &pt) in pts.iter().enumerate() {
            let q = if mask >> i & 1 == 1 { c.neg(pt) } else { pt };
            s = c.add(s, q);
        }
        if s == r {
            return Outcome::Relation(xs);
        }
    }
    Outcome::NoSigns
}

// ── Semaev: numeric S_k, interpolation, symmetrisation ─────────────

/// `S_3(x1, x2, V)` as a polynomial in `V` (low degree first).
fn s3_in_v(f: &Field, a: Fq, b: Fq, x1: Fq, x2: Fq) -> Vec<Fq> {
    let d = f.sub(x1, x2);
    let s = f.add(x1, x2);
    let pr = f.mul(x1, x2);
    let c2 = f.sqr(d);
    let c1 = f.neg(f.scale(f.add(f.mul(s, f.add(pr, a)), f.scale(b, 2)), 2));
    let c0 = f.sub(f.sqr(f.sub(pr, a)), f.scale(f.mul(b, s), 4));
    vec![c0, c1, c2]
}

fn horner(f: &Field, c: &[Fq], x: Fq) -> Fq {
    c.iter()
        .rev()
        .fold(f.zero(), |acc, &ci| f.add(f.mul(acc, x), ci))
}

/// Determinant by cofactor expansion with minors memoised on column masks
/// (division-free, `O(k·2^k)` for `k ≤ 10`).
fn det(f: &Field, m: &[Vec<Fq>]) -> Fq {
    let k = m.len();
    let mut memo: HashMap<u32, Fq> = HashMap::new();
    fn rec(f: &Field, m: &[Vec<Fq>], row: usize, cols: u32, memo: &mut HashMap<u32, Fq>) -> Fq {
        if row == m.len() {
            return f.one();
        }
        if let Some(v) = memo.get(&cols) {
            return *v;
        }
        let mut acc = f.zero();
        let mut sign_pos = true;
        for j in 0..m.len() {
            if cols >> j & 1 == 1 {
                continue;
            }
            let e = m[row][j];
            if !f.is_zero(e) {
                let minor = rec(f, m, row + 1, cols | 1 << j, memo);
                let t = f.mul(e, minor);
                acc = if sign_pos {
                    f.add(acc, t)
                } else {
                    f.sub(acc, t)
                };
            }
            sign_pos = !sign_pos;
        }
        memo.insert(cols, acc);
        acc
    }
    let _ = k;
    rec(f, m, 0, 0, &mut memo)
}

/// Resultant of two polynomials with formal degrees `len − 1` (Sylvester).
fn resultant(f: &Field, p: &[Fq], q: &[Fq]) -> Fq {
    let (dp, dq) = (p.len() - 1, q.len() - 1);
    let k = dp + dq;
    let mut m = vec![vec![f.zero(); k]; k];
    for i in 0..dq {
        for (j, &c) in p.iter().rev().enumerate() {
            m[i][i + j] = c;
        }
    }
    for i in 0..dp {
        for (j, &c) in q.iter().rev().enumerate() {
            m[dq + i][i + j] = c;
        }
    }
    det(f, &m)
}

/// Interpolate a univariate polynomial of degree `≤ d` from its values at
/// `0, 1, …, d` (inverse Vandermonde over `F_p`, applied to `F_q` values).
fn vandermonde_inv(d: usize, p: u64) -> Vec<Vec<u64>> {
    let k = d + 1;
    let mut m: Vec<Vec<u64>> = (0..k)
        .map(|i| {
            let mut row: Vec<u64> = (0..k).map(|j| powp(i as u64, j as u64, p)).collect();
            row.extend((0..k).map(|j| (i == j) as u64));
            row
        })
        .collect();
    for col in 0..k {
        let piv = (col..k)
            .find(|&r| m[r][col] != 0)
            .expect("singular Vandermonde");
        m.swap(col, piv);
        let inv = invp(m[col][col], p);
        for v in m[col].iter_mut() {
            *v = mulp(*v, inv, p);
        }
        for r in 0..k {
            if r != col && m[r][col] != 0 {
                let fac = m[r][col];
                let pivot_row = m[col].clone();
                for (v, pv) in m[r].iter_mut().zip(&pivot_row) {
                    *v = (*v + p - mulp(fac, *pv, p)) % p;
                }
            }
        }
    }
    // V · c = values  ⇒  c = V⁻¹ · values; rows of V are evaluation points.
    m.into_iter().map(|row| row[k..].to_vec()).collect()
}

fn apply_vinv(f: &Field, vinv: &[Vec<u64>], vals: &[Fq]) -> Vec<Fq> {
    vinv.iter()
        .map(|row| {
            row.iter()
                .zip(vals)
                .fold(f.zero(), |acc, (&w, &v)| f.add(acc, f.scale(v, w)))
        })
        .collect()
}

/// `S_k(xs..., V)` as a polynomial in `V`, for `xs.len() = k − 1 ≥ 2`.
fn semaev_in_last(
    f: &Field,
    a: Fq,
    b: Fq,
    xs: &[Fq],
    vinvs: &HashMap<usize, Vec<Vec<u64>>>,
) -> Vec<Fq> {
    let k = xs.len() + 1;
    if k == 3 {
        return s3_in_v(f, a, b, xs[0], xs[1]);
    }
    // S_k(x_1..x_{k-1}, V) = Res_U(S_{k-1}(x_1..x_{k-3}, U), S_3(x_{k-2}, x_{k-1}, U))
    // only gives the value at fixed V; with V free we use
    // S_k(x_1..x_{k-1}, V) = Res_U(S_{k-1}(x_1..x_{k-2}, U), S_3(x_{k-1}, V, U)),
    // evaluated at V = 0..d and interpolated.
    let d = 1usize << (k - 2);
    let inner = semaev_in_last(f, a, b, &xs[..k - 2], vinvs);
    let vals: Vec<Fq> = (0..=d)
        .map(|v| {
            let s3 = s3_in_v(f, a, b, xs[k - 2], f.c(v as u64));
            resultant(f, &inner, &s3)
        })
        .collect();
    apply_vinv(f, &vinvs[&d], &vals)
}

fn semaev_value(f: &Field, a: Fq, b: Fq, xs: &[Fq], vinvs: &HashMap<usize, Vec<Vec<u64>>>) -> Fq {
    let last = *xs.last().unwrap();
    let poly = semaev_in_last(f, a, b, &xs[..xs.len() - 1], vinvs);
    horner(f, &poly, last)
}

/// The symmetrised summation polynomial `T(e_1..e_m, X)` with
/// `T(e(x_1..x_m), X) = S_{m+1}(x_1, …, x_m, X)`.  Keys are
/// `[b_1, …, b_m, d]` meaning `e_1^{b_1}⋯e_m^{b_m} X^d`.
struct SymSemaev {
    m: usize,
    terms: BTreeMap<Vec<u32>, Fq>,
    build_ms: f64,
    grid_points: usize,
}

fn build_sym_semaev(f: &Field, a: Fq, b: Fq, m: usize) -> SymSemaev {
    let t0 = Instant::now();
    let d = 1usize << (m - 1); // degree in each variable of S_{m+1}
    let mut vinvs = HashMap::new();
    for e in 0..=m + 1 {
        let dd = 1usize << e;
        if dd <= d.max(2) {
            vinvs.insert(dd, vandermonde_inv(dd, f.p));
        }
    }
    vinvs.insert(d, vandermonde_inv(d, f.p));
    let nv = m + 1;
    let side = d + 1;
    let total = side.pow(nv as u32);
    // Values on the grid {0..d}^{m+1}; index = Σ idx_i · side^i.
    let mut vals = vec![f.zero(); total];
    for (lin, v) in vals.iter_mut().enumerate() {
        let mut r = lin;
        let xs: Vec<Fq> = (0..nv)
            .map(|_| {
                let c = r % side;
                r /= side;
                f.c(c as u64)
            })
            .collect();
        *v = semaev_value(f, a, b, &xs, &vinvs);
    }
    // Axis-by-axis inverse Vandermonde → monomial coefficients.
    let vinv = &vinvs[&d];
    let mut stride = 1usize;
    for _axis in 0..nv {
        let mut out = vals.clone();
        for base in 0..total {
            if (base / stride) % side != 0 {
                continue;
            }
            let fiber: Vec<Fq> = (0..side).map(|j| vals[base + j * stride]).collect();
            let coeffs = apply_vinv(f, vinv, &fiber);
            for (j, c) in coeffs.into_iter().enumerate() {
                out[base + j * stride] = c;
            }
        }
        vals = out;
        stride *= side;
    }
    // Dense coefficients → sparse polynomial in (x_1..x_m, X).
    let mut s: BTreeMap<Vec<u32>, Fq> = BTreeMap::new();
    for (lin, c) in vals.iter().enumerate() {
        if f.is_zero(*c) {
            continue;
        }
        let mut r = lin;
        let e: Vec<u32> = (0..nv)
            .map(|_| {
                let v = (r % side) as u32;
                r /= side;
                v
            })
            .collect();
        s.insert(e, *c);
    }
    // Fundamental theorem: peel the lex-leading x-monomial with e-products.
    let elem: Vec<Mp> = (1..=m).map(|k| elementary(f, m, k)).collect();
    let mut eprod_cache: HashMap<Vec<u32>, Mp> = HashMap::new();
    let mut out: BTreeMap<Vec<u32>, Fq> = BTreeMap::new();
    while let Some((lead, _)) = s.last_key_value() {
        let ax: Vec<u32> = lead[..m].to_vec();
        for w in ax.windows(2) {
            assert!(
                w[0] >= w[1],
                "summation polynomial not symmetric: lead {:?}",
                ax
            );
        }
        let bexp: Vec<u32> = (0..m)
            .map(|i| ax[i] - if i + 1 < m { ax[i + 1] } else { 0 })
            .collect();
        let lo = {
            let mut k = ax.clone();
            k.push(0);
            k
        };
        let hi = {
            let mut k = ax.clone();
            k.push(u32::MAX);
            k
        };
        let cx: Vec<(u32, Fq)> = s.range(lo..=hi).map(|(k, c)| (k[m], *c)).collect();
        let ep = eprod_cache
            .entry(bexp.clone())
            .or_insert_with(|| {
                let mut acc = mp_const(f, m, f.one());
                for (i, &bi) in bexp.iter().enumerate() {
                    for _ in 0..bi {
                        acc = mp_mul(f, &acc, &elem[i]);
                    }
                }
                acc
            })
            .clone();
        for (dx, c) in cx {
            let mut key = bexp.clone();
            key.push(dx);
            out.insert(key, c);
            for (e, ec) in &ep {
                let mut k = e.clone();
                k.push(dx);
                let v = s.entry(k.clone()).or_insert(f.zero());
                *v = f.sub(*v, f.mul(c, *ec));
                if f.is_zero(*v) {
                    s.remove(&k);
                }
            }
        }
    }
    SymSemaev {
        m,
        terms: out,
        build_ms: t0.elapsed().as_secs_f64() * 1e3,
        grid_points: total,
    }
}

fn elementary(f: &Field, m: usize, k: usize) -> Mp {
    let mut r = Mp::new();
    for mask in 0u32..(1 << m) {
        if mask.count_ones() as usize == k {
            let e: Vec<u32> = (0..m).map(|i| mask >> i & 1).collect();
            r.insert(e, f.one());
        }
    }
    r
}

impl SymSemaev {
    /// Specialise `X = x_R`: a polynomial in `e_1..e_m` over `F_q`.
    fn at(&self, f: &Field, xr: Fq) -> Mp {
        let mut r = Mp::new();
        let maxd = self.terms.keys().map(|k| k[self.m]).max().unwrap_or(0);
        let pw: Vec<Fq> = (0..=maxd)
            .scan(f.one(), |acc, _| {
                let v = *acc;
                *acc = f.mul(*acc, xr);
                Some(v)
            })
            .collect();
        for (k, c) in &self.terms {
            let e = k[..self.m].to_vec();
            let v = r.entry(e).or_insert(f.zero());
            *v = f.add(*v, f.mul(*c, pw[k[self.m] as usize]));
        }
        r.retain(|_, c| !f.is_zero(*c));
        r
    }
    fn total_degree(&self) -> u32 {
        self.terms
            .keys()
            .map(|k| k[..self.m].iter().sum::<u32>())
            .max()
            .unwrap_or(0)
    }
}

// ── The two systems ─────────────────────────────────────────────────

struct System {
    polys: Vec<Poly>,
    n_vars: usize,
    build_ms: f64,
}

fn sem_system(f: &Field, t: &SymSemaev, xr: Fq) -> System {
    let t0 = Instant::now();
    let spec = t.at(f, xr);
    let polys: Vec<Poly> = (0..f.n)
        .map(|k| component(f, &spec, k))
        .filter(|p| !p.is_empty())
        .collect();
    System {
        polys,
        n_vars: t.m,
        build_ms: t0.elapsed().as_secs_f64() * 1e3,
    }
}

/// `g` (monic, low degree first) from `e` with `g = ∏(x − x_i)`.
fn g_from_e(e: &[u64], p: u64) -> Vec<u64> {
    let m = e.len();
    let mut g = vec![0u64; m + 1];
    g[m] = 1;
    for (j, &ej) in e.iter().enumerate() {
        let k = j + 1;
        g[m - k] = if k % 2 == 1 { (p - ej % p) % p } else { ej % p };
    }
    g
}

/// The Riemann–Roch basis of `L(k·O)` by pole order: `(order, is_y, i)` for
/// `x^i` (order `2i`) or `x^i·y` (order `2i + 3`).
fn rr_basis(k: usize) -> Vec<(usize, bool, usize)> {
    let mut v = Vec::new();
    for ord in 0..=k {
        if ord == 1 {
            continue;
        }
        if ord % 2 == 0 {
            v.push((ord, false, ord / 2));
        } else {
            v.push((ord, true, (ord - 3) / 2));
        }
    }
    v
}

struct RrSystem {
    sys: System,
    /// `g_0..g_{m-1}` as polynomials in the unknowns (component 0 gives `F_p` values).
    g: Vec<Mp>,
}

fn rr_system(c: &Curve, m: usize, xr: Fq, yr: Fq) -> RrSystem {
    let t0 = Instant::now();
    let f = c.f;
    let n = f.n;
    let k = m + 1;
    let basis = rr_basis(k);
    let free: Vec<(usize, bool, usize)> = basis
        .iter()
        .copied()
        .filter(|&(o, _, _)| o != 0 && o != k)
        .collect();
    assert_eq!(free.len(), m - 1);
    let nv = n * (m - 1);
    // Value of a basis function at −R = (x_R, −y_R).
    let at_neg_r = |is_y: bool, i: usize| {
        let xi = f.pow(xr, i as u128);
        if is_y {
            f.mul(xi, f.neg(yr))
        } else {
            xi
        }
    };
    let lead = *basis.iter().find(|&&(o, _, _)| o == k).unwrap();
    // c_j = Σ_s θ^s z_{j,s}
    let coef: Vec<Mp> = (0..free.len())
        .map(|j| {
            let mut mp = Mp::new();
            for s in 0..n {
                let mut e = vec![0u32; nv];
                e[j * n + s] = 1;
                mp.insert(e, f.theta_pow(s));
            }
            mp
        })
        .collect();
    // c_0 = −(lead(−R) + Σ c_j φ_j(−R))
    let mut c0 = mp_const(f, nv, f.neg(at_neg_r(lead.1, lead.2)));
    for (j, &(_, is_y, i)) in free.iter().enumerate() {
        c0 = mp_add(f, &c0, &mp_scale(f, &coef[j], f.neg(at_neg_r(is_y, i))));
    }
    let deg_a = k / 2;
    let deg_b = if k >= 3 { (k - 3) / 2 } else { 0 };
    let mut a_poly: XPoly = vec![Mp::new(); deg_a + 1];
    let mut b_poly: XPoly = vec![Mp::new(); deg_b + 1];
    let mut place = |is_y: bool, i: usize, mp: Mp| {
        let slot = if is_y { &mut b_poly[i] } else { &mut a_poly[i] };
        *slot = mp_add(f, slot, &mp);
    };
    place(lead.1, lead.2, mp_const(f, nv, f.one()));
    place(false, 0, c0);
    for (j, &(_, is_y, i)) in free.iter().enumerate() {
        place(is_y, i, coef[j].clone());
    }
    let h: XPoly = vec![
        mp_const(f, nv, c.b),
        mp_const(f, nv, c.a),
        Mp::new(),
        mp_const(f, nv, f.one()),
    ];
    let a2 = xp_mul(f, &a_poly, &a_poly);
    let b2h = xp_mul(f, &xp_mul(f, &b_poly, &b_poly), &h);
    let mut norm = xp_sub(f, &a2, &b2h);
    while norm.last().is_some_and(|c| c.is_empty()) {
        norm.pop();
    }
    assert_eq!(
        norm.len(),
        k + 1,
        "norm has degree {} not {}",
        norm.len() - 1,
        k
    );
    let lc = norm[k]
        .get(&vec![0u32; nv])
        .copied()
        .expect("constant leading coefficient");
    assert_eq!(norm[k].len(), 1);
    let lci = f.inv(lc);
    let norm: XPoly = norm.iter().map(|c| mp_scale(f, c, lci)).collect();
    // Synthetic division by (x − x_R).
    let mut q: XPoly = vec![Mp::new(); k];
    q[k - 1] = norm[k].clone();
    for j in (1..k).rev() {
        q[j - 1] = mp_add(f, &norm[j], &mp_scale(f, &q[j], xr));
    }
    let rem = mp_add(f, &norm[0], &mp_scale(f, &q[0], xr));
    assert!(rem.is_empty(), "x_R is not a root of the norm");
    let g: Vec<Mp> = q[..m].to_vec();
    let mut polys = Vec::new();
    for gi in &g {
        for s in 1..n {
            let pc = component(f, gi, s);
            if !pc.is_empty() {
                polys.push(pc);
            }
        }
    }
    RrSystem {
        sys: System {
            polys,
            n_vars: nv,
            build_ms: t0.elapsed().as_secs_f64() * 1e3,
        },
        g,
    }
}

// ── Solving and scoring ─────────────────────────────────────────────

#[derive(Debug)]
struct ArmResult {
    status: &'static str,
    algebraic: usize,
    outcomes: BTreeMap<String, usize>,
    relations: BTreeSet<Vec<u64>>,
    report: Option<ZeroDimStats>,
    max_input_degree: u32,
    n_eq: usize,
    n_vars: usize,
    build_ms: f64,
}

fn max_degree(polys: &[Poly]) -> u32 {
    polys
        .iter()
        .flat_map(|p| p.iter().map(|(e, _)| e.iter().sum::<u32>()))
        .max()
        .unwrap_or(0)
}

/// One instrument for both arms: grevlex F4 to a Gröbner basis, then
/// multiplication matrices, a Krylov minimal polynomial and eigenvectors
/// (`jv_cover::solve_zero_dim`, the Joux–Vitse reproduction's solver).  The
/// deadline is set here, per solve.
fn run_arm(
    c: &Curve,
    r: Pt,
    sys: &System,
    to_g: &dyn Fn(&[u64]) -> Vec<u64>,
    budget: (u64, u32),
    rng: &mut StdRng,
) -> ArmResult {
    let opts =
        F4Options::new(Ordering::Grevlex, budget.1).with_budget(Duration::from_secs(budget.0));
    let (sols, st) = solve_zero_dim(&sys.polys, sys.n_vars, c.f.p, &opts, rng);
    let verdict = if st.timed_out || st.incomplete {
        Verdict::Undetermined
    } else if st.inconsistent || sols.is_empty() {
        Verdict::Inconsistent
    } else {
        Verdict::Solutions(sols)
    };
    let rep = st;
    let mut res = ArmResult {
        status: "",
        algebraic: 0,
        outcomes: BTreeMap::new(),
        relations: BTreeSet::new(),
        report: None,
        max_input_degree: max_degree(&sys.polys),
        n_eq: sys.polys.len(),
        n_vars: sys.n_vars,
        build_ms: sys.build_ms,
    };
    match &verdict {
        Verdict::Inconsistent => res.status = "inconsistent",
        Verdict::Undetermined => res.status = "undetermined",
        Verdict::Solutions(sols) => {
            res.status = "solved";
            res.algebraic = sols.len();
            for s in sols {
                let g = to_g(s);
                let o = verify(c, r, &g, rng);
                let key = match &o {
                    Outcome::NoSplit => "nosplit",
                    Outcome::Twist => "twist",
                    Outcome::NoSigns => "nosigns",
                    Outcome::Relation(_) => "relation",
                };
                *res.outcomes.entry(key.to_string()).or_default() += 1;
                if let Outcome::Relation(xs) = o {
                    res.relations.insert(xs);
                }
            }
        }
    }
    res.report = Some(rep);
    res
}

fn sem_to_g(p: u64) -> impl Fn(&[u64]) -> Vec<u64> {
    move |e: &[u64]| g_from_e(e, p)
}

fn rr_to_g<'a>(f: &'a Field, g: &'a [Mp]) -> impl Fn(&[u64]) -> Vec<u64> + 'a {
    move |z: &[u64]| {
        let mut out: Vec<u64> = g
            .iter()
            .map(|gi| {
                let v = mp_eval_fp(f, gi, z);
                debug_assert!(f.in_fp(v));
                v.0[0]
            })
            .collect();
        out.push(1);
        out
    }
}

// ── Targets ─────────────────────────────────────────────────────────

fn random_fb_point(c: &Curve, rng: &mut StdRng) -> Pt {
    let f = c.f;
    loop {
        let x = f.c(rng.gen_range(0..f.p));
        if let Some(y) = f.sqrt(c.h(x), rng) {
            let y = if rng.gen::<bool>() { y } else { f.neg(y) };
            return Pt::A(x, y);
        }
    }
}

fn random_point(c: &Curve, rng: &mut StdRng) -> Pt {
    let f = c.f;
    loop {
        let x = f.random(rng);
        if let Some(y) = f.sqrt(c.h(x), rng) {
            return Pt::A(x, y);
        }
    }
}

/// A target with `x_R ∉ F_p`, `y_R ≠ 0`.
fn draw_target(c: &Curve, m: usize, planted: bool, rng: &mut StdRng) -> (Pt, Option<Vec<u64>>) {
    loop {
        let (r, plant) = if planted {
            let pts: Vec<Pt> = (0..m).map(|_| random_fb_point(c, rng)).collect();
            let mut s = Pt::Inf;
            for &p in &pts {
                s = c.add(s, p);
            }
            let mut xs: Vec<u64> = pts
                .iter()
                .map(|p| match p {
                    Pt::A(x, _) => x.0[0],
                    Pt::Inf => unreachable!(),
                })
                .collect();
            xs.sort();
            (s, Some(xs))
        } else {
            (random_point(c, rng), None)
        };
        if let Pt::A(x, y) = r {
            if !c.f.in_fp(x) && !c.f.is_zero(y) {
                return (r, plant);
            }
        }
    }
}

fn random_curve<'a>(f: &'a Field, rng: &mut StdRng) -> Curve<'a> {
    loop {
        let a = f.random(rng);
        let b = f.random(rng);
        if f.in_fp(a) && f.in_fp(b) {
            continue;
        }
        let disc = f.add(f.scale(f.mul(f.sqr(a), a), 4), f.scale(f.sqr(b), 27));
        if !f.is_zero(disc) {
            return Curve { f, a, b };
        }
    }
}

// ── Brute force (validation only) ───────────────────────────────────

/// Positive control for the `rr` construction without F4: interpolate the
/// function `f ∈ L((m+1)O)` through a planted decomposition and check that its
/// Weil coordinates satisfy every equation of `rr_system`.  `None` when the
/// planted points are not distinct (the interpolation is then singular).
fn rr_plant_check(c: &Curve, m: usize, r: Pt, xs: &[u64], rng: &mut StdRng) -> Option<bool> {
    let f = c.f;
    let Pt::A(xr, yr) = r else { return None };
    let lifted: Vec<(Fq, Fq)> = xs
        .iter()
        .map(|&x| (f.c(x), f.sqrt(c.h(f.c(x)), rng).unwrap()))
        .collect();
    let mut pts = None;
    for mask in 0u32..(1 << m) {
        let ps: Vec<Pt> = lifted
            .iter()
            .enumerate()
            .map(|(i, &(x, y))| Pt::A(x, if mask >> i & 1 == 1 { f.neg(y) } else { y }))
            .collect();
        let s = ps.iter().fold(Pt::Inf, |acc, &p| c.add(acc, p));
        if s == r {
            pts = Some(ps);
            break;
        }
    }
    let pts = pts?;
    let k = m + 1;
    let basis = rr_basis(k);
    let lead = *basis.iter().find(|&&(o, _, _)| o == k).unwrap();
    let rest: Vec<(usize, bool, usize)> =
        basis.iter().copied().filter(|&(o, _, _)| o != k).collect();
    let ev = |(_, is_y, i): (usize, bool, usize), x: Fq, y: Fq| {
        let xi = f.pow(x, i as u128);
        if is_y {
            f.mul(xi, y)
        } else {
            xi
        }
    };
    // rows: Σ_j coef_j φ_j(P) = −lead(P)
    let mut a: Vec<Vec<Fq>> = pts
        .iter()
        .map(|p| {
            let Pt::A(x, y) = *p else { unreachable!() };
            let mut row: Vec<Fq> = rest.iter().map(|&b| ev(b, x, y)).collect();
            row.push(f.neg(ev(lead, x, y)));
            row
        })
        .collect();
    let w = rest.len();
    for col in 0..w {
        let piv = (col..a.len()).find(|&i| !f.is_zero(a[i][col]))?;
        a.swap(col, piv);
        let inv = f.inv(a[col][col]);
        for v in a[col].iter_mut() {
            *v = f.mul(*v, inv);
        }
        for i in 0..a.len() {
            if i != col && !f.is_zero(a[i][col]) {
                let fac = a[i][col];
                let pr = a[col].clone();
                for (v, pv) in a[i].iter_mut().zip(&pr) {
                    *v = f.sub(*v, f.mul(fac, *pv));
                }
            }
        }
    }
    let sol: Vec<Fq> = (0..w).map(|j| a[j][w]).collect();
    // f(−R) = 0 must hold for the interpolated f.
    let at = |x: Fq, y: Fq| {
        rest.iter()
            .zip(&sol)
            .fold(ev(lead, x, y), |acc, (&b, &cj)| {
                f.add(acc, f.mul(cj, ev(b, x, y)))
            })
    };
    if !f.is_zero(at(xr, f.neg(yr))) {
        return Some(false);
    }
    // Free coefficients in rr_system's order: every basis function but lead and 1.
    let free: Vec<Fq> = rest
        .iter()
        .zip(&sol)
        .filter(|(&(o, _, _), _)| o != 0)
        .map(|(_, &cj)| cj)
        .collect();
    let z: Vec<u64> = free.iter().flat_map(|cj| cj.0[..f.n].to_vec()).collect();
    let rs = rr_system(c, m, xr, yr);
    let holds = rs.sys.polys.iter().all(|p| f4_fp::eval(p, &z, f.p) == 0);
    // Negative control: one coordinate moved by 1 must break some equation.
    let mut z1 = z.clone();
    z1[0] = (z1[0] + 1) % f.p;
    let broken = rs.sys.polys.iter().any(|p| f4_fp::eval(p, &z1, f.p) != 0);
    Some(holds && broken)
}

fn brute_force(c: &Curve, r: Pt, m: usize, rng: &mut StdRng) -> BTreeSet<Vec<u64>> {
    let f = c.f;
    let fb: Vec<u64> = (0..f.p).filter(|&x| f.is_square(c.h(f.c(x)))).collect();
    let mut out = BTreeSet::new();
    let mut idx = vec![0usize; m];
    loop {
        let xs: Vec<u64> = idx.iter().map(|&i| fb[i]).collect();
        // g = ∏(x − x_i)
        let mut g = vec![1u64];
        for &x in &xs {
            let mut ng = vec![0u64; g.len() + 1];
            for (i, &gi) in g.iter().enumerate() {
                ng[i + 1] = (ng[i + 1] + gi) % f.p;
                ng[i] = (ng[i] + f.p - mulp(gi, x, f.p)) % f.p;
            }
            g = ng;
        }
        if let Outcome::Relation(v) = verify(c, r, &g, rng) {
            out.insert(v);
        }
        // next non-decreasing index tuple
        let mut pos = m;
        loop {
            if pos == 0 {
                return out;
            }
            pos -= 1;
            if idx[pos] + 1 < fb.len() {
                idx[pos] += 1;
                for j in pos + 1..m {
                    idx[j] = idx[pos];
                }
                break;
            }
        }
    }
}

// ── CLI ─────────────────────────────────────────────────────────────

fn arg<T: std::str::FromStr>(args: &[String], key: &str, default: T) -> T {
    args.iter()
        .position(|a| a == key)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn arm_json(name: &str, a: &ArmResult) -> serde_json::Value {
    let rep = a.report.as_ref().unwrap();
    serde_json::json!({
        "arm": name,
        "status": a.status,
        "n_vars": a.n_vars,
        "n_eq": a.n_eq,
        "input_degree": a.max_input_degree,
        "build_ms": a.build_ms,
        "solving_degree": rep.solving_degree,
        "degree_reached": rep.degree_reached,
        "max_rows": rep.max_rows,
        "max_cols": rep.max_cols,
        "f4_muls": rep.f4_muls,
        "lin_muls": rep.lin_muls,
        "delta": rep.delta,
        "f4_ms": rep.f4_ms,
        "lin_ms": rep.lin_ms,
        "timed_out": rep.timed_out,
        "incomplete": rep.incomplete,
        "nonseparating": rep.nonseparating,
        "algebraic_solutions": a.algebraic,
        "outcomes": a.outcomes,
        "relations": a.relations.iter().collect::<Vec<_>>(),
    })
}

fn validate(p: u64, seed: u64, targets: usize, cells: &[(usize, usize)], arms: &[&str]) {
    let mut ok = true;
    for &(n, m) in cells {
        let mut rng = StdRng::seed_from_u64(seed ^ ((n as u64) << 8 | m as u64));
        let f = Field::new(p, n, &mut rng);
        let c = random_curve(&f, &mut rng);
        let t = build_sym_semaev(&f, c.a, c.b, m);
        // Identity check: T(e(x), X) = S_{m+1}(x, X) at random points of F_q.
        let mut vinvs = HashMap::new();
        for e in 1..=m {
            let dd = 1usize << e;
            vinvs.insert(dd, vandermonde_inv(dd, p));
        }
        for _ in 0..4 {
            let xs: Vec<Fq> = (0..=m).map(|_| f.random(&mut rng)).collect();
            let lhs = semaev_value(&f, c.a, c.b, &xs, &vinvs);
            // e_k(x_1..x_m) over F_q
            let mut ev = vec![f.one()];
            for &x in &xs[..m] {
                let mut ne = vec![f.zero(); ev.len() + 1];
                for (i, &v) in ev.iter().enumerate() {
                    ne[i] = f.add(ne[i], v);
                    ne[i + 1] = f.add(ne[i + 1], f.mul(v, x));
                }
                ev = ne;
            }
            let mut rhs = f.zero();
            for (k, cf) in &t.terms {
                let mut term = *cf;
                for i in 0..m {
                    term = f.mul(term, f.pow(ev[i + 1], k[i] as u128));
                }
                term = f.mul(term, f.pow(xs[m], k[m] as u128));
                rhs = f.add(rhs, term);
            }
            if lhs != rhs {
                println!("FAIL n={n} m={m}: symmetrisation identity");
                ok = false;
            }
        }
        let o = (600u64, 60u32);
        for ti in 0..targets {
            let planted = ti % 2 == 0;
            let (r, plant) = draw_target(&c, m, planted, &mut rng);
            let Pt::A(xr, yr) = r else { unreachable!() };
            let truth = brute_force(&c, r, m, &mut rng);
            let planted_found = plant.as_ref().map(|x| truth.contains(x));
            let mut good = planted_found != Some(false);
            let mut line = format!("n={n} m={m} t={ti} planted={planted} truth={}", truth.len());
            for arm in arms {
                let a = match *arm {
                    "sem" => {
                        let ss = sem_system(&f, &t, xr);
                        run_arm(&c, r, &ss, &sem_to_g(p), o, &mut rng)
                    }
                    _ => {
                        let rs = rr_system(&c, m, xr, yr);
                        let to_g = rr_to_g(&f, &rs.g);
                        run_arm(&c, r, &rs.sys, &to_g, o, &mut rng)
                    }
                };
                good &= a.relations == truth;
                line += &format!(
                    " {arm}={}({},{:?})",
                    a.relations.len(),
                    a.status,
                    a.outcomes
                );
            }
            ok &= good;
            println!("{line} {}", if good { "ok" } else { "MISMATCH" });
        }
        println!(
            "n={n} m={m}: |T| = {} terms, deg_e = {}, build {:.1} ms",
            t.terms.len(),
            t.total_degree(),
            t.build_ms
        );
    }
    println!(
        "{}",
        if ok {
            "VALIDATION PASSED"
        } else {
            "VALIDATION FAILED"
        }
    );
    if !ok {
        std::process::exit(1);
    }
}

/// The 12- and 15-unknown `rr` systems overflow the default 8 MB main stack
/// inside F4; the whole run goes on one thread with a 1 GB stack.
fn main() {
    std::thread::Builder::new()
        .stack_size(1 << 30)
        .spawn(real_main)
        .expect("spawn")
        .join()
        .expect("run");
}

fn real_main() {
    let args: Vec<String> = std::env::args().collect();
    let p: u64 = arg(&args, "--p", 65521);
    let seed: u64 = arg(&args, "--seed", 20261016);
    if args.iter().any(|a| a == "--rr-plant-check") {
        let targets: usize = arg(&args, "--targets", 8);
        let mut ok = true;
        for &(n, m) in &[(2usize, 2usize), (3, 2), (3, 3), (4, 3), (4, 4), (5, 4)] {
            let mut rng = StdRng::seed_from_u64(seed ^ ((n as u64) << 8 | m as u64));
            let f = Field::new(p, n, &mut rng);
            let c = random_curve(&f, &mut rng);
            let (mut pass, mut skip) = (0, 0);
            for _ in 0..targets {
                let (r, plant) = draw_target(&c, m, true, &mut rng);
                match rr_plant_check(&c, m, r, plant.as_ref().unwrap(), &mut rng) {
                    Some(true) => pass += 1,
                    Some(false) => ok = false,
                    None => skip += 1,
                }
            }
            println!("n={n} m={m} p={p}: planted f satisfies every rr equation on {pass}/{targets} targets ({skip} skipped: repeated summands)");
        }
        println!(
            "{}",
            if ok {
                "RR PLANT CHECK PASSED"
            } else {
                "RR PLANT CHECK FAILED"
            }
        );
        if !ok {
            std::process::exit(1);
        }
        return;
    }
    if args.iter().any(|a| a == "--validate") {
        let targets: usize = arg(&args, "--targets", 6);
        let cells: Vec<(usize, usize)> =
            match arg::<String>(&args, "--cells", String::new()).as_str() {
                "" => vec![(2, 2), (3, 2), (3, 3), (4, 3), (4, 4)],
                s => s
                    .split(';')
                    .map(|c| {
                        let v: Vec<usize> = c.split(',').map(|x| x.parse().unwrap()).collect();
                        (v[0], v[1])
                    })
                    .collect(),
            };
        let arms: String = arg(&args, "--arms", "sem,rr".to_string());
        let arms: Vec<&str> = arms.split(',').collect();
        validate(p, seed, targets, &cells, &arms);
        return;
    }
    let n: usize = arg(&args, "--n", 3);
    let m: usize = arg(&args, "--m", 3);
    let targets: usize = arg(&args, "--targets", 4);
    let kind: String = arg(&args, "--kind", "planted".to_string());
    let budget_s: u64 = arg(&args, "--budget-s", 300);
    let dmax: u32 = arg(&args, "--dmax", 60);
    let arms: String = arg(&args, "--arms", "sem,rr".to_string());
    let out_path: String = arg(&args, "--out", String::new());
    // `--only K`: draw targets 0..=K from the same stream, solve only K (one
    // process per target, so an out-of-memory kill loses one row, not a cell).
    let only: Option<usize> = args
        .iter()
        .position(|a| a == "--only")
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok());
    assert!(n >= 2 && n <= NMAX && m >= 2 && m <= 4 && n * (m - 1) <= 16);
    let mut out: Box<dyn Write> = if out_path.is_empty() {
        Box::new(std::io::stdout())
    } else {
        Box::new(
            std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(&out_path)
                .unwrap(),
        )
    };
    let mut rng = StdRng::seed_from_u64(seed ^ ((n as u64) << 8 | m as u64));
    let f = Field::new(p, n, &mut rng);
    let c = random_curve(&f, &mut rng);
    let t = build_sym_semaev(&f, c.a, c.b, m);
    let setup = serde_json::json!({
        "row": "setup", "p": p, "n": n, "m": m, "seed": seed,
        "modulus": f.mu, "a": c.a.0[..n], "b": c.b.0[..n],
        "sem_terms": t.terms.len(), "sem_total_degree": t.total_degree(),
        "sem_build_ms": t.build_ms, "sem_grid_points": t.grid_points,
        "budget_s": budget_s, "dmax": dmax,
    });
    writeln!(out, "{}", setup).unwrap();
    let o = (budget_s, dmax);
    let planted = kind == "planted";
    // Targets are drawn from their own stream so arms and budgets never shift them.
    let mut trng = StdRng::seed_from_u64(
        seed ^ 0x7a72_6765 ^ ((n as u64) << 8 | m as u64) ^ (planted as u64) << 20,
    );
    for ti in 0..targets {
        let (r, plant) = draw_target(&c, m, planted, &mut trng);
        if only.is_some_and(|k| k != ti) {
            continue;
        }
        let Pt::A(xr, yr) = r else { unreachable!() };
        let mut row = serde_json::json!({
            "row": "target", "n": n, "m": m, "p": p, "seed": seed, "kind": kind, "target": ti,
            "x_r": xr.0[..n], "planted": plant,
        });
        let mut arms_out = Vec::new();
        for arm in arms.split(',') {
            let a = match arm {
                "sem" => {
                    let ss = sem_system(&f, &t, xr);
                    run_arm(&c, r, &ss, &sem_to_g(p), o, &mut rng)
                }
                "rr" => {
                    let rs = rr_system(&c, m, xr, yr);
                    let to_g = rr_to_g(&f, &rs.g);
                    run_arm(&c, r, &rs.sys, &to_g, o, &mut rng)
                }
                other => panic!("unknown arm {other}"),
            };
            let mut j = arm_json(arm, &a);
            if let Some(px) = &plant {
                j["planted_found"] = serde_json::json!(a.relations.contains(px));
            }
            arms_out.push(j);
        }
        row["arms"] = serde_json::json!(arms_out);
        writeln!(out, "{}", row).unwrap();
        out.flush().unwrap();
    }
}
