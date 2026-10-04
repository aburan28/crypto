//! **The sieving variant of the cover route** ([JV12] §3.2) on the genus-3
//! GHS cover `H / F_q`, `q = p²`, with the factor base of `F_p`-abscissae.
//!
//! Companion to `research/notes/index-calculus/RESEARCH_COVER_DECOMPOSITION_LEDGER.md`
//! §11, which registered the construction and the predictions before this
//! file existed.  The route's instance, cover, factor base, Jacobian and the
//! six-point descent test are [`super::jv_cover`]'s; this module adds the
//! relation search among factor-base points alone.
//!
//! A function `f = A(x) + B(x)·y ∈ L(m·∞)` on `H : y² = h(x)` has `m` zeros,
//! with abscissae the roots of `F = A² − B²h`.  Writing `A = A₀ + tA₁` over
//! `F_q = F_p(t)`, `t² = ω`, the condition `F ∈ F_p[x]` is the one identity
//! `2A₀A₁ = Im(B²h)`; for a fixed `B` and a monic divisor `A₀ | Im(B²h)` the
//! one-parameter family `A = A₀ + t·s·A₁′`, `B ↦ √s·B` satisfies it for every
//! `s ∈ F_p`, and along it
//!
//! ```text
//! F(x, s) = A₀(x)² + ω s² A₁′(x)² − s·Re(B²h)(x)      (quadratic in s).
//! ```
//!
//! The sieve walks `x` over the factor base's abscissae, computes the two
//! roots in `s` (the discriminant is `N(B²h)(x)`, a square on every base
//! abscissa) and increments a counter at each; a counter at `m` is a
//! relation `Σ_i (Q_i) − m·∞ ∼ 0` among `m` factor-base points, verified in
//! the Jacobian before it is used.  Everything is counted in `F_p`
//! multiplications (the ledger's unit), and the additions and table lookups
//! of the sieve are counted beside them (`S⁺`).

use std::cell::Cell;
use std::collections::HashMap;
use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::Serialize;

use super::f4_fp;
use super::gaudry_cubic::{am, mm, sm, wiedemann_u64, SparseRel};
use super::jv_cover::{
    factor_base, generate_spec, nagao_decompose, opts_for, peval, pmul, rho_e, unit_costs,
    verify_dec, BaseEl, Ctx, Div, Fld, Fq, RhoRunE, Spec, E2,
};
use super::residual_walk::inv_mod;

/// The harness's price of one `F_p` inversion, as in `jv_cover`.
const INV_COST: u64 = 16;

// ── F_p[x], counted ──────────────────────────────────────────────────────

/// Polynomials over `F_p`, coefficients low to high with no trailing zero,
/// every multiplication counted (an inversion as [`INV_COST`]).
pub struct FpRing {
    pub p: u64,
    muls: Cell<u64>,
    adds: Cell<u64>,
}

pub type FpPoly = Vec<u64>;

fn trim(mut v: FpPoly) -> FpPoly {
    while v.last() == Some(&0) {
        v.pop();
    }
    v
}

pub fn fdeg(a: &[u64]) -> isize {
    a.len() as isize - 1
}

impl FpRing {
    pub fn new(p: u64) -> FpRing {
        FpRing {
            p,
            muls: Cell::new(0),
            adds: Cell::new(0),
        }
    }
    fn count(&self, k: u64) {
        self.muls.set(self.muls.get() + k);
    }
    pub fn muls(&self) -> u64 {
        self.muls.get()
    }
    /// Additions of the evaluation-based root search ([`FpRing::roots_by_table`]).
    pub fn adds(&self) -> u64 {
        self.adds.get()
    }
    pub fn reset(&self) {
        self.muls.set(0);
        self.adds.set(0);
    }
    pub fn inv(&self, a: u64) -> u64 {
        self.count(INV_COST);
        inv_mod(a, self.p)
    }
    pub fn add(&self, a: &[u64], b: &[u64]) -> FpPoly {
        let n = a.len().max(b.len());
        trim(
            (0..n)
                .map(|i| {
                    am(
                        a.get(i).copied().unwrap_or(0),
                        b.get(i).copied().unwrap_or(0),
                        self.p,
                    )
                })
                .collect(),
        )
    }
    pub fn sub(&self, a: &[u64], b: &[u64]) -> FpPoly {
        let n = a.len().max(b.len());
        trim(
            (0..n)
                .map(|i| {
                    sm(
                        a.get(i).copied().unwrap_or(0),
                        b.get(i).copied().unwrap_or(0),
                        self.p,
                    )
                })
                .collect(),
        )
    }
    pub fn scale(&self, a: &[u64], k: u64) -> FpPoly {
        self.count(a.len() as u64);
        trim(a.iter().map(|&c| mm(c, k, self.p)).collect())
    }
    pub fn mul(&self, a: &[u64], b: &[u64]) -> FpPoly {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        self.count((a.len() * b.len()) as u64);
        let mut out = vec![0u64; a.len() + b.len() - 1];
        for (i, &x) in a.iter().enumerate() {
            if x == 0 {
                continue;
            }
            for (j, &y) in b.iter().enumerate() {
                out[i + j] = am(out[i + j], mm(x, y, self.p), self.p);
            }
        }
        trim(out)
    }
    pub fn monic(&self, a: &[u64]) -> FpPoly {
        match a.last() {
            None => Vec::new(),
            Some(&1) => a.to_vec(),
            Some(&l) => {
                let li = self.inv(l);
                self.scale(a, li)
            }
        }
    }
    /// `(q, r)` with `a = q·b + r`, `deg r < deg b`.
    pub fn divrem(&self, a: &[u64], b: &[u64]) -> (FpPoly, FpPoly) {
        let db = fdeg(b);
        assert!(db >= 0, "division by zero polynomial");
        let mut r = a.to_vec();
        if fdeg(&r) < db {
            return (Vec::new(), r);
        }
        let li = if *b.last().unwrap() == 1 {
            1
        } else {
            self.inv(*b.last().unwrap())
        };
        let mut q = vec![0u64; r.len() - b.len() + 1];
        while fdeg(&r) >= db {
            let shift = r.len() - b.len();
            let coef = if li == 1 {
                *r.last().unwrap()
            } else {
                self.count(1);
                mm(*r.last().unwrap(), li, self.p)
            };
            q[shift] = coef;
            self.count(b.len() as u64);
            for (i, &bc) in b.iter().enumerate() {
                r[shift + i] = sm(r[shift + i], mm(coef, bc, self.p), self.p);
            }
            r = trim(r);
        }
        (trim(q), r)
    }
    pub fn rem(&self, a: &[u64], b: &[u64]) -> FpPoly {
        self.divrem(a, b).1
    }
    pub fn exact_div(&self, a: &[u64], b: &[u64]) -> FpPoly {
        let (q, r) = self.divrem(a, b);
        debug_assert!(r.is_empty(), "inexact division");
        q
    }
    /// Monic gcd.
    pub fn gcd(&self, a: &[u64], b: &[u64]) -> FpPoly {
        let (mut x, mut y) = (a.to_vec(), b.to_vec());
        while !y.is_empty() {
            let r = self.rem(&x, &y);
            x = y;
            y = r;
        }
        self.monic(&x)
    }
    pub fn deriv(&self, a: &[u64]) -> FpPoly {
        if a.len() < 2 {
            return Vec::new();
        }
        self.count(a.len() as u64 - 1);
        trim(
            a.iter()
                .enumerate()
                .skip(1)
                .map(|(i, &c)| mm(c, (i as u64) % self.p, self.p))
                .collect(),
        )
    }
    pub fn eval(&self, a: &[u64], x: u64) -> u64 {
        if a.is_empty() {
            return 0;
        }
        self.count(a.len() as u64 - 1);
        a.iter()
            .rev()
            .fold(0u64, |acc, &c| am(mm(acc, x, self.p), c, self.p))
    }
    pub fn mulmod(&self, a: &[u64], b: &[u64], m: &[u64]) -> FpPoly {
        self.rem(&self.mul(a, b), m)
    }
    pub fn powmod(&self, a: &[u64], mut e: u128, m: &[u64]) -> FpPoly {
        let mut base = self.rem(a, m);
        let mut acc = vec![1u64];
        while e > 0 {
            if e & 1 == 1 {
                acc = self.mulmod(&acc, &base, m);
            }
            e >>= 1;
            if e > 0 {
                base = self.mulmod(&base, &base, m);
            }
        }
        acc
    }
    /// `f(g) mod m` by Horner.
    pub fn compose_mod(&self, f: &[u64], g: &[u64], m: &[u64]) -> FpPoly {
        let mut acc: FpPoly = Vec::new();
        for &c in f.iter().rev() {
            acc = self.mulmod(&acc, g, m);
            acc = self.add(&acc, &[c]);
        }
        acc
    }
    /// Yun's squarefree decomposition of a monic `a` (`deg a < p`): pairs
    /// `(a_i, i)` with `a = Π a_i^i`, each `a_i` monic and squarefree.
    pub fn squarefree(&self, a: &[u64]) -> Vec<(FpPoly, usize)> {
        let a = self.monic(a);
        let mut out = Vec::new();
        if fdeg(&a) <= 0 {
            return out;
        }
        let da = self.deriv(&a);
        let a1 = self.gcd(&a, &da);
        let mut b = self.exact_div(&a, &a1);
        let c = self.exact_div(&da, &a1);
        let mut d = self.sub(&c, &self.deriv(&b));
        let mut i = 1;
        while fdeg(&b) > 0 {
            let g = self.gcd(&b, &d);
            if fdeg(&g) > 0 {
                out.push((g.clone(), i));
            }
            b = self.exact_div(&b, &g);
            let cc = self.exact_div(&d, &g);
            d = self.sub(&cc, &self.deriv(&b));
            i += 1;
        }
        out
    }
    /// `x^{p^k} mod z` for `k = 1, 2, …`: `x^p` by repeated squaring, then
    /// each next power by the Frobenius matrix (`(Σ c_i x^i)^p = Σ c_i x^{ip}`,
    /// `n²` multiplications a step instead of a composition).
    fn frobenius_powers(&self, z: &[u64], upto: usize) -> Vec<FpPoly> {
        let n = fdeg(z) as usize;
        let x = vec![0u64, 1];
        let xp = self.powmod(&x, self.p as u128, z);
        let mut out = vec![xp.clone()];
        if upto < 2 || n < 2 {
            return out;
        }
        // column i: x^{ip} mod z
        let mut cols: Vec<FpPoly> = Vec::with_capacity(n);
        let mut cur = vec![1u64];
        for _ in 0..n {
            cols.push(cur.clone());
            cur = self.mulmod(&cur, &xp, z);
        }
        for _ in 1..upto {
            let prev = out.last().unwrap();
            let mut next = vec![0u64; n];
            for (i, &c) in prev.iter().enumerate() {
                if c == 0 {
                    continue;
                }
                for (j, &v) in cols[i].iter().enumerate() {
                    next[j] = am(next[j], mm(c, v, self.p), self.p);
                }
            }
            self.count((n * n) as u64);
            out.push(trim(next));
        }
        out
    }
    /// Distinct-degree factorisation of a monic squarefree `z`: pairs
    /// `(g_d, d)`, `g_d` the product of the irreducible factors of degree `d`.
    pub fn ddf(&self, z: &[u64]) -> Vec<(FpPoly, usize)> {
        self.ddf_from(z, 1)
    }
    /// [`FpRing::ddf`] for a `z` known to have no irreducible factor of
    /// degree below `from`.
    pub fn ddf_from(&self, z: &[u64], from: usize) -> Vec<(FpPoly, usize)> {
        let z0 = self.monic(z);
        let mut z = z0.clone();
        let mut out = Vec::new();
        if fdeg(&z) <= 0 {
            return out;
        }
        let x = vec![0u64, 1];
        let n = fdeg(&z) as usize;
        let powers = self.frobenius_powers(&z0, n / 2);
        for (d, xpd) in powers.iter().enumerate().map(|(i, v)| (i + 1, v)) {
            if d < from {
                continue;
            }
            if fdeg(&z) < 2 * d as isize {
                break;
            }
            let g = self.gcd(&z, &self.sub(&self.rem(xpd, &z), &x));
            if fdeg(&g) > 0 {
                out.push((g.clone(), d));
                z = self.exact_div(&z, &g);
                if fdeg(&z) <= 0 {
                    break;
                }
            }
        }
        if fdeg(&z) > 0 {
            out.push((z.clone(), fdeg(&z) as usize));
        }
        out
    }
    /// Cantor–Zassenhaus: the irreducible factors of a monic `g` that is a
    /// product of distinct irreducibles of degree `d`.
    pub fn edf(&self, g: &[u64], d: usize, rng: &mut StdRng) -> Vec<FpPoly> {
        let g = self.monic(g);
        let n = fdeg(&g);
        if n <= 0 {
            return Vec::new();
        }
        if n as usize == d {
            return vec![g];
        }
        let e = ((self.p as u128).pow(d as u32) - 1) / 2;
        loop {
            let r: FpPoly = trim((0..n as usize).map(|_| rng.gen_range(0..self.p)).collect());
            if fdeg(&r) < 1 {
                continue;
            }
            let w = self.sub(&self.powmod(&r, e, &g), &[1]);
            let h = self.gcd(&g, &w);
            if fdeg(&h) > 0 && fdeg(&h) < n {
                let other = self.exact_div(&g, &h);
                let mut out = self.edf(&h, d, rng);
                out.extend(self.edf(&other, d, rng));
                return out;
            }
        }
    }
    /// The distinct roots in `F_p` of `a` by evaluating it at every `x`
    /// with forward differences: `deg a · p` additions and no
    /// multiplication beyond the table's `deg²`.
    pub fn roots_by_table(&self, a: &[u64]) -> Vec<u64> {
        let p = self.p;
        let d = a.len().max(1) - 1;
        if d == 0 {
            return Vec::new();
        }
        let mut v: Vec<u64> = (0..=d as u64).map(|x| self.eval(a, x)).collect();
        let mut t = Vec::with_capacity(d + 1);
        for k in 0..=d {
            t.push(v[0]);
            for i in 0..(d - k) {
                v[i] = sm(v[i + 1], v[i], p);
            }
        }
        self.adds.set(self.adds.get() + (d * (d + 1) / 2) as u64);
        let mut out = Vec::new();
        for x in 0..p {
            if t[0] == 0 {
                out.push(x);
            }
            for k in 0..d {
                t[k] = am(t[k], t[k + 1], p);
            }
        }
        self.adds.set(self.adds.get() + d as u64 * p);
        out
    }
    /// Monic irreducible factors with multiplicities of a non-zero `a`:
    /// squarefree decomposition, the linear factors by evaluation, then
    /// distinct-degree and Cantor–Zassenhaus on the root-free cofactor.
    pub fn factor(&self, a: &[u64], rng: &mut StdRng) -> Vec<(FpPoly, usize)> {
        let mut out = Vec::new();
        for (sq, mult) in self.squarefree(a) {
            let mut rest = sq;
            for r in self.roots_by_table(&rest) {
                let lin = vec![(self.p - r) % self.p, 1];
                rest = self.exact_div(&rest, &lin);
                out.push((lin, mult));
            }
            let dr = fdeg(&rest);
            if dr <= 0 {
                continue;
            }
            if dr <= 3 {
                // no root and degree 2 or 3: irreducible
                out.push((rest, mult));
                continue;
            }
            for (gd, d) in self.ddf_from(&rest, 2) {
                for f in self.edf(&gd, d, rng) {
                    out.push((f, mult));
                }
            }
        }
        out
    }
    /// The distinct roots in `F_p` of a non-zero `a`.
    pub fn roots(&self, a: &[u64], rng: &mut StdRng) -> Vec<u64> {
        let a = self.monic(a);
        if fdeg(&a) <= 0 {
            return Vec::new();
        }
        let x = vec![0u64, 1];
        let xp = self.powmod(&x, self.p as u128, &a);
        let g = self.gcd(&a, &self.sub(&xp, &x));
        let mut out: Vec<u64> = self
            .edf(&g, 1, rng)
            .into_iter()
            .map(|lin| (self.p - lin[0]) % self.p)
            .collect();
        out.sort_unstable();
        out
    }
    /// Every monic divisor of degree `k` of `Π f_i^{e_i}` (duplicates when
    /// the factorisation has them, none otherwise).
    pub fn divisors_of_degree(&self, factors: &[(FpPoly, usize)], k: usize) -> Vec<FpPoly> {
        let mut out = Vec::new();
        fn rec(
            ring: &FpRing,
            factors: &[(FpPoly, usize)],
            i: usize,
            k: usize,
            acc: FpPoly,
            out: &mut Vec<FpPoly>,
        ) {
            if k == 0 {
                out.push(acc);
                return;
            }
            if i == factors.len() {
                return;
            }
            let (f, e) = &factors[i];
            let d = fdeg(f) as usize;
            let mut cur = acc.clone();
            for j in 0..=*e {
                if j * d > k {
                    break;
                }
                rec(ring, factors, i + 1, k - j * d, cur.clone(), out);
                if j < *e {
                    cur = ring.mul(&cur, f);
                }
            }
        }
        rec(self, factors, 0, k, vec![1], &mut out);
        out
    }
}

// ── The sieve's tables and lines ─────────────────────────────────────────

/// Per-run tables: `h = h₀ + t·h₁`, the factor base's abscissae, the square
/// roots of `F_p`, `ω⁻¹`.
pub struct SieveSetup {
    pub p: u64,
    pub w: u64,
    pub winv: u64,
    pub m: usize,
    /// `⌊m/2⌋`, the degree of `A`.
    pub m1: usize,
    /// `⌊(m − 7)/2⌋`, the degree of `B`.
    pub m2: usize,
    pub h: Vec<E2>,
    pub h0: FpPoly,
    pub h1: FpPoly,
    /// `in_base[x]`: `x` is the abscissa of a factor-base class.
    pub in_base: Vec<bool>,
    /// `sqrt[a]`: a square root of `a` in `F_p`, or `u64::MAX`.
    pub sqrt: Vec<u64>,
    /// `rho[x] = √N(h(x))` for a base abscissa `x` (`N(h(x))` is a square in
    /// `F_p` exactly when `h(x)` is a square in `F_q`), `0` elsewhere.
    pub rho: Vec<u64>,
    /// `F_p` multiplications spent on the tables.
    pub table_muls: u64,
}

pub fn sieve_setup(fq: &Fq, h: &[E2], base: &[BaseEl], m: usize) -> SieveSetup {
    let p = fq.p;
    assert!(m >= 8, "m ≥ 8");
    let mut in_base = vec![false; p as usize];
    for b in base {
        in_base[b.x as usize] = true;
    }
    let mut sqrt = vec![u64::MAX; p as usize];
    for x in 0..p {
        let sq = mm(x, x, p) as usize;
        if sqrt[sq] == u64::MAX || x < sqrt[sq] {
            sqrt[sq] = x;
        }
    }
    let h0: FpPoly = trim(h.iter().map(|c| c.0[0]).collect());
    let h1: FpPoly = trim(h.iter().map(|c| c.0[1]).collect());
    let ring = FpRing::new(p);
    let mut rho = vec![0u64; p as usize];
    for b in base {
        let x = b.x;
        let (a, c) = (ring.eval(&h0, x), ring.eval(&h1, x));
        let nrm = sm(mm(a, a, p), mm(fq.w % p, mm(c, c, p), p), p);
        ring.count(4);
        let r = sqrt[nrm as usize];
        assert!(r != u64::MAX, "N(h(x)) is a square on a base abscissa");
        rho[x as usize] = r;
    }
    SieveSetup {
        p,
        w: fq.w,
        winv: inv_mod(fq.w, p),
        m,
        m1: m / 2,
        m2: (m - 7) / 2,
        h: h.to_vec(),
        h0,
        h1,
        in_base,
        sqrt,
        rho,
        table_muls: p + INV_COST + ring.muls(),
    }
}

/// The smallest `m ≥ 8` with `p^{m−6}/m! ≥ margin·|F|` (`|F|` the factor
/// base's size): the rule of the note's §11.1, fixed from `p` alone.
pub fn choose_m(p: u64, base_len: usize, margin: f64) -> usize {
    let mut m = 8usize;
    loop {
        let avail = (p as f64).powi(m as i32 - 6) / (1..=m).map(|k| k as f64).product::<f64>();
        if avail >= margin * base_len as f64 || m >= 16 {
            return m;
        }
        m += 1;
    }
}

/// One line: `A₀` (monic), `A₁′`, `B`, and the three polynomials of
/// `F(x, s) = P₁(x) + s² P₂(x) − s P₃(x)`.
#[derive(Clone, Debug)]
pub struct Line {
    pub a0: FpPoly,
    pub a1: FpPoly,
    pub b: Vec<E2>,
    /// `A₀²`.
    pub p1: FpPoly,
    /// `ω·A₁′²`.
    pub p2: FpPoly,
    /// `Re(B²h)`.
    pub p3: FpPoly,
    /// `N(B) = B₀² − ωB₁²`, whose value at a base abscissa `x` times
    /// `√N(h(x))` is the square root of the discriminant in `s`.
    pub p4: FpPoly,
}

/// The number of `B`'s the run can enumerate for this `m`.  For `m` odd, `B`
/// is monic of degree `m₂` with `m₂` free coefficients in `F_q`.  For `m`
/// even, `B` has degree `m₂` exactly and its leading coefficient runs over
/// representatives of `F_q^× / (F_p^× ∪ t·F_p^×)` (`(p + 1)/2` of them):
/// `B` and `c·B` with `c² ∈ F_p` give the same line, so this is the one `B`
/// per line.
pub fn b_space(setup: &SieveSetup) -> u128 {
    let p = setup.p as u128;
    if !setup.m.is_multiple_of(2) {
        p.pow(2 * setup.m2 as u32)
    } else {
        p.div_ceil(2) * p.pow(2 * setup.m2 as u32)
    }
}

/// The representatives of `F_q^× / (F_p^× ∪ t·F_p^×)`: `1`, and `u + t` for
/// the `u ∈ F_p^×` that are the smaller of the pair `{u, ω/u}`.
pub fn lead_reps(p: u64, w: u64) -> Vec<E2> {
    let mut out = vec![E2::ONE];
    for u in 1..p {
        let v = mm(w % p, inv_mod(u, p), p);
        if u < v {
            out.push(E2([u, 1]));
        }
    }
    out
}

/// A multiplier coprime to `n` for the enumeration's affine scramble (so that
/// `k ↦ (k·c + c/7) mod n` is a bijection).
pub fn coprime_scramble(n: u128, rng: &mut StdRng) -> u128 {
    fn gcd(mut a: u128, mut b: u128) -> u128 {
        while b != 0 {
            let t = a % b;
            a = b;
            b = t;
        }
        a
    }
    if n <= 2 {
        return 1;
    }
    loop {
        let c = rng.gen_range(2..n.min(1 << 40));
        if gcd(c, n) == 1 {
            return c;
        }
    }
}

/// The `k`-th `B` in the run's order: a bijection of `[0, b_space)` onto the
/// coefficient vectors, scrambled by an affine map coprime to `p`.
pub fn b_from_index(setup: &SieveSetup, k: u128, scramble: u128) -> Vec<E2> {
    let p = setup.p as u128;
    let n = b_space(setup);
    let mut idx = (k.wrapping_mul(scramble) + scramble / 7) % n;
    let c = setup.m2;
    let mut b: Vec<E2> = Vec::with_capacity(c + 1);
    for _ in 0..c {
        let lo = (idx % p) as u64;
        idx /= p;
        let hi = (idx % p) as u64;
        idx /= p;
        b.push(E2([lo, hi]));
    }
    if !setup.m.is_multiple_of(2) {
        b.push(E2::ONE);
    } else {
        let reps = lead_reps(setup.p, setup.w);
        b.push(reps[(idx % reps.len() as u128) as usize]);
    }
    b
}

/// The lines of one `B`: `A₀` over the monic divisors of `Im(B²h)` of the
/// admissible degrees, `A₁′ = Im(B²h)/(2A₀)`.  Counts its `F_q` work on
/// `fq` and its `F_p` work on `ring`.
pub fn lines_for_b(
    ring: &FpRing,
    fq: &Fq,
    setup: &SieveSetup,
    b: &[E2],
    rng: &mut StdRng,
) -> Vec<Line> {
    let p = setup.p;
    if b.is_empty() {
        return Vec::new();
    }
    let b2h = pmul(fq, &pmul(fq, b, b), &setup.h);
    let g: FpPoly = trim(b2h.iter().map(|c| c.0[1]).collect());
    let re: FpPoly = trim(b2h.iter().map(|c| c.0[0]).collect());
    if g.is_empty() {
        return Vec::new();
    }
    let dg = fdeg(&g) as usize;
    let (lo, hi) = if setup.m.is_multiple_of(2) {
        (setup.m1, setup.m1)
    } else {
        (dg.saturating_sub(setup.m1), setup.m1)
    };
    if lo > hi || dg < lo {
        return Vec::new();
    }
    let factors = ring.factor(&g, rng);
    let half = ring.inv(2);
    let b0: FpPoly = trim(b.iter().map(|c| c.0[0]).collect());
    let b1: FpPoly = trim(b.iter().map(|c| c.0[1]).collect());
    let nb = ring.sub(
        &ring.mul(&b0, &b0),
        &ring.scale(&ring.mul(&b1, &b1), setup.w % p),
    );
    let mut out = Vec::new();
    for k in lo..=hi.min(dg) {
        for a0 in ring.divisors_of_degree(&factors, k) {
            let a1 = ring.scale(&ring.exact_div(&g, &a0), half);
            if setup.m.is_multiple_of(2) && fdeg(&a1) >= setup.m1 as isize {
                continue;
            }
            if !setup.m.is_multiple_of(2) {
                // the identity 2A₀A₁′ = G is symmetric: (A₁′/lc, lc·A₀) is the
                // same line up to the scalar t/(ωsc), so keep one of the pair
                let a1m = ring.monic(&a1);
                if (a0.len(), &a0) > (a1m.len(), &a1m) {
                    continue;
                }
            }
            let p1 = ring.mul(&a0, &a0);
            let p2 = ring.scale(&ring.mul(&a1, &a1), setup.w % p);
            out.push(Line {
                a0,
                a1,
                b: b.to_vec(),
                p1,
                p2,
                p3: re.clone(),
                p4: nb.clone(),
            });
        }
    }
    out
}

/// Counters indexed by `s`, tagged by line so that no reset is needed.
pub struct Counters {
    tag: Vec<u32>,
    cnt: Vec<u8>,
    id: u32,
}

impl Counters {
    pub fn new(p: u64) -> Counters {
        Counters {
            tag: vec![u32::MAX; p as usize],
            cnt: vec![0; p as usize],
            id: 0,
        }
    }
}

/// What one line's sieve cost.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct SieveCost {
    /// `x` visited (the whole of `F_p`, the differences must advance).
    pub steps: u64,
    /// `x` in the factor base (the two roots in `s` computed).
    pub base_steps: u64,
    pub muls: u64,
    pub adds: u64,
    /// Counter updates (one per root).
    pub lookups: u64,
    pub roots: u64,
}

/// Forward-difference table of `P` at `x = 0`: `Δ^k P(0)`, `k = 0..deg`.
fn difference_table(ring: &FpRing, poly: &[u64], adds: &mut u64) -> Vec<u64> {
    let p = ring.p;
    let d = poly.len().max(1) - 1;
    let mut v: Vec<u64> = (0..=d as u64).map(|x| ring.eval(poly, x)).collect();
    // v[k] = Δ^k P(0): k-th forward difference at 0
    let mut table = Vec::with_capacity(d + 1);
    for k in 0..=d {
        table.push(v[0]);
        for i in 0..(d - k) {
            v[i] = sm(v[i + 1], v[i], p);
            *adds += 1;
        }
    }
    table
}

/// Sieve one line: the `s` whose counter reached `m`, and the cost.
/// `all_counts`, when given, receives the final counter of every `s` (for
/// the tests' brute-force comparison).
pub fn sieve_line(
    ring: &FpRing,
    setup: &SieveSetup,
    line: &Line,
    ctr: &mut Counters,
    all_counts: Option<&mut Vec<u8>>,
) -> (Vec<u64>, SieveCost) {
    let p = setup.p;
    let m = setup.m as u8;
    let mut cost = SieveCost::default();
    ctr.id = ctr.id.wrapping_add(1);
    if ctr.id == u32::MAX {
        ctr.tag.iter_mut().for_each(|t| *t = u32::MAX);
        ctr.id = 0;
    }
    let id = ctr.id;
    let m0 = ring.muls();
    // F(x, s) = a s² − Re s + c with a = ω A₁′(x)², c = A₀(x)²: on a base
    // abscissa the discriminant is Re² − 4ac = N(B²h)(x) = (N(B)(x)·√N(h(x)))²,
    // a square always, so every base step has its two roots
    //   s = (Re(x) ± N(B)(x)·ρ(x)) / (2a)
    // and the tables advanced are P₂ = ωA₁′², P₃ = Re(B²h), P₄ = N(B).
    let mut t2 = difference_table(ring, &line.p2, &mut cost.adds);
    let mut t3 = difference_table(ring, &line.p3, &mut cost.adds);
    let mut t4 = difference_table(ring, &line.p4, &mut cost.adds);
    let mut hits: Vec<u64> = Vec::new();
    // (Re, r, 2a) of the base steps, inverted in batches
    let mut pending: Vec<(u64, u64, u64)> = Vec::with_capacity(256);
    let bump = |s: u64, ctr: &mut Counters, hits: &mut Vec<u64>, cost: &mut SieveCost| {
        if s == 0 {
            return;
        }
        let i = s as usize;
        cost.lookups += 1;
        cost.roots += 1;
        if ctr.tag[i] != id {
            ctr.tag[i] = id;
            ctr.cnt[i] = 0;
        }
        ctr.cnt[i] += 1;
        if ctr.cnt[i] == m {
            hits.push(s);
        }
    };
    let flush = |pending: &mut Vec<(u64, u64, u64)>,
                 ctr: &mut Counters,
                 hits: &mut Vec<u64>,
                 cost: &mut SieveCost| {
        if pending.is_empty() {
            return;
        }
        // Montgomery's batch inversion of the 2a's: 3(n−1) products, one inverse
        let n = pending.len();
        let mut pre = Vec::with_capacity(n);
        let mut acc = 1u64;
        for &(_, _, a2) in pending.iter() {
            pre.push(acc);
            acc = mm(acc, a2, p);
        }
        let mut inv = inv_mod(acc, p);
        cost.muls += 3 * (n as u64 - 1) + INV_COST;
        for i in (0..n).rev() {
            let (re, r, a2) = pending[i];
            let inv_a2 = mm(inv, pre[i], p);
            inv = mm(inv, a2, p);
            let s1 = mm(am(re, r, p), inv_a2, p);
            cost.muls += 1;
            bump(s1, ctr, hits, cost);
            if r != 0 {
                let s2 = mm(sm(re, r, p), inv_a2, p);
                cost.muls += 1;
                bump(s2, ctr, hits, cost);
            }
        }
        pending.clear();
    };
    for x in 0..p {
        cost.steps += 1;
        if setup.in_base[x as usize] {
            cost.base_steps += 1;
            let a = t2[0];
            let re = t3[0];
            let r = mm(t4[0], setup.rho[x as usize], p);
            cost.muls += 1;
            if a == 0 {
                if re != 0 {
                    // A₁′(x) = 0: F is linear in s, s = A₀(x)²/Re(x)
                    let c = ring.eval(&line.a0, x);
                    let s = mm(mm(c, c, p), inv_mod(re, p), p);
                    cost.muls += INV_COST + 2;
                    bump(s, &mut *ctr, &mut hits, &mut cost);
                }
            } else {
                debug_assert_eq!(
                    mm(r, r, p),
                    sm(
                        mm(re, re, p),
                        mm(4 % p, mm(a, ring.eval(&line.p1, x), p), p),
                        p
                    ),
                    "the discriminant is N(B²h)(x)"
                );
                pending.push((re, r, am(a, a, p)));
                if pending.len() == 256 {
                    flush(&mut pending, &mut *ctr, &mut hits, &mut cost);
                }
            }
        }
        // advance the difference tables to x + 1
        for t in [&mut t2, &mut t3, &mut t4] {
            let d = t.len() - 1;
            for k in 0..d {
                t[k] = am(t[k], t[k + 1], p);
            }
            cost.adds += d as u64;
        }
    }
    flush(&mut pending, ctr, &mut hits, &mut cost);
    cost.muls += ring.muls() - m0;
    if let Some(out) = all_counts {
        out.clear();
        out.resize(p as usize, 0);
        for s in 0..p as usize {
            if ctr.tag[s] == id {
                out[s] = ctr.cnt[s];
            }
        }
    }
    (hits, cost)
}

/// A relation `Σ ε_i B_i = 0` over `m` distinct factor-base classes.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct SieveRel {
    pub terms: Vec<(usize, i8)>,
}

/// The polynomial `F(·, s)` of a line, in `F_p[x]`.
pub fn line_poly(ring: &FpRing, line: &Line, s: u64) -> FpPoly {
    let s2 = mm(s, s, ring.p);
    ring.count(1);
    let a = ring.add(&line.p1, &ring.scale(&line.p2, s2));
    ring.sub(&a, &ring.scale(&line.p3, s))
}

/// Turn a hit `s` on a line into the relation: the `m` roots of `F(·, s)`,
/// the ordinate `y_i = −A(x_i)/(√s·B(x_i))` of each zero of `f`, the base
/// index and sign.  `None` if the hit does not describe `m` distinct base
/// points (counted by the caller).
pub fn relation_from_hit(
    ring: &FpRing,
    fq: &Fq,
    setup: &SieveSetup,
    base: &[BaseEl],
    by_x: &HashMap<u64, usize>,
    line: &Line,
    s: u64,
    rng: &mut StdRng,
) -> Option<SieveRel> {
    let p = setup.p;
    let f = line_poly(ring, line, s);
    if fdeg(&f) != setup.m as isize {
        return None;
    }
    let roots = ring.roots(&f, rng);
    if roots.len() != setup.m {
        return None;
    }
    // √s ∈ F_q: in F_p if s is a square, else t·√(s/ω)
    let sq = if setup.sqrt[s as usize] != u64::MAX {
        E2([setup.sqrt[s as usize], 0])
    } else {
        let r = setup.sqrt[mm(s, setup.winv, p) as usize];
        ring.count(1);
        if r == u64::MAX {
            return None;
        }
        E2([0, r])
    };
    let mut terms = Vec::with_capacity(setup.m);
    for &x in &roots {
        let &idx = by_x.get(&x)?;
        let a = E2([ring.eval(&line.a0, x), mm(s, ring.eval(&line.a1, x), p)]);
        ring.count(1);
        let bx = peval(fq, &line.b, &fq.from_fp(x));
        let den = fq.mul(&sq, &bx);
        if den == E2::ZERO {
            return None;
        }
        let y = fq.neg(&fq.mul(&a, &fq.inv(&den)));
        let sign = if y == base[idx].y {
            1i8
        } else if y == fq.neg(&base[idx].y) {
            -1i8
        } else {
            return None;
        };
        terms.push((idx, sign));
    }
    terms.sort_unstable();
    Some(SieveRel { terms })
}

/// `Σ ε_i B_i = 0` in the Jacobian?
pub fn verify_rel(jac: &super::jv_cover::Hyp<Fq>, base: &[BaseEl], rel: &SieveRel) -> bool {
    let mut acc = jac.identity();
    for &(i, e) in &rel.terms {
        let d = if e == 1 {
            base[i].d.clone()
        } else {
            jac.neg(&base[i].d)
        };
        acc = jac.add(&acc, &d);
    }
    jac.is_identity(&acc)
}

/// The square core of homogeneous rows over `ncols` columns: singleton
/// columns eliminated to a fixed point, then the surplus trimmed by the
/// rule of [`square_core`] (the relation touching the fewest weight-2
/// columns goes first).  `None` while the rows do not cover their columns.
pub fn sieve_core(rels: &[SparseRel], ncols: usize) -> Option<super::gaudry_cubic::SquareCore> {
    let mut alive = vec![true; rels.len()];
    let weights = |alive: &[bool]| -> Vec<usize> {
        let mut w = vec![0usize; ncols];
        for (i, r) in rels.iter().enumerate() {
            if alive[i] {
                for &(c, _) in &r.cols {
                    w[c] += 1;
                }
            }
        }
        w
    };
    loop {
        loop {
            let w = weights(&alive);
            let mut dropped = false;
            for (i, r) in rels.iter().enumerate() {
                if alive[i] && r.cols.iter().any(|&(c, _)| w[c] == 1) {
                    alive[i] = false;
                    dropped = true;
                }
            }
            if !dropped {
                break;
            }
        }
        let w = weights(&alive);
        let n_rows = alive.iter().filter(|&&a| a).count();
        let n_cols = w.iter().filter(|&&x| x > 0).count();
        if n_rows == 0 || n_rows < n_cols {
            return None;
        }
        if n_rows == n_cols {
            return Some(super::gaudry_cubic::SquareCore {
                rows: (0..rels.len()).filter(|&i| alive[i]).collect(),
                columns: (0..ncols).map(|c| w[c] > 0).collect(),
            });
        }
        let victim = (0..rels.len())
            .filter(|&i| alive[i])
            .min_by_key(|&i| rels[i].cols.iter().filter(|&&(c, _)| w[c] == 2).count())?;
        alive[victim] = false;
    }
}

// ── The method end to end ────────────────────────────────────────────────

#[derive(Clone, Debug, Default, Serialize)]
pub struct SieveDlpReport {
    pub p: u64,
    pub seed: u64,
    pub l: u64,
    pub bits: f64,
    pub base: usize,
    pub c_add_e: f64,
    pub c_add_j: f64,
    /// `m` as the rule chose it from `p`, and the `m` the run ended on.
    pub m_rule: usize,
    pub m_final: usize,
    pub margin: f64,
    pub relations_available: f64,
    /// Relation phase.
    pub bs: u64,
    pub lines: u64,
    pub sieve_steps: u64,
    pub base_steps: u64,
    pub roots: u64,
    pub hits: u64,
    pub false_hits: u64,
    pub relations: usize,
    pub rels_verified: u64,
    pub rels_failed_verify: u64,
    /// Relations found again (or as their negation), not counted in `relations`.
    pub duplicates: u64,
    /// Per `m` the run went through: `(m, B's, lines, relations, hits)`.
    pub per_m: Vec<(usize, u64, u64, u64, u64)>,
    /// Measured relations per line against `p/m!`.
    pub rels_per_line: f64,
    pub expected_rels_per_line: f64,
    pub rate_ratio: f64,
    /// Phases, in `F_p` multiplications.
    pub setup_muls: u64,
    pub base_muls: u64,
    pub table_muls: u64,
    pub enum_muls: u64,
    /// Additions of the enumeration's root search by evaluation.
    pub enum_adds: u64,
    pub sieve_muls: u64,
    pub sieve_adds: u64,
    pub sieve_lookups: u64,
    pub extract_muls: u64,
    pub verify_muls: u64,
    pub relation_muls: u64,
    pub c_rel: f64,
    /// The descent: two residuals decomposed by the six-point test.
    pub descent_residuals: u64,
    pub descent_tests_per_success: f64,
    pub descent_successes: u64,
    pub descent_weil_muls: u64,
    pub descent_f4_muls: u64,
    pub descent_lin_muls: u64,
    pub descent_post_muls: u64,
    pub descent_verify_muls: u64,
    pub descent_stream_muls: u64,
    pub descent_muls: u64,
    pub descent_c_cov: f64,
    pub descent_stopped: u64,
    pub descent_incomplete: u64,
    pub descent_timed_out: u64,
    /// Descent rows found but outside the core's columns at the solve.
    pub descent_rejected: u64,
    pub unknowns: usize,
    pub filtered_out: usize,
    pub row_weight: f64,
    pub la_attempts: u64,
    pub la_ops: u64,
    pub la_muls: u64,
    pub la_ops_last: u64,
    pub total_muls: u64,
    pub total_plus: u64,
    pub s: f64,
    pub s_plus: f64,
    pub solved: bool,
    pub correct: bool,
    pub exhausted: bool,
    pub rho: Vec<RhoRunE>,
    pub rho_s_mean: f64,
    pub rho_s_sd: f64,
    pub s_over_rho: f64,
    pub s_plus_over_rho: f64,
    pub relation_over_rho: f64,
    pub descent_over_rho: f64,
    pub la_over_rho: f64,
    pub wall_ms: f64,
}

/// One batch of `B` indices sieved on one thread.
struct BatchOut {
    rels: Vec<SieveRel>,
    /// `(B index, line index within B, s)` of each relation, in order.
    prov: Vec<(u128, usize, u64)>,
    bs: u64,
    lines: u64,
    cost: SieveCost,
    hits: u64,
    false_hits: u64,
    enum_muls: u64,
    enum_adds: u64,
    extract_muls: u64,
    verify_muls: u64,
    verified: u64,
    failed: u64,
}

#[allow(clippy::too_many_arguments)]
fn sieve_batch(
    spec: &Spec,
    base: &[BaseEl],
    by_x: &HashMap<u64, usize>,
    m: usize,
    scramble: u128,
    ks: std::ops::Range<u128>,
    seed: u64,
) -> BatchOut {
    let ctx = Ctx::new(spec);
    let fq = &ctx.f.f;
    let jac = ctx.jac();
    let ring = FpRing::new(spec.p);
    let setup = sieve_setup(fq, &ctx.cov.hx, base, m);
    let mut ctr = Counters::new(spec.p);
    let mut out = BatchOut {
        rels: Vec::new(),
        prov: Vec::new(),
        bs: 0,
        lines: 0,
        cost: SieveCost::default(),
        hits: 0,
        false_hits: 0,
        enum_muls: 0,
        enum_adds: 0,
        extract_muls: 0,
        verify_muls: 0,
        verified: 0,
        failed: 0,
    };
    let mut rng =
        StdRng::seed_from_u64(seed ^ (ks.start as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    for k in ks {
        out.bs += 1;
        let b = b_from_index(&setup, k, scramble);
        fq.reset_muls();
        ring.reset();
        let lines = lines_for_b(&ring, fq, &setup, &b, &mut rng);
        out.enum_muls += fq.muls() + ring.muls();
        out.enum_adds += ring.adds();
        for (li, line) in lines.iter().enumerate() {
            out.lines += 1;
            let (hits, cost) = sieve_line(&ring, &setup, line, &mut ctr, None);
            out.cost.steps += cost.steps;
            out.cost.base_steps += cost.base_steps;
            out.cost.muls += cost.muls;
            out.cost.adds += cost.adds;
            out.cost.lookups += cost.lookups;
            out.cost.roots += cost.roots;
            for s in hits {
                out.hits += 1;
                fq.reset_muls();
                ring.reset();
                let rel = relation_from_hit(&ring, fq, &setup, base, by_x, line, s, &mut rng);
                out.extract_muls += fq.muls() + ring.muls();
                match rel {
                    None => out.false_hits += 1,
                    Some(rel) => {
                        fq.reset_muls();
                        let ok = verify_rel(&jac, base, &rel);
                        out.verify_muls += fq.muls();
                        if ok {
                            out.verified += 1;
                            out.rels.push(rel);
                            out.prov.push((k, li, s));
                        } else {
                            out.failed += 1;
                        }
                    }
                }
            }
        }
    }
    out
}

/// The route with the sieve: relations among factor-base points by the
/// sieve of the note's §11, the two residuals' descent by the six-point
/// test (F4 stopped at the Bézout staircase when `stop` says so), Wiedemann
/// mod `ℓ`, and rho on the same group (`rho_runs = 0` takes `rho_s_ref`).
#[allow(clippy::too_many_arguments)]
pub fn run_cover_sieve_dlp(
    p: u64,
    seed: u64,
    rho_runs: usize,
    rho_s_ref: f64,
    stop: Option<usize>,
    margin: f64,
    m_override: Option<usize>,
) -> SieveDlpReport {
    let start = Instant::now();
    let spec = generate_spec(p, seed);
    let ctx = Ctx::new(&spec);
    let jac = ctx.jac();
    let l = spec.l;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x51E7E);
    let (c_e, c_j) = unit_costs(&spec, &ctx);
    ctx.f.reset_muls();
    let base = factor_base(&ctx, &mut rng);
    let base_muls = ctx.f.muls();
    let by_x: HashMap<u64, usize> = base.iter().enumerate().map(|(i, b)| (b.x, i)).collect();
    let small = base.len();
    let unknowns = small + 1;
    let m_rule = choose_m(p, small, margin);
    let mut m = m_override.unwrap_or(m_rule);
    let fact = |m: usize| (1..=m).map(|k| k as f64).product::<f64>();
    let mut rep = SieveDlpReport {
        p,
        seed,
        l,
        bits: (l as f64).log2(),
        base: small,
        c_add_e: c_e,
        c_add_j: c_j,
        m_rule,
        m_final: m,
        margin,
        relations_available: (p as f64).powi(m as i32 - 6) / fact(m),
        base_muls,
        ..Default::default()
    };
    // ── the residual stream for the descent: R = aG′ + bQ′, R₀ + i·M
    ctx.f.reset_muls();
    let (al, be) = (rng.gen_range(1..l), rng.gen_range(1..l));
    let step = jac.add(
        &jac.mul(&spec.gj, al as u128),
        &jac.mul(&spec.qj, be as u128),
    );
    let mut a = rng.gen_range(0..l);
    let mut b = rng.gen_range(1..l);
    let mut r = jac.add(&jac.mul(&spec.gj, a as u128), &jac.mul(&spec.qj, b as u128));
    rep.setup_muls = ctx.f.muls() + spec.transfer_muls;
    ctx.f.reset_muls();
    // ── the relation phase: the sieve over the lines, in the run's order of B
    let setup = sieve_setup(&ctx.f.f, &ctx.cov.hx, &base, m);
    rep.table_muls = setup.table_muls;
    let scramble_for = |m: usize| -> u128 {
        let mut s = StdRng::seed_from_u64(seed ^ ((m as u64) << 32));
        let n = b_space(&sieve_setup(&ctx.f.f, &ctx.cov.hx, &base, m));
        coprime_scramble(n, &mut s)
    };
    let mut scramble = scramble_for(m);
    let mut space = b_space(&setup);
    let mut next_k: u128 = 0;
    let mut sieve_rows: Vec<SparseRel> = Vec::new();
    let mut sieve_rels: Vec<SieveRel> = Vec::new();
    let mut seen: std::collections::HashSet<Vec<(usize, i8)>> = std::collections::HashSet::new();
    let threads = rayon::current_num_threads().max(1);
    // about 4·10⁶ sieve steps per thread and round
    let chunk: u128 = (4_000_000 / p as u128).max(64);
    // sieve one round of batches; false when every m up to 16 is exhausted
    let mut sieve_round = |rep: &mut SieveDlpReport,
                           sieve_rows: &mut Vec<SparseRel>,
                           sieve_rels: &mut Vec<SieveRel>|
     -> bool {
        if next_k >= space {
            // the lines of this m are exhausted: continue with m + 1 (the rule
            // registered in §11.3), the relations found so far kept
            if m >= 16 {
                return false;
            }
            m += 1;
            rep.m_final = m;
            scramble = scramble_for(m);
            space = b_space(&sieve_setup(&ctx.f.f, &ctx.cov.hx, &base, m));
            next_k = 0;
        }
        let ranges: Vec<std::ops::Range<u128>> = (0..threads as u128)
            .map(|t| {
                let lo = (next_k + t * chunk).min(space);
                let hi = (lo + chunk).min(space);
                lo..hi
            })
            .filter(|r| r.start < r.end)
            .collect();
        next_k = ranges.last().map(|r| r.end).unwrap_or(space);
        let outs: Vec<BatchOut> = ranges
            .into_par_iter()
            .map(|ks| sieve_batch(&spec, &base, &by_x, m, scramble, ks, seed))
            .collect();
        if rep.per_m.last().map(|e| e.0) != Some(m) {
            rep.per_m.push((m, 0, 0, 0, 0));
        }
        for o in outs {
            let e = rep.per_m.last_mut().unwrap();
            e.1 += o.bs;
            e.2 += o.lines;
            e.4 += o.hits;
            rep.bs += o.bs;
            rep.lines += o.lines;
            rep.sieve_steps += o.cost.steps;
            rep.base_steps += o.cost.base_steps;
            rep.sieve_muls += o.cost.muls;
            rep.sieve_adds += o.cost.adds;
            rep.sieve_lookups += o.cost.lookups;
            rep.roots += o.cost.roots;
            rep.hits += o.hits;
            rep.false_hits += o.false_hits;
            rep.enum_muls += o.enum_muls;
            rep.enum_adds += o.enum_adds;
            rep.extract_muls += o.extract_muls;
            rep.verify_muls += o.verify_muls;
            rep.rels_verified += o.verified;
            rep.rels_failed_verify += o.failed;
            for rel in o.rels {
                // a relation found twice (or as its negation) is one relation
                let neg: Vec<(usize, i8)> = rel.terms.iter().map(|&(i, e)| (i, -e)).collect();
                if seen.contains(&rel.terms) || seen.contains(&neg) {
                    rep.duplicates += 1;
                    continue;
                }
                seen.insert(rel.terms.clone());
                rep.per_m.last_mut().unwrap().3 += 1;
                let cols: Vec<(usize, u64)> = rel
                    .terms
                    .iter()
                    .map(|&(i, e)| (i, if e == 1 { 1 } else { l - 1 }))
                    .collect();
                sieve_rows.push(SparseRel { cols, rhs: 0 });
                sieve_rels.push(rel);
            }
        }
        rep.relations = sieve_rels.len();
        if std::env::var_os("JV_SIEVE_PROGRESS").is_some() {
            eprintln!(
                "  [sieve p={p} m={m}] B's {} lines {} hits {} relations {} (need ≈ {}) enum {:.2e} sieve {:.2e} muls",
                rep.bs, rep.lines, rep.hits, rep.relations, small, rep.enum_muls as f64, rep.sieve_muls as f64
            );
        }
        true
    };
    // ── the descent: residuals decomposed by the six-point test until two
    // rows lie inside the core's columns
    let mut descent_rows: Vec<SparseRel> = Vec::new();
    let mut descend_round = |rep: &mut SieveDlpReport, descent_rows: &mut Vec<SparseRel>| {
        let mut batch: Vec<(u64, u64, Div<E2>, u64)> = Vec::with_capacity(64);
        for _ in 0..64 {
            let idx = rep.descent_residuals + batch.len() as u64;
            batch.push((a, b, r.clone(), idx));
            r = jac.add(&r, &step);
            a = (a + al) % l;
            b = (b + be) % l;
        }
        rep.descent_stream_muls += ctx.f.muls();
        ctx.f.reset_muls();
        let f4_before = f4_fp::field_ops_total();
        let found: Vec<_> = batch
            .par_iter()
            .map_init(
                || Ctx::new(&spec),
                |c, (a, b, r, idx)| {
                    let opts = opts_for(24, 600.0, stop);
                    let mut lrng =
                        StdRng::seed_from_u64(seed ^ idx.wrapping_mul(0x9E37_79B9_7F4A_7C15));
                    let (decs, cost) = nagao_decompose(c, &base, &by_x, r, &opts, &mut lrng);
                    let jc = c.jac();
                    let v0 = c.f.muls();
                    let ok: Vec<_> = decs
                        .into_iter()
                        .filter(|d| verify_dec(&jc, &base, r, d))
                        .collect();
                    let vm = c.f.muls() - v0;
                    (*a, *b, ok, cost, vm)
                },
            )
            .collect();
        rep.descent_f4_muls += f4_fp::field_ops_total() - f4_before;
        if std::env::var_os("JV_SIEVE_PROGRESS").is_some() {
            eprintln!(
                "  [descent p={p}] residuals {} successes {}",
                rep.descent_residuals + 64,
                rep.descent_successes
            );
        }
        for (a, b, decs, cost, vm) in found {
            rep.descent_residuals += 1;
            rep.descent_weil_muls += cost.weil_muls;
            rep.descent_lin_muls += cost.lin_muls;
            rep.descent_post_muls += cost.post_muls;
            rep.descent_verify_muls += vm;
            if cost.incomplete {
                rep.descent_incomplete += 1;
            }
            if cost.timed_out {
                rep.descent_timed_out += 1;
            }
            if cost.stopped_at.is_some() {
                rep.descent_stopped += 1;
            }
            if let Some(d) = decs.first() {
                rep.descent_successes += 1;
                let mut cols: Vec<(usize, u64)> = d
                    .terms
                    .iter()
                    .map(|&(i, e)| (i, if e == 1 { 1 } else { l - 1 }))
                    .collect();
                cols.push((small, (l - b) % l));
                cols.sort_unstable();
                descent_rows.push(SparseRel { cols, rhs: a });
            }
        }
    };
    let mut la_ops = 0u64;
    let mut next_attempt = small;
    let descent_cap = 400u64 * 720 * 4;
    'collect: loop {
        if !sieve_round(&mut rep, &mut sieve_rows, &mut sieve_rels) {
            rep.exhausted = true;
            break 'collect;
        }
        if sieve_rows.len() < next_attempt {
            continue;
        }
        let Some(core) = sieve_core(&sieve_rows, unknowns) else {
            next_attempt = sieve_rows.len() + (small / 20).max(1);
            continue;
        };
        // two descent rows inside the core's columns, decomposing more
        // residuals until there are (a row outside the core is counted and
        // kept for a later, larger core)
        let fits = |row: &SparseRel| row.cols.iter().all(|&(c, _)| c == small || core.columns[c]);
        while descent_rows.iter().filter(|r| fits(r)).count() < 2
            && rep.descent_residuals < descent_cap
        {
            descend_round(&mut rep, &mut descent_rows);
        }
        let usable: Vec<&SparseRel> = descent_rows.iter().filter(|r| fits(r)).take(2).collect();
        if usable.len() < 2 {
            break 'collect;
        }
        rep.descent_rejected = descent_rows.len() as u64 - 2;
        rep.la_attempts += 1;
        let mut map = vec![usize::MAX; unknowns];
        let mut k = 0;
        for c in 0..unknowns {
            if core.columns[c] {
                map[c] = k;
                k += 1;
            }
        }
        map[small] = k;
        let mapped = |row: &SparseRel| SparseRel {
            cols: row.cols.iter().map(|&(c, v)| (map[c], v)).collect(),
            rhs: row.rhs,
        };
        // the core minus one sieve row, plus the two descent rows: square
        let mut sel: Vec<SparseRel> = core.rows[..core.rows.len() - 1]
            .iter()
            .map(|&i| mapped(&sieve_rows[i]))
            .collect();
        sel.extend(usable.iter().map(|r| mapped(r)));
        let before = la_ops;
        let x = wiedemann_u64(&sel, sel.len(), l, &mut rng, &mut la_ops);
        let dd = x.map(|x| x[k]);
        let ec = super::jv_cover::EllE::new(&ctx.f, &spec.alpha);
        if dd.is_some_and(|dd| ec.mul(&spec.g, dd) == spec.q) {
            rep.solved = true;
            rep.correct = dd == Some(spec.d);
            rep.unknowns = sel.len();
            rep.filtered_out = sieve_rows.len() + 1 - core.rows.len();
            rep.la_ops_last = la_ops - before;
            break 'collect;
        }
        next_attempt = sieve_rows.len() + (small / 20).max(1);
    }
    let full_rels = sieve_rows;
    rep.descent_muls = rep.descent_weil_muls
        + rep.descent_f4_muls
        + rep.descent_lin_muls
        + rep.descent_post_muls
        + rep.descent_verify_muls
        + rep.descent_stream_muls;
    rep.descent_c_cov = (rep.descent_weil_muls
        + rep.descent_f4_muls
        + rep.descent_lin_muls
        + rep.descent_post_muls) as f64
        / rep.descent_residuals.max(1) as f64;
    rep.descent_tests_per_success =
        rep.descent_residuals as f64 / rep.descent_successes.max(1) as f64;
    rep.relations = sieve_rels.len();
    rep.expected_rels_per_line = p as f64 / fact(rep.m_final);
    rep.rels_per_line = rep.relations as f64 / rep.lines.max(1) as f64;
    rep.rate_ratio = rep.rels_per_line / rep.expected_rels_per_line;
    rep.relation_muls =
        rep.table_muls + rep.enum_muls + rep.sieve_muls + rep.extract_muls + rep.verify_muls;
    rep.c_rel = rep.relation_muls as f64 / rep.relations.max(1) as f64;
    rep.la_ops = la_ops;
    rep.la_muls = la_ops * 16;
    rep.row_weight = full_rels.iter().map(|r| r.cols.len()).sum::<usize>() as f64
        / full_rels.len().max(1) as f64;
    rep.total_muls =
        rep.setup_muls + rep.base_muls + rep.relation_muls + rep.descent_muls + rep.la_muls;
    rep.total_plus = rep.total_muls + rep.sieve_adds + rep.sieve_lookups + rep.enum_adds;
    let sqrt_l = (l as f64).sqrt();
    rep.s = rep.total_muls as f64 / (c_e * sqrt_l);
    rep.s_plus = rep.total_plus as f64 / (c_e * sqrt_l);
    rep.rho = (0..rho_runs as u64)
        .into_par_iter()
        .map(|k| rho_e(&spec, seed * 1000 + k))
        .collect();
    let ss: Vec<f64> = rep.rho.iter().map(|r| r.s).collect();
    rep.rho_s_mean = if ss.is_empty() {
        rho_s_ref
    } else {
        ss.iter().sum::<f64>() / ss.len() as f64
    };
    rep.rho_s_sd = if ss.len() < 2 {
        0.0
    } else {
        (ss.iter().map(|s| (s - rep.rho_s_mean).powi(2)).sum::<f64>() / (ss.len() - 1) as f64)
            .sqrt()
    };
    let rho_muls = rep.rho_s_mean * sqrt_l * c_e;
    rep.s_over_rho = rep.total_muls as f64 / rho_muls;
    rep.s_plus_over_rho = rep.total_plus as f64 / rho_muls;
    rep.relation_over_rho = (rep.setup_muls + rep.base_muls + rep.relation_muls) as f64 / rho_muls;
    rep.descent_over_rho = rep.descent_muls as f64 / rho_muls;
    rep.la_over_rho = rep.la_muls as f64 / rho_muls;
    rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
    rep
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::jv_cover::Hyp;

    fn ring_and_rng(p: u64) -> (FpRing, StdRng) {
        (FpRing::new(p), StdRng::seed_from_u64(7))
    }

    #[test]
    fn factoring_recovers_a_planted_product() {
        let (ring, mut rng) = ring_and_rng(1009);
        for trial in 0..20 {
            // random monic factors of degrees 1..4 with multiplicities
            let mut prod = vec![1u64];
            let mut planted: Vec<(FpPoly, usize)> = Vec::new();
            let nf = 2 + trial % 3;
            for _ in 0..nf {
                let d = 1 + rng.gen_range(0..4usize);
                let mut f: FpPoly = (0..d).map(|_| rng.gen_range(0..1009)).collect();
                f.push(1);
                let e = 1 + rng.gen_range(0..2usize);
                for _ in 0..e {
                    prod = ring.mul(&prod, &f);
                }
                planted.push((f, e));
            }
            let got = ring.factor(&prod, &mut rng);
            let mut back = vec![1u64];
            for (f, e) in &got {
                assert!(fdeg(f) >= 1);
                for _ in 0..*e {
                    back = ring.mul(&back, f);
                }
            }
            assert_eq!(back, prod, "trial {trial}");
            // every reported factor is irreducible: no root unless linear, and
            // for degree 2..4 no factor of smaller degree
            for (f, _) in &got {
                let d = fdeg(f) as usize;
                if d > 1 {
                    assert!(ring.roots(f, &mut rng).is_empty(), "{f:?}");
                    let dd = ring.ddf(f);
                    assert_eq!(dd, vec![(f.clone(), d)], "{f:?}");
                }
            }
            // the divisor enumeration multiplies out
            let d2 = ring.divisors_of_degree(&got, 2);
            for d in &d2 {
                assert_eq!(fdeg(d), 2);
                assert!(ring.rem(&prod, d).is_empty());
            }
        }
    }

    #[test]
    fn roots_are_exact() {
        let (ring, mut rng) = ring_and_rng(503);
        for _ in 0..20 {
            let k = 1 + rng.gen_range(0..6usize);
            let mut rs: Vec<u64> = (0..k).map(|_| rng.gen_range(0..503)).collect();
            rs.sort_unstable();
            rs.dedup();
            let mut f = vec![1u64];
            for &r in &rs {
                f = ring.mul(&f, &[(503 - r) % 503, 1]);
            }
            // times an irreducible quadratic: x² − ω
            let w = Fq::new(503).w;
            f = ring.mul(&f, &[(503 - w) % 503, 0, 1]);
            assert_eq!(ring.roots(&f, &mut rng), rs);
        }
    }

    fn small_instance(p: u64, seed: u64) -> (Spec, Ctx, Vec<BaseEl>, HashMap<u64, usize>) {
        let spec = generate_spec(p, seed);
        let ctx = Ctx::new(&spec);
        let mut rng = StdRng::seed_from_u64(seed);
        let base = factor_base(&ctx, &mut rng);
        let by_x = base.iter().enumerate().map(|(i, b)| (b.x, i)).collect();
        (spec, ctx, base, by_x)
    }

    #[test]
    fn the_line_identity_holds_in_fq() {
        // F(x, s) computed in F_q from f = A + By is in F_p[x] and equals the
        // line's P₁ + s²P₂ − sP₃, for m = 9 and m = 10
        let (spec, ctx, base, _) = small_instance(101, 3);
        let fq = &ctx.f.f;
        let ring = FpRing::new(spec.p);
        let mut rng = StdRng::seed_from_u64(11);
        for m in [9usize, 10] {
            let setup = sieve_setup(fq, &ctx.cov.hx, &base, m);
            let scramble = coprime_scramble(b_space(&setup), &mut rng);
            let mut lines_seen = 0;
            for k in 0..200u128 {
                let b = b_from_index(&setup, k, scramble);
                for line in lines_for_b(&ring, fq, &setup, &b, &mut rng) {
                    lines_seen += 1;
                    for _ in 0..3 {
                        let s = rng.gen_range(1..spec.p);
                        // A = A₀ + t·s·A₁′, B′ = √s·B
                        let a: Vec<E2> = (0..line.a0.len().max(line.a1.len()))
                            .map(|i| {
                                E2([
                                    line.a0.get(i).copied().unwrap_or(0),
                                    mm(s, line.a1.get(i).copied().unwrap_or(0), spec.p),
                                ])
                            })
                            .collect();
                        let sq = fq.sqrt(&fq.from_fp(s), &mut rng).expect("√s ∈ F_q");
                        let bb: Vec<E2> = line.b.iter().map(|c| fq.mul(c, &sq)).collect();
                        let a2 = pmul(fq, &a, &a);
                        let b2h = pmul(fq, &pmul(fq, &bb, &bb), &ctx.cov.hx);
                        let f: Vec<E2> = (0..a2.len().max(b2h.len()))
                            .map(|i| {
                                fq.sub(
                                    &a2.get(i).copied().unwrap_or(E2::ZERO),
                                    &b2h.get(i).copied().unwrap_or(E2::ZERO),
                                )
                            })
                            .collect();
                        assert!(f.iter().all(|c| c.in_fp()), "F ∉ F_p[x] at m = {m}");
                        let f_p: FpPoly = trim(f.iter().map(|c| c.0[0]).collect());
                        assert_eq!(f_p, line_poly(&ring, &line, s), "m = {m}");
                        assert_eq!(fdeg(&f_p), m as isize);
                    }
                }
            }
            assert!(lines_seen >= 20, "m = {m}: {lines_seen} lines in 200 B's");
        }
    }

    #[test]
    fn the_sieve_counts_what_brute_force_counts() {
        // every counter equals the number of base abscissae that are roots of
        // F(·, s), for every s, on lines at p = 53 and p = 101
        for (p, m) in [(53u64, 12usize), (101, 11), (101, 10)] {
            let (_spec, ctx, base, _) = small_instance(p, 5);
            let fq = &ctx.f.f;
            let ring = FpRing::new(p);
            let mut rng = StdRng::seed_from_u64(9);
            let setup = sieve_setup(fq, &ctx.cov.hx, &base, m);
            let mut ctr = Counters::new(p);
            let scr = coprime_scramble(b_space(&setup), &mut rng);
            let mut lines_seen = 0;
            for k in 0..60u128 {
                let b = b_from_index(&setup, k, scr);
                for line in lines_for_b(&ring, fq, &setup, &b, &mut rng) {
                    lines_seen += 1;
                    let mut counts = Vec::new();
                    let (hits, _) = sieve_line(&ring, &setup, &line, &mut ctr, Some(&mut counts));
                    for s in 1..p {
                        let f = line_poly(&ring, &line, s);
                        // an x where A₀ and B share a root is a root of F(·, s) for
                        // every s (a double one); the sieve skips it, as the
                        // relation would contain P + ι(P)
                        let brute = (0..p)
                            .filter(|&x| {
                                setup.in_base[x as usize]
                                    && ring.eval(&f, x) == 0
                                    && !(ring.eval(&line.p1, x) == 0
                                        && ring.eval(&line.p2, x) == 0
                                        && ring.eval(&line.p3, x) == 0)
                            })
                            .count() as u8;
                        assert_eq!(counts[s as usize], brute, "p = {p}, m = {m}, s = {s}");
                        assert_eq!(hits.contains(&s), brute as usize == m);
                    }
                }
            }
            assert!(lines_seen >= 5, "p = {p}: {lines_seen} lines");
        }
    }

    #[test]
    fn relations_sum_to_zero_in_the_jacobian() {
        // p = 251, m = 9 (p/9! ≈ 7·10⁻⁴ relations per line): a few relations,
        // each verified by Cantor arithmetic
        let (spec, ctx, base, by_x) = small_instance(251, 2);
        let fq = &ctx.f.f;
        let jac: Hyp<Fq> = ctx.jac();
        let ring = FpRing::new(spec.p);
        let mut rng = StdRng::seed_from_u64(4);
        let setup = sieve_setup(fq, &ctx.cov.hx, &base, 9);
        let scr = coprime_scramble(b_space(&setup), &mut rng);
        let mut ctr = Counters::new(spec.p);
        let (mut rels, mut false_hits, mut lines) = (0, 0, 0u64);
        let mut k = 0u128;
        while rels < 4 && k < 200_000 {
            let b = b_from_index(&setup, k, scr);
            k += 1;
            for line in lines_for_b(&ring, fq, &setup, &b, &mut rng) {
                lines += 1;
                let (hits, _) = sieve_line(&ring, &setup, &line, &mut ctr, None);
                for s in hits {
                    match relation_from_hit(&ring, fq, &setup, &base, &by_x, &line, s, &mut rng) {
                        None => false_hits += 1,
                        Some(rel) => {
                            assert_eq!(rel.terms.len(), 9);
                            assert!(verify_rel(&jac, &base, &rel), "{rel:?}");
                            rels += 1;
                        }
                    }
                }
            }
        }
        assert!(rels >= 4, "{rels} relations in {lines} lines");
        assert_eq!(false_hits, 0);
    }

    #[test]
    fn m_rule_matches_the_registration() {
        // the note's §11.1: 12 at 53, 11 at 101, 10 at 251, 9 at 503, 1009, 1511
        for (p, m) in [
            (53u64, 12usize),
            (101, 11),
            (251, 10),
            (503, 9),
            (1009, 9),
            (1511, 9),
        ] {
            assert_eq!(choose_m(p, (p / 2) as usize, 1.25), m, "p = {p}");
        }
    }

    #[test]
    #[ignore]
    fn end_to_end_at_p101_solves_the_logarithm() {
        let r = run_cover_sieve_dlp(101, 1, 4, 1.3, Some(64), 1.25, None);
        assert!(r.solved && r.correct, "{r:?}");
        assert_eq!(r.m_final, 11);
        assert_eq!(r.rels_failed_verify, 0);
    }
}
