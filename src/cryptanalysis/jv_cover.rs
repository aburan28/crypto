//! **The cover-and-decomposition route on `E(F_{p⁶})`: a genus-3 cover over
//! `F_{p²}` and Nagao decompositions on its Jacobian, measured end to end.**
//!
//! Companion to `research/notes/index-calculus/RESEARCH_COVER_DECOMPOSITION_LEDGER.md`
//! (registered before this module existed) and a reproduction of Joux and
//! Vitse, *Cover and decomposition index calculus on elliptic curves made
//! practical* (Eurocrypt 2012).  The curves are the elliptic curves over
//! `F_{q³}`, `q = p²`, of the form `y² = x(x − α)(x − σ(α))`, `α ∈ F_{q³} ∖ F_q`,
//! whose GHS cover is the genus-3 hyperelliptic curve `H : y² = F(x)N(x)`
//! over `F_q` (`N` the minimal polynomial of `α`).  The conorm–norm map
//! `E(F_{q³}) → Jac_H(F_q)` carries the DLP over at unit cost; on `Jac_H(F_q)`
//! the factor base is `{(x, y) : x ∈ F_p}` (`≈ p` points, `≈ p/2` classes under
//! `y ↦ −y`) and a Nagao decomposition writes a divisor as a sum of `ng = 6`
//! of them, the test being a quadratic system of six equations in six
//! unknowns over `F_p` (the Weil restriction of the six coefficients of a
//! monic sextic over `F_q` that must lie in `F_p`).  A residual decomposes
//! with probability `1/720`.
//!
//! Everything is counted in `F_p` multiplications through the counter of the
//! base field `F_{p²}`, which every operation of the tower reaches.

use std::cell::Cell;
use std::collections::{BTreeSet, HashMap};
use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::Serialize;

use super::f4_fp::{self, F4Options, F4Trace, Ordering as F4Ordering};
use super::gaudry_cubic::{am, mm, sm, square_core, wiedemann_u64, PolyRing, SparseRel, UPoly};
use super::residual_walk::{inv_mod, is_prime_u64, mix64, pow_mod};
use std::sync::{Arc, Mutex, OnceLock};

// ── A field, abstractly ──────────────────────────────────────────────────

pub trait Fld {
    type E: Copy + Eq + std::fmt::Debug;
    fn zero(&self) -> Self::E;
    fn one(&self) -> Self::E;
    fn add(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn sub(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn neg(&self, a: &Self::E) -> Self::E;
    fn mul(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn inv(&self, a: &Self::E) -> Self::E;
    fn is_zero(&self, a: &Self::E) -> bool {
        *a == self.zero()
    }
}

// ── F_{p²} = F_p[t] / (t² − ω) ───────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default, PartialOrd, Ord)]
pub struct E2(pub [u64; 2]);

impl E2 {
    pub const ZERO: E2 = E2([0, 0]);
    pub const ONE: E2 = E2([1, 0]);
    pub fn in_fp(&self) -> bool {
        self.0[1] == 0
    }
}

/// `F_{p²}` with a counter of `F_p` multiplications: a product is `5` (four
/// products and the scaling by `ω`), an inverse `21` (the norm, one `F_p`
/// inversion at the harness's `16`, two scalings).
#[derive(Clone)]
pub struct Fq {
    pub p: u64,
    pub w: u64,
    muls: Cell<u64>,
}

impl Fq {
    pub fn new(p: u64) -> Fq {
        assert!(p % 2 == 1 && p > 3);
        let w = (2..p)
            .find(|&w| pow_mod(w, (p - 1) / 2, p) != 1)
            .expect("a non-residue");
        Fq {
            p,
            w,
            muls: Cell::new(0),
        }
    }
    fn count(&self, k: u64) {
        self.muls.set(self.muls.get() + k);
    }
    pub fn count_public(&self, k: u64) {
        self.count(k);
    }
    pub fn muls(&self) -> u64 {
        self.muls.get()
    }
    pub fn reset_muls(&self) {
        self.muls.set(0);
    }
    pub fn from_fp(&self, x: u64) -> E2 {
        E2([x % self.p, 0])
    }
    /// Multiplication by an element of `F_p`: two products.
    pub fn scale(&self, a: &E2, k: u64) -> E2 {
        self.count(2);
        E2([mm(a.0[0], k, self.p), mm(a.0[1], k, self.p)])
    }
    pub fn sq(&self, a: &E2) -> E2 {
        self.mul(a, a)
    }
    pub fn conj(&self, a: &E2) -> E2 {
        E2([a.0[0], (self.p - a.0[1]) % self.p])
    }
    pub fn pow(&self, a: &E2, mut e: u128) -> E2 {
        let mut base = *a;
        let mut acc = E2::ONE;
        while e > 0 {
            if e & 1 == 1 {
                acc = self.mul(&acc, &base);
            }
            base = self.sq(&base);
            e >>= 1;
        }
        acc
    }
    fn order(&self) -> u128 {
        (self.p as u128) * (self.p as u128) - 1
    }
    pub fn is_square(&self, a: &E2) -> bool {
        *a == E2::ZERO || self.pow(a, self.order() / 2) == E2::ONE
    }
    pub fn random(&self, rng: &mut StdRng) -> E2 {
        E2([rng.gen_range(0..self.p), rng.gen_range(0..self.p)])
    }
    /// Tonelli–Shanks in `F_{p²}^×`.
    pub fn sqrt(&self, a: &E2, rng: &mut StdRng) -> Option<E2> {
        if *a == E2::ZERO {
            return Some(E2::ZERO);
        }
        if !self.is_square(a) {
            return None;
        }
        let q1 = self.order();
        let s = q1.trailing_zeros();
        let odd = q1 >> s;
        let z = loop {
            let z = self.random(rng);
            if z != E2::ZERO && !self.is_square(&z) {
                break z;
            }
        };
        let mut m = s;
        let mut cc = self.pow(&z, odd);
        let mut t = self.pow(a, odd);
        let mut r = self.pow(a, odd.div_ceil(2));
        while t != E2::ONE {
            let mut i = 0;
            let mut tt = t;
            while tt != E2::ONE {
                tt = self.sq(&tt);
                i += 1;
            }
            let mut b = cc;
            for _ in 0..(m - i - 1) {
                b = self.sq(&b);
            }
            m = i;
            cc = self.sq(&b);
            t = self.mul(&t, &cc);
            r = self.mul(&r, &b);
        }
        Some(r)
    }
}

impl Fld for Fq {
    type E = E2;
    fn zero(&self) -> E2 {
        E2::ZERO
    }
    fn one(&self) -> E2 {
        E2::ONE
    }
    fn add(&self, a: &E2, b: &E2) -> E2 {
        E2([am(a.0[0], b.0[0], self.p), am(a.0[1], b.0[1], self.p)])
    }
    fn sub(&self, a: &E2, b: &E2) -> E2 {
        E2([sm(a.0[0], b.0[0], self.p), sm(a.0[1], b.0[1], self.p)])
    }
    fn neg(&self, a: &E2) -> E2 {
        E2([(self.p - a.0[0]) % self.p, (self.p - a.0[1]) % self.p])
    }
    fn mul(&self, a: &E2, b: &E2) -> E2 {
        let p = self.p;
        self.count(5);
        let m00 = mm(a.0[0], b.0[0], p);
        let m11 = mm(a.0[1], b.0[1], p);
        E2([
            am(m00, mm(self.w, m11, p), p),
            am(mm(a.0[0], b.0[1], p), mm(a.0[1], b.0[0], p), p),
        ])
    }
    fn inv(&self, a: &E2) -> E2 {
        assert!(*a != E2::ZERO, "inverse of zero in F_p²");
        let p = self.p;
        self.count(3 + 16 + 2);
        let norm = sm(
            mm(a.0[0], a.0[0], p),
            mm(self.w, mm(a.0[1], a.0[1], p), p),
            p,
        );
        let ni = inv_mod(norm, p);
        E2([mm(a.0[0], ni, p), mm((p - a.0[1]) % p, ni, p)])
    }
}

// ── F_{p⁶} = F_{p²}[θ] / (θ³ − s) ────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default, PartialOrd, Ord)]
pub struct E6(pub [E2; 3]);

impl E6 {
    pub const ZERO: E6 = E6([E2::ZERO; 3]);
    pub const ONE: E6 = E6([E2::ONE, E2::ZERO, E2::ZERO]);
    pub fn in_fq(&self) -> bool {
        self.0[1] == E2::ZERO && self.0[2] == E2::ZERO
    }
}

/// `F_{q³}`, `q = p²`, over `F_q`: `s ∈ F_q` a non-cube (so `θ³ − s` is
/// irreducible; it needs `3 | q − 1`, i.e. `p ≢ 0 (mod 3)`), `ζ = s^{(q−1)/3}`
/// a primitive cube root of unity with `θ^q = ζθ`, so the `q`-Frobenius is
/// `σ(a₀ + a₁θ + a₂θ²) = a₀ + ζa₁θ + ζ²a₂θ²`.
#[derive(Clone)]
pub struct Fq3 {
    pub f: Fq,
    pub s: E2,
    pub zeta: E2,
    pub zeta2: E2,
}

impl Fq3 {
    pub fn new(p: u64) -> Fq3 {
        assert!(!p.is_multiple_of(3));
        let f = Fq::new(p);
        let q1 = f.order();
        let s = (1..p * p)
            .map(|k| E2([k % p, (k / p) % p]))
            .find(|s| *s != E2::ZERO && f.pow(s, q1 / 3) != E2::ONE)
            .expect("a non-cube");
        let zeta = f.pow(&s, q1 / 3);
        let zeta2 = f.mul(&zeta, &zeta);
        f.reset_muls();
        Fq3 { f, s, zeta, zeta2 }
    }
    pub fn muls(&self) -> u64 {
        self.f.muls()
    }
    pub fn reset_muls(&self) {
        self.f.reset_muls()
    }
    pub fn from_fq(&self, x: E2) -> E6 {
        E6([x, E2::ZERO, E2::ZERO])
    }
    /// Scale by an element of `F_q`: three products.
    pub fn scale(&self, a: &E6, k: &E2) -> E6 {
        E6(core::array::from_fn(|i| self.f.mul(&a.0[i], k)))
    }
    pub fn sq(&self, a: &E6) -> E6 {
        self.mul(a, a)
    }
    /// The `q`-Frobenius, two products.
    pub fn sigma(&self, a: &E6) -> E6 {
        E6([
            a.0[0],
            self.f.mul(&a.0[1], &self.zeta),
            self.f.mul(&a.0[2], &self.zeta2),
        ])
    }
    pub fn pow(&self, a: &E6, mut e: u128) -> E6 {
        let mut base = *a;
        let mut acc = E6::ONE;
        while e > 0 {
            if e & 1 == 1 {
                acc = self.mul(&acc, &base);
            }
            base = self.sq(&base);
            e >>= 1;
        }
        acc
    }
    fn order(&self) -> u128 {
        (self.f.p as u128).pow(6) - 1
    }
    pub fn is_square(&self, a: &E6) -> bool {
        *a == E6::ZERO || self.pow(a, self.order() / 2) == E6::ONE
    }
    pub fn random(&self, rng: &mut StdRng) -> E6 {
        E6(core::array::from_fn(|_| self.f.random(rng)))
    }
    pub fn sqrt(&self, a: &E6, rng: &mut StdRng) -> Option<E6> {
        if *a == E6::ZERO {
            return Some(E6::ZERO);
        }
        if !self.is_square(a) {
            return None;
        }
        let q1 = self.order();
        let s = q1.trailing_zeros();
        let odd = q1 >> s;
        let z = loop {
            let z = self.random(rng);
            if z != E6::ZERO && !self.is_square(&z) {
                break z;
            }
        };
        let mut m = s;
        let mut cc = self.pow(&z, odd);
        let mut t = self.pow(a, odd);
        let mut r = self.pow(a, odd.div_ceil(2));
        while t != E6::ONE {
            let mut i = 0;
            let mut tt = t;
            while tt != E6::ONE {
                tt = self.sq(&tt);
                i += 1;
            }
            let mut b = cc;
            for _ in 0..(m - i - 1) {
                b = self.sq(&b);
            }
            m = i;
            cc = self.sq(&b);
            t = self.mul(&t, &cc);
            r = self.mul(&r, &b);
        }
        Some(r)
    }
}

impl Fld for Fq3 {
    type E = E6;
    fn zero(&self) -> E6 {
        E6::ZERO
    }
    fn one(&self) -> E6 {
        E6::ONE
    }
    fn add(&self, a: &E6, b: &E6) -> E6 {
        E6(core::array::from_fn(|i| self.f.add(&a.0[i], &b.0[i])))
    }
    fn sub(&self, a: &E6, b: &E6) -> E6 {
        E6(core::array::from_fn(|i| self.f.sub(&a.0[i], &b.0[i])))
    }
    fn neg(&self, a: &E6) -> E6 {
        E6(core::array::from_fn(|i| self.f.neg(&a.0[i])))
    }
    /// Schoolbook with `θ³ = s`: nine products and two reductions.
    fn mul(&self, a: &E6, b: &E6) -> E6 {
        let f = &self.f;
        let mut d = [E2::ZERO; 5];
        for i in 0..3 {
            for j in 0..3 {
                d[i + j] = f.add(&d[i + j], &f.mul(&a.0[i], &b.0[j]));
            }
        }
        E6([
            f.add(&d[0], &f.mul(&self.s, &d[3])),
            f.add(&d[1], &f.mul(&self.s, &d[4])),
            d[2],
        ])
    }
    /// Through the norm to `F_q`: `a⁻¹ = σ(a)σ²(a) / N(a)`.
    fn inv(&self, a: &E6) -> E6 {
        assert!(*a != E6::ZERO, "inverse of zero in F_p⁶");
        let s1 = self.sigma(a);
        let s2 = self.sigma(&s1);
        let b = self.mul(&s1, &s2);
        let n = self.mul(a, &b);
        debug_assert!(n.in_fq());
        let ni = self.f.inv(&n.0[0]);
        self.scale(&b, &ni)
    }
}

// ── Polynomials over a field, low degree first, no trailing zeros ─────────

pub type Poly<E> = Vec<E>;

pub(crate) fn ptrim<F: Fld>(f: &F, p: &mut Poly<F::E>) {
    while p.last().is_some_and(|c| f.is_zero(c)) {
        p.pop();
    }
}
pub(crate) fn padd<F: Fld>(f: &F, a: &[F::E], b: &[F::E]) -> Poly<F::E> {
    let n = a.len().max(b.len());
    let mut out: Poly<F::E> = (0..n)
        .map(|i| match (a.get(i), b.get(i)) {
            (Some(x), Some(y)) => f.add(x, y),
            (Some(x), None) => *x,
            (None, Some(y)) => *y,
            _ => f.zero(),
        })
        .collect();
    ptrim(f, &mut out);
    out
}
pub(crate) fn psub<F: Fld>(f: &F, a: &[F::E], b: &[F::E]) -> Poly<F::E> {
    let n = a.len().max(b.len());
    let mut out: Poly<F::E> = (0..n)
        .map(|i| match (a.get(i), b.get(i)) {
            (Some(x), Some(y)) => f.sub(x, y),
            (Some(x), None) => *x,
            (None, Some(y)) => f.neg(y),
            _ => f.zero(),
        })
        .collect();
    ptrim(f, &mut out);
    out
}
pub(crate) fn pmul<F: Fld>(f: &F, a: &[F::E], b: &[F::E]) -> Poly<F::E> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![f.zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if f.is_zero(x) {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            out[i + j] = f.add(&out[i + j], &f.mul(x, y));
        }
    }
    ptrim(f, &mut out);
    out
}
pub(crate) fn pscale<F: Fld>(f: &F, a: &[F::E], k: &F::E) -> Poly<F::E> {
    let mut out: Poly<F::E> = a.iter().map(|x| f.mul(x, k)).collect();
    ptrim(f, &mut out);
    out
}
/// `a = q b + r`, `deg r < deg b`.
pub(crate) fn pdivrem<F: Fld>(f: &F, a: &[F::E], b: &[F::E]) -> (Poly<F::E>, Poly<F::E>) {
    assert!(!b.is_empty(), "division by the zero polynomial");
    let mut r: Poly<F::E> = a.to_vec();
    ptrim(f, &mut r);
    if r.len() < b.len() {
        return (Vec::new(), r);
    }
    let lci = f.inv(b.last().unwrap());
    let mut q = vec![f.zero(); r.len() - b.len() + 1];
    while r.len() >= b.len() {
        let k = r.len() - b.len();
        let c = f.mul(r.last().unwrap(), &lci);
        q[k] = c;
        for (i, x) in b.iter().enumerate() {
            let t = f.mul(&c, x);
            r[k + i] = f.sub(&r[k + i], &t);
        }
        // the leading term is now exactly zero
        r.pop();
        ptrim(f, &mut r);
    }
    ptrim(f, &mut q);
    (q, r)
}
pub(crate) fn pmonic<F: Fld>(f: &F, a: &[F::E]) -> Poly<F::E> {
    match a.last() {
        None => Vec::new(),
        Some(lc) => {
            let i = f.inv(lc);
            let mut out: Poly<F::E> = a.iter().map(|x| f.mul(x, &i)).collect();
            *out.last_mut().unwrap() = f.one();
            out
        }
    }
}
/// `(g, s, t)` with `g = s a + t b` monic (or zero).
pub(crate) fn pxgcd<F: Fld>(f: &F, a: &[F::E], b: &[F::E]) -> (Poly<F::E>, Poly<F::E>, Poly<F::E>) {
    let (mut r0, mut r1) = (a.to_vec(), b.to_vec());
    let (mut s0, mut s1): (Poly<F::E>, Poly<F::E>) = (vec![f.one()], Vec::new());
    let (mut t0, mut t1): (Poly<F::E>, Poly<F::E>) = (Vec::new(), vec![f.one()]);
    ptrim(f, &mut r0);
    ptrim(f, &mut r1);
    while !r1.is_empty() {
        let (q, r) = pdivrem(f, &r0, &r1);
        let s = psub(f, &s0, &pmul(f, &q, &s1));
        let t = psub(f, &t0, &pmul(f, &q, &t1));
        r0 = std::mem::replace(&mut r1, r);
        s0 = std::mem::replace(&mut s1, s);
        t0 = std::mem::replace(&mut t1, t);
    }
    match r0.last() {
        None => (r0, s0, t0),
        Some(lc) => {
            let i = f.inv(lc);
            (pscale(f, &r0, &i), pscale(f, &s0, &i), pscale(f, &t0, &i))
        }
    }
}
pub(crate) fn peval<F: Fld>(f: &F, a: &[F::E], x: &F::E) -> F::E {
    a.iter()
        .rev()
        .fold(f.zero(), |acc, c| f.add(&f.mul(&acc, x), c))
}
/// `a mod m` (`m ≠ 0`).
pub(crate) fn pmod<F: Fld>(f: &F, a: &[F::E], m: &[F::E]) -> Poly<F::E> {
    pdivrem(f, a, m).1
}

// ── Jacobians of y² = f(x), f monic of odd degree, by Cantor ──────────────

/// A divisor class in Mumford form `(u, v)`, `u` monic, `v² ≡ f (mod u)`;
/// reduced when `deg u ≤ g`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Div<E> {
    pub u: Poly<E>,
    pub v: Poly<E>,
}

pub struct Hyp<'a, F: Fld> {
    pub f: &'a F,
    /// `f(x)`, monic of degree `2g + 1`.
    pub h: Poly<F::E>,
    pub g: usize,
}

impl<'a, F: Fld> Hyp<'a, F> {
    pub fn identity(&self) -> Div<F::E> {
        Div {
            u: vec![self.f.one()],
            v: Vec::new(),
        }
    }
    pub fn is_identity(&self, d: &Div<F::E>) -> bool {
        d.u.len() == 1
    }
    pub fn neg(&self, d: &Div<F::E>) -> Div<F::E> {
        Div {
            u: d.u.clone(),
            v: d.v.iter().map(|c| self.f.neg(c)).collect(),
        }
    }
    /// Is `(u, v)` a valid Mumford pair for this curve (`u` monic, `u | f − v²`)?
    pub fn is_valid(&self, d: &Div<F::E>) -> bool {
        let f = self.f;
        if d.u.last() != Some(&f.one()) {
            return false;
        }
        let r = pmod(f, &psub(f, &self.h, &pmul(f, &d.v, &d.v)), &d.u);
        r.is_empty() && d.v.len() < d.u.len().max(1)
    }
    /// Cantor composition, not reduced.
    pub fn compose(&self, a: &Div<F::E>, b: &Div<F::E>) -> Div<F::E> {
        let f = self.f;
        let (d1, e1, e2) = pxgcd(f, &a.u, &b.u);
        let vs = padd(f, &a.v, &b.v);
        let (d, c1, c2) = pxgcd(f, &d1, &vs);
        let s1 = pmul(f, &c1, &e1);
        let s2 = pmul(f, &c1, &e2);
        let s3 = c2;
        let d2 = pmul(f, &d, &d);
        let (u, r) = pdivrem(f, &pmul(f, &a.u, &b.u), &d2);
        debug_assert!(r.is_empty());
        let t1 = pmul(f, &pmul(f, &s1, &a.u), &b.v);
        let t2 = pmul(f, &pmul(f, &s2, &b.u), &a.v);
        let t3 = pmul(f, &s3, &padd(f, &pmul(f, &a.v, &b.v), &self.h));
        let num = padd(f, &padd(f, &t1, &t2), &t3);
        let (vq, vr) = pdivrem(f, &num, &d);
        debug_assert!(vr.is_empty());
        let u = pmonic(f, &u);
        let v = pmod(f, &vq, &u);
        Div { u, v }
    }
    pub fn reduce(&self, mut d: Div<F::E>) -> Div<F::E> {
        let f = self.f;
        while d.u.len() as isize - 1 > self.g as isize {
            let num = psub(f, &self.h, &pmul(f, &d.v, &d.v));
            let (u2, r) = pdivrem(f, &num, &d.u);
            debug_assert!(r.is_empty());
            let u2 = pmonic(f, &u2);
            let nv: Poly<F::E> = d.v.iter().map(|c| f.neg(c)).collect();
            let v2 = pmod(f, &nv, &u2);
            d = Div { u: u2, v: v2 };
        }
        d.v = pmod(f, &d.v, &d.u);
        d
    }
    pub fn add(&self, a: &Div<F::E>, b: &Div<F::E>) -> Div<F::E> {
        if self.is_identity(a) {
            return b.clone();
        }
        if self.is_identity(b) {
            return a.clone();
        }
        self.reduce(self.compose(a, b))
    }
    pub fn sub(&self, a: &Div<F::E>, b: &Div<F::E>) -> Div<F::E> {
        self.add(a, &self.neg(b))
    }
    pub fn mul(&self, a: &Div<F::E>, mut k: u128) -> Div<F::E> {
        let mut acc = self.identity();
        let mut base = a.clone();
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(&acc, &base);
            }
            k >>= 1;
            if k > 0 {
                base = self.add(&base, &base);
            }
        }
        acc
    }
}

// ── The weak curve y² = x³ + a₂x² + a₄x over F_{q³}, a₂ = −(α+σα), a₄ = α·σα

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PtE6 {
    pub x: E6,
    pub y: E6,
    pub inf: bool,
}

impl PtE6 {
    pub const INF: PtE6 = PtE6 {
        x: E6::ZERO,
        y: E6::ZERO,
        inf: true,
    };
}

pub struct EllE<'a> {
    pub f: &'a Fq3,
    pub a2: E6,
    pub a4: E6,
    ops: Cell<u64>,
}

impl<'a> EllE<'a> {
    pub fn new(f: &'a Fq3, alpha: &E6) -> EllE<'a> {
        let sa = f.sigma(alpha);
        EllE {
            f,
            a2: f.neg(&f.add(alpha, &sa)),
            a4: f.mul(alpha, &sa),
            ops: Cell::new(0),
        }
    }
    /// `y² = x³ + a₂x² + a₄x` from its coefficients (the model of a curve
    /// with full 2-torsion after one root is moved to `0`; §17's walk).
    pub fn from_a2_a4(f: &'a Fq3, a2: E6, a4: E6) -> EllE<'a> {
        EllE {
            f,
            a2,
            a4,
            ops: Cell::new(0),
        }
    }
    pub fn ops(&self) -> u64 {
        self.ops.get()
    }
    pub fn reset_ops(&self) {
        self.ops.set(0);
    }
    pub fn rhs(&self, x: &E6) -> E6 {
        let f = self.f;
        let x2 = f.sq(x);
        f.add(
            &f.add(&f.mul(&x2, x), &f.mul(&self.a2, &x2)),
            &f.mul(&self.a4, x),
        )
    }
    pub fn on_curve(&self, p: &PtE6) -> bool {
        p.inf || self.f.sq(&p.y) == self.rhs(&p.x)
    }
    pub fn neg(&self, p: &PtE6) -> PtE6 {
        if p.inf {
            *p
        } else {
            PtE6 {
                x: p.x,
                y: self.f.neg(&p.y),
                inf: false,
            }
        }
    }
    /// Affine addition: one inversion, three products.
    pub fn add(&self, p: &PtE6, q: &PtE6) -> PtE6 {
        let f = self.f;
        if p.inf {
            return *q;
        }
        if q.inf {
            return *p;
        }
        self.ops.set(self.ops.get() + 1);
        let lam = if p.x == q.x {
            if p.y == f.neg(&q.y) || p.y == E6::ZERO {
                return PtE6::INF;
            }
            // (3x² + 2a₂x + a₄) / 2y
            let x2 = f.sq(&p.x);
            let num = f.add(
                &f.add(
                    &f.add(&x2, &f.add(&x2, &x2)),
                    &f.mul(&f.add(&self.a2, &self.a2), &p.x),
                ),
                &self.a4,
            );
            f.mul(&num, &f.inv(&f.add(&p.y, &p.y)))
        } else {
            f.mul(&f.sub(&q.y, &p.y), &f.inv(&f.sub(&q.x, &p.x)))
        };
        let x3 = f.sub(&f.sub(&f.sub(&f.sq(&lam), &self.a2), &p.x), &q.x);
        let y3 = f.sub(&f.mul(&lam, &f.sub(&p.x, &x3)), &p.y);
        PtE6 {
            x: x3,
            y: y3,
            inf: false,
        }
    }
    pub fn sub(&self, p: &PtE6, q: &PtE6) -> PtE6 {
        self.add(p, &self.neg(q))
    }
    /// [`EllE::mul`] for a scalar above 2⁶⁴ (the group order itself).
    pub fn mul_u128(&self, pt: &PtE6, mut k: u128) -> PtE6 {
        let mut acc = PtE6::INF;
        let mut base = *pt;
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(&acc, &base);
            }
            k >>= 1;
            if k > 0 {
                base = self.add(&base, &base);
            }
        }
        acc
    }
    pub fn mul(&self, pt: &PtE6, mut k: u64) -> PtE6 {
        let mut acc = PtE6::INF;
        let mut base = *pt;
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(&acc, &base);
            }
            k >>= 1;
            if k > 0 {
                base = self.add(&base, &base);
            }
        }
        acc
    }
    pub fn random_point(&self, rng: &mut StdRng) -> PtE6 {
        loop {
            let x = self.f.random(rng);
            if let Some(y) = self.f.sqrt(&self.rhs(&x), rng) {
                if y != E6::ZERO {
                    return PtE6 { x, y, inf: false };
                }
            }
        }
    }
}

/// `#E(F_{q³})` by baby-step giant-step on a random point, accepted when it
/// is four times a prime (the group is then `Z/2 × Z/2ℓ` and the DLP lives
/// in the subgroup of order `ℓ`).
fn group_order_4_prime(ec: &EllE, rng: &mut StdRng) -> Option<u64> {
    // returns ℓ = #E/4 when #E = 4ℓ with ℓ prime;
    // the order is near p⁶, which passes 2⁶⁴ at p = 1,622 while ℓ = #E/4 stays
    // below it up to p = 2,039: the search runs in u128 (note §16)
    let p = ec.f.f.p as u128;
    let q3 = p.pow(6);
    let two_sqrt = 2 * isqrt128(q3) + 2;
    let lo = q3 + 1 - two_sqrt;
    let width = 2 * two_sqrt;
    let pt = ec.random_point(rng);
    let steps = isqrt128(width) + 1;
    let mut table: HashMap<PtE6, u128> = HashMap::new();
    let mut jp = PtE6::INF;
    for j in 0..steps {
        if !jp.inf {
            table.entry(jp).or_insert(j);
        }
        jp = ec.add(&jp, &pt);
    }
    let giant = ec.mul_u128(&pt, steps);
    let mut t = ec.mul_u128(&pt, lo);
    let mut i = 0u128;
    let mut found = None;
    while i * steps <= width + steps {
        let m = if t.inf {
            Some(lo + i * steps)
        } else {
            table.get(&ec.neg(&t)).map(|&j| lo + i * steps + j)
        };
        if let Some(m) = m {
            found = Some(m);
            break;
        }
        t = ec.add(&t, &giant);
        i += 1;
    }
    let m = found?;
    if m % 4 != 0 || m / 4 > u64::MAX as u128 || !is_prime_u64((m / 4) as u64) {
        return None;
    }
    let other = ec.random_point(rng);
    ec.mul_u128(&other, m).inf.then_some((m / 4) as u64)
}

fn isqrt128(n: u128) -> u128 {
    let mut r = (n as f64).sqrt() as u128;
    while r * r > n {
        r -= 1;
    }
    while (r + 1) * (r + 1) <= n {
        r += 1;
    }
    r
}

// ── The genus-3 cover H : y² = F(x)·N(x) over F_q ─────────────────────────

/// Everything about `α` the cover needs, over `F_q` and `F_{q³}`.
pub struct Cover {
    pub alpha: E6,
    pub sal: E6,
    pub s2al: E6,
    /// `N(x) = (x−α)(x−σα)(x−σ²α) ∈ F_q[x]`.
    pub nx: Poly<E2>,
    /// `F(x) = N(x)(x + φ + φ^σ + φ^{σ²})` (quartic, monic) and `f_H = F·N`.
    pub fx: Poly<E2>,
    pub hx: Poly<E2>,
    /// `D₀ = (α − σ²α)(σα − σ²α)` and its conjugates `D₁ = σD₀`, `D₂ = σ²D₀`.
    pub d: [E6; 3],
}

impl Cover {
    /// `None` when `α ∈ F_q` (no cover) or the polynomials do not land in `F_q`.
    pub fn new(f: &Fq3, alpha: &E6) -> Option<Cover> {
        if alpha.in_fq() {
            return None;
        }
        let sal = f.sigma(alpha);
        let s2al = f.sigma(&sal);
        let lin = |r: &E6| -> Poly<E6> { vec![f.neg(r), E6::ONE] };
        let n6 = pmul(f, &pmul(f, &lin(alpha), &lin(&sal)), &lin(&s2al));
        let d0 = f.mul(&f.sub(alpha, &s2al), &f.sub(&sal, &s2al));
        let d1 = f.sigma(&d0);
        let d2 = f.sigma(&d1);
        let trace = f.add(&f.add(alpha, &sal), &s2al);
        let xt: Poly<E6> = vec![trace, E6::ONE];
        let mut fx6 = pmul(f, &n6, &xt);
        let terms = [
            (&d0, pmul(f, &lin(alpha), &lin(&sal))),
            (&d1, pmul(f, &lin(&sal), &lin(&s2al))),
            (&d2, pmul(f, &lin(alpha), &lin(&s2al))),
        ];
        for (di, poly) in terms {
            fx6 = padd(f, &fx6, &pscale(f, &poly, di));
        }
        let to_fq = |p: &Poly<E6>| -> Option<Poly<E2>> {
            p.iter().map(|c| c.in_fq().then_some(c.0[0])).collect()
        };
        let nx = to_fq(&n6)?;
        let fx = to_fq(&fx6)?;
        let hx = pmul(&f.f, &fx, &nx);
        Some(Cover {
            alpha: *alpha,
            sal,
            s2al,
            nx,
            fx,
            hx,
            d: [d0, d1, d2],
        })
    }
    pub fn embed(&self, f: &Fq3, p: &[E2]) -> Poly<E6> {
        p.iter().map(|c| f.from_fq(*c)).collect()
    }
    /// The cover map `π(x, y) = (F(x)/(4N(x)), y·A₁A₂/(8N(x)²))` with
    /// `A₁ = (x−α)² − D₁` and `A₂ = (x−σα)² − D₂` (the paper's
    /// `y(x−φ^σ(x))(x−φ^{σ²}(x))/(8N(x)(x−σ²α))` with the factors of `N`
    /// cancelled); `None` at a pole.
    pub fn map_point(&self, f: &Fq3, x: &E6, y: &E6) -> Option<(E6, E6)> {
        let n6 = self.embed(f, &self.nx);
        let f6 = self.embed(f, &self.fx);
        let nx = peval(f, &n6, x);
        if nx == E6::ZERO {
            return None;
        }
        let fx = peval(f, &f6, x);
        let four = f.from_fq(f.f.from_fp(4));
        let eight = f.from_fq(f.f.from_fp(8));
        let xx = f.mul(&fx, &f.inv(&f.mul(&four, &nx)));
        let a1 = f.sub(&f.sq(&f.sub(x, &self.alpha)), &self.d[1]);
        let a2 = f.sub(&f.sq(&f.sub(x, &self.sal)), &self.d[2]);
        let yy = f.mul(
            &f.mul(y, &f.mul(&a1, &a2)),
            &f.inv(&f.mul(&eight, &f.sq(&nx))),
        );
        Some((xx, yy))
    }
    /// `π*(P)` for an affine `P = (X₀, Y₀)` of `E`: the four points of `H`
    /// over it as the Mumford pair `(G, v)` over `F_{q³}`, `G = F − 4X₀N`
    /// (monic quartic) and `v ≡ 8Y₀N²/(A₁A₂) (mod G)`.  `None` when
    /// `A₁A₂` is not invertible mod `G`.
    pub fn pullback(&self, f: &Fq3, pt: &PtE6) -> Option<Div<E6>> {
        let n6 = self.embed(f, &self.nx);
        let f6 = self.embed(f, &self.fx);
        let four = f.from_fq(f.f.from_fp(4));
        let eight = f.from_fq(f.f.from_fp(8));
        let g = psub(f, &f6, &pscale(f, &n6, &f.mul(&four, &pt.x)));
        let lin = |r: &E6| -> Poly<E6> { vec![f.neg(r), E6::ONE] };
        let a1 = psub(
            f,
            &pmul(f, &lin(&self.alpha), &lin(&self.alpha)),
            &[self.d[1]],
        );
        let a2 = psub(f, &pmul(f, &lin(&self.sal), &lin(&self.sal)), &[self.d[2]]);
        let a12 = pmod(f, &pmul(f, &a1, &a2), &g);
        let (gc, s, _) = pxgcd(f, &a12, &g);
        if gc.len() != 1 {
            return None;
        }
        let num = pscale(f, &pmul(f, &n6, &n6), &f.mul(&eight, &pt.y));
        let v = pmod(f, &pmul(f, &num, &s), &g);
        Some(Div { u: g, v })
    }
}

impl Cover {
    /// Frobenius on a Mumford pair over `F_{q³}`.
    pub fn sigma_div(&self, f: &Fq3, d: &Div<E6>) -> Div<E6> {
        Div {
            u: d.u.iter().map(|c| f.sigma(c)).collect(),
            v: d.v.iter().map(|c| f.sigma(c)).collect(),
        }
    }
    /// The conorm–norm transfer `Φ : E(F_{q³}) → Jac_H(F_q)`:
    /// `Φ(P) = Σ_{i<3} σ^i(π*P − 4∞) + W`, `W = (N, 0)` the class of
    /// `Σ(β_i, 0) − 3∞ = π*(O) − 4∞`, of order two (`−3W = W`).
    /// `None` on the point at infinity's degenerate fibre or a non-invertible
    /// denominator.
    pub fn transfer(&self, f: &Fq3, pt: &PtE6) -> Option<Div<E2>> {
        let h6 = self.embed(f, &self.hx);
        let jac3 = Hyp { f, h: h6, g: 3 };
        let jacq = Hyp {
            f: &f.f,
            h: self.hx.clone(),
            g: 3,
        };
        if pt.inf {
            return Some(jacq.identity());
        }
        let d0 = self.pullback(f, pt)?;
        let d1 = self.sigma_div(f, &d0);
        let d2 = self.sigma_div(f, &d1);
        debug_assert!(jac3.is_valid(&d0) && jac3.is_valid(&d1) && jac3.is_valid(&d2));
        let s = jac3.reduce(jac3.compose(&jac3.compose(&d0, &d1), &d2));
        let to_fq = |p: &Poly<E6>| -> Option<Poly<E2>> {
            p.iter().map(|c| c.in_fq().then_some(c.0[0])).collect()
        };
        let sq = Div {
            u: to_fq(&s.u)?,
            v: to_fq(&s.v)?,
        };
        let w = Div {
            u: self.nx.clone(),
            v: Vec::new(),
        };
        Some(jacq.add(&sq, &w))
    }
}

// ── The Weil restriction of the Nagao system ──────────────────────────────
//
// For a reduced divisor `D = (u, v)` of degree 3 the functions
// `f = (λ₀ + λ₁x)·u + (µ₀ + x)(y − v) = A(x) + µ(x)·y`, `A = (λ₀+λ₁x)u − µ v`,
// vanish on `D` and on six further points (pole order 9 at ∞); their
// abscissae are the roots of `F̃ = (µ² f_H − A²) / u`, monic of degree 6, whose
// coefficients are quadratic forms in `(λ₀, λ₁, µ₀) ∈ F_q³`.  Writing
// `λ₀ = a₀ + t b₀`, `λ₁ = a₁ + t b₁`, `µ₀ = c₀ + t d₀` over `F_p` and asking the
// six lower coefficients to lie in `F_p` gives six quadrics in six unknowns.

/// Monomials of a quadratic form in `(λ₀, λ₁, µ₀)`: `1, λ₀, λ₁, µ₀, λ₀², λ₀λ₁,
/// λ₀µ₀, λ₁², λ₁µ₀, µ₀²`.
type Q10 = [E2; 10];
const MONO: [&[usize]; 10] = [
    &[],
    &[0],
    &[1],
    &[2],
    &[0, 0],
    &[0, 1],
    &[0, 2],
    &[1, 1],
    &[1, 2],
    &[2, 2],
];

/// An affine form `c + a·λ₀ + b·λ₁ + d·µ₀`.
type Aff = [E2; 4];

fn aff_mul(f: &Fq, a: &Aff, b: &Aff) -> Q10 {
    let m = |x: usize, y: usize| f.mul(&a[x], &b[y]);
    let sum = |x: E2, y: E2| f.add(&x, &y);
    [
        m(0, 0),
        sum(m(0, 1), m(1, 0)),
        sum(m(0, 2), m(2, 0)),
        sum(m(0, 3), m(3, 0)),
        m(1, 1),
        sum(m(1, 2), m(2, 1)),
        sum(m(1, 3), m(3, 1)),
        m(2, 2),
        sum(m(2, 3), m(3, 2)),
        m(3, 3),
    ]
}
fn q_add(f: &Fq, a: &Q10, b: &Q10) -> Q10 {
    core::array::from_fn(|i| f.add(&a[i], &b[i]))
}
fn q_sub(f: &Fq, a: &Q10, b: &Q10) -> Q10 {
    core::array::from_fn(|i| f.sub(&a[i], &b[i]))
}
fn q_scale(f: &Fq, a: &Q10, k: &E2) -> Q10 {
    core::array::from_fn(|i| f.mul(&a[i], k))
}

/// The `F_p`-polynomials `Re` and `Im` of each of the ten monomials in the six
/// unknowns `(a₀, b₀, a₁, b₁, c₀, d₀)`: `(X + tY)(X′ + tY′) = XX′ + ωYY′ + t(XY′ + YX′)`.
pub struct ReIm {
    re: Vec<Vec<([u32; 6], u64)>>,
    im: Vec<Vec<([u32; 6], u64)>>,
}

impl ReIm {
    pub fn new(p: u64, w: u64) -> ReIm {
        let var = |l: usize, y: bool| -> [u32; 6] {
            let mut e = [0u32; 6];
            e[2 * l + y as usize] = 1;
            e
        };
        let add = |a: [u32; 6], b: [u32; 6]| -> [u32; 6] { core::array::from_fn(|i| a[i] + b[i]) };
        let mut re = Vec::new();
        let mut im = Vec::new();
        for m in MONO {
            match m.len() {
                0 => {
                    re.push(vec![([0u32; 6], 1u64)]);
                    im.push(vec![]);
                }
                1 => {
                    re.push(vec![(var(m[0], false), 1)]);
                    im.push(vec![(var(m[0], true), 1)]);
                }
                _ => {
                    let (i, j) = (m[0], m[1]);
                    re.push(vec![
                        (add(var(i, false), var(j, false)), 1),
                        (add(var(i, true), var(j, true)), w % p),
                    ]);
                    im.push(vec![
                        (add(var(i, false), var(j, true)), 1),
                        (add(var(i, true), var(j, false)), 1),
                    ]);
                }
            }
        }
        ReIm { re, im }
    }
}

/// `q₀, …, q₆`: the quotient `(µ² f_H − A²)/u` as quadratic forms, `q₆ = 1`.
pub struct WeilForms {
    pub q: [Q10; 7],
}

/// Build the quotient forms for the divisor `(u, v)` (`u` monic of degree 3,
/// `v` of degree ≤ 2) on `y² = f_H`.
pub fn weil_forms(f: &Fq, hx: &[E2], u: &[E2], v: &[E2]) -> WeilForms {
    let z = E2::ZERO;
    let uc = |k: isize| -> E2 {
        if k >= 0 && (k as usize) < u.len() {
            u[k as usize]
        } else {
            z
        }
    };
    let vc = |k: isize| -> E2 {
        if k >= 0 && (k as usize) < v.len() {
            v[k as usize]
        } else {
            z
        }
    };
    let fc = |k: isize| -> E2 {
        if k >= 0 && (k as usize) < hx.len() {
            hx[k as usize]
        } else {
            z
        }
    };
    // A_k = −v_{k−1} + λ₀ u_k + λ₁ u_{k−1} − µ₀ v_k, k = 0..4
    let a: Vec<Aff> = (0..5isize)
        .map(|k| [f.neg(&vc(k - 1)), uc(k), uc(k - 1), f.neg(&vc(k))])
        .collect();
    let mut n: Vec<Q10> = vec![[z; 10]; 10];
    // µ² f_H: µ₀² f_m + 2µ₀ f_{m−1} + f_{m−2}
    for (m, nm) in n.iter_mut().enumerate() {
        let m = m as isize;
        nm[0] = fc(m - 2);
        let f1 = fc(m - 1);
        nm[3] = f.add(&f1, &f1);
        nm[9] = fc(m);
    }
    // − A²
    for i in 0..5 {
        for j in i..5 {
            let mut prod = aff_mul(f, &a[i], &a[j]);
            if i != j {
                prod = q_add(f, &prod, &prod);
            }
            n[i + j] = q_sub(f, &n[i + j], &prod);
        }
    }
    // divide by u (monic, degree 3): quotient of degree 6
    let mut q: [Q10; 7] = [[z; 10]; 7];
    for m in (3..=9).rev() {
        let c = n[m];
        q[m - 3] = c;
        for j in 0..3 {
            let t = q_scale(f, &c, &uc(j as isize));
            n[m - 3 + j] = q_sub(f, &n[m - 3 + j], &t);
        }
    }
    WeilForms { q }
}

impl WeilForms {
    /// The six quadrics `Im q_k = 0`, `k = 0..5`, in `(a₀, b₀, a₁, b₁, c₀, d₀)`.
    pub fn polys(&self, ri: &ReIm, p: u64) -> Vec<f4_fp::Poly> {
        (0..6)
            .map(|k| {
                let mut terms: Vec<(Vec<u32>, u64)> = Vec::new();
                for (mi, c) in self.q[k].iter().enumerate() {
                    // c·m: Im = c_r·Im(m) + c_i·Re(m)
                    for (e, k0) in &ri.im[mi] {
                        terms.push((e.to_vec(), mm(*k0, c.0[0], p)));
                    }
                    for (e, k0) in &ri.re[mi] {
                        terms.push((e.to_vec(), mm(*k0, c.0[1], p)));
                    }
                }
                f4_fp::normalise(&terms, p, F4Ordering::Grevlex)
            })
            .filter(|f| !f.is_empty())
            .collect()
    }
    /// `q_k` at `(λ₀, λ₁, µ₀)`.
    pub fn eval(&self, f: &Fq, k: usize, lam: &[E2; 3]) -> E2 {
        let mut acc = E2::ZERO;
        for (mi, c) in self.q[k].iter().enumerate() {
            let mut t = *c;
            for &i in MONO[mi] {
                t = f.mul(&t, &lam[i]);
            }
            acc = f.add(&acc, &t);
        }
        acc
    }
}

// ── A zero-dimensional solver on top of F4: multiplication matrices ───────

#[derive(Clone, Debug, Default, Serialize)]
pub struct ZeroDimStats {
    pub degree_reached: u32,
    pub solving_degree: u32,
    pub max_rows: usize,
    pub max_cols: usize,
    pub f4_muls: u64,
    /// `F_p` multiplications of the normal forms, the Krylov minimal
    /// polynomial, the roots and the eigenvectors.
    pub lin_muls: u64,
    /// Standard monomials: the degree of the ideal.
    pub delta: usize,
    pub inconsistent: bool,
    /// Not zero-dimensional, truncated by the degree bound, or too large.
    pub incomplete: bool,
    pub timed_out: bool,
    /// A root of the minimal polynomial whose eigenspace was not a line.
    pub nonseparating: usize,
    pub f4_ms: f64,
    /// Wall time of everything after F4 in the solver.
    pub lin_ms: f64,
    /// F4 stopped early at this staircase (`F4Options::stop_staircase`).
    pub stopped_at: Option<usize>,
    /// F4 replayed a recorded trace (the useful rows of an earlier system
    /// of the same shape) instead of selecting pairs.
    pub replayed: bool,
    /// The replay found other leading monomials and the full run was done.
    pub trace_mismatch: bool,
}

/// A recorded trace with its record: replays that held, replays that
/// diverged.  A trace recorded on a system of an uncommon shape diverges on
/// most of the rest; once it has diverged three times and more often than
/// it held, it is dropped and the next full run records a new one.
pub struct TraceEntry {
    pub trace: Arc<F4Trace>,
    pub held: u64,
    pub diverged: u64,
}

/// Traces of stopped F4 runs on the six-quadric systems, one per `p`
/// (every residual's system has the same shape up to rare variants),
/// recorded by a full run and replayed on the rest when a context asks.
fn trace_cache() -> &'static Mutex<HashMap<u64, TraceEntry>> {
    static CACHE: OnceLock<Mutex<HashMap<u64, TraceEntry>>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

/// Forget the recorded traces (tests).
pub fn clear_traces() {
    trace_cache().lock().unwrap().clear();
}

/// All solutions over `F_p` of a system whose ideal is zero-dimensional: F4 to
/// a Gröbner basis, the multiplication matrices on the standard monomials, the
/// minimal polynomial of a random linear form by Krylov, its roots in `F_p`,
/// and the eigenvectors of the transposed matrices at those roots.
pub fn solve_zero_dim(
    input: &[f4_fp::Poly],
    n: usize,
    p: u64,
    opts: &F4Options,
    rng: &mut StdRng,
) -> (Vec<Vec<u64>>, ZeroDimStats) {
    solve_zero_dim_traced(input, n, p, opts, rng, false)
}

/// [`solve_zero_dim`] with F4's trace recorded on the first system of a
/// `p` and replayed on the later ones when `trace` is set; a replay whose
/// shape diverges falls back to the full run (both counted).
pub fn solve_zero_dim_traced(
    input: &[f4_fp::Poly],
    n: usize,
    p: u64,
    opts: &F4Options,
    rng: &mut StdRng,
    trace: bool,
) -> (Vec<Vec<u64>>, ZeroDimStats) {
    let mut st = ZeroDimStats::default();
    let r = if !trace {
        f4_fp::f4(input, n, p, opts)
    } else {
        let have = trace_cache()
            .lock()
            .unwrap()
            .get(&p)
            .map(|e| e.trace.clone());
        match have {
            Some(t) => {
                // a diverging replay finishes as a full run inside the engine
                let r = f4_fp::f4_replay(input, n, p, opts, &t);
                let mut cache = trace_cache().lock().unwrap();
                if let Some(e) = cache.get_mut(&p) {
                    if Arc::ptr_eq(&e.trace, &t) {
                        if r.trace_mismatch {
                            e.diverged += 1;
                            if e.diverged >= 3 && e.diverged > e.held {
                                cache.remove(&p);
                            }
                        } else {
                            e.held += 1;
                        }
                    }
                }
                if r.trace_mismatch {
                    st.trace_mismatch = true;
                } else {
                    st.replayed = true;
                }
                r
            }
            None => {
                let (r, t) = f4_fp::f4_record(input, n, p, opts);
                if !r.timed_out && !r.inconsistent && r.pairs_above_bound == 0 {
                    trace_cache()
                        .lock()
                        .unwrap()
                        .entry(p)
                        .or_insert_with(|| TraceEntry {
                            trace: Arc::new(t),
                            held: 0,
                            diverged: 0,
                        });
                }
                r
            }
        }
    };
    let t_lin = Instant::now();
    st.f4_muls = r.field_ops;
    st.degree_reached = r.degree_reached;
    st.solving_degree = r.solving_degree;
    st.max_rows = r.max_rows;
    st.max_cols = r.max_cols;
    st.f4_ms = r.ms;
    st.stopped_at = r.staircase_at_stop;
    if r.timed_out {
        st.timed_out = true;
        return (Vec::new(), st);
    }
    if r.inconsistent {
        st.inconsistent = true;
        return (Vec::new(), st);
    }
    if r.pairs_above_bound > 0 {
        st.incomplete = true;
        return (Vec::new(), st);
    }
    let lms: Vec<&Vec<u32>> = r.basis.iter().map(|g| &g[0].0).collect();
    // the smallest pure power of each variable
    let mut cap = vec![u32::MAX; n];
    for lm in &lms {
        if let Some(k) = (0..n).find(|&k| lm[k] > 0 && (0..n).all(|i| i == k || lm[i] == 0)) {
            cap[k] = cap[k].min(lm[k]);
        }
    }
    if cap.contains(&u32::MAX) {
        st.incomplete = true;
        return (Vec::new(), st);
    }
    // standard monomials by pruned depth-first search
    let mut std_mons: Vec<Vec<u32>> = Vec::new();
    fn walk(
        k: usize,
        mono: &mut Vec<u32>,
        cap: &[u32],
        lms: &[&Vec<u32>],
        out: &mut Vec<Vec<u32>>,
    ) -> bool {
        if lms
            .iter()
            .any(|l| l.iter().zip(mono.iter()).all(|(a, b)| a <= b))
        {
            return true;
        }
        if k == cap.len() {
            out.push(mono.clone());
            return out.len() <= 4096;
        }
        for e in 0..cap[k] {
            mono[k] = e;
            if !walk(k + 1, mono, cap, lms, out) {
                mono[k] = 0;
                return false;
            }
        }
        mono[k] = 0;
        true
    }
    let mut mono = vec![0u32; n];
    if !walk(0, &mut mono, &cap, &lms, &mut std_mons) {
        st.incomplete = true;
        return (Vec::new(), st);
    }
    let delta = std_mons.len();
    st.delta = delta;
    if delta == 0 {
        st.inconsistent = true;
        return (Vec::new(), st);
    }
    // Monomials packed eight bits per variable; divisibility by the guard-bit
    // trick (exponents stay below 128).
    let pack = |e: &[u32]| -> u64 {
        e.iter()
            .enumerate()
            .fold(0u64, |a, (k, &x)| a | ((x as u64) << (8 * k)))
    };
    const GUARD: u64 = 0x8080_8080_8080_8080;
    let divides_packed = |a: u64, b: u64| -> bool { ((b | GUARD) - a) & GUARD == GUARD };
    let mask: u64 = if n >= 8 {
        u64::MAX
    } else {
        (1u64 << (8 * n)) - 1
    };
    let guard = GUARD & mask;
    let divides_p =
        |a: u64, b: u64| -> bool { (((b & mask) | guard) - (a & mask)) & guard == guard };
    let _ = divides_packed;
    let index: crate::cryptanalysis::fx_hash::FxMap<u64, usize> = std_mons
        .iter()
        .enumerate()
        .map(|(i, m)| (pack(m), i))
        .collect();
    let bas: Vec<(u64, Vec<(u64, u64)>)> = r
        .basis
        .iter()
        .map(|g| {
            (
                pack(&g[0].0),
                g[1..].iter().map(|(e, c)| (pack(e), *c)).collect(),
            )
        })
        .collect();
    // NF(m) as a dense vector over the standard monomials, memoised: a standard
    // monomial is itself; otherwise `m = c·lm(g)` and `NF(m) = −Σ g_t NF(c·t)`.
    let mut memo: crate::cryptanalysis::fx_hash::FxMap<u64, Vec<u64>> = Default::default();
    let mut cnt = 0u64;
    #[allow(clippy::too_many_arguments)]
    fn nf_mono(
        m: u64,
        delta: usize,
        p: u64,
        index: &crate::cryptanalysis::fx_hash::FxMap<u64, usize>,
        bas: &[(u64, Vec<(u64, u64)>)],
        divides_p: &dyn Fn(u64, u64) -> bool,
        memo: &mut crate::cryptanalysis::fx_hash::FxMap<u64, Vec<u64>>,
        cnt: &mut u64,
    ) -> Option<Vec<u64>> {
        if let Some(&k) = index.get(&m) {
            let mut v = vec![0u64; delta];
            v[k] = 1;
            return Some(v);
        }
        if let Some(v) = memo.get(&m) {
            return Some(v.clone());
        }
        let (lm, tail) = bas.iter().find(|(lm, _)| divides_p(*lm, m))?;
        let cof = m - lm;
        let mut acc = vec![0u64; delta];
        for (t, c) in tail {
            let sub = nf_mono(t + cof, delta, p, index, bas, divides_p, memo, cnt)?;
            let neg = p - c;
            for (a, x) in acc.iter_mut().zip(&sub) {
                if *x != 0 {
                    *a = am(*a, mm(neg, *x, p), p);
                    *cnt += 1;
                }
            }
        }
        memo.insert(m, acc.clone());
        Some(acc)
    }
    // multiplication matrices: column j of M_i is NF(x_i · b_j)
    let mut mats: Vec<Vec<Vec<u64>>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut mi = vec![vec![0u64; delta]; delta];
        for (j, b) in std_mons.iter().enumerate() {
            let xb = pack(b) + (1u64 << (8 * i));
            let Some(col) = nf_mono(xb, delta, p, &index, &bas, &divides_p, &mut memo, &mut cnt)
            else {
                st.incomplete = true;
                return (Vec::new(), st);
            };
            for (k, &c) in col.iter().enumerate() {
                mi[k][j] = c;
            }
        }
        mats.push(mi);
    }
    // A random linear form `ℓ = Σ c_i x_i` and the characteristic polynomial of
    // its matrix (exact, by Hessenberg reduction: a Krylov minimal polynomial
    // of one random vector misses a root with probability `1/p` per root).  If
    // some root's eigenspace is not a line, `ℓ` does not separate the points
    // and a new form is drawn.
    let ring = PolyRing::new(p);
    let mut sols: Vec<Vec<u64>> = Vec::new();
    for _attempt in 0..16 {
        let coef: Vec<u64> = (0..n).map(|_| rng.gen_range(1..p)).collect();
        let mut m = vec![vec![0u64; delta]; delta];
        for i in 0..n {
            for a in 0..delta {
                for b in 0..delta {
                    m[a][b] = am(m[a][b], mm(coef[i], mats[i][a][b], p), p);
                }
            }
        }
        cnt += (n * delta * delta) as u64;
        let cp = charpoly_mod_p(&m, p, &mut cnt);
        let roots = ring.roots(&UPoly(cp), rng);
        let mut round: Vec<Vec<u64>> = Vec::new();
        let mut clean = true;
        for rt in roots {
            // null space of Mᵀ − rt·I
            let mut a: Vec<Vec<u64>> = (0..delta)
                .map(|i| {
                    (0..delta)
                        .map(|j| if i == j { sm(m[j][i], rt, p) } else { m[j][i] })
                        .collect()
                })
                .collect();
            let mut mc = 0u64;
            let piv = rref_mod_p_local(&mut a, p, &mut mc);
            cnt += mc;
            if piv.len() + 1 != delta {
                clean = false;
                break;
            }
            let free = (0..delta).find(|c| !piv.contains(c)).unwrap();
            let mut w = vec![0u64; delta];
            w[free] = 1;
            for (ri, &pc) in piv.iter().enumerate() {
                w[pc] = (p - a[ri][free]) % p;
            }
            let Some(kk) = w.iter().position(|&x| x != 0) else {
                continue;
            };
            let wi = inv_mod(w[kk], p);
            let mut sol = vec![0u64; n];
            for i in 0..n {
                // (M_iᵀ w)_kk = Σ_a M_i[a][kk] w_a
                let mut acc = 0u64;
                for aa in 0..delta {
                    acc = am(acc, mm(mats[i][aa][kk], w[aa], p), p);
                }
                sol[i] = mm(acc, wi, p);
            }
            cnt += (n * (delta + 1)) as u64;
            // a genuine solution of the input
            if input.iter().all(|f| f4_fp::eval(f, &sol, p) == 0) {
                round.push(sol);
            }
        }
        if clean {
            sols = round;
            break;
        }
        st.nonseparating += 1;
        if _attempt == 15 {
            st.incomplete = true;
        }
    }
    cnt += ring.muls.get();
    sols.sort();
    sols.dedup();
    st.lin_muls = cnt;
    st.lin_ms = t_lin.elapsed().as_secs_f64() * 1e3;
    (sols, st)
}

/// Characteristic polynomial (monic, coefficients low to high) of a square
/// matrix over `F_p` by reduction to upper Hessenberg form and the standard
/// recurrence; multiplications counted into `cnt`.
pub fn charpoly_mod_p(m: &[Vec<u64>], p: u64, cnt: &mut u64) -> Vec<u64> {
    let n = m.len();
    let mut h: Vec<Vec<u64>> = m.to_vec();
    for mm_ in 1..n.saturating_sub(1) {
        // zero out h[i][mm_-1] for i > mm_ by similarity with the pivot h[mm_][mm_-1]
        if h[mm_][mm_ - 1] == 0 {
            let Some(i) = ((mm_ + 1)..n).find(|&i| h[i][mm_ - 1] != 0) else {
                continue;
            };
            h.swap(i, mm_);
            for row in h.iter_mut() {
                row.swap(i, mm_);
            }
        }
        let t = h[mm_][mm_ - 1];
        let tinv = inv_mod(t, p);
        for i in (mm_ + 1)..n {
            let u = mm(h[i][mm_ - 1], tinv, p);
            *cnt += 1;
            if u == 0 {
                continue;
            }
            for j in (mm_ - 1)..n {
                let v = mm(u, h[mm_][j], p);
                h[i][j] = sm(h[i][j], v, p);
            }
            *cnt += (n - mm_ + 1) as u64;
            for row in h.iter_mut() {
                let v = mm(u, row[i], p);
                row[mm_] = am(row[mm_], v, p);
            }
            *cnt += n as u64;
        }
    }
    // A_0 = 1; A_k = (X − h_kk) A_{k−1} − Σ_{i=1}^{k−1} h_{k−i,k} (Π_{j=k−i+1}^{k} h_{j,j−1}) A_{k−i−1}
    let mut a: Vec<Vec<u64>> = vec![vec![1]];
    for k in 1..=n {
        let hk = h[k - 1][k - 1];
        let prev = &a[k - 1];
        let mut cur = vec![0u64; k + 1];
        for (d, &c) in prev.iter().enumerate() {
            cur[d + 1] = am(cur[d + 1], c, p);
            cur[d] = sm(cur[d], mm(hk, c, p), p);
        }
        *cnt += prev.len() as u64;
        // product of sub-diagonal entries h_{j,j-1}, j = k−i+1..k (1-indexed)
        let mut sub = 1u64;
        for i in 1..k {
            // the new sub-diagonal factor h_{k−i+1, k−i} (1-indexed)
            sub = mm(sub, h[k - i][k - i - 1], p);
            *cnt += 1;
            let coef = mm(h[k - i - 1][k - 1], sub, p);
            *cnt += 1;
            if coef != 0 {
                let ap = &a[k - i - 1];
                for (d, &c) in ap.iter().enumerate() {
                    cur[d] = sm(cur[d], mm(coef, c, p), p);
                }
                *cnt += ap.len() as u64;
            }
        }
        a.push(cur);
    }
    a.pop().unwrap()
}

fn rref_mod_p_local(m: &mut [Vec<u64>], p: u64, muls: &mut u64) -> Vec<usize> {
    super::gaudry_cubic::rref_mod_p(m, p, muls)
}

// ── An instance, its factor base, and the decomposition test ──────────────

/// The data of an instance, free of the field contexts (which hold counters
/// and so cannot be shared across threads).
#[derive(Clone, Debug)]
pub struct Spec {
    pub p: u64,
    pub alpha: E6,
    /// The prime order of the subgroup the logarithm lives in.
    pub l: u64,
    pub g: PtE6,
    pub q: PtE6,
    pub d: u64,
    /// `Φ(G)` and `Φ(Q)` on `Jac_H(F_q)`.
    pub gj: Div<E2>,
    pub qj: Div<E2>,
    /// `F_p` multiplications spent building the instance's transfer.
    pub transfer_muls: u64,
}

/// Per-thread contexts: the tower, the cover, and the `Re`/`Im` tables.
pub struct Ctx {
    pub f: Fq3,
    pub cov: Cover,
    pub reim: ReIm,
    /// Replay F4's trace on the six-quadric systems (§12 of the note).
    pub trace: bool,
}

impl Ctx {
    pub fn new(spec: &Spec) -> Ctx {
        Self::with_trace(spec, false)
    }
    pub fn with_trace(spec: &Spec, trace: bool) -> Ctx {
        let f = Fq3::new(spec.p);
        let cov = Cover::new(&f, &spec.alpha).expect("cover");
        let reim = ReIm::new(spec.p, f.f.w);
        f.reset_muls();
        Ctx {
            f,
            cov,
            reim,
            trace,
        }
    }
    pub fn jac(&self) -> Hyp<'_, Fq> {
        Hyp {
            f: &self.f.f,
            h: self.cov.hx.clone(),
            g: 3,
        }
    }
}

pub fn generate_spec(p: u64, seed: u64) -> Spec {
    let f = Fq3::new(p);
    let mut rng = StdRng::seed_from_u64(seed ^ 0xC07E_4);
    loop {
        let alpha = f.random(&mut rng);
        if alpha.in_fq() {
            continue;
        }
        let ec = EllE::new(&f, &alpha);
        let Some(l) = group_order_4_prime(&ec, &mut rng) else {
            continue;
        };
        let Some(cov) = Cover::new(&f, &alpha) else {
            continue;
        };
        let g = loop {
            let g = ec.mul(&ec.random_point(&mut rng), 4);
            if !g.inf && ec.mul(&g, l).inf {
                break g;
            }
        };
        let d = rng.gen_range(1..l);
        let q = ec.mul(&g, d);
        f.reset_muls();
        let (Some(gj), Some(qj)) = (cov.transfer(&f, &g), cov.transfer(&f, &q)) else {
            continue;
        };
        let transfer_muls = f.muls();
        let jac = Hyp {
            f: &f.f,
            h: cov.hx.clone(),
            g: 3,
        };
        if jac.is_identity(&gj)
            || !jac.is_identity(&jac.mul(&gj, l as u128))
            || jac.mul(&gj, d as u128) != qj
        {
            continue;
        }
        return Spec {
            p,
            alpha,
            l,
            g,
            q,
            d,
            gj,
            qj,
            transfer_muls,
        };
    }
}

/// An instance on a given weak curve `y² = x(x − α)(x − σα)` with given
/// points: `G` of prime order `l`, `Q = [d]G` (`d` is the planted answer, used
/// only to report correctness; the route never reads it).  `None` when the
/// cover refuses `α` or a transfer degenerates.  §18 builds this from the
/// curve the walk reached and the transported points.
pub fn spec_from_curve(p: u64, alpha: E6, l: u64, g: PtE6, q: PtE6, d: u64) -> Option<Spec> {
    let f = Fq3::new(p);
    let cov = Cover::new(&f, &alpha)?;
    f.reset_muls();
    let gj = cov.transfer(&f, &g)?;
    let qj = cov.transfer(&f, &q)?;
    let transfer_muls = f.muls();
    let jac = Hyp {
        f: &f.f,
        h: cov.hx.clone(),
        g: 3,
    };
    if jac.is_identity(&gj) || !jac.is_identity(&jac.mul(&gj, l as u128)) {
        return None;
    }
    Some(Spec {
        p,
        alpha,
        l,
        g,
        q,
        d,
        gj,
        qj,
        transfer_muls,
    })
}

#[derive(Clone, Debug)]
pub struct BaseEl {
    pub x: u64,
    pub y: E2,
    pub d: Div<E2>,
}

/// One element per `x ∈ F_p` with `f_H(x)` a non-zero square in `F_q`: the
/// class of `(x, ±y)`, `y` the lexicographically smaller root.
pub fn factor_base(ctx: &Ctx, rng: &mut StdRng) -> Vec<BaseEl> {
    let f = &ctx.f.f;
    let p = f.p;
    let mut out = Vec::new();
    for x in 0..p {
        let xe = f.from_fp(x);
        let rhs = peval(f, &ctx.cov.hx, &xe);
        let Some(y) = f.sqrt(&rhs, rng) else { continue };
        if y == E2::ZERO {
            continue;
        }
        let ny = f.neg(&y);
        let y = if y.0 <= ny.0 { y } else { ny };
        out.push(BaseEl {
            x,
            y,
            d: Div {
                u: vec![f.neg(&xe), E2::ONE],
                v: vec![y],
            },
        });
    }
    out
}

/// `R = Σ ε_i B_i` over six distinct base elements, `ε_i = ±1`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub struct Dec6 {
    pub terms: [(usize, i8); 6],
}

pub fn verify_dec(jac: &Hyp<Fq>, base: &[BaseEl], r: &Div<E2>, d: &Dec6) -> bool {
    let mut acc = jac.identity();
    for &(i, e) in &d.terms {
        let b = if e == 1 {
            base[i].d.clone()
        } else {
            jac.neg(&base[i].d)
        };
        acc = jac.add(&acc, &b);
    }
    acc == *r
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct NagaoCost {
    pub weil_muls: u64,
    pub f4_muls: u64,
    pub lin_muls: u64,
    /// Roots of the sextic, the ordinates, the soundness checks.
    pub post_muls: u64,
    pub delta: usize,
    pub degree_reached: u32,
    pub max_rows: usize,
    pub max_cols: usize,
    pub f4_ms: f64,
    pub lin_ms: f64,
    pub total_ms: f64,
    pub inconsistent: bool,
    pub incomplete: bool,
    pub timed_out: bool,
    pub degenerate: bool,
    pub nonseparating: usize,
    pub stopped_at: Option<usize>,
    pub replayed: bool,
    pub trace_mismatch: bool,
    /// Solutions over `F_p` of the six quadrics, and those whose sextic
    /// splits into six distinct roots with admissible ordinates.
    pub fp_solutions: usize,
    pub split: usize,
    pub decompositions: usize,
}

impl NagaoCost {
    pub fn total(&self) -> u64 {
        self.weil_muls + self.f4_muls + self.lin_muls + self.post_muls
    }
}

pub fn jv_options(max_degree: u32, budget_secs: f64) -> F4Options {
    F4Options::new(F4Ordering::Grevlex, max_degree)
        .with_budget(std::time::Duration::from_secs_f64(budget_secs))
}

/// [`jv_options`] with F4 stopped as soon as the leading monomials found so
/// far leave at most `staircase` standard monomials.  Six quadrics in six
/// unknowns have Bézout degree `2⁶ = 64`: once a partial basis's staircase
/// is `64`, it generates a subideal `J ⊆ I` with `dim R/J ≤ 64 ≤ dim R/I`,
/// so `J = I` and the partial basis is a Gröbner basis of `I`; every
/// further F4 step would only certify that.  A system whose ideal has
/// degree below `64` stops at a strict superideal's staircase, which the
/// solver's final check of every candidate against the input covers.
pub fn jv_options_stopped(max_degree: u32, budget_secs: f64, staircase: usize) -> F4Options {
    jv_options(max_degree, budget_secs).stopping_below(staircase)
}

pub(crate) fn opts_for(max_degree: u32, budget_secs: f64, stop: Option<usize>) -> F4Options {
    match stop {
        Some(k) => jv_options_stopped(max_degree, budget_secs, k),
        None => jv_options(max_degree, budget_secs),
    }
}

/// Every six-point decomposition of the reduced divisor `r` over the base, by
/// the Weil-restricted Nagao system.
pub fn nagao_decompose(
    ctx: &Ctx,
    base: &[BaseEl],
    by_x: &HashMap<u64, usize>,
    r: &Div<E2>,
    opts: &F4Options,
    rng: &mut StdRng,
) -> (Vec<Dec6>, NagaoCost) {
    let f = &ctx.f.f;
    let p = f.p;
    let mut cost = NagaoCost::default();
    if r.u.len() != 4 {
        cost.degenerate = true;
        return (Vec::new(), cost);
    }
    let t_all = Instant::now();
    let m0 = f.muls();
    let mut v = r.v.clone();
    v.resize(3, E2::ZERO);
    let forms = weil_forms(f, &ctx.cov.hx, &r.u, &v);
    let polys = forms.polys(&ctx.reim, p);
    f.count_public(6 * 10 * 8);
    cost.weil_muls = f.muls() - m0;
    let (sols, st) = solve_zero_dim_traced(&polys, 6, p, opts, rng, ctx.trace);
    cost.f4_muls = st.f4_muls;
    cost.lin_muls = st.lin_muls;
    cost.delta = st.delta;
    cost.degree_reached = st.degree_reached;
    cost.max_rows = st.max_rows;
    cost.max_cols = st.max_cols;
    cost.f4_ms = st.f4_ms;
    cost.lin_ms = st.lin_ms;
    cost.inconsistent = st.inconsistent;
    cost.incomplete = st.incomplete;
    cost.timed_out = st.timed_out;
    cost.nonseparating = st.nonseparating;
    cost.stopped_at = st.stopped_at;
    cost.replayed = st.replayed;
    cost.trace_mismatch = st.trace_mismatch;
    cost.fp_solutions = sols.len();
    let m1 = f.muls();
    let ring = PolyRing::new(p);
    let mut out: BTreeSet<Dec6> = BTreeSet::new();
    let uc = |k: isize| -> E2 {
        if k >= 0 && (k as usize) < r.u.len() {
            r.u[k as usize]
        } else {
            E2::ZERO
        }
    };
    let vc = |k: isize| -> E2 {
        if k >= 0 && (k as usize) < v.len() {
            v[k as usize]
        } else {
            E2::ZERO
        }
    };
    'sol: for s in &sols {
        let lam = [E2([s[0], s[1]]), E2([s[2], s[3]]), E2([s[4], s[5]])];
        let mut coeffs = Vec::with_capacity(7);
        for k in 0..6 {
            let c = forms.eval(f, k, &lam);
            if c.0[1] != 0 {
                continue 'sol;
            }
            coeffs.push(c.0[0]);
        }
        coeffs.push(1);
        let roots = ring.roots(&UPoly(coeffs), rng);
        if roots.len() != 6 {
            continue;
        }
        // A_k = −v_{k−1} + λ₀u_k + λ₁u_{k−1} − µ₀v_k
        let ak: Vec<E2> = (0..5isize)
            .map(|k| {
                let mut a = f.neg(&vc(k - 1));
                a = f.add(&a, &f.mul(&lam[0], &uc(k)));
                a = f.add(&a, &f.mul(&lam[1], &uc(k - 1)));
                f.sub(&a, &f.mul(&lam[2], &vc(k)))
            })
            .collect();
        // Per root, the admissible signs: one when `y = −A/µ` is determined, both
        // when `µ(x) = 0` (then `A(x) = 0` too, and `u(x) = 0`: the reduced
        // divisor contains a point over this abscissa, `y` is not determined by
        // the formula and the group decides).
        let jc = ctx.jac();
        let mut options: Vec<(usize, Vec<i8>)> = Vec::with_capacity(6);
        for &x in &roots {
            let xe = f.from_fp(x);
            let Some(&i) = by_x.get(&x) else {
                continue 'sol;
            };
            let mu = f.add(&lam[2], &xe);
            if mu == E2::ZERO {
                options.push((i, vec![1, -1]));
                continue;
            }
            let y = f.neg(&f.mul(&peval(f, &ak, &xe), &f.inv(&mu)));
            if y == E2::ZERO || f.sq(&y) != peval(f, &ctx.cov.hx, &xe) {
                continue 'sol;
            }
            // R = −Σ (Q_i − ∞): the sign of the base element is the opposite
            let eps: i8 = if y == base[i].y {
                -1
            } else if y == f.neg(&base[i].y) {
                1
            } else {
                continue 'sol;
            };
            options.push((i, vec![eps]));
        }
        if options.len() != 6 {
            continue;
        }
        let free: Vec<usize> = (0..6).filter(|&k| options[k].1.len() > 1).collect();
        if free.len() > 4 {
            continue;
        }
        cost.split += 1;
        for combo in 0..(1u32 << free.len()) {
            let mut terms: Vec<(usize, i8)> = options.iter().map(|(i, e)| (*i, e[0])).collect();
            for (bit, &k) in free.iter().enumerate() {
                terms[k].1 = if combo & (1 << bit) == 0 { 1 } else { -1 };
            }
            terms.sort_unstable();
            let t: [(usize, i8); 6] = core::array::from_fn(|k| terms[k]);
            let d = Dec6 { terms: t };
            if free.is_empty() || verify_dec(&jc, base, r, &d) {
                out.insert(d);
            }
        }
    }
    cost.post_muls = f.muls() - m1 + ring.muls.get();
    cost.decompositions = out.len();
    cost.total_ms = t_all.elapsed().as_secs_f64() * 1e3;
    (out.into_iter().collect(), cost)
}

// ── The oracle: meet in the middle over the three-point sums ──────────────

type JKey = ([E2; 4], [E2; 3]);

fn jkey(d: &Div<E2>) -> JKey {
    let mut u = [E2::ZERO; 4];
    for (i, c) in d.u.iter().enumerate().take(4) {
        u[i] = *c;
    }
    let mut v = [E2::ZERO; 3];
    for (i, c) in d.v.iter().enumerate().take(3) {
        v[i] = *c;
    }
    (u, v)
}

pub struct Table3 {
    entries: Vec<(Div<E2>, [u32; 3], u8)>,
    map: HashMap<JKey, Vec<u32>>,
}

pub fn table3(jac: &Hyp<Fq>, base: &[BaseEl]) -> Table3 {
    let m = base.len();
    let neg: Vec<Div<E2>> = base.iter().map(|b| jac.neg(&b.d)).collect();
    let mut entries = Vec::new();
    let mut map: HashMap<JKey, Vec<u32>> = HashMap::new();
    for i in 0..m {
        for j in (i + 1)..m {
            for k in (j + 1)..m {
                for mask in 0..8u8 {
                    let pick = |t: usize, idx: usize| -> &Div<E2> {
                        if mask & (1 << t) == 0 {
                            &base[idx].d
                        } else {
                            &neg[idx]
                        }
                    };
                    let s = jac.add(&jac.add(pick(0, i), pick(1, j)), pick(2, k));
                    map.entry(jkey(&s)).or_default().push(entries.len() as u32);
                    entries.push((s, [i as u32, j as u32, k as u32], mask));
                }
            }
        }
    }
    Table3 { entries, map }
}

pub fn mitm6(jac: &Hyp<Fq>, t: &Table3, r: &Div<E2>) -> Vec<Dec6> {
    let mut out: BTreeSet<Dec6> = BTreeSet::new();
    for (s1, idx1, mask1) in &t.entries {
        let rem = jac.sub(r, s1);
        let Some(hits) = t.map.get(&jkey(&rem)) else {
            continue;
        };
        for &h in hits {
            let (_, idx2, mask2) = &t.entries[h as usize];
            let mut terms: Vec<(usize, i8)> = Vec::with_capacity(6);
            for (idx, mask) in [(idx1, mask1), (idx2, mask2)] {
                for (tt, &i) in idx.iter().enumerate() {
                    terms.push((i as usize, if mask & (1 << tt) == 0 { 1 } else { -1 }));
                }
            }
            terms.sort_unstable();
            if terms.windows(2).any(|w| w[0].0 == w[1].0) {
                continue;
            }
            out.insert(Dec6 {
                terms: core::array::from_fn(|k| terms[k]),
            });
        }
    }
    out.into_iter().collect()
}

// ── Pollard rho on the same group E(F_{q³}) ───────────────────────────────

#[derive(Clone, Debug, Default, Serialize)]
pub struct RhoRunE {
    pub seed: u64,
    pub steps: u64,
    pub group_ops: u64,
    pub s: f64,
    pub correct: bool,
}

/// r-adding Pollard rho with exhaustive storage on `⟨G, Q⟩ ⊂ E(F_{q³})`,
/// charging every group operation including the walk's setup, as
/// `gaudry_quartic::rho4` does.
pub fn rho_e(spec: &Spec, seed: u64) -> RhoRunE {
    let f = Fq3::new(spec.p);
    let ec = EllE::new(&f, &spec.alpha);
    let n = spec.l;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x8D06);
    let r = 32usize;
    let mults: Vec<(u64, u64, PtE6)> = (0..r)
        .map(|_| {
            let al = rng.gen_range(0..n);
            let be = rng.gen_range(1..n);
            (al, be, ec.add(&ec.mul(&spec.g, al), &ec.mul(&spec.q, be)))
        })
        .collect();
    let mut table: HashMap<PtE6, (u64, u64)> = HashMap::new();
    let mut steps = 0u64;
    let mut found = None;
    'outer: loop {
        let mut a = rng.gen_range(0..n);
        let mut b = rng.gen_range(1..n);
        let mut l = ec.add(&ec.mul(&spec.g, a), &ec.mul(&spec.q, b));
        loop {
            steps += 1;
            if let Some(&(a2, b2)) = table.get(&l) {
                let db = (b + n - b2) % n;
                if db != 0 {
                    let da = (a2 + n - a) % n;
                    found = Some(mm(da, inv_mod(db, n), n));
                    break 'outer;
                }
                break;
            }
            table.insert(l, (a, b));
            let h = l.x.0[0].0[0]
                ^ l.x.0[0].0[1].wrapping_mul(0x9E37_79B9)
                ^ l.x.0[1].0[0].rotate_left(17)
                ^ l.x.0[2].0[1].rotate_left(33);
            let j = (mix64(h) % r as u64) as usize;
            l = ec.add(&l, &mults[j].2);
            a = (a + mults[j].0) % n;
            b = (b + mults[j].1) % n;
            if steps > 64 * (n as f64).sqrt() as u64 + 1_000_000 {
                break 'outer;
            }
        }
    }
    let ops = ec.ops();
    RhoRunE {
        seed,
        steps,
        group_ops: ops,
        s: ops as f64 / (n as f64).sqrt(),
        correct: found == Some(spec.d),
    }
}

/// [`rho_e`] without the stored walk (note §15): the same r-adding walk,
/// `threads` walkers from random starts, a point distinguished when `dp_bits`
/// bits of its mixed hash vanish, the distinguished points shared between
/// the walkers; a walker that meets no distinguished point in `40·2^dp_bits`
/// steps restarts from a new random point.  Every group operation of every
/// walker is charged, the walk's multipliers once per walker.
pub fn rho_e_dp(spec: &Spec, seed: u64, dp_bits: u32, threads: usize) -> RhoRunE {
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
    use std::sync::Mutex;
    let n = spec.l;
    let r = 32usize;
    let scalars: Vec<(u64, u64)> = {
        let mut rng = StdRng::seed_from_u64(seed ^ 0x8D06);
        (0..r)
            .map(|_| (rng.gen_range(0..n), rng.gen_range(1..n)))
            .collect()
    };
    let table: Mutex<HashMap<PtE6, (u64, u64)>> = Mutex::new(HashMap::new());
    let found: Mutex<Option<u64>> = Mutex::new(None);
    let stop = AtomicBool::new(false);
    let steps = AtomicU64::new(0);
    let ops = AtomicU64::new(0);
    let cap = 64 * (n as f64).sqrt() as u64 + 1_000_000;
    let mask = (1u64 << dp_bits) - 1;
    let tail_cap = 40u64 << dp_bits;
    let progress = std::env::var("JV_RHO_PROGRESS").is_ok();
    let start = Instant::now();
    std::thread::scope(|sc| {
        for w in 0..threads.max(1) {
            let (table, found, stop, steps, ops, scalars) =
                (&table, &found, &stop, &steps, &ops, &scalars);
            sc.spawn(move || {
                let f = Fq3::new(spec.p);
                let ec = EllE::new(&f, &spec.alpha);
                let mults: Vec<PtE6> = scalars
                    .iter()
                    .map(|&(al, be)| ec.add(&ec.mul(&spec.g, al), &ec.mul(&spec.q, be)))
                    .collect();
                let mut rng = StdRng::seed_from_u64(
                    seed ^ 0x8D06 ^ (w as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15),
                );
                let mut local = 0u64;
                'outer: while !stop.load(Ordering::Relaxed) {
                    let mut a = rng.gen_range(0..n);
                    let mut b = rng.gen_range(1..n);
                    let mut l = ec.add(&ec.mul(&spec.g, a), &ec.mul(&spec.q, b));
                    let mut tail = 0u64;
                    loop {
                        let h = l.x.0[0].0[0]
                            ^ l.x.0[0].0[1].wrapping_mul(0x9E37_79B9)
                            ^ l.x.0[1].0[0].rotate_left(17)
                            ^ l.x.0[2].0[1].rotate_left(33);
                        let idx = mix64(h);
                        if (idx >> 20) & mask == 0 {
                            let mut t = table.lock().unwrap();
                            if let Some(&(a2, b2)) = t.get(&l) {
                                let db = (b + n - b2) % n;
                                if db != 0 {
                                    let da = (a2 + n - a) % n;
                                    *found.lock().unwrap() = Some(mm(da, inv_mod(db, n), n));
                                    stop.store(true, Ordering::Relaxed);
                                    break 'outer;
                                }
                                drop(t);
                                break; // the same walk met its own point: restart
                            }
                            t.insert(l, (a, b));
                            drop(t);
                            tail = 0;
                        }
                        let j = (idx % r as u64) as usize;
                        l = ec.add(&l, &mults[j]);
                        a = (a + scalars[j].0) % n;
                        b = (b + scalars[j].1) % n;
                        local += 1;
                        tail += 1;
                        if local & ((1 << 20) - 1) == 0 {
                            let tot = steps.fetch_add(1 << 20, Ordering::Relaxed) + (1 << 20);
                            if progress && w == 0 && tot & ((1 << 26) - 1) < (1 << 20) {
                                eprintln!(
                                    "  [rho-dp p={} seed={}] steps {:.3e} ({:.2} sqrt l) dps {} {:.0} s",
                                    spec.p,
                                    seed,
                                    tot as f64,
                                    tot as f64 / (n as f64).sqrt(),
                                    table.lock().unwrap().len(),
                                    start.elapsed().as_secs_f64()
                                );
                            }
                            if tot > cap || stop.load(Ordering::Relaxed) {
                                stop.store(true, Ordering::Relaxed);
                                break 'outer;
                            }
                        }
                        if tail > tail_cap {
                            break;
                        }
                    }
                }
                steps.fetch_add(local & ((1 << 20) - 1), Ordering::Relaxed);
                ops.fetch_add(ec.ops(), Ordering::Relaxed);
            });
        }
    });
    let ops = ops.load(Ordering::Relaxed);
    let correct = *found.lock().unwrap() == Some(spec.d);
    RhoRunE {
        seed,
        steps: steps.load(Ordering::Relaxed),
        group_ops: ops,
        s: ops as f64 / (n as f64).sqrt(),
        correct,
    }
}

/// One distinguished-point rho run on the instance `(p, seed)` (note §15).
#[derive(Clone, Debug, Default, Serialize)]
pub struct RhoDpReport {
    pub p: u64,
    pub seed: u64,
    pub run: u64,
    pub l: u64,
    pub bits: f64,
    pub dp_bits: u32,
    pub threads: usize,
    pub c_add_e: f64,
    pub steps: u64,
    pub group_ops: u64,
    pub s: f64,
    pub correct: bool,
    pub wall_ms: f64,
}

pub fn run_rho_dp(p: u64, seed: u64, run: u64, dp_bits: u32, threads: usize) -> RhoDpReport {
    let spec = generate_spec(p, seed);
    let ctx = Ctx::new(&spec);
    let (c_e, _) = unit_costs(&spec, &ctx);
    let start = Instant::now();
    let r = rho_e_dp(&spec, seed.wrapping_mul(1000) + run, dp_bits, threads);
    RhoDpReport {
        p,
        seed,
        run,
        l: spec.l,
        bits: (spec.l as f64).log2(),
        dp_bits,
        threads,
        c_add_e: c_e,
        steps: r.steps,
        group_ops: r.group_ops,
        s: r.s,
        correct: r.correct,
        wall_ms: start.elapsed().as_secs_f64() * 1e3,
    }
}

/// `F_p` multiplications per affine addition on `E(F_{q³})` and per addition
/// on `Jac_H(F_q)`, measured.
pub fn unit_costs(spec: &Spec, ctx: &Ctx) -> (f64, f64) {
    let ec = EllE::new(&ctx.f, &spec.alpha);
    ctx.f.reset_muls();
    let mut acc = spec.g;
    for _ in 0..64 {
        acc = ec.add(&acc, &spec.q);
    }
    let c_e = ctx.f.muls() as f64 / 64.0;
    let jac = ctx.jac();
    ctx.f.reset_muls();
    let mut jacc = spec.gj.clone();
    for _ in 0..64 {
        jacc = jac.add(&jacc, &spec.qj);
    }
    let c_j = ctx.f.muls() as f64 / 64.0;
    ctx.f.reset_muls();
    (c_e, c_j)
}

// ── Experiment: C_cov, the cost of one decomposition test ─────────────────

#[derive(Clone, Debug, Default, Serialize)]
pub struct CostStatsC {
    pub mean: f64,
    pub min: u64,
    pub max: u64,
}

fn stats_c(v: &[u64]) -> CostStatsC {
    if v.is_empty() {
        return CostStatsC::default();
    }
    CostStatsC {
        mean: v.iter().sum::<u64>() as f64 / v.len() as f64,
        min: *v.iter().min().unwrap(),
        max: *v.iter().max().unwrap(),
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct CcovReport {
    pub p: u64,
    pub seed: u64,
    pub l: u64,
    pub bits: f64,
    pub base: usize,
    pub c_add_e: f64,
    pub c_add_j: f64,
    pub max_degree: u32,
    pub budget_secs: f64,
    /// F4 stopped at this staircase, or run to a certified basis.
    pub stop_staircase: Option<usize>,
    pub stopped: usize,
    /// F4's trace replayed (the first system of a `p` records it).
    pub trace: bool,
    pub replayed: usize,
    pub trace_mismatches: usize,
    pub random_residuals: usize,
    pub constructed_residuals: usize,
    pub random_decomposable: usize,
    pub expected_rate: f64,
    pub planted_found: usize,
    pub oracle_checked: usize,
    pub mismatches: usize,
    pub unverified: usize,
    pub incomplete: usize,
    pub timed_out: usize,
    pub nonseparating: usize,
    pub c_cov: CostStatsC,
    pub weil_muls: CostStatsC,
    pub f4_muls: CostStatsC,
    pub lin_muls: CostStatsC,
    pub post_muls: CostStatsC,
    pub delta: CostStatsC,
    pub degree_reached: CostStatsC,
    pub max_rows: CostStatsC,
    pub max_cols: CostStatsC,
    pub f4_ms: CostStatsC,
    pub lin_ms: CostStatsC,
    pub total_ms: CostStatsC,
    pub fp_solutions: CostStatsC,
    pub c_cov_decomposable: CostStatsC,
    pub per_test: Vec<NagaoCost>,
    pub wall_ms: f64,
}

/// A random six-sum and the decomposition it plants.
fn planted(jac: &Hyp<Fq>, base: &[BaseEl], rng: &mut StdRng) -> (Div<E2>, Dec6) {
    let mut idx: Vec<usize> = Vec::new();
    while idx.len() < 6 {
        let i = rng.gen_range(0..base.len());
        if !idx.contains(&i) {
            idx.push(i);
        }
    }
    let mut r = jac.identity();
    let mut terms: Vec<(usize, i8)> = Vec::new();
    for &i in &idx {
        let e: i8 = if rng.gen_bool(0.5) { 1 } else { -1 };
        let b = if e == 1 {
            base[i].d.clone()
        } else {
            jac.neg(&base[i].d)
        };
        r = jac.add(&r, &b);
        terms.push((i, e));
    }
    terms.sort_unstable();
    (
        r,
        Dec6 {
            terms: core::array::from_fn(|k| terms[k]),
        },
    )
}

#[allow(clippy::too_many_arguments)]
pub fn run_cover_ccov(
    p: u64,
    seed: u64,
    random: usize,
    constructed: usize,
    max_degree: u32,
    budget_secs: f64,
    oracle_up_to_base: usize,
    stop: Option<usize>,
    trace: bool,
) -> CcovReport {
    let start = Instant::now();
    let spec = generate_spec(p, seed);
    let ctx = Ctx::with_trace(&spec, trace);
    let mut rng = StdRng::seed_from_u64(seed ^ 0xC0C0);
    let base = factor_base(&ctx, &mut rng);
    let by_x: HashMap<u64, usize> = base.iter().enumerate().map(|(i, b)| (b.x, i)).collect();
    let (c_e, c_j) = unit_costs(&spec, &ctx);
    let jac = ctx.jac();
    let table = (base.len() <= oracle_up_to_base).then(|| table3(&jac, &base));
    let mut rep = CcovReport {
        p,
        seed,
        l: spec.l,
        bits: (spec.l as f64).log2(),
        base: base.len(),
        c_add_e: c_e,
        c_add_j: c_j,
        max_degree,
        budget_secs,
        stop_staircase: stop,
        trace,
        random_residuals: random,
        constructed_residuals: constructed,
        expected_rate: 1.0 / 720.0,
        ..Default::default()
    };
    let mut jobs: Vec<(Div<E2>, Option<Dec6>)> = Vec::new();
    for _ in 0..random {
        let a = rng.gen_range(0..spec.l);
        let b = rng.gen_range(1..spec.l);
        jobs.push((
            jac.add(&jac.mul(&spec.gj, a as u128), &jac.mul(&spec.qj, b as u128)),
            None,
        ));
    }
    for _ in 0..constructed {
        let (r, w) = planted(&jac, &base, &mut rng);
        jobs.push((r, Some(w)));
    }
    // Sequential: F4's per-run operation counter is a difference of a
    // process-wide counter and interleaves under concurrency.
    let results: Vec<(Vec<Dec6>, NagaoCost, bool, bool, Option<bool>)> = jobs
        .iter()
        .enumerate()
        .map(|(k, (r, want))| {
            let mut lrng =
                StdRng::seed_from_u64(seed ^ (k as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            let opts = opts_for(max_degree, budget_secs, stop);
            let (found, cost) = nagao_decompose(&ctx, &base, &by_x, r, &opts, &mut lrng);
            let verified = found.iter().all(|d| verify_dec(&jac, &base, r, d));
            let planted_ok = want.as_ref().map(|w| found.contains(w));
            let oracle_ok = table.as_ref().map(|t| {
                let o = mitm6(&jac, t, r);
                cost.incomplete || cost.timed_out || o == found
            });
            (found, cost, verified, oracle_ok.unwrap_or(true), planted_ok)
        })
        .collect();
    let mut totals = Vec::new();
    let mut totals_dec = Vec::new();
    let mut costs = Vec::new();
    for (k, (found, cost, verified, oracle_ok, planted_ok)) in results.into_iter().enumerate() {
        if !verified {
            rep.unverified += 1;
        }
        if table.is_some() {
            rep.oracle_checked += 1;
            if !oracle_ok {
                rep.mismatches += 1;
            }
        }
        if cost.incomplete {
            rep.incomplete += 1;
        }
        if cost.timed_out {
            rep.timed_out += 1;
        }
        rep.nonseparating += cost.nonseparating;
        if cost.stopped_at.is_some() {
            rep.stopped += 1;
        }
        if cost.replayed {
            rep.replayed += 1;
        }
        if cost.trace_mismatch {
            rep.trace_mismatches += 1;
        }
        if k < random {
            if !found.is_empty() {
                rep.random_decomposable += 1;
            }
        } else if planted_ok == Some(true) {
            rep.planted_found += 1;
        }
        if !cost.incomplete && !cost.timed_out && !cost.degenerate {
            totals.push(cost.total());
            if !found.is_empty() {
                totals_dec.push(cost.total());
            }
            costs.push(cost);
        }
    }
    let col = |g: &dyn Fn(&NagaoCost) -> u64| stats_c(&costs.iter().map(g).collect::<Vec<u64>>());
    rep.c_cov = stats_c(&totals);
    rep.c_cov_decomposable = stats_c(&totals_dec);
    rep.weil_muls = col(&|c| c.weil_muls);
    rep.f4_muls = col(&|c| c.f4_muls);
    rep.lin_muls = col(&|c| c.lin_muls);
    rep.post_muls = col(&|c| c.post_muls);
    rep.delta = col(&|c| c.delta as u64);
    rep.degree_reached = col(&|c| c.degree_reached as u64);
    rep.max_rows = col(&|c| c.max_rows as u64);
    rep.max_cols = col(&|c| c.max_cols as u64);
    rep.f4_ms = col(&|c| c.f4_ms as u64);
    rep.lin_ms = col(&|c| c.lin_ms as u64);
    rep.total_ms = col(&|c| c.total_ms as u64);
    rep.fp_solutions = col(&|c| c.fp_solutions as u64);
    rep.per_test = costs;
    rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
    rep
}

// ── The method end to end ─────────────────────────────────────────────────

#[derive(Clone, Debug, Default, Serialize)]
pub struct CoverDlpReport {
    pub p: u64,
    pub seed: u64,
    pub l: u64,
    pub bits: f64,
    pub base: usize,
    pub c_add_e: f64,
    pub c_add_j: f64,
    pub stop_staircase: Option<usize>,
    pub stopped: u64,
    pub trace: bool,
    pub replayed: u64,
    pub trace_mismatches: u64,
    pub residuals: u64,
    pub decompositions: u64,
    pub decomposition_rate: f64,
    pub expected_rate: f64,
    pub relations: usize,
    pub floor_residuals: f64,
    pub residual_ratio: f64,
    pub incomplete: u64,
    pub timed_out: u64,
    pub cross_checked: u64,
    pub mismatches: u64,
    /// Phases, in `F_p` multiplications.
    pub setup_muls: u64,
    pub base_muls: u64,
    pub stream_muls: u64,
    pub weil_muls: u64,
    pub f4_muls: u64,
    pub lin_muls: u64,
    pub post_muls: u64,
    pub verify_muls: u64,
    pub oracle_muls: u64,
    /// `C_cov` as the run paid it, per residual.
    pub c_cov: f64,
    pub unknowns: usize,
    pub filtered_out: usize,
    pub phi: f64,
    pub row_weight: f64,
    pub la_attempts: u64,
    pub la_ops: u64,
    pub la_muls: u64,
    pub wiedemann_constant: f64,
    pub la_ops_last: u64,
    pub total_muls: u64,
    pub s: f64,
    pub solved: bool,
    pub correct: bool,
    pub rho: Vec<RhoRunE>,
    pub rho_s_mean: f64,
    pub rho_s_sd: f64,
    pub s_over_rho: f64,
    pub relation_over_rho: f64,
    pub r: f64,
    pub r_last: f64,
    /// `360p·C_cov / (ρ_S · (√ℓ) · c_add) + r`, the registered formula with
    /// this run's own constants.
    pub predicted_s_over_rho: f64,
    pub wall_ms: f64,
}

/// `rho_runs = 0` takes rho's `S` from `rho_s_ref` (rho's table does not fit in
/// memory above `ℓ ≈ 2^48`); the report then carries no rho runs.
pub fn run_cover_dlp(
    p: u64,
    seed: u64,
    rho_runs: usize,
    check_every: u64,
    rho_s_ref: f64,
    stop: Option<usize>,
    trace: bool,
) -> CoverDlpReport {
    let start = Instant::now();
    let spec = generate_spec(p, seed);
    let ctx = Ctx::with_trace(&spec, trace);
    let jac = ctx.jac();
    let l = spec.l;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x3D1F6);
    let (c_e, c_j) = unit_costs(&spec, &ctx);
    ctx.f.reset_muls();
    let base = factor_base(&ctx, &mut rng);
    let base_muls = ctx.f.muls();
    let by_x: HashMap<u64, usize> = base.iter().enumerate().map(|(i, b)| (b.x, i)).collect();
    let small = base.len();
    let unknowns = small + 1;
    let table = (check_every > 0 && small <= 40).then(|| table3(&jac, &base));
    ctx.f.reset_muls();
    // The residual stream R_i = R_0 + i·M over ⟨G′, Q′⟩, one Jacobian
    // addition each; its start and step are setup.
    let (al, be) = (rng.gen_range(1..l), rng.gen_range(1..l));
    let step = jac.add(
        &jac.mul(&spec.gj, al as u128),
        &jac.mul(&spec.qj, be as u128),
    );
    let mut a = rng.gen_range(0..l);
    let mut b = rng.gen_range(1..l);
    let mut r = jac.add(&jac.mul(&spec.gj, a as u128), &jac.mul(&spec.qj, b as u128));
    let setup_muls = ctx.f.muls() + spec.transfer_muls;
    ctx.f.reset_muls();
    let mut rep = CoverDlpReport {
        p,
        seed,
        l,
        bits: (l as f64).log2(),
        base: small,
        c_add_e: c_e,
        c_add_j: c_j,
        stop_staircase: stop,
        trace,
        expected_rate: 1.0 / 720.0,
        floor_residuals: unknowns as f64 * 720.0,
        base_muls,
        setup_muls,
        ..Default::default()
    };
    let mut full_rels: Vec<SparseRel> = Vec::new();
    let mut next_attempt = unknowns;
    let mut la_ops = 0u64;
    let cap = 400 * (unknowns as u64) * 720;
    'collect: while rep.residuals < cap {
        let mut batch: Vec<(u64, u64, Div<E2>, u64)> = Vec::with_capacity(64);
        for _ in 0..64 {
            let idx = rep.residuals + batch.len() as u64;
            batch.push((a, b, r.clone(), idx));
            r = jac.add(&r, &step);
            a = (a + al) % l;
            b = (b + be) % l;
        }
        rep.stream_muls += ctx.f.muls();
        ctx.f.reset_muls();
        let f4_before = f4_fp::field_ops_total();
        let found: Vec<(u64, u64, Div<E2>, Vec<Dec6>, NagaoCost, u64, bool, bool)> = batch
            .par_iter()
            .map_init(
                || Ctx::with_trace(&spec, trace),
                |c, (a, b, r, idx)| {
                    // the deadline is absolute: one budget per test, not per run
                    let opts = opts_for(24, 600.0, stop);
                    let mut lrng =
                        StdRng::seed_from_u64(seed ^ idx.wrapping_mul(0x9E37_79B9_7F4A_7C15));
                    let (decs, cost) = nagao_decompose(c, &base, &by_x, r, &opts, &mut lrng);
                    let jc = c.jac();
                    let v0 = c.f.muls();
                    let ok: Vec<Dec6> = decs
                        .into_iter()
                        .filter(|d| verify_dec(&jc, &base, r, d))
                        .collect();
                    let vm = c.f.muls() - v0;
                    let (checked, was_checked) =
                        match (&table, check_every > 0 && idx % check_every == 0) {
                            (Some(t), true) => {
                                let o = mitm6(&jc, t, r);
                                (cost.incomplete || cost.timed_out || o == ok, true)
                            }
                            _ => (true, false),
                        };
                    (*a, *b, r.clone(), ok, cost, vm, checked, was_checked)
                },
            )
            .collect();
        // the F4 work is counted once, globally, by the batch's difference
        rep.f4_muls += f4_fp::field_ops_total() - f4_before;
        for (a, b, _r, decs, cost, vm, checked, was_checked) in found {
            rep.residuals += 1;
            rep.weil_muls += cost.weil_muls;
            rep.lin_muls += cost.lin_muls;
            rep.post_muls += cost.post_muls;
            rep.verify_muls += vm;
            if cost.incomplete {
                rep.incomplete += 1;
            }
            if cost.timed_out {
                rep.timed_out += 1;
            }
            if cost.stopped_at.is_some() {
                rep.stopped += 1;
            }
            if cost.replayed {
                rep.replayed += 1;
            }
            if cost.trace_mismatch {
                rep.trace_mismatches += 1;
            }
            if was_checked {
                rep.cross_checked += 1;
            }
            if !checked {
                rep.mismatches += 1;
            }
            let mut any = false;
            for d in decs {
                any = true;
                rep.decompositions += 1;
                let mut cols: Vec<(usize, u64)> = d
                    .terms
                    .iter()
                    .map(|&(i, e)| (i, if e == 1 { 1 } else { l - 1 }))
                    .collect();
                cols.push((small, (l - b) % l));
                cols.sort_unstable();
                full_rels.push(SparseRel { cols, rhs: a });
            }
            if any {
                rep.relations += 1;
            }
            if full_rels.len() >= next_attempt {
                if let Some(core) = square_core(&full_rels, small) {
                    rep.la_attempts += 1;
                    let mut map = vec![usize::MAX; unknowns];
                    let mut k = 0;
                    for c in 0..unknowns {
                        if core.columns.get(c).copied().unwrap_or(false) {
                            map[c] = k;
                            k += 1;
                        }
                    }
                    let sel: Vec<SparseRel> = core
                        .rows
                        .iter()
                        .map(|&i| SparseRel {
                            cols: full_rels[i]
                                .cols
                                .iter()
                                .map(|&(c, v)| (map[c], v))
                                .collect(),
                            rhs: full_rels[i].rhs,
                        })
                        .collect();
                    let before = la_ops;
                    let x = wiedemann_u64(&sel, sel.len(), l, &mut rng, &mut la_ops);
                    let dd = x.map(|x| x[map[small]]);
                    let ec = EllE::new(&ctx.f, &spec.alpha);
                    if dd.is_some_and(|dd| ec.mul(&spec.g, dd) == spec.q) {
                        rep.solved = true;
                        rep.correct = dd == Some(spec.d);
                        rep.unknowns = sel.len();
                        rep.filtered_out = full_rels.len() - sel.len();
                        rep.la_ops_last = la_ops - before;
                        break 'collect;
                    }
                    next_attempt = full_rels.len() + (sel.len() / 20).max(1);
                } else {
                    next_attempt = full_rels.len() + (unknowns / 20).max(1);
                }
            }
        }
    }
    rep.stream_muls += ctx.f.muls();
    rep.oracle_muls = rep.weil_muls + rep.f4_muls + rep.lin_muls + rep.post_muls;
    rep.c_cov = rep.oracle_muls as f64 / rep.residuals.max(1) as f64;
    rep.la_ops = la_ops;
    rep.la_muls = la_ops * 16;
    rep.decomposition_rate = rep.relations as f64 / rep.residuals.max(1) as f64;
    rep.residual_ratio = rep.residuals as f64 / rep.floor_residuals;
    rep.phi = rep.unknowns as f64 / small as f64;
    rep.row_weight = full_rels.iter().map(|r| r.cols.len()).sum::<usize>() as f64
        / full_rels.len().max(1) as f64;
    rep.wiedemann_constant = la_ops as f64 / (rep.unknowns.max(1) as f64).powi(2);
    rep.total_muls = rep.setup_muls
        + rep.base_muls
        + rep.stream_muls
        + rep.oracle_muls
        + rep.verify_muls
        + rep.la_muls;
    let sqrt_l = (l as f64).sqrt();
    rep.s = rep.total_muls as f64 / (c_e * sqrt_l);
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
    rep.relation_over_rho =
        (rep.setup_muls + rep.base_muls + rep.stream_muls + rep.oracle_muls + rep.verify_muls)
            as f64
            / rho_muls;
    rep.r = rep.la_muls as f64 / rho_muls;
    rep.r_last = rep.la_ops_last as f64 * 16.0 / rho_muls;
    rep.predicted_s_over_rho =
        360.0 * p as f64 * rep.c_cov / (rep.rho_s_mean * sqrt_l * c_e) + rep.r;
    rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
    rep
}

#[cfg(test)]
mod tests {
    use super::*;

    fn weak(p: u64, seed: u64) -> (Fq3, E6, StdRng) {
        let f = Fq3::new(p);
        let mut rng = StdRng::seed_from_u64(seed);
        let alpha = loop {
            let a = f.random(&mut rng);
            if !a.in_fq() {
                break a;
            }
        };
        (f, alpha, rng)
    }

    #[test]
    fn tower_arithmetic() {
        let (f, _, mut rng) = weak(53, 1);
        for _ in 0..20 {
            let a = f.random(&mut rng);
            let b = f.random(&mut rng);
            let c = f.random(&mut rng);
            assert_eq!(
                f.mul(&a, &f.add(&b, &c)),
                f.add(&f.mul(&a, &b), &f.mul(&a, &c))
            );
            assert_eq!(f.mul(&f.mul(&a, &b), &c), f.mul(&a, &f.mul(&b, &c)));
            if a != E6::ZERO {
                assert_eq!(f.mul(&a, &f.inv(&a)), E6::ONE);
            }
            // σ is a field automorphism of order 3 fixing F_q
            assert_eq!(f.sigma(&f.mul(&a, &b)), f.mul(&f.sigma(&a), &f.sigma(&b)));
            assert_eq!(f.sigma(&f.sigma(&f.sigma(&a))), a);
            // σ is the q-power: check on θ
        }
        let theta = E6([E2::ZERO, E2::ONE, E2::ZERO]);
        let q = (53u128) * 53;
        assert_eq!(f.pow(&theta, q), f.sigma(&theta));
        // F_{p²}
        for _ in 0..20 {
            let a = f.f.random(&mut rng);
            if a != E2::ZERO {
                assert_eq!(f.f.mul(&a, &f.f.inv(&a)), E2::ONE);
            }
        }
        // square roots
        let a = f.random(&mut rng);
        let s = f.sq(&a);
        let r = f.sqrt(&s, &mut rng).unwrap();
        assert_eq!(f.sq(&r), s);
    }

    #[test]
    fn polynomial_gcd_and_division() {
        let (f, _, mut rng) = weak(53, 2);
        let fq = &f.f;
        let rp = |rng: &mut StdRng, d: usize| -> Poly<E2> {
            let mut p: Poly<E2> = (0..=d).map(|_| fq.random(rng)).collect();
            ptrim(fq, &mut p);
            p
        };
        let a = rp(&mut rng, 6);
        let b = rp(&mut rng, 4);
        let (q, r) = pdivrem(fq, &a, &b);
        assert_eq!(padd(fq, &pmul(fq, &q, &b), &r), a);
        let c = rp(&mut rng, 2);
        let (g, s, t) = pxgcd(fq, &pmul(fq, &a, &c), &pmul(fq, &b, &c));
        assert_eq!(
            padd(
                fq,
                &pmul(fq, &s, &pmul(fq, &a, &c)),
                &pmul(fq, &t, &pmul(fq, &b, &c))
            ),
            g
        );
        assert!(g.len() >= c.len());
    }

    #[test]
    fn cover_map_lands_on_the_curve() {
        // Y² = X(X − α)(X − σα) on H, tested through y² = F·N at random x.
        let (f, alpha, mut rng) = weak(53, 3);
        let cov = Cover::new(&f, &alpha).expect("cover");
        let ec = EllE::new(&f, &alpha);
        let n6 = cov.embed(&f, &cov.nx);
        let h6 = cov.embed(&f, &cov.hx);
        let mut checked = 0;
        for _ in 0..30 {
            let x = f.random(&mut rng);
            let y2 = peval(&f, &h6, &x);
            let Some(y) = f.sqrt(&y2, &mut rng) else {
                continue;
            };
            let Some((bx, by)) = cov.map_point(&f, &x, &y) else {
                continue;
            };
            assert!(
                ec.on_curve(&PtE6 {
                    x: bx,
                    y: by,
                    inf: false
                }),
                "π(x, y) off E"
            );
            assert!(peval(&f, &n6, &x) != E6::ZERO);
            checked += 1;
        }
        assert!(checked >= 5, "{checked}");
    }

    #[test]
    fn jacobian_group_law() {
        let (f, alpha, mut rng) = weak(53, 4);
        let cov = Cover::new(&f, &alpha).expect("cover");
        let jac = Hyp {
            f: &f.f,
            h: cov.hx.clone(),
            g: 3,
        };
        let fq = &f.f;
        // random reduced divisors from three random points of H(F_q)
        let rand_div = |rng: &mut StdRng| -> Div<E2> {
            let mut d = jac.identity();
            let mut k = 0;
            while k < 3 {
                let x = fq.random(rng);
                if let Some(y) = fq.sqrt(&peval(fq, &cov.hx, &x), rng) {
                    if y != E2::ZERO {
                        let pt = Div {
                            u: vec![fq.neg(&x), E2::ONE],
                            v: vec![y],
                        };
                        assert!(jac.is_valid(&pt));
                        d = jac.add(&d, &pt);
                        k += 1;
                    }
                }
            }
            d
        };
        for _ in 0..6 {
            let a = rand_div(&mut rng);
            let b = rand_div(&mut rng);
            let c = rand_div(&mut rng);
            assert!(jac.is_valid(&a) && jac.is_valid(&b));
            assert_eq!(jac.add(&a, &b), jac.add(&b, &a));
            assert_eq!(jac.add(&jac.add(&a, &b), &c), jac.add(&a, &jac.add(&b, &c)));
            assert!(jac.is_identity(&jac.sub(&a, &a)));
            assert_eq!(jac.mul(&a, 5), jac.add(&jac.mul(&a, 2), &jac.mul(&a, 3)));
        }
    }

    #[test]
    fn transfer_is_a_homomorphism_into_the_jacobian_over_fq() {
        let (f, alpha, mut rng) = weak(53, 5);
        let cov = Cover::new(&f, &alpha).expect("cover");
        let ec = EllE::new(&f, &alpha);
        let jacq = Hyp {
            f: &f.f,
            h: cov.hx.clone(),
            g: 3,
        };
        let mut done = 0;
        while done < 4 {
            let p1 = ec.random_point(&mut rng);
            let p2 = ec.random_point(&mut rng);
            let p3 = ec.add(&p1, &p2);
            let (Some(a), Some(b), Some(c)) = (
                cov.transfer(&f, &p1),
                cov.transfer(&f, &p2),
                cov.transfer(&f, &p3),
            ) else {
                continue;
            };
            assert!(jacq.is_valid(&a) && jacq.is_valid(&b) && jacq.is_valid(&c));
            assert_eq!(jacq.add(&a, &b), c, "Φ(P+Q) ≠ Φ(P) + Φ(Q)");
            assert!(!jacq.is_identity(&a));
            done += 1;
        }
    }

    fn small_ctx(p: u64, seed: u64) -> (Spec, Ctx, Vec<BaseEl>, HashMap<u64, usize>, StdRng) {
        let spec = generate_spec(p, seed);
        let ctx = Ctx::new(&spec);
        let mut rng = StdRng::seed_from_u64(seed ^ 77);
        let base = factor_base(&ctx, &mut rng);
        let by_x: HashMap<u64, usize> = base.iter().enumerate().map(|(i, b)| (b.x, i)).collect();
        (spec, ctx, base, by_x, rng)
    }

    #[test]
    fn instance_transfers_the_dlp() {
        let (spec, ctx, base, _, _) = small_ctx(53, 1);
        let jac = ctx.jac();
        assert!(is_prime_u64(spec.l));
        assert!(jac.is_identity(&jac.mul(&spec.gj, spec.l as u128)));
        assert_eq!(jac.mul(&spec.gj, spec.d as u128), spec.qj);
        assert!(base.len() > 10, "{}", base.len());
        for b in &base {
            assert!(jac.is_valid(&b.d));
        }
    }

    #[test]
    fn nagao_finds_planted_six_sums_and_agrees_with_the_oracle() {
        let (_, ctx, base, by_x, mut rng) = small_ctx(53, 2);
        let jac = ctx.jac();
        let t = table3(&jac, &base);
        let opts = jv_options(24, 120.0);
        let mut found_planted = 0;
        for trial in 0..6 {
            let mut idx: Vec<usize> = Vec::new();
            while idx.len() < 6 {
                let i = rng.gen_range(0..base.len());
                if !idx.contains(&i) {
                    idx.push(i);
                }
            }
            let mut r = jac.identity();
            let mut terms: Vec<(usize, i8)> = Vec::new();
            for &i in &idx {
                let e: i8 = if rng.gen_bool(0.5) { 1 } else { -1 };
                let b = if e == 1 {
                    base[i].d.clone()
                } else {
                    jac.neg(&base[i].d)
                };
                r = jac.add(&r, &b);
                terms.push((i, e));
            }
            terms.sort_unstable();
            let want = Dec6 {
                terms: core::array::from_fn(|k| terms[k]),
            };
            let (found, cost) = nagao_decompose(&ctx, &base, &by_x, &r, &opts, &mut rng);
            let oracle = mitm6(&jac, &t, &r);
            assert!(oracle.contains(&want));
            assert_eq!(found, oracle, "trial {trial}: {cost:?}");
            assert!(found.iter().all(|d| verify_dec(&jac, &base, &r, d)));
            if found.contains(&want) {
                found_planted += 1;
            }
            eprintln!("trial {trial}: {cost:?}");
        }
        assert_eq!(found_planted, 6);
    }

    #[test]
    fn charpoly_satisfies_cayley_hamilton_and_has_the_right_roots() {
        let mut rng = StdRng::seed_from_u64(31);
        for &p in &[53u64, 1009] {
            for n in [1usize, 2, 5, 12, 20] {
                let m: Vec<Vec<u64>> = (0..n)
                    .map(|_| (0..n).map(|_| rng.gen_range(0..p)).collect())
                    .collect();
                let mut cnt = 0;
                let cp = charpoly_mod_p(&m, p, &mut cnt);
                assert_eq!(cp.len(), n + 1);
                assert_eq!(cp[n], 1);
                // cp(M) = 0: Horner on matrices
                let mut acc = vec![vec![0u64; n]; n];
                for c in cp.iter().rev() {
                    let mut next = vec![vec![0u64; n]; n];
                    for i in 0..n {
                        for j in 0..n {
                            let mut v = 0u64;
                            for k in 0..n {
                                v = am(v, mm(acc[i][k], m[k][j], p), p);
                            }
                            next[i][j] = v;
                        }
                        next[i][i] = am(next[i][i], *c, p);
                    }
                    acc = next;
                }
                assert!(acc.iter().all(|r| r.iter().all(|&x| x == 0)), "p={p} n={n}");
            }
        }
        // a diagonal matrix: the eigenvalues are the roots
        let p = 53;
        let d = [3u64, 3, 17, 40];
        let m: Vec<Vec<u64>> = (0..4)
            .map(|i| (0..4).map(|j| if i == j { d[i] } else { 0 }).collect())
            .collect();
        let mut cnt = 0;
        let cp = charpoly_mod_p(&m, p, &mut cnt);
        for &r in &d {
            let v = cp.iter().rev().fold(0u64, |a, &c| am(mm(a, r, p), c, p));
            assert_eq!(v, 0);
        }
    }

    #[test]
    fn the_solver_loses_no_rational_solution_at_small_p() {
        // p = 53: a Krylov minimal polynomial of one vector loses a root with
        // probability 1/p, a non-separating form likewise.  Plant many sums.
        let (_, ctx, base, by_x, mut rng) = small_ctx(53, 6);
        let jac = ctx.jac();
        let opts = jv_options(24, 120.0);
        for trial in 0..120 {
            let (r, want) = planted(&jac, &base, &mut rng);
            let (found, cost) = nagao_decompose(&ctx, &base, &by_x, &r, &opts, &mut rng);
            assert!(
                !cost.incomplete && !cost.timed_out,
                "trial {trial}: {cost:?}"
            );
            assert!(
                found.contains(&want),
                "trial {trial}: missed {want:?}, {cost:?}"
            );
        }
    }

    /// Slow (40 s): 1,500 planted six-sums at `p = 61`.  The only admissible
    /// misses are residuals the solver flags `incomplete` (a positive-dimensional
    /// ideal, about 1 in 700), and there must be no other.  Two defects were
    /// found by this test: a Krylov vector missing a root's eigenspace, and
    /// the reduced divisor containing a point over a planted abscissa, where
    /// `µ(x) = u(x) = A(x) = 0` leaves `y` undetermined (both signs are now
    /// tried and the group decides).
    #[test]
    #[ignore]
    fn solver_is_exact_on_1500_planted_sums_at_p61() {
        let (_, ctx, base, by_x, mut rng) = small_ctx(61, 1);
        let jac = ctx.jac();
        let opts = jv_options(24, 120.0);
        let (mut flagged, mut missed) = (0, 0);
        for _ in 0..1500 {
            let (r, want) = planted(&jac, &base, &mut rng);
            let (found, cost) = nagao_decompose(&ctx, &base, &by_x, &r, &opts, &mut rng);
            if found.contains(&want) {
                continue;
            }
            if cost.incomplete {
                flagged += 1;
            } else {
                missed += 1;
            }
        }
        assert_eq!(missed, 0, "unflagged misses (flagged: {flagged})");
    }

    #[test]
    fn stopped_f4_agrees_with_the_full_solver_and_the_oracle() {
        let (spec, ctx, base, by_x, mut rng) = small_ctx(53, 7);
        let jac = ctx.jac();
        let t = table3(&jac, &base);
        let full = jv_options(24, 120.0);
        let stop = jv_options_stopped(24, 120.0, 64);
        let (mut stopped, mut cheaper, mut n) = (0, 0, 0);
        for k in 0..30 {
            let r = if k % 2 == 0 {
                planted(&jac, &base, &mut rng).0
            } else {
                let a = rng.gen_range(0..spec.l);
                jac.add(&jac.mul(&spec.gj, a as u128), &spec.qj)
            };
            let (fa, ca) = nagao_decompose(&ctx, &base, &by_x, &r, &full, &mut rng);
            let (fb, cb) = nagao_decompose(&ctx, &base, &by_x, &r, &stop, &mut rng);
            if ca.degenerate || ca.incomplete || cb.incomplete || ca.timed_out || cb.timed_out {
                continue;
            }
            n += 1;
            let oracle = mitm6(&jac, &t, &r);
            assert_eq!(fa, oracle, "full solver vs oracle, residual {k}");
            assert_eq!(fb, oracle, "stopped solver vs oracle, residual {k}");
            assert!(fb.iter().all(|d| verify_dec(&jac, &base, &r, d)));
            if cb.stopped_at == Some(64) {
                stopped += 1;
            }
            // "cheaper" by the F4 degree and matrix the stop saves, which are
            // per-run facts; `f4_muls` is a difference of the process-wide
            // counter and interleaves with other tests' F4 runs under
            // `cargo test`'s parallelism.
            if cb.degree_reached < ca.degree_reached && cb.max_rows < ca.max_rows {
                cheaper += 1;
            }
        }
        assert!(n >= 25, "{n}");
        assert_eq!(stopped, n, "every system stops at the Bézout staircase");
        assert_eq!(
            cheaper, n,
            "every stopped run ends at a lower degree with a smaller matrix"
        );
    }

    #[test]
    fn distinguished_point_rho_finds_the_planted_logarithm() {
        // p = 53: l ≈ 2^32, sqrt l ≈ 6.5·10⁴ steps; three runs, two walkers
        let spec = generate_spec(53, 1);
        for run in 0..3u64 {
            let r = rho_e_dp(&spec, 100 + run, 6, 2);
            assert!(r.correct, "run {run}: {r:?}");
            assert!(r.s > 0.2 && r.s < 8.0, "run {run}: S = {}", r.s);
        }
    }

    #[test]
    fn replayed_f4_agrees_with_the_full_solver_and_the_oracle() {
        let (spec, ctx, base, by_x, mut rng) = small_ctx(53, 7);
        let traced = Ctx::with_trace(&spec, true);
        clear_traces();
        let jac = ctx.jac();
        let t = table3(&jac, &base);
        let stop = jv_options_stopped(24, 120.0, 64);
        let (mut replayed, mut cheaper, mut n, mut mismatches) = (0, 0, 0, 0);
        for k in 0..40 {
            let r = if k % 2 == 0 {
                planted(&jac, &base, &mut rng).0
            } else {
                let a = rng.gen_range(0..spec.l);
                jac.add(&jac.mul(&spec.gj, a as u128), &spec.qj)
            };
            let (fa, ca) = nagao_decompose(&ctx, &base, &by_x, &r, &stop, &mut rng);
            let (fb, cb) = nagao_decompose(&traced, &base, &by_x, &r, &stop, &mut rng);
            if ca.degenerate || ca.incomplete || cb.incomplete || ca.timed_out || cb.timed_out {
                continue;
            }
            n += 1;
            let oracle = mitm6(&jac, &t, &r);
            assert_eq!(fa, oracle, "stopped solver vs oracle, residual {k}");
            assert_eq!(fb, oracle, "replayed solver vs oracle, residual {k}");
            assert!(fb.iter().all(|d| verify_dec(&jac, &base, &r, d)));
            assert_eq!(cb.stopped_at, Some(64));
            if cb.replayed {
                replayed += 1;
                // the replay keeps only the rows that produced a pivot, so its
                // matrices are smaller at the same degree (per-run facts; the
                // multiplication counter interleaves with other tests' runs)
                if cb.max_rows < ca.max_rows && cb.degree_reached <= ca.degree_reached {
                    cheaper += 1;
                }
                if replayed <= 3 {
                    eprintln!(
                        "full f4 {} ({} × {}, degree {}) vs replay {} ({} × {}, degree {})",
                        ca.f4_muls,
                        ca.max_rows,
                        ca.max_cols,
                        ca.degree_reached,
                        cb.f4_muls,
                        cb.max_rows,
                        cb.max_cols,
                        cb.degree_reached
                    );
                }
            }
            if cb.trace_mismatch {
                mismatches += 1;
            }
        }
        assert!(n >= 30, "{n}");
        assert!(
            replayed >= n - 1 - mismatches,
            "replayed {replayed} of {n} ({mismatches} mismatches)"
        );
        eprintln!("replayed {replayed} of {n}, cheaper {cheaper}, mismatches {mismatches}");
        let t = trace_cache()
            .lock()
            .unwrap()
            .get(&53)
            .map(|e| e.trace.clone())
            .unwrap();
        for (i, st) in t.steps.iter().enumerate() {
            eprintln!(
                "  trace step {i}: degree {} useful rows {} new lms {}",
                st.degree,
                st.rows.len(),
                st.new_lms.len()
            );
        }
        assert!(
            cheaper * 10 >= replayed * 9,
            "cheaper {cheaper} of {replayed} replays"
        );
    }
}
