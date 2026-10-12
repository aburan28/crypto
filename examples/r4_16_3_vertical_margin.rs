//! Round 4, lead 16.3: the vertical security margin.
//!
//! How much does it cost to *find* and *evaluate* a vertical `ℓ`-isogeny
//! (a transfer across a conductor gap containing the prime `ℓ`), counted in
//! field multiplications over a 61-bit toy prime field, as `ℓ` grows?
//!
//! Instances.  Every instance is a crater curve of an `ℓ`-volcano built by
//! the CM method with `D = −7` (`j = −3375`, class number one): a prime
//! `p = (t² + 7v²)/4` in `[2^60, 2^61)` with `ℓ | v`, so `ℓ²` divides the
//! Frobenius discriminant, and `t ≡ 2λ (mod ℓ)`, so Frobenius acts on `E[ℓ]`
//! as the scalar `λ`.  Every one of the `ℓ + 1` lines in `E[ℓ]` is then a
//! rational kernel; the `ℓ − (−7/ℓ)` non-horizontal ones are descending
//! (vertical) isogenies to curves with conductor `ℓ`.  Because the crater has
//! a single vertex, a codomain with `j ≠ −3375` certifies a vertical edge.
//!
//! * `λ = 1`: `E[ℓ] ⊂ E(F_p)`, the rational-kernel case.  Kernel finding is a
//!   cofactor ladder; evaluation is timed for Vélu (`O(ℓ)`) and √élu
//!   (Bernstein–De Feo–Leroux–Smith 2020) separately, for `ℓ` up to `2^20`.
//! * `λ` a primitive root mod `ℓ`: the kernel's x-coordinates generate
//!   `F_{p^k}`, `k = (ℓ − 1)/2`, the generic case.  Kernel finding is a
//!   cofactor ladder in `F_{p^k}` plus a minimal polynomial (the "torsion
//!   point over an extension" route); evaluation then uses the kernel
//!   polynomial over `F_p` (`O(ℓ)` per point).  The library's own Elkies
//!   route (`Φ_ℓ mod p` from the `q`-series, a root, Elkies' kernel) is timed
//!   alongside it in contended wall time, because it is not instrumented.
//!
//! Unit.  `M` counts every multiplication or squaring in `F_p` (including
//! the `a_i·b_j` products inside polynomial and extension-field arithmetic);
//! `I` counts `F_p` inversions (extended Euclid), reported separately and
//! converted at `I = 91 M` (the Fermat cost at 61 bits) where a single
//! number is needed.  Additions are not counted.
//!
//! Validation.  Vélu and √élu must agree on the codomain and on every image;
//! images are checked to have order dividing `#E` on the codomain; kernel
//! polynomials are certified by the library's `verify_kernel` (division
//! polynomial test, independent of how they were found) and the library's
//! Vélu codomain must have the same `j`; for small `ℓ` the codomain `j` must
//! be a root of the library's `Φ_ℓ(−3375, Y)`.
//!
//! Usage:
//!   r4_16_3_vertical_margin rational <max_log2_ell> <seeds> <out.json>
//!   r4_16_3_vertical_margin generic <max_ell> <time_cap_s> <out.json>
//!   r4_16_3_vertical_margin elkies <max_ell> <time_cap_s> <out.json>

use crypto_lib::cryptanalysis::isogeny_walk::curve::Model;
use crypto_lib::cryptanalysis::isogeny_walk::field::Field as LField;
use crypto_lib::cryptanalysis::isogeny_walk::{kernel, modpoly, poly as lpoly};
use num_bigint::{BigInt, BigUint};
use num_traits::{Signed, ToPrimitive, Zero};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::cell::Cell;
use std::time::Instant;

// ---------------------------------------------------------------------------
// Operation counters.

thread_local! {
    static MULS: Cell<u64> = const { Cell::new(0) };
    static INVS: Cell<u64> = const { Cell::new(0) };
}

fn add_m(n: u64) {
    MULS.with(|c| c.set(c.get() + n));
}

fn add_i() {
    INVS.with(|c| c.set(c.get() + 1));
}

#[derive(Clone, Copy, Debug, Default)]
struct Cnt {
    m: u64,
    i: u64,
}

impl Cnt {
    fn now() -> Cnt {
        Cnt {
            m: MULS.with(|c| c.get()),
            i: INVS.with(|c| c.get()),
        }
    }
    fn since(start: Cnt) -> Cnt {
        let n = Cnt::now();
        Cnt {
            m: n.m - start.m,
            i: n.i - start.i,
        }
    }
    fn plus(self, o: Cnt) -> Cnt {
        Cnt {
            m: self.m + o.m,
            i: self.i + o.i,
        }
    }
    /// Multiplication-equivalents with `I = 91 M`.
    fn meq(&self) -> u64 {
        self.m + 91 * self.i
    }
    fn js(&self) -> Value {
        json!({"M": self.m, "I": self.i, "Meq": self.meq()})
    }
}

// ---------------------------------------------------------------------------
// F_p, p < 2^62.

#[derive(Clone, Copy, Debug)]
struct Fp {
    p: u64,
}

impl Fp {
    #[inline]
    fn mul(&self, a: u64, b: u64) -> u64 {
        add_m(1);
        ((a as u128 * b as u128) % self.p as u128) as u64
    }
    #[inline]
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    #[inline]
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }
    #[inline]
    fn neg(&self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }
    fn red(&self, v: i128) -> u64 {
        v.rem_euclid(self.p as i128) as u64
    }
    fn inv(&self, a: u64) -> u64 {
        add_i();
        assert!(a != 0, "inverse of zero");
        let (mut r0, mut r1) = (self.p as i128, a as i128);
        let (mut s0, mut s1) = (0i128, 1i128);
        while r1 != 0 {
            let q = r0 / r1;
            (r0, r1) = (r1, r0 - q * r1);
            (s0, s1) = (s1, s0 - q * s1);
        }
        self.red(s0)
    }
    fn pow(&self, a: u64, mut e: u128) -> u64 {
        let mut r = 1u64;
        let mut b = a;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, b);
            }
            b = self.mul(b, b);
            e >>= 1;
        }
        r
    }
    fn legendre(&self, a: u64) -> i32 {
        if a == 0 {
            return 0;
        }
        if self.pow(a, ((self.p - 1) / 2) as u128) == 1 {
            1
        } else {
            -1
        }
    }
    fn sqrt(&self, a: u64) -> Option<u64> {
        if a == 0 {
            return Some(0);
        }
        if self.legendre(a) != 1 {
            return None;
        }
        // Tonelli–Shanks.
        let mut q = self.p - 1;
        let mut s = 0;
        while q.is_multiple_of(2) {
            q /= 2;
            s += 1;
        }
        let mut z = 2u64;
        while self.legendre(z) != -1 {
            z += 1;
        }
        let mut m = s;
        let mut c = self.pow(z, q as u128);
        let mut t = self.pow(a, q as u128);
        let mut r = self.pow(a, q.div_ceil(2) as u128);
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
}

// ---------------------------------------------------------------------------
// Polynomials over F_p (little-endian coefficient vectors), Karatsuba.

const KARA: usize = 24;

fn trim(mut a: Vec<u64>) -> Vec<u64> {
    while a.last() == Some(&0) {
        a.pop();
    }
    a
}

fn school(f: &Fp, a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    add_m((a.len() * b.len()) as u64);
    let p = f.p as u128;
    let mut acc = vec![0u128; a.len() + b.len() - 1];
    // p < 2^61: each product < 2^122, so 24 terms fit; reduce every 24 rows.
    for (i, &ai) in a.iter().enumerate() {
        for (j, &bj) in b.iter().enumerate() {
            acc[i + j] += ai as u128 * bj as u128;
        }
        if i % 24 == 23 {
            for v in acc.iter_mut() {
                *v %= p;
            }
        }
    }
    acc.into_iter().map(|v| (v % p) as u64).collect()
}

fn padd(f: &Fp, a: &[u64], b: &[u64]) -> Vec<u64> {
    let n = a.len().max(b.len());
    (0..n)
        .map(|i| f.add(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0)))
        .collect()
}

fn psub(f: &Fp, a: &[u64], b: &[u64]) -> Vec<u64> {
    let n = a.len().max(b.len());
    (0..n)
        .map(|i| f.sub(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0)))
        .collect()
}

fn add_into(f: &Fp, acc: &mut [u64], src: &[u64], off: usize) {
    for (i, &v) in src.iter().enumerate() {
        acc[off + i] = f.add(acc[off + i], v);
    }
}

fn kmul(f: &Fp, a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let (a, b) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    if b.len() <= KARA {
        return school(f, a, b);
    }
    let mut out = vec![0u64; a.len() + b.len() - 1];
    if a.len() >= 2 * b.len() {
        // Unbalanced: chunks of a against b.
        let mut off = 0;
        while off < a.len() {
            let end = (off + b.len()).min(a.len());
            let part = kmul(f, &a[off..end], b);
            add_into(f, &mut out, &part, off);
            off = end;
        }
        return out;
    }
    let m = a.len() / 2;
    let (a0, a1) = a.split_at(m);
    if b.len() <= m {
        let z0 = kmul(f, a0, b);
        let z1 = kmul(f, a1, b);
        add_into(f, &mut out, &z0, 0);
        add_into(f, &mut out, &z1, m);
        return out;
    }
    let (b0, b1) = b.split_at(m);
    let z0 = kmul(f, a0, b0);
    let z2 = kmul(f, a1, b1);
    let z1 = kmul(f, &padd(f, a0, a1), &padd(f, b0, b1));
    let mid = psub(f, &psub(f, &z1, &z0), &z2);
    add_into(f, &mut out, &z0, 0);
    add_into(f, &mut out, &z2, 2 * m);
    let mid = trim(mid);
    add_into(f, &mut out, &mid, m);
    out
}

/// `1/a mod x^n` (`a[0] ≠ 0`) by Newton iteration.
fn inv_series(f: &Fp, a: &[u64], n: usize) -> Vec<u64> {
    let mut g = vec![f.inv(a[0])];
    let mut m = 1;
    while m < n {
        let m2 = (2 * m).min(n);
        let at: Vec<u64> = a.iter().take(m2).copied().collect();
        let mut e = kmul(f, &at, &g);
        e.truncate(m2);
        e.resize(m2, 0);
        // e = 2 − a·g
        for v in e.iter_mut() {
            *v = f.neg(*v);
        }
        e[0] = f.add(e[0], 2);
        let mut ng = kmul(f, &g, &e);
        ng.truncate(m2);
        g = ng;
        m = m2;
    }
    g.resize(n, 0);
    g
}

/// Remainder of `a` by monic `b`; `binv` caches `1/rev(b)` to some precision.
fn prem(f: &Fp, a: &[u64], b: &[u64], binv: &mut Vec<u64>) -> Vec<u64> {
    let a = trim(a.to_vec());
    let db = b.len() - 1;
    if a.len() <= db {
        return a;
    }
    let dq = a.len() - 1 - db; // quotient degree
    if db < 32 || dq < 32 {
        // Classical division, b monic.
        let mut r = a;
        for i in (db..r.len()).rev() {
            let q = r[i];
            if q == 0 {
                continue;
            }
            for j in 0..db {
                let t = f.mul(q, b[j]);
                r[i - db + j] = f.sub(r[i - db + j], t);
            }
            r[i] = 0;
        }
        r.truncate(db);
        return trim(r);
    }
    let need = dq + 1;
    if binv.len() < need {
        let rb: Vec<u64> = b.iter().rev().copied().collect();
        *binv = inv_series(f, &rb, need);
    }
    let ra: Vec<u64> = a.iter().rev().take(need).copied().collect();
    let mut rq = kmul(f, &ra, &binv[..need]);
    rq.truncate(need);
    rq.resize(need, 0);
    let q: Vec<u64> = rq.into_iter().rev().collect();
    let qb = kmul(f, &q, b);
    let mut r = psub(f, &a, &qb);
    r.truncate(db);
    trim(r)
}

/// Product tree whose leaves are `Z − x_i`.
struct Tree {
    levels: Vec<Vec<Vec<u64>>>,
    inv: Vec<Vec<Vec<u64>>>,
}

impl Tree {
    fn new(f: &Fp, xs: &[u64]) -> Tree {
        let leaves: Vec<Vec<u64>> = xs.iter().map(|&x| vec![f.neg(x), 1]).collect();
        let mut levels = vec![leaves];
        while levels.last().unwrap().len() > 1 {
            let prev = levels.last().unwrap();
            let mut next = Vec::with_capacity(prev.len().div_ceil(2));
            for pair in prev.chunks(2) {
                if pair.len() == 2 {
                    next.push(kmul(f, &pair[0], &pair[1]));
                } else {
                    next.push(pair[0].clone());
                }
            }
            levels.push(next);
        }
        let inv = levels.iter().map(|l| vec![Vec::new(); l.len()]).collect();
        Tree { levels, inv }
    }
    /// `g(x_i)` for every leaf, by a remainder tree.
    fn multipoint(&mut self, f: &Fp, g: &[u64]) -> Vec<u64> {
        let top = self.levels.len() - 1;
        let r = prem(f, g, &self.levels[top][0].clone(), &mut self.inv[top][0]);
        let mut cur = vec![r];
        for lv in (0..top).rev() {
            let n = self.levels[lv].len();
            let mut next = Vec::with_capacity(n);
            for (j, parent) in cur.iter().enumerate() {
                for c in [2 * j, 2 * j + 1] {
                    if c < n {
                        let node = self.levels[lv][c].clone();
                        next.push(prem(f, parent, &node, &mut self.inv[lv][c]));
                    }
                }
            }
            cur = next;
        }
        cur.into_iter().map(|r| *r.first().unwrap_or(&0)).collect()
    }
}

// ---------------------------------------------------------------------------
// A field abstraction so the Montgomery ladder runs over F_p and F_{p^k}.

trait Fld {
    type E: Clone + PartialEq + std::fmt::Debug;
    fn mul(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn add(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn sub(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn zero(&self) -> Self::E;
    fn is_zero(&self, a: &Self::E) -> bool;
}

impl Fld for Fp {
    type E = u64;
    fn mul(&self, a: &u64, b: &u64) -> u64 {
        Fp::mul(self, *a, *b)
    }
    fn add(&self, a: &u64, b: &u64) -> u64 {
        Fp::add(self, *a, *b)
    }
    fn sub(&self, a: &u64, b: &u64) -> u64 {
        Fp::sub(self, *a, *b)
    }
    fn zero(&self) -> u64 {
        0
    }
    fn is_zero(&self, a: &u64) -> bool {
        *a == 0
    }
}

/// `F_{p^k} = F_p[X]/(X^k + X^m + c)`.
#[derive(Clone, Debug)]
struct Fpk {
    f: Fp,
    k: usize,
    m: usize,
    c: u64,
}

impl Fpk {
    fn reduce(&self, mut r: Vec<u64>) -> Vec<u64> {
        let k = self.k;
        r.resize((2 * k).max(k), 0);
        for i in (k..r.len()).rev() {
            let q = r[i];
            if q == 0 {
                continue;
            }
            r[i] = 0;
            let j = i - k + self.m;
            r[j] = self.f.sub(r[j], q);
            let t = self.f.mul(q, self.c);
            r[i - k] = self.f.sub(r[i - k], t);
        }
        r.truncate(k);
        r
    }
    fn modulus(&self) -> Vec<u64> {
        let mut v = vec![0u64; self.k + 1];
        v[0] = self.c;
        v[self.m] = self.f.add(v[self.m], 1);
        v[self.k] = 1;
        v
    }
    fn x(&self) -> Vec<u64> {
        let mut v = vec![0u64; self.k];
        if self.k > 1 {
            v[1] = 1;
        } else {
            v[0] = self.f.neg(self.c);
        }
        v
    }
    fn one(&self) -> Vec<u64> {
        let mut v = vec![0u64; self.k];
        v[0] = 1;
        v
    }
    fn pow_big(&self, a: &[u64], e: &BigUint) -> Vec<u64> {
        let mut r = self.one();
        for i in (0..e.bits()).rev() {
            r = Fld::mul(self, &r, &r);
            if e.bit(i) {
                r = Fld::mul(self, &r, &a.to_vec());
            }
        }
        r
    }
    fn inv(&self, a: &[u64]) -> Vec<u64> {
        let g = poly_xgcd_inv(&self.f, &trim(a.to_vec()), &self.modulus())
            .expect("non-invertible extension element");
        let mut g = g;
        g.resize(self.k, 0);
        g
    }
}

impl Fld for Fpk {
    type E = Vec<u64>;
    fn mul(&self, a: &Vec<u64>, b: &Vec<u64>) -> Vec<u64> {
        self.reduce(kmul(&self.f, a, b))
    }
    fn add(&self, a: &Vec<u64>, b: &Vec<u64>) -> Vec<u64> {
        a.iter().zip(b).map(|(x, y)| self.f.add(*x, *y)).collect()
    }
    fn sub(&self, a: &Vec<u64>, b: &Vec<u64>) -> Vec<u64> {
        a.iter().zip(b).map(|(x, y)| self.f.sub(*x, *y)).collect()
    }
    fn zero(&self) -> Vec<u64> {
        vec![0u64; self.k]
    }
    fn is_zero(&self, a: &Vec<u64>) -> bool {
        a.iter().all(|v| *v == 0)
    }
}

/// Classical polynomial remainder (any nonzero `b`).
fn poly_rem(f: &Fp, a: &[u64], b: &[u64]) -> Vec<u64> {
    let b = trim(b.to_vec());
    let db = b.len() - 1;
    let lc_inv = f.inv(b[db]);
    let mut r = trim(a.to_vec());
    while r.len() > db {
        let i = r.len() - 1;
        let q = f.mul(r[i], lc_inv);
        for j in 0..=db {
            let t = f.mul(q, b[j]);
            r[i - db + j] = f.sub(r[i - db + j], t);
        }
        r = trim(r);
    }
    r
}

fn poly_gcd(f: &Fp, a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut a = trim(a.to_vec());
    let mut b = trim(b.to_vec());
    while !b.is_empty() {
        let r = poly_rem(f, &a, &b);
        a = b;
        b = r;
    }
    a
}

/// `a^{-1} mod m` by the extended Euclidean algorithm.
fn poly_xgcd_inv(f: &Fp, a: &[u64], m: &[u64]) -> Option<Vec<u64>> {
    let (mut r0, mut r1) = (trim(m.to_vec()), trim(a.to_vec()));
    let (mut s0, mut s1): (Vec<u64>, Vec<u64>) = (Vec::new(), vec![1]);
    while !r1.is_empty() {
        // q, r = divmod(r0, r1)
        let d1 = r1.len() - 1;
        let lc_inv = f.inv(r1[d1]);
        let mut r = r0.clone();
        let mut q = vec![0u64; r0.len().saturating_sub(d1).max(1)];
        while r.len() > d1 {
            let i = r.len() - 1;
            let c = f.mul(r[i], lc_inv);
            q[i - d1] = c;
            for j in 0..=d1 {
                let t = f.mul(c, r1[j]);
                r[i - d1 + j] = f.sub(r[i - d1 + j], t);
            }
            r = trim(r);
        }
        let ns = psub(f, &s0, &trim(kmul(f, &trim(q), &s1)));
        (r0, r1) = (r1, r);
        (s0, s1) = (s1, trim(ns));
    }
    if r0.len() != 1 {
        return None;
    }
    let c = f.inv(r0[0]);
    Some(s0.iter().map(|v| f.mul(*v, c)).collect())
}

// ---------------------------------------------------------------------------
// x-only Montgomery arithmetic, `By² = x³ + Ax² + x`, points as `(X : Z)`.

#[derive(Clone, Debug, PartialEq)]
struct Xz<E> {
    x: E,
    z: E,
}

fn xdbl<F: Fld>(f: &F, p: &Xz<F::E>, a24: &F::E) -> Xz<F::E> {
    let t0 = f.add(&p.x, &p.z);
    let t1 = f.sub(&p.x, &p.z);
    let t0 = f.mul(&t0, &t0);
    let t1 = f.mul(&t1, &t1);
    let x = f.mul(&t0, &t1);
    let t = f.sub(&t0, &t1);
    let z = f.mul(&t, &f.add(&t1, &f.mul(a24, &t)));
    Xz { x, z }
}

fn xadd<F: Fld>(f: &F, p: &Xz<F::E>, q: &Xz<F::E>, d: &Xz<F::E>) -> Xz<F::E> {
    let u = f.mul(&f.sub(&p.x, &p.z), &f.add(&q.x, &q.z));
    let v = f.mul(&f.add(&p.x, &p.z), &f.sub(&q.x, &q.z));
    let s = f.add(&u, &v);
    let t = f.sub(&u, &v);
    let x = f.mul(&d.z, &f.mul(&s, &s));
    let z = f.mul(&d.x, &f.mul(&t, &t));
    Xz { x, z }
}

/// `[n]P` by the Montgomery ladder (`n ≥ 1`), base point projective.
fn ladder<F: Fld>(f: &F, p: &Xz<F::E>, n: &BigUint, a24: &F::E) -> Xz<F::E> {
    if n.is_zero() {
        return Xz {
            x: p.x.clone(),
            z: f.zero(),
        };
    }
    let mut r0 = p.clone();
    let mut r1 = xdbl(f, p, a24);
    for i in (0..n.bits() - 1).rev() {
        if n.bit(i) {
            r0 = xadd(f, &r1, &r0, p);
            r1 = xdbl(f, &r1, a24);
        } else {
            r1 = xadd(f, &r1, &r0, p);
            r0 = xdbl(f, &r0, a24);
        }
    }
    r0
}

fn affine(f: &Fp, p: &Xz<u64>) -> u64 {
    f.mul(p.x, f.inv(p.z))
}

/// `j(A) = 256(A² − 3)³/(A² − 4)`.
fn j_of_a(f: &Fp, a: u64) -> u64 {
    let a2 = f.mul(a, a);
    let t = f.sub(a2, 3);
    let num = f.mul(256, f.mul(t, f.mul(t, t)));
    f.mul(num, f.inv(f.sub(a2, 4)))
}

// ---------------------------------------------------------------------------
// Number theory helpers.

fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for sp in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n.is_multiple_of(sp) {
            return n == sp;
        }
    }
    let mut d = n - 1;
    let mut s = 0;
    while d.is_multiple_of(2) {
        d /= 2;
        s += 1;
    }
    let mulm = |a: u64, b: u64| ((a as u128 * b as u128) % n as u128) as u64;
    let powm = |mut a: u64, mut e: u64| {
        let mut r = 1u64;
        while e > 0 {
            if e & 1 == 1 {
                r = mulm(r, a);
            }
            a = mulm(a, a);
            e >>= 1;
        }
        r
    };
    'w: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powm(a, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulm(x, x);
            if x == n - 1 {
                continue 'w;
            }
        }
        return false;
    }
    true
}

fn next_prime(mut n: u64) -> u64 {
    while !is_prime_u64(n) {
        n += 1;
    }
    n
}

fn primitive_root(l: u64) -> u64 {
    let mut fac = Vec::new();
    let mut m = l - 1;
    let mut d = 2;
    while d * d <= m {
        if m.is_multiple_of(d) {
            fac.push(d);
            while m.is_multiple_of(d) {
                m /= d;
            }
        }
        d += 1;
    }
    if m > 1 {
        fac.push(m);
    }
    let pw = |a: u64, mut e: u64| {
        let mut r = 1u128;
        let mut b = a as u128;
        while e > 0 {
            if e & 1 == 1 {
                r = r * b % l as u128;
            }
            b = b * b % l as u128;
            e >>= 1;
        }
        r as u64
    };
    (2..l)
        .find(|&g| fac.iter().all(|&q| pw(g, (l - 1) / q) != 1))
        .unwrap()
}

fn val(n: &BigUint, l: u64) -> (u32, BigUint) {
    let lb = BigUint::from(l);
    let mut n = n.clone();
    let mut e = 0;
    while !n.is_zero() && (&n % &lb).is_zero() {
        n /= &lb;
        e += 1;
    }
    (e, n)
}

// ---------------------------------------------------------------------------
// Instances: crater curves of an ℓ-volcano with CM by Z[(1+√−7)/2].

#[derive(Clone, Debug)]
struct Inst {
    ell: u64,
    lambda: u64,
    p: u64,
    t: i64,
    v: u64,
    n: u128,
    /// Montgomery coefficient and `(A + 2)/4`.
    a: u64,
    a24: u64,
    /// `χ(x³ + Ax² + x)` of the x-coordinates of the order-`N` curve.
    chi_n: i32,
    /// Short Weierstrass `y² = x³ − 35x − 98`; `x_w = s·x_m + α`.
    alpha: u64,
    s: u64,
    jcrater: u64,
}

fn make_instance(ell: u64, lambda: u64, rng: &mut StdRng) -> Inst {
    let lo: u128 = 1 << 62;
    let hi: u128 = 1 << 63;
    loop {
        let wmax = (((hi - 1) / 7) as f64).sqrt() as u128 / ell as u128;
        let w = rng.gen_range(1..=wmax.max(1)) as u64;
        let v = ell * w;
        let sv = 7u128 * (v as u128) * (v as u128);
        if sv >= hi {
            continue;
        }
        let tmax = ((hi - 1 - sv) as f64).sqrt() as i128;
        let tmin = if sv >= lo {
            0
        } else {
            ((lo - sv) as f64).sqrt().ceil() as i128
        };
        if tmin > tmax {
            continue;
        }
        // t ≡ 2λ (mod ℓ), t ≡ v (mod 2): one class mod 2ℓ.
        let l = ell as i128;
        let r_l = (2 * lambda as i128) % l;
        let r2 = (v % 2) as i128;
        let r = if r_l % 2 == r2 { r_l } else { r_l + l };
        let span = (tmax - tmin) / (2 * l);
        if span < 1 {
            continue;
        }
        let j = rng.gen_range(0..=span);
        let mut t = tmin - tmin.rem_euclid(2 * l) + r + 2 * l * j;
        if t < tmin {
            t += 2 * l;
        }
        if t > tmax {
            continue;
        }
        let four_p = (t * t) as u128 + sv;
        if !four_p.is_multiple_of(4) || four_p < lo || four_p >= hi {
            continue;
        }
        let p = (four_p / 4) as u64;
        if !is_prime_u64(p) {
            continue;
        }
        let n = (p as i128 + 1 - t) as u128;
        if let Some(inst) = finish_instance(ell, lambda, p, t as i64, v, n, rng) {
            return inst;
        }
    }
}

fn finish_instance(
    ell: u64,
    lambda: u64,
    p: u64,
    t: i64,
    v: u64,
    n: u128,
    rng: &mut StdRng,
) -> Option<Inst> {
    let f = Fp { p };
    let aw = f.red(-35);
    // Roots of x³ − 35x − 98 = (x − 7)(x² + 7x + 14).
    let mut roots = vec![7u64];
    let disc = f.red(-7);
    if let Some(sq) = f.sqrt(disc) {
        let inv2 = f.inv(2);
        roots.push(f.mul(f.sub(sq, 7), inv2));
        roots.push(f.mul(f.sub(f.neg(sq), 7), inv2));
    }
    for alpha in roots {
        let w = f.add(f.mul(3, f.mul(alpha, alpha)), aw);
        let Some(s) = f.sqrt(w) else { continue };
        if s == 0 {
            continue;
        }
        let a = f.mul(f.mul(3, alpha), f.inv(s));
        let a24 = f.mul(f.add(a, 2), f.inv(4));
        let jcrater = j_of_a(&f, a);
        assert_eq!(jcrater, f.red(-3375), "CM curve must have j = −3375");
        // which side of the x-line has order N?
        let nb = BigUint::from(n);
        let mut chi_n = 0;
        for _ in 0..8 {
            let x = rng.gen_range(2..p);
            let rhs = f.mul(x, f.add(f.mul(x, x), f.add(f.mul(a, x), 1)));
            let ch = f.legendre(rhs);
            if ch == 0 {
                continue;
            }
            let q = ladder(&f, &Xz { x, z: 1 }, &nb, &a24);
            if q.z == 0 {
                if chi_n != 0 && chi_n != ch {
                    return None;
                }
                chi_n = ch;
            } else if chi_n == ch {
                return None;
            }
        }
        if chi_n == 0 {
            return None;
        }
        return Some(Inst {
            ell,
            lambda,
            p,
            t,
            v,
            n,
            a,
            a24,
            chi_n,
            alpha,
            s,
            jcrater,
        });
    }
    None
}

// ---------------------------------------------------------------------------
// Kernel finding, rational case: a point of order ℓ in E(F_p).

fn find_rational_kernel(inst: &Inst, rng: &mut StdRng) -> (u64, u32) {
    let f = Fp { p: inst.p };
    let (e, m) = val(&BigUint::from(inst.n), inst.ell);
    assert!(e >= 2, "vertical λ=1 instance must have ℓ² | N");
    let lb = BigUint::from(inst.ell);
    let mut tries = 0;
    loop {
        tries += 1;
        let x = rng.gen_range(2..inst.p);
        let rhs = f.mul(x, f.add(f.mul(x, x), f.add(f.mul(inst.a, x), 1)));
        if f.legendre(rhs) != inst.chi_n {
            continue;
        }
        let mut q = ladder(&f, &Xz { x, z: 1 }, &m, &inst.a24);
        if q.z == 0 {
            continue;
        }
        for _ in 0..e {
            let r = ladder(&f, &q, &lb, &inst.a24);
            if r.z == 0 {
                return (affine(&f, &q), tries);
            }
            q = r;
        }
        panic!("ℓ-part larger than v_ℓ(N)");
    }
}

// ---------------------------------------------------------------------------
// Vélu, O(ℓ), x-only Montgomery (Costello–Hisil images, Meyer–Reith codomain).

struct IsoOut {
    a_new: u64,
    images: Vec<u64>,
}

fn edwards_codomain(f: &Fp, a: u64, ell: u64, prod_plus: u64, prod_minus: u64) -> u64 {
    // a_E = A + 2, d_E = A − 2; a' = a^ℓ Π(x+1)^8, d' = d^ℓ Π(x−1)^8.
    let pow8 = |v: u64| {
        let v2 = f.mul(v, v);
        let v4 = f.mul(v2, v2);
        f.mul(v4, v4)
    };
    let ae = f.mul(f.pow(f.add(a, 2), ell as u128), pow8(prod_plus));
    let de = f.mul(f.pow(f.sub(a, 2), ell as u128), pow8(prod_minus));
    // A' = 2(a' + d')/(a' − d')
    f.mul(f.mul(2, f.add(ae, de)), f.inv(f.sub(ae, de)))
}

fn velu(f: &Fp, a: u64, a24: u64, ell: u64, xp: u64, pts: &[u64]) -> IsoOut {
    let d = ((ell - 1) / 2) as usize;
    let p1 = Xz { x: xp, z: 1 };
    let mut mult = vec![p1.clone()];
    if d >= 2 {
        mult.push(xdbl(f, &p1, &a24));
    }
    for i in 2..d {
        let nx = xadd(f, &mult[i - 1], &p1, &mult[i - 2]);
        mult.push(nx);
    }
    let mut pp = 1u64; // Π (X_i + Z_i)
    let mut pm = 1u64; // Π (X_i − Z_i)
    let sums: Vec<(u64, u64)> = mult
        .iter()
        .map(|q| (f.sub(q.x, q.z), f.add(q.x, q.z)))
        .collect();
    for (m_, p_) in &sums {
        pm = f.mul(pm, *m_);
        pp = f.mul(pp, *p_);
    }
    let a_new = edwards_codomain(f, a, ell, pp, pm);
    let mut images = Vec::new();
    for &x in pts {
        let (xm, xp_) = (f.sub(x, 1), f.add(x, 1));
        let mut num = 1u64;
        let mut den = 1u64;
        for (m_, p_) in &sums {
            let t0 = f.mul(*m_, xp_);
            let t1 = f.mul(*p_, xm);
            num = f.mul(num, f.add(t0, t1));
            den = f.mul(den, f.sub(t0, t1));
        }
        let xn = f.mul(x, f.mul(num, num));
        let zn = f.mul(den, den);
        images.push(f.mul(xn, f.inv(zn)));
    }
    IsoOut { a_new, images }
}

// ---------------------------------------------------------------------------
// √élu (Bernstein–De Feo–Leroux–Smith 2020), x-only Montgomery.

fn odd_multiples(f: &Fp, a24: u64, p1: &Xz<u64>, start: u64, count: usize) -> Vec<Xz<u64>> {
    // x([start + 2i]P), i < count, start odd.
    if count == 0 {
        return Vec::new();
    }
    let p2 = xdbl(f, p1, &a24);
    let first = ladder(f, p1, &BigUint::from(start), &a24);
    let mut out = vec![first];
    if count > 1 {
        let second = ladder(f, p1, &BigUint::from(start + 2), &a24);
        out.push(second);
    }
    for i in 2..count {
        let nx = xadd(f, &out[i - 1], &p2, &out[i - 2]);
        out.push(nx);
    }
    out
}

fn batch_affine(f: &Fp, pts: &[Xz<u64>]) -> Vec<u64> {
    let n = pts.len();
    if n == 0 {
        return Vec::new();
    }
    let mut pre = vec![1u64; n];
    let mut acc = 1u64;
    for (i, q) in pts.iter().enumerate() {
        pre[i] = acc;
        acc = f.mul(acc, q.z);
    }
    let mut inv = f.inv(acc);
    let mut out = vec![0u64; n];
    for i in (0..n).rev() {
        out[i] = f.mul(pts[i].x, f.mul(inv, pre[i]));
        inv = f.mul(inv, pts[i].z);
    }
    out
}

fn sqrt_velu(f: &Fp, a: u64, a24: u64, ell: u64, xp: u64, pts: &[u64]) -> (IsoOut, Value) {
    let b = (((ell - 1) as f64).sqrt() / 2.0).floor() as u64;
    assert!(b >= 1, "√élu needs ℓ ≥ 5");
    let bp = (ell - 1) / (4 * b);
    let p1 = Xz { x: xp, z: 1 };
    // I = {2b(2i+1)}, via Q_{i+1} = Q_i + [4b]P with difference Q_{i−1}.
    let q0 = ladder(f, &p1, &BigUint::from(2 * b), &a24);
    let q4 = xdbl(f, &q0, &a24);
    let mut ipts = vec![q0.clone()];
    if bp > 1 {
        ipts.push(xadd(f, &q0, &q4, &q0));
    }
    for i in 2..bp as usize {
        let nx = xadd(f, &ipts[i - 1], &q4, &ipts[i - 2]);
        ipts.push(nx);
    }
    let jpts = odd_multiples(f, a24, &p1, 1, b as usize);
    let kstart = 4 * b * bp + 1;
    let kcount = if kstart <= ell - 2 {
        ((ell - 2 - kstart) / 2 + 1) as usize
    } else {
        0
    };
    let kpts = odd_multiples(f, a24, &p1, kstart, kcount);
    let xi = batch_affine(f, &ipts);
    let xj = batch_affine(f, &jpts);
    let xk = batch_affine(f, &kpts);
    // J-side constants: x, x², x² + 1 + 2Ax.
    let jc: Vec<(u64, u64)> = xj
        .iter()
        .map(|&x| {
            let x2 = f.mul(x, x);
            let mid = f.add(f.add(x2, 1), f.mul(f.mul(2, a), x));
            (x, mid)
        })
        .collect();
    let mut tree = Tree::new(f, &xi);
    // Evaluate Δ·Π_{s∈S}(α − x_s) (den) and Δ·Π_{s∈S}(1 − α x_s) (num).
    // `need_num = false` skips the reversed product (not needed at α = ±1).
    let eval = |tree: &mut Tree, alpha: u64, need_num: bool| -> (u64, u64) {
        let a2 = f.mul(alpha, alpha);
        let mut quads_den = Vec::with_capacity(jc.len());
        let mut quads_num = Vec::with_capacity(jc.len());
        for &(x, mid) in &jc {
            // c2 = (α − x)², c1 = −2(α²x + α·mid + x), c0 = (αx − 1)².
            let amx = f.sub(alpha, x);
            let c2 = f.mul(amx, amx);
            let ax = f.mul(alpha, x);
            let axm1 = f.sub(ax, 1);
            let c0 = f.mul(axm1, axm1);
            let c1 = f.neg(f.mul(2, f.add(f.add(f.mul(a2, x), f.mul(alpha, mid)), x)));
            quads_den.push(vec![c0, c1, c2]);
            quads_num.push(vec![c2, c1, c0]);
        }
        let prod = |qs: Vec<Vec<u64>>| -> Vec<u64> {
            let mut cur = qs;
            while cur.len() > 1 {
                let mut nxt = Vec::with_capacity(cur.len().div_ceil(2));
                for pair in cur.chunks(2) {
                    if pair.len() == 2 {
                        nxt.push(kmul(f, &pair[0], &pair[1]));
                    } else {
                        nxt.push(pair[0].clone());
                    }
                }
                cur = nxt;
            }
            cur.pop().unwrap()
        };
        let ej_den = prod(quads_den);
        let rd = tree.multipoint(f, &ej_den);
        let mut den = rd.iter().fold(1u64, |acc, v| f.mul(acc, *v));
        for &x in &xk {
            den = f.mul(den, f.sub(alpha, x));
        }
        let mut num = 0u64;
        if need_num {
            let ej_num = prod(quads_num);
            let rn = tree.multipoint(f, &ej_num);
            num = rn.iter().fold(1u64, |acc, v| f.mul(acc, *v));
            for &x in &xk {
                num = f.mul(num, f.sub(1, f.mul(alpha, x)));
            }
        }
        (den, num)
    };
    // Codomain: Π(x_s + 1) and Π(x_s − 1) up to the common Δ and signs that
    // the eighth powers remove: den(−1) = Δ Π(−1 − x_s), den(1) = Δ Π(1 − x_s).
    let (dm1, _) = eval(&mut tree, f.neg(1), false);
    let (d1, _) = eval(&mut tree, 1, false);
    let a_new = edwards_codomain(f, a, ell, dm1, d1);
    let mut images = Vec::new();
    for &x in pts {
        let (den, num) = eval(&mut tree, x, true);
        let xn = f.mul(x, f.mul(num, num));
        let zn = f.mul(den, den);
        images.push(f.mul(xn, f.inv(zn)));
    }
    let meta = json!({"b": b, "b_prime": bp, "I": xi.len(), "J": xj.len(), "K": xk.len()});
    (IsoOut { a_new, images }, meta)
}

// ---------------------------------------------------------------------------
// Library cross-checks (uncounted).

fn lib_field(p: u64) -> LField {
    LField::new(&BigUint::from(p)).expect("prime field")
}

/// Short-Weierstrass kernel polynomial from Montgomery-x kernel polynomial:
/// `h_w(X) = Σ c_j s^{d−j} (X − α)^j`.
fn to_weierstrass_kernel(f: &Fp, hm: &[u64], alpha: u64, s: u64) -> Vec<u64> {
    let d = hm.len() - 1;
    let mut spow = vec![1u64; d + 1];
    for i in 1..=d {
        spow[i] = f.mul(spow[i - 1], s);
    }
    let mut h: Vec<u64> = Vec::new();
    for j in (0..=d).rev() {
        // h = h·(X − α) + c_j s^{d−j}
        let mut nh = vec![0u64; h.len() + 1];
        for (i, &c) in h.iter().enumerate() {
            nh[i + 1] = f.add(nh[i + 1], c);
            nh[i] = f.sub(nh[i], f.mul(c, alpha));
        }
        nh[0] = f.add(nh[0], f.mul(hm[j], spow[d - j]));
        h = nh;
    }
    trim(h)
}

/// Certify `h_w` with the library and return the library's codomain `j`.
fn lib_certify(p: u64, hw: &[u64], ell: u64) -> Result<u64, String> {
    let lf = lib_field(p);
    let e = Model {
        a: lf.from_i64(-35),
        b: lf.from_i64(-98),
    };
    let h: Vec<_> = hw.iter().map(|&c| lf.from_u64(c)).collect();
    let cod = kernel::verify_kernel(&lf, &e, ell, &h).map_err(|e| format!("{e:?}"))?;
    let j = cod.j(&lf).ok_or("singular codomain")?;
    Ok(lf.to_big(&j).to_u64().unwrap())
}

fn kernel_poly_from_points(f: &Fp, a24: u64, xp: u64, ell: u64) -> Vec<u64> {
    let d = ((ell - 1) / 2) as usize;
    let p1 = Xz { x: xp, z: 1 };
    let mut mult = vec![p1.clone()];
    if d >= 2 {
        mult.push(xdbl(f, &p1, &a24));
    }
    for i in 2..d {
        let nx = xadd(f, &mult[i - 1], &p1, &mult[i - 2]);
        mult.push(nx);
    }
    let xs = batch_affine(f, &mult);
    let mut h = vec![1u64];
    for x in xs {
        h = kmul(f, &h, &[f.neg(x), 1]);
    }
    h
}

fn order_check(f: &Fp, a_new: u64, x: u64, n: u128) -> bool {
    let a24 = f.mul(f.add(a_new, 2), f.inv(4));
    ladder(f, &Xz { x, z: 1 }, &BigUint::from(n), &a24).z == 0
}

// ---------------------------------------------------------------------------
// Mode: rational kernel (λ = 1), Vélu vs √élu.

fn mode_rational(max_log2: f64, seeds: u64) -> Value {
    let mut ells = Vec::new();
    let mut e = 3.0;
    while e <= max_log2 + 1e-9 {
        let l = next_prime(2f64.powf(e).round() as u64);
        if l != 7 && !ells.contains(&l) {
            ells.push(l);
        }
        e += 0.5;
    }
    let mut rows = Vec::new();
    for &ell in &ells {
        for seed in 0..seeds {
            let mut rng = StdRng::seed_from_u64(0x163 ^ (ell << 8) ^ seed);
            // Restrict to ℓ² ∥ N so that E(F_p)[ℓ^∞] = E[ℓ] (uniform lines).
            let mut rejected = 0u64;
            let inst = loop {
                let inst = make_instance(ell, 1, &mut rng);
                if val(&BigUint::from(inst.n), ell).0 == 2 {
                    break inst;
                }
                rejected += 1;
            };
            let f = Fp { p: inst.p };
            // Evaluation points: two random points on the order-N side.
            let mut pts = Vec::new();
            while pts.len() < 2 {
                let x = rng.gen_range(2..inst.p);
                let rhs = f.mul(x, f.add(f.mul(x, x), f.add(f.mul(inst.a, x), 1)));
                if f.legendre(rhs) == inst.chi_n {
                    pts.push(x);
                }
            }
            // Kernel finding + Vélu, retried if the line is horizontal.
            let c0 = Cnt::now();
            let w0 = Instant::now();
            let mut find_cnt = Cnt::default();
            let mut wasted = Cnt::default();
            let mut horizontal_retries = 0;
            let (xp, tries, velu_out, velu_cnt, velu_wall) = loop {
                let cf = Cnt::now();
                let (xp, tries) = find_rational_kernel(&inst, &mut rng);
                find_cnt = find_cnt.plus(Cnt::since(cf));
                let cv = Cnt::now();
                let wv = Instant::now();
                let out = velu(&f, inst.a, inst.a24, ell, xp, &pts[..1]);
                let vc = Cnt::since(cv);
                let vw = wv.elapsed().as_secs_f64();
                let jn = j_of_a(&f, out.a_new);
                if jn == inst.jcrater {
                    horizontal_retries += 1;
                    wasted = wasted.plus(vc);
                    continue;
                }
                break (xp, tries, out, vc, vw);
            };
            let total_find_velu = Cnt::since(c0);
            let _ = w0;
            // Marginal cost of a second image point with Vélu.
            let cv2 = Cnt::now();
            let velu2 = velu(&f, inst.a, inst.a24, ell, xp, &pts);
            let velu2_cnt = Cnt::since(cv2);
            // √élu.
            let (sq_out, sq_cnt, sq_wall, sq_meta, sq2_cnt) = if ell >= 11 {
                let cs = Cnt::now();
                let ws = Instant::now();
                let (o, meta) = sqrt_velu(&f, inst.a, inst.a24, ell, xp, &pts[..1]);
                let sc = Cnt::since(cs);
                let sw = ws.elapsed().as_secs_f64();
                let cs2 = Cnt::now();
                let (o2, _) = sqrt_velu(&f, inst.a, inst.a24, ell, xp, &pts);
                let sc2 = Cnt::since(cs2);
                assert_eq!(o2.images, velu2.images, "√élu images ≠ Vélu images");
                (Some(o), Some(sc), Some(sw), meta, Some(sc2))
            } else {
                (None, None, None, Value::Null, None)
            };
            // Validation (uncounted).
            let jn = j_of_a(&f, velu_out.a_new);
            let mut checks = vec![];
            if let Some(o) = &sq_out {
                assert_eq!(o.a_new, velu_out.a_new, "√élu codomain ≠ Vélu codomain");
                assert_eq!(o.images, velu_out.images, "√élu image ≠ Vélu image");
                checks.push("sqrtvelu_eq_velu");
            }
            for (k, &x) in velu2.images.iter().enumerate() {
                assert!(
                    order_check(&f, velu2.a_new, x, inst.n),
                    "image {k} has wrong order"
                );
            }
            checks.push("images_order_N_on_codomain");
            let mut lib_j = Value::Null;
            if ell <= 151 {
                let hm = kernel_poly_from_points(&f, inst.a24, xp, ell);
                let hw = to_weierstrass_kernel(&f, &hm, inst.alpha, inst.s);
                let lj = lib_certify(inst.p, &hw, ell).expect("library rejects kernel");
                assert_eq!(lj, jn, "library Vélu j ≠ Montgomery j");
                lib_j = json!(lj);
                checks.push("lib_verify_kernel_and_j");
            }
            let row = json!({
                "ell": ell, "log2_ell": (ell as f64).log2(), "seed": seed, "lambda": inst.lambda, "instances_rejected": rejected,
                "p": inst.p, "t": inst.t, "v": inst.v, "N": inst.n.to_string(),
                "v_ell_N": val(&BigUint::from(inst.n), ell).0,
                "find": find_cnt.js(), "find_point_tries": tries,
                "horizontal_retries": horizontal_retries, "wasted_velu": wasted.js(),
                "velu_1pt": velu_cnt.js(), "velu_2pt": velu2_cnt.js(),
                "velu_wall_s_contended": velu_wall,
                "find_plus_velu_total": total_find_velu.js(),
                "sqrtvelu_1pt": sq_cnt.map(|c| c.js()), "sqrtvelu_2pt": sq2_cnt.map(|c| c.js()),
                "sqrtvelu_wall_s_contended": sq_wall, "sqrtvelu_sets": sq_meta,
                "codomain_j": jn, "crater_j": inst.jcrater, "lib_codomain_j": lib_j,
                "checks": checks,
            });
            eprintln!(
                "rational ℓ={ell} seed={seed} find={} velu={} sqrtvelu={:?}",
                find_cnt.meq(),
                velu_cnt.meq(),
                sq_cnt.map(|c| c.meq())
            );
            rows.push(row);
        }
    }
    json!({"mode": "rational", "rows": rows})
}

// ---------------------------------------------------------------------------
// Mode: generic kernel (λ primitive), torsion point over F_{p^k}.

fn find_irreducible(f: &Fp, k: usize, rng: &mut StdRng) -> (Fpk, u64) {
    // Ben-Or test on X^k + X^m + c.
    let mut cands = 0;
    loop {
        cands += 1;
        let m = if k > 1 { rng.gen_range(1..k) } else { 0 };
        let c = rng.gen_range(1..f.p);
        let fk = Fpk { f: *f, k, m, c };
        if k == 1 {
            return (fk, cands);
        }
        let modp = fk.modulus();
        let x = fk.x();
        let pb = BigUint::from(f.p);
        let mut xq = x.clone();
        let mut ok = true;
        for _ in 1..=k / 2 {
            xq = fk.pow_big(&xq, &pb);
            let g = poly_gcd(f, &modp, &trim(psub(f, &xq, &x)));
            if g.len() > 1 {
                ok = false;
                break;
            }
        }
        if ok {
            return (fk, cands);
        }
    }
}

/// Berlekamp–Massey: connection polynomial of `s`.
fn berlekamp_massey(f: &Fp, s: &[u64]) -> Vec<u64> {
    let mut c = vec![1u64];
    let mut b = vec![1u64];
    let mut l = 0usize;
    let mut m = 1usize;
    let mut bb = 1u64;
    for n in 0..s.len() {
        let mut d = s[n];
        for i in 1..=l {
            d = f.add(d, f.mul(*c.get(i).unwrap_or(&0), s[n - i]));
        }
        if d == 0 {
            m += 1;
            continue;
        }
        let coef = f.mul(d, f.inv(bb));
        let mut nc = c.clone();
        nc.resize(nc.len().max(b.len() + m), 0);
        for (i, &bi) in b.iter().enumerate() {
            nc[i + m] = f.sub(nc[i + m], f.mul(coef, bi));
        }
        if 2 * l <= n {
            l = n + 1 - l;
            b = c;
            bb = d;
            m = 1;
        } else {
            m += 1;
        }
        c = nc;
    }
    c.resize(l + 1, 0);
    c
}

fn mode_generic(max_ell: u64, cap_s: f64) -> Value {
    let mut rows = Vec::new();
    let t_start = Instant::now();
    let mut ell = 5u64;
    while ell <= max_ell {
        if !is_prime_u64(ell) || ell == 7 {
            ell += 2;
            continue;
        }
        if t_start.elapsed().as_secs_f64() > cap_s {
            eprintln!("generic: time cap reached before ℓ={ell}");
            break;
        }
        let mut rng = StdRng::seed_from_u64(0x1630 ^ ell ^ seed_salt());
        let lambda = primitive_root(ell);
        let k = ((ell - 1) / 2) as usize;
        // Restrict to instances whose ℓ-Sylow over F_(p^k) is exactly E[ℓ]
        // (v_ℓ = 2 on the twist, 0 on E): then a cofactor multiple is a
        // uniform point of E[ℓ] and every line is equally likely.
        let mut rejected = 0u64;
        let (inst, n_e, n_t, ve, vt) = loop {
            let inst = make_instance(ell, lambda, &mut rng);
            let (n_e, n_t) = orders_over_ext(&inst, k);
            let (ve, _) = val(&n_e, ell);
            let (vt, _) = val(&n_t, ell);
            if ve == 0 && vt == 2 {
                break (inst, n_e, n_t, ve, vt);
            }
            rejected += 1;
        };
        let f = Fp { p: inst.p };
        let wall = Instant::now();
        // Extension setup.
        let c_setup = Cnt::now();
        let (fk, cands) = find_irreducible(&f, k, &mut rng);
        let setup = Cnt::since(c_setup);
        let (ntarget, side) = if vt >= 2 { (n_t, "twist") } else { (n_e, "E") };
        let (e, m) = val(&ntarget, ell);
        assert!(e >= 2, "ℓ² must divide the order over F_(p^k)");
        // Lift A, a24 into F_{p^k}.
        let mut a24k = fk.zero();
        a24k[0] = inst.a24;
        let lb = BigUint::from(ell);
        // Evaluation with the kernel polynomial over F_p (O(ℓ) per point).
        let mut pts = Vec::new();
        while pts.len() < 2 {
            let x = rng.gen_range(2..inst.p);
            let rhs = f.mul(x, f.add(f.mul(x, x), f.add(f.mul(inst.a, x), 1)));
            if f.legendre(rhs) == inst.chi_n {
                pts.push(x);
            }
        }
        let horner =
            |h: &[u64], x: u64| h.iter().rev().fold(0u64, |acc, c| f.add(f.mul(acc, x), *c));
        let mut horizontal_retries = 0u64;
        let mut wasted = Cnt::default();
        let mut attempts = 0;
        let mut att_final: u64;
        let (ladder_cnt, minpoly_cnt, find_cnt, find_wall, minpoly, _x1, a_new, codomain_cnt) = loop {
            // Kernel finding: ladder by the cofactor, then down to order ℓ.
            let c_find = Cnt::now();
            att_final = 0;
            let q = 'outer: loop {
                attempts += 1;
                att_final += 1;
                let x: Vec<u64> = (0..k).map(|_| rng.gen_range(0..inst.p)).collect();
                let mut q = ladder(&fk, &Xz { x, z: fk.one() }, &m, &a24k);
                if fk.is_zero(&q.z) {
                    continue;
                }
                for _ in 0..e {
                    let r = ladder(&fk, &q, &lb, &a24k);
                    if fk.is_zero(&r.z) {
                        break 'outer q;
                    }
                    q = r;
                }
            };
            let ladder_cnt = Cnt::since(c_find);
            let c_mp = Cnt::now();
            let x1 = Fld::mul(&fk, &q.x, &fk.inv(&q.z));
            // Minimal polynomial of x1 by Berlekamp–Massey on a random projection.
            let proj: Vec<u64> = (0..k).map(|_| rng.gen_range(0..inst.p)).collect();
            let mut pw = fk.one();
            let mut seq = Vec::with_capacity(2 * k);
            for _ in 0..2 * k {
                let s = pw
                    .iter()
                    .zip(&proj)
                    .fold(0u64, |acc, (a, b)| f.add(acc, f.mul(*a, *b)));
                seq.push(s);
                pw = Fld::mul(&fk, &pw, &x1);
            }
            let conn = berlekamp_massey(&f, &seq);
            let minpoly: Vec<u64> = conn.iter().rev().copied().collect(); // monic, deg L
            let minpoly_cnt = Cnt::since(c_mp);
            let find_cnt = Cnt::since(c_find);
            let find_wall = wall.elapsed().as_secs_f64();
            assert_eq!(minpoly.len() - 1, k, "minimal polynomial degree ≠ (ℓ−1)/2");
            // Validate: minpoly(x1) = 0 in F_{p^k}.
            let mut acc = fk.zero();
            for &c in minpoly.iter().rev() {
                acc = Fld::mul(&fk, &acc, &x1);
                acc[0] = f.add(acc[0], c);
            }
            assert!(fk.is_zero(&acc), "minpoly does not vanish at x1");
            let c_ev = Cnt::now();
            let h1 = horner(&minpoly, 1);
            let hm1 = horner(&minpoly, f.neg(1));
            // Π(x_s + 1) = (−1)^k h(−1), Π(x_s − 1) = (−1)^k h(1); signs vanish.
            let a_new = edwards_codomain(&f, inst.a, ell, hm1, h1);
            let codomain_cnt = Cnt::since(c_ev);
            // A horizontal line (codomain j = crater j) is not a vertical
            // transfer: count its cost as wasted and draw another point.
            if j_of_a(&f, a_new) == inst.jcrater {
                horizontal_retries += 1;
                wasted = wasted.plus(Cnt::since(c_find));
                continue;
            }
            break (
                ladder_cnt,
                minpoly_cnt,
                find_cnt,
                find_wall,
                minpoly,
                x1,
                a_new,
                codomain_cnt,
            );
        };
        let c_pt = Cnt::now();
        let rev: Vec<u64> = minpoly.iter().rev().copied().collect();
        let mut images = Vec::new();
        for &x in &pts {
            let den = horner(&minpoly, x);
            let num = horner(&rev, x);
            images.push(f.mul(f.mul(x, f.mul(num, num)), f.inv(f.mul(den, den))));
        }
        let eval_cnt = Cnt::since(c_pt);
        // Checks (uncounted).
        let jn = j_of_a(&f, a_new);
        for &x in &images {
            assert!(order_check(&f, a_new, x, inst.n), "image has wrong order");
        }
        let mut lib_j = Value::Null;
        let mut checks = vec!["minpoly_vanishes", "images_order_N_on_codomain"];
        if ell <= 211 {
            let hw = to_weierstrass_kernel(&f, &minpoly, inst.alpha, inst.s);
            let lj = lib_certify(inst.p, &hw, ell).expect("library rejects kernel");
            assert_eq!(lj, jn, "library j ≠ Montgomery j");
            lib_j = json!(lj);
            checks.push("lib_verify_kernel_and_j");
        }
        let vertical = jn != inst.jcrater;
        eprintln!(
            "generic ℓ={ell} k={k} setup={} ladder={} minpoly={} eval/pt={} vertical={vertical} wall={find_wall:.2}s",
            setup.meq(),
            ladder_cnt.meq(),
            minpoly_cnt.meq(),
            eval_cnt.meq() / 2
        );
        rows.push(json!({
            "ell": ell, "k": k, "lambda": lambda, "p": inst.p, "t": inst.t, "v": inst.v,
            "side": side, "v_ell_E": ve, "v_ell_twist": vt, "instances_rejected": rejected,
            "irreducible_candidates": cands, "ext_modulus": [k, fk.m, fk.c],
            "setup": setup.js(), "ladder": ladder_cnt.js(), "ladder_attempts_total": attempts, "ladder_attempts_final_find": att_final, "find_wall_includes": "setup and horizontal retries",
            "cofactor_bits": m.bits(), "minpoly": minpoly_cnt.js(), "find": find_cnt.js(),
            "find_incl_setup": find_cnt.plus(setup).js(),
            "find_wall_s_contended": find_wall,
            "codomain_from_h": codomain_cnt.js(), "eval_2pt_from_h": eval_cnt.js(),
            "codomain_j": jn, "crater_j": inst.jcrater, "vertical": vertical,
            "horizontal_retries": horizontal_retries, "wasted_on_horizontal": wasted.js(),
            "lib_codomain_j": lib_j, "checks": checks,
        }));
        ell += 2;
    }
    json!({"mode": "generic", "rows": rows})
}

// ---------------------------------------------------------------------------
// Mode: the library's Elkies route (contended wall time only).

fn mode_elkies(max_ell: u64, cap_s: f64) -> Value {
    let mut rows = Vec::new();
    let t_start = Instant::now();
    let mut ell = 3u64;
    while ell <= max_ell {
        if !is_prime_u64(ell) || ell == 7 {
            ell += 2;
            continue;
        }
        if t_start.elapsed().as_secs_f64() > cap_s {
            eprintln!("elkies: time cap reached before ℓ={ell}");
            break;
        }
        let mut rng = StdRng::seed_from_u64(0x1631 ^ ell);
        let lambda = primitive_root(ell);
        let inst = make_instance(ell, lambda, &mut rng);
        let lf = lib_field(inst.p);
        let e = Model {
            a: lf.from_i64(-35),
            b: lf.from_i64(-98),
        };
        let j = e.j(&lf).unwrap();
        let w0 = Instant::now();
        let phi = modpoly::modular_polynomial(&lf, ell).expect("Φ_ℓ");
        let t_phi = w0.elapsed().as_secs_f64();
        let w1 = Instant::now();
        let py = phi.at_x(&lf, &j);
        let roots = lpoly::roots(&lf, &py, 1);
        let t_roots = w1.elapsed().as_secs_f64();
        let others: Vec<_> = roots.iter().filter(|r| **r != j).cloned().collect();
        let w2 = Instant::now();
        let mut ok = 0;
        let mut errs = Vec::new();
        let mut t_first = None;
        for r in &others {
            match kernel::elkies_isogeny(&lf, &e, &phi, r) {
                Ok(_) => {
                    ok += 1;
                    if t_first.is_none() {
                        t_first = Some(w2.elapsed().as_secs_f64());
                    }
                }
                Err(er) => errs.push(format!("{er:?}")),
            }
        }
        let t_all = w2.elapsed().as_secs_f64();
        eprintln!(
            "elkies ℓ={ell} phi={t_phi:.3}s roots={t_roots:.3}s n_roots={} elkies_ok={ok}/{}",
            roots.len(),
            others.len()
        );
        rows.push(json!({
            "ell": ell, "p": inst.p, "t": inst.t, "v": inst.v,
            "phi_wall_s_contended": t_phi, "roots_wall_s_contended": t_roots,
            "n_distinct_rational_roots": roots.len(), "n_non_crater_roots": others.len(),
            "elkies_first_kernel_wall_s_contended": t_first,
            "elkies_all_kernels_wall_s_contended": t_all,
            "elkies_ok": ok, "elkies_errors": errs,
            "phi_checks_passed": phi.checks_passed,
        }));
        ell += 2;
    }
    json!({"mode": "elkies", "rows": rows})
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        eprintln!(
            "usage: {} rational|generic|elkies <max> <seeds|cap_s> <out.json>",
            args[0]
        );
        std::process::exit(2);
    }
    let t0 = Instant::now();
    let mut out = match args[1].as_str() {
        "rational" => mode_rational(args[2].parse().unwrap(), args[3].parse().unwrap()),
        "generic" => mode_generic(args[2].parse().unwrap(), args[3].parse().unwrap()),
        "elkies" => mode_elkies(args[2].parse().unwrap(), args[3].parse().unwrap()),
        "lines" => mode_lines(args[2].parse().unwrap(), args[3].parse().unwrap()),
        "analyze" => mode_analyze(&args[2], &args[3]),
        m => panic!("unknown mode {m}"),
    };
    out["argv"] = json!(args);
    out["unit"] = json!("M = F_p multiplications+squarings (incl. inside polynomial/extension arithmetic); I = F_p inversions; Meq = M + 91 I");
    out["field_bits"] = json!(61);
    out["total_wall_s_contended"] = json!(t0.elapsed().as_secs_f64());
    let s = serde_json::to_string_pretty(&out).unwrap();
    std::fs::write(&args[4], s).expect("write output");
}

/// Optional seed salt (`R4_SEED_SALT`) for repeated generic runs.
fn seed_salt() -> u64 {
    std::env::var("R4_SEED_SALT")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .map(|s| s << 32)
        .unwrap_or(0)
}

/// Diagnostic: distribution of codomain j over sampled kernel points.
fn mode_lines(ell: u64, samples: u64) -> Value {
    let mut rng = StdRng::seed_from_u64(0x1632 ^ ell ^ seed_salt());
    let inst = make_instance(ell, 1, &mut rng);
    let v_ell_n = val(&BigUint::from(inst.n), ell).0;
    let f = Fp { p: inst.p };
    let mut hist: std::collections::BTreeMap<u64, u64> = Default::default();
    for _ in 0..samples {
        let (xp, _) = find_rational_kernel(&inst, &mut rng);
        let out = velu(&f, inst.a, inst.a24, ell, xp, &[]);
        *hist.entry(j_of_a(&f, out.a_new)).or_default() += 1;
    }
    let crater = hist.get(&inst.jcrater).copied().unwrap_or(0);
    eprintln!(
        "lines ℓ={ell} p={} v_ell_N={v_ell_n} distinct_j={} crater_hits={crater}/{samples}",
        inst.p,
        hist.len()
    );
    json!({"mode": "lines", "ell": ell, "p": inst.p, "v_ell_N": v_ell_n, "crater_hits": crater,
           "samples": samples, "hist": hist.values().collect::<Vec<_>>()})
}

/// `(#E(F_(p^k)), #E^twist(F_(p^k)))` from the trace recurrence.
fn orders_over_ext(inst: &Inst, k: usize) -> (BigUint, BigUint) {
    let pb = BigInt::from(inst.p);
    let (mut t0, mut t1) = (BigInt::from(2), BigInt::from(inst.t));
    for _ in 1..k {
        let t2 = BigInt::from(inst.t) * &t1 - &pb * &t0;
        t0 = t1;
        t1 = t2;
    }
    let pk = pb.pow(k as u32);
    let one = BigInt::from(1);
    let n_e = (&pk + &one - &t1).abs().to_biguint().unwrap();
    let n_t = (&pk + &one + &t1).abs().to_biguint().unwrap();
    (n_e, n_t)
}

// ---------------------------------------------------------------------------
// Mode: fits and the derived crossover table (reads the recorded JSON).

/// Two-sided 95% Student-t quantile (Cornish–Fisher expansion in 1/df).
fn t975(df: f64) -> f64 {
    let z: f64 = 1.959964;
    z + (z.powi(3) + z) / (4.0 * df)
        + (5.0 * z.powi(5) + 16.0 * z.powi(3) + 3.0 * z) / (96.0 * df * df)
        + (3.0 * z.powi(7) + 19.0 * z.powi(5) + 17.0 * z.powi(3) - 15.0 * z) / (384.0 * df.powi(3))
}

#[derive(Clone, Copy, Debug)]
struct Fit {
    slope: f64,
    slope_lo: f64,
    slope_hi: f64,
    /// log2 of the constant: log2 cost ≈ icpt + slope·log2 ℓ.
    icpt: f64,
    n: usize,
    r2: f64,
}

impl Fit {
    fn eval_log2(&self, log2_l: f64) -> f64 {
        self.icpt + self.slope * log2_l
    }
    fn js(&self, name: &str, range: &str) -> Value {
        json!({"series": name, "range": range, "exponent": self.slope,
               "exponent_ci95": [self.slope_lo, self.slope_hi],
               "log2_constant": self.icpt, "n": self.n, "r2": self.r2})
    }
}

/// OLS of log2 y on log2 x.
fn fit_loglog(pts: &[(f64, f64)]) -> Fit {
    let n = pts.len();
    assert!(n >= 3, "need at least three points to fit");
    let xs: Vec<f64> = pts.iter().map(|p| p.0.log2()).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.1.log2()).collect();
    let mx = xs.iter().sum::<f64>() / n as f64;
    let my = ys.iter().sum::<f64>() / n as f64;
    let sxx: f64 = xs.iter().map(|x| (x - mx).powi(2)).sum();
    let sxy: f64 = xs.iter().zip(&ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let slope = sxy / sxx;
    let icpt = my - slope * mx;
    let sse: f64 = xs
        .iter()
        .zip(&ys)
        .map(|(x, y)| (y - icpt - slope * x).powi(2))
        .sum();
    let syy: f64 = ys.iter().map(|y| (y - my).powi(2)).sum();
    let se = (sse / (n as f64 - 2.0) / sxx).sqrt();
    let t = t975(n as f64 - 2.0);
    Fit {
        slope,
        slope_lo: slope - t * se,
        slope_hi: slope + t * se,
        icpt,
        n,
        r2: 1.0 - sse / syy,
    }
}

fn meq(v: &Value) -> f64 {
    v["Meq"].as_f64().unwrap()
}

/// Mean of `get(row)` per distinct ℓ.
fn per_ell(rows: &[Value], get: impl Fn(&Value) -> Option<f64>) -> Vec<(f64, f64)> {
    let mut m: std::collections::BTreeMap<u64, (f64, f64)> = Default::default();
    for r in rows {
        if let Some(y) = get(r) {
            let e = m.entry(r["ell"].as_u64().unwrap()).or_insert((0.0, 0.0));
            e.0 += y;
            e.1 += 1.0;
        }
    }
    m.into_iter().map(|(l, (s, c))| (l as f64, s / c)).collect()
}

/// Smallest log2 ℓ in [lo, hi] with cost(log2 ℓ) ≥ target (bisection on a
/// non-decreasing model); None if cost(hi) < target.
fn crossover(cost: impl Fn(f64) -> f64, target: f64, lo: f64, hi: f64) -> Option<f64> {
    if cost(hi) < target {
        return None;
    }
    if cost(lo) >= target {
        return Some(lo);
    }
    let (mut a, mut b) = (lo, hi);
    for _ in 0..200 {
        let m = 0.5 * (a + b);
        if cost(m) >= target {
            b = m;
        } else {
            a = m;
        }
    }
    Some(b)
}

fn mode_analyze(dir: &str, out: &str) -> Value {
    let load = |name: &str| -> Value {
        let s = std::fs::read_to_string(format!("{dir}/{name}")).expect("read json");
        serde_json::from_str(&s).expect("parse json")
    };
    let rat = load("rational.json");
    let gen = load("generic.json");
    let elk = load("elkies.json");
    let rr = rat["rows"].as_array().unwrap().clone();
    let gr = gen["rows"].as_array().unwrap().clone();
    let er = elk["rows"].as_array().unwrap().clone();
    let mut text = String::new();
    let mut fits = Vec::new();

    // Rational case.
    let find_r = per_ell(&rr, |r| Some(meq(&r["find"])));
    let velu = per_ell(&rr, |r| Some(meq(&r["velu_1pt"])));
    let velu_marg = per_ell(&rr, |r| Some(meq(&r["velu_2pt"]) - meq(&r["velu_1pt"])));
    let sq = per_ell(&rr, |r| {
        r["sqrtvelu_1pt"].get("Meq").and_then(|v| v.as_f64())
    });
    let sq_marg = per_ell(&rr, |r| {
        let a = r["sqrtvelu_2pt"].get("Meq")?.as_f64()?;
        let b = r["sqrtvelu_1pt"].get("Meq")?.as_f64()?;
        Some(a - b)
    });
    let sel = |v: &[(f64, f64)], lo: f64| -> Vec<(f64, f64)> {
        v.iter().filter(|p| p.0 >= lo).copied().collect()
    };
    let f_find_r = fit_loglog(&find_r);
    let f_velu = fit_loglog(&sel(&velu, 1024.0));
    let f_velu_m = fit_loglog(&sel(&velu_marg, 1024.0));
    let f_sq_all = fit_loglog(&sel(&sq, 1024.0));
    let f_sq_top = fit_loglog(&sel(&sq, 65536.0));
    let f_sq_m = fit_loglog(&sel(&sq_marg, 1024.0));
    fits.push(f_find_r.js("rational_find (cofactor ladder)", "all ell"));
    fits.push(f_velu.js("velu_1pt (codomain + 1 image)", "ell >= 2^10"));
    fits.push(f_velu_m.js("velu_marginal_per_extra_point", "ell >= 2^10"));
    fits.push(f_sq_all.js("sqrtvelu_1pt (codomain + 1 image)", "ell >= 2^10"));
    fits.push(f_sq_top.js("sqrtvelu_1pt (codomain + 1 image)", "ell >= 2^16"));
    fits.push(f_sq_m.js("sqrtvelu_marginal_per_extra_point", "ell >= 2^10"));
    let find_r_mean = find_r.iter().map(|p| p.1).sum::<f64>() / find_r.len() as f64;
    let find_r_max = find_r.iter().map(|p| p.1).fold(0.0, f64::max);
    let (l_top, sq_top) = *sq.last().unwrap();

    // Generic case.
    let ladder_pa = per_ell(&gr, |r| {
        Some(meq(&r["ladder"]) / r["ladder_attempts_final_find"].as_f64().unwrap())
    });
    let minp = per_ell(&gr, |r| Some(meq(&r["minpoly"])));
    let setup = per_ell(&gr, |r| Some(meq(&r["setup"])));
    // M only: each point also costs one inversion (Kohel/Costello–Hisil form).
    let evalh = per_ell(&gr, |r| r["eval_2pt_from_h"]["M"].as_f64().map(|m| m / 2.0));
    let gfind: Vec<(f64, f64)> = ladder_pa
        .iter()
        .zip(&minp)
        .map(|(a, b)| (a.0, 2.0 * a.1 + b.1))
        .collect();
    let gfind_meas = per_ell(&gr, |r| Some(meq(&r["find"])));
    let f_lpa = fit_loglog(&sel(&ladder_pa, 30.0));
    let f_minp = fit_loglog(&sel(&minp, 30.0));
    let f_setup = fit_loglog(&sel(&setup, 30.0));
    let f_evalh = fit_loglog(&sel(&evalh, 30.0));
    let f_gfind = fit_loglog(&sel(&gfind, 30.0));
    let f_gfind_meas = fit_loglog(&sel(&gfind_meas, 30.0));
    let f_gfind_top = fit_loglog(&sel(&gfind, 150.0));
    fits.push(f_lpa.js(
        "generic_ladder_per_attempt (F_p^k, k=(ell-1)/2)",
        "ell >= 30",
    ));
    fits.push(f_minp.js("generic_minpoly", "ell >= 30"));
    fits.push(f_gfind.js(
        "generic_find_expected (2 x ladder/attempt + minpoly)",
        "ell >= 30",
    ));
    fits.push(f_gfind_top.js(
        "generic_find_expected (2 x ladder/attempt + minpoly)",
        "ell >= 150",
    ));
    fits.push(f_gfind_meas.js("generic_find_measured (realised attempts)", "ell >= 30"));
    fits.push(f_setup.js(
        "generic_ext_setup (irreducible search, Ben-Or)",
        "ell >= 30",
    ));
    fits.push(f_evalh.js(
        "generic_eval_per_point_from_kernel_poly (M; +1 I)",
        "ell >= 30",
    ));
    let (lg_top, gfind_top) = *gfind.last().unwrap();

    // Elkies route (library, contended wall time).
    let phi_w: Vec<(f64, f64)> = er
        .iter()
        .filter(|r| r["ell"].as_f64().unwrap() >= 11.0)
        .map(|r| {
            (
                r["ell"].as_f64().unwrap(),
                r["phi_wall_s_contended"].as_f64().unwrap(),
            )
        })
        .collect();
    let f_phi = fit_loglog(&phi_w);
    fits.push(f_phi.js("library_Phi_ell_mod_p_wall_seconds_contended", "ell >= 11"));

    text.push_str("Fits: log2(cost) = log2_constant + exponent * log2(ell), OLS, 95% t-interval\n");
    for f in &fits {
        text.push_str(&format!(
            "  {:<62} {:<11} exponent {:.3} [{:.3}, {:.3}]  log2 c {:.2}  n={} r2={:.4}\n",
            f["series"].as_str().unwrap(),
            f["range"].as_str().unwrap(),
            f["exponent"].as_f64().unwrap(),
            f["exponent_ci95"][0].as_f64().unwrap(),
            f["exponent_ci95"][1].as_f64().unwrap(),
            f["log2_constant"].as_f64().unwrap(),
            f["n"].as_u64().unwrap(),
            f["r2"].as_f64().unwrap()
        ));
    }

    // Derived crossover table.
    let mut table = Vec::new();
    text.push_str("\nDerived (extrapolated) crossover: smallest ell with find+evaluate >= g * sqrt(pi r / 4) F_p mults\n");
    for &g in &[6.0f64, 12.0] {
        for &rb in &[64.0f64, 128.0, 256.0] {
            let log2_rho = (g * (std::f64::consts::PI / 4.0).sqrt()).log2() + rb / 2.0;
            let s = rb / 61.0; // log p scaling for ladder-type phases (p ≈ r)
            let log2_lmax = rb / 2.0 + 1.0; // ℓ² | t² − 4p ⇒ ℓ ≤ 2√p
            let hi = 400.0;
            let rho = 2f64.powf(log2_rho);
            // Rational kernel, Vélu (measured fit).
            let c_velu = |x: f64| s * find_r_max + 2f64.powf(f_velu.eval_log2(x));
            // Rational kernel, √élu measured power law (ℓ ≥ 2^16 fit).
            let c_sq_fit = |x: f64| s * find_r_max + 2f64.powf(f_sq_top.eval_log2(x));
            // Rational kernel, √élu asymptotic model √ℓ·(log2 ℓ)², calibrated
            // at the largest measured ℓ.
            let k_sq = sq_top / (l_top.sqrt() * l_top.log2().powi(2));
            let c_sq_asym = |x: f64| s * find_r_max + k_sq * 2f64.powf(x / 2.0) * x * x;
            // Generic kernel: measured fit (find, scaled by log p) + Kohel eval.
            let c_gen = |x: f64| {
                s * 2f64.powf(f_gfind.eval_log2(x)) + 2f64.powf(f_evalh.eval_log2(x)) + 91.0
            };
            let c_gen_lo = |x: f64| {
                // the same, with the lower end of the exponent interval
                let lg = f_gfind.icpt
                    + f_gfind.slope_lo * x
                    + (f_gfind.slope - f_gfind.slope_lo) * lg_top.log2();
                s * 2f64.powf(lg) + 2f64.powf(f_evalh.eval_log2(x)) + 91.0
            };
            // Generic kernel, asymptotic model ℓ²·log2 ℓ (fast F_(p^k)
            // arithmetic), calibrated at the largest measured ℓ.
            let k_g = gfind_top / (lg_top * lg_top * lg_top.log2());
            let c_gen_asym =
                |x: f64| s * k_g * 2f64.powf(2.0 * x) * x + 2f64.powf(f_evalh.eval_log2(x)) + 91.0;
            let row = |name: &str, c: &dyn Fn(f64) -> f64| -> Value {
                let x = crossover(c, rho, 1.0, hi);
                let note = match x {
                    Some(v) if v > log2_lmax => {
                        "beyond ell_max = 2 sqrt(p): no vertical gap that large exists"
                    }
                    Some(_) => "crossover inside the possible range",
                    None => "no crossover below 2^400",
                };
                json!({"model": name, "log2_ell_star": x, "note": note})
            };
            let rows = vec![
                row("rational kernel, Velu (fit ell>=2^10)", &c_velu),
                row(
                    "rational kernel, sqrt-elu (power-law fit ell>=2^16)",
                    &c_sq_fit,
                ),
                row(
                    "rational kernel, sqrt-elu (model sqrt(l) log^2 l, calibrated)",
                    &c_sq_asym,
                ),
                row("generic kernel k=(l-1)/2 (fit ell>=30, central)", &c_gen),
                row(
                    "generic kernel k=(l-1)/2 (fit, exponent at CI low end)",
                    &c_gen_lo,
                ),
                row(
                    "generic kernel k=(l-1)/2 (model l^2 log l, calibrated)",
                    &c_gen_asym,
                ),
            ];
            text.push_str(&format!(
                "  g={g} M/op, r=2^{rb}: rho = 2^{log2_rho:.2} M, ell_max = 2^{log2_lmax:.0}\n"
            ));
            for r in &rows {
                text.push_str(&format!(
                    "    {:<66} log2 ell* = {:>8}  ({})\n",
                    r["model"].as_str().unwrap(),
                    r["log2_ell_star"]
                        .as_f64()
                        .map(|v| format!("{v:.2}"))
                        .unwrap_or("-".into()),
                    r["note"].as_str().unwrap()
                ));
            }
            table.push(
                json!({"g_mults_per_group_op": g, "log2_r": rb, "log2_rho_mults": log2_rho,
                              "log2_ell_max": log2_lmax, "logp_scale": s, "rows": rows}),
            );
        }
    }

    // Intermediate k (kernel over F_(p^k) with k ≪ ℓ): derived formula.
    // find(k) from the generic fit evaluated at ℓ = 2k + 1; evaluation by
    // √élu over F_(p^k): √élu(ℓ) × m(k), m(k) = F_(p^k) mult cost in F_p mults
    // (Karatsuba with threshold 24 + trinomial reduction, measured below).
    let mut mk = Vec::new();
    for &k in &[1usize, 2, 4, 8, 16, 32, 64] {
        let fk = Fpk {
            f: Fp {
                p: 2305843009213693951,
            },
            k,
            m: if k > 1 { 1 } else { 0 },
            c: 3,
        };
        let a: Vec<u64> = (0..k as u64).map(|i| i + 2).collect();
        let c0 = Cnt::now();
        let _ = Fld::mul(&fk, &a, &a);
        mk.push((k, Cnt::since(c0).m as f64));
    }
    let mut kt = Vec::new();
    text.push_str(
        "\nDerived: kernel over F_(p^k), r = 2^128, g = 6 (find(k) + sqrt-elu(ell) x m(k))\n",
    );
    let log2_rho128 = (6.0 * (std::f64::consts::PI / 4.0).sqrt()).log2() + 64.0;
    for &(k, mkv) in &mk {
        let findk = if k == 1 {
            find_r_max
        } else {
            2f64.powf(f_gfind.eval_log2(((2 * k + 1) as f64).log2()))
        };
        let s = 128.0 / 61.0;
        let c = |x: f64| s * findk + mkv * k_sq_model(sq_top, l_top, x);
        let x = crossover(c, 2f64.powf(log2_rho128), 1.0, 400.0);
        text.push_str(&format!(
            "  k={k:<3} m(k)={mkv:<6} find(k)=2^{:.1} M  log2 ell* = {}\n",
            (s * findk).log2(),
            x.map(|v| format!("{v:.1}")).unwrap_or("-".into())
        ));
        kt.push(
            json!({"k": k, "ext_mult_cost_M": mkv, "log2_find_M_scaled": (s*findk).log2(),
                       "log2_ell_star": x}),
        );
    }
    text.push_str(&format!(
        "\nrational find: mean {find_r_mean:.0} Meq, max {find_r_max:.0} Meq per kernel point (61-bit p)\n"
    ));
    print!("{text}");
    let res = json!({
        "mode": "analyze", "fits": fits, "crossover": table, "intermediate_k_r128_g6": kt,
        "rational_find_mean_Meq": find_r_mean, "rational_find_max_Meq": find_r_max,
        "sqrtvelu_calibration": {"ell": l_top, "Meq": sq_top},
        "generic_calibration": {"ell": lg_top, "Meq": gfind_top},
        "assumptions": {
            "rho_group_ops": "sqrt(pi r / 4)",
            "g": "field mults per group op: 6 (affine add with batched inversion), 12 as sensitivity",
            "logp_scale": "ladder-type phases scale with log2 p / 61, p ~ r; Velu/sqrt-elu counts taken as p-independent",
            "ell_max": "ell^2 | t^2 - 4p implies ell <= 2 sqrt(p)",
        },
        "text": text,
    });
    std::fs::write(out, serde_json::to_string_pretty(&res).unwrap()).expect("write");
    res
}

/// √élu asymptotic model √ℓ·(log2 ℓ)² calibrated at `(l_top, sq_top)`.
fn k_sq_model(sq_top: f64, l_top: f64, log2_l: f64) -> f64 {
    let k = sq_top / (l_top.sqrt() * l_top.log2().powi(2));
    k * 2f64.powf(log2_l / 2.0) * log2_l * log2_l
}
