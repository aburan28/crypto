//! B5a's generators, natively (`rounds/B5a-extension-fields/PROTOCOL.md`):
//! - `instances`: design §4's curves over `GF(p^k)`, searched in this
//!   module's own arithmetic, which shares nothing with the tool;
//! - `cases`: `conformance/v2-b5a/`, every instance re-checked first
//!   (design §4.5's method 3) and every document written from it;
//! - the ICV1 identity of a curve over `GF(p^k)`, computed as
//!   `scripts/curve_id.py`'s reference computes it, for the slugs the
//!   instances record and the cases expect.
//!
//! The declaration's first generators were Python (#1178, closed before
//! anything ran).  This module replaces them: its instance records and its
//! parameter files are theirs byte for byte; only the generator each file
//! names is this one.

use std::cell::OnceCell;
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use num_bigint::BigUint;
use num_traits::ToPrimitive;

use super::bench::Bench;
use super::identity::sha256_hex;
use super::json::{self, obj, J};
use super::runs;
use super::stats;

/// The instances' label: every choice is SHAKE-256 of it, extended.
pub const LABEL: &str = "ic-tool-programme/B5a/instances";
/// A subgroup the index calculus can use lies below this (design §4.2).
const R_MAX: u128 = 1 << 44;
/// This module, as the records name their generator.
pub const GENERATOR: &str = "src/bin/icprog/b5a.rs";

// ── GF(p) and GF(p)[t] ──────────────────────────────────────────────
//
// Coefficients are below `p < 2^63`, so a product of two fits a `u128`.

fn powm(mut b: u128, mut e: u128, p: u128) -> u128 {
    let mut r = 1 % p;
    b %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % p;
        }
        b = b * b % p;
        e >>= 1;
    }
    r
}

/// A polynomial over GF(p), lowest coefficient first, with no zero on top.
type Poly = Vec<u128>;

fn trim(mut a: Poly) -> Poly {
    while a.last() == Some(&0) {
        a.pop();
    }
    a
}

fn pmul(a: &[u128], b: &[u128], p: u128) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0u128; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x != 0 {
            for (j, &y) in b.iter().enumerate() {
                out[i + j] = (out[i + j] + x * y % p) % p;
            }
        }
    }
    trim(out)
}

fn psub(a: &[u128], b: &[u128], p: u128) -> Poly {
    let n = a.len().max(b.len());
    trim(
        (0..n)
            .map(|i| {
                let (x, y) = (
                    a.get(i).copied().unwrap_or(0),
                    b.get(i).copied().unwrap_or(0),
                );
                (x + p - y % p) % p
            })
            .collect(),
    )
}

/// Quotient and remainder by a nonzero `b`.
fn pdivmod(a: &[u128], b: &[u128], p: u128) -> (Poly, Poly) {
    let mut a = trim(a.to_vec());
    let mut q = vec![0u128; (a.len() + 1).saturating_sub(b.len()).max(1)];
    let lead = powm(*b.last().expect("a nonzero divisor"), p - 2, p);
    while a.len() >= b.len() && !a.is_empty() {
        let c = a[a.len() - 1] * lead % p;
        let s = a.len() - b.len();
        q[s] = c;
        for (i, &y) in b.iter().enumerate() {
            a[s + i] = (a[s + i] + p - c * y % p) % p;
        }
        a = trim(a);
    }
    (trim(q), a)
}

fn pgcd(a: &[u128], b: &[u128], p: u128) -> Poly {
    let (mut a, mut b) = (trim(a.to_vec()), trim(b.to_vec()));
    while !b.is_empty() {
        let r = pdivmod(&a, &b, p).1;
        (a, b) = (b, r);
    }
    a
}

/// Rabin's test for the monic `t^k + c_{k-1} t^{k-1} + ... + c_0`
/// (`modc = [c_0, ..., c_{k-1}]`): `t^(p^k) = t` modulo it, and
/// `gcd(t^(p^(k/l)) - t, f) = 1` for every prime `l` dividing `k`.
pub fn rabin_irreducible(p: u128, modc: &[u128]) -> bool {
    let k = modc.len();
    let mut f: Poly = modc.iter().map(|c| c % p).collect();
    f.push(1);
    let t_pow = |mut e: u128| -> Poly {
        let (mut out, mut base): (Poly, Poly) = (vec![1], vec![0, 1]);
        while e > 0 {
            if e & 1 == 1 {
                out = pdivmod(&pmul(&out, &base, p), &f, p).1;
            }
            base = pdivmod(&pmul(&base, &base, p), &f, p).1;
            e >>= 1;
        }
        out
    };
    let Some(pk) = p.checked_pow(k as u32) else {
        return false;
    };
    if !psub(&t_pow(pk), &[0, 1], p).is_empty() {
        return false;
    }
    for ell in (2..=k).filter(|&d| k.is_multiple_of(d) && (2..d).all(|e| !d.is_multiple_of(e))) {
        let x = psub(&t_pow(p.pow((k / ell) as u32)), &[0, 1], p);
        if pgcd(&f, &x, p).len() != 1 {
            return false;
        }
    }
    true
}

fn shake128bits(label: &str) -> u128 {
    let d = crypto_lib::hash::sha3::shake256(label.as_bytes(), 16);
    u128::from_be_bytes(d.try_into().expect("sixteen bytes"))
}

/// `label` with `:i` appended for every `i > 0`, as B2's generator labels.
pub fn labelled(label: &str, i: u64) -> String {
    if i == 0 {
        label.to_string()
    } else {
        format!("{label}:{i}")
    }
}

// ── GF(p^k) ─────────────────────────────────────────────────────────

/// An element: its `k` coefficients in the basis `1, t, ..., t^(k-1)`.
pub type El = Vec<u128>;

/// `GF(p)[t] / (t^k + c_{k-1} t^{k-1} + ... + c_0)`.
pub struct Fpk {
    pub p: u128,
    pub k: usize,
    pub modc: Vec<u128>,
    pub q: u128,
    f: Poly,
    z: OnceCell<El>,
}

impl Fpk {
    pub fn new(p: u128, modc: &[u128]) -> Result<Self, String> {
        let k = modc.len();
        let q = p
            .checked_pow(k as u32)
            .filter(|q| *q < 1 << 126)
            .ok_or_else(|| format!("GF({p}^{k}) is past this generator's range"))?;
        let modc: Vec<u128> = modc.iter().map(|c| c % p).collect();
        let mut f = modc.clone();
        f.push(1);
        Ok(Fpk {
            p,
            k,
            modc,
            q,
            f,
            z: OnceCell::new(),
        })
    }

    pub fn zero(&self) -> El {
        vec![0; self.k]
    }

    pub fn one(&self) -> El {
        let mut e = self.zero();
        e[0] = 1 % self.p;
        e
    }

    /// `u`, of any length, reduced modulo the field's polynomial.
    pub fn reduce(&self, u: &[u128]) -> El {
        let (p, k) = (self.p, self.k);
        let mut u: Vec<u128> = u.iter().map(|x| x % p).collect();
        for i in (k..u.len()).rev() {
            let c = u[i];
            if c != 0 {
                u[i] = 0;
                for j in 0..k {
                    u[i - k + j] = (u[i - k + j] + p - c * self.modc[j] % p) % p;
                }
            }
        }
        u.resize(k.max(u.len()), 0);
        u.truncate(k);
        u
    }

    pub fn add(&self, u: &[u128], v: &[u128]) -> El {
        u.iter().zip(v).map(|(x, y)| (x + y) % self.p).collect()
    }

    pub fn sub(&self, u: &[u128], v: &[u128]) -> El {
        u.iter()
            .zip(v)
            .map(|(x, y)| (x + self.p - y) % self.p)
            .collect()
    }

    pub fn mul(&self, u: &[u128], v: &[u128]) -> El {
        let p = self.p;
        let mut prod = vec![0u128; u.len() + v.len() - 1];
        for (i, &x) in u.iter().enumerate() {
            if x != 0 {
                for (j, &y) in v.iter().enumerate() {
                    prod[i + j] = (prod[i + j] + x % p * (y % p)) % p;
                }
            }
        }
        self.reduce(&prod)
    }

    pub fn scalar(&self, c: u128) -> El {
        let mut e = self.zero();
        e[0] = c % self.p;
        e
    }

    pub fn pow(&self, u: &[u128], mut e: u128) -> El {
        let (mut out, mut base) = (self.one(), u.to_vec());
        while e > 0 {
            if e & 1 == 1 {
                out = self.mul(&out, &base);
            }
            base = self.mul(&base, &base);
            e >>= 1;
        }
        out
    }

    /// The inverse of a nonzero `u`, by the extended Euclidean algorithm in
    /// GF(p)[t].
    pub fn inv(&self, u: &[u128]) -> Result<El, String> {
        let p = self.p;
        let (mut r0, mut r1) = (self.f.clone(), trim(u.to_vec()));
        let (mut s0, mut s1): (Poly, Poly) = (Vec::new(), vec![1]);
        if r1.is_empty() {
            return Err("the inverse of zero".into());
        }
        while !r1.is_empty() {
            let (quo, rem) = pdivmod(&r0, &r1, p);
            (r0, r1) = (r1, rem);
            let next = psub(&s0, &pmul(&quo, &s1, p), p);
            (s0, s1) = (s1, next);
        }
        if r0.len() != 1 {
            return Err("the modulus is reducible".into());
        }
        let c = powm(r0[0], p - 2, p);
        Ok(self.reduce(&s0.iter().map(|x| x * c % p).collect::<Vec<_>>()))
    }

    pub fn is_square(&self, u: &[u128]) -> bool {
        u.iter().all(|&x| x == 0) || self.pow(u, (self.q - 1) / 2) == self.one()
    }

    /// The element whose coefficient `i` is SHAKE-256 of `label/i`, 16
    /// bytes, modulo `p`.
    pub fn element(&self, label: &str) -> El {
        (0..self.k)
            .map(|i| shake128bits(&format!("{label}/{i}")) % self.p)
            .collect()
    }

    /// The first labelled non-square, the same for every root taken.
    fn non_residue(&self) -> &El {
        self.z.get_or_init(|| {
            let mut i = 0;
            loop {
                let z = self.element(&labelled(&format!("{LABEL}/non-residue"), i));
                if !self.is_square(&z) {
                    return z;
                }
                i += 1;
            }
        })
    }

    /// A square root by Tonelli and Shanks, or `None`: the root the
    /// declaration's generator took, step for step.
    pub fn sqrt(&self, u: &[u128]) -> Option<El> {
        if u.iter().all(|&x| x == 0) {
            return Some(self.zero());
        }
        if !self.is_square(u) {
            return None;
        }
        let (mut m, mut s) = (self.q - 1, 0u32);
        while m.is_multiple_of(2) {
            m /= 2;
            s += 1;
        }
        let z = self.non_residue().clone();
        let one = self.one();
        let (mut c, mut t, mut r) = (self.pow(&z, m), self.pow(u, m), self.pow(u, m.div_ceil(2)));
        while t != one {
            let (mut j, mut t2) = (0u32, t.clone());
            while t2 != one {
                t2 = self.mul(&t2, &t2);
                j += 1;
            }
            let bb = self.pow(&c, 1u128 << (s - j - 1));
            let bb2 = self.mul(&bb, &bb);
            (s, c, t, r) = (j, bb2.clone(), self.mul(&t, &bb2), self.mul(&r, &bb));
        }
        Some(r)
    }

    /// `Σ e_i p^i`, an element as one integer.
    fn pack(&self, u: &[u128]) -> u128 {
        u.iter().rev().fold(0, |acc, &x| acc * self.p + x)
    }
}

// ── y^2 = x^3 + ax + b over GF(p^k) ──────────────────────────────────

/// A point, `None` the identity.
pub type Pt = Option<(El, El)>;

pub struct Curve<'f> {
    pub f: &'f Fpk,
    pub a: El,
    pub b: El,
}

impl Curve<'_> {
    pub fn singular(&self) -> bool {
        let f = self.f;
        let a3 = f.mul(&self.a, &f.mul(&self.a, &self.a));
        let lhs = f.add(
            &f.mul(&f.scalar(4), &a3),
            &f.mul(&f.scalar(27), &f.mul(&self.b, &self.b)),
        );
        lhs == f.zero()
    }

    pub fn rhs(&self, x: &[u128]) -> El {
        let f = self.f;
        f.add(&f.add(&f.mul(x, &f.mul(x, x)), &f.mul(&self.a, x)), &self.b)
    }

    pub fn on_curve(&self, pt: &Pt) -> bool {
        match pt {
            None => true,
            Some((x, y)) => self.f.mul(y, y) == self.rhs(x),
        }
    }

    pub fn neg(&self, pt: &Pt) -> Pt {
        pt.as_ref()
            .map(|(x, y)| (x.clone(), self.f.sub(&self.f.zero(), y)))
    }

    pub fn add(&self, pt: &Pt, qt: &Pt) -> Pt {
        let f = self.f;
        let ((x1, y1), (x2, y2)) = match (pt, qt) {
            (None, _) => return qt.clone(),
            (_, None) => return pt.clone(),
            (Some(a), Some(b)) => (a, b),
        };
        let lam = if x1 == x2 {
            if f.add(y1, y2) == f.zero() {
                return None;
            }
            let num = f.add(&f.mul(&f.scalar(3), &f.mul(x1, x1)), &self.a);
            f.mul(
                &num,
                &f.inv(&f.mul(&f.scalar(2), y1)).expect("2y is nonzero"),
            )
        } else {
            f.mul(
                &f.sub(y2, y1),
                &f.inv(&f.sub(x2, x1)).expect("x2 - x1 is nonzero"),
            )
        };
        let x3 = f.sub(&f.sub(&f.mul(&lam, &lam), x1), x2);
        let y3 = f.sub(&f.mul(&lam, &f.sub(x1, &x3)), y1);
        Some((x3, y3))
    }

    /// `[k]P`, left to right.
    pub fn mul(&self, k: u128, pt: &Pt) -> Pt {
        let mut r: Pt = None;
        for i in (0..128 - k.leading_zeros()).rev() {
            r = self.add(&r, &r);
            if k >> i & 1 == 1 {
                r = self.add(&r, pt);
            }
        }
        r
    }

    /// The first point whose abscissa is labelled `label`, `label:1`, ...
    pub fn point(&self, label: &str) -> (El, El) {
        let mut i = 0;
        loop {
            let x = self.f.element(&labelled(label, i));
            if let Some(y) = self.f.sqrt(&self.rhs(&x)) {
                return (x, y);
            }
            i += 1;
        }
    }

    fn key(&self, pt: &Pt) -> (u128, u128) {
        match pt {
            None => (u128::MAX, u128::MAX),
            Some((x, y)) => (self.f.pack(x), self.f.pack(y)),
        }
    }
}

// ── the group's order, from its points' orders ──────────────────────

/// The `m` in `[lo, hi]` with `[m]P = O`: a set, or every multiple of a
/// point's small order.
#[derive(Clone, Debug, PartialEq)]
enum Multiples {
    Of(u128),
    Set(BTreeSet<u128>),
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

impl Multiples {
    fn len(&self, lo: u128, hi: u128) -> u128 {
        match self {
            Multiples::Of(j) => hi / j - (lo - 1) / j,
            Multiples::Set(s) => s.len() as u128,
        }
    }

    fn first(&self, lo: u128) -> Option<u128> {
        match self {
            Multiples::Of(j) => Some(lo.div_ceil(*j) * j),
            Multiples::Set(s) => s.iter().next().copied(),
        }
    }

    fn and(self, other: Multiples) -> Multiples {
        match (self, other) {
            (Multiples::Of(a), Multiples::Of(b)) => Multiples::Of(a / gcd(a, b) * b),
            (Multiples::Of(j), Multiples::Set(s)) | (Multiples::Set(s), Multiples::Of(j)) => {
                Multiples::Set(s.into_iter().filter(|m| m.is_multiple_of(j)).collect())
            }
            (Multiples::Set(a), Multiples::Set(b)) => Multiples::Set(&a & &b),
        }
    }
}

fn isqrt(n: u128) -> u128 {
    if n < 2 {
        return n;
    }
    let mut x = (n as f64).sqrt() as u128;
    while x * x > n {
        x -= 1;
    }
    while (x + 1) * (x + 1) <= n {
        x += 1;
    }
    x
}

/// Every `m` in `[lo, hi]` with `[m]P = O`, by baby steps and giant steps:
/// B2's generator's `multiples_in`.
fn multiples_in(e: &Curve, pt: &Pt, lo: u128, hi: u128) -> Multiples {
    let width = hi - lo;
    let b = isqrt(width) + 1;
    let mut table: HashMap<(u128, u128), u128> = HashMap::new();
    let mut r: Pt = None;
    for j in 0..b {
        if r.is_none() && j > 0 {
            return Multiples::Of(j);
        }
        table.entry(e.key(&r)).or_insert(j);
        r = e.add(&r, pt);
    }
    let step = e.neg(&e.mul(b, pt));
    let mut q = e.neg(&e.mul(lo, pt));
    let mut hits = BTreeSet::new();
    for g in 0..width / b + 2 {
        if let Some(&j) = table.get(&e.key(&q)) {
            if g * b + j <= width {
                hits.insert(lo + g * b + j);
            }
        }
        q = e.add(&q, &step);
    }
    Multiples::Set(hits)
}

/// `#E` over a field of `q` elements: the Hasse interval narrowed by the
/// orders of the points `point(0)`, `point(1)`, ... (64 at most) until one
/// value is left, or `None` if it never is.
fn group_order(e: &Curve, q: u128, point: impl Fn(u64) -> Pt) -> Option<u128> {
    let s = isqrt(q);
    let (lo, hi) = (q + 1 - 2 * s - 2, q + 1 + 2 * s + 2);
    let mut cands: Option<Multiples> = None;
    for j in 0..64 {
        let hits = multiples_in(e, &point(j), lo, hi);
        let now = match cands.take() {
            None => hits,
            Some(c) => c.and(hits),
        };
        if now.len(lo, hi) == 1 {
            return now.first(lo);
        }
        cands = Some(now);
    }
    None
}

// ── integers: exact primality and factors ───────────────────────────

/// Below this, Miller–Rabin with the first thirteen primes as bases is
/// exact (`ψ_13`).
const PSI13: u128 = 3_317_044_064_679_887_385_961_981;

fn mulmod(a: u128, b: u128, m: u128) -> u128 {
    if m <= u64::MAX as u128 {
        return a % m * (b % m) % m;
    }
    let (mut a, mut b, mut r) = (a % m, b % m, 0u128);
    while b > 0 {
        if b & 1 == 1 {
            r = (r + a) % m;
        }
        a = (a + a) % m;
        b >>= 1;
    }
    r
}

fn powmod(mut b: u128, mut e: u128, m: u128) -> u128 {
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

/// B1's exact test: Miller–Rabin on thirteen bases, below `ψ_13`.
pub fn prime_exact(n: u128) -> Result<bool, String> {
    const BASES: [u128; 13] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41];
    if n >= PSI13 {
        return Err(format!("{n} is beyond the exact range"));
    }
    if n < 2 {
        return Ok(false);
    }
    for p in BASES {
        if n.is_multiple_of(p) {
            return Ok(n == p);
        }
    }
    let (mut d, mut s) = (n - 1, 0);
    while d.is_multiple_of(2) {
        d /= 2;
        s += 1;
    }
    'bases: for b in BASES {
        let mut x = powmod(b, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'bases;
            }
        }
        return Ok(false);
    }
    Ok(true)
}

/// A nontrivial factor of an odd composite, by Brent's variant of
/// Pollard's rho.
fn split(n: u128) -> u128 {
    for c in 1u128.. {
        let f = |x: u128| (mulmod(x, x, n) + c) % n;
        let (mut x, mut y, mut g) = (2u128, 2u128, 1u128);
        while g == 1 {
            x = f(x);
            y = f(f(y));
            g = gcd(x.abs_diff(y), n);
        }
        if g != n {
            return g;
        }
    }
    unreachable!("some constant splits a composite")
}

/// The prime factors of `n`, with multiplicity, in increasing order.
fn factor(n: u128) -> Result<Vec<u128>, String> {
    let (mut out, mut m, mut d) = (Vec::new(), n, 2u128);
    while d <= 1 << 16 && d * d <= m {
        while m.is_multiple_of(d) {
            out.push(d);
            m /= d;
        }
        d += if d == 2 { 1 } else { 2 };
    }
    let mut todo = if m > 1 { vec![m] } else { Vec::new() };
    while let Some(m) = todo.pop() {
        if prime_exact(m)? {
            out.push(m);
        } else {
            let d = if m.is_multiple_of(2) { 2 } else { split(m) };
            todo.push(d);
            todo.push(m / d);
        }
    }
    out.sort_unstable();
    Ok(out)
}

/// `(r, h)` by the instance's rule, or `None`.
fn subgroup(order: u128, q: u128, rule: &str) -> Result<Option<(u128, u128)>, String> {
    let floor = 4 * isqrt(q) + 4;
    if rule == "prime" {
        return Ok((prime_exact(order)? && order > floor).then_some((order, 1)));
    }
    let r = *factor(order)?.last().ok_or("an order of one")?;
    if !(floor < r && r < R_MAX) || order.is_multiple_of(r * r) {
        return Ok(None);
    }
    let h = order / r;
    if rule == "cofactor" && h == 1 {
        return Ok(None);
    }
    Ok(Some((r, h)))
}

// ── ICV1 for GF(p^k) (`docs/curves/ICV1.md`; `scripts/curve_id.py`) ─────

/// An ICV1 identity: the canonical string and the slug.
#[derive(Clone, Debug, PartialEq)]
pub struct Icv1 {
    pub icv1: String,
    pub slug: String,
}

fn strs(u: &[u128]) -> J {
    J::Arr(u.iter().map(|c| J::Str(c.to_string())).collect())
}

/// `fpk-<p>-<k>-<modhash8>`: an extension's field.
pub fn extension_field(p: u128, modc: &[u128]) -> String {
    let coeffs: Vec<String> = modc.iter().map(|c| (c % p).to_string()).collect();
    let hash = sha256_hex(format!("fpk-modulus:{p}:{}", coeffs.join(",")).as_bytes());
    format!("fpk-{p}-{}-{}", modc.len(), &hash[..8])
}

/// The identity of `y^2 = x^3 + ax + b` over `GF(p)[t]/(f)`, with `#E` the
/// whole group's order; the caller has checked that `f` is irreducible.
pub fn extension_id(
    p: u128,
    modc: &[u128],
    a: &[u128],
    b: &[u128],
    order: u128,
) -> Result<Icv1, String> {
    let k = modc.len();
    if p < 5 {
        return Err("an extension's characteristic is a prime above 3".into());
    }
    if k < 2 || a.len() != k || b.len() != k {
        return Err("an extension has degree k >= 2, and a and b k coefficients".into());
    }
    let f = Fpk::new(p, modc)?;
    let (a, b): (El, El) = (
        a.iter().map(|x| x % p).collect(),
        b.iter().map(|x| x % p).collect(),
    );
    let field = extension_field(p, &f.modc);
    let a3 = f.mul(&f.mul(&a, &a), &a);
    let bb = f.mul(&b, &b);
    let disc = f.reduce(
        &a3.iter()
            .zip(&bb)
            .map(|(x, y)| (4 * x + 27 * y) % p)
            .collect::<Vec<_>>(),
    );
    if disc.iter().all(|&x| x == 0) {
        return Err("singular curve".into());
    }
    let model_json = json::dumps_compact(
        &obj([
            ("a", strs(&a)),
            ("b", strs(&b)),
            ("field", J::Str(field.clone())),
            ("form", J::Str("y^2=x^3+a*x+b".into())),
            ("k", J::Str(k.to_string())),
            ("modulus", strs(&f.modc)),
            ("p", J::Str(p.to_string())),
            ("v", J::Str("1".into())),
        ]),
        true,
        true,
    );
    let trace = f.q as i128 + 1 - order as i128;
    if trace * trace > 4 * f.q as i128 {
        return Err(format!(
            "order {order} violates the Hasse bound over a field of {} elements",
            f.q
        ));
    }
    let scaled: El = a3.iter().map(|x| 6912 % p * x % p).collect();
    let j = f.mul(&scaled, &f.inv(&disc)?);
    let j: Vec<String> = j.iter().map(u128::to_string).collect();
    let model = sha256_hex(model_json.as_bytes());
    let bits = 128 - p.leading_zeros();
    let t = if trace < 0 {
        format!("tm{}", -trace)
    } else {
        format!("t{trace}")
    };
    Ok(Icv1 {
        icv1: format!(
            "ICV1:{field}:{trace}:{order}:{}:unk:unk:r:{}",
            j.join(","),
            &model[..12]
        ),
        slug: format!("icv1-fp{bits}k{k}-{t}-{}", &model[..8]),
    })
}

// ── the instances (design §4) ───────────────────────────────────────

/// id, `p`, `k`, the modulus rule, the subgroup rule, the purpose.
pub const SPECS: [(&str, u128, usize, &str, &str, &str); 8] = [
    (
        "G1",
        271,
        3,
        "binomial",
        "prime",
        "ic-gaudry-cubic at n ≈ 2^24 (§11.2's first prime)",
    ),
    (
        "G2",
        523,
        3,
        "binomial",
        "prime",
        "ic-gaudry-cubic at n ≈ 2^27",
    ),
    (
        "G3",
        1039,
        3,
        "binomial",
        "prime",
        "ic-gaudry-cubic at n ≈ 2^30",
    ),
    (
        "H1",
        271,
        3,
        "binomial",
        "cofactor",
        "a cofactor: refused by ic-gaudry-cubic",
    ),
    (
        "E2",
        (1 << 31) - 1,
        2,
        "labelled",
        "largest",
        "the one-word arithmetic at its widest, q ≈ 2^62",
    ),
    (
        "E5",
        2053,
        5,
        "labelled",
        "largest",
        "a general degree, q ≈ 2^55",
    ),
    (
        "E11",
        37,
        11,
        "labelled",
        "largest",
        "past eight coefficients, q ≈ 2^57",
    ),
    (
        "B2",
        34359738421,
        2,
        "labelled",
        "largest",
        "q ≈ 2^70, past one word: rho-bignum",
    ),
];

fn least_non_cube(p: u128) -> u128 {
    assert_eq!(p % 3, 1, "a binomial cubic needs p = 1 mod 3");
    (2..p)
        .find(|&c| powm(c, (p - 1) / 3, p) != 1)
        .expect("a non-cube")
}

/// The first labelled monic modulus of degree `k` that Rabin's test finds
/// irreducible, with its index.
fn labelled_modulus(label: &str, p: u128, k: usize) -> (Vec<u128>, u64) {
    let mut i = 0;
    loop {
        let lab = labelled(&format!("{label}/modulus"), i);
        let modc: Vec<u128> = (0..k)
            .map(|j| shake128bits(&format!("{lab}/{j}")) % p)
            .collect();
        if modc[0] != 0 && rabin_irreducible(p, &modc) {
            return (modc, i);
        }
        i += 1;
    }
}

/// `round(log2(x), 2)` as Python computes it.
fn log2_2(x: u128) -> f64 {
    format!("{:.2}", (x as f64).log2())
        .parse()
        .expect("a decimal")
}

/// One instance, searched as the declaration's generator searched it.
pub fn search(id: &str) -> Result<J, String> {
    let &(_, p, k, mod_rule, sub_rule, purpose) = SPECS
        .iter()
        .find(|s| s.0 == id)
        .ok_or_else(|| format!("no instance {id}"))?;
    if !prime_exact(p)? {
        return Err(format!("{p} is not prime"));
    }
    let label = format!("{LABEL}/{id}");
    let (modc, mod_index) = if mod_rule == "binomial" {
        let mut m = vec![0u128; k];
        m[0] = p - least_non_cube(p);
        (m, J::Null)
    } else {
        let (m, i) = labelled_modulus(&label, p, k);
        (m, J::Int(i.into()))
    };
    if !rabin_irreducible(p, &modc) {
        return Err(format!("{id}: the modulus is reducible"));
    }
    let f = Fpk::new(p, &modc)?;
    let mut i = 0u64;
    loop {
        let lab = labelled(&format!("{label}/curve"), i);
        let e = Curve {
            f: &f,
            a: f.element(&format!("{lab}/a")),
            b: f.element(&format!("{lab}/b")),
        };
        i += 1;
        if e.singular() {
            continue;
        }
        let Some(order) = group_order(&e, f.q, |j| Some(e.point(&format!("{lab}/order/{j}"))))
        else {
            continue;
        };
        let Some((r, h)) = subgroup(order, f.q, sub_rule)? else {
            continue;
        };
        let ident = extension_id(p, &modc, &e.a, &e.b, order)?;
        return Ok(obj([
            ("id", J::Str(id.into())),
            ("purpose", J::Str(purpose.into())),
            ("p", J::Str(p.to_string())),
            ("k", J::Int(k as i128)),
            ("modulus", strs(&modc)),
            ("modulus_rule", J::Str(mod_rule.into())),
            ("modulus_index", mod_index),
            ("curve_index", J::Int((i - 1).into())),
            ("a", strs(&e.a)),
            ("b", strs(&e.b)),
            ("order", J::Str(order.to_string())),
            ("r", J::Str(r.to_string())),
            ("h", J::Str(h.to_string())),
            ("log2_q", J::Float(log2_2(f.q))),
            ("log2_r", J::Float(log2_2(r))),
            ("slug", J::Str(ident.slug)),
            ("icv1", J::Str(ident.icv1)),
        ]));
    }
}

/// `instances.json`'s text, from the records in `SPECS`'s order.
pub fn instances_text(records: Vec<J>) -> String {
    json::dumps_utf8(
        &obj([
            ("label", J::Str(LABEL.into())),
            ("generator", J::Str(GENERATOR.into())),
            ("instances", J::Arr(records)),
        ]),
        1,
    ) + "\n"
}

pub fn instances(ids: &[String]) -> Result<String, String> {
    let ids: Vec<&str> = if ids.is_empty() {
        SPECS.iter().map(|s| s.0).collect()
    } else {
        ids.iter().map(String::as_str).collect()
    };
    let mut records = Vec::new();
    for id in ids {
        let rec = search(id)?;
        eprintln!("{id}: {}", json::dumps_line(&rec, false));
        records.push(rec);
    }
    Ok(instances_text(records))
}

/// `instances.json` written, refusing to overwrite it; with `check`, every
/// instance searched again and the file's bytes compared.
pub fn write_instances(root: &Path, check: bool) -> Result<(J, bool), String> {
    let path = root.join(INSTANCES);
    if check {
        let want = std::fs::read_to_string(&path).map_err(|e| format!("{INSTANCES}: {e}"))?;
        let ok = instances(&[])? == want;
        return Ok((
            obj([
                ("instances", J::Int(SPECS.len() as i128)),
                ("same", J::Bool(ok)),
            ]),
            ok,
        ));
    }
    if path.exists() {
        return Err("instances.json exists; B5a's instances are frozen (use --check)".into());
    }
    let text = instances(&[])?;
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    std::fs::write(&path, text).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok((obj([("instances", J::Int(SPECS.len() as i128))]), true))
}

/// The big-endian integer of `label`'s SHA-256, modulo `r - 1`, plus one.
pub fn known_log(label: &str, r: u128) -> u128 {
    let d = BigUint::from_bytes_be(&crypto_lib::hash::sha256::sha256(label.as_bytes()));
    (d % BigUint::from(r - 1)).to_u128().expect("below r") + 1
}

// ── measurement 5's replay ──────────────────────────────────────────

fn coeffs(v: &J) -> Result<El, String> {
    v.as_arr()
        .ok_or("a coefficient list is not a list")?
        .iter()
        .map(|c| {
            c.as_str()
                .ok_or("a coefficient is not a string")?
                .parse::<u128>()
                .map_err(|e| e.to_string())
        })
        .collect()
}

fn doc_point(v: &J) -> Result<(El, El), String> {
    Ok((coeffs(v.at("x")?)?, coeffs(v.at("y")?)?))
}

/// `[scalar]G = Q` on a `prime_extension` document's short Weierstrass
/// curve, in this module's arithmetic, and the known logarithm: the same
/// record as `f0::replay` writes for a binary document.
pub fn replay(doc: &J, scalar: Option<&BigUint>) -> Result<J, String> {
    let field = doc.at("field")?;
    let p: u128 = field
        .at("p")?
        .as_str()
        .ok_or("`field.p`")?
        .parse()
        .map_err(|e| format!("`field.p`: {e}"))?;
    let f = Fpk::new(p, &coeffs(field.at("modulus")?)?)?;
    let curve = doc.at("curve")?;
    if curve.get("form").and_then(J::as_str) != Some("short_weierstrass") {
        return Err("the replay reads short Weierstrass documents".into());
    }
    let e = Curve {
        f: &f,
        a: f.reduce(&coeffs(curve.at("a")?)?),
        b: f.reduce(&coeffs(curve.at("b")?)?),
    };
    let sub = doc.at("subgroup")?;
    let g = Some(doc_point(sub.at("generator")?)?);
    let order: u128 = sub
        .at("order")?
        .as_str()
        .ok_or("`subgroup.order`")?
        .parse()
        .map_err(|e| format!("`subgroup.order`: {e}"))?;
    let target = doc.at("target")?;
    let known: Option<u128> = match target.get("known_log") {
        Some(k) => Some(
            k.as_str()
                .ok_or("`target.known_log`")?
                .parse()
                .map_err(|e| format!("`target.known_log`: {e}"))?,
        ),
        None => None,
    };
    let q = match known {
        Some(d) => e.mul(d, &g),
        None => Some(doc_point(target.at("point")?)?),
    };
    let scalar = scalar
        .map(|s| s.to_u128().ok_or("a scalar past 128 bits"))
        .transpose()?;
    // A point off the curve replays nothing: the document is not what it says.
    let on_curve = e.on_curve(&g) && e.on_curve(&q);
    let replays = on_curve && scalar.is_some_and(|d| e.mul(d, &g) == q);
    let known_matches = match (known, scalar) {
        (None, _) => J::Null,
        (Some(d), Some(s)) => J::Bool(s == d % order),
        (Some(_), None) => J::Bool(false),
    };
    Ok(obj([
        ("replays", J::Bool(replays)),
        ("points_on_curve", J::Bool(on_curve)),
        (
            "known_log",
            known.map_or(J::Null, |d| J::Str(d.to_string())),
        ),
        ("known_matches", known_matches),
    ]))
}

// ── measurement 6: the estimates' constants ─────────────────────────

/// The calibration's instances: `rho-negation` on E2, E5 and E11,
/// `rho-bignum` on B2, and `ic-gaudry-cubic` paired with `rho-negation` on
/// G1–G3, each on its known-answer document.
pub const CALIBRATE: [&str; 7] = ["E2", "E5", "E11", "B2", "G1", "G2", "G3"];
/// Runs of each, round by round.
pub const CALIBRATION_RUNS: u32 = 5;
/// Measurement 5's hour, for each run.
const HOUR_S: u64 = 3600;

fn calibration_out(runs_dir: &Path, id: &str, k: u32) -> PathBuf {
    runs_dir
        .join("calibrate")
        .join(id)
        .join(format!("r{k}.price.json"))
}

/// Five isolated runs of each instance, every instance once a round; a
/// contended or failed attempt is kept and the run made again, at most
/// twice, and the first clean, complete attempt is the run.
pub fn calibrate(root: &Path, b: &Bench, ic: &Path, runs_dir: &Path) -> Result<(), String> {
    let params = root.join(CASES_DIR).join("params");
    let tree = runs_dir.join("calibrate");
    for k in 1..=CALIBRATION_RUNS {
        for id in CALIBRATE {
            let doc = params.join(format!("{id}-known.json"));
            let paired = json::read(&doc)?.at("method")?.at("solve")?.as_str() == Some("paired");
            let out = calibration_out(runs_dir, id, k);
            let mut rep = J::Null;
            for a in 0..=runs::RETRIES {
                let attempt = runs::attempt(&out, a);
                if !attempt.exists() {
                    if let Some(parent) = attempt.parent() {
                        std::fs::create_dir_all(parent)
                            .map_err(|e| format!("{}: {e}", parent.display()))?;
                    }
                    let mut cmd: Vec<String> = vec![
                        "timeout".into(),
                        HOUR_S.to_string(),
                        ic.to_string_lossy().into_owned(),
                        "price".into(),
                        "--params".into(),
                        doc.to_string_lossy().into_owned(),
                        "--json".into(),
                        "--out".into(),
                        attempt.to_string_lossy().into_owned(),
                    ];
                    if paired {
                        cmd.extend(["--repeats", "1", "--repeats-fast", "1"].map(String::from));
                    }
                    b.launch(&cmd, &attempt, &tree)?;
                }
                rep = runs::load(&attempt);
                if runs::clean(&attempt)?
                    && rep.get("status").and_then(J::as_str) == Some("complete")
                {
                    break;
                }
            }
            println!(
                "r{k} {id}: {}",
                rep.get("status").and_then(J::as_str).unwrap_or("None")
            );
        }
    }
    Ok(())
}

fn num_at(v: &J, path: &[&str]) -> Option<f64> {
    let mut cur = v;
    for k in path {
        cur = cur.get(k)?;
    }
    cur.as_f64()
}

/// A run's figure: its first clean, complete attempt, if any.
fn calibration_figure(runs_dir: &Path, id: &str, k: u32) -> Result<Option<J>, String> {
    let out = calibration_out(runs_dir, id, k);
    for a in 0..=runs::RETRIES {
        let attempt = runs::attempt(&out, a);
        if attempt.exists() && runs::clean(&attempt)? {
            let rep = runs::load(&attempt);
            if rep.get("status").and_then(J::as_str) == Some("complete") {
                return Ok(Some(rep));
            }
        }
    }
    Ok(None)
}

fn median_or_null(xs: &[f64]) -> J {
    if xs.is_empty() {
        J::Null
    } else {
        J::Float(stats::median(xs))
    }
}

/// The constants measurement 6 sets, from its runs: each instance's
/// medians over its runs, then
/// - `rho-negation`'s step cost: the median of E2's, E5's and E11's;
/// - `rho-bignum`'s: B2's, at the width of B2's field;
/// - `ic-gaudry-cubic`'s exponent: the least-squares slope of the log of
///   its median operations on the log of `r` over G1–G3, and its anchor
///   G3's medians, the largest rung's.
///
/// They are reported here; B5a's results pull request writes them into
/// `estimates.json`, whose provisional constants they replace.
pub fn estimates(root: &Path, runs_dir: &Path) -> Result<J, String> {
    let inst = json::read(&root.join(INSTANCES))?;
    let recs: Vec<Rec> = inst
        .at("instances")?
        .as_arr()
        .ok_or("`instances`")?
        .iter()
        .map(Rec::read)
        .collect::<Result<_, _>>()?;
    let rec = |id: &str| recs.iter().find(|r| r.id == id).ok_or(format!("no {id}"));
    let mut per = Vec::new();
    let mut rho_ns: HashMap<&str, J> = HashMap::new();
    let mut gaudry: Vec<(f64, f64, f64, f64)> = Vec::new(); // log2 r, ops, wall s, ns/op
    for id in CALIBRATE {
        let (mut step_ns, mut ops, mut wall_s) = (Vec::new(), Vec::new(), Vec::new());
        let mut missing = 0;
        for k in 1..=CALIBRATION_RUNS {
            let Some(rep) = calibration_figure(runs_dir, id, k)? else {
                missing += 1;
                continue;
            };
            if let (Some(ns), Some(steps)) = (
                num_at(&rep, &["rho", "online_ns"]),
                num_at(&rep, &["rho_counts", "steps"]),
            ) {
                step_ns.push(ns / steps);
            }
            if let (Some(o), Some(set), Some(on)) = (
                num_at(&rep, &["counts", "total_ops"]),
                num_at(&rep, &["ic", "setup_ns"]),
                num_at(&rep, &["ic", "online_ns"]),
            ) {
                ops.push(o);
                wall_s.push((set + on) / 1e9);
            }
        }
        let r = rec(id)?;
        let row = obj([
            ("id", s(id)),
            ("slug", s(&r.slug)),
            ("runs", J::Int((step_ns.len().max(ops.len())) as i128)),
            ("missing", J::Int(missing)),
            ("rho_ns_per_step", median_or_null(&step_ns)),
            ("ic_total_ops", median_or_null(&ops)),
            ("ic_wall_s", median_or_null(&wall_s)),
        ]);
        rho_ns.insert(id, row.get("rho_ns_per_step").cloned().unwrap_or(J::Null));
        if !ops.is_empty() {
            let (o, w) = (stats::median(&ops), stats::median(&wall_s));
            gaudry.push(((r.r as f64).log2(), o, w, w * 1e9 / o));
        }
        per.push(row);
    }
    let negation: Vec<f64> = ["E2", "E5", "E11"]
        .iter()
        .filter_map(|id| rho_ns.get(id).and_then(J::as_f64))
        .collect();
    let b2 = rec("B2")?;
    let width = 128 - b2.field()?.q.leading_zeros();
    let fit = stats::fit(
        &gaudry.iter().map(|g| g.0).collect::<Vec<_>>(),
        &gaudry.iter().map(|g| g.1.log2()).collect::<Vec<_>>(),
        0.95,
    );
    let g3 = rec("G3")?;
    let anchor = gaudry
        .last()
        .filter(|_| gaudry.len() == 3)
        .map(|&(lr, o, w, nso)| {
            obj([
                ("p", J::Int(g3.p as i128)),
                ("log2_r", J::Float(lr)),
                ("total_ops", J::Float(o)),
                ("wall_s", J::Float(w)),
                ("ns_per_op", J::Float(nso)),
            ])
        });
    Ok(obj([
        (
            "what",
            s(
                "B5a's measurement 6: the constants its F2 estimates read for extension fields, \
               from five isolated runs of each instance; B5a's results pull request writes them \
               into estimates.json",
            ),
        ),
        ("instances", J::Arr(per)),
        (
            "rho-negation",
            obj([(
                "extension",
                obj([
                    ("ns_per_step", median_or_null(&negation)),
                    ("from", s("the median of E2's, E5's and E11's medians")),
                ]),
            )]),
        ),
        (
            "rho-bignum",
            obj([(
                "extension",
                obj([(
                    "points",
                    J::Arr(vec![obj([
                        ("width", J::Int(width.into())),
                        ("ns_per_step", rho_ns.get("B2").cloned().unwrap_or(J::Null)),
                    ])]),
                )]),
            )]),
        ),
        (
            "ic-gaudry-cubic",
            obj([
                (
                    "alpha",
                    fit.as_ref()
                        .and_then(|f| f.get("beta").cloned())
                        .unwrap_or(J::Null),
                ),
                ("fit", fit.unwrap_or(J::Null)),
                ("anchor", anchor.unwrap_or(J::Null)),
            ]),
        ),
    ]))
}

// ── the conformance cases (`conformance/v2-b5a/`) ───────────────────

const STEP: &str = "B5a";
/// The cases' label: generators, known logarithms and public points.
pub const CASES_LABEL: &str = "ic-conformance-v2-b5a";
/// The suite's rho seed for target `T01`, as B1 to B4 use.
const RHO_SEED: i128 = 0x230000 + 1;
const TIMEOUT: i128 = 1800;
/// B2's document, copied byte for byte (C103, C107).
const C050_DOC: &str = "C050-cubic-extension.json";
const C050_KNOWN_LOG: &str = "554837";
const CASES_DIR: &str = "research/ic_tool_program/conformance/v2-b5a";
const INSTANCES: &str = "research/ic_tool_program/rounds/B5a-extension-fields/instances.json";

fn s(v: &str) -> J {
    J::Str(v.into())
}

fn rho_alone(seed: i128) -> J {
    obj([
        ("solve", s("rho")),
        (
            "rho",
            obj([("pipeline", s("auto")), ("seed", J::Int(seed))]),
        ),
    ])
}

fn paired_auto(seed: i128) -> J {
    obj([
        ("solve", s("paired")),
        (
            "index_calculus",
            obj([("pipeline", s("auto")), ("recipe", s("auto"))]),
        ),
        (
            "rho",
            obj([("pipeline", s("auto")), ("seed", J::Int(seed))]),
        ),
    ])
}

/// An instance's record, as numbers.
struct Rec {
    id: String,
    p: u128,
    modc: Vec<u128>,
    a: El,
    b: El,
    order: u128,
    r: u128,
    h: u128,
    slug: String,
    icv1: String,
    log2_r: J,
    json: J,
}

fn big(v: &J, key: &str) -> Result<u128, String> {
    v.at(key)?
        .as_str()
        .ok_or_else(|| format!("`{key}` is not a string"))?
        .parse()
        .map_err(|e| format!("`{key}`: {e}"))
}

fn bigs(v: &J, key: &str) -> Result<Vec<u128>, String> {
    v.at(key)?
        .as_arr()
        .ok_or_else(|| format!("`{key}` is not a list"))?
        .iter()
        .map(|c| {
            c.as_str()
                .ok_or("a coefficient is not a string")?
                .parse()
                .map_err(|e| format!("`{key}`: {e}"))
        })
        .collect()
}

impl Rec {
    fn read(v: &J) -> Result<Rec, String> {
        let text = |k: &str| -> Result<String, String> {
            Ok(v.at(k)?
                .as_str()
                .ok_or_else(|| format!("`{k}`"))?
                .to_string())
        };
        Ok(Rec {
            id: text("id")?,
            p: big(v, "p")?,
            modc: bigs(v, "modulus")?,
            a: bigs(v, "a")?,
            b: bigs(v, "b")?,
            order: big(v, "order")?,
            r: big(v, "r")?,
            h: big(v, "h")?,
            slug: text("slug")?,
            icv1: text("icv1")?,
            log2_r: v.at("log2_r")?.clone(),
            json: v.clone(),
        })
    }

    fn field(&self) -> Result<Fpk, String> {
        Fpk::new(self.p, &self.modc)
    }
}

/// `[h]P` for the first labelled `P` with `[h]P` not the identity, checked
/// to lie on the curve with order `r`.
fn subgroup_point(e: &Curve, label: &str, h: u128, r: u128) -> Result<(El, El), String> {
    let mut j = 0;
    loop {
        let g = e.mul(h, &Some(e.point(&labelled(label, j))));
        if let Some(pt) = g.clone() {
            if !e.on_curve(&g) || e.mul(r, &g).is_some() {
                return Err(format!("{label}: [h]P does not have order r"));
            }
            return Ok(pt);
        }
        j += 1;
    }
}

/// An instance re-checked (design §4.5's method 3): the modulus
/// irreducible, the curve non-singular, `r` prime with `r > 4√q` (so `#E`
/// is the only multiple of `r` in the Hasse interval) and `r²` not
/// dividing `#E`, the generator of order `r`, and the slug the reference's.
fn checked(rec: &Rec, f: &Fpk) -> Result<(El, El), String> {
    let id = &rec.id;
    if !rabin_irreducible(rec.p, &rec.modc) {
        return Err(format!("{id}: the modulus is reducible"));
    }
    let e = Curve {
        f,
        a: rec.a.clone(),
        b: rec.b.clone(),
    };
    if e.singular() {
        return Err(format!("{id}: singular"));
    }
    if rec.order != rec.h * rec.r || !prime_exact(rec.r)? || rec.order.is_multiple_of(rec.r * rec.r)
    {
        return Err(format!(
            "{id}: the order is not h·r with r prime and r² not dividing it"
        ));
    }
    let g = subgroup_point(&e, &format!("{CASES_LABEL}/{id}/generator"), rec.h, rec.r)?;
    let q = f.q;
    let sq = isqrt(q);
    let (lo, hi) = (q + 1 - 2 * sq - 2, q + 1 + 2 * sq + 2);
    let multiples: Vec<u128> = (lo - lo % rec.r..=hi)
        .step_by(usize::try_from(rec.r).map_err(|_| "r is past usize")?)
        .filter(|&m| m >= lo)
        .collect();
    if rec.r <= 4 * sq + 4 || multiples != [rec.order] {
        return Err(format!(
            "{id}: #E is not the one multiple of r in the Hasse interval"
        ));
    }
    let ident = extension_id(rec.p, &rec.modc, &rec.a, &rec.b, rec.order)?;
    if ident.slug != rec.slug || ident.icv1 != rec.icv1 {
        return Err(format!("{id}: the identity is not the record's"));
    }
    Ok(g)
}

fn point_doc(pt: &(El, El)) -> J {
    obj([("x", strs(&pt.0)), ("y", strs(&pt.1))])
}

fn doc(name: &str, rec: &Rec, g: &(El, El), target: J, method: J, curve: Option<J>) -> J {
    assert!(name.chars().count() <= 120, "{name}");
    let j = &rec.json;
    obj([
        ("schema_version", J::Int(2)),
        ("name", s(name)),
        (
            "field",
            obj([
                ("kind", s("prime_extension")),
                ("p", j.get("p").cloned().unwrap_or(J::Null)),
                ("degree", j.get("k").cloned().unwrap_or(J::Null)),
                ("modulus", j.get("modulus").cloned().unwrap_or(J::Null)),
            ]),
        ),
        (
            "curve",
            curve.unwrap_or_else(|| {
                obj([
                    ("form", s("short_weierstrass")),
                    ("a", j.get("a").cloned().unwrap_or(J::Null)),
                    ("b", j.get("b").cloned().unwrap_or(J::Null)),
                ])
            }),
        ),
        (
            "subgroup",
            obj([
                ("order", j.get("r").cloned().unwrap_or(J::Null)),
                ("cofactor", j.get("h").cloned().unwrap_or(J::Null)),
                ("generator", point_doc(g)),
            ]),
        ),
        ("target", target),
        ("method", method),
    ])
}

fn strs_of(v: &[&str]) -> J {
    J::Arr(v.iter().map(|x| s(x)).collect())
}

fn kv(pairs: Vec<(&str, J)>) -> J {
    J::Obj(pairs.into_iter().map(|(k, v)| (k.to_string(), v)).collect())
}

#[allow(clippy::too_many_arguments)]
fn case(
    cid: &str,
    purpose: &str,
    files: J,
    argv: J,
    expect: J,
    timeout: i128,
    extra: Vec<(&str, J)>,
) -> J {
    let mut out = vec![
        ("id", s(cid)),
        ("step", s(STEP)),
        ("checks", s(purpose)),
        ("files", files),
        ("argv", argv),
        ("env", obj([("RAYON_NUM_THREADS", s("1"))])),
        ("expect", expect),
        ("timeout_s", J::Int(timeout)),
    ];
    out.extend(extra);
    kv(out)
}

fn price_argv() -> J {
    strs_of(&[
        "price",
        "--params",
        "{tmp}/params.json",
        "--json",
        "--out",
        "{tmp}/report.json",
        "--repeats",
        "1",
        "--repeats-fast",
        "1",
    ])
}

fn plain_price_argv() -> J {
    strs_of(&[
        "price",
        "--params",
        "{tmp}/params.json",
        "--json",
        "--out",
        "{tmp}/report.json",
    ])
}

fn check_argv() -> J {
    strs_of(&[
        "check",
        "--params",
        "{tmp}/params.json",
        "--json",
        "--out",
        "{tmp}/report.json",
    ])
}

fn copy(path: &str) -> J {
    obj([("params.json", obj([("copy", s(path))]))])
}

fn copy_set(path: &str, set: J) -> J {
    obj([("params.json", obj([("copy", s(path)), ("set", set)]))])
}

fn considered(rows: &[(&str, bool, Option<&str>)]) -> J {
    J::Arr(
        rows.iter()
            .map(|(p, a, g)| {
                let mut kv = vec![
                    ("pipeline".to_string(), s(p)),
                    ("admitted".to_string(), J::Bool(*a)),
                ];
                if let Some(g) = g {
                    kv.push(("gate".to_string(), s(g)));
                }
                J::Obj(kv)
            })
            .collect(),
    )
}

/// `complete` followed by `rest`, in that order.
fn complete_and(rest: Vec<(&str, J)>) -> J {
    let mut out = vec![
        ("status", s("complete")),
        ("result.verified", J::Bool(true)),
    ];
    out.extend(rest);
    kv(out)
}

/// `E` written as `y² + a1xy + a3y = x³ + a2x² + a4x + a6` by
/// `x = X + s`, `y = Y + λX + μ`, with `s`, `λ` and `μ` labelled, and the
/// map of points to it.
fn general_form<'f>(
    f: &'f Fpk,
    e: &Curve,
    label: &str,
) -> (J, impl Fn(&(El, El)) -> (El, El) + 'f) {
    let [sv, lam, mu] = ["s", "lambda", "mu"].map(|v| f.element(&format!("{label}/{v}")));
    let (two, three) = (f.scalar(2), f.scalar(3));
    let a1 = f.mul(&two, &lam);
    let a3 = f.mul(&two, &mu);
    let a2 = f.sub(&f.mul(&three, &sv), &f.mul(&lam, &lam));
    let a4 = f.sub(
        &f.add(&f.mul(&three, &f.mul(&sv, &sv)), &e.a),
        &f.mul(&two, &f.mul(&lam, &mu)),
    );
    let a6 = f.sub(
        &f.add(
            &f.add(&f.mul(&sv, &f.mul(&sv, &sv)), &f.mul(&e.a, &sv)),
            &e.b,
        ),
        &f.mul(&mu, &mu),
    );
    let curve = obj([
        ("form", s("general_weierstrass")),
        ("a1", strs(&a1)),
        ("a2", strs(&a2)),
        ("a3", strs(&a3)),
        ("a4", strs(&a4)),
        ("a6", strs(&a6)),
    ]);
    let to_general = move |pt: &(El, El)| {
        let x = f.sub(&pt.0, &sv);
        let y = f.sub(&f.sub(&pt.1, &f.mul(&lam, &x)), &mu);
        let lhs = f.add(
            &f.add(&f.mul(&y, &y), &f.mul(&a1, &f.mul(&x, &y))),
            &f.mul(&a3, &y),
        );
        let rhs = f.add(
            &f.add(
                &f.add(&f.mul(&x, &f.mul(&x, &x)), &f.mul(&a2, &f.mul(&x, &x))),
                &f.mul(&a4, &x),
            ),
            &a6,
        );
        assert_eq!(lhs, rhs, "the general form's equation");
        (x, y)
    };
    (curve, to_general)
}

/// What the cases' files hold: the parameter documents (by name, in the
/// order written) and the cases, sorted by id.
fn build(root: &Path) -> Result<(Vec<(String, String)>, Vec<J>), String> {
    let inst = json::read(&root.join(INSTANCES))?;
    let recs: Vec<Rec> = inst
        .at("instances")?
        .as_arr()
        .ok_or("`instances`")?
        .iter()
        .map(Rec::read)
        .collect::<Result<_, _>>()?;
    let rec = |id: &str| recs.iter().find(|r| r.id == id).ok_or(format!("no {id}"));
    let mut files: Vec<(String, String)> = Vec::new();
    let mut cases: Vec<J> = Vec::new();
    let put = |files: &mut Vec<(String, String)>, name: String, d: &J| {
        files.push((name, json::dumps(d, 1) + "\n"));
    };
    let (paired, rho) = (paired_auto(RHO_SEED), rho_alone(RHO_SEED));
    let mut known: HashMap<String, u128> = HashMap::new();

    // Every instance: its known-answer document, and measurement 5's public point.
    let fields: Vec<Fpk> = recs.iter().map(Rec::field).collect::<Result<_, _>>()?;
    let mut gens: HashMap<String, (El, El)> = HashMap::new();
    for (r, f) in recs.iter().zip(&fields) {
        let g = checked(r, f)?;
        let e = Curve {
            f,
            a: r.a.clone(),
            b: r.b.clone(),
        };
        let k = known_log(&format!("{CASES_LABEL}/{}/known_log", r.id), r.r);
        let method = if r.id.starts_with('G') {
            paired.clone()
        } else {
            rho.clone()
        };
        put(
            &mut files,
            format!("{}-known.json", r.id),
            &doc(
                &format!("B5a {}: {}, a known logarithm's multiple", r.id, r.slug),
                r,
                &g,
                obj([("known_log", s(&k.to_string()))]),
                method.clone(),
                None,
            ),
        );
        let t = subgroup_point(&e, &format!("{CASES_LABEL}/{}/T001", r.id), r.h, r.r)?;
        put(
            &mut files,
            format!("{}-T001.json", r.id),
            &doc(
                &format!(
                    "B5a {}-T001: {}, a public point for measurement 5",
                    r.id, r.slug
                ),
                r,
                &g,
                obj([("point", point_doc(&t))]),
                method,
                None,
            ),
        );
        known.insert(r.id.clone(), k);
        gens.insert(r.id.clone(), g);
    }
    let gaudry_rho = considered(&[
        ("ic-gaudry-cubic", true, None),
        ("rho-negation", true, None),
    ]);

    // C103: C050's successor; C050 names B5 as its `until`, and B5 is now two steps.
    let c050 = std::fs::read_to_string(
        root.join("research/ic_tool_program/conformance/v2-b2/params")
            .join(C050_DOC),
    )
    .map_err(|e| format!("{C050_DOC}: {e}"))?;
    files.push((C050_DOC.into(), c050));
    cases.push(case(
        "C103-cubic-extension-rho",
        "C050's successor from B5a: its GF(1009^3) document, with a general modulus and h = 1524, \
         recovers its known logarithm with rho-negation",
        copy(&format!("{{here}}/{C050_DOC}")),
        plain_price_argv(),
        obj([
            ("exit", J::Int(0)),
            ("json_file", s("{tmp}/report.json")),
            (
                "json_paths",
                complete_and(vec![
                    ("route.rho.pipeline", s("rho-negation")),
                    ("result.scalar", s(C050_KNOWN_LOG)),
                    ("result.known_answer", J::Bool(true)),
                ]),
            ),
        ]),
        300,
        vec![("supersedes", s("C050-cubic-extension"))],
    ));

    // C104-C106: Gaudry's index calculus at F0, paired with rho-negation.
    for (cid, id) in [("C104", "G1"), ("C105", "G2"), ("C106", "G3")] {
        let r = rec(id)?;
        let log2_r = json::dumps_line(&r.log2_r, false);
        cases.push(case(
            &format!("{cid}-gaudry-cubic-{}", id.to_lowercase()),
            &format!(
                "ic-gaudry-cubic paired with rho-negation at F0 on {id}, GF({}^3) with t^3 - c and a \
                 prime group order of 2^{log2_r}",
                r.p
            ),
            copy(&format!("{{here}}/{id}-known.json")),
            price_argv(),
            obj([
                ("exit", J::Int(0)),
                ("json_file", s("{tmp}/report.json")),
                (
                    "json_paths",
                    complete_and(vec![
                        ("fidelity", s("F0")),
                        ("result.scalar", s(&known[id].to_string())),
                        ("result.known_answer", J::Bool(true)),
                        ("all_verified", J::Bool(true)),
                        ("ic_and_rho_agree", J::Bool(true)),
                        ("route.ic.pipeline", s("ic-gaudry-cubic")),
                        ("route.rho.pipeline", s("rho-negation")),
                        ("speedup_eligible", J::Bool(false)),
                        ("curve_id.slug", s(&r.slug)),
                    ]),
                ),
                (
                    "json_contains",
                    obj([
                        ("route.considered", gaudry_rho.clone()),
                        (
                            "disclosures",
                            J::Arr(vec![obj([("code", s("study-pipeline"))])]),
                        ),
                    ]),
                ),
            ]),
            TIMEOUT,
            Vec::new(),
        ));
    }

    // C107-C109: the index calculus's gates, each with rho admitted.
    for (cid, purpose, files_, gate) in [
        (
            "C107-cubic-modulus-not-binomial",
            "C050's document paired: its general cubic modulus is refused by ic-gaudry-cubic as \
             modulus-not-binomial, and rho-negation admits it",
            copy_set(
                &format!("{{here}}/{C050_DOC}"),
                obj([("method", paired_auto(7))]),
            ),
            "modulus-not-binomial",
        ),
        (
            "C108-gaudry-cofactor-not-one",
            "H1 paired: a cofactor of 395 is refused by ic-gaudry-cubic as cofactor-not-one, and \
             rho-negation admits it",
            copy_set("{here}/H1-known.json", obj([("method", paired.clone())])),
            "cofactor-not-one",
        ),
        (
            "C109-quadratic-extension-paired",
            "E2 paired: GF(p^2) is refused by ic-gaudry-cubic as extension-degree-not-three, and \
             rho-negation admits it",
            copy_set("{here}/E2-known.json", obj([("method", paired.clone())])),
            "extension-degree-not-three",
        ),
    ] {
        cases.push(case(
            cid,
            purpose,
            files_,
            price_argv(),
            obj([
                ("exit", J::Int(3)),
                ("json_file", s("{tmp}/report.json")),
                (
                    "json_paths",
                    obj([
                        ("refusal.code", s("no-ic-route")),
                        ("refusal.class", s("unsupported")),
                    ]),
                ),
                (
                    "json_contains",
                    obj([(
                        "route.considered",
                        considered(&[
                            ("ic-gaudry-cubic", false, Some(gate)),
                            ("rho-negation", true, None),
                        ]),
                    )]),
                ),
            ]),
            300,
            Vec::new(),
        ));
    }

    // C110-C113: rho alone on every width and degree.
    type Extra<'a> = &'a [(&'a str, bool, Option<&'a str>)];
    let b2_extra: Extra = &[
        ("rho-negation", false, Some("field-wider-than-one-word")),
        ("rho-bignum", true, None),
    ];
    for (cid, id, purpose, pipe, extra) in [
        (
            "C110",
            "E2",
            "the one-word arithmetic at its widest, q = (2^31 - 1)^2",
            "rho-negation",
            &[][..],
        ),
        (
            "C111",
            "E5",
            "a degree-5 field, q ≈ 2^55",
            "rho-negation",
            &[][..],
        ),
        (
            "C112",
            "E11",
            "a degree-11 field, past eight coefficients",
            "rho-negation",
            &[][..],
        ),
        (
            "C113",
            "B2",
            "q ≈ 2^70, past one word: rho-negation is refused and rho-bignum runs",
            "rho-bignum",
            b2_extra,
        ),
    ] {
        let r = rec(id)?;
        let mut expect = vec![
            ("exit", J::Int(0)),
            ("json_file", s("{tmp}/report.json")),
            (
                "json_paths",
                complete_and(vec![
                    ("route.rho.pipeline", s(pipe)),
                    ("result.scalar", s(&known[id].to_string())),
                    ("result.known_answer", J::Bool(true)),
                    ("curve_id.slug", s(&r.slug)),
                ]),
            ),
        ];
        if !extra.is_empty() {
            expect.push((
                "json_contains",
                obj([("route.considered", considered(extra))]),
            ));
        }
        cases.push(case(
            &format!("{cid}-rho-{}", id.to_lowercase()),
            &format!("{id} under solve: rho: {purpose}, the logarithm recovered and verified"),
            copy(&format!("{{here}}/{id}-known.json")),
            plain_price_argv(),
            kv(expect),
            TIMEOUT,
            Vec::new(),
        ));
    }

    // C114: G1's curve in general Weierstrass form, with both points mapped to it.
    let g1 = rec("G1")?;
    let f1 = &fields[recs.iter().position(|r| r.id == "G1").expect("G1")];
    let e1 = Curve {
        f: f1,
        a: g1.a.clone(),
        b: g1.b.clone(),
    };
    let (g, k) = (&gens["G1"], known["G1"]);
    let (curve, to_general) = general_form(f1, &e1, &format!("{CASES_LABEL}/G1/general"));
    let gg = to_general(g);
    let q = e1.mul(k, &Some(g.clone())).ok_or("[k]G is the identity")?;
    let qg = to_general(&q);
    put(
        &mut files,
        "C114-general-weierstrass-g1.json".into(),
        &doc(
            "C114: G1's curve in general Weierstrass form, its points mapped",
            g1,
            &gg,
            obj([("point", point_doc(&qg))]),
            rho.clone(),
            Some(curve),
        ),
    );
    cases.push(case(
        "C114-general-weierstrass-extension",
        "a general Weierstrass curve over GF(271^3), isomorphic to G1's: converted, the conversion \
         recorded, and the logarithm of the mapped target recovered by rho",
        copy("{here}/C114-general-weierstrass-g1.json"),
        plain_price_argv(),
        obj([
            ("exit", J::Int(0)),
            ("json_file", s("{tmp}/report.json")),
            (
                "json_paths",
                complete_and(vec![
                    ("result.scalar", s(&k.to_string())),
                    ("route.rho.pipeline", s("rho-negation")),
                    ("conversion.from", s("general_weierstrass")),
                    ("conversion.to", s("short_weierstrass")),
                ]),
            ),
        ]),
        TIMEOUT,
        Vec::new(),
    ));

    // C115, C116: G1 under check, and at F2.
    cases.push(case(
        "C115-gaudry-check",
        "G1 under check: valid, with ic-gaudry-cubic and rho-negation admitted",
        copy("{here}/G1-known.json"),
        check_argv(),
        obj([
            ("exit", J::Int(0)),
            ("json_file", s("{tmp}/report.json")),
            ("json_paths", obj([("status", s("checks_passed"))])),
            (
                "json_contains",
                obj([("route.considered", gaudry_rho.clone())]),
            ),
        ]),
        120,
        Vec::new(),
    ));
    cases.push(case(
        "C116-gaudry-estimates",
        "G1 at fidelity F2: estimated, with an estimate for each arm",
        copy_set("{here}/G1-known.json", obj([("method.fidelity", s("F2"))])),
        price_argv(),
        obj([
            ("exit", J::Int(0)),
            ("json_file", s("{tmp}/report.json")),
            (
                "json_paths",
                obj([
                    ("status", s("estimated")),
                    ("estimate.ic.pipeline", s("ic-gaudry-cubic")),
                    ("estimate.rho.pipeline", s("rho-negation")),
                ]),
            ),
        ]),
        120,
        Vec::new(),
    ));

    // C117: the pipeline named on a binary document.
    cases.push(case(
        "C117-gaudry-on-a-binary-field",
        "ic-gaudry-cubic named on the smoke row's binary document: refused with its gate, \
         no-pipeline-for-field",
        copy_set(
            "{cases}/C009-translation-of-smoke-row.json",
            obj([(
                "method.index_calculus",
                obj([("pipeline", s("ic-gaudry-cubic")), ("recipe", s("auto"))]),
            )]),
        ),
        price_argv(),
        obj([
            ("exit", J::Int(3)),
            ("json_file", s("{tmp}/report.json")),
            (
                "json_paths",
                obj([
                    ("refusal.code", s("no-pipeline-for-field")),
                    ("refusal.class", s("unsupported")),
                ]),
            ),
        ]),
        120,
        Vec::new(),
    ));

    // C118: Gaudry's S4 solve by recipe.
    cases.push(case(
        "C118-gaudry-groebner-oracle",
        "G1 with recipe {oracle: groebner}: Gaudry's symmetrised S4 solve recovers the logarithm, \
         verified with the matched rho's",
        copy_set(
            "{here}/G1-known.json",
            obj([(
                "method.index_calculus",
                obj([
                    ("pipeline", s("ic-gaudry-cubic")),
                    ("recipe", obj([("oracle", s("groebner"))])),
                ]),
            )]),
        ),
        price_argv(),
        obj([
            ("exit", J::Int(0)),
            ("json_file", s("{tmp}/report.json")),
            (
                "json_paths",
                complete_and(vec![
                    ("result.scalar", s(&k.to_string())),
                    ("all_verified", J::Bool(true)),
                    ("route.ic.pipeline", s("ic-gaudry-cubic")),
                    ("ic.recipe.oracle", s("groebner")),
                ]),
            ),
        ]),
        TIMEOUT,
        Vec::new(),
    ));
    cases.sort_by(|x, y| {
        let id = |c: &J| c.get("id").and_then(J::as_str).unwrap_or("").to_string();
        id(x).cmp(&id(y))
    });
    Ok((files, cases))
}

/// Every file the cases' directory holds, by its path in the directory.  `cases.json`
/// names its generator by path and command; what pins the generator is
/// that it writes every file again byte for byte (`--check`, and this
/// module's test), not a hash of its source, which any edit would move.
pub fn case_texts(root: &Path) -> Result<Vec<(String, String)>, String> {
    let (files, cases) = build(root)?;
    let instances = std::fs::read(root.join(INSTANCES)).map_err(|e| format!("{INSTANCES}: {e}"))?;
    let mut out: Vec<(String, String)> = files
        .into_iter()
        .map(|(name, text)| (format!("params/{name}"), text))
        .collect();
    let rules = [
        "The rules of ../v2/cases.json hold. {here} is this directory's params/; {cases} is still \
         ../v2/params/. `icprog conformance` runs these after the earlier steps' cases.",
        "The supersedes rule: C050 names B5 as its `until`, and B5 is now two steps, B5a and B5b, \
         so C103 supersedes C050 from B5a on.",
        "params/C050-cubic-extension.json is B2's file, copied byte for byte for C103 and C107.",
        "params/*-T001.json are no case's: they are B5a's measurement 5's public points, frozen \
         here with the instances' known-answer documents. Every instance comes from \
         ../../rounds/B5a-extension-fields/instances.json, re-checked by this generator.",
    ];
    let text = json::dumps_utf8(
        &obj([
            (
                "suite",
                s("ic tool programme conformance suite v2, B5a's cases"),
            ),
            (
                "design",
                s("research/ic_tool_program/design/extension-fields.md; \
                   research/ic_tool_program/rounds/B5a-extension-fields/PROTOCOL.md"),
            ),
            (
                "includes",
                s(
                    "../v1/cases.json (step B0), ../v2/cases.json (B1) and the later steps' sets, \
                   run first, with the until and supersedes rules",
                ),
            ),
            ("rules", J::Arr(rules.iter().map(|r| s(r)).collect())),
            (
                "generator",
                obj([("path", s(GENERATOR)), ("command", s("icprog b5a cases"))]),
            ),
            (
                "instances",
                obj([
                    ("path", s(INSTANCES)),
                    ("sha256", s(&sha256_hex(&instances))),
                ]),
            ),
            ("cases", J::Arr(cases)),
        ]),
        1,
    ) + "\n";
    out.push(("cases.json".into(), text));
    let mut sorted = out.clone();
    sorted.sort_by(|x, y| x.0.cmp(&y.0));
    let sums: String = sorted
        .iter()
        .map(|(rel, t)| format!("{}  {rel}\n", sha256_hex(t.as_bytes())))
        .collect();
    out.push(("SHA256SUMS".into(), sums));
    Ok(out)
}

/// Writes the cases' directory, refusing to overwrite it; with `check`,
/// compares every file's bytes instead.
pub fn cases(root: &Path, check: bool) -> Result<(J, bool), String> {
    let dir = root.join(CASES_DIR);
    let out = case_texts(root)?;
    if check {
        let bad: Vec<J> = out
            .iter()
            .filter(|(rel, t)| std::fs::read_to_string(dir.join(rel)).ok().as_deref() != Some(t))
            .map(|(rel, _)| s(rel))
            .collect();
        let ok = bad.is_empty();
        return Ok((
            obj([
                ("files", J::Int(out.len() as i128)),
                ("mismatches", J::Arr(bad)),
            ]),
            ok,
        ));
    }
    if dir.join("cases.json").exists() {
        return Err("cases.json exists; B5a's cases are frozen (use --check)".into());
    }
    for (rel, t) in &out {
        let path = dir.join(rel);
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
        }
        std::fs::write(&path, t).map_err(|e| format!("{}: {e}", path.display()))?;
    }
    Ok((obj([("files", J::Int(out.len() as i128))]), true))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `#E` by counting, one abscissa at a time (curve_id.py's selftest).
    fn count(f: &Fpk, a: &[u128], b: &[u128]) -> u128 {
        let e = Curve {
            f,
            a: a.to_vec(),
            b: b.to_vec(),
        };
        let mut total = 1;
        for idx in 0..f.q {
            let x: El = (0..f.k).map(|i| idx / f.p.pow(i as u32) % f.p).collect();
            let rhs = e.rhs(&x);
            total += if rhs == f.zero() {
                1
            } else if f.is_square(&rhs) {
                2
            } else {
                0
            };
        }
        total
    }

    /// `#E(GF(p^k))` of a curve over GF(p), from `#E(GF(p))` by the trace
    /// recurrence.
    fn lifted(p: u128, a: u128, b: u128, k: u32) -> u128 {
        let f = Fpk::new(p, &[0, 1]).unwrap(); // only its prime field is used
        let mut n = 1i128;
        for x in 0..p {
            let rhs = (x * x % p * x + a * x + b) % p;
            n += if rhs == 0 {
                1
            } else if powm(rhs, (p - 1) / 2, p) == 1 {
                2
            } else {
                0
            };
        }
        let _ = f;
        let t1 = p as i128 + 1 - n;
        let (mut prev, mut cur) = (2i128, t1);
        for _ in 1..k {
            (prev, cur) = (cur, t1 * cur - p as i128 * prev);
        }
        (p as i128).pow(k) as u128 + 1 - cur as u128
    }

    #[test]
    fn the_field_counts_and_inverts_as_the_reference_does() {
        let f = Fpk::new(5, &[2, 0]).unwrap();
        assert!(rabin_irreducible(5, &[2, 0]));
        assert_eq!(count(&f, &[1, 0], &[1, 0]), 27);
        assert_eq!(lifted(5, 1, 1, 2), 27);
        let g = Fpk::new(7, &[4, 0, 0]).unwrap();
        assert!(rabin_irreducible(7, &[4, 0, 0]));
        assert_eq!(count(&g, &[2, 0, 0], &[3, 0, 0]), lifted(7, 2, 3, 3));
        for v in [[1u128, 2, 3], [0, 0, 5], [6, 0, 1]] {
            assert_eq!(g.mul(&v, &g.inv(&v).unwrap()), vec![1, 0, 0]);
        }
        // t^2 + 1 is reducible mod 5 (2^2 = -1), and so is t^3 - 1 mod 7;
        // t^3 - 2 is not, 2 being no cube mod 7.
        assert!(!rabin_irreducible(5, &[1, 0]));
        assert!(!rabin_irreducible(7, &[6, 0, 0]));
        assert!(rabin_irreducible(7, &[5, 0, 0]));
        // Square roots square back, at both 2-adic depths.
        for (f, x) in [(&f, vec![3u128, 4]), (&g, vec![1, 5, 2])] {
            let s = f.mul(&x, &x);
            let r = f.sqrt(&s).unwrap();
            assert_eq!(f.mul(&r, &r), s);
        }
    }

    #[test]
    fn the_identity_matches_the_references_vectors() {
        let e = extension_id(5, &[2, 0], &[1, 0], &[1, 0], 27).unwrap();
        assert!(e.icv1.starts_with("ICV1:fpk-5-2-"), "{}", e.icv1);
        assert!(e.slug.starts_with("icv1-fp3k2-tm1-"), "{}", e.slug);
        // A curve defined over GF(p) keeps its j-invariant: prime_id(5, 1,
        // 1, 9) has j = 2.
        assert_eq!(e.icv1.split(':').nth(4), Some("2,0"));
        let g = Fpk::new(7, &[4, 0, 0]).unwrap();
        let order = count(&g, &[1, 2, 0], &[3, 0, 5]);
        let e3 = extension_id(7, &[4, 0, 0], &[1, 2, 0], &[3, 0, 5], order).unwrap();
        assert!(e3.slug.starts_with("icv1-fp3k3-t"), "{}", e3.slug);
        assert!(e3.icv1.starts_with("ICV1:fpk-7-3-"), "{}", e3.icv1);
        assert!(extension_id(5, &[2, 0], &[0, 0], &[0, 0], 27).is_err());
        assert!(extension_id(5, &[2, 0], &[1, 0], &[1, 0], 40).is_err());
    }

    #[test]
    fn group_orders_match_counting() {
        let f = Fpk::new(7, &[4, 0, 0]).unwrap();
        for (i, (a, b)) in [([1u128, 2, 0], [3u128, 0, 5]), ([0, 1, 1], [2, 2, 0])]
            .into_iter()
            .enumerate()
        {
            let e = Curve {
                f: &f,
                a: a.to_vec(),
                b: b.to_vec(),
            };
            let got = group_order(&e, f.q, |j| Some(e.point(&format!("test/{i}/{j}"))));
            assert_eq!(got, Some(count(&f, &a, &b)), "{a:?} {b:?}");
        }
    }

    #[test]
    fn factors_and_primes_are_exact() {
        assert_eq!(
            factor(2 * 3 * 3 * 65537 * 1_000_003).unwrap(),
            vec![2, 3, 3, 65537, 1_000_003]
        );
        let (p1, p2) = (4_294_967_311u128, 1_099_511_627_791u128); // 2^32 + 15, 2^40 + 15
        assert!(prime_exact(p1).unwrap() && prime_exact(p2).unwrap());
        assert_eq!(factor(p1 * p2 * 4).unwrap(), vec![2, 2, p1, p2]);
        assert!(!prime_exact(3_215_031_751).unwrap()); // a strong pseudoprime to 2, 3, 5, 7
        assert!(prime_exact(PSI13).is_err());
        assert_eq!(subgroup(4 * p1, 1 << 62, "largest").unwrap(), None); // r below 4 sqrt q
        let m61 = (1u128 << 61) - 1; // a Mersenne prime, past R_MAX
        assert_eq!(subgroup(4 * m61, 1 << 62, "largest").unwrap(), None);
        assert_eq!(subgroup(4 * p2, 1 << 60, "largest").unwrap(), Some((p2, 4)));
        assert_eq!(subgroup(p2, 1 << 60, "cofactor").unwrap(), None);
    }

    fn root() -> &'static Path {
        Path::new(env!("CARGO_MANIFEST_DIR"))
    }

    #[test]
    fn the_frozen_cases_are_written_again_byte_for_byte() {
        let (summary, ok) = cases(root(), true).unwrap();
        assert!(ok, "{}", json::dumps(&summary, 1));
    }

    #[test]
    fn known_logarithms_replay_and_others_do_not() {
        let params = root().join(CASES_DIR).join("params");
        for id in ["G1", "E2", "E11", "B2"] {
            let doc = json::read(&params.join(format!("{id}-known.json"))).unwrap();
            let r: u128 = doc
                .at("subgroup")
                .unwrap()
                .at("order")
                .unwrap()
                .as_str()
                .unwrap()
                .parse()
                .unwrap();
            let k = known_log(&format!("{CASES_LABEL}/{id}/known_log"), r);
            let yes = replay(&doc, Some(&BigUint::from(k))).unwrap();
            assert_eq!(yes.get("replays"), Some(&J::Bool(true)), "{id}");
            assert_eq!(yes.get("known_matches"), Some(&J::Bool(true)), "{id}");
            let no = replay(&doc, Some(&BigUint::from(k + 1))).unwrap();
            assert_eq!(no.get("replays"), Some(&J::Bool(false)), "{id}");
        }
        // C050, B2's document, under a general cubic modulus.
        let c050 = json::read(
            &root()
                .join("research/ic_tool_program/conformance/v2-b2/params")
                .join(C050_DOC),
        )
        .unwrap();
        let yes = replay(&c050, Some(&BigUint::from(554_837u32))).unwrap();
        assert_eq!(yes.get("replays"), Some(&J::Bool(true)));
    }
}
