//! Smith's two rationality conditions for genus-3 hyperelliptic tractable
//! subgroups, computed natively, and the per-curve `full_success` they give.
//!
//! Source: B. Smith, "Isogenies and the discrete logarithm problem in
//! Jacobians of genus 3 hyperelliptic curves", J. Cryptology 22 (2009),
//! arXiv:0806.2995: §3 (tractable subgroups, Lemma 1), §5 (the matrix `M` of
//! eq. (2), the determinant (3), Proposition 2), §6 (`s(t)` of eq. (7),
//! Lemma 4, Proposition 4), §8 (Theorem 2, Examples 1 and 2).
//!
//! For `H: y^2 = F(x)` over `F_p` (`deg F` 7 or 8, `p` odd, `p < 2^31`):
//!
//! 1. Factor `F`.  A rational tractable subgroup is a Frobenius-stable
//!    partition of the eight Weierstrass points into pairs.  One pair per
//!    Galois orbit of pairs is built with both points in `K = F_p[z]/(g)` for
//!    an irreducible factor `g` (so `z` is a root); a pair joining two factors
//!    of equal degree needs one Cantor-Zassenhaus root in `K`.
//! 2. A pair's row of `M` (eq. (2)) lives over `K`.  The traces
//!    `Tr_{K/F_p}(z^j row)` span the `F_p`-points of the (Galois-stable) row
//!    space, so their `F_p`-kernel is `ker M` with an `F_p` basis `alpha, beta`
//!    (Galois descent).
//! 3. Proposition 2: `disc = (sum a_i b_{i+3})^2 - (sum a_i a_{i+3})(sum b_i b_{i+3})`
//!    (indices mod 6, `i = 0..5`) is a square in `F_p`.  Zero counts as a
//!    square and is tallied separately.
//! 4. The rational roots of `det(M_alpha + lambda M_beta) = 0` (projective in
//!    `lambda`) give lines `L`; Gaussian elimination gives `N` and `D`.
//! 5. Proposition 4: `f0 + f1 x + f2 x^2 = F mod G(t, x)`, `s(t)` from eq. (7);
//!    the isogeny is rational iff the leading coefficient of `s` is a square.
//!
//! Independent checks beside the formulas: the pair condition `g(W') = g(W'')`
//! on every pair, `deg g = 3`, Lemma 4 (`s / lc(s)` is a square in `F_p[t]`),
//! the square class of `lc(s)` against `s(t_0)` computed as a norm (which
//! also covers maps not in Smith's normal form), and, for small `p`, a point
//! count of Recillas' curve `X` from the fibres of `g`: `#X(F_p) = #H(F_p)`
//! when the isogeny is rational and `2p + 2 - #H(F_p)` when `J_X` is the
//! quadratic twist of `J_H / S` (Smith, Remark 4).
//!
//! Scope: these are Smith's conditions for his explicit construction; this
//! module builds no isogeny and no codomain beyond the point count of `X`.

use super::smith_kernel_coverage::{chebotarev_report, Fpk};
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};

pub type Poly = Vec<u64>;

// ---------------------------------------------------------------------------
// F_p with p odd and below 2^31, so products fit in u64.

#[inline]
fn md(a: u64, b: u64, p: u64) -> u64 {
    a * b % p
}

#[inline]
fn ad(a: u64, b: u64, p: u64) -> u64 {
    let s = a + b;
    if s >= p {
        s - p
    } else {
        s
    }
}

#[inline]
fn sb(a: u64, b: u64, p: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a + p - b
    }
}

#[inline]
fn ng(a: u64, p: u64) -> u64 {
    if a == 0 {
        0
    } else {
        p - a
    }
}

pub fn pow_mod(b: u64, mut e: u128, p: u64) -> u64 {
    let mut r = 1 % p;
    let mut b = b % p;
    while e > 0 {
        if e & 1 == 1 {
            r = md(r, b, p);
        }
        b = md(b, b, p);
        e >>= 1;
    }
    r
}

pub fn inv_mod(a: u64, p: u64) -> u64 {
    assert!(!a.is_multiple_of(p), "inverse of zero");
    pow_mod(a, (p - 2) as u128, p)
}

/// Quadratic character of `a` mod `p`: 0, 1 or -1.
pub fn chi(a: u64, p: u64) -> i8 {
    let a = a % p;
    if a == 0 {
        0
    } else if pow_mod(a, ((p - 1) / 2) as u128, p) == 1 {
        1
    } else {
        -1
    }
}

/// Square root mod `p` (the smaller representative), Tonelli-Shanks.
pub fn sqrt_mod(a: u64, p: u64) -> Option<u64> {
    let a = a % p;
    if a == 0 {
        return Some(0);
    }
    if chi(a, p) != 1 {
        return None;
    }
    let (mut q, mut s) = (p - 1, 0u32);
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }
    let mut z = 2;
    while chi(z, p) != -1 {
        z += 1;
    }
    let mut m = s;
    let mut c = pow_mod(z, q as u128, p);
    let mut t = pow_mod(a, q as u128, p);
    let mut r = pow_mod(a, q.div_ceil(2) as u128, p);
    while t != 1 {
        let mut i = 0;
        let mut tt = t;
        while tt != 1 {
            tt = md(tt, tt, p);
            i += 1;
        }
        let b = pow_mod(c, 1u128 << (m - i - 1), p);
        m = i;
        c = md(b, b, p);
        t = md(t, c, p);
        r = md(r, b, p);
    }
    Some(r.min(p - r))
}

/// The smallest quadratic non-residue mod `p`.
pub fn least_nonsquare(p: u64) -> u64 {
    (2..p)
        .find(|&c| chi(c, p) == -1)
        .expect("odd p has a non-residue")
}

/// SplitMix64: a fixed, documented generator so every run replays from its seed.
#[derive(Clone, Debug)]
pub struct SplitMix64(pub u64);

impl SplitMix64 {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `0..n` (rejection sampling, no modulo bias).
    pub fn below(&mut self, n: u64) -> u64 {
        let zone = u64::MAX - u64::MAX % n;
        loop {
            let x = self.next_u64();
            if x < zone {
                return x % n;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Polynomials over F_p, little-endian, trimmed.

fn trim(a: &mut Poly) {
    while a.last() == Some(&0) {
        a.pop();
    }
}

/// Degree, -1 for the zero polynomial.
pub fn deg(a: &[u64]) -> isize {
    a.iter().rposition(|&c| c != 0).map_or(-1, |i| i as isize)
}

fn padd(a: &[u64], b: &[u64], p: u64) -> Poly {
    let n = a.len().max(b.len());
    let mut r: Poly = (0..n)
        .map(|i| ad(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0), p))
        .collect();
    trim(&mut r);
    r
}

fn psub(a: &[u64], b: &[u64], p: u64) -> Poly {
    let n = a.len().max(b.len());
    let mut r: Poly = (0..n)
        .map(|i| sb(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0), p))
        .collect();
    trim(&mut r);
    r
}

fn pscale(a: &[u64], c: u64, p: u64) -> Poly {
    let mut r: Poly = a.iter().map(|&x| md(x, c, p)).collect();
    trim(&mut r);
    r
}

fn pmul(a: &[u64], b: &[u64], p: u64) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = ad(r[i + j], md(x, y, p), p);
        }
    }
    trim(&mut r);
    r
}

fn pdivrem(a: &[u64], b: &[u64], p: u64) -> (Poly, Poly) {
    let mut r = a.to_vec();
    trim(&mut r);
    let mut b = b.to_vec();
    trim(&mut b);
    assert!(!b.is_empty(), "division by zero polynomial");
    let db = b.len() - 1;
    if r.len() <= db {
        return (Vec::new(), r);
    }
    let li = inv_mod(b[db], p);
    let mut q = vec![0u64; r.len() - db];
    while r.len() > db {
        let k = r.len() - 1 - db;
        let c = md(r[r.len() - 1], li, p);
        q[k] = c;
        for (i, &bi) in b.iter().enumerate() {
            r[k + i] = sb(r[k + i], md(c, bi, p), p);
        }
        trim(&mut r);
    }
    trim(&mut q);
    (q, r)
}

fn prem(a: &[u64], b: &[u64], p: u64) -> Poly {
    pdivrem(a, b, p).1
}

fn pmonic(a: &[u64], p: u64) -> Poly {
    let mut a = a.to_vec();
    trim(&mut a);
    match a.last() {
        None => a,
        Some(&l) => pscale(&a, inv_mod(l, p), p),
    }
}

fn pgcd(a: &[u64], b: &[u64], p: u64) -> Poly {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    trim(&mut a);
    trim(&mut b);
    while !b.is_empty() {
        let r = prem(&a, &b, p);
        a = b;
        b = r;
    }
    pmonic(&a, p)
}

fn pmulmod(a: &[u64], b: &[u64], m: &[u64], p: u64) -> Poly {
    prem(&pmul(a, b, p), m, p)
}

fn ppowmod(a: &[u64], mut e: u128, m: &[u64], p: u64) -> Poly {
    let mut r = prem(&[1], m, p);
    let mut b = prem(a, m, p);
    while e > 0 {
        if e & 1 == 1 {
            r = pmulmod(&r, &b, m, p);
        }
        b = pmulmod(&b, &b, m, p);
        e >>= 1;
    }
    r
}

fn pderiv(a: &[u64], p: u64) -> Poly {
    let mut r: Poly = a
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, &c)| md(i as u64 % p, c, p))
        .collect();
    trim(&mut r);
    r
}

pub fn peval(a: &[u64], x: u64, p: u64) -> u64 {
    a.iter().rev().fold(0, |acc, &c| ad(md(acc, x, p), c, p))
}

/// Determinant over `F_p` by Gaussian elimination.
fn det_fp(mut a: Vec<Vec<u64>>, p: u64) -> u64 {
    let n = a.len();
    let mut det = 1u64;
    for c in 0..n {
        let Some(k) = (c..n).find(|&i| a[i][c] != 0) else {
            return 0;
        };
        if k != c {
            a.swap(k, c);
            det = ng(det, p);
        }
        det = md(det, a[c][c], p);
        let iv = inv_mod(a[c][c], p);
        let pivot = a[c].clone();
        for row in a.iter_mut().skip(c + 1) {
            let f = md(row[c], iv, p);
            if f == 0 {
                continue;
            }
            for j in c..n {
                row[j] = sb(row[j], md(f, pivot[j], p), p);
            }
        }
    }
    det
}

/// Norm of `h` in `F_p[x]/(m)` for monic `m`: the determinant of
/// multiplication by `h`, equal to the product of `h` over the roots of `m`.
fn norm_mod_monic(h: &[u64], m: &[u64], p: u64) -> u64 {
    let k = deg(m) as usize;
    let mut mat = vec![vec![0u64; k]; k];
    let mut cur = prem(h, m, p);
    for j in 0..k {
        for (i, &c) in cur.iter().enumerate() {
            mat[i][j] = c;
        }
        cur = prem(&pmul(&cur, &[0, 1], p), m, p);
    }
    det_fp(mat, p)
}

/// Reduced row echelon form in place (zero rows dropped); returns pivots.
fn rref(mat: &mut Vec<Vec<u64>>, p: u64) -> Vec<usize> {
    let rows = mat.len();
    if rows == 0 {
        return Vec::new();
    }
    let cols = mat[0].len();
    let mut piv = Vec::new();
    let mut r = 0;
    for c in 0..cols {
        if r == rows {
            break;
        }
        let Some(k) = (r..rows).find(|&i| mat[i][c] != 0) else {
            continue;
        };
        mat.swap(r, k);
        let iv = inv_mod(mat[r][c], p);
        for x in mat[r].iter_mut() {
            *x = md(*x, iv, p);
        }
        let pivot = mat[r].clone();
        for (i, row) in mat.iter_mut().enumerate() {
            if i == r || row[c] == 0 {
                continue;
            }
            let f = row[c];
            for j in 0..cols {
                row[j] = sb(row[j], md(f, pivot[j], p), p);
            }
        }
        piv.push(c);
        r += 1;
    }
    mat.truncate(r);
    piv
}

/// Rank and a basis of `{x : A x = 0}` for the rows of `A`.
fn nullspace(rows: &[Vec<u64>], cols: usize, p: u64) -> (usize, Vec<Vec<u64>>) {
    let mut m = rows.to_vec();
    let piv = rref(&mut m, p);
    let mut basis = Vec::new();
    for f in (0..cols).filter(|c| !piv.contains(c)) {
        let mut x = vec![0u64; cols];
        x[f] = 1;
        for (ri, &pc) in piv.iter().enumerate() {
            x[pc] = ng(m[ri][f], p);
        }
        basis.push(x);
    }
    (piv.len(), basis)
}

// ---------------------------------------------------------------------------
// Factorisation over F_p: distinct-degree, then Cantor-Zassenhaus.

fn edf(g: &[u64], d: usize, p: u64, rng: &mut SplitMix64, out: &mut Vec<Poly>) {
    let n = deg(g) as usize;
    if n == d {
        out.push(g.to_vec());
        return;
    }
    let e = ((p as u128).pow(d as u32) - 1) / 2;
    loop {
        let mut a: Poly = (0..n).map(|_| rng.below(p)).collect();
        trim(&mut a);
        if deg(&a) < 1 {
            continue;
        }
        let b = psub(&ppowmod(&a, e, g, p), &[1], p);
        let h = pgcd(&b, g, p);
        let dh = deg(&h);
        if dh > 0 && (dh as usize) < n {
            let q = pmonic(&pdivrem(g, &h, p).0, p);
            edf(&h, d, p, rng, out);
            edf(&q, d, p, rng, out);
            return;
        }
    }
}

/// Monic irreducible factors of `f` (any nonzero leading coefficient),
/// sorted by degree then coefficients; `None` if `f` is not squarefree.
pub fn factor_squarefree(f: &[u64], p: u64, rng: &mut SplitMix64) -> Option<Vec<Poly>> {
    let f = pmonic(f, p);
    if deg(&f) <= 0 {
        return Some(Vec::new());
    }
    if deg(&pgcd(&f, &pderiv(&f, p), p)) != 0 {
        return None;
    }
    let mut out = Vec::new();
    let mut rem = f;
    let mut h = prem(&[0, 1], &rem, p);
    let mut d = 0usize;
    while deg(&rem) > 0 {
        d += 1;
        let n = deg(&rem) as usize;
        if 2 * d > n {
            out.push(rem.clone());
            break;
        }
        h = ppowmod(&h, p as u128, &rem, p);
        let g = pgcd(&psub(&h, &[0, 1], p), &rem, p);
        if deg(&g) > 0 {
            edf(&g, d, p, rng, &mut out);
            rem = pmonic(&pdivrem(&rem, &g, p).0, p);
            h = prem(&h, &rem, p);
        }
    }
    out.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    Some(out)
}

// ---------------------------------------------------------------------------
// Extension fields K = F_p[z]/(m), deg m <= 8.

pub type El = [u64; 8];

#[derive(Clone, Debug)]
pub struct Ext {
    pub p: u64,
    pub d: usize,
    m: [u64; 9],
    tr: [u64; 8],
}

impl Ext {
    /// `F_p[z]/(m)` for monic irreducible `m` of degree 1..=8.
    pub fn new(m: &[u64], p: u64) -> Ext {
        let d = m.len() - 1;
        assert!(
            (1..=8).contains(&d) && m[d] == 1,
            "monic modulus of degree 1..=8"
        );
        let mut mm = [0u64; 9];
        mm[..=d].copy_from_slice(m);
        let mut k = Ext {
            p,
            d,
            m: mm,
            tr: [0; 8],
        };
        // Tr(z^j) = trace of multiplication by z^j on the basis z^i.
        let z = k.z();
        let mut zj = k.one();
        for j in 0..d {
            let mut t = 0;
            let mut zi = k.one();
            for i in 0..d {
                t = ad(t, k.mul(&zj, &zi)[i], p);
                zi = k.mul(&zi, &z);
            }
            k.tr[j] = t;
            zj = k.mul(&zj, &z);
        }
        k
    }

    /// `F_p` itself (`z = 0`).
    pub fn prime(p: u64) -> Ext {
        Ext::new(&[0, 1], p)
    }

    pub fn zero(&self) -> El {
        [0; 8]
    }

    pub fn c(&self, x: u64) -> El {
        let mut a = [0; 8];
        a[0] = x % self.p;
        a
    }

    pub fn one(&self) -> El {
        self.c(1)
    }

    /// The class of `z`: a root of the modulus.
    pub fn z(&self) -> El {
        if self.d == 1 {
            self.c(ng(self.m[0], self.p))
        } else {
            let mut a = [0; 8];
            a[1] = 1;
            a
        }
    }

    pub fn is_zero(&self, a: &El) -> bool {
        a[..self.d].iter().all(|&x| x == 0)
    }

    pub fn add(&self, a: &El, b: &El) -> El {
        let mut r = [0; 8];
        for i in 0..self.d {
            r[i] = ad(a[i], b[i], self.p);
        }
        r
    }

    pub fn sub(&self, a: &El, b: &El) -> El {
        let mut r = [0; 8];
        for i in 0..self.d {
            r[i] = sb(a[i], b[i], self.p);
        }
        r
    }

    pub fn neg(&self, a: &El) -> El {
        let mut r = [0; 8];
        for i in 0..self.d {
            r[i] = ng(a[i], self.p);
        }
        r
    }

    pub fn mul(&self, a: &El, b: &El) -> El {
        let (d, p) = (self.d, self.p);
        let mut r = [0u128; 15];
        for i in 0..d {
            if a[i] == 0 {
                continue;
            }
            for j in 0..d {
                r[i + j] += a[i] as u128 * b[j] as u128;
            }
        }
        let mut s = [0u64; 15];
        for k in 0..2 * d - 1 {
            s[k] = (r[k] % p as u128) as u64;
        }
        for k in (d..2 * d - 1).rev() {
            let c = s[k];
            if c == 0 {
                continue;
            }
            for i in 0..d {
                s[k - d + i] = sb(s[k - d + i], md(c, self.m[i], p), p);
            }
        }
        let mut out = [0u64; 8];
        out[..d].copy_from_slice(&s[..d]);
        out
    }

    pub fn pow(&self, a: &El, mut e: u128) -> El {
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

    pub fn frob(&self, a: &El) -> El {
        self.pow(a, self.p as u128)
    }

    pub fn trace(&self, a: &El) -> u64 {
        (0..self.d).fold(0, |s, j| ad(s, md(a[j], self.tr[j], self.p), self.p))
    }

    /// Inverse by the extended Euclidean algorithm in `F_p[z]`.
    pub fn inv(&self, a: &El) -> El {
        let p = self.p;
        let m: Poly = self.m[..=self.d].to_vec();
        let mut r0 = m.clone();
        let mut r1: Poly = a[..self.d].to_vec();
        trim(&mut r1);
        assert!(!r1.is_empty(), "inverse of zero");
        let (mut s0, mut s1): (Poly, Poly) = (Vec::new(), vec![1]);
        while deg(&r1) > 0 {
            let (q, r) = pdivrem(&r0, &r1, p);
            let s2 = psub(&s0, &pmul(&q, &s1, p), p);
            r0 = r1;
            r1 = r;
            s0 = s1;
            s1 = s2;
        }
        let s = prem(&pscale(&s1, inv_mod(r1[0], p), p), &m, p);
        let mut out = [0; 8];
        out[..s.len()].copy_from_slice(&s);
        out
    }

    /// `|K| = p^d`, if it fits in a `u128`.
    pub fn size(&self) -> Option<u128> {
        (self.p as u128).checked_pow(self.d as u32)
    }
}

// Polynomials over K, used only to find a root of a factor that splits in K.

fn kp_trim(k: &Ext, a: &mut Vec<El>) {
    while a.last().is_some_and(|x| k.is_zero(x)) {
        a.pop();
    }
}

fn kp_divrem(k: &Ext, a: &[El], b: &[El]) -> (Vec<El>, Vec<El>) {
    let mut r = a.to_vec();
    kp_trim(k, &mut r);
    let mut b = b.to_vec();
    kp_trim(k, &mut b);
    let db = b.len() - 1;
    if r.len() <= db {
        return (Vec::new(), r);
    }
    let li = k.inv(&b[db]);
    let mut q = vec![k.zero(); r.len() - db];
    while r.len() > db {
        let s = r.len() - 1 - db;
        let c = k.mul(&r[r.len() - 1], &li);
        q[s] = c;
        for (i, bi) in b.iter().enumerate() {
            r[s + i] = k.sub(&r[s + i], &k.mul(&c, bi));
        }
        kp_trim(k, &mut r);
    }
    (q, r)
}

fn kp_mulmod(k: &Ext, a: &[El], b: &[El], m: &[El]) -> Vec<El> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut r = vec![k.zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            r[i + j] = k.add(&r[i + j], &k.mul(x, y));
        }
    }
    kp_divrem(k, &r, m).1
}

fn kp_monic(k: &Ext, a: &[El]) -> Vec<El> {
    let mut a = a.to_vec();
    kp_trim(k, &mut a);
    if let Some(l) = a.last().copied() {
        let li = k.inv(&l);
        for x in a.iter_mut() {
            *x = k.mul(x, &li);
        }
    }
    a
}

fn kp_gcd(k: &Ext, a: &[El], b: &[El]) -> Vec<El> {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    kp_trim(k, &mut a);
    kp_trim(k, &mut b);
    while !b.is_empty() {
        let r = kp_divrem(k, &a, &b).1;
        a = b;
        b = r;
    }
    kp_monic(k, &a)
}

/// One root in `K` of a monic polynomial over `K` with distinct roots, all in `K`.
fn kp_one_root(k: &Ext, f: &[El], rng: &mut SplitMix64) -> El {
    let mut f = kp_monic(k, f);
    let e = (k.size().expect("|K| fits in u128") - 1) / 2;
    loop {
        let n = f.len() - 1;
        if n == 1 {
            return k.neg(&f[0]);
        }
        let mut a: Vec<El> = (0..n)
            .map(|_| {
                let mut x = [0; 8];
                for c in x.iter_mut().take(k.d) {
                    *c = rng.below(k.p);
                }
                x
            })
            .collect();
        kp_trim(k, &mut a);
        if a.len() < 2 {
            continue;
        }
        let mut r = vec![k.one()];
        let mut b = kp_divrem(k, &a, &f).1;
        let mut ee = e;
        while ee > 0 {
            if ee & 1 == 1 {
                r = kp_mulmod(k, &r, &b, &f);
            }
            b = kp_mulmod(k, &b, &b, &f);
            ee >>= 1;
        }
        if r.is_empty() {
            r.push(k.zero());
        }
        r[0] = k.sub(&r[0], &k.one());
        let h = kp_gcd(k, &r, &f);
        let dh = h.len() - 1;
        if dh > 0 && dh < n {
            let q = kp_monic(k, &kp_divrem(k, &f, &h).0);
            f = if dh < q.len() { h } else { q };
        }
    }
}

// ---------------------------------------------------------------------------
// Weierstrass points and tractable subgroups (Smith §3).

/// A point `(u : v)` of `P^1(K)`; finite `x` is `(x : 1)`, infinity `(1 : 0)`.
#[derive(Clone, Debug)]
pub struct KPt {
    pub u: El,
    pub v: El,
}

/// One pair of Weierstrass points, both over `k`, standing for its Galois orbit.
#[derive(Clone, Debug)]
pub struct PairRep {
    pub k: Ext,
    pub a: KPt,
    pub b: KPt,
}

/// A rational tractable subgroup: one representative pair per Galois orbit.
#[derive(Clone, Debug)]
pub struct Subgroup {
    pub pairs: Vec<PairRep>,
    pub label: String,
}

#[derive(Clone, Debug)]
pub struct Weierstrass {
    pub p: u64,
    /// `F(x) = F~(x, 1)`, degree 7 or 8.
    pub f: Poly,
    /// Monic irreducible factors of `F(x)`.
    pub factors: Vec<Poly>,
    /// `deg F = 7`: `(1 : 0)` is a Weierstrass point.
    pub infinite: bool,
}

impl Weierstrass {
    pub fn new(f: &[u64], p: u64, rng: &mut SplitMix64) -> Option<Weierstrass> {
        let mut f = f.to_vec();
        trim(&mut f);
        let d = deg(&f);
        assert!(d == 7 || d == 8, "genus 3 needs deg F in {{7, 8}}");
        let factors = factor_squarefree(&f, p, rng)?;
        Some(Weierstrass {
            p,
            f,
            factors,
            infinite: d == 7,
        })
    }

    /// Degrees of the factors of the octic form, sorted (infinity counts as 1).
    pub fn cycle_type(&self) -> Vec<u32> {
        let mut ct: Vec<u32> = self.factors.iter().map(|g| (g.len() - 1) as u32).collect();
        if self.infinite {
            ct.push(1);
        }
        ct.sort_unstable();
        ct
    }
}

fn finite_pt(k: &Ext, x: El) -> KPt {
    KPt { u: x, v: k.one() }
}

fn rational_matchings(
    rat: &[KPt],
    fp: &Ext,
    used: &mut Vec<bool>,
    cur: &mut Vec<(usize, usize)>,
    out: &mut Vec<(Vec<PairRep>, String)>,
) {
    let Some(i) = used.iter().position(|&u| !u) else {
        let pairs = cur
            .iter()
            .map(|&(i, j)| PairRep {
                k: fp.clone(),
                a: rat[i].clone(),
                b: rat[j].clone(),
            })
            .collect();
        let label = cur
            .iter()
            .map(|(i, j)| format!("r{i}{j}"))
            .collect::<Vec<_>>()
            .join(".");
        out.push((pairs, label));
        return;
    };
    used[i] = true;
    for j in i + 1..rat.len() {
        if used[j] {
            continue;
        }
        used[j] = true;
        cur.push((i, j));
        rational_matchings(rat, fp, used, cur, out);
        cur.pop();
        used[j] = false;
    }
    used[i] = false;
}

struct DegreeClass<'a> {
    d: usize,
    gs: Vec<&'a Poly>,
    exts: Vec<Ext>,
    roots: HashMap<(usize, usize), El>,
}

fn partial_matchings(
    dc: &mut DegreeClass,
    done: &mut Vec<bool>,
    cur: &mut Vec<(PairRep, String)>,
    rng: &mut SplitMix64,
    out: &mut Vec<(Vec<PairRep>, String)>,
) {
    let Some(i) = done.iter().position(|&u| !u) else {
        out.push((
            cur.iter().map(|c| c.0.clone()).collect(),
            cur.iter()
                .map(|c| c.1.clone())
                .collect::<Vec<_>>()
                .join("."),
        ));
        return;
    };
    let d = dc.d;
    done[i] = true;
    let k = dc.exts[i].clone();
    let z = k.z();
    if d.is_multiple_of(2) {
        // The pairs {r, r^(p^(d/2))} inside the d-cycle of g_i.
        let mut w = z;
        for _ in 0..d / 2 {
            w = k.frob(&w);
        }
        cur.push((
            PairRep {
                k: k.clone(),
                a: finite_pt(&k, z),
                b: finite_pt(&k, w),
            },
            format!("s{d}:{i}"),
        ));
        partial_matchings(dc, done, cur, rng, out);
        cur.pop();
    }
    for j in i + 1..dc.gs.len() {
        if done[j] {
            continue;
        }
        let s0 = *dc.roots.entry((i, j)).or_insert_with(|| {
            let f: Vec<El> = dc.gs[j].iter().map(|&c| k.c(c)).collect();
            kp_one_root(&k, &f, rng)
        });
        done[j] = true;
        let mut s = s0;
        for shift in 0..d {
            cur.push((
                PairRep {
                    k: k.clone(),
                    a: finite_pt(&k, z),
                    b: finite_pt(&k, s),
                },
                format!("x{d}:{i}{j}+{shift}"),
            ));
            partial_matchings(dc, done, cur, rng, out);
            cur.pop();
            s = k.frob(&s);
        }
        done[j] = false;
    }
    done[i] = false;
}

/// Every rational tractable subgroup `S(H)` (Smith §3, Lemma 1).
pub fn tractable_subgroups(w: &Weierstrass, rng: &mut SplitMix64) -> Vec<Subgroup> {
    let p = w.p;
    let fp = Ext::prime(p);
    let mut rat: Vec<KPt> = Vec::new();
    let mut by_deg: BTreeMap<usize, Vec<&Poly>> = BTreeMap::new();
    for g in &w.factors {
        let d = g.len() - 1;
        if d == 1 {
            rat.push(finite_pt(&fp, fp.c(ng(g[0], p))));
        } else {
            by_deg.entry(d).or_default().push(g);
        }
    }
    if w.infinite {
        rat.push(KPt {
            u: fp.one(),
            v: fp.zero(),
        });
    }
    if rat.len() % 2 == 1 {
        return Vec::new();
    }
    let mut classes: Vec<Vec<(Vec<PairRep>, String)>> = Vec::new();
    let mut opts = Vec::new();
    rational_matchings(
        &rat,
        &fp,
        &mut vec![false; rat.len()],
        &mut Vec::new(),
        &mut opts,
    );
    classes.push(opts);
    for (&d, gs) in &by_deg {
        let mut dc = DegreeClass {
            d,
            exts: gs.iter().map(|g| Ext::new(g, p)).collect(),
            gs: gs.clone(),
            roots: HashMap::new(),
        };
        let mut opts = Vec::new();
        let n = dc.gs.len();
        partial_matchings(
            &mut dc,
            &mut vec![false; n],
            &mut Vec::new(),
            rng,
            &mut opts,
        );
        if opts.is_empty() {
            return Vec::new();
        }
        classes.push(opts);
    }
    let mut acc: Vec<(Vec<PairRep>, String)> = vec![(Vec::new(), String::new())];
    for opts in classes {
        let mut next = Vec::with_capacity(acc.len() * opts.len());
        for (base, lab) in &acc {
            for (o, l) in &opts {
                let mut v = base.clone();
                v.extend(o.iter().cloned());
                let label = match (lab.is_empty(), l.is_empty()) {
                    (true, _) => l.clone(),
                    (_, true) => lab.clone(),
                    _ => format!("{lab}.{l}"),
                };
                next.push((v, label));
            }
        }
        acc = next;
    }
    acc.into_iter()
        .map(|(pairs, label)| Subgroup { pairs, label })
        .collect()
}

// ---------------------------------------------------------------------------
// Smith §5: the matrix M, its kernel, Proposition 2, the trigonal map.

/// `(a, b, c)` of `F_i = a u^2 + b uv + c v^2` vanishing at the pair.
pub fn pair_quadratic(pr: &PairRep) -> [El; 3] {
    let k = &pr.k;
    let (x1, y1, x2, y2) = (&pr.a.u, &pr.a.v, &pr.b.u, &pr.b.v);
    // (y1 u - x1 v)(y2 u - x2 v)
    let a = k.mul(y1, y2);
    let b = k.neg(&k.add(&k.mul(y1, x2), &k.mul(x1, y2)));
    let c = k.mul(x1, x2);
    [a, b, c]
}

/// The row `(a^2, ab, ac, c^2, -cb, b^2 - ac)` of eq. (2).
pub fn m_row(k: &Ext, q: &[El; 3]) -> [El; 6] {
    let [a, b, c] = q;
    let ac = k.mul(a, c);
    [
        k.mul(a, a),
        k.mul(a, b),
        ac,
        k.mul(c, c),
        k.neg(&k.mul(c, b)),
        k.sub(&k.mul(b, b), &ac),
    ]
}

/// `F_p` spanning set of the row space of `M` (traces of `z^j` times each row).
pub fn m_rows_fp(sg: &Subgroup) -> Vec<Vec<u64>> {
    let mut rows = Vec::new();
    for pr in &sg.pairs {
        let k = &pr.k;
        let row = m_row(k, &pair_quadratic(pr));
        let z = k.z();
        let mut zj = k.one();
        for _ in 0..k.d {
            rows.push(row.iter().map(|x| k.trace(&k.mul(&zj, x))).collect());
            zj = k.mul(&zj, &z);
        }
    }
    rows
}

/// `sum_{i=0..5} a_i b_{i+3}` (indices mod 6).
pub fn qform(a: &[u64], b: &[u64], p: u64) -> u64 {
    (0..6).fold(0, |s, i| ad(s, md(a[i], b[(i + 3) % 6], p), p))
}

/// The matrix `M_gamma` of Lemma 2 (2) whose kernel is the line of `gamma`.
pub fn m_gamma(g: &[u64], p: u64) -> Vec<Vec<u64>> {
    let n = |x: u64| ng(x, p);
    vec![
        vec![0, n(g[3]), n(g[4]), n(g[5])],
        vec![g[3], 0, n(g[2]), g[1]],
        vec![g[4], g[2], 0, n(g[0])],
        vec![g[5], n(g[1]), g[0], 0],
    ]
}

/// A trigonal map `g = (N~ : D~)` with `N~, D~` binary cubic forms given by
/// their coefficients of `(u^3, u^2 v, u v^2, v^3)`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct TrigMap {
    pub n: [u64; 4],
    pub d: [u64; 4],
    /// Smith's form `N = x^3 + n1 x + n0`, `D = x^2 + d1 x + d0`.
    pub normal_form: bool,
}

impl TrigMap {
    /// `(n1, n0, d1, d0)` when in normal form.
    pub fn smith_coeffs(&self) -> Option<(u64, u64, u64, u64)> {
        self.normal_form
            .then_some((self.n[2], self.n[3], self.d[2], self.d[3]))
    }
}

/// The line `L` of a point `gamma` of `Gr(1,3)` and the projection from it.
pub fn trig_map_of(gamma: &[u64], p: u64) -> Option<TrigMap> {
    let mut m = m_gamma(gamma, p);
    let piv = rref(&mut m, p);
    if piv.len() != 2 {
        return None;
    }
    let mut n = [0u64; 4];
    let mut d = [0u64; 4];
    n.copy_from_slice(&m[0]);
    d.copy_from_slice(&m[1]);
    Some(TrigMap {
        n,
        d,
        normal_form: piv == [0, 1],
    })
}

fn cubic_at(k: &Ext, c: &[u64; 4], pt: &KPt) -> El {
    let (u, v) = (&pt.u, &pt.v);
    let u2 = k.mul(u, u);
    let v2 = k.mul(v, v);
    let t0 = k.mul(&k.c(c[0]), &k.mul(&u2, u));
    let t1 = k.mul(&k.c(c[1]), &k.mul(&u2, v));
    let t2 = k.mul(&k.c(c[2]), &k.mul(u, &v2));
    let t3 = k.mul(&k.c(c[3]), &k.mul(&v2, v));
    k.add(&k.add(&t0, &t1), &k.add(&t2, &t3))
}

/// `g(P) = g(Q)` for the pair, as `N~(P) D~(Q) - N~(Q) D~(P) = 0` in `K`.
pub fn pair_condition(map: &TrigMap, pr: &PairRep) -> bool {
    let k = &pr.k;
    let l = k.mul(&cubic_at(k, &map.n, &pr.a), &cubic_at(k, &map.d, &pr.b));
    let r = k.mul(&cubic_at(k, &map.n, &pr.b), &cubic_at(k, &map.d, &pr.a));
    k.is_zero(&k.sub(&l, &r))
}

fn cubic_poly(c: &[u64; 4]) -> Poly {
    let mut q = vec![c[3], c[2], c[1], c[0]];
    trim(&mut q);
    q
}

/// `N~` and `D~` have no common zero on `P^1`, so `g` has degree 3.
pub fn degree_three(map: &TrigMap, p: u64) -> bool {
    if map.n[0] == 0 && map.d[0] == 0 {
        return false;
    }
    let (a, b) = (cubic_poly(&map.n), cubic_poly(&map.d));
    if a.is_empty() || b.is_empty() {
        return false;
    }
    deg(&pgcd(&a, &b, p)) == 0
}

// ---------------------------------------------------------------------------
// Smith §6: f0, f1, f2, s(t) (eq. (7)) and Proposition 4.

/// `[f0, f1, f2]` in `F_p[t]` with `f0 + f1 x + f2 x^2 = F(x) mod G(t, x)`,
/// `G = x^3 - t x^2 + (n1 - d1 t) x + (n0 - d0 t)`.
pub fn f_mod_g(f: &[u64], n1: u64, n0: u64, d1: u64, d0: u64, p: u64) -> [Poly; 3] {
    let g2: Poly = vec![0, ng(1, p)];
    let g1: Poly = vec![n1, ng(d1, p)];
    let g0: Poly = vec![n0, ng(d0, p)];
    let mut c: Vec<Poly> = f
        .iter()
        .map(|&a| {
            let mut v = vec![a];
            trim(&mut v);
            v
        })
        .collect();
    while c.len() < 3 {
        c.push(Vec::new());
    }
    for k in (3..c.len()).rev() {
        let ck = std::mem::take(&mut c[k]);
        // x^3 = -g2 x^2 - g1 x - g0
        c[k - 1] = psub(&c[k - 1], &pmul(&ck, &g2, p), p);
        c[k - 2] = psub(&c[k - 2], &pmul(&ck, &g1, p), p);
        c[k - 3] = psub(&c[k - 3], &pmul(&ck, &g0, p), p);
    }
    [c[0].clone(), c[1].clone(), c[2].clone()]
}

/// `(g0, g1, g2)` of `G(t, x)` in normal form.
pub fn g_coeffs(n1: u64, n0: u64, d1: u64, d0: u64, p: u64) -> [Poly; 3] {
    let mut g0 = vec![n0, ng(d0, p)];
    let mut g1 = vec![n1, ng(d1, p)];
    trim(&mut g0);
    trim(&mut g1);
    [g0, g1, vec![0, ng(1, p)]]
}

fn term(coef: i64, fs: &[&Poly], p: u64) -> Poly {
    let mut acc: Poly = vec![(coef.rem_euclid(p as i64)) as u64];
    trim(&mut acc);
    for f in fs {
        acc = pmul(&acc, f, p);
    }
    acc
}

/// `s(t)` of Lemma 4, eq. (7).
pub fn smith_s(fs: &[Poly; 3], gs: &[Poly; 3], p: u64) -> Poly {
    let (f0, f1, f2) = (&fs[0], &fs[1], &fs[2]);
    let (g0, g1, g2) = (&gs[0], &gs[1], &gs[2]);
    let terms = [
        term(1, &[f0, f0, f0], p),
        term(-1, &[f0, f0, f1, g2], p),
        term(-2, &[f0, f0, f2, g1], p),
        term(1, &[f0, f0, f2, g2, g2], p),
        term(1, &[f0, f1, f1, g1], p),
        term(3, &[f0, f1, f2, g0], p),
        term(-1, &[f0, f1, f2, g1, g2], p),
        term(-2, &[f0, f2, f2, g0, g2], p),
        term(1, &[f0, f2, f2, g1, g1], p),
        term(-1, &[f1, f1, f1, g0], p),
        term(1, &[f1, f1, f2, g0, g2], p),
        term(-1, &[f1, f2, f2, g0, g1], p),
        term(1, &[f2, f2, f2, g0, g0], p),
    ];
    terms.iter().fold(Vec::new(), |acc, t| padd(&acc, t, p))
}

/// Proposition 3's `delta_0, delta_2, delta_4` (for the Example 1 control).
pub fn smith_deltas(fs: &[Poly; 3], gs: &[Poly; 3], p: u64) -> [Poly; 3] {
    let (f0, f1, f2) = (&fs[0], &fs[1], &fs[2]);
    let (g0, g1, g2) = (&gs[0], &gs[1], &gs[2]);
    let sum = |ts: Vec<Poly>| ts.iter().fold(Vec::new(), |acc, t| padd(&acc, t, p));
    let d4 = sum(vec![
        term(-27, &[g0, g0], p),
        term(18, &[g0, g1, g2], p),
        term(-4, &[g0, g2, g2, g2], p),
        term(-4, &[g1, g1, g1], p),
        term(1, &[g1, g1, g2, g2], p),
    ]);
    let d2 = sum(vec![
        term(12, &[f0, g1], p),
        term(-4, &[f0, g2, g2], p),
        term(-18, &[f1, g0], p),
        term(2, &[f1, g1, g2], p),
        term(12, &[f2, g0, g2], p),
        term(-4, &[f2, g1, g1], p),
    ]);
    let d0 = sum(vec![term(-4, &[f0, f2], p), term(1, &[f1, f1], p)]);
    [d0, d2, d4]
}

/// Lemma 4: `s / lc(s)` is the square of a monic polynomial over `F_p`.
pub fn is_lc_times_square(s: &[u64], p: u64) -> bool {
    let n = deg(s);
    if n < 0 || n % 2 == 1 {
        return false;
    }
    let sm = pmonic(s, p);
    let m = (n / 2) as usize;
    let mut h = vec![0u64; m + 1];
    h[m] = 1;
    let inv2 = inv_mod(2, p);
    for k in 1..=m {
        let mut acc = 0;
        for i in (m - k + 1)..=m {
            let j = 2 * m - k - i;
            if j > m - k && j <= m {
                acc = ad(acc, md(h[i], h[j], p), p);
            }
        }
        h[m - k] = md(sb(sm[2 * m - k], acc, p), inv2, p);
    }
    pmul(&h, &h, p) == sm
}

/// Square class of `lc(s)` from a value: `s(t_0) = prod F(x)` over the fibre
/// `g(x) = t_0`, a norm, for the first `t_0 in F_p` where it is defined and
/// nonzero.  Works for maps not in normal form.  0 if no `t_0` qualifies.
pub fn eval_alpha_chi(f: &[u64], map: &TrigMap, p: u64) -> i8 {
    for t0 in 0..p {
        let mut gc: Poly = (0..4)
            .map(|i| sb(map.n[3 - i], md(t0, map.d[3 - i], p), p))
            .collect();
        trim(&mut gc);
        if deg(&gc) != 3 {
            continue;
        }
        let gm = pmonic(&gc, p);
        let v = norm_mod_monic(&prem(f, &gm, p), &gm, p);
        if v != 0 {
            return chi(v, p);
        }
    }
    0
}

// ---------------------------------------------------------------------------
// Point count of Recillas' curve X from the fibres of g (small p).

/// `#X(F_p)` for the curve of the trigonal construction on the normal-form
/// map `map`, for `deg F = 8`; `None` where a fibre is degenerate (a triple
/// point, a Weierstrass point at a ramification point of `g`, ...).
pub fn x_point_count(f: &[u64], map: &TrigMap, p: u64) -> Option<u64> {
    let (n1, n0, d1, d0) = map.smith_coeffs()?;
    if deg(f) != 8 {
        return None;
    }
    let chi_f: Vec<i8> = (0..p).map(|x| chi(peval(f, x, p), p)).collect();
    let chi_inf = chi(f[8], p);
    let fd = pderiv(f, p);
    let np: Poly = {
        let mut v = vec![n0, n1, 0, 1];
        trim(&mut v);
        v
    };
    let dp: Poly = {
        let mut v = vec![d0, d1, 1];
        trim(&mut v);
        v
    };
    // Numerator of g' = (N' D - N D') / D^2.
    let gnum = psub(
        &pmul(&pderiv(&np, p), &dp, p),
        &pmul(&np, &pderiv(&dp, p), p),
        p,
    );
    let roots = |g: &[u64]| -> Vec<u64> { (0..p).filter(|&x| peval(g, x, p) == 0).collect() };
    let mut total = 0u64;
    for t0 in 0..p {
        let g: Poly = vec![
            sb(n0, md(d0, t0, p), p),
            sb(n1, md(d1, t0, p), p),
            ng(t0, p),
            1,
        ];
        let w = pgcd(&g, f, p);
        let gd = pgcd(&g, &pderiv(&g, p), p);
        match deg(&w) {
            0 => match deg(&gd) {
                0 => {
                    let rs = roots(&g);
                    total += match rs.len() {
                        3 => {
                            let c: Vec<i8> = rs.iter().map(|&x| chi_f[x as usize]).collect();
                            if c[0] == c[1] && c[1] == c[2] {
                                4
                            } else {
                                0
                            }
                        }
                        1 => {
                            let q2 = pdivrem(&g, &[ng(rs[0], p), 1], p).0;
                            if chi(norm_mod_monic(&prem(f, &q2, p), &q2, p), p) == 1 {
                                2
                            } else {
                                0
                            }
                        }
                        0 => 1,
                        _ => return None,
                    };
                }
                1 => {
                    let xb = ng(gd[0], p);
                    let lin = pdivrem(&g, &pmul(&gd, &gd, p), p).0;
                    let x3 = ng(lin[0], p);
                    let c = chi_f[xb as usize] * chi_f[x3 as usize];
                    total += 1 + if c == 1 { 2 } else { 0 };
                }
                _ => return None,
            },
            2 => {
                if deg(&gd) != 0 {
                    return None;
                }
                let c = md(
                    norm_mod_monic(&prem(&fd, &w, p), &w, p),
                    norm_mod_monic(&prem(&gnum, &w, p), &w, p),
                    p,
                );
                match chi(c, p) {
                    1 => total += 2,
                    -1 => {}
                    _ => return None,
                }
            }
            _ => return None,
        }
    }
    // The fibre over t = infinity: x = infinity and the roots of D.
    let w = pgcd(&dp, f, p);
    let dd = pgcd(&dp, &pderiv(&dp, p), p);
    match deg(&w) {
        0 => {
            if deg(&dd) == 0 {
                let rs = roots(&dp);
                total += match rs.len() {
                    2 => {
                        let (a, b) = (chi_f[rs[0] as usize], chi_f[rs[1] as usize]);
                        if a == b && b == chi_inf {
                            4
                        } else {
                            0
                        }
                    }
                    0 => {
                        if chi(norm_mod_monic(&prem(f, &dp, p), &dp, p), p) == 1 {
                            2
                        } else {
                            0
                        }
                    }
                    _ => return None,
                };
            } else {
                let delta = ng(dd[0], p);
                let c = chi_f[delta as usize] * chi_inf;
                total += 1 + if c == 1 { 2 } else { 0 };
            }
        }
        2 => {
            if deg(&dd) != 0 {
                return None;
            }
            let c = md(
                md(
                    norm_mod_monic(&prem(&fd, &dp, p), &dp, p),
                    norm_mod_monic(&prem(&np, &dp, p), &dp, p),
                    p,
                ),
                norm_mod_monic(&prem(&pderiv(&dp, p), &dp, p), &dp, p),
                p,
            );
            match chi(c, p) {
                1 => total += 2,
                -1 => {}
                _ => return None,
            }
        }
        _ => return None,
    }
    Some(total)
}

/// `#H(F_p)` for `y^2 = F(x)`, `deg F` 7 or 8 (small `p`).
pub fn h_point_count(f: &[u64], p: u64) -> u64 {
    let affine: i64 = (0..p).map(|x| 1 + chi(peval(f, x, p), p) as i64).sum();
    let inf: i64 = if deg(f) == 8 {
        1 + chi(f[8], p) as i64
    } else {
        1
    };
    (affine + inf) as u64
}

// ---------------------------------------------------------------------------
// Per-subgroup and per-curve outcomes.

#[derive(Clone, Debug, Default, Serialize)]
pub struct MapOutcome {
    pub map: TrigMap,
    pub degree3: bool,
    pub pairs_ok: bool,
    /// Leading coefficient of `s(t)` (normal form only).
    pub alpha: Option<u64>,
    /// Square class used for Proposition 4: eq. (7) in normal form, the
    /// evaluated norm otherwise.
    pub alpha_chi: i8,
    /// Square class of `s(t_0)` (independent of eq. (7)).
    pub eval_chi: i8,
    /// Lemma 4 holds (normal form only).
    pub lemma4_ok: Option<bool>,
}

impl MapOutcome {
    pub fn valid(&self) -> bool {
        self.degree3 && self.pairs_ok && self.alpha_chi != 0
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct SubgroupOutcome {
    pub label: String,
    /// Rank of the `F_p` row space of `M` (Smith: 4).
    pub m_rank: usize,
    pub disc: u64,
    pub disc_chi: i8,
    /// `det(M_alpha + lambda M_beta)` vanishes identically.
    pub pencil_in_grassmannian: bool,
    pub maps: Vec<MapOutcome>,
}

impl SubgroupOutcome {
    /// Proposition 2: the discriminant (4) is a square (zero included).
    pub fn trig_rational(&self) -> bool {
        self.m_rank == 4 && self.disc_chi >= 0
    }

    /// Proposition 4 for some rational trigonal map of this subgroup, on the
    /// curve (`twist = false`) or its quadratic twist (`lc(s)` times `c^3`).
    pub fn iso_rational(&self, twist: bool) -> bool {
        let want = if twist { -1 } else { 1 };
        self.trig_rational() && self.maps.iter().any(|m| m.valid() && m.alpha_chi == want)
    }

    /// The two rational trigonal maps give different square classes.
    pub fn lambda_disagree(&self) -> bool {
        let v: Vec<i8> = self
            .maps
            .iter()
            .filter(|m| m.valid())
            .map(|m| m.alpha_chi)
            .collect();
        v.len() == 2 && v[0] != v[1]
    }

    /// Proposition 2 holds but no valid map could be built.
    pub fn construction_failed(&self) -> bool {
        self.trig_rational() && !self.maps.iter().any(|m| m.valid())
    }
}

/// Run Propositions 2 and 4 on one tractable subgroup.
pub fn analyze_subgroup(f: &[u64], sg: &Subgroup, p: u64) -> SubgroupOutcome {
    let rows = m_rows_fp(sg);
    let (rank, basis) = nullspace(&rows, 6, p);
    let mut out = SubgroupOutcome {
        label: sg.label.clone(),
        m_rank: rank,
        ..Default::default()
    };
    if rank != 4 {
        return out;
    }
    let (al, be) = (&basis[0], &basis[1]);
    let (qab, qaa, qbb) = (qform(al, be, p), qform(al, al, p), qform(be, be, p));
    let disc = sb(md(qab, qab, p), md(qaa, qbb, p), p);
    out.disc = disc;
    out.disc_chi = chi(disc, p);
    if out.disc_chi < 0 {
        return out;
    }
    // Roots (lambda : mu) of mu^2 Q(alpha) + lambda mu q(alpha, beta) + lambda^2 Q(beta),
    // Q(x) = q(x, x) / 2; gamma = mu alpha + lambda beta.
    let comb = |mu: u64, la: u64| -> Vec<u64> {
        (0..6)
            .map(|i| ad(md(mu, al[i], p), md(la, be[i], p), p))
            .collect()
    };
    let mut gammas: Vec<Vec<u64>> = Vec::new();
    if qbb != 0 {
        let r = sqrt_mod(disc, p).expect("square discriminant");
        let iq = inv_mod(qbb, p);
        let l1 = md(sb(r, qab, p), iq, p);
        gammas.push(comb(1, l1));
        if r != 0 {
            let l2 = md(sb(ng(r, p), qab, p), iq, p);
            gammas.push(comb(1, l2));
        }
    } else if qab != 0 {
        gammas.push(be.clone());
        let l = ng(md(qaa, inv_mod(md(2, qab, p), p), p), p);
        gammas.push(comb(1, l));
    } else if qaa != 0 {
        gammas.push(be.clone());
    } else {
        out.pencil_in_grassmannian = true;
        return out;
    }
    for g in &gammas {
        let Some(map) = trig_map_of(g, p) else {
            out.maps.push(MapOutcome::default());
            continue;
        };
        let degree3 = degree_three(&map, p);
        let pairs_ok = sg.pairs.iter().all(|pr| pair_condition(&map, pr));
        let eval_chi = if degree3 {
            eval_alpha_chi(f, &map, p)
        } else {
            0
        };
        // A line meeting the twisted cubic gives a map of degree < 3: no
        // trigonal map, so Proposition 4 is not evaluated for it.
        let (alpha, alpha_chi, lemma4_ok) = match map.smith_coeffs() {
            _ if !degree3 => (None, 0, None),
            Some((n1, n0, d1, d0)) => {
                let fs = f_mod_g(f, n1, n0, d1, d0, p);
                let gs = g_coeffs(n1, n0, d1, d0, p);
                let s = smith_s(&fs, &gs, p);
                match s.last() {
                    Some(&lc) => (Some(lc), chi(lc, p), Some(is_lc_times_square(&s, p))),
                    None => (None, 0, Some(false)),
                }
            }
            None => (None, eval_chi, None),
        };
        out.maps.push(MapOutcome {
            map,
            degree3,
            pairs_ok,
            alpha,
            alpha_chi: if degree3 { alpha_chi } else { 0 },
            eval_chi,
            lemma4_ok,
        });
    }
    out
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct CurveOutcome {
    pub cycle_type: Vec<u32>,
    pub subgroups: Vec<SubgroupOutcome>,
}

impl CurveOutcome {
    pub fn n_tractable(&self) -> usize {
        self.subgroups.len()
    }
    pub fn n_trig(&self) -> usize {
        self.subgroups.iter().filter(|s| s.trig_rational()).count()
    }
    pub fn n_iso(&self, twist: bool) -> usize {
        self.subgroups
            .iter()
            .filter(|s| s.iso_rational(twist))
            .count()
    }
    /// Some rational tractable subgroup passes Propositions 2 and 4.
    pub fn full_success(&self, twist: bool) -> bool {
        self.subgroups.iter().any(|s| s.iso_rational(twist))
    }
    pub fn diag(&self) -> Diag {
        let mut d = Diag::default();
        for s in &self.subgroups {
            d.m_rank_deficient += (s.m_rank != 4) as u64;
            d.disc_zero += (s.m_rank == 4 && s.disc_chi == 0) as u64;
            d.pencil_in_grassmannian += s.pencil_in_grassmannian as u64;
            d.lambda_disagree += s.lambda_disagree() as u64;
            d.construction_failed += s.construction_failed() as u64;
            for m in &s.maps {
                d.maps += 1;
                d.non_normal_form += (!m.map.normal_form) as u64;
                d.degree_drop += (!m.degree3) as u64;
                d.pair_condition_fail += (m.degree3 && !m.pairs_ok) as u64;
                d.eval_mismatch += (m.map.normal_form
                    && m.degree3
                    && m.eval_chi != 0
                    && m.alpha_chi != m.eval_chi) as u64;
                d.lemma4_fail += (m.lemma4_ok == Some(false)) as u64;
            }
        }
        d
    }
}

/// Diagnostics of the instrument; every "fail" or "mismatch" must stay 0.
#[derive(Clone, Copy, Debug, Default, Serialize, PartialEq, Eq)]
pub struct Diag {
    pub maps: u64,
    pub m_rank_deficient: u64,
    pub disc_zero: u64,
    pub pencil_in_grassmannian: u64,
    pub non_normal_form: u64,
    pub degree_drop: u64,
    pub pair_condition_fail: u64,
    pub eval_mismatch: u64,
    pub lemma4_fail: u64,
    pub lambda_disagree: u64,
    pub construction_failed: u64,
    pub x_agree: u64,
    pub x_disagree: u64,
    pub x_not_applicable: u64,
    pub x_count_error: u64,
}

impl Diag {
    pub fn merge(&mut self, o: &Diag) {
        self.maps += o.maps;
        self.m_rank_deficient += o.m_rank_deficient;
        self.disc_zero += o.disc_zero;
        self.pencil_in_grassmannian += o.pencil_in_grassmannian;
        self.non_normal_form += o.non_normal_form;
        self.degree_drop += o.degree_drop;
        self.pair_condition_fail += o.pair_condition_fail;
        self.eval_mismatch += o.eval_mismatch;
        self.lemma4_fail += o.lemma4_fail;
        self.lambda_disagree += o.lambda_disagree;
        self.construction_failed += o.construction_failed;
        self.x_agree += o.x_agree;
        self.x_disagree += o.x_disagree;
        self.x_not_applicable += o.x_not_applicable;
        self.x_count_error += o.x_count_error;
    }
}

/// Factor `F`, enumerate `S(H)`, run Propositions 2 and 4 on each subgroup.
/// `None` if `F` is not squarefree.
pub fn analyze_curve(f: &[u64], p: u64, rng: &mut SplitMix64) -> Option<CurveOutcome> {
    assert!(
        p % 2 == 1 && (7..(1 << 31)).contains(&p),
        "odd prime 7 <= p < 2^31"
    );
    let w = Weierstrass::new(f, p, rng)?;
    let subgroups = tractable_subgroups(&w, rng)
        .iter()
        .map(|sg| analyze_subgroup(&w.f, sg, p))
        .collect();
    Some(CurveOutcome {
        cycle_type: w.cycle_type(),
        subgroups,
    })
}

/// Point-count check of Proposition 4 on every valid normal-form map:
/// `#X = #H` must go with a square `lc(s)`, `#X = 2p + 2 - #H` with a
/// non-square.  `twist` flips the square class (the model is then `c F`).
pub fn x_check(f: &[u64], p: u64, out: &CurveOutcome, twist: bool) -> Diag {
    let mut d = Diag::default();
    let h = h_point_count(f, p) as i64;
    let a1 = h - p as i64 - 1;
    for s in &out.subgroups {
        if !s.trig_rational() {
            continue;
        }
        for m in s.maps.iter().filter(|m| m.valid()) {
            let want = if twist { -m.alpha_chi } else { m.alpha_chi };
            match (a1 != 0).then(|| x_point_count(f, &m.map, p)).flatten() {
                None => d.x_not_applicable += 1,
                Some(x) => {
                    let x = x as i64;
                    let says = if x == h {
                        1
                    } else if x == 2 * (p as i64 + 1) - h {
                        -1
                    } else {
                        0
                    };
                    if says == 0 {
                        d.x_count_error += 1;
                    } else if says == want {
                        d.x_agree += 1;
                    } else {
                        d.x_disagree += 1;
                    }
                }
            }
        }
    }
    d
}

// ---------------------------------------------------------------------------
// Clopper-Pearson intervals (regularised incomplete beta, bisection).

fn ln_gamma(x: f64) -> f64 {
    // Lanczos, g = 7, n = 9.
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        let pi = std::f64::consts::PI;
        return (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + 7.5;
    for (i, &c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

fn beta_cf(a: f64, b: f64, x: f64) -> f64 {
    let tiny = 1e-300;
    let (qab, qap, qam) = (a + b, a + 1.0, a - 1.0);
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < tiny {
        d = tiny;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..100_000 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < tiny {
            d = tiny;
        }
        c = 1.0 + aa / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < tiny {
            d = tiny;
        }
        c = 1.0 + aa / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-15 {
            break;
        }
    }
    h
}

/// Regularised incomplete beta `I_x(a, b)`.
pub fn inc_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let lbt = ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln();
    let bt = lbt.exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        bt * beta_cf(a, b, x) / a
    } else {
        1.0 - bt * beta_cf(b, a, 1.0 - x) / b
    }
}

fn beta_quantile(a: f64, b: f64, q: f64) -> f64 {
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if inc_beta(a, b, mid) < q {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Two-sided Clopper-Pearson interval for `k` successes in `n` trials.
pub fn clopper_pearson(k: u64, n: u64, conf: f64) -> (f64, f64) {
    assert!(n > 0 && k <= n);
    let a = (1.0 - conf) / 2.0;
    let lo = if k == 0 {
        0.0
    } else {
        beta_quantile(k as f64, (n - k + 1) as f64, a)
    };
    let hi = if k == n {
        1.0
    } else {
        beta_quantile((k + 1) as f64, (n - k) as f64, 1.0 - a)
    };
    (lo, hi)
}

// ---------------------------------------------------------------------------
// Predictions.

/// Lemma 1: `#S(H)` from the factor-degree type (independent of the code above).
pub fn lemma1_count(ct: &[u32]) -> u32 {
    let table: &[(&[u32], u32)] = &[
        (&[8], 1),
        (&[2, 6], 1),
        (&[1, 1, 6], 1),
        (&[1, 1, 2, 4], 1),
        (&[2, 2, 4], 3),
        (&[1, 1, 1, 1, 4], 3),
        (&[2, 3, 3], 3),
        (&[1, 1, 3, 3], 3),
        (&[4, 4], 5),
        (&[1, 1, 2, 2, 2], 7),
        (&[1, 1, 1, 1, 2, 2], 9),
        (&[1, 1, 1, 1, 1, 1, 2], 15),
        (&[2, 2, 2, 2], 25),
        (&[1; 8], 105),
    ];
    table.iter().find(|(t, _)| *t == ct).map_or(0, |t| t.1)
}

/// Theorem 2 (eq. (12)) with per-subgroup success `q`: the `S_8`-weighted
/// mean of `1 - (1 - q)^{#S}`; `q = 1/4` gives Smith's 0.1857.
pub fn theorem2(q: f64) -> f64 {
    let rep = chebotarev_report();
    rep.rows
        .iter()
        .map(|r| r.class_size as f64 * (1.0 - (1.0 - q).powi(r.stable_pair_kernels as i32)))
        .sum::<f64>()
        / 40320.0
}

/// Chebotarev density of a rational tractable subgroup (`P(#S(H) > 0)`).
pub fn chebotarev_tractable() -> f64 {
    chebotarev_report().s8_pair_fraction
}

// ---------------------------------------------------------------------------
// Question 1: random curves over a large prime.

/// One random squarefree `F` of degree 8 (all coefficients uniform,
/// leading coefficient nonzero) for curve `index` of a run with `seed`.
pub fn random_curve(p: u64, seed: u64, index: u64) -> (Poly, SplitMix64) {
    let mut rng = SplitMix64(seed ^ index.wrapping_mul(0xD1B5_4A32_D192_ED03));
    rng.next_u64();
    loop {
        let mut f: Poly = (0..8).map(|_| rng.below(p)).collect();
        f.push(1 + rng.below(p - 1));
        let fm = pmonic(&f, p);
        if deg(&pgcd(&fm, &pderiv(&fm, p), p)) == 0 {
            return (f, rng);
        }
    }
}

/// Per-curve row of a random run.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct RandomRow {
    pub index: u64,
    pub cycle_type: String,
    pub n_tractable: u32,
    pub n_trig: u32,
    pub n_iso: u32,
    pub full_success: bool,
    pub diag: Diag,
}

pub const RANDOM_HEADER: &str = "index,cycle_type,n_tractable,n_trig,n_iso,full_success,maps,m_rank_deficient,disc_zero,pencil_in_grassmannian,non_normal_form,degree_drop,pair_condition_fail,eval_mismatch,lemma4_fail,lambda_disagree,construction_failed";

impl RandomRow {
    pub fn csv(&self) -> String {
        let d = &self.diag;
        format!(
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            self.index,
            self.cycle_type,
            self.n_tractable,
            self.n_trig,
            self.n_iso,
            self.full_success as u8,
            d.maps,
            d.m_rank_deficient,
            d.disc_zero,
            d.pencil_in_grassmannian,
            d.non_normal_form,
            d.degree_drop,
            d.pair_condition_fail,
            d.eval_mismatch,
            d.lemma4_fail,
            d.lambda_disagree,
            d.construction_failed
        )
    }

    pub fn parse(line: &str) -> Option<RandomRow> {
        let v: Vec<&str> = line.trim().split(',').collect();
        if v.len() != 17 {
            return None;
        }
        let u = |i: usize| v[i].parse::<u64>().ok();
        Some(RandomRow {
            index: u(0)?,
            cycle_type: v[1].to_string(),
            n_tractable: u(2)? as u32,
            n_trig: u(3)? as u32,
            n_iso: u(4)? as u32,
            full_success: u(5)? == 1,
            diag: Diag {
                maps: u(6)?,
                m_rank_deficient: u(7)?,
                disc_zero: u(8)?,
                pencil_in_grassmannian: u(9)?,
                non_normal_form: u(10)?,
                degree_drop: u(11)?,
                pair_condition_fail: u(12)?,
                eval_mismatch: u(13)?,
                lemma4_fail: u(14)?,
                lambda_disagree: u(15)?,
                construction_failed: u(16)?,
                ..Default::default()
            },
        })
    }
}

pub fn ct_string(ct: &[u32]) -> String {
    ct.iter()
        .map(|d| d.to_string())
        .collect::<Vec<_>>()
        .join("+")
}

pub fn random_row(p: u64, seed: u64, index: u64) -> RandomRow {
    let (f, mut rng) = random_curve(p, seed, index);
    let out = analyze_curve(&f, p, &mut rng).expect("squarefree by construction");
    RandomRow {
        index,
        cycle_type: ct_string(&out.cycle_type),
        n_tractable: out.n_tractable() as u32,
        n_trig: out.n_trig() as u32,
        n_iso: out.n_iso(false) as u32,
        full_success: out.full_success(false),
        diag: out.diag(),
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct FractionRow {
    pub name: String,
    pub unit: String,
    pub k: u64,
    pub n: u64,
    pub value: f64,
    /// Clopper-Pearson at the Bonferroni level `1 - 0.05/4`.
    pub cp_lo: f64,
    pub cp_hi: f64,
    /// Percentile curve-cluster bootstrap at the same level (subgroup units).
    pub boot_lo: Option<f64>,
    pub boot_hi: Option<f64>,
    pub prediction: f64,
    pub prediction_source: String,
    pub smith_example2: f64,
    pub margin: f64,
    pub ci_contains_prediction: bool,
    pub verdict: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct RandomSummary {
    pub p: u64,
    pub seed: u64,
    pub n_curves: u64,
    pub confidence: f64,
    pub bootstrap_resamples: u64,
    pub subgroups: u64,
    pub subgroups_per_curve: f64,
    pub fractions: Vec<FractionRow>,
    /// Protocol Section 5: "Smith §8 reproduced" or the fractions that are not.
    pub overall_verdict: String,
    /// Every instrument invariant of protocol Section 4 is 0.
    pub instrument_valid: bool,
    pub diag: Diag,
    pub cycle_type_counts: BTreeMap<String, u64>,
}

impl Diag {
    /// The instrument invariants of the protocol (Section 4) all hold.
    pub fn instrument_valid(&self) -> bool {
        self.pair_condition_fail == 0
            && self.eval_mismatch == 0
            && self.lemma4_fail == 0
            && self.x_disagree == 0
            && self.x_count_error == 0
    }
}

/// Fixed verdict wording of the protocol's Question 1 decision rule.
pub fn equivalence_verdict(lo: f64, hi: f64, pred: f64, margin: f64) -> String {
    if lo >= pred - margin && hi <= pred + margin {
        format!("equivalent to the prediction within ±{margin}")
    } else if hi < pred - margin || lo > pred + margin {
        format!("differs from the prediction by more than {margin}")
    } else {
        format!("inconclusive at ±{margin}")
    }
}

/// Summary of a random run, computed from its rows only.
pub fn summarize_random(
    rows: &[RandomRow],
    p: u64,
    seed: u64,
    boot: u64,
    boot_seed: u64,
) -> RandomSummary {
    let conf = 1.0 - 0.05 / 4.0;
    let margin = 0.01;
    let n = rows.len() as u64;
    let tract = rows.iter().filter(|r| r.n_tractable > 0).count() as u64;
    let subs: u64 = rows.iter().map(|r| r.n_tractable as u64).sum();
    let trig: u64 = rows.iter().map(|r| r.n_trig as u64).sum();
    let iso: u64 = rows.iter().map(|r| r.n_iso as u64).sum();
    let full = rows.iter().filter(|r| r.full_success).count() as u64;
    let mut diag = Diag::default();
    let mut cts: BTreeMap<String, u64> = BTreeMap::new();
    for r in rows {
        diag.merge(&r.diag);
        *cts.entry(r.cycle_type.clone()).or_default() += 1;
    }
    // Curve-cluster bootstrap for the two subgroup-level ratios.
    let compact: Vec<[u8; 3]> = rows
        .iter()
        .map(|r| [r.n_tractable as u8, r.n_trig as u8, r.n_iso as u8])
        .collect();
    let mut rng = SplitMix64(boot_seed);
    let mut b2 = Vec::with_capacity(boot as usize);
    let mut b3 = Vec::with_capacity(boot as usize);
    for _ in 0..boot {
        let (mut s, mut t, mut i) = (0u64, 0u64, 0u64);
        for _ in 0..n {
            let r = compact[rng.below(n) as usize];
            s += r[0] as u64;
            t += r[1] as u64;
            i += r[2] as u64;
        }
        b2.push(t as f64 / s.max(1) as f64);
        b3.push(i as f64 / t.max(1) as f64);
    }
    let pct = |v: &mut Vec<f64>| -> Option<(f64, f64)> {
        if v.is_empty() {
            return None;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let a = (1.0 - conf) / 2.0;
        let lo = v[((a * v.len() as f64).floor() as usize).min(v.len() - 1)];
        let hi = v[(((1.0 - a) * v.len() as f64).ceil() as usize).min(v.len()) - 1];
        Some((lo, hi))
    };
    let (bt2, bt3) = (pct(&mut b2), pct(&mut b3));
    let mk = |name: &str,
              unit: &str,
              k: u64,
              nn: u64,
              pred: f64,
              src: &str,
              smith: f64,
              bt: Option<(f64, f64)>| {
        let (lo, hi) = clopper_pearson(k, nn, conf);
        FractionRow {
            name: name.into(),
            unit: unit.into(),
            k,
            n: nn,
            value: k as f64 / nn as f64,
            cp_lo: lo,
            cp_hi: hi,
            boot_lo: bt.map(|b| b.0),
            boot_hi: bt.map(|b| b.1),
            prediction: pred,
            prediction_source: src.into(),
            smith_example2: smith,
            margin,
            ci_contains_prediction: lo <= pred && pred <= hi,
            verdict: equivalence_verdict(lo, hi, pred, margin),
        }
    };
    let fractions = vec![
        mk(
            "curves with a rational tractable subgroup",
            "curves",
            tract,
            n,
            chebotarev_tractable(),
            "Chebotarev density of #S(H) > 0 over S_8 (Lemma 1)",
            502_005.0 / 1e6,
            None,
        ),
        mk(
            "tractable subgroups with a rational trigonal map (Prop. 2)",
            "subgroups",
            trig,
            subs,
            0.5,
            "Hypothesis 2",
            501_629.0 / 1_002_244.0,
            bt2,
        ),
        mk(
            "of those, subgroups with a rational isogeny (Prop. 4)",
            "subgroups",
            iso,
            trig,
            0.5,
            "Hypothesis 3",
            250_560.0 / 501_629.0,
            bt3,
        ),
        mk(
            "curves with full_success",
            "curves",
            full,
            n,
            theorem2(0.25),
            "Theorem 2, eq. (12)",
            185_814.0 / 1e6,
            None,
        ),
    ];
    let not_eq: Vec<String> = fractions
        .iter()
        .filter(|f| !f.verdict.starts_with("equivalent"))
        .map(|f| format!("{} ({})", f.name, f.verdict))
        .collect();
    let overall_verdict = if !diag.instrument_valid() {
        "instrument invariant failed — run invalid".to_string()
    } else if not_eq.is_empty() {
        "Smith §8 reproduced (all four fractions equivalent within ±0.01)".to_string()
    } else {
        format!("not reproduced: {}", not_eq.join("; "))
    };
    RandomSummary {
        p,
        seed,
        n_curves: n,
        confidence: conf,
        bootstrap_resamples: boot,
        subgroups: subs,
        subgroups_per_curve: subs as f64 / n as f64,
        fractions,
        overall_verdict,
        instrument_valid: diag.instrument_valid(),
        diag,
        cycle_type_counts: cts,
    }
}

// ---------------------------------------------------------------------------
// Question 2: exhaustive census over a small prime, per isogeny class.

/// Strata of the independence null: the value of `#S(H)` (Lemma 1).
pub const STRATA: [u32; 9] = [0, 1, 3, 5, 7, 9, 15, 25, 105];

fn stratum_index(ns: u32) -> usize {
    STRATA
        .iter()
        .position(|&s| s == ns)
        .expect("#S(H) is a Lemma 1 value")
}

#[derive(Clone, Copy, Debug)]
struct OcticRec {
    key: [i16; 3],
    stratum: u8,
    n_trig: u8,
    full_curve: bool,
    full_twist: bool,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct TypeRow {
    pub cycle_type: String,
    pub curves: u64,
    pub subgroups: u64,
    pub trig_subgroups: u64,
    pub iso_subgroups: u64,
    pub full_success_curves: u64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct ClassRow {
    pub a1: i64,
    pub a2: i64,
    pub a3: i64,
    pub curves: u64,
    pub tractable_curves: u64,
    pub trig_curves: u64,
    pub full_success_curves: u64,
    /// Curves per `#S(H)` stratum (`STRATA` order).
    pub by_stratum: [u64; 9],
}

pub const CLASS_HEADER: &str = "a1,a2,a3,curves,tractable_curves,trig_curves,full_success_curves,s0,s1,s3,s5,s7,s9,s15,s25,s105";
pub const TYPE_HEADER: &str =
    "cycle_type,curves,subgroups,trig_subgroups,iso_subgroups,full_success_curves";

impl ClassRow {
    pub fn csv(&self) -> String {
        let s: Vec<String> = self.by_stratum.iter().map(|x| x.to_string()).collect();
        format!(
            "{},{},{},{},{},{},{},{}",
            self.a1,
            self.a2,
            self.a3,
            self.curves,
            self.tractable_curves,
            self.trig_curves,
            self.full_success_curves,
            s.join(",")
        )
    }

    pub fn parse(line: &str) -> Option<ClassRow> {
        let v: Vec<&str> = line.trim().split(',').collect();
        if v.len() != 16 {
            return None;
        }
        let i = |k: usize| v[k].parse::<i64>().ok();
        let mut by = [0u64; 9];
        for (j, b) in by.iter_mut().enumerate() {
            *b = v[7 + j].parse().ok()?;
        }
        Some(ClassRow {
            a1: i(0)?,
            a2: i(1)?,
            a3: i(2)?,
            curves: i(3)? as u64,
            tractable_curves: i(4)? as u64,
            trig_curves: i(5)? as u64,
            full_success_curves: i(6)? as u64,
            by_stratum: by,
        })
    }
}

impl TypeRow {
    pub fn csv(&self) -> String {
        format!(
            "{},{},{},{},{},{}",
            self.cycle_type,
            self.curves,
            self.subgroups,
            self.trig_subgroups,
            self.iso_subgroups,
            self.full_success_curves
        )
    }

    pub fn parse(line: &str) -> Option<TypeRow> {
        let v: Vec<&str> = line.trim().split(',').collect();
        if v.len() != 6 {
            return None;
        }
        let u = |k: usize| v[k].parse::<u64>().ok();
        Some(TypeRow {
            cycle_type: v[0].to_string(),
            curves: u(1)?,
            subgroups: u(2)?,
            trig_subgroups: u(3)?,
            iso_subgroups: u(4)?,
            full_success_curves: u(5)?,
        })
    }
}

/// Raw output of the census: per-class rows, per-type rows, diagnostics and
/// the per-stratum success counts the permutation null draws from.
#[derive(Clone, Debug, Default, Serialize)]
pub struct CensusData {
    pub p: u64,
    pub octics_enumerated: u64,
    pub classes: Vec<ClassRow>,
    pub types: Vec<TypeRow>,
    /// Successes (`full_success` curves) per `#S(H)` stratum.
    pub stratum_successes: [u64; 9],
    pub diag: Diag,
    pub twist_constant: u64,
}

#[derive(Default)]
struct CensusAcc {
    classes: HashMap<[i16; 3], ClassRow>,
    types: BTreeMap<String, TypeRow>,
    diag: Diag,
    strat: [u64; 9],
}

impl CensusAcc {
    fn add(&mut self, rec: &OcticRec, ct: String, out: &CurveOutcome, diag: &Diag) {
        self.diag.merge(diag);
        let t = self.types.entry(ct.clone()).or_insert_with(|| TypeRow {
            cycle_type: ct,
            ..Default::default()
        });
        t.curves += 2;
        t.subgroups += 2 * out.n_tractable() as u64;
        t.trig_subgroups += 2 * out.n_trig() as u64;
        t.iso_subgroups += (out.n_iso(false) + out.n_iso(true)) as u64;
        t.full_success_curves += rec.full_curve as u64 + rec.full_twist as u64;
        let [a1, a2, a3] = rec.key;
        for (key, full) in [
            ([a1, a2, a3], rec.full_curve),
            ([-a1, a2, -a3], rec.full_twist),
        ] {
            let e = self.classes.entry(key).or_insert_with(|| ClassRow {
                a1: key[0] as i64,
                a2: key[1] as i64,
                a3: key[2] as i64,
                ..Default::default()
            });
            e.curves += 1;
            e.tractable_curves += (rec.stratum != 0) as u64;
            e.trig_curves += (rec.n_trig > 0) as u64;
            e.full_success_curves += full as u64;
            e.by_stratum[rec.stratum as usize] += 1;
            self.strat[rec.stratum as usize] += full as u64;
        }
    }

    fn merge(mut self, o: CensusAcc) -> CensusAcc {
        self.diag.merge(&o.diag);
        for (k, t) in o.types {
            let e = self.types.entry(k).or_insert_with(|| TypeRow {
                cycle_type: t.cycle_type.clone(),
                ..Default::default()
            });
            e.curves += t.curves;
            e.subgroups += t.subgroups;
            e.trig_subgroups += t.trig_subgroups;
            e.iso_subgroups += t.iso_subgroups;
            e.full_success_curves += t.full_success_curves;
        }
        for (k, c) in o.classes {
            let e = self.classes.entry(k).or_insert_with(|| ClassRow {
                a1: c.a1,
                a2: c.a2,
                a3: c.a3,
                ..Default::default()
            });
            e.curves += c.curves;
            e.tractable_curves += c.tractable_curves;
            e.trig_curves += c.trig_curves;
            e.full_success_curves += c.full_success_curves;
            for s in 0..9 {
                e.by_stratum[s] += c.by_stratum[s];
            }
        }
        for s in 0..9 {
            self.strat[s] += o.strat[s];
        }
        self
    }
}

/// Exhaustive census over the monic octics with zero `x^7` coefficient
/// (`p^7` of them, or the first `limit` codes for a dry run).  Every
/// squarefree one stands for two curves, `y^2 = f` and its twist `y^2 = c f`
/// (`c` the least non-residue), as in `smith_kernel_coverage::prime_report`.
/// Classes are `L`-polynomials from point counts over `F_p, F_{p^2}, F_{p^3}`.
/// The point-count check of Proposition 4 runs on the codes divisible by
/// `xcheck_every` (0: never).
pub fn census(p: u64, xcheck_every: u64, limit: Option<u64>) -> CensusData {
    use rayon::prelude::*;
    assert!(p % 2 == 1 && (7..=37).contains(&p));
    let pu = p as u32;
    let fields = [Fpk::new(pu, 1), Fpk::new(pu, 2), Fpk::new(pu, 3)];
    let c_tw = least_nonsquare(p);
    let total = limit.unwrap_or(p.pow(7)).min(p.pow(7));
    let acc = (0..total)
        .into_par_iter()
        .fold(CensusAcc::default, |mut acc, code| {
            let mut c = [0u32; 8];
            let mut v = code;
            for cd in c.iter_mut().take(7) {
                *cd = (v % p) as u32;
                v /= p;
            }
            let mut f: Poly = c.iter().map(|&x| x as u64).collect();
            f.push(1);
            let mut rng = SplitMix64(code);
            let Some(out) = analyze_curve(&f, p, &mut rng) else {
                return acc;
            };
            let s: Vec<i64> = fields.iter().map(|fk| fk.char_sum(&c)).collect();
            let pk = [-s[0] - 1, -s[1] - 1, -s[2] - 1];
            let e1 = pk[0];
            let e2 = (e1 * pk[0] - pk[1]) / 2;
            let e3 = (e2 * pk[0] - e1 * pk[1] + pk[2]) / 3;
            let key = [(-e1) as i16, e2 as i16, (-e3) as i16];
            let mut diag = out.diag();
            if xcheck_every > 0 && code % xcheck_every == 0 {
                diag.merge(&x_check(&f, p, &out, false));
                let ft = pscale(&f, c_tw, p);
                diag.merge(&x_check(&ft, p, &out, true));
            }
            let ns = out.n_tractable() as u32;
            let rec = OcticRec {
                key,
                stratum: stratum_index(ns) as u8,
                n_trig: out.n_trig() as u8,
                full_curve: out.full_success(false),
                full_twist: out.full_success(true),
            };
            acc.add(&rec, ct_string(&out.cycle_type), &out, &diag);
            acc
        })
        .reduce(CensusAcc::default, CensusAcc::merge);
    let mut classes: Vec<ClassRow> = acc.classes.into_values().collect();
    classes.sort_by_key(|c| (c.a1, c.a2, c.a3));
    CensusData {
        p,
        octics_enumerated: total,
        classes,
        types: acc.types.into_values().collect(),
        stratum_successes: acc.strat,
        diag: acc.diag,
        twist_constant: c_tw,
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct NullSummary {
    pub permutations: u64,
    pub seed: u64,
    pub class_fraction_median: f64,
    pub class_fraction_lo: f64,
    pub class_fraction_hi: f64,
    pub r_median: f64,
    pub r_lo: f64,
    pub r_hi: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct CensusSummary {
    pub p: u64,
    pub octics_enumerated: u64,
    pub curves: u64,
    pub classes: u64,
    /// Direct fraction of curves with a rational tractable subgroup.
    pub direct_tractable: f64,
    /// U: curves whose class has a member with a rational tractable subgroup.
    pub class_tractable: f64,
    /// D: direct fraction of curves with `full_success`.
    pub direct_full: f64,
    /// C: curves whose class has a member with `full_success`.
    pub class_full: f64,
    /// Classes with at least one tractable member, and how many of them
    /// have a `full_success` member.
    pub classes_with_tractable: u64,
    pub classes_with_full: u64,
    /// `r = (C - D) / (U - D)`.
    pub r: f64,
    /// The band of `r` alone (Section 5 thresholds).
    pub r_band: String,
    /// The verdict after the control gates of protocol Sections 4 and 5.
    pub verdict: String,
    /// Every instrument invariant of protocol Section 4 is 0.
    pub instrument_valid: bool,
    /// Curves, classes and tractable curves equal the frozen 2026-10-13
    /// counts (`None` for a prime that study did not run).
    pub control_frozen_counts_match: Option<bool>,
    /// Control: the class closure of "has a rational tractable subgroup"
    /// equals its direct fraction (a property fixed by `chi_pi mod 2`).
    pub control_tractable_closure_equal: bool,
    pub null: NullSummary,
    pub p2_measured: f64,
    pub p3_measured: f64,
    pub p2_heuristic: f64,
    pub smallp_direct_prediction_heuristic: f64,
    pub smallp_direct_prediction_measured: f64,
    pub asymptotic_direct_prediction: f64,
    pub asymptotic_ceiling: f64,
    pub diag: Diag,
}

/// Fixed verdict wording of the protocol's Question 2 decision rule.
pub fn class_verdict(r: f64) -> &'static str {
    if r >= 0.75 {
        "near the ceiling"
    } else if r <= 0.25 {
        "near the direct fraction"
    } else {
        "intermediate"
    }
}

/// `(curves, classes, tractable curves)` frozen by the 2026-10-13 study
/// (`research/weak_isogeny_classes_20261013/data/smith_p<p>.json`).
pub fn frozen_20261013(p: u64) -> Option<(u64, u64, u64)> {
    match p {
        7 => Some((1_411_788, 5231, 691_308)),
        11 => Some((35_431_220, 21_701, 17_513_540)),
        _ => None,
    }
}

/// Overall Question 2 verdict: the common per-prime verdict, else
/// "prime-dependent — inconclusive".
pub fn overall_q2(per_prime: &[CensusSummary]) -> String {
    let v: Vec<&str> = per_prime.iter().map(|s| s.verdict.as_str()).collect();
    match v.first() {
        Some(first) if v.iter().all(|x| x == first) => first.to_string(),
        Some(_) => "prime-dependent — inconclusive".to_string(),
        None => "no census".to_string(),
    }
}

/// Summary of a census, from its class and type rows only.
pub fn summarize_census(
    p: u64,
    octics: u64,
    classes: &[ClassRow],
    types: &[TypeRow],
    stratum_successes: &[u64; 9],
    diag: Diag,
    perms: u64,
    seed: u64,
) -> CensusSummary {
    let curves: u64 = classes.iter().map(|c| c.curves).sum();
    let tf = curves as f64;
    let direct_tract: u64 = classes.iter().map(|c| c.tractable_curves).sum();
    let class_tract: u64 = classes
        .iter()
        .filter(|c| c.tractable_curves > 0)
        .map(|c| c.curves)
        .sum();
    let direct_full: u64 = classes.iter().map(|c| c.full_success_curves).sum();
    let class_full: u64 = classes
        .iter()
        .filter(|c| c.full_success_curves > 0)
        .map(|c| c.curves)
        .sum();
    let (u, d, c) = (
        class_tract as f64 / tf,
        direct_full as f64 / tf,
        class_full as f64 / tf,
    );
    let r = (c - d) / (u - d);
    // Independence null: within each #S(H) stratum, scatter that stratum's
    // successes uniformly over its curves, then recompute C.
    let n_cls = classes.len();
    let mut slots: Vec<Vec<u32>> = vec![Vec::new(); 9];
    for (ci, cl) in classes.iter().enumerate() {
        for s in 0..9 {
            for _ in 0..cl.by_stratum[s] {
                slots[s].push(ci as u32);
            }
        }
    }
    let mut rng = SplitMix64(seed);
    let mut cs = Vec::with_capacity(perms as usize);
    for _ in 0..perms {
        let mut hit = vec![false; n_cls];
        for s in 0..9 {
            let k = stratum_successes[s] as usize;
            let v = &mut slots[s];
            let n = v.len();
            for i in 0..k.min(n) {
                let j = i + rng.below((n - i) as u64) as usize;
                v.swap(i, j);
                hit[v[i] as usize] = true;
            }
        }
        let cf: u64 = classes
            .iter()
            .zip(&hit)
            .filter(|(_, &h)| h)
            .map(|(c, _)| c.curves)
            .sum();
        cs.push(cf as f64 / tf);
    }
    cs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |x: f64| -> f64 {
        if cs.is_empty() {
            return f64::NAN;
        }
        cs[((x * (cs.len() - 1) as f64).round() as usize).min(cs.len() - 1)]
    };
    let rr = |cv: f64| (cv - d) / (u - d);
    let null = NullSummary {
        permutations: perms,
        seed,
        class_fraction_median: q(0.5),
        class_fraction_lo: q(0.025),
        class_fraction_hi: q(0.975),
        r_median: rr(q(0.5)),
        r_lo: rr(q(0.025)),
        r_hi: rr(q(0.975)),
    };
    let frozen_ok = frozen_20261013(p)
        .map(|(cv, cl, tr)| cv == curves && cl == classes.len() as u64 && tr == direct_tract);
    let verdict = if !diag.instrument_valid() {
        "instrument invariant failed — run invalid".to_string()
    } else if direct_tract != class_tract || frozen_ok == Some(false) || null.r_median < 0.75 {
        "control failed — inconclusive".to_string()
    } else {
        class_verdict(r).to_string()
    };
    // Smith's model at this p: exact type counts, per-subgroup success q.
    let subs: u64 = types.iter().map(|t| t.subgroups).sum();
    let trig: u64 = types.iter().map(|t| t.trig_subgroups).sum();
    let iso: u64 = types.iter().map(|t| t.iso_subgroups).sum();
    let p2m = trig as f64 / subs as f64;
    let p3m = iso as f64 / trig as f64;
    let p2h = (p as f64 + 1.0) / (2.0 * p as f64);
    let model = |qq: f64| -> f64 {
        types
            .iter()
            .map(|t| {
                let ns = t.subgroups.checked_div(t.curves).unwrap_or(0);
                t.curves as f64 * (1.0 - (1.0 - qq).powi(ns as i32))
            })
            .sum::<f64>()
            / tf
    };
    CensusSummary {
        p,
        octics_enumerated: octics,
        curves,
        classes: classes.len() as u64,
        direct_tractable: direct_tract as f64 / tf,
        class_tractable: u,
        direct_full: d,
        class_full: c,
        classes_with_tractable: classes.iter().filter(|c| c.tractable_curves > 0).count() as u64,
        classes_with_full: classes.iter().filter(|c| c.full_success_curves > 0).count() as u64,
        r,
        r_band: class_verdict(r).to_string(),
        verdict,
        instrument_valid: diag.instrument_valid(),
        control_frozen_counts_match: frozen_ok,
        control_tractable_closure_equal: direct_tract == class_tract,
        null,
        p2_measured: p2m,
        p3_measured: p3m,
        p2_heuristic: p2h,
        smallp_direct_prediction_heuristic: model(p2h * 0.5),
        smallp_direct_prediction_measured: model(p2m * p3m),
        asymptotic_direct_prediction: theorem2(0.25),
        asymptotic_ceiling: chebotarev_tractable(),
        diag,
    }
}

// ---------------------------------------------------------------------------
// Smith's Example 1 (positive control).

/// `H: y^2 = x^7 + 28x^6 + 15x^5 + 20x^4 + 33x^3 + 12x^2 + 29x + 2` over `F_37`.
pub const EXAMPLE1_P: u64 = 37;
pub const EXAMPLE1_F: [u64; 8] = [2, 29, 12, 33, 20, 15, 28, 1];
/// `N = x^3 + 16x + 22`, `D = x^2 + 32x + 18` as `(n1, n0, d1, d0)`.
pub const EXAMPLE1_ND: (u64, u64, u64, u64) = (16, 22, 32, 18);

#[derive(Clone, Debug, Serialize)]
pub struct Example1Report {
    pub p: u64,
    pub f: Poly,
    pub cycle_type: Vec<u32>,
    pub tractable_subgroups: usize,
    pub trig_rational: bool,
    pub maps: Vec<MapOutcome>,
    pub smith_map_found: bool,
    pub iso_rational: bool,
    pub h_points: u64,
    pub h_points_expected: u64,
    pub delta1_squared_matches_64s: bool,
    pub delta0_matches: bool,
    pub delta2_matches: bool,
    pub delta4_matches: bool,
    pub f_coeffs_match: bool,
    pub reversed_model_iso_rational: bool,
    pub reversed_model_x_check: Diag,
    pub passed: bool,
}

fn poly_from_desc(desc: &[u64]) -> Poly {
    // Coefficients written highest degree first, as in the paper.
    let mut v: Poly = desc.iter().rev().copied().collect();
    trim(&mut v);
    v
}

/// Reproduce Smith's Example 1 through the same code path as every run.
pub fn example1() -> Example1Report {
    let p = EXAMPLE1_P;
    let f = EXAMPLE1_F.to_vec();
    let mut rng = SplitMix64(1);
    let out = analyze_curve(&f, p, &mut rng).expect("squarefree");
    let sg = out.subgroups.first();
    let maps = sg.map(|s| s.maps.clone()).unwrap_or_default();
    let (n1, n0, d1, d0) = EXAMPLE1_ND;
    let smith_map_found = maps
        .iter()
        .any(|m| m.map.smith_coeffs() == Some((n1, n0, d1, d0)));
    let fs = f_mod_g(&f, n1, n0, d1, d0, p);
    let gs = g_coeffs(n1, n0, d1, d0, p);
    let s = smith_s(&fs, &gs, p);
    let [dl0, dl2, dl4] = smith_deltas(&fs, &gs, p);
    let delta1 = poly_from_desc(&[35, 8, 33, 3]);
    let delta1_ok = pmul(&delta1, &delta1, p) == pscale(&s, 64, p);
    let delta0_ok = dl0 == poly_from_desc(&[27, 20, 33, 6, 16, 8, 9, 2, 31, 15, 16]);
    let delta2_ok = dl2 == poly_from_desc(&[20, 18, 29, 14, 6, 20, 12, 16]);
    let delta4_ok = dl4 == poly_from_desc(&[27, 36, 13, 21, 0]);
    // The equations of X: constant terms -f0, -f1, -f2 (in b's: 2b01, 2b02 + b11).
    let mf0 = poly_from_desc(&[19, 10, 12, 7, 1, 30]);
    let mf1 = poly_from_desc(&[5, 26, 15, 23, 19, 17]);
    let mf2 = poly_from_desc(&[36, 29, 7, 13, 21, 18]);
    let f_ok =
        psub(&[], &fs[0], p) == mf0 && psub(&[], &fs[1], p) == mf1 && psub(&[], &fs[2], p) == mf2;
    // Z(H; T) numerator 1 + 4T + ...: #H(F_37) = 37 + 1 + 4.
    let h_points = h_point_count(&f, p);
    // The same curve with x -> 1/x (a degree-8 model): same verdict, and the
    // point count of X applies.
    let mut fr: Poly = vec![0u64; 9];
    for (i, &c) in f.iter().enumerate() {
        fr[8 - i] = c;
    }
    let mut rng2 = SplitMix64(2);
    let outr = analyze_curve(&fr, p, &mut rng2).expect("squarefree");
    let xr = x_check(&fr, p, &outr, false);
    let iso = out.full_success(false);
    let passed = out.n_tractable() == 1
        && out.n_trig() == 1
        && smith_map_found
        && iso
        && h_points == 42
        && delta1_ok
        && delta0_ok
        && delta2_ok
        && delta4_ok
        && f_ok
        && outr.full_success(false)
        && xr.x_disagree == 0
        && xr.x_count_error == 0;
    Example1Report {
        p,
        f,
        cycle_type: out.cycle_type.clone(),
        tractable_subgroups: out.n_tractable(),
        trig_rational: out.n_trig() == 1,
        maps,
        smith_map_found,
        iso_rational: iso,
        h_points,
        h_points_expected: 42,
        delta1_squared_matches_64s: delta1_ok,
        delta0_matches: delta0_ok,
        delta2_matches: delta2_ok,
        delta4_matches: delta4_ok,
        f_coeffs_match: f_ok,
        reversed_model_iso_rational: outr.full_success(false),
        reversed_model_x_check: xr,
        passed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::smith_kernel_coverage::{
        cycle_type as ct_reference, lagrangians, perm_of_type, stable_kernels,
    };

    fn rand_octic(p: u64, rng: &mut SplitMix64) -> Poly {
        let mut f: Poly = (0..8).map(|_| rng.below(p)).collect();
        f.push(1 + rng.below(p - 1));
        f
    }

    fn shift(f: &[u64], b: u64, p: u64) -> Poly {
        // f(x + b) by Horner in F_p[x].
        let mut r: Poly = Vec::new();
        for &c in f.iter().rev() {
            r = padd(&pmul(&r, &[b, 1], p), &[c], p);
        }
        r
    }

    #[test]
    fn field_helpers() {
        for p in [7u64, 11, 13, 37, 2_147_483_647] {
            for a in [2u64, 3, 5, 6, 10] {
                if chi(a, p) == 1 {
                    let r = sqrt_mod(a, p).unwrap();
                    assert_eq!(md(r, r, p), a % p);
                }
            }
        }
        let k = Ext::new(&[3, 0, 0, 1], 7); // x^3 + 3, irreducible over F_7
        let a = [1, 2, 3, 0, 0, 0, 0, 0];
        assert_eq!(k.mul(&a, &k.inv(&a)), k.one());
        assert_eq!(k.trace(&k.one()), 3);
        // Tr(z) = -(coefficient of z^2) = 0.
        assert_eq!(k.trace(&k.z()), 0);
        // Frobenius has order 3.
        let z = k.z();
        assert_eq!(k.frob(&k.frob(&k.frob(&z))), z);
        assert_ne!(k.frob(&z), z);
    }

    #[test]
    fn clopper_pearson_reference_values() {
        let (lo, hi) = clopper_pearson(5, 10, 0.95);
        assert!((lo - 0.187_086).abs() < 1e-5 && (hi - 0.812_914).abs() < 1e-5);
        let (lo, hi) = clopper_pearson(0, 10, 0.95);
        assert_eq!(lo, 0.0);
        assert!((hi - (1.0 - 0.025f64.powf(0.1))).abs() < 1e-9);
        let (lo, hi) = clopper_pearson(18_574, 100_000, 0.95);
        assert!(lo < 0.18574 && hi > 0.18574 && hi - lo < 0.0055 && hi - lo > 0.0045);
    }

    #[test]
    fn factorisation_matches_reference() {
        let mut rng = SplitMix64(11);
        for p in [7u64, 11, 13] {
            for _ in 0..400 {
                let f = rand_octic(p, &mut rng);
                let fu: Vec<u32> = f.iter().map(|&x| x as u32).collect();
                let reference = ct_reference(
                    &pmonic(&f, p).iter().map(|&x| x as u32).collect::<Vec<_>>(),
                    p as u32,
                );
                let mine = factor_squarefree(&f, p, &mut rng);
                match (reference, mine) {
                    (None, None) => {}
                    (Some(ct), Some(fs)) => {
                        let mut d: Vec<u32> = fs.iter().map(|g| (g.len() - 1) as u32).collect();
                        d.sort_unstable();
                        assert_eq!(d, ct, "{fu:?}");
                        let prod = fs.iter().fold(vec![1u64], |a, g| pmul(&a, g, p));
                        assert_eq!(prod, pmonic(&f, p));
                    }
                    (a, b) => panic!("squarefree disagreement {a:?} {b:?}"),
                }
            }
        }
    }

    #[test]
    fn lemma1_counts_and_descent_rank() {
        let lags = lagrangians();
        let mut rng = SplitMix64(12);
        let mut seen = std::collections::BTreeSet::new();
        for p in [7u64, 13, 17] {
            for _ in 0..1500 {
                let f = rand_octic(p, &mut rng);
                let Some(w) = Weierstrass::new(&f, p, &mut rng) else {
                    continue;
                };
                let ct = w.cycle_type();
                let sgs = tractable_subgroups(&w, &mut rng);
                let (pair, _) = stable_kernels(&perm_of_type(&ct), &lags);
                assert_eq!(sgs.len() as u32, pair, "{ct:?}");
                assert_eq!(sgs.len() as u32, lemma1_count(&ct), "{ct:?}");
                for sg in &sgs {
                    // Every pair is a pair of distinct Weierstrass points.
                    for pr in &sg.pairs {
                        let k = &pr.k;
                        for pt in [&pr.a, &pr.b] {
                            let x = pt.u;
                            let mut acc = k.zero();
                            for &c in w.f.iter().rev() {
                                acc = k.add(&k.mul(&acc, &x), &k.c(c));
                            }
                            assert!(k.is_zero(&acc));
                        }
                        assert!(!k.is_zero(&k.sub(&pr.a.u, &pr.b.u)));
                    }
                    let (rank, _) = nullspace(&m_rows_fp(sg), 6, p);
                    assert!(rank <= 4);
                }
                if !sgs.is_empty() {
                    seen.insert(ct);
                }
            }
        }
        assert!(seen.len() >= 10, "types covered: {}", seen.len());
    }

    #[test]
    fn example1_positive_control() {
        let r = example1();
        assert_eq!(r.cycle_type, vec![1, 1, 6]);
        assert!(r.smith_map_found, "{:?}", r.maps);
        assert!(r.delta1_squared_matches_64s);
        assert!(r.delta0_matches && r.delta2_matches && r.delta4_matches);
        assert!(r.f_coeffs_match);
        assert_eq!(r.h_points, 42);
        assert!(r.iso_rational);
        assert!(r.reversed_model_iso_rational);
        assert!(r.passed, "{r:?}");
    }

    #[test]
    fn eq7_is_the_norm() {
        let mut rng = SplitMix64(13);
        let p = 101;
        let mut checked = 0;
        for _ in 0..300 {
            let f = rand_octic(p, &mut rng);
            let Some(out) = analyze_curve(&f, p, &mut rng) else {
                continue;
            };
            for m in out.subgroups.iter().flat_map(|s| s.maps.iter()) {
                let Some((n1, n0, d1, d0)) = m.map.smith_coeffs() else {
                    continue;
                };
                let fs = f_mod_g(&f, n1, n0, d1, d0, p);
                let gs = g_coeffs(n1, n0, d1, d0, p);
                let s = smith_s(&fs, &gs, p);
                for t0 in 0..20 {
                    let g = vec![
                        sb(n0, md(d0, t0, p), p),
                        sb(n1, md(d1, t0, p), p),
                        ng(t0, p),
                        1,
                    ];
                    let norm = norm_mod_monic(&prem(&f, &g, p), &g, p);
                    assert_eq!(peval(&s, t0, p), norm);
                    // f0 + f1 x + f2 x^2 = F mod G at t0.
                    let at: Poly = fs.iter().map(|q| peval(q, t0, p)).collect();
                    let mut at = at;
                    trim(&mut at);
                    assert_eq!(at, prem(&f, &g, p));
                }
                checked += 1;
            }
        }
        assert!(checked > 50, "{checked}");
    }

    #[test]
    fn prop2_brute_force_lines_and_lambdas() {
        // Count every rational line of P^3 meeting the four secants directly
        // (a 4x4 determinant over K), and every lambda in P^1(F_p) where
        // det(M_alpha + lambda M_beta) = 0 (a generic determinant).
        let p = 7u64;
        let mut lines: Vec<([u64; 4], [u64; 4])> = Vec::new();
        for i in 0..4 {
            for j in i + 1..4 {
                let free1: Vec<usize> = (i + 1..4).filter(|&c| c != j).collect();
                let free2: Vec<usize> = (j + 1..4).collect();
                let n1 = p.pow(free1.len() as u32);
                let n2 = p.pow(free2.len() as u32);
                for a in 0..n1 {
                    for b in 0..n2 {
                        let mut r1 = [0u64; 4];
                        let mut r2 = [0u64; 4];
                        r1[i] = 1;
                        r2[j] = 1;
                        let mut x = a;
                        for &c in &free1 {
                            r1[c] = x % p;
                            x /= p;
                        }
                        let mut y = b;
                        for &c in &free2 {
                            r2[c] = y % p;
                            y /= p;
                        }
                        lines.push((r1, r2));
                    }
                }
            }
        }
        assert_eq!(lines.len() as u64, (p * p + 1) * (p * p + p + 1));
        let e = |k: &Ext, pt: &KPt| -> [El; 4] {
            let (u, v) = (&pt.u, &pt.v);
            [
                k.mul(&k.mul(u, u), u),
                k.mul(&k.mul(u, u), v),
                k.mul(&k.mul(u, v), v),
                k.mul(&k.mul(v, v), v),
            ]
        };
        let kdet = |k: &Ext, mut a: Vec<[El; 4]>| -> El {
            let mut det = k.one();
            for c in 0..4 {
                let Some(r) = (c..4).find(|&i| !k.is_zero(&a[i][c])) else {
                    return k.zero();
                };
                if r != c {
                    a.swap(r, c);
                    det = k.neg(&det);
                }
                det = k.mul(&det, &a[c][c]);
                let iv = k.inv(&a[c][c]);
                for i in c + 1..4 {
                    let f = k.mul(&a[i][c], &iv);
                    for j in c..4 {
                        a[i][j] = k.sub(&a[i][j], &k.mul(&f, &a[c][j]));
                    }
                }
            }
            det
        };
        let mut rng = SplitMix64(14);
        let (mut tested, mut by_class) = (0, [0u32; 3]);
        for _ in 0..200 {
            let f = rand_octic(p, &mut rng);
            let Some(w) = Weierstrass::new(&f, p, &mut rng) else {
                continue;
            };
            for sg in tractable_subgroups(&w, &mut rng).iter().take(3) {
                let o = analyze_subgroup(&w.f, sg, p);
                if o.m_rank != 4 || o.pencil_in_grassmannian {
                    continue;
                }
                let count = lines
                    .iter()
                    .filter(|(r1, r2)| {
                        sg.pairs.iter().all(|pr| {
                            let k = &pr.k;
                            let rows = vec![
                                r1.map(|x| k.c(x)),
                                r2.map(|x| k.c(x)),
                                e(k, &pr.a),
                                e(k, &pr.b),
                            ];
                            k.is_zero(&kdet(k, rows))
                        })
                    })
                    .count();
                let want = match o.disc_chi {
                    1 => 2,
                    0 => 1,
                    _ => 0,
                };
                assert_eq!(count, want, "{} disc_chi {}", sg.label, o.disc_chi);
                // lambda search with a generic 4x4 determinant over F_p.
                let (_, basis) = nullspace(&m_rows_fp(sg), 6, p);
                let mut roots = 0;
                for lam in 0..=p {
                    let g: Vec<u64> = (0..6)
                        .map(|i| {
                            if lam == p {
                                basis[1][i]
                            } else {
                                ad(basis[0][i], md(lam, basis[1][i], p), p)
                            }
                        })
                        .collect();
                    if det_fp(m_gamma(&g, p), p) == 0 {
                        roots += 1;
                    }
                }
                assert_eq!(roots, want);
                by_class[(o.disc_chi + 1) as usize] += 1;
                tested += 1;
            }
        }
        assert!(
            tested > 40 && by_class[0] > 5 && by_class[2] > 5,
            "{tested} {by_class:?}"
        );
    }

    #[test]
    fn normal_form_brute_force() {
        let p = 7u64;
        let mut rng = SplitMix64(15);
        let mut tested = 0;
        for _ in 0..150 {
            let f = rand_octic(p, &mut rng);
            let Some(w) = Weierstrass::new(&f, p, &mut rng) else {
                continue;
            };
            for sg in tractable_subgroups(&w, &mut rng).iter().take(2) {
                let o = analyze_subgroup(&w.f, sg, p);
                if o.m_rank != 4 || o.pencil_in_grassmannian {
                    continue;
                }
                let mut brute = std::collections::BTreeSet::new();
                for code in 0..p.pow(4) {
                    let (n1, n0, d1, d0) = (code % p, code / p % p, code / 49 % p, code / 343);
                    let map = TrigMap {
                        n: [1, 0, n1, n0],
                        d: [0, 1, d1, d0],
                        normal_form: true,
                    };
                    if degree_three(&map, p) && sg.pairs.iter().all(|pr| pair_condition(&map, pr)) {
                        brute.insert((n1, n0, d1, d0));
                    }
                }
                let mine: std::collections::BTreeSet<_> = o
                    .maps
                    .iter()
                    .filter(|m| m.degree3 && m.map.normal_form)
                    .filter_map(|m| m.map.smith_coeffs())
                    .collect();
                assert_eq!(mine, brute, "{}", sg.label);
                if !brute.is_empty() {
                    assert!(o.trig_rational());
                }
                tested += 1;
            }
        }
        assert!(tested > 30, "{tested}");
    }

    #[test]
    fn prop4_against_point_count_of_x() {
        let mut rng = SplitMix64(16);
        let mut total = Diag::default();
        for p in [101u64, 103] {
            for _ in 0..250 {
                let f = rand_octic(p, &mut rng);
                let Some(out) = analyze_curve(&f, p, &mut rng) else {
                    continue;
                };
                total.merge(&out.diag());
                total.merge(&x_check(&f, p, &out, false));
                let c = least_nonsquare(p);
                total.merge(&x_check(&pscale(&f, c, p), p, &out, true));
            }
        }
        assert_eq!(total.x_disagree, 0, "{total:?}");
        assert_eq!(total.x_count_error, 0, "{total:?}");
        assert!(total.x_agree > 100, "{total:?}");
        assert_eq!(total.pair_condition_fail, 0);
        assert_eq!(total.eval_mismatch, 0);
        assert_eq!(total.lemma4_fail, 0);
    }

    #[test]
    fn model_independence_and_twist() {
        let mut rng = SplitMix64(17);
        let p = 13u64;
        let c = least_nonsquare(p);
        let mut tested = 0;
        for _ in 0..300 {
            let f = rand_octic(p, &mut rng);
            let Some(out) = analyze_curve(&f, p, &mut rng) else {
                continue;
            };
            if out.n_tractable() == 0 {
                continue;
            }
            let key = |o: &CurveOutcome, tw: bool| {
                (o.n_tractable(), o.n_trig(), o.n_iso(tw), o.full_success(tw))
            };
            // Translate x -> x + 3 and scale by a square: the same curve.
            let g = pscale(&shift(&f, 3, p), 4, p);
            let og = analyze_curve(&g, p, &mut rng).unwrap();
            if og.diag().lambda_disagree == 0 && out.diag().lambda_disagree == 0 {
                assert_eq!(key(&out, false), key(&og, false));
            }
            // The twist: Proposition 4 flips, Proposition 2 does not.
            let ot = analyze_curve(&pscale(&f, c, p), p, &mut rng).unwrap();
            assert_eq!(ot.n_trig(), out.n_trig());
            assert_eq!(key(&ot, false), key(&out, true));
            tested += 1;
        }
        assert!(tested > 50);
    }

    #[test]
    fn predictions_reproduce_smith() {
        assert!((theorem2(0.25) - 0.18574).abs() < 1e-4);
        assert!((chebotarev_tractable() - 0.50159).abs() < 1e-4);
    }
}
