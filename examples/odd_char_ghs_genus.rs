//! Odd-characteristic GHS Weil descent: the genus of the cover for every
//! branch type and extension degree `n <= 8`, and a check of the `n = 3`
//! hyperelliptic family of Joux–Vitse.
//!
//! The file uses `std` only, so it also builds without the crate:
//!
//! ```text
//! rustc --edition 2021 -O examples/odd_char_ghs_genus.rs -o odd_char_ghs_genus
//! ./odd_char_ghs_genus research/odd_char_descent_genus3_20261014/results
//! # or: cargo run --release --example odd_char_ghs_genus -- <out-dir>
//! # tests: rustc --edition 2021 -O --test examples/odd_char_ghs_genus.rs -o t && ./t
//! ```
//!
//! **Part A (types).** `E/F_{q^n}`, `q` odd, and `x: E -> P^1` of degree 2
//! with branch set `B` (four points of `P^1`, `σ^n B = B`, where `σ` is the
//! `q`-Frobenius of the `F_q`-structure that `x` is taken over). Let
//! `S = ∪ σ^i B`, `r = |S|` and `V = span_{F_2}{1_{σ^i B}}`, with
//! `m = dim V`. The compositum of the `n` conjugate double covers is a
//! `(Z/2)^m` cover of `P^1`, tamely branched with index 2 over the `r`
//! points of `S`. Riemann–Hurwitz gives `g = 1 + 2^{m−2}(r−4)`. The tool
//! checks this against the sum of the genera of the quadratic subcovers
//! (Kani–Rosen, the identity Diem uses in Compositio 2011, Prop. 4.10). The
//! genus depends only on the type: the `σ`-orbits that `S` is made of and
//! the positions `B` takes in them. `q` enters only through how many branch
//! sets of each type exist; those counts are exact here and must add up to
//! `q^{4n} − q^{2n}`. That total is the number of reduced degree-4 divisors
//! on `P^1` over `F_{q^n}`.
//!
//! **Part B (`n = 3` hyperelliptic family).** Joux–Vitse (Eurocrypt 2012,
//! Appendix A.1) take `E: y^2 = h(x)(x−α)(x−σα)` with `α ∈ F_{q^3} \ F_q`
//! and `h ∈ F_q[x]` of degree 1 or 2. Their `H: y^2 = 4N^2 h(F/(4N))` lives
//! over `F_q`, with `N` the minimal polynomial of `α` and
//! `F = N·(x + φ + φ^σ + φ^{σ^2})`. The tool counts points on `H` over
//! `F_p`, `F_{p^2}` and `F_{p^3}`, and on `E` over `F_{p^3}`. It checks
//! `L_H(T) = 1 − tT^3 + p^3T^6`, which is `L` of `Res_{F_{p^3}/F_p} E`. It
//! records the Frobenius cycle type on the eight Weierstrass points of `H`
//! and the number of Smith's rational tractable subgroups that this cycle
//! type allows (Smith, J. Cryptology 2009, Lemma 1).

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

// ===========================================================================
// Small unsigned big integers (base 10^9), enough for exact type counts.

const BASE: u64 = 1_000_000_000;

#[derive(Clone, Debug, PartialEq, Eq)]
struct Big(Vec<u64>); // little-endian limbs, base 10^9, no trailing zero limbs

impl Big {
    fn from_u64(mut v: u64) -> Big {
        let mut limbs = Vec::new();
        while v > 0 {
            limbs.push(v % BASE);
            v /= BASE;
        }
        Big(limbs)
    }
    fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
    fn trim(mut self) -> Big {
        while self.0.last() == Some(&0) {
            self.0.pop();
        }
        self
    }
    fn add(&self, o: &Big) -> Big {
        let n = self.0.len().max(o.0.len());
        let mut out = Vec::with_capacity(n + 1);
        let mut carry = 0;
        for i in 0..n {
            let s = self.0.get(i).copied().unwrap_or(0) + o.0.get(i).copied().unwrap_or(0) + carry;
            out.push(s % BASE);
            carry = s / BASE;
        }
        if carry > 0 {
            out.push(carry);
        }
        Big(out).trim()
    }
    /// `self − o`; panics if `o > self`.
    fn sub(&self, o: &Big) -> Big {
        assert!(self.cmp_big(o) != std::cmp::Ordering::Less, "negative Big");
        let mut out = Vec::with_capacity(self.0.len());
        let mut borrow = 0i64;
        for i in 0..self.0.len() {
            let mut d = self.0[i] as i64 - o.0.get(i).copied().unwrap_or(0) as i64 - borrow;
            if d < 0 {
                d += BASE as i64;
                borrow = 1;
            } else {
                borrow = 0;
            }
            out.push(d as u64);
        }
        Big(out).trim()
    }
    fn mul_small(&self, k: u64) -> Big {
        assert!(k < BASE);
        let mut out = Vec::with_capacity(self.0.len() + 1);
        let mut carry = 0u64;
        for &l in &self.0 {
            let v = l * k + carry;
            out.push(v % BASE);
            carry = v / BASE;
        }
        while carry > 0 {
            out.push(carry % BASE);
            carry /= BASE;
        }
        Big(out).trim()
    }
    fn mul(&self, o: &Big) -> Big {
        if self.is_zero() || o.is_zero() {
            return Big(Vec::new());
        }
        let mut acc = vec![0u128; self.0.len() + o.0.len() + 1];
        for (i, &a) in self.0.iter().enumerate() {
            for (j, &b) in o.0.iter().enumerate() {
                acc[i + j] += a as u128 * b as u128;
            }
        }
        let mut out = Vec::with_capacity(acc.len());
        let mut carry = 0u128;
        for v in acc {
            let s = v + carry;
            out.push((s % BASE as u128) as u64);
            carry = s / BASE as u128;
        }
        assert_eq!(carry, 0);
        Big(out).trim()
    }
    /// Exact division by a small integer; panics on a non-zero remainder.
    fn div_exact(&self, k: u64) -> Big {
        assert!(k > 0 && k < BASE);
        let mut out = vec![0u64; self.0.len()];
        let mut rem = 0u64;
        for i in (0..self.0.len()).rev() {
            let cur = rem * BASE + self.0[i];
            out[i] = cur / k;
            rem = cur % k;
        }
        assert_eq!(rem, 0, "inexact division");
        Big(out).trim()
    }
    fn pow_small(base: u64, e: u32) -> Big {
        let mut r = Big::from_u64(1);
        for _ in 0..e {
            r = r.mul_small(base);
        }
        r
    }
    fn cmp_big(&self, o: &Big) -> std::cmp::Ordering {
        self.0
            .len()
            .cmp(&o.0.len())
            .then_with(|| self.0.iter().rev().cmp(o.0.iter().rev()))
    }
    fn log10(&self) -> f64 {
        match self.0.len() {
            0 => f64::NEG_INFINITY,
            1 => (self.0[0] as f64).log10(),
            n => {
                let top = self.0[n - 1] as f64 * BASE as f64 + self.0[n - 2] as f64;
                top.log10() + 9.0 * (n - 2) as f64
            }
        }
    }
}

impl std::fmt::Display for Big {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.0.split_last() {
            None => write!(f, "0"),
            Some((top, rest)) => {
                write!(f, "{top}")?;
                for l in rest.iter().rev() {
                    write!(f, "{l:09}")?;
                }
                Ok(())
            }
        }
    }
}

// ===========================================================================
// Part A: branch types.

fn gcd(a: u32, b: u32) -> u32 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn full(k: u32) -> u64 {
    if k == 64 {
        u64::MAX
    } else {
        (1u64 << k) - 1
    }
}

/// Rotate a subset of `Z/k` (bitmask) by `s`: position `i` goes to `i + s`.
fn rot(mask: u64, k: u32, s: u32) -> u64 {
    let s = s % k;
    if s == 0 {
        return mask;
    }
    ((mask << s) | (mask >> (k - s))) & full(k)
}

fn canon(mask: u64, k: u32) -> u64 {
    (0..k).map(|s| rot(mask, k, s)).min().unwrap()
}

fn stabilizer(mask: u64, k: u32) -> u32 {
    (0..k).filter(|&s| rot(mask, k, s) == mask).count() as u32
}

/// One `σ`-orbit of size `k` and the positions (canonical up to rotation)
/// that the branch set occupies in it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Part {
    k: u32,
    mask: u64,
}

/// All parts that can occur for extension degree `n`. `σ^n` acts on an
/// orbit of size `k` as a shift by `n`, so `B ∩ O` is a union of cosets of
/// `gcd(n, k)·Z/k`, each of size `k / gcd(n, k) <= 4`.
fn parts(n: u32) -> Vec<Part> {
    let mut out = BTreeSet::new();
    for k in 1..=4 * n {
        let g = gcd(n, k);
        let l = k / g;
        if l > 4 {
            continue;
        }
        for jmask in 1u64..(1u64 << g) {
            if jmask.count_ones() * l > 4 {
                continue;
            }
            let mut mask = 0u64;
            for j in 0..g {
                if (jmask >> j) & 1 == 1 {
                    for t in 0..l {
                        mask |= 1u64 << (j + t * g);
                    }
                }
            }
            out.insert(Part {
                k,
                mask: canon(mask, k),
            });
        }
    }
    out.into_iter().collect()
}

/// Every type for degree `n`: a multiset of parts (distinct orbits) holding
/// four points in total.
fn types(n: u32) -> Vec<Vec<Part>> {
    fn rec(ps: &[Part], start: usize, rem: u32, cur: &mut Vec<Part>, out: &mut Vec<Vec<Part>>) {
        if rem == 0 {
            out.push(cur.clone());
            return;
        }
        for i in start..ps.len() {
            let c = ps[i].mask.count_ones();
            if c <= rem {
                cur.push(ps[i]);
                rec(ps, i, rem - c, cur, out);
                cur.pop();
            }
        }
    }
    let ps = parts(n);
    let mut out = Vec::new();
    rec(&ps, 0, 4, &mut Vec::new(), &mut out);
    out
}

#[derive(Clone, Debug)]
struct Analysis {
    r: u32,
    m: u32,
    genus: i64,
    period: u32,
    hyper_via_g: bool,
}

fn xor_basis(vs: &[u64]) -> Vec<u64> {
    let mut basis: Vec<u64> = Vec::new();
    for &v in vs {
        let mut x = v;
        for &b in &basis {
            x = x.min(x ^ b);
        }
        if x != 0 {
            basis.push(x);
            basis.sort_unstable_by(|a, b| b.cmp(a));
        }
    }
    basis
}

fn span(basis: &[u64]) -> Vec<u64> {
    (0u64..(1u64 << basis.len()))
        .map(|c| {
            basis
                .iter()
                .enumerate()
                .filter(|(i, _)| (c >> i) & 1 == 1)
                .fold(0u64, |a, (_, &b)| a ^ b)
        })
        .collect()
}

/// `4g − 4 = 2^w (r − 4)` for a `(Z/2)^w` cover of `P^1` branched with
/// index 2 over `r` points (`w = 0`, `r = 0`: `P^1` itself).
fn genus_of(w: u32, r: u32) -> i64 {
    let num = (1i64 << w) * (r as i64 - 4) + 4;
    assert_eq!(num % 4, 0, "non-integral genus w={w} r={r}");
    num / 4
}

fn analyse(n: u32, ty: &[Part]) -> Analysis {
    let mut offs = Vec::new();
    let mut r = 0u32;
    for p in ty {
        offs.push(r);
        r += p.k;
    }
    assert!(r <= 64);
    let shift = |v: u64| -> u64 {
        let mut out = 0u64;
        for (p, &o) in ty.iter().zip(&offs) {
            out |= rot((v >> o) & full(p.k), p.k, 1) << o;
        }
        out
    };
    let b = ty
        .iter()
        .zip(&offs)
        .fold(0u64, |a, (p, &o)| a | (p.mask << o));
    let mut vecs = vec![b];
    for i in 1..n as usize {
        vecs.push(shift(vecs[i - 1]));
    }
    assert_eq!(shift(vecs[n as usize - 1]), b, "σ^n B != B");
    let union = vecs.iter().fold(0u64, |a, &v| a | v);
    assert_eq!(union, full(r), "S is not the union of the σ^i B");
    let period = (1..=n).find(|&i| vecs[(i % n) as usize] == b).unwrap();
    let basis = xor_basis(&vecs);
    let m = basis.len() as u32;
    let all = span(&basis);
    // Kani–Rosen: g(C) is the sum of the genera of the quadratic subcovers.
    let kr: i64 = all
        .iter()
        .filter(|&&v| v != 0)
        .map(|&v| {
            assert_eq!(v.count_ones() % 2, 0);
            v.count_ones() as i64 / 2 - 1
        })
        .sum();
    let genus = genus_of(m, r);
    assert_eq!(genus, kr, "Riemann–Hurwitz and Kani–Rosen disagree");
    // Hyperelliptic through an involution of the cover group G: some
    // hyperplane W of V whose (Z/2)^{m−1} cover has genus 0. Genus 0 needs
    // 2^{m−1}(4 − r_W) = 4, so m <= 3; larger m is skipped (exactly, not as a
    // heuristic), which keeps n up to 14 cheap.
    let mut hyper_via_g = false;
    if genus >= 2 && m <= 3 {
        for lam in 1u64..(1u64 << m) {
            let or_w = all
                .iter()
                .enumerate()
                .filter(|(c, _)| (lam & *c as u64).count_ones().is_multiple_of(2))
                .fold(0u64, |a, (_, &v)| a | v);
            if genus_of(m - 1, or_w.count_ones()) == 0 {
                hyper_via_g = true;
                break;
            }
        }
    }
    Analysis {
        r,
        m,
        genus,
        period,
        hyper_via_g,
    }
}

fn mobius(n: u32) -> i32 {
    let mut n = n;
    let mut mu = 1;
    let mut d = 2;
    while d * d <= n {
        if n.is_multiple_of(d) {
            n /= d;
            if n.is_multiple_of(d) {
                return 0;
            }
            mu = -mu;
        }
        d += 1;
    }
    if n > 1 {
        mu = -mu;
    }
    mu
}

/// Number of `σ`-orbits of size `k` on `P^1(F_q-bar)`.
fn orbits(q: u64, k: u32) -> Big {
    if k == 1 {
        return Big::from_u64(q + 1);
    }
    let (mut pos, mut neg) = (Big::from_u64(0), Big::from_u64(0));
    for d in 1..=k {
        if k.is_multiple_of(d) {
            let t = Big::pow_small(q, d);
            match mobius(k / d) {
                1 => pos = pos.add(&t),
                -1 => neg = neg.add(&t),
                _ => {}
            }
        }
    }
    pos.sub(&neg).div_exact(k as u64)
}

/// Number of branch sets `B` (4-subsets of `P^1(F_q-bar)` with
/// `σ^n B = B`) of the given type.
fn count_type(q: u64, ty: &[Part]) -> Big {
    let mut by_k: BTreeMap<u32, Vec<Part>> = BTreeMap::new();
    for p in ty {
        by_k.entry(p.k).or_default().push(*p);
    }
    let mut total = Big::from_u64(1);
    for (k, ps) in &by_k {
        let n_orb = orbits(q, *k);
        for i in 0..ps.len() {
            let avail = n_orb.sub(&Big::from_u64(i as u64));
            total = total.mul(&avail);
            total = total.mul_small((*k / stabilizer(ps[i].mask, *k)) as u64);
        }
        let mut mult: BTreeMap<Part, u64> = BTreeMap::new();
        for p in ps {
            *mult.entry(*p).or_default() += 1;
        }
        for &c in mult.values() {
            for f in 2..=c {
                total = total.div_exact(f);
            }
        }
    }
    total
}

fn describe(ty: &[Part]) -> (String, String) {
    let sizes: Vec<String> = ty.iter().map(|p| p.k.to_string()).collect();
    let pats: Vec<String> = ty
        .iter()
        .map(|p| {
            let pos: Vec<String> = (0..p.k)
                .filter(|i| (p.mask >> i) & 1 == 1)
                .map(|i| i.to_string())
                .collect();
            format!("{}:{{{}}}", p.k, pos.join(","))
        })
        .collect();
    (sizes.join("+"), pats.join(";"))
}

struct TypeRow {
    n: u32,
    ty: Vec<Part>,
    a: Analysis,
}

fn all_types(nmax: u32) -> Vec<TypeRow> {
    let mut rows = Vec::new();
    for n in 1..=nmax {
        for ty in types(n) {
            let a = analyse(n, &ty);
            rows.push(TypeRow { n, ty, a });
        }
    }
    rows
}

// ===========================================================================
// Part B: finite fields F_{p^k} (k <= 3) and polynomials.

type El = [u64; 3];

#[derive(Clone)]
struct Gf {
    p: u64,
    k: usize,
    modulus: Vec<u64>, // monic, length k + 1, ascending
    chi: Vec<i8>,      // quadratic character by element index
}

fn poly_has_root_mod_p(m: &[u64], p: u64) -> bool {
    (0..p).any(|x| m.iter().rev().fold(0u64, |acc, &c| (acc * x + c) % p) == 0)
}

impl Gf {
    fn new(p: u64, k: usize) -> Gf {
        assert!((1..=3).contains(&k));
        let modulus = if k == 1 {
            vec![0, 1]
        } else {
            let mut found = None;
            for code in 0..p.pow(k as u32) {
                let mut m: Vec<u64> = (0..k).map(|i| (code / p.pow(i as u32)) % p).collect();
                m.push(1);
                if !poly_has_root_mod_p(&m, p) {
                    found = Some(m);
                    break;
                }
            }
            found.expect("irreducible modulus")
        };
        Gf::with_modulus(p, modulus)
    }
    fn with_modulus(p: u64, modulus: Vec<u64>) -> Gf {
        let k = modulus.len() - 1;
        assert!(k == 1 || !poly_has_root_mod_p(&modulus, p));
        let mut f = Gf {
            p,
            k,
            modulus,
            chi: Vec::new(),
        };
        let size = f.size();
        let mut chi = vec![-1i8; size];
        chi[0] = 0;
        for i in 1..size {
            let a = f.el(i);
            chi[f.idx(&f.mul(&a, &a))] = 1;
        }
        f.chi = chi;
        f
    }
    fn size(&self) -> usize {
        self.p.pow(self.k as u32) as usize
    }
    fn el(&self, i: usize) -> El {
        let mut e = [0u64; 3];
        let mut v = i as u64;
        for c in e.iter_mut().take(self.k) {
            *c = v % self.p;
            v /= self.p;
        }
        e
    }
    fn idx(&self, a: &El) -> usize {
        (0..self.k).rev().fold(0u64, |acc, i| acc * self.p + a[i]) as usize
    }
    fn int(&self, v: i64) -> El {
        [v.rem_euclid(self.p as i64) as u64, 0, 0]
    }
    fn add(&self, a: &El, b: &El) -> El {
        let mut c = [0u64; 3];
        for i in 0..self.k {
            c[i] = (a[i] + b[i]) % self.p;
        }
        c
    }
    fn sub(&self, a: &El, b: &El) -> El {
        let mut c = [0u64; 3];
        for i in 0..self.k {
            c[i] = (a[i] + self.p - b[i]) % self.p;
        }
        c
    }
    fn mul(&self, a: &El, b: &El) -> El {
        let (p, k) = (self.p, self.k);
        let mut t = [0u64; 5];
        for i in 0..k {
            for j in 0..k {
                t[i + j] = (t[i + j] + a[i] * b[j]) % p;
            }
        }
        for d in (k..2 * k - 1).rev() {
            let c = t[d];
            if c == 0 {
                continue;
            }
            t[d] = 0;
            for i in 0..k {
                t[d - k + i] = (t[d - k + i] + (p - c) * self.modulus[i]) % p;
            }
        }
        [t[0], t[1], t[2]]
    }
    fn pow(&self, a: &El, mut e: u64) -> El {
        let mut r = self.int(1);
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
    fn frob(&self, a: &El) -> El {
        self.pow(a, self.p)
    }
    fn chi_of(&self, a: &El) -> i64 {
        self.chi[self.idx(a)] as i64
    }
    fn in_prime_field(&self, a: &El) -> bool {
        a[1] == 0 && a[2] == 0
    }
}

type Poly = Vec<El>; // ascending coefficients over a Gf

fn pmul(f: &Gf, a: &Poly, b: &Poly) -> Poly {
    let mut out = vec![[0u64; 3]; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = f.add(&out[i + j], &f.mul(x, y));
        }
    }
    out
}

fn padd(f: &Gf, a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    (0..n)
        .map(|i| {
            let x = a.get(i).copied().unwrap_or([0; 3]);
            let y = b.get(i).copied().unwrap_or([0; 3]);
            f.add(&x, &y)
        })
        .collect()
}

fn pscale(f: &Gf, c: &El, a: &Poly) -> Poly {
    a.iter().map(|x| f.mul(c, x)).collect()
}

fn peval(f: &Gf, a: &Poly, x: &El) -> El {
    a.iter()
        .rev()
        .fold([0u64; 3], |acc, c| f.add(&f.mul(&acc, x), c))
}

fn ptrim(a: &mut Poly) {
    while a.len() > 1 && a.last() == Some(&[0u64; 3]) {
        a.pop();
    }
}

/// Polynomials over F_p as plain `u64` coefficient vectors (ascending).
mod fp {
    pub fn trim(a: &mut Vec<u64>) {
        while a.len() > 1 && *a.last().unwrap() == 0 {
            a.pop();
        }
    }
    pub fn inv(a: u64, p: u64) -> u64 {
        let mut r = 1u64;
        let mut b = a % p;
        let mut e = p - 2;
        while e > 0 {
            if e & 1 == 1 {
                r = r * b % p;
            }
            b = b * b % p;
            e >>= 1;
        }
        r
    }
    pub fn rem(a: &[u64], m: &[u64], p: u64) -> Vec<u64> {
        let mut r = a.to_vec();
        trim(&mut r);
        let dm = m.len() - 1;
        let li = inv(m[dm], p);
        while r.len() > dm && !(r.len() == 1 && r[0] == 0) {
            let d = r.len() - 1;
            let c = r[d] * li % p;
            for i in 0..=dm {
                r[d - dm + i] = (r[d - dm + i] + (p - c) * m[i]) % p;
            }
            r.pop();
            trim(&mut r);
            if r.len() <= dm {
                break;
            }
        }
        if r.is_empty() {
            r.push(0);
        }
        r
    }
    pub fn mulmod(a: &[u64], b: &[u64], m: &[u64], p: u64) -> Vec<u64> {
        let mut t = vec![0u64; a.len() + b.len() - 1];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                t[i + j] = (t[i + j] + x * y) % p;
            }
        }
        rem(&t, m, p)
    }
    pub fn is_zero(a: &[u64]) -> bool {
        a.iter().all(|&c| c == 0)
    }
    pub fn gcd(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
        let (mut x, mut y) = (a.to_vec(), b.to_vec());
        trim(&mut x);
        trim(&mut y);
        while !is_zero(&y) {
            let r = rem(&x, &y, p);
            x = y;
            y = r;
        }
        let li = inv(*x.last().unwrap(), p);
        x.iter().map(|c| c * li % p).collect()
    }
    pub fn div(a: &[u64], m: &[u64], p: u64) -> Vec<u64> {
        let mut r = a.to_vec();
        trim(&mut r);
        let dm = m.len() - 1;
        let li = inv(m[dm], p);
        if r.len() <= dm {
            return vec![0];
        }
        let mut q = vec![0u64; r.len() - dm];
        while r.len() > dm {
            let d = r.len() - 1;
            let c = r[d] * li % p;
            q[d - dm] = c;
            for i in 0..=dm {
                r[d - dm + i] = (r[d - dm + i] + (p - c) * m[i]) % p;
            }
            r.pop();
        }
        q
    }
    /// Sorted degrees of the irreducible factors of a squarefree `f`
    /// (distinct-degree factorisation); `None` if `f` is not squarefree.
    pub fn cycle_type(f: &[u64], p: u64) -> Option<Vec<u32>> {
        let mut g = f.to_vec();
        trim(&mut g);
        let deriv: Vec<u64> = g
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, &c)| (i as u64 % p) * c % p)
            .collect();
        if is_zero(&deriv) || gcd(&g, &deriv, p).len() != 1 {
            return None;
        }
        let mut out = Vec::new();
        let mut xp = vec![0u64, 1];
        let mut d = 0u32;
        while g.len() > 1 {
            d += 1;
            if 2 * d > (g.len() - 1) as u32 {
                out.push((g.len() - 1) as u32);
                break;
            }
            // xp <- xp^p mod g
            let mut acc = vec![1u64];
            let mut base = rem(&xp, &g, p);
            let mut e = p;
            while e > 0 {
                if e & 1 == 1 {
                    acc = mulmod(&acc, &base, &g, p);
                }
                base = mulmod(&base, &base, &g, p);
                e >>= 1;
            }
            xp = acc;
            let mut h = xp.clone();
            h.resize(h.len().max(2), 0);
            h[1] = (h[1] + p - 1) % p;
            trim(&mut h);
            let gd = gcd(&g, &h, p);
            let deg = (gd.len() - 1) as u32;
            if deg > 0 {
                for _ in 0..deg / d {
                    out.push(d);
                }
                g = div(&g, &gd, p);
                xp = rem(&xp, &g, p);
            }
        }
        out.sort_unstable();
        Some(out)
    }
}

/// Number of perfect matchings of eight points into four pairs that a
/// permutation with the given cycle type maps to itself: the number of
/// `F_q`-rational tractable subgroups `#S(H)` (Smith 2009, Prop. 1, Def. 1).
fn stable_matchings(ct: &[u32]) -> u32 {
    assert_eq!(ct.iter().sum::<u32>(), 8);
    let mut perm = [0usize; 8];
    let mut start = 0usize;
    for &c in ct {
        let c = c as usize;
        for i in 0..c {
            perm[start + i] = start + (i + 1) % c;
        }
        start += c;
    }
    fn rec(perm: &[usize; 8], used: u8, pairs: &mut Vec<(usize, usize)>, count: &mut u32) {
        if used == 0xff {
            let key = |a: usize, b: usize| (a.min(b), a.max(b));
            let set: BTreeSet<(usize, usize)> = pairs.iter().map(|&(a, b)| key(a, b)).collect();
            if pairs
                .iter()
                .all(|&(a, b)| set.contains(&key(perm[a], perm[b])))
            {
                *count += 1;
            }
            return;
        }
        let a = (0..8).find(|&i| (used >> i) & 1 == 0).unwrap();
        for b in a + 1..8 {
            if (used >> b) & 1 == 0 {
                pairs.push((a, b));
                rec(perm, used | (1 << a) | (1 << b), pairs, count);
                pairs.pop();
            }
        }
    }
    let mut count = 0;
    rec(&perm, 0, &mut Vec::new(), &mut count);
    count
}

/// Smith 2009, Lemma 1: `#S(H)` by the degrees of the `F_q`-irreducible
/// factors of the binary octic. Transcribed from the source.
fn smith_lemma1(ct: &[u32]) -> u32 {
    let mut s = ct.to_vec();
    s.sort_unstable_by(|a, b| b.cmp(a));
    match s.as_slice() {
        [8] | [6, 2] | [6, 1, 1] | [4, 2, 1, 1] => 1,
        [4, 2, 2] | [4, 1, 1, 1, 1] | [3, 3, 2] | [3, 3, 1, 1] => 3,
        [4, 4] => 5,
        [2, 2, 2, 1, 1] => 7,
        [2, 2, 1, 1, 1, 1] => 9,
        [2, 1, 1, 1, 1, 1, 1] => 15,
        [2, 2, 2, 2] => 25,
        [1, 1, 1, 1, 1, 1, 1, 1] => 105,
        _ => 0,
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Arm {
    Linear,    // h = x            (deg H = 7)
    Split,     // h = x(x − 1)     (deg H = 8)
    Conjugate, // h = x^2 − c, c the least non-residue (deg H = 8)
}

impl Arm {
    fn name(&self) -> &'static str {
        match self {
            Arm::Linear => "linear",
            Arm::Split => "split",
            Arm::Conjugate => "conjugate",
        }
    }
}

#[derive(Clone, Debug)]
struct JvRow {
    p: u64,
    arm: Arm,
    h: Vec<u64>,
    alpha: El,
    modulus: Vec<u64>,
    n_e: u64,
    t: i64,
    e_k_2torsion: u64,
    f_h: Vec<u64>,
    n: [u64; 3],
    l_ok: bool,
    jac_order: i64,
    cycle_type: Option<Vec<u32>>,
    stable: Option<u32>,
    lemma1: Option<u32>,
}

fn least_nonresidue(p: u64) -> u64 {
    (2..p)
        .find(|&c| (1..p).all(|y| y * y % p != c))
        .expect("p odd")
}

/// The Joux–Vitse model `H: y^2 = 4N^2 h(F/(4N))` and its checks.
fn jv_instance(k3: &Gf, k2: &Gf, k1: &Gf, alpha: El, arm: Arm) -> JvRow {
    let p = k3.p;
    let x: Poly = vec![[0; 3], k3.int(1)];
    let lin = |c: &El| -> Poly { vec![k3.sub(&[0; 3], c), k3.int(1)] };
    let a0 = alpha;
    let a1 = k3.frob(&a0);
    let a2 = k3.frob(&a1);
    assert_eq!(k3.frob(&a2), a0);
    assert!(a0 != a1, "α in F_p");
    let n_poly = pmul(k3, &pmul(k3, &lin(&a0), &lin(&a1)), &lin(&a2));
    let d0 = k3.mul(&k3.sub(&a0, &a2), &k3.sub(&a1, &a2));
    let d1 = k3.frob(&d0);
    let d2 = k3.frob(&d1);
    let tr = k3.add(&k3.add(&a0, &a1), &a2);
    // F = xN + Tr(α)N + D(x−α)(x−σα) + σD(x−σα)(x−σ²α) + σ²D(x−α)(x−σ²α)
    let mut f_poly = padd(k3, &pmul(k3, &x, &n_poly), &pscale(k3, &tr, &n_poly));
    f_poly = padd(
        k3,
        &f_poly,
        &pscale(k3, &d0, &pmul(k3, &lin(&a0), &lin(&a1))),
    );
    f_poly = padd(
        k3,
        &f_poly,
        &pscale(k3, &d1, &pmul(k3, &lin(&a1), &lin(&a2))),
    );
    f_poly = padd(
        k3,
        &f_poly,
        &pscale(k3, &d2, &pmul(k3, &lin(&a0), &lin(&a2))),
    );
    for c in n_poly.iter().chain(f_poly.iter()) {
        assert!(k3.in_prime_field(c), "N or F not over F_p");
    }
    let h: Vec<u64> = match arm {
        Arm::Linear => vec![0, 1],
        Arm::Split => vec![0, p - 1, 1],
        Arm::Conjugate => vec![p - least_nonresidue(p), 0, 1],
    };
    // 4N^2 h(F/(4N)) = 4 h0 N^2 + h1 F N + (h2/4) F^2
    let four = k3.int(4);
    let n2 = pmul(k3, &n_poly, &n_poly);
    let mut f_h = pscale(k3, &k3.mul(&four, &k3.int(h[0] as i64)), &n2);
    f_h = padd(
        k3,
        &f_h,
        &pscale(k3, &k3.int(h[1] as i64), &pmul(k3, &f_poly, &n_poly)),
    );
    if h.len() == 3 {
        let quarter = k3.int(fp::inv(4, p) as i64);
        let c = k3.mul(&quarter, &k3.int(h[2] as i64));
        f_h = padd(k3, &f_h, &pscale(k3, &c, &pmul(k3, &f_poly, &f_poly)));
    }
    ptrim(&mut f_h);
    let f_h_fp: Vec<u64> = f_h.iter().map(|c| c[0]).collect();
    let deg_h = f_h_fp.len() - 1;
    // E over K: y^2 = h(x)(x − α)(x − σα).
    let h_k: Poly = h.iter().map(|&c| k3.int(c as i64)).collect();
    let e_poly = pmul(k3, &h_k, &pmul(k3, &lin(&a0), &lin(&a1)));
    let deg_e = e_poly.len() - 1;
    let mut n_e: i64 = 0;
    for i in 0..k3.size() {
        n_e += 1 + k3.chi_of(&peval(k3, &e_poly, &k3.el(i)));
    }
    n_e += if deg_e == 3 {
        1
    } else {
        1 + k3.chi_of(&e_poly[deg_e])
    };
    let q3 = p.pow(3) as i64;
    let t = q3 + 1 - n_e;
    // #E(K)[2]: K-roots of the cubic, or of the resolvent cubic of the quartic.
    let e_k_2torsion = {
        let cubic: Poly = if deg_e == 3 {
            e_poly.clone()
        } else {
            let li = k3.pow(&e_poly[4], k3.size() as u64 - 2);
            let c: Vec<El> = e_poly.iter().map(|v| k3.mul(v, &li)).collect();
            let (b, cc, d, e) = (c[3], c[2], c[1], c[0]);
            // y^3 − c y^2 + (bd − 4e) y − (b^2 e − 4ce + d^2)
            let t1 = k3.sub(&k3.mul(&b, &d), &k3.mul(&four, &e));
            let t0 = k3.add(
                &k3.sub(
                    &k3.mul(&k3.mul(&b, &b), &e),
                    &k3.mul(&four, &k3.mul(&cc, &e)),
                ),
                &k3.mul(&d, &d),
            );
            vec![k3.sub(&[0; 3], &t0), t1, k3.sub(&[0; 3], &cc), k3.int(1)]
        };
        1 + (0..k3.size())
            .filter(|&i| peval(k3, &cubic, &k3.el(i)) == [0; 3])
            .count() as u64
    };
    // Points on H over F_p, F_{p^2}, F_{p^3}.
    let count_h = |f: &Gf| -> u64 {
        let poly: Poly = f_h_fp.iter().map(|&c| f.int(c as i64)).collect();
        let mut s: i64 = 0;
        for i in 0..f.size() {
            s += 1 + f.chi_of(&peval(f, &poly, &f.el(i)));
        }
        s += if deg_h == 7 {
            1
        } else {
            1 + f.chi_of(&f.int(f_h_fp[deg_h] as i64))
        };
        s as u64
    };
    let n = [count_h(k1), count_h(k2), count_h(k3)];
    let pi = p as i64;
    let l_ok = n[0] as i64 == pi + 1 && n[1] as i64 == pi * pi + 1 && n[2] as i64 == q3 + 1 - 3 * t;
    // L(1) from the point counts (genus 3, Newton's identities).
    let s1 = pi + 1 - n[0] as i64;
    let s2 = pi * pi + 1 - n[1] as i64;
    let s3 = q3 + 1 - n[2] as i64;
    let e1 = s1;
    let e2 = (s1 * s1 - s2) / 2;
    let e3 = (s1 * s1 * s1 - 3 * s1 * s2 + 2 * s3) / 6;
    let (c1, c2, c3) = (-e1, e2, -e3);
    let jac_order = 1 + c1 + c2 + c3 + pi * c2 + pi * pi * c1 + q3;
    let mut cycle_type = fp::cycle_type(&f_h_fp, p);
    if let Some(ct) = cycle_type.as_mut() {
        if deg_h == 7 {
            ct.push(1); // the Weierstrass point at infinity
            ct.sort_unstable();
        }
    }
    let stable = cycle_type.as_ref().map(|ct| stable_matchings(ct));
    let lemma1 = cycle_type.as_ref().map(|ct| smith_lemma1(ct));
    JvRow {
        p,
        arm,
        h,
        alpha,
        modulus: k3.modulus.clone(),
        n_e: n_e as u64,
        t,
        e_k_2torsion,
        f_h: f_h_fp,
        n,
        l_ok,
        jac_order,
        cycle_type,
        stable,
        lemma1,
    }
}

fn join<T: std::fmt::Display>(v: &[T], sep: &str) -> String {
    v.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(sep)
}

// ===========================================================================
// Driver.

const NMAX: u32 = 8;
/// Degrees for the minimal-genus table (genus only, no counts).
const NMAX_MIN: u32 = 14;
const PRIMES: [u64; 9] = [5, 7, 11, 13, 17, 19, 23, 29, 31];

/// Per `n`: the smallest genus over types of period exactly `n` (the
/// branch set is not fixed by any `σ^d`, `d < n`; this holds for every `x`
/// when `F_q(j(E)) = F_{q^n}`), and whether any type at all has `g = n`.
fn min_genus_table(nmax: u32) -> String {
    let mut csv = String::from(
        "n,types,min_genus_period_n,m_r_at_min,types_at_min,any_type_genus_eq_n,generic_genus,\
         check_total_q5\n",
    );
    for n in 1..=nmax {
        let ts = types(n);
        let an: Vec<Analysis> = ts.iter().map(|ty| analyse(n, ty)).collect();
        let per_n: Vec<&Analysis> = an.iter().filter(|a| a.period == n).collect();
        let gmin = per_n.iter().map(|a| a.genus).min().unwrap();
        let at_min: BTreeSet<(u32, u32)> = per_n
            .iter()
            .filter(|a| a.genus == gmin)
            .map(|a| (a.m, a.r))
            .collect();
        let n_at_min = per_n.iter().filter(|a| a.genus == gmin).count();
        let mr: Vec<String> = at_min.iter().map(|(m, r)| format!("m{m}r{r}")).collect();
        let any_eq = an.iter().any(|a| a.genus == n as i64);
        // Completeness at q = 5: the type counts must add up to 5^{4n} − 5^{2n}.
        let total = ts
            .iter()
            .fold(Big::from_u64(0), |acc, ty| acc.add(&count_type(5, ty)));
        let q2 = Big::pow_small(5, 2 * n);
        let ok = total == q2.mul(&q2).sub(&q2);
        writeln!(
            csv,
            "{n},{},{gmin},{},{n_at_min},{any_eq},{},{}",
            ts.len(),
            mr.join(" "),
            (n as i64 - 1) * (1i64 << n) + 1,
            if ok { "ok" } else { "FAIL" }
        )
        .unwrap();
    }
    csv
}

fn main() {
    let out = std::env::args()
        .nth(1)
        .expect("usage: odd_char_ghs_genus <out-dir>");
    let out = Path::new(&out);
    fs::create_dir_all(out).expect("create out dir");
    let mut failures = 0u32;

    // ---- Part A: types.
    let rows = all_types(NMAX);
    let mut csv = String::from(
        "n,type,orbit_sizes,pattern,r,m,genus,period,hyperelliptic_via_G,genus3,g_equals_n\n",
    );
    let mut idx_by_n: BTreeMap<u32, u32> = BTreeMap::new();
    for row in &rows {
        let i = idx_by_n.entry(row.n).or_default();
        *i += 1;
        let (sizes, pat) = describe(&row.ty);
        let a = &row.a;
        writeln!(
            csv,
            "{},{},{},{},{},{},{},{},{},{},{}",
            row.n,
            i,
            sizes,
            pat,
            a.r,
            a.m,
            a.genus,
            a.period,
            a.hyper_via_g,
            a.genus == 3,
            a.genus == row.n as i64
        )
        .unwrap();
    }
    fs::write(out.join("types.csv"), csv).unwrap();

    // genus_by_n: (n, period, genus) aggregates, as plotted.
    let mut agg: BTreeMap<(u32, u32, i64, u32, u32), (u32, u32)> = BTreeMap::new();
    for row in &rows {
        let e = agg
            .entry((row.n, row.a.period, row.a.genus, row.a.m, row.a.r))
            .or_default();
        e.0 += 1;
        if row.a.hyper_via_g {
            e.1 += 1;
        }
    }
    let mut csv = String::from("n,period,genus,m,r,types,hyperelliptic_via_G_types\n");
    for ((n, per, g, m, r), (c, hy)) in &agg {
        writeln!(csv, "{n},{per},{g},{m},{r},{c},{hy}").unwrap();
    }
    fs::write(out.join("genus_by_n.csv"), csv).unwrap();
    let min_table = min_genus_table(NMAX_MIN);
    if min_table.contains("FAIL") {
        failures += 1;
    }
    fs::write(out.join("min_genus.csv"), min_table).unwrap();

    // counts: exact branch-set counts per (q, n).
    let mut csv = String::from(
        "p,d,q,n,total_branch_sets,check_total,check_sigma_stable,min_genus_period_n,\
         genus3_branch_sets,genus3_hyperelliptic_via_G,genus3_not_hyperelliptic_via_G,\
         log10_share_genus3\n",
    );
    for &p in &PRIMES {
        for d in 1..=2u32 {
            let q = p.pow(d);
            for n in 1..=NMAX {
                let mut total = Big::from_u64(0);
                let mut stable = Big::from_u64(0);
                let mut g3_h = Big::from_u64(0);
                let mut g3_n = Big::from_u64(0);
                let mut min_g: Option<i64> = None;
                for row in rows.iter().filter(|r| r.n == n) {
                    let c = count_type(q, &row.ty);
                    total = total.add(&c);
                    if row.a.period == 1 {
                        stable = stable.add(&c);
                    }
                    if row.a.genus == 3 {
                        if row.a.hyper_via_g {
                            g3_h = g3_h.add(&c);
                        } else {
                            g3_n = g3_n.add(&c);
                        }
                    }
                    if row.a.period == n && !c.is_zero() {
                        min_g = Some(min_g.map_or(row.a.genus, |g: i64| g.min(row.a.genus)));
                    }
                }
                let big_q = Big::pow_small(q, n);
                let q2 = big_q.mul(&big_q);
                let expect = q2.mul(&q2).sub(&q2);
                let ok_total = total == expect;
                let small_q2 = Big::pow_small(q, 2);
                let ok_stable = stable == small_q2.mul(&small_q2).sub(&small_q2);
                if !ok_total || !ok_stable {
                    failures += 1;
                }
                let g3 = g3_h.add(&g3_n);
                let share = g3.log10() - total.log10();
                writeln!(
                    csv,
                    "{p},{d},{q},{n},{total},{},{},{},{g3},{g3_h},{g3_n},{}",
                    if ok_total { "ok" } else { "FAIL" },
                    if ok_stable { "ok" } else { "FAIL" },
                    min_g.map_or("none".to_string(), |g| g.to_string()),
                    if g3.is_zero() {
                        "-inf".to_string()
                    } else {
                        format!("{share:.4}")
                    }
                )
                .unwrap();
            }
        }
    }
    fs::write(out.join("counts.csv"), csv).unwrap();

    // ---- Part B: Joux–Vitse family over F_p.
    let mut csv = String::from(
        "p,arm,h_coeffs_asc,k3_modulus_asc,alpha,num_E_Fp3,t,E_Fp3_2torsion,deg_H,\
         H_coeffs_asc,N1,N2,N3,L_check,jac_order_L1,jac_eq_numE,weierstrass_cycle_type,\
         stable_pair_partitions,smith_lemma1,lemma1_match\n",
    );
    for &p in &PRIMES {
        let k1 = Gf::new(p, 1);
        let k2 = Gf::new(p, 2);
        let k3 = Gf::new(p, 3);
        let alphas = [p as usize, p as usize + 1, (p * p + p + 2) as usize];
        for &ai in &alphas {
            let alpha = k3.el(ai);
            for arm in [Arm::Linear, Arm::Split, Arm::Conjugate] {
                let r = jv_instance(&k3, &k2, &k1, alpha, arm);
                let jac_eq = r.jac_order == r.n_e as i64;
                let match1 = r.stable.is_some() && r.stable == r.lemma1;
                if !r.l_ok || !jac_eq || !match1 {
                    failures += 1;
                }
                writeln!(
                    csv,
                    "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                    r.p,
                    r.arm.name(),
                    join(&r.h, " "),
                    join(&r.modulus, " "),
                    join(&r.alpha[..3], " "),
                    r.n_e,
                    r.t,
                    r.e_k_2torsion,
                    r.f_h.len() - 1,
                    join(&r.f_h, " "),
                    r.n[0],
                    r.n[1],
                    r.n[2],
                    if r.l_ok { "ok" } else { "FAIL" },
                    r.jac_order,
                    jac_eq,
                    r.cycle_type
                        .as_ref()
                        .map_or("not_squarefree".to_string(), |c| join(c, " ")),
                    r.stable.map_or("na".to_string(), |s| s.to_string()),
                    r.lemma1.map_or("na".to_string(), |s| s.to_string()),
                    match1
                )
                .unwrap();
            }
        }
    }
    fs::write(out.join("jv_family.csv"), csv).unwrap();

    println!(
        "types: {} over n = 1..{NMAX}; failures: {failures}",
        rows.len()
    );
    if failures > 0 {
        eprintln!("FAILED checks: {failures}");
        std::process::exit(1);
    }
}

// ===========================================================================
// Tests against values the sources state.

#[cfg(test)]
mod tests {
    use super::*;

    fn rows_for(n: u32) -> Vec<TypeRow> {
        types(n)
            .into_iter()
            .map(|ty| {
                let a = analyse(n, &ty);
                TypeRow { n, ty, a }
            })
            .collect()
    }

    /// Smith 2009, Lemma 1, against brute force over every cycle type of S_8.
    #[test]
    fn smith_lemma1_table_matches_bruteforce() {
        fn partitions(rem: u32, max: u32, cur: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
            if rem == 0 {
                out.push(cur.clone());
                return;
            }
            for v in (1..=rem.min(max)).rev() {
                cur.push(v);
                partitions(rem - v, v, cur, out);
                cur.pop();
            }
        }
        let mut all = Vec::new();
        partitions(8, 8, &mut Vec::new(), &mut all);
        assert_eq!(all.len(), 22);
        for ct in all {
            assert_eq!(stable_matchings(&ct), smith_lemma1(&ct), "{ct:?}");
        }
    }

    /// Momose–Chao (ECC 2007, slides 24–25): over a cubic extension the
    /// genus-3 covers are one hyperelliptic shape of covering degree 8
    /// (`y^2 = e g(x)(x−α)(x−α^q)`, deg g = 1 or 2) and the non-hyperelliptic
    /// Types 1 and 2 of covering degree 4.
    #[test]
    fn n3_genus3_types_are_the_published_shapes() {
        let mut g3: Vec<(String, String, u32, bool)> = rows_for(3)
            .into_iter()
            .filter(|r| r.a.genus == 3)
            .map(|r| {
                let (s, p) = describe(&r.ty);
                (s, p, r.a.m, r.a.hyper_via_g)
            })
            .collect();
        g3.sort();
        let expect = vec![
            (
                "1+1+3".to_string(),
                "1:{0};1:{0};3:{0,1}".to_string(),
                3,
                true,
            ),
            ("2+3".to_string(), "2:{0,1};3:{0,1}".to_string(), 3, true),
            ("3+3".to_string(), "3:{0,1};3:{0,1}".to_string(), 2, false),
            ("6".to_string(), "6:{0,1,3,4}".to_string(), 2, false),
        ];
        assert_eq!(g3, expect);
    }

    /// Momose–Chao (ECC 2007, slides 20–21): with `g = n` (their Condition (C),
    /// `g = d·g0` for `g0 = 1`) the odd-characteristic cases are
    /// `(d, n) = (2,2)`, `(3,2)`, `(3,3)`, `(5,4)` and `d = 2^n − 1`.
    /// Here `d` is our `n` and their `n` is our `m`.
    #[test]
    fn minimal_covers_match_momose_chao() {
        let mut found = BTreeSet::new();
        for n in 1..=NMAX {
            for r in rows_for(n) {
                if r.a.period == n && r.a.genus == n as i64 && n >= 2 {
                    found.insert((n, r.a.m, r.a.hyper_via_g));
                }
            }
        }
        let expect: BTreeSet<(u32, u32, bool)> = [
            (2, 2, true),
            (3, 2, false),
            (3, 3, true),
            (5, 4, false),
            (7, 3, false),
        ]
        .into_iter()
        .collect();
        assert_eq!(found, expect);
    }

    /// Joux–Vitse, Appendix C: the two `n = 6` shapes (12) and (13) have genus 9.
    #[test]
    fn jv_appendix_c_genus9() {
        let rows = rows_for(6);
        let find = |sizes: &str, pat: &str| {
            rows.iter()
                .find(|r| describe(&r.ty) == (sizes.to_string(), pat.to_string()))
                .map(|r| r.a.genus)
        };
        // (12): {α, σα, σ²α, σ³α}, α of degree 6.
        assert_eq!(find("6", "6:{0,1,2,3}"), Some(9));
        // (13): {∞, β, α, σα}, β of degree 2, α of degree 3.
        assert_eq!(find("1+2+3", "1:{0};2:{0};3:{0,1}"), Some(9));
    }

    /// Diem, Compositio 2011, Prop. 4.10(a): `g <= (2^n − 1)(2n − 1)`.
    #[test]
    fn diem_upper_bound() {
        for n in 1..=NMAX {
            let bound = ((1i64 << n) - 1) * (2 * n as i64 - 1);
            for r in rows_for(n) {
                assert!(r.a.genus <= bound.max(1), "n={n} g={}", r.a.genus);
            }
        }
    }

    /// Derived here: genus 3 forces `(m, r) ∈ {(2,6), (3,5)}` and period 2 or 3;
    /// period >= 4 forces genus >= 5. The generic genus is `(n−1)2^n + 1`.
    #[test]
    fn genus3_needs_period_two_or_three() {
        for n in 1..=NMAX {
            let rows = rows_for(n);
            for r in &rows {
                if r.a.genus == 3 {
                    assert!(matches!((r.a.m, r.a.r), (2, 6) | (3, 5)));
                    assert!(r.a.period == 2 || r.a.period == 3);
                }
                if r.a.period >= 4 {
                    assert!(r.a.genus >= 5);
                }
            }
            let gmax = rows.iter().map(|r| r.a.genus).max().unwrap();
            assert_eq!(gmax, (n as i64 - 1) * (1i64 << n) + 1);
        }
    }

    /// Derived here, and Momose–Chao slide 20 ("2|d, char(k) ≠ 2: (d,n) = (2,2)"):
    /// for even n >= 6 no branch type of any period has genus exactly n
    /// (g − 1 odd forces m = 2, hence r = n + 3 <= 8). At n = 4 genus 4
    /// occurs only on period-2 types (m = 2, r = 7).
    #[test]
    fn no_genus_n_cover_for_even_n() {
        for n in [6u32, 8, 10] {
            for ty in types(n) {
                assert_ne!(analyse(n, &ty).genus, n as i64, "n={n}");
            }
        }
        for ty in types(4) {
            let a = analyse(4, &ty);
            if a.genus == 4 {
                assert_eq!((a.period, a.m, a.r), (2, 2, 7));
            }
        }
    }

    /// Repository census at n = 5 (research/ghs_cover_census_n5_20261013):
    /// `m ∈ {1, 4, 5}` and no genus in 6..=16.
    #[test]
    fn n5_census_statements() {
        for r in rows_for(5) {
            assert!(matches!(r.a.m, 1 | 4 | 5), "m={}", r.a.m);
            assert!(!(6..=16).contains(&r.a.genus), "g={}", r.a.genus);
        }
    }

    #[test]
    fn counts_are_complete() {
        for q in [5u64, 7, 9] {
            for n in 1..=4 {
                let total = rows_for(n)
                    .iter()
                    .fold(Big::from_u64(0), |a, r| a.add(&count_type(q, &r.ty)));
                let big_q = Big::pow_small(q, n);
                let q2 = big_q.mul(&big_q);
                assert_eq!(total, q2.mul(&q2).sub(&q2), "q={q} n={n}");
            }
        }
    }

    /// The repository's explicit F_7 instance (research/ghs_cubic_hyperelliptic_genus3_20261013):
    /// K = F_7[t]/(t^3 + 3), E0: y^2 = x(x−1)(x−t)(x−2t), #E0(K) = 328,
    /// t_E = 16, N1, N2, N3 = 8, 50, 296.
    #[test]
    fn jv_model_matches_repository_f7_instance() {
        let k3 = Gf::with_modulus(7, vec![3, 0, 0, 1]);
        let (k1, k2) = (Gf::new(7, 1), Gf::new(7, 2));
        let r = jv_instance(&k3, &k2, &k1, [0, 1, 0], Arm::Split);
        assert_eq!(r.n_e, 328);
        assert_eq!(r.t, 16);
        assert_eq!(r.n, [8, 50, 296]);
        assert!(r.l_ok);
        assert_eq!(r.jac_order, 328);
    }

    /// Joux–Vitse family: L-polynomial and Weierstrass cycle types.
    #[test]
    fn jv_family_small_primes() {
        for p in [5u64, 7, 11] {
            let (k1, k2, k3) = (Gf::new(p, 1), Gf::new(p, 2), Gf::new(p, 3));
            for arm in [Arm::Linear, Arm::Split, Arm::Conjugate] {
                let r = jv_instance(&k3, &k2, &k1, k3.el(p as usize), arm);
                assert!(r.l_ok, "p={p} {arm:?}");
                let ct = r.cycle_type.clone().expect("squarefree");
                let expect = if arm == Arm::Conjugate {
                    vec![2, 6]
                } else {
                    vec![1, 1, 3, 3]
                };
                assert_eq!(ct, expect, "p={p} {arm:?}");
            }
        }
    }
}
