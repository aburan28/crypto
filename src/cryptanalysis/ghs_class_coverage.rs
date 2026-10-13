//! Isogeny-class coverage of classical binary GHS Weil descent.
//!
//! The GHS magic number `m(b)` of `E_b: y^2 + xy = x^3 + ax^2 + b` over
//! `K = F_{2^N}`, relative to `N = n*l`, is not an isogeny invariant, while
//! the ECDLP transfers along any explicit isogeny of degree coprime to the
//! subgroup order.  By Tate, two ordinary curves over `K` are `K`-isogenous
//! exactly when their traces agree, so the set a GHS attacker can reach
//! (given the isogeny) is the union of the trace classes that contain one
//! small-`m` curve.  This module measures, exhaustively for small `N`:
//!
//! * the **direct** fraction of curves with `2 <= m <= M`;
//! * the **class** fraction: curves whose trace class contains a curve with
//!   `2 <= m <= M` (Tate: isogenous over `K`);
//! * the **3-walk** fraction: curves joined to such a curve by a path of
//!   `3`-isogenies and Frobenius conjugations (`j -> j^2`), with the path
//!   length.  This is the part of the class that a small-degree walk
//!   reaches without leaving the `3`-isogeny volcano component; the rest of
//!   the class needs other primes `l` or a vertical/horizontal class-group
//!   walk, which is not executed here.
//!
//! ## Exact definitions used
//!
//! * Magic number (Maurer–Menezes–Teske, LMS JCM 5 (2002), eq. (5)):
//!   `m(b) = dim_F2 Span{(1, sqrt(b_i)) : 0 <= i < n}`, `b_i = sigma^i(b)`,
//!   `sigma = x -> x^(2^l)`.  The leading `1` matters: it adds one to the
//!   plain orbit dimension exactly when `(x - 1)` does not divide `Ord_b`
//!   (MMT Theorem 5(i)); [`tests::plain_orbit_undercounts`] records the
//!   difference.  [`crate::cryptanalysis::ec_trapdoor::magic_number`] uses
//!   the same MMT definition ([`tests::mmt_matches_ec_trapdoor`]).
//! * Applicability to the twist `Tr_{K/F2}(a) = 1`: MMT Lemma 6, condition
//!   (8): `n` odd, or `j_0 = 2^e`.  With `P(x) = (x^n - 1)/(x - 1)`,
//!   `j_0 = 2^e` iff `P(sigma) b != 0`, and `P(sigma) b = Tr_{K/k}(b)`, so
//!   the twist is admitted iff `n` is odd or `Tr_{K/k}(b) != 0`.  A twist
//!   failing (8) is counted as *not* covered (GHS Theorem 1 does not apply);
//!   that is conservative.
//! * Traces: for `x != 0`, `y = xz` gives `z^2 + z = x + a + b/x^2`, so
//!   `#E = 2^N + 1 + (-1)^{Tr a} W(sqrt b)` with the Kloosterman sum
//!   `W(c) = sum_{u != 0} (-1)^{Tr(1/u) + Tr(c u)}`, i.e. `t = -(-1)^{Tr a} W`.
//!   All `W(c)` come from one fast Walsh–Hadamard transform of `Tr(1/u)`.
//!   [`tests::traces_match_naive_count`] checks this against direct counts.
//!
//! The genus of the GHS curve is `2^{m-1}` or `2^{m-1} - 1` (GHS Theorem 1);
//! the threshold here is on `m`.  Curves with `m = 1` are subfield curves
//! (`b` in `k = F_{2^l}`); MMT exclude them and so does this module: a weak
//! source has `2 <= m <= M`.
//!
//! Scope: these are exact structural counts over toy fields.  They say
//! nothing about descent cost, Jacobian order, relation collection, or
//! transfer of a prime-order subgroup.  Class membership is an existence
//! statement; building the isogeny is the cost the 3-walk column partly
//! measures.

use crate::cryptanalysis::binary_isogeny::phi_l_mod2_monomials;
use serde::Serialize;
use std::collections::VecDeque;

/// Largest supported field degree (log tables of `2 * 2^N` `u32`s).
pub const MAX_N: u32 = 26;

/// `F_{2^N}` with discrete-log tables, `N <= MAX_N`.
pub struct Gf2Log {
    pub n: u32,
    /// Monic modulus including the `z^N` bit.
    pub modulus: u64,
    q1: u32,
    exp: Vec<u32>,
    log: Vec<u32>,
    tr_mask: u32,
}

fn clmul_mod(mut a: u64, mut b: u64, modulus: u64, n: u32) -> u64 {
    let mut r = 0u64;
    while b != 0 {
        if b & 1 == 1 {
            r ^= a;
        }
        b >>= 1;
        a <<= 1;
        if (a >> n) & 1 == 1 {
            a ^= modulus;
        }
    }
    r
}

fn pow_mod(mut a: u64, mut e: u64, modulus: u64, n: u32) -> u64 {
    let mut r = 1u64;
    while e > 0 {
        if e & 1 == 1 {
            r = clmul_mod(r, a, modulus, n);
        }
        a = clmul_mod(a, a, modulus, n);
        e >>= 1;
    }
    r
}

/// Binary polynomial gcd on bit-packed `u64`s.
fn gf2x_gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let db = 63 - b.leading_zeros();
        while a != 0 && 63 - a.leading_zeros() >= db {
            a ^= b << (63 - a.leading_zeros() - db);
        }
        std::mem::swap(&mut a, &mut b);
    }
    a
}

/// `x^(2^k) mod f`, with `x^(2^k)` reduced by repeated squaring.
fn x_pow_2k(k: u32, modulus: u64, n: u32) -> u64 {
    let mut r = 2u64; // x
    if n == 1 {
        r = clmul_mod(1, 2 ^ modulus, modulus, n);
    }
    for _ in 0..k {
        r = clmul_mod(r, r, modulus, n);
    }
    r
}

/// Ben-Or irreducibility test for a degree-`n` binary polynomial.
fn is_irreducible(modulus: u64, n: u32) -> bool {
    for i in 1..=n / 2 {
        let t = x_pow_2k(i, modulus, n) ^ 2;
        if gf2x_gcd(modulus, t) != 1 {
            return false;
        }
    }
    true
}

fn prime_factors(mut v: u64) -> Vec<u64> {
    let mut out = Vec::new();
    let mut p = 2u64;
    while p * p <= v {
        if v % p == 0 {
            out.push(p);
            while v % p == 0 {
                v /= p;
            }
        }
        p += 1;
    }
    if v > 1 {
        out.push(v);
    }
    out
}

impl Gf2Log {
    /// Build the field with the lexicographically first irreducible
    /// trinomial or pentanomial and the smallest primitive element.
    pub fn new(n: u32) -> Self {
        assert!((2..=MAX_N).contains(&n), "N must be in 2..={MAX_N}");
        let top = 1u64 << n;
        let mut modulus = 0u64;
        'search: for weight_two in [true, false] {
            if weight_two {
                for a in 1..n {
                    let f = top | (1 << a) | 1;
                    if is_irreducible(f, n) {
                        modulus = f;
                        break 'search;
                    }
                }
            } else {
                for a in 3..n {
                    for b in 2..a {
                        for c in 1..b {
                            let f = top | (1 << a) | (1 << b) | (1 << c) | 1;
                            if is_irreducible(f, n) {
                                modulus = f;
                                break 'search;
                            }
                        }
                    }
                }
            }
        }
        assert!(modulus != 0, "no irreducible tri/pentanomial of degree {n}");
        let q1 = (1u64 << n) - 1;
        let factors = prime_factors(q1);
        let mut g = 2u64;
        loop {
            if factors.iter().all(|p| pow_mod(g, q1 / p, modulus, n) != 1) {
                break;
            }
            g += 1;
        }
        let q1u = q1 as u32;
        let mut exp = vec![0u32; 2 * q1 as usize];
        let mut log = vec![u32::MAX; 1usize << n];
        let mut x = 1u64;
        for i in 0..q1 as usize {
            exp[i] = x as u32;
            log[x as usize] = i as u32;
            x = clmul_mod(x, g, modulus, n);
        }
        for i in 0..q1 as usize {
            exp[q1 as usize + i] = exp[i];
        }
        let mut f = Gf2Log { n, modulus, q1: q1u, exp, log, tr_mask: 0 };
        let mut mask = 0u32;
        for i in 0..n {
            if f.trace_slow(1 << i) == 1 {
                mask |= 1 << i;
            }
        }
        f.tr_mask = mask;
        f
    }

    pub fn size(&self) -> usize {
        1usize << self.n
    }

    #[inline]
    pub fn mul(&self, a: u32, b: u32) -> u32 {
        if a == 0 || b == 0 {
            return 0;
        }
        self.exp[(self.log[a as usize] + self.log[b as usize]) as usize]
    }

    #[inline]
    pub fn inv(&self, a: u32) -> u32 {
        debug_assert!(a != 0);
        let la = self.log[a as usize];
        self.exp[((self.q1 - la) % self.q1) as usize]
    }

    /// `a^(2^k)` (the `k`-th power of absolute Frobenius).
    #[inline]
    pub fn frob(&self, a: u32, k: u32) -> u32 {
        if a == 0 {
            return 0;
        }
        let la = self.log[a as usize] as u64;
        let e = (la << (k % self.n)) % self.q1 as u64;
        self.exp[e as usize]
    }

    #[inline]
    pub fn sqrt(&self, a: u32) -> u32 {
        self.frob(a, self.n - 1)
    }

    fn trace_slow(&self, a: u32) -> u32 {
        let mut s = 0u32;
        let mut x = a;
        for _ in 0..self.n {
            s ^= x;
            x = self.mul(x, x);
        }
        debug_assert!(s <= 1);
        s
    }

    /// Absolute trace `Tr_{K/F2}`.
    #[inline]
    pub fn trace(&self, a: u32) -> u32 {
        (a & self.tr_mask).count_ones() & 1
    }

    /// Relative trace `Tr_{K/k}(a) = sum_{i<n} a^(2^(l i))`.
    pub fn rel_trace(&self, a: u32, l: u32) -> u32 {
        let n = self.n / l;
        let mut s = 0u32;
        for i in 0..n {
            s ^= self.frob(a, l * i);
        }
        s
    }
}

/// Traces `t(b)` of `y^2 + xy = x^3 + b` (the `a = 0` family) for every
/// `b`, index `b`; entry `0` is unused.  The twist (`Tr a = 1`) has `-t(b)`.
pub fn traces_a0(f: &Gf2Log) -> Vec<i32> {
    let q = f.size();
    // F[u] = (-1)^{Tr(1/u)} for u != 0; F[0] = +1 (subtracted below).
    let mut h: Vec<i32> = vec![1; q];
    for u in 1..q as u32 {
        if f.trace(f.inv(u)) == 1 {
            h[u as usize] = -1;
        }
    }
    let mut len = 1;
    while len < q {
        for i in (0..q).step_by(2 * len) {
            for j in i..i + len {
                let (x, y) = (h[j], h[j + len]);
                h[j] = x + y;
                h[j + len] = x - y;
            }
        }
        len <<= 1;
    }
    // v(c): bit i is Tr(c z^i), so that v(c).u = Tr(c u) for all u.
    let n = f.n;
    let w: Vec<u32> = (0..n)
        .map(|j| {
            let mut v = 0u32;
            for i in 0..n {
                if f.trace(f.mul(1 << j, 1 << i)) == 1 {
                    v |= 1 << i;
                }
            }
            v
        })
        .collect();
    let mut vc = vec![0u32; q];
    for c in 1..q {
        let low = c.trailing_zeros();
        vc[c] = vc[c & (c - 1)] ^ w[low as usize];
    }
    let mut t = vec![0i32; q];
    for b in 1..q as u32 {
        let c = f.sqrt(b);
        let kloosterman = h[vc[c as usize] as usize] - 1;
        t[b as usize] = -kloosterman;
    }
    t
}

/// Direct point count of `y^2 + xy = x^3 + ax^2 + b` (test reference).
pub fn naive_count(f: &Gf2Log, a: u32, b: u32) -> u64 {
    let mut count = 1u64; // infinity
    for x in 0..f.size() as u32 {
        let x2 = f.mul(x, x);
        let rhs = f.mul(x2, x) ^ f.mul(a, x2) ^ b;
        for y in 0..f.size() as u32 {
            if f.mul(y, y) ^ f.mul(x, y) == rhs {
                count += 1;
            }
        }
    }
    count
}

/// MMT magic number `dim Span{(1, sqrt(sigma^i b))}` relative to `n = N/l`.
pub fn magic_mmt(f: &Gf2Log, b: u32, l: u32) -> u32 {
    let n = f.n / l;
    let s = f.sqrt(b);
    let one = 1u64 << f.n;
    let rows: Vec<u64> = (0..n).map(|i| one | f.frob(s, l * i) as u64).collect();
    rank_u64(rows)
}

/// Plain orbit dimension `dim Span{sqrt(sigma^i b)}` (the definition
/// `ec_trapdoor::magic_number` used before it was corrected to MMT's).
pub fn magic_plain(f: &Gf2Log, b: u32, l: u32) -> u32 {
    let n = f.n / l;
    let s = f.sqrt(b);
    rank_u64((0..n).map(|i| f.frob(s, l * i) as u64).collect())
}

fn rank_u64(mut rows: Vec<u64>) -> u32 {
    let mut rank = 0u32;
    for i in 0..rows.len() {
        let r = rows[i];
        if r == 0 {
            continue;
        }
        rank += 1;
        let pivot = 63 - r.leading_zeros();
        for row in rows.iter_mut().skip(i + 1) {
            if (*row >> pivot) & 1 == 1 {
                *row ^= r;
            }
        }
    }
    rank
}

/// MMT Lemma 6 condition (8) for the `Tr a = 1` twist.
pub fn twist_admitted(f: &Gf2Log, b: u32, l: u32) -> bool {
    (f.n / l) % 2 == 1 || f.rel_trace(b, l) != 0
}

// ---------------------------------------------------------------------------
// Polynomials over F_{2^N} (small degree) for Phi_3 root finding.

fn ptrim(p: &mut Vec<u32>) {
    while p.last() == Some(&0) {
        p.pop();
    }
}

fn pmod(f: &Gf2Log, a: &[u32], m: &[u32]) -> Vec<u32> {
    let mut r = a.to_vec();
    ptrim(&mut r);
    let dm = m.len() - 1;
    let lead_inv = f.inv(m[dm]);
    while r.len() > dm {
        let c = f.mul(*r.last().unwrap(), lead_inv);
        let shift = r.len() - 1 - dm;
        for (i, &mi) in m.iter().enumerate() {
            r[shift + i] ^= f.mul(c, mi);
        }
        ptrim(&mut r);
    }
    r
}

fn pmulmod(f: &Gf2Log, a: &[u32], b: &[u32], m: &[u32]) -> Vec<u32> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut r = vec![0u32; a.len() + b.len() - 1];
    for (i, &ai) in a.iter().enumerate() {
        if ai == 0 {
            continue;
        }
        for (j, &bj) in b.iter().enumerate() {
            r[i + j] ^= f.mul(ai, bj);
        }
    }
    pmod(f, &r, m)
}

fn pgcd(f: &Gf2Log, a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    ptrim(&mut a);
    ptrim(&mut b);
    while !b.is_empty() {
        let r = pmod(f, &a, &b);
        a = b;
        b = r;
    }
    if let Some(&lead) = a.last() {
        let li = f.inv(lead);
        for c in a.iter_mut() {
            *c = f.mul(*c, li);
        }
    }
    a
}

fn pdiv_exact(f: &Gf2Log, a: &[u32], m: &[u32]) -> Vec<u32> {
    let mut r = a.to_vec();
    ptrim(&mut r);
    let dm = m.len() - 1;
    let lead_inv = f.inv(m[dm]);
    let mut qt = vec![0u32; r.len().saturating_sub(dm)];
    while r.len() > dm {
        let c = f.mul(*r.last().unwrap(), lead_inv);
        let shift = r.len() - 1 - dm;
        qt[shift] = c;
        for (i, &mi) in m.iter().enumerate() {
            r[shift + i] ^= f.mul(c, mi);
        }
        ptrim(&mut r);
    }
    qt
}

/// Distinct roots in `F_{2^N}` of a nonzero polynomial (ascending
/// coefficients).
pub fn roots(f: &Gf2Log, poly: &[u32]) -> Vec<u32> {
    let mut g = poly.to_vec();
    ptrim(&mut g);
    if g.len() <= 1 {
        return Vec::new();
    }
    // h = gcd(g, X^q - X): the product of distinct rational linear factors.
    let mut r = vec![0u32, 1u32];
    r = pmod(f, &r, &g);
    for _ in 0..f.n {
        r = pmulmod(f, &r, &r, &g);
    }
    if r.len() < 2 {
        r.resize(2, 0);
    }
    r[1] ^= 1;
    ptrim(&mut r);
    let h = if r.is_empty() { pgcd(f, &g, &[]) } else { pgcd(f, &g, &r) };
    let mut out = Vec::new();
    split_roots(f, h, &mut out, 1);
    out.sort_unstable();
    out
}

fn split_roots(f: &Gf2Log, h: Vec<u32>, out: &mut Vec<u32>, mut seed: u32) {
    match h.len() {
        0 | 1 => {}
        2 => out.push(f.mul(h[0], f.inv(h[1]))),
        _ => loop {
            // Tr(delta X) mod h for delta = generator powers.
            let delta = f.exp[(seed as usize) % f.q1 as usize];
            seed += 1;
            let mut x = pmod(f, &[0, delta], &h);
            let mut acc = x.clone();
            for _ in 1..f.n {
                x = pmulmod(f, &x, &x, &h);
                let len = acc.len().max(x.len());
                acc.resize(len, 0);
                for (i, &c) in x.iter().enumerate() {
                    acc[i] ^= c;
                }
                ptrim(&mut acc);
            }
            let d = pgcd(f, &h, &acc);
            if d.len() > 1 && d.len() < h.len() {
                let e = pdiv_exact(f, &h, &d);
                split_roots(f, d, out, seed);
                split_roots(f, e, out, seed + 7);
                return;
            }
        },
    }
}

/// `Phi_l(X, j) mod 2` as a polynomial in `X` (ascending), `l in {2, 3}`.
pub fn phi_mod2_at(f: &Gf2Log, monos: &[(u32, u32)], j: u32) -> Vec<u32> {
    let max_x = monos.iter().map(|m| m.0).max().unwrap_or(0) as usize;
    let max_y = monos.iter().map(|m| m.1).max().unwrap_or(0) as usize;
    let mut jp = vec![1u32; max_y + 1];
    for k in 1..=max_y {
        jp[k] = f.mul(jp[k - 1], j);
    }
    let mut p = vec![0u32; max_x + 1];
    for &(i, k) in monos {
        p[i as usize] ^= jp[k as usize];
    }
    ptrim(&mut p);
    p
}

// ---------------------------------------------------------------------------
// Coverage experiment.

#[derive(Clone, Debug, Serialize)]
pub struct ThresholdRow {
    /// Threshold: a curve is "weak" when `2 <= m <= max_m` (and condition
    /// (8) holds for the twist family).  Subfield curves (`m = 1`) are never
    /// weak sources: the trace kills the large subgroup (MMT exclude them).
    pub max_m: u32,
    pub curves_total: u64,
    pub classes_total: u64,
    pub direct_weak_curves: u64,
    pub classes_with_weak_member: u64,
    /// Curves whose `K`-isogeny class (trace) contains a weak curve.
    pub class_covered_curves: u64,
    /// `None` when the 3-isogeny graph was not built for this `N`.
    pub walk3_covered_curves: Option<u64>,
    /// Histogram of 0-1 BFS distance (3-isogeny steps; Frobenius free) to
    /// the nearest weak curve, over curves in weak classes.  Index = steps.
    pub walk3_distance_histogram: Option<Vec<u64>>,
    /// Median over weak classes of `log2(sqrt(class size))`, a proxy for a
    /// generic meet-in-the-middle class walk (curves, not group operations).
    pub median_log2_sqrt_weak_class_size: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct FactorRow {
    pub big_n: u32,
    pub n: u32,
    pub l: u32,
    /// Number of `b in K*` (a = 0 family) with each MMT magic number.
    pub magic_histogram_a0: Vec<u64>,
    /// `b` with `magic_plain + 1 == magic_mmt` (plain orbit undercounts).
    pub plain_undercount_b: u64,
    /// Subfield curves (`m = 1`, both families): excluded as weak sources.
    pub subfield_curves: u64,
    pub twist_failing_condition8: u64,
    pub thresholds: Vec<ThresholdRow>,
}

#[derive(Clone, Debug, Serialize)]
pub struct FieldReport {
    pub big_n: u32,
    pub modulus_hex: String,
    pub curves_total: u64,
    pub classes_total: u64,
    pub walk3_built: bool,
    pub walk3_edges: Option<u64>,
    pub factorisations: Vec<FactorRow>,
}

/// Union-free 3-isogeny + Frobenius graph on the `2(q-1)` ordinary curves
/// `(b, family)`; vertex id `2 b + family`.
pub struct Walk3Graph {
    pub adj3: Vec<Vec<u32>>,
    pub edges: u64,
}

pub fn build_walk3(f: &Gf2Log, t: &[i32]) -> Walk3Graph {
    let monos = phi_l_mod2_monomials(3).expect("Phi_3 tabled");
    let q = f.size();
    let mut adj3: Vec<Vec<u32>> = vec![Vec::new(); 2 * q];
    let mut edges = 0u64;
    for b in 1..q as u32 {
        let j = f.inv(b);
        let poly = phi_mod2_at(f, &monos, j);
        for jp in roots(f, &poly) {
            if jp == 0 {
                continue;
            }
            let bp = f.inv(jp);
            let same = t[bp as usize] == t[b as usize];
            debug_assert!(same || t[bp as usize] == -t[b as usize]);
            for fam in 0..2u32 {
                let other = if same { fam } else { 1 - fam };
                adj3[(2 * b + fam) as usize].push(2 * bp + other);
                edges += 1;
            }
        }
    }
    Walk3Graph { adj3, edges }
}

fn frob_neighbours(f: &Gf2Log, v: u32) -> [u32; 2] {
    let (b, fam) = (v / 2, v % 2);
    [2 * f.frob(b, 1) + fam, 2 * f.sqrt(b) + fam]
}

/// 0-1 BFS from the weak vertices: 3-isogeny steps cost 1, Frobenius 0.
fn walk_distances(f: &Gf2Log, g: &Walk3Graph, weak: &[bool]) -> Vec<u32> {
    let mut dist = vec![u32::MAX; weak.len()];
    let mut dq = VecDeque::new();
    for (v, &w) in weak.iter().enumerate() {
        if w {
            dist[v] = 0;
            dq.push_back(v as u32);
        }
    }
    while let Some(v) = dq.pop_front() {
        let d = dist[v as usize];
        for u in frob_neighbours(f, v) {
            if dist[u as usize] > d {
                dist[u as usize] = d;
                dq.push_front(u);
            }
        }
        for &u in &g.adj3[v as usize] {
            if dist[u as usize] > d + 1 {
                dist[u as usize] = d + 1;
                dq.push_back(u);
            }
        }
    }
    dist
}

/// Run the exhaustive coverage experiment for `F_{2^N}`.
pub fn field_report(big_n: u32, build_walk: bool) -> FieldReport {
    let f = Gf2Log::new(big_n);
    let q = f.size();
    let t = traces_a0(&f);
    let bound = 2 * ((q as f64).sqrt().ceil() as i32) + 2;
    let off = bound;
    let width = (2 * bound + 1) as usize;
    let idx = |tr: i32| (tr + off) as usize;
    let mut class_size = vec![0u64; width];
    for b in 1..q {
        class_size[idx(t[b])] += 1;
        class_size[idx(-t[b])] += 1;
    }
    let classes_total = class_size.iter().filter(|&&c| c > 0).count() as u64;
    let curves_total = 2 * (q as u64 - 1);
    let graph = if build_walk { Some(build_walk3(&f, &t)) } else { None };

    let mut factorisations = Vec::new();
    for n in 2..=big_n {
        if big_n % n != 0 {
            continue;
        }
        let l = big_n / n;
        let mut mvals = vec![0u8; q];
        let mut adm = vec![false; q];
        let mut hist = vec![0u64; n as usize + 1];
        let mut undercount = 0u64;
        let mut fail8 = 0u64;
        for b in 1..q as u32 {
            let m = magic_mmt(&f, b, l);
            mvals[b as usize] = m as u8;
            hist[m as usize] += 1;
            if magic_plain(&f, b, l) + 1 == m {
                undercount += 1;
            }
            adm[b as usize] = twist_admitted(&f, b, l);
            if !adm[b as usize] {
                fail8 += 1;
            }
        }
        let subfield = 2 * hist[1];
        let mut thresholds = Vec::new();
        for max_m in 2..=n {
            let weak_v = |b: usize, fam: usize| {
                let m = mvals[b] as u32;
                (2..=max_m).contains(&m) && (fam == 0 || adm[b])
            };
            let mut weak_class = vec![false; width];
            let mut direct = 0u64;
            for b in 1..q {
                for fam in 0..2 {
                    if weak_v(b, fam) {
                        direct += 1;
                        let tr = if fam == 0 { t[b] } else { -t[b] };
                        weak_class[idx(tr)] = true;
                    }
                }
            }
            let classes_weak = weak_class.iter().filter(|&&w| w).count() as u64;
            let class_cov: u64 = (0..width).filter(|&i| weak_class[i]).map(|i| class_size[i]).sum();
            let mut sizes: Vec<f64> = (0..width)
                .filter(|&i| weak_class[i])
                .map(|i| (class_size[i] as f64).log2() / 2.0)
                .collect();
            sizes.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let median = if sizes.is_empty() { None } else { Some(sizes[sizes.len() / 2]) };
            let (walk_cov, walk_hist) = if let Some(g) = &graph {
                let mut weak = vec![false; 2 * q];
                for b in 1..q {
                    for fam in 0..2 {
                        weak[2 * b + fam] = weak_v(b, fam);
                    }
                }
                let dist = walk_distances(&f, g, &weak);
                let mut covered = 0u64;
                let mut hist_d: Vec<u64> = Vec::new();
                for b in 1..q {
                    for fam in 0..2 {
                        let tr = if fam == 0 { t[b] } else { -t[b] };
                        if !weak_class[idx(tr)] {
                            continue;
                        }
                        let d = dist[2 * b + fam];
                        if d != u32::MAX {
                            covered += 1;
                            if hist_d.len() <= d as usize {
                                hist_d.resize(d as usize + 1, 0);
                            }
                            hist_d[d as usize] += 1;
                        }
                    }
                }
                (Some(covered), Some(hist_d))
            } else {
                (None, None)
            };
            thresholds.push(ThresholdRow {
                max_m,
                curves_total,
                classes_total,
                direct_weak_curves: direct,
                classes_with_weak_member: classes_weak,
                class_covered_curves: class_cov,
                walk3_covered_curves: walk_cov,
                walk3_distance_histogram: walk_hist,
                median_log2_sqrt_weak_class_size: median,
            });
        }
        factorisations.push(FactorRow {
            big_n,
            n,
            l,
            magic_histogram_a0: hist,
            plain_undercount_b: undercount,
            subfield_curves: subfield,
            twist_failing_condition8: fail8,
            thresholds,
        });
    }
    FieldReport {
        big_n,
        modulus_hex: format!("{:#x}", f.modulus),
        curves_total,
        classes_total,
        walk3_built: graph.is_some(),
        walk3_edges: graph.as_ref().map(|g| g.edges),
        factorisations,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn traces_match_naive_count() {
        for n in [5u32, 6, 7, 8] {
            let f = Gf2Log::new(n);
            let t = traces_a0(&f);
            // a with Tr(a) = 1 for the twist.
            let a1 = (1..f.size() as u32).find(|&a| f.trace(a) == 1).unwrap();
            let q = f.size() as i64;
            for b in 1..f.size() as u32 {
                let c0 = naive_count(&f, 0, b) as i64;
                assert_eq!(c0, q + 1 - t[b as usize] as i64, "N={n} b={b}");
                assert_eq!(c0 % 4, 0, "a = 0 curves have #E = 0 mod 4");
                let c1 = naive_count(&f, a1, b) as i64;
                assert_eq!(c1, q + 1 + t[b as usize] as i64, "twist N={n} b={b}");
            }
        }
    }

    #[test]
    fn mmt_theorem5_counts() {
        // N = 12, n = 3 (q = 2^4): x^3 - 1 = (x - 1)(x^2 + x + 1), so
        // #{b != 0 : m = 1} = q - 1 and #{m = 3} = q^3 - q.
        let f = Gf2Log::new(12);
        let mut h = [0u64; 4];
        for b in 1..f.size() as u32 {
            h[magic_mmt(&f, b, 4) as usize] += 1;
        }
        assert_eq!(h, [0, 15, 0, 4096 - 16]);
        // n = 2 (q = 2^6): m = 1 on k*, m = 2 elsewhere.
        let mut h2 = [0u64; 3];
        for b in 1..f.size() as u32 {
            h2[magic_mmt(&f, b, 6) as usize] += 1;
        }
        assert_eq!(h2, [0, 63, 4096 - 64]);
    }

    #[test]
    fn plain_orbit_undercounts() {
        // MMT Theorem 5(i): m = plain + 1 exactly when (x - 1) does not
        // divide Ord_b, which for n = 3 is the b with Tr_{K/k} b = 0 and
        // b not in k.
        let f = Gf2Log::new(12);
        for b in 1..f.size() as u32 {
            let (mm, mp) = (magic_mmt(&f, b, 4), magic_plain(&f, b, 4));
            let expect_extra = f.rel_trace(b, 4) == 0;
            assert_eq!(mm, mp + expect_extra as u32, "b={b}");
        }
    }

    #[test]
    fn mmt_matches_ec_trapdoor() {
        use crate::binary_ecc::{F2mElement, IrreduciblePoly};
        use num_bigint::BigUint;
        use crate::cryptanalysis::ec_trapdoor::{magic_number, FieldTower};
        let f = Gf2Log::new(12);
        let irr = IrreduciblePoly {
            degree: 12,
            low_terms: (0..12).filter(|i| (f.modulus >> i) & 1 == 1).collect(),
        };
        let tower = FieldTower::new(12, 3, 4, irr);
        for b in (1..f.size() as u32).step_by(7) {
            let be = F2mElement::from_biguint(&BigUint::from(b), 12);
            let old = magic_number(&tower, &F2mElement::zero(12), &be);
            assert_eq!(old, magic_mmt(&f, b, 4), "b={b}");
        }
    }

    #[test]
    fn walk3_neighbours_are_isogenous() {
        let f = Gf2Log::new(10);
        let t = traces_a0(&f);
        let g = build_walk3(&f, &t);
        assert!(g.edges > 0);
        for b in 1..f.size() as u32 {
            for fam in 0..2u32 {
                let tv = if fam == 0 { t[b as usize] } else { -t[b as usize] };
                for &u in &g.adj3[(2 * b + fam) as usize] {
                    let (bu, fu) = (u / 2, u % 2);
                    let tu = if fu == 0 { t[bu as usize] } else { -t[bu as usize] };
                    assert_eq!(tu, tv, "3-isogeny preserves the trace");
                }
            }
        }
    }

    #[test]
    fn roots_finds_planted_roots() {
        let f = Gf2Log::new(9);
        // (X - 3)(X - 77)(X^2 + X + c) with c chosen so the quadratic is
        // irreducible (Tr(c) = 1).
        let c = (1..f.size() as u32).find(|&c| f.trace(c) == 1).unwrap();
        let lin = |r: u32| vec![r, 1u32];
        let mul = |a: &[u32], b: &[u32]| {
            let mut out = vec![0u32; a.len() + b.len() - 1];
            for (i, &x) in a.iter().enumerate() {
                for (j, &y) in b.iter().enumerate() {
                    out[i + j] ^= f.mul(x, y);
                }
            }
            out
        };
        let p = mul(&mul(&lin(3), &lin(77)), &[c, 1, 1]);
        assert_eq!(roots(&f, &p), vec![3, 77]);
    }
}
