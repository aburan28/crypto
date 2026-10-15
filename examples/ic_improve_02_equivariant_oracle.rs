//! IC improvement item 02: the G2 Frobenius-equivariant Macaulay pass as a
//! Koblitz decomposition oracle, measured against the incumbent oracle.
//!
//! Helpers (barrel system, 192-bit masks, F_{2^d}, orbits) are copied from
//! `examples/ic_next10_06_g2_full_pass.rs`.  New here: the gauge-fixed
//! system (the target-fixed system the pipeline solves), the incumbent step
//! on the crate kernel, specialisation of the equivariant rows to the gauge,
//! and the whole-oracle arms through the pipeline's public solver.
//!
//! See `research/ic_improve_20261014/02-equivariant-oracle/PROTOCOL.md`.
//!
//! ```bash
//! cargo run --release --example ic_improve_02_equivariant_oracle -- --mode step --n 7,11 --targets 8
//! cargo run --release --example ic_improve_02_equivariant_oracle -- --mode oracle --n 5,7 --targets 12
//! ```
#![allow(clippy::too_many_arguments, clippy::type_complexity)]

use std::collections::HashMap;
use std::time::Instant;

use crypto_lib::binary_ecc::f2m::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crypto_lib::cryptanalysis::gf2_elim;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    f4_word_ops_thread, solve_boolean_system, split_rule_default, SolveOptions, SolverEngine,
};
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};

const WORDS: usize = 3;
type Mask = [u64; WORDS];

fn mask_var(v: usize) -> Mask {
    let mut m = [0u64; WORDS];
    m[v / 64] |= 1u64 << (v % 64);
    m
}
fn mask_or(a: &Mask, b: &Mask) -> Mask {
    let mut m = [0u64; WORDS];
    for i in 0..WORDS {
        m[i] = a[i] | b[i];
    }
    m
}
fn mask_deg(a: &Mask) -> u32 {
    a.iter().map(|w| w.count_ones()).sum()
}

/// Boolean polynomial = sorted set of monomials (F_2 coefficients).
type Poly = Vec<Mask>;

fn poly_add(a: &Poly, b: &Poly) -> Poly {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                out.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                out.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}
fn poly_mul(a: &Poly, b: &Poly) -> Poly {
    let mut all: Vec<Mask> = Vec::with_capacity(a.len() * b.len());
    for x in a {
        for y in b {
            all.push(mask_or(x, y));
        }
    }
    normalize(all)
}
fn normalize(mut all: Vec<Mask>) -> Poly {
    all.sort_unstable();
    let mut out = Vec::with_capacity(all.len());
    let mut i = 0;
    while i < all.len() {
        let mut j = i;
        while j < all.len() && all[j] == all[i] {
            j += 1;
        }
        if (j - i) % 2 == 1 {
            out.push(all[i]);
        }
        i = j;
    }
    out
}
fn poly_mul_mono(p: &Poly, m: &Mask) -> Poly {
    normalize(p.iter().map(|t| mask_or(t, m)).collect())
}

/// Field structure constants for F_2[z]/(f): reduced[i][j] = z^{i+j} mod f, squares[k] = z^{2k} mod f.
struct Field {
    n: usize,
    reduced: Vec<Vec<u64>>,
    squares: Vec<u64>,
}
fn bits_of(e: &F2mElement) -> u64 {
    e.raw_bits().first().copied().unwrap_or(0)
}
fn elem_from_bits(bits: u64, n: u32) -> F2mElement {
    let positions: Vec<u32> = (0..n).filter(|k| (bits >> k) & 1 == 1).collect();
    F2mElement::from_bit_positions(&positions, n)
}
impl Field {
    fn new(n: u32, irr: &IrreduciblePoly) -> Self {
        let mono = |k: u32| F2mElement::from_bit_positions(&[k], n);
        let nn = n as usize;
        let mut reduced = vec![vec![0u64; nn]; nn];
        for i in 0..nn {
            for j in 0..nn {
                reduced[i][j] = bits_of(&mono(i as u32).mul(&mono(j as u32), irr));
            }
        }
        let squares = (0..nn).map(|k| reduced[k][k]).collect();
        Field {
            n: nn,
            reduced,
            squares,
        }
    }
}

/// Symbolic field element: n coordinate polynomials in the poly basis.
#[derive(Clone)]
struct Sym {
    coords: Vec<Poly>,
}
impl Sym {
    fn zero(n: usize) -> Self {
        Sym {
            coords: vec![Vec::new(); n],
        }
    }
    fn constant(v: u64, n: usize) -> Self {
        Sym {
            coords: (0..n)
                .map(|k| {
                    if (v >> k) & 1 == 1 {
                        vec![[0u64; WORDS]]
                    } else {
                        Vec::new()
                    }
                })
                .collect(),
        }
    }
    fn from_subspace_vars(basis: &[u64], offset: usize, n: usize) -> Self {
        let mut coords: Vec<Vec<Mask>> = vec![Vec::new(); n];
        for (t, &b) in basis.iter().enumerate() {
            for k in 0..n {
                if (b >> k) & 1 == 1 {
                    coords[k].push(mask_var(offset + t));
                }
            }
        }
        Sym {
            coords: coords.into_iter().map(normalize).collect(),
        }
    }
    fn add(&self, o: &Self) -> Self {
        Sym {
            coords: self
                .coords
                .iter()
                .zip(o.coords.iter())
                .map(|(a, b)| poly_add(a, b))
                .collect(),
        }
    }
    fn mul(&self, o: &Self, f: &Field) -> Self {
        let n = f.n;
        let mut acc: Vec<Vec<Mask>> = vec![Vec::new(); n];
        for i in 0..n {
            if self.coords[i].is_empty() {
                continue;
            }
            for j in 0..n {
                if o.coords[j].is_empty() {
                    continue;
                }
                let prod = poly_mul(&self.coords[i], &o.coords[j]);
                let red = f.reduced[i][j];
                for k in 0..n {
                    if (red >> k) & 1 == 1 {
                        acc[k].extend_from_slice(&prod);
                    }
                }
            }
        }
        Sym {
            coords: acc.into_iter().map(normalize).collect(),
        }
    }
    fn square(&self, f: &Field) -> Self {
        let n = f.n;
        let mut acc: Vec<Vec<Mask>> = vec![Vec::new(); n];
        for i in 0..n {
            let sq = f.squares[i];
            for k in 0..n {
                if (sq >> k) & 1 == 1 {
                    acc[k].extend_from_slice(&self.coords[i]);
                }
            }
        }
        Sym {
            coords: acc.into_iter().map(normalize).collect(),
        }
    }
    /// Multiply every coordinate by a Boolean variable.
    fn scale_var(&self, v: usize) -> Self {
        let m = mask_var(v);
        Sym {
            coords: self.coords.iter().map(|p| poly_mul_mono(p, &m)).collect(),
        }
    }
}

/// Binary S3 for y² + xy = x³ + ax² + b: (x1 + x2)² x3² + x1 x2 x3 + (x1 x2)² + b.
fn sym_s3(x1: &Sym, x2: &Sym, x3: &Sym, b: u64, f: &Field) -> Vec<Poly> {
    let s = x1.add(x2).square(f);
    let t1 = s.mul(&x3.square(f), f);
    let p12 = x1.mul(x2, f);
    let t2 = p12.mul(x3, f);
    let t3 = p12.square(f);
    let total = t1.add(&t2).add(&t3).add(&Sym::constant(b, f.n));
    total.coords
}

fn rank_u64(vectors: &[u64]) -> usize {
    let mut rows = vectors.to_vec();
    let mut rank = 0;
    for bit in (0..64).rev() {
        if let Some(p) = (rank..rows.len()).find(|&r| (rows[r] >> bit) & 1 == 1) {
            rows.swap(rank, p);
            let pv = rows[rank];
            for r in 0..rows.len() {
                if r != rank && (rows[r] >> bit) & 1 == 1 {
                    rows[r] ^= pv;
                }
            }
            rank += 1;
        }
    }
    rank
}

/// Normal basis (conjugates of a normal element) and P^{-1} rows.
fn normal_basis(n: u32, irr: &IrreduciblePoly, rng: &mut StdRng) -> (Vec<u64>, Vec<u64>) {
    loop {
        let alpha = elem_from_bits(rng.gen_range(1..(1u64 << n)), n);
        let conj: Vec<u64> = (0..n)
            .map(|j| bits_of(&alpha.square_k_times(j, irr)))
            .collect();
        if rank_u64(&conj) as u32 != n {
            continue;
        }
        let nn = n as usize;
        let mut rows: Vec<(u64, u64)> = (0..nn)
            .map(|i| {
                let mut p_row = 0u64;
                for (j, &c) in conj.iter().enumerate() {
                    if (c >> i) & 1 == 1 {
                        p_row |= 1 << j;
                    }
                }
                (p_row, 1u64 << i)
            })
            .collect();
        for col in 0..nn {
            let piv = (col..nn).find(|&r| (rows[r].0 >> col) & 1 == 1).unwrap();
            rows.swap(col, piv);
            let pv = rows[col];
            for r in 0..nn {
                if r != col && (rows[r].0 >> col) & 1 == 1 {
                    rows[r].0 ^= pv.0;
                    rows[r].1 ^= pv.1;
                }
            }
        }
        return (conj, rows.iter().map(|r| r.1).collect());
    }
}
fn to_normal(coords: &[Poly], pinv: &[u64]) -> Vec<Poly> {
    pinv.iter()
        .map(|row| {
            let mut acc: Poly = Vec::new();
            for (i, c) in coords.iter().enumerate() {
                if (row >> i) & 1 == 1 {
                    acc = poly_add(&acc, c);
                }
            }
            acc
        })
        .collect()
}

struct Layout {
    n: usize,
    s0: usize,
    s1: usize,
    s2: usize,
    c1: usize,
    c2: usize,
    u1: usize,
    u2: usize,
    n_vars: usize,
}
impl Layout {
    fn new(n: usize, l: usize) -> Self {
        let (s0, s1, s2) = (0, n, 2 * n);
        let (c1, c2) = (3 * n, 3 * n + l);
        let (u1, u2) = (3 * n + 2 * l, 4 * n + 2 * l);
        Layout {
            n,
            s0,
            s1,
            s2,
            c1,
            c2,
            u1,
            u2,
            n_vars: 5 * n + 2 * l,
        }
    }
    fn shift_mask(&self, m: &Mask) -> Mask {
        let mut out = *m;
        for base in [self.s0, self.s1, self.s2, self.u1, self.u2] {
            // rotate the n-bit block starting at `base` by one position upward
            let mut bits = Vec::with_capacity(self.n);
            for k in 0..self.n {
                let v = base + k;
                bits.push((out[v / 64] >> (v % 64)) & 1);
            }
            for k in 0..self.n {
                let v = base + ((k + 1) % self.n);
                if bits[k] == 1 {
                    out[v / 64] |= 1 << (v % 64);
                } else {
                    out[v / 64] &= !(1 << (v % 64));
                }
            }
        }
        out
    }
}

#[derive(Clone, Copy)]
enum Fam {
    S3(usize),
    U(u8, usize),
    OneLin(u8),
    OneQuad(u8, usize, usize),
}

fn build_symmetric(
    n: u32,
    l: u32,
    irr: &IrreduciblePoly,
    f: &Field,
    b: u64,
    r: &F2mElement,
    normal: &[u64],
    pinv: &[u64],
) -> (Layout, Vec<Poly>, Vec<Fam>) {
    let nn = n as usize;
    let lay = Layout::new(nn, l as usize);
    let window: Vec<F2mElement> = (0..l)
        .map(|k| F2mElement::from_bit_positions(&[k], n))
        .collect();
    let u1 = Sym::from_subspace_vars(normal, lay.u1, nn);
    let u2 = Sym::from_subspace_vars(normal, lay.u2, nn);
    let barrel = |s_off: usize, c_off: usize| -> Sym {
        let mut acc = Sym::zero(nn);
        for k in 0..nn {
            let shifted: Vec<u64> = window
                .iter()
                .map(|w| bits_of(&w.square_k_times(k as u32, irr)))
                .collect();
            acc = acc.add(&Sym::from_subspace_vars(&shifted, c_off, nn).scale_var(s_off + k));
        }
        acc
    };
    let mut xr = Sym::zero(nn);
    for k in 0..nn {
        xr = xr.add(
            &Sym::constant(bits_of(&r.square_k_times(k as u32, irr)), nn).scale_var(lay.s0 + k),
        );
    }
    let mut eqs = Vec::new();
    let mut fam = Vec::new();
    for (j, e) in to_normal(&sym_s3(&u1, &u2, &xr, b, f), pinv)
        .into_iter()
        .enumerate()
    {
        eqs.push(e);
        fam.push(Fam::S3(j));
    }
    for (i, (u, s_off, c_off)) in [(&u1, lay.s1, lay.c1), (&u2, lay.s2, lay.c2)]
        .into_iter()
        .enumerate()
    {
        let diff = u.add(&barrel(s_off, c_off));
        for (j, e) in to_normal(&diff.coords, pinv).into_iter().enumerate() {
            eqs.push(e);
            fam.push(Fam::U(i as u8, j));
        }
    }
    for (blk, off) in [(0u8, lay.s0), (1u8, lay.s1), (2u8, lay.s2)] {
        let mut monos: Vec<Mask> = (0..nn).map(|k| mask_var(off + k)).collect();
        monos.push([0u64; WORDS]);
        eqs.push(normalize(monos));
        fam.push(Fam::OneLin(blk));
        for k in 0..nn {
            for k2 in (k + 1)..nn {
                eqs.push(vec![mask_or(&mask_var(off + k), &mask_var(off + k2))]);
                fam.push(Fam::OneQuad(blk, k, k2));
            }
        }
    }
    (lay, eqs, fam)
}

fn fam_key(f: Fam) -> (u8, usize, usize, u8) {
    match f {
        Fam::S3(j) => (0, j, 0, 0),
        Fam::U(i, j) => (1, j, 0, i),
        Fam::OneLin(b) => (2, 0, 0, b),
        Fam::OneQuad(b, k, k2) => (3, k.min(k2), k.max(k2), b),
    }
}
fn shift_fam(f: Fam, n: usize) -> Fam {
    match f {
        Fam::S3(j) => Fam::S3((j + 1) % n),
        Fam::U(i, j) => Fam::U(i, (j + 1) % n),
        Fam::OneLin(b) => Fam::OneLin(b),
        Fam::OneQuad(b, k, k2) => Fam::OneQuad(b, (k + 1) % n, (k2 + 1) % n),
    }
}

/// Multipliers of degree ≤ d over nv variables.
fn multipliers(nv: usize, d: u32, out: &mut Vec<Mask>) {
    fn rec(start: usize, nv: usize, left: u32, cur: Mask, out: &mut Vec<Mask>) {
        out.push(cur);
        if left == 0 {
            return;
        }
        for v in start..nv {
            rec(v + 1, nv, left - 1, mask_or(&cur, &mask_var(v)), out);
        }
    }
    rec(0, nv, d, [0u64; WORDS], out);
}

// ── F_{2^d} ───────────────────────────────────────────────────────────
struct Gf {
    d: u32,
    poly: u64,
    exp: Vec<u32>,
    log: Vec<u32>,
}
impl Gf {
    fn new(d: u32) -> Self {
        for low in (1u64..(1u64 << d)).step_by(2) {
            let f = (1u64 << d) | low;
            if Self::is_irreducible(f, d) {
                let mut g = Gf {
                    d,
                    poly: f,
                    exp: Vec::new(),
                    log: Vec::new(),
                };
                if d <= 16 {
                    g.build_tables();
                }
                return g;
            }
        }
        panic!("no irreducible of degree {d}");
    }
    fn is_irreducible(f: u64, d: u32) -> bool {
        let g = Gf {
            d,
            poly: f,
            exp: Vec::new(),
            log: Vec::new(),
        };
        let mut x = 2u64;
        for _ in 1..=(d / 2) {
            x = g.mul_slow(x, x);
            if Self::poly_gcd(x ^ 2, f) != 1 {
                return false;
            }
        }
        true
    }
    fn poly_gcd(mut a: u64, mut b: u64) -> u64 {
        while b != 0 {
            if a == 0 {
                return b;
            }
            let (da, db) = (63 - a.leading_zeros(), 63 - b.leading_zeros());
            if da < db {
                std::mem::swap(&mut a, &mut b);
                continue;
            }
            a ^= b << (da - db);
        }
        a
    }
    fn mul_slow(&self, mut a: u64, mut b: u64) -> u64 {
        let mut r = 0u64;
        while b != 0 {
            if b & 1 == 1 {
                r ^= a;
            }
            b >>= 1;
            a <<= 1;
            if (a >> self.d) & 1 == 1 {
                a ^= self.poly;
            }
        }
        r
    }
    #[inline]
    fn mul(&self, a: u64, b: u64) -> u64 {
        if !self.exp.is_empty() {
            if a == 0 || b == 0 {
                return 0;
            }
            return self.exp[(self.log[a as usize] + self.log[b as usize]) as usize] as u64;
        }
        self.mul_slow(a, b)
    }
    fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1u64;
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
        self.pow(a, (1u64 << self.d) - 2)
    }
    fn build_tables(&mut self) {
        let order = (1u64 << self.d) - 1;
        let mut fs = Vec::new();
        let (mut m, mut q) = (order, 2u64);
        while q * q <= m {
            if m % q == 0 {
                fs.push(q);
                while m % q == 0 {
                    m /= q;
                }
            }
            q += 1;
        }
        if m > 1 {
            fs.push(m);
        }
        let mut g = 2u64;
        while !fs.iter().all(|&f| self.pow_slow(g, order / f) != 1) {
            g += 1;
        }
        let size = 1usize << self.d;
        let mut exp = vec![0u32; 2 * size];
        let mut log = vec![0u32; size];
        let mut x = 1u64;
        for i in 0..(order as usize) {
            exp[i] = x as u32;
            log[x as usize] = i as u32;
            x = self.mul_slow(x, g);
        }
        for i in (order as usize)..(2 * size) {
            exp[i] = exp[i - order as usize];
        }
        self.exp = exp;
        self.log = log;
    }
    fn pow_slow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1u64;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul_slow(r, a);
            }
            a = self.mul_slow(a, a);
            e >>= 1;
        }
        r
    }
    fn root_of_unity(&self, n: u64, rng: &mut StdRng) -> u64 {
        let order = (1u64 << self.d) - 1;
        assert_eq!(order % n, 0);
        loop {
            let g = rng.gen_range(2..(1u64 << self.d));
            let z = self.pow(g, order / n);
            if z != 1 {
                return z;
            }
        }
    }
}
fn character_orbits(n: usize) -> Vec<(usize, usize)> {
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for j in 1..n {
        if seen[j] {
            continue;
        }
        let (mut k, mut size) = (j, 0);
        while !seen[k] {
            seen[k] = true;
            size += 1;
            k = (2 * k) % n;
        }
        out.push((j, size));
    }
    out
}

// ── new for the full pass ─────────────────────────────────────────────

/// Reduced row echelon form over `F_{2^d}`: the nonzero reduced rows and
/// their pivot columns, in pivot order.
#[allow(dead_code)]
fn rref_gf(gf: &Gf, mut m: Vec<Vec<u64>>, cols: usize) -> (Vec<Vec<u64>>, Vec<usize>) {
    let mut rank = 0;
    let mut piv = Vec::new();
    for c in 0..cols {
        if rank == m.len() {
            break;
        }
        let Some(p) = (rank..m.len()).find(|&r| m[r][c] != 0) else {
            continue;
        };
        m.swap(rank, p);
        let inv = gf.inv(m[rank][c]);
        for v in m[rank].iter_mut() {
            if *v != 0 {
                *v = gf.mul(*v, inv);
            }
        }
        let pivot = m[rank].clone();
        let nz: Vec<usize> = (c..cols).filter(|&k| pivot[k] != 0).collect();
        for r in 0..m.len() {
            if r != rank && m[r][c] != 0 {
                let f = m[r][c];
                for &k in &nz {
                    m[r][k] ^= gf.mul(f, pivot[k]);
                }
            }
        }
        piv.push(c);
        rank += 1;
    }
    m.truncate(rank);
    (m, piv)
}

/// Reduced row echelon form over `F₂` (bit-packed rows).
#[allow(dead_code)]
fn rref_bits(mut m: Vec<Vec<u64>>, cols: usize) -> (Vec<Vec<u64>>, Vec<usize>) {
    let mut rank = 0;
    let mut piv = Vec::new();
    for c in 0..cols {
        if rank == m.len() {
            break;
        }
        let (w, bit) = (c / 64, 1u64 << (c % 64));
        let Some(p) = (rank..m.len()).find(|&r| m[r][w] & bit != 0) else {
            continue;
        };
        m.swap(rank, p);
        let pivot = m[rank].clone();
        for r in 0..m.len() {
            if r != rank && m[r][w] & bit != 0 {
                for (x, pv) in m[r][w..].iter_mut().zip(pivot[w..].iter()) {
                    *x ^= pv;
                }
            }
        }
        piv.push(c);
        rank += 1;
    }
    m.truncate(rank);
    (m, piv)
}

fn trace_table(gf: &Gf) -> Vec<u8> {
    (0..(1u64 << gf.d)).map(|x| trace(gf, x) as u8).collect()
}

fn trace(gf: &Gf, x: u64) -> u64 {
    let (mut acc, mut y) = (0u64, x);
    for _ in 0..gf.d {
        acc ^= y;
        y = gf.mul(y, y);
    }
    debug_assert!(acc <= 1);
    acc
}

// ════════════════════════════════════════════════════════════════════════
// New for item 02
// ════════════════════════════════════════════════════════════════════════

/// Thread CPU seconds (`CLOCK_THREAD_CPUTIME_ID`).
fn tcpu() -> f64 {
    let mut ts = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    unsafe {
        libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut ts);
    }
    ts.tv_sec as f64 + ts.tv_nsec as f64 * 1e-9
}

fn col_key(m: &Mask) -> (std::cmp::Reverse<u32>, Mask) {
    (std::cmp::Reverse(mask_deg(m)), *m)
}

/// Gauge `s₀ = e₀` on one monomial: `None` if it contains `s₀[k]`, `k ≠ 0`;
/// otherwise `s₀[0]` is dropped and every index is shifted down by `n`.
fn spec_mask(m: &Mask, n: usize) -> Option<Mask> {
    let low = m[0] & ((1u64 << n) - 1);
    if low & !1 != 0 {
        return None;
    }
    let mut out = [0u64; WORDS];
    for i in 0..WORDS {
        out[i] = m[i] >> n;
        if i + 1 < WORDS {
            out[i] |= m[i + 1] << (64 - n);
        }
    }
    Some(out)
}

fn spec_poly(p: &Poly, n: usize) -> Poly {
    normalize(p.iter().filter_map(|m| spec_mask(m, n)).collect())
}

/// Degree-`degree` Macaulay rows `t·f` (as monomial lists), as `build_macaulay`
/// enumerates them: every multiplier of degree `≤ degree − deg f`.
fn macaulay_rows(eqs: &[Poly], nv: usize, degree: u32) -> Vec<Poly> {
    let mut rows = Vec::new();
    for p in eqs {
        let pdeg = p.iter().map(mask_deg).max().unwrap_or(0);
        if p.is_empty() || pdeg > degree {
            continue;
        }
        let mut mults = Vec::new();
        multipliers(nv, degree - pdeg, &mut mults);
        for m in mults {
            let r = poly_mul_mono(p, &m);
            if !r.is_empty() {
                rows.push(r);
            }
        }
    }
    rows
}

/// A reduced row space in a degree-descending column order.
struct Reduced {
    cols: Vec<Mask>,
    rows: Vec<Vec<u64>>,
    piv: Vec<usize>,
    word_ops: u64,
}
impl Reduced {
    fn rank(&self) -> usize {
        self.piv.len()
    }
    fn refutes(&self) -> bool {
        self.piv.iter().any(|&p| mask_deg(&self.cols[p]) == 0)
    }
    fn count_below(&self, deg: u32) -> usize {
        self.piv
            .iter()
            .filter(|&&p| mask_deg(&self.cols[p]) < deg)
            .count()
    }
    fn row_poly(&self, i: usize) -> Poly {
        let mut out = Vec::new();
        for (w, &word) in self.rows[i].iter().enumerate() {
            let mut x = word;
            while x != 0 {
                let b = x.trailing_zeros() as usize;
                out.push(self.cols[w * 64 + b]);
                x &= x - 1;
            }
        }
        out.sort_unstable();
        out
    }
    /// The reduced rows of degree `< deg` (pivot degree `< deg`), as polynomials.
    fn polys_below(&self, deg: u32) -> Vec<Poly> {
        (0..self.piv.len())
            .filter(|&i| mask_deg(&self.cols[self.piv[i]]) < deg)
            .map(|i| self.row_poly(i))
            .collect()
    }
}

/// RREF with the crate kernel (`gf2_elim::rref_counted`), columns sorted by
/// descending degree, then mask.
fn reduce_rows(rows: &[Poly]) -> Reduced {
    let mut cols: Vec<Mask> = rows.iter().flatten().copied().collect();
    cols.sort_unstable_by_key(col_key);
    cols.dedup();
    let idx: HashMap<Mask, usize> = cols.iter().enumerate().map(|(i, c)| (*c, i)).collect();
    let words = cols.len().div_ceil(64).max(1);
    let mut mat: Vec<Vec<u64>> = rows
        .iter()
        .map(|r| {
            let mut row = vec![0u64; words];
            for c in r {
                let k = idx[c];
                row[k / 64] ^= 1u64 << (k % 64);
            }
            row
        })
        .collect();
    let mut ops = 0u64;
    let rank = gf2_elim::rref_counted(&mut mat, cols.len(), &mut ops);
    mat.truncate(rank);
    let piv: Vec<usize> = mat
        .iter()
        .map(|r| {
            let w = r
                .iter()
                .position(|&x| x != 0)
                .expect("pivot row is nonzero");
            w * 64 + r[w].trailing_zeros() as usize
        })
        .collect();
    assert!(piv.windows(2).all(|w| w[0] < w[1]), "pivot convention");
    Reduced {
        cols,
        rows: mat,
        piv,
        word_ops: ops,
    }
}

fn rref_gf_counted(
    gf: &Gf,
    mut m: Vec<Vec<u64>>,
    cols: usize,
    muls: &mut u64,
) -> (Vec<Vec<u64>>, Vec<usize>) {
    let mut rank = 0;
    let mut piv = Vec::new();
    for c in 0..cols {
        if rank == m.len() {
            break;
        }
        let Some(p) = (rank..m.len()).find(|&r| m[r][c] != 0) else {
            continue;
        };
        m.swap(rank, p);
        let inv = gf.inv(m[rank][c]);
        for v in m[rank].iter_mut() {
            if *v != 0 {
                *v = gf.mul(*v, inv);
                *muls += 1;
            }
        }
        let pivot = m[rank].clone();
        let nz: Vec<usize> = (c..cols).filter(|&k| pivot[k] != 0).collect();
        for r in 0..m.len() {
            if r != rank && m[r][c] != 0 {
                let f = m[r][c];
                for &k in &nz {
                    m[r][k] ^= gf.mul(f, pivot[k]);
                }
                *muls += nz.len() as u64;
            }
        }
        piv.push(c);
        rank += 1;
    }
    m.truncate(rank);
    (m, piv)
}

/// Per-size context shared by every target.
struct Ctx {
    n: u32,
    nn: usize,
    l: u32,
    d: u32,
    irr: IrreduciblePoly,
    f: Field,
    normal: Vec<u64>,
    pinv: Vec<u64>,
    gf: Gf,
    zeta: u64,
    ttab: Vec<u8>,
}
impl Ctx {
    fn new(n: u32, l: u32) -> Self {
        let nn = n as usize;
        let d = {
            let mut k = 1u32;
            let mut v = 2u64 % n as u64;
            while v != 1 {
                v = v * 2 % n as u64;
                k += 1;
            }
            k
        };
        for (_, size) in character_orbits(nn) {
            assert_eq!(
                size as u32, d,
                "n = {n}: a character orbit is not of size d"
            );
        }
        let irr = find_irreducible(n).expect("irreducible");
        let f = Field::new(n, &irr);
        let mut rng = StdRng::seed_from_u64(0x62 + n as u64);
        let (normal, pinv) = normal_basis(n, &irr, &mut rng);
        let gf = Gf::new(d);
        let zeta = gf.root_of_unity(n as u64, &mut rng);
        let ttab = trace_table(&gf);
        Ctx {
            n,
            nn,
            l,
            d,
            irr,
            f,
            normal,
            pinv,
            gf,
            zeta,
            ttab,
        }
    }
    fn target(&self, t: u64) -> u64 {
        let mut rng = StdRng::seed_from_u64(20261014 ^ (self.n as u64 * 1000 + t));
        rng.gen_range(1..(1u64 << self.n))
    }
    fn fmul(&self, a: u64, b: u64) -> u64 {
        let mut acc = 0u64;
        for i in 0..self.nn {
            if (a >> i) & 1 == 0 {
                continue;
            }
            for j in 0..self.nn {
                if (b >> j) & 1 == 1 {
                    acc ^= self.f.reduced[i][j];
                }
            }
        }
        acc
    }
    fn normal_coords(&self, x: u64) -> u64 {
        let mut u = 0u64;
        for (j, row) in self.pinv.iter().enumerate() {
            if (row & x).count_ones() & 1 == 1 {
                u |= 1 << j;
            }
        }
        u
    }
    /// Every root of the gauge-fixed system, by enumerating
    /// `(k₁, c₁, k₂, c₂)` and testing `S₃(x₁, x₂, x_R) = 0` directly.
    fn brute_roots(&self, r: u64) -> Vec<u64> {
        let (nn, l) = (self.nn, self.l as usize);
        let xs: Vec<Vec<u64>> = (0..nn)
            .map(|k| {
                (0..(1u64 << l))
                    .map(|c| {
                        bits_of(&elem_from_bits(c, self.n).square_k_times(k as u32, &self.irr))
                    })
                    .collect()
            })
            .collect();
        let r2 = self.fmul(r, r);
        let mut roots = Vec::new();
        for k1 in 0..nn {
            for c1 in 0..(1u64 << l) {
                let x1 = xs[k1][c1 as usize];
                for k2 in 0..nn {
                    for c2 in 0..(1u64 << l) {
                        let x2 = xs[k2][c2 as usize];
                        let s = x1 ^ x2;
                        let p12 = self.fmul(x1, x2);
                        let v = self.fmul(self.fmul(s, s), r2)
                            ^ self.fmul(p12, r)
                            ^ self.fmul(p12, p12)
                            ^ 1;
                        if v == 0 {
                            let mut mask = 0u64;
                            mask |= 1 << k1;
                            mask |= 1 << (nn + k2);
                            mask |= c1 << (2 * nn);
                            mask |= c2 << (2 * nn + l);
                            mask |= self.normal_coords(x1) << (2 * nn + 2 * l);
                            mask |= self.normal_coords(x2) << (3 * nn + 2 * l);
                            roots.push(mask);
                        }
                    }
                }
            }
        }
        roots.sort_unstable();
        roots
    }
}

/// Output of the candidate's equivariant pass.
struct PassOut {
    /// Inverse-transformed degree-fall rows, gauge-specialised.
    spec_rows: Vec<Poly>,
    sym_rank: usize,
    sym_fall: usize,
    triv_word_ops: u64,
    gf_muls: u64,
    row_orbits: usize,
    sym_rows: usize,
}

/// The G2 equivariant degree-`degree` pass on the symmetric system, with
/// the inverse transform of the degree-fall rows and the gauge `s₀ = e₀`.
fn equivariant_pass(cx: &Ctx, lay: &Layout, eqs: &[Poly], fam: &[Fam], degree: u32) -> PassOut {
    let nn = cx.nn;
    let index: HashMap<(u8, usize, usize, u8), usize> = fam
        .iter()
        .enumerate()
        .map(|(i, &f)| (fam_key(f), i))
        .collect();
    // Row orbit representatives.
    let mut reps: Vec<(Mask, usize)> = Vec::new();
    let mut total_rows = 0usize;
    for (ei, p) in eqs.iter().enumerate() {
        let pdeg = p.iter().map(mask_deg).max().unwrap_or(0);
        if pdeg > degree {
            continue;
        }
        let mut mults = Vec::new();
        multipliers(lay.n_vars, degree - pdeg, &mut mults);
        for m in mults {
            if poly_mul_mono(p, &m).is_empty() {
                continue;
            }
            total_rows += 1;
            let mut canon = (m, ei);
            let (mut mm, mut ee) = (m, fam[ei]);
            for _ in 1..nn {
                mm = lay.shift_mask(&mm);
                ee = shift_fam(ee, nn);
                canon = canon.min((mm, index[&fam_key(ee)]));
            }
            if canon == (m, ei) {
                reps.push((m, ei));
            }
        }
    }
    // Column orbits.
    let mut col_info: HashMap<Mask, (Mask, u32, bool)> = HashMap::new();
    let mut canon_set: Vec<Mask> = Vec::new();
    let mut rep_cols: Vec<Vec<Mask>> = Vec::with_capacity(reps.len());
    for &(m, ei) in &reps {
        let cs = poly_mul_mono(&eqs[ei], &m);
        for c in &cs {
            if col_info.contains_key(c) {
                continue;
            }
            let mut canon = *c;
            let mut mm = *c;
            for _ in 1..nn {
                mm = lay.shift_mask(&mm);
                canon = canon.min(mm);
            }
            let fixed = lay.shift_mask(&canon) == canon;
            let mut mm = canon;
            for t in 0..(if fixed { 1 } else { nn }) {
                col_info.insert(mm, (canon, t as u32, fixed));
                mm = lay.shift_mask(&mm);
            }
            canon_set.push(canon);
        }
        rep_cols.push(cs);
    }
    canon_set.sort_unstable_by_key(col_key);
    let triv_index: HashMap<Mask, usize> =
        canon_set.iter().enumerate().map(|(i, c)| (*c, i)).collect();
    let free_canon: Vec<Mask> = canon_set
        .iter()
        .copied()
        .filter(|c| !col_info[c].2)
        .collect();
    let free_index: HashMap<Mask, usize> = free_canon
        .iter()
        .enumerate()
        .map(|(i, c)| (*c, i))
        .collect();
    let n_triv = canon_set.len();
    let orbit = |canon: &Mask, fixed: bool| -> Vec<Mask> {
        let mut out = Vec::with_capacity(nn);
        let mut mm = *canon;
        for _ in 0..(if fixed { 1 } else { nn }) {
            out.push(mm);
            mm = lay.shift_mask(&mm);
        }
        out
    };

    // Trivial block over F₂ with the crate kernel.
    let words = n_triv.div_ceil(64).max(1);
    let mut t_rows: Vec<Vec<u64>> = Vec::with_capacity(reps.len());
    for cs in &rep_cols {
        let mut row = vec![0u64; words];
        for c in cs {
            let k = triv_index[&col_info[c].0];
            row[k / 64] ^= 1u64 << (k % 64);
        }
        t_rows.push(row);
    }
    let mut triv_ops = 0u64;
    let triv_rank = gf2_elim::rref_counted(&mut t_rows, n_triv, &mut triv_ops);
    t_rows.truncate(triv_rank);
    let mut sym_rank = triv_rank;
    let mut sym_fall = 0usize;
    let mut fall_sym: Vec<Poly> = Vec::new();
    for row in &t_rows {
        let w = row.iter().position(|&x| x != 0).unwrap();
        let p = w * 64 + row[w].trailing_zeros() as usize;
        if mask_deg(&canon_set[p]) >= degree {
            continue;
        }
        sym_fall += 1;
        let mut u: Poly = Vec::new();
        for (wi, &word) in row.iter().enumerate() {
            let mut x = word;
            while x != 0 {
                let k = wi * 64 + x.trailing_zeros() as usize;
                u.extend(orbit(&canon_set[k], col_info[&canon_set[k]].2));
                x &= x - 1;
            }
        }
        u.sort_unstable();
        fall_sym.push(u);
    }

    // Character blocks over F_{2^d}.
    let gf = &cx.gf;
    let zeta_inv = gf.inv(cx.zeta);
    let zeta_pows: Vec<u64> = (0..nn).map(|t| gf.pow(cx.zeta, t as u64)).collect();
    let mut gf_muls = 0u64;
    for (j, size) in character_orbits(nn) {
        let zj = gf.pow(zeta_inv, j as u64);
        let zpow: Vec<u64> = (0..nn).map(|t| gf.pow(zj, t as u64)).collect();
        let mut rows_k: Vec<Vec<u64>> = Vec::new();
        for (ri, &(m, ei)) in reps.iter().enumerate() {
            if lay.shift_mask(&m) == m && index[&fam_key(shift_fam(fam[ei], nn))] == ei {
                continue;
            }
            let mut row = vec![0u64; free_canon.len()];
            let mut any = false;
            for c in &rep_cols[ri] {
                let (canon, t, fixed) = col_info[c];
                if fixed {
                    continue;
                }
                row[free_index[&canon]] ^= zpow[t as usize];
                any = true;
            }
            if any {
                rows_k.push(row);
            }
        }
        let (rr, pv) = rref_gf_counted(gf, rows_k, free_canon.len(), &mut gf_muls);
        sym_rank += rr.len() * size;
        let zjt: Vec<u64> = (0..nn).map(|t| zeta_pows[(j * t) % nn]).collect();
        for (row, &p) in rr.iter().zip(&pv) {
            if mask_deg(&free_canon[p]) >= degree {
                continue;
            }
            sym_fall += size;
            for b in 0..cx.d {
                let lambda = 1u64 << b;
                let mut u: Poly = Vec::new();
                for (k, &wk) in row.iter().enumerate() {
                    if wk == 0 {
                        continue;
                    }
                    let lw = gf.mul(lambda, wk);
                    let members = orbit(&free_canon[k], false);
                    for (t, mm) in members.iter().enumerate() {
                        if cx.ttab[gf.mul(lw, zjt[t]) as usize] == 1 {
                            u.push(*mm);
                        }
                    }
                }
                u.sort_unstable();
                fall_sym.push(u);
            }
        }
    }
    let spec_rows: Vec<Poly> = fall_sym
        .iter()
        .map(|p| spec_poly(p, nn))
        .filter(|p| !p.is_empty())
        .collect();
    PassOut {
        spec_rows,
        sym_rank,
        sym_fall,
        triv_word_ops: triv_ops,
        gf_muls,
        row_orbits: reps.len(),
        sym_rows: total_rows,
    }
}

fn gauge_system(eqs: &[Poly], n: usize) -> Vec<Poly> {
    eqs.iter()
        .map(|p| spec_poly(p, n))
        .filter(|p| !p.is_empty())
        .collect()
}

fn to_u64(p: &Poly, nv: usize) -> F2BoolPoly {
    let monos = p
        .iter()
        .map(|m| {
            assert!(m[1] == 0 && m[2] == 0, "gauge system exceeds 64 unknowns");
            F2BoolMono::from_mask(m[0])
        })
        .collect();
    F2BoolPoly::from_monos(monos, nv)
}

fn sorted_polys(mut v: Vec<Poly>) -> Vec<Poly> {
    v.sort();
    v
}

/// Incumbent step: gauge-fixed degree-`degree` Macaulay matrix, crate kernel.
fn incumbent_step(geqs: &[Poly], gnv: usize, degree: u32) -> (Reduced, usize, f64) {
    let c0 = tcpu();
    let rows = macaulay_rows(geqs, gnv, degree);
    let nrows = rows.len();
    let red = reduce_rows(&rows);
    (red, nrows, tcpu() - c0)
}

fn run_step(cx: &Ctx, t: u64, degree: u32, check: bool, posctl: bool) {
    let r = cx.target(t);
    let rel = elem_from_bits(r, cx.n);
    let (lay, eqs, fam) =
        build_symmetric(cx.n, cx.l, &cx.irr, &cx.f, 1, &rel, &cx.normal, &cx.pinv);
    let geqs = gauge_system(&eqs, cx.nn);
    let gnv = lay.n_vars - cx.nn;
    let wall = Instant::now();
    // Arms interleaved: I, C, I (A/A).
    let (ired, irows, i_cpu) = incumbent_step(&geqs, gnv, degree);
    let c0 = tcpu();
    let pass = equivariant_pass(cx, &lay, &eqs, &fam, degree);
    let pass_cpu = tcpu() - c0;
    let c1 = tcpu();
    let cred = reduce_rows(&pass.spec_rows);
    let c_cpu = pass_cpu + (tcpu() - c1);
    let (ired2, _, i2_cpu) = incumbent_step(&geqs, gnv, degree);
    assert_eq!(ired2.piv, ired.piv);
    drop(ired2);
    // Gate: refutation bit and the linear-consequence space.
    let i_lin = sorted_polys(ired.polys_below(2));
    let c_lin = sorted_polys(cred.polys_below(2));
    let gate = ired.refutes() == cred.refutes() && i_lin == c_lin;
    // Containment span(spec rows) ⊆ span(I rows) (not charged).
    let contained = if check {
        let mut rows: Vec<Poly> = (0..ired.rank()).map(|i| ired.row_poly(i)).collect();
        rows.extend((0..cred.rank()).map(|i| cred.row_poly(i)));
        Some(reduce_rows(&rows).rank() == ired.rank())
    } else {
        None
    };
    // Positive control: crate kernel on the full symmetric matrix.
    let pos = if posctl {
        let c0 = tcpu();
        let rows = macaulay_rows(&eqs, lay.n_vars, degree);
        let red = reduce_rows(&rows);
        let cpu = tcpu() - c0;
        Some((cpu, red.rank(), red.word_ops))
    } else {
        None
    };
    let opt = |x: Option<String>| x.unwrap_or_else(|| "null".into());
    println!(
        r#"{{"mode":"step","n":{},"l":{},"d":{},"degree":{degree},"target":{t},"x_r":{r},"gauge_vars":{gnv},"gauge_eqs":{},"gauge_rows":{irows},"gauge_cols":{},"gauge_rank":{},"gauge_fall":{},"gauge_linear":{},"gauge_refutes":{},"i_word_ops":{},"i_cpu":{i_cpu:.4},"i2_cpu":{i2_cpu:.4},"sym_vars":{},"sym_rows":{},"row_orbits":{},"sym_rank":{},"sym_fall":{},"spec_rows":{},"spec_rank":{},"spec_linear":{},"spec_refutes":{},"c_triv_word_ops":{},"c_gf_muls":{},"c_spec_word_ops":{},"c_pass_cpu":{pass_cpu:.4},"c_cpu":{c_cpu:.4},"gate":{gate},"contained":{},"pos_cpu":{},"pos_rank":{},"pos_word_ops":{},"wall_s":{:.2}}}"#,
        cx.n,
        cx.l,
        cx.d,
        geqs.len(),
        ired.cols.len(),
        ired.rank(),
        ired.count_below(degree),
        i_lin.len(),
        ired.refutes(),
        ired.word_ops,
        lay.n_vars,
        pass.sym_rows,
        pass.row_orbits,
        pass.sym_rank,
        pass.sym_fall,
        pass.spec_rows.len(),
        cred.rank(),
        c_lin.len(),
        cred.refutes(),
        pass.triv_word_ops,
        pass.gf_muls,
        cred.word_ops,
        opt(contained.map(|b| b.to_string())),
        opt(pos.map(|p| format!("{:.4}", p.0))),
        opt(pos.map(|p| p.1.to_string())),
        opt(pos.map(|p| p.2.to_string())),
        wall.elapsed().as_secs_f64(),
    );
}

fn solve(sys: &[F2BoolPoly], nv: usize) -> (Vec<u64>, f64, u64, usize, bool) {
    let opts = SolveOptions {
        engine: SolverEngine::default(),
        max_solutions: usize::MAX,
        node_budget: usize::MAX,
        split_rule: split_rule_default(),
    };
    let w0 = f4_word_ops_thread();
    let c0 = tcpu();
    let (mut roots, stats) = solve_boolean_system(sys, nv, &opts);
    let cpu = tcpu() - c0;
    roots.sort_unstable();
    (
        roots,
        cpu,
        f4_word_ops_thread() - w0,
        stats.splits,
        stats.exhausted,
    )
}

fn run_oracle(cx: &Ctx, t: u64, degree: u32) {
    let r = cx.target(t);
    let rel = elem_from_bits(r, cx.n);
    let (lay, eqs, fam) =
        build_symmetric(cx.n, cx.l, &cx.irr, &cx.f, 1, &rel, &cx.normal, &cx.pinv);
    let gnv = lay.n_vars - cx.nn;
    assert!(
        gnv <= 64,
        "gauge system needs ≤ 64 unknowns for the pipeline solver"
    );
    let truth = cx.brute_roots(r);
    let wall = Instant::now();
    // I: the incumbent oracle (system construction charged to both arms).
    let c0 = tcpu();
    let geqs = gauge_system(&eqs, cx.nn);
    let gsys: Vec<F2BoolPoly> = geqs.iter().map(|p| to_u64(p, gnv)).collect();
    let build_cpu = tcpu() - c0;
    let (i_roots, i_solve, i_ops, i_splits, i_exh) = solve(&gsys, gnv);
    let i_cpu = build_cpu + i_solve;
    // C: equivariant root pass, specialised fall rows appended, same solver.
    let c0 = tcpu();
    let pass = equivariant_pass(cx, &lay, &eqs, &fam, degree);
    let cred = reduce_rows(&pass.spec_rows);
    let mut csys = gsys.clone();
    csys.extend((0..cred.rank()).map(|i| to_u64(&cred.row_poly(i), gnv)));
    let pass_cpu = tcpu() - c0;
    let (c_roots, c_solve, c_ops, c_splits, c_exh) = solve(&csys, gnv);
    let c_cpu = build_cpu + pass_cpu + c_solve;
    // F (attribution control): the incumbent's own degree-3 Macaulay step on
    // the gauge-fixed system, its degree-fall rows appended, same solver.
    let c0 = tcpu();
    let (fred, _, _) = incumbent_step(&geqs, gnv, degree);
    let mut fsys = gsys.clone();
    fsys.extend(fred.polys_below(degree).iter().map(|p| to_u64(p, gnv)));
    let f_extra = fsys.len() - gsys.len();
    let f_step = tcpu() - c0;
    let (f_roots, f_solve, f_ops, f_splits, f_exh) = solve(&fsys, gnv);
    let f_cpu = build_cpu + f_step + f_solve;
    // I again (A/A).
    let (i2_roots, i2_solve, _, _, _) = solve(&gsys, gnv);
    let i2_cpu = build_cpu + i2_solve;
    let gate = !i_exh && !c_exh && i_roots == truth && c_roots == truth && i2_roots == truth;
    let f_gate = !f_exh && f_roots == truth;
    println!(
        r#"{{"mode":"oracle","n":{},"l":{},"d":{},"degree":{degree},"target":{t},"x_r":{r},"gauge_vars":{gnv},"roots":{},"decomposes":{},"i_cpu":{i_cpu:.5},"i2_cpu":{i2_cpu:.5},"i_word_ops":{i_ops},"i_splits":{i_splits},"c_cpu":{c_cpu:.5},"c_pass_cpu":{pass_cpu:.5},"c_solve_cpu":{c_solve:.5},"c_word_ops":{c_ops},"c_triv_word_ops":{},"c_gf_muls":{},"c_spec_word_ops":{},"c_extra_rows":{},"c_splits":{c_splits},"f_cpu":{f_cpu:.5},"f_step_cpu":{f_step:.5},"f_solve_cpu":{f_solve:.5},"f_word_ops":{f_ops},"f_extra_rows":{f_extra},"f_splits":{f_splits},"f_gate":{f_gate},"gate":{gate},"wall_s":{:.2}}}"#,
        cx.n,
        cx.l,
        cx.d,
        truth.len(),
        !truth.is_empty(),
        pass.triv_word_ops,
        pass.gf_muls,
        cred.word_ops,
        cred.rank(),
        wall.elapsed().as_secs_f64(),
    );
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let mut ns = vec![7u32];
    let (mut l, mut degree, mut targets, mut first) = (3u32, 3u32, 4u64, 0u64);
    let mut mode = "step".to_string();
    let (mut check, mut posctl) = (false, false);
    let mut i = 1;
    while i < argv.len() {
        match argv[i].as_str() {
            "--n" => {
                i += 1;
                ns = argv[i].split(',').filter_map(|s| s.parse().ok()).collect();
            }
            "--l" => {
                i += 1;
                l = argv[i].parse().unwrap();
            }
            "--degree" => {
                i += 1;
                degree = argv[i].parse().unwrap();
            }
            "--targets" => {
                i += 1;
                targets = argv[i].parse().unwrap();
            }
            "--first" => {
                i += 1;
                first = argv[i].parse().unwrap();
            }
            "--mode" => {
                i += 1;
                mode = argv[i].clone();
            }
            "--check" => check = true,
            "--posctl" => posctl = true,
            _ => {}
        }
        i += 1;
    }
    for &n in &ns {
        let cx = Ctx::new(n, l);
        for t in first..first + targets {
            match mode.as_str() {
                "step" => run_step(&cx, t, degree, check, posctl),
                "oracle" => run_oracle(&cx, t, degree),
                _ => panic!("--mode step|oracle"),
            }
        }
    }
}
