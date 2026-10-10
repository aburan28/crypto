//! G2 wide-mask probe: the equivariant block reduction of the barrel
//! decomposition system at `n = 13 … 31`, beyond the 64-variable engine.
//!
//! Same formulation as `g2_equivariant_orbit_probe` (normal-basis Weil
//! restriction, auxiliary `u_i = Σ_k s_i[k] τ^k(y_i)`, one-hot selectors for
//! both summands and the target), re-implemented over 192-bit monomial masks
//! so that `5n + 2l ≤ 192`. It builds the degree-`D` Macaulay matrix of the
//! symmetric system, splits it into the trivial block over `F₂` and one
//! character block per Galois orbit over `F_{2^d}`, computes the block ranks,
//! and optionally (`--control`) the full `F₂` rank by dense elimination for
//! the identity check and the wall-clock ratio.
//!
//! ```bash
//! cargo run --release --example g2_wide_probe -- --n 13,17 --l 3 --degree 3 [--control]
//! ```

use std::collections::HashMap;
use std::time::Instant;

use crypto_lib::binary_ecc::f2m::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

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
        Field { n: nn, reduced, squares }
    }
}

/// Symbolic field element: n coordinate polynomials in the poly basis.
#[derive(Clone)]
struct Sym {
    coords: Vec<Poly>,
}
impl Sym {
    fn zero(n: usize) -> Self {
        Sym { coords: vec![Vec::new(); n] }
    }
    fn constant(v: u64, n: usize) -> Self {
        Sym { coords: (0..n).map(|k| if (v >> k) & 1 == 1 { vec![[0u64; WORDS]] } else { Vec::new() }).collect() }
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
        Sym { coords: coords.into_iter().map(normalize).collect() }
    }
    fn add(&self, o: &Self) -> Self {
        Sym { coords: self.coords.iter().zip(o.coords.iter()).map(|(a, b)| poly_add(a, b)).collect() }
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
        Sym { coords: acc.into_iter().map(normalize).collect() }
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
        Sym { coords: acc.into_iter().map(normalize).collect() }
    }
    /// Multiply every coordinate by a Boolean variable.
    fn scale_var(&self, v: usize) -> Self {
        let m = mask_var(v);
        Sym { coords: self.coords.iter().map(|p| poly_mul_mono(p, &m)).collect() }
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
        let conj: Vec<u64> = (0..n).map(|j| bits_of(&alpha.square_k_times(j, irr))).collect();
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
        Layout { n, s0, s1, s2, c1, c2, u1, u2, n_vars: 5 * n + 2 * l }
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

fn build_symmetric(n: u32, l: u32, irr: &IrreduciblePoly, f: &Field, b: u64, r: &F2mElement, normal: &[u64], pinv: &[u64]) -> (Layout, Vec<Poly>, Vec<Fam>) {
    let nn = n as usize;
    let lay = Layout::new(nn, l as usize);
    let window: Vec<F2mElement> = (0..l).map(|k| F2mElement::from_bit_positions(&[k], n)).collect();
    let u1 = Sym::from_subspace_vars(normal, lay.u1, nn);
    let u2 = Sym::from_subspace_vars(normal, lay.u2, nn);
    let barrel = |s_off: usize, c_off: usize| -> Sym {
        let mut acc = Sym::zero(nn);
        for k in 0..nn {
            let shifted: Vec<u64> = window.iter().map(|w| bits_of(&w.square_k_times(k as u32, irr))).collect();
            acc = acc.add(&Sym::from_subspace_vars(&shifted, c_off, nn).scale_var(s_off + k));
        }
        acc
    };
    let mut xr = Sym::zero(nn);
    for k in 0..nn {
        xr = xr.add(&Sym::constant(bits_of(&r.square_k_times(k as u32, irr)), nn).scale_var(lay.s0 + k));
    }
    let mut eqs = Vec::new();
    let mut fam = Vec::new();
    for (j, e) in to_normal(&sym_s3(&u1, &u2, &xr, b, f), pinv).into_iter().enumerate() {
        eqs.push(e);
        fam.push(Fam::S3(j));
    }
    for (i, (u, s_off, c_off)) in [(&u1, lay.s1, lay.c1), (&u2, lay.s2, lay.c2)].into_iter().enumerate() {
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
                let mut g = Gf { d, poly: f, exp: Vec::new(), log: Vec::new() };
                if d <= 16 {
                    g.build_tables();
                }
                return g;
            }
        }
        panic!("no irreducible of degree {d}");
    }
    fn is_irreducible(f: u64, d: u32) -> bool {
        let g = Gf { d, poly: f, exp: Vec::new(), log: Vec::new() };
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
fn rank_gf(gf: &Gf, mut m: Vec<Vec<u64>>) -> usize {
    if m.is_empty() {
        return 0;
    }
    let cols = m[0].len();
    let mut rank = 0;
    for c in 0..cols {
        let Some(p) = (rank..m.len()).find(|&r| m[r][c] != 0) else { continue };
        m.swap(rank, p);
        let inv = gf.inv(m[rank][c]);
        for v in m[rank].iter_mut() {
            *v = gf.mul(*v, inv);
        }
        let pivot = m[rank].clone();
        for r in 0..m.len() {
            if r != rank && m[r][c] != 0 {
                let f = m[r][c];
                for (x, pv) in m[r].iter_mut().zip(pivot.iter()) {
                    if *pv != 0 {
                        *x ^= gf.mul(f, *pv);
                    }
                }
            }
        }
        rank += 1;
        if rank == m.len() {
            break;
        }
    }
    rank
}
fn rank_bits(mut m: Vec<Vec<u64>>, cols: usize) -> usize {
    let mut rank = 0;
    for c in 0..cols {
        let (w, bit) = (c / 64, 1u64 << (c % 64));
        let Some(p) = (rank..m.len()).find(|&r| m[r][w] & bit != 0) else { continue };
        m.swap(rank, p);
        let pivot = m[rank].clone();
        for r in 0..m.len() {
            if r != rank && m[r][w] & bit != 0 {
                for (x, pv) in m[r].iter_mut().zip(pivot.iter()) {
                    *x ^= pv;
                }
            }
        }
        rank += 1;
        if rank == m.len() {
            break;
        }
    }
    rank
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

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let mut ns = vec![13u32, 17];
    let mut l = 3u32;
    let mut degree = 3u32;
    let mut control = false;
    let mut i = 1;
    while i < argv.len() {
        match argv[i].as_str() {
            "--n" => {
                i += 1;
                ns = argv[i].split(',').filter_map(|s| s.parse().ok()).collect();
            }
            "--l" => {
                i += 1;
                l = argv[i].parse().unwrap_or(3);
            }
            "--degree" => {
                i += 1;
                degree = argv[i].parse().unwrap_or(3);
            }
            "--control" => control = true,
            _ => {}
        }
        i += 1;
    }
    println!("# G2 wide probe (192-variable masks)\n");
    println!("| n | l | d | vars | eqs | degree | sym rows | sym cols | block rows | block cols (free) | trivial rank | character block ranks | blocks s | full rank | full s | ratio | identity |");
    println!("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---|");
    for &n in &ns {
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
        let irr = find_irreducible(n).expect("irreducible");
        let f = Field::new(n, &irr);
        let mut rng = StdRng::seed_from_u64(0x62 + n as u64);
        let (normal, pinv) = normal_basis(n, &irr, &mut rng);
        let r = elem_from_bits(rng.gen_range(1..(1u64 << n)), n);
        let (lay, eqs, fam) = build_symmetric(n, l, &irr, &f, 1, &r, &normal, &pinv);
        if lay.n_vars > 64 * WORDS {
            println!("| {n} | {l} | {d} | {} | — | — | too many variables | | | | | | | | | | |", lay.n_vars);
            continue;
        }
        let index: HashMap<(u8, usize, usize, u8), usize> = fam.iter().enumerate().map(|(i, &f)| (fam_key(f), i)).collect();
        // σ-invariance check of the equation set
        let mut invariant = true;
        for (ei, p) in eqs.iter().enumerate() {
            let tgt = &eqs[index[&fam_key(shift_fam(fam[ei], nn))]];
            let mut sh: Vec<Mask> = p.iter().map(|m| lay.shift_mask(m)).collect();
            sh.sort_unstable();
            if &sh != tgt {
                invariant = false;
                break;
            }
        }
        if !invariant {
            println!("<!-- n={n}: equation set NOT σ-invariant — bug -->");
            continue;
        }
        let t_build = Instant::now();
        // Macaulay rows (representatives only) at the given degree.
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
                let cols = poly_mul_mono(p, &m);
                if cols.is_empty() {
                    continue;
                }
                total_rows += 1;
                let mut canon = (m, ei);
                let (mut mm, mut ee) = (m, fam[ei]);
                for _ in 1..nn {
                    mm = lay.shift_mask(&mm);
                    ee = shift_fam(ee, nn);
                    let eidx = index[&fam_key(ee)];
                    canon = canon.min((mm, eidx));
                }
                if canon == (m, ei) {
                    reps.push((m, ei));
                }
            }
        }
        // columns of representatives, orbit classification
        let mut col_info: HashMap<Mask, (Mask, u32, bool)> = HashMap::new();
        let mut col_index: HashMap<Mask, usize> = HashMap::new();
        let mut total_cols: std::collections::HashSet<Mask> = std::collections::HashSet::new();
        let mut rep_cols: Vec<Vec<Mask>> = Vec::with_capacity(reps.len());
        for &(m, ei) in &reps {
            let cs = poly_mul_mono(&eqs[ei], &m);
            for c in &cs {
                if !col_info.contains_key(c) {
                    let mut canon = *c;
                    let mut mm = *c;
                    for _ in 1..nn {
                        mm = lay.shift_mask(&mm);
                        canon = canon.min(mm);
                    }
                    let fixed = lay.shift_mask(&canon) == canon;
                    let mut t = 0u32;
                    let mut mm = canon;
                    while mm != *c {
                        mm = lay.shift_mask(&mm);
                        t += 1;
                    }
                    let next = col_index.len();
                    col_index.entry(canon).or_insert(next);
                    col_info.insert(*c, (canon, t, fixed));
                }
                // full column set (for reporting): all shifts of c
                let mut mm = *c;
                for _ in 0..nn {
                    total_cols.insert(mm);
                    mm = lay.shift_mask(&mm);
                }
            }
            rep_cols.push(cs);
        }
        let build_s = t_build.elapsed().as_secs_f64();
        let ncols = col_index.len();
        let t0 = Instant::now();
        // trivial block
        let words = ncols.div_ceil(64);
        let mut t_rows: Vec<Vec<u64>> = Vec::with_capacity(reps.len());
        for cs in &rep_cols {
            let mut row = vec![0u64; words];
            for c in cs {
                let (canon, _, _) = col_info[c];
                let k = col_index[&canon];
                row[k / 64] ^= 1u64 << (k % 64);
            }
            t_rows.push(row);
        }
        let rank_triv = rank_bits(t_rows, ncols);
        // character blocks
        let gf = Gf::new(d);
        let zeta_inv = gf.inv(gf.root_of_unity(n as u64, &mut rng));
        let free_cols: Vec<Mask> = {
            let mut v: Vec<Mask> = col_info.values().filter(|(_, _, fx)| !fx).map(|(c, _, _)| *c).collect();
            v.sort_unstable();
            v.dedup();
            v
        };
        let free_index: HashMap<Mask, usize> = free_cols.iter().enumerate().map(|(i, c)| (*c, i)).collect();
        let mut block_ranks = Vec::new();
        for (j, size) in character_orbits(nn) {
            let zj = gf.pow(zeta_inv, j as u64);
            let zpow: Vec<u64> = (0..nn).map(|t| gf.pow(zj, t as u64)).collect();
            let mut rows_k: Vec<Vec<u64>> = Vec::new();
            for (ri, &(m, ei)) in reps.iter().enumerate() {
                let fixed_row = lay.shift_mask(&m) == m && index[&fam_key(shift_fam(fam[ei], nn))] == ei;
                if fixed_row {
                    continue;
                }
                let mut row = vec![0u64; free_cols.len()];
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
            block_ranks.push((j, size, rank_gf(&gf, rows_k)));
        }
        let blocks_s = t0.elapsed().as_secs_f64();
        let sum: usize = rank_triv + block_ranks.iter().map(|(_, s, r)| s * r).sum::<usize>();
        // optional full control: dense F_2 rank of the whole symmetric Macaulay matrix
        let (full_rank, full_s) = if control {
            let t1 = Instant::now();
            let all_cols: Vec<Mask> = {
                let mut v: Vec<Mask> = total_cols.iter().copied().collect();
                v.sort_unstable();
                v
            };
            let cidx: HashMap<Mask, usize> = all_cols.iter().enumerate().map(|(i, c)| (*c, i)).collect();
            let w = all_cols.len().div_ceil(64);
            let mut rows: Vec<Vec<u64>> = Vec::with_capacity(total_rows);
            for (ei, p) in eqs.iter().enumerate() {
                let pdeg = p.iter().map(mask_deg).max().unwrap_or(0);
                if pdeg > degree {
                    continue;
                }
                let mut mults = Vec::new();
                multipliers(lay.n_vars, degree - pdeg, &mut mults);
                for m in mults {
                    let cs = poly_mul_mono(p, &m);
                    if cs.is_empty() {
                        continue;
                    }
                    let mut row = vec![0u64; w];
                    for c in &cs {
                        let k = cidx[c];
                        row[k / 64] ^= 1u64 << (k % 64);
                    }
                    rows.push(row);
                }
            }
            let rk = rank_bits(rows, all_cols.len());
            (Some(rk), t1.elapsed().as_secs_f64())
        } else {
            (None, 0.0)
        };
        let detail: Vec<String> = block_ranks.iter().map(|(j, s, r)| format!("χ{j}×{s}: {r}")).collect();
        println!(
            "| {n} | {l} | {d} | {} | {} | {degree} | {} | {} | {} | {} | {} | {} | {:.1} (+{:.1} build) | {} | {:.1} | {} | {} |",
            lay.n_vars,
            eqs.len(),
            total_rows,
            total_cols.len(),
            reps.len(),
            free_cols.len(),
            rank_triv,
            detail.join(", "),
            blocks_s,
            build_s,
            full_rank.map(|r| r.to_string()).unwrap_or("—".into()),
            full_s,
            full_rank.map(|_| format!("{:.1}×", full_s / blocks_s.max(1e-9))).unwrap_or("—".into()),
            full_rank.map(|r| if r == sum { "HOLDS".to_string() } else { format!("FAILS (Σ={sum})") }).unwrap_or(format!("Σ={sum}"))
        );
    }
}
