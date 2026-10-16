//! G2 scoping probe: does the Frobenius-equivariant reduction of the barrel
//! decomposition system buy anything over gauge-fixing the target's shift?
//!
//! The barrel formulation puts summand `i` in `τ^{k_i}(W)` with a one-hot
//! selector `s_i[k]`, and — to make the system invariant under the cyclic
//! shift of all selectors — also gives the target a selector `s_0[k]` so
//! that `X_R = Σ_k s_0[k] R^{2^k}`. The shift `σ: k ↦ k + 1` acts freely on
//! the selector blocks, so the symmetric system is `n` copies of the
//! gauge-fixed system (`k_0 = 0`) glued by `σ`. Faugère–Svartz splits the
//! symmetric Macaulay matrix into blocks of `1/n` the size — but if those
//! blocks are no smaller than the gauge-fixed matrix, the reduction only
//! recovers what gauge-fixing gives for free, at extension-field cost.
//!
//! This probe measures exactly that, by pure counting, on toy Koblitz
//! instances (`K_0`, `n ∈ {7, 11, 13}`, window `W = ⟨1, z, …, z^{l−1}⟩`,
//! two summands): it builds both systems in a **normal basis** (so `σ` is a
//! permutation of variables and of equations), enumerates the Macaulay
//! rows and columns at degrees 2–4 exactly as `koblitz_groebner::build_macaulay`
//! does, counts `σ`-orbits for the symmetric system, and compares the orbit
//! counts (= equivariant block dimensions) with the gauge-fixed row and
//! column counts. It also reports `rank` for both systems from
//! `macaulay_profile`.
//!
//! Formulation (degree kept at 3 by auxiliary variables):
//!   u_i   = Σ_k s_i[k] · τ^k(y_i),   y_i = Σ_t c_i[t] w_t        (n quadratic constraints per summand)
//!   S₃(u_1, u_2, X_R) = 0                                         (n cubic equations)
//!   Σ_k s_i[k] = 1,  s_i[k] s_i[k'] = 0                            (one-hot, per block)
//!
//! ```bash
//! cargo run --release --example g2_equivariant_orbit_probe -- [--n 7,11,13] [--l 3] [--dmax 4]
//! ```

use std::collections::{HashMap, HashSet};
use std::time::Instant;

use crypto_lib::binary_ecc::f2m::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::koblitz_groebner::{
    macaulay_profile, sym_semaev_s3, FieldStructure, SymElement,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible;
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn bits_of(e: &F2mElement) -> u64 {
    e.raw_bits().first().copied().unwrap_or(0)
}

fn elem_from_bits(bits: u64, n: u32) -> F2mElement {
    let positions: Vec<u32> = (0..n).filter(|k| (bits >> k) & 1 == 1).collect();
    F2mElement::from_bit_positions(&positions, n)
}

/// Rank of a set of n-bit vectors over F_2.
fn rank_f2(vectors: &[u64]) -> usize {
    let mut rows: Vec<u64> = vectors.to_vec();
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

/// Find a normal element α (its conjugates α^{2^j} form a basis) and return
/// the conjugates in poly-basis coordinates plus the inverse change of basis
/// (rows: normal coordinate j as an F_2-linear form on poly coordinates).
fn normal_basis(n: u32, irr: &IrreduciblePoly, rng: &mut StdRng) -> (Vec<F2mElement>, Vec<u64>) {
    loop {
        let alpha = elem_from_bits(rng.gen_range(1..(1u64 << n)), n);
        let conj: Vec<F2mElement> = (0..n).map(|j| alpha.square_k_times(j, irr)).collect();
        let cols: Vec<u64> = conj.iter().map(bits_of).collect();
        if rank_f2(&cols) as u32 != n {
            continue;
        }
        // Invert the n×n matrix P whose column j is cols[j]: solve P a = v for v = e_i.
        // Gauss-Jordan on the augmented [P | I] with rows indexed by poly coordinate.
        let nn = n as usize;
        let mut rows: Vec<(u64, u64)> = (0..nn)
            .map(|i| {
                let mut p_row = 0u64;
                for (j, &c) in cols.iter().enumerate() {
                    if (c >> i) & 1 == 1 {
                        p_row |= 1 << j;
                    }
                }
                (p_row, 1u64 << i)
            })
            .collect();
        for col in 0..nn {
            let piv = (col..nn)
                .find(|&r| (rows[r].0 >> col) & 1 == 1)
                .expect("invertible");
            rows.swap(col, piv);
            let pv = rows[col];
            for r in 0..nn {
                if r != col && (rows[r].0 >> col) & 1 == 1 {
                    rows[r].0 ^= pv.0;
                    rows[r].1 ^= pv.1;
                }
            }
        }
        // rows[j].1 is now row j of P^{-1}: normal coordinate j = Σ_i (P^{-1})_{ji} v_i.
        let pinv: Vec<u64> = rows.iter().map(|r| r.1).collect();
        return (conj, pinv);
    }
}

/// Convert poly-basis coordinate polynomials into normal-basis coordinates.
fn to_normal_coords(coords: &[F2BoolPoly], pinv: &[u64], n_vars: usize) -> Vec<F2BoolPoly> {
    pinv.iter()
        .map(|row| {
            let mut acc = F2BoolPoly::zero(n_vars);
            for (i, c) in coords.iter().enumerate() {
                if (row >> i) & 1 == 1 {
                    acc = acc.add(c);
                }
            }
            acc
        })
        .collect()
}

/// Variable layout for one system.
struct Layout {
    n: usize,
    with_target_selector: bool,
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
    fn new(n: usize, l: usize, with_target_selector: bool) -> Self {
        let mut off = 0;
        let s0 = off;
        if with_target_selector {
            off += n;
        }
        let s1 = off;
        off += n;
        let s2 = off;
        off += n;
        let c1 = off;
        off += l;
        let c2 = off;
        off += l;
        let u1 = off;
        off += n;
        let u2 = off;
        off += n;
        Layout {
            n,
            with_target_selector,
            s0,
            s1,
            s2,
            c1,
            c2,
            u1,
            u2,
            n_vars: off,
        }
    }

    /// The cyclic shift σ on a monomial mask: selector blocks and u blocks rotate by one.
    fn shift_mask(&self, mask: u64) -> u64 {
        let n = self.n;
        let rot = |m: u64, base: usize| -> u64 {
            let block = (m >> base) & ((1u64 << n) - 1);
            let rotated = ((block << 1) | (block >> (n - 1))) & ((1u64 << n) - 1);
            (m & !(((1u64 << n) - 1) << base)) | (rotated << base)
        };
        let mut m = mask;
        if self.with_target_selector {
            m = rot(m, self.s0);
        }
        m = rot(m, self.s1);
        m = rot(m, self.s2);
        m = rot(m, self.u1);
        m = rot(m, self.u2);
        m
    }
}

/// Equation families, each with its own σ-action on the index.
#[derive(Clone, Copy, Debug)]
enum Fam {
    S3Coord,          // n, cyclic
    UConstraint(u8),  // n per summand, cyclic
    OneHotLinear(u8), // 1 per block, fixed
    OneHotQuad(u8),   // C(n,2) per block: (k, k') ↦ (k+1, k'+1)
}

struct System {
    eqs: Vec<F2BoolPoly>,
    fam: Vec<(Fam, usize, usize)>, // (family, index a, index b) for σ-action
}

fn build_system(
    n: u32,
    l: u32,
    irr: &IrreduciblePoly,
    st: &FieldStructure,
    b: &F2mElement,
    r: &F2mElement,
    normal: &[F2mElement],
    pinv: &[u64],
    with_target_selector: bool,
) -> (Layout, System) {
    let lay = Layout::new(n as usize, l as usize, with_target_selector);
    let nv = lay.n_vars;
    let nn = n as usize;
    let window: Vec<F2mElement> = (0..l)
        .map(|k| F2mElement::from_bit_positions(&[k], n))
        .collect();
    let s_elem = |var: usize| SymElement::from_subspace_vars(&[F2mElement::one(n)], var, n, nv);
    // u_i as symbolic elements in the normal basis (coordinates are the u variables).
    let u1 = SymElement::from_subspace_vars(normal, lay.u1, n, nv);
    let u2 = SymElement::from_subspace_vars(normal, lay.u2, n, nv);
    // Σ_k s_i[k] τ^k(y_i)
    let barrel = |s_off: usize, c_off: usize| -> SymElement {
        let mut acc = SymElement::zero(n, nv);
        for k in 0..nn {
            let shifted: Vec<F2mElement> = window
                .iter()
                .map(|w| w.square_k_times(k as u32, irr))
                .collect();
            let y_k = SymElement::from_subspace_vars(&shifted, c_off, n, nv);
            acc = acc.add(&s_elem(s_off + k).mul(&y_k, st));
        }
        acc
    };
    let x_r = if with_target_selector {
        let mut acc = SymElement::zero(n, nv);
        for k in 0..nn {
            let rk = SymElement::constant(&r.square_k_times(k as u32, irr), n, nv);
            acc = acc.add(&s_elem(lay.s0 + k).mul(&rk, st));
        }
        acc
    } else {
        SymElement::constant(r, n, nv)
    };
    let mut eqs = Vec::new();
    let mut fam = Vec::new();
    // S3 in normal coordinates
    for (j, e) in to_normal_coords(&sym_semaev_s3(&u1, &u2, &x_r, b, st), pinv, nv)
        .into_iter()
        .enumerate()
    {
        eqs.push(e);
        fam.push((Fam::S3Coord, j, 0));
    }
    // u-constraints in normal coordinates
    for (i, (u, s_off, c_off)) in [(&u1, lay.s1, lay.c1), (&u2, lay.s2, lay.c2)]
        .into_iter()
        .enumerate()
    {
        let diff = u.add(&barrel(s_off, c_off));
        for (j, e) in to_normal_coords(&diff.coords, pinv, nv)
            .into_iter()
            .enumerate()
        {
            eqs.push(e);
            fam.push((Fam::UConstraint(i as u8), j, 0));
        }
    }
    // one-hot constraints
    let mut blocks = vec![(1u8, lay.s1), (2u8, lay.s2)];
    if with_target_selector {
        blocks.push((0u8, lay.s0));
    }
    for (blk, off) in blocks {
        let mut monos: Vec<F2BoolMono> =
            (0..nn).map(|k| F2BoolMono::var((off + k) as u32)).collect();
        monos.push(F2BoolMono::one());
        eqs.push(F2BoolPoly::from_monos(monos, nv));
        fam.push((Fam::OneHotLinear(blk), 0, 0));
        for k in 0..nn {
            for k2 in (k + 1)..nn {
                let m = F2BoolMono::var((off + k) as u32).mul(F2BoolMono::var((off + k2) as u32));
                eqs.push(F2BoolPoly::from_monos(vec![m], nv));
                fam.push((Fam::OneHotQuad(blk), k, k2));
            }
        }
    }
    (lay, System { eqs, fam })
}

/// Monomials of degree ≤ d over nv variables, as masks.
fn monomials_up_to(nv: usize, d: u32, out: &mut Vec<u64>) {
    fn rec(start: usize, nv: usize, left: u32, cur: u64, out: &mut Vec<u64>) {
        out.push(cur);
        if left == 0 {
            return;
        }
        for v in start..nv {
            rec(v + 1, nv, left - 1, cur | (1u64 << v), out);
        }
    }
    rec(0, nv, d, 0, out);
}

/// Macaulay rows (multiplier mask, equation index) and the column set at degree d,
/// built exactly like `koblitz_groebner::build_macaulay` (parity cancellation).
fn macaulay_shape(sys: &System, nv: usize, d: u32) -> (Vec<(u64, usize)>, HashSet<u64>) {
    let mut rows = Vec::new();
    let mut cols = HashSet::new();
    for (ei, p) in sys.eqs.iter().enumerate() {
        let pdeg = p
            .terms
            .iter()
            .map(|t| t.mask.count_ones())
            .max()
            .unwrap_or(0);
        if pdeg > d {
            continue;
        }
        let mut mults = Vec::new();
        monomials_up_to(nv, d - pdeg, &mut mults);
        for m in mults {
            let mut all: Vec<u64> = p.terms.iter().map(|t| t.mask | m).collect();
            all.sort_unstable();
            let mut row_nonempty = false;
            let mut i = 0;
            while i < all.len() {
                let mut j = i;
                while j < all.len() && all[j] == all[i] {
                    j += 1;
                }
                if (j - i) % 2 == 1 {
                    cols.insert(all[i]);
                    row_nonempty = true;
                }
                i = j;
            }
            if row_nonempty {
                rows.push((m, ei));
            }
        }
    }
    (rows, cols)
}

/// σ on an equation index.
fn shift_eq(
    sys: &System,
    lay: &Layout,
    ei: usize,
    index: &HashMap<(u8, usize, usize, u8), usize>,
) -> usize {
    let n = lay.n;
    let (f, a, b2) = sys.fam[ei];
    let key = match f {
        Fam::S3Coord => (0u8, (a + 1) % n, 0usize, 0u8),
        Fam::UConstraint(i) => (1u8, (a + 1) % n, 0usize, i),
        Fam::OneHotLinear(blk) => (2u8, 0usize, 0usize, blk),
        Fam::OneHotQuad(blk) => {
            let (k, k2) = ((a + 1) % n, (b2 + 1) % n);
            (3u8, k.min(k2), k.max(k2), blk)
        }
    };
    index[&key]
}

fn fam_key(f: Fam, a: usize, b2: usize) -> (u8, usize, usize, u8) {
    match f {
        Fam::S3Coord => (0, a, 0, 0),
        Fam::UConstraint(i) => (1, a, 0, i),
        Fam::OneHotLinear(blk) => (2, 0, 0, blk),
        Fam::OneHotQuad(blk) => (3, a.min(b2), a.max(b2), blk),
    }
}

// ── Block-rank milestone: F_{2^d} arithmetic and the character blocks ──

/// Small binary field F_{2^d}, elements as u64 bit-polynomials, d ≤ 40.
#[derive(Clone)]
struct Gf {
    d: u32,
    /// Irreducible polynomial of degree d (bit d set).
    poly: u64,
    /// log/antilog tables (d ≤ 16): exp[i] = g^i, log[x] = i.
    exp: Vec<u32>,
    log: Vec<u32>,
}

impl Gf {
    fn new(d: u32) -> Self {
        // Search for an irreducible polynomial of degree d by Ben-Or's test.
        for low in 1u64..(1u64 << d) {
            let f = (1u64 << d) | low;
            if low & 1 == 0 {
                continue;
            }
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
        panic!("no irreducible polynomial of degree {d}");
    }
    /// Find a generator of the multiplicative group and tabulate it.
    fn build_tables(&mut self) {
        let order = (1u64 << self.d) - 1;
        // factor the order for the generator test
        let mut fs = Vec::new();
        let mut m = order;
        let mut q = 2u64;
        while q * q <= m {
            if m.is_multiple_of(q) {
                fs.push(q);
                while m.is_multiple_of(q) {
                    m /= q;
                }
            }
            q += 1;
        }
        if m > 1 {
            fs.push(m);
        }
        let mut g = 2u64;
        loop {
            if fs.iter().all(|&f| self.pow_slow(g, order / f) != 1) {
                break;
            }
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
    fn is_irreducible(f: u64, d: u32) -> bool {
        // gcd(x^{2^i} − x, f) = 1 for i = 1..d/2
        let g = Gf {
            d,
            poly: f,
            exp: Vec::new(),
            log: Vec::new(),
        };
        let mut x = 2u64; // x
        for _ in 1..=(d / 2) {
            x = g.mul_slow(x, x);
            let t = x ^ 2;
            if Self::poly_gcd(t, f) != 1 {
                return false;
            }
        }
        true
    }
    fn poly_gcd(mut a: u64, mut b: u64) -> u64 {
        while b != 0 {
            let (da, db) = (63 - a.leading_zeros(), 63 - b.leading_zeros());
            if a == 0 {
                return b;
            }
            if da < db {
                std::mem::swap(&mut a, &mut b);
                continue;
            }
            a ^= b << (da - db);
        }
        a
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
    #[inline]
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
    /// A primitive n-th root of unity (n | 2^d − 1).
    fn root_of_unity(&self, n: u64, rng: &mut StdRng) -> u64 {
        let order = (1u64 << self.d) - 1;
        assert_eq!(order % n, 0, "n must divide 2^d − 1");
        loop {
            let g = rng.gen_range(2..(1u64 << self.d));
            let z = self.pow(g, order / n);
            if z != 1 {
                return z;
            }
        }
    }
}

/// Rank of a dense matrix over F_{2^d}.
fn rank_gf(gf: &Gf, mut m: Vec<Vec<u64>>) -> usize {
    if m.is_empty() {
        return 0;
    }
    let cols = m[0].len();
    let mut rank = 0;
    for c in 0..cols {
        let Some(p) = (rank..m.len()).find(|&r| m[r][c] != 0) else {
            continue;
        };
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
                    *x ^= gf.mul(f, *pv);
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

/// Rank of a dense bit matrix over F_2 (rows as Vec<u64> words).
fn rank_bits(mut m: Vec<Vec<u64>>, cols: usize) -> usize {
    let mut rank = 0;
    for c in 0..cols {
        let (w, bit) = (c / 64, 1u64 << (c % 64));
        let Some(p) = (rank..m.len()).find(|&r| m[r][w] & bit != 0) else {
            continue;
        };
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

/// The columns (with parity cancellation) of one Macaulay row.
fn row_columns(p: &F2BoolPoly, m: u64) -> Vec<u64> {
    let mut all: Vec<u64> = p.terms.iter().map(|t| t.mask | m).collect();
    all.sort_unstable();
    let mut out = Vec::new();
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

/// Galois orbits of the nonzero characters j ∈ (Z/n)^* under j ↦ 2j.
fn character_orbits(n: usize) -> Vec<(usize, usize)> {
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for j in 1..n {
        if seen[j] {
            continue;
        }
        let mut k = j;
        let mut size = 0;
        while !seen[k] {
            seen[k] = true;
            size += 1;
            k = (2 * k) % n;
        }
        out.push((j, size));
    }
    out
}

/// Build the trivial block over F_2 and one character block per Galois orbit,
/// and return (rank_full, rank_trivial, Vec<(j, orbit size, rank_j)>, wall seconds for blocks).
fn block_ranks(
    sys: &System,
    lay: &Layout,
    d: u32,
    field_d: u32,
    index: &HashMap<(u8, usize, usize, u8), usize>,
    rng: &mut StdRng,
) -> (usize, usize, Vec<(usize, usize, usize)>, f64) {
    let n = lay.n;
    let nv = lay.n_vars;
    // Row representatives and their column sets.
    let (rows, _cols) = macaulay_shape(sys, nv, d);
    let mut reps: Vec<(u64, usize)> = Vec::new();
    for &(m, ei) in &rows {
        let mut canon = (m, ei);
        let (mut mm, mut ee) = (m, ei);
        for _ in 1..n {
            mm = lay.shift_mask(mm);
            ee = shift_eq(sys, lay, ee, index);
            canon = canon.min((mm, ee));
        }
        if canon == (m, ei) {
            reps.push((m, ei));
        }
    }
    // Column orbit bookkeeping: canonical representative, shift t, fixed?
    let mut col_info: HashMap<u64, (u64, u32, bool)> = HashMap::new();
    let mut col_index: HashMap<u64, usize> = HashMap::new(); // canonical → block column index
    let classify = |c: u64,
                    col_info: &mut HashMap<u64, (u64, u32, bool)>,
                    col_index: &mut HashMap<u64, usize>| {
        if col_info.contains_key(&c) {
            return;
        }
        let mut canon = c;
        let mut m = c;
        for _ in 1..n {
            m = lay.shift_mask(m);
            canon = canon.min(m);
        }
        let fixed = lay.shift_mask(canon) == canon;
        // shift t with σ^t(canon) = c
        let mut t = 0u32;
        let mut m = canon;
        while m != c {
            m = lay.shift_mask(m);
            t += 1;
        }
        let next = col_index.len();
        col_index.entry(canon).or_insert(next);
        col_info.insert(c, (canon, t, fixed));
    };
    let mut rep_cols: Vec<Vec<u64>> = Vec::with_capacity(reps.len());
    for &(m, ei) in &reps {
        let cs = row_columns(&sys.eqs[ei], m);
        for &c in &cs {
            classify(c, &mut col_info, &mut col_index);
        }
        rep_cols.push(cs);
    }
    let ncols = col_index.len();
    let t0 = Instant::now();
    // Trivial block over F_2: T[r̄][c̄] = parity of columns of r̄ in the orbit of c̄.
    let words = ncols.div_ceil(64);
    let mut t_rows: Vec<Vec<u64>> = Vec::with_capacity(reps.len());
    for cs in &rep_cols {
        let mut row = vec![0u64; words];
        for &c in cs {
            let (canon, _, _) = col_info[&c];
            let k = col_index[&canon];
            row[k / 64] ^= 1u64 << (k % 64);
        }
        t_rows.push(row);
    }
    let rank_triv = rank_bits(t_rows, ncols);
    // Character blocks over F_{2^d}: N_j[r̄][c̄] = Σ_t P[r̄, σ^t c̄] ζ^{−jt}, free orbits only.
    let gf = Gf::new(field_d);
    let zeta = gf.root_of_unity(n as u64, rng);
    let zeta_inv = gf.inv(zeta);
    // free column representatives
    let free_cols: Vec<u64> = {
        let mut v: Vec<u64> = col_info
            .values()
            .filter(|(_, _, fixed)| !fixed)
            .map(|(canon, _, _)| *canon)
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    let free_index: HashMap<u64, usize> =
        free_cols.iter().enumerate().map(|(i, &c)| (c, i)).collect();
    let mut block_results = Vec::new();
    for (j, size) in character_orbits(n) {
        let zj = gf.pow(zeta_inv, j as u64);
        let mut rows_k: Vec<Vec<u64>> = Vec::new();
        for (ri, &(m, _ei)) in reps.iter().enumerate() {
            // fixed rows contribute nothing to nontrivial blocks
            if lay.shift_mask(m) == m {
                let (_, ei) = reps[ri];
                if shift_eq(sys, lay, ei, index) == ei {
                    continue;
                }
            }
            let mut row = vec![0u64; free_cols.len()];
            for &c in &rep_cols[ri] {
                let (canon, t, fixed) = col_info[&c];
                if fixed {
                    continue;
                }
                let k = free_index[&canon];
                row[k] ^= gf.pow(zj, t as u64);
            }
            rows_k.push(row);
        }
        let r = rank_gf(&gf, rows_k);
        block_results.push((j, size, r));
    }
    let secs = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let rank_full = macaulay_profile(&sys.eqs, nv, d)
        .map(|p| p.rank)
        .unwrap_or(0);
    let full_secs = t1.elapsed().as_secs_f64();
    eprintln!(
        "timing n={} degree={}: full F_2 rank {:.3}s, equivariant blocks {:.3}s, ratio {:.1}x",
        n,
        d,
        full_secs,
        secs,
        full_secs / secs.max(1e-9)
    );
    (rank_full, rank_triv, block_results, secs)
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let mut ns = vec![7u32, 11, 13];
    let mut l = 3u32;
    let mut dmax = 4u32;
    let mut blocks = false;
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
            "--dmax" => {
                i += 1;
                dmax = argv[i].parse().unwrap_or(4);
            }
            "--blocks" => blocks = true,
            _ => {}
        }
        i += 1;
    }
    println!("# G2 orbit probe: equivariant block dimensions vs the gauge-fixed system\n");
    println!("| n | l | system | vars | eqs | degree | rows | cols | row σ-orbits | col σ-orbits | rank | seconds |");
    println!("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|");
    for &n in &ns {
        let irr = find_irreducible(n).expect("irreducible");
        let st = FieldStructure::new(n, &irr);
        let b = F2mElement::one(n);
        let mut rng = StdRng::seed_from_u64(0x62 + n as u64);
        let (normal, pinv) = normal_basis(n, &irr, &mut rng);
        let r = elem_from_bits(rng.gen_range(1..(1u64 << n)), n);
        for with_sel in [true, false] {
            let (lay, sys) = build_system(n, l, &irr, &st, &b, &r, &normal, &pinv, with_sel);
            if lay.n_vars > 64 {
                println!(
                    "| {n} | {l} | {} | {} | — | — | too many variables | | | | | |",
                    if with_sel { "symmetric" } else { "gauge-fixed" },
                    lay.n_vars
                );
                continue;
            }
            // Sanity: the symmetric system is σ-invariant as a set of equations.
            let index: HashMap<(u8, usize, usize, u8), usize> = sys
                .fam
                .iter()
                .enumerate()
                .map(|(i, &(f, a, b2))| (fam_key(f, a, b2), i))
                .collect();
            if with_sel {
                let mut invariant = true;
                for (ei, p) in sys.eqs.iter().enumerate() {
                    let target = &sys.eqs[shift_eq(&sys, &lay, ei, &index)];
                    let mut shifted: Vec<u64> =
                        p.terms.iter().map(|t| lay.shift_mask(t.mask)).collect();
                    shifted.sort_unstable();
                    let mut tgt: Vec<u64> = target.terms.iter().map(|t| t.mask).collect();
                    tgt.sort_unstable();
                    if shifted != tgt {
                        invariant = false;
                        break;
                    }
                }
                println!("<!-- n={n}: symmetric equation set is σ-invariant: {invariant} -->");
            }
            for d in 2..=dmax {
                let t0 = Instant::now();
                let (rows, cols) = macaulay_shape(&sys, lay.n_vars, d);
                if cols.len() > 3_000_000 {
                    println!(
                        "| {n} | {l} | {} | {} | {} | {d} | {} | {} | too large | | | |",
                        if with_sel { "symmetric" } else { "gauge-fixed" },
                        lay.n_vars,
                        sys.eqs.len(),
                        rows.len(),
                        cols.len()
                    );
                    break;
                }
                let (row_orbits, col_orbits) = if with_sel {
                    let mut seen_cols: HashSet<u64> = HashSet::new();
                    for &c in &cols {
                        let mut canon = c;
                        let mut m = c;
                        for _ in 1..n {
                            m = lay.shift_mask(m);
                            canon = canon.min(m);
                        }
                        seen_cols.insert(canon);
                    }
                    let mut seen_rows: HashSet<(u64, usize)> = HashSet::new();
                    for &(m, ei) in &rows {
                        let mut canon = (m, ei);
                        let (mut mm, mut ee) = (m, ei);
                        for _ in 1..n {
                            mm = lay.shift_mask(mm);
                            ee = shift_eq(&sys, &lay, ee, &index);
                            canon = canon.min((mm, ee));
                        }
                        seen_rows.insert(canon);
                    }
                    (seen_rows.len(), seen_cols.len())
                } else {
                    (rows.len(), cols.len())
                };
                let rank = macaulay_profile(&sys.eqs, lay.n_vars, d)
                    .map(|p| p.rank.to_string())
                    .unwrap_or("—".into());
                if blocks && with_sel && d <= 3 {
                    let dd = {
                        let mut k = 1u32;
                        let mut v = 2u64 % n as u64;
                        while v != 1 {
                            v = v * 2 % n as u64;
                            k += 1;
                        }
                        k
                    };
                    let (rf, rt, bl, secs) = block_ranks(&sys, &lay, d, dd, &index, &mut rng);
                    let sum: usize = rt + bl.iter().map(|(_, size, r)| size * r).sum::<usize>();
                    let detail: Vec<String> = bl
                        .iter()
                        .map(|(j, size, r)| format!("χ_{j} (orbit {size}): rank {r}"))
                        .collect();
                    println!(
                        "<!-- blocks n={n} d={d}: full rank {rf} | trivial rank {rt} | {} | Σ = {sum} | identity {} | field F_2^{dd} | {:.1}s -->",
                        detail.join(", "),
                        if sum == rf { "HOLDS" } else { "FAILS" },
                        secs
                    );
                }
                println!(
                    "| {n} | {l} | {} | {} | {} | {d} | {} | {} | {} | {} | {} | {:.1} |",
                    if with_sel { "symmetric" } else { "gauge-fixed" },
                    lay.n_vars,
                    sys.eqs.len(),
                    rows.len(),
                    cols.len(),
                    row_orbits,
                    col_orbits,
                    rank,
                    t0.elapsed().as_secs_f64()
                );
            }
        }
    }
}
