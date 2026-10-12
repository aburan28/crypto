//! E5 of `research/fall_degree_bounds_20261008/PREREGISTRATION.md`: a native
//! engine for the Huang–Kosters–Yeo last fall degree `d_F` of the Weil descent
//! of Semaev's `S₃` over an `F₂`-subspace `V ⊂ F_{2^n}` (`m = 2`).
//!
//! It replaces `lfd_fast.sage` and keeps its conventions exactly:
//!
//! - curve `y² + xy = x³ + b`, target `x₃`, factor base `X₁, X₂ ∈ V`,
//!   `dim V = n' = ⌊n/2⌋`, `N = 2n'` Boolean unknowns, the `n` coordinate
//!   equations of `S₃(X₁, X₂, x₃) = (X₁X₂ + X₁x₃ + X₂x₃)² + X₁X₂x₃ + b`;
//! - the Boolean ring with the rule "multiply `g` by a monomial `m` only if
//!   `deg m + deg g ≤ c`" (the ledger's §0 convention);
//! - `V_{F,c}` = the mutant closure of the degree-`c` Macaulay rows; `e₀` =
//!   the first `c ≥ 2` with `dim V_{F,c} = dim (I ∩ R_{≤c})` and `c ≥ DG`;
//!   `d_last` = the least `c ≤ e₀` from which every `c' ∈ [c, e₀]` agrees.
//!
//! What it does differently, and why it is exact:
//!
//! - **No Gröbner basis.** With the field equations the Boolean ideal is
//!   radical, so `I = I(S)` for the solution set `S ⊂ F₂^N`, and
//!   `dim (I ∩ R_{≤c}) = #monomials(≤ c) − rank(evaluation of those
//!   monomials on S)`. `S` is found exactly by solving, for each `X₁ ∈ V`,
//!   the quadratic `S₃ = 0` in `X₂` over `F_{2^n}` and testing `X₂ ∈ V`
//!   (`2^{n'}` quadratic solves).
//! - **`DG`** (the top degree of the reduced deglex Gröbner basis of `I(S)`,
//!   `x₀ > x₁ > …`, as in Sage's `BooleanPolynomialRing(order='deglex')`) is
//!   read off the standard monomials of `S` (Buchberger–Möller order scan):
//!   the reduced basis's leading monomials are the minimal non-standard
//!   monomials. `S = ∅` gives `I = (1)` and `DG = 0`.
//! - **Closure** by a work queue: every row that enters the echelon basis
//!   with degree `< c` has its `N` variable multiples queued. Rows keep their
//!   leading (deglex-largest) monomial as pivot, so `V ∩ R_{≤d}` is spanned
//!   by the rows of degree `≤ d` and the queue reaches the same fixpoint as
//!   the Sage loop. The closure stops early once it reaches `dim (I ∩ R_{≤c})`,
//!   since `V_{F,c} ⊆ I`.
//!
//! Field: `F_{2^n}` with Sage's `modulus='minimal_weight'`, which for
//! `GF(2)` is NTL's sparse irreducible (the first trinomial `t^n + t^k + 1`
//! by increasing `k`, else the first pentanomial by increasing
//! `(k₃, k₂, k₁)`); `--modulus HEX` overrides it.
//!
//! ```bash
//! # reproduce a committed Sage record (subfield cells: V is fixed by n)
//! hky_last_fall_degree replay research/fall_degree_bounds_20261008/runs/n10_subfield_rand.jsonl
//! # a new cell: n, family, target, draws, seed, cmax
//! hky_last_fall_degree cell --n 24 --family random --target rand --draws 8 --seed 20261012 --cmax 6
//! ```
use rayon::prelude::*;
use serde_json::{json, Value};
use std::collections::HashMap;
use std::time::Instant;

// ------------------------------------------------------------------ F_{2^n}

#[derive(Clone, Copy)]
pub struct Field {
    pub n: u32,
    pub poly: u64,
}

impl Field {
    pub fn mul(&self, mut a: u64, mut b: u64) -> u64 {
        let mut r = 0;
        while b != 0 {
            if b & 1 == 1 {
                r ^= a;
            }
            b >>= 1;
            a <<= 1;
            if a >> self.n & 1 == 1 {
                a ^= self.poly;
            }
        }
        r
    }
    pub fn sqr(&self, a: u64) -> u64 {
        self.mul(a, a)
    }
    fn inv(&self, a: u64) -> u64 {
        let (mut r, mut e, mut base) = (1u64, (1u64 << self.n) - 2, a);
        while e != 0 {
            if e & 1 == 1 {
                r = self.mul(r, base);
            }
            base = self.sqr(base);
            e >>= 1;
        }
        r
    }
    fn mask(&self) -> u64 {
        (1u64 << self.n) - 1
    }
}

fn poly_mod(mut a: u64, b: u64) -> u64 {
    let db = 63 - b.leading_zeros();
    while a != 0 && 63 - a.leading_zeros() >= db {
        a ^= b << (63 - a.leading_zeros() - db);
    }
    a
}

fn poly_gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = poly_mod(a, b);
        a = b;
        b = r;
    }
    a
}

fn irreducible(p: u64, n: u32) -> bool {
    let f = Field { n, poly: p };
    let mut t = 2u64;
    for _ in 1..=n / 2 {
        t = f.sqr(t);
        if poly_gcd(p, t ^ 2) != 1 {
            return false;
        }
    }
    true
}

/// NTL `BuildSparseIrred`, which Sage's `minimal_weight` uses over `GF(2)`.
fn sparse_irreducible(n: u32) -> u64 {
    for k in 1..n {
        let p = (1u64 << n) | (1 << k) | 1;
        if irreducible(p, n) {
            return p;
        }
    }
    for k3 in 3..n {
        for k2 in 2..k3 {
            for k1 in 1..k2 {
                let p = (1u64 << n) | (1 << k3) | (1 << k2) | (1 << k1) | 1;
                if irreducible(p, n) {
                    return p;
                }
            }
        }
    }
    panic!("no sparse irreducible of degree {n}")
}

/// Parse Sage's `z^4 + z^3 + 1` notation.
fn parse_elem(s: &str) -> u64 {
    let mut v = 0u64;
    for term in s.split('+').map(str::trim) {
        v ^= match term {
            "1" => 1,
            "0" => 0,
            "z" => 2,
            t => 1u64 << t.strip_prefix("z^").expect("term").parse::<u32>().unwrap(),
        };
    }
    v
}

/// Solutions `t` of `t² + t = c` (all of them), by linear algebra on the
/// `F₂`-linear map `t ↦ t² + t`.
pub struct ArtinSchreier {
    /// rows (image, source) in echelon form on the image
    rows: Vec<(u64, u64)>,
    kernel: Vec<u64>,
}

impl ArtinSchreier {
    pub fn new(f: &Field) -> Self {
        let n = f.n;
        let mut rows: Vec<(u64, u64)> = (0..n)
            .map(|i| {
                let t = 1u64 << i;
                (f.sqr(t) ^ t, t)
            })
            .collect();
        let mut basis: Vec<(u64, u64)> = vec![];
        let mut kernel = vec![];
        for (mut img, mut src) in rows.drain(..) {
            for &(bi, bs) in &basis {
                if img & (1u64 << (63 - bi.leading_zeros())) != 0 {
                    img ^= bi;
                    src ^= bs;
                }
            }
            if img == 0 {
                kernel.push(src);
            } else {
                basis.push((img, src));
                basis.sort_by(|a, b| b.0.cmp(&a.0));
            }
        }
        ArtinSchreier { rows: basis, kernel }
    }
    fn solve(&self, mut c: u64) -> Vec<u64> {
        let mut t = 0u64;
        for &(bi, bs) in &self.rows {
            if c & (1u64 << (63 - bi.leading_zeros())) != 0 {
                c ^= bi;
                t ^= bs;
            }
        }
        if c != 0 {
            return vec![];
        }
        // kernel of t² + t is F₂ = {0, 1}
        let mut out = vec![t];
        for &k in &self.kernel {
            let more: Vec<u64> = out.iter().map(|&x| x ^ k).collect();
            out.extend(more);
        }
        out
    }
}

/// All roots in `F_{2^n}` of `α T² + β T + γ`.
pub fn quad_roots(f: &Field, asl: &ArtinSchreier, a: u64, b: u64, c: u64) -> Vec<u64> {
    match (a, b) {
        (0, 0) => {
            if c == 0 {
                panic!("degenerate quadratic: every T is a root")
            } else {
                vec![]
            }
        }
        (0, b) => vec![f.mul(c, f.inv(b))],
        (a, 0) => {
            // T² = c/a; square root = x^(2^(n−1))
            let mut r = f.mul(c, f.inv(a));
            for _ in 0..f.n - 1 {
                r = f.sqr(r);
            }
            vec![r]
        }
        (a, b) => {
            // T = (b/a)·t, t² + t = a·c / b²
            let s = f.mul(b, f.inv(a));
            let rhs = f.mul(f.mul(a, c), f.inv(f.sqr(b)));
            asl.solve(rhs).into_iter().map(|t| f.mul(s, t)).collect()
        }
    }
}

// ------------------------------------------------------------------ subspace

struct Subspace {
    basis: Vec<u64>,
    /// echelon rows (vector, coordinate combination) for membership
    ech: Vec<(u64, u64)>,
}

impl Subspace {
    fn new(basis: Vec<u64>) -> Self {
        let mut ech: Vec<(u64, u64)> = vec![];
        for (j, &v) in basis.iter().enumerate() {
            let (mut x, mut c) = (v, 1u64 << j);
            for &(e, ec) in &ech {
                if x & (1u64 << (63 - e.leading_zeros())) != 0 {
                    x ^= e;
                    c ^= ec;
                }
            }
            assert!(x != 0, "basis is dependent");
            ech.push((x, c));
            ech.sort_by(|a, b| b.0.cmp(&a.0));
        }
        Subspace { basis, ech }
    }
    /// Coordinates of `x` in `basis`, or `None` if `x ∉ V`.
    fn coords(&self, mut x: u64) -> Option<u64> {
        let mut c = 0u64;
        for &(e, ec) in &self.ech {
            if x & (1u64 << (63 - e.leading_zeros())) != 0 {
                x ^= e;
                c ^= ec;
            }
        }
        (x == 0).then_some(c)
    }
    fn elem(&self, coords: u64) -> u64 {
        self.basis
            .iter()
            .enumerate()
            .filter(|(j, _)| coords >> j & 1 == 1)
            .fold(0, |acc, (_, &v)| acc ^ v)
    }
}

// ------------------------------------------------------------------ the system

/// The `n` coordinate equations, each a list of multilinear monomials (bit
/// masks over `N = 2n'` variables: `x_j` is bit `j`, `y_k` is bit `n' + k`).
fn descend(f: &Field, v: &Subspace, x3: u64, b: u64) -> Vec<Vec<u64>> {
    let np = v.basis.len();
    let n = f.n as usize;
    let mut eqs: Vec<HashMap<u64, ()>> = vec![HashMap::new(); n];
    let toggle = |coeff: u64, mono: u64, eqs: &mut Vec<HashMap<u64, ()>>| {
        for (i, eq) in eqs.iter_mut().enumerate() {
            if coeff >> i & 1 == 1 {
                if eq.remove(&mono).is_none() {
                    eq.insert(mono, ());
                }
            }
        }
    };
    // S₃ = Σ_{j,k} [(v_j v_k)² + v_j v_k x₃] x_j y_k + Σ_j (v_j x₃)² (x_j + y_j) + b
    for j in 0..np {
        for k in 0..np {
            let p = f.mul(v.basis[j], v.basis[k]);
            let coeff = f.sqr(p) ^ f.mul(p, x3);
            toggle(coeff, (1u64 << j) | (1u64 << (np + k)), &mut eqs);
        }
        let c = f.sqr(f.mul(v.basis[j], x3));
        toggle(c, 1u64 << j, &mut eqs);
        toggle(c, 1u64 << (np + j), &mut eqs);
    }
    toggle(b, 0, &mut eqs);
    eqs.into_iter()
        .map(|e| e.into_keys().collect::<Vec<u64>>())
        .collect()
}

/// The solution set, as variable assignments (bit masks).
fn solutions(f: &Field, asl: &ArtinSchreier, v: &Subspace, x3: u64, b: u64) -> Vec<u64> {
    let np = v.basis.len();
    let mut out = vec![];
    for c1 in 0..(1u64 << np) {
        let x1 = v.elem(c1);
        // (A X₂ + C)² + X₁x₃ X₂ + b with A = X₁ + x₃, C = X₁x₃
        let a = x1 ^ x3;
        let c = f.mul(x1, x3);
        for x2 in quad_roots(f, asl, f.sqr(a), c, f.sqr(c) ^ b) {
            if let Some(c2) = v.coords(x2) {
                out.push(c1 | (c2 << np));
            }
        }
    }
    out.sort_unstable();
    out
}

// ------------------------------------------------------------------ monomials

/// Multilinear monomials of degree `≤ c` in `nv` variables, in Sage's deglex
/// order (`x₀ > x₁ > …`), ascending.
pub fn monomials_deglex(nv: u32, c: u32) -> Vec<u64> {
    let mut out = vec![];
    for d in 0..=c {
        let mut level = vec![];
        fn rec(start: u32, nv: u32, left: u32, acc: u64, out: &mut Vec<u64>) {
            if left == 0 {
                out.push(acc);
                return;
            }
            for i in start..nv {
                rec(i + 1, nv, left - 1, acc | 1 << i, out);
            }
        }
        rec(0, nv, d, 0, &mut level);
        // ascending deglex within a degree: m < m' iff at the first variable
        // (x₀ first) where they differ, m' contains it. As bit masks with
        // x₀ = bit 0 this is the order of the bit-reversed masks.
        level.sort_by_key(|&m| m.reverse_bits());
        out.extend(level);
    }
    out
}

pub fn eval(mono: u64, point: u64) -> bool {
    point & mono == mono
}

/// `dim (I(S) ∩ R_{≤c})` for each `c ≤ cmax`, and `DG`.
pub fn ideal_side(nv: u32, sols: &[u64], cmax: u32) -> (Vec<usize>, u32) {
    let monos = monomials_deglex(nv, cmax);
    let s = sols.len();
    if s == 0 {
        let dims = (0..=cmax)
            .map(|c| monos.iter().filter(|m| m.count_ones() <= c).count())
            .collect();
        return (dims, 0);
    }
    // Buchberger–Möller scan in ascending deglex: a monomial is standard iff
    // its evaluation vector on S is independent of the earlier standard ones.
    assert!(s <= 128, "solution set too large for the u128 evaluation basis");
    let mut basis: Vec<u128> = vec![];
    let mut standard: Vec<u64> = vec![];
    for &m in &monos {
        let mut v: u128 = sols
            .iter()
            .enumerate()
            .fold(0, |acc, (i, &p)| acc | ((eval(m, p) as u128) << i));
        for &b in &basis {
            v = v.min(v ^ b);
        }
        if v != 0 {
            basis.push(v);
            basis.sort_unstable_by(|a, b| b.cmp(a));
            standard.push(m);
        }
    }
    let std_set: std::collections::HashSet<u64> = standard.iter().copied().collect();
    let dims = (0..=cmax)
        .map(|c| {
            monos.iter().filter(|m| m.count_ones() <= c).count()
                - standard.iter().filter(|m| m.count_ones() <= c).count()
        })
        .collect();
    // Minimal non-standard monomials (all proper divisors standard). Their
    // degree is at most max standard degree + 1 ≤ ⌊log₂ s⌋ + 1 ≤ cmax here.
    let max_std = standard.iter().map(|m| m.count_ones()).max().unwrap_or(0);
    let dg = monomials_deglex(nv, max_std + 1)
        .into_iter()
        .filter(|m| !std_set.contains(m))
        .filter(|&m| (0..nv).filter(|i| m >> i & 1 == 1).all(|i| std_set.contains(&(m & !(1 << i)))))
        .map(|m| m.count_ones())
        .max()
        .unwrap_or(0);
    (dims, dg)
}

// ------------------------------------------------------------------ closure

/// An echelon basis keyed by each row's leading (deglex-largest) monomial.
struct Echelon {
    words: usize,
    pivots: HashMap<usize, Vec<u64>>,
}

impl Echelon {
    fn new(cols: usize) -> Self {
        Echelon {
            words: cols.div_ceil(64),
            pivots: HashMap::new(),
        }
    }
    fn top(v: &[u64]) -> Option<usize> {
        v.iter()
            .enumerate()
            .rev()
            .find(|(_, &w)| w != 0)
            .map(|(i, &w)| i * 64 + 63 - w.leading_zeros() as usize)
    }
    /// Reduce and insert; returns the new row's leading column if it was new.
    fn insert(&mut self, mut v: Vec<u64>) -> Option<usize> {
        while let Some(t) = Self::top(&v) {
            match self.pivots.get(&t) {
                Some(row) => {
                    let w = t / 64 + 1;
                    for (x, &y) in v[..w].iter_mut().zip(&row[..w]) {
                        *x ^= y;
                    }
                }
                None => {
                    self.pivots.insert(t, v);
                    return Some(t);
                }
            }
        }
        None
    }
    fn dim(&self) -> usize {
        self.pivots.len()
    }
}

struct Ring {
    nv: u32,
    monos: Vec<u64>,
    index: HashMap<u64, usize>,
    words: usize,
}

impl Ring {
    fn new(nv: u32, c: u32) -> Self {
        let monos = monomials_deglex(nv, c);
        let index = monos.iter().enumerate().map(|(i, &m)| (m, i)).collect();
        let words = monos.len().div_ceil(64);
        Ring {
            nv,
            monos,
            index,
            words,
        }
    }
    fn row(&self, poly: &[u64]) -> Vec<u64> {
        let mut r = vec![0u64; self.words];
        for &m in poly {
            let i = self.index[&m];
            r[i / 64] ^= 1 << (i % 64);
        }
        r
    }
    fn degree_of_col(&self, col: usize) -> u32 {
        self.monos[col].count_ones()
    }
    /// `x_i · row` in the Boolean ring; `None` if it leaves degree `≤ c`.
    fn times_var(&self, row: &[u64], i: u32) -> Option<Vec<u64>> {
        let mut out = vec![0u64; self.words];
        for (w, &word) in row.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                let b = bits.trailing_zeros() as usize;
                bits &= bits - 1;
                let m = self.monos[w * 64 + b] | (1u64 << i);
                let j = *self.index.get(&m)?;
                out[j / 64] ^= 1 << (j % 64);
            }
        }
        Some(out)
    }
}

/// `(dim V_{F,c} without mutants, dim V_{F,c})`, the second capped at `target`.
pub fn closure_dims(eqs: &[Vec<u64>], nv: u32, c: u32, target: usize) -> (usize, usize) {
    let ring = Ring::new(nv, c);
    let mut ech = Echelon::new(ring.monos.len());
    let _ = ech.words;
    // Macaulay rows m·f with deg m ≤ c − deg f
    let mut low_rows: Vec<usize> = vec![];
    for f in eqs {
        let df = f.iter().map(|m| m.count_ones()).max().unwrap_or(0);
        if df > c {
            continue;
        }
        let base = ring.row(f);
        for &m in ring.monos.iter().filter(|m| m.count_ones() <= c - df) {
            let mut r = base.clone();
            // multiply by the monomial m, one variable at a time
            let mut ok = true;
            for i in 0..nv {
                if m >> i & 1 == 1 {
                    match ring.times_var(&r, i) {
                        Some(x) => r = x,
                        None => {
                            ok = false;
                            break;
                        }
                    }
                }
            }
            if ok {
                if let Some(t) = ech.insert(r) {
                    if ring.degree_of_col(t) < c {
                        low_rows.push(t);
                    }
                }
            }
        }
    }
    let plain = ech.dim();
    // mutant closure: queue variable multiples of every row of degree < c
    let mut queue = low_rows;
    while let Some(t) = queue.pop() {
        if ech.dim() >= target {
            break;
        }
        // 1 ∈ V (pivot at the constant monomial, column 0): V_{F,c} is everything
        if ech.pivots.contains_key(&0) {
            return (plain, ring.monos.len().min(target));
        }
        let row = ech.pivots[&t].clone();
        let products: Vec<Vec<u64>> = (0..nv)
            .into_par_iter()
            .filter_map(|i| ring.times_var(&row, i))
            .collect();
        for p in products {
            if let Some(t2) = ech.insert(p) {
                if ring.degree_of_col(t2) < c {
                    queue.push(t2);
                }
            }
        }
    }
    if ech.pivots.contains_key(&0) {
        return (plain, ring.monos.len().min(target));
    }
    (plain, ech.dim())
}

// ------------------------------------------------------------------ measurement

struct Draw {
    n: u32,
    family: String,
    target: String,
    seed: u64,
    field: Field,
    v: Subspace,
    x3: u64,
    b: u64,
}

fn measure(d: &Draw, cmax: u32) -> Value {
    let t0 = Instant::now();
    let f = &d.field;
    let np = d.v.basis.len() as u32;
    let nv = 2 * np;
    let asl = ArtinSchreier::new(f);
    let eqs = descend(f, &d.v, d.x3, d.b);
    let sols = solutions(f, &asl, &d.v, d.x3, d.b);
    // cross-check: every solution satisfies every descended equation
    for &p in &sols {
        for e in &eqs {
            assert!(!e.iter().fold(false, |acc, &m| acc ^ eval(m, p)), "solution check");
        }
    }
    let (idim, dg) = ideal_side(nv, &sols, cmax);
    let mut plain_dims = serde_json::Map::new();
    let mut clos_dims = serde_json::Map::new();
    let mut ideal_dims = serde_json::Map::new();
    let (mut ffd, mut prev_plain, mut e0) = (None, 0usize, None);
    let mut clo = HashMap::new();
    for c in 2..=cmax {
        let (plain, closed) = closure_dims(&eqs, nv, c, idim[c as usize]);
        // operational first fall: plain Macaulay space gains a row of degree < c
        // (Sage counts rows of degree < c; equivalently dim(plain ∩ R_{<c}) grows
        // past the previous plain dimension)
        let low = plain_low(&eqs, nv, c);
        if ffd.is_none() && low > prev_plain {
            ffd = Some(c);
        }
        prev_plain = plain;
        plain_dims.insert(c.to_string(), json!(plain));
        clos_dims.insert(c.to_string(), json!(closed));
        ideal_dims.insert(c.to_string(), json!(idim[c as usize]));
        clo.insert(c, closed);
        if e0.is_none() && closed == idim[c as usize] && c >= dg {
            e0 = Some(c);
            break;
        }
    }
    let d_last = e0.map(|e| {
        let mut dl = e;
        for c in (2..e).rev() {
            if clo.get(&c) == Some(&idim[c as usize]) {
                dl = c;
            } else {
                break;
            }
        }
        dl
    });
    json!({
        "schema": "hky_last_fall_degree/v1",
        "n": d.n, "nprime": np, "family": d.family, "target": d.target, "seed": d.seed,
        "modulus": format!("{:#x}", f.poly),
        "v_basis": d.v.basis.iter().map(|x| format!("{x:#x}")).collect::<Vec<_>>(),
        "x3": format!("{:#x}", d.x3), "b": format!("{:#x}", d.b),
        "nsol": sols.len(), "unsat": sols.is_empty(), "DG": dg,
        "ffd": ffd, "d_last": d_last, "e0": e0,
        "plain_dims": plain_dims, "closure_dims": clos_dims, "ideal_dims": ideal_dims,
        "secs": (t0.elapsed().as_secs_f64() * 10.0).round() / 10.0,
    })
}

/// Rows of degree `< c` in the plain (no-mutant) degree-`c` Macaulay space,
/// as Sage's `low` count: `dim(plain ∩ R_{≤c−1})`.
fn plain_low(eqs: &[Vec<u64>], nv: u32, c: u32) -> usize {
    let ring = Ring::new(nv, c);
    let mut ech = Echelon::new(ring.monos.len());
    for f in eqs {
        let df = f.iter().map(|m| m.count_ones()).max().unwrap_or(0);
        if df > c {
            continue;
        }
        for &m in ring.monos.iter().filter(|m| m.count_ones() <= c - df) {
            let mut r = ring.row(f);
            let mut ok = true;
            for i in 0..nv {
                if m >> i & 1 == 1 {
                    match ring.times_var(&r, i) {
                        Some(x) => r = x,
                        None => {
                            ok = false;
                            break;
                        }
                    }
                }
            }
            if ok {
                ech.insert(r);
            }
        }
    }
    ech.pivots.keys().filter(|&&t| ring.degree_of_col(t) < c).count()
}

// ------------------------------------------------------------------ draws

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

fn subfield_basis(f: &Field, np: u32) -> Vec<u64> {
    // Sage enumerates the subfield; any basis spans the same V, and every
    // dimension measured here is invariant under a change of basis.
    let mut basis = vec![];
    let mut sp = Subspace::new(vec![]);
    for e in 1..=f.mask() {
        if (0..np).fold(e, |y, _| f.sqr(y)) == e && sp.coords(e).is_none() {
            basis.push(e);
            sp = Subspace::new(basis.clone());
            if basis.len() == np as usize {
                break;
            }
        }
    }
    basis
}

fn new_draw(n: u32, family: &str, target: &str, seed: u64, poly: u64) -> Draw {
    let f = Field { n, poly };
    let np = n / 2;
    let mut rng = Rng(seed ^ (n as u64) << 40);
    let v = match family {
        "subfield" => Subspace::new(subfield_basis(&f, np)),
        _ => loop {
            let b: Vec<u64> = (0..np).map(|_| rng.next() & f.mask()).collect();
            let mut sp = Subspace::new(vec![]);
            let mut ok = true;
            for (j, &x) in b.iter().enumerate() {
                if x == 0 || sp.coords(x).is_some() {
                    ok = false;
                    break;
                }
                sp = Subspace::new(b[..=j].to_vec());
            }
            if ok {
                break Subspace::new(b);
            }
        },
    };
    let b = loop {
        let x = rng.next() & f.mask();
        if x != 0 {
            break x;
        }
    };
    let x3 = if target == "sat" {
        let asl = ArtinSchreier::new(&f);
        loop {
            let v1 = v.elem(rng.next() & ((1 << np) - 1));
            let v2 = v.elem(rng.next() & ((1 << np) - 1));
            // (v₁v₂ + (v₁+v₂)T)² + v₁v₂T + b = (v₁+v₂)²T² + v₁v₂T + (v₁v₂)² + b
            let p = f.mul(v1, v2);
            let s = v1 ^ v2;
            if s == 0 && p == 0 {
                continue;
            }
            let roots = quad_roots(&f, &asl, f.sqr(s), p, f.sqr(p) ^ b);
            if !roots.is_empty() {
                break roots[(rng.next() % roots.len() as u64) as usize];
            }
        }
    } else {
        rng.next() & f.mask()
    };
    Draw {
        n,
        family: family.into(),
        target: target.into(),
        seed,
        field: f,
        v,
        x3,
        b,
    }
}

fn arg<T: std::str::FromStr>(args: &[String], name: &str, default: Option<T>) -> T {
    match args.iter().position(|a| a == name) {
        Some(i) => args[i + 1].parse().ok().expect("bad argument"),
        None => default.unwrap_or_else(|| panic!("missing {name}")),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("replay") => {
            // Reproduce subfield records of lfd_fast.sage: V is fixed by n.
            let text = std::fs::read_to_string(&args[2]).expect("read records");
            for line in text.lines().filter(|l| !l.trim().is_empty()) {
                let rec: Value = serde_json::from_str(line).unwrap();
                let n = rec["n"].as_u64().unwrap() as u32;
                let family = rec["family"].as_str().unwrap();
                if family != "subfield" {
                    println!("{}", json!({"skip": "random V not stored", "seed": rec["seed"]}));
                    continue;
                }
                let poly = args
                    .iter()
                    .position(|a| a == "--modulus")
                    .map(|i| u64::from_str_radix(args[i + 1].trim_start_matches("0x"), 16).unwrap())
                    .unwrap_or_else(|| sparse_irreducible(n));
                let f = Field { n, poly };
                let d = Draw {
                    n,
                    family: family.into(),
                    target: if rec["unsat"].as_bool() == Some(true) { "rand".into() } else { "?".into() },
                    seed: rec["seed"].as_u64().unwrap(),
                    field: f,
                    v: Subspace::new(subfield_basis(&f, n / 2)),
                    x3: parse_elem(rec["x3"].as_str().unwrap()),
                    b: parse_elem(rec["b"].as_str().unwrap()),
                };
                let cmax = rec["closure_dims"]
                    .as_object()
                    .unwrap()
                    .keys()
                    .map(|k| k.parse::<u32>().unwrap())
                    .max()
                    .unwrap();
                let mine = measure(&d, cmax);
                let mut agree = true;
                for key in ["DG", "d_last", "e0", "unsat"] {
                    agree &= mine[key] == rec[key];
                }
                for key in ["closure_dims", "ideal_dims", "plain_dims"] {
                    for (c, v) in mine[key].as_object().unwrap() {
                        agree &= rec[key].get(c).map_or(true, |r| r == v);
                    }
                }
                println!(
                    "{}",
                    json!({"seed": rec["seed"], "n": n, "agree": agree,
                           "sage": {"DG": rec["DG"], "d_last": rec["d_last"], "e0": rec["e0"], "closure_dims": rec["closure_dims"], "ideal_dims": rec["ideal_dims"]},
                           "native": {"DG": mine["DG"], "d_last": mine["d_last"], "e0": mine["e0"], "closure_dims": mine["closure_dims"], "ideal_dims": mine["ideal_dims"], "nsol": mine["nsol"]}})
                );
            }
        }
        Some("cell") => {
            let n: u32 = arg(&args, "--n", None);
            let family: String = arg(&args, "--family", Some("random".to_string()));
            let target: String = arg(&args, "--target", Some("rand".to_string()));
            let draws: u64 = arg(&args, "--draws", Some(8));
            let seed: u64 = arg(&args, "--seed", None);
            let cmax: u32 = arg(&args, "--cmax", Some(6));
            let poly = sparse_irreducible(n);
            for s in 0..draws {
                let d = new_draw(n, &family, &target, seed.wrapping_mul(1000).wrapping_add(s), poly);
                println!("{}", measure(&d, cmax));
            }
        }
        Some("ladder") => {
            // E4: d_F (= PC refutation degree) on committed unsat dreg_ladder draws.
            use crypto_lib::binary_ecc::F2mElement;
            use crypto_lib::cryptanalysis::koblitz_groebner::{build_decomposition_system, FieldStructure};
            use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
            let m: usize = arg(&args, "--m", Some(3));
            let cmax: u32 = arg(&args, "--cmax", Some(6));
            let text = std::fs::read_to_string(&args[2]).expect("read ladder records");
            for line in text.lines().filter(|l| !l.trim().is_empty()) {
                let rec: Value = serde_json::from_str(line).unwrap();
                if rec.get("draw").is_none() || rec["solutions"].as_u64() != Some(0) {
                    continue;
                }
                let t0 = Instant::now();
                let n = rec["n"].as_u64().unwrap() as u32;
                let irr = find_irreducible_sparse(n).unwrap();
                let st = FieldStructure::new(n, &irr);
                let el = |w: u64| F2mElement::from_biguint(&num_bigint::BigUint::from(w), n);
                let basis: Vec<F2mElement> = rec["v_basis"].as_array().unwrap().iter().map(|w| el(w.as_u64().unwrap())).collect();
                let sys = build_decomposition_system(&basis, &el(rec["x_r"].as_u64().unwrap()), &F2mElement::one(n), m, &st)
                    .expect("system builds");
                assert_eq!(sys.n_vars as u64, rec["n_vars"].as_u64().unwrap(), "rebuilt shape");
                let eqs: Vec<Vec<u64>> = sys.equations.iter().map(|p| p.terms.iter().map(|t| t.mask).collect()).collect();
                let nv = sys.n_vars as u32;
                let mut d_f = None;
                let mut dims = serde_json::Map::new();
                for c in 2..=cmax {
                    let full = monomials_deglex(nv, c).len();
                    let (_, closed) = closure_dims(&eqs, nv, c, full);
                    dims.insert(c.to_string(), json!(closed));
                    if closed == full {
                        d_f = Some(c);
                        break;
                    }
                }
                println!("{}", json!({
                    "schema": "hky_ladder_d_f/v1", "m": m, "cell": rec["cell"], "n": n, "ell": rec["ell"], "draw": rec["draw"],
                    "n_vars": nv, "n_eqs": eqs.len(), "d_star": rec["outcome"], "ffd_recorded": rec["ffd"],
                    "d_f": d_f, "censored_above": if d_f.is_none() { Some(cmax) } else { None },
                    "closure_dims": dims, "secs": (t0.elapsed().as_secs_f64() * 10.0).round() / 10.0,
                }));
            }
        }
        _ => eprintln!("usage: hky_last_fall_degree replay RECORDS.jsonl [--modulus HEX] | cell --n N [--family random|subfield] [--target rand|sat] [--draws K] --seed S [--cmax C]"),
    }
}
