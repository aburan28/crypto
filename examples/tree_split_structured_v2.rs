//! Tree-split cut-out screen, round 2: the follow-ups that
//! `research/ic_tree_split_structured_v_20261011/` (round 1) left open.
//!
//! Protocol: `research/ic_tree_split_structured_v2_20261011/PREREGISTRATION.md`.
//! Round 1's harness, `examples/tree_split_structured.rs`, stays frozen with
//! its records; this file copies its field, curve, subspace and screen code
//! unchanged and adds:
//!
//! - `frob2`: a Frobenius-stable `V = ker g(σ)` with `g | t^n − 1`,
//!   `deg g = ℓ`, and `g ∤ t^d − 1` for every proper divisor `d` of `n`, so
//!   `V` lies in no proper subfield (round 1's `frob` excluded only the
//!   subfields themselves, and six of its cells sat inside `GF(2^7)`);
//! - `nbw`: a normal-basis window `V = span{θ, θ², …, θ^{2^{ℓ−1}}}` for a
//!   seeded normal element `θ`;
//! - the object `Z = {(e + f, e·f) : e, f ∈ T2} ⊂ GF(2)^{2n}` with round 1's
//!   predecessor's controls `ZR` (`Sym²` of a random set of `|T2|` field
//!   elements) and `R`;
//! - `poly`, `rand` and `subfield` exactly as in round 1.
//!
//! ```bash
//! cargo run --release --example tree_split_structured_v2 -- selftest
//! cargo run --release --example tree_split_structured_v2 -- cell --object T2 --vtype nbw --n 19 --a 0 --ell 6 --seed 20261012
//! cargo run --release --example tree_split_structured_v2 -- readout runs/registered/cells.jsonl
//! cargo run --release --example tree_split_structured_v2 -- lean runs/registered/cells.jsonl lean/TreeSplit2/Data.lean
//! ```
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashSet};
use std::io::{BufRead, Write};
use std::time::Instant;

// ------------------------------------------------------------------ RNG

/// `splitmix64`, chosen because it is four lines in any language.
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

const TAG_V: u64 = 1;
const TAG_CONTROL: u64 = 2;
const TAG_SAMPLES_CURVE: u64 = 3;
const TAG_SAMPLES_R: u64 = 4;
const TAG_ZR: u64 = 5;
const TAG_SAMPLES_ZR: u64 = 6;

fn stream(seed: u64, n: u32, a: u32, ell: u32, tag: u64) -> Rng {
    Rng(seed ^ ((a as u64) << 40) ^ ((n as u64) << 32) ^ ((ell as u64) << 24) ^ (tag << 56))
}

// ------------------------------------------------------------------ GF(2)[t] and GF(2^n)

fn deg(p: u64) -> i32 {
    63 - p.leading_zeros() as i32
}

fn poly_mod(mut a: u64, b: u64) -> u64 {
    let db = deg(b);
    while a != 0 && deg(a) >= db {
        a ^= b << (deg(a) - db);
    }
    a
}

fn poly_div(mut a: u64, b: u64) -> (u64, u64) {
    let db = deg(b);
    let mut q = 0;
    while a != 0 && deg(a) >= db {
        let s = deg(a) - db;
        q |= 1 << s;
        a ^= b << s;
    }
    (q, a)
}

fn poly_gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = poly_mod(a, b);
        a = b;
        b = r;
    }
    a
}

fn poly_mul(a: u64, b: u64) -> u64 {
    let mut r = 0;
    for i in 0..64 {
        if b >> i & 1 == 1 {
            r ^= a << i;
        }
    }
    r
}

#[derive(Clone, Copy)]
struct Field {
    n: u32,
    poly: u64,
}

impl Field {
    /// The same choice as `ic_tree_split_cutout_20260930/cutout.py`: the
    /// first irreducible trinomial `t^n + t^k + 1`, else the first pentanomial
    /// in lexicographic `(a, b, c)` order.
    fn new(n: u32) -> Field {
        for k in 1..n {
            let p = (1u64 << n) | (1 << k) | 1;
            if is_irreducible(p, n) {
                return Field { n, poly: p };
            }
        }
        for a in 1..n {
            for b in a + 1..n {
                for c in b + 1..n {
                    let p = (1u64 << n) | (1 << a) | (1 << b) | (1 << c) | 1;
                    if is_irreducible(p, n) {
                        return Field { n, poly: p };
                    }
                }
            }
        }
        panic!("no irreducible of degree {n}");
    }

    fn mul(&self, mut a: u64, mut b: u64) -> u64 {
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

    fn sqr(&self, a: u64) -> u64 {
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

    fn trace(&self, a: u64) -> u64 {
        let (mut t, mut x) = (0, a);
        for _ in 0..self.n {
            t ^= x;
            x = self.sqr(x);
        }
        assert!(t <= 1, "trace outside GF(2)");
        t
    }

    /// `z` with `z² + z = c` when `Tr(c) = 0`; `n` odd.
    fn half_trace(&self, c: u64) -> u64 {
        let (mut z, mut x) = (0, c);
        for _ in (0..self.n).step_by(2) {
            z ^= x;
            x = self.sqr(self.sqr(x));
        }
        z
    }

    fn mask(&self) -> u64 {
        (1u64 << self.n) - 1
    }
}

fn is_irreducible(p: u64, n: u32) -> bool {
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

// ------------------------------------------------------------------ Koblitz curve

type Pt = Option<(u64, u64)>;

struct Curve {
    f: Field,
    a: u64,
}

impl Curve {
    fn points_with_x(&self, x: u64) -> Vec<(u64, u64)> {
        let f = &self.f;
        if x == 0 {
            return vec![];
        }
        let c = x ^ self.a ^ f.inv(f.sqr(x));
        if f.trace(c) != 0 {
            return vec![];
        }
        let y = f.mul(x, f.half_trace(c));
        vec![(x, y), (x, y ^ x)]
    }

    fn add(&self, p: Pt, q: Pt) -> Pt {
        let f = &self.f;
        let ((x1, y1), (x2, y2)) = match (p, q) {
            (None, _) => return q,
            (_, None) => return p,
            (Some(p), Some(q)) => (p, q),
        };
        if x1 == x2 {
            if y1 ^ y2 == x1 || x1 == 0 {
                return None;
            }
            let lam = x1 ^ f.mul(y1, f.inv(x1));
            let x3 = f.sqr(lam) ^ lam ^ self.a;
            let y3 = f.sqr(x1) ^ f.mul(lam ^ 1, x3);
            return Some((x3, y3));
        }
        let lam = f.mul(y1 ^ y2, f.inv(x1 ^ x2));
        let x3 = f.sqr(lam) ^ lam ^ x1 ^ x2 ^ self.a;
        let y3 = f.mul(lam, x1 ^ x3) ^ x3 ^ y1;
        Some((x3, y3))
    }

    fn on_curve(&self, (x, y): (u64, u64)) -> bool {
        let f = &self.f;
        f.sqr(y) ^ f.mul(x, y) == f.mul(f.sqr(x), x) ^ f.mul(self.a, f.sqr(x)) ^ 1
    }
}

// ------------------------------------------------------------------ factor-base subspaces

/// Irreducible factors of `t^n − 1` over GF(2), sorted by `(degree, bits)`.
fn factor_xn1(n: u32) -> Vec<u64> {
    let mut rem = (1u64 << n) | 1;
    let mut out = vec![];
    for d in 1..=n {
        if deg(rem) < d as i32 {
            break;
        }
        for p in (1u64 << d)..(1u64 << (d + 1)) {
            while deg(rem) >= d as i32 {
                let (q, r) = poly_div(rem, p);
                if r != 0 {
                    break;
                }
                out.push(p);
                rem = q;
            }
        }
    }
    assert_eq!(rem, 1, "t^{n} - 1 did not factor completely");
    out
}

/// Kernel of the GF(2)-linear map `x ↦ Σ_k g_k x^{2^k}` on `GF(2^n)`.
fn frobenius_kernel(f: &Field, g: u64) -> Vec<u64> {
    let n = f.n;
    // rows (image, source) for the basis z^i; eliminate on the image.
    let mut rows: Vec<(u64, u64)> = (0..n)
        .map(|i| {
            let mut x = 1u64 << i;
            let mut img = 0;
            for k in 0..=deg(g) {
                if g >> k & 1 == 1 {
                    img ^= x;
                }
                x = f.sqr(x);
            }
            (img, 1u64 << i)
        })
        .collect();
    let mut r = 0;
    for col in 0..n {
        let bit = 1u64 << col;
        if let Some(p) = (r..rows.len()).find(|&i| rows[i].0 & bit != 0) {
            rows.swap(r, p);
            let piv = rows[r];
            for (i, row) in rows.iter_mut().enumerate() {
                if i != r && row.0 & bit != 0 {
                    row.0 ^= piv.0;
                    row.1 ^= piv.1;
                }
            }
            r += 1;
        }
    }
    let basis: Vec<u64> = rows[r..].iter().map(|row| row.1).collect();
    span(&basis)
}

fn span(basis: &[u64]) -> Vec<u64> {
    let mut out = vec![0u64];
    for &b in basis {
        let ext: Vec<u64> = out.iter().map(|&s| s ^ b).collect();
        out.extend(ext);
    }
    out.sort_unstable();
    out.dedup();
    out
}

#[derive(Serialize, Deserialize, Clone, Debug)]
struct VInfo {
    /// For `frob`: the chosen `g` (bit `k` is the coefficient of `t^k`), the
    /// number of admissible `g`, and the index chosen.
    g: Option<u64>,
    candidates: Option<usize>,
    choice: Option<usize>,
}

fn build_v(vtype: &str, f: &Field, ell: u32, seed: u64, rng: &mut Rng) -> (Vec<u64>, VInfo) {
    let n = f.n;
    let none = VInfo {
        g: None,
        candidates: None,
        choice: None,
    };
    match vtype {
        "rand" => {
            let mut basis: Vec<u64> = vec![];
            let mut sp: HashSet<u64> = [0].into();
            while basis.len() < ell as usize {
                let v = rng.next() & f.mask();
                if !sp.contains(&v) {
                    let ext: Vec<u64> = sp.iter().map(|&s| s ^ v).collect();
                    sp.extend(ext);
                    basis.push(v);
                }
            }
            (span(&basis), none)
        }
        "poly" => {
            assert!(2 * ell - 1 <= n, "poly V needs 2ℓ − 1 ≤ n");
            ((0..1u64 << ell).collect(), none)
        }
        "subfield" => {
            assert!(n % ell == 0, "subfield V needs ℓ | n");
            let g = (1u64 << ell) | 1;
            (
                frobenius_kernel(f, g),
                VInfo {
                    g: Some(g),
                    candidates: Some(1),
                    choice: Some(0),
                },
            )
        }
        "frob" => {
            let fac = factor_xn1(n);
            let subfields: Vec<u64> = (1..n)
                .filter(|d| n % d == 0)
                .map(|d| (1u64 << d) | 1)
                .collect();
            let mut cands = vec![];
            for mask in 1u64..(1 << fac.len()) {
                let g = (0..fac.len())
                    .filter(|&i| mask >> i & 1 == 1)
                    .fold(1u64, |acc, i| poly_mul(acc, fac[i]));
                if deg(g) == ell as i32 && !subfields.contains(&g) {
                    cands.push(g);
                }
            }
            assert!(
                !cands.is_empty(),
                "no Frobenius-stable V of dim {ell} in GF(2^{n})"
            );
            // `seed mod #candidates`: consecutive seeds pick distinct `g`.
            let choice = (seed % cands.len() as u64) as usize;
            let g = cands[choice];
            (
                frobenius_kernel(f, g),
                VInfo {
                    g: Some(g),
                    candidates: Some(cands.len()),
                    choice: Some(choice),
                },
            )
        }
        "frob2" => {
            let fac = factor_xn1(n);
            let proper: Vec<u64> = (1..n)
                .filter(|d| n % d == 0)
                .map(|d| (1u64 << d) | 1)
                .collect();
            let mut cands = vec![];
            for mask in 1u64..(1 << fac.len()) {
                let g = (0..fac.len())
                    .filter(|&i| mask >> i & 1 == 1)
                    .fold(1u64, |acc, i| poly_mul(acc, fac[i]));
                if deg(g) == ell as i32 && proper.iter().all(|&t| poly_mod(t, g) != 0) {
                    cands.push(g);
                }
            }
            assert!(
                !cands.is_empty(),
                "no subfield-free Frobenius-stable V of dim {ell} in GF(2^{n})"
            );
            let choice = (seed % cands.len() as u64) as usize;
            let g = cands[choice];
            (
                frobenius_kernel(f, g),
                VInfo {
                    g: Some(g),
                    candidates: Some(cands.len()),
                    choice: Some(choice),
                },
            )
        }
        "nbw" => {
            // First seeded θ whose n conjugates are linearly independent.
            let theta = loop {
                let t = rng.next() & f.mask();
                let conj: Vec<u64> = (0..n)
                    .scan(t, |x, _| {
                        let y = *x;
                        *x = f.sqr(*x);
                        Some(y)
                    })
                    .collect();
                if gf2_rank(&conj) == n as usize {
                    break t;
                }
            };
            let basis: Vec<u64> = (0..ell)
                .scan(theta, |x, _| {
                    let y = *x;
                    *x = f.sqr(*x);
                    Some(y)
                })
                .collect();
            (
                span(&basis),
                VInfo {
                    g: Some(theta),
                    candidates: None,
                    choice: None,
                },
            )
        }
        _ => panic!("unknown vtype {vtype}"),
    }
}

/// GF(2) rank of a list of bit vectors.
fn gf2_rank(vs: &[u64]) -> usize {
    let mut basis: Vec<u64> = vec![];
    for &v in vs {
        let mut v = v;
        for &b in &basis {
            v = v.min(v ^ b);
        }
        if v != 0 {
            basis.push(v);
            basis.sort_unstable_by(|a, b| b.cmp(a));
        }
    }
    basis.len()
}

/// `{(e + f) << n | e·f : e, f ∈ t}`, unordered with `e = f` allowed.
fn sym2(f: &Field, t: &[u64]) -> Vec<u64> {
    let mut out = BTreeSet::new();
    for i in 0..t.len() {
        for j in i..t.len() {
            out.insert(((t[i] ^ t[j]) << f.n) | f.mul(t[i], t[j]));
        }
    }
    out.into_iter().collect()
}

fn random_bits(nbits: u32, size: usize, rng: &mut Rng) -> Vec<u64> {
    let mask = (1u64 << nbits) - 1;
    let mut s = BTreeSet::new();
    while s.len() < size {
        let x = rng.next() & mask;
        if x != 0 {
            s.insert(x);
        }
    }
    s.into_iter().collect()
}

// ------------------------------------------------------------------ the sets

fn tree_split_t2(curve: &Curve, v: &[u64]) -> (usize, Vec<u64>) {
    let base: Vec<(u64, u64)> = v.iter().flat_map(|&x| curve.points_with_x(x)).collect();
    assert!(base.iter().all(|&p| curve.on_curve(p)));
    let mut t2 = BTreeSet::new();
    for i in 0..base.len() {
        for j in i..base.len() {
            if let Some((x, _)) = curve.add(Some(base[i]), Some(base[j])) {
                t2.insert(x);
            }
        }
    }
    (base.len(), t2.into_iter().collect())
}

fn random_control(f: &Field, size: usize, rng: &mut Rng) -> Vec<u64> {
    let mut s = BTreeSet::new();
    while s.len() < size {
        let x = rng.next() & f.mask();
        if x != 0 {
            s.insert(x);
        }
    }
    s.into_iter().collect()
}

// ------------------------------------------------------------------ cut-out screen

/// Multilinear monomials of degree `≤ d` in `nvars` variables, as bitmasks,
/// by degree and then lexicographic order of their index sets.
fn monomials(nvars: u32, d: u32) -> Vec<u64> {
    fn rec(start: u32, nvars: u32, left: u32, acc: u64, out: &mut Vec<u64>) {
        if left == 0 {
            out.push(acc);
            return;
        }
        for i in start..nvars {
            rec(i + 1, nvars, left - 1, acc | 1 << i, out);
        }
    }
    let mut out = vec![];
    for k in 0..=d {
        rec(0, nvars, k, 0, &mut out);
    }
    out
}

fn eval_row(p: u64, monos: &[u64], w: usize) -> Vec<u64> {
    let mut row = vec![0u64; w];
    for (j, &m) in monos.iter().enumerate() {
        if p & m == m {
            row[j / 64] |= 1 << (j % 64);
        }
    }
    row
}

/// GF(2) row echelon form of `rows`: returns `(rows', pivots)`, where row
/// `r` of the first `pivots.len()` has its leading bit at `pivots[r]` and
/// every later row is zero there.
fn echelon(mut rows: Vec<Vec<u64>>, w: usize) -> (Vec<Vec<u64>>, Vec<usize>) {
    let s = rows.len();
    let mut pivots = vec![];
    let mut r = 0;
    for col in 0..w * 64 {
        if r == s {
            break;
        }
        let (cw, bit) = (col / 64, 1u64 << (col % 64));
        let Some(p) = (r..s).find(|&i| rows[i][cw] & bit != 0) else {
            continue;
        };
        rows.swap(r, p);
        let (head, tail) = rows.split_at_mut(r + 1);
        let pivot = &head[r];
        tail.par_iter_mut().for_each(|row| {
            if row[cw] & bit != 0 {
                for (x, &y) in row[cw..].iter_mut().zip(&pivot[cw..]) {
                    *x ^= y;
                }
            }
        });
        pivots.push(col);
        r += 1;
    }
    rows.truncate(r);
    (rows, pivots)
}

fn in_span(rows: &[Vec<u64>], pivots: &[usize], mut v: Vec<u64>) -> bool {
    for (row, &col) in rows.iter().zip(pivots) {
        let (cw, bit) = (col / 64, 1u64 << (col % 64));
        if v[cw] & bit != 0 {
            for (x, &y) in v[cw..].iter_mut().zip(&row[cw..]) {
                *x ^= y;
            }
        }
    }
    v.iter().all(|&x| x == 0)
}

#[derive(Serialize, Deserialize, Clone, Debug)]
struct DegreeRow {
    d: u32,
    m: usize,
    rank: Option<usize>,
    deficit: Option<usize>,
    false_zeros: Option<usize>,
    skipped: Option<String>,
    secs: Option<f64>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
struct Profile {
    size: usize,
    nvars: u32,
    /// `None` when censored: no `D ≤ d_max` within the column cap had zero
    /// false zeros.
    d_s: Option<u32>,
    censored: bool,
    rows: Vec<DegreeRow>,
}

fn profile(
    points: &[u64],
    nvars: u32,
    d_max: u32,
    k: usize,
    max_cols: usize,
    rng: &mut Rng,
) -> Profile {
    let set: HashSet<u64> = points.iter().copied().collect();
    let mask = (1u64 << nvars) - 1;
    let mut samples = Vec::with_capacity(k);
    while samples.len() < k {
        let p = rng.next() & mask;
        if !set.contains(&p) {
            samples.push(p);
        }
    }
    let mut rows = vec![];
    let mut d_s = None;
    for d in 1..=d_max {
        let monos = monomials(nvars, d);
        if monos.len() > max_cols {
            rows.push(DegreeRow {
                d,
                m: monos.len(),
                rank: None,
                deficit: None,
                false_zeros: None,
                skipped: Some(format!("{} monomials > max_cols {}", monos.len(), max_cols)),
                secs: None,
            });
            break;
        }
        let t0 = Instant::now();
        let w = monos.len().div_ceil(64);
        let mat: Vec<Vec<u64>> = points.par_iter().map(|&p| eval_row(p, &monos, w)).collect();
        let (ech, piv) = echelon(mat, w);
        let rank = piv.len();
        // Full column rank: no nonzero polynomial of degree ≤ d vanishes on
        // the set, so every sample is a common zero of the (empty) ideal part.
        let false_zeros = if rank == monos.len() {
            k
        } else {
            samples
                .par_iter()
                .filter(|&&p| in_span(&ech, &piv, eval_row(p, &monos, w)))
                .count()
        };
        rows.push(DegreeRow {
            d,
            m: monos.len(),
            rank: Some(rank),
            deficit: Some(points.len().min(monos.len()) - rank),
            false_zeros: Some(false_zeros),
            skipped: None,
            secs: Some((t0.elapsed().as_secs_f64() * 100.0).round() / 100.0),
        });
        if false_zeros == 0 {
            d_s = Some(d);
            break;
        }
    }
    Profile {
        size: points.len(),
        nvars,
        d_s,
        censored: d_s.is_none(),
        rows,
    }
}

// ------------------------------------------------------------------ cell

#[derive(Serialize, Deserialize, Clone, Debug)]
struct Cell {
    object: String,
    vtype: String,
    n: u32,
    a: u32,
    ell: u32,
    seed: u64,
    field_poly: u64,
    v_size: usize,
    v_info: VInfo,
    factor_base_points: usize,
    t2_size: usize,
    /// `|Z|` for `object = Z`.
    z_size: Option<usize>,
    /// For `subfield` only: every element of `T2` lies in `V`.
    t2_in_v: Option<bool>,
    curve: Profile,
    r: Profile,
    /// `object = Z` only: `Sym²` of a random set of `|T2|` field elements.
    zr: Option<Profile>,
    d_max: u32,
    samples: usize,
    max_cols: usize,
}

#[allow(clippy::too_many_arguments)]
fn run_cell(
    object: &str,
    vtype: &str,
    n: u32,
    a: u32,
    ell: u32,
    seed: u64,
    d_max: u32,
    k: usize,
    max_cols: usize,
) -> Cell {
    assert!(n % 2 == 1, "half-trace solver needs n odd");
    let f = Field::new(n);
    let curve = Curve { f, a: a as u64 };
    let (v, v_info) = build_v(vtype, &f, ell, seed, &mut stream(seed, n, a, ell, TAG_V));
    assert_eq!(v.len(), 1 << ell, "dim V != ℓ");
    if matches!(vtype, "frob" | "frob2" | "subfield") {
        let vs: HashSet<u64> = v.iter().copied().collect();
        assert!(
            v.iter().all(|&x| vs.contains(&f.sqr(x))),
            "V not Frobenius-stable"
        );
    }
    let (fb, t2) = tree_split_t2(&curve, &v);
    let t2_in_v = (vtype == "subfield").then(|| {
        let vs: HashSet<u64> = v.iter().copied().collect();
        t2.iter().all(|x| vs.contains(x))
    });
    let (set, nvars, zr_set) = match object {
        "T2" => (t2.clone(), n, None),
        "Z" => {
            let rt = random_control(&f, t2.len(), &mut stream(seed, n, a, ell, TAG_ZR));
            (sym2(&f, &t2), 2 * n, Some(sym2(&f, &rt)))
        }
        _ => panic!("unknown object {object}"),
    };
    let control = random_bits(nvars, set.len(), &mut stream(seed, n, a, ell, TAG_CONTROL));
    let cp = profile(
        &set,
        nvars,
        d_max,
        k,
        max_cols,
        &mut stream(seed, n, a, ell, TAG_SAMPLES_CURVE),
    );
    let rp = profile(
        &control,
        nvars,
        d_max,
        k,
        max_cols,
        &mut stream(seed, n, a, ell, TAG_SAMPLES_R),
    );
    let zp = zr_set.map(|z| {
        profile(
            &z,
            nvars,
            d_max,
            k,
            max_cols,
            &mut stream(seed, n, a, ell, TAG_SAMPLES_ZR),
        )
    });
    Cell {
        object: object.into(),
        vtype: vtype.into(),
        n,
        a,
        ell,
        seed,
        field_poly: f.poly,
        v_size: v.len(),
        v_info,
        factor_base_points: fb,
        t2_size: t2.len(),
        z_size: (object == "Z").then_some(set.len()),
        t2_in_v,
        curve: cp,
        r: rp,
        zr: zp,
        d_max,
        samples: k,
        max_cols,
    }
}

// ------------------------------------------------------------------ readout

/// `D_s(curve) − min over controls of D_s`, or `None` if any arm is censored.
fn gap(c: &Cell) -> Option<i64> {
    let dc = c.curve.d_s? as i64;
    let mut m = c.r.d_s? as i64;
    if let Some(z) = &c.zr {
        m = m.min(z.d_s? as i64);
    }
    Some(dc - m)
}

/// The registered decision rule (PREREGISTRATION.md §4), mirrored in
/// `lean/TreeSplit2/Rule.lean`, which is what the Lean theorem checks.
/// `slope_clause = false` for the `Z` groups.
fn verdict(cells: &[&Cell], slope_clause: bool) -> (String, Vec<String>) {
    let mut series: Vec<(u32, u32, u64)> = cells.iter().map(|c| (c.n, c.a, c.seed)).collect();
    series.sort();
    series.dedup();
    let mut notes = vec![];
    let (mut all_closed, mut all_alive_gap, mut slopes_ok, mut nonempty) =
        (true, true, 0usize, 0usize);
    for &(n, a, seed) in &series {
        let mut pts: Vec<(i64, i64, i64)> = vec![];
        for c in cells.iter().filter(|c| (c.n, c.a, c.seed) == (n, a, seed)) {
            match gap(c) {
                Some(g) => pts.push((c.ell as i64, c.curve.d_s.unwrap() as i64, g)),
                None => notes.push(format!(
                    "n={n} a={a} seed={seed} ell={}: censored, excluded",
                    c.ell
                )),
            }
        }
        pts.sort();
        if pts.is_empty() {
            continue;
        }
        nonempty += 1;
        if pts.iter().any(|&(_, _, g)| g < -1) {
            all_closed = false;
        }
        if pts.last().unwrap().2 > -2 {
            all_alive_gap = false;
        }
        let k = pts.len() as i64;
        let (sl, sd, sll, sld) = pts.iter().fold((0, 0, 0, 0), |(a, b, c, d), &(l, ds, _)| {
            (a + l, b + ds, c + l * l, d + l * ds)
        });
        let (num, den) = (k * sld - sl * sd, k * sll - sl * sl);
        slopes_ok += (den > 0 && 100 * num <= 15 * den) as usize;
        notes.push(format!(
            "n={n} a={a} seed={seed}: gaps {:?}, slope {:.3}",
            pts.iter().map(|p| p.2).collect::<Vec<_>>(),
            if den > 0 {
                num as f64 / den as f64
            } else {
                f64::NAN
            }
        ));
    }
    let slope_ok = !slope_clause || 2 * slopes_ok >= nonempty;
    let v = if nonempty > 0 && all_alive_gap && slope_ok {
        "ALIVE"
    } else if all_closed {
        "CLOSED"
    } else {
        "INCONCLUSIVE"
    };
    (v.into(), notes)
}

fn read_cells(path: &str) -> Vec<Cell> {
    let f = std::fs::File::open(path).expect("open cells");
    std::io::BufReader::new(f)
        .lines()
        .map(|l| serde_json::from_str(&l.unwrap()).expect("parse cell"))
        .collect()
}

/// The registered groups: `(name, object, vtype, n filter, slope clause)`.
type Group = (
    &'static str,
    &'static str,
    &'static str,
    fn(u32) -> bool,
    bool,
);
const GROUPS: [Group; 6] = [
    ("T2/frob2 (n = 21)", "T2", "frob2", |n| n == 21, true),
    ("T2/nbw", "T2", "nbw", |_| true, true),
    (
        "T2/poly (n = 23, 25)",
        "T2",
        "poly",
        |n| n == 23 || n == 25,
        true,
    ),
    ("Z/poly", "Z", "poly", |_| true, false),
    ("Z/nbw", "Z", "nbw", |_| true, false),
    ("Z/frob2", "Z", "frob2", |_| true, false),
];

fn readout(path: &str) {
    let cells = read_cells(path);
    println!("| object | vtype | n | a | seed | ℓ | |F| | |T2| | |Z| | D_s(curve) | D_s(R) | D_s(ZR) | gap |");
    println!("|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|");
    let ds = |p: &Profile| p.d_s.map_or("cens".to_string(), |d| d.to_string());
    for c in &cells {
        println!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
            c.object,
            c.vtype,
            c.n,
            c.a,
            c.seed,
            c.ell,
            c.factor_base_points,
            c.t2_size,
            c.z_size.map_or("—".into(), |z| z.to_string()),
            ds(&c.curve),
            ds(&c.r),
            c.zr.as_ref().map_or("—".into(), ds),
            gap(c).map_or("—".into(), |g| g.to_string())
        );
    }
    let gate: Vec<&Cell> = cells.iter().filter(|c| c.vtype == "subfield").collect();
    let gate_ok = !gate.is_empty()
        && gate.iter().all(|c| {
            c.t2_in_v == Some(true)
                && c.curve.rows.first().and_then(|r| r.deficit).unwrap_or(0)
                    >= (c.n - c.ell) as usize
        });
    println!(
        "\npositive control (subfield, deficit(1) ≥ n − ℓ and T2 ⊆ V): {}",
        if gate_ok { "PASS" } else { "FAIL" }
    );
    for (name, obj, vt, keep, slope) in GROUPS {
        let sel: Vec<&Cell> = cells
            .iter()
            .filter(|c| c.object == obj && c.vtype == vt && keep(c.n))
            .collect();
        if sel.is_empty() {
            continue;
        }
        let (v, notes) = verdict(&sel, slope);
        println!("\n{name}: {v}");
        for n in notes {
            println!("  {n}");
        }
    }
}

// ------------------------------------------------------------------ Lean data

fn lean_opt(x: Option<u32>) -> String {
    x.map_or("none".into(), |d| format!("some {d}"))
}

fn emit_lean(path: &str, out: &str) {
    let cells = read_cells(path);
    let mut s = String::new();
    s.push_str("import TreeSplit2.Rule\n\n");
    s.push_str(
        "/-! Generated by `examples/tree_split_structured_v2.rs lean` from the registered\n",
    );
    s.push_str(&format!(
        "cell records `{}`. Do not edit: regenerate. -/\n",
        path.replace("/-", "/ -")
    ));
    s.push_str("\nnamespace TreeSplit2\n\n/-- The registered cells, one row each. -/\n");
    s.push_str("def cells : List Row := [\n");
    for (i, c) in cells.iter().enumerate() {
        let def1 = c.curve.rows.first().and_then(|r| r.deficit).unwrap_or(0);
        s.push_str(&format!(
            "  {{ obj := .{}, vtype := .{}, n := {}, a := {}, seed := {}, ell := {}, size := {}, dsCurve := {}, dsR := {}, dsZR := {}, deficit1 := {}, t2InV := {} }}{}\n",
            c.object.to_lowercase(),
            c.vtype,
            c.n,
            c.a,
            c.seed,
            c.ell,
            c.z_size.unwrap_or(c.t2_size),
            lean_opt(c.curve.d_s),
            lean_opt(c.r.d_s),
            c.zr.as_ref().map_or("none".into(), |z| format!("some ({})", lean_opt(z.d_s))),
            def1,
            c.t2_in_v.unwrap_or(false),
            if i + 1 < cells.len() { "," } else { "" }
        ));
    }
    s.push_str("]\n\nend TreeSplit2\n");
    std::fs::write(out, s).expect("write lean");
}

// ------------------------------------------------------------------ self-test

fn selftest() {
    let mut rng = Rng(20261011);
    for (a, n) in [(0u32, 13u32), (1, 17), (0, 7), (1, 21)] {
        let f = Field::new(n);
        for _ in 0..50 {
            let x = (rng.next() & f.mask()).max(1);
            assert_eq!(f.mul(x, f.inv(x)), 1, "inverse");
        }
        let curve = Curve { f, a: a as u64 };
        let mut pts = vec![];
        while pts.len() < 30 {
            let x = rng.next() & f.mask();
            pts.extend(curve.points_with_x(x));
        }
        for i in 0..pts.len() {
            let (p, q, r) = (
                Some(pts[i]),
                Some(pts[(i + 1) % 30]),
                Some(pts[(i + 7) % 30]),
            );
            assert_eq!(curve.add(p, q), curve.add(q, p), "commutativity");
            assert_eq!(
                curve.add(curve.add(p, q), r),
                curve.add(p, curve.add(q, r)),
                "associativity"
            );
            assert!(curve.add(p, q).map_or(true, |s| curve.on_curve(s)));
        }
    }
    // #K_0(GF(2^7)) = 116, and 116·P = O.
    let f7 = Field::new(7);
    let k0 = Curve { f: f7, a: 0 };
    let count = 1
        + (0..128u64)
            .map(|x| if x == 0 { 1 } else { k0.points_with_x(x).len() })
            .sum::<usize>();
    assert_eq!(count, 116, "#K_0(GF(2^7))");
    for p in (1..128u64).flat_map(|x| k0.points_with_x(x)).take(5) {
        let mut acc: Pt = None;
        for _ in 0..116 {
            acc = k0.add(acc, Some(p));
        }
        assert_eq!(acc, None, "116·P");
    }
    // Frobenius-stable subspaces and subfields.
    let f21 = Field::new(21);
    assert_eq!(
        factor_xn1(21).iter().map(|&p| deg(p)).collect::<Vec<_>>(),
        vec![1, 2, 3, 3, 6, 6]
    );
    assert_eq!(
        factor_xn1(17).iter().map(|&p| deg(p)).collect::<Vec<_>>(),
        vec![1, 8, 8]
    );
    let sub7 = frobenius_kernel(&f21, (1 << 7) | 1);
    assert_eq!(sub7.len(), 128);
    assert!(
        sub7.iter().all(|&x| (0..7).fold(x, |y, _| f21.sqr(y)) == x),
        "GF(2^7) = ker(σ^7 − 1)"
    );
    for ell in 4..=8 {
        let (v, info) = build_v("frob", &f21, ell, ell as u64, &mut Rng(ell as u64));
        assert_eq!(v.len(), 1 << ell);
        assert!(v.iter().all(|x| v.binary_search(&f21.sqr(*x)).is_ok()));
        println!(
            "frob n=21 ℓ={ell}: {} candidate g",
            info.candidates.unwrap()
        );
    }
    // Linear-algebra screen checks (as in the 2026-09-30 study §5).
    let mut r = Rng(7);
    let basis: Vec<u64> = (0..8).map(|_| r.next() & 0xFFFFF).collect();
    let shift = r.next() & 0xFFFFF;
    let aff: Vec<u64> = span(&basis).into_iter().map(|x| x ^ shift).collect();
    assert_eq!(aff.len(), 256);
    assert_eq!(
        profile(&aff, 20, 6, 2000, 90000, &mut Rng(1)).d_s,
        Some(1),
        "affine subspace"
    );
    let quad: Vec<u64> = (0..1u64 << 14)
        .filter(|x| (x & 1) & (x >> 1 & 1) == (x >> 2 & 1))
        .collect();
    let qp = profile(&quad, 14, 6, 2000, 90000, &mut Rng(2));
    assert_eq!(qp.d_s, Some(2), "quadric D_s");
    assert_eq!(qp.rows[1].deficit, Some(3), "quadric deficit at D = 2");
    let rnd = random_control(&Field { n: 20, poly: 0 }, 1000, &mut Rng(3));
    assert_eq!(
        profile(&rnd, 20, 6, 2000, 90000, &mut Rng(4)).d_s,
        Some(3),
        "random 1000 in GF(2)^20"
    );
    // round 2: subfield-free Frobenius kernels, normal-basis windows, Sym².
    let counts: Vec<usize> = (4..=8)
        .map(|ell| {
            let fac = factor_xn1(21);
            (1u64..(1 << fac.len()))
                .map(|m| {
                    (0..fac.len())
                        .filter(|&i| m >> i & 1 == 1)
                        .fold(1u64, |acc, i| poly_mul(acc, fac[i]))
                })
                .filter(|&g| {
                    deg(g) == ell
                        && [(1u64 << 3) | 1, (1 << 7) | 1]
                            .iter()
                            .all(|&t| poly_mod(t, g) != 0)
                })
                .count()
        })
        .collect();
    assert_eq!(counts, vec![0, 2, 4, 2, 3], "frob2 candidates at n = 21");
    for ell in 5..=8 {
        let (v, _) = build_v("frob2", &f21, ell, ell as u64, &mut Rng(1));
        assert_eq!(v.len(), 1 << ell);
        assert!(v.iter().all(|x| v.binary_search(&f21.sqr(*x)).is_ok()));
        // in no proper subfield: some element is outside GF(2^3) and GF(2^7)
        for d in [3u32, 7] {
            assert!(v.iter().any(|&x| (0..d).fold(x, |y, _| f21.sqr(y)) != x));
        }
    }
    for n in [17u32, 19, 21, 23] {
        let f = Field::new(n);
        for ell in 4..=8 {
            let (v, info) = build_v("nbw", &f, ell, 0, &mut Rng(ell as u64));
            assert_eq!(v.len(), 1 << ell, "nbw dim");
            let th = info.g.unwrap();
            let conj: Vec<u64> = (0..n)
                .scan(th, |x, _| {
                    let y = *x;
                    *x = f.sqr(*x);
                    Some(y)
                })
                .collect();
            assert_eq!(gf2_rank(&conj), n as usize, "θ normal");
        }
    }
    let f13 = Field::new(13);
    let t: Vec<u64> = (1..40u64).collect();
    assert_eq!(
        sym2(&f13, &t).len(),
        39 * 40 / 2,
        "Sym² injective on unordered pairs"
    );
    println!("selftest: ok");
}

// ------------------------------------------------------------------ main

fn arg<T: std::str::FromStr>(args: &[String], name: &str, default: Option<T>) -> T {
    match args.iter().position(|a| a == name) {
        Some(i) => args[i + 1].parse().ok().expect("bad argument"),
        None => default.unwrap_or_else(|| panic!("missing {name}")),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("selftest") => selftest(),
        Some("cell") => {
            let vtype: String = arg(&args, "--vtype", None);
            let object: String = arg(&args, "--object", Some("T2".to_string()));
            let cell = run_cell(
                &object,
                &vtype,
                arg(&args, "--n", None),
                arg(&args, "--a", None),
                arg(&args, "--ell", None),
                arg(&args, "--seed", None),
                arg(&args, "--d-max", Some(6)),
                arg(&args, "--samples", Some(2000)),
                arg(&args, "--max-cols", Some(90000)),
            );
            let mut out = std::io::stdout();
            writeln!(out, "{}", serde_json::to_string(&cell).unwrap()).unwrap();
        }
        Some("readout") => readout(&args[2]),
        Some("lean") => emit_lean(&args[2], &args[3]),
        _ => eprintln!(
            "usage: tree_split_structured selftest | cell … | readout CELLS | lean CELLS OUT"
        ),
    }
}
