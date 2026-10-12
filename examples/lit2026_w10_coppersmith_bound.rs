//! lit2026 W10: can a Coppersmith lattice be the prime-field `S₅` small-root
//! oracle that index calculus at arity `m = 4` needs (`δ ≈ 1/5`)?
//!
//! Everything here is native Rust; Python only plots the recorded JSON.
//!
//! **Problem.** Fix a prime `p`, a curve `y² = x³ + ax + b` over `F_p` and a
//! target x-coordinate `z`. The arity-`m` decomposition oracle must find
//! `(x₁, …, x_m)` with `|x_i| < X = p^δ` and
//! `f(x₁, …, x_m) := S_{m+1}(x₁, …, x_m, z) ≡ 0 (mod p)`. For `S₃, S₄, S₅`
//! (`m = 2, 3, 4`) this example computes the support of `f` exactly (multivariate
//! interpolation of Semaev's resultant recurrence, checked against planted
//! decompositions), its Newton polytope, and the symmetric rewriting in
//! elementary symmetric functions `e₁, …, e_m`.
//!
//! **Strategy family analysed** ("full-rank triangular single-polynomial
//! Coppersmith", which contains Jochemsz–May basic and extended strategies and
//! the Meers–Nowakowski/cuso shift sets for a single polynomial and a known
//! modulus): choose a vertex `a` of `Newt(f)` as leading monomial, a finite
//! monomial set `M`, and for each `β ∈ M` the shift `x^{β−k a} f^k p^{ℓ−k}`
//! with `k = k(β) ∈ [0, ℓ]` such that its whole support lies in `M`. As `p → ∞`
//! at fixed `ℓ`, Howgrave-Graham succeeds iff `Σ_{β∈M} k(β) > Σ_{β∈M} λ·β`,
//! `λ_i = log_p X_i`.
//!
//! **Exact optimiser (the "LP").** For fixed `(a, ℓ)` and a bounding box, the
//! choice of `M` and `k(·)` is a maximum-weight-closure problem (its LP is
//! totally unimodular, so the LP optimum is integral); we solve it exactly by
//! min-cut (Dinic) inside Dinkelbach iteration for the ratio
//! `δ_ℓ = max (Σk − Σe·β)/Σc·β`. That is the exact optimum of the family at
//! multiplicity `ℓ`, over *all* monomial sets in the box — not only boxes,
//! polytopes or JM-shaped sets.
//!
//! **Upper bound for every `ℓ`.** For any vertex `c` of `Newt(f)` with
//! `(k−j)a + jc ∈ supp(f^k)` (always true for `c = a`; true for generic
//! coefficients otherwise) the family can only succeed if `λ·c < 1` (charging
//! argument in the study README; checked numerically on every optimum here).
//! For `S₅` (box `[0,8]^4`, corner `(8,8,8,8)`) this is `δ < 1/32`.
//!
//! ```bash
//! cargo run --release --example lit2026_w10_coppersmith_bound -- \
//!     [--part known|support|families|closure|lattice|all] [--json out.json] [--max-ell N]
//! ```

use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

use crypto_lib::cryptanalysis::lattice::{lll_reduce, lll_reduce_hp};
use crypto_lib::cryptanalysis::prime_fast::{find_a3_curve, FastCurve};
use num_bigint::BigInt;
use num_traits::{ToPrimitive, Zero};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

// ───────────────────────────── F_p (p < 2^62) ─────────────────────────────

#[derive(Clone, Copy, Debug)]
struct Fp {
    p: u64,
}

impl Fp {
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }
    fn neg(&self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }
    fn mul(&self, a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % self.p as u128) as u64
    }
    fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1u64 % self.p;
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
        assert!(a != 0, "inverse of zero");
        self.pow(a, self.p - 2)
    }
    fn c(&self, v: i64) -> u64 {
        let p = self.p as i128;
        (((v as i128 % p) + p) % p) as u64
    }
}

// ─────────────────────────── Semaev polynomials ───────────────────────────

#[derive(Clone, Copy, Debug)]
struct Weier {
    f: Fp,
    a: u64,
    b: u64,
}

impl Weier {
    /// Coefficients `[c0, c1, c2]` of `S₃(u, v, Y)` in `Y`.
    fn s3_y(&self, u: u64, v: u64) -> [u64; 3] {
        let f = self.f;
        let (a, b) = (self.a, self.b);
        let uv = f.mul(u, v);
        let upv = f.add(u, v);
        let d = f.sub(u, v);
        let c2 = f.mul(d, d);
        let inner = f.add(f.mul(upv, f.add(uv, a)), f.mul(f.c(2), b));
        let c1 = f.neg(f.mul(f.c(2), inner));
        let t = f.sub(uv, a);
        let c0 = f.sub(f.mul(t, t), f.mul(f.mul(f.c(4), b), upv));
        [c0, c1, c2]
    }
    fn s3(&self, x1: u64, x2: u64, x3: u64) -> u64 {
        let c = self.s3_y(x1, x2);
        let f = self.f;
        f.add(f.add(f.mul(c[2], f.mul(x3, x3)), f.mul(c[1], x3)), c[0])
    }
    /// `S₄(x₁,x₂,x₃,x₄) = Res_Y(S₃(x₁,x₂,Y), S₃(x₃,x₄,Y))`.
    fn s4(&self, x1: u64, x2: u64, x3: u64, x4: u64) -> u64 {
        res22(self.f, &self.s3_y(x1, x2), &self.s3_y(x3, x4))
    }
    /// Coefficients of `S₄(x₁,x₂,x₃,Y)` in `Y` (degree ≤ 4).
    fn s4_y(&self, x1: u64, x2: u64, x3: u64) -> Vec<u64> {
        let vals: Vec<u64> = (0..5u64).map(|t| self.s4(x1, x2, x3, t)).collect();
        let vinv = inv_vandermonde(self.f, 5);
        apply_mat(self.f, &vinv, &vals)
    }
    /// `S₅(x₁..x₅) = Res_Y(S₄(x₁,x₂,x₃,Y), S₃(x₄,x₅,Y))` (6×6 Sylvester).
    fn s5(&self, x: &[u64; 5]) -> u64 {
        let a = self.s4_y(x[0], x[1], x[2]);
        let b = self.s3_y(x[3], x[4]);
        sylvester(self.f, &a, &b)
    }
    fn s_n(&self, n: usize, x: &[u64]) -> u64 {
        match n {
            3 => self.s3(x[0], x[1], x[2]),
            4 => self.s4(x[0], x[1], x[2], x[3]),
            5 => self.s5(&[x[0], x[1], x[2], x[3], x[4]]),
            _ => unreachable!(),
        }
    }
}

fn res22(f: Fp, a: &[u64; 3], b: &[u64; 3]) -> u64 {
    // Res(a2Y²+a1Y+a0, b2Y²+b1Y+b0) = (a2b0−a0b2)² − (a2b1−a1b2)(a1b0−a0b1)
    let t1 = f.sub(f.mul(a[2], b[0]), f.mul(a[0], b[2]));
    let t2 = f.sub(f.mul(a[2], b[1]), f.mul(a[1], b[2]));
    let t3 = f.sub(f.mul(a[1], b[0]), f.mul(a[0], b[1]));
    f.sub(f.mul(t1, t1), f.mul(t2, t3))
}

/// Sylvester resultant of `A` (coeffs low→high, formal degree len−1) and `B`.
fn sylvester(f: Fp, a: &[u64], b: &[u64]) -> u64 {
    let da = a.len() - 1;
    let db = b.len() - 1;
    let n = da + db;
    let mut m = vec![vec![0u64; n]; n];
    for r in 0..db {
        for (i, &c) in a.iter().rev().enumerate() {
            m[r][r + i] = c;
        }
    }
    for r in 0..da {
        for (i, &c) in b.iter().rev().enumerate() {
            m[db + r][r + i] = c;
        }
    }
    det(f, m)
}

fn det(f: Fp, mut m: Vec<Vec<u64>>) -> u64 {
    let n = m.len();
    let mut d = 1u64;
    for col in 0..n {
        let Some(piv) = (col..n).find(|&r| m[r][col] != 0) else {
            return 0;
        };
        if piv != col {
            m.swap(piv, col);
            d = f.neg(d);
        }
        d = f.mul(d, m[col][col]);
        let inv = f.inv(m[col][col]);
        for r in col + 1..n {
            if m[r][col] == 0 {
                continue;
            }
            let fac = f.mul(m[r][col], inv);
            for c in col..n {
                let v = f.mul(fac, m[col][c]);
                m[r][c] = f.sub(m[r][c], v);
            }
        }
    }
    d
}

/// Inverse of the Vandermonde matrix on nodes `0, 1, …, n−1`.
fn inv_vandermonde(f: Fp, n: usize) -> Vec<Vec<u64>> {
    let mut a: Vec<Vec<u64>> = (0..n)
        .map(|k| {
            let mut row: Vec<u64> = (0..n).map(|e| f.pow(k as u64 % f.p, e as u64)).collect();
            row.extend((0..n).map(|j| u64::from(j == k)));
            row
        })
        .collect();
    for col in 0..n {
        let piv = (col..n).find(|&r| a[r][col] != 0).expect("vandermonde singular");
        a.swap(piv, col);
        let inv = f.inv(a[col][col]);
        for c in 0..2 * n {
            a[col][c] = f.mul(a[col][c], inv);
        }
        for r in 0..n {
            if r != col && a[r][col] != 0 {
                let fac = a[r][col];
                for c in 0..2 * n {
                    let v = f.mul(fac, a[col][c]);
                    a[r][c] = f.sub(a[r][c], v);
                }
            }
        }
    }
    a.into_iter().map(|row| row[n..].to_vec()).collect()
}

fn apply_mat(f: Fp, m: &[Vec<u64>], v: &[u64]) -> Vec<u64> {
    m.iter()
        .map(|row| row.iter().zip(v).fold(0u64, |acc, (&x, &y)| f.add(acc, f.mul(x, y))))
        .collect()
}

/// Dense polynomial in `m` variables, degree ≤ `deg` in each; index
/// `Σ e_i (deg+1)^i`.
#[derive(Clone, Debug)]
struct Dense {
    m: usize,
    deg: usize,
    c: Vec<u64>,
}

impl Dense {
    fn exps(&self, idx: usize) -> Vec<u8> {
        let mut e = vec![0u8; self.m];
        let mut r = idx;
        for v in e.iter_mut() {
            *v = (r % (self.deg + 1)) as u8;
            r /= self.deg + 1;
        }
        e
    }
    fn eval(&self, f: Fp, x: &[u64]) -> u64 {
        let mut acc = 0u64;
        for (idx, &c) in self.c.iter().enumerate() {
            if c == 0 {
                continue;
            }
            let e = self.exps(idx);
            let mut t = c;
            for i in 0..self.m {
                t = f.mul(t, f.pow(x[i], e[i] as u64));
            }
            acc = f.add(acc, t);
        }
        acc
    }
    fn support(&self) -> Vec<Vec<u8>> {
        (0..self.c.len()).filter(|&i| self.c[i] != 0).map(|i| self.exps(i)).collect()
    }
    fn to_sparse(&self) -> BTreeMap<Vec<u8>, u64> {
        (0..self.c.len())
            .filter(|&i| self.c[i] != 0)
            .map(|i| (self.exps(i), self.c[i]))
            .collect()
    }
}

/// Tensor-product interpolation of `g` (degree ≤ `deg` per variable) on the
/// grid `{0..deg}^m`.
fn interpolate(f: Fp, m: usize, deg: usize, g: &dyn Fn(&[u64]) -> u64) -> Dense {
    let n = deg + 1;
    let total = n.pow(m as u32);
    let mut vals = vec![0u64; total];
    let mut pt = vec![0u64; m];
    for (idx, v) in vals.iter_mut().enumerate() {
        let mut r = idx;
        for x in pt.iter_mut() {
            *x = (r % n) as u64;
            r /= n;
        }
        *v = g(&pt);
    }
    let vinv = inv_vandermonde(f, n);
    for axis in 0..m {
        let stride = n.pow(axis as u32);
        let mut out = vals.clone();
        for idx in 0..total {
            if (idx / stride) % n != 0 {
                continue;
            }
            let line: Vec<u64> = (0..n).map(|k| vals[idx + k * stride]).collect();
            let co = apply_mat(f, &vinv, &line);
            for k in 0..n {
                out[idx + k * stride] = co[k];
            }
        }
        vals = out;
    }
    Dense { m, deg, c: vals }
}

/// `f(x₁..x_m) = S_{m+1}(x₁..x_m, z)` as a dense polynomial.
fn s_with_target(w: &Weier, m: usize, z: u64) -> Dense {
    let deg = 1usize << (m - 1);
    interpolate(w.f, m, deg, &|x: &[u64]| {
        let mut full = x.to_vec();
        full.push(z);
        w.s_n(m + 1, &full)
    })
}

// ──────────────── symmetric rewriting in e₁..e_m (exact, mod p) ────────────────

type Sparse = BTreeMap<Vec<u8>, u64>;

fn sp_mul(f: Fp, a: &Sparse, b: &Sparse) -> Sparse {
    let mut out: Sparse = BTreeMap::new();
    for (ea, &ca) in a {
        for (eb, &cb) in b {
            let e: Vec<u8> = ea.iter().zip(eb).map(|(x, y)| x + y).collect();
            let v = out.entry(e).or_insert(0);
            *v = f.add(*v, f.mul(ca, cb));
        }
    }
    out.retain(|_, c| *c != 0);
    out
}

fn elementary(m: usize, i: usize) -> Sparse {
    let mut out = BTreeMap::new();
    for mask in 0u32..(1 << m) {
        if mask.count_ones() as usize == i {
            let e: Vec<u8> = (0..m).map(|j| ((mask >> j) & 1) as u8).collect();
            out.insert(e, 1u64);
        }
    }
    out
}

/// Write a symmetric polynomial as `G(e₁, …, e_m)`; returns `G` (exponents of
/// `e_i`) or `Err` if the input is not symmetric. The remainder reaching zero
/// is the exact identity check.
fn symmetric_reduce(f: Fp, m: usize, poly: &Sparse) -> Result<Sparse, String> {
    let es: Vec<Sparse> = (1..=m).map(|i| elementary(m, i)).collect();
    let mut memo: HashMap<Vec<u8>, Sparse> = HashMap::new();
    let one: Sparse = BTreeMap::from([(vec![0u8; m], 1u64)]);
    memo.insert(vec![0u8; m], one);
    fn epow(
        f: Fp,
        g: &[u8],
        es: &[Sparse],
        memo: &mut HashMap<Vec<u8>, Sparse>,
    ) -> Sparse {
        if let Some(v) = memo.get(g) {
            return v.clone();
        }
        let i = g.iter().position(|&x| x > 0).unwrap();
        let mut h = g.to_vec();
        h[i] -= 1;
        let prev = epow(f, &h, es, memo);
        let v = sp_mul(f, &prev, &es[i]);
        memo.insert(g.to_vec(), v.clone());
        v
    }
    let mut rem = poly.clone();
    let mut out: Sparse = BTreeMap::new();
    let mut guard = 0usize;
    while let Some((lead, &coef)) = rem.iter().next_back() {
        guard += 1;
        if guard > 1_000_000 {
            return Err("no termination".into());
        }
        let lead = lead.clone();
        if lead.windows(2).any(|w| w[0] < w[1]) {
            return Err(format!("not symmetric: leading monomial {lead:?}"));
        }
        let g: Vec<u8> = (0..m)
            .map(|i| lead[i] - if i + 1 < m { lead[i + 1] } else { 0 })
            .collect();
        let ex = epow(f, &g, &es, &mut memo);
        for (e, &c) in &ex {
            let v = rem.entry(e.clone()).or_insert(0);
            *v = f.sub(*v, f.mul(coef, c));
        }
        rem.retain(|_, c| *c != 0);
        out.insert(g, coef);
    }
    Ok(out)
}

// ──────────────────────── max-weight closure (Dinic) ────────────────────────

struct Dinic {
    n: usize,
    head: Vec<u32>,
    to: Vec<u32>,
    cap: Vec<i64>,
    level: Vec<i32>,
    it: Vec<u32>,
    start: Vec<u32>,
}

const INF: i64 = i64::MAX / 4;

impl Dinic {
    /// Build from an edge list `(u, v, cap)`; reverse edges are added.
    fn new(n: usize, edges: &[(u32, u32, i64)]) -> Self {
        let mut deg = vec![0u32; n + 1];
        for &(u, v, _) in edges {
            deg[u as usize] += 1;
            deg[v as usize] += 1;
        }
        let mut start = vec![0u32; n + 1];
        for i in 0..n {
            start[i + 1] = start[i] + deg[i];
        }
        let total = start[n] as usize;
        let mut to = vec![0u32; total];
        let mut cap = vec![0i64; total];
        let mut rev = vec![0u32; total];
        let mut fill = start.clone();
        for &(u, v, c) in edges {
            let iu = fill[u as usize];
            fill[u as usize] += 1;
            let iv = fill[v as usize];
            fill[v as usize] += 1;
            to[iu as usize] = v;
            cap[iu as usize] = c;
            rev[iu as usize] = iv;
            to[iv as usize] = u;
            cap[iv as usize] = 0;
            rev[iv as usize] = iu;
        }
        // store rev in head (reuse a vector name to keep the struct small)
        Dinic { n, head: rev, to, cap, level: vec![0; n], it: vec![0; n], start }
    }

    fn bfs(&mut self, s: usize, t: usize) -> bool {
        self.level.iter_mut().for_each(|l| *l = -1);
        let mut q = std::collections::VecDeque::new();
        self.level[s] = 0;
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            for e in self.start[u]..self.start[u + 1] {
                let v = self.to[e as usize] as usize;
                if self.cap[e as usize] > 0 && self.level[v] < 0 {
                    self.level[v] = self.level[u] + 1;
                    q.push_back(v);
                }
            }
        }
        self.level[t] >= 0
    }

    fn dfs(&mut self, s: usize, t: usize) -> i64 {
        // iterative augmenting-path search on the level graph
        let mut total = 0i64;
        loop {
            let mut path: Vec<u32> = Vec::new(); // edge indices
            let mut u = s;
            loop {
                if u == t {
                    break;
                }
                let mut advanced = false;
                while self.it[u] < self.start[u + 1] {
                    let e = self.it[u] as usize;
                    let v = self.to[e] as usize;
                    if self.cap[e] > 0 && self.level[v] == self.level[u] + 1 {
                        path.push(e as u32);
                        u = v;
                        advanced = true;
                        break;
                    }
                    self.it[u] += 1;
                }
                if !advanced {
                    if u == s {
                        return total;
                    }
                    // retreat: dead end
                    self.level[u] = -1;
                    let e = path.pop().unwrap() as usize;
                    let back = self.head[e] as usize;
                    u = self.to[back] as usize;
                    self.it[u] += 1;
                }
            }
            let mut push = INF;
            for &e in &path {
                push = push.min(self.cap[e as usize]);
            }
            for &e in &path {
                self.cap[e as usize] -= push;
                let r = self.head[e as usize] as usize;
                self.cap[r] += push;
            }
            total += push;
        }
    }

    fn max_flow(&mut self, s: usize, t: usize) -> i64 {
        let mut flow = 0i64;
        while self.bfs(s, t) {
            for i in 0..self.n {
                self.it[i] = self.start[i];
            }
            flow += self.dfs(s, t);
        }
        flow
    }

    /// Source side of the min cut after `max_flow`.
    fn source_side(&self, s: usize) -> Vec<bool> {
        let mut seen = vec![false; self.n];
        let mut stack = vec![s];
        seen[s] = true;
        while let Some(u) = stack.pop() {
            for e in self.start[u]..self.start[u + 1] {
                let v = self.to[e as usize] as usize;
                if self.cap[e as usize] > 0 && !seen[v] {
                    seen[v] = true;
                    stack.push(v);
                }
            }
        }
        seen
    }
}

/// A closure problem: node `i` carries `pos[i]` unit gains (each "k ≥ j"
/// level), and an x-node carries costs `cc[i]` (∝ δ) and `ce[i]` (fixed part);
/// `edges (u,v)`: selecting `u` forces `v`.
struct Closure {
    pos: Vec<i64>,
    cc: Vec<i64>,
    ce: Vec<i64>,
    is_x: Vec<bool>,
    orbit: Vec<i64>,
    /// For the theorem check: `J(u)·orbit` on x-nodes (max charge received).
    charge_cap: Vec<i64>,
    edges: Vec<(u32, u32)>,
    den: i64,
}

#[derive(Clone, Debug)]
struct ClosureResult {
    delta_num: i128,
    delta_den: i128,
    delta: f64,
    dim: i64,
    sum_k: i64,
    iterations: usize,
    nodes: usize,
    edges: usize,
    /// `Σ_{u∈M} J(u)` — the charging upper bound on `Σk`.
    charge_bound: i64,
}

fn gcd(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

impl Closure {
    fn solve_at(&self, p: i128, q: i128) -> (i128, Vec<bool>) {
        let n = self.pos.len();
        let s = n;
        let t = n + 1;
        let mut e: Vec<(u32, u32, i64)> = Vec::with_capacity(self.edges.len() + n);
        let mut pos_sum: i128 = 0;
        for i in 0..n {
            let w: i128 = q * self.den as i128 * self.pos[i] as i128
                - p * self.cc[i] as i128
                - q * self.ce[i] as i128;
            if w > 0 {
                pos_sum += w;
                e.push((s as u32, i as u32, i64::try_from(w).expect("cap overflow")));
            } else if w < 0 {
                e.push((i as u32, t as u32, i64::try_from(-w).expect("cap overflow")));
            }
        }
        assert!(pos_sum < (INF / 2) as i128, "capacity sum overflow");
        for &(u, v) in &self.edges {
            e.push((u, v, INF));
        }
        let mut g = Dinic::new(n + 2, &e);
        let flow = g.max_flow(s, t) as i128;
        let side = g.source_side(s);
        (pos_sum - flow, side[..n].to_vec())
    }

    fn dinkelbach(&self) -> Option<ClosureResult> {
        let (mut p, mut q) = (0i128, 1i128);
        let mut iters = 0;
        let mut best: Option<ClosureResult> = None;
        loop {
            iters += 1;
            let (val, sel) = self.solve_at(p, q);
            if val <= 0 && best.is_some() {
                break;
            }
            let (mut k, mut csum, mut esum, mut dim, mut jb) = (0i128, 0i128, 0i128, 0i64, 0i64);
            for i in 0..self.pos.len() {
                if sel[i] {
                    k += self.pos[i] as i128;
                    csum += self.cc[i] as i128;
                    esum += self.ce[i] as i128;
                    if self.is_x[i] {
                        dim += self.orbit[i];
                        jb += self.charge_cap[i];
                    }
                }
            }
            if csum == 0 {
                return None;
            }
            let num = k * self.den as i128 - esum;
            let g = gcd(num, csum).max(1);
            let (np, nq) = (num / g, csum / g);
            if num <= 0 {
                return best;
            }
            if let Some(b) = &best {
                if np * b.delta_den <= b.delta_num * nq {
                    break;
                }
            }
            best = Some(ClosureResult {
                delta_num: np,
                delta_den: nq,
                delta: np as f64 / nq as f64,
                dim,
                sum_k: k as i64,
                iterations: iters,
                nodes: self.pos.len(),
                edges: self.edges.len(),
                charge_bound: jb,
            });
            p = np;
            q = nq;
            if iters > 60 {
                break;
            }
        }
        best
    }
}

/// Box-supported `f` (Newton polytope `[0,d]^m`, every lattice point of
/// `k·[0,d]^m` in `supp f^k`), leading vertex `a = d·1_top`, multiplicity
/// `ell`, bounding box `[0,L]^m`, bound exponents `λ_i = (δ cnum_i + enm_i)/den`.
/// Nodes `(β, r)`: "the side-`r` box anchored at `β` (towards `−` on top
/// coordinates, `+` otherwise) lies in `M`". Optional quotient by the
/// permutations preserving `(top_i, cnum_i, enm_i)` blocks.
#[derive(Clone, Debug)]
struct BoxSpec {
    m: usize,
    d: i64,
    ell: i64,
    l: i64,
    top: Vec<bool>,
    cnum: Vec<i64>,
    enm: Vec<i64>,
    den: i64,
    symmetric: bool,
}

fn multinomial_orbit(beta: &[i64], blocks: &[Vec<usize>]) -> i64 {
    let fact = |n: usize| (1..=n as i64).product::<i64>();
    let mut o = 1i64;
    for b in blocks {
        let mut counts: HashMap<i64, usize> = HashMap::new();
        for &i in b {
            *counts.entry(beta[i]).or_default() += 1;
        }
        let mut v = fact(b.len());
        for (_, c) in counts {
            v /= fact(c);
        }
        o *= v;
    }
    o
}

impl BoxSpec {
    fn blocks(&self) -> Vec<Vec<usize>> {
        if !self.symmetric {
            return (0..self.m).map(|i| vec![i]).collect();
        }
        let mut map: BTreeMap<(bool, i64, i64), Vec<usize>> = BTreeMap::new();
        for i in 0..self.m {
            map.entry((self.top[i], self.cnum[i], self.enm[i])).or_default().push(i);
        }
        map.into_values().collect()
    }

    fn build(&self) -> Closure {
        let blocks = self.blocks();
        let canon = |b: &mut Vec<i64>| {
            for bl in &blocks {
                let mut vals: Vec<i64> = bl.iter().map(|&i| b[i]).collect();
                vals.sort_unstable();
                for (j, &i) in bl.iter().enumerate() {
                    b[i] = vals[j];
                }
            }
        };
        let rmax = self.d * self.ell;
        let base = self.l + 1;
        let key = |b: &[i64], r: i64| -> u64 {
            let mut k = 0u64;
            for &x in b.iter().rev() {
                k = k * base as u64 + x as u64;
            }
            k * (rmax as u64 + 1) + r as u64
        };
        // enumerate canonical β
        let mut betas: Vec<Vec<i64>> = Vec::new();
        let mut cur = vec![0i64; self.m];
        loop {
            let mut c = cur.clone();
            canon(&mut c);
            if c == cur {
                betas.push(cur.clone());
            }
            let mut i = 0;
            loop {
                if i == self.m {
                    break;
                }
                cur[i] += 1;
                if cur[i] <= self.l {
                    break;
                }
                cur[i] = 0;
                i += 1;
            }
            if i == self.m {
                break;
            }
        }
        let valid = |b: &[i64], r: i64| -> bool {
            (0..self.m).all(|i| if self.top[i] { b[i] - r >= 0 } else { b[i] + r <= self.l })
        };
        let mut index: HashMap<u64, u32> = HashMap::new();
        let mut nodes: Vec<(usize, i64)> = Vec::new();
        for (bi, b) in betas.iter().enumerate() {
            for r in 0..=rmax {
                if valid(b, r) {
                    index.insert(key(b, r), nodes.len() as u32);
                    nodes.push((bi, r));
                }
            }
        }
        let n = nodes.len();
        let mut cl = Closure {
            pos: vec![0; n],
            cc: vec![0; n],
            ce: vec![0; n],
            is_x: vec![false; n],
            orbit: vec![0; n],
            charge_cap: vec![0; n],
            edges: Vec::new(),
            den: self.den,
        };
        let mut tmp = vec![0i64; self.m];
        for (id, &(bi, r)) in nodes.iter().enumerate() {
            let b = &betas[bi];
            let o = multinomial_orbit(b, &blocks);
            cl.orbit[id] = o;
            if r == 0 {
                cl.is_x[id] = true;
                cl.cc[id] = o * (0..self.m).map(|i| self.cnum[i] * b[i]).sum::<i64>();
                cl.ce[id] = o * (0..self.m).map(|i| self.enm[i] * b[i]).sum::<i64>();
                // charging bound with c = the all-d corner
                let j = (0..self.m).map(|i| b[i] / self.d).min().unwrap().min(self.ell);
                cl.charge_cap[id] = o * j;
            } else {
                if r % self.d == 0 {
                    cl.pos[id] = o;
                }
                let mut targets: Vec<u32> = Vec::with_capacity(1 << self.m);
                for mask in 0u32..(1 << self.m) {
                    for i in 0..self.m {
                        tmp[i] = b[i];
                        if (mask >> i) & 1 == 1 {
                            tmp[i] += if self.top[i] { -1 } else { 1 };
                        }
                    }
                    let mut c = tmp.clone();
                    canon(&mut c);
                    let t = *index.get(&key(&c, r - 1)).expect("sub-box node missing");
                    targets.push(t);
                }
                targets.sort_unstable();
                targets.dedup();
                for t in targets {
                    cl.edges.push((id as u32, t));
                }
            }
        }
        cl
    }
}

/// General support (exact sumset recursion): nodes `(β, j)` = "the shift
/// `x^{β−ja} f^j` has its support in `M`".
#[derive(Clone, Debug)]
struct GenSpec {
    supp: Vec<Vec<i64>>,
    a: Vec<i64>,
    ell: i64,
    l: Vec<i64>,
    cnum: Vec<i64>,
    enm: Vec<i64>,
    den: i64,
}

impl GenSpec {
    fn build(&self) -> Closure {
        let m = self.a.len();
        let mn: Vec<i64> = (0..m).map(|i| self.supp.iter().map(|s| s[i]).min().unwrap()).collect();
        let mx: Vec<i64> = (0..m).map(|i| self.supp.iter().map(|s| s[i]).max().unwrap()).collect();
        let valid = |b: &[i64], j: i64| -> bool {
            (0..m).all(|i| {
                b[i] + j * (mn[i] - self.a[i]) >= 0
                    && b[i] + j * (mx[i] - self.a[i]) <= self.l[i]
                    && b[i] >= j * self.a[i]
            })
        };
        let mut index: HashMap<(Vec<i64>, i64), u32> = HashMap::new();
        let mut nodes: Vec<(Vec<i64>, i64)> = Vec::new();
        let mut cur = vec![0i64; m];
        loop {
            for j in 0..=self.ell {
                if valid(&cur, j) {
                    index.insert((cur.clone(), j), nodes.len() as u32);
                    nodes.push((cur.clone(), j));
                }
            }
            let mut i = 0;
            loop {
                if i == m {
                    break;
                }
                cur[i] += 1;
                if cur[i] <= self.l[i] {
                    break;
                }
                cur[i] = 0;
                i += 1;
            }
            if i == m {
                break;
            }
        }
        let n = nodes.len();
        // charging corner: the support point maximising Σ cnum·s + enm·s at δ≈0.3
        let mut cl = Closure {
            pos: vec![0; n],
            cc: vec![0; n],
            ce: vec![0; n],
            is_x: vec![false; n],
            orbit: vec![1; n],
            charge_cap: vec![0; n],
            edges: Vec::new(),
            den: self.den,
        };
        for (id, (b, j)) in nodes.iter().enumerate() {
            if *j == 0 {
                cl.is_x[id] = true;
                cl.cc[id] = (0..m).map(|i| self.cnum[i] * b[i]).sum();
                cl.ce[id] = (0..m).map(|i| self.enm[i] * b[i]).sum();
                // J(u) for c = a (unconditional charging point)
                let jj = (0..m)
                    .filter(|&i| self.a[i] > 0)
                    .map(|i| b[i] / self.a[i])
                    .min()
                    .unwrap_or(self.ell)
                    .min(self.ell);
                cl.charge_cap[id] = jj;
            } else {
                cl.pos[id] = 1;
                let mut targets: Vec<u32> = self
                    .supp
                    .iter()
                    .map(|s| {
                        let t: Vec<i64> = (0..m).map(|i| b[i] + s[i] - self.a[i]).collect();
                        *index.get(&(t, j - 1)).expect("sumset node missing")
                    })
                    .collect();
                targets.sort_unstable();
                targets.dedup();
                for t in targets {
                    cl.edges.push((id as u32, t));
                }
            }
        }
        cl
    }
}

// ───────────────────────── closed-form box families ─────────────────────────

/// Box region `[0,L]^m`, box support `[0,d]^m`, any vertex: `Σk`, `Σ|β|`
/// (equal bounds). `k(β) = min(ℓ, min_i ⌊β_i/d⌋)` (or reflected).
fn box_region(m: u32, d: u64, ell: u64, l: u64) -> (f64, f64) {
    let mut sk = 0f64;
    for j in 1..=ell {
        let nj = (l + 1).saturating_sub(j * d) as f64;
        sk += nj.powi(m as i32);
    }
    let sb = m as f64 * ((l + 1) as f64).powi(m as i32 - 1) * (l as f64) * ((l + 1) as f64) / 2.0;
    (sk, sb)
}

fn jm_basic(m: u32, d: u64, ell: u64) -> f64 {
    let (k, b) = box_region(m, d, ell, d * ell);
    k / b
}

/// Best box extension `L = dℓ + t`, `t ∈ [0, 4dℓ]`.
fn jm_extended(m: u32, d: u64, ell: u64) -> (f64, u64) {
    let mut best = (0f64, 0u64);
    for t in 0..=(4 * d * ell) {
        let (k, b) = box_region(m, d, ell, d * ell + t);
        if k / b > best.0 {
            best = (k / b, t);
        }
    }
    best
}

/// Boneh–Durfee x+y-shift lattice (full rank, triangular), `ℓ` and `t`:
/// x-shifts `x^i f^k e^{ℓ−k}` (`i+k ≤ ℓ`), y-shifts `y^j f^k e^{ℓ−k}` (`1 ≤ j ≤ t`).
/// Returns the asymptotic-in-`e` bound on δ (`X = e^δ`, `Y = e^{1/2}`).
fn bd_xy(ell: i64, t: i64) -> f64 {
    let (mut sk, mut sx, mut sy) = (0f64, 0f64, 0f64);
    for k in 0..=ell {
        for i in 0..=(ell - k) {
            // leading monomial x^{i+k} y^k
            sk += k as f64;
            sx += (i + k) as f64;
            sy += k as f64;
        }
        for j in 1..=t {
            sk += k as f64;
            sx += k as f64;
            sy += (k + j) as f64;
        }
    }
    (sk - sy / 2.0) / sx
}

// ─────────────────────────────── parts ───────────────────────────────

fn rat(r: &ClosureResult) -> String {
    format!("{}/{}", r.delta_num, r.delta_den)
}

fn part_known(max_ell: i64, out: &mut Vec<Value>) {
    println!("## 1. Known bounds reproduced\n");
    // (a) univariate degree d
    println!("### 1a. Univariate `f` of degree `d` (Coppersmith 1/d)\n");
    println!("| d | ℓ | leading | closure-optimal δ_ℓ | dim | classic ℓ/(dℓ+d−1) |");
    println!("|---:|---:|:--|---:|---:|---:|");
    for d in [2i64, 3, 5] {
        for ell in [1i64, 2, 4, 8, 16] {
            for (lab, a) in [("x^d", d), ("1", 0)] {
                let spec = GenSpec {
                    supp: (0..=d).map(|e| vec![e]).collect(),
                    a: vec![a],
                    ell,
                    l: vec![d * ell + 2 * d],
                    cnum: vec![1],
                    enm: vec![0],
                    den: 1,
                };
                let r = spec.build().dinkelbach().unwrap();
                let classic = ell as f64 / (d * ell + d - 1) as f64;
                println!(
                    "| {d} | {ell} | {lab} | {:.5} ({}) | {} | {:.5} |",
                    r.delta,
                    rat(&r),
                    r.dim,
                    classic
                );
                out.push(json!({"part":"known","case":"univariate","d":d,"ell":ell,"leading":lab,
                    "delta":r.delta,"delta_rational":rat(&r),"dim":r.dim,"classic":classic,"limit":1.0/d as f64}));
            }
        }
    }
    // (b) Boneh–Durfee
    println!("\n### 1b. Boneh–Durfee `f(x,y) = 1 + x(A+y) mod e`, `X = e^δ`, `Y = e^{{1/2}}`\n");
    println!("Closed form x+y-shift lattice (best `t`) and the closure optimum over all monomial sets (leading `xy`).\n");
    println!("| ℓ | BD x+y shifts: best δ (t) | closure optimum δ_ℓ | dim | Σk ≤ charge bound |");
    println!("|---:|---:|---:|---:|:--|");
    let mut ells: Vec<i64> = vec![1, 2, 3, 4, 6, 8, 12, 16, 24, 32];
    ells.retain(|&e| e <= max_ell.max(8) * 4);
    for &ell in &ells {
        let mut best = (0f64, 0i64);
        for t in 0..=(2 * ell) {
            let v = bd_xy(ell, t);
            if v > best.0 {
                best = (v, t);
            }
        }
        let spec = GenSpec {
            supp: vec![vec![0, 0], vec![1, 0], vec![1, 1]],
            a: vec![1, 1],
            ell,
            l: vec![ell + 1, 2 * ell + 2],
            cnum: vec![2, 0],
            enm: vec![0, 1],
            den: 2,
        };
        let r = spec.build().dinkelbach().unwrap();
        println!(
            "| {ell} | {:.5} (t={}) | {:.5} | {} | {} ≤ {} |",
            best.0, best.1, r.delta, r.dim, r.sum_k, r.charge_bound
        );
        out.push(json!({"part":"known","case":"boneh_durfee","ell":ell,"bd_xy_best":best.0,"bd_xy_t":best.1,
            "closure_delta":r.delta,"closure_rational":rat(&r),"dim":r.dim}));
    }
    let mut asym = (0f64, 0f64);
    for tau in 0..=200 {
        let tau = tau as f64 / 200.0;
        let ell = 2000i64;
        let v = bd_xy(ell, (tau * ell as f64).round() as i64);
        if v > asym.0 {
            asym = (v, tau);
        }
    }
    let lit = (7.0 - 2.0 * 7f64.sqrt()) / 6.0;
    println!(
        "\nBD x+y shifts at ℓ = 2000, best τ = {:.3}: δ = {:.5}; literature (7−2√7)/6 = {:.5}. \
         The 0.2929 = 1−1/√2 bound needs a non-full-rank sublattice (Boneh–Durfee 1999) or the \
         unravelled linearisation `u = 1 + xy` (Herrmann–May 2010), outside this family.\n",
        asym.1, asym.0, lit
    );
    out.push(json!({"part":"known","case":"boneh_durfee_asymptotic","ell":2000,"tau":asym.1,"delta":asym.0,"literature":lit}));
    // (c) repo's JM hypercube numbers
    println!("### 1c. Jochemsz–May basic lattice on the box `[0,d]^m` (repo RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md §2.2)\n");
    println!("| m | d | ℓ=1 | ℓ=2 | ℓ=3 | ℓ=10 | ℓ=10⁴ | 2/(m(m+1)d) |");
    println!("|---:|---:|---:|---:|---:|---:|---:|---:|");
    for (m, d) in [(2u32, 2u64), (3, 4), (4, 8)] {
        let v: Vec<f64> = [1u64, 2, 3, 10, 10_000].iter().map(|&l| jm_basic(m, d, l)).collect();
        let lim = 2.0 / (m as f64 * (m as f64 + 1.0) * d as f64);
        println!(
            "| {m} | {d} | {:.5} | {:.5} | {:.5} | {:.5} | {:.5} | {:.5} |",
            v[0], v[1], v[2], v[3], v[4], lim
        );
        out.push(json!({"part":"known","case":"jm_basic_box","m":m,"d":d,"ell":[1,2,3,10,10000],"delta":v,"limit":lim}));
    }
    println!();
}

fn generic_weier(p: u64, seed: u64) -> (Weier, StdRng) {
    let mut rng = StdRng::seed_from_u64(seed);
    let f = Fp { p };
    let a = rng.gen_range(1..p);
    let b = rng.gen_range(1..p);
    (Weier { f, a, b }, rng)
}

fn largest_prime_below(n: u64) -> u64 {
    let mut p = n - 1;
    loop {
        if p % 2 == 1 && is_prime(p) {
            return p;
        }
        p -= 1;
    }
}

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for q in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % q == 0 {
            return n == q;
        }
    }
    let f = Fp { p: n };
    let mut d = n - 1;
    let mut s = 0;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'outer: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = f.pow(a, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = f.mul(x, x);
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

/// Square root mod p (Tonelli–Shanks), `None` for non-residues.
fn sqrt_mod(f: Fp, a: u64) -> Option<u64> {
    if a == 0 {
        return Some(0);
    }
    let p = f.p;
    if f.pow(a, (p - 1) / 2) != 1 {
        return None;
    }
    if p % 4 == 3 {
        return Some(f.pow(a, (p + 1) / 4));
    }
    let (mut q, mut s) = (p - 1, 0u32);
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }
    let mut z = 2u64;
    while f.pow(z, (p - 1) / 2) != p - 1 {
        z += 1;
    }
    let (mut m, mut c, mut t, mut r) = (s, f.pow(z, q), f.pow(a, q), f.pow(a, (q + 1) / 2));
    while t != 1 {
        let mut i = 0;
        let mut tt = t;
        while tt != 1 {
            tt = f.mul(tt, tt);
            i += 1;
        }
        let b = f.pow(c, 1u64 << (m - i - 1));
        m = i;
        c = f.mul(b, b);
        t = f.mul(t, c);
        r = f.mul(r, b);
    }
    Some(r)
}

/// Affine point arithmetic on a `Weier` curve (None = infinity).
fn pt_add(w: &Weier, p1: Option<(u64, u64)>, p2: Option<(u64, u64)>) -> Option<(u64, u64)> {
    let f = w.f;
    match (p1, p2) {
        (None, q) | (q, None) => q,
        (Some((x1, y1)), Some((x2, y2))) => {
            let lam = if x1 == x2 {
                if f.add(y1, y2) == 0 {
                    return None;
                }
                let num = f.add(f.mul(f.c(3), f.mul(x1, x1)), w.a);
                f.mul(num, f.inv(f.add(y1, y1)))
            } else {
                f.mul(f.sub(y2, y1), f.inv(f.sub(x2, x1)))
            };
            let x3 = f.sub(f.sub(f.mul(lam, lam), x1), x2);
            let y3 = f.sub(f.mul(lam, f.sub(x1, x3)), y1);
            Some((x3, y3))
        }
    }
}

fn random_point(w: &Weier, rng: &mut StdRng, bound: u64) -> (u64, u64) {
    loop {
        let x = rng.gen_range(1..=bound.min(w.f.p - 1));
        let rhs = w.f.add(w.f.add(w.f.mul(w.f.mul(x, x), x), w.f.mul(w.a, x)), w.b);
        if let Some(y) = sqrt_mod(w.f, rhs) {
            let y = if rng.gen_bool(0.5) { y } else { w.f.neg(y) };
            return (x, y);
        }
    }
}

struct SupportInfo {
    m: usize,
    dense: Dense,
    support: Vec<Vec<u8>>,
    e_support: Vec<Vec<u8>>,
}

fn analyse_support(m: usize, out: &mut Vec<Value>) -> SupportInfo {
    let p = largest_prime_below(1u64 << 61);
    let (w, mut rng) = generic_weier(p, 0x5EED_0000 + m as u64);
    let deg = 1usize << (m - 1);
    // target: x-coordinate of a random point
    let z = random_point(&w, &mut rng, p - 1).0;
    let t0 = Instant::now();
    let dense = s_with_target(&w, m, z);
    let interp_ms = t0.elapsed().as_secs_f64() * 1e3;
    // 1) the interpolant agrees with the direct evaluation off-grid (degree bound holds)
    let mut offgrid_ok = 0;
    for _ in 0..8 {
        let x: Vec<u64> = (0..m).map(|_| rng.gen_range(0..p)).collect();
        let mut full = x.clone();
        full.push(z);
        if dense.eval(w.f, &x) == w.s_n(m + 1, &full) {
            offgrid_ok += 1;
        }
    }
    // 2) vanishes at genuine decompositions R = P1 + … + Pm (fresh R each time)
    let mut vanish_ok = 0;
    for _ in 0..8 {
        let pts: Vec<(u64, u64)> = (0..m).map(|_| random_point(&w, &mut rng, p - 1)).collect();
        let mut r = None;
        for &q in &pts {
            r = pt_add(&w, r, Some(q));
        }
        let Some((zr, _)) = r else { continue };
        let xs: Vec<u64> = pts.iter().map(|q| q.0).collect();
        let mut full = xs.clone();
        full.push(zr);
        if w.s_n(m + 1, &full) == 0 {
            vanish_ok += 1;
        }
    }
    // and the interpolated f at the planted decomposition of z itself
    let mut planted_ok = 0;
    for _ in 0..4 {
        // decompose z: pick P1..P_{m−1}, set P_m = R − ΣP_i where R has x = z
        let ry = sqrt_mod(w.f, w.f.add(w.f.add(w.f.mul(w.f.mul(z, z), z), w.f.mul(w.a, z)), w.b))
            .unwrap();
        let rpt = (z, ry);
        let pts: Vec<(u64, u64)> = (0..m - 1).map(|_| random_point(&w, &mut rng, p - 1)).collect();
        let mut acc = Some(rpt);
        for &q in &pts {
            acc = pt_add(&w, acc, Some((q.0, w.f.neg(q.1))));
        }
        let Some(last) = acc else { continue };
        let mut xs: Vec<u64> = pts.iter().map(|q| q.0).collect();
        xs.push(last.0);
        if dense.eval(w.f, &xs) == 0 {
            planted_ok += 1;
        }
    }
    let support = dense.support();
    let sp = dense.to_sparse();
    // symmetry
    let symmetric = sp.keys().all(|e| {
        let mut s = e.clone();
        s.sort_unstable_by(|a, b| b.cmp(a));
        sp.get(&s) == sp.get(e)
    });
    let corners_present = (0u32..(1 << m))
        .filter(|mask| {
            let e: Vec<u8> = (0..m).map(|i| if (mask >> i) & 1 == 1 { deg as u8 } else { 0 }).collect();
            sp.contains_key(&e)
        })
        .count();
    let max_total = support.iter().map(|e| e.iter().map(|&x| x as usize).sum::<usize>()).max().unwrap();
    let units_present = (0..m).all(|i| {
        let mut e = vec![0u8; m];
        e[i] = 1;
        sp.contains_key(&e)
    }) && sp.contains_key(&vec![0u8; m]);
    // symmetric rewriting
    let t1 = Instant::now();
    let g = symmetric_reduce(w.f, m, &sp).expect("symmetric reduction");
    let sym_ms = t1.elapsed().as_secs_f64() * 1e3;
    let e_support: Vec<Vec<u8>> = g.keys().cloned().collect();
    let e_wdeg = e_support
        .iter()
        .map(|e| e.iter().enumerate().map(|(i, &x)| (i + 1) * x as usize).sum::<usize>())
        .max()
        .unwrap();
    let e_maxdeg: Vec<u8> = (0..m).map(|i| e_support.iter().map(|e| e[i]).max().unwrap()).collect();
    println!(
        "| S{} | {} | {}^{} = {} | {} | {}/{} | {} | {} | {} | {}/8 | {}/8 | {}/4 | {} | {:?} | {} | {:.0} / {:.0} |",
        m + 1,
        m,
        deg + 1,
        m,
        (deg + 1).pow(m as u32),
        support.len(),
        corners_present,
        1 << m,
        max_total,
        if symmetric { "yes" } else { "NO" },
        if units_present { "yes" } else { "no" },
        offgrid_ok,
        vanish_ok,
        planted_ok,
        e_support.len(),
        e_maxdeg,
        e_wdeg,
        interp_ms,
        sym_ms
    );
    out.push(json!({"part":"support","poly":format!("S{}",m+1),"m":m,"deg_per_var":deg,"p":p,"curve_a":w.a,"curve_b":w.b,"z":z,
        "box_size":(deg+1).pow(m as u32),"support_size":support.len(),"box_vertices_present":corners_present,
        "max_total_degree":max_total,"symmetric":symmetric,"has_1_and_units":units_present,
        "offgrid_agree_of_8":offgrid_ok,"semaev_vanishes_of_8":vanish_ok,"interpolant_vanishes_at_planted_of_4":planted_ok,
        "e_support_size":e_support.len(),"e_max_degree":e_maxdeg,"e_weighted_degree":e_wdeg}));
    SupportInfo { m, dense, support, e_support }
}

fn part_families(out: &mut Vec<Value>) {
    println!("## 3. Closed-form box strategies (any leading vertex; equal bounds)\n");
    println!("JM basic = `M = [0,dℓ]^m`; JM extended = best `M = [0,dℓ+t]^m`. Upper bound for the whole family: 1/(md).\n");
    println!("| poly | m | d | ℓ | JM basic δ_ℓ | JM extended δ_ℓ (t) | 1/(md) |");
    println!("|:--|---:|---:|---:|---:|---:|---:|");
    for (name, m, d) in [("S3", 2u32, 2u64), ("S4", 3, 4), ("S5", 4, 8)] {
        for ell in [1u64, 2, 3, 4, 6, 8, 12, 16, 24, 32, 64, 128] {
            let b = jm_basic(m, d, ell);
            let (e, t) = jm_extended(m, d, ell);
            println!(
                "| {name} | {m} | {d} | {ell} | {b:.6} | {e:.6} ({t}) | {:.6} |",
                1.0 / (m as f64 * d as f64)
            );
            out.push(json!({"part":"families","poly":name,"m":m,"d":d,"ell":ell,"jm_basic":b,"jm_extended":e,"t":t,
                "upper_bound":1.0/(m as f64 * d as f64)}));
        }
    }
    println!();
}

fn part_closure(max_ell: i64, infos: &[SupportInfo], out: &mut Vec<Value>) {
    println!("## 4. Exact closure optimum (all monomial sets in the bounding box)\n");
    println!("Box support (Newton polytope `[0,d]^m`, verified in §2), every vertex class `a = d·1_S`, `|S| = s`. ");
    println!("Bounding box `[0,L]^m`. `best box ≤ L` = max over `L' ≤ L` of the box strategy `M = [0,L']^m` (closed form); ");
    println!("`bounded` = that best box is `[0,L]^m` itself (the optimum would grow with a larger bounding box). ");
    println!("`Σk ≤ charge` is the machine check of the 1/(md) argument on the optimum.\n");
    println!("| poly | ℓ | vertex s | box L | δ_ℓ (exact) | δ_ℓ | best box ≤ L (L') | closure = best box | dim | Σk ≤ charge | nodes | edges | iters | ms |");
    println!("|:--|---:|---:|---:|:--|---:|---:|:--|---:|:--|---:|---:|---:|---:|");
    let mut jobs: Vec<(&str, usize, i64, i64, i64, Vec<usize>)> = Vec::new();
    for ell in 1..=(2 * max_ell.max(6)) {
        jobs.push(("S3", 2, 2, ell, 4 * ell, vec![0, 1, 2]));
    }
    for ell in 1..=max_ell.max(5) {
        let t = jm_extended(3, 4, ell as u64).1 as i64;
        jobs.push(("S4", 3, 4, ell, 4 * ell + t + 4, vec![0, 1, 2, 3]));
    }
    if max_ell >= 1 {
        let t1 = jm_extended(4, 8, 1).1 as i64;
        jobs.push(("S5", 4, 8, 1, 8 + t1 + 6, vec![0, 4]));
        jobs.push(("S5", 4, 8, 1, 24, vec![0, 1, 2, 3, 4]));
    }
    if max_ell >= 2 {
        jobs.push(("S5", 4, 8, 2, 32, vec![0, 4]));
        jobs.push(("S5", 4, 8, 2, 24, vec![1, 2, 3, 4]));
    }
    if max_ell >= 3 {
        jobs.push(("S5", 4, 8, 3, 32, vec![4]));
    }
    for (name, m, d, ell, l, classes) in jobs {
        let (mut bb, mut bl) = (0f64, 0i64);
        for lp in 0..=l {
            let (k, b) = box_region(m as u32, d as u64, ell as u64, lp as u64);
            if b > 0.0 && k / b > bb {
                bb = k / b;
                bl = lp;
            }
        }
        for s in classes {
            let top: Vec<bool> = (0..m).map(|i| i < s).collect();
            let spec = BoxSpec { m, d, ell, l, top, cnum: vec![1; m], enm: vec![0; m], den: 1, symmetric: true };
            let t0 = Instant::now();
            let cl = spec.build();
            let r = cl.dinkelbach().unwrap();
            let ms = t0.elapsed().as_secs_f64() * 1e3;
            assert!(r.sum_k <= r.charge_bound, "charging bound violated");
            let same = (r.delta - bb).abs() <= 1e-12 * bb.max(1e-300);
            let bounded = bl == l;
            println!(
                "| {name} | {ell} | {s} | {l} | {} | {:.6} | {:.6} ({bl}{}) | {} | {} | {} ≤ {} | {} | {} | {} | {:.0} |",
                rat(&r), r.delta, bb, if bounded { ", bounded" } else { "" }, if same { "yes" } else { "NO" },
                r.dim, r.sum_k, r.charge_bound, r.nodes, r.edges, r.iterations, ms
            );
            out.push(json!({"part":"closure","encoding":"x","poly":name,"m":m,"d":d,"ell":ell,"vertex_s":s,"L":l,
                "delta":r.delta,"delta_rational":rat(&r),"best_box_delta":bb,"best_box_L":bl,"bounded":bounded,
                "closure_equals_best_box":same,"dim":r.dim,"sum_k":r.sum_k,"charge_bound":r.charge_bound,
                "nodes":r.nodes,"edges":r.edges,"iterations":r.iterations,"ms":ms,"jm_basic":jm_basic(m as u32,d as u64,ell as u64)}));
        }
    }
    // symmetric (e₁..e_m) encoding, general sumset engine, weights δ·i
    println!("\n### 4b. Symmetric encoding `G(e₁..e_m)`, bounds `|e_i| ≲ X^i` (weights `i·δ`)\n");
    println!("Leading monomial: every vertex of `Newt(G)` found by 64 random weight directions. Upper bound 1/(weighted degree).\n");
    println!("| poly | ℓ | #vertices tried | best vertex | δ_ℓ | dim | nodes | edges | ms |");
    println!("|:--|---:|---:|:--|---:|---:|---:|---:|---:|");
    for info in infos {
        let m = info.m;
        let supp: Vec<Vec<i64>> = info.e_support.iter().map(|e| e.iter().map(|&x| x as i64).collect()).collect();
        let mx: Vec<i64> = (0..m).map(|i| supp.iter().map(|s| s[i]).max().unwrap()).collect();
        // vertices via random directions (any unique maximiser is a vertex)
        let mut rng = StdRng::seed_from_u64(0xE5 + m as u64);
        let mut verts: Vec<Vec<i64>> = Vec::new();
        for _ in 0..64 {
            let wv: Vec<f64> = (0..m).map(|_| rng.gen_range(-1.0..1.0)).collect();
            let score = |s: &Vec<i64>| s.iter().zip(&wv).map(|(&a, &b)| a as f64 * b).sum::<f64>();
            let best = supp.iter().max_by(|a, b| score(a).partial_cmp(&score(b)).unwrap()).unwrap();
            if !verts.contains(best) {
                verts.push(best.clone());
            }
        }
        verts.sort();
        let ell_cap = match m {
            2 => 4,
            3 => 2,
            _ => 1,
        };
        for ell in 1..=ell_cap.min(max_ell.max(1)) {
            let t0 = Instant::now();
            let mut best: Option<(ClosureResult, Vec<i64>)> = None;
            let mut nodes = 0;
            let mut edges = 0;
            for a in &verts {
                let l: Vec<i64> = mx.iter().map(|&x| x * ell + x).collect();
                // keep the node count bounded for S5
                let size: f64 = l.iter().map(|&x| (x + 1) as f64).product();
                if size > 3.0e6 {
                    continue;
                }
                let spec = GenSpec {
                    supp: supp.clone(),
                    a: a.clone(),
                    ell,
                    l,
                    cnum: (1..=m as i64).collect(),
                    enm: vec![0; m],
                    den: 1,
                };
                let cl = spec.build();
                nodes = nodes.max(cl.pos.len());
                edges = edges.max(cl.edges.len());
                if let Some(r) = cl.dinkelbach() {
                    if best.as_ref().map_or(true, |b| r.delta > b.0.delta) {
                        best = Some((r, a.clone()));
                    }
                }
            }
            let ms = t0.elapsed().as_secs_f64() * 1e3;
            if let Some((r, a)) = best {
                println!(
                    "| S{} | {ell} | {} | {:?} | {:.6} | {} | {} | {} | {:.0} |",
                    m + 1, verts.len(), a, r.delta, r.dim, nodes, edges, ms
                );
                out.push(json!({"part":"closure","encoding":"e","poly":format!("S{}",m+1),"m":m,"ell":ell,
                    "vertices":verts.len(),"best_vertex":a,"delta":r.delta,"delta_rational":rat(&r),"dim":r.dim,"ms":ms}));
            } else {
                println!("| S{} | {ell} | {} | – | skipped (box > 3e6 nodes) | | | | |", m + 1, verts.len());
                out.push(json!({"part":"closure","encoding":"e","poly":format!("S{}",m+1),"m":m,"ell":ell,"skipped":true}));
            }
        }
    }
    println!();
}

/// Finite lattice: plant a root of `S_{m+1}(·, z)` with coordinates ≤ B on a
/// 32-bit curve, build the JM basic lattice at multiplicity `ell`, LLL, count
/// reduced rows that vanish at the root over Z.
fn lattice_trial(
    curve: &FastCurve,
    m: usize,
    ell: i64,
    bound: u64,
    rng: &mut StdRng,
) -> Option<(usize, bool, usize, f64)> {
    let f = Fp { p: curve.f.p };
    let w = Weier { f, a: curve.f.from_mont(curve.a), b: curve.f.from_mont(curve.b) };
    let d = 1i64 << (m - 1);
    let pts: Vec<(u64, u64)> = (0..m).map(|_| random_point(&w, rng, bound)).collect();
    let mut r = None;
    for &q in &pts {
        r = pt_add(&w, r, Some(q));
    }
    let (z, _) = r?;
    let root: Vec<BigInt> = pts.iter().map(|q| BigInt::from(q.0)).collect();
    let poly = s_with_target(&w, m, z);
    // normalise so the corner coefficient is 1
    let corner_idx: usize = (0..m).map(|i| (d as usize) * (d as usize + 1).pow(i as u32)).sum();
    let lc = poly.c[corner_idx];
    if lc == 0 {
        return None;
    }
    let inv = f.inv(lc);
    let fpoly: HashMap<Vec<i64>, BigInt> = (0..poly.c.len())
        .filter(|&i| poly.c[i] != 0)
        .map(|i| {
            let e: Vec<i64> = poly.exps(i).iter().map(|&x| x as i64).collect();
            (e, BigInt::from(f.mul(poly.c[i], inv)))
        })
        .collect();
    let eval = |pp: &HashMap<Vec<i64>, BigInt>| -> BigInt {
        let mut acc = BigInt::zero();
        for (e, c) in pp {
            let mut t = c.clone();
            for i in 0..m {
                t *= root[i].pow(e[i] as u32);
            }
            acc += t;
        }
        acc
    };
    let pz = BigInt::from(f.p);
    debug_assert!((eval(&fpoly) % &pz).is_zero());
    let mul = |a: &HashMap<Vec<i64>, BigInt>, b: &HashMap<Vec<i64>, BigInt>| {
        let mut o: HashMap<Vec<i64>, BigInt> = HashMap::new();
        for (ea, ca) in a {
            for (eb, cb) in b {
                let e: Vec<i64> = ea.iter().zip(eb).map(|(x, y)| x + y).collect();
                *o.entry(e).or_insert_with(BigInt::zero) += ca * cb;
            }
        }
        o.retain(|_, c| !c.is_zero());
        o
    };
    let mut fpow = vec![HashMap::from([(vec![0i64; m], BigInt::from(1))])];
    for k in 1..=ell {
        let next = mul(&fpow[k as usize - 1], &fpoly);
        fpow.push(next);
    }
    let side = d * ell + 1;
    let dim = (side as usize).pow(m as u32);
    let mono = |e: &[i64]| -> usize {
        let mut k = 0usize;
        for &x in e.iter().rev() {
            k = k * side as usize + x as usize;
        }
        k
    };
    let bb = BigInt::from(bound);
    let mut basis: Vec<Vec<BigInt>> = Vec::with_capacity(dim);
    for idx in 0..dim {
        let mut beta = vec![0i64; m];
        let mut rr = idx;
        for x in beta.iter_mut() {
            *x = (rr % side as usize) as i64;
            rr /= side as usize;
        }
        let k = beta.iter().map(|&x| x / d).min().unwrap().min(ell);
        let shift: Vec<i64> = beta.iter().map(|&x| x - k * d).collect();
        let pk = pz.pow((ell - k) as u32);
        let mut row = vec![BigInt::zero(); dim];
        for (e, c) in &fpow[k as usize] {
            let ee: Vec<i64> = e.iter().zip(&shift).map(|(x, y)| x + y).collect();
            let tot: i64 = ee.iter().sum();
            row[mono(&ee)] = c * &pk * bb.pow(tot as u32);
        }
        basis.push(row);
    }
    let t0 = Instant::now();
    // f64 Gram–Schmidt is documented sound to ~30 dimensions; use the
    // high-precision variant above that.
    let reduced = if dim <= 30 { lll_reduce(&mut basis, 0.99) } else { lll_reduce_hp(&mut basis, 0.99) };
    if reduced.is_err() {
        return None;
    }
    let ms = t0.elapsed().as_secs_f64() * 1e3;
    let pm = pz.pow(ell as u32);
    let mut vanishing = 0usize;
    let mut hg_first = true;
    for (ri, row) in basis.iter().enumerate() {
        let mut pp: HashMap<Vec<i64>, BigInt> = HashMap::new();
        let mut norm2 = BigInt::zero();
        for (idx, c) in row.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            norm2 += c * c;
            let mut beta = vec![0i64; m];
            let mut rr = idx;
            for x in beta.iter_mut() {
                *x = (rr % side as usize) as i64;
                rr /= side as usize;
            }
            let tot: i64 = beta.iter().sum();
            pp.insert(beta, c / bb.pow(tot as u32));
        }
        if pp.is_empty() {
            continue;
        }
        if eval(&pp).is_zero() {
            vanishing += 1;
        }
        if ri < m {
            let norm = norm2.to_f64().unwrap_or(f64::INFINITY).sqrt();
            hg_first &= norm * (dim as f64).sqrt() < pm.to_f64().unwrap_or(f64::INFINITY);
        }
    }
    Some((vanishing, hg_first, dim, ms))
}

fn part_lattice(out: &mut Vec<Value>) {
    println!("## 5. Finite lattices on a 32-bit curve (repo LLL, δ_LLL = 0.99)\n");
    let (curve, _) = find_a3_curve(32, 3).expect("32-bit curve");
    println!(
        "Curve `{}`: p = {}, a = −3, b = {}, prime order n = {}. Planted `R = P₁ + … + P_m`, \
         `x(P_i) ≤ B = ⌊p^δ⌋`; success = at least `m` reduced rows vanish at the root over Z \
         (the count needed to eliminate `m` unknowns). JM basic lattice `M = [0, dℓ]^m`.\n",
        curve.name,
        curve.f.p,
        curve.f.from_mont(curve.b),
        curve.n
    );
    println!("| poly | ℓ | dim | δ | B | ≥m vanishing / trials | HG on first m rows | mean vanishing rows | ms per LLL |");
    println!("|:--|---:|---:|---:|---:|---:|---:|---:|---:|");
    let cells: Vec<(usize, i64, Vec<f64>, usize)> = vec![
        (2, 1, vec![0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.20], 8),
        (2, 2, vec![0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.20], 8),
        (3, 1, vec![0.035, 0.05, 0.065, 0.08, 0.10, 0.125], 6),
    ];
    for (m, ell, deltas, trials) in cells {
        let mut zero_streak = 0;
        for &delta in &deltas {
            let bound = ((curve.f.p as f64).powf(delta)).floor().max(2.0) as u64;
            let mut rng = StdRng::seed_from_u64(0x10_0000 + (m as u64) * 1000 + ell as u64 * 100 + (delta * 1000.0) as u64);
            let (mut succ, mut hg, mut vsum, mut done, mut dim, mut ms) = (0, 0, 0, 0, 0, 0f64);
            for _ in 0..trials {
                if let Some((v, h, dd, t)) = lattice_trial(&curve, m, ell, bound, &mut rng) {
                    done += 1;
                    vsum += v;
                    dim = dd;
                    ms += t;
                    if v >= m {
                        succ += 1;
                    }
                    if h {
                        hg += 1;
                    }
                }
            }
            let mean_ms = ms / done.max(1) as f64;
            println!(
                "| S{} | {ell} | {dim} | {delta:.3} | {bound} | {succ}/{done} | {hg}/{done} | {:.1} | {mean_ms:.0} |",
                m + 1,
                vsum as f64 / done.max(1) as f64
            );
            out.push(json!({"part":"lattice","poly":format!("S{}",m+1),"m":m,"ell":ell,"dim":dim,"delta":delta,"bound":bound,
                "trials":done,"successes":succ,"hg_first_m":hg,"mean_vanishing":vsum as f64/done.max(1) as f64,"ms_per_lll":mean_ms,
                "p":curve.f.p,"curve":curve.name}));
            if succ == 0 && m == 2 {
                zero_streak += 1;
                if zero_streak >= 2 {
                    break;
                }
            } else {
                zero_streak = 0;
            }
        }
    }
    println!();
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).cloned();
    let part = get("--part").unwrap_or_else(|| "all".into());
    let max_ell: i64 = get("--max-ell").and_then(|s| s.parse().ok()).unwrap_or(3);
    let mut out: Vec<Value> = Vec::new();
    let t0 = Instant::now();
    println!("# lit2026 W10: Coppersmith reach for the prime-field S₃/S₄/S₅ small-root oracle\n");
    let all = part == "all";
    if all || part == "known" {
        part_known(max_ell, &mut out);
    }
    let mut infos = Vec::new();
    if all || part == "support" || part == "closure" {
        println!("## 2. Supports of `f = S_{{m+1}}(x₁..x_m, z)` over a 61-bit prime field\n");
        println!("| poly | m | box | #monomials | box vertices | max total deg | symmetric | 1, x_i present | off-grid = direct | S vanishes at decompositions | f vanishes at planted | #e-monomials | e max degrees | e weighted deg | ms interp / sym |");
        println!("|:--|---:|---:|---:|---:|---:|:--|:--|---:|---:|---:|---:|:--|---:|:--|");
        for m in 2..=4 {
            infos.push(analyse_support(m, &mut out));
        }
        for info in &infos {
            let deg = 1usize << (info.m - 1);
            assert!(info.support.iter().all(|e| e.iter().all(|&x| x as usize <= deg)));
            let _ = &info.dense;
        }
        println!();
    }
    if all || part == "families" {
        part_families(&mut out);
    }
    if all || part == "closure" {
        part_closure(max_ell, &infos, &mut out);
    }
    if all || part == "lattice" {
        part_lattice(&mut out);
    }
    println!("total wall {:.1} s", t0.elapsed().as_secs_f64());
    if let Some(path) = get("--json") {
        std::fs::write(&path, serde_json::to_string_pretty(&json!({
            "schema":"lit2026-w10-coppersmith/v1","part":part,"max_ell":max_ell,"records":out}))
        .unwrap())
        .unwrap();
        eprintln!("wrote {path}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn univariate_closure_reaches_classic_and_tends_to_one_over_d() {
        for d in [2i64, 3] {
            let mut last = 0.0;
            for ell in [1i64, 2, 4, 8] {
                let spec = GenSpec {
                    supp: (0..=d).map(|e| vec![e]).collect(),
                    a: vec![d],
                    ell,
                    l: vec![d * ell + 2 * d],
                    cnum: vec![1],
                    enm: vec![0],
                    den: 1,
                };
                let r = spec.build().dinkelbach().unwrap();
                let classic = ell as f64 / (d * ell + d - 1) as f64;
                assert!(r.delta >= classic - 1e-12, "d={d} ell={ell}: {} < {classic}", r.delta);
                assert!(r.delta < 1.0 / d as f64);
                assert!(r.delta >= last - 1e-12);
                last = r.delta;
            }
            assert!(last > 0.9 / d as f64);
        }
    }

    #[test]
    fn boneh_durfee_xy_shifts_reproduce_0_2847() {
        let mut best = 0f64;
        for tau in 0..=100 {
            let ell = 1000;
            best = best.max(bd_xy(ell, tau * ell / 100));
        }
        let lit = (7.0 - 2.0 * 7f64.sqrt()) / 6.0;
        assert!((best - lit).abs() < 2e-3, "{best} vs {lit}");
    }

    #[test]
    fn jm_basic_matches_repo_numbers_and_limit() {
        assert!((jm_basic(2, 2, 1) - 1.0 / 18.0).abs() < 1e-12);
        assert!((jm_basic(2, 2, 2) - 1.0 / 10.0).abs() < 1e-12);
        assert!((jm_basic(2, 2, 3) - 0.119).abs() < 1e-3);
        for (m, d) in [(2u32, 2u64), (3, 4), (4, 8)] {
            let lim = 2.0 / (m as f64 * (m as f64 + 1.0) * d as f64);
            assert!((jm_basic(m, d, 20_000) - lim).abs() < 1e-3 * lim + 1e-6);
        }
    }

    #[test]
    fn box_engine_equals_sumset_engine_on_s3_box() {
        for ell in 1..=3i64 {
            for top in [vec![true, true], vec![false, false], vec![true, false]] {
                let l = 2 * ell + 2;
                let b = BoxSpec {
                    m: 2,
                    d: 2,
                    ell,
                    l,
                    top: top.clone(),
                    cnum: vec![1, 1],
                    enm: vec![0, 0],
                    den: 1,
                    symmetric: false,
                }
                .build()
                .dinkelbach()
                .unwrap();
                let bs = BoxSpec { symmetric: true, ..BoxSpec { m: 2, d: 2, ell, l, top: top.clone(), cnum: vec![1, 1], enm: vec![0, 0], den: 1, symmetric: false } }
                    .build()
                    .dinkelbach()
                    .unwrap();
                let a: Vec<i64> = top.iter().map(|&t| if t { 2 } else { 0 }).collect();
                let supp: Vec<Vec<i64>> = (0..3).flat_map(|i| (0..3).map(move |j| vec![i, j])).collect();
                let g = GenSpec { supp, a, ell, l: vec![l, l], cnum: vec![1, 1], enm: vec![0, 0], den: 1 }
                    .build()
                    .dinkelbach()
                    .unwrap();
                assert_eq!((b.delta_num, b.delta_den), (g.delta_num, g.delta_den), "ell={ell} top={top:?}");
                if top[0] == top[1] {
                    assert_eq!((bs.delta_num, bs.delta_den), (b.delta_num, b.delta_den));
                }
            }
        }
    }

    #[test]
    fn closure_never_beats_jm_box_from_below_and_respects_charge_bound() {
        // the closure optimum contains the JM basic set, so δ_ℓ ≥ JM basic
        for ell in 1..=3i64 {
            let r = BoxSpec {
                m: 3,
                d: 4,
                ell,
                l: 4 * ell + 4,
                top: vec![true; 3],
                cnum: vec![1; 3],
                enm: vec![0; 3],
                den: 1,
                symmetric: true,
            }
            .build()
            .dinkelbach()
            .unwrap();
            assert!(r.delta >= jm_basic(3, 4, ell as u64) - 1e-12);
            assert!(r.delta < 1.0 / 12.0);
            assert!(r.sum_k <= r.charge_bound);
        }
    }

    #[test]
    fn semaev_s4_s5_vanish_at_decompositions_and_have_box_support() {
        let p = largest_prime_below(1u64 << 40);
        let (w, mut rng) = generic_weier(p, 7);
        for m in [2usize, 3] {
            let pts: Vec<(u64, u64)> = (0..m).map(|_| random_point(&w, &mut rng, p - 1)).collect();
            let mut r = None;
            for &q in &pts {
                r = pt_add(&w, r, Some(q));
            }
            let z = r.unwrap().0;
            let f = s_with_target(&w, m, z);
            let xs: Vec<u64> = pts.iter().map(|q| q.0).collect();
            assert_eq!(f.eval(w.f, &xs), 0);
            let deg = 1u8 << (m - 1);
            let s = f.to_sparse();
            assert!(s.contains_key(&vec![deg; m]));
            assert!(s.contains_key(&vec![0u8; m]));
            // symmetric rewriting re-expands exactly (remainder reached zero)
            let g = symmetric_reduce(w.f, m, &s).unwrap();
            assert!(!g.is_empty());
        }
        // S5 vanishing at a genuine 5-term relation
        let pts: Vec<(u64, u64)> = (0..4).map(|_| random_point(&w, &mut rng, p - 1)).collect();
        let mut r = None;
        for &q in &pts {
            r = pt_add(&w, r, Some(q));
        }
        let z = r.unwrap().0;
        assert_eq!(w.s5(&[pts[0].0, pts[1].0, pts[2].0, pts[3].0, z]), 0);
        let mut bad = [pts[0].0, pts[1].0, pts[2].0, pts[3].0, z];
        bad[0] = w.f.add(bad[0], 1);
        assert_ne!(w.s5(&bad), 0);
    }
}
