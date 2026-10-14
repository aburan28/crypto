//! IC next-10 item 05: Coppersmith reach for several polynomials and shared
//! targets (extends `lit2026_w10_coppersmith_bound`).
//!
//! **Family.** Unknowns `x_v` with bounds `p^{λ_v}` (`λ_v = δ` small, `1` large).
//! Polynomials `f_1..f_r` mod p with box supports (degree `d_{jv}` in each of
//! their variables, dense, as Semaev's `S_{m+1}(·, z)` is). Leading monomial of
//! each `f_j` = its top corner `a_j`; a common graded order makes every product
//! triangular. Shifts `p^{ℓ−|k|} x^γ Π f_j^{k_j}` with leading monomial
//! `β = γ + Σ k_j a_j`; for a down-closed monomial set `M` every such shift has
//! support in `M`, so the best gain at `β` is
//! `k(β) = max{|k| : |k| ≤ ℓ, Σ_{j∋v} k_j d_{jv} ≤ β_v ∀v}`.
//! Howgrave-Graham at `p → ∞`: `Σ_M k(β) > Σ_M λ·β`.
//!
//! **Exact optimiser.** Maximum-weight down-closed `M` = min-cut (Dinic) on
//! symmetry orbits, inside Dinkelbach iteration on the ratio
//! `δ = Σ(k − Σ_large β)/Σ_small β`; returns the best `M` (an exact rational)
//! and certifies that no `M` beats it by `2^-30` or more.
//!
//! **LLL part.** Builds the closure-optimal lattice for `S₃` alone and for the
//! shared-target pair `S₃(x₁,x₂,z₁), S₃(x₁,x₃,z₂)` on a toy prime and counts
//! reduced rows vanishing over Z at planted roots.
//!
//! ```bash
//! cargo run --release --example ic_next10_05_multipoly_lattice -- --part lp --json lp.json
//! cargo run --release --example ic_next10_05_multipoly_lattice -- --part lll --bits 40 --ell 1 --json lll.json
//! ```

use std::collections::HashMap;
use std::time::Instant;

use crypto_lib::cryptanalysis::lattice::{lll_reduce, lll_reduce_hp};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

// ───────────────────────────── Dinic (as in W10) ─────────────────────────────

struct Dinic {
    n: usize,
    rev: Vec<u32>,
    to: Vec<u32>,
    cap: Vec<i64>,
    level: Vec<i32>,
    it: Vec<u32>,
    start: Vec<u32>,
}

const INF: i64 = i64::MAX / 4;

impl Dinic {
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
        Dinic { n, rev, to, cap, level: vec![0; n], it: vec![0; n], start }
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
        let mut total = 0i64;
        loop {
            let mut path: Vec<u32> = Vec::new();
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
                    self.level[u] = -1;
                    let e = path.pop().unwrap() as usize;
                    let back = self.rev[e] as usize;
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
                let r = self.rev[e as usize] as usize;
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

// ───────────────────────────── system description ────────────────────────────

#[derive(Clone)]
struct Sys {
    name: String,
    nvars: usize,
    large: Vec<bool>,
    /// interchangeable variable groups (the system is invariant under
    /// permuting the variables of a group, together with the polynomials)
    groups: Vec<Vec<usize>>,
    /// per polynomial: (variable, degree) of its box support
    polys: Vec<Vec<(usize, i64)>>,
}

impl Sys {
    /// all k-vectors with 1 ≤ |k| ≤ ℓ
    fn kvecs(&self, ell: i64) -> Vec<Vec<i64>> {
        let r = self.polys.len();
        let mut out = Vec::new();
        let mut k = vec![0i64; r];
        loop {
            let s: i64 = k.iter().sum();
            if s >= 1 && s <= ell {
                out.push(k.clone());
            }
            let mut i = 0;
            loop {
                if i == r {
                    return out;
                }
                k[i] += 1;
                if k[i] <= ell {
                    break;
                }
                k[i] = 0;
                i += 1;
            }
        }
    }
    /// top-corner exponent vector of poly j
    fn corner(&self, j: usize) -> Vec<i64> {
        let mut a = vec![0i64; self.nvars];
        for &(v, d) in &self.polys[j] {
            a[v] = d;
        }
        a
    }
    /// best k-vector at β (max |k|, then lexicographically largest)
    fn best_k(&self, beta: &[i64], kv: &[Vec<i64>]) -> Vec<i64> {
        let mut best = vec![0i64; self.polys.len()];
        let mut bs = 0i64;
        let mut need = vec![0i64; self.nvars];
        for k in kv {
            need.iter_mut().for_each(|x| *x = 0);
            for (j, p) in self.polys.iter().enumerate() {
                for &(v, d) in p {
                    need[v] += k[j] * d;
                }
            }
            if need.iter().zip(beta).all(|(n, b)| n <= b) {
                let s: i64 = k.iter().sum();
                if s > bs {
                    bs = s;
                    best = k.clone();
                }
            }
        }
        best
    }
    /// Pointwise ceiling for every M, every ℓ and every sublattice (top-corner
    /// leading monomials): success needs some j with λ·a_j < 1.
    fn ceiling(&self) -> Option<f64> {
        let mut best: Option<f64> = None;
        for p in &self.polys {
            if p.iter().any(|&(v, d)| self.large[v] && d > 0) {
                continue;
            }
            let s: i64 = p.iter().map(|&(_, d)| d).sum();
            let c = 1.0 / s as f64;
            best = Some(best.map_or(c, |b: f64| b.max(c)));
        }
        best
    }
}

fn gcd(a: i128, b: i128) -> i128 {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

fn multisets(len: usize, max: i64) -> Vec<Vec<i64>> {
    // non-increasing tuples of length len with entries in [0, max]
    let mut out = Vec::new();
    let mut cur = Vec::with_capacity(len);
    fn rec(len: usize, hi: i64, cur: &mut Vec<i64>, out: &mut Vec<Vec<i64>>) {
        if cur.len() == len {
            out.push(cur.clone());
            return;
        }
        for v in (0..=hi).rev() {
            cur.push(v);
            rec(len, v, cur, out);
            cur.pop();
        }
    }
    rec(len, max, &mut cur, &mut out);
    out
}

fn orbit_size(t: &[i64]) -> i64 {
    let mut f = 1i64;
    for i in 1..=t.len() as i64 {
        f *= i;
    }
    let mut run = 1i64;
    for i in 1..t.len() {
        if t[i] == t[i - 1] {
            run += 1;
            f /= run;
        } else {
            run = 1;
        }
    }
    f
}

#[derive(Clone, Debug)]
struct LpResult {
    num: i128,
    den: i128,
    dim: i64,
    sum_k: i64,
    group_max: Vec<i64>,
    is_box: bool,
    bounded: bool,
    nodes: usize,
    iters: usize,
    /// full monomial list of the optimal M with its k-vector (only when asked)
    monos: Vec<(Vec<i64>, Vec<i64>)>,
}

/// Exact optimum over down-closed M ⊂ Π_g [0, ext_g]^{|g|}.
fn solve(sys: &Sys, ell: i64, ext: &[i64], want_monos: bool) -> LpResult {
    let kv = sys.kvecs(ell);
    let gl: Vec<Vec<Vec<i64>>> = sys
        .groups
        .iter()
        .zip(ext)
        .map(|(g, &e)| multisets(g.len(), e))
        .collect();
    let gidx: Vec<HashMap<Vec<i64>, u32>> = gl
        .iter()
        .map(|l| l.iter().enumerate().map(|(i, t)| (t.clone(), i as u32)).collect())
        .collect();
    let ng = sys.groups.len();
    let mut stride = vec![1usize; ng];
    for g in 1..ng {
        stride[g] = stride[g - 1] * gl[g - 1].len();
    }
    let n = stride[ng - 1] * gl[ng - 1].len();
    let decode = |mut i: usize| -> Vec<usize> {
        let mut r = vec![0usize; ng];
        for g in 0..ng {
            r[g] = i % gl[g].len();
            i /= gl[g].len();
        }
        r
    };
    // per node: orbit, k, small sum, large sum
    let mut orb = vec![0i64; n];
    let mut kk = vec![0i64; n];
    let mut ss = vec![0i64; n];
    let mut ls = vec![0i64; n];
    let mut beta = vec![0i64; sys.nvars];
    for i in 0..n {
        let r = decode(i);
        let mut o = 1i64;
        for g in 0..ng {
            let t = &gl[g][r[g]];
            o *= orbit_size(t);
            for (pos, &v) in sys.groups[g].iter().enumerate() {
                beta[v] = t[pos];
            }
        }
        orb[i] = o;
        let k = sys.best_k(&beta, &kv);
        kk[i] = k.iter().sum();
        for v in 0..sys.nvars {
            if sys.large[v] {
                ls[i] += beta[v];
            } else {
                ss[i] += beta[v];
            }
        }
    }
    // down edges
    let mut down: Vec<(u32, u32)> = Vec::new();
    for i in 0..n {
        let r = decode(i);
        for g in 0..ng {
            let t = &gl[g][r[g]];
            for pos in 0..t.len() {
                if t[pos] == 0 || (pos + 1 < t.len() && t[pos + 1] == t[pos]) {
                    continue;
                }
                let mut t2 = t.clone();
                t2[pos] -= 1;
                let j = i - r[g] * stride[g] + gidx[g][&t2] as usize * stride[g];
                down.push((i as u32, j as u32));
            }
        }
    }
    const SC: i128 = 1 << 30;
    let (mut pn, mut pd): (i128, i128) = (0, 1); // current δ = pn/pd
    let mut best: Option<(Vec<bool>, i128, i128)> = None;
    let mut iters = 0;
    loop {
        iters += 1;
        let s = n;
        let t = n + 1;
        let mut edges: Vec<(u32, u32, i64)> = Vec::with_capacity(down.len() + n);
        for &(u, v) in &down {
            edges.push((u, v, INF));
        }
        let mut posw: i128 = 0;
        for i in 0..n {
            let w = orb[i] as i128 * ((kk[i] - ls[i]) as i128 * pd - ss[i] as i128 * pn);
            assert!(w.abs() < (INF as i128) / 4, "weight overflow");
            if w > 0 {
                edges.push((s as u32, i as u32, w as i64));
                posw += w;
            } else if w < 0 {
                edges.push((i as u32, t as u32, (-w) as i64));
            }
        }
        assert!(posw < INF as i128);
        let mut d = Dinic::new(n + 2, &edges);
        let flow = d.max_flow(s, t) as i128;
        let gain = posw - flow; // max closure weight
        if gain <= 0 {
            break;
        }
        let side = d.source_side(s);
        let sel: Vec<bool> = side[..n].to_vec();
        let mut a: i128 = 0;
        let mut b: i128 = 0;
        for i in 0..n {
            if sel[i] {
                a += orb[i] as i128 * (kk[i] - ls[i]) as i128;
                b += orb[i] as i128 * ss[i] as i128;
            }
        }
        if b == 0 {
            // infinite ratio would mean success at every δ; cannot happen for p | stuff
            panic!("degenerate closure");
        }
        let g = gcd(a, b);
        let (a, b) = (a / g, b / g);
        // ratio must strictly beat previous best
        if let Some((_, ba, bb)) = &best {
            assert!(a * bb > ba * b, "Dinkelbach did not improve");
        }
        best = Some((sel, a, b));
        // next δ: just above a/b on the 2^-30 grid
        pn = (a * SC) / b + 1;
        pd = SC;
        if iters > 60 {
            panic!("no convergence");
        }
    }
    let Some((sel, a, b)) = best else {
        return LpResult {
            num: 0,
            den: 1,
            dim: 0,
            sum_k: 0,
            group_max: vec![0; ng],
            is_box: false,
            bounded: false,
            nodes: n,
            iters,
            monos: vec![],
        };
    };
    let mut dim = 0i64;
    let mut sk = 0i64;
    let mut gmax = vec![-1i64; ng];
    let mut monos = Vec::new();
    for i in 0..n {
        if !sel[i] {
            continue;
        }
        dim += orb[i];
        sk += orb[i] * kk[i];
        let r = decode(i);
        for g in 0..ng {
            gmax[g] = gmax[g].max(gl[g][r[g]][0]);
        }
        if want_monos {
            // expand the orbit: all distinct permutations within each group
            let mut partial: Vec<Vec<i64>> = vec![vec![0i64; sys.nvars]];
            for g in 0..ng {
                let t = &gl[g][r[g]];
                let perms = distinct_perms(t);
                let mut next = Vec::new();
                for b0 in &partial {
                    for pm in &perms {
                        let mut b1 = b0.clone();
                        for (pos, &v) in sys.groups[g].iter().enumerate() {
                            b1[v] = pm[pos];
                        }
                        next.push(b1);
                    }
                }
                partial = next;
            }
            for b1 in partial {
                let k = sys.best_k(&b1, &kv);
                monos.push((b1, k));
            }
        }
    }
    let boxsize: i64 = sys
        .groups
        .iter()
        .zip(&gmax)
        .map(|(g, &m)| (m + 1).pow(g.len() as u32))
        .product();
    let bounded = gmax.iter().zip(ext).any(|(m, e)| m == e);
    LpResult {
        num: a,
        den: b,
        dim,
        sum_k: sk,
        group_max: gmax,
        is_box: boxsize == dim,
        bounded,
        nodes: n,
        iters,
        monos,
    }
}

fn distinct_perms(t: &[i64]) -> Vec<Vec<i64>> {
    let mut v = t.to_vec();
    v.sort();
    let mut out = vec![v.clone()];
    // next_permutation
    loop {
        let n = v.len();
        if n < 2 {
            return out;
        }
        let mut i = n - 1;
        while i > 0 && v[i - 1] >= v[i] {
            i -= 1;
        }
        if i == 0 {
            return out;
        }
        let mut j = n - 1;
        while v[j] <= v[i - 1] {
            j -= 1;
        }
        v.swap(i - 1, j);
        v[i..].reverse();
        out.push(v.clone());
    }
}

/// grow the bounding box until the optimum no longer touches it (or a node cap)
fn solve_grow(sys: &Sys, ell: i64, ext0: &[i64], node_cap: usize, want: bool) -> (LpResult, Vec<i64>) {
    let mut ext = ext0.to_vec();
    loop {
        let r = solve(sys, ell, &ext, want);
        if !r.bounded {
            return (r, ext);
        }
        let mut next = ext.clone();
        for (g, e) in next.iter_mut().enumerate() {
            if r.group_max[g] == ext[g] {
                *e += 2;
            }
        }
        // node estimate
        let est: f64 = sys
            .groups
            .iter()
            .zip(&next)
            .map(|(g, &e)| {
                let mut c = 1f64;
                for i in 0..g.len() {
                    c = c * (e as f64 + 1.0 + i as f64) / (i as f64 + 1.0);
                }
                c
            })
            .product();
        if est > node_cap as f64 {
            return (r, ext);
        }
        ext = next;
    }
}

fn s3(nv: usize, a: usize, b: usize) -> Vec<(usize, i64)> {
    let _ = nv;
    vec![(a, 2), (b, 2)]
}

fn part_lp(max_ell_s3: i64, only: Option<String>, ells_override: Option<Vec<i64>>, out: &mut Vec<Value>) {
    println!("## LP: exact down-closed optimum, top-corner leading monomials\n");
    println!("| system | n small | n large | r | ℓ | δ_ℓ exact | δ_ℓ | dim | Σk | group max | box? | bounded? | nodes | iters | ceiling | ms |");
    println!("|:--|--:|--:|--:|--:|:--|--:|--:|--:|:--|:--|:--|--:|--:|--:|--:|");
    let mut cases: Vec<(Sys, Vec<i64>, Vec<i64>)> = Vec::new(); // sys, ells, ext per unit ℓ (+const)
    let s3single = Sys {
        name: "C1 S3".into(),
        nvars: 2,
        large: vec![false; 2],
        groups: vec![vec![0, 1]],
        polys: vec![s3(2, 0, 1)],
    };
    cases.push((s3single, (1..=max_ell_s3).collect(), vec![]));
    let s4single = Sys {
        name: "C1 S4".into(),
        nvars: 3,
        large: vec![false; 3],
        groups: vec![vec![0, 1, 2]],
        polys: vec![vec![(0, 4), (1, 4), (2, 4)]],
    };
    cases.push((s4single, vec![1, 2], vec![]));
    let n1 = Sys {
        name: "N1 S3 x S3 disjoint".into(),
        nvars: 4,
        large: vec![false; 4],
        groups: vec![vec![0, 1], vec![2, 3]],
        polys: vec![s3(4, 0, 1), s3(4, 2, 3)],
    };
    cases.push((n1, (1..=6).collect(), vec![]));
    let p3 = Sys {
        name: "V3 S3 pair s=1".into(),
        nvars: 3,
        large: vec![false; 3],
        groups: vec![vec![0], vec![1, 2]],
        polys: vec![s3(3, 0, 1), s3(3, 0, 2)],
    };
    cases.push((p3, (1..=max_ell_s3).collect(), vec![]));
    let t3 = Sys {
        name: "V3 S3 triple s=1".into(),
        nvars: 4,
        large: vec![false; 4],
        groups: vec![vec![0], vec![1, 2, 3]],
        polys: vec![s3(4, 0, 1), s3(4, 0, 2), s3(4, 0, 3)],
    };
    cases.push((t3, (1..=6).collect(), vec![]));
    let p4s1 = Sys {
        name: "V3 S4 pair s=1".into(),
        nvars: 5,
        large: vec![false; 5],
        groups: vec![vec![0], vec![1, 2], vec![3, 4]],
        polys: vec![vec![(0, 4), (1, 4), (2, 4)], vec![(0, 4), (3, 4), (4, 4)]],
    };
    cases.push((p4s1, vec![1, 2], vec![]));
    let p4s2 = Sys {
        name: "V3 S4 pair s=2".into(),
        nvars: 4,
        large: vec![false; 4],
        groups: vec![vec![0, 1], vec![2, 3]],
        polys: vec![vec![(0, 4), (1, 4), (2, 4)], vec![(0, 4), (1, 4), (3, 4)]],
    };
    cases.push((p4s2, vec![1, 2, 3], vec![]));
    let ch4 = Sys {
        name: "V2 S4 chain (t large)".into(),
        nvars: 4,
        large: vec![false, false, false, true],
        groups: vec![vec![0, 1], vec![2], vec![3]],
        polys: vec![vec![(0, 2), (1, 2), (3, 2)], vec![(2, 2), (3, 2)]],
    };
    cases.push((ch4, vec![1, 2, 3], vec![]));
    for (sys, ells, _) in cases {
        if let Some(o) = &only {
            if !sys.name.contains(o.as_str()) {
                continue;
            }
        }
        let ells = ells_override.clone().unwrap_or(ells);
        for ell in ells {
            let per: Vec<i64> = sys
                .groups
                .iter()
                .map(|g| {
                    let dmax = sys
                        .polys
                        .iter()
                        .flat_map(|p| p.iter())
                        .filter(|(v, _)| g.contains(v))
                        .map(|&(_, d)| d)
                        .max()
                        .unwrap_or(1);
                    let nl = sys.polys.len() as i64;
                    dmax * ell * nl.min(2) + 2
                })
                .collect();
            let cap = if sys.nvars >= 5 { 1_300_000 } else { 2_000_000 };
            let t0 = Instant::now();
            let (r, ext) = solve_grow(&sys, ell, &per, cap, false);
            let ms = t0.elapsed().as_millis();
            let nsmall = sys.large.iter().filter(|l| !**l).count();
            let ceil = sys.ceiling();
            println!(
                "| {} | {} | {} | {} | {} | {}/{} | {:.6} | {} | {} | {:?} (ext {:?}) | {} | {} | {} | {} | {} | {} |",
                sys.name,
                nsmall,
                sys.nvars - nsmall,
                sys.polys.len(),
                ell,
                r.num,
                r.den,
                r.num as f64 / r.den as f64,
                r.dim,
                r.sum_k,
                r.group_max,
                ext,
                if r.is_box { "box" } else { "non-box" },
                if r.bounded { "BINDS (lower bound)" } else { "no" },
                r.nodes,
                r.iters,
                ceil.map_or("none (0)".to_string(), |c| format!("{c:.6}")),
                ms
            );
            out.push(json!({"part":"lp","system":sys.name,"n_small":nsmall,"n_large":sys.nvars-nsmall,
                "r":sys.polys.len(),"ell":ell,"num":r.num.to_string(),"den":r.den.to_string(),
                "delta":r.num as f64/r.den as f64,"dim":r.dim,"sum_k":r.sum_k,"group_max":r.group_max,
                "ext":ext,"is_box":r.is_box,"bounded":r.bounded,"nodes":r.nodes,"iters":r.iters,
                "ceiling":ceil,"ms":ms}));
        }
    }
    println!();
}

// ───────────────────────────── LLL validation ────────────────────────────────

fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}
fn powmod(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64;
    a %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, p);
        }
        a = mulmod(a, a, p);
        e >>= 1;
    }
    r
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
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'w: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'w;
            }
        }
        return false;
    }
    true
}
fn inv(a: u64, p: u64) -> u64 {
    powmod(a, p - 2, p)
}
type Pt = Option<(u64, u64)>;
fn add(p: u64, a: u64, p1: Pt, p2: Pt) -> Pt {
    let (x1, y1) = match p1 {
        None => return p2,
        Some(v) => v,
    };
    let (x2, y2) = match p2 {
        None => return p1,
        Some(v) => v,
    };
    let lam = if x1 == x2 {
        if (y1 + y2) % p == 0 {
            return None;
        }
        let num = (mulmod(3, mulmod(x1, x1, p), p) + a) % p;
        mulmod(num, inv(mulmod(2, y1, p), p), p)
    } else {
        mulmod((y2 + p - y1) % p, inv((x2 + p - x1) % p, p), p)
    };
    let x3 = (mulmod(lam, lam, p) + 2 * p - x1 - x2) % p;
    let y3 = (mulmod(lam, (x1 + p - x3) % p, p) + p - y1) % p;
    Some((x3, y3))
}
fn small_point(p: u64, a: u64, b: u64, bound: u64, rng: &mut StdRng) -> (u64, u64) {
    loop {
        let x = rng.gen_range(0..bound.max(1));
        let r = (mulmod(mulmod(x, x, p), x, p) + mulmod(a, x, p) + b) % p;
        if r == 0 {
            continue;
        }
        if powmod(r, (p - 1) / 2, p) == 1 {
            let y = powmod(r, (p + 1) / 4, p);
            return (x, y);
        }
    }
}

type Poly = HashMap<Vec<i64>, BigInt>;

/// S₃(x_i, x_j, z) over Z with coefficients reduced mod p into [0,p)
fn s3_poly(nv: usize, i: usize, j: usize, z: u64, p: u64, a: u64, b: u64) -> Poly {
    // S3 = (x1−x2)² z² − 2((x1+x2)(x1x2+a) + 2b) z + (x1x2 − a)² − 4b(x1+x2)
    let mut c: HashMap<(i64, i64), i128> = HashMap::new();
    let (z, a, b) = (z as i128, a as i128, b as i128);
    let mut addc = |e: (i64, i64), v: i128| *c.entry(e).or_insert(0) += v;
    // (x1² − 2x1x2 + x2²) z²
    addc((2, 0), z * z);
    addc((1, 1), -2 * z * z);
    addc((0, 2), z * z);
    // −2z( x1²x2 + x1x2² + a x1 + a x2 + 2b )
    addc((2, 1), -2 * z);
    addc((1, 2), -2 * z);
    addc((1, 0), -2 * z * a);
    addc((0, 1), -2 * z * a);
    addc((0, 0), -4 * z * b);
    // x1²x2² − 2a x1x2 + a²
    addc((2, 2), 1);
    addc((1, 1), -2 * a);
    addc((0, 0), a * a);
    // −4b x1 − 4b x2
    addc((1, 0), -4 * b);
    addc((0, 1), -4 * b);
    let pp = p as i128;
    let mut out = Poly::new();
    for ((e1, e2), v) in c {
        let r = ((v % pp) + pp) % pp;
        if r != 0 {
            let mut e = vec![0i64; nv];
            e[i] = e1;
            e[j] = e2;
            out.insert(e, BigInt::from(r));
        }
    }
    out
}

fn pmul(a: &Poly, b: &Poly) -> Poly {
    let mut o = Poly::new();
    for (ea, ca) in a {
        for (eb, cb) in b {
            let e: Vec<i64> = ea.iter().zip(eb).map(|(x, y)| x + y).collect();
            *o.entry(e).or_insert_with(BigInt::zero) += ca * cb;
        }
    }
    o.retain(|_, c| !c.is_zero());
    o
}

fn peval(f: &Poly, root: &[BigInt]) -> BigInt {
    let mut acc = BigInt::zero();
    for (e, c) in f {
        let mut t = c.clone();
        for (i, &x) in e.iter().enumerate() {
            if x > 0 {
                t *= root[i].pow(x as u32);
            }
        }
        acc += t;
    }
    acc
}

struct Trial {
    vanish: usize,
    dim: usize,
    ms: f64,
    recovered: bool,
}

/// One planted instance. `pair`: shared-target system on (x1,x2,x3); else single S₃ on (x1,x2).
fn lll_trial(
    pair: bool,
    ell: i64,
    monos: &[(Vec<i64>, Vec<i64>)],
    p: u64,
    a: u64,
    b: u64,
    bound: u64,
    rng: &mut StdRng,
) -> Trial {
    let nv = if pair { 3 } else { 2 };
    let p1 = small_point(p, a, b, bound, rng);
    let p2 = small_point(p, a, b, bound, rng);
    let p3 = small_point(p, a, b, bound, rng);
    let r1 = add(p, a, Some(p1), Some(p2));
    let r2 = add(p, a, Some(p1), Some(p3));
    let (Some((z1, _)), Some((z2, _))) = (r1, r2) else {
        return Trial { vanish: 0, dim: 0, ms: 0.0, recovered: false };
    };
    let mut polys = vec![s3_poly(nv, 0, 1, z1, p, a, b)];
    let mut root = vec![BigInt::from(p1.0), BigInt::from(p2.0)];
    if pair {
        polys.push(s3_poly(nv, 0, 2, z2, p, a, b));
        root.push(BigInt::from(p3.0));
    }
    let pz = BigInt::from(p);
    for f in &polys {
        assert!((peval(f, &root) % &pz).is_zero(), "planted root not a root");
    }
    // powers
    let mut pw: Vec<Vec<Poly>> = Vec::new();
    for f in &polys {
        let mut v = vec![Poly::from([(vec![0i64; nv], BigInt::one())])];
        for k in 1..=ell as usize {
            let nx = pmul(&v[k - 1], f);
            v.push(nx);
        }
        pw.push(v);
    }
    let corners: Vec<Vec<i64>> = if pair {
        vec![vec![2, 2, 0], vec![2, 0, 2]]
    } else {
        vec![vec![2, 2]]
    };
    let col: HashMap<Vec<i64>, usize> =
        monos.iter().enumerate().map(|(i, (m, _))| (m.clone(), i)).collect();
    let dim = monos.len();
    let bb = BigInt::from(bound.max(1));
    let mut basis: Vec<Vec<BigInt>> = Vec::with_capacity(dim);
    for (beta, k) in monos {
        let mut gamma = beta.clone();
        for (j, &kj) in k.iter().enumerate() {
            for v in 0..nv {
                gamma[v] -= kj * corners[j][v];
            }
        }
        assert!(gamma.iter().all(|&g| g >= 0));
        let ks: i64 = k.iter().sum();
        let mut poly = Poly::from([(gamma.clone(), pz.pow((ell - ks) as u32))]);
        for (j, &kj) in k.iter().enumerate() {
            poly = pmul(&poly, &pw[j][kj as usize]);
        }
        let mut row = vec![BigInt::zero(); dim];
        for (e, c) in &poly {
            let ci = *col.get(e).expect("shift support outside M");
            let tot: i64 = e.iter().sum();
            row[ci] = c * bb.pow(tot as u32);
        }
        basis.push(row);
    }
    let t0 = Instant::now();
    let ok = if dim <= 30 {
        lll_reduce(&mut basis, 0.99)
    } else {
        lll_reduce_hp(&mut basis, 0.99)
    };
    let ms = t0.elapsed().as_secs_f64() * 1e3;
    if ok.is_err() {
        return Trial { vanish: 0, dim, ms, recovered: false };
    }
    let mut vanish = 0;
    for row in &basis {
        let mut pp = Poly::new();
        for (ci, c) in row.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let e = &monos[ci].0;
            let tot: i64 = e.iter().sum();
            pp.insert(e.clone(), c / bb.pow(tot as u32));
        }
        if !pp.is_empty() && peval(&pp, &root).is_zero() {
            vanish += 1;
        }
    }
    Trial { vanish, dim, ms, recovered: vanish >= nv }
}

fn part_lll(bits: u32, ells: &[i64], trials: usize, deltas: &[f64], out: &mut Vec<Value>) {
    let mut p = (1u64 << bits) - 1;
    while !(p % 4 == 3 && is_prime(p)) {
        p -= 1;
    }
    let a = p - 3;
    let b = 31u64;
    println!("## LLL validation: p = {p} ({bits}-bit), a = -3, b = {b}, δ_LLL = 0.99\n");
    println!("Success = at least n reduced rows vanish over Z at the planted root (n = 2 single, 3 pair).\n");
    println!("| system | ℓ | M (closure-optimal) | dim | LP δ_ℓ | δ | B = ⌊p^δ⌋ | success | mean vanishing rows | mean LLL ms |");
    println!("|:--|--:|:--|--:|--:|--:|--:|:--|--:|--:|");
    let single = Sys {
        name: "C1 S3".into(),
        nvars: 2,
        large: vec![false; 2],
        groups: vec![vec![0, 1]],
        polys: vec![s3(2, 0, 1)],
    };
    let pair = Sys {
        name: "V3 S3 pair s=1".into(),
        nvars: 3,
        large: vec![false; 3],
        groups: vec![vec![0], vec![1, 2]],
        polys: vec![s3(3, 0, 1), s3(3, 0, 2)],
    };
    for &ell in ells {
        for (is_pair, sys) in [(false, &single), (true, &pair)] {
            let per: Vec<i64> = sys.groups.iter().map(|_| 4 * ell + 2).collect();
            let (r, _) = solve_grow(sys, ell, &per, 2_000_000, true);
            let lpd = r.num as f64 / r.den as f64;
            let mut rng = StdRng::seed_from_u64(0x05_1a77 + ell as u64 * 1000 + is_pair as u64);
            for &delta in deltas {
                let bound = (p as f64).powf(delta).floor() as u64;
                let mut succ = 0;
                let mut vs = 0usize;
                let mut msum = 0.0;
                let mut dim = 0;
                let mut done = 0;
                for _ in 0..trials {
                    let t = lll_trial(is_pair, ell, &r.monos, p, a, b, bound, &mut rng);
                    if t.dim == 0 {
                        continue;
                    }
                    done += 1;
                    dim = t.dim;
                    succ += t.recovered as usize;
                    vs += t.vanish;
                    msum += t.ms;
                }
                let dn = done.max(1) as f64;
                println!(
                    "| {} | {} | max {:?} | {} | {:.5} | {:.3} | {} | {}/{} | {:.2} | {:.0} |",
                    sys.name,
                    ell,
                    r.group_max,
                    dim,
                    lpd,
                    delta,
                    bound,
                    succ,
                    done,
                    vs as f64 / dn,
                    msum / dn
                );
                out.push(json!({"part":"lll","system":sys.name,"ell":ell,"bits":bits,"p":p.to_string(),
                    "dim":dim,"lp_delta":lpd,"delta":delta,"bound":bound,"success":succ,"trials":done,
                    "mean_vanish":vs as f64/dn,"mean_ms":msum/dn}));
            }
        }
    }
    println!();
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).cloned();
    let part = get("--part").unwrap_or_else(|| "lp".into());
    let mut out = Vec::new();
    let t0 = Instant::now();
    println!("# IC next-10 item 05: multi-polynomial / shared-target Coppersmith reach\n");
    if part == "lp" || part == "all" {
        let me: i64 = get("--max-ell").and_then(|s| s.parse().ok()).unwrap_or(8);
        let only = get("--only");
        let eo: Option<Vec<i64>> = get("--ells").map(|s| s.split(',').map(|x| x.parse().unwrap()).collect());
        part_lp(me, only, eo, &mut out);
    }
    if part == "lll" || part == "all" {
        let bits: u32 = get("--bits").and_then(|s| s.parse().ok()).unwrap_or(40);
        let ells: Vec<i64> = get("--ell")
            .unwrap_or_else(|| "1".into())
            .split(',')
            .map(|s| s.parse().unwrap())
            .collect();
        let trials: usize = get("--trials").and_then(|s| s.parse().ok()).unwrap_or(8);
        let deltas: Vec<f64> = get("--deltas")
            .unwrap_or_else(|| "0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24".into())
            .split(',')
            .map(|s| s.parse().unwrap())
            .collect();
        part_lll(bits, &ells, trials, &deltas, &mut out);
    }
    println!("total wall {:.1} s", t0.elapsed().as_secs_f64());
    if let Some(path) = get("--json") {
        std::fs::write(&path, serde_json::to_string_pretty(&json!({"schema":"ic-next10-05/v1","part":part,"records":out})).unwrap()).unwrap();
        eprintln!("wrote {path}");
    }
}
