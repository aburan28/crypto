//! **Factor-base designs: the solving-degree slope per design.**
//!
//! Runs the arms preregistered in
//! `research/factor_base_designs_20261015/PREREGISTRATION.md` §3 on
//! `p₁ = 786433`, `m ∈ {2, 3}`, `t` rising from 2 until a limit stops a
//! system, two planted targets per cell, and scores each (arm, `m`) series by
//! the least-squares slope `σ = ΔD/Δb` of `D = solving_degree_max` (the
//! repository's degree-bounded F4, `f4_fp::f4`, grevlex, exactly as
//! `examples/pkm_tower_pilot.rs` records it) against `b = log₂|V|`.
//!
//! Arms, each one presentation of the per-summand set `V`:
//!
//! - `A0`: the Kummer tower of the pilot, `x = g·y_0`, `y_{j+1} = y_j²`,
//!   `y_{t−1}² = 1` (`t` unknowns per summand).
//! - `N`: as `A0`, with every `S₃` of the summation given random
//!   coefficients on its own support and a planted zero.
//! - `A1z`, `A1g`: the exponent interval `x = h·∏_j(1 + β_j b_j)`,
//!   `β_j = g^{2^j} − 1`, bits `b_j² = b_j` and partial products
//!   `π_1 = h(1 + β_0 b_0)`, `π_{j+1} = π_j(1 + β_j b_j)`, `x` substituted
//!   into `S₃` (`2t − 1` unknowns). `A1z`: `g = ζ` of order `2^t`, the set of
//!   `A0`. `A1g`: `ord(g)` not a power of 2.
//! - `A4`: the parallel tower `x = h·y_0 + h'·z_0`, Kummer towers of depths
//!   `⌈t/2⌉` on `y` and `⌊t/2⌋` on `z` (`t + 1` unknowns, `x` one of them);
//!   `|V|` is counted exactly and `b = log₂|V|`.
//! - `A11`: the norm chain. A function `f = y − P(X)`, `deg P = m − 1`,
//!   through `−R`; its norm `N(X) = X³ + AX + B − P(X)²` is `(X − x_R)·G(X)`
//!   with `deg G = m`, and `G(gX)` must divide `X^{2^t} − 1`: `t` squaring
//!   steps modulo the monic `G_s(X) = G(gX)/(lead·g^m)`.
//!
//! Every set, every presentation and every planted solution is checked before
//! anything is measured; a failed check panics and is never a degree.
//!
//! ```bash
//! cargo build --release --example fb_designs_slope
//! B=target/release/examples/fb_designs_slope
//! # The study: one child process per system, under `ulimit -v` and a budget.
//! $B run --arms A0,N,A1z,A1g,A4,A11 --m 2,3 --out-dir research/factor_base_designs_20261015/runs
//! # One system, its row on stdout.
//! $B one --arm A1z --m 2 --t 3 --target 0
//! # The instrument check: A0 rebuilt with the pilot's seed rule, against its rows.
//! $B repro --pilot research/pkm_tower_pilot_20260924/runs --m2-max-t 6 --m3-max-t 3
//! # σ per (arm, m) and the preregistered verdict.
//! $B analyze research/factor_base_designs_20261015/runs --pilot research/pkm_tower_pilot_20260924/runs
//! ```
//!
//! Amendment 2 (F1–F4) adds the prime index to the seed (`p₁ = 0`, so every
//! `p₁` seed is unchanged; `p₂ = 2013265921` is 1), the arm `N11` (A11 with
//! the equations defining `s_i` given random coefficients on their support
//! and a planted zero) and the norm chain at `m = 4`
//! (`f = (X + q_0)·y − P(X)`, `deg P = 2`, quartic `G`):
//!
//! ```bash
//! R=research/factor_base_designs_20261015/runs
//! $B check-p2                                  # F1's instrument check at p₂
//! $B regress --rows $R/A11-m3.jsonl --arm A11 --m 3 --t 2   # p₁ rows unchanged
//! $B run --prime p2 --arms A0,A11 --m 2,3 --file 'F1-{arm}-m{m}-p2' --out-dir $R
//! $B run --arms N11 --m 2,3 --file 'F2-{arm}-m{m}' --out-dir $R
//! $B run --arms A0,A11 --m 4 --file 'F4-{arm}-m{m}' --out-dir $R
//! $B run --arms A0,A11 --m 3 --t-min 5 --t-max 5 --budget 3600 --grace 300 --file 'F3-{arm}-m3-t5' --out-dir $R
//! $B analyze $R --amendment2                   # every Amendment 2 decision rule
//! $B compare-rows A.jsonl B.jsonl               # two runs of the same systems, row by row
//! ```
//!
//! A stage diagnostic of the solver axis (`AGENTS.md` §2): it prices no
//! pipeline and reports no `S`.

use std::collections::{BTreeMap, HashMap};
use std::fs::{File, OpenOptions};
use std::io::{Read, Write};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use crypto_lib::cryptanalysis::f4_fp::{self, F4Options, Ordering, Poly};
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps, PrimeCurve, PrimePoint};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

/// `p₁ = 3·2^18 + 1`, the pilot's prime.
const P1: u64 = 786_433;
/// `p₂ = 15·2^27 + 1`, the replication prime (§5, Amendment 2).
const P2: u64 = 2_013_265_921;

/// A prime by its index in the seed rule (Amendment 2): `p₁ = 0`, `p₂ = 1`.
fn prime_of(idx: usize) -> u64 {
    match idx {
        0 => P1,
        1 => P2,
        _ => panic!("unknown prime index {idx}"),
    }
}

fn prime_label(idx: usize) -> &'static str {
    ["p1", "p2"][idx]
}

fn parse_prime(s: &str) -> usize {
    match s {
        "p1" => 0,
        "p2" => 1,
        _ => panic!("unknown prime `{s}` (p1 or p2)"),
    }
}
/// "FBDS", the preregistered seed base.
const SEED: u64 = 0x4642_4453;
/// The pilot's seed base ("PKMT"), for the instrument check only.
const PILOT_SEED: u64 = 0x504B_4D54;

// ── Field arithmetic (copied from the pilot) ───────────────────────

fn mulm(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

fn addm(a: u64, b: u64, p: u64) -> u64 {
    (a + b) % p
}

fn subm(a: u64, b: u64, p: u64) -> u64 {
    (a + p - b % p) % p
}

fn powm(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64;
    a %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mulm(r, a, p);
        }
        a = mulm(a, a, p);
        e >>= 1;
    }
    r
}

fn invm(a: u64, p: u64) -> u64 {
    assert!(!a.is_multiple_of(p), "inverse of zero");
    powm(a, p - 2, p)
}

/// The multiplicative order of `a` in `F_p^*`.
fn order(a: u64, p: u64) -> u64 {
    let mut n = p - 1;
    let mut primes = Vec::new();
    let mut q = 2;
    let mut r = n;
    while q * q <= r {
        if r.is_multiple_of(q) {
            primes.push(q);
            while r.is_multiple_of(q) {
                r /= q;
            }
        }
        q += 1;
    }
    if r > 1 {
        primes.push(r);
    }
    for q in primes {
        while n.is_multiple_of(q) && powm(a, n / q, p) == 1 {
            n /= q;
        }
    }
    n
}

// ── Polynomials while building a system (copied from the pilot) ────

/// A polynomial over `F_p` in `n` variables, as a map from exponent
/// vectors to non-zero coefficients.
#[derive(Clone, Debug)]
struct Pol {
    n: usize,
    terms: BTreeMap<Vec<u32>, u64>,
}

impl Pol {
    fn zero(n: usize) -> Self {
        Pol {
            n,
            terms: BTreeMap::new(),
        }
    }

    fn constant(n: usize, c: u64, p: u64) -> Self {
        let mut f = Pol::zero(n);
        if !c.is_multiple_of(p) {
            f.terms.insert(vec![0; n], c % p);
        }
        f
    }

    fn var(n: usize, i: usize) -> Self {
        let mut e = vec![0; n];
        e[i] = 1;
        let mut f = Pol::zero(n);
        f.terms.insert(e, 1);
        f
    }

    fn add(&self, o: &Pol, p: u64) -> Pol {
        let mut f = self.clone();
        for (e, c) in &o.terms {
            let v = f.terms.entry(e.clone()).or_insert(0);
            *v = addm(*v, *c, p);
        }
        f.terms.retain(|_, c| *c != 0);
        f
    }

    fn scale(&self, k: u64, p: u64) -> Pol {
        let mut f = Pol::zero(self.n);
        for (e, c) in &self.terms {
            let v = mulm(*c, k, p);
            if v != 0 {
                f.terms.insert(e.clone(), v);
            }
        }
        f
    }

    fn sub(&self, o: &Pol, p: u64) -> Pol {
        self.add(&o.scale(p - 1, p), p)
    }

    fn mul(&self, o: &Pol, p: u64) -> Pol {
        let mut f = Pol::zero(self.n);
        for (e1, c1) in &self.terms {
            for (e2, c2) in &o.terms {
                let e: Vec<u32> = e1.iter().zip(e2).map(|(a, b)| a + b).collect();
                let v = f.terms.entry(e).or_insert(0);
                *v = addm(*v, mulm(*c1, *c2, p), p);
            }
        }
        f.terms.retain(|_, c| *c != 0);
        f
    }

    fn eval(&self, x: &[u64], p: u64) -> u64 {
        let mut acc = 0u64;
        for (e, c) in &self.terms {
            let mut t = *c;
            for (i, k) in e.iter().enumerate() {
                if *k > 0 {
                    t = mulm(t, powm(x[i], *k as u64, p), p);
                }
            }
            acc = addm(acc, t, p);
        }
        acc
    }

    fn degree(&self) -> u32 {
        self.terms
            .keys()
            .map(|e| e.iter().sum::<u32>())
            .max()
            .unwrap_or(0)
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    fn to_poly(&self, p: u64) -> Poly {
        let terms: Vec<(Vec<u32>, u64)> = self.terms.iter().map(|(e, c)| (e.clone(), *c)).collect();
        f4_fp::normalise(&terms, p, Ordering::Grevlex)
    }
}

/// Semaev's third summation polynomial for `y² = x³ + ax + b`, on three
/// polynomial arguments (the pilot's `s3`).
fn s3(x1: &Pol, x2: &Pol, x3: &Pol, a: u64, b: u64, p: u64) -> Pol {
    let n = x1.n;
    let d = x1.sub(x2, p);
    let t1 = d.mul(&d, p).mul(x3, p).mul(x3, p);
    let s = x1.add(x2, p);
    let pr = x1.mul(x2, p);
    let inner = s
        .mul(&pr.add(&Pol::constant(n, a, p), p), p)
        .add(&Pol::constant(n, mulm(2, b, p), p), p);
    let t2 = inner.mul(x3, p).scale(p - 2, p);
    let q = pr.sub(&Pol::constant(n, a, p), p);
    let t3 = q.mul(&q, p).sub(&s.scale(mulm(4, b, p), p), p);
    t1.add(&t2, p).add(&t3, p)
}

// ── Numeric univariate polynomials (coefficients low to high) ──────

fn npmul(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = addm(r[i + j], mulm(x, y, p), p);
        }
    }
    r
}

/// `a = q·d + r` for a monic `d` of degree `k`: `(q, r)` with `r` of length
/// `k`.
fn npdivmod_monic(a: &[u64], d: &[u64], p: u64) -> (Vec<u64>, Vec<u64>) {
    let k = d.len() - 1;
    assert_eq!(d[k], 1, "divisor not monic");
    let mut r = a.to_vec();
    if r.len() < k {
        r.resize(k, 0);
        return (Vec::new(), r);
    }
    let mut q = vec![0u64; r.len() - k];
    for i in (k..r.len()).rev() {
        let c = r[i];
        q[i - k] = c;
        if c != 0 {
            for (j, &dj) in d.iter().enumerate() {
                r[i - k + j] = subm(r[i - k + j], mulm(c, dj, p), p);
            }
        }
    }
    r.truncate(k);
    (q, r)
}

/// `X^{2^t} mod d` for a monic `d`, by `t` squarings from `X`.
fn x_pow2_mod(d: &[u64], t: usize, p: u64) -> Vec<u64> {
    let k = d.len() - 1;
    let mut r = vec![0u64; k];
    if k == 1 {
        r[0] = subm(0, d[0], p);
    } else {
        r[1] = 1;
    }
    for _ in 0..t {
        let sq = npmul(&r, &r, p);
        r = npdivmod_monic(&sq, d, p).1;
    }
    r
}

// ── Arms and sets ──────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Arm {
    A0,
    N,
    A1z,
    A1g,
    A4,
    A11,
    N11,
}

impl Arm {
    const ALL: [Arm; 7] = [
        Arm::A0,
        Arm::N,
        Arm::A1z,
        Arm::A1g,
        Arm::A4,
        Arm::A11,
        Arm::N11,
    ];

    fn name(self) -> &'static str {
        match self {
            Arm::A0 => "A0",
            Arm::N => "N",
            Arm::A1z => "A1z",
            Arm::A1g => "A1g",
            Arm::A4 => "A4",
            Arm::A11 => "A11",
            Arm::N11 => "N11",
        }
    }

    /// The arm's index in the seed rule.
    fn id(self) -> u64 {
        self as u64
    }

    fn parse(s: &str) -> Arm {
        Arm::ALL
            .into_iter()
            .find(|a| a.name() == s)
            .unwrap_or_else(|| panic!("unknown arm `{s}`"))
    }
}

/// How one summand's abscissa is presented.
enum Shape {
    /// `x = g·y_0`, `y_{j+1} = y_j²`, `y_{t−1}² = 1`.
    Kummer { g: u64 },
    /// `x = h·∏_j(1 + β_j b_j)` with bits and partial products; `k_of` maps
    /// an element `h·gen^k` to `k`.
    Interval {
        h: u64,
        beta: Vec<u64>,
        k_of: HashMap<u64, u64>,
    },
    /// `x = h·y_0 + h'·z_0` with Kummer towers of depths `a` and `b`;
    /// `ij_of` maps an element to the first `(i, j)` with
    /// `x = h·ζ_a^i + h'·ζ_b^j`.
    Parallel {
        h: u64,
        h2: u64,
        za: u64,
        zb: u64,
        a: usize,
        b: usize,
        ij_of: HashMap<u64, (u64, u64)>,
    },
}

/// One summand's factor base and its membership presentation.
struct Design {
    t: usize,
    /// Unknowns per summand.
    per: usize,
    /// The distinct elements of `V`, in construction order.
    v: Vec<u64>,
    shape: Shape,
    info: Value,
}

/// The pilot's Kummer equations on `depth` unknowns from `base`:
/// `y_{j+1} − y_j²`, and `1 − y_{depth−1}²`.
fn kummer_eqs(n: usize, base: usize, depth: usize, p: u64) -> Vec<Pol> {
    (0..depth)
        .map(|j| {
            let yj = Pol::var(n, base + j);
            let next = if j + 1 < depth {
                Pol::var(n, base + j + 1)
            } else {
                Pol::constant(n, 1, p)
            };
            next.sub(&yj.mul(&yj, p), p)
        })
        .collect()
}

fn kummer_values(y0: u64, depth: usize, p: u64) -> Vec<u64> {
    let mut y = vec![y0; depth];
    for j in 1..depth {
        y[j] = mulm(y[j - 1], y[j - 1], p);
    }
    y
}

impl Design {
    /// The abscissa of summand `i` (unknowns from `base`) as a polynomial.
    fn x_expr(&self, n: usize, base: usize, p: u64) -> Pol {
        let t = self.t;
        match &self.shape {
            Shape::Kummer { g } => Pol::var(n, base).scale(*g, p),
            Shape::Interval { beta, .. } => {
                // b_j at base + j, π_k (k = 1 … t−1) at base + t + k − 1.
                let pi_last = Pol::var(n, base + 2 * t - 2);
                let factor =
                    Pol::constant(n, 1, p).add(&Pol::var(n, base + t - 1).scale(beta[t - 1], p), p);
                pi_last.mul(&factor, p)
            }
            Shape::Parallel { .. } => Pol::var(n, base),
        }
    }

    /// The membership equations of one summand.
    fn membership(&self, n: usize, base: usize, p: u64) -> Vec<Pol> {
        let t = self.t;
        match &self.shape {
            Shape::Kummer { .. } => kummer_eqs(n, base, t, p),
            Shape::Interval { h, beta, .. } => {
                let one = Pol::constant(n, 1, p);
                let bit = |j: usize| Pol::var(n, base + j);
                let pi = |k: usize| Pol::var(n, base + t + k - 1);
                let factor = |j: usize| one.add(&bit(j).scale(beta[j], p), p);
                let mut eqs: Vec<Pol> = (0..t)
                    .map(|j| bit(j).mul(&bit(j), p).sub(&bit(j), p))
                    .collect();
                eqs.push(pi(1).sub(&factor(0).scale(*h, p), p));
                for k in 1..t - 1 {
                    eqs.push(pi(k + 1).sub(&pi(k).mul(&factor(k), p), p));
                }
                eqs
            }
            Shape::Parallel { h, h2, a, b, .. } => {
                let x = Pol::var(n, base);
                let y0 = Pol::var(n, base + 1);
                let z0 = Pol::var(n, base + 1 + a);
                let mut eqs = vec![x.sub(&y0.scale(*h, p), p).sub(&z0.scale(*h2, p), p)];
                eqs.extend(kummer_eqs(n, base + 1, *a, p));
                eqs.extend(kummer_eqs(n, base + 1 + a, *b, p));
                eqs
            }
        }
    }

    /// The values of one summand's unknowns at an element of `V`.
    fn values(&self, x: u64, p: u64) -> Vec<u64> {
        let t = self.t;
        match &self.shape {
            Shape::Kummer { g } => kummer_values(mulm(x, invm(*g, p), p), t, p),
            Shape::Interval { h, beta, k_of } => {
                let k = k_of[&x];
                let bits: Vec<u64> = (0..t).map(|j| (k >> j) & 1).collect();
                let factor = |j: usize| addm(1, mulm(beta[j], bits[j], p), p);
                let mut vals = bits.clone();
                let mut pi = mulm(*h, factor(0), p);
                vals.push(pi);
                for k in 1..t - 1 {
                    pi = mulm(pi, factor(k), p);
                    vals.push(pi);
                }
                vals
            }
            Shape::Parallel {
                za,
                zb,
                a,
                b,
                ij_of,
                ..
            } => {
                let (i, j) = ij_of[&x];
                let mut vals = vec![x];
                vals.extend(kummer_values(powm(*za, i, p), *a, p));
                vals.extend(kummer_values(powm(*zb, j, p), *b, p));
                vals
            }
        }
    }

    /// `V` has no repeats, and every element of `V` satisfies the
    /// presentation with its abscissa equal to the element. Panics
    /// otherwise.
    fn check(&self, p: u64) {
        let mut s = self.v.clone();
        s.sort_unstable();
        s.dedup();
        assert_eq!(s.len(), self.v.len(), "V has repeats");
        let n = self.per;
        let eqs = self.membership(n, 0, p);
        let x = self.x_expr(n, 0, p);
        for &e in &self.v {
            let vals = self.values(e, p);
            assert_eq!(vals.len(), n, "wrong number of values");
            for (k, f) in eqs.iter().enumerate() {
                assert_eq!(f.eval(&vals, p), 0, "membership equation {k} fails at {e}");
            }
            assert_eq!(
                x.eval(&vals, p),
                e,
                "the presentation's abscissa is not {e}"
            );
        }
    }

    /// `b = log₂|V|`.
    fn bits(&self) -> f64 {
        (self.v.len() as f64).log2()
    }
}

/// An element of exact order `2^t` in `F_p^*` (the pilot's draw).
fn root_of_unity(p: u64, t: usize, rng: &mut StdRng) -> u64 {
    assert!(
        (p - 1).is_multiple_of(1u64 << t),
        "2^{t} does not divide p − 1"
    );
    loop {
        let h = rng.gen_range(2..p);
        if powm(h, (p - 1) / 2, p) == p - 1 {
            return powm(h, (p - 1) >> t, p);
        }
    }
}

/// `A0`, `N`, `A11`: the pilot's Kummer coset `g·μ_{2^t}`, drawn as the
/// pilot draws it (`ζ`, then `g`).
fn kummer_design(p: u64, t: usize, rng: &mut StdRng) -> Design {
    let zeta = root_of_unity(p, t, rng);
    let g = rng.gen_range(2..p);
    let v: Vec<u64> = (0..1u64 << t)
        .map(|k| mulm(g, powm(zeta, k, p), p))
        .collect();
    Design {
        t,
        per: t,
        v,
        shape: Shape::Kummer { g },
        info: json!({"g": g, "zeta": zeta}),
    }
}

/// `A1z` (`generic = false`): `h·μ_{2^t}` with `gen = ζ`. `A1g`: `h·gen^k`,
/// `0 ≤ k < 2^t`, with `ord(gen)` not a power of 2.
fn interval_design(p: u64, t: usize, generic: bool, rng: &mut StdRng) -> Design {
    assert!(t >= 2, "the interval presentation needs t ≥ 2");
    let two = (p - 1).trailing_zeros() as u64;
    loop {
        let (gen, h) = if generic {
            let h = rng.gen_range(2..p);
            let gen = rng.gen_range(2..p);
            // ord(gen) is a power of 2 exactly when gen^{2^two} = 1.
            if powm(gen, 1u64 << two, p) == 1 {
                continue;
            }
            (gen, h)
        } else {
            let zeta = root_of_unity(p, t, rng);
            (zeta, rng.gen_range(2..p))
        };
        let mut k_of = HashMap::new();
        let mut v = Vec::with_capacity(1 << t);
        let mut e = h;
        for k in 0..1u64 << t {
            if k_of.insert(e, k).is_some() {
                break;
            }
            v.push(e);
            e = mulm(e, gen, p);
        }
        if v.len() != 1 << t {
            continue;
        }
        let beta: Vec<u64> = (0..t)
            .map(|j| subm(powm(gen, 1u64 << j, p), 1, p))
            .collect();
        let ord = order(gen, p);
        let info = json!({"h": h, "gen": gen, "ord_gen": ord, "beta": &beta});
        return Design {
            t,
            per: 2 * t - 1,
            v,
            shape: Shape::Interval { h, beta, k_of },
            info,
        };
    }
}

/// `A4`: `{h·ζ_a^i + h'·ζ_b^j}`, depths `a = ⌈t/2⌉`, `b = ⌊t/2⌋`.
fn parallel_design(p: u64, t: usize, rng: &mut StdRng) -> Design {
    assert!(t >= 2, "the parallel tower needs t ≥ 2");
    let (a, b) = (t.div_ceil(2), t / 2);
    let za = root_of_unity(p, a, rng);
    let zb = root_of_unity(p, b, rng);
    let h = rng.gen_range(2..p);
    let h2 = rng.gen_range(2..p);
    let mut ij_of = HashMap::new();
    let mut v = Vec::new();
    for i in 0..1u64 << a {
        for j in 0..1u64 << b {
            let x = addm(mulm(h, powm(za, i, p), p), mulm(h2, powm(zb, j, p), p), p);
            if let std::collections::hash_map::Entry::Vacant(e) = ij_of.entry(x) {
                e.insert((i, j));
                v.push(x);
            }
        }
    }
    let collisions = (1usize << t) - v.len();
    let info = json!({"h": h, "h2": h2, "zeta_a": za, "zeta_b": zb, "a": a, "b": b,
                      "pairs": 1usize << t, "collisions": collisions});
    Design {
        t,
        per: t + 1,
        v,
        shape: Shape::Parallel {
            h,
            h2,
            za,
            zb,
            a,
            b,
            ij_of,
        },
        info,
    }
}

fn build_design(arm: Arm, p: u64, t: usize, rng: &mut StdRng) -> Design {
    let d = match arm {
        Arm::A0 | Arm::N | Arm::A11 | Arm::N11 => kummer_design(p, t, rng),
        Arm::A1z => interval_design(p, t, false, rng),
        Arm::A1g => interval_design(p, t, true, rng),
        Arm::A4 => parallel_design(p, t, rng),
    };
    d.check(p);
    d
}

// ── Curves and targets (copied from the pilot) ─────────────────────

fn nonsingular(a: u64, b: u64, p: u64) -> bool {
    let d = addm(
        mulm(4, mulm(mulm(a, a, p), a, p), p),
        mulm(27, mulm(b, b, p), p),
        p,
    );
    d != 0
}

/// A target curve on which at least `need` elements of `V` are abscissae.
fn target_curve(p: u64, v: &[u64], need: usize, rng: &mut StdRng) -> (PrimeCurve, Vec<u64>) {
    loop {
        let a = rng.gen_range(0..p);
        let b = rng.gen_range(1..p);
        if !nonsingular(a, b, p) {
            continue;
        }
        let curve = PrimeCurve { p, a, b };
        let on: Vec<u64> = v
            .iter()
            .copied()
            .filter(|&x| curve.legendre(curve.rhs(x)) == 1)
            .collect();
        if on.len() >= need {
            return (curve, on);
        }
    }
}

/// A planted target `R = P_1 + … + P_m` over distinct base abscissae.
struct Target {
    x_r: u64,
    r: PrimePoint,
    pts: Vec<PrimePoint>,
    /// The abscissae of the partial sums `P_1 + … + P_{k+1}`, `k < m − 1`.
    partial: Vec<u64>,
}

/// The pilot's `plant`: `m` distinct abscissae of `on`, random signs; `None`
/// when a partial sum is the identity.
fn plant(curve: &PrimeCurve, m: usize, on: &[u64], rng: &mut StdRng) -> Option<Target> {
    let p = curve.p;
    let mut ops = GroupOps::default();
    let mut pick: Vec<u64> = Vec::with_capacity(m);
    while pick.len() < m {
        let x = on[rng.gen_range(0..on.len())];
        if !pick.contains(&x) {
            pick.push(x);
        }
    }
    let pts: Vec<PrimePoint> = pick
        .iter()
        .map(|&x| {
            let y = curve.sqrt(curve.rhs(x)).expect("abscissa on the curve");
            PrimePoint::affine(x, if rng.gen::<bool>() { y } else { p - y })
        })
        .collect();
    let mut acc = pts[0];
    let mut partial = Vec::new();
    for pt in &pts[1..] {
        acc = curve.add(&mut ops, acc, *pt);
        if acc.infinity {
            return None;
        }
        partial.push(acc.x);
    }
    Some(Target {
        x_r: acc.x,
        r: acc,
        pts,
        partial,
    })
}

// ── Systems ────────────────────────────────────────────────────────

/// `S₃(x_1, x_2, x_R)` at `m = 2`; at `m = 3` the pilot's chain
/// `S₃(x_1, x_2, u_0)`, `S₃(u_0, x_3, x_R)`.
fn summation(d: &Design, curve: &PrimeCurve, m: usize, n: usize, x_r: u64) -> Vec<Pol> {
    let p = curve.p;
    let (a, b) = (curve.a, curve.b);
    let x = |i: usize| d.x_expr(n, i * d.per, p);
    let xr = Pol::constant(n, x_r, p);
    if m == 2 {
        return vec![s3(&x(0), &x(1), &xr, a, b, p)];
    }
    let u = |k: usize| Pol::var(n, m * d.per + k);
    let mut eqs = vec![s3(&x(0), &x(1), &u(0), a, b, p)];
    for k in 1..m - 2 {
        eqs.push(s3(&u(k - 1), &x(k + 1), &u(k), a, b, p));
    }
    eqs.push(s3(&u(m - 3), &x(m - 1), &xr, a, b, p));
    eqs
}

/// A tower-style system (`A0`, `N`, `A1z`, `A1g`, `A4`): the summation, then
/// each summand's membership. `null` replaces every summation polynomial by
/// random coefficients on its support, shifted to vanish at the planted
/// solution.
fn tower_system(
    d: &Design,
    curve: &PrimeCurve,
    m: usize,
    target: &Target,
    null: Option<&mut StdRng>,
) -> (Vec<Pol>, usize, Vec<u64>) {
    let p = curve.p;
    let n = m * d.per + (m - 2);
    let mut sol = Vec::with_capacity(n);
    for pt in &target.pts {
        sol.extend(d.values(pt.x, p));
    }
    sol.extend_from_slice(&target.partial[..m - 2]);
    assert_eq!(sol.len(), n);
    let mut eqs = summation(d, curve, m, n, target.x_r);
    if let Some(rng) = null {
        for f in eqs.iter_mut() {
            let mut g = Pol::zero(n);
            for e in f.terms.keys() {
                g.terms.insert(e.clone(), rng.gen_range(1..p));
            }
            let c = g.eval(&sol, p);
            *f = g.sub(&Pol::constant(n, c, p), p);
        }
    }
    for i in 0..m {
        eqs.extend(d.membership(n, i * d.per, p));
    }
    (eqs, n, sol)
}

/// Gaussian elimination mod `p` on an augmented `k × (k + 1)` system;
/// `None` when it is singular.
fn solve_linear(mut a: Vec<Vec<u64>>, p: u64) -> Option<Vec<u64>> {
    let k = a.len();
    for c in 0..k {
        let r = (c..k).find(|&r| a[r][c] != 0)?;
        a.swap(c, r);
        let inv = invm(a[c][c], p);
        for x in a[c].iter_mut() {
            *x = mulm(*x, inv, p);
        }
        let prow = a[c].clone();
        for (r, row) in a.iter_mut().enumerate() {
            if r != c && row[c] != 0 {
                let f = row[c];
                for (x, &y) in row.iter_mut().zip(&prow) {
                    *x = subm(*x, mulm(f, y, p), p);
                }
            }
        }
    }
    Some(a.into_iter().map(|row| row[k]).collect())
}

/// The Lagrange interpolant through `(xs, ys)`, coefficients low to high.
fn interpolate(xs: &[u64], ys: &[u64], p: u64) -> Vec<u64> {
    let m = xs.len();
    let mut pnum = vec![0u64; m];
    for i in 0..m {
        let mut basis = vec![1u64];
        let mut denom = 1u64;
        for k in 0..m {
            if k != i {
                basis = npmul(&basis, &[subm(0, xs[k], p), 1], p);
                denom = mulm(denom, subm(xs[i], xs[k], p), p);
            }
        }
        let w = mulm(ys[i], invm(denom, p), p);
        for (c, bk) in pnum.iter_mut().zip(&basis) {
            *c = addm(*c, mulm(w, *bk, p), p);
        }
    }
    pnum
}

/// The norm-chain system (`A11`) for `m ∈ {2, 3, 4}` on the coset
/// `g·μ_{2^t}`.
///
/// Unknowns: `c_1 … c_{m−1}` (the function `f = y − P(X)` through `−R`,
/// `P(X) = y(−R) + Σ c_i (X^i − x_R^i)`; `c_1 = λ` at `m = 2`), then
/// `s_0 … s_{m−1}` (`G_s = X^m + Σ s_i X^i`), then per step `j < t` the
/// quotient `q_j` (`m − 1` coefficients) and `r_{j+1}` (`m` coefficients).
/// At `m = 4` (Amendment 2, F4) the function is `f = Q(X)·y − P(X)` with
/// `Q = X + q_0` and `P(X) = Q(x_R)·y(−R) + c_1(X − x_R) + c_2(X² − x_R²)`;
/// the free coefficients are `q_0, c_1, c_2` in that order, and
/// `N(X) = Q²(X³ + AX + B) − P²` has degree 5, so `G` is a quartic.
/// Returns `None` when the planted function has a vanishing leading
/// coefficient, or at `m = 4` when no such `f` with `Q` monic passes through
/// the summands (the caller replants).
fn norm_chain_system(
    g: u64,
    t: usize,
    curve: &PrimeCurve,
    m: usize,
    target: &Target,
) -> Option<(Vec<Pol>, usize, Vec<u64>, Value)> {
    assert!(
        (2..=4).contains(&m),
        "the norm chain is written for m ∈ {{2, 3, 4}}"
    );
    let p = curve.p;
    let (ca, cb) = (curve.a, curve.b);
    let x_r = target.x_r;
    let y0 = subm(0, target.r.y, p); // y(−R)
    let n_c = m - 1;
    let n_s = m;
    let per_step = 2 * m - 1;
    let n = n_c + n_s + per_step * t;
    let cvar = |i: usize| Pol::var(n, i - 1); // c_i, i = 1 … m−1
    let svar = |i: usize| Pol::var(n, n_c + i);
    let step_base = |j: usize| n_c + n_s + per_step * j;

    // P(X) and N(X) = X³ + AX + B − P(X)², as coefficient lists of Pol.
    let (pc, deg_n) = if m == 4 {
        // P(X) = Q(x_R)·y(−R) + c_1(X − x_R) + c_2(X² − x_R²), unknowns
        // q_0 = var 0, c_1 = var 1, c_2 = var 2.
        let q0 = Pol::var(n, 0);
        let (c1, c2) = (Pol::var(n, 1), Pol::var(n, 2));
        let p0 = Pol::constant(n, mulm(x_r, y0, p), p)
            .add(&q0.scale(y0, p), p)
            .sub(&c1.scale(x_r, p), p)
            .sub(&c2.scale(mulm(x_r, x_r, p), p), p);
        (vec![p0, c1, c2], 5)
    } else {
        let mut pc: Vec<Pol> = vec![Pol::constant(n, y0, p)];
        for i in 1..m {
            let xri = powm(x_r, i as u64, p);
            pc[0] = pc[0].sub(&cvar(i).scale(xri, p), p);
            pc.push(cvar(i));
        }
        (pc, 3.max(2 * (m - 1)))
    };
    let mut nc: Vec<Pol> = vec![Pol::zero(n); deg_n + 1];
    if m == 4 {
        // Q² · (X³ + AX + B), Q = X + q_0.
        let q0 = Pol::var(n, 0);
        let q2 = [q0.mul(&q0, p), q0.scale(2, p), Pol::constant(n, 1, p)];
        let cubic = [cb, ca, 0, 1];
        for (i, qi) in q2.iter().enumerate() {
            for (j, &cj) in cubic.iter().enumerate() {
                nc[i + j] = nc[i + j].add(&qi.scale(cj, p), p);
            }
        }
    } else {
        nc[3] = nc[3].add(&Pol::constant(n, 1, p), p);
        nc[1] = nc[1].add(&Pol::constant(n, ca, p), p);
        nc[0] = nc[0].add(&Pol::constant(n, cb, p), p);
    }
    for i in 0..pc.len() {
        for j in 0..pc.len() {
            nc[i + j] = nc[i + j].sub(&pc[i].mul(&pc[j], p), p);
        }
    }
    // Exact division by (X − x_R), verified: the remainder is the zero
    // polynomial in the c's.
    let mut gq: Vec<Pol> = vec![Pol::zero(n); deg_n];
    gq[deg_n - 1] = nc[deg_n].clone();
    for k in (1..deg_n).rev() {
        gq[k - 1] = nc[k].add(&gq[k].scale(x_r, p), p);
    }
    let rem = nc[0].add(&gq[0].scale(x_r, p), p);
    assert!(rem.is_zero(), "N(X) is not divisible by X − x_R: {rem:?}");
    assert_eq!(gq.len(), m + 1, "deg G ≠ m");
    let lead = gq[m].clone();
    let gpow = |e: usize| powm(g, e as u64, p);
    // lead·g^{m−i}·s_i − G_i = 0, i = 0 … m−1 (at m = 2, lead = 1:
    // g·s_1 = e_1, g²·s_0 = e_0).
    let mut eqs: Vec<Pol> = (0..m)
        .map(|i| lead.mul(&svar(i), p).scale(gpow(m - i), p).sub(&gq[i], p))
        .collect();
    // The chain: r_0 = X; r_j² = q_j·G_s + r_{j+1}; r_t = 1.
    let mut gs: Vec<Pol> = (0..m).map(svar).collect();
    gs.push(Pol::constant(n, 1, p));
    let mut r: Vec<Pol> = (0..m)
        .map(|k| Pol::constant(n, u64::from(k == 1), p))
        .collect();
    for j in 0..t {
        let b = step_base(j);
        let q: Vec<Pol> = (0..m - 1).map(|k| Pol::var(n, b + k)).collect();
        let rn: Vec<Pol> = (0..m).map(|k| Pol::var(n, b + m - 1 + k)).collect();
        let mut lhs = vec![Pol::zero(n); 2 * m - 1];
        for (i, ri) in r.iter().enumerate() {
            for (k, rk) in r.iter().enumerate() {
                lhs[i + k] = lhs[i + k].add(&ri.mul(rk, p), p);
            }
        }
        let mut rhs = vec![Pol::zero(n); 2 * m - 1];
        for (i, qi) in q.iter().enumerate() {
            for (k, gk) in gs.iter().enumerate() {
                rhs[i + k] = rhs[i + k].add(&qi.mul(gk, p), p);
            }
        }
        for (k, rk) in rn.iter().enumerate() {
            rhs[k] = rhs[k].add(rk, p);
        }
        for k in (0..2 * m - 1).rev() {
            eqs.push(lhs[k].sub(&rhs[k], p));
        }
        r = rn;
    }
    for (k, rk) in r.iter().enumerate() {
        eqs.push(rk.sub(&Pol::constant(n, u64::from(k == 0), p), p));
    }

    // The planted solution: P through the m summand points, then the c's,
    // the s's and the chain's values.
    let xs: Vec<u64> = target.pts.iter().map(|pt| pt.x).collect();
    let ys: Vec<u64> = target.pts.iter().map(|pt| pt.y).collect();
    let mut sol = vec![0u64; n];
    if m == 4 {
        // (x_i + q_0)·y_i = P(x_i) at the four summands: a linear system in
        // (q_0, P_0, P_1, P_2); singular means no f with Q monic.
        let rows: Vec<Vec<u64>> = (0..4)
            .map(|i| {
                let x2 = mulm(xs[i], xs[i], p);
                let rhs = subm(0, mulm(xs[i], ys[i], p), p);
                vec![ys[i], p - 1, subm(0, xs[i], p), subm(0, x2, p), rhs]
            })
            .collect();
        let v = solve_linear(rows, p)?;
        let (q0, p0, p1, p2) = (v[0], v[1], v[2], v[3]);
        let at_r = addm(
            addm(p0, mulm(p1, x_r, p), p),
            mulm(p2, mulm(x_r, x_r, p), p),
            p,
        );
        assert_eq!(
            mulm(addm(x_r, q0, p), y0, p),
            at_r,
            "the planted function misses −R"
        );
        sol[0] = q0;
        sol[1] = p1;
        sol[2] = p2;
        assert_eq!(pc[0].eval(&sol, p), p0, "P_0 is not Q(x_R)·y(−R) + …");
        assert!(
            lead.terms.len() == 1 && lead.terms.get(&vec![0; n]) == Some(&1),
            "lead(G) is not the constant 1 at m = 4"
        );
    } else {
        let pnum = interpolate(&xs, &ys, p);
        let p_at = |x: u64| {
            pnum.iter()
                .rev()
                .fold(0u64, |acc, &c| addm(mulm(acc, x, p), c, p))
        };
        assert_eq!(p_at(x_r), y0, "the planted function misses −R");
        sol[..(m - 1)].copy_from_slice(&pnum[1..m]);
    }
    let lead_v = lead.eval(&sol, p);
    if lead_v == 0 {
        return None;
    }
    let mut gsn = vec![0u64; m + 1];
    for i in 0..m {
        let s = mulm(
            gq[i].eval(&sol, p),
            invm(mulm(lead_v, gpow(m - i), p), p),
            p,
        );
        sol[n_c + i] = s;
        gsn[i] = s;
    }
    gsn[m] = 1;
    for &x in &xs {
        let z = mulm(x, invm(g, p), p);
        let val = gsn
            .iter()
            .rev()
            .fold(0u64, |acc, &c| addm(mulm(acc, z, p), c, p));
        assert_eq!(val, 0, "x/g is not a root of G_s");
    }
    let mut rv: Vec<u64> = (0..m).map(|k| u64::from(k == 1)).collect();
    for j in 0..t {
        let b = step_base(j);
        let (q, rn) = npdivmod_monic(&npmul(&rv, &rv, p), &gsn, p);
        sol[b..b + m - 1].copy_from_slice(&q);
        sol[b + m - 1..b + 2 * m - 1].copy_from_slice(&rn);
        rv = rn;
    }
    assert!(
        rv[0] == 1 && rv[1..].iter().all(|&c| c == 0),
        "X^(2^t) ≢ 1 mod G_s at the planted solution"
    );
    // norm_chain_excess (Lean): one equation more than unknowns.
    assert_eq!(eqs.len(), n + 1, "the norm chain's excess is not 1");
    let mut info = json!({"y_neg_r": y0, "lead_at_planted": lead_v, "n_c": n_c, "n_s": n_s,
                      "per_step": per_step});
    if m == 4 {
        merge(
            &mut info,
            &json!({"function": "Q(X)·y − P(X), Q = X + q0, deg P = 2", "q0": sol[0],
                    "lead_G": "1 (checked symbolically)", "excess": eqs.len() - n}),
        );
    }
    Some((eqs, n, sol, info))
}

/// The norm chain's presentation on `V`: for random `m`-subsets of `V`,
/// `X^{2^t} ≡ 1 mod ∏(X − v/g)`; with one element outside `V`, not.
/// Every element of `V` is `g` times a `2^t`-th root of unity.
fn check_norm_chain(d: &Design, m: usize, p: u64, rng: &mut StdRng) {
    let Shape::Kummer { g } = d.shape else {
        panic!("the norm chain is on a Kummer coset")
    };
    let ginv = invm(g, p);
    for &v in &d.v {
        assert_eq!(powm(mulm(v, ginv, p), 1u64 << d.t, p), 1, "{v} ∉ g·μ");
    }
    let set: std::collections::HashSet<u64> = d.v.iter().copied().collect();
    for trial in 0..128 {
        let mut zs: Vec<u64> = Vec::new();
        while zs.len() < m {
            let z = if trial % 2 == 1 && zs.is_empty() {
                let x = rng.gen_range(1..p);
                if set.contains(&x) {
                    continue;
                }
                mulm(x, ginv, p)
            } else {
                mulm(d.v[rng.gen_range(0..d.v.len())], ginv, p)
            };
            if !zs.contains(&z) {
                zs.push(z);
            }
        }
        let mut poly = vec![1u64];
        for &z in &zs {
            poly = npmul(&poly, &[subm(0, z, p), 1], p);
        }
        let r = x_pow2_mod(&poly, d.t, p);
        let one = r[0] == 1 && r[1..].iter().all(|&c| c == 0);
        assert_eq!(
            one,
            trial % 2 == 0,
            "norm-chain membership wrong on trial {trial}"
        );
    }
}

// ── Measuring ──────────────────────────────────────────────────────

/// The planted solution is checked, then F4 (grevlex, the pilot's degree
/// bound `n + d_max + 6`) runs once.
fn measure(eqs: &[Pol], n: usize, p: u64, sol: &[u64], budget: Duration) -> Value {
    for (k, e) in eqs.iter().enumerate() {
        assert_eq!(
            e.eval(sol, p),
            0,
            "equation {k} does not vanish at the planted solution"
        );
    }
    let d_in: Vec<u32> = eqs.iter().map(Pol::degree).collect();
    let d_max = d_in.iter().copied().max().unwrap_or(1);
    let max_degree = n as u32 + d_max + 6;
    let input: Vec<Poly> = eqs
        .iter()
        .map(|e| e.to_poly(p))
        .filter(|f| !f.is_empty())
        .collect();
    let opts = F4Options::new(Ordering::Grevlex, max_degree).with_budget(budget);
    let started = Instant::now();
    let r = f4_fp::f4(&input, n, p, &opts);
    let wall = started.elapsed().as_secs_f64();
    let planted_ok = !r.inconsistent && r.basis.iter().all(|f| f4_fp::eval(f, sol, p) == 0);
    let status = if r.timed_out {
        "timeout"
    } else if !planted_ok {
        "planted_violation"
    } else {
        "ok"
    };
    let reason = match status {
        "timeout" => Value::from(format!("F4 budget of {} s", budget.as_secs())),
        "planted_violation" => Value::from("the planted solution fails F4's basis"),
        _ => Value::Null,
    };
    json!({
        "input_degrees": d_in,
        "max_degree_bound": max_degree,
        "D": r.solving_degree_max,
        "solving_degree_max": r.solving_degree_max,
        "last_productive_degree": r.solving_degree,
        "degree_reached": r.degree_reached,
        "max_cols_to_solution": r.max_cols_to_solution,
        "steps_to_solution": r.steps_to_solution,
        "inconsistent": r.inconsistent,
        "basis_len": r.basis.len(),
        "steps": r.steps,
        "max_rows": r.max_rows,
        "max_cols": r.max_cols,
        "f4_ms": r.ms,
        "time": wall,
        "timed_out": r.timed_out,
        "pairs_above_bound": r.pairs_above_bound,
        "planted_ok": planted_ok,
        "status": status,
        "failure_reason": reason,
    })
}

fn merge(into: &mut Value, from: &Value) {
    if let (Some(a), Some(b)) = (into.as_object_mut(), from.as_object()) {
        for (k, v) in b {
            a.insert(k.clone(), v.clone());
        }
    }
}

/// Amendment 1 item 3, with Amendment 2's prime index (`p₁ = 0`, so every
/// `p₁` seed is unchanged).
fn cell_seed(prime: usize, arm: Arm, m: usize, t: usize) -> u64 {
    SEED ^ ((prime as u64) << 48) ^ (arm.id() << 40) ^ ((m as u64) << 32) ^ ((t as u64) << 16)
}

fn target_seed(cell: u64, k: usize) -> u64 {
    cell ^ ((k as u64 + 1) << 8)
}

/// Build and measure one system; prints a description line (`"phase":
/// "built"`) before F4 runs, and the row after.
fn one(prime: usize, arm: Arm, m: usize, t: usize, k: usize, budget: Duration) -> Value {
    let (mut row, eqs, n, sol) = build_system(prime, arm, m, t, k);
    let mut built = row.clone();
    merge(&mut built, &json!({"phase": "built"}));
    println!("{built}");
    std::io::stdout().flush().ok();
    let r = measure(&eqs, n, prime_of(prime), &sol, budget);
    merge(&mut row, &r);
    row
}

/// One system with every check before measuring: its row so far (no F4
/// fields), equations, unknown count and planted solution.
fn build_system(
    prime: usize,
    arm: Arm,
    m: usize,
    t: usize,
    k: usize,
) -> (Value, Vec<Pol>, usize, Vec<u64>) {
    let p = prime_of(prime);
    let cs = cell_seed(prime, arm, m, t);
    let mut rng = StdRng::seed_from_u64(cs);
    let d = build_design(arm, p, t, &mut rng);
    let (curve, on) = target_curve(p, &d.v, m.max(2), &mut rng);
    let ts = target_seed(cs, k);
    let mut trng = StdRng::seed_from_u64(ts);
    let mut extra = Value::Null;
    let (eqs, n, sol, target) = loop {
        let Some(tg) = plant(&curve, m, &on, &mut trng) else {
            continue;
        };
        if arm == Arm::A11 || arm == Arm::N11 {
            let Shape::Kummer { g } = d.shape else {
                unreachable!()
            };
            let Some((mut eqs, n, sol, mut info)) = norm_chain_system(g, t, &curve, m, &tg) else {
                continue;
            };
            if arm == Arm::N11 {
                // Amendment 2, F2: the m equations defining s_i from the
                // function's coefficients get random coefficients on their
                // own support, shifted to vanish at the planted solution.
                for f in eqs[..m].iter_mut() {
                    let mut g = Pol::zero(n);
                    for e in f.terms.keys() {
                        g.terms.insert(e.clone(), trng.gen_range(1..p));
                    }
                    let c = g.eval(&sol, p);
                    *f = g.sub(&Pol::constant(n, c, p), p);
                }
                merge(
                    &mut info,
                    &json!({"null": "s-defining equations randomised", "null_eqs": m}),
                );
            }
            extra = info;
            break (eqs, n, sol, tg);
        }
        let null = (arm == Arm::N).then_some(&mut trng);
        let (eqs, n, sol) = tower_system(&d, &curve, m, &tg, null);
        break (eqs, n, sol, tg);
    };
    if arm == Arm::A11 || arm == Arm::N11 {
        check_norm_chain(&d, m, p, &mut StdRng::seed_from_u64(ts ^ 0xC4EC));
    }
    let row = json!({
        "arm": arm.name(),
        "m": m,
        "t": t,
        "b": d.bits(),
        "V_size": d.v.len(),
        "target": k,
        "target_kind": "planted",
        "n_vars": n,
        "n_eqs": eqs.len(),
        "p": p,
        "prime": prime_label(prime),
        "seed_base": SEED,
        "cell_seed": cs,
        "target_seed": ts,
        "curve": {"a": curve.a, "b": curve.b},
        "x_r": target.x_r,
        "planted_x": target.pts.iter().map(|pt| pt.x).collect::<Vec<_>>(),
        "set": d.info,
        "presentation": extra,
        "engine": "f4_fp",
        "order": "grevlex",
    });
    (row, eqs, n, sol)
}

// ── The driver: one child process per system ───────────────────────

fn git_commit() -> String {
    let head = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default();
    let dirty = Command::new("git")
        .args([
            "diff",
            "--quiet",
            "HEAD",
            "--",
            "examples",
            "src",
            "Cargo.toml",
            "Cargo.lock",
        ])
        .status()
        .map(|s| !s.success())
        .unwrap_or(true);
    if dirty {
        format!("{head}-dirty")
    } else {
        head
    }
}

fn host() -> Value {
    let name = std::fs::read_to_string("/proc/sys/kernel/hostname")
        .map(|s| s.trim().to_string())
        .unwrap_or_default();
    let cpu = std::fs::read_to_string("/proc/cpuinfo")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("model name"))
                .and_then(|l| l.split(':').nth(1))
                .map(|x| x.trim().to_string())
        })
        .unwrap_or_default();
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(0);
    json!({"hostname": name, "cpu": cpu, "cores": cores})
}

struct RunArgs {
    arms: Vec<Arm>,
    ms: Vec<usize>,
    t_min: usize,
    t_max: usize,
    targets: usize,
    budget: u64,
    vmem_kb: u64,
    grace: u64,
    out_dir: String,
    /// Prime index (`p1` = 0, `p2` = 1).
    prime: usize,
    /// Row file name without extension; `{arm}` and `{m}` are substituted.
    file: String,
}

/// Run one system in a child under `ulimit -v` and a hard wall limit of
/// budget + grace; the row, or a failure row with its reason.
fn run_child(a: &RunArgs, arm: Arm, m: usize, t: usize, k: usize, log: &mut File) -> Value {
    let exe = std::env::current_exe().expect("current exe");
    let mut child = Command::new("sh")
        .arg("-c")
        .arg(format!("ulimit -v {} && exec \"$0\" \"$@\"", a.vmem_kb))
        .arg(&exe)
        .args([
            "one",
            "--prime",
            prime_label(a.prime),
            "--arm",
            arm.name(),
            "--m",
            &m.to_string(),
            "--t",
            &t.to_string(),
            "--target",
            &k.to_string(),
            "--budget",
            &a.budget.to_string(),
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn a child");
    let mut out = child.stdout.take().expect("stdout");
    let mut err = child.stderr.take().expect("stderr");
    let out_t = std::thread::spawn(move || {
        let mut s = String::new();
        out.read_to_string(&mut s).ok();
        s
    });
    let err_t = std::thread::spawn(move || {
        let mut s = String::new();
        err.read_to_string(&mut s).ok();
        s
    });
    let started = Instant::now();
    let limit = Duration::from_secs(a.budget + a.grace);
    let mut killed = false;
    let status = loop {
        if let Some(st) = child.try_wait().expect("wait") {
            break st;
        }
        if started.elapsed() > limit {
            child.kill().ok();
            killed = true;
            break child.wait().expect("wait");
        }
        std::thread::sleep(Duration::from_millis(200));
    };
    let wall = started.elapsed().as_secs_f64();
    let stdout = out_t.join().unwrap_or_default();
    let stderr = err_t.join().unwrap_or_default();
    if !stderr.is_empty() {
        write!(log, "{stderr}").ok();
    }
    let rows: Vec<Value> = stdout
        .lines()
        .filter_map(|l| serde_json::from_str(l).ok())
        .collect();
    let last = rows.last().cloned();
    let mut row = match last {
        Some(r) if r.get("phase").is_none() && status.success() => r,
        other => {
            let mut r = other.unwrap_or_else(|| {
                json!({"arm": arm.name(), "m": m, "t": t, "target": k, "target_kind": "planted",
                       "p": prime_of(a.prime), "prime": prime_label(a.prime)})
            });
            if let Some(o) = r.as_object_mut() {
                o.remove("phase");
            }
            let tail: Vec<&str> = stderr.lines().rev().take(3).collect();
            let (st, why) = if killed {
                (
                    "timeout",
                    format!(
                        "killed at the hard wall limit of {} s (budget {} s + grace {} s): F4 did not return to its deadline check",
                        a.budget + a.grace,
                        a.budget,
                        a.grace
                    ),
                )
            } else if stderr.contains("memory allocation") {
                (
                    "aborted",
                    format!(
                        "memory allocation failed under ulimit -v {} ({status})",
                        a.vmem_kb
                    ),
                )
            } else {
                (
                    "aborted",
                    format!("child exited: {status}; stderr: {}", tail.join(" | ")),
                )
            };
            merge(
                &mut r,
                &json!({"status": st, "failure_reason": why, "planted_ok": null, "D": null,
                        "solving_degree_max": null, "degree_reached": null, "time": wall}),
            );
            r
        }
    };
    merge(
        &mut row,
        &json!({"process_wall_s": wall, "budget_s": a.budget, "ulimit_v_kb": a.vmem_kb}),
    );
    row
}

fn run(a: &RunArgs) {
    std::fs::create_dir_all(&a.out_dir).expect("create --out-dir");
    let commit = git_commit();
    let host = host();
    eprintln!("commit {commit}, host {host}");
    for &arm in &a.arms {
        for &m in &a.ms {
            let name = a
                .file
                .replace("{arm}", arm.name())
                .replace("{m}", &m.to_string());
            let base = format!("{}/{name}", a.out_dir);
            let mut out = OpenOptions::new()
                .create(true)
                .append(true)
                .open(format!("{base}.jsonl"))
                .expect("open rows");
            let mut log = OpenOptions::new()
                .create(true)
                .append(true)
                .open(format!("{base}.log"))
                .expect("open log");
            writeln!(
                log,
                "== {} m={m} {}: commit {commit}, budget {} s, grace {} s, ulimit -v {}, started {}",
                arm.name(),
                prime_label(a.prime),
                a.budget,
                a.grace,
                a.vmem_kb,
                utc_now()
            )
            .ok();
            for t in a.t_min..=a.t_max {
                let mut all_ok = true;
                for k in 0..a.targets {
                    let mut row = run_child(a, arm, m, t, k, &mut log);
                    merge(&mut row, &json!({"git_commit": commit, "host": host}));
                    let line = format!(
                        "{:4} m={m} t={t:2} target={k} b={:.3} |V|={} n_vars={} n_eqs={} D={} reached={} {:.1} s {}{}",
                        arm.name(),
                        row["b"].as_f64().unwrap_or(f64::NAN),
                        row["V_size"],
                        row["n_vars"],
                        row["n_eqs"],
                        row["D"],
                        row["degree_reached"],
                        row["time"].as_f64().unwrap_or(f64::NAN),
                        row["status"].as_str().unwrap_or("?"),
                        row["failure_reason"]
                            .as_str()
                            .map(|s| format!(" ({s})"))
                            .unwrap_or_default(),
                    );
                    eprintln!("{line}");
                    writeln!(log, "{line}").ok();
                    writeln!(out, "{row}").expect("write a row");
                    out.flush().ok();
                    log.flush().ok();
                    if row["status"] != "ok" {
                        all_ok = false;
                    }
                }
                if !all_ok {
                    let msg = format!(
                        "-- {} m={m}: a system at t = {t} hit a limit or failed; the series stops",
                        arm.name()
                    );
                    eprintln!("{msg}");
                    writeln!(log, "{msg}").ok();
                    break;
                }
            }
            writeln!(log, "== {} m={m}: ended {}", arm.name(), utc_now()).ok();
        }
    }
}

fn utc_now() -> String {
    Command::new("date")
        .args(["-u", "+%Y-%m-%dT%H:%M:%SZ"])
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default()
}

// ── The instrument check ───────────────────────────────────────────

fn read_rows(dir: &str) -> Vec<(String, Value)> {
    let mut files: Vec<_> = std::fs::read_dir(dir)
        .expect("read a runs directory")
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "jsonl"))
        .collect();
    files.sort();
    let mut rows = Vec::new();
    for f in files {
        let name = f.file_name().unwrap().to_string_lossy().to_string();
        let text = std::fs::read_to_string(&f).expect("read rows");
        for l in text.lines().filter(|l| !l.trim().is_empty()) {
            if let Ok(v) = serde_json::from_str::<Value>(l) {
                if v.get("summary").is_none() {
                    rows.push((name.clone(), v));
                }
            }
        }
    }
    rows
}

/// The pilot's committed planted Kummer `tower` rows (f4_fp engine, with
/// `solving_degree_max`) at `(m, t, target_index)`.
fn pilot_rows(
    rows: &[(String, Value)],
    m: usize,
    t: usize,
    k: Option<usize>,
) -> Vec<&(String, Value)> {
    rows.iter()
        .filter(|(_, r)| {
            r["kind"] == "kummer"
                && r["control"] == "tower"
                && r["g"] == 1
                && r["p"] == P1
                && r["m"] == m
                && r["t"] == t
                && r["target"] == "planted"
                && r.get("solving_degree_max").is_some_and(|x| !x.is_null())
                && r.get("engine").is_none_or(|e| e == "f4_fp")
                && k.is_none_or(|k| r["target_index"] == k)
        })
        .collect()
}

/// A0 with the pilot's seed rule (`0x504B4D54 ^ kind<<40 ^ m<<32 ^ t<<16 ^
/// g`, kind = Kummer = 0, g = 1) and the pilot's draw order: the same
/// tower, curve and planted targets, measured by this example's code.
fn repro(pilot: &str, m2_max: usize, m3_max: usize, budget: Duration) -> bool {
    let rows = read_rows(pilot);
    let p = P1;
    let mut all = true;
    let mut reached_differs = 0usize;
    let mut systems = 0usize;
    let mut cells = Vec::new();
    for t in 1..=m2_max {
        cells.push((2usize, t));
    }
    for t in 2..=m3_max {
        cells.push((3usize, t));
    }
    let gm = std::env::var("F4_FP_GM").unwrap_or_else(|_| "unset (Gebauer–Möller on)".into());
    println!(
        "A0 instrument check: the pilot's seed rule and draw order, against its committed rows"
    );
    println!("F4_FP_GM = {gm}");
    println!("criterion: x_r, curve, tower, n_vars, n_equations, D = solving_degree_max, inconsistent and planted_ok equal on every committed row;");
    println!("degree_reached (the certification tail, not D) is compared and reported, not scored");
    println!("m  t  target  x_r      ours: D reached inconsistent planted_ok | pilot rows (file: D reached)  verdict");
    for (m, t) in cells {
        let seed = PILOT_SEED ^ ((m as u64) << 32) ^ ((t as u64) << 16) ^ 1;
        let mut rng = StdRng::seed_from_u64(seed);
        let d = build_design(Arm::A0, p, t, &mut rng);
        let (curve, on) = target_curve(p, &d.v, m.max(2), &mut rng);
        let mut targets = Vec::new();
        while targets.len() < 2 {
            if let Some(tg) = plant(&curve, m, &on, &mut rng) {
                targets.push(tg);
            }
        }
        for (k, tg) in targets.iter().enumerate() {
            let (eqs, n, sol) = tower_system(&d, &curve, m, tg, None);
            let r = measure(&eqs, n, p, &sol, budget);
            systems += 1;
            let mine = pilot_rows(&rows, m, t, Some(k));
            let mut ok = !mine.is_empty();
            let mut reached_same = true;
            let mut cites = Vec::new();
            for (f, pr) in &mine {
                let same = pr["x_r"] == tg.x_r
                    && pr["curve"]["a"] == curve.a
                    && pr["curve"]["b"] == curve.b
                    && pr["tower"] == d.info
                    && pr["n_vars"] == n
                    && pr["n_equations"] == eqs.len()
                    && pr["solving_degree_max"] == r["D"]
                    && pr["inconsistent"] == r["inconsistent"]
                    && pr["planted_ok"] == r["planted_ok"];
                let reached = pr["degree_reached"] == r["degree_reached"];
                ok &= same;
                reached_same &= reached;
                cites.push(format!(
                    "{f}: {} {}{}",
                    pr["solving_degree_max"],
                    pr["degree_reached"],
                    if !same {
                        " DIFFERS"
                    } else if !reached {
                        " (reached differs)"
                    } else {
                        ""
                    }
                ));
            }
            all &= ok;
            if !reached_same {
                reached_differs += 1;
            }
            println!(
                "{m}  {t}  {k}       {:<8} {:>2} {:>2} {} {} | {}  {}",
                tg.x_r,
                r["D"],
                r["degree_reached"],
                r["inconsistent"],
                r["planted_ok"],
                if cites.is_empty() {
                    "no pilot row".to_string()
                } else {
                    cites.join("; ")
                },
                if ok { "MATCH" } else { "MISMATCH" }
            );
        }
    }
    println!(
        "instrument check: {} ({systems} systems; degree_reached differs from the pilot on {reached_differs})",
        if all {
            "PASS: every system matches every committed pilot row on the criterion"
        } else {
            "FAIL: MISMATCH"
        }
    );
    all
}

// ── Analysis ───────────────────────────────────────────────────────

/// Least squares `D = α + σ·b`: `(σ, α)`.
fn lsq(pts: &[(f64, f64)]) -> Option<(f64, f64)> {
    if pts.len() < 2 {
        return None;
    }
    let n = pts.len() as f64;
    let mb = pts.iter().map(|x| x.0).sum::<f64>() / n;
    let md = pts.iter().map(|x| x.1).sum::<f64>() / n;
    let sxx: f64 = pts.iter().map(|x| (x.0 - mb).powi(2)).sum();
    if sxx == 0.0 {
        return None;
    }
    let sxy: f64 = pts.iter().map(|x| (x.0 - mb) * (x.1 - md)).sum();
    let s = sxy / sxx;
    Some((s, md - s * mb))
}

fn verdict(cells: usize, sigma: Option<f64>) -> &'static str {
    match sigma {
        _ if cells < 3 => "inconclusive (fewer than 3 completed cells)",
        Some(s) if s <= 0.25 => "candidate (σ ≤ 0.25 on ≥ 3 cells)",
        Some(s) if s < 0.5 => "weak (0.25 < σ < 0.5)",
        Some(_) => "no middle ground at this arity (σ ≥ 0.5)",
        None => "inconclusive (no slope)",
    }
}

/// A cell: `(t, b, D values of its completed rows, statuses)`.
type Cell = (usize, f64, Vec<u64>, Vec<String>);

fn cells_of(rows: &[(String, Value)], arm: &str, m: usize) -> Vec<Cell> {
    let mut by_t: BTreeMap<usize, Cell> = BTreeMap::new();
    for (_, r) in rows {
        if r["arm"] != arm || r["m"] != m {
            continue;
        }
        let t = r["t"].as_u64().unwrap_or(0) as usize;
        let e = by_t.entry(t).or_insert((
            t,
            r["b"].as_f64().unwrap_or(f64::NAN),
            Vec::new(),
            Vec::new(),
        ));
        if e.1.is_nan() {
            e.1 = r["b"].as_f64().unwrap_or(f64::NAN);
        }
        let st = r["status"].as_str().unwrap_or("?").to_string();
        if st == "ok" {
            e.2.push(r["D"].as_u64().expect("an ok row has D"));
        }
        e.3.push(st);
    }
    by_t.into_values().collect()
}

fn complete(c: &Cell) -> bool {
    !c.3.is_empty() && c.3.iter().all(|s| s == "ok")
}

fn mean(xs: &[u64]) -> f64 {
    xs.iter().sum::<u64>() as f64 / xs.len() as f64
}

fn analyze(dir: &str, pilot: Option<&str>) {
    // The first run's rows only: Amendment 2's files (F1-* … F4-*) have their
    // own analysis (`--amendment2`).
    let mut rows = read_rows(dir);
    rows.retain(|(f, _)| group_of(f) == "p1");
    println!("Factor-base designs on p1 = {P1}: σ = ΔD/Δb by least squares over completed cells");
    println!("(D = solving_degree_max of f4_fp grevlex; a cell is completed when every planted system in it is `ok`;");
    println!(" a cell's D is the mean over its systems; rows from {dir})");
    let mut summary = Vec::new();
    for arm in Arm::ALL {
        for m in [2usize, 3] {
            let cells = cells_of(&rows, arm.name(), m);
            if cells.is_empty() {
                continue;
            }
            println!();
            println!("── {} m = {m}", arm.name());
            println!("  t   b       D per system   status");
            for c in &cells {
                println!(
                    "  {:<3} {:<7.3} {:<14} {}",
                    c.0,
                    c.1,
                    format!("{:?}", c.2),
                    c.3.join(",")
                );
            }
            let pts: Vec<(f64, f64)> = cells
                .iter()
                .filter(|c| complete(c))
                .map(|c| (c.1, mean(&c.2)))
                .collect();
            let fit = lsq(&pts);
            let v = verdict(pts.len(), fit.map(|f| f.0));
            match fit {
                Some((s, a)) => println!(
                    "  completed cells {}: σ = {s:.4}, intercept {a:.3}  →  {v}",
                    pts.len()
                ),
                None => println!("  completed cells {}: no slope  →  {v}", pts.len()),
            }
            summary.push((arm.name(), m, pts.len(), fit, v));
        }
    }
    println!();
    println!("Summary");
    println!("  arm  m  cells  σ        verdict");
    for (arm, m, n, fit, v) in &summary {
        println!(
            "  {:<4} {m}  {n:<5}  {:<8} {v}",
            arm,
            fit.map_or("-".to_string(), |f| format!("{:.4}", f.0))
        );
    }

    // The preregistered predictions, as differences at equal b.
    let d_at = |arm: &str, m: usize, t: usize| -> Option<f64> {
        cells_of(&rows, arm, m)
            .into_iter()
            .find(|c| c.0 == t && complete(c))
            .map(|c| mean(&c.2))
    };
    println!();
    println!("Differences at equal t (completed cells only; prediction in brackets)");
    for m in [2usize, 3] {
        for t in 2..=18 {
            let a0 = d_at("A0", m, t);
            let z = d_at("A1z", m, t);
            let g = d_at("A1g", m, t);
            let a4 = d_at("A4", m, t);
            let a0h = d_at("A0", m, t.div_ceil(2));
            let a11 = d_at("A11", m, t);
            if [a0, z, g, a4, a11].iter().all(Option::is_none) {
                continue;
            }
            let f = |x: Option<f64>| x.map_or("-".to_string(), |v| format!("{v:.1}"));
            let df = |x: Option<f64>, y: Option<f64>| match (x, y) {
                (Some(x), Some(y)) => format!("{:+.1}", x - y),
                _ => "-".to_string(),
            };
            println!(
                "  m={m} t={t:<2} A0 {:<5} A1z {:<5} A1g {:<5} A4 {:<5} A11 {:<5} | A1z−A0 {:<5} [±1]  A1g−A1z {:<5} [±1]  A4−A0(⌈t/2⌉) {:<5} [≤ +1 if depth]  A4−A0(t) {:<5} [0 if Bezout]",
                f(a0),
                f(z),
                f(g),
                f(a4),
                f(a11),
                df(z, a0),
                df(g, z),
                df(a4, a0h),
                df(a4, a0),
            );
        }
    }

    if let Some(pdir) = pilot {
        let prow = read_rows(pdir);
        println!();
        println!(
            "A0 (FBDS seeds) against the pilot's committed planted Kummer rows at the same (m, t)"
        );
        println!("(different seeds, so different towers, curves and targets: a cell-by-cell comparison of D)");
        for m in [2usize, 3] {
            for c in cells_of(&rows, "A0", m) {
                let pr = pilot_rows(&prow, m, c.0, None);
                let mut pd: Vec<u64> = pr
                    .iter()
                    .filter(|(_, r)| r["timed_out"] != true)
                    .filter_map(|(_, r)| r["solving_degree_max"].as_u64())
                    .collect();
                pd.sort_unstable();
                pd.dedup();
                let ours: Vec<u64> = {
                    let mut x = c.2.clone();
                    x.sort_unstable();
                    x.dedup();
                    x
                };
                let same = !pd.is_empty() && ours == pd;
                println!(
                    "  m={m} t={:<2} A0 D {:?} ({})  pilot D {:?} ({} rows)  {}",
                    c.0,
                    ours,
                    c.3.join(","),
                    pd,
                    pr.len(),
                    if pr.is_empty() {
                        "no pilot cell"
                    } else if pd.is_empty() {
                        "every pilot row at this cell timed out"
                    } else if same {
                        "same"
                    } else {
                        "DIFFERENT"
                    }
                );
            }
        }
    }
}

// ── Amendment 2: the p₂ instrument check and the p₁ regression ─────

/// F1's instrument check (Amendment 2): `p₂` itself, the arithmetic path
/// `f4_fp` takes there, and one A0 and one A11 planted system at `t = 2`
/// (`m = 2` and `m = 3`, target 0), each checked at its planted solution,
/// measured, and its basis checked to reduce every input to zero.
fn check_p2(budget: Duration) -> bool {
    let p = P2;
    let mut all = true;
    let mut say = |what: &str, ok: bool| {
        println!("  [{}] {what}", if ok { "ok" } else { "FAIL" });
        all &= ok;
    };
    println!("Amendment 2, F1: instrument check at p2 = {p}, before any p2 row");
    println!("the prime");
    say("p2 = 15·2^27 + 1", p == 15 * (1 << 27) + 1);
    say(
        "p2 is prime (trial division by every q ≤ √p2)",
        (2u64..)
            .take_while(|q| q * q <= p)
            .all(|q| !p.is_multiple_of(q)),
    );
    say(
        &format!(
            "2^27 ∥ p2 − 1 (v2 = {}), so t ≤ 27",
            (p - 1).trailing_zeros()
        ),
        (p - 1).trailing_zeros() == 27,
    );
    println!("field arithmetic");
    say(
        "p2 < 2^31 < 2^32, the bound f4_fp's Barrett::new asserts",
        p < 1 << 31,
    );
    say(
        "p2 > 2^16: f4_fp takes the wide path (one Barrett reduction per update), as at p1; the 32-bit and lazy-accumulation paths are for p ≤ 2^16 only",
        p > 1 << 16,
    );
    let worst = u128::from(p - 1) + u128::from(p - 1) * u128::from(p - 1);
    say(
        &format!("largest wide-path operand x + (p − f)·y ≤ (p−1) + (p−1)² = {worst} < 2^64"),
        worst < 1u128 << 64,
    );
    let mut rng = StdRng::seed_from_u64(SEED ^ 0xC4EC);
    let ok = (0..100_000).all(|_| {
        let a = rng.gen_range(1..p);
        let b = rng.gen_range(0..p);
        let exact = (u128::from(a) * u128::from(b) % u128::from(p)) as u64;
        mulm(a, b, p) == exact && mulm(a, invm(a, p), p) == 1 && powm(a, p - 1, p) == 1
    });
    say(
        "harness arithmetic on 100 000 random pairs: a·b against a u128 product, a·a⁻¹ = 1, a^(p−1) = 1",
        ok,
    );
    println!("planted systems at t = 2 (the p2 seeds of F1, target 0)");
    for (arm, m) in [(Arm::A0, 2), (Arm::A11, 2), (Arm::A0, 3), (Arm::A11, 3)] {
        let (row, eqs, n, sol) = build_system(1, arm, m, 2, 0);
        let r = measure(&eqs, n, p, &sol, budget);
        let input: Vec<Poly> = eqs
            .iter()
            .map(|e| e.to_poly(p))
            .filter(|f| !f.is_empty())
            .collect();
        let d_max = eqs.iter().map(Pol::degree).max().unwrap_or(1);
        let opts = F4Options::new(Ordering::Grevlex, n as u32 + d_max + 6).with_budget(budget);
        let f = f4_fp::f4(&input, n, p, &opts);
        let complete = !f.timed_out && f.pairs_above_bound == 0;
        let reduce_zero = input
            .iter()
            .all(|e| f4_fp::reduce(e, &f.basis, p, Ordering::Grevlex).is_empty());
        println!(
            "  {} m={m} t=2: n_vars {n}, n_eqs {}, cell_seed {}, curve {}, x_r {}, D {}, max_cols {}, status {}, planted_ok {}",
            arm.name(),
            eqs.len(),
            row["cell_seed"],
            row["curve"],
            row["x_r"],
            r["D"],
            r["max_cols"],
            r["status"],
            r["planted_ok"]
        );
        say(
            "every equation vanishes at the planted solution (checked before F4)",
            true,
        );
        say(
            "F4 finished, consistent, the planted solution satisfies its basis",
            r["status"] == "ok" && r["planted_ok"] == true,
        );
        say(
            "no pair above the degree bound, and every input reduces to 0 modulo F4's basis",
            complete && reduce_zero,
        );
        say(
            "a second F4 run gives the same D",
            r["D"] == f.solving_degree_max,
        );
    }
    println!(
        "instrument check at p2: {}",
        if all {
            "PASS: F1 may run"
        } else {
            "FAIL: F1 stops"
        }
    );
    all
}

/// The fields a regenerated row must reproduce exactly: D and the
/// structural fields (seeds, set, curve, target, system shape, F4 trace).
const REGRESS_FIELDS: [&str; 27] = [
    "D",
    "solving_degree_max",
    "last_productive_degree",
    "degree_reached",
    "max_cols",
    "max_rows",
    "max_cols_to_solution",
    "steps",
    "steps_to_solution",
    "basis_len",
    "inconsistent",
    "planted_ok",
    "status",
    "pairs_above_bound",
    "input_degrees",
    "max_degree_bound",
    "n_vars",
    "n_eqs",
    "V_size",
    "b",
    "cell_seed",
    "target_seed",
    "curve",
    "x_r",
    "planted_x",
    "set",
    "presentation",
];

/// Rebuild and measure the committed `p₁` rows of one cell with the current
/// code; compare each field of [`REGRESS_FIELDS`] as serialised JSON.
fn regress(file: &str, arm: Arm, m: usize, t: usize, budget: Duration) -> bool {
    let text = std::fs::read_to_string(file).expect("read the committed rows");
    let committed: Vec<Value> = text
        .lines()
        .filter_map(|l| serde_json::from_str::<Value>(l).ok())
        .filter(|r| r["arm"] == arm.name() && r["m"] == m && r["t"] == t)
        .collect();
    println!(
        "Regression: {} m={m} t={t} on p1, rebuilt with the current code (prime index 0), against {file}",
        arm.name()
    );
    println!(
        "fields compared as serialised JSON: {}",
        REGRESS_FIELDS.join(", ")
    );
    let mut all = !committed.is_empty();
    for c in &committed {
        let k = c["target"].as_u64().expect("a target") as usize;
        let (mut row, eqs, n, sol) = build_system(0, arm, m, t, k);
        merge(&mut row, &measure(&eqs, n, P1, &sol, budget));
        let differ: Vec<&str> = REGRESS_FIELDS
            .iter()
            .copied()
            .filter(|f| serde_json::to_vec(&c[*f]).ok() != serde_json::to_vec(&row[*f]).ok())
            .collect();
        all &= differ.is_empty();
        println!(
            "  target {k}: D {} (committed {}), max_cols {} (committed {}), cell_seed {}, target_seed {}: {}",
            row["D"],
            c["D"],
            row["max_cols"],
            c["max_cols"],
            row["cell_seed"],
            row["target_seed"],
            if differ.is_empty() {
                "IDENTICAL on every field".to_string()
            } else {
                format!("DIFFERS on {}", differ.join(", "))
            }
        );
    }
    println!(
        "regression: {}",
        if all {
            "PASS (byte-identical on D and the structural fields)"
        } else {
            "FAIL"
        }
    );
    all
}

// ── Amendment 2: analysis ──────────────────────────────────────────

/// A cell of an Amendment 2 series.
struct Cell2 {
    t: usize,
    b: f64,
    d: Vec<u64>,
    cols: Vec<Option<u64>>,
    time: Vec<f64>,
    status: Vec<String>,
    shape: (u64, u64),
}

impl Cell2 {
    fn complete(&self) -> bool {
        !self.status.is_empty() && self.status.iter().all(|s| s == "ok")
    }

    fn d_mean(&self) -> f64 {
        mean(&self.d)
    }

    /// The mean of its systems' `max_cols` (completed cells only).
    fn cols_mean(&self) -> f64 {
        let c: Vec<u64> = self.cols.iter().flatten().copied().collect();
        mean(&c)
    }
}

/// Which study group a row file belongs to: `F1` … `F4`, or `p1` for the
/// committed p₁ rows of the first run.
fn group_of(file: &str) -> &str {
    match file.get(..3) {
        Some(g @ ("F1-" | "F2-" | "F3-" | "F4-")) => &g[..2],
        _ => "p1",
    }
}

fn cells2(rows: &[(String, Value)], group: &str, arm: &str, m: usize) -> Vec<Cell2> {
    let mut by_t: BTreeMap<usize, Cell2> = BTreeMap::new();
    for (f, r) in rows {
        if group_of(f) != group || r["arm"] != arm || r["m"] != m {
            continue;
        }
        let t = r["t"].as_u64().unwrap_or(0) as usize;
        let e = by_t.entry(t).or_insert(Cell2 {
            t,
            b: f64::NAN,
            d: Vec::new(),
            cols: Vec::new(),
            time: Vec::new(),
            status: Vec::new(),
            shape: (0, 0),
        });
        if e.b.is_nan() {
            e.b = r["b"].as_f64().unwrap_or(f64::NAN);
        }
        e.shape = (
            r["n_vars"].as_u64().unwrap_or(0),
            r["n_eqs"].as_u64().unwrap_or(0),
        );
        let st = r["status"].as_str().unwrap_or("?").to_string();
        if st == "ok" {
            e.d.push(r["D"].as_u64().expect("an ok row has D"));
        }
        e.cols.push(r["max_cols"].as_u64());
        e.time.push(r["time"].as_f64().unwrap_or(f64::NAN));
        e.status.push(st);
    }
    by_t.into_values().collect()
}

/// `(cells, σ fit, Δlog₂ max_cols/Δb fit)` over the completed cells.
type Fits = (usize, Option<(f64, f64)>, Option<(f64, f64)>);

fn fits(cells: &[&Cell2]) -> Fits {
    let done: Vec<&&Cell2> = cells.iter().filter(|c| c.complete()).collect();
    let d: Vec<(f64, f64)> = done.iter().map(|c| (c.b, c.d_mean())).collect();
    let w: Vec<(f64, f64)> = done.iter().map(|c| (c.b, c.cols_mean().log2())).collect();
    (done.len(), lsq(&d), lsq(&w))
}

fn show_series(title: &str, cells: &[Cell2]) -> Fits {
    println!();
    println!("── {title}");
    println!("  t   b      n_vars n_eqs  D per system  max_cols per system   time per system (s)   status");
    for c in cells {
        let cols: Vec<String> = c
            .cols
            .iter()
            .map(|x| x.map_or("-".into(), |v| v.to_string()))
            .collect();
        let time: Vec<String> = c.time.iter().map(|x| format!("{x:.1}")).collect();
        println!(
            "  {:<3} {:<6.3} {:<6} {:<6} {:<13} {:<20} {:<21} {}",
            c.t,
            c.b,
            c.shape.0,
            c.shape.1,
            format!("{:?}", c.d),
            cols.join(","),
            time.join(","),
            c.status.join(",")
        );
    }
    let refs: Vec<&Cell2> = cells.iter().collect();
    let f = fits(&refs);
    print_fits(&f);
    f
}

fn print_fits(f: &Fits) {
    let (n, s, w) = f;
    println!(
        "  completed cells {n}: σ = {}  →  {};  Δlog₂ max_cols/Δb = {} (secondary, descriptive)",
        s.map_or("-".into(), |x| format!("{:.4}", x.0)),
        verdict(*n, s.map(|x| x.0)),
        w.map_or("-".into(), |x| format!("{:.3}", x.0)),
    );
}

fn find(cells: &[Cell2], t: usize) -> Option<&Cell2> {
    cells.iter().find(|c| c.t == t && c.complete())
}

fn analyze2(dir: &str) {
    let rows = read_rows(dir);
    println!("Factor-base designs, Amendment 2 (F1–F4): rows from {dir}");
    println!("(D = solving_degree_max of f4_fp grevlex; a cell is completed when both of its planted systems are `ok`;");
    println!(" a cell's D is the mean over its systems, and so is its max_cols (Amendment 2 does not say how a cell's");
    println!(" width is formed; the mean mirrors the D rule). σ and Δlog₂ max_cols/Δb are least-squares slopes on b");
    println!(" over completed cells. Groups: p1 = the committed rows of the first run; F1–F4 = files F1-*, …, F4-*.)");

    let s = |g: &str, a: &str, m: usize| cells2(&rows, g, a, m);
    let mut secondary: Vec<(String, Fits)> = Vec::new();

    println!();
    println!("════ F1: replication on p2 = {P2} (A0, A11; m = 2, 3)");
    for (a, m) in [("A0", 2), ("A11", 2), ("A0", 3), ("A11", 3)] {
        let f = show_series(&format!("F1 {a} m = {m} p2"), &s("F1", a, m));
        secondary.push((format!("{a} m={m} p2"), f));
    }
    println!();
    println!("════ F2: the A11-shape null N11 on p1 (m = 2, 3); A11 p1 rows (committed) beside it");
    for m in [2usize, 3] {
        let f = show_series(&format!("F2 N11 m = {m} p1"), &s("F2", "N11", m));
        secondary.push((format!("N11 m={m} p1"), f));
        show_series(
            &format!("A11 m = {m} p1 (committed, for the rule)"),
            &s("p1", "A11", m),
        );
    }
    println!();
    println!("════ F3: the m = 3, t = 5 cell on p1, F4 budget 3600 s, hard kill 3900 s");
    for a in ["A0", "A11"] {
        show_series(&format!("F3 {a} m = 3 t = 5 p1"), &s("F3", a, 3));
    }
    println!();
    println!("════ F4: m = 4 on p1 (A0, A11)");
    for a in ["A0", "A11"] {
        let f = show_series(&format!("F4 {a} m = 4 p1"), &s("F4", a, 4));
        secondary.push((format!("{a} m=4 p1"), f));
    }
    let a11_m4: Vec<(u64, u64)> = rows
        .iter()
        .filter(|(f, r)| group_of(f) == "F4" && r["arm"] == "A11")
        .filter_map(|(_, r)| Some((r["n_vars"].as_u64()?, r["n_eqs"].as_u64()?)))
        .collect();
    let excess_ok = !a11_m4.is_empty() && a11_m4.iter().all(|(n, e)| *e == n + 1);
    println!(
        "  A11 m = 4 excess n_eqs − n_vars = 1 on every built row ({} rows): {}",
        a11_m4.len(),
        if excess_ok { "yes" } else { "NO" }
    );

    // ── Decision rules
    println!();
    println!("════ Decision rules (PREREGISTRATION.md, Amendment 2)");

    println!();
    println!("F1. A11 replicates on p2 at m = 3?  Replicates: D(A11) ≤ D(A0) − 1 and max_cols(A11) < max_cols(A0) at every");
    println!("    shared completed t, with ≥ 3 shared completed cells. Fails: any shared completed cell with D(A11) ≥ D(A0).");
    let (a0, a11) = (s("F1", "A0", 3), s("F1", "A11", 3));
    let mut shared = 0usize;
    let (mut every, mut any_fail) = (true, false);
    println!("  t   D A0   D A11   max_cols A0   max_cols A11   D ≤ A0 − 1   cols lower");
    for t in 2..=27 {
        let (Some(x), Some(y)) = (find(&a0, t), find(&a11, t)) else {
            continue;
        };
        shared += 1;
        let dcond = y.d_mean() <= x.d_mean() - 1.0;
        let ccond = y.cols_mean() < x.cols_mean();
        every &= dcond && ccond;
        any_fail |= y.d_mean() >= x.d_mean();
        println!(
            "  {t:<3} {:<6.1} {:<7.1} {:<13.0} {:<14.0} {:<12} {}",
            x.d_mean(),
            y.d_mean(),
            x.cols_mean(),
            y.cols_mean(),
            dcond,
            ccond
        );
    }
    let f1 = if any_fail {
        "FAILS TO REPLICATE"
    } else if shared >= 3 && every {
        "REPLICATES"
    } else {
        "INCONCLUSIVE"
    };
    println!("  shared completed cells: {shared}  →  A11 {f1}");
    println!("  (descriptive, no rule: the same comparison at m = 2 on p2)");
    let (b0, b11) = (s("F1", "A0", 2), s("F1", "A11", 2));
    for t in 2..=27 {
        if let (Some(x), Some(y)) = (find(&b0, t), find(&b11, t)) {
            println!(
                "    m=2 t={t:<2} D A0 {:.1}  D A11 {:.1}  max_cols A0 {:.0}  A11 {:.0}",
                x.d_mean(),
                y.d_mean(),
                x.cols_mean(),
                y.cols_mean()
            );
        }
    }

    println!();
    println!("F2. Shape or structure?  Shape, not structure: D(N11) = D(A11) in every completed cell of F2.");
    println!("    Structure beyond shape: D(N11) ≥ D(A11) + 1 in at least 2 completed cells. Otherwise mixed, cell by cell.");
    println!(
        "    (A11 = the committed p1 rows; a cell enters when it is completed in both series.)"
    );
    let mut f2_cells = 0usize;
    let (mut f2_equal, mut f2_above) = (true, 0usize);
    let mut f2_n11_complete = 0usize;
    for m in [2usize, 3] {
        let (n11, a11) = (s("F2", "N11", m), s("p1", "A11", m));
        f2_n11_complete += n11.iter().filter(|c| c.complete()).count();
        for c in n11.iter().filter(|c| c.complete()) {
            match find(&a11, c.t) {
                Some(a) => {
                    f2_cells += 1;
                    let diff = c.d_mean() - a.d_mean();
                    f2_equal &= diff == 0.0;
                    if diff >= 1.0 {
                        f2_above += 1;
                    }
                    println!(
                        "  m={m} t={:<2} D N11 {:.1}  D A11 {:.1}  N11 − A11 {diff:+.1}   max_cols N11 {:.0}  A11 {:.0}",
                        c.t,
                        c.d_mean(),
                        a.d_mean(),
                        c.cols_mean(),
                        a.cols_mean()
                    );
                }
                None => println!(
                    "  m={m} t={:<2} D N11 {:.1}  A11 has no completed cell here (not compared)",
                    c.t,
                    c.d_mean()
                ),
            }
        }
    }
    let f2 = if f2_cells == 0 {
        "NO SHARED COMPLETED CELL (no verdict)"
    } else if f2_equal {
        "SHAPE, NOT STRUCTURE"
    } else if f2_above >= 2 {
        "STRUCTURE BEYOND SHAPE"
    } else {
        "MIXED (cell by cell above)"
    };
    println!("  completed N11 cells {f2_n11_complete}, compared {f2_cells}, N11 above A11 by ≥ 1 in {f2_above}  →  {f2}");

    println!();
    println!("F3. The m = 3, t = 5 cell on p1, reported as recorded; it joins the m = 3 p1 σ fit only if both systems complete.");
    let mut f3_both = BTreeMap::new();
    for a in ["A0", "A11"] {
        let f3 = s("F3", a, 3);
        let c = f3.iter().find(|c| c.t == 5);
        let both = c.is_some_and(|c| c.complete() && c.status.len() == 2);
        // None until both systems of the cell are recorded.
        f3_both.insert(a, c.is_some_and(|c| c.status.len() == 2).then_some(both));
        let old = s("p1", a, 3);
        let base: Vec<&Cell2> = old.iter().filter(|c| c.t != 5).collect();
        print!("  {a}: ");
        match c {
            Some(c) => println!(
                "D {:?}, max_cols {:?}, time {:?} s, status {:?}",
                c.d, c.cols, c.time, c.status
            ),
            None => println!("no F3 row"),
        }
        print!("    m = 3 p1 fit without F3 (t ≠ 5): ");
        print_fits(&fits(&base));
        if both {
            let mut with = base.clone();
            with.extend(c);
            print!("    m = 3 p1 fit with the F3 cell as a fourth point: ");
            print_fits(&fits(&with));
        } else {
            println!("    the F3 cell is not completed: it does not join the fit");
        }
    }

    println!();
    println!("F4. σ verdicts at m = 4 on p1 (thresholds of §1, derived at m = 4)");
    let mut f4_sigma = BTreeMap::new();
    for a in ["A0", "A11"] {
        let c = s("F4", a, 4);
        let refs: Vec<&Cell2> = c.iter().collect();
        let f = fits(&refs);
        print!("  {a} m = 4: ");
        print_fits(&f);
        f4_sigma.insert(a, f);
    }

    println!();
    println!("════ Secondary (descriptive, not a verdict): Δlog₂ max_cols/Δb per (arm, m, prime), completed cells");
    for a in ["A0", "N", "A1z", "A1g", "A4", "A11"] {
        for m in [2usize, 3] {
            let c = s("p1", a, m);
            if c.is_empty() {
                continue;
            }
            let refs: Vec<&Cell2> = c.iter().collect();
            secondary.push((format!("{a} m={m} p1 (committed)"), fits(&refs)));
        }
    }
    secondary.sort_by(|x, y| x.0.cmp(&y.0));
    println!("  series                       cells  Δlog₂ max_cols/Δb   σ");
    for (name, (n, s, w)) in &secondary {
        println!(
            "  {name:<28} {n:<6} {:<19} {}",
            w.map_or("-".into(), |x| format!("{:.3}", x.0)),
            s.map_or("-".into(), |x| format!("{:.4}", x.0)),
        );
    }

    println!();
    println!("════ Predictions scorecard (Amendment 2)");
    let held = |b: Option<bool>| match b {
        Some(true) => "held",
        Some(false) => "did not hold",
        None => "not decided",
    };
    let f1_held = match f1 {
        "REPLICATES" => Some(true),
        "FAILS TO REPLICATE" => Some(false),
        _ => None,
    };
    println!(
        "  F1: A11 replicates                          observed: {f1:<34} {}",
        held(f1_held)
    );
    let f2_held = match f2 {
        "SHAPE, NOT STRUCTURE" => Some(true),
        x if x.starts_with("NO SHARED") => None,
        _ => Some(false),
    };
    println!(
        "  F2: shape, not structure                    observed: {f2:<34} {}",
        held(f2_held)
    );
    let shown = |b: Option<bool>| b.map_or("not recorded".to_string(), |x| x.to_string());
    println!(
        "  F3: both A11 systems complete within 3600 s observed: A11 both complete = {:<13} {}  (A0 both complete = {}; \"A0's may not\" is not scored)",
        shown(f3_both["A11"]),
        held(f3_both["A11"]),
        shown(f3_both["A0"])
    );
    let (n0, s0, _) = f4_sigma["A0"];
    let a0_held = if n0 < 3 { None } else { s0.map(|x| x.0 >= 0.5) };
    println!(
        "  F4: A0's σ ≥ 0.5 at m = 4                   observed: σ = {} on {n0} cells{:<14} {}",
        s0.map_or("-".into(), |x| format!("{:.4}", x.0)),
        "",
        held(a0_held)
    );
    let (n11, s11, _) = f4_sigma["A11"];
    println!(
        "  F4: A11 at m = 4, no prediction              observed: σ = {} on {n11} cells, {}",
        s11.map_or("-".into(), |x| format!("{:.4}", x.0)),
        verdict(n11, s11.map(|x| x.0))
    );
}

// ── Arguments ──────────────────────────────────────────────────────

/// Compare two row files system by system, keyed by (arm, m, t, target):
/// the construction (seeds, `x_R`, curve, counts) must agree on every
/// shared row; where both systems completed, so must F4's deterministic
/// output (D, `degree_reached`, `max_cols`, `inconsistent`, `planted_ok`).
/// A status that differs (one host finished within the budget, the other
/// did not) is listed, not counted as a disagreement. Returns false on any
/// disagreement.
fn compare_rows(a: &str, b: &str) -> bool {
    let load = |f: &str| -> BTreeMap<(String, u64, u64, u64), Value> {
        std::fs::read_to_string(f)
            .unwrap_or_else(|e| panic!("read {f}: {e}"))
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| serde_json::from_str::<Value>(l).expect("a JSON row"))
            .filter(|r| r.get("summary").is_none())
            .map(|r| {
                let key = (
                    r["arm"].as_str().unwrap_or("?").to_string(),
                    r["m"].as_u64().unwrap_or(0),
                    r["t"].as_u64().unwrap_or(0),
                    r["target"].as_u64().unwrap_or(0),
                );
                (key, r)
            })
            .collect()
    };
    let (ra, rb) = (load(a), load(b));
    const BUILD: [&str; 8] = [
        "p",
        "cell_seed",
        "target_seed",
        "x_r",
        "curve",
        "planted_x",
        "n_vars",
        "n_eqs",
    ];
    const F4: [&str; 5] = [
        "D",
        "degree_reached",
        "max_cols",
        "inconsistent",
        "planted_ok",
    ];
    let (mut same, mut differ, mut status) = (0usize, 0usize, 0usize);
    println!("compare-rows\n  A = {a}\n  B = {b}");
    for (k, x) in &ra {
        let Some(y) = rb.get(k) else { continue };
        let mut bad: Vec<&str> = BUILD.iter().copied().filter(|f| x[*f] != y[*f]).collect();
        let both_ok = x["status"] == "ok" && y["status"] == "ok";
        if both_ok {
            bad.extend(F4.iter().copied().filter(|f| x[*f] != y[*f]));
        }
        let label = format!("{} m={} t={} target={}", k.0, k.1, k.2, k.3);
        if !bad.is_empty() {
            differ += 1;
            println!("  DIFFER {label}: {}", bad.join(", "));
        } else if x["status"] != y["status"] {
            status += 1;
            println!(
                "  status {label}: A {} (D {}, {} s), B {} (D {}, {} s); construction identical",
                x["status"], x["D"], x["time"], y["status"], y["D"], y["time"]
            );
        } else {
            same += 1;
        }
    }
    let only = |p: &BTreeMap<_, Value>, q: &BTreeMap<_, Value>| -> Vec<String> {
        p.keys()
            .filter(|k| !q.contains_key(*k))
            .map(|k: &(String, u64, u64, u64)| {
                format!("{} m={} t={} target={}", k.0, k.1, k.2, k.3)
            })
            .collect()
    };
    let (oa, ob) = (only(&ra, &rb), only(&rb, &ra));
    println!(
        "  shared {}: identical {same}, status differs {status}, disagree {differ}",
        same + status + differ
    );
    println!(
        "  only in A: {}",
        if oa.is_empty() {
            "none".into()
        } else {
            oa.join("; ")
        }
    );
    println!(
        "  only in B: {}",
        if ob.is_empty() {
            "none".into()
        } else {
            ob.join("; ")
        }
    );
    differ == 0
}

fn flag(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .map(|i| args.get(i + 1).cloned().expect("a flag value"))
}

fn list(s: &str) -> Vec<usize> {
    s.split(',').map(|x| x.parse().expect("a number")).collect()
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cmd = args.first().map(String::as_str).unwrap_or("");
    match cmd {
        "one" => {
            let arm = Arm::parse(&flag(&args, "--arm").expect("--arm"));
            let m: usize = flag(&args, "--m").expect("--m").parse().unwrap();
            let t: usize = flag(&args, "--t").expect("--t").parse().unwrap();
            let k: usize = flag(&args, "--target").map_or(0, |s| s.parse().unwrap());
            let budget = flag(&args, "--budget").map_or(300, |s| s.parse().unwrap());
            let prime = flag(&args, "--prime").map_or(0, |s| parse_prime(&s));
            let row = one(prime, arm, m, t, k, Duration::from_secs(budget));
            println!("{row}");
        }
        "run" => {
            let a = RunArgs {
                arms: flag(&args, "--arms").map_or(Arm::ALL.to_vec(), |s| {
                    s.split(',').map(Arm::parse).collect()
                }),
                ms: flag(&args, "--m").map_or(vec![2, 3], |s| list(&s)),
                t_min: flag(&args, "--t-min").map_or(2, |s| s.parse().unwrap()),
                t_max: flag(&args, "--t-max").map_or_else(
                    || {
                        let prime = flag(&args, "--prime").map_or(0, |s| parse_prime(&s));
                        (prime_of(prime) - 1).trailing_zeros() as usize
                    },
                    |s| s.parse().unwrap(),
                ),
                targets: flag(&args, "--targets").map_or(2, |s| s.parse().unwrap()),
                budget: flag(&args, "--budget").map_or(300, |s| s.parse().unwrap()),
                vmem_kb: flag(&args, "--vmem-kb").map_or(13_000_000, |s| s.parse().unwrap()),
                grace: flag(&args, "--grace").map_or(120, |s| s.parse().unwrap()),
                out_dir: flag(&args, "--out-dir").expect("--out-dir"),
                prime: flag(&args, "--prime").map_or(0, |s| parse_prime(&s)),
                file: flag(&args, "--file").unwrap_or_else(|| "{arm}-m{m}".into()),
            };
            assert!(a.t_min >= 2, "t starts at 2");
            let two = (prime_of(a.prime) - 1).trailing_zeros() as usize;
            assert!(a.t_max <= two, "2^t must divide p − 1 (t ≤ {two})");
            run(&a);
        }
        "repro" => {
            let pilot = flag(&args, "--pilot").expect("--pilot");
            let m2 = flag(&args, "--m2-max-t").map_or(6, |s| s.parse().unwrap());
            let m3 = flag(&args, "--m3-max-t").map_or(3, |s| s.parse().unwrap());
            let budget = flag(&args, "--budget").map_or(300, |s| s.parse().unwrap());
            if !repro(&pilot, m2, m3, Duration::from_secs(budget)) {
                std::process::exit(1);
            }
        }
        "analyze" => {
            let dir = args.get(1).expect("a runs directory");
            if args.iter().any(|a| a == "--amendment2") {
                analyze2(dir);
            } else {
                analyze(dir, flag(&args, "--pilot").as_deref());
            }
        }
        "check-p2" => {
            let budget = flag(&args, "--budget").map_or(300, |s| s.parse().unwrap());
            if !check_p2(Duration::from_secs(budget)) {
                std::process::exit(1);
            }
        }
        "compare-rows" => {
            let (a, b) = (args.get(1).expect("rows A"), args.get(2).expect("rows B"));
            if !compare_rows(a, b) {
                std::process::exit(1);
            }
        }
        "regress" => {
            let file = flag(&args, "--rows").expect("--rows");
            let arm = Arm::parse(&flag(&args, "--arm").expect("--arm"));
            let m: usize = flag(&args, "--m").expect("--m").parse().unwrap();
            let t: usize = flag(&args, "--t").expect("--t").parse().unwrap();
            let budget = flag(&args, "--budget").map_or(300, |s| s.parse().unwrap());
            if !regress(&file, arm, m, t, Duration::from_secs(budget)) {
                std::process::exit(1);
            }
        }
        _ => {
            eprintln!(
                "usage: fb_designs_slope run|one|repro|analyze|check-p2|regress|compare-rows …  (see the module docs)"
            );
            std::process::exit(2);
        }
    }
}
