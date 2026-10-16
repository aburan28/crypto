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
/// `p₂ = 15·2^27 + 1`, the replication prime (§5, Amendment 2 F1).
const P2: u64 = 2_013_265_921;

/// The prime index of the seed rule (Amendment 2): `p₁ = 0`, `p₂ = 1`.
fn prime_index(p: u64) -> u64 {
    match p {
        P1 => 0,
        P2 => 1,
        _ => panic!("no prime index registered for p = {p}"),
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
}

impl Arm {
    const ALL: [Arm; 6] = [Arm::A0, Arm::N, Arm::A1z, Arm::A1g, Arm::A4, Arm::A11];

    fn name(self) -> &'static str {
        match self {
            Arm::A0 => "A0",
            Arm::N => "N",
            Arm::A1z => "A1z",
            Arm::A1g => "A1g",
            Arm::A4 => "A4",
            Arm::A11 => "A11",
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
        Arm::A0 | Arm::N | Arm::A11 => kummer_design(p, t, rng),
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

/// The norm-chain system (`A11`) for `m ∈ {2, 3}` on the coset `g·μ_{2^t}`.
///
/// Unknowns: `c_1 … c_{m−1}` (the function `f = y − P(X)` through `−R`,
/// `P(X) = y(−R) + Σ c_i (X^i − x_R^i)`; `c_1 = λ` at `m = 2`), then
/// `s_0 … s_{m−1}` (`G_s = X^m + Σ s_i X^i`), then per step `j < t` the
/// quotient `q_j` (`m − 1` coefficients) and `r_{j+1}` (`m` coefficients).
/// Returns `None` when the planted function has a vanishing leading
/// coefficient (the caller replants).
fn norm_chain_system(
    g: u64,
    t: usize,
    curve: &PrimeCurve,
    m: usize,
    target: &Target,
) -> Option<(Vec<Pol>, usize, Vec<u64>, Value)> {
    assert!(
        m == 2 || m == 3,
        "the norm chain is written for m ∈ {{2, 3}}"
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
    let mut pc: Vec<Pol> = vec![Pol::constant(n, y0, p)];
    for i in 1..m {
        let xri = powm(x_r, i as u64, p);
        pc[0] = pc[0].sub(&cvar(i).scale(xri, p), p);
        pc.push(cvar(i));
    }
    let deg_n = 3.max(2 * (m - 1));
    let mut nc: Vec<Pol> = vec![Pol::zero(n); deg_n + 1];
    nc[3] = nc[3].add(&Pol::constant(n, 1, p), p);
    nc[1] = nc[1].add(&Pol::constant(n, ca, p), p);
    nc[0] = nc[0].add(&Pol::constant(n, cb, p), p);
    for i in 0..m {
        for j in 0..m {
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
    let p_at = |x: u64| {
        pnum.iter()
            .rev()
            .fold(0u64, |acc, &c| addm(mulm(acc, x, p), c, p))
    };
    assert_eq!(p_at(x_r), y0, "the planted function misses −R");
    let mut sol = vec![0u64; n];
    sol[..(m - 1)].copy_from_slice(&pnum[1..m]);
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
    let info = json!({"y_neg_r": y0, "lead_at_planted": lead_v, "n_c": n_c, "n_s": n_s,
                      "per_step": per_step});
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

/// `0x46424453 ⊕ prime≪48 ⊕ arm≪40 ⊕ m≪32 ⊕ t≪16` (Amendments 1 and 2;
/// every p₁ seed is unchanged).
fn cell_seed(p: u64, arm: Arm, m: usize, t: usize) -> u64 {
    SEED ^ (prime_index(p) << 48) ^ (arm.id() << 40) ^ ((m as u64) << 32) ^ ((t as u64) << 16)
}

fn target_seed(cell: u64, k: usize) -> u64 {
    cell ^ ((k as u64 + 1) << 8)
}

/// Build and measure one system; prints a description line (`"phase":
/// "built"`) before F4 runs, and the row after.
fn one(p: u64, arm: Arm, m: usize, t: usize, k: usize, budget: Duration) -> Value {
    let cs = cell_seed(p, arm, m, t);
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
        if arm == Arm::A11 {
            let Shape::Kummer { g } = d.shape else {
                unreachable!()
            };
            let Some((eqs, n, sol, info)) = norm_chain_system(g, t, &curve, m, &tg) else {
                continue;
            };
            extra = info;
            break (eqs, n, sol, tg);
        }
        let null = (arm == Arm::N).then_some(&mut trng);
        let (eqs, n, sol) = tower_system(&d, &curve, m, &tg, null);
        break (eqs, n, sol, tg);
    };
    if arm == Arm::A11 {
        check_norm_chain(&d, m, p, &mut StdRng::seed_from_u64(ts ^ 0xC4EC));
    }
    let mut row = json!({
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
        "prime_index": prime_index(p),
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
    let mut built = row.clone();
    merge(&mut built, &json!({"phase": "built"}));
    println!("{built}");
    std::io::stdout().flush().ok();
    let r = measure(&eqs, n, p, &sol, budget);
    merge(&mut row, &r);
    row
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
    p: u64,
    arms: Vec<Arm>,
    ms: Vec<usize>,
    t_min: usize,
    t_max: usize,
    targets: usize,
    budget: u64,
    vmem_kb: u64,
    grace: u64,
    out_dir: String,
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
            "--p",
            &a.p.to_string(),
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
                json!({"arm": arm.name(), "m": m, "t": t, "target": k, "target_kind": "planted"})
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
            let base = format!("{}/{}-m{m}", a.out_dir, arm.name());
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
                "== {} m={m}: commit {commit}, budget {} s, ulimit -v {}, started {}",
                arm.name(),
                a.budget,
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

/// The mean `max_cols` of a completed cell's systems, per `t`.
fn widths_of(rows: &[(String, Value)], arm: &str, m: usize) -> BTreeMap<usize, f64> {
    let mut by_t: BTreeMap<usize, Vec<u64>> = BTreeMap::new();
    for (_, r) in rows {
        if r["arm"] == arm && r["m"] == m && r["status"] == "ok" {
            let t = r["t"].as_u64().unwrap_or(0) as usize;
            by_t.entry(t)
                .or_default()
                .push(r["max_cols"].as_u64().expect("an ok row has max_cols"));
        }
    }
    by_t.into_iter().map(|(t, w)| (t, mean(&w))).collect()
}

fn analyze(dir: &str, pilot: Option<&str>) {
    let rows = read_rows(dir);
    let mut primes: Vec<u64> = rows.iter().filter_map(|(_, r)| r["p"].as_u64()).collect();
    primes.sort_unstable();
    primes.dedup();
    println!(
        "Factor-base designs on p = {primes:?}: σ = ΔD/Δb by least squares over completed cells"
    );
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
            // Secondary, descriptive (Amendment 2): Δ log₂ max_cols / Δb.
            let widths = widths_of(&rows, arm.name(), m);
            let wpts: Vec<(f64, f64)> = cells
                .iter()
                .filter(|c| complete(c))
                .filter_map(|c| widths.get(&c.0).map(|w| (c.1, w.log2())))
                .collect();
            if let Some((ws, _)) = lsq(&wpts) {
                println!(
                    "  width growth Δlog2(max_cols)/Δb = {ws:.4} over {} cells (descriptive)",
                    wpts.len()
                );
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

    // Amendment 2 F1: A11 against A0 at shared completed t (cell D and cell
    // mean max_cols); the rule is registered for m = 3 on p₂.
    println!();
    println!(
        "A11 against A0 at shared completed t (Amendment 2 F1 rule; registered for m = 3 on p2)"
    );
    for m in [2usize, 3] {
        let a11: BTreeMap<usize, f64> = cells_of(&rows, "A11", m)
            .into_iter()
            .filter(complete)
            .map(|c| (c.0, mean(&c.2)))
            .collect();
        let a0: BTreeMap<usize, f64> = cells_of(&rows, "A0", m)
            .into_iter()
            .filter(complete)
            .map(|c| (c.0, mean(&c.2)))
            .collect();
        let (w11, w0) = (widths_of(&rows, "A11", m), widths_of(&rows, "A0", m));
        let shared: Vec<usize> = a11.keys().filter(|t| a0.contains_key(t)).copied().collect();
        if shared.is_empty() {
            continue;
        }
        let mut all_pass = true;
        let mut any_fail = false;
        for t in &shared {
            let (d11, d0) = (a11[t], a0[t]);
            let (c11, c0) = (w11[t], w0[t]);
            let pass = d11 <= d0 - 1.0 && c11 < c0;
            all_pass &= pass;
            any_fail |= d11 >= d0;
            println!(
                "  m={m} t={t:<2} A11 D {d11:.1} cols {c11:.0} | A0 D {d0:.1} cols {c0:.0} | {}",
                if pass {
                    "A11 D ≤ A0 D − 1 and narrower"
                } else {
                    "condition not met"
                }
            );
        }
        let v = if any_fail {
            "fails to replicate (a shared cell has A11 D ≥ A0 D)"
        } else if all_pass && shared.len() >= 3 {
            "replicates"
        } else {
            "inconclusive"
        };
        if m == 3 {
            println!("  m={m}: {} shared completed cells → {v}", shared.len());
        } else {
            println!(
                "  m={m}: {} shared completed cells (descriptive: the rule is registered for m = 3)",
                shared.len()
            );
        }
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

// ── Arguments ──────────────────────────────────────────────────────

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
            let p = flag(&args, "--p").map_or(P1, |s| s.parse().unwrap());
            let row = one(p, arm, m, t, k, Duration::from_secs(budget));
            println!("{row}");
        }
        "run" => {
            let p = flag(&args, "--p").map_or(P1, |s| s.parse().unwrap());
            let two_adic = (p - 1).trailing_zeros() as usize;
            let a = RunArgs {
                p,
                arms: flag(&args, "--arms").map_or(Arm::ALL.to_vec(), |s| {
                    s.split(',').map(Arm::parse).collect()
                }),
                ms: flag(&args, "--m").map_or(vec![2, 3], |s| list(&s)),
                t_min: flag(&args, "--t-min").map_or(2, |s| s.parse().unwrap()),
                t_max: flag(&args, "--t-max").map_or(two_adic, |s| s.parse().unwrap()),
                targets: flag(&args, "--targets").map_or(2, |s| s.parse().unwrap()),
                budget: flag(&args, "--budget").map_or(300, |s| s.parse().unwrap()),
                vmem_kb: flag(&args, "--vmem-kb").map_or(13_000_000, |s| s.parse().unwrap()),
                grace: flag(&args, "--grace").map_or(120, |s| s.parse().unwrap()),
                out_dir: flag(&args, "--out-dir").expect("--out-dir"),
            };
            assert!(a.t_min >= 2, "t starts at 2");
            assert!(a.t_max <= two_adic, "2^t must divide p − 1");
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
            analyze(dir, flag(&args, "--pilot").as_deref());
        }
        _ => {
            eprintln!("usage: fb_designs_slope run|one|repro|analyze …  (see the module docs)");
            std::process::exit(2);
        }
    }
}
