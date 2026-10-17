//! `volcano-rho` — one fixed ECDLP instance carried down every stratum of
//! an `ℓ`-isogeny volcano, and Pollard rho timed on each curve.
//!
//! 1. Search a random prime-field curve `E/F_p` (`p ≡ 3 mod 4` for ℓ = 3, `p ≡ 1 mod 8` for ℓ = 2) whose
//!    Frobenius discriminant `D = t² − 4p = f²·D_K` has `v_ℓ(f) ≥ depth`,
//!    with `#E = h·n`, `n` prime and `h` small.  The `ℓ`-volcano of `E`
//!    then has `v_ℓ(f)` levels below its crater.
//! 2. Walk the `ℓ`-isogeny graph breadth first from `E`.  Each edge is an
//!    explicit Vélu isogeny from a rational root of the `ℓ`-division
//!    polynomial (`ℓ = 2`: the cubic; `ℓ = 3`: `ψ₃`), keyed by `j`.
//!    A vertex's level is `depth − dist(v, floor)`, the floor being the
//!    vertices with exactly one rational `ℓ`-kernel.
//! 3. Fix `P` of order `n` and `Q = kP` on `E`, and push `(P, Q)` through
//!    the edge maps to every curve, so each curve holds the same discrete
//!    log `k` (`ℓ ∤ n`, so every map is injective on `⟨P⟩`).  Every edge is
//!    checked: image on the codomain, `φ(P + Q) = φ(P) + φ(Q)`, `n·φ(P) = O`.
//! 4. Run single-walker distinguished-point Pollard rho (32-adding walk,
//!    no negation map) on every curve with the same trial seeds, plus
//!    isomorphic re-scalings of the root as a representation control, and
//!    record group operations and wall time per solve.
//!
//! Isogenous curves share `#E(F_p)`, so the expected rho cost is the same
//! on every curve; this tool measures how the observed costs spread.

use std::collections::{HashMap, HashSet, VecDeque};
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use crypto_lib::cryptanalysis::isogeny_walk::field::Field;
use crypto_lib::cryptanalysis::isogeny_walk::poly;
use num_bigint::BigUint;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde_json::json;

// ---------------------------------------------------------------- F_p, p < 2^62

#[derive(Clone, Copy)]
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
    fn sm(&self, k: u64, a: u64) -> u64 {
        self.mul(k % self.p, a)
    }
    fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1;
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
        inv_mod(a, self.p)
    }
    /// Tonelli–Shanks (a single power when `p ≡ 3 mod 4`); `None` for a
    /// non-residue.
    fn sqrt(&self, a: u64) -> Option<u64> {
        if a == 0 {
            return Some(0);
        }
        if self.pow(a, (self.p - 1) / 2) != 1 {
            return None;
        }
        if self.p % 4 == 3 {
            return Some(self.pow(a, (self.p + 1) / 4));
        }
        let s = (self.p - 1).trailing_zeros();
        let q = (self.p - 1) >> s;
        let mut z = 2;
        while self.pow(z, (self.p - 1) / 2) != self.p - 1 {
            z += 1;
        }
        let (mut m, mut c) = (s, self.pow(z, q));
        let (mut t, mut r) = (self.pow(a, q), self.pow(a, q.div_ceil(2)));
        while t != 1 {
            let (mut i, mut t2) = (0, t);
            while t2 != 1 {
                t2 = self.mul(t2, t2);
                i += 1;
            }
            let b = self.pow(c, 1 << (m - i - 1));
            m = i;
            c = self.mul(b, b);
            t = self.mul(t, c);
            r = self.mul(r, b);
        }
        Some(r)
    }
}

fn inv_mod(a: u64, m: u64) -> u64 {
    let (mut r0, mut r1) = (m as i128, a as i128);
    let (mut s0, mut s1) = (0i128, 1i128);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (s0, s1) = (s1, s0 - q * s1);
    }
    assert_eq!(r0, 1, "not invertible");
    s0.rem_euclid(m as i128) as u64
}

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

fn powmod(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1 % m;
    a %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, m);
        }
        a = mulmod(a, a, m);
        e >>= 1;
    }
    r
}

/// Deterministic Miller–Rabin for `n < 2^64`.
fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for q in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n.is_multiple_of(q) {
            return n == q;
        }
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
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

fn factor(mut n: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut d = 2u64;
    while d * d <= n {
        if n.is_multiple_of(d) {
            let mut e = 0;
            while n.is_multiple_of(d) {
                n /= d;
                e += 1;
            }
            out.push((d, e));
        }
        d += if d == 2 { 1 } else { 2 };
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

/// The largest prime factor of `n < 2^64`: trial division to 2^16, then
/// Pollard–Brent on what is left.
fn largest_prime_factor(mut n: u64) -> u64 {
    let mut best = 1;
    let mut d = 2u64;
    while d < 1 << 16 && d * d <= n {
        while n.is_multiple_of(d) {
            n /= d;
            best = d;
        }
        d += if d == 2 { 1 } else { 2 };
    }
    let mut stack = vec![n];
    while let Some(m) = stack.pop() {
        if m == 1 {
            continue;
        }
        if is_prime(m) {
            best = best.max(m);
            continue;
        }
        let g = brent(m);
        stack.push(g);
        stack.push(m / g);
    }
    best
}

/// A nontrivial factor of the composite `n` (Pollard–Brent).
fn brent(n: u64) -> u64 {
    for c in 1u64.. {
        let f = |x: u64| (mulmod(x, x, n) + c) % n;
        let (mut x, mut y, mut g) = (2u64, 2u64, 1u64);
        while g == 1 {
            x = f(x);
            y = f(f(y));
            g = gcd(x.abs_diff(y), n);
        }
        if g != n {
            return g;
        }
    }
    unreachable!()
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

// ---------------------------------------------------------------- curves

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Curve {
    a: u64,
    b: u64,
}

type Pt = Option<(u64, u64)>;

impl Curve {
    fn rhs(&self, f: &Fp, x: u64) -> u64 {
        f.add(f.mul(f.add(f.mul(x, x), self.a), x), self.b)
    }
    fn on(&self, f: &Fp, p: Pt) -> bool {
        p.is_none_or(|(x, y)| f.mul(y, y) == self.rhs(f, x))
    }
    fn j(&self, f: &Fp) -> Option<u64> {
        let a3 = f.sm(4, f.mul(f.mul(self.a, self.a), self.a));
        let den = f.add(a3, f.sm(27, f.mul(self.b, self.b)));
        (den != 0).then(|| f.mul(f.sm(1728, a3), f.inv(den)))
    }
    fn add(&self, f: &Fp, p: Pt, q: Pt) -> Pt {
        let (Some((x1, y1)), Some((x2, y2))) = (p, q) else {
            return p.or(q);
        };
        let l = if x1 == x2 {
            if f.add(y1, y2) == 0 {
                return None;
            }
            f.mul(f.add(f.sm(3, f.mul(x1, x1)), self.a), f.inv(f.add(y1, y1)))
        } else {
            f.mul(f.sub(y2, y1), f.inv(f.sub(x2, x1)))
        };
        let x3 = f.sub(f.sub(f.mul(l, l), x1), x2);
        Some((x3, f.sub(f.mul(l, f.sub(x1, x3)), y1)))
    }
    fn mul(&self, f: &Fp, p: Pt, k: u64) -> Pt {
        let (mut r, mut b, mut k) = (None, p, k);
        while k > 0 {
            if k & 1 == 1 {
                r = self.add(f, r, b);
            }
            b = self.add(f, b, b);
            k >>= 1;
        }
        r
    }
    fn neg(f: &Fp, p: Pt) -> Pt {
        p.map(|(x, y)| (x, f.neg(y)))
    }
    fn random_point(&self, f: &Fp, rng: &mut SmallRng) -> (u64, u64) {
        loop {
            let x = rng.gen_range(0..f.p);
            if let Some(y) = f.sqrt(self.rhs(f, x)) {
                if y != 0 {
                    return (x, y);
                }
            }
        }
    }
}

/// The order of `E(F_p)`: every `M` in the Hasse interval with `M·P = O`
/// by baby-step giant-step, intersected over random points until one is
/// left.  `None` if ten points do not pin it down.
fn group_order(f: &Fp, e: &Curve, rng: &mut SmallRng) -> Option<u64> {
    let s = (4.0 * (f.p as f64).sqrt()).ceil() as u64 + 2;
    let lo = f.p + 1 - s / 2;
    let m = (s as f64).sqrt().ceil() as u64 + 1;
    let mut cands: Option<HashSet<u64>> = None;
    for _ in 0..10 {
        let p = Some(e.random_point(f, rng));
        let mut baby: HashMap<Pt, u64> = HashMap::with_capacity(m as usize);
        let mut jp: Pt = None;
        for j in 0..m {
            baby.entry(Curve::neg(f, jp)).or_insert(j);
            jp = e.add(f, jp, p);
        }
        let step = e.mul(f, p, m);
        let mut g = e.mul(f, p, lo);
        let mut found = HashSet::new();
        for i in 0..=m {
            if let Some(&j) = baby.get(&g) {
                // Every j with −jP = g differs by a multiple of ord(P) < m
                // only when ord(P) is tiny; the intersection handles that.
                found.insert(lo + i * m + j);
            }
            g = e.add(f, g, step);
        }
        cands = Some(match cands {
            None => found,
            Some(c) => c.intersection(&found).copied().collect(),
        });
        let c = cands.as_ref().unwrap();
        let hasse: Vec<u64> = c
            .iter()
            .copied()
            .filter(|&v| (v as i128 - (f.p as i128 + 1)).pow(2) <= 4 * f.p as i128)
            .collect();
        if hasse.len() == 1 {
            return Some(hasse[0]);
        }
    }
    None
}

// ---------------------------------------------------------------- isogenies

/// The rational `x`-coordinates generating `ℓ`-kernels (`ℓ ∈ {2, 3}`).
fn kernel_xs(f: &Fp, field: &Field, e: &Curve, ell: u64) -> Vec<u64> {
    let fe = |v: u64| field.from_u64(v);
    let h = if ell == 2 {
        // x³ + ax + b
        vec![fe(e.b), fe(e.a), fe(0), fe(1)]
    } else {
        // ψ₃ = 3x⁴ + 6ax² + 12bx − a²
        vec![
            fe(f.neg(f.mul(e.a, e.a))),
            fe(f.sm(12, e.b)),
            fe(f.sm(6, e.a)),
            fe(0),
            fe(3),
        ]
    };
    let h = poly::trim(field, h);
    poly::roots(field, &h, 0x5eed)
        .iter()
        .map(|r| field.to_big(r).iter_u64_digits().next().unwrap_or(0))
        .collect()
}

/// Vélu's isogeny with kernel `⟨(x0, ·)⟩` of order `ℓ ∈ {2, 3}`.
#[derive(Clone, Copy, Debug)]
struct Velu {
    ell: u64,
    x0: u64,
    tq: u64,
    uq: u64,
    codomain: Curve,
}

fn velu(f: &Fp, e: &Curve, ell: u64, x0: u64) -> Velu {
    let gx = f.add(f.sm(3, f.mul(x0, x0)), e.a);
    let (tq, uq) = if ell == 2 {
        (gx, 0)
    } else {
        (f.add(gx, gx), f.sm(4, e.rhs(f, x0)))
    };
    let w = f.add(uq, f.mul(x0, tq));
    Velu {
        ell,
        x0,
        tq,
        uq,
        codomain: Curve {
            a: f.sub(e.a, f.sm(5, tq)),
            b: f.sub(e.b, f.sm(7, w)),
        },
    }
}

impl Velu {
    fn map(&self, f: &Fp, p: Pt) -> Pt {
        let (x, y) = p?;
        if x == self.x0 {
            return None;
        }
        let i = f.inv(f.sub(x, self.x0));
        let i2 = f.mul(i, i);
        let i3 = f.mul(i2, i);
        let xx = f.add(x, f.add(f.mul(self.tq, i), f.mul(self.uq, i2)));
        let yy = f.sub(
            y,
            f.mul(
                y,
                f.add(f.mul(f.add(self.uq, self.uq), i3), f.mul(self.tq, i2)),
            ),
        );
        Some((xx, yy))
    }
}

// ---------------------------------------------------------------- pull-back walk

/// One step back toward the root curve: a dual isogeny followed by the
/// isomorphism `(x, y) ↦ (u²x, u³y)` onto the previous curve's model, or
/// that isomorphism alone (an isomorphic control).
#[derive(Clone, Copy)]
enum Pull {
    Dual(Velu, u64),
    Scale(u64),
}

/// How a walk reads a point: the point itself (`native`), or its canonical
/// image `[ℓ^lift] · φ̂_1 ∘ … ∘ φ̂_m(R)` on the root curve (`pullback`), which
/// is `±ℓ^M φ⁻¹(R)` for every curve, so every curve steers by the same
/// numbers.
#[derive(Clone)]
struct Lens {
    chain: Vec<Pull>,
    root: Curve,
    lift: u64,
    native: bool,
}

impl Lens {
    fn native() -> Self {
        Lens {
            chain: Vec::new(),
            root: Curve { a: 0, b: 0 },
            lift: 1,
            native: true,
        }
    }
    fn read(&self, f: &Fp, pt: (u64, u64)) -> (u64, u64) {
        if self.native {
            return pt;
        }
        let mut q = pt;
        for step in &self.chain {
            let (v, u) = match step {
                Pull::Dual(v, u) => (Some(v), *u),
                Pull::Scale(u) => (None, *u),
            };
            if let Some(v) = v {
                q = v
                    .map(f, Some(q))
                    .expect("a point of order n is never in a dual kernel");
            }
            let u2 = f.mul(u, u);
            q = (f.mul(u2, q.0), f.mul(f.mul(u2, u), q.1));
        }
        self.root
            .mul(f, Some(q), self.lift)
            .expect("l^lift is invertible mod n")
    }
}

// ---------------------------------------------------------------- rho

struct Instance {
    curve: Curve,
    p: (u64, u64),
    q: (u64, u64),
}

struct Solve {
    k: Option<u64>,
    iters: u64,
    restarts: u32,
    ns: u128,
}

fn rho(f: &Fp, inst: &Instance, n: u64, seed: u64, dp_bits: u32, lens: &Lens) -> Solve {
    let e = &inst.curve;
    let (p, q) = (Some(inst.p), Some(inst.q));
    let mut rng = SmallRng::seed_from_u64(seed);
    let t0 = Instant::now();
    let steps: Vec<(u64, u64, (u64, u64))> = (0..32)
        .map(|_| loop {
            let (a, b) = (rng.gen_range(0..n), rng.gen_range(0..n));
            if let Some(r) = e.add(f, e.mul(f, p, a), e.mul(f, q, b)) {
                break (a, b, r);
            }
        })
        .collect();
    let mut table: HashMap<u64, (u64, u64, u64)> = HashMap::new();
    let (mut iters, mut restarts) = (0u64, 0u32);
    let dp_mask = (1u64 << dp_bits) - 1;
    let max_run = 40u64 << dp_bits;
    'restart: loop {
        let (mut a, mut b) = (rng.gen_range(0..n), rng.gen_range(0..n));
        let Some(mut x) = e.add(f, e.mul(f, p, a), e.mul(f, q, b)) else {
            restarts += 1;
            continue;
        };
        let mut run = 0u64;
        let mut c = lens.read(f, x);
        loop {
            let s = &steps[(c.0.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 59) as usize];
            iters += 1;
            run += 1;
            a = (a + s.0) % n;
            b = (b + s.1) % n;
            match e.add(f, Some(x), Some(s.2)) {
                Some(nx) => {
                    x = nx;
                    c = lens.read(f, x);
                }
                None => {
                    // aP + bQ = O.
                    if b != 0 {
                        let k = mulmod(n - a % n, inv_mod(b, n), n) % n;
                        return Solve {
                            k: Some(k),
                            iters,
                            restarts,
                            ns: t0.elapsed().as_nanos(),
                        };
                    }
                    restarts += 1;
                    continue 'restart;
                }
            }
            if c.0.wrapping_mul(0xC2B2_AE3D_27D4_EB4F) >> 40 & dp_mask == 0 {
                if let Some(&(y2, a2, b2)) = table.get(&c.0) {
                    // Same point: (a − a2)P = (b2 − b)Q; opposite: (a + a2)P = −(b + b2)Q.
                    let (num, den) = if y2 == c.1 {
                        ((a + n - a2) % n, (b2 + n - b) % n)
                    } else {
                        ((n - (a + a2) % n) % n, (b + b2) % n)
                    };
                    if den != 0 {
                        let k = mulmod(num, inv_mod(den, n), n);
                        return Solve {
                            k: Some(k),
                            iters,
                            restarts,
                            ns: t0.elapsed().as_nanos(),
                        };
                    }
                    restarts += 1;
                    continue 'restart;
                }
                table.insert(c.0, (c.1, a, b));
                run = 0;
            }
            if run > max_run {
                restarts += 1;
                continue 'restart;
            }
        }
    }
}

/// One solve by the cryptanalysis repository's `ca solve --alg rho`:
/// group operations and time as `ca` reports them.
fn rho_ca(ca: &str, extra: &[String], f: &Fp, inst: &Instance, n: u64, seed: u64) -> Solve {
    let out = std::process::Command::new(ca)
        .args(["solve", "--alg", "rho", "--group", "ec"])
        .args(["--p", &f.p.to_string()])
        .args([
            "--a",
            &inst.curve.a.to_string(),
            "--b",
            &inst.curve.b.to_string(),
        ])
        .args(["--order", &n.to_string()])
        .args(["--g", &format!("{},{}", inst.p.0, inst.p.1)])
        .args(["--h", &format!("{},{}", inst.q.0, inst.q.1)])
        .args(["--seed", &seed.to_string()])
        .args(extra)
        .output()
        .unwrap_or_else(|e| die(&format!("cannot run {ca}: {e}")));
    let v: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap_or_else(|_| {
        die(&format!(
            "{ca} printed no JSON: {}",
            String::from_utf8_lossy(&out.stderr)
        ))
    });
    Solve {
        k: v["x"].as_u64(),
        iters: v["ops"].as_u64().unwrap_or(0),
        restarts: 0,
        ns: (v["seconds"].as_f64().unwrap_or(0.0) * 1e9) as u128,
    }
}

// ---------------------------------------------------------------- statistics

fn stats(v: &[f64]) -> serde_json::Value {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len() as f64;
    let mean = s.iter().sum::<f64>() / n;
    let sd = (s.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0).max(1.0)).sqrt();
    let q = |p: f64| s[((p * (n - 1.0)).round() as usize).min(s.len() - 1)];
    json!({
        "n": s.len(), "min": s[0], "max": s[s.len() - 1], "mean": mean,
        "median": q(0.5), "p10": q(0.1), "p90": q(0.9), "sd": sd,
        "ci95_half": 1.96 * sd / n.sqrt(),
    })
}

// ---------------------------------------------------------------- ecbench export

/// Exponent n₁ of `E(F_p)` (the lcm of random points' orders); the group is
/// `Z/n₁ × Z/n₂` with `n₂ = #E / n₁`. A point's order divides `#E`; 24
/// random points give the exponent with overwhelming probability.
fn group_structure(f: &Fp, e: &Curve, order: u64, rng: &mut SmallRng) -> (u64, u64) {
    let primes: Vec<u64> = factor(order).into_iter().map(|(q, _)| q).collect();
    let mut expo: u128 = 1;
    for _ in 0..24 {
        let pt = e.random_point(f, rng);
        let mut ord = order;
        for &q in &primes {
            while ord % q == 0 && e.mul(f, Some(pt), ord / q).is_none() {
                ord /= q;
            }
        }
        let g = gcd(expo as u64, ord) as u128;
        expo = expo / g * ord as u128;
    }
    let n1 = expo as u64;
    (n1, order / n1)
}

/// ecbench input for every walked curve (`--ecbench`): see `--help`.
#[allow(clippy::too_many_arguments)]
fn write_ecbench(
    args: &Args,
    dir: &PathBuf,
    f: &Fp,
    vs: &[Vertex],
    levels: &[Option<u32>],
    walked: &[usize],
    order: u64,
    n: u64,
    rng: &mut SmallRng,
) -> (usize, std::collections::BTreeMap<String, usize>) {
    fs::create_dir_all(dir).unwrap();
    let mut curves = Vec::new();
    // (vertex, model, level) of each entry of `curves`.
    let mut owners: Vec<(usize, usize, u32)> = Vec::new();
    let mut side = String::new();
    let mut structures = std::collections::BTreeMap::new();
    for &v in walked {
        let e = vs[v].curve;
        let (n1, n2) = group_structure(f, &e, order, rng);
        *structures.entry(format!("{n1}x{n2}")).or_insert(0usize) += 1;
        let lv = levels[v].expect("walked curves have a level");
        let horizontal: Vec<usize> = vs[v].nbrs.iter().copied().filter(|&w| levels[w] == Some(lv)).collect();
        let models = if lv == 0 { args.crater_models.unwrap_or(args.models) } else { args.models }.max(1);
        for m in 0..models {
            // Model 0 is the walked model; others are (u⁴a, u⁶b) with the generator mapped.
            let u = if m == 0 { 1 } else { rng.gen_range(2..f.p) };
            let (u2, u3) = (f.mul(u, u), f.mul(f.mul(u, u), u));
            let cv = Curve {
                a: f.mul(f.mul(u2, u2), e.a),
                b: f.mul(f.mul(u3, u3), e.b),
            };
            let g = (f.mul(u2, vs[v].p.0), f.mul(u3, vs[v].p.1));
            assert!(cv.on(f, Some(g)));
            curves.push(json!({
                "kind": "prime_explicit", "p": f.p, "a": cv.a, "b": cv.b,
                "group_order": order, "r": n, "gx": g.0, "gy": g.1,
            }));
            owners.push((v, m, lv));
            side += &json!({
                "vertex": v, "model": m, "u": u, "level": lv, "j": vs[v].j,
                "p": f.p, "a": cv.a, "b": cv.b, "gx": g.0, "gy": g.1,
                "kernels": vs[v].kernels.len(), "n1": n1, "n2": n2,
                "dist_from_root": vs[v].dist, "parent": vs[v].parent,
                "nbrs": vs[v].nbrs, "horizontal_nbrs": horizontal,
            })
            .to_string();
            side.push('\n');
        }
    }
    let ic = |size: u32| json!({"id": "ic.pipeline", "params": {
        "factor_base": format!("prime-abscissa:size={size}"), "oracle": "mitm:negation_folded=1"}});
    let spec = |label: &str, arms: serde_json::Value, seed: u64| json!({
        "schema": "ecbench.spec/v1",
        "label": label,
        "description": format!(
            "volcano-rho --ecbench: every curve of an {}-isogeny volcano over F_{} ({} curves x {} models), \
             isogeny-transported generators, planted_shared targets; written by {}",
            args.ell, f.p, walked.len(), args.models.max(1),
            std::env::args().collect::<Vec<_>>().join(" ")),
        "workloads": {"curves": curves, "targets_per_curve": args.targets,
                      "target_seed": args.seed, "target_kind": "planted_shared"},
        "arms": arms,
        "measurement": {"rounds": 1, "warmup": 1, "order": "shuffle", "seed": seed,
                        "isolation_required": "L0", "timeout_seconds": 600},
    });
    let pick = |keep: &dyn Fn(&(usize, usize, u32)) -> bool| -> Vec<serde_json::Value> {
        curves.iter().zip(&owners).filter(|(_, o)| keep(o)).map(|(c, _)| c.clone()).collect()
    };
    let main = spec(
        "isogeny volcano: rho and IC on every curve",
        json!([
            {"name": "rho-neg", "role": "reference", "method": {"id": "rho.negation"}},
            {"name": "ic-mitm", "role": "candidate", "method": ic(args.ic_size)},
        ]),
        args.seed,
    );
    // BSGS's cost is fixed by its target: two paired targets check exactness.
    let mut bsgs = spec(
        "isogeny volcano: BSGS exactness",
        json!([{"name": "bsgs-neg", "role": "baseline", "method": {"id": "bsgs.negation"}}]),
        args.seed,
    );
    bsgs["workloads"]["targets_per_curve"] = json!(2);
    // The planted IC control (one more base point) on each curve's model 0,
    // same workloads as `spec.json`, so runs pair by workload.
    let mut icctl = spec(
        "isogeny volcano: planted IC control (factor base + 1)",
        json!([{"name": "ic-mitm-plus1", "role": "control", "method": ic(args.ic_size + 1)}]),
        args.seed,
    );
    icctl["workloads"]["curves"] = json!(pick(&|o| o.1 == 0));
    let aa = spec(
        "isogeny volcano: independent A/A rho",
        json!([{"name": "rho-neg-aa", "role": "control", "method": {"id": "rho.negation"}}]),
        args.seed ^ 0xA_A000_0000,
    );
    let mut aa = aa;
    if args.aa_model0 {
        aa["workloads"]["curves"] = json!(pick(&|o| o.1 == 0));
    }
    // Positive control on a j = 0 / j = 1728 crater: the automorphism-folded
    // walk (which refuses every other curve) beside rho.negation, crater only.
    let crater_j = walked.iter().find(|&&v| levels[v] == Some(0)).map(|&v| vs[v].j);
    let aut = match crater_j.unwrap_or(u64::MAX) {
        0 if f.p % 3 == 1 && (n - 1) % 3 == 0 => Some(("rho-cubic", "rho.cubic")),
        1728 if f.p % 4 == 1 && (n - 1) % 4 == 0 => Some(("rho-quartic", "rho.quartic")),
        _ => None,
    };
    if let Some((name, id)) = aut {
        let mut crater = spec(
            "isogeny volcano: automorphism positive control on the crater",
            json!([
                // Named apart from spec.json's rho-neg: same seeds, so pooling
                // the two would count each crater solve twice.
                {"name": "rho-neg-crater", "role": "reference", "method": {"id": "rho.negation"}},
                {"name": name, "role": "control", "method": {"id": id}},
            ]),
            args.seed,
        );
        let crater_curves = pick(&|o| o.2 == 0);
        crater["workloads"]["curves"] = json!(crater_curves);
        fs::write(dir.join("spec-crater-aut.json"), serde_json::to_string_pretty(&crater).unwrap()).unwrap();
    }
    fs::write(dir.join("spec.json"), serde_json::to_string_pretty(&main).unwrap()).unwrap();
    fs::write(dir.join("spec-aa.json"), serde_json::to_string_pretty(&aa).unwrap()).unwrap();
    fs::write(dir.join("spec-bsgs.json"), serde_json::to_string_pretty(&bsgs).unwrap()).unwrap();
    fs::write(dir.join("spec-icctl.json"), serde_json::to_string_pretty(&icctl).unwrap()).unwrap();
    fs::write(dir.join("curves.jsonl"), side).unwrap();
    (curves.len(), structures)
}

// ---------------------------------------------------------------- instance study

/// Two-sample Kolmogorov–Smirnov statistic and its asymptotic p-value.
fn ks(a: &[f64], b: &[f64]) -> (f64, f64) {
    let (mut x, mut y) = (a.to_vec(), b.to_vec());
    x.sort_by(|p, q| p.partial_cmp(q).unwrap());
    y.sort_by(|p, q| p.partial_cmp(q).unwrap());
    let (n1, n2) = (x.len() as f64, y.len() as f64);
    let (mut i, mut j, mut d) = (0, 0, 0f64);
    while i < x.len() && j < y.len() {
        let v = x[i].min(y[j]);
        while i < x.len() && x[i] <= v {
            i += 1;
        }
        while j < y.len() && y[j] <= v {
            j += 1;
        }
        d = d.max((i as f64 / n1 - j as f64 / n2).abs());
    }
    let lam = ((n1 * n2 / (n1 + n2)).sqrt() + 0.12 + 0.11 / (n1 * n2 / (n1 + n2)).sqrt()) * d;
    let p = (1..=100)
        .map(|j| {
            let j = j as f64;
            2.0 * (-1f64).powf(j - 1.0) * (-2.0 * j * j * lam * lam).exp()
        })
        .sum::<f64>()
        .clamp(0.0, 1.0);
    (d, p)
}

/// Pollard rho on one curve, five arms of `--trials` solves each. Trial t of
/// an arm draws what it varies from its own seed, so every arm replays alone.
#[allow(clippy::too_many_arguments)]
fn instance_study(
    args: &Args,
    f: &Fp,
    e0: Curve,
    n: u64,
    cof: u64,
    p0: (u64, u64),
    k0: u64,
    dp_bits: u32,
) {
    // (name, walk seed varies, k varies, generator varies, small k)
    let arms: [(&str, bool, bool, bool, bool); 5] = [
        ("walk", true, false, false, false),
        ("k", false, true, false, false),
        ("generator", false, false, true, false),
        ("all", true, true, true, false),
        ("small_k", false, true, false, true),
    ];
    let fixed_walk = args.seed.wrapping_mul(0x2545_F491_4F6C_DD1D) ^ 0x5eed;
    let lens = Lens::native();
    let mut out_arms = Vec::new();
    let mut rows = String::from("arm,trial,walk_seed,k,px,py,iters,restarts,ns\n");
    let mut samples: Vec<Vec<f64>> = Vec::new();
    let t0 = Instant::now();
    for (ai, &(name, vw, vk, vg, small)) in arms.iter().enumerate() {
        let res: Vec<(u64, u64, (u64, u64), Solve)> = (0..args.trials)
            .into_par_iter()
            .map(|t| {
                let mut r = SmallRng::seed_from_u64(
                    args.seed ^ ((ai as u64 + 1) << 48) ^ t.wrapping_mul(0x9E37_79B9_7F4A_7C15),
                );
                let walk = if vw { r.gen::<u64>() } else { fixed_walk };
                let k = match (vk, small) {
                    (true, true) => r.gen_range(1..=(1u64 << 16).min(n - 1)),
                    (true, false) => r.gen_range(1..n),
                    _ => k0,
                };
                let p = if vg {
                    loop {
                        if let Some(pt) = e0.mul(f, Some(e0.random_point(f, &mut r)), cof) {
                            break pt;
                        }
                    }
                } else {
                    p0
                };
                let q = e0.mul(f, Some(p), k).expect("1 <= k < n");
                let inst = Instance { curve: e0, p, q };
                let s = rho(f, &inst, n, walk, dp_bits, &lens);
                (walk, k, p, s)
            })
            .collect();
        let wrong = res.iter().filter(|r| r.3.k != Some(r.1)).count();
        if wrong > 0 {
            die(&format!("arm {name}: {wrong} solves returned the wrong k"));
        }
        for (t, (w, k, p, s)) in res.iter().enumerate() {
            rows += &format!(
                "{name},{t},{w},{k},{},{},{},{},{}\n",
                p.0, p.1, s.iters, s.restarts, s.ns
            );
        }
        let it: Vec<f64> = res.iter().map(|r| r.3.iters as f64).collect();
        let ns: Vec<f64> = res
            .iter()
            .map(|r| r.3.ns as f64 / r.3.iters as f64)
            .collect();
        let distinct_k = res.iter().map(|r| r.1).collect::<HashSet<_>>().len();
        let distinct_p = res.iter().map(|r| r.2).collect::<HashSet<_>>().len();
        let distinct_w = res.iter().map(|r| r.0).collect::<HashSet<_>>().len();
        eprintln!(
            "{name:>10}: mean {:>10.0} iters  ({} solves, {distinct_w} walk seeds, {distinct_k} k, {distinct_p} generators)",
            it.iter().sum::<f64>() / it.len() as f64,
            it.len()
        );
        out_arms.push(json!({
            "arm": name, "varies": {"walk_seed": vw, "k": vk, "generator": vg, "small_k": small},
            "distinct": {"walk_seeds": distinct_w, "k": distinct_k, "generators": distinct_p},
            "iters": stats(&it), "ns_per_iter": stats(&ns), "wrong": 0,
        }));
        samples.push(it);
    }
    let base = &samples[0];
    let bm = base.iter().sum::<f64>() / base.len() as f64;
    let bsd =
        (base.iter().map(|x| (x - bm).powi(2)).sum::<f64>() / (base.len() as f64 - 1.0)).sqrt();
    for (i, a) in out_arms.iter_mut().enumerate().skip(1) {
        let s = &samples[i];
        let m = s.iter().sum::<f64>() / s.len() as f64;
        let sd = (s.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (s.len() as f64 - 1.0)).sqrt();
        let se = (bsd * bsd / base.len() as f64 + sd * sd / s.len() as f64).sqrt();
        let (d, p) = ks(base, s);
        a["vs_walk"] = json!({
            "mean_ratio": m / bm, "diff": m - bm, "diff_ci95_half": 1.96 * se,
            "welch_z": (m - bm) / se, "ks_d": d, "ks_p": p,
        });
    }
    let expected = (std::f64::consts::PI * n as f64 / 2.0).sqrt();
    let result = json!({
        "schema": "volcano-rho.instance-study/v1",
        "curve": {"p": f.p, "a": e0.a, "b": e0.b, "order": n * cof, "n": n},
        "P": [p0.0, p0.1], "k": k0, "seed": args.seed, "trials_per_arm": args.trials,
        "fixed_walk_seed": fixed_walk, "dp_bits": dp_bits,
        "rho": {"walk": "native single-walker 32-adding DP, no negation map",
                "expected_iters": expected,
                "expected_iters_r32": expected / (1.0 - 1.0 / 32.0f64).sqrt()},
        "arms": out_arms, "seconds": t0.elapsed().as_secs_f64(),
        "argv": std::env::args().collect::<Vec<_>>(),
    });
    fs::create_dir_all(&args.out).unwrap();
    fs::write(
        args.out.join("instance_study.json"),
        serde_json::to_string_pretty(&result).unwrap(),
    )
    .unwrap();
    fs::write(args.out.join("instance_trials.csv"), rows).unwrap();
    eprintln!("wrote {}", args.out.display());
}

// ---------------------------------------------------------------- main

struct Args {
    bits: u32,
    ell: u64,
    depth: u32,
    seed: u64,
    trials: u64,
    per_level: usize,
    /// Distinguished-point bits; `None` picks `min(8, log2(n)/2 − 6)`.
    dp_bits: Option<u32>,
    max_cofactor: u64,
    /// `any`, `split`, `inert` or `ramified`: the crater shape the search requires.
    crater: String,
    /// Search j = 0 curves y² = x³ + b only (class D_K = −3, with the
    /// order-6 automorphism on its crater).
    j0: bool,
    /// Search j = 1728 curves y² = x³ + ax only (class D_K = −4, with the
    /// order-4 automorphism on its crater).
    j1728: bool,
    /// `p,a,b`: measure this curve instead of searching for one.
    curve: Option<(u64, u64, u64)>,
    /// `#E(F_p)` for `--curve`; counted when absent.
    order: Option<u64>,
    /// `native` (this binary's walk) or `ca` (`ca solve --alg rho`).
    engine: String,
    /// `native` (each curve hashes its own x) or `pullback` (every curve
    /// hashes the root's canonical point, so walks coincide across curves).
    walk: String,
    ca: String,
    ca_args: Vec<String>,
    controls: usize,
    /// Measure rho on the root curve only, redrawing the walk seed, k, the
    /// generator, or all three per trial (`instance_study`).
    instance_study: bool,
    /// Write ecbench specs for every walked curve to this directory
    /// instead of running rho (`write_ecbench`).
    ecbench: Option<PathBuf>,
    /// Isomorphic models per curve in the ecbench export (1 = as walked).
    models: usize,
    /// Models per crater curve, when crater curves need more replicates
    /// than the rest (`--crater-models`; default `models`).
    crater_models: Option<usize>,
    /// Curves kept per level, crater first (`--level-caps 16,4,12`); a
    /// level past the list uses `--per-level`.
    level_caps: Vec<usize>,
    /// Restrict the A/A spec to each curve's model 0.
    aa_model0: bool,
    /// Paired `planted_shared` targets per curve in the ecbench export.
    targets: u64,
    /// `prime-abscissa` factor-base size of the IC arm (the planted
    /// control arm uses one more).
    ic_size: u32,
    out: PathBuf,
}

fn parse_args() -> Args {
    let mut a = Args {
        bits: 40,
        ell: 3,
        depth: 3,
        seed: 1,
        trials: 100,
        per_level: 10,
        dp_bits: None,
        max_cofactor: 64,
        crater: "any".into(),
        j0: false,
        j1728: false,
        curve: None,
        order: None,
        engine: "native".into(),
        walk: "native".into(),
        ca: "ca".into(),
        ca_args: Vec::new(),
        controls: 3,
        instance_study: false,
        ecbench: None,
        models: 1,
        crater_models: None,
        level_caps: Vec::new(),
        aa_model0: false,
        targets: 64,
        ic_size: 32,
        out: PathBuf::from("volcano-rho-out"),
    };
    let v: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < v.len() {
        let val = v.get(i + 1).cloned().unwrap_or_default();
        match v[i].as_str() {
            "--bits" => a.bits = val.parse().unwrap(),
            "--ell" => a.ell = val.parse().unwrap(),
            "--depth" => a.depth = val.parse().unwrap(),
            "--seed" => a.seed = val.parse().unwrap(),
            "--trials" => a.trials = val.parse().unwrap(),
            "--per-level" => a.per_level = val.parse().unwrap(),
            "--dp-bits" => {
                a.dp_bits = Some(
                    val.parse()
                        .unwrap_or_else(|_| die("--dp-bits takes an integer")),
                )
            }
            "--max-cofactor" => a.max_cofactor = val.parse().unwrap(),
            "--crater" => a.crater = val.clone(),
            "--j1728" => {
                a.j1728 = true;
                i += 1; // a flag: no value
                continue;
            }
            "--ecbench" => a.ecbench = Some(PathBuf::from(val)),
            "--models" => a.models = val.parse().unwrap_or_else(|_| die("--models takes an integer")),
            "--crater-models" => {
                a.crater_models = Some(val.parse().unwrap_or_else(|_| die("--crater-models takes an integer")))
            }
            "--level-caps" => {
                a.level_caps = val
                    .split(',')
                    .map(|x| x.trim().parse().unwrap_or_else(|_| die("--level-caps takes N,N,...")))
                    .collect()
            }
            "--aa-model0" => {
                a.aa_model0 = true;
                i += 1; // a flag: no value
                continue;
            }
            "--targets" => a.targets = val.parse().unwrap_or_else(|_| die("--targets takes an integer")),
            "--ic-size" => a.ic_size = val.parse().unwrap_or_else(|_| die("--ic-size takes an integer")),
            "--instance-study" => {
                a.instance_study = true;
                i += 1; // a flag: no value
                continue;
            }
            "--j0" => {
                a.j0 = true;
                i += 1; // a flag: no value
                continue;
            }
            "--curve" => {
                let v: Vec<u64> = val
                    .split(',')
                    .map(|x| {
                        x.trim()
                            .parse()
                            .unwrap_or_else(|_| die("--curve takes P,A,B in decimal"))
                    })
                    .collect();
                if v.len() != 3 {
                    die("--curve takes P,A,B");
                }
                a.curve = Some((v[0], v[1], v[2]));
            }
            "--order" => {
                a.order = Some(
                    val.parse()
                        .unwrap_or_else(|_| die("--order takes a decimal integer")),
                )
            }
            "--engine" => a.engine = val.clone(),
            "--walk" => a.walk = val.clone(),
            "--ca" => a.ca = val.clone(),
            "--ca-arg" => a.ca_args.push(val.clone()),
            "--controls" => a.controls = val.parse().unwrap(),
            "--out" => a.out = PathBuf::from(val),
            "-h" | "--help" => {
                println!(
                    "volcano-rho: one ECDLP instance carried to every level of an ell-isogeny volcano,
Pollard rho measured on each curve.

  your curve:   volcano-rho --curve P,A,B [--order N] [options]
  search one:   volcano-rho [--bits 40] [--depth 3] [--crater any|split|inert|ramified] [--j0|--j1728]
                            [--max-cofactor 64] [options]

options: [--ell 3] [--seed 1] [--trials 100] [--per-level 10] [--controls 3]
         [--dp-bits auto] [--engine native|ca] [--walk native|pullback] [--ca PATH] [--ca-arg ARG]... [--out DIR]
         [--instance-study] [--ecbench DIR [--models 1] [--crater-models M] [--level-caps C0,C1,..]
                        [--aa-model0] [--targets 64] [--ic-size 32]]

--curve takes y^2 = x^3 + Ax + B over F_P, P prime < 2^62, in decimal. The
generator (cofactor times a seeded random point) and the secret k are drawn
from --seed. --engine ca runs `ca solve --alg rho` from the cryptanalysis
repository on every curve (--ca-arg passes extra options, e.g. --no-negation).
Writes results.json, trials.csv and instances.jsonl to --out.

--instance-study skips the walk and measures the root curve alone, --trials
solves per arm: the walk seed redrawn (P, k fixed), k redrawn, the generator
redrawn, all three redrawn, and k drawn from [1, 2^16]. Writes
instance_study.json and instance_trials.csv.

--ecbench DIR skips rho and writes ecbench input for up to --per-level
walked curves per level: spec.json (rho.negation and ic.pipeline at
--ic-size), spec-bsgs.json (bsgs.negation, 2 targets: an exactness check),
spec-icctl.json (ic.pipeline at --ic-size + 1 on model 0: the planted
control), spec-aa.json (rho.negation under another measurement seed: the
independent A/A control), spec-crater-aut.json on a j = 0 or 1728 crater
(rho.cubic / rho.quartic) and curves.jsonl (level, kernels, group
structure, model, neighbours of each workload curve). Each curve appears in
--models isomorphic models; every model carries the isogeny-transported
generator, and planted_shared targets give the same k everywhere."
                );
                std::process::exit(0);
            }
            o => die(&format!("unknown argument {o} (see --help)")),
        }
        i += 2;
    }
    if a.j0 && a.j1728 {
        die("--j0 and --j1728 exclude each other");
    }
    if a.ell != 2 && a.ell != 3 {
        die("--ell must be 2 or 3");
    }
    if !["any", "split", "inert", "ramified"].contains(&a.crater.as_str()) {
        die("--crater must be any, split, inert or ramified");
    }
    if !["native", "pullback"].contains(&a.walk.as_str()) {
        die("--walk must be native or pullback");
    }
    if a.walk == "pullback" && a.engine != "native" {
        die("--walk pullback needs --engine native");
    }
    if !["native", "ca"].contains(&a.engine.as_str()) {
        die("--engine must be native or ca");
    }
    if a.trials == 0 {
        die("--trials must be at least 1");
    }
    if a.curve.is_none() && (a.depth < 1 || !(20..=60).contains(&a.bits)) {
        die("search mode needs --depth >= 1 and 20 <= --bits <= 60");
    }
    a
}

fn vl(mut x: u128, l: u128) -> u32 {
    let mut v = 0;
    while x > 0 && x.is_multiple_of(l) {
        x /= l;
        v += 1;
    }
    v
}

struct Vertex {
    j: u64,
    curve: Curve,
    dist: u32,
    parent: Option<usize>,
    edge: Option<Velu>,
    kernels: Vec<u64>,
    nbrs: Vec<usize>,
    p: (u64, u64),
    q: (u64, u64),
}

struct Found {
    f: Fp,
    root: Curve,
    order: u64,
    n: u64,
    t: i64,
    d_abs: u64,
    depth: u32,
    tried: u64,
}

/// Random curves over random primes until one has an `ℓ`-volcano of the
/// requested depth and `#E = h·n` with `n` prime, `h ≤ max_cofactor`.
fn search(args: &Args, rng: &mut SmallRng) -> Found {
    let mut tried = 0u64;
    loop {
        let p = loop {
            let c = rng.gen_range((1u64 << (args.bits - 1))..(1u64 << args.bits));
            // p ≡ 3 mod 4 caps v₂(t² − 4p) at 3, so a 2-volcano deeper than
            // one level needs p ≡ 1 mod 8; ℓ = 3 keeps p ≡ 3 mod 4.
            let c = if args.ell == 2 { (c & !7) | 1 } else { c | 3 };
            // An ordinary j = 0 curve needs p ≡ 1 mod 3 (with the above: 1 mod 24
            // for ℓ = 2, 7 mod 12 for ℓ = 3).
            let c = if args.j0 {
                c - c % 24 + if args.ell == 2 { 1 } else { 7 }
            } else if args.j1728 {
                // An ordinary j = 1728 curve needs p ≡ 1 mod 4 (1 mod 8 for ℓ = 2 above).
                (c & !3) | 1
            } else {
                c
            };
            if is_prime(c) {
                break c;
            }
        };
        let f = Fp { p };
        for _ in 0..200 {
            tried += 1;
            let e = Curve {
                a: if args.j0 { 0 } else { rng.gen_range(1..p) },
                b: if args.j1728 { 0 } else { rng.gen_range(1..p) },
            };
            match e.j(&f) {
                Some(0) if args.j0 => {}
                Some(j) if args.j1728 && j == 1728 % p => {}
                Some(j) if !args.j0 && !args.j1728 && j != 0 && j != 1728 % p => {}
                _ => continue,
            }
            let Some(order) = group_order(&f, &e, rng) else {
                continue;
            };
            let t = p as i64 + 1 - order as i64;
            let d_abs = (4 * p as i128 - (t as i128) * (t as i128)) as u64;
            let (depth, kind) = ell_depth(d_abs, args.ell);
            if depth != args.depth {
                continue;
            }
            let n = factor(order).last().unwrap().0;
            if order / n > args.max_cofactor || n.is_multiple_of(args.ell) || n < (1 << 20) {
                continue;
            }
            if args.crater != "any" && !kind.starts_with(args.crater.as_str()) {
                continue;
            }
            return Found {
                f,
                root: e,
                order,
                n,
                t,
                d_abs,
                depth,
                tried,
            };
        }
    }
}

fn die(msg: &str) -> ! {
    eprintln!("volcano-rho: {msg}");
    std::process::exit(2)
}

/// The caller's curve `y² = x³ + ax + b` over `F_p`: count its points
/// unless `--order` gives `#E`, check the order, and take `n` the largest
/// prime factor.
fn from_curve(args: &Args, (p, a, b): (u64, u64, u64), rng: &mut SmallRng) -> Found {
    if !(5..1 << 62).contains(&p) || !is_prime(p) {
        die("--curve needs a prime 5 <= p < 2^62");
    }
    let f = Fp { p };
    let e = Curve { a: a % p, b: b % p };
    if e.j(&f).is_none() {
        die("singular curve: 4a^3 + 27b^2 = 0 mod p");
    }
    let order = match args.order {
        Some(o) => o,
        None => group_order(&f, &e, rng)
            .unwrap_or_else(|| die("point counting did not converge; pass --order")),
    };
    let t = p as i64 + 1 - order as i64;
    if (t as i128) * (t as i128) > 4 * p as i128 {
        die("--order is outside the Hasse interval");
    }
    for _ in 0..8 {
        if e.mul(&f, Some(e.random_point(&f, rng)), order).is_some() {
            die("--order is not #E(F_p): a random point does not vanish under it");
        }
    }
    let n = largest_prime_factor(order);
    if n.is_multiple_of(args.ell) {
        die(
            "ell divides the subgroup order n, so an ell-isogeny can kill <P>; use the other --ell",
        );
    }
    if t == 0 {
        eprintln!("volcano-rho: warning: t = 0, the curve is supersingular; levels may not follow the ordinary volcano");
    }
    let d_abs = (4 * p as i128 - (t as i128) * (t as i128)) as u64;
    Found {
        f,
        root: e,
        order,
        n,
        t,
        d_abs,
        depth: ell_depth(d_abs, args.ell).0,
        tried: 0,
    }
}

/// `(D_K, f)` with `−d_abs = f²·D_K`, `D_K` fundamental.
fn fundamental(d_abs: u64) -> (i64, u64) {
    let (mut core, mut f) = (1u64, 1u64);
    for (q, e) in factor(d_abs) {
        core *= q.pow(e % 2);
        f *= q.pow(e / 2);
    }
    let d0 = -(core as i64);
    if d0.rem_euclid(4) == 1 {
        (d0, f)
    } else {
        (4 * d0, f / 2)
    }
}

/// The `ℓ`-volcano's depth `v_ℓ(f)` and crater shape from `|D|` alone,
/// `D = t² − 4p = f²·D_K`, with no factorisation: the shape is the
/// Kronecker symbol `(D_K/ℓ)`, which `D / ℓ^{2·depth}` shares.
fn ell_depth(d_abs: u64, ell: u64) -> (u32, &'static str) {
    if ell == 2 {
        let e = d_abs.trailing_zeros();
        let m = (d_abs >> e) as i128;
        if e % 2 == 1 {
            return ((e - 3) / 2, crater_kind(0));
        }
        let d = -m; // D / 4^{e/2}
        if d.rem_euclid(4) == 1 {
            let r = match d.rem_euclid(8) {
                1 => 1,
                _ => -1,
            };
            (e / 2, crater_kind(r))
        } else {
            (e / 2 - 1, crater_kind(0))
        }
    } else {
        let e = vl(d_abs as u128, ell as u128);
        let d = -((d_abs / ell.pow(2 * (e / 2))) as i128);
        let r = if d.rem_euclid(ell as i128) == 0 {
            0
        } else if powmod(d.rem_euclid(ell as i128) as u64, (ell - 1) / 2, ell) == 1 {
            1
        } else {
            -1
        };
        (e / 2, crater_kind(r))
    }
}

fn crater_kind(r: i32) -> &'static str {
    match r {
        1 => "split (crater is a cycle; two horizontal edges per crater vertex)",
        -1 => "inert (crater is one vertex; no horizontal edges)",
        _ => "ramified (one horizontal edge per crater vertex)",
    }
}

fn main() {
    let args = parse_args();
    let mut rng = SmallRng::seed_from_u64(args.seed);
    let ell = args.ell;
    let t_search = Instant::now();
    let found = match args.curve {
        Some(c) => from_curve(&args, c, &mut rng),
        None => search(&args, &mut rng),
    };
    let search_s = t_search.elapsed().as_secs_f64();
    let (f, n, depth) = (found.f, found.n, found.depth);
    let crater = ell_depth(found.d_abs, ell).1;
    // Keep 1/θ well below the √n-step walk, so a cycle always holds a DP.
    let dp_bits = args
        .dp_bits
        .unwrap_or(((63 - n.leading_zeros()) / 2).saturating_sub(6).min(8));
    // Full D_K needs |D| factored; skip it when trial division would be slow.
    let (dk, cond) = if found.d_abs < 1 << 52 {
        let (dk, c) = fundamental(found.d_abs);
        (Some(dk), Some(c))
    } else {
        (None, None)
    };
    let field = Field::new(&BigUint::from(f.p)).expect("p prime");
    eprintln!(
        "curve: p={} a={} b={} #E={} n={} h={} t={} D=-{} depth={} crater: {} ({} curves tried, {:.1}s)",
        f.p,
        found.root.a,
        found.root.b,
        found.order,
        n,
        found.order / n,
        found.t,
        found.d_abs,
        depth,
        crater,
        found.tried,
        search_s
    );

    // The instance: P of order n, Q = kP.
    let e0 = found.root;
    let cof = found.order / n;
    let p0 = loop {
        if let Some(pt) = e0.mul(&f, Some(e0.random_point(&f, &mut rng)), cof) {
            break pt;
        }
    };
    assert!(e0.mul(&f, Some(p0), n).is_none());
    let k_secret = rng.gen_range(1..n);
    let q0 = e0.mul(&f, Some(p0), k_secret).unwrap();
    if args.instance_study {
        instance_study(&args, &f, e0, n, cof, p0, k_secret, dp_bits);
        return;
    }

    // 2. Breadth-first walk with explicit Vélu maps.
    let r_meas = 2 * depth + 2;
    let r_explore = r_meas + depth;
    let mut vs: Vec<Vertex> = vec![Vertex {
        j: e0.j(&f).unwrap(),
        curve: e0,
        dist: 0,
        parent: None,
        edge: None,
        kernels: Vec::new(),
        nbrs: Vec::new(),
        p: p0,
        q: q0,
    }];
    let mut by_j: HashMap<u64, usize> = HashMap::from([(vs[0].j, 0)]);
    let mut queue = VecDeque::from([0usize]);
    let mut edges_checked = 0u64;
    let mut skipped_special = 0u64;
    while let Some(v) = queue.pop_front() {
        let e = vs[v].curve;
        let kernels = kernel_xs(&f, &field, &e, ell);
        vs[v].kernels = kernels.clone();
        for x0 in kernels {
            let iso = velu(&f, &e, ell, x0);
            let c = iso.codomain;
            let Some(j2) = c.j(&f) else { continue };
            if j2 == 0 || j2 == 1728 % f.p {
                skipped_special += 1;
                continue;
            }
            let w = match by_j.get(&j2) {
                Some(&w) => w,
                None => {
                    if vs[v].dist >= r_explore || vs.len() >= 20_000 {
                        continue;
                    }
                    let (pp, qq) = (Some(vs[v].p), Some(vs[v].q));
                    let (ip, iq) = (iso.map(&f, pp), iso.map(&f, qq));
                    let isum = iso.map(&f, e.add(&f, pp, qq));
                    assert!(c.on(&f, ip) && c.on(&f, iq), "image off the codomain");
                    assert_eq!(c.add(&f, ip, iq), isum, "Vélu map is not a homomorphism");
                    assert!(
                        ip.is_some() && c.mul(&f, ip, n).is_none(),
                        "order not preserved"
                    );
                    edges_checked += 1;
                    let w = vs.len();
                    vs.push(Vertex {
                        j: j2,
                        curve: c,
                        dist: vs[v].dist + 1,
                        parent: Some(v),
                        edge: Some(iso),
                        kernels: Vec::new(),
                        nbrs: Vec::new(),
                        p: ip.unwrap(),
                        q: iq.unwrap(),
                    });
                    by_j.insert(j2, w);
                    queue.push_back(w);
                    w
                }
            };
            if w != v && !vs[v].nbrs.contains(&w) {
                vs[v].nbrs.push(w);
            }
            if w != v && !vs[w].nbrs.contains(&v) {
                vs[w].nbrs.push(v);
            }
        }
    }

    // Levels: depth − distance to the floor (vertices with one ℓ-kernel).
    let level_of = |v: usize| -> Option<u32> {
        if depth == 0 {
            return Some(0); // no ℓ in the conductor: the whole component is crater
        }
        let mut seen = HashMap::from([(v, 0u32)]);
        let mut q = VecDeque::from([v]);
        while let Some(u) = q.pop_front() {
            let du = seen[&u];
            if vs[u].kernels.len() == 1 {
                return Some(depth.saturating_sub(du));
            }
            if du == depth {
                continue;
            }
            for &w in &vs[u].nbrs {
                if let std::collections::hash_map::Entry::Vacant(s) = seen.entry(w) {
                    s.insert(du + 1);
                    q.push_back(w);
                }
            }
        }
        None
    };
    let levels: Vec<Option<u32>> = (0..vs.len())
        .map(|v| (vs[v].dist <= r_meas).then(|| level_of(v)).flatten())
        .collect();
    for (v, lv) in levels.iter().enumerate() {
        let Some(lv) = lv else { continue };
        for &w in &vs[v].nbrs {
            if let Some(lw) = levels[w] {
                assert!(lv.abs_diff(lw) <= 1, "edge skips a level");
            }
        }
    }
    let mut per_level: Vec<Vec<usize>> = vec![Vec::new(); depth as usize + 1];
    for (v, lv) in levels.iter().enumerate() {
        if let Some(lv) = lv {
            per_level[*lv as usize].push(v);
        }
    }
    eprintln!(
        "walk: {} curves, {} verified edges; per level (0 = crater): {:?}",
        vs.len(),
        edges_checked,
        per_level.iter().map(|l| l.len()).collect::<Vec<_>>()
    );

    // 3. Curves to measure: the root, then up to `per_level` per level in BFS order.
    let mut measured: Vec<usize> = vec![0];
    for l in &per_level {
        for &v in l.iter().filter(|&&v| v != 0).take(args.per_level) {
            measured.push(v);
        }
    }
    measured.sort();
    struct Arm {
        label: String,
        vertex: Option<usize>,
        level: Option<u32>,
        inst: Instance,
        lens: Lens,
    }
    let pullback = args.walk == "pullback";
    // Each measured vertex's dual edge back to its parent: the rational
    // kernel whose codomain has the parent's j, the model isomorphism onto
    // the parent's model, and a check that φ̂(φ(P)) = ±[ℓ]P.
    let lift_to = measured.iter().map(|&v| vs[v].dist).max().unwrap_or(0);
    let mut duals: HashMap<usize, Pull> = HashMap::new();
    if pullback {
        for &v0 in &measured {
            let mut v = v0;
            while let Some(par) = vs[v].parent {
                if duals.contains_key(&v) {
                    break;
                }
                let (ev, ep) = (vs[v].curve, vs[par].curve);
                let want = ep.mul(&f, Some(vs[par].p), ell).expect("[l]P is not O").0;
                let found = kernel_xs(&f, &field, &ev, ell).into_iter().find_map(|x0| {
                    let d = velu(&f, &ev, ell, x0);
                    let c = d.codomain;
                    if c.j(&f) != Some(vs[par].j) {
                        return None;
                    }
                    let u2 = f.mul(f.mul(ep.b, c.a), f.inv(f.mul(c.b, ep.a)));
                    let u = f.sqrt(u2)?;
                    let q = d.map(&f, Some(vs[v].p))?;
                    (f.mul(u2, q.0) == want).then_some(Pull::Dual(d, u))
                });
                duals.insert(
                    v,
                    found.unwrap_or_else(|| {
                        die("internal: no dual isogeny found for a volcano edge")
                    }),
                );
                v = par;
            }
        }
    }
    let lens_for = |v: usize| -> Lens {
        if !pullback {
            return Lens::native();
        }
        let mut chain = Vec::new();
        let mut cur = v;
        while vs[cur].parent.is_some() {
            chain.push(duals[&cur]);
            cur = vs[cur].parent.unwrap();
        }
        Lens {
            chain,
            root: e0,
            lift: ell.pow(lift_to - vs[v].dist),
            native: false,
        }
    };
    let mut arms: Vec<Arm> = measured
        .iter()
        .map(|&v| Arm {
            label: format!("v{v}"),
            vertex: Some(v),
            level: levels[v],
            inst: Instance {
                curve: vs[v].curve,
                p: vs[v].p,
                q: vs[v].q,
            },
            lens: lens_for(v),
        })
        .collect();
    for c in 0..args.controls {
        // Isomorphic model (u⁴a, u⁶b), (x, y) ↦ (u²x, u³y): same curve, new coordinates.
        let u = rng.gen_range(2..f.p);
        let (u2, u3) = (f.mul(u, u), f.mul(f.mul(u, u), u));
        let cv = Curve {
            a: f.mul(f.mul(u2, u2), e0.a),
            b: f.mul(f.mul(u3, u3), e0.b),
        };
        let mv = |pt: (u64, u64)| (f.mul(u2, pt.0), f.mul(u3, pt.1));
        let (cp, cq) = (mv(p0), mv(q0));
        assert!(cv.on(&f, Some(cp)) && cv.on(&f, Some(cq)));
        arms.push(Arm {
            label: format!("iso{c}"),
            vertex: None,
            level: levels[0],
            inst: Instance {
                curve: cv,
                p: cp,
                q: cq,
            },
            lens: if pullback {
                Lens {
                    chain: vec![Pull::Scale(f.inv(u))],
                    root: e0,
                    lift: ell.pow(lift_to),
                    native: false,
                }
            } else {
                Lens::native()
            },
        });
    }

    // Every arm's instance in `ca`'s format, for any external solver.
    fs::create_dir_all(&args.out).unwrap();
    let mut jl = String::new();
    for a in &arms {
        jl += &json!({
            "arm": a.label, "vertex": a.vertex, "level": a.level, "group": "ec",
            "p": f.p, "a": a.inst.curve.a, "b": a.inst.curve.b, "order": n,
            "g": format!("{},{}", a.inst.p.0, a.inst.p.1),
            "h": format!("{},{}", a.inst.q.0, a.inst.q.1), "x": k_secret,
        })
        .to_string();
        jl.push('\n');
    }
    fs::write(args.out.join("instances.jsonl"), jl).unwrap();
    if let Some(dir) = &args.ecbench {
        // Up to --per-level curves per level, in BFS order; the root always.
        let mut taken = vec![0usize; depth as usize + 1];
        let walked: Vec<usize> = (0..vs.len())
            .filter(|&v| {
                let Some(l) = levels[v] else { return false };
                let cap = args.level_caps.get(l as usize).copied().unwrap_or(args.per_level);
                let keep = v == 0 || taken[l as usize] < cap;
                if keep {
                    taken[l as usize] += 1;
                }
                keep
            })
            .collect();
        let (curves, n1n2) = write_ecbench(&args, dir, &f, &vs, &levels, &walked, found.order, n, &mut rng);
        eprintln!(
            "ecbench: {} curves x {} models = {} workload curves, {} targets each; group structures {:?} -> {}",
            walked.len(),
            args.models.max(1),
            curves,
            args.targets,
            n1n2,
            dir.display()
        );
        return;
    }

    // 4. Pollard rho: every arm × the same trial seeds, interleaved.
    let mut jobs: Vec<(usize, u64)> = (0..arms.len())
        .flat_map(|a| (0..args.trials).map(move |s| (a, s)))
        .collect();
    jobs.shuffle(&mut rng);
    let seed_base = args.seed.wrapping_mul(0x1000_0000_01b3);
    eprintln!(
        "rho: {} arms x {} trials = {} solves",
        arms.len(),
        args.trials,
        jobs.len()
    );
    let t_rho = Instant::now();
    let results: Vec<(usize, u64, Solve)> = jobs
        .par_iter()
        .map(|&(a, s)| {
            let inst = &arms[a].inst;
            let r = if args.engine == "ca" {
                rho_ca(&args.ca, &args.ca_args, &f, inst, n, seed_base ^ s)
            } else {
                rho(&f, inst, n, seed_base ^ s, dp_bits, &arms[a].lens)
            };
            (a, s, r)
        })
        .collect();
    let rho_s = t_rho.elapsed().as_secs_f64();
    let wrong = results.iter().filter(|r| r.2.k != Some(k_secret)).count();
    if wrong > 0 {
        die(&format!(
            "{wrong} rho solves returned no or the wrong logarithm"
        ));
    }

    // In pullback mode every curve should take exactly the same walk.
    let mut seed_iters: HashMap<u64, Vec<u64>> = HashMap::new();
    for (_, s, r) in &results {
        seed_iters.entry(*s).or_default().push(r.iters);
    }
    let identical_seeds = seed_iters
        .values()
        .filter(|v| v.iter().all(|&x| x == v[0]))
        .count();
    if pullback {
        eprintln!(
            "pullback: {identical_seeds}/{} seeds took the same number of steps on every curve",
            seed_iters.len()
        );
    }

    // 5. Records and summaries.
    fs::create_dir_all(&args.out).unwrap();
    let mut csv = String::from("arm,vertex,level,seed,iters,restarts,ns\n");
    let mut by_arm: Vec<Vec<&Solve>> = vec![Vec::new(); arms.len()];
    let mut sorted: Vec<&(usize, u64, Solve)> = results.iter().collect();
    sorted.sort_by_key(|r| (r.0, r.1));
    for (a, s, r) in &sorted {
        by_arm[*a].push(r);
        csv += &format!(
            "{},{},{},{},{},{},{}\n",
            arms[*a].label,
            arms[*a].vertex.map_or(String::new(), |v| v.to_string()),
            arms[*a].level.map_or(String::new(), |l| l.to_string()),
            seed_base ^ s,
            r.iters,
            r.restarts,
            r.ns
        );
    }
    fs::write(args.out.join("trials.csv"), csv).unwrap();

    let summarize = |rs: &[&Solve]| {
        let it: Vec<f64> = rs.iter().map(|r| r.iters as f64).collect();
        let ms: Vec<f64> = rs.iter().map(|r| r.ns as f64 / 1e6).collect();
        let nsop: Vec<f64> = rs.iter().map(|r| r.ns as f64 / r.iters as f64).collect();
        json!({"iters": stats(&it), "ms": stats(&ms), "ns_per_op": stats(&nsop)})
    };
    let arm_json: Vec<serde_json::Value> = arms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            json!({
                "arm": a.label, "vertex": a.vertex, "level": a.level,
                "kind": if a.vertex.is_some() { "isogenous" } else { "isomorphic-control" },
                "curve": {"a": a.inst.curve.a, "b": a.inst.curve.b, "j": a.inst.curve.j(&f)},
                "P": [a.inst.p.0, a.inst.p.1], "Q": [a.inst.q.0, a.inst.q.1],
                "summary": summarize(&by_arm[i]),
            })
        })
        .collect();
    let mut level_json = Vec::new();
    for l in 0..=depth {
        let rs: Vec<&Solve> = arms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.vertex.is_some() && a.level == Some(l))
            .flat_map(|(i, _)| by_arm[i].iter().copied())
            .collect();
        if !rs.is_empty() {
            let curves = arms
                .iter()
                .filter(|a| a.vertex.is_some() && a.level == Some(l))
                .count();
            level_json.push(json!({"level": l, "curves": curves, "summary": summarize(&rs)}));
        }
    }
    let all_iso: Vec<&Solve> = arms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.vertex.is_some())
        .flat_map(|(i, _)| by_arm[i].iter().copied())
        .collect();
    let ctrl: Vec<&Solve> = arms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.vertex.is_none())
        .flat_map(|(i, _)| by_arm[i].iter().copied())
        .collect();
    let vertices: Vec<serde_json::Value> = vs
        .iter()
        .enumerate()
        .map(|(i, v)| {
            json!({
                "id": i, "j": v.j, "a": v.curve.a, "b": v.curve.b, "dist": v.dist,
                "level": levels[i], "kernels": v.kernels.len(), "nbrs": v.nbrs,
                "parent": v.parent,
                "edge": v.edge.map(|e| json!({"ell": e.ell, "kernel_x": e.x0})),
            })
        })
        .collect();
    let expected = (std::f64::consts::PI * n as f64 / 2.0).sqrt();
    let host = fs::read_to_string("/proc/cpuinfo").ok().and_then(|c| {
        c.lines()
            .find(|l| l.starts_with("model name"))
            .map(|l| l.split(':').nth(1).unwrap().trim().to_string())
    });
    let out = json!({
        "schema": "volcano-rho/v1",
        "command": std::env::args().collect::<Vec<_>>().join(" "),
        "params": {"bits": args.bits, "ell": ell, "depth": args.depth, "seed": args.seed, "crater": args.crater,
                   "curve": args.curve.map(|(p, a, b)| [p, a, b]), "order": args.order,
                   "engine": args.engine, "walk": args.walk, "ca": (args.engine == "ca").then_some(&args.ca), "ca_args": args.ca_args,
                   "trials": args.trials, "per_level": args.per_level, "dp_bits": dp_bits,
                   "max_cofactor": args.max_cofactor, "controls": args.controls},
        "host": {"cpu": host, "threads": rayon::current_num_threads()},
        "problem": {"p": f.p, "root": {"a": e0.a, "b": e0.b}, "order": found.order, "n": n,
                    "cofactor": found.order / n, "trace": found.t, "disc": -(found.d_abs as i64),
                    "fundamental_disc": dk, "conductor": cond, "volcano_depth": depth,
                    "crater": crater, "given": args.curve.is_some(), "P": [p0.0, p0.1], "Q": [q0.0, q0.1], "k": k_secret,
                    "search_curves_tried": found.tried, "search_s": search_s},
        "walk": {"radius_measured": r_meas, "radius_explored": r_explore, "vertices": vs.len(),
                 "edges_verified": edges_checked, "special_j_skipped": skipped_special,
                 "per_level_found": per_level.iter().map(|l| l.len()).collect::<Vec<_>>()},
        "pullback": {"enabled": pullback, "lift_to_dist": lift_to,
                     "seeds_identical_on_every_curve": identical_seeds, "seeds": seed_iters.len()},
        "rho": {"walk": if args.engine == "ca" {
                    "ca solve --alg rho (cryptanalysis repository; its defaults plus ca_args); iters = its reported ops, time = its reported seconds"
                } else {
                    "32-adding (Teske), single walker, distinguished points, no negation map"
                },
                "expected_iters_sqrt_pi_n_over_2": expected,
                "expected_iters_negation_sqrt_pi_n_over_4": (std::f64::consts::PI * n as f64 / 4.0).sqrt(), "solves": results.len(),
                "all_correct": wrong == 0, "wall_s": rho_s},
        "arms": arm_json,
        "levels": level_json,
        "all_isogenous": summarize(&all_iso),
        "controls": if ctrl.is_empty() { serde_json::Value::Null } else { summarize(&ctrl) },
        "vertices": vertices,
    });
    fs::write(
        args.out.join("results.json"),
        serde_json::to_string_pretty(&out).unwrap(),
    )
    .unwrap();

    println!("\nPollard rho on one instance (k = {k_secret}, n = {n}) across the {ell}-volcano");
    println!(
        "engine {}: expected iterations ~ sqrt(pi n / 2) = {expected:.0}, with a negation map sqrt(pi n / 4) = {:.0}",
        args.engine,
        (std::f64::consts::PI * n as f64 / 4.0).sqrt()
    );
    println!(
        "{:<8}{:>6}{:>12}{:>12}{:>12}{:>10}{:>10}{:>10}{:>9}",
        "arm", "level", "min it", "mean it", "max it", "min ms", "mean ms", "max ms", "ns/op"
    );
    let row = |name: &str, lv: String, s: &serde_json::Value| {
        println!(
            "{:<8}{:>6}{:>12.0}{:>12.0}{:>12.0}{:>10.2}{:>10.2}{:>10.2}{:>9.1}",
            name,
            lv,
            s["iters"]["min"].as_f64().unwrap(),
            s["iters"]["mean"].as_f64().unwrap(),
            s["iters"]["max"].as_f64().unwrap(),
            s["ms"]["min"].as_f64().unwrap(),
            s["ms"]["mean"].as_f64().unwrap(),
            s["ms"]["max"].as_f64().unwrap(),
            s["ns_per_op"]["mean"].as_f64().unwrap()
        );
    };
    for a in &arm_json {
        row(
            a["arm"].as_str().unwrap(),
            a["level"].to_string(),
            &a["summary"],
        );
    }
    println!("-- pooled by level --");
    for l in &level_json {
        row(
            &format!("L{}", l["level"]),
            format!("{}c", l["curves"]),
            &l["summary"],
        );
    }
    row("ALL", "iso".into(), &out["all_isogenous"]);
    if !ctrl.is_empty() {
        row("CTRL", "iso~".into(), &out["controls"]);
    }
    eprintln!("wrote {}", args.out.display());
}
