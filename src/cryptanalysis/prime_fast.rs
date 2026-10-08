//! Native-width prime-field ECDLP engine — the speed-up prototype for
//! the prime regime of the index-calculus ladder.
//!
//! Every prime-field rho and index-calculus path in this crate runs on
//! heap-allocated `BigUint` through [`crate::ecc::field::FieldElement`],
//! whose inverse is a Fermat exponentiation on a fixed-iteration ladder
//! (about `2·log₂p` big-integer multiplications per inverse) and whose
//! affine point addition pays one such inverse per step.  This module
//! replaces that stack for `p < 2^62` with:
//!
//! 1. **Montgomery arithmetic on `u64`** ([`Fp64`]): one `u128` product
//!    and one reduction per multiplication, no allocation, and a
//!    binary extended-Euclid inverse.
//! 2. **Simultaneous inversion** (Montgomery's trick, [`Fp64::batch_inv`]):
//!    `W` affine additions share one inverse at the price of `3(W−1)`
//!    multiplications.
//! 3. **Parallel distinguished-point Pollard rho** ([`rho_parallel`]):
//!    `r`-adding walk with a 1024-entry table, the negation map with a
//!    deterministic fruitless-2-cycle escape, lock-step walkers per
//!    thread so every step is batched, and a shared DP table.
//! 4. **Index-calculus relation collection without summation
//!    polynomials** ([`collect_relations`]): the Semaev `S₃` root
//!    `X` of `S₃(x_R, x_i, X) = 0` *is* `x(R ± F_i)`, so a
//!    2-decomposition sweep is one batched subtraction per factor-base
//!    element plus an array lookup — no square root, no inversion per
//!    element.  The `S₄` quartic roots are likewise `x(R ∓ F_i)` matched
//!    against a precomputed table of `x(F_j ± F_k)` ([`PairTable`]),
//!    turning the `B²/2` root-findings per target into `2B` lookups.
//! 5. **Linear algebra on `u64`** mod the group order (Montgomery form,
//!    Gauss–Jordan) and a **Hasse-interval BSGS point counter**
//!    ([`point_order_hasse`]) so the `a = −3` bench ladder can be
//!    extended past the `u64` brute-force counter's reach.
//!
//! The relation sets produced by (4) are *identical* to the ones the
//! `S₃` / `S₄` sweeps in [`crate::cryptanalysis::ec_index_calculus`]
//! produce for the same factor base (the summation polynomial vanishes
//! exactly when the x-coordinates come from points summing to zero), so
//! the staged timers remain comparable; only the cost per target
//! changes.  Everything is public-synthetic known-answer only, and every
//! recovered logarithm is re-checked in the group before it is returned.
//!
//! The benchmark driver is `examples/prime_ecdlp_fast_bench.rs`.

use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Mutex;
use std::time::Instant;

use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::cryptanalysis::residual_walk::is_prime_u64;
use crate::ecc::curve::CurveParams;
use crate::ecc::field::FieldElement;
use crate::ecc::point::Point;

/// Largest field size the Montgomery layer supports (`p < 2^62` keeps
/// `a + b` and `t + m·p` inside `u64` / `u128`).
pub const MAX_BITS: u32 = 62;

// ── Montgomery F_p on u64 ─────────────────────────────────────────────

/// Montgomery field context for an odd prime `p < 2^62`, `R = 2^64`.
/// Values are held as `a·R mod p`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fp64 {
    pub p: u64,
    /// `−p⁻¹ mod 2^64`.
    n_inv: u64,
    /// `R mod p` — the Montgomery form of 1.
    pub one: u64,
    /// `R² mod p`.
    r2: u64,
    /// `R³ mod p` (turns a canonical inverse of a Montgomery value back
    /// into Montgomery form in one multiplication).
    r3: u64,
    /// `(p − 1) / 2`, the sign boundary for the negation map.
    pub half: u64,
}

impl Fp64 {
    pub fn new(p: u64) -> Self {
        assert!(p & 1 == 1 && p > 3, "Fp64 needs an odd prime > 3");
        assert!(p < (1u64 << MAX_BITS), "Fp64 supports p < 2^{MAX_BITS}");
        // Newton iteration for p⁻¹ mod 2^64 (each step doubles the
        // number of correct bits; 6 steps from 1 correct bit is enough).
        let mut inv: u64 = 1;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        }
        let n_inv = inv.wrapping_neg();
        let r = ((1u128 << 64) % p as u128) as u64;
        let r2 = ((r as u128 * r as u128) % p as u128) as u64;
        let mut f = Fp64 {
            p,
            n_inv,
            one: r,
            r2,
            r3: 0,
            half: (p - 1) / 2,
        };
        f.r3 = f.mul(r2, r2);
        f
    }

    /// Montgomery reduction of `t < p·2^64`: returns `t·R⁻¹ mod p`.
    #[inline(always)]
    pub fn reduce(&self, t: u128) -> u64 {
        let m = (t as u64).wrapping_mul(self.n_inv);
        let u = ((t + (m as u128) * (self.p as u128)) >> 64) as u64;
        if u >= self.p {
            u - self.p
        } else {
            u
        }
    }

    #[inline(always)]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        self.reduce(a as u128 * b as u128)
    }

    #[inline(always)]
    pub fn sqr(&self, a: u64) -> u64 {
        self.mul(a, a)
    }

    #[inline(always)]
    pub fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }

    #[inline(always)]
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }

    #[inline(always)]
    pub fn neg(&self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }

    /// Canonical → Montgomery.
    #[inline]
    pub fn to_mont(&self, a: u64) -> u64 {
        self.mul(a % self.p, self.r2)
    }

    /// Montgomery → canonical.
    #[inline(always)]
    pub fn from_mont(&self, a: u64) -> u64 {
        self.reduce(a as u128)
    }

    /// Inverse of a non-zero Montgomery value, in Montgomery form.
    #[inline]
    pub fn inv(&self, a: u64) -> u64 {
        // inv_mod gives (a·R)⁻¹ = a⁻¹·R⁻¹; times R³ through one
        // Montgomery multiplication (which divides by R) is a⁻¹·R.
        let c = inv_mod_u64(a, self.p);
        self.mul(c, self.r3)
    }

    /// `a^e` for a Montgomery value `a` (result in Montgomery form).
    pub fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut acc = self.one;
        while e > 0 {
            if e & 1 == 1 {
                acc = self.mul(acc, a);
            }
            a = self.sqr(a);
            e >>= 1;
        }
        acc
    }

    /// Montgomery's trick: invert every (non-zero) element of `xs` in
    /// place with a single field inversion and `3(k − 1)` multiplications.
    pub fn batch_inv(&self, xs: &mut [u64], scratch: &mut Vec<u64>) {
        let k = xs.len();
        if k == 0 {
            return;
        }
        scratch.clear();
        scratch.reserve(k);
        let mut acc = self.one;
        for &x in xs.iter() {
            scratch.push(acc);
            acc = self.mul(acc, x);
        }
        let mut inv_acc = self.inv(acc);
        for i in (0..k).rev() {
            let xi = xs[i];
            xs[i] = self.mul(inv_acc, scratch[i]);
            inv_acc = self.mul(inv_acc, xi);
        }
    }

    /// Square root of a Montgomery value (Tonelli–Shanks), in Montgomery
    /// form.  `None` if `a` is a non-residue.
    pub fn sqrt(&self, a: u64) -> Option<u64> {
        if a == 0 {
            return Some(0);
        }
        let p = self.p;
        if self.pow(a, (p - 1) / 2) != self.one {
            return None;
        }
        if p % 4 == 3 {
            return Some(self.pow(a, (p + 1) / 4));
        }
        let mut q = p - 1;
        let mut s = 0u32;
        while q & 1 == 0 {
            q >>= 1;
            s += 1;
        }
        let minus_one = self.neg(self.one);
        let mut z = 2u64;
        while self.pow(self.to_mont(z), (p - 1) / 2) != minus_one {
            z += 1;
        }
        let mut m = s;
        let mut c = self.pow(self.to_mont(z), q);
        let mut t = self.pow(a, q);
        let mut r = self.pow(a, (q + 1) / 2);
        loop {
            if t == self.one {
                return Some(r);
            }
            let mut i = 0u32;
            let mut tt = t;
            while tt != self.one {
                tt = self.sqr(tt);
                i += 1;
                if i == m {
                    return None;
                }
            }
            let mut b = c;
            for _ in 0..(m - i - 1) {
                b = self.sqr(b);
            }
            m = i;
            c = self.sqr(b);
            t = self.mul(t, c);
            r = self.mul(r, b);
        }
    }
}

/// `a⁻¹ mod p` for odd `p < 2^63` and `0 < a < p`, by the binary
/// extended Euclidean algorithm (HAC 14.61): no divisions.
pub fn inv_mod_u64(a: u64, p: u64) -> u64 {
    debug_assert!(a != 0 && a < p && p & 1 == 1);
    let (mut u, mut v) = (a, p);
    // Invariants: x1·a ≡ u, x2·a ≡ v (mod p).
    let (mut x1, mut x2) = (1u64, 0u64);
    #[inline(always)]
    fn halve(x: u64, p: u64) -> u64 {
        if x & 1 == 0 {
            x >> 1
        } else {
            ((x as u128 + p as u128) >> 1) as u64
        }
    }
    while u != 1 && v != 1 {
        while u & 1 == 0 {
            u >>= 1;
            x1 = halve(x1, p);
        }
        while v & 1 == 0 {
            v >>= 1;
            x2 = halve(x2, p);
        }
        if u >= v {
            u -= v;
            x1 = if x1 >= x2 { x1 - x2 } else { x1 + p - x2 };
        } else {
            v -= u;
            x2 = if x2 >= x1 { x2 - x1 } else { x2 + p - x1 };
        }
    }
    if u == 1 {
        x1
    } else {
        x2
    }
}

#[inline(always)]
fn add_mod(a: u64, b: u64, n: u64) -> u64 {
    let s = a + b;
    if s >= n {
        s - n
    } else {
        s
    }
}

#[inline(always)]
fn neg_mod(a: u64, n: u64) -> u64 {
    if a == 0 {
        0
    } else {
        n - a
    }
}

// ── Curve and points ──────────────────────────────────────────────────

/// Finite affine point with Montgomery-form coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Pt {
    pub x: u64,
    pub y: u64,
}

/// Short-Weierstrass curve `y² = x³ + ax + b` over `F_p`, `p < 2^62`,
/// with a prime-order generator subgroup of order `n`.
#[derive(Clone, Debug)]
pub struct FastCurve {
    pub f: Fp64,
    /// `a`, `b` in Montgomery form.
    pub a: u64,
    pub b: u64,
    /// Subgroup order (prime, `< 2^62`).
    pub n: u64,
    pub name: String,
}

impl FastCurve {
    pub fn new(name: &str, p: u64, a: u64, b: u64, n: u64) -> Self {
        let f = Fp64::new(p);
        FastCurve {
            f,
            a: f.to_mont(a),
            b: f.to_mont(b),
            n,
            name: name.to_string(),
        }
    }

    /// Adopt a crate [`CurveParams`] whose `p` and `n` fit in 62 bits.
    pub fn from_params(c: &CurveParams) -> Option<Self> {
        let p = biguint_to_u64(&c.p)?;
        let n = biguint_to_u64(&c.n)?;
        if p >= (1u64 << MAX_BITS) || n >= (1u64 << MAX_BITS) {
            return None;
        }
        Some(Self::new(
            c.name,
            p,
            biguint_to_u64(&c.a)?,
            biguint_to_u64(&c.b)?,
            n,
        ))
    }

    /// Group-order field context (Montgomery mod `n`) for the linear
    /// algebra.
    pub fn scalar_field(&self) -> Fp64 {
        Fp64::new(self.n)
    }

    pub fn point(&self, x: u64, y: u64) -> Pt {
        Pt {
            x: self.f.to_mont(x),
            y: self.f.to_mont(y),
        }
    }

    pub fn canonical(&self, p: Pt) -> (u64, u64) {
        (self.f.from_mont(p.x), self.f.from_mont(p.y))
    }

    pub fn to_textbook(&self, p: Option<Pt>) -> Point {
        match p {
            None => Point::Infinity,
            Some(pt) => {
                let (x, y) = self.canonical(pt);
                let m = BigUint::from(self.f.p);
                Point::Affine {
                    x: FieldElement::new(BigUint::from(x), m.clone()),
                    y: FieldElement::new(BigUint::from(y), m),
                }
            }
        }
    }

    pub fn from_textbook(&self, p: &Point) -> Option<Pt> {
        match p {
            Point::Infinity => None,
            Point::Affine { x, y } => Some(self.point(
                biguint_to_u64(&x.value).expect("coordinate fits u64"),
                biguint_to_u64(&y.value).expect("coordinate fits u64"),
            )),
        }
    }

    pub fn is_on_curve(&self, p: Pt) -> bool {
        let f = &self.f;
        let lhs = f.sqr(p.y);
        let rhs = f.add(f.add(f.mul(f.sqr(p.x), p.x), f.mul(self.a, p.x)), self.b);
        lhs == rhs
    }

    /// `p1 + p2` for `x1 ≠ x2`, given `inv_dx = (x2 − x1)⁻¹`.
    #[inline(always)]
    pub fn add_with_inv(&self, p1: Pt, p2: Pt, inv_dx: u64) -> Pt {
        let f = &self.f;
        let lam = f.mul(f.sub(p2.y, p1.y), inv_dx);
        let x3 = f.sub(f.sub(f.sqr(lam), p1.x), p2.x);
        let y3 = f.sub(f.mul(lam, f.sub(p1.x, x3)), p1.y);
        Pt { x: x3, y: y3 }
    }

    #[inline]
    pub fn neg(&self, p: Pt) -> Pt {
        Pt {
            x: p.x,
            y: self.f.neg(p.y),
        }
    }

    pub fn double(&self, p: Pt) -> Option<Pt> {
        if p.y == 0 {
            return None;
        }
        let f = &self.f;
        let x2 = f.sqr(p.x);
        let num = f.add(f.add(f.add(x2, x2), x2), self.a);
        let den = f.inv(f.add(p.y, p.y));
        let lam = f.mul(num, den);
        let x3 = f.sub(f.sqr(lam), f.add(p.x, p.x));
        let y3 = f.sub(f.mul(lam, f.sub(p.x, x3)), p.y);
        Some(Pt { x: x3, y: y3 })
    }

    /// General addition with the point at infinity as `None`.
    pub fn add(&self, p1: Option<Pt>, p2: Option<Pt>) -> Option<Pt> {
        match (p1, p2) {
            (None, q) | (q, None) => q,
            (Some(a), Some(b)) => {
                if a.x == b.x {
                    if a.y == b.y {
                        self.double(a)
                    } else {
                        None
                    }
                } else {
                    let inv = self.f.inv(self.f.sub(b.x, a.x));
                    Some(self.add_with_inv(a, b, inv))
                }
            }
        }
    }

    /// `[k]P` by double-and-add (variable time; everything here is public).
    pub fn mul(&self, p: Option<Pt>, mut k: u64) -> Option<Pt> {
        let mut acc: Option<Pt> = None;
        let mut base = p;
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(acc, base);
            }
            base = self.add(base, base);
            k >>= 1;
        }
        acc
    }

    /// `[u]G + [v]Q`.
    pub fn combine(&self, g: Pt, q: Pt, u: u64, v: u64) -> Option<Pt> {
        self.add(self.mul(Some(g), u), self.mul(Some(q), v))
    }

    /// Lift a canonical x-coordinate to the curve point whose
    /// Montgomery-form `y ≤ (p−1)/2` (the sign convention every table
    /// in this module uses).  `None` if `x³ + ax + b` is a non-residue.
    pub fn lift_x(&self, x_canonical: u64) -> Option<Pt> {
        let f = &self.f;
        let x = f.to_mont(x_canonical);
        let rhs = f.add(f.add(f.mul(f.sqr(x), x), f.mul(self.a, x)), self.b);
        let y = f.sqrt(rhs)?;
        let y = if y > f.half { f.neg(y) } else { y };
        Some(Pt { x, y })
    }

    /// Negation-map representative: the one of `±P` with `y ≤ (p−1)/2`
    /// (in Montgomery form; any fixed choice works, this one is free).
    #[inline(always)]
    pub fn canonical_sign(&self, p: Pt) -> (Pt, bool) {
        if p.y > self.f.half {
            (self.neg(p), true)
        } else {
            (p, false)
        }
    }
}

fn biguint_to_u64(v: &BigUint) -> Option<u64> {
    let digits = v.to_u64_digits();
    match digits.len() {
        0 => Some(0),
        1 => Some(digits[0]),
        _ => None,
    }
}

// ── Parallel distinguished-point rho ──────────────────────────────────

#[derive(Clone, Debug)]
pub struct RhoConfig {
    pub threads: usize,
    /// Lock-step walkers per thread (the Montgomery-trick batch width).
    pub walkers: usize,
    /// `log₂` of the r-adding table size.
    pub table_bits: u32,
    /// A point is distinguished when this many hashed bits of `x` vanish.
    pub dp_bits: u32,
    pub negation: bool,
    pub seed: u64,
    /// Give up after this many total steps (0 = unlimited).
    pub max_steps: u64,
}

impl RhoConfig {
    /// Heuristic sizing for a prime-order group of `bits` bits: enough
    /// walkers to batch well without the DP start-up overhead
    /// (`walkers · 2^dp_bits`) exceeding a few percent of the expected
    /// `√(πn/4)` steps.
    pub fn for_bits(bits: u32, negation: bool) -> Self {
        let avail = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        let expected = (std::f64::consts::PI * 2f64.powi(bits as i32) / if negation { 4.0 } else { 2.0 })
            .sqrt();
        let walkers_total = (expected / 64.0).clamp(32.0, (avail * 1024) as f64) as usize;
        let threads = (walkers_total / 64).clamp(1, avail);
        let walkers = (walkers_total / threads).max(16);
        let dp_bits = (expected / (8.0 * (walkers * threads) as f64))
            .log2()
            .floor()
            .clamp(0.0, 24.0) as u32;
        RhoConfig {
            threads,
            walkers,
            table_bits: 10,
            dp_bits,
            negation,
            seed: 1,
            max_steps: 0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct RhoResult {
    pub log: u64,
    pub total_steps: u64,
    pub distinguished_points: u64,
    pub wall_secs: f64,
    pub steps_per_sec: f64,
    pub threads: usize,
    pub walkers_per_thread: usize,
    pub dp_bits: u32,
}

#[derive(Clone, Copy)]
struct TableEntry {
    pt: Pt,
    u: u64,
    v: u64,
}

struct Shared {
    dps: Mutex<HashMap<u64, (u64, u64, u64)>>, // x → (y, a, b)
    done: AtomicBool,
    answer: Mutex<Option<u64>>,
    total_steps: AtomicU64,
    dp_count: AtomicU64,
}

#[inline(always)]
fn branch_index(x: u64, table_bits: u32) -> usize {
    (x.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - table_bits)) as usize
}

#[inline(always)]
fn is_distinguished(x: u64, dp_bits: u32) -> bool {
    if dp_bits == 0 {
        return true;
    }
    (x.rotate_left(23).wrapping_mul(0xD6E8_FEB8_6659_FD93) >> (64 - dp_bits)) == 0
}

/// Solve `aG + bQ` collisions: `a1 + b1·k ≡ ±(a2 + b2·k)`.
fn solve_collision(
    curve: &FastCurve,
    g: Pt,
    q: Pt,
    (y1, a1, b1): (u64, u64, u64),
    (y2, a2, b2): (u64, u64, u64),
) -> Option<u64> {
    let fn_ = curve.scalar_field();
    let n = curve.n;
    let candidates: [(u64, u64); 1] = if y1 == y2 {
        // a1 + b1 k = a2 + b2 k  ⟹  k = (a1 − a2) / (b2 − b1)
        [(fn_.sub(a1, a2), fn_.sub(b2, b1))]
    } else {
        // a1 + b1 k = −(a2 + b2 k)  ⟹  k = −(a1 + a2) / (b1 + b2)
        [(neg_mod(add_mod(a1, a2, n), n), add_mod(b1, b2, n))]
    };
    for (num, den) in candidates {
        if den == 0 {
            continue;
        }
        let k = ((num as u128 * inv_mod_u64(den, n) as u128) % n as u128) as u64;
        if curve.mul(Some(g), k) == Some(q) {
            return Some(k);
        }
    }
    None
}

/// Parallel distinguished-point Pollard rho for `Q = [k]G` in a
/// prime-order subgroup of order `curve.n`.
pub fn rho_parallel(curve: &FastCurve, g: Pt, q: Pt, cfg: &RhoConfig) -> Option<RhoResult> {
    let started = Instant::now();
    let n = curve.n;
    // Shared r-adding table T_j = u_j G + v_j Q.
    let mut rng = StdRng::seed_from_u64(cfg.seed ^ 0xA5A5_5A5A_1234_5678);
    let table_len = 1usize << cfg.table_bits;
    let mut table = Vec::with_capacity(table_len);
    while table.len() < table_len {
        let u = rng.gen_range(1..n);
        let v = rng.gen_range(1..n);
        if let Some(pt) = curve.combine(g, q, u, v) {
            table.push(TableEntry { pt, u, v });
        }
    }
    let shared = Shared {
        dps: Mutex::new(HashMap::new()),
        done: AtomicBool::new(false),
        answer: Mutex::new(None),
        total_steps: AtomicU64::new(0),
        dp_count: AtomicU64::new(0),
    };
    std::thread::scope(|s| {
        for tid in 0..cfg.threads {
            let table = &table;
            let shared = &shared;
            s.spawn(move || rho_worker(curve, g, q, table, cfg, tid, shared));
        }
    });
    let log = (*shared.answer.lock().unwrap())?;
    let wall = started.elapsed().as_secs_f64();
    let total_steps = shared.total_steps.load(Ordering::Relaxed);
    Some(RhoResult {
        log,
        total_steps,
        distinguished_points: shared.dp_count.load(Ordering::Relaxed),
        wall_secs: wall,
        steps_per_sec: total_steps as f64 / wall.max(1e-9),
        threads: cfg.threads,
        walkers_per_thread: cfg.walkers,
        dp_bits: cfg.dp_bits,
    })
}

fn rho_worker(
    curve: &FastCurve,
    g: Pt,
    q: Pt,
    table: &[TableEntry],
    cfg: &RhoConfig,
    tid: usize,
    shared: &Shared,
) {
    let f = curve.f;
    let n = curve.n;
    let w = cfg.walkers;
    let mut rng = StdRng::seed_from_u64(
        cfg.seed ^ (tid as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15),
    );
    let step_cap: u32 = (40u64 << cfg.dp_bits).max(256).min(u32::MAX as u64) as u32;

    let mut xs = vec![0u64; w];
    let mut ys = vec![0u64; w];
    let mut as_ = vec![0u64; w];
    let mut bs = vec![0u64; w];
    let mut prev = vec![u64::MAX; w];
    let mut steps = vec![0u32; w];
    let mut d = vec![0u64; w];
    let mut js = vec![0usize; w];
    let mut reset = vec![false; w];
    let mut scratch = Vec::with_capacity(w);

    let restart = |i: usize,
                   rng: &mut StdRng,
                   xs: &mut [u64],
                   ys: &mut [u64],
                   as_: &mut [u64],
                   bs: &mut [u64],
                   prev: &mut [u64],
                   steps: &mut [u32]| {
        loop {
            let a = rng.gen_range(1..n);
            let b = rng.gen_range(1..n);
            if let Some(pt) = curve.combine(g, q, a, b) {
                let (pt, flipped) = if cfg.negation {
                    curve.canonical_sign(pt)
                } else {
                    (pt, false)
                };
                xs[i] = pt.x;
                ys[i] = pt.y;
                as_[i] = if flipped { neg_mod(a, n) } else { a };
                bs[i] = if flipped { neg_mod(b, n) } else { b };
                prev[i] = u64::MAX;
                steps[i] = 0;
                return;
            }
        }
    };
    for i in 0..w {
        restart(i, &mut rng, &mut xs, &mut ys, &mut as_, &mut bs, &mut prev, &mut steps);
    }

    let mut local_steps = 0u64;
    let mut local_dps = 0u64;
    'outer: loop {
        if shared.done.load(Ordering::Relaxed) {
            break;
        }
        if cfg.max_steps > 0 && shared.total_steps.load(Ordering::Relaxed) > cfg.max_steps {
            break;
        }
        for i in 0..w {
            let j = branch_index(xs[i], cfg.table_bits);
            js[i] = j;
            let dd = f.sub(table[j].pt.x, xs[i]);
            if dd == 0 {
                reset[i] = true;
                d[i] = f.one;
            } else {
                reset[i] = false;
                d[i] = dd;
            }
        }
        f.batch_inv(&mut d, &mut scratch);
        for i in 0..w {
            if reset[i] {
                restart(i, &mut rng, &mut xs, &mut ys, &mut as_, &mut bs, &mut prev, &mut steps);
                continue;
            }
            let t = &table[js[i]];
            let lam = f.mul(f.sub(t.pt.y, ys[i]), d[i]);
            let x3 = f.sub(f.sub(f.sqr(lam), xs[i]), t.pt.x);
            let mut y3 = f.sub(f.mul(lam, f.sub(xs[i], x3)), ys[i]);
            let mut a = add_mod(as_[i], t.u, n);
            let mut b = add_mod(bs[i], t.v, n);
            if cfg.negation && y3 > f.half {
                y3 = f.p - y3;
                a = neg_mod(a, n);
                b = neg_mod(b, n);
            }
            if x3 == prev[i] {
                // Fruitless 2-cycle {prev, cur}: both members are known
                // (prev == new).  Escape deterministically by doubling the
                // member with the smaller x, so merged walks stay merged.
                let (cx, cy, ca, cb) = if xs[i] < x3 {
                    (xs[i], ys[i], as_[i], bs[i])
                } else {
                    (x3, y3, a, b)
                };
                match curve.double(Pt { x: cx, y: cy }) {
                    Some(pt) => {
                        let (pt, flipped) = if cfg.negation {
                            curve.canonical_sign(pt)
                        } else {
                            (pt, false)
                        };
                        let a2 = add_mod(ca, ca, n);
                        let b2 = add_mod(cb, cb, n);
                        xs[i] = pt.x;
                        ys[i] = pt.y;
                        as_[i] = if flipped { neg_mod(a2, n) } else { a2 };
                        bs[i] = if flipped { neg_mod(b2, n) } else { b2 };
                        prev[i] = u64::MAX;
                    }
                    None => {
                        restart(i, &mut rng, &mut xs, &mut ys, &mut as_, &mut bs, &mut prev, &mut steps);
                        continue;
                    }
                }
            } else {
                prev[i] = xs[i];
                xs[i] = x3;
                ys[i] = y3;
                as_[i] = a;
                bs[i] = b;
            }
            steps[i] += 1;
            if is_distinguished(xs[i], cfg.dp_bits) {
                local_dps += 1;
                let mut dps = shared.dps.lock().unwrap();
                match dps.get(&xs[i]) {
                    Some(&other) => {
                        if (other.1, other.2) != (as_[i], bs[i]) {
                            if let Some(k) =
                                solve_collision(curve, g, q, (ys[i], as_[i], bs[i]), other)
                            {
                                *shared.answer.lock().unwrap() = Some(k);
                                shared.done.store(true, Ordering::Relaxed);
                                drop(dps);
                                break 'outer;
                            }
                        }
                    }
                    None => {
                        dps.insert(xs[i], (ys[i], as_[i], bs[i]));
                    }
                }
                drop(dps);
                steps[i] = 0;
            } else if steps[i] > step_cap {
                // Probable fruitless cycle of length ≥ 4: abandon the trail.
                restart(i, &mut rng, &mut xs, &mut ys, &mut as_, &mut bs, &mut prev, &mut steps);
            }
        }
        local_steps += w as u64;
        if local_steps >= 1 << 16 {
            shared.total_steps.fetch_add(local_steps, Ordering::Relaxed);
            shared.dp_count.fetch_add(local_dps, Ordering::Relaxed);
            local_steps = 0;
            local_dps = 0;
        }
    }
    shared.total_steps.fetch_add(local_steps, Ordering::Relaxed);
    shared.dp_count.fetch_add(local_dps, Ordering::Relaxed);
}

// ── Index calculus: factor base, pair table, relations, linear algebra ─

/// Small-x factor base: the points with the smallest x-coordinates
/// `x = 1, 2, 3, …` that lie on the curve, with the `y ≤ (p−1)/2` sign.
/// Same policy as [`crate::cryptanalysis::ec_index_calculus::build_factor_base`].
#[derive(Clone, Debug)]
pub struct FastFactorBase {
    pub pts: Vec<Pt>,
    /// Canonical x of each entry.
    pub xs: Vec<u64>,
    pub max_x: u64,
    /// `index[x] = i` for `x = xs[i]`, else `u32::MAX`.
    index: Vec<u32>,
}

impl FastFactorBase {
    pub fn build(curve: &FastCurve, size: usize) -> Self {
        let mut pts = Vec::with_capacity(size);
        let mut xs = Vec::with_capacity(size);
        let mut x = 1u64;
        while pts.len() < size && x < curve.f.p {
            if let Some(pt) = curve.lift_x(x) {
                pts.push(pt);
                xs.push(x);
            }
            x += 1;
        }
        let max_x = xs.last().copied().unwrap_or(0);
        let mut index = vec![u32::MAX; max_x as usize + 1];
        for (i, &x) in xs.iter().enumerate() {
            index[x as usize] = i as u32;
        }
        FastFactorBase {
            pts,
            xs,
            max_x,
            index,
        }
    }

    pub fn len(&self) -> usize {
        self.pts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pts.is_empty()
    }

    #[inline(always)]
    pub fn lookup(&self, x_canonical: u64) -> Option<usize> {
        if x_canonical > self.max_x {
            return None;
        }
        let i = self.index[x_canonical as usize];
        if i == u32::MAX {
            None
        } else {
            Some(i as usize)
        }
    }
}

/// Meet-in-the-middle table for 3-decompositions: every
/// `x(F_j + F_k)` (`j ≤ k`) and `x(F_j − F_k)` (`j < k`), sorted by
/// canonical x.  `B²` entries of 16 bytes.
#[derive(Clone, Debug)]
pub struct PairTable {
    /// `(x_canonical, j, k_with_sign)`; the top bit of the third field
    /// is set for a difference `F_j − F_k`.
    entries: Vec<(u64, u32, u32)>,
}

const PAIR_DIFF_BIT: u32 = 1 << 31;

impl PairTable {
    pub fn build(curve: &FastCurve, fb: &FastFactorBase) -> Self {
        let f = curve.f;
        let b = fb.len();
        let mut entries = Vec::with_capacity(b * b);
        let mut d = Vec::with_capacity(b);
        let mut scratch = Vec::with_capacity(b);
        for j in 0..b {
            let fj = fb.pts[j];
            if let Some(dbl) = curve.double(fj) {
                entries.push((f.from_mont(dbl.x), j as u32, j as u32));
            }
            // One inverse per (j, k) pair serves both F_j + F_k and F_j − F_k.
            d.clear();
            for k in (j + 1)..b {
                d.push(f.sub(fb.pts[k].x, fj.x));
            }
            f.batch_inv(&mut d, &mut scratch);
            for (idx, k) in ((j + 1)..b).enumerate() {
                let fk = fb.pts[k];
                let inv = d[idx];
                let sum = curve.add_with_inv(fj, fk, inv);
                let diff = curve.add_with_inv(fj, curve.neg(fk), inv);
                entries.push((f.from_mont(sum.x), j as u32, k as u32));
                entries.push((f.from_mont(diff.x), j as u32, k as u32 | PAIR_DIFF_BIT));
            }
        }
        entries.sort_unstable();
        PairTable { entries }
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// All `(j, k, sigma)` with `x(F_j + sigma·F_k) == x`.
    #[inline]
    pub fn lookup(&self, x: u64) -> &[(u64, u32, u32)] {
        let lo = self.entries.partition_point(|e| e.0 < x);
        let mut hi = lo;
        while hi < self.entries.len() && self.entries[hi].0 == x {
            hi += 1;
        }
        &self.entries[lo..hi]
    }
}

/// One relation `a·G + b·Q = Σ entries[i].1 · F_{entries[i].0}`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastRelation {
    pub a: u64,
    pub b: u64,
    pub entries: Vec<(u32, i8)>,
}

fn merge_entries(mut entries: Vec<(u32, i8)>) -> Vec<(u32, i8)> {
    entries.sort_unstable_by_key(|e| e.0);
    let mut merged: Vec<(u32, i8)> = Vec::with_capacity(entries.len());
    for (idx, s) in entries {
        if let Some(last) = merged.last_mut() {
            if last.0 == idx {
                last.1 += s;
                continue;
            }
        }
        merged.push((idx, s));
    }
    merged.retain(|e| e.1 != 0);
    merged
}

/// Fail-closed check: re-add the decomposition in the group.
fn relation_holds(curve: &FastCurve, fb: &FastFactorBase, r: Pt, entries: &[(u32, i8)]) -> bool {
    let mut acc: Option<Pt> = None;
    for &(idx, s) in entries {
        let pt = fb.pts[idx as usize];
        let pt = if s < 0 { curve.neg(pt) } else { pt };
        for _ in 0..s.unsigned_abs() {
            acc = curve.add(acc, Some(pt));
        }
    }
    acc == Some(r)
}

#[derive(Clone, Debug, Default)]
pub struct CollectStats {
    pub targets: u64,
    pub relations: usize,
    pub rows_one: usize,
    pub rows_two: usize,
    pub rows_three: usize,
    pub rejected: usize,
}

/// Collect `needed` relations by walking targets `R = aG + bQ`
/// (fresh random `(a, b)` after every productive target, `R ← R + G`
/// otherwise) and decomposing each by direct subtraction:
///
/// - 2-decomposition: `x(R ∓ F_i) ∈ FB` — one batched inverse per
///   target, one array lookup per `(i, sign)`;
/// - 3-decomposition (when `pairs` is given): `x(R ∓ F_i)` matched
///   against the pair table.
pub fn collect_relations(
    curve: &FastCurve,
    g: Pt,
    q: Pt,
    fb: &FastFactorBase,
    pairs: Option<&PairTable>,
    needed: usize,
    max_targets: u64,
    seed: u64,
) -> (Vec<FastRelation>, CollectStats) {
    let f = curve.f;
    let n = curve.n;
    let b_len = fb.len();
    let mut rng = StdRng::seed_from_u64(seed);
    let mut rels: Vec<FastRelation> = Vec::with_capacity(needed);
    let mut stats = CollectStats::default();
    let mut d = vec![0u64; b_len];
    let mut degenerate = vec![false; b_len];
    let mut scratch = Vec::with_capacity(b_len);
    let mut seen: HashSet<Vec<(u32, i8)>> = HashSet::new();

    let fresh = |rng: &mut StdRng| -> (u64, u64, Pt) {
        loop {
            let a = rng.gen_range(1..n);
            let b = rng.gen_range(1..n);
            if let Some(r) = curve.combine(g, q, a, b) {
                return (a, b, r);
            }
        }
    };
    let (mut a, mut b, mut r) = fresh(&mut rng);

    while rels.len() < needed && stats.targets < max_targets {
        stats.targets += 1;
        for i in 0..b_len {
            let dd = f.sub(r.x, fb.pts[i].x);
            degenerate[i] = dd == 0;
            d[i] = if dd == 0 { f.one } else { dd };
        }
        f.batch_inv(&mut d, &mut scratch);
        seen.clear();
        let mut found: Vec<Vec<(u32, i8)>> = Vec::new();
        for i in 0..b_len {
            let fi = fb.pts[i];
            if degenerate[i] {
                // R = ±F_i outright.
                let s = if r.y == fi.y { 1i8 } else { -1 };
                found.push(vec![(i as u32, s)]);
                continue;
            }
            let inv = d[i];
            // T = R − ε·F_i for ε = +1 (lam uses y_R + y_i) and ε = −1
            // (lam uses y_R − y_i); both share the inverse of x_R − x_i.
            let lam_minus = f.mul(f.add(r.y, fi.y), inv);
            let lam_plus = f.mul(f.sub(r.y, fi.y), inv);
            for (eps, lam) in [(1i8, lam_minus), (-1i8, lam_plus)] {
                let x3 = f.sub(f.sub(f.sqr(lam), r.x), fi.x);
                let xc = f.from_mont(x3);
                let mut y3: Option<u64> = None;
                if let Some(j) = fb.lookup(xc) {
                    if j >= i {
                        let y = f.sub(f.mul(lam, f.sub(r.x, x3)), r.y);
                        y3 = Some(y);
                        let tau = if y == fb.pts[j].y { 1i8 } else { -1 };
                        found.push(vec![(i as u32, eps), (j as u32, tau)]);
                    }
                }
                if let Some(pairs) = pairs {
                    if xc <= pairs.entries.last().map(|e| e.0).unwrap_or(0) {
                        for &(_, j, k_raw) in pairs.lookup(xc) {
                            if (j as usize) < i {
                                continue;
                            }
                            let sigma = if k_raw & PAIR_DIFF_BIT != 0 { -1i8 } else { 1 };
                            let k = (k_raw & !PAIR_DIFF_BIT) as usize;
                            let fk = fb.pts[k];
                            let fk = if sigma < 0 { curve.neg(fk) } else { fk };
                            let Some(pjk) = curve.add(Some(fb.pts[j as usize]), Some(fk)) else {
                                continue;
                            };
                            let y = *y3.get_or_insert_with(|| f.sub(f.mul(lam, f.sub(r.x, x3)), r.y));
                            let tau = if y == pjk.y { 1i8 } else { -1 };
                            found.push(vec![(i as u32, eps), (j, tau), (k as u32, tau * sigma)]);
                        }
                    }
                }
            }
        }
        let mut productive = false;
        for entries in found {
            let merged = merge_entries(entries);
            if merged.is_empty() || !seen.insert(merged.clone()) {
                continue;
            }
            if !relation_holds(curve, fb, r, &merged) {
                stats.rejected += 1;
                continue;
            }
            match merged.len() {
                1 => stats.rows_one += 1,
                2 => stats.rows_two += 1,
                _ => stats.rows_three += 1,
            }
            rels.push(FastRelation {
                a,
                b,
                entries: merged,
            });
            productive = true;
            if rels.len() >= needed {
                break;
            }
        }
        if productive {
            let (na, nb, nr) = fresh(&mut rng);
            a = na;
            b = nb;
            r = nr;
        } else {
            match curve.add(Some(r), Some(g)) {
                Some(nr) => {
                    r = nr;
                    a = add_mod(a, 1, n);
                }
                None => {
                    let (na, nb, nr) = fresh(&mut rng);
                    a = na;
                    b = nb;
                    r = nr;
                }
            }
        }
    }
    stats.relations = rels.len();
    (rels, stats)
}

/// Solve the relation system for `x = log_G Q` by Gauss–Jordan
/// elimination mod `n` on `u64` (Montgomery form).  Unknowns are
/// `(y_0 … y_{m−1}, x)`; row `r` reads
/// `Σ e_j y_j − b·x ≡ a (mod n)`.  Returns `None` unless `x` is pinned.
pub fn solve_relations(curve: &FastCurve, m: usize, rels: &[FastRelation]) -> Option<u64> {
    let fnum = curve.scalar_field();
    let n = curve.n;
    let cols = m + 1;
    let mut rows: Vec<Vec<u64>> = Vec::with_capacity(rels.len());
    let mut rhs: Vec<u64> = Vec::with_capacity(rels.len());
    for rel in rels {
        let mut row = vec![0u64; cols];
        for &(j, s) in &rel.entries {
            let v = if s < 0 {
                neg_mod(s.unsigned_abs() as u64 % n, n)
            } else {
                s as u64 % n
            };
            row[j as usize] = fnum.to_mont(v);
        }
        row[m] = fnum.to_mont(neg_mod(rel.b % n, n));
        rows.push(row);
        rhs.push(fnum.to_mont(rel.a % n));
    }
    let nrows = rows.len();
    let mut pivot_row_of_col = vec![usize::MAX; cols];
    let mut row = 0usize;
    for col in 0..cols {
        if row >= nrows {
            break;
        }
        let Some(piv) = (row..nrows).find(|&r| rows[r][col] != 0) else {
            continue;
        };
        rows.swap(row, piv);
        rhs.swap(row, piv);
        let inv = fnum.inv(rows[row][col]);
        for c in 0..cols {
            rows[row][c] = fnum.mul(rows[row][c], inv);
        }
        rhs[row] = fnum.mul(rhs[row], inv);
        let (pivot_vals, pivot_rhs) = (rows[row].clone(), rhs[row]);
        for r in 0..nrows {
            if r == row || rows[r][col] == 0 {
                continue;
            }
            let factor = rows[r][col];
            for c in col..cols {
                if pivot_vals[c] != 0 {
                    rows[r][c] = fnum.sub(rows[r][c], fnum.mul(factor, pivot_vals[c]));
                }
            }
            rhs[r] = fnum.sub(rhs[r], fnum.mul(factor, pivot_rhs));
        }
        pivot_row_of_col[col] = row;
        row += 1;
    }
    // Consistency: zero rows must have zero rhs.
    for r in row..nrows {
        if rhs[r] != 0 && rows[r].iter().all(|&v| v == 0) {
            return None;
        }
    }
    let pr = pivot_row_of_col[m];
    if pr == usize::MAX {
        return None;
    }
    Some(fnum.from_mont(rhs[pr]))
}

#[derive(Clone, Debug, Default)]
pub struct IcReport {
    pub factor_base_ms: f64,
    pub pair_table_ms: f64,
    pub relations_ms: f64,
    pub linear_algebra_ms: f64,
    pub verify_ms: f64,
    pub total_ms: f64,
    pub factor_base_size: usize,
    pub pair_table_entries: usize,
    pub relations: usize,
    pub targets: u64,
    pub rows_one: usize,
    pub rows_two: usize,
    pub rows_three: usize,
}

/// End-to-end index calculus: factor base → (pair table) → relations →
/// linear algebra → group verification.  `summands` is 2 or 3.
pub fn ic_solve(
    curve: &FastCurve,
    g: Pt,
    q: Pt,
    fb_size: usize,
    extra_relations: usize,
    summands: u32,
    max_targets: u64,
    seed: u64,
) -> Option<(u64, IcReport)> {
    let t0 = Instant::now();
    let fb = FastFactorBase::build(curve, fb_size);
    let factor_base_ms = t0.elapsed().as_secs_f64() * 1e3;
    if fb.is_empty() {
        return None;
    }
    let t1 = Instant::now();
    let pairs = if summands >= 3 {
        Some(PairTable::build(curve, &fb))
    } else {
        None
    };
    let pair_table_ms = t1.elapsed().as_secs_f64() * 1e3;
    let needed = fb.len() + extra_relations;
    let t2 = Instant::now();
    let (rels, stats) =
        collect_relations(curve, g, q, &fb, pairs.as_ref(), needed, max_targets, seed);
    let relations_ms = t2.elapsed().as_secs_f64() * 1e3;
    if rels.len() < needed {
        return None;
    }
    let t3 = Instant::now();
    let x = solve_relations(curve, fb.len(), &rels)?;
    let linear_algebra_ms = t3.elapsed().as_secs_f64() * 1e3;
    let t4 = Instant::now();
    if curve.mul(Some(g), x) != Some(q) {
        return None;
    }
    let verify_ms = t4.elapsed().as_secs_f64() * 1e3;
    Some((
        x,
        IcReport {
            factor_base_ms,
            pair_table_ms,
            relations_ms,
            linear_algebra_ms,
            verify_ms,
            total_ms: t0.elapsed().as_secs_f64() * 1e3,
            factor_base_size: fb.len(),
            pair_table_entries: pairs.map(|p| p.len()).unwrap_or(0),
            relations: rels.len(),
            targets: stats.targets,
            rows_one: stats.rows_one,
            rows_two: stats.rows_two,
            rows_three: stats.rows_three,
        },
    ))
}

// ── Point counting and bench-curve generation ─────────────────────────

/// Order of `g` found as the unique `N` in the Hasse interval with
/// `[N]g = O`, by baby-step giant-step on `t = p + 1 − N`.  Returns
/// `None` when the interval holds several candidates (small-order `g`)
/// or the match is not unique — callers only accept prime results.
pub fn point_order_hasse(curve: &FastCurve, g: Pt) -> Option<u64> {
    let p = curve.f.p;
    let s = (p as f64).sqrt().ceil() as u64 + 1;
    let lo: i128 = -(2 * s as i128);
    let hi: i128 = 2 * s as i128;
    let width = (hi - lo + 1) as u64;
    let m = ((width as f64).sqrt().ceil() as u64).max(2);
    // Baby steps: x([j]G) → j, with y for sign resolution.
    let mut baby: HashMap<u64, (u64, u64)> = HashMap::with_capacity(m as usize);
    let mut cur: Option<Pt> = None;
    for j in 0..m {
        if let Some(pt) = cur {
            baby.entry(pt.x).or_insert((j, pt.y));
        }
        cur = curve.add(cur, Some(g));
    }
    let q0 = curve.mul(Some(g), p + 1); // [p+1]G = [t]G
    let m_g = curve.mul(Some(g), m);
    let neg_m_g = m_g.map(|pt| curve.neg(pt));
    // Y_i = Q0 − [lo + i·m]G, i = 0 ..= width/m.
    let lo_g = {
        let mag = lo.unsigned_abs() as u64;
        curve.mul(Some(g), mag).map(|pt| curve.neg(pt)) // [lo]G with lo < 0
    };
    let mut y = curve.add(q0, lo_g.map(|pt| curve.neg(pt))); // Q0 − [lo]G
    let mut candidates: Vec<i128> = Vec::new();
    let iters = width / m + 1;
    for i in 0..=iters {
        match y {
            None => candidates.push(lo + (i as i128) * m as i128),
            Some(pt) => {
                if let Some(&(j, by)) = baby.get(&pt.x) {
                    let t = if by == pt.y {
                        lo + (i as i128) * m as i128 + j as i128
                    } else {
                        lo + (i as i128) * m as i128 - j as i128
                    };
                    candidates.push(t);
                }
            }
        }
        y = curve.add(y, neg_m_g);
    }
    candidates.sort_unstable();
    candidates.dedup();
    candidates.retain(|&t| t >= lo && t <= hi && curve.mul(Some(g), (p as i128 + 1 - t) as u64).is_none());
    if candidates.len() != 1 {
        return None;
    }
    Some((p as i128 + 1 - candidates[0]) as u64)
}

/// Deterministic `a = −3` prime-order bench curve at `bits` bits in the
/// requested residue class of `p mod 4`: `p` the largest such prime
/// below `2^bits`, `b` the smallest positive value giving prime order,
/// `G` the smallest-x point.  Reproduces
/// [`crate::cryptanalysis::research_bench::bench_curves_a_minus_3`].
pub fn find_a3_curve(bits: u32, residue_mod_4: u64) -> Option<(FastCurve, Pt)> {
    assert!((8..=MAX_BITS).contains(&bits));
    let mut p = (1u64 << bits) - 1;
    while !(p % 4 == residue_mod_4 && is_prime_u64(p)) {
        p -= 1;
    }
    let a = p - 3;
    for b in 1u64..1_000_000 {
        // Singular if 4a³ + 27b² ≡ 0.
        let f = Fp64::new(p);
        let am = f.to_mont(a);
        let bm = f.to_mont(b);
        let disc = f.add(
            f.mul(f.to_mont(4), f.mul(f.sqr(am), am)),
            f.mul(f.to_mont(27), f.sqr(bm)),
        );
        if disc == 0 {
            continue;
        }
        let probe = FastCurve::new("probe", p, a, b, 1);
        let Some(g) = (1u64..10_000).find_map(|x| probe.lift_x(x)) else {
            continue;
        };
        let Some(order) = point_order_hasse(&probe, g) else {
            continue;
        };
        if !is_prime_u64(order) {
            continue;
        }
        let name = format!(
            "a3-fast-{bits}bit-{}",
            if residue_mod_4 == 3 { "p256class" } else { "cryptoproclass" }
        );
        let curve = FastCurve::new(&name, p, a, b, order);
        let g = curve.lift_x(probe.canonical(g).0)?;
        if curve.mul(Some(g), order).is_some() {
            continue;
        }
        return Some((curve, g));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::research_bench::bench_curves_a_minus_3;

    fn ladder(bits: u32) -> (FastCurve, Pt) {
        let (_, c) = bench_curves_a_minus_3()
            .into_iter()
            .find(|(b, c)| *b == bits && c.name.contains("p256class"))
            .expect("ladder rung");
        let fc = FastCurve::from_params(&c).unwrap();
        let g = fc.from_textbook(&c.generator()).unwrap();
        (fc, g)
    }

    #[test]
    fn montgomery_matches_biguint() {
        let mut rng = StdRng::seed_from_u64(7);
        for &p in &[65_519u64, 1_048_571, 16_777_199, (1u64 << 61) - 1, 4_611_686_018_427_387_847] {
            let f = Fp64::new(p);
            let pb = BigUint::from(p);
            for _ in 0..200 {
                let a = rng.gen_range(0..p);
                let b = rng.gen_range(0..p);
                let (am, bm) = (f.to_mont(a), f.to_mont(b));
                assert_eq!(f.from_mont(am), a);
                let prod = ((BigUint::from(a) * BigUint::from(b)) % &pb).to_u64_digits();
                assert_eq!(f.from_mont(f.mul(am, bm)), prod.first().copied().unwrap_or(0));
                assert_eq!(f.from_mont(f.add(am, bm)), ((a as u128 + b as u128) % p as u128) as u64);
                assert_eq!(f.from_mont(f.sub(am, bm)), ((a as u128 + p as u128 - b as u128) % p as u128) as u64);
                if a != 0 {
                    let inv = f.inv(am);
                    assert_eq!(f.mul(inv, am), f.one);
                    let e = BigUint::from(a).modpow(&(&pb - BigUint::from(2u32)), &pb);
                    assert_eq!(f.from_mont(inv), e.to_u64_digits()[0]);
                }
                if let Some(r) = f.sqrt(am) {
                    assert_eq!(f.sqr(r), am);
                }
            }
            // Batch inversion agrees with pointwise.
            let xs: Vec<u64> = (0..37).map(|_| f.to_mont(rng.gen_range(1..p))).collect();
            let mut batched = xs.clone();
            f.batch_inv(&mut batched, &mut Vec::new());
            for (x, b) in xs.iter().zip(batched.iter()) {
                assert_eq!(f.inv(*x), *b);
            }
        }
    }

    #[test]
    fn fast_arithmetic_matches_textbook_points() {
        for bits in [16u32, 20, 24] {
            let (fc, g) = ladder(bits);
            let (_, c) = bench_curves_a_minus_3()
                .into_iter()
                .find(|(b, c)| *b == bits && c.name.contains("p256class"))
                .unwrap();
            let a_fe = c.a_fe();
            let tg = c.generator();
            assert!(fc.is_on_curve(g));
            let mut rng = StdRng::seed_from_u64(bits as u64);
            for _ in 0..20 {
                let k = rng.gen_range(1..fc.n);
                let l = rng.gen_range(1..fc.n);
                let kg = fc.mul(Some(g), k);
                let lg = fc.mul(Some(g), l);
                let sum = fc.add(kg, lg);
                let textbook = tg
                    .scalar_mul(&BigUint::from(k), &a_fe)
                    .add(&tg.scalar_mul(&BigUint::from(l), &a_fe), &a_fe);
                assert_eq!(fc.to_textbook(sum), textbook);
            }
            assert_eq!(fc.mul(Some(g), fc.n), None);
        }
    }

    #[test]
    fn rho_recovers_known_logs() {
        for bits in [16u32, 20, 24] {
            let (fc, g) = ladder(bits);
            let mut rng = StdRng::seed_from_u64(99 + bits as u64);
            for negation in [true, false] {
                let k = rng.gen_range(1..fc.n);
                let q = fc.mul(Some(g), k).unwrap();
                let mut cfg = RhoConfig::for_bits(bits, negation);
                cfg.threads = 2;
                cfg.seed = 5;
                let res = rho_parallel(&fc, g, q, &cfg).expect("rho solves");
                assert_eq!(res.log, k, "bits={bits} negation={negation}");
            }
        }
    }

    #[test]
    fn ic_two_and_three_decomposition_recover_known_logs() {
        for bits in [16u32, 20] {
            let (fc, g) = ladder(bits);
            let mut rng = StdRng::seed_from_u64(3 + bits as u64);
            let k = rng.gen_range(1..fc.n);
            let q = fc.mul(Some(g), k).unwrap();
            for summands in [2u32, 3] {
                let (x, rep) = ic_solve(&fc, g, q, 60, 6, summands, 50_000_000, 11)
                    .unwrap_or_else(|| panic!("ic bits={bits} summands={summands}"));
                assert_eq!(x, k);
                assert_eq!(rep.relations, rep.factor_base_size + 6);
                if summands == 3 {
                    assert!(rep.rows_three > 0);
                }
            }
        }
    }

    #[test]
    fn fast_relations_match_the_s3_sweep() {
        // Every 2-decomposition the direct subtraction finds satisfies
        // S₃(x_R, x_i, x_j) = 0, and the staged BigUint driver's relation
        // for the same target is among them.
        use crate::cryptanalysis::ec_index_calculus::{build_factor_base, semaev_s3};
        let (fc, g) = ladder(16);
        let (_, c) = bench_curves_a_minus_3()
            .into_iter()
            .find(|(b, c)| *b == 16 && c.name.contains("p256class"))
            .unwrap();
        let fb = FastFactorBase::build(&fc, 40);
        let fb_text = build_factor_base(&c, 40);
        assert_eq!(fb.len(), fb_text.len());
        let q = fc.mul(Some(g), 12345).unwrap();
        let (rels, _) = collect_relations(&fc, g, q, &fb, None, 5, 1_000_000, 1);
        assert_eq!(rels.len(), 5);
        let a_fe = c.a_fe();
        let b_fe = c.fe(c.b.clone());
        for rel in &rels {
            if rel.entries.len() != 2 {
                continue;
            }
            let r = fc.combine(g, q, rel.a, rel.b).unwrap();
            let (xr, _) = fc.canonical(r);
            let xi = fb.xs[rel.entries[0].0 as usize];
            let xj = fb.xs[rel.entries[1].0 as usize];
            let s3 = semaev_s3(&c.fe(BigUint::from(xr)), &c.fe(BigUint::from(xi)), &c.fe(BigUint::from(xj)), &a_fe, &b_fe);
            assert!(s3.is_zero(), "S3 must vanish on a found relation");
        }
    }

    #[test]
    fn hasse_bsgs_reproduces_the_ladder_rungs() {
        // The committed 16/20/24-bit a=−3 rungs, and the 28-bit P-256-class
        // rung found independently by `examples/find_a3_bench_curves`
        // (p = 268435399, b = 3, G = (1, 1), n = 268407199).
        for bits in [16u32, 20, 24] {
            let (fc, g) = ladder(bits);
            assert_eq!(point_order_hasse(&fc, g), Some(fc.n), "bits={bits}");
            let (gen, gg) = find_a3_curve(bits, 3).unwrap();
            assert_eq!((gen.f.p, gen.n, gen.canonical(gg)), (fc.f.p, fc.n, fc.canonical(g)));
        }
        let (gen, gg) = find_a3_curve(28, 3).unwrap();
        assert_eq!(gen.f.p, 268_435_399);
        assert_eq!(gen.f.from_mont(gen.b), 3);
        assert_eq!(gen.canonical(gg), (1, 1));
        assert_eq!(gen.n, 268_407_199);
    }

    #[test]
    fn generated_curves_have_prime_order_and_solve() {
        let (fc, g) = find_a3_curve(36, 3).unwrap();
        assert!(is_prime_u64(fc.n));
        assert_eq!(fc.mul(Some(g), fc.n), None);
        let k = 0xDEAD_BEEFu64 % fc.n;
        let q = fc.mul(Some(g), k).unwrap();
        let mut cfg = RhoConfig::for_bits(36, true);
        cfg.threads = cfg.threads.min(4);
        let res = rho_parallel(&fc, g, q, &cfg).unwrap();
        assert_eq!(res.log, k);
    }
}
