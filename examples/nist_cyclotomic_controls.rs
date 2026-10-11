//! Controlled re-test of the comparative claims in
//! `research/notes/cm-isogeny/RESEARCH_NIST_CYCLOTOMIC_POINT_HEIGHTS.md`,
//! following the frozen protocol
//! `research/nist_cyclotomic_controls_20261011/PROTOCOL.md`
//! and `.claude/skills/controlled-experiment/SKILL.md`.
//!
//! E1  small-factor content of n−1, n+1, twist vs random Hasse-interval primes
//!     (+ planted-smooth positive control)
//! E2  small-point census power: planted integer / rational points (positive),
//!     random b (negative), NIST b (target)
//! E3  ĥ(G̃)/ln p vs lifted random points (+ small-x positive control)
//! E4  two-point-lift pairing vs random pairs (+ dependent-pair positive control)
//!
//! ```bash
//! cargo run --release --example nist_cyclotomic_controls -- \
//!     [--out research/nist_cyclotomic_controls_20261011] [--protocol v1|v2|quick]
//! ```
//! `v1` runs `PROTOCOL.md` as frozen.  `v2` runs `AMENDMENT_1.md`: the same
//! instruments, metrics and rules with the larger N the power check needs,
//! and E2 is not re-run because it has no rank test.  `quick` is a smoke run
//! at small N, and its verdicts mean nothing.  Outputs carry the mode as a
//! suffix: `rows_<mode>.jsonl`, `RESULTS_<mode>.md`, `results_<mode>.json`.

use std::fmt::Write as _;
use std::fs;
use std::io::Write as _;
use std::path::PathBuf;
use std::time::Instant;

use crypto_lib::cryptanalysis::curve_traits::arith::{self, is_prime, log2_abs};
use crypto_lib::cryptanalysis::ecm::ecm_factor;
use num_bigint::{BigInt, BigUint, RandBigInt};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde_json::json;

const LN2: f64 = std::f64::consts::LN_2;
const SEED: u64 = 20261011;
const VN: [&str; 3] = ["n-1", "n+1", "twist"];

struct NistCurve {
    name: &'static str,
    p: &'static str,
    n: &'static str,
    b_hex: &'static str,
    gx_hex: &'static str,
    gy_hex: &'static str,
}

const CURVES: &[NistCurve] = &[
    NistCurve {
        name: "P-192",
        p: "6277101735386680763835789423207666416083908700390324961279",
        n: "6277101735386680763835789423176059013767194773182842284081",
        b_hex: "64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1",
        gx_hex: "188da80eb03090f67cbf20eb43a18800f4ff0afd82ff1012",
        gy_hex: "07192b95ffc8da78631011ed6b24cdd573f977a11e794811",
    },
    NistCurve {
        name: "P-224",
        p: "26959946667150639794667015087019630673557916260026308143510066298881",
        n: "26959946667150639794667015087019625940457807714424391721682722368061",
        b_hex: "b4050a850c04b3abf54132565044b0b7d7bfd8ba270b39432355ffb4",
        gx_hex: "b70e0cbd6bb4bf7f321390b94a03c1d356c21122343280d6115c1d21",
        gy_hex: "bd376388b5f723fb4c22dfe6cd4375a05a07476444d5819985007e34",
    },
    NistCurve {
        name: "P-256",
        p: "115792089210356248762697446949407573530086143415290314195533631308867097853951",
        n: "115792089210356248762697446949407573529996955224135760342422259061068512044369",
        b_hex: "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
        gx_hex: "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296",
        gy_hex: "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5",
    },
];

fn big(s: &str) -> BigInt {
    BigInt::parse_bytes(s.as_bytes(), 10).expect("decimal")
}
fn bighex(s: &str) -> BigInt {
    BigInt::parse_bytes(s.as_bytes(), 16).expect("hex")
}
fn ubig(x: &BigInt) -> BigUint {
    x.to_biguint().expect("non-negative")
}
fn ln_abs(x: &BigInt) -> f64 {
    log2_abs(x) * LN2
}

// ────────────────────────────── rationals and E(Q) ──────────────────────────────

#[derive(Clone, Debug, PartialEq, Eq)]
struct Rat {
    n: BigInt,
    d: BigInt,
}
impl Rat {
    fn new(n: BigInt, d: BigInt) -> Rat {
        assert!(!d.is_zero());
        let g = n.gcd(&d);
        let (mut n, mut d) = (n / &g, d / &g);
        if d.is_negative() {
            n = -n;
            d = -d;
        }
        Rat { n, d }
    }
    fn int(x: &BigInt) -> Rat {
        Rat {
            n: x.clone(),
            d: BigInt::one(),
        }
    }
    fn small(x: i64) -> Rat {
        Rat::int(&BigInt::from(x))
    }
    fn add(&self, o: &Rat) -> Rat {
        Rat::new(&self.n * &o.d + &o.n * &self.d, &self.d * &o.d)
    }
    fn sub(&self, o: &Rat) -> Rat {
        Rat::new(&self.n * &o.d - &o.n * &self.d, &self.d * &o.d)
    }
    fn mul(&self, o: &Rat) -> Rat {
        Rat::new(&self.n * &o.n, &self.d * &o.d)
    }
    fn div(&self, o: &Rat) -> Rat {
        Rat::new(&self.n * &o.d, &self.d * &o.n)
    }
    fn is_zero(&self) -> bool {
        self.n.is_zero()
    }
}

/// `y² = x³ + A x + B` over `Q`.
struct CurveQ {
    a: Rat,
    b: Rat,
}
type Pt = Option<(Rat, Rat)>;

impl CurveQ {
    fn on_curve(&self, p: &Pt) -> bool {
        match p {
            None => true,
            Some((x, y)) => {
                let lhs = y.mul(y);
                let rhs = x.mul(x).mul(x).add(&self.a.mul(x)).add(&self.b);
                lhs == rhs
            }
        }
    }
    fn add(&self, p: &Pt, q: &Pt) -> Pt {
        let (x1, y1) = match p {
            None => return q.clone(),
            Some(v) => v,
        };
        let (x2, y2) = match q {
            None => return p.clone(),
            Some(v) => v,
        };
        let lambda = if x1 == x2 {
            if y1.add(y2).is_zero() {
                return None;
            }
            let num = Rat::small(3).mul(&x1.mul(x1)).add(&self.a);
            num.div(&Rat::small(2).mul(y1))
        } else {
            y2.sub(y1).div(&x2.sub(x1))
        };
        let x3 = lambda.mul(&lambda).sub(x1).sub(x2);
        let y3 = lambda.mul(&x1.sub(&x3)).sub(y1);
        Some((x3, y3))
    }
    /// Discriminant `−16(4A³ + 27B²)`.
    /// `ĥ(P) ≈ h_x([2^k]P) / 4^k` for `k = 0..=doublings`; returns the
    /// sequence so the reader sees the convergence.  The curve and point
    /// are first moved to an integral model `(x, y) ↦ (u²x, u³y)` (the
    /// canonical height is an isomorphism invariant and the naive-height
    /// limit does not depend on the model), then doubled in Jacobian
    /// coordinates over `Z`, with one gcd per evaluated height.
    fn height_sequence(&self, p: &Pt, doublings: u32) -> Vec<f64> {
        self.height_sequence_on(p, doublings, false)
    }

    /// As [`Self::height_sequence`], but with `rational = true` the naive
    /// height is taken of the reduced rational `x = X/(Z²u²)` on the
    /// original model, not of `X/Z²` on the denominator-cleared one.  The
    /// cleared model inflates `h_x` by about `2 log u` (u includes the
    /// point's own denominators), which survives the `4^{-k}` as a bias.
    /// Amendment 2 (E4r) uses `rational = true`; v1/v2 keep `false`.
    fn height_sequence_on(&self, p: &Pt, doublings: u32, rational: bool) -> Vec<f64> {
        let (x, y) = match p {
            None => return vec![0.0; doublings as usize + 1],
            Some(v) => v,
        };
        let u: BigInt = self.a.d.lcm(&self.b.d).lcm(&x.d).lcm(&y.d);
        let u2 = &u * &u;
        let a_int = self.a.mul(&Rat::int(&(&u2 * &u2)));
        assert!(a_int.d.is_one());
        let a = a_int.n;
        let xs = x.mul(&Rat::int(&u2));
        let ys = y.mul(&Rat::int(&(&u2 * &u)));
        assert!(xs.d.is_one() && ys.d.is_one());
        let (mut xj, mut yj, mut zj) = (xs.n, ys.n, BigInt::one());
        let mut seq = Vec::new();
        for k in 0..=doublings {
            let z2 = &zj * &zj;
            let den = if rational { &z2 * &u2 } else { z2.clone() };
            let g = xj.gcd(&den);
            let h = ln_abs(&(&xj / &g)).max(ln_abs(&(&den / &g)));
            seq.push(h / 4f64.powi(k as i32));
            if k < doublings {
                // Jacobian doubling: S = 4XY², M = 3X² + aZ⁴, X' = M² − 2S, Y' = M(S − X') − 8Y⁴, Z' = 2YZ.
                let y2 = &yj * &yj;
                let s = BigInt::from(4) * &xj * &y2;
                let m = BigInt::from(3) * &xj * &xj + &a * &z2 * &z2;
                let x3 = &m * &m - BigInt::from(2) * &s;
                let y3 = &m * (&s - &x3) - BigInt::from(8) * &y2 * &y2;
                let z3 = BigInt::from(2) * &yj * &zj;
                xj = x3;
                yj = y3;
                zj = z3;
            }
        }
        seq
    }
}

/// Affine arithmetic on `y² = x³ − 3x + b` over `F_p` for the lifts.
fn ec_add_fp(
    p: &BigInt,
    a: &BigInt,
    p1: &(BigInt, BigInt),
    p2: &(BigInt, BigInt),
) -> (BigInt, BigInt) {
    let modp = |x: BigInt| x.mod_floor(p);
    let inv = |x: &BigInt| modp(x.clone()).modpow(&(p - 2), p);
    let lambda = if p1.0 == p2.0 {
        modp((BigInt::from(3) * &p1.0 * &p1.0 + a) * inv(&(BigInt::from(2) * &p1.1)))
    } else {
        modp((&p2.1 - &p1.1) * inv(&(&p2.0 - &p1.0)))
    };
    let x3 = modp(&lambda * &lambda - &p1.0 - &p2.0);
    let y3 = modp(&lambda * (&p1.0 - &x3) - &p1.1);
    (x3, y3)
}
fn ec_mul_fp(p: &BigInt, a: &BigInt, k: &BigInt, g: &(BigInt, BigInt)) -> (BigInt, BigInt) {
    let mut r: Option<(BigInt, BigInt)> = None;
    for i in (0..k.bits()).rev() {
        if let Some(cur) = &r {
            r = Some(ec_add_fp(p, a, cur, cur));
        }
        if k.bit(i) {
            r = Some(match &r {
                None => g.clone(),
                Some(cur) => ec_add_fp(p, a, cur, g),
            });
        }
    }
    r.expect("k > 0")
}

// ────────────────────────────── statistics ──────────────────────────────

/// Two-sided empirical p-value of `t` against `ctrl` (add-one smoothing),
/// plus the target's rank (number of controls strictly below it).
fn rank_p(t: f64, ctrl: &[f64]) -> (usize, f64) {
    let n = ctrl.len() as f64;
    let ge = ctrl.iter().filter(|&&c| c >= t).count() as f64;
    let le = ctrl.iter().filter(|&&c| c <= t).count() as f64;
    let below = ctrl.iter().filter(|&&c| c < t).count();
    let p = (2.0 * ((1.0 + ge) / (n + 1.0)).min((1.0 + le) / (n + 1.0))).min(1.0);
    (below, p)
}
fn quantile(v: &[f64], q: f64) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let i = ((s.len() - 1) as f64 * q).round() as usize;
    s[i]
}

// ────────────────────────────── E1 ──────────────────────────────

/// M1 = log2 of the product (with multiplicity) of returned prime factors < 2^32;
/// M2 = number of distinct returned prime factors < 2^32. Instrument: arith::factor(x, 2^18).
fn small_factor_metrics(x: &BigInt) -> (f64, f64) {
    let f = arith::factor(&ubig(x), 1 << 18);
    let lim = BigUint::one() << 32;
    let mut m1 = 0.0;
    let mut m2 = 0.0;
    for (q, e) in &f.primes {
        if *q < lim {
            m1 += *e as f64 * q.to_f64().unwrap().log2();
            m2 += 1.0;
        }
    }
    (m1, m2)
}

fn random_prime_in(rng: &mut StdRng, lo: &BigInt, hi: &BigInt) -> BigInt {
    loop {
        let c = rng.gen_bigint_range(lo, hi);
        let c = if c.is_even() { c + 1 } else { c };
        if &c < hi && is_prime(&ubig(&c)) {
            return c;
        }
    }
}

/// Planted positive control: random 2^20-smooth part of ~120 bits times a random
/// prime bringing the total to `bits` bits.
fn planted_smooth(rng: &mut StdRng, bits: u64) -> BigInt {
    let mut s = BigInt::one();
    while s.bits() < 120 {
        let q = loop {
            let c: u64 = rng.gen_range(3..(1u64 << 20));
            if is_prime(&BigUint::from(c)) {
                break c;
            }
        };
        s *= q;
    }
    let rest = bits - s.bits();
    let lo = BigInt::one() << (rest - 1);
    let hi = BigInt::one() << rest;
    s * random_prime_in(rng, &lo, &hi)
}

// ────────────────────────────── E2 census ──────────────────────────────

/// The census of `nist_cyclotomic_point_heights.rs`, verbatim: integer |x| ≤ sx,
/// rational x = u/d², 2 ≤ d ≤ sd, |u| ≤ su, residue filter mod 64, 63, 65, 11.
/// Returns (integer points, rational points).
fn census(b: &BigInt, sx: i64, sd: i64, su: i64) -> (Vec<i64>, Vec<(i64, i64)>) {
    let sq_table = |m: u64| -> Vec<bool> {
        let mut v = vec![false; m as usize];
        for i in 0..m {
            v[((i * i) % m) as usize] = true;
        }
        v
    };
    let mods = [64u64, 63, 65, 11];
    let tables: Vec<Vec<bool>> = mods.iter().map(|&m| sq_table(m)).collect();
    let bmods: Vec<i128> = mods
        .iter()
        .map(|&m| b.mod_floor(&BigInt::from(m)).to_i128().unwrap())
        .collect();
    let passes = |u: i128, dpow4: i128, dpow6: i128| -> bool {
        mods.iter()
            .zip(&tables)
            .zip(&bmods)
            .all(|((&m, tab), &bm)| {
                let m = m as i128;
                let v = ((u % m) * (u % m) % m * (u % m) % m - 3 * (u % m) % m * (dpow4 % m) % m
                    + bm * (dpow6 % m) % m)
                    .rem_euclid(m);
                tab[v as usize]
            })
    };
    let mut ints = Vec::new();
    for x in -sx..=sx {
        if !passes(x as i128, 1, 1) {
            continue;
        }
        let xb = BigInt::from(x);
        let v = &xb * &xb * &xb - BigInt::from(3) * &xb + b;
        if v.is_negative() {
            continue;
        }
        let s = v.sqrt();
        if &s * &s == v {
            ints.push(x);
        }
    }
    let mut rats = Vec::new();
    for d in 2..=sd {
        let (d4, d6) = ((d as i128).pow(4), (d as i128).pow(6));
        let (bd4, bd6) = (BigInt::from(d4), BigInt::from(d6));
        for u in -su..=su {
            if u.gcd(&d) != 1 || !passes(u as i128, d4, d6) {
                continue;
            }
            let ub = BigInt::from(u);
            let w = &ub * &ub * &ub - BigInt::from(3) * &ub * &bd4 + b * &bd6;
            if w.is_negative() {
                continue;
            }
            let s = w.sqrt();
            if &s * &s == w {
                rats.push((u, d));
            }
        }
    }
    (ints, rats)
}

// ────────────────────────────── E3/E4 helpers ──────────────────────────────

fn sqrt_mod(a: &BigInt, p: &BigInt) -> Option<BigInt> {
    // p ≡ 3 (mod 4) for P-192/P-256; Tonelli–Shanks for P-224.
    let a = a.mod_floor(p);
    if a.is_zero() {
        return Some(BigInt::zero());
    }
    let one = BigInt::one();
    let pm1 = p - &one;
    if a.modpow(&(&pm1 >> 1), p) != one {
        return None;
    }
    if (p % 4u32) == BigInt::from(3) {
        return Some(a.modpow(&((p + 1u32) >> 2), p));
    }
    let mut q = pm1.clone();
    let mut s = 0u32;
    while q.is_even() {
        q >>= 1;
        s += 1;
    }
    let mut z = BigInt::from(2);
    while z.modpow(&(&pm1 >> 1), p) == one {
        z += 1;
    }
    let (mut m, mut c, mut t, mut r) = (
        s,
        z.modpow(&q, p),
        a.modpow(&q, p),
        a.modpow(&((&q + 1u32) >> 1), p),
    );
    while t != one {
        let mut i = 0u32;
        let mut t2 = t.clone();
        while t2 != one {
            t2 = (&t2 * &t2) % p;
            i += 1;
        }
        let b = c.modpow(&(BigInt::one() << (m - i - 1)), p);
        m = i;
        c = (&b * &b) % p;
        t = (&t * &c) % p;
        r = (&r * &b) % p;
    }
    Some(r)
}

fn random_point(
    rng: &mut StdRng,
    p: &BigInt,
    b: &BigInt,
    xmax: Option<&BigInt>,
) -> (BigInt, BigInt) {
    loop {
        let x = match xmax {
            Some(m) => rng.gen_bigint_range(&BigInt::one(), m),
            None => rng.gen_bigint_range(&BigInt::zero(), p),
        };
        let rhs = (&x * &x * &x - BigInt::from(3) * &x + b).mod_floor(p);
        if let Some(y) = sqrt_mod(&rhs, p) {
            return (x, y);
        }
    }
}

/// One-point lift through (x, y) ∈ Z², a = −3: returns (ĥ, h_x) in nats.
fn one_point_lift_height(x: &BigInt, y: &BigInt) -> (f64, f64) {
    let bq = y * y - x * x * x + BigInt::from(3) * x;
    let e = CurveQ {
        a: Rat::small(-3),
        b: Rat::int(&bq),
    };
    let pt: Pt = Some((Rat::int(x), Rat::int(y)));
    let seq = e.height_sequence(&pt, 5);
    (*seq.last().unwrap(), seq[0])
}

/// Two-point lift through (x1,y1), (x2,y2) ∈ Z²: returns (|cos|, det/(h1 h2)).
fn two_point_pairing(p1: &(BigInt, BigInt), p2: &(BigInt, BigInt)) -> (f64, f64) {
    let (x1, y1, x2, y2) = (
        Rat::int(&p1.0),
        Rat::int(&p1.1),
        Rat::int(&p2.0),
        Rat::int(&p2.1),
    );
    let c1 = y1.mul(&y1).sub(&x1.mul(&x1).mul(&x1));
    let c2 = y2.mul(&y2).sub(&x2.mul(&x2).mul(&x2));
    let a = c1.sub(&c2).div(&x1.sub(&x2));
    let b = c1.sub(&a.mul(&x1));
    let e = CurveQ { a, b };
    let pg: Pt = Some((x1, y1));
    let pr: Pt = Some((x2, y2));
    assert!(e.on_curve(&pg) && e.on_curve(&pr));
    pairing_on(&e, &pg, &pr)
}

fn pairing_on(e: &CurveQ, pg: &Pt, pr: &Pt) -> (f64, f64) {
    // The three heights are independent; run them concurrently (same values).
    let h = |q: &Pt| *e.height_sequence(q, 5).last().unwrap();
    let ((hg, hr), hgr) = rayon::join(|| rayon::join(|| h(pg), || h(pr)), || h(&e.add(pg, pr)));
    let pair = (hgr - hg - hr) / 2.0;
    (
        (pair / (hg * hr).sqrt()).abs(),
        (hg * hr - pair * pair) / (hg * hr),
    )
}

/// Amendment 2 (E4r): the pairing from rational-model heights, at k = 5
/// and k = 4 doublings: `(|cos|, det/(ĥĥ), |cos|_k4, det_k4)`.
fn pairing_on_rational(e: &CurveQ, pg: &Pt, pr: &Pt) -> [f64; 4] {
    let h = |q: &Pt| e.height_sequence_on(q, 5, true);
    let ((sg, sr), sgr) = rayon::join(|| rayon::join(|| h(pg), || h(pr)), || h(&e.add(pg, pr)));
    let at = |k: usize| {
        let (hg, hr, hgr) = (sg[k], sr[k], sgr[k]);
        let pair = (hgr - hg - hr) / 2.0;
        (
            (pair / (hg * hr).sqrt()).abs(),
            (hg * hr - pair * pair) / (hg * hr),
        )
    };
    let ((c5, r5), (c4, r4)) = (at(5), at(4));
    [c5, r5, c4, r4]
}

fn two_point_curve(p1: &(BigInt, BigInt), p2: &(BigInt, BigInt)) -> (CurveQ, Pt, Pt) {
    let (x1, y1, x2, y2) = (
        Rat::int(&p1.0),
        Rat::int(&p1.1),
        Rat::int(&p2.0),
        Rat::int(&p2.1),
    );
    let c1 = y1.mul(&y1).sub(&x1.mul(&x1).mul(&x1));
    let c2 = y2.mul(&y2).sub(&x2.mul(&x2).mul(&x2));
    let a = c1.sub(&c2).div(&x1.sub(&x2));
    let b = c1.sub(&a.mul(&x1));
    let e = CurveQ { a, b };
    let (pg, pr): (Pt, Pt) = (Some((x1, y1)), Some((x2, y2)));
    assert!(e.on_curve(&pg) && e.on_curve(&pr));
    (e, pg, pr)
}

/// Amendment 2 (E1b): `factor(x, 2^18)` as in E1, then ECM (B1 = 2000,
/// 20 curves) on every composite it leaves, recursively.  Returns
/// `(M1, M2, M1_rho, M2_rho, unsplit_composites)`; the `_rho` pair is the
/// v1/v2 instrument's value on the same input.
fn small_factor_metrics_ecm(x: &BigInt, seed: u64) -> (f64, f64, f64, f64, usize) {
    let f = arith::factor(&ubig(x), 1 << 18);
    let lim = BigUint::one() << 32;
    let mut primes: Vec<(BigUint, u32)> = f.primes.clone();
    let (mut m1r, mut m2r) = (0.0, 0.0);
    for (q, e) in &f.primes {
        if *q < lim {
            m1r += *e as f64 * q.to_f64().unwrap().log2();
            m2r += 1.0;
        }
    }
    let mut stack: Vec<(BigUint, u32)> = f.composites.clone();
    let mut unsplit = 0usize;
    let mut s = seed;
    while let Some((c, e)) = stack.pop() {
        if is_prime(&c) {
            primes.push((c, e));
            continue;
        }
        s = s.wrapping_add(1000);
        match ecm_factor(&c, 2000, 20, s) {
            Some(d) => {
                let other = &c / &d;
                stack.push((d, e));
                stack.push((other, e));
            }
            None => unsplit += 1,
        }
    }
    // Merge repeated primes (ECM may split c = q·q·r into pieces sharing q).
    primes.sort();
    let mut merged: Vec<(BigUint, u32)> = Vec::new();
    for (q, e) in primes {
        match merged.last_mut() {
            Some((lq, le)) if *lq == q => *le += e,
            _ => merged.push((q, e)),
        }
    }
    let (mut m1, mut m2) = (0.0, 0.0);
    for (q, e) in &merged {
        if *q < lim {
            m1 += *e as f64 * q.to_f64().unwrap().log2();
            m2 += 1.0;
        }
    }
    (m1, m2, m1r, m2r, unsplit)
}

/// Amendment 2: E1b (ECM-completed small-factor content) and E4r
/// (rational-model heights for the two-point pairing).  `quick` shrinks N.
fn run_amendment2(out: &std::path::Path, mode: &str) {
    let quick = mode == "quick3";
    let (n_e1b, n_pos1b, n_e4) = if quick { (40, 5, 10) } else { (800, 50, 500) };
    let mut rows = fs::File::create(out.join(format!("rows_{mode}.jsonl"))).unwrap();
    let mut md = String::new();
    let mut summary = Vec::new();
    let started = Instant::now();
    writeln!(md, "# Controlled re-test: results ({mode})\n\n{} Seed base {SEED}. Protocol file hashes are recorded in `results_{mode}.json`. Every verdict below is the decision rule applied mechanically.\n", if quick { "QUICK3 smoke run of Amendment 2 at small N. Verdicts are meaningless." } else { "v3, `AMENDMENT_2.md`: E1b and E4r only." }).unwrap();

    // ── E1b ──
    let alpha = 0.05 / 9.0;
    let pmin = 2.0 / (n_e1b as f64 + 1.0);
    writeln!(md, "## E1b — M1 with ECM-completed factoring (rho 2^18, then ECM B1 = 2000 × 20 curves per composite)\n\nNull: the first N = {n_e1b} of v2's Hasse-interval primes per curve (same seeds). Family: 9 targets × M1, Bonferroni α = 0.05/9 = {alpha:.4}. Power check: min attainable p = {pmin:.4}, {}.\n", if pmin < alpha { "passes" } else { "**FAILS**" }).unwrap();
    writeln!(md, "| curve | value | M1 target (rho → ECM) | control median [2.5%, 97.5%] | rank / N | p (two-sided) | verdict | controls whose M1 changed under ECM | unsplit composites (controls) |\n|:--|:--|--:|:--|--:|--:|:--|--:|--:|").unwrap();
    let mut pos_lines = Vec::new();
    for (ci, c) in CURVES.iter().enumerate() {
        let p = big(c.p);
        let n = big(c.n);
        let mut rng = StdRng::seed_from_u64(SEED + ci as u64);
        let w = (BigInt::from(4) * &p).sqrt();
        let lo = &p + 1 - &w;
        let hi = &p + 1 + &w;
        let nps: Vec<BigInt> = (0..n_e1b)
            .map(|_| random_prime_in(&mut rng, &lo, &hi))
            .collect();
        let vals: Vec<BigInt> = nps
            .iter()
            .flat_map(|np| [np - 1, np + 1, BigInt::from(2) * &p + 2 - np])
            .collect();
        let base = 0xE1B0_0000u64 + (ci as u64) * 0x10_0000;
        let mets: Vec<_> = vals
            .par_iter()
            .enumerate()
            .map(|(i, v)| small_factor_metrics_ecm(v, base + i as u64))
            .collect();
        let mut ctrl: [Vec<f64>; 3] = Default::default();
        let mut changed = [0usize; 3];
        let mut unsplit = [0usize; 3];
        for (i, &(m1, m2, m1r, m2r, us)) in mets.iter().enumerate() {
            let k = i % 3;
            ctrl[k].push(m1);
            if (m1 - m1r).abs() > 1e-9 {
                changed[k] += 1;
            }
            unsplit[k] += us;
            writeln!(rows, "{}", json!({"exp":"E1b","arm":"control","curve":c.name,"i":i/3,"value":VN[k],"M1":m1,"M2":m2,"M1_rho":m1r,"M2_rho":m2r,"unsplit":us})).unwrap();
        }
        let targets = [&n - 1, &n + 1, BigInt::from(2) * &p + 2 - &n];
        let tm: Vec<_> = targets
            .par_iter()
            .enumerate()
            .map(|(k, v)| small_factor_metrics_ecm(v, base + 0xF_0000 + k as u64))
            .collect();
        for (k, &(m1, m2, m1r, m2r, us)) in tm.iter().enumerate() {
            writeln!(rows, "{}", json!({"exp":"E1b","arm":"target","curve":c.name,"value":VN[k],"M1":m1,"M2":m2,"M1_rho":m1r,"M2_rho":m2r,"unsplit":us})).unwrap();
            let cv = &ctrl[k];
            let (r, pv) = rank_p(m1, cv);
            let verdict = if pv < alpha {
                "**unusual**"
            } else {
                "within control distribution"
            };
            writeln!(md, "| {} | {} | {:.1} → {:.1} | {:.1} [{:.1}, {:.1}] | {}/{} | {:.4} | {} | {}/{} | {} |", c.name, ["n − 1", "n + 1", "twist"][k], m1r, m1, quantile(cv, 0.5), quantile(cv, 0.025), quantile(cv, 0.975), r, cv.len(), pv, verdict, changed[k], cv.len(), unsplit[k]).unwrap();
            summary.push(json!({"exp":"E1b","curve":c.name,"value":VN[k],"metric":"M1","target":m1,"target_rho":m1r,"median":quantile(cv,0.5),"rank":r,"N":cv.len(),"p":pv,"alpha":alpha,"unusual":pv<alpha,"controls_changed":changed[k],"unsplit":unsplit[k]}));
        }
        // Positive control: a prime in [2^31, 2^32) hidden in a value of n's size.
        let mut prng = StdRng::seed_from_u64(SEED + 400 + ci as u64);
        let planted: Vec<(BigInt, BigInt)> = (0..n_pos1b)
            .map(|_| {
                let q32 =
                    random_prime_in(&mut prng, &(BigInt::one() << 31), &(BigInt::one() << 32));
                let q100 =
                    random_prime_in(&mut prng, &(BigInt::one() << 99), &(BigInt::one() << 100));
                let rest = n.bits() - 132;
                let r = random_prime_in(
                    &mut prng,
                    &(BigInt::one() << (rest - 1)),
                    &(BigInt::one() << rest),
                );
                (&q32 * &q100 * &r, q32)
            })
            .collect();
        let pm: Vec<_> = planted
            .par_iter()
            .enumerate()
            .map(|(i, (v, _))| small_factor_metrics_ecm(v, base + 0xE_0000 + i as u64))
            .collect();
        let mut hits = 0;
        for ((_, q32), &(m1, _, m1r, _, _)) in planted.iter().zip(&pm) {
            let want = q32.to_f64().unwrap().log2();
            let hit = (m1 - want).abs() < 1e-6;
            hits += hit as usize;
            writeln!(rows, "{}", json!({"exp":"E1b","arm":"positive","curve":c.name,"M1":m1,"M1_rho":m1r,"planted_log2":want,"found":hit})).unwrap();
        }
        pos_lines.push(format!(
            "| {} | {}/{} planted 32-bit primes recovered exactly |",
            c.name, hits, n_pos1b
        ));
        summary.push(
            json!({"exp":"E1b","arm":"positive","curve":c.name,"detected":hits,"of":n_pos1b}),
        );
        eprintln!(
            "E1b {} done ({:.0}s)",
            c.name,
            started.elapsed().as_secs_f64()
        );
    }
    writeln!(md, "\nPositive control (`q32 · q100 · r`, q32 a random prime in [2^31, 2^32)); pass = M1 equals log2(q32) exactly:\n\n| curve | detection |\n|:--|:--|\n{}\n", pos_lines.join("\n")).unwrap();

    // ── E4r ──
    let alpha4 = 0.05 / 6.0;
    let pmin4 = 2.0 / (n_e4 as f64 + 1.0);
    writeln!(md, "## E4r — two-point-lift pairing from rational-model heights (k = 5; k = 4 for the accuracy check)\n\nNull: the same N = {n_e4} random pairs per curve as v2 (same seeds). Bonferroni α = 0.05/6 = {alpha4:.4}. Power check: min attainable p = {pmin4:.4}, {}. Accuracy gate (pre-registered): on the dependent pairs, abs(|cos| − 1) ≤ 1e-4 and abs(det/(ĥĥ)) ≤ 1e-4 at k = 5 on every curve, else no orthogonality verdict.\n", if pmin4 < alpha4 { "passes" } else { "**FAILS**" }).unwrap();
    writeln!(md, "| curve | metric | target | random pairs median [2.5%, 97.5%] (N) | rank | p | verdict | dependent pairs ([2]G̃; [3]G̃) | median abs(k5 − k4) over controls |\n|:--|:--|--:|:--|--:|--:|:--|:--|--:|").unwrap();
    let mut gate = true;
    for (ci, c) in CURVES.iter().enumerate() {
        let p = big(c.p);
        let b = bighex(c.b_hex);
        let mut rng = StdRng::seed_from_u64(SEED + 300 + ci as u64);
        let pairs: Vec<_> = (0..n_e4)
            .map(|_| {
                (
                    random_point(&mut rng, &p, &b, None),
                    random_point(&mut rng, &p, &b, None),
                )
            })
            .collect();
        let prs: Vec<[f64; 4]> = pairs
            .par_iter()
            .map(|(p1, p2)| {
                let (e, a, r) = two_point_curve(p1, p2);
                pairing_on_rational(&e, &a, &r)
            })
            .collect();
        eprintln!(
            "E4r {} random pairs done ({:.0}s)",
            c.name,
            started.elapsed().as_secs_f64()
        );
        for v in &prs {
            writeln!(rows, "{}", json!({"exp":"E4r","arm":"control","curve":c.name,"abs_cos":v[0],"norm_reg":v[1],"abs_cos_k4":v[2],"norm_reg_k4":v[3]})).unwrap();
        }
        let g = (bighex(c.gx_hex), bighex(c.gy_hex));
        let r = ec_mul_fp(&p, &BigInt::from(-3), &big("1311768467294899695"), &g);
        let bq = &g.1 * &g.1 - &g.0 * &g.0 * &g.0 + BigInt::from(3) * &g.0;
        let e1 = CurveQ {
            a: Rat::small(-3),
            b: Rat::int(&bq),
        };
        let gt: Pt = Some((Rat::int(&g.0), Rat::int(&g.1)));
        let g2 = e1.add(&gt, &gt);
        let g3 = e1.add(&g2, &gt);
        assert!(e1.on_curve(&g2) && e1.on_curve(&g3));
        let (t, (d2, d3)) = rayon::join(
            || {
                let (e, a, rr) = two_point_curve(&g, &r);
                pairing_on_rational(&e, &a, &rr)
            },
            || {
                rayon::join(
                    || pairing_on_rational(&e1, &gt, &g2),
                    || pairing_on_rational(&e1, &gt, &g3),
                )
            },
        );
        writeln!(rows, "{}", json!({"exp":"E4r","arm":"target","curve":c.name,"abs_cos":t[0],"norm_reg":t[1],"abs_cos_k4":t[2],"norm_reg_k4":t[3]})).unwrap();
        writeln!(rows, "{}", json!({"exp":"E4r","arm":"positive","curve":c.name,"g2":[d2[0],d2[1],d2[2],d2[3]],"g3":[d3[0],d3[1],d3[2],d3[3]]})).unwrap();
        let ok = [d2, d3]
            .iter()
            .all(|d| (d[0] - 1.0).abs() <= 1e-4 && d[1].abs() <= 1e-4);
        gate &= ok;
        for (mi, name) in [(0usize, "|cos|"), (1, "det/(ĥĥ)")] {
            let cv: Vec<f64> = prs.iter().map(|v| v[mi]).collect();
            let dk: Vec<f64> = prs.iter().map(|v| (v[mi] - v[mi + 2]).abs()).collect();
            let (rk, pv) = rank_p(t[mi], &cv);
            let verdict = if pv < alpha4 {
                "**unusual**"
            } else {
                "within control distribution"
            };
            writeln!(md, "| {} | {} | {:.3e} | {:.3e} [{:.3e}, {:.3e}] ({}) | {}/{} | {:.3} | {} | {:.3e}; {:.3e} | {:.2e} |", c.name, name, t[mi], quantile(&cv, 0.5), quantile(&cv, 0.025), quantile(&cv, 0.975), cv.len(), rk, cv.len(), pv, verdict, d2[mi], d3[mi], quantile(&dk, 0.5)).unwrap();
            summary.push(json!({"exp":"E4r","curve":c.name,"metric":name,"target":t[mi],"median":quantile(&cv,0.5),"q025":quantile(&cv,0.025),"q975":quantile(&cv,0.975),"rank":rk,"N":cv.len(),"p":pv,"alpha":alpha4,"unusual":pv<alpha4,"dependent":[d2[mi],d3[mi]],"median_abs_k5_minus_k4":quantile(&dk,0.5)}));
        }
        eprintln!(
            "E4r {} done ({:.0}s)",
            c.name,
            started.elapsed().as_secs_f64()
        );
    }
    writeln!(
        md,
        "\nAccuracy gate: **{}**.\n",
        if gate {
            "passed"
        } else {
            "FAILED — no orthogonality verdict is licensed"
        }
    )
    .unwrap();
    summary.push(json!({"exp":"E4r","accuracy_gate_passed":gate}));
    writeln!(
        md,
        "\n---\nRun time {:.0} s. Raw per-sample rows: `rows_{mode}.jsonl`.",
        started.elapsed().as_secs_f64()
    )
    .unwrap();
    let hash_of = |name: &str| -> String {
        let bytes = fs::read(format!("research/nist_cyclotomic_controls_20261011/{name}"))
            .unwrap_or_default();
        crypto_lib::hash::sha256::sha256(&bytes)
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect()
    };
    let hashes = json!({"PROTOCOL.md": hash_of("PROTOCOL.md"), "AMENDMENT_1.md": hash_of("AMENDMENT_1.md"), "AMENDMENT_2.md": hash_of("AMENDMENT_2.md")});
    fs::write(out.join(format!("RESULTS_{mode}.md")), &md).unwrap();
    fs::write(
        out.join(format!("results_{mode}.json")),
        serde_json::to_string_pretty(&json!({"protocol_sha256": hashes, "seed": SEED, "mode": mode, "N": {"E1b": n_e1b, "E1b_pos": n_pos1b, "E4r": n_e4}, "summary": summary})).unwrap(),
    )
    .unwrap();
    eprintln!("wrote {}", out.display());
}

// ────────────────────────────── figure ──────────────────────────────

/// One figure row: control values, the target, and positive-control values.
struct PlotRow {
    label: String,
    ctrl: Vec<f64>,
    target: f64,
    pos: Vec<f64>,
}

fn plot_rows(
    rows: &[serde_json::Value],
    exp: &str,
    key: &str,
    pos_key: Option<&str>,
) -> Vec<PlotRow> {
    let f = |r: &serde_json::Value, k: &str| r[k].as_f64();
    let mut out = Vec::new();
    for c in CURVES {
        let values: Vec<Option<&str>> = if exp.starts_with("E1") {
            VN.iter().map(|v| Some(*v)).collect()
        } else {
            vec![None]
        };
        for v in values {
            let sel = |arm: &'static str| {
                rows.iter().filter(move |r| {
                    r["exp"] == exp
                        && r["curve"] == c.name
                        && r["arm"] == arm
                        && (arm.starts_with("positive") || v.is_none_or(|v| r["value"] == v))
                })
            };
            let ctrl: Vec<f64> = sel("control").filter_map(|r| f(r, key)).collect();
            if ctrl.is_empty() {
                continue;
            }
            let target = sel("target").find_map(|r| f(r, key)).unwrap_or(f64::NAN);
            let pos: Vec<f64> = match (exp, pos_key) {
                ("E1" | "E1b", Some(pk)) if v == Some("n-1") => {
                    sel("positive").filter_map(|r| f(r, pk)).collect()
                }
                ("E3", Some(pk)) => sel("positive_small_x").filter_map(|r| f(r, pk)).collect(),
                ("E4" | "E4r", Some(pk)) => sel("positive")
                    .flat_map(|r| {
                        let i = if pk == "abs_cos" { 0 } else { 1 };
                        [r["g2"][i].as_f64(), r["g3"][i].as_f64()]
                    })
                    .flatten()
                    .collect(),
                _ => Vec::new(),
            };
            let label = match v {
                Some(v) => format!("{} {}", c.name, v),
                None => c.name.to_string(),
            };
            out.push(PlotRow {
                label,
                ctrl,
                target,
                pos,
            });
        }
    }
    out
}

/// Six panels (E1 M1, E1 M2, E3 two metrics, E4 two metrics): per row a
/// histogram of the control arm, its 2.5–97.5% band, the target (red
/// diamond) and positive controls (green ticks).
fn controls_svg(rows: &[serde_json::Value], mode: &str) -> String {
    let v3 = mode == "v3" || mode == "quick3";
    let panels: Vec<(String, Vec<PlotRow>)> = if v3 {
        vec![
            (
                "E1b · M1, ECM-completed (bits)".into(),
                plot_rows(rows, "E1b", "M1", Some("M1")),
            ),
            (
                "E1b · M1 before ECM (rho only, same inputs)".into(),
                plot_rows(rows, "E1b", "M1_rho", None),
            ),
            (
                "E4r · |cos ∠|, rational-model heights".into(),
                plot_rows(rows, "E4r", "abs_cos", Some("abs_cos")),
            ),
            (
                "E4r · det/(ĥĥ), rational-model heights".into(),
                plot_rows(rows, "E4r", "norm_reg", Some("norm_reg")),
            ),
        ]
    } else {
        vec![
            (
                "E1 · M1 = log2 of the 2^32-smooth part found (bits)".into(),
                plot_rows(rows, "E1", "M1", Some("M1")),
            ),
            (
                "E1 · M2 = # distinct primes below 2^32".into(),
                plot_rows(rows, "E1", "M2", None),
            ),
            (
                "E3 · ĥ/ln p of one-point lift".into(),
                plot_rows(rows, "E3", "hhat_over_lnp", Some("hhat_over_lnp")),
            ),
            (
                "E3 · h_x − ĥ (nats)".into(),
                plot_rows(rows, "E3", "hx_minus_hhat", Some("hx_minus_hhat")),
            ),
            (
                "E4 · |cos ∠| of two-point lift".into(),
                plot_rows(rows, "E4", "abs_cos", Some("abs_cos")),
            ),
            (
                "E4 · det/(ĥĥ) of two-point lift".into(),
                plot_rows(rows, "E4", "norm_reg", Some("norm_reg")),
            ),
        ]
    };
    let (pw, rh, top, left, gap) = (430.0, 34.0, 46.0, 104.0, 40.0);
    let col_h: Vec<f64> = panels
        .iter()
        .map(|(_, r)| r.len() as f64 * rh + 70.0)
        .collect();
    let h_row: Vec<f64> = col_h
        .chunks(2)
        .map(|c| c.iter().fold(0.0f64, |a, &b| a.max(b)))
        .collect();
    let w = 2.0 * (left + pw) + gap + 70.0;
    let h = top + h_row.iter().sum::<f64>() + 40.0;
    let mut s = String::new();
    writeln!(s, r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" font-family="Helvetica, Arial, sans-serif" font-size="11">"##).unwrap();
    writeln!(s, r##"<rect width="100%" height="100%" fill="white"/><text x="10" y="20" font-size="14" font-weight="bold">Controlled re-test ({mode}): control distributions (grey), target (red ◆), positive controls (green |)</text>"##).unwrap();
    writeln!(s, r##"<text x="10" y="36" fill="#555">Recorded values from rows_{mode}.jsonl. Band = control 2.5–97.5%. N per row is printed at the right. Seed base {SEED}.</text>"##).unwrap();
    for (pi, (title, prow)) in panels.iter().enumerate() {
        let x0 = left + (pi % 2) as f64 * (left + pw + gap);
        let y0 = top + h_row[..pi / 2].iter().sum::<f64>() + 16.0;
        writeln!(
            s,
            r##"<text x="{}" y="{}" font-weight="bold">{}</text>"##,
            x0 - left + 10.0,
            y0,
            title
        )
        .unwrap();
        // Axis spans controls and target; positive controls outside it are
        // marked at the edge with their range.
        let all: Vec<f64> = prow
            .iter()
            .flat_map(|r| r.ctrl.iter().chain(std::iter::once(&r.target)))
            .copied()
            .filter(|v| v.is_finite())
            .collect();
        if all.is_empty() {
            continue;
        }
        let (mut lo, mut hi) = all
            .iter()
            .fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
        if hi - lo < 1e-12 {
            hi = lo + 1.0;
        }
        let pad = 0.03 * (hi - lo);
        lo -= pad;
        hi += pad;
        let sx = |v: f64| x0 + (v - lo) / (hi - lo) * pw;
        for (ri, r) in prow.iter().enumerate() {
            let yc = y0 + 14.0 + ri as f64 * rh + rh / 2.0;
            writeln!(
                s,
                r##"<text x="{}" y="{}" text-anchor="end">{}</text>"##,
                x0 - 6.0,
                yc + 4.0,
                r.label
            )
            .unwrap();
            let (q1, q2) = (quantile(&r.ctrl, 0.025), quantile(&r.ctrl, 0.975));
            writeln!(
                s,
                r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="#e3e3e3"/>"##,
                sx(q1),
                yc - rh / 2.0 + 3.0,
                (sx(q2) - sx(q1)).max(1.0),
                rh - 6.0
            )
            .unwrap();
            let bins = 60usize;
            let mut cnt = vec![0usize; bins];
            for &v in &r.ctrl {
                let b = (((v - lo) / (hi - lo)) * bins as f64)
                    .floor()
                    .clamp(0.0, bins as f64 - 1.0) as usize;
                cnt[b] += 1;
            }
            let cmax = *cnt.iter().max().unwrap() as f64;
            for (b, &c) in cnt.iter().enumerate() {
                if c == 0 {
                    continue;
                }
                let bh = (rh - 8.0) * c as f64 / cmax;
                writeln!(
                    s,
                    r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="#888"/>"##,
                    x0 + b as f64 * pw / bins as f64,
                    yc + rh / 2.0 - 4.0 - bh,
                    pw / bins as f64 - 0.5,
                    bh
                )
                .unwrap();
            }
            let (mut below, mut above) = (Vec::new(), Vec::new());
            for &v in &r.pos {
                if v < lo {
                    below.push(v);
                } else if v > hi {
                    above.push(v);
                } else {
                    writeln!(s, r##"<line x1="{0:.1}" x2="{0:.1}" y1="{1:.1}" y2="{2:.1}" stroke="#1a9850" stroke-width="1.5"/>"##, sx(v), yc - rh / 2.0 + 2.0, yc + rh / 2.0 - 2.0).unwrap();
                }
            }
            let fmt = |v: f64| {
                if v.abs() >= 10.0 {
                    format!("{v:.0}")
                } else {
                    format!("{v:.3}")
                }
            };
            let span = |v: &[f64]| {
                v.iter()
                    .fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)))
            };
            if !below.is_empty() {
                let (a, b) = span(&below);
                writeln!(s, r##"<path d="M{:.1},{:.1} l8,-5 l0,10 z" fill="#1a9850"/><text x="{:.1}" y="{:.1}" fill="#1a9850" font-size="9">positive: {} to {} (n={})</text>"##, x0, yc, x0 + 10.0, yc - 6.0, fmt(a), fmt(b), below.len()).unwrap();
            }
            if !above.is_empty() {
                let (a, b) = span(&above);
                writeln!(s, r##"<path d="M{:.1},{:.1} l-8,-5 l0,10 z" fill="#1a9850"/><text x="{:.1}" y="{:.1}" fill="#1a9850" font-size="9" text-anchor="end">positive: {} to {} (n={})</text>"##, x0 + pw, yc, x0 + pw - 10.0, yc - 6.0, fmt(a), fmt(b), above.len()).unwrap();
            }
            if r.target.is_finite() {
                let tx = sx(r.target);
                writeln!(
                    s,
                    r##"<path d="M{:.1},{:.1} l5,6 l-5,6 l-5,-6 z" fill="#d73027"/>"##,
                    tx,
                    yc - 6.0
                )
                .unwrap();
            }
            writeln!(
                s,
                r##"<text x="{}" y="{}" fill="#555" font-size="9">N={}</text>"##,
                x0 + pw + 3.0,
                yc + 3.0,
                r.ctrl.len()
            )
            .unwrap();
        }
        let yb = y0 + 14.0 + prow.len() as f64 * rh;
        writeln!(
            s,
            r##"<line x1="{x0}" x2="{}" y1="{yb}" y2="{yb}" stroke="#333"/>"##,
            x0 + pw
        )
        .unwrap();
        for k in 0..=4 {
            let v = lo + (hi - lo) * k as f64 / 4.0;
            let lab = if (hi - lo).abs() < 0.001 {
                format!("{v:.6}")
            } else if (hi - lo).abs() < 0.1 {
                format!("{v:.4}")
            } else if (hi - lo).abs() < 10.0 {
                format!("{v:.2}")
            } else {
                format!("{v:.0}")
            };
            writeln!(s, r##"<line x1="{0:.1}" x2="{0:.1}" y1="{yb}" y2="{1}" stroke="#333"/><text x="{0:.1}" y="{2}" text-anchor="middle" font-size="9">{3}</text>"##, sx(v), yb + 4.0, yb + 15.0, lab).unwrap();
        }
    }
    s.push_str("</svg>\n");
    s
}

// ────────────────────────────── main ──────────────────────────────

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mode = args
        .iter()
        .position(|a| a == "--protocol")
        .map(|i| args[i + 1].clone())
        .unwrap_or_else(|| "v1".into());
    let out = args
        .iter()
        .position(|a| a == "--out")
        .map(|i| PathBuf::from(&args[i + 1]))
        .unwrap_or_else(|| PathBuf::from("research/nist_cyclotomic_controls_20261011"));
    fs::create_dir_all(&out).unwrap();
    if args.iter().any(|a| a == "--plot") {
        // Lay out recorded rows only; nothing is recomputed.
        let text = fs::read_to_string(out.join(format!("rows_{mode}.jsonl"))).expect("rows file");
        let rows: Vec<serde_json::Value> = text
            .lines()
            .map(|l| serde_json::from_str(l).unwrap())
            .collect();
        fs::write(
            out.join(format!("controls_{mode}.svg")),
            controls_svg(&rows, &mode),
        )
        .unwrap();
        eprintln!("wrote controls_{mode}.svg");
        return;
    }
    if mode == "v3" || mode == "quick3" {
        run_amendment2(&out, &mode);
        return;
    }
    let (n_e1, n_pos1, n_e2, n_e3, n_pos3, n_e4) = match mode.as_str() {
        "quick" => (40, 5, 3, 20, 5, 10),
        "v1" => (400, 20, 10, 200, 20, 100),
        "v2" => (2000, 20, 0, 1000, 20, 500),
        m => panic!("--protocol {m}: expected v1, v2, v3, quick or quick3"),
    };
    let run_e2 = n_e2 > 0;
    let mut rows = fs::File::create(out.join(format!("rows_{mode}.jsonl"))).unwrap();
    let mut md = String::new();
    let mut summary = Vec::new();
    let started = Instant::now();
    let mode_line = match mode.as_str() {
        "quick" => "QUICK smoke run at small N. Verdicts are meaningless (see the power column).",
        "v1" => "v1, `PROTOCOL.md` exactly as frozen. Its N fail the power check, so no rank test here can reach its α (see the power column).",
        _ => "v2, `PROTOCOL.md` as amended by `AMENDMENT_1.md`: larger N only. E2 is not re-run.",
    };
    writeln!(md, "# Controlled re-test: results ({mode})\n\n{mode_line} Seed base {SEED}. Protocol file hashes are recorded in `results_{mode}.json`. Every verdict below is the decision rule applied mechanically.\n").unwrap();

    // ── E1 ──
    let alpha1 = 0.05 / 18.0;
    writeln!(md, "## E1 — small-factor content (M1 = log2 of 2^32-smooth part found by factor(x, 2^18); M2 = # distinct such primes)\n").unwrap();
    writeln!(md, "Null: N = {n_e1} random primes n' in the Hasse interval per curve. Bonferroni α = 0.05/18 = {:.4}. Power check: min attainable p = 2/(N+1) = {:.4}, {}.\n", alpha1, 2.0 / (n_e1 as f64 + 1.0), if 2.0 / (n_e1 as f64 + 1.0) < alpha1 { "passes" } else { "**FAILS (no target can be declared unusual)**" }).unwrap();
    writeln!(md, "| curve | value | metric | target | control median [2.5%, 97.5%] | rank / N | p (two-sided) | verdict |\n|:--|:--|:--|--:|:--|--:|--:|:--|").unwrap();
    let mut pos_lines = Vec::new();
    for (ci, c) in CURVES.iter().enumerate() {
        let p = big(c.p);
        let n = big(c.n);
        let mut rng = StdRng::seed_from_u64(SEED + ci as u64);
        let w = (BigInt::from(4) * &p).sqrt();
        let lo = &p + 1 - &w;
        let hi = &p + 1 + &w;
        let mut ctrl: [Vec<f64>; 6] = Default::default();
        // Draw the whole null sample first (sequential, seeded), then run the
        // instrument in parallel; collect() keeps sample order.
        let nps: Vec<BigInt> = (0..n_e1)
            .map(|_| random_prime_in(&mut rng, &lo, &hi))
            .collect();
        let vals: Vec<BigInt> = nps
            .iter()
            .flat_map(|np| [np - 1, np + 1, BigInt::from(2) * &p + 2 - np])
            .collect();
        let mets: Vec<(f64, f64)> = vals.par_iter().map(small_factor_metrics).collect();
        for i in 0..n_e1 {
            for k in 0..3 {
                let (m1, m2) = mets[3 * i + k];
                ctrl[2 * k].push(m1);
                ctrl[2 * k + 1].push(m2);
                writeln!(rows, "{}", json!({"exp":"E1","arm":"control","curve":c.name,"i":i,"value":VN[k],"M1":m1,"M2":m2})).unwrap();
            }
        }
        let targets = [&n - 1, &n + 1, BigInt::from(2) * &p + 2 - &n];
        let tmets: Vec<(f64, f64)> = targets.par_iter().map(small_factor_metrics).collect();
        for (k, &(m1, m2)) in tmets.iter().enumerate() {
            writeln!(
                rows,
                "{}",
                json!({"exp":"E1","arm":"target","curve":c.name,"value":VN[k],"M1":m1,"M2":m2})
            )
            .unwrap();
            for (mi, (name, t)) in [("M1", m1), ("M2", m2)].iter().enumerate() {
                let cv = &ctrl[2 * k + mi];
                let (r, pv) = rank_p(*t, cv);
                let verdict = if pv < alpha1 {
                    "**unusual**"
                } else {
                    "within control distribution"
                };
                writeln!(
                    md,
                    "| {} | {} | {} | {:.1} | {:.1} [{:.1}, {:.1}] | {}/{} | {:.4} | {} |",
                    c.name,
                    ["n − 1", "n + 1", "twist"][k],
                    name,
                    t,
                    quantile(cv, 0.5),
                    quantile(cv, 0.025),
                    quantile(cv, 0.975),
                    r,
                    cv.len(),
                    pv,
                    verdict
                )
                .unwrap();
                summary.push(json!({"exp":"E1","curve":c.name,"value":VN[k],"metric":name,"target":t,"median":quantile(cv,0.5),"q025":quantile(cv,0.025),"q975":quantile(cv,0.975),"rank":r,"N":cv.len(),"p":pv,"alpha":alpha1,"unusual":pv<alpha1}));
            }
        }
        // positive control
        let q975 = quantile(&ctrl[0], 0.975);
        let mut hits = 0;
        let planted: Vec<BigInt> = (0..n_pos1)
            .map(|_| planted_smooth(&mut rng, n.bits()))
            .collect();
        let pmets: Vec<(f64, f64)> = planted.par_iter().map(small_factor_metrics).collect();
        for &(m1, _) in &pmets {
            writeln!(
                rows,
                "{}",
                json!({"exp":"E1","arm":"positive","curve":c.name,"M1":m1})
            )
            .unwrap();
            if m1 > q975 {
                hits += 1;
            }
        }
        pos_lines.push(format!("| {} | {}/{} planted-smooth values above the control 97.5th percentile of M1 ({:.1}) |", c.name, hits, n_pos1, q975));
        summary
            .push(json!({"exp":"E1","arm":"positive","curve":c.name,"detected":hits,"of":n_pos1}));
        eprintln!(
            "E1 {} done ({:.0}s)",
            c.name,
            started.elapsed().as_secs_f64()
        );
    }
    writeln!(md, "\nPositive control (planted 2^20-smooth part of ≈120 bits):\n\n| curve | detection |\n|:--|:--|\n{}\n", pos_lines.join("\n")).unwrap();

    // ── E2 ──
    if run_e2 {
        writeln!(
            md,
            "## E2 — census power (|x| ≤ 10^6; x = u/d², 2 ≤ d ≤ 40, |u| ≤ 10^5)\n"
        )
        .unwrap();
        writeln!(md, "| curve | planted integer found | planted d=2 found | planted d=3 found | random-b curves with any find | NIST b finds |\n|:--|--:|--:|--:|--:|--:|").unwrap();
        for (ci, c) in CURVES.iter().enumerate() {
            let p = big(c.p);
            let mut rng = StdRng::seed_from_u64(SEED + 100 + ci as u64);
            let sp = p.sqrt();
            // Draw every curve first (sequential, seeded, same order as the protocol:
            // planted integer, planted d = 2, planted d = 3, then random b), then run
            // the census in parallel; collect() keeps order.
            let mut planted: Vec<(usize, BigInt, (i64, i64))> = Vec::new();
            for kind in 0..3 {
                for _ in 0..n_e2 {
                    let (b, want) = loop {
                        if kind == 0 {
                            let x0: i64 = rng.gen_range(-1_000_000..=1_000_000);
                            let y0 = rng.gen_bigint_range(&BigInt::zero(), &sp);
                            let xb = BigInt::from(x0);
                            break (
                                &y0 * &y0 - &xb * &xb * &xb + BigInt::from(3) * &xb,
                                (x0, 1i64),
                            );
                        }
                        let d: i64 = if kind == 1 { 2 } else { 3 };
                        let u: i64 = rng.gen_range(-100_000..=100_000);
                        if u.gcd(&d) != 1 {
                            continue;
                        }
                        let d4 = BigInt::from(d.pow(4));
                        let d6 = BigInt::from(d.pow(6));
                        let ub = BigInt::from(u);
                        let vmax = (&p * &d6).sqrt();
                        let v = rng.gen_bigint_range(&BigInt::zero(), &vmax);
                        let num = &v * &v - &ub * &ub * &ub + BigInt::from(3) * &ub * &d4;
                        if num.mod_floor(&d6).is_zero() {
                            break (num / &d6, (u, d));
                        }
                    };
                    planted.push((kind, b, want));
                }
            }
            let negs: Vec<BigInt> = (0..n_e2)
                .map(|_| rng.gen_bigint_range(&BigInt::zero(), &p))
                .collect();
            let pos_res: Vec<_> = planted
                .par_iter()
                .map(|(_, b, _)| census(b, 1_000_000, 40, 100_000))
                .collect();
            let neg_res: Vec<_> = negs
                .par_iter()
                .map(|b| census(b, 1_000_000, 40, 100_000))
                .collect();
            let mut found = [0usize; 3];
            for ((kind, _, want), (ints, rats)) in planted.iter().zip(&pos_res) {
                let hit = if want.1 == 1 {
                    ints.contains(&want.0)
                } else {
                    rats.contains(want)
                };
                if hit {
                    found[*kind] += 1;
                }
                writeln!(rows, "{}", json!({"exp":"E2","arm":"positive","curve":c.name,"kind":kind,"planted":[want.0,want.1],"found":hit})).unwrap();
            }
            let mut neg_hits = 0;
            for (ints, rats) in &neg_res {
                if !ints.is_empty() || !rats.is_empty() {
                    neg_hits += 1;
                }
                writeln!(
                rows,
                "{}",
                json!({"exp":"E2","arm":"negative","curve":c.name,"finds":ints.len()+rats.len()})
            )
            .unwrap();
            }
            let (ints, rats) = census(&bighex(c.b_hex), 1_000_000, 40, 100_000);
            writeln!(
                md,
                "| {} | {}/{} | {}/{} | {}/{} | {}/{} | {} |",
                c.name,
                found[0],
                n_e2,
                found[1],
                n_e2,
                found[2],
                n_e2,
                neg_hits,
                n_e2,
                ints.len() + rats.len()
            )
            .unwrap();
            summary.push(json!({"exp":"E2","curve":c.name,"pos_int":found[0],"pos_d2":found[1],"pos_d3":found[2],"of":n_e2,"neg_hits":neg_hits,"nist_finds":ints.len()+rats.len()}));
            eprintln!(
                "E2 {} done ({:.0}s)",
                c.name,
                started.elapsed().as_secs_f64()
            );
        }
    }

    // ── E3 ──
    let alpha3 = 0.05 / 6.0;
    writeln!(md, "\n## E3 — ĥ of the one-point lift: generator vs random points (nats; ĥ := lim h_x([2^k]P)/4^k)\n").unwrap();
    writeln!(md, "Bonferroni α = 0.05/6 = {:.4}. Power check: min attainable p = 2/(N+1) = {:.4}, {}.\n\n| curve | metric | G̃ | random points median [2.5%, 97.5%] (N) | rank | p | verdict | small-x positive control median [min, max] |\n|:--|:--|--:|:--|--:|--:|:--|:--|", alpha3, 2.0 / (n_e3 as f64 + 1.0), if 2.0 / (n_e3 as f64 + 1.0) < alpha3 { "passes" } else { "**FAILS (no target can be declared unusual)**" }).unwrap();
    for (ci, c) in CURVES.iter().enumerate() {
        let p = big(c.p);
        let b = bighex(c.b_hex);
        let lp = ln_abs(&p);
        let mut rng = StdRng::seed_from_u64(SEED + 200 + ci as u64);
        let (mut r_h, mut r_d) = (Vec::new(), Vec::new());
        let pts: Vec<(BigInt, BigInt)> = (0..n_e3)
            .map(|_| random_point(&mut rng, &p, &b, None))
            .collect();
        let hs: Vec<(f64, f64)> = pts
            .par_iter()
            .map(|(x, y)| one_point_lift_height(x, y))
            .collect();
        for &(h, hx) in &hs {
            r_h.push(h / lp);
            r_d.push(hx - h);
            writeln!(rows, "{}", json!({"exp":"E3","arm":"control","curve":c.name,"hhat_over_lnp":h/lp,"hx_minus_hhat":hx-h})).unwrap();
        }
        let (mut s_h, mut s_d) = (Vec::new(), Vec::new());
        let xmax = BigInt::one() << 32;
        let spts: Vec<(BigInt, BigInt)> = (0..n_pos3)
            .map(|_| random_point(&mut rng, &p, &b, Some(&xmax)))
            .collect();
        let shs: Vec<(f64, f64)> = spts
            .par_iter()
            .map(|(x, y)| one_point_lift_height(x, y))
            .collect();
        for &(h, hx) in &shs {
            s_h.push(h / lp);
            s_d.push(hx - h);
            writeln!(rows, "{}", json!({"exp":"E3","arm":"positive_small_x","curve":c.name,"hhat_over_lnp":h/lp,"hx_minus_hhat":hx-h,"hx_over_lnp":hx/lp})).unwrap();
        }
        let (h, hx) = one_point_lift_height(&bighex(c.gx_hex), &bighex(c.gy_hex));
        writeln!(rows, "{}", json!({"exp":"E3","arm":"target","curve":c.name,"hhat_over_lnp":h/lp,"hx_minus_hhat":hx-h})).unwrap();
        for (name, t, cv, sv) in [
            ("ĥ/ln p", h / lp, &r_h, &s_h),
            ("h_x − ĥ", hx - h, &r_d, &s_d),
        ] {
            let (r, pv) = rank_p(t, cv);
            let verdict = if pv < alpha3 {
                "**unusual**"
            } else {
                "within control distribution"
            };
            writeln!(md, "| {} | {} | {:.4} | {:.4} [{:.4}, {:.4}] ({}) | {}/{} | {:.3} | {} | {:.4} [{:.4}, {:.4}] |", c.name, name, t, quantile(cv, 0.5), quantile(cv, 0.025), quantile(cv, 0.975), cv.len(), r, cv.len(), pv, verdict, quantile(sv, 0.5), quantile(sv, 0.0), quantile(sv, 1.0)).unwrap();
            summary.push(json!({"exp":"E3","curve":c.name,"metric":name,"target":t,"median":quantile(cv,0.5),"q025":quantile(cv,0.025),"q975":quantile(cv,0.975),"rank":r,"N":cv.len(),"p":pv,"alpha":alpha3,"pos_median":quantile(sv,0.5)}));
        }
        eprintln!(
            "E3 {} done ({:.0}s)",
            c.name,
            started.elapsed().as_secs_f64()
        );
    }

    // ── E4 ──
    writeln!(
        md,
        "\n## E4 — two-point-lift pairing: (G̃, R̃) vs random pairs\n"
    )
    .unwrap();
    writeln!(md, "Bonferroni α = 0.05/6 = {:.4}. Power check: min attainable p = 2/(N+1) = {:.4}, {}.\n\n| curve | metric | target | random pairs median [2.5%, 97.5%] (N) | rank | p | verdict | dependent-pair positive control ([2]G̃; [3]G̃) |\n|:--|:--|--:|:--|--:|--:|:--|:--|", alpha3, 2.0 / (n_e4 as f64 + 1.0), if 2.0 / (n_e4 as f64 + 1.0) < alpha3 { "passes" } else { "**FAILS (no target can be declared unusual)**" }).unwrap();
    for (ci, c) in CURVES.iter().enumerate() {
        let p = big(c.p);
        let b = bighex(c.b_hex);
        let mut rng = StdRng::seed_from_u64(SEED + 300 + ci as u64);
        let (mut cs, mut rs) = (Vec::new(), Vec::new());
        let pairs: Vec<_> = (0..n_e4)
            .map(|_| {
                (
                    random_point(&mut rng, &p, &b, None),
                    random_point(&mut rng, &p, &b, None),
                )
            })
            .collect();
        let prs: Vec<(f64, f64)> = pairs
            .par_iter()
            .map(|(p1, p2)| two_point_pairing(p1, p2))
            .collect();
        eprintln!(
            "E4 {} random pairs done ({:.0}s)",
            c.name,
            started.elapsed().as_secs_f64()
        );
        for &(cos, reg) in &prs {
            cs.push(cos);
            rs.push(reg);
            writeln!(
                rows,
                "{}",
                json!({"exp":"E4","arm":"control","curve":c.name,"abs_cos":cos,"norm_reg":reg})
            )
            .unwrap();
        }
        let g = (bighex(c.gx_hex), bighex(c.gy_hex));
        let r = ec_mul_fp(&p, &BigInt::from(-3), &big("1311768467294899695"), &g);
        // positive control: dependent rational points on E_1 through G̃
        let bq = &g.1 * &g.1 - &g.0 * &g.0 * &g.0 + BigInt::from(3) * &g.0;
        let e1 = CurveQ {
            a: Rat::small(-3),
            b: Rat::int(&bq),
        };
        let gt: Pt = Some((Rat::int(&g.0), Rat::int(&g.1)));
        let g2 = e1.add(&gt, &gt);
        let g3 = e1.add(&g2, &gt);
        assert!(e1.on_curve(&g2) && e1.on_curve(&g3));
        let ((tc, tr), ((pc2, pr2), (pc3, pr3))) = rayon::join(
            || two_point_pairing(&g, &r),
            || rayon::join(|| pairing_on(&e1, &gt, &g2), || pairing_on(&e1, &gt, &g3)),
        );
        writeln!(
            rows,
            "{}",
            json!({"exp":"E4","arm":"target","curve":c.name,"abs_cos":tc,"norm_reg":tr})
        )
        .unwrap();
        writeln!(
            rows,
            "{}",
            json!({"exp":"E4","arm":"positive","curve":c.name,"g2":[pc2,pr2],"g3":[pc3,pr3]})
        )
        .unwrap();
        for (name, t, cv, pos) in [
            ("|cos|", tc, &cs, format!("{:.6}; {:.6}", pc2, pc3)),
            ("det/(ĥĥ)", tr, &rs, format!("{:.2e}; {:.2e}", pr2, pr3)),
        ] {
            let (rk, pv) = rank_p(t, cv);
            let verdict = if pv < alpha3 {
                "**unusual**"
            } else {
                "within control distribution"
            };
            writeln!(
                md,
                "| {} | {} | {:.5} | {:.5} [{:.5}, {:.5}] ({}) | {}/{} | {:.3} | {} | {} |",
                c.name,
                name,
                t,
                quantile(cv, 0.5),
                quantile(cv, 0.025),
                quantile(cv, 0.975),
                cv.len(),
                rk,
                cv.len(),
                pv,
                verdict,
                pos
            )
            .unwrap();
            summary.push(json!({"exp":"E4","curve":c.name,"metric":name,"target":t,"median":quantile(cv,0.5),"q025":quantile(cv,0.025),"q975":quantile(cv,0.975),"rank":rk,"N":cv.len(),"p":pv,"alpha":alpha3,"positive":pos}));
        }
        eprintln!(
            "E4 {} done ({:.0}s)",
            c.name,
            started.elapsed().as_secs_f64()
        );
    }

    writeln!(
        md,
        "\n---\nRun time {:.0} s. Raw per-sample rows: `rows_{mode}.jsonl`.",
        started.elapsed().as_secs_f64()
    )
    .unwrap();
    let hash_of = |name: &str| -> String {
        let bytes = fs::read(format!("research/nist_cyclotomic_controls_20261011/{name}"))
            .expect("protocol file");
        crypto_lib::hash::sha256::sha256(&bytes)
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect()
    };
    let mut hashes = json!({"PROTOCOL.md": hash_of("PROTOCOL.md")});
    if mode == "v2" {
        hashes["AMENDMENT_1.md"] = json!(hash_of("AMENDMENT_1.md"));
    }
    fs::write(out.join(format!("RESULTS_{mode}.md")), &md).unwrap();
    fs::write(
        out.join(format!("results_{mode}.json")),
        serde_json::to_string_pretty(&json!({
            "protocol_sha256": hashes, "seed": SEED, "mode": mode,
            "N": {"E1": n_e1, "E1_pos": n_pos1, "E2": n_e2, "E3": n_e3, "E3_pos": n_pos3, "E4": n_e4},
            "summary": summary,
        }))
        .unwrap(),
    )
    .unwrap();
    eprintln!("wrote {}", out.display());
}
