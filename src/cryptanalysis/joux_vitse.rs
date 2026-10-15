//! **Joux–Vitse `(n − 1)`-decomposition at `n = 5`: one trial, measured
//! and checked, on curves over `F_{p⁵}`.**
//!
//! Companion to `research/joux_vitse_n1_20261015/`.  Gaudry and Diem
//! decompose a point of `E(F_{qⁿ})` into `n` points of the factor base
//! `F = {P : x(P) ∈ F_q}`; the Weil-descended system is square and its
//! solve cost grows as `2^{O(n²)}`.  Joux and Vitse (Vitse, PhD thesis 2011,
//! §7.3; J. Cryptology 2013) decompose into `n − 1` points instead:
//! `S_n(x₁, …, x_{n−1}, x_R)`, symmetrised in the elementary symmetric
//! functions `e₁, …, e_{n−1}` of the abscissae and Weil-restricted, is `n`
//! equations of total degree `2^{n−2}` in `n − 1` unknowns over `F_q`.  The
//! system is overdetermined: for most `R` it is inconsistent and its
//! Gröbner basis is `{1}`, at the price of a success probability of about
//! `1 / ((n − 1)! q)` per trial.
//!
//! At `n = 5` this module builds, for a curve `y² = x³ + ax + b` over
//! `F_{p⁵} = F_p[t]/(t⁵ − c)` (`p ≡ 1 mod 5`):
//!
//! * `H(e₁, …, e₄; x_R) = S₅(x₁, …, x₄, x_R)` on the 495 monomials of total
//!   degree `≤ 8` in the `e`'s, by interpolation from resultant evaluations
//!   of `S₅` (as [`super::gaudry_quartic::SymmetrisedS5`] does over
//!   `F_{p⁴}`), checked on fresh points ([`SymmetrisedS5x5`]);
//! * per trial, the five `F_p`-components of `H(e; x_R)`: five equations of
//!   degree 8 in four unknowns ([`SymmetrisedS5x5::weil_restrict`]);
//! * the trial itself ([`trial`]): degree-bounded F4 over `F_p`
//!   ([`super::f4_fp`]) at the Macaulay bound `5·7 + 1 = 36`, optionally
//!   replaying a recorded trace (Joux–Vitse's "F4Remake",
//!   [`super::f4_fp::f4_replay`]); a consistent system is solved, every
//!   `F_p`-solution `e` is split as `X⁴ − e₁X³ + e₂X² − e₃X + e₄`, and
//!   every split solution whose roots lift to base points is turned into
//!   a relation by trying the sign patterns and **re-adding the points in
//!   the group**;
//! * the independent check ([`mitm_decompositions`]): every decomposition
//!   of `R` into four distinct base points, by meet in the middle over
//!   pairs, without any polynomial system.  On every trial the two sets of
//!   decompositions must be equal.
//!
//! Costs are counted in `F_p` multiplications: the F4 row reductions by
//! [`super::f4_fp`]'s counter, everything else by [`Fp5`]'s.  The module
//! builds no linear algebra of its own beyond the small dense solve that a
//! complete DLP at toy size needs ([`solve_logs`]).

use std::cell::Cell;
use std::collections::{BTreeSet, HashMap};
use std::time::Instant;

use rand::rngs::StdRng;
use rand::Rng;
use rayon::prelude::*;
use serde::Serialize;

use super::f4_fp::{self, F4Options, F4Trace, Ordering, Poly, Verdict};

#[inline]
fn mm(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}
#[inline]
fn am(a: u64, b: u64, p: u64) -> u64 {
    let s = a + b;
    if s >= p {
        s - p
    } else {
        s
    }
}
#[inline]
fn sm(a: u64, b: u64, p: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a + p - b
    }
}
fn pow_mod(mut b: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1 % p;
    b %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mm(r, b, p);
        }
        b = mm(b, b, p);
        e >>= 1;
    }
    r
}
fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    let mut d = 2;
    while d * d <= n {
        if n % d == 0 {
            return false;
        }
        d += 1;
    }
    true
}

// ── F_{p⁵} = F_p[t] / (t⁵ − c) ───────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default, PartialOrd, Ord, Serialize)]
pub struct E5(pub [u64; 5]);

impl E5 {
    pub const ZERO: E5 = E5([0; 5]);
    pub const ONE: E5 = E5([1, 0, 0, 0, 0]);
    pub fn is_zero(&self) -> bool {
        self.0 == [0; 5]
    }
    pub fn in_fp(&self) -> bool {
        self.0[1..] == [0; 4]
    }
}

/// `F_{p⁵}` with a counter of `F_p` multiplications.  Charges: `29` per
/// product (25 schoolbook products and 4 reductions by `c`), `5` per
/// multiplication by an element of `F_p`, `4` per Frobenius (the basis is
/// an eigenbasis: `t^p = w t`, `w = c^{(p−1)/5}`), and an inversion by
/// Itoh–Tsujii through the norm (two products, three Frobenius maps, the
/// norm, one `F_p` inversion charged at its square-and-multiply count).
#[derive(Clone)]
pub struct Fp5 {
    pub p: u64,
    pub c: u64,
    /// `w^i`, `i = 0..5`: the Frobenius acts as `a_i ↦ w^i a_i`.
    frob: [u64; 5],
    muls: Cell<u64>,
}

impl Fp5 {
    pub fn new(p: u64) -> Fp5 {
        assert!(
            p % 5 == 1 && is_prime(p) && p < 1 << 12,
            "t⁵ − c is irreducible over F_p for some c only when p ≡ 1 (mod 5); p < 2^12 keeps p⁵ in u64"
        );
        // t⁵ − c is irreducible iff c is not a fifth power (p ≡ 1 mod 5).
        let c = (2..p)
            .find(|&c| pow_mod(c, (p - 1) / 5, p) != 1)
            .expect("a non-fifth-power");
        let w = pow_mod(c, (p - 1) / 5, p);
        let frob = core::array::from_fn(|i| pow_mod(w, i as u64, p));
        Fp5 {
            p,
            c,
            frob,
            muls: Cell::new(0),
        }
    }
    #[inline]
    fn count(&self, k: u64) {
        self.muls.set(self.muls.get() + k);
    }
    pub fn muls(&self) -> u64 {
        self.muls.get()
    }
    pub fn from_fp(&self, x: u64) -> E5 {
        E5([x % self.p, 0, 0, 0, 0])
    }
    pub fn add(&self, a: &E5, b: &E5) -> E5 {
        E5(core::array::from_fn(|i| am(a.0[i], b.0[i], self.p)))
    }
    pub fn sub(&self, a: &E5, b: &E5) -> E5 {
        E5(core::array::from_fn(|i| sm(a.0[i], b.0[i], self.p)))
    }
    pub fn neg(&self, a: &E5) -> E5 {
        E5(core::array::from_fn(|i| (self.p - a.0[i]) % self.p))
    }
    pub fn scale(&self, a: &E5, k: u64) -> E5 {
        self.count(5);
        E5(core::array::from_fn(|i| mm(a.0[i], k, self.p)))
    }
    pub fn mul(&self, a: &E5, b: &E5) -> E5 {
        let p = self.p;
        self.count(29);
        let mut d = [0u128; 9];
        for i in 0..5 {
            for j in 0..5 {
                d[i + j] += a.0[i] as u128 * b.0[j] as u128;
            }
        }
        let r: [u64; 9] = core::array::from_fn(|i| (d[i] % p as u128) as u64);
        E5([
            am(r[0], mm(self.c, r[5], p), p),
            am(r[1], mm(self.c, r[6], p), p),
            am(r[2], mm(self.c, r[7], p), p),
            am(r[3], mm(self.c, r[8], p), p),
            r[4],
        ])
    }
    pub fn sq(&self, a: &E5) -> E5 {
        self.mul(a, a)
    }
    pub fn frobenius(&self, a: &E5) -> E5 {
        self.count(4);
        E5(core::array::from_fn(|i| mm(a.0[i], self.frob[i], self.p)))
    }
    fn frob2(&self, a: &E5) -> E5 {
        self.frobenius(&self.frobenius(a))
    }
    /// `N(a) = a^{1 + p + p² + p³ + p⁴} ∈ F_p` and `a^{p + p² + p³ + p⁴}`.
    fn norm_and_cofactor(&self, a: &E5) -> (u64, E5) {
        let u = self.frobenius(a); // a^p
        let u2 = self.mul(&u, &self.frobenius(&u)); // a^{p+p²}
        let u4 = self.mul(&u2, &self.frob2(&u2)); // a^{p+p²+p³+p⁴}
        let p = self.p;
        // the constant coefficient of a·u4: 5 products, 4 of them by c
        self.count(9);
        let mut n = mm(a.0[0], u4.0[0], p);
        for i in 1..5 {
            n = am(n, mm(self.c, mm(a.0[i], u4.0[5 - i], p), p), p);
        }
        (n, u4)
    }
    pub fn norm(&self, a: &E5) -> u64 {
        self.norm_and_cofactor(a).0
    }
    fn inv_fp(&self, x: u64) -> u64 {
        let p = self.p;
        let e = p - 2;
        self.count(2 * (64 - e.leading_zeros()) as u64);
        pow_mod(x, e, p)
    }
    pub fn inv(&self, a: &E5) -> E5 {
        let (n, u4) = self.norm_and_cofactor(a);
        assert!(n != 0, "inverse of zero in F_p⁵");
        self.scale(&u4, self.inv_fp(n))
    }
    pub fn pow(&self, a: &E5, mut e: u128) -> E5 {
        let mut base = *a;
        let mut acc = E5::ONE;
        while e > 0 {
            if e & 1 == 1 {
                acc = self.mul(&acc, &base);
            }
            base = self.sq(&base);
            e >>= 1;
        }
        acc
    }
    fn order(&self) -> u128 {
        (self.p as u128).pow(5) - 1
    }
    /// Quadratic character through the norm: for an extension of odd
    /// degree, `a` is a square in `F_{p⁵}` iff `N(a)` is one in `F_p`.
    pub fn is_square(&self, a: &E5) -> bool {
        if a.is_zero() {
            return true;
        }
        let n = self.norm(a);
        self.count(2 * 64);
        pow_mod(n, (self.p - 1) / 2, self.p) == 1
    }
    pub fn random(&self, rng: &mut StdRng) -> E5 {
        E5(core::array::from_fn(|_| rng.gen_range(0..self.p)))
    }
    /// Tonelli–Shanks over `F_{p⁵}^×`, with a fixed non-residue so that the
    /// root returned for a given square is deterministic.
    pub fn sqrt(&self, a: &E5) -> Option<E5> {
        if a.is_zero() {
            return Some(E5::ZERO);
        }
        if !self.is_square(a) {
            return None;
        }
        let q1 = self.order();
        let s = q1.trailing_zeros();
        let odd = q1 >> s;
        // a non-square of F_p stays one in F_{p⁵} (odd degree)
        let z = self.from_fp(
            (2..self.p)
                .find(|&z| pow_mod(z, (self.p - 1) / 2, self.p) == self.p - 1)
                .unwrap(),
        );
        let mut m = s;
        let mut cc = self.pow(&z, odd);
        let mut t = self.pow(a, odd);
        let mut r = self.pow(a, odd.div_ceil(2));
        while t != E5::ONE {
            let mut i = 0;
            let mut tt = t;
            while tt != E5::ONE {
                tt = self.sq(&tt);
                i += 1;
            }
            let mut b = cc;
            for _ in 0..(m - i - 1) {
                b = self.sq(&b);
            }
            m = i;
            cc = self.sq(&b);
            t = self.mul(&t, &cc);
            r = self.mul(&r, &b);
        }
        Some(r)
    }
}

// ── The curve `y² = x³ + a x + b` over `F_{p⁵}` ─────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize)]
pub struct Pt5 {
    pub x: E5,
    pub y: E5,
    pub inf: bool,
}

impl Pt5 {
    pub const INF: Pt5 = Pt5 {
        x: E5::ZERO,
        y: E5::ZERO,
        inf: true,
    };
}

#[derive(Clone)]
pub struct Curve5 {
    pub f: Fp5,
    pub a: E5,
    pub b: E5,
    ops: Cell<u64>,
}

impl Curve5 {
    /// Coefficients with `j(E) ∉ F_p`, so that the curve is not a twist of
    /// one defined over the base field (`F_{p⁵}` has no other subfield).
    pub fn random(p: u64, rng: &mut StdRng) -> Curve5 {
        let f = Fp5::new(p);
        loop {
            let a = f.random(rng);
            let b = f.random(rng);
            let a3 = f.mul(&f.sq(&a), &a);
            let disc = f.add(&f.scale(&a3, 4), &f.scale(&f.sq(&b), 27));
            if disc.is_zero() || a.is_zero() {
                continue;
            }
            // j = 1728 · 4a³ / disc
            let j = f.mul(&f.scale(&a3, 4), &f.inv(&disc));
            if j.in_fp() {
                continue;
            }
            return Curve5::with_params(p, a, b);
        }
    }
    pub fn with_params(p: u64, a: E5, b: E5) -> Curve5 {
        Curve5 {
            f: Fp5::new(p),
            a,
            b,
            ops: Cell::new(0),
        }
    }
    pub fn ops(&self) -> u64 {
        self.ops.get()
    }
    pub fn rhs(&self, x: &E5) -> E5 {
        let f = &self.f;
        f.add(&f.mul(&f.add(&f.sq(x), &self.a), x), &self.b)
    }
    pub fn lift_x(&self, x: &E5) -> Option<Pt5> {
        let y = self.f.sqrt(&self.rhs(x))?;
        Some(Pt5 {
            x: *x,
            y,
            inf: false,
        })
    }
    pub fn on_curve(&self, pt: &Pt5) -> bool {
        pt.inf || self.f.sq(&pt.y) == self.rhs(&pt.x)
    }
    pub fn neg(&self, pt: &Pt5) -> Pt5 {
        if pt.inf {
            return *pt;
        }
        Pt5 {
            x: pt.x,
            y: self.f.neg(&pt.y),
            inf: false,
        }
    }
    /// Affine addition: one inversion, two products and a squaring.
    pub fn add(&self, p1: &Pt5, p2: &Pt5) -> Pt5 {
        let f = &self.f;
        self.ops.set(self.ops.get() + 1);
        if p1.inf {
            return *p2;
        }
        if p2.inf {
            return *p1;
        }
        let lambda = if p1.x == p2.x {
            if p1.y != p2.y || p1.y.is_zero() {
                return Pt5::INF;
            }
            let x2 = f.sq(&p1.x);
            let num = f.add(&f.add(&f.add(&x2, &x2), &x2), &self.a);
            f.mul(&num, &f.inv(&f.add(&p1.y, &p1.y)))
        } else {
            f.mul(&f.sub(&p2.y, &p1.y), &f.inv(&f.sub(&p2.x, &p1.x)))
        };
        let x3 = f.sub(&f.sub(&f.sq(&lambda), &p1.x), &p2.x);
        let y3 = f.sub(&f.mul(&lambda, &f.sub(&p1.x, &x3)), &p1.y);
        Pt5 {
            x: x3,
            y: y3,
            inf: false,
        }
    }
    pub fn sub(&self, p1: &Pt5, p2: &Pt5) -> Pt5 {
        self.add(p1, &self.neg(p2))
    }
    pub fn mul(&self, pt: &Pt5, mut k: u128) -> Pt5 {
        let mut acc = Pt5::INF;
        let mut base = *pt;
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(&acc, &base);
            }
            base = self.add(&base, &base);
            k >>= 1;
        }
        acc
    }
    /// `#E(F_{p⁵})` by summing the quadratic character over every `x`
    /// (`p⁵` norms; feasible for `p ≤ 41` or so).
    pub fn order_by_count(&self) -> u128 {
        let f = &self.f;
        let p = f.p;
        let q = (p as u128).pow(5);
        let (ca, cb) = (self.a, self.b);
        // enumerate x by its five coordinates, in parallel over the first
        let total: i128 = (0..p)
            .into_par_iter()
            .map(|x0| {
                let c = Curve5::with_params(p, ca, cb);
                let mut s: i128 = 0;
                for idx in 0..p.pow(4) {
                    let mut v = idx;
                    let mut x = [x0, 0, 0, 0, 0];
                    for k in 1..5 {
                        x[k] = v % p;
                        v /= p;
                    }
                    let r = c.rhs(&E5(x));
                    if r.is_zero() {
                        continue;
                    }
                    let n = c.f.norm(&r);
                    s += if pow_mod(n, (p - 1) / 2, p) == 1 { 1 } else { -1 };
                }
                s
            })
            .sum();
        (q as i128 + 1 + total) as u128
    }
}

// ── Summation polynomials, numerically ───────────────────────────────────

fn poly_mul(f: &Fp5, a: &[E5], b: &[E5]) -> Vec<E5> {
    let mut out = vec![E5::ZERO; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = f.add(&out[i + j], &f.mul(x, y));
        }
    }
    out
}
fn poly_sub(f: &Fp5, a: &[E5], b: &[E5]) -> Vec<E5> {
    let n = a.len().max(b.len());
    (0..n)
        .map(|i| f.sub(a.get(i).unwrap_or(&E5::ZERO), b.get(i).unwrap_or(&E5::ZERO)))
        .collect()
}
fn poly_scale(f: &Fp5, a: &[E5], k: &E5) -> Vec<E5> {
    a.iter().map(|x| f.mul(x, k)).collect()
}

fn det(f: &Fp5, mut m: Vec<Vec<E5>>) -> E5 {
    let n = m.len();
    let mut d = E5::ONE;
    for c in 0..n {
        let Some(r) = (c..n).find(|&r| !m[r][c].is_zero()) else {
            return E5::ZERO;
        };
        if r != c {
            m.swap(r, c);
            d = f.neg(&d);
        }
        d = f.mul(&d, &m[c][c]);
        let inv = f.inv(&m[c][c]);
        for r in (c + 1)..n {
            if m[r][c].is_zero() {
                continue;
            }
            let k = f.mul(&m[r][c], &inv);
            for j in c..n {
                let v = f.mul(&k, &m[c][j]);
                m[r][j] = f.sub(&m[r][j], &v);
            }
        }
    }
    d
}

impl Curve5 {
    /// `S₃(x₁, x₂, X) = (x₁ − x₂)² X² − 2((x₁ + x₂)(x₁x₂ + a) + 2b) X +
    /// (x₁x₂ − a)² − 4b(x₁ + x₂)`, low to high in `X`.
    fn s3_in_last(&self, x1: &E5, x2: &E5) -> [E5; 3] {
        let f = &self.f;
        let s = f.add(x1, x2);
        let pr = f.mul(x1, x2);
        let two_b = f.add(&self.b, &self.b);
        let t = f.add(&f.mul(&s, &f.add(&pr, &self.a)), &two_b);
        let c1 = f.neg(&f.add(&t, &t));
        let u = f.sub(&pr, &self.a);
        let c0 = f.sub(&f.sq(&u), &f.mul(&f.add(&two_b, &two_b), &s));
        [c0, c1, f.sq(&f.sub(x1, x2))]
    }

    /// `S₃(x₁, x₂, x₃)`.
    pub fn s3(&self, x1: &E5, x2: &E5, x3: &E5) -> E5 {
        let f = &self.f;
        let [c0, c1, c2] = self.s3_in_last(x1, x2);
        f.add(&f.add(&f.mul(&f.mul(&c2, x3), x3), &f.mul(&c1, x3)), &c0)
    }

    /// `S₄(x₃, x₄, x₅, Y) = Res_Z(S₃(x₃, x₄, Z), S₃(x₅, Y, Z))` in `Y`.
    fn s4_in_last(&self, x3: &E5, x4: &E5, x5: &E5) -> Vec<E5> {
        let f = &self.f;
        let [a0, a1, a2] = self.s3_in_last(x3, x4);
        let (a, b) = (&self.a, &self.b);
        let x5s = f.sq(x5);
        let two = |v: &E5| f.add(v, v);
        let four = |v: &E5| two(&two(v));
        let b2 = vec![x5s, f.neg(&two(x5)), E5::ONE];
        let b1 = vec![
            f.neg(&two(&f.add(&f.mul(a, x5), &two(b)))),
            f.neg(&two(&f.add(&x5s, a))),
            f.neg(&two(x5)),
        ];
        let b0 = vec![
            f.sub(&f.sq(a), &four(&f.mul(b, x5))),
            f.neg(&f.add(&two(&f.mul(a, x5)), &four(b))),
            x5s,
        ];
        let p1 = poly_sub(f, &poly_scale(f, &b0, &a2), &poly_scale(f, &b2, &a0));
        let p2 = poly_sub(f, &poly_scale(f, &b1, &a2), &poly_scale(f, &b2, &a1));
        let p3 = poly_sub(f, &poly_scale(f, &b0, &a1), &poly_scale(f, &b1, &a0));
        poly_sub(f, &poly_mul(f, &p1, &p1), &poly_mul(f, &p2, &p3))
    }

    /// `S₅(x₁, …, x₅) = Res_Y(S₃(x₁, x₂, Y), S₄(x₃, x₄, x₅, Y))`, the
    /// `6 × 6` Sylvester determinant at formal degrees 2 and 4.
    pub fn s5(&self, x: &[E5; 5]) -> E5 {
        let f = &self.f;
        let a = self.s3_in_last(&x[0], &x[1]);
        let mut b = self.s4_in_last(&x[2], &x[3], &x[4]);
        b.resize(5, E5::ZERO);
        let mut m = vec![vec![E5::ZERO; 6]; 6];
        for i in 0..4 {
            for k in 0..3 {
                m[i][i + k] = a[2 - k];
            }
        }
        for i in 0..2 {
            for k in 0..5 {
                m[4 + i][i + k] = b[4 - k];
            }
        }
        det(f, m)
    }
}

// ── The symmetrised S₅ over F_{p⁵}, by interpolation ─────────────────────

/// Monomials `e₁^a e₂^b e₃^c e₄^d` of total degree `≤ d`.
pub fn exponents4(d: u8) -> Vec<[u8; 4]> {
    let mut v = Vec::new();
    for a in 0..=d {
        for b in 0..=(d - a) {
            for c in 0..=(d - a - b) {
                for e in 0..=(d - a - b - c) {
                    v.push([a, b, c, e]);
                }
            }
        }
    }
    v
}

fn elementary(f: &Fp5, x: &[E5; 4]) -> [E5; 4] {
    let e1 = f.add(&f.add(&x[0], &x[1]), &f.add(&x[2], &x[3]));
    let mut e2 = E5::ZERO;
    let mut e3 = E5::ZERO;
    for i in 0..4 {
        for j in (i + 1)..4 {
            let xij = f.mul(&x[i], &x[j]);
            e2 = f.add(&e2, &xij);
            for k in (j + 1)..4 {
                e3 = f.add(&e3, &f.mul(&xij, &x[k]));
            }
        }
    }
    let e4 = f.mul(&f.mul(&x[0], &x[1]), &f.mul(&x[2], &x[3]));
    [e1, e2, e3, e4]
}

fn mono_values(f: &Fp5, e: &[E5; 4], monos: &[[u8; 4]]) -> Vec<E5> {
    let mut pows = [[E5::ONE; 9]; 4];
    for i in 0..4 {
        for k in 1..9 {
            pows[i][k] = f.mul(&pows[i][k - 1], &e[i]);
        }
    }
    monos
        .iter()
        .map(|m| {
            f.mul(
                &f.mul(&pows[0][m[0] as usize], &pows[1][m[1] as usize]),
                &f.mul(&pows[2][m[2] as usize], &pows[3][m[3] as usize]),
            )
        })
        .collect()
}

/// `S₅(x₁, x₂, x₃, x₄, x_R) = H(e₁, …, e₄, x_R)`: `495 · 9` coefficients
/// in `F_{p⁵}` (total degree `≤ 8` in the `e`'s, degree `≤ 8` in `x_R`).
#[derive(Clone)]
pub struct SymmetrisedS5x5 {
    pub monos: Vec<[u8; 4]>,
    pub coef: Vec<[E5; 9]>,
    /// `F_p` multiplications spent building it (once per curve).
    pub precompute_muls: u64,
}

impl SymmetrisedS5x5 {
    fn coefficients_in_last(curve: &Curve5, x: &[E5; 4], vinv: &[[u64; 9]; 9]) -> [E5; 9] {
        let f = &curve.f;
        let vals: Vec<E5> = (0..9u64)
            .map(|k| curve.s5(&[x[0], x[1], x[2], x[3], f.from_fp(k)]))
            .collect();
        core::array::from_fn(|j| {
            vals.iter()
                .enumerate()
                .fold(E5::ZERO, |acc, (k, v)| f.add(&acc, &f.scale(v, vinv[j][k])))
        })
    }

    fn vandermonde_inverse(p: u64) -> [[u64; 9]; 9] {
        // Gauss–Jordan on [V | I] over F_p, V[k][j] = k^j, nodes 0..8.
        let mut m: Vec<Vec<u64>> = (0..9u64)
            .map(|k| {
                let mut row: Vec<u64> = (0..9).map(|j| pow_mod(k, j, p)).collect();
                row.extend((0..9).map(|i| u64::from(i == k)));
                row
            })
            .collect();
        for c in 0..9 {
            let r = (c..9).find(|&r| m[r][c] != 0).expect("p > 8");
            m.swap(r, c);
            let inv = pow_mod(m[c][c], p - 2, p);
            for v in m[c].iter_mut() {
                *v = mm(*v, inv, p);
            }
            let piv = m[c].clone();
            for (r, row) in m.iter_mut().enumerate() {
                if r != c && row[c] != 0 {
                    let k = row[c];
                    for (x, &y) in row.iter_mut().zip(&piv) {
                        *x = sm(*x, mm(k, y, p), p);
                    }
                }
            }
        }
        core::array::from_fn(|j| core::array::from_fn(|k| m[j][9 + k]))
    }

    pub fn precompute(curve: &Curve5, rng: &mut StdRng) -> SymmetrisedS5x5 {
        let f = &curve.f;
        let p = f.p;
        let before = f.muls();
        let vinv = Self::vandermonde_inverse(p);
        let monos = exponents4(8);
        let n = monos.len();
        loop {
            let mut sys: Vec<Vec<E5>> = Vec::with_capacity(n);
            for _ in 0..n {
                let x: [E5; 4] = core::array::from_fn(|_| f.random(rng));
                let e = elementary(f, &x);
                let mut row = mono_values(f, &e, &monos);
                row.extend(Self::coefficients_in_last(curve, &x, &vinv));
                sys.push(row);
            }
            let mut ok = true;
            for c in 0..n {
                let Some(r) = (c..n).find(|&r| !sys[r][c].is_zero()) else {
                    ok = false;
                    break;
                };
                sys.swap(r, c);
                let inv = f.inv(&sys[c][c]);
                let pivot: Vec<E5> = sys[c].iter().map(|v| f.mul(v, &inv)).collect();
                sys[c] = pivot.clone();
                let work: u64 = sys
                    .par_iter_mut()
                    .enumerate()
                    .map(|(r, row)| {
                        if r == c || row[c].is_zero() {
                            return 0;
                        }
                        let fl = Fp5::new(p);
                        let k = row[c];
                        for j in c..row.len() {
                            if !pivot[j].is_zero() {
                                let v = fl.mul(&k, &pivot[j]);
                                row[j] = fl.sub(&row[j], &v);
                            }
                        }
                        fl.muls()
                    })
                    .sum();
                f.count(work);
            }
            if !ok {
                continue;
            }
            let coef: Vec<[E5; 9]> = (0..n)
                .map(|r| core::array::from_fn(|j| sys[r][n + j]))
                .collect();
            let pre = SymmetrisedS5x5 {
                monos: monos.clone(),
                coef,
                precompute_muls: 0,
            };
            let good = (0..8).all(|_| {
                let x: [E5; 4] = core::array::from_fn(|_| f.random(rng));
                let want = Self::coefficients_in_last(curve, &x, &vinv);
                let e = elementary(f, &x);
                let mv = mono_values(f, &e, &pre.monos);
                (0..9).all(|j| {
                    let got = mv
                        .iter()
                        .zip(&pre.coef)
                        .fold(E5::ZERO, |acc, (m, c)| f.add(&acc, &f.mul(m, &c[j])));
                    got == want[j]
                })
            });
            assert!(good, "the symmetrised S₅ does not reproduce the resultant");
            return SymmetrisedS5x5 {
                precompute_muls: f.muls() - before,
                ..pre
            };
        }
    }

    /// `H(e; x_R)` evaluated at `e ∈ F_{p⁵}⁴` (used by the tests).
    pub fn eval(&self, f: &Fp5, e: &[E5; 4], x_r: &E5) -> E5 {
        let mv = mono_values(f, e, &self.monos);
        let mut pows = [E5::ONE; 9];
        for k in 1..9 {
            pows[k] = f.mul(&pows[k - 1], x_r);
        }
        mv.iter().zip(&self.coef).fold(E5::ZERO, |acc, (m, c)| {
            let v = (0..9).fold(E5::ZERO, |a, j| f.add(&a, &f.mul(&c[j], &pows[j])));
            f.add(&acc, &f.mul(m, &v))
        })
    }

    /// The five `F_p`-components of `H(e₁, …, e₄; x_R)`: the descended
    /// system, five polynomials of degree `≤ 8` in four unknowns, as
    /// [`f4_fp`] input (`ϕ₀, …, ϕ₄` of Vitse's (7.9)).
    pub fn weil_restrict(&self, f: &Fp5, x_r: &E5) -> Vec<Poly> {
        let mut pows = [E5::ONE; 9];
        for k in 1..9 {
            pows[k] = f.mul(&pows[k - 1], x_r);
        }
        let mut terms: Vec<Vec<(Vec<u32>, u64)>> = vec![Vec::new(); 5];
        for (m, c) in self.monos.iter().zip(&self.coef) {
            let v = (0..9).fold(E5::ZERO, |acc, j| f.add(&acc, &f.mul(&c[j], &pows[j])));
            let e: Vec<u32> = m.iter().map(|&x| x as u32).collect();
            for (i, comp) in terms.iter_mut().enumerate() {
                if v.0[i] != 0 {
                    comp.push((e.clone(), v.0[i]));
                }
            }
        }
        terms
            .iter()
            .map(|t| f4_fp::normalise(t, f.p, Ordering::Grevlex))
            .filter(|t| !t.is_empty())
            .collect()
    }
}

// ── Factor base, meet in the middle ──────────────────────────────────────

/// `F = {P : x(P) ∈ F_p}` modulo `±`: one point per abscissa whose
/// `x³ + ax + b` is a square in `F_{p⁵}` (the root [`Fp5::sqrt`] returns).
pub fn factor_base(curve: &Curve5) -> Vec<Pt5> {
    (0..curve.f.p)
        .filter_map(|x| curve.lift_x(&curve.f.from_fp(x)))
        .filter(|pt| !pt.y.is_zero())
        .collect()
}

/// The target-independent table of the meet-in-the-middle oracle:
/// `x(P_i ± P_j)` for every pair of distinct base points.
pub struct PairTable {
    table: HashMap<E5, Vec<(usize, usize)>>,
}

impl PairTable {
    pub fn new(curve: &Curve5, base: &[Pt5]) -> PairTable {
        let n = base.len();
        let mut table: HashMap<E5, Vec<(usize, usize)>> = HashMap::new();
        for i in 0..n {
            for j in (i + 1)..n {
                for q in [curve.add(&base[i], &base[j]), curve.sub(&base[i], &base[j])] {
                    if !q.inf {
                        table.entry(q.x).or_default().push((i, j));
                    }
                }
            }
        }
        PairTable { table }
    }
}

/// Every set of four distinct base abscissae `{x_i, x_j, x_k, x_l}` with
/// `±P_i ± P_j ± P_k ± P_l = R`: `x(R − s_k P_k − s_l P_l)` looked up among
/// the `x(P_i ± P_j)`.  No polynomial system is involved.
pub fn mitm_decompositions(
    curve: &Curve5,
    base: &[Pt5],
    pairs: &PairTable,
    r: &Pt5,
) -> BTreeSet<[u64; 4]> {
    let n = base.len();
    let mut out = BTreeSet::new();
    for k in 0..n {
        for l in (k + 1)..n {
            for (sk, sl) in [(1, 1), (1, -1), (-1, 1), (-1, -1)] {
                let pk = if sk == 1 { base[k] } else { curve.neg(&base[k]) };
                let pl = if sl == 1 { base[l] } else { curve.neg(&base[l]) };
                let v = curve.sub(&curve.sub(r, &pk), &pl);
                if v.inf {
                    continue;
                }
                if let Some(ps) = pairs.table.get(&v.x) {
                    for &(i, j) in ps {
                        if i == k || i == l || j == k || j == l {
                            continue;
                        }
                        let mut s = [base[i].x.0[0], base[j].x.0[0], base[k].x.0[0], base[l].x.0[0]];
                        s.sort_unstable();
                        out.insert(s);
                    }
                }
            }
        }
    }
    out
}

// ── One trial ────────────────────────────────────────────────────────────

/// A relation `R = Σ s_i P_{idx_i}`, checked by re-adding in the group.
#[derive(Clone, Debug, Serialize)]
pub struct Relation {
    /// `(factor-base index, sign)`, four terms (indices may repeat).
    pub terms: Vec<(usize, i8)>,
}

/// How a consistent system's `F_p`-solutions split.
#[derive(Clone, Debug, Default, Serialize)]
pub struct SolutionSplit {
    /// `e ∈ F_p⁴` solving the descended system.
    pub solutions: usize,
    /// ... whose quartic splits over `F_p` into four distinct roots.
    pub split_distinct: usize,
    /// ... which splits with a repeated root.
    pub split_repeated: usize,
    /// ... which does not split over `F_p` (roots in `F_{p²}`, `F_{p³}`,
    /// `F_{p⁴}`): a decomposition over an extension, not a relation.
    pub not_split: usize,
    /// Split solutions with a root whose point is not over `F_{p⁵}` (on
    /// the quadratic twist), or with no sign pattern re-adding to `R`.
    pub split_not_relation: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct TrialReport {
    pub consistent: bool,
    /// F4 (or replay) verdict could not be reached under the degree bound.
    pub undetermined: bool,
    /// The run was a trace replay that matched the recorded run.
    pub replayed: bool,
    /// A replay diverged and the full F4 was run after it.
    pub replay_mismatch: bool,
    /// `F_p` multiplications: building the descended system.
    pub muls_system: u64,
    /// `F_p` multiplications: F4 row reductions of the consistency test
    /// (the replay, plus the full run after a mismatch).
    pub muls_f4: u64,
    /// `F_p` multiplications: solving a consistent system (substitution
    /// tree of `f4_fp::solve`), splitting, lifting and re-adding.
    pub muls_solve: u64,
    pub solving_degree: u32,
    pub max_rows: usize,
    pub max_cols: usize,
    pub steps: usize,
    pub ms: f64,
    pub split: SolutionSplit,
    pub relations: Vec<Relation>,
    /// Sorted abscissae of the relations with four distinct base points
    /// (the sets [`mitm_decompositions`] must reproduce).
    pub decompositions: BTreeSet<[u64; 4]>,
}

/// Roots in `F_p` of `X⁴ − e₁X³ + e₂X² − e₃X + e₄`, with multiplicity.
fn quartic_roots(e: &[u64], p: u64) -> Vec<u64> {
    let mut c = vec![e[3] % p, (p - e[2] % p) % p, e[1] % p, (p - e[0] % p) % p, 1]; // low→high
    let mut roots = Vec::new();
    for x in 0..p {
        loop {
            if c.len() < 2 {
                break;
            }
            // synthetic division by (X − x)
            let deg = c.len() - 1;
            let mut q = vec![0u64; deg];
            let mut acc = 0u64;
            for k in (0..=deg).rev() {
                acc = am(mm(acc, x, p), c[k], p);
                if k > 0 {
                    q[k - 1] = acc;
                }
            }
            if acc == 0 {
                roots.push(x);
                c = q;
            } else {
                break;
            }
        }
    }
    roots
}

/// One decomposition trial on `R`: the descended system, its consistency
/// test by F4 (replaying `trace` when given), and on a consistent system
/// the solutions, their split, and the relations re-added in the group.
pub fn trial(
    curve: &Curve5,
    sym: &SymmetrisedS5x5,
    base: &[Pt5],
    base_index: &HashMap<u64, usize>,
    r: &Pt5,
    opts: &F4Options,
    trace: Option<&F4Trace>,
) -> TrialReport {
    let t0 = Instant::now();
    let f = &curve.f;
    let p = f.p;
    let m0 = f.muls();
    let system = sym.weil_restrict(f, &r.x);
    let muls_system = f.muls() - m0;
    let (mut replayed, mut replay_mismatch) = (false, false);
    let mut muls_f4 = 0u64;
    let rep = match trace {
        Some(tr) => {
            let rr = f4_fp::f4_replay(&system, 4, p, opts, tr);
            muls_f4 += rr.field_ops;
            if rr.trace_mismatch {
                replay_mismatch = true;
                let full = f4_fp::f4(&system, 4, p, opts);
                muls_f4 += full.field_ops;
                full
            } else {
                replayed = true;
                rr
            }
        }
        None => {
            let full = f4_fp::f4(&system, 4, p, opts);
            muls_f4 += full.field_ops;
            full
        }
    };
    let mut out = TrialReport {
        consistent: !rep.inconsistent,
        undetermined: false,
        replayed,
        replay_mismatch,
        muls_system,
        muls_f4,
        muls_solve: 0,
        solving_degree: rep.solving_degree,
        max_rows: rep.max_rows,
        max_cols: rep.max_cols,
        steps: rep.steps,
        ms: 0.0,
        split: SolutionSplit::default(),
        relations: Vec::new(),
        decompositions: BTreeSet::new(),
    };
    if rep.inconsistent {
        out.ms = t0.elapsed().as_secs_f64() * 1e3;
        return out;
    }
    // A consistent system: solve it (few solutions), split, lift, re-add.
    let m1 = f.muls();
    let sol = f4_fp::solve(&system, 4, p, opts);
    let mut solve_muls = sol.field_ops;
    match sol.verdict {
        Verdict::Inconsistent => {
            out.consistent = false;
        }
        Verdict::Undetermined => {
            out.undetermined = true;
        }
        Verdict::Solutions(sols) => {
            out.split.solutions = sols.len();
            for e in &sols {
                let roots = quartic_roots(e, p);
                solve_muls += 8 * p;
                if roots.len() < 4 {
                    out.split.not_split += 1;
                    continue;
                }
                let distinct: BTreeSet<u64> = roots.iter().copied().collect();
                if distinct.len() == 4 {
                    out.split.split_distinct += 1;
                } else {
                    out.split.split_repeated += 1;
                }
                let Some(idx): Option<Vec<usize>> =
                    roots.iter().map(|x| base_index.get(x).copied()).collect()
                else {
                    out.split.split_not_relation += 1;
                    continue;
                };
                let mut found = None;
                for signs in 0u8..16 {
                    let mut acc = Pt5::INF;
                    for (k, &i) in idx.iter().enumerate() {
                        let pt = if signs >> k & 1 == 0 {
                            base[i]
                        } else {
                            curve.neg(&base[i])
                        };
                        acc = curve.add(&acc, &pt);
                    }
                    if acc == *r {
                        found = Some(signs);
                        break;
                    }
                }
                match found {
                    Some(signs) => {
                        let terms: Vec<(usize, i8)> = idx
                            .iter()
                            .enumerate()
                            .map(|(k, &i)| (i, if signs >> k & 1 == 0 { 1 } else { -1 }))
                            .collect();
                        if distinct.len() == 4 {
                            let mut s = [roots[0], roots[1], roots[2], roots[3]];
                            s.sort_unstable();
                            out.decompositions.insert(s);
                        }
                        out.relations.push(Relation { terms });
                    }
                    None => out.split.split_not_relation += 1,
                }
            }
        }
    }
    solve_muls += f.muls() - m1;
    out.muls_solve = solve_muls;
    out.ms = t0.elapsed().as_secs_f64() * 1e3;
    out
}

/// Re-add a relation in the group: `Σ s_i P_i == R`.
pub fn verify_relation(curve: &Curve5, base: &[Pt5], rel: &Relation, r: &Pt5) -> bool {
    let mut acc = Pt5::INF;
    for &(i, s) in &rel.terms {
        let pt = if s == 1 { base[i] } else { curve.neg(&base[i]) };
        acc = curve.add(&acc, &pt);
    }
    acc == *r
}

/// The F4 options of a trial: grevlex, the Macaulay bound
/// `Σ (deg ϕ_i − 1) + 1 = 5·7 + 1 = 36` of Vitse's Hypothesis 7.3.4.
pub fn trial_options() -> F4Options {
    F4Options::new(Ordering::Grevlex, 36)
}

// ── A complete toy DLP: dense linear algebra modulo ℓ ────────────────────

fn inv_mod_u128(a: u128, m: u128) -> Option<u128> {
    let (mut t, mut nt) = (0i128, 1i128);
    let (mut r, mut nr) = (m as i128, (a % m) as i128);
    while nr != 0 {
        let q = r / nr;
        (t, nt) = (nt, t - q * nt);
        (r, nr) = (nr, r - q * nr);
    }
    if r != 1 {
        return None;
    }
    Some(t.rem_euclid(m as i128) as u128)
}
fn mulm(a: u128, b: u128, m: u128) -> u128 {
    // m < 2^64 in every use here
    a % m * (b % m) % m
}

/// Solve `Σ_i s_{ji} L_i − h b_j k ≡ h a_j (mod ℓ)` for the unknowns
/// `(L_1, …, L_m, k)` by Gauss–Jordan; returns `k` when it is determined.
/// The cost (multiplications modulo `ℓ`) is returned with it.
pub fn solve_logs(
    rows: &[(Vec<(usize, i8)>, u128, u128)],
    m: usize,
    h: u128,
    ell: u128,
) -> (Option<u128>, u64) {
    let cols = m + 1;
    let mut a: Vec<Vec<u128>> = rows
        .iter()
        .map(|(terms, aj, bj)| {
            let mut row = vec![0u128; cols + 1];
            for &(i, s) in terms {
                row[i] = if s == 1 {
                    (row[i] + 1) % ell
                } else {
                    (row[i] + ell - 1) % ell
                };
            }
            row[m] = (ell - mulm(h, *bj, ell)) % ell;
            row[cols] = mulm(h, *aj, ell);
            row
        })
        .collect();
    let mut ops = 0u64;
    let mut piv_col = Vec::new();
    let mut r = 0;
    for c in 0..cols {
        let Some(pr) = (r..a.len()).find(|&i| a[i][c] != 0) else {
            continue;
        };
        a.swap(r, pr);
        let inv = inv_mod_u128(a[r][c], ell).expect("ℓ prime");
        for v in a[r].iter_mut() {
            *v = mulm(*v, inv, ell);
        }
        ops += (cols + 1) as u64;
        let pivot = a[r].clone();
        for i in 0..a.len() {
            if i != r && a[i][c] != 0 {
                let fct = a[i][c];
                for (x, &y) in a[i].iter_mut().zip(&pivot) {
                    *x = (*x + ell - mulm(fct, y, ell)) % ell;
                }
                ops += (cols + 1) as u64;
            }
        }
        piv_col.push(c);
        r += 1;
    }
    // k is determined iff column m is a pivot whose row has no other free
    // (non-pivot) column with a non-zero entry.
    let pivots: BTreeSet<usize> = piv_col.iter().copied().collect();
    for (row, &c) in a.iter().zip(&piv_col) {
        if c == m {
            let free_dep = (0..cols).any(|j| j != m && !pivots.contains(&j) && row[j] != 0);
            return (if free_dep { None } else { Some(row[cols]) }, ops);
        }
    }
    (None, ops)
}

// ── Pollard rho on the same group, for the paired reference ──────────────

/// Pollard rho with an `r`-adding walk (`r = 20`) and Brent's cycle
/// detection on `⟨P⟩` of prime order `ℓ`, without the negation map.
/// Returns `(k, group operations)` with `kP = Q`; all group operations are
/// counted, the `r` multipliers' set-up included.
pub fn rho(curve: &Curve5, pp: &Pt5, qq: &Pt5, ell: u128, rng: &mut StdRng) -> (Option<u128>, u64) {
    const R: usize = 20;
    let ops0 = curve.ops();
    loop {
        let mults: Vec<(u128, u128, Pt5)> = (0..R)
            .map(|_| {
                let a = rng.gen_range(0..ell);
                let b = rng.gen_range(0..ell);
                (a, b, curve.add(&curve.mul(pp, a), &curve.mul(qq, b)))
            })
            .collect();
        let idx = |pt: &Pt5| -> usize {
            if pt.inf {
                0
            } else {
                ((pt.x.0[0] * 7 + pt.x.0[1] * 3 + pt.x.0[2]) % R as u64) as usize
            }
        };
        let (mut a, mut b) = (rng.gen_range(0..ell), rng.gen_range(0..ell));
        let mut x = curve.add(&curve.mul(pp, a), &curve.mul(qq, b));
        let (mut sa, mut sb, mut sx) = (a, b, x);
        let mut power = 1u64;
        let mut lam = 0u64;
        loop {
            let j = idx(&x);
            x = curve.add(&x, &mults[j].2);
            a = (a + mults[j].0) % ell;
            b = (b + mults[j].1) % ell;
            lam += 1;
            if x == sx {
                // aP + bQ = saP + sbQ  ⇒  (a − sa) = k (sb − b)
                let db = (sb + ell - b) % ell;
                if db == 0 {
                    break;
                }
                let k = mulm((a + ell - sa) % ell, inv_mod_u128(db, ell).unwrap(), ell);
                return (Some(k), curve.ops() - ops0);
            }
            if lam == power {
                (sa, sb, sx) = (a, b, x);
                power *= 2;
                lam = 0;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn fp5_inverse_frobenius_and_sqrt() {
        let f = Fp5::new(11);
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..50 {
            let a = f.random(&mut rng);
            if a.is_zero() {
                continue;
            }
            assert_eq!(f.mul(&a, &f.inv(&a)), E5::ONE);
            assert_eq!(f.frobenius(&a), f.pow(&a, 11));
            let s = f.sq(&a);
            let r = f.sqrt(&s).unwrap();
            assert_eq!(f.sq(&r), s);
            assert_eq!(f.is_square(&a), f.pow(&a, (f.order() / 2) as u128) == E5::ONE);
        }
    }

    #[test]
    fn s5_vanishes_on_four_point_sums() {
        let mut rng = StdRng::seed_from_u64(2);
        let c = Curve5::random(11, &mut rng);
        let f = &c.f;
        let pts: Vec<Pt5> = std::iter::repeat_with(|| c.lift_x(&f.random(&mut rng)))
            .flatten()
            .take(4)
            .collect();
        let s = c.add(&c.add(&pts[0], &pts[1]), &c.add(&pts[2], &pts[3]));
        let v = c.s5(&[pts[0].x, pts[1].x, pts[2].x, pts[3].x, s.x]);
        assert!(v.is_zero());
        // and S₃ on two-point sums
        let t = c.add(&pts[0], &pts[1]);
        assert!(c.s3(&pts[0].x, &pts[1].x, &t.x).is_zero());
    }

    #[test]
    fn symmetrised_s5_and_descended_system_find_planted_decomposition() {
        let mut rng = StdRng::seed_from_u64(3);
        let c = Curve5::random(11, &mut rng);
        let sym = SymmetrisedS5x5::precompute(&c, &mut rng);
        let base = factor_base(&c);
        assert!(base.len() >= 4);
        let index: HashMap<u64, usize> = base.iter().enumerate().map(|(i, p)| (p.x.0[0], i)).collect();
        let r = c.add(&c.add(&base[0], &base[1]), &c.sub(&base[2], &base[3]));
        let rep = trial(&c, &sym, &base, &index, &r, &trial_options(), None);
        assert!(rep.consistent);
        assert!(!rep.relations.is_empty());
        for rel in &rep.relations {
            assert!(verify_relation(&c, &base, rel, &r));
        }
        let pairs = PairTable::new(&c, &base);
        assert_eq!(mitm_decompositions(&c, &base, &pairs, &r), rep.decompositions);
    }
}
