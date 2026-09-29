//! **Joux–Vitse `(k − 1)`-point decompositions at `k = 5`: the stage the
//! rho-parity programme registered.**
//!
//! Companion to `research/notes/index-calculus/RESEARCH_RHO_PARITY_PROGRAMME.md`
//! §6.  On `E(F_{p⁵})` with the subspace base `F = {P : x(P) ∈ F_p}`, a
//! residual is a sum of *four* base points about `1/(24p)` of the time, so
//! the relation phase tries `≈ 12p²` residuals — `∝ n^{2/5}` against rho's
//! `n^{1/2}` — and the linear algebra over `|F| ≈ p/2` unknowns is `∝ n^{2/5}`
//! too.  The method's `S / rho` therefore *closes*, as `n^{−1/10}`, with the
//! constant set by `C″`, the cost of one four-point test: the symmetrised
//! `S₅(x₁, …, x₄, x_R)` Weil-restricted to five equations of degree `≤ 8` in
//! four unknowns, overdetermined by one.  This module builds the field, the
//! curves, the test and its independent oracle, and measures `C″`.
//!
//! `F_{p⁵} = F_p[t]/(t⁵ − c)` needs `p ≡ 1 (mod 5)` (then a non-fifth-power
//! `c` exists and the binomial is irreducible); the Frobenius `t ↦ ζt` with
//! `ζ = c^{(p−1)/5}` gives the inverse through the norm.

use std::cell::Cell;
use std::collections::{BTreeSet, HashMap};
use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::Serialize;

use super::f4_fp::{self, F4Options, Ordering as F4Ordering, Verdict};
use super::gaudry_cubic::{am, mm, sm, PolyRing, UPoly};
use super::gaudry_quartic::exponents4;
use super::residual_walk::{inv_mod, is_prime_u64, pow_mod};

// ── F_{p⁵} = F_p[t] / (t⁵ − c) ───────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default, PartialOrd, Ord)]
pub struct E5(pub [u64; 5]);

impl E5 {
    pub const ZERO: E5 = E5([0; 5]);
    pub const ONE: E5 = E5([1, 0, 0, 0, 0]);
    pub fn is_zero(&self) -> bool {
        self.0 == [0; 5]
    }
    pub fn in_fp(&self) -> bool {
        self.0[1..].iter().all(|&v| v == 0)
    }
}

/// `F_{p⁵}` with a multiplication counter in `F_p` multiplications, by the
/// convention of [`super::gaudry_quartic::Fp4`]: a product is `25`
/// multiplications and `4` reductions by `c`, `29`; an inverse is charged
/// as performed (four Frobenius twists, three products, the norm's constant
/// term, one `F_p` inversion at `16`, one scaling), `133`.
#[derive(Clone)]
pub struct Fp5 {
    pub p: u64,
    pub c: u64,
    /// `c^{(p−1)/5}`, a primitive fifth root of unity: `t^p = zeta·t`.
    pub zeta: u64,
    muls: Cell<u64>,
}

impl Fp5 {
    pub fn new(p: u64) -> Fp5 {
        assert!(
            p % 5 == 1,
            "t⁵ − c is irreducible over F_p only when p ≡ 1 (mod 5)"
        );
        let c = (2..p)
            .find(|&c| pow_mod(c, (p - 1) / 5, p) != 1)
            .expect("a non-fifth-power");
        let zeta = pow_mod(c, (p - 1) / 5, p);
        Fp5 {
            p,
            c,
            zeta,
            muls: Cell::new(0),
        }
    }
    fn count(&self, k: u64) {
        self.muls.set(self.muls.get() + k);
    }
    pub fn muls(&self) -> u64 {
        self.muls.get()
    }
    pub fn reset_muls(&self) {
        self.muls.set(0);
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
    /// Multiplication by an element of `F_p`: five multiplications.
    pub fn scale(&self, a: &E5, k: u64) -> E5 {
        self.count(5);
        E5(core::array::from_fn(|i| mm(a.0[i], k, self.p)))
    }
    /// Schoolbook with `t⁵ = c`: 25 products and 4 reductions.
    pub fn mul(&self, a: &E5, b: &E5) -> E5 {
        let p = self.p;
        self.count(29);
        let mut d = [0u64; 9];
        for i in 0..5 {
            for j in 0..5 {
                d[i + j] = am(d[i + j], mm(a.0[i], b.0[j], p), p);
            }
        }
        E5([
            am(d[0], mm(self.c, d[5], p), p),
            am(d[1], mm(self.c, d[6], p), p),
            am(d[2], mm(self.c, d[7], p), p),
            am(d[3], mm(self.c, d[8], p), p),
            d[4],
        ])
    }
    pub fn sq(&self, a: &E5) -> E5 {
        self.mul(a, a)
    }
    /// The Frobenius `a ↦ a^p`: `a_i t^i ↦ a_i ζ^i t^i`, four multiplications.
    pub fn frobenius(&self, a: &E5) -> E5 {
        let p = self.p;
        self.count(4);
        let mut z = 1u64;
        let mut out = [0u64; 5];
        for i in 0..5 {
            out[i] = mm(a.0[i], z, p);
            z = mm(z, self.zeta, p);
        }
        E5(out)
    }
    /// Inverse through the norm: `a⁻¹ = B / N(a)` with `B = a^p a^{p²} a^{p³}
    /// a^{p⁴}` and `N(a) = a·B ∈ F_p`.
    pub fn inv(&self, a: &E5) -> E5 {
        let p = self.p;
        assert!(!a.is_zero(), "inverse of zero in F_p⁵");
        let f1 = self.frobenius(a);
        let f2 = self.frobenius(&f1);
        let f3 = self.frobenius(&f2);
        let f4 = self.frobenius(&f3);
        let b = self.mul(&self.mul(&f1, &f2), &self.mul(&f3, &f4));
        // Constant term of a·B (the norm; the other coordinates vanish).
        let mut n0 = mm(a.0[0], b.0[0], p);
        for i in 1..5 {
            n0 = am(n0, mm(self.c, mm(a.0[i], b.0[5 - i], p), p), p);
        }
        self.count(9 + 16);
        let ninv = inv_mod(n0, p);
        self.scale(&b, ninv)
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
    pub fn is_square(&self, a: &E5) -> bool {
        a.is_zero() || self.pow(a, self.order() / 2) == E5::ONE
    }
    pub fn random(&self, rng: &mut StdRng) -> E5 {
        E5(core::array::from_fn(|_| rng.gen_range(0..self.p)))
    }
    /// Tonelli–Shanks over `F_{p⁵}^×`.
    pub fn sqrt(&self, a: &E5, rng: &mut StdRng) -> Option<E5> {
        if a.is_zero() {
            return Some(E5::ZERO);
        }
        if !self.is_square(a) {
            return None;
        }
        let q1 = self.order();
        let s = q1.trailing_zeros();
        let odd = q1 >> s;
        let z = loop {
            let z = self.random(rng);
            if !z.is_zero() && !self.is_square(&z) {
                break z;
            }
        };
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
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

pub struct Curve5 {
    pub f: Fp5,
    pub a: E5,
    pub b: E5,
    ops: Cell<u64>,
}

impl Curve5 {
    /// Coefficients outside `F_p` (the only proper subfield of `F_{p⁵}`),
    /// so the curve is not a twist of one defined over a subfield.
    pub fn random(p: u64, rng: &mut StdRng) -> Curve5 {
        let f = Fp5::new(p);
        loop {
            let a = f.random(rng);
            let b = f.random(rng);
            if a.in_fp() || b.in_fp() {
                continue;
            }
            let a3 = f.mul(&f.sq(&a), &a);
            let disc = f.add(&f.scale(&a3, 4), &f.scale(&f.sq(&b), 27));
            if !disc.is_zero() {
                return Curve5 {
                    f,
                    a,
                    b,
                    ops: Cell::new(0),
                };
            }
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
    pub fn rhs(&self, x: &E5) -> E5 {
        let f = &self.f;
        f.add(&f.mul(&f.add(&f.sq(x), &self.a), x), &self.b)
    }
    pub fn lift_x(&self, x: &E5, rng: &mut StdRng) -> Option<Pt5> {
        let y = self.f.sqrt(&self.rhs(x), rng)?;
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
    pub fn ops(&self) -> u64 {
        self.ops.get()
    }
    pub fn reset_ops(&self) {
        self.ops.set(0);
    }
    /// Affine addition: one inversion and three multiplications for
    /// distinct points.
    pub fn add(&self, p1: &Pt5, p2: &Pt5) -> Pt5 {
        let f = &self.f;
        self.ops.set(self.ops.get() + 1);
        if p1.inf {
            return *p2;
        }
        if p2.inf {
            return *p1;
        }
        if p1.x == p2.x {
            if p1.y != p2.y || p1.y.is_zero() {
                return Pt5::INF;
            }
            let num = f.add(&f.scale(&f.sq(&p1.x), 3), &self.a);
            let lam = f.mul(&num, &f.inv(&f.scale(&p1.y, 2)));
            let x3 = f.sub(&f.sq(&lam), &f.add(&p1.x, &p1.x));
            let y3 = f.sub(&f.mul(&lam, &f.sub(&p1.x, &x3)), &p1.y);
            return Pt5 {
                x: x3,
                y: y3,
                inf: false,
            };
        }
        let lam = f.mul(&f.sub(&p2.y, &p1.y), &f.inv(&f.sub(&p2.x, &p1.x)));
        let x3 = f.sub(&f.sub(&f.sq(&lam), &p1.x), &p2.x);
        let y3 = f.sub(&f.mul(&lam, &f.sub(&p1.x, &x3)), &p1.y);
        Pt5 {
            x: x3,
            y: y3,
            inf: false,
        }
    }
    pub fn sub(&self, p1: &Pt5, p2: &Pt5) -> Pt5 {
        self.add(p1, &self.neg(p2))
    }
    pub fn mul(&self, pt: &Pt5, mut k: u64) -> Pt5 {
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

fn det5(f: &Fp5, mut m: Vec<Vec<E5>>) -> E5 {
    let n = m.len();
    let mut det = E5::ONE;
    for c in 0..n {
        let Some(r) = (c..n).find(|&r| !m[r][c].is_zero()) else {
            return E5::ZERO;
        };
        if r != c {
            m.swap(r, c);
            det = f.neg(&det);
        }
        det = f.mul(&det, &m[c][c]);
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
    det
}

impl Curve5 {
    /// `S₃(x₁, x₂, X)` as a polynomial in `X`.
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
    /// `S₄(x₃, x₄, x₅, Y)` as a polynomial in `Y`, of degree ≤ 4.
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
    /// `S₅(x₁, …, x₅) = Res_Y(S₃(x₁, x₂, Y), S₄(x₃, x₄, x₅, Y))`.
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
        det5(f, m)
    }
}

// ── The symmetrised S₅ in four points, by interpolation ───────────────────

fn elementary4(f: &Fp5, x: &[E5; 4]) -> [E5; 4] {
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

fn mono_values4(f: &Fp5, e: &[E5; 4], monos: &[[u8; 4]]) -> Vec<E5> {
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

/// `S₅(x₁, …, x₄, X) = H(e₁, …, e₄, X)`, `H` of total degree `≤ 8` in the
/// `e`s and degree `≤ 8` in `X`: `495 · 9` coefficients in `F_{p⁵}`, found by
/// interpolation as [`super::gaudry_quartic::SymmetrisedS5`] does over
/// `F_{p⁴}`, and checked against fresh evaluations.
pub struct SymmetrisedS5Q {
    pub monos: Vec<[u8; 4]>,
    pub coef: Vec<[E5; 9]>,
    pub precompute_muls: u64,
}

impl SymmetrisedS5Q {
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

    pub fn precompute(curve: &Curve5, rng: &mut StdRng) -> SymmetrisedS5Q {
        let f = &curve.f;
        let p = f.p;
        let before = f.muls();
        let vinv: [[u64; 9]; 9] = {
            let mut m: Vec<Vec<u64>> = (0..9u64)
                .map(|k| {
                    let mut row: Vec<u64> = (0..9).map(|j| pow_mod(k, j, p)).collect();
                    row.extend((0..9).map(|i| u64::from(i == k)));
                    row
                })
                .collect();
            let mut scratch = 0u64;
            super::gaudry_cubic::rref_mod_p(&mut m, p, &mut scratch);
            core::array::from_fn(|j| core::array::from_fn(|k| m[j][9 + k]))
        };
        let monos = exponents4(8);
        let n = monos.len();
        loop {
            let mut sys: Vec<Vec<E5>> = Vec::with_capacity(n);
            for _ in 0..n {
                let x: [E5; 4] = core::array::from_fn(|_| f.random(rng));
                let e = elementary4(f, &x);
                let mut row = mono_values4(f, &e, &monos);
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
                sys.par_iter_mut().enumerate().for_each(|(r, row)| {
                    if r == c || row[c].is_zero() {
                        return;
                    }
                    let fl = Fp5::new(p);
                    let k = row[c];
                    for j in c..row.len() {
                        if !pivot[j].is_zero() {
                            let v = fl.mul(&k, &pivot[j]);
                            row[j] = fl.sub(&row[j], &v);
                        }
                    }
                });
                f.count(29 * (n as u64) * (n + 9 - c) as u64);
            }
            if !ok {
                continue;
            }
            let coef: Vec<[E5; 9]> = (0..n)
                .map(|r| core::array::from_fn(|j| sys[r][n + j]))
                .collect();
            let pre = SymmetrisedS5Q {
                monos: monos.clone(),
                coef,
                precompute_muls: 0,
            };
            let good = (0..6).all(|_| {
                let x: [E5; 4] = core::array::from_fn(|_| f.random(rng));
                let want = Self::coefficients_in_last(curve, &x, &vinv);
                let e = elementary4(f, &x);
                let mv = mono_values4(f, &e, &pre.monos);
                (0..9).all(|j| {
                    let got = mv
                        .iter()
                        .zip(&pre.coef)
                        .fold(E5::ZERO, |acc, (m, c)| f.add(&acc, &f.mul(m, &c[j])));
                    got == want[j]
                })
            });
            assert!(good, "the symmetrised S₅ does not reproduce the resultant");
            return SymmetrisedS5Q {
                precompute_muls: f.muls() - before,
                ..pre
            };
        }
    }

    pub fn eval(&self, f: &Fp5, e: &[E5; 4], x_r: &E5) -> E5 {
        let mut pows = [E5::ONE; 9];
        for k in 1..9 {
            pows[k] = f.mul(&pows[k - 1], x_r);
        }
        let mv = mono_values4(f, e, &self.monos);
        mv.iter().zip(&self.coef).fold(E5::ZERO, |acc, (m, c)| {
            let v = (0..9).fold(E5::ZERO, |a, j| f.add(&a, &f.mul(&c[j], &pows[j])));
            f.add(&acc, &f.mul(m, &v))
        })
    }

    /// The five `F_p`-components of `H(e₁, …, e₄, x_R)`: five polynomials of
    /// total degree `≤ 8` in four `F_p` unknowns.
    pub fn weil_restrict(&self, f: &Fp5, x_r: &E5) -> [HashMap<[u8; 4], u64>; 5] {
        let mut pows = [E5::ONE; 9];
        for k in 1..9 {
            pows[k] = f.mul(&pows[k - 1], x_r);
        }
        let mut out: [HashMap<[u8; 4], u64>; 5] = Default::default();
        for (m, c) in self.monos.iter().zip(&self.coef) {
            let v = (0..9).fold(E5::ZERO, |acc, j| f.add(&acc, &f.mul(&c[j], &pows[j])));
            for (i, comp) in out.iter_mut().enumerate() {
                if v.0[i] != 0 {
                    comp.insert(*m, v.0[i]);
                }
            }
        }
        out
    }
}

// ── Instances, the base, the oracle ──────────────────────────────────────

pub struct Instance5 {
    pub curve: Curve5,
    pub n: u64,
    pub g: Pt5,
    pub d: u64,
    pub q: Pt5,
}

fn random_point(curve: &Curve5, rng: &mut StdRng) -> Pt5 {
    loop {
        let x = curve.f.random(rng);
        if let Some(pt) = curve.lift_x(&x, rng) {
            return pt;
        }
    }
}

fn isqrt(n: u64) -> u64 {
    let mut r = (n as f64).sqrt() as u64;
    while r * r > n {
        r -= 1;
    }
    while (r + 1) * (r + 1) <= n {
        r += 1;
    }
    r
}

/// `#E(F_{p⁵})` when it is prime: baby-step giant-step for the multiple of a
/// random point's order inside the Hasse interval `q + 1 ± 2√q`.
fn prime_group_order(curve: &Curve5, rng: &mut StdRng) -> Option<u64> {
    let p = curve.f.p;
    let q = p.pow(5);
    let two_sqrt = 2 * isqrt(q) + 2;
    let lo = q + 1 - two_sqrt;
    let width = 2 * two_sqrt;
    let pt = random_point(curve, rng);
    let steps = isqrt(width) + 1;
    let mut table: HashMap<(E5, E5), u64> = HashMap::new();
    let mut jp = Pt5::INF;
    for j in 0..steps {
        if !jp.inf {
            table.entry((jp.x, jp.y)).or_insert(j);
        }
        jp = curve.add(&jp, &pt);
    }
    let giant = curve.mul(&pt, steps);
    let mut t = curve.mul(&pt, lo);
    let mut i = 0u64;
    while i * steps <= width + steps {
        let m = if t.inf {
            Some(lo + i * steps)
        } else {
            let neg = curve.neg(&t);
            table.get(&(neg.x, neg.y)).map(|&j| lo + i * steps + j)
        };
        if let Some(m) = m {
            return is_prime_u64(m).then_some(m);
        }
        t = curve.add(&t, &giant);
        i += 1;
    }
    None
}

pub fn generate_instance5(p: u64, seed: u64) -> Instance5 {
    let mut rng = StdRng::seed_from_u64(seed ^ 0x5A11);
    loop {
        let curve = Curve5::random(p, &mut rng);
        let Some(n) = prime_group_order(&curve, &mut rng) else {
            continue;
        };
        let g = random_point(&curve, &mut rng);
        if !curve.mul(&g, n).inf {
            continue;
        }
        let d = rng.gen_range(1..n);
        let q = curve.mul(&g, d);
        curve.reset_ops();
        curve.f.reset_muls();
        return Instance5 { curve, n, g, d, q };
    }
}

pub fn factor_base(curve: &Curve5, rng: &mut StdRng) -> Vec<Pt5> {
    (0..curve.f.p)
        .filter_map(|x| curve.lift_x(&curve.f.from_fp(x), rng))
        .collect()
}

pub struct PairTable5 {
    map: HashMap<E5, Vec<(u32, u32, i8, E5)>>,
}

pub fn pair_table(curve: &Curve5, base: &[Pt5]) -> PairTable5 {
    let mut map: HashMap<E5, Vec<(u32, u32, i8, E5)>> = HashMap::new();
    for i in 0..base.len() {
        for j in (i + 1)..base.len() {
            for (s, v) in [
                (1i8, curve.add(&base[i], &base[j])),
                (-1, curve.sub(&base[i], &base[j])),
            ] {
                if !v.inf {
                    map.entry(v.x)
                        .or_default()
                        .push((i as u32, j as u32, s, v.y));
                }
            }
        }
    }
    PairTable5 { map }
}

/// One decomposition `R = Σ s_t P_{i_t}` over four distinct base points.
pub type Quad = [(usize, i64); 4];

/// Every `R = s_k P_k + s_l P_l + t (P_i + s P_j)` by meet in the middle
/// against the pair table: `4·C(|F|, 2)` group operations.
pub fn mitm4_signed(curve: &Curve5, base: &[Pt5], table: &PairTable5, r: &Pt5) -> Vec<Quad> {
    let n = base.len();
    let mut out: BTreeSet<Quad> = BTreeSet::new();
    for k in 0..n {
        for sk in [1i64, -1] {
            let pk = if sk == 1 {
                base[k]
            } else {
                curve.neg(&base[k])
            };
            let tk = curve.sub(r, &pk);
            for l in (k + 1)..n {
                for sl in [1i64, -1] {
                    let pl = if sl == 1 {
                        base[l]
                    } else {
                        curve.neg(&base[l])
                    };
                    let v = curve.sub(&tk, &pl);
                    if v.inf {
                        continue;
                    }
                    let Some(cands) = table.map.get(&v.x) else {
                        continue;
                    };
                    for &(i, j, s, y) in cands {
                        let (i, j) = (i as usize, j as usize);
                        if i == k || i == l || j == k || j == l {
                            continue;
                        }
                        let t: i64 = if v.y == y { 1 } else { -1 };
                        let mut terms = [(k, sk), (l, sl), (i, t), (j, t * s as i64)];
                        terms.sort_unstable();
                        out.insert(terms);
                    }
                }
            }
        }
    }
    out.into_iter().collect()
}

fn verify_quad(curve: &Curve5, base: &[Pt5], r: &Pt5, t: &Quad) -> bool {
    let mut acc = *r;
    for &(i, s) in t {
        acc = if s == 1 {
            curve.sub(&acc, &base[i])
        } else {
            curve.add(&acc, &base[i])
        };
    }
    acc.inf
}

// ── The four-point test ───────────────────────────────────────────────────

#[derive(Clone, Debug, Default, Serialize)]
pub struct Jv4Cost {
    pub weil_muls: u64,
    pub f4_muls: u64,
    pub roots_muls: u64,
    pub group_ops: u64,
    pub solving_degree: u32,
    pub degree_reached: u32,
    pub max_rows: usize,
    pub max_cols: usize,
    pub f4_ms: f64,
    pub inconsistent: bool,
    pub undetermined: bool,
    pub timed_out: bool,
    pub e_solutions: usize,
    pub decompositions: usize,
}

fn polys_of(comps: &[HashMap<[u8; 4], u64>; 5], p: u64) -> Vec<f4_fp::Poly> {
    comps
        .iter()
        .map(|c| {
            let terms: Vec<(Vec<u32>, u64)> = c
                .iter()
                .map(|(e, v)| (e.iter().map(|&x| x as u32).collect(), *v))
                .collect();
            f4_fp::normalise(&terms, p, F4Ordering::Grevlex)
        })
        .filter(|f| !f.is_empty())
        .collect()
}

/// Grevlex F4 with pairs bounded at `max_degree` and a wall-clock budget.
pub fn jv4_options(max_degree: u32, budget_secs: f64) -> F4Options {
    F4Options::new(F4Ordering::Grevlex, max_degree)
        .with_budget(std::time::Duration::from_secs_f64(budget_secs))
}

fn resolve_signs(curve: &Curve5, base: &[Pt5], r: &Pt5, idx: [usize; 4]) -> Option<Quad> {
    let [i1, i2, i3, i4] = idx;
    for s2 in [1i64, -1] {
        for s3 in [1i64, -1] {
            for s4 in [1i64, -1] {
                let pick = |i: usize, s: i64| if s == 1 { base[i] } else { curve.neg(&base[i]) };
                let t = curve.sub(
                    &curve.sub(&curve.sub(r, &pick(i2, s2)), &pick(i3, s3)),
                    &pick(i4, s4),
                );
                if t.inf || t.x != base[i1].x {
                    continue;
                }
                let s1 = if t.y == base[i1].y { 1 } else { -1 };
                let mut terms = [(i1, s1), (i2, s2), (i3, s3), (i4, s4)];
                terms.sort_unstable();
                return Some(terms);
            }
        }
    }
    None
}

/// Every four-point decomposition of `r` over the base, by the symmetrised
/// `S₅` and F4, with the cost split.
pub fn jv4_decompose(
    curve: &Curve5,
    base: &[Pt5],
    by_x: &HashMap<u64, usize>,
    pre: &SymmetrisedS5Q,
    r: &Pt5,
    opts: &F4Options,
    rng: &mut StdRng,
) -> (Vec<Quad>, Jv4Cost) {
    let f = &curve.f;
    let p = f.p;
    let mut cost = Jv4Cost::default();
    let mut out: BTreeSet<Quad> = BTreeSet::new();
    if r.inf {
        return (Vec::new(), cost);
    }
    let m0 = f.muls();
    let comps = pre.weil_restrict(f, &r.x);
    cost.weil_muls = f.muls() - m0;
    let polys = polys_of(&comps, p);
    let rep = f4_fp::solve(&polys, 4, p, opts);
    cost.f4_muls = rep.field_ops;
    cost.solving_degree = rep.solving_degree;
    cost.degree_reached = rep.degree_reached;
    cost.max_rows = rep.max_rows;
    cost.max_cols = rep.max_cols;
    cost.f4_ms = rep.ms;
    cost.timed_out = rep.timed_out;
    let es = match rep.verdict {
        Verdict::Inconsistent => {
            cost.inconsistent = true;
            Vec::new()
        }
        Verdict::Undetermined => {
            cost.undetermined = true;
            Vec::new()
        }
        Verdict::Solutions(s) => s,
    };
    cost.e_solutions = es.len();
    let ring = PolyRing::new(p);
    let g0 = curve.ops();
    for e in &es {
        // T⁴ − e₁T³ + e₂T² − e₃T + e₄, low to high.
        let quartic = UPoly(vec![e[3], (p - e[2]) % p, e[1], (p - e[0]) % p, 1]);
        let roots = ring.roots(&quartic, rng);
        if roots.len() != 4 {
            continue;
        }
        let Some(idx) = roots
            .iter()
            .map(|x| by_x.get(x).copied())
            .collect::<Option<Vec<usize>>>()
        else {
            continue;
        };
        let mut sorted = idx.clone();
        sorted.sort_unstable();
        sorted.dedup();
        if sorted.len() != 4 {
            continue;
        }
        if let Some(t) = resolve_signs(curve, base, r, [idx[0], idx[1], idx[2], idx[3]]) {
            out.insert(t);
        }
    }
    cost.roots_muls = ring.muls.get();
    cost.group_ops = curve.ops() - g0;
    cost.decompositions = out.len();
    (out.into_iter().collect(), cost)
}

// ── Experiment: C″ ────────────────────────────────────────────────────────

#[derive(Clone, Debug, Default, Serialize)]
pub struct CostStats5 {
    pub mean: f64,
    pub min: u64,
    pub max: u64,
}

fn stats(v: &[u64]) -> CostStats5 {
    if v.is_empty() {
        return CostStats5::default();
    }
    CostStats5 {
        mean: v.iter().sum::<u64>() as f64 / v.len() as f64,
        min: *v.iter().min().unwrap(),
        max: *v.iter().max().unwrap(),
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct CsecondReport {
    pub p: u64,
    pub seed: u64,
    pub n: u64,
    pub bits: f64,
    pub base: usize,
    pub fp_muls_per_add: f64,
    pub precompute_muls: u64,
    pub precompute_ms: f64,
    pub max_degree: u32,
    pub budget_secs: f64,
    pub random_residuals: usize,
    pub constructed_residuals: usize,
    pub random_decomposable: usize,
    pub expected_rate: f64,
    pub planted_found: usize,
    /// Decompositions the test returned that fail the group check (none expected).
    pub unverified: usize,
    pub mismatches: usize,
    pub undetermined: usize,
    pub timed_out: usize,
    /// Per-test cost over the random residuals that finished: `C″` and its split.
    pub c_second: CostStats5,
    pub weil_muls: CostStats5,
    pub f4_muls: CostStats5,
    pub roots_muls: CostStats5,
    pub sign_group_ops: CostStats5,
    pub f4_ms: CostStats5,
    pub solving_degree: CostStats5,
    pub degree_reached: CostStats5,
    pub max_rows: CostStats5,
    pub max_cols: CostStats5,
    /// The same over the constructed residuals that finished.
    pub c_second_decomposable: CostStats5,
    pub f4_ms_decomposable: CostStats5,
    pub degree_reached_decomposable: CostStats5,
    pub max_cols_decomposable: CostStats5,
    pub per_test: Vec<Jv4Cost>,
    pub mitm4_group_ops: f64,
    pub wall_ms: f64,
}

/// `C″` measured: `random` residuals `aG + bQ` and `constructed` sums of four
/// base points, every one cross-checked against [`mitm4_signed`].
pub fn run_jv5_csecond(
    p: u64,
    seed: u64,
    random: usize,
    constructed: usize,
    max_degree: u32,
    budget_secs: f64,
) -> CsecondReport {
    let start = Instant::now();
    let inst = generate_instance5(p, seed);
    let curve = &inst.curve;
    let f = &curve.f;
    let n = inst.n;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x5C02);
    let base = factor_base(curve, &mut rng);
    let by_x: HashMap<u64, usize> = base
        .iter()
        .enumerate()
        .map(|(i, q)| (q.x.0[0], i))
        .collect();
    let fp_per_add = {
        f.reset_muls();
        let mut acc = base[0];
        for _ in 0..64 {
            acc = curve.add(&acc, &base[1]);
        }
        f.muls() as f64 / 64.0
    };
    let t0 = Instant::now();
    let pre = SymmetrisedS5Q::precompute(curve, &mut rng);
    let precompute_ms = t0.elapsed().as_secs_f64() * 1e3;
    let table = pair_table(curve, &base);
    let mut rep = CsecondReport {
        p,
        seed,
        n,
        bits: (n as f64).log2(),
        base: base.len(),
        fp_muls_per_add: fp_per_add,
        precompute_muls: pre.precompute_muls,
        precompute_ms,
        max_degree,
        budget_secs,
        random_residuals: random,
        constructed_residuals: constructed,
        expected_rate: 1.0 / (24.0 * p as f64),
        ..Default::default()
    };
    let mut costs: Vec<Jv4Cost> = Vec::new();
    let mut costs_dec: Vec<Jv4Cost> = Vec::new();
    let mut mitm_ops = 0u64;
    let mut check_rng = StdRng::seed_from_u64(seed ^ 0x5EED);
    let mut check = |r: &Pt5, want_planted: Option<Quad>, rep: &mut CsecondReport| -> Jv4Cost {
        let opts = jv4_options(max_degree, budget_secs);
        let (found, cost) = jv4_decompose(curve, &base, &by_x, &pre, r, &opts, &mut check_rng);
        rep.unverified += found
            .iter()
            .filter(|t| !verify_quad(curve, &base, r, t))
            .count();
        let g0 = curve.ops();
        let oracle = mitm4_signed(curve, &base, &table, r);
        mitm_ops += curve.ops() - g0;
        if cost.undetermined {
            rep.undetermined += 1;
        } else if found != oracle {
            rep.mismatches += 1;
        }
        if cost.timed_out {
            rep.timed_out += 1;
        }
        if let Some(t) = want_planted {
            if found.contains(&t) {
                rep.planted_found += 1;
            }
        }
        cost
    };
    for _ in 0..random {
        let a = rng.gen_range(0..n);
        let b = rng.gen_range(1..n);
        let r = curve.add(&curve.mul(&inst.g, a), &curve.mul(&inst.q, b));
        let c = check(&r, None, &mut rep);
        if c.decompositions > 0 {
            rep.random_decomposable += 1;
        }
        costs.push(c);
    }
    for _ in 0..constructed {
        let idx: [usize; 4] = loop {
            let idx: [usize; 4] = core::array::from_fn(|_| rng.gen_range(0..base.len()));
            let mut s = idx.to_vec();
            s.sort_unstable();
            s.dedup();
            if s.len() == 4 {
                break idx;
            }
        };
        let signs: [i64; 4] = core::array::from_fn(|_| if rng.gen_bool(0.5) { 1 } else { -1 });
        let mut r = Pt5::INF;
        for t in 0..4 {
            let q = if signs[t] == 1 {
                base[idx[t]]
            } else {
                curve.neg(&base[idx[t]])
            };
            r = curve.add(&r, &q);
        }
        let mut terms = [
            (idx[0], signs[0]),
            (idx[1], signs[1]),
            (idx[2], signs[2]),
            (idx[3], signs[3]),
        ];
        terms.sort_unstable();
        let c = check(&r, Some(terms), &mut rep);
        costs_dec.push(c);
    }
    let total = |c: &Jv4Cost| {
        c.weil_muls + c.f4_muls + c.roots_muls + (c.group_ops as f64 * fp_per_add) as u64
    };
    let done =
        |v: &[Jv4Cost]| -> Vec<Jv4Cost> { v.iter().filter(|c| !c.undetermined).cloned().collect() };
    let col = |v: &[Jv4Cost], g: &dyn Fn(&Jv4Cost) -> u64| -> CostStats5 {
        stats(&v.iter().map(g).collect::<Vec<u64>>())
    };
    let cd = done(&costs);
    let dd = done(&costs_dec);
    rep.c_second = col(&cd, &total);
    rep.weil_muls = col(&cd, &|c| c.weil_muls);
    rep.f4_muls = col(&cd, &|c| c.f4_muls);
    rep.roots_muls = col(&cd, &|c| c.roots_muls);
    rep.sign_group_ops = col(&cd, &|c| c.group_ops);
    rep.f4_ms = col(&cd, &|c| c.f4_ms as u64);
    rep.solving_degree = col(&cd, &|c| c.solving_degree as u64);
    rep.degree_reached = col(&cd, &|c| c.degree_reached as u64);
    rep.max_rows = col(&cd, &|c| c.max_rows as u64);
    rep.max_cols = col(&cd, &|c| c.max_cols as u64);
    rep.c_second_decomposable = col(&dd, &total);
    rep.f4_ms_decomposable = col(&dd, &|c| c.f4_ms as u64);
    rep.degree_reached_decomposable = col(&dd, &|c| c.degree_reached as u64);
    rep.max_cols_decomposable = col(&dd, &|c| c.max_cols as u64);
    rep.per_test = costs.iter().chain(costs_dec.iter()).cloned().collect();
    rep.mitm4_group_ops = mitm_ops as f64 / (random + constructed).max(1) as f64;
    rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
    rep
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fp5_is_a_field_with_a_norm_inverse() {
        let f = Fp5::new(31);
        let mut rng = StdRng::seed_from_u64(1);
        assert_eq!(pow_mod(f.zeta, 5, 31), 1);
        assert_ne!(f.zeta, 1);
        for _ in 0..50 {
            let a = f.random(&mut rng);
            let b = f.random(&mut rng);
            assert_eq!(f.mul(&a, &b), f.mul(&b, &a));
            // Frobenius is a^p.
            assert_eq!(f.frobenius(&a), f.pow(&a, 31));
            if !a.is_zero() {
                assert_eq!(f.mul(&a, &f.inv(&a)), E5::ONE);
            }
            let s = f.sq(&a);
            let r = f.sqrt(&s, &mut rng).expect("a square has a root");
            assert_eq!(f.sq(&r), s);
        }
        f.reset_muls();
        f.inv(&E5([3, 1, 4, 1, 5]));
        assert_eq!(f.muls(), 133);
    }

    #[test]
    fn s5_vanishes_on_four_point_sums_and_the_symmetrisation_reproduces_it() {
        let mut rng = StdRng::seed_from_u64(2);
        let curve = Curve5::random(31, &mut rng);
        let f = &curve.f;
        let mut pts: Vec<Pt5> = Vec::new();
        while pts.len() < 4 {
            let x = f.random(&mut rng);
            if let Some(q) = curve.lift_x(&x, &mut rng) {
                pts.push(q);
            }
        }
        let sum = pts.iter().fold(Pt5::INF, |acc, q| curve.add(&acc, q));
        assert!(curve.on_curve(&sum));
        let x = [pts[0].x, pts[1].x, pts[2].x, pts[3].x, sum.x];
        assert!(curve.s5(&x).is_zero());
        let pre = SymmetrisedS5Q::precompute(&curve, &mut rng);
        assert_eq!(pre.monos.len(), 495);
        let e = elementary4(f, &[pts[0].x, pts[1].x, pts[2].x, pts[3].x]);
        assert!(pre.eval(f, &e, &sum.x).is_zero());
        let y: [E5; 4] = core::array::from_fn(|_| f.random(&mut rng));
        assert!(!pre.eval(f, &elementary4(f, &y), &sum.x).is_zero());
    }

    #[test]
    fn prime_order_instances_and_the_oracle_find_planted_quadruples() {
        let inst = generate_instance5(31, 3);
        assert!(is_prime_u64(inst.n));
        assert!(inst.curve.mul(&inst.g, inst.n).inf);
        let mut rng = StdRng::seed_from_u64(4);
        let base = factor_base(&inst.curve, &mut rng);
        let table = pair_table(&inst.curve, &base);
        let curve = &inst.curve;
        for _ in 0..4 {
            let idx: [usize; 4] = loop {
                let idx: [usize; 4] = core::array::from_fn(|_| rng.gen_range(0..base.len()));
                let mut s = idx.to_vec();
                s.sort_unstable();
                s.dedup();
                if s.len() == 4 {
                    break idx;
                }
            };
            let r = idx
                .iter()
                .fold(Pt5::INF, |acc, &i| curve.add(&acc, &base[i]));
            let found = mitm4_signed(curve, &base, &table, &r);
            let mut want = [(idx[0], 1i64), (idx[1], 1), (idx[2], 1), (idx[3], 1)];
            want.sort_unstable();
            assert!(found.contains(&want), "{want:?} not in {found:?}");
            assert!(found.iter().all(|t| verify_quad(curve, &base, &r, t)));
        }
    }
}
