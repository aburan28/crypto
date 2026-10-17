//! **Gaudry's index calculus at `k = 4`: the `S₅` solve, measured before
//! the pipeline.**
//!
//! Companion to `research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md`
//! §11.16–11.17.  §11.16 derived everything the `k = 4` crossover with rho
//! depends on except `C₄`, the cost of one solve of the symmetrised `S₅`
//! system over `F_{p⁴}`, and bounded it below at `7.1·10¹⁰` `F_p`
//! multiplications from the characteristic polynomial alone.  This module
//! measures it, in the design of [`super::gaudry_cubic`] (§11.4–11.6):
//! Weil restriction of the symmetrised summation polynomial, a Macaulay
//! matrix at the regularity degree with Macaulay's row selection, forward
//! elimination, memoised normal forms, the multiplication matrix of `e₁`,
//! its characteristic polynomial and eigenvectors, and the split of the
//! quartic whose roots are the four abscissae.
//!
//! It builds nothing else of the pipeline — no relation collection, no
//! linear algebra — and every solve is checked against an independent
//! meet-in-the-middle decomposition oracle ([`mitm_decompositions`]).
//!
//! `p ≡ 1 (mod 4)` throughout: `F_{p⁴} = F_p[t]/(t⁴ − c)` needs it, since for
//! `p ≡ 3 (mod 4)` every binomial `t⁴ − c` is reducible.

use std::cell::Cell;
use std::collections::{BTreeSet, HashMap};
use std::time::Instant;

use rand::rngs::StdRng;
use rand::Rng;
use rayon::prelude::*;
use serde::Serialize;

use super::gaudry_cubic::{
    am, charpoly_mod_p, eigenvectors_mod_p, mm, rref_mod_p, sm, split_eigenspace, PolyRing, UPoly,
};
use super::residual_walk::{inv_mod, pow_mod};

// ── F_{p⁴} = F_p[t] / (t⁴ − c) ───────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default, PartialOrd, Ord)]
pub struct E4(pub [u64; 4]);

impl E4 {
    pub const ZERO: E4 = E4([0; 4]);
    pub const ONE: E4 = E4([1, 0, 0, 0]);
    pub fn is_zero(&self) -> bool {
        self.0 == [0; 4]
    }
    pub fn in_fp(&self) -> bool {
        self.0[1] == 0 && self.0[2] == 0 && self.0[3] == 0
    }
}

/// `F_{p⁴}` with a multiplication counter in `F_p` multiplications, charged
/// by the convention [`super::gaudry_cubic::Fp3`] uses at `k = 3`.
#[derive(Clone)]
pub struct Fp4 {
    pub p: u64,
    pub c: u64,
    muls: Cell<u64>,
}

impl Fp4 {
    pub fn new(p: u64) -> Fp4 {
        assert!(
            p % 4 == 1,
            "t⁴ − c is irreducible over F_p only when p ≡ 1 (mod 4)"
        );
        let c = (2..p)
            .find(|&c| pow_mod(c, (p - 1) / 2, p) == p - 1)
            .expect("a non-square");
        Fp4 {
            p,
            c,
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
    pub fn from_fp(&self, x: u64) -> E4 {
        E4([x % self.p, 0, 0, 0])
    }
    pub fn add(&self, a: &E4, b: &E4) -> E4 {
        E4(core::array::from_fn(|i| am(a.0[i], b.0[i], self.p)))
    }
    pub fn sub(&self, a: &E4, b: &E4) -> E4 {
        E4(core::array::from_fn(|i| sm(a.0[i], b.0[i], self.p)))
    }
    pub fn neg(&self, a: &E4) -> E4 {
        E4(core::array::from_fn(|i| (self.p - a.0[i]) % self.p))
    }
    /// Multiplication by an element of `F_p`: four multiplications.
    pub fn scale(&self, a: &E4, k: u64) -> E4 {
        self.count(4);
        E4(core::array::from_fn(|i| mm(a.0[i], k, self.p)))
    }
    /// Schoolbook with `t⁴ = c`: 16 products and 3 reductions.
    pub fn mul(&self, a: &E4, b: &E4) -> E4 {
        let p = self.p;
        self.count(19);
        let mut d = [0u64; 7];
        for i in 0..4 {
            for j in 0..4 {
                d[i + j] = am(d[i + j], mm(a.0[i], b.0[j], p), p);
            }
        }
        E4([
            am(d[0], mm(self.c, d[4], p), p),
            am(d[1], mm(self.c, d[5], p), p),
            am(d[2], mm(self.c, d[6], p), p),
            d[3],
        ])
    }
    pub fn sq(&self, a: &E4) -> E4 {
        self.mul(a, a)
    }
    /// Inverse through the norm to `F_{p²} = F_p[s]/(s² − c)`, `s = t²`:
    /// with `a = A + B t`, `a⁻¹ = (A − B t) / (A² − s B²)`.  Twenty-four
    /// multiplications and one inversion in `F_p`, charged `16` as the count
    /// of `30` for `Fp3::inv` implies (it performs `14`): `40` in all, the
    /// figure §11.16's `97` per addition rests on.
    pub fn inv(&self, a: &E4) -> E4 {
        let p = self.p;
        let c = self.c;
        let (a0, a1, a2, a3) = (a.0[0], a.0[1], a.0[2], a.0[3]);
        let sqr = |x0: u64, x1: u64| -> (u64, u64) {
            let x01 = mm(x0, x1, p);
            (
                am(mm(x0, x0, p), mm(c, mm(x1, x1, p), p), p),
                am(x01, x01, p),
            )
        };
        let m2 = |x0: u64, x1: u64, y0: u64, y1: u64| -> (u64, u64) {
            (
                am(mm(x0, y0, p), mm(c, mm(x1, y1, p), p), p),
                am(mm(x0, y1, p), mm(x1, y0, p), p),
            )
        };
        let (aa0, aa1) = sqr(a0, a2);
        let (bb0, bb1) = sqr(a1, a3);
        // s·B² = c·bb1 + bb0·s.
        let (n0, n1) = (sm(aa0, mm(c, bb1, p), p), sm(aa1, bb0, p));
        let d = sm(mm(n0, n0, p), mm(c, mm(n1, n1, p), p), p);
        assert!(d != 0, "inverse of zero in F_p⁴");
        let dinv = inv_mod(d, p);
        let (i0, i1) = (mm(n0, dinv, p), (p - mm(n1, dinv, p)) % p);
        let (x0, x2) = m2(a0, a2, i0, i1);
        let (y0, y2) = m2(a1, a3, i0, i1);
        self.count(40);
        E4([x0, (p - y0) % p, x2, (p - y2) % p])
    }
    pub fn pow(&self, a: &E4, mut e: u128) -> E4 {
        let mut base = *a;
        let mut acc = E4::ONE;
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
        (self.p as u128).pow(4) - 1
    }
    pub fn is_square(&self, a: &E4) -> bool {
        a.is_zero() || self.pow(a, self.order() / 2) == E4::ONE
    }
    pub fn random(&self, rng: &mut StdRng) -> E4 {
        E4(core::array::from_fn(|_| rng.gen_range(0..self.p)))
    }
    /// Tonelli–Shanks over `F_{p⁴}^×`.
    pub fn sqrt(&self, a: &E4, rng: &mut StdRng) -> Option<E4> {
        if a.is_zero() {
            return Some(E4::ZERO);
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
        while t != E4::ONE {
            let mut i = 0;
            let mut tt = t;
            while tt != E4::ONE {
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

// ── The curve `y² = x³ + a x + b` over `F_{p⁴}` ─────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Pt4 {
    pub x: E4,
    pub y: E4,
    pub inf: bool,
}

impl Pt4 {
    pub const INF: Pt4 = Pt4 {
        x: E4::ZERO,
        y: E4::ZERO,
        inf: true,
    };
}

#[derive(Clone)]
pub struct Curve4 {
    pub f: Fp4,
    pub a: E4,
    pub b: E4,
    /// Group operations (additions and doublings) performed.
    ops: Cell<u64>,
}

impl Curve4 {
    /// Coefficients drawn outside `F_{p²}` so that the curve is not a
    /// twist of one defined over a subfield.
    pub fn random(p: u64, rng: &mut StdRng) -> Curve4 {
        let f = Fp4::new(p);
        loop {
            let a = f.random(rng);
            let b = f.random(rng);
            if a.0[1] == 0 && a.0[3] == 0 || b.0[1] == 0 && b.0[3] == 0 {
                continue;
            }
            let a3 = f.mul(&f.sq(&a), &a);
            let disc = f.add(&f.scale(&a3, 4), &f.scale(&f.sq(&b), 27));
            if !disc.is_zero() {
                return Curve4 {
                    f,
                    a,
                    b,
                    ops: Cell::new(0),
                };
            }
        }
    }
    pub fn rhs(&self, x: &E4) -> E4 {
        let f = &self.f;
        f.add(&f.mul(&f.add(&f.sq(x), &self.a), x), &self.b)
    }
    pub fn lift_x(&self, x: &E4, rng: &mut StdRng) -> Option<Pt4> {
        let y = self.f.sqrt(&self.rhs(x), rng)?;
        Some(Pt4 {
            x: *x,
            y,
            inf: false,
        })
    }
    pub fn on_curve(&self, pt: &Pt4) -> bool {
        pt.inf || self.f.sq(&pt.y) == self.rhs(&pt.x)
    }
    pub fn neg(&self, pt: &Pt4) -> Pt4 {
        if pt.inf {
            return *pt;
        }
        Pt4 {
            x: pt.x,
            y: self.f.neg(&pt.y),
            inf: false,
        }
    }
    /// Affine addition: one inversion and three multiplications for
    /// distinct points, the decomposition behind §11.16's `c_add`.
    pub fn add(&self, p1: &Pt4, p2: &Pt4) -> Pt4 {
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
                return Pt4::INF;
            }
            // (3x² + a) / 2y
            let x2 = f.sq(&p1.x);
            let num = f.add(&f.add(&f.add(&x2, &x2), &x2), &self.a);
            f.mul(&num, &f.inv(&f.add(&p1.y, &p1.y)))
        } else {
            f.mul(&f.sub(&p2.y, &p1.y), &f.inv(&f.sub(&p2.x, &p1.x)))
        };
        let x3 = f.sub(&f.sub(&f.sq(&lambda), &p1.x), &p2.x);
        let y3 = f.sub(&f.mul(&lambda, &f.sub(&p1.x, &x3)), &p1.y);
        Pt4 {
            x: x3,
            y: y3,
            inf: false,
        }
    }
    pub fn sub(&self, p1: &Pt4, p2: &Pt4) -> Pt4 {
        self.add(p1, &self.neg(p2))
    }
    pub fn ops(&self) -> u64 {
        self.ops.get()
    }
    /// A fresh curve with its own counters, for a worker thread.
    pub fn with_params(p: u64, a: E4, b: E4) -> Curve4 {
        Curve4 {
            f: Fp4::new(p),
            a,
            b,
            ops: Cell::new(0),
        }
    }
    pub fn reset_ops(&self) {
        self.ops.set(0);
    }
    /// Double-and-add.
    pub fn mul(&self, pt: &Pt4, mut k: u64) -> Pt4 {
        let mut acc = Pt4::INF;
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

/// Univariate polynomials over `F_{p⁴}`, low to high.
fn poly_mul(f: &Fp4, a: &[E4], b: &[E4]) -> Vec<E4> {
    let mut out = vec![E4::ZERO; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = f.add(&out[i + j], &f.mul(x, y));
        }
    }
    out
}
fn poly_sub(f: &Fp4, a: &[E4], b: &[E4]) -> Vec<E4> {
    let n = a.len().max(b.len());
    (0..n)
        .map(|i| f.sub(a.get(i).unwrap_or(&E4::ZERO), b.get(i).unwrap_or(&E4::ZERO)))
        .collect()
}
fn poly_scale(f: &Fp4, a: &[E4], k: &E4) -> Vec<E4> {
    a.iter().map(|x| f.mul(x, k)).collect()
}

/// Determinant over `F_{p⁴}` by Gaussian elimination.
fn det4(f: &Fp4, mut m: Vec<Vec<E4>>) -> E4 {
    let n = m.len();
    let mut det = E4::ONE;
    for c in 0..n {
        let Some(r) = (c..n).find(|&r| !m[r][c].is_zero()) else {
            return E4::ZERO;
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

impl Curve4 {
    /// `S₃(x₁, x₂, X)` as a polynomial in `X` (Semaev):
    /// `(x₁ − x₂)² X² − 2((x₁ + x₂)(x₁x₂ + a) + 2b) X + (x₁x₂ − a)² − 4b(x₁ + x₂)`.
    fn s3_in_last(&self, x1: &E4, x2: &E4) -> [E4; 3] {
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

    /// `S₄(x₃, x₄, x₅, Y)` as a polynomial in `Y`, of degree ≤ 4:
    /// `Res_Z(S₃(x₃, x₄, Z), S₃(x₅, Y, Z))`, with the resultant of two
    /// quadratics `(a₂b₀ − a₀b₂)² − (a₂b₁ − a₁b₂)(a₁b₀ − a₀b₁)`.
    pub(crate) fn s4_in_last(&self, x3: &E4, x4: &E4, x5: &E4) -> Vec<E4> {
        let f = &self.f;
        let [a0, a1, a2] = self.s3_in_last(x3, x4);
        let (a, b) = (&self.a, &self.b);
        let x5s = f.sq(x5);
        let two = |v: &E4| f.add(v, v);
        let four = |v: &E4| two(&two(v));
        // S₃(x₅, Y, Z) = b₂(Y) Z² + b₁(Y) Z + b₀(Y).
        let b2 = vec![x5s, f.neg(&two(x5)), E4::ONE];
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
    /// recursive definition, as the determinant of the `6 × 6` Sylvester
    /// matrix at formal degrees `2` and `4`.
    pub fn s5(&self, x: &[E4; 5]) -> E4 {
        let f = &self.f;
        let a = self.s3_in_last(&x[0], &x[1]);
        let mut b = self.s4_in_last(&x[2], &x[3], &x[4]);
        b.resize(5, E4::ZERO);
        let mut m = vec![vec![E4::ZERO; 6]; 6];
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
        det4(f, m)
    }
}

// ── The symmetrised S₅, by interpolation ──────────────────────────────────

/// Monomials `e₁^a e₂^b e₃^c e₄^d` of total degree `≤ d`.
pub(crate) fn exponents4(d: u8) -> Vec<[u8; 4]> {
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

/// `S₅(x₁, x₂, x₃, x₄, x_R) = H(e₁, e₂, e₃, e₄, x_R)`, with `H` of total
/// degree `≤ 8` in the elementary symmetric `e`'s and degree `≤ 8` in `x_R`:
/// `495 · 9` coefficients in `F_{p⁴}`.  Found by interpolation rather than
/// symbolic symmetrisation — `S₅` evaluated as a resultant at `495` random
/// points and the `495 × 495` system solved once for all nine powers of
/// `x_R` — and then checked against fresh evaluations.
pub struct SymmetrisedS5 {
    pub monos: Vec<[u8; 4]>,
    pub coef: Vec<[E4; 9]>,
    /// `F_p` multiplications spent building it (once per curve).
    pub precompute_muls: u64,
}

fn elementary(f: &Fp4, x: &[E4; 4]) -> [E4; 4] {
    let e1 = f.add(&f.add(&x[0], &x[1]), &f.add(&x[2], &x[3]));
    let mut e2 = E4::ZERO;
    let mut e3 = E4::ZERO;
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

fn mono_values(f: &Fp4, e: &[E4; 4], monos: &[[u8; 4]]) -> Vec<E4> {
    let mut pows = [[E4::ONE; 9]; 4];
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

impl SymmetrisedS5 {
    /// The nine coefficients of `S₅(x₁, …, x₄, X)` in `X`, from values at
    /// `X = 0, …, 8` through the inverse Vandermonde matrix over `F_p`.
    fn coefficients_in_last(curve: &Curve4, x: &[E4; 4], vinv: &[[u64; 9]; 9]) -> [E4; 9] {
        let f = &curve.f;
        let vals: Vec<E4> = (0..9u64)
            .map(|k| curve.s5(&[x[0], x[1], x[2], x[3], f.from_fp(k)]))
            .collect();
        core::array::from_fn(|j| {
            vals.iter()
                .enumerate()
                .fold(E4::ZERO, |acc, (k, v)| f.add(&acc, &f.scale(v, vinv[j][k])))
        })
    }

    pub fn precompute(curve: &Curve4, rng: &mut StdRng) -> SymmetrisedS5 {
        let f = &curve.f;
        let p = f.p;
        let before = f.muls();
        // Inverse Vandermonde at the nodes 0..8, over F_p.
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
            // Row j of the reduced matrix is [I | V⁻¹]; V maps coefficients to
            // values, so coefficient j = Σ_k V⁻¹[j][k] · value_k.
            core::array::from_fn(|j| core::array::from_fn(|k| m[j][9 + k]))
        };
        let monos = exponents4(8);
        let n = monos.len();
        loop {
            // Augmented system [V | C] over F_{p⁴}: one row per sample.
            let mut sys: Vec<Vec<E4>> = Vec::with_capacity(n);
            for _ in 0..n {
                let x: [E4; 4] = core::array::from_fn(|_| f.random(rng));
                let e = elementary(f, &x);
                let mut row = mono_values(f, &e, &monos);
                row.extend(Self::coefficients_in_last(curve, &x, &vinv));
                sys.push(row);
            }
            // Gauss–Jordan.
            let mut ok = true;
            for c in 0..n {
                let Some(r) = (c..n).find(|&r| !sys[r][c].is_zero()) else {
                    ok = false;
                    break;
                };
                sys.swap(r, c);
                let inv = f.inv(&sys[c][c]);
                let pivot: Vec<E4> = sys[c].iter().map(|v| f.mul(v, &inv)).collect();
                sys[c] = pivot.clone();
                sys.par_iter_mut().enumerate().for_each(|(r, row)| {
                    if r == c || row[c].is_zero() {
                        return;
                    }
                    // A local field for the parallel rows: the counter is
                    // charged below, from the operation count.
                    let fl = Fp4::new(p);
                    let k = row[c];
                    for j in c..row.len() {
                        if !pivot[j].is_zero() {
                            let v = fl.mul(&k, &pivot[j]);
                            row[j] = fl.sub(&row[j], &v);
                        }
                    }
                });
                f.count(19 * (n as u64) * (n + 9 - c) as u64);
            }
            if !ok {
                continue;
            }
            let coef: Vec<[E4; 9]> = (0..n)
                .map(|r| core::array::from_fn(|j| sys[r][n + j]))
                .collect();
            let pre = SymmetrisedS5 {
                monos: monos.clone(),
                coef,
                precompute_muls: 0,
            };
            // Check against fresh evaluations of the resultant.
            let good = (0..8).all(|_| {
                let x: [E4; 4] = core::array::from_fn(|_| f.random(rng));
                let want = Self::coefficients_in_last(curve, &x, &vinv);
                let e = elementary(f, &x);
                let mv = mono_values(f, &e, &pre.monos);
                (0..9).all(|j| {
                    let got = mv
                        .iter()
                        .zip(&pre.coef)
                        .fold(E4::ZERO, |acc, (m, c)| f.add(&acc, &f.mul(m, &c[j])));
                    got == want[j]
                })
            });
            assert!(good, "the symmetrised S₅ does not reproduce the resultant");
            return SymmetrisedS5 {
                precompute_muls: f.muls() - before,
                ..pre
            };
        }
    }

    /// The four `F_p`-components of `H(e₁, …, e₄, x_R)`, as sparse maps
    /// over the monomials of total degree `≤ 8`.
    pub fn weil_restrict(&self, f: &Fp4, x_r: &E4) -> [HashMap<[u8; 4], u64>; 4] {
        let mut pows = [E4::ONE; 9];
        for k in 1..9 {
            pows[k] = f.mul(&pows[k - 1], x_r);
        }
        let mut out: [HashMap<[u8; 4], u64>; 4] = Default::default();
        for (m, c) in self.monos.iter().zip(&self.coef) {
            let v = (0..9).fold(E4::ZERO, |acc, j| f.add(&acc, &f.mul(&c[j], &pows[j])));
            for (i, comp) in out.iter_mut().enumerate() {
                if v.0[i] != 0 {
                    comp.insert(*m, v.0[i]);
                }
            }
        }
        out
    }
}

// ── Monomials in four variables ───────────────────────────────────────────

/// The degree of `m` when unknown `k` has weight `w[k]`: the total degree
/// for weights all 1.
fn wdeg<const N: usize>(m: &[u8; N], w: &[u8; N]) -> u16 {
    m.iter().zip(w).map(|(&x, &k)| x as u16 * k as u16).sum()
}

/// Weighted degree, then reverse lexicographic: grevlex for weights all 1.
fn wgrevlex_cmp<const N: usize>(a: &[u8; N], b: &[u8; N], w: &[u8; N]) -> std::cmp::Ordering {
    let (da, db) = (wdeg(a, w), wdeg(b, w));
    if da != db {
        return da.cmp(&db);
    }
    for i in (0..N).rev() {
        if a[i] != b[i] {
            return b[i].cmp(&a[i]);
        }
    }
    std::cmp::Ordering::Equal
}

/// Monomials of weighted degree `≤ d` in `N` variables (any order).
fn exponents_w<const N: usize>(d: u8, w: &[u8; N]) -> Vec<[u8; N]> {
    let mut out = Vec::new();
    let mut cur = [0u8; N];
    fn rec<const N: usize>(
        i: usize,
        left: u8,
        w: &[u8; N],
        cur: &mut [u8; N],
        out: &mut Vec<[u8; N]>,
    ) {
        if i == N {
            out.push(*cur);
            return;
        }
        for k in 0..=left / w[i] {
            cur[i] = k;
            rec(i + 1, left - k * w[i], w, cur, out);
        }
        cur[i] = 0;
    }
    rec(0, d, w, &mut cur, &mut out);
    out
}

/// Monomials of weighted degree `≤ d`, weighted grevlex descending.
fn monomials_w_desc<const N: usize>(d: u8, w: &[u8; N]) -> Vec<[u8; N]> {
    let mut v = exponents_w::<N>(d, w);
    v.sort_by(|x, y| wgrevlex_cmp(y, x, w));
    v
}

fn leading_monomial_w<const N: usize>(f: &HashMap<[u8; N], u64>, w: &[u8; N]) -> Option<[u8; N]> {
    f.iter()
        .filter(|(_, &c)| c != 0)
        .map(|(m, _)| *m)
        .max_by(|a, b| wgrevlex_cmp(a, b, w))
}

/// As `gaudry_cubic`'s `triangularise_leads`, for four components.
fn triangularise_leads4<const N: usize>(
    comps: &[HashMap<[u8; N], u64>],
    w: &[u8; N],
    p: u64,
    muls: &mut u64,
) -> (Vec<HashMap<[u8; N], u64>>, Vec<Option<[u8; N]>>) {
    let mut fs: Vec<HashMap<[u8; N], u64>> = comps.to_vec();
    for _round in 0..16 {
        fs.sort_by(
            |a, b| match (leading_monomial_w(a, w), leading_monomial_w(b, w)) {
                (Some(x), Some(y)) => wgrevlex_cmp(&y, &x, w),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            },
        );
        let mut changed = false;
        for i in 1..fs.len() {
            let (Some(li), Some(lj)) = (
                leading_monomial_w(&fs[i], w),
                leading_monomial_w(&fs[i - 1], w),
            ) else {
                continue;
            };
            if li == lj {
                let k = mm(fs[i][&li], inv_mod(fs[i - 1][&lj], p), p);
                let prev = fs[i - 1].clone();
                for (&m, &c) in &prev {
                    let e = fs[i].entry(m).or_insert(0);
                    *e = sm(*e, mm(k, c, p), p);
                    *muls += 1;
                }
                fs[i].retain(|_, c| *c != 0);
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    let leads = fs.iter().map(|f| leading_monomial_w(f, w)).collect();
    (fs, leads)
}

// ── Arithmetic on u16 rows ────────────────────────────────────────────────

/// `row ← row − f·piv` elementwise mod `P`, vectorisable when `P` is a
/// compile-time constant; the solver's matrices are dense `u16` rows.
#[inline(always)]
fn axpy_const<const P: u32>(row: &mut [u16], piv: &[u16], f: u32) {
    let nf = P - f;
    for (r, &v) in row.iter_mut().zip(piv) {
        *r = ((*r as u32 + nf * v as u32) % P) as u16;
    }
}

#[inline(always)]
fn axpy_dyn(row: &mut [u16], piv: &[u16], f: u32, p: u32) {
    let nf = p - f;
    for (r, &v) in row.iter_mut().zip(piv) {
        *r = ((*r as u32 + nf * v as u32) % p) as u16;
    }
}

#[inline(always)]
/// Every prime `p ≡ 1 (mod 4)` below 300 gets the constant-modulus loop
/// (the `k = 4` solves run there); others take the run-time remainder.
fn axpy_portable(row: &mut [u16], piv: &[u16], f: u32, p: u32) {
    match p {
        5 => axpy_const::<5>(row, piv, f),
        13 => axpy_const::<13>(row, piv, f),
        17 => axpy_const::<17>(row, piv, f),
        29 => axpy_const::<29>(row, piv, f),
        37 => axpy_const::<37>(row, piv, f),
        41 => axpy_const::<41>(row, piv, f),
        53 => axpy_const::<53>(row, piv, f),
        61 => axpy_const::<61>(row, piv, f),
        73 => axpy_const::<73>(row, piv, f),
        89 => axpy_const::<89>(row, piv, f),
        97 => axpy_const::<97>(row, piv, f),
        101 => axpy_const::<101>(row, piv, f),
        109 => axpy_const::<109>(row, piv, f),
        113 => axpy_const::<113>(row, piv, f),
        137 => axpy_const::<137>(row, piv, f),
        149 => axpy_const::<149>(row, piv, f),
        157 => axpy_const::<157>(row, piv, f),
        173 => axpy_const::<173>(row, piv, f),
        181 => axpy_const::<181>(row, piv, f),
        193 => axpy_const::<193>(row, piv, f),
        197 => axpy_const::<197>(row, piv, f),
        229 => axpy_const::<229>(row, piv, f),
        233 => axpy_const::<233>(row, piv, f),
        241 => axpy_const::<241>(row, piv, f),
        257 => axpy_const::<257>(row, piv, f),
        269 => axpy_const::<269>(row, piv, f),
        277 => axpy_const::<277>(row, piv, f),
        281 => axpy_const::<281>(row, piv, f),
        293 => axpy_const::<293>(row, piv, f),
        _ => axpy_dyn(row, piv, f, p),
    }
}

/// The same loop compiled for wider vectors, chosen at run time: the
/// release build targets baseline x86-64, where the multiply and the
/// constant remainder barely vectorise.  Integer arithmetic, so the result
/// is the same whichever runs.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vl,avx2")]
unsafe fn axpy_avx512(row: &mut [u16], piv: &[u16], f: u32, p: u32) {
    axpy_portable(row, piv, f, p)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn axpy_avx2(row: &mut [u16], piv: &[u16], f: u32, p: u32) {
    axpy_portable(row, piv, f, p)
}

#[inline]
pub(crate) fn axpy(row: &mut [u16], piv: &[u16], f: u32, p: u32) {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx512bw") && is_x86_feature_detected!("avx512vl") {
            // SAFETY: the features were detected on this CPU.
            return unsafe { axpy_avx512(row, piv, f, p) };
        }
        if is_x86_feature_detected!("avx2") {
            // SAFETY: as above.
            return unsafe { axpy_avx2(row, piv, f, p) };
        }
    }
    axpy_portable(row, piv, f, p)
}

// ── Dense kernels on u16 rows ─────────────────────────────────────────────
//
// The quotient-side routines (characteristic polynomial, eigenvectors) do
// exactly the operations of `gaudry_cubic`'s `charpoly_mod_p` and
// `eigenvectors_mod_p`, charge the same count and return the same values;
// they run on `u16` rows with the vectorised `axpy` and in parallel over
// rows instead of through `u128` remainders one entry at a time.

/// As `gaudry_cubic::rref_mod_p`, on `u16` rows, rows eliminated in parallel.
fn rref_u16(m: &mut [Vec<u16>], p: u32, muls: &mut u64) -> Vec<usize> {
    let rows = m.len();
    let cols = if rows > 0 { m[0].len() } else { 0 };
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..cols {
        if r >= rows {
            break;
        }
        let Some(pr) = (r..rows).find(|&i| m[i][c] != 0) else {
            continue;
        };
        m.swap(r, pr);
        let inv = inv_mod(m[r][c] as u64, p as u64) as u32;
        for v in &mut m[r][c..] {
            *v = ((*v as u32 * inv) % p) as u16;
        }
        *muls += (cols - c) as u64;
        let pivot = m[r].clone();
        let hit = m
            .par_iter_mut()
            .enumerate()
            .filter(|(i, row)| *i != r && row[c] != 0)
            .map(|(_, row)| {
                let f = row[c] as u32;
                axpy(&mut row[c..], &pivot[c..], f, p);
            })
            .count();
        *muls += (hit * (cols - c)) as u64;
        pivots.push(c);
        r += 1;
    }
    pivots
}

/// As `gaudry_cubic::charpoly_mod_p`.  The Hessenberg step at column `j`
/// is a similarity by `I + Σ fᵢ E_{i,j+1}`, whose elementary factors
/// commute (their products vanish), so its row operations can all run
/// before its column operations: the same matrix as the interleaved order.
pub(crate) fn charpoly_u16(a: &[Vec<u64>], p: u32, muls: &mut u64) -> UPoly {
    let n = a.len();
    let mut h: Vec<Vec<u16>> = a
        .iter()
        .map(|r| r.iter().map(|&v| v as u16).collect())
        .collect();
    for j in 0..n.saturating_sub(2) {
        let Some(piv) = ((j + 1)..n).find(|&i| h[i][j] != 0) else {
            continue;
        };
        if piv != j + 1 {
            h.swap(piv, j + 1);
            for row in h.iter_mut() {
                row.swap(piv, j + 1);
            }
        }
        let inv = inv_mod(h[j + 1][j] as u64, p as u64) as u32;
        let fs: Vec<(usize, u32)> = ((j + 2)..n)
            .filter(|&i| h[i][j] != 0)
            .map(|i| (i, (h[i][j] as u32 * inv) % p))
            .collect();
        if fs.is_empty() {
            continue;
        }
        *muls += 2 * (n * fs.len()) as u64;
        let mut fv = vec![0u32; n];
        for &(i, f) in &fs {
            fv[i] = f;
        }
        // row_i −= fᵢ · row_{j+1}
        let pivot = h[j + 1].clone();
        h.par_iter_mut().enumerate().for_each(|(i, row)| {
            if fv[i] != 0 {
                axpy(row, &pivot, fv[i], p);
            }
        });
        // col_{j+1} += Σ fᵢ · col_i
        h.par_iter_mut().for_each(|row| {
            let s: u64 = fs.iter().map(|&(i, f)| f as u64 * row[i] as u64).sum();
            row[j + 1] = ((row[j + 1] as u64 + s) % p as u64) as u16;
        });
    }
    // p_0 = 1; p_k(x) = (x − h[k-1][k-1]) p_{k-1}
    //                   − Σ_{i=1}^{k-1} h[i-1][k-1] (Π_{m=i}^{k-1} h[m][m-1]) p_{i-1},
    // charged as `PolyRing` charges its product and scalings.
    let pp = p as u64;
    let mut polys: Vec<Vec<u32>> = vec![vec![1]];
    let mut acc: Vec<u64> = Vec::with_capacity(n);
    for k in 1..=n {
        let prev = &polys[k - 1];
        let hk = (pp - h[k - 1][k - 1] as u64) % pp;
        let mut pk = vec![0u32; k + 1];
        for t in 0..=k {
            let hi = if t >= 1 { prev[t - 1] as u64 } else { 0 };
            let lo = if t < k { hk * prev[t] as u64 } else { 0 };
            pk[t] = ((hi + lo) % pp) as u32;
        }
        *muls += 2 * k as u64;
        acc.clear();
        acc.resize(k, 0);
        let mut prod = 1u64;
        for i in (1..k).rev() {
            prod = prod * h[i][i - 1] as u64 % pp;
            let coef = h[i - 1][k - 1] as u64 * prod % pp;
            if coef != 0 {
                for (s, &q) in acc.iter_mut().zip(&polys[i - 1]) {
                    *s += coef * q as u64;
                }
                *muls += i as u64;
            }
        }
        for (v, &s) in pk.iter_mut().zip(&acc) {
            *v = ((*v as u64 + pp - s % pp) % pp) as u32;
        }
        polys.push(pk);
    }
    UPoly(polys.pop().unwrap().into_iter().map(u64::from).collect())
}

/// As `gaudry_cubic::eigenvectors_mod_p`: a kernel basis of `m − λI`.
pub(crate) fn eigenvectors_u16(
    m: &[Vec<u16>],
    lambda: u64,
    p: u32,
    muls: &mut u64,
) -> Vec<Vec<u64>> {
    let n = m.len();
    let mut a: Vec<Vec<u16>> = m.to_vec();
    for (i, row) in a.iter_mut().enumerate() {
        row[i] = ((row[i] as u64 + p as u64 - lambda) % p as u64) as u16;
    }
    let pivots = rref_u16(&mut a, p, muls);
    let mut is_piv = vec![false; n];
    for &c in &pivots {
        is_piv[c] = true;
    }
    let mut basis = Vec::new();
    for fc in (0..n).filter(|&c| !is_piv[c]) {
        let mut v = vec![0u64; n];
        v[fc] = 1;
        for (r, &pc) in pivots.iter().enumerate() {
            v[pc] = (p as u64 - a[r][fc] as u64) % p as u64;
        }
        basis.push(v);
    }
    basis
}

// ── The solve ─────────────────────────────────────────────────────────────

/// Everything one solve cost, by phase, in `F_p` multiplications counted as
/// `gaudry_cubic` counts them: only multiplications actually performed in
/// the elimination and the normal forms, the characteristic polynomial and
/// eigenvectors through the same routines.
#[derive(Clone, Debug, Default, Serialize)]
pub struct QuarticSolve {
    pub equation_degree: u8,
    pub macaulay_degree: u8,
    pub rows: usize,
    pub cols: usize,
    pub pivots: usize,
    pub zero_rows: usize,
    pub dim: usize,
    pub reached_pivots: usize,
    pub weil_muls: u64,
    pub triangularise_muls: u64,
    pub echelon_muls: u64,
    pub nf_muls: u64,
    pub charpoly_muls: u64,
    pub roots_muls: u64,
    pub eigenvector_muls: u64,
    pub verify_muls: u64,
    pub split_muls: u64,
    pub total_muls: u64,
    pub rational_eigenvalues: usize,
    pub degenerate_eigenspaces: usize,
    pub e_solutions: usize,
    pub build_ms: f64,
    pub echelon_ms: f64,
    pub nf_ms: f64,
    pub charpoly_ms: f64,
    pub eigen_ms: f64,
    pub total_ms: f64,
    pub outcome: String,
    pub replay: ReplayOutcome,
    /// Held rows a replay reduced into its echelon to reach the recorded
    /// rank (`top_up`): 0 when the replayed rows had it.
    #[serde(skip_serializing_if = "is_zero")]
    pub topped_up: u64,
}

fn is_zero(v: &u64) -> bool {
    *v == 0
}

impl QuarticSolve {
    /// Sum the phases into `total_muls` (for solves built outside this
    /// module).
    pub fn close_public(&mut self) {
        self.close();
    }
    fn close(&mut self) {
        self.total_muls = self.weil_muls
            + self.triangularise_muls
            + self.echelon_muls
            + self.nf_muls
            + self.charpoly_muls
            + self.roots_muls
            + self.eigenvector_muls
            + self.verify_muls
            + self.split_muls;
    }
}

pub(crate) struct Echelon {
    pub(crate) rows: Vec<Vec<u16>>,
    pub(crate) pivot_row: Vec<u32>,
    pub(crate) zero_rows: usize,
    /// Input rows that reduced to zero (or were zero), by index.
    pub(crate) zero_ids: Vec<u32>,
}

/// Forward elimination by leading column.  Every unused row has zeros
/// left of its lead, so the rows with a nonzero entry in column `c` are
/// exactly those whose lead is `c`: the lowest-indexed becomes the pivot,
/// is normalised, and is subtracted from the rest, which then move to the
/// bucket of their new lead.  The same operations as a column-by-column
/// elimination, without scanning rows that are already zero there;
/// `muls` counts the nonzero entries of each pivot row once per row it is
/// subtracted from, and once for its normalisation.
pub(crate) fn echelon_by_lead(
    mut rows: Vec<Vec<u16>>,
    p: u32,
    muls: &mut u64,
    debug: bool,
) -> Echelon {
    let ncols = rows.first().map_or(0, |r| r.len());
    let mut buckets: Vec<Vec<u32>> = vec![Vec::new(); ncols];
    let mut zero_rows = 0usize;
    let mut zero_ids: Vec<u32> = Vec::new();
    for (r, row) in rows.iter().enumerate() {
        match row.iter().position(|&x| x != 0) {
            Some(c) => buckets[c].push(r as u32),
            None => {
                zero_rows += 1;
                zero_ids.push(r as u32);
            }
        }
    }
    // Nonzero entries of each row (kept only for `sparsest`).
    let mut row_nnz: Vec<u32> = if sparsest {
        rows.par_iter()
            .map(|r| r.iter().filter(|&&x| x != 0).count() as u32)
            .collect()
    } else {
        Vec::new()
    };
    let mut pivot_row = vec![u32::MAX; ncols];
    let start = Instant::now();
    let mut npiv = 0usize;
    for c in 0..ncols {
        let mut b = std::mem::take(&mut buckets[c]);
        if b.is_empty() {
            continue;
        }
        b.sort_unstable();
        if sparsest && b.len() > 1 {
            let best = (0..b.len())
                .min_by_key(|&i| (row_nnz[b[i] as usize], i))
                .unwrap_or(0);
            b.swap(0, best);
        }
        let pr = b[0] as usize;
        let mut piv = std::mem::take(&mut rows[pr]);
        let inv = inv_mod(piv[c] as u64, p as u64) as u32;
        let mut nz: Vec<u32> = Vec::new();
        for j in c..ncols {
            if piv[j] != 0 {
                piv[j] = ((piv[j] as u32 * inv) % p) as u16;
                nz.push(j as u32);
            }
        }
        *muls += nz.len() as u64;
        let dense = nz.len() * 8 > ncols - c;
        let mut others: Vec<(u32, Vec<u16>, Option<usize>, u32)> = b[1..]
            .iter()
            .map(|&r| (r, std::mem::take(&mut rows[r as usize]), None, 0))
            .collect();
        *muls += (others.len() * nz.len()) as u64;
        others.par_iter_mut().for_each(|(_, row, lead, nnz)| {
            let f = row[c] as u32;
            if dense {
                axpy(&mut row[c..], &piv[c..], f, p);
            } else {
                let nf = p - f;
                for &j in &nz {
                    let j = j as usize;
                    row[j] = ((row[j] as u32 + nf * piv[j] as u32) % p) as u16;
                }
            }
            *lead = row[c + 1..].iter().position(|&x| x != 0).map(|k| k + c + 1);
            if sparsest {
                if let Some(l) = *lead {
                    *nnz = row[l..].iter().filter(|&&x| x != 0).count() as u32;
                }
            }
        });
        for (r, row, lead, nnz) in others {
            if sparsest {
                row_nnz[r as usize] = nnz;
            }
            rows[r as usize] = row;
            match lead {
                Some(l) => buckets[l].push(r),
                None => {
                    zero_rows += 1;
                    zero_ids.push(r);
                }
            }
        }
        rows[pr] = piv;
        pivot_row[c] = pr as u32;
        npiv += 1;
        if debug && npiv.is_multiple_of(2000) {
            eprintln!(
                "    echelon: {npiv} pivots at column {c}/{ncols}, {:.3e} muls, {:.0} s",
                *muls as f64,
                start.elapsed().as_secs_f64()
            );
        }
    }
    zero_ids.sort_unstable();
    Echelon {
        rows,
        pivot_row,
        zero_rows,
        zero_ids,
    }
}

/// Reduce `held` rows, in order, into a forward echelon until `need` new
/// pivots are found (or the rows run out): each is reduced left to right
/// by the pivot rows (one multiplication per nonzero entry of a pivot row
/// used) and, if it survives, normalised (one per nonzero) as the pivot of
/// its lead.  The echelon then spans what the full row set spans up to the
/// rows not reached.  Returns the number of held rows reduced.
fn top_up_echelon(
    ech: &mut Echelon,
    held: Vec<Vec<u16>>,
    mut need: usize,
    p: u32,
    muls: &mut u64,
) -> u64 {
    let ncols = ech.pivot_row.len();
    let mut used = 0u64;
    for mut v in held {
        if need == 0 {
            break;
        }
        used += 1;
        let mut lead = None;
        for c in 0..ncols {
            if v[c] == 0 {
                continue;
            }
            let pr = ech.pivot_row[c];
            if pr == u32::MAX {
                lead = Some(c);
                break;
            }
            let f = v[c] as u32;
            let piv = &ech.rows[pr as usize];
            *muls += piv[c..].iter().filter(|&&x| x != 0).count() as u64;
            axpy(&mut v[c..], &piv[c..], f, p);
        }
        let Some(c) = lead else {
            ech.zero_rows += 1;
            continue;
        };
        let inv = inv_mod(v[c] as u64, p as u64) as u32;
        for x in &mut v[c..] {
            if *x != 0 {
                *x = ((*x as u32 * inv) % p) as u16;
                *muls += 1;
            }
        }
        ech.pivot_row[c] = ech.rows.len() as u32;
        ech.rows.push(v);
        need -= 1;
    }
    used
}

/// Memoised normal forms over the standard monomials, computed for a set
/// of target columns by one sweep from the right over the pivots they
/// reach — iteratively, since a chain of pivots can be tens of thousands
/// deep.
pub(crate) struct NormalForms4<'a> {
    pub(crate) ech: &'a Echelon,
    pub(crate) is_std: Vec<bool>,
    pub(crate) std_index: Vec<u32>,
    pub(crate) dim: usize,
    pub(crate) p: u32,
    pub(crate) memo: Vec<Option<(Vec<u16>, u64)>>,
    pub(crate) reached: usize,
}

impl<'a> NormalForms4<'a> {
    /// Make every column in `targets` resolvable.  `Err(c)` names a column
    /// reached that is neither standard nor a pivot.
    pub(crate) fn ensure(&mut self, targets: &[usize], muls: &mut u64) -> Result<(), usize> {
        let ncols = self.is_std.len();
        let mut need: Vec<usize> = Vec::new();
        let mut seen = vec![false; ncols];
        let mut stack: Vec<usize> = Vec::new();
        for &t in targets {
            if !self.is_std[t] && self.memo[t].is_none() && !seen[t] {
                seen[t] = true;
                stack.push(t);
            }
        }
        while let Some(c) = stack.pop() {
            let pr = self.ech.pivot_row[c];
            if pr == u32::MAX {
                return Err(c);
            }
            need.push(c);
            let row = &self.ech.rows[pr as usize];
            for j in (c + 1)..ncols {
                if row[j] != 0 && !self.is_std[j] && self.memo[j].is_none() && !seen[j] {
                    seen[j] = true;
                    stack.push(j);
                }
            }
        }
        need.sort_unstable_by(|a, b| b.cmp(a));
        // Right to left in blocks: within a block, the terms on columns
        // already resolved are summed in parallel; the terms on columns of
        // the block itself are added afterwards, right to left.  The same
        // products as one column at a time, so the same vectors and count.
        const BLOCK: usize = 128;
        let p = self.p;
        let mut in_block = vec![false; ncols];
        for block in need.chunks(BLOCK) {
            for &c in block {
                in_block[c] = true;
            }
            let this = &*self;
            let in_block_ref = &in_block;
            let partial: Vec<(Vec<u16>, u64, Vec<(usize, u32)>)> = block
                .par_iter()
                .map(|&c| {
                    let row = &this.ech.rows[this.ech.pivot_row[c] as usize];
                    let mut v = vec![0u16; this.dim];
                    let mut m = 0u64;
                    let mut later = Vec::new();
                    for j in (c + 1)..ncols {
                        let a = row[j] as u32;
                        if a == 0 {
                            continue;
                        }
                        if this.is_std[j] {
                            let i = this.std_index[j] as usize;
                            v[i] = ((v[i] as u32 + p - a) % p) as u16;
                        } else if in_block_ref[j] {
                            later.push((j, a));
                        } else {
                            let (w, nnz) = this.memo[j].as_ref().expect("computed right to left");
                            axpy(&mut v, w, a, p);
                            m += *nnz;
                        }
                    }
                    (v, m, later)
                })
                .collect();
            for (&c, (mut v, m, later)) in block.iter().zip(partial) {
                *muls += m;
                for (j, a) in later {
                    let (w, nnz) = self.memo[j].as_ref().expect("earlier in the block");
                    axpy(&mut v, w, a, p);
                    *muls += *nnz;
                }
                let nnz = v.iter().filter(|&&x| x != 0).count() as u64;
                self.memo[c] = Some((v, nnz));
                self.reached += 1;
            }
            for &c in block {
                in_block[c] = false;
            }
        }
        Ok(())
    }

    pub(crate) fn vector(&self, c: usize) -> Vec<u64> {
        if self.is_std[c] {
            let mut v = vec![0u64; self.dim];
            v[self.std_index[c] as usize] = 1;
            v
        } else {
            let (w, _) = self.memo[c].as_ref().expect("ensured");
            w.iter().map(|&x| x as u64).collect()
        }
    }
}

fn eval4<const N: usize>(
    comp: &HashMap<[u8; N], u64>,
    e: &[u64; N],
    p: u64,
    muls: &mut u64,
) -> u64 {
    let mut acc = 0u64;
    for (m, &c) in comp {
        let mut t = c;
        for k in 0..N {
            for _ in 0..m[k] {
                t = mm(t, e[k], p);
                *muls += 1;
            }
        }
        acc = am(acc, t, p);
    }
    acc
}

/// All `F_p`-rational solutions of four polynomial equations in four
/// unknowns, each of total degree `equation_degree`, through a Macaulay
/// matrix at `macaulay_degree` (the regularity bound `4(d − 1) + 1` for a
/// generic system) and the eigenvectors of the multiplication matrix of
/// the first unknown.  `None` when the degree cut does not close the
/// quotient; every returned solution satisfies all four equations.
pub fn solve_system4(
    comps: &[HashMap<[u8; 4], u64>],
    equation_degree: u8,
    macaulay_degree: u8,
    p: u64,
    rng: &mut StdRng,
    st: &mut QuarticSolve,
    debug: bool,
) -> Option<Vec<[u64; 4]>> {
    solve_system4_traced(
        comps,
        equation_degree,
        macaulay_degree,
        p,
        rng,
        st,
        debug,
        None,
        None,
    )
}

/// [`solve_system4`] in `N` unknowns, each equation shifted by its own
/// degree: the `D₄` system's solve, for callers outside the crate.
pub fn solve_system_n<const N: usize>(
    comps: &[HashMap<[u8; N], u64>],
    equation_degree: u8,
    macaulay_degree: u8,
    p: u64,
    rng: &mut StdRng,
    st: &mut QuarticSolve,
) -> Option<Vec<[u64; N]>> {
    solve_system_n_inner::<N>(
        comps,
        equation_degree,
        macaulay_degree,
        p,
        rng,
        st,
        false,
        None,
        None,
    )
    .ok()
    .flatten()
}

/// The shape of one Macaulay solve, recorded on a first system and replayed
/// on later systems with the same monomial support (F4Remake's idea, for a
/// Macaulay matrix): the leading monomials after triangularisation, the
/// number of rows built, the input rows that reduced to zero, and the number
/// of pivots.  A replay builds only the rows that did not reduce to zero.
#[derive(Clone, Debug, Default, Serialize)]
pub struct ShapeTrace {
    pub leads: Vec<Option<Vec<u8>>>,
    pub rows_built: usize,
    pub zero_ids: Vec<u32>,
    pub pivots: usize,
    pub dim: usize,
}

/// How a traced solve went.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub enum ReplayOutcome {
    /// No trace was given: a full solve, recorded.
    #[default]
    Recorded,
    /// The trace matched: zero rows were not built.
    Replayed,
    /// The trace did not match this system (leads, rows, pivots or
    /// dimension differ): the full solve ran and its cost is in `st`.
    FellBack,
}

/// [`solve_system4`] with a shape trace.  With `replay = Some(t)`, the rows
/// `t` records as reducing to zero are not built; if the replayed echelon
/// does not reach `t`'s pivots and dimension, or the leads or row count
/// differ, the full solve runs instead (charged in `st` as well: a fallback
/// is never free).  With `record = Some(out)`, the full solve's shape is
/// written to `out`.
#[allow(clippy::too_many_arguments)]
pub fn solve_system4_traced(
    comps: &[HashMap<[u8; 4], u64>],
    equation_degree: u8,
    macaulay_degree: u8,
    p: u64,
    rng: &mut StdRng,
    st: &mut QuarticSolve,
    debug: bool,
    replay: Option<&ShapeTrace>,
    record: Option<&mut ShapeTrace>,
) -> Option<Vec<[u64; 4]>> {
    let mut outcome = ReplayOutcome::Recorded;
    if let Some(t) = replay {
        match solve_system4_inner(
            comps,
            equation_degree,
            macaulay_degree,
            p,
            rng,
            st,
            debug,
            Some(t),
            None,
        ) {
            Ok(sols) => {
                st.replay = ReplayOutcome::Replayed;
                return sols;
            }
            Err(()) => outcome = ReplayOutcome::FellBack,
        }
    }
    let r = solve_system4_inner(
        comps,
        equation_degree,
        macaulay_degree,
        p,
        rng,
        st,
        debug,
        None,
        record,
    )
    .unwrap_or(None);
    st.replay = outcome;
    r
}

/// The solve.  `Err(())` only from a replay whose shape did not match.
#[allow(clippy::too_many_arguments)]
fn solve_system4_inner(
    comps: &[HashMap<[u8; 4], u64>],
    equation_degree: u8,
    macaulay_degree: u8,
    p: u64,
    rng: &mut StdRng,
    st: &mut QuarticSolve,
    debug: bool,
    replay: Option<&ShapeTrace>,
    record: Option<&mut ShapeTrace>,
) -> Result<Option<Vec<[u64; 4]>>, ()> {
    solve_system_n_inner::<4>(
        comps,
        equation_degree,
        macaulay_degree,
        p,
        rng,
        st,
        debug,
        replay,
        record,
    )
}

/// The solve in `N` unknowns.  `equation_degree > 0`: every equation is
/// shifted as one of that degree (the `k = 4` systems); `0`: each equation
/// by its own total degree after triangularisation.  `Err(())` only from a
/// replay whose shape did not match.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_system_n_inner<const N: usize>(
    comps: &[HashMap<[u8; N], u64>],
    equation_degree: u8,
    macaulay_degree: u8,
    p: u64,
    rng: &mut StdRng,
    st: &mut QuarticSolve,
    debug: bool,
    replay: Option<&ShapeTrace>,
    record: Option<&mut ShapeTrace>,
) -> Result<Option<Vec<[u64; N]>>, ()> {
    solve_system_w_inner::<N>(
        comps,
        &[1; N],
        false,
        false,
        false,
        equation_degree,
        macaulay_degree,
        p,
        rng,
        st,
        debug,
        replay,
        record,
    )
}

/// [`solve_system_n_inner`] in a weighted grading: unknown `k` has weight
/// `w[k]`, the Macaulay matrix holds the monomials of weighted degree
/// `≤ macaulay_degree`, each equation is shifted by the monomials of
/// weighted degree `≤ macaulay_degree − (its weighted degree)`, and the
/// multiplication matrix is that of the first unknown of weight 1 (so its
/// border stays inside the matrix).  Weights all 1 are the unweighted solve,
/// operation for operation.  A rational eigenvalue whose eigenspace holds
/// more than one solution is split by the multiplication matrices of the
/// other unknowns when their border is inside the matrix; otherwise its
/// solutions are lost, which the caller's exact check sees (they are never
/// returned wrong: every returned solution satisfies all equations).
/// `chains`: the eigen step along the chains of the multiplication
/// unknown ([`solve_along_chains`]) instead of the dense
/// multiplication matrix; `sparsest`: the elimination's pivot in each
/// column is the sparsest candidate row ([`echelon_by_lead_with`]).
/// Both off is the unweighted solve's arithmetic.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_system_w_inner<const N: usize>(
    comps: &[HashMap<[u8; N], u64>],
    w: &[u8; N],
    chains: bool,
    sparsest: bool,
    top_up: bool,
    equation_degree: u8,
    macaulay_degree: u8,
    p: u64,
    rng: &mut StdRng,
    st: &mut QuarticSolve,
    debug: bool,
    replay: Option<&ShapeTrace>,
    record: Option<&mut ShapeTrace>,
) -> Result<Option<Vec<[u64; N]>>, ()> {
    assert!(p < 1 << 15, "rows are u16 and products must fit u32");
    let Some(mvar) = w.iter().position(|&k| k == 1) else {
        st.outcome = "no unknown of weight 1".into();
        return Ok(None);
    };
    let t0 = Instant::now();
    st.equation_degree = equation_degree;
    st.macaulay_degree = macaulay_degree;
    let cols = monomials_w_desc(macaulay_degree, w);
    let col_index: HashMap<[u8; N], usize> =
        cols.iter().enumerate().map(|(i, m)| (*m, i)).collect();
    let mut tri = 0u64;
    let (comps, leads) = triangularise_leads4(comps, w, p, &mut tri);
    st.triangularise_muls += tri;
    let leads_v: Vec<Option<Vec<u8>>> = leads.iter().map(|l| l.map(|m| m.to_vec())).collect();
    if let Some(t) = replay {
        if t.leads != leads_v {
            st.outcome = "replay: leads differ".into();
            return Err(());
        }
    }
    let mut shifts_by: HashMap<u8, Vec<[u8; N]>> = HashMap::new();
    let degree_of = |c: &HashMap<[u8; N], u64>| -> u8 {
        if equation_degree > 0 {
            equation_degree
        } else {
            c.iter()
                .filter(|(_, &v)| v != 0)
                .map(|(m, _)| wdeg(m, w) as u8)
                .max()
                .unwrap_or(0)
        }
    };
    let divides = |lm: &[u8; N], m: &[u8; N]| (0..N).all(|k| lm[k] <= m[k]);
    let mut rows: Vec<Vec<u16>> = Vec::new();
    // With `top_up`, the rows a replay skips, kept for a rank shortfall.
    let mut held: Vec<Vec<u16>> = Vec::new();
    let mut built = 0usize;
    let mut skip_at = 0usize;
    for (i, comp) in comps.iter().enumerate() {
        let dg = degree_of(comp);
        if dg > macaulay_degree {
            st.outcome = "an equation above the Macaulay degree".into();
            return Ok(None);
        }
        let shifts = shifts_by
            .entry(dg)
            .or_insert_with(|| monomials_w_desc::<N>(macaulay_degree - dg, w));
        for sh in shifts.iter() {
            if leads
                .iter()
                .take(i)
                .any(|lm| lm.is_some_and(|lm| divides(&lm, sh)))
            {
                continue;
            }
            let index = built;
            built += 1;
            let mut skip = false;
            if let Some(t) = replay {
                if skip_at < t.zero_ids.len() && t.zero_ids[skip_at] as usize == index {
                    skip_at += 1;
                    if !top_up {
                        continue;
                    }
                    skip = true;
                }
            }
            let mut row = vec![0u16; cols.len()];
            for (&e, &c) in comp {
                let m: [u8; N] = core::array::from_fn(|k| e[k] + sh[k]);
                let ci = col_index[&m];
                row[ci] = am(row[ci] as u64, c, p) as u16;
            }
            if skip {
                held.push(row);
            } else {
                rows.push(row);
            }
        }
    }
    if let Some(t) = replay {
        if t.rows_built != built {
            st.outcome = "replay: row count differs".into();
            return Err(());
        }
    }
    st.rows = rows.len();
    st.cols = cols.len();
    st.build_ms = t0.elapsed().as_secs_f64() * 1e3;
    if debug {
        eprintln!(
            "  macaulay degree {macaulay_degree}: {} rows × {} columns ({:.2} GB), leads {leads:?}",
            st.rows,
            st.cols,
            (st.rows * st.cols * 2) as f64 / 1e9
        );
    }
    let t1 = Instant::now();
    let mut ech_muls = 0u64;
    let mut ech = echelon_by_lead_with(rows, p as u32, &mut ech_muls, debug, sparsest);
    if let Some(t) = replay {
        if top_up {
            let have = ech.pivot_row.iter().filter(|&&r| r != u32::MAX).count();
            if have < t.pivots {
                st.topped_up +=
                    top_up_echelon(&mut ech, held, t.pivots - have, p as u32, &mut ech_muls);
            }
        }
    }
    st.echelon_muls += ech_muls;
    st.echelon_ms = t1.elapsed().as_secs_f64() * 1e3;
    st.zero_rows = ech.zero_rows;
    st.pivots = ech.pivot_row.iter().filter(|&&r| r != u32::MAX).count();
    let is_std: Vec<bool> = (0..cols.len())
        .map(|c| ech.pivot_row[c] == u32::MAX && wdeg(&cols[c], w) < macaulay_degree as u16)
        .collect();
    let standard: Vec<usize> = (0..cols.len()).filter(|&c| is_std[c]).collect();
    let dim = standard.len();
    st.dim = dim;
    // With `top_up` a replay needs no rank check: it stops at the recorded
    // rank, which a specialised system cannot exceed when the recorded one
    // had the generic rank (Macaulay rank is lower semicontinuous), or it
    // has used every row.  Either way it spans what the full matrix spans.
    if let (Some(t), false) = (replay, top_up) {
        if t.pivots != st.pivots || t.dim != dim {
            st.outcome = "replay: pivots or dimension differ".into();
            return Err(());
        }
    }
    if let Some(out) = record {
        *out = ShapeTrace {
            leads: leads_v.clone(),
            rows_built: built,
            zero_ids: ech.zero_ids.clone(),
            pivots: st.pivots,
            dim,
        };
    }
    let mut std_index = vec![u32::MAX; cols.len()];
    for (i, &c) in standard.iter().enumerate() {
        std_index[c] = i as u32;
    }
    if debug {
        eprintln!(
            "  echelon: {} pivots, {} zero rows, dim {dim}, {:.3e} muls, {:.0} s",
            st.pivots,
            st.zero_rows,
            st.echelon_muls as f64,
            st.echelon_ms / 1e3
        );
    }
    // Columns of 1 and of each unknown.
    let mut low: Vec<usize> = Vec::with_capacity(N + 1);
    for v in 0..=N {
        let m: [u8; N] = core::array::from_fn(|k| u8::from(v > 0 && k == v - 1));
        match col_index.get(&m) {
            Some(&c) => low.push(c),
            None => {
                st.outcome = "low-degree monomials missing".into();
                return Ok(None);
            }
        }
    }
    if low.iter().any(|&c| !is_std[c]) {
        st.outcome = "1 or a variable is not standard".into();
        return Ok(None);
    }
    let idx: Vec<usize> = low.iter().map(|&c| std_index[c] as usize).collect();
    let t2 = Instant::now();
    if chains {
        return Ok(solve_along_chains::<N>(
            &ech, &cols, &col_index, &is_std, &standard, &std_index, &idx, mvar, &comps, p, rng,
            st, t0, t2,
        ));
    }
    let mut nfs = NormalForms4 {
        ech: &ech,
        is_std,
        std_index,
        dim,
        p: p as u32,
        memo: vec![None; cols.len()],
        reached: 0,
    };
    // Multiplication matrix of a variable: column b ↦ NF(var · b).
    let mult_matrix = |nfs: &mut NormalForms4,
                       var: usize,
                       muls: &mut u64|
     -> Option<Vec<Vec<u64>>> {
        let mut targets = Vec::with_capacity(dim);
        for &bc in &standard {
            let mut m = cols[bc];
            m[var] += 1;
            targets.push(*col_index.get(&m)?);
        }
        if debug {
            let border = targets.iter().filter(|&&t| !nfs.is_std[t]).count();
            eprintln!("  multiplication by unknown {var}: {border} of {dim} products outside the standard monomials");
        }
        nfs.ensure(&targets, muls).ok()?;
        let mut out = vec![vec![0u64; dim]; dim];
        for (bi, &tc) in targets.iter().enumerate() {
            for (ri, v) in nfs.vector(tc).into_iter().enumerate() {
                out[ri][bi] = v;
            }
        }
        Some(out)
    };
    let mut nf_muls = 0u64;
    let m_e1 = mult_matrix(&mut nfs, mvar, &mut nf_muls);
    st.nf_muls += nf_muls;
    st.nf_ms = t2.elapsed().as_secs_f64() * 1e3;
    st.reached_pivots = nfs.reached;
    let Some(m_e1) = m_e1 else {
        st.outcome = "border unreachable at this degree".into();
        return Ok(None);
    };
    if debug {
        eprintln!(
            "  normal forms: {} pivots reached, {:.3e} muls, {:.0} s",
            st.reached_pivots,
            st.nf_muls as f64,
            st.nf_ms / 1e3
        );
    }
    let t3 = Instant::now();
    let mut cp_muls = 0u64;
    let cp = charpoly_u16(&m_e1, p as u32, &mut cp_muls);
    st.charpoly_muls += cp_muls;
    st.charpoly_ms = t3.elapsed().as_secs_f64() * 1e3;
    let t4 = Instant::now();
    let ring = PolyRing::new(p);
    let lambdas = ring.roots(&cp, rng);
    st.roots_muls += ring.muls.get();
    st.rational_eigenvalues = lambdas.len();
    let mut mt = vec![vec![0u16; dim]; dim];
    for i in 0..dim {
        for j in 0..dim {
            mt[j][i] = m_e1[i][j] as u16;
        }
    }
    drop(m_e1);
    // The other unknowns, in order, for splitting a degenerate eigenspace.
    let other_vars: Vec<usize> = (0..N).filter(|&v| v != mvar).collect();
    let mut others: Vec<Option<Vec<Vec<u64>>>> = vec![None; N - 1];
    let mut out = Vec::new();
    for lam in lambdas {
        let mut ev = 0u64;
        let mut ker = eigenvectors_u16(&mt, lam, p as u32, &mut ev);
        if ker.len() > 1 {
            st.degenerate_eigenspaces += 1;
            for (k, &var) in other_vars.iter().enumerate() {
                if ker.len() <= 1 {
                    break;
                }
                if others[k].is_none() {
                    let mut nm = 0u64;
                    others[k] = mult_matrix(&mut nfs, var, &mut nm);
                    st.nf_muls += nm;
                }
                // An unknown whose border leaves the matrix (weight above
                // 1) cannot split it; with weights all 1 none does.
                let Some(mv) = &others[k] else {
                    continue;
                };
                ker = split_eigenspace(&ker, mv, p, rng, &mut ev);
            }
        }
        st.eigenvector_muls += ev;
        for v in &ker {
            if v[idx[0]] == 0 {
                continue;
            }
            let inv = inv_mod(v[idx[0]], p);
            let e: [u64; N] = core::array::from_fn(|k| mm(v[idx[k + 1]], inv, p));
            if e[mvar] != lam {
                continue;
            }
            let mut vm = 0u64;
            let ok = comps.iter().all(|c| eval4(c, &e, p, &mut vm) == 0);
            st.verify_muls += vm;
            if ok {
                out.push(e);
            }
        }
    }
    st.eigen_ms = t4.elapsed().as_secs_f64() * 1e3;
    out.sort_unstable();
    out.dedup();
    st.e_solutions = out.len();
    st.outcome = "solved".into();
    st.total_ms = t0.elapsed().as_secs_f64() * 1e3;
    st.close();
    Ok(Some(out))
}

// ── The quotient along the chains of the multiplication unknown ───────────
//
// With `x` the multiplication unknown (weight 1) and `B` the standard
// monomials, an order ideal, `B` splits into chains `h, hx, …, hx^{ℓ−1}`
// from the heads `h ∈ B` not divisible by `x`; `M_x` maps each chain
// element to the next, and only a chain's last element leaves `B`, to the
// border monomial `hx^ℓ`.  A left eigenvector `w` of `M_x` for `λ` (the
// evaluation at a solution, up to scale) is then fixed by its values `u`
// on the heads, `w(hx^t) = λ^t u_h`, and the border relations
// `w · NF(h_j x^{ℓ_j}) = λ^{ℓ_j} u_j` are `A(λ) u = 0` for the `k × k`
// matrix `A(λ)_{ji} = Σ_t NF(h_j x^{ℓ_j})[h_i x^t] λ^t − δ_{ji} λ^{ℓ_j}`,
// `k` the number of chains.  `det A(λ) = ± det(λI − M_x)`, of degree
// `Σ ℓ_j = dim`.  So only the `k` border normal forms are computed, the
// characteristic polynomial is `det A` (at every point of `F_p` when
// `p ≤ dim + 1`, whose zeros are then the rational eigenvalues, otherwise
// at `dim + 1` points and interpolated), and an eigenvector is a kernel
// vector of a `k × k` matrix.

/// What the chain step found, or why it could not run (the caller then
/// reports the solve unsolved, its work charged).
struct Chains {
    /// For each standard index: (chain, position).
    pos: Vec<(u32, u32)>,
    /// Chain lengths.
    len: Vec<u32>,
    /// For each chain, its border normal form as (chain, position, coef).
    nf: Vec<Vec<(u32, u32, u64)>>,
}

/// The heads, chains and border normal forms.  `Err` names the failure.
#[allow(clippy::too_many_arguments)]
fn chains_of<const N: usize>(
    ech: &Echelon,
    cols: &[[u8; N]],
    col_index: &HashMap<[u8; N], usize>,
    is_std: &[bool],
    standard: &[usize],
    std_index: &[u32],
    mvar: usize,
    p: u32,
    muls: &mut u64,
) -> Result<Chains, String> {
    let dim = standard.len();
    let ncols = cols.len();
    // Successor of each standard monomial under x, or its border column.
    let mut next: Vec<Result<usize, usize>> = Vec::with_capacity(dim);
    for &bc in standard {
        let mut m = cols[bc];
        m[mvar] += 1;
        let Some(&t) = col_index.get(&m) else {
            return Err("border unreachable at this degree".into());
        };
        next.push(if is_std[t] {
            Ok(std_index[t] as usize)
        } else {
            Err(t)
        });
    }
    let mut pos = vec![(u32::MAX, u32::MAX); dim];
    let mut len = Vec::new();
    let mut ends = Vec::new();
    for s in 0..dim {
        if cols[standard[s]][mvar] != 0 {
            continue;
        }
        let j = len.len() as u32;
        let (mut cur, mut t) = (s, 0u32);
        loop {
            if pos[cur].0 != u32::MAX {
                return Err("chain structure: a standard monomial on two chains".into());
            }
            pos[cur] = (j, t);
            t += 1;
            match next[cur] {
                Ok(n) => cur = n,
                Err(b) => {
                    ends.push(b);
                    break;
                }
            }
        }
        len.push(t);
    }
    if pos.iter().any(|&(j, _)| j == u32::MAX) {
        return Err("chain structure: standard monomials are not an order ideal in x".into());
    }
    // nnz of each pivot row from its lead: bookkeeping, not arithmetic.
    let nnz: Vec<u64> = (0..ncols)
        .map(|c| match ech.pivot_row[c] {
            u32::MAX => 0,
            r => ech.rows[r as usize][c..]
                .iter()
                .filter(|&&v| v != 0)
                .count() as u64,
        })
        .collect();
    // Each border normal form by reducing e_b left to right by the pivot
    // rows (normalised, lead 1): one multiplication per nonzero entry of a
    // pivot row per use.
    let reduced: Vec<Result<(Vec<u16>, u64), String>> = ends
        .par_iter()
        .map(|&b| {
            let mut v = vec![0u16; ncols];
            v[b] = 1;
            let mut m = 0u64;
            for c in b..ncols {
                if v[c] == 0 || is_std[c] {
                    continue;
                }
                let pr = ech.pivot_row[c];
                if pr == u32::MAX {
                    return Err("border unreachable at this degree".into());
                }
                let f = v[c] as u32;
                axpy(&mut v[c..], &ech.rows[pr as usize][c..], f, p);
                m += nnz[c];
            }
            Ok((v, m))
        })
        .collect();
    let mut nf = Vec::with_capacity(ends.len());
    for r in reduced {
        let (v, m) = r?;
        *muls += m;
        nf.push(
            standard
                .iter()
                .enumerate()
                .filter(|&(_, &c)| v[c] != 0)
                .map(|(s, &c)| (pos[s].0, pos[s].1, v[c] as u64))
                .collect(),
        );
    }
    Ok(Chains { pos, len, nf })
}

impl Chains {
    fn k(&self) -> usize {
        self.len.len()
    }

    /// `A(λ)`, charged one multiplication per term and per power.
    fn matrix(&self, lam: u64, p: u64, muls: &mut u64) -> Vec<Vec<u64>> {
        let top = *self.len.iter().max().unwrap_or(&0) as usize;
        let mut pw = vec![1u64; top + 1];
        for t in 1..=top {
            pw[t] = mm(pw[t - 1], lam, p);
        }
        *muls += top as u64;
        let k = self.k();
        let mut a = vec![vec![0u64; k]; k];
        for (j, row) in self.nf.iter().enumerate() {
            for &(i, t, c) in row {
                let e = &mut a[j][i as usize];
                *e = am(*e, mm(c, pw[t as usize], p), p);
            }
            *muls += row.len() as u64;
            let d = &mut a[j][j];
            *d = sm(*d, pw[self.len[j] as usize], p);
        }
        a
    }

    /// `det A(λ)` by elimination, charged as `rref_mod_p` charges.
    fn det(&self, lam: u64, p: u64, muls: &mut u64) -> u64 {
        let mut a = self.matrix(lam, p, muls);
        let k = a.len();
        let mut det = 1u64;
        for c in 0..k {
            let Some(r) = (c..k).find(|&r| a[r][c] != 0) else {
                return 0;
            };
            if r != c {
                a.swap(r, c);
                det = (p - det) % p;
            }
            det = mm(det, a[c][c], p);
            let inv = inv_mod(a[c][c], p);
            *muls += 1;
            for r in (c + 1)..k {
                if a[r][c] == 0 {
                    continue;
                }
                let f = mm(a[r][c], inv, p);
                for j in (c + 1)..k {
                    a[r][j] = sm(a[r][j], mm(f, a[c][j], p), p);
                }
                *muls += (k - c) as u64;
            }
        }
        det
    }

    /// The left eigenvectors of `M_x` for `λ` (as functionals on the
    /// standard monomials): the kernel of `A(λ)`, along the chains.
    fn eigenvectors(&self, lam: u64, p: u64, muls: &mut u64) -> Vec<Vec<u64>> {
        let a = self.matrix(lam, p, muls);
        let ker = eigenvectors_mod_p(&a, 0, p, muls);
        let top = *self.len.iter().max().unwrap_or(&0) as usize;
        let mut pw = vec![1u64; top + 1];
        for t in 1..=top {
            pw[t] = mm(pw[t - 1], lam, p);
        }
        *muls += top as u64;
        ker.iter()
            .map(|u| {
                *muls += self.pos.len() as u64;
                self.pos
                    .iter()
                    .map(|&(j, t)| mm(u[j as usize], pw[t as usize], p))
                    .collect()
            })
            .collect()
    }
}

/// The rational roots of `det A`: its zeros on `F_p` when `p ≤ dim + 1`,
/// else the roots of its interpolation on `0, …, dim`.  `Err` if the
/// interpolated determinant does not have degree `dim` and leading
/// coefficient `±1`.
fn chain_eigenvalues(
    ch: &Chains,
    dim: usize,
    p: u64,
    rng: &mut StdRng,
    cp_muls: &mut u64,
    roots_muls: &mut u64,
) -> Result<Vec<u64>, String> {
    let npts = if p <= dim as u64 + 1 {
        p
    } else {
        dim as u64 + 1
    };
    let vals: Vec<(u64, u64)> = (0..npts)
        .into_par_iter()
        .map(|x| {
            let mut m = 0u64;
            (ch.det(x, p, &mut m), m)
        })
        .collect();
    *cp_muls += vals.iter().map(|v| v.1).sum::<u64>();
    if npts == p {
        return Ok((0..p).filter(|&x| vals[x as usize].0 == 0).collect());
    }
    // Newton's divided differences on 0, …, dim (inverses of 1, …, dim
    // once), then the Newton form expanded: O(dim²) each.
    let n = dim + 1;
    let inv: Vec<u64> = (0..n as u64)
        .map(|d| if d == 0 { 0 } else { inv_mod(d, p) })
        .collect();
    let mut dd: Vec<u64> = vals.iter().map(|v| v.0).collect();
    for lvl in 1..n {
        for i in (lvl..n).rev() {
            dd[i] = mm(sm(dd[i], dd[i - 1], p), inv[lvl], p);
        }
        *cp_muls += (n - lvl) as u64;
    }
    // Σ dd[i] Π_{j<i} (X − j), Horner from the top.
    let mut poly = vec![0u64; n];
    for i in (0..n).rev() {
        // poly ← poly · (X − i) + dd[i]
        let mut nxt = vec![0u64; n];
        for d in 0..n - 1 {
            if poly[d] != 0 {
                nxt[d + 1] = am(nxt[d + 1], poly[d], p);
                nxt[d] = sm(nxt[d], mm(poly[d], i as u64 % p, p), p);
            }
        }
        *cp_muls += (n - 1 - i) as u64;
        nxt[0] = am(nxt[0], dd[i], p);
        poly = nxt;
    }
    let lead = poly[dim];
    if lead != 1 && lead != p - 1 {
        return Err(format!("det A(λ) has leading coefficient {lead}, not ±1"));
    }
    let ring = PolyRing::new(p);
    let r = ring.roots(&UPoly(poly), rng);
    *roots_muls += ring.muls.get();
    Ok(r)
}

/// [`split_eigenspace`] for a multiplication unknown whose border may
/// leave the matrix: `wᵢ · M_v` is known only at the standard monomials
/// `b` with `v·b` standard (`known`: `(b, v·b)` as standard indices), which
/// determine `a` when the `w` restricted to them are independent.
/// Returns `w` unchanged if they are not (the eigenspace stays unsplit).
fn split_eigenspace_partial(
    w: &[Vec<u64>],
    known: &[(usize, usize)],
    p: u64,
    rng: &mut StdRng,
    muls: &mut u64,
) -> Vec<Vec<u64>> {
    let k = w.len();
    let mut a = vec![vec![0u64; k]; k];
    for i in 0..k {
        let mut sys: Vec<Vec<u64>> = known
            .iter()
            .map(|&(l, t)| {
                let mut row: Vec<u64> = (0..k).map(|j| w[j][l]).collect();
                row.push(w[i][t]);
                row
            })
            .collect();
        let pivots = rref_mod_p(&mut sys, p, muls);
        if pivots.len() != k || pivots.iter().any(|&c| c >= k) {
            return w.to_vec();
        }
        for (r, &pc) in pivots.iter().enumerate() {
            a[i][pc] = sys[r][k];
        }
    }
    let mut at = vec![vec![0u64; k]; k];
    for i in 0..k {
        for j in 0..k {
            at[j][i] = a[i][j];
        }
    }
    let cp = charpoly_mod_p(&a, p, muls);
    let ring = PolyRing::new(p);
    let mus = ring.roots(&cp, rng);
    *muls += ring.muls.get();
    let dim = w[0].len();
    let mut out = Vec::new();
    for mu in mus {
        for v in eigenvectors_mod_p(&at, mu, p, muls) {
            let mut u = vec![0u64; dim];
            for j in 0..k {
                if v[j] == 0 {
                    continue;
                }
                for r in 0..dim {
                    if w[j][r] != 0 {
                        u[r] = am(u[r], mm(v[j], w[j][r], p), p);
                        *muls += 1;
                    }
                }
            }
            out.push(u);
        }
    }
    out
}

/// The eigen step of [`solve_system_w_inner`] along the chains: border
/// normal forms (`nf_muls`), `det A` at points (`charpoly_muls`), its
/// roots (`roots_muls`), kernels of `A(λ)` and splits (`eigenvector_muls`),
/// each solution checked in every equation (`verify_muls`).  A degenerate
/// eigenspace is split by the other unknowns where their products stay
/// standard.  `None` (with `st.outcome`) when the chain step cannot run.
#[allow(clippy::too_many_arguments)]
fn solve_along_chains<const N: usize>(
    ech: &Echelon,
    cols: &[[u8; N]],
    col_index: &HashMap<[u8; N], usize>,
    is_std: &[bool],
    standard: &[usize],
    std_index: &[u32],
    idx: &[usize],
    mvar: usize,
    comps: &[HashMap<[u8; N], u64>],
    p: u64,
    rng: &mut StdRng,
    st: &mut QuarticSolve,
    t0: Instant,
    t2: Instant,
) -> Option<Vec<[u64; N]>> {
    let dim = standard.len();
    let mut nf_muls = 0u64;
    let ch = chains_of::<N>(
        ech,
        cols,
        col_index,
        is_std,
        standard,
        std_index,
        mvar,
        p as u32,
        &mut nf_muls,
    );
    st.nf_muls += nf_muls;
    st.nf_ms = t2.elapsed().as_secs_f64() * 1e3;
    let ch = match ch {
        Ok(ch) => ch,
        Err(why) => {
            st.outcome = why;
            return None;
        }
    };
    st.reached_pivots = ch.k();
    let t3 = Instant::now();
    let (mut cp, mut rm) = (0u64, 0u64);
    let lambdas = chain_eigenvalues(&ch, dim, p, rng, &mut cp, &mut rm);
    st.charpoly_muls += cp;
    st.roots_muls += rm;
    st.charpoly_ms = t3.elapsed().as_secs_f64() * 1e3;
    let lambdas = match lambdas {
        Ok(l) => l,
        Err(why) => {
            st.outcome = why;
            return None;
        }
    };
    st.rational_eigenvalues = lambdas.len();
    let t4 = Instant::now();
    // For each other unknown, the standard monomials it keeps standard.
    let known: Vec<Vec<(usize, usize)>> = (0..N)
        .filter(|&v| v != mvar)
        .map(|v| {
            (0..dim)
                .filter_map(|s| {
                    let mut m = cols[standard[s]];
                    m[v] += 1;
                    let &t = col_index.get(&m)?;
                    is_std[t].then(|| (s, std_index[t] as usize))
                })
                .collect()
        })
        .collect();
    let mut out = Vec::new();
    for lam in lambdas {
        let mut ev = 0u64;
        let mut ker = ch.eigenvectors(lam, p, &mut ev);
        if ker.len() > 1 {
            st.degenerate_eigenspaces += 1;
            for kn in &known {
                if ker.len() <= 1 {
                    break;
                }
                ker = split_eigenspace_partial(&ker, kn, p, rng, &mut ev);
            }
        }
        st.eigenvector_muls += ev;
        for v in &ker {
            if v[idx[0]] == 0 {
                continue;
            }
            let inv = inv_mod(v[idx[0]], p);
            let e: [u64; N] = core::array::from_fn(|k| mm(v[idx[k + 1]], inv, p));
            if e[mvar] != lam {
                continue;
            }
            let mut vm = 0u64;
            let ok = comps.iter().all(|c| eval4(c, &e, p, &mut vm) == 0);
            st.verify_muls += vm;
            if ok {
                out.push(e);
            }
        }
    }
    st.eigen_ms = t4.elapsed().as_secs_f64() * 1e3;
    out.sort_unstable();
    out.dedup();
    st.e_solutions = out.len();
    st.outcome = "solved".into();
    st.total_ms = t0.elapsed().as_secs_f64() * 1e3;
    st.close();
    Some(out)
}

/// The multisets `{x₁, …, x₄} ⊂ F_p` whose elementary symmetric functions
/// are `e` — the roots of `T⁴ − e₁T³ + e₂T² − e₃T + e₄`, when it splits.
pub(crate) fn split_quartic(
    e: &[u64; 4],
    p: u64,
    rng: &mut StdRng,
    muls: &mut u64,
) -> Option<[u64; 4]> {
    let ring = PolyRing::new(p);
    let f = UPoly(vec![e[3], (p - e[2]) % p, e[1], (p - e[0]) % p, 1]);
    let roots = ring.roots(&f, rng);
    *muls += ring.muls.get();
    for &r1 in &roots {
        for &r2 in &roots {
            for &r3 in &roots {
                for &r4 in &roots {
                    if !(r1 <= r2 && r2 <= r3 && r3 <= r4) {
                        continue;
                    }
                    let x = [r1, r2, r3, r4];
                    let mut sy = [0u64; 5];
                    sy[0] = 1;
                    for &v in &x {
                        for k in (1..5).rev() {
                            sy[k] = am(sy[k], mm(sy[k - 1], v, p), p);
                        }
                    }
                    if sy[1] == e[0] && sy[2] == e[1] && sy[3] == e[2] && sy[4] == e[3] {
                        return Some(x);
                    }
                }
            }
        }
    }
    None
}

/// Gaudry's solve at `k = 4`: the abscissae `x₁ ≤ … ≤ x₄` in `F_p` with
/// `S₅(x₁, …, x₄, x_R) = 0`.  `C₄` is `st.total_muls`.
pub fn solve_s5_subspace(
    curve: &Curve4,
    pre: &SymmetrisedS5,
    x_r: &E4,
    rng: &mut StdRng,
    st: &mut QuarticSolve,
    debug: bool,
) -> Option<Vec<[u64; 4]>> {
    let f = &curve.f;
    let before = f.muls();
    let comps = pre.weil_restrict(f, x_r);
    st.weil_muls += f.muls() - before;
    let es = solve_system4(&comps, 8, 29, f.p, rng, st, debug)?;
    let mut out = Vec::new();
    for e in es {
        let mut sm_ = 0u64;
        if let Some(x) = split_quartic(&e, f.p, rng, &mut sm_) {
            out.push(x);
        }
        st.split_muls += sm_;
    }
    st.close();
    out.sort_unstable();
    Some(out)
}

// ── The independent check ─────────────────────────────────────────────────

/// The subspace factor base: one point per `x ∈ F_p` whose right-hand side
/// is a square in `F_{p⁴}`.
pub fn factor_base(curve: &Curve4, rng: &mut StdRng) -> Vec<Pt4> {
    (0..curve.f.p)
        .filter_map(|x| curve.lift_x(&curve.f.from_fp(x), rng))
        .collect()
}

/// Every set of four distinct base abscissae `{x_i, x_j, x_k, x_l}` with
/// `±P_i ± P_j ± P_k ± P_l = R`, by meet in the middle over pairs:
/// `x(R − s_k P_k − s_l P_l)` looked up among the `x(P_i ± P_j)`.
pub fn mitm_decompositions(curve: &Curve4, base: &[Pt4], r: &Pt4) -> BTreeSet<[u64; 4]> {
    let n = base.len();
    let mut table: HashMap<E4, Vec<(usize, usize)>> = HashMap::new();
    for i in 0..n {
        for j in (i + 1)..n {
            for q in [curve.add(&base[i], &base[j]), curve.sub(&base[i], &base[j])] {
                if !q.inf {
                    table.entry(q.x).or_default().push((i, j));
                }
            }
        }
    }
    let mut out = BTreeSet::new();
    for k in 0..n {
        for l in (k + 1)..n {
            for (sk, sl) in [(1, 1), (1, -1), (-1, 1), (-1, -1)] {
                let pk = if sk == 1 {
                    base[k]
                } else {
                    curve.neg(&base[k])
                };
                let pl = if sl == 1 {
                    base[l]
                } else {
                    curve.neg(&base[l])
                };
                let v = curve.sub(&curve.sub(r, &pk), &pl);
                if v.inf {
                    continue;
                }
                if let Some(pairs) = table.get(&v.x) {
                    for &(i, j) in pairs {
                        if i == k || i == l || j == k || j == l {
                            continue;
                        }
                        let mut xs = [
                            base[i].x.0[0],
                            base[j].x.0[0],
                            base[k].x.0[0],
                            base[l].x.0[0],
                        ];
                        xs.sort_unstable();
                        out.insert(xs);
                    }
                }
            }
        }
    }
    out
}

/// The solver's decompositions in the oracle's terms: four distinct
/// abscissae, all in the base.
pub fn base_quadruples(sols: &[[u64; 4]], base: &[Pt4]) -> BTreeSet<[u64; 4]> {
    let xs: std::collections::HashSet<u64> = base.iter().map(|q| q.x.0[0]).collect();
    sols.iter()
        .filter(|x| x[0] < x[1] && x[1] < x[2] && x[2] < x[3])
        .filter(|x| x.iter().all(|v| xs.contains(v)))
        .copied()
        .collect()
}

// ── r∞: the linear algebra and rho at k = 4, measured (§11.18) ───────────

/// A prime-order curve over `F_{p⁴}` with a planted logarithm `Q = [d]G`.
pub struct Instance4 {
    pub curve: Curve4,
    pub n: u64,
    pub g: Pt4,
    pub d: u64,
    pub q: Pt4,
}

/// The same instance as plain data, shareable across threads.
#[derive(Clone, Copy)]
pub struct Spec4 {
    pub p: u64,
    pub a: E4,
    pub b: E4,
    pub n: u64,
    pub g: Pt4,
    pub d: u64,
    pub q: Pt4,
}

impl Instance4 {
    pub fn spec(&self) -> Spec4 {
        Spec4 {
            p: self.curve.f.p,
            a: self.curve.a,
            b: self.curve.b,
            n: self.n,
            g: self.g,
            d: self.d,
            q: self.q,
        }
    }
}

fn random_point(curve: &Curve4, rng: &mut StdRng) -> Pt4 {
    loop {
        let x = curve.f.random(rng);
        if let Some(pt) = curve.lift_x(&x, rng) {
            return pt;
        }
    }
}

/// `#E(F_{p⁴})` when it is prime: baby-step giant-step for the multiple of a
/// random point's order inside the Hasse interval `q + 1 ± 2p²`, accepted
/// only when that multiple is prime (then it is the point's order, and the
/// group's, since the interval is shorter than the order).
fn prime_group_order(curve: &Curve4, rng: &mut StdRng) -> Option<u64> {
    let p = curve.f.p;
    let q = p.pow(4);
    let lo = q + 1 - 2 * p * p;
    let width = 4 * p * p;
    let pt = random_point(curve, rng);
    let steps = (width as f64).sqrt() as u64 + 1;
    let mut table: HashMap<(E4, E4), u64> = HashMap::new();
    let mut jp = Pt4::INF;
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
            return super::residual_walk::is_prime_u64(m).then_some(m);
        }
        t = curve.add(&t, &giant);
        i += 1;
    }
    None
}

pub fn generate_instance4(p: u64, seed: u64) -> Instance4 {
    let mut rng = <StdRng as rand::SeedableRng>::seed_from_u64(seed ^ 0x4A11);
    loop {
        let curve = Curve4::random(p, &mut rng);
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
        return Instance4 { curve, n, g, d, q };
    }
}

/// `x(P_i + s·P_j) ↦ (i, j, s, y)` over the base, built once per curve.
pub struct PairTable {
    pub(crate) map: HashMap<E4, Vec<(u32, u32, i8, E4)>>,
}

pub fn pair_table(curve: &Curve4, base: &[Pt4]) -> PairTable {
    let mut map: HashMap<E4, Vec<(u32, u32, i8, E4)>> = HashMap::new();
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
    PairTable { map }
}

/// Every decomposition `R = Σ s_t P_{i_t}` over four distinct base points,
/// with its signs, by meet in the middle against the pair table.
pub fn mitm_signed(
    curve: &Curve4,
    base: &[Pt4],
    table: &PairTable,
    r: &Pt4,
) -> Vec<Vec<(usize, i64)>> {
    let n = base.len();
    let mut out: BTreeSet<Vec<(usize, i64)>> = BTreeSet::new();
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
                        let mut terms = vec![(k, sk), (l, sl), (i, t), (j, t * s as i64)];
                        terms.sort_unstable();
                        out.insert(terms);
                    }
                }
            }
        }
    }
    out.into_iter().collect()
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct RhoRun4 {
    pub seed: u64,
    pub steps: u64,
    pub group_ops: u64,
    pub s: f64,
    pub correct: bool,
}

/// r-adding Pollard rho with exhaustive storage, as `gaudry_cubic::run_rho3`
/// runs it, charging every group operation including the walk's setup.
pub fn rho4(inst: &Spec4, seed: u64) -> RhoRun4 {
    let curve = Curve4::with_params(inst.p, inst.a, inst.b);
    let n = inst.n;
    let mut rng = <StdRng as rand::SeedableRng>::seed_from_u64(seed ^ 0x8D04);
    let r = 32usize;
    let mults: Vec<(u64, u64, Pt4)> = (0..r)
        .map(|_| {
            let al = rng.gen_range(0..n);
            let be = rng.gen_range(1..n);
            (
                al,
                be,
                curve.add(&curve.mul(&inst.g, al), &curve.mul(&inst.q, be)),
            )
        })
        .collect();
    let mut table: HashMap<Pt4, (u64, u64)> = HashMap::new();
    let mut steps = 0u64;
    let mut found = None;
    'outer: loop {
        let mut a = rng.gen_range(0..n);
        let mut b = rng.gen_range(1..n);
        let mut l = curve.add(&curve.mul(&inst.g, a), &curve.mul(&inst.q, b));
        loop {
            steps += 1;
            if let Some(&(a2, b2)) = table.get(&l) {
                let db = (b + n - b2) % n;
                if db != 0 {
                    let da = (a2 + n - a) % n;
                    found = Some(mm(da, inv_mod(db, n), n));
                    break 'outer;
                }
                break;
            }
            table.insert(l, (a, b));
            let h = l.x.0[0] ^ l.x.0[1].wrapping_mul(0x9E37_79B9) ^ l.x.0[2].rotate_left(17);
            let j = (super::residual_walk::mix64(h) % r as u64) as usize;
            l = curve.add(&l, &mults[j].2);
            a = (a + mults[j].0) % n;
            b = (b + mults[j].1) % n;
            if steps > 64 * (n as f64).sqrt() as u64 + 1_000_000 {
                break 'outer;
            }
        }
    }
    let ops = curve.ops();
    RhoRun4 {
        seed,
        steps,
        group_ops: ops,
        s: ops as f64 / (n as f64).sqrt(),
        correct: found == Some(inst.d),
    }
}

/// One `k = 4` run of §11.18: relations from [`mitm_signed`], filtered and
/// solved as §11.7 does, the logarithm verified, and rho on the same group.
#[derive(Clone, Debug, Default, Serialize)]
pub struct K4LaReport {
    pub p: u64,
    pub seed: u64,
    pub n: u64,
    pub bits: f64,
    pub base: usize,
    pub residuals: u64,
    pub decompositions: u64,
    pub decomposition_rate: f64,
    pub relations: usize,
    /// Unknowns in the square core the successful solve used.
    pub unknowns: usize,
    pub filtered_out: usize,
    pub phi: f64,
    pub row_weight: f64,
    pub la_attempts: u64,
    /// Multiplications mod `n`, every attempt.
    pub la_ops: u64,
    pub la_ops_last: u64,
    /// `la_ops / unknowns²`, against §11.16's `20`.
    pub wiedemann_constant: f64,
    pub solved: bool,
    pub correct: bool,
    /// Group operations the relation phase spent (residuals and the oracle),
    /// kept out of `r`.
    pub relation_group_ops: u64,
    pub fp_muls_per_add: f64,
    pub rho: Vec<RhoRun4>,
    pub rho_s_mean: f64,
    pub rho_s_sd: f64,
    /// `LA / rho` in `F_p` multiplications: `la_ops·16 / (S_rho·√n·c_add)`.
    pub r: f64,
    pub wall_ms: f64,
}

pub fn run_k4_la(p: u64, seed: u64, rho_runs: usize) -> K4LaReport {
    use super::gaudry_cubic::{square_core, wiedemann_u64, SparseRel};
    let start = Instant::now();
    let inst = generate_instance4(p, seed);
    let curve = &inst.curve;
    let n = inst.n;
    let mut rng = <StdRng as rand::SeedableRng>::seed_from_u64(seed ^ 0x1A4);
    let base = factor_base(curve, &mut rng);
    let fp_per_add = {
        curve.f.reset_muls();
        let mut acc = base[0];
        for _ in 0..64 {
            acc = curve.add(&acc, &base[1]);
        }
        curve.f.muls() as f64 / 64.0
    };
    curve.reset_ops();
    let table = pair_table(curve, &base);
    let spec = inst.spec();
    let small = base.len();
    let unknowns = small + 1;
    let mut rep = K4LaReport {
        p,
        seed,
        n,
        bits: (n as f64).log2(),
        base: small,
        fp_muls_per_add: fp_per_add,
        ..Default::default()
    };
    let mut rel_ops = curve.ops();
    let mut full_rels: Vec<SparseRel> = Vec::new();
    let mut next_attempt = unknowns;
    let mut la_ops = 0u64;
    'collect: while rep.residuals < 50_000_000 {
        let batch: Vec<(u64, u64)> = (0..64)
            .map(|_| (rng.gen_range(0..n), rng.gen_range(1..n)))
            .collect();
        let found: Vec<(u64, u64, Pt4, Vec<Vec<(usize, i64)>>, u64)> = batch
            .par_iter()
            .map(|&(a, b)| {
                let c = Curve4::with_params(spec.p, spec.a, spec.b);
                let r = c.add(&c.mul(&spec.g, a), &c.mul(&spec.q, b));
                let decs = mitm_signed(&c, &base, &table, &r);
                (a, b, r, decs, c.ops())
            })
            .collect();
        for (a, b, r, decs, ops) in found {
            rep.residuals += 1;
            rel_ops += ops;
            for terms in decs {
                // Verify aG + bQ = Σ s P before using it.
                let mut acc = r;
                for &(i, s) in &terms {
                    acc = if s == 1 {
                        curve.sub(&acc, &base[i])
                    } else {
                        curve.add(&acc, &base[i])
                    };
                }
                if !acc.inf {
                    continue;
                }
                rep.decompositions += 1;
                let mut cols: Vec<(usize, u64)> = terms
                    .iter()
                    .map(|&(i, s)| (i, if s == 1 { 1 } else { n - 1 }))
                    .collect();
                cols.push((small, (n - b) % n));
                cols.sort_unstable();
                full_rels.push(SparseRel { cols, rhs: a });
            }
            if full_rels.len() >= next_attempt {
                if let Some(core) = square_core(&full_rels, small) {
                    rep.la_attempts += 1;
                    let mut map = vec![usize::MAX; unknowns];
                    let mut k = 0;
                    for c in 0..unknowns {
                        if core.columns.get(c).copied().unwrap_or(false) {
                            map[c] = k;
                            k += 1;
                        }
                    }
                    let sel: Vec<SparseRel> = core
                        .rows
                        .iter()
                        .map(|&i| SparseRel {
                            cols: full_rels[i]
                                .cols
                                .iter()
                                .map(|&(c, v)| (map[c], v))
                                .collect(),
                            rhs: full_rels[i].rhs,
                        })
                        .collect();
                    let before = la_ops;
                    let x = wiedemann_u64(&sel, sel.len(), n, &mut rng, &mut la_ops);
                    let dd = x.map(|x| x[map[small]]);
                    if dd.is_some_and(|dd| curve.mul(&inst.g, dd) == inst.q) {
                        rep.solved = true;
                        rep.correct = dd == Some(inst.d);
                        rep.unknowns = sel.len();
                        rep.filtered_out = full_rels.len() - sel.len();
                        rep.la_ops_last = la_ops - before;
                        break 'collect;
                    }
                    next_attempt = full_rels.len() + (sel.len() / 20).max(1);
                } else {
                    next_attempt = full_rels.len() + (unknowns / 20).max(1);
                }
            }
        }
    }
    rep.relations = full_rels.len();
    rep.relation_group_ops = rel_ops;
    rep.la_ops = la_ops;
    rep.decomposition_rate = rep.decompositions as f64 / rep.residuals.max(1) as f64;
    rep.phi = rep.unknowns as f64 / small as f64;
    rep.row_weight = full_rels.iter().map(|r| r.cols.len()).sum::<usize>() as f64
        / full_rels.len().max(1) as f64;
    rep.wiedemann_constant = la_ops as f64 / (rep.unknowns.max(1) as f64).powi(2);
    rep.rho = (0..rho_runs as u64)
        .into_par_iter()
        .map(|k| rho4(&spec, seed * 1000 + k))
        .collect();
    let ss: Vec<f64> = rep.rho.iter().map(|r| r.s).collect();
    rep.rho_s_mean = ss.iter().sum::<f64>() / ss.len().max(1) as f64;
    rep.rho_s_sd = (ss.iter().map(|s| (s - rep.rho_s_mean).powi(2)).sum::<f64>()
        / (ss.len().max(2) - 1) as f64)
        .sqrt();
    rep.r = la_ops as f64 * 16.0 / (rep.rho_s_mean * (n as f64).sqrt() * fp_per_add);
    rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
    rep
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn u16_kernels_match_the_reference_routines() {
        use super::super::gaudry_cubic::{charpoly_mod_p, eigenvectors_mod_p};
        let mut rng = StdRng::seed_from_u64(7);
        for (p, n, zeros) in [
            (269u64, 40usize, 0.0),
            (269, 33, 0.6),
            (257, 25, 0.9),
            (13, 30, 0.3),
        ] {
            let mut m: Vec<Vec<u64>> = (0..n)
                .map(|_| {
                    (0..n)
                        .map(|_| {
                            if rng.gen::<f64>() < zeros {
                                0
                            } else {
                                rng.gen_range(0..p)
                            }
                        })
                        .collect()
                })
                .collect();
            // A planted eigenvalue with a two-dimensional eigenspace.
            for row in m.iter_mut().take(2) {
                row.iter_mut().for_each(|v| *v = 0);
            }
            m[0][0] = 5 % p;
            m[1][1] = 5 % p;
            let (mut c_ref, mut c_new) = (0u64, 0u64);
            let cp_ref = charpoly_mod_p(&m, p, &mut c_ref);
            let cp_new = charpoly_u16(&m, p as u32, &mut c_new);
            assert_eq!(cp_ref.0, cp_new.0, "p {p} n {n}");
            assert_eq!(c_ref, c_new, "charpoly count, p {p} n {n}");
            let m16: Vec<Vec<u16>> = m
                .iter()
                .map(|r| r.iter().map(|&v| v as u16).collect())
                .collect();
            for lam in [5 % p, 0, 1] {
                let (mut e_ref, mut e_new) = (0u64, 0u64);
                assert_eq!(
                    eigenvectors_mod_p(&m, lam, p, &mut e_ref),
                    eigenvectors_u16(&m16, lam, p as u32, &mut e_new)
                );
                assert_eq!(e_ref, e_new, "eigenvector count, p {p} n {n} λ {lam}");
            }
        }
    }

    #[test]
    fn fp4_is_a_field() {
        let f = Fp4::new(269);
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..50 {
            let a = f.random(&mut rng);
            let b = f.random(&mut rng);
            assert_eq!(f.mul(&a, &b), f.mul(&b, &a));
            if !a.is_zero() {
                assert_eq!(f.mul(&a, &f.inv(&a)), E4::ONE);
            }
            let s = f.sq(&a);
            let r = f.sqrt(&s, &mut rng).expect("a square has a root");
            assert_eq!(f.sq(&r), s);
        }
        // The inversion is charged as §11.16 derives it: 40.
        f.reset_muls();
        f.inv(&E4([3, 1, 4, 1]));
        assert_eq!(f.muls(), 40);
    }

    #[test]
    fn s5_vanishes_on_four_point_sums_and_is_symmetric() {
        let mut rng = StdRng::seed_from_u64(2);
        let curve = Curve4::random(269, &mut rng);
        let f = &curve.f;
        let mut pts: Vec<Pt4> = Vec::new();
        while pts.len() < 4 {
            let x = f.random(&mut rng);
            if let Some(q) = curve.lift_x(&x, &mut rng) {
                pts.push(q);
            }
        }
        let sum = pts.iter().fold(Pt4::INF, |acc, q| curve.add(&acc, q));
        assert!(curve.on_curve(&sum));
        let x = [pts[0].x, pts[1].x, pts[2].x, pts[3].x, sum.x];
        assert!(curve.s5(&x).is_zero(), "S5 must vanish on a four-point sum");
        // Signs do not matter: x(P1 − P2 + P3 + P4) is a root too.
        let mixed = curve.add(&curve.add(&curve.sub(&pts[0], &pts[1]), &pts[2]), &pts[3]);
        let xm = [pts[0].x, pts[1].x, pts[2].x, pts[3].x, mixed.x];
        assert!(curve.s5(&xm).is_zero());
        // Symmetric in the first four arguments at a random point.
        let r: [E4; 5] = core::array::from_fn(|_| f.random(&mut rng));
        let v = curve.s5(&r);
        assert!(!v.is_zero());
        for perm in [[1, 0, 2, 3], [2, 3, 0, 1], [3, 1, 2, 0], [0, 2, 1, 3]] {
            let q = [r[perm[0]], r[perm[1]], r[perm[2]], r[perm[3]], r[4]];
            assert_eq!(curve.s5(&q), v, "not symmetric under {perm:?}");
        }
    }

    /// Four random equations of total degree `d` in four unknowns with a
    /// planted `F_p`-rational root: the solver must find it, every answer
    /// must satisfy all four, and the quotient must have the Bézout
    /// dimension `d⁴` at the Macaulay bound `4(d − 1) + 1`.
    fn planted(d: u8, seed: u64) {
        let p = 269u64;
        let mut rng = StdRng::seed_from_u64(seed);
        let root: [u64; 4] = core::array::from_fn(|_| rng.gen_range(0..p));
        let monos = exponents4(d);
        let comps: Vec<HashMap<[u8; 4], u64>> = (0..4)
            .map(|_| {
                let mut c: HashMap<[u8; 4], u64> =
                    monos.iter().map(|m| (*m, rng.gen_range(0..p))).collect();
                let mut dummy = 0u64;
                let v = eval4(&c, &root, p, &mut dummy);
                let e = c.entry([0, 0, 0, 0]).or_insert(0);
                *e = sm(*e, v, p);
                c
            })
            .collect();
        let mut st = QuarticSolve::default();
        let sols = solve_system4(&comps, d, 4 * (d - 1) + 1, p, &mut rng, &mut st, false)
            .expect("a generic system closes at the Macaulay bound");
        assert_eq!(st.dim, (d as usize).pow(4), "{st:?}");
        assert!(sols.contains(&root), "planted {root:?} not in {sols:?}");
        for s in &sols {
            let mut dummy = 0u64;
            assert!(comps.iter().all(|c| eval4(c, s, p, &mut dummy) == 0));
        }
        assert!(st.total_muls > 0 && st.echelon_muls > 0 && st.charpoly_muls > 0);
    }

    #[test]
    fn the_solver_finds_planted_roots_of_quadrics() {
        for seed in 0..4 {
            planted(2, seed);
        }
    }

    #[test]
    fn the_solver_finds_planted_roots_of_cubics() {
        planted(3, 7);
    }

    /// §11.18's pipeline end to end at a toy size: a prime-order curve,
    /// relations from the meet-in-the-middle oracle, filtering, Wiedemann,
    /// the planted logarithm recovered — and rho recovering it too.
    #[test]
    fn k4_linear_algebra_and_rho_recover_the_planted_logarithm() {
        let rep = run_k4_la(29, 3, 4);
        assert!(rep.solved && rep.correct, "{rep:?}");
        assert!(super::super::residual_walk::is_prime_u64(rep.n));
        assert!(rep.rho.iter().all(|r| r.correct), "{:?}", rep.rho);
        assert!(rep.la_ops > 0 && rep.r > 0.0);
        assert!(
            (rep.row_weight - 5.0).abs() < 1e-9,
            "four points and d: {}",
            rep.row_weight
        );
        assert_eq!(rep.fp_muls_per_add, 97.0);
    }
}
