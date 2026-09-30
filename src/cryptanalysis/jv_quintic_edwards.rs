//! **The `k = 5` stage with the 2-torsion symmetry: Joux–Vitse four-point
//! decompositions on an Edwards curve over `F_{p⁵}`, in the `y`-coordinate.**
//!
//! Companion to `research/notes/index-calculus/RESEARCH_RHO_PARITY_PROGRAMME.md`
//! §6 and to [`super::jv_quintic`], which measured the plain symmetrisation.
//! On `x² + y² = 1 + d x² y²` the point `T = (0, −1)` has order two and
//! `P + T = (−x, −y)`, `−P = (−x, y)`: the `y`-coordinate is invariant under
//! negation and changes sign under `T`.  The relation `P₁ + ⋯ + P₄ = R` is
//! preserved by translating an even number of its points by `T`, so the
//! summation polynomial in `(y₁, …, y₄, y_R)` has every monomial all-even
//! or all-odd in its exponents (Faugère–Gaudry–Huot–Renault):
//! `S₅ = A(y²; y_R²) + (y₁y₂y₃y₄ y_R) · B(y²; y_R²)`.  Symmetrised over `S₄`
//! it is a polynomial in the elementary symmetric functions `e_i` of the
//! `y_i²` (total degree `≤ 4` in `A`, `≤ 3` in `B`) and one product
//! variable `π = y₁y₂y₃y₄` with `π² = e₄` — degrees `4` and `3` where the
//! Weierstrass symmetrisation has `8`.  The factor base `{P : y(P) ∈ F_p}`
//! is stable under `⟨−1, T⟩`, its columns are the `y²` classes (`≈ p/4` of
//! them), and a relation `R = Σ s_i P_i + εT` is used through `[4]`, which
//! kills `T`.
//!
//! `S₃` in the `y`-coordinate, derived from the addition law with the
//! spurious factor `1 − d y₁²y₂²` removed:
//! `S₃ = (y₁² + y₂² + y₃²) − d(y₁²y₂² + y₁²y₃² + y₂²y₃²) + d y₁²y₂²y₃²
//!       − 2(1 − d) y₁y₂y₃ − 1`.

use std::cell::Cell;
use std::collections::{BTreeSet, HashMap};
use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::Serialize;

use super::f4_fp::{self, F4Options, Ordering as F4Ordering, Verdict};
use super::gaudry_cubic::{mm, PolyRing, UPoly};
use super::gaudry_quartic::exponents4;
use super::jv_quintic::{Fp5, E5};
use super::residual_walk::{is_prime_u64, pow_mod};

// ── The Edwards curve x² + y² = 1 + d x² y² over F_{p⁵} ─────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PtE {
    pub x: E5,
    pub y: E5,
}

impl PtE {
    pub const ID: PtE = PtE {
        x: E5::ZERO,
        y: E5::ONE,
    };
    pub fn is_identity(&self) -> bool {
        self.x.is_zero() && self.y == E5::ONE
    }
}

pub struct Ed5 {
    pub f: Fp5,
    pub d: E5,
    ops: Cell<u64>,
}

impl Ed5 {
    /// `d` a non-square of `F_{p⁵}` outside `F_p`: the addition law is then
    /// complete (`a = 1` is a square), and the curve is not defined over a
    /// subfield.
    pub fn random(p: u64, rng: &mut StdRng) -> Ed5 {
        let f = Fp5::new(p);
        loop {
            let d = f.random(rng);
            if d.in_fp() || d.is_zero() || f.is_square(&d) {
                continue;
            }
            return Ed5 {
                f,
                d,
                ops: Cell::new(0),
            };
        }
    }
    pub fn with_params(p: u64, d: E5) -> Ed5 {
        Ed5 {
            f: Fp5::new(p),
            d,
            ops: Cell::new(0),
        }
    }
    pub fn t(&self) -> PtE {
        PtE {
            x: E5::ZERO,
            y: self.f.neg(&E5::ONE),
        }
    }
    pub fn ops(&self) -> u64 {
        self.ops.get()
    }
    pub fn reset_ops(&self) {
        self.ops.set(0);
    }
    pub fn on_curve(&self, pt: &PtE) -> bool {
        let f = &self.f;
        let x2 = f.sq(&pt.x);
        let y2 = f.sq(&pt.y);
        f.add(&x2, &y2) == f.add(&E5::ONE, &f.mul(&self.d, &f.mul(&x2, &y2)))
    }
    pub fn neg(&self, pt: &PtE) -> PtE {
        PtE {
            x: self.f.neg(&pt.x),
            y: pt.y,
        }
    }
    /// `x² = (1 − y²) / (1 − d y²)`; the point with that `y` and the
    /// square root the field returns, if any.
    pub fn lift_y(&self, y: &E5, rng: &mut StdRng) -> Option<PtE> {
        let f = &self.f;
        let y2 = f.sq(y);
        let den = f.sub(&E5::ONE, &f.mul(&self.d, &y2));
        if den.is_zero() {
            return None;
        }
        let x2 = f.mul(&f.sub(&E5::ONE, &y2), &f.inv(&den));
        let x = f.sqrt(&x2, rng)?;
        Some(PtE { x, y: *y })
    }
    /// The unified addition law, one operation: two products for the
    /// denominators, one inversion shared through their product.
    pub fn add(&self, p1: &PtE, p2: &PtE) -> PtE {
        let f = &self.f;
        self.ops.set(self.ops.get() + 1);
        let x1x2 = f.mul(&p1.x, &p2.x);
        let y1y2 = f.mul(&p1.y, &p2.y);
        let x1y2 = f.mul(&p1.x, &p2.y);
        let y1x2 = f.mul(&p1.y, &p2.x);
        let k = f.mul(&self.d, &f.mul(&x1x2, &y1y2));
        let den_x = f.add(&E5::ONE, &k);
        let den_y = f.sub(&E5::ONE, &k);
        let inv = f.inv(&f.mul(&den_x, &den_y));
        let x3 = f.mul(&f.add(&x1y2, &y1x2), &f.mul(&inv, &den_y));
        let y3 = f.mul(&f.sub(&y1y2, &x1x2), &f.mul(&inv, &den_x));
        PtE { x: x3, y: y3 }
    }
    pub fn sub(&self, p1: &PtE, p2: &PtE) -> PtE {
        self.add(p1, &self.neg(p2))
    }
    pub fn mul(&self, pt: &PtE, mut k: u64) -> PtE {
        let mut acc = PtE::ID;
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

fn random_point(curve: &Ed5, rng: &mut StdRng) -> PtE {
    loop {
        let y = curve.f.random(rng);
        if let Some(pt) = curve.lift_y(&y, rng) {
            if !pt.is_identity() && pt != curve.t() {
                return pt;
            }
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

/// `#E(F_{p⁵})` by baby-step giant-step on a random point, accepted when it
/// is four times a prime; the prime-order subgroup is then `[4]E`.
fn group_order_4_prime(curve: &Ed5, rng: &mut StdRng) -> Option<u64> {
    let p = curve.f.p;
    let q = p.pow(5);
    let two_sqrt = 2 * isqrt(q) + 2;
    let lo = q + 1 - two_sqrt;
    let width = 2 * two_sqrt;
    let pt = random_point(curve, rng);
    let steps = isqrt(width) + 1;
    let mut table: HashMap<PtE, u64> = HashMap::new();
    let mut jp = PtE::ID;
    for j in 0..steps {
        if !jp.is_identity() {
            table.entry(jp).or_insert(j);
        }
        jp = curve.add(&jp, &pt);
    }
    let giant = curve.mul(&pt, steps);
    let mut t = curve.mul(&pt, lo);
    let mut i = 0u64;
    let mut found = None;
    while i * steps <= width + steps {
        let m = if t.is_identity() {
            Some(lo + i * steps)
        } else {
            table.get(&curve.neg(&t)).map(|&j| lo + i * steps + j)
        };
        if let Some(m) = m {
            found = Some(m);
            break;
        }
        t = curve.add(&t, &giant);
        i += 1;
    }
    let m = found?;
    if m % 4 != 0 || !is_prime_u64(m / 4) {
        return None;
    }
    // The order must annihilate a second random point too.
    let other = random_point(curve, rng);
    curve.mul(&other, m).is_identity().then_some(m)
}

pub struct InstanceE {
    pub curve: Ed5,
    /// The prime order of `[4]E`.
    pub n: u64,
    pub g: PtE,
    pub d: u64,
    pub q: PtE,
}

pub fn generate_instance_edwards(p: u64, seed: u64) -> InstanceE {
    let mut rng = StdRng::seed_from_u64(seed ^ 0xED05);
    loop {
        let curve = Ed5::random(p, &mut rng);
        let Some(m) = group_order_4_prime(&curve, &mut rng) else {
            continue;
        };
        let n = m / 4;
        let g = loop {
            let g = curve.mul(&random_point(&curve, &mut rng), 4);
            if !g.is_identity() && curve.mul(&g, n).is_identity() {
                break g;
            }
        };
        let d = rng.gen_range(1..n);
        let q = curve.mul(&g, d);
        curve.reset_ops();
        curve.f.reset_muls();
        return InstanceE { curve, n, g, d, q };
    }
}

/// One point per `y²` class: `y ∈ [2, (p − 1)/2]` with `(1 − y²)/(1 − d y²)`
/// a square, `x` the root the field returns.  `y = 0` gives the 4-torsion
/// points and `y = ±1` the identity and `T`; `−y` gives `±P + T`.
pub fn factor_base(curve: &Ed5, rng: &mut StdRng) -> Vec<PtE> {
    let p = curve.f.p;
    (2..=(p - 1) / 2)
        .filter_map(|y| curve.lift_y(&curve.f.from_fp(y), rng))
        .collect()
}

/// `y(P)² ∈ F_p` of a base point, the column key.
pub fn y2_key(curve: &Ed5, pt: &PtE) -> u64 {
    mm(pt.y.0[0], pt.y.0[0], curve.f.p)
}

// ── Summation polynomials in y ───────────────────────────────────────────

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

impl Ed5 {
    /// `S₃(y₁, y₂, Z)` as a quadratic in `Z`: `[c₀, c₁, c₂]` with
    /// `c₂ = 1 − d(y₁² + y₂²) + d y₁²y₂²`, `c₁ = −2(1 − d) y₁y₂`,
    /// `c₀ = y₁² + y₂² − d y₁²y₂² − 1`.
    fn s3_in_last(&self, y1: &E5, y2: &E5) -> [E5; 3] {
        let f = &self.f;
        let d = &self.d;
        let (a, b) = (f.sq(y1), f.sq(y2));
        let ab = f.mul(&a, &b);
        let s = f.add(&a, &b);
        let c2 = f.add(&f.sub(&E5::ONE, &f.mul(d, &s)), &f.mul(d, &ab));
        let one_minus_d = f.sub(&E5::ONE, d);
        let u = f.mul(y1, y2);
        let c1 = f.neg(&f.scale(&f.mul(&one_minus_d, &u), 2));
        let c0 = f.sub(&f.sub(&s, &f.mul(d, &ab)), &E5::ONE);
        [c0, c1, c2]
    }
    /// `S₃(y₁, y₂, y₃)`.
    pub fn s3(&self, y: &[E5; 3]) -> E5 {
        let f = &self.f;
        let [c0, c1, c2] = self.s3_in_last(&y[0], &y[1]);
        f.add(&f.add(&f.mul(&c2, &f.sq(&y[2])), &f.mul(&c1, &y[2])), &c0)
    }
    /// `S₃(y₅, Y, Z)` with `Y` symbolic: `Z`-coefficients as polynomials in `Y`.
    fn s3_in_last_symbolic(&self, y5: &E5) -> [Vec<E5>; 3] {
        let f = &self.f;
        let d = &self.d;
        let a = f.sq(y5);
        // c₂(Y) = (1 − d a) − d(1 − a) Y²
        let c2 = vec![
            f.sub(&E5::ONE, &f.mul(d, &a)),
            E5::ZERO,
            f.neg(&f.mul(d, &f.sub(&E5::ONE, &a))),
        ];
        // c₁(Y) = −2(1 − d) y₅ Y
        let c1 = vec![
            E5::ZERO,
            f.neg(&f.scale(&f.mul(&f.sub(&E5::ONE, d), y5), 2)),
        ];
        // c₀(Y) = (a − 1) + (1 − d a) Y²
        let c0 = vec![
            f.sub(&a, &E5::ONE),
            E5::ZERO,
            f.sub(&E5::ONE, &f.mul(d, &a)),
        ];
        [c0, c1, c2]
    }
    /// `S₄(y₃, y₄, y₅, Y) = Res_Z(S₃(y₃, y₄, Z), S₃(y₅, Y, Z))`, degree ≤ 4 in `Y`.
    fn s4_in_last(&self, y3: &E5, y4: &E5, y5: &E5) -> Vec<E5> {
        let f = &self.f;
        let [a0, a1, a2] = self.s3_in_last(y3, y4);
        let [b0, b1, b2] = self.s3_in_last_symbolic(y5);
        let p1 = poly_sub(f, &poly_scale(f, &b0, &a2), &poly_scale(f, &b2, &a0));
        let p2 = poly_sub(f, &poly_scale(f, &b1, &a2), &poly_scale(f, &b2, &a1));
        let p3 = poly_sub(f, &poly_scale(f, &b0, &a1), &poly_scale(f, &b1, &a0));
        poly_sub(f, &poly_mul(f, &p1, &p1), &poly_mul(f, &p2, &p3))
    }
    /// `S₅(y₁, …, y₅) = Res_Y(S₃(y₁, y₂, Y), S₄(y₃, y₄, y₅, Y))`.
    pub fn s5(&self, y: &[E5; 5]) -> E5 {
        let f = &self.f;
        let a = self.s3_in_last(&y[0], &y[1]);
        let mut b = self.s4_in_last(&y[2], &y[3], &y[4]);
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

// ── The symmetrised S₅ in the squares, by interpolation ──────────────────

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
    let mut pows = [[E5::ONE; 5]; 4];
    for i in 0..4 {
        for k in 1..5 {
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

/// Solve `[M | C]` over `F_{p⁵}` for the coefficient block; `None` if singular.
fn solve_block(f: &Fp5, mut sys: Vec<Vec<E5>>, n: usize) -> Option<Vec<Vec<E5>>> {
    let p = f.p;
    for c in 0..n {
        let r = (c..n).find(|&r| !sys[r][c].is_zero())?;
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
        f.count_public(29 * (n as u64) * (sys[0].len() - c) as u64);
    }
    Some(sys.into_iter().map(|row| row[n..].to_vec()).collect())
}

/// `S₅(y₁, …, y₄, Y) = A(e; Y) + π · B(e; Y)` with `e = e(y²)`,
/// `π = y₁y₂y₃y₄`, `A` of total degree `≤ 4` in `e` and even degree `≤ 8`
/// in `Y` (five coefficients), `B` of total degree `≤ 3` and odd degree
/// `≤ 7` in `Y` (four coefficients).
pub struct SymmetrisedS5Ed {
    pub monos_a: Vec<[u8; 4]>,
    pub coef_a: Vec<[E5; 5]>,
    pub monos_b: Vec<[u8; 4]>,
    pub coef_b: Vec<[E5; 4]>,
    pub precompute_muls: u64,
}

impl SymmetrisedS5Ed {
    /// The nine `Y`-coefficients of `S₅(y, Y)` from values at `Y = 0..8`.
    fn coefficients_in_last(curve: &Ed5, y: &[E5; 4], vinv: &[[u64; 9]; 9]) -> [E5; 9] {
        let f = &curve.f;
        let vals: Vec<E5> = (0..9u64)
            .map(|k| curve.s5(&[y[0], y[1], y[2], y[3], f.from_fp(k)]))
            .collect();
        core::array::from_fn(|j| {
            vals.iter()
                .enumerate()
                .fold(E5::ZERO, |acc, (k, v)| f.add(&acc, &f.scale(v, vinv[j][k])))
        })
    }

    pub fn precompute(curve: &Ed5, rng: &mut StdRng) -> SymmetrisedS5Ed {
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
        let monos_a = exponents4(4);
        let monos_b = exponents4(3);
        let (na, nb) = (monos_a.len(), monos_b.len());
        loop {
            let mut sys_a: Vec<Vec<E5>> = Vec::with_capacity(na);
            let mut sys_b: Vec<Vec<E5>> = Vec::with_capacity(nb);
            let mut samples = 0;
            while sys_a.len() < na || sys_b.len() < nb {
                let y: [E5; 4] = core::array::from_fn(|_| f.random(rng));
                let pi = f.mul(&f.mul(&y[0], &y[1]), &f.mul(&y[2], &y[3]));
                if pi.is_zero() {
                    continue;
                }
                let sq: [E5; 4] = core::array::from_fn(|i| f.sq(&y[i]));
                let e = elementary4(f, &sq);
                let c = Self::coefficients_in_last(curve, &y, &vinv);
                if sys_a.len() < na {
                    let mut row = mono_values4(f, &e, &monos_a);
                    row.extend([c[0], c[2], c[4], c[6], c[8]]);
                    sys_a.push(row);
                }
                if sys_b.len() < nb {
                    let inv = f.inv(&pi);
                    let mut row = mono_values4(f, &e, &monos_b);
                    row.extend([c[1], c[3], c[5], c[7]].map(|v| f.mul(&v, &inv)));
                    sys_b.push(row);
                }
                samples += 1;
            }
            let _ = samples;
            let (Some(ca), Some(cb)) = (solve_block(f, sys_a, na), solve_block(f, sys_b, nb))
            else {
                continue;
            };
            let pre = SymmetrisedS5Ed {
                monos_a: monos_a.clone(),
                coef_a: ca.iter().map(|r| core::array::from_fn(|j| r[j])).collect(),
                monos_b: monos_b.clone(),
                coef_b: cb.iter().map(|r| core::array::from_fn(|j| r[j])).collect(),
                precompute_muls: 0,
            };
            let good = (0..6).all(|_| {
                let y: [E5; 4] = core::array::from_fn(|_| f.random(rng));
                let yr = f.random(rng);
                let want = curve.s5(&[y[0], y[1], y[2], y[3], yr]);
                pre.eval(f, &y, &yr) == want
            });
            assert!(
                good,
                "the symmetrised S₅ (Edwards) does not reproduce the resultant"
            );
            return SymmetrisedS5Ed {
                precompute_muls: f.muls() - before,
                ..pre
            };
        }
    }

    /// `A(e(y²); Y) + π B(e(y²); Y)` at a point.
    pub fn eval(&self, f: &Fp5, y: &[E5; 4], y_r: &E5) -> E5 {
        let sq: [E5; 4] = core::array::from_fn(|i| f.sq(&y[i]));
        let e = elementary4(f, &sq);
        let pi = f.mul(&f.mul(&y[0], &y[1]), &f.mul(&y[2], &y[3]));
        let mut pows = [E5::ONE; 9];
        for k in 1..9 {
            pows[k] = f.mul(&pows[k - 1], y_r);
        }
        let ma = mono_values4(f, &e, &self.monos_a);
        let a = ma.iter().zip(&self.coef_a).fold(E5::ZERO, |acc, (m, c)| {
            let v = (0..5).fold(E5::ZERO, |s, j| f.add(&s, &f.mul(&c[j], &pows[2 * j])));
            f.add(&acc, &f.mul(m, &v))
        });
        let mb = mono_values4(f, &e, &self.monos_b);
        let b = mb.iter().zip(&self.coef_b).fold(E5::ZERO, |acc, (m, c)| {
            let v = (0..4).fold(E5::ZERO, |s, j| f.add(&s, &f.mul(&c[j], &pows[2 * j + 1])));
            f.add(&acc, &f.mul(m, &v))
        });
        f.add(&a, &f.mul(&pi, &b))
    }

    /// The five `F_p`-components of `A(e; y_R) + π B(e; y_R)` as polynomials
    /// in `(e₁, e₂, e₃, e₄, π)`: exponents `[a, b, c, d, 0]` and `[a, b, c, d, 1]`.
    pub fn weil_restrict(&self, f: &Fp5, y_r: &E5) -> [HashMap<[u8; 5], u64>; 5] {
        let mut pows = [E5::ONE; 9];
        for k in 1..9 {
            pows[k] = f.mul(&pows[k - 1], y_r);
        }
        let mut out: [HashMap<[u8; 5], u64>; 5] = Default::default();
        for (m, c) in self.monos_a.iter().zip(&self.coef_a) {
            let v = (0..5).fold(E5::ZERO, |s, j| f.add(&s, &f.mul(&c[j], &pows[2 * j])));
            for (i, comp) in out.iter_mut().enumerate() {
                if v.0[i] != 0 {
                    comp.insert([m[0], m[1], m[2], m[3], 0], v.0[i]);
                }
            }
        }
        for (m, c) in self.monos_b.iter().zip(&self.coef_b) {
            let v = (0..4).fold(E5::ZERO, |s, j| f.add(&s, &f.mul(&c[j], &pows[2 * j + 1])));
            for (i, comp) in out.iter_mut().enumerate() {
                if v.0[i] != 0 {
                    comp.insert([m[0], m[1], m[2], m[3], 1], v.0[i]);
                }
            }
        }
        out
    }
}

// ── The oracle ────────────────────────────────────────────────────────────

/// One decomposition `R = Σ s_t P_{i_t} + ε T` over four distinct base points.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub struct QuadE {
    pub terms: [(usize, i64); 4],
    pub eps: u8,
}

pub struct PairTableE {
    map: HashMap<u64, Vec<(u32, u32, i8, PtE)>>,
}

/// `y(P_i + s P_j)² ↦ (i, j, s, the point)`, built once per curve.
pub fn pair_table(curve: &Ed5, base: &[PtE]) -> PairTableE {
    let mut map: HashMap<u64, Vec<(u32, u32, i8, PtE)>> = HashMap::new();
    for i in 0..base.len() {
        for j in (i + 1)..base.len() {
            for (s, v) in [
                (1i8, curve.add(&base[i], &base[j])),
                (-1, curve.sub(&base[i], &base[j])),
            ] {
                map.entry(y2_key(curve, &v))
                    .or_default()
                    .push((i as u32, j as u32, s, v));
            }
        }
    }
    PairTableE { map }
}

/// `v` against a table point `W`: `W`, `−W`, `W + T` or `−W + T`.
fn classify(curve: &Ed5, v: &PtE, w: &PtE) -> Option<(i64, u8)> {
    let f = &curve.f;
    if v.y == w.y {
        if v.x == w.x {
            return Some((1, 0));
        }
        if v.x == f.neg(&w.x) {
            return Some((-1, 0));
        }
    } else if v.y == f.neg(&w.y) {
        if v.x == f.neg(&w.x) {
            return Some((1, 1));
        }
        if v.x == w.x {
            return Some((-1, 1));
        }
    }
    None
}

/// Every `R = s_k P_k + s_l P_l ± (P_i + s P_j) + εT` by meet in the middle.
pub fn mitm4_signed(curve: &Ed5, base: &[PtE], table: &PairTableE, r: &PtE) -> Vec<QuadE> {
    let n = base.len();
    let mut out: BTreeSet<QuadE> = BTreeSet::new();
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
                    let Some(cands) = table.map.get(&y2_key(curve, &v)) else {
                        continue;
                    };
                    for &(i, j, s, w) in cands {
                        let (i, j) = (i as usize, j as usize);
                        if i == k || i == l || j == k || j == l {
                            continue;
                        }
                        if let Some((t, eps)) = classify(curve, &v, &w) {
                            let mut terms = [(k, sk), (l, sl), (i, t), (j, t * s as i64)];
                            terms.sort_unstable();
                            out.insert(QuadE { terms, eps });
                        }
                    }
                }
            }
        }
    }
    out.into_iter().collect()
}

pub fn verify_quad(curve: &Ed5, base: &[PtE], r: &PtE, q: &QuadE) -> bool {
    let mut acc = *r;
    for &(i, s) in &q.terms {
        acc = if s == 1 {
            curve.sub(&acc, &base[i])
        } else {
            curve.add(&acc, &base[i])
        };
    }
    if q.eps == 1 {
        acc = curve.sub(&acc, &curve.t());
    }
    acc.is_identity()
}

// ── The four-point test ───────────────────────────────────────────────────

#[derive(Clone, Debug, Default, Serialize)]
pub struct Jv4CostE {
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

fn polys_of(comps: &[HashMap<[u8; 5], u64>; 5], p: u64) -> Vec<f4_fp::Poly> {
    let mut polys: Vec<f4_fp::Poly> = comps
        .iter()
        .map(|c| {
            let terms: Vec<(Vec<u32>, u64)> = c
                .iter()
                .map(|(e, v)| (e.iter().map(|&x| x as u32).collect(), *v))
                .collect();
            f4_fp::normalise(&terms, p, F4Ordering::Grevlex)
        })
        .filter(|f| !f.is_empty())
        .collect();
    // π² − e₄ = 0
    polys.push(f4_fp::normalise(
        &[(vec![0, 0, 0, 0, 2], 1), (vec![0, 0, 0, 1, 0], p - 1)],
        p,
        F4Ordering::Grevlex,
    ));
    polys
}

pub fn jv4_options(max_degree: u32, budget_secs: f64) -> F4Options {
    F4Options::new(F4Ordering::Grevlex, max_degree)
        .with_budget(std::time::Duration::from_secs_f64(budget_secs))
}

/// Signs and the `T`-component of a candidate index quadruple by group
/// arithmetic: `R − Σ s_i P_i ∈ {O, T}` for one of the sixteen sign vectors.
fn resolve(curve: &Ed5, base: &[PtE], r: &PtE, idx: [usize; 4]) -> Option<QuadE> {
    let t = curve.t();
    for mask in 0..16u32 {
        let mut acc = *r;
        let mut terms = [(0usize, 0i64); 4];
        for (k, &i) in idx.iter().enumerate() {
            let s: i64 = if mask & (1 << k) == 0 { 1 } else { -1 };
            acc = if s == 1 {
                curve.sub(&acc, &base[i])
            } else {
                curve.add(&acc, &base[i])
            };
            terms[k] = (i, s);
        }
        let eps = if acc.is_identity() {
            0
        } else if acc == t {
            1
        } else {
            continue;
        };
        terms.sort_unstable();
        return Some(QuadE { terms, eps });
    }
    None
}

/// Every four-point decomposition of `r` over the base, by the symmetrised
/// Edwards `S₅` and F4, with the cost split.
pub fn jv4_decompose(
    curve: &Ed5,
    base: &[PtE],
    by_y2: &HashMap<u64, usize>,
    pre: &SymmetrisedS5Ed,
    r: &PtE,
    opts: &F4Options,
    rng: &mut StdRng,
) -> (Vec<QuadE>, Jv4CostE) {
    let f = &curve.f;
    let p = f.p;
    let mut cost = Jv4CostE::default();
    let mut out: BTreeSet<QuadE> = BTreeSet::new();
    let m0 = f.muls();
    let comps = pre.weil_restrict(f, &r.y);
    cost.weil_muls = f.muls() - m0;
    let polys = polys_of(&comps, p);
    let rep = f4_fp::solve(&polys, 5, p, opts);
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
        // T⁴ − e₁T³ + e₂T² − e₃T + e₄: the roots are the y_i².
        let quartic = UPoly(vec![e[3], (p - e[2]) % p, e[1], (p - e[0]) % p, 1]);
        let roots = ring.roots(&quartic, rng);
        if roots.len() != 4 {
            continue;
        }
        let Some(idx) = roots
            .iter()
            .map(|x| by_y2.get(x).copied())
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
        if let Some(q) = resolve(curve, base, r, [idx[0], idx[1], idx[2], idx[3]]) {
            out.insert(q);
        }
    }
    cost.roots_muls = ring.muls.get();
    cost.group_ops = curve.ops() - g0;
    cost.decompositions = out.len();
    (out.into_iter().collect(), cost)
}

// ── Experiment: C″ with the symmetry ─────────────────────────────────────

#[derive(Clone, Debug, Default, Serialize)]
pub struct CostStatsE {
    pub mean: f64,
    pub min: u64,
    pub max: u64,
}

fn stats(v: &[u64]) -> CostStatsE {
    if v.is_empty() {
        return CostStatsE::default();
    }
    CostStatsE {
        mean: v.iter().sum::<u64>() as f64 / v.len() as f64,
        min: *v.iter().min().unwrap(),
        max: *v.iter().max().unwrap(),
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct CsecondEdReport {
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
    /// `1/(192p)`: eight subgroup points per class-quadruple of the `p/4`
    /// classes (`rate_census`).  The frozen run `27_jv_quintic_edwards_csecond`
    /// carries this field as `1/(24p)`, the Weierstrass rate, an error the
    /// census corrected; the ledger recomputes the rate from `p`.
    pub expected_rate: f64,
    pub planted_found: usize,
    pub unverified: usize,
    pub mismatches: usize,
    pub undetermined: usize,
    pub timed_out: usize,
    pub c_second: CostStatsE,
    pub weil_muls: CostStatsE,
    pub f4_muls: CostStatsE,
    pub roots_muls: CostStatsE,
    pub sign_group_ops: CostStatsE,
    pub f4_ms: CostStatsE,
    pub solving_degree: CostStatsE,
    pub degree_reached: CostStatsE,
    pub max_rows: CostStatsE,
    pub max_cols: CostStatsE,
    pub c_second_decomposable: CostStatsE,
    pub f4_ms_decomposable: CostStatsE,
    pub degree_reached_decomposable: CostStatsE,
    pub max_cols_decomposable: CostStatsE,
    pub per_test: Vec<Jv4CostE>,
    pub mitm4_group_ops: f64,
    pub wall_ms: f64,
}

pub fn run_jv5_edwards_csecond(
    p: u64,
    seed: u64,
    random: usize,
    constructed: usize,
    max_degree: u32,
    budget_secs: f64,
) -> CsecondEdReport {
    let start = Instant::now();
    let inst = generate_instance_edwards(p, seed);
    let curve = &inst.curve;
    let f = &curve.f;
    let n = inst.n;
    let mut rng = StdRng::seed_from_u64(seed ^ 0xED02);
    let base = factor_base(curve, &mut rng);
    let by_y2: HashMap<u64, usize> = base
        .iter()
        .enumerate()
        .map(|(i, q)| (y2_key(curve, q), i))
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
    let pre = SymmetrisedS5Ed::precompute(curve, &mut rng);
    let precompute_ms = t0.elapsed().as_secs_f64() * 1e3;
    let table = pair_table(curve, &base);
    let mut rep = CsecondEdReport {
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
        expected_rate: 1.0 / (192.0 * p as f64),
        ..Default::default()
    };
    let mut costs: Vec<Jv4CostE> = Vec::new();
    let mut costs_dec: Vec<Jv4CostE> = Vec::new();
    let mut mitm_ops = 0u64;
    let mut check_rng = StdRng::seed_from_u64(seed ^ 0x5EED);
    let mut check = |r: &PtE, want: Option<QuadE>, rep: &mut CsecondEdReport| -> Jv4CostE {
        let opts = jv4_options(max_degree, budget_secs);
        let (found, cost) = jv4_decompose(curve, &base, &by_y2, &pre, r, &opts, &mut check_rng);
        rep.unverified += found
            .iter()
            .filter(|q| !verify_quad(curve, &base, r, q))
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
        if let Some(w) = want {
            if found.contains(&w) {
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
        let eps: u8 = if rng.gen_bool(0.5) { 1 } else { 0 };
        let mut r = PtE::ID;
        for t in 0..4 {
            let q = if signs[t] == 1 {
                base[idx[t]]
            } else {
                curve.neg(&base[idx[t]])
            };
            r = curve.add(&r, &q);
        }
        if eps == 1 {
            r = curve.add(&r, &curve.t());
        }
        let mut terms = [
            (idx[0], signs[0]),
            (idx[1], signs[1]),
            (idx[2], signs[2]),
            (idx[3], signs[3]),
        ];
        terms.sort_unstable();
        let c = check(&r, Some(QuadE { terms, eps }), &mut rep);
        costs_dec.push(c);
    }
    let total = |c: &Jv4CostE| {
        c.weil_muls + c.f4_muls + c.roots_muls + (c.group_ops as f64 * fp_per_add) as u64
    };
    let done = |v: &[Jv4CostE]| -> Vec<Jv4CostE> {
        v.iter().filter(|c| !c.undetermined).cloned().collect()
    };
    let col = |v: &[Jv4CostE], g: &dyn Fn(&Jv4CostE) -> u64| -> CostStatsE {
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

// ── The rational 4-torsion point: saturation, not symmetry ───────────────
//
// `Q₄ = (1, 0)` has order four (`2Q₄ = T`) and translates by
// `P + Q₄ = (y, −x)`, `P − Q₄ = (−y, x)`: it swaps the coordinates.  No
// function of `y` alone is `Q₄`-invariant, so — unlike `T`, which the
// `y²` classes absorb — `Q₄` cannot halve the degree of the symmetrised
// `S₅`.  The second degree-halving symmetry of the Edwards form is the
// *other* 2-torsion point, `y ↦ −1/(√d · y)` (the identity
// `d y₁²y₂² · S₃(−1/(√d y₁), −1/(√d y₂), y₃) = S₃(y₁, y₂, y₃)`, tested
// below); it keeps `y ∈ F_p` only when `√d ∈ F_p`, i.e. on a curve
// defined over `F_p`, which the programme excludes.
//
// What `Q₄` does buy is saturation of the residual: `R` is decomposable
// through `Q₄` when `R − t_q Q₄` is a four-point sum for some `t_q ∈ Z/4`,
// and since `y(R − Q₄) = x_R` and `y(R + Q₄) = −x_R` the classes `t_q`
// odd share one `y²` value, `x_R²`.  A relation `R = Σ s_i P_i + t_q Q₄`
// is used through `[4]`, as before.  The rate doubles (`1/(96p)` against
// `1/(192p)`, §`rate_census`) and the test runs F4 twice (`y_R`, `x_R`):
// a wash, measured by `run_jv5_edwards4_csecond`.

impl Ed5 {
    /// `Q₄ = (1, 0)`: `2Q₄ = T`, `P + Q₄ = (y, −x)`, `P − Q₄ = (−y, x)`.
    pub fn q4(&self) -> PtE {
        PtE {
            x: E5::ONE,
            y: E5::ZERO,
        }
    }
    /// `P − Q₄ = (−y, x)` without a group operation.
    pub fn sub_q4(&self, pt: &PtE) -> PtE {
        PtE {
            x: self.f.neg(&pt.y),
            y: pt.x,
        }
    }
    /// The homomorphism `E → Z/4` with kernel the prime-order subgroup:
    /// `φ(P) = c` with `nP = c · (nQ₄)`.  `φ(Q₄) = 1`, `φ(T) = 2`.
    pub fn torsion_class(&self, n: u64, pt: &PtE) -> u8 {
        let np = self.mul(pt, n);
        let nq = self.mul(&self.q4(), n);
        let mut acc = PtE::ID;
        for c in 0..4u8 {
            if acc == np {
                return c;
            }
            acc = self.add(&acc, &nq);
        }
        panic!("nP outside <nQ₄>: n is not the odd part of the order");
    }
}

/// One decomposition `R = Σ s_t P_{i_t} + t_q Q₄`, `t_q ∈ Z/4`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub struct QuadE4 {
    pub terms: [(usize, i64); 4],
    pub tq: u8,
}

fn saturate(quads: &[QuadE], offset: u8) -> impl Iterator<Item = QuadE4> + '_ {
    quads.iter().map(move |q| QuadE4 {
        terms: q.terms,
        tq: (offset + 2 * q.eps) % 4,
    })
}

/// Every `R = Σ s_i P_i + t_q Q₄`: the oracle on `R` (`t_q = 2ε`) and on
/// `R − Q₄` (`t_q = 1 + 2ε`).
pub fn mitm4_saturated(curve: &Ed5, base: &[PtE], table: &PairTableE, r: &PtE) -> Vec<QuadE4> {
    let mut out: BTreeSet<QuadE4> = BTreeSet::new();
    out.extend(saturate(&mitm4_signed(curve, base, table, r), 0));
    out.extend(saturate(
        &mitm4_signed(curve, base, table, &curve.sub_q4(r)),
        1,
    ));
    out.into_iter().collect()
}

pub fn verify_quad4(curve: &Ed5, base: &[PtE], r: &PtE, q: &QuadE4) -> bool {
    let mut acc = *r;
    for &(i, s) in &q.terms {
        acc = if s == 1 {
            curve.sub(&acc, &base[i])
        } else {
            curve.add(&acc, &base[i])
        };
    }
    for _ in 0..q.tq {
        acc = curve.sub_q4(&acc);
    }
    acc.is_identity()
}

/// The four-point test on both translates: `jv4_decompose` at `y_R` and at
/// `y(R − Q₄) = x_R`, with the cost of each run.
pub fn jv4_decompose_saturated(
    curve: &Ed5,
    base: &[PtE],
    by_y2: &HashMap<u64, usize>,
    pre: &SymmetrisedS5Ed,
    r: &PtE,
    opts: &F4Options,
    rng: &mut StdRng,
) -> (Vec<QuadE4>, [Jv4CostE; 2]) {
    let mut out: BTreeSet<QuadE4> = BTreeSet::new();
    let (q0, c0) = jv4_decompose(curve, base, by_y2, pre, r, opts, rng);
    out.extend(saturate(&q0, 0));
    let (q1, c1) = jv4_decompose(curve, base, by_y2, pre, &curve.sub_q4(r), opts, rng);
    out.extend(saturate(&q1, 1));
    (out.into_iter().collect(), [c0, c1])
}

// ── The rate census ───────────────────────────────────────────────────────

/// Every signed four-point sum over the base, counted exactly: how many
/// distinct points of the prime-order subgroup are `Σ s_i P_i + εT`
/// (2-torsion) and `Σ s_i P_i + t_q Q₄` (saturated).  Per class-quadruple
/// the sixteen sign vectors give sixteen points.  With `Q₄` free each
/// lands in the subgroup for exactly one of the four translates: `16` per
/// quadruple.  With `T` free a point lands in the subgroup for one of the
/// two translates when its `Z/4` class is even, and the parity of
/// `Σ s_i c_i` does not depend on the signs: a quadruple whose classes
/// sum to an even number gives `16`, an odd one gives none, `8` on
/// average.  With `|F| = p/4` this is `1/(192p)` and `1/(96p)` of the
/// subgroup, not the `1/(24p)` the 2-torsion report field
/// `expected_rate` carried before this census.
#[derive(Clone, Debug, Default, Serialize)]
pub struct RateCensus {
    pub p: u64,
    pub seed: u64,
    pub n: u64,
    pub base: usize,
    pub quadruples: u64,
    /// Sign vectors whose sum lies in the subgroup for some `ε`: `16` on
    /// the quadruples of even class parity, `8` per quadruple on average.
    pub two_torsion_in_subgroup: u64,
    pub two_torsion_distinct: u64,
    pub two_torsion_predicted: u64,
    pub two_torsion_rate: f64,
    /// `rate · 192 p`: `1` when `|F| = p/4` exactly.
    pub two_torsion_rate_x_192p: f64,
    pub saturated_in_subgroup: u64,
    pub saturated_distinct: u64,
    pub saturated_predicted: u64,
    pub saturated_rate: f64,
    pub saturated_rate_x_96p: f64,
    pub group_ops: u64,
    pub wall_ms: f64,
}

fn fingerprint(pt: &PtE) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    pt.hash(&mut h);
    h.finish()
}

pub fn rate_census(p: u64, seed: u64) -> RateCensus {
    use std::collections::HashSet;
    let start = Instant::now();
    let inst = generate_instance_edwards(p, seed);
    let curve = &inst.curve;
    let n = inst.n;
    let mut rng = StdRng::seed_from_u64(seed ^ 0xED02);
    let base = factor_base(curve, &mut rng);
    let m = base.len();
    curve.reset_ops();
    let class: Vec<u8> = base.iter().map(|q| curve.torsion_class(n, q)).collect();
    // Pair sums P_i + s P_j, s = ±1, and their classes.
    let mut pair: HashMap<(usize, usize, i8), (PtE, u8)> = HashMap::new();
    for i in 0..m {
        for j in (i + 1)..m {
            pair.insert(
                (i, j, 1),
                (curve.add(&base[i], &base[j]), (class[i] + class[j]) % 4),
            );
            pair.insert(
                (i, j, -1),
                (curve.sub(&base[i], &base[j]), (class[i] + 4 - class[j]) % 4),
            );
        }
    }
    let mut two: HashSet<u64> = HashSet::new();
    let mut sat: HashSet<u64> = HashSet::new();
    let mut two_hits = 0u64;
    let mut sat_hits = 0u64;
    let mut quadruples = 0u64;
    for i in 0..m {
        for j in (i + 1)..m {
            for k in (j + 1)..m {
                for l in (k + 1)..m {
                    quadruples += 1;
                    for s in [1i8, -1] {
                        let (a, ca) = pair[&(i, j, s)];
                        for u in [1i8, -1] {
                            let (b, cb) = pair[&(k, l, u)];
                            for v in [1i8, -1] {
                                let (x, cx) = if v == 1 {
                                    (curve.add(&a, &b), (ca + cb) % 4)
                                } else {
                                    (curve.sub(&a, &b), (ca + 4 - cb) % 4)
                                };
                                // x and −x: the sixteen sign vectors.
                                for (pt, c) in [(x, cx), (curve.neg(&x), (4 - cx) % 4)] {
                                    // Saturated: the one t_q with c + t_q ≡ 0.
                                    let tq = (4 - c) % 4;
                                    let mut w = pt;
                                    for _ in 0..tq {
                                        w = PtE {
                                            x: w.y,
                                            y: curve.f.neg(&w.x),
                                        };
                                    }
                                    sat_hits += 1;
                                    sat.insert(fingerprint(&w));
                                    if c % 2 == 0 {
                                        two_hits += 1;
                                        two.insert(fingerprint(&w));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    RateCensus {
        p,
        seed,
        n,
        base: m,
        quadruples,
        two_torsion_in_subgroup: two_hits,
        two_torsion_distinct: two.len() as u64,
        two_torsion_predicted: 8 * quadruples,
        two_torsion_rate: two.len() as f64 / n as f64,
        two_torsion_rate_x_192p: two.len() as f64 / n as f64 * 192.0 * p as f64,
        saturated_in_subgroup: sat_hits,
        saturated_distinct: sat.len() as u64,
        saturated_predicted: 16 * quadruples,
        saturated_rate: sat.len() as f64 / n as f64,
        saturated_rate_x_96p: sat.len() as f64 / n as f64 * 96.0 * p as f64,
        group_ops: curve.ops(),
        wall_ms: start.elapsed().as_secs_f64() * 1e3,
    }
}

// ── Experiment: C″ with the residual saturated by Q₄ ─────────────────────

#[derive(Clone, Debug, Default, Serialize)]
pub struct CsecondEd4Report {
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
    /// `1/(96p)` on `|F| = p/4` classes: sixteen subgroup points per
    /// class-quadruple (`rate_census`).
    pub expected_rate: f64,
    pub planted_found: usize,
    pub planted_by_tq: [usize; 4],
    pub constructed_by_tq: [usize; 4],
    pub unverified: usize,
    pub mismatches: usize,
    pub undetermined: usize,
    pub timed_out: usize,
    /// Both runs together: the cost of one saturated test.
    pub c_second: CostStatsE,
    /// The run at `y_R` and the run at `x_R`, separately.
    pub c_second_at_y: CostStatsE,
    pub c_second_at_x: CostStatsE,
    pub weil_muls: CostStatsE,
    pub f4_muls: CostStatsE,
    pub roots_muls: CostStatsE,
    pub sign_group_ops: CostStatsE,
    pub f4_ms: CostStatsE,
    pub degree_reached: CostStatsE,
    pub max_rows: CostStatsE,
    pub max_cols: CostStatsE,
    pub c_second_decomposable: CostStatsE,
    pub f4_ms_decomposable: CostStatsE,
    pub degree_reached_decomposable: CostStatsE,
    pub per_test: Vec<[Jv4CostE; 2]>,
    pub mitm4_group_ops: f64,
    pub wall_ms: f64,
}

pub fn run_jv5_edwards4_csecond(
    p: u64,
    seed: u64,
    random: usize,
    constructed: usize,
    max_degree: u32,
    budget_secs: f64,
) -> CsecondEd4Report {
    let start = Instant::now();
    let inst = generate_instance_edwards(p, seed);
    let curve = &inst.curve;
    let f = &curve.f;
    let n = inst.n;
    let mut rng = StdRng::seed_from_u64(seed ^ 0xED04);
    let base = factor_base(curve, &mut rng);
    let by_y2: HashMap<u64, usize> = base
        .iter()
        .enumerate()
        .map(|(i, q)| (y2_key(curve, q), i))
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
    let pre = SymmetrisedS5Ed::precompute(curve, &mut rng);
    let precompute_ms = t0.elapsed().as_secs_f64() * 1e3;
    let table = pair_table(curve, &base);
    let mut rep = CsecondEd4Report {
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
        expected_rate: 1.0 / (96.0 * p as f64),
        ..Default::default()
    };
    let mut costs: Vec<[Jv4CostE; 2]> = Vec::new();
    let mut costs_dec: Vec<[Jv4CostE; 2]> = Vec::new();
    let mut mitm_ops = 0u64;
    let mut check_rng = StdRng::seed_from_u64(seed ^ 0x5EED);
    let mut check = |r: &PtE, want: Option<QuadE4>, rep: &mut CsecondEd4Report| -> [Jv4CostE; 2] {
        let opts = jv4_options(max_degree, budget_secs);
        let (found, cost) =
            jv4_decompose_saturated(curve, &base, &by_y2, &pre, r, &opts, &mut check_rng);
        rep.unverified += found
            .iter()
            .filter(|q| !verify_quad4(curve, &base, r, q))
            .count();
        let g0 = curve.ops();
        let oracle = mitm4_saturated(curve, &base, &table, r);
        mitm_ops += curve.ops() - g0;
        if cost[0].undetermined || cost[1].undetermined {
            rep.undetermined += 1;
        } else if found != oracle {
            rep.mismatches += 1;
        }
        if cost[0].timed_out || cost[1].timed_out {
            rep.timed_out += 1;
        }
        if let Some(w) = want {
            rep.constructed_by_tq[w.tq as usize] += 1;
            if found.contains(&w) {
                rep.planted_found += 1;
                rep.planted_by_tq[w.tq as usize] += 1;
            }
        }
        cost
    };
    for _ in 0..random {
        let a = rng.gen_range(0..n);
        let b = rng.gen_range(1..n);
        let r = curve.add(&curve.mul(&inst.g, a), &curve.mul(&inst.q, b));
        let c = check(&r, None, &mut rep);
        if c[0].decompositions + c[1].decompositions > 0 {
            rep.random_decomposable += 1;
        }
        costs.push(c);
    }
    for t in 0..constructed {
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
        // Every t_q in turn, so each translate is planted.
        let tq = (t % 4) as u8;
        let mut r = PtE::ID;
        for k in 0..4 {
            let q = if signs[k] == 1 {
                base[idx[k]]
            } else {
                curve.neg(&base[idx[k]])
            };
            r = curve.add(&r, &q);
        }
        for _ in 0..tq {
            r = curve.add(&r, &curve.q4());
        }
        let mut terms = [
            (idx[0], signs[0]),
            (idx[1], signs[1]),
            (idx[2], signs[2]),
            (idx[3], signs[3]),
        ];
        terms.sort_unstable();
        let c = check(&r, Some(QuadE4 { terms, tq }), &mut rep);
        costs_dec.push(c);
    }
    let one = |c: &Jv4CostE| {
        c.weil_muls + c.f4_muls + c.roots_muls + (c.group_ops as f64 * fp_per_add) as u64
    };
    let done = |v: &[[Jv4CostE; 2]]| -> Vec<[Jv4CostE; 2]> {
        v.iter()
            .filter(|c| !c[0].undetermined && !c[1].undetermined)
            .cloned()
            .collect()
    };
    let both = |v: &[[Jv4CostE; 2]], g: &dyn Fn(&Jv4CostE) -> u64| -> CostStatsE {
        stats(&v.iter().map(|c| g(&c[0]) + g(&c[1])).collect::<Vec<u64>>())
    };
    let pooled = |v: &[[Jv4CostE; 2]], g: &dyn Fn(&Jv4CostE) -> u64| -> CostStatsE {
        stats(
            &v.iter()
                .flat_map(|c| [g(&c[0]), g(&c[1])])
                .collect::<Vec<u64>>(),
        )
    };
    let cd = done(&costs);
    let dd = done(&costs_dec);
    rep.c_second = both(&cd, &one);
    rep.c_second_at_y = stats(&cd.iter().map(|c| one(&c[0])).collect::<Vec<u64>>());
    rep.c_second_at_x = stats(&cd.iter().map(|c| one(&c[1])).collect::<Vec<u64>>());
    rep.weil_muls = both(&cd, &|c| c.weil_muls);
    rep.f4_muls = both(&cd, &|c| c.f4_muls);
    rep.roots_muls = both(&cd, &|c| c.roots_muls);
    rep.sign_group_ops = both(&cd, &|c| c.group_ops);
    rep.f4_ms = both(&cd, &|c| c.f4_ms as u64);
    rep.degree_reached = pooled(&cd, &|c| c.degree_reached as u64);
    rep.max_rows = pooled(&cd, &|c| c.max_rows as u64);
    rep.max_cols = pooled(&cd, &|c| c.max_cols as u64);
    rep.c_second_decomposable = both(&dd, &one);
    rep.f4_ms_decomposable = both(&dd, &|c| c.f4_ms as u64);
    rep.degree_reached_decomposable = pooled(&dd, &|c| c.degree_reached as u64);
    rep.per_test = costs.iter().chain(costs_dec.iter()).cloned().collect();
    rep.mitm4_group_ops = mitm_ops as f64 / (random + constructed).max(1) as f64;
    rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
    rep
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edwards_group_law_and_torsion() {
        let mut rng = StdRng::seed_from_u64(7);
        let curve = Ed5::random(31, &mut rng);
        let p1 = random_point(&curve, &mut rng);
        let p2 = random_point(&curve, &mut rng);
        assert!(curve.on_curve(&p1) && curve.on_curve(&p2));
        let s = curve.add(&p1, &p2);
        assert!(curve.on_curve(&s));
        assert_eq!(curve.add(&p2, &p1), s);
        assert!(curve.sub(&s, &p2) == p1);
        let t = curve.t();
        assert!(curve.add(&t, &t).is_identity());
        let pt = curve.add(&p1, &t);
        assert_eq!(pt.x, curve.f.neg(&p1.x));
        assert_eq!(pt.y, curve.f.neg(&p1.y));
        // Associativity on a triple.
        let p3 = random_point(&curve, &mut rng);
        assert_eq!(
            curve.add(&curve.add(&p1, &p2), &p3),
            curve.add(&p1, &curve.add(&p2, &p3))
        );
    }

    #[test]
    fn s3_is_symmetric_and_vanishes_on_sums() {
        let mut rng = StdRng::seed_from_u64(8);
        let curve = Ed5::random(31, &mut rng);
        let f = &curve.f;
        let p1 = random_point(&curve, &mut rng);
        let p2 = random_point(&curve, &mut rng);
        let s = curve.add(&p1, &p2);
        let d = curve.sub(&p1, &p2);
        assert!(curve.s3(&[p1.y, p2.y, s.y]).is_zero());
        assert!(curve.s3(&[p1.y, p2.y, d.y]).is_zero());
        // Even flips preserve it (P₁ + T, (P₁ + P₂) + T); an odd flip does not.
        assert!(curve.s3(&[f.neg(&p1.y), p2.y, f.neg(&s.y)]).is_zero());
        let y: [E5; 3] = core::array::from_fn(|_| f.random(&mut rng));
        let v = curve.s3(&y);
        assert!(!v.is_zero());
        assert_eq!(curve.s3(&[y[1], y[2], y[0]]), v);
        assert_eq!(curve.s3(&[y[2], y[0], y[1]]), v);
        assert_eq!(curve.s3(&[y[1], y[0], y[2]]), v);
    }

    #[test]
    fn s5_vanishes_on_four_point_sums_and_the_symmetrisation_reproduces_it() {
        let mut rng = StdRng::seed_from_u64(9);
        let curve = Ed5::random(31, &mut rng);
        let f = &curve.f;
        let pts: Vec<PtE> = (0..4).map(|_| random_point(&curve, &mut rng)).collect();
        let sum = pts.iter().fold(PtE::ID, |acc, q| curve.add(&acc, q));
        assert!(curve
            .s5(&[pts[0].y, pts[1].y, pts[2].y, pts[3].y, sum.y])
            .is_zero());
        // With T on one point and on the sum.
        let st = curve.add(&sum, &curve.t());
        assert!(curve
            .s5(&[f.neg(&pts[0].y), pts[1].y, pts[2].y, pts[3].y, st.y])
            .is_zero());
        let pre = SymmetrisedS5Ed::precompute(&curve, &mut rng);
        assert_eq!(pre.monos_a.len(), 70);
        assert_eq!(pre.monos_b.len(), 35);
        let y = [pts[0].y, pts[1].y, pts[2].y, pts[3].y];
        assert!(pre.eval(f, &y, &sum.y).is_zero());
        let z: [E5; 4] = core::array::from_fn(|_| f.random(&mut rng));
        assert!(!pre.eval(f, &z, &sum.y).is_zero());
    }

    #[test]
    fn instances_base_and_oracle() {
        let inst = generate_instance_edwards(31, 3);
        let curve = &inst.curve;
        assert!(is_prime_u64(inst.n));
        assert!(curve.mul(&inst.g, inst.n).is_identity());
        let mut rng = StdRng::seed_from_u64(4);
        let base = factor_base(curve, &mut rng);
        assert!(base.len() >= 4, "{}", base.len());
        let table = pair_table(curve, &base);
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
            let mut r = idx
                .iter()
                .fold(PtE::ID, |acc, &i| curve.add(&acc, &base[i]));
            let eps = rng.gen_bool(0.5) as u8;
            if eps == 1 {
                r = curve.add(&r, &curve.t());
            }
            let found = mitm4_signed(curve, &base, &table, &r);
            let mut terms = [(idx[0], 1i64), (idx[1], 1), (idx[2], 1), (idx[3], 1)];
            terms.sort_unstable();
            let want = QuadE { terms, eps };
            assert!(found.contains(&want), "{want:?} not in {found:?}");
            assert!(found.iter().all(|q| verify_quad(curve, &base, &r, q)));
        }
    }

    #[test]
    fn q4_has_order_four_and_swaps_the_coordinates() {
        let mut rng = StdRng::seed_from_u64(11);
        let curve = Ed5::random(31, &mut rng);
        let f = &curve.f;
        let q4 = curve.q4();
        assert!(curve.on_curve(&q4));
        assert_eq!(curve.add(&q4, &q4), curve.t());
        assert!(curve.mul(&q4, 4).is_identity());
        let p1 = random_point(&curve, &mut rng);
        let plus = curve.add(&p1, &q4);
        assert_eq!(
            plus,
            PtE {
                x: p1.y,
                y: f.neg(&p1.x)
            }
        );
        let minus = curve.sub(&p1, &q4);
        assert_eq!(
            minus,
            PtE {
                x: f.neg(&p1.y),
                y: p1.x
            }
        );
        assert_eq!(curve.sub_q4(&p1), minus);
        // y(P − Q₄) = x_P: no function of y alone is Q₄-invariant.
        assert_eq!(minus.y, p1.x);
        // The class homomorphism on an instance.
        let inst = generate_instance_edwards(31, 3);
        let c = &inst.curve;
        assert_eq!(c.torsion_class(inst.n, &c.q4()), 1);
        assert_eq!(c.torsion_class(inst.n, &c.t()), 2);
        assert_eq!(c.torsion_class(inst.n, &inst.g), 0);
        let s = random_point(c, &mut rng);
        let cs = c.torsion_class(inst.n, &s);
        assert_eq!(c.torsion_class(inst.n, &c.add(&s, &c.q4())), (cs + 1) % 4);
        assert_eq!(c.torsion_class(inst.n, &c.neg(&s)), (4 - cs) % 4);
    }

    #[test]
    fn the_other_two_torsion_needs_a_square_root_of_d() {
        // d y₁²y₂² · S₃(−1/(δy₁), −1/(δy₂), y₃) = S₃(y₁, y₂, y₃) for δ² = d:
        // the degree-halving involution y ↦ −1/(δy) exists exactly when d is
        // a square, and keeps y ∈ F_p exactly when δ ∈ F_p.
        let mut rng = StdRng::seed_from_u64(12);
        let f = Fp5::new(31);
        for _ in 0..8 {
            let delta = loop {
                let v = f.random(&mut rng);
                if !v.is_zero() {
                    break v;
                }
            };
            let d = f.sq(&delta);
            let curve = Ed5::with_params(31, d);
            let y: [E5; 3] = loop {
                let y: [E5; 3] = core::array::from_fn(|_| f.random(&mut rng));
                if !y[0].is_zero() && !y[1].is_zero() {
                    break y;
                }
            };
            let flip = |v: &E5| f.neg(&f.inv(&f.mul(&delta, v)));
            let lhs = f.mul(
                &f.mul(&d, &f.mul(&f.sq(&y[0]), &f.sq(&y[1]))),
                &curve.s3(&[flip(&y[0]), flip(&y[1]), y[2]]),
            );
            assert_eq!(lhs, curve.s3(&y));
            // δ ∈ F_p only when d ∈ F_p: on the programme's curves (d ∉ F_p)
            // the involution leaves the factor base.
            let in_fp = f.from_fp(3);
            assert!(delta.in_fp() == flip(&in_fp).in_fp() || !delta.in_fp());
        }
    }

    #[test]
    fn saturation_finds_every_translate() {
        let inst = generate_instance_edwards(31, 3);
        let curve = &inst.curve;
        let mut rng = StdRng::seed_from_u64(13);
        let base = factor_base(curve, &mut rng);
        let table = pair_table(curve, &base);
        for tq in 0..4u8 {
            let idx: [usize; 4] = loop {
                let idx: [usize; 4] = core::array::from_fn(|_| rng.gen_range(0..base.len()));
                let mut s = idx.to_vec();
                s.sort_unstable();
                s.dedup();
                if s.len() == 4 {
                    break idx;
                }
            };
            let mut r = idx
                .iter()
                .fold(PtE::ID, |acc, &i| curve.add(&acc, &base[i]));
            for _ in 0..tq {
                r = curve.add(&r, &curve.q4());
            }
            let found = mitm4_saturated(curve, &base, &table, &r);
            let mut terms = [(idx[0], 1i64), (idx[1], 1), (idx[2], 1), (idx[3], 1)];
            terms.sort_unstable();
            let want = QuadE4 { terms, tq };
            assert!(found.contains(&want), "{want:?} not in {found:?}");
            assert!(found.iter().all(|q| verify_quad4(curve, &base, &r, q)));
            // The unsaturated oracle sees it only when t_q is even.
            let plain = mitm4_signed(curve, &base, &table, &r);
            assert_eq!(plain.contains(&QuadE { terms, eps: tq / 2 }), tq % 2 == 0);
        }
    }

    #[test]
    fn saturated_four_point_test_agrees_with_the_oracle() {
        let rep = run_jv5_edwards4_csecond(271, 1, 1, 4, 24, 120.0);
        assert_eq!(rep.planted_found, 4, "{:?}", rep.planted_by_tq);
        assert_eq!(rep.planted_by_tq, [1, 1, 1, 1]);
        assert_eq!(rep.mismatches, 0);
        assert_eq!(rep.unverified, 0);
        assert_eq!(rep.undetermined, 0);
    }

    #[test]
    fn census_counts_eight_and_sixteen_per_quadruple() {
        let c = rate_census(31, 3);
        assert!(c.quadruples > 0);
        // Collisions are O(1/p); at p = 31 allow a few.
        let two = c.two_torsion_distinct as f64 / c.two_torsion_predicted as f64;
        let sat = c.saturated_distinct as f64 / c.saturated_predicted as f64;
        assert!(two > 0.4 && two <= 2.0, "{c:?}");
        assert!(sat > 0.8 && sat <= 1.0, "{c:?}");
        // Even-parity quadruples contribute all sixteen sign vectors, odd
        // ones none; saturation takes every sign vector of every quadruple.
        assert_eq!(c.two_torsion_in_subgroup % 16, 0);
        assert!(c.two_torsion_in_subgroup <= 16 * c.quadruples);
        assert_eq!(c.saturated_in_subgroup, 16 * c.quadruples);
        assert_eq!(c.saturated_distinct % 16, 0);
    }
}
