//! (2^e, 2^e)-isogeny chains between products of elliptic curves through
//! Jacobians of genus-2 curves, in Mumford coordinates.
//!
//! Formulas: gluing, Richelot correspondence and splitting follow the
//! Castryck–Decru SageMath implementation (Pope, Oudompheng, Panny; MIT,
//! notice in ../LICENSE.castryck-decru),
//! re-derived here for fixed-degree polynomials so that every routine is a
//! straight-line program: no branch depends on field values.  Degenerate
//! inputs (which do not occur for honest data) produce garbage instead of
//! an early exit; FO re-encryption rejects the result.

use crate::ct::{Choice, ConditionallySelectable};
use crate::field::Field;
use crate::weier::{WCurve, WPoint};

/// Reduced divisor with monic u = x^2 + u1 x + u0, v = v1 x + v0.
#[derive(Clone, Copy, Debug)]
pub struct Div<F: Field> {
    pub u1: F,
    pub u0: F,
    pub v1: F,
    pub v0: F,
}

impl<F: Field> ConditionallySelectable for Div<F> {
    fn conditional_select(a: &Self, b: &Self, c: Choice) -> Self {
        Div {
            u1: F::conditional_select(&a.u1, &b.u1, c),
            u0: F::conditional_select(&a.u0, &b.u0, c),
            v1: F::conditional_select(&a.v1, &b.v1, c),
            v0: F::conditional_select(&a.v0, &b.v0, c),
        }
    }
}

/// Polynomial helpers on little-endian coefficient vectors (fixed lengths).
fn pmul<F: Field>(a: &[F], b: &[F]) -> Vec<F> {
    let mut r = vec![F::zero(); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            r[i + j] += a[i] * b[j];
        }
    }
    r
}

fn padd<F: Field>(a: &[F], b: &[F]) -> Vec<F> {
    let n = a.len().max(b.len());
    (0..n)
        .map(|i| {
            let x = if i < a.len() { a[i] } else { F::zero() };
            let y = if i < b.len() { b[i] } else { F::zero() };
            x + y
        })
        .collect()
}

fn pscale<F: Field>(a: &[F], s: F) -> Vec<F> {
    a.iter().map(|x| *x * s).collect()
}

fn pneg<F: Field>(a: &[F]) -> Vec<F> {
    a.iter().map(|x| -*x).collect()
}

/// a mod (x^2 + m1 x + m0) -> (r1, r0).
fn rem_quad<F: Field>(a: &[F], m1: F, m0: F) -> (F, F) {
    let mut c: Vec<F> = a.to_vec();
    let mut k = c.len();
    while k > 2 {
        k -= 1;
        let t = c[k];
        c[k] = F::zero();
        c[k - 1] -= t * m1;
        c[k - 2] -= t * m0;
    }
    let r1 = if c.len() > 1 { c[1] } else { F::zero() };
    (r1, c[0])
}

/// Quotient of a by the monic polynomial m (exact or not; remainder dropped).
fn quo_monic<F: Field>(a: &[F], m: &[F]) -> Vec<F> {
    let dm = m.len() - 1;
    let da = a.len() - 1;
    let mut c: Vec<F> = a.to_vec();
    let mut q = vec![F::zero(); da - dm + 1];
    for k in (0..=(da - dm)).rev() {
        let t = c[k + dm];
        q[k] = t;
        for j in 0..=dm {
            c[k + j] -= t * m[j];
        }
    }
    q
}

/// a mod m for monic m (deg m = 4): returns 4 coefficients.
fn rem_monic4<F: Field>(a: &[F], m: &[F; 5]) -> [F; 4] {
    let mut c: Vec<F> = a.to_vec();
    while c.len() < 5 {
        c.push(F::zero());
    }
    for k in (4..c.len()).rev() {
        let t = c[k];
        for j in 0..5 {
            c[k - 4 + j] -= t * m[j];
        }
    }
    [c[0], c[1], c[2], c[3]]
}

/// Inverse of (r1 x + r0) modulo the monic x^2 + p1 x + p0 (generic case).
fn inv_lin_mod_quad<F: Field>(r1: F, r0: F, p1: F, p0: F) -> (F, F) {
    let res = r0.square() - r0 * r1 * p1 + r1.square() * p0;
    let ri = res.invert();
    (-r1 * ri, (r0 - r1 * p1) * ri)
}

/// (a1 x + a0)(b1 x + b0) mod (x^2 + p1 x + p0).
fn mul_lin_mod_quad<F: Field>(a1: F, a0: F, b1: F, b0: F, p1: F, p0: F) -> (F, F) {
    let c2 = a1 * b1;
    let c1 = a1 * b0 + a0 * b1;
    let c0 = a0 * b0;
    (c1 - c2 * p1, c0 - c2 * p0)
}

/// 4x4 determinant (row-major).
fn det4<F: Field>(m: &[[F; 4]; 4]) -> F {
    let mut d = F::zero();
    for c in 0..4 {
        let minor = minor3(m, 0, c);
        let s = if c % 2 == 0 { m[0][c] * minor } else { -(m[0][c] * minor) };
        d += s;
    }
    d
}

fn minor3<F: Field>(m: &[[F; 4]; 4], r: usize, c: usize) -> F {
    let rows: Vec<usize> = (0..4).filter(|&i| i != r).collect();
    let cols: Vec<usize> = (0..4).filter(|&j| j != c).collect();
    let a = |i: usize, j: usize| m[rows[i]][cols[j]];
    a(0, 0) * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) - a(0, 1) * (a(1, 0) * a(2, 2) - a(1, 2) * a(2, 0))
        + a(0, 2) * (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0))
}

/// Inverse of a (deg <= 3) modulo the monic quartic m, via Cramer's rule on
/// the multiplication matrix.  Constant time.
fn inv_mod_quartic<F: Field>(a: &[F; 4], m: &[F; 5]) -> [F; 4] {
    // columns: a * x^j mod m
    let mut cols = [[F::zero(); 4]; 4];
    let mut cur: Vec<F> = a.to_vec();
    for j in 0..4 {
        let r = rem_monic4(&cur, m);
        cols[j] = r;
        let mut nxt = vec![F::zero()];
        nxt.extend_from_slice(&r);
        cur = nxt;
    }
    // matrix M[i][j] = coefficient i of column j; solve M c = e0
    let mut mm = [[F::zero(); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            mm[i][j] = cols[j][i];
        }
    }
    let det = det4(&mm);
    let di = det.invert();
    // c_j = cofactor(0, j)/det  (column of the adjugate)  -> adj[j][0] = C[0][j]
    let mut c = [F::zero(); 4];
    for j in 0..4 {
        let minor = minor3(&mm, 0, j);
        let cof = if j % 2 == 0 { minor } else { -minor };
        c[j] = cof * di;
    }
    c
}

// ---------------------------------------------------------------------------
// Jacobian arithmetic on y^2 = h(x), deg h = 6
// ---------------------------------------------------------------------------

/// Generic Cantor addition of divisors with coprime u's.
pub fn div_add<F: Field>(h: &[F; 7], a: &Div<F>, b: &Div<F>) -> Div<F> {
    // w = (vb - va) * (ua mod ub)^-1 mod ub
    let (r1, r0) = (a.u1 - b.u1, a.u0 - b.u0);
    let (i1, i0) = inv_lin_mod_quad(r1, r0, b.u1, b.u0);
    let (w1, w0) = mul_lin_mod_quad(b.v1 - a.v1, b.v0 - a.v0, i1, i0, b.u1, b.u0);
    // v3 = va + ua * w  (cubic)
    let ua = [a.u0, a.u1, F::one()];
    let uw = pmul(&ua, &[w0, w1]);
    let v3 = padd(&[a.v0, a.v1], &uw);
    let ub = [b.u0, b.u1, F::one()];
    let uu = pmul(&ua, &ub);
    reduce_from(h, &v3, &uu)
}

/// u'' = (h - v3^2) / U (U monic quartic), monic; v'' = -v3 mod u''.
fn reduce_from<F: Field>(h: &[F; 7], v3: &[F], uu: &[F]) -> Div<F> {
    let v3sq = pmul(v3, v3);
    let num = padd(h, &pneg(&v3sq));
    let num = &num[..7];
    let q = quo_monic(num, uu); // degree 2
    let li = q[2].invert();
    let u1 = q[1] * li;
    let u0 = q[0] * li;
    let (r1, r0) = rem_quad(v3, u1, u0);
    Div { u1, u0, v1: -r1, v0: -r0 }
}

/// Generic doubling.
pub fn div_dbl<F: Field>(h: &[F; 7], a: &Div<F>) -> Div<F> {
    let u = [a.u0, a.u1, F::one()];
    let v = [a.v0, a.v1];
    let vsq = pmul(&v, &v);
    let num = padd(h, &pneg(&vsq));
    let t = quo_monic(&num, &u); // (h - v^2)/u, degree 4
    let (t1, t0) = rem_quad(&t, a.u1, a.u0);
    let (i1, i0) = inv_lin_mod_quad(a.v1.double(), a.v0.double(), a.u1, a.u0);
    let (k1, k0) = mul_lin_mod_quad(t1, t0, i1, i0, a.u1, a.u0);
    let uk = pmul(&u, &[k0, k1]);
    let v3 = padd(&v, &uk);
    let uu = pmul(&u, &u);
    reduce_from(h, &v3, &uu)
}

pub fn div_dbl_iter<F: Field>(h: &[F; 7], a: &Div<F>, n: usize) -> Div<F> {
    let mut r = *a;
    for _ in 0..n {
        r = div_dbl(h, &r);
    }
    r
}

// ---------------------------------------------------------------------------
// Gluing  E1 x E2 -> Jac(H)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
pub struct Glue<F: Field> {
    pub h: [F; 7],
    #[allow(dead_code)] // kept for debugging the gluing map
    s1: F,
    s2: F,
    t1: F,
    t2: F,
    s1i: F,
}

impl<F: Field> Glue<F> {
    /// `a` = x-coordinates of the kernel's 2-torsion on E1: x(P1), x(Q1),
    /// x(P1+Q1); `b` the same on E2.
    pub fn new(a: [F; 3], b: [F; 3]) -> Self {
        let m = [[a[0] * b[0], a[0], b[0]], [a[1] * b[1], a[1], b[1]], [a[2] * b[2], a[2], b[2]]];
        // adjugate and determinant of the 3x3 matrix
        let cof = |i: usize, j: usize| {
            let r: Vec<usize> = (0..3).filter(|&x| x != i).collect();
            let c: Vec<usize> = (0..3).filter(|&x| x != j).collect();
            let v = m[r[0]][c[0]] * m[r[1]][c[1]] - m[r[0]][c[1]] * m[r[1]][c[0]];
            if (i + j) % 2 == 0 {
                v
            } else {
                -v
            }
        };
        let det = m[0][0] * cof(0, 0) + m[0][1] * cof(0, 1) + m[0][2] * cof(0, 2);
        // (R, S, T) = M^-1 (1,1,1) = adj(M) (1,1,1) / det, adj = cof^T
        let adj_row = |k: usize| cof(0, k) + cof(1, k) + cof(2, k); // row k of adj = column k of cof
        let rr_det = adj_row(0); // R * det
        let ss_det = adj_row(1);
        let tt_det = adj_row(2);
        let da = (a[0] - a[1]) * (a[1] - a[2]) * (a[2] - a[0]);
        let db = (b[0] - b[1]) * (b[1] - b[2]) * (b[2] - b[0]);
        // RD = R * det = rr_det; s1 = -da / RD, t1 = db / RD, s2 = -T/R, t2 = -S/R

        let rdi = rr_det.invert();
        let s1 = -da * rdi;
        let t1 = db * rdi;
        let s2 = -tt_det * rdi;
        let t2 = -ss_det * rdi;
        let _ = det;
        let s1i = s1.invert();
        let at: Vec<F> = a.iter().map(|x| (*x - s2) * s1i).collect();
        // h = s1 (x^2 - at0)(x^2 - at1)(x^2 - at2)
        let q0 = [-at[0], F::zero(), F::one()];
        let q1 = [-at[1], F::zero(), F::one()];
        let q2 = [-at[2], F::zero(), F::one()];
        let hh = pscale(&pmul(&pmul(&q0, &q1), &q2), s1);
        let mut h = [F::zero(); 7];
        h.copy_from_slice(&hh[..7]);
        Glue { h, s1, s2, t1, t2, s1i }
    }

    fn map_e1(&self, p: &WPoint<F>) -> Div<F> {
        // u = x^2 + (s2 - xP)/s1, v = yP/s1
        Div { u1: F::zero(), u0: (self.s2 - p.x) * self.s1i, v1: F::zero(), v0: p.y * self.s1i }
    }
    fn map_e2(&self, p: &WPoint<F>) -> Div<F> {
        // u = x^2 - t1/(xP - t2), v = yP/(xP - t2) * x
        let di = (p.x - self.t2).invert();
        Div { u1: F::zero(), u0: -(self.t1 * di), v1: p.y * di, v0: F::zero() }
    }
    /// Image of (P1, P2); either component may be infinity (selected in
    /// constant time).
    pub fn map(&self, p1: &WPoint<F>, p2: &WPoint<F>) -> Div<F> {
        let d1 = self.map_e1(p1);
        let d2 = self.map_e2(p2);
        let s = div_add(&self.h, &d1, &d2);
        let s = Div::conditional_select(&s, &d2, p1.inf);
        Div::conditional_select(&s, &d1, p2.inf & !p1.inf)
    }
}

// ---------------------------------------------------------------------------
// Richelot isogeny Jac(H) -> Jac(H')
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct Richelot<F: Field> {
    g1: [F; 3],
    g2: [F; 3],
    h1: Vec<F>,
    h11: Vec<F>,
    h12: Vec<F>,
    h22: Vec<F>,
    pub hnew: [F; 7],
}

fn inv3x3<F: Field>(m: &[[F; 3]; 3]) -> [[F; 3]; 3] {
    let c = |i: usize, j: usize| {
        let r: Vec<usize> = (0..3).filter(|&x| x != i).collect();
        let k: Vec<usize> = (0..3).filter(|&x| x != j).collect();
        let v = m[r[0]][k[0]] * m[r[1]][k[1]] - m[r[0]][k[1]] * m[r[1]][k[0]];
        if (i + j) % 2 == 0 {
            v
        } else {
            -v
        }
    };
    let det = m[0][0] * c(0, 0) + m[0][1] * c(0, 1) + m[0][2] * c(0, 2);
    let di = det.invert();
    let mut r = [[F::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = c(j, i) * di;
        }
    }
    r
}

fn kernel_quadratics<F: Field>(h: &[F; 7], d1: &Div<F>, d2: &Div<F>) -> ([F; 3], [F; 3], [F; 3]) {
    let g1 = [d1.u0, d1.u1, F::one()];
    let g2 = [d2.u0, d2.u1, F::one()];
    let g12 = pmul(&g1, &g2);
    let g3v = quo_monic(&h[..], &g12);
    let g3 = [g3v[0], g3v[1], g3v[2]];
    (g1, g2, g3)
}

impl<F: Field> Richelot<F> {
    /// `d1`, `d2`: the kernel's two generators (order 2).
    pub fn new(h: &[F; 7], d1: &Div<F>, d2: &Div<F>) -> Self {
        let (g1, g2, g3) = kernel_quadratics(h, d1, d2);
        let delta = [g1, g2, g3];
        let dinv = inv3x3(&delta);
        let mk = |col: usize| vec![-dinv[2][col], dinv[1][col].double(), -dinv[0][col]];
        let h1 = mk(0);
        let h2 = mk(1);
        let h3 = mk(2);
        let hn = pmul(&pmul(&h1, &h2), &h3);
        let mut hnew = [F::zero(); 7];
        hnew.copy_from_slice(&hn[..7]);
        let h11 = pmul(&h1, &h1);
        let h12 = pmul(&h1, &h2);
        let h22 = pmul(&h2, &h2);
        Richelot { g1, g2, h1, h11, h12, h22, hnew }
    }

    pub fn map(&self, d: &Div<F>) -> Div<F> {
        let s = -d.u1;
        let p = d.u0;
        let (g11, g10) = (self.g1[1] - d.u1, self.g1[0] - d.u0);
        let (g21, g20) = (self.g2[1] - d.u1, self.g2[0] - d.u0);
        let (v1, v0) = (d.v1, d.v0);
        let two = F::from_i64(2);
        let c11 = g11.square() * p + g11 * g10 * s + g10.square();
        let c12 = two * g11 * g21 * p + (g11 * g20 + g21 * g10) * s + two * g10 * g20;
        let c22 = g21.square() * p + g21 * g20 * s + g20.square();
        let px = padd(&padd(&pscale(&self.h11, c11), &pscale(&self.h12, c12)), &pscale(&self.h22, c22));
        let py2 = v1.square() * p + v1 * v0 * s + v0.square();
        let a = two * v1 * g11 * p + v1 * g10 * s + v0 * g11 * s + two * v0 * g10;
        let b = v1 * g11 * s * p + two * v1 * g10 * p + v0 * g11 * (s.square() - two * p) + v0 * g10 * s;
        let py1 = pmul(&[-b, a], &self.h1); // degree 3
        let u = [d.u0, d.u1, F::one()];
        let py0 = pscale(&pmul(&self.h11, &u), c11); // degree 6
        // make Px monic
        let lc_i = px[4].invert();
        let pxm = [px[0] * lc_i, px[1] * lc_i, px[2] * lc_i, px[3] * lc_i, F::one()];
        let py1r = rem_monic4(&py1, &pxm);
        let inv = inv_mod_quartic(&py1r, &pxm);
        let rhs = padd(&pscale(&self.hnew, py2), &py0);
        let rhsr = rem_monic4(&rhs, &pxm);
        let prod = pmul(&inv, &rhsr);
        let pyr = rem_monic4(&prod, &pxm);
        let py: Vec<F> = pyr.iter().map(|x| -*x).collect();
        // Dx = (hnew - Py^2) / Pxm, made monic; Dy = -Py mod Dx
        let pysq = pmul(&py, &py);
        let num = padd(&self.hnew, &pneg(&pysq));
        let q = quo_monic(&num[..7], &pxm);
        let li = q[2].invert();
        let u1 = q[1] * li;
        let u0 = q[0] * li;
        let (r1, r0) = rem_quad(&py, u1, u0);
        Div { u1, u0, v1: -r1, v0: -r0 }
    }
}

// ---------------------------------------------------------------------------
// Splitting Jac(H) -> E1 x E2
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct Split<F: Field> {
    a: F,
    b: F,
    d: F,
    p1: Vec<F>,
    p2: Vec<F>,
    k1: F,
    k2: F,
    pub e1: WCurve<F>,
    pub e2: WCurve<F>,
}

impl<F: Field> Split<F> {
    pub fn new(h: &[F; 7], d1: &Div<F>, d2: &Div<F>) -> Self {
        let (g1, g2, g3) = kernel_quadratics(h, d1, d2);
        let m = [g1, g2, g3];
        // right kernel vector of the singular matrix: cross product of rows 0, 1
        let r0 = m[0];
        let r1 = m[1];
        let u = r0[1] * r1[2] - r0[2] * r1[1];
        let v = r0[2] * r1[0] - r0[0] * r1[2];
        let w = r0[0] * r1[1] - r0[1] * r1[0];
        let d = u * F::from_i64(2).invert();
        // roots of X^2 - v X + u w / 4
        let disc = v.square() - u * w;
        let (sq, _) = disc.sqrt();
        let half = F::from_i64(2).invert();
        let ad = (v + sq) * half;
        let b = (v - sq) * half;
        let a = ad * d.invert();
        let ev = |g: &[F; 3], x0: F, x1: F, x2: F| g[0] * x0 + g[1] * x1 + g[2] * x2;
        let (h11, h21, h31) = (ev(&g1, F::one(), a, a.square()), ev(&g2, F::one(), a, a.square()), ev(&g3, F::one(), a, a.square()));
        let (h10, h20, h30) = (ev(&g1, d.square(), b * d, b.square()), ev(&g2, d.square(), b * d, b.square()), ev(&g3, d.square(), b * d, b.square()));
        let p1 = pmul(&pmul(&[h10, h11], &[h20, h21]), &[h30, h31]);
        let p2 = pmul(&pmul(&[h11, h10], &[h21, h20]), &[h31, h30]);
        let p1n = pmul(&pmul(&[h10 * h21 * h31, F::one()], &[h20 * h11 * h31, F::one()]), &[h30 * h11 * h21, F::one()]);
        let p2n = pmul(&pmul(&[h11 * h20 * h30, F::one()], &[h21 * h10 * h30, F::one()]), &[h31 * h10 * h20, F::one()]);
        let e1 = WCurve::new(p1n[2], p1n[1], p1n[0]);
        let e2 = WCurve::new(p2n[2], p2n[1], p2n[0]);
        Split { a, b, d, p1, p2, k1: h11 * h21 * h31, k2: h10 * h20 * h30, e1, e2 }
    }

    pub fn map(&self, dv: &Div<F>) -> (WPoint<F>, WPoint<F>) {
        let (a, b, d) = (self.a, self.b, self.d);
        let xd = [d, F::one()]; // x + d
        let axb = [b, a]; // a x + b
        let xd2 = pmul(&xd, &xd);
        let u_ = padd(&padd(&pscale(&xd2, dv.u0), &pscale(&pmul(&axb, &xd), dv.u1)), &pmul(&axb, &axb));
        let vv = padd(&pscale(&pmul(&xd2, &xd), dv.v0), &pscale(&pmul(&axb, &xd2), dv.v1));
        // reduce vv mod u_ (u_ has leading coefficient u_[2])
        let lci = u_[2].invert();
        let m1 = u_[1] * lci;
        let m0 = u_[0] * lci;
        let (v1, v0) = rem_quad(&vv, m1, m0);
        let s = -m1;
        let p = m0;
        let two = F::from_i64(2);
        // E1 side
        let u1p1 = -(s.square() - two * p);
        let u1p0 = p.square();
        let num = padd(&self.p1, &[v0.square(), -v1.square()]);
        let num = pscale(&num, (two * v0).invert());
        let (w1, w0) = rem_quad(&num, u1p1, u1p0);
        let wsq = pmul(&[w0, w1], &[w0, w1]);
        let red = quo_monic(&padd(&self.p1, &pneg(&wsq)), &[u1p0, u1p1, F::one()]); // linear
        let xp1 = -red[0] * red[1].invert();
        let yp1 = w1 * xp1 + w0;
        // E2 side
        let pi = p.invert();
        let pi2 = pi.square();
        let u2p1 = -(s.square() - two * p) * pi2;
        let u2p0 = pi2;
        let c21 = v1 * (s.square() - two * p) + v0 * s;
        let (q1, q0) = rem_quad(&[F::zero(), F::zero(), c21], u2p1, u2p0);
        let c4 = p * (v1.square() * p + v1 * v0 * s + v0.square());
        let v20 = padd(&self.p2, &[F::zero(), F::zero(), F::zero(), F::zero(), c4]);
        let (r1, r0) = rem_quad(&v20, u2p1, u2p0);
        let (i1, i0) = inv_lin_mod_quad(q1, q0, u2p1, u2p0);
        let (z1, z0) = mul_lin_mod_quad(i1, i0, r1, r0, u2p1, u2p0);
        let zsq = pmul(&[z0, z1], &[z0, z1]);
        let red2 = quo_monic(&padd(&self.p2, &pneg(&zsq)), &[u2p0, u2p1, F::one()]);
        let xp2 = -red2[0] * red2[1].invert();
        let yp2 = z1 * xp2 + z0;
        (WPoint::new(self.k1 * xp1, self.k1 * yp1), WPoint::new(self.k2 * xp2, self.k2 * yp2))
    }
}

/// Check v^2 = h mod u (testing aid).
pub fn div_on_curve<F: Field>(h: &[F; 7], d: &Div<F>) -> bool {
    let v = [d.v0, d.v1];
    let num = padd(h, &pneg(&pmul(&v, &v)));
    let (r1, r0) = rem_quad(&num, d.u1, d.u0);
    bool::from(r1.is_zero() & r0.is_zero())
}

pub(crate) fn _unused<F: Field>(x: F) -> Choice {
    x.ct_eq(&x)
}

// ---------------------------------------------------------------------------
// Batched versions (Montgomery's trick across independent divisors)
// ---------------------------------------------------------------------------

/// Batch inversion over any Field; zeros map to zero.  Constant time.
pub fn batch_inv<F: Field>(xs: &[F]) -> Vec<F> {
    let n = xs.len();
    if n == 0 {
        return vec![];
    }
    let mut pre = Vec::with_capacity(n);
    let mut zs = Vec::with_capacity(n);
    let mut acc = F::one();
    for x in xs {
        let z = x.is_zero();
        zs.push(z);
        pre.push(acc);
        acc *= F::conditional_select(x, &F::one(), z);
    }
    let mut inv = acc.invert();
    let mut out = vec![F::zero(); n];
    for k in (0..n).rev() {
        let xx = F::conditional_select(&xs[k], &F::one(), zs[k]);
        out[k] = F::conditional_select(&(inv * pre[k]), &F::zero(), zs[k]);
        inv *= xx;
    }
    out
}

/// Doubling of several divisors in lockstep (2 batch inversions in total).
pub fn div_dbl_batch<F: Field>(h: &[F; 7], ds: &[Div<F>]) -> Vec<Div<F>> {
    // stage 1: t mod u and the resultant of 2v mod u
    let mut st = Vec::with_capacity(ds.len());
    let mut res = Vec::with_capacity(ds.len());
    for a in ds {
        let u = [a.u0, a.u1, F::one()];
        let v = [a.v0, a.v1];
        let num = padd(h, &pneg(&pmul(&v, &v)));
        let t = quo_monic(&num, &u);
        let (t1, t0) = rem_quad(&t, a.u1, a.u0);
        let (r1, r0) = (a.v1.double(), a.v0.double());
        res.push(r0.square() - r0 * r1 * a.u1 + r1.square() * a.u0);
        st.push((t1, t0, r1, r0));
    }
    let ri = batch_inv(&res);
    // stage 2: v3 and the quadratic quotient
    let mut qs = Vec::with_capacity(ds.len());
    let mut lcs = Vec::with_capacity(ds.len());
    for (k, a) in ds.iter().enumerate() {
        let (t1, t0, r1, r0) = st[k];
        let (i1, i0) = (-r1 * ri[k], (r0 - r1 * a.u1) * ri[k]);
        let (k1, k0) = mul_lin_mod_quad(t1, t0, i1, i0, a.u1, a.u0);
        let u = [a.u0, a.u1, F::one()];
        let v3 = padd(&[a.v0, a.v1], &pmul(&u, &[k0, k1]));
        let uu = pmul(&u, &u);
        let num = padd(h, &pneg(&pmul(&v3, &v3)));
        let q = quo_monic(&num[..7], &uu);
        lcs.push(q[2]);
        qs.push((q, v3));
    }
    let li = batch_inv(&lcs);
    qs.iter()
        .enumerate()
        .map(|(k, (q, v3))| {
            let u1 = q[1] * li[k];
            let u0 = q[0] * li[k];
            let (r1, r0) = rem_quad(v3, u1, u0);
            Div { u1, u0, v1: -r1, v0: -r0 }
        })
        .collect()
}

pub fn div_dbl_iter_batch<F: Field>(h: &[F; 7], ds: &[Div<F>], n: usize) -> Vec<Div<F>> {
    let mut r = ds.to_vec();
    for _ in 0..n {
        r = div_dbl_batch(h, &r);
    }
    r
}

impl<F: Field> Richelot<F> {
    /// Batched map (3 batch inversions for the whole slice).
    pub fn map_batch(&self, ds: &[Div<F>]) -> Vec<Div<F>> {
        let two = F::from_i64(2);
        struct S1<F> {
            px: Vec<F>,
            py2: F,
            py1: Vec<F>,
            py0: Vec<F>,
        }
        let mut s1 = Vec::with_capacity(ds.len());
        let mut lcs = Vec::with_capacity(ds.len());
        for d in ds {
            let s = -d.u1;
            let p = d.u0;
            let (g11, g10) = (self.g1[1] - d.u1, self.g1[0] - d.u0);
            let (g21, g20) = (self.g2[1] - d.u1, self.g2[0] - d.u0);
            let (v1, v0) = (d.v1, d.v0);
            let c11 = g11.square() * p + g11 * g10 * s + g10.square();
            let c12 = two * g11 * g21 * p + (g11 * g20 + g21 * g10) * s + two * g10 * g20;
            let c22 = g21.square() * p + g21 * g20 * s + g20.square();
            let px = padd(&padd(&pscale(&self.h11, c11), &pscale(&self.h12, c12)), &pscale(&self.h22, c22));
            let py2 = v1.square() * p + v1 * v0 * s + v0.square();
            let a = two * v1 * g11 * p + v1 * g10 * s + v0 * g11 * s + two * v0 * g10;
            let b = v1 * g11 * s * p + two * v1 * g10 * p + v0 * g11 * (s.square() - two * p) + v0 * g10 * s;
            let py1 = pmul(&[-b, a], &self.h1);
            let u = [d.u0, d.u1, F::one()];
            let py0 = pscale(&pmul(&self.h11, &u), c11);
            lcs.push(px[4]);
            s1.push(S1 { px, py2, py1, py0 });
        }
        let lci = batch_inv(&lcs);
        // stage 2: monic Px, multiplication matrix and its determinant
        let mut s2 = Vec::with_capacity(ds.len());
        let mut dets = Vec::with_capacity(ds.len());
        for (k, x) in s1.iter().enumerate() {
            let li = lci[k];
            let pxm = [x.px[0] * li, x.px[1] * li, x.px[2] * li, x.px[3] * li, F::one()];
            let py1r = rem_monic4(&x.py1, &pxm);
            let mut cols = [[F::zero(); 4]; 4];
            let mut cur: Vec<F> = py1r.to_vec();
            for j in 0..4 {
                let r = rem_monic4(&cur, &pxm);
                cols[j] = r;
                let mut nxt = vec![F::zero()];
                nxt.extend_from_slice(&r);
                cur = nxt;
            }
            let mut mm = [[F::zero(); 4]; 4];
            for i in 0..4 {
                for j in 0..4 {
                    mm[i][j] = cols[j][i];
                }
            }
            dets.push(det4(&mm));
            s2.push((pxm, mm));
        }
        let deti = batch_inv(&dets);
        // stage 3: Py and the new u (before normalisation)
        let mut s3 = Vec::with_capacity(ds.len());
        let mut q2s = Vec::with_capacity(ds.len());
        for k in 0..ds.len() {
            let (pxm, mm) = &s2[k];
            let mut inv = [F::zero(); 4];
            for j in 0..4 {
                let minor = minor3(mm, 0, j);
                let cof = if j % 2 == 0 { minor } else { -minor };
                inv[j] = cof * deti[k];
            }
            let rhs = padd(&pscale(&self.hnew, s1[k].py2), &s1[k].py0);
            let rhsr = rem_monic4(&rhs, pxm);
            let pyr = rem_monic4(&pmul(&inv, &rhsr), pxm);
            let py: Vec<F> = pyr.iter().map(|x| -*x).collect();
            let num = padd(&self.hnew, &pneg(&pmul(&py, &py)));
            let q = quo_monic(&num[..7], pxm);
            q2s.push(q[2]);
            s3.push((q, py));
        }
        let q2i = batch_inv(&q2s);
        s3.iter()
            .enumerate()
            .map(|(k, (q, py))| {
                let u1 = q[1] * q2i[k];
                let u0 = q[0] * q2i[k];
                let (r1, r0) = rem_quad(py, u1, u0);
                Div { u1, u0, v1: -r1, v0: -r0 }
            })
            .collect()
    }
}
