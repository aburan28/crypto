//! Genus-2 curves C: y^2 = f(x) (deg f = 5 or 6, odd characteristic) and (2,2)-isogenies of their
//! Jacobians: Richelot's construction, its degenerate case (codomain E1 x E2, "splitting"), the
//! inverse "gluing" E1 x E2 -> J(C) (Howe–Leprévost–Poonen), Igusa–Clebsch invariants, and the
//! superspecial Richelot graph over F_{p^2} (Jacobians and products of supersingular curves).
//! Correctness is checked through L-polynomials: isogenous abelian surfaces have equal
//! L-polynomials, and L_{E1 x E2} = L_{E1} L_{E2}.
//!
//! Richelot: f = G1 G2 G3 with G_i = g_i2 x^2 + g_i1 x + g_i0 (a degree-1 G_i means a root at
//! infinity). Delta = det(g_ij); H1 = G2' G3 - G2 G3', H2 = G3' G1 - G3 G1', H3 = G1' G2 - G1 G2'.
//! If Delta != 0 the codomain is the Jacobian of y^2 = Delta H1 H2 H3; if Delta = 0 there is an
//! involution of P^1 swapping the roots of every G_i, and in a coordinate z where it is z -> -z,
//! (1-z)^6 f = c6 z^6 + c4 z^4 + c2 z^2 + c0 and J(C) ~ E1 x E2 with
//! E1: y^2 = c6 u^3 + c4 u^2 + c2 u + c0, E2: y^2 = c0 u^3 + c2 u^2 + c4 u + c6.
use crate::curve::Curve;
use crate::field::{Field, Rng, Zp, Zp2};
use crate::poly::{self, Poly};
use std::collections::{HashMap, HashSet};

/// L-polynomial data of a genus-2 curve over F_p: L(T) = 1 + c1 T + c2 T^2 + p c1 T^3 + p^2 T^4.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LPoly {
    pub c1: i64,
    pub c2: i64,
}

impl LPoly {
    /// #J(F_p) = L(1).
    pub fn jac_order(&self, p: i64) -> i64 {
        1 + self.c1 + self.c2 + p * self.c1 + p * p
    }
    /// L of E1 x E2 from the traces a1, a2: (1 - a1 T + p T^2)(1 - a2 T + p T^2).
    pub fn product(a1: i64, a2: i64, p: i64) -> LPoly {
        LPoly {
            c1: -(a1 + a2),
            c2: 2 * p + a1 * a2,
        }
    }
}

/// Quadratic character table of F_p (p small).
fn legendre_table(p: u64) -> Vec<i8> {
    let mut t = vec![-1i8; p as usize];
    t[0] = 0;
    for x in 1..p {
        t[((x * x) % p) as usize] = 1;
    }
    t
}

/// #C(F_p) and #C(F_{p^2}) for C: y^2 = f(x) over F_p (naive, O(p^2); p up to ~10^4).
pub fn point_counts(fp: &Zp, f: &Poly<Zp>) -> (i64, i64) {
    let p = fp.p;
    let leg = legendre_table(p);
    let d = f.len() - 1;
    let lc = f[d];
    let inf1 = if d == 5 {
        1
    } else {
        1 + leg[lc as usize] as i64
    };
    let mut n1 = inf1;
    for x in 0..p {
        n1 += 1 + leg[poly::eval(fp, f, x) as usize] as i64;
    }
    // over F_{p^2}: chi(a) = legendre(Norm(a)); lc in F_p is a square in F_{p^2}
    let f2 = Zp2::new(p);
    let fl: Vec<(u64, u64)> = f.iter().map(|&c| (c, 0)).collect();
    let inf2 = if d == 5 { 1 } else { 2 };
    let mut n2 = inf2;
    for a in 0..p {
        for b in 0..p {
            let v = poly::eval(&f2, &fl, (a, b));
            let nrm = fp.sub(fp.mul(v.0, v.0), fp.mul(f2.nr, fp.mul(v.1, v.1)));
            n2 += 1 + leg[nrm as usize] as i64;
        }
    }
    (n1, n2)
}

pub fn lpoly(fp: &Zp, f: &Poly<Zp>) -> LPoly {
    let p = fp.p as i64;
    let (n1, n2) = point_counts(fp, f);
    let c1 = n1 - p - 1;
    // N2 = p^2 + 1 + 2 c2 - c1^2
    let c2 = (n2 - p * p - 1 + c1 * c1) / 2;
    LPoly { c1, c2 }
}

/// Trace of Frobenius of y^2 = cubic(u) over F_p (naive).
pub fn elliptic_trace(fp: &Zp, cubic: &Poly<Zp>) -> i64 {
    let p = fp.p;
    let leg = legendre_table(p);
    let mut n = 1i64;
    for x in 0..p {
        n += 1 + leg[poly::eval(fp, cubic, x) as usize] as i64;
    }
    p as i64 + 1 - n
}

pub struct Richelot<F: Field> {
    pub delta: F::E,
    pub h: [Poly<F>; 3],
}

fn coeffs<F: Field>(f: &F, g: &Poly<F>) -> [F::E; 3] {
    [0, 1, 2].map(|k| if k < g.len() { g[k] } else { f.zero() })
}

fn det3<F: Field>(f: &F, m: &[[F::E; 3]; 3]) -> F::E {
    let t = |a, b, c, d| f.sub(f.mul(a, d), f.mul(b, c));
    let x = f.mul(m[0][0], t(m[1][1], m[1][2], m[2][1], m[2][2]));
    let y = f.mul(m[0][1], t(m[1][0], m[1][2], m[2][0], m[2][2]));
    let z = f.mul(m[0][2], t(m[1][0], m[1][1], m[2][0], m[2][1]));
    f.add(f.sub(x, y), z)
}

fn bracket<F: Field>(f: &F, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
    poly::sub(
        f,
        &poly::mul(f, &poly::derivative(f, a), b),
        &poly::mul(f, a, &poly::derivative(f, b)),
    )
}

/// Richelot data for the (2,2)-subgroup given by the quadratic splitting f = G1 G2 G3.
pub fn richelot<F: Field>(f: &F, g: &[Poly<F>; 3]) -> Richelot<F> {
    let m = [coeffs(f, &g[0]), coeffs(f, &g[1]), coeffs(f, &g[2])];
    let delta = det3(f, &m);
    let h = [
        bracket(f, &g[1], &g[2]),
        bracket(f, &g[2], &g[0]),
        bracket(f, &g[0], &g[1]),
    ];
    Richelot { delta, h }
}

impl<F: Field> Richelot<F> {
    /// Codomain sextic Delta H1 H2 H3 (Delta != 0).
    pub fn codomain(&self, f: &F) -> Poly<F> {
        let p = poly::mul(f, &poly::mul(f, &self.h[0], &self.h[1]), &self.h[2]);
        poly::scale(f, &p, self.delta)
    }
}

impl<F: Field> Richelot<F> {
    /// Image of a point (x, y) of C (y != 0) under the Richelot correspondence (Smith 2005, §8.4):
    /// the two points (x', y') with G1(x) H1(x') + G2(x) H2(x') = 0 and
    /// y y' = G1(x) H1(x') (x - x'), on the codomain written as Delta y'^2 = H1 H2 H3.
    /// The x' are the roots of a quadratic over F, so the points are returned over the quadratic
    /// extension `ext` (with an embedding of F). Returns (x', y') pairs scaled to the model
    /// Y^2 = Delta H1 H2 H3 (Y = Delta y').
    pub fn image_point<G: Field>(
        &self,
        f: &F,
        g: &[Poly<F>; 3],
        ext: &G,
        emb: impl Fn(F::E) -> G::E,
        x: F::E,
        y: F::E,
        rng: &mut Rng,
    ) -> Vec<(G::E, G::E)> {
        let g1x = poly::eval(f, &g[0], x);
        let g2x = poly::eval(f, &g[1], x);
        let q = poly::add(
            f,
            &poly::scale(f, &self.h[0], g1x),
            &poly::scale(f, &self.h[1], g2x),
        );
        let qe: Vec<G::E> = q.iter().map(|&c| emb(c)).collect();
        let h1e: Vec<G::E> = self.h[0].iter().map(|&c| emb(c)).collect();
        // the x' are the roots of a quadratic: closed form (one square root), not Cantor-Zassenhaus
        let mut qe = qe;
        poly::trim(ext, &mut qe);
        let roots = match qe.len() {
            3 => {
                let (c, b, a) = (qe[0], qe[1], qe[2]);
                let disc = ext.sub(ext.mul(b, b), ext.mul(ext.from_u64(4), ext.mul(a, c)));
                match ext.sqrt(disc) {
                    Some(s) => {
                        let i2a = ext.inv(ext.add(a, a));
                        let r1 = ext.mul(ext.sub(s, b), i2a);
                        let r2 = ext.mul(ext.sub(ext.neg(s), b), i2a);
                        if r1 == r2 {
                            vec![r1]
                        } else {
                            vec![r1, r2]
                        }
                    }
                    None => vec![],
                }
            }
            2 => vec![ext.neg(ext.div(qe[0], qe[1]))],
            _ => poly::roots(ext, &qe, rng),
        };
        let (xe, ye, g1e, de) = (emb(x), emb(y), emb(g1x), emb(self.delta));
        roots
            .into_iter()
            .map(|xp| {
                let yp = ext.div(
                    ext.mul(ext.mul(g1e, poly::eval(ext, &h1e, xp)), ext.sub(xe, xp)),
                    ye,
                );
                (xp, ext.mul(de, yp))
            })
            .collect()
    }
}

/// Degenerate Richelot (Delta = 0): the two elliptic cubics (E1, E2) with J(C) ~ E1 x E2, when
/// the fixed points of the involution swapping the roots of each G_i are F-rational.
pub fn split<F: Field>(f: &F, fsex: &Poly<F>, g: &[Poly<F>; 3]) -> Option<(Poly<F>, Poly<F>)> {
    // sigma(x) = (alpha x + beta)/(gamma x - alpha) swaps the roots u, v of G_i iff
    // -alpha g_i1 + beta g_i2 - gamma g_i0 = 0
    let m: Vec<[F::E; 3]> = g
        .iter()
        .map(|gi| {
            let c = coeffs(f, gi);
            [f.neg(c[1]), c[2], f.neg(c[0])]
        })
        .collect();
    let (alpha, beta, gamma) = kernel3(f, &m)?;
    // coordinate z with sigma: z -> -z, x = (r - s z)/(1 - z) (s = infinity: x = r + z)
    let six = 6usize;
    let mut fz: Poly<F> = vec![f.zero(); 7];
    let fc: Vec<F::E> = (0..=six)
        .map(|k| if k < fsex.len() { fsex[k] } else { f.zero() })
        .collect();
    if f.is_zero(gamma) {
        // fixed points: r = beta/(2 alpha) and infinity; z = x - r... sigma(x) = -x - beta/alpha
        let r = f.neg(f.div(beta, f.add(alpha, alpha)));
        // f(r + z)
        let lin = vec![r, f.one()];
        let mut pw = vec![f.one()];
        for k in 0..=six {
            fz = poly::add(f, &fz, &poly::scale(f, &pw, fc[k]));
            pw = poly::mul(f, &pw, &lin);
        }
    } else {
        // gamma x^2 - 2 alpha x - beta = 0
        let disc = f.add(f.mul(alpha, alpha), f.mul(gamma, beta));
        let sq = f.sqrt(disc)?;
        let r = f.div(f.add(alpha, sq), gamma);
        let s = f.div(f.sub(alpha, sq), gamma);
        // sum_k a_k (r - s z)^k (1 - z)^(6-k)
        let num = vec![r, f.neg(s)];
        let den = vec![f.one(), f.neg(f.one())];
        for k in 0..=six {
            let mut t = vec![fc[k]];
            for _ in 0..k {
                t = poly::mul(f, &t, &num);
            }
            for _ in k..six {
                t = poly::mul(f, &t, &den);
            }
            fz = poly::add(f, &fz, &t);
        }
    }
    fz.resize(7, f.zero());
    if ![1, 3, 5].iter().all(|&k| f.is_zero(fz[k])) {
        return None;
    }
    let e1 = vec![fz[0], fz[2], fz[4], fz[6]];
    let e2 = vec![fz[6], fz[4], fz[2], fz[0]];
    Some((e1, e2))
}

/// A non-zero vector of the kernel of a 3x3 matrix (rows), if singular.
fn kernel3<F: Field>(f: &F, m: &[[F::E; 3]]) -> Option<(F::E, F::E, F::E)> {
    // cross products of row pairs
    for (a, b) in [(0, 1), (0, 2), (1, 2)] {
        let (u, v) = (m[a], m[b]);
        let c = (
            f.sub(f.mul(u[1], v[2]), f.mul(u[2], v[1])),
            f.sub(f.mul(u[2], v[0]), f.mul(u[0], v[2])),
            f.sub(f.mul(u[0], v[1]), f.mul(u[1], v[0])),
        );
        if !(f.is_zero(c.0) && f.is_zero(c.1) && f.is_zero(c.2)) {
            return Some(c);
        }
    }
    None
}

/// Short Weierstrass model y^2 = x^3 + A x + B of y^2 = a u^3 + b u^2 + c u + d.
pub fn cubic_to_short<F: Field>(f: &F, cub: &Poly<F>) -> Curve<F::E> {
    let (d, c, b, a) = (cub[0], cub[1], cub[2], cub[3]);
    // X = a u, Y = a y: Y^2 = X^3 + b X^2 + a c X + a^2 d; then X = t - b/3
    let b2 = b;
    let c2 = f.mul(a, c);
    let d2 = f.mul(f.mul(a, a), d);
    let three = f.from_u64(3);
    let s = f.div(b2, three);
    // t^3 + (c2 - b2^2/3) t + (2 b2^3/27 - b2 c2/3 + d2)
    let aa = f.sub(c2, f.mul(b2, s));
    let bb = f.add(
        f.sub(f.mul(f.from_u64(2), f.mul(s, f.mul(s, s))), f.mul(s, c2)),
        d2,
    );
    Curve::new(aa, bb)
}

/// Gluing: the genus-2 curve C with J(C) (2,2)-isogenous to E1 x E2 for E1: y^2 = prod(x - a_i),
/// E2: y^2 = prod(x - b_i), along the 2-torsion identification a_i <-> b_i (Howe–Leprévost–Poonen
/// 2000, in the normalisation used by Castryck–Decru 2022). None if the gluing degenerates.
pub fn glue<F: Field>(f: &F, a: [F::E; 3], b: [F::E; 3]) -> Option<Poly<F>> {
    let d = |x: F::E, y: F::E| f.sub(x, y);
    let (a1, a2, a3) = (a[0], a[1], a[2]);
    let (b1, b2, b3) = (b[0], b[1], b[2]);
    let sq = |x: F::E| f.mul(x, x);
    let inv_or = |x: F::E| if f.is_zero(x) { None } else { Some(f.inv(x)) };
    let aa1 = f.add(
        f.add(
            f.mul(sq(d(a1, a2)), inv_or(d(b1, b2))?),
            f.mul(sq(d(a2, a3)), inv_or(d(b2, b3))?),
        ),
        f.mul(sq(d(a3, a1)), inv_or(d(b3, b1))?),
    );
    let bb1 = f.add(
        f.add(
            f.mul(sq(d(b1, b2)), inv_or(d(a1, a2))?),
            f.mul(sq(d(b2, b3)), inv_or(d(a2, a3))?),
        ),
        f.mul(sq(d(b3, b1)), inv_or(d(a3, a1))?),
    );
    let aa2 = f.add(
        f.add(f.mul(a1, d(b3, b2)), f.mul(a2, d(b1, b3))),
        f.mul(a3, d(b2, b1)),
    );
    let bb2 = f.add(
        f.add(f.mul(b1, d(a3, a2)), f.mul(b2, d(a1, a3))),
        f.mul(b3, d(a2, a1)),
    );
    if f.is_zero(aa1) || f.is_zero(bb1) || f.is_zero(aa2) || f.is_zero(bb2) {
        return None;
    }
    let da = f.mul(f.mul(sq(d(a1, a2)), sq(d(a2, a3))), sq(d(a3, a1)));
    let db = f.mul(f.mul(sq(d(b1, b2)), sq(d(b2, b3))), sq(d(b3, b1)));
    let big_a = f.mul(db, f.div(aa1, aa2));
    let big_b = f.mul(da, f.div(bb1, bb2));
    let quad = |x1: F::E, x2: F::E, x3: F::E, y1: F::E, y2: F::E, y3: F::E| -> Poly<F> {
        // A (x2 - x1)(x1 - x3) x^2 + B (y2 - y1)(y1 - y3)
        vec![
            f.mul(big_b, f.mul(d(y2, y1), d(y1, y3))),
            f.zero(),
            f.mul(big_a, f.mul(d(x2, x1), d(x1, x3))),
        ]
    };
    let q1 = quad(a1, a2, a3, b1, b2, b3);
    let q2 = quad(a2, a3, a1, b2, b3, b1);
    let q3 = quad(a3, a1, a2, b3, b1, b2);
    // (the published formula carries a factor -1, a quadratic twist relative to these models of
    // E1 and E2; without it the L-polynomial matches L_{E1} L_{E2}, tested for p = 1 and 3 mod 4)
    let c = poly::mul(f, &poly::mul(f, &q1, &q2), &q3);
    if c.len() < 6 {
        return None;
    }
    Some(c)
}

/// Igusa–Clebsch invariants (I2, I4, I6, I10) of y^2 = lc prod (x - r_i) from its six roots.
pub fn igusa_clebsch<F: Field>(f: &F, lc: F::E, r: &[F::E; 6]) -> [F::E; 4] {
    let d = |i: usize, j: usize| f.sq(f.sub(r[i], r[j]));
    let mut i2 = f.zero();
    let mut i4 = f.zero();
    let mut i6 = f.zero();
    let mut i10 = f.one();
    for i in 0..6 {
        for j in i + 1..6 {
            i10 = f.mul(i10, d(i, j));
        }
    }
    // perfect matchings {a,b},{c,e},{g,h}
    for b in 1..6 {
        let rest: Vec<usize> = (1..6).filter(|&x| x != b).collect();
        let c = rest[0];
        for &e in &rest[1..] {
            let last: Vec<usize> = rest.iter().copied().filter(|&x| x != c && x != e).collect();
            i2 = f.add(i2, f.mul(f.mul(d(0, b), d(c, e)), d(last[0], last[1])));
        }
    }
    // splits into two triples {0,a,b} | {c,e,g}
    for a in 1..6 {
        for b in a + 1..6 {
            let t: Vec<usize> = (1..6).filter(|&x| x != a && x != b).collect();
            let s1 = [0, a, b];
            let tri = |s: [usize; 3]| f.mul(f.mul(d(s[0], s[1]), d(s[1], s[2])), d(s[2], s[0]));
            let s2 = [t[0], t[1], t[2]];
            let base = f.mul(tri(s1), tri(s2));
            i4 = f.add(i4, base);
            // the six bijections s1 -> s2
            for perm in [
                [0, 1, 2],
                [0, 2, 1],
                [1, 0, 2],
                [1, 2, 0],
                [2, 0, 1],
                [2, 1, 0],
            ] {
                let m = f.mul(
                    f.mul(d(s1[0], s2[perm[0]]), d(s1[1], s2[perm[1]])),
                    d(s1[2], s2[perm[2]]),
                );
                i6 = f.add(i6, f.mul(base, m));
            }
        }
    }
    let l2 = f.sq(lc);
    let l4 = f.sq(l2);
    let l6 = f.mul(l4, l2);
    let l10 = f.mul(l6, l4);
    [f.mul(l2, i2), f.mul(l4, i4), f.mul(l6, i6), f.mul(l10, i10)]
}

/// A key for the F-bar isomorphism class from (I2 : I4 : I6 : I10) in weighted projective space.
pub fn ic_key<F: Field>(f: &F, ic: &[F::E; 4]) -> (u8, Vec<F::E>) {
    let [i2, i4, i6, i10] = *ic;
    if !f.is_zero(i2) {
        let s = f.inv(i2);
        let s2 = f.sq(s);
        let s3 = f.mul(s2, s);
        let s5 = f.mul(s3, s2);
        (0, vec![f.mul(i4, s2), f.mul(i6, s3), f.mul(i10, s5)])
    } else if !f.is_zero(i4) {
        let s = f.inv(i4);
        let s3 = f.mul(f.sq(s), s);
        let s5 = f.mul(s3, f.sq(s));
        (1, vec![f.mul(f.sq(i6), s3), f.mul(f.sq(i10), s5)])
    } else if !f.is_zero(i6) {
        let s = f.inv(i6);
        let s5 = f.mul(f.sq(f.sq(s)), s);
        (2, vec![f.mul(f.mul(f.sq(i10), i10), s5)])
    } else {
        (3, vec![])
    }
}

/// A vertex of the superspecial (2,2)-graph: a Jacobian (sextic with its six roots) or a product
/// of two supersingular elliptic curves (with the roots of their 2-torsion cubics).
#[derive(Clone, Debug)]
pub enum Vertex<E> {
    Jac { lc: E, roots: [E; 6] },
    Prod { a: [E; 3], b: [E; 3] },
}

pub type VKey<E> = (u8, Vec<E>);

/// Isomorphism-class key of a vertex.
pub fn vertex_key(f: &Zp2, v: &Vertex<(u64, u64)>) -> VKey<(u64, u64)> {
    match v {
        Vertex::Jac { lc, roots } => ic_key(f, &igusa_clebsch(f, *lc, roots)),
        Vertex::Prod { a, b } => {
            let j = |r: &[(u64, u64); 3]| {
                let c = poly::from_roots(f, r);
                crate::curve::jinv(f, &cubic_to_short(f, &c))
            };
            let (j1, j2) = (j(a), j(b));
            (9, if j1 <= j2 { vec![j1, j2] } else { vec![j2, j1] })
        }
    }
}

/// The 15 (2,2)-neighbours of a vertex over F_{p^2} (all 2-torsion rational).
pub fn neighbours(f: &Zp2, v: &Vertex<(u64, u64)>, rng: &mut Rng) -> Vec<Vertex<(u64, u64)>> {
    let mut out = vec![];
    match v {
        Vertex::Jac { lc, roots } => {
            let fsex = poly::scale(f, &poly::from_roots(f, roots), *lc);
            for (pa, pb, pc) in matchings6() {
                let g = [
                    poly::scale(f, &poly::from_roots(f, &[roots[pa.0], roots[pa.1]]), *lc),
                    poly::from_roots(f, &[roots[pb.0], roots[pb.1]]),
                    poly::from_roots(f, &[roots[pc.0], roots[pc.1]]),
                ];
                let rl = richelot(f, &g);
                if f.is_zero(rl.delta) {
                    if let Some((e1, e2)) = split(f, &fsex, &g) {
                        let (ra, rb) = (poly::roots(f, &e1, rng), poly::roots(f, &e2, rng));
                        if ra.len() == 3 && rb.len() == 3 {
                            out.push(Vertex::Prod {
                                a: [ra[0], ra[1], ra[2]],
                                b: [rb[0], rb[1], rb[2]],
                            });
                        }
                    }
                } else {
                    let c = rl.codomain(f);
                    let rs = poly::roots(f, &c, rng);
                    if rs.len() == 6 {
                        out.push(Vertex::Jac {
                            lc: *c.last().unwrap(),
                            roots: [rs[0], rs[1], rs[2], rs[3], rs[4], rs[5]],
                        });
                    }
                }
            }
        }
        Vertex::Prod { a, b } => {
            // 9 products of 2-isogenies (kernel (a_i, 0) on E1, (b_j, 0) on E2)
            let two_iso =
                |r: &[(u64, u64); 3], k: usize, rng: &mut Rng| -> Option<[(u64, u64); 3]> {
                    let c = poly::from_roots(f, r);
                    let e = cubic_to_short(f, &c);
                    // short model shifts x by -b/3: the kernel abscissa is r_k + (sum r)/3
                    let shift = f.div(f.add(f.add(r[0], r[1]), r[2]), f.from_u64(3));
                    let x0 = f.sub(r[k], shift);
                    let iso = crate::kernel::velu::velu_cyclic(
                        f,
                        &e,
                        &crate::curve::Pt::Aff(x0, f.zero()),
                        2,
                    );
                    let rr = poly::roots(f, &vec![iso.cod.b, iso.cod.a, f.zero(), f.one()], rng);
                    if rr.len() == 3 {
                        Some([rr[0], rr[1], rr[2]])
                    } else {
                        None
                    }
                };
            for i in 0..3 {
                for j in 0..3 {
                    if let (Some(na), Some(nb)) = (two_iso(a, i, rng), two_iso(b, j, rng)) {
                        out.push(Vertex::Prod { a: na, b: nb });
                    }
                }
            }
            // 6 gluings
            for perm in [
                [0, 1, 2],
                [0, 2, 1],
                [1, 0, 2],
                [1, 2, 0],
                [2, 0, 1],
                [2, 1, 0],
            ] {
                let bp = [b[perm[0]], b[perm[1]], b[perm[2]]];
                match glue(f, *a, bp) {
                    Some(c) => {
                        let rs = poly::roots(f, &c, rng);
                        if rs.len() == 6 {
                            out.push(Vertex::Jac {
                                lc: *c.last().unwrap(),
                                roots: [rs[0], rs[1], rs[2], rs[3], rs[4], rs[5]],
                            });
                        }
                    }
                    None => {
                        // degenerate gluing (E1 = E2 along an isomorphism): the neighbour is a
                        // product again; omitted (counted by the caller as missing degree)
                    }
                }
            }
        }
    }
    out
}

/// The 15 perfect matchings of {0..5}.
pub fn matchings6() -> Vec<((usize, usize), (usize, usize), (usize, usize))> {
    let mut out = vec![];
    for b in 1..6 {
        let rest: Vec<usize> = (1..6).filter(|&x| x != b).collect();
        let c = rest[0];
        for &e in &rest[1..] {
            let last: Vec<usize> = rest.iter().copied().filter(|&x| x != c && x != e).collect();
            out.push(((0, b), (c, e), (last[0], last[1])));
        }
    }
    out
}

pub struct SsGraph {
    pub jacobians: usize,
    pub products: usize,
    /// sum over vertices of the number of neighbours found
    pub edges: usize,
}

/// Breadth-first search of the superspecial Richelot graph over F_{p^2} from E0 x E0
/// (E0: y^2 = x^3 + x, p = 3 mod 4), identifying vertices by Igusa–Clebsch invariants (Jacobians)
/// or by the pair of j-invariants (products).
pub fn superspecial_graph(p: u64, max_vertices: usize, rng: &mut Rng) -> SsGraph {
    let f = Zp2::new(p);
    let i = (0u64, 1u64); // sqrt(-1)
    let e0 = [f.zero(), i, f.neg(i)]; // x^3 + x = x (x - i)(x + i)
    let start = Vertex::Prod { a: e0, b: e0 };
    let mut seen: HashSet<VKey<(u64, u64)>> = HashSet::new();
    let mut queue = vec![start.clone()];
    seen.insert(vertex_key(&f, &start));
    let (mut jac, mut prod, mut edges) = (0, 1, 0);
    let mut k = 0;
    let mut kinds: HashMap<u8, usize> = HashMap::new();
    while k < queue.len() && seen.len() < max_vertices {
        let v = queue[k].clone();
        k += 1;
        let ns = neighbours(&f, &v, rng);
        edges += ns.len();
        for n in ns {
            let key = vertex_key(&f, &n);
            if seen.insert(key.clone()) {
                *kinds.entry(key.0).or_default() += 1;
                match n {
                    Vertex::Jac { .. } => jac += 1,
                    Vertex::Prod { .. } => prod += 1,
                }
                queue.push(n);
            }
        }
    }
    SsGraph {
        jacobians: jac,
        products: prod,
        edges,
    }
}
