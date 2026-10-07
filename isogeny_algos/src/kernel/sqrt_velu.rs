//! sqrt-Velu (Bernstein-De Feo-Leroux-Smith 2020), adapted to short Weierstrass form.
//! The kernel indices {1..n}, n=(l-1)/2, are split as I (block centres) u (I +- J) u K with
//! J={1..b}; products over I+-J are evaluated as a resultant-style product
//! prod_{c in I} G(x_c), G(Z) = prod_{j in J} F_j(X,Z), in a truncated power-series ring.
//! Baseline: naive polynomial arithmetic (no multipoint evaluation), so the asymptotic
//! saving is NOT yet realised; this is the structural baseline for later optimisation.
use crate::curve::*;
use crate::field::Field;

type Ser<F> = Vec<<F as Field>::E>; // truncated power series, length m

fn smul<F: Field>(f: &F, a: &Ser<F>, b: &Ser<F>) -> Ser<F> {
    let m = a.len();
    let mut r = vec![f.zero(); m];
    for i in 0..m {
        if f.is_zero(a[i]) {
            continue;
        }
        for j in 0..m - i {
            r[i + j] = f.add(r[i + j], f.mul(a[i], b[j]));
        }
    }
    r
}
fn sadd<F: Field>(f: &F, a: &Ser<F>, b: &Ser<F>) -> Ser<F> {
    a.iter().zip(b).map(|(&x, &y)| f.add(x, y)).collect()
}
fn sconst<F: Field>(f: &F, c: F::E, m: usize) -> Ser<F> {
    let mut v = vec![f.zero(); m];
    v[0] = c;
    v
}
fn sinv<F: Field>(f: &F, a: &Ser<F>) -> Ser<F> {
    let m = a.len();
    let a0i = f.inv(a[0]);
    let mut r = vec![f.zero(); m];
    r[0] = a0i;
    for k in 1..m {
        let mut s = f.zero();
        for i in 1..=k {
            s = f.add(s, f.mul(a[i], r[k - i]));
        }
        r[k] = f.neg(f.mul(s, a0i));
    }
    r
}

pub struct SqrtVeluIso<F: Field> {
    pub dom: Curve<F::E>,
    pub cod: Curve<F::E>,
    pub deg: u64,
    n: usize,
    p1: F::E,
    xs_j: Vec<F::E>,
    xs_i: Vec<F::E>,
    xs_k: Vec<F::E>,
    d0: F::E,
}

fn pt_x<F: Field>(p: &Pt<F::E>) -> F::E {
    match p {
        Pt::Aff(x, _) => *x,
        Pt::Inf => panic!("kernel point of small order"),
    }
}

/// Choose block parameters and compute x-coordinates of I, J, K points from a generator P.
fn index_sets<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    p: &Pt<F::E>,
    ell: u64,
) -> (usize, Vec<F::E>, Vec<F::E>, Vec<F::E>) {
    let n = ((ell - 1) / 2) as usize;
    let mut b = ((n as f64).sqrt() / 2.0).floor() as usize;
    if b < 1 {
        b = 1;
    }
    let bp = n / (2 * b + 1);
    if bp == 0 {
        let mut ks = vec![];
        let mut cur = *p;
        for _ in 0..n {
            ks.push(pt_x::<F>(&cur));
            cur = padd(f, e, &cur, p);
        }
        return (0, vec![], vec![], ks);
    }
    let klen = n - (2 * b + 1) * bp;
    let msmall = (2 * b + 1).max(b + klen);
    let mut small = vec![Pt::Inf; msmall + 1];
    small[1] = *p;
    for k in 2..=msmall {
        small[k] = padd(f, e, &small[k - 1], p);
    }
    let step = small[2 * b + 1];
    let mut xs_i = vec![];
    let mut c = small[b + 1];
    let mut centres = vec![];
    for _ in 0..bp {
        xs_i.push(pt_x::<F>(&c));
        centres.push(c);
        c = padd(f, e, &c, &step);
    }
    let xs_j: Vec<F::E> = (1..=b).map(|j| pt_x::<F>(&small[j])).collect();
    let last = centres[bp - 1];
    let mut xs_k = vec![];
    for r in 1..=klen {
        // index (2b+1) bp + r = c_{bp-1} + (b + r)
        let q = padd(f, e, &last, &small[b + r]);
        xs_k.push(pt_x::<F>(&q));
    }
    (b, xs_j, xs_i, xs_k)
}

/// Polynomial in Z with series coefficients: coeffs[k] is coefficient of Z^k.
type PolyS<F> = Vec<Ser<F>>;

fn pmul_s<F: Field>(f: &F, a: &PolyS<F>, b: &PolyS<F>) -> PolyS<F> {
    let m = a[0].len();
    let mut r: PolyS<F> = vec![vec![f.zero(); m]; a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            let t = smul(f, &a[i], &b[j]);
            r[i + j] = sadd(f, &r[i + j], &t);
        }
    }
    r
}

/// prod_{k in I+-J} (X - x_k) * D0 as a series, where X is given as the series `xs` ("alpha + delta"),
/// or, when `inf` is true, the series in eps of prod (1 - eps x_k) * D0.
fn ij_product<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    xs_j: &[F::E],
    xs_i: &[F::E],
    m: usize,
    inf: bool,
    alpha: F::E,
) -> Ser<F> {
    let c = |v: u64| f.from_u64(v);
    // build G(Z) = prod_j F_j
    let mut g: PolyS<F> = vec![sconst(f, f.one(), m)];
    for &xj in xs_j {
        // F_j(X,Z) = (Z-xj)^2 X^2 - (2(Z+xj)(Z xj + A)+4B) X + ((Z xj - A)^2 - 4B(Z+xj))
        // coefficients in Z: c2 Z^2 + c1 Z + c0 where each is a polynomial in X
        let xj2 = f.mul(xj, xj);
        // X^2 coeff: Z^2 - 2xj Z + xj^2
        let a2 = [xj2, f.neg(f.mul(c(2), xj)), f.one()];
        // X^1 coeff (negated): 2 xj Z^2 + 2(A + xj^2) Z + 2 A xj + 4B
        let a1 = [
            f.add(f.mul(c(2), f.mul(e.a, xj)), f.mul(c(4), e.b)),
            f.mul(c(2), f.add(e.a, xj2)),
            f.mul(c(2), xj),
        ];
        // X^0 coeff: xj^2 Z^2 - (2 A xj + 4B) Z + A^2 - 4B xj
        let a0 = [
            f.sub(f.mul(e.a, e.a), f.mul(c(4), f.mul(e.b, xj))),
            f.neg(f.add(f.mul(c(2), f.mul(e.a, xj)), f.mul(c(4), e.b))),
            xj2,
        ];
        let mut fj: PolyS<F> = vec![];
        for k in 0..3 {
            let mut s = vec![f.zero(); m];
            if inf {
                // eps^2 F(1/eps) = a2 - eps*a1 + eps^2*a0 (a1 stored as negated X coeff)
                s[0] = a2[k];
                if m > 1 {
                    s[1] = f.neg(a1[k]);
                }
                if m > 2 {
                    s[2] = a0[k];
                }
            } else {
                // X = alpha + delta: a2 (alpha+d)^2 - a1 (alpha+d) + a0
                let al = alpha;
                s[0] = f.add(f.sub(f.mul(a2[k], f.mul(al, al)), f.mul(a1[k], al)), a0[k]);
                if m > 1 {
                    s[1] = f.sub(f.mul(f.mul(c(2), al), a2[k]), a1[k]);
                }
                if m > 2 {
                    s[2] = a2[k];
                }
            }
            fj.push(s);
        }
        g = pmul_s(f, &g, &fj);
    }
    // evaluate G at each x_c by Horner, multiply
    let mut acc = sconst(f, f.one(), m);
    for &xc in xs_i {
        let mut r = vec![f.zero(); m];
        for k in (0..g.len()).rev() {
            r = smul(f, &r, &sconst(f, xc, m));
            r = sadd(f, &r, &g[k]);
        }
        acc = smul(f, &acc, &r);
    }
    acc
}

fn d0_of<F: Field>(f: &F, xs_j: &[F::E], xs_i: &[F::E]) -> F::E {
    let mut d0 = f.one();
    for &xc in xs_i {
        for &xj in xs_j {
            let d = f.sub(xc, xj);
            d0 = f.mul(d0, f.mul(d, d));
        }
    }
    d0
}

/// sqrt-Velu codomain from an x-only generator... here the generator point P.
pub fn sqrt_velu<F: Field>(f: &F, e: &Curve<F::E>, p: &Pt<F::E>, ell: u64) -> SqrtVeluIso<F> {
    let n = ((ell - 1) / 2) as usize;
    let (_b, xs_j, xs_i, xs_k) = index_sets::<F>(f, e, p, ell);
    let m = 4;
    let d0 = d0_of(f, &xs_j, &xs_i);
    // series of prod (1 - eps x_k)
    let mut prod = sconst(f, f.one(), m);
    if !xs_i.is_empty() {
        let mut ij = ij_product(f, e, &xs_j, &xs_i, m, true, f.zero());
        let inv0 = f.inv(ij[0]);
        ij = ij.iter().map(|&c| f.mul(c, inv0)).collect();
        prod = smul(f, &prod, &ij);
    }
    for &x in xs_i.iter().chain(xs_k.iter()) {
        let lin = vec![f.one(), f.neg(x), f.zero(), f.zero()];
        prod = smul(f, &prod, &lin);
    }
    let e1 = f.neg(prod[1]);
    let e2 = prod[2];
    let e3 = f.neg(prod[3]);
    let p1 = e1;
    let p2 = f.sub(f.mul(e1, e1), f.mul(f.from_u64(2), e2));
    let p3 = f.add(
        f.sub(
            f.mul(e1, f.mul(e1, e1)),
            f.mul(f.from_u64(3), f.mul(e1, e2)),
        ),
        f.mul(f.from_u64(3), e3),
    );
    let cod = super::kohel::codomain_from_sums(f, e, n, (p1, p2, p3));
    SqrtVeluIso {
        dom: *e,
        cod,
        deg: ell,
        n,
        p1,
        xs_j,
        xs_i,
        xs_k,
        d0,
    }
}

// ------------------------------------------------------------------ fast version
// R = F[delta]/delta^m;  R[Z] polynomials are packed into F[Y] with blocks of width 2m-1
// (Kronecker substitution), so products use the fast univariate multiplication.

fn pack<F: Field>(f: &F, a: &PolyS<F>, m: usize) -> Vec<F::E> {
    let w = 2 * m - 1;
    let mut v = vec![f.zero(); a.len() * w];
    for (k, c) in a.iter().enumerate() {
        v[k * w..k * w + m].copy_from_slice(&c[..m]);
    }
    v
}
fn unpack<F: Field>(f: &F, v: &[F::E], m: usize, len: usize) -> PolyS<F> {
    let w = 2 * m - 1;
    (0..len)
        .map(|k| {
            (0..m)
                .map(|t| {
                    if k * w + t < v.len() {
                        v[k * w + t]
                    } else {
                        f.zero()
                    }
                })
                .collect()
        })
        .collect()
}
fn mul_rz<F: Field>(f: &F, a: &PolyS<F>, b: &PolyS<F>, m: usize) -> PolyS<F> {
    let prod = crate::poly::mul_raw(f, &pack(f, a, m), &pack(f, b, m));
    unpack(f, &prod, m, a.len() + b.len() - 1)
}

/// Coefficients (in Z, each a series in delta) of F_j(X, Z) for X = 1/eps (inf) or X = alpha + delta.
fn fj_coeffs<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    xj: F::E,
    m: usize,
    inf: bool,
    alpha: F::E,
) -> PolyS<F> {
    let c = |v: u64| f.from_u64(v);
    let xj2 = f.mul(xj, xj);
    let a2 = [xj2, f.neg(f.mul(c(2), xj)), f.one()];
    let a1 = [
        f.add(f.mul(c(2), f.mul(e.a, xj)), f.mul(c(4), e.b)),
        f.mul(c(2), f.add(e.a, xj2)),
        f.mul(c(2), xj),
    ];
    let a0 = [
        f.sub(f.mul(e.a, e.a), f.mul(c(4), f.mul(e.b, xj))),
        f.neg(f.add(f.mul(c(2), f.mul(e.a, xj)), f.mul(c(4), e.b))),
        xj2,
    ];
    (0..3)
        .map(|k| {
            let mut s = vec![f.zero(); m];
            if inf {
                s[0] = a2[k];
                if m > 1 {
                    s[1] = f.neg(a1[k]);
                }
                if m > 2 {
                    s[2] = a0[k];
                }
            } else {
                let al = alpha;
                s[0] = f.add(f.sub(f.mul(a2[k], f.mul(al, al)), f.mul(a1[k], al)), a0[k]);
                if m > 1 {
                    s[1] = f.sub(f.mul(f.mul(c(2), al), a2[k]), a1[k]);
                }
                if m > 2 {
                    s[2] = a2[k];
                }
            }
            s
        })
        .collect()
}

/// Same output as `ij_product`, in O~(sqrt l): product tree over J in R[Z], multipoint evaluation at I.
fn ij_product_fast<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    xs_j: &[F::E],
    xs_i: &[F::E],
    m: usize,
    inf: bool,
    alpha: F::E,
    tree_i: &[Vec<crate::poly::Poly<F>>],
) -> Ser<F> {
    let mut layer: Vec<PolyS<F>> = xs_j
        .iter()
        .map(|&xj| fj_coeffs(f, e, xj, m, inf, alpha))
        .collect();
    while layer.len() > 1 {
        layer = layer
            .chunks(2)
            .map(|c| {
                if c.len() == 2 {
                    mul_rz(f, &c[0], &c[1], m)
                } else {
                    c[0].clone()
                }
            })
            .collect();
    }
    let g = &layer[0];
    // evaluate each delta-component at the centres
    let mut vals: Vec<Ser<F>> = vec![vec![f.zero(); m]; xs_i.len()];
    for t in 0..m {
        let mut comp: crate::poly::Poly<F> = g.iter().map(|c| c[t]).collect();
        crate::poly::trim(f, &mut comp);
        let ev = if xs_i.len() <= 8 {
            xs_i.iter()
                .map(|&x| crate::poly::eval(f, &comp, x))
                .collect()
        } else {
            crate::poly::multipoint_eval_tree(f, &comp, tree_i)
        };
        for (i, v) in ev.into_iter().enumerate() {
            vals[i][t] = v;
        }
    }
    // product of the values in R (balanced)
    while vals.len() > 1 {
        vals = vals
            .chunks(2)
            .map(|c| {
                if c.len() == 2 {
                    smul(f, &c[0], &c[1])
                } else {
                    c[0].clone()
                }
            })
            .collect();
    }
    vals.pop().unwrap_or_else(|| sconst(f, f.one(), m))
}

/// x-only Weierstrass ladder: ((X_m : Z_m), (X_{m+1} : Z_{m+1})) for x(P) = x0.
fn ladder_pair<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    x0: F::E,
    m: usize,
) -> ((F::E, F::E), (F::E, F::E)) {
    use crate::kernel::xonly::{xadd_c, xdbl_c, XConst};
    let c = XConst::new(f, e);
    let xadd_w = |_f: &F, _e: &Curve<F::E>, p: (F::E, F::E), q: (F::E, F::E), d: (F::E, F::E)| {
        xadd_c(f, &c, p, q, d)
    };
    let xdbl_w = |_f: &F, _e: &Curve<F::E>, p: (F::E, F::E)| xdbl_c(f, &c, p);
    let p = (x0, f.one());
    let (mut r0, mut r1) = (p, xdbl_w(f, e, p));
    let bits = usize::BITS - m.leading_zeros();
    for i in (0..bits.saturating_sub(1)).rev() {
        if (m >> i) & 1 == 0 {
            r1 = xadd_w(f, e, r0, r1, p);
            r0 = xdbl_w(f, e, r0);
        } else {
            r0 = xadd_w(f, e, r0, r1, p);
            r1 = xdbl_w(f, e, r1);
        }
    }
    (r0, r1)
}

/// sqrt-Velu from x(P) only: index sets from projective x-only chains + one batch inversion,
/// products in R[Z] by product trees and evaluation at the centres by remainder trees.
pub fn sqrt_velu_fast<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    x0: F::E,
    ell: u64,
) -> Option<SqrtVeluIso<F>> {
    use crate::kernel::xonly::{xadd_c, xdbl_c, XConst};
    let xc = XConst::new(f, e);
    let xadd_w = |_f: &F, _e: &Curve<F::E>, p: (F::E, F::E), q: (F::E, F::E), d: (F::E, F::E)| {
        xadd_c(f, &xc, p, q, d)
    };
    let n = ((ell - 1) / 2) as usize;
    let mut b = ((n as f64 / 2.0).sqrt()).floor() as usize;
    if b < 1 {
        b = 1;
    }
    let bp = n / (2 * b + 1);
    if bp < 2 {
        // tiny degree: plain x-only Velu data in the same structure
        let xs = crate::kernel::xonly::multiples_x_fast(f, e, x0, n)?;
        let mut prod = sconst(f, f.one(), 4);
        for &x in &xs {
            prod = smul(f, &prod, &vec![f.one(), f.neg(x), f.zero(), f.zero()]);
        }
        let (e1, e2, e3) = (f.neg(prod[1]), prod[2], f.neg(prod[3]));
        let p1 = e1;
        let p2 = f.sub(f.mul(e1, e1), f.mul(f.from_u64(2), e2));
        let p3 = f.add(
            f.sub(
                f.mul(e1, f.mul(e1, e1)),
                f.mul(f.from_u64(3), f.mul(e1, e2)),
            ),
            f.mul(f.from_u64(3), e3),
        );
        let cod = super::kohel::codomain_from_sums(f, e, n, (p1, p2, p3));
        return Some(SqrtVeluIso {
            dom: *e,
            cod,
            deg: ell,
            n,
            p1,
            xs_j: vec![],
            xs_i: vec![],
            xs_k: xs,
            d0: f.one(),
        });
    }
    let p = (x0, f.one());
    // J chain: x_1 .. x_{2b+1}
    let mut proj: Vec<(F::E, F::E)> = vec![p, xdbl_c(f, &xc, p)];
    for k in 2..=2 * b {
        let d = proj[k - 2];
        if f.is_zero(d.0) {
            return None;
        }
        proj.push(xadd_w(f, e, proj[k - 1], p, d));
    }
    let s = proj[2 * b]; // (2b+1)P
                         // I chain: c_0 = (b+1)P, c_{-1} = -bP, c_{i+1} = c_i + S
    let mut centres: Vec<(F::E, F::E)> = vec![proj[b]];
    let mut prev = proj[b - 1];
    for _ in 1..bp {
        if f.is_zero(prev.0) {
            return None;
        }
        let next = xadd_w(f, e, *centres.last().unwrap(), s, prev);
        prev = *centres.last().unwrap();
        centres.push(next);
    }
    // K: indices t0+1 .. n with t0 = (2b+1) bp, from a ladder pair at t0 and a chain
    let t0 = (2 * b + 1) * bp;
    let klen = n - t0;
    let mut kpts: Vec<(F::E, F::E)> = vec![];
    if klen > 0 {
        let (rt0, rt1) = ladder_pair(f, e, x0, t0);
        let (mut a0, mut a1) = (rt0, rt1);
        kpts.push(a1);
        for _ in 1..klen {
            if f.is_zero(a0.0) {
                return None;
            }
            let a2 = xadd_w(f, e, a1, p, a0);
            a0 = a1;
            a1 = a2;
            kpts.push(a1);
        }
    }
    let all: Vec<(F::E, F::E)> = proj[..b]
        .iter()
        .chain(centres.iter())
        .chain(kpts.iter())
        .copied()
        .collect();
    let zs: Vec<F::E> = all.iter().map(|q| q.1).collect();
    if zs.iter().any(|&z| f.is_zero(z)) {
        return None;
    }
    let zi = f.batch_inv(&zs);
    let xs: Vec<F::E> = all.iter().zip(zi).map(|(q, iz)| f.mul(q.0, iz)).collect();
    let xs_j = xs[..b].to_vec();
    let xs_i = xs[b..b + bp].to_vec();
    let xs_k = xs[b + bp..].to_vec();
    let m = 4;
    let tree_i = crate::poly::subproduct_tree(f, &xs_i);
    // D0 = prod_c h_J(x_c)^2
    let hj = crate::poly::from_roots(f, &xs_j);
    let hv = if xs_i.len() <= 8 {
        xs_i.iter()
            .map(|&x| crate::poly::eval(f, &hj, x))
            .collect::<Vec<_>>()
    } else {
        crate::poly::multipoint_eval_tree(f, &hj, &tree_i)
    };
    let mut d0 = f.one();
    for v in hv {
        d0 = f.mul(d0, f.sq(v));
    }
    let mut prod = sconst(f, f.one(), m);
    let ij = ij_product_fast(f, e, &xs_j, &xs_i, m, true, f.zero(), &tree_i);
    let inv0 = f.inv(ij[0]);
    prod = smul(
        f,
        &prod,
        &ij.iter().map(|&c| f.mul(c, inv0)).collect::<Vec<_>>(),
    );
    for &x in xs_i.iter().chain(xs_k.iter()) {
        prod = smul(f, &prod, &vec![f.one(), f.neg(x), f.zero(), f.zero()]);
    }
    let (e1, e2, e3) = (f.neg(prod[1]), prod[2], f.neg(prod[3]));
    let p1 = e1;
    let p2 = f.sub(f.mul(e1, e1), f.mul(f.from_u64(2), e2));
    let p3 = f.add(
        f.sub(
            f.mul(e1, f.mul(e1, e1)),
            f.mul(f.from_u64(3), f.mul(e1, e2)),
        ),
        f.mul(f.from_u64(3), e3),
    );
    let cod = super::kohel::codomain_from_sums(f, e, n, (p1, p2, p3));
    let _ = d0;
    Some(SqrtVeluIso {
        dom: *e,
        cod,
        deg: ell,
        n,
        p1,
        xs_j,
        xs_i,
        xs_k,
        d0,
    })
}

impl<F: Field> SqrtVeluIso<F> {
    /// (f(x), f'(x)) for x outside the kernel.
    fn eval_xd(&self, f: &F, alpha: F::E) -> Option<(F::E, F::E)> {
        let e = &self.dom;
        let m = 4;
        // h(alpha+delta) = prod_{k in S}(alpha + delta - x_k)
        let mut h = sconst(f, f.one(), m);
        if !self.xs_i.is_empty() {
            let ij = if self.xs_i.len() > 8 {
                let tree_i = crate::poly::subproduct_tree(f, &self.xs_i);
                ij_product_fast(f, e, &self.xs_j, &self.xs_i, m, false, alpha, &tree_i)
            } else {
                ij_product(f, e, &self.xs_j, &self.xs_i, m, false, alpha)
            };
            let d0i = f.inv(self.d0);
            let ij: Ser<F> = ij.iter().map(|&c| f.mul(c, d0i)).collect();
            h = smul(f, &h, &ij);
        }
        for &x in self.xs_i.iter().chain(self.xs_k.iter()) {
            let lin = vec![f.sub(alpha, x), f.one(), f.zero(), f.zero()];
            h = smul(f, &h, &lin);
        }
        if f.is_zero(h[0]) {
            return None;
        }
        // L = h'/h mod delta^3: L0 = s1, L1 = -s2, L2 = s3
        let hp = vec![
            h[1],
            f.mul(f.from_u64(2), h[2]),
            f.mul(f.from_u64(3), h[3]),
            f.zero(),
        ];
        let l = smul(f, &hp, &sinv(f, &h));
        let (s1, s2, s3) = (l[0], f.neg(l[1]), l[2]);
        let c = |v: u64| f.from_u64(v);
        let al2 = f.mul(alpha, alpha);
        let fa = f.add(f.mul(f.add(al2, e.a), alpha), e.b);
        let nn = c(self.n as u64);
        let g1 = f.add(f.mul(c(6), al2), f.mul(c(2), e.a)); // 6a^2+2A
                                                            // f = a - g1 s1 + 4F s2 + 2(n a - p1)
        let fx = f.add(
            f.add(alpha, f.mul(c(2), f.sub(f.mul(nn, alpha), self.p1))),
            f.add(f.neg(f.mul(g1, s1)), f.mul(c(4), f.mul(fa, s2))),
        );
        // f' = 1 - 12 a s1 - g1 s1' + 4F' s2 + 4F s2' + 2n ; s1' = -s2, s2' = -2 s3
        let fap = f.add(f.mul(c(3), al2), e.a);
        let fpx = f.add(
            f.add(
                f.add(f.one(), f.mul(c(2), nn)),
                f.neg(f.mul(f.mul(c(12), alpha), s1)),
            ),
            f.add(
                f.add(f.mul(g1, s2), f.mul(c(4), f.mul(fap, s2))),
                f.neg(f.mul(c(8), f.mul(fa, s3))),
            ),
        );
        Some((fx, fpx))
    }
}

impl<F: Field> Isogeny<F> for SqrtVeluIso<F> {
    fn domain(&self) -> &Curve<F::E> {
        &self.dom
    }
    fn codomain(&self) -> &Curve<F::E> {
        &self.cod
    }
    fn degree(&self) -> u64 {
        self.deg
    }
    fn eval_x(&self, f: &F, x: F::E) -> Option<F::E> {
        self.eval_xd(f, x).map(|v| v.0)
    }
    fn eval(&self, f: &F, p: &Pt<F::E>) -> Pt<F::E> {
        match *p {
            Pt::Inf => Pt::Inf,
            Pt::Aff(x, y) => match self.eval_xd(f, x) {
                None => Pt::Inf,
                Some((fx, fpx)) => Pt::Aff(fx, f.mul(y, fpx)),
            },
        }
    }
}
