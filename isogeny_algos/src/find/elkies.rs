//! Elkies' normalised codomain from Phi_l (via derivatives of Phi), then the isogeny itself by
//! Bostan-Morain-Salvy-Schost: the x-map f satisfies f(wp_E(z)) = wp_E'(z); expand both in z,
//! and recover f as a Pade approximant of degree (l-1, l). Baseline: quadratic series
//! arithmetic and dense Gaussian elimination for the Pade step.
use super::modpoly::Phi;
use crate::curve::*;
use crate::field::Field;
use crate::poly::{self, Poly};

/// Coefficients c_k (k>=2) of wp(z) = z^-2 (1 + sum_{k>=2} c_k z^{2k}) as U(v) = 1 + sum c_k v^k.
pub fn u_series<F: Field>(f: &F, e: &Curve<F::E>, n: usize) -> Vec<F::E> {
    let mut u = vec![f.zero(); n + 1];
    u[0] = f.one();
    if n >= 2 {
        u[2] = f.neg(f.div(e.a, f.from_u64(5)));
    }
    if n >= 3 {
        u[3] = f.neg(f.div(e.b, f.from_u64(7)));
    }
    for k in 4..=n {
        let mut s = f.zero();
        for m in 2..=k - 2 {
            s = f.add(s, f.mul(u[m], u[k - m]));
        }
        let den = f.from_u64(((2 * k + 1) * (k - 3)) as u64);
        u[k] = f.div(f.mul(f.from_u64(3), s), den);
    }
    u
}

fn ser_mul<F: Field>(f: &F, a: &[F::E], b: &[F::E], n: usize) -> Vec<F::E> {
    let mut r = vec![f.zero(); n];
    for i in 0..a.len().min(n) {
        for j in 0..b.len().min(n - i) {
            r[i + j] = f.add(r[i + j], f.mul(a[i], b[j]));
        }
    }
    r
}
/// a(b(w)) with b(0)=0, to n terms
fn ser_compose<F: Field>(f: &F, a: &[F::E], b: &[F::E], n: usize) -> Vec<F::E> {
    let mut res = vec![f.zero(); n];
    for k in (0..a.len()).rev() {
        res = ser_mul(f, &res, b, n);
        res[0] = f.add(res[0], a[k]);
    }
    res
}
fn ser_inv<F: Field>(f: &F, a: &[F::E], n: usize) -> Vec<F::E> {
    let a0i = f.inv(a[0]);
    let mut r = vec![f.zero(); n];
    r[0] = a0i;
    for k in 1..n {
        let mut s = f.zero();
        for i in 1..=k.min(a.len() - 1) {
            s = f.add(s, f.mul(a[i], r[k - i]));
        }
        r[k] = f.neg(f.mul(s, a0i));
    }
    r
}

pub fn solve_linear<F: Field>(f: &F, mut m: Vec<Vec<F::E>>, mut b: Vec<F::E>) -> Option<Vec<F::E>> {
    let n = b.len();
    for col in 0..n {
        let pr = (col..n).find(|&r| !f.is_zero(m[r][col]))?;
        m.swap(col, pr);
        b.swap(col, pr);
        let inv = f.inv(m[col][col]);
        for c in col..n {
            m[col][c] = f.mul(m[col][c], inv);
        }
        b[col] = f.mul(b[col], inv);
        for r in 0..n {
            if r != col && !f.is_zero(m[r][col]) {
                let fct = m[r][col];
                for c in col..n {
                    let t = f.mul(fct, m[col][c]);
                    m[r][c] = f.sub(m[r][c], t);
                }
                let t = f.mul(fct, b[col]);
                b[r] = f.sub(b[r], t);
            }
        }
    }
    Some(b)
}

/// Normalised codomain curve E' for the root `jt` of Phi_l(j, Y) (Elkies' formulas).
pub fn elkies_codomain<F: Field>(
    f: &F,
    phi: &Phi<F>,
    e: &Curve<F::E>,
    jt: F::E,
) -> Option<Curve<F::E>> {
    let c = |v: u64| f.from_u64(v);
    let j = jinv(f, e);
    if j == c(0) || j == c(1728) || jt == c(0) || jt == c(1728) {
        return None;
    }
    let e4 = f.mul(f.neg(c(48)), e.a);
    let e6 = f.mul(f.neg(c(864)), e.b);
    let (px, py) = phi.grad(f, j, jt);
    if f.is_zero(py) || f.is_zero(px) {
        return None;
    }
    let dj = f.neg(f.div(f.mul(e6, j), e4));
    let r = f.div(f.mul(f.mul(c(phi.ell as u64), px), dj), f.mul(py, jt));
    let e4p = f.div(f.mul(jt, f.mul(r, r)), f.sub(jt, c(1728)));
    let e6p = f.mul(r, e4p);
    Some(Curve::new(
        f.neg(f.div(e4p, c(48))),
        f.neg(f.div(e6p, c(864))),
    ))
}

/// Normalised l-isogeny E -> E' (E' from `elkies_codomain`), via BMSS series + Pade.
pub fn bmss_isogeny<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    ep: &Curve<F::E>,
    ell: usize,
) -> Option<RatIsogeny<F>> {
    let n = 2 * ell + 2;
    let u = u_series(f, e, n);
    let up = u_series(f, ep, n);
    // U'(v) for wp_{E'} = U'(v)/v ; named `up`. v(w) solves v = w U(v)
    let mut v = vec![f.zero(); n];
    v[1] = f.one();
    for _ in 0..n {
        let uv = ser_compose(f, &u, &v, n);
        let mut nv = vec![f.zero(); n];
        for i in 0..n - 1 {
            nv[i + 1] = uv[i];
        }
        v = nv;
    }
    // sigma(w) = U(v(w)) / U'(v(w))
    let uv = ser_compose(f, &u, &v, n);
    let upv = ser_compose(f, &up, &v, n);
    let sigma = ser_mul(f, &uv, &ser_inv(f, &upv, n), n);
    // Pade: sum_{i=1}^{l} nt_i sigma_{k-i} = -sigma_k, k = l..2l-1
    let mut m = vec![vec![f.zero(); ell]; ell];
    let mut rhs = vec![f.zero(); ell];
    for k in ell..2 * ell {
        for i in 1..=ell {
            m[k - ell][i - 1] = sigma[k - i];
        }
        rhs[k - ell] = f.neg(sigma[k]);
    }
    let sol = solve_linear(f, m, rhs)?;
    let mut nt = vec![f.one()];
    nt.extend(sol);
    let mut dt = vec![f.zero(); ell];
    for k in 0..ell {
        let mut s = f.zero();
        for i in 0..=k {
            s = f.add(s, f.mul(nt[i], sigma[k - i]));
        }
        dt[k] = s;
    }
    let num: Poly<F> = (0..=ell).map(|i| nt[ell - i]).collect();
    let den: Poly<F> = (0..ell).map(|i| dt[ell - 1 - i]).collect();
    if f.is_zero(*num.last()?) || f.is_zero(*den.last()?) {
        return None;
    }
    // h = sqrt(den)
    let d = (ell - 1) / 2;
    let mut h = vec![f.zero(); d + 1];
    h[d] = f.one();
    let inv2 = f.inv(f.from_u64(2));
    for k in 1..=d {
        let mut s = den[2 * d - k];
        for i in 1..k {
            s = f.sub(s, f.mul(h[d - i], h[d - k + i]));
        }
        h[d - k] = f.mul(s, inv2);
    }
    if poly::mul(f, &h, &h) != den {
        return None;
    }
    Some(RatIsogeny {
        dom: *e,
        cod: *ep,
        deg: ell as u64,
        ker: h,
        num,
        den,
    })
}

/// All F-rational normalised l-isogenies out of `e` found through Phi_l roots.
pub fn elkies_isogenies<F: Field>(
    f: &F,
    phi: &Phi<F>,
    e: &Curve<F::E>,
    rng: &mut crate::field::Rng,
) -> Vec<RatIsogeny<F>> {
    let j = jinv(f, e);
    let mut out = vec![];
    for jt in phi.neighbors(f, j, rng) {
        if let Some(ep) = elkies_codomain(f, phi, e, jt) {
            if let Some(iso) = bmss_isogeny(f, e, &ep, phi.ell) {
                out.push(iso);
            }
        }
    }
    out
}
