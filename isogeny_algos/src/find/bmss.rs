//! The family of algorithms that turn a pair of l-isogenous curves E, E~ (and, for some of them, the
//! sum sigma of the kernel abscissas) into the normalised isogeny, as surveyed in
//! Bostan–Morain–Salvy–Schost, "Fast algorithms for computing isogenies between elliptic curves",
//! Math. Comp. 77 (2008), Table 1:
//!
//! | method          | complexity (paper) | needs sigma |
//! |-----------------|--------------------|-------------|
//! | LinearAlgebra   | O(l^omega)         | no          |
//! | Stark1972       | O(l M(l))          | no          |
//! | Atkin1992       | O(l M(l))          | yes         |
//! | AtkinModComp    | O(M(l) sqrt(l)+..) | yes         |
//! | Elkies1992      | O(l^2)             | yes         |
//! | Elkies1998      | O(l^2)             | yes         |
//! | FastElkies      | O(M(l))            | yes         |
//! | FastElkies'     | O(M(l) log l)      | no          |
//!
//! Everything here is the quadratic/cubic *baseline* (schoolbook series, dense linear algebra): the
//! asymptotic columns are the targets for later optimisation, not what this code achieves.
//! Odd l only; characteristic must exceed ~8l. Each method returns the kernel polynomial g
//! (D = g^2), from which the isogeny is assembled with Kohel's formula and cross-checked: the
//! codomain must equal E~.
use super::elkies::{solve_linear, u_series};
use super::modpoly::Phi;
use crate::curve::*;
use crate::field::Field;
use crate::kernel::kohel::kohel;
use crate::poly::Poly;
use crate::series;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    LinearAlgebra,
    Stark1972,
    Atkin1992,
    AtkinModComp,
    Elkies1992,
    Elkies1998,
    FastElkies,
    FastElkiesPrime,
}

impl Method {
    pub const ALL: [Method; 8] = [
        Method::LinearAlgebra,
        Method::Stark1972,
        Method::Atkin1992,
        Method::AtkinModComp,
        Method::Elkies1992,
        Method::Elkies1998,
        Method::FastElkies,
        Method::FastElkiesPrime,
    ];
    pub fn needs_sigma(self) -> bool {
        !matches!(
            self,
            Method::LinearAlgebra | Method::Stark1972 | Method::FastElkiesPrime
        )
    }
    pub fn name(self) -> &'static str {
        match self {
            Method::LinearAlgebra => "linear_algebra",
            Method::Stark1972 => "stark1972",
            Method::Atkin1992 => "atkin1992",
            Method::AtkinModComp => "atkin_modcomp",
            Method::Elkies1992 => "elkies1992",
            Method::Elkies1998 => "elkies1998",
            Method::FastElkies => "fast_elkies",
            Method::FastElkiesPrime => "fast_elkies_prime",
        }
    }
}

/// Coefficients c_k of wp(z) = z^-2 + sum_{k>=1} c_k z^{2k}, k = 0..=n (c_0 unused = 0).
fn c_coeffs<F: Field>(f: &F, e: &Curve<F::E>, n: usize) -> Vec<F::E> {
    let u = u_series(f, e, n + 1);
    let mut c = vec![f.zero(); n + 1];
    for k in 1..=n {
        c[k] = u[k + 1];
    }
    c
}

/// V(Z) = 1 + sum_{k>=2} u_k Z^k with wp = Z^-1 V(Z), Z = z^2, to length n.
fn v_series<F: Field>(f: &F, e: &Curve<F::E>, n: usize) -> Vec<F::E> {
    let u = u_series(f, e, n);
    u[..n].to_vec()
}

/// Monic polynomial of degree d whose first d power sums are q[1..=d] (Newton identities).
fn poly_from_power_sums<F: Field>(f: &F, q: &[F::E], d: usize) -> Poly<F> {
    let mut e = vec![f.one()];
    for k in 1..=d {
        let mut s = f.zero();
        for i in 1..=k {
            let t = f.mul(e[k - i], q[i]);
            s = if i % 2 == 1 { f.add(s, t) } else { f.sub(s, t) };
        }
        e.push(f.div(s, f.from_u64(k as u64)));
    }
    let mut g = vec![f.zero(); d + 1];
    for k in 0..=d {
        g[d - k] = if k % 2 == 0 { e[k] } else { f.neg(e[k]) };
    }
    g
}

/// g with g^2 = D (D monic of even degree 2d); None if D is not a square.
fn sqrt_monic<F: Field>(f: &F, dd: &Poly<F>) -> Option<Poly<F>> {
    if dd.len() % 2 == 0 {
        return None;
    }
    let d = (dd.len() - 1) / 2;
    let lead = *dd.last()?;
    let li = f.inv(lead);
    let dm: Vec<F::E> = dd.iter().map(|&c| f.mul(c, li)).collect();
    let mut g = vec![f.zero(); d + 1];
    g[d] = f.one();
    let inv2 = f.inv(f.from_u64(2));
    for k in 1..=d {
        let mut s = dm[2 * d - k];
        for i in 1..k {
            s = f.sub(s, f.mul(g[d - i], g[d - k + i]));
        }
        g[d - k] = f.mul(s, inv2);
    }
    if crate::poly::mul(f, &g, &g) != dm {
        return None;
    }
    Some(g)
}

/// h_1..h_count of N/D = x + sum h_i x^-i by the recurrence (16) of BMSS (Elkies 1998).
fn h_recurrence<F: Field>(f: &F, e: &Curve<F::E>, et: &Curve<F::E>, count: usize) -> Vec<F::E> {
    let c = |v: u64| f.from_u64(v);
    let mut h = vec![f.zero(); count.max(3) + 1];
    h[1] = f.div(f.sub(e.a, et.a), c(5));
    h[2] = f.div(f.sub(e.b, et.b), c(7));
    for k in 3..=count {
        let mut s = f.zero();
        for i in 1..=k - 2 {
            s = f.add(s, f.mul(h[i], h[k - 1 - i]));
        }
        let den = c(((k - 2) * (2 * k + 3)) as u64);
        let t1 = f.div(f.mul(c(3), s), den);
        let t2 = f.mul(
            f.div(c((2 * k - 3) as u64), c((2 * k + 3) as u64)),
            f.mul(e.a, h[k - 2]),
        );
        let t3 = if k >= 4 {
            f.mul(
                f.div(c((2 * (k - 3)) as u64), c((2 * k + 3) as u64)),
                f.mul(e.b, h[k - 3]),
            )
        } else {
            f.zero()
        };
        h[k] = f.sub(f.sub(t1, t2), t3);
    }
    h
}

/// Power sums q_0..=q_d of g from h_1..h_{d-1} and sigma, via (18):
/// h_i = (4i+2) q_{i+1} + (4i-2) A q_{i-1} + (4i-4) B q_{i-2}.
fn q_from_h<F: Field>(f: &F, e: &Curve<F::E>, d: usize, sigma: F::E, h: &[F::E]) -> Vec<F::E> {
    let c = |v: u64| f.from_u64(v);
    let mut q = vec![f.zero(); d + 2];
    q[0] = c(d as u64);
    q[1] = f.div(sigma, c(2));
    for i in 1..d {
        let mut s = f.sub(h[i], f.mul(c((4 * i - 2) as u64), f.mul(e.a, q[i - 1])));
        if i >= 2 {
            s = f.sub(s, f.mul(c((4 * i - 4) as u64), f.mul(e.b, q[i - 2])));
        }
        q[i + 1] = f.div(s, c((4 * i + 2) as u64));
    }
    q
}

fn elkies1998<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    et: &Curve<F::E>,
    ell: usize,
    sigma: F::E,
) -> Option<Poly<F>> {
    let d = (ell - 1) / 2;
    let h = h_recurrence(f, e, et, d.max(2));
    let q = q_from_h(f, e, d, sigma, &h);
    Some(poly_from_power_sums(f, &q, d))
}

fn elkies1992<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    et: &Curve<F::E>,
    ell: usize,
    sigma: F::E,
) -> Option<Poly<F>> {
    let c = |v: u64| f.from_u64(v);
    let d = (ell - 1) / 2;
    let (ck, ct) = (c_coeffs(f, e, d + 1), c_coeffs(f, et, d + 1));
    // mu[k][j]: d^{2k} wp / dz^{2k} = sum_j mu[k][j] wp^j
    let mut mu: Vec<Vec<F::E>> = vec![vec![f.zero(), f.one(), f.zero()]];
    for k in 0..d.saturating_sub(1) {
        let prev = &mu[k];
        let at = |j: i64| {
            if j >= 0 && (j as usize) < prev.len() {
                prev[j as usize]
            } else {
                f.zero()
            }
        };
        let mut row = vec![f.zero(); k + 3];
        for (j, slot) in row.iter_mut().enumerate() {
            let ji = j as i64;
            let a1 = f.mul(c(((2 * ji - 2) * (2 * ji - 1)).unsigned_abs()), at(ji - 1));
            let a1 = if (2 * ji - 2) * (2 * ji - 1) < 0 {
                f.neg(a1)
            } else {
                a1
            };
            let a2 = f.mul(
                c(((2 * ji + 1) * (2 * ji + 2)) as u64),
                f.mul(e.a, at(ji + 1)),
            );
            let a3 = f.mul(
                c(((2 * ji + 2) * (2 * ji + 4)) as u64),
                f.mul(e.b, at(ji + 2)),
            );
            *slot = f.add(f.add(a1, a2), a3);
        }
        mu.push(row);
    }
    let mut q = vec![f.zero(); d + 2];
    q[0] = c(d as u64);
    q[1] = f.div(sigma, c(2));
    let mut fact = f.one(); // (2k)!
    let mut fact_next = f.one();
    let _ = &mut fact_next;
    for k in 1..d {
        fact = f.mul(fact, f.mul(c((2 * k - 1) as u64), c((2 * k) as u64)));
        let mut s = f.div(f.mul(fact, f.sub(ct[k], ck[k])), c(2));
        for j in 0..=k {
            s = f.sub(s, f.mul(mu[k][j], q[j]));
        }
        q[k + 1] = f.div(s, mu[k][k + 1]);
    }
    Some(poly_from_power_sums(f, &q, d))
}

/// F(Z) = -sigma Z + 2 sum_{k>=1} (l c_k - c~_k) Z^{k+1} / ((2k+1)(2k+2)), length n.
fn atkin_f<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    et: &Curve<F::E>,
    ell: usize,
    sigma: F::E,
    n: usize,
) -> Vec<F::E> {
    let c = |v: u64| f.from_u64(v);
    let (ck, ct) = (c_coeffs(f, e, n), c_coeffs(f, et, n));
    let mut fs = vec![f.zero(); n];
    if n > 1 {
        fs[1] = f.neg(sigma);
    }
    for k in 1..n.saturating_sub(1) {
        let num = f.mul(c(2), f.sub(f.mul(c(ell as u64), ck[k]), ct[k]));
        fs[k + 1] = f.add(fs[k + 1], f.div(num, c(((2 * k + 1) * (2 * k + 2)) as u64)));
    }
    fs
}

fn atkin1992<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    et: &Curve<F::E>,
    ell: usize,
    sigma: F::E,
) -> Option<Poly<F>> {
    let gser = series::exp(f, &atkin_f(f, e, et, ell, sigma, ell), ell);
    let v = v_series(f, e, ell + 1);
    // t[idx] = coefficient of Z^{idx-(l-1)} in T = Z^{1-l} G(Z)
    let mut t = gser.clone();
    let mut d = vec![f.zero(); ell];
    // powers of V
    let mut vp: Vec<Vec<F::E>> = vec![vec![f.one()]];
    let mut cur = {
        let mut x = vec![f.zero(); ell + 1];
        x[0] = f.one();
        x
    };
    for _ in 1..ell {
        cur = series::mul(f, &cur, &v, ell + 1);
        vp.push(cur.clone());
    }
    for i in (0..ell).rev() {
        let ti = t[ell - 1 - i];
        d[i] = ti;
        if !f.is_zero(ti) {
            // wp^i = Z^{-i} V^i : exponent -i + j  <->  idx = (ell-1) - i + j
            for j in 0..=i {
                let idx = ell - 1 - i + j;
                t[idx] = f.sub(t[idx], f.mul(ti, vp[i][j]));
            }
        }
    }
    if t.iter().any(|&x| !f.is_zero(x)) {
        return None;
    }
    sqrt_monic(f, &d)
}

fn atkin_modcomp<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    et: &Curve<F::E>,
    ell: usize,
    sigma: F::E,
) -> Option<Poly<F>> {
    let c = |v: u64| f.from_u64(v);
    let n = ell + 2;
    let gser = series::exp(f, &atkin_f(f, e, et, ell, sigma, n), n);
    // J(x) = x^{-1/2} I(x), recurrence (23)
    let mut a = vec![f.zero(); n + 1];
    a[0] = f.one();
    a[2] = f.neg(f.div(e.a, c(10)));
    for i in 2..n {
        let mut s = f.mul(c((2 * i - 3) as u64), f.mul(e.b, a[i - 2]));
        s = f.add(s, f.mul(f.mul(c(2), f.mul(e.a, c(i as u64))), a[i - 1]));
        let coef = f.div(c((2 * i - 1) as u64), c((2 * (i + 1) * (2 * i + 3)) as u64));
        a[i + 1] = f.neg(f.mul(coef, s));
    }
    let jsq = series::mul(f, &a, &a, n);
    // I^2 = x J^2
    let mut z = vec![f.zero(); n];
    z[1..n].copy_from_slice(&jsq[..n - 1]);
    let gc = series::compose(f, &gser, &z, n);
    let jinvsq = series::inv(f, &jsq, n);
    let jpow = series::pow(f, &jinvsq, ell - 1, n);
    let r = series::mul(f, &gc, &jpow, n);
    // D(1/x) = x^{1-l} R(x)
    let mut dd = vec![f.zero(); ell];
    for m in 0..ell {
        dd[ell - 1 - m] = r[m];
    }
    sqrt_monic(f, &dd)
}

/// S(x) = x T(x^2) solving (B x^6 + A x^4 + 1) S'^2 = 1 + A~ S^4 + B~ S^6, S = x + O(x^5), to length n.
fn s_ode<F: Field>(f: &F, e: &Curve<F::E>, et: &Curve<F::E>, n: usize) -> Vec<F::E> {
    let c = |v: u64| f.from_u64(v);
    let mut den = vec![f.zero(); n];
    den[0] = f.one();
    if n > 4 {
        den[4] = e.a;
    }
    if n > 6 {
        den[6] = e.b;
    }
    let cc = series::inv(f, &den, n);
    let mut s = vec![f.zero(); n];
    s[1] = f.one();
    for _ in 0..n {
        let s4 = series::pow(f, &s, 4, n);
        let s6 = series::pow(f, &s, 6, n);
        let mut inner = vec![f.zero(); n];
        inner[0] = f.one();
        for k in 0..n {
            inner[k] = f.add(inner[k], f.add(f.mul(et.a, s4[k]), f.mul(et.b, s6[k])));
        }
        let rr = series::mul(f, &cc, &inner, n);
        // sqrt(rr), rr[0] = 1
        let mut sp = vec![f.zero(); n];
        sp[0] = f.one();
        let inv2 = f.inv(c(2));
        for k in 1..n {
            let mut t = rr[k];
            for i in 1..k {
                t = f.sub(t, f.mul(sp[i], sp[k - i]));
            }
            sp[k] = f.mul(t, inv2);
        }
        let mut ns = vec![f.zero(); n];
        for k in 0..n - 1 {
            ns[k + 1] = f.div(sp[k], c((k + 1) as u64));
        }
        if ns == s {
            break;
        }
        s = ns;
    }
    s
}

/// U(w) = 1/T(w)^2 where S(x) = x T(x^2); N/D = x U(1/x) so h_i = U[i+1]; length m.
fn u_from_ode<F: Field>(f: &F, e: &Curve<F::E>, et: &Curve<F::E>, m: usize) -> Vec<F::E> {
    let s = s_ode(f, e, et, 2 * m + 2);
    let t: Vec<F::E> = (0..m).map(|i| s[2 * i + 1]).collect();
    let t2 = series::mul(f, &t, &t, m);
    series::inv(f, &t2, m)
}

fn fast_elkies<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    et: &Curve<F::E>,
    ell: usize,
    sigma: F::E,
) -> Option<Poly<F>> {
    let d = (ell - 1) / 2;
    let u = u_from_ode(f, e, et, d + 3);
    let mut h = vec![f.zero(); d.max(2) + 1];
    for i in 1..h.len() {
        h[i] = u[i + 1];
    }
    let q = q_from_h(f, e, d, sigma, &h);
    Some(poly_from_power_sums(f, &q, d))
}

fn fast_elkies_prime<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    et: &Curve<F::E>,
    ell: usize,
) -> Option<Poly<F>> {
    let m = 2 * ell + 4;
    let u = u_from_ode(f, e, et, m);
    // U = Nh/Dh with deg Nh <= l, deg Dh <= l-1, Dh(0) = 1: sum_{i=1}^{l-1} dh_i U_{k-i} = -U_k, k = l+1..2l-1
    let r = ell - 1;
    let mut mat = vec![vec![f.zero(); r]; r];
    let mut rhs = vec![f.zero(); r];
    for k in ell + 1..2 * ell {
        for i in 1..=r {
            mat[k - ell - 1][i - 1] = u[k - i];
        }
        rhs[k - ell - 1] = f.neg(u[k]);
    }
    let sol = solve_linear(f, mat, rhs)?;
    // D(x) = x^{l-1} Dh(1/x)
    let mut dd = vec![f.zero(); ell];
    dd[ell - 1] = f.one();
    for i in 1..=r {
        dd[ell - 1 - i] = sol[i - 1];
    }
    sqrt_monic(f, &dd)
}

fn linear_algebra<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    et: &Curve<F::E>,
    ell: usize,
) -> Option<Poly<F>> {
    // Z^l N(wp) = sum_i n_i Z^{l-i} V^i ; Z^l wp~ D(wp) = V~ sum_i d_i Z^{l-1-i} V^i.   (n_l = d_{l-1} = 1)
    let n = 2 * ell + 1;
    let v = v_series(f, e, n);
    let vt = v_series(f, et, n);
    let mut vp: Vec<Vec<F::E>> = vec![{
        let mut x = vec![f.zero(); n];
        x[0] = f.one();
        x
    }];
    for i in 1..=ell {
        let nx = series::mul(f, &vp[i - 1], &v, n);
        vp.push(nx);
    }
    // unknowns: n_0..n_{l-1} (indices 0..l), d_0..d_{l-2} (indices l..2l-1)
    let unk = 2 * ell - 1;
    let mut mat = vec![vec![f.zero(); unk]; unk];
    let mut rhs = vec![f.zero(); unk];
    // residual coefficient of Z^t (t = 1..=2l-1) of: Z^0 V^l [n_l] - V~ V^{l-1} [d_{l-1}]  + sum unknown terms
    let vt_vl1 = series::mul(f, &vt, &vp[ell - 1], n);
    for t in 1..=unk {
        let row = t - 1;
        // known part: V^l - V~ V^{l-1}
        rhs[row] = f.neg(f.sub(vp[ell][t], vt_vl1[t]));
        for i in 0..ell {
            // n_i Z^{l-i} V^i  -> coefficient at Z^t is V^i[t-(l-i)]
            if t + i >= ell {
                mat[row][i] = vp[i][t + i - ell];
            }
        }
        for i in 0..ell - 1 {
            // - d_i V~ Z^{l-1-i} V^i
            let prod = series::mul(f, &vt, &vp[i], n);
            if t + 1 + i >= ell {
                mat[row][ell + i] = f.neg(prod[t + 1 + i - ell]);
            }
        }
    }
    let sol = solve_linear(f, mat, rhs)?;
    let mut dd = vec![f.zero(); ell];
    dd[ell - 1] = f.one();
    for i in 0..ell - 1 {
        dd[i] = sol[ell + i];
    }
    sqrt_monic(f, &dd)
}

/// Laurent series helper for Stark: coefficients for exponents base.. (absolute exponent bound kmax).
fn stark<F: Field>(f: &F, e: &Curve<F::E>, et: &Curve<F::E>, ell: usize) -> Option<Poly<F>> {
    let k0 = 4 * ell + 16;
    let v = v_series(f, e, k0 + ell + 2);
    let vt = v_series(f, et, k0 + ell + 2);
    // powers of V up to degree l
    let mut vp: Vec<Vec<F::E>> = vec![{
        let mut x = vec![f.zero(); k0 + ell + 2];
        x[0] = f.one();
        x
    }];
    for i in 1..=ell + 1 {
        let nx = series::mul(f, &vp[i - 1], &v, k0 + ell + 2);
        vp.push(nx);
    }
    // T = wp~ = Z^-1 V~ : base = -1, exponents -1 .. k0-1
    let mut base: i64 = -1;
    let mut kmax: i64 = k0 as i64;
    let mut t: Vec<F::E> = vt[..(kmax - base) as usize].to_vec();
    let mut q_prev2: Poly<F> = vec![f.one()];
    let mut q_prev1: Poly<F> = vec![];
    let mut iterations = 0;
    loop {
        if q_prev1.len() >= ell {
            break; // deg q_n >= l-1
        }
        iterations += 1;
        if iterations > 4 * ell + 10 || kmax <= 2 {
            return None;
        }
        // polynomial part a(x) of T in terms of wp
        let lowest = base.min(0);
        let r = (-lowest).max(0) as usize;
        let mut a = vec![f.zero(); r + 1];
        for m in (0..=r).rev() {
            let idx = (-(m as i64) - base) as usize;
            let c = if idx < t.len() { t[idx] } else { f.zero() };
            if f.is_zero(c) {
                continue;
            }
            a[m] = c;
            // subtract c * wp^m = c Z^{-m} V^m
            for j in 0..vp[m].len() {
                let ex = -(m as i64) + j as i64;
                if ex >= kmax {
                    break;
                }
                let idx = (ex - base) as usize;
                t[idx] = f.sub(t[idx], f.mul(c, vp[m][j]));
            }
        }
        crate::poly::trim(f, &mut a);
        // q_n = a q_{n-1} + q_{n-2}
        let q_new = crate::poly::add(f, &crate::poly::mul(f, &a, &q_prev1), &q_prev2);
        q_prev2 = q_prev1;
        q_prev1 = q_new;
        if q_prev1.len() >= ell {
            break;
        }
        // T <- 1 / T  (T now has valuation s >= 1)
        let first = t.iter().position(|&c| !f.is_zero(c))?;
        let s = base + first as i64;
        let w: Vec<F::E> = t[first..].to_vec();
        if w.is_empty() {
            return None;
        }
        let rel = w.len();
        t = series::inv(f, &w, rel);
        base = -s;
        kmax = base + rel as i64;
    }
    if q_prev1.len() != ell {
        return None;
    }
    let li = f.inv(*q_prev1.last()?);
    let dm: Vec<F::E> = q_prev1.iter().map(|&c| f.mul(c, li)).collect();
    sqrt_monic(f, &dm)
}

/// Run one algorithm of the family; None if it fails (degenerate input, or sigma missing).
pub fn isogeny<F: Field>(
    f: &F,
    m: Method,
    e: &Curve<F::E>,
    et: &Curve<F::E>,
    ell: usize,
    sigma: Option<F::E>,
) -> Option<RatIsogeny<F>> {
    if ell < 3 || ell % 2 == 0 {
        return None;
    }
    let g = match m {
        Method::LinearAlgebra => linear_algebra(f, e, et, ell),
        Method::Stark1972 => stark(f, e, et, ell),
        Method::Atkin1992 => atkin1992(f, e, et, ell, sigma?),
        Method::AtkinModComp => atkin_modcomp(f, e, et, ell, sigma?),
        Method::Elkies1992 => elkies1992(f, e, et, ell, sigma?),
        Method::Elkies1998 => elkies1998(f, e, et, ell, sigma?),
        Method::FastElkies => fast_elkies(f, e, et, ell, sigma?),
        Method::FastElkiesPrime => fast_elkies_prime(f, e, et, ell),
    }?;
    let iso = kohel(f, e, &g, ell as u64);
    if iso.cod != *et {
        return None;
    }
    Some(iso)
}

/// sigma = sum of the abscissas of the non-zero kernel points, from the second-order relation of
/// Phi_l (Elkies): differentiating Phi(j(tau), j(tau')) = 0 twice and using Ramanujan's identities
/// isolates E2 - l E2', and  sigma = (E2' - l E2)/12  (E~ must be the *normalised* codomain).
pub fn sigma_from_phi<F: Field>(
    f: &F,
    phi: &Phi,
    e: &Curve<F::E>,
    et: &Curve<F::E>,
) -> Option<F::E> {
    let c = |v: u64| f.from_u64(v);
    let ell = c(phi.ell as u64);
    let (j, jt) = (jinv(f, e), jinv(f, et));
    if f.is_zero(e.a) || f.is_zero(et.a) || f.is_zero(e.b) || f.is_zero(et.b) {
        return None;
    }
    let (e4, e6) = (f.mul(f.neg(c(48)), e.a), f.mul(f.neg(c(864)), e.b));
    let (e4t, e6t) = (f.mul(f.neg(c(48)), et.a), f.mul(f.neg(c(864)), et.b));
    let (r, rt) = (f.div(e6, e4), f.div(e6t, e4t));
    let (px, py) = phi.grad(f, j, jt);
    let (pxx, pxy, pyy) = phi.hessian(f, j, jt);
    let dj = f.neg(f.mul(j, r));
    let djt = f.neg(f.div(f.mul(jt, rt), ell));
    let w = f.mul(px, f.mul(j, r));
    if f.is_zero(w) {
        return None;
    }
    let q = f.add(
        f.add(
            f.mul(pxx, f.mul(dj, dj)),
            f.mul(f.mul(c(2), pxy), f.mul(dj, djt)),
        ),
        f.mul(pyy, f.mul(djt, djt)),
    );
    let half = f.inv(c(2));
    let third2 = f.div(c(2), c(3));
    let t1 = f.mul(
        px,
        f.mul(j, f.add(f.mul(e4, half), f.mul(third2, f.mul(r, r)))),
    );
    let t2 = f.div(
        f.mul(
            py,
            f.mul(jt, f.add(f.mul(e4t, half), f.mul(third2, f.mul(rt, rt)))),
        ),
        f.mul(ell, ell),
    );
    let la_minus_ap = f.div(f.mul(f.mul(c(6), ell), f.add(q, f.add(t1, t2))), w);
    Some(f.neg(f.div(la_minus_ap, c(12))))
}

/// All F-rational normalised l-isogenies out of `e`, end to end: roots of Phi_l(j, Y), Elkies'
/// normalised codomain for each, sigma from the second-order relation of Phi (when the method
/// needs it), then the chosen algorithm of the family.
pub fn isogenies_via_phi<F: Field>(
    f: &F,
    phi: &Phi,
    e: &Curve<F::E>,
    m: Method,
    rng: &mut crate::field::Rng,
) -> Vec<RatIsogeny<F>> {
    let mut out = vec![];
    for jt in phi.neighbors(f, jinv(f, e), rng) {
        let Some(et) = super::elkies::elkies_codomain(f, phi, e, jt) else {
            continue;
        };
        let sigma = if m.needs_sigma() {
            match sigma_from_phi(f, phi, e, &et) {
                Some(s) => Some(s),
                None => continue,
            }
        } else {
            None
        };
        if let Some(iso) = isogeny(f, m, e, &et, phi.ell, sigma) {
            out.push(iso);
        }
    }
    out
}
