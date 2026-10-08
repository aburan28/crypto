//! Schoof–Elkies–Atkin point counting (Schoof 1985; Elkies 1991/1998; Atkin 1988–92), built on the
//! isogeny machinery of this crate:
//!  * l = 2: t is even iff x^3 + a x + b has a root in F_q;
//!  * Elkies primes (Phi_l(j, Y) has a root in F_q): the kernel polynomial h of the rational
//!    l-isogeny (Elkies codomain + a BMSS algorithm), then the Frobenius eigenvalue lambda on its
//!    kernel by comparing (x^q, y^q) with [lambda](x, y) in F_q[x, y]/(h(x), y^2 - f(x)), and
//!    t = lambda + q / lambda mod l;
//!  * Atkin primes (no root): the common degree r of the irreducible factors of Phi_l(j, Y) is
//!    the order of the eigenvalue ratio of Frobenius, which leaves a set of possible t mod l;
//!  * recombination: t = t_E mod M_E in the Hasse interval, candidates filtered by the Atkin sets
//!    and tested by walking [q + 1 - t] P in steps of [M_E] P (or baby-step giant-step on the
//!    progression when it is long).
use crate::bigint::Big;
use crate::curve::*;
use crate::field::{Field, Rng};
use crate::find::bmss;
use crate::find::modpoly::Phi;
use crate::int::Int;
use crate::poly::{self, Poly};
use std::collections::HashMap;

#[derive(Clone, Debug, Default)]
pub struct SeaStats {
    pub elkies: Vec<(u64, u64)>, // (l^k, t mod l^k); k > 1 from isogeny cycles
    pub atkin: Vec<(u64, usize, usize)>, // (l, r, number of t candidates mod l)
    pub m_elkies: u128,
    pub candidates: u128,
    pub tested: u64,
}

/// t mod 2 (0 if the 2-division polynomial has an F_q-root).
pub fn trace_mod2<F: Field>(f: &F, e: &Curve<F::E>) -> u64 {
    let fx = vec![e.b, e.a, f.zero(), f.one()];
    let xq = poly::powmod_big(f, &poly::x_poly(f), &f.q(), &fx);
    let g = poly::gcd(f, &fx, &poly::sub(f, &xq, &poly::x_poly(f)));
    if poly::deg(f, &g) > 0 {
        0
    } else {
        1
    }
}

/// Frobenius eigenvalue on the kernel with kernel polynomial h (Elkies prime l).
pub fn elkies_eigenvalue<F: Field>(f: &F, e: &Curve<F::E>, h: &Poly<F>, ell: u64) -> Option<u64> {
    let x = poly::x_poly(f);
    let fx = vec![e.b, e.a, f.zero(), f.one()];
    let fxh = poly::rem(f, &fx, h);
    let q = f.q();
    let xq = poly::powmod_big(f, &x, &q, h);
    let yq = poly::powmod_big(f, &fxh, &q.sub_small(1).shr(1), h);
    let one = poly::constant(f, f.one());
    let reduce = |p: &Poly<F>| poly::rem(f, p, h);
    let mulh = |a: &Poly<F>, b: &Poly<F>| poly::mulmod(f, a, b, h);
    // P = (x, y): multiples [k]P = (A_k, y B_k)
    let (mut a, mut b) = (reduce(&x), one.clone());
    let xr = reduce(&x);
    let d = (ell - 1) / 2;
    for lam in 1..=d {
        if a == xq {
            return Some(if b == yq { lam } else { ell - lam });
        }
        if lam == d {
            break;
        }
        // next multiple
        let (s, xsum) = if lam == 1 {
            // doubling: slope = y S with S = (3x^2 + a)/(2 f(x))
            let num = reduce(&poly::add(f, &poly::scale(f, &mulh(&xr, &xr), f.from_u64(3)), &poly::constant(f, e.a)));
            let den = poly::scale(f, &fxh, f.from_u64(2));
            (mulh(&num, &poly::invmod(f, &den, h)?), poly::scale(f, &xr, f.from_u64(2)))
        } else {
            // [lam]P + P: slope = y (B - 1)/(A - x)
            let num = poly::sub(f, &b, &one);
            let den = poly::sub(f, &a, &xr);
            (mulh(&num, &poly::invmod(f, &den, h)?), poly::add(f, &xr, &a))
        };
        let a3 = poly::sub(f, &mulh(&fxh, &mulh(&s, &s)), &xsum);
        let b3 = poly::sub(f, &mulh(&s, &poly::sub(f, &xr, &a3)), &one);
        a = a3;
        b = b3;
    }
    None
}

/// [k](x, y) in F_q[x, y]/(g(x), y^2 - f(x)) as (A, B) with [k](x, y) = (A, y B) (k >= 1).
pub fn ring_mul<F: Field>(f: &F, e: &Curve<F::E>, g: &Poly<F>, k: u64) -> Option<(Poly<F>, Poly<F>)> {
    let fx = poly::rem(f, &vec![e.b, e.a, f.zero(), f.one()], g);
    let one = poly::constant(f, f.one());
    let mulg = |a: &Poly<F>, b: &Poly<F>| poly::mulmod(f, a, b, g);
    let dbl = |p: &(Poly<F>, Poly<F>)| -> Option<(Poly<F>, Poly<F>)> {
        let (a, b) = p;
        // slope (3 A^2 + a)/(2 y B) = y S with S = (3 A^2 + a) / (2 f B);
        // A3 = f S^2 - 2 A; B3 = S (A - A3) - B
        let num = poly::add(f, &poly::scale(f, &mulg(a, a), f.from_u64(3)), &poly::constant(f, e.a));
        let den = poly::scale(f, &mulg(&fx, b), f.from_u64(2));
        let s = mulg(&num, &poly::invmod(f, &den, g)?);
        let a3 = poly::sub(f, &mulg(&fx, &mulg(&s, &s)), &poly::scale(f, a, f.from_u64(2)));
        let b3 = poly::sub(f, &mulg(&s, &poly::sub(f, a, &a3)), b);
        Some((a3, b3))
    };
    let add = |p: &(Poly<F>, Poly<F>), q: &(Poly<F>, Poly<F>)| -> Option<(Poly<F>, Poly<F>)> {
        if p.0 == q.0 {
            return if p.1 == q.1 { dbl(p) } else { None };
        }
        // S = (B2 - B1)/(A2 - A1); A3 = f S^2 - A1 - A2; B3 = S (A1 - A3) - B1
        let s = mulg(&poly::sub(f, &q.1, &p.1), &poly::invmod(f, &poly::sub(f, &q.0, &p.0), g)?);
        let a3 = poly::sub(f, &poly::sub(f, &mulg(&fx, &mulg(&s, &s)), &p.0), &q.0);
        let b3 = poly::sub(f, &mulg(&s, &poly::sub(f, &p.0, &a3)), &p.1);
        Some((a3, b3))
    };
    let base = (poly::rem(f, &poly::x_poly(f), g), one);
    let mut acc: Option<(Poly<F>, Poly<F>)> = None;
    let mut b = base;
    let mut k = k;
    while k > 0 {
        if k & 1 == 1 {
            acc = Some(match acc {
                None => b.clone(),
                Some(a) => add(&a, &b)?,
            });
        }
        k >>= 1;
        if k > 0 {
            b = dbl(&b)?;
        }
    }
    acc
}

/// Numerator of h(N/D) for h of degree m: sum_i h_i N^i D^(m - i).
pub fn compose_rational<F: Field>(f: &F, h: &Poly<F>, n: &Poly<F>, d: &Poly<F>) -> Poly<F> {
    let m = h.len() - 1;
    let mut npow = vec![poly::constant(f, f.one())];
    for i in 1..=m {
        npow.push(poly::mul(f, &npow[i - 1], n));
    }
    let mut dpow = vec![poly::constant(f, f.one())];
    for i in 1..=m {
        dpow.push(poly::mul(f, &dpow[i - 1], d));
    }
    let mut out: Poly<F> = vec![f.zero()];
    for i in 0..=m {
        out = poly::add(f, &out, &poly::scale(f, &poly::mul(f, &npow[i], &dpow[m - i]), h[i]));
    }
    out
}

/// Isogeny cycles (Couveignes–Morain 1994) for an Elkies prime: the Frobenius eigenvalue mod
/// l^k. phi1: E -> E1 is the rational l-isogeny on whose kernel Frobenius acts by lam1 (mod l).
/// Following the non-backtracking chain of rational l-isogenies E1 -> E2 -> ..., the preimage
/// on E of the next kernel is a Frobenius-stable cyclic subgroup of order l^i, cut out by
/// g_i = numerator of h_{psi}(x-map of the chain so far) (degree l^(i-1)(l-1)/2); Frobenius acts
/// on it by lam_i = lam_{i-1} + c l^(i-1), found by comparing x^q with x([lam_i] P) mod g_i.
/// Returns (lambda mod l^k, k) for the largest k reached (k <= k_max).
pub fn eigenvalue_cycle<F: Field>(
    f: &F,
    phi: &Phi<F>,
    e: &Curve<F::E>,
    iso1: &RatIsogeny<F>,
    ell: u64,
    lam1: u64,
    k_max: u32,
    rng: &mut Rng,
) -> (u64, u32) {
    let mut lam = lam1;
    let mut modulus = ell;
    let mut k = 1u32;
    let mut prev_j = jinv(f, e);
    let mut cur = iso1.cod;
    let (mut n, mut d) = (iso1.num.clone(), iso1.den.clone());
    let q = f.q();
    while k < k_max {
        let jc = jinv(f, &cur);
        let nbrs = phi.neighbors(f, jc, rng);
        let Some(&jn) = nbrs.iter().find(|&&r| r != prev_j) else { break };
        let Some(et) = crate::find::elkies::elkies_codomain(f, phi, &cur, jn) else { break };
        let Some(psi) = bmss::isogeny(f, bmss::Method::FastElkiesPrime, &cur, &et, ell as usize, None) else { break };
        let g = poly::monic(f, &compose_rational(f, &psi.ker, &n, &d));
        let xq = poly::powmod_big(f, &poly::x_poly(f), &q, &g);
        // [lam + c l^k] P for c = 0..l-1, by repeated addition of [l^k] P
        let Some(mut cur_pt) = ring_mul(f, e, &g, lam) else { break };
        let Some(step) = ring_mul(f, e, &g, modulus) else { break };
        let mut found = None;
        for c in 0..ell {
            if cur_pt.0 == xq {
                found = Some(lam + c * modulus);
                break;
            }
            if c + 1 < ell {
                // (A1, y B1) + (A2, y B2)
                let mulg = |a: &Poly<F>, b: &Poly<F>| poly::mulmod(f, a, b, &g);
                let fx = poly::rem(f, &vec![e.b, e.a, f.zero(), f.one()], &g);
                let Some(inv) = poly::invmod(f, &poly::sub(f, &step.0, &cur_pt.0), &g) else { break };
                let s = mulg(&poly::sub(f, &step.1, &cur_pt.1), &inv);
                let a3 = poly::sub(f, &poly::sub(f, &mulg(&fx, &mulg(&s, &s)), &cur_pt.0), &step.0);
                let b3 = poly::sub(f, &mulg(&s, &poly::sub(f, &cur_pt.0, &a3)), &cur_pt.1);
                cur_pt = (a3, b3);
            }
        }
        let Some(l2) = found else { break };
        lam = l2;
        modulus *= ell;
        k += 1;
        // extend the chain's x-map: psi_x(N/D) = psi.num(N/D) / psi.den(N/D)
        let nn = compose_rational(f, &psi.num, &n, &d);
        let dd = poly::mul(f, &compose_rational(f, &psi.den, &n, &d), &d);
        n = nn;
        d = dd;
        prev_j = jc;
        cur = psi.cod;
    }
    (lam % modulus, k)
}

/// Degree r of the irreducible factors of Phi_l(j, Y) when it has no root (Atkin prime): the
/// least k with gcd(Y^(q^k) - Y, Phi_l(j, Y)) != 1. One exponentiation gives Y^q; the further
/// Frobenius powers are modular compositions with it (`poly::Composer`).
pub fn atkin_degree<F: Field>(f: &F, g: &Poly<F>) -> Option<usize> {
    let g = poly::monic(f, g);
    let xi = poly::powmod_big(f, &poly::x_poly(f), &f.q(), &g);
    atkin_degree_with(f, &g, &xi)
}

/// `atkin_degree` for g = Phi_l(j, Y) monic, with xi = Y^q mod g already known (SEA has it from
/// the root test). The Frobenius powers Y^(q^k) come from the Q-matrix of xi (one
/// matrix-vector product each). When g is squarefree and has no root, Frobenius permutes its
/// l+1 roots like an element of a non-split torus of PGL_2(F_l), which acts freely on P^1(F_l):
/// all factors then have one degree r, r | deg g, and Y^(q^k) = Y mod g exactly when r | k, so
/// an equality at the divisors of deg g replaces the gcd at every k. Otherwise (repeated roots,
/// or no equality up to deg g) the gcd scan decides.
pub fn atkin_degree_with<F: Field>(f: &F, g: &Poly<F>, xi: &Poly<F>) -> Option<usize> {
    let y = poly::x_poly(f);
    let n = g.len() - 1;
    let has_gcd = |yk: &Poly<F>| poly::deg(f, &poly::gcd(f, g, &poly::sub(f, yk, &y))) > 0;
    if has_gcd(xi) {
        return Some(1);
    }
    let sqfree = poly::deg(f, &poly::gcd(f, g, &poly::derivative(f, g))) == 0;
    let qm = poly::Composer::with_baby_steps(f, xi, g, n);
    if sqfree {
        let mut yk = xi.clone();
        for k in 2..=n {
            yk = qm.compose(f, &yk);
            if n % k == 0 && yk == y {
                return Some(k);
            }
        }
    }
    let mut yk = xi.clone();
    for k in 2..=n {
        yk = qm.compose(f, &yk);
        if has_gcd(&yk) {
            return Some(k);
        }
    }
    None
}

/// The values t mod l for which the eigenvalue ratio of x^2 - t x + q has order r in F_{l^2}*
/// (small l: plain modular arithmetic, F_{l^2} = F_l[s]/(s^2 - nr)).
pub fn atkin_candidates(ell: u64, q_mod_l: u64, r: usize) -> Vec<u64> {
    let l = ell;
    let md = |a: u64| a % l;
    let sqrt_l = |a: u64| (0..l).find(|&x| x * x % l == md(a));
    let nr = (2..l).find(|&a| sqrt_l(a).is_none()).unwrap_or(2);
    // F_{l^2} elements (a, b) = a + b s
    let mul = |x: (u64, u64), y: (u64, u64)| (md(x.0 * y.0 + nr * md(x.1 * y.1)), md(x.0 * y.1 + x.1 * y.0));
    let pow = |mut x: (u64, u64), mut e: u64| {
        let mut r = (1u64, 0u64);
        while e > 0 {
            if e & 1 == 1 {
                r = mul(r, x);
            }
            x = mul(x, x);
            e >>= 1;
        }
        r
    };
    let inv = |x: (u64, u64)| pow(x, l * l - 2);
    let inv2 = (l + 1) / 2;
    let mut out = vec![];
    for t in 0..l {
        let disc = md(t * t + l * l * 4 - 4 * q_mod_l);
        let s = match sqrt_l(disc) {
            Some(v) => (v, 0),
            None => (0, sqrt_l(md(disc * pow((nr, 0), l * l - 2).0)).unwrap()),
        };
        let l1 = mul((md(t + s.0), s.1), (inv2, 0));
        let l2 = mul((md(t + l - s.0), md(l - s.1)), (inv2, 0));
        if l1 == (0, 0) || l2 == (0, 0) {
            continue;
        }
        let g = mul(l1, inv(l2));
        let mut x = g;
        let mut ord = 1usize;
        while x != (1, 0) && ord <= (l * l) as usize {
            x = mul(x, g);
            ord += 1;
        }
        if ord == r {
            out.push(t);
        }
    }
    out
}

/// Atkin-prime isogeny over the F_{p^d} tower. E/F_p has no rational l-isogeny, but over a large
/// enough extension F_{p^d} Frobenius^d acquires an F_l-rational eigenvalue on E[l] (the smallest
/// such d is the order of lambda over F_l, <= ord(lambda/mu)), so there E has a rational
/// l-isogeny. We scan d = 2, 3, ...: over each F_{p^d} we run the ordinary Elkies codomain +
/// kernel-polynomial + eigenvalue computation, and return the first (nu, d) for which the
/// eigenvalue step succeeds. Because elkies_eigenvalue only returns a value once Frobenius^d
/// genuinely acts as [nu] on the kernel, that nu is a true eigenvalue of pi^d, i.e. nu = lambda^d
/// or mu^d mod l. Returns None if no degree <= RMAX works. Curve over the prime field F_p.
pub fn atkin_eigenvalue_tower(zp: &crate::field::Zp, e: &Curve<u64>, ell: u64, rng: &mut Rng) -> Option<(u64, usize)> {
    let phi0 = Phi::compute(zp, ell as usize);
    atkin_eigenvalue_tower_with(zp, &phi0, e, ell, rng)
}

/// `atkin_eigenvalue_tower` with the base-field Phi_l supplied (as cached by SEA). Phi_l has
/// integer coefficients, so over F_{p^d} it is the F_p polynomial with its coefficients embedded:
/// no modular-polynomial computation happens in the extension. The roots come from factoring
/// Phi_l(j, Y) over F_p: for an irreducible factor f of degree d, F_{p^d} is built as
/// F_p[z]/(f) and z is itself a root, so no root finding happens in the extension either (its
/// conjugates z^(p^i) give the same eigenvalue). Should every such factor fail, the earlier
/// route (a fixed F_{p^d} and Cantor-Zassenhaus there) is tried for each d that has roots.
pub fn atkin_eigenvalue_tower_with(
    zp: &crate::field::Zp,
    phi0: &Phi<crate::field::Zp>,
    e: &Curve<u64>,
    ell: u64,
    rng: &mut Rng,
) -> Option<(u64, usize)> {
    use crate::fpr::{FpR, RMAX};
    let dmax = ATKIN_TOWER_MAX_DEGREE.min(RMAX);
    let g = poly::monic(zp, &phi0.y_poly(zp, jinv(zp, e)));
    // the eigenvalue over F_{p^d}, from a root jt of Phi_l(j, Y) there
    let try_root = |fr: &FpR, jt: crate::fpr::ER| -> Option<u64> {
        let e_ext = Curve::new(fr.embed(e.a), fr.embed(e.b));
        let phi = Phi { ell: phi0.ell, c: phi0.c.iter().map(|row| row.iter().map(|&x| fr.embed(x)).collect()).collect() };
        let et = crate::find::elkies::elkies_codomain(fr, &phi, &e_ext, jt)?;
        let iso = bmss::isogeny(fr, bmss::Method::FastElkiesPrime, &e_ext, &et, ell as usize, None)?;
        elkies_eigenvalue(fr, &e_ext, &iso.ker, ell)
    };
    // squarefree part (repeated roots only for special j), then distinct-degree factorisation
    let dg = poly::derivative(zp, &g);
    let sq = poly::gcd(zp, &g, &dg);
    let gs = if poly::deg(zp, &sq) > 0 { poly::monic(zp, &poly::divrem(zp, &g, &sq).0) } else { g.clone() };
    let parts = poly::ddf(zp, &gs);
    for (k, part) in &parts {
        if *k < 2 || *k > dmax {
            continue;
        }
        let mut facs = vec![];
        poly::edf(zp, part, *k, rng, &mut facs);
        for f in facs.iter().take(5) {
            let fr = FpR::from_modulus(zp.p, f);
            if let Some(lam) = try_root(&fr, fr.gen()) {
                return Some((lam, *k));
            }
        }
    }
    // fallback: roots by Cantor-Zassenhaus in a fixed F_{p^d}, for each d a factor degree divides
    for d in 2..=dmax {
        if !parts.iter().any(|(k, _)| *k >= 2 && d % *k == 0) {
            continue;
        }
        let fr = FpR::new(zp.p, d);
        let e_ext = Curve::new(fr.embed(e.a), fr.embed(e.b));
        let phi = Phi { ell: phi0.ell, c: phi0.c.iter().map(|row| row.iter().map(|&x| fr.embed(x)).collect()).collect() };
        let roots = phi.neighbors(&fr, jinv(&fr, &e_ext), rng);
        for &jt in roots.iter().take(5) {
            if let Some(lam) = try_root(&fr, jt) {
                return Some((lam, d));
            }
        }
    }
    None
}

/// t mod l for an Atkin prime, refined by the F_{p^d} eigenvalue nu = lambda^d (or mu^d) mod l:
/// of the Atkin candidates keep those whose roots lambda, mu of x^2 - t x + q (in F_{l^2})
/// satisfy lambda^d = nu or mu^d = nu. Usually a single value, i.e. an Atkin prime pinned to
/// Elkies strength.
pub fn atkin_candidates_tower(ell: u64, q_mod_l: u64, d: usize, nu: u64) -> Vec<u64> {
    let l = ell;
    let md = |a: u64| a % l;
    let sqrt_l = |a: u64| (0..l).find(|&x| x * x % l == md(a));
    let nr = (2..l).find(|&a| sqrt_l(a).is_none()).unwrap_or(2);
    let mul = |x: (u64, u64), y: (u64, u64)| (md(x.0 * y.0 + nr * md(x.1 * y.1)), md(x.0 * y.1 + x.1 * y.0));
    let pow = |mut x: (u64, u64), mut e: u64| {
        let mut rr = (1u64, 0u64);
        while e > 0 {
            if e & 1 == 1 {
                rr = mul(rr, x);
            }
            x = mul(x, x);
            e >>= 1;
        }
        rr
    };
    let mut out = vec![];
    for t in 0..l {
        let disc = md(t * t + l * l * 4 - 4 * q_mod_l);
        // genuine Atkin: discriminant a non-residue (lambda, mu not in F_l)
        if sqrt_l(disc).is_some() {
            continue;
        }
        let s = (0, sqrt_l(md(disc * pow((nr, 0), l * l - 2).0)).unwrap());
        let l1 = mul((md(t + s.0), s.1), (md((l + 1) / 2), 0));
        let l2 = mul((md(t + l - s.0), md(l - s.1)), (md((l + 1) / 2), 0));
        if pow(l1, d as u64) == (md(nu), 0) || pow(l2, d as u64) == (md(nu), 0) {
            out.push(t);
        }
    }
    out
}

/// #E(F_q) by SEA with primes l <= max_ell (the Phi_l are computed on demand into `phis`),
/// with isogeny cycles for Elkies primes up to kernel-polynomial degree `CYCLE_DEGREE`.
pub fn sea<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    max_ell: usize,
    phis: &mut HashMap<usize, Phi<F>>,
    rng: &mut Rng,
) -> Option<(Big, SeaStats)> {
    sea_opts(f, e, max_ell, phis, rng, CYCLE_DEGREE)
}

/// Unique t mod l for an Atkin prime via the F_{p^d} tower, or None if the refined candidate set
/// is not a singleton. This turns an Atkin prime into an Elkies-strength congruence for SEA.
pub fn atkin_trace_tower(fp: &crate::field::Zp, e: &Curve<u64>, ell: u64, rng: &mut Rng) -> Option<u64> {
    let (nu, d) = atkin_eigenvalue_tower(fp, e, ell, rng)?;
    let cands = atkin_candidates_tower(ell, fp.p % ell, d, nu);
    if cands.len() == 1 {
        Some(cands[0])
    } else {
        None
    }
}

/// SEA over F_p that additionally resolves Atkin primes to unique congruences through F_{p^d}
/// towers (when the tower pins t mod l), so both Elkies and Atkin primes contribute congruences.
pub fn sea_atkin_tower(
    fp: &crate::field::Zp,
    e: &Curve<u64>,
    max_ell: usize,
    phis: &mut HashMap<usize, Phi<crate::field::Zp>>,
    rng: &mut Rng,
) -> Option<(Big, SeaStats)> {
    let ee = *e;
    let fpc = fp.clone();
    let mut rng2 = Rng::new(0xA7C0 ^ fp.p);
    // the resolver cannot borrow `phis` (sea_resolved holds it), so it keeps its own base-field
    // Phi_l cache; each is computed once per prime and only embedded into the towers
    let mut cache: HashMap<u64, Phi<crate::field::Zp>> = HashMap::new();
    let mut resolver = move |l: u64| {
        let phi0 = cache.entry(l).or_insert_with(|| Phi::compute(&fpc, l as usize));
        let (nu, d) = atkin_eigenvalue_tower_with(&fpc, phi0, &ee, l, &mut rng2)?;
        let cands = atkin_candidates_tower(l, fpc.p % l, d, nu);
        if cands.len() == 1 {
            Some(cands[0])
        } else {
            None
        }
    };
    sea_resolved(fp, e, max_ell, phis, rng, CYCLE_DEGREE, Some(&mut resolver))
}

/// Stop adding primes once at most this many candidates for t remain.
pub const SEA_STOP_CANDIDATES: f64 = 4_194_304.0;
/// Candidate counts up to this are walked one point addition each; above, baby-step giant-step.
pub const SEA_WALK_MAX: u128 = 1024;

/// Default bound on the degree l^(k-1)(l-1)/2 of the polynomials used by isogeny cycles.
pub const CYCLE_DEGREE: usize = 40;

/// Largest tower degree F_{p^d} scanned for the Atkin-prime eigenvalue (Atkin primes whose
/// lambda/mu has larger order are left as candidate sets). Bounds the cost per Atkin prime.
pub const ATKIN_TOWER_MAX_DEGREE: usize = 6;

/// SEA with an explicit bound on the isogeny-cycle polynomial degree (0: no cycles; t mod l
/// only). An Elkies prime with two rational l-isogenies contributes t mod l^k for the largest k
/// with l^(k-1)(l-1)/2 <= cycle_degree.
pub fn sea_opts<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    max_ell: usize,
    phis: &mut HashMap<usize, Phi<F>>,
    rng: &mut Rng,
    cycle_degree: usize,
) -> Option<(Big, SeaStats)> {
    sea_resolved(f, e, max_ell, phis, rng, cycle_degree, None)
}

/// SEA with an optional Atkin resolver: for an Atkin prime l, `atkin_resolver(l)` may return a
/// unique t mod l (e.g. from the F_{p^d} tower), which is CRT'd in as a congruence instead of
/// being kept only as a filtering set. Generic in F; the resolver (when given) is specific to the
/// base field it closes over.
#[allow(clippy::too_many_arguments)]
pub fn sea_resolved<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    max_ell: usize,
    phis: &mut HashMap<usize, Phi<F>>,
    rng: &mut Rng,
    cycle_degree: usize,
    mut atkin_resolver: Option<&mut dyn FnMut(u64) -> Option<u64>>,
) -> Option<(Big, SeaStats)> {
    let q = Int::from_big(&f.q());
    let j = jinv(f, e);
    assert!(j != f.zero() && j != f.from_u64(1728), "SEA here needs j != 0, 1728");
    let mut st = SeaStats::default();
    // Elkies congruences (t mod m) and Atkin sets
    let mut te = Int::from(trace_mod2(f, e) as i64);
    let mut m = Int::from(2i64);
    let mut atkin_sets: Vec<(u64, Vec<u64>)> = vec![];
    let hasse = (&q * &Int::from(4i64)).isqrt(); // |t| <= 2 sqrt(q)
    let width = &hasse * &Int::from(2i64) + Int::one();
    let mut ell = 3usize;
    while ell <= max_ell {
        if !crate::field::is_prime(ell as u64) {
            ell += 2;
            continue;
        }
        let phi = phis.entry(ell).or_insert_with(|| Phi::compute(f, ell));
        // roots of Phi_l(j, Y) as in Phi::neighbors (gcd with Y^q - Y, then splitting), keeping
        // Y^q mod Phi_l(j, Y) for the Atkin degree below
        let g = poly::monic(f, &phi.y_poly(f, j));
        let xi = poly::powmod_big(f, &poly::x_poly(f), &f.q(), &g);
        let mut roots = vec![];
        poly::split_roots(f, &poly::gcd(f, &g, &poly::sub(f, &xi, &poly::x_poly(f))), rng, &mut roots);
        let l = ell as u64;
        let ql = q.mod_u64(l);
        if !roots.is_empty() {
            // Elkies: kernel polynomial of a rational l-isogeny, eigenvalue
            let mut done = false;
            for jt in roots.iter().take(2) {
                let Some(et) = crate::find::elkies::elkies_codomain(f, phi, e, *jt) else { continue };
                let Some(iso) = bmss::isogeny(f, bmss::Method::FastElkiesPrime, e, &et, ell, None) else { continue };
                if let Some(lam1) = elkies_eigenvalue(f, e, &iso.ker, l) {
                    // isogeny cycles: eigenvalue mod l^k when the eigenvalues are distinct
                    let mut kmax = 1u32;
                    while roots.len() == 2 && (l.pow(kmax) * (l - 1) / 2) as usize <= cycle_degree {
                        kmax += 1;
                    }
                    let (lam, k) = if kmax > 1 { eigenvalue_cycle(f, phi, e, &iso, l, lam1, kmax, rng) } else { (lam1, 1) };
                    let lk = l.pow(k);
                    let lam_inv = Int::from(lam as i64).inv_mod(&Int::from(lk)).unwrap().mod_u64(lk);
                    let tl = ((lam as u128 + q.mod_u64(lk) as u128 * lam_inv as u128) % lk as u128) as u64;
                    st.elkies.push((lk, tl));
                    // CRT
                    let (_, u, _) = Int::xgcd(&m, &Int::from(lk));
                    let diff = (&Int::from(tl as i64) - &te).modulo(&Int::from(lk));
                    te = &te + &(&m * &(&diff * &u).modulo(&Int::from(lk)));
                    m = &m * &Int::from(lk);
                    te = te.modulo(&m);
                    done = true;
                    break;
                }
            }
            if !done && roots.len() >= 2 {
                // fall back to the repeated-eigenvalue constraint is not valid here; skip prime
            }
        } else if let Some(r) = atkin_degree_with(f, &g, &xi) {
            // tower resolution to a unique congruence (Elkies-strength from an Atkin prime)
            let resolved = atkin_resolver.as_deref_mut().and_then(|res| res(l));
            if let Some(tl) = resolved {
                st.elkies.push((l, tl));
                let (_, u, _) = Int::xgcd(&m, &Int::from(l as i64));
                let diff = (&Int::from(tl as i64) - &te).modulo(&Int::from(l as i64));
                te = &te + &(&m * &(&diff * &u).modulo(&Int::from(l as i64)));
                m = &m * &Int::from(l as i64);
                te = te.modulo(&m);
            } else {
                let cands = atkin_candidates(l, ql, r);
                st.atkin.push((l, r, cands.len()));
                if !cands.is_empty() && (cands.len() as u64) < l {
                    atkin_sets.push((l, cands));
                }
            }
        }
        // stop once the search over the candidates (Hasse width / M_E of them, ~2 sqrt(.) point
        // additions by baby-step giant-step; the Atkin sets only filter the hits) is cheaper
        // than another prime (roots of Phi_l(j, Y), a kernel polynomial, a Frobenius power)
        if width.to_f64() / m.to_f64() < SEA_STOP_CANDIDATES {
            break;
        }
        ell += 2;
    }
    st.m_elkies = m.to_i128().map(|v| v as u128).unwrap_or(u128::MAX);
    // candidates t = te + k m in [-hasse, hasse]: k from ceil((-hasse - te)/m) to floor((hasse - te)/m)
    let k_lo = -(&hasse + &te).div_floor(&m);
    let k_hi = (&hasse - &te).div_floor(&m);
    if k_hi < k_lo {
        return None;
    }
    let count = (&(&k_hi - &k_lo) + &Int::one()).to_i128()? as u128;
    st.candidates = count;
    let q1 = &q + &Int::one();
    let to_big = |v: &Int| -> Big { v.mag().clone() };
    let ok_atkin = |t: &Int| atkin_sets.iter().all(|(l, s)| s.contains(&t.mod_u64(*l)));
    let neg = |p: Pt<F::E>| match p {
        Pt::Inf => Pt::Inf,
        Pt::Aff(x, y) => Pt::Aff(x, f.neg(y)),
    };
    // t with [q + 1 - t] P = 0 for every sampled P (intersection over points)
    let mut cands: Option<Vec<Int>> = None;
    for _attempt in 0..12 {
        let p = random_point_f(f, e, rng);
        let mut hits = vec![];
        let t0 = &te + &(&k_lo * &m);
        if count <= SEA_WALK_MAX {
            // walk R_k = [q + 1 - t_k] P, t_k = t0 + k m: R_{k+1} = R_k - [m] P
            let mut r = pmul_big(f, e, &p, &to_big(&(&q1 - &t0)));
            let neg_mp = neg(pmul_big(f, e, &p, &to_big(&m)));
            let mut t = t0;
            for _ in 0..count {
                if r == Pt::Inf && ok_atkin(&t) {
                    hits.push(t.clone());
                }
                st.tested += 1;
                r = padd(f, e, &r, &neg_mp);
                t = &t + &m;
            }
        } else {
            // baby-step giant-step: [q+1-t0] P = [k m] P, k in [0, count)
            let base = pmul_big(f, e, &p, &to_big(&(&q1 - &t0)));
            let mp = pmul_big(f, e, &p, &to_big(&m));
            let bs = (count as f64).sqrt().ceil() as u128;
            let mut baby: HashMap<Pt<F::E>, u128> = HashMap::new();
            let mut cur = Pt::Inf;
            for jj in 0..bs {
                baby.entry(cur).or_insert(jj);
                cur = padd(f, e, &cur, &mp);
            }
            let ngiant = neg(cur);
            let mut g = base;
            let mut i = 0u128;
            while i * bs < count {
                if let Some(&jj) = baby.get(&g) {
                    let kp = i * bs + jj;
                    if kp < count {
                        let t = &t0 + &(&Int::from(kp as i128) * &m);
                        if ok_atkin(&t) {
                            hits.push(t);
                        }
                    }
                }
                g = padd(f, e, &g, &ngiant);
                i += 1;
                st.tested += 1;
            }
        }
        let next: Vec<Int> = match cands {
            None => hits,
            Some(c) => c.into_iter().filter(|t| hits.contains(t)).collect(),
        };
        if next.len() == 1 {
            cands = Some(next);
            break;
        }
        cands = Some(next);
    }
    let c = cands?;
    let found = if c.len() == 1 { Some(c[0].clone()) } else { None };
    let t = found?;
    Some(((&q1 - &t).mag().clone(), st))
}
