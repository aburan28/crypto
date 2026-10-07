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
    pub elkies: Vec<(u64, u64)>, // (l, t mod l)
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

/// Degree r of the irreducible factors of Phi_l(j, Y) when it has no root (Atkin prime).
pub fn atkin_degree<F: Field>(f: &F, g: &Poly<F>) -> Option<usize> {
    let g = poly::monic(f, g);
    let y = poly::x_poly(f);
    let q = f.q();
    let mut yk = y.clone();
    for k in 1..=g.len() {
        yk = poly::powmod_big(f, &yk, &q, &g);
        let d = poly::gcd(f, &g, &poly::sub(f, &yk, &y));
        if poly::deg(f, &d) > 0 {
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

/// #E(F_q) by SEA with primes l <= max_ell (the Phi_l are computed on demand into `phis`).
pub fn sea<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    max_ell: usize,
    phis: &mut HashMap<usize, Phi<F>>,
    rng: &mut Rng,
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
        let g = phi.y_poly(f, j);
        let roots = phi.neighbors(f, j, rng);
        let l = ell as u64;
        let ql = q.mod_u64(l);
        if !roots.is_empty() {
            // Elkies: kernel polynomial of a rational l-isogeny, eigenvalue
            let mut done = false;
            for jt in roots.iter().take(2) {
                let Some(et) = crate::find::elkies::elkies_codomain(f, phi, e, *jt) else { continue };
                let Some(iso) = bmss::isogeny(f, bmss::Method::FastElkiesPrime, e, &et, ell, None) else { continue };
                if let Some(lam) = elkies_eigenvalue(f, e, &iso.ker, l) {
                    let lam_inv = Int::from(lam as i64).inv_mod(&Int::from(l)).unwrap().mod_u64(l);
                    let tl = (lam + ql * lam_inv) % l;
                    st.elkies.push((l, tl));
                    // CRT
                    let (_, u, _) = Int::xgcd(&m, &Int::from(l));
                    let diff = (&Int::from(tl as i64) - &te).modulo(&Int::from(l));
                    te = &te + &(&m * &(&diff * &u).modulo(&Int::from(l)));
                    m = &m * &Int::from(l);
                    te = te.modulo(&m);
                    done = true;
                    break;
                }
            }
            if !done && roots.len() >= 2 {
                // fall back to the repeated-eigenvalue constraint is not valid here; skip prime
            }
        } else if let Some(r) = atkin_degree(f, &g) {
            let cands = atkin_candidates(l, ql, r);
            st.atkin.push((l, r, cands.len()));
            if !cands.is_empty() && (cands.len() as u64) < l {
                atkin_sets.push((l, cands));
            }
        }
        // stop once the walk over the candidates (Hasse width / M_E point additions; the Atkin
        // sets only filter the hits, they do not shorten the walk) is cheaper than another prime
        // (a modular polynomial, a kernel polynomial, a Frobenius power)
        if width.to_f64() / m.to_f64() < 65_536.0 {
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
        if count <= 1 << 22 {
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
