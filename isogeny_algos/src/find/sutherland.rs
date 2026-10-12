//! Sutherland-style modular polynomials (Sutherland 2011, "Computing modular polynomials with the
//! Chinese remainder theorem"): compute Phi_l mod p from the l-isogeny graph, then CRT to the
//! integer polynomial. The mod-p step is independent of the q-expansion methods in `modpoly.rs`
//! (Hecke/Newton, linear algebra): for a curve E/F_p (p = 1 mod l) whose l+1 l-isogenies are all
//! rational, every subgroup's codomain j-invariant is computed by Velu, so
//!     Phi_l(X, j(E)) = prod_{C} (X - j(E/C))
//! is read off directly; interpolating over l+1 such curves (Phi_l is monic of degree l+1 in Y)
//! recovers the bivariate Phi_l(X, Y) mod p. The result is checked (in tests) for exact equality
//! against the q-expansion Phi_l.
use crate::curve::{is_smooth, jinv, rhs, Curve};
use crate::field::{is_prime, Field, Rng, Zp};
use crate::find::divpoly::division_poly;
use crate::int::Int;
use crate::kernel::radical::{model3, model7, tate_model, to_short};
use crate::kernel::xonly::{multiples_x_fast, velu_xonly_fast, xadd_c, xdbl_c, XConst};
use crate::poly;
use std::collections::{HashSet, VecDeque};

/// Lagrange interpolation over F_p: the unique polynomial (coeffs low->high) of degree < xs.len()
/// through (xs[k], ys[k]). xs must be distinct.
fn lagrange(fp: &Zp, xs: &[u64], ys: &[u64]) -> Vec<u64> {
    let n = xs.len();
    let mut acc = vec![0u64];
    for k in 0..n {
        // basis L_k(Y) = prod_{m!=k} (Y - xs[m]) / (xs[k] - xs[m])
        let mut num = vec![fp.one()];
        let mut den = fp.one();
        for (m, &xm) in xs.iter().enumerate() {
            if m == k {
                continue;
            }
            num = poly::mul(fp, &num, &vec![fp.neg(xm), fp.one()]);
            den = fp.mul(den, fp.sub(xs[k], xm));
        }
        let scale = fp.mul(ys[k], fp.inv(den));
        let term = poly::scale(fp, &num, scale);
        acc = poly::add(fp, &acc, &term);
    }
    // pad/truncate to length n
    acc.resize(n, 0);
    acc
}

/// The l+1 cyclic subgroups C of E[l], each as the x-coordinates of (C \ O) / +-1, when all of
/// them are rational; None otherwise. With p = 1 mod l, all l+1 are rational iff Frobenius is
/// +-1 on E[l] (its determinant is p = 1), iff every root of psi_l lies in F_p, iff
/// x^p = x mod psi_l: one modular power decides, for any curve. The (l^2-1)/2 roots are then
/// grouped by x-only multiples (x-only, so the twist case needs no square root).
fn full_groups(fp: &Zp, e: &Curve<u64>, ell: usize, rng: &mut Rng) -> Option<Vec<Vec<u64>>> {
    let psi = poly::monic(fp, &division_poly(fp, e, ell));
    let x = poly::x_poly(fp);
    if poly::deg(
        fp,
        &poly::sub(fp, &poly::powmod(fp, &x, fp.p as u128, &psi), &x),
    ) >= 0
    {
        return None;
    }
    let mut xs = vec![];
    poly::split_roots(fp, &psi, rng, &mut xs);
    let roots: HashSet<u64> = xs.iter().copied().collect();
    let mut taken = HashSet::new();
    let mut groups = vec![];
    for &x0 in &xs {
        if taken.contains(&x0) {
            continue;
        }
        let g = multiples_x_fast(fp, e, x0, (ell - 1) / 2)?;
        for &xk in &g {
            // x(kP) is a root of psi_l not yet in any subgroup (the subgroups partition the roots)
            if !roots.contains(&xk) || !taken.insert(xk) {
                return None;
            }
        }
        groups.push(g);
    }
    (groups.len() == ell + 1).then_some(groups)
}

/// Velu's codomain for the odd-order kernel whose (C \ O) / +-1 has these x-coordinates:
/// v = sum (6x^2 + 2a), w = sum (4(x^3 + ax + b) + x (6x^2 + 2a)), (a - 5v, b - 7w).
fn velu_cod(fp: &Zp, e: &Curve<u64>, xs: &[u64]) -> Curve<u64> {
    let (mut v, mut w) = (0, 0);
    for &x in xs {
        let x2 = fp.sq(x);
        let vq = fp.add(fp.mul(6, x2), fp.mul(2, e.a));
        let uq = fp.mul(4, fp.add(fp.mul(fp.add(x2, e.a), x), e.b));
        v = fp.add(v, vq);
        w = fp.add(w, fp.add(uq, fp.mul(x, vq)));
    }
    Curve::new(fp.sub(e.a, fp.mul(5, v)), fp.sub(e.b, fp.mul(7, w)))
}

/// Images under a Velu isogeny with representatives (x_Q, v_Q, u_Q) of the x-coordinates xs:
/// x + sum v_Q / (x - x_Q) + u_Q / (x - x_Q)^2 (one batch inversion; x != x_Q because the
/// kernel's order is prime to the order of the points mapped).
fn push_x(fp: &Zp, reps: &[(u64, u64, u64)], xs: &[u64]) -> Vec<u64> {
    let dens: Vec<u64> = xs
        .iter()
        .flat_map(|&x| reps.iter().map(move |r| fp.sub(x, r.0)))
        .collect();
    let inv = fp.batch_inv(&dens);
    xs.iter()
        .enumerate()
        .map(|(i, &x)| {
            reps.iter().enumerate().fold(x, |acc, (k, r)| {
                let t = inv[i * reps.len() + k];
                fp.add(acc, fp.add(fp.mul(r.1, t), fp.mul(r.2, fp.sq(t))))
            })
        })
        .collect()
}

/// ([m]P, [m+1]P) by the x-only Montgomery ladder (m >= 1, projective (X, Z)).
fn ladder(fp: &Zp, c: &XConst<Zp>, p: (u64, u64), m: u64) -> ((u64, u64), (u64, u64)) {
    let (mut r0, mut r1) = (p, xdbl_c(fp, c, p));
    for i in (0..63 - m.leading_zeros()).rev() {
        if m >> i & 1 == 1 {
            r0 = xadd_c(fp, c, r1, r0, p);
            r1 = xdbl_c(fp, c, r1);
        } else {
            r1 = xadd_c(fp, c, r1, r0, p);
            r0 = xdbl_c(fp, c, r0);
        }
    }
    (r0, r1)
}

/// The multiples of l^2 in the Hasse interval of F_p, as (first multiplier, number of further
/// steps): N = l^2 (m_lo + k), k = 0..=steps.
fn hasse_multiples(p: u64, l2: u64) -> (u64, u64) {
    let s = (p as f64).sqrt() as u64 + 1;
    let (lo, hi) = (p + 1 - 2 * s, p + 1 + 2 * s);
    let m_lo = lo.div_ceil(l2);
    (m_lo, (hi / l2).saturating_sub(m_lo))
}

/// A necessary condition for qualifying that is far cheaper than x^p = x mod psi_l: E or its
/// quadratic twist then has all of E[l] rational, so l^2 divides its point count N, which lies
/// in the Hasse interval. For a random x on E and (if `twist_too`) one on the twist (x-only
/// arithmetic is shared by both), test whether [N]P = O for some multiple N of l^2 in the
/// interval: one ladder to the first multiple, then differential additions in steps of [l^2]P.
/// It only saves work: a false positive costs the exact test, a false negative a fresh curve.
/// A curve with a rational point of order l can only qualify on its own side (Frobenius is 1,
/// not -1, on that point), so it needs no twist point.
fn may_qualify(
    fp: &Zp,
    e: &Curve<u64>,
    l2: u64,
    (m_lo, steps): (u64, u64),
    twist_too: bool,
    rng: &mut Rng,
) -> bool {
    let c = XConst::new(fp, e);
    let scan = |x: u64| {
        let s = ladder(fp, &c, (x, 1), l2).0; // [l^2]P
        if s.1 == 0 {
            return true;
        }
        let (mut prev, mut cur) = ladder(fp, &c, s, m_lo - 1);
        for _ in 0..=steps {
            if cur.1 == 0 {
                return true;
            }
            let next = xadd_c(fp, &c, cur, s, prev);
            (prev, cur) = (cur, next);
        }
        false
    };
    let (mut on_e, mut on_twist) = (false, !twist_too);
    while !(on_e && on_twist) {
        let x = fp.random(rng);
        let r = rhs(fp, e, x);
        if r == 0 {
            continue;
        }
        let side = if fp.legendre(r) == 1 {
            &mut on_e
        } else {
            &mut on_twist
        };
        if *side {
            continue;
        }
        *side = true;
        if scan(x) {
            return true;
        }
    }
    false
}

/// The F_p-rational m-isogenies of E for a small prime m, as (codomain, Velu representatives
/// (x_Q, v_Q, u_Q)): for m = 2 one per rational root x0 of x^3 + a x + b (v = 3 x0^2 + a,
/// u = 0), for odd m x-only Velu from each rational root of psi_m (a rational x(P) makes the
/// kernel's x-coordinates rational).
fn small_isogenies(
    fp: &Zp,
    e: &Curve<u64>,
    m: usize,
    rng: &mut Rng,
) -> Vec<(Curve<u64>, Vec<(u64, u64, u64)>)> {
    if m == 2 {
        let cubic = vec![e.b, e.a, 0, 1];
        return poly::roots(fp, &cubic, rng)
            .into_iter()
            .map(|x0| {
                let v = fp.add(fp.mul(3, fp.sq(x0)), e.a);
                let cod = Curve::new(
                    fp.sub(e.a, fp.mul(5, v)),
                    fp.sub(e.b, fp.mul(7, fp.mul(x0, v))),
                );
                (cod, vec![(x0, v, 0)])
            })
            .collect();
    }
    poly::roots(fp, &division_poly(fp, e, m), rng)
        .into_iter()
        .filter_map(|x0| velu_xonly_fast(fp, e, x0, m as u64).map(|iso| (iso.cod, iso.reps)))
        .collect()
}

/// A random curve with a rational point of order l, for the l whose modular curve X_1(l) has a
/// one-parameter model at hand (l = 3, 5, 7; `kernel::radical`); None for other l. Frobenius
/// then fixes that point and has determinant p = 1 mod l, so it is unipotent on E[l] and is the
/// identity (all l+1 l-isogenies rational) about one time in l, against about one in
/// l (l^2 - 1) / 2 for a uniformly random curve.
fn x1_curve(fp: &Zp, ell: usize, rng: &mut Rng) -> Option<Curve<u64>> {
    let a = match ell {
        3 => model3(fp, fp.random(rng), fp.random(rng)),
        5 => {
            let b = fp.random(rng);
            tate_model(fp, b, b)
        }
        7 => model7(fp, fp.random(rng)),
        _ => return None,
    };
    Some(to_short(fp, a))
}

/// Phi_l mod p by the isogeny-graph method. Returns the (l+2) x (l+2) coefficient matrix c[i][j]
/// (coefficient of X^i Y^j) mod p, or None if fewer than l+1 curves with all l+1 l-isogenies
/// rational were found within the search budget. Requires p prime and p = 1 mod l.
///
/// Random curves qualify rarely (about 2 / |SL_2(F_l)|), but qualifying curves come in families:
/// an isogeny of degree m prime to l maps E[l] isomorphically onto E'[l] and commutes with
/// Frobenius, so every rational 2-, 3- or 5-isogenous curve (m != l) of a qualifying curve
/// qualifies too, and Velu's x-map carries its l-torsion subgroups over: no division polynomial,
/// no root finding, just the roots of a small polynomial for the m-isogenies. Those are used
/// first; then the l-codomains (qualifying when above the floor of the l-volcano); random curves
/// only to find the first one: from X_1(l) where available, otherwise uniform behind a cheap
/// necessary test. The interpolation is unique, so the order in which the
/// l+1 curves turn up does not change the result.
pub fn phi_mod_p(p: u64, ell: usize, rng: &mut Rng) -> Option<Vec<Vec<u64>>> {
    assert!(is_prime(p), "p must be prime");
    assert!(
        p % ell as u64 == 1,
        "need p = 1 mod l so the full l-torsion can be rational"
    );
    let fp = Zp::new(p);
    let l1 = ell + 1;
    // Phi_l is monic of degree l+1 in Y (c[i][l+1] = [i = 0]), so each X^i coefficient minus
    // its known Y^(l+1) term has degree <= l in Y: l + 1 distinct j0 determine it
    let need = l1;
    let mut j0s: Vec<u64> = vec![];
    let mut cols: Vec<Vec<u64>> = vec![]; // per curve: coeffs of Phi(X, j0) in X, length l+2
                                          // j-invariants already tested (the property depends only on j: twists flip the sign of
                                          // Frobenius, and j != 0, 1728 has no other twists)
    let mut tested = HashSet::new();
    // curves known to qualify (prime-to-l isogenous to one that does) with their l-torsion
    // subgroups, qualifying curves whose prime-to-l neighbours are not generated yet, and
    // l-codomains (qualifying above the floor)
    let mut sure: VecDeque<(Curve<u64>, Vec<Vec<u64>>)> = VecDeque::new();
    let mut pending: VecDeque<(Curve<u64>, Vec<Vec<u64>>)> = VecDeque::new();
    let mut walk: VecDeque<Curve<u64>> = VecDeque::new();
    let mut tries = 0u64;
    // the l^2 filter pays when the Hasse interval holds few multiples of l^2. Measured (medians
    // over 8-16 primes, exact against Hecke): uniform curves l = 11 / 13 at 17 bits 93 -> 5.8 ms
    // and 395 -> 33 ms; X_1 samples l = 5 / 7 at 17 bits 0.38 -> 0.26 ms and 0.65 -> 0.44 ms; with
    // 164-455 multiples (l = 3, 5 at 20 bits) it was 2-5x slower, hence the bound
    let l2 = (ell * ell) as u64;
    let hm = hasse_multiples(p, l2);
    let filter = hm.0 >= 2 && hm.1 <= 4 * l2;
    while j0s.len() < need && tries < 400_000 {
        let mut from_x1 = false;
        let (e, known) = if let Some((e, g)) = sure.pop_front() {
            (e, Some(g))
        } else if let Some((h, hg)) = pending.pop_front() {
            for m in [2usize, 3, 5].into_iter().filter(|&m| m != ell) {
                for (cod, reps) in small_isogenies(&fp, &h, m, rng) {
                    let groups = hg.iter().map(|g| push_x(&fp, &reps, g)).collect();
                    sure.push_back((cod, groups));
                }
            }
            continue;
        } else if let Some(e) = walk.pop_front() {
            (e, None)
        } else {
            tries += 1;
            let e = match x1_curve(&fp, ell, rng) {
                Some(e) => {
                    from_x1 = true;
                    e
                }
                None => Curve::new(fp.random(rng), fp.random(rng)),
            };
            (e, None)
        };
        if !is_smooth(&fp, &e) {
            continue;
        }
        let j0 = jinv(&fp, &e);
        if j0 == 0 || j0 == 1728 || !tested.insert(j0) {
            continue;
        }
        let groups = match known {
            Some(g) => g,
            None => {
                if filter && !may_qualify(&fp, &e, l2, hm, !from_x1, rng) {
                    continue;
                }
                let Some(g) = full_groups(&fp, &e, ell, rng) else {
                    continue;
                };
                g
            }
        };
        let cods: Vec<Curve<u64>> = groups.iter().map(|g| velu_cod(&fp, &e, g)).collect();
        // Phi(X, j0) = prod over the l+1 subgroups of (X - j(E/C))
        let mut px = vec![fp.one()];
        for cod in &cods {
            px = poly::mul(&fp, &px, &vec![fp.neg(jinv(&fp, cod)), fp.one()]);
        }
        px.resize(l1 + 1, 0); // monic degree l+1: coeffs c0..c_{l+1}
        j0s.push(j0);
        cols.push(px);
        pending.push_back((e, groups));
        walk.extend(cods);
    }
    if j0s.len() < need {
        return None;
    }
    // interpolate each X-power's coefficient as a polynomial in Y
    let mut c = vec![vec![0u64; l1 + 1]; l1 + 1];
    c[0][l1] = 1;
    for i in 0..=l1 {
        let ys: Vec<u64> = cols
            .iter()
            .zip(&j0s)
            .map(|(col, &j0)| {
                if i == 0 {
                    fp.sub(col[0], fp.pow(j0, l1 as u128))
                } else {
                    col[i]
                }
            })
            .collect();
        let py = lagrange(&fp, &j0s, &ys);
        c[i][..l1].copy_from_slice(&py);
    }
    Some(c)
}

/// Integer Phi_l by the Sutherland method: CRT of `phi_mod_p` over enough primes p = 1 mod l to
/// exceed the coefficient height bound. Coefficients are returned centered (symmetric residues).
pub fn phi_crt(ell: usize) -> Option<Vec<Vec<Int>>> {
    use crate::bigint::Big;
    let l = ell as f64;
    let nats = 6.0 * l * l.ln() + 16.0 * l + 14.0 * l.sqrt() * l.ln();
    let bits = (nats / std::f64::consts::LN_2) as usize + 24;
    // primes p = 1 mod l near 2^20 (small enough that full-structure curves are easy to find,
    // large enough that few primes are needed)
    let mut primes = vec![];
    let mut cand = (1u64 << 20) + 1;
    while primes.len() * 20 < bits {
        if cand % ell as u64 == 1 && is_prime(cand) {
            primes.push(cand);
        }
        cand += 1;
        if cand > (1u64 << 21) {
            return None;
        }
    }
    let mut mats = vec![];
    for (idx, &p) in primes.iter().enumerate() {
        let mut rng = Rng::new(0x5107 ^ p ^ (idx as u64));
        mats.push(phi_mod_p(p, ell, &mut rng)?);
    }
    let rows = mats[0].len();
    let cols = mats[0][0].len();
    let mut m = Big::from_u64(1);
    for &p in &primes {
        m = m.mul_small(p);
    }
    let mi = Int::from_big(&m);
    let out = (0..rows)
        .map(|i| {
            (0..cols)
                .map(|k| {
                    let mut x = Big::zero();
                    let mut mm = Big::from_u64(1);
                    for (t, &p) in primes.iter().enumerate() {
                        let fp = Zp::new(p);
                        let r = mats[t][i][k];
                        let delta = fp.mul(fp.sub(r, x.rem_small(p)), fp.inv(mm.rem_small(p)));
                        x = x.add(&mm.mul_small(delta));
                        mm = mm.mul_small(p);
                    }
                    let xi = Int::from_big(&x);
                    if x.add(&x) > m {
                        &xi - &mi
                    } else {
                        xi
                    }
                })
                .collect()
        })
        .collect();
    Some(out)
}
