//! Sutherland-style modular polynomials (Sutherland 2011, "Computing modular polynomials with the
//! Chinese remainder theorem"): compute Phi_l mod p from the l-isogeny graph, then CRT to the
//! integer polynomial. The mod-p step is independent of the q-expansion methods in `modpoly.rs`
//! (Hecke/Newton, linear algebra): for a curve E/F_p (p = 1 mod l) whose l+1 l-isogenies are all
//! rational, every subgroup's codomain j-invariant is computed by Velu, so
//!     Phi_l(X, j(E)) = prod_{C} (X - j(E/C))
//! is read off directly; interpolating over l+1 such curves (Phi_l is monic of degree l+1 in Y)
//! recovers the bivariate Phi_l(X, Y) mod p. The result is checked (in tests) for exact equality
//! against the q-expansion Phi_l.
use crate::curve::{is_smooth, jinv, Curve};
use crate::field::{is_prime, Field, Rng, Zp};
use crate::find::divpoly::division_poly;
use crate::int::Int;
use crate::kernel::xonly::velu_xonly_fast;
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

/// The l+1 codomains of the l-isogenies of E, when all of them are rational; None otherwise.
/// With p = 1 mod l, all l+1 are rational iff Frobenius is +-1 on E[l] (its determinant is
/// p = 1), iff every root of psi_l lies in F_p, iff x^p = x mod psi_l: one modular power decides,
/// for any curve. The (l^2-1)/2 roots then split into the l+1 cyclic subgroups, each read off by
/// x-only Velu from one of its x-coordinates (x-only, so the twist case needs no square root).
fn full_codomains(fp: &Zp, e: &Curve<u64>, ell: usize, rng: &mut Rng) -> Option<Vec<Curve<u64>>> {
    let psi = poly::monic(fp, &division_poly(fp, e, ell));
    let x = poly::x_poly(fp);
    if poly::deg(fp, &poly::sub(fp, &poly::powmod(fp, &x, fp.p as u128, &psi), &x)) >= 0 {
        return None;
    }
    let mut xs = vec![];
    poly::split_roots(fp, &psi, rng, &mut xs);
    let roots: HashSet<u64> = xs.iter().copied().collect();
    let mut taken = HashSet::new();
    let mut cods = vec![];
    for &x0 in &xs {
        if taken.contains(&x0) {
            continue;
        }
        let iso = velu_xonly_fast(fp, e, x0, ell as u64)?;
        for r in &iso.reps {
            // x(kP) is a root of psi_l not yet in any subgroup (the subgroups partition the roots)
            if !roots.contains(&r.0) || !taken.insert(r.0) {
                return None;
            }
        }
        cods.push(iso.cod);
    }
    (cods.len() == ell + 1).then_some(cods)
}

/// Phi_l mod p by the isogeny-graph method. Returns the (l+2) x (l+2) coefficient matrix c[i][j]
/// (coefficient of X^i Y^j) mod p, or None if fewer than l+1 curves with all l+1 l-isogenies
/// rational were found within the search budget. Requires p prime and p = 1 mod l.
///
/// Such curves cluster: they are the vertices of the l-isogeny volcano above its floor, so the
/// codomains of each one found are tried next (a walk through the volcano), and random curves
/// are drawn only when that runs dry. The interpolation is unique, so the order in which the
/// l+1 curves turn up does not change the result.
pub fn phi_mod_p(p: u64, ell: usize, rng: &mut Rng) -> Option<Vec<Vec<u64>>> {
    assert!(is_prime(p), "p must be prime");
    assert!(p % ell as u64 == 1, "need p = 1 mod l so the full l-torsion can be rational");
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
    let mut walk: VecDeque<Curve<u64>> = VecDeque::new();
    let mut tries = 0u64;
    while j0s.len() < need && tries < 400_000 {
        let e = match walk.pop_front() {
            Some(e) => e,
            None => {
                tries += 1;
                Curve::new(fp.random(rng), fp.random(rng))
            }
        };
        if !is_smooth(&fp, &e) {
            continue;
        }
        let j0 = jinv(&fp, &e);
        if j0 == 0 || j0 == 1728 || !tested.insert(j0) {
            continue;
        }
        let Some(cods) = full_codomains(&fp, &e, ell, rng) else { continue };
        // Phi(X, j0) = prod over the l+1 subgroups of (X - j(E/C))
        let mut px = vec![fp.one()];
        for cod in &cods {
            px = poly::mul(&fp, &px, &vec![fp.neg(jinv(&fp, cod)), fp.one()]);
        }
        px.resize(l1 + 1, 0); // monic degree l+1: coeffs c0..c_{l+1}
        j0s.push(j0);
        cols.push(px);
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
            .map(|(col, &j0)| if i == 0 { fp.sub(col[0], fp.pow(j0, l1 as u128)) } else { col[i] })
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
