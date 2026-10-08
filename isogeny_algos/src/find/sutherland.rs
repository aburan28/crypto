//! Sutherland-style modular polynomials (Sutherland 2011, "Computing modular polynomials with the
//! Chinese remainder theorem"): compute Phi_l mod p from the l-isogeny graph, then CRT to the
//! integer polynomial. The mod-p step is independent of the q-expansion methods in `modpoly.rs`
//! (Hecke/Newton, linear algebra): for a curve E/F_p (p = 1 mod l) whose l+1 l-isogenies are all
//! rational, every subgroup's codomain j-invariant is computed by Velu/Kohel, so
//!     Phi_l(X, j(E)) = prod_{C} (X - j(E/C))
//! is read off directly; interpolating over l+2 such curves recovers the bivariate Phi_l(X, Y)
//! mod p. The result is checked (in tests) for exact equality against the q-expansion Phi_l.
use crate::curve::{is_smooth, jinv, Curve};
use crate::field::{is_prime, Field, Rng, Zp};
use crate::find::divpoly::kernel_polys;
use crate::int::Int;
use crate::kernel::kohel::kohel;
use crate::poly;
use std::collections::HashSet;

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

/// Phi_l mod p by the isogeny-graph method. Returns the (l+2) x (l+2) coefficient matrix c[i][j]
/// (coefficient of X^i Y^j) mod p, or None if fewer than l+2 curves with all l+1 l-isogenies
/// rational were found within the search budget. Requires p prime and p = 1 mod l.
pub fn phi_mod_p(p: u64, ell: usize, rng: &mut Rng) -> Option<Vec<Vec<u64>>> {
    assert!(is_prime(p), "p must be prime");
    assert!(p % ell as u64 == 1, "need p = 1 mod l so the full l-torsion can be rational");
    let fp = Zp::new(p);
    let l1 = ell + 1;
    let need = l1 + 1; // l + 2 distinct j0 for a degree <= l+1 interpolation in Y
    let mut j0s: Vec<u64> = vec![];
    let mut cols: Vec<Vec<u64>> = vec![]; // per curve: coeffs of Phi(X, j0) in X, length l+2
    let mut seen = HashSet::new();
    let mut tries = 0u64;
    while j0s.len() < need && tries < 400_000 {
        tries += 1;
        let e = Curve::new(fp.random(rng), fp.random(rng));
        if !is_smooth(&fp, &e) {
            continue;
        }
        let j0 = jinv(&fp, &e);
        if j0 == 0 || j0 == 1728 || seen.contains(&j0) {
            continue;
        }
        let kps = kernel_polys(&fp, &e, ell as u64, rng);
        if kps.len() != l1 {
            continue; // need every l+1 l-isogeny rational (pi scalar on E[l])
        }
        // Phi(X, j0) = prod over the l+1 subgroups of (X - j(E/C))
        let mut px = vec![fp.one()];
        for h in &kps {
            let iso = kohel(&fp, &e, h, ell as u64);
            let jc = jinv(&fp, &iso.cod);
            px = poly::mul(&fp, &px, &vec![fp.neg(jc), fp.one()]);
        }
        px.resize(l1 + 1, 0); // monic degree l+1: coeffs c0..c_{l+1}
        seen.insert(j0);
        j0s.push(j0);
        cols.push(px);
    }
    if j0s.len() < need {
        return None;
    }
    // interpolate each X-power's coefficient as a polynomial in Y
    let mut c = vec![vec![0u64; l1 + 1]; l1 + 1];
    for i in 0..=l1 {
        let ys: Vec<u64> = cols.iter().map(|col| col[i]).collect();
        let py = lagrange(&fp, &j0s, &ys);
        for (j, &v) in py.iter().enumerate() {
            if j <= l1 {
                c[i][j] = v;
            }
        }
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
