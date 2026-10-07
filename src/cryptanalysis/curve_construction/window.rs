//! Index calculus over `F_2` priced from a zeta function — the native port of
//! the legacy window model in `scripts/ecc2k130_hyperelliptic_cover_boundary.py`.
//!
//! From point counts `#C(F_{2^k})`: place counts `N_d` by Möbius inversion; the
//! number of `b`-smooth effective divisors of degree `g` is the `T^g`
//! coefficient of `∏_{d≤b} (1 − T^d)^{−N_d}`, and all effective divisors of
//! degree `g` the same product over every `d ≤ g`.  One operation is one
//! divisor-class step plus one smoothness test; linear algebra is charged
//! `|FB|²·g` multiplications modulo `r`.  The figure is a **model** (a stage
//! diagnostic of a hypothetical curve), not a measurement.

use super::weil::log2_big;
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Signed, ToPrimitive, Zero};
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct Cell {
    pub genus: usize,
    pub model: String,
    pub smoothness_bound_b: usize,
    pub log2_factor_base: f64,
    pub log2_smooth_probability: f64,
    pub log2_relations: f64,
    pub log2_linear_algebra: f64,
    pub log2_total: f64,
    pub log2_ratio_to_rho: f64,
    pub beats_rho: bool,
}

fn mobius(n: usize) -> i64 {
    let (mut m, mut res, mut p) = (n, 1i64, 2usize);
    while p * p <= m {
        if m % p == 0 {
            m /= p;
            if m % p == 0 {
                return 0;
            }
            res = -res;
        }
        p += 1;
    }
    if m > 1 {
        -res
    } else {
        res
    }
}

/// `N_d` for `d = 1..=dmax` from `points[k] = #C(F_{2^k})` (index 0 unused).
pub fn place_counts(points: &[BigInt], dmax: usize) -> Vec<BigUint> {
    let mut out = vec![BigUint::zero(); dmax + 1];
    for d in 1..=dmax {
        let mut acc = BigInt::zero();
        for e in 1..=d {
            if d % e == 0 {
                acc += &points[e] * mobius(d / e);
            }
        }
        let (q, r) = (&acc / d as i64, &acc % d as i64);
        assert!(
            r.is_zero() && !q.is_negative(),
            "place count N_{d} must be a non-negative integer"
        );
        out[d] = q.to_biguint().expect("non-negative");
    }
    out
}

/// `#C(F_{2^k}) = 2^k + 1`: the random-polynomial model behind Enge–Gaudry.
pub fn rational_model_points(g: usize) -> Vec<BigInt> {
    let mut v = vec![BigInt::zero()];
    v.extend((1..=g).map(|k| (BigInt::one() << k) + 1));
    v
}

/// Price index calculus on a genus-`g` curve with the given point counts,
/// minimised over smoothness bounds `2 ≤ b ≤ bmax`.
pub fn price(points: &[BigInt], g: usize, bmax: usize, model: &str, log2_rho: f64) -> Cell {
    let places = place_counts(points, g);
    let mut ser = vec![BigUint::zero(); g + 1];
    ser[0] = BigUint::one();
    let mut smooth_at = vec![BigUint::zero(); bmax.min(g) + 1];
    let mut fb_at = vec![BigUint::zero(); bmax.min(g) + 1];
    let mut fb = BigUint::zero();
    for d in 1..=g {
        let nd = &places[d];
        if !nd.is_zero() {
            // binom(nd + j − 1, j), j = 0..=g/d
            let jmax = g / d;
            let mut binom = Vec::with_capacity(jmax + 1);
            let mut b = BigUint::one();
            binom.push(b.clone());
            for j in 1..=jmax {
                b = b * (nd + BigUint::from(j - 1)) / BigUint::from(j);
                binom.push(b.clone());
            }
            for t in (d..=g).rev() {
                let mut add = BigUint::zero();
                let mut j = 1;
                while j * d <= t {
                    if !ser[t - j * d].is_zero() {
                        add += &ser[t - j * d] * &binom[j];
                    }
                    j += 1;
                }
                ser[t] += add;
            }
        }
        if d <= bmax.min(g) {
            fb += nd;
            smooth_at[d] = ser[g].clone();
            fb_at[d] = fb.clone();
        }
    }
    let log2_total_divisors = log2_big(&ser[g]);
    let mut best: Option<Cell> = None;
    for b in 2..=bmax.min(g) {
        let fbv = &fb_at[b];
        if fbv < &BigUint::from(2u32) || smooth_at[b].is_zero() {
            continue;
        }
        let lfb = log2_big(fbv);
        let lprob = log2_big(&smooth_at[b]) - log2_total_divisors;
        let rel = lfb - lprob;
        let lin = 2.0 * lfb + (g as f64).log2();
        let tot = rel.max(lin) + 1.0;
        if best.as_ref().is_none_or(|c| tot < c.log2_total) {
            best = Some(Cell {
                genus: g,
                model: model.to_string(),
                smoothness_bound_b: b,
                log2_factor_base: lfb,
                log2_smooth_probability: lprob,
                log2_relations: rel,
                log2_linear_algebra: lin,
                log2_total: tot,
                log2_ratio_to_rho: tot - log2_rho,
                beats_rho: tot < log2_rho,
            });
        }
    }
    best.expect("some smoothness bound admits relations")
}

/// `L_{2^g}(1/2, √2)` in bits — the Enge–Gaudry heuristic, used only beyond
/// the genus where [`price`] can be evaluated; an extrapolation.
pub fn l_half_bits(g: f64) -> f64 {
    let ln_n = g * std::f64::consts::LN_2;
    std::f64::consts::SQRT_2 * (ln_n * ln_n.ln()).sqrt() / std::f64::consts::LN_2
}

pub fn to_f64_bits(x: &BigUint) -> f64 {
    x.to_f64().map(|v| v.log2()).unwrap_or_else(|| log2_big(x))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rational_model_divisor_counts() {
        // P¹ over F_2 has 2^{n+1} − 1 effective divisors of degree ≤ ... and
        // exactly 2^n monic polynomials of degree n: N_d = #irreducibles.
        let pts = rational_model_points(6);
        let places = place_counts(&pts, 6);
        let irr: Vec<u64> = places.iter().skip(1).map(|x| x.to_u64().unwrap()).collect();
        assert_eq!(irr, vec![3, 1, 2, 3, 6, 9]);
    }
}
