//! The Klein quartic `x³y + y³z + z³x = 0`: the `n = 3` instance of the curve
//! the ECC2K-130 question asks for, for the *other* Koblitz sign.
//!
//! Over `F_2` its Jacobian is isogenous to `Res_{F_8/F_2}(E_1)` (checked from
//! point counts over `F_{2^k}`, `k ≤ 8`); over `Q` its point counts match the
//! Weil restriction of the conductor-49 CM curve from `Q(ζ_7)^+` at every
//! prime checked; and its `F_2`-rational automorphism `(x:y:z) ↦ (y:z:x)`
//! presents it as a cyclic triple cover of an elliptic curve branched at the
//! two points of one degree-2 place.

use super::enumerate::{
    form_value, gl3_f2, plane_count, plane_smooth, projective_points, quartic_monomials,
    substitute_quartic,
};
use super::gf2k::SmallField;
use super::weil::frobenius_trace;
use num_traits::ToPrimitive;
use serde::Serialize;
use std::collections::BTreeSet;

pub fn klein_form(monos: &[Vec<usize>]) -> u64 {
    [[3usize, 1, 0], [0, 3, 1], [1, 0, 3]]
        .iter()
        .map(|e| {
            1u64 << monos
                .iter()
                .position(|m| m.as_slice() == e)
                .expect("quartic monomial")
        })
        .fold(0, |a, b| a | b)
}

#[derive(Clone, Debug, Serialize)]
pub struct ModPCheck {
    pub primes_checked: usize,
    pub largest_prime: u64,
    pub mismatches: usize,
    /// `(p, #C(F_p), a_p of the conductor-49 curve, predicted)` for the first primes.
    pub sample: Vec<(u64, i64, i64, i64)>,
}

#[derive(Clone, Debug, Serialize)]
pub struct KleinReport {
    pub counts_over_f2k: Vec<i64>,
    pub predicted_res_f8_f2_of_e1: Vec<i64>,
    pub matches_weil_restriction: bool,
    pub smooth_checked_through_k: u32,
    pub smooth: bool,
    pub gl3_orbit_size: usize,
    pub stabiliser_order: usize,
    pub cyclic_automorphism_preserves_form: bool,
    pub cyclic_automorphism_fixed_points: usize,
    pub fixed_points_first_field_degree: u32,
    pub quotient_genus_by_riemann_hurwitz: i64,
    pub over_q: ModPCheck,
}

fn primes_up_to(n: u64) -> Vec<u64> {
    let mut sieve = vec![true; n as usize + 1];
    let mut out = Vec::new();
    for i in 2..=n as usize {
        if sieve[i] {
            out.push(i as u64);
            let mut j = i * i;
            while j <= n as usize {
                sieve[j] = false;
                j += i;
            }
        }
    }
    out
}

fn legendre(a: i64, p: u64) -> i64 {
    let a = a.rem_euclid(p as i64) as u64;
    if a == 0 {
        return 0;
    }
    let mut r = 1u64;
    let mut b = a;
    let mut e = (p - 1) / 2;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % p;
        }
        b = b * b % p;
        e >>= 1;
    }
    if r == 1 {
        1
    } else {
        -1
    }
}

/// `a_p` of `y² + xy = x³ − x² − 2x − 1` (conductor 49, CM by `Q(√−7)`), odd `p ≠ 7`:
/// completing the square gives `a_p = −Σ_x (4x³ − 3x² − 8x − 4 | p)`.
pub fn a_p_conductor49(p: u64) -> i64 {
    let pi = p as i64;
    -(0..pi)
        .map(|x| legendre(4 * x * x * x - 3 * x * x - 8 * x - 4, p))
        .sum::<i64>()
}

fn klein_count_mod_p(p: u64) -> i64 {
    // z = 1 chart: x³y + y³ + x = 0; at z = 0 the points (1:0:0), (0:1:0).
    let cubes: Vec<u64> = (0..p).map(|v| v * v % p * v % p).collect();
    let mut c = 2i64;
    for x in 0..p {
        for y in 0..p {
            if (cubes[x as usize] * y + cubes[y as usize] + x).is_multiple_of(p) {
                c += 1;
            }
        }
    }
    c
}

pub fn run(max_prime: u64) -> KleinReport {
    let monos = quartic_monomials();
    let klein = klein_form(&monos);
    let counts: Vec<i64> = (1..=8)
        .map(|k| plane_count(&SmallField::new(k), &monos, klein))
        .collect();
    let predicted: Vec<i64> = (1..=8u32)
        .map(|k| {
            let res_trace = if k % 3 == 0 {
                3 * frobenius_trace(1, k).to_i64().expect("small")
            } else {
                0
            };
            (1i64 << k) + 1 - res_trace
        })
        .collect();
    let smooth = plane_smooth(&monos, klein, 8);
    let orbit: BTreeSet<u64> = gl3_f2()
        .into_iter()
        .map(|a| substitute_quartic(&monos, klein, a))
        .collect();
    let sigma = [[0u8, 1, 0], [0, 0, 1], [1, 0, 0]];
    let invariant = substitute_quartic(&monos, klein, sigma) == klein;

    // fixed points of (x:y:z) ↦ (y:z:x) on the curve over F_{2^k}, k ≤ 6
    let mut fixed = 0usize;
    let mut first_k = 0u32;
    for k in 1..=6u32 {
        let field = SmallField::new(k);
        let mut here = 0usize;
        for pt in projective_points(&field, 2) {
            if form_value(&field, &monos, klein, &pt) != 0 {
                continue;
            }
            let img = [pt[1], pt[2], pt[0]];
            // projectively equal: all 2×2 minors vanish
            let eq = (0..3)
                .all(|i| (0..3).all(|j| field.mul(pt[i], img[j]) == field.mul(pt[j], img[i])));
            if eq {
                here += 1;
            }
        }
        if here > fixed {
            if fixed == 0 {
                first_k = k;
            }
            fixed = here;
        }
    }
    // Riemann–Hurwitz for a degree-3 cyclic quotient of a genus-3 curve:
    // 2·3 − 2 = 3(2h − 2) + 2·fixed.
    let quotient_genus = ((4 - 2 * fixed as i64) / 3 + 2) / 2;

    let mut sample = Vec::new();
    let mut mismatches = 0;
    let primes: Vec<u64> = primes_up_to(max_prime)
        .into_iter()
        .filter(|&p| p > 2 && p != 7)
        .collect();
    for &p in &primes {
        let c = klein_count_mod_p(p);
        let ap = a_p_conductor49(p);
        let pred = p as i64 + 1 - if p % 7 == 1 || p % 7 == 6 { 3 * ap } else { 0 };
        if c != pred {
            mismatches += 1;
        }
        if sample.len() < 12 {
            sample.push((p, c, ap, pred));
        }
    }

    KleinReport {
        counts_over_f2k: counts.clone(),
        predicted_res_f8_f2_of_e1: predicted.clone(),
        matches_weil_restriction: counts == predicted,
        smooth_checked_through_k: 8,
        smooth,
        gl3_orbit_size: orbit.len(),
        stabiliser_order: 168 / orbit.len(),
        cyclic_automorphism_preserves_form: invariant,
        cyclic_automorphism_fixed_points: fixed,
        fixed_points_first_field_degree: first_k,
        quotient_genus_by_riemann_hurwitz: quotient_genus,
        over_q: ModPCheck {
            primes_checked: primes.len(),
            largest_prime: *primes.last().unwrap_or(&0),
            mismatches,
            sample,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conductor49_traces() {
        // CM by Q(√−7): a_p = 0 at primes inert in Q(√−7) (p ≡ 3, 5, 6 mod 7).
        for p in [3u64, 5, 13, 17, 19, 31] {
            assert_eq!(a_p_conductor49(p), 0, "p = {p}");
        }
        assert_ne!(a_p_conductor49(11), 0);
    }

    #[test]
    fn klein_over_small_primes() {
        let r = run(60);
        assert!(r.matches_weil_restriction);
        assert_eq!(r.over_q.mismatches, 0);
        assert_eq!(r.cyclic_automorphism_fixed_points, 2);
        assert_eq!(r.quotient_genus_by_riemann_hurwitz, 1);
    }
}
