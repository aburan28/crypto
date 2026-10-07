//! Route 2: geometrically cyclic covers of prime degree `n` of a curve `B` of
//! genus ≤ 1 over `F_2` (`B = P¹`, or an elliptic curve of trace `t`).
//!
//! The Klein quartic is a cyclic triple cover of `E_1` whose automorphism is
//! defined over `F_2`.  For `C → B` geometrically cyclic of prime degree `n`
//! with `r` branch points (all totally ramified),
//! `g(C) = n(g_B − 1) + 1 + (n − 1)r/2` and `Jac(C) ~ Jac(B) × P`.  `A` is simple
//! of dimension `n − 1 > g_B`, so `A ⊆ Jac(C)` forces `A ⊆ P`.  If Frobenius acts
//! on the group by a unit of order `d` (`d | n − 1`), two derived constraints
//! bound `r`:
//!
//! * **class field theory** — over `F_{2^d}` the group is a quotient of a ray
//!   class group; its `n`-part comes either from `Pic⁰(B)`, which needs Frobenius
//!   to have an eigenvalue of order `d` on `Jac(B)[n]` (for elliptic `B`: a root
//!   `u` of `x² − tx + 2 mod n` with `ord_n(u) = d`), or from a branch place whose
//!   residue field contains `μ_n`, so `r ≥ (n − 1)/d` (`r ≥ 2` when
//!   `d = n − 1`, since the constants then absorb one place's `n`-part);
//! * **`μ_d`-stability** — Frobenius permutes the `n − 1` non-trivial characters
//!   in orbits of length `d`, so the multiset of `P`'s Weil numbers is stable
//!   under `μ_d`; it then contains `μ_d·{ζα, ζᾱ}`, with `2d(n − 1)` distinct
//!   elements: `dim P ≥ d(n − 1)`, i.e. `r ≥ 2(d + 1 − g_B)`.

use super::weil::{curve_order, mult_order};
use num_bigint::BigInt;
use num_traits::{ToPrimitive, Zero};
use serde::Serialize;

#[derive(Clone, Copy, Debug, Serialize)]
pub enum Base {
    ProjectiveLine,
    /// An elliptic curve over `F_2` with this Frobenius trace (−2 ≤ t ≤ 2).
    Elliptic(i64),
}

impl Base {
    pub fn genus(self) -> u64 {
        match self {
            Base::ProjectiveLine => 0,
            Base::Elliptic(_) => 1,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct CyclicRow {
    pub d: u32,
    pub picard_part_available: bool,
    pub cft_min_branch_points: u32,
    pub mu_d_min_branch_points: u32,
    pub min_branch_points: u32,
    pub min_genus: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct CyclicBound {
    pub n: u32,
    pub base: Base,
    pub rows: Vec<CyclicRow>,
    pub min_genus: u64,
    pub argmin_d: u32,
    /// For an elliptic base: the family with two branch points (genus `n`) is
    /// excluded for every twist.
    pub two_branch_points_excluded: bool,
}

fn picard_available(n: u32, base: Base, d: u32) -> bool {
    match base {
        Base::ProjectiveLine => false,
        Base::Elliptic(t) => (1..n as i64).any(|u| {
            (u * u - t * u + 2).rem_euclid(n as i64) == 0
                && mult_order(u as u64, n as u64) == d as u64
        }),
    }
}

pub fn bound_for(n: u32, base: Base) -> CyclicBound {
    let gb = base.genus();
    let mut rows = Vec::new();
    for d in 1..n {
        if !(n - 1).is_multiple_of(d) {
            continue;
        }
        let pic = picard_available(n, base, d);
        let cft = if pic {
            0
        } else if d == n - 1 {
            2
        } else {
            (n - 1) / d
        };
        let mu = 2 * (d + 1 - gb as u32);
        let r = cft.max(mu).max(2);
        let genus = (n as i64 * (gb as i64 - 1) + 1 + (n as i64 - 1) * r as i64 / 2) as u64;
        rows.push(CyclicRow {
            d,
            picard_part_available: pic,
            cft_min_branch_points: cft,
            mu_d_min_branch_points: mu,
            min_branch_points: r,
            min_genus: genus,
        });
    }
    let best = rows
        .iter()
        .min_by_key(|r| r.min_genus)
        .expect("d = 1 always present")
        .clone();
    CyclicBound {
        n,
        base,
        two_branch_points_excluded: matches!(base, Base::Elliptic(_))
            && rows.iter().all(|r| r.min_branch_points > 2),
        min_genus: best.min_genus,
        argmin_d: best.d,
        rows,
    }
}

/// `E_0` as the base: the case the Klein quartic generalises.
pub fn bound(n: u32) -> CyclicBound {
    bound_for(n, Base::Elliptic(-1))
}

/// Every base of genus ≤ 1 over `F_2`: `P¹` and the five elliptic isogeny classes.
pub fn all_bases(n: u32) -> Vec<CyclicBound> {
    let mut out = vec![bound_for(n, Base::ProjectiveLine)];
    out.extend((-2..=2).map(|t| bound_for(n, Base::Elliptic(t))));
    out
}

/// `#E_0(F_{2^d}) mod n` for the record.
pub fn e0_order_mod(n: u32, d: u32) -> u64 {
    let m = curve_order(-1, d) % BigInt::from(n);
    debug_assert!(!m.is_zero() || n != 131);
    m.to_u64().expect("small")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn toy_and_challenge_bounds() {
        // n = 3: the Klein-type d = 1 cover reaches genus 3.
        assert_eq!(bound(3).min_genus, 3);
        // n = 5: genus 9 at best.
        assert_eq!(bound(5).min_genus, 9);
        // n = 131: d = 10 gives r ≥ max(20, 13) = 20, genus 1301.
        let b = bound(131);
        assert_eq!((b.min_genus, b.argmin_d), (1301, 10));
        assert!(b.two_branch_points_excluded);
        // over P¹: r ≥ max(2d + 2, 130/d) = 22 at d = 10, genus 1300.
        assert_eq!(bound_for(131, Base::ProjectiveLine).min_genus, 1300);
        assert!(all_bases(131).iter().all(|b| b.min_genus >= 1300));
    }
}
