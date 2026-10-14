//! The composite-degree counterfactual of
//! `research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md` §7, priced with
//! the Maurer–Menezes–Teske (2002) magic number.
//!
//! For a field degree `N` and each factorisation `N = n·l`, the GHS descent to
//! `F_{2^l}` reaches genus `2^(m−1)` (type II) or `2^(m−1) − 1` (type I) for the
//! `(m, type)` pairs [`mmt_attainable`] lists.  A transfer can be injective only
//! if the target group is large enough, `l·genus ≥ N − 2` (the note's
//! admissibility rule, kept verbatim); [`counterfactual`] takes the least
//! admissible genus and prices it with the note's models:
//!
//! * `l = 1`, `m ≤ 12`: [`window::price`] on the random-polynomial point counts
//!   `#C(F_{2^k}) = 2^k + 1` (boundary E of the note);
//! * `l = 1`, larger `m`: the description-length floor, `m − 1` bits;
//! * `l ≥ 2`, `m ≤ 24`: the fixed-genus model `g!·q^(2−2/g)` (pessimistic, not
//!   tight, once `g` is large against `log q`);
//! * `l ≥ 2`, larger `m`: the floor `l·(m − 1)` bits.
//!
//! [`Rule::Legacy`] reproduces the note's table, which used the plain orbit
//! dimension and genus `2^(m−1)` (both superseded, see
//! `docs/audits/GHS_MAGIC_NUMBER_MMT_CORRECTION_20261013.md`).  Every cost is a
//! model figure for a hypothetical curve, not a measurement.

use super::weil::log2_big;
use super::window::{price, rational_model_points};
use crate::cryptanalysis::isogeny_degree_search::ball::cyclotomic_coset_sizes;
use num_bigint::BigUint;
use num_traits::One;
use serde::Serialize;
use std::collections::BTreeSet;

/// Which magic-number rule feeds the table.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Rule {
    /// MMT eq. (5): `dim span{(1, σ^i(√b))}`; genus `2^(m−1)` or `2^(m−1) − 1`.
    Mmt,
    /// Plain orbit dimension with genus `2^(m−1)` (the note's pre-2026-10-13 table).
    Legacy,
}

/// One attainable GHS outcome over `F_{2^N}/F_{2^l}`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Outcome {
    pub m: u32,
    pub type_i: bool,
}

impl Outcome {
    pub fn genus(&self) -> BigUint {
        let g = BigUint::one() << (self.m - 1);
        if self.type_i {
            g - BigUint::one()
        } else {
            g
        }
    }
}

/// Every `(m, type)` attained by some `b ≠ 0` with `n = N/l` the relative degree.
///
/// `Ord_b = (x − 1)^e·D` divides `x^n − 1 = (x − 1)^{2^v}·∏ f_i^{2^v}`
/// (`n = 2^v·n'`, `n'` odd), so `deg D` runs over the subset sums of the
/// non-`(x − 1)` factor degrees (each repeated `2^v` times) and `0 ≤ e ≤ 2^v`.
/// `e ≥ 1` gives `m = deg D + e`, type II; `e = 0` gives `m = deg D + 1`,
/// type I.  `(deg D, e) = (0, 0)` is `b = 0`.
pub fn mmt_attainable(n: u32) -> Vec<Outcome> {
    let mut odd = n;
    let mut mult = 1u32;
    while odd.is_multiple_of(2) {
        odd /= 2;
        mult *= 2;
    }
    let sizes: Vec<u32> = cyclotomic_coset_sizes(n)
        .into_iter()
        .filter(|s| *s != 1)
        .collect();
    let mut reach = vec![false; n as usize + 1];
    reach[0] = true;
    for s in sizes {
        for d in (0..=n as usize).rev() {
            if reach[d] && d + s as usize <= n as usize {
                reach[d + s as usize] = true;
            }
        }
    }
    let mut out = BTreeSet::new();
    for d in (0..=n).filter(|d| reach[*d as usize]) {
        for e in 0..=mult {
            match (d, e) {
                (0, 0) => {}
                (_, 0) => {
                    out.insert((d + 1, true));
                }
                _ => {
                    out.insert((d + e, false));
                }
            }
        }
    }
    out.into_iter()
        .map(|(m, type_i)| Outcome { m, type_i })
        .collect()
}

/// The legacy outcomes, as the note's script computed them: the dimensions of
/// the Frobenius-stable subspaces of `F_q^n` with `q = 2^l`, i.e. subset sums of
/// the factor degrees of `x^n − 1` over **`F_q`**, genus `2^(m−1)`.
///
/// That is two errors at once.  The magic number is an `F_2`-dimension, so the
/// factorisation must be over `F_2` whatever `l` is
/// (`attainable_pairs_match_an_exhaustive_scan_for_l_above_1` scans `F_{2^12}`
/// at `l = 2, 3, 4, 6` and finds the `F_2` set); and the `+1` and the type were
/// missing.  At `l = 1` the first error vanishes.
fn legacy_attainable(n: u32, l: u32) -> Vec<Outcome> {
    let mut odd = n;
    let mut mult = 1u32;
    while odd.is_multiple_of(2) {
        odd /= 2;
        mult *= 2;
    }
    let q = (0..l).fold(1u64, |acc, _| acc * 2 % u64::from(odd));
    let mut seen = vec![false; odd as usize];
    let mut sizes = Vec::new();
    for s in 0..odd {
        if seen[s as usize] {
            continue;
        }
        let (mut c, mut size) = (s, 0u32);
        loop {
            seen[c as usize] = true;
            size += 1;
            c = (u64::from(c) * q % u64::from(odd)) as u32;
            if c == s {
                break;
            }
        }
        sizes.extend(std::iter::repeat_n(size, mult as usize));
    }
    let mut reach = vec![false; n as usize + 1];
    reach[0] = true;
    for sz in sizes {
        for d in (0..=n as usize).rev() {
            if reach[d] && d + sz as usize <= n as usize {
                reach[d + sz as usize] = true;
            }
        }
    }
    (1..=n)
        .filter(|m| reach[*m as usize])
        .map(|m| Outcome { m, type_i: false })
        .collect()
}

/// `log2` of Pollard rho on a group of `2^(N−2)` elements:
/// `sqrt(π·2^(N−2) / (4N))`, the note's unit.
pub fn log2_rho(n_total: u32) -> f64 {
    0.5 * (std::f64::consts::PI.log2() + f64::from(n_total)
        - 2.0
        - (4.0 * f64::from(n_total)).log2())
}

fn log2_factorial(g: u64) -> f64 {
    (2..=g).map(|i| (i as f64).log2()).sum()
}

/// One row of the counterfactual table.
#[derive(Clone, Debug, Serialize)]
pub struct Row {
    pub field_degree_n_total: u32,
    pub subfield_degree_l: u32,
    pub extension_degree_n: u32,
    pub rule: Rule,
    /// Smallest attainable magic numbers, for the table's "available" column.
    pub available_magic_numbers: Vec<u32>,
    /// `None` when no attainable genus reaches `l·genus ≥ N − 2`.
    pub least_admissible_magic_number: Option<u32>,
    pub type_i: Option<bool>,
    pub genus: Option<String>,
    pub log2_genus: Option<f64>,
    pub cost_model: Option<String>,
    pub log2_cost: Option<f64>,
    pub smoothness_bound_b: Option<usize>,
    pub log2_rho_on_that_group: f64,
    pub breaks_rho: Option<bool>,
}

/// The table for `N`, every factorisation `N = n·l` with `n ≥ 2`.
pub fn counterfactual(n_total: u32, rule: Rule) -> Vec<Row> {
    let rho = log2_rho(n_total);
    let need = BigUint::from(n_total - 2);
    let mut rows = Vec::new();
    for l in 1..=n_total {
        if !n_total.is_multiple_of(l) || n_total / l < 2 {
            continue;
        }
        let n = n_total / l;
        let outcomes = match rule {
            Rule::Mmt => mmt_attainable(n),
            Rule::Legacy => legacy_attainable(n, l),
        };
        let mut available: Vec<u32> = outcomes.iter().map(|o| o.m).collect();
        available.dedup();
        available.truncate(12);
        let best = outcomes
            .iter()
            .filter(|o| BigUint::from(l) * o.genus() >= need)
            .min_by(|x, y| x.genus().cmp(&y.genus()));
        let mut row = Row {
            field_degree_n_total: n_total,
            subfield_degree_l: l,
            extension_degree_n: n,
            rule,
            available_magic_numbers: available,
            least_admissible_magic_number: best.map(|o| o.m),
            type_i: best.map(|o| o.type_i),
            genus: best.map(|o| o.genus().to_string()),
            log2_genus: best.map(|o| log2_big(&o.genus())),
            cost_model: None,
            log2_cost: None,
            smoothness_bound_b: None,
            log2_rho_on_that_group: rho,
            breaks_rho: None,
        };
        if let Some(o) = best {
            let g = o.genus();
            let m = o.m;
            if l == 1 && m <= 12 {
                let gg = usize::try_from(g.iter_u64_digits().next().unwrap_or(0)).unwrap();
                let cell = price(
                    &rational_model_points(gg),
                    gg,
                    200,
                    "random-polynomial",
                    rho,
                );
                row.cost_model = Some("large genus over F_2 (boundary E)".into());
                row.log2_cost = Some(cell.log2_total);
                row.smoothness_bound_b = Some(cell.smoothness_bound_b);
            } else if l == 1 {
                row.cost_model = Some("large genus over F_2, floor = bits per divisor".into());
                row.log2_cost = Some(f64::from(m - 1));
            } else if m <= 24 {
                let gg = g.iter_u64_digits().next().unwrap_or(0);
                row.cost_model = Some("fixed genus over F_2^l: g! q^(2-2/g)".into());
                row.log2_cost = Some(log2_factorial(gg) + f64::from(l) * (2.0 - 2.0 / gg as f64));
            } else {
                row.cost_model = Some("fixed genus over F_2^l: floor l(m-1)".into());
                row.log2_cost = Some(f64::from(l) * f64::from(m - 1));
            }
            row.breaks_rho = row.log2_cost.map(|c| c < rho);
        }
        rows.push(row);
    }
    rows
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_ecc::{F2mElement, IrreduciblePoly};
    use crate::cryptanalysis::ec_trapdoor::{
        ghs_genus_with_type, ghs_type_i, magic_number, mmt_condition_8, FieldTower,
    };
    use crate::cryptanalysis::ghs_screen::{screen_curve, GhsCurveInput};

    fn row(rows: &[Row], l: u32) -> &Row {
        rows.iter().find(|r| r.subfield_degree_l == l).unwrap()
    }

    /// The pre-fix rule reproduces the note's §7 table (`l = 1`: genus 128 at
    /// `2^36.41`; the dagger rows), so the port prices as the note did.
    #[test]
    fn legacy_rule_reproduces_the_note_table() {
        let rows = counterfactual(130, Rule::Legacy);
        let r1 = row(&rows, 1);
        assert_eq!(
            (r1.least_admissible_magic_number, r1.genus.as_deref()),
            (Some(8), Some("128"))
        );
        assert!(
            (r1.log2_cost.unwrap() - 36.41).abs() < 0.01,
            "{:?}",
            r1.log2_cost
        );
        assert!((r1.log2_rho_on_that_group - 60.31).abs() < 0.005);
        assert_eq!(r1.breaks_rho, Some(true));
        for (l, m, cost) in [(2, 7, 299.93), (13, 5, 68.63), (26, 4, 60.80)] {
            let r = row(&rows, l);
            assert_eq!(r.least_admissible_magic_number, Some(m), "l = {l}");
            assert!(
                (r.log2_cost.unwrap() - cost).abs() < 0.01,
                "l = {l}: {:?}",
                r.log2_cost
            );
        }
        let r131 = counterfactual(131, Rule::Legacy);
        assert_eq!(r131[0].least_admissible_magic_number, Some(130));
        assert_eq!(r131[0].log2_cost, Some(129.0));
    }

    /// `mmt_attainable` equals the `(m, type)` pairs an exhaustive scan of
    /// `ec_trapdoor` finds over `F_{2^N}/F_2`.
    #[test]
    fn attainable_pairs_match_an_exhaustive_scan() {
        for (n, low) in [(10u32, vec![0u32, 3]), (12, vec![0, 3])] {
            let irr = IrreduciblePoly {
                degree: n,
                low_terms: low,
            };
            let tower = FieldTower::new(n, n, 1, irr);
            let zero = F2mElement::zero(n);
            let mut seen = BTreeSet::new();
            for v in 1u64..(1 << n) {
                let b = F2mElement::from_biguint(&BigUint::from(v), n);
                seen.insert((magic_number(&tower, &zero, &b), ghs_type_i(&tower, &b)));
            }
            let want: Vec<(u32, bool)> =
                mmt_attainable(n).iter().map(|o| (o.m, o.type_i)).collect();
            assert_eq!(seen.into_iter().collect::<Vec<_>>(), want, "n = {n}");
        }
    }

    /// The attainable set does not depend on `l`: it is the `F_2` set for `n`
    /// at every tower of `F_{2^12}`.  (The note's script used the `F_{2^l}`
    /// factorisation of `x^n − 1`, which would allow `m = 2` at `n = 3`,
    /// `l = 4`; no `b` has it.)
    #[test]
    fn attainable_pairs_match_an_exhaustive_scan_for_l_above_1() {
        let irr = IrreduciblePoly {
            degree: 12,
            low_terms: vec![0, 3],
        };
        let zero = F2mElement::zero(12);
        for (n, l) in [(6u32, 2u32), (4, 3), (3, 4), (2, 6)] {
            let tower = FieldTower::new(12, n, l, irr.clone());
            let mut seen = BTreeSet::new();
            for v in 1u64..4096 {
                let b = F2mElement::from_biguint(&BigUint::from(v), 12);
                seen.insert((magic_number(&tower, &zero, &b), ghs_type_i(&tower, &b)));
            }
            let want: Vec<(u32, bool)> =
                mmt_attainable(n).iter().map(|o| (o.m, o.type_i)).collect();
            assert_eq!(
                seen.into_iter().collect::<Vec<_>>(),
                want,
                "n = {n}, l = {l}"
            );
        }
        // the F_q factorisation would have admitted m = 2 for n = 3, l = 4
        assert!(legacy_attainable(3, 4).iter().any(|o| o.m == 2));
        assert!(!mmt_attainable(3).iter().any(|o| o.m == 2));
    }

    /// Over `F_{2^130}/F_2` no `b` has `m = 8`; the least admissible outcome is
    /// `m = 9`, type I, genus 255, and that curve is priced below rho in the
    /// note's model.  The figure itself is pinned from the run recorded in
    /// `research/ghs_counterfactual_n130_20261013/`.
    #[test]
    fn n130_l1_least_admissible_is_genus_255() {
        let out = mmt_attainable(130);
        assert!(!out.iter().any(|o| o.m == 8));
        assert!(out.contains(&Outcome { m: 9, type_i: true }));
        let rows = counterfactual(130, Rule::Mmt);
        let r = row(&rows, 1);
        assert_eq!(r.least_admissible_magic_number, Some(9));
        assert_eq!((r.type_i, r.genus.as_deref()), (Some(true), Some("255")));
        let cost = r.log2_cost.unwrap();
        assert!((cost - EXPECTED_LOG2_COST_G255).abs() < 0.01, "cost {cost}");
        assert_eq!(r.breaks_rho, Some(cost < r.log2_rho_on_that_group));
    }

    /// Model value from the recorded run; see the test above.
    const EXPECTED_LOG2_COST_G255: f64 = 55.60;

    /// `h = (x^n − 1)/f^2` as a bit vector (index = degree), `f = Φ_5`.
    fn gf2_div_exact(mut num: Vec<u8>, den: &[u8]) -> Vec<u8> {
        let dd = den.len() - 1;
        let mut q = vec![0u8; num.len() - dd];
        for i in (dd..num.len()).rev() {
            if num[i] == 1 {
                q[i - dd] = 1;
                for (j, c) in den.iter().enumerate() {
                    num[i - dd + j] ^= c;
                }
            }
        }
        assert!(num.iter().all(|c| *c == 0), "not exact");
        q
    }

    /// A concrete `b ∈ F_{2^130}` with `Ord_b = Φ_5(x)^2`: `(x − 1) ∤ Ord_b`,
    /// `deg Ord_b = 8`, so `m = 9`, type I, genus 255, and condition (8) holds
    /// for `a = 1` (`Tr(1) = 0` at even `N`).
    #[test]
    fn genus_255_witness_exists_over_f2_130() {
        let big_n = 130u32;
        // the first irreducible trinomial z^130 + z^k + 1
        let (k, irr) = (1u32..65)
            .find_map(|k| {
                let modulus =
                    (BigUint::one() << 130usize) | (BigUint::one() << k as usize) | BigUint::one();
                let input = GhsCurveInput {
                    absolute_degree: 130,
                    modulus,
                    a: BigUint::one(),
                    b: BigUint::one(),
                };
                screen_curve(&input, &BigUint::from(4u32)).ok().map(|_| {
                    (
                        k,
                        IrreduciblePoly {
                            degree: 130,
                            low_terms: vec![0, k],
                        },
                    )
                })
            })
            .expect("an irreducible trinomial of degree 130");
        let _ = k;
        let tower = FieldTower::new(big_n, big_n, 1, irr.clone());
        let phi5 = [1u8, 1, 1, 1, 1];
        let mut den = vec![0u8; 9];
        for (i, a) in phi5.iter().enumerate() {
            for (j, b) in phi5.iter().enumerate() {
                den[i + j] ^= a & b;
            }
        }
        let mut num = vec![0u8; 131];
        num[0] = 1;
        num[130] = 1;
        let h = gf2_div_exact(num, &den);
        // s = h(σ) t, so Φ_5(σ)^2 s = (σ^130 − 1) t = 0
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut t = F2mElement::zero(big_n);
        for i in 0..big_n {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            if (state >> 33) & 1 == 1 {
                t = t.add(&F2mElement::from_bit_positions(&[i], big_n));
            }
        }
        let mut s = F2mElement::zero(big_n);
        for (i, c) in h.iter().enumerate() {
            if *c == 1 {
                s = s.add(&tower.frobenius_iter(&t, i as u32));
            }
        }
        let b = s.square(&irr);
        let a = F2mElement::one(big_n);
        assert_eq!(magic_number(&tower, &a, &b), 9);
        assert_eq!(
            ghs_genus_with_type(&tower, &b),
            (BigUint::from(255u32), true)
        );
        assert!(mmt_condition_8(&tower, &a, &b));
    }
}
