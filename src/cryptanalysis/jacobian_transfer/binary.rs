//! Binary-field routes: classical GHS Weil descent for every factorisation of
//! the field degree, and the Galbraith–Hess–Smart isogeny extension priced as
//! a heuristic.
//!
//! The magic number and genus come from [`crate::cryptanalysis::ec_trapdoor`]
//! and are exact invariants.  The isogeny-extension rows are heuristic: they
//! count how many curves over the same field have a small magic number for a
//! given tower, estimate how many of those lie in the target's isogeny class,
//! and price the walk needed to reach one.

use super::cost::{log2_add, CostModel};
use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use crate::cryptanalysis::ec_trapdoor::{audit_curve, DescentRow};
use num_bigint::BigUint;
use num_traits::Zero;

/// The GHS descent data for one factorisation `N = n * l`, as the
/// repository's `ec_trapdoor::audit_curve` reports it.  Its magic number
/// adds the Frobenius differences of `a` as generators without reducing `a`
/// modulo `{α² + α}`, so it can exceed the GHS value; the screen records it
/// beside [`ghs_magic`] and prices with the latter.
pub fn ghs_rows(
    n_total: u32,
    irr: &IrreduciblePoly,
    a: &F2mElement,
    b: &F2mElement,
) -> Vec<DescentRow> {
    audit_curve(n_total, irr, a, b)
}

/// The GHS magic number for the tower `F_2 <= F_(2^l) <= F_(2^N)`, `n = N/l`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GhsMagic {
    /// `dim_F2 Span{(1, σ^i(√b)) : 0 <= i < n}`, `σ` the `2^l`-Frobenius.
    pub magic_b: u32,
    /// `Tr_(F_(2^N)/F_2)(a)`.  `a` matters only modulo `{α² + α}`, i.e.
    /// through this trace: `E_(a,b) ≅ E_(a',b)` when the traces agree.
    pub trace_a: u32,
    /// Whether the reduced `a` may add one generator: its trace is 1 and
    /// `n` is even, so `a` cannot be taken in `F_2`.
    pub a_may_add_one: bool,
    /// `magic_b`.  When `a_may_add_one`, the GHS value is `magic_b` or
    /// `magic_b + 1` and the screen prices both.
    pub magic: u32,
    /// The `σ`-orbit length of `√b`; type I (genus `2^(m-1) - 1`) when it equals `magic`.
    pub orbit_length: u32,
}

fn f2_rank(mut rows: Vec<BigUint>) -> u32 {
    let mut rank = 0u32;
    let mut i = 0;
    while i < rows.len() {
        rows.sort_unstable_by(|x, y| y.cmp(x));
        rows.truncate(rows.iter().take_while(|r| !r.is_zero()).count());
        if i >= rows.len() {
            break;
        }
        let pivot = rows[i].clone();
        let top = pivot.bits() - 1;
        for r in rows.iter_mut().skip(i + 1) {
            if r.bit(top) {
                *r ^= &pivot;
            }
        }
        rank += 1;
        i += 1;
    }
    rank
}

/// The GHS magic number of `E_(a,b)` over `F_(2^N)` for base degree `l`, computed
/// as in `research/cover_route_prime_binary_20261008/ghs_magic.rs`.
pub fn ghs_magic(
    n_total: u32,
    l: u32,
    irr: &IrreduciblePoly,
    a: &F2mElement,
    b: &F2mElement,
) -> GhsMagic {
    let n_rel = n_total / l;
    let gamma = b.square_k_times(n_total - 1, irr);
    let mut rows = Vec::with_capacity(n_rel as usize);
    let mut cur = gamma.clone();
    for _ in 0..n_rel {
        rows.push((cur.to_biguint() << 1u32) | BigUint::from(1u32));
        cur = cur.square_k_times(l, irr);
    }
    let magic_b = f2_rank(rows);
    let mut trace = a.clone();
    let mut acc = a.clone();
    for _ in 1..n_total {
        acc = acc.square(irr);
        trace = trace.add(&acc);
    }
    let trace_a = u32::from(!trace.is_zero());
    let a_may_add_one = trace_a == 1 && n_rel.is_multiple_of(2);
    let mut orbit_length = n_rel;
    for d in 1..=n_rel {
        if n_rel.is_multiple_of(d) && gamma.square_k_times(d * l, irr) == gamma {
            orbit_length = d;
            break;
        }
    }
    GhsMagic {
        magic_b,
        trace_a,
        a_may_add_one,
        magic: magic_b,
        orbit_length,
    }
}

/// Degrees (with multiplicity) of the irreducible factors of `x^n - 1` over
/// `F_2`: for `n = 2^e * n'` with `n'` odd, each 2-cyclotomic coset modulo
/// `n'` of size `d` contributes a factor of degree `d` with multiplicity `2^e`.
pub fn cyclotomic_factor_degrees(n: u32) -> Vec<(u32, u32)> {
    let mut e = 0u32;
    let mut odd = n;
    while odd.is_multiple_of(2) {
        odd /= 2;
        e += 1;
    }
    let mult = 1u32 << e;
    let mut seen = vec![false; odd as usize];
    let mut out = Vec::new();
    for start in 0..odd {
        if seen[start as usize] {
            continue;
        }
        let mut size = 0u32;
        let mut x = start;
        loop {
            seen[x as usize] = true;
            size += 1;
            x = (2 * x) % odd;
            if x == start {
                break;
            }
        }
        out.push((size, mult));
    }
    out.sort_unstable();
    out
}

/// `log2` of the number of monic divisors of `x^n - 1` over `F_2` of each
/// degree `0..=n`, `None` where there is none.
pub fn log2_divisor_counts(n: u32) -> Vec<Option<f64>> {
    let mut counts: Vec<f64> = vec![0.0; n as usize + 1];
    counts[0] = 1.0;
    for (degree, mult) in cyclotomic_factor_degrees(n) {
        let mut next = vec![0.0; n as usize + 1];
        for (deg, &c) in counts.iter().enumerate() {
            if c == 0.0 {
                continue;
            }
            for copies in 0..=mult {
                let d = deg + (copies * degree) as usize;
                if d > n as usize {
                    break;
                }
                next[d] += c;
            }
        }
        counts = next;
    }
    counts
        .into_iter()
        .map(|c| if c > 0.0 { Some(c.log2()) } else { None })
        .collect()
}

/// One heuristic isogeny-extension row for the tower with base degree `l`
/// and relative degree `n`.
#[derive(Clone, Debug, PartialEq)]
pub struct IsogenyExtensionEstimate {
    /// Annihilator degree `d` of `sqrt(b')` aimed at; the magic number of such a
    /// curve is at most `d + 1`.
    pub annihilator_degree: u32,
    /// Magic number bound `d + 1` used for the genus.
    pub magic_bound: u32,
    /// `log2` of the number of `b'` over the top field whose `sqrt(b')` has
    /// annihilator of degree exactly `d` for this tower (an upper bound:
    /// divisors of degree `d` times `2^(l d)`).
    pub log2_candidate_curves: f64,
    /// `log2` of the expected number of those curves in the target's isogeny
    /// class, with the class size taken as `sqrt(|t^2 - 4q|) / pi`.
    pub log2_expected_isogenous: f64,
    /// `log2` cost of walking the isogeny class until such a curve is hit:
    /// the reciprocal of the fraction of weak curves, times one isogeny step.
    pub log2_search_cost: f64,
}

/// Isogeny-extension estimates for one tower, for every achievable
/// annihilator degree `d` below `direct_magic` (a curve already reaching the
/// direct magic number gains nothing from a neighbour with the same one).
pub fn isogeny_extension_estimates(
    model: &CostModel,
    n_total: u32,
    l: u32,
    n_rel: u32,
    direct_magic: u32,
    log2_abs_disc: f64,
) -> Vec<IsogenyExtensionEstimate> {
    let counts = log2_divisor_counts(n_rel);
    let log2_class_size = 0.5 * log2_abs_disc - std::f64::consts::PI.log2();
    let log2_all_curves = f64::from(n_total) + 1.0;
    let step = model.log2_field_mult(f64::from(n_total)) + f64::from(n_total).log2();
    let mut out = Vec::new();
    for d in 1..n_rel.min(direct_magic.saturating_sub(1)) + 1 {
        let Some(log2_divisors) = counts.get(d as usize).copied().flatten() else {
            continue;
        };
        if d + 1 >= direct_magic {
            break;
        }
        let log2_candidate_curves = log2_divisors + f64::from(l) * f64::from(d);
        let log2_expected_isogenous = log2_candidate_curves + log2_class_size - log2_all_curves;
        let log2_search_cost = (log2_all_curves - log2_candidate_curves).max(0.0) + step;
        out.push(IsogenyExtensionEstimate {
            annihilator_degree: d,
            magic_bound: d + 1,
            log2_candidate_curves,
            log2_expected_isogenous,
            log2_search_cost,
        });
    }
    out
}

/// `2^(m-1)` as a decimal string, the generic GHS genus for magic number `m`.
pub fn generic_genus(magic: u32) -> BigUint {
    if magic == 0 {
        BigUint::from(0u32)
    } else {
        BigUint::from(1u32) << (magic - 1)
    }
}

/// `log2` of the total of a search phase and a solve phase.
pub fn log2_total(parts: &[f64]) -> f64 {
    parts.iter().copied().fold(f64::NEG_INFINITY, log2_add)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cyclotomic_degrees_match_known_factorisations() {
        // x^31 - 1 = (x + 1) * six quintics.
        assert_eq!(
            cyclotomic_factor_degrees(31),
            vec![(1, 1), (5, 1), (5, 1), (5, 1), (5, 1), (5, 1), (5, 1)]
        );
        // x^7 - 1 = (x + 1)(x^3 + x + 1)(x^3 + x^2 + 1).
        assert_eq!(cyclotomic_factor_degrees(7), vec![(1, 1), (3, 1), (3, 1)]);
        // x^6 - 1 = (x^3 - 1)^2.
        assert_eq!(cyclotomic_factor_degrees(6), vec![(1, 2), (2, 2)]);
        // Prime n with 2 a primitive root: x + 1 times one factor of degree n - 1.
        assert_eq!(cyclotomic_factor_degrees(131), vec![(1, 1), (130, 1)]);
        let total: u32 = cyclotomic_factor_degrees(155)
            .iter()
            .map(|(d, m)| d * m)
            .sum();
        assert_eq!(total, 155);
    }

    #[test]
    fn f2_rank_counts_independent_rows() {
        let r = |v: &[u32]| f2_rank(v.iter().map(|&x| BigUint::from(x)).collect());
        assert_eq!(r(&[0b101, 0b011, 0b110]), 2);
        assert_eq!(r(&[0b1, 0b10, 0b100]), 3);
        assert_eq!(r(&[0, 0b111, 0b111]), 1);
    }

    #[test]
    fn divisor_counts_are_binomial_for_distinct_factors() {
        // x^31 - 1: degrees 0,1,5,6,10,11,...; C(6, k) quintics.
        let c = log2_divisor_counts(31);
        assert!((c[5].unwrap() - 6f64.log2()).abs() < 1e-9);
        assert!((c[10].unwrap() - 15f64.log2()).abs() < 1e-9);
        assert!(c[2].is_none() && c[7].is_none());
        assert!((c[0].unwrap()).abs() < 1e-12);
        assert!((c[31].unwrap()).abs() < 1e-12);
    }

    #[test]
    fn isogeny_extension_for_oakley_group_three_is_hopeless_at_small_magic() {
        // N = 155, l = 5, n = 31: curves with annihilator degree 5 number
        // 6 * 2^25, the class has about 2^76 members, 2^156 curves in all.
        let model = CostModel::default();
        let rows = isogeny_extension_estimates(&model, 155, 5, 31, 26, 155.0);
        let d5 = rows.iter().find(|r| r.annihilator_degree == 5).unwrap();
        assert!((d5.log2_candidate_curves - (6f64.log2() + 25.0)).abs() < 1e-9);
        assert!(d5.log2_expected_isogenous < -40.0);
        assert!(d5.log2_search_cost > 120.0);
        assert!(rows.iter().all(|r| r.magic_bound < 26));
    }

    #[test]
    fn no_extension_rows_when_direct_magic_is_already_minimal() {
        let model = CostModel::default();
        assert!(isogeny_extension_estimates(&model, 155, 31, 5, 5, 155.0)
            .iter()
            .all(|r| r.magic_bound < 5));
        assert!(isogeny_extension_estimates(&model, 163, 1, 163, 1, 163.0).is_empty());
    }
}
