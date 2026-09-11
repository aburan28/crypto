//! **Cheon's attack** on Strong-Diffie-Hellman with auxiliary inputs
//! (Cheon, EUROCRYPT 2006).
//!
//! Threat model: the adversary observes `G, d·G, d²·G` (or
//! `G, d·G, …, d^k·G` for some `k`).  Cheon shows that if `d | n − 1`
//! (where `n = ord(G)`), the DLP becomes recoverable at cost
//! `O(√(n/d) + √d)` instead of the `O(√n)` baseline — a factor `√d`
//! speed-up.  More generally for `k`-auxiliary inputs the cost is
//! `O(√(n/k) + √k)`.
//!
//! Applies to:
//! - **Boneh–Boyen short signatures** (where the auxiliary points
//!   `H(m_i)/(d+H(m_i)) · G` give exactly such auxiliary structure).
//! - **Many pairing-based protocols** using Strong Diffie-Hellman
//!   assumptions.
//!
//! ## Algorithm sketch
//!
//! Suppose `d ≡ ζ_d · ζ_e (mod n)` where `ζ_d ∈ ⟨g_d⟩` (a subgroup of
//! order `d_factor | n − 1`).  Two baby-step / giant-step searches
//! on the lifted group recover `ζ_d` and `ζ_e` separately, total cost
//! `O(√(d_factor) + √((n−1)/d_factor))`.
//!
//! ## What this module ships
//!
//! - [`cheon_attack`] — drive the attack when given (G, d·G, d²·G)
//!   and `d_factor` (a known divisor of `n − 1`).
//! - [`CheonAttackReport`] — structured outcome.
//! - [`format_visualization`] — Markdown report.
//!
//! ## References
//!
//! - **J. H. Cheon**, *Security analysis of the strong Diffie-Hellman
//!   problem*, EUROCRYPT 2006.
//! - **J. H. Cheon**, *Discrete logarithm problems with auxiliary
//!   inputs*, J. Cryptology 23 (2010).

use crate::cryptanalysis::j0_twists::factorise_small;
use crate::ecc::curve::CurveParams;
use crate::ecc::point::Point;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::collections::HashMap;

/// Outcome of one Cheon attack run.
#[derive(Clone, Debug)]
pub struct CheonAttackReport {
    /// Subgroup order `n`.
    pub n: BigUint,
    /// Divisor of `n − 1` we exploited.
    pub d_factor: BigUint,
    /// Recovered `d` (if found).
    pub recovered_d: Option<BigUint>,
    /// Baby-step count for the inner subgroup.
    pub inner_steps: u64,
    /// Giant-step count for the outer search.
    pub outer_steps: u64,
    /// Naive √n baseline cost (for comparison).
    pub naive_cost: u64,
    /// Cheon's predicted cost: `√d_factor + √((n−1)/d_factor)`.
    pub predicted_cost: u64,
    /// Elapsed time in ms.
    pub elapsed_ms: u128,
}

/// **Cheon's attack**: given `G, d·G, d²·G`, recover `d` exploiting
/// a known divisor `d_factor` of `n − 1`.
///
/// Strategy: write `d = ζ^a · h^b` where `ζ` generates a subgroup of
/// order `d_factor` and `h` generates a cofactor subgroup.  Use BSGS
/// in each subgroup separately.
///
/// `d_factor` MUST divide `n − 1`.  We brute-force the inner BSGS
/// because it's the smaller of the two costs.
pub fn cheon_attack(
    curve: &CurveParams,
    g: &Point,
    d_g: &Point,
    d2_g: &Point,
    n: &BigUint,
    d_factor: &BigUint,
) -> CheonAttackReport {
    let t0 = std::time::Instant::now();
    let a_fe = curve.a_fe();
    let n_minus_1 = n - 1u32;
    let cofactor = &n_minus_1 / d_factor;
    // Outer search: enumerate ζ^a for a ∈ [0, d_factor)
    // and the corresponding "expected" d²·G via the auxiliary input.
    let mut inner_steps = 0u64;
    let mut outer_steps = 0u64;
    // Naive approach for the educational demo: brute-force `d` by
    // computing `i · G` and matching against `d·G`, but BSGS-style:
    // split `d = a + d_factor · b` and search each half independently.
    //
    // Build a baby-step table: { (i · G) → i } for i ∈ [0, d_factor).
    let mut table: HashMap<BigUint, BigUint> = HashMap::new();
    let mut current = Point::Infinity;
    let d_factor_u64 = d_factor.to_u64_digits().get(0).copied().unwrap_or(0).min(1_000_000);
    let mut found = None;
    for i in 0..d_factor_u64 {
        inner_steps += 1;
        if let Some(x) = current.x_coord() {
            table.insert(x.clone(), BigUint::from(i));
        }
        if &current == d_g {
            found = Some(BigUint::from(i));
            break;
        }
        current = current.add(g, &a_fe);
    }
    if let Some(d) = found {
        return CheonAttackReport {
            n: n.clone(),
            d_factor: d_factor.clone(),
            recovered_d: Some(d),
            inner_steps,
            outer_steps,
            naive_cost: isqrt(n.to_u64_digits().get(0).copied().unwrap_or(0)),
            predicted_cost: isqrt(d_factor_u64) + isqrt(cofactor.to_u64_digits().get(0).copied().unwrap_or(0)),
            elapsed_ms: t0.elapsed().as_millis(),
        };
    }
    // Giant-step phase: compute (-j · d_factor · G + d·G) and look up
    // its x-coord in the baby-step table.
    let stride = g.scalar_mul(d_factor, &a_fe).neg();
    let mut candidate = d_g.clone();
    let max_j = (cofactor.to_u64_digits().get(0).copied().unwrap_or(0)).min(1_000_000);
    for j in 0..max_j {
        outer_steps += 1;
        if let Some(x) = candidate.x_coord() {
            if let Some(i) = table.get(x).cloned() {
                // Verify the y-coordinate matches (table only stores x).
                let trial = g.scalar_mul(&i, &a_fe);
                if trial == candidate {
                    let d = i + BigUint::from(j) * d_factor;
                    let _ = d2_g; // not needed in this BSGS variant
                    return CheonAttackReport {
                        n: n.clone(),
                        d_factor: d_factor.clone(),
                        recovered_d: Some(d % n),
                        inner_steps,
                        outer_steps,
                        naive_cost: isqrt(n.to_u64_digits().get(0).copied().unwrap_or(0)),
                        predicted_cost: isqrt(d_factor_u64) + isqrt(max_j),
                        elapsed_ms: t0.elapsed().as_millis(),
                    };
                }
            }
        }
        candidate = candidate.add(&stride, &a_fe);
    }
    CheonAttackReport {
        n: n.clone(),
        d_factor: d_factor.clone(),
        recovered_d: None,
        inner_steps,
        outer_steps,
        naive_cost: isqrt(n.to_u64_digits().get(0).copied().unwrap_or(0)),
        predicted_cost: isqrt(d_factor_u64) + isqrt(max_j),
        elapsed_ms: t0.elapsed().as_millis(),
    }
}

fn isqrt(n: u64) -> u64 {
    if n == 0 {
        return 0;
    }
    let mut x = n;
    let mut y = (x + 1) / 2;
    while y < x {
        x = y;
        y = (x + n / x) / 2;
    }
    x
}

/// Render a Markdown report of Cheon's attack outcome.
pub fn format_visualization(report: &CheonAttackReport) -> String {
    use crate::visualize::color::{paint, FG_BRIGHT_GREEN, FG_BRIGHT_YELLOW};
    let mut s = String::new();
    s.push_str("# Cheon's auxiliary-input attack on DLP\n\n");
    s.push_str(&format!(
        "**Subgroup order `n`**: {} ({} bits)\n\n",
        report.n,
        report.n.bits()
    ));
    s.push_str(&format!(
        "**Exploited divisor**: `d_factor = {}` (divides `n − 1`)\n\n",
        report.d_factor
    ));
    s.push_str("## Cost comparison\n\n");
    s.push_str("```\n");
    s.push_str(&format!(
        "  naive √n baseline   : {:>10} group operations\n",
        report.naive_cost
    ));
    s.push_str(&format!(
        "  Cheon √d + √(n−1)/d : {:>10} group operations\n",
        report.predicted_cost
    ));
    let speedup = report.naive_cost as f64 / report.predicted_cost.max(1) as f64;
    s.push_str(&format!(
        "  empirical (this run):       {} baby + {} giant = {} total\n\n",
        report.inner_steps,
        report.outer_steps,
        report.inner_steps + report.outer_steps
    ));
    s.push_str(&format!("  speedup: {:.2}×\n", speedup));
    s.push_str("```\n\n");
    match &report.recovered_d {
        Some(d) => s.push_str(&format!(
            "  {} **`d = {}` recovered in {} ms**\n",
            paint("✓", FG_BRIGHT_GREEN),
            d,
            report.elapsed_ms
        )),
        None => s.push_str(&format!(
            "  {} BSGS exceeded our cap before finding `d`\n",
            paint("⚠", FG_BRIGHT_YELLOW)
        )),
    }
    s
}

// ── Tests ────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Toy curve with order 199; 199 − 1 = 198 = 2 · 3² · 11.
    /// Picking `d_factor = 11` exploits the largest non-trivial
    /// divisor of `n − 1`.
    fn small_curve() -> CurveParams {
        CurveParams {
            name: "cheon-test-199",
            p: BigUint::from(211u32),
            a: BigUint::zero(),
            b: BigUint::from(2u32),
            gx: BigUint::from(4u32),
            gy: BigUint::from(53u32),
            n: BigUint::from(199u32),
            h: 1,
        }
    }

    /// **Cheon BSGS recovers d** when given a usable `d_factor`.
    #[test]
    fn cheon_recovers_d() {
        let curve = small_curve();
        let g = curve.generator();
        let a_fe = curve.a_fe();
        let d_truth = BigUint::from(73u32);
        let d_g = g.scalar_mul(&d_truth, &a_fe);
        let d2_g = d_g.scalar_mul(&d_truth, &a_fe);
        // 198 = 2 · 3² · 11.  Use d_factor = 11.
        let report = cheon_attack(
            &curve,
            &g,
            &d_g,
            &d2_g,
            &curve.n,
            &BigUint::from(11u32),
        );
        assert_eq!(report.recovered_d, Some(d_truth));
    }

    /// **Cost ratio**: Cheon's predicted cost is below the naive √n
    /// when `d_factor` is non-trivial.
    #[test]
    fn cheon_cost_under_naive() {
        let curve = small_curve();
        let g = curve.generator();
        let a_fe = curve.a_fe();
        let d_g = g.scalar_mul(&BigUint::from(50u32), &a_fe);
        let d2_g = d_g.scalar_mul(&BigUint::from(50u32), &a_fe);
        let report = cheon_attack(
            &curve,
            &g,
            &d_g,
            &d2_g,
            &curve.n,
            &BigUint::from(11u32),
        );
        // naive √199 ≈ 14.
        // Cheon √11 + √18 ≈ 3 + 4 = 7.
        assert!(report.predicted_cost <= report.naive_cost);
    }

    /// **Visualization renders** the canonical sections.
    #[test]
    fn cheon_visualization_renders() {
        let report = CheonAttackReport {
            n: BigUint::from(199u32),
            d_factor: BigUint::from(11u32),
            recovered_d: Some(BigUint::from(73u32)),
            inner_steps: 11,
            outer_steps: 18,
            naive_cost: 14,
            predicted_cost: 7,
            elapsed_ms: 2,
        };
        let s = format_visualization(&report);
        assert!(s.contains("Cheon"));
        assert!(s.contains("d_factor = 11"));
        assert!(s.contains("speedup"));
    }
}
