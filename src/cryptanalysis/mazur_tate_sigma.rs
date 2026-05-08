//! **Mazur-Tate `p`-adic sigma function** for ECDLP recovery via
//! `p`-adic heights (Phase 3 of Kim-program / nonabelian Chabauty).
//!
//! Mazur-Tate 1991 introduced the `p`-adic analog of the Weierstrass
//! σ-function for elliptic curves with good ordinary reduction.  Key
//! properties:
//!
//! 1. `σ_p(P) ∈ Q_p` is defined for `P ∈ E(Q_p) − {O}`.
//! 2. `σ_p` has a simple zero at `O`: `σ_p(z) = z + O(z³)` near `O`,
//!    where `z = -x/y` is the formal-group parameter.
//! 3. **The defining functional equation** (Mazur-Tate):
//!
//!    ```text
//!    σ_p(P + Q) · σ_p(P − Q) / σ_p(P)² σ_p(Q)² = (x(P) − x(Q))^{-1}
//!    ```
//!
//!    (or similar; conventions vary).
//!
//! 4. **The bilinear height pairing**:
//!
//!    ```text
//!    h_p(P, Q) := log_p[σ_p(P + Q) σ_p(P − Q) / σ_p(P)² σ_p(Q)²] - log_p[(x(P) − x(Q))^{-1}]
//!    ```
//!
//!    is bilinear and symmetric in `P, Q`.
//!
//! 5. **The doubling-quadratic relation**:
//!
//!    ```text
//!    log_p σ_p(d · P) - d · log_p σ_p(P) = h_p(P, P) · d(d-1)/2
//!    ```
//!
//!    which is **quadratic in d**.
//!
//! # The proposed ECDLP attack
//!
//! Given `P, Q ∈ E(F_p)` with `Q = d · P` (planted), lift to
//! `P̂, Q̂ ∈ E(Z_p)` via Hensel.  Compute:
//!
//! ```text
//! L_P := log_p σ_p(P̂)
//! L_Q := log_p σ_p(Q̂)
//! H   := h_p(P̂, P̂)     (computable from P̂ alone via height formula)
//! ```
//!
//! Then:
//!
//! ```text
//! L_Q - d · L_P = H · d(d-1)/2
//!         d² · H/2 - d · (L_P + H/2) + L_Q = 0
//! ```
//!
//! Solve quadratic for `d` mod `p`.  Two roots; need a second
//! constraint to pick the correct one (e.g., reduction mod a small
//! prime or trial verification).
//!
//! # Why this might work where the canonical-lift Smart attack failed
//!
//! The Smart attack uses only `log_F` (the abelian formal-group
//! logarithm), which doesn't capture the height information.
//! Mazur-Tate `σ_p` adds the **height pairing** `h_p`, which is a
//! genuinely non-abelian invariant: it depends on **products** of
//! coordinates (`x(P) - x(Q)` term in the functional equation) in
//! a way that doesn't factor through the formal group.
//!
//! As far as the agent's literature survey can determine, **this
//! attack idea has never been published**.  It would be a new
//! application of Mazur-Tate to ECDLP.
//!
//! # Implementation status — HONEST CAVEATS
//!
//! This module ships a **simplified σ_p stub**, NOT the full
//! Mazur-Tate σ-function.  The differences:
//!
//! - The actual Mazur-Tate σ is constructed via `p`-adic Eisenstein
//!   series `E_2, E_4, E_6` and satisfies the functional equation
//!   `σ(P+Q)σ(P-Q) = σ(P)²σ(Q)² · (x(P)-x(Q))^{-1}` (or similar).
//! - This module's `σ_p_series` uses the simple expansion
//!   `σ ≈ z + (a/5)z⁵ + (b/7)z⁷ + ...`, which is the formal-group
//!   parameter, NOT the Mazur-Tate σ.
//!
//! As a consequence: the `mazur_tate_quadratic_recovery_probe` test
//! finds that estimated `h(P, P)` values are NOT consistent across
//! `d`, indicating that the simplified σ does NOT satisfy the
//! Mazur-Tate height-pairing relation.
//!
//! **Empirical conclusion**: the simplified σ does not yield ECDLP
//! recovery.  This is an **honest negative result** at the level
//! of the simplified implementation.
//!
//! **What's needed for a full attempt**:
//! 1. Implement `p`-adic Eisenstein series `E_2(τ)`, `E_4(τ)`, `E_6(τ)`
//!    with `τ = log_p` of the formal-group parameter.
//! 2. Construct the actual Mazur-Tate σ using these.
//! 3. Verify the functional equation empirically.
//! 4. Run the recovery test.
//!
//! Each step is several hundred lines of code and substantial
//! mathematics.  This is left as future work.
//!
//! # The conceptual framework — still valid
//!
//! Even with the simplified implementation's negative result, the
//! **conceptual idea** of using Mazur-Tate σ for ECDLP recovery
//! via the quadratic relation `log σ(d·P) - d·log σ(P) = h·d(d-1)/2`
//! remains the most-promising direction.  The negative result here
//! says only that the simplified stub doesn't work; it does NOT
//! say the full attack doesn't work.
//!
//! As far as the agent's literature survey can determine, this
//! attack idea has never been published as a serious ECDLP-attack
//! proposal.  The infrastructure to attempt it is the natural
//! continuation of this codebase.

use crate::cryptanalysis::canonical_lift::{ZpCurve, ZpInt};
use crate::cryptanalysis::coleman_integration::PSeries;
use num_bigint::BigInt;
use num_traits::{One, Zero};

/// `σ_p(z)` as a `z`-power-series for short-Weierstrass `E: y² = x³
/// + ax + b`.
///
/// **Derivation** (Cohen "Number Theory I" §7.4.2 / Silverman III.5):
/// for short Weierstrass with invariants `g_2 = -4a`, `g_3 = -4b`,
/// the Weierstrass `℘`-function expands as:
///
/// ```text
/// ℘(z) = z⁻² - (a/5) z² - (b/7) z⁴ - …
/// ```
///
/// The Weierstrass `ζ`-function `ζ' = -℘ + z⁻²` integrates to:
///
/// ```text
/// ζ(z) - 1/z = ∫₀ᶻ -℘(t) dt - (-1/z + 1/0)  (regularised)
///            = (a/15) z³ + (b/35) z⁵ + …
/// ```
///
/// Then `σ(z) = z · exp[∫(ζ - 1/z) dz]`:
///
/// ```text
/// σ(z) = z · exp[(a/60) z⁴ + (b/210) z⁶ + …]
///      = z + (a/60) z⁵ + (b/210) z⁷ + (a²/7200) z⁹ + …
/// ```
///
/// **This is the corrected formula** (the previous version used
/// `a/5`, `b/7`, which were the `℘`-coefficients, not the `σ`-
/// coefficients).
///
/// Returns the power series in `z` with `n` terms.
pub fn sigma_p_series(curve: &ZpCurve, n: usize) -> PSeries {
    let p = &curve.p;
    let prec = curve.precision;
    let mut coefs = vec![ZpInt::zero(p, prec); n];
    if n >= 2 {
        // σ_p(z) = z + ...
        coefs[1] = ZpInt::one(p, prec);
    }
    if n >= 6 {
        // Coefficient of z⁵: a / 60.
        let sixty = ZpInt::new(BigInt::from(60), p, prec);
        if let Some(inv60) = sixty.inverse() {
            coefs[5] = curve.a.mul(&inv60);
        }
    }
    if n >= 8 {
        // Coefficient of z⁷: b / 210.
        let two_ten = ZpInt::new(BigInt::from(210), p, prec);
        if let Some(inv) = two_ten.inverse() {
            coefs[7] = curve.b.mul(&inv);
        }
    }
    if n >= 10 {
        // Coefficient of z⁹: a² / 7200.
        // (From exp expansion: (a/60)²/2 z⁸·z = a²/7200 z⁹.)
        let seventy_two_hundred = ZpInt::new(BigInt::from(7200), p, prec);
        if let Some(inv) = seventy_two_hundred.inverse() {
            coefs[9] = curve.a.mul(&curve.a).mul(&inv);
        }
    }
    PSeries::from_coef_vec(p, prec, coefs)
}

/// `log_p σ_p(z)` as a power-series in `z`.  Defined as the natural
/// `p`-adic logarithm of `σ_p(z) / z` (since `σ_p` has a simple zero
/// at `O`, so `σ_p(z) / z = 1 + O(z²)` and `log` converges).
///
/// `log_p(1 + u) = u - u²/2 + u³/3 - …` truncated to `n` terms.
pub fn log_sigma_minus_log_z(curve: &ZpCurve, n: usize) -> Option<PSeries> {
    let p = &curve.p;
    let prec = curve.precision;
    // u(z) = σ_p(z) / z - 1.  Compute as power series.
    let sigma = sigma_p_series(curve, n + 2);
    // Strip the leading z: divide by z (= shift coefficients).
    let mut u_coefs = vec![ZpInt::zero(p, prec); n];
    for i in 0..n {
        if i + 1 < sigma.coefs.len() {
            u_coefs[i] = sigma.coefs[i + 1].clone();
        }
    }
    // Subtract 1 (the constant term of σ_p/z is 1, since σ_p starts
    // at z).
    u_coefs[0] = u_coefs[0].sub(&ZpInt::one(p, prec));
    let u = PSeries::from_coef_vec(p, prec, u_coefs);
    // log(1 + u) = u - u²/2 + u³/3 - …
    let mut result = PSeries::new(p, prec, n);
    let mut u_pow = u.clone();
    let mut sign: i64 = 1;
    for k in 1..n.min(20) {
        let k_zp = ZpInt::new(BigInt::from(k as i64), p, prec);
        let k_inv = k_zp.inverse()?;
        let term = u_pow.mul(&PSeries::from_coef_vec(p, prec, vec![
            ZpInt::new(BigInt::from(sign), p, prec).mul(&k_inv),
        ]));
        result = result.add(&term);
        u_pow = u_pow.mul(&u);
        // Truncate to keep length manageable.
        if u_pow.n() > n + 3 {
            u_pow.coefs.truncate(n + 3);
        }
        sign = -sign;
    }
    Some(result)
}

/// **The headline ECDLP-recovery test**: given `P, Q = d·P` in the
/// formal group, attempt to recover `d` via the quadratic relation:
///
/// ```text
/// log σ(d·P) - d · log σ(P) = h(P, P) · d(d-1)/2
/// ```
///
/// At leading order in `z`, `σ(z) = z + O(z⁵)` for short Weierstrass,
/// so `log σ(z) - log z = O(z⁴) / z = O(z⁴)` → vanishes at first
/// order.  The height correction `h(P, P)` enters at higher order.
///
/// For toy curves at low precision, this might or might not have
/// the d-recovery signal.  The experiment: tabulate empirically.
pub fn run_mazur_tate_recovery_test(
    curve: &ZpCurve,
    z_p: &ZpInt,
    d: u64,
) -> Option<MazurTateResult> {
    // Series length must stay below p (so we don't try to invert
    // p in any antiderivative).  For p = 11, use n = 10.
    let p_u: u64 = curve.p.clone().try_into().unwrap_or(11);
    let series_n = (p_u - 1) as usize;
    let log_sigma_minus_log_z = log_sigma_minus_log_z(curve, series_n)?;

    // For Q = d·P in formal group: z_Q = d · z_P (additive at leading
    // order for short Weierstrass).
    let z_q = ZpInt::new(BigInt::from(d as i64), &curve.p, curve.precision).mul(z_p);

    // log σ(P) = log z_P + log_sigma_minus_log_z(z_P)
    // log σ(Q) = log z_Q + log_sigma_minus_log_z(z_Q)
    //
    // Subtraction: log σ(Q) - d · log σ(P)
    //            = log z_Q - d · log z_P + (correction)
    //            = log(d · z_P) - d · log z_P + (correction)
    //            = log d + log z_P - d · log z_P + (correction)
    //            = log d + (1 - d) · log z_P + (correction)
    //
    // For toy: log z_P doesn't make sense in p-adic for v_p(z_P) ≥ 1
    // unless we treat it formally.  The interesting part is the
    // (correction) term, which is determined by σ_p's structure.

    // Pragmatic: just compute log_sigma_minus_log_z at z_P and z_Q,
    // and look at the relation log σ(Q) - d · log σ(P) under the
    // formal-group additive law.
    let val_p = log_sigma_minus_log_z.evaluate(z_p);
    let val_q = log_sigma_minus_log_z.evaluate(&z_q);

    let d_zp = ZpInt::new(BigInt::from(d as i64), &curve.p, curve.precision);
    let predicted_via_d = d_zp.mul(&val_p);
    let actual = val_q;
    let diff = actual.sub(&predicted_via_d);

    // For Mazur-Tate quadratic: diff should equal h(P, P) · d(d-1)/2
    // (some constant times d(d-1)/2).
    // The relationship: diff / [d(d-1)/2] should be independent of d.
    let d_d_minus_1 = ZpInt::new(BigInt::from((d * (d - 1)) as i64), &curve.p, curve.precision);
    let two = ZpInt::new(BigInt::from(2), &curve.p, curve.precision);
    let half_d_d_minus_1 = d_d_minus_1.mul(&two.inverse()?);
    let estimated_h = if half_d_d_minus_1.value.is_zero() {
        None
    } else {
        Some(diff.mul(&half_d_d_minus_1.inverse()?))
    };

    Some(MazurTateResult {
        d,
        z_p: z_p.value.clone(),
        z_q: z_q.value,
        val_p: val_p.value,
        val_q: actual.value,
        diff: diff.value,
        estimated_h: estimated_h.map(|h| h.value),
    })
}

#[derive(Clone, Debug)]
pub struct MazurTateResult {
    pub d: u64,
    pub z_p: BigInt,
    pub z_q: BigInt,
    pub val_p: BigInt,
    pub val_q: BigInt,
    pub diff: BigInt,
    pub estimated_h: Option<BigInt>,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `σ_p(z)` has expected leading coefficients.
    #[test]
    fn sigma_p_series_leading_terms() {
        let p = BigInt::from(11);
        let prec = 6u32;
        let curve = ZpCurve::new(BigInt::from(1), BigInt::from(1), &p, prec);
        let sigma = sigma_p_series(&curve, 10);
        // σ_p(z) = z + 0·z² + 0·z³ + 0·z⁴ + (a/5)·z⁵ + 0·z⁶ + (b/7)·z⁷ + ...
        assert_eq!(sigma.coefs[0].value, BigInt::zero());
        assert_eq!(sigma.coefs[1].value, BigInt::one());
        assert_eq!(sigma.coefs[2].value, BigInt::zero());
        assert_eq!(sigma.coefs[3].value, BigInt::zero());
        assert_eq!(sigma.coefs[4].value, BigInt::zero());
        // a/5 mod 11^6: a=1, 5^{-1} mod 11^6.
        // We'll just check it's non-zero.
        assert!(!sigma.coefs[5].value.is_zero());
    }

    /// **The Mazur-Tate ECDLP recovery probe**: tabulate the
    /// quadratic-correction signal across many `d`.  See whether
    /// the estimated `h(P, P)` is consistent across `d` (which
    /// would confirm the relation; if so, `d` is recoverable).
    #[test]
    fn mazur_tate_quadratic_recovery_probe() {
        let p = BigInt::from(11);
        let prec = 6u32;
        let curve = ZpCurve::new(BigInt::from(1), BigInt::from(1), &p, prec);

        // Pick z_P with v_p(z_P) = 1.
        let z_p = ZpInt::new(BigInt::from(11), &p, prec);

        println!();
        println!("=== Mazur-Tate p-adic σ ECDLP recovery probe ===");
        println!();
        println!("Hypothesis: log σ(d·P) - d·log σ(P) = h(P,P) · d(d-1)/2");
        println!("If estimated_h is constant across d, the relation is exact");
        println!("and d is recoverable as (positive root of quadratic).");
        println!();
        println!("{:>3} {:>12} {:>15} {:>15} {:>15}", "d", "diff", "estimated h", "z_q", "(d-1)·val_p");

        let mut estimated_hs: Vec<BigInt> = Vec::new();
        for d in 2u64..=8 {
            if let Some(r) = run_mazur_tate_recovery_test(&curve, &z_p, d) {
                let h_str = r.estimated_h.as_ref()
                    .map(|h| h.to_string())
                    .unwrap_or_else(|| "N/A".into());
                println!(
                    "{:>3} {:>12} {:>15} {:>15} {:>15}",
                    r.d,
                    r.diff.to_string(),
                    h_str,
                    r.z_q.to_string(),
                    "—"
                );
                if let Some(h) = r.estimated_h {
                    estimated_hs.push(h);
                }
            }
        }

        // Check consistency: are all estimated h's equal?
        let consistent = if estimated_hs.is_empty() {
            false
        } else {
            estimated_hs.iter().all(|h| h == &estimated_hs[0])
        };

        println!();
        if consistent && !estimated_hs.is_empty() {
            println!("✓ ALL estimated h(P, P) values are CONSISTENT.");
            println!("  Common h = {}", estimated_hs[0]);
            println!("  This confirms the Mazur-Tate quadratic relation at toy scale.");
            println!("  ECDLP recovery via the quadratic formula is feasible HERE.");
            println!();
            println!("  Caveat: this is for points GENERATED in formal group via");
            println!("  z_Q = d · z_P (additive law).  Real attack on F_p points");
            println!("  requires lifting and computing σ at non-formal-group points.");
        } else {
            println!("✗ Estimated h(P, P) varies with d: {:?}", estimated_hs);
            println!("  Either (a) the simplified formula is missing terms, or");
            println!("  (b) precision artifacts dominate, or");
            println!("  (c) the relation requires more careful regularization.");
        }
    }
}
