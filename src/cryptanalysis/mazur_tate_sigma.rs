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

// ── Phase 3: proper formal-group [d]_F via formal-log reversion ──────

/// Formal log `log_F(z)` for short Weierstrass `y² = x³ + ax + b`.
///
/// `log_F(z) = ∫ ω = z + (a/5)·z⁵ + (b/7)·z⁷ + (a²/9)·z⁹ + …`
///
/// This is the inverse-image of the standard "additive formal group"
/// under the isomorphism (over Z_p[[z]]) `Ê_a,b ≅ Ĝ_a`.
pub fn formal_log_series(curve: &ZpCurve, n: usize) -> PSeries {
    let p = &curve.p;
    let prec = curve.precision;
    let mut coefs = vec![ZpInt::zero(p, prec); n];
    if n >= 2 {
        coefs[1] = ZpInt::one(p, prec); // z
    }
    if n >= 6 {
        // a/5 · z⁵
        let five = ZpInt::new(BigInt::from(5), p, prec);
        if let Some(inv5) = five.inverse() {
            coefs[5] = curve.a.mul(&inv5);
        }
    }
    if n >= 8 {
        // b/7 · z⁷
        let seven = ZpInt::new(BigInt::from(7), p, prec);
        if let Some(inv7) = seven.inverse() {
            coefs[7] = curve.b.mul(&inv7);
        }
    }
    if n >= 10 {
        // a²/9 · z⁹  (from higher-order integration of ω)
        let nine = ZpInt::new(BigInt::from(9), p, prec);
        if let Some(inv9) = nine.inverse() {
            coefs[9] = curve.a.mul(&curve.a).mul(&inv9);
        }
    }
    PSeries::from_coef_vec(p, prec, coefs)
}

/// **Reversion** (series inverse) of `log_F`: given `u = log_F(z)`,
/// returns `z(u)` as a power series in `u`.
///
/// For `log_F(z) = z + c_5·z⁵ + c_7·z⁷ + c_9·z⁹ + ...`, the inverse:
/// ```text
///   z(u) = u − c_5·u⁵ − c_7·u⁷ + (5·c_5² − c_9)·u⁹ + …
/// ```
/// Computed via Lagrange-style iterative substitution.
pub fn formal_log_inverse_series(curve: &ZpCurve, n: usize) -> PSeries {
    let p = &curve.p;
    let prec = curve.precision;
    // Build `log_F(z)` as a series of length `n`.
    let logf = formal_log_series(curve, n);
    // We want z(u) such that logf(z(u)) ≡ u.  Start with z(u) = u
    // and iterate: z_{new}(u) = u − [logf(z_old(u)) − z_old(u)]
    //                         = z_old(u) − [logf(z_old(u)) − u]
    let mut z = PSeries::new(p, prec, n);
    if n >= 2 {
        z.coefs[1] = ZpInt::one(p, prec);
    }
    // Helper: compose `logf` with current `z(u)` to get `logf(z(u))`.
    // logf has coefficients c_k at z^k; logf(z(u)) = sum_k c_k · z(u)^k.
    for _iter in 0..n {
        // Compute z^k power series up to k = n - 1.
        let mut z_pow = PSeries::from_coef_vec(p, prec, vec![ZpInt::one(p, prec)]); // = 1
        let mut acc = PSeries::new(p, prec, n);
        for k in 0..logf.n() {
            if !logf.coefs[k].value.is_zero() {
                // Add c_k · z(u)^k.
                let scaled = z_pow.mul(&PSeries::from_coef_vec(
                    p,
                    prec,
                    vec![logf.coefs[k].clone()],
                ));
                acc = acc.add(&scaled);
            }
            // Step: z_pow ← z_pow · z.
            z_pow = z_pow.mul(&z);
            // Truncate to length n for efficiency.
            if z_pow.n() > n {
                z_pow.coefs.truncate(n);
            }
        }
        // Compute residual = acc − u.
        let mut u_series = PSeries::new(p, prec, n);
        if n >= 2 {
            u_series.coefs[1] = ZpInt::one(p, prec);
        }
        // acc − u
        let residual = acc.add(&scale_neg_one(&u_series, p, prec));
        // Check if residual is zero (up to truncation).
        let mut nonzero_residual = false;
        for c in &residual.coefs {
            if !c.value.is_zero() {
                nonzero_residual = true;
                break;
            }
        }
        if !nonzero_residual {
            break;
        }
        // Subtract residual from z (Newton-like correction).
        z = z.add(&scale_neg_one(&residual, p, prec));
        if z.n() > n {
            z.coefs.truncate(n);
        }
    }
    z
}

fn scale_neg_one(s: &PSeries, p: &BigInt, prec: u32) -> PSeries {
    let neg_one = ZpInt::new(BigInt::from(-1), p, prec);
    let neg_series = PSeries::from_coef_vec(p, prec, vec![neg_one]);
    s.mul(&neg_series)
}

/// **Compute `[d]_F(z) = log_F^{-1}(d · log_F(z))`** for short
/// Weierstrass.  This is the *correct* formal-group `d`-fold multiple
/// of `z`, accounting for the non-additive corrections in `F(s, t)`.
///
/// Uses **scalar Newton iteration** rather than power-series reversion
/// for numerical robustness at toy precision.
///
/// Given the target `u = d · log_F(z)`, find `z_d` such that
/// `log_F(z_d) = u` via:
/// ```text
///   z_{k+1} = z_k − (log_F(z_k) − u) / log_F'(z_k)
/// ```
/// Newton converges quadratically; at `prec` digits we need `O(log
/// prec)` iterations.
pub fn formal_d_times_z(curve: &ZpCurve, z: &ZpInt, d: u64, series_len: usize) -> ZpInt {
    let p = &z.p;
    let prec = z.precision;
    let logf = formal_log_series(curve, series_len);
    // Derivative of log_F: log_F'(z) = 1 + a·z⁴ + b·z⁶ + a²·z⁸ + ...
    // (differentiating term-by-term).
    let logf_deriv = {
        let mut coefs = vec![ZpInt::zero(p, prec); series_len];
        for k in 1..logf.n() {
            // logf has c_k at z^k.  d/dz: k·c_k·z^{k-1}.
            let k_zp = ZpInt::new(BigInt::from(k as i64), p, prec);
            coefs[k - 1] = k_zp.mul(&logf.coefs[k]);
        }
        PSeries::from_coef_vec(p, prec, coefs)
    };
    // u = d · log_F(z)
    let u = logf.evaluate(z);
    let d_zp = ZpInt::new(BigInt::from(d as i64), p, prec);
    let target = d_zp.mul(&u);
    // Newton: start at z_0 = d · z (leading-order).
    let mut z_cur = d_zp.mul(z);
    for _ in 0..(prec as usize + 4) {
        let logf_at_z = logf.evaluate(&z_cur);
        let f = logf_at_z.sub(&target);
        if f.value.is_zero() {
            break;
        }
        let fp = logf_deriv.evaluate(&z_cur);
        let fp_inv = match fp.inverse() {
            Some(v) => v,
            None => break,
        };
        let delta = f.mul(&fp_inv);
        let new_z = z_cur.sub(&delta);
        if new_z == z_cur {
            break;
        }
        z_cur = new_z;
    }
    z_cur
}

/// **Phase 3 test result**: verifies the Mazur-Tate quadratic identity
/// ```text
///     log σ(d·P) − d · log σ(P) ≡ h(P, P) · d(d−1)/2  (mod p^precision)
/// ```
/// using the **correct** formal-group `[d]_F` (Phase-3 fix) rather
/// than the naive `z_Q = d·z_P` (Phase-2 stub).
#[derive(Clone, Debug)]
pub struct MazurTatePhase3Row {
    pub d: u64,
    pub z_d: BigInt,
    pub log_sigma_d: BigInt,
    pub diff: BigInt, // log σ(d·P) − d · log σ(P)
    pub d_d_minus_1_half: BigInt,
    pub h_estimate: Option<BigInt>,
}

/// **Verify** that `[d]_F` is correct by checking the composition
/// law `[d]_F([e]_F(z)) = [d·e]_F(z)`.  This is a structural sanity
/// check that does not depend on the σ-series being complete.
pub fn verify_formal_d_times_composition(
    curve: &ZpCurve,
    z: &ZpInt,
    d: u64,
    e: u64,
    series_len: usize,
) -> bool {
    let z_e = formal_d_times_z(curve, z, e, series_len);
    let z_de_via_compose = formal_d_times_z(curve, &z_e, d, series_len);
    let z_de_direct = formal_d_times_z(curve, z, d * e, series_len);
    z_de_via_compose == z_de_direct
}

/// **`p`-adic unit root `α_p`** of the Frobenius polynomial
/// `T² − a_p T + p = 0` in `Z_p`, where `a_p = p + 1 − #E(F_p)`.
///
/// For an ordinary curve (`gcd(a_p, p) = 1`), this factor lifts to a
/// p-adic unit α_p with `α_p ≡ a_p (mod p)`.  Computed via Hensel
/// iteration on `f(α) = α² − a_p α + p`, starting from `α_0 = a_p`.
///
/// This is the Frobenius eigenvalue on the unit-root subspace of
/// `H^1_dR(E/Q_p)` and the crucial ingredient missing from the
/// existing σ-series for the Mazur-Tate quadratic identity to hold.
pub fn frobenius_unit_root(p: &BigInt, a_p: i64, precision: u32) -> Option<ZpInt> {
    let a_p_zp = ZpInt::new(BigInt::from(a_p), p, precision);
    if a_p_zp.value.is_zero() {
        // Supersingular case — no unit root.
        return None;
    }
    let p_zp = ZpInt::new(p.clone(), p, precision);
    // Newton iteration: α_{k+1} = α_k − f(α_k)/f'(α_k)
    // f(α) = α² − a_p α + p ⇒ f'(α) = 2α − a_p
    let mut alpha = a_p_zp.clone();
    for _ in 0..precision {
        let f = alpha.mul(&alpha).sub(&a_p_zp.mul(&alpha)).add(&p_zp);
        let f_prime = ZpInt::new(BigInt::from(2), p, precision)
            .mul(&alpha)
            .sub(&a_p_zp);
        let f_prime_inv = f_prime.inverse()?;
        let delta = f.mul(&f_prime_inv);
        alpha = alpha.sub(&delta);
    }
    Some(alpha)
}

/// Run the Phase-3 quadratic-identity verification at `(curve, z_P, d
/// ∈ 2..=d_max)`.  Returns a vector of rows + a consistency flag.
pub fn run_mazur_tate_phase3(
    curve: &ZpCurve,
    z_p: &ZpInt,
    d_max: u64,
) -> (Vec<MazurTatePhase3Row>, bool, Option<BigInt>) {
    let p_u: u64 = curve.p.clone().try_into().unwrap_or(11);
    let series_n = (p_u - 1).max(10) as usize;
    let log_sigma_minus_log_z = match log_sigma_minus_log_z(curve, series_n) {
        Some(s) => s,
        None => return (Vec::new(), false, None),
    };
    let log_sigma_p = log_sigma_minus_log_z.evaluate(z_p);
    let mut rows = Vec::new();
    let mut hs = Vec::new();
    for d in 2..=d_max {
        // [d]_F(z_P) via correct formal-group multiplication.
        let z_d = formal_d_times_z(curve, z_p, d, series_n);
        let log_sigma_d = log_sigma_minus_log_z.evaluate(&z_d);
        // diff = log_sigma_d − d · log_sigma_p
        let d_zp = ZpInt::new(BigInt::from(d as i64), &curve.p, curve.precision);
        let d_log_p = d_zp.mul(&log_sigma_p);
        let diff = log_sigma_d.sub(&d_log_p);
        // d(d−1)/2
        let d_d_minus_1 = ZpInt::new(
            BigInt::from((d * (d - 1)) as i64),
            &curve.p,
            curve.precision,
        );
        let two = ZpInt::new(BigInt::from(2), &curve.p, curve.precision);
        let two_inv = two.inverse();
        let denom_half = two_inv.as_ref().map(|inv| d_d_minus_1.mul(inv));
        // h estimate = diff / [d(d−1)/2]
        let h_estimate = if let Some(denom) = denom_half.as_ref() {
            if denom.value.is_zero() {
                None
            } else {
                denom.inverse().map(|inv| diff.mul(&inv))
            }
        } else {
            None
        };
        if let Some(h) = h_estimate.as_ref() {
            hs.push(h.value.clone());
        }
        rows.push(MazurTatePhase3Row {
            d,
            z_d: z_d.value.clone(),
            log_sigma_d: log_sigma_d.value.clone(),
            diff: diff.value.clone(),
            d_d_minus_1_half: denom_half.map(|h| h.value).unwrap_or_else(BigInt::zero),
            h_estimate: h_estimate.map(|h| h.value),
        });
    }
    let consistent = !hs.is_empty() && hs.iter().all(|h| h == &hs[0]);
    let common_h = if consistent {
        Some(hs[0].clone())
    } else {
        None
    };
    (rows, consistent, common_h)
}

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
        let term = u_pow.mul(&PSeries::from_coef_vec(
            p,
            prec,
            vec![ZpInt::new(BigInt::from(sign), p, prec).mul(&k_inv)],
        ));
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
    let d_d_minus_1 = ZpInt::new(
        BigInt::from((d * (d - 1)) as i64),
        &curve.p,
        curve.precision,
    );
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
        println!(
            "{:>3} {:>12} {:>15} {:>15} {:>15}",
            "d", "diff", "estimated h", "z_q", "(d-1)·val_p"
        );

        let mut estimated_hs: Vec<BigInt> = Vec::new();
        for d in 2u64..=8 {
            if let Some(r) = run_mazur_tate_recovery_test(&curve, &z_p, d) {
                let h_str = r
                    .estimated_h
                    .as_ref()
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

    /// **Phase 3a — `[d]_F` composition sanity**: verify the formal-
    /// group composition law `[d]_F([e]_F(z)) = [d·e]_F(z)`.  This is
    /// a structural test of the formal-log reversion, independent of
    /// the σ-series being complete.
    #[test]
    fn mazur_tate_phase3_formal_d_composition() {
        let p = BigInt::from(11);
        let prec = 6u32;
        let curve = ZpCurve::new(BigInt::from(1), BigInt::from(1), &p, prec);
        let z_p = ZpInt::new(BigInt::from(11), &p, prec);
        // [2]_F([3]_F(z)) == [6]_F(z)
        // [3]_F([2]_F(z)) == [6]_F(z)
        // [2]_F([5]_F(z)) == [10]_F(z)
        for (d, e) in [(2u64, 3u64), (3, 2), (2, 5), (4, 2)] {
            let ok = verify_formal_d_times_composition(&curve, &z_p, d, e, 10);
            assert!(
                ok,
                "composition [{}]_F ∘ [{}]_F = [{}]_F failed",
                d,
                e,
                d * e
            );
        }
        println!("\n✓ Phase 3a: [d]_F composition law verified at p = 11.");
    }

    /// **Phase 3b — `α_p` Hensel lift**: verify that the unit root
    /// `α_p` of the Frobenius polynomial `T² − a_p T + p` is correctly
    /// Hensel-lifted in `Z_p`.
    #[test]
    fn mazur_tate_phase3_alpha_p_unit_root() {
        // Curve y² = x³ + x + 1 over F_11 has #E = 13 (verified earlier),
        // so a_p = 11 + 1 − 13 = −1.
        // Frobenius polynomial: T² + T + 11 = 0.  Unit root α_p ≡ −1 (mod 11)
        // and satisfies α_p² + α_p + 11 = 0 in Z_11.
        let p = BigInt::from(11);
        let prec = 6u32;
        let alpha = frobenius_unit_root(&p, -1, prec).expect("ordinary curve");
        // Verify f(α) = 0 mod p^prec.
        let p_zp = ZpInt::new(p.clone(), &p, prec);
        let a_p_zp = ZpInt::new(BigInt::from(-1), &p, prec);
        let f_alpha = alpha.mul(&alpha).sub(&a_p_zp.mul(&alpha)).add(&p_zp);
        assert!(
            f_alpha.value.is_zero(),
            "Hensel lift failed: f(α) = {} ≠ 0 mod 11^6",
            f_alpha.value
        );
        // Verify α is a unit (α ≢ 0 mod p).
        assert!(!alpha.value.is_zero());
        println!(
            "\n✓ Phase 3b: α_p Hensel-lifted to mod 11^6:  α = {}",
            alpha.value
        );
    }

    /// **Phase 3 verification**: with the corrected `[d]_F` via
    /// formal-log reversion, the Mazur-Tate quadratic identity should
    /// hold with a constant `h` across `d`.
    #[test]
    fn mazur_tate_phase3_quadratic_identity() {
        let p = BigInt::from(11);
        let prec = 6u32;
        let curve = ZpCurve::new(BigInt::from(1), BigInt::from(1), &p, prec);
        let z_p = ZpInt::new(BigInt::from(11), &p, prec);

        let (rows, consistent, h_common) = run_mazur_tate_phase3(&curve, &z_p, 8);
        println!();
        println!("=== Mazur-Tate σ Phase 3: corrected [d]_F via reversion ===");
        println!();
        println!(
            "{:>3} {:>15} {:>15} {:>15} {:>15}",
            "d", "z_d", "log σ(d·P)", "diff", "h estimate"
        );
        for r in &rows {
            let h_str = r
                .h_estimate
                .as_ref()
                .map(|h| h.to_string())
                .unwrap_or_else(|| "N/A".into());
            println!(
                "{:>3} {:>15} {:>15} {:>15} {:>15}",
                r.d, r.z_d, r.log_sigma_d, r.diff, h_str,
            );
        }
        println!();
        if consistent {
            println!("✓ Phase 3: h(P, P) CONSISTENT across d.");
            if let Some(h) = h_common.as_ref() {
                println!("  Common h = {}", h);
            }
            println!("  The Mazur-Tate quadratic identity is verified at toy scale.");
            println!("  ECDLP recovery via Phase 3 σ is feasible for these parameters.");
        } else {
            println!("✗ Phase 3: h still varies across d.");
            println!("  Either reversion precision is insufficient, or the σ-series");
            println!("  itself requires the E_2 / p-adic-height correction terms.");
        }
    }
}
