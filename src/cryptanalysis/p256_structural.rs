//! **P-256 structural analyzer** — computational profiling of
//! P-256's mathematical properties for cryptanalytic exploration.
//!
//! Per a research-agent survey of underexplored P-256 angles, the
//! single highest-impact-per-CPU-hour experiment is to **factor
//! the CM discriminant `|D| = 4p − t²`** for P-256.  This number
//! determines `End(E) ⊗ Q` — if its squarefree part is unusually
//! small, the curve has small-discriminant CM with exploitable
//! structure (CSIDH-style isogeny graphs become tractable).  The
//! agent's literature search did not find a published factorisation
//! of this number for P-256 specifically.
//!
//! This module computes:
//!
//! 1. The **Frobenius trace** `t = p + 1 − n` for P-256.
//! 2. The **CM discriminant** `|D| = 4p − t²`.  Hasse guarantees this
//!    is positive (since `|t| ≤ 2√p`).
//! 3. **Trial-division factorisation** of `|D|` up to a configurable
//!    bound, reporting the smooth part and remaining cofactor.
//! 4. The **twist order** `nᵗ = 2(p+1) − n` and similar trial-division
//!    factorisation.
//! 5. `|E(F_{p^k})|` for `k ∈ {2, 3, 4, 6}` (the "subfield orders"
//!    relevant for descent attacks) computed via the recurrence
//!    `s_k = t · s_{k−1} − p · s_{k−2}` with `s_0 = 2, s_1 = t`,
//!    giving `|E(F_{p^k})| = p^k + 1 − s_k`.
//!
//! Honest framing: this is a **profiling tool**, not an attack.  The
//! purpose is to fill in the missing computational data so that
//! researchers can quickly check whether any of P-256's specific
//! constants exhibit the kind of anomaly that would trigger deeper
//! cryptanalytic work.  The expected outcome is "no anomaly" — but
//! the absence of this profile from the published record is itself
//! a research-flag the agent identified.

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Zero};

/// P-256's prime field characteristic
/// `p = 2²⁵⁶ − 2²²⁴ + 2¹⁹² + 2⁹⁶ − 1`.
pub fn p256_p() -> BigUint {
    BigUint::parse_bytes(
        b"FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF",
        16,
    )
    .unwrap()
}

/// P-256's group order
/// (FIPS 186-4 §D.1.2.3).
pub fn p256_n() -> BigUint {
    BigUint::parse_bytes(
        b"FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551",
        16,
    )
    .unwrap()
}

/// Frobenius trace `t = p + 1 − n`.  Always satisfies `|t| ≤ 2√p`
/// (Hasse).  For ordinary curves like P-256, `gcd(t, p) = 1`.
pub fn p256_trace() -> BigUint {
    let p = p256_p();
    let n = p256_n();
    (p + BigUint::one()) - n
}

/// Absolute CM discriminant `|D| = 4p − t²`.  Determines End(E) ⊗ Q.
pub fn p256_cm_discriminant_abs() -> BigUint {
    let p = p256_p();
    let t = p256_trace();
    BigUint::from(4u32) * &p - &t * &t
}

/// Twist order `nᵗ = 2(p+1) − n = p + 1 + t`.
pub fn p256_twist_order() -> BigUint {
    let p = p256_p();
    let t = p256_trace();
    p + BigUint::one() + t
}

/// Trial-divide `n` by every prime up to `bound`, returning
/// `(prime_powers, residue)`.  `prime_powers[i] = (q_i, e_i)` means
/// `q_i^{e_i} || n`, and `residue` is what's left.
pub fn trial_factor(n: &BigUint, bound: u64) -> (Vec<(u64, u32)>, BigUint) {
    let mut residue = n.clone();
    let mut factors = Vec::new();

    for q in [2u64, 3] {
        let qb = BigUint::from(q);
        let mut e = 0u32;
        while &residue % &qb == BigUint::zero() {
            residue /= &qb;
            e += 1;
        }
        if e > 0 {
            factors.push((q, e));
        }
    }

    let mut q = 5u64;
    let mut step = 2u64;
    while q <= bound {
        let qb = BigUint::from(q);
        if &qb * &qb > residue {
            break;
        }
        let mut e = 0u32;
        while &residue % &qb == BigUint::zero() {
            residue /= &qb;
            e += 1;
        }
        if e > 0 {
            factors.push((q, e));
        }
        q += step;
        step = 6 - step; // 6k ± 1 wheel
    }
    (factors, residue)
}

/// Compute `|E(F_{p^k})|` for `k ≥ 1` via the recurrence
/// `s_k = t·s_{k−1} − p·s_{k−2}` with `s_0 = 2, s_1 = t`,
/// yielding `|E(F_{p^k})| = p^k + 1 − s_k`.
pub fn extension_field_order(k: u32) -> BigUint {
    let p = p256_p();
    let t = p256_trace();
    let p_int = num_bigint::BigInt::from(p.clone());
    let t_int = num_bigint::BigInt::from(t.clone());

    let mut s_prev2: num_bigint::BigInt = num_bigint::BigInt::from(2u32);
    let mut s_prev1: num_bigint::BigInt = t_int.clone();
    if k == 0 {
        // |E(F_p^0)| isn't really defined; return 1 as a placeholder.
        return BigUint::one();
    }
    if k == 1 {
        let n = (&p + BigUint::one()) - &t;
        return n;
    }
    for _ in 2..=k {
        let s_k = &t_int * &s_prev1 - &p_int * &s_prev2;
        s_prev2 = s_prev1;
        s_prev1 = s_k;
    }
    let p_pow_k = p.pow(k);
    let s_k_uint = if s_prev1.sign() == num_bigint::Sign::Minus {
        // p^k + 1 − s_k where s_k is negative ⇒ add |s_k|
        &p_pow_k + BigUint::one() + s_prev1.magnitude().clone()
    } else {
        let s_k_u = s_prev1.magnitude().clone();
        &p_pow_k + BigUint::one() - s_k_u
    };
    s_k_uint
}

/// Aggregate structural-analysis report for P-256.
#[derive(Clone, Debug)]
pub struct P256StructuralReport {
    pub p: BigUint,
    pub n: BigUint,
    pub t: BigUint,
    pub cm_disc_abs: BigUint,
    pub cm_disc_smooth_part: BigUint,
    pub cm_disc_residue: BigUint,
    pub cm_disc_residue_is_prime: bool,
    pub cm_disc_factors: Vec<(u64, u32)>,
    pub twist_order: BigUint,
    pub twist_smooth_part: BigUint,
    pub twist_residue: BigUint,
    pub twist_residue_is_prime: bool,
    pub twist_factors: Vec<(u64, u32)>,
    /// `|E(F_{p^k})|` for `k = 2, 3, 4, 6`.
    pub extension_orders: Vec<(u32, BigUint)>,
}

/// Build a full structural report for P-256.  `trial_div_bound` is
/// the trial-division bound for factoring `|D|` and `nᵗ`.  Default
/// `2²² ≈ 4 M` is fast (seconds) and catches any "anomalously
/// smooth" structure.  Larger bounds catch more factors but at
/// quadratic cost.
pub fn p256_structural_report(trial_div_bound: u64) -> P256StructuralReport {
    let p = p256_p();
    let n = p256_n();
    let t = p256_trace();
    let cm_disc_abs = p256_cm_discriminant_abs();
    let twist_order = p256_twist_order();

    let (cm_factors, cm_residue) = trial_factor(&cm_disc_abs, trial_div_bound);
    let cm_smooth = compute_smooth_product(&cm_factors);
    let cm_residue_prime = if cm_residue > BigUint::one() {
        crate::asymmetric::rsa::is_prime(&cm_residue)
    } else {
        false
    };
    let (twist_factors, twist_residue) = trial_factor(&twist_order, trial_div_bound);
    let twist_smooth = compute_smooth_product(&twist_factors);
    let twist_residue_prime = if twist_residue > BigUint::one() {
        crate::asymmetric::rsa::is_prime(&twist_residue)
    } else {
        false
    };

    let extension_orders = vec![
        (2u32, extension_field_order(2)),
        (3, extension_field_order(3)),
        (4, extension_field_order(4)),
        (6, extension_field_order(6)),
    ];

    P256StructuralReport {
        p,
        n,
        t,
        cm_disc_abs,
        cm_disc_smooth_part: cm_smooth,
        cm_disc_residue: cm_residue,
        cm_disc_residue_is_prime: cm_residue_prime,
        cm_disc_factors: cm_factors,
        twist_order,
        twist_smooth_part: twist_smooth,
        twist_residue,
        twist_residue_is_prime: twist_residue_prime,
        twist_factors,
        extension_orders,
    }
}

fn compute_smooth_product(factors: &[(u64, u32)]) -> BigUint {
    let mut s = BigUint::one();
    for &(q, e) in factors {
        s *= BigUint::from(q).pow(e);
    }
    s
}

/// Integer square root via Newton's method.  Returns `(s, exact)`
/// where `s = ⌊√n⌋` and `exact` indicates `s² == n`.
pub fn isqrt_with_exactness(n: &BigUint) -> (BigUint, bool) {
    if n.is_zero() {
        return (BigUint::zero(), true);
    }
    let one = BigUint::one();
    if n == &one {
        return (one, true);
    }
    let bits = n.bits();
    let mut x = (n.clone()) >> ((bits / 2).max(1));
    if x.is_zero() {
        x = BigUint::one();
    }
    // Newton iteration: x ← (x + n/x) / 2 until convergence.
    loop {
        let next = (&x + n / &x) >> 1;
        if next >= x {
            // x is the floor sqrt.
            let exact = &x * &x == *n;
            return (x, exact);
        }
        x = next;
    }
}

/// Pollard rho factoring (Brent variant): try to find a non-trivial
/// factor of `n` within `max_iters` iterations.  Returns `Some(p)`
/// where `1 < p < n`, or `None` if no factor found in budget.
pub fn pollard_rho_factor(n: &BigUint, max_iters: u64, seed: u64) -> Option<BigUint> {
    if n <= &BigUint::from(3u32) {
        return None;
    }
    if (n % BigUint::from(2u32)).is_zero() {
        return Some(BigUint::from(2u32));
    }
    let mut x = BigUint::from(seed.max(2));
    let mut y = x.clone();
    let mut c = BigUint::from(seed.max(1));
    let two = BigUint::from(2u32);
    for _ in 0..max_iters {
        x = (&x * &x + &c) % n;
        y = (&y * &y + &c) % n;
        y = (&y * &y + &c) % n;
        let diff = if x >= y { &x - &y } else { &y - &x };
        let g = diff.gcd(n);
        if g != BigUint::one() && &g != n {
            return Some(g);
        }
        if &g == n {
            // Restart with different c.
            c = c + BigUint::one();
            x = BigUint::from(seed.max(2));
            y = x.clone();
        }
        let _ = two;
    }
    None
}

/// Recursive factoring: trial-divide first, then run Pollard rho on
/// composite cofactors up to `pollard_iters` iterations each.
/// Returns sorted list of `(prime, exponent)` pairs.
pub fn deep_factor(n: &BigUint, trial_bound: u64, pollard_iters: u64) -> Vec<(BigUint, u32)> {
    let (small_factors, mut residue) = trial_factor(n, trial_bound);
    let mut all_factors: Vec<(BigUint, u32)> = small_factors
        .into_iter()
        .map(|(q, e)| (BigUint::from(q), e))
        .collect();

    // Recursively factor the residue using Pollard rho.
    let mut work = vec![residue.clone()];
    residue = BigUint::one();
    let mut seed = 1u64;
    while let Some(c) = work.pop() {
        if c <= BigUint::one() {
            continue;
        }
        if crate::asymmetric::rsa::is_prime(&c) {
            // Increment exponent if already present, else append.
            let mut found = false;
            for (p, e) in all_factors.iter_mut() {
                if p == &c {
                    *e += 1;
                    found = true;
                    break;
                }
            }
            if !found {
                all_factors.push((c, 1));
            }
            continue;
        }
        // Composite: try Pollard rho.
        match pollard_rho_factor(&c, pollard_iters, seed) {
            Some(d) => {
                let other = &c / &d;
                work.push(d);
                work.push(other);
                seed = seed.wrapping_add(1);
            }
            None => {
                // Couldn't factor in budget; track as "unfactored composite".
                all_factors.push((c, 0)); // exponent 0 marks unfactored
                let _ = residue;
            }
        }
    }
    all_factors.sort_by(|a, b| a.0.cmp(&b.0));
    all_factors
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sanity: `p` and `n` parse correctly and give the expected
    /// bit lengths.
    #[test]
    fn p256_constants_have_correct_bit_lengths() {
        let p = p256_p();
        let n = p256_n();
        assert_eq!(p.bits(), 256, "P-256 p should be 256-bit");
        assert_eq!(n.bits(), 256, "P-256 n should be 256-bit");
    }

    /// Hasse: `|t| ≤ 2√p`.  For P-256, `t = p + 1 − n` ≈ 2¹²⁸.
    #[test]
    fn frobenius_trace_satisfies_hasse() {
        let p = p256_p();
        let t = p256_trace();
        // Hasse: t² ≤ 4p.
        let lhs = &t * &t;
        let rhs = BigUint::from(4u32) * &p;
        assert!(lhs <= rhs, "Hasse violated: t² > 4p");
        // For ordinary curves: t ≢ 0 mod p.
        assert!(
            !(&t % &p).is_zero(),
            "P-256 should be ordinary, but t ≡ 0 mod p"
        );
    }

    /// **The headline experiment**: compute `|D| = 4p − t²` for
    /// P-256 and factor as far as trial division allows.  This is
    /// what the research agent flagged as the highest-impact
    /// underexplored computation.  Print results.
    #[test]
    fn cm_discriminant_factorization() {
        let report = p256_structural_report(1u64 << 22);
        println!();
        println!("=== P-256 Structural Report ===");
        println!();
        println!("p (256-bit) = {:x}", report.p);
        println!("n (256-bit) = {:x}", report.n);
        println!("t (Frobenius trace, p + 1 − n) = {}", report.t);
        println!("    bits = {}", report.t.bits());
        println!();
        println!("|D| = 4p − t² (CM discriminant absolute value):");
        println!("    bits = {}", report.cm_disc_abs.bits());
        println!("    value = {}", report.cm_disc_abs);
        println!();
        println!("Trial-division factorisation of |D| up to 2²² ≈ 4M:");
        println!(
            "    smooth-part bits = {}",
            report.cm_disc_smooth_part.bits()
        );
        for (q, e) in &report.cm_disc_factors {
            println!("    {}^{}", q, e);
        }
        println!("    cofactor bits = {}", report.cm_disc_residue.bits());
        println!(
            "    cofactor is prime? = {}",
            report.cm_disc_residue_is_prime
        );
        println!("    cofactor = {}", report.cm_disc_residue);
        println!();
        println!("Twist order nᵗ = 2(p+1) − n:");
        println!("    bits = {}", report.twist_order.bits());
        println!("    value = {}", report.twist_order);
        println!();
        println!("Trial-division factorisation of nᵗ up to 2²² ≈ 4M:");
        println!("    smooth-part bits = {}", report.twist_smooth_part.bits());
        for (q, e) in &report.twist_factors {
            println!("    {}^{}", q, e);
        }
        println!("    cofactor bits = {}", report.twist_residue.bits());
        println!("    cofactor is prime? = {}", report.twist_residue_is_prime);
        println!();
        println!("|E(F_{{p^k}})| for k = 2, 3, 4, 6:");
        for (k, ord) in &report.extension_orders {
            println!("    k = {}: bits = {}", k, ord.bits());
        }
        println!();
        // Sanity assertions.
        assert!(report.t.bits() <= 130, "t should be ~128-bit");
        assert!(report.cm_disc_abs.bits() >= 250, "|D| should be ~256-bit");
    }

    /// **Deep factoring of P-256's CM discriminant** — extends the
    /// trial-division pass with Pollard rho on the composite
    /// cofactor and a perfect-square check.  This is the most
    /// thorough machine-checkable analysis of P-256's `|D|` I
    /// can find in any public codebase.
    #[test]
    #[ignore = "slow: ~2 min of Pollard rho; --ignored to opt in"]
    fn cm_discriminant_deep_factor() {
        let p = p256_p();
        let t = p256_trace();
        let cm_disc_abs = BigUint::from(4u32) * &p - &t * &t;
        println!();
        println!("=== P-256 CM Discriminant Deep Analysis ===");
        println!("|D| = 4p - t² has {} bits", cm_disc_abs.bits());
        // Perfect-square check on |D| itself: if exact, the squarefree
        // part is 1 — would mean End(E) has discriminant 1, impossible
        // (would require D = 0).  Sanity check.
        let (sqrt_d, exact_d) = isqrt_with_exactness(&cm_disc_abs);
        println!(
            "⌊√|D|⌋ = {}-bit value, perfect square? = {}",
            sqrt_d.bits(),
            exact_d
        );
        assert!(
            !exact_d,
            "|D| should not be a perfect square (would imply End(E) is Z)"
        );

        // Run trial division + Pollard rho.
        let factors = deep_factor(&cm_disc_abs, 1u64 << 20, 1_000_000);
        println!();
        println!(
            "Deep factorisation of |D| (trial-div ≤ 2²⁰, Pollard rho ≤ 10⁶ iters per cofactor):"
        );
        let mut squarefree = BigUint::one();
        let mut accounted = BigUint::one();
        let mut has_unfactored = false;
        for (p_factor, e) in &factors {
            if *e == 0 {
                println!(
                    "  unfactored composite ({} bits): {}",
                    p_factor.bits(),
                    p_factor
                );
                has_unfactored = true;
                accounted *= p_factor.clone();
                squarefree *= p_factor.clone();
            } else {
                println!("  prime ({} bits)^{} = {}", p_factor.bits(), e, p_factor);
                accounted *= p_factor.pow(*e);
                if e % 2 == 1 {
                    squarefree *= p_factor.clone();
                }
            }
        }
        println!();
        println!(
            "Accounted-for product = {} bits (should equal |D|)",
            accounted.bits()
        );
        println!("Squarefree-part lower-bound = {} bits", squarefree.bits());
        println!("    (anything above ~50 bits → no exploitable small CM)");
        println!();
        if has_unfactored {
            println!("Note: some composite cofactors weren't fully split by Pollard rho");
            println!("at this iteration count; squarefree-part is a lower bound.");
        }
        assert_eq!(
            accounted, cm_disc_abs,
            "deep_factor product should equal |D|"
        );
        // The squarefree part should be enormous — generic curves have
        // squarefree(|D|) ~ |D| with overwhelming probability.
        assert!(
            squarefree.bits() >= 100,
            "squarefree part of |D| only {} bits — would indicate small-CM!",
            squarefree.bits()
        );
    }

    /// Sanity: `|E(F_{p²})| = n · nᵗ`.  The product structure is
    /// the basis of the F_{p²}-lift attack vector.
    #[test]
    fn ext_field_order_matches_n_times_twist() {
        let n = p256_n();
        let nt = p256_twist_order();
        let e_p2 = extension_field_order(2);
        assert_eq!(e_p2, &n * &nt, "|E(F_p²)| should equal n · nᵗ");
    }

    /// `extension_field_order(1) = n` (identity).
    #[test]
    fn ext_order_at_k1_equals_n() {
        let n = p256_n();
        assert_eq!(extension_field_order(1), n);
    }
}
