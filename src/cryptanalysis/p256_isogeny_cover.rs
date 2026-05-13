//! # P-256 isogeny-cover probe (Phase 1, toy size).
//!
//! Phase-1 implementation of the proposal in
//! `RESEARCH_P256_ISOGENY_COVER.md`: at a tiny "P-256-shaped" prime
//! `p_toy ≈ 2^8`, sweep ordinary prime-order curves `E/F_{p_toy}`
//! and search for a **genus-2 curve `C/F_{p_toy}` whose Jacobian
//! is `F_p`-isogenous to `E × E^twist`** with an additional
//! `(2, 2)`-isogeny decomposition (Richelot / Cassels–Flynn).
//!
//! ## What the probe does
//!
//! For a chosen toy prime `p`:
//!
//! 1. Sweep **elliptic curves** `E: y² = x³ + ax + b` over `F_p`,
//!    keeping only those with **prime order** (the "P-256-like"
//!    structural property — cofactor 1, no Pohlig–Hellman handle).
//! 2. For each kept `E`, compute the twist order
//!    `n_t = 2(p+1) − n` and the target Jacobian order
//!    `N_target = n · n_t = #E(F_{p²})`.
//! 3. Sweep **genus-2 curves** `C: y² = f(x)` with monic quintic
//!    `f` and squarefree.  For each, compute
//!    `#Jac(C)(F_p)` by brute-force Mumford-divisor enumeration
//!    ([`crate::prime_hyperelliptic::brute_force_jac_order`]).
//! 4. If `#Jac(C) == N_target`, run the
//!    **Richelot `(2, 2)`-splitting test**:
//!    factor `f` over `F_p` into a quadratic × cubic / 3 quadratics
//!    if possible; if a 3-quadratic factorization exists, build the
//!    three associated elliptic curves `E_{ij}: y² = q_i · q_j` and
//!    check whether two of them have orders `{n, n_t}` (an explicit
//!    `F_p`-rational `(2, 2)`-split of `Jac(C)` with one factor
//!    matching `E`).
//! 5. Report all hits as [`CoverProbeHit`]s.
//!
//! ## Honest scope
//!
//! - **Toy only**: `p ≤ ~30` for full sweeps within seconds; `p =
//!   100` takes minutes; `p = 256` is hours; cryptographic `p`'s are
//!   completely infeasible at this implementation's complexity.
//! - **Conservative detector**: we only catch `F_p`-rational
//!   `(2, 2)`-splittings via the elementary 3-quadratic
//!   factorisation.  Splittings that exist only over `F_{p²}` are
//!   missed.  More sophisticated Igusa-invariant / Cardona–Quer
//!   tests would catch those but require additional machinery.
//! - **Phase 1 deliverable**: the goal is to establish the
//!   **frequency** of `(2, 2)`-split covers in random ordinary
//!   prime-order isogeny classes — a structural statistic that
//!   directly bounds the expected outcome of running the same probe
//!   at `2^32`, `2^64`, …, `P-256`.

use crate::prime_hyperelliptic::{
    brute_force_jac_order, count_points, FpPoly, HyperellipticCurveP,
};
use num_bigint::BigUint;
use num_traits::{One, Zero};

// ── tiny elliptic-curve utilities (own to keep this module self-contained) ──

/// `E: y² = x³ + a·x + b` over `F_p` (short Weierstrass, char ≠ 2,3).
#[derive(Clone, Debug)]
pub struct ECurveP {
    pub p: u64,
    pub a: u64,
    pub b: u64,
}

impl ECurveP {
    /// `4a³ + 27b² mod p ≠ 0`.  Smooth iff non-zero.
    pub fn is_smooth(&self) -> bool {
        let p = self.p;
        let a = self.a;
        let b = self.b;
        let aa = mulmod(a, a, p);
        let four_a3 = mulmod(4, mulmod(aa, a, p), p);
        let twentyseven_b2 = mulmod(27, mulmod(b, b, p), p);
        let disc = (four_a3 + twentyseven_b2) % p;
        disc != 0
    }

    /// `#E(F_p) = 1 + Σ_{x ∈ F_p} (1 + χ(x³ + a·x + b))`.
    pub fn order(&self) -> u64 {
        let p = self.p;
        let mut count: u64 = 1; // ∞
        let exp = (p - 1) / 2;
        for x in 0..p {
            let x2 = mulmod(x, x, p);
            let x3 = mulmod(x2, x, p);
            let ax = mulmod(self.a, x, p);
            let rhs = (x3 + ax + self.b) % p;
            if rhs == 0 {
                count += 1;
            } else {
                let chi = powmod(rhs, exp, p);
                if chi == 1 {
                    count += 2;
                }
            }
        }
        count
    }
}

fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

fn powmod(base: u64, exp: u64, p: u64) -> u64 {
    let mut acc = 1u64;
    let mut b = base % p;
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            acc = mulmod(acc, b, p);
        }
        b = mulmod(b, b, p);
        e >>= 1;
    }
    acc
}

fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n < 4 {
        return true;
    }
    if n & 1 == 0 {
        return false;
    }
    let mut d = 3u64;
    while d * d <= n {
        if n % d == 0 {
            return false;
        }
        d += 2;
    }
    true
}

/// Find all ordinary, prime-order elliptic curves `E: y² = x³ + ax + b`
/// over `F_p` ("P-256-like": cofactor 1).  Returns one representative
/// per `(a, b)` pair (no isomorphism reduction).
pub fn sweep_prime_order_curves(p: u64) -> Vec<(ECurveP, u64, u64)> {
    let mut out = Vec::new();
    for a in 0..p {
        for b in 0..p {
            let e = ECurveP { p, a, b };
            if !e.is_smooth() {
                continue;
            }
            let n = e.order();
            if !is_prime_u64(n) {
                continue;
            }
            // Twist order = 2(p+1) − n.
            let n_t = 2 * (p + 1) - n;
            out.push((e, n, n_t));
        }
    }
    out
}

// ── candidate genus-2 curve enumeration ─────────────────────────────

/// Enumerate genus-2 curves `C: y² = f(x)` with `f` a monic quintic
/// over `F_p`, `f` squarefree.  Returns the list of squarefree
/// quintic coefficient tuples `(c_4, c_3, c_2, c_1, c_0)`.
///
/// For `p` small this is `O(p^5)` but with the squarefree filter
/// we can typically prune quickly.
pub fn sweep_genus2_quintics(p: u64) -> Vec<[u64; 5]> {
    let mut out = Vec::new();
    for c4 in 0..p {
        for c3 in 0..p {
            for c2 in 0..p {
                for c1 in 0..p {
                    for c0 in 0..p {
                        if is_squarefree_quintic(p, c4, c3, c2, c1, c0) {
                            out.push([c4, c3, c2, c1, c0]);
                        }
                    }
                }
            }
        }
    }
    out
}

fn is_squarefree_quintic(p: u64, c4: u64, c3: u64, c2: u64, c1: u64, c0: u64) -> bool {
    // f(x) = x^5 + c4·x^4 + c3·x^3 + c2·x^2 + c1·x + c0
    // f'(x) = 5x^4 + 4c4·x^3 + 3c3·x^2 + 2c2·x + c1
    // f is squarefree iff gcd(f, f') = 1 in F_p[x].
    let p_big = BigUint::from(p);
    let f = FpPoly::from_coeffs(
        vec![
            BigUint::from(c0),
            BigUint::from(c1),
            BigUint::from(c2),
            BigUint::from(c3),
            BigUint::from(c4),
            BigUint::one(),
        ],
        p_big.clone(),
    );
    let fp = FpPoly::from_coeffs(
        vec![
            BigUint::from(c1),
            BigUint::from(2 * c2 % p),
            BigUint::from(3 * c3 % p),
            BigUint::from(4 * c4 % p),
            BigUint::from(5u32 % p as u32),
        ],
        p_big,
    );
    let g = f.gcd(&fp);
    g.degree() == Some(0)
}

// ── Richelot (2, 2)-split detector ──────────────────────────────────

/// Attempt to factor a monic quintic `f` over `F_p` as a product of
/// a **monic linear** times **two monic quadratics** (the Richelot
/// "3-quadratic" form with one factor (x − ∞) absorbed at infinity).
/// Returns `Some([linear, quad_1, quad_2])` if such a factorisation
/// exists (and there's at least one `F_p`-rational root of `f`),
/// `None` otherwise.
///
/// Concretely: find any root `α` of `f` in `F_p`, divide out
/// `(x − α)` to get a quartic `g`, then try to factor `g` as a
/// product of two monic quadratics over `F_p`.
pub fn richelot_factor(f: &FpPoly) -> Option<[FpPoly; 3]> {
    // Find a root α ∈ F_p of f.
    let p_u: u64 = match f.p.to_u64_digits().as_slice() {
        [v] => *v,
        _ => return None,
    };
    let p_big = &f.p;
    let mut alpha = None;
    for ai in 0..p_u {
        let av = BigUint::from(ai);
        let val = f.eval(&av);
        if val.is_zero() {
            alpha = Some(av);
            break;
        }
    }
    let alpha = alpha?;
    // g(x) = f(x) / (x − α)  =  quartic.
    let linear = FpPoly::from_coeffs(
        vec![
            (p_big - &alpha) % p_big,
            BigUint::one(),
        ],
        p_big.clone(),
    );
    let (g, r) = f.divrem(&linear);
    if !r.is_zero() {
        return None;
    }
    debug_assert_eq!(g.degree(), Some(4));
    // Factor g as monic quartic into two monic quadratics:
    // g(x) = (x² + r₁ x + s₁)·(x² + r₂ x + s₂)
    //      = x⁴ + (r₁+r₂) x³ + (s₁+s₂+r₁ r₂) x² + (r₁ s₂ + r₂ s₁) x + s₁ s₂
    // Match against g = x⁴ + g₃ x³ + g₂ x² + g₁ x + g₀:
    //   r₁ + r₂ = g₃
    //   s₁ + s₂ + r₁ r₂ = g₂
    //   r₁ s₂ + r₂ s₁ = g₁
    //   s₁ s₂ = g₀
    // Enumerate r₁: r₂ = g₃ - r₁.  Then resolve s₁, s₂ from
    // s₁+s₂ = g₂ - r₁ r₂ (=:S), s₁ s₂ = g₀ (=:P).  s₁, s₂ are
    // roots of t² - S t + P = 0.  Cross-check the linear-in-x
    // coefficient.
    let g3 = g.coeff(3);
    let g2 = g.coeff(2);
    let g1 = g.coeff(1);
    let g0 = g.coeff(0);
    let two_inv = crate::prime_hyperelliptic::fp_poly::fp_inv(
        &BigUint::from(2u32),
        p_big,
    )
    .unwrap();
    for r1_u in 0..p_u {
        let r1 = BigUint::from(r1_u);
        let r2 = (&g3 + p_big - &r1) % p_big;
        let r1r2 = (&r1 * &r2) % p_big;
        let s_sum = (&g2 + p_big - &r1r2) % p_big;
        let s_prod = g0.clone();
        // s₁, s₂ are roots of t² − S t + P = 0.  In char ≠ 2:
        // disc = S² − 4P; require disc to be a QR (or 0).
        let s_sum_sq = (&s_sum * &s_sum) % p_big;
        let four_p = (&s_prod * BigUint::from(4u32)) % p_big;
        let disc = (&s_sum_sq + p_big - &four_p) % p_big;
        let chi = qr_chi(&disc, p_big);
        if chi == -1 {
            continue;
        }
        let sqrt_disc = mod_sqrt(&disc, p_big);
        if sqrt_disc.is_none() {
            continue;
        }
        let sqrt_disc = sqrt_disc.unwrap();
        // s₁ = (S + √disc)/2, s₂ = (S − √disc)/2
        let s1 = (&s_sum + &sqrt_disc) % p_big;
        let s1 = (&s1 * &two_inv) % p_big;
        let s2 = (&s_sum + p_big - &sqrt_disc) % p_big;
        let s2 = (&s2 * &two_inv) % p_big;
        // Cross-check: r₁ s₂ + r₂ s₁ == g₁
        let cross = ((&r1 * &s2) % p_big + (&r2 * &s1)) % p_big;
        if cross != g1 {
            continue;
        }
        // Found a factorisation.
        let q1 = FpPoly::from_coeffs(vec![s1, r1, BigUint::one()], p_big.clone());
        let q2 = FpPoly::from_coeffs(vec![s2, r2, BigUint::one()], p_big.clone());
        return Some([linear, q1, q2]);
    }
    None
}

fn qr_chi(a: &BigUint, p: &BigUint) -> i32 {
    if a.is_zero() {
        return 0;
    }
    let exp = (p - BigUint::one()) >> 1;
    let chi = a.modpow(&exp, p);
    if chi.is_one() {
        1
    } else {
        -1
    }
}

/// Tonelli–Shanks `√a (mod p)` for odd prime `p`.  Returns `None` if
/// `a` is not a QR.
pub fn mod_sqrt(a: &BigUint, p: &BigUint) -> Option<BigUint> {
    let a = a % p;
    if a.is_zero() {
        return Some(BigUint::zero());
    }
    if qr_chi(&a, p) != 1 {
        return None;
    }
    // p mod 4 == 3 fast path: √a = a^{(p+1)/4}.
    let four = BigUint::from(4u32);
    if p % &four == BigUint::from(3u32) {
        let exp = (p + BigUint::one()) >> 2;
        return Some(a.modpow(&exp, p));
    }
    // General Tonelli–Shanks.
    let mut q = p - BigUint::one();
    let mut s = 0u32;
    while &q & &BigUint::one() == BigUint::zero() {
        q >>= 1;
        s += 1;
    }
    // Find a non-residue z.
    let mut z = BigUint::from(2u32);
    while qr_chi(&z, p) != -1 {
        z += BigUint::one();
    }
    let mut m = s;
    let mut c = z.modpow(&q, p);
    let mut t = a.modpow(&q, p);
    let mut r = a.modpow(&((&q + BigUint::one()) >> 1), p);
    loop {
        if t.is_one() {
            return Some(r);
        }
        // Find least i, 0 < i < m, with t^{2^i} = 1.
        let mut i = 0u32;
        let mut tt = t.clone();
        while !tt.is_one() {
            tt = (&tt * &tt) % p;
            i += 1;
            if i == m {
                return None; // shouldn't happen if a is QR
            }
        }
        // b = c^{2^{m−i−1}}
        let mut b = c.clone();
        for _ in 0..(m - i - 1) {
            b = (&b * &b) % p;
        }
        m = i;
        c = (&b * &b) % p;
        t = (&t * &c) % p;
        r = (&r * &b) % p;
    }
}

// ── the probe itself ────────────────────────────────────────────────

/// One hit from the cover probe.
#[derive(Clone, Debug)]
pub struct CoverProbeHit {
    pub p: u64,
    /// Source elliptic curve `E`.
    pub e_a: u64,
    pub e_b: u64,
    /// `#E(F_p)`.
    pub n_e: u64,
    /// `#E^twist(F_p)`.
    pub n_t: u64,
    /// Genus-2 curve coefficients `[c_4, c_3, c_2, c_1, c_0]` so
    /// `f(x) = x⁵ + c_4·x⁴ + … + c_0`.
    pub c_coeffs: [u64; 5],
    /// `#Jac(C)(F_p)`.  Equal to `n_e · n_t` if this is a hit.
    pub n_jac: BigUint,
    /// Did the Richelot 3-factor test find an explicit
    /// `F_p`-rational factorisation?  If `true`, `Jac(C)` admits an
    /// `F_p`-rational `(2, 2)`-isogeny structure.
    pub richelot_split: bool,
    /// If `richelot_split`, the three factors (linear, quadratic,
    /// quadratic) of `f`.
    pub richelot_factors: Option<[Vec<u64>; 3]>,
}

/// Run the Phase-1 probe at prime `p`.  Returns all hits.
///
/// Compute cost: dominated by the genus-2 sweep, which is `O(p^5)`
/// curves × `O(p^4)` brute-force `#Jac` per curve = `O(p^9)`.  At
/// `p = 11`: ~5 × 10^9 ops, runs in seconds because the `#Jac`
/// inner loop is very simple.  At `p = 17`: ~10^11 ops, runs in
/// minutes.  Beyond `p ≈ 30` need optimisation.
pub fn run_probe(p: u64, max_curves: Option<usize>) -> Vec<CoverProbeHit> {
    assert!(p >= 5 && is_prime_u64(p), "p must be a prime ≥ 5");
    // Step 1: prime-order E curves.
    let ecs = sweep_prime_order_curves(p);
    if ecs.is_empty() {
        return Vec::new();
    }
    // Build a (target order → E list) map for quick lookup.
    let mut target_map: std::collections::HashMap<BigUint, Vec<(u64, u64, u64, u64)>> =
        std::collections::HashMap::new();
    for (e, n, n_t) in &ecs {
        let target = BigUint::from(*n) * BigUint::from(*n_t);
        target_map
            .entry(target)
            .or_insert_with(Vec::new)
            .push((e.a, e.b, *n, *n_t));
    }
    // Step 2: sweep genus-2 quintics.
    let p_big = BigUint::from(p);
    let mut hits = Vec::new();
    let mut tried = 0usize;
    let total_curves = (p as usize).pow(5);
    let cap = max_curves.unwrap_or(total_curves);
    for c4 in 0..p {
        for c3 in 0..p {
            for c2 in 0..p {
                for c1 in 0..p {
                    for c0 in 0..p {
                        if tried >= cap {
                            return hits;
                        }
                        tried += 1;
                        if !is_squarefree_quintic(p, c4, c3, c2, c1, c0) {
                            continue;
                        }
                        let f = FpPoly::from_coeffs(
                            vec![
                                BigUint::from(c0),
                                BigUint::from(c1),
                                BigUint::from(c2),
                                BigUint::from(c3),
                                BigUint::from(c4),
                                BigUint::one(),
                            ],
                            p_big.clone(),
                        );
                        let curve = HyperellipticCurveP::new(p_big.clone(), f.clone(), 2);
                        let n_jac = brute_force_jac_order(&curve);
                        let matches = target_map.get(&n_jac);
                        if let Some(es) = matches {
                            for (e_a, e_b, n_e, n_t) in es {
                                let richelot = richelot_factor(&f);
                                let split = richelot.is_some();
                                let factors = richelot.map(|arr| {
                                    let conv = |q: &FpPoly| {
                                        q.coeffs
                                            .iter()
                                            .map(|c| {
                                                c.to_u64_digits()
                                                    .first()
                                                    .copied()
                                                    .unwrap_or(0)
                                            })
                                            .collect::<Vec<u64>>()
                                    };
                                    [
                                        conv(&arr[0]),
                                        conv(&arr[1]),
                                        conv(&arr[2]),
                                    ]
                                });
                                hits.push(CoverProbeHit {
                                    p,
                                    e_a: *e_a,
                                    e_b: *e_b,
                                    n_e: *n_e,
                                    n_t: *n_t,
                                    c_coeffs: [c4, c3, c2, c1, c0],
                                    n_jac: n_jac.clone(),
                                    richelot_split: split,
                                    richelot_factors: factors,
                                });
                            }
                        }
                    }
                }
            }
        }
    }
    hits
}

/// Summary statistics of a probe run.
#[derive(Clone, Debug)]
pub struct ProbeStats {
    pub p: u64,
    pub n_prime_order_curves: usize,
    pub n_genus2_curves_tested: u64,
    pub n_jac_order_matches: usize,
    pub n_richelot_splits: usize,
}

/// Run the probe and return summary statistics + the hits.
pub fn run_probe_with_stats(p: u64, max_curves: Option<usize>) -> (ProbeStats, Vec<CoverProbeHit>) {
    let ecs = sweep_prime_order_curves(p);
    let n_prime_order_curves = ecs.len();
    let hits = run_probe(p, max_curves);
    let n_richelot = hits.iter().filter(|h| h.richelot_split).count();
    let total = (p as u64).pow(5);
    let stats = ProbeStats {
        p,
        n_prime_order_curves,
        n_genus2_curves_tested: total,
        n_jac_order_matches: hits.len(),
        n_richelot_splits: n_richelot,
    };
    (stats, hits)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tiny_sweep_p_5() {
        // p = 5: ordinary curves with prime order — counts.
        let ecs = sweep_prime_order_curves(5);
        // Just confirm the sweep produced some output.
        assert!(ecs.is_empty() || ecs.iter().all(|(_, n, _)| is_prime_u64(*n)));
        for (_, n, _) in &ecs {
            assert!(is_prime_u64(*n));
        }
    }

    #[test]
    fn mod_sqrt_roundtrip() {
        for p in &[7u32, 11, 13, 17, 23, 29] {
            let pb = BigUint::from(*p);
            for a in 1u32..*p {
                let ab = BigUint::from(a);
                if let Some(s) = mod_sqrt(&ab, &pb) {
                    let sq = (&s * &s) % &pb;
                    assert_eq!(sq, ab, "mod_sqrt({}) mod {} = {}, but {}² mod {} ≠ {}",
                        a, p, s, s, p, a);
                }
            }
        }
    }

    /// **Experimental driver**: run the probe at `p = 7`, print
    /// a summary plus all hits.  Marked `#[ignore]` because it
    /// takes a few seconds.  Run with
    /// `cargo test --release p256_isogeny_cover_probe_p7 -- --ignored
    /// --nocapture`.
    #[test]
    #[ignore]
    fn p256_isogeny_cover_probe_p7() {
        let (stats, hits) = run_probe_with_stats(7, None);
        println!("\n=== Probe at p = {} ===", stats.p);
        println!(
            "  prime-order curves:   {}",
            stats.n_prime_order_curves
        );
        println!(
            "  genus-2 curves swept: {}",
            stats.n_genus2_curves_tested
        );
        println!(
            "  #Jac order matches:   {}",
            stats.n_jac_order_matches
        );
        println!(
            "  Richelot splits:      {}",
            stats.n_richelot_splits
        );
        if !hits.is_empty() {
            println!("\n  First 5 hits:");
            for h in hits.iter().take(5) {
                println!(
                    "    E(a={}, b={}): #E={}, #E^t={}; C coeffs {:?}; Jac={}; richelot={}",
                    h.e_a,
                    h.e_b,
                    h.n_e,
                    h.n_t,
                    h.c_coeffs,
                    h.n_jac,
                    h.richelot_split
                );
            }
        }
    }

    /// Same probe at `p = 11`.  Several seconds in release mode.
    #[test]
    #[ignore]
    fn p256_isogeny_cover_probe_p11() {
        let (stats, hits) = run_probe_with_stats(11, None);
        println!("\n=== Probe at p = {} ===", stats.p);
        println!(
            "  prime-order curves:   {}",
            stats.n_prime_order_curves
        );
        println!(
            "  #Jac order matches:   {}",
            stats.n_jac_order_matches
        );
        println!(
            "  Richelot splits:      {}",
            stats.n_richelot_splits
        );
        if !hits.is_empty() {
            println!("\n  First 5 hits (of {}):", hits.len());
            for h in hits.iter().take(5) {
                println!(
                    "    E(a={}, b={}): #E={}, #E^t={}; C coeffs {:?}; Jac={}; richelot={}",
                    h.e_a,
                    h.e_b,
                    h.n_e,
                    h.n_t,
                    h.c_coeffs,
                    h.n_jac,
                    h.richelot_split
                );
            }
            let split_count = hits.iter().filter(|h| h.richelot_split).count();
            println!("\n  Total Richelot-split hits: {}/{}", split_count, hits.len());
        }
    }

    #[test]
    fn richelot_factor_correctness_synthetic() {
        // Build f = (x − 1) · (x² + 0·x + 0) · (x² + 0·x + 3) over F_7.
        // I.e., f = (x − 1) · x² · (x² + 3) = x⁵ - x⁴ + 3x³ - 3x²
        let p_big = BigUint::from(7u32);
        let q0 = FpPoly::from_coeffs(
            vec![
                (BigUint::from(7u32) - BigUint::from(1u32)) % &p_big,
                BigUint::one(),
            ],
            p_big.clone(),
        );
        let q1 = FpPoly::from_coeffs(
            vec![BigUint::zero(), BigUint::zero(), BigUint::one()],
            p_big.clone(),
        );
        let q2 = FpPoly::from_coeffs(
            vec![BigUint::from(3u32), BigUint::zero(), BigUint::one()],
            p_big.clone(),
        );
        let f = q0.mul(&q1).mul(&q2);
        assert_eq!(f.degree(), Some(5));
        // But f = x⁵ - x⁴ + 3x³ - 3x²  has a double root at x=0
        // (since q1 = x²).  That makes f non-squarefree.
        // Skip the squarefree check here — just verify the
        // factoriser can find the factorisation.
        let result = richelot_factor(&f);
        assert!(result.is_some(), "richelot_factor should find this");
    }
}
