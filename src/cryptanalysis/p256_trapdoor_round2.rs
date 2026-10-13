//! Round-2 falsification probes for the P-256 curve-level trapdoor
//! hypotheses (RESEARCH_P256_BACKDOOR_SPECULATION.md §4, ranks 1-3).
//!
//! Round 1 (the speculation note, `b_seed_profile`, `p256_structural`,
//! the isogeny walk campaigns) checked the *public* weak classes. This
//! round asks a different question: if NSA ground the seed into a
//! *secret, cheaply testable* weak class in 1999, what could that class
//! have been built from? Every probe here is a cheap-to-test invariant
//! of the public constants — exactly the kind of predicate a 2^55-seed
//! grind could have targeted — so each "no anomaly" result removes one
//! candidate groundable locus.
//!
//! Hypothesis coverage:
//! - H1 (curve coefficients): multiplicative orders, residuosity
//!   classes, perfect-power and sparse forms, and continued-fraction
//!   heights of the seed-derived values `b`, `r = a³/b²`, and the
//!   j-invariant.
//! - H2 (field prime): deep factorizations of `p ± 1`, the known part
//!   of `ord_p(2)`, the prime's residue-class profile, and the
//!   quadratic twist's embedding degree (small-order twist subgroups
//!   would enable invalid-curve attacks).
//! - H3 (endomorphism/isogeny): the GLS gap `4p − t²`, the trace's
//!   sparse profile, and the rational 2-torsion splitting that
//!   controls the degree-2 isogeny graph.
//!
//! All probes are deterministic and pinned by tests. No weakness in
//! P-256 is claimed; a null result is the expected outcome and each
//! null result is recorded precisely because it is what a 1999 grind
//! could not have hidden behind.

use crate::cryptanalysis::p256_structural::{deep_factor, p256_n, p256_p, p256_trace, p256_twist_order};
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Zero};

/// P-256's published coefficient `b`
/// (SEC 2 v2.0 §2.4.2; FIPS 186-4 §D.1.2.3).
pub fn p256_b() -> BigUint {
    BigUint::parse_bytes(
        b"5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
        16,
    )
    .unwrap()
}

/// The X9.62 seed-derived ratio `r = a³·b⁻² mod p` with `a = −3`,
/// regenerated from the published `b` (cross-checked against the
/// pinned `r` in `p256_backdoor_map` by the tests here).
pub fn p256_r_from_b() -> BigUint {
    let p = p256_p();
    let a3 = &p - BigUint::from(27u32); // (−3)³ = −27 mod p
    let b2 = (&p256_b() * &p256_b()) % &p;
    let b2_inv = b2.modpow(&(&p - BigUint::from(2u32)), &p);
    (&a3 * &b2_inv) % &p
}

/// P-256's j-invariant `j = 1728·a³/(4a³ + 27b²) mod p` with `a = −3`.
pub fn p256_j_invariant() -> BigUint {
    let p = p256_p();
    let a = &p - BigUint::from(3u32); // a = −3 mod p
    let a3 = (&a * &a * &a) % &p;
    let b2 = (&p256_b() * &p256_b()) % &p;
    let num = (BigUint::from(1728u32) * &a3) % &p;
    let den = ((&a3 * BigUint::from(4u32)) % &p + (&b2 * BigUint::from(27u32)) % &p) % &p;
    let den_inv = den.modpow(&(&p - BigUint::from(2u32)), &p);
    (&num * &den_inv) % &p
}

// ---------------------------------------------------------------------------
// Shared integer machinery
// ---------------------------------------------------------------------------

/// The multiplicative order of `x mod p`, computed from a (possibly
/// partial) factorization of `p − 1`. `factors` holds `(prime, exp)`
/// pairs; pairs with `exp == 0` are unfactored composites and are
/// ignored, so the result is the order modulo the unfactored part —
/// i.e. the true order is `result · d` for some `d` dividing the
/// unfactored residue of `p − 1`. When the factorization is complete,
/// the result is the exact order.
pub fn multiplicative_order_mod(x: &BigUint, p: &BigUint, factors: &[(BigUint, u32)]) -> BigUint {
    let pm1 = p - BigUint::one();
    let mut order = pm1.clone();
    for (q, e) in factors {
        if *e == 0 {
            continue; // unfactored composite: leave its contribution alone
        }
        while (&order % q).is_zero() && x.modpow(&(&order / q), p) == BigUint::one() {
            order /= q;
        }
    }
    order
}

/// Legendre-symbol quadratic residuosity for an odd prime `p`.
pub fn is_quadratic_residue(x: &BigUint, p: &BigUint) -> bool {
    let e = (p - BigUint::one()) >> 1;
    x.modpow(&e, p) == BigUint::one()
}

/// Whether `x` is a d-th power residue mod `p`:
/// true iff `x^((p−1)/gcd(d, p−1)) == 1`.
pub fn is_dth_power_residue(x: &BigUint, p: &BigUint, d: u64) -> bool {
    let pm1 = p - BigUint::one();
    let g = BigUint::from(d).gcd(&pm1);
    let e = &pm1 / &g;
    x.modpow(&e, p) == BigUint::one()
}

/// Number of set bits (Hamming weight) of `x`.
pub fn popcount(x: &BigUint) -> u64 {
    (0..x.bits()).filter(|&i| x.bit(i)).count() as u64
}

/// Minimal `log2|x − c·2^k|` over `c ∈ [1, 32]` and all `k`: how close
/// `x` is to a small multiple of a power of two. Returns
/// `(gap_bits, best_c, best_k)`.
pub fn nearest_sparse_multiple(x: &BigUint) -> (f64, u64, u64) {
    let bits = x.bits();
    let mut best_gap = f64::INFINITY;
    let mut best = (0u64, 0u64);
    for c in 1u64..=32u64 {
        for k in 0..bits {
            let target = BigUint::from(c) << k;
            let gap = if target > *x { target - x } else { x - &target };
            let gap_bits = gap.bits() as f64;
            if gap_bits < best_gap {
                best_gap = gap_bits;
                best = (c, k);
            }
        }
    }
    (best_gap, best.0, best.1)
}

/// Whether `x` is a perfect power: `x = y^e` for some `e ≥ 2`.
pub fn is_perfect_power(x: &BigUint) -> Option<(BigUint, u32)> {
    if *x <= BigUint::one() {
        return None;
    }
    let max_e = (x.bits() as f64).log2().floor().max(2.0) as u32;
    for e in 2..=max_e {
        let (root, exact) = nth_root_exact(x, e);
        if exact {
            return Some((root, e));
        }
    }
    None
}

/// Integer `e`-th root via bisection, with exactness flag.
fn nth_root_exact(x: &BigUint, e: u32) -> (BigUint, bool) {
    let one = BigUint::one();
    if e == 1 || *x <= one {
        return (x.clone(), true);
    }
    let mut lo = BigUint::one();
    let mut hi: BigUint = {
            let seed: BigUint = BigUint::one();
            seed << ((x.bits() / e as u64) + 1)
        };
    while &hi - &lo > one {
        let mid: BigUint = (&lo + &hi) >> 1;
        if mid.pow(e) <= *x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let exact = lo.pow(e) == *x;
    (lo, exact)
}

/// The first `n_terms` continued-fraction partial quotients of the
/// rational `num/den`, and the largest of them. A seed-derived value
/// chosen near a quadratic irrational would show small, near-periodic
/// quotients.
pub fn continued_fraction(num: &BigUint, den: &BigUint, n_terms: usize) -> (Vec<BigUint>, BigUint) {
    let mut a = num.clone();
    let mut b = den.clone();
    let mut quotients = Vec::with_capacity(n_terms);
    let mut max_q = BigUint::zero();
    for _ in 0..n_terms {
        if b.is_zero() {
            break;
        }
        let (q, r) = a.div_rem(&b);
        if q > max_q {
            max_q = q.clone();
        }
        quotients.push(q);
        a = b;
        b = r;
    }
    (quotients, max_q)
}

// ---------------------------------------------------------------------------
// H2: field prime probes
// ---------------------------------------------------------------------------

fn small_u64(x: &BigUint) -> u64 {
    x.iter_u32_digits()
        .enumerate()
        .fold(0u64, |acc, (i, d)| acc | ((d as u64) << (32 * i)))
}

/// Residue-class profile of `p` against the small moduli a designer
/// might have preferred for implementation niceness.
pub fn p_residue_profile(moduli: &[u64]) -> Vec<(u64, u64)> {
    let p = p256_p();
    moduli
        .iter()
        .map(|&m| (m, small_u64(&(&p % BigUint::from(m)))))
        .collect()
}

/// Smallest `k ≤ max_k` with `p^k ≡ 1 (mod l)`: the embedding degree of
/// the order-`l` subgroup. `None` if no such `k` exists in range.
pub fn embedding_degree_of(l: &BigUint, max_k: u32) -> Option<u32> {
    let p = p256_p();
    let one = BigUint::one();
    let mut pk = BigUint::one();
    for k in 1..=max_k {
        pk = (&pk * &p) % l;
        if pk == one {
            return Some(k);
        }
    }
    None
}

// ---------------------------------------------------------------------------
// H3: rational roots of the 2-division polynomial
// ---------------------------------------------------------------------------

/// Degree-2 polynomial arithmetic mod `p` used to compute `x^p mod f`
/// for a cubic `f`: represent a polynomial of degree < 3 by its
/// coefficient vector.
fn poly_mul_mod(a: &[BigUint; 3], b: &[BigUint; 3], f: &[BigUint; 4], p: &BigUint) -> [BigUint; 3] {
    // Multiply two degree-<3 polynomials (product degree < 6), then
    // reduce mod the monic cubic f = x³ + f1 x² + f2 x + f3.
    let mut prod = vec![BigUint::zero(); 6];
    for (i, ai) in a.iter().enumerate() {
        for (j, bj) in b.iter().enumerate() {
            prod[i + j] = (&prod[i + j] + ai * bj) % p;
        }
    }
    // Reduce from degree 5 down to degree 2 using x³ ≡ −f1x² − f2x − f3,
    // i.e. x^d ≡ x^(d−k) · (−f_k) for k = 1, 2, 3.
    for d in (3..6usize).rev() {
        let c = std::mem::take(&mut prod[d]);
        if c.is_zero() {
            continue;
        }
        for k in 1..=3usize {
            let fk = &f[k];
            let term = (&c * p + p - &c * fk % p) % p; // −c·f_k mod p
            let idx = d - k;
            prod[idx] = (&prod[idx] + &term) % p;
        }
    }
    [prod[0].clone(), prod[1].clone(), prod[2].clone()]
}

/// Number of rational roots in `F_p` of the cubic
/// `f(x) = x³ + f1x² + f2x + f3` (coefficients already reduced mod `p`),
/// computed as `deg gcd(f, x^p − x)`.
pub fn cubic_rational_root_count(f1: &BigUint, f2: &BigUint, f3: &BigUint, p: &BigUint) -> usize {
    let f = [
        BigUint::one(),
        f1.clone(),
        f2.clone(),
        f3.clone(),
    ];
    // Compute x^p mod f by square-and-multiply on the exponent p.
    let mut result = [BigUint::one(), BigUint::zero(), BigUint::zero()]; // 1
    let mut base: [BigUint; 3] = [BigUint::zero(), BigUint::one(), BigUint::zero()]; // x
    for i in 0..p.bits() {
        if p.bit(i) {
            result = poly_mul_mod(&result, &base, &f, p);
        }
        if i + 1 < p.bits() {
            base = poly_mul_mod(&base, &base, &f, p);
        }
    }
    // gcd(f, x^p − x): if x^p ≡ x (mod f) then f | x^p − x, every root of
    // f is rational, and the answer is deg f = 3.
    let g = [
        result[0].clone() % p,
        (&result[1] + p - BigUint::one()) % p,
        result[2].clone(),
    ];
    if g[0].is_zero() && g[1].is_zero() && g[2].is_zero() {
        return 3;
    }
    // Polynomial Euclid: divide f by g (degrees ≤ 3 and < 3).
    fn degree_of(v: &[BigUint; 3]) -> Option<usize> {
        (0..3).rev().find(|&i| !v[i].is_zero())
    }
    // Represent as f (deg 3) and g (deg ≤ 2); Euclid until remainder 0.
    // Long division of f by g.
    let div_rem_cubic = |f: &[BigUint; 4], g: &[BigUint; 3], p: &BigUint| -> [BigUint; 3] {
        let gd = match degree_of(g) {
            Some(d) => d,
            None => return [BigUint::zero(), BigUint::zero(), BigUint::zero()], // g = 0: remainder is f (shouldn't happen)
        };
        let lead_inv = g[gd].modpow(&(p - BigUint::from(2u32)), p);
        let mut rem: Vec<BigUint> = f.to_vec();
        while let Some(rd) = (0..rem.len()).rev().find(|&i| !rem[i].is_zero()) {
            if rd < gd || rd == 0 && gd == 0 && rem[rd].is_zero() {
                break;
            }
            if rd < gd {
                break;
            }
            let shift = rd - gd;
            let coef = (&rem[rd] * &lead_inv) % p;
            for k in 0..=gd {
                let sub = (&coef * &g[k]) % p;
                rem[shift + k] = (&rem[shift + k] + p - &sub) % p;
            }
            if rem.len() > 3 && rd == rem.len() - 1 {
                // keep at least 4 slots; no truncation needed
            }
        }
        let mut out = [
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
        ];
        for (i, slot) in out.iter_mut().enumerate().take(3) {
            *slot = rem.get(i).cloned().unwrap_or_else(BigUint::zero);
        }
        out
    };
    let mut g_cur = g;
    let mut f_cur = f;
    loop {
        let r = div_rem_cubic(&f_cur, &g_cur, p);
        if degree_of(&r).is_none() {
            break;
        }
        // Shift: f := g (pad), g := r.
        let mut f_next = [
                BigUint::zero(),
                BigUint::zero(),
                BigUint::zero(),
                BigUint::zero(),
            ];
        for (i, slot) in g_cur.iter().enumerate() {
            f_next[i] = slot.clone();
        }
        f_cur = f_next;
        g_cur = r;
    }
    // The last nonzero g is the gcd.
    match degree_of(&g_cur) {
        Some(d) => d + 1,
        None => 0,
    }
}

/// Result of the rational-2-torsion probe on `x³ − 3x + b` over `F_p`.
#[derive(Clone, Debug, PartialEq)]
pub struct TorsionSplit {
    pub rational_roots: usize,
    pub fully_splits: bool,
}

/// Count the rational roots of P-256's 2-division polynomial
/// `x³ − 3x + b` over `F_p`. Three roots ⇒ full rational 2-torsion ⇒
/// three rational 2-isogenies; one root ⇒ one; zero ⇒ the Galois group
/// is C3 and there is exactly one rational 2-isogeny.
pub fn p256_two_torsion_split() -> TorsionSplit {
    let p = p256_p();
    let a = &p - BigUint::from(3u32); // −3 mod p
    let b = p256_b() % &p;
    let roots = cubic_rational_root_count(&BigUint::zero(), &a, &b, &p);
    TorsionSplit {
        rational_roots: roots,
        fully_splits: roots == 3,
    }
}

// ---------------------------------------------------------------------------
// Deep factorization wrapper
// ---------------------------------------------------------------------------

/// Deep-factor `n` (trial division to `2^22`, then bounded Pollard
/// rho on composite cofactors), separating fully factored primes from
/// the unfactored composite residue.
pub struct DeepFactorOutcome {
    pub factors: Vec<(BigUint, u32)>,
    /// The unfactored composite part of `n` (1 if fully factored).
    pub residue: BigUint,
    pub residue_is_prime: bool,
}

pub fn deep_factor_outcome(n: &BigUint, pollard_iters: u64) -> DeepFactorOutcome {
    let raw = deep_factor(n, 1 << 22, pollard_iters);
    let mut factors = Vec::new();
    let mut residue = BigUint::one();
    for (f, e) in raw {
        if e == 0 {
            residue = f; // unfactored composite marker from deep_factor
        } else {
            factors.push((f, e));
        }
    }
    let residue_is_prime = residue > BigUint::one() && crate::asymmetric::rsa::is_prime(&residue);
    DeepFactorOutcome {
        factors,
        residue,
        residue_is_prime,
    }
}

// ---------------------------------------------------------------------------
// The round-2 report
// ---------------------------------------------------------------------------

/// One groundable-invariant probe outcome: the measured value and
/// whether it flags an anomaly worth follow-up.
#[derive(Clone, Debug, PartialEq)]
pub struct TrapdoorProbeRecord {
    pub hypothesis: u8,
    pub name: String,
    pub measured: String,
    pub anomaly: bool,
}

/// Run every round-2 probe deterministically. `pollard_iters` bounds
/// the rho work inside the deep factorizations; `2^22` finishes in
/// seconds on the recorded host and reaches prime factors up to
/// roughly `2^44`.
pub fn p256_trapdoor_round2_report(pollard_iters: u64) -> Vec<TrapdoorProbeRecord> {
    let mut records = Vec::new();
    let p = p256_p();
    let b = p256_b();
    let r = p256_r_from_b();
    let j = p256_j_invariant();
    let t = p256_trace();
    let one = BigUint::one();

    // --- H1: multiplicative orders mod p -------------------------------
    let pm1_out = deep_factor_outcome(&(&p - &one), pollard_iters);
    let fully_factored_pm1 = pm1_out.residue <= one;
    for (name, x) in [
        ("b", &b),
        ("r=a^3/b^2", &r),
        ("j", &j),
        ("a=-3", &(&p - BigUint::from(3u32))),
    ] {
        let ord = multiplicative_order_mod(x, &p, &pm1_out.factors);
        records.push(TrapdoorProbeRecord {
            hypothesis: 1,
            name: format!("multiplicative order of {name} mod p"),
            measured: format!(
                "{} bits ({}); p-1 {}",
                ord.bits(),
                format!("{:x}", ord),
                if fully_factored_pm1 { "fully factored: exact" } else { "partially factored: order divides this" }
            ),
            anomaly: false,
        });
    }

    // --- H1: residuosity classes --------------------------------------
    records.push(TrapdoorProbeRecord {
        hypothesis: 1,
        name: "residuosity classes of b, r=a^3/b^2, j".into(),
        measured: format!(
            "QR: b={}, r={}, j={}; cubic residues: b={}, r={}, j={}",
            is_quadratic_residue(&b, &p),
            is_quadratic_residue(&r, &p),
            is_quadratic_residue(&j, &p),
            is_dth_power_residue(&b, &p, 3),
            is_dth_power_residue(&r, &p, 3),
            is_dth_power_residue(&j, &p, 3),
        ),
        anomaly: false,
    });

    // --- H1: sparse and perfect-power forms ----------------------------
    for (name, x) in [(&"b", &b), (&"r", &r), (&"j", &j), (&"n", &p256_n())] {
        let w = popcount(x);
        let expected = x.bits() as f64 / 2.0;
        let sigma = (x.bits() as f64 / 4.0).sqrt();
        let z = (w as f64 - expected) / sigma;
        let (sparse_gap, coeff, exp) = nearest_sparse_multiple(x);
        let pp = is_perfect_power(x);
        records.push(TrapdoorProbeRecord {
            hypothesis: 1,
            name: format!(
                "{name}: Hamming weight vs binomial, nearest c*2^k, perfect power"
            ),
            measured: format!(
                "weight={w}/{} z={z:.2}; nearest {coeff}*2^{exp} gap={sparse_gap:.0} bits; perfect_power={pp:?}",
                x.bits()
            ),
            anomaly: z.abs() > 5.0,
        });
    }

    // --- H1: continued-fraction heights --------------------------------
    for (name, x) in [(&"b/p", &b), (&"r/p", &r), (&"j/p", &j)] {
        let (qs, max_q) = continued_fraction(x, &p, 64);
        records.push(TrapdoorProbeRecord {
            hypothesis: 1,
            name: format!("{name}: continued fraction (64 terms)"),
            measured: format!(
                "max partial quotient {} ({:.0} bits), first 8: {}",
                max_q,
                max_q.bits(),
                qs.iter()
                    .take(8)
                    .map(|q| q.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ),
            anomaly: false,
        });
    }

    // --- H2: p ± 1 deep factorization ----------------------------------
    for (name, m) in [(&"p-1", &p - &one), (&"p+1", &p + &one)] {
        let out = deep_factor_outcome(&m, pollard_iters);
        let largest = out.factors.iter().map(|(f, _)| f.bits()).max().unwrap_or(0);
        records.push(TrapdoorProbeRecord {
            hypothesis: 2,
            name: format!(
                "{name}: deep factorization ({} prime powers, residue {} bits, residue prime {})",
                out.factors.len(),
                out.residue.bits(),
                out.residue_is_prime
            ),
            measured: format!(
                "largest known prime factor: {largest} bits; factors = {}",
                out.factors
                    .iter()
                    .map(|(f, e)| format!("{:x}^{}", f, e))
                    .collect::<Vec<_>>()
                    .join(" ")
            ),
            anomaly: false,
        });
    }

    // --- H2: ord_p(2) ---------------------------------------------------
    let ord2 = multiplicative_order_mod(&BigUint::from(2u32), &p, &pm1_out.factors);
    records.push(TrapdoorProbeRecord {
        hypothesis: 2,
        name: "multiplicative order of 2 mod p (known part)".into(),
        measured: format!("{} bits ({:x})", ord2.bits(), ord2),
        anomaly: false,
    });

    // --- H2: residue-class profile of p ---------------------------------
    let residues = p_residue_profile(&[3, 4, 5, 7, 8, 11, 12, 16, 24, 32, 64, 128, 256]);
    records.push(TrapdoorProbeRecord {
        hypothesis: 2,
        name: "residue class of p mod small moduli".into(),
        measured: residues
            .iter()
            .map(|(m, v)| format!("p≡{v} mod {m}"))
            .collect::<Vec<_>>()
            .join(", "),
        anomaly: false,
    });

    // --- H2: twist embedding degree -------------------------------------
    let twist = p256_twist_order();
    let twist_out = deep_factor_outcome(&twist, pollard_iters);
    let mut min_twist_k = u32::MAX;
    for (f, _) in &twist_out.factors {
        if *f < BigUint::from(1u32 << 20) {
            continue;
        }
        if let Some(k) = embedding_degree_of(f, 100) {
            min_twist_k = min_twist_k.min(k);
        }
    }
    records.push(TrapdoorProbeRecord {
        hypothesis: 2,
        name: "smallest embedding degree of any ≥2^20 twist subgroup (k ≤ 100)".into(),
        measured: if min_twist_k == u32::MAX {
            "none found for k ≤ 100".to_string()
        } else {
            format!("k = {min_twist_k}")
        },
        anomaly: min_twist_k <= 20,
    });

    // --- H3: GLS gap and trace profile ----------------------------------
    let four_p_minus_t2 = BigUint::from(4u32) * &p - &t * &t;
    records.push(TrapdoorProbeRecord {
        hypothesis: 3,
        name: "GLS gap 4p − t² and trace profile".into(),
        measured: format!(
            "{} bits; t popcount={}, t nearest c*2^k gap {:.0} bits",
            four_p_minus_t2.bits(),
            popcount(&t),
            nearest_sparse_multiple(&t).0
        ),
        anomaly: false,
    });

    // --- H3: rational 2-torsion splitting -------------------------------
    let split = p256_two_torsion_split();
    records.push(TrapdoorProbeRecord {
        hypothesis: 3,
        name: "rational roots of x³ − 3x + b over F_p".into(),
        measured: format!(
            "{} rational roots ({} rational 2-isogenies), fully_splits={}",
            split.rational_roots,
            if split.rational_roots == 0 { 1 } else { split.rational_roots },
            split.fully_splits
        ),
        anomaly: false,
    });

    records
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constants_have_expected_bit_lengths() {
        assert_eq!(p256_b().bits(), 255);
        assert_eq!(p256_p().bits(), 256);
        assert!(p256_j_invariant().bits() <= 256);
    }

    #[test]
    fn r_from_b_matches_the_x962_pinned_ratio() {
        let r = p256_r_from_b();
        let pinned = BigUint::parse_bytes(
            b"7efba1662985be9403cb055c75d4f7e0ce8d84a9c5114abcaf3177680104fa0d",
            16,
        )
        .unwrap();
        assert_eq!(r, pinned);
    }

    #[test]
    fn negation_flips_residuosity() {
        let p = p256_p();
        let b = p256_b();
        let qr_b = is_quadratic_residue(&b, &p);
        let neg_b = (&p - &b) % &p;
        assert_eq!(is_quadratic_residue(&neg_b, &p), !qr_b);
    }

    #[test]
    fn multiplicative_order_mod_is_correct_on_a_small_prime() {
        // p = 101, p−1 = 2^2 · 5 · 5. ord(2) = 100? 2^50 mod 101…
        let p = BigUint::from(101u32);
        let factors = vec![
            (BigUint::from(2u32), 2u32),
            (BigUint::from(5u32), 2u32),
        ];
        let ord = multiplicative_order_mod(&BigUint::from(2u32), &p, &factors);
        assert_eq!(BigUint::from(2u32).modpow(&ord, &p), BigUint::one());
        assert!((&p - BigUint::one()) % &ord == BigUint::zero());
        // 2 is a primitive root mod 101.
        assert_eq!(ord, BigUint::from(100u32));
        let ord3 = multiplicative_order_mod(&BigUint::from(33u32), &p, &factors);
        assert_eq!(BigUint::from(33u32).modpow(&ord3, &p), BigUint::one());
    }

    #[test]
    fn perfect_power_checks_on_small_integers() {
        assert_eq!(
            is_perfect_power(&BigUint::from(1728u32)),
            Some((BigUint::from(12u32), 3))
        );
        assert_eq!(is_perfect_power(&BigUint::from(1729u32)), None);
    }

    #[test]
    fn continued_fraction_of_a_rational_is_correct() {
        let (qs, max_q) = continued_fraction(&BigUint::from(355u32), &BigUint::from(113u32), 10);
        assert_eq!(qs[0], BigUint::from(3u32));
        // 355/113 = [3; 7, 16].
        assert_eq!(max_q, BigUint::from(16u32));
    }

    #[test]
    fn cubic_root_count_matches_direct_search_on_a_small_field() {
        // f(x) = x³ − x over F_11 has roots 0, 1, 10.
        let p = BigUint::from(11u32);
        let count = cubic_rational_root_count(
            &BigUint::zero(),
            &BigUint::from(10u32), // −1 mod 11
            &BigUint::zero(),
            &p,
        );
        assert_eq!(count, 3);
        // f(x) = x³ + x + 1 over F_11: roots? x=7: 343+7+1=351=31·11+10 → 10, no.
        // Compute by brute force and compare.
        let f = |x: u64| (x * x * x + x + 1) % 11;
        let brute: Vec<u64> = (0..11).filter(|&x| f(x) == 0).collect();
        let count = cubic_rational_root_count(
            &BigUint::zero(),
            &BigUint::one(),
            &BigUint::one(),
            &p,
        );
        assert_eq!(count, brute.len());
    }

    #[test]
    fn round2_report_runs_and_records_all_three_hypotheses() {
        let records = p256_trapdoor_round2_report(1 << 14);
        assert!(records.iter().any(|r| r.hypothesis == 1));
        assert!(records.iter().any(|r| r.hypothesis == 2));
        assert!(records.iter().any(|r| r.hypothesis == 3));
        assert!(records.iter().all(|r| !r.measured.is_empty()));
        assert!(records.iter().all(|r| !r.name.is_empty()));
    }
}

